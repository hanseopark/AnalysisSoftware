// provided by Gamma Conversion Group, $ALICE_PHYSICS/PWGGA/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion
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


void PrepareCocktail(   TString     nameFileCocktail            = "",
                        TString     nameFilePi0                 = "",
                        TString     suffix                      = "eps",
                        TString     cutSelection                = "",
                        TString     option                      = "",
                        TString     directphotonPlots           = "",
                        Double_t    rapidity                    = 0.85,
                        TString     period                      = "",
                        Int_t       numberOfBins                = 30,
                        Int_t       mode                        = 0,
                        Bool_t      producePlotsInOrPtRange     = kFALSE,
                        Bool_t      producePlotsForThesis       = kFALSE
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
    fdirectphoton                                               = directphotonPlots;
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

    //***************************** Calculate scaling factor ********************************************************
    eventNormScalingFactor                                      = ReturnCocktailNormalization(fEnergyFlag, fEventCutSelection);
    
    //**************************** Determine Centrality *************************************************************
    TString centrality                                          = GetCentralityString(fEventCutSelection);
    TString collisionSystem                                     = ReturnFullCollisionsSystem(option);
    TString cent                                                = "";
    TString textMeasurement                                     = ""; //"#gamma";
    if(option.Contains("PbPb")) cent                            = Form("%s %s", centrality.Data(), collisionSystem.Data());
    else                        cent                            = collisionSystem;
    
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
            cocktailInputParametrizationsMtScaled[i]            = (TF1*)MtScaledParam(cocktailInputParametrizations[0], motherParticlesPDG[i].Atoi(), motherParticlesPDG[0].Atoi(), mtScaleFactor[i], kFALSE, kTRUE);
    }
  
    if (cocktailSettingsList->FindObject("2212_pt"))
        cocktailInputParametrizationProton                      = (TF1*)cocktailSettingsList->FindObject(Form("2212_pt"));
    
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

    //***************************** ranges and scaling factors ******************************************************
    Double_t deltaPtGen                                         = ptGenMax-ptGenMin;
    Double_t deltaRap                                           = 2*rapidity;
    Double_t deltaEta                                           = ReturnDeltaEta(fGammaCutSelection);   // 2*0.9
    Double_t eta                                                = deltaEta*0.5;                         // 0.9
    Double_t deltaPhi                                           = 2*TMath::Pi();
    Double_t deltaEtaCalo                                       = 0;
    Double_t deltaPhiCalo                                       = 0;
    Double_t minPhiCalo                                         = 0;
    Double_t maxPhiCalo                                         = 0;
    Double_t etaCalo                                            = 0;
    if (isCalo){
        deltaEtaCalo                                            = ReturnDeltaEtaCalo(fClusterCutSelection);
        etaCalo                                                 = deltaEtaCalo*0.5;
        deltaPhiCalo                                            = ReturnDeltaPhiCalo(fClusterCutSelection);
        TString phiMinCut(fClusterCutSelection(                 GetClusterPhiMinCutPosition(fClusterCutSelection),1));
        TString phiMaxCut(fClusterCutSelection(                 GetClusterPhiMaxCutPosition(fClusterCutSelection),1));
        minPhiCalo                                              = AnalyseClusterMinPhiCut(phiMinCut.Atoi());
        maxPhiCalo                                              = AnalyseClusterMaxPhiCut(phiMaxCut.Atoi());
    }
    
    cout << "========================================"  << endl;
    cout << "deltaRap         = "     << deltaRap       << endl;
    cout << "deltaPtGen       = "     << deltaPtGen     << endl;
    cout << "deltaPt          = "     << ptMax - ptMin  << endl;
    cout << "deltaEta         = "     << deltaEta       << endl;
    if (isCalo)
        cout << "deltaEtaCalo     = " << deltaEtaCalo   << endl;
    cout << "deltaPhi         = "     << deltaPhi       << endl;
    if (isCalo)
        cout << "deltaPhiCalo     = " << deltaPhiCalo   << endl;
    cout << "========================================"  << endl;

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
    if (histoPi0YieldData) cout << "Found the corrected pi0 yield in " << nameFilePi0.Data() << ", will produce cocktail and data pi0 plot." << endl;

    //***************************** Eta file ************************************************************************
    TString nameFileEta                                         = nameFilePi0;
    nameFileEta.ReplaceAll("Pi0_", "Eta_");
    TFile* fileEta                                              = new TFile(nameFileEta);
    if (fileEta && !fileEta->IsZombie()) {
        histoEtaYieldData                                       = (TH1D*)fileEta->Get("CorrectedYieldTrueEff");
    } else {
        histoEtaYieldData                                       = NULL;
    }
    if (histoEtaYieldData) cout << "Found the corrected eta yield in " << nameFileEta.Data() << ", will produce cocktail and data eta plot." << endl;

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
        histoGammaSumPtY->Scale(eventNormScalingFactor);
        histoGammaSumPtPhi->Sumw2();
        histoGammaSumPtPhi->Scale(eventNormScalingFactor);
    } else {
        cout << "ERROR: Gamma histo not found!" << endl;
        return;
    }
    
    histoDecayChannels                                          = new TH1F*[nMotherParticles];
    histoGammaPtY                                               = new TH2F*[nMotherParticles];
    histoGammaPtPhi                                             = new TH2F*[nMotherParticles];
    histoGammaMotherPtY                                         = new TH2F*[nMotherParticles];
    histoGammaMotherPtPhi                                       = new TH2F*[nMotherParticles];
    histoGammaMotherPtGammaPt                                   = new TH2F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (hasMother[i]) {
            histoDecayChannels[i]                               = (TH1F*)histoListCocktail->FindObject(Form("DecayChannels_%s", motherParticles[i].Data()));
            
            histoGammaPtY[i]                                    = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_Gamma_From_%s", motherParticles[i].Data()));
            histoGammaPtY[i]->SetName(Form("Gamma_From_%s_Pt_Y_OrBin", motherParticles[i].Data()));
            histoGammaPtY[i]->Sumw2();
            histoGammaPtY[i]->Scale(eventNormScalingFactor);
            
            histoGammaPtPhi[i]                                  = (TH2F*)histoListCocktail->FindObject(Form("Pt_Phi_Gamma_From_%s", motherParticles[i].Data()));
            histoGammaPtPhi[i]->SetName(Form("Gamma_From_%s_Pt_Phi_OrBin", motherParticles[i].Data()));
            histoGammaPtPhi[i]->Sumw2();
            histoGammaPtPhi[i]->Scale(eventNormScalingFactor);
            
            histoGammaMotherPtY[i]                              = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_%s", motherParticles[i].Data()));
            histoGammaMotherPtY[i]->SetName(Form("%s_Pt_Y_OrBin", motherParticles[i].Data()));
            histoGammaMotherPtY[i]->Sumw2();
            histoGammaMotherPtY[i]->Scale(eventNormScalingFactor);
            
            histoGammaMotherPtPhi[i]                            = (TH2F*)histoListCocktail->FindObject(Form("Pt_Phi_%s", motherParticles[i].Data()));
            histoGammaMotherPtPhi[i]->SetName(Form("%s_Pt_Phi_OrBin", motherParticles[i].Data()));
            histoGammaMotherPtPhi[i]->Sumw2();
            histoGammaMotherPtPhi[i]->Scale(eventNormScalingFactor);

            histoGammaMotherPtGammaPt[i]                        = (TH2F*)histoListCocktail->FindObject(Form("PtGamma_PtMother_%s", motherParticles[i].Data()));
            histoGammaMotherPtGammaPt[i]->SetName(Form("%s_PtGamma_PtMother_OrBin", motherParticles[i].Data()));
            histoGammaMotherPtGammaPt[i]->Sumw2();
            histoGammaMotherPtGammaPt[i]->Scale(eventNormScalingFactor);
        } else {
            histoDecayChannels[i]                               = NULL;
            histoGammaPtY[i]                                    = NULL;
            histoGammaPtPhi[i]                                  = NULL;
            histoGammaMotherPtY[i]                              = NULL;
            histoGammaMotherPtPhi[i]                            = NULL;
            histoGammaMotherPtGammaPt[i]                        = NULL;
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
    histoGammaSumPtOrBin->Scale(deltaPhi);
    histoGammaSumPtOrBin->GetXaxis()->SetRangeUser(ptMin, ptMax);
    histoGammaSumYOrBin                                         = (TH1F*)histoGammaSumPtY->ProjectionY("Gamma_Y_OrBin", 1, histoGammaSumPtY->GetNbinsX(), "e");
    SetHistogramTitles(histoGammaSumYOrBin,"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
    histoGammaSumYOrBin->Sumw2();
    histoGammaSumYOrBin->Scale(deltaPhi);
    histoGammaSumPhiOrBin                                       = (TH1F*)histoGammaSumPtPhi->ProjectionY("Gamma_Phi_OrBin", 1, histoGammaSumPtPhi->GetNbinsX(), "e");
    SetHistogramTitles(histoGammaSumPhiOrBin,"","#phi","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
    histoGammaSumPhiOrBin->Sumw2();
    histoGammaSumPhiOrBin->Scale(deltaPhi);
    if(doGammaPtOrBinFromPtPhi) {
        histoGammaSumPtOrBin2                                   = (TH1F*)histoGammaSumPtPhi->ProjectionX("Gamma_Pt_OrBin_2", 1, histoGammaSumPtY->GetNbinsY(), "e");
        SetHistogramTitles(histoGammaSumPtOrBin2,"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
        histoGammaSumPtOrBin2->Sumw2();
        histoGammaSumPtOrBin2->Scale(deltaPhi);
        histoGammaSumPtOrBin2->GetXaxis()->SetRangeUser(ptMin, ptMax);
    } else {
        histoGammaSumPtOrBin2                                   = NULL;
    }
    
    histoGammaPtOrBin                                           = new TH1F*[nMotherParticles];
    histoGammaPtOrBin2                                          = new TH1F*[nMotherParticles];
    histoGammaYOrBin                                            = new TH1F*[nMotherParticles];
    histoGammaPhiOrBin                                          = new TH1F*[nMotherParticles];
    histoGammaMotherPtOrBin                                     = new TH1F*[nMotherParticles];
    histoGammaMotherYOrBin                                      = new TH1F*[nMotherParticles];
    histoGammaMotherPhiOrBin                                    = new TH1F*[nMotherParticles];
    histoGammaMotherPtGammaOrBin                                = new TH1F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        // gamma
        if (histoGammaPtY[i]) {
            histoGammaPtOrBin[i]                                = (TH1F*)histoGammaPtY[i]->ProjectionX(Form("Gamma_From_%s_Pt_OrBin", motherParticles[i].Data()),1,histoGammaPtY[i]->GetNbinsY(),"e");
            SetHistogramTitles(histoGammaPtOrBin[i],"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoGammaPtOrBin[i]->Sumw2();
            histoGammaPtOrBin[i]->Scale(deltaPhi);
            histoGammaPtOrBin[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
            histoGammaYOrBin[i]                                 = (TH1F*)histoGammaPtY[i]->ProjectionY(Form("Gamma_From_%s_Y_OrBin", motherParticles[i].Data()),1,histoGammaPtY[i]->GetNbinsX(),"e");
            SetHistogramTitles(histoGammaYOrBin[i],"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoGammaYOrBin[i]->Sumw2();
            histoGammaYOrBin[i]->Scale(deltaPhi);
        } else {
            histoGammaPtOrBin[i]                                = NULL;
            histoGammaYOrBin[i]                                 = NULL;
        }
        if (histoGammaPtPhi[i]) {
            histoGammaPhiOrBin[i]                               = (TH1F*)histoGammaPtPhi[i]->ProjectionY(Form("Gamma_From_%s_Phi_OrBin",motherParticles[i].Data()),1,histoGammaPtPhi[i]->GetNbinsX(),"e");
            SetHistogramTitles(histoGammaPhiOrBin[i],"","#phi","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoGammaPhiOrBin[i]->Sumw2();
            histoGammaPhiOrBin[i]->Scale(deltaPhi);
            if (doGammaPtOrBinFromPtPhi) {
                histoGammaPtOrBin2[i]                           = (TH1F*)histoGammaPtPhi[i]->ProjectionX(Form("Gamma_From_%s_Pt_OrBin_2",motherParticles[i].Data()),1,histoGammaPtY[i]->GetNbinsY(),"e");
                SetHistogramTitles(histoGammaPtOrBin2[i],"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
                histoGammaPtOrBin2[i]->Sumw2();
                histoGammaPtOrBin2[i]->Scale(deltaPhi);
                histoGammaPtOrBin2[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
            } else {
                histoGammaPtOrBin2[i]                           = NULL;
            }
        } else {
            histoGammaPhiOrBin[i]                               = NULL;
            histoGammaPtOrBin2[i]                               = NULL;
        }
        
        //mother
        if (histoGammaMotherPtY[i]) {
            histoGammaMotherPtOrBin[i]                          = (TH1F*)histoGammaMotherPtY[i]->ProjectionX(Form("%s_Pt_OrBin",motherParticles[i].Data()),1,histoGammaMotherPtY[i]->GetNbinsY(),"e");
            SetHistogramTitles(histoGammaMotherPtOrBin[i],"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoGammaMotherPtOrBin[i]->Sumw2();
            histoGammaMotherPtOrBin[i]->Scale(deltaPhi);
            histoGammaMotherPtOrBin[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
            histoGammaMotherYOrBin[i]                           = (TH1F*)histoGammaMotherPtY[i]->ProjectionY(Form("%s_Y_OrBin",motherParticles[i].Data()),1,histoGammaMotherPtY[i]->GetNbinsX(),"e");
            SetHistogramTitles(histoGammaMotherYOrBin[i],"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoGammaMotherYOrBin[i]->Sumw2();
            histoGammaMotherYOrBin[i]->Scale(deltaPhi);
        } else {
            histoGammaMotherPtOrBin[i]                          = NULL;
            histoGammaMotherYOrBin[i]                           = NULL;
        }
        if (histoGammaMotherPtPhi[i]) {
            histoGammaMotherPhiOrBin[i]                         = (TH1F*)histoGammaMotherPtPhi[i]->ProjectionY(Form("%s_Phi_OrBin",motherParticles[i].Data()),1,histoGammaMotherPtPhi[i]->GetNbinsX(),"e");
            SetHistogramTitles(histoGammaMotherPhiOrBin[i],"","#phi","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoGammaMotherPhiOrBin[i]->Sumw2();
            histoGammaMotherPhiOrBin[i]->Scale(deltaPhi);
        } else
            histoGammaMotherPhiOrBin[i]                         = NULL;

        if (histoGammaMotherPtGammaPt[i]) {
            histoGammaMotherPtGammaOrBin[i]                     = (TH1F*)histoGammaMotherPtGammaPt[i]->ProjectionX(Form("%s_PtGamma_OrBin",motherParticles[i].Data()),1,histoGammaMotherPtGammaPt[i]->GetNbinsY(),"e");
            SetHistogramTitles(histoGammaMotherPtGammaOrBin[i],"","#it{p}_{T, #gamma} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoGammaMotherPtGammaOrBin[i]->Sumw2();
            histoGammaMotherPtGammaOrBin[i]->Scale(deltaPhi/2); // histogram filled once for each gamma from pi0 -> gamma gamma
            histoGammaMotherPtGammaOrBin[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
        } else
            histoGammaMotherPtGammaOrBin[i]                     = NULL;
    }
    
    if (!histoGammaMotherPtOrBin[0]) {
        cout << "ERROR: Didn't get pi0 pt spectrum, returning!" << endl;
        return;
    }

    //***************************** Scale spectra *******************************************************************
    histoGammaSumPtOrBin->Scale(                                1./nEvents);
    if (histoGammaSumPtOrBin2) histoGammaSumPtOrBin2->Scale(    1./nEvents);
    histoGammaSumYOrBin->Scale(                                 1./nEvents);
    histoGammaSumPhiOrBin->Scale(                               1./nEvents);
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPtOrBin[i])               histoGammaPtOrBin[i]->Scale(            1./nEvents);
        if (histoGammaPtOrBin2[i])              histoGammaPtOrBin2[i]->Scale(           1./nEvents);
        if (histoGammaYOrBin[i])                histoGammaYOrBin[i]->Scale(             1./nEvents);
        if (histoGammaPhiOrBin[i])              histoGammaPhiOrBin[i]->Scale(           1./nEvents);
        if (histoGammaMotherPtOrBin[i])         histoGammaMotherPtOrBin[i]->Scale(      1./nEvents);
        if (histoGammaMotherYOrBin[i])          histoGammaMotherYOrBin[i]->Scale(       1./nEvents);
        if (histoGammaMotherPhiOrBin[i])        histoGammaMotherPhiOrBin[i]->Scale(     1./nEvents);
        if (histoGammaMotherPtGammaOrBin[i])    histoGammaMotherPtGammaOrBin[i]->Scale( 1./nEvents);
    }

    //***************************** Rebin pt spectra ****************************************************************
    histoGammaSumPt                                             = (TH1F*)histoGammaSumPtOrBin->Clone("Gamma_Pt");
    histoGammaSumPt->Sumw2();
    histoGammaSumPt->GetXaxis()->SetRangeUser(ptMin, ptMax);
    RebinSpectrum(histoGammaSumPt,"");
    
    histoGammaPt                                                = new TH1F*[nMotherParticles];
    histoGammaMotherPt                                          = new TH1F*[nMotherParticles];
    histoGammaMotherPtGamma                                     = new TH1F*[nMotherParticles];
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
        if (histoGammaMotherPtGammaOrBin[i]) {
            histoGammaMotherPtGamma[i]                          = (TH1F*)histoGammaMotherPtGammaOrBin[i]->Clone(Form("%s_PtGamma", motherParticles[i].Data()));
            histoGammaMotherPtGamma[i]->Sumw2();
            histoGammaMotherPtGamma[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
            RebinSpectrum(histoGammaMotherPtGamma[i],"");
        } else {
            histoGammaMotherPtGamma[i]                          = NULL;
        }
    }

    //***************************** Transform yields ****************************************************************
    histoGammaSumPt                                             = ConvertYieldHisto(histoGammaSumPt);
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPt[i])        histoGammaPt[i]             = ConvertYieldHisto(histoGammaPt[i]);
        if (histoGammaMotherPt[i])  histoGammaMotherPt[i]       = ConvertYieldHisto(histoGammaMotherPt[i]);
        if (histoGammaMotherPtGamma[i])
            histoGammaMotherPtGamma[i]                          = ConvertYieldHisto(histoGammaMotherPtGamma[i]);
    }
    
    //***************************** generated eta yield in analyzed eta binning *************************************
    if (histoGammaMotherPtOrBin[1]) {
        histoGeneratedEtaPt                                     = (TH1F*)histoGammaMotherPtOrBin[1]->Clone(Form("%s_inEtaBinning",histoGammaMotherPtOrBin[1]->GetName()));
        histoGeneratedEtaPt->Sumw2();
        RebinSpectrum(histoGeneratedEtaPt, (TH1F*)histoEtaYieldData, "");
        histoGeneratedEtaPt                                     = ConvertYieldHisto(histoGeneratedEtaPt);
    }

    // adapt plotting range for original binned histograms
    if (!producePlotsInOrPtRange) {
        ptPlotMin = ptMin;
        ptPlotMax = ptMax;
    } else {
        ptPlotMin = ptGenMin;
        ptPlotMax = ptGenMax;
    }
    if (ptPlotMin == 0 || ptPlotMin < 3e-1)
        ptPlotMin = 3e-1;

    //***************************** Sanity check ********************************************************************
    TH1D*   histoGammaResummedPt                                = (TH1D*)histoGammaPt[0]->Clone("histoGammaResummedPt");
    histoGammaResummedPt->Sumw2();
    TH1D*   histoGammaResummedPtOrBin                           = (TH1D*)histoGammaPtOrBin[0]->Clone("histoGammaResummedPtOrBin");
    histoGammaResummedPtOrBin->Sumw2();
    for (Int_t i=1; i<nMotherParticles; i++) {
        if (histoGammaPt[i])        histoGammaResummedPt->Add(      histoGammaPt[i]);
        if (histoGammaPtOrBin[i])   histoGammaResummedPtOrBin->Add( histoGammaPtOrBin[i]);
    }

    TH1D*   differenceGammaSumPt                                = (TH1D*)histoGammaSumPt->Clone("differenceGammaSumPt");
    differenceGammaSumPt->Sumw2();
    differenceGammaSumPt->Add(histoGammaResummedPt, -1);
    TH1D*   differenceGammaSumPtOrBin                           = (TH1D*)histoGammaSumPtOrBin->Clone("differenceGammaSumPtOrBin");
    differenceGammaSumPtOrBin->Sumw2();
    differenceGammaSumPtOrBin->Add(histoGammaResummedPtOrBin, -1);

    Bool_t  isSane                                              = kTRUE;
    for (Int_t i=1; i<differenceGammaSumPt->GetNbinsX()+1; i++) {
        if ( !((0 <= (differenceGammaSumPt->GetBinContent(i) + differenceGammaSumPt->GetBinError(i))) && (0 >= (differenceGammaSumPt->GetBinContent(i) - differenceGammaSumPt->GetBinError(i)))) ) {
            cout << "WARNING: " << differenceGammaSumPt->GetName() << " at p_T[" << i << "] = " << differenceGammaSumPt->GetBinCenter(i) << " GeV/c is " << differenceGammaSumPt->GetBinContent(i) << " +/- " << differenceGammaSumPt->GetBinError(i) << endl;
            isSane                                              = kFALSE;
        }
    }
    for (Int_t i=1; i<differenceGammaSumPtOrBin->GetNbinsX()+1; i++) {
        if ( !((0 <= (differenceGammaSumPtOrBin->GetBinContent(i) + differenceGammaSumPtOrBin->GetBinError(i))) && (0 >= (differenceGammaSumPtOrBin->GetBinContent(i) - differenceGammaSumPtOrBin->GetBinError(i)))) ) {
            cout << "WARNING: " << differenceGammaSumPtOrBin->GetName() << " at p_T[" << i << "] = " << differenceGammaSumPtOrBin->GetBinCenter(i) << " GeV/c is " << differenceGammaSumPtOrBin->GetBinContent(i) << " +/- " << differenceGammaSumPtOrBin->GetBinError(i) << endl;
            isSane                                              = kFALSE;
        }
    }

    if (!isSane) {
        cout << "ERROR: summed contributions to cocktail are not equal to the cocktail, check the normalizations! Plots are produced but histograms are not written!" << endl;
    }

    TH1D* dummyHist                                             = NULL;
    //***************************** Plot cocktail mothers (pt) ******************************************************
    TCanvas *canvasMothers                                      = new TCanvas("canvasMothers","",1100,1200);
    DrawGammaCanvasSettings(canvasMothers, 0.15, 0.01, 0.01, 0.075);
    canvasMothers->SetLogy();
    canvasMothers->SetLogx();
    
    TLegend* legendMothers                                      = GetAndSetLegend2(0.2, 0.98-(0.04*nRows), 0.95, 0.98, 40, 6);
    legendMothers->SetBorderSize(0);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, ptPlotMin, ptPlotMax);
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", 1e-7, 6e2, 1.0, 1.5);
    dummyHist->SetLabelOffset(-0.015, "X");
    dummyHist->SetTitleOffset(0.8, "X");
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
    
    PutProcessLabelAndEnergyOnPlot(                 0.22, 0.22, 0.032, cent, textMeasurement, "", 42, 0.03);
    if (producePlotsForThesis) PutThisThesisLabel(  0.22, 0.17, 0.032, 0.03, 1.25, 42);
    else PutALICESimulationLabel(                   0.22, 0.17, 0.032, 0.03, 1.25, 42);
    
    canvasMothers->SaveAs(Form("%s/CocktailMothers_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendMothers;
    delete canvasMothers;

    //***************************** Plot cocktail mothers + input param (pt) ****************************************
    TCanvas *canvasMothersParam                                 = new TCanvas("canvasMothersParam","",1100,1200);
    DrawGammaCanvasSettings(canvasMothersParam, 0.15, 0.01, 0.01, 0.075);
    canvasMothersParam->SetLogy();
    canvasMothersParam->SetLogx();
    
    TLegend* legendMothersParam                                 = GetAndSetLegend2(0.2, 0.98-(0.04*nRows), 0.95, 0.98, 40, 6);
    legendMothersParam->SetBorderSize(0);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, ptPlotMin, ptPlotMax);
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", 1e-7, 6e2, 1.0, 1.5);
    dummyHist->SetLabelOffset(-0.015, "X");
    dummyHist->SetTitleOffset(0.8, "X");
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
    
    PutProcessLabelAndEnergyOnPlot(                 0.22, 0.22, 0.032, cent, textMeasurement, "", 42, 0.03);
    if (producePlotsForThesis) PutThisThesisLabel(  0.22, 0.17, 0.032, 0.03, 1.25, 42);
    else PutALICESimulationLabel(                   0.22, 0.17, 0.032, 0.03, 1.25, 42);

    canvasMothersParam->SaveAs(Form("%s/CocktailMothersInclParam_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendMothersParam;
    delete canvasMothersParam;

    //***************************** Plot integrated mother ratios to pi0 ********************************************

    cout << "PLOTTING: Integrated particle ratios" << endl;
    Double_t textSizeLabelsPixel                 = 54;
    TCanvas* canvasIntPartRatios       = new TCanvas("canvasIntPartRatios","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasIntPartRatios, 0.1, 0.01, 0.01, 0.125);

        Double_t textsizeLabels = 0;
        Double_t textsizeFac    = 0;
        if (canvasIntPartRatios->XtoPixel(canvasIntPartRatios->GetX2()) <canvasIntPartRatios->YtoPixel(canvasIntPartRatios->GetY1()) ){
            textsizeLabels      = (Double_t)textSizeLabelsPixel/canvasIntPartRatios->XtoPixel(canvasIntPartRatios->GetX2()) ;
            textsizeFac         = (Double_t)1./canvasIntPartRatios->XtoPixel(canvasIntPartRatios->GetX2()) ;
        } else {
            textsizeLabels      = (Double_t)textSizeLabelsPixel/canvasIntPartRatios->YtoPixel(canvasIntPartRatios->GetY1());
            textsizeFac         = (Double_t)1./canvasIntPartRatios->YtoPixel(canvasIntPartRatios->GetY1());
        }

    TH2F * histomTscalingPoints = (TH2F*) histMtScalingFactors->Clone("MtScalingClone1");
    DrawGammaSetMarker(histomTscalingPoints, 23, 1.5, kBlue+2 , kBlue+2);
    
    TH2F * histoIntPartRatios = (TH2F*) histMtScalingFactors->Clone("MtScalingClone2");
    TH2F * histoTheoPartRatios = (TH2F*) histMtScalingFactors->Clone("MtScalingClone3");
    Double_t theoRatios[nMotherParticles] = {1,0.2131528511,0.0293771568,0.2121738246,0.2488197353,
                                               0.2488197353,0.2488197353,0.0641642177,4.69567604091221e-7,
                                               0.028698447,0.0286394262,0.028698447,0.0285501625,0.2267748773,0.2659074433,0.,0.};
    Double_t intpi0 = 0;
    Double_t intpi0err = 0;
    Double_t intratio = 0;
    Double_t intratioerr = 0;
    Double_t intpart[nMotherParticles];
    Double_t intparterr[nMotherParticles];
    intpi0 = histoGammaMotherPtOrBin[0]->IntegralAndError(histoGammaMotherPtOrBin[0]->FindBin( 0.1), histoGammaMotherPtOrBin[0]->FindBin(16), intpi0err, "width");
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaMotherPtOrBin[i]) {
            intpart[i] = histoGammaMotherPtOrBin[i]->IntegralAndError(histoGammaMotherPtOrBin[i]->FindBin( ptGenMin), histoGammaMotherPtOrBin[i]->FindBin(ptGenMax), intparterr[i], "width");
            intratio = intpart[i]/intpi0;
            intratioerr = intratio*pow(pow(intparterr[i]/intpart[i],2)+pow(intpi0err/intpi0,2),0.5);
            histoIntPartRatios->SetBinContent(i+1,intratio);
            histomTscalingPoints->GetXaxis()->SetBinLabel(i+1,Form("%s/#pi^{0}",motherParticlesLatex[i].Data()));
            histomTscalingPoints->SetBinContent(i+1,mtScaleFactor[i]);
            histomTscalingPoints->SetBinError(i+1,0);
            histoIntPartRatios->SetBinError(i+1,intratioerr);
            histoTheoPartRatios->SetBinContent(i+1,theoRatios[i]);
            histoTheoPartRatios->SetBinError(i+1,0.001);
        }else {
            histomTscalingPoints->SetBinContent(i+1,-1);
            histomTscalingPoints->SetBinError(i+1,0);
            histoIntPartRatios->SetBinContent(i+1,-1);
            histoIntPartRatios->SetBinError(i+1,0);
            histoTheoPartRatios->SetBinContent(i+1,-1);
            histoTheoPartRatios->SetBinError(i+1,0);
        }
    }
    DrawGammaSetMarker(histoIntPartRatios, 21, 1.5, kRed+2 , kRed+2);
    DrawGammaSetMarker(histoTheoPartRatios, 24, 1.5, kGreen+2 , kGreen+2);
    histomTscalingPoints->GetXaxis()->SetLabelSize(0.8*textsizeLabels);
    histomTscalingPoints->GetYaxis()->SetRangeUser(0.,1.05);
    histomTscalingPoints->GetXaxis()->SetRangeUser(1.,17);
    histomTscalingPoints->Draw("p");
    histoIntPartRatios->Draw("pesame");
    histoTheoPartRatios->Draw("psame");

    TLegend* legendIntRatio = GetAndSetLegend2(0.67, 0.65, 0.9, 0.65+(textsizeLabels*4*0.9), textSizeLabelsPixel);
    legendIntRatio->AddEntry(histomTscalingPoints,"mT scaling","p");
    legendIntRatio->AddEntry(histoIntPartRatios,"int. ratios","p");
    legendIntRatio->AddEntry(histoTheoPartRatios,"theory ratios","p");
        legendIntRatio->Draw();

    canvasIntPartRatios->Update();
    canvasIntPartRatios->SaveAs(Form("%s/IntegratedRatios.%s",outputDir.Data(),suffix.Data()));

    //***************************** Plot ratio cocktail mothers to pi0 (pt) *****************************************
    TCanvas *canvasMothersRatio                                 = new TCanvas("canvasMothersRatio","",1100,1200);
    DrawGammaCanvasSettings(canvasMothersRatio, 0.08, 0.01, 0.01, 0.075);
    canvasMothersRatio->SetLogy();
    canvasMothersRatio->SetLogx();
    
    TLegend* legendMothersRatio                                 = GetAndSetLegend2(0.13, 0.98-(0.04*nRows), 0.95, 0.98, 40, 5);
    legendMothersRatio->SetBorderSize(0);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, ptPlotMin, ptPlotMax);
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "particle ratio", 1e-3, 10, 1.0, 0.9);
    dummyHist->SetLabelOffset(-0.015, "X");
    dummyHist->SetTitleOffset(0.8, "X");
    dummyHist->Draw();

    TH1F* tempRatio                                             = NULL;
    Bool_t isMtScaled                                           = kFALSE;
    for (Int_t i=1; i<nMotherParticles; i++) {
        //if (!cocktailInputParametrizations[i]) isMtScaled       = kTRUE;
        
        // baryons
        if (cocktailInputParametrizationProton && ((i>= 9 && i<=13) || i==16) && histoGammaMotherPtOrBin[i]) {
            tempRatio                                           = (TH1F*)histoGammaMotherPtOrBin[i]->Clone("tempRatio");
            tempRatio->Sumw2();
            tempRatio->Divide(cocktailInputParametrizationProton);
            DrawGammaSetMarker(                             tempRatio, cocktailMarker[i], 1,  cocktailColor[i],  cocktailColor[i]);
            if(isMtScaled)  legendMothersRatio->AddEntry(   tempRatio, Form("%s / p *",       motherParticlesLatex[i].Data()), "l"); //from m_{T} scaling
            else            legendMothersRatio->AddEntry(   tempRatio, Form("%s / p",         motherParticlesLatex[i].Data()), "l");
            tempRatio->SetLineWidth(2);
            tempRatio->Draw("csamehist");
        }
      
        // mesons
        else if (histoGammaMotherPtOrBin[0] && histoGammaMotherPtOrBin[i]) {
            tempRatio                                           = (TH1F*)histoGammaMotherPtOrBin[i]->Clone("tempRatio");
            tempRatio->Sumw2();
            tempRatio->Divide(histoGammaMotherPtOrBin[0]);
            DrawGammaSetMarker(                             tempRatio, cocktailMarker[i], 1,  cocktailColor[i],  cocktailColor[i]);
            if(isMtScaled)  legendMothersRatio->AddEntry(   tempRatio, Form("%s / #pi^{0} *", motherParticlesLatex[i].Data()), "l");
            else            legendMothersRatio->AddEntry(   tempRatio, Form("%s / #pi^{0}",   motherParticlesLatex[i].Data()), "l");
            tempRatio->SetLineWidth(2);
            if(i!=0 && i%2==0) tempRatio->SetLineStyle(9);
            tempRatio->Draw("csamehist");
        }
      
        //isMtScaled                                              = kFALSE;
    }
    legendMothersRatio->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(                 0.7, 0.22, 0.032, cent, textMeasurement, "", 42, 0.03);
    if (producePlotsForThesis) PutThisThesisLabel(  0.7, 0.17, 0.032, 0.03, 1.25, 42);
    else PutALICESimulationLabel(                   0.7, 0.17, 0.032, 0.03, 1.25, 42);

    canvasMothersRatio->SaveAs(Form("%s/CocktailMothersRatio_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete tempRatio;
    delete legendMothersRatio;
    delete canvasMothersRatio;

    //***************************** Plot cocktail mothers (y) *******************************************************
    TCanvas *canvasMothersY                                     = new TCanvas("canvasMothersY","",1100,1200);
    DrawGammaCanvasSettings(canvasMothersY, 0.12, 0.025, 0.01, 0.075);
    canvasMothersY->SetLogy();

    TLegend* legendMothersY                                     = GetAndSetLegend2(0.2, 0.98-(0.04*nRows), 0.95, 0.98, 40, 6);
    legendMothersY->SetBorderSize(0);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, -fRapidity, fRapidity);
    SetHistogramm(dummyHist, "y", "#frac{1}{N_{ev}} #frac{d#it{N}}{dy}", 1e-4, 20, 0.9, 1.25);
    dummyHist->Draw();
    
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaMotherYOrBin[i]) {
            DrawGammaSetMarker(         histoGammaMotherYOrBin[i], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendMothersY->AddEntry(   histoGammaMotherYOrBin[i], Form("%s", motherParticlesLatex[i].Data()), "p");
            histoGammaMotherYOrBin[i]->Draw("same");
        }
    }
    legendMothersY->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(                 0.2, 0.22, 0.032, cent, textMeasurement, "", 42, 0.03);
    if (producePlotsForThesis) PutThisThesisLabel(  0.2, 0.17, 0.032, 0.03, 1.25, 42);
    else PutALICESimulationLabel(                   0.2, 0.17, 0.032, 0.03, 1.25, 42);
    
    canvasMothersY->SaveAs(Form("%s/CocktailMothersY_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendMothersY;
    delete canvasMothersY;

    //***************************** Plot cocktail mothers (phi) *****************************************************
    TCanvas *canvasMothersPhi                                   = new TCanvas("canvasMothersPhi","",1100,1200);
    DrawGammaCanvasSettings(canvasMothersPhi, 0.12, 0.025, 0.01, 0.075);
    canvasMothersPhi->SetLogy();

    TLegend* legendMothersPhi                                   = GetAndSetLegend2(0.2, 0.98-(0.04*nRows), 0.95, 0.98, 40, 6);
    legendMothersPhi->SetBorderSize(0);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, 0, 7.);
    SetHistogramm(dummyHist, "#phi", "#frac{1}{N_{ev}} #frac{d#it{N}}{d#phi}", 1e-4, 20, 0.9, 1.25);
    dummyHist->Draw();
    
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaMotherPhiOrBin[i]) {
            DrawGammaSetMarker(         histoGammaMotherPhiOrBin[i], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendMothersPhi->AddEntry( histoGammaMotherPhiOrBin[i], Form("%s", motherParticlesLatex[i].Data()), "p");
            histoGammaMotherPhiOrBin[i]->Draw("same");
        }
    }
    legendMothersPhi->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(                 0.2, 0.22, 0.032, cent, textMeasurement, "", 42, 0.03);
    if (producePlotsForThesis) PutThisThesisLabel(  0.2, 0.17, 0.032, 0.03, 1.25, 42);
    else PutALICESimulationLabel(                   0.2, 0.17, 0.032, 0.03, 1.25, 42);

    canvasMothersPhi->SaveAs(Form("%s/CocktailMothersPhi_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendMothersPhi;
    delete canvasMothersPhi;

    //***************************** Plot cocktail gammas (pt) *******************************************************
    TCanvas *canvasGammas                                       = new TCanvas("canvasGammas","",1100,1200);
    DrawGammaCanvasSettings(canvasGammas, 0.15, 0.01, 0.01, 0.075);
    canvasGammas->SetLogy();
    canvasGammas->SetLogx();
    
    TLegend* legendGammas                                       = GetAndSetLegend2(0.2, 0.98-(0.05*nRows), 0.95, 0.98, 40, 6);
    legendGammas->SetBorderSize(0);
    legendGammas->SetHeader("#gamma from");
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, ptPlotMin, ptPlotMax);
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", 1e-9, histoGammaSumPtOrBin->GetMaximum()*100, 1.0, 1.5);
    dummyHist->SetLabelOffset(-0.015, "X");
    dummyHist->SetTitleOffset(0.8, "X");
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
    
    PutProcessLabelAndEnergyOnPlot(                 0.22, 0.22, 0.032, cent, textMeasurement, "", 42, 0.03);
    if (producePlotsForThesis) PutThisThesisLabel(  0.22, 0.17, 0.032, 0.03, 1.25, 42);
    else PutALICESimulationLabel(                   0.22, 0.17, 0.032, 0.03, 1.25, 42);
    dummyHist->Draw("same,axis");
    
    canvasGammas->SaveAs(Form("%s/CocktailGammas_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendGammas;
    delete canvasGammas;

    //***************************** Plot cocktail gammas to pi0 ratio ***********************************************
    TCanvas *canvasGammasRatio                                  = new TCanvas("canvasGammasRatio","",1100,1200);
    DrawGammaCanvasSettings(canvasGammasRatio, 0.12, 0.01, 0.01, 0.075);
    canvasGammasRatio->SetLogy();
    canvasGammasRatio->SetLogx();
    
    TLegend* legendGammasRatio                                  = GetAndSetLegend2(0.2, 0.98-(0.05*nRows), 0.95, 0.98, 40, 6);
    legendGammasRatio->SetBorderSize(0);
    legendGammasRatio->SetHeader("#gamma from");
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, ptPlotMin, ptPlotMax);
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#gamma_{decay} / #pi^{0}", 1e-6, 5e2, 1.0, 1.3);
    dummyHist->SetLabelOffset(-0.015, "X");
    dummyHist->SetTitleOffset(0.8, "X");
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
            if(i!=0 && i%2==0) tempRatioGammas->SetLineStyle(9);
            tempRatioGammas->Draw("csamehist");
        }
    }
    legendGammasRatio->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(                 0.18, 0.22, 0.032, cent, textMeasurement, "", 42, 0.03);
    if (producePlotsForThesis) PutThisThesisLabel(  0.18, 0.17, 0.032, 0.03, 1.25, 42);
    else PutALICESimulationLabel(                   0.18, 0.17, 0.032, 0.03, 1.25, 42);
    dummyHist->Draw("same,axis");
    
    canvasGammasRatio->SaveAs(Form("%s/CocktailGammasRatioToPi0_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete tempRatioGammas;
    delete legendGammasRatio;
    delete canvasGammasRatio;

    //***************************** Plot cocktail gammas to all gammas ratio ****************************************
    TCanvas *canvasGammasRatio2                                 = new TCanvas("canvasGammasRatio2","",1100,1200);
    DrawGammaCanvasSettings(canvasGammasRatio2, 0.12, 0.01, 0.01, 0.075);
    canvasGammasRatio2->SetLogy();
    canvasGammasRatio2->SetLogx();
    
    TLegend* legendGammasRatio2                                 = GetAndSetLegend2(0.2, 0.98-(0.05*nRows), 0.95, 0.98, 40, 6);
    legendGammasRatio2->SetBorderSize(0);
    legendGammasRatio2->SetHeader("#gamma from");
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, ptPlotMin, ptPlotMax);
    SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#gamma_{source} / #gamma_{decay}", 1e-6, 5e2, 1.0, 1.3);
    dummyHist->SetLabelOffset(-0.015, "X");
    dummyHist->SetTitleOffset(0.8, "X");
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
            if(i!=0 && i%2==0) tempRatioGammas2->SetLineStyle(9);
            tempRatioGammas2->Draw("csamehist");
        }
    }
    legendGammasRatio2->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(                 0.18, 0.22, 0.032, cent, textMeasurement, "", 42, 0.03);
    if (producePlotsForThesis) PutThisThesisLabel(  0.18, 0.17, 0.032, 0.03, 1.25, 42);
    else PutALICESimulationLabel(                   0.18, 0.17, 0.032, 0.03, 1.25, 42);
    dummyHist->Draw("same,axis");
    
    canvasGammasRatio2->SaveAs(Form("%s/CocktailGammasRatioToAll_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete tempRatioGammas2;
    delete legendGammasRatio2;
    delete canvasGammasRatio2;

    //***************************** Plot cocktail gammas (y) ********************************************************
    TCanvas *canvasGammasY                                      = new TCanvas("canvasGammasY","",1100,1200);
    DrawGammaCanvasSettings(canvasGammasY, 0.12, 0.025, 0.01, 0.075);
    canvasGammasY->SetLogy();
    
    TLegend* legendGammasY                                      = GetAndSetLegend2(0.2, 0.98-(0.05*nRows), 0.95, 0.98, 40, 6);
    legendGammasY->SetBorderSize(0);
    legendGammasY->SetHeader("#gamma from");
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, -fRapidity, fRapidity);
    SetHistogramm(dummyHist, "y", "#frac{1}{N_{ev}} #frac{d#it{N}}{dy}", 5e-7, 50, 0.9, 1.25);
    dummyHist->Draw();
    
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaYOrBin[i]) {
            DrawGammaSetMarker(         histoGammaYOrBin[i], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendGammasY->AddEntry(    histoGammaYOrBin[i], motherParticlesLatex[i].Data(), "p");
            histoGammaYOrBin[i]->Draw("same");
        }
    }
    legendGammasY->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(                 0.2, 0.22, 0.032, cent, textMeasurement, "", 42, 0.03);
    if (producePlotsForThesis) PutThisThesisLabel(  0.2, 0.17, 0.032, 0.03, 1.25, 42);
    else PutALICESimulationLabel(                   0.2, 0.17, 0.032, 0.03, 1.25, 42);

    canvasGammasY->SaveAs(Form("%s/CocktailGammasY_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendGammasY;
    delete canvasGammasY;

    //***************************** Plot cocktail mothers (phi) *****************************************************
    TCanvas *canvasGammasPhi                                    = new TCanvas("canvasGammasPhi","",1100,1200);
    DrawGammaCanvasSettings(canvasGammasPhi, 0.12, 0.025, 0.01, 0.075);
    canvasGammasPhi->SetLogy();
    
    TLegend* legendGammasPhi                                    = GetAndSetLegend2(0.2, 0.98-(0.05*nRows), 0.95, 0.98, 40, 6);
    legendGammasPhi->SetBorderSize(0);
    legendGammasPhi->SetHeader("#gamma from");
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, 0, 7.);
    SetHistogramm(dummyHist, "#phi", "#frac{1}{N_{ev}} #frac{d#it{N}}{d#phi}", 5e-7, 50, 0.9, 1.25);
    dummyHist->Draw();
    
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaMotherPhiOrBin[i]) {
            DrawGammaSetMarker(         histoGammaPhiOrBin[i], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendGammasPhi->AddEntry(  histoGammaPhiOrBin[i], motherParticlesLatex[i].Data(), "p");
            histoGammaPhiOrBin[i]->Draw("same");
        }
    }
    legendGammasPhi->Draw("same");
    
    PutProcessLabelAndEnergyOnPlot(                 0.2, 0.22, 0.032, cent, textMeasurement, "", 42, 0.03);
    if (producePlotsForThesis) PutThisThesisLabel(  0.2, 0.17, 0.032, 0.03, 1.25, 42);
    else PutALICESimulationLabel(                   0.2, 0.17, 0.032, 0.03, 1.25, 42);

    canvasGammasPhi->SaveAs(Form("%s/CocktailGammasPhi_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendGammasPhi;
    delete canvasGammasPhi;

    TH1D*   dummyHistRatio                                          = NULL;
    //***************************** Plot mT scaling cross check *****************************************************
    TCanvas* canvasMtCrossCheck                                 = NULL;
    TLegend* legendMtCrossCheck                                 = NULL;
    TPad* padMtCrossCheck                                       = NULL;
    TPad* padMtCrossCheckRatio                                  = NULL;
    dummyHist                                                   = NULL;
    TH1D* tempRatio1                                            = NULL;
    TH1D* tempRatio2                                            = NULL;
    for (Int_t particle=0; particle<nMotherParticles; particle++) {
        if (histoGammaMotherPtOrBin[particle] && cocktailInputParametrizations[particle] && (particle==0 || cocktailInputParametrizationsMtScaled[particle])) {
            canvasMtCrossCheck                                  = new TCanvas("canvasMtCrossCheck","",1100,1200);
            padMtCrossCheck                                     = new TPad("padMtCrossCheck", "", 0., 0.25, 1., 1.,-1, -1, -2);
            padMtCrossCheckRatio                                = new TPad("padMtCrossCheckRatio", "", 0., 0., 1., 0.25,-1, -1, -2);
            legendMtCrossCheck                                  = GetAndSetLegend2(0.55, 0.87-(0.048*3), 0.9, 0.87, 40);
            legendMtCrossCheck->SetBorderSize(0);
            
            DrawGammaCanvasSettings(canvasMtCrossCheck, 0.165, 0.015, 0.025, 0.25);
            DrawGammaPadSettings(padMtCrossCheck,       0.165, 0.015, 0.025, 0.);
            DrawGammaPadSettings(padMtCrossCheckRatio,  0.165, 0.015, 0.0, 0.25);
            
            padMtCrossCheck->Draw();
            padMtCrossCheck->SetLogy();
            padMtCrossCheck->SetLogx();
            
            padMtCrossCheckRatio->Draw();
            padMtCrossCheckRatio->SetLogx();

            // dummy hist
            dummyHist                                           = new TH1D("dummyHist", "", 1000, ptPlotMin, ptPlotMax);
            dummyHistRatio                                      = new TH1D("dummyHist", "", 1000, ptPlotMin, ptPlotMax);
            
            // spectrum + parametrizations
            padMtCrossCheck->cd();
            SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", 5e-8, 1e2, 1.0, 1.8);
            dummyHist->Draw();

            legendMtCrossCheck->AddEntry(histoGammaMotherPtOrBin[particle], Form("%s", motherParticlesLatex[particle].Data()), "l");
            legendMtCrossCheck->AddEntry(cocktailInputParametrizations[particle], Form("%s param.", motherParticlesLatex[particle].Data()), "l");
            if(particle!=0){
                cocktailInputParametrizationsMtScaled[particle]->SetLineColor(kBlack);
                cocktailInputParametrizationsMtScaled[particle]->SetLineStyle(4);
                cocktailInputParametrizationsMtScaled[particle]->SetLineWidth(2);
                legendMtCrossCheck->AddEntry(cocktailInputParametrizationsMtScaled[particle], Form("%s m_{T} scaled param.", motherParticlesLatex[particle].Data()), "l");
                cocktailInputParametrizationsMtScaled[particle]->Draw("same");
            }
            histoGammaMotherPtOrBin[particle]->Draw("e1same");
            cocktailInputParametrizations[particle]->Draw("same");
            legendMtCrossCheck->Draw("same");
            
            PutProcessLabelAndEnergyOnPlot(                 0.22, 0.30, 0.032, cent, textMeasurement, "", 42, 0.03);
            if (producePlotsForThesis) PutThisThesisLabel(  0.22, 0.25, 0.032, 0.03, 1.25, 42);
            else PutALICESimulationLabel(                   0.22, 0.25, 0.032, 0.03, 1.25, 42);

            // ratios of parametrizations to spectrum
            padMtCrossCheckRatio->cd();
            SetStyleHistoTH1ForGraphs(dummyHistRatio, "#it{p}_{T} (GeV/#it{c})","#frac{spec}{param}", 0.12, 0.1, 0.12, 0.1, 1.1, 0.6, 510, 505);
            dummyHistRatio->GetXaxis()->SetLabelOffset(-0.025);
            dummyHistRatio->GetYaxis()->SetRangeUser(0,2.3);
            
            tempRatio1                                          = (TH1D*)CalculateRatioToTF1((TH1D*)histoGammaMotherPtOrBin[particle], cocktailInputParametrizations[particle]);
            tempRatio1->SetLineColor(cocktailColor[particle]);
            if(particle!=0){
              tempRatio2                                          = (TH1D*)CalculateRatioToTF1((TH1D*)histoGammaMotherPtOrBin[particle], cocktailInputParametrizationsMtScaled[particle]);
              tempRatio2->SetLineColor(kBlack);
              tempRatio2->SetMarkerColor(kBlack);
            }
            dummyHistRatio->Draw();
            tempRatio1->Draw("e1same");
            if(particle!=0)tempRatio2->Draw("e1same");
            
            canvasMtCrossCheck->SaveAs(Form("%s/MtScaling%s_%.2f_%s.%s",outputDir.Data(), motherParticles[particle].Data(),fRapidity,cutSelection.Data(),suffix.Data()));
        }
    }

    //***************************** Plot pi0 from data vs. cocktail *************************************************
    if (histoPi0YieldData) {
        // get proper pt range for plotting
        for (Int_t i=1; i<histoPi0YieldData->GetNbinsX()+1; i++) {
            if (histoPi0YieldData->GetBinContent(i)) {
                ptPlotMin                                           = histoPi0YieldData->GetXaxis()->GetBinLowEdge(i);
                break;
            }
        }
        ptPlotMin                                                   = ptPlotMin/2;
        if(ptPlotMin == 0 || ptPlotMin < 1e-3) ptPlotMin            = 1e-3;
        ptPlotMax                                                   = ptMax*2;
        
        
        TCanvas *canvasPi0                                          = new TCanvas("canvasPi0","",1100,1200);
        DrawGammaCanvasSettings(canvasPi0, 0.165, 0.015, 0.025, 0.25);
        canvasPi0->SetLogy();
        canvasPi0->SetLogx();

        TPad *padSpectrum = new TPad("padSpectrum", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings(padSpectrum,       0.165, 0.015, 0.025, 0.);
        padSpectrum->SetBorderSize(0);
        padSpectrum->Draw();
        padSpectrum->SetLogy();
        padSpectrum->SetLogx();

        TPad *padRatio = new TPad("padRatio","", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings(padRatio,  0.165, 0.015, 0.0, 0.25);
        padRatio->Draw();
        padRatio->SetLogx();

        padSpectrum->cd();

        TLegend* legendPi0                                          = GetAndSetLegend2(0.7, 0.95-(0.045*2), 0.85, 0.95, 40);
        legendPi0->SetBorderSize(0);
        dummyHist                                                   = new TH1D("dummyHist", "", 1000, ptPlotMin, ptPlotMax);
        SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", histoPi0YieldData->GetMinimum(0)*0.1, histoPi0YieldData->GetMaximum()*2, 1.0, 1.7);
        dummyHist->SetLabelOffset(-0.015, "X");
        dummyHist->SetTitleOffset(0.8, "X");
        dummyHist->Draw();
        
        DrawGammaSetMarker(histoPi0YieldData, 24, 1, kBlack,  kBlack);
        DrawGammaSetMarker(histoGammaMotherPt[0], 20, 1, kBlue,  kBlue);

        legendPi0->AddEntry(histoPi0YieldData,      Form("%s data", motherParticlesLatex[0].Data()), "p");
        legendPi0->AddEntry(histoGammaMotherPt[0],  Form("%s cocktail", motherParticlesLatex[0].Data()), "p");

        histoPi0YieldData->Draw("same");
        histoGammaMotherPt[0]->Draw("same");
        legendPi0->Draw("same");

        PutProcessLabelAndEnergyOnPlot(                 0.22, 0.22, 0.03, cent, textMeasurement, "", 42, 0.03);
        if (producePlotsForThesis) PutThisThesisLabel(  0.22, 0.17, 0.032, 0.03, 1.25, 42);
        else PutALICESimulationLabel(                   0.22, 0.17, 0.032, 0.03, 1.25, 42);

        padRatio->cd();
        dummyHistRatio                                              = new TH1D("dummyHistRatio", "", 1000, ptPlotMin, ptPlotMax);
        SetStyleHistoTH1ForGraphs(dummyHistRatio, "#it{p}_{T} (GeV/#it{c})","#frac{data}{gen.}", 0.12, 0.1, 0.12, 0.1, 1.1, 0.6, 510, 505);
        dummyHistRatio->GetXaxis()->SetLabelOffset(-0.025);
        dummyHistRatio->GetYaxis()->SetRangeUser(0.65,1.55);
        dummyHistRatio->Draw();
        DrawGammaLines(ptPlotMin,ptPlotMax,1,1,0.1, kGray+1, 1);
        DrawGammaLines(ptPlotMin,ptPlotMax,0.9,0.9,0.1, kGray+1, 7);
        DrawGammaLines(ptPlotMin,ptPlotMax,1.1,1.1,0.1, kGray+1, 7);
        DrawGammaLines(ptPlotMin,ptPlotMax,1.2,1.2,0.1, kGray+1, 8);
        
        ratioPi0DataCocktail = (TH1D*)histoPi0YieldData->Clone("ratioPi0DataCocktail");
        ratioPi0DataCocktail->Divide(histoPi0YieldData,histoGammaMotherPt[0],1.,1.,"");
        ratioPi0DataCocktail->SetLineColor(cocktailColor[0]);
        ratioPi0DataCocktail->SetLineColor(kBlack);
        ratioPi0DataCocktail->SetMarkerColor(kBlack);
        ratioPi0DataCocktail->Draw("same");

        canvasPi0->SaveAs(Form("%s/Pi0DataCocktail_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
        delete legendPi0;
        delete canvasPi0;
    }

    //***************************** Plot pi0 vs. gamma pt + vs. pi0 pt **********************************************
    if (histoGammaMotherPtGamma) {
        TCanvas *canvasPi02                                         = new TCanvas("canvasPi02","",1100,1200);
        TPad *padSpectrum2                                          = new TPad("padSpectrum2", "", 0., 0.25, 1., 1.,-1, -1, -2);
        TPad *padRatio2                                             = new TPad("padRatio2","", 0., 0., 1., 0.25,-1, -1, -2);

        DrawGammaCanvasSettings(canvasPi02, 0.165, 0.015, 0.025,    0.25);
        DrawGammaPadSettings(padSpectrum2,  0.165, 0.015, 0.025,    0.);
        DrawGammaPadSettings(padRatio2,     0.165, 0.015, 0.0,      0.25);

        canvasPi02->SetLogy();
        canvasPi02->SetLogx();

        padSpectrum2->SetBorderSize(0);
        padSpectrum2->Draw();
        padSpectrum2->SetLogy();
        padSpectrum2->SetLogx();

        padRatio2->Draw();
        padRatio2->SetLogx();

        padSpectrum2->cd();

        TLegend* legendPi02                                         = GetAndSetLegend2(0.7, 0.95-(0.045*2), 0.85, 0.95, 40);
        legendPi02->SetBorderSize(0);

        dummyHist                                                   = new TH1D("dummyHist", "", 1000, ptPlotMin, ptPlotMax);
        SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", histoGammaMotherPt[0]->GetMinimum(0)*0.1, histoGammaMotherPt[0]->GetMaximum()*2, 1.0, 1.7);
        dummyHist->SetLabelOffset(-0.015, "X");
        dummyHist->SetTitleOffset(0.8, "X");
        dummyHist->Draw();

        DrawGammaSetMarker(histoGammaMotherPt[0],       24, 1, kBlack,  kBlack);
        DrawGammaSetMarker(histoGammaMotherPtGamma[0],  20, 1, kBlue,  kBlue);

        legendPi02->AddEntry(histoGammaMotherPt[0],         Form("%s(#it{p}_{T, #pi^{0}}) cocktail", motherParticlesLatex[0].Data()), "p");
        legendPi02->AddEntry(histoGammaMotherPtGamma[0],    Form("%s(#it{p}_{T, #gamma}) cocktail", motherParticlesLatex[0].Data()), "p");

        histoGammaMotherPt[0]->Draw("same");
        histoGammaMotherPtGamma[0]->Draw("same");
        legendPi02->Draw("same");

        PutProcessLabelAndEnergyOnPlot(                 0.22, 0.22, 0.03, cent, textMeasurement, "", 42, 0.03);
        if (producePlotsForThesis) PutThisThesisLabel(  0.22, 0.17, 0.032, 0.03, 1.25, 42);
        else PutALICESimulationLabel(                   0.22, 0.17, 0.032, 0.03, 1.25, 42);

        padRatio2->cd();
        padRatio2->SetLogy();
        dummyHistRatio                                              = new TH1D("dummyHistRatio", "", 1000, ptPlotMin, ptPlotMax);
        SetStyleHistoTH1ForGraphs(dummyHistRatio, "#it{p}_{T} (GeV/#it{c})","#pi^{}(#it{p}_{T, #pi^{0}}) / #pi^{0}(#it{p}_{T, #gamma})", 0.12, 0.1, 0.12, 0.1, 1.1, 0.6, 510, 505);
        dummyHistRatio->GetXaxis()->SetLabelOffset(-0.025);
        dummyHistRatio->GetYaxis()->SetRangeUser(.9,20);
        dummyHistRatio->Draw();
//        DrawGammaLines(ptPlotMin,ptPlotMax,1,1,0.1, kGray+1, 1);
//        DrawGammaLines(ptPlotMin,ptPlotMax,0.9,0.9,0.1, kGray+1, 7);
//        DrawGammaLines(ptPlotMin,ptPlotMax,1.1,1.1,0.1, kGray+1, 7);
//        DrawGammaLines(ptPlotMin,ptPlotMax,1.2,1.2,0.1, kGray+1, 8);

        TH1D* tempRatio = (TH1D*)histoGammaMotherPt[0]->Clone("tempRatio");
        tempRatio->Divide(histoGammaMotherPt[0],histoGammaMotherPtGamma[0],1.,1.,"");
        tempRatio->SetLineColor(cocktailColor[0]);
        tempRatio->SetLineColor(kBlack);
        tempRatio->SetMarkerColor(kBlack);
        tempRatio->Draw("same");

        canvasPi02->SaveAs(Form("%s/Pi0VersusGammaPt_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
        delete legendPi02;
        delete canvasPi02;
    }

    //***************************** Plot eta from data vs. cocktail *************************************************
    if (histoEtaYieldData && histoGeneratedEtaPt) {

        // get proper pt range for plotting
        for (Int_t i=1; i<histoEtaYieldData->GetNbinsX()+1; i++) {
            if (histoEtaYieldData->GetBinContent(i)) {
                ptPlotMin                                           = histoEtaYieldData->GetXaxis()->GetBinLowEdge(i);
                break;
            }
        }
        ptPlotMin                                                   = ptPlotMin/2;
        if(ptPlotMin == 0 || ptPlotMin < 1e-3) ptPlotMin            = 1e-3;
        ptPlotMax                                                   = ptMax*2;

        TCanvas *canvasEta                                          = new TCanvas("canvasEta","",1100,1200);
        DrawGammaCanvasSettings(canvasEta, 0.165, 0.015, 0.025, 0.25);
        canvasEta->SetLogy();
        canvasEta->SetLogx();

        TPad *padSpectrumEta                                        = new TPad("padSpectrum", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings(padSpectrumEta,       0.165, 0.015, 0.025, 0.);
        padSpectrumEta->SetBorderSize(0);
        padSpectrumEta->Draw();
        padSpectrumEta->SetLogy();
        padSpectrumEta->SetLogx();

        TPad *padRatioEta                                           = new TPad("padRatioEta","", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings(padRatioEta,  0.165, 0.015, 0.0, 0.25);
        padRatioEta->Draw();
        padRatioEta->SetLogx();

        padSpectrumEta->cd();

        TLegend* legendEta                                          = GetAndSetLegend2(0.7, 0.95-(0.045*2), 0.85, 0.95, 40);
        legendEta->SetBorderSize(0);
        dummyHist                                                   = new TH1D("dummyHist", "", 1000, ptPlotMin, ptPlotMax);
        SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{1}{2#pi#it{p}_{T}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", histoEtaYieldData->GetMinimum(0)*0.1, histoEtaYieldData->GetMaximum()*2, 1.0, 1.7);
        dummyHist->SetLabelOffset(-0.015, "X");
        dummyHist->SetTitleOffset(0.8, "X");
        dummyHist->Draw();

        DrawGammaSetMarker(histoEtaYieldData,   24, 1, kBlack,  kBlack);
        DrawGammaSetMarker(histoGeneratedEtaPt, 20, 1, kBlue,  kBlue);

        legendEta->AddEntry(histoEtaYieldData,  Form("%s data", motherParticlesLatex[1].Data()), "p");
        legendEta->AddEntry(histoGeneratedEtaPt,Form("%s cocktail", motherParticlesLatex[1].Data()), "p");

        histoEtaYieldData->Draw("same");
        histoGeneratedEtaPt->Draw("same");
        legendEta->Draw("same");

        PutProcessLabelAndEnergyOnPlot(                 0.22, 0.22, 0.03, cent, textMeasurement, "", 42, 0.03);
        if (producePlotsForThesis) PutThisThesisLabel(  0.22, 0.17, 0.032, 0.03, 1.25, 42);
        else PutALICESimulationLabel(                   0.22, 0.17, 0.032, 0.03, 1.25, 42);

        padRatioEta->cd();
        dummyHistRatio                                              = new TH1D("dummyHistRatio", "", 1000, ptPlotMin, ptPlotMax);
        SetStyleHistoTH1ForGraphs(dummyHistRatio, "#it{p}_{T} (GeV/#it{c})","#frac{data}{gen.}", 0.12, 0.1, 0.12, 0.1, 1.1, 0.6, 510, 505);
        dummyHistRatio->GetXaxis()->SetLabelOffset(-0.025);
        dummyHistRatio->GetYaxis()->SetRangeUser(0.65,1.55);
        dummyHistRatio->Draw();
        DrawGammaLines(ptPlotMin,ptPlotMax,1,1,0.1, kGray+1, 1);
        DrawGammaLines(ptPlotMin,ptPlotMax,0.9,0.9,0.1, kGray+1, 7);
        DrawGammaLines(ptPlotMin,ptPlotMax,1.1,1.1,0.1, kGray+1, 7);
        DrawGammaLines(ptPlotMin,ptPlotMax,1.2,1.2,0.1, kGray+1, 8);

        ratioEtaDataCocktail = (TH1D*)histoEtaYieldData->Clone("ratioEtaDataCocktail");
        ratioEtaDataCocktail->Divide(histoEtaYieldData,histoGeneratedEtaPt,1.,1.,"");
        ratioEtaDataCocktail->SetLineColor(cocktailColor[0]);
        ratioEtaDataCocktail->SetLineColor(kBlack);
        ratioEtaDataCocktail->SetMarkerColor(kBlack);
        ratioEtaDataCocktail->Draw("same");

        canvasEta->SaveAs(Form("%s/EtaDataCocktail_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
        delete legendEta;
        delete canvasEta;
    }
    delete dummyHist;
    delete dummyHistRatio;

    //***************************** Save histograms *****************************************************************
    if (isSane) {
        SaveHistos();
        CreateBRTableLatex();
    } else {
        TString nameOutput                                      = Form("%s/%s/GammaCocktail%s_%.2f_%s.root", fCutSelection.Data(), fEnergyFlag.Data(), fPeriodFlag.Data(), fRapidity, fCutSelection.Data());
        cout << "ERROR: Possible problem with cocktail normalization, cocktail unequal to summed contributions. If " << nameOutput.Data() << " was already created it will be deleated!" << endl;
        gSystem->Exec("rm "+nameOutput);
    }

    //***************************** Delete objects ******************************************************************
    DeleteObjects();
}

//************************** Initialize binning *********************************************************************
void Initialize(TString energy, Int_t numberOfBins){
    
    InitializeBinning("Pi0", numberOfBins, energy, fdirectphoton, fMode, fEventCutSelection, fClusterCutSelection);
    
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

void RebinSpectrum(TH1F* Spectrum, TH1F* SpectrumForBinning, TString NewName){

    if (!SpectrumForBinning) return;

    if(NewName.CompareTo(""))
        NewName                     = Spectrum->GetName();

    Double_t newBinContent          = 0.;
    for (Int_t i=1; i<Spectrum->GetNbinsX()+1; i++) {
        newBinContent               = Spectrum->GetBinContent(i) * Spectrum->GetBinWidth(i);
        Spectrum->SetBinContent(i, newBinContent);
    }

    Int_t       nBins               = SpectrumForBinning->GetNbinsX();
    Double_t*   binsPt              = new Double_t[nBins+1];
    for (Int_t i=0; i<nBins+1; i++) {
        if (i<nBins)    binsPt[i]   = SpectrumForBinning->GetXaxis()->GetBinLowEdge(i+1);
        else            binsPt[i]   = SpectrumForBinning->GetXaxis()->GetBinUpEdge(i);
    }

    TH1D* deltaPt                   = new TH1D("deltaPt","",nBins,binsPt);
    for(Int_t iPt=1;iPt<nBins+1;iPt++){
        deltaPt->SetBinContent(iPt, binsPt[iPt]-binsPt[iPt-1]);
        deltaPt->SetBinError(iPt,   0);
    }

    *Spectrum                       = *((TH1F*)Spectrum->Rebin(nBins,NewName,binsPt));
    Spectrum->Divide(deltaPt);
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
    } else if (particleName.CompareTo("K0s") == 0) {
        particle                            = pdg->GetParticle(310);
        mass                                = particle->Mass();
    } else if (particleName.CompareTo("K0l") == 0) {
        particle                            = pdg->GetParticle(130);
        mass                                = particle->Mass();
    } else if (particleName.CompareTo("Lambda") == 0) {
        particle                            = pdg->GetParticle(3122);
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
    delete[] histoGammaMotherPtGammaPt;
    delete[] histoGammaMotherPtGammaOrBin;
    delete[] histoGammaMotherPtGamma;
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
    histoGammaSumPtOrBin->Write(                            histoGammaSumPtOrBin->GetName(),    TObject::kOverwrite);
    if (histoGammaSumPtOrBin2) histoGammaSumPtOrBin2->Write(histoGammaSumPtOrBin2->GetName(),   TObject::kOverwrite);
    histoGammaSumYOrBin->Write(                             histoGammaSumYOrBin->GetName(),     TObject::kOverwrite);
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPtOrBin[i])               histoGammaPtOrBin[i]->Write(            histoGammaPtOrBin[i]->GetName(),            TObject::kOverwrite);
        if (histoGammaPtOrBin2[i])              histoGammaPtOrBin2[i]->Write(           histoGammaPtOrBin[i]->GetName(),            TObject::kOverwrite);
        if (histoGammaYOrBin[i])                histoGammaYOrBin[i]->Write(             histoGammaYOrBin[i]->GetName(),             TObject::kOverwrite);
        if (histoGammaPhiOrBin[i])              histoGammaPhiOrBin[i]->Write(           histoGammaPhiOrBin[i]->GetName(),           TObject::kOverwrite);
        if (histoGammaMotherPtOrBin[i])         histoGammaMotherPtOrBin[i]->Write(      histoGammaMotherPtOrBin[i]->GetName(),      TObject::kOverwrite);
        if (histoGammaMotherYOrBin[i])          histoGammaMotherYOrBin[i]->Write(       histoGammaMotherYOrBin[i]->GetName(),       TObject::kOverwrite);
        if (histoGammaMotherPhiOrBin[i])        histoGammaMotherPhiOrBin[i]->Write(     histoGammaMotherPhiOrBin[i]->GetName(),     TObject::kOverwrite);
        if (histoGammaMotherPtGammaOrBin[i])    histoGammaMotherPtGammaOrBin[i]->Write( histoGammaMotherPtGammaOrBin[i]->GetName(), TObject::kOverwrite);
    }
    
    // write rebinned histograms
    histoGammaSumPt->Write(histoGammaSumPt->GetName(), TObject::kOverwrite);
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaPt[i])                histoGammaPt[i]->Write(             histoGammaPt[i]->GetName(),             TObject::kOverwrite);
        if (histoGammaMotherPt[i])          histoGammaMotherPt[i]->Write(       histoGammaMotherPt[i]->GetName(),       TObject::kOverwrite);
        if (histoGammaMotherPtGamma[i])     histoGammaMotherPtGamma[i]->Write(  histoGammaMotherPtGamma[i]->GetName(),  TObject::kOverwrite);
    }
    
    // write input parametrizations and mt scaled ones
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (cocktailInputParametrizations[i])           cocktailInputParametrizations[i]->Write(        cocktailInputParametrizations[i]->GetName(),            TObject::kOverwrite);
        if (cocktailInputParametrizationsMtScaled[i])   cocktailInputParametrizationsMtScaled[i]->Write(cocktailInputParametrizationsMtScaled[i]->GetName(),    TObject::kOverwrite);
    }
    
    if (histoPi0YieldData)          histoPi0YieldData->Write(   "Pi0_invYield",                     TObject::kOverwrite);
    if (ratioPi0DataCocktail)       ratioPi0DataCocktail->Write(ratioPi0DataCocktail->GetName(),    TObject::kOverwrite);

    if (histoEtaYieldData)          histoEtaYieldData->Write(   "Eta_invYield",                     TObject::kOverwrite);
    if (histoGeneratedEtaPt)        histoGeneratedEtaPt->Write( histoGeneratedEtaPt->GetName(),     TObject::kOverwrite);
    if (ratioEtaDataCocktail)       ratioEtaDataCocktail->Write(ratioEtaDataCocktail->GetName(),    TObject::kOverwrite);
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

