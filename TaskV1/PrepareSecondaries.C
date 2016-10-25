// provided by Gamma Conversion Group, $ALICE_PHYSICS/PWGGA/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion
// ***************************************************************************************************************
// **   Friederike Bock,    friederike.bock@cern.ch                                                             **
// **   Nicolas Schmidt,    nicolas.schmidt@cern.ch                                                             **
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
#include "PrepareSecondaries.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
//#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/ExtractSignalBinning.h"

void PrepareSecondaries(    TString     meson                       = "",
                            TString     nameFileCocktail            = "",
                            TString     suffix                      = "eps",
                            TString     cutSelection                = "",
                            TString     option                      = "",
                            TString     directphotonPlots           = "",
                            Double_t    rapidity                    = 0.85,
                            TString     period                      = "",
                            Int_t       numberOfBins                = 30,
                            Int_t       mode                        = 0,
                            Bool_t      producePlotsInOrPtRange     = kFALSE
                     ) {
    
    gROOT->Reset();
    
    //************************** Set general style settings *********************************************************
    StyleSettingsThesis();
    SetPlotStyle();
    
    //************************** Set output directory ***************************************************************
    TString outputDir                                           = Form("%s/%s/%s/PrepareSecondaries",cutSelection.Data(),option.Data(),suffix.Data());
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
    
    //************************** direct photon ana ******************************************************************
    Bool_t doSecondaryGamma                                     = kFALSE;
    if (directphotonPlots.CompareTo("directPhoton")==0) {
        cout << "secondary photon histos will be produced" << endl;
        doSecondaryGamma                                        = kTRUE;
    }
    
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
    fAnalyzedMeson                                              = "";
    if(meson.CompareTo("Pi0") == 0){
        fAnalyzedMeson                                          = "Pi0";
        Initialize("Pi0", fEnergyFlag, numberOfBins);
    } else if (meson.CompareTo("Eta") == 0) {
        fAnalyzedMeson                                          = "Eta";
        Initialize("Eta", fEnergyFlag, numberOfBins);
        cout << "ERROR: Eta not yet fully implemented, returning" << endl;
        return;
    } else if(meson.CompareTo("Pi0EtaBinning") == 0) {
        fAnalyzedMeson                                          = "Pi0";
        Initialize("Pi0EtaBinning", fEnergyFlag, numberOfBins);
    } else   {
        cout << "ERROR: Meson not specified correctly, returning" << endl;
        return;
    }

    //***************************** Cocktail file *******************************************************************
    TFile fileCocktail(nameFileCocktail.Data());
    TDirectoryFile* topDirCocktail                              = (TDirectoryFile*)fileCocktail.Get("HadronicCocktailMC");
    if (!topDirCocktail) {
        cout << "ERROR: TopDirCocktail not found!" << endl;
        return;
    }
    TList* histoListCocktail                                    = NULL;
    if (fAnalyzedMeson.CompareTo("Pi0") == 0) {
        cout << "searching for " << Form("HadronicCocktailMC_pi0_%.2f", rapidity) << endl;
        histoListCocktail                                       = (TList*)topDirCocktail->Get(Form("HadronicCocktailMC_pi0_%.2f", rapidity));
        if (!histoListCocktail) {
            cout << "ERROR: Folder with rapidity " << rapidity << " not contained in cocktail file!" << endl;
            return;
        }
    } else if (fAnalyzedMeson.CompareTo("Eta") == 0) {
        cout << "searching for " << Form("HadronicCocktailMC_eta_%.2f", rapidity) << endl;
        histoListCocktail                                       = (TList*)topDirCocktail->Get(Form("HadronicCocktailMC_eta_%.2f", rapidity));
        if (!histoListCocktail) {
            cout << "ERROR: Folder with rapidity " << rapidity << " not contained in cocktail file!" << endl;
            return;
        }
    } else {
        cout << "ERROR: Meson not specified correctly, returning" << endl;
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
                    TH2F* tempHist                              = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_Pi0_From_%s", motherParticles[j].Data()));
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
    Double_t deltaEta                                           = 2*0.9;
    Double_t deltaPtGen                                         = ptGenMax-ptGenMin;
    Double_t deltaPhi                                           = 2*TMath::Pi();
    cout << "================================"  << endl;
    cout << "deltaRap   "   << deltaRap         << endl;
    cout << "deltaEta   "   << deltaEta         << endl;
    cout << "deltaPtGen "   << deltaPtGen       << endl;
    cout << "deltaPt    "   << ptMax - ptMin    << endl;
    cout << "deltaPhi   "   << deltaPhi         << endl;
    cout << "================================"  << endl;

    //***************************** Get number of spectra ***********************************************************
    Int_t nSpectra                                              = 0;
    for (Int_t i=0; i<nMotherParticles; i++)
        if (hasMother[i]) nSpectra++;
    Int_t nRows                                                 = 0;
    if (nSpectra%6 == 0) nRows                                  = nSpectra/6 + 1;
    else nRows                                                  = (nSpectra+1)/6 + 1;

    //***************************** Get number of events (cocktail) *************************************************
    histoNEvents                                                = (TH1F*)histoListCocktail->FindObject("NEvents");
    nEvents                                                     = histoNEvents->GetEntries();
    cout << nEvents << " events" << endl;

    //***************************** Read histograms from cocktail file **********************************************
    // mesons
    histoDecayChannels                                          = new TH1F*[nMotherParticles];
    histoMesonDaughterPtY                                       = new TH2F*[nMotherParticles];
    histoMesonDaughterPtPhi                                     = new TH2F*[nMotherParticles];
    histoMesonMotherPtY                                         = new TH2F*[nMotherParticles];
    histoMesonMotherPtPhi                                       = new TH2F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (hasMother[i]) {
            histoDecayChannels[i]                               = (TH1F*)histoListCocktail->FindObject(Form("DecayChannels_%s",motherParticles[i].Data()));
            
            histoMesonDaughterPtY[i]                            = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_%s_From_%s",fAnalyzedMeson.Data(),motherParticles[i].Data()));
            histoMesonDaughterPtY[i]->SetName(Form("%s_From_%s_Pt_Y_OrBin",fAnalyzedMeson.Data(),motherParticles[i].Data()));
            histoMesonDaughterPtY[i]->Sumw2();
            
            histoMesonDaughterPtPhi[i]                          = (TH2F*)histoListCocktail->FindObject(Form("Pt_Phi_%s_From_%s",fAnalyzedMeson.Data(),motherParticles[i].Data()));
            histoMesonDaughterPtPhi[i]->SetName(Form("%s_From_%s_Pt_Phi_OrBin",fAnalyzedMeson.Data(),motherParticles[i].Data()));
            histoMesonDaughterPtPhi[i]->Sumw2();
            
            histoMesonMotherPtY[i]                              = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_%s",motherParticles[i].Data()));
            histoMesonMotherPtY[i]->SetName(Form("%s_Pt_Y_OrBin", motherParticles[i].Data()));
            histoMesonMotherPtY[i]->Sumw2();
            
            histoMesonMotherPtPhi[i]                            = (TH2F*)histoListCocktail->FindObject(Form("Pt_Phi_%s",motherParticles[i].Data()));
            histoMesonMotherPtPhi[i]->SetName(Form("%s_Pt_Phi_OrBin", motherParticles[i].Data()));
            histoMesonMotherPtPhi[i]->Sumw2();
        } else {
            histoDecayChannels[i]                               = NULL;
            histoMesonDaughterPtY[i]                            = NULL;
            histoMesonDaughterPtPhi[i]                          = NULL;
            histoMesonMotherPtY[i]                              = NULL;
            histoMesonMotherPtPhi[i]                            = NULL;
        }
    }
    
    // gammas
    histoGammaFromXFromMotherPtY                                = new TH2F*[nMotherParticles];
    histoGammaFromXFromMotherPtPhi                              = new TH2F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (doSecondaryGamma) {
            if (hasMother[i]) {
                histoGammaFromXFromMotherPtY[i]                 = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_Gamma_From_X_From_%s",motherParticles[i].Data()));
                if (histoGammaFromXFromMotherPtY[i]) {
                    histoGammaFromXFromMotherPtY[i]->SetName(Form("Gamma_From_X_From_%s_Pt_Y_OrBin", motherParticles[i].Data()));
                    histoGammaFromXFromMotherPtY[i]->Sumw2();
                }
                
                histoGammaFromXFromMotherPtPhi[i]               = (TH2F*)histoListCocktail->FindObject(Form("Pt_Phi_Gamma_From_X_From_%s",motherParticles[i].Data()));
                if (histoGammaFromXFromMotherPtPhi[i]) {
                    histoGammaFromXFromMotherPtPhi[i]->SetName(Form("Gamma_From_X_From_%s_Pt_Phi_OrBin", motherParticles[i].Data()));
                    histoGammaFromXFromMotherPtPhi[i]->Sumw2();
                }
            } else {
                histoGammaFromXFromMotherPtY[i]                 = NULL;
                histoGammaFromXFromMotherPtPhi[i]               = NULL;
            }
        } else {
            histoGammaFromXFromMotherPtY[i]                     = NULL;
            histoGammaFromXFromMotherPtPhi[i]                   = NULL;
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
                    if (tempBinLabel.Contains("#pi^{0}"))
                        decayChannelsLatex[i][bin-2]            = tempBinLabel;
                }
            }
        }
    }

    //***************************** Project from 2D histograms ******************************************************
    // mesons
    histoMesonDaughterPtOrBin                                   = new TH1F*[nMotherParticles];
    histoMesonDaughterYOrBin                                    = new TH1F*[nMotherParticles];
    histoMesonDaughterPhiOrBin                                  = new TH1F*[nMotherParticles];
    histoMesonMotherPtOrBin                                     = new TH1F*[nMotherParticles];
    histoMesonMotherYOrBin                                      = new TH1F*[nMotherParticles];
    histoMesonMotherPhiOrBin                                    = new TH1F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoMesonDaughterPtY[i]) {
            histoMesonDaughterPtOrBin[i]                        = (TH1F*)histoMesonDaughterPtY[i]->ProjectionX(Form("%s_From_%s_Pt_OrBin",fAnalyzedMeson.Data(),motherParticles[i].Data()),1,histoMesonDaughterPtY[i]->GetNbinsY(),"e");
            SetHistogramTitles(histoMesonDaughterPtOrBin[i],"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoMesonDaughterPtOrBin[i]->Sumw2();
            histoMesonDaughterPtOrBin[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
            histoMesonDaughterYOrBin[i]                         = (TH1F*)histoMesonDaughterPtY[i]->ProjectionY(Form("%s_From_%s_Y_OrBin",fAnalyzedMeson.Data(),motherParticles[i].Data()),1,histoMesonDaughterPtY[i]->GetNbinsX(),"e");
            SetHistogramTitles(histoMesonDaughterYOrBin[i],"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoMesonDaughterYOrBin[i]->Sumw2();
        } else {
            histoMesonDaughterPtOrBin[i]                        = NULL;
            histoMesonDaughterYOrBin[i]                         = NULL;
        }
        if (histoMesonDaughterPtPhi[i]) {
            histoMesonDaughterPhiOrBin[i]                       = (TH1F*)histoMesonDaughterPtPhi[i]->ProjectionY(Form("%s_From_%s_Phi_OrBin",fAnalyzedMeson.Data(),motherParticles[i].Data()),1,histoMesonDaughterPtPhi[i]->GetNbinsX(),"e");
            SetHistogramTitles(histoMesonDaughterPhiOrBin[i],"","#phi","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoMesonDaughterPhiOrBin[i]->Sumw2();
        } else
            histoMesonDaughterPhiOrBin[i]                       = NULL;
        
        if (histoMesonMotherPtY[i]) {
            histoMesonMotherPtOrBin[i]                          = (TH1F*)histoMesonMotherPtY[i]->ProjectionX(Form("%s_Pt_OrBin",motherParticles[i].Data()),1,histoMesonMotherPtY[i]->GetNbinsY(),"e");
            SetHistogramTitles(histoMesonMotherPtOrBin[i],"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoMesonMotherPtOrBin[i]->Sumw2();
            histoMesonMotherPtOrBin[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
            histoMesonMotherPtOrBin[i]->Scale(deltaPhi);
            histoMesonMotherYOrBin[i]                           = (TH1F*)histoMesonMotherPtY[i]->ProjectionY(Form("%s_Y_OrBin",motherParticles[i].Data()),1,histoMesonMotherPtY[i]->GetNbinsX(),"e");
            SetHistogramTitles(histoMesonMotherYOrBin[i],"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoMesonMotherYOrBin[i]->Sumw2();
        } else {
            histoMesonMotherPtOrBin[i]                          = NULL;
            histoMesonMotherYOrBin[i]                           = NULL;
        }
        if (histoMesonMotherPtPhi[i]) {
            histoMesonMotherPhiOrBin[i]                         = (TH1F*)histoMesonMotherPtPhi[i]->ProjectionY(Form("%s_Phi_OrBin",motherParticles[i].Data()),1,histoMesonMotherPtPhi[i]->GetNbinsX(),"e");
            SetHistogramTitles(histoMesonMotherPhiOrBin[i],"","#phi","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoMesonMotherPhiOrBin[i]->Sumw2();
        } else
            histoMesonMotherPhiOrBin[i]                         = NULL;
    }

    if (!histoMesonMotherPtOrBin[1]) {
        cout << "ERROR: Didn't get K0s pt spectrum, returning!" << endl;
        return;
    }
    
    // gammas
    histoGammaFromXFromMotherPtOrBin                            = new TH1F*[nMotherParticles];
    histoGammaFromXFromMotherYOrBin                             = new TH1F*[nMotherParticles];
    histoGammaFromXFromMotherPhiOrBin                           = new TH1F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        if(doSecondaryGamma) {
            if (histoGammaFromXFromMotherPtY[i]) {
                histoGammaFromXFromMotherPtOrBin[i]             = (TH1F*)histoGammaFromXFromMotherPtY[i]->ProjectionX(Form("Gamma_From_X_From_%s_Pt_OrBin",motherParticles[i].Data()),1,histoGammaFromXFromMotherPtY[i]->GetNbinsY(),"e");
                SetHistogramTitles(histoGammaFromXFromMotherPtOrBin[i],"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
                histoGammaFromXFromMotherPtOrBin[i]->Sumw2();
                histoGammaFromXFromMotherPtOrBin[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
                histoGammaFromXFromMotherYOrBin[i]              = (TH1F*)histoGammaFromXFromMotherPtY[i]->ProjectionY(Form("Gamma_From_X_From_%s_Y_OrBin",motherParticles[i].Data()),1,histoGammaFromXFromMotherPtY[i]->GetNbinsX(),"e");
                SetHistogramTitles(histoGammaFromXFromMotherYOrBin[i],"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
                histoGammaFromXFromMotherYOrBin[i]->Sumw2();
            } else {
                histoGammaFromXFromMotherPtOrBin[i]             = NULL;
                histoGammaFromXFromMotherYOrBin[i]              = NULL;
            }
            
            if (histoGammaFromXFromMotherPtPhi[i]) {
                histoGammaFromXFromMotherPhiOrBin[i]            = (TH1F*)histoGammaFromXFromMotherPtPhi[i]->ProjectionY(Form("Gamma_From_X_From_%s_Phi_OrBin",motherParticles[i].Data()),1,histoGammaFromXFromMotherPtPhi[i]->GetNbinsX(),"e");
                SetHistogramTitles(histoGammaFromXFromMotherPhiOrBin[i],"","#phi","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
                histoGammaFromXFromMotherPhiOrBin[i]->Sumw2();
            } else {
                histoGammaFromXFromMotherPhiOrBin[i]            = NULL;
            }
        } else {
            histoGammaFromXFromMotherPtOrBin[i]                 = NULL;
            histoGammaFromXFromMotherYOrBin[i]                  = NULL;
            histoGammaFromXFromMotherPhiOrBin[i]                = NULL;
        }
    }
    
    //***************************** Scale spectra *******************************************************************
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoMesonDaughterPtOrBin[i])           histoMesonDaughterPtOrBin[i]->Scale(        1./nEvents);
        if (histoMesonDaughterYOrBin[i])            histoMesonDaughterYOrBin[i]->Scale(         1./nEvents);
        if (histoMesonDaughterPhiOrBin[i])          histoMesonDaughterPhiOrBin[i]->Scale(       1./nEvents);
        if (histoMesonMotherPtOrBin[i])             histoMesonMotherPtOrBin[i]->Scale(          1./nEvents);
        if (histoMesonMotherYOrBin[i])              histoMesonMotherYOrBin[i]->Scale(           1./nEvents);
        if (histoMesonMotherPhiOrBin[i])            histoMesonMotherPhiOrBin[i]->Scale(         1./nEvents);
        if (histoGammaFromXFromMotherPtOrBin[i])    histoGammaFromXFromMotherPtOrBin[i]->Scale( 1./nEvents);
        if (histoGammaFromXFromMotherYOrBin[i])     histoGammaFromXFromMotherYOrBin[i]->Scale(  1./nEvents);
        if (histoGammaFromXFromMotherPhiOrBin[i])   histoGammaFromXFromMotherPhiOrBin[i]->Scale(1./nEvents);
    }

    
    //
    // ******************************************
    // * PPPPP   LL       OOOO   TTTTTT   SSSSS *
    // * PP  PP  LL      OO  OO    TT    SS     *
    // * PP  PP  LL      OO  OO    TT     SSSS  *
    // * PPPPP   LL      OO  OO    TT        SS *
    // * PP      LL      OO  OO    TT        SS *
    // * PP      LLLLLL   OOOO     TT    SSSSS  *
    // ******************************************
    //
    
    // adapt plotting range for original binned histograms
    ///if (!producePlotsInOrPtRange) {
    //    ptGenMin = ptMin;
    //    ptGenMax = ptMax;
    //}
    //
    // plotting
    //
    
    
    //***************************** Save histograms *****************************************************************
    SaveMesonHistos();
    if (doSecondaryGamma) SavePhotonHistos();
    CreateBRTableLatex();
}

//************************** Initialize binning *********************************************************************
void Initialize(TString meson, TString energy, Int_t numberOfBins){
    
    InitializeBinning(meson, numberOfBins, energy, fdirectphoton, fMode, fEventCutSelection, fClusterCutSelection);
    
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


//************************** Routine for saving histograms **********************************************************
void SaveMesonHistos() {
    TString nameOutput                  = Form("%s/%s/Secondary%s%s_%.2f_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fAnalyzedMeson.Data(),fPeriodFlag.Data(),fRapidity,fCutSelection.Data());
    cout << "INFO: writing into: " << nameOutput << endl;
    TFile *outputFile                   = new TFile(nameOutput,"UPDATE");
    
    // write number of events
    histoNEvents->Write("NEvents", TObject::kOverwrite);
    
    // write projections
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoMesonDaughterPtOrBin[i])   histoMesonDaughterPtOrBin[i]->Write(    histoMesonDaughterPtOrBin[i]->GetName(),    TObject::kOverwrite);
        if (histoMesonDaughterYOrBin[i])    histoMesonDaughterYOrBin[i]->Write(     histoMesonDaughterYOrBin[i]->GetName(),     TObject::kOverwrite);
        if (histoMesonDaughterPhiOrBin[i])  histoMesonDaughterPhiOrBin[i]->Write(   histoMesonDaughterPhiOrBin[i]->GetName(),   TObject::kOverwrite);
        if (histoMesonMotherPtOrBin[i])     histoMesonMotherPtOrBin[i]->Write(      histoMesonMotherPtOrBin[i]->GetName(),      TObject::kOverwrite);
        if (histoMesonMotherYOrBin[i])      histoMesonMotherYOrBin[i]->Write(       histoMesonMotherYOrBin[i]->GetName(),       TObject::kOverwrite);
        if (histoMesonMotherPhiOrBin[i])    histoMesonMotherPhiOrBin[i]->Write(     histoMesonMotherPhiOrBin[i]->GetName(),     TObject::kOverwrite);
    }
}

//************************** Routine for saving histograms **********************************************************
void SavePhotonHistos() {
    TString nameOutput                  = Form("%s/%s/SecondaryGamma%s_%.2f_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fPeriodFlag.Data(),fRapidity,fCutSelection.Data());
    cout << "INFO: writing into: " << nameOutput << endl;
    TFile *outputFile                   = new TFile(nameOutput,"UPDATE");
    
    // write number of events
    histoNEvents->Write("NEvents", TObject::kOverwrite);
    
    // write projections
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaFromXFromMotherPtOrBin[i])    histoGammaFromXFromMotherPtOrBin[i]->Write( histoGammaFromXFromMotherPtOrBin[i]->GetName(),     TObject::kOverwrite);
        if (histoGammaFromXFromMotherYOrBin[i])     histoGammaFromXFromMotherYOrBin[i]->Write(  histoGammaFromXFromMotherYOrBin[i]->GetName(),      TObject::kOverwrite);
        if (histoGammaFromXFromMotherPhiOrBin[i])   histoGammaFromXFromMotherPhiOrBin[i]->Write(histoGammaFromXFromMotherPhiOrBin[i]->GetName(),    TObject::kOverwrite);
    }
}

//************************** Create tex file containing BR table ****************************************************
void CreateBRTableLatex() {
    
    // tex file
    TString texFileName                  = Form("%s/%s/HadronicCocktail%s_%.2f_%s_BranchingRatioTable.tex", fCutSelection.Data(), fEnergyFlag.Data(), fPeriodFlag.Data(), fRapidity, fCutSelection.Data());
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