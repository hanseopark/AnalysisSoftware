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
#include "TList.h"
#include "TFitResultPtr.h"
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
                            Bool_t      producePlotsInOrPtRange     = kFALSE,
                            Bool_t      doRapidityCorrection        = kFALSE
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
    if (option.Contains("PbPb")) cent                           = Form("%s %s", centrality.Data(), collisionSystem.Data());
    else                         cent                           = collisionSystem;
    
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
    
    cocktailInputParametrizationPi0                             = (TF1*)cocktailSettingsList->FindObject("111_pt");
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
        else {
            if (cocktailInputParametrizationPi0)
                cocktailInputParametrizationsMtScaled[i]        = (TF1*)MtScaledParam(cocktailInputParametrizationPi0, i);
            else
                cocktailInputParametrizationsMtScaled[i]        = (TF1*)MtScaledParam(cocktailInputParametrizations[0], i);
        }
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
    
    Double_t scalingEta                                         = 1.;
    Double_t scalingPhi                                         = 1.;
    if (isCalo) {
        scalingEta                                              = deltaEtaCalo/deltaEta;
        scalingPhi                                              = deltaPhiCalo/deltaPhi;
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
    cout << "add. scaling eta = "     << scalingEta     << endl;
    cout << "add. scaling phi = "     << scalingPhi     << endl;
    cout << "========================================"  << endl;

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
    
    //***************************** Correct for non-flat y distributions ********************************************
    histoMesonDaughterPtYCorr                                   = new TH2F*[nMotherParticles];
    histoGammaFromXFromMotherPtYCorr                            = new TH2F*[nMotherParticles];
    listYSlicesMesonDaughter                                    = new TList*[nMotherParticles];
    listYSlicesGammaFromXFromMother                             = new TList*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (doRapidityCorrection) {
            if (histoMesonDaughterPtY[i]){
                histoMesonDaughterPtYCorr[i]                    = (TH2F*)histoMesonDaughterPtY[i]->Clone(Form("%s_yCorr", histoMesonDaughterPtY[i]->GetName()));
                listYSlicesMesonDaughter[i]                     = new TList;
                listYSlicesMesonDaughter[i]->SetName(Form("%s_ySlices", histoMesonDaughterPtY[i]->GetName()));
                CorrectForNonFlatRapidity(histoMesonDaughterPtYCorr[i], histoMesonDaughterPtY[i], listYSlicesMesonDaughter[i]);
            } else {
                histoMesonDaughterPtYCorr[i]                    = NULL;
                listYSlicesMesonDaughter[i]                     = NULL;
            }
            
            if (histoGammaFromXFromMotherPtY[i]) {
                histoGammaFromXFromMotherPtYCorr[i]             = (TH2F*)histoGammaFromXFromMotherPtY[i]->Clone(Form("%s_yCorr", histoGammaFromXFromMotherPtY[i]->GetName()));
                listYSlicesGammaFromXFromMother[i]              = new TList;
                listYSlicesGammaFromXFromMother[i]->SetName(Form("%s_ySlices",  histoGammaFromXFromMotherPtY[i]->GetName()));
                CorrectForNonFlatRapidity(histoGammaFromXFromMotherPtYCorr[i],  histoGammaFromXFromMotherPtY[i], listYSlicesGammaFromXFromMother[i]);
            } else {
                histoGammaFromXFromMotherPtYCorr[i]             = NULL;
                listYSlicesGammaFromXFromMother[i]              = NULL;
            }
        } else {
            if (histoMesonDaughterPtY[i])
                histoMesonDaughterPtYCorr[i]                    = (TH2F*)histoMesonDaughterPtY[i]->Clone(histoMesonDaughterPtY[i]->GetName());
            else
                histoMesonDaughterPtYCorr[i]                    = NULL;
            listYSlicesMesonDaughter[i]                         = NULL;
            if (histoGammaFromXFromMotherPtY[i])
                histoGammaFromXFromMotherPtYCorr[i]             = (TH2F*)histoGammaFromXFromMotherPtY[i]->Clone(histoGammaFromXFromMotherPtY[i]->GetName());
            else
                histoGammaFromXFromMotherPtYCorr[i]             = NULL;
            listYSlicesGammaFromXFromMother[i]                  = NULL;
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
        
        if (histoMesonDaughterPtYCorr[i]) {
            // project pt distributions (from corrected pt-y histograms)
            histoMesonDaughterPtOrBin[i]                        = (TH1F*)histoMesonDaughterPtYCorr[i]->ProjectionX(Form("%s_From_%s_Pt_OrBin",fAnalyzedMeson.Data(),motherParticles[i].Data()),1,histoMesonDaughterPtY[i]->GetNbinsY(),"e");
            SetHistogramTitles(histoMesonDaughterPtOrBin[i],"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoMesonDaughterPtOrBin[i]->Sumw2();
            histoMesonDaughterPtOrBin[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
            
            // project y distributions (from corrected pt-y histograms)
            histoMesonDaughterYOrBin[i]                         = (TH1F*)histoMesonDaughterPtYCorr[i]->ProjectionY(Form("%s_From_%s_Y_OrBin",fAnalyzedMeson.Data(),motherParticles[i].Data()),1,histoMesonDaughterPtY[i]->GetNbinsX(),"e");
            SetHistogramTitles(histoMesonDaughterYOrBin[i],"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoMesonDaughterYOrBin[i]->Sumw2();
        } else {
            histoMesonDaughterPtOrBin[i]                        = NULL;
            histoMesonDaughterYOrBin[i]                         = NULL;
        }
        
        // project phi distributions
        if (histoMesonDaughterPtPhi[i]) {
            histoMesonDaughterPhiOrBin[i]                       = (TH1F*)histoMesonDaughterPtPhi[i]->ProjectionY(Form("%s_From_%s_Phi_OrBin",fAnalyzedMeson.Data(),motherParticles[i].Data()),1,histoMesonDaughterPtPhi[i]->GetNbinsX(),"e");
            SetHistogramTitles(histoMesonDaughterPhiOrBin[i],"","#phi","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoMesonDaughterPhiOrBin[i]->Sumw2();
        } else
            histoMesonDaughterPhiOrBin[i]                       = NULL;
        
        // project pt and y distributions
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
        
        // project phi distributions
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
        if (doSecondaryGamma) {
            
            if (histoGammaFromXFromMotherPtYCorr[i]) {
                // project pt distributions (from corrected pt-y histograms)
                histoGammaFromXFromMotherPtOrBin[i]             = (TH1F*)histoGammaFromXFromMotherPtYCorr[i]->ProjectionX(Form("Gamma_From_X_From_%s_Pt_OrBin",motherParticles[i].Data()),1,histoGammaFromXFromMotherPtY[i]->GetNbinsY(),"e");
                SetHistogramTitles(histoGammaFromXFromMotherPtOrBin[i],"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
                histoGammaFromXFromMotherPtOrBin[i]->Sumw2();
                histoGammaFromXFromMotherPtOrBin[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
                
                // project y distributions (from corrected pt-y histograms)
                histoGammaFromXFromMotherYOrBin[i]              = (TH1F*)histoGammaFromXFromMotherPtYCorr[i]->ProjectionY(Form("Gamma_From_X_From_%s_Y_OrBin",motherParticles[i].Data()),1,histoGammaFromXFromMotherPtY[i]->GetNbinsX(),"e");
                SetHistogramTitles(histoGammaFromXFromMotherYOrBin[i],"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
                histoGammaFromXFromMotherYOrBin[i]->Sumw2();
            } else {
                histoGammaFromXFromMotherPtOrBin[i]             = NULL;
                histoGammaFromXFromMotherYOrBin[i]              = NULL;
            }
            
            // project phi distributions
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
    Double_t factorNEvents                                 = 1./nEvents;
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoMesonDaughterPtOrBin[i])           histoMesonDaughterPtOrBin[i]->Scale(        factorNEvents);
        if (histoMesonDaughterYOrBin[i])            histoMesonDaughterYOrBin[i]->Scale(         factorNEvents);
        if (histoMesonDaughterPhiOrBin[i])          histoMesonDaughterPhiOrBin[i]->Scale(       factorNEvents);
        if (histoMesonMotherPtOrBin[i])             histoMesonMotherPtOrBin[i]->Scale(          factorNEvents);
        if (histoMesonMotherYOrBin[i])              histoMesonMotherYOrBin[i]->Scale(           factorNEvents);
        if (histoMesonMotherPhiOrBin[i])            histoMesonMotherPhiOrBin[i]->Scale(         factorNEvents);
        if (histoGammaFromXFromMotherPtOrBin[i])    histoGammaFromXFromMotherPtOrBin[i]->Scale( factorNEvents*scalingEta*scalingPhi);
        if (histoGammaFromXFromMotherYOrBin[i])     histoGammaFromXFromMotherYOrBin[i]->Scale(  factorNEvents*scalingEta*scalingPhi);
        if (histoGammaFromXFromMotherPhiOrBin[i])   histoGammaFromXFromMotherPhiOrBin[i]->Scale(factorNEvents*scalingEta*scalingPhi);
    }

    //***************************** calculate (pi0 from X)/pi0_param ************************************************
    histoRatioPi0FromXToPi0Param                                = new TH1F*[nMotherParticles];
    if (fAnalyzedMeson.CompareTo("Pi0") == 0) {
        if (cocktailInputParametrizationPi0) {
            for (Int_t i=0; i<nMotherParticles; i++) {
                if(histoMesonDaughterPtOrBin[i]) {
                    histoRatioPi0FromXToPi0Param[i]             = (TH1F*)histoMesonDaughterPtOrBin[i]->Clone(Form("Ratio_Pi0_From_%s_To_Pi0Param",motherParticles[i].Data()));
                    SetHistogramTitles(histoRatioPi0FromXToPi0Param[i],"","#it{p}_{T} (GeV/#it{c})",Form("(#pi^{0} from %s)/#pi^{0} param.",motherParticlesLatex[i].Data()));
                    Double_t binContent                         = 0.;
                    Double_t binError                           = 0.;
                    Double_t binIntegralParam                   = 0.;
                    for(Int_t bin=1; bin<histoRatioPi0FromXToPi0Param[i]->GetNbinsX()+1; bin++) {
                        // get bin content & error
                        binContent                              = histoRatioPi0FromXToPi0Param[i]->GetBinContent(i) * histoRatioPi0FromXToPi0Param[i]->GetBinWidth(i);
                        binError                                = histoRatioPi0FromXToPi0Param[i]->GetBinError(i);
                        binIntegralParam                        = cocktailInputParametrizationPi0->Integral(histoRatioPi0FromXToPi0Param[i]->GetXaxis()->GetBinLowEdge(i),histoRatioPi0FromXToPi0Param[i]->GetXaxis()->GetBinUpEdge(i));
                        
                        // set new bin content & error
                        if(binIntegralParam) {
                            histoRatioPi0FromXToPi0Param[i]->SetBinContent( i, binContent/binIntegralParam);
                            histoRatioPi0FromXToPi0Param[i]->SetBinError(   i, binError/binIntegralParam);
                        } else {
                            histoRatioPi0FromXToPi0Param[i]->SetBinContent( i, 0);
                            histoRatioPi0FromXToPi0Param[i]->SetBinError(   i, 0);
                        }
                    }
                } else {
                    histoRatioPi0FromXToPi0Param[i]             = NULL;
                }
            }
        } else {
            cout << "Pi0 param missing, can't calculate histoRatioPi0FromXToPi0Param." << endl;
            for (Int_t i=0; i<nMotherParticles; i++)
                histoRatioPi0FromXToPi0Param[i]                 = NULL;
        }
    } else {
        for (Int_t i=0; i<nMotherParticles; i++)
            histoRatioPi0FromXToPi0Param[i]                     = NULL;
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
    TFile *outputFile                   = new TFile(nameOutput,"RECREATE");
    
    // write number of events
    histoNEvents->Write("NEvents", TObject::kOverwrite);
    
    // write y vs pt histograms
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoMesonDaughterPtY[i])       histoMesonDaughterPtY[i]->Write(        histoMesonDaughterPtY[i]->GetName(),        TObject::kOverwrite);
        if (histoMesonDaughterPtYCorr[i])   histoMesonDaughterPtYCorr[i]->Write(    histoMesonDaughterPtYCorr[i]->GetName(),    TObject::kOverwrite);
        if (listYSlicesMesonDaughter[i])    listYSlicesMesonDaughter[i]->Write(     listYSlicesMesonDaughter[i]->GetName(),     TObject::kSingleKey);
    }
    
    // write projections
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoMesonDaughterPtOrBin[i])   histoMesonDaughterPtOrBin[i]->Write(    histoMesonDaughterPtOrBin[i]->GetName(),    TObject::kOverwrite);
        if (histoMesonDaughterYOrBin[i])    histoMesonDaughterYOrBin[i]->Write(     histoMesonDaughterYOrBin[i]->GetName(),     TObject::kOverwrite);
        if (histoMesonDaughterPhiOrBin[i])  histoMesonDaughterPhiOrBin[i]->Write(   histoMesonDaughterPhiOrBin[i]->GetName(),   TObject::kOverwrite);
        if (histoMesonMotherPtOrBin[i])     histoMesonMotherPtOrBin[i]->Write(      histoMesonMotherPtOrBin[i]->GetName(),      TObject::kOverwrite);
        if (histoMesonMotherYOrBin[i])      histoMesonMotherYOrBin[i]->Write(       histoMesonMotherYOrBin[i]->GetName(),       TObject::kOverwrite);
        if (histoMesonMotherPhiOrBin[i])    histoMesonMotherPhiOrBin[i]->Write(     histoMesonMotherPhiOrBin[i]->GetName(),     TObject::kOverwrite);
    }

    // write input params
    if (cocktailInputParametrizationPi0)                cocktailInputParametrizationPi0->Write(         cocktailInputParametrizationPi0->GetName(),         TObject::kOverwrite);
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (cocktailInputParametrizations[i])           cocktailInputParametrizations[i]->Write(        cocktailInputParametrizations[i]->GetName(),        TObject::kOverwrite);
        if (cocktailInputParametrizationsMtScaled[i])   cocktailInputParametrizationsMtScaled[i]->Write(cocktailInputParametrizationsMtScaled[i]->GetName(),TObject::kOverwrite);
    }

    // write ratio pi0 from X to pi0 param
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoRatioPi0FromXToPi0Param[i]) histoRatioPi0FromXToPi0Param[i]->Write(histoRatioPi0FromXToPi0Param[i]->GetName(), TObject::kOverwrite);
    }
}

//************************** Routine for saving histograms **********************************************************
void SavePhotonHistos() {
    TString nameOutput                  = Form("%s/%s/SecondaryGamma%s_%.2f_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fPeriodFlag.Data(),fRapidity,fCutSelection.Data());
    cout << "INFO: writing into: " << nameOutput << endl;
    TFile *outputFile                   = new TFile(nameOutput,"RECREATE");
    
    // write number of events
    histoNEvents->Write("NEvents", TObject::kOverwrite);
    
    // write y vs pt histograms
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaFromXFromMotherPtYCorr[i])    histoGammaFromXFromMotherPtYCorr[i]->Write( histoGammaFromXFromMotherPtYCorr[i]->GetName(), TObject::kOverwrite);
        if (histoGammaFromXFromMotherPtY[i])        histoGammaFromXFromMotherPtY[i]->Write(     histoGammaFromXFromMotherPtY[i]->GetName(),     TObject::kOverwrite);
        if (listYSlicesGammaFromXFromMother[i])     listYSlicesGammaFromXFromMother[i]->Write(  listYSlicesGammaFromXFromMother[i]->GetName(),  TObject::kSingleKey);
    }
    
    // write projections
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoGammaFromXFromMotherPtOrBin[i])    histoGammaFromXFromMotherPtOrBin[i]->Write( histoGammaFromXFromMotherPtOrBin[i]->GetName(), TObject::kOverwrite);
        if (histoGammaFromXFromMotherYOrBin[i])     histoGammaFromXFromMotherYOrBin[i]->Write(  histoGammaFromXFromMotherYOrBin[i]->GetName(),  TObject::kOverwrite);
        if (histoGammaFromXFromMotherPhiOrBin[i])   histoGammaFromXFromMotherPhiOrBin[i]->Write(histoGammaFromXFromMotherPhiOrBin[i]->GetName(),TObject::kOverwrite);
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

//************************** Function to correct for non-flat rapidity distribution *********************************
void CorrectForNonFlatRapidity(TH2F* histCorr, TH2F* histOr, TList* list) {
    
    Int_t       nBinsPt                     = histOr->GetNbinsX();
    TH1F**      ySlice                      = new TH1F*[nBinsPt];
    TH1F*       ySliceTemp                  = NULL;
    TF1**       yFit                        = new TF1*[nBinsPt];
    Double_t    nEntriesFirstSlice          = 0.;
    Double_t    fitdNdy                     = 0.;
    Double_t    fitdNdyErr                  = 0.;
    Double_t    histoMindNdy                = 0.;
    Double_t    histoMindNdyErr             = 0.;
    Double_t    histoMaxdNdy                = 0.;
    Double_t    histoMaxdNdyErr             = 0.;
    Double_t    corrFactor                  = 0.;
    Int_t       fitStatus                   = 0;
    
    for (Int_t ptBin=1; ptBin<nBinsPt+1; ptBin++) {

        // initialize
        ySlice[ptBin-1]                     = NULL;
        yFit[ptBin-1]                       = NULL;
        
        // project y distribution
        ySlice[ptBin-1]                     = (TH1F*)histOr->ProjectionY(Form("%s_ptBin%d", histOr->GetName(), ptBin),ptBin,ptBin,"e");
        ySlice[ptBin-1]->Sumw2();
        
        // test for entries
        if (ptBin==1) nEntriesFirstSlice    = ySlice[ptBin-1]->GetEntries();
        if (ySlice[ptBin-1]->GetEntries()<nEntriesFirstSlice*0.1) continue;

        // fit y distribution with constant
        yFit[ptBin-1]                       = new TF1(Form("%s_ptBin%d_fit", histOr->GetName(), ptBin), "[0]", -fRapidity, fRapidity);
        TFitResultPtr fitResult             = ySlice[ptBin-1]->Fit(yFit[ptBin-1], "SQNRIME");
        fitStatus                           = fitResult;
        if (!(fitStatus==0 || fitStatus>=1000)) continue;
        
        // get fit value
        fitdNdy                             = yFit[ptBin-1]->GetParameter(0);
        fitdNdyErr                          = yFit[ptBin-1]->GetParError(0);
        
        // test for maximum deviation from fit value in histogram
        histoMindNdy                        = ySlice[ptBin-1]->GetMinimum(0.);
        histoMindNdyErr                     = ySlice[ptBin-1]->GetBinError(GetMinimumBinAboveThreshold(ySlice[ptBin-1],0.));
        histoMaxdNdy                        = ySlice[ptBin-1]->GetMaximum();
        histoMaxdNdyErr                     = ySlice[ptBin-1]->GetBinError(ySlice[ptBin-1]->GetMaximumBin());
        if ( ((fitdNdy-fitdNdyErr)<(histoMindNdy+histoMindNdyErr)) || ((fitdNdy+fitdNdyErr)>(histoMaxdNdy-histoMaxdNdyErr)) ) continue;
        
        // correct y vs. pt hist for non-flat y dist.
        for (Int_t yBin=1; yBin<histOr->GetNbinsY()+1; yBin++) {
            
            if (!histOr->GetBinContent(ptBin, yBin)) continue;
            
            corrFactor                      = fitdNdy / histOr->GetBinContent(ptBin, yBin);
            histCorr->SetBinContent(ptBin,  yBin, histOr->GetBinContent(ptBin, yBin)*corrFactor);
            histCorr->SetBinError(ptBin,    yBin, histOr->GetBinError(ptBin, yBin)*corrFactor);
        }
        
        // add to list
        list->Add(ySlice[ptBin-1]);
        list->Add(yFit[ptBin-1]);
    }
}

//************************** Function to correct for non-flat rapidity distribution *********************************
Int_t GetMinimumBinAboveThreshold(TH1F* hist, Double_t thres) {

    TH1F*       histTemp        = (TH1F*)hist->Clone("histTemp");
    Double_t    tempBinContent  = 0.;
    
    for (Int_t i=1; i<hist->GetNbinsX()+1; i++) {
        tempBinContent          = histTemp->GetBinContent(i);
        if (tempBinContent>thres)   histTemp->SetBinContent(i, 1/tempBinContent);
        else                        histTemp->SetBinContent(i, 0);
    }
    
    return histTemp->GetMaximumBin();
}


