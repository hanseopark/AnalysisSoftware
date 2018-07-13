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
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/ExtractSignalBinning.h"

void PrepareSecondaries(    TString     meson                       = "",
                            TString     nameFileCocktail            = "",
                            TString     suffix                      = "eps",
                            TString     cutSelection                = "",
                            TString     option                      = "",
                            TString     directphotonPlots           = "",
                            TString     rapidityAndeta              = "0.85",
                            TString     period                      = "",
                            Int_t       numberOfBins                = 30,
                            Int_t       mode                        = 0,
                            Bool_t      producePlotsInOrPtRange     = kFALSE,
                            Bool_t      doRapidityCorrection        = kFALSE,
                            Bool_t      producePlotsForThesis       = kFALSE,
                            Int_t       nMotherParticleToAnalyseCur = 16
                     ) {

    gROOT->Reset();

    //************************** Set number of particles to be analyzed *********************************************
    nMotherParticleToAnalyse                                    = nMotherParticleToAnalyseCur;

    //************************** Set general style settings *********************************************************
    StyleSettingsThesis();
    SetPlotStyle();

    //************************** Set output directory ***************************************************************
    TString outputDir                                           = Form("%s/%s/%s/PrepareSecondaries",cutSelection.Data(),option.Data(),suffix.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    //************************** Set global variables ***************************************************************
    fEnergyFlag                                                 = option;
    fPeriodFlag                                                 = period;
    fdirectphoton                                               = directphotonPlots;
    fSuffix                                                     = suffix;
    fMode                                                       = mode;
    frapidityAndeta                                             = rapidityAndeta;
    cout << "Pictures are saved as " << suffix.Data() << endl;
    TObjArray *tempArr  = rapidityAndeta.Tokenize("_");
    if(tempArr->GetEntries()<1){
        cout << "rapidity not set" << endl;
        delete tempArr;
    } else if(tempArr->GetEntries()==1){
        fRapidity           = ((TString)((TObjString*)tempArr->At(0))->GetString()).Atof();
    } else if(tempArr->GetEntries()==2){
        fRapidity           = ((TString)((TObjString*)tempArr->At(0))->GetString()).Atof();
        fPseudoRapidity     = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atof();
    }
    Double_t rapidity = fRapidity;
    Double_t pseudorapidity = fPseudoRapidity;
    //************************** direct photon ana ******************************************************************
    Bool_t doSecondaryGamma                                     = kTRUE;

    //***************************** Separate cutstrings *************************************************************
    if(cutSelection.Length() == 0){
        cout<<"ERROR: Cut selection is not set, please do!"<<endl;
        return;
    }
    fCutSelection                                               = cutSelection;
    ReturnSeparatedCutNumberAdvanced(fCutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, fMode);

    Bool_t isPCM                                                = 0;
    Bool_t isCalo                                               = 0;
    if ( mode == 0 || mode == 2 || mode == 3 || mode == 13)
        isPCM                                                   = 1;
    if ( mode == 2 || mode == 3 || mode == 4 || mode == 5 || mode == 10 || mode == 12 )
        isCalo                                                  = 1;

    //***************************** Calculate scaling factor ********************************************************
    eventNormScalingFactor                                      = ReturnCocktailNormalization(fEnergyFlag, fEventCutSelection);

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
    if(nameFileCocktail.Contains("2050")) nameFileCocktail.ReplaceAll("2050","2040");
    TFile* fileCocktail     = new TFile(nameFileCocktail.Data());
    TDirectoryFile* topDirCocktail                              = (TDirectoryFile*)fileCocktail->Get("HadronicCocktailMC");
    if (!topDirCocktail) {
        cout << "ERROR: TopDirCocktail not found!" << endl;
        return;
    }
    TList* histoListCocktail                                    = NULL;
    if (fAnalyzedMeson.CompareTo("Pi0") == 0) {
        cout << "searching for " << Form("HadronicCocktailMC_pi0_%s", rapidityAndeta.Data()) << endl;
        histoListCocktail                                       = (TList*)topDirCocktail->Get(Form("HadronicCocktailMC_pi0_%s", rapidityAndeta.Data()));
        if (!histoListCocktail) {
            cout << "ERROR: Folder with rapidity " << rapidityAndeta.Data() << " not contained in cocktail file!" << endl;
            return;
        }
    } else if (fAnalyzedMeson.CompareTo("Eta") == 0) {
        cout << "searching for " << Form("HadronicCocktailMC_eta_%s", rapidityAndeta.Data()) << endl;
        histoListCocktail                                       = (TList*)topDirCocktail->Get(Form("HadronicCocktailMC_eta_%s", rapidityAndeta.Data()));
        if (!histoListCocktail) {
            cout << "ERROR: Folder with rapidity " << rapidityAndeta.Data() << " not contained in cocktail file!" << endl;
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
    // initialize mt scaling factors
    for (Int_t i=0; i<nMotherParticles; i++) mtScaleFactor[i]   = -9999.;
    // read in mt scaling factors
    TString tempBinLabel                                        = "";
    for (Int_t i=1; i<histMtScalingFactors->GetNbinsX()+1; i++) {
        tempBinLabel                                            = (TString)histMtScalingFactors->GetXaxis()->GetBinLabel(i);
        for (Int_t j=0; j<nMotherParticleToAnalyse; j++) {
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
            for (Int_t j=0; j<nMotherParticleToAnalyse; j++) {
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
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        paramTemp                                               = (TF1*)cocktailSettingsList->FindObject(Form("%s_pt", motherParticlesPDG[i].Data()));
        paramMtTemp                                             = (TF1*)cocktailSettingsList->FindObject(Form("%s_pt_mtScaled", motherParticlesPDG[i].Data()));
        if (paramTemp)  cocktailInputParametrizations[i]        = new TF1(*paramTemp);
        else            cocktailInputParametrizations[i]        = NULL;
        if (paramMtTemp)
            cocktailInputParametrizationsMtScaled[i]            = new TF1(*paramMtTemp);
        else {
            if (cocktailInputParametrizationPi0) {
                cout << motherParticlesPDG[i].Atoi() << endl;
                cocktailInputParametrizationsMtScaled[i]        = (TF1*)MtScaledParam(cocktailInputParametrizationPi0, motherParticlesPDG[i].Atoi(), 111, mtScaleFactor[i], kFALSE, kTRUE);
                if (cocktailInputParametrizationsMtScaled[i])
                    cocktailInputParametrizationsMtScaled[i]->SetName(Form("%s_pt_mtScaled", motherParticlesPDG[i].Data()));
            } else
                cocktailInputParametrizationsMtScaled[i]        = NULL;
        }
    }
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        for (Int_t j=0; j<18; j++) {
            decayChannelsBR[i][j]                               = 0.;
        }
    }

    histoDecayChannelsBR                                        = new TH1F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        if (!hasMother[i]) continue;
        histoDecayChannelsBR[i]                                 = (TH1F*)cocktailSettingsList->FindObject(Form("PythiaBR_%s", motherParticles[i].Data()));
        for (Int_t j=0; j<18; j++) {
            if (histoDecayChannelsBR[i]->GetBinContent(j+2))
                decayChannelsBR[i][j]                           = histoDecayChannelsBR[i]->GetBinContent(j+2);
        }
    }
    cout << __LINE__ << endl;
    //***************************** ranges and scaling factors ******************************************************
    Double_t deltaPtGen                                         = ptGenMax-ptGenMin;
    Double_t deltaRap                                           = 2*rapidity;
    Double_t deltaRapGen                                        = 2.0;
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
        minPhiCalo                                              = AnalyseClusterMinPhiCut(CutNumberToInteger(phiMinCut));
        maxPhiCalo                                              = AnalyseClusterMaxPhiCut(CutNumberToInteger(phiMaxCut));
    }

    Double_t scalingEta                                         = 1.;
    Double_t scalingPhi                                         = 1.;
    if (isCalo && !isPCM) {
        scalingEta                                              = deltaEtaCalo/deltaEta;
        scalingPhi                                              = deltaPhiCalo/deltaPhi;
    }

    cout << "========================================"  << endl;
    cout << "deltaRapGen      = "     << deltaRapGen    << endl;
    cout << "deltaRap         = "     << deltaRap       << endl;
    cout << "deltaPtGen       = "     << deltaPtGen     << endl;
    cout << "deltaPt          = "     << ptMax - ptMin  << endl;
    if (isPCM)
        cout << "deltaEta         = "     << deltaEta       << endl;
    if (isCalo && !isPCM)
        cout << "deltaEtaCalo     = " << deltaEtaCalo   << endl;
    if (isPCM)
        cout << "deltaPhi         = "     << deltaPhi       << endl;
    if (isCalo && !isPCM)
        cout << "deltaPhiCalo     = " << deltaPhiCalo   << endl;
    cout << "add. scaling eta = "     << scalingEta     << " (for photons only)" << endl;
    cout << "add. scaling phi = "     << scalingPhi     << " (for photons only)" << endl;
    cout << "========================================"  << endl;

    //***************************** pt range for plots **************************************************************
    if (!producePlotsInOrPtRange) {
        ptPlotMin                                               = ptMin;
        ptPlotMax                                               = ptMax;
    } else {
        ptPlotMin                                               = ptGenMin;
        ptPlotMax                                               = ptGenMax;
    }
    if (ptPlotMin == 0 || ptPlotMin < 3e-1) ptPlotMin           = 3e-1;

    //***************************** Get number of spectra ***********************************************************
    Int_t nSpectra                                              = 0;
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++)
        if (hasMother[i]) nSpectra++;
    Int_t nRows                                                 = 0;
    if (nSpectra%6 == 0) nRows                                  = nSpectra/6 + 1;
    else nRows                                                  = (nSpectra+1)/6 + 1;
    cout << __LINE__ << endl;
    //***************************** Get number of events (cocktail) *************************************************
    histoNEvents                                                = (TH1F*)histoListCocktail->FindObject("NEvents");
    nEvents                                                     = histoNEvents->GetEntries();
    cout << nEvents << " events" << endl;
    cout << __LINE__ << endl;
    //***************************** Read histograms from cocktail file **********************************************
    // mesons
    histoDecayChannels                                          = new TH1F*[nMotherParticles];
    histoMesonDaughterPtY                                       = new TH2F*[nMotherParticles];
    histoMesonDaughterPtPhi                                     = new TH2F*[nMotherParticles];
    histoMesonMotherPtY                                         = new TH2F*[nMotherParticles];
    histoMesonMotherPtPhi                                       = new TH2F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        if (hasMother[i]) {
            histoDecayChannels[i]                               = (TH1F*)histoListCocktail->FindObject(Form("DecayChannels_%s",motherParticles[i].Data()));

            histoMesonDaughterPtY[i]                            = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_%s_From_%s",fAnalyzedMeson.Data(),motherParticles[i].Data()));
            histoMesonDaughterPtY[i]->SetName(Form("%s_From_%s_Pt_Y_OrBin",fAnalyzedMeson.Data(),motherParticles[i].Data()));
            histoMesonDaughterPtY[i]->Sumw2();
            histoMesonDaughterPtY[i]->Scale(eventNormScalingFactor);

            histoMesonDaughterPtPhi[i]                          = (TH2F*)histoListCocktail->FindObject(Form("Pt_Phi_%s_From_%s",fAnalyzedMeson.Data(),motherParticles[i].Data()));
            histoMesonDaughterPtPhi[i]->SetName(Form("%s_From_%s_Pt_Phi_OrBin",fAnalyzedMeson.Data(),motherParticles[i].Data()));
            histoMesonDaughterPtPhi[i]->Sumw2();
            histoMesonDaughterPtPhi[i]->Scale(eventNormScalingFactor);

            histoMesonMotherPtY[i]                              = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_%s",motherParticles[i].Data()));
            histoMesonMotherPtY[i]->SetName(Form("%s_Pt_Y_OrBin", motherParticles[i].Data()));
            histoMesonMotherPtY[i]->Sumw2();
            histoMesonMotherPtY[i]->Scale(eventNormScalingFactor);

            histoMesonMotherPtPhi[i]                            = (TH2F*)histoListCocktail->FindObject(Form("Pt_Phi_%s",motherParticles[i].Data()));
            histoMesonMotherPtPhi[i]->SetName(Form("%s_Pt_Phi_OrBin", motherParticles[i].Data()));
            histoMesonMotherPtPhi[i]->Sumw2();
            histoMesonMotherPtPhi[i]->Scale(eventNormScalingFactor);
        } else {
            histoDecayChannels[i]                               = NULL;
            histoMesonDaughterPtY[i]                            = NULL;
            histoMesonDaughterPtPhi[i]                          = NULL;
            histoMesonMotherPtY[i]                              = NULL;
            histoMesonMotherPtPhi[i]                            = NULL;
        }
    }
    cout << __LINE__ << endl;
    // gammas
    histoGammaFromXFromMotherPtY                                = new TH2F*[nMotherParticles];
    histoGammaFromXFromMotherPtPhi                              = new TH2F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        if (doSecondaryGamma) {
            if (hasMother[i]) {
                histoGammaFromXFromMotherPtY[i]                 = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_Gamma_From_X_From_%s",motherParticles[i].Data()));
                if (histoGammaFromXFromMotherPtY[i]) {
                    histoGammaFromXFromMotherPtY[i]->SetName(Form("Gamma_From_X_From_%s_Pt_Y_OrBin", motherParticles[i].Data()));
                    histoGammaFromXFromMotherPtY[i]->Sumw2();
                    histoGammaFromXFromMotherPtY[i]->Scale(eventNormScalingFactor);
                }

                histoGammaFromXFromMotherPtPhi[i]               = (TH2F*)histoListCocktail->FindObject(Form("Pt_Phi_Gamma_From_X_From_%s",motherParticles[i].Data()));
                if (histoGammaFromXFromMotherPtPhi[i]) {
                    histoGammaFromXFromMotherPtPhi[i]->SetName(Form("Gamma_From_X_From_%s_Pt_Phi_OrBin", motherParticles[i].Data()));
                    histoGammaFromXFromMotherPtPhi[i]->Sumw2();
                    histoGammaFromXFromMotherPtPhi[i]->Scale(eventNormScalingFactor);
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
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        for (Int_t j=0; j<18; j++) {
            decayChannelsLatex[i][j]                            = "";
        }
    }

    tempBinLabel                                                = "";
    Int_t counter                                               = 0;
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
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
    cout << __LINE__ << endl;
    //***************************** Correct for non-flat y distributions ********************************************
    histoMesonDaughterPtYCorr                                   = new TH2F*[nMotherParticles];
    histoGammaFromXFromMotherPtYCorr                            = new TH2F*[nMotherParticles];
    listYSlicesMesonDaughter                                    = new TList*[nMotherParticles];
    listYSlicesGammaFromXFromMother                             = new TList*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
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
    cout << __LINE__ << endl;
    //***************************** Project from 2D histograms ******************************************************
    // mesons
    histoMesonDaughterPtOrBin                                   = new TH1F*[nMotherParticles];
    histoMesonDaughterYOrBin                                    = new TH1F*[nMotherParticles];
    histoMesonDaughterPhiOrBin                                  = new TH1F*[nMotherParticles];
    histoMesonMotherPtOrBin                                     = new TH1F*[nMotherParticles];
    histoMesonMotherYOrBin                                      = new TH1F*[nMotherParticles];
    histoMesonMotherPhiOrBin                                    = new TH1F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {

        if (histoMesonDaughterPtYCorr[i]) {
            // project pt distributions (from corrected pt-y histograms, if doRapidityCorrection = true, otherwise from uncorrected one)
            histoMesonDaughterPtOrBin[i]                        = (TH1F*)histoMesonDaughterPtYCorr[i]->ProjectionX(Form("%s_From_%s_Pt_OrBin",fAnalyzedMeson.Data(),motherParticles[i].Data()),1,histoMesonDaughterPtY[i]->GetNbinsY(),"e");
            histoMesonDaughterPtOrBin[i]->Sumw2();
            histoMesonDaughterPtOrBin[i]->Scale(1/(histoMesonDaughterPtYCorr[i]->GetYaxis()->GetBinWidth(1)));
            histoMesonDaughterPtOrBin[i]->Scale(1/(histoMesonDaughterPtYCorr[i]->GetXaxis()->GetBinWidth(1)));
            histoMesonDaughterPtOrBin[i]->Scale(1/deltaPtGen);
            histoMesonDaughterPtOrBin[i]->Scale(1/deltaRapGen);
            histoMesonDaughterPtOrBin[i]->Scale(1/deltaRap);
            histoMesonDaughterPtOrBin[i]->Scale(1/deltaPhi);
            histoMesonDaughterPtOrBin[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
            SetHistogramTitles(histoMesonDaughterPtOrBin[i],"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");

            // project y distributions (from corrected pt-y histograms, if doRapidityCorrection = true, otherwise from uncorrected one)
            histoMesonDaughterYOrBin[i]                         = (TH1F*)histoMesonDaughterPtYCorr[i]->ProjectionY(Form("%s_From_%s_Y_OrBin",fAnalyzedMeson.Data(),motherParticles[i].Data()),1,histoMesonDaughterPtY[i]->GetNbinsX(),"e");
            histoMesonDaughterYOrBin[i]->Sumw2();
            histoMesonDaughterYOrBin[i]->Scale(1/(histoMesonDaughterPtYCorr[i]->GetXaxis()->GetBinWidth(1)));
            histoMesonDaughterYOrBin[i]->Scale(1/(histoMesonDaughterPtYCorr[i]->GetYaxis()->GetBinWidth(1)));
            histoMesonDaughterYOrBin[i]->Scale(1/deltaPtGen);
            histoMesonDaughterYOrBin[i]->Scale(1/deltaRapGen);
            histoMesonDaughterYOrBin[i]->Scale(1/deltaPhi);
            SetHistogramTitles(histoMesonDaughterYOrBin[i],"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
        } else {
            histoMesonDaughterPtOrBin[i]                        = NULL;
            histoMesonDaughterYOrBin[i]                         = NULL;
        }

        // project phi distributions
        if (histoMesonDaughterPtPhi[i]) {
            histoMesonDaughterPhiOrBin[i]                       = (TH1F*)histoMesonDaughterPtPhi[i]->ProjectionY(Form("%s_From_%s_Phi_OrBin",fAnalyzedMeson.Data(),motherParticles[i].Data()),1,histoMesonDaughterPtPhi[i]->GetNbinsX(),"e");
            histoMesonDaughterPhiOrBin[i]->Sumw2();
            histoMesonDaughterPhiOrBin[i]->Scale(1/(histoMesonDaughterPtPhi[i]->GetXaxis()->GetBinWidth(1)));
            histoMesonDaughterPhiOrBin[i]->Scale(1/(histoMesonDaughterPtPhi[i]->GetYaxis()->GetBinWidth(1)));
            histoMesonDaughterPhiOrBin[i]->Scale(1/deltaPtGen);
            histoMesonDaughterPhiOrBin[i]->Scale(1/deltaRapGen);
            histoMesonDaughterPhiOrBin[i]->Scale(1/deltaPhi);
            SetHistogramTitles(histoMesonDaughterPhiOrBin[i],"","#phi","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
        } else
            histoMesonDaughterPhiOrBin[i]                       = NULL;

        if (histoMesonMotherPtY[i]) {
            // project pt distributions
            histoMesonMotherPtOrBin[i]                          = (TH1F*)histoMesonMotherPtY[i]->ProjectionX(Form("%s_Pt_OrBin",motherParticles[i].Data()),1,histoMesonMotherPtY[i]->GetNbinsY(),"e");
            histoMesonMotherPtOrBin[i]->Sumw2();
            histoMesonMotherPtOrBin[i]->Scale(1/(histoMesonMotherPtY[i]->GetYaxis()->GetBinWidth(1)));
            histoMesonMotherPtOrBin[i]->Scale(1/(histoMesonMotherPtY[i]->GetXaxis()->GetBinWidth(1)));
            histoMesonMotherPtOrBin[i]->Scale(1/deltaPtGen);
            histoMesonMotherPtOrBin[i]->Scale(1/deltaRapGen);
            histoMesonMotherPtOrBin[i]->Scale(1/deltaRap);
            histoMesonMotherPtOrBin[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
            SetHistogramTitles(histoMesonMotherPtOrBin[i],"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");

            // project y distributions
            histoMesonMotherYOrBin[i]                           = (TH1F*)histoMesonMotherPtY[i]->ProjectionY(Form("%s_Y_OrBin",motherParticles[i].Data()),1,histoMesonMotherPtY[i]->GetNbinsX(),"e");
            histoMesonMotherYOrBin[i]->Sumw2();
            histoMesonMotherYOrBin[i]->Scale(1/(histoMesonMotherPtY[i]->GetXaxis()->GetBinWidth(1)));
            histoMesonMotherYOrBin[i]->Scale(1/(histoMesonMotherPtY[i]->GetYaxis()->GetBinWidth(1)));
            histoMesonMotherYOrBin[i]->Scale(1/deltaPtGen);
            histoMesonMotherYOrBin[i]->Scale(1/deltaRapGen);
            SetHistogramTitles(histoMesonMotherYOrBin[i],"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
        } else {
            histoMesonMotherPtOrBin[i]                          = NULL;
            histoMesonMotherYOrBin[i]                           = NULL;
        }

        // project phi distributions
        if (histoMesonMotherPtPhi[i]) {
            histoMesonMotherPhiOrBin[i]                         = (TH1F*)histoMesonMotherPtPhi[i]->ProjectionY(Form("%s_Phi_OrBin",motherParticles[i].Data()),1,histoMesonMotherPtPhi[i]->GetNbinsX(),"e");
            histoMesonMotherPhiOrBin[i]->Sumw2();
            histoMesonMotherPhiOrBin[i]->Scale(1/(histoMesonMotherPtPhi[i]->GetXaxis()->GetBinWidth(1)));
            histoMesonMotherPhiOrBin[i]->Scale(1/(histoMesonMotherPtPhi[i]->GetYaxis()->GetBinWidth(1)));
            histoMesonMotherPhiOrBin[i]->Scale(1/deltaPtGen);
            histoMesonMotherPhiOrBin[i]->Scale(1/deltaRapGen);
            SetHistogramTitles(histoMesonMotherPhiOrBin[i],"","#phi","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
        } else
            histoMesonMotherPhiOrBin[i]                         = NULL;
    }

    if (!histoMesonMotherPtOrBin[1]) {
        cout << "ERROR: Didn't get K0s pt spectrum, returning!" << endl;
        return;
    }
    cout << __LINE__ << endl;
    // gammas
    histoGammaFromXFromMotherPtOrBin                            = new TH1F*[nMotherParticles];
    histoGammaFromXFromMotherYOrBin                             = new TH1F*[nMotherParticles];
    histoGammaFromXFromMotherPhiOrBin                           = new TH1F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        if (doSecondaryGamma) {

            if (histoGammaFromXFromMotherPtYCorr[i]) {
                // project pt distributions (from corrected pt-y histograms, if doRapidityCorrection = true, otherwise from uncorrected one)
                histoGammaFromXFromMotherPtOrBin[i]             = (TH1F*)histoGammaFromXFromMotherPtYCorr[i]->ProjectionX(Form("Gamma_From_X_From_%s_Pt_OrBin",motherParticles[i].Data()),1,histoGammaFromXFromMotherPtY[i]->GetNbinsY(),"e");
                histoGammaFromXFromMotherPtOrBin[i]->Sumw2();
                histoGammaFromXFromMotherPtOrBin[i]->Scale(1/(histoGammaFromXFromMotherPtYCorr[i]->GetYaxis()->GetBinWidth(1)));
                histoGammaFromXFromMotherPtOrBin[i]->Scale(1/(histoGammaFromXFromMotherPtYCorr[i]->GetXaxis()->GetBinWidth(1)));
                histoGammaFromXFromMotherPtOrBin[i]->Scale(1/deltaPtGen);
                histoGammaFromXFromMotherPtOrBin[i]->Scale(1/deltaRapGen);
                histoGammaFromXFromMotherPtOrBin[i]->Scale(1/deltaRap);
                histoGammaFromXFromMotherPtOrBin[i]->Scale(1/deltaPhi);
                histoGammaFromXFromMotherPtOrBin[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
                SetHistogramTitles(histoGammaFromXFromMotherPtOrBin[i],"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");

                // project y distributions (from corrected pt-y histograms, if doRapidityCorrection = true, otherwise from uncorrected one)
                histoGammaFromXFromMotherYOrBin[i]              = (TH1F*)histoGammaFromXFromMotherPtYCorr[i]->ProjectionY(Form("Gamma_From_X_From_%s_Y_OrBin",motherParticles[i].Data()),1,histoGammaFromXFromMotherPtY[i]->GetNbinsX(),"e");
                histoGammaFromXFromMotherYOrBin[i]->Sumw2();
                histoGammaFromXFromMotherYOrBin[i]->Scale(1/(histoGammaFromXFromMotherPtYCorr[i]->GetXaxis()->GetBinWidth(1)));
                histoGammaFromXFromMotherYOrBin[i]->Scale(1/(histoGammaFromXFromMotherPtYCorr[i]->GetYaxis()->GetBinWidth(1)));
                histoGammaFromXFromMotherYOrBin[i]->Scale(1/deltaPtGen);
                histoGammaFromXFromMotherYOrBin[i]->Scale(1/deltaRapGen);
                histoGammaFromXFromMotherYOrBin[i]->Scale(1/deltaPhi);
                SetHistogramTitles(histoGammaFromXFromMotherYOrBin[i],"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            } else {
                histoGammaFromXFromMotherPtOrBin[i]             = NULL;
                histoGammaFromXFromMotherYOrBin[i]              = NULL;
            }

            // project phi distributions
            if (histoGammaFromXFromMotherPtPhi[i]) {
                histoGammaFromXFromMotherPhiOrBin[i]            = (TH1F*)histoGammaFromXFromMotherPtPhi[i]->ProjectionY(Form("Gamma_From_X_From_%s_Phi_OrBin",motherParticles[i].Data()),1,histoGammaFromXFromMotherPtPhi[i]->GetNbinsX(),"e");
                histoGammaFromXFromMotherPhiOrBin[i]->Sumw2();
                histoGammaFromXFromMotherPhiOrBin[i]->Scale(1/(histoGammaFromXFromMotherPtPhi[i]->GetXaxis()->GetBinWidth(1)));
                histoGammaFromXFromMotherPhiOrBin[i]->Scale(1/(histoGammaFromXFromMotherPtPhi[i]->GetYaxis()->GetBinWidth(1)));
                histoGammaFromXFromMotherPhiOrBin[i]->Scale(1/deltaPtGen);
                histoGammaFromXFromMotherPhiOrBin[i]->Scale(1/deltaRapGen);
                histoGammaFromXFromMotherPhiOrBin[i]->Scale(1/deltaPhi);
                SetHistogramTitles(histoGammaFromXFromMotherPhiOrBin[i],"","#phi","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            } else {
                histoGammaFromXFromMotherPhiOrBin[i]            = NULL;
            }
        } else {
            histoGammaFromXFromMotherPtOrBin[i]                 = NULL;
            histoGammaFromXFromMotherYOrBin[i]                  = NULL;
            histoGammaFromXFromMotherPhiOrBin[i]                = NULL;
        }
    }

    //***************************** Scale spectra and params ********************************************************
    Double_t factorNEvents                                 = 1./nEvents;
    Double_t factorDecay                                   = 1.;
    TF1* corrFactorFromDecayLength = new TF1("decayLaw","1-exp(-x/[0])",0,500);
    cout << "========================================"  << endl;
    cout << "correction factors from decay length (fraction of particles decayed up to R = '" << ReturnMeanR(fMode) << "cm'): " << endl;
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        // calculate fraction of given particles which are decayed up to mean R of reconstruction
        if(motherParticles_ctau[i]<0.01) factorDecay = 1.;
        else{
          corrFactorFromDecayLength->SetParameter(0,motherParticles_ctau[i]);
          factorDecay = corrFactorFromDecayLength->Eval(ReturnMeanR(fMode));
        }
        cout << motherParticles[i] << "\t--> " << factorDecay << endl;
        motherFactorDecayLength[i]                          = factorDecay;

        // spectra
        if (histoMesonDaughterPtOrBin[i])           histoMesonDaughterPtOrBin[i]->Scale(        factorNEvents*factorDecay);
        if (histoMesonDaughterYOrBin[i])            histoMesonDaughterYOrBin[i]->Scale(         factorNEvents*factorDecay);
        if (histoMesonDaughterPhiOrBin[i])          histoMesonDaughterPhiOrBin[i]->Scale(       factorNEvents*factorDecay);
        if (histoMesonMotherPtOrBin[i])             histoMesonMotherPtOrBin[i]->Scale(          factorNEvents*factorDecay);
        if (histoMesonMotherYOrBin[i])              histoMesonMotherYOrBin[i]->Scale(           factorNEvents*factorDecay);
        if (histoMesonMotherPhiOrBin[i])            histoMesonMotherPhiOrBin[i]->Scale(         factorNEvents*factorDecay);
        if (histoGammaFromXFromMotherPtOrBin[i])    histoGammaFromXFromMotherPtOrBin[i]->Scale( factorNEvents*scalingEta*scalingPhi*factorDecay);
        if (histoGammaFromXFromMotherYOrBin[i])     histoGammaFromXFromMotherYOrBin[i]->Scale(  factorNEvents*scalingEta*scalingPhi*factorDecay);
        if (histoGammaFromXFromMotherPhiOrBin[i])   histoGammaFromXFromMotherPhiOrBin[i]->Scale(factorNEvents*scalingEta*scalingPhi*factorDecay);

        // parametrizations
        if (cocktailInputParametrizations[i])
            cocktailInputParametrizations[i]                    = (TF1*)ScaleTF1(cocktailInputParametrizations[i], eventNormScalingFactor*factorDecay, cocktailInputParametrizations[i]->GetName());
        if (cocktailInputParametrizationsMtScaled[i])
            cocktailInputParametrizationsMtScaled[i]            = (TF1*)ScaleTF1(cocktailInputParametrizationsMtScaled[i], eventNormScalingFactor*factorDecay, cocktailInputParametrizationsMtScaled[i]->GetName());
    }
    cout << "========================================"  << endl;
    delete corrFactorFromDecayLength;
    cout << __LINE__ << endl;
    //***************************** Read histograms from cocktail input file ****************************************
    TList* cocktailInputList                                    = GetCocktailInputList(fEnergyFlag, centrality);
    if(centrality.CompareTo("20-50%")==0) cocktailInputList     = GetCocktailInputList(fEnergyFlag, "20-40%");
    histoMesonMotherCocktailInputPtMeasBin                      = new TH1F*[nCocktailInputParticles];
    histoMesonMotherPtMeasBin                                   = new TH1F*[nCocktailInputParticles];

    // load cocktail input spectra
    for (Int_t i=0; i<nCocktailInputParticles; i++) {
        histoMesonMotherCocktailInputPtMeasBin[i]               = GetCocktailInputSpectrum(cocktailInputList, i, Form("%s_CocktailInput_Pt_MeasBin",motherParticles[i+1].Data()));

        // K0l not included in cocktail input -> clone K0s
        if (i>0 && !histoMesonMotherCocktailInputPtMeasBin[i]) {
            if (cocktailInputParticles[i].CompareTo(cocktailInputParticles[i-1]) == 0 && histoMesonMotherCocktailInputPtMeasBin[i-1]) {
                histoMesonMotherCocktailInputPtMeasBin[i]       = (TH1F*)histoMesonMotherCocktailInputPtMeasBin[i-1]->Clone(Form("%s_CocktailInput_Pt_MeasBin",motherParticles[i+1].Data()));
            }
            else
                histoMesonMotherCocktailInputPtMeasBin[i] = NULL;
        }
    }

    // scale cocktail input spectra and rebin generated spectra
    for (Int_t i=0; i<nCocktailInputParticles; i++) {

        if (histoMesonMotherCocktailInputPtMeasBin[i]) {
            histoMesonMotherCocktailInputPtMeasBin[i]->Sumw2();
            histoMesonMotherCocktailInputPtMeasBin[i]->Scale(eventNormScalingFactor*motherFactorDecayLength[i+1]);

            if (histoMesonMotherPtOrBin[i+1]) {
                histoMesonMotherPtMeasBin[i]                    = (TH1F*)histoMesonMotherPtOrBin[i+1]->Clone(Form("%s_Pt_MeasBin",motherParticles[i+1].Data()));
                histoMesonMotherPtMeasBin[i]->Sumw2();
                RebinSpectrum(histoMesonMotherPtMeasBin[i], histoMesonMotherCocktailInputPtMeasBin[i], "");
            } else {
                histoMesonMotherPtMeasBin[i]                    = NULL;
            }
        } else {
            histoMesonMotherPtMeasBin[i]                        = NULL;
        }
    }
    cout << __LINE__ << endl;
    //***************************** calculate (pi0 from X)/pi0_param ************************************************
    histoRatioPi0FromXToPi0Param                                = new TH1F*[nMotherParticles];
    if (fAnalyzedMeson.CompareTo("Pi0") == 0) {
        if (cocktailInputParametrizationPi0) {
            for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
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
            for (Int_t i=0; i<nMotherParticleToAnalyse; i++)
                histoRatioPi0FromXToPi0Param[i]                 = NULL;
        }
    } else {
        for (Int_t i=0; i<nMotherParticleToAnalyse; i++)
            histoRatioPi0FromXToPi0Param[i]                     = NULL;
    }

    TH1D* dummyHist                                             = NULL;
    //***************************** Plot cocktail mothers (pt) ******************************************************
    cout << __LINE__ << endl;
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

    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        if (histoMesonMotherPtOrBin[i]) {
            DrawGammaSetMarker(         histoMesonMotherPtOrBin[i], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendMothers->AddEntry(    histoMesonMotherPtOrBin[i], Form("%s", motherParticlesLatex[i].Data()), "l");
            histoMesonMotherPtOrBin[i]->SetLineWidth(2);
            histoMesonMotherPtOrBin[i]->Draw("csamehist");
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

    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        if (histoMesonMotherPtOrBin[i]) {
            if (cocktailInputParametrizations[i]) {
                cocktailInputParametrizations[i]->SetLineColor(cocktailColor[i]);
                cocktailInputParametrizations[i]->SetLineStyle(2);
                cocktailInputParametrizations[i]->SetLineWidth(2);
                cocktailInputParametrizations[i]->Draw("same");
            }
            legendMothersParam->AddEntry(histoMesonMotherPtOrBin[i], Form("%s", motherParticlesLatex[i].Data()), "l");
            histoMesonMotherPtOrBin[i]->Draw("csamehist");
        }
    }
    legendMothersParam->Draw("same");

    PutProcessLabelAndEnergyOnPlot(                 0.22, 0.22, 0.032, cent, textMeasurement, "", 42, 0.03);
    if (producePlotsForThesis) PutThisThesisLabel(  0.22, 0.17, 0.032, 0.03, 1.25, 42);
    else PutALICESimulationLabel(                   0.22, 0.17, 0.032, 0.03, 1.25, 42);

    canvasMothersParam->SaveAs(Form("%s/CocktailMothersInclParam_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendMothersParam;
    delete canvasMothersParam;

    //***************************** Plot cocktail mothers (y) *******************************************************
    TCanvas *canvasMothersY                                     = new TCanvas("canvasMothersY","",1100,1200);
    DrawGammaCanvasSettings(canvasMothersY, 0.12, 0.025, 0.01, 0.075);
    canvasMothersY->SetLogy();

    TLegend* legendMothersY                                     = GetAndSetLegend2(0.2, 0.98-(0.04*nRows), 0.95, 0.98, 40, 6);
    legendMothersY->SetBorderSize(0);
    dummyHist                                                   = new TH1D("dummyHist", "", 1000, -fRapidity, fRapidity);
    SetHistogramm(dummyHist, "y", "#frac{1}{N_{ev}} #frac{d#it{N}}{dy}", 1e-5, 1, 0.9, 1.25);
    dummyHist->Draw();

    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        if (histoMesonMotherYOrBin[i]) {
            DrawGammaSetMarker(         histoMesonMotherYOrBin[i], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendMothersY->AddEntry(   histoMesonMotherYOrBin[i], Form("%s", motherParticlesLatex[i].Data()), "p");
            histoMesonMotherYOrBin[i]->Draw("same");
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
    SetHistogramm(dummyHist, "#phi", "#frac{1}{N_{ev}} #frac{d#it{N}}{d#phi}", 1e-5, 1, 0.9, 1.25);
    dummyHist->Draw();

    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        if (histoMesonMotherPhiOrBin[i]) {
            DrawGammaSetMarker(         histoMesonMotherPhiOrBin[i], cocktailMarker[i], 1, cocktailColor[i],  cocktailColor[i]);
            legendMothersPhi->AddEntry( histoMesonMotherPhiOrBin[i], Form("%s", motherParticlesLatex[i].Data()), "p");
            histoMesonMotherPhiOrBin[i]->Draw("same");
        }
    }
    legendMothersPhi->Draw("same");

    PutProcessLabelAndEnergyOnPlot(                 0.2, 0.22, 0.032, cent, textMeasurement, "", 42, 0.03);
    if (producePlotsForThesis) PutThisThesisLabel(  0.2, 0.17, 0.032, 0.03, 1.25, 42);
    else PutALICESimulationLabel(                   0.2, 0.17, 0.032, 0.03, 1.25, 42);

    canvasMothersPhi->SaveAs(Form("%s/CocktailMothersPhi_%.2f_%s.%s",outputDir.Data(),fRapidity,cutSelection.Data(),suffix.Data()));
    delete legendMothersPhi;
    delete canvasMothersPhi;

    //***************************** Plot K0s, K0l and Lambda cross check ********************************************
    TCanvas* canvasInputCrossCheck                              = NULL;
    TLegend* legendInputCrossCheck                              = NULL;
    TPad* padInputCrossCheck                                    = NULL;
    TPad* padInputCrossCheckRatio                               = NULL;
    dummyHist                                                   = NULL;
    TH1D* dummyHistRatio                                        = NULL;
    TH1D* tempRatio1                                            = NULL;
    TH1D* tempRatio2                                            = NULL;
    TH1D* tempRatio3                                            = NULL;
    TH1D* tempRatio4                                            = NULL;
    for (Int_t i=0; i<nCocktailInputParticles; i++) {
        if (histoMesonMotherCocktailInputPtMeasBin[i] && histoMesonMotherPtMeasBin[i] && (cocktailInputParametrizations[i+1] || cocktailInputParametrizationsMtScaled[i+1])) {

            // calc. ratios
            if (cocktailInputParametrizations[i+1]) {
                tempRatio1                                      = (TH1D*)CalculateRatioToTF1((TH1D*)histoMesonMotherCocktailInputPtMeasBin[i],  cocktailInputParametrizations[i+1]);
                tempRatio2                                      = (TH1D*)CalculateRatioToTF1((TH1D*)histoMesonMotherPtMeasBin[i],               cocktailInputParametrizations[i+1]);
            }
            if (cocktailInputParametrizationsMtScaled[i+1]) {
                tempRatio3                                      = (TH1D*)CalculateRatioToTF1((TH1D*)histoMesonMotherCocktailInputPtMeasBin[i],  cocktailInputParametrizationsMtScaled[i+1]);
                tempRatio4                                      = (TH1D*)CalculateRatioToTF1((TH1D*)histoMesonMotherPtMeasBin[i],               cocktailInputParametrizationsMtScaled[i+1]);
            }

            // start plotting
            canvasInputCrossCheck                               = new TCanvas("canvasInputCrossCheck","",1100,1200);
            padInputCrossCheck                                  = new TPad("padInputCrossCheck", "", 0., 0.25, 1., 1.,-1, -1, -2);
            padInputCrossCheckRatio                             = new TPad("padInputCrossCheckRatio", "", 0., 0., 1., 0.25,-1, -1, -2);
            legendInputCrossCheck                               = GetAndSetLegend2(0.4, 0.9-(0.048*5), 0.9, 0.9, 40);
            legendInputCrossCheck->SetBorderSize(0);

            DrawGammaCanvasSettings(canvasInputCrossCheck, 0.165, 0.015, 0.025, 0.25);
            DrawGammaPadSettings(padInputCrossCheck,       0.165, 0.015, 0.025, 0.);
            DrawGammaPadSettings(padInputCrossCheckRatio,  0.165, 0.015, 0.0, 0.25);

            padInputCrossCheck->Draw();
            padInputCrossCheck->SetLogy();
            padInputCrossCheck->SetLogx();

            padInputCrossCheckRatio->Draw();
            padInputCrossCheckRatio->SetLogx();

            dummyHist                                           = new TH1D("dummyHist", "", 1000, ptPlotMin, ptPlotMax);
            dummyHistRatio                                      = new TH1D("dummyHist", "", 1000, ptPlotMin, ptPlotMax);

            // spectra and params.
            padInputCrossCheck->cd();

            SetHistogramm(dummyHist, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})", 5e-8, 1e2, 1.0, 1.8);
            dummyHist->Draw();

            DrawGammaSetMarker(histoMesonMotherCocktailInputPtMeasBin[i],   24, 2.0, kBlack,  kBlack);
            DrawGammaSetMarker(histoMesonMotherPtMeasBin[i],                27, 2.0, kBlack,  kBlack);

            histoMesonMotherCocktailInputPtMeasBin[i]->Draw("e1,same");
            histoMesonMotherPtMeasBin[i]->Draw("e1,same");

            legendInputCrossCheck->AddEntry(histoMesonMotherCocktailInputPtMeasBin[i], Form("%s cocktail input (#times %.3e)", motherParticlesLatex[i+1].Data(), motherFactorDecayLength[i+1]), "p");
            legendInputCrossCheck->AddEntry(histoMesonMotherPtMeasBin[i], Form("%s generated (#times %.3e)", motherParticlesLatex[i+1].Data(), motherFactorDecayLength[i+1]), "p");

            if (cocktailInputParametrizations[i+1]) {
                cocktailInputParametrizations[i+1]->SetLineColor(kBlue-2);
                cocktailInputParametrizations[i+1]->SetLineStyle(1);
                cocktailInputParametrizations[i+1]->SetLineWidth(2);
                cocktailInputParametrizations[i+1]->Draw("same");
                legendInputCrossCheck->AddEntry(cocktailInputParametrizations[i+1], Form("%s param. (#times %.3e)", motherParticlesLatex[i+1].Data(), motherFactorDecayLength[i+1]), "l");
            }

            if (cocktailInputParametrizationsMtScaled[i+1]) {
                cocktailInputParametrizationsMtScaled[i+1]->SetLineColor(kOrange+2);
                cocktailInputParametrizationsMtScaled[i+1]->SetLineStyle(1);
                cocktailInputParametrizationsMtScaled[i+1]->SetLineWidth(2);
                cocktailInputParametrizationsMtScaled[i+1]->Draw("same");
                legendInputCrossCheck->AddEntry(cocktailInputParametrizationsMtScaled[i+1], Form("%s param. m_{T} scaled (#times %.3e)", motherParticlesLatex[i+1].Data(), motherFactorDecayLength[i+1]), "l");
            }

            legendInputCrossCheck->Draw("same");

            PutProcessLabelAndEnergyOnPlot(                 0.22, 0.30, 0.032, cent, textMeasurement, "", 42, 0.03);
            if (producePlotsForThesis) PutThisThesisLabel(  0.22, 0.25, 0.032, 0.03, 1.25, 42);
            else PutALICESimulationLabel(                   0.22, 0.25, 0.032, 0.03, 1.25, 42);

            // ratios of spectra to params.
            padInputCrossCheckRatio->cd();
            SetStyleHistoTH1ForGraphs(dummyHistRatio, "#it{p}_{T} (GeV/#it{c})","#frac{spec.}{param.}", 0.12, 0.1, 0.12, 0.1, 1.1, 0.6, 510, 505);
            dummyHistRatio->GetXaxis()->SetLabelOffset(-0.025);
            dummyHistRatio->GetYaxis()->SetRangeUser(0.5,1.5);
            dummyHistRatio->Draw();

            DrawGammaLines(ptPlotMin, ptPlotMax, 1.0, 1.0, 2.0, kGray+2, 2);
            DrawGammaLines(ptPlotMin, ptPlotMax, 1.2, 1.2, 2.0, kGray+2, 8);
            DrawGammaLines(ptPlotMin, ptPlotMax, 0.8, 0.8, 2.0, kGray+2, 8);

            if (cocktailInputParametrizations[i+1]) {
                DrawGammaSetMarker(tempRatio1, 24, 2.0, kBlue-2,  kBlue-2);
                DrawGammaSetMarker(tempRatio2, 27, 2.0, kBlue-2,  kBlue-2);
                tempRatio1->Draw("e1,same");
                tempRatio2->Draw("e1,same");
            }

            if (cocktailInputParametrizationsMtScaled[i+1]) {
                DrawGammaSetMarker(tempRatio3, 24, 2.0, kOrange+2,  kOrange+2);
                DrawGammaSetMarker(tempRatio4, 27, 2.0, kOrange+2,  kOrange+2);
                tempRatio3->Draw("e1,same");
                tempRatio4->Draw("e1,same");
            }

            canvasInputCrossCheck->SaveAs(Form("%s/InputCrossCheck%s_%.2f_%s.%s",outputDir.Data(), motherParticles[i+1].Data(),fRapidity,cutSelection.Data(),suffix.Data()));
        }
    }
    delete dummyHist;
    delete dummyHistRatio;

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
void RebinSpectrum(TH1F* Spectrum, TH1F* SpectrumForBinning, TString NewName){

    if (!SpectrumForBinning) return;

    if(NewName.CompareTo("")==0)
        NewName                     = Spectrum->GetName();

    Double_t newBinContent          = 0.;
    Double_t newBinError            = 0.;
    for (Int_t i=1; i<Spectrum->GetNbinsX()+1; i++) {
        newBinContent               = Spectrum->GetBinContent(i) * Spectrum->GetBinWidth(i);
        newBinError                 = Spectrum->GetBinError(i) * Spectrum->GetBinWidth(i);
        Spectrum->SetBinContent(i,  newBinContent);
        Spectrum->SetBinError(i,    newBinError);
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

//************************** Calculate ratio of TF1 to TH1D *********************************************************
TH1D* CalculateRatioToTF1(TH1D* hist, TF1* func) {

    if (!hist) return NULL;

    TH1D* resultHist            = (TH1D*)hist->Clone(Form("%s_to_%s", hist->GetName(), func->GetName()));

    Double_t tempValHist        = 0.;
    Double_t tempValFunc        = 0.;

    Double_t tempVal            = 0.;
    Double_t tempErr            = 0.;

    for (Int_t i=1; i<hist->GetNbinsX()+1; i++) {
        if (hist->GetBinContent(i)) {

            tempValHist         = hist->GetBinContent(i)*hist->GetBinWidth(i);
            tempValFunc         = func->Integral(hist->GetXaxis()->GetBinLowEdge(i), hist->GetXaxis()->GetBinUpEdge(i));

            tempVal             = tempValHist/tempValFunc;
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
    } else if (particleName.CompareTo("EtaPrime") == 0 || particleName.CompareTo("EtaPrim") == 0) {
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
    TString nameOutput                  = Form("%s/%s/Secondary%s%s_%s_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fAnalyzedMeson.Data(),fPeriodFlag.Data(),frapidityAndeta.Data(),fCutSelection.Data());
    cout << "INFO: writing into: " << nameOutput << endl;
    TFile *outputFile                   = new TFile(nameOutput,"RECREATE");

    // write number of events
    histoNEvents->Write("NEvents", TObject::kOverwrite);

    // write y vs pt histograms
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        if (histoMesonDaughterPtY[i])       histoMesonDaughterPtY[i]->Write(        histoMesonDaughterPtY[i]->GetName(),        TObject::kOverwrite);
        if (histoMesonDaughterPtYCorr[i])   histoMesonDaughterPtYCorr[i]->Write(    histoMesonDaughterPtYCorr[i]->GetName(),    TObject::kOverwrite);
        if (listYSlicesMesonDaughter[i])    listYSlicesMesonDaughter[i]->Write(     listYSlicesMesonDaughter[i]->GetName(),     TObject::kSingleKey);
    }

    // write projections
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        if (histoMesonDaughterPtOrBin[i])   histoMesonDaughterPtOrBin[i]->Write(    histoMesonDaughterPtOrBin[i]->GetName(),    TObject::kOverwrite);
        if (histoMesonDaughterYOrBin[i])    histoMesonDaughterYOrBin[i]->Write(     histoMesonDaughterYOrBin[i]->GetName(),     TObject::kOverwrite);
        if (histoMesonDaughterPhiOrBin[i])  histoMesonDaughterPhiOrBin[i]->Write(   histoMesonDaughterPhiOrBin[i]->GetName(),   TObject::kOverwrite);
        if (histoMesonMotherPtOrBin[i])     histoMesonMotherPtOrBin[i]->Write(      histoMesonMotherPtOrBin[i]->GetName(),      TObject::kOverwrite);
        if (histoMesonMotherYOrBin[i])      histoMesonMotherYOrBin[i]->Write(       histoMesonMotherYOrBin[i]->GetName(),       TObject::kOverwrite);
        if (histoMesonMotherPhiOrBin[i])    histoMesonMotherPhiOrBin[i]->Write(     histoMesonMotherPhiOrBin[i]->GetName(),     TObject::kOverwrite);
    }

    // write input params
    if (cocktailInputParametrizationPi0)                cocktailInputParametrizationPi0->Write(         cocktailInputParametrizationPi0->GetName(),         TObject::kOverwrite);
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        if (cocktailInputParametrizations[i])           cocktailInputParametrizations[i]->Write(        cocktailInputParametrizations[i]->GetName(),        TObject::kOverwrite);
        if (cocktailInputParametrizationsMtScaled[i])   cocktailInputParametrizationsMtScaled[i]->Write(cocktailInputParametrizationsMtScaled[i]->GetName(),TObject::kOverwrite);
    }

    // write input spectra and rebinned generate spectra
    for (Int_t i=0; i<nCocktailInputParticles; i++) {
        if (histoMesonMotherCocktailInputPtMeasBin[i])  histoMesonMotherCocktailInputPtMeasBin[i]->Write(   histoMesonMotherCocktailInputPtMeasBin[i]->GetName(),   TObject::kOverwrite);
        if (histoMesonMotherPtMeasBin[i])               histoMesonMotherPtMeasBin[i]->Write(                histoMesonMotherPtMeasBin[i]->GetName(),                TObject::kOverwrite);
    }

    // write ratio pi0 from X to pi0 param
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        if (histoRatioPi0FromXToPi0Param[i]) histoRatioPi0FromXToPi0Param[i]->Write(histoRatioPi0FromXToPi0Param[i]->GetName(), TObject::kOverwrite);
    }

    outputFile->Write();
    outputFile->Close();
}

//************************** Routine for saving histograms **********************************************************
void SavePhotonHistos() {
    TString nameOutput                  = Form("%s/%s/SecondaryGamma%s_%s_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fPeriodFlag.Data(),frapidityAndeta.Data(),fCutSelection.Data());
    cout << "INFO: writing into: " << nameOutput << endl;
    TFile *outputFile                   = new TFile(nameOutput,"RECREATE");

    // write number of events
    histoNEvents->Write("NEvents", TObject::kOverwrite);

    // write y vs pt histograms
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        if (histoGammaFromXFromMotherPtYCorr[i])    histoGammaFromXFromMotherPtYCorr[i]->Write( histoGammaFromXFromMotherPtYCorr[i]->GetName(), TObject::kOverwrite);
        if (histoGammaFromXFromMotherPtY[i])        histoGammaFromXFromMotherPtY[i]->Write(     histoGammaFromXFromMotherPtY[i]->GetName(),     TObject::kOverwrite);
        if (listYSlicesGammaFromXFromMother[i])     listYSlicesGammaFromXFromMother[i]->Write(  listYSlicesGammaFromXFromMother[i]->GetName(),  TObject::kSingleKey);
    }

    // write projections
    for (Int_t i=0; i<nMotherParticleToAnalyse; i++) {
        if (histoGammaFromXFromMotherPtOrBin[i])    histoGammaFromXFromMotherPtOrBin[i]->Write( histoGammaFromXFromMotherPtOrBin[i]->GetName(), TObject::kOverwrite);
        if (histoGammaFromXFromMotherYOrBin[i])     histoGammaFromXFromMotherYOrBin[i]->Write(  histoGammaFromXFromMotherYOrBin[i]->GetName(),  TObject::kOverwrite);
        if (histoGammaFromXFromMotherPhiOrBin[i])   histoGammaFromXFromMotherPhiOrBin[i]->Write(histoGammaFromXFromMotherPhiOrBin[i]->GetName(),TObject::kOverwrite);
    }

    // write input spectra and rebinned generate spectra
    for (Int_t i=0; i<nCocktailInputParticles; i++) {
        if (histoMesonMotherCocktailInputPtMeasBin[i])  histoMesonMotherCocktailInputPtMeasBin[i]->Write(   histoMesonMotherCocktailInputPtMeasBin[i]->GetName(),   TObject::kOverwrite);
        if (histoMesonMotherPtMeasBin[i])               histoMesonMotherPtMeasBin[i]->Write(                histoMesonMotherPtMeasBin[i]->GetName(),                TObject::kOverwrite);
    }

    outputFile->Write();
    outputFile->Close();
}

//************************** Create tex file containing BR table ****************************************************
void CreateBRTableLatex() {

    // tex file
    TString texFileName                  = Form("%s/%s/HadronicCocktail%s_%s_%s_BranchingRatioTable.tex", fCutSelection.Data(), fEnergyFlag.Data(), fPeriodFlag.Data(), frapidityAndeta.Data(), fCutSelection.Data());
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
    for (Int_t i=nMotherParticleToAnalyse-1; i>=0; i--) {
        if(hasMother[i]) {
            lastMother                  = i;
            break;
        }
    }

    // loop over particles
    for(Int_t particle=0; particle<nMotherParticleToAnalyse; particle++) {
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

//************************** Function to load correct list cocktail input from comp. file ***************************
TList* GetCocktailInputList(TString energy, TString centrality) {

    if (energy.CompareTo("") == 0) {
        cout << "no energy given!" << endl;
        return NULL;
    }

    // get file of cocktail input objects
    TString fCollSysFile                                        = "";
    TString fCollSys                                            = "";
    if (energy.Contains("PbPb")) {
        fCollSysFile                                            = "PbPb";
        fCollSys                                                = "PbPb";
    } else if (energy.Contains("pPb")) {
        fCollSysFile                                            = "PPb";
        fCollSys                                                = "pPb";
    } else  {
        fCollSysFile                                            = "PP";
        fCollSys                                                = "pp";
    }

    TString fileName                                            = Form("CocktailInput/CocktailInput%s.root", fCollSysFile.Data());
    TFile file(fileName.Data());
    if (file.IsZombie()) {
        cout << "file " << fileName.Data() << " not found!" << endl;
        return NULL;
    }

    // get list of cocktail input objects in file
    TString                         fEnergy                     = "";
    if (energy.Contains("9"))       fEnergy                     = "0.9TeV";
    else if (energy.Contains("2."))  fEnergy                     = "2.76TeV";
    else if (energy.Contains("5"))  fEnergy                     = "5TeV";
    else if (energy.Contains("7"))  fEnergy                     = "7TeV";
    else if (energy.Contains("8"))  fEnergy                     = "8TeV";
    else if (energy.Contains("13")) fEnergy                     = "13TeV";
    else {
        cout << "energy " << energy.Data() << " not recognized, failed to load correct list of cocktail input objects!" << endl;
        return NULL;
    }

    TString                                         fCentrality = "";
    if (fCollSys.CompareTo("pp")) {
        if (centrality.CompareTo("0-100%")==0)      fCentrality = "MB";
        else if (centrality.CompareTo("0-5%")==0)   fCentrality = "0005";
        else if (centrality.CompareTo("5-10%")==0)  fCentrality = "0510";
        else if (centrality.CompareTo("0-10%")==0)  fCentrality = "0010";
        else if (centrality.CompareTo("10-20%")==0) fCentrality = "1020";
        else if (centrality.CompareTo("0-20%")==0)  fCentrality = "0020";
        else if (centrality.CompareTo("20-40%")==0) fCentrality = "2040";
        else if (centrality.CompareTo("20-50%")==0) fCentrality = "2050";
        else if (centrality.CompareTo("40-60%")==0) fCentrality = "4060";
        else if (centrality.CompareTo("60-80%")==0) fCentrality = "6080";
        else if (centrality.CompareTo("80-90%")==0) fCentrality = "8090";
        else if (centrality.CompareTo("80-100%")==0)fCentrality = "80100";
        else if (centrality.CompareTo("20-30%")==0) fCentrality = "2030";
        else if (centrality.CompareTo("30-40%")==0) fCentrality = "3040";
        else {
            cout << "centrality " << centrality.Data() << " not recognized, failed to load correct list of cocktail input objects!" << endl;
            return NULL;
        }
    }

    TString                             listName                = "";
    if (fCollSys.CompareTo("pp")==0)          listName                = Form("%s_%s",  fCollSys.Data(), fEnergy.Data());
    else if (fCollSys.CompareTo("pPb")==0)    listName                = Form("%s_%s",  fEnergy.Data(), fCentrality.Data());
    else                                      listName                = Form("%s_%s_%s",  fCollSys.Data(), fEnergy.Data(), fCentrality.Data());
    TList* list                                                 = (TList*)file.Get(listName.Data());
    if (!list) {
        cout << "list " << listName.Data() << " not contained in file " << fileName.Data() << "!" << endl;
        return NULL;
    }
    list->SetOwner(kTRUE);

    return list;
}

//************************** Function to get spectrum from cocktail input list **************************************
TH1F* GetCocktailInputSpectrum(TList* list, Int_t particle, TString name) {

    if (!list) {
        cout << "no cocktail input list given!" << endl;
        return NULL;
    }

    // read cocktail input spectra
    TH1F*               histo                                   = NULL;
    TObject*            tempObject                              = NULL;
    TGraphErrors*       tempGraphErrs                           = NULL;
    TGraphAsymmErrors*  tempGraphAsymmErrs                      = NULL;
    TString             tempObjectName                          = Form("%sStat", cocktailInputParticles[particle].Data());
    TString             tempObjectClassName                     = "";
    if ( list->FindObject(tempObjectName.Data()) ) {

        tempObject                                              = (TObject*)list->FindObject(tempObjectName.Data());
        tempObjectClassName                                     = (TString)tempObject->ClassName();

        if (tempObjectClassName.Contains("TH1")) {
            histo                                               = (TH1F*)list->FindObject(tempObjectName.Data());
        } else if (tempObjectClassName.Contains("TGraphErrors")) {
            tempGraphErrs                                       = (TGraphErrors*)list->FindObject(tempObjectName.Data());
            histo                                               = TransformGraphToTH1F(tempGraphErrs);
        } else if (tempObjectClassName.Contains("TGraphAsymmErrors")) {
            tempGraphAsymmErrs                                  = (TGraphAsymmErrors*)list->FindObject(tempObjectName.Data());
            histo                                               = TransformGraphToTH1F(tempGraphAsymmErrs);
        } else {
            cout << tempObjectName.Data() << " object type (" << tempObjectClassName.Data() << ") not recognized!" << endl;
            histo                                               = NULL;
        }
    }

    if (histo) histo->SetName(name);
    return histo;
}

//************************** Function to transform graph to TH1D ****************************************************
TH1F* TransformGraphToTH1F(TGraphErrors* graph) {

    if (!graph) return NULL;

    Double_t*   xValue                                          = graph->GetX();
    Double_t*   yValue                                          = graph->GetY();
    Double_t*   xError                                          = graph->GetEX();
    Double_t*   yError                                          = graph->GetEY();
    Int_t       nPoints                                         = graph->GetN();

    Double_t*   newBinningX                                     = new Double_t[nPoints+1];
    for(Int_t i = 0; i < nPoints; i++)  newBinningX[i]          = xValue[i] - xError[i];
    newBinningX[nPoints]                                        = xValue[nPoints-1] + xError[nPoints-1];

    TH1F* histo                                                 = new TH1F("histo","",nPoints,newBinningX);
    for(Int_t i = 1; i <= nPoints; i++){
        histo->SetBinContent(i,   yValue[i-1]);
        histo->SetBinError(i,     yError[i-1]);
    }
    return histo;

    delete[] newBinningX;
}

TH1F* TransformGraphToTH1F(TGraphAsymmErrors* graph) {

    if (!graph) return                                          NULL;

    Double_t*   xValue                                          = graph->GetX();
    Double_t*   yValue                                          = graph->GetY();
    Double_t*   xErrorLow                                       = graph->GetEXlow();
    Double_t*   xErrorHigh                                      = graph->GetEXhigh();
    Double_t*   yErrorLow                                       = graph->GetEYlow();
    Double_t*   yErrorHigh                                      = graph->GetEYhigh();
    Int_t       nPoints                                         = graph->GetN();

    Double_t*   newBinningX                                     = new Double_t[nPoints+1];
    for(Int_t i = 0; i < nPoints; i++) newBinningX[i]           = xValue[i] - xErrorLow[i];
    newBinningX[nPoints]                                        = xValue[nPoints-1] + xErrorHigh[nPoints-1];

    TH1F* histo                                                 = new TH1F("histo","",nPoints,newBinningX);
    for(Int_t i = 1; i <= nPoints; i++){
        histo->SetBinContent(i,   yValue[i-1]);
        histo->SetBinError(i,     (yErrorLow[i-1] + yErrorHigh[i-1]) / 2);
    }

    return histo;
    delete[] newBinningX;
}
