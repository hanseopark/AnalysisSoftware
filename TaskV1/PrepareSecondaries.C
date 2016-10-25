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

void PrepareSecondaries(    TString nameFileCocktail        = "",
                            TString suffix                  = "eps",
                            TString cutSelection            = "",
                            TString option                  = "",
                            TString directphotonPlots       = "",
                            Double_t rapidity               = 0.85,
                            TString period                  = "",
                            Int_t numberOfBins              = 30,
                            Int_t mode                      = 0,
                            Bool_t producePlotsInOrPtRange  = kFALSE
                     ) {
    
    gROOT->Reset();
    
    //************************** Set general style settings *********************************************************
    StyleSettingsThesis();
    SetPlotStyle();
    
    //************************** Set output directory ***************************************************************
    TString outputDir                                           = Form("%s/%s/%s/PrepareSecondaryPi0",cutSelection.Data(),option.Data(),suffix.Data());
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
    TDirectoryFile* topDirCocktail                              = (TDirectoryFile*)fileCocktail.Get("HadronicCocktailMC");
    if (!topDirCocktail) {
        cout << "ERROR: TopDirCocktail not found!" << endl;
        return;
    }
    TList* histoListCocktail                                    = (TList*)topDirCocktail->Get(Form("HadronicCocktailMC_pi0_%.2f", rapidity));
    cout << "searching for " << Form("HadronicCocktailMC_pi0_%.2f", rapidity) << endl;
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
    histoDecayChannels                                          = new TH1F*[nMotherParticles];
    histoPi0PtY                                                 = new TH2F*[nMotherParticles];
    histoPi0PtPhi                                               = new TH2F*[nMotherParticles];
    histoPi0MotherPtY                                           = new TH2F*[nMotherParticles];
    histoPi0MotherPtPhi                                         = new TH2F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (hasMother[i]) {
            cout << "searching for mother[" << i << "] = " << motherParticles[i].Data() << endl;
            
            histoDecayChannels[i]                               = (TH1F*)histoListCocktail->FindObject(Form("DecayChannels_%s", motherParticles[i].Data()));
            
            histoPi0PtY[i]                                      = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_Pi0_From_%s", motherParticles[i].Data()));
            histoPi0PtY[i]->SetName(Form("Pi0_From_%s_Pt_Y_OrBin", motherParticles[i].Data()));
            histoPi0PtY[i]->Sumw2();
            
            histoPi0PtPhi[i]                                    = (TH2F*)histoListCocktail->FindObject(Form("Pt_Phi_Pi0_From_%s", motherParticles[i].Data()));
            histoPi0PtPhi[i]->SetName(Form("Pi0_From_%s_Pt_Phi_OrBin", motherParticles[i].Data()));
            histoPi0PtPhi[i]->Sumw2();
            
            histoPi0MotherPtY[i]                                = (TH2F*)histoListCocktail->FindObject(Form("Pt_Y_%s", motherParticles[i].Data()));
            histoPi0MotherPtY[i]->SetName(Form("%s_Pt_Y_OrBin", motherParticles[i].Data()));
            histoPi0MotherPtY[i]->Sumw2();
            
            histoPi0MotherPtPhi[i]                              = (TH2F*)histoListCocktail->FindObject(Form("Pt_Phi_%s", motherParticles[i].Data()));
            histoPi0MotherPtPhi[i]->SetName(Form("%s_Pt_Phi_OrBin", motherParticles[i].Data()));
            histoPi0MotherPtPhi[i]->Sumw2();
        } else {
            histoDecayChannels[i]                               = NULL;
            histoPi0PtY[i]                                      = NULL;
            histoPi0PtPhi[i]                                    = NULL;
            histoPi0MotherPtY[i]                                = NULL;
            histoPi0MotherPtPhi[i]                              = NULL;
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
    histoPi0PtOrBin                                             = new TH1F*[nMotherParticles];
    histoPi0YOrBin                                              = new TH1F*[nMotherParticles];
    histoPi0PhiOrBin                                            = new TH1F*[nMotherParticles];
    histoPi0MotherPtOrBin                                       = new TH1F*[nMotherParticles];
    histoPi0MotherYOrBin                                        = new TH1F*[nMotherParticles];
    histoPi0MotherPhiOrBin                                      = new TH1F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoPi0PtY[i]) {
            histoPi0PtOrBin[i]                                  = (TH1F*)histoPi0PtY[i]->ProjectionX(Form("Pi0_From_%s_Pt_OrBin", motherParticles[i].Data()), 1, histoPi0PtY[i]->GetNbinsY(), "e");
            SetHistogramTitles(histoPi0PtOrBin[i],"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoPi0PtOrBin[i]->Sumw2();
            histoPi0PtOrBin[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
//             histoPi0PtOrBin[i]->Scale(deltaPhi);
            histoPi0YOrBin[i]                                   = (TH1F*)histoPi0PtY[i]->ProjectionY(Form("Pi0_From_%s_Y_OrBin", motherParticles[i].Data()), 1, histoPi0PtY[i]->GetNbinsX(), "e");
            SetHistogramTitles(histoPi0YOrBin[i],"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoPi0YOrBin[i]->Sumw2();
//             histoPi0YOrBin[i]->Scale(deltaPhi);
        } else {
            histoPi0PtOrBin[i]                                  = NULL;
            histoPi0YOrBin[i]                                   = NULL;
        }
        if (histoPi0PtPhi[i]) {
            histoPi0PhiOrBin[i]                                 = (TH1F*)histoPi0PtPhi[i]->ProjectionY(Form("Pi0_From_%s_Phi_OrBin", motherParticles[i].Data()), 1, histoPi0PtPhi[i]->GetNbinsX(), "e");
            SetHistogramTitles(histoPi0PhiOrBin[i],"","#phi","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoPi0PhiOrBin[i]->Sumw2();
//             histoPi0PhiOrBin[i]->Scale(deltaRap);
        } else
            histoPi0PhiOrBin[i]                                 = NULL;
        
        if (histoPi0MotherPtY[i]) {
            histoPi0MotherPtOrBin[i]                            = (TH1F*)histoPi0MotherPtY[i]->ProjectionX(Form("%s_Pt_OrBin", motherParticles[i].Data()), 1, histoPi0MotherPtY[i]->GetNbinsY(), "e");
            SetHistogramTitles(histoPi0MotherPtOrBin[i],"","#it{p}_{T} (GeV/#it{c})","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoPi0MotherPtOrBin[i]->Sumw2();
            histoPi0MotherPtOrBin[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
            histoPi0MotherPtOrBin[i]->Scale(deltaPhi);
            histoPi0MotherYOrBin[i]                             = (TH1F*)histoPi0MotherPtY[i]->ProjectionY(Form("%s_Y_OrBin", motherParticles[i].Data()), 1, histoPi0MotherPtY[i]->GetNbinsX(), "e");
            SetHistogramTitles(histoPi0MotherYOrBin[i],"","y","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoPi0MotherYOrBin[i]->Sumw2();
//             histoPi0MotherYOrBin[i]->Scale(deltaPhi);
        } else {
            histoPi0MotherPtOrBin[i]                            = NULL;
            histoPi0MotherYOrBin[i]                             = NULL;
        }
        if (histoPi0MotherPtPhi[i]) {
            histoPi0MotherPhiOrBin[i]                           = (TH1F*)histoPi0MotherPtPhi[i]->ProjectionY(Form("%s_Phi_OrBin", motherParticles[i].Data()), 1, histoPi0MotherPtPhi[i]->GetNbinsX(), "e");
            SetHistogramTitles(histoPi0MotherPhiOrBin[i],"","#phi","#frac{1}{N_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}dy} ((GeV/#it{c})^{-1})");
            histoPi0MotherPhiOrBin[i]->Sumw2();
//             histoPi0MotherPhiOrBin[i]->Scale(deltaRap);
        } else
            histoPi0MotherPhiOrBin[i]                           = NULL;
    }

    if (!histoPi0MotherPtOrBin[1]) {
        cout << "ERROR: Didn't get K0s pt spectrum, returning!" << endl;
        return;
    }
    
    //***************************** Scale spectra *******************************************************************
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoPi0PtOrBin[i])           histoPi0PtOrBin[i]->Scale(        1./nEvents);
        if (histoPi0YOrBin[i])            histoPi0YOrBin[i]->Scale(         1./nEvents);
        if (histoPi0PhiOrBin[i])          histoPi0PhiOrBin[i]->Scale(       1./nEvents);
        if (histoPi0MotherPtOrBin[i])     histoPi0MotherPtOrBin[i]->Scale(  1./nEvents);
        if (histoPi0MotherYOrBin[i])      histoPi0MotherYOrBin[i]->Scale(   1./nEvents);
        if (histoPi0MotherPhiOrBin[i])    histoPi0MotherPhiOrBin[i]->Scale( 1./nEvents);
    }

    //***************************** Rebin pt spectra ****************************************************************
    histoPi0Pt                                                  = new TH1F*[nMotherParticles];
    histoPi0MotherPt                                            = new TH1F*[nMotherParticles];
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoPi0PtOrBin[i]) {
            histoPi0Pt[i]                                       = (TH1F*)histoPi0PtOrBin[i]->Clone(Form("Pi0_From_%s_Pt", motherParticles[i].Data()));
            histoPi0Pt[i]->Sumw2();
            histoPi0Pt[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
            RebinSpectrum(histoPi0Pt[i],"");
        } else {
            histoPi0Pt[i]                                       = NULL;
        }
        if (histoPi0MotherPtOrBin[i]) {
            histoPi0MotherPt[i]                                 = (TH1F*)histoPi0MotherPtOrBin[i]->Clone(Form("%s_Pt", motherParticles[i].Data()));
            histoPi0MotherPt[i]->Sumw2();
            histoPi0MotherPt[i]->GetXaxis()->SetRangeUser(ptMin, ptMax);
            RebinSpectrum(histoPi0MotherPt[i],"");
        } else {
            histoPi0MotherPt[i]                                 = NULL;
        }
    }
    
    //***************************** Transform yields ****************************************************************
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoPi0Pt[i])        histoPi0Pt[i]                 = ConvertYieldHisto(histoPi0Pt[i]);
        if (histoPi0MotherPt[i])  histoPi0MotherPt[i]           = ConvertYieldHisto(histoPi0MotherPt[i]);
    }

    
    // adapt plotting range for original binned histograms
    if (!producePlotsInOrPtRange) {
        ptGenMin = ptMin;
        ptGenMax = ptMax;
    }
    //
    // plotting
    //
    
    
    //***************************** Save histograms *****************************************************************
    SaveHistos();
    CreateBRTableLatex();
    
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
    delete[] histoPi0PtPhi;
    delete[] histoPi0PtOrBin;
    delete[] histoPi0Pt;
    delete[] histoPi0YOrBin;
    delete[] histoPi0PhiOrBin;
    delete[] histoPi0MotherPtY;
    delete[] histoPi0MotherPtOrBin;
    delete[] histoPi0MotherPt;
    delete[] histoPi0MotherYOrBin;
    delete[] histoPi0MotherPtPhi;
    delete[] histoPi0MotherPhiOrBin;
    delete[] cocktailInputParametrizations;
    delete[] cocktailInputParametrizationsMtScaled;
}

//************************** Routine for saving histograms **********************************************************
void SaveHistos() {
    TString nameOutput                  = Form("%s/%s/HadronicCocktail%s_%.2f_%s.root", fCutSelection.Data(), fEnergyFlag.Data(), fPeriodFlag.Data(), fRapidity, fCutSelection.Data());
    cout << "INFO: writing into: " << nameOutput << endl;
    TFile *outputFile                   = new TFile(nameOutput,"UPDATE");
    
    // write number of events
    histoNEvents->Write("NEvents", TObject::kOverwrite);
    
    // write binning histogram
    fDeltaPt->Write("deltaPt", TObject::kOverwrite);
    
    // write projections
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoPi0PtOrBin[i])        histoPi0PtOrBin[i]->Write(        histoPi0PtOrBin[i]->GetName(),        TObject::kOverwrite);
        if (histoPi0YOrBin[i])         histoPi0YOrBin[i]->Write(         histoPi0YOrBin[i]->GetName(),         TObject::kOverwrite);
        if (histoPi0PhiOrBin[i])       histoPi0PhiOrBin[i]->Write(       histoPi0PhiOrBin[i]->GetName(),       TObject::kOverwrite);
        if (histoPi0MotherPtOrBin[i])  histoPi0MotherPtOrBin[i]->Write(  histoPi0MotherPtOrBin[i]->GetName(),  TObject::kOverwrite);
        if (histoPi0MotherYOrBin[i])   histoPi0MotherYOrBin[i]->Write(   histoPi0MotherYOrBin[i]->GetName(),   TObject::kOverwrite);
        if (histoPi0MotherPhiOrBin[i]) histoPi0MotherPhiOrBin[i]->Write( histoPi0MotherPhiOrBin[i]->GetName(), TObject::kOverwrite);
    }
    
    // write rebinned histograms
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (histoPi0Pt[i])       histoPi0Pt[i]->Write(histoPi0Pt[i]->GetName(),             TObject::kOverwrite);
        if (histoPi0MotherPt[i]) histoPi0MotherPt[i]->Write(histoPi0MotherPt[i]->GetName(), TObject::kOverwrite);
    }
    
    // write input parametrizations and mt scaled ones
    for (Int_t i=0; i<nMotherParticles; i++) {
        if (cocktailInputParametrizations[i])         cocktailInputParametrizations[i]->Write(          cocktailInputParametrizations[i]->GetName(),         TObject::kOverwrite);
        if (cocktailInputParametrizationsMtScaled[i]) cocktailInputParametrizationsMtScaled[i]->Write(  cocktailInputParametrizationsMtScaled[i]->GetName(), TObject::kOverwrite);
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