// provided by Gamma Conversion Group, $ALICE_ROOT/PWGGA/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion
// ***************************************************************************************************************
// **************** this macro has been written by                                                          ******
// ***************** Friederike Bock, fbock@cern.ch  on the basis of work by Martin Wilde                   ******
// ***************** as of 29th Oct 2015 it superceeds CorrectGamma.C ********************************************
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
#include "TH1F.h"
#include "TH1D.h"
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
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
//#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"

void CorrectGammaEffiResol(TH1D* histoGammaCorr, TH1D *histoAllSec, TH1D *histoK0sSec, TH1D* histoPurity, TH1D* histoConvProb, TH1D* histoRecoEff, Double_t deltaEta, Double_t scaling, Double_t nEvt){
    histoGammaCorr->Scale(1./nEvt);
    histoGammaCorr->Add(histoAllSec,-1);
    histoGammaCorr->Add(histoK0sSec,-1);
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
    histoGammaCorr->Divide(histoGammaCorr,histoConvProb,1.,1.,"");
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    histoGammaCorr->Scale(1./deltaEta);
    histoGammaCorr->Scale(scaling);
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoGammaCorr->GetBinContent(i)/histoGammaCorr->GetBinCenter(i);
        Double_t newBinError = histoGammaCorr->GetBinError(i)/histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,newBinError);
    }
}

void CorrectGammaEffiResolCocktail(TH1D* histoGammaCorr, TH1D *histoK0sSec, TH1D *histoK0lSec, TH1D *histoLambdaSec, TH1D *histoRestSec, TH1D* histoPurity, TH1D* histoConvProb, TH1D* histoRecoEff, Double_t deltaEta, Double_t scaling, Double_t nEvt){
    histoGammaCorr->Scale(1./nEvt);
    histoGammaCorr->Add(histoK0sSec,-1);
    histoGammaCorr->Add(histoK0lSec,-1);
    histoGammaCorr->Add(histoLambdaSec,-1);
    histoGammaCorr->Add(histoRestSec,-1);
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
    histoGammaCorr->Divide(histoGammaCorr,histoConvProb,1.,1.,"");
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    histoGammaCorr->Scale(1./deltaEta);
    histoGammaCorr->Scale(scaling);
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoGammaCorr->GetBinContent(i)/histoGammaCorr->GetBinCenter(i);
        Double_t newBinError = histoGammaCorr->GetBinError(i)/histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,newBinError);
    }
}

void CorrectGammaEffiResol(TH1D* histoGammaCorr, TH1D *histoAllSec, TH1D *histoK0sSec, TH1D* histoPurity, TH1D* histoRecoEff, Double_t deltaEta, Double_t scaling, Double_t nEvt){
    histoGammaCorr->Scale(1./nEvt);
    histoGammaCorr->Add(histoAllSec,-1);
    histoGammaCorr->Add(histoK0sSec,-1);
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    histoGammaCorr->Scale(1./deltaEta);
    histoGammaCorr->Scale(scaling);
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoGammaCorr->GetBinContent(i)/histoGammaCorr->GetBinCenter(i);
        Double_t newBinError = histoGammaCorr->GetBinError(i)/histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,newBinError);
    }
}

void CorrectGammaEffiResolCocktail(TH1D* histoGammaCorr, TH1D *histoK0sSec, TH1D *histoK0lSec, TH1D *histoLambdaSec, TH1D *histoRestSec, TH1D* histoPurity, TH1D* histoRecoEff, Double_t deltaEta, Double_t scaling, Double_t nEvt){
    histoGammaCorr->Scale(1./nEvt);
    histoGammaCorr->Add(histoK0sSec,-1);
    histoGammaCorr->Add(histoK0lSec,-1);
    histoGammaCorr->Add(histoLambdaSec,-1);
    histoGammaCorr->Add(histoRestSec,-1);
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    histoGammaCorr->Scale(1./deltaEta);
    histoGammaCorr->Scale(scaling);
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoGammaCorr->GetBinContent(i)/histoGammaCorr->GetBinCenter(i);
        Double_t newBinError = histoGammaCorr->GetBinError(i)/histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,newBinError);
    }
}

void CorrectGammaSecAndPurity(TH1D* histoGammaCorr, TH1D* histoSecGamma, TH1D* histoSecGammaAddK0s, TH1D* histoPurity){
    histoGammaCorr->Sumw2();
    histoGammaCorr->Add(histoSecGamma,-1);
    histoGammaCorr->Add(histoSecGammaAddK0s,-1);
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
}

void CorrectGammaSecAndPurityCocktail(TH1D* histoGammaCorr, TH1D* histoSecGammaK0s,TH1D* histoSecGammaK0l,TH1D* histoSecGammaLambda, TH1D* histoSecGammaRest, TH1D* histoPurity){
    histoGammaCorr->Sumw2();
    histoGammaCorr->Add(histoSecGammaK0s,-1);
    histoGammaCorr->Add(histoSecGammaK0l,-1);
    histoGammaCorr->Add(histoSecGammaLambda,-1);
    histoGammaCorr->Add(histoSecGammaRest,-1);
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
}

void CorrectGammaUnfoldResol(TH1D* histoGammaCorr, TH1D* histoConvProb, TH1D* histoRecoEff, Double_t deltaEta, Double_t scaling, Double_t nEvt){
    histoGammaCorr->Divide(histoGammaCorr,histoConvProb,1.,1.,"");
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    histoGammaCorr->Scale(1./deltaEta);
    histoGammaCorr->Scale(scaling);
    histoGammaCorr->Scale(1./nEvt);
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoGammaCorr->GetBinContent(i)/histoGammaCorr->GetBinCenter(i);
        Double_t newBinError = histoGammaCorr->GetBinError(i)/histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,newBinError);
    }
}

void CorrectGammaUnfoldResol(TH1D* histoGammaCorr, TH1D* histoRecoEff, Double_t deltaEta, Double_t scaling, Double_t nEvt){
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    histoGammaCorr->Scale(1./deltaEta);
    histoGammaCorr->Scale(scaling);
    histoGammaCorr->Scale(1./nEvt);
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoGammaCorr->GetBinContent(i)/histoGammaCorr->GetBinCenter(i);
        Double_t newBinError = histoGammaCorr->GetBinError(i)/histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,newBinError);
    }
}

Bool_t ConvertCocktailSecondaryToRaw(TH1D* histoGammaSec, TH1D* histoConvProb, TH1D* histoRecoEff, TH2D* responseMatrix, Double_t nEvt, Bool_t useResponseMatrix, Int_t nIterationsUnfolding = 5){
    histoGammaSec->Scale(nEvt);
    
    if (!histoConvProb || !histoRecoEff) return kFALSE;
    histoGammaSec->Multiply(histoConvProb);
    histoGammaSec->Multiply(histoRecoEff);
    
    if (useResponseMatrix) {
        if (!responseMatrix) return kFALSE;
        cout << "will use unfolding for " << histoGammaSec->GetName() << endl;
        RooUnfoldResponse response(0,0,responseMatrix);
        RooUnfoldBayes unfold_SpectrumCocktail (&response, histoGammaSec, nIterationsUnfolding);
        histoGammaSec = (TH1D*)unfold_SpectrumCocktail.Hreco();
    }
    
    return kTRUE;
}

Bool_t ConvertCocktailSecondaryToRaw(TH1D* histoGammaSec, TH1D* histoRecoEff, TH2D* responseMatrix, Double_t nEvt, Bool_t useResponseMatrix, Int_t nIterationsUnfolding = 5){
    histoGammaSec->Scale(nEvt);
    
    if (!histoRecoEff) return kFALSE;
    histoGammaSec->Multiply(histoRecoEff);
    
    if (useResponseMatrix) {
        if (!responseMatrix) return kFALSE;
        cout << "will use unfolding for " << histoGammaSec->GetName() << endl;
        RooUnfoldResponse response(0,0,responseMatrix);
        RooUnfoldBayes unfold_SpectrumCocktail (&response, histoGammaSec, nIterationsUnfolding);
        histoGammaSec = (TH1D*)unfold_SpectrumCocktail.Hreco();
    }

    return kTRUE;
}

void CorrectGammaMC(TH1D* histoMCGammaSpec_MCPtCorr,  Double_t deltaEta, Double_t scaling, Double_t nEvtMC){
    histoMCGammaSpec_MCPtCorr->Scale(1./deltaEta);
    histoMCGammaSpec_MCPtCorr->Scale(scaling);
    histoMCGammaSpec_MCPtCorr->Scale(1./nEvtMC);
    for (Int_t i = 1; i < histoMCGammaSpec_MCPtCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoMCGammaSpec_MCPtCorr->GetBinContent(i)/histoMCGammaSpec_MCPtCorr->GetBinCenter(i);
        Double_t newBinError = histoMCGammaSpec_MCPtCorr->GetBinError(i)/histoMCGammaSpec_MCPtCorr->GetBinCenter(i);
        histoMCGammaSpec_MCPtCorr->SetBinContent(i,newBinContent);
        histoMCGammaSpec_MCPtCorr->SetBinError(i,newBinError);
    }
}


//*********************************************************************************************************
//***************************** Main routine to correct inclusive photons *********************************
//*********************************************************************************************************
void  CorrectGammaV2(   const char *nameUnCorrectedFile     = "myOutput", 
                        const char *nameCorrectionFile      = "", 
                        TString cutSelection                ="", 
                        TString suffix                      = "eps", 
                        TString nameMeson                   = "Pi0", 
                        TString isMC                        = "kFALSE",
                        TString option                      = "", 
                        TString optionPeriod                = "", 
                        TString fEstimatePileup             = "",
                        Int_t mode                          = 9
                    ){
    
    gROOT->Reset();
    
    Int_t   nIterationsUnfolding                = 4;
    Bool_t  useUnfoldingForCocktailSecondaries  = kFALSE;
    
    //****************************************************************************************
    //*************************** Catch old versions *****************************************
    //****************************************************************************************
    if (mode == 1 || mode > 5){
        cout << "This macro isn't designed to run for Dalitz or the old software version" << endl;
        return;
    } else if ( mode == 4 || mode == 5){
        cout << "WARNING: this macro is still under construction for this mode" << endl;
    } else if ( mode == 2 || mode == 3){
        cout << "WARNING: running hybrid mode, this macro is still under construction for this mode" << endl;    
    }
    //*****************************************************************************************
    //************************** Set general style settings ***********************************
    //*****************************************************************************************
    StyleSettingsThesis();  
    SetPlotStyle();

    //*****************************************************************************************
    //*************************** Determine K0s scaling factor ********************************
    //*****************************************************************************************
    Double_t doubleAddFactorK0s     = ReturnCorrectK0ScalingFactor( option,  cutSelection);
    if(isMC.CompareTo("kTRUE") == 0){
        doubleAddFactorK0s = 0.;
        cout<<" MC --> doubleAddFactorK0s = 0"<<endl;
    } else {
        cout << "The additional K0 correction factor is: "  << doubleAddFactorK0s<<endl;
    }

    //*****************************************************************************************
    //*************************** Find out cutSelections for different things *****************
    //*****************************************************************************************    
    TString fEventCutSelection      = "";
    TString fGammaCutSelection      = "";
    TString fClusterCutSelection    = "";
    TString fElectronCutSelection   = "";
    TString fMesonCutSelection      = "";
    TString fCutSelection           = cutSelection;
    ReturnSeparatedCutNumberAdvanced(cutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);

    Bool_t isPCM                    = 0;
    Bool_t isCalo                   = 0;
    if ( mode == 0 || mode == 2 || mode == 3)
        isPCM                       = 1;
    if ( mode == 2 || mode == 3 || mode == 4 || mode == 5)
        isCalo                      = 1;
 
    
    //*****************************************************************************************
    //*************************** Determine centrality from CutString *************************
    //*****************************************************************************************        
    TString centrality              = GetCentralityString(fEventCutSelection);
    TString collisionSystem         = ReturnFullCollisionsSystem(option);
    TString cent                    = "";
    TString textMeasurement         = "#gamma";
    TString detectionProcess        = ReturnFullTextReconstructionProcess(mode);
    TString detectionProcess1       = "";
    TString detectionProcess2       = "";
    if (isPCM && isCalo){
        detectionProcess1           = ReturnFullTextReconstructionProcess(mode,1);
        detectionProcess            = detectionProcess1;
        detectionProcess2           = ReturnFullTextReconstructionProcess(mode,2);
    }
    if(option.Contains("PbPb")){
        cent                        = Form("%s %s", centrality.Data(), collisionSystem.Data());
    } else {
        cent                        = collisionSystem;
    }    
    
    //******************************************************************************************
    //***************************** Set output directory ***************************************
    //******************************************************************************************
    TString outputDir               = Form("%s/%s/%s/CorrectGammaV2",cutSelection.Data(),option.Data(),suffix.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    //******************************************************************************************
    //************************ Define additional global variables ******************************
    //******************************************************************************************
    TString textPi0New              = Form("Gamma_%s",nameMeson.Data());
    
    TString textPrefix2;
    if (isMC.CompareTo("kTRUE")==0){
        textPrefix2                 = "MC";
    } else {
        textPrefix2                 = "data";
    }
    
    Double_t deltaEta               = ReturnDeltaEta(fGammaCutSelection);
    Double_t eta                    = deltaEta*0.5;
    Double_t deltaEtaCalo           = 0;
    Double_t deltaPhiCalo           = 0;
    Double_t minPhiCalo             = 0;
    Double_t maxPhiCalo             = 0;
    Double_t etaCalo                = 0;
    if (isCalo){
        deltaEtaCalo                = ReturnDeltaEtaCalo(fClusterCutSelection);
        etaCalo                     = deltaEtaCalo*0.5;
        deltaPhiCalo                = ReturnDeltaPhiCalo(fClusterCutSelection);
        TString phiMinCut(fClusterCutSelection(GetClusterPhiMinCutPosition(fClusterCutSelection),1));
        TString phiMaxCut(fClusterCutSelection(GetClusterPhiMaxCutPosition(fClusterCutSelection),1));
        minPhiCalo                  = AnalyseClusterMinPhiCut(phiMinCut.Atoi());
        maxPhiCalo                  = AnalyseClusterMaxPhiCut(phiMaxCut.Atoi());
    }
    
    Double_t scaling                = 1./(2.*TMath::Pi());
    Double_t scalingCalo            = 0.;
    if (isCalo) scalingCalo         = 1./deltaPhiCalo;
    
    Bool_t kDoPileup                = kFALSE;
    if (fEstimatePileup.CompareTo("EstimateTrainPileUp") == 0) {
        kDoPileup                   = kTRUE;
        if (mode == 4 || mode == 5) {
            cout << "requested out-of-bunch pileup correction for calo mode, skipping" << endl;
            kDoPileup               = kFALSE;
        }
    }
    
    // Read True Combinatorial Background
    TString combinatorics[17]                               = { "Elec+Elec","Elec+Pion","Elec+Kaon","Elec+Proton","Elec+Muon","Pion+Pion","Pion+Kaon","Pion+Proton",
                                                                "Pion+Muon","Kaon+Kaon","Kaon+Proton","Kaon+Muon","Proton+Proton","Proton+Muon","Muon+Muon","Rest","All"};
    TString combinatoricsCalo[11]                           = { "Electron","Pion","Proton","Kaon","Neutron","K0s","Lambda","Muon","K0l","Rest","All"};
    
    Color_t colorsCombinatorics[17]                         = { kAzure-1, kRed+1, kOrange+7, kMagenta+2, kRed-9,
                                                                kBlue-6, kAzure+5, kPink, kCyan-3, kGreen-3,
                                                                kSpring+9, kGreen+2, kBlue+2, kMagenta-6, kSpring+4,
                                                                kCyan+2, 809};
    Color_t colorsCombinatoricsCalo[11]                     = { kAzure-1, kRed+1, kOrange+7, kMagenta+2, kRed-9, kBlue,
                                                                kBlue-6, kAzure+5, kPink, kCyan-3, kGreen-3};
    
    Style_t markersCombinatorics[17]                        = { 20, 21, 24, 25, 27,
                                                                28, 29, 30, 33, 34, 
                                                                20, 21, 24, 25, 27,
                                                                28, 29};
    Style_t markersCombinatoricsCalo[11]                    = { 20, 21, 24, 25, 27,
                                                                28, 29, 30, 33, 34, 20};
    
    TString decays[9]                                       = { "Pi0","Eta","Etap","Omega","Rho",
                                                                "Phi","Sigma", "All decays","direct #gamma"
                                                              };
    TString decaysLabels[9]                                  = { "#gamma from #pi^{0}","#gamma from #eta","#gamma from #eta'","#gamma from #omega","#gamma from #rho^{0}",
                                                                "#gamma from #phi", "#gamma from #Sigma^{0}", "All decays","direct #gamma"
                                                              };

    Color_t colorsDecay[9]                                  = { kRed+1, kAzure-1, kOrange+7, kGreen+2, kCyan-3, 809, 
                                                                kBlue+2, kBlack, kRed-1
                                                              };
    
    
    //****************************************************************************************** 
    //*******************File definitions and reading histograms from files ********************
    //******************************************************************************************
    TFile*  fileUnCorrected                                             = new TFile(nameUnCorrectedFile);
    TFile*  fileCorrections                                             = new TFile(nameCorrectionFile);
    if (!fileUnCorrected) {
        cout << "file " << nameUnCorrectedFile << " not found!" << endl;
        return;
    }
    if (!nameCorrectionFile) {
        cout << "file " << nameCorrectionFile << " not found!" << endl;
        return;
    }

    //******************* Conversion gamma spectra *********************************************
    TH1D*   histoESDConvGammaPt                                         = NULL;
    TH1D*   histoESDConvGammaPt_OriginalBin                             = NULL;
    if (isPCM) {
        histoESDConvGammaPt                                             = (TH1D*)fileUnCorrected->Get("ESD_ConvGamma_Pt");
        histoESDConvGammaPt_OriginalBin                                 = (TH1D*)fileUnCorrected->Get("ESD_ConvGamma_Pt_OriginalBinning");
    }

    //******************* Calo gamma spectra ***************************************************
    TH1D*   histoESDCaloGammaPt                                         = NULL;
    TH1D*   histoESDCaloGammaPt_OriginalBin                             = NULL;
    if (isCalo) {
        histoESDCaloGammaPt                                             = (TH1D*)fileUnCorrected->Get("ESD_CaloGamma_Pt");
        histoESDCaloGammaPt_OriginalBin                                 = (TH1D*)fileUnCorrected->Get("ESD_CaloGamma_Pt_OriginalBinning");
    }
    
    //******************* Calculate number of events for normalization *************************
    TH1D*   histoEventQuality                                           = (TH1D*)fileUnCorrected->Get("NEvents");
    Float_t nEvt                                                        = 0.;
    if (option.CompareTo("PbPb_2.76TeV") == 0)  nEvt                    = histoEventQuality->GetBinContent(1);
    else                                        nEvt                    = GetNEvents(histoEventQuality);
    
    TH1F*   histoEventQualityMC                                         = (TH1F*)fileCorrections->Get("NEvents");
    Float_t nEvtMC                                                      = 0.;
    if (option.CompareTo("PbPb_2.76TeV") == 0)  nEvtMC                  = histoEventQualityMC->GetBinContent(1);
    else                                        nEvtMC                  = GetNEvents(histoEventQualityMC);
    
    //******************* Binning histogram ****************************************************
    TH1D*   deltaPt                                                     = (TH1D*)fileUnCorrected->Get("deltaPt");
    for (Int_t i = 0; i < deltaPt->GetNbinsX() +1; i++) deltaPt->SetBinError(i, 0);

    //******************* Primary gamma correction factors (PCM and Calo) **********************
    TH1D*   histoGammaPurity_Pt                                         = NULL;
    TH1D*   histoGammaTruePurity_Pt                                     = NULL;
    TH1D*   histoGammaTruePurity_OriginalBin_Pt                         = NULL;
    TH1D*   histoGammaPrimaryRecoEff_Pt                                 = NULL;
    TH1D*   histoGammaPrimaryRecoEff_MCPt                               = NULL;
    if ( isPCM || isCalo ) {
        // reconstructed validated photons / all reconstructed photon candidates in MC
        histoGammaPurity_Pt                                             = (TH1D*)fileCorrections->Get("GammaPurity_Pt");
        // reconstructed validated primary photons / (all reconstructed photon candidates in MC - validated reconstructed secondary photons)
        // rescaling of histoGammaPurity_Pt for primary photons
        histoGammaTruePurity_Pt                                         = (TH1D*)fileCorrections->Get("GammaTruePurity_Pt");
        // reconstructed validated primary photons / (all reconstructed photon candidates in MC - validated reconstructed secondary photons)
        // rescaling of histoGammaPurity_Pt for primary photons
        // rebinned
        histoGammaTruePurity_OriginalBin_Pt                             = (TH1D*)fileCorrections->Get("GammaTruePurity_OriginalBinning_Pt");
        // reconstructed validated primary photons vs rec pt / converted MC photons vs MCpt
        // this reconstruction efficiency includes the transverse momentum resolution correction and the acceptance correction
        histoGammaPrimaryRecoEff_Pt                                     = (TH1D*)fileCorrections->Get("GammaPrimaryRecoEff_Pt");
        // reconstructed validated primary photons vs MC pt / converted MC photons vs MCpt
        // this reconstruction efficiency includes the acceptance correction
        histoGammaPrimaryRecoEff_MCPt                                   = (TH1D*)fileCorrections->Get("GammaPrimaryRecoEff_MCPt");
    }
    TH1D*   histoGammaConvProb_MCPt                                     = NULL;
    if ( isPCM ) {
        // converted photons MC/ all generated photons (always based on primaries only)
        histoGammaConvProb_MCPt                                         = (TH1D*)fileCorrections->Get("GammaConvProb_MCPt");
    }
    
    //******************* Primary gamma correction factors (PCM-Calo hybrid) *******************
    TH1D*   histoGammaCaloPurity_Pt                                     = NULL;
    TH1D*   histoGammaCaloTruePurity_Pt                                 = NULL;
    TH1D*   histoGammaCaloTruePurity_OrigBin_Pt                         = NULL;
    TH1D*   histoGammaCaloPrimaryRecoEff_Pt                             = NULL;
    TH1D*   histoGammaCaloPrimaryRecoEff_MCPt                           = NULL;
    // Calorimeter photons 
    if ( isPCM && isCalo ) {
        // reconstructed validated photons / all reconstructed photon candidates in MC
        histoGammaCaloPurity_Pt                                         = (TH1D*)fileCorrections->Get("GammaCaloPurity_Pt");
        // reconstructed validated primary photons / (all reconstructed photon candidates in MC - validated reconstructed secondary photons)
        // rescaling of histoGammaPurity_Pt for primary photons
        histoGammaCaloTruePurity_Pt                                     = (TH1D*)fileCorrections->Get("GammaCaloTruePurity_Pt");
        // reconstructed validated primary photons / (all reconstructed photon candidates in MC - validated reconstructed secondary photons)
        // rescaling of histoGammaPurity_Pt for primary photons
        // rebinned
        histoGammaCaloTruePurity_OrigBin_Pt                             = (TH1D*)fileCorrections->Get("GammaCaloTruePurity_OriginalBinning_Pt");
        // reconstructed validated primary photons vs rec pt / converted MC photons vs MCpt 
        // this reconstruction efficiency includes the transverse momentum resolution correction and the acceptance correction    
        histoGammaCaloPrimaryRecoEff_Pt                                 = (TH1D*)fileCorrections->Get("GammaCaloPrimaryRecoEff_Pt");
        // reconstructed validated primary photons vs MC pt / converted MC photons vs MCpt 
        // this reconstruction efficiency includes the acceptance correction        
        histoGammaCaloPrimaryRecoEff_MCPt                               = (TH1D*)fileCorrections->Get("GammaCaloPrimaryRecoEff_MCPt");
    }

    //******************* MC histograms (PCM and Calo) *****************************************
    TH1D*   histoMCrecGamma_Pt                                          = NULL;
    TH1D*   histoMCrecGamma_OriginalBin_Pt                              = NULL;
    TH1D*   histoMCrecPrimaryGamma_Pt                                   = NULL;
    TH1D*   histoMCGammaConv_MCPt                                       = NULL;
    TH1D*   histoGammaTrueConv_Pt                                       = NULL;
    TH1D*   histoGammaTruePrimaryConv_Pt                                = NULL;
    TH2D*   histoGammaTruePrimaryConv_recPt_MCPt                        = NULL;
    if (isPCM) {
        // reconstructed photon candidates in MC
        histoMCrecGamma_Pt                                              = (TH1D*)fileCorrections->Get("MCrec_ConvGamma_Pt");
        histoMCrecGamma_OriginalBin_Pt                                  = (TH1D*)fileCorrections->Get("MCrec_ConvGamma_OriginalBinning_Pt");
        // reconstructed photon candidates in MC - validated secondary photons
        histoMCrecPrimaryGamma_Pt                                       = (TH1D*)fileCorrections->Get("MCrec_PrimaryConvGamma_Pt");
        // MC converted photons regardless of source
        histoMCGammaConv_MCPt                                           = (TH1D*)fileCorrections->Get("MC_ConvGamma_MCPt");
        // validated reconstructed photons
        histoGammaTrueConv_Pt                                           = (TH1D*)fileCorrections->Get("TrueConvGamma_Pt");
        // validated reconstructed primary photons
        histoGammaTruePrimaryConv_Pt                                    = (TH1D*)fileCorrections->Get("TruePrimaryConvGamma_Pt");
        // response matrix for photons recPt vs MC pt
        histoGammaTruePrimaryConv_recPt_MCPt                            = (TH2D*)fileCorrections->Get("TruePrimaryConvGamma_recPt_MCPt");
        // response matrix for photons recPt vs MC pt rebinned
    }
    TH1D*   histoMCrecBackground_Pt                                     = NULL;
    TH1D*   histoMCAllGamma_MCPt                                        = NULL;
    TH1D*   histoMCAllGamma_OriginalBin_MCPt                            = NULL;
    if (isPCM ||isCalo) {
        // reconstructed background in MC
        histoMCrecBackground_Pt                                         = (TH1D*)fileCorrections->Get("MCrec_Background");
        // MC input for all gamma, regardless of source
        histoMCAllGamma_MCPt                                            = (TH1D*)fileCorrections->Get("MC_AllGamma_MCPt");
        histoMCAllGamma_OriginalBin_MCPt                                = (TH1D*)fileCorrections->Get("MC_AllGamma_OriginalBinning_MCPt");
    }

    //******************* MC histograms (Calo) *************************************************
    TH1D*   histoMCrecCaloBackground_Pt                                 = NULL;
    TH1D*   histoMCAllGammaCalo_MCPt                                    = NULL;
    TH1D*   histoMCAllGammaCalo_OriginalBin_MCPt                        = NULL;
    if (isPCM && isCalo) {
        // reconstructed background in MC
        histoMCrecCaloBackground_Pt                                     = (TH1D*)fileCorrections->Get("MCrec_Calo_Background");
        // MC input for all gamma, regardless of source
        histoMCAllGammaCalo_MCPt                                        = (TH1D*)fileCorrections->Get("MC_AllGammaEMCAcc_MCPt");
        histoMCAllGammaCalo_OriginalBin_MCPt                            = (TH1D*)fileCorrections->Get("MC_AllGammaEMCAcc_OriginalBinning_MCPt");
    }
    TH1D*   histoMCrecGammaCalo_Pt                                      = NULL;
    TH1D*   histoMCrecGammaCalo_OriginalBin_Pt                          = NULL;
    TH1D*   histoMCrecPrimaryGammaCalo_Pt                               = NULL;
    TH1D*   histoGammaTrueCalo_Pt                                       = NULL;
    TH1D*   histoGammaTruePrimaryCalo_Pt                                = NULL;
    TH2D*   histoGammaTruePrimaryCalo_recPt_MCPt                        = NULL;
    if (isCalo) {
        // reconstructed photon candidates in MC
        histoMCrecGammaCalo_Pt                                          = (TH1D*)fileCorrections->Get("MCrec_CaloGamma_Pt");
        histoMCrecGammaCalo_OriginalBin_Pt                              = (TH1D*)fileCorrections->Get("MCrec_CaloGamma_OriginalBinning_Pt");
        // reconstructed photon candidates in MC - validated secondary photons
        histoMCrecPrimaryGammaCalo_Pt                                   = (TH1D*)fileCorrections->Get("MCrec_PrimaryCaloGamma_Pt");
        // validated reconstructed photons
        histoGammaTrueCalo_Pt                                           = (TH1D*)fileCorrections->Get("TrueCaloGamma_Pt");
        // validated reconstructed primary photons
        histoGammaTruePrimaryCalo_Pt                                    = (TH1D*)fileCorrections->Get("TruePrimaryCaloGamma_Pt");
        // response matrix for photons recPt vs MC pt (check for binning)
        histoGammaTruePrimaryCalo_recPt_MCPt                            = (TH2D*)fileCorrections->Get("TruePrimaryCaloGamma_recPt_MCPt");
    }

    //******************* MC decay gammas ******************************************************
    TH1D**  histoPhotonSource_MCPt                                      = new TH1D*[9];
    for (Int_t i = 0;i<7;i++) {
        histoPhotonSource_MCPt[i]                                       = (TH1D*)fileCorrections->Get(Form("MC_DecayGamma%s_Pt",decays[i].Data()));
        histoPhotonSource_MCPt[i]->Sumw2();
        if (i == 0) histoPhotonSource_MCPt[7]                           = (TH1D*)histoPhotonSource_MCPt[0]->Clone("MC_DecayGammaAll_Pt");
        else        histoPhotonSource_MCPt[7]->Add(histoPhotonSource_MCPt[i]);
    }
    histoPhotonSource_MCPt[8]                                           = (TH1D*)histoMCAllGamma_OriginalBin_MCPt->Clone("MC_DirectPhotons");
    histoPhotonSource_MCPt[8]->Sumw2();
    histoPhotonSource_MCPt[8]->Add(histoPhotonSource_MCPt[7],-1);

    //******************* MC combinatorial gammas (PCM) ****************************************
    TH1D**  histoCombinatorialSpecies_Pt                                = NULL;
    if (isPCM) {
        histoCombinatorialSpecies_Pt                                    = new TH1D*[17];
        for(Int_t i = 0;i<17;i++){
            histoCombinatorialSpecies_Pt[i]                             = (TH1D*)fileCorrections->Get(Form("ESD_TrueComb%s_Pt",combinatorics[i].Data()));
            histoCombinatorialSpecies_Pt[i]->SetMinimum(1e-10);
        }
    }

    //******************* MC combinatorial gammas (Calo) ***************************************
    TH1D**  histoCombinatorialSpeciesCalo_Pt                            = NULL;
    if (isCalo) {
        histoCombinatorialSpeciesCalo_Pt                                = new TH1D*[11];
        for(Int_t i = 0;i<11;i++){
            histoCombinatorialSpeciesCalo_Pt[i]                         = (TH1D*)fileCorrections->Get(Form("ESD_TrueComb%s_Pt",combinatoricsCalo[i].Data()));
            histoCombinatorialSpeciesCalo_Pt[i]->SetMinimum(1e-10);
        }
    }

    //******************* MC true secondary gammas (PCM) ***************************************
    TH1D*   histoGammaTrueSecConv_Pt                                    = NULL;
    TH1D*   histoGammaTrueSecConvGammaFromXFromK0s_Pt                   = NULL;
    TH1D*   histoGammaTrueSecConvGammaFromXFromK0l_Pt                   = NULL;
    TH1D*   histoGammaTrueSecConvGammaFromXFromLambda_Pt                = NULL;
    if (isPCM) {
        histoGammaTrueSecConv_Pt                                        = (TH1D*)fileCorrections->Get("TrueSecondaryConvGamma_Pt");
        histoGammaTrueSecConvGammaFromXFromK0s_Pt                       = (TH1D*)fileCorrections->Get("TrueSecondaryConvGammaFromXFromK0s_Pt");
        histoGammaTrueSecConvGammaFromXFromK0l_Pt                       = (TH1D*)fileCorrections->Get("TrueSecondaryConvGammaFromXFromK0l_Pt");
        histoGammaTrueSecConvGammaFromXFromLambda_Pt                    = (TH1D*)fileCorrections->Get("TrueSecondaryConvGammaFromXFromLambda_Pt");
    }
    
    //******************* MC true secondary gammas (Calo) **************************************
    TH1D*   histoGammaTrueSecCalo_Pt                                    = NULL;
    TH1D*   histoGammaTrueSecCaloGammaFromXFromK0s_Pt                   = NULL;
    TH1D*   histoGammaTrueSecCaloGammaFromXFromK0l_Pt                   = NULL;
    TH1D*   histoGammaTrueSecCaloGammaFromXFromLambda_Pt                = NULL;
    if (isCalo) {
        histoGammaTrueSecCalo_Pt                                        = (TH1D*)fileCorrections->Get("TrueSecondaryCaloGamma_Pt");
        histoGammaTrueSecCaloGammaFromXFromK0s_Pt                       = (TH1D*)fileCorrections->Get("TrueSecondaryCaloGammaFromXFromK0s_Pt");
        histoGammaTrueSecCaloGammaFromXFromK0l_Pt                       = (TH1D*)fileCorrections->Get("TrueSecondaryCaloGammaFromXFromK0l_Pt");
        histoGammaTrueSecCaloGammaFromXFromLambda_Pt                    = (TH1D*)fileCorrections->Get("TrueSecondaryCaloGammaFromXFromLambda_Pt");
    }
    
    //******************* MC true secondary gamma fractions (PCM or Calo) **********************
    TH1D*   histoFracAllGammaToSec_Pt                                   = NULL;
    TH1D*   histoFracAllGammaToSec_OriginalBin_Pt                       = NULL;
    TH1D*   histoFracAllGammaToSecFromXFromK0s_Pt                       = NULL;
    TH1D*   histoFracAllGammaToSecFromXFromK0s_OriginalBin_Pt           = NULL;
    TH1D*   histoFracAllGammaToSecFromXFromK0l_Pt                       = NULL;
    TH1D*   histoFracAllGammaToSecFromXFromK0l_OriginalBin_Pt           = NULL;
    TH1D*   histoFracAllGammaToSecFromXFromLambda_Pt                    = NULL;
    TH1D*   histoFracAllGammaToSecFromXFromLambda_OriginalBin_Pt        = NULL;
    if ( isPCM || isCalo ) {
        histoFracAllGammaToSec_Pt                                       = (TH1D*)fileCorrections->Get("FracAllGammaToSec");
        histoFracAllGammaToSec_OriginalBin_Pt                           = (TH1D*)fileCorrections->Get("FracAllGammaToSecOriginalBinning");
        histoFracAllGammaToSecFromXFromK0s_Pt                           = (TH1D*)fileCorrections->Get("FracAllGammaToSecFromXFromK0s");
        histoFracAllGammaToSecFromXFromK0s_OriginalBin_Pt               = (TH1D*)fileCorrections->Get("FracAllGammaToSecFromXFromK0sOriginalBinning");
        histoFracAllGammaToSecFromXFromK0l_Pt                           = (TH1D*)fileCorrections->Get("FracAllGammaToSecFromXFromK0l");
        histoFracAllGammaToSecFromXFromK0l_OriginalBin_Pt               = (TH1D*)fileCorrections->Get("FracAllGammaToSecFromXFromK0lOriginalBinning");
        histoFracAllGammaToSecFromXFromLambda_Pt                        = (TH1D*)fileCorrections->Get("FracAllGammaToSecFromXFromLambda");
        histoFracAllGammaToSecFromXFromLambda_OriginalBin_Pt            = (TH1D*)fileCorrections->Get("FracAllGammaToSecFromXFromLambdaOriginalBinning");
    }
    
    //******************* Secondary gamma spectra from cocktail ********************************
    Bool_t  hasCocktailInput                                            = kTRUE;
    TH1D*   histoGammaTrueSecCocktailGammaFromXFromK0s_Pt               = NULL;
    TH1D*   histoGammaTrueSecCocktailGammaFromXFromK0s_PtOrBin          = NULL;
    TH1D*   histoGammaTrueSecCocktailGammaFromXFromK0l_Pt               = NULL;
    TH1D*   histoGammaTrueSecCocktailGammaFromXFromK0l_PtOrBin          = NULL;
    TH1D*   histoGammaTrueSecCocktailGammaFromXFromLambda_Pt            = NULL;
    TH1D*   histoGammaTrueSecCocktailGammaFromXFromLambda_PtOrBin       = NULL;
    TH1D*   histoGammaTrueSecCocktailGammaRest_Pt                       = NULL;
    TH1D*   histoGammaTrueSecCocktailGammaRest_PtOrBin                  = NULL;
    if( isPCM || isCalo ){
        histoGammaTrueSecCocktailGammaFromXFromK0s_Pt                   = (TH1D*)fileUnCorrected->Get("CocktailSecondaryGammaFromXFromK0s_Pt");
        histoGammaTrueSecCocktailGammaFromXFromK0s_PtOrBin              = (TH1D*)fileUnCorrected->Get("CocktailSecondaryGammaFromXFromK0s_PtOrBin");
        histoGammaTrueSecCocktailGammaFromXFromK0l_Pt                   = (TH1D*)fileUnCorrected->Get("CocktailSecondaryGammaFromXFromK0l_Pt");
        histoGammaTrueSecCocktailGammaFromXFromK0l_PtOrBin              = (TH1D*)fileUnCorrected->Get("CocktailSecondaryGammaFromXFromK0l_PtOrBin");
        histoGammaTrueSecCocktailGammaFromXFromLambda_Pt                = (TH1D*)fileUnCorrected->Get("CocktailSecondaryGammaFromXFromLambda_Pt");
        histoGammaTrueSecCocktailGammaFromXFromLambda_PtOrBin           = (TH1D*)fileUnCorrected->Get("CocktailSecondaryGammaFromXFromLambda_PtOrBin");

        if ( isPCM && !isCalo ) {
            histoGammaTrueSecCocktailGammaRest_Pt                       = (TH1D*)fileCorrections->Get("fHistoGammaTrueSecondaryConvGammaRestPt");
            histoGammaTrueSecCocktailGammaRest_PtOrBin                  = (TH1D*)fileCorrections->Get("fHistoGammaTrueSecondaryConvGammaRestPtOrBin");
        }
        if ( isCalo && !isPCM ) {
            histoGammaTrueSecCocktailGammaRest_Pt                       = (TH1D*)fileCorrections->Get("fHistoGammaTrueSecondaryCaloRestPt");
            histoGammaTrueSecCocktailGammaRest_PtOrBin                  = (TH1D*)fileCorrections->Get("fHistoGammaTrueSecondaryCaloRestPtOrBin");
        }
        
        // test for secondary spectra from cocktail
        if (!histoGammaTrueSecCocktailGammaFromXFromK0s_Pt || !histoGammaTrueSecCocktailGammaFromXFromK0l_Pt
            || !histoGammaTrueSecCocktailGammaFromXFromLambda_Pt || !histoGammaTrueSecCocktailGammaRest_Pt)
            hasCocktailInput                                            = kFALSE;
        if (!histoGammaTrueSecCocktailGammaFromXFromK0s_PtOrBin || !histoGammaTrueSecCocktailGammaFromXFromK0l_PtOrBin
            || !histoGammaTrueSecCocktailGammaFromXFromLambda_PtOrBin || !histoGammaTrueSecCocktailGammaRest_PtOrBin)
            hasCocktailInput                                            = kFALSE;
        if (!hasCocktailInput) cout << "secondary spectra from cocktail not found, will not use" << endl;
    }

    //******************* Secondary gamma reco eff *********************************************
    TH1D*   histoGammaSecondaryFromXFromK0sRecoEff_MCPt                 = NULL;
    TH1D*   histoGammaSecondaryFromXFromK0sRecoEff_MCPtOrBin            = NULL;
    TH1D*   histoGammaSecondaryFromXFromK0sRecoEff_Pt                   = NULL;
    TH1D*   histoGammaSecondaryFromXFromK0sRecoEff_PtOrBin              = NULL;
    
    TH1D*   histoGammaSecondaryFromXFromK0lRecoEff_MCPt                 = NULL;
    TH1D*   histoGammaSecondaryFromXFromK0lRecoEff_MCPtOrBin            = NULL;
    TH1D*   histoGammaSecondaryFromXFromK0lRecoEff_Pt                   = NULL;
    TH1D*   histoGammaSecondaryFromXFromK0lRecoEff_PtOrBin              = NULL;
    
    TH1D*   histoGammaSecondaryFromXFromLambdaRecoEff_MCPt              = NULL;
    TH1D*   histoGammaSecondaryFromXFromLambdaRecoEff_MCPtOrBin         = NULL;
    TH1D*   histoGammaSecondaryFromXFromLambdaRecoEff_Pt                = NULL;
    TH1D*   histoGammaSecondaryFromXFromLambdaRecoEff_PtOrBin           = NULL;
    if ( hasCocktailInput && (isPCM || isCalo) ) {
        histoGammaSecondaryFromXFromK0sRecoEff_MCPt                     = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromK0sRecoEff_MCPt");
        histoGammaSecondaryFromXFromK0sRecoEff_MCPtOrBin                = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromK0sRecoEff_MCPtOrBin");
        histoGammaSecondaryFromXFromK0sRecoEff_Pt                       = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromK0sRecoEff_Pt");
        histoGammaSecondaryFromXFromK0sRecoEff_PtOrBin                  = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromK0sRecoEff_PtOrBin");
        
        histoGammaSecondaryFromXFromK0lRecoEff_MCPt                     = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromK0lRecoEff_MCPt");
        histoGammaSecondaryFromXFromK0lRecoEff_MCPtOrBin                = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromK0lRecoEff_MCPtOrBin");
        histoGammaSecondaryFromXFromK0lRecoEff_Pt                       = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromK0lRecoEff_Pt");
        histoGammaSecondaryFromXFromK0lRecoEff_PtOrBin                  = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromK0lRecoEff_PtOrBin");

        histoGammaSecondaryFromXFromLambdaRecoEff_MCPt                  = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromLambdaRecoEff_MCPt");
        histoGammaSecondaryFromXFromLambdaRecoEff_MCPtOrBin             = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromLambdaRecoEff_MCPtOrBin");
        histoGammaSecondaryFromXFromLambdaRecoEff_Pt                    = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromLambdaRecoEff_Pt");
        histoGammaSecondaryFromXFromLambdaRecoEff_PtOrBin               = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromLambdaRecoEff_PtOrBin");
    }
    
    //******************* Secondary gamma conv prob ********************************************
    TH1D*   histoGammaSecondaryFromXFromK0sConvProb_MCPt                = NULL;
    TH1D*   histoGammaSecondaryFromXFromK0sConvProb_MCPtOrBin           = NULL;
    TH1D*   histoGammaSecondaryFromXFromK0lConvProb_MCPt                = NULL;
    TH1D*   histoGammaSecondaryFromXFromK0lConvProb_MCPtOrBin           = NULL;
    TH1D*   histoGammaSecondaryFromXFromLambdaConvProb_MCPt             = NULL;
    TH1D*   histoGammaSecondaryFromXFromLambdaConvProb_MCPtOrBin        = NULL;
    if ( hasCocktailInput && isPCM ) {
        histoGammaSecondaryFromXFromK0sConvProb_MCPt                    = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromK0sConvProb_MCPt");
        histoGammaSecondaryFromXFromK0sConvProb_MCPtOrBin               = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromK0sConvProb_MCPtOrBin");
        histoGammaSecondaryFromXFromK0lConvProb_MCPt                    = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromK0lConvProb_MCPt");
        histoGammaSecondaryFromXFromK0lConvProb_MCPtOrBin               = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromK0lConvProb_MCPtOrBin");
        histoGammaSecondaryFromXFromLambdaConvProb_MCPt                 = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromLambdaConvProb_MCPt");
        histoGammaSecondaryFromXFromLambdaConvProb_MCPtOrBin            = (TH1D*)fileCorrections->Get("SecondaryGammaFromXFromLambdaConvProb_MCPtOrBin");

    }

    //******************* Secondary gamma response matrices ************************************
    TH2D*   histoGammaTrueSecondaryFromXFromK0s_MCPt_recPt              = NULL;
    TH2D*   histoGammaTrueSecondaryFromXFromK0s_MCPt_recPtOrBin         = NULL;
    TH2D*   histoGammaTrueSecondaryFromXFromK0l_MCPt_recPt              = NULL;
    TH2D*   histoGammaTrueSecondaryFromXFromK0l_MCPt_recPtOrBin         = NULL;
    TH2D*   histoGammaTrueSecondaryFromXFromLambda_MCPt_recPt           = NULL;
    TH2D*   histoGammaTrueSecondaryFromXFromLambda_MCPt_recPtOrBin      = NULL;
    if ( hasCocktailInput && isPCM && !isCalo ) {
        // PCM
        histoGammaTrueSecondaryFromXFromK0s_MCPt_recPt                  = (TH2D*)fileCorrections->Get("TrueSecondaryConvGammaFromXFromK0s_MCPt_recPt");
        histoGammaTrueSecondaryFromXFromK0s_MCPt_recPtOrBin             = (TH2D*)fileCorrections->Get("TrueSecondaryConvGammaFromXFromK0s_MCPt_recPt_orBin");
        histoGammaTrueSecondaryFromXFromK0l_MCPt_recPt                  = (TH2D*)fileCorrections->Get("TrueSecondaryConvGammaFromXFromK0l_MCPt_recPt");
        histoGammaTrueSecondaryFromXFromK0l_MCPt_recPtOrBin             = (TH2D*)fileCorrections->Get("TrueSecondaryConvGammaFromXFromK0l_MCPt_recPt_orBin");
        histoGammaTrueSecondaryFromXFromLambda_MCPt_recPt               = (TH2D*)fileCorrections->Get("TrueSecondaryConvGammaFromXFromLambda_MCPt_recPt");
        histoGammaTrueSecondaryFromXFromLambda_MCPt_recPtOrBin          = (TH2D*)fileCorrections->Get("TrueSecondaryConvGammaFromXFromLambda_MCPt_recPt_orBin");
    }
    if ( hasCocktailInput && isCalo && !isPCM ) {
        // Calo
        histoGammaTrueSecondaryFromXFromK0s_MCPt_recPt                  = (TH2D*)fileCorrections->Get("TrueSecondaryCaloGammaFromXFromK0s_MCPt_recPt");
        histoGammaTrueSecondaryFromXFromK0s_MCPt_recPtOrBin             = (TH2D*)fileCorrections->Get("TrueSecondaryCaloGammaFromXFromK0s_MCPt_recPt_orBin");
        histoGammaTrueSecondaryFromXFromK0l_MCPt_recPt                  = (TH2D*)fileCorrections->Get("TrueSecondaryCaloGammaFromXFromK0l_MCPt_recPt");
        histoGammaTrueSecondaryFromXFromK0l_MCPt_recPtOrBin             = (TH2D*)fileCorrections->Get("TrueSecondaryCaloGammaFromXFromK0l_MCPt_recPt_orBin");
        histoGammaTrueSecondaryFromXFromLambda_MCPt_recPt               = (TH2D*)fileCorrections->Get("TrueSecondaryCaloGammaFromXFromLambda_MCPt_recPt");
        histoGammaTrueSecondaryFromXFromLambda_MCPt_recPtOrBin          = (TH2D*)fileCorrections->Get("TrueSecondaryCaloGammaFromXFromLambda_MCPt_recPt_orBin");
    }
    
    //******************* Calculate raw secondary spectra from cocktail input ******************
    TH1D*   histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt              = NULL;
    TH1D*   histoGammaSecGammaFromXFromK0s_Cocktail_Raw_PtOrBin         = NULL;
    TH1D*   histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt              = NULL;
    TH1D*   histoGammaSecGammaFromXFromK0l_Cocktail_Raw_PtOrBin         = NULL;
    TH1D*   histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt           = NULL;
    TH1D*   histoGammaSecGammaFromXFromLambda_Cocktail_Raw_PtOrBin      = NULL;
    if ( hasCocktailInput && (isPCM || isCalo) ) {
        cout << "calculating raw secondary spectra from cocktail" << endl;

        // K0s: clone cocktail spectra
        histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt                  = (TH1D*) histoGammaTrueSecCocktailGammaFromXFromK0s_Pt->Clone("histoGammaTrueSecConvGammaFromXFromK0s_Cocktail_Raw_Pt");
        histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt->Sumw2();
        histoGammaSecGammaFromXFromK0s_Cocktail_Raw_PtOrBin             = (TH1D*) histoGammaTrueSecCocktailGammaFromXFromK0s_PtOrBin->Clone("histoGammaTrueSecConvGammaFromXFromK0s_Cocktail_Raw_PtOrBin");
        histoGammaSecGammaFromXFromK0s_Cocktail_Raw_PtOrBin->Sumw2();
        
        // K0l: clone cocktail spectra
        histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt                  = (TH1D*) histoGammaTrueSecCocktailGammaFromXFromK0l_Pt->Clone("histoGammaTrueSecConvGammaFromXFromK0l_Cocktail_Raw_Pt");
        histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt->Sumw2();
        histoGammaSecGammaFromXFromK0l_Cocktail_Raw_PtOrBin             = (TH1D*) histoGammaTrueSecCocktailGammaFromXFromK0l_PtOrBin->Clone("histoGammaTrueSecConvGammaFromXFromK0l_Cocktail_Raw_PtOrBin");
        histoGammaSecGammaFromXFromK0l_Cocktail_Raw_PtOrBin->Sumw2();

        // Lambda: clone cocktail spectra
        histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt               = (TH1D*) histoGammaTrueSecCocktailGammaFromXFromLambda_Pt->Clone("histoGammaTrueSecConvGammaFromXFromLambda_Cocktail_Raw_Pt");
        histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt->Sumw2();
        histoGammaSecGammaFromXFromLambda_Cocktail_Raw_PtOrBin          = (TH1D*) histoGammaTrueSecCocktailGammaFromXFromLambda_PtOrBin->Clone("histoGammaTrueSecConvGammaFromXFromLambda_Cocktail_Raw_PtOrBin");
        histoGammaSecGammaFromXFromLambda_Cocktail_Raw_PtOrBin->Sumw2();
        
        if ( isPCM && !isCalo ) {

            // K0s: calculate raw yield
            if (useUnfoldingForCocktailSecondaries) {
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt, histoGammaSecondaryFromXFromK0sConvProb_MCPt,histoGammaSecondaryFromXFromK0sRecoEff_MCPt,histoGammaTrueSecondaryFromXFromK0s_MCPt_recPt,nEvt,kTRUE,nIterationsUnfolding);
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_PtOrBin,histoGammaSecondaryFromXFromK0sConvProb_MCPtOrBin,histoGammaSecondaryFromXFromK0sRecoEff_MCPtOrBin,histoGammaTrueSecondaryFromXFromK0s_MCPt_recPtOrBin,nEvt,kTRUE,nIterationsUnfolding);
            } else {
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt, histoGammaSecondaryFromXFromK0sConvProb_MCPt,histoGammaSecondaryFromXFromK0sRecoEff_Pt,histoGammaTrueSecondaryFromXFromK0s_MCPt_recPt,nEvt,kFALSE);
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_PtOrBin,histoGammaSecondaryFromXFromK0sConvProb_MCPtOrBin,histoGammaSecondaryFromXFromK0sRecoEff_PtOrBin,histoGammaTrueSecondaryFromXFromK0s_MCPt_recPtOrBin,nEvt,kFALSE);
            }
   
            // K0l: calculate raw yield
            hasCocktailInput                                            = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt,histoGammaSecondaryFromXFromK0lConvProb_MCPt,histoGammaSecondaryFromXFromK0lRecoEff_Pt,histoGammaTrueSecondaryFromXFromK0l_MCPt_recPt,nEvt,kFALSE);
            hasCocktailInput                                            = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromK0l_Cocktail_Raw_PtOrBin,histoGammaSecondaryFromXFromK0lConvProb_MCPtOrBin,histoGammaSecondaryFromXFromK0lRecoEff_PtOrBin,histoGammaTrueSecondaryFromXFromK0l_MCPt_recPtOrBin,nEvt,kFALSE);

            // Lambda: calculate raw yield
            hasCocktailInput                                            = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt,histoGammaSecondaryFromXFromLambdaConvProb_MCPt,histoGammaSecondaryFromXFromLambdaRecoEff_Pt,histoGammaTrueSecondaryFromXFromLambda_MCPt_recPt,nEvt,kFALSE);
            hasCocktailInput                                            = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromLambda_Cocktail_Raw_PtOrBin,histoGammaSecondaryFromXFromLambdaConvProb_MCPtOrBin,histoGammaSecondaryFromXFromLambdaRecoEff_PtOrBin,histoGammaTrueSecondaryFromXFromLambda_MCPt_recPtOrBin,nEvt,kFALSE);
        }
        
        if ( isCalo && !isPCM ) {
            // K0s: calculate raw yield
            if (useUnfoldingForCocktailSecondaries) {
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt,histoGammaSecondaryFromXFromK0sRecoEff_MCPt,histoGammaTrueSecondaryFromXFromK0s_MCPt_recPt,nEvt,kTRUE,nIterationsUnfolding);
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_PtOrBin,histoGammaSecondaryFromXFromK0sRecoEff_MCPtOrBin,histoGammaTrueSecondaryFromXFromK0s_MCPt_recPtOrBin,nEvt,kTRUE,nIterationsUnfolding);
            } else {
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt,histoGammaSecondaryFromXFromK0sRecoEff_Pt,histoGammaTrueSecondaryFromXFromK0s_MCPt_recPt,nEvt,kFALSE);
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_PtOrBin,histoGammaSecondaryFromXFromK0sRecoEff_PtOrBin,histoGammaTrueSecondaryFromXFromK0s_MCPt_recPtOrBin,nEvt,kFALSE);
            }
            
            // K0l: calculate raw yield
            hasCocktailInput                                            = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt,histoGammaSecondaryFromXFromK0lRecoEff_Pt,histoGammaTrueSecondaryFromXFromK0l_MCPt_recPt,nEvt,kFALSE);
            hasCocktailInput                                            = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromK0l_Cocktail_Raw_PtOrBin,histoGammaSecondaryFromXFromK0lRecoEff_PtOrBin,histoGammaTrueSecondaryFromXFromK0l_MCPt_recPtOrBin,nEvt,kFALSE);
            
            // Lambda: calculate raw yield
            hasCocktailInput                                            = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt,histoGammaSecondaryFromXFromLambdaRecoEff_Pt,histoGammaTrueSecondaryFromXFromLambda_MCPt_recPt,nEvt,kFALSE);
            hasCocktailInput                                            = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromXFromLambda_Cocktail_Raw_PtOrBin,histoGammaSecondaryFromXFromLambdaRecoEff_PtOrBin,histoGammaTrueSecondaryFromXFromLambda_MCPt_recPtOrBin,nEvt,kFALSE);
        }
    }

    //******************* Calculate secondary fractions from cocktail input ********************
    TH1D*   histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt              = NULL;
    TH1D*   histoFracAllGammaToSecFromXFromK0s_Cocktail_PtOrBin         = NULL;
    TH1D*   histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt              = NULL;
    TH1D*   histoFracAllGammaToSecFromXFromK0l_Cocktail_PtOrBin         = NULL;
    TH1D*   histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt           = NULL;
    TH1D*   histoFracAllGammaToSecFromXFromLambda_Cocktail_PtOrBin      = NULL;
    TH1D*   histoFracAllGammaToSecRest_Cocktail_Pt                      = NULL;
    //TH1D*   histoFracAllGammaToSecRest_Cocktail_PtOrBin                 = NULL;
    if ( hasCocktailInput && isPCM && !isCalo ) {
        cout << "calculating secondary fractions from cocktail" << endl;
        
        histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt                  = (TH1D*)histoESDConvGammaPt->Clone("FracAllGammaToSecFromXFromK0s");
        histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt->Divide(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt,histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt,1,1,"B");
        histoFracAllGammaToSecFromXFromK0s_Cocktail_PtOrBin             = (TH1D*)histoESDConvGammaPt_OriginalBin->Clone("FracAllGammaToSecFromXFromK0sOriginalBinning");
        histoFracAllGammaToSecFromXFromK0s_Cocktail_PtOrBin->Divide(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_PtOrBin,histoFracAllGammaToSecFromXFromK0s_Cocktail_PtOrBin,1,1,"B");
        
        histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt                  = (TH1D*)histoESDConvGammaPt->Clone("FracAllGammaToSecFromXFromK0l");
        histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt->Divide(histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt,histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt,1,1,"B");
        histoFracAllGammaToSecFromXFromK0l_Cocktail_PtOrBin             = (TH1D*)histoESDConvGammaPt_OriginalBin->Clone("FracAllGammaToSecFromXFromK0lOriginalBinning");
        histoFracAllGammaToSecFromXFromK0l_Cocktail_PtOrBin->Divide(histoGammaSecGammaFromXFromK0l_Cocktail_Raw_PtOrBin,histoFracAllGammaToSecFromXFromK0l_Cocktail_PtOrBin,1,1,"B");
        
        histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt               = (TH1D*)histoESDConvGammaPt->Clone("FracAllGammaToSecFromXFromLambda");
        histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt->Divide(histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt,histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt,1,1,"B");
        histoFracAllGammaToSecFromXFromLambda_Cocktail_PtOrBin          = (TH1D*)histoESDConvGammaPt_OriginalBin->Clone("FracAllGammaToSecFromXFromLambdaOriginalBinning");
        histoFracAllGammaToSecFromXFromLambda_Cocktail_PtOrBin->Divide(histoGammaSecGammaFromXFromLambda_Cocktail_Raw_PtOrBin,histoFracAllGammaToSecFromXFromLambda_Cocktail_PtOrBin,1,1,"B");
        
        histoFracAllGammaToSecRest_Cocktail_Pt                          = (TH1D*)histoESDConvGammaPt->Clone("FracAllGammaToSecRest");
        histoFracAllGammaToSecRest_Cocktail_Pt->Divide(histoGammaTrueSecCocktailGammaRest_Pt,histoFracAllGammaToSecRest_Cocktail_Pt,1,1,"B");
    }
    if ( hasCocktailInput && isCalo && !isPCM ) {
        cout << "calculating secondary fractions from cocktail" << endl;
        
        histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt                  = (TH1D*)histoESDCaloGammaPt->Clone("FracAllGammaToSecFromXFromK0s");
        histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt->Divide(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt,histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt,1,1,"B");
        histoFracAllGammaToSecFromXFromK0s_Cocktail_PtOrBin             = (TH1D*)histoESDCaloGammaPt_OriginalBin->Clone("FracAllGammaToSecFromXFromK0sOriginalBinning");
        histoFracAllGammaToSecFromXFromK0s_Cocktail_PtOrBin->Divide(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_PtOrBin,histoFracAllGammaToSecFromXFromK0s_Cocktail_PtOrBin,1,1,"B");
        
        histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt                  = (TH1D*)histoESDCaloGammaPt->Clone("FracAllGammaToSecFromXFromK0l");
        histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt->Divide(histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt,histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt,1,1,"B");
        histoFracAllGammaToSecFromXFromK0l_Cocktail_PtOrBin             = (TH1D*)histoESDCaloGammaPt_OriginalBin->Clone("FracAllGammaToSecFromXFromK0lOriginalBinning");
        histoFracAllGammaToSecFromXFromK0l_Cocktail_PtOrBin->Divide(histoGammaSecGammaFromXFromK0l_Cocktail_Raw_PtOrBin,histoFracAllGammaToSecFromXFromK0l_Cocktail_PtOrBin,1,1,"B");
        
        histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt               = (TH1D*)histoESDCaloGammaPt->Clone("FracAllGammaToSecFromXFromLambda");
        histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt->Divide(histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt,histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt,1,1,"B");
        histoFracAllGammaToSecFromXFromLambda_Cocktail_PtOrBin          = (TH1D*)histoESDCaloGammaPt_OriginalBin->Clone("FracAllGammaToSecFromXFromLambdaOriginalBinning");
        histoFracAllGammaToSecFromXFromLambda_Cocktail_PtOrBin->Divide(histoGammaSecGammaFromXFromLambda_Cocktail_Raw_PtOrBin,histoFracAllGammaToSecFromXFromLambda_Cocktail_PtOrBin,1,1,"B");
        
        histoFracAllGammaToSecRest_Cocktail_Pt                          = (TH1D*)histoESDCaloGammaPt->Clone("FracAllGammaToSecRest");
        histoFracAllGammaToSecRest_Cocktail_Pt->Divide(histoGammaTrueSecCocktailGammaRest_Pt,histoFracAllGammaToSecRest_Cocktail_Pt,1,1,"B");
    }

    //******************* Pileup correction factors ********************************************
    TFile*  doPileUpCorr                                                = NULL;
    if (kDoPileup){
        doPileUpCorr                                                    = new  TFile(Form("%s/%s/%s_%s_GammaConvV1DCAHistogramms%s_%s.root", cutSelection.Data(), option.Data(), nameMeson.Data(),textPrefix2.Data(), optionPeriod.Data(), cutSelection.Data()));
        if (doPileUpCorr->IsZombie())   doPileUpCorr                    = 0;
    } else                              doPileUpCorr                    = 0;
    
    TH1D*   histoESDConvGammaPtPileUp                                   = NULL;
    TH1D*   histoPhotonDCAzFullPt                                       = NULL;
    TH1D*   histoRatioWithWithoutPileUp                                 = NULL;
    TF1*    histoRatioWithWithoutPileUpFit                              = NULL;
    TH1D*   histoPileUpCorrectionFactor                                 = NULL;
    TH1D*   histoPileUpCorrectionFactorOrBin                            = NULL;
    if(doPileUpCorr && isPCM){
        histoPhotonDCAzFullPt                                           = (TH1D*)fileUnCorrected->Get("ESD_GammaPtDCAzBin_Full");
        histoESDConvGammaPtPileUp                                       = (TH1D*)fileUnCorrected->Get("ESD_ConvGamma_Pt_PileUp");
        histoPileUpCorrectionFactor                                     = (TH1D*)histoESDConvGammaPtPileUp->Clone("PileUpCorrectionFactor");
        histoPileUpCorrectionFactor->Divide(histoPileUpCorrectionFactor,histoESDConvGammaPt,1,1,"B");
        
        // fit correction factor to get back to original binning
        cout << "fitting ratio with to without pileup to extract correction factor in original binning" << endl;
        
        histoRatioWithWithoutPileUp                                     = (TH1D*)histoESDConvGammaPt->Clone("histoRatioWithWithoutPileUp");
        histoRatioWithWithoutPileUp->Divide(histoRatioWithWithoutPileUp,histoESDConvGammaPtPileUp,1,1,"B");
        
        Int_t   fitStatus                                               = 0;
                histoRatioWithWithoutPileUpFit                          = new TF1("histoRatioWithWithoutPileUpFit", "1+[0]/TMath::Power((x-[1]), [2])",
                                                                                  histoESDConvGammaPt_OriginalBin->GetXaxis()->GetXmin(),
                                                                                  histoESDConvGammaPt_OriginalBin->GetXaxis()->GetXmax());
        histoRatioWithWithoutPileUpFit->SetParameters(1, 0, 1);
        histoRatioWithWithoutPileUpFit->SetName(Form("%s_fit", histoRatioWithWithoutPileUp->GetName()));
        TFitResultPtr histoRatioWithWithoutPileUpFitResult              = histoRatioWithWithoutPileUp->Fit(histoRatioWithWithoutPileUpFit, "SMNRE+","",
                                                                                                           histoRatioWithWithoutPileUp->GetXaxis()->GetXmin(),
                                                                                                           histoRatioWithWithoutPileUp->GetXaxis()->GetXmax());
        fitStatus                                                       = histoRatioWithWithoutPileUpFitResult;
        
        if (fitStatus == 0 || fitStatus >= 1000) {
            // accepting fits, if everything went fine (i.e. fitstatus = 0) of if only improve had problems (i.e. fitstatus >= 1000)
            cout << "fit status: " << fitStatus << endl;
            
            histoPileUpCorrectionFactorOrBin                            = (TH1D*)histoESDConvGammaPt_OriginalBin->Clone(Form("%s_OriginalBinning", histoPileUpCorrectionFactor->GetName()));
            histoPileUpCorrectionFactorOrBin->Reset("ICES");
            histoPileUpCorrectionFactorOrBin->Sumw2();
            
            Double_t binContent, binError;
            for (Int_t i=1; i<histoPileUpCorrectionFactorOrBin->GetNbinsX()+1; i++) {
                binContent                                              = histoRatioWithWithoutPileUpFit->Integral( histoPileUpCorrectionFactorOrBin->GetXaxis()->GetBinLowEdge(i),
                                                                                                                    histoPileUpCorrectionFactorOrBin->GetXaxis()->GetBinUpEdge(i),
                                                                                                                    histoRatioWithWithoutPileUpFitResult->GetParams()) / histoPileUpCorrectionFactorOrBin->GetBinWidth(i);
                
                binError                                                = histoRatioWithWithoutPileUpFit->IntegralError(    histoPileUpCorrectionFactorOrBin->GetXaxis()->GetBinLowEdge(i),
                                                                                                                            histoPileUpCorrectionFactorOrBin->GetXaxis()->GetBinUpEdge(i),
                                                                                                                            histoRatioWithWithoutPileUpFitResult->GetParams(),
                                                                                                                            histoRatioWithWithoutPileUpFitResult->GetCovarianceMatrix().GetMatrixArray()) / histoPileUpCorrectionFactorOrBin->GetBinWidth(i);
                
                if (binContent > 0) {
                    histoPileUpCorrectionFactorOrBin->SetBinContent(i, 1/binContent);
                    histoPileUpCorrectionFactorOrBin->SetBinError(  i, binError/binContent/binContent);
                } else {
                    histoPileUpCorrectionFactorOrBin->SetBinContent(i, 1);
                    histoPileUpCorrectionFactorOrBin->SetBinError(  i, 0);
                }
            }

        } else {
            cout << "fit failed, no pileup correction will be applied to unfolded spectra" << endl;
        }
    }

    //******************* MC pileup histograms *************************************************
    TH1D*   histoGammaPurity_PileUp_Pt                                  = NULL;
    TH1D*   histoGammaTruePurity_PileUp_Pt                              = NULL;
    TH1D*   histoGammaRecoEff_PileUp_Pt                                 = NULL;
    TH1D*   histoGammaPrimaryRecoEff_PileUp_Pt                          = NULL;
    if( doPileUpCorr && isPCM ){
        histoGammaPurity_PileUp_Pt                                      = (TH1D*)fileCorrections->Get("GammaPurity_PileUp_Pt");
        histoGammaTruePurity_PileUp_Pt                                  = (TH1D*)fileCorrections->Get("GammaTruePurity_PileUp_Pt");
        histoGammaRecoEff_PileUp_Pt                                     = (TH1D*)fileCorrections->Get("GammaRecoEff_PileUp_Pt");
        histoGammaPrimaryRecoEff_PileUp_Pt                              = (TH1D*)fileCorrections->Get("GammaPrimaryRecoEff_PileUp_Pt");
    }
    TH1D*   histoMCrecGamma_PileUp_Pt                                   = NULL;
    TH1D*   histoPileUpCorrectionFactorMC                               = NULL;
    if(doPileUpCorr && isPCM){
        histoMCrecGamma_PileUp_Pt                                       = (TH1D*)fileCorrections->Get("MCrec_ConvGamma_Pt_PileUp");
        if (histoMCrecGamma_PileUp_Pt){
            histoPileUpCorrectionFactorMC                               = (TH1D*)histoMCrecGamma_PileUp_Pt->Clone("PileUpCorrectionFactorMC");
            histoPileUpCorrectionFactorMC->Divide(histoPileUpCorrectionFactorMC,histoMCrecGamma_Pt,1,1,"B");
        }
    }
    TH1D*   histoFracAllGammaToSec_PileUp_Pt                            = NULL;
    TH1D*   histoFracAllGammaToSecFromXFromK0s_PileUp_Pt                = NULL;
    TH1D*   histoMCrecPhotonDCAzFullPt                                  = NULL;
    if(doPileUpCorr&&isPCM){
        histoFracAllGammaToSec_PileUp_Pt                                = (TH1D*)fileCorrections->Get("FracAllGammaToSecPileUp");
        histoFracAllGammaToSecFromXFromK0s_PileUp_Pt                    = (TH1D*)fileCorrections->Get("FracAllGammaToSecFromXFromK0sPileUp");
        histoMCrecPhotonDCAzFullPt                                      = (TH1D*)fileCorrections->Get("MCrec_GammaPtDCAzBin_Full");     // category must be adapted
    }

    //******************* Determine max Pt *****************************************************
    Double_t maxPtGamma                                                 = 0.;
    if (isPCM && !isCalo) maxPtGamma                                    = histoESDConvGammaPt->GetXaxis()->GetBinUpEdge(histoESDConvGammaPt->GetNbinsX());
    if (isCalo && !isPCM) maxPtGamma                                    = histoESDCaloGammaPt->GetXaxis()->GetBinUpEdge(histoESDCaloGammaPt->GetNbinsX());

    //******************* Proper scaling of background *****************************************
    TH1D *ScalingGammaBackground_Pt                                     = NULL;
    if (isPCM && !isCalo) {
        ScalingGammaBackground_Pt                                       = (TH1D*) histoESDConvGammaPt->Clone("ScalingGammaBackground_Pt");
        ScalingGammaBackground_Pt->Divide(ScalingGammaBackground_Pt, histoMCrecGamma_Pt, 1., 1, "");
    }
    if (isCalo && !isPCM) {
        ScalingGammaBackground_Pt                                       = (TH1D*) histoESDCaloGammaPt->Clone("ScalingGammaBackground_Pt");
        ScalingGammaBackground_Pt->Divide(ScalingGammaBackground_Pt, histoMCrecGammaCalo_Pt, 1., 1, "");
    }
    ScalingGammaBackground_Pt->Scale(nEvtMC/nEvt);
    histoMCrecBackground_Pt->Scale(1./nEvtMC);
    TH1D *histoGammaMCBackground_Pt                                     = (TH1D*)histoMCrecBackground_Pt->Clone("histoGammaMCBackground_Pt");
    histoMCrecBackground_Pt->Multiply(ScalingGammaBackground_Pt);

    //******************* Calculate secondary spectra from data ********************************
    TH1D *histoSecondaryGammaSpecPt                                     = NULL;
    TH1D *histoSecondaryGammaFromXFromK0sSpecPt                         = NULL;
    TH1D *histoSecondaryGammaFromXFromK0lSpecPt                         = NULL;
    TH1D *histoSecondaryGammaFromXFromLambdaSpecPt                      = NULL;
    if (isPCM && !isCalo) {
        histoSecondaryGammaSpecPt                                       = (TH1D*)histoESDConvGammaPt->Clone("SecondaryGammaSpecPt");
        histoSecondaryGammaFromXFromK0sSpecPt                           = (TH1D*)histoESDConvGammaPt->Clone("SecondaryGammaSpecFromXFromK0sPt");
        histoSecondaryGammaFromXFromK0lSpecPt                           = (TH1D*)histoESDConvGammaPt->Clone("SecondaryGammaSpecFromXFromK0lPt");
        histoSecondaryGammaFromXFromLambdaSpecPt                        = (TH1D*)histoESDConvGammaPt->Clone("SecondaryGammaSpecFromXFromLambdaPt");
    }
    if (isCalo && !isPCM) {
        histoSecondaryGammaSpecPt                                       = (TH1D*)histoESDCaloGammaPt->Clone("SecondaryGammaSpecPt");
        histoSecondaryGammaFromXFromK0sSpecPt                           = (TH1D*)histoESDCaloGammaPt->Clone("SecondaryGammaSpecFromXFromK0sPt");
        histoSecondaryGammaFromXFromK0lSpecPt                           = (TH1D*)histoESDCaloGammaPt->Clone("SecondaryGammaSpecFromXFromK0lPt");
        histoSecondaryGammaFromXFromLambdaSpecPt                        = (TH1D*)histoESDCaloGammaPt->Clone("SecondaryGammaSpecFromXFromLambdaPt");
    }
    
    histoSecondaryGammaSpecPt->Multiply(histoFracAllGammaToSec_Pt);
    histoSecondaryGammaFromXFromK0sSpecPt->Multiply(histoFracAllGammaToSecFromXFromK0s_Pt);
    if (histoFracAllGammaToSecFromXFromK0l_Pt)  histoSecondaryGammaFromXFromK0lSpecPt->Multiply(histoFracAllGammaToSecFromXFromK0l_Pt);
    else                                        histoSecondaryGammaFromXFromK0lSpecPt = NULL;
    histoSecondaryGammaFromXFromLambdaSpecPt->Multiply(histoFracAllGammaToSecFromXFromLambda_Pt);
    
    histoSecondaryGammaSpecPt->Scale(1./nEvt);
    histoSecondaryGammaFromXFromK0sSpecPt->Scale(1./nEvt);
    histoSecondaryGammaFromXFromK0sSpecPt->Scale(doubleAddFactorK0s);
    if (histoSecondaryGammaFromXFromK0lSpecPt) histoSecondaryGammaFromXFromK0lSpecPt->Scale(1./nEvt);
    histoSecondaryGammaFromXFromLambdaSpecPt->Scale(1./nEvt);

    //******************* Scale MC/cocktail secondary spectra **********************************
    if (isPCM && !isCalo) {
        histoGammaTrueSecConv_Pt->Scale(1./nEvtMC);
        histoGammaTrueSecConvGammaFromXFromK0s_Pt->Scale(1./nEvtMC);
        if (histoGammaTrueSecConvGammaFromXFromK0l_Pt) histoGammaTrueSecConvGammaFromXFromK0l_Pt->Scale(1./nEvtMC);
        histoGammaTrueSecConvGammaFromXFromLambda_Pt->Scale(1./nEvtMC);
    }
    if (isCalo && !isPCM) {
        histoGammaTrueSecCalo_Pt->Scale(1./nEvtMC);
        histoGammaTrueSecCaloGammaFromXFromK0s_Pt->Scale(1./nEvtMC);
        if (histoGammaTrueSecCaloGammaFromXFromK0l_Pt) histoGammaTrueSecCaloGammaFromXFromK0l_Pt->Scale(1./nEvtMC);
        histoGammaTrueSecCaloGammaFromXFromLambda_Pt->Scale(1./nEvtMC);
    }
    if (hasCocktailInput) {
        histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt->Scale(1./nEvt);
        histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt->Scale(1./nEvt);
        histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt->Scale(1./nEvt);
        histoGammaTrueSecCocktailGammaRest_Pt->Scale(1./nEvtMC);
    }

    //******************* Calculate pileup corr. secondary spectra from data *******************
    TH1D *histoSecondaryGammaSpecPtPileUp                               = NULL;
    TH1D *histoSecondaryGammaFromXFromK0sSpecPtPileUp                   = NULL;
    if(doPileUpCorr && isPCM){
        histoSecondaryGammaSpecPtPileUp                                 = (TH1D*)histoESDConvGammaPtPileUp->Clone("SecondaryGammaSpecPtPileUp");
        histoSecondaryGammaSpecPtPileUp->Multiply(histoFracAllGammaToSec_PileUp_Pt);
        histoSecondaryGammaSpecPtPileUp->Scale(1./nEvt);
        histoSecondaryGammaFromXFromK0sSpecPtPileUp                     = (TH1D*)histoESDConvGammaPtPileUp->Clone("SecondaryGammaSpecFromXFromK0sPtPileUp");
        histoSecondaryGammaFromXFromK0sSpecPtPileUp->Multiply(histoFracAllGammaToSecFromXFromK0s_PileUp_Pt);
        histoSecondaryGammaFromXFromK0sSpecPtPileUp->Scale(1./nEvt);
        histoSecondaryGammaFromXFromK0sSpecPtPileUp->Scale(doubleAddFactorK0s);
    }

    //**********************************************************************************
    //******************** PrimVtx DCA Plot ********************************************
    //**********************************************************************************
    if( doPileUpCorr && isPCM ){
        TCanvas *canvasPileUpCorrFactor     = GetAndSetCanvas("canvasPileUpCorrFactor");

            DrawGammaSetMarker(histoPileUpCorrectionFactor, 20, 3, 1, 1);
            if (histoMCrecGamma_PileUp_Pt)DrawGammaSetMarker(histoPileUpCorrectionFactorMC, 24, 3, 2, 2);

            SetHistogramm(histoPileUpCorrectionFactor,"#it{p}_{T} (GeV/#it{c})","Correction Factor (%)",0.84,1.02);
            if (histoMCrecGamma_PileUp_Pt)SetHistogramm(histoPileUpCorrectionFactorMC,"#it{p}_{T} (GeV/#it{c})","Correction Factor (%)",0.84,1.02);

            histoPileUpCorrectionFactor->DrawCopy("");
            if (histoMCrecGamma_PileUp_Pt)histoPileUpCorrectionFactorMC->Draw("same");

            TLegend* legendPileUpCorrFactor = GetAndSetLegend(0.3,0.2,2.2,1,cent);
            legendPileUpCorrFactor->AddEntry(histoPileUpCorrectionFactor,"Correction Factor Data","lp");
            if (histoMCrecGamma_PileUp_Pt)legendPileUpCorrFactor->AddEntry(histoPileUpCorrectionFactorMC,"Correction Factor MC","lp");
            legendPileUpCorrFactor->Draw();
        
        canvasPileUpCorrFactor->SaveAs(Form("%s/%s_PileUpCorrFactor_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasPileUpCorrFactor;
        
        Bool_t doPileUpCorrSimplePlot       = kTRUE;
        if (doPileUpCorrSimplePlot && !textPrefix2.CompareTo("data")) {
            TCanvas *canvasPileUpCorrFactor2= GetAndSetCanvas("canvasPileUpCorrFactor2");
            
            SetHistogramm(histoPileUpCorrectionFactor,"#it{p}_{T} (GeV/#it{c})","Correction Factor (%)",0.90,1.02);
            DrawGammaSetMarker(histoPileUpCorrectionFactor, 24, 1.5, kBlack, kBlack);
            histoPileUpCorrectionFactor->DrawCopy("");
            DrawGammaLines(0., maxPtGamma,1.0, 1.0, 1, kGray+2, 2);
            histoPileUpCorrectionFactor->DrawCopy("same");

            PutProcessLabelAndEnergyOnPlot( 0.15, 0.3, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
            
            canvasPileUpCorrFactor2->SaveAs(Form("%s/%s_PileUpCorrFactorDataOnly_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data(),suffix.Data()));
            delete canvasPileUpCorrFactor2;
        }
        
        if (histoPileUpCorrectionFactorOrBin) {
            
            // rebin pileup correction factor from original binning
            TH1D* histoPileUpCorrectionFactorTemp           = (TH1D*)histoPileUpCorrectionFactorOrBin->Clone("histoPileUpCorrectionFactorTemp");
            histoPileUpCorrectionFactorTemp                 = RebinTH1D(histoPileUpCorrectionFactorTemp,histoPileUpCorrectionFactor,kFALSE);
            histoPileUpCorrectionFactorTemp->SetBinContent( 1, 0);
            histoPileUpCorrectionFactorTemp->SetBinError(   1, 0);
            
            // plot
            TCanvas *canvasPileUpCorrFactor3                = GetAndSetCanvas("canvasPileUpCorrFactor3");
            
            SetHistogramm(histoPileUpCorrectionFactorTemp,"#it{p}_{T} (GeV/#it{c})","Correction Factor (%)",0.90,1.02);
            DrawGammaSetMarker(histoPileUpCorrectionFactorTemp, 24, 1.5, kBlack, kBlack);
            histoPileUpCorrectionFactorTemp->DrawCopy("");
            DrawGammaLines(0., maxPtGamma,1.0, 1.0, 1, kGray+2, 2);
            histoPileUpCorrectionFactorTemp->DrawCopy("same");
            
            PutProcessLabelAndEnergyOnPlot( 0.15, 0.3, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
            
            canvasPileUpCorrFactor3->SaveAs(Form("%s/%s_PileUpCorrFactorDataOnlyRebinned_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data(),suffix.Data()));
            delete canvasPileUpCorrFactor3;
        }
    }

    //**********************************************************************************
    //******************** Background Plot *********************************************
    //**********************************************************************************
    TCanvas *canvasBackground               = GetAndSetCanvas("canvasBackground");
    canvasBackground->SetLogy();

        TH1D* histoGammaRawSpectrum_PtMC    = NULL;
        TH1D* histoGammaRawSpectrum_Pt      = NULL;
        TLegend* legendBackground           = GetAndSetLegend(0.6,0.7,2,1);
        if (isPCM && !isCalo) {
            histoGammaRawSpectrum_PtMC      = (TH1D*) histoMCrecGamma_Pt->Clone("histoGammaRawSpectrum_PtMC");
            histoGammaRawSpectrum_PtMC->Scale(1./nEvtMC);
            histoGammaRawSpectrum_Pt        = (TH1D*) histoESDConvGammaPt->Clone("histoGammaRawSpectrum_Pt");
            histoGammaRawSpectrum_Pt->Scale(1./nEvt);
            
            DrawGammaSetMarker(histoGammaRawSpectrum_Pt, 20, 1.0, 1, 1);
            DrawGammaSetMarker(histoMCrecBackground_Pt, 20, 1.0, 2, 2);
            
            SetHistogramm(histoMCrecBackground_Pt,"#it{p}_{T} (GeV/#it{c})",Form("Background for #gamma in |#eta| < %g",eta));
            SetHistogramm(histoGammaRawSpectrum_Pt,"#it{p}_{T} (GeV/#it{c})","Raw yield",1e-8,5);
            
            histoGammaRawSpectrum_Pt->DrawCopy("");
            histoMCrecBackground_Pt->Draw("same");
        }
        if (isCalo && !isPCM) {
            histoGammaRawSpectrum_PtMC      = (TH1D*) histoMCrecGammaCalo_Pt->Clone("histoGammaRawSpectrum_PtMC");
            histoGammaRawSpectrum_PtMC->Scale(1./nEvtMC);
            histoGammaRawSpectrum_Pt        = (TH1D*) histoESDCaloGammaPt->Clone("histoGammaRawSpectrum_Pt");
            histoGammaRawSpectrum_Pt->Scale(1./nEvt);
            
            DrawGammaSetMarker(histoGammaRawSpectrum_Pt, 20, 1.0, 1, 1);
            DrawGammaSetMarker(histoMCrecBackground_Pt, 20, 1.0, 2, 2);
            
            SetHistogramm(histoMCrecBackground_Pt,"#it{p}_{T} (GeV/#it{c})",Form("Background for #gamma in |#eta| < %g",etaCalo));
            SetHistogramm(histoGammaRawSpectrum_Pt,"#it{p}_{T} (GeV/#it{c})","Raw yield",1e-10,1e-2);
            
            histoGammaRawSpectrum_Pt->DrawCopy("");
            histoMCrecBackground_Pt->Draw("same");
        }

        legendBackground->AddEntry(histoGammaRawSpectrum_Pt,"raw #gamma spectrum","pl");
        legendBackground->AddEntry(histoMCrecBackground_Pt,"#gamma background","pl");
        legendBackground->Draw();
    
        PutProcessLabelAndEnergyOnPlot( 0.6, 0.95, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);

    canvasBackground->SaveAs(Form("%s/%s_Background_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasBackground;

    //**********************************************************************************
    //******************** Secondary Spectra Plot **************************************
    //**********************************************************************************
    TCanvas *canvasSecSpec                      = GetAndSetCanvas("canvasSecSpec");
    canvasSecSpec->SetLogy();
    
        TLegend* legendSecSpec                  = NULL;
        if (!hasCocktailInput)  legendSecSpec   = GetAndSetLegend(0.45,0.7,5);
        else                    legendSecSpec   = GetAndSetLegend(0.45,0.5,8);

        // Plotting MC spectra
        if (isPCM && !isCalo) {
            SetHistogramm(histoGammaTrueSecConv_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma",1e-10,1e-1);
            SetHistogramm(histoGammaTrueSecConvGammaFromXFromK0s_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
            if (histoGammaTrueSecConvGammaFromXFromK0l_Pt) SetHistogramm(histoGammaTrueSecConvGammaFromXFromK0l_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
            SetHistogramm(histoGammaTrueSecConvGammaFromXFromLambda_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
            
            DrawGammaSetMarker(histoGammaTrueSecConv_Pt, 20, 1.0, kRed+1, kRed+1);
            DrawGammaSetMarker(histoGammaTrueSecConvGammaFromXFromK0s_Pt, 24, 1.0, kRed-6, kRed-6);
            if (histoGammaTrueSecConvGammaFromXFromK0l_Pt) DrawGammaSetMarker(histoGammaTrueSecConvGammaFromXFromK0l_Pt, 22, 1.0, kRed-6, kRed-6);
            DrawGammaSetMarker(histoGammaTrueSecConvGammaFromXFromLambda_Pt, 21, 1.0, kRed-6, kRed-6);
            
            histoGammaTrueSecConv_Pt->Draw();
            histoGammaTrueSecConvGammaFromXFromK0s_Pt->Draw("same");
            if (histoGammaTrueSecConvGammaFromXFromK0l_Pt) histoGammaTrueSecConvGammaFromXFromK0l_Pt->Draw("same");
            histoGammaTrueSecConvGammaFromXFromLambda_Pt->Draw("same");
            
            legendSecSpec->AddEntry(histoGammaTrueSecConv_Pt,"Raw MC Sec #gamma Spectrum","pl");
            legendSecSpec->AddEntry(histoGammaTrueSecConvGammaFromXFromK0s_Pt,"Raw MC Sec #gamma from X from K^{0}_{S}","pl");
            if (histoGammaTrueSecConvGammaFromXFromK0l_Pt) legendSecSpec->AddEntry(histoGammaTrueSecConvGammaFromXFromK0l_Pt,"Raw MC Sec #gamma from X from K^{0}_{L}","pl");
            legendSecSpec->AddEntry(histoGammaTrueSecConvGammaFromXFromLambda_Pt,"Raw MC Sec #gamma from X from #Lambda","pl");
        }

        if (isCalo && !isPCM) {
            SetHistogramm(histoGammaTrueSecCalo_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma",1e-10,1e-1);
            SetHistogramm(histoGammaTrueSecCaloGammaFromXFromK0s_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary #gamma");
            if (histoGammaTrueSecCaloGammaFromXFromK0l_Pt) SetHistogramm(histoGammaTrueSecCaloGammaFromXFromK0s_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary #gamma");
            SetHistogramm(histoGammaTrueSecCaloGammaFromXFromLambda_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary #gamma");
            
            DrawGammaSetMarker(histoGammaTrueSecCalo_Pt, 20, 1.0, kRed+1, kRed+1);
            DrawGammaSetMarker(histoGammaTrueSecCaloGammaFromXFromK0s_Pt, 24, 1.0, kRed-6, kRed-6);
            if (histoGammaTrueSecCaloGammaFromXFromK0l_Pt) DrawGammaSetMarker(histoGammaTrueSecCaloGammaFromXFromK0l_Pt, 22, 1.0, kRed-6, kRed-6);
            DrawGammaSetMarker(histoGammaTrueSecCaloGammaFromXFromLambda_Pt, 21, 1.0, kRed-6, kRed-6);
            
            histoGammaTrueSecCalo_Pt->Draw();
            histoGammaTrueSecCaloGammaFromXFromK0s_Pt->Draw("same");
            if (histoGammaTrueSecCaloGammaFromXFromK0l_Pt) histoGammaTrueSecCaloGammaFromXFromK0l_Pt->Draw("same");
            histoGammaTrueSecCaloGammaFromXFromLambda_Pt->Draw("same");
            
            legendSecSpec->AddEntry(histoGammaTrueSecCalo_Pt,"Raw MC Sec #gamma Spectrum","pl");
            legendSecSpec->AddEntry(histoGammaTrueSecCaloGammaFromXFromK0s_Pt,"Raw MC Sec #gamma from X from K^{0}_{S}","pl");
            if (histoGammaTrueSecCaloGammaFromXFromK0l_Pt) legendSecSpec->AddEntry(histoGammaTrueSecCaloGammaFromXFromK0l_Pt,"Raw MC Sec #gamma from X from K^{0}_{L}","pl");
            legendSecSpec->AddEntry(histoGammaTrueSecCaloGammaFromXFromLambda_Pt,"Raw MC Sec #gamma from X from #Lambda","pl");
        }

        if (hasCocktailInput) {
            if (isPCM && !isCalo) {
                SetHistogramm(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
                SetHistogramm(histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
                SetHistogramm(histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
            }
            
            if (isCalo && !isPCM) {
                SetHistogramm(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary #gamma");
                SetHistogramm(histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary #gamma");
                SetHistogramm(histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary #gamma");
            }
            
            DrawGammaSetMarker(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt, 20, 1.0, kBlue-6, kBlue-6);
            DrawGammaSetMarker(histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt, 24, 1.0, kBlue-4, kBlue-4);
            DrawGammaSetMarker(histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt, 34, 1.0, kBlue+2, kBlue+2);
            
            histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt->DrawCopy("same");
            histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt->DrawCopy("same");
            histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt->DrawCopy("same");
            
            legendSecSpec->AddEntry(histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt,"Raw cocktail Sec #gamma from X from K^{0}_{s}","pl");
            legendSecSpec->AddEntry(histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt,"Raw cocktail Sec #gamma from X from K^{0}_{l}","pl");
            legendSecSpec->AddEntry(histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt,"Raw cocktail Sec #gamma from X from #lambda","pl");
            
        } else {
            if (isPCM && !isCalo) {
                SetHistogramm(histoSecondaryGammaSpecPt, "#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
                SetHistogramm(histoSecondaryGammaFromXFromK0sSpecPt,"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
                if (histoSecondaryGammaFromXFromK0lSpecPt) SetHistogramm(histoSecondaryGammaFromXFromK0lSpecPt,"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
                SetHistogramm(histoSecondaryGammaFromXFromLambdaSpecPt,"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
            }
            
            if (isCalo && !isPCM) {
                SetHistogramm(histoSecondaryGammaSpecPt, "#it{p}_{T} (GeV/#it{c})","Secondary #gamma");
                SetHistogramm(histoSecondaryGammaFromXFromK0sSpecPt,"#it{p}_{T} (GeV/#it{c})","Secondary #gamma");
                if (histoSecondaryGammaFromXFromK0lSpecPt) SetHistogramm(histoSecondaryGammaFromXFromK0lSpecPt,"#it{p}_{T} (GeV/#it{c})","Secondary #gamma");
                SetHistogramm(histoSecondaryGammaFromXFromLambdaSpecPt,"#it{p}_{T} (GeV/#it{c})","Secondary #gamma");
            }
            
            DrawGammaSetMarker(histoSecondaryGammaSpecPt, 21, 1.0, kBlue+1, kBlue+1);
            DrawGammaSetMarker(histoSecondaryGammaFromXFromK0sSpecPt, 25, 1.0, kBlue-6, kBlue-6);
            if (histoSecondaryGammaFromXFromK0lSpecPt) DrawGammaSetMarker(histoSecondaryGammaFromXFromK0lSpecPt, 22, 1.0, kBlue-2, kBlue-2);
            DrawGammaSetMarker(histoSecondaryGammaFromXFromLambdaSpecPt, 20, 1.0, kBlue-4, kBlue-4);
            
            
            histoSecondaryGammaSpecPt->Draw("same");
            histoSecondaryGammaFromXFromK0sSpecPt->DrawCopy("same");
            if (histoSecondaryGammaFromXFromK0lSpecPt) histoSecondaryGammaFromXFromK0lSpecPt->DrawCopy("same");
            histoSecondaryGammaFromXFromLambdaSpecPt->DrawCopy("same");
            
            legendSecSpec->AddEntry(histoSecondaryGammaSpecPt,"Raw data Sec #gamma Spectrum ","pl");
            legendSecSpec->AddEntry(histoSecondaryGammaFromXFromK0sSpecPt,"Raw data Sec #gamma from X from K^{0}_{S}","pl");
            legendSecSpec->AddEntry((TObject*)0, "(scaled with K^{0}_{s} factor)","");
            if (histoSecondaryGammaFromXFromK0lSpecPt) legendSecSpec->AddEntry(histoSecondaryGammaFromXFromK0lSpecPt,"Raw data Sec #gamma from X from K^{0}_{L}","pl");
            legendSecSpec->AddEntry(histoSecondaryGammaFromXFromLambdaSpecPt,"Raw data Sec #gamma from X from #Lambda","pl");
        }
        legendSecSpec->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.18, 0.3, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
    canvasSecSpec->SaveAs(Form("%s/%s_SecondarySpectra_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasSecSpec;

    //**********************************************************************************
    //******************** Secondary Fractions Plot ************************************
    //**********************************************************************************
    TCanvas *canvasSecFrac                  = GetAndSetCanvas("canvasSecFrac");
    canvasSecFrac->SetTopMargin(0.035);
    
        TLegend* legendSecFrac              = GetAndSetLegend(0.6,0.70,5,1);

        if (!hasCocktailInput) {
            if (isPCM && !isCalo) {
                SetHistogramm(histoFracAllGammaToSec_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary Converted #gamma");
                SetHistogramm(histoFracAllGammaToSecFromXFromK0s_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary Converted #gamma");
                if (histoFracAllGammaToSecFromXFromK0l_Pt) SetHistogramm(histoFracAllGammaToSecFromXFromK0l_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary Converted #gamma");
                SetHistogramm(histoFracAllGammaToSecFromXFromLambda_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary Converted #gamma");
            }
            
            if (isCalo && !isPCM) {
                SetHistogramm(histoFracAllGammaToSec_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary #gamma");
                SetHistogramm(histoFracAllGammaToSecFromXFromK0s_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary #gamma");
                if (histoFracAllGammaToSecFromXFromK0l_Pt) SetHistogramm(histoFracAllGammaToSecFromXFromK0l_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary #gamma");
                SetHistogramm(histoFracAllGammaToSecFromXFromLambda_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary #gamma");
            }
        
            DrawGammaSetMarker(histoFracAllGammaToSec_Pt, 20, 1.0, 1, 1);
            DrawGammaSetMarker(histoFracAllGammaToSecFromXFromK0s_Pt, 24, 1.0, kBlue-6, kBlue-6);
            if (histoFracAllGammaToSecFromXFromK0l_Pt) DrawGammaSetMarker(histoFracAllGammaToSecFromXFromK0l_Pt, 22, 1.0, kBlue-6, kBlue-6);
            DrawGammaSetMarker(histoFracAllGammaToSecFromXFromLambda_Pt, 21, 1.0, kRed-6, kRed-6);

            histoFracAllGammaToSec_Pt->Draw();
            histoFracAllGammaToSecFromXFromK0s_Pt->Draw("same");
            if (histoFracAllGammaToSecFromXFromK0l_Pt) histoFracAllGammaToSecFromXFromK0l_Pt->Draw("same");
            histoFracAllGammaToSecFromXFromLambda_Pt->Draw("same");

            legendSecFrac->AddEntry(histoFracAllGammaToSec_Pt,"Fraction Sec #gamma","pl");
            legendSecFrac->AddEntry(histoFracAllGammaToSecFromXFromK0s_Pt,Form("Fraction Sec #gamma from X from K^{0}_{S}"),"pl");
            if (histoFracAllGammaToSecFromXFromK0l_Pt) legendSecFrac->AddEntry(histoFracAllGammaToSecFromXFromK0l_Pt,Form("Fraction Sec #gamma from X from K^{0}_{L}"),"pl");
            legendSecFrac->AddEntry(histoFracAllGammaToSecFromXFromLambda_Pt,Form("Fraction Sec #gamma from X from #Lambda"),"pl");
        } else {
            if (isPCM && !isCalo) {
                SetHistogramm(histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary Converted #gamma");
                SetHistogramm(histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary Converted #gamma");
                SetHistogramm(histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary Converted #gamma");
                SetHistogramm(histoFracAllGammaToSecRest_Cocktail_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary Converted #gamma");
            }
            
            if (isCalo && !isPCM) {
                SetHistogramm(histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary #gamma");
                SetHistogramm(histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary #gamma");
                SetHistogramm(histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary #gamma");
                SetHistogramm(histoFracAllGammaToSecRest_Cocktail_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary #gamma");
            }
            
            DrawGammaSetMarker(histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt,      20, 1.0, kBlue-6, kBlue-6);
            DrawGammaSetMarker(histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt,      21, 1.0, kRed-4, kRed-4);
            DrawGammaSetMarker(histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt,   20, 1.0, kGreen-1, kGreen-1);
            DrawGammaSetMarker(histoFracAllGammaToSecRest_Cocktail_Pt,              21, 1.0, kCyan+3, kCyan+3);

            histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt->Draw("same");
            histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt->Draw("same");
            histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt->Draw("same");
            histoFracAllGammaToSecRest_Cocktail_Pt->Draw("same");

            legendSecFrac->AddEntry(histoFracAllGammaToSecFromXFromK0s_Cocktail_Pt,     Form("sec. #gamma from X from K^{0}_{s}"),"pl");
            legendSecFrac->AddEntry(histoFracAllGammaToSecFromXFromK0l_Cocktail_Pt,     Form("sec. #gamma from X from K^{0}_{l}"),"pl");
            legendSecFrac->AddEntry(histoFracAllGammaToSecFromXFromLambda_Cocktail_Pt,  Form("sec. #gamma from X from #Lambda"),"pl");
            legendSecFrac->AddEntry(histoFracAllGammaToSecRest_Cocktail_Pt,             Form("sec. #gamma rest"),"pl");
        }
        legendSecFrac->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.75, 0.65, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
    canvasSecFrac->SaveAs(Form("%s/%s_SecondaryGammaFraction_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasSecFrac;

    //**********************************************************************************
    //******************** Purity Plot *************************************************
    //**********************************************************************************
    TCanvas *canvasPurity                       = GetAndSetCanvas("canvasPurity");

        TLegend* legendPurity                   = NULL;

        if ( (isPCM && !isCalo) || (isCalo && !isPCM) ){
            DrawGammaSetMarker(histoGammaPurity_Pt, 24, 1.0, kRed+2, kRed+2);
            DrawGammaSetMarker(histoGammaTruePurity_Pt, 20, 1.0, 1, 1);
            if(doPileUpCorr)DrawGammaSetMarker(histoGammaTruePurity_PileUp_Pt, 20, 1.0, kBlue+2, kBlue+2);
            
            if (isPCM && !isCalo) SetHistogramm(histoGammaTruePurity_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{pur,#gamma} in |#eta| < %g",eta),0.6, 1.);
            if (isCalo && !isPCM) SetHistogramm(histoGammaTruePurity_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{pur,#gamma} in |#eta| < %g",etaCalo),0.6, 1.);
            histoGammaTruePurity_Pt->Draw();
            histoGammaPurity_Pt->Draw("same");
            if (doPileUpCorr) histoGammaTruePurity_PileUp_Pt->Draw("same");
            
            if (doPileUpCorr)   legendPurity    = GetAndSetLegend(0.18,0.15,3);
            else                legendPurity    = GetAndSetLegend(0.18,0.15,2);
            legendPurity->AddEntry(histoGammaPurity_Pt,"purity");
            legendPurity->AddEntry(histoGammaTruePurity_Pt,"rescaled purity, secondaries removed");
            if (doPileUpCorr) legendPurity->AddEntry(histoGammaTruePurity_PileUp_Pt,"rescaled purity, secondaries removed, pileup");
            legendPurity->Draw();

            PutProcessLabelAndEnergyOnPlot( 0.18, 0.5, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
            canvasPurity->SaveAs(Form("%s/%s_Purity_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        }
        if ( isPCM && isCalo ){
            DrawGammaSetMarker(histoGammaCaloPurity_Pt, 24, 1.0, kRed+2, kRed+2);
            DrawGammaSetMarker(histoGammaCaloTruePurity_Pt, 20, 1.0, 1, 1);
            
            SetHistogramm(histoGammaCaloTruePurity_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{pur,#gamma} in |#eta| < %g",etaCalo),0.5, 1.);
            histoGammaCaloTruePurity_Pt->Draw();
            histoGammaCaloPurity_Pt->Draw("same");
            
            legendPurity                        = GetAndSetLegend(0.45,0.15,2);
            legendPurity->AddEntry(histoGammaPurity_Pt,"purity");
            legendPurity->AddEntry(histoGammaTruePurity_Pt,"rescaled purity, secondaries removed");
            legendPurity->Draw();

            PutProcessLabelAndEnergyOnPlot( 0.7, 0.4, 0.035, cent, textMeasurement, detectionProcess2, 42, 0.03);
        
            canvasPurity->SaveAs(Form("%s/%s_PurityCalo_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        }
    
    delete canvasPurity;

    TCanvas*    canvasPurity2                   = GetAndSetCanvas("canvasPurity2");
    Bool_t      doPurityPlotSimple              = kTRUE;
    if (doPurityPlotSimple && ((isPCM && !isCalo) || (isCalo && !isPCM))) {
        
        if (isPCM && !isCalo) SetHistogramm(histoGammaTruePurity_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{pur,#gamma} in |#eta| < %g",eta),0.8, 1.1);
        if (isCalo && !isPCM) SetHistogramm(histoGammaTruePurity_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{pur,#gamma} in |#eta| < %g",etaCalo),0.8, 1.1);
        
        histoGammaTruePurity_Pt->Draw();
        DrawGammaSetMarker(histoGammaTruePurity_Pt, 24, 1.5, 1, 1);
        DrawGammaLines(0., maxPtGamma,1.0, 1.0, 1, kGray+2, 2);
        histoGammaTruePurity_Pt->Draw("same");
        
        PutProcessLabelAndEnergyOnPlot( 0.18, 0.3, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
        canvasPurity2->SaveAs(Form("%s/%s_TruePurity_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    }
    delete canvasPurity2;

    //**********************************************************************************
    //******************** Conversion Prob Plot ****************************************
    //**********************************************************************************
    if ( isPCM ){
        TCanvas *canvasConvProb         = GetAndSetCanvas("canvasConvProb");
            DrawGammaSetMarker(histoGammaConvProb_MCPt, 24, 2.0, 1, 1);
            SetHistogramm(histoGammaConvProb_MCPt, "#it{p}_{T} (GeV/#it{c})",Form("#it{P}_{conv} in |#eta| < %g",eta), 0.04, 0.10);
            histoGammaConvProb_MCPt->Draw();

            TF1 *fConv                  = new TF1("line","[0]",2.5,25.);
            histoGammaConvProb_MCPt->Fit(fConv,"QRME0");
            Double_t parameterProb[1];
            fConv->GetParameters(parameterProb);

            DrawGammaLines(0., maxPtGamma,parameterProb[0], parameterProb[0], 1, kGray+2, 2);
            histoGammaConvProb_MCPt->Draw("same");

            PutProcessLabelAndEnergyOnPlot( 0.75, 0.3, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);

        canvasConvProb->SaveAs(Form("%s/%s_ConversionProb_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasConvProb;

        // secondary reco eff plot (only used for cocktail sec corr)
        if (hasCocktailInput) {
            TCanvas *canvasConvProbSec  = GetAndSetCanvas("canvasConvProbSec");
            canvasConvProbSec->SetTopMargin(0.035);
            
            DrawGammaSetMarker(histoGammaSecondaryFromXFromK0sConvProb_MCPt, 20, 1.0, 1, 1);
            DrawGammaSetMarker(histoGammaSecondaryFromXFromK0lConvProb_MCPt, 24, 1.0, kBlue-4, kBlue-4);
            DrawGammaSetMarker(histoGammaSecondaryFromXFromLambdaConvProb_MCPt, 34, 1.0, kCyan+2, kCyan+2);
            
            SetHistogramm(histoGammaSecondaryFromXFromK0sConvProb_MCPt, "#it{p}_{T} (GeV/#it{c})",Form("#it{P}_{conv} in |#eta| < %g",eta), 0., 0.09);
            histoGammaSecondaryFromXFromK0sConvProb_MCPt->Draw();
            histoGammaSecondaryFromXFromK0lConvProb_MCPt->Draw("same");
            histoGammaSecondaryFromXFromLambdaConvProb_MCPt->Draw("same");
            
            TF1 *fConvK0s               = new TF1("line","[0]",2.5,25.);
            TF1 *fConvK0l               = new TF1("line","[0]",2.5,25.);
            TF1 *fConvLambda            = new TF1("line","[0]",2.5,25.);
            histoGammaSecondaryFromXFromK0sConvProb_MCPt->Fit(fConvK0s,"QRME0");
            histoGammaSecondaryFromXFromK0lConvProb_MCPt->Fit(fConvK0l,"QRME0");
            histoGammaSecondaryFromXFromLambdaConvProb_MCPt->Fit(fConvLambda,"QRME0");
            Double_t parameterProbK0s[1];
            Double_t parameterProbK0l[1];
            Double_t parameterProbLambda[1];
            fConvK0s->GetParameters(parameterProbK0s);
            fConvK0l->GetParameters(parameterProbK0l);
            fConvLambda->GetParameters(parameterProbLambda);
            
            // not flat (K0s) or not enough stat (K0l + Lambda)
            //DrawGammaLines(0., maxPtGamma,parameterProbK0s[0], parameterProbK0s[0],1,1,2);
            //DrawGammaLines(0., maxPtGamma,parameterProbK0l[0], parameterProbK0l[0],1,kBlue-4,2);
            //DrawGammaLines(0., maxPtGamma,parameterProbLambda[0], parameterProbLambda[0],1,kCyan+2,2);
            
            PutProcessLabelAndEnergyOnPlot( 0.75, 0.3, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
            TLegend* legendSecConvProb = GetAndSetLegend(0.55,0.75,4);
            legendSecConvProb->AddEntry(histoGammaSecondaryFromXFromK0sConvProb_MCPt,"Secondary #gamma from K^{0}_{S}","p");
            legendSecConvProb->AddEntry(histoGammaSecondaryFromXFromK0lConvProb_MCPt,"Secondary #gamma from K^{0}_{L}","p");
            legendSecConvProb->AddEntry(histoGammaSecondaryFromXFromLambdaConvProb_MCPt,"Secondary #gamma from #Lambda","p");
            legendSecConvProb->Draw();
            
            canvasConvProbSec->SaveAs(Form("%s/%s_ConversionProbSecondaries_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
            delete canvasConvProbSec;
        }
    }

    //**********************************************************************************
    //******************** Reconstruction Eff Plot *************************************
    //**********************************************************************************
    TCanvas *canvasRecoEff = GetAndSetCanvas("canvasRecoEff");
        if ( isPCM && !isCalo ){
            DrawGammaSetMarker(histoGammaPrimaryRecoEff_MCPt, 20, 1.0, 1, 1);
            SetHistogramm(histoGammaPrimaryRecoEff_MCPt,"#it{p}_{T,MC} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);
            histoGammaPrimaryRecoEff_MCPt->Draw();
        
            PutProcessLabelAndEnergyOnPlot( 0.75, 0.95, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
            canvasRecoEff->SaveAs(Form("%s/%s_ReconstructionEff_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        }
        if ( isCalo && !isPCM ){
            DrawGammaSetMarker(histoGammaPrimaryRecoEff_MCPt, 20, 1.0, 1, 1);
            SetHistogramm(histoGammaPrimaryRecoEff_MCPt,"#it{p}_{T,MC} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
            histoGammaPrimaryRecoEff_MCPt->Draw();
        
            PutProcessLabelAndEnergyOnPlot( 0.75, 0.95, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
            canvasRecoEff->SaveAs(Form("%s/%s_ReconstructionEff_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        }
        if ( isPCM && isCalo ){
            DrawGammaSetMarker(histoGammaCaloPrimaryRecoEff_MCPt, 20, 1.0, 1, 1);
            SetHistogramm(histoGammaCaloPrimaryRecoEff_MCPt,"#it{p}_{T,MC} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
            histoGammaCaloPrimaryRecoEff_MCPt->Draw();
        
            PutProcessLabelAndEnergyOnPlot( 0.75, 0.95, 0.035, cent, textMeasurement, detectionProcess2, 42, 0.03);
        
            canvasRecoEff->SaveAs(Form("%s/%s_ReconstructionEffCalo_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        
        }    
    delete canvasRecoEff;

    // secondary reco eff plot (only used for cocktail sec corr)
    if (hasCocktailInput) {
        TCanvas *canvasRecoEffSec       = GetAndSetCanvas("canvasRecoEffSec");
        if ( isPCM || isCalo ){
            DrawGammaSetMarker(histoGammaSecondaryFromXFromK0sRecoEff_Pt, 20, 1.0, 1, 1);
            DrawGammaSetMarker(histoGammaSecondaryFromXFromK0lRecoEff_Pt, 24, 1.0, kBlue-4, kBlue-4);
            DrawGammaSetMarker(histoGammaSecondaryFromXFromLambdaRecoEff_Pt, 34, 1.0, kCyan+2, kCyan+2);
            
            if (isPCM && !isCalo) SetHistogramm(histoGammaSecondaryFromXFromK0lRecoEff_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);
            if (isCalo && !isPCM) SetHistogramm(histoGammaSecondaryFromXFromK0lRecoEff_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
            
            histoGammaSecondaryFromXFromK0lRecoEff_Pt->Draw();
            histoGammaSecondaryFromXFromK0sRecoEff_Pt->Draw("same");
            histoGammaSecondaryFromXFromLambdaRecoEff_Pt->Draw("same");
            
            PutProcessLabelAndEnergyOnPlot( 0.75, 0.95, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
            TLegend* legendSecRecoEff   = GetAndSetLegend(0.40,0.76,4);
            legendSecRecoEff->AddEntry(histoGammaSecondaryFromXFromK0sRecoEff_Pt,"Secondary #gamma from K^{0}_{S}","p");
            legendSecRecoEff->AddEntry(histoGammaSecondaryFromXFromK0lRecoEff_Pt,"Secondary #gamma from K^{0}_{L}","p");
            legendSecRecoEff->AddEntry(histoGammaSecondaryFromXFromLambdaRecoEff_Pt,"Secondary #gamma from #Lambda","p");
            legendSecRecoEff->Draw();
        }
        
        canvasRecoEffSec->SaveAs(Form("%s/%s_ReconstructionEffSec_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasRecoEffSec;
    }

    //**********************************************************************************
    //******************** Reconstruction Eff Comparison Plot **************************
    //**********************************************************************************
    TCanvas* canvasCompRecoEff = GetAndSetCanvas("canvasCompRecoEff");
        DrawGammaCanvasSettings( canvasCompRecoEff,0.1, 0.015, 0.02, 0.09);
        TPad* padCompRecoEff = new TPad("padCompRecoEff", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padCompRecoEff, 0.10, 0.015, 0.02, 0.);
        padCompRecoEff->Draw();
        TPad* padBinCompRecoEffRatio = new TPad("padBinCompRecoEffRatio", "", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings( padBinCompRecoEffRatio,  0.1, 0.015, 0.0, 0.3);
        padBinCompRecoEffRatio->Draw();

        if ( (isPCM && !isCalo) || (isCalo && !isPCM) ){
            padCompRecoEff->cd();
            
            if (isPCM && !isCalo) {
                SetHistogramm(histoGammaPrimaryRecoEff_MCPt, "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);
                SetHistogramm(histoGammaPrimaryRecoEff_Pt, "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0.05, 1.0);
                if(doPileUpCorr)SetHistogramm(histoGammaPrimaryRecoEff_PileUp_Pt, "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);
            }
            if (isCalo && !isPCM) {
                SetHistogramm(histoGammaPrimaryRecoEff_MCPt, "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
                SetHistogramm(histoGammaPrimaryRecoEff_Pt, "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0.05, 1.0);
                if(doPileUpCorr)SetHistogramm(histoGammaPrimaryRecoEff_PileUp_Pt, "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
            }

            DrawGammaSetMarker(histoGammaPrimaryRecoEff_MCPt, 20, 1.0, kBlack, kBlack);
            DrawGammaSetMarker(histoGammaPrimaryRecoEff_Pt, 24, 1.0, kBlue+1, kBlue+1);
            if(doPileUpCorr)DrawGammaSetMarker(histoGammaPrimaryRecoEff_PileUp_Pt, 20, 1.0, 5, 5);

            histoGammaPrimaryRecoEff_Pt->Draw("e1");
            if(doPileUpCorr)histoGammaPrimaryRecoEff_PileUp_Pt->Draw("e1,same");
            histoGammaPrimaryRecoEff_MCPt->Draw("same,e1");

            TLegend* legendCompRecoEff;
            if(doPileUpCorr)
                legendCompRecoEff = GetAndSetLegend(0.15 ,0.05,3,1);
            else
                legendCompRecoEff = GetAndSetLegend(0.15 ,0.05,2,1);
            legendCompRecoEff->SetTextSize(0.04);
            legendCompRecoEff->AddEntry(histoGammaPrimaryRecoEff_MCPt,"MC pT Reconstruction Efficiency (unfolding)");
            legendCompRecoEff->AddEntry(histoGammaPrimaryRecoEff_Pt,"Reconstruction Efficiency for primary  #gamma");
            if(doPileUpCorr)legendCompRecoEff->AddEntry(histoGammaPrimaryRecoEff_PileUp_Pt,"Reconstruction Efficiency for primary  #gamma Pile Up");
            legendCompRecoEff->Draw();
        
            PutProcessLabelAndEnergyOnPlot( 0.75, 0.95, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        }

        TH1D* histoGammaResolCorrEff_Pt = NULL;
        if ( (isPCM && !isCalo) || (isCalo && !isPCM) ) {
            histoGammaResolCorrEff_Pt = (TH1D*)histoGammaPrimaryRecoEff_Pt->Clone("histoGammaResolCorrEff_Pt");
            histoGammaResolCorrEff_Pt->Divide(histoGammaPrimaryRecoEff_MCPt);

            padBinCompRecoEffRatio->cd();
            TH1D* histoGammaPrimaryRecoEff_PileUp_PtRatio;
            if(doPileUpCorr){
                histoGammaPrimaryRecoEff_PileUp_PtRatio = (TH1D*)histoGammaPrimaryRecoEff_PileUp_Pt->Clone("GammaPrimaryRecoEffRatioPileUp");
                histoGammaPrimaryRecoEff_PileUp_PtRatio->Divide(histoGammaPrimaryRecoEff_MCPt);
            }
            SetStyleHistoTH1ForGraphs(  histoGammaResolCorrEff_Pt, "#it{p}_{T} (GeV/#it{c})","#frac{Reco Eff x}{Primary MC Reco Eff}" , 0.14, 0.15,
                                        0.12, 0.10,  0.85, 0.4, 510, 505);
            histoGammaResolCorrEff_Pt->GetYaxis()->SetRangeUser(0.7,1.25);

            histoGammaResolCorrEff_Pt->DrawCopy("e1");
            if(doPileUpCorr)histoGammaPrimaryRecoEff_PileUp_PtRatio->DrawCopy("e1,same");
        
            DrawGammaLines(0., maxPtGamma,1, 1,0.5,kGray+2,2);
            
            canvasCompRecoEff ->SaveAs(Form("%s/%s_CompRecEff_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        }
        if (isPCM && isCalo){
            padCompRecoEff->cd();
            SetHistogramm(histoGammaCaloPrimaryRecoEff_MCPt, "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
            SetHistogramm(histoGammaCaloPrimaryRecoEff_Pt, "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0.05, 1.0);
        
            DrawGammaSetMarker(histoGammaCaloPrimaryRecoEff_MCPt, 20, 1.0, kBlack, kBlack);
            DrawGammaSetMarker(histoGammaCaloPrimaryRecoEff_Pt, 24, 1.0, kBlue+1, kBlue+1);
        
            histoGammaCaloPrimaryRecoEff_Pt->Draw("e1");
            histoGammaCaloPrimaryRecoEff_MCPt->Draw("same,e1");

            TLegend* legendCompRecoEff = GetAndSetLegend(0.15 ,0.05,3,1);
            legendCompRecoEff->SetTextSize(0.04);
            legendCompRecoEff->AddEntry(histoGammaCaloPrimaryRecoEff_MCPt,"MC pT Reconstruction Efficiency (unfolding)");
            legendCompRecoEff->AddEntry(histoGammaCaloPrimaryRecoEff_Pt,"Reconstruction Efficiency for primary  #gamma");
            legendCompRecoEff->Draw();
        
            PutProcessLabelAndEnergyOnPlot( 0.75, 0.95, 0.035, cent, textMeasurement, detectionProcess2, 42, 0.03);
        }

        TH1D* histoGammaCaloResolCorrEff_Pt = NULL;
        if (isPCM && isCalo){
            histoGammaCaloResolCorrEff_Pt = (TH1D*)histoGammaCaloPrimaryRecoEff_Pt->Clone("histoGammaCaloResolCorrEff_Pt");
            histoGammaCaloResolCorrEff_Pt->Divide(histoGammaCaloPrimaryRecoEff_MCPt);

            padBinCompRecoEffRatio->cd();
            SetStyleHistoTH1ForGraphs(  histoGammaCaloResolCorrEff_Pt, "#it{p}_{T} (GeV/#it{c})","#frac{Reco Eff x}{Primary MC Reco Eff}" , 0.14, 0.15,
                                        0.12, 0.10,  0.85, 0.4, 510, 505);
            histoGammaCaloResolCorrEff_Pt->GetYaxis()->SetRangeUser(0.7,1.25);

            histoGammaCaloResolCorrEff_Pt->DrawCopy("e1");
        
            DrawGammaLines(0., maxPtGamma,1, 1,0.5,kGray+2,2);
    
            canvasCompRecoEff ->SaveAs(Form("%s/%s_CompCaloRecEff_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        }
        
    delete padCompRecoEff;
    delete padBinCompRecoEffRatio;
    delete canvasCompRecoEff;

    //**********************************************************************************
    //******************** Response Matrix for detector resolution *********************
    //**********************************************************************************
    if (isPCM){
        SetStyleHistoTH2ForGraphs(  histoGammaTruePrimaryConv_recPt_MCPt, "Reconstructed #it{p}_{T} (GeV/#it{c})","MC #it{p}_{T} (GeV/#it{c})", 0.035, 0.04,  
                                    0.035, 0.04, 0.9, 1.0, 510, 510);
        histoGammaTruePrimaryConv_recPt_MCPt->GetYaxis()->SetRangeUser(0,25);
        histoGammaTruePrimaryConv_recPt_MCPt->GetXaxis()->SetRangeUser(0,25);
    }
    if (isCalo){
        SetStyleHistoTH2ForGraphs(  histoGammaTruePrimaryCalo_recPt_MCPt, "Reconstructed #it{p}_{T} (GeV/#it{c})","MC #it{p}_{T} (GeV/#it{c})", 0.035, 0.04,  
                                    0.035, 0.04, 0.9, 1.0, 510, 510);
        histoGammaTruePrimaryCalo_recPt_MCPt->GetYaxis()->SetRangeUser(0,25);
        histoGammaTruePrimaryCalo_recPt_MCPt->GetXaxis()->SetRangeUser(0,25);
    }

    //**********************************************************************************
    //******************** Response Matrix Plot ****************************************
    //**********************************************************************************
    TCanvas * canvasResponseMatrix = new TCanvas("canvasResponseMatrix","",480,440);  // gives the page size
    DrawGammaCanvasSettings( canvasResponseMatrix, 0.09, 0.105, 0.02, 0.085);    
    canvasResponseMatrix->SetLogz(1);
    canvasResponseMatrix->cd();

        if (isPCM){
            histoGammaTruePrimaryConv_recPt_MCPt->Draw("colz");
            PutProcessLabelAndEnergyOnPlot( 0.18, 0.95, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
            canvasResponseMatrix->SaveAs(Form("%s/%s_ResponseMatrix_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        } 
        
        if (isCalo){
            histoGammaTruePrimaryCalo_recPt_MCPt->Draw("colz");
            PutProcessLabelAndEnergyOnPlot( 0.18, 0.95, 0.035, cent, textMeasurement, detectionProcess2, 42, 0.03);
            canvasResponseMatrix->SaveAs(Form("%s/%s_ResponseMatrixCalo_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
            
        }
    //     delete canvasResponseMatrix;

    // response matrix for cocktail secondary corr
    if (hasCocktailInput) {
        TCanvas* canvasResponseMatrixSec = new TCanvas("canvasResponseMatrixSec","",480,440);  // gives the page size
        DrawGammaCanvasSettings( canvasResponseMatrixSec, 0.09, 0.105, 0.02, 0.085);
        canvasResponseMatrixSec->SetLogz(1);
        canvasResponseMatrixSec->cd();
        
        if (isPCM || isCalo){
            SetStyleHistoTH2ForGraphs(  histoGammaTrueSecondaryFromXFromK0s_MCPt_recPtOrBin, "Reconstructed #it{p}_{T} (GeV/#it{c})","MC #it{p}_{T} (GeV/#it{c})", 0.035, 0.04,
                                      0.035, 0.04, 0.9, 1.0, 510, 510);
            histoGammaTrueSecondaryFromXFromK0s_MCPt_recPtOrBin->GetYaxis()->SetRangeUser(0,25);
            histoGammaTrueSecondaryFromXFromK0s_MCPt_recPtOrBin->GetXaxis()->SetRangeUser(0,25);
            histoGammaTrueSecondaryFromXFromK0s_MCPt_recPtOrBin->Draw("colz");
            PutProcessLabelAndEnergyOnPlot( 0.18, 0.95, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        }

        canvasResponseMatrixSec->SaveAs(Form("%s/%s_ResponseMatrixSecK0s_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    }

    //**********************************************************************************
    //************************ Unfolding of inclusive gamma spectrum *******************
    //**********************************************************************************
    TH1D* histoGammaCorrUnfoldReso_Pt                           = NULL;
    TH1D* histoGammaCorrUnfoldReso_BinByBin_Pt                  = NULL;
    TH1D* histoGammaResolCorrUnfold_Pt                          = NULL;
    TH1D* histoGammaCorrUnfoldReso_PtNotCorrected               = NULL;
    TH1D* histoGammaResolCorrUnfold_BinByBin_Pt                 = NULL;
    TH1D* histoSecondaryGammaSpecPtOriginalBin                  = NULL;
    TH1D* histoSecondaryGammaFromXFromK0sSpecPtOriginalBin      = NULL;
    // do the same with MC rec gammas as a sanity check
    TH1D* histoMCrecGammaCorr_Pt                                = NULL;

    if (isPCM && !isCalo){
        // Correct inclusive photon spectrum for pileup contribution
        if(doPileUpCorr) histoESDConvGammaPt_OriginalBin->Multiply(histoPileUpCorrectionFactorOrBin);
        
        // correct measured gamma spectra for secondaries and purity
        if(!hasCocktailInput){
            // determine secondary contribution in orginal binning
            histoSecondaryGammaSpecPtOriginalBin                = (TH1D*)histoESDConvGammaPt_OriginalBin->Clone("SecondaryGammaSpecPt");
            histoSecondaryGammaFromXFromK0sSpecPtOriginalBin    = (TH1D*)histoESDConvGammaPt_OriginalBin->Clone("SecondaryGammaSpecFromXFromK0sPt");
            histoSecondaryGammaSpecPtOriginalBin->Multiply(histoFracAllGammaToSec_OriginalBin_Pt);
            histoSecondaryGammaFromXFromK0sSpecPtOriginalBin->Multiply(histoFracAllGammaToSecFromXFromK0s_OriginalBin_Pt);
            histoSecondaryGammaFromXFromK0sSpecPtOriginalBin->Scale(doubleAddFactorK0s);

            //subtract secondary contribution from inclusive spectrum
            CorrectGammaSecAndPurity(histoESDConvGammaPt_OriginalBin, histoSecondaryGammaSpecPtOriginalBin, histoSecondaryGammaFromXFromK0sSpecPtOriginalBin, histoGammaTruePurity_OriginalBin_Pt );
            CorrectGammaSecAndPurity(histoMCrecGamma_OriginalBin_Pt, histoSecondaryGammaSpecPtOriginalBin, histoSecondaryGammaFromXFromK0sSpecPtOriginalBin, histoGammaTruePurity_OriginalBin_Pt );
        } else {
            //subtract secondary contribution from cocktail from inclusive spectrum
            cout << "will use cocktail for secondary correction" << endl;
            CorrectGammaSecAndPurityCocktail(
                histoESDConvGammaPt_OriginalBin,
                histoGammaSecGammaFromXFromK0s_Cocktail_Raw_PtOrBin,
                histoGammaSecGammaFromXFromK0l_Cocktail_Raw_PtOrBin,
                histoGammaSecGammaFromXFromLambda_Cocktail_Raw_PtOrBin,
                histoGammaTrueSecCocktailGammaRest_PtOrBin,
                histoGammaTruePurity_OriginalBin_Pt );
            CorrectGammaSecAndPurityCocktail(
                histoMCrecGamma_OriginalBin_Pt,
                histoGammaSecGammaFromXFromK0s_Cocktail_Raw_PtOrBin,
                histoGammaSecGammaFromXFromK0l_Cocktail_Raw_PtOrBin,
                histoGammaSecGammaFromXFromLambda_Cocktail_Raw_PtOrBin,
                histoGammaTrueSecCocktailGammaRest_PtOrBin,
                histoGammaTruePurity_OriginalBin_Pt );
        }
        
        // create histograms for unfolding for different techniques
        histoGammaCorrUnfoldReso_Pt                             = (TH1D*)histoESDConvGammaPt_OriginalBin->Clone("histoGammaCorrUnfoldReso_Pt");
        histoMCrecGammaCorr_Pt                                  = (TH1D*)histoMCrecGamma_OriginalBin_Pt->Clone("GammaSpecCorrESDMC");
        histoGammaCorrUnfoldReso_BinByBin_Pt                    = (TH1D*)histoESDConvGammaPt_OriginalBin->Clone("histoGammaCorrUnfoldReso_BinByBin_Pt");
        // TH1D *histoGammaCorrUnfoldReso_SvD_Pt                = (TH1D*)histoESDConvGammaPt_OriginalBin->Clone("histoGammaCorrUnfoldReso_SvD_Pt");
        // TH1D *histoGammaCorrUnfoldReso_TUnfold_Pt            = (TH1D*)histoESDConvGammaPt_OriginalBin->Clone("histoGammaCorrUnfoldReso_TUnfold_Pt");

        // Unfolding setup:
        // RooUnfoldResponse constructor - create from already-filled histograms
        // "response" gives the response matrix, measured X truth.
        // "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
        // but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
        // in "truth" for unmeasured events (inefficiency).
        // "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
        // to indicate, respectively, no fakes and/or no inefficiency.
        RooUnfoldResponse   response(0,0,histoGammaTruePrimaryConv_recPt_MCPt);
        
        // unfold with different techniques, but same response matrix
        RooUnfoldBayes      unfold_Spectrum (&response,histoGammaCorrUnfoldReso_Pt, nIterationsUnfolding);
        RooUnfoldBayes      unfold_SpectrumMCrec (&response,histoMCrecGammaCorr_Pt, nIterationsUnfolding);
        RooUnfoldBinByBin   unfold_SpectrumBinByBin (&response,histoGammaCorrUnfoldReso_BinByBin_Pt);
        //RooUnfoldSvd      unfold_SpectrumSvD (&response, histoGammaCorrUnfoldReso_SvD_Pt, 20);
        //RooUnfoldTUnfold  unfold_SpectrumTUnfold (&response,histoGammaCorrUnfoldReso_TUnfold_Pt);
        
        // get histograms from RooUnfold and rebin them 
        histoGammaCorrUnfoldReso_Pt                             = (TH1D*)unfold_Spectrum.Hreco();
        histoMCrecGammaCorr_Pt                                  = (TH1D*)unfold_SpectrumMCrec.Hreco();
        histoMCrecGammaCorr_Pt                                  = RebinTH1D(histoMCrecGammaCorr_Pt,histoESDConvGammaPt,kTRUE);
        histoGammaResolCorrUnfold_Pt                            = (TH1D*)histoESDConvGammaPt_OriginalBin->Clone("histoGammaResolCorrUnfold_Pt");
        histoGammaResolCorrUnfold_Pt->Divide(histoGammaCorrUnfoldReso_Pt);    
        histoGammaCorrUnfoldReso_Pt                             = RebinTH1D(histoGammaCorrUnfoldReso_Pt,histoESDConvGammaPt,kTRUE);

        histoGammaCorrUnfoldReso_BinByBin_Pt                    = (TH1D*)unfold_SpectrumBinByBin.Hreco();
        histoGammaResolCorrUnfold_BinByBin_Pt                   = (TH1D*)histoESDConvGammaPt_OriginalBin->Clone("histoGammaResolCorrUnfold_BinByBin_Pt");
        histoGammaResolCorrUnfold_BinByBin_Pt->Divide(histoGammaCorrUnfoldReso_BinByBin_Pt);

        histoGammaCorrUnfoldReso_BinByBin_Pt                    = RebinTH1D(histoGammaCorrUnfoldReso_BinByBin_Pt,histoESDConvGammaPt,kTRUE);
        // histoGammaCorrUnfoldReso_SvD_Pt                      = (TH1D*)unfold_SpectrumSvD.Hreco();
        // histoGammaCorrUnfoldReso_SvD_Pt                      = RebinTH1D(histoGammaCorrUnfoldReso_SvD_Pt,histoESDConvGammaPt,kTRUE);
        // histoGammaCorrUnfoldReso_TUnfold_Pt                  = (TH1D*)unfold_SpectrumTUnfold.Hreco();
        // histoGammaCorrUnfoldReso_TUnfold_Pt                  = RebinTH1D(histoGammaCorrUnfoldReso_TUnfold_Pt,histoESDConvGammaPt,kTRUE);

        // Correct inclusive photon spectrum with conversion probability & reco efficiency (both versus MC pt)
        histoGammaCorrUnfoldReso_PtNotCorrected                 = (TH1D*)histoGammaCorrUnfoldReso_Pt->Clone("GammaUnfoldNotCorrected");
        CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
        CorrectGammaUnfoldResol(histoMCrecGammaCorr_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
        CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_BinByBin_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
        //CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_SvD_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
        //CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_TUnfold_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
        
        SetHistogramm(histoGammaCorrUnfoldReso_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        SetHistogramm(histoGammaCorrUnfoldReso_BinByBin_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        SetHistogramm(histoMCrecGammaCorr_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        DrawGammaSetMarker(histoGammaCorrUnfoldReso_Pt, 20, 1.0, 1, 1);
        DrawGammaSetMarker(histoGammaCorrUnfoldReso_BinByBin_Pt, 24, 1.0, kBlue, kBlue);
        DrawGammaSetMarker(histoMCrecGammaCorr_Pt, 20, 1.0, kGreen-1, kGreen-1);
        //SetHistogramm(histoGammaCorrUnfoldReso_SvD_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        //SetHistogramm(histoGammaCorrUnfoldReso_TUnfold_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
    }

    if (isCalo && !isPCM) {
        // correct measured gamma spectra for secondaries and purity
        if(!hasCocktailInput){
            // determine secondary contribution in orginal binning
            TH1D *histoSecondaryGammaSpecPtOriginalBin              = (TH1D*)histoESDCaloGammaPt_OriginalBin->Clone("SecondaryGammaSpecPt");
            TH1D *histoSecondaryGammaFromXFromK0sSpecPtOriginalBin  = (TH1D*)histoESDCaloGammaPt_OriginalBin->Clone("SecondaryGammaSpecFromXFromK0sPt");
            histoSecondaryGammaSpecPtOriginalBin->Multiply(histoFracAllGammaToSec_OriginalBin_Pt);
            histoSecondaryGammaFromXFromK0sSpecPtOriginalBin->Multiply(histoFracAllGammaToSecFromXFromK0s_OriginalBin_Pt);
            histoSecondaryGammaFromXFromK0sSpecPtOriginalBin->Scale(doubleAddFactorK0s);
            
            //subtract secondary contribution from inclusive spectrum
            CorrectGammaSecAndPurity(histoESDCaloGammaPt_OriginalBin, histoSecondaryGammaSpecPtOriginalBin, histoSecondaryGammaFromXFromK0sSpecPtOriginalBin, histoGammaTruePurity_OriginalBin_Pt );
            CorrectGammaSecAndPurity(histoMCrecGammaCalo_OriginalBin_Pt, histoSecondaryGammaSpecPtOriginalBin, histoSecondaryGammaFromXFromK0sSpecPtOriginalBin, histoGammaTruePurity_OriginalBin_Pt );
        } else {
            //subtract secondary contribution from cocktail from inclusive spectrum
            cout << "will use cocktail for secondary correction" << endl;
            CorrectGammaSecAndPurityCocktail(
                                             histoESDCaloGammaPt_OriginalBin,
                                             histoGammaSecGammaFromXFromK0s_Cocktail_Raw_PtOrBin,
                                             histoGammaSecGammaFromXFromK0l_Cocktail_Raw_PtOrBin,
                                             histoGammaSecGammaFromXFromLambda_Cocktail_Raw_PtOrBin,
                                             histoGammaTrueSecCocktailGammaRest_PtOrBin,
                                             histoGammaTruePurity_OriginalBin_Pt );
            CorrectGammaSecAndPurityCocktail(
                                             histoMCrecGammaCalo_OriginalBin_Pt,
                                             histoGammaSecGammaFromXFromK0s_Cocktail_Raw_PtOrBin,
                                             histoGammaSecGammaFromXFromK0l_Cocktail_Raw_PtOrBin,
                                             histoGammaSecGammaFromXFromLambda_Cocktail_Raw_PtOrBin,
                                             histoGammaTrueSecCocktailGammaRest_PtOrBin,
                                             histoGammaTruePurity_OriginalBin_Pt );
        }

        // create histograms for unfolding for different techniques
        histoGammaCorrUnfoldReso_Pt                                 = (TH1D*)histoESDCaloGammaPt_OriginalBin->Clone("histoGammaCorrUnfoldReso_Pt");
        histoGammaCorrUnfoldReso_BinByBin_Pt                        = (TH1D*)histoESDCaloGammaPt_OriginalBin->Clone("histoGammaCorrUnfoldReso_BinByBin_Pt");
        histoMCrecGammaCorr_Pt                                      = (TH1D*)histoMCrecGammaCalo_OriginalBin_Pt->Clone("GammaSpecCorrESDMC");
        
        // Unfolding setup:
        // RooUnfoldResponse constructor - create from already-filled histograms
        // "response" gives the response matrix, measured X truth.
        // "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
        // but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
        // in "truth" for unmeasured events (inefficiency).
        // "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
        // to indicate, respectively, no fakes and/or no inefficiency.
        RooUnfoldResponse   response(0,0,histoGammaTruePrimaryCalo_recPt_MCPt);

        // unfold with different techniques, but same response matrix
        RooUnfoldBayes      unfold_Spectrum (&response,histoGammaCorrUnfoldReso_Pt, nIterationsUnfolding);
        RooUnfoldBayes      unfold_SpectrumMCrec (&response,histoMCrecGammaCorr_Pt, nIterationsUnfolding);
        RooUnfoldBinByBin   unfold_SpectrumBinByBin (&response,histoGammaCorrUnfoldReso_BinByBin_Pt);

        // get histograms from RooUnfold and rebin them
        histoGammaCorrUnfoldReso_Pt                                 = (TH1D*)unfold_Spectrum.Hreco();
        histoMCrecGammaCorr_Pt                                      = (TH1D*)unfold_SpectrumMCrec.Hreco();

        histoMCrecGammaCorr_Pt                                      = RebinTH1D(histoMCrecGammaCorr_Pt,histoESDCaloGammaPt,kTRUE);
        histoGammaResolCorrUnfold_Pt                                = (TH1D*)histoESDCaloGammaPt_OriginalBin->Clone("histoGammaResolCorrUnfold_Pt");
        histoGammaResolCorrUnfold_Pt->Divide(histoGammaCorrUnfoldReso_Pt);
        histoGammaCorrUnfoldReso_Pt                                 = RebinTH1D(histoGammaCorrUnfoldReso_Pt,histoESDCaloGammaPt,kTRUE);

        histoGammaCorrUnfoldReso_BinByBin_Pt                        = (TH1D*)unfold_SpectrumBinByBin.Hreco();
        histoGammaResolCorrUnfold_BinByBin_Pt                       = (TH1D*)histoESDCaloGammaPt_OriginalBin->Clone("histoGammaResolCorrUnfold_BinByBin_Pt");
        histoGammaResolCorrUnfold_BinByBin_Pt->Divide(histoGammaCorrUnfoldReso_BinByBin_Pt);
        histoGammaCorrUnfoldReso_BinByBin_Pt                        = RebinTH1D(histoGammaCorrUnfoldReso_BinByBin_Pt,histoESDCaloGammaPt,kTRUE);

        // Correct inclusive photon spectrum with conversion probability & reco efficiency (both versus MC pt)
        histoGammaCorrUnfoldReso_PtNotCorrected                     = (TH1D*)histoGammaCorrUnfoldReso_Pt->Clone("GammaUnfoldNotCorrected");
        CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_Pt, histoGammaPrimaryRecoEff_MCPt, deltaEtaCalo, scalingCalo, nEvt);
        CorrectGammaUnfoldResol(histoMCrecGammaCorr_Pt, histoGammaPrimaryRecoEff_MCPt, deltaEtaCalo, scalingCalo, nEvt);
        CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_BinByBin_Pt, histoGammaPrimaryRecoEff_MCPt, deltaEtaCalo, scalingCalo, nEvt);

        SetHistogramm(histoGammaCorrUnfoldReso_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        SetHistogramm(histoGammaCorrUnfoldReso_BinByBin_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        SetHistogramm(histoMCrecGammaCorr_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        DrawGammaSetMarker(histoGammaCorrUnfoldReso_Pt, 20, 1.0, 1, 1);
        DrawGammaSetMarker(histoGammaCorrUnfoldReso_BinByBin_Pt, 24, 1.0, kBlue, kBlue);
        DrawGammaSetMarker(histoMCrecGammaCorr_Pt, 20, 1.0, kGreen-1, kGreen-1);
    }

    TH1D* histoGammaCaloCorrUnfoldReso_Pt                           = NULL;
    TH1D* histoGammaCaloCorrUnfoldReso_BinByBin_Pt                  = NULL;
    TH1D* histoGammaCaloCorrUnfoldResoPileUp_Pt                     = NULL;
    TH1D* histoGammaCaloResolCorrUnfold_Pt                          = NULL;
    TH1D* histoGammaCaloCorrUnfoldReso_PtNotCorrected               = NULL;
    TH1D* histoGammaCaloResolCorrUnfold_BinByBin_Pt                 = NULL;
    
//     if (isPCM && isCalo){
//         // determine secondary contribution in orginal binning
//         TH1D *histoSecondaryGammaSpecPtOriginalBin             = (TH1D*) histoESDConvGammaPt_OriginalBin->Clone("SecondaryGammaSpecPt");
//         TH1D *histoSecondaryGammaFromXFromK0sSpecPtOriginalBin = (TH1D*) histoESDConvGammaPt_OriginalBin->Clone("SecondaryGammaSpecFromXFromK0sPt");
//         histoSecondaryGammaSpecPtOriginalBin->Multiply(histoFracAllGammaToSec_OriginalBin_Pt);
//         histoSecondaryGammaFromXFromK0sSpecPtOriginalBin->Multiply(histoFracAllGammaToSecFromXFromK0s_OriginalBin_Pt);
//         histoSecondaryGammaFromXFromK0sSpecPtOriginalBin->Scale(doubleAddFactorK0s);
// 
//         //subtract secondary contribution from inclusive spectrum
//         CorrectGammaSecAndPurity(histoESDConvGammaPt_OriginalBin, histoSecondaryGammaSpecPtOriginalBin, histoSecondaryGammaFromXFromK0sSpecPtOriginalBin, histoGammaTruePurity_OriginalBin_Pt );
//         
//         
//         // create histograms for unfolding for different techniques
//         histoGammaCorrUnfoldReso_Pt             = (TH1D*) histoESDConvGammaPt_OriginalBin->Clone("histoGammaCorrUnfoldReso_Pt");
//         histoGammaCorrUnfoldReso_BinByBin_Pt    = (TH1D*) histoESDConvGammaPt_OriginalBin->Clone("histoGammaCorrUnfoldReso_BinByBin_Pt");
//         // TH1D *histoGammaCorrUnfoldReso_SvD_Pt       = (TH1D*) histoESDConvGammaPt_OriginalBin->Clone("histoGammaCorrUnfoldReso_SvD_Pt");
//         // TH1D *histoGammaCorrUnfoldReso_TUnfold_Pt   = (TH1D*) histoESDConvGammaPt_OriginalBin->Clone("histoGammaCorrUnfoldReso_TUnfold_Pt");
// 
//         // Unfolding setup:
//         // RooUnfoldResponse constructor - create from already-filled histograms
//         // "response" gives the response matrix, measured X truth.
//         // "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
//         // but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
//         // in "truth" for unmeasured events (inefficiency).
//         // "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
//         // to indicate, respectively, no fakes and/or no inefficiency.
//         RooUnfoldResponse   response(0,0,histoGammaTruePrimaryConv_recPt_MCPt);
//         
//         // unfold with different techniques, but same response matrix
//         RooUnfoldBayes      unfold_Spectrum (&response,histoGammaCorrUnfoldReso_Pt, 4);
//         RooUnfoldBinByBin   unfold_SpectrumBinByBin (&response,histoGammaCorrUnfoldReso_BinByBin_Pt);
//         //RooUnfoldSvd      unfold_SpectrumSvD (&response, histoGammaCorrUnfoldReso_SvD_Pt, 20);   
//         //RooUnfoldTUnfold  unfold_SpectrumTUnfold (&response,histoGammaCorrUnfoldReso_TUnfold_Pt);
//         
//         // get histograms from RooUnfold and rebin them 
//         histoGammaCorrUnfoldReso_Pt             = (TH1D*) unfold_Spectrum.Hreco();
//         histoGammaResolCorrUnfold_Pt            = (TH1D*) histoESDConvGammaPt_OriginalBin->Clone("histoGammaResolCorrUnfold_Pt");
//         histoGammaResolCorrUnfold_Pt->Divide(histoGammaCorrUnfoldReso_Pt);    
//         histoGammaCorrUnfoldReso_Pt             = RebinTH1D(histoGammaCorrUnfoldReso_Pt,histoESDConvGammaPt,kTRUE);
// 
//         histoGammaCorrUnfoldReso_BinByBin_Pt    = (TH1D*) unfold_SpectrumBinByBin.Hreco();
//         histoGammaResolCorrUnfold_BinByBin_Pt   = (TH1D*) histoESDConvGammaPt_OriginalBin->Clone("histoGammaResolCorrUnfold_BinByBin_Pt");
//         histoGammaResolCorrUnfold_BinByBin_Pt->Divide(histoGammaCorrUnfoldReso_BinByBin_Pt);
// 
//         histoGammaCorrUnfoldReso_BinByBin_Pt    = RebinTH1D(histoGammaCorrUnfoldReso_BinByBin_Pt,histoESDConvGammaPt,kTRUE);
//         // histoGammaCorrUnfoldReso_SvD_Pt      = (TH1D*) unfold_SpectrumSvD.Hreco();
//         // histoGammaCorrUnfoldReso_SvD_Pt      = RebinTH1D(histoGammaCorrUnfoldReso_SvD_Pt,histoESDConvGammaPt,kTRUE);
//         // histoGammaCorrUnfoldReso_TUnfold_Pt  = (TH1D*) unfold_SpectrumTUnfold.Hreco();
//         // histoGammaCorrUnfoldReso_TUnfold_Pt  = RebinTH1D(histoGammaCorrUnfoldReso_TUnfold_Pt,histoESDConvGammaPt,kTRUE);
// 
//         // Correct inclusive photon spectrum with conversion probability & reco efficiency (both versus MC pt)
//         histoGammaCorrUnfoldReso_PtNotCorrected = (TH1D*) histoGammaCorrUnfoldReso_Pt->Clone("GammaUnfoldNotCorrected");
//         CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
//         CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_BinByBin_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
//         //CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_SvD_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
//         //CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_TUnfold_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
//         
//         // Correct inclusive photon spectrum for pileup contribution
//         
//         SetHistogramm(histoGammaCorrUnfoldReso_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
//         DrawGammaSetMarker(histoGammaCorrUnfoldReso_Pt, 20, 1.0, 1, 1);
//         SetHistogramm(histoGammaCorrUnfoldReso_BinByBin_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
//         DrawGammaSetMarker(histoGammaCorrUnfoldReso_BinByBin_Pt, 24, 1.0, kBlue, kBlue);
//         //SetHistogramm(histoGammaCorrUnfoldReso_SvD_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
//         //SetHistogramm(histoGammaCorrUnfoldReso_TUnfold_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
//     }    
    
    //**********************************************************************************
    //******************** Corrected Photon Spectrum Plot ******************************
    //**********************************************************************************
    TCanvas *canvasCorrGammaSpecUnfold          = GetAndSetCanvas("canvasCorrGammaSpecUnfold");
    canvasCorrGammaSpecUnfold->SetLogy();
        TH1D* Dummy                             = NULL;
        if (isPCM && !isCalo) Dummy             = (TH1D*)histoESDConvGammaPt->Clone("Dummy");
        if (isCalo && !isPCM) Dummy             = (TH1D*)histoESDCaloGammaPt->Clone("Dummy");
    
        if (isPCM && !isCalo) SetHistogramm(Dummy,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",1e-8,10);
        if (isCalo && !isPCM) SetHistogramm(Dummy,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",1e-10,1e-2);

        Dummy->Reset();
        Dummy->Draw();
        histoGammaCorrUnfoldReso_Pt->DrawCopy("same");
        histoGammaCorrUnfoldReso_BinByBin_Pt->DrawCopy("same");
        //histoGammaCorrUnfoldReso_SvD_Pt->DrawCopy("same");
        //histoGammaCorrUnfoldReso_TUnfold_Pt->DrawCopy("same");

        TLegend* legendGammaSpecUnfoldComp                = GetAndSetLegend(0.6,0.8,2);
        legendGammaSpecUnfoldComp->AddEntry(histoGammaCorrUnfoldReso_Pt,"Bayesian unfolding");
        legendGammaSpecUnfoldComp->AddEntry(histoGammaCorrUnfoldReso_BinByBin_Pt,"Bin-by-Bin unfolding");
        legendGammaSpecUnfoldComp->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.68, 0.78, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
    canvasCorrGammaSpecUnfold->SaveAs(Form("%s/%s_%s_CorrGammaUnfoldSpectrumPurityMinusSec_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasCorrGammaSpecUnfold;

    //**********************************************************************************
    //******************** Resolution Correction Plot **********************************
    //**********************************************************************************
    TCanvas *canvasResolutionCorr               = GetAndSetCanvas("canvasResolutionCorr");
        TH1D* Dummy2                            = NULL;
        if (isPCM && !isCalo) Dummy2            = (TH1D*)histoESDConvGammaPt->Clone("Dummy2");
        if (isCalo && !isPCM) Dummy2            = (TH1D*)histoESDCaloGammaPt->Clone("Dummy2");
        SetHistogramm(Dummy2,"#it{p}_{T} (GeV/#it{c})", "resolution correction",0,2);
        Dummy2->Reset();
        Dummy2->Draw();

        DrawGammaSetMarker(histoGammaResolCorrUnfold_Pt, 20, 1.0, kBlue+1, kBlue+1);
        DrawGammaSetMarker(histoGammaResolCorrUnfold_BinByBin_Pt, 21, 1.0, kGray+2, kGray+2);
        DrawGammaSetMarker(histoGammaResolCorrEff_Pt, 24, 1.5, kBlack, kBlack);

        DrawGammaLines(0., maxPtGamma,1, 1,1, kGray);
        DrawGammaLines(0., maxPtGamma,0.8, 0.8,1, kGray, 7);
        DrawGammaLines(0., maxPtGamma,1.2, 1.2,1, kGray, 7);
        
        histoGammaResolCorrUnfold_BinByBin_Pt->DrawCopy("same");
        histoGammaResolCorrUnfold_Pt->DrawCopy("same");
        histoGammaResolCorrEff_Pt->DrawCopy("same");
        //histoGammaCorrUnfoldReso_SvD_Pt->DrawCopy("same");
        //histoGammaCorrUnfoldReso_TUnfold_Pt->DrawCopy("same");

        TLegend* legendResolutionCorr                = GetAndSetLegend(0.15,0.80,3);
        legendResolutionCorr->AddEntry(histoGammaResolCorrEff_Pt,"from effi ");
        legendResolutionCorr->AddEntry(histoGammaResolCorrUnfold_Pt,"from Bayesian unfolding");
        legendResolutionCorr->AddEntry(histoGammaResolCorrUnfold_BinByBin_Pt,"from bin-by-bin unfolding");
        legendResolutionCorr->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.18, 0.3, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
 
    canvasResolutionCorr->SaveAs(Form("%s/%s_%s_ResolutionCorrection_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasResolutionCorr;
    
    //******************************************************************************************
    //****************************************************************************************** 
    //****************************************************************************************** 
    // copy raw inc gamma spectrum and correct for secondary contamination, purity, conversion probability & reco effi (vs rec pt with sec correction)
    TH1D* histoGammaCorrEffiReso_Pt                 = NULL;
    if (isPCM && !isCalo) {
        histoGammaCorrEffiReso_Pt                   = (TH1D*)histoESDConvGammaPt->Clone("CorrGammaSpecPurityMinusSec");
        
        if (!hasCocktailInput)
            CorrectGammaEffiResol(histoGammaCorrEffiReso_Pt, histoSecondaryGammaSpecPt, histoSecondaryGammaFromXFromK0sSpecPt, histoGammaTruePurity_Pt, histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_Pt, deltaEta, scaling, nEvt);
        else
            CorrectGammaEffiResolCocktail(histoGammaCorrEffiReso_Pt, histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt, histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt, histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt, histoGammaTrueSecCocktailGammaRest_Pt, histoGammaTruePurity_Pt, histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_Pt, deltaEta, scaling, nEvt);
    }
    if (isCalo && !isPCM) {
        histoGammaCorrEffiReso_Pt                   = (TH1D*)histoESDCaloGammaPt->Clone("CorrGammaSpecPurityMinusSec");
        
        if (!hasCocktailInput)
            CorrectGammaEffiResol( histoGammaCorrEffiReso_Pt, histoSecondaryGammaSpecPt, histoSecondaryGammaFromXFromK0sSpecPt,
                                  histoGammaTruePurity_Pt, histoGammaPrimaryRecoEff_Pt, deltaEtaCalo, scalingCalo, nEvt);
        else
            CorrectGammaEffiResolCocktail(histoGammaCorrEffiReso_Pt, histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt, histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt, histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt, histoGammaTrueSecCocktailGammaRest_Pt, histoGammaTruePurity_Pt, histoGammaPrimaryRecoEff_Pt, deltaEtaCalo, scalingCalo, nEvt);
    }
    DrawGammaSetMarker(histoGammaCorrEffiReso_Pt, 20, 1.0, kGreen+2, kGreen+2);
    SetHistogramm(histoGammaCorrEffiReso_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");

    // Do same as before with secondary corrections applied
    TH1D* histoGammaCorrEffiReso_PileUp_Pt              = NULL;
    TH1D* histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt    = NULL;
    if(doPileUpCorr){
        // copy raw inc gamma spectrum after pure data driven pileup correction, correct with pileup corrected secondary histograms, pileup corrected purity, conversion prob, reco effi
        // (pileup corrected, vs rec pT with sec correction) all contain MC update of pileup correction
        histoGammaCorrEffiReso_PileUp_Pt                = (TH1D*)histoESDConvGammaPtPileUp->Clone("CorrGammaSpecPurityMinusSecPileUp");
        
        if (!hasCocktailInput) {
            CorrectGammaEffiResol( histoGammaCorrEffiReso_PileUp_Pt, histoSecondaryGammaSpecPtPileUp, histoSecondaryGammaFromXFromK0sSpecPtPileUp,
                                  histoGammaTruePurity_PileUp_Pt, histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_PileUp_Pt, deltaEta, scaling, nEvt);
        } else {
            cout << __LINE__ << ": MC updated bin-by-bin correction for pileup corrected spectra not implemented yet for use of cocktail-based secondary correction, will use MC-based secondary correction." << endl;
            CorrectGammaEffiResol( histoGammaCorrEffiReso_PileUp_Pt, histoSecondaryGammaSpecPtPileUp, histoSecondaryGammaFromXFromK0sSpecPtPileUp,
                                  histoGammaTruePurity_PileUp_Pt, histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_PileUp_Pt, deltaEta, scaling, nEvt);
        }
        
        DrawGammaSetMarker(histoGammaCorrEffiReso_PileUp_Pt, 20, 1.0, kMagenta, kMagenta);
        SetHistogramm(histoGammaCorrEffiReso_PileUp_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");

        // copy raw inc gamma spectrum after pure data driven pileup correction, correct with pileup corrected secondary histograms, pileup corrected purity, conversion prob, reco effi
        // (pileup corrected, vs rec pT with sec correction)
        histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt    = (TH1D*)histoESDConvGammaPtPileUp->Clone("CorrGammaSpecPurityMinusSecPileUpNoMCUpdate");
        
        if (!hasCocktailInput) {
            CorrectGammaEffiResol( histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt, histoSecondaryGammaSpecPt, histoSecondaryGammaFromXFromK0sSpecPt,
                                  histoGammaTruePurity_Pt, histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_Pt, deltaEta, scaling, nEvt);
        } else {
            CorrectGammaEffiResolCocktail( histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt, histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt, histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt, histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt, histoGammaTrueSecCocktailGammaRest_Pt, histoGammaTruePurity_Pt, histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_Pt, deltaEta, scaling, nEvt);
        }
        
        DrawGammaSetMarker(histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt, 20, 1.0, kBlue+2, kBlue+2);
        SetHistogramm(histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt->Draw("same");
    }

    //*************************************************************************************************
    //********** Sanity check for photon spectrum (do we get back what we put in) *********************
    //*************************************************************************************************
    // correct input inc gamma spectrum for delta Eta, delta phi & number of events
    TH1D* histoMCGammaSpec_MCPt                      = (TH1D*)histoMCAllGamma_MCPt->Clone("GammaSpecMC");
    if (isPCM && !isCalo) CorrectGammaMC(histoMCGammaSpec_MCPt, deltaEta,       scaling,        nEvtMC);
    if (isCalo && !isPCM) CorrectGammaMC(histoMCGammaSpec_MCPt, deltaEtaCalo,   scalingCalo,    nEvtMC);
    DrawGammaSetMarker(histoMCGammaSpec_MCPt, 20, 1.0, 2, 2);
    SetHistogramm(histoMCGammaSpec_MCPt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");

    //*************************************************************************************************
    //********************************* Comparison Gamma Spec *****************************************
    //*************************************************************************************************
    TH1D* histoRatioGammaCorrMinusSecGammaUnfold            = (TH1D*)histoGammaCorrEffiReso_Pt->Clone("histoRatioGammaCorrMinusSecGammaUnfold");
    histoRatioGammaCorrMinusSecGammaUnfold->Divide(histoGammaCorrUnfoldReso_Pt);
    TH1D* histoRatioGammaCorrMinusSecGammaPileUp            = NULL;
    TH1D* histoRatioGammaCorrMinusSecGammaPileUpNoMCUpdate  = NULL;
    if(doPileUpCorr){
        histoRatioGammaCorrMinusSecGammaPileUp              = (TH1D*)histoGammaCorrEffiReso_Pt->Clone("histoRatioGammaCorrMinusSecGammaPileUp");
        histoRatioGammaCorrMinusSecGammaPileUp->Divide(histoGammaCorrEffiReso_PileUp_Pt);

        histoRatioGammaCorrMinusSecGammaPileUpNoMCUpdate    = (TH1D*)histoGammaCorrEffiReso_Pt->Clone("histoRatioGammaCorrMinusSecGammaPileUp");
        histoRatioGammaCorrMinusSecGammaPileUpNoMCUpdate->Divide(histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt);
    }
    
    //*************************************************************************************************
    //********************************* Decay gamma to all gamma **************************************
    //*************************************************************************************************
    TCanvas *canvasRatioAllDiffDecay    = GetAndSetCanvas("canvasRatioAllDiffDecay");
    canvasRatioAllDiffDecay->SetLogy(0);
    canvasRatioAllDiffDecay->cd();
    
        TH1D* histoAllDiffDecayGamma    = (TH1D*)histoMCAllGamma_OriginalBin_MCPt->Clone("RatioAllToDecay");
        histoAllDiffDecayGamma->Scale(histoAllDiffDecayGamma->GetBinWidth(5));
        histoAllDiffDecayGamma->Rebin(5);
        histoAllDiffDecayGamma->Scale(1./histoAllDiffDecayGamma->GetBinWidth(5));
        
        TH1D* histoAllDecayGamma        = (TH1D*)histoPhotonSource_MCPt[7]->Clone("AllDecay");
        histoAllDecayGamma->Scale(histoAllDecayGamma->GetBinWidth(5));
        histoAllDecayGamma->Rebin(5);
        histoAllDecayGamma->Scale(1./histoAllDecayGamma->GetBinWidth(5));
        
        histoAllDiffDecayGamma->Divide(histoAllDiffDecayGamma,histoAllDecayGamma,1,1,"B");
        histoAllDiffDecayGamma->GetYaxis()->SetRangeUser(0.8,1.5);
        
        SetHistogramm(histoAllDiffDecayGamma,"#it{p}_{T} (GeV/#it{c})", "#frac{N_{#gamma_{incl}}}{N_{#gamma_{decay}}}");
        histoAllDiffDecayGamma->Draw("p,e1");

        DrawGammaLines(0., 50,1, 1,0.5, kGray+2, 7);
 
        PutProcessLabelAndEnergyOnPlot( 0.7, 0.95, 0.035, cent, "#gamma", "", 42, 0.03);
        
    canvasRatioAllDiffDecay->SaveAs(Form("%s/%s_AllGammaDivDecayGammaSpectrumMC_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasRatioAllDiffDecay;
   
    //*************************************************************************************************
    //********************************* Gamma from Decay **********************************************
    //*************************************************************************************************
    // correct input spectra of different decay channels with delta eta, delta phi & events number
    for (Int_t i = 0; i< 9; i++){
        CorrectGammaMC(histoPhotonSource_MCPt[i], deltaEta, scaling, nEvtMC);
        SetHistogramm(histoPhotonSource_MCPt[i],"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        
        DrawGammaSetMarker(histoPhotonSource_MCPt[i], 22, 1.0, colorsDecay[i] , colorsDecay[i]);
    }
        
    TCanvas *canvasDecayGammaSpecMC = GetAndSetCanvas("canvasDecayGammaSpecMC");
    canvasDecayGammaSpecMC->SetLogy();
    
        TLegend* legendDecaySpectra = GetAndSetLegend(0.5,0.7,5,2);
        histoPhotonSource_MCPt[7]->DrawCopy("chist");
        for (Int_t i = 0; i< 9; i++){   
            histoPhotonSource_MCPt[i]->DrawCopy("chist,same");
            legendDecaySpectra->AddEntry(histoPhotonSource_MCPt[i],decaysLabels[i].Data());
        }
        legendDecaySpectra->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.7, 0.68, 0.035, cent, "#gamma", detectionProcess, 42, 0.03);
        
    canvasDecayGammaSpecMC->SaveAs(Form("%s/%s_DecayGammaSpectrumMC_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasDecayGammaSpecMC;
    
    //*************************************************************************************************
    //**************** Compare cocktail secondaries to MC secondaries *********************************
    //*************************************************************************************************
    if(hasCocktailInput){
        TH1D* histoRatioSecondariesCocktailMCK0s_Raw_Pt     = NULL;
        TH1D* histoRatioSecondariesCocktailMCK0l_Raw_Pt     = NULL;
        TH1D* histoRatioSecondariesCocktailMCLambda_Raw_Pt  = NULL;
        TCanvas *canvasSecondaryComparison                  = GetAndSetCanvas("canvasSecondaryComparison");
        
        histoRatioSecondariesCocktailMCK0s_Raw_Pt           = (TH1D*)histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt->Clone("histoRatioSecondariesCocktailMCK0s_Raw_Pt");
        if (isPCM && !isCalo) histoRatioSecondariesCocktailMCK0s_Raw_Pt->Divide(histoGammaTrueSecConvGammaFromXFromK0s_Pt);
        if (isCalo && !isPCM) histoRatioSecondariesCocktailMCK0s_Raw_Pt->Divide(histoGammaTrueSecCaloGammaFromXFromK0s_Pt);
        SetHistogramm(histoRatioSecondariesCocktailMCK0s_Raw_Pt,"#it{p}_{T} (GeV/#it{c})", "Cocktail/MC");
        DrawGammaSetMarker(histoRatioSecondariesCocktailMCK0s_Raw_Pt, 20, 2.0, 1 , 1);

        if(0){
            histoRatioSecondariesCocktailMCK0l_Raw_Pt       = (TH1D*) histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt -> Clone("histoRatioSecondariesCocktailMCK0l_Raw_Pt");
            histoRatioSecondariesCocktailMCK0l_Raw_Pt->Divide(histoGammaTrueSecConvGammaFromXFromK0l_Pt);
            SetHistogramm(histoRatioSecondariesCocktailMCK0l_Raw_Pt,"#it{p}_{T} (GeV/#it{c})", "Cocktail/MC");
            DrawGammaSetMarker(histoRatioSecondariesCocktailMCK0l_Raw_Pt, 22, 2.0, kBlue-4 , kBlue-4);
            
            histoRatioSecondariesCocktailMCLambda_Raw_Pt    = (TH1D*) histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt -> Clone("histoRatioSecondariesCocktailMCLambda_Raw_Pt");
            histoRatioSecondariesCocktailMCLambda_Raw_Pt->Divide(histoGammaTrueSecConvGammaFromXFromLambda_Pt);
            SetHistogramm(histoRatioSecondariesCocktailMCLambda_Raw_Pt,"#it{p}_{T} (GeV/#it{c})", "Cocktail/MC");
            DrawGammaSetMarker(histoRatioSecondariesCocktailMCLambda_Raw_Pt, 24, 2.0, kCyan+2 , kCyan+2);
        }
        histoRatioSecondariesCocktailMCK0s_Raw_Pt->Draw();
        DrawGammaLines(0,histoRatioSecondariesCocktailMCK0s_Raw_Pt->GetXaxis()->GetXmax(),1,1,1,kGray+2,2);
        histoRatioSecondariesCocktailMCK0s_Raw_Pt->Draw("same");
        if(0){
            histoRatioSecondariesCocktailMCK0l_Raw_Pt->Draw("same");
            histoRatioSecondariesCocktailMCLambda_Raw_Pt->Draw("same");
        }
        TLegend* legendCompareSecCocktailMC                 = GetAndSetLegend(0.2,0.75,5,2);
        legendCompareSecCocktailMC->SetBorderSize(0);
        legendCompareSecCocktailMC->AddEntry(histoRatioSecondariesCocktailMCK0s_Raw_Pt,"Secondary #gamma from K^{0}_{S}","p");
        legendCompareSecCocktailMC->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.7, 0.28, 0.035, cent, "#gamma", detectionProcess, 42, 0.03);
        
        canvasSecondaryComparison->SaveAs(Form("%s/%s_SecondaryComparisonCocktailMC_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasSecondaryComparison;
    }
    
    //*************************************************************************************************
    //******************************* Gamma Combinatorial Background **********************************
    //*************************************************************************************************
    // pure BG plotting
    TCanvas *canvasCombBackSpecMC = GetAndSetCanvas("canvasCombBackSpecMC");
    canvasCombBackSpecMC->SetLogy();

        TLegend* legendCombSpectra = GetAndSetLegend(0.25,0.75,4,3);
        if (isPCM && !isCalo) {
            for(Int_t i = 0;i<17;i++){
                SetHistogramm(histoCombinatorialSpecies_Pt[i],"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
                DrawGammaSetMarker(histoCombinatorialSpecies_Pt[i], markersCombinatorics[i], 1., colorsCombinatorics[i], colorsCombinatorics[i]);
                histoCombinatorialSpecies_Pt[i]->Scale(1./nEvtMC);
            
                if(i==0) histoCombinatorialSpecies_Pt[i]->DrawCopy("");
                else histoCombinatorialSpecies_Pt[i]->DrawCopy("same");
            
                legendCombSpectra->AddEntry(histoCombinatorialSpecies_Pt[i],combinatorics[i]);
            }
        }
        if (isCalo && !isPCM) {
            for(Int_t i = 0;i<11;i++){
                SetHistogramm(histoCombinatorialSpeciesCalo_Pt[i],"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
                DrawGammaSetMarker(histoCombinatorialSpeciesCalo_Pt[i], markersCombinatoricsCalo[i], 1., colorsCombinatoricsCalo[i], colorsCombinatoricsCalo[i]);
                histoCombinatorialSpeciesCalo_Pt[i]->Scale(1./nEvtMC);
                
                if(i==0) {
                    histoCombinatorialSpeciesCalo_Pt[i]->GetYaxis()->SetRangeUser(1e-10, 1e-3);
                    histoCombinatorialSpeciesCalo_Pt[i]->DrawCopy("");
                }
                else histoCombinatorialSpeciesCalo_Pt[i]->DrawCopy("same");
                
                legendCombSpectra->AddEntry(histoCombinatorialSpeciesCalo_Pt[i],combinatoricsCalo[i]);
            }
        }

        legendCombSpectra->Draw();
        PutProcessLabelAndEnergyOnPlot( 0.75, 0.68, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
    
    canvasCombBackSpecMC->SaveAs(Form("%s/%s_CombBackSpectrumMC_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasCombBackSpecMC;
    
    // plot ratio to real primary photons
    TCanvas *canvasSignalToCombBackgroundRatio      = GetAndSetCanvas("canvasSignalToCombBackgroundRatio");
    canvasSignalToCombBackgroundRatio->SetLogy();

        TLegend* legendSignalToCombBackgroundRatio  = GetAndSetLegend(0.15,0.76,4,3);
        TH1D **histoSignalToCombBackgroundRatio     = NULL;
        TH1D *SummedSmallContributionsCombBack      = NULL; //10,11,12,13,14,15
        if (isPCM && !isCalo) {
            histoSignalToCombBackgroundRatio        = new TH1D*[17];
            
            for(Int_t i = 0;i<16;i++){
                histoSignalToCombBackgroundRatio[i] = (TH1D*) histoCombinatorialSpecies_Pt[i]->Clone(Form("ESD_TrueCombRatioSignal_%s_Pt",combinatorics[i].Data()));
                histoSignalToCombBackgroundRatio[i]->Scale(nEvtMC);
                
                if(i==9){
                    SummedSmallContributionsCombBack = (TH1D*)histoSignalToCombBackgroundRatio[i]->Clone("SummedSmallContributions");
                    SetHistogramm(SummedSmallContributionsCombBack,"#it{p}_{T} (GeV/#it{c})","SummedSmallContributions",10,5e7);
                    SummedSmallContributionsCombBack->SetMinimum(1e-5);
                }
                
                if(i>9){
                    SummedSmallContributionsCombBack->Add(histoSignalToCombBackgroundRatio[i]);
                }
                histoSignalToCombBackgroundRatio[i]->Divide(histoSignalToCombBackgroundRatio[i],histoGammaTruePrimaryConv_Pt,1,1,"");
                SetHistogramm(histoSignalToCombBackgroundRatio[i],"#it{p}_{T} (GeV/#it{c})","identified BG/ real primary photons",10,5e7);
                histoSignalToCombBackgroundRatio[i]->SetMinimum(1e-5);
                
                if(i==0){
                    histoSignalToCombBackgroundRatio[i]->GetYaxis()->SetRangeUser(1e-5,100);
                    histoSignalToCombBackgroundRatio[i]->DrawCopy("e1");
                }
                if(i<9){
                    DrawGammaSetMarker(histoSignalToCombBackgroundRatio[i], markersCombinatorics[i], 1., colorsCombinatorics[i], colorsCombinatorics[i]);
                    histoSignalToCombBackgroundRatio[i]->DrawCopy("e1same");
                }
                else continue;
                
                legendSignalToCombBackgroundRatio->AddEntry(histoSignalToCombBackgroundRatio[i],combinatorics[i]);
            }
            
            SummedSmallContributionsCombBack->Divide(SummedSmallContributionsCombBack,histoGammaTruePrimaryConv_Pt,1,1,"");
            legendSignalToCombBackgroundRatio->AddEntry(SummedSmallContributionsCombBack,"Protos+Kaons+Muons");
        }
        if (isCalo && !isPCM) {
            histoSignalToCombBackgroundRatio        = new TH1D*[10];
        
            for(Int_t i = 0;i<10;i++){
                histoSignalToCombBackgroundRatio[i] = (TH1D*)histoCombinatorialSpeciesCalo_Pt[i]->Clone(Form("ESD_TrueCombRatioSignal_%s_Pt",combinatoricsCalo[i].Data()));
                histoSignalToCombBackgroundRatio[i]->Scale(nEvtMC);
            
                if(i==2){
                    SummedSmallContributionsCombBack = (TH1D*)histoSignalToCombBackgroundRatio[i]->Clone("SummedSmallContributions");
                    SetHistogramm(SummedSmallContributionsCombBack,"#it{p}_{T} (GeV/#it{c})","SummedSmallContributions",10,5e7);
                    SummedSmallContributionsCombBack->SetMinimum(1e-5);
                }
                if(i==5 || i==6 || i==7 || i==9){
                    SummedSmallContributionsCombBack->Add(histoSignalToCombBackgroundRatio[i]);
                }
                
                histoSignalToCombBackgroundRatio[i]->Divide(histoSignalToCombBackgroundRatio[i],histoGammaTruePrimaryCalo_Pt,1,1,"");
                SetHistogramm(histoSignalToCombBackgroundRatio[i],"#it{p}_{T} (GeV/#it{c})","identified BG/real primary photons",10,5e7);
                histoSignalToCombBackgroundRatio[i]->SetMinimum(1e-5);
            
                if(i==0){
                    histoSignalToCombBackgroundRatio[i]->GetYaxis()->SetRangeUser(1e-5,100);
                    histoSignalToCombBackgroundRatio[i]->DrawCopy("e1");
                }
                if(i<2 || (i>2 && i<5) || i==8){
                    DrawGammaSetMarker(histoSignalToCombBackgroundRatio[i], markersCombinatoricsCalo[i], 1., colorsCombinatoricsCalo[i], colorsCombinatoricsCalo[i]);
                    histoSignalToCombBackgroundRatio[i]->DrawCopy("e1same");
                } else continue;
            
                legendSignalToCombBackgroundRatio->AddEntry(histoSignalToCombBackgroundRatio[i],combinatoricsCalo[i]);
            }
            
            SummedSmallContributionsCombBack->Divide(SummedSmallContributionsCombBack,histoGammaTruePrimaryCalo_Pt,1,1,"");
            legendSignalToCombBackgroundRatio->AddEntry(SummedSmallContributionsCombBack,"p+K0s+Lambda+mu+rest");
        }
   
        DrawGammaSetMarker(SummedSmallContributionsCombBack, markersCombinatorics[10], 1., colorsCombinatorics[10], colorsCombinatorics[10]);
        SummedSmallContributionsCombBack->DrawCopy("e1same");
        legendSignalToCombBackgroundRatio->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.2, 0.76, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
    canvasSignalToCombBackgroundRatio->SaveAs(Form("%s/%s_CombBackgroundRatioToSignal_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasSignalToCombBackgroundRatio;

    // plot ratio to summed total MC BG
    TCanvas *canvasRatioCombBackToBack              = GetAndSetCanvas("canvasRatioCombBackToBack");

    canvasRatioCombBackToBack->SetLogy();
        TLegend* legendRatioCombBackToBack          = GetAndSetLegend(0.15,0.76,4,3);
        
        TH1D **histoRatioCombBackToBack                 = NULL;
        TH1D *SummedSmallContributionsCombBackToBack    = NULL; //10,11,12,13,14,15
        if (isPCM && !isCalo) {
            histoRatioCombBackToBack                    = new TH1D*[17];
   
            for(Int_t i = 0;i<16;i++){
            
                histoRatioCombBackToBack[i] = (TH1D*) histoCombinatorialSpecies_Pt[i]->Clone(Form("ESD_TrueCombRatioSignal_%s_Pt",combinatorics[i].Data()));
            
                if(i==9){
                    SummedSmallContributionsCombBackToBack = (TH1D*)histoRatioCombBackToBack[i]->Clone("SummedSmallContributions");
                    SetHistogramm(SummedSmallContributionsCombBackToBack,"#it{p}_{T} (GeV/#it{c})","SummedSmallContributions",10,5e7);
                    SummedSmallContributionsCombBackToBack->SetMinimum(1e-5);
                }
            
                if(i>9){
                    SummedSmallContributionsCombBackToBack->Add(histoRatioCombBackToBack[i]);
                }
                histoRatioCombBackToBack[i]->Divide(histoRatioCombBackToBack[i],histoGammaMCBackground_Pt,1,1,"B");
                SetHistogramm(histoRatioCombBackToBack[i],"#it{p}_{T} (GeV/#it{c})","identified BG/Total BG",10,5e7);
                histoRatioCombBackToBack[i]->SetMinimum(1e-5);
            
                if(i==0){
                    histoRatioCombBackToBack[i]->GetYaxis()->SetRangeUser(1e-4,400);
                    histoRatioCombBackToBack[i]->DrawCopy("e1");
                }
                if(i<9){
                    DrawGammaSetMarker(histoRatioCombBackToBack[i], markersCombinatorics[i], 1., colorsCombinatorics[i], colorsCombinatorics[i]);
                    histoRatioCombBackToBack[i]->DrawCopy("e1same");
                }
                else continue;
            
                legendRatioCombBackToBack->AddEntry(histoRatioCombBackToBack[i],combinatorics[i]);
            }
            
            legendRatioCombBackToBack->AddEntry(SummedSmallContributionsCombBackToBack,"Protos+Kaons+Muons");
        }
        if (isCalo && !isPCM) {
            histoRatioCombBackToBack                    = new TH1D*[10];
            
            for(Int_t i = 0;i<10;i++){
                histoRatioCombBackToBack[i] = (TH1D*)histoCombinatorialSpeciesCalo_Pt[i]->Clone(Form("ESD_TrueCombRatioSignal_%s_Pt",combinatoricsCalo[i].Data()));
                
                if(i==2){
                    SummedSmallContributionsCombBackToBack = (TH1D*)histoRatioCombBackToBack[i]->Clone("SummedSmallContributions");
                    SetHistogramm(SummedSmallContributionsCombBackToBack,"#it{p}_{T} (GeV/#it{c})","SummedSmallContributions",10,5e7);
                    SummedSmallContributionsCombBackToBack->SetMinimum(1e-5);
                }
                if(i==5 || i==6 || i==7 || i==9){
                    SummedSmallContributionsCombBackToBack->Add(histoRatioCombBackToBack[i]);
                }
                
                histoRatioCombBackToBack[i]->Divide(histoRatioCombBackToBack[i],histoGammaMCBackground_Pt,1,1,"B");
                SetHistogramm(histoRatioCombBackToBack[i],"#it{p}_{T} (GeV/#it{c})","identified BG/Total BG",10,5e7);
                histoRatioCombBackToBack[i]->SetMinimum(1e-5);
                
                if(i==0){
                    histoRatioCombBackToBack[i]->GetYaxis()->SetRangeUser(1e-4,400);
                    histoRatioCombBackToBack[i]->DrawCopy("e1");
                }
                if(i<2 || (i>2 && i<5) || i==8){
                    DrawGammaSetMarker(histoRatioCombBackToBack[i], markersCombinatoricsCalo[i], 1., colorsCombinatoricsCalo[i], colorsCombinatoricsCalo[i]);
                    histoRatioCombBackToBack[i]->DrawCopy("e1same");
                } else continue;
                
                legendRatioCombBackToBack->AddEntry(histoRatioCombBackToBack[i],combinatoricsCalo[i]);
            }
            
            legendRatioCombBackToBack->AddEntry(SummedSmallContributionsCombBackToBack,"p+K0s+Lambda+mu+rest");
        }
    
        SummedSmallContributionsCombBackToBack->Divide(SummedSmallContributionsCombBackToBack,histoGammaMCBackground_Pt,1,1,"B");
        DrawGammaSetMarker(SummedSmallContributionsCombBackToBack, markersCombinatorics[10], 1., colorsCombinatorics[10], colorsCombinatorics[10]);
        SummedSmallContributionsCombBackToBack->DrawCopy("e1same");
        legendRatioCombBackToBack->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.75, 0.3, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
    
    canvasRatioCombBackToBack->SaveAs(Form("%s/%s_RatioCombBackToBack_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasRatioCombBackToBack;
    
    //*************************************************************************************************
    //*********************** Compare Gamma Yields MC/Reconstructed ***********************************
    //*************************************************************************************************
    TCanvas *canvasPurBackGammaSpec = GetAndSetCanvas("canvasPurBackGammaSpec");
    DrawGammaCanvasSettings( canvasPurBackGammaSpec, 0.1, 0.015, 0.02, 0.09);
        TPad *padPurBackSpec = new TPad("padNLOHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padPurBackSpec, 0.10, 0.015, 0.02, 0.);
        padPurBackSpec->Draw();
        TPad *padPurBackRatio = new TPad("padNLORatios", "", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings( padPurBackRatio, 0.1, 0.015, 0.0, 0.35);
        padPurBackRatio->Draw();
    
        padPurBackSpec->cd();
        padPurBackSpec->SetLogy();

        DrawGammaSetMarker(histoGammaCorrEffiReso_Pt, 20, 1.0, kBlue+2,  kBlue+2);
        DrawGammaSetMarker(histoGammaCorrUnfoldReso_Pt, 24, 1.0, kAzure+7,  kAzure+7);

        histoGammaCorrEffiReso_Pt->DrawCopy("e1");
        histoGammaCorrUnfoldReso_Pt->DrawCopy("e1,same");
        
        TLegend* legendPurBackSpectra = GetAndSetLegend(0.5,0.85,2);
        legendPurBackSpectra->AddEntry(histoGammaCorrEffiReso_Pt,"#gamma with Purity and corrected for Secondaries","p");
        legendPurBackSpectra->AddEntry(histoGammaCorrUnfoldReso_Pt,"#gamma unfolded with Secondary correction","p");
        legendPurBackSpectra->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.18, 0.2, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
        padPurBackRatio->cd();
        TH1D *histGammaRatioUnfoldMinusSec = (TH1D*) histoGammaCorrEffiReso_Pt->Clone("histGammaRatioMinusSec");
        histGammaRatioUnfoldMinusSec->Divide(histGammaRatioUnfoldMinusSec,histoGammaCorrUnfoldReso_Pt,1,1,"");
        histGammaRatioUnfoldMinusSec->GetYaxis()->SetRangeUser(0.88,1.12);
        SetStyleHistoTH1ForGraphs(  histGammaRatioUnfoldMinusSec, "#it{p}_{T} (GeV/#it{c})","#frac{#gamma corr legacy}{#gamma corr unfold}" , 0.14, 0.15,  
                                    0.12, 0.10,  0.85, 0.4, 510, 505);
        
        DrawGammaSetMarker(histGammaRatioUnfoldMinusSec, 24, 1.0, kAzure+7,  kAzure+7);
        histGammaRatioUnfoldMinusSec->DrawCopy("");
        
        DrawGammaLines(0., maxPtGamma,1, 1,0.5);
    
    canvasPurBackGammaSpec->SaveAs(Form("%s/%s_%s_GammaSpectraComparison_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasPurBackGammaSpec;
    
    if (isPCM && !isCalo) {
        CorrectGammaMC(histoMCGammaConv_MCPt, deltaEta, scaling, nEvtMC);
        CorrectGammaMC(histoGammaTrueConv_Pt, deltaEta, scaling, nEvtMC);
    }
    if (isCalo && !isPCM) {
        CorrectGammaMC(histoGammaTrueCalo_Pt, deltaEtaCalo, scalingCalo, nEvtMC);
    }

    TCanvas *canvasGammaPurityConvBin = GetAndSetCanvas("canvasGammaPurityConvBin");
    DrawGammaCanvasSettings(canvasGammaPurityConvBin,0.1, 0.015, 0.02, 0.09);
    TPad* padGammaPurityConvBin       = new TPad("padGammaPurityConvBin", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings(padGammaPurityConvBin, 0.10, 0.015, 0.02, 0.);
    padGammaPurityConvBin->Draw();
    TPad* padGammaPurityConvBinRatio  = new TPad("padGammaPurityConvBinRatio", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings(padGammaPurityConvBinRatio,  0.1, 0.015, 0.0, 0.3);
    padGammaPurityConvBinRatio->Draw();

    padGammaPurityConvBin->cd();
    padGammaPurityConvBin->SetLogy();
    
        DrawGammaSetMarker(histoGammaCorrEffiReso_Pt, 24, 1.0, 1, 1);
        histoMCGammaSpec_MCPt->GetYaxis()->SetRangeUser(5e-10,100);
        histoMCGammaSpec_MCPt->DrawCopy("e1");
    
        DrawGammaSetMarker(histoGammaCorrUnfoldReso_Pt, 24, 1.0, 1, 1);
        histoGammaCorrUnfoldReso_Pt->DrawCopy("e1,same");
    
        TLegend *legendGammaSpectraConvBin = GetAndSetLegend(0.45,0.8,2);
        legendGammaSpectraConvBin->AddEntry(histoGammaCorrUnfoldReso_Pt,"corrected data #gamma spectrum");
        legendGammaSpectraConvBin->AddEntry(histoMCGammaSpec_MCPt,"input MC #gamma spectrum");
        legendGammaSpectraConvBin->Draw();
    
        PutProcessLabelAndEnergyOnPlot( 0.18, 0.3, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
    
    padGammaPurityConvBinRatio->cd();

        TH1D* tempRatioDataMC               = NULL;
        tempRatioDataMC                     = (TH1D*)histoGammaCorrUnfoldReso_Pt->Clone("tempRatioDataMC");
        tempRatioDataMC->Divide(tempRatioDataMC,histoMCGammaSpec_MCPt,1.,1.,"");
    
        SetStyleHistoTH1ForGraphs(tempRatioDataMC, "#it{p}_{T} (GeV/#it{c})", "#frac{data}{MC}", 0.14, 0.15, 0.12, 0.10, 0.85, 0.4, 510, 505);
        DrawGammaSetMarker(tempRatioDataMC, 24, 1.0, kBlack, kBlack);
        tempRatioDataMC->GetYaxis()->SetRangeUser(0.5, 1.5);
    
        tempRatioDataMC->DrawCopy("e1");
        DrawGammaLines(0., tempRatioDataMC->GetXaxis()->GetBinUpEdge(tempRatioDataMC->GetNbinsX()), 1., 1., 0.5, kBlue+1);
    
    canvasGammaPurityConvBin->SaveAs(Form("%s/%s_%s_GammaSpec_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasGammaPurityConvBin;

    TH1D* histoGammaCorrFac_Pt           = (TH1D*) histoGammaPrimaryRecoEff_Pt->Clone("GammaCorrFac_Pt");
    if (isPCM && !isCalo) histoGammaCorrFac_Pt->Multiply(histoGammaConvProb_MCPt);
    
    //***********************************************************************************************************
    //******************************* Write output file for photons *********************************************
    //***********************************************************************************************************
    TString nameOutput = Form("%s/%s/%s_%s_GammaConvV1Correction_%s.root",cutSelection.Data(),option.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data());
    TFile* fileCorrectedOutput = new TFile(nameOutput,"RECREATE");

        //________________________ writing MC quantities to file _____________________________
        // input spectrum corrected for deta, dPhi, Nevt
        if (histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt)     histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt->Write("histoGammaSecGammaFromXFromK0s_Cocktail_Raw_Pt", TObject::kOverwrite);
        if (histoGammaTrueSecCocktailGammaFromXFromK0s_Pt)      histoGammaTrueSecCocktailGammaFromXFromK0s_Pt->Write("histoGammaTrueSecCocktailGammaFromXFromK0s_Pt", TObject::kOverwrite);
        if (histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt)     histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt->Write("histoGammaSecGammaFromXFromK0l_Cocktail_Raw_Pt", TObject::kOverwrite);
        if (histoGammaTrueSecCocktailGammaFromXFromK0l_Pt)      histoGammaTrueSecCocktailGammaFromXFromK0l_Pt->Write("histoGammaTrueSecCocktailGammaFromXFromK0l_Pt", TObject::kOverwrite);
        if (histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt)  histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt->Write("histoGammaSecGammaFromXFromLambda_Cocktail_Raw_Pt", TObject::kOverwrite);
        if (histoGammaTrueSecCocktailGammaFromXFromLambda_Pt)   histoGammaTrueSecCocktailGammaFromXFromLambda_Pt->Write("histoGammaTrueSecCocktailGammaFromXFromLambda_Pt", TObject::kOverwrite);
        if (histoMCGammaSpec_MCPt)                              histoMCGammaSpec_MCPt->Write("GammaSpecMC", TObject::kOverwrite);
        if (histoSecondaryGammaSpecPt)                          histoSecondaryGammaSpecPt->Write("histoSecondaryGammaSpecPt", TObject::kOverwrite);
        if (histoSecondaryGammaFromXFromK0sSpecPt)              histoSecondaryGammaFromXFromK0sSpecPt->Write("histoSecondaryGammaFromXFromK0sSpecPt", TObject::kOverwrite);
        if (histoGammaTrueSecCocktailGammaRest_Pt)              histoGammaTrueSecCocktailGammaRest_Pt->Write("histoGammaTrueSecCocktailGammaRest_Pt", TObject::kOverwrite);
    
        // reconstructed MC gamma spectrum corrected deta, dPhi, Nevt
        if (histoMCrecGammaCorr_Pt)                             histoMCrecGammaCorr_Pt->Write("GammaSpecCorrESDMC", TObject::kOverwrite);
    
        // split in different source (pi0,eta,...)
        for (Int_t i= 0; i< 8; i++)
            if (histoPhotonSource_MCPt[i])                      histoPhotonSource_MCPt[i]->Write(histoPhotonSource_MCPt[i]->GetName(), TObject::kOverwrite);
    
        // pure direct photons
        if (histoPhotonSource_MCPt[8])                          histoPhotonSource_MCPt[8]->Write("MC_DirectPhoton_Pt", TObject::kOverwrite);
        // input spectrum of converted photons corrected for deta, dphi, Nevt
        if (histoMCGammaConv_MCPt)                              histoMCGammaConv_MCPt->Write("MC_ConvGamma_MCPt", TObject::kOverwrite);
        // reconstructed MC true photons corrected for deta, dphi, Nevt
        if (histoGammaTrueConv_Pt)                              histoGammaTrueConv_Pt->Write("TrueConvGamma_Pt", TObject::kOverwrite);
        if (histoGammaTrueCalo_Pt)                              histoGammaTrueCalo_Pt->Write("TrueCaloGamma_Pt", TObject::kOverwrite);
        // recontructed MC photon candidates, uncorrected
        if (histoMCrecGamma_Pt)                                 histoMCrecGamma_Pt->Write("MCrec_ConvGamma_Pt", TObject::kOverwrite);
        if (histoMCrecGammaCalo_Pt)                             histoMCrecGammaCalo_Pt->Write("MCrec_CaloGamma_Pt", TObject::kOverwrite);

        // summed reconstructed BG in MC
        if (histoMCrecBackground_Pt)                            histoMCrecBackground_Pt->Write("MCrec_Background", TObject::kOverwrite);
        // combinatorial histos split
        if (histoCombinatorialSpecies_Pt)
            for (Int_t i = 0;i<17;i++)                          histoCombinatorialSpecies_Pt[i]->Write(histoCombinatorialSpecies_Pt[i]->GetName(), TObject::kOverwrite);

        if (histoCombinatorialSpeciesCalo_Pt)
            for (Int_t i=0; i<11; i++)                          histoCombinatorialSpeciesCalo_Pt[i]->Write(histoCombinatorialSpeciesCalo_Pt[i]->GetName(), TObject::kOverwrite);


        //_________________________ writing correction factors to file _________________________
        // photon purity without secondary subtraction
        if (histoGammaPurity_Pt)                                histoGammaPurity_Pt->Write("GammaPurityWSec_Pt", TObject::kOverwrite);
        if (histoGammaCaloPurity_Pt)                            histoGammaCaloPurity_Pt->Write("histoGammaCaloPurityWSec_Pt", TObject::kOverwrite);
        // photon purity with secondary subtraction
        if (histoGammaTruePurity_Pt)                            histoGammaTruePurity_Pt->Write("GammaPurityWOSec_Pt", TObject::kOverwrite);
        if (histoGammaCaloTruePurity_Pt)                        histoGammaCaloTruePurity_Pt->Write("GammaCaloPurityWOSec_Pt", TObject::kOverwrite);
        // photon reconstruction efficiency including resolution correction
        if (histoGammaPrimaryRecoEff_Pt)                        histoGammaPrimaryRecoEff_Pt->Write("GammaRecoEff_WithResolCorr_Pt", TObject::kOverwrite);
        if (histoGammaCaloPrimaryRecoEff_Pt)                    histoGammaCaloPrimaryRecoEff_Pt->Write("GammaCaloRecoEff_WithResolCorr_Pt", TObject::kOverwrite);
        // photon reconstruction efficiency without resolution correction
        if (histoGammaPrimaryRecoEff_MCPt)                      histoGammaPrimaryRecoEff_MCPt->Write("GammaRecoEff_MCPt", TObject::kOverwrite);
        if (histoGammaCaloPrimaryRecoEff_MCPt)                  histoGammaCaloPrimaryRecoEff_MCPt->Write("GammaCaloRecoEff_MCPt", TObject::kOverwrite);
        // photon conversion probability
        if (histoGammaConvProb_MCPt)                            histoGammaConvProb_MCPt->Write("GammaConvProb_Pt", TObject::kOverwrite);
        // resolution correction in case of unfolding
        if (histoGammaResolCorrUnfold_Pt)                       histoGammaResolCorrUnfold_Pt->Write("GammaResolCorrUnfold_Pt", TObject::kOverwrite);
        // photon correction factors (conv Prob, efficiency incl. resolution correction)
        if (histoGammaCorrFac_Pt)                               histoGammaCorrFac_Pt->Write("GammaCorrFac_Pt", TObject::kOverwrite);

        // ________________________ writing data quantities to file
        // uncorrected spectrum (scaled by 1/Nevt)
        if (histoGammaRawSpectrum_Pt)                           histoGammaRawSpectrum_Pt->Write("GammaRaw_Pt", TObject::kOverwrite);
        // corrected spectrum (legacy corrections: purity, effi incl resolution correction, secondaries, conv prob, plus trivial factors  )
        if (histoGammaCorrEffiReso_Pt)                          histoGammaCorrEffiReso_Pt->Write("GammaCorrEffiResol_Pt", TObject::kOverwrite);
        // corrected spectrum (unfolding corrections: purity, secondaries, unfolding resolution correction, effi without unfolding corr, conv prob, plus trivial factors )
        if (histoGammaCorrUnfoldReso_Pt)                        histoGammaCorrUnfoldReso_Pt->Write("GammaCorrUnfold_Pt", TObject::kOverwrite);
        if(doPileUpCorr){
            // pileup correction factor
            if (histoPileUpCorrectionFactor)                    histoPileUpCorrectionFactor->Write("PileUpCorrectionFactor", TObject::kOverwrite);
            // same as histoGammaCorrEffiReso_Pt with additional pileup correction
            if (histoGammaCorrEffiReso_PileUp_Pt)               histoGammaCorrEffiReso_PileUp_Pt->Write("GammaCorrEffiResolPileup_Pt", TObject::kOverwrite);
            // -> original binning
            if (histoPileUpCorrectionFactorOrBin)               histoPileUpCorrectionFactorOrBin->Write("PileUpCorrectionFactorOrBin",  TObject::kOverwrite);
            if (histoRatioWithWithoutPileUp)                    histoRatioWithWithoutPileUp->Write(     "RatioWithWithoutPileUp",       TObject::kOverwrite);
            if (histoRatioWithWithoutPileUpFit)                 histoRatioWithWithoutPileUpFit->Write(  "RatioWithWithoutPileUpFit",    TObject::kOverwrite);
        }
        
    fileCorrectedOutput->Close();
}
