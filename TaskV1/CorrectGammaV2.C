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
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"

//*********************************************************************************************************
//*********************************************************************************************************
//*********************************************************************************************************
void CorrectGammaEffiResol( TH1D* histoGammaCorr, 
                            TH1D *histoAllSec, 
                            TH1D *histoK0sSec, 
                            TH1D* histoPurity, 
                            TH1D* histoConvProb, 
                            TH1D* histoRecoEff,
                            Double_t deltaEta, 
                            Double_t scaling, 
                            Double_t nEvt
                          ){
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

//*********************************************************************************************************
//*** Correct inclusive raw gamma candite spectrum from conversions with bin-by-bin efficiencies **********
// - scale with 1/nEvt first
// - subtract secodary gamma's 
// - take out impurities
// - conversion probability (histoConvProb)
// - reco efficiency (histoRecoEff)
// - scale with 1/deltaEta, scaling, 1/nEvt
// - convert to inv yield (scale by bin 1/bin center)
//*********************************************************************************************************
void CorrectGammaEffiResolCocktail( TH1D* histoGammaCorr, 
                                    TH1D** histoHadronSec, 
                                    TH1D* histoRestSec, 
                                    TH1D* histoPurity, 
                                    TH1D* histoConvProb, 
                                    TH1D* histoRecoEff, 
                                    Double_t deltaEta, 
                                    Double_t scaling,
                                    Double_t nEvt
                                  ){
    // scale with 1/nEvt as secondaries are normalized per event
    histoGammaCorr->Scale(1./nEvt);
    // subtract secondary gamma's from decays
    for (Int_t k = 0; k < 3; k++){
        if (histoHadronSec[k]) histoGammaCorr->Add(histoHadronSec[k],-1);
    }
    // subtract secondary gamma's from mat interactions
    histoGammaCorr->Add(histoRestSec,-1);
    // Multiply with purity for primary particles
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
    // divide by P_{conv}
    histoGammaCorr->Divide(histoGammaCorr,histoConvProb,1.,1.,"");
    // divide by reconstruction efficiency
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    // scale to with 1/deltaEta and 1/deltaPhi
    histoGammaCorr->Scale(1./deltaEta);
    histoGammaCorr->Scale(scaling);
    // convert to inv yield by dividing by pT at bin center
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoGammaCorr->GetBinContent(i)/histoGammaCorr->GetBinCenter(i);
        Double_t newBinError = histoGammaCorr->GetBinError(i)/histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,newBinError);
    }
}

//*********************************************************************************************************
//*********************************************************************************************************
//*********************************************************************************************************
void CorrectGammaEffiResol( TH1D* histoGammaCorr, 
                            TH1D *histoAllSec, 
                            TH1D *histoK0sSec, 
                            TH1D* histoPurity, 
                            TH1D* histoRecoEff, 
                            Double_t deltaEta, 
                            Double_t scaling, 
                            Double_t nEvt
                          ){
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

//*********************************************************************************************************
//*********************************************************************************************************
//*********************************************************************************************************
void CorrectGammaEffiResolCocktail( TH1D* histoGammaCorr, 
                                    TH1D** histoHadronSec, 
                                    TH1D* histoRestSec, 
                                    TH1D* histoPurity, 
                                    TH1D* histoRecoEff, 
                                    Double_t deltaEta, 
                                    Double_t scaling, 
                                    Double_t nEvt
                                  ){
    histoGammaCorr->Scale(1./nEvt);
    for (Int_t k = 0; k < 3; k++){
        if (histoHadronSec[k]) histoGammaCorr->Add(histoHadronSec[k],-1);
    }    
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

//*********************************************************************************************************
//*********************************************************************************************************
//*********************************************************************************************************
void CorrectGammaSecAndPurity(  TH1D* histoGammaCorr, 
                                TH1D* histoSecGamma, 
                                TH1D* histoSecGammaAddK0s, 
                                TH1D* histoPurity
                             ){
    histoGammaCorr->Sumw2();
    histoGammaCorr->Add(histoSecGamma,-1);
    histoGammaCorr->Add(histoSecGammaAddK0s,-1);
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
}

//*********************************************************************************************************
//*********************************************************************************************************
//*********************************************************************************************************
void CorrectGammaSecAndPurityCocktail(  TH1D* histoGammaCorr, 
                                        TH1D** histoSecGammaHadrons,
                                        TH1D* histoSecGammaRest, 
                                        TH1D* histoPurity
                                     ){
    histoGammaCorr->Sumw2();
    for (Int_t k = 0; k < 3; k++){
        if (histoSecGammaHadrons[k])histoGammaCorr->Add(histoSecGammaHadrons[k],-1);
    }    
    histoGammaCorr->Add(histoSecGammaRest,-1);
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
}

//*********************************************************************************************************
//*********** Correct inclusive raw gamma spectrum from conversions using unfolded as inputs **************
// - conversion probability (histoConvProb)
// - reco efficiency (histoRecoEff)
// - scale with 1/deltaEta, scaling, 1/nEvt
// - convert to inv yield (scale by bin 1/bin center)
//*********************************************************************************************************
void CorrectGammaUnfoldResol(   TH1D* histoGammaCorr, 
                                TH1D* histoConvProb, 
                                TH1D* histoRecoEff, 
                                Double_t deltaEta, 
                                Double_t scaling, 
                                Double_t nEvt
                            ){
    // scale with 1/reco eff
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    // scale with 1/ P_{conv}
    histoGammaCorr->Divide(histoGammaCorr,histoConvProb,1.,1.,"");
    // normalize to rapidity window
    histoGammaCorr->Scale(1./deltaEta);
    // normalize to eta window
    histoGammaCorr->Scale(scaling);
    // normalize per event
    histoGammaCorr->Scale(1./nEvt);
    // convert to inv yield (divide each yield by 1/pT)
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoGammaCorr->GetBinContent(i)/histoGammaCorr->GetBinCenter(i);
        Double_t newBinError = histoGammaCorr->GetBinError(i)/histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,newBinError);
    }
}

//*********************************************************************************************************
//*********************************************************************************************************
//*********************************************************************************************************
void CorrectGammaUnfoldResol(   TH1D* histoGammaCorr, 
                                TH1D* histoRecoEff, 
                                Double_t deltaEta, 
                                Double_t scaling,
                                Double_t nEvt
                            ){
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

//*********************************************************************************************************
//*********************************************************************************************************
//*********************************************************************************************************
Bool_t ConvertCocktailSecondaryToRaw(   TH1D* histoGammaSec, 
                                        TH1D* histoConvProb, 
                                        TH1D* histoRecoEff, 
                                        TH2D* responseMatrix, 
                                        Double_t nEvt, 
                                        Bool_t useResponseMatrix, 
                                        Int_t nIterationsUnfolding = 5
                                    ){
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

//*********************************************************************************************************
//*********************************************************************************************************
//*********************************************************************************************************
Bool_t ConvertCocktailSecondaryToRaw(   TH1D* histoGammaSec, 
                                        TH1D* histoRecoEff, 
                                        TH2D* responseMatrix, 
                                        Double_t nEvt, 
                                        Bool_t useResponseMatrix, 
                                        Int_t nIterationsUnfolding = 5
                                    ){
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

//*********************************************************************************************************
//*********************************************************************************************************
//*********************************************************************************************************
void CorrectGammaMC(    TH1D* histoMCGammaSpec_MCPtCorr,  
                        Double_t deltaEta, 
                        Double_t scaling, 
                        Double_t nEvtMC
                   ){
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
                        TString cutSelection                = "", 
                        TString suffix                      = "eps", 
                        TString nameMeson                   = "Pi0", 
                        TString isMC                        = "kFALSE",
                        TString energy                      = "", 
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

    TString nameRec;
    Bool_t isRunMC                  = kFALSE;
    if (isMC.CompareTo("kTRUE")==0){
        nameRec                     = "MC";
        isRunMC                     = kTRUE;
    } else {
        nameRec                     = "data";
    }

    //*****************************************************************************************
    //*************************** Determine K0s scaling factor ********************************
    //*****************************************************************************************
    Double_t doubleAddFactorK0s     = ReturnCorrectK0ScalingFactor( energy,  cutSelection);
    if(isRunMC){
        doubleAddFactorK0s = 0.;
        cout<<" MC --> doubleAddFactorK0s = 0"<<endl;
    } else {
        cout << "The additional K0 correction factor is: "  << doubleAddFactorK0s<<endl;
    }
    Double_t scaleFactorsSec[4]     = {1.+doubleAddFactorK0s, 1., 1., 1.};
    
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
    TString collisionSystem         = ReturnFullCollisionsSystem(energy);
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
    if(energy.Contains("PbPb")){
        cent                        = Form("%s %s", centrality.Data(), collisionSystem.Data());
    } else {
        cent                        = collisionSystem;
    }
    
    //******************************************************************************************
    //***************************** Set output directory ***************************************
    //******************************************************************************************
    TString outputDir               = Form("%s/%s/%s/CorrectGammaV2",cutSelection.Data(),energy.Data(),suffix.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    //******************************************************************************************
    //************************ Define additional global variables ******************************
    //******************************************************************************************
    TString textPi0New              = Form("Gamma_%s",nameMeson.Data());
    
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
    
    //****************************************************************************************** 
    //********************************** secondary labels **************************************
    //****************************************************************************************** 
    TString nameSecondaries[4]                              = { "K0s", "Lambda", "K0l", "Rest" };
    TString nameLabelSecondaries[4]                         = { "K^{0}_{s}", "#Lambda", "K^{0}_{l}", "rest" };
    Style_t markerStyleSec[4]                               = { 21, 33, 29, 34};
    Style_t markerStyleSecWithToy[4]                        = { 25, 27, 30, 28};
    Size_t markerSizeSec[4]                                 = { 1.5, 1.75, 2., 1.5};
    Color_t colorSecFromToy[4]                              = { kRed-2, kOrange+1, kCyan-2, kBlue+2};

    //****************************************************************************************** 
    //************************ True Combinatorial Background labels ****************************
    //****************************************************************************************** 
    TString combinatorics[17]                               = { "Elec+Elec","Elec+Pion","Elec+Kaon","Elec+Proton","Elec+Muon","Pion+Pion","Pion+Kaon","Pion+Proton",
                                                                "Pion+Muon","Kaon+Kaon","Kaon+Proton","Kaon+Muon","Proton+Proton","Proton+Muon","Muon+Muon","Rest","All"};
    TString combinatoricsLabels[17]                         = { "e^{#pm}e^{#mp}","e^{#pm}#pi^{#mp}","e^{#pm}K^{#mp}","e^{#pm}#bar{p}(p)","e^{#pm}#mu^{#mp}","#pi^{#pm}#pi^{#mp}","#pi^{#pm}K^{#mp}",
                                                                "#pi^{#pm}#bar{p}(p)","#pi^{#pm}#mu^{#mp}","K^{#pm}K^{#mp}","K^{#pm}#bar{p}(p)","K^{#pm}#mu^{#mp}","p(#bar{p})#bar{p}(p)","p(#bar{p})#mu^{#mp}","#mu^{#pm}#mu^{#mp}","rest","all"};
    Color_t colorsCombinatorics[17]                         = { kAzure-1, kRed+1, kOrange+7, kMagenta+2, kRed-9,
                                                                kBlue-6, kAzure+5, kPink, kCyan-3, kGreen-3,
                                                                kSpring+9, kGreen+2, kBlue+2, kMagenta-6, kSpring+4,
                                                                kCyan+2, 809};
    Style_t markersCombinatorics[17]                        = { 20, 21, 24, 25, 27,
                                                                28, 29, 30, 33, 34, 
                                                                20, 21, 24, 25, 27,
                                                                28, 29};
    TString combinatoricsCalo[11]                           = { "Electron","Pion","Proton","Kaon","Neutron","K0s","Lambda","Muon","K0l","Rest","All"};
    TString combinatoricsLabelsCalo[11]                     = { "e^{#pm}","#pi^{#pm}","p(#bar{p})","K^{#pm}","n","K^{0}_{s}","#Lambda","#mu^{#pm}","K^{0}_{l}","rest","all"};
    Color_t colorsCombinatoricsCalo[11]                     = { kAzure-1, kRed+1, kOrange+7, kMagenta+2, kRed-9, kBlue,
                                                                kBlue-6, kAzure+5, kPink, kCyan-3, kGreen-3};
    Style_t markersCombinatoricsCalo[11]                    = { 20, 21, 24, 25, 27,
                                                                28, 29, 30, 33, 34, 20};    
                                                                
    //****************************************************************************************** 
    //********************************** decay labels ******************************************
    //******************************************************************************************                                                                 
    TString decays[9]                                       = { "Pi0","Eta","Etap","Omega","Rho",
                                                                "Phi","Sigma", "All decays","direct #gamma"};
    TString decaysLabels[9]                                 = { "#gamma from #pi^{0}","#gamma from #eta","#gamma from #eta'","#gamma from #omega","#gamma from #rho^{0}",
                                                                "#gamma from #phi", "#gamma from #Sigma^{0}", "All decays","direct #gamma"};
    Color_t colorsDecay[9]                                  = { kRed+1, kAzure-1, kOrange+7, kGreen+2, kCyan-3, 809, 
                                                                kBlue+2, kBlack, kRed-1};
    
    //****************************************************************************************** 
    //*******************File definitions and reading histograms from files ********************
    //******************************************************************************************
    TFile*  fileUnCorrected                                         = new TFile(nameUnCorrectedFile);
    TFile*  fileCorrections                                         = new TFile(nameCorrectionFile);
    if (!fileUnCorrected) {
        cout << "file " << nameUnCorrectedFile << " not found!" << endl;
        return;
    }
    if (!nameCorrectionFile) {
        cout << "file " << nameCorrectionFile << " not found!" << endl;
        return;
    }

    //******************* Conversion gamma spectra *********************************************
    TH1D*   histoESDConvGammaPt                                     = NULL;
    TH1D*   histoESDConvGammaPt_OriginalBin                         = NULL;
    if (isPCM) {
        histoESDConvGammaPt                                         = (TH1D*)fileUnCorrected->Get("ESD_ConvGamma_Pt");
        histoESDConvGammaPt_OriginalBin                             = (TH1D*)fileUnCorrected->Get("ESD_ConvGamma_Pt_OriginalBinning");
    }

    //******************* Calo gamma spectra ***************************************************
    TH1D*   histoESDCaloGammaPt                                     = NULL;
    TH1D*   histoESDCaloGammaPt_OriginalBin                         = NULL;
    if (isCalo) {
        histoESDCaloGammaPt                                         = (TH1D*)fileUnCorrected->Get("ESD_CaloGamma_Pt");
        histoESDCaloGammaPt_OriginalBin                             = (TH1D*)fileUnCorrected->Get("ESD_CaloGamma_Pt_OriginalBinning");
    }

    //******************* Determine max Pt *****************************************************
    Double_t maxPtGamma                                             = 0.;
    if (isPCM && !isCalo) maxPtGamma                                = histoESDConvGammaPt->GetXaxis()->GetBinUpEdge(histoESDConvGammaPt->GetNbinsX());
    if (isCalo && !isPCM) maxPtGamma                                = histoESDCaloGammaPt->GetXaxis()->GetBinUpEdge(histoESDCaloGammaPt->GetNbinsX());

    //******************* Calculate number of events for normalization *************************
    TH1D*   histoEventQuality                                       = (TH1D*)fileUnCorrected->Get("NEvents");
    Float_t nEvt                                                    = 0.;
    if (energy.CompareTo("PbPb_2.76TeV") == 0)  nEvt                = histoEventQuality->GetBinContent(1);
    else                                        nEvt                = GetNEvents(histoEventQuality);
    
    TH1F*   histoEventQualityMC                                     = (TH1F*)fileCorrections->Get("NEvents");
    Float_t nEvtMC                                                  = 0.;
    if (energy.CompareTo("PbPb_2.76TeV") == 0)  nEvtMC              = histoEventQualityMC->GetBinContent(1);
    else                                        nEvtMC              = GetNEvents(histoEventQualityMC);
    
    //******************* Binning histogram ****************************************************
    TH1D*   deltaPt                                                 = (TH1D*)fileUnCorrected->Get("deltaPt");
    for (Int_t i = 0; i < deltaPt->GetNbinsX() +1; i++) deltaPt->SetBinError(i, 0);

    //******************* Primary gamma correction factors (PCM and Calo) **********************
    TH1D*   histoGammaPurity_Pt                                     = NULL;
    TH1D*   histoGammaTruePurity_Pt                                 = NULL;
    TH1D*   histoGammaTruePurity_OriginalBin_Pt                     = NULL;
    TH1D*   histoGammaPrimaryRecoEff_Pt                             = NULL;
    TH1D*   histoGammaPrimaryRecoEff_Pt_OrBin                       = NULL;
    TH1D*   histoGammaPrimaryRecoEff_MCPt                           = NULL;
    TH1D*   histoGammaPrimaryRecoEff_MCPt_OrBin                     = NULL;
    if ( isPCM || isCalo ) {
        // reconstructed validated photons / all reconstructed photon candidates in MC
        histoGammaPurity_Pt                                         = (TH1D*)fileCorrections->Get("GammaPurity_Pt");
        // reconstructed validated primary photons / (all reconstructed photon candidates in MC - validated reconstructed secondary photons)
        // rescaling of histoGammaPurity_Pt for primary photons
        histoGammaTruePurity_Pt                                     = (TH1D*)fileCorrections->Get("GammaTruePurity_Pt");
        // reconstructed validated primary photons / (all reconstructed photon candidates in MC - validated reconstructed secondary photons)
        // rescaling of histoGammaPurity_Pt for primary photons
        // rebinned
        histoGammaTruePurity_OriginalBin_Pt                         = (TH1D*)fileCorrections->Get("GammaTruePurity_OriginalBinning_Pt");
        // reconstructed validated primary photons vs rec pt / converted MC photons vs MCpt
        // this reconstruction efficiency includes the transverse momentum resolution correction and the acceptance correction
        histoGammaPrimaryRecoEff_Pt                                 = (TH1D*)fileCorrections->Get("GammaPrimaryRecoEff_Pt");
        histoGammaPrimaryRecoEff_Pt_OrBin                           = (TH1D*)fileCorrections->Get("GammaPrimaryRecoEff_Pt_OriginalBinning");
        // reconstructed validated primary photons vs MC pt / converted MC photons vs MCpt
        // this reconstruction efficiency includes the acceptance correction
        histoGammaPrimaryRecoEff_MCPt                               = (TH1D*)fileCorrections->Get("GammaPrimaryRecoEff_MCPt");
        histoGammaPrimaryRecoEff_MCPt_OrBin                         = (TH1D*)fileCorrections->Get("GammaPrimaryRecoEff_MCPt_OriginalBinning");
    }
    TH1D*   histoGammaConvProb_MCPt                                 = NULL;
    TH1D*   histoGammaConvProb_MCPt_OrBin                           = NULL;
    if ( isPCM ) {
        // converted photons MC/ all generated photons (always based on primaries only)
        histoGammaConvProb_MCPt                                     = (TH1D*)fileCorrections->Get("GammaConvProb_MCPt");
        histoGammaConvProb_MCPt_OrBin                               = (TH1D*)fileCorrections->Get("GammaConvProb_MCPt_OriginalBinning");
    }
    
    //******************* Primary gamma correction factors (PCM-Calo hybrid) *******************
    TH1D*   histoGammaCaloPurity_Pt                                 = NULL;
    TH1D*   histoGammaCaloTruePurity_Pt                             = NULL;
    TH1D*   histoGammaCaloTruePurity_OrigBin_Pt                     = NULL;
    TH1D*   histoGammaCaloPrimaryRecoEff_Pt                         = NULL;
    TH1D*   histoGammaCaloPrimaryRecoEff_MCPt                       = NULL;
    // Calorimeter photons 
    if ( isPCM && isCalo ) {
        // reconstructed validated photons / all reconstructed photon candidates in MC
        histoGammaCaloPurity_Pt                                     = (TH1D*)fileCorrections->Get("GammaCaloPurity_Pt");
        // reconstructed validated primary photons / (all reconstructed photon candidates in MC - validated reconstructed secondary photons)
        // rescaling of histoGammaPurity_Pt for primary photons
        histoGammaCaloTruePurity_Pt                                 = (TH1D*)fileCorrections->Get("GammaCaloTruePurity_Pt");
        // reconstructed validated primary photons / (all reconstructed photon candidates in MC - validated reconstructed secondary photons)
        // rescaling of histoGammaPurity_Pt for primary photons
        // rebinned
        histoGammaCaloTruePurity_OrigBin_Pt                         = (TH1D*)fileCorrections->Get("GammaCaloTruePurity_OriginalBinning_Pt");
        // reconstructed validated primary photons vs rec pt / converted MC photons vs MCpt 
        // this reconstruction efficiency includes the transverse momentum resolution correction and the acceptance correction    
        histoGammaCaloPrimaryRecoEff_Pt                             = (TH1D*)fileCorrections->Get("GammaCaloPrimaryRecoEff_Pt");
        // reconstructed validated primary photons vs MC pt / converted MC photons vs MCpt 
        // this reconstruction efficiency includes the acceptance correction        
        histoGammaCaloPrimaryRecoEff_MCPt                           = (TH1D*)fileCorrections->Get("GammaCaloPrimaryRecoEff_MCPt");
    }

    //******************* MC histograms (PCM and Calo) *****************************************
    TH1D*   histoMCrecGamma_Pt                                      = NULL;
    TH1D*   histoMCrecGamma_OriginalBin_Pt                          = NULL;
    TH1D*   histoMCrecPrimaryGamma_Pt                               = NULL;
    TH1D*   histoMCGammaConv_MCPt                                   = NULL;
    TH1D*   histoGammaTrueConv_Pt                                   = NULL;
    TH1D*   histoGammaTruePrimaryConv_Pt                            = NULL;
    TH2D*   histoGammaTruePrimaryConv_recPt_MCPt                    = NULL;
    if (isPCM) {
        // reconstructed photon candidates in MC
        histoMCrecGamma_Pt                                          = (TH1D*)fileCorrections->Get("MCrec_ConvGamma_Pt");
        histoMCrecGamma_OriginalBin_Pt                              = (TH1D*)fileCorrections->Get("MCrec_ConvGamma_OriginalBinning_Pt");
        // reconstructed photon candidates in MC - validated secondary photons
        histoMCrecPrimaryGamma_Pt                                   = (TH1D*)fileCorrections->Get("MCrec_PrimaryConvGamma_Pt");
        // MC converted photons regardless of source
        histoMCGammaConv_MCPt                                       = (TH1D*)fileCorrections->Get("MC_ConvGamma_MCPt");
        // validated reconstructed photons
        histoGammaTrueConv_Pt                                       = (TH1D*)fileCorrections->Get("TrueConvGamma_Pt");
        // validated reconstructed primary photons
        histoGammaTruePrimaryConv_Pt                                = (TH1D*)fileCorrections->Get("TruePrimaryConvGamma_Pt");
        // response matrix for photons recPt vs MC pt
        histoGammaTruePrimaryConv_recPt_MCPt                        = (TH2D*)fileCorrections->Get("TruePrimaryConvGamma_recPt_MCPt");
        // response matrix for photons recPt vs MC pt rebinned
    }
    TH1D*   histoMCrecBackground_Pt                                 = NULL;
    TH1D*   histoMCAllGamma_MCPt                                    = NULL;
    TH1D*   histoMCAllGamma_OriginalBin_MCPt                        = NULL;
    if (isPCM ||isCalo) {
        // reconstructed background in MC
        histoMCrecBackground_Pt                                     = (TH1D*)fileCorrections->Get("MCrec_Background");
        // MC input for all gamma, regardless of source
        histoMCAllGamma_MCPt                                        = (TH1D*)fileCorrections->Get("MC_AllGamma_MCPt");
        histoMCAllGamma_OriginalBin_MCPt                            = (TH1D*)fileCorrections->Get("MC_AllGamma_OriginalBinning_MCPt");
    }

    //******************* MC histograms (Calo) *************************************************
    TH1D*   histoMCrecCaloBackground_Pt                             = NULL;
    TH1D*   histoMCAllGammaCalo_MCPt                                = NULL;
    TH1D*   histoMCAllGammaCalo_OriginalBin_MCPt                    = NULL;
    if (isPCM && isCalo) {
        // reconstructed background in MC
        histoMCrecCaloBackground_Pt                                 = (TH1D*)fileCorrections->Get("MCrec_Calo_Background");
        // MC input for all gamma, regardless of source
        histoMCAllGammaCalo_MCPt                                    = (TH1D*)fileCorrections->Get("MC_AllGammaEMCAcc_MCPt");
        histoMCAllGammaCalo_OriginalBin_MCPt                        = (TH1D*)fileCorrections->Get("MC_AllGammaEMCAcc_OriginalBinning_MCPt");
    }
    TH1D*   histoMCrecGammaCalo_Pt                                  = NULL;
    TH1D*   histoMCrecGammaCalo_OriginalBin_Pt                      = NULL;
    TH1D*   histoMCrecPrimaryGammaCalo_Pt                           = NULL;
    TH1D*   histoGammaTrueCalo_Pt                                   = NULL;
    TH1D*   histoGammaTruePrimaryCalo_Pt                            = NULL;
    TH2D*   histoGammaTruePrimaryCalo_recPt_MCPt                    = NULL;
    if (isCalo) {
        // reconstructed photon candidates in MC
        histoMCrecGammaCalo_Pt                                      = (TH1D*)fileCorrections->Get("MCrec_CaloGamma_Pt");
        histoMCrecGammaCalo_OriginalBin_Pt                          = (TH1D*)fileCorrections->Get("MCrec_CaloGamma_OriginalBinning_Pt");
        // reconstructed photon candidates in MC - validated secondary photons
        histoMCrecPrimaryGammaCalo_Pt                               = (TH1D*)fileCorrections->Get("MCrec_PrimaryCaloGamma_Pt");
        // validated reconstructed photons
        histoGammaTrueCalo_Pt                                       = (TH1D*)fileCorrections->Get("TrueCaloGamma_Pt");
        // validated reconstructed primary photons
        histoGammaTruePrimaryCalo_Pt                                = (TH1D*)fileCorrections->Get("TruePrimaryCaloGamma_Pt");
        // response matrix for photons recPt vs MC pt (check for binning)
        histoGammaTruePrimaryCalo_recPt_MCPt                        = (TH2D*)fileCorrections->Get("TruePrimaryCaloGamma_recPt_MCPt");
    }

    //******************* MC decay gammas ******************************************************
    TH1D**  histoPhotonSource_MCPt                                  = new TH1D*[9];
    for (Int_t i = 0;i<7;i++) {
        histoPhotonSource_MCPt[i]                                   = (TH1D*)fileCorrections->Get(Form("MC_DecayGamma%s_Pt",decays[i].Data()));
        histoPhotonSource_MCPt[i]->Sumw2();
        if (i == 0) histoPhotonSource_MCPt[7]                       = (TH1D*)histoPhotonSource_MCPt[0]->Clone("MC_DecayGammaAll_Pt");
        else        histoPhotonSource_MCPt[7]->Add(histoPhotonSource_MCPt[i]);
    }
    histoPhotonSource_MCPt[8]                                       = (TH1D*)histoMCAllGamma_OriginalBin_MCPt->Clone("MC_DirectPhotons");
    histoPhotonSource_MCPt[8]->Sumw2();
    histoPhotonSource_MCPt[8]->Add(histoPhotonSource_MCPt[7],-1);

    //******************* MC combinatorial gammas (PCM) ****************************************
    TH1D**  histoCombinatorialSpecies_Pt                            = NULL;
    if (isPCM) {
        histoCombinatorialSpecies_Pt                                = new TH1D*[17];
        for(Int_t i = 0;i<17;i++){
            histoCombinatorialSpecies_Pt[i]                         = (TH1D*)fileCorrections->Get(Form("ESD_TrueComb%s_Pt",combinatorics[i].Data()));
            histoCombinatorialSpecies_Pt[i]->SetMinimum(1e-10);
        }
    }

    //******************* MC combinatorial gammas (Calo) ***************************************
    TH1D**  histoCombinatorialSpeciesCalo_Pt                        = NULL;
    if (isCalo) {
        histoCombinatorialSpeciesCalo_Pt                            = new TH1D*[11];
        for(Int_t i = 0;i<11;i++){
            histoCombinatorialSpeciesCalo_Pt[i]                     = (TH1D*)fileCorrections->Get(Form("ESD_TrueComb%s_Pt",combinatoricsCalo[i].Data()));
            histoCombinatorialSpeciesCalo_Pt[i]->SetMinimum(1e-10);
        }
    }

    //******************* MC true secondary gammas (PCM) ***************************************
    TH1D*   histoGammaTrueSecConvGammaFromX_Pt[4]                   = {NULL, NULL, NULL, NULL};
    TH1D*   histoGammaTrueSecConvGammaFromX_PtOrgBin[4]             = {NULL, NULL, NULL, NULL};
    if (isPCM) {
        for (Int_t k = 0; k < 4; k++){
            histoGammaTrueSecConvGammaFromX_Pt[k]                   = (TH1D*)fileCorrections->Get(Form("TrueSecondaryConvGammaFromXFrom%s_Pt", nameSecondaries[k].Data()));
            histoGammaTrueSecConvGammaFromX_PtOrgBin[k]             = (TH1D*)fileCorrections->Get(Form("TrueSecondaryConvGammaFromXFrom%s_Pt_OriginalBinning", nameSecondaries[k].Data()));
        }
    }
    
    //******************* MC true secondary gammas (Calo) **************************************
    TH1D*   histoGammaTrueSecCaloGammaFromX_Pt[4]                   = {NULL, NULL, NULL, NULL};
    TH1D*   histoGammaTrueSecCaloGammaFromX_PtOrgBin[4]             = {NULL, NULL, NULL, NULL};
    if (isCalo) {
        for (Int_t k = 0; k < 4; k++){
            histoGammaTrueSecCaloGammaFromX_Pt[k]                   = (TH1D*)fileCorrections->Get(Form("TrueSecondaryCaloGammaFromXFrom%s_Pt", nameSecondaries[k].Data()));
            histoGammaTrueSecCaloGammaFromX_PtOrgBin[k]             = (TH1D*)fileCorrections->Get(Form("TrueSecondaryCaloGammaFromXFrom%s_Pt_OriginalBinning", nameSecondaries[k].Data()));
        }    
    }
    
    //******************* MC true secondary gamma fractions (PCM or Calo) **********************
    TH1D* histoFracAllGammaToSecFromX_Pt[4]                         = {NULL, NULL, NULL, NULL};
    TH1D* histoFracAllGammaToSecFromX_OrBin_Pt[4]                   = {NULL, NULL, NULL, NULL};
    if ( isPCM || isCalo ) {
        for (Int_t k = 0; k<4; k++){
            histoFracAllGammaToSecFromX_Pt[k]                       = (TH1D*)fileCorrections->Get(Form("FracAllGammaToSecFromXFrom%s",nameSecondaries[k].Data()));
            histoFracAllGammaToSecFromX_OrBin_Pt[k]                 = (TH1D*)fileCorrections->Get(Form("FracAllGammaToSecFromXFrom%sOriginalBinning",nameSecondaries[k].Data()));
        }
    }
    
    //****************************************************************************************** 
    //******************* Load Secondary gamma spectra from cocktail ***************************
    //****************************************************************************************** 
    Bool_t  hasCocktailInput                                       = kTRUE;
    TH1D*   histoGammaTrueSecCocktailGammaFromX_Pt[4]              = { NULL, NULL, NULL, NULL};
    TH1D*   histoGammaTrueSecCocktailGammaFromX_PtOrBin[4]         = { NULL, NULL, NULL, NULL};
    
    if( (isPCM || isCalo) && !isRunMC ){
        for (Int_t k = 0; k < 4; k++){
            if (k < 3){
                histoGammaTrueSecCocktailGammaFromX_Pt[k]          = (TH1D*)fileUnCorrected->Get(Form("CocktailSecondaryGammaFromX%s_Pt",nameSecondaries[k].Data()));
                histoGammaTrueSecCocktailGammaFromX_PtOrBin[k]     = (TH1D*)fileUnCorrected->Get(Form("CocktailSecondaryGammaFromX%s_PtOrBin",nameSecondaries[k].Data()));
            } else if ( isPCM && !isCalo ) {
                histoGammaTrueSecCocktailGammaFromX_Pt[k]          = (TH1D*)fileCorrections->Get("TrueSecondaryConvGammaFromXFromRest_Pt");
                histoGammaTrueSecCocktailGammaFromX_PtOrBin[k]     = (TH1D*)fileCorrections->Get("TrueSecondaryConvGammaFromXFromRest_Pt_OriginalBinning");            
            } else if ( isCalo && !isPCM ) {
                histoGammaTrueSecCocktailGammaFromX_Pt[k]          = (TH1D*)fileCorrections->Get("TrueSecondaryCaloGammaFromXFromRest_Pt");
                histoGammaTrueSecCocktailGammaFromX_PtOrBin[k]     = (TH1D*)fileCorrections->Get("TrueSecondaryCaloGammaFromXFromRest_Pt_OriginalBinning");            
            }
            if (!histoGammaTrueSecCocktailGammaFromX_Pt[k] || !histoGammaTrueSecCocktailGammaFromX_PtOrBin[k])
                hasCocktailInput                                   = kFALSE;
        }    
        if (!hasCocktailInput) cout << "secondary spectra from cocktail not found, will not use" << endl;
    } else {
        hasCocktailInput                                           = kFALSE;
    }    

    //****************************************************************************************** 
    //******************* Secondary gamma reco eff *********************************************
    //****************************************************************************************** 
    TH1D* histoGammaSecFromXRecoEff_MCPt[3]                         = { NULL, NULL, NULL };
    TH1D* histoGammaSecFromXRecoEff_MCPt_Unscaled[3]                = { NULL, NULL, NULL };
    TH1D* histoGammaSecFromXRecoEff_MCPtOrBin[3]                    = { NULL, NULL, NULL };
    TH1D* histoGammaSecFromXRecoEff_RecPt[3]                        = { NULL, NULL, NULL };
    TH1D* histoGammaSecFromXRecoEff_RecPt_Unscaled[3]               = { NULL, NULL, NULL };
    TH1D* histoGammaSecFromXRecoEff_RecPtOrBin_Unscaled[3]          = { NULL, NULL, NULL };
    TH1D* histoGammaSecFromXRecoEff_RecPtOrBin[3]                   = { NULL, NULL, NULL };
    TH1D* ratioGammaSecEffMCPt[3]                                   = { NULL, NULL, NULL };
    TH1D* ratioGammaSecEffRecPt[3]                                  = { NULL, NULL, NULL };
    Double_t constOffsetEffMCPt[3]                                  = { 1, 1, 1};
    Double_t constOffsetEffRecPt[3]                                 = { 1, 1, 1};
    Double_t minPtFitSec[3]                                         = { 1.5, 1, 0.};
    
    Double_t maxPtFitSec = 50;
    if (energy.Contains("pPb")){
        maxPtFitSec = maxPtGamma;
    }

    if ( hasCocktailInput && (isPCM || isCalo) ) {        
        TF1* constant                                               = new TF1("constant", "[0]", 0, 50);
        // taken directly from MC (statistics sufficient)
        for (Int_t k = 0; k < 3; k++ ){
            histoGammaSecFromXRecoEff_MCPt[k]                       = (TH1D*)fileCorrections->Get(Form("SecondaryGammaFromXFrom%sRecoEff_MCPt", nameSecondaries[k].Data()));
            histoGammaSecFromXRecoEff_MCPt_Unscaled[k]              = (TH1D*)histoGammaSecFromXRecoEff_MCPt[k]->Clone(Form("SecondaryGammaFromXFrom%sRecoEff_MCPt_Unscaled", nameSecondaries[k].Data()));
            histoGammaSecFromXRecoEff_MCPtOrBin[k]                  = (TH1D*)fileCorrections->Get(Form("SecondaryGammaFromXFrom%sRecoEff_MCPtOrBin", nameSecondaries[k].Data()));
            histoGammaSecFromXRecoEff_RecPt[k]                      = (TH1D*)fileCorrections->Get(Form("SecondaryGammaFromXFrom%sRecoEff_Pt", nameSecondaries[k].Data()));
            histoGammaSecFromXRecoEff_RecPt_Unscaled[k]             = (TH1D*)histoGammaSecFromXRecoEff_RecPt[k]->Clone(Form("SecondaryGammaFromXFrom%sRecoEff_Pt_Unscaled", nameSecondaries[k].Data()));
            histoGammaSecFromXRecoEff_RecPtOrBin[k]                 = (TH1D*)fileCorrections->Get(Form("SecondaryGammaFromXFrom%sRecoEff_PtOrBin", nameSecondaries[k].Data()));
            histoGammaSecFromXRecoEff_RecPtOrBin_Unscaled[k]        = (TH1D*)histoGammaSecFromXRecoEff_RecPtOrBin[k]->Clone(Form("SecondaryGammaFromXFrom%sRecoEff_PtOrBin_Unscaled", nameSecondaries[k].Data()));
            ratioGammaSecEffMCPt[k]                                 = (TH1D*)histoGammaSecFromXRecoEff_MCPt[k]->Clone(Form("RatioSecEffFrom%sToPrim_MCPt", nameSecondaries[k].Data()));
            ratioGammaSecEffMCPt[k]->Divide(histoGammaPrimaryRecoEff_MCPt);
            ratioGammaSecEffRecPt[k]                                = (TH1D*)histoGammaSecFromXRecoEff_RecPt[k]->Clone(Form("RatioSecEffFrom%sToPrim_RecPt", nameSecondaries[k].Data()));
            ratioGammaSecEffRecPt[k]->Divide(histoGammaPrimaryRecoEff_Pt);
            
            // constant fit to scale primary reco eff (statistics insufficient) MC pt
            ratioGammaSecEffMCPt[k]->Fit(constant,"SMNR0E+","", minPtFitSec[k], maxPtFitSec);
            constOffsetEffMCPt[k]                                   = constant->GetParameter(0);
            // constant fit to scale primary reco eff (statistics insufficient) MC pt
            ratioGammaSecEffRecPt[k]->Fit(constant,"SMNR0E+","", minPtFitSec[k], maxPtFitSec);
            constOffsetEffRecPt[k]                                  = constant->GetParameter(0);

            for (Int_t i = histoGammaSecFromXRecoEff_RecPt[k]->FindBin(minPtFitSec[k]); i < histoGammaSecFromXRecoEff_RecPt[k]->GetNbinsX()+1; i++){
                histoGammaSecFromXRecoEff_RecPt[k]->SetBinContent(i, histoGammaPrimaryRecoEff_Pt->GetBinContent(i)*constOffsetEffRecPt[k]);
                histoGammaSecFromXRecoEff_RecPt[k]->SetBinError(i, histoGammaPrimaryRecoEff_Pt->GetBinError(i)*constOffsetEffRecPt[k]);
                histoGammaSecFromXRecoEff_MCPt[k]->SetBinContent(i, histoGammaPrimaryRecoEff_MCPt->GetBinContent(i)*constOffsetEffMCPt[k]);
                histoGammaSecFromXRecoEff_MCPt[k]->SetBinError(i, histoGammaPrimaryRecoEff_MCPt->GetBinError(i)*constOffsetEffMCPt[k]);
            }    
            for (Int_t i = histoGammaSecFromXRecoEff_RecPtOrBin[k]->FindBin(minPtFitSec[k]); i < histoGammaSecFromXRecoEff_RecPtOrBin[k]->GetNbinsX()+1; i++){
                histoGammaSecFromXRecoEff_RecPtOrBin[k]->SetBinContent(i, histoGammaPrimaryRecoEff_Pt_OrBin->GetBinContent(i)*constOffsetEffRecPt[k]);
                histoGammaSecFromXRecoEff_RecPtOrBin[k]->SetBinError(i, histoGammaPrimaryRecoEff_Pt_OrBin->GetBinError(i)*constOffsetEffRecPt[k]);
                histoGammaSecFromXRecoEff_MCPtOrBin[k]->SetBinContent(i, histoGammaPrimaryRecoEff_MCPt_OrBin->GetBinContent(i)*constOffsetEffMCPt[k]);
                histoGammaSecFromXRecoEff_MCPtOrBin[k]->SetBinError(i, histoGammaPrimaryRecoEff_MCPt_OrBin->GetBinError(i)*constOffsetEffMCPt[k]);
            }    

            //****************************************************************************************** 
            //********************************* plot ratio efficiencies ********************************
            //****************************************************************************************** 

            TCanvas *canvasSecEffiRatio                               = GetAndSetCanvas("canvasSecEffiRatio");
            
                SetHistogramm(ratioGammaSecEffRecPt[k],"#it{p}_{T} (GeV/#it{c})","#epsilon_{eff,#gamma, sec}/#epsilon_{eff,#gamma, prim}",0.0,3.0);
                DrawGammaSetMarker(ratioGammaSecEffRecPt[k], 20, 1, kRed+2, kRed+2);
                ratioGammaSecEffRecPt[k]->Draw();
                DrawGammaSetMarkerTF1( constant, 9, 2, kRed-6); 
                constant->Draw("same");
                
                TLegend* legendSecEffRatio = GetAndSetLegend2(0.15,0.93-3*1.1*0.035, 0.4,0.93, 0.035, 1, cent, 42, 0.1);
                legendSecEffRatio->AddEntry(ratioGammaSecEffRecPt[k], Form("sec %s reco. eff/prim",nameLabelSecondaries[k].Data()),"lp");
                legendSecEffRatio->AddEntry(constant, Form("const fit: %2.2f", constOffsetEffRecPt[k]),"l");
                legendSecEffRatio->Draw();
                
            canvasSecEffiRatio->SaveAs(Form("%s/%s_RatioSecEffiToPrim%sPt_%s_%s.%s",outputDir.Data(),textPi0New.Data(),nameSecondaries[k].Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));        
            delete canvasSecEffiRatio;
            
            //****************************************************************************************** 
            //************************* plot reco effs and fits for secondaries ************************
            //****************************************************************************************** 
            TCanvas *canvasSecEffiFit                               = GetAndSetCanvas("canvasSecEffiFit");
            
                SetHistogramm(histoGammaPrimaryRecoEff_Pt,"#it{p}_{T} (GeV/#it{c})","#epsilon_{eff,#gamma} in |#eta| < 0.9",0.0,1.02);
                DrawGammaSetMarker(histoGammaPrimaryRecoEff_Pt, 20, 1, kGray+2, kGray+2);
                histoGammaPrimaryRecoEff_Pt->Draw();
                DrawGammaSetMarker(histoGammaPrimaryRecoEff_Pt_OrBin, 24, 1, kBlack, kBlack);
                histoGammaPrimaryRecoEff_Pt_OrBin->Draw("same");
                            
                DrawGammaSetMarker(histoGammaSecFromXRecoEff_RecPt_Unscaled[k], 20, 1, kBlue-8, kBlue-8);
                histoGammaSecFromXRecoEff_RecPt_Unscaled[k]->Draw("same");
                DrawGammaSetMarker(histoGammaSecFromXRecoEff_RecPtOrBin_Unscaled[k], 24, 1, kBlue+2, kBlue+2);
                histoGammaSecFromXRecoEff_RecPtOrBin_Unscaled[k]->Draw("same");
                DrawGammaSetMarker(histoGammaSecFromXRecoEff_RecPt[k], 20, 1, kRed-8, kRed-8);
                histoGammaSecFromXRecoEff_RecPt[k]->Draw("same");
                DrawGammaSetMarker(histoGammaSecFromXRecoEff_RecPtOrBin[k], 20, 1, kRed+2, kRed+2);
                histoGammaSecFromXRecoEff_RecPtOrBin[k]->Draw("same");
                
                TLegend* legendSecEffFits = GetAndSetLegend2(0.15,0.93-7*1.1*0.035, 0.4,0.93, 0.035, 1, cent, 42, 0.1);
                legendSecEffFits->AddEntry(histoGammaPrimaryRecoEff_Pt,"primary reco. eff.","lp");
                legendSecEffFits->AddEntry(histoGammaPrimaryRecoEff_Pt_OrBin,"rebin primary reco. eff.","lp");
                legendSecEffFits->AddEntry(histoGammaSecFromXRecoEff_RecPtOrBin_Unscaled[k], Form("%s reco. eff",nameLabelSecondaries[k].Data()),"lp");
                legendSecEffFits->AddEntry(histoGammaSecFromXRecoEff_RecPt_Unscaled[k], Form("rebin %s reco. eff",nameLabelSecondaries[k].Data()),"lp");
                legendSecEffFits->AddEntry(histoGammaSecFromXRecoEff_RecPt[k], Form("new %s reco. eff",nameLabelSecondaries[k].Data()),"lp");
                legendSecEffFits->AddEntry(histoGammaSecFromXRecoEff_RecPtOrBin[k], Form("new org.bin. %s reco. eff",nameLabelSecondaries[k].Data()),"lp");
                legendSecEffFits->Draw();
            
            canvasSecEffiFit->SaveAs(Form("%s/%s_SecEffiFits%sPt_%s_%s.%s",outputDir.Data(),textPi0New.Data(),nameSecondaries[k].Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));        
            delete canvasSecEffiFit;
        }
        delete constant;
    }
    
    //****************************************************************************************** 
    //******************* Secondary gamma conv prob ********************************************
    //****************************************************************************************** 
    TH1D* histoGammaSecondaryFromXConvProb_MCPt[3]                  = { NULL, NULL, NULL };
    TH1D* histoGammaSecondaryFromXConvProb_MCPt_Unscaled[3]         = { NULL, NULL, NULL };
    TH1D* histoGammaSecondaryFromXConvProb_MCPtOrBin[3]             = { NULL, NULL, NULL };
    TH1D* histoGammaSecondaryFromXConvProb_MCPtOrBin_Unscaled[3]    = { NULL, NULL, NULL };
    TH1D* ratioGammaConvProbMCPt[3]                                 = { NULL, NULL, NULL };
    Double_t constOffsetConvProbMCPt[3]                             = { 1, 1, 1};
    
    if ( hasCocktailInput && isPCM ) {
        TF1* linear                                                   = new TF1("linear", "[0]/TMath::Power((x-[1]), [2])", 0, 50);        
        for (Int_t k = 0; k < 3; k++){
            histoGammaSecondaryFromXConvProb_MCPt[k]                = (TH1D*)fileCorrections->Get(Form("SecondaryGammaFromXFrom%sConvProb_MCPt", nameSecondaries[k].Data()));
            histoGammaSecondaryFromXConvProb_MCPtOrBin[k]           = (TH1D*)fileCorrections->Get(Form("SecondaryGammaFromXFrom%sConvProb_MCPtOrBin", nameSecondaries[k].Data()));
            histoGammaSecondaryFromXConvProb_MCPt_Unscaled[k]       = (TH1D*)histoGammaSecondaryFromXConvProb_MCPt[k]->Clone(Form("SecondaryGammaFromXFrom%sConvProb_MCPt_Unscaled", nameSecondaries[k].Data()));
            histoGammaSecondaryFromXConvProb_MCPtOrBin_Unscaled[k]  = (TH1D*)histoGammaSecondaryFromXConvProb_MCPtOrBin[k]->Clone(Form("SecondaryGammaFromXFrom%sConvProb_MCPtOrBin_Unscaled",
                                                                                                                                       nameSecondaries[k].Data()));
            ratioGammaConvProbMCPt[k]                               = (TH1D*)histoGammaSecondaryFromXConvProb_MCPt[k]->Clone(Form("RatioConvProbFrom%sToPrim_MCPt", nameSecondaries[k].Data()));
            ratioGammaConvProbMCPt[k]->Divide(histoGammaConvProb_MCPt);
            ratioGammaConvProbMCPt[k]->Fit(linear,"SMNR0E+", "", minPtFitSec[k], maxPtFitSec);
            constOffsetConvProbMCPt[k]                              = linear->GetParameter(0);
            
            for (Int_t i=histoGammaSecondaryFromXConvProb_MCPt[k]->FindBin(minPtFitSec[k]); i<histoGammaSecondaryFromXConvProb_MCPt[k]->GetNbinsX()+1; i++) {    
                histoGammaSecondaryFromXConvProb_MCPt[k]->SetBinContent(i, histoGammaConvProb_MCPt->GetBinContent(i)*linear->Eval(histoGammaSecondaryFromXConvProb_MCPt[k]->GetBinCenter(i)));
                histoGammaSecondaryFromXConvProb_MCPt[k]->SetBinError(  i, histoGammaConvProb_MCPt->GetBinError(i)*linear->Eval(histoGammaSecondaryFromXConvProb_MCPt[k]->GetBinCenter(i)));
            }
            for (Int_t i=histoGammaSecondaryFromXConvProb_MCPtOrBin[k]->FindBin(minPtFitSec[k]); i<histoGammaSecondaryFromXConvProb_MCPtOrBin[k]->GetNbinsX()+1; i++) {    
                histoGammaSecondaryFromXConvProb_MCPtOrBin[k]->SetBinContent(i, histoGammaConvProb_MCPt_OrBin->GetBinContent(i)*linear->Eval(histoGammaSecondaryFromXConvProb_MCPtOrBin[k]->GetBinCenter(i)));
                histoGammaSecondaryFromXConvProb_MCPtOrBin[k]->SetBinError(  i, histoGammaConvProb_MCPt_OrBin->GetBinError(i)*linear->Eval(histoGammaSecondaryFromXConvProb_MCPtOrBin[k]->GetBinCenter(i)));
            }
            //****************************************************************************************** 
            //************************* plot ratio conversion probabilies ****************************** 
            //****************************************************************************************** 
            TCanvas *canvasSecConvProbRatio                               = GetAndSetCanvas("canvasSecConvProbRatio");
                SetHistogramm(ratioGammaConvProbMCPt[k],"#it{p}_{T} (GeV/#it{c})","P_{conv, sec}/P_{conv, prim}",0.0,3.0);
                DrawGammaSetMarker(ratioGammaConvProbMCPt[k], 20, 1, kRed+2, kRed+2);
                ratioGammaConvProbMCPt[k]->Draw();
                DrawGammaSetMarkerTF1( linear, 9, 2, kRed-6); 
                linear->Draw("same");
                
                TLegend* legendSecConvProbRatio = GetAndSetLegend2(0.15,0.93-3*1.1*0.035, 0.4,0.93, 0.035, 1, cent, 42, 0.1);
                legendSecConvProbRatio->AddEntry(ratioGammaConvProbMCPt[k], Form("sec %s P_{conv}/prim",nameLabelSecondaries[k].Data()),"lp");
                legendSecConvProbRatio->AddEntry(linear, "linear fit","l");
                legendSecConvProbRatio->Draw();
            canvasSecConvProbRatio->SaveAs(Form("%s/%s_RatioSecConvProbToPrim%sPt_%s_%s.%s",outputDir.Data(),textPi0New.Data(),nameSecondaries[k].Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));        
            delete canvasSecConvProbRatio;
            
            //****************************************************************************************** 
            //*************** plot conversion probabilities for secondaries and fits *******************
            //****************************************************************************************** 
            TCanvas *canvasSecConvProbFit     = GetAndSetCanvas("canvasSecConvProbFit");
            canvasSecConvProbFit->SetTopMargin(0.035);
                SetHistogramm(histoGammaConvProb_MCPt,"#it{p}_{T} (GeV/#it{c})","P_{conv} in |#eta| < 0.9",0.0,1.5e-1);
                DrawGammaSetMarker(histoGammaConvProb_MCPt, 20, 1, kGray+2, kGray+2);
                histoGammaConvProb_MCPt->Draw();

                DrawGammaSetMarker(histoGammaConvProb_MCPt_OrBin, 24, 1, kBlack, kBlack);
                histoGammaConvProb_MCPt_OrBin->Draw("same");
                            
                DrawGammaSetMarker(histoGammaSecondaryFromXConvProb_MCPt_Unscaled[k], 20, 1, kBlue-8, kBlue-8);
                histoGammaSecondaryFromXConvProb_MCPt_Unscaled[k]->Draw("same");
                DrawGammaSetMarker(histoGammaSecondaryFromXConvProb_MCPtOrBin_Unscaled[k], 24, 1, kBlue+2, kBlue+2);
                histoGammaSecondaryFromXConvProb_MCPtOrBin_Unscaled[k]->Draw("same");
                DrawGammaSetMarker(histoGammaSecondaryFromXConvProb_MCPt[k], 20, 1, kRed-8, kRed-8);
                histoGammaSecondaryFromXConvProb_MCPt[k]->Draw("same");
                DrawGammaSetMarker(histoGammaSecondaryFromXConvProb_MCPtOrBin[k], 24, 1, kRed+2, kRed+2);
                histoGammaSecondaryFromXConvProb_MCPtOrBin[k]->Draw("same");
                
                TLegend* legendSecConvProbFits = GetAndSetLegend2(0.15,0.93-4*1.1*0.035, 0.7,0.93, 0.035, 2, cent, 42, 0.1);
                legendSecConvProbFits->AddEntry(histoGammaConvProb_MCPt,"primary","lp");
                legendSecConvProbFits->AddEntry(histoGammaConvProb_MCPt_OrBin,"rebin primary","lp");
                legendSecConvProbFits->AddEntry(histoGammaSecondaryFromXConvProb_MCPtOrBin_Unscaled[k], Form("from %s ",nameLabelSecondaries[k].Data()),"lp");
                legendSecConvProbFits->AddEntry(histoGammaSecondaryFromXConvProb_MCPt_Unscaled[k], Form("rebin from %s",nameLabelSecondaries[k].Data()),"lp");
                legendSecConvProbFits->AddEntry(histoGammaSecondaryFromXConvProb_MCPtOrBin[k], Form("new from  %s reco. eff",nameLabelSecondaries[k].Data()),"lp");
                legendSecConvProbFits->AddEntry(histoGammaSecondaryFromXConvProb_MCPt[k], Form("new rebin from %s ",nameLabelSecondaries[k].Data()),"lp");
                legendSecConvProbFits->Draw();

            canvasSecConvProbFit->SaveAs(Form("%s/%s_ConversionProbFits%sMCPt_%s_%s.%s",outputDir.Data(),textPi0New.Data(),nameSecondaries[k].Data(), nameRec.Data(),cutSelection.Data(),suffix.Data()));        
            delete canvasSecConvProbFit;
        }
        delete linear;
    }

    //****************************************************************************************** 
    //******************* Secondary gamma response matrices ************************************
    //****************************************************************************************** 
    TH2D* histoGammaTrueSecondaryFromX_MCPt_recPt[3]                = { NULL, NULL, NULL };
    TH2D* histoGammaTrueSecondaryFromX_MCPt_recPtOrBin[3]           = { NULL, NULL, NULL };
    // PCM
    if ( hasCocktailInput && isPCM && !isCalo ) {
        for (Int_t k = 0; k <3; k++){
            histoGammaTrueSecondaryFromX_MCPt_recPt[k]              = (TH2D*)fileCorrections->Get(Form("TrueSecondaryConvGammaFromXFrom%s_MCPt_recPt",nameSecondaries[k].Data()));
            histoGammaTrueSecondaryFromX_MCPt_recPtOrBin[k]         = (TH2D*)fileCorrections->Get(Form("TrueSecondaryConvGammaFromXFrom%s_MCPt_recPt_orBin",nameSecondaries[k].Data()));
        }    
    }
    // Calo
    if ( hasCocktailInput && isCalo && !isPCM ) {
        for (Int_t k = 0; k <3; k++){
            histoGammaTrueSecondaryFromX_MCPt_recPt[k]              = (TH2D*)fileCorrections->Get(Form("TrueSecondaryCaloGammaFromXFrom%s_MCPt_recPt",nameSecondaries[k].Data()));
            histoGammaTrueSecondaryFromX_MCPt_recPtOrBin[k]         = (TH2D*)fileCorrections->Get(Form("TrueSecondaryCaloGammaFromXFrom%s_MCPt_recPt_orBin",nameSecondaries[k].Data()));
        }    
    }
    
    //****************************************************************************************** 
    //******************* Calculate raw secondary spectra from cocktail input ******************
    //****************************************************************************************** 
    TH1D* histoGammaSecGammaFromX_Cocktail_Raw_Pt[3]                = { NULL, NULL, NULL };
    TH1D* histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin[3]           = { NULL, NULL, NULL };
    Double_t scaleCocktailDecayLength[3]                            = {1.0, 1.0, 0.5/15.34};
    if ( hasCocktailInput && (isPCM || isCalo) ) {
        cout << "calculating raw secondary spectra from cocktail" << endl;
        // K0s: clone cocktail spectra
        for (Int_t k = 0; k < 3; k++){
            histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]              = (TH1D*) histoGammaTrueSecCocktailGammaFromX_Pt[k]->Clone(Form("histoGammaTrueSecConvGammaFromXFrom%s_Cocktail_Raw_Pt", 
                                                                                                                              nameSecondaries[k].Data()));
            histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]->Sumw2();
            histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin[k]         = (TH1D*) histoGammaTrueSecCocktailGammaFromX_PtOrBin[k]->Clone(Form("histoGammaTrueSecConvGammaFromXFrom%s_Cocktail_Raw_PtOrBin", 
                                                                                                                                   nameSecondaries[k].Data()));
            histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin[k]->Sumw2();
        }
        
        if ( isPCM && !isCalo ) {
            // K0s: calculate raw yield
            if (useUnfoldingForCocktailSecondaries) {
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt[0], histoGammaSecondaryFromXConvProb_MCPt[0],
                                                                                                        histoGammaSecFromXRecoEff_MCPt[0], histoGammaTrueSecondaryFromX_MCPt_recPt[0], nEvt, kTRUE,
                                                                                                        nIterationsUnfolding);
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin[0], histoGammaSecondaryFromXConvProb_MCPtOrBin[0],
                                                                                                        histoGammaSecFromXRecoEff_MCPtOrBin[0], histoGammaTrueSecondaryFromX_MCPt_recPtOrBin[0], nEvt, kTRUE,
                                                                                                        nIterationsUnfolding);
            } else {
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt[0], histoGammaSecondaryFromXConvProb_MCPt[0],
                                                                                                        histoGammaSecFromXRecoEff_RecPt[0], histoGammaTrueSecondaryFromX_MCPt_recPt[0], nEvt, kFALSE);
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin[0], histoGammaSecondaryFromXConvProb_MCPtOrBin[0],
                                                                                                        histoGammaSecFromXRecoEff_RecPtOrBin[0], histoGammaTrueSecondaryFromX_MCPt_recPtOrBin[0], nEvt, kFALSE);
            }
            // K0l and Lambda
            for (Int_t k = 1; k < 3; k++){
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k], histoGammaSecondaryFromXConvProb_MCPt[k],
                                                                                                        histoGammaSecFromXRecoEff_RecPt[k], histoGammaTrueSecondaryFromX_MCPt_recPt[k], nEvt, kFALSE);
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin[k], histoGammaSecondaryFromXConvProb_MCPtOrBin[k],
                                                                                                        histoGammaSecFromXRecoEff_RecPtOrBin[k], histoGammaTrueSecondaryFromX_MCPt_recPtOrBin[k], nEvt,                                                                                            
                                                                                                        kFALSE);
            }
        }
        
        if ( isCalo && !isPCM ) {
            // K0s: calculate raw yield
            if (useUnfoldingForCocktailSecondaries) {
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt[0], histoGammaSecFromXRecoEff_MCPt[0],
                                                                                                        histoGammaTrueSecondaryFromX_MCPt_recPt[0], nEvt, kTRUE, nIterationsUnfolding);
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin[0], histoGammaSecFromXRecoEff_MCPtOrBin[0],
                                                                                                        histoGammaTrueSecondaryFromX_MCPt_recPtOrBin[0], nEvt, kTRUE, nIterationsUnfolding);
            } else {
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt[0], histoGammaSecFromXRecoEff_RecPt[0],
                                                                                                        histoGammaTrueSecondaryFromX_MCPt_recPt[0], nEvt, kFALSE);
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin[0], histoGammaSecFromXRecoEff_RecPtOrBin[0],
                                                                                                        histoGammaTrueSecondaryFromX_MCPt_recPtOrBin[0], nEvt, kFALSE);
            }
            // K0l and Lambda
            for (Int_t k = 1; k < 3; k++){            
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k], histoGammaSecFromXRecoEff_RecPt[k],
                                                                                                        histoGammaTrueSecondaryFromX_MCPt_recPt[k], nEvt, kFALSE);
                hasCocktailInput                                        = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin[k], histoGammaSecFromXRecoEff_RecPtOrBin[k], 
                                                                                                        histoGammaTrueSecondaryFromX_MCPt_recPtOrBin[k], nEvt, kFALSE);
            }
        }
        if (hasCocktailInput){
            for (Int_t k = 0; k<3; k++){
                histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]->Scale(scaleCocktailDecayLength[k]);
                histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin[k]->Scale(scaleCocktailDecayLength[k]);
            }    
        }        
    }

    //****************************************************************************************** 
    //******************* Calculate secondary fractions from cocktail input ********************
    //****************************************************************************************** 
    TH1D* histoFracAllGammaToSecFromX_Cocktail_Pt[4]                    = { NULL, NULL, NULL, NULL };
    TH1D* histoFracAllGammaToSecFromX_Cocktail_PtOrBin[3]               = { NULL, NULL, NULL };
    if ( hasCocktailInput && isPCM && !isCalo ) {
        cout << "calculating secondary fractions from cocktail" << endl;
        for (Int_t k = 0; k < 3; k++){
            histoFracAllGammaToSecFromX_Cocktail_Pt[k]                  = (TH1D*)histoESDConvGammaPt->Clone(Form("FracAllGammaToSecFromXFrom%s", nameSecondaries[k].Data()));
            histoFracAllGammaToSecFromX_Cocktail_Pt[k]->Divide(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k],histoFracAllGammaToSecFromX_Cocktail_Pt[k],1,1,"B");
            histoFracAllGammaToSecFromX_Cocktail_PtOrBin[k]             = (TH1D*)histoESDConvGammaPt_OriginalBin->Clone(Form("FracAllGammaToSecFromXFrom%sOriginalBinning", nameSecondaries[k].Data()));
            histoFracAllGammaToSecFromX_Cocktail_PtOrBin[k]->Divide(histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin[k],histoFracAllGammaToSecFromX_Cocktail_PtOrBin[k],1,1,"B");
        }
        
        histoFracAllGammaToSecFromX_Cocktail_Pt[3]                      = (TH1D*)histoESDConvGammaPt->Clone("FracAllGammaToSecRest");
        histoFracAllGammaToSecFromX_Cocktail_Pt[3]->Divide(histoGammaTrueSecCocktailGammaFromX_Pt[3],histoFracAllGammaToSecFromX_Cocktail_Pt[3],1,1,"B");
    }
    if ( hasCocktailInput && isCalo && !isPCM ) {
        cout << "calculating secondary fractions from cocktail" << endl;
        for (Int_t k = 0; k < 3; k++){
            histoFracAllGammaToSecFromX_Cocktail_Pt[k]                  = (TH1D*)histoESDCaloGammaPt->Clone(Form("FracAllGammaToSecFromXFrom%s", nameSecondaries[k].Data()));
            histoFracAllGammaToSecFromX_Cocktail_Pt[k]->Divide(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k],histoFracAllGammaToSecFromX_Cocktail_Pt[k],1,1,"B");
            histoFracAllGammaToSecFromX_Cocktail_PtOrBin[k]             = (TH1D*)histoESDCaloGammaPt_OriginalBin->Clone(Form("FracAllGammaToSecFromXFrom%sOriginalBinning", nameSecondaries[k].Data()));
            histoFracAllGammaToSecFromX_Cocktail_PtOrBin[k]->Divide(histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin[k],histoFracAllGammaToSecFromX_Cocktail_PtOrBin[k],1,1,"B");
        }
        
        histoFracAllGammaToSecFromX_Cocktail_Pt[3]                      = (TH1D*)histoESDCaloGammaPt->Clone("FracAllGammaToSecRest");
        histoFracAllGammaToSecFromX_Cocktail_Pt[3]->Divide(histoGammaTrueSecCocktailGammaFromX_Pt[3],histoFracAllGammaToSecFromX_Cocktail_Pt[3],1,1,"B");
    }

    //****************************************************************************************** 
    //******************* Pileup correction factors ********************************************
    //****************************************************************************************** 
    TFile*  doPileUpCorr                                                = NULL;
    if (kDoPileup){
        doPileUpCorr                                                    = new  TFile(Form("%s/%s/%s_%s_GammaConvV1DCAHistogramms%s_%s.root", cutSelection.Data(), energy.Data(), nameMeson.Data(),nameRec.Data(), optionPeriod.Data(), cutSelection.Data()));
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
                                                                                                                    histoRatioWithWithoutPileUpFitResult->GetParams()) / 
                                                                                                                    histoPileUpCorrectionFactorOrBin->GetBinWidth(i);
                binError                                                = histoRatioWithWithoutPileUpFit->IntegralError(    histoPileUpCorrectionFactorOrBin->GetXaxis()->GetBinLowEdge(i),
                                                                                                                            histoPileUpCorrectionFactorOrBin->GetXaxis()->GetBinUpEdge(i),
                                                                                                                            histoRatioWithWithoutPileUpFitResult->GetParams(),
                                                                                                                            histoRatioWithWithoutPileUpFitResult->GetCovarianceMatrix().GetMatrixArray()) /
                                                                                                                            histoPileUpCorrectionFactorOrBin->GetBinWidth(i);
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

    //****************************************************************************************** 
    //******************* MC pileup histograms *************************************************
    //****************************************************************************************** 
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

    //****************************************************************************************** 
    //******************* Proper scaling of background *****************************************
    //****************************************************************************************** 
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

    //****************************************************************************************** 
    //******************* Calculate secondary spectra from data ********************************
    //****************************************************************************************** 
    TH1D *histoSecondaryGammaFromXSpecPt[4]                             = { NULL, NULL, NULL, NULL };
    if (isPCM && !isCalo) {
        for (Int_t k = 0; k < 4; k++){
            histoSecondaryGammaFromXSpecPt[k]                           = (TH1D*)histoESDConvGammaPt->Clone(Form("SecondaryGammaSpecFromXFrom%sPt",nameSecondaries[k].Data()));
            // overwrite spectra loaded from pure MC if running data
            if(!isRunMC) histoGammaTrueSecConvGammaFromX_Pt[k]          = (TH1D*)histoESDConvGammaPt->Clone(Form("SecondaryMCGammaSpecFromXFrom%sPt",nameSecondaries[k].Data()));
        }
    }
    if (isCalo && !isPCM) {
        for (Int_t k = 0; k < 4; k++){
            histoSecondaryGammaFromXSpecPt[k]                           = (TH1D*)histoESDCaloGammaPt->Clone(Form("SecondaryGammaSpecFromXFrom%sPt",nameSecondaries[k].Data()));
            // overwrite spectra loaded from pure MC if running data
            if (!isRunMC)histoGammaTrueSecCaloGammaFromX_Pt[k]          = (TH1D*)histoESDCaloGammaPt->Clone(Form("SecondaryMCGammaSpecFromXFrom%sPt",nameSecondaries[k].Data()));
        }
    }
    for (Int_t k = 0; k < 4; k++){
        if (histoFracAllGammaToSecFromX_Pt[k]) {
            histoSecondaryGammaFromXSpecPt[k]->Multiply(histoFracAllGammaToSecFromX_Pt[k]);
            histoSecondaryGammaFromXSpecPt[k]->Scale(1./nEvt);
            histoSecondaryGammaFromXSpecPt[k]->Scale(scaleFactorsSec[k]);
            // scale reconstructed gamma with secondaries
            if (!isRunMC && isPCM && !isCalo) {
                histoGammaTrueSecConvGammaFromX_Pt[k]->Multiply(histoFracAllGammaToSecFromX_Pt[k]);
            } else if (!isRunMC && isCalo && !isPCM) {
                histoGammaTrueSecCaloGammaFromX_Pt[k]->Multiply(histoFracAllGammaToSecFromX_Pt[k]);
            }    
        } else {
            histoSecondaryGammaFromXSpecPt[k]                           = NULL;
            if (!isRunMC && isPCM && !isCalo) {
                histoGammaTrueSecConvGammaFromX_Pt[k]                   = NULL;
            } else if (!isRunMC && isCalo && !isPCM) {
                histoGammaTrueSecCaloGammaFromX_Pt[k]                   = NULL;
            }    
        }
    }
    
    //****************************************************************************************** 
    //******************* Scale MC/cocktail secondary spectra **********************************
    //****************************************************************************************** 
    if (isPCM && !isCalo) {
        for (Int_t k = 0; k<4; k++){
            if(histoGammaTrueSecConvGammaFromX_Pt[k]) histoGammaTrueSecConvGammaFromX_Pt[k]->Scale(1./nEvtMC);
        }
    }
    if (isCalo && !isPCM) {
        for (Int_t k = 0; k<4; k++){
            if (histoGammaTrueSecCaloGammaFromX_Pt[k]) histoGammaTrueSecCaloGammaFromX_Pt[k]->Scale(1./nEvtMC);       
        }
    }    
    if (hasCocktailInput) {
        for (Int_t k = 0; k<3; k++){
            histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]->Scale(1./nEvt);
        }
        histoGammaTrueSecCocktailGammaFromX_Pt[3]->Scale(1./nEvtMC);
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
    if( doPileUpCorr && isPCM && !isRunMC){
        TCanvas *canvasPileUpCorrFactor     = GetAndSetCanvas("canvasPileUpCorrFactor");

            DrawGammaSetMarker(histoPileUpCorrectionFactor, 20, 3, 1, 1);
            if (histoMCrecGamma_PileUp_Pt)DrawGammaSetMarker(histoPileUpCorrectionFactorMC, 24, 3, 2, 2);

            SetHistogramm(histoPileUpCorrectionFactor,"#it{p}_{T} (GeV/#it{c})","Correction Factor (%)",0.84,1.02);
            if (histoMCrecGamma_PileUp_Pt)SetHistogramm(histoPileUpCorrectionFactorMC,"#it{p}_{T} (GeV/#it{c})","Correction Factor (%)",0.84,1.02);

            histoPileUpCorrectionFactor->DrawCopy("");
            if (histoMCrecGamma_PileUp_Pt)histoPileUpCorrectionFactorMC->Draw("same");

            TLegend* legendPileUpCorrFactor = GetAndSetLegend(0.55,0.2,2.2,1,cent);
            legendPileUpCorrFactor->AddEntry(histoPileUpCorrectionFactor,"Correction Factor Data","lp");
            if (histoMCrecGamma_PileUp_Pt)legendPileUpCorrFactor->AddEntry(histoPileUpCorrectionFactorMC,"Correction Factor MC","lp");
            legendPileUpCorrFactor->Draw();
        
        canvasPileUpCorrFactor->SaveAs(Form("%s/%s_PileUpCorrFactor_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasPileUpCorrFactor;
        
        Bool_t doPileUpCorrSimplePlot       = kTRUE;
        if (doPileUpCorrSimplePlot && !isRunMC) {
            TCanvas *canvasPileUpCorrFactor2= GetAndSetCanvas("canvasPileUpCorrFactor2");
            
            SetHistogramm(histoPileUpCorrectionFactorOrBin,"#it{p}_{T} (GeV/#it{c})","Correction Factor (%)",0.90,1.02);
            DrawGammaSetMarker(histoPileUpCorrectionFactorOrBin, 24, 1.5, kBlack, kBlack);
            histoPileUpCorrectionFactorOrBin->DrawCopy("");
            DrawGammaLines(0., histoPileUpCorrectionFactorOrBin->GetXaxis()->GetBinUpEdge(histoPileUpCorrectionFactorOrBin->GetNbinsX()),1.0, 1.0, 1, kGray+2, 2);
            histoPileUpCorrectionFactorOrBin->DrawCopy("same");

            PutProcessLabelAndEnergyOnPlot( 0.15, 0.25, 0.035, cent, detectionProcess, "", 42, 0.03);
            
            canvasPileUpCorrFactor2->SaveAs(Form("%s/%s_%s_PileUpCorrFactor_%s.%s",outputDir.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
            delete canvasPileUpCorrFactor2;
        }
        
        if (histoPileUpCorrectionFactorOrBin && !isRunMC) {
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
            
            PutProcessLabelAndEnergyOnPlot( 0.15, 0.25, 0.035, cent, detectionProcess, "", 42, 0.03);
            
            canvasPileUpCorrFactor3->SaveAs(Form("%s/%s_%s_PileUpCorrFactorRebinned_%s.%s",outputDir.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
            delete canvasPileUpCorrFactor3;
        }
    }

    //**********************************************************************************
    //******************** Background Plot *********************************************
    //**********************************************************************************
    TCanvas* canvasBackground                      = new TCanvas("canvasBackground","",200,10,1000*1.25,1100*1.25);  // gives the page size
    DrawGammaCanvasSettings( canvasBackground, 0.1, 0.02, 0.02, 0.09);
    canvasBackground->SetLogy();
    
        TH1D* histoGammaRawSpectrum_PtMC    = NULL;
        TH1D* histoGammaRawSpectrum_Pt      = NULL;
        TLegend* legendBackground           = GetAndSetLegend(0.6,0.75,2,1);
        if (isPCM && !isCalo) {
            histoGammaRawSpectrum_PtMC      = (TH1D*) histoMCrecGamma_Pt->Clone("histoGammaRawSpectrum_PtMC");
            histoGammaRawSpectrum_PtMC->Scale(1./nEvtMC);
            histoGammaRawSpectrum_Pt        = (TH1D*) histoESDConvGammaPt->Clone("histoGammaRawSpectrum_Pt");
            histoGammaRawSpectrum_Pt->Scale(1./nEvt);
            
            DrawGammaSetMarker(histoGammaRawSpectrum_Pt, 20, 1.0, 1, 1);
            DrawGammaSetMarker(histoMCrecBackground_Pt, 20, 1.0, 2, 2);
            
            SetHistogramm(histoMCrecBackground_Pt,"#it{p}_{T} (GeV/#it{c})",Form("Background for #gamma in |#eta| < %g",eta),-99,-99,1.0,1.7);
            SetHistogramm(histoGammaRawSpectrum_Pt,"#it{p}_{T} (GeV/#it{c})","Raw yield",1e-8,5,1.0,1.2);
            
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
            
            SetHistogramm(histoMCrecBackground_Pt,"#it{p}_{T} (GeV/#it{c})",Form("Background for #gamma in |#eta| < %g",etaCalo),-99,-99,1.0,1.7);
            SetHistogramm(histoGammaRawSpectrum_Pt,"#it{p}_{T} (GeV/#it{c})","Raw yield",1e-10,1e-2,1.0,1.2);
            
            histoGammaRawSpectrum_Pt->DrawCopy("");
            histoMCrecBackground_Pt->Draw("same");
        }

        legendBackground->AddEntry(histoGammaRawSpectrum_Pt,"raw #gamma spectrum","pl");
        legendBackground->AddEntry(histoMCrecBackground_Pt,"#gamma background","pl");
        legendBackground->Draw();
    
        PutProcessLabelAndEnergyOnPlot( 0.935, 0.95, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,31);

    canvasBackground->SaveAs(Form("%s/%s_Background_%s_%s.%s",outputDir.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasBackground;

    //**********************************************************************************
    //******************** Secondary Spectra Plot **************************************
    //**********************************************************************************
    TCanvas* canvasSecSpec                      = new TCanvas("canvasSecSpec","",200,10,1000*1.25,1100*1.25);  // gives the page size
    DrawGammaCanvasSettings( canvasSecSpec, 0.15, 0.02, 0.02, 0.09);
    canvasSecSpec->SetLogy();
    
        TLegend* legendSecSpec                  = GetAndSetLegend2(0.35, 0.935-0.035*1.1*4, 0.95, 0.935, 0.035, 3, "", 42, 0.1);

        // secondary gamma spectra plotting as would they would be subtracted using MC fractions
        for (Int_t k = 0; k < 4; k++){
            legendSecSpec->AddEntry((TObject*)0, Form("sec #gamma from %s:", nameLabelSecondaries[k].Data()),"");
            // plotting MC unscaled secondaries
            if (isPCM && !isCalo) {
                if (histoGammaTrueSecConvGammaFromX_Pt[k]){
                    SetHistogramm(histoGammaTrueSecConvGammaFromX_Pt[k],"#it{p}_{T} (GeV/#it{c})","#frac{1}{#it{N}_{ev.}} #frac{d#it{N}}{d#it{p}_{T}} (#it{c}/GeV)",1e-10,1e-1,1.0,1.6);
                    DrawGammaSetMarker(histoGammaTrueSecConvGammaFromX_Pt[k], markerStyleSecWithToy[k], markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                    histoGammaTrueSecConvGammaFromX_Pt[k]->Draw("same");
                    legendSecSpec->AddEntry(histoGammaTrueSecConvGammaFromX_Pt[k], "MC     ","pl");
                } else {
                    legendSecSpec->AddEntry((TObject*)0, "","");
                }    
            } else if (isCalo && !isPCM) {
                if (histoGammaTrueSecCaloGammaFromX_Pt[k]){
                    SetHistogramm(histoGammaTrueSecCaloGammaFromX_Pt[k],"#it{p}_{T} (GeV/#it{c})","#frac{1}{#it{N}_{ev.}} #frac{d#it{N}}{d#it{p}_{T}} (#it{c}/GeV)",1e-10,1e-1,1.0,1.6);
                    DrawGammaSetMarker(histoGammaTrueSecCaloGammaFromX_Pt[k], markerStyleSecWithToy[k], markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                    histoGammaTrueSecCaloGammaFromX_Pt[k]->Draw("same");
                    legendSecSpec->AddEntry(histoGammaTrueSecCaloGammaFromX_Pt[k], "MC    ","pl");
                } else {
                    legendSecSpec->AddEntry((TObject*)0, "","");
                }
            }
            // plotting cocktail on top of the MC unscaled secondaries
            if (hasCocktailInput) {
                if (histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]){
                    if (isPCM && !isCalo) {
                        SetHistogramm(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k],"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
                    } else if (isCalo && !isPCM) {
                        SetHistogramm(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k],"#it{p}_{T} (GeV/#it{c})","Secondary #gamma");
                    }    
                    DrawGammaSetMarker(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k], markerStyleSec[k] , markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                    histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]->DrawCopy("same");
                    legendSecSpec->AddEntry(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k],"cocktail","pl");
                } else {
                    legendSecSpec->AddEntry((TObject*)0, "","");
                }    
            } else {
                if (histoSecondaryGammaFromXSpecPt[k]){
                    if (isPCM && !isCalo) {
                        SetHistogramm(histoSecondaryGammaFromXSpecPt[k], "#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
                    } else if (isCalo && !isPCM) {
                        SetHistogramm(histoSecondaryGammaFromXSpecPt[k], "#it{p}_{T} (GeV/#it{c})","Secondary #gamma");
                    }
                    DrawGammaSetMarker(histoSecondaryGammaFromXSpecPt[k], markerStyleSec[k], markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                    histoSecondaryGammaFromXSpecPt[k]->Draw("same");
                    legendSecSpec->AddEntry(histoSecondaryGammaFromXSpecPt[k],"data scaled","pl");
                } else {
                    legendSecSpec->AddEntry((TObject*)0, "","");
                }
            }
        }
        legendSecSpec->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.95, 0.935-0.035*1.1*4, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,31);
        
    canvasSecSpec->SaveAs(Form("%s/%s_%s_SecondarySpectra_%s.%s", outputDir.Data(), textPi0New.Data(), nameRec.Data(), cutSelection.Data(), suffix.Data()));
    delete canvasSecSpec;

    //**********************************************************************************
    //******************** Secondary Fractions Plot ************************************
    //**********************************************************************************
    TCanvas* canvasSecFrac = new TCanvas("canvasSecFrac","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasSecFrac, 0.072, 0.02, 0.035, 0.09);
        Int_t nColumnsSec                   = 1;
        Double_t startXSecLegend            = 0.65;
        Int_t nRowsSec                      = 4;
        Double_t marginSecLegend            = 0.15;
        if (hasCocktailInput){
            nColumnsSec                     = 3;
            startXSecLegend                 = 0.55;
            marginSecLegend                 = 0.15 ;
        }
        TLegend* legendSecFrac              = GetAndSetLegend2(startXSecLegend, 0.935-0.035*1.1*nRowsSec, 0.95, 0.935, 0.035, nColumnsSec, "", 42, marginSecLegend);
        if (!hasCocktailInput) {
            for (Int_t k = 0; k < 4; k++){
                if (histoFracAllGammaToSecFromX_Pt[k]){ 
                    SetHistogramm(histoFracAllGammaToSecFromX_Pt[k],"#it{p}_{T} (GeV/#it{c})","C_{sec, #gamma}", -99., -99., 1.0, 0.9);
                    DrawGammaSetMarker(histoFracAllGammaToSecFromX_Pt[k], markerStyleSec[k], markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                    histoFracAllGammaToSecFromX_Pt[k]->Draw("same");
                    legendSecFrac->AddEntry(histoFracAllGammaToSecFromX_Pt[k],Form("#gamma from %s", nameLabelSecondaries[k].Data()),"pl");
                }
            }
        } else {
            for (Int_t k = 0; k < 4; k++){
                if (histoFracAllGammaToSecFromX_Cocktail_Pt[k]){ 
                    SetHistogramm(histoFracAllGammaToSecFromX_Cocktail_Pt[k],"#it{p}_{T} (GeV/#it{c})","C_{sec, #gamma}", -99., -99., 1.0, 0.9);
                    SetHistogramm(histoFracAllGammaToSecFromX_Pt[k],"#it{p}_{T} (GeV/#it{c})","C_{sec, #gamma}", -99., -99., 1.0, 0.9);
                    DrawGammaSetMarker(histoFracAllGammaToSecFromX_Cocktail_Pt[k], markerStyleSec[k], markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                    DrawGammaSetMarker(histoFracAllGammaToSecFromX_Pt[k], markerStyleSecWithToy[k], markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                    if (k == 0){
                        histoFracAllGammaToSecFromX_Pt[k]->Draw();
                        histoFracAllGammaToSecFromX_Cocktail_Pt[k]->Draw("same");
                    } else {
                        histoFracAllGammaToSecFromX_Pt[k]->Draw("same");
                        histoFracAllGammaToSecFromX_Cocktail_Pt[k]->Draw("same");
                    }   
                    legendSecFrac->AddEntry((TObject*)0,Form("#gamma from %s", nameLabelSecondaries[k].Data()),"");   
                    legendSecFrac->AddEntry(histoFracAllGammaToSecFromX_Cocktail_Pt[k],"cocktail","pl");   
                    legendSecFrac->AddEntry(histoFracAllGammaToSecFromX_Pt[k],"MC","pl");   
                }
            }
        }
        legendSecFrac->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.94, 0.94-0.035*1.1*nRowsSec, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,31);
        
    canvasSecFrac->SaveAs(Form("%s/%s_%s_SecondaryGammaFraction_%s.%s", outputDir.Data(), textPi0New.Data(), nameRec.Data(), cutSelection.Data(), suffix.Data()));
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
            
            if (doPileUpCorr)   legendPurity    = GetAndSetLegend(0.16,0.15,3);
            else                legendPurity    = GetAndSetLegend(0.16,0.15,2);
            legendPurity->AddEntry(histoGammaPurity_Pt,"purity");
            legendPurity->AddEntry(histoGammaTruePurity_Pt,"rescaled purity, secondaries removed");
            if (doPileUpCorr) legendPurity->AddEntry(histoGammaTruePurity_PileUp_Pt,"rescaled purity, secondaries removed, pileup");
            legendPurity->Draw();

            PutProcessLabelAndEnergyOnPlot( 0.18, 0.4, 0.035, cent, detectionProcess, "", 42, 0.03);
        
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
        
        PutProcessLabelAndEnergyOnPlot( 0.18, 0.3, 0.035, cent, detectionProcess, "", 42, 0.03);
        
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

            PutProcessLabelAndEnergyOnPlot( 0.95, 0.23, 0.035, cent, detectionProcess, "", 42, 0.03,"", 1, 1.25,31);

        canvasConvProb->SaveAs(Form("%s/%s_ConversionProb_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasConvProb;

        // secondary reco eff plot (only used for cocktail sec corr)
        if (hasCocktailInput) {
            TCanvas *canvasConvProbSec  = GetAndSetCanvas("canvasConvProbSec");
            canvasConvProbSec->SetTopMargin(0.035);
            TLegend* legendSecConvProb       = GetAndSetLegend2(0.25, 0.935-0.035*1.4*2, 0.65, 0.935, 0.035, 2, "", 42, 0.12);
        
            DrawGammaSetMarker(histoGammaConvProb_MCPt, 20, 2.0, 1, 1);
            SetHistogramm(histoGammaConvProb_MCPt, "#it{p}_{T} (GeV/#it{c})",Form("#it{P}_{conv} in |#eta| < %g",eta), 0.0, 0.15);
            histoGammaConvProb_MCPt->Draw();
            legendSecConvProb->AddEntry(histoGammaConvProb_MCPt,"Prim. #gamma","p");
                
            for (Int_t k = 0; k < 3; k++){
                SetHistogramm(histoGammaSecondaryFromXConvProb_MCPt[k],"#it{p}_{T} (GeV/#it{c})",Form("#it{P}_{conv} in |#eta| < %g",eta), 0., 1.0);
                DrawGammaSetMarker(histoGammaSecondaryFromXConvProb_MCPt[k], markerStyleSecWithToy[k], markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                histoGammaSecondaryFromXConvProb_MCPt[k]->Draw("same");
                legendSecConvProb->AddEntry(histoGammaSecondaryFromXConvProb_MCPt[k],Form("Sec. #gamma from %s",nameLabelSecondaries[k].Data()),"p");
            }
        
            PutProcessLabelAndEnergyOnPlot( 0.95, 0.95, 0.035, cent, detectionProcess,"", 42, 0.03,"", 1, 1.25,31);
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
        
            PutProcessLabelAndEnergyOnPlot( 0.95, 0.95, 0.035, cent, detectionProcess,"", 42, 0.03, "", 1, 1.25,31);
        
            canvasRecoEff->SaveAs(Form("%s/%s_ReconstructionEff_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        }
        if ( isCalo && !isPCM ){
            DrawGammaSetMarker(histoGammaPrimaryRecoEff_MCPt, 20, 1.0, 1, 1);
            SetHistogramm(histoGammaPrimaryRecoEff_MCPt,"#it{p}_{T,MC} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
            histoGammaPrimaryRecoEff_MCPt->Draw();
        
            PutProcessLabelAndEnergyOnPlot( 0.95, 0.95, 0.035, cent, detectionProcess,"", 42, 0.03, "", 1, 1.25,31);
        
            canvasRecoEff->SaveAs(Form("%s/%s_ReconstructionEff_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        }
        if ( isPCM && isCalo ){
            DrawGammaSetMarker(histoGammaCaloPrimaryRecoEff_MCPt, 20, 1.0, 1, 1);
            SetHistogramm(histoGammaCaloPrimaryRecoEff_MCPt,"#it{p}_{T,MC} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
            histoGammaCaloPrimaryRecoEff_MCPt->Draw();
        
            PutProcessLabelAndEnergyOnPlot( 0.95, 0.95, 0.035, cent, detectionProcess,"", 42, 0.03, "", 1, 1.25,31);
    
            canvasRecoEff->SaveAs(Form("%s/%s_ReconstructionEffCalo_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        
        }    
    delete canvasRecoEff;

    // secondary reco eff plot (only used for cocktail sec corr)
    if (hasCocktailInput) {
        TCanvas *canvasRecoEffSec       = GetAndSetCanvas("canvasRecoEffSec");
        TLegend* legendSecRecoEff       = GetAndSetLegend2(0.25, 0.935-0.035*1.4*2, 0.65, 0.935, 0.035, 2, "", 42, 0.12);
        
        DrawGammaSetMarker(histoGammaPrimaryRecoEff_Pt, 20, 1.0, 1, 1);
        if ( isPCM && !isCalo ){
            SetHistogramm(histoGammaPrimaryRecoEff_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);
        } else if ( isCalo && !isPCM ){
            SetHistogramm(histoGammaPrimaryRecoEff_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
        } else {
            SetHistogramm(histoGammaPrimaryRecoEff_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);    
        }
        legendSecRecoEff->AddEntry(histoGammaPrimaryRecoEff_Pt,"Prim. #gamma","p");
        
        histoGammaPrimaryRecoEff_Pt->Draw();
        
        if ( isPCM || isCalo ){
            for (Int_t k = 0; k < 3; k++){
                if (isPCM && !isCalo) SetHistogramm(histoGammaSecFromXRecoEff_RecPt[k],"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);
                if (isCalo && !isPCM) SetHistogramm(histoGammaSecFromXRecoEff_RecPt[k],"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
                DrawGammaSetMarker(histoGammaSecFromXRecoEff_RecPt[k], markerStyleSecWithToy[k], markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                histoGammaSecFromXRecoEff_RecPt[k]->Draw("same");
                legendSecRecoEff->AddEntry(histoGammaSecFromXRecoEff_RecPt[k],Form("Sec. #gamma from %s",nameLabelSecondaries[k].Data()),"p");
            }
        
        }
        PutProcessLabelAndEnergyOnPlot( 0.95, 0.95, 0.035, cent, detectionProcess, "", 42, 0.03,"", 1, 1.25,31);
        legendSecRecoEff->Draw();
        
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
            if(doPileUpCorr)DrawGammaSetMarker(histoGammaPrimaryRecoEff_PileUp_Pt, 25, 1.0,807, 807);

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
        
            PutProcessLabelAndEnergyOnPlot( 0.935, 0.95, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,31);
        }

        TH1D* histoGammaResolCorrEff_Pt = NULL;
        if ( (isPCM && !isCalo) || (isCalo && !isPCM) ) {
            histoGammaResolCorrEff_Pt = (TH1D*)histoGammaPrimaryRecoEff_Pt->Clone("histoGammaResolCorrEff_Pt");
            histoGammaResolCorrEff_Pt->Divide(histoGammaPrimaryRecoEff_MCPt);

            padBinCompRecoEffRatio->cd();
            TH1D* histoGammaPrimaryRecoEff_PileUp_PtRatio       = NULL;
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
        
            PutProcessLabelAndEnergyOnPlot( 0.935, 0.95, 0.035, cent, detectionProcess2, "", 42, 0.03,"",1,1.25,31);
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
            PutProcessLabelAndEnergyOnPlot( 0.15, 0.95, 0.035, cent, detectionProcess, "", 42, 0.03);
            canvasResponseMatrix->SaveAs(Form("%s/%s_ResponseMatrix_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        } 
        
        if (isCalo){
            histoGammaTruePrimaryCalo_recPt_MCPt->Draw("colz");
            PutProcessLabelAndEnergyOnPlot( 0.15, 0.95, 0.035, cent,  detectionProcess2, "", 42, 0.03);
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
            SetStyleHistoTH2ForGraphs(  histoGammaTrueSecondaryFromX_MCPt_recPtOrBin[0], "Reconstructed #it{p}_{T} (GeV/#it{c})","MC #it{p}_{T} (GeV/#it{c})", 0.035, 0.04,
                                      0.035, 0.04, 0.9, 1.0, 510, 510);
            histoGammaTrueSecondaryFromX_MCPt_recPtOrBin[0]->GetYaxis()->SetRangeUser(0,25);
            histoGammaTrueSecondaryFromX_MCPt_recPtOrBin[0]->GetXaxis()->SetRangeUser(0,25);
            histoGammaTrueSecondaryFromX_MCPt_recPtOrBin[0]->Draw("colz");
            PutProcessLabelAndEnergyOnPlot( 0.15, 0.95, 0.035, cent, detectionProcess,"", 42, 0.03);
        }

        canvasResponseMatrixSec->SaveAs(Form("%s/%s_ResponseMatrixSecK0s_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasResponseMatrixSec;
    }

    //**********************************************************************************
    //************************ Unfolding of inclusive gamma spectrum *******************
    //**********************************************************************************
    TH1D* histoGammaCorrUnfoldReso_Pt                           = NULL;
    TH1D* histoGammaCorrUnfoldReso_BinByBin_Pt                  = NULL;
    TH1D* histoGammaResolCorrUnfold_Pt                          = NULL;
    TH1D* histoGammaCorrUnfoldReso_PtNotCorrected               = NULL;
    TH1D* histoGammaResolCorrUnfold_BinByBin_Pt                 = NULL;
    TH1D* histoSecondaryGammaFromXSpecPtOrBin[4]                = { NULL, NULL, NULL, NULL };
    // do the same with MC rec gammas as a sanity check
    TH1D* histoMCrecGammaCorr_Pt                                = NULL;

    if (isPCM && !isCalo){
        // Correct inclusive photon spectrum for pileup contribution
        if(doPileUpCorr) histoESDConvGammaPt_OriginalBin->Multiply(histoPileUpCorrectionFactorOrBin);
        
        // correct measured gamma spectra for secondaries and purity
        if(!hasCocktailInput){
            // determine secondary contribution in orginal binning
            histoSecondaryGammaFromXSpecPtOrBin[3]              = (TH1D*)histoESDConvGammaPt_OriginalBin->Clone("SecondaryGammaSpecPt");
            histoSecondaryGammaFromXSpecPtOrBin[3]->Multiply(histoFracAllGammaToSecFromX_OrBin_Pt[3]);
            histoSecondaryGammaFromXSpecPtOrBin[0]              = (TH1D*)histoESDConvGammaPt_OriginalBin->Clone("SecondaryGammaSpecFromXFromK0sPt");
            histoSecondaryGammaFromXSpecPtOrBin[0]->Multiply(histoFracAllGammaToSecFromX_OrBin_Pt[0]);
            histoSecondaryGammaFromXSpecPtOrBin[0]->Scale(doubleAddFactorK0s);

            //subtract secondary contribution from inclusive spectrum
            CorrectGammaSecAndPurity(   histoESDConvGammaPt_OriginalBin, 
                                        histoSecondaryGammaFromXSpecPtOrBin[3], 
                                        histoSecondaryGammaFromXSpecPtOrBin[0], 
                                        histoGammaTruePurity_OriginalBin_Pt 
                                    );
            CorrectGammaSecAndPurity(   histoMCrecGamma_OriginalBin_Pt, 
                                        histoSecondaryGammaFromXSpecPtOrBin[3],
                                        histoSecondaryGammaFromXSpecPtOrBin[0],
                                        histoGammaTruePurity_OriginalBin_Pt 
                                    );
        } else {
            //subtract secondary contribution from cocktail from inclusive spectrum
            cout << "will use cocktail for secondary correction" << endl;
            CorrectGammaSecAndPurityCocktail(   histoESDConvGammaPt_OriginalBin,
                                                histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin,
                                                histoGammaTrueSecCocktailGammaFromX_PtOrBin[3],
                                                histoGammaTruePurity_OriginalBin_Pt 
                                            );
            CorrectGammaSecAndPurityCocktail(   histoMCrecGamma_OriginalBin_Pt,
                                                histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin,
                                                histoGammaTrueSecCocktailGammaFromX_PtOrBin[3],
                                                histoGammaTruePurity_OriginalBin_Pt 
                                            );
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
        CorrectGammaUnfoldResol( histoGammaCorrUnfoldReso_Pt, 
                                 histoGammaConvProb_MCPt, 
                                 histoGammaPrimaryRecoEff_MCPt, 
                                 deltaEta, scaling, nEvt);
        CorrectGammaUnfoldResol( histoMCrecGammaCorr_Pt, 
                                 histoGammaConvProb_MCPt, 
                                 histoGammaPrimaryRecoEff_MCPt, 
                                 deltaEta, scaling, nEvt);
        CorrectGammaUnfoldResol( histoGammaCorrUnfoldReso_BinByBin_Pt,
                                 histoGammaConvProb_MCPt, 
                                 histoGammaPrimaryRecoEff_MCPt, 
                                 deltaEta, scaling, nEvt);
        //CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_SvD_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
        //CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_TUnfold_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
        
        SetHistogramm(histoGammaCorrUnfoldReso_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", -99, -99, 1.0, 1.7);
        SetHistogramm(histoGammaCorrUnfoldReso_BinByBin_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", -99, -99, 1.0, 1.7);
        SetHistogramm(histoMCrecGammaCorr_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", -99, -99, 1.0, 1.7);
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
            histoSecondaryGammaFromXSpecPtOrBin[3]              = (TH1D*)histoESDCaloGammaPt_OriginalBin->Clone("SecondaryGammaSpecPt");
            histoSecondaryGammaFromXSpecPtOrBin[0]              = (TH1D*)histoESDCaloGammaPt_OriginalBin->Clone("SecondaryGammaSpecFromXFromK0sPt");
            histoSecondaryGammaFromXSpecPtOrBin[3]->Multiply(histoFracAllGammaToSecFromX_OrBin_Pt[3]);
            histoSecondaryGammaFromXSpecPtOrBin[0]->Multiply(histoFracAllGammaToSecFromX_OrBin_Pt[0]);
            histoSecondaryGammaFromXSpecPtOrBin[0]->Scale(doubleAddFactorK0s);
            
            //subtract secondary contribution from inclusive spectrum
            CorrectGammaSecAndPurity(   histoESDCaloGammaPt_OriginalBin, 
                                        histoSecondaryGammaFromXSpecPtOrBin[3], 
                                        histoSecondaryGammaFromXSpecPtOrBin[0],
                                        histoGammaTruePurity_OriginalBin_Pt
                                    );
            CorrectGammaSecAndPurity(   histoMCrecGammaCalo_OriginalBin_Pt, 
                                        histoSecondaryGammaFromXSpecPtOrBin[3],
                                        histoSecondaryGammaFromXSpecPtOrBin[0], 
                                        histoGammaTruePurity_OriginalBin_Pt 
                                    );
        } else {
            //subtract secondary contribution from cocktail from inclusive spectrum
            cout << "will use cocktail for secondary correction" << endl;
            CorrectGammaSecAndPurityCocktail(histoESDCaloGammaPt_OriginalBin,
                                             histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin,
                                             histoGammaTrueSecCocktailGammaFromX_PtOrBin[3],
                                             histoGammaTruePurity_OriginalBin_Pt );
            CorrectGammaSecAndPurityCocktail(histoMCrecGammaCalo_OriginalBin_Pt,
                                             histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin,
                                             histoGammaTrueSecCocktailGammaFromX_PtOrBin[3],
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
        CorrectGammaUnfoldResol(    histoGammaCorrUnfoldReso_Pt, 
                                    histoGammaPrimaryRecoEff_MCPt, 
                                    deltaEtaCalo, 
                                    scalingCalo, 
                                    nEvt
                               );
        CorrectGammaUnfoldResol(    histoMCrecGammaCorr_Pt, 
                                    histoGammaPrimaryRecoEff_MCPt, 
                                    deltaEtaCalo, 
                                    scalingCalo,
                                    nEvt
                               );
        CorrectGammaUnfoldResol(    histoGammaCorrUnfoldReso_BinByBin_Pt, 
                                    histoGammaPrimaryRecoEff_MCPt,
                                    deltaEtaCalo, 
                                    scalingCalo, 
                                    nEvt
                               );

        SetHistogramm(histoGammaCorrUnfoldReso_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", -99, -99, 1.0, 1.7);
        SetHistogramm(histoGammaCorrUnfoldReso_BinByBin_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", -99, -99, 1.0, 1.7);
        SetHistogramm(histoMCrecGammaCorr_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", -99, -99, 1.0, 1.7);
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
    
    //**********************************************************************************
    //******************** Corrected Photon Spectrum Plot ******************************
    //**********************************************************************************
    TCanvas *canvasCorrGammaSpecUnfold     = new TCanvas("canvasDecayGammaSpecMC","",200,10,1000*1.25,1100*1.25);  // gives the page size
    DrawGammaCanvasSettings( canvasCorrGammaSpecUnfold, 0.16, 0.02, 0.02, 0.09);
    canvasCorrGammaSpecUnfold->SetLogy();
    
        TH1D* Dummy                             = NULL;
        if (isPCM && !isCalo) Dummy             = (TH1D*)histoESDConvGammaPt->Clone("Dummy");
        if (isCalo && !isPCM) Dummy             = (TH1D*)histoESDCaloGammaPt->Clone("Dummy");
    
        if (isPCM && !isCalo) SetHistogramm(Dummy,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",1e-8,10,1.0,1.7);
        if (isCalo && !isPCM) SetHistogramm(Dummy,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",1e-10,1e-2,1.0,1.7);
        Dummy->Draw();
        histoGammaCorrUnfoldReso_Pt->DrawCopy("same");
        histoGammaCorrUnfoldReso_BinByBin_Pt->DrawCopy("same");
        //histoGammaCorrUnfoldReso_SvD_Pt->DrawCopy("same");
        //histoGammaCorrUnfoldReso_TUnfold_Pt->DrawCopy("same");

        TLegend* legendGammaSpecUnfoldComp      = GetAndSetLegend2(0.6,0.945-2*1.1*0.035, 0.95,0.945, 0.035, 1, "", 42, 0.1); 
        legendGammaSpecUnfoldComp->AddEntry(histoGammaCorrUnfoldReso_Pt,"Bayesian unfolding");
        legendGammaSpecUnfoldComp->AddEntry(histoGammaCorrUnfoldReso_BinByBin_Pt,"Bin-by-Bin unfolding");
        legendGammaSpecUnfoldComp->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.935, 0.945-2*1.1*0.035, 0.035, cent, detectionProcess,"", 42, 0.03, "", 1, 1.25,31);
        
    canvasCorrGammaSpecUnfold->SaveAs(Form("%s/%s_%s_CorrGammaUnfoldSpectrumPurityMinusSec_%s.%s",outputDir.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasCorrGammaSpecUnfold;

    //**********************************************************************************
    //******************** Resolution Correction Plot **********************************
    //**********************************************************************************
    TCanvas *canvasResolutionCorr               = GetAndSetCanvas("canvasResolutionCorr");
        TH1D* Dummy2                            = NULL;
        if (isPCM && !isCalo) Dummy2            = (TH1D*)histoESDConvGammaPt->Clone("Dummy2");
        if (isCalo && !isPCM) Dummy2            = (TH1D*)histoESDCaloGammaPt->Clone("Dummy2");
        SetHistogramm(Dummy2,"#it{p}_{T} (GeV/#it{c})", "resolution correction",0,2);
        Dummy2->Draw();

        DrawGammaSetMarker(histoGammaResolCorrUnfold_Pt, 20, 1.0, kBlue+1, kBlue+1);
        DrawGammaSetMarker(histoGammaResolCorrUnfold_BinByBin_Pt, 21, 1.0, kGray+2, kGray+2);
        DrawGammaSetMarker(histoGammaResolCorrEff_Pt, 24, 1.5, kRed+2, kRed+2);

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

        PutProcessLabelAndEnergyOnPlot( 0.18, 0.25, 0.035, cent, detectionProcess,"", 42, 0.03);
 
    canvasResolutionCorr->SaveAs(Form("%s/%s_%s_ResolutionCorrection_%s.%s",outputDir.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasResolutionCorr;
    
    //******************************************************************************************
    // copy raw inc gamma spectrum and correct for;
    // - secondary contamination, 
    // - purity, 
    // - conversion probability & 
    // - reco effi (vs rec pt with sec correction)
    //****************************************************************************************** 
    TH1D* histoGammaCorrEffiReso_Pt                 = NULL;
    // correct alternate way for conversion reco
    if (isPCM && !isCalo) {
        histoGammaCorrEffiReso_Pt                   = (TH1D*)histoESDConvGammaPt->Clone("CorrGammaSpecPurityMinusSec");
        if (!hasCocktailInput)
            CorrectGammaEffiResol(  histoGammaCorrEffiReso_Pt, 
                                    histoSecondaryGammaFromXSpecPt[3], 
                                    histoSecondaryGammaFromXSpecPt[0], 
                                    histoGammaTruePurity_Pt, 
                                    histoGammaConvProb_MCPt, 
                                    histoGammaPrimaryRecoEff_Pt, 
                                    deltaEta, scaling, nEvt
                                 );
        else
            CorrectGammaEffiResolCocktail(  histoGammaCorrEffiReso_Pt, 
                                            histoGammaSecGammaFromX_Cocktail_Raw_Pt, 
                                            histoGammaTrueSecCocktailGammaFromX_Pt[3], 
                                            histoGammaTruePurity_Pt, 
                                            histoGammaConvProb_MCPt,
                                            histoGammaPrimaryRecoEff_Pt, 
                                            deltaEta, scaling, nEvt
                                         );
    }
    // correct alternate way for calo reco
    if (isCalo && !isPCM) {
        histoGammaCorrEffiReso_Pt                   = (TH1D*)histoESDCaloGammaPt->Clone("CorrGammaSpecPurityMinusSec");
        
        if (!hasCocktailInput)
            CorrectGammaEffiResol(  histoGammaCorrEffiReso_Pt, 
                                    histoSecondaryGammaFromXSpecPt[3], 
                                    histoSecondaryGammaFromXSpecPt[0],
                                    histoGammaTruePurity_Pt,
                                    histoGammaPrimaryRecoEff_Pt, 
                                    deltaEtaCalo, scalingCalo, nEvt
                                 );
        else
            CorrectGammaEffiResolCocktail(histoGammaCorrEffiReso_Pt,
                                          histoGammaSecGammaFromX_Cocktail_Raw_Pt, 
                                          histoGammaTrueSecCocktailGammaFromX_Pt[3], 
                                          histoGammaTruePurity_Pt, 
                                          histoGammaPrimaryRecoEff_Pt, 
                                          deltaEtaCalo, scalingCalo, nEvt
                                         );
    }
    DrawGammaSetMarker(histoGammaCorrEffiReso_Pt, 20, 1.0, kGreen+2, kGreen+2);
    SetHistogramm(histoGammaCorrEffiReso_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",-99,-99, 1.0,1.7);

    //****************************************************************************************** 
    TH1D* histoGammaCorrEffiReso_PileUp_Pt              = NULL;
    TH1D* histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt    = NULL;
    if(doPileUpCorr){
        //****************************************************************************************** 
        // copy raw inc gamma spectrum after pure data driven pileup correction,
        // - secondary contamination, 
        // - pileup corrected purity, (pileup corrected, vs rec pT with sec correction) 
        // - conversion probability 
        // - reco effi (vs rec pt with sec correction)
        //****************************************************************************************** 
        histoGammaCorrEffiReso_PileUp_Pt                = (TH1D*)histoESDConvGammaPtPileUp->Clone("CorrGammaSpecPurityMinusSecPileUp");
        
        if (!hasCocktailInput) {
            CorrectGammaEffiResol(  histoGammaCorrEffiReso_PileUp_Pt,
                                    histoSecondaryGammaSpecPtPileUp, 
                                    histoSecondaryGammaFromXFromK0sSpecPtPileUp,
                                    histoGammaTruePurity_PileUp_Pt, 
                                    histoGammaConvProb_MCPt, 
                                    histoGammaPrimaryRecoEff_PileUp_Pt,
                                    deltaEta, scaling, nEvt
                                 );
        } else {
            cout << __LINE__ << ": MC updated bin-by-bin correction for pileup corrected spectra not implemented yet for use of cocktail-based secondary correction, will use MC-based secondary correction." << endl;
            CorrectGammaEffiResol(  histoGammaCorrEffiReso_PileUp_Pt,
                                    histoSecondaryGammaSpecPtPileUp, 
                                    histoSecondaryGammaFromXFromK0sSpecPtPileUp,
                                    histoGammaTruePurity_PileUp_Pt,
                                    histoGammaConvProb_MCPt, 
                                    histoGammaPrimaryRecoEff_PileUp_Pt, 
                                    deltaEta, scaling, nEvt
                                 );
        }
        
        DrawGammaSetMarker(histoGammaCorrEffiReso_PileUp_Pt, 20, 1.0, kMagenta, kMagenta);
        SetHistogramm(histoGammaCorrEffiReso_PileUp_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", -99, -99, 1.0, 1.7);

        //****************************************************************************************** 
        // copy raw inc gamma spectrum after pure data driven pileup correction,
        // - secondary contamination, 
        // - purity, 
        // - conversion probability 
        // - reco effi (vs rec pt with sec correction)
        //****************************************************************************************** 
        histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt    = (TH1D*)histoESDConvGammaPtPileUp->Clone("CorrGammaSpecPurityMinusSecPileUpNoMCUpdate");
        
        if (!hasCocktailInput) {
            CorrectGammaEffiResol(  histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt,
                                    histoSecondaryGammaFromXSpecPt[3], 
                                    histoSecondaryGammaFromXSpecPt[0],
                                    histoGammaTruePurity_Pt, 
                                    histoGammaConvProb_MCPt, 
                                    histoGammaPrimaryRecoEff_Pt, 
                                    deltaEta, scaling, nEvt
                                 );
        } else {
            CorrectGammaEffiResolCocktail( histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt, 
                                           histoGammaSecGammaFromX_Cocktail_Raw_Pt, 
                                           histoGammaTrueSecCocktailGammaFromX_Pt[3], 
                                           histoGammaTruePurity_Pt, 
                                           histoGammaConvProb_MCPt, 
                                           histoGammaPrimaryRecoEff_Pt, 
                                           deltaEta, scaling, nEvt
                                         );
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
    SetHistogramm(histoMCGammaSpec_MCPt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", -99, -99, 1.0, 1.7);

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
 
        PutProcessLabelAndEnergyOnPlot( 0.95, 0.95, 0.035, cent, "#gamma", "", 42, 0.03,"", 1, 1.25,31);
        
    canvasRatioAllDiffDecay->SaveAs(Form("%s/%s_AllGammaDivDecayGammaSpectrumMC_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasRatioAllDiffDecay;
   
    //*************************************************************************************************
    //********************************* Gamma from Decay **********************************************
    //*************************************************************************************************
    // correct input spectra of different decay channels with delta eta, delta phi & events number
    for (Int_t i = 0; i< 9; i++){
        CorrectGammaMC(histoPhotonSource_MCPt[i], deltaEta, scaling, nEvtMC);
        SetHistogramm(histoPhotonSource_MCPt[i],"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",-99., -99., 1.0, 1.65);        
        DrawGammaSetMarker(histoPhotonSource_MCPt[i], 22, 1.0, colorsDecay[i] , colorsDecay[i]);
    }
        
    TCanvas *canvasDecayGammaSpecMC     = new TCanvas("canvasDecayGammaSpecMC","",200,10,1000*1.25,1100*1.25);  // gives the page size
    DrawGammaCanvasSettings( canvasDecayGammaSpecMC, 0.16, 0.02, 0.02, 0.09);
    canvasDecayGammaSpecMC->SetLogy();
    
        TLegend* legendDecaySpectra     = GetAndSetLegend2(0.5, 0.935-0.035*1.1*5, 0.95, 0.935, 0.035, 2, "", 42, 0.2); 
        histoPhotonSource_MCPt[7]->DrawCopy("chist");
        for (Int_t i = 0; i< 9; i++){   
            histoPhotonSource_MCPt[i]->DrawCopy("chist,same");
            legendDecaySpectra->AddEntry(histoPhotonSource_MCPt[i],decaysLabels[i].Data(),"l");
        }
        legendDecaySpectra->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.94, 0.935-0.035*1.1*5, 0.035, cent, detectionProcess,"", 42, 0.03,"", 1, 1.25,31);
        
    canvasDecayGammaSpecMC->SaveAs(Form("%s/%s_DecayGammaSpectrumMC_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasDecayGammaSpecMC;
    
    //*************************************************************************************************
    //**************** Compare cocktail secondaries to MC secondaries *********************************
    //*************************************************************************************************
    if(hasCocktailInput){
        TH1D* histoRatioSecondariesCocktailMC_Raw_Pt[3]     = {NULL, NULL, NULL};
        TCanvas *canvasSecondaryComparison                  = GetAndSetCanvas("canvasSecondaryComparison");
        DrawGammaCanvasSettings( canvasDecayGammaSpecMC, 0.07, 0.02, 0.02, 0.085);
        TLegend* legendCompareSecCocktailMC                 = GetAndSetLegend2(0.75, 0.935-0.035*1.1*3, 0.95, 0.935, 0.035, 1, "", 42, 0.15); 
        legendCompareSecCocktailMC->SetBorderSize(0);
        
        for (Int_t k = 0; k< 3; k++){
            histoRatioSecondariesCocktailMC_Raw_Pt[k]       = (TH1D*)histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]->Clone(Form("histoRatioSecondariesCocktailMC%s_Raw_Pt", nameSecondaries[k].Data()));
            if (isPCM && !isCalo) histoRatioSecondariesCocktailMC_Raw_Pt[k]->Divide(histoGammaTrueSecConvGammaFromX_Pt[k]);
            if (isCalo && !isPCM) histoRatioSecondariesCocktailMC_Raw_Pt[k]->Divide(histoGammaTrueSecCaloGammaFromX_Pt[k]);
            SetHistogramm(histoRatioSecondariesCocktailMC_Raw_Pt[k],"#it{p}_{T} (GeV/#it{c})", "Cocktail/MC",-99,-99,1.0,0.9);
            DrawGammaSetMarker(histoRatioSecondariesCocktailMC_Raw_Pt[k], markerStyleSec[k], markerSizeSec[k]*2, colorSecFromToy[k] , colorSecFromToy[k]);
            if (k == 0)
                histoRatioSecondariesCocktailMC_Raw_Pt[k]->Draw();
            else 
                histoRatioSecondariesCocktailMC_Raw_Pt[k]->Draw("same");
            legendCompareSecCocktailMC->AddEntry(histoRatioSecondariesCocktailMC_Raw_Pt[k],Form("Sec. #gamma from %s", nameLabelSecondaries[k].Data()),"p");
        }
        legendCompareSecCocktailMC->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.935, 0.24, 0.035, cent, detectionProcess,"", 42, 0.03,"", 1, 1.25,31);
        
        canvasSecondaryComparison->SaveAs(Form("%s/%s_SecondaryComparisonCocktailMC_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasSecondaryComparison;
    }
    
    //*************************************************************************************************
    //******************************* Gamma Combinatorial Background **********************************
    //*************************************************************************************************
    // pure BG plotting
    TCanvas *canvasCombBackSpecMC     = new TCanvas("canvasCombBackSpecMC","",200,10,1000*1.25,1100*1.25);  // gives the page size
    DrawGammaCanvasSettings( canvasCombBackSpecMC, 0.16, 0.02, 0.02, 0.09);
    canvasCombBackSpecMC->SetLogy();
    
        TLegend* legendCombSpectra = GetAndSetLegend2(0.35,0.945-4*1.1*0.035, 0.95,0.945, 0.035, 5, "", 42, 0.1);
        if (isPCM && !isCalo) {
            for(Int_t i = 0;i<17;i++){
                SetHistogramm(histoCombinatorialSpecies_Pt[i],"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",-99., -99., 1.0, 1.7);        
                DrawGammaSetMarker(histoCombinatorialSpecies_Pt[i], markersCombinatorics[i], 1., colorsCombinatorics[i], colorsCombinatorics[i]);
                histoCombinatorialSpecies_Pt[i]->Scale(1./nEvtMC);
            
                if(i==0) histoCombinatorialSpecies_Pt[i]->DrawCopy("");
                else histoCombinatorialSpecies_Pt[i]->DrawCopy("same");
            
                legendCombSpectra->AddEntry(histoCombinatorialSpecies_Pt[i],combinatoricsLabels[i]);
            }
        }
        if (isCalo && !isPCM) {
            for(Int_t i = 0;i<11;i++){
                SetHistogramm(histoCombinatorialSpeciesCalo_Pt[i],"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",-99., -99., 1.0, 1.7);
                DrawGammaSetMarker(histoCombinatorialSpeciesCalo_Pt[i], markersCombinatoricsCalo[i], 1., colorsCombinatoricsCalo[i], colorsCombinatoricsCalo[i]);
                histoCombinatorialSpeciesCalo_Pt[i]->Scale(1./nEvtMC);
                
                if(i==0) {
                    histoCombinatorialSpeciesCalo_Pt[i]->GetYaxis()->SetRangeUser(1e-10, 1e-3);
                    histoCombinatorialSpeciesCalo_Pt[i]->DrawCopy("");
                }
                else histoCombinatorialSpeciesCalo_Pt[i]->DrawCopy("same");
                
                legendCombSpectra->AddEntry(histoCombinatorialSpeciesCalo_Pt[i],combinatoricsLabelsCalo[i]);
            }
        }

        legendCombSpectra->Draw();
        PutProcessLabelAndEnergyOnPlot( 0.935, 0.93-4*1.1*0.035, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,31);
    
    canvasCombBackSpecMC->SaveAs(Form("%s/%s_CombBackSpectrumMC_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasCombBackSpecMC;
    
    // plot ratio to real primary photons
    TCanvas *canvasSignalToCombBackgroundRatio      = GetAndSetCanvas("canvasSignalToCombBackgroundRatio");
    canvasSignalToCombBackgroundRatio->SetLogy();

        TLegend* legendSignalToCombBackgroundRatio  = GetAndSetLegend2(0.15,0.93-2*1.1*0.035, 0.95,0.93, 0.035, 6, "", 42, 0.1);
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
                
                legendSignalToCombBackgroundRatio->AddEntry(histoSignalToCombBackgroundRatio[i],combinatoricsLabels[i]);
            }
            
            SummedSmallContributionsCombBack->Divide(SummedSmallContributionsCombBack,histoGammaTruePrimaryConv_Pt,1,1,"");
            legendSignalToCombBackgroundRatio->AddEntry(SummedSmallContributionsCombBack,"p(#bar{p})K^{#pm}#mu^{#pm}");
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
            
                legendSignalToCombBackgroundRatio->AddEntry(histoSignalToCombBackgroundRatio[i],combinatoricsLabelsCalo[i]);
            }
            
            SummedSmallContributionsCombBack->Divide(SummedSmallContributionsCombBack,histoGammaTruePrimaryCalo_Pt,1,1,"");
            legendSignalToCombBackgroundRatio->AddEntry(SummedSmallContributionsCombBack,"p+K0s+Lambda+mu+rest");
        }
   
        DrawGammaSetMarker(SummedSmallContributionsCombBack, markersCombinatorics[10], 1., colorsCombinatorics[10], colorsCombinatorics[10]);
        SummedSmallContributionsCombBack->DrawCopy("e1same");
        legendSignalToCombBackgroundRatio->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.15, 0.93-0.035*1.05*2, 0.035, cent, detectionProcess, "" , 42, 0.03,"",1,1.25,11);
        
    canvasSignalToCombBackgroundRatio->SaveAs(Form("%s/%s_CombBackgroundRatioToSignal_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasSignalToCombBackgroundRatio;

    // plot ratio to summed total MC BG
    TCanvas *canvasRatioCombBackToBack              = GetAndSetCanvas("canvasRatioCombBackToBack");

    canvasRatioCombBackToBack->SetLogy();
        TLegend* legendRatioCombBackToBack          = GetAndSetLegend2(0.15,0.93-2*1.1*0.035, 0.95,0.93, 0.035, 6, "", 42, 0.1);
        
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
                    histoRatioCombBackToBack[i]->GetYaxis()->SetRangeUser(1e-4,40);
                    histoRatioCombBackToBack[i]->DrawCopy("e1");
                }
                if(i<9){
                    DrawGammaSetMarker(histoRatioCombBackToBack[i], markersCombinatorics[i], 1., colorsCombinatorics[i], colorsCombinatorics[i]);
                    histoRatioCombBackToBack[i]->DrawCopy("e1same");
                }
                else continue;
            
                legendRatioCombBackToBack->AddEntry(histoRatioCombBackToBack[i],combinatoricsLabels[i]);
            }
            
            legendRatioCombBackToBack->AddEntry(SummedSmallContributionsCombBackToBack,"p(#bar{p})K^{#pm}#mu^{#pm}");
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
                    histoRatioCombBackToBack[i]->GetYaxis()->SetRangeUser(1e-4,40);
                    histoRatioCombBackToBack[i]->DrawCopy("e1");
                }
                if(i<2 || (i>2 && i<5) || i==8){
                    DrawGammaSetMarker(histoRatioCombBackToBack[i], markersCombinatoricsCalo[i], 1., colorsCombinatoricsCalo[i], colorsCombinatoricsCalo[i]);
                    histoRatioCombBackToBack[i]->DrawCopy("e1same");
                } else continue;
                
                legendRatioCombBackToBack->AddEntry(histoRatioCombBackToBack[i],combinatoricsLabelsCalo[i]);
            }
            
            legendRatioCombBackToBack->AddEntry(SummedSmallContributionsCombBackToBack,"p+K0s+Lambda+mu+rest");
        }
    
        SummedSmallContributionsCombBackToBack->Divide(SummedSmallContributionsCombBackToBack,histoGammaMCBackground_Pt,1,1,"B");
        DrawGammaSetMarker(SummedSmallContributionsCombBackToBack, markersCombinatorics[10], 1., colorsCombinatorics[10], colorsCombinatorics[10]);
        SummedSmallContributionsCombBackToBack->DrawCopy("e1same");
        legendRatioCombBackToBack->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.15, 0.93-0.035*1.05*2, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,11);
    
    canvasRatioCombBackToBack->SaveAs(Form("%s/%s_RatioCombBackToBack_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasRatioCombBackToBack;
    
    //*************************************************************************************************
    //***************** Compare Gamma Yields alternate corr method vs unfolded ************************
    //*************************************************************************************************
    TCanvas *canvaskGammaSpecAlternateCorrMethods   = new TCanvas("canvaskGammaSpecAlternateCorrMethods","",200,10,1000*1.25,1300*1.25);  // gives the page size
    DrawGammaCanvasSettings(canvaskGammaSpecAlternateCorrMethods,0.16, 0.015, 0.02, 0.09);
    TPad* padAlternateCorrMethods       = new TPad("padAlternateCorrMethods", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings(padAlternateCorrMethods, 0.16, 0.015, 0.02, 0.);
    padAlternateCorrMethods->Draw();
    TPad* padAlternateCorrMethodsRatio  = new TPad("padAlternateCorrMethodsRatio", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings(padAlternateCorrMethodsRatio,  0.16, 0.015, 0.0, 0.27);
    padAlternateCorrMethodsRatio->Draw();

    padAlternateCorrMethodsRatio->Draw();
    padAlternateCorrMethods->cd();
    padAlternateCorrMethods->SetLogy();

        TLegend* legendAlternateCorrMethods = GetAndSetLegend(0.25,0.85,2);
        // plot unfolded spectrum
        histoGammaCorrUnfoldReso_Pt->DrawCopy("e1");
        legendAlternateCorrMethods->AddEntry(histoGammaCorrUnfoldReso_Pt,"#gamma with Bayes unfolding","p");
        // copy pileup correct hist if pileup correction enabled
        if (doPileUpCorr){
            DrawGammaSetMarker(histoGammaCorrEffiReso_PileUp_Pt,  24, 1.0, kAzure+7,  kAzure+7);
            histoGammaCorrEffiReso_PileUp_Pt->DrawCopy("e1,same");
            legendAlternateCorrMethods->AddEntry(histoGammaCorrEffiReso_PileUp_Pt,"#gamma with bin-by-bin MC based corrections","p");
        // take hist without pileup corr if disabled
        } else {
            DrawGammaSetMarker(histoGammaCorrEffiReso_Pt,  24, 1.0, kAzure+7,  kAzure+7);
            histoGammaCorrEffiReso_Pt->DrawCopy("e1,same");
            legendAlternateCorrMethods->AddEntry(histoGammaCorrEffiReso_Pt,"#gamma with bin-by-bin MC based corrections","p");
        }    
        legendAlternateCorrMethods->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.20, 0.15, 0.035, cent, detectionProcess,"", 42, 0.03);
        
    padAlternateCorrMethodsRatio->cd();
        TH1D *histoRatioGammaLegacyCorrVsUnfold = NULL;
        if (doPileUpCorr)
            histoRatioGammaLegacyCorrVsUnfold = (TH1D*) histoGammaCorrEffiReso_PileUp_Pt->Clone("histoRatioGammaLegacyCorrVsUnfold");
        else     
            histoRatioGammaLegacyCorrVsUnfold = (TH1D*) histoGammaCorrEffiReso_Pt->Clone("histoRatioGammaLegacyCorrVsUnfold");
        histoRatioGammaLegacyCorrVsUnfold->Divide(histoRatioGammaLegacyCorrVsUnfold,histoGammaCorrUnfoldReso_Pt,1,1,"");
        histoRatioGammaLegacyCorrVsUnfold->GetYaxis()->SetRangeUser(0.88,1.12);
        SetStyleHistoTH1ForGraphs(  histoRatioGammaLegacyCorrVsUnfold, "#it{p}_{T} (GeV/#it{c})","#frac{#gamma corr legacy}{#gamma corr unfold}" , 0.10, 0.12,  
                                    0.10, 0.12,  0.85, 0.6, 510, 505);
        
        DrawGammaSetMarker(histoRatioGammaLegacyCorrVsUnfold, 24, 1.0, kAzure+7,  kAzure+7);
        histoRatioGammaLegacyCorrVsUnfold->DrawCopy("");
        
        DrawGammaLines(0., maxPtGamma,1, 1,0.5);
    
    canvaskGammaSpecAlternateCorrMethods->SaveAs(Form("%s/%s_%s_GammaSpectraComparison_%s.%s",outputDir.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
    delete canvaskGammaSpecAlternateCorrMethods;


    //*************************************************************************************************
    //***************** Compare Gamma reconstructed gamma yields to MC input **************************
    //*************************************************************************************************    
    if (isPCM && !isCalo) {
        CorrectGammaMC(histoMCGammaConv_MCPt, deltaEta, scaling, nEvtMC);
        CorrectGammaMC(histoGammaTrueConv_Pt, deltaEta, scaling, nEvtMC);
    }
    if (isCalo && !isPCM) {
        CorrectGammaMC(histoGammaTrueCalo_Pt, deltaEtaCalo, scalingCalo, nEvtMC);
    }

    TCanvas *canvasGammaPurityConvBin   = new TCanvas("canvasGammaPurityConvBin","",200,10,1000*1.25,1300*1.25);  // gives the page size
    DrawGammaCanvasSettings(canvasGammaPurityConvBin,0.16, 0.015, 0.02, 0.09);
    TPad* padGammaPurityConvBin       = new TPad("padGammaPurityConvBin", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings(padGammaPurityConvBin, 0.16, 0.015, 0.02, 0.);
    padGammaPurityConvBin->Draw();
    TPad* padGammaPurityConvBinRatio  = new TPad("padGammaPurityConvBinRatio", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings(padGammaPurityConvBinRatio,  0.16, 0.015, 0.0, 0.27);
    padGammaPurityConvBinRatio->Draw();

    padGammaPurityConvBin->cd();
    padGammaPurityConvBin->SetLogy();
    
        DrawGammaSetMarker(histoGammaCorrEffiReso_Pt, 24, 1.0, 1, 1);
        histoMCGammaSpec_MCPt->GetYaxis()->SetRangeUser(5e-10,100);
        histoMCGammaSpec_MCPt->DrawCopy("e1");
    
        DrawGammaSetMarker(histoGammaCorrUnfoldReso_Pt, 24, 1.0, 1, 1);
        histoGammaCorrUnfoldReso_Pt->DrawCopy("e1,same");
    
        TLegend *legendGammaSpectraConvBin = GetAndSetLegend(0.45,0.85,2);
        legendGammaSpectraConvBin->AddEntry(histoGammaCorrUnfoldReso_Pt,"corrected data #gamma spectrum");
        legendGammaSpectraConvBin->AddEntry(histoMCGammaSpec_MCPt,"input MC #gamma spectrum");
        legendGammaSpectraConvBin->Draw();
    
        PutProcessLabelAndEnergyOnPlot( 0.20, 0.15, 0.035, cent, detectionProcess, "", 42, 0.03);
    
    padGammaPurityConvBinRatio->cd();

        TH1D* tempRatioDataMC               = NULL;
        tempRatioDataMC                     = (TH1D*)histoGammaCorrUnfoldReso_Pt->Clone("tempRatioDataMC");
        tempRatioDataMC->Divide(tempRatioDataMC,histoMCGammaSpec_MCPt,1.,1.,"");
    
        SetStyleHistoTH1ForGraphs(tempRatioDataMC, "#it{p}_{T} (GeV/#it{c})", "#frac{data}{MC}", 0.10, 0.14, 0.10, 0.14, 0.85, 0.5, 510, 505);
        DrawGammaSetMarker(tempRatioDataMC, 24, 1.0, kBlack, kBlack);
        if (isRunMC){
            tempRatioDataMC->GetYaxis()->SetRangeUser(0.85, 1.15);
        } else if (energy.Contains("pPb")) {
            tempRatioDataMC->GetYaxis()->SetRangeUser(0.5, 2.7);
        } else {
            tempRatioDataMC->GetYaxis()->SetRangeUser(0.5, 1.5);
        }
        tempRatioDataMC->DrawCopy("e1");
        DrawGammaLines(0., tempRatioDataMC->GetXaxis()->GetBinUpEdge(tempRatioDataMC->GetNbinsX()), 1., 1., 0.5, kBlue+1);
    
    canvasGammaPurityConvBin->SaveAs(Form("%s/%s_%s_GammaSpec_%s.%s",outputDir.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasGammaPurityConvBin;

    TH1D* histoGammaCorrFac_Pt           = (TH1D*) histoGammaPrimaryRecoEff_Pt->Clone("GammaCorrFac_Pt");
    if (isPCM && !isCalo) histoGammaCorrFac_Pt->Multiply(histoGammaConvProb_MCPt);
    
    //***********************************************************************************************************
    //******************************* Write output file for photons *********************************************
    //***********************************************************************************************************
    TString nameOutput = Form("%s/%s/%s_%s_GammaConvV1Correction_%s.root",cutSelection.Data(),energy.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data());
    TFile* fileCorrectedOutput = new TFile(nameOutput,"RECREATE");

        //________________________ writing MC quantities to file _____________________________
        // input spectrum corrected for deta, dPhi, Nevt
        if (histoMCGammaSpec_MCPt)                              histoMCGammaSpec_MCPt->Write("GammaSpecMC", TObject::kOverwrite);
        if (histoSecondaryGammaFromXSpecPt[3])                  histoSecondaryGammaFromXSpecPt[3]->Write("histoSecondaryGammaFromXSpecPtRest", TObject::kOverwrite);
        if (histoGammaTrueSecCocktailGammaFromX_Pt[3])          histoGammaTrueSecCocktailGammaFromX_Pt[3]->Write("histoGammaTrueSecCocktailGammaRest_Pt", TObject::kOverwrite);
    
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

    
        //_________________________ writing secondary correction factors to file _________________________
        for (Int_t k = 0; k < 3; k++){
            if (histoGammaSecFromXRecoEff_MCPt[k])                  histoGammaSecFromXRecoEff_MCPt[k]->Write(histoGammaSecFromXRecoEff_MCPt[k]->GetName(), TObject::kOverwrite);
            if (histoGammaSecFromXRecoEff_MCPtOrBin[k])             histoGammaSecFromXRecoEff_MCPtOrBin[k]->Write(histoGammaSecFromXRecoEff_MCPtOrBin[k]->GetName(), TObject::kOverwrite);
            if (histoGammaSecFromXRecoEff_RecPt[k])                 histoGammaSecFromXRecoEff_RecPt[k]->Write(histoGammaSecFromXRecoEff_RecPt[k]->GetName(), TObject::kOverwrite);
            if (histoGammaSecFromXRecoEff_RecPtOrBin[k])            histoGammaSecFromXRecoEff_RecPtOrBin[k]->Write(histoGammaSecFromXRecoEff_RecPtOrBin[k]->GetName(), TObject::kOverwrite);
            if (histoGammaSecondaryFromXConvProb_MCPt[k])           histoGammaSecondaryFromXConvProb_MCPt[k]->Write(histoGammaSecondaryFromXConvProb_MCPt[k]->GetName(), TObject::kOverwrite);
            if (histoGammaSecondaryFromXConvProb_MCPtOrBin[k])      histoGammaSecondaryFromXConvProb_MCPtOrBin[k]->Write(histoGammaSecondaryFromXConvProb_MCPtOrBin[k]->GetName(), TObject::kOverwrite);
            if (histoGammaSecGammaFromX_Cocktail_Raw_Pt[k])         histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]->Write(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]->GetName(), TObject::kOverwrite);
            if (histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin[k])    histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin[k]->Write(histoGammaSecGammaFromX_Cocktail_Raw_PtOrBin[k]->GetName(), TObject::kOverwrite);
            if (histoGammaTrueSecCocktailGammaFromX_Pt[k])          histoGammaTrueSecCocktailGammaFromX_Pt[k]->Write(Form("histoGammaTrueSecCocktailGammaFromXFrom%s_Pt", nameSecondaries[k].Data()),
                                                                                                                     TObject::kOverwrite);
            if (histoSecondaryGammaFromXSpecPt[k])                  histoSecondaryGammaFromXSpecPt[k]->Write(Form("histoSecondaryGammaFromXFrom%sSpecPt", nameSecondaries[k].Data()), TObject::kOverwrite);
        }
        
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
