//  **********************************************************************************
//  ******     provided by Gamma Conversion Group, PWGGA,                        *****
//  ******     Friederike Bock, friederike.bock@cern.ch                          *****
//  **********************************************************************************

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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
// #include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

Bool_t kIsMC = kFALSE;

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    // TString name;
};

Double_t FindLargestBin1DHist(TH1* hist ){
    Double_t largestContent     = 0;
    for (Int_t i= 0; i < hist->GetNbinsX(); i++){
        if (largestContent < hist->GetBinContent(i)){
            largestContent = hist->GetBinContent(i);
        }    
    }    
    return largestContent;
}

Double_t FindSmallestBin1DHist(TH1* hist, Double_t maxStart = 1e6 ){
    Double_t smallesContent     = maxStart;
    for (Int_t i= 0; i < hist->GetNbinsX(); i++){
        if (hist->GetBinContent(i) != 0 && smallesContent > hist->GetBinContent(i)){
            smallesContent = hist->GetBinContent(i);
        }    
    }    
    return smallesContent;
}


void CorrectYieldDalitz(TH1D* histoCorrectedYield,TH1D* histoRawGGYield, TH1D* histoEffiPt, TH1D* histoAcceptance, Double_t deltaRapid, Double_t scaling, Double_t nEvt, TString nameMeson){
    histoCorrectedYield->Sumw2();
    histoCorrectedYield->Add(histoRawGGYield,-1.);
    histoCorrectedYield->Scale(1./nEvt);
    histoCorrectedYield->Divide(histoCorrectedYield,histoEffiPt,1.,1.,"");
    histoCorrectedYield->Divide(histoCorrectedYield,histoAcceptance,1.,1.,"");
    histoCorrectedYield->Scale(1./deltaRapid);
    histoCorrectedYield->Scale(scaling);
    for (Int_t i = 1; i < histoCorrectedYield->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoCorrectedYield->GetBinContent(i)/histoCorrectedYield->GetBinCenter(i);
        Double_t newBinError = histoCorrectedYield->GetBinError(i)/histoCorrectedYield->GetBinCenter(i);
        histoCorrectedYield->SetBinContent(i,newBinContent);
        histoCorrectedYield->SetBinError(i,newBinError);
    }
    if (nameMeson.CompareTo("Pi0") == 0 ||nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
        histoCorrectedYield->Scale(1./0.01198);
    }else{
        histoCorrectedYield->Scale(1./0.000068);
    }
}


void CorrectYield(TH1D* histoCorrectedYield, TH1D** histoRawSecYield, TH1D** histoRawSecYieldFromToy, TH1D* histoEffiPt, TH1D* histoAcceptance, Double_t deltaRapid, Double_t scaling, Double_t nEvt, TString nameMeson){
    histoCorrectedYield->Sumw2();
    
    // subtract standard secondary yield before event scaling
    for (Int_t j = 0; j < 4; j++){
        if (kIsMC){
            if (histoRawSecYield[j]) histoCorrectedYield->Add(histoRawSecYield[j],-1.);
        } else {
            if (histoRawSecYieldFromToy[j] && j < 3){
                cout << "will take secondary yield from toy approach for component: " << j << endl;
            } else {
                if (histoRawSecYield[j]) histoCorrectedYield->Add(histoRawSecYield[j],-1.);
            }    
        }    
    }    
    // scale with number of events
    histoCorrectedYield->Scale(1./nEvt);
    // subtract secondary yield from toy approach after scaling with nuber of events
    if (!kIsMC){
        for (Int_t j = 0; j < 3; j++){
            if (histoRawSecYieldFromToy[j]) histoCorrectedYield->Add(histoRawSecYieldFromToy[j],-1.);
        }    
    }    
    // correct with acceptance and efficiency
    histoCorrectedYield->Divide(histoCorrectedYield,histoAcceptance,1.,1.,"");
    histoCorrectedYield->Divide(histoCorrectedYield,histoEffiPt,1.,1.,"");
    // scale with 1/ (Delta rapidiy)
    histoCorrectedYield->Scale(1./deltaRapid);
    // scale with 1/(2 pi)
    histoCorrectedYield->Scale(scaling);
    // scale with 1/pT
    for (Int_t i = 1; i < histoCorrectedYield->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoCorrectedYield->GetBinContent(i)/histoCorrectedYield->GetBinCenter(i);
        Double_t newBinError = histoCorrectedYield->GetBinError(i)/histoCorrectedYield->GetBinCenter(i);
        histoCorrectedYield->SetBinContent(i,newBinContent);
        histoCorrectedYield->SetBinError(i,newBinError);
    }
    // scale with 1/ BR
    if (nameMeson.CompareTo("Pi0") == 0 ||nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
        histoCorrectedYield->Scale(1./0.98798);
    }else{
        histoCorrectedYield->Scale(1./0.3931);
    }
}

void CorrectYieldInclResonanceFeedDown(TH1D* histoCorrectedYield, TH1D** histoRawSecYield, TH1D** histoRawSecYieldFromToy, TH1D** histoResonanceFeedDownPi0, TH1D* histoEffiPt, TH1D* histoAcceptance, Double_t deltaRapid, Double_t scaling, Double_t nEvt, TString nameMeson){
    histoCorrectedYield->Sumw2();
    
    // subtract standard secondary yield before event scaling
    for (Int_t j = 0; j < 4; j++){
        if (kIsMC){
            if (histoRawSecYield[j]) histoCorrectedYield->Add(histoRawSecYield[j],-1.);
        } else {
            if (histoRawSecYieldFromToy[j] && j < 3){
                cout << "will take secondary yield from toy approach for component: " << j << endl;
            } else {
                if (histoRawSecYield[j]) histoCorrectedYield->Add(histoRawSecYield[j],-1.);
            }
        }
    }
    // scale with number of events
    histoCorrectedYield->Scale(1./nEvt);
    // subtract secondary yield from toy approach after scaling with nuber of events
    if (!kIsMC){
        for (Int_t j = 0; j < 3; j++){
            if (histoRawSecYieldFromToy[j]) histoCorrectedYield->Add(histoRawSecYieldFromToy[j],-1.);
        }
    }
    // correct with acceptance and efficiency
    histoCorrectedYield->Divide(histoCorrectedYield,histoAcceptance,1.,1.,"");
    histoCorrectedYield->Divide(histoCorrectedYield,histoEffiPt,1.,1.,"");
    
    // subtract contribution from resonance feed down
    if (!kIsMC) {
        for (Int_t j = 0; j < 8; j++) {
            if (histoResonanceFeedDownPi0[j]) histoCorrectedYield->Add(histoResonanceFeedDownPi0[j],-1.);
        }
    }
    
    // scale with 1/ (Delta rapidiy)
    histoCorrectedYield->Scale(1./deltaRapid);
    
    // scale with 1/(2 pi)
    histoCorrectedYield->Scale(scaling);
    // scale with 1/pT
    for (Int_t i = 1; i < histoCorrectedYield->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoCorrectedYield->GetBinContent(i)/histoCorrectedYield->GetBinCenter(i);
        Double_t newBinError = histoCorrectedYield->GetBinError(i)/histoCorrectedYield->GetBinCenter(i);
        histoCorrectedYield->SetBinContent(i,newBinContent);
        histoCorrectedYield->SetBinError(i,newBinError);
    }
    // scale with 1/ BR
    if (nameMeson.CompareTo("Pi0") == 0 ||nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
        histoCorrectedYield->Scale(1./0.98798);
    } else {
        histoCorrectedYield->Scale(1./0.3931);
    }
}

void CompileFullCorrectionFactor(TH1D* histoEffiPt, TH1D* histoAcceptance, Double_t deltaRapid){
    histoEffiPt->Sumw2();
    histoEffiPt->Multiply(histoEffiPt,histoAcceptance,1.,1.,"");
    histoEffiPt->Scale(deltaRapid);
}


void ScaleMCYield(TH1D* histoCorrectedToBeScaled, Double_t deltaRapid, Double_t scaling, Double_t nEvtMC, TString nameMeson, Bool_t optionDalitz ){
    histoCorrectedToBeScaled->Sumw2();
    histoCorrectedToBeScaled->Scale(1./deltaRapid);
    histoCorrectedToBeScaled->Scale(scaling);
    histoCorrectedToBeScaled->Scale(1./nEvtMC);
    for (Int_t i = 1; i < histoCorrectedToBeScaled->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoCorrectedToBeScaled->GetBinContent(i)/histoCorrectedToBeScaled->GetBinCenter(i);
        Double_t newBinError = histoCorrectedToBeScaled->GetBinError(i)/histoCorrectedToBeScaled->GetBinCenter(i);
        histoCorrectedToBeScaled->SetBinContent(i,newBinContent);
        histoCorrectedToBeScaled->SetBinError(i,newBinError);
    }
    if (nameMeson.CompareTo("Pi0") == 0 ||nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
        if (!optionDalitz){
            histoCorrectedToBeScaled->Scale(1./0.98798);
        } else {
            histoCorrectedToBeScaled->Scale(1./0.01198);
        }
    }else{
        if (!optionDalitz){
            histoCorrectedToBeScaled->Scale(1./0.3931);
        } else {
            histoCorrectedToBeScaled->Scale(1./6.8e-5);
        }
        
    }
}

// *******************************************************************************************************************************************
// ******************************** Main function ********************************************************************************************
// *******************************************************************************************************************************************
void  CorrectSignalV2(  TString fileNameUnCorrectedFile = "myOutput", 
                        TString fileNameCorrectionFile  = "", 
                        TString fCutSelection           = "", 
                        TString suffix                  = "gif", 
                        TString nameMeson               = "", 
                        TString isMC                    = "", 
                        TString optionEnergy            = "", 
                        TString optionPeriod            = "",
                        TString fEstimatePileup         = "", 
                        Bool_t optDalitz                = kFALSE, 
                        Int_t mode                      = 9,
                        Int_t triggerSet                = -1
                        
                     ){  
    
    // ******************************************************************************************
    // ********************** general style settings ********************************************
    // ******************************************************************************************
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    StyleSettingsThesis(suffix);  
    SetPlotStyle();
    
    TString fDalitz     ="";
    Bool_t kDalitz      = optDalitz;
    if (kDalitz){
        fDalitz         = "Dalitz";
    } else {
        fDalitz         = "";
    }
    
    // *******************************************************************************************
    // *********************** setting global variables ******************************************
    // *******************************************************************************************
    TString outputDir           = Form("%s/%s/%s/%s/CorrectSignalV2%s",fCutSelection.Data(),optionEnergy.Data(),optionPeriod.Data(),suffix.Data(), fDalitz.Data());
    gSystem->Exec("mkdir -p "+outputDir);
        
    TString date                = ReturnDateString();
    TString prefix2             = "";

    // plot labeling
    TString textProcess         = ReturnMesonString (nameMeson);
    if(textProcess.CompareTo("") == 0 ){
        cout << "Meson unknown" << endl;
        return ;
    }
    TString fTextMeasurement    = Form("%s #rightarrow #gamma#gamma", textProcess.Data());
    TString fDetectionProcess   = ReturnFullTextReconstructionProcess(mode);
    TString collisionSystem     = ReturnFullCollisionsSystem(optionEnergy);
    
    // cut strings
    TString fEventCutSelection      = "";
    TString fGammaCutSelection      = "";
    TString fClusterCutSelection    = "";
    TString fElectronCutSelection   = "";
    TString fMesonCutSelection      = "";
    
    if (mode == 9){
        if (kDalitz){
            ReturnSeparatedCutNumber(fCutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection,kTRUE);
        } else {
            ReturnSeparatedCutNumber(fCutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection);
        }
        fEventCutSelection          = fGammaCutSelection(0,7);
        fGammaCutSelection          = fGammaCutSelection(7,fGammaCutSelection.Length()-7);
        cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
    } else {
        ReturnSeparatedCutNumberAdvanced(fCutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
    }
    
    // scaling factors
    Double_t energy                 = ReturnCollisionEnergy( optionEnergy);
    Double_t doubleAddFactorK0s     = ReturnCorrectK0ScalingFactor( optionEnergy,  fEventCutSelection);
    
    // Set flags for MC case
    if (isMC.CompareTo("kTRUE") ==0){ 
        prefix2             = "MC";
        doubleAddFactorK0s  = 0.;
        kIsMC               = kTRUE;
    } else { 
        prefix2             = "data";
    }
    
    // Set flag for pi0/eta
    Bool_t kIsEta           = kFALSE;
    if (!nameMeson.Contains("Pi0") ) 
        kIsEta              = kTRUE;
    
    cout << "The additional K0 correction factor is: "  << doubleAddFactorK0s<<endl;
    TString textMeson=ReturnMesonString ( nameMeson);
    if (textMeson.CompareTo("") == 0) return;

    // flags for collisions sytem
    if (collisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;  
    }
    Int_t kCollisionSystem = 0; // 0 : pp, 1: PbPb, 2: pPb
    if (optionEnergy.CompareTo("PbPb_2.76TeV") == 0) kCollisionSystem = 1;
    if (optionEnergy.CompareTo("pPb_5.023TeV") == 0) kCollisionSystem = 2;
    
    TString centralityString    = GetCentralityString(fEventCutSelection);
    TString centralityString2   = GetCentralityString(fEventCutSelection);
    TString fTextCent       ="";
    if (centralityString.CompareTo("pp")==0){
        fTextCent           = "MinBias";  
        centralityString    = "";
    } else {
        fTextCent = Form("%s central", centralityString.Data());
        if ( !centralityString.Contains("0-100%") )
            collisionSystem = Form("%s", collisionSystem.Data());
    }
    if (optionPeriod.CompareTo("") != 0 ){
        collisionSystem     = Form("%s, %s",collisionSystem.Data(),optionPeriod.Data());
    }   
        
    TString rapidityRange   = "";
    Double_t deltaRapid     = ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
    
    TString trigger         = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
    
    // Initialize bin for single invariant mass plot
    Int_t fExampleBin       = 2;
    Double_t scaleFacSinBin = 1.0;
    fExampleBin             = ReturnSingleInvariantMassBinPlotting (nameMeson, optionEnergy, mode, trigger.Atoi(), scaleFacSinBin, triggerSet);
    
    //Variable defintion
    Double_t scaling        = 1./(2.*TMath::Pi());

    // Naming conventions
    TString nameSecMeson[4]                         = {"K0S", "Lambda", "K0L", "Rest"};
    TString nameSecMesonPlot[4]                     = {"K_{s}^{0}", "#Lambda", "K_{l}^{0}", "Rest"};
    TString nameResMeson[8]                         = {"Eta", "Rho0", "EtaPrime", "Omega", "Phi", "Delta0", "Delta+", "JPsi"};
    TString nameResMesonPlot[8]                     = {"#eta", "#rho^{0}", "#eta'", "#omega", "#phi", "#Delta^{0}", "#Delta^{+}", "J/#psi"};
    TString nameIntRange[6]                         = {"", "Wide", "Narrow", "Left", "LeftWide", "LeftNarrow"};
    TString nameIntRangePlot[6]                     = {"right/ normal int", "right/wide int", "right/narrow int", "left/ normal int", "left/wide int", "left/narrow int"};
    TString nameIntBckResult[3]                     = {"pol2_normal","exp_normal","exp2_normal"};

    Int_t pdgSecMeson[3]                            = {310, 3122, 130};
    Double_t decayLength[3]                         = {0.026844, 0.0789, 15.34};
    Double_t stacklength                            = 4.5;
    
    // Color and style setting for plotting
    Color_t colorIntRanges[6]                       = {kBlack, kGray+3, kGray+1, kBlue, kBlue+2, kBlue-5};
    Color_t colorIntRangesSec[6]                    = {kRed+1, kRed+2, kRed-5, kBlue, kBlue+2, kBlue-5};
    Style_t markerStyleIntRanges[6]                 = {20, 24, 24, 20, 24, 24};
    Size_t markerSizeIntRanges[6]                   = {1, 1, 1, 1, 1, 1};
    Color_t colorCat[6]                             = { kRed+1, 807, 800, kGreen+2, kCyan+2, kBlue+1};
    Color_t colorCatMC[6]                           = { kRed+3, 807+2, 800+2, kGreen+4, kCyan+4, kBlue+3};
    Style_t styleCat[6]                             = { 20, 21, 29, 33, 20, 21};
    Style_t styleCatMC[6]                           = { 24, 25, 30, 27, 24, 25};
    Color_t colorSec[4]                             = {kRed+2, 807, kCyan+2, kBlue};
    Style_t markerStyleSec[4]                       = {21, 33, 29, 34};
    Style_t markerStyleSecWithToy[4]                = {25, 27, 30, 28};
    Size_t markerSizeSec[4]                         = {1.5, 1.75, 2., 1.5};
    Color_t colorSecFromToy[3]                      = {kRed-2, kOrange+1, kCyan-2};
    Style_t markerStyleSecFromToy[3]                = {20, 33, 29};
    Size_t markerSizeSecFromToy[3]                  = {1.7, 2, 2.2};
    Style_t markerStyleFeedDown[8]                  = {20, 25, 34, 27, 24, 21, 28, 33};
    Size_t  markerSizeFeedDown[8]                   = {1.7, 1.7, 2.2, 2.2, 1.7, 1.7, 2.2, 2.2};
    Color_t colorFeedDown[8]                        = {kRed-2, kOrange+1, kCyan-2, kBlue-2, kRed-2, kOrange+1, kCyan-2, kBlue-2};
    

    // set correct meson mass
    Double_t mesonMassExpect    = 0;
    if( !kIsEta )    
        mesonMassExpect         = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    if( kIsEta )     
        mesonMassExpect         = TDatabasePDG::Instance()->GetParticle(221)->Mass();
    
    
    //*******************************************************************************************************
    //***********************************Reading data file **************************************************
    //*******************************************************************************************************
    // File definitions
    TFile fileUncorrected(fileNameUnCorrectedFile.Data());
    if (fileUncorrected.IsZombie()) return;
    TH1F *histoNumberOfGoodESDTracksVtx             = (TH1F*)fileUncorrected.Get("GoodESDTracks");
    TH1D *histoEventQuality                         = (TH1D*)fileUncorrected.Get("NEvents");
    TH1D *histoUnCorrectedYield[6]                  = {NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D *RatioRaw[6]                               = {NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D *fHistoYieldDiffBckResult[3]               = {NULL, NULL, NULL};
    for (Int_t k= 0; k < 6; k++){
        histoUnCorrectedYield[k]                    = (TH1D*)fileUncorrected.Get(Form("histoYieldMeson%s", nameIntRange[k].Data()));
        RatioRaw[k]                                 = (TH1D*) histoUnCorrectedYield[k]->Clone(Form("RatioRaw%s_%s",nameIntRange[k].Data(),nameIntRange[0].Data())); 
        RatioRaw[k]->Divide(RatioRaw[k],histoUnCorrectedYield[0],1.,1.,"B");
    }
    for (Int_t k= 0; k < 3; k++){
        fHistoYieldDiffBckResult[k]                 = (TH1D*)fileUncorrected.Get(Form("histoYieldMesonDiffBckResult%s", nameIntBckResult[k].Data()));
        fHistoYieldDiffBckResult[k]->Scale(100.);
    }
    TH1D *histoFWHMMeson                            = (TH1D*)fileUncorrected.Get("histoFWHMMeson"); 
    TH1D *histoMassMeson                            = (TH1D*)fileUncorrected.Get("histoMassMeson");
    TH1D *histoWidthGaussianMeson                   = (TH1D*)fileUncorrected.Get("histoWidthGaussianMeson"); 
    TH1D *histoMassGaussianMeson                    = (TH1D*)fileUncorrected.Get("histoMassGaussianMeson");
    TH1D* histoInvMassSignalPlusBG                  = (TH1D*)fileUncorrected.Get(Form("Mapping_GG_InvMass_in_Pt_Bin%02d",fExampleBin));
    TH1D* histoInvMassBG                            = (TH1D*)fileUncorrected.Get(Form("Mapping_BckNorm_InvMass_in_Pt_Bin%02d",fExampleBin));
    TH1D* histoInvMassSignal                        = (TH1D*)fileUncorrected.Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",fExampleBin));
    TF1* fitInvMassSignal                           = (TF1*)fileUncorrected.Get(Form("Signal_InvMassFit_in_Pt_Bin%02d",fExampleBin));
    TH1D *deltaPt                                   = (TH1D*)fileUncorrected.Get("deltaPt");
    for (Int_t i = 0; i < deltaPt->GetNbinsX() +1; i++){
        deltaPt->SetBinError(i, 0);
    }
    // calculate number of events
    Float_t nEvt    = 0;
    if (kCollisionSystem == 1){
        nEvt        = histoEventQuality->GetBinContent(1);
    } else {
        nEvt        = GetNEvents(histoEventQuality);
        // BinContent 5 - Zvertex-position, BinContent 4 - no Trigger Bit, BinContent 7 - PileUp 
    }

    Double_t scaleFactorMeasXSecForExternalInput              = 1;
    if ( kCollisionSystem == 0){
        // obtain effective xSection
        Int_t isV0AND           = 0;
        if (histoEventQuality->GetNbinsX() > 7){
            if (histoEventQuality->GetBinContent(9) > 0){
                isV0AND         = 1;
            }
        }
        if (optionEnergy.CompareTo("8TeV") == 0){
            isV0AND             = 1;
        }
        if (trigger.Atoi() == 10 || trigger.Atoi() == 52 || trigger.Atoi() == 83  || trigger.Atoi() == 85 || trigger.Atoi() == 81 ){
            isV0AND             = 1;
        }

        Double_t xSectionEff    = ReturnCorrectXSection( optionEnergy, isV0AND);
        Double_t xSectionINEL   = ReturnCorrectXSection( optionEnergy, 3);
        if (xSectionINEL != 0){
            scaleFactorMeasXSecForExternalInput               = xSectionINEL/xSectionEff; // was: xSectionEff/xSectionINEL;
        }
    }
    // The trigger rejection factor has to be given externally 
    Double_t triggerRejection                       = ReturnTriggerRejectionFactor(optionEnergy, trigger.Atoi());
    cout << "trigger rejection factor set to: " << triggerRejection << endl;
    if (triggerRejection != 1.)
        scaleFactorMeasXSecForExternalInput                   = scaleFactorMeasXSecForExternalInput*triggerRejection;
    if (scaleFactorMeasXSecForExternalInput != 1)
        cout << "The secondary correction from the toy has to be scaled with " << scaleFactorMeasXSecForExternalInput << endl;

    // read cocktail input if available
    TString strExternalInputName                            = "";
    TH1D* histoExternalInputSecPi0[3]                       = {NULL, NULL, NULL};
    TH1D* histoExternalInputFeedDownPi0[8]                  = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    Bool_t foundCocktailInput                               = kFALSE;
    if (!kIsMC){
        for (Int_t j = 0; j < 3; j++){
            histoExternalInputSecPi0[j]                     = (TH1D*)fileUncorrected.Get(Form("histoSecPi0YieldFrom%s_FromCocktail",nameSecMeson[j].Data()));
            if (histoExternalInputSecPi0[j]){
                foundCocktailInput                          = kTRUE;
                cout << "Using the cocktail input for secondary correction" << endl;
                if (j==0) strExternalInputName              = "Cocktail";
            }
        }
        
        // contribution from resonance feed down
        for (Int_t j = 0; j < 8; j++) {
            histoExternalInputFeedDownPi0[j]                = (TH1D*)fileUncorrected.Get(Form("histoResonanceFeedDownPi0YieldFrom%s_FromCocktail",nameResMeson[j].Data()));
        }
    }
    // read toy MC input if available
    Bool_t foundToyMCInput                          = kFALSE;
    // only use toyMC if cocktail input is not available
    if (!kIsMC && !foundCocktailInput){
        for (Int_t j = 0; j < 3; j++){
            histoExternalInputSecPi0[j]                    = (TH1D*)fileUncorrected.Get(Form("histoSecPi0YieldFrom%s_FromToy",nameSecMeson[j].Data()));
            if (histoExternalInputSecPi0[j]){
                histoExternalInputSecPi0[j]->Scale(scaleFactorMeasXSecForExternalInput);
                foundToyMCInput                         = kTRUE;
                cout << "Using the toyMC input for secondary correction" << endl;
                if (j==0) strExternalInputName = "Toy";
            }
        }
    }

    // set min and max values for pT
    Double_t maxPtMeson     = histoUnCorrectedYield[0]->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield[0]->GetNbinsX());
    Double_t minPtMeson     = 0;
    Int_t ptBin             = 1;
    while (histoUnCorrectedYield[0]->GetBinContent(ptBin) == 0.){
        ptBin++;
        minPtMeson          = histoUnCorrectedYield[0]->GetXaxis()->GetBinLowEdge(ptBin);
    }

    Double_t minPtMesonSec  = 0.3;
    if (mode == 0 && kCollisionSystem==1)
        minPtMesonSec       = minPtMeson;
    else if (mode == 2)
        minPtMesonSec       = minPtMeson;
    else if (mode == 4)
        minPtMesonSec       = minPtMeson;
    cout << "minimum pT: " << minPtMeson << ", maximum pT: " << maxPtMeson << endl;


    //*******************************************************************************************************
    //***********************************Reading MC correction file *****************************************
    //*******************************************************************************************************    
    TFile* fileCorrections =         new TFile(fileNameCorrectionFile.Data());
    if (fileCorrections->IsZombie()) return;
    TH1F *histoEventQualityMC                       = (TH1F*)fileCorrections->Get("NEvents");
    TH1D* histoEffiPt[6]                            = {NULL, NULL, NULL, NULL, NULL, NULL};
    for (Int_t k= 0; k < 6; k++){
        histoEffiPt[k]                              = (TH1D*)fileCorrections->Get(Form("Meson%sEffiPt", nameIntRange[k].Data())); //not yet correct MesonEffiPt
    }
    if (kIsMC){
        for (Int_t j = 0; j < 3; j++){
            histoExternalInputSecPi0[j]                    = (TH1D*)fileCorrections->Get(Form("MCSecPi0From%s_Rebinned",nameSecMeson[j].Data()));
            if (histoExternalInputSecPi0[j]){
                histoExternalInputSecPi0[j]->Scale(1/nEvt);
                foundToyMCInput                         = kTRUE;
            }
        }
    }

    TH1D* histoAcceptance                           = (TH1D*)fileCorrections->Get("fMCMesonAccepPt");
    TH1D* histoSecAcceptance[4]                     = {NULL, NULL, NULL, NULL};
    Bool_t hasNewSecQuantities                      = kFALSE;
    Int_t nAccHistSec                               = 0;
    Bool_t modifiedSecAcc[4]                        = {kFALSE, kFALSE, kFALSE, kFALSE};
    for (Int_t j = 0; j< 4; j++){
        histoSecAcceptance[j]                       = (TH1D*)fileCorrections->Get(Form("fMCSecPi0From%sAccepPt",nameSecMeson[j].Data()));
        if (histoSecAcceptance[j]){
            if ( mode == 4 &&  ( j == 1 || j == 2) ){
                Double_t    accSec                  = ReturnDeltaEtaCalo(fClusterCutSelection)/deltaRapid*ReturnDeltaPhiCalo(fClusterCutSelection)/(2*TMath::Pi());
                for (Int_t iPt = histoSecAcceptance[j]->FindBin(minPtMeson); iPt< histoSecAcceptance[j]->GetNbinsX()+1; iPt++ ){
                    histoSecAcceptance[j]->SetBinContent(iPt, accSec);
                    histoSecAcceptance[j]->SetBinError(iPt, accSec*0.1);
                }    
                modifiedSecAcc[j]                   = kTRUE;
            }    
            hasNewSecQuantities                     = kTRUE;
            nAccHistSec++;
        }
    }

    TH1D *histoAcceptanceWOEvtWeights               = (TH1D*)fileCorrections->Get("fMCMesonAccepPtWOEvtWeights");

    // load histograms without weighting for comparison of efficiencies & enable the comparison if those are present in the input file
    Bool_t containsWOWeights                        = kFALSE;
    TH1D* histoAcceptanceWOWeights                  = NULL;
    histoAcceptanceWOWeights                        = (TH1D*)fileCorrections->Get("fMCMesonAccepPtWOWeights"); 
    TH1D* histoTrueEffiPtWOWeights[3]               = {NULL, NULL, NULL};
    TH1D* histoTrueEffiPt[3]                        = {NULL, NULL, NULL};
    TH1D* histoTrueEffiPtFixed[3]                   = {NULL, NULL, NULL};
    TH1D* histoEffiPtFixed[3]                       = {NULL, NULL, NULL};

    for (Int_t k = 0; k < 3; k++){
        histoTrueEffiPtWOWeights[k]                 = (TH1D*)fileCorrections->Get(Form("TrueMeson%sEffiPtUnweighted", nameIntRange[k].Data())); 
        histoTrueEffiPt[k]                          = (TH1D*)fileCorrections->Get(Form("TrueMeson%sEffiPt", nameIntRange[k].Data())); 
        histoTrueEffiPtFixed[k]                     = (TH1D*)histoTrueEffiPt[k]->Clone(Form("TrueMeson%sEffiPtFixed", nameIntRange[k].Data())); 
        histoTrueEffiPtFixed[k]                     = FixEfficiency(histoTrueEffiPtFixed[k],histoTrueEffiPt[k],optionEnergy,centralityString);
        histoEffiPtFixed[k]                         = (TH1D*)histoEffiPt[k]->Clone(Form("Meson%sEffiPtFixed", nameIntRange[k].Data())); 
        histoEffiPtFixed[k]                         = FixEfficiency(histoEffiPtFixed[k],histoEffiPt[k],optionEnergy,centralityString);
    }    
    // read secondary pi0 efficiencies
    TH1D* histoSecTrueEffi[3][4]                    = {{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL}};
    TH1D* histoSecTrueEffiFromFit[3][4]             = {{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL}};
    TH1D* histoRatioSecEffDivTrueEff[3][4]          = {{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL}};
    TF1*  fithistoRatioSecEffDivTrueEff[3][4]       = {{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL},{NULL, NULL, NULL, NULL}};

    Bool_t modifiedSecTrueEffi[3][4]                = {{kFALSE, kFALSE, kFALSE, kFALSE}, {kFALSE, kFALSE, kFALSE, kFALSE}, {kFALSE, kFALSE, kFALSE, kFALSE}};
    Int_t nEffHistSec                               = 0;
    for (Int_t j = 0; j< 4; j++){
        for (Int_t k = 0; k< 3; k++){
            histoSecTrueEffi[k][j]                  = (TH1D*)fileCorrections->Get(Form("TrueSecFrom%s%sEffiPt",nameSecMeson[j].Data(), nameIntRange[k].Data()));
            if (histoSecTrueEffi[k][j]){
                histoRatioSecEffDivTrueEff[k][j]    = (TH1D*)histoSecTrueEffi[k][j]->Clone(Form("ratioSecEffDivTrueEff%s%s",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                histoRatioSecEffDivTrueEff[k][j]->Divide(histoRatioSecEffDivTrueEff[k][j],histoTrueEffiPt[k]);
                // fit the K0s efficieny ratio with an exponential if cocktail or toyMC input is used
                TF1*  fitConst                      = new TF1("fitConst","[0]");
                fitConst->SetLineColor(colorSec[j]);
                histoRatioSecEffDivTrueEff[k][j]->Fit(fitConst);
                cout << fitConst->GetParameter(0) << "\t +-" << fitConst->GetParError(0) << endl;

                if((j==0) && (foundCocktailInput || foundToyMCInput)){
//                    fithistoRatioSecEffDivTrueEff[k][j] = FitObject("h",Form("fitexpEffi%s_%s",nameSecMeson[j].Data(),nameIntRange[k].Data()),"Pi0",histoRatioSecEffDivTrueEff[k][0],minPtMesonSec,maxPtMeson,NULL,"QNRME+");
                    if(kCollisionSystem==1) histoRatioSecEffDivTrueEff[k][j]->Fit(fitConst,"QNRME+","",4,maxPtMeson);
                    else histoRatioSecEffDivTrueEff[k][j]->Fit(fitConst,"QNRME+","",2.5,maxPtMeson);
                    
                    fithistoRatioSecEffDivTrueEff[k][j] = new TF1(Form("fitexpEffi%s_%s",nameSecMeson[j].Data(),nameIntRange[k].Data()),"[0]/pow(x,[1])+[2]");
                    fithistoRatioSecEffDivTrueEff[k][j]->SetRange(minPtMesonSec,maxPtMeson);
                    fithistoRatioSecEffDivTrueEff[k][j]->FixParameter(2,fitConst->GetParameter(0));
                    histoRatioSecEffDivTrueEff[k][j]->Fit(fithistoRatioSecEffDivTrueEff[k][j],"QNRME+","",minPtMesonSec,maxPtMeson);
                }   

                
                cout << "rel stat. err sec effi: " << k << "\t"<< j << endl;
                Int_t nBinsActive                   = 0;
                Int_t nBinsTot                      = 0;
                for (Int_t iPt = 1; iPt< histoSecTrueEffi[k][j]->GetNbinsX()+1; iPt++){
                    if (histoSecTrueEffi[k][j]->GetBinCenter(iPt) > minPtMeson && histoSecTrueEffi[k][j]->GetBinContent(iPt) != 0){
                        nBinsActive++;
                        nBinsTot++;
                    }  else if (histoSecTrueEffi[k][j]->GetBinCenter(iPt) > minPtMeson){
                        nBinsTot++;
                    }
                }

                if (optionEnergy.CompareTo("2.76TeV") == 0){
                    histoSecTrueEffi[k][j]              = (TH1D*)histoTrueEffiPt[k]->Clone(Form("TrueSecFrom%s%sEffiPt",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                    // crude assumptions // need to be validated for other energies
                    if (mode == 4 ){
                        if (j == 0 ){
                            if ( fitConst->GetParameter(0) > 0.3 && fitConst->GetParameter(0) < 0.7)
                                histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                            else 
                                histoSecTrueEffi[k][j]->Scale(0.5);     
                        } else if ( (j == 1 || j == 2) ){
                            if ( (fitConst->GetParameter(0) > 0.7 && fitConst->GetParameter(0) < 1.3) && !(Double_t)nBinsActive/nBinsTot < 0.5 )
                               histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                            else     
                                histoSecTrueEffi[k][j]->Scale(1.);
                        } else if (j == 3) {
                            if ( (fitConst->GetParameter(0) > 0.6 && fitConst->GetParameter(0) < 1.2) && !(Double_t)nBinsActive/nBinsTot < 0.5 )
                                histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                            else     
                                histoSecTrueEffi[k][j]->Scale(0.75);   
                        }    
                    } else if (mode == 2 ) {
                        cout << "entered here" << endl;
                        if (j == 0){
                            if ( fitConst->GetParameter(0) > 0.2 && fitConst->GetParameter(0) < 0.3)
                                histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                            else 
                                histoSecTrueEffi[k][j]->Scale(0.25);     
                        } else if (j == 1 || j == 2){
                            if ( (fitConst->GetParameter(0) > 0. && fitConst->GetParameter(0) < 0.15) && !(Double_t)nBinsActive/nBinsTot < 0.5 )
                                histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                            else 
                                histoSecTrueEffi[k][j]->Scale(0.05);
                        } else if (j == 3){
                            if ( (fitConst->GetParameter(0) > 0.1 && fitConst->GetParameter(0) < 0.3) && !(Double_t)nBinsActive/nBinsTot < 0.5 )
                                histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                            else 
                                histoSecTrueEffi[k][j]->Scale(0.15);   
                        }    
                    } else if  ( mode == 0) {
                        if (j == 0){
                            if ( fitConst->GetParameter(0) > 0.2 && fitConst->GetParameter(0) < 0.35)
                                histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                            else 
                                histoSecTrueEffi[k][j]->Scale(0.25);     
                        } else if (j == 1 || j == 2){
                            if ( (fitConst->GetParameter(0) > 0.1 && fitConst->GetParameter(0) < 0.3) && !(Double_t)nBinsActive/nBinsTot < 0.5 )
                                histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                            else 
                                histoSecTrueEffi[k][j]->Scale(0.2);
                        } else if (j == 3){
                            if ( (fitConst->GetParameter(0) > 0.0 && fitConst->GetParameter(0) < 0.15) && !(Double_t)nBinsActive/nBinsTot < 0.5 )
                                histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                            else 
                                histoSecTrueEffi[k][j]->Scale(0.1);
                        }    
                    }
                    if ( mode == 0 || mode == 2 || mode == 4){
                        modifiedSecTrueEffi[k][j]   = kTRUE;
                        cout << "adjusted sec effi, due to to little stat" << endl;
                    }
                } else if (optionEnergy.CompareTo("8TeV") == 0 || optionEnergy.CompareTo("7TeV") == 0 || optionEnergy.CompareTo("900GeV") == 0){
                    if (mode == 4 ){
                        modifiedSecTrueEffi[k][j]   = kTRUE;
                        if (j == 0 ){
                          if(optionEnergy.CompareTo("900GeV") == 0){
                            histoSecTrueEffi[k][j]              = (TH1D*)histoTrueEffiPt[k]->Clone(Form("TrueSecFrom%s%sEffiPt",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                            cout << Form("SECONDARIES: Fixed %s ",nameSecMeson[j].Data()) << "efficiency" << endl;
                            histoSecTrueEffi[k][j]->Scale(0.5);
                          }else{
                            histoSecTrueEffi[k][j]              = (TH1D*)histoTrueEffiPt[k]->Clone(Form("TrueSecFrom%s%sEffiPt",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                            cout << Form("SECONDARIES: Calculated the %s ",nameSecMeson[j].Data()) << "efficiency from the fit" << endl;
                            histoSecTrueEffi[k][j] ->Multiply(fithistoRatioSecEffDivTrueEff[k][j]);
                          }
                        } else if ( j == 1 ){
                            histoSecTrueEffi[k][j]              = (TH1D*)histoTrueEffiPt[k]->Clone(Form("TrueSecFrom%s%sEffiPt",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                            cout << Form("SECONDARIES: Fixed %s ",nameSecMeson[j].Data()) << "efficiency" << endl;
                            histoSecTrueEffi[k][j]->Scale(0.5);
                        } else if ( j == 2 ){
                            histoSecTrueEffi[k][j]              = (TH1D*)histoTrueEffiPt[k]->Clone(Form("TrueSecFrom%s%sEffiPt",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                            cout << Form("SECONDARIES: Fixed %s ",nameSecMeson[j].Data()) << "efficiency" << endl;
                            histoSecTrueEffi[k][j]->Scale(1.0);
                        } else if ( j == 3 ){
                            histoSecTrueEffi[k][j]              = (TH1D*)histoTrueEffiPt[k]->Clone(Form("TrueSecFrom%s%sEffiPt",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                            cout << Form("SECONDARIES: Calculated the %s ",nameSecMeson[j].Data()) << "efficiency from the fit" << endl;
                            histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                        }
                    } else if (mode == 2 ) {
                        modifiedSecTrueEffi[k][j]   = kTRUE;
                        if (j == 0 ){
                          if(optionEnergy.CompareTo("900GeV") == 0){
                            histoSecTrueEffi[k][j]              = (TH1D*)histoTrueEffiPt[k]->Clone(Form("TrueSecFrom%s%sEffiPt",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                            cout << Form("SECONDARIES: Fixed %s ",nameSecMeson[j].Data()) << "efficiency" << endl;
                            histoSecTrueEffi[k][j]->Scale(0.5);
                          }else{
                            histoSecTrueEffi[k][j]              = (TH1D*)histoTrueEffiPt[k]->Clone(Form("TrueSecFrom%s%sEffiPt",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                            cout << Form("SECONDARIES: Calculated the %s ",nameSecMeson[j].Data()) << "efficiency from the fit" << endl;
                            histoSecTrueEffi[k][j] ->Multiply(fithistoRatioSecEffDivTrueEff[k][j]);
                          }
                        } else if ( (j == 1 || j == 2) ){
                            histoSecTrueEffi[k][j]              = (TH1D*)histoTrueEffiPt[k]->Clone(Form("TrueSecFrom%s%sEffiPt",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                            cout << Form("SECONDARIES: Fixed %s ",nameSecMeson[j].Data()) << "efficiency" << endl;
                            histoSecTrueEffi[k][j]->Scale(0.1);
                        } else if ( j == 3 ){
                            histoSecTrueEffi[k][j]              = (TH1D*)histoTrueEffiPt[k]->Clone(Form("TrueSecFrom%s%sEffiPt",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                            cout << Form("SECONDARIES: Calculated the %s ",nameSecMeson[j].Data()) << "efficiency from the fit" << endl;
                            histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                        }
                    }
                } else if (optionEnergy.CompareTo("pPb_5.023TeV") == 0){
                    histoSecTrueEffi[k][j]              = (TH1D*)histoTrueEffiPt[k]->Clone(Form("TrueSecFrom%s%sEffiPt",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                    // crude assumptions // need to be validated for other energies
//                     if (mode == 4 ){
//                         if (j == 0 ){
//                             if ( fitConst->GetParameter(0) > 0.3 && fitConst->GetParameter(0) < 0.7)
//                                 histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
//                             else 
//                                 histoSecTrueEffi[k][j]->Scale(0.5);     
//                         } else if ( (j == 1 || j == 2) ){
//                             if ( (fitConst->GetParameter(0) > 0.7 && fitConst->GetParameter(0) < 1.3) && !(Double_t)nBinsActive/nBinsTot < 0.5 )
//                                histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
//                             else     
//                                 histoSecTrueEffi[k][j]->Scale(1.);
//                         } else if (j == 3) {
//                             if ( (fitConst->GetParameter(0) > 0.6 && fitConst->GetParameter(0) < 1.2) && !(Double_t)nBinsActive/nBinsTot < 0.5 )
//                                 histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
//                             else     
//                                 histoSecTrueEffi[k][j]->Scale(0.75);   
//                         }    
//                     } else 
                    if (mode == 2 ) {
                        cout << "entered here" << endl;
                        if (j == 0){
                            if ( fitConst->GetParameter(0) > 0.2 && fitConst->GetParameter(0) < 0.3)
                                histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                            else 
                                histoSecTrueEffi[k][j]->Scale(0.25);     
                        } else if (j == 1 || j == 2){
                            if ( (fitConst->GetParameter(0) > 0. && fitConst->GetParameter(0) < 0.15) && !(Double_t)nBinsActive/nBinsTot < 0.5 )
                                histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                            else 
                                histoSecTrueEffi[k][j]->Scale(0.05);
                        } else if (j == 3){
                            if ( (fitConst->GetParameter(0) > 0.1 && fitConst->GetParameter(0) < 0.3) && !(Double_t)nBinsActive/nBinsTot < 0.5 )
                                histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                            else 
                                histoSecTrueEffi[k][j]->Scale(0.15);   
                        }    
                    } else if  ( mode == 0) {
                        if (j == 1 || j == 2){
                            if ( (fitConst->GetParameter(0) > 0.1 && fitConst->GetParameter(0) < 0.3) && !(Double_t)nBinsActive/nBinsTot < 0.5 )
                                histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                            else 
                                histoSecTrueEffi[k][j]->Scale(0.0);
                            modifiedSecTrueEffi[k][j]   = kTRUE;
                            
                        } else if (j == 3){
                            if ( (fitConst->GetParameter(0) > 0.0 && fitConst->GetParameter(0) < 0.15) && !(Double_t)nBinsActive/nBinsTot < 0.5 )
                                histoSecTrueEffi[k][j]->Scale(fitConst->GetParameter(0));
                            else 
                                histoSecTrueEffi[k][j]->Scale(0.1);
                            modifiedSecTrueEffi[k][j]   = kTRUE;
                        }    
                    }
                    if ( mode == 2 ){
                        modifiedSecTrueEffi[k][j]   = kTRUE;
                        cout << "adjusted sec effi, due to to little stat" << endl;
                    }
                }

                // use the fits from the MC efficiency ratio to get the secondary efficiencies for the raw yield calculation
                if (!modifiedSecTrueEffi[k][j] && (foundCocktailInput||foundToyMCInput)){
                    histoSecTrueEffi[k][j]              = (TH1D*)histoTrueEffiPt[k]->Clone(Form("TrueSecFrom%s%sEffiPt",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                    // for the K0s (j==0) use the exponential fit and for the other particles use the constant fit
                    if(j==0){
                        histoSecTrueEffi[k][j] ->Multiply(fithistoRatioSecEffDivTrueEff[k][j]);
                    } else {
                        histoSecTrueEffi[k][j] ->Scale(fitConst->GetParameter(0));
                    }
                    cout << Form("SECONDARIES: Calculated the %s ",nameSecMeson[j].Data()) << "efficiency from the fit" << endl;
                }
            }
        }

        if (histoSecTrueEffi[0][j]){
            nEffHistSec++;
        }
    }
//     return;
    
    if (histoAcceptanceWOWeights && histoTrueEffiPtWOWeights[0] && histoTrueEffiPtWOWeights[1] && histoTrueEffiPtWOWeights[2]) 
        containsWOWeights                           = kTRUE;
    

    TH1D* histoInputMesonPt                         = (TH1D*)fileCorrections->Get("MC_Meson_genPt");
    TH1D* histoInputMesonOldBinPt                   = (TH1D*)fileCorrections->Get("MC_Meson_genPt_oldBin");
    TH1D* histoInputMesonOldBinPtWOWeights          = NULL;
    TH1D* histoInputMesonOldBinPtWeights            = NULL;
    histoInputMesonOldBinPtWOWeights                = (TH1D*)fileCorrections->Get("MC_Meson_genPt_WOWeights");
    histoInputMesonOldBinPtWeights                  = (TH1D*)fileCorrections->Get("MC_Meson_genPt_Weights");
    TH1D* histoMCInputAddedSig                      = NULL;
    TH1D* histoMCInputWOWeightingAddedSig           = NULL;
    TH1D* histoMCInputWeightsAddedSig               = NULL;
    TH1F *histoEventQualityMCAddedSig               = NULL;
    histoMCInputAddedSig                            = (TH1D*)fileCorrections->Get("MC_Meson_genPt_oldBin_AddedSig");
    histoMCInputWOWeightingAddedSig                 = (TH1D*)fileCorrections->Get("MC_Meson_genPt_WOWeights_AddedSig");
    histoMCInputWeightsAddedSig                     = (TH1D*)fileCorrections->Get("MC_Meson_genPt_Weights_AddedSig");
    histoEventQualityMCAddedSig                     = (TH1F*)fileCorrections->Get("NEvents_AddedSig");

    TH1D* histoMCInputJetJetMC                      = NULL;
    TH1F *histoEventQualityMCJetJet                 = NULL;
    histoMCInputJetJetMC                            = (TH1D*)fileCorrections->Get("MC_Meson_genPt_oldBin_JetJetMC");
    histoEventQualityMCJetJet                       = (TH1F*)fileCorrections->Get("NEvents_JetJetMC");
        
    TH1D* histoTrueMassMeson                        = (TH1D*)fileCorrections->Get("histoTrueMassMeson");
    TH1D* histoTrueFWHMMeson                        = (TH1D*)fileCorrections->Get("histoTrueFWHMMeson");
    TH1D* histoTrueMassGaussianMeson                = (TH1D*)fileCorrections->Get("histoTrueMassGaussianMeson");
    TH1D* histoTrueWidthGaussianMeson               = (TH1D*)fileCorrections->Get("histoTrueWidthGaussianMeson");
    TH1D* histoMCrecMassMeson                       = (TH1D*)fileCorrections->Get("histoMassMeson");
    TH1D* histoMCrecMassGaussianMeson               = (TH1D*)fileCorrections->Get("histoMassGaussianMeson");
    TH1D* histoMCrecFWHMMeson                       = (TH1D*)fileCorrections->Get("histoFWHMMeson");
    TH1D* histoMCrecWidthGaussMeson                 = (TH1D*)fileCorrections->Get("histoWidthGaussianMeson");
    
    Float_t nEvtMC = 0;
    if (kCollisionSystem == 1){
        nEvtMC = histoEventQualityMC->GetBinContent(1);
    } else {
        nEvtMC = GetNEvents(histoEventQualityMC);
        // BinContent 5 - Zvertex-position, BinContent 4 - no Trigger Bit, BinContent 7 - PileUp 
    }
    Float_t nEvtMCAddSig = 0;
    if (histoEventQualityMCAddedSig){
        if (kCollisionSystem == 1){
            nEvtMCAddSig = histoEventQualityMCAddedSig->GetBinContent(1);
        } else {
            nEvtMCAddSig = GetNEvents(histoEventQualityMCAddedSig);
            // BinContent 5 - Zvertex-position, BinContent 4 - no Trigger Bit, BinContent 7 - PileUp 
        }
    }
    Float_t nEvtMCJetJet = 0;
    if (histoEventQualityMCJetJet){
        if (kCollisionSystem == 1){
            nEvtMCJetJet = histoEventQualityMCJetJet->GetBinContent(1);
        } else {
            nEvtMCJetJet = GetNEvents(histoEventQualityMCJetJet);
            // BinContent 5 - Zvertex-position, BinContent 4 - no Trigger Bit, BinContent 7 - PileUp 
        }
    }

    
    //*******************************************************************************************************
    //***********************************Read additional MC file for secondaries ****************************
    //*******************************************************************************************************
    TH1D* histoDefaultTrueSecFracMeson[3][4]        = {{ NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}};
    TF1* fitDefaultSecFrac[3][4]                    = {{ NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}};
    TFitResultPtr resultSecFrac[3][2];
    
    TFile* fileCorrectionsSecPi0 = new TFile("ExternalInput/PCM/SecondaryFractionHistogramms7TeV.root");
    if (fileCorrectionsSecPi0->IsZombie()) return;
    
    for (Int_t k = 0; k < 3; k++){
        histoDefaultTrueSecFracMeson[k][0]          = (TH1D*)fileCorrectionsSecPi0->Get(Form("TrueSecFracFromK0S%s", nameIntRange[k].Data()));
        histoDefaultTrueSecFracMeson[k][1]          = (TH1D*)histoDefaultTrueSecFracMeson[k][0]->Clone(Form("TrueSecFracFromLambda%s", nameIntRange[k].Data()));
        histoDefaultTrueSecFracMeson[k][2]          = (TH1D*)histoDefaultTrueSecFracMeson[k][0]->Clone(Form("TrueSecFracFromK0L%s", nameIntRange[k].Data()));
        
        // set bins for lambda and K0L to 0
        for (Int_t i = 1; i < histoDefaultTrueSecFracMeson[k][0]->GetNbinsX()+1; i++){
            histoDefaultTrueSecFracMeson[k][1]->SetBinContent(i,0);
            histoDefaultTrueSecFracMeson[k][1]->SetBinError(i,0);
            histoDefaultTrueSecFracMeson[k][2]->SetBinContent(i,0);
            histoDefaultTrueSecFracMeson[k][2]->SetBinError(i,0);
        }    
        
        // take out the K0s from the full histo
        histoDefaultTrueSecFracMeson[k][3]          = (TH1D*)fileCorrectionsSecPi0->Get(Form("TrueSecFrac%s", nameIntRange[k].Data()));
        histoDefaultTrueSecFracMeson[k][3]->Sumw2();
        histoDefaultTrueSecFracMeson[k][3]->Add(histoDefaultTrueSecFracMeson[k][0],-1);
    }    
    
    for (Int_t k = 0; k < 3; k++){
        cout << "fitting default hist for "<< nameIntRange[k].Data() << " K0s" << endl;
        fitDefaultSecFrac[k][0]                     = new TF1(Form("fitDefaultSecFromK0SFrac%s",nameIntRange[k].Data()),"[0]/pow(x,[1])");
        fitDefaultSecFrac[k][0]->SetRange(minPtMesonSec, maxPtMeson);
        resultSecFrac[k][0]                         = histoDefaultTrueSecFracMeson[k][0]->Fit(fitDefaultSecFrac[k][0],"SNRME+","",minPtMesonSec, maxPtMeson);
        cout << "fitting default hist for "<< nameIntRange[k].Data() << " Rest" << endl;
        fitDefaultSecFrac[k][3]                     = new TF1(Form("fitDefaultSecFrac%s",nameIntRange[k].Data()),"[0]/pow(x,[1])");
        fitDefaultSecFrac[k][3]->SetRange(minPtMesonSec, maxPtMeson);
        resultSecFrac[k][1]                         = histoDefaultTrueSecFracMeson[k][3]->Fit(fitDefaultSecFrac[k][3],"SNRME+","",minPtMesonSec, maxPtMeson);
    }
    cout << "survived fitting" << endl;

    //*******************************************************************************************************
    //******************** Read secondary histograms from current MC corr file ******************************
    //*******************************************************************************************************
    TH1D *histoYieldTrueSecFracMeson[3][4]          = { { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}};
    TH1D *histoYieldTrueSecFracMeson_orig[3][4]     = { { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}};
    TF1* fitSecFracPurePowerlaw[3][4]               = { { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}};
    TF1* fitSecFracPLWithConst[3][4]                = { { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}};
    for (Int_t k = 0; k < 3; k++){
        for (Int_t j = 0; j < 4; j++){
            fitSecFracPurePowerlaw[k][j]            = new TF1( Form("fitSecFracPurePowerlawFrom%s%s",nameSecMeson[k].Data(),nameIntRange[k].Data()) ,"[0]/pow(x,[1])");
            fitSecFracPLWithConst[k][j]             = new TF1( Form("fitSecFracPLWithConstFrom%s%s",nameSecMeson[k].Data(),nameIntRange[k].Data()) ,"[0]/pow(x,[1])+[2]");
            if ( optionEnergy.CompareTo("8TeV") == 0 ){
               if ( mode != 2 && mode != 4 ) fitSecFracPLWithConst[k][j]->SetParLimits(1,0,1000);
            }else if ( (mode == 2 || mode == 4) ) fitSecFracPLWithConst[k][j]->SetParLimits(1,0,1000);
        }    
    }
    TH1D* histoYieldSecMeson[6][4]                  = { { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL},
                                                        { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL} };
    TH1D* histoRatioYieldSecMeson[6][4]             = { { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL},
                                                        { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL} };
    TH1D* histoYieldSecMesonFromExternalInput[6][3] = { { NULL, NULL, NULL}, { NULL, NULL, NULL}, { NULL, NULL, NULL},
                                                        { NULL, NULL, NULL}, { NULL, NULL, NULL}, { NULL, NULL, NULL} };
    TH1D* histoRatioYieldSecMesonFromExtInput[6][3] = { { NULL, NULL, NULL}, { NULL, NULL, NULL}, { NULL, NULL, NULL},
                                                        { NULL, NULL, NULL}, { NULL, NULL, NULL}, { NULL, NULL, NULL} };
    TH1D* histoYieldResonanceFeedDownPi0FromExternalInput[6][8] = { {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                                    {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                                    {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                                    {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                                    {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                                    {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}  };
    TH1D* histoRatioYieldResonanceFeedDownPi0FromExternalInput[6][8] = {    {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                                            {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                                            {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                                            {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                                            {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                                            {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}  };
    TH1D* histoYieldTrueGGFracMeson[3]              = { NULL, NULL, NULL };
    TH1D* histoYieldGGMeson[6]                      = { NULL, NULL, NULL, NULL, NULL, NULL };
    Bool_t haveSec[4]                               = { kTRUE, kTRUE, kTRUE, kTRUE };
    Bool_t haveSecUsed[4]                           = { kTRUE, kTRUE, kTRUE, kTRUE };
    Double_t scalingFacSec[4]                       = { 1+doubleAddFactorK0s, 1, 1, 1};
    
    Bool_t doK0SecCorrection                        = kFALSE;
    Int_t doK0SecCorrectionWithDefaultHisto         = 0;
    Bool_t doGGCorrection                           = kFALSE;
    
    if ( !kIsEta ) {
        doK0SecCorrection = kTRUE;
        if ( (mode == 0 || mode == 9) && !kIsMC && kCollisionSystem != 2 )                  doK0SecCorrectionWithDefaultHisto = 1;
        if ( (mode == 0 || mode == 9 || mode == 1) && !kIsMC && kCollisionSystem == 2 )     doK0SecCorrectionWithDefaultHisto = 2;
    }
    if (optDalitz) {
        doGGCorrection                          = kTRUE;
        doK0SecCorrection                       = kFALSE;
        doK0SecCorrectionWithDefaultHisto       = kFALSE;
    }    

    if (doK0SecCorrection){
        // create secondary yield for different integration windows: normal, wide, narrow
        for (Int_t k = 0; k < 3; k++){
            // initialize histos for this integration range for different secondary components: K0S, Lambda, K0L, Rest
            for (Int_t j = 0; j < 4; j++){
                histoYieldTrueSecFracMeson[k][j]            = (TH1D*)fileCorrections->Get(Form("TrueSecFrom%sFrac%s",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                if (histoYieldTrueSecFracMeson[k][j]){
                    histoYieldTrueSecFracMeson_orig[k][j]   = (TH1D*)histoYieldTrueSecFracMeson[k][j]->Clone(Form("TrueSecFrom%sFrac%s_orig",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                    // set fractions to 0 if none of the pt bins is above 1e-3 => 0.1%
                    if (FindLargestBin1DHist(histoYieldTrueSecFracMeson[k][j]) < 1e-3 ){
                        haveSecUsed[j]                      = kFALSE;
                        for (Int_t i = 1; i < histoYieldTrueSecFracMeson[k][j]->GetNbinsX()+1; i++){
                            histoYieldTrueSecFracMeson[k][j]->SetBinContent(i, 0);
                            histoYieldTrueSecFracMeson[k][j]->SetBinError(i, 0);                            
                        }
                    }
                    // exception for PCM mode set K0L to 0
                    if ( mode == 0 && j == 2){
                        haveSecUsed[j]                      = kFALSE;
                        for (Int_t i = 1; i < histoYieldTrueSecFracMeson[k][j]->GetNbinsX()+1; i++){
                            histoYieldTrueSecFracMeson[k][j]->SetBinContent(i, 0);
                            histoYieldTrueSecFracMeson[k][j]->SetBinError(i, 0);                            
                        }
                        histoYieldTrueSecFracMeson_orig[k][j]   = (TH1D*)histoYieldTrueSecFracMeson_orig[0][0]->Clone(Form("TrueSecFrom%sFrac%s_orig",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                    }    
                } else {
                    haveSec[j]                              = kFALSE;
                    haveSecUsed[j]                          = kFALSE;
                    histoYieldTrueSecFracMeson[k][j]        = (TH1D*)histoYieldTrueSecFracMeson[0][0]->Clone(Form("TrueSecFrom%sFrac%s",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                    histoYieldTrueSecFracMeson_orig[k][j]   = (TH1D*)histoYieldTrueSecFracMeson_orig[0][0]->Clone(Form("TrueSecFrom%sFrac%s_orig",nameSecMeson[j].Data(), nameIntRange[k].Data()));
                }
            }
        }
        // modify respective histos
        for (Int_t k = 0; k < 3; k++){
            // this is a really really ugly hack!
            if (doK0SecCorrectionWithDefaultHisto == 1){
                for (Int_t i = 1; i < histoYieldTrueSecFracMeson[k][0]->GetNbinsX()+1; i++){

                    Double_t ptStart        = histoYieldTrueSecFracMeson[k][3]->GetXaxis()->GetBinLowEdge(i);
                    Double_t ptEnd          = histoYieldTrueSecFracMeson[k][3]->GetXaxis()->GetBinUpEdge(i);
                    Double_t binWidth       = ptEnd-ptStart;
                    Double_t secFrac        = fitDefaultSecFrac[k][3]->Integral(ptStart, ptEnd, resultSecFrac[k][1]->GetParams()) / binWidth;
                    Double_t errorSecFrac   = fitDefaultSecFrac[k][3]->IntegralError(ptStart, ptEnd, resultSecFrac[k][1]->GetParams(), resultSecFrac[k][1]->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                    histoYieldTrueSecFracMeson[k][3]->SetBinContent(i, secFrac);
                    histoYieldTrueSecFracMeson[k][3]->SetBinError(i, errorSecFrac);

                    secFrac                 = fitDefaultSecFrac[k][0]->Integral(ptStart, ptEnd, resultSecFrac[k][0]->GetParams()) / binWidth;
                    errorSecFrac            = fitDefaultSecFrac[k][0]->IntegralError(ptStart, ptEnd, resultSecFrac[k][0]->GetParams(), resultSecFrac[k][0]->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                    histoYieldTrueSecFracMeson[k][0]->SetBinContent(i, secFrac);
                    histoYieldTrueSecFracMeson[k][0]->SetBinError(i, errorSecFrac);
                    histoYieldTrueSecFracMeson[k][1]->SetBinContent(i, 0);
                    histoYieldTrueSecFracMeson[k][1]->SetBinError(i, 0);
                    histoYieldTrueSecFracMeson[k][2]->SetBinContent(i, 0);
                    histoYieldTrueSecFracMeson[k][2]->SetBinError(i, 0);
                }
            // fit current histo with pure powerlaw fit (goes to 0 at high pt)    
            } else if ( doK0SecCorrectionWithDefaultHisto == 2 ){
                for (Int_t j = 0; j < 4; j++){
                    // exception for PCM mode set K0L to 0
                    if (haveSec[j] && haveSecUsed[j] && j != 2){
                        cout << "fitting current hist for "<< nameIntRange[k].Data() << " " << nameSecMeson[j].Data() << endl;
                        fitSecFracPurePowerlaw[k][j]->SetRange(minPtMesonSec, maxPtMeson);
                        TFitResultPtr resultCurr    = histoYieldTrueSecFracMeson[k][j]->Fit(fitSecFracPurePowerlaw[k][j],"SNRME+","",minPtMesonSec, maxPtMeson);
                        for (Int_t i = 1; i < histoYieldTrueSecFracMeson[k][j]->GetNbinsX()+1; i++){
                            Double_t ptStart        = histoYieldTrueSecFracMeson[k][j]->GetXaxis()->GetBinLowEdge(i);
                            Double_t ptEnd          = histoYieldTrueSecFracMeson[k][j]->GetXaxis()->GetBinUpEdge(i);
                            Double_t binWidth       = ptEnd-ptStart;
                            Double_t secFrac        = fitSecFracPurePowerlaw[k][j]->Integral(ptStart, ptEnd, resultCurr->GetParams()) / binWidth;
                            Double_t errorSecFrac   = fitSecFracPurePowerlaw[k][j]->IntegralError(ptStart, ptEnd, resultCurr->GetParams(), resultCurr->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                            histoYieldTrueSecFracMeson[k][j]->SetBinContent(i, secFrac);
                            histoYieldTrueSecFracMeson[k][j]->SetBinError(i, errorSecFrac);
                        }
                    // if secondary fraction for particular component not contained, put fraction to 0    
                    } else {
                        for (Int_t i = 1; i < histoYieldTrueSecFracMeson[k][j]->GetNbinsX()+1; i++){
                            histoYieldTrueSecFracMeson[k][j]->SetBinContent(i, 0);
                            histoYieldTrueSecFracMeson[k][j]->SetBinError(i, 0);                            
                        }    
                    }    
                }
            // fit current histo with powerlaw fit + constant (goes to constant  at high pt)        
            } else if ( doK0SecCorrectionWithDefaultHisto == 0 ){
                for (Int_t j = 0; j < 4; j++){
                    if (haveSec[j] && haveSecUsed[j]){
                        cout << "fitting current hist for "<< nameIntRange[k].Data() << " " << nameSecMeson[j].Data() << endl;
                        fitSecFracPLWithConst[k][j]->SetRange(minPtMesonSec, maxPtMeson);
                        TFitResultPtr resultCurr  = histoYieldTrueSecFracMeson[k][j]->Fit(fitSecFracPLWithConst[k][j],"SNRME+","",minPtMesonSec, maxPtMeson);
                        histoYieldTrueSecFracMeson[k][j]->Fit(fitSecFracPLWithConst[k][j],"NRME+","",minPtMesonSec, maxPtMeson);
                        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(histoYieldTrueSecFracMeson[k][j]);
                    // if secondary fraction for particular component not contained, put fraction to 0    
                    } else {
                        for (Int_t i = 1; i < histoYieldTrueSecFracMeson[k][j]->GetNbinsX()+1; i++){
                            histoYieldTrueSecFracMeson[k][j]->SetBinContent(i, 0);
                            histoYieldTrueSecFracMeson[k][j]->SetBinError(i, 0);                            
                        }    
                    }    
                }
            }   
        }    
        for (Int_t k = 0; k < 6; k++){
            // create secondary yield from MC fractions and reconstructed yield
            for (Int_t j = 0; j < 4; j++){
                histoYieldSecMeson[k][j]            = (TH1D*)histoUnCorrectedYield[k]->Clone(Form("SecYieldFrom%sMeson%s", nameSecMeson[j].Data(), nameIntRange[k].Data()));
                histoYieldSecMeson[k][j]->Sumw2();
                histoYieldSecMeson[k][j]->Multiply(histoYieldTrueSecFracMeson[k%3][j]);
                histoYieldSecMeson[k][j]->Scale(scalingFacSec[j]);

                histoRatioYieldSecMeson[k][j]       = (TH1D*)histoYieldSecMeson[k][j]->Clone(Form("RatioSecYieldFrom%sMeson%sToRaw", nameSecMeson[j].Data(), nameIntRange[k].Data()));
                histoRatioYieldSecMeson[k][j]->Sumw2();
                histoRatioYieldSecMeson[k][j]->Divide(histoRatioYieldSecMeson[k][j],histoUnCorrectedYield[k]);
            }
            
            // create secondary yield from external input (Cocktail or Toy) and sec effi and acc
            for (Int_t j = 0; j < 3; j++){
                if (histoSecAcceptance[j] && histoSecTrueEffi[k%3][j] && histoExternalInputSecPi0[j]){
                    histoYieldSecMesonFromExternalInput[k][j] = (TH1D*)histoExternalInputSecPi0[j]->Clone(Form("SecYieldFrom%sMeson%sFrom%s", nameSecMeson[j].Data(), nameIntRange[k].Data(),strExternalInputName.Data()));
                    histoYieldSecMesonFromExternalInput[k][j]->Sumw2();
                    histoYieldSecMesonFromExternalInput[k][j]->Multiply(histoSecAcceptance[j]);
                    histoYieldSecMesonFromExternalInput[k][j]->Multiply(histoSecTrueEffi[k%3][j]);

                    histoRatioYieldSecMesonFromExtInput[k][j]    = (TH1D*)histoYieldSecMesonFromExternalInput[k][j]->Clone(Form("RatioSecYieldFrom%sMeson%sFrom%sToRaw", nameSecMeson[j].Data(), nameIntRange[k].Data(),strExternalInputName.Data()));
                    histoRatioYieldSecMesonFromExtInput[k][j]->Sumw2();
                    histoRatioYieldSecMesonFromExtInput[k][j]->Scale(nEvt);
                    histoRatioYieldSecMesonFromExtInput[k][j]->Divide(histoRatioYieldSecMesonFromExtInput[k][j],histoUnCorrectedYield[k]);

                }
            }
            
            // resonance feed down contributions
            for (Int_t j = 0; j < 8; j++) {
                if (histoExternalInputFeedDownPi0[j]) {
                    histoYieldResonanceFeedDownPi0FromExternalInput[k][j] = (TH1D*)histoExternalInputFeedDownPi0[j]->Clone(Form("ResonanceFeedDownYieldFrom%sMeson%sFrom%s", nameSecMeson[j].Data(), nameIntRange[k].Data(),strExternalInputName.Data()));
                    histoYieldResonanceFeedDownPi0FromExternalInput[k][j]->Sumw2();
                    
                    histoRatioYieldResonanceFeedDownPi0FromExternalInput[k][j]    = (TH1D*)histoYieldResonanceFeedDownPi0FromExternalInput[k][j]->Clone(Form("RatioResonanceFeedDownYieldFrom%sMeson%sFrom%sToRaw", nameSecMeson[j].Data(), nameIntRange[k].Data(),strExternalInputName.Data()));
                    histoRatioYieldResonanceFeedDownPi0FromExternalInput[k][j]->Sumw2();
                    histoRatioYieldResonanceFeedDownPi0FromExternalInput[k][j]->Multiply(histoAcceptance);
                    histoRatioYieldResonanceFeedDownPi0FromExternalInput[k][j]->Multiply(histoTrueEffiPt[k%3]);
                    histoRatioYieldResonanceFeedDownPi0FromExternalInput[k][j]->Scale(nEvt);
                    histoRatioYieldResonanceFeedDownPi0FromExternalInput[k][j]->Divide(histoYieldResonanceFeedDownPi0FromExternalInput[k][j],histoUnCorrectedYield[k]);
                }
            }
        }
    // Do gamma-gamma correction for Dalitz channel
    } else if (doGGCorrection){
        for (Int_t k = 0; k < 3; k++){
            histoYieldTrueGGFracMeson[k]        = (TH1D*)fileCorrections->Get(Form("TrueGGFrac%s", nameIntRange[k].Data()));
            histoYieldGGMeson[k]                = (TH1D*)histoUnCorrectedYield[k]->Clone(Form("GGFracMeson%s", nameIntRange[k].Data()));
            histoYieldGGMeson[k]->Sumw2();
            histoYieldGGMeson[k]->Multiply(histoYieldTrueGGFracMeson[k]);
            
            histoYieldGGMeson[k+3]              = (TH1D*)histoUnCorrectedYield[k+3]->Clone(Form("GGFracMeson%s", nameIntRange[k+3].Data()));
            histoYieldGGMeson[k+3]->Sumw2();
            histoYieldGGMeson[k+3]->Multiply(histoYieldTrueGGFracMeson[k]);
        }
    }    

    Int_t nSecComp                  = 0;
    Int_t nSecCompUsed              = 0;
    for (Int_t j = 0; j < 4; j++){
        if (haveSec[j]) nSecComp++;
        if (haveSecUsed[j]) nSecCompUsed++;
    }    
            
//     return;

    //*******************************************************************************************************
    //******************************* Read pileup correction file data **************************************
    //*******************************************************************************************************    
    TString fileNameDCAData                         = Form("%s/%s/%s_Data_GammaConvV1DCATestAnalysed%s.root",fCutSelection.Data(),optionEnergy.Data(),nameMeson.Data(),optionPeriod.Data());
    
    TFile* fileDCAAnalysisData                      = new TFile(fileNameDCAData.Data());
    Bool_t kDCAFileDataExists                       = kTRUE;
    if (fileDCAAnalysisData->IsZombie())         
        kDCAFileDataExists                          = kFALSE;
    cout << kDCAFileDataExists << endl;
    TH1D *histoFracCatvsPt[6];
    TH1D *histoFracIntHistBGvsPt[6];
    for (Int_t i = 0; i < 6 ; i++){
        histoFracCatvsPt[i]                         = NULL;
        histoFracIntHistBGvsPt[i]                   = NULL;
    }   
    TH1D* histoCorrectionFactorsHistvsPt            = NULL;
    TH1D* histoCorrectionFactorsFitvsPt             = NULL;
    TH1D* histoCorrectionFactorsHistvsPtCatA        = NULL;
    TH1D* histoCorrectionFactorsHistvsPtCatC        = NULL;
    TH1D* histoCorrectionFactorsHistvsPtCatD        = NULL;
    TH1D* histoDCAZUnderMesonAllCat_AllPt           = NULL;
    if (kDCAFileDataExists){
        histoCorrectionFactorsHistvsPt              = (TH1D*)fileDCAAnalysisData->Get("fHistCorrectionFactorsHistAllCat_vsPt");
        histoCorrectionFactorsFitvsPt               = (TH1D*)fileDCAAnalysisData->Get("fHistCorrectionFactorsFitAllCat_vsPt");
        histoCorrectionFactorsHistvsPtCatA          = (TH1D*)fileDCAAnalysisData->Get("fHistCorrectionFactorsHistvsPt_0");
        histoCorrectionFactorsHistvsPtCatC          = (TH1D*)fileDCAAnalysisData->Get("fHistCorrectionFactorsHistvsPt_1");
        histoCorrectionFactorsHistvsPtCatD          = (TH1D*)fileDCAAnalysisData->Get("fHistCorrectionFactorsHistvsPt_2");
        histoDCAZUnderMesonAllCat_AllPt             = (TH1D*)fileDCAAnalysisData->Get("HistDCAZUnderMesonAllCat_AllPt");
        for (Int_t i = 0; i < 6 ; i++){
            histoFracCatvsPt[i]                     = (TH1D*)fileDCAAnalysisData->Get(Form("fHistFracCat_%i_vsPt",i+1));
            histoFracIntHistBGvsPt[i]               = (TH1D*)fileDCAAnalysisData->Get(Form("fHistFracIntHistBGvsPt_Cat_%i_Variant_1",i+1));
        }    
    }

    //*******************************************************************************************************
    //********************************** Read pileup correction file MC *************************************
    //*******************************************************************************************************        
    TString fileNameDCAMonteCarlo                   = Form("%s/%s/%s_MC_GammaConvV1DCATestAnalysed%s.root",fCutSelection.Data(),optionEnergy.Data(),nameMeson.Data(),optionPeriod.Data());

    TFile* fileDCAAnalysisMonteCarlo                = new TFile(fileNameDCAMonteCarlo.Data());
    Bool_t kDCAFileMCExists                         = kTRUE;
    if (fileDCAAnalysisMonteCarlo->IsZombie())     
        kDCAFileMCExists                            = kFALSE;
    cout << kDCAFileMCExists << endl;
    
    TH1D *histoMCFracCatvsPt[6];
    TH1D *histoMCFracIntHistBGvsPt[6];
    for (Int_t i = 0; i < 6 ; i++){
        histoMCFracCatvsPt[i]                       = NULL;
        histoMCFracIntHistBGvsPt[i]                 = NULL;
    }
    TH1D* histoMCDCAZUnderMesonAllCatDecomp_AllPt[7]    = { NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TString nameDiffCatDecomp[7]                        = { "HistDCAZUnderMesonAllCat_AllPt", "HistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt", "HistDCAZTruePrimaryMesonDalitzAllCat_AllPt",
                                                            "HistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt", "HistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt", "HistDCAZGarbageAllCat_AllPt", 
                                                            "HistDCAZTrueBackgroundAllCat_AllPt" };
    TString nameDiffCatDecompPlot[7]                    = { "Total MC", "Prim", "Dalitz", "Sec. #pi^{0} from K^{0}_{s}", "Sec. #pi^{0} from X", "garbage", "#gamma#gamma BG" };
                                                                        
    Color_t colorDiffCatDecomp[7]                       = {kGray, kRed+2, kGreen+2, 807, kViolet+2, kCyan+2, kPink+2};
    Style_t markerStyleCatDecomp[7]                     = {24, 20, 20, 20, 20, 20, 20};
    Size_t markerSizeCatDecomp[7]                       = {1, 1, 1, 1, 1, 1, 1};

    if (kDCAFileMCExists){
        for (Int_t l = 0; l < 7; l++){
            histoMCDCAZUnderMesonAllCatDecomp_AllPt[l]              = (TH1D*)fileDCAAnalysisMonteCarlo->Get(nameDiffCatDecomp[l].Data());
        }    
        for (Int_t i = 0; i < 6 ; i++){
            histoMCFracCatvsPt[i]                                   = (TH1D*)fileDCAAnalysisMonteCarlo->Get(Form("fHistFracCat_%i_vsPt",i+1));
            histoMCFracIntHistBGvsPt[i]                             = (TH1D*)fileDCAAnalysisMonteCarlo->Get(Form("fHistFracIntHistBGvsPt_Cat_%i_Variant_1",i+1));
        }
    }
    
    TH1D *histoBGEstimateA                          = (TH1D*)histoUnCorrectedYield[0]->Clone("histoBGEstimateA");
    TH1D *histoBGEstimateB                          = (TH1D*)histoUnCorrectedYield[0]->Clone("histoBGEstimateB");
    TH1D *histoBGEstimateCatA                       = (TH1D*)histoUnCorrectedYield[0]->Clone("histoBGEstimateCatA");
    TH1D *histoBGEstimateCatC                       = (TH1D*)histoUnCorrectedYield[0]->Clone("histoBGEstimateCatC");
    TH1D *histoBGEstimateCatD                       = (TH1D*)histoUnCorrectedYield[0]->Clone("histoBGEstimateCatD");
    
    // ************************************************************************************************
    // ********************** Plot dca distribution with MC component for all pT **********************
    // ************************************************************************************************
    if (kDCAFileDataExists && kDCAFileMCExists){
        cout << "Plotting dca distribution with MC component for all pT" << endl;
        TCanvas* canvasDCAMCComponents = new TCanvas("canvasDCAMCComponents","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasDCAMCComponents, 0.08, 0.02, 0.02, 0.09);
        canvasDCAMCComponents->SetLogy();
        
            TLegend* legendDCAMCComponents0 = GetAndSetLegend2(0.7,0.7,0.85,0.95, 0.04*900,1); 
            if (histoDCAZUnderMesonAllCat_AllPt){
                DrawAutoGammaMesonHistos( histoDCAZUnderMesonAllCat_AllPt, 
                                "","dca_{z} #gamma (cm)", "d(dca_{z})/#it{N}_{evt}", 
                                kFALSE, 2.,1, kFALSE,
                                kTRUE,1e-8,10*histoMCDCAZUnderMesonAllCatDecomp_AllPt[0]->GetMaximum(), 
                                kTRUE, -6., 6.);
                DrawGammaSetMarker(histoDCAZUnderMesonAllCat_AllPt, 20, 1., kBlack, kBlack);
                histoDCAZUnderMesonAllCat_AllPt->GetYaxis()->SetTitleOffset(0.8);
                histoDCAZUnderMesonAllCat_AllPt->DrawCopy("p,e1"); 
                legendDCAMCComponents0->AddEntry(histoDCAZUnderMesonAllCat_AllPt,"Data","p");
            }
            for (Int_t l = 0; l < 7; l++){
                if (histoMCDCAZUnderMesonAllCatDecomp_AllPt[l] && histoMCDCAZUnderMesonAllCatDecomp_AllPt[l]->GetEntries() > 0){
                    DrawGammaSetMarker(histoMCDCAZUnderMesonAllCatDecomp_AllPt[l], markerStyleCatDecomp[l], markerSizeCatDecomp[l], colorDiffCatDecomp[l], colorDiffCatDecomp[l]);
                    histoMCDCAZUnderMesonAllCatDecomp_AllPt[l]->DrawCopy("same,p,e1"); 
                    legendDCAMCComponents0->AddEntry(histoMCDCAZUnderMesonAllCatDecomp_AllPt[l],nameDiffCatDecompPlot[l].Data(),"p");
                }   
            }
                
            TLatex *labelEnergy = new TLatex(0.11,0.9,Form("%s",collisionSystem.Data()));
            SetStyleTLatex( labelEnergy, 0.04,4);
            labelEnergy->Draw();
            legendDCAMCComponents0->Draw();
            
        canvasDCAMCComponents->Update(); 
        canvasDCAMCComponents->SaveAs(Form("%s/%s_MC_DCAzDecomposition.%s",outputDir.Data(),nameMeson.Data(),suffix.Data()));
    }   
    
    // ************************************************************************************************
    // ********************** Plot fraction of mesons in different categories *************************
    // ************************************************************************************************    
    if (kDCAFileDataExists){
        cout << "Plotting fraction of mesons in different categories" << endl;
        TCanvas* canvasCorrFrac = new TCanvas("canvasCorrFrac","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasCorrFrac, 0.08, 0.02, 0.02, 0.09);

        canvasCorrFrac->cd();
        
        DrawAutoGammaMesonHistos( histoFracCatvsPt[0], 
                        "", Form("#it{p}_{T,%s} (GeV/#it{c})",textProcess.Data()), Form("#it{N}_{%s per cat}/(#it{N}_{%s}) (%s)", textProcess.Data(), textProcess.Data(), "%"), 
                        kFALSE, 2.,1e-8, kFALSE,
                        kTRUE, 0, 100., 
                        kFALSE, 0., 10.);  
        DrawGammaSetMarker(histoFracCatvsPt[0], styleCat[0], 1., colorCat[0], colorCat[0]);
        histoFracCatvsPt[0]->GetYaxis()->SetTitleOffset(0.8);
        histoFracCatvsPt[0]->DrawCopy("p,e1"); 
    
        if(histoMCFracCatvsPt[0]){
            DrawGammaSetMarker(histoMCFracCatvsPt[0], styleCatMC[0], 1., colorCatMC[0], colorCatMC[0]);
            histoMCFracCatvsPt[0]->GetYaxis()->SetTitleOffset(0.8);
            histoMCFracCatvsPt[0]->DrawCopy("same,p,e1"); 
        }
        TLegend* legendFractionCat = GetAndSetLegend2(0.65,0.7,0.95,0.95, 0.04*900,1); 
        if (kDCAFileMCExists) legendFractionCat->SetNColumns(3);
            else legendFractionCat->SetNColumns(2);
        legendFractionCat->AddEntry((TObject*)0,"Cat 1","");
        legendFractionCat->AddEntry(histoFracCatvsPt[0],"Data","p");
        if(histoMCFracCatvsPt[0])legendFractionCat->AddEntry(histoMCFracCatvsPt[0],"MC","p");
        
        for (Int_t i = 1; i< 6; i++){
            DrawGammaSetMarker(histoFracCatvsPt[i], styleCat[i], 1., colorCat[i], colorCat[i]);
            histoFracCatvsPt[i]->DrawCopy("same,p,e1"); 
            if(histoMCFracCatvsPt[i]){
                DrawGammaSetMarker(histoMCFracCatvsPt[i], styleCatMC[i], 1., colorCatMC[i], colorCatMC[i]);
                histoMCFracCatvsPt[i]->DrawCopy("same,p,e1"); 
            }
            legendFractionCat->AddEntry((TObject*)0,Form("Cat %i",i+1),"");
            legendFractionCat->AddEntry(histoFracCatvsPt[i],"Data","p");
            if(histoMCFracCatvsPt[i])legendFractionCat->AddEntry(histoMCFracCatvsPt[i],"MC","p");
        }   
        legendFractionCat->Draw();
        canvasCorrFrac->Update(); 
        canvasCorrFrac->SaveAs(Form("%s/%s_FractionPerCatVsPt_ComparedToMC.%s",outputDir.Data(),nameMeson.Data(),suffix.Data()));
    }   
    
    // ************************************************************************************************
    // ***** Plot fraction of BG from out of bunch pileup in different cateogries compared to MC ******
    // ************************************************************************************************
    if (kDCAFileDataExists){
        cout << "Ploting fraction of BG from out of bunch pileup in different cateogries compared to MC" << endl;
        TCanvas* canvasCorrFrac = new TCanvas("canvasCorrFrac","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasCorrFrac, 0.08, 0.02, 0.02, 0.09);

        canvasCorrFrac->cd();
        
        DrawAutoGammaMesonHistos( histoFracIntHistBGvsPt[0], 
                    "", Form("#it{p}_{T,%s} (GeV/#it{c})", textMeson.Data()), "BG/Total (%)", 
                    kFALSE, 2.,1e-8, kFALSE,
                    kTRUE, 0, 20, 
                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoFracIntHistBGvsPt[0], styleCat[0], 1., colorCat[0], colorCat[0]);
        histoFracIntHistBGvsPt[0]->GetYaxis()->SetTitleOffset(0.8);
        histoFracIntHistBGvsPt[0]->DrawCopy("p,e1"); 
    
        if(histoMCFracIntHistBGvsPt[0]){
            DrawGammaSetMarker(histoMCFracIntHistBGvsPt[0], styleCatMC[0], 1., colorCatMC[0], colorCatMC[0]);
            histoMCFracIntHistBGvsPt[0]->DrawCopy("same,p,e1"); 
        }
        TLegend* legendFractionCat = GetAndSetLegend2(0.65,0.7,0.95,0.95, 0.04*900,1); 
        if (kDCAFileMCExists) legendFractionCat->SetNColumns(3);
            else legendFractionCat->SetNColumns(2);
        legendFractionCat->AddEntry((TObject*)0,"Cat 1","");
        legendFractionCat->AddEntry(histoFracIntHistBGvsPt[0],"Data","p");
        if(histoMCFracIntHistBGvsPt[0])legendFractionCat->AddEntry(histoMCFracIntHistBGvsPt[0],"MC","p");
        
        for (Int_t i = 1; i< 6; i++){
            DrawGammaSetMarker(histoFracIntHistBGvsPt[i], styleCat[i], 1., colorCat[i], colorCat[i]);
            histoFracIntHistBGvsPt[i]->DrawCopy("same,p,e1"); 
            if(histoMCFracIntHistBGvsPt[i]){
                DrawGammaSetMarker(histoMCFracIntHistBGvsPt[i], styleCatMC[i], 1., colorCatMC[i], colorCatMC[i]);
                histoMCFracIntHistBGvsPt[i]->DrawCopy("same,p,e1"); 
            }
            legendFractionCat->AddEntry((TObject*)0,Form("Cat %i",i+1),"");
            legendFractionCat->AddEntry(histoFracIntHistBGvsPt[i],"Data","p");
            if(histoMCFracIntHistBGvsPt[i])legendFractionCat->AddEntry(histoMCFracIntHistBGvsPt[i],"MC","p");
        }   
        legendFractionCat->Draw();
        canvasCorrFrac->Update(); 
        canvasCorrFrac->SaveAs(Form("%s/%s_FracBGOverIntHist_ComparedToMC.%s",outputDir.Data(),nameMeson.Data(),suffix.Data()));
    }   
    
    Color_t colorMethod[5] = {kBlack, kCyan+2, kRed+2, kGreen+2, kBlue+2};
    Style_t styleMethod[5] = {20,24,21,29,33};

    // ************************************************************************************************
    // ********** Plot total contamination from  out of bunch pileup with different methods ***********
    // ************************************************************************************************
    Double_t maxFracBG  = 8;
    if (optionEnergy.CompareTo("2.76TeV") == 0) 
        maxFracBG       = 20;
    if (optionEnergy.CompareTo("8TeV") == 0) 
        maxFracBG       = 25;
    if (kDCAFileDataExists){
        cout << "Plotting total contamination from  out of bunch pileup with different methods" << endl;
        TCanvas* canvasCorrFrac = new TCanvas("canvasCorrFrac","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasCorrFrac, 0.08, 0.02, 0.02, 0.09);

        canvasCorrFrac->cd();
        if(nameMeson.Contains("Eta")){
            DrawAutoGammaMesonHistos( histoCorrectionFactorsHistvsPt, 
                        "", "#it{p}_{T,#eta} (GeV/#it{c})", "Contamination from Pileup (%)", 
                        kFALSE, 2.,1e-8, kFALSE,
                        kTRUE, 0, maxFracBG, 
                        kFALSE, 0., 10.);
        } else {
            DrawAutoGammaMesonHistos( histoCorrectionFactorsHistvsPt, 
                                    "", "#it{p}_{T,#pi^{0}} (GeV/#it{c})", "Contamination from Pileup (%)", 
                                    kFALSE, 2.,1e-8, kFALSE,
                                    kTRUE, 0, maxFracBG, 
                                    kFALSE, 0., 10.);
        }
        DrawGammaSetMarker(histoCorrectionFactorsHistvsPt, styleMethod[0], 1.2, colorMethod[0], colorMethod[0]);
        histoCorrectionFactorsHistvsPt->GetYaxis()->SetTitleOffset(0.8);
        histoCorrectionFactorsHistvsPt->DrawCopy("p,e1"); 
        TF1* fitCorrectionFactorsHistvsPt = new TF1("fitCorrectionFactorsHistvsPt","[0]/pow(x,[1])+[2]");
        fitCorrectionFactorsHistvsPt->SetRange(0.4, maxPtMeson);
        TFitResultPtr resultCorrectionFactorsHistvsPt = histoCorrectionFactorsHistvsPt->Fit(fitCorrectionFactorsHistvsPt,"SNRME+","",0.4, maxPtMeson);
        TString bla= WriteParameterToFile(fitCorrectionFactorsHistvsPt);
        cout << bla.Data()<< endl;
        fitCorrectionFactorsHistvsPt->SetLineColor(colorMethod[0]);
        fitCorrectionFactorsHistvsPt->Draw("same");
        
        TF1* fitCorrectionFactorsFitvsPt = NULL;
        TFitResultPtr resultCorrectionFactorsFitvsPt ;
        if (histoCorrectionFactorsFitvsPt) {
            DrawGammaSetMarker(histoCorrectionFactorsFitvsPt, styleMethod[1], 1.2, colorMethod[1], colorMethod[1]);
            histoCorrectionFactorsFitvsPt->DrawCopy("same,p,e1"); 
            fitCorrectionFactorsFitvsPt = new TF1("fitCorrectionFactorsFitvsPt","[0]/pow(x,[1])+[2]");
            fitCorrectionFactorsFitvsPt->SetRange(0.4, maxPtMeson);
            resultCorrectionFactorsFitvsPt = histoCorrectionFactorsFitvsPt->Fit(fitCorrectionFactorsFitvsPt,"SNRME+","",0.4, maxPtMeson);
            fitCorrectionFactorsFitvsPt->SetLineColor(colorMethod[1]);
            fitCorrectionFactorsFitvsPt->Draw("same");
        }
        DrawGammaSetMarker(histoCorrectionFactorsHistvsPtCatA, styleMethod[2], 1.2, colorMethod[2], colorMethod[2]);
        histoCorrectionFactorsHistvsPtCatA->DrawCopy("same,p,e1"); 
        TF1* fitCorrectionFactorsHistvsPtCatA = new TF1("fitCorrectionFactorsHistvsPtCatA","[0]/pow(x,[1])+[2]");
        fitCorrectionFactorsHistvsPtCatA->SetRange(0.4, maxPtMeson);
        
        TFitResultPtr resultCorrectionFactorsHistvsPtCatA = histoCorrectionFactorsHistvsPtCatA->Fit(fitCorrectionFactorsHistvsPtCatA,"SNRME+","",0.4, maxPtMeson);
        fitCorrectionFactorsHistvsPtCatA->SetLineColor(colorMethod[2]);
        fitCorrectionFactorsHistvsPtCatA->Draw("same");

        DrawGammaSetMarker(histoCorrectionFactorsHistvsPtCatC, styleMethod[3], 1.2, colorMethod[3], colorMethod[3]);
        histoCorrectionFactorsHistvsPtCatC->DrawCopy("same,p,e1"); 
        TF1* fitCorrectionFactorsHistvsPtCatC = new TF1("fitCorrectionFactorsHistvsPtCatC","[0]/pow(x,[1])+[2]");
        fitCorrectionFactorsHistvsPtCatC->SetRange(0.4, maxPtMeson);
        TFitResultPtr resultCorrectionFactorsHistvsPtCatC = histoCorrectionFactorsHistvsPtCatC->Fit(fitCorrectionFactorsHistvsPtCatC,"SNRME+","",0.4, maxPtMeson);
        fitCorrectionFactorsHistvsPtCatC->SetLineColor(colorMethod[3]);
        fitCorrectionFactorsHistvsPtCatC->Draw("same");

        DrawGammaSetMarker(histoCorrectionFactorsHistvsPtCatD, styleMethod[4], 1.2, colorMethod[4], colorMethod[4]);
        histoCorrectionFactorsHistvsPtCatD->DrawCopy("same,p,e1"); 
        TF1* fitCorrectionFactorsHistvsPtCatD = new TF1("fitCorrectionFactorsHistvsPtCatD","[0]/pow(x,[1])+[2]");
        fitCorrectionFactorsHistvsPtCatD->SetRange(0.4, maxPtMeson);
        TFitResultPtr resultCorrectionFactorsHistvsPtCatD = histoCorrectionFactorsHistvsPtCatD->Fit(fitCorrectionFactorsHistvsPtCatD,"SNRME+","",0.4, maxPtMeson);
        fitCorrectionFactorsHistvsPtCatD->SetLineColor(colorMethod[4]);
        fitCorrectionFactorsHistvsPtCatD->Draw("same");
                
        TLegend* legendFractionCat = new TLegend(0.5,0.8,0.95,0.95);
        legendFractionCat->SetTextSize(0.03);
        legendFractionCat->SetFillColor(0);
        legendFractionCat->SetLineColor(0);
        legendFractionCat->SetNColumns(2);
        legendFractionCat->AddEntry((TObject*)0,"Method A","");
        legendFractionCat->AddEntry(histoCorrectionFactorsHistvsPt,"Data","p");
        if (histoCorrectionFactorsFitvsPt) legendFractionCat->AddEntry((TObject*)0,"Method B","");
        if (histoCorrectionFactorsFitvsPt) legendFractionCat->AddEntry(histoCorrectionFactorsFitvsPt,"Data","p");
        legendFractionCat->AddEntry((TObject*)0,"Method A sep Cat","");
        legendFractionCat->AddEntry(histoCorrectionFactorsHistvsPtCatA,"Data","p");
        legendFractionCat->AddEntry((TObject*)0,"Method C sep Cat","");
        legendFractionCat->AddEntry(histoCorrectionFactorsHistvsPtCatC,"Data","p");
        legendFractionCat->AddEntry((TObject*)0,"Method D sep Cat","");
        legendFractionCat->AddEntry(histoCorrectionFactorsHistvsPtCatD,"Data","p");
    
        TLatex *labelEnergy = new TLatex(0.11,0.9,Form("%s", collisionSystem.Data()));
        SetStyleTLatex( labelEnergy, 0.04,4);
        labelEnergy->Draw();
        
        legendFractionCat->Draw();
        canvasCorrFrac->Update(); 
        canvasCorrFrac->SaveAs(Form("%s/%s_FinalBGEstimate_AllMethods.%s",outputDir.Data(),nameMeson.Data(),suffix.Data()));

        // ************************************************************************************************
        // ** Calculate and plot final correction factor for  out of bunch pileup with different methods **
        // ************************************************************************************************
        
        for (Int_t i = 2; i < histoBGEstimateA->GetNbinsX()+1; i++){
            Double_t ptStart = histoBGEstimateA->GetXaxis()->GetBinLowEdge(i);
            Double_t ptEnd = histoBGEstimateA->GetXaxis()->GetBinUpEdge(i);
            Double_t binWidth = ptEnd-ptStart;
            Double_t bgEstimate = (100-fitCorrectionFactorsHistvsPt->Integral(ptStart, ptEnd, resultCorrectionFactorsHistvsPt->GetParams()) / binWidth )/100.;
            Double_t errorBGEstimate = (fitCorrectionFactorsHistvsPt->IntegralError(ptStart, ptEnd, resultCorrectionFactorsHistvsPt->GetParams(), 
                                                                                    resultCorrectionFactorsHistvsPt->GetCovarianceMatrix().GetMatrixArray() ) / binWidth )/100.;
            histoBGEstimateA->SetBinContent(i, bgEstimate);
            histoBGEstimateA->SetBinError(i, errorBGEstimate);
            
            if (fitCorrectionFactorsFitvsPt){
                bgEstimate = (100-fitCorrectionFactorsFitvsPt->Integral(ptStart, ptEnd, resultCorrectionFactorsFitvsPt->GetParams()) / binWidth )/100.;
                errorBGEstimate = (fitCorrectionFactorsFitvsPt->IntegralError(ptStart, ptEnd, resultCorrectionFactorsFitvsPt->GetParams(), 
                                                                              resultCorrectionFactorsFitvsPt->GetCovarianceMatrix().GetMatrixArray() ) / binWidth )/100.;
                histoBGEstimateB->SetBinContent(i, bgEstimate);
                histoBGEstimateB->SetBinError(i, errorBGEstimate);
            }
            
            bgEstimate = (100-fitCorrectionFactorsHistvsPtCatA->Integral(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatA->GetParams()) / binWidth )/100.;
            errorBGEstimate = (fitCorrectionFactorsHistvsPtCatA->IntegralError(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatA->GetParams(), 
                                                                               resultCorrectionFactorsHistvsPtCatA->GetCovarianceMatrix().GetMatrixArray() ) / binWidth )/100.;
            histoBGEstimateCatA->SetBinContent(i, bgEstimate);
            histoBGEstimateCatA->SetBinError(i, errorBGEstimate);
            
            bgEstimate = (100-fitCorrectionFactorsHistvsPtCatC->Integral(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatC->GetParams()) / binWidth )/100.;
            errorBGEstimate = (fitCorrectionFactorsHistvsPtCatC->IntegralError(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatC->GetParams(), 
                                                                               resultCorrectionFactorsHistvsPtCatC->GetCovarianceMatrix().GetMatrixArray() ) / binWidth )/100.;
            histoBGEstimateCatC->SetBinContent(i, bgEstimate);
            histoBGEstimateCatC->SetBinError(i, errorBGEstimate);
            
            bgEstimate = (100-fitCorrectionFactorsHistvsPtCatD->Integral(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatD->GetParams()) / binWidth )/100.;
            errorBGEstimate = (fitCorrectionFactorsHistvsPtCatD->IntegralError(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatD->GetParams(), 
                                                                               resultCorrectionFactorsHistvsPtCatD->GetCovarianceMatrix().GetMatrixArray() ) / binWidth )/100.;
            histoBGEstimateCatD->SetBinContent(i, bgEstimate);
            histoBGEstimateCatD->SetBinError(i, errorBGEstimate);
        }

        canvasCorrFrac->cd();
        if(nameMeson.Contains("Eta")){
            DrawAutoGammaMesonHistos( histoBGEstimateA, 
                                "", "#it{p}_{T,#eta} (GeV/#it{c})", "Correction factor", 
                                kFALSE, 2.,1e-8, kFALSE,
                                kTRUE, 0.8, 1,
                                kFALSE, 0., 10.);
        } else {
            DrawAutoGammaMesonHistos( histoBGEstimateA, 
                                "", "#it{p}_{T,#pi^{0}} (GeV/#it{c})", "Correction factor", 
                                kFALSE, 2.,1e-8, kFALSE,
                                kTRUE, 0.7, 1,
                                kFALSE, 0., 10.);
        }
        DrawGammaSetMarker(histoBGEstimateA, styleMethod[0], 1.2, colorMethod[0], colorMethod[0]);
        histoBGEstimateA->GetYaxis()->SetTitleOffset(0.8);
        histoBGEstimateA->DrawCopy("p,e1"); 
        
        if ( fitCorrectionFactorsFitvsPt) {
            DrawGammaSetMarker(histoBGEstimateB, styleMethod[1], 1.2, colorMethod[1], colorMethod[1]);
            histoBGEstimateB->DrawCopy("same,p,e1"); 
        }
        
        DrawGammaSetMarker(histoBGEstimateCatA, styleMethod[2], 1.2, colorMethod[2], colorMethod[2]);
        histoBGEstimateCatA->DrawCopy("same,p,e1"); 
        
        DrawGammaSetMarker(histoBGEstimateCatC, styleMethod[3], 1.2, colorMethod[3], colorMethod[3]);
        histoBGEstimateCatC->DrawCopy("same,p,e1"); 
        
        DrawGammaSetMarker(histoBGEstimateCatD, styleMethod[4], 1.2, colorMethod[4], colorMethod[4]);
        histoBGEstimateCatD->DrawCopy("same,p,e1"); 
                
        TLegend* legendDiffMethods = new TLegend(0.55,0.15,0.95,0.3);
        legendDiffMethods->SetTextSize(0.03);
        legendDiffMethods->SetFillColor(0);
        legendDiffMethods->SetLineColor(0);
        legendDiffMethods->SetNColumns(2);
        legendDiffMethods->AddEntry((TObject*)0,"Method A","");
        legendDiffMethods->AddEntry(histoBGEstimateA,"Data","p");
        if (fitCorrectionFactorsFitvsPt){
            legendDiffMethods->AddEntry((TObject*)0,"Method B","");
            legendDiffMethods->AddEntry(histoBGEstimateB,"Data","p");
        }
        legendDiffMethods->AddEntry((TObject*)0,"Method A sep Cat","");
        legendDiffMethods->AddEntry(histoBGEstimateCatA,"Data","p");
        legendDiffMethods->AddEntry((TObject*)0,"Method C sep Cat","");
        legendDiffMethods->AddEntry(histoBGEstimateCatC,"Data","p");
        legendDiffMethods->AddEntry((TObject*)0,"Method D sep Cat","");
        legendDiffMethods->AddEntry(histoBGEstimateCatD,"Data","p");
        legendDiffMethods->Draw();
        
        TLatex *labelEnergy2 = new TLatex(0.11,0.15,Form("%s", collisionSystem.Data()));
        SetStyleTLatex( labelEnergy2, 0.04,4);
        labelEnergy2->Draw();

        canvasCorrFrac->Update(); 
        canvasCorrFrac->SaveAs(Form("%s/%s_FinalCorrectionFactor_AllMethods.%s",outputDir.Data(),nameMeson.Data(),suffix.Data()));
    }   
    
    cout << "made it!!" << endl;

    //********************************************************************************************************
    //************************ Subtracting BG from out of bunch pileup ***************************************
    //************************ and correcting secondary yield accordingly ************************************
    //********************************************************************************************************
    if (!kIsMC && !optDalitz && kDCAFileDataExists){
        for (Int_t k = 0; k < 6; k++){
            cout << "subtracting yield from out of bunch collisions for int range: " << nameIntRange[k].Data() << endl;
            histoUnCorrectedYield[k]->Multiply(histoBGEstimateCatA);
            // recalculate secondary yield in full MC approach with raw yield corrected for pileup
            if (doK0SecCorrection){
                for (Int_t j = 0; j < 4; j++){
                    cout << "recalculating sec yield for "<< nameSecMeson[j].Data() << " in int range: " << nameIntRange[k].Data() << endl;
                    histoYieldSecMeson[k][j]            = (TH1D*)histoUnCorrectedYield[k]->Clone(Form("SecYieldFrom%sMeson%s", nameSecMeson[j].Data(), nameIntRange[k].Data()));
                    histoYieldSecMeson[k][j]->Sumw2();
                    histoYieldSecMeson[k][j]->Multiply(histoYieldTrueSecFracMeson[k%3][j]); // was [k]
                    histoYieldSecMeson[k][j]->Scale(scalingFacSec[j]);  
                    
                    histoRatioYieldSecMeson[k][j]       = (TH1D*)histoYieldSecMeson[k][j]->Clone(Form("RatioSecYieldFrom%sMeson%sToRaw", nameSecMeson[j].Data(), nameIntRange[k].Data()));
                    histoRatioYieldSecMeson[k][j]->Sumw2();
                    histoRatioYieldSecMeson[k][j]->Divide(histoRatioYieldSecMeson[k][j],histoUnCorrectedYield[k]);
                }
            }    
        }    
    }



    
    //**********************************************************************************
    //******************** Mass Plot *********************************************
    //**********************************************************************************
    
    TH1D* histoRatioRecMass         = NULL;
    TH1D* histoRatioValRecMass      = NULL;
    TH1D* histoRatioRecMassGauss    = NULL; 
    TH1D* histoRatioValRecMassGauss = NULL;
    if (!kIsMC){ 
        TCanvas* canvasMass = new TCanvas("canvasMass","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasMass, 0.092, 0.01, 0.02, 0.082);
        
        if ( !kIsEta ){
            histoMassMeson->GetYaxis()->SetRangeUser(0.130,0.140);
            if (mode == 4 ) histoMassMeson->GetYaxis()->SetRangeUser(0.122,0.140);
            if (mode == 4 && optionEnergy.CompareTo("8TeV")==0 && trigger.CompareTo("81")==0) histoMassMeson->GetYaxis()->SetRangeUser(0.13,0.180);
            if (mode == 2 ) histoMassMeson->GetYaxis()->SetRangeUser(0.128,0.140);
            if (kCollisionSystem == 1 && mode > 1) histoMassMeson->GetYaxis()->SetRangeUser(0.130,0.155);
        } else {
            histoMassMeson->GetYaxis()->SetRangeUser(0.54,0.56);
            if (mode == 2 || mode == 4 ) histoMassMeson->GetYaxis()->SetRangeUser(0.50,0.57);
        }               
        histoMassMeson->GetYaxis()->SetNdivisions(510); 
        
        DrawAutoGammaMesonHistos( histoMassMeson, 
                                    "", "#it{p}_{T} (GeV/#it{c})", Form("Mass for %s in |#it{y}| < %s (GeV/#it{c}^{2})",textMeson.Data(), rapidityRange.Data()), 
                                    kFALSE, 0., 0.7, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoMassMeson, 20, 0.8, kBlack, kBlack); 
        histoMassMeson->DrawCopy("same,e1,p");
        DrawGammaSetMarker(histoTrueMassMeson, 24, 0.8, kRed+2, kRed+2);
        histoTrueMassMeson->DrawCopy("same,e1,p"); 
        
        DrawGammaLines(0., maxPtMeson,mesonMassExpect, mesonMassExpect,0.1);
        
        TLegend* legendMass = GetAndSetLegend2(0.65, 0.13, 0.95, 0.13+(0.035*3), 0.035, 1, "", 42, 0.15);
        legendMass->AddEntry(histoMassMeson,"reconstructed Data");
        if (histoMCrecMassMeson){
            DrawGammaSetMarker(histoMCrecMassMeson, 21, 0.8, kRed-4, kRed-4);
            histoMCrecMassMeson->DrawCopy("same,e1,p"); 
            legendMass->AddEntry(histoMCrecMassMeson,"reconstructed MC");
        }
        if( !kIsEta ) legendMass->AddEntry(histoTrueMassMeson,"True reconstructed #pi^{0}");
        if( kIsEta ) legendMass->AddEntry(histoTrueMassMeson,"True reconstructed #eta");

        
        legendMass->Draw();
        PutProcessLabelAndEnergyOnPlot(0.15, 0.97, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
        
        canvasMass->Update();
        canvasMass->SaveAs(Form("%s/%s_Mass_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));

        if (histoMassMeson && histoTrueMassMeson && histoMCrecMassMeson){

            TCanvas* canvasMassRatio = new TCanvas("canvasMassRatio","",200,10,1350,900);  // gives the page size
            DrawGammaCanvasSettings( canvasMassRatio, 0.10, 0.01, 0.02, 0.10);

            histoRatioRecMass           = (TH1D*)histoMCrecMassMeson->Clone("histoRatioRecMass");
            histoRatioValRecMass        = (TH1D*)histoTrueMassMeson->Clone("histoRatioValRecMass");
            histoRatioRecMass->Divide(histoRatioRecMass, histoMassMeson, 1., 1., "");
            histoRatioValRecMass->Divide(histoRatioValRecMass, histoMassMeson, 1., 1., "");

            TF1* fitPol0                = new TF1("fitPol0","[0]",1.5,maxPtMeson);
            histoRatioRecMass->Fit(fitPol0,"NRME+","",1.5,maxPtMeson);
            Double_t recMassRatio       = fitPol0->GetParameter(0);
            Double_t recMassRatioError  = fitPol0->GetParError(0);
            histoRatioValRecMass->Fit(fitPol0,"NRME+","",1.5,maxPtMeson);
            Double_t valMassRatio       = fitPol0->GetParameter(0);
            Double_t valMassRatioError  = fitPol0->GetParError(0);
            
            DrawGammaSetMarker(histoRatioRecMass, 20, 0.8, kBlack, kBlack); 
            DrawAutoGammaMesonHistos( histoRatioRecMass, 
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("Ratio m_{MC}/m_{data} for %s in |#it{y}| < %s (GeV/#it{c}^{2})",textMeson.Data(), rapidityRange.Data()), 
                                        kFALSE, 0., 0.7, kFALSE,
                                        kTRUE, 0.95, 1.05, 
                                        kFALSE, 0., 10.);
            DrawGammaSetMarker(histoRatioRecMass, 20, 0.8, kBlack, kBlack); 
            histoRatioRecMass->DrawCopy("e1,p");
            DrawGammaSetMarker(histoRatioValRecMass, 24, 0.8, kRed+2, kRed+2);
            histoRatioValRecMass->DrawCopy("same,e1,p"); 
            
            DrawGammaLines(0., maxPtMeson,1, 1,0.1);
            
            TLegend* legendMassRatio = GetAndSetLegend2(0.15, 0.12, 0.6, 0.12+(0.035*4), 0.035, 2, "", 42, 0.15);
            legendMassRatio->AddEntry(histoRatioRecMass,"rec MC/ data");
            legendMassRatio->AddEntry((TObject*)0,Form("%0.4f #pm %0.4f", recMassRatio, recMassRatioError),"");
            legendMassRatio->AddEntry(histoRatioValRecMass,"val. rec MC/ data");
            legendMassRatio->AddEntry((TObject*)0,Form("%0.4f #pm %0.4f", valMassRatio, valMassRatioError),"");
            
            if(histoMassGaussianMeson && histoMCrecMassGaussianMeson && histoTrueMassGaussianMeson){
                histoRatioRecMassGauss          = (TH1D*)histoMCrecMassGaussianMeson->Clone("histoRatioRecMassGauss");
                histoRatioValRecMassGauss       = (TH1D*)histoTrueMassGaussianMeson->Clone("histoRatioValRecMassGauss");
                histoRatioRecMassGauss->Divide(histoRatioRecMassGauss, histoMassGaussianMeson, 1., 1., "");
                histoRatioValRecMassGauss->Divide(histoRatioValRecMassGauss, histoMassGaussianMeson, 1., 1., "");
                histoRatioRecMassGauss->Fit(fitPol0,"NRME+","",1.5,maxPtMeson);
                Double_t recMassGaussRatio      = fitPol0->GetParameter(0);
                Double_t recMassGaussRatioError = fitPol0->GetParError(0);
                histoRatioValRecMassGauss->Fit(fitPol0,"NRME+","",1.5,maxPtMeson);
                Double_t valMassGaussRatio      = fitPol0->GetParameter(0);
                Double_t valMassGaussRatioError = fitPol0->GetParError(0);

                
                DrawGammaSetMarker(histoRatioRecMassGauss, 20, 0.8, kGray+2, kGray+2);
                histoRatioRecMassGauss->DrawCopy("same,e1,p"); 
            
                DrawGammaSetMarker(histoRatioValRecMassGauss, 24, 0.8, kGreen+4, kGreen+4);
                histoRatioValRecMassGauss->DrawCopy("same,e1,p"); 
                
                legendMassRatio->AddEntry(histoRatioRecMassGauss,"rec MC/ data Gauss"); 
                legendMassRatio->AddEntry((TObject*)0,Form("%0.4f #pm %0.4f", recMassGaussRatio, recMassGaussRatioError),"");
    
                legendMassRatio->AddEntry(histoRatioValRecMassGauss,"val. rec MC/ data Gauss");
                legendMassRatio->AddEntry((TObject*)0,Form("%0.4f #pm %0.4f", valMassGaussRatio, valMassGaussRatioError),"");
            }
                        
            legendMassRatio->Draw();
            
            PutProcessLabelAndEnergyOnPlot(0.15, 0.97, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            
            canvasMassRatio->Update();
            canvasMassRatio->SaveAs(Form("%s/%s_RatioMass_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
        }    
        
        //**********************************************************************************
        //******************** Mass Plot compared to pure Gaussian fit *********************
        //**********************************************************************************                
        if (histoMassGaussianMeson){
            canvasMass->cd();
            
            if ( !kIsEta ){
                histoMassMeson->GetYaxis()->SetRangeUser(0.125,0.150);
                if (mode == 4 ) histoMassMeson->GetYaxis()->SetRangeUser(0.122,0.140);
                if (mode == 2 ) histoMassMeson->GetYaxis()->SetRangeUser(0.125,0.140);
                if (kCollisionSystem == 1 && mode > 1) histoMassMeson->GetYaxis()->SetRangeUser(0.130,0.155);
            } else {
                histoMassMeson->GetYaxis()->SetRangeUser(0.52,0.58);
                if (mode == 2 || mode == 4 ) histoMassMeson->GetYaxis()->SetRangeUser(0.48,0.57);
            }               

            histoMassMeson->DrawCopy("e1,p"); 
            histoTrueMassMeson->DrawCopy("same,e1,p"); 
            
            // count lines to be put for legend
            Int_t nLinesLegend = 2;
            if (histoMCrecMassMeson) nLinesLegend++;
            if (histoMassGaussianMeson) nLinesLegend++;
            if (histoMCrecMassGaussianMeson) nLinesLegend++;
            if (histoTrueMassGaussianMeson) nLinesLegend++;
            
            TLegend* legendMass4 = GetAndSetLegend2(0.55, 0.12, 0.95, 0.12+(0.035*nLinesLegend), 0.035, 1, "", 42, 0.1); 
            legendMass4->AddEntry(histoMassMeson,"reconstructed Data");
            // put reconstructed MC if available
            if (histoMCrecMassMeson){
                histoMCrecMassMeson->DrawCopy("same,e1,p"); 
                legendMass4->AddEntry(histoMCrecMassMeson,"reconstructed MC");
            }
            if(!kIsEta ) legendMass4->AddEntry(histoTrueMassMeson,"True reconstructed #pi^{0}");
            if(kIsEta ) legendMass4->AddEntry(histoTrueMassMeson,"True reconstructed #eta");
            
            // put Mass obtained from pure gaussian if available
            if (histoMassGaussianMeson){
                DrawGammaSetMarker(histoMassGaussianMeson, 20, 0.8, kGray+2, kGray+2);
                histoMassGaussianMeson->DrawCopy("same,e1,p"); 
                legendMass4->AddEntry(histoMassGaussianMeson,"reconstructed Data, pure Gauss"); 
            }
            // put MC Mass obtained from pure gaussian if available
            if (histoMCrecMassGaussianMeson){
                DrawGammaSetMarker(histoMCrecMassGaussianMeson, 20, 0.8, kGreen-2, kGreen-2);
                histoMCrecMassGaussianMeson->DrawCopy("same,e1,p"); 
                legendMass4->AddEntry(histoMCrecMassGaussianMeson,"reconstructed MC, pure Gauss");
            }
            // put True MC Mass obtained from pure gaussian if available
            if (histoTrueMassGaussianMeson){
                DrawGammaSetMarker(histoTrueMassGaussianMeson, 24, 0.8, kGreen+4, kGreen+4);
                histoTrueMassGaussianMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendMass4->AddEntry(histoTrueMassGaussianMeson,"True reconstructed #pi^{0}, pure Gauss");
                if(kIsEta ) legendMass4->AddEntry(histoTrueMassGaussianMeson,"True reconstructed #eta, pure Gauss");
            }
            
            DrawGammaLines(0., maxPtMeson,mesonMassExpect, mesonMassExpect,0.1);
            PutProcessLabelAndEnergyOnPlot(0.15, 0.95, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            legendMass4->Draw();
            canvasMass->Update(); 
            canvasMass->SaveAs(Form("%s/%s_MassComparisonPureGaussian_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
        }
        
        //**********************************************************************************
        //******************** Mass Plot further decomposed for PCM + Calo *****************
        //**********************************************************************************        
        canvasMass->cd();
        if (mode==2 || mode == 3){
            // read additional histos
            TH1D* histoTrueMassCaloPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueMassMesonCaloPhoton");            
            TH1D* histoTrueMassCaloConvPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueMassMesonCaloConvPhoton");
            TH1D* histoTrueMassCaloMergedClusterMeson =          (TH1D*)fileCorrections->Get("histoTrueMassMesonCaloMergedCluster");

            // start plotting
            if ( !kIsEta ){
                histoMassMeson->GetYaxis()->SetRangeUser(0.120,0.140);
                if (kCollisionSystem == 1 && mode > 1) histoMassMeson->GetYaxis()->SetRangeUser(0.130,0.155);
                if (optionEnergy.CompareTo("8TeV")==0) histoMassMeson->GetYaxis()->SetRangeUser(0.120,0.155);
            } else {
                histoMassMeson->GetYaxis()->SetRangeUser(0.52,0.58);
                if (mode == 2 || mode == 4 ) histoMassMeson->GetYaxis()->SetRangeUser(0.48,0.57);
            }               

            histoMassMeson->DrawCopy("e1,p"); 
            histoTrueMassMeson->DrawCopy("same,e1,p"); 

            // count lines to be put for legend
            Int_t nLinesLegend = 2;
            if (histoMCrecMassMeson) nLinesLegend++;
            if (histoTrueMassCaloPhotonMeson) nLinesLegend++;
            if (histoTrueMassCaloConvPhotonMeson) nLinesLegend++;
            if (histoTrueMassCaloMergedClusterMeson) nLinesLegend++;
                        
            TLegend* legendMass2 = GetAndSetLegend2(0.55, 0.12, 0.95, 0.12+(0.035*nLinesLegend), 0.035, 1, "", 42, 0.1); 
            legendMass2->AddEntry(histoMassMeson,"reconstructed Data");            
            // put reconstructed MC if available
            if (histoMCrecMassMeson){
                histoMCrecMassMeson->DrawCopy("same,e1,p"); 
                legendMass2->AddEntry(histoMCrecMassMeson,"reconstructed MC");
            }
            if(!kIsEta ) legendMass2->AddEntry(histoTrueMassMeson,"True reconstructed #pi^{0}");
            if(kIsEta ) legendMass2->AddEntry(histoTrueMassMeson,"True reconstructed #eta");
            // put reconstructed validated real gammagamma Mass position
            if (histoTrueMassCaloPhotonMeson){
                DrawGammaSetMarker(histoTrueMassCaloPhotonMeson, 25, 0.8, kGreen+2, kGreen+2);
                histoTrueMassCaloPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendMass2->AddEntry(histoTrueMassCaloPhotonMeson,"True reconstructed #pi^{0}, cluster real #gamma");
                if(kIsEta ) legendMass2->AddEntry(histoTrueMassCaloPhotonMeson,"True reconstructed #eta, cluster real #gamma");
            }
            // put reconstructed validated real gamma gamma_{conv} Mass position
            if (histoTrueMassCaloConvPhotonMeson){
                DrawGammaSetMarker(histoTrueMassCaloConvPhotonMeson, 25, 0.8, kCyan+2, kCyan+2);
                histoTrueMassCaloConvPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendMass2->AddEntry(histoTrueMassCaloConvPhotonMeson,"True reconstructed #pi^{0}, cluster conv #gamma");
                if(kIsEta ) legendMass2->AddEntry(histoTrueMassCaloConvPhotonMeson,"True reconstructed #eta, cluster conv #gamma");
            }
            // put reconstructed validated real gamma gamma_{merged} Mass position
            if (histoTrueMassCaloMergedClusterMeson){
                DrawGammaSetMarker(histoTrueMassCaloMergedClusterMeson, 25, 0.8, kViolet+2, kViolet+2);
                histoTrueMassCaloMergedClusterMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendMass2->AddEntry(histoTrueMassCaloMergedClusterMeson,"True reconstructed #pi^{0}, merged cluster #gamma");
                if(kIsEta ) legendMass2->AddEntry(histoTrueMassCaloMergedClusterMeson,"True reconstructed #eta, merged cluster #gamma");
            }
            
            DrawGammaLines(0., maxPtMeson,mesonMassExpect, mesonMassExpect,0.1);
            PutProcessLabelAndEnergyOnPlot(0.15, 0.95, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            legendMass2->Draw();
            
            canvasMass->Update(); 
            canvasMass->SaveAs(Form("%s/%s_MassAddedInfos_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
        }

        //**********************************************************************************
        //******************** Mass Plot further decomposed for Calo + Calo *****************
        //**********************************************************************************        
        if (mode==4 || mode == 5){
            // read additional histos
            TH1D* histoTrueMassCaloPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueMassMesonCaloPhoton");
            TH1D* histoTrueMassCaloConvPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueMassMesonCaloConvPhoton");
            TH1D* histoTrueMassMixedCaloConvPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueMassMesonMixedCaloConvPhoton");
            if (!kIsEta){
                histoMassMeson->GetYaxis()->SetRangeUser(0.120,0.140);
                if (kCollisionSystem == 1 && mode > 1) histoMassMeson->GetYaxis()->SetRangeUser(0.130,0.155);
                if (optionEnergy.CompareTo("8TeV")==0) histoMassMeson->GetYaxis()->SetRangeUser(0.120,0.155);
            } else {
                histoMassMeson->GetYaxis()->SetRangeUser(0.52,0.58);
                if (mode == 2 || mode == 4 ) histoMassMeson->GetYaxis()->SetRangeUser(0.48,0.57);
            }               

            histoMassMeson->DrawCopy("e1,p"); 
            histoTrueMassMeson->DrawCopy("same,e1,p"); 

            // count lines to be put for legend
            Int_t nLinesLegend = 2;
            if (histoMCrecMassMeson) nLinesLegend++;
            if (histoTrueMassCaloPhotonMeson) nLinesLegend++;
            if (histoTrueMassCaloConvPhotonMeson) nLinesLegend++;
            if (histoTrueMassMixedCaloConvPhotonMeson) nLinesLegend++;
                        
            TLegend* legendMass3 = GetAndSetLegend2(0.55, 0.12, 0.95, 0.12+(0.035*nLinesLegend), 0.035, 1, "", 42, 0.1); 
            legendMass3->AddEntry(histoMassMeson,"reconstructed Data");
            // put reconstructed MC if available
            if (histoMCrecMassMeson){
                histoMCrecMassMeson->DrawCopy("same,e1,p"); 
                legendMass3->AddEntry(histoMCrecMassMeson,"reconstructed MC");
            }

            if(!kIsEta ) legendMass3->AddEntry(histoTrueMassMeson,"True reconstructed #pi^{0}");
            if(kIsEta ) legendMass3->AddEntry(histoTrueMassMeson,"True reconstructed #eta");
            // put reconstructed validated real gammagamma Mass position
            if (histoTrueMassCaloPhotonMeson){
                DrawGammaSetMarker(histoTrueMassCaloPhotonMeson, 25, 0.8, kGreen+2, kGreen+2);
                histoTrueMassCaloPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendMass3->AddEntry(histoTrueMassCaloPhotonMeson,"True reconstructed #pi^{0}, #gamma#gamma");
                if(kIsEta ) legendMass3->AddEntry(histoTrueMassCaloPhotonMeson,"True reconstructed #eta, #gamma#gamma");
            }
            // put reconstructed validated real gamma_{conv} gamma_{conv} Mass position
            if (histoTrueMassCaloConvPhotonMeson){
                DrawGammaSetMarker(histoTrueMassCaloConvPhotonMeson, 25, 0.8, kCyan+2, kCyan+2);
                histoTrueMassCaloConvPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendMass3->AddEntry(histoTrueMassCaloConvPhotonMeson,"True reconstructed #pi^{0}, #gamma_{conv}#gamma_{conv}");
                if(kIsEta ) legendMass3->AddEntry(histoTrueMassCaloConvPhotonMeson,"True reconstructed #eta, #gamma_{conv}#gamma_{conv}");
            }
            // put reconstructed validated real gamma gamma_{conv} Mass position
            if (histoTrueMassMixedCaloConvPhotonMeson){
                DrawGammaSetMarker(histoTrueMassMixedCaloConvPhotonMeson, 25, 0.8, kBlue+2, kBlue+2);
                histoTrueMassMixedCaloConvPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendMass3->AddEntry(histoTrueMassMixedCaloConvPhotonMeson,"True reconstructed #pi^{0}, #gamma#gamma_{conv}");
                if(kIsEta ) legendMass3->AddEntry(histoTrueMassMixedCaloConvPhotonMeson,"True reconstructed #eta, #gamma#gamma_{conv}");
            }
            
            DrawGammaLines(0., maxPtMeson,mesonMassExpect, mesonMassExpect,0.1);
            
            PutProcessLabelAndEnergyOnPlot(0.15, 0.95, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            legendMass3->Draw();
            canvasMass->Update(); 
            canvasMass->SaveAs(Form("%s/%s_MassAddedInfos_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
            delete legendMass3;
        }
        delete canvasMass;
        delete legendMass;
    }    
    //**********************************************************************************
    //******************** FWHM Plot *********************************************
    //**********************************************************************************
    if (!kIsMC){ 
        
        TCanvas* canvasFWHM = new TCanvas("canvasFWHM","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasFWHM, 0.07, 0.01, 0.031, 0.082);
        
        Double_t maxFWHM        = 0.030;
        if (mode == 2)
          maxFWHM               = 0.015;
        if (kIsEta) 
          maxFWHM               = 0.022;
        if (kIsEta && mode == 4)
          maxFWHM               = 0.060;
        if (kIsEta && mode == 2)
          maxFWHM               = 0.060;

        if (mode == 4 && optionEnergy.CompareTo("8TeV") == 0) maxFWHM = 0.05;

        Double_t minFWHM        = -0.004;
        if (kIsEta) minFWHM     = 0.00;
        if (mode == 4) minFWHM  = 0.00; 
        //if (mode == 2 && optionEnergy.CompareTo("8TeV") == 0) minFWHM = 0.003;

        
        histoFWHMMeson->Sumw2();
        histoFWHMMeson->Scale(1./2.35);
        DrawAutoGammaMesonHistos( histoFWHMMeson, 
                                    "", "#it{p}_{T} (GeV/#it{c})", Form("FWHM/2.35 for %s in |#it{y}| < %s (GeV/#it{c}^{2})",textMeson.Data(), rapidityRange.Data()), 
                                    kFALSE, 1.5,-20., kFALSE,
                                    kTRUE, minFWHM, maxFWHM, 
                                    kFALSE, 0., 10., 62, 
                                    0.04, 42, 0.03, 0.9, 0.8);
                                        
        TLegend* legendFWHM = GetAndSetLegend2(0.1, 0.12, 0.45, 0.12+(0.035*3), 0.035, 1, "", 42, 0.1); 
        legendFWHM->AddEntry(histoFWHMMeson,"reconstructed Data");
        
        DrawGammaSetMarker(histoFWHMMeson, 20, 0.8, kBlack, kBlack); 
        histoFWHMMeson->DrawCopy("same,e1,p"); 
        
        if (histoMCrecFWHMMeson){
            histoMCrecFWHMMeson->Scale(1./2.35);
            DrawGammaSetMarker(histoMCrecFWHMMeson, 21, 0.8, kRed-4, kRed-4);
            histoMCrecFWHMMeson->DrawCopy("same,e1,p"); 
            legendFWHM->AddEntry(histoMCrecFWHMMeson,"reconstructed MC");
        }

        histoTrueFWHMMeson->Sumw2();
        histoTrueFWHMMeson->Scale(1./2.35);  
        DrawGammaSetMarker(histoTrueFWHMMeson, 24, 0.8, kRed+2, kRed+2);
        histoTrueFWHMMeson->DrawCopy("same,e1,p"); 
        if(!kIsEta) legendFWHM->AddEntry(histoTrueFWHMMeson,"True reconstructed #pi^{0}");
        if(kIsEta ) legendFWHM->AddEntry(histoTrueFWHMMeson,"True reconstructed #eta");

        legendFWHM->Draw();
        PutProcessLabelAndEnergyOnPlot(0.7, 0.94, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
        canvasFWHM->Update();
        canvasFWHM->SaveAs(Form("%s/%s_FWHM_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));

        //**********************************************************************************
        //******************** FWHM Plot comparison to pure Gaussian width *****************
        //**********************************************************************************        
        
        if (histoWidthGaussianMeson){
            Double_t maxFWHMGaus    = 0.030;
            if (mode == 2)
                maxFWHMGaus         = 0.015;
            if (kIsEta) 
                maxFWHMGaus         = 0.022;
            if (kIsEta && mode == 4)
                maxFWHMGaus         = 0.060;
            if (kIsEta && mode == 2)
                maxFWHMGaus         = 0.060;
          
            histoFWHMMeson->GetYaxis()->SetRangeUser(minFWHM, maxFWHMGaus);
            histoFWHMMeson->DrawCopy("e1,p"); 
            histoTrueFWHMMeson->DrawCopy("same,e1,p"); 

            // count lines to be put for legend
            Int_t nLinesLegend = 2;
            if (histoMCrecFWHMMeson) nLinesLegend++;
            if (histoWidthGaussianMeson) nLinesLegend++;
            if (histoMCrecWidthGaussMeson) nLinesLegend++;
            if (histoTrueWidthGaussianMeson) nLinesLegend++;

            TLegend* legendFWHM4 = GetAndSetLegend2(0.1, 0.12, 0.45, 0.12+(0.035*nLinesLegend), 0.035, 1, "", 42, 0.1); 
            legendFWHM4->AddEntry(histoFWHMMeson,"reconstructed Data");
            //  put rec MC width if available
            if (histoMCrecFWHMMeson){
                histoMCrecFWHMMeson->DrawCopy("same,e1,p"); 
                legendFWHM4->AddEntry(histoMCrecFWHMMeson,"reconstructed MC");
            }
            if(!kIsEta ) legendFWHM4->AddEntry(histoTrueFWHMMeson,"True reconstructed #pi^{0}");
            if(kIsEta ) legendFWHM4->AddEntry(histoTrueFWHMMeson,"True reconstructed #eta");
            
            //  put gaussian data width if available
            if (histoWidthGaussianMeson){
                DrawGammaSetMarker(histoWidthGaussianMeson, 20, 1.0, kGray+2, kGray+2);
                histoWidthGaussianMeson->DrawCopy("same,e1,p"); 
                legendFWHM4->AddEntry(histoWidthGaussianMeson,"reconstructed Data, #sigma pure Gauss");
            }
            //  put gaussian MC width if available
            if (histoMCrecWidthGaussMeson){
                DrawGammaSetMarker(histoMCrecWidthGaussMeson, 20, 1.0, kGreen-2, kGreen-2);
                histoMCrecWidthGaussMeson->DrawCopy("same,e1,p"); 
                legendFWHM4->AddEntry(histoMCrecWidthGaussMeson,"reconstructed MC, #sigma pure Gauss");
            }
            //  put gaussian validated MC width if available
            if (histoTrueWidthGaussianMeson){
                DrawGammaSetMarker(histoTrueWidthGaussianMeson, 24, 1.0, kGreen+4, kGreen+4);
                histoTrueWidthGaussianMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendFWHM4->AddEntry(histoTrueWidthGaussianMeson,"True reconstructed #pi^{0}, #sigma pure Gauss");
                if(kIsEta ) legendFWHM4->AddEntry(histoTrueWidthGaussianMeson,"True reconstructed #eta, #sigma pure Gauss");
            }
            
            PutProcessLabelAndEnergyOnPlot(0.7, 0.94, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            legendFWHM4->Draw();
            canvasFWHM->Update(); 
            canvasFWHM->SaveAs(Form("%s/%s_FWHMComparisonPureGauss_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
        }
        //**********************************************************************************
        //******************** FWHM Plot further decomposed for PCM + Calo *****************
        //**********************************************************************************        
        Double_t maxFWHMAdd         = 0.030;
        if (mode == 2)      
            maxFWHMAdd              = 0.015;
        if ( (mode == 2 || mode == 4) && optionEnergy.CompareTo("8TeV")==0)
            maxFWHMAdd              = 0.025;
        if (kIsEta) 
            maxFWHMAdd              = 0.030;
        if (kIsEta && mode == 4)
            maxFWHMAdd              = 0.070;
        if (kIsEta && mode == 2)
            maxFWHMAdd              = 0.060;

        if (mode==2 || mode == 3){

            // read additional histos   
            TH1D* histoTrueFWHMCaloPhotonMeson          = (TH1D*)fileCorrections->Get("histoTrueFWHMMesonCaloPhoton");
            TH1D* histoTrueFWHMCaloConvPhotonMeson      = (TH1D*)fileCorrections->Get("histoTrueFWHMMesonCaloConvPhoton");
            TH1D* histoTrueFWHMCaloMergedClusterMeson   = (TH1D*)fileCorrections->Get("histoTrueFWHMMesonCaloMergedCluster");
            
            // start plotting
            histoFWHMMeson->GetYaxis()->SetRangeUser(minFWHM, maxFWHMAdd); 
            histoFWHMMeson->DrawCopy("e1,p"); 
            histoTrueFWHMMeson->DrawCopy("same,e1,p"); 

            // count lines to be put for legend
            Int_t nLinesLegend = 2;
            if (histoMCrecFWHMMeson) nLinesLegend++;
            if (histoTrueFWHMCaloPhotonMeson) nLinesLegend++;
            if (histoTrueFWHMCaloConvPhotonMeson) nLinesLegend++;
            if (histoTrueFWHMCaloMergedClusterMeson) nLinesLegend++;
            
            TLegend* legendFWHM2 = GetAndSetLegend2(0.1, 0.12, 0.45, 0.12+(0.035*nLinesLegend), 0.035, 1, "", 42, 0.1); 
            legendFWHM2->AddEntry(histoFWHMMeson,"reconstructed Data");
            // plot MC rec FWHM if available
            if (histoMCrecFWHMMeson){
                histoMCrecFWHMMeson->DrawCopy("same,e1,p"); 
                legendFWHM2->AddEntry(histoMCrecFWHMMeson,"reconstructed MC");
            }

            if(!kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMMeson,"True reconstructed #pi^{0}");
            if(kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMMeson,"True reconstructed #eta");
            
            // plot validated gamma gamma FWHM if available
            if (histoTrueFWHMCaloPhotonMeson){
                histoTrueFWHMCaloPhotonMeson->Scale(1./2.35);
                DrawGammaSetMarker(histoTrueFWHMCaloPhotonMeson, 25, 0.8, kGreen+2, kGreen+2);
                histoTrueFWHMCaloPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMCaloPhotonMeson,"True reconstructed #pi^{0}, cluster real #gamma");
                if(kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMCaloPhotonMeson,"True reconstructed #eta, cluster real #gamma");
            }
            // plot validated gamma gamma_{conv} FWHM if available
            if (histoTrueFWHMCaloConvPhotonMeson){
                histoTrueFWHMCaloConvPhotonMeson->Scale(1./2.35);
                DrawGammaSetMarker(histoTrueFWHMCaloConvPhotonMeson, 25, 0.8, kCyan+2, kCyan+2);
                histoTrueFWHMCaloConvPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMCaloConvPhotonMeson,"True reconstructed #pi^{0}, cluster conv #gamma");
                if(kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMCaloConvPhotonMeson,"True reconstructed #eta, cluster conv #gamma");
            }
            // plot validated gamma gamma_{merged} FWHM if available
            if (histoTrueFWHMCaloMergedClusterMeson){
                histoTrueFWHMCaloMergedClusterMeson->Scale(1./2.35);
                DrawGammaSetMarker(histoTrueFWHMCaloMergedClusterMeson, 25, 0.8, kViolet+2, kViolet+2);
                histoTrueFWHMCaloMergedClusterMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMCaloMergedClusterMeson,"True reconstructed #pi^{0}, merged cluster #gamma");
                if(kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMCaloMergedClusterMeson,"True reconstructed #eta, merged cluster #gamma");
            }
            
            PutProcessLabelAndEnergyOnPlot(0.7, 0.94, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            legendFWHM2->Draw();
            canvasFWHM->Update(); 
            canvasFWHM->SaveAs(Form("%s/%s_FWHMAddedInfos_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
        }

        //**********************************************************************************
        //******************** FWHM Plot further decomposed for Calo + Calo *****************
        //**********************************************************************************        

        if (mode==4 || mode == 5){
            // read additional histos
            TH1D* histoTrueFWHMCaloPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueFWHMMesonCaloPhoton");
            TH1D* histoTrueFWHMCaloConvPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueFWHMMesonCaloConvPhoton");
            TH1D* histoTrueFWHMMixedCaloConvPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueFWHMMesonMixedCaloConvPhoton");

            histoFWHMMeson->GetYaxis()->SetRangeUser(minFWHM, maxFWHMAdd); 
            histoFWHMMeson->DrawCopy("e1,p"); 
            histoTrueFWHMMeson->DrawCopy("same,e1,p"); 

            // count lines to be put for legend
            Int_t nLinesLegend = 2;
            if (histoMCrecFWHMMeson) nLinesLegend++;
            if (histoTrueFWHMCaloPhotonMeson) nLinesLegend++;
            if (histoTrueFWHMCaloConvPhotonMeson) nLinesLegend++;
            if (histoTrueFWHMMixedCaloConvPhotonMeson) nLinesLegend++;
            
            Int_t nColumns                  = 1;
            if (nLinesLegend > 5) nColumns  = 2;
            
            TLegend* legendFWHM3 = GetAndSetLegend2(0.1, 0.12, 0.1+nColumns*0.35, 0.12+(0.035*nLinesLegend/(Double_t)nColumns), 0.035, nColumns, "", 42, 0.1); 
            legendFWHM3->AddEntry(histoFWHMMeson,"reconstructed Data");
            // plot MC rec FWHM if available
            if (histoMCrecFWHMMeson){
                histoMCrecFWHMMeson->DrawCopy("same,e1,p"); 
                legendFWHM3->AddEntry(histoMCrecFWHMMeson,"reconstructed MC");
            }

            if(!kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMMeson,"True reconstructed #pi^{0}");
            if(kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMMeson,"True reconstructed #eta");
            // plot validated gamma gamma FWHM if available
            if (histoTrueFWHMCaloPhotonMeson){
                histoTrueFWHMCaloPhotonMeson->Scale(1./2.35);
                DrawGammaSetMarker(histoTrueFWHMCaloPhotonMeson, 25, 0.8, kGreen+2, kGreen+2);
                histoTrueFWHMCaloPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMCaloPhotonMeson,"True reconstructed #pi^{0}, #gamma#gamma");
                if(kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMCaloPhotonMeson,"True reconstructed #eta, #gamma#gamma");
            }
            // plot validated gamma_{conv} gamma_{conv} FWHM if available
            if (histoTrueFWHMCaloConvPhotonMeson){
                histoTrueFWHMCaloConvPhotonMeson->Scale(1./2.35);
                DrawGammaSetMarker(histoTrueFWHMCaloConvPhotonMeson, 25, 0.8, kCyan+2, kCyan+2);
                histoTrueFWHMCaloConvPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMCaloConvPhotonMeson,"True reconstructed #pi^{0}, #gamma_{conv}#gamma_{conv}");
                if(kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMCaloConvPhotonMeson,"True reconstructed #eta, #gamma_{conv}#gamma_{conv}");
            }
            // plot validated gamma gamma_{conv} FWHM if available
            if (histoTrueFWHMMixedCaloConvPhotonMeson){
                histoTrueFWHMMixedCaloConvPhotonMeson->Scale(1./2.35);
                DrawGammaSetMarker(histoTrueFWHMMixedCaloConvPhotonMeson, 25, 0.8, kViolet+2, kViolet+2);
                histoTrueFWHMMixedCaloConvPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMMixedCaloConvPhotonMeson,"True reconstructed #pi^{0}, #gamma#gamma_{conv}");
                if(kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMMixedCaloConvPhotonMeson,"True reconstructed #eta, #gamma#gamma_{conv}");
            }
            
            legendFWHM3->Draw();
            PutProcessLabelAndEnergyOnPlot(0.7, 0.94, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());

            canvasFWHM->Update(); 
            canvasFWHM->SaveAs(Form("%s/%s_FWHMAddedInfos_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
            delete legendFWHM3;
        }

        delete canvasFWHM;
        delete legendFWHM;
    }    

    TH1D* histoTrueEffiPtUnmod[3]               = { NULL, NULL, NULL};
    TF1* fitEffiBiasWOWeightsPol0[3]            = { NULL, NULL, NULL}; 
    TF1* fitEffiBiasWOWeightsPol1[3]            = { NULL, NULL, NULL};
    Color_t colorEffiRatio[3]                   = { 807, kGreen+2, kCyan+2 };
    Color_t colorEffiShadePol0[3]               = { 806, kGreen-7, kCyan-7 };
    Color_t colorEffiShadePol1[3]               = { 791, kGreen-6, kCyan-6 };
    
    for (Int_t k = 0; k < 3; k++){
        histoTrueEffiPtUnmod[k]                 = (TH1D*) histoTrueEffiPt[k]->Clone(Form("histoTrueEffi%sPtUnmod", nameIntRange[k].Data())); 
        fitEffiBiasWOWeightsPol0[k]             = new TF1(Form("fitEffiBiasWOWeights%sPol0", nameIntRange[k].Data()),"[0]",0.4,maxPtMeson);
        fitEffiBiasWOWeightsPol1[k]             = new TF1(Form("fitEffiBiasWOWeights%sPol1", nameIntRange[k].Data()),"[0]/pow(x,[1])+[2]",0.4,maxPtMeson);
    }    

    TCanvas* canvasCompEffSimple = new TCanvas("canvasCompEffSimple","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasCompEffSimple, 0.10, 0.01, 0.035, 0.09);

    if (containsWOWeights &&  nameMeson.Contains("Pi0")){ //(mode!=0 || mode!=9) &&

//         if ( !(mode == 0 && optionEnergy.CompareTo("2.76TeV") == 0 )){
            TH1D* histoRatioEffWOWeightingEff[3]        = {NULL, NULL, NULL};
            TH1D* histoRatioEffWOWeightingEffCFPol0[3]  = {NULL, NULL, NULL};
            TH1D* histoRatioEffWOWeightingEffCFPol1[3]  = {NULL, NULL, NULL};
            for (Int_t k = 0; k < 3; k++){
                histoRatioEffWOWeightingEff[k]          = (TH1D*) histoEffiPt[k]->Clone();
                histoRatioEffWOWeightingEff[k]->Divide(histoRatioEffWOWeightingEff[k], histoTrueEffiPtWOWeights[k], 1., 1., "B");

                histoRatioEffWOWeightingEff[k]->Fit(fitEffiBiasWOWeightsPol0[k],"NRME+","",0.4,maxPtMeson);
                cout << WriteParameterToFile(fitEffiBiasWOWeightsPol0[k]) << endl;
                fitEffiBiasWOWeightsPol0[k]->SetLineColor(colorEffiShadePol0[k]);
                fitEffiBiasWOWeightsPol0[k]->SetLineStyle(1);

                histoRatioEffWOWeightingEffCFPol0[k]    = (TH1D*)histoRatioEffWOWeightingEff[k]->Clone(Form("histoRatioEffWOWeighting%sEffCFPol0",nameIntRange[k].Data()));
                (TVirtualFitter::GetFitter())->GetConfidenceIntervals(histoRatioEffWOWeightingEffCFPol0[k]);
                histoRatioEffWOWeightingEffCFPol0[k]->SetStats(kFALSE);
                histoRatioEffWOWeightingEffCFPol0[k]->SetFillColor(colorEffiShadePol0[k]);
                histoRatioEffWOWeightingEffCFPol0[k]->SetMarkerSize(0);

                if( !((optionEnergy.CompareTo("8TeV")==0 || optionEnergy.CompareTo("900GeV") == 0) && mode == 2 && k == 0) ) fitEffiBiasWOWeightsPol1[k]->SetParLimits(2,0.5,1.5);
                histoRatioEffWOWeightingEff[k]->Fit(fitEffiBiasWOWeightsPol1[k],"NRME+","",0.4,maxPtMeson    );
                cout << WriteParameterToFile(fitEffiBiasWOWeightsPol1[k]) << endl;
                fitEffiBiasWOWeightsPol1[k]->SetLineColor(colorEffiShadePol1[k]);
                fitEffiBiasWOWeightsPol1[k]->SetLineStyle(7);

                histoRatioEffWOWeightingEffCFPol1[k]    = (TH1D*)histoRatioEffWOWeightingEff[k]->Clone("histoRatioEffWOWeightingNormalEffCFPol1");
                (TVirtualFitter::GetFitter())->GetConfidenceIntervals(histoRatioEffWOWeightingEffCFPol1[k]);
                histoRatioEffWOWeightingEffCFPol1[k]->SetStats(kFALSE);
                histoRatioEffWOWeightingEffCFPol1[k]->SetFillColor(colorEffiShadePol1[k]);
                histoRatioEffWOWeightingEffCFPol1[k]->SetFillStyle(3003);
                histoRatioEffWOWeightingEffCFPol1[k]->SetMarkerSize(0);
                for (Int_t i=1; i< histoTrueEffiPt[k]->GetNbinsX(); i++){
                    if (histoTrueEffiPt[k]->GetBinContent(i) == 0){
                        histoRatioEffWOWeightingEffCFPol1[k]->SetBinContent(i,1);
                        histoRatioEffWOWeightingEffCFPol1[k]->SetBinContent(i,0);
                        histoRatioEffWOWeightingEffCFPol0[k]->SetBinContent(i,1);
                        histoRatioEffWOWeightingEffCFPol0[k]->SetBinContent(i,0);
                    }
                }


                    DrawAutoGammaMesonHistos( histoRatioEffWOWeightingEff[k],
                                            "", "#it{p}_{T} (GeV/#it{c})", Form("#epsilon_{eff,%s, rec}/#epsilon_{eff,%s, true wo weights} ", textMeson.Data(), textMeson.Data()),
                                            kFALSE, 1.3, 3e-6, kFALSE,
                                            kTRUE, 0.8, 1.5,
                                            kFALSE, 0., 10.);
                    DrawGammaSetMarker(histoRatioEffWOWeightingEff[k], 24, 1., colorEffiRatio[k], colorEffiRatio[k]);
                    histoRatioEffWOWeightingEff[k]->Draw("e1");
                    histoRatioEffWOWeightingEffCFPol0[k]->Draw("e3,same");
                    histoRatioEffWOWeightingEffCFPol1[k]->Draw("e3,same");
                    fitEffiBiasWOWeightsPol0[k]->Draw("same");
                    fitEffiBiasWOWeightsPol1[k]->Draw("same");
                    histoRatioEffWOWeightingEff[k]->Draw("e1,same");

                    PutProcessLabelAndEnergyOnPlot(0.72, 0.25, 28, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data(), 63, 0.03);

                canvasCompEffSimple->Update();
                canvasCompEffSimple->SaveAs(Form("%s/%s_EffiCompW0Weighting%sRatio_%s.%s",outputDir.Data(),nameMeson.Data(),nameIntRange[k].Data(),fCutSelection.Data(),suffix.Data()));

                // correct true effi to normal effi
                histoTrueEffiPt[k]->Multiply(histoTrueEffiPt[k],histoRatioEffWOWeightingEffCFPol1[k]);
            }

            // plotting of final comparison
            TH1D* histoRatioEffWOWeightingTrueEffCorr        = (TH1D*) histoEffiPt[0]->Clone();
            histoRatioEffWOWeightingTrueEffCorr->Divide(histoRatioEffWOWeightingTrueEffCorr, histoTrueEffiPt[0], 1., 1., "B");

            histoRatioEffWOWeightingEff[0]->Draw("e1");
            DrawGammaSetMarker(histoRatioEffWOWeightingTrueEffCorr, 20, 1.5, kAzure-6, kAzure-6);
            histoRatioEffWOWeightingTrueEffCorr->Draw("same,e1");

            canvasCompEffSimple->Update();
            canvasCompEffSimple->SaveAs(Form("%s/%s_EffiCompW0WeightingNormalRatioAfterFix_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
//         }
    }

    //**********************************************************************************
    //******************** Acceptance Plot *********************************************
    //**********************************************************************************
    if (kIsMC){
        // Plot simple primary acceptance
        TCanvas* canvasAcceptance = new TCanvas("canvasAcceptance2","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasAcceptance, 0.1, 0.01, 0.02, 0.10);

        Double_t rangeAcc[2]    = {0.5, 1.02};
        if (mode == 0){
            if (optionEnergy.CompareTo("pPb_5.023TeV")==0){
                rangeAcc[0]     = 0.2; 
                rangeAcc[1]     = 1.02;
                if (kIsEta) 
                    rangeAcc[0] = 0.1; 
            } else {    
                rangeAcc[0]     = 0.7; 
                rangeAcc[1]     = 1.02; 
                if (kIsEta) 
                    rangeAcc[0] = 0.5; 
            }    
        } else if (mode == 2 || mode == 4){
            rangeAcc[0]         = 0.; 
            rangeAcc[1]         = 0.3;             
        }
        
        DrawAutoGammaMesonHistos( histoAcceptance, 
                                    "", "#it{p}_{T} (GeV/#it{c})", Form("A_{%s} in |#it{y}| < %s",textMeson.Data(),rapidityRange.Data()), 
                                    kFALSE, 1.3, 3e-6, kFALSE,
                                    kTRUE, rangeAcc[0], rangeAcc[1], 
                                    kFALSE, 0., 10.);
                
        DrawGammaSetMarker(histoAcceptance, 20, 1.5, kAzure-6, kAzure-6);
        histoAcceptance->DrawCopy("e1"); 

        PutProcessLabelAndEnergyOnPlot(0.72, 0.25, 28, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data(), 63, 0.03);

        canvasAcceptance->Update();
        canvasAcceptance->SaveAs(Form("%s/%s_Acceptance_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
        
        // Plot acceptance without weights
        if (containsWOWeights){
            canvasAcceptance->cd();
            DrawGammaSetMarker(histoAcceptanceWOWeights, 24, 1., kBlack, kBlack);  
            histoAcceptanceWOWeights->DrawCopy("e1, same"); 

            TLegend* legendAccComp = GetAndSetLegend2(0.45, 0.12, 0.65, 0.22, 28);
            legendAccComp->AddEntry(histoAcceptance,"with weights");
            legendAccComp->AddEntry(histoAcceptanceWOWeights,"without weights");
            legendAccComp->Draw();
        
            
            canvasAcceptance->Update();
            canvasAcceptance->SaveAs(Form("%s/%s_AcceptanceCompWAndW0Weighting_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));

            // plot ratio of acceptances with and without weight
            TH1D* histoRatioAccWWOWeighting        =     (TH1D*) histoAcceptance->Clone(); 
            histoRatioAccWWOWeighting->Divide(histoRatioAccWWOWeighting, histoAcceptanceWOWeights, 1., 1., "B");

            DrawAutoGammaMesonHistos( histoRatioAccWWOWeighting, 
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("A_{%s, wo weights}/A_{%s, w weights} ", textMeson.Data(), textMeson.Data()), 
                                        kFALSE, 1.3, 3e-6, kFALSE,
                                        kTRUE, 0.99, 1.01, 
                                        kFALSE, 0., 10.);
            histoRatioAccWWOWeighting->Draw("e1");
            
            PutProcessLabelAndEnergyOnPlot(0.72, 0.25, 28, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data(), 63, 0.03);
            
            canvasAcceptance->Update();
            canvasAcceptance->SaveAs(Form("%s/%s_AcceptanceCompWAndW0WeightingRatio_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
    
        }

        // Plot acceptance with secondary acceptances
        if (hasNewSecQuantities){
            canvasAcceptance->cd();
            
            TLegend* legendAccSec = GetAndSetLegend2(0.45, 0.12, 0.65, 0.12+(nAccHistSec+1)*0.03, 28);
            histoAcceptance->GetYaxis()->SetRangeUser(rangeAcc[0], rangeAcc[1]*1.5);
            histoAcceptance->DrawCopy("e1");
            legendAccSec->AddEntry(histoAcceptance,"prim.");
            for (Int_t j = 0; j < 4; j++){
                if (histoSecAcceptance[j]){
                    DrawGammaSetMarker(histoSecAcceptance[j], markerStyleSec[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                    histoSecAcceptance[j]->Draw("e1,same");
                    if (!modifiedSecAcc[j]){
                        legendAccSec->AddEntry(histoSecAcceptance[j],Form("sec. from %s", nameSecMesonPlot[j].Data()));
                    } else {
                        legendAccSec->AddEntry(histoSecAcceptance[j],Form("sec. from %s, adj.", nameSecMesonPlot[j].Data()));
                    }    
                }    
            
            }
            legendAccSec->Draw();
            PutProcessLabelAndEnergyOnPlot(0.72, 0.25, 28, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data(), 63, 0.03);
            
            canvasAcceptance->Update();
            canvasAcceptance->SaveAs(Form("%s/%s_AcceptanceWithSec_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
        }    

        
        
        delete canvasAcceptance;
        
        //**********************************************************************************
        //******************** Secondary Fraction     **************************************
        //**********************************************************************************                        
        if (!kDalitz && nameMeson.Contains("Pi0")){
            cout << "Plotting fractions for secondaries" << endl;
            TCanvas* canvasSecFrac = new TCanvas("canvasSecFrac","",200,10,1350,900);  // gives the page size
            DrawGammaCanvasSettings( canvasSecFrac, 0.09, 0.018, 0.04, 0.08);

            Double_t rangeSecRatio[2]   = {0, 0.10};    
            if (mode == 0){
                rangeSecRatio[1]        = 0.025;
            } else if (mode == 2){
                rangeSecRatio[1]        = 0.05;
            } else if (mode == 4){
                rangeSecRatio[1]        = 0.05;
            }    
            
            TH2F* histo2DDummySecHad;
            histo2DDummySecHad         = new TH2F("histo2DDummySecHad","histo2DDummySecHad",1000,0, maxPtMeson,
                                                                                    100000, rangeSecRatio[0], rangeSecRatio[1]*100);
            SetStyleHistoTH2ForGraphs(histo2DDummySecHad, "#it{p}_{T} (GeV/#it{c})", "#it{r}_{X} = #frac{X->#pi^{0}}{#pi^{0}}", 0.035   ,0.04, 0.035,0.04, 0.9,1.,510,505);
            histo2DDummySecHad->GetYaxis()->SetRangeUser(0,rangeSecRatio[1]);
            histo2DDummySecHad->DrawCopy();         

            TLegend* legendSecFrac  = GetAndSetLegend2(0.65, 0.93-nSecComp*0.035*1.15, 0.94, 0.93, 0.035, 2, "", 42, 0.125);            
            for (Int_t j = 0; j < 4; j++){
                if (haveSec[j]){
                    if (histoYieldTrueSecFracMeson_orig[0][j]->GetMaximum() > 1e-4){
                        DrawGammaSetMarker(histoYieldTrueSecFracMeson_orig[0][j], markerStyleSec[j] , markerSizeSec[j], colorSec[j], colorSec[j]);
                        histoYieldTrueSecFracMeson_orig[0][j]->DrawCopy("p,e1,same");
                        legendSecFrac->AddEntry(histoYieldTrueSecFracMeson_orig[0][j],Form("#it{r}_{%s}",nameSecMesonPlot[j].Data()));
                        if (doK0SecCorrectionWithDefaultHisto == 0 && fitSecFracPLWithConst[0][j]){
                            fitSecFracPLWithConst[0][j]->SetLineColor(colorSec[j]);  
                            fitSecFracPLWithConst[0][j]->Draw("same");
                            legendSecFrac->AddEntry(fitSecFracPLWithConst[0][j],Form("fit to #it{r}_{%s}", nameSecMesonPlot[j].Data()),"l");
                        } else if (doK0SecCorrectionWithDefaultHisto == 2 && fitSecFracPurePowerlaw[0][j]) {
                            fitSecFracPurePowerlaw[0][j]->SetLineColor(colorSec[j]);  
                            fitSecFracPurePowerlaw[0][j]->Draw("same");
                            legendSecFrac->AddEntry(fitSecFracPurePowerlaw[0][j],Form("fit to #it{r}_{%s}", nameSecMesonPlot[j].Data()),"l");
                        } else if (doK0SecCorrectionWithDefaultHisto == 1 && fitDefaultSecFrac[0][j]) {    
                            fitDefaultSecFrac[0][j]->SetLineColor(colorSec[j]);  
                            fitDefaultSecFrac[0][j]->Draw("same");
                            legendSecFrac->AddEntry(fitDefaultSecFrac[0][j],Form("fit to #it{r}_{%s}", nameSecMesonPlot[j].Data()),"l");
                        } else {
                            legendSecFrac->AddEntry((TObject*)0, "","");
                        }
                    }
                }    
            }    
            legendSecFrac->Draw();
            
            PutProcessLabelAndEnergyOnPlot(0.15, 0.93, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            canvasSecFrac->Update();
            canvasSecFrac->SaveAs(Form("%s/%s_FracSecondaries_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
            
            histo2DDummySecHad->DrawCopy();         
            canvasSecFrac->cd();
            TLegend* legendSecFracUsed  = GetAndSetLegend2(0.65, 0.93-nSecCompUsed*0.035*1.15, 0.94, 0.93, 0.035, 2, "", 42, 0.125);            
            for (Int_t j = 0; j < 4; j++){
                if (haveSecUsed[j]){
                    DrawGammaSetMarker(histoYieldTrueSecFracMeson[0][j], markerStyleSec[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                    histoYieldTrueSecFracMeson[0][j]->DrawCopy("p,e1,same");  
                    legendSecFracUsed->AddEntry(histoYieldTrueSecFracMeson[0][j],Form("#it{r}_{%s}",nameSecMesonPlot[j].Data()));
                    if (doK0SecCorrectionWithDefaultHisto == 0 && fitSecFracPLWithConst[0][j]){
                        fitSecFracPLWithConst[0][j]->SetLineColor(colorSec[j]);  
                        fitSecFracPLWithConst[0][j]->Draw("same");
                        legendSecFracUsed->AddEntry(fitSecFracPLWithConst[0][j],Form("fit to #it{r}_{%s}", nameSecMesonPlot[j].Data()),"l");
                    } else if (doK0SecCorrectionWithDefaultHisto == 2 && fitSecFracPurePowerlaw[0][j]) {
                        fitSecFracPurePowerlaw[0][j]->SetLineColor(colorSec[j]);  
                        fitSecFracPurePowerlaw[0][j]->Draw("same");
                        legendSecFracUsed->AddEntry(fitSecFracPurePowerlaw[0][j],Form("fit to #it{r}_{%s}", nameSecMesonPlot[j].Data()),"l");
                    } else if (doK0SecCorrectionWithDefaultHisto == 1 && fitDefaultSecFrac[0][j]) {    
                        fitDefaultSecFrac[0][j]->SetLineColor(colorSec[j]);  
                        fitDefaultSecFrac[0][j]->Draw("same");
                        legendSecFracUsed->AddEntry(fitDefaultSecFrac[0][j],Form("fit to #it{r}_{%s}", nameSecMesonPlot[j].Data()),"l");
                    } else {
                        legendSecFracUsed->AddEntry((TObject*)0, "","");
                    }    
                }    
            }    
            legendSecFracUsed->Draw();
            
            PutProcessLabelAndEnergyOnPlot(0.15, 0.93, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            canvasSecFrac->Update();
            canvasSecFrac->SaveAs(Form("%s/%s_FracSecondariesSmoothed_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));

            // plot K0s component for different integration ranges
            canvasSecFrac->cd();
            histo2DDummySecHad->DrawCopy();        
            TLegend* legendSecFracK0s  = GetAndSetLegend2(0.65, 0.93-3*0.035*1.15, 0.94, 0.93, 0.035, 1, "", 42, 0.125);            
            for (Int_t k = 0; k < 3; k++){
                DrawGammaSetMarker(histoYieldTrueSecFracMeson[k][0], markerStyleIntRanges[k] , markerSizeIntRanges[k], colorIntRanges[k], colorIntRanges[k]);  
                histoYieldTrueSecFracMeson[k][0]->DrawCopy("p,e1,same");  
                legendSecFracK0s->AddEntry(histoYieldTrueSecFracMeson[k][0],Form("#it{r}_{K^{0}_{s}} %s",nameIntRangePlot[k].Data()));
            }    
            legendSecFracK0s->Draw();
            
            PutProcessLabelAndEnergyOnPlot(0.15, 0.93, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            canvasSecFrac->Update();
            canvasSecFrac->SaveAs(Form("%s/%s_FracSecondariesFromK0sDiffIntRanges_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));

            
            // plot all fractions in Logy scale
            canvasSecFrac->cd();
            canvasSecFrac->SetTopMargin(0.01);
            canvasSecFrac->SetLogy(1);
            histo2DDummySecHad->GetYaxis()->SetRangeUser(1e-6,rangeSecRatio[1]*10);
            histo2DDummySecHad->DrawCopy();         

            for (Int_t j = 0; j < 4; j++){
                if (haveSec[j]){
                    if (histoYieldTrueSecFracMeson_orig[0][j]->GetMaximum() > 1e-4){
                        histoYieldTrueSecFracMeson_orig[0][j]->DrawCopy("e1,same");
                        if (doK0SecCorrectionWithDefaultHisto == 0 && fitSecFracPLWithConst[0][j]){
                            fitSecFracPLWithConst[0][j]->Draw("same");
                        } else if (doK0SecCorrectionWithDefaultHisto == 2 && fitSecFracPurePowerlaw[0][j]) {
                            fitSecFracPurePowerlaw[0][j]->Draw("same");
                        } else if (doK0SecCorrectionWithDefaultHisto == 1 && fitDefaultSecFrac[0][j]) {    
                            fitDefaultSecFrac[0][j]->Draw("same");
                        } 
                    }
                }
                legendSecFrac->Draw();
            }    
            PutProcessLabelAndEnergyOnPlot(0.15, 0.93, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            canvasSecFrac->Update();
            canvasSecFrac->SaveAs(Form("%s/%s_FracSecondariesLogy_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));

            histo2DDummySecHad->DrawCopy();
            for (Int_t j = 0; j < 4; j++){
                if (haveSecUsed[j]){
                    DrawGammaSetMarker(histoYieldTrueSecFracMeson[0][j], markerStyleSec[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                    histoYieldTrueSecFracMeson[0][j]->DrawCopy("e1,same");        
                    if (doK0SecCorrectionWithDefaultHisto == 0 && fitSecFracPLWithConst[0][j]){
                        fitSecFracPLWithConst[0][j]->Draw("same");
                    } else if (doK0SecCorrectionWithDefaultHisto == 2 && fitSecFracPurePowerlaw[0][j]) {
                        fitSecFracPurePowerlaw[0][j]->Draw("same");
                    } else if (doK0SecCorrectionWithDefaultHisto == 1 && fitDefaultSecFrac[0][j]) {    
                        fitDefaultSecFrac[0][j]->Draw("same");
                    }
                }
            }
            legendSecFracUsed->Draw();
            
            PutProcessLabelAndEnergyOnPlot(0.15, 0.93, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            canvasSecFrac->Update();
            canvasSecFrac->SaveAs(Form("%s/%s_FracSecondariesSmoothedLogy_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));            
            delete canvasSecFrac;
            
        }
        
        //**********************************************************************************
        //******************** Efficiency Simple Plot **************************************
        //**********************************************************************************
        TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasEffSimple, 0.065, 0.012, 0.035, 0.09);

        TH2F* histo2DDummyEffi;
        histo2DDummyEffi         = new TH2F("histo2DDummyEffi","histo2DDummyEffi",1000,0, histoTrueEffiPtUnmod[0]->GetXaxis()->GetBinUpEdge(histoTrueEffiPtUnmod[0]->GetNbinsX()),
                                                                                1000, 0, 1.3*histoTrueEffiPtUnmod[0]->GetMaximum());
        SetStyleHistoTH2ForGraphs(histo2DDummyEffi, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{eff}", 0.035,0.04, 0.035,0.04, 0.9,0.8, 510, 505);
        histo2DDummyEffi->DrawCopy();         

            DrawGammaSetMarker(histoTrueEffiPtUnmod[0], 20, 1., kBlack, kBlack);
            histoTrueEffiPtUnmod[0]->DrawCopy("same,e1,p");
            if (containsWOWeights){
                DrawGammaSetMarker(histoEffiPt[0], 25, 1., kGreen+2, kGreen+2);
                histoEffiPt[0]->DrawCopy("same,e1,p");
                DrawGammaSetMarker(histoTrueEffiPtWOWeights[0], 24, 1., 807, 807);
                histoTrueEffiPtWOWeights[0]->DrawCopy("same,e1,p");
                DrawGammaSetMarker(histoTrueEffiPt[0], 21, 1., kBlue+1, kBlue+1);
                histoTrueEffiPt[0]->DrawCopy("same,e1,p");
            } else if (mode == 4 || mode == 2 || mode == 0){
                DrawGammaSetMarker(histoEffiPt[0], 25, 1., kGreen+2, kGreen+2);
                histoEffiPt[0]->DrawCopy("same,e1,p");            
            }    
            TLegend* legendEff = GetAndSetLegend2(0.5,0.13,0.95,0.24, 0.035, 1, "", 42, 0.1);
            if (containsWOWeights){
                legendEff->AddEntry(histoTrueEffiPtUnmod[0],"validated efficiency, w weights");
                legendEff->AddEntry(histoTrueEffiPtWOWeights[0],"validated efficiency, w/o weights"); 
                legendEff->AddEntry(histoEffiPt[0],"reconstructed efficiency, as in Data"); 
                legendEff->AddEntry(histoTrueEffiPt[0],"corr validated efficiency"); 
            } else if (mode == 4 || mode == 2 || mode == 0 ) {
                legendEff->AddEntry(histoEffiPt[0],"reconstructed efficiency, as in Data"); 
                legendEff->AddEntry(histoTrueEffiPtUnmod[0],"validated efficiency");            
            } else {
                legendEff->AddEntry(histoTrueEffiPtUnmod[0],"validated efficiency");
            }
            legendEff->Draw();

        canvasEffSimple->Update();
        PutProcessLabelAndEnergyOnPlot(0.12, 0.95, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());

        canvasEffSimple->SaveAs(Form("%s/%s_TrueEffSimple_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data())); 

        //**********************************************************************************
        //**************** Plot primary efficiency together with sec effis *****************
        //**********************************************************************************        
        if (hasNewSecQuantities){
            canvasEffSimple->cd();
            histo2DDummyEffi->DrawCopy();

            TLegend* legendEffWithSec = GetAndSetLegend2(0.12,0.8-(nEffHistSec+2)*0.035,0.4,0.8 , 0.035, 1, "", 42, 0.1);
            if (containsWOWeights){
                DrawGammaSetMarker(histoTrueEffiPt[0], 20, 1., kBlack, kBlack);
                histoTrueEffiPt[0]->DrawCopy("same,e1,p");
                DrawGammaSetMarker(histoEffiPt[0], 25, 1., kGreen+2, kGreen+2);
                histoEffiPt[0]->DrawCopy("same,e1,p");
                legendEffWithSec->AddEntry(histoTrueEffiPt[0],"val. prim");
                legendEffWithSec->AddEntry(histoEffiPt[0],"rec. prim, as in Data"); 
            } else if (mode == 4){
                DrawGammaSetMarker(histoTrueEffiPtUnmod[0], 20, 1., kBlack, kBlack);
                histoTrueEffiPtUnmod[0]->DrawCopy("same,e1,p");
                DrawGammaSetMarker(histoEffiPt[0], 25, 1., kGreen+2, kGreen+2);
                histoEffiPt[0]->DrawCopy("same,e1,p");
                legendEffWithSec->AddEntry(histoTrueEffiPtUnmod[0],"val. prim");
                legendEffWithSec->AddEntry(histoEffiPt[0],"rec. prim, as in Data"); 
            } else {
                DrawGammaSetMarker(histoTrueEffiPtUnmod[0], 20, 1., kBlack, kBlack);
                histoTrueEffiPtUnmod[0]->DrawCopy("same,e1,p");
                legendEffWithSec->AddEntry(histoTrueEffiPtUnmod[0],"val. prim");
            }
            for (Int_t j = 0; j < 4; j++){
                if (histoSecTrueEffi[0][j]){
                    DrawGammaSetMarker(histoSecTrueEffi[0][j],  markerStyleSec[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                    histoSecTrueEffi[0][j]->DrawCopy("same,e1");  
                    if (!modifiedSecTrueEffi[0][j])
                        legendEffWithSec->AddEntry(histoSecTrueEffi[0][j],Form("val. #pi^{0} from %s",nameSecMesonPlot[j].Data()),"p");
                    else 
                        legendEffWithSec->AddEntry(histoSecTrueEffi[0][j],Form("val. #pi^{0} from %s adj.",nameSecMesonPlot[j].Data()),"p");
                }
            }
            legendEffWithSec->Draw();

            canvasEffSimple->Update();
            PutProcessLabelAndEnergyOnPlot(0.12, 0.95, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());

            canvasEffSimple->SaveAs(Form("%s/%s_TrueEffWithSecEffi_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));     

            canvasEffSimple->cd();
            TH2F* histo2DDummyEffiRatio;
            histo2DDummyEffiRatio       = new TH2F("histo2DDummyEffiRatio","histo2DDummyEffiRatio",1000,0, histoTrueEffiPtUnmod[0]->GetXaxis()->GetBinUpEdge(histoTrueEffiPtUnmod[0]->GetNbinsX()),
                                                                                    1000, 0, 2);
            SetStyleHistoTH2ForGraphs(histo2DDummyEffiRatio, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{sec,eff}/#epsilon_{eff}", 0.035,0.04, 0.035,0.04, 0.9,0.8, 510, 505);
            histo2DDummyEffiRatio->DrawCopy();         

            TLegend* legendEffWithSecRatio = GetAndSetLegend2(0.12,0.8-(nEffHistSec)*0.035,0.4,0.8 , 0.035, 1, "", 42, 0.1);
            Bool_t  plotted             = kFALSE;
            for (Int_t j = 0; j < 4; j++){
                if (histoSecTrueEffi[0][j]){
                    plotted                             = kTRUE;
                    DrawGammaSetMarker(histoRatioSecEffDivTrueEff[0][j],  markerStyleSec[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                    histoRatioSecEffDivTrueEff[0][j]->DrawCopy("same,e1");  
                    legendEffWithSecRatio->AddEntry(histoRatioSecEffDivTrueEff[0][j],Form("val. #pi^{0} from %s",nameSecMesonPlot[j].Data()),"p");
                }
            }
            fithistoRatioSecEffDivTrueEff[0][0]->SetLineColor(kRed);
            fithistoRatioSecEffDivTrueEff[0][0]->SetLineWidth(5);
            fithistoRatioSecEffDivTrueEff[0][0]->Draw("same");

            legendEffWithSecRatio->Draw();
            canvasEffSimple->Update();
            PutProcessLabelAndEnergyOnPlot(0.12, 0.95, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());

            if (plotted)
                canvasEffSimple->SaveAs(Form("%s/%s_RatioSecEffiToTrueEff_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));

        }
        delete canvasEffSimple;

        //*********************************************************************************
        //********************** True Efficiency Plot ******************************************
        //*********************************************************************************
        
        TCanvas* canvasEffi = new TCanvas("canvasEffi","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasEffi, 0.065, 0.01, 0.035, 0.09);

            histo2DDummyEffi->DrawCopy();         

            TLegend* legendTrueEff = GetAndSetLegend2(0.5,0.13,0.95,0.13+3*0.035, 0.035, 1, "", 42, 0.1);
            for (Int_t k = 0; k < 3; k++){
                DrawGammaSetMarker(histoTrueEffiPtUnmod[k], markerStyleIntRanges[k], markerSizeIntRanges[k], colorIntRanges[k], colorIntRanges[k]);
                histoTrueEffiPtUnmod[k]->DrawCopy("same,e1,p");
                legendTrueEff->AddEntry(histoTrueEffiPt[k],nameIntRangePlot[k].Data());
            }
                    
            PutProcessLabelAndEnergyOnPlot(0.12, 0.95, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            
        canvasEffi->Update();
        canvasEffi->SaveAs(Form("%s/%s_%s_TrueEfficiency_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        delete legendTrueEff;
        delete canvasEffi;
    }
    
    
    //**********************************************************************************
    //*************************** MC Yield *********************************************
    //***** need to do it for MC and data in order to have it in the output file  ******
    //**********************************************************************************

    TCanvas* canvasMCYieldMeson = new TCanvas("canvasMCYieldMeson","",1350,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasMCYieldMeson, 0.13, 0.02, 0.02, 0.09);
    canvasMCYieldMeson->SetLogy();

        TH1D *histoMCYieldMesonOldBin = (TH1D*)histoInputMesonOldBinPt->Clone();
        histoMCYieldMesonOldBin->SetName("MCYield_Meson_oldBin");
        ScaleMCYield(histoMCYieldMesonOldBin,  deltaRapid,  scaling,  nEvtMC,  nameMeson ,optDalitz);
        Float_t integralMB = 0;
        integralMB = histoMCYieldMesonOldBin->Integral(histoMCYieldMesonOldBin->FindBin(8),histoMCYieldMesonOldBin->FindBin(10));
        DrawAutoGammaMesonHistos( histoMCYieldMesonOldBin, 
                                "", "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                                kFALSE, 3., 4e-10, kTRUE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
        if (histoInputMesonOldBinPtWOWeights){
            ScaleMCYield(histoInputMesonOldBinPtWOWeights,  deltaRapid,  scaling,  nEvtMC,  nameMeson ,optDalitz);
            histoInputMesonOldBinPtWOWeights->SetName("MCYield_Meson_oldBinWOWeights");
        }   
        if (histoMCInputAddedSig){
            ScaleMCYield(histoMCInputAddedSig,  deltaRapid,  scaling,  nEvtMCAddSig,  nameMeson ,optDalitz);
            histoMCInputAddedSig->SetName("MCYield_Meson_oldBin_AddedSig");
        }
        Float_t integralJetJet = 0;
        if (histoMCInputJetJetMC){
            ScaleMCYield(histoMCInputJetJetMC,  deltaRapid,  scaling,  nEvtMCJetJet,  nameMeson ,optDalitz);
            integralJetJet = histoMCInputJetJetMC->Integral(histoMCInputJetJetMC->FindBin(8),histoMCInputJetJetMC->FindBin(10));
            histoMCInputJetJetMC->SetName("MCYield_Meson_oldBin_JetJetMC");
            if (integralJetJet > 0){
                histoMCInputJetJetMC->Scale(integralMB/integralJetJet);
                cout << endl << endl << "Scaled Jet Jet MC with: " << integralMB/integralJetJet << endl << endl;
            }    
        }   

        
        if (histoMCInputWOWeightingAddedSig){
            ScaleMCYield(histoMCInputWOWeightingAddedSig,  deltaRapid,  scaling,  nEvtMCAddSig,  nameMeson ,optDalitz);
            histoMCInputWOWeightingAddedSig->SetName("MCYield_Meson_oldBinWOWeights_AddedSig");
        }   

        TH1D *histoMCYieldMeson = (TH1D*)histoInputMesonPt->Clone();
        histoMCYieldMeson->SetName("MCYield_Meson");
        ScaleMCYield(histoMCYieldMeson,  deltaRapid,  scaling,  nEvtMC,  nameMeson ,optDalitz);
        DrawAutoGammaMesonHistos(histoMCYieldMeson , 
                                "", "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                                kFALSE, 3., 4e-10, kTRUE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., histoUnCorrectedYield[0]->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield[0]->GetNbinsX()));

        TF1* fitTsallisMC;
        if (nameMeson.CompareTo("Pi0")==0 || nameMeson.CompareTo("Pi0EtaBinning")==0 ){
            fitTsallisMC= FitObject("l","fitTsallisMC","Pi0",histoMCYieldMesonOldBin,0.3,histoUnCorrectedYield[0]->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield[0]->GetNbinsX()),NULL,"QNRME+");
        } else { 
            fitTsallisMC= FitObject("l","fitTsallisMC","Eta",histoMCYieldMesonOldBin,0.3,histoUnCorrectedYield[0]->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield[0]->GetNbinsX()),NULL,"QNRME+");
        }
        DrawGammaSetMarkerTF1(fitTsallisMC, 1, 1.5, kBlue);
        TString forOutput= WriteParameterToFile(fitTsallisMC);
        cout << forOutput.Data()<< endl;
        histoMCYieldMesonOldBin->Draw("l,hist,same");
        if (histoMCInputJetJetMC){
            histoMCInputJetJetMC->SetLineColor(kRed+1);
            histoMCInputJetJetMC->Draw("l,hist,same");
        }    
        fitTsallisMC->Draw("same");

    canvasMCYieldMeson->SaveAs(Form("%s/%s_histoMCYieldMeson_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
    delete canvasMCYieldMeson; 
        
    //**********************************************************************************
    //******************** RAW Yield spectrum ******************************************
    //**********************************************************************************
        
    TCanvas* canvasRAWYield = new TCanvas("canvasRAWYield","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRAWYield, 0.10, 0.01, 0.02, 0.10); 
    canvasRAWYield->SetLogy(1);

        TH1D* histoUnCorrectedYieldDrawing = (TH1D*)histoUnCorrectedYield[0]->Clone();
        histoUnCorrectedYieldDrawing->Scale(1./nEvt);
        DrawAutoGammaMesonHistos( histoUnCorrectedYieldDrawing, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "RAW Yield/ #it{N}_{Evt}", 
                                    kTRUE, 3., 4e-10, kTRUE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
        histoUnCorrectedYieldDrawing->SetLineWidth(0.5); 
        DrawGammaSetMarker(histoUnCorrectedYieldDrawing, 20, 0.5, kBlack, kBlack);  
        histoUnCorrectedYieldDrawing->DrawCopy("e1");

        PutProcessLabelAndEnergyOnPlot(0.62, 0.95, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
        
    canvasRAWYield->Update();
    canvasRAWYield->SaveAs(Form("%s/%s_%s_RAWYieldPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
    delete canvasRAWYield;

    
     // **************************************************************************************
    // ************** Plot raw yield with differnt yield extraction methods ***********
    // **************************************************************************************
    TCanvas* canvasRawYield = new TCanvas("canvasRawYield","",1350,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasRawYield, 0.13, 0.02, 0.02, 0.09);
    canvasRawYield->SetLogy();

    TPad* padRawYieldHistos = new TPad("padRawYieldHistos", "", 0., 0.3, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padRawYieldHistos, 0.12, 0.02, 0.02, 0.);
    padRawYieldHistos->Draw();

    TPad* padRawYieldRatios = new TPad("padRawYieldRatios", "", 0., 0., 1., 0.3,-1, -1, -2);
    DrawGammaPadSettings( padRawYieldRatios, 0.12, 0.02, 0., 0.18);
    padRawYieldRatios->Draw();

    padRawYieldHistos->cd();
    padRawYieldHistos->SetLogy(); 
         
        TH2F* histo2DDummyPtRaw;
        histo2DDummyPtRaw               = new TH2F("histo2DDummyPtRaw","histo2DDummyPtRaw",1000,0, histoUnCorrectedYield[0]->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield[0]->GetNbinsX()),
                                                                                10000, 0.01*FindSmallestBin1DHist(histoUnCorrectedYield[0],1e6 ), 3*histoUnCorrectedYield[0]->GetMaximum());
        SetStyleHistoTH2ForGraphs(histo2DDummyPtRaw, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.033,0.04, 0.033,0.04, 1,1.35);
        histo2DDummyPtRaw->DrawCopy();         
        
        TLegend* legendYieldRaw = GetAndSetLegend2(0.15,0.03,0.66,0.03+6*0.035, 0.035, 1, "", 42, 0.1);
        for (Int_t k = 0; k < 6; k++){
            DrawGammaSetMarker(histoUnCorrectedYield[k], markerStyleIntRanges[k], markerSizeIntRanges[k], colorIntRanges[k], colorIntRanges[k]);
            histoUnCorrectedYield[k]->DrawCopy("e1,same"); 
            legendYieldRaw->AddEntry(histoUnCorrectedYield[k],nameIntRangePlot[k].Data());            
        }    
        legendYieldRaw->Draw();    
        PutProcessLabelAndEnergyOnPlot(0.6, 0.95, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
    
    padRawYieldRatios->cd();
    padRawYieldRatios->SetTickx();
    padRawYieldRatios->SetTicky();
    padRawYieldRatios->SetLogy(0);
    
        Double_t rangeRatioPtRaw[2]    = {0.5, 1.53};
        TH2F* histo2DDummyRatioPtRaw;
        histo2DDummyRatioPtRaw         = new TH2F("histo2DDummyRatioPtRaw","histo2DDummyRatioPtRaw",1000,0, histoUnCorrectedYield[0]->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield[0]->GetNbinsX()),
                                                                                1000, rangeRatioPtRaw[0], rangeRatioPtRaw[1]);
        SetStyleHistoTH2ForGraphs(histo2DDummyRatioPtRaw, "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.07,0.1, 0.07,0.1, 0.8,0.55, 510, 505);
        histo2DDummyRatioPtRaw->DrawCopy();         

        for (Int_t k = 0; k < 6; k++){
            DrawGammaSetMarker(RatioRaw[k], markerStyleIntRanges[k], markerSizeIntRanges[k], colorIntRanges[k], colorIntRanges[k]);
            RatioRaw[k]->DrawCopy("e1,same"); 

        }    
        DrawGammaLines(0., histoUnCorrectedYield[0]->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield[0]->GetNbinsX()),1., 1.,0.1);
        
    canvasRawYield->Update();
    canvasRawYield->SaveAs(Form("%s/%s_%s_RawYieldDiffIntRanges_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(),  fCutSelection.Data(), suffix.Data()));
       
    //***********************************************************************************************
    //*********************************** correction for yield **************************************
    //***********************************************************************************************
    TH1D* histoCorrectedYieldNorm[6]        = {NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoCorrectedYieldTrue[6]        = {NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoCorrectedYieldTrueFixed[3]   = {NULL, NULL, NULL};
    TH1D* histoCorrectedYieldFixed[3]       = {NULL, NULL, NULL};
    TH1D* histoCompleteCorr                 = (TH1D*)histoTrueEffiPt[0]->Clone();
    TH1D *RatioTrue[6]                      = {NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D *RatioNormal[6]                    = {NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D *RatioNormalToTrue                 = NULL;
    TH1D *RatioTrueMCInput                  = NULL;
    
    // corrected for resonance feed down
    TH1D* histoFeedDownCorrectedYieldNorm[6]        = {NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoFeedDownCorrectedYieldTrue[6]        = {NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoFeedDownCorrectedYieldTrueFixed[3]   = {NULL, NULL, NULL};
    TH1D* histoFeedDownCorrectedYieldFixed[3]       = {NULL, NULL, NULL};
    TH1D* FeedDownCorrectedRatioTrue[6]             = {NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* FeedDownCorrectedRatioNormal[6]           = {NULL, NULL, NULL, NULL, NULL, NULL};

    for (Int_t k = 0; k < 6; k++){
        Int_t m = k;
        histoCorrectedYieldNorm[k]      = (TH1D*)histoUnCorrectedYield[k]->Clone(Form("CorrectedYieldNormEff%s",nameIntRange[k].Data()));
        histoCorrectedYieldTrue[k]      = (TH1D*)histoUnCorrectedYield[k]->Clone(Form("CorrectedYieldTrueEff%s",nameIntRange[k].Data()));
        
        // resonance feed down correction
        histoFeedDownCorrectedYieldNorm[k]      = (TH1D*)histoUnCorrectedYield[k]->Clone(Form("FeedDownCorrectedYieldNormEff%s",nameIntRange[k].Data()));
        histoFeedDownCorrectedYieldTrue[k]      = (TH1D*)histoUnCorrectedYield[k]->Clone(Form("FeedDownCorrectedYieldTrueEff%s",nameIntRange[k].Data()));
        
        if (k < 3){
            histoCorrectedYieldFixed[k]     = (TH1D*)histoUnCorrectedYield[k]->Clone(Form("CorrectedYieldEff%sFixed",nameIntRange[k].Data()));
            histoCorrectedYieldTrueFixed[k] = (TH1D*)histoUnCorrectedYield[k]->Clone(Form("CorrectedYieldTrueEff%sFixed",nameIntRange[k].Data()));
            
            // resonance feed down correction
            histoFeedDownCorrectedYieldFixed[k]     = (TH1D*)histoUnCorrectedYield[k]->Clone(Form("FeedDownCorrectedYieldEff%sFixed",nameIntRange[k].Data()));
            histoFeedDownCorrectedYieldTrueFixed[k] = (TH1D*)histoUnCorrectedYield[k]->Clone(Form("FeedDownCorrectedYieldTrueEff%sFixed",nameIntRange[k].Data()));
        } else {    
            m = k-3;
        }
        
        if (!optDalitz){
            cout << "correcting spectra in " << nameIntRange[k].Data() << endl;
            cout << k << "\t" << k << "\t" << m << endl;
            CorrectYield(histoCorrectedYieldNorm[k], histoYieldSecMeson[k], histoYieldSecMesonFromExternalInput[k] ,histoEffiPt[m], histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
            CorrectYield(histoCorrectedYieldTrue[k], histoYieldSecMeson[k], histoYieldSecMesonFromExternalInput[k], histoTrueEffiPt[m], histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
            
            // resonance feed down correction
            CorrectYieldInclResonanceFeedDown(histoFeedDownCorrectedYieldNorm[k], histoYieldSecMeson[k], histoYieldSecMesonFromExternalInput[k], histoYieldResonanceFeedDownPi0FromExternalInput[k], histoEffiPt[m], histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
            CorrectYieldInclResonanceFeedDown(histoFeedDownCorrectedYieldTrue[k], histoYieldSecMeson[k], histoYieldSecMesonFromExternalInput[k], histoYieldResonanceFeedDownPi0FromExternalInput[k], histoTrueEffiPt[m], histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
            
            if (k < 3){
                CorrectYield(histoCorrectedYieldTrueFixed[k], histoYieldSecMeson[k], histoYieldSecMesonFromExternalInput[k], histoTrueEffiPtFixed[k], histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
                CorrectYield(histoCorrectedYieldFixed[k], histoYieldSecMeson[k], histoYieldSecMesonFromExternalInput[k], histoEffiPtFixed[k], histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
                
                // resonance feed down correction
                CorrectYieldInclResonanceFeedDown(histoFeedDownCorrectedYieldTrueFixed[k], histoYieldSecMeson[k], histoYieldSecMesonFromExternalInput[k], histoYieldResonanceFeedDownPi0FromExternalInput[k], histoTrueEffiPtFixed[k], histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
                CorrectYieldInclResonanceFeedDown(histoFeedDownCorrectedYieldFixed[k], histoYieldSecMeson[k], histoYieldSecMesonFromExternalInput[k], histoYieldResonanceFeedDownPi0FromExternalInput[k], histoEffiPtFixed[k], histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
            }    
        } else {
            CorrectYieldDalitz(histoCorrectedYieldNorm[k], histoYieldGGMeson[k], histoEffiPt[k], histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
            CorrectYieldDalitz(histoCorrectedYieldTrue[k], histoYieldGGMeson[k], histoTrueEffiPt[k], histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
        }    

        RatioTrue[k]                    = (TH1D*) histoCorrectedYieldTrue[k]->Clone(); 
        RatioTrue[k]->Divide(RatioTrue[k],histoCorrectedYieldTrue[0],1.,1.,"B");
        RatioNormal[k]                  = (TH1D*) histoCorrectedYieldNorm[k]->Clone(); 
        RatioNormal[k]->Divide(RatioNormal[k],histoCorrectedYieldNorm[0],1.,1.,"B");
        
        // resonance feed down corrected
        FeedDownCorrectedRatioTrue[k]   = (TH1D*) histoFeedDownCorrectedYieldTrue[k]->Clone();
        FeedDownCorrectedRatioTrue[k]->Divide(FeedDownCorrectedRatioTrue[k],histoFeedDownCorrectedYieldTrue[0],1.,1.,"B");
        FeedDownCorrectedRatioNormal[k] = (TH1D*) histoFeedDownCorrectedYieldNorm[k]->Clone();
        FeedDownCorrectedRatioNormal[k]->Divide(FeedDownCorrectedRatioNormal[k],histoFeedDownCorrectedYieldNorm[0],1.,1.,"B");
        
        if (k == 0){
            RatioTrueMCInput            = (TH1D*) histoMCYieldMeson->Clone();
            RatioTrueMCInput->Divide(RatioTrueMCInput,histoCorrectedYieldTrue[k],1.,1.,"B");
            RatioNormalToTrue           = (TH1D*) histoCorrectedYieldNorm[k]->Clone();  
            RatioNormalToTrue->Divide(RatioNormalToTrue,histoCorrectedYieldTrue[k],1.,1.,"B");
        }    
    }
    
    if (!optDalitz){
        CompileFullCorrectionFactor( histoCompleteCorr, histoAcceptance, deltaRapid);    
    }
    
//     return;
    // **************************************************************************************
    // ************** Plot corrected yield with differnt yield extraction methods ***********
    // **************************************************************************************
    TCanvas* canvasCorrectedYield = new TCanvas("canvasCorrectedYield","",1350,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasCorrectedYield, 0.13, 0.02, 0.02, 0.09);
    canvasCorrectedYield->SetLogy();

    TPad* padCorrectedYieldHistos = new TPad("padCorrectedYieldHistos", "", 0., 0.3, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padCorrectedYieldHistos, 0.12, 0.02, 0.02, 0.);
    padCorrectedYieldHistos->Draw();

    TPad* padCorrectedYieldRatios = new TPad("padCorrectedYieldRatios", "", 0., 0., 1., 0.3,-1, -1, -2);
    DrawGammaPadSettings( padCorrectedYieldRatios, 0.12, 0.02, 0., 0.18);
    padCorrectedYieldRatios->Draw();

    padCorrectedYieldHistos->cd();
    padCorrectedYieldHistos->SetLogy(); 
    
        TH2F* histo2DDummyPt;
        histo2DDummyPt               = new TH2F("histo2DDummyPt","histo2DDummyPt",1000,0, histoCorrectedYieldTrue[0]->GetXaxis()->GetBinUpEdge(histoCorrectedYieldTrue[0]->GetNbinsX()),
                                                                                10000, 0.01*FindSmallestBin1DHist(histoCorrectedYieldTrue[0],1e6 ), 3*histoCorrectedYieldTrue[0]->GetMaximum());
        SetStyleHistoTH2ForGraphs(histo2DDummyPt, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.033,0.04, 0.033,0.04, 1,1.35);
        histo2DDummyPt->DrawCopy();         
        
        TLegend* legendYield3 = GetAndSetLegend2(0.15,0.03,0.66,0.03+6*0.035, 0.035, 1, "", 42, 0.1);
        for (Int_t k = 0; k < 6; k++){
            DrawGammaSetMarker(histoCorrectedYieldTrue[k], markerStyleIntRanges[k], markerSizeIntRanges[k], colorIntRanges[k], colorIntRanges[k]);
            histoCorrectedYieldTrue[k]->DrawCopy("e1,same"); 
            legendYield3->AddEntry(histoCorrectedYieldTrue[k],nameIntRangePlot[k].Data());
        }    
        legendYield3->Draw();    
        PutProcessLabelAndEnergyOnPlot(0.6, 0.95, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
    
    padCorrectedYieldRatios->cd();
    padCorrectedYieldRatios->SetTickx();
    padCorrectedYieldRatios->SetTicky();
    padCorrectedYieldRatios->SetLogy(0);
    
        Double_t rangeRatioPt[2]    = {0.8, 1.23};
//         if (kIsEta){
//             rangeRatioPt[0]         = 0.4
//             rangeRatioPt[1]         = 1.63;
//         }    
        TH2F* histo2DDummyRatioPt;
        histo2DDummyRatioPt         = new TH2F("histo2DDummyRatioPt","histo2DDummyRatioPt",1000,0, histoCorrectedYieldTrue[0]->GetXaxis()->GetBinUpEdge(histoCorrectedYieldTrue[0]->GetNbinsX()),
                                                                                1000, rangeRatioPt[0], rangeRatioPt[1]);
        SetStyleHistoTH2ForGraphs(histo2DDummyRatioPt, "#it{p}_{T} (GeV/#it{c})", "#frac{modified}{standard}", 0.07,0.1, 0.07,0.1, 0.8,0.55, 510, 505);
        histo2DDummyRatioPt->DrawCopy();         

        for (Int_t k = 0; k < 6; k++){
            DrawGammaSetMarker(RatioTrue[k], markerStyleIntRanges[k], markerSizeIntRanges[k], colorIntRanges[k], colorIntRanges[k]);
            RatioTrue[k]->DrawCopy("e1,same"); 
        }    
        DrawGammaLines(0., histoCorrectedYieldTrue[0]->GetXaxis()->GetBinUpEdge(histoCorrectedYieldTrue[0]->GetNbinsX()),1., 1.,0.1);
        
    canvasCorrectedYield->Update();
    canvasCorrectedYield->SaveAs(Form("%s/%s_%s_CorrectedYieldTrueEff_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(),  fCutSelection.Data(), suffix.Data()));

    padCorrectedYieldHistos->cd();
    
        histo2DDummyPt->DrawCopy();         
        
        for (Int_t k = 0; k < 6; k++){
            DrawGammaSetMarker(histoCorrectedYieldNorm[k], markerStyleIntRanges[k], markerSizeIntRanges[k], colorIntRanges[k], colorIntRanges[k]);
            histoCorrectedYieldNorm[k]->DrawCopy("e1,same"); 
        }    
        legendYield3->Draw();    
        PutProcessLabelAndEnergyOnPlot(0.6, 0.95, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
    
    padCorrectedYieldRatios->cd();
    
        histo2DDummyRatioPt->DrawCopy();         

        for (Int_t k = 0; k < 6; k++){
            DrawGammaSetMarker(RatioNormal[k], markerStyleIntRanges[k], markerSizeIntRanges[k], colorIntRanges[k], colorIntRanges[k]);
            RatioNormal[k]->DrawCopy("e1,same"); 
        }    
        DrawGammaLines(0., histoCorrectedYieldNorm[0]->GetXaxis()->GetBinUpEdge(histoCorrectedYieldNorm[0]->GetNbinsX()),1., 1.,0.1);
        
    canvasCorrectedYield->Update();
    canvasCorrectedYield->SaveAs(Form("%s/%s_%s_CorrectedYieldNormalEff_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(),  fCutSelection.Data(), suffix.Data()));
    
    //***********************************************************************************************
    //*************** Plot comparison to resonance feed down corrected yield ************************
    //***********************************************************************************************
    TH1D* ratioFeedDownCorrectedYieldTrueToStandard         = NULL;
    if (histoFeedDownCorrectedYieldTrue[0]) {
        ratioFeedDownCorrectedYieldTrueToStandard           = (TH1D*)histoFeedDownCorrectedYieldTrue[0]->Clone("RatioFeedDownCorrectedYieldTrueEffToCorrectedYieldTrueEff");
        ratioFeedDownCorrectedYieldTrueToStandard->Sumw2();
        ratioFeedDownCorrectedYieldTrueToStandard->Divide(ratioFeedDownCorrectedYieldTrueToStandard,histoCorrectedYieldTrue[0],1.,1.,"B");

        TCanvas* canvasFeedDownComparison                   = new TCanvas("canvasFeedDownComparison","",1350,1500);
        DrawGammaCanvasSettings( canvasFeedDownComparison, 0.13, 0.02, 0.02, 0.09);
        canvasCorrectedYield->SetLogy();

        TPad* padComparisonHistos = new TPad("padComparisonHistos", "", 0., 0.3, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padComparisonHistos, 0.12, 0.02, 0.02, 0.);
        padComparisonHistos->Draw();

        TPad* padComparisonRatio = new TPad("padComparisonRatio", "", 0., 0., 1., 0.3,-1, -1, -2);
        DrawGammaPadSettings( padComparisonRatio, 0.12, 0.02, 0., 0.18);
        padComparisonRatio->Draw();

        padComparisonHistos->cd();
        padComparisonHistos->SetLogy();

        TH2F* histo2DComparisonDummy                        = new TH2F("histo2DComparisonDummy", "histo2DComparisonDummy", 1000, 0, histoCorrectedYieldTrue[0]->GetXaxis()->GetBinUpEdge(histoCorrectedYieldTrue[0]->GetNbinsX()), 10000, 0.01*FindSmallestBin1DHist(histoCorrectedYieldTrue[0],1e6), 3*histoCorrectedYieldTrue[0]->GetMaximum());
        SetStyleHistoTH2ForGraphs(histo2DComparisonDummy, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.033, 0.04, 0.033,0.04, 1, 1.35);
        histo2DComparisonDummy->DrawCopy();

        TLegend* legendComparison                           = GetAndSetLegend2(0.15,0.03,0.66,0.03+2*0.035, 0.035, 1, "", 42, 0.1);

        DrawGammaSetMarker(histoCorrectedYieldTrue[0], 20, 2.0, kBlack, kBlack);
        histoCorrectedYieldTrue[0]->DrawCopy("e1,same");
        legendComparison->AddEntry(histoCorrectedYieldTrue[0],"corr. yield true eff");

        DrawGammaSetMarker(histoFeedDownCorrectedYieldTrue[0], 34, 2.0, kBlue-3, kBlue-3);
        histoFeedDownCorrectedYieldTrue[0]->DrawCopy("e1,same");
        legendComparison->AddEntry(histoFeedDownCorrectedYieldTrue[0],"corr. yield true eff, res. feed down corr.");

        legendComparison->Draw();
        PutProcessLabelAndEnergyOnPlot(0.7, 0.95, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());

        padComparisonRatio->cd();
        padComparisonRatio->SetTickx();
        padComparisonRatio->SetTicky();
        padComparisonRatio->SetLogy(0);

        Double_t minYComparisonRatio                        = 0.65;
        Double_t maxYComparisonRatio                        = 1.05;

        TH2F* histo2DComparisonDummyRatio                   = new TH2F("histo2DComparisonDummyRatio","histo2DComparisonDummyRatio",1000,0, histoCorrectedYieldTrue[0]->GetXaxis()->GetBinUpEdge(histoCorrectedYieldTrue[0]->GetNbinsX()), 1000, minYComparisonRatio, maxYComparisonRatio);
        SetStyleHistoTH2ForGraphs(histo2DComparisonDummyRatio, "#it{p}_{T} (GeV/#it{c})", "#frac{feed down. corr.}{standard}", 0.07,0.1, 0.07,0.1, 0.8,0.55, 510, 505);
        histo2DComparisonDummyRatio->DrawCopy();

        DrawGammaLines(0., histoCorrectedYieldTrue[0]->GetXaxis()->GetBinUpEdge(histoCorrectedYieldTrue[0]->GetNbinsX()), 1., 1., 2.0, kGray+2, 2);
        DrawGammaLines(0., histoCorrectedYieldTrue[0]->GetXaxis()->GetBinUpEdge(histoCorrectedYieldTrue[0]->GetNbinsX()), 0.9, 0.9, 2.0, kGray+2, 8);
        DrawGammaSetMarker(ratioFeedDownCorrectedYieldTrueToStandard, 34, 2.0, kBlue-3, kBlue-3);
        ratioFeedDownCorrectedYieldTrueToStandard->DrawCopy("e1,same");

        canvasFeedDownComparison->SaveAs(Form("%s/%s_%s_CorrectedYieldTrueEffComparisonToResonanceFeedDownCorrected_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(),  fCutSelection.Data(), suffix.Data()));

        delete histo2DComparisonDummy;
        delete histo2DComparisonDummyRatio;

        delete padComparisonHistos;
        delete padComparisonRatio;
        delete canvasFeedDownComparison;
    }

    // **************************************************************************************
    // ************** Plot corrected yield with differnt efficiencies & MC yield ************
    // **************************** Sanity check for MC *************************************
    // **************************************************************************************
    
    if (kIsMC){
        canvasCorrectedYield->cd();

        padCorrectedYieldHistos->cd();
        padCorrectedYieldHistos->SetLogy(); 

            histo2DDummyPt->DrawCopy();         
            DrawGammaSetMarker(histoCorrectedYieldTrue[0], 20, 1., kBlack, kBlack);  
            histoCorrectedYieldTrue[0]->DrawCopy("e1");  
            DrawGammaSetMarker(histoMCYieldMeson, 24, 1., kRed+2, kRed+2);  
            histoMCYieldMeson->DrawCopy("e1,same"); 
            DrawGammaSetMarker(histoCorrectedYieldNorm[0], 24, 1., kGreen+2, kGreen+2); 
            histoCorrectedYieldNorm[0]->DrawCopy("e1,same"); 
        
            cout << "here" << endl; 
            TLegend* legendYield4 = GetAndSetLegend2(0.15,0.03,0.66,0.03+3*0.035, 0.035, 1, "", 42, 0.1);
            legendYield4->AddEntry(histoCorrectedYieldTrue[0],"corr true eff");
            legendYield4->AddEntry(histoMCYieldMeson,"MC input (possibly weighted)");
            legendYield4->AddEntry(histoCorrectedYieldNorm[0],"normal eff");
            legendYield4->Draw();

            PutProcessLabelAndEnergyOnPlot(0.6, 0.95, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
        padCorrectedYieldRatios->cd();
        padCorrectedYieldRatios->SetTickx();
        padCorrectedYieldRatios->SetTicky();
        padCorrectedYieldRatios->SetLogy(0);
        
            histo2DDummyRatioPt->DrawCopy();         
            
            DrawGammaSetMarker(RatioTrue[0], 20, 1., kBlack, kBlack);
            for(Int_t b = 0; b< RatioTrue[0]->GetNbinsX(); b++){
                RatioTrue[0]->SetBinError(b+1,histoCorrectedYieldTrue[0]->GetBinError(b+1)/histoCorrectedYieldTrue[0]->GetBinContent(b+1));
            }
            RatioTrue[0]->SetFillColor(kGray+2);
            RatioTrue[0]->SetFillStyle(1);
            
            RatioTrue[0]->DrawCopy("p,e2");  
            DrawGammaSetMarker(RatioTrueMCInput, 24, 1., kRed+2, kRed+2);
            RatioTrueMCInput->DrawCopy("e1,same"); 
            DrawGammaSetMarker(RatioNormalToTrue, 25, 1., kGreen+2, kGreen+2); 
            RatioNormalToTrue->DrawCopy("e1,same"); 
            DrawGammaLines(0., histoCorrectedYieldTrue[0]->GetXaxis()->GetBinUpEdge(histoCorrectedYieldTrue[0]->GetNbinsX()),1., 1.,0.1);
            
        canvasCorrectedYield->Update();
        canvasCorrectedYield->SaveAs(Form("%s/%s_%s_CorrectedYield_SanityCheck_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(),  fCutSelection.Data(), suffix.Data()));
    }

    delete canvasCorrectedYield;
    delete legendYield3;

    //***********************************************************************************************
    //***************************  Secondary RAW Yield  *********************************************
    //***********************************************************************************************
    if (!optDalitz && nameMeson.Contains("Pi0") ){ //&& !kIsMC
        TCanvas* canvasRAWYieldSec  = new TCanvas("canvasRAWYieldSec","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRAWYieldSec, 0.10, 0.01, 0.02, 0.08); 
        canvasRAWYieldSec->SetLogy(1);
        
        Int_t nColumnsSec           = 1;
        Int_t nSecPlot              = nSecCompUsed;
//         if (!kIsMC){
            nColumnsSec             = 2;
            nSecPlot                = 4;
//         }
        
        TLegend* legendSecRAWYield  = GetAndSetLegend2(0.65,0.93-(nSecPlot+1)*0.035,0.93,0.93, 0.035, nColumnsSec, "", 42, 0.12);
        
        Double_t minRawYieldDrawing = 1e6;
        for (Int_t iPt = 1; iPt < histoUnCorrectedYieldDrawing->GetNbinsX()+1; iPt++){
             if ( histoUnCorrectedYieldDrawing->GetBinContent(iPt) > 0 && minRawYieldDrawing > histoUnCorrectedYieldDrawing->GetBinContent(iPt)){
                minRawYieldDrawing  = histoUnCorrectedYieldDrawing->GetBinContent(iPt);
            }        
        }    
        minRawYieldDrawing          = minRawYieldDrawing*1e-4;
        
        histoUnCorrectedYieldDrawing->GetYaxis()->SetRangeUser(minRawYieldDrawing, histoUnCorrectedYieldDrawing->GetMaximum());
        DrawGammaSetMarker(histoUnCorrectedYieldDrawing, 20, 1.7, kBlack, kBlack);
        histoUnCorrectedYieldDrawing->Draw("e1");
        legendSecRAWYield->AddEntry(histoUnCorrectedYieldDrawing,"RAW yield");
//         if (!kIsMC) legendSecRAWYield->AddEntry((TObject*)0,"","");
        legendSecRAWYield->AddEntry((TObject*)0,"","");
        
        for (Int_t j = 0; j < 4; j++){
            if (histoYieldSecMeson[0][j] && haveSecUsed[j]){
                histoYieldSecMeson[0][j]->Scale(1./nEvt);
                DrawGammaSetMarker(histoYieldSecMeson[0][j],  markerStyleSecWithToy[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                histoYieldSecMeson[0][j]->DrawCopy("same,e1");  
                legendSecRAWYield->AddEntry(histoYieldSecMeson[0][j],Form("#pi^{0} from %s",nameSecMesonPlot[j].Data()),"p");
            } else if ( j < 3){ //!kIsMC &&
                if (histoYieldSecMesonFromExternalInput[0][j] && histoYieldSecMesonFromExternalInput[0][j]->GetEntries()){
//                     if (!kIsMC) legendSecRAWYield->AddEntry((TObject*)0,Form("#pi^{0} from %s",nameSecMesonPlot[j].Data()),"");
                    legendSecRAWYield->AddEntry((TObject*)0,Form("#pi^{0} from %s",nameSecMesonPlot[j].Data()),"");
                }    
            }    
            if (j < 3){
                if (histoYieldSecMesonFromExternalInput[0][j] ){ //&& !kIsMC 
                    if (histoYieldSecMesonFromExternalInput[0][j]->GetEntries() > 0){
                        cout << "plotting toy approach" << endl;
                        DrawGammaSetMarker(histoYieldSecMesonFromExternalInput[0][j],  markerStyleSecFromToy[j] , markerSizeSecFromToy[j], colorSecFromToy[j], colorSecFromToy[j]);  
                        histoYieldSecMesonFromExternalInput[0][j]->DrawCopy("same,e1");  
                        legendSecRAWYield->AddEntry(histoYieldSecMesonFromExternalInput[0][j],"Toy appr.","p");
                    } else if (histoYieldSecMeson[0][j] &&  haveSecUsed[j]){ //!kIsMC &&
                        legendSecRAWYield->AddEntry((TObject*)0,"","");
                    } 
                } else if (histoYieldSecMeson[0][j] && haveSecUsed[j]){ //&& !kIsMC 
                    legendSecRAWYield->AddEntry((TObject*)0,"","");
                } 
            }     
        }
        legendSecRAWYield->Draw();

        PutProcessLabelAndEnergyOnPlot(0.15, 0.23, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
        
        canvasRAWYieldSec->Update();
        canvasRAWYieldSec->SaveAs(Form("%s/%s_%s_RAWYieldSecPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        

        if ( (histoYieldSecMesonFromExternalInput[0][0] && !kIsMC) ){
            TCanvas* canvasSecFrac2 = new TCanvas("canvasSecFrac","",200,10,1350,900);  // gives the page size
            DrawGammaCanvasSettings( canvasSecFrac2, 0.09, 0.018, 0.04, 0.08);

            Double_t rangeSecRatio[2]   = {0, 0.30};    
            if (mode == 0){
                rangeSecRatio[1]        = 0.05;
            } else if (mode == 2){
                rangeSecRatio[1]        = 0.05;
            } else if (mode == 4){
                rangeSecRatio[1]        = 0.075;
            }    
            
            TH2F* histo2DDummySecHad2;
            histo2DDummySecHad2         = new TH2F("histo2DDummySecHad2","histo2DDummySecHad2",1000,0, maxPtMeson,
                                                                                    100000, rangeSecRatio[0], rangeSecRatio[1]*100);
            SetStyleHistoTH2ForGraphs(histo2DDummySecHad2, "#it{p}_{T} (GeV/#it{c})", "#it{r}_{X} = #frac{X->#pi^{0}}{#pi^{0}}", 0.035   ,0.04, 0.035,0.04, 0.9,1.,510,505);
            histo2DDummySecHad2->GetYaxis()->SetRangeUser(0,rangeSecRatio[1]);
            histo2DDummySecHad2->DrawCopy();         

            TLegend* legendSecRAWRatio  = GetAndSetLegend2(0.65,0.93-(nSecPlot)*0.035,0.93,0.93, 0.035, nColumnsSec, "", 42, 0.12);
            for (Int_t j = 0; j < 4; j++){
                if (histoRatioYieldSecMeson[0][j] && haveSecUsed[j]){
                    DrawGammaSetMarker(histoRatioYieldSecMeson[0][j],  markerStyleSecWithToy[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                    histoRatioYieldSecMeson[0][j]->DrawCopy("same,e1");  
                    legendSecRAWRatio->AddEntry(histoRatioYieldSecMeson[0][j],Form("#pi^{0} from %s",nameSecMesonPlot[j].Data()),"p");
                } else if (!kIsMC && j < 3){
                    if (histoRatioYieldSecMesonFromExtInput[0][j] && histoRatioYieldSecMesonFromExtInput[0][j]->GetEntries()){
                        if (!kIsMC) legendSecRAWRatio->AddEntry((TObject*)0,Form("#pi^{0} from %s",nameSecMesonPlot[j].Data()),"");
                    }    
                }    
                if (j < 3){
                    if (histoRatioYieldSecMesonFromExtInput[0][j] && !kIsMC ){
                        if (histoRatioYieldSecMesonFromExtInput[0][j]->GetEntries() > 0){
                            DrawGammaSetMarker(histoRatioYieldSecMesonFromExtInput[0][j],  markerStyleSecFromToy[j] , markerSizeSecFromToy[j], colorSecFromToy[j], colorSecFromToy[j]);  
                            histoRatioYieldSecMesonFromExtInput[0][j]->DrawCopy("same,e1");  
                            legendSecRAWRatio->AddEntry(histoRatioYieldSecMesonFromExtInput[0][j],"Toy appr.","p");
                        } else if (histoRatioYieldSecMeson[0][j] && !kIsMC && haveSecUsed[j]){
                            legendSecRAWRatio->AddEntry((TObject*)0,"","");
                        } 
                    } else if (histoRatioYieldSecMeson[0][j] && !kIsMC && haveSecUsed[j]){
                        legendSecRAWRatio->AddEntry((TObject*)0,"","");
                    } 
                }     
            }
            legendSecRAWRatio->Draw();
            PutProcessLabelAndEnergyOnPlot(0.15, 0.93, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            
            canvasSecFrac2->Update();
            canvasSecFrac2->SaveAs(Form("%s/%s_%s_EffectiveSecCorrPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        }

        delete canvasRAWYieldSec;
        
        // resonance feed down yield plot
        TCanvas* canvasFeedDownYield  = new TCanvas("canvasFeedDownYield","",200,10,1350,900);
        DrawGammaCanvasSettings( canvasFeedDownYield, 0.10, 0.01, 0.02, 0.08);
        canvasFeedDownYield->SetLogy(1);
        
        nColumnsSec                 = 2;
        nSecPlot                    = 4;
        
        TLegend* legendFeedDownYield  = GetAndSetLegend2(0.65,0.93-(nSecPlot+1)*0.035,0.93,0.93, 0.035, nColumnsSec, "", 42, 0.12);
        
        minRawYieldDrawing          = 1e6;
        for (Int_t iPt = 1; iPt < histoUnCorrectedYieldDrawing->GetNbinsX()+1; iPt++){
            if ( histoUnCorrectedYieldDrawing->GetBinContent(iPt) > 0 && minRawYieldDrawing > histoUnCorrectedYieldDrawing->GetBinContent(iPt)){
                minRawYieldDrawing  = histoUnCorrectedYieldDrawing->GetBinContent(iPt);
            }
        }
        minRawYieldDrawing          = minRawYieldDrawing*1e-4;
        
        histoUnCorrectedYieldDrawing->GetYaxis()->SetRangeUser(minRawYieldDrawing, histoUnCorrectedYieldDrawing->GetMaximum());
        DrawGammaSetMarker(histoUnCorrectedYieldDrawing, 20, 1.7, kBlack, kBlack);
        histoUnCorrectedYieldDrawing->Draw("e1");
        legendFeedDownYield->AddEntry(histoUnCorrectedYieldDrawing,"RAW yield");
        legendFeedDownYield->AddEntry((TObject*)0,"","");
        
        for (Int_t j = 0; j < 8; j++){
            if (histoYieldResonanceFeedDownPi0FromExternalInput[0][j] && histoYieldResonanceFeedDownPi0FromExternalInput[0][j]->GetEntries()) {
                histoYieldResonanceFeedDownPi0FromExternalInput[0][j]->Multiply(histoAcceptance);
                histoYieldResonanceFeedDownPi0FromExternalInput[0][j]->Multiply(histoTrueEffiPt[0]);
                DrawGammaSetMarker(histoYieldResonanceFeedDownPi0FromExternalInput[0][j],  markerStyleFeedDown[j] , markerSizeFeedDown[j], colorFeedDown[j], colorFeedDown[j]);
                histoYieldResonanceFeedDownPi0FromExternalInput[0][j]->DrawCopy("same,e1");
                legendFeedDownYield->AddEntry(histoYieldResonanceFeedDownPi0FromExternalInput[0][j],Form("#pi^{0} from %s",nameResMesonPlot[j].Data()),"p");
            }
        }
        legendFeedDownYield->Draw();
        
        PutProcessLabelAndEnergyOnPlot(0.15, 0.23, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
        
        canvasFeedDownYield->Update();
        canvasFeedDownYield->SaveAs(Form("%s/%s_%s_ResonanceFeedDownYieldPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        
    } else if (optDalitz){
        TCanvas* canvasRAWYieldSec  = new TCanvas("canvasRAWYieldSec","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRAWYieldSec, 0.10, 0.01, 0.02, 0.08); 
        canvasRAWYieldSec->SetLogy(1);

            DrawGammaSetMarker(histoUnCorrectedYieldDrawing, 20, 1., kBlack, kBlack);
            histoUnCorrectedYieldDrawing->Draw("e1");
            histoYieldGGMeson[0]->Scale(1./nEvt);
            DrawGammaSetMarker(histoYieldGGMeson[0], 20, 1., kBlue, kBlue); 
            histoYieldGGMeson[0]->DrawCopy("same,e1");
                
            TLegend* legendSecRAWYield  = GetAndSetLegend2(0.6,0.93-2*0.035,0.93,0.93, 0.035, 1, "", 42, 0.1);
            legendSecRAWYield->AddEntry(histoUnCorrectedYieldDrawing,"RAW yield");
            legendSecRAWYield->AddEntry(histoYieldGGMeson[0],"Contamination from #gamma#gamma");
            legendSecRAWYield->Draw();

            PutProcessLabelAndEnergyOnPlot(0.62, 0.78, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
        
        canvasRAWYieldSec->Update();
        canvasRAWYieldSec->SaveAs(Form("%s/%s_%s_RAWYieldContGGPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        delete canvasRAWYieldSec;
    }

    // *******************************************************************************************
    // ****** Show fractions of cluster origin in MC to total for Calo related analysis path *****
    // *******************************************************************************************
    if ((mode==2 || mode == 3) && kIsMC ){
        TH1D* histoTrueTotalRecYield            = (TH1D*)fileCorrections->Get("histoYieldTrueMesonFixedWindow");
        if (histoTrueTotalRecYield){
            TH1D* histoTrueTotalRecYieldGamma       = (TH1D*)fileCorrections->Get("histoYieldTrueMesonGammaFixedWindow");
            TH1D* histoTrueTotalRecYieldConvGamma   = (TH1D*)fileCorrections->Get("histoYieldTrueMesonGammaConvGammaFixedWindow");
            
            TH1D* ratioTrueGammaDivTotal            = NULL;
            TH1D* ratioTrueConvGammaDivTotal        = NULL;
            if (histoTrueTotalRecYieldGamma){
                ratioTrueGammaDivTotal              = (TH1D*)histoTrueTotalRecYieldGamma->Clone("ratioTrueGammaDivTotal");
                ratioTrueGammaDivTotal->Divide(ratioTrueGammaDivTotal,histoTrueTotalRecYield,1.,1.,"B");
            }
            if (histoTrueTotalRecYieldConvGamma){
                ratioTrueConvGammaDivTotal          = (TH1D*)histoTrueTotalRecYieldConvGamma->Clone("ratioTrueConvGammaDivTotal");
                ratioTrueConvGammaDivTotal->Divide(ratioTrueConvGammaDivTotal,histoTrueTotalRecYield,1.,1.,"B");
            }
            
            if (histoTrueTotalRecYieldGamma){
                TCanvas* canvasFracDifferentContrib = new TCanvas("canvasFracDifferentContrib","",200,10,1350,900);  // gives the page size
                DrawGammaCanvasSettings( canvasFracDifferentContrib, 0.09, 0.02, 0.04, 0.09);
                            
                DrawAutoGammaMesonHistos( ratioTrueGammaDivTotal, 
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("#frac{%s->XX}{%s->ALL}", textProcess.Data(), textProcess.Data() ), 
                                        kFALSE, 1.5, 0, kFALSE,
                                        kTRUE, 0., 1.2, 
                                        kFALSE, 0., 10.);
                ratioTrueGammaDivTotal->GetYaxis()->SetTitleOffset(0.9);
                DrawGammaSetMarker(ratioTrueGammaDivTotal, 20, 1., kBlack, kBlack);  
                DrawGammaSetMarker(ratioTrueConvGammaDivTotal, 20, 1., kBlue, kBlue);
                ratioTrueGammaDivTotal->DrawCopy("e1");  
                ratioTrueConvGammaDivTotal->DrawCopy("e1,same");  
                
                TLegend* legendFracDifferentContrib = GetAndSetLegend2(0.6,0.93-2*0.035*1.1,0.93,0.93, 0.035, 1, "", 42, 0.1);
                legendFracDifferentContrib->AddEntry(ratioTrueGammaDivTotal,"XX= #gamma_{PCM}, #gamma_{cluster}");
                legendFracDifferentContrib->AddEntry(ratioTrueConvGammaDivTotal,"XX= #gamma_{PCM}, #gamma_{conv,cluster}");
                legendFracDifferentContrib->Draw();

                PutProcessLabelAndEnergyOnPlot(0.15, 0.92, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
                canvasFracDifferentContrib->Update();
                canvasFracDifferentContrib->SaveAs(Form("%s/%s_RelativeContributionsToTruePeak_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
                delete canvasFracDifferentContrib;
            }
        }        
    }    

    // *******************************************************************************************
    // ****** Show fractions of cluster origin in MC to total for Calo related analysis path *****
    // *******************************************************************************************    
    if ((mode==4 || mode == 5) && kIsMC ){
        TH1D* histoTrueTotalRecYield                = (TH1D*)fileCorrections->Get("histoYieldTrueMesonFixedWindow");
        if (histoTrueTotalRecYield){
            TH1D* histoTrueTotalRecYieldGamma           = (TH1D*)fileCorrections->Get("histoYieldTrueMesonGammaFixedWindow");
            TH1D* histoTrueTotalRecYieldConvGamma       = (TH1D*)fileCorrections->Get("histoYieldTrueMesonGammaConvGammaFixedWindow");
            TH1D* histoTrueTotalRecYieldConvGamma2      = (TH1D*)fileCorrections->Get("histoYieldTrueMesonConvGammaConvGammaFixedWindow");
            TH1D* ratioTrueGammaDivTotal                = NULL;
            TH1D* ratioTrueConvGammaDivTotal            = NULL;
            TH1D* ratioTrueConvGamma2DivTotal           = NULL;
            if (histoTrueTotalRecYieldGamma){
                ratioTrueGammaDivTotal                  = (TH1D*)histoTrueTotalRecYieldGamma->Clone("ratioTrueGammaDivTotal");
                ratioTrueGammaDivTotal->Divide(ratioTrueGammaDivTotal,histoTrueTotalRecYield,1.,1.,"B");
            }
            if (histoTrueTotalRecYieldConvGamma){
                ratioTrueConvGammaDivTotal              = (TH1D*)histoTrueTotalRecYieldConvGamma->Clone("ratioTrueConvGammaDivTotal");
                ratioTrueConvGammaDivTotal->Divide(ratioTrueConvGammaDivTotal,histoTrueTotalRecYield,1.,1.,"B");
            }
            if (histoTrueTotalRecYieldConvGamma2){
                ratioTrueConvGamma2DivTotal             = (TH1D*)histoTrueTotalRecYieldConvGamma2->Clone("ratioTrueConvGamma2DivTotal");
                ratioTrueConvGamma2DivTotal->Divide(ratioTrueConvGamma2DivTotal,histoTrueTotalRecYield,1.,1.,"B");
            }
            
            if (histoTrueTotalRecYieldGamma){
                TCanvas* canvasFracDifferentContrib     = new TCanvas("canvasFracDifferentContrib","",200,10,1350,900);  // gives the page size
                DrawGammaCanvasSettings( canvasFracDifferentContrib, 0.09, 0.02, 0.04, 0.09);
                            
                DrawAutoGammaMesonHistos( ratioTrueGammaDivTotal, 
                                        "", "#it{p}_{T} (GeV/#it{c})", Form("#frac{%s->XX}{%s->ALL}", textProcess.Data(), textProcess.Data() ), 
                                        kFALSE, 1.5, 0, kFALSE,
                                        kTRUE, 0., 1.2, 
                                        kFALSE, 0., 10.);
                ratioTrueGammaDivTotal->GetYaxis()->SetTitleOffset(0.9);
                DrawGammaSetMarker(ratioTrueGammaDivTotal, 20, 1., kBlack, kBlack);  
                DrawGammaSetMarker(ratioTrueConvGammaDivTotal, 20, 1., kBlue, kBlue);
                DrawGammaSetMarker(ratioTrueConvGamma2DivTotal, 20, 1., kRed+2, kRed+2);
                ratioTrueGammaDivTotal->DrawCopy("e1");  
                ratioTrueConvGammaDivTotal->DrawCopy("e1,same");  
                ratioTrueConvGamma2DivTotal->DrawCopy("e1,same");  
                
                TLegend* legendFracDifferentContrib     = GetAndSetLegend2(0.6,0.93-3*0.035*1.1,0.93,0.93, 0.035, 1, "", 42, 0.1);
                legendFracDifferentContrib->AddEntry(ratioTrueGammaDivTotal,"XX= #gamma_{cluster}, #gamma_{cluster}");
                legendFracDifferentContrib->AddEntry(ratioTrueConvGammaDivTotal,"XX= #gamma_{cluster}, #gamma_{conv,cluster}");
                legendFracDifferentContrib->AddEntry(ratioTrueConvGamma2DivTotal,"XX= #gamma_{conv,cluster}, #gamma_{conv,cluster}");
                legendFracDifferentContrib->Draw();

                PutProcessLabelAndEnergyOnPlot(0.15, 0.92, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
                canvasFracDifferentContrib->Update();
                canvasFracDifferentContrib->SaveAs(Form("%s/%s_RelativeContributionsToTruePeak_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
                delete canvasFracDifferentContrib;
            }
        }        
    }    

    //***********************************************************************************************
    //***************************  correction for yield in mt bins **********************************
    //***********************************************************************************************    
    Int_t nBinsPt                       = histoCorrectedYieldTrue[0]->GetNbinsX();
    
    //*************************************************************************************************
    //******************** Output of the systematic Error due to Signal extraction ********************
    //*************************************************************************************************
    Double_t  binsXCenter[100];
    Double_t  binsXWidth[100];
    binsXWidth[0]=       0.;
    for (Int_t i = 1; i < nBinsPt +1; i++){
        binsXCenter[i]  = histoCorrectedYieldTrue[0]->GetBinCenter(i);
        binsXWidth[i]   = histoCorrectedYieldTrue[0]->GetBinWidth(i)/2.;
    }


    SysErrorConversion sysErr[6][100];
    for (Int_t k = 0; k < 6; k++){
        for (Int_t i = 1; i < nBinsPt +1; i++){
            if (mode == 9 || mode == 0 || mode == 1){
    //          binYValue[i] = histoCorrectedYieldTrue->GetBinContent(i);
                sysErr[k][i].value  = histoCorrectedYieldTrue[k]->GetBinContent(i);
                sysErr[k][i].error  = histoCorrectedYieldTrue[k]->GetBinError(i);
            } else {
                sysErr[k][i].value  = histoCorrectedYieldNorm[k]->GetBinContent(i);
                sysErr[k][i].error  = histoCorrectedYieldNorm[k]->GetBinError(i);
            }           
        }
    }    
    Double_t differencesToStandard[6][100];
    Double_t differencesToStandardErr[6][100];
    Double_t relDifferencesToStandard[6][100];
    Double_t relDifferencesToStandardErr[6][100];

    Double_t largestDifferenceNeg[100];
    Double_t largestDifferencePos[100];
    Double_t largestDifferenceNegError[100];
    Double_t largestDifferencePosError[100];
    Double_t relLargestDifferenceNeg[100];
    Double_t relLargestDifferencePos[100];
    Double_t relLargestDifferenceNegError[100];
    Double_t relLargestDifferencePosError[100];

    for (Int_t i = 1; i < nBinsPt +1; i++){
        largestDifferenceNeg[i]         = 0;
        largestDifferencePos[i]         = 0;
        largestDifferenceNegError[i]    = 0;
        largestDifferencePosError[i]    = 0;
        relLargestDifferenceNeg[i]      = 0;
        relLargestDifferencePos[i]      = 0;
        relLargestDifferenceNegError[i] = 0;
        relLargestDifferencePosError[i] = 0;
        //Calculate differences
        for (Int_t k = 0; k < 6; k++){
            differencesToStandard[k][i]     = sysErr[k][i].value - sysErr[0][i].value;
            differencesToStandardErr[k][i]  = TMath::Sqrt(TMath::Abs(TMath::Power(sysErr[k][i].error,2)-TMath::Power(sysErr[0][i].error,2)));
            if (sysErr[0][i].value != 0){
                relDifferencesToStandard[k][i]      = (differencesToStandard[k][i]/sysErr[0][i].value)*100.;
                relDifferencesToStandardErr[k][i]   = (differencesToStandardErr[k][i]/sysErr[0][i].value)*100.;;
                
            } else {
                relDifferencesToStandard[k][i]      = 0;
                relDifferencesToStandardErr[k][i]   = 0;
            }
            if (TMath::Abs(relDifferencesToStandard[k][i]) < 75. ){
                if( differencesToStandard[k][i] < 0 && differencesToStandard[k][i] < largestDifferenceNeg[i]){
                    largestDifferenceNeg[i]         = differencesToStandard[k][i];
                    largestDifferenceNegError[i]    = differencesToStandardErr[k][i];
                    relLargestDifferenceNeg[i]      = relDifferencesToStandard[k][i];
                    relLargestDifferenceNegError[i] = relDifferencesToStandardErr[k][i];
                } else if (differencesToStandard[k][i] > 0 && differencesToStandard[k][i] > largestDifferencePos[i]){
                    largestDifferencePos[i]         = differencesToStandard[k][i];
                    largestDifferencePosError[i]    = differencesToStandardErr[k][i];
                    relLargestDifferencePos[i]      = relDifferencesToStandard[k][i];
                    relLargestDifferencePosError[i] = relDifferencesToStandardErr[k][i];                    
                }    
            } 
        }
    }  

    cout << "done filling" << endl;
    const char *nameFileSysErrDat = Form("%s/%s/%s_%s_SystematicErrorYieldExtraction_%s.dat",fCutSelection.Data(),optionEnergy.Data() ,nameMeson.Data(),prefix2.Data(),fCutSelection.Data());
    fstream fileSysErrDat;
    fileSysErrDat.open(nameFileSysErrDat, ios::out);
    fileSysErrDat << "Calculation of the systematic error due to the yield extraction" << endl;
    fileSysErrDat << "Bin" << "\t" << "Normal value" << "\t" << "Normal error" << endl;
    for(Int_t i = 1; i < (nBinsPt +1); i++){
        fileSysErrDat << i << "\t" << sysErr[0][i].value << "\t" << sysErr[0][i].error << endl;
    }
    fileSysErrDat << endl;
    fileSysErrDat << "Bin" << "\t" << "Wide value" << "\t" << "Wide error" << "\t" << "Diff to Norm" << endl;
    for(Int_t i = 1; i < (nBinsPt +1); i++){
        fileSysErrDat << i  << "\t" << sysErr[1][i].value << "\t" << sysErr[1][i].error << "\t" << differencesToStandard[1][i] 
                            << "\t" << differencesToStandardErr[1][i] << "\t" << relDifferencesToStandard[1][i] << "\t" << relDifferencesToStandardErr[1][i] <<endl;
    }
    fileSysErrDat << endl;
    fileSysErrDat << "Bin" << "\t" << "Narrow value" << "\t" << "Narrow error" << "\t" << "Diff to Norm" << endl;
    for(Int_t i = 1; i < (nBinsPt +1); i++){
        fileSysErrDat << i  << "\t" << sysErr[2][i].value << "\t" << sysErr[2][i].error << "\t" << differencesToStandard[2][i] 
                            << "\t" << differencesToStandardErr[2][i] << "\t" << relDifferencesToStandard[2][i] << "\t" << relDifferencesToStandardErr[2][i] <<endl;
    }
    fileSysErrDat << endl;
    fileSysErrDat << "Bin" << "\t" << "Left value" << "\t" << "Left error" << "\t" << "Diff to Norm" << endl;
    for(Int_t i = 1; i < (nBinsPt +1); i++){
        fileSysErrDat << i  << "\t" << sysErr[3][i].value << "\t" << sysErr[3][i].error << "\t" << differencesToStandard[3][i] 
                            << "\t" << differencesToStandardErr[3][i] << "\t" << relDifferencesToStandard[3][i] << "\t" << relDifferencesToStandardErr[3][i] <<endl;
    }
    fileSysErrDat << endl;
    fileSysErrDat << "Bin" << "\t" << "Left Wide value" << "\t" << "Left Wide error" << "\t" << "Diff to Norm" << endl;
    for(Int_t i = 1; i < (nBinsPt +1); i++){
        fileSysErrDat << i  << "\t" << sysErr[4][i].value << "\t" << sysErr[4][i].error << "\t" << differencesToStandard[4][i] 
                            << "\t" << differencesToStandardErr[4][i] << "\t" << relDifferencesToStandard[4][i] << "\t" << relDifferencesToStandardErr[4][i] <<endl;
    }
    fileSysErrDat << endl;
    fileSysErrDat << "Bin" << "\t" << "Left Narrow value" << "\t" << "Left Narrow error" << "\t" << "Diff to Norm" << endl;
    for(Int_t i = 1; i < (nBinsPt +1); i++){
        fileSysErrDat << i  << "\t" << sysErr[5][i].value << "\t" << sysErr[5][i].error << "\t" << differencesToStandard[5][i] 
                            << "\t" << differencesToStandardErr[5][i] << "\t" << relDifferencesToStandard[5][i] << "\t" << relDifferencesToStandardErr[5][i] <<endl;
    }
    fileSysErrDat << endl;

    fileSysErrDat << "Bin" << "\t" << "Largest Dev Neg" << "\t" << "Largest Dev Pos"  << endl;
    for(Int_t i = 1; i < (nBinsPt +1); i++){
        fileSysErrDat << i << "\t" << largestDifferenceNeg[i] <<  "\t" << largestDifferenceNegError[i] << "\t" << largestDifferencePos[i] << "\t" << largestDifferencePosError[i]<<endl;
    }
    fileSysErrDat << "Bin" << "\t" << "Largest Rel Dev Neg" << "\t" << "Largest Rel Dev Pos"  << endl;
    for(Int_t i = 1; i < (nBinsPt +1); i++){
        fileSysErrDat << i << "\t" << relLargestDifferenceNeg[i] <<  "\t" << relLargestDifferenceNegError[i] << "\t" << relLargestDifferencePos[i] << "\t" << relLargestDifferencePosError[i]<<endl;
    }

    fileSysErrDat.close();
    
    TGraphAsymmErrors* SystErrGraphNeg = new TGraphAsymmErrors(nBinsPt+1, binsXCenter, relLargestDifferenceNeg, binsXWidth, binsXWidth, relLargestDifferenceNegError, relLargestDifferenceNegError);
    TGraphAsymmErrors* SystErrGraphPos = new TGraphAsymmErrors(nBinsPt+1, binsXCenter, relLargestDifferencePos, binsXWidth, binsXWidth, relLargestDifferencePosError, relLargestDifferencePosError);

    // ************************************** Plot sys error as calculated ********************************************************************
    TCanvas* canvasSysYieldExtraction = new TCanvas("canvasSysYieldExtraction","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasSysYieldExtraction, 0.065, 0.012, 0.035, 0.09);

    TH2F* histo2DDummySys;
    histo2DDummySys         = new TH2F("histo2DDummySys","histo2DDummySys",1000,0, histoTrueEffiPtUnmod[0]->GetXaxis()->GetBinUpEdge(histoTrueEffiPtUnmod[0]->GetNbinsX()),
                                                                            1000, -20, 20);
    SetStyleHistoTH2ForGraphs(histo2DDummySys, "#it{p}_{T} (GeV/#it{c})", "max Sys err %", 0.035,0.04, 0.035,0.04, 0.9,0.8, 510, 505);
    histo2DDummySys->DrawCopy();         

        DrawGammaSetMarkerTGraphAsym(SystErrGraphPos, 20, 1.,kBlue+1,kBlue+1);
        SystErrGraphPos->Draw("pX0,csame");
        DrawGammaSetMarkerTGraphAsym(SystErrGraphNeg, 21, 1.,kCyan+1,kCyan+1);
        SystErrGraphNeg->Draw("pX0,csame");
        
        DrawGammaSetMarker(fHistoYieldDiffBckResult[0], 31, 1.,kBlack,kBlack);
        fHistoYieldDiffBckResult[0]->DrawCopy("e1,p,SAME");
        DrawGammaSetMarker(fHistoYieldDiffBckResult[1], 33, 1.,kRed+2,kRed+2);
        fHistoYieldDiffBckResult[1]->DrawCopy("e1,p,SAME");
        DrawGammaSetMarker(fHistoYieldDiffBckResult[2], 34, 1.,kGreen+2,kGreen+2);
        fHistoYieldDiffBckResult[2]->DrawCopy("e1,p,SAME");

        TLegend* legendSys = GetAndSetLegend2(0.5,0.13,0.95,0.33, 0.035, 1, "", 42, 0.1);
        legendSys->AddEntry(SystErrGraphPos,"pos max","p"); 
        legendSys->AddEntry(SystErrGraphNeg,"neg max","p");
        legendSys->AddEntry(fHistoYieldDiffBckResult[0],"var pol2","p");
        legendSys->AddEntry(fHistoYieldDiffBckResult[1],"var exp","p");
        legendSys->AddEntry(fHistoYieldDiffBckResult[2],"var exp2","p");
        legendSys->Draw();

    canvasSysYieldExtraction->Update();
    PutProcessLabelAndEnergyOnPlot(0.12, 0.95, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());

    canvasSysYieldExtraction->SaveAs(Form("%s/%s_%s_SysYieldExtraction_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data())); 

    
    Double_t relBGEstimate[100];
    Double_t relBGEstimateError[100];
    for (Int_t i = 1; i < nBinsPt +1; i++){
        relBGEstimateError[i] = 0.;
        if ( TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateCatC->GetBinContent(i)) > TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateCatD->GetBinContent(i)) &&
             TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateA->GetBinContent(i)) > TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateA->GetBinContent(i))){
            relBGEstimate[i] = TMath::Abs(histoBGEstimateCatA->GetBinContent(i)- histoBGEstimateCatC->GetBinContent(i)) *100;
        } else if ( TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateCatD->GetBinContent(i)) > TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateCatC->GetBinContent(i)) &&
                    TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateA->GetBinContent(i)) > TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateA->GetBinContent(i))) {
            relBGEstimate[i] = TMath::Abs(histoBGEstimateCatA->GetBinContent(i)- histoBGEstimateCatC->GetBinContent(i)) *100;
        } else {
            relBGEstimate[i] = TMath::Abs(histoBGEstimateCatA->GetBinContent(i)- histoBGEstimateA->GetBinContent(i)) *100;
        }
    }   
    TGraphAsymmErrors* SystErrGraphBGEstimate = new TGraphAsymmErrors(nBinsPt+1, binsXCenter, relBGEstimate, binsXWidth, binsXWidth, relBGEstimateError, relBGEstimateError);

    
    // ********************************************************************************************************************************
    // ****************************** Write file with corrections only ****************************************************************
    // ********************************************************************************************************************************
    const char* nameOutput2 = Form("%s/%s/%s_%s_GammaConv_OnlyCorrectionFactor%s_%s.root",fCutSelection.Data(),optionEnergy.Data(),nameMeson.Data(),prefix2.Data(),optionPeriod.Data(),fCutSelection.Data());
    TFile* correctedOutput2 = new TFile(nameOutput2,"RECREATE");  
        if (histoAcceptance)                histoAcceptance->Write();
        if (histoTrueEffiPt[0])             histoTrueEffiPt[0]->Write("TrueMesonEffiPt");
        if (histoCompleteCorr)              histoCompleteCorr->Write("EffiTimesAcceptanceTimesDeltaY");
    correctedOutput2->Write();
    correctedOutput2->Close();
    
    // ********************************************************************************************************************************
    // ****************************** Write file with all further needed histograms ***************************************************
    // ********************************************************************************************************************************    
    const char* nameOutput = Form("%s/%s/%s_%s_GammaConvV1%sCorrection%s_%s.root",fCutSelection.Data(),optionEnergy.Data(),nameMeson.Data(),prefix2.Data(),fDalitz.Data(),optionPeriod.Data(),fCutSelection.Data());
    TFile* correctedOutput = new TFile(nameOutput,"RECREATE");  

        if (SystErrGraphPos)                    SystErrGraphPos->Write(Form("%s_SystErrorRelPos_YieldExtraction_%s",nameMeson.Data(),centralityString2.Data()),TObject::kOverwrite);
        if (SystErrGraphNeg)                    SystErrGraphNeg->Write(Form("%s_SystErrorRelNeg_YieldExtraction_%s",nameMeson.Data(),centralityString2.Data()),TObject::kOverwrite);
        if (SystErrGraphBGEstimate)             SystErrGraphBGEstimate->Write(Form("%s_SystErrorRel_BGEstimate_%s",nameMeson.Data(),centralityString2.Data()),TObject::kOverwrite);
        if (histoBGEstimateCatA)                histoBGEstimateCatA->Write("BGEstimateFromPileup");
        if (histoCorrectionFactorsHistvsPtCatA) histoCorrectionFactorsHistvsPtCatA->Write("PileupContamination");
        for (Int_t k = 0; k < 6; k++){
            if (histoCorrectedYieldTrue[k])     histoCorrectedYieldTrue[k]->Write();
            if (histoCorrectedYieldNorm[k])     histoCorrectedYieldNorm[k]->Write();
            
            // resonance feed down corrected
            if (!kIsEta) {
                if (histoFeedDownCorrectedYieldTrue[k])     histoFeedDownCorrectedYieldTrue[k]->Write();
                if (histoFeedDownCorrectedYieldNorm[k])     histoFeedDownCorrectedYieldNorm[k]->Write();
            }
        }
        if(!kIsEta){
            if (ratioFeedDownCorrectedYieldTrueToStandard)  ratioFeedDownCorrectedYieldTrueToStandard->Write();
        }
        for (Int_t k = 0; k < 3; k++){
            if (histoTrueEffiPt[k])             histoTrueEffiPt[k]->Write(Form("TrueMesonEffi%sPt",nameIntRange[k].Data()));
            if (histoEffiPt[k])                 histoEffiPt[k]->Write(Form("MesonEffi%sPt",nameIntRange[k].Data()));
            if (histoCorrectedYieldFixed[k])    histoCorrectedYieldFixed[k]->Write();
            if (histoCorrectedYieldTrueFixed[k])histoCorrectedYieldTrueFixed[k]->Write();
            
            // resonance feed down corrected
            if (!kIsEta) {
                if (histoFeedDownCorrectedYieldFixed[k])    histoFeedDownCorrectedYieldFixed[k]->Write();
                if (histoFeedDownCorrectedYieldTrueFixed[k])histoFeedDownCorrectedYieldTrueFixed[k]->Write();
            }
            
            for (Int_t j = 0; j < 4; j++){
                if (histoYieldSecMeson[k][j])           histoYieldSecMeson[k][j]->Write();
                if (histoRatioYieldSecMeson[k][j])      histoRatioYieldSecMeson[k][j]->Write();
                if (histoSecTrueEffi[k][j])             histoSecTrueEffi[k][j]->Write();
                if (k == 0){
                    if (histoSecAcceptance[j])              histoSecAcceptance[j]->Write();
                }    
                if (j < 3){
                    if (histoYieldSecMesonFromExternalInput[k][j])        histoYieldSecMesonFromExternalInput[k][j]->Write();
                    if (histoRatioYieldSecMesonFromExtInput[k][j])   histoRatioYieldSecMesonFromExtInput[k][j]->Write();
                }
            }
        }    
        
        if (histoUnCorrectedYield[0])           histoUnCorrectedYield[0]->Write();
        if (histoFWHMMeson)                     histoFWHMMeson->Write();
        if (histoMassMeson)                     histoMassMeson->Write();
        if (histoTrueFWHMMeson)                 histoTrueFWHMMeson->Write();
        if (histoTrueMassMeson)                 histoTrueMassMeson->Write();
        if (histoMassGaussianMeson)             histoMassGaussianMeson->Write("histoMassGaussianMeson");
        if (histoTrueMassGaussianMeson)         histoTrueMassGaussianMeson->Write("histoTrueMassGaussianMeson");
        if (histoMCrecMassMeson)                histoMCrecMassMeson->Write("histoMassMesonRecMC");
        if (histoMCrecMassGaussianMeson)        histoMCrecMassGaussianMeson->Write("histoMassGaussianMesonRecMC");
        if (histoWidthGaussianMeson)            histoWidthGaussianMeson->Write("histoWidthGaussianMeson");
        if (histoTrueWidthGaussianMeson)        histoTrueWidthGaussianMeson->Write("histoTrueWidthGaussianMeson");
        if (histoMCrecFWHMMeson)                histoMCrecFWHMMeson->Write("histoFWHMMesonRecMC");
        if (histoMCrecWidthGaussMeson)          histoMCrecWidthGaussMeson->Write("histoWidthGaussianMesonRecMC");
        if (histoRatioRecMass)                  histoRatioRecMass->Write("histoRatioRecMass");
        if (histoRatioValRecMass)               histoRatioValRecMass->Write("histoRatioValRecMass");
        if (histoRatioRecMassGauss)             histoRatioRecMassGauss->Write("histoRatioRecMassGauss");
        if (histoRatioValRecMassGauss)          histoRatioValRecMassGauss->Write("histoRatioValRecMassGauss");

        if (histoAcceptance)                    histoAcceptance->Write();
        if (histoAcceptanceWOEvtWeights)        histoAcceptanceWOEvtWeights->Write();
                                                                                                
        if (histoEventQuality)                  histoEventQuality->Write();
        if (histoNumberOfGoodESDTracksVtx)      histoNumberOfGoodESDTracksVtx->Write();
        if (histoInputMesonPt)                  histoInputMesonPt->Write();
        if (histoMCYieldMeson)                  histoMCYieldMeson->Write();
        if (histoMCYieldMesonOldBin)            histoMCYieldMesonOldBin->Write();
        if (histoInputMesonOldBinPtWOWeights)   histoInputMesonOldBinPtWOWeights->Write();
        if (histoInputMesonOldBinPtWeights)     histoInputMesonOldBinPtWeights->Write("WeightsMeson");
        if (histoMCInputAddedSig)               histoMCInputAddedSig->Write();
        if (histoMCInputWOWeightingAddedSig)    histoMCInputWOWeightingAddedSig->Write();
        if (histoMCInputWeightsAddedSig)        histoMCInputWeightsAddedSig->Write("WeightsMeson_AddedSig");
        if (histoUnCorrectedYieldDrawing){
            histoUnCorrectedYieldDrawing->SetName("histoYieldMesonPerEvent");
            histoUnCorrectedYieldDrawing->Write();
        }
        if (deltaPt)                            deltaPt->Write("deltaPt");
        if (histoInvMassSignalPlusBG)           histoInvMassSignalPlusBG->Write(Form("InvMassSigPlusBG_PtBin%02d",fExampleBin));
            else cout << "couldn't find SG+BG Inv Mass" << endl;
        if (histoInvMassBG)                     histoInvMassBG->Write(Form("InvMassBG_PtBin%02d",fExampleBin));
            else cout << "couldn't find BG Inv Mass" << endl;
        if (histoInvMassSignal)                 histoInvMassSignal->Write(Form("InvMassSig_PtBin%02d",fExampleBin));
            else cout << "couldn't find SG Inv Mass" << endl;
        if (fitInvMassSignal)                   fitInvMassSignal->Write(Form("FitInvMassSig_PtBin%02d",fExampleBin));
            else cout << "couldn't find Fit Inv Mass" << endl;
        
    correctedOutput->Write();
    correctedOutput->Close();
    
}
