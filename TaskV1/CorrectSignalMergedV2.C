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

struct SysErrorConversion {
    Double_t value;
    Double_t error;
};

void CorrectYield( TH1D* histoCorrectedYield, 
                   TH1D* histoPurity,
                   TH1D** histoRawSecYield, 
                   TH1D** histoRawSecYieldFromToy, 
                   TH1D* histoEffiPt, 
                   TH1D* histoAcceptance, 
                   Double_t deltaRapid, 
                   Double_t scaling, 
                   Double_t nEvt, 
                   TString nameMeson,
                   Bool_t isMC
                 ){
    histoCorrectedYield->Sumw2();
    histoCorrectedYield->Multiply(histoPurity);
    histoCorrectedYield->Scale(1./nEvt);
    if (histoRawSecYield){
        for (Int_t j = 0; j< 4;  j++){
            if (!isMC && j < 3 && histoRawSecYieldFromToy){
                if (histoRawSecYieldFromToy[j]){ 
                    histoCorrectedYield->Add(histoRawSecYieldFromToy[j],-1.);
                    cout << "will take secondary yield from toy approach for component: " << j << endl;  
                } else if (histoRawSecYield[j]) histoCorrectedYield->Add(histoRawSecYield[j],-1.);    
            } else {    
                if (histoRawSecYield[j]) histoCorrectedYield->Add(histoRawSecYield[j],-1.);
            }    
        }
    }        
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
    // no BR-ration correction needed
}

void CompileFullCorrectionFactor( TH1D* histoEffiPt, TH1D* histoAcceptance, Double_t deltaRapid){
    histoEffiPt->Sumw2();
    histoEffiPt->Multiply(histoEffiPt,histoAcceptance,1.,1.,"");
    histoEffiPt->Scale(deltaRapid);
}


void ScaleMCYield( TH1D* histoCorrectedToBeScaled, Double_t deltaRapid, Double_t scaling, Double_t nEvtMC, TString nameMeson){
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
    // no BR-ration correction needed
}


void  CorrectSignalMergedV2(    TString fileNameUnCorrectedFile = "myOutput", 
                                TString fileNameCorrectionFile  = "", 
                                TString fCutSelection           = "", 
                                TString suffix                  = "gif", 
                                TString nameMeson               = "", 
                                Bool_t  kIsMC                   = kFALSE, 
                                TString optionEnergy            = "", 
                                TString optionPeriod            = "",
                                Int_t mode                      = 10, 
                                TString fileNameTheory          = "ExternalInput/Theory/TheoryCompilationPP.root"
                            ){  
    
    //*************************************************************************************************
    //******************************** Set Style settings globally ************************************
    //*************************************************************************************************
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis(suffix);  
    SetPlotStyle();
        
    TString outputDir                           = "";
    if (!(optionPeriod.CompareTo("") == 0 || optionPeriod.CompareTo("No") == 0)){
        outputDir                               = Form("%s/%s/%s/%s/CorrectSignalMerged", fCutSelection.Data(), optionEnergy.Data(), optionPeriod.Data(), suffix.Data());
    } else {
        optionPeriod                            = "";
        outputDir                               = Form("%s/%s/%s/CorrectSignalMerged", fCutSelection.Data(), optionEnergy.Data(), suffix.Data());
    }    
    gSystem->Exec("mkdir -p "+outputDir);
        
    TString date                                = ReturnDateString();
    TString prefix2                             = "";
    TString textProcess                         = ReturnMesonString (nameMeson);
    if(textProcess.CompareTo("") == 0 ){
        cout << "Meson unknown" << endl;
        return ;
    }
    TString fTextMeasurement                    = Form("%s #rightarrow #gamma#gamma", textProcess.Data());
    TString fDetectionProcess                   = ReturnFullTextReconstructionProcess(mode, 0, textProcess);

    //*************************************************************************************************
    //******************************* Disentangle cut selections **************************************
    //*************************************************************************************************
    if (!(mode == 10 || mode == 11 )){
        cout << "incorrect mode: " << mode << endl;
        return ;
    }    
    TString fEventCutSelection                  = "";
    TString fClusterCutSelection                = "";
    TString fClusterMergedCutSelection          = "";
    TString dummyString                         = "";
    TString fMesonCutSelection                  = "";
    ReturnSeparatedCutNumberAdvanced( fCutSelection, fEventCutSelection, fClusterCutSelection, fClusterMergedCutSelection, dummyString, fMesonCutSelection, mode);
   
    //*************************************************************************************************
    //********************************** Set global options *******************************************
    //*************************************************************************************************
    TString collisionSystem                     = ReturnFullCollisionsSystem(optionEnergy);
    
    Double_t energy                             = ReturnCollisionEnergy( optionEnergy);
    Double_t doubleAddFactorK0s                 = ReturnCorrectK0ScalingFactor( optionEnergy,  fEventCutSelection);
    cout << "The additional K0 correction factor is: "  << doubleAddFactorK0s << endl;

    TString fNLMString                          = "";
    Int_t fNLMmin                               = ReturnClusterNLM(fClusterMergedCutSelection);
    if (ReturnClusterNLM(fClusterMergedCutSelection) == 1)
        fNLMString                  = Form("%i local maximum", ReturnClusterNLM(fClusterMergedCutSelection));        
    else 
        fNLMString                  = Form("%i local maxima", ReturnClusterNLM(fClusterMergedCutSelection));

    // Set flag for pi0/eta
    Bool_t kIsEta                               = kFALSE;
    if (!nameMeson.Contains("Pi0") ) 
        kIsEta                                  = kTRUE;
    TString textMeson                           = ReturnMesonString ( nameMeson);
    if (textMeson.CompareTo("") == 0) return;
    
    // Set flags for MC case
    if (kIsMC){ 
        prefix2                                 = "MC";
        doubleAddFactorK0s                      = 0.;
    } else { 
        prefix2                                 = "data";
    }
   
    if (collisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;  
    }
    
    Int_t kCollisionSystem                      = 0; // 0 : pp, 1: PbPb, 2: pPb
    if (optionEnergy.CompareTo("PbPb_2.76TeV") == 0) 
        kCollisionSystem                        = 1;
    if (optionEnergy.CompareTo("pPb_5.023TeV") == 0) 
        kCollisionSystem                        = 2;
    
    TString intermediate                        = GetCentralityString(fEventCutSelection);
    TString fTextCent                           = "";
    if (intermediate.CompareTo("pp")==0){
        fTextCent                               = "MinBias";  
        intermediate                            = "";
    } else {
        fTextCent                               = Form("%s central", intermediate.Data());
        if ( !intermediate.Contains("0-100%") )
            collisionSystem                     = Form("%s", collisionSystem.Data());
    }
    
    if (optionPeriod.CompareTo("") != 0 && optionPeriod.CompareTo("No") != 0 ){
        collisionSystem                         = Form("%s, %s",collisionSystem.Data(),optionPeriod.Data());
    }   
        
    TString rapidityRange                       = "";
    Double_t deltaRapid                         = ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
    TString trigger                             = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
    TString nameTrigger                         = ReturnTriggerName(trigger.Atoi());
    
    //Variable defintion
    Double_t scaling                            = 1./(2.*TMath::Pi());
    
    
    TString labelsBG[10]                        = { "PiCh", "P", "KCh", "N", "K0s", 
                                                    "Lambda", "Mu", "K0l", "Rest", "AllNonEM"};
    TString labelsBGPlot[10]                    = { "#pi^{#pm}", "p/#bar{p}", "K^{#pm}", "n", "K^{0}_{s}", 
                                                    "#Lambda", "#mu^{#pm}", "K^{0}_{L}", "rest", "all BG non em."};
    Color_t colorBGPlot[10]                     = { kRed+1, kBlue+1, kCyan-2, kTeal+4, 800,
                                                    kOrange+7, kGreen+2, kViolet+2, kAzure+5, kGray+2 };
    Style_t markerStyleBGPlot[10]               = { 20, 21, 24, 25, 33, 
                                                    27, 34, 24, 28, 20 };
    
    TString nameSecMeson[4]                     = {"K0S", "Lambda", "K0L", "Rest"};
    TString nameSecMesonPlot[4]                 = {"K_{s}^{0}", "#Lambda", "K_{l}^{0}", "Rest"};
    Color_t colorSec[4]                         = {kRed+2, 807, kCyan+2, kBlue};
    Style_t markerStyleSec[4]                   = {21, 33, 29, 34};
    Style_t markerStyleSecWithToy[4]            = {25, 27, 30, 28};
    Size_t markerSizeSec[4]                     = {1.5, 1.75, 2., 1.5};
    Color_t colorSecFromToy[3]                  = {kRed-2, kOrange+1, kCyan-2};
    Style_t markerStyleSecFromToy[3]            = {20, 33, 29};
    Size_t markerSizeSecFromToy[3]              = {1.7, 2, 2.2};
    Int_t pdgSecMeson[3]                        = {310, 3122, 130};
    Double_t decayLength[3]                     = {0.026844, 0.0789, 15.34};
    Double_t stacklength                        = 4.5;
    
    
    //************************************************************************************************
    //*********************** Loading uncorrected file ***********************************************
    //************************************************************************************************
    TFile fileUncorrected(fileNameUnCorrectedFile.Data());  
    if (fileUncorrected.IsZombie()) return;
    TH1D *histoEventQuality                     = (TH1D*)fileUncorrected.Get("NEvents");
    TH1D *histoUnCorrectedYield                 = (TH1D*)fileUncorrected.Get("histoYieldMesonM02");
    TH1D *deltaPt                               = (TH1D*)fileUncorrected.Get("deltaPt");
    for (Int_t i = 0; i < deltaPt->GetNbinsX() +1; i++){
        deltaPt->SetBinError(i, 0);
    }  
    
    Float_t nEvt                                = 0;
    if (kCollisionSystem == 1){
        nEvt                                    = histoEventQuality->GetBinContent(1);
    } else {
        nEvt                                    = GetNEvents(histoEventQuality);
        // BinContent 5 - Zvertex-position, BinContent 4 - no Trigger Bit, BinContent 7 - PileUp 
    }
    
    Double_t scaleFactorMeasXSecForToy              = 1;
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
            scaleFactorMeasXSecForToy               = xSectionEff/xSectionINEL;
        }    
    } 
    // The trigger rejection factor has to be given externally 
    Double_t triggerRejection                       = ReturnTriggerRejectionFactor(optionEnergy, trigger.Atoi());
    cout << "trigger rejection factor set to: " << triggerRejection << endl;
    if (triggerRejection != 1.)
        scaleFactorMeasXSecForToy                   = scaleFactorMeasXSecForToy*triggerRejection;
    if (scaleFactorMeasXSecForToy != 1)
        cout << "The secondary correction from the toy has to be scaled with " << scaleFactorMeasXSecForToy << endl;
//     return;

    // read pt distribution for different NLM bins (0-10, single NLM bins, 11 (>3), 12 (all))
    TH1D* histoPtInNLMBins[13]                      = { NULL, NULL, NULL, NULL, NULL,
                                                        NULL, NULL, NULL, NULL, NULL,
                                                        NULL, NULL, NULL};                                        
    TH1D* histoRatioPtInNLMBins[12]                 = { NULL, NULL, NULL, NULL, NULL,
                                                        NULL, NULL, NULL, NULL, NULL,
                                                        NULL, NULL};                                        
    Int_t nNLMBinsFilled                            = 0;
    Bool_t enableNLMBinsPlotting                    = kFALSE;
    for (Int_t i = 0; i< 11; i++){
        histoPtInNLMBins[i]                         = (TH1D*)fileUncorrected.Get(Form("ClusterMergedNLM%d",i+1));
        if (i > 2){
            if (!histoPtInNLMBins[11] && histoPtInNLMBins[i]){
                histoPtInNLMBins[11]                    = (TH1D*)histoPtInNLMBins[i]->Clone("ClusterMergedNLMBigger3");
                histoPtInNLMBins[11]->Sumw2();
            } else if (histoPtInNLMBins[i]){
                histoPtInNLMBins[11]->Add(histoPtInNLMBins[i]);            
            }                
        }    
        if (!histoPtInNLMBins[12] && histoPtInNLMBins[i]){
            cout << "generated summed histo" << endl; 
            nNLMBinsFilled++;
            histoPtInNLMBins[12]                    = (TH1D*)histoPtInNLMBins[i]->Clone("ClusterMergedNLMALL");
            histoPtInNLMBins[12]->Sumw2();
        } else if (histoPtInNLMBins[i]){
            nNLMBinsFilled++;
            histoPtInNLMBins[12]->Add(histoPtInNLMBins[i]);            
        }    
    }
    cout << "filled " << nNLMBinsFilled << endl;
    
    if (nNLMBinsFilled > 1){
        cout << "I am correct" << endl;
        enableNLMBinsPlotting                       = kTRUE;
        for (Int_t i = 0; i< 12; i++ ){
            cout << i << " found " << histoPtInNLMBins[i] << "\t" << histoPtInNLMBins[12] << endl;
            if (histoPtInNLMBins[i] && histoPtInNLMBins[12]){
                cout << "producing ratio of NLM bins for " << i+1 << endl;
                histoRatioPtInNLMBins[i]            = (TH1D*)histoPtInNLMBins[i]->Clone(Form("RatioClusterMergedNLM%d",i+1));
                histoRatioPtInNLMBins[i]->Divide(histoRatioPtInNLMBins[i],histoPtInNLMBins[12],1.,1.,"B");
            }
        }    
    }
    
    
    // read toy MC input if available
    TH1D* histoToyMCInputSecPi0[3]                  = {NULL, NULL, NULL};
    Bool_t foundToyMCInput                          = kFALSE;
    for (Int_t j = 0; j < 3; j++){
        histoToyMCInputSecPi0[j]                    = (TH1D*)fileUncorrected.Get(Form("histoSecPi0YieldFrom%s_FromToy",nameSecMeson[j].Data()));
        if (histoToyMCInputSecPi0[j]){
            histoToyMCInputSecPi0[j]->Scale(scaleFactorMeasXSecForToy);
            foundToyMCInput                         = kTRUE;
        }
    }

    // load theory prediction for direct to fragmentation photons
    TFile* fileTheory                           = new TFile(fileNameTheory.Data());
    TString nameTheoryGraph                     = Form("graphPromptPhotonDivFragementationNLOVogelsang_%s",((TString)ReturnCollisionEnergyStringForTheory(optionEnergy)).Data());
    TString nameTheoryFit                       = Form("ratioFitNLOPromptDivFragGamma%s",((TString)ReturnCollisionEnergyStringForTheory(optionEnergy)).Data());
    TDirectory* directoryTheoDirGamma           = (TDirectory*)fileTheory->Get("DirectPhoton");
    TGraphAsymmErrors* graphPromptDivFragTheo   = (TGraphAsymmErrors*)directoryTheoDirGamma->Get(nameTheoryGraph.Data());
    TF1* fitPromptdivFragTheo                   = (TF1*)directoryTheoDirGamma->Get(nameTheoryFit.Data());
    
    // set min and max values for pT
    Double_t maxPtMeson     = histoUnCorrectedYield->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield->GetNbinsX());
    Double_t minPtMeson     = 0;
    Int_t ptBin             = 1;
    while (histoUnCorrectedYield->GetBinContent(ptBin) == 0. && ptBin < histoUnCorrectedYield->GetNbinsX()){
        ptBin++;
        minPtMeson          = histoUnCorrectedYield->GetXaxis()->GetBinLowEdge(ptBin);
    }
    
    Double_t minPtMesonSec  = minPtMeson;
    cout << "minimum pT: " << minPtMeson << ", maximum pT: " << maxPtMeson << endl;
     
    //************************************************************************************************
    //************************** Loading corrections file ********************************************
    //************************************************************************************************
    TFile* fileCorrections                      = new TFile(fileNameCorrectionFile.Data());
    if (fileCorrections->IsZombie()) return;
    TH1F *histoEventQualityMC                   = (TH1F*)fileCorrections->Get("NEvents");
    TString namePurity                          = "TruePurityPi0Pt";
    if (kIsEta) namePurity                      = "TruePurityEtaPt";
    TH1D* histoMesonPurityPt                    = (TH1D*)fileCorrections->Get(namePurity.Data());

    TH1D *histoAcceptance                       = (TH1D*)fileCorrections->Get("AcceptancePt");
    
    // read pt distribution for different NLM bins (0-10, single NLM bins, 11 (>3), 12 (all))
    TH1D* histoMCPtInNLMBins[13]                = { NULL, NULL, NULL, NULL, NULL,
                                                    NULL, NULL, NULL, NULL, NULL,
                                                    NULL, NULL, NULL};                                        
    TH1D* histoMCRatioPtInNLMBins[12]           = { NULL, NULL, NULL, NULL, NULL,
                                                    NULL, NULL, NULL, NULL, NULL,
                                                    NULL, NULL};                                        
    if (enableNLMBinsPlotting){
        for (Int_t i = 0; i< 11; i++){
            histoMCPtInNLMBins[i]                           = (TH1D*)fileCorrections->Get(Form("ClusterMergedNLM%d",i+1));
            if (histoMCPtInNLMBins[i]){
                histoMCPtInNLMBins[i]->SetName(Form("ClusterMergedNLM%dMC",i+1));
                if (i > 2){
                    if (!histoMCPtInNLMBins[11]){
                        histoMCPtInNLMBins[11]              = (TH1D*)histoMCPtInNLMBins[i]->Clone("ClusterMergedNLMBigger3MC");
                        histoMCPtInNLMBins[11]->Sumw2();
                    } else {
                        histoMCPtInNLMBins[11]->Add(histoMCPtInNLMBins[i]);            
                    }                
                }    
                if (!histoMCPtInNLMBins[12] ){
                    histoMCPtInNLMBins[12]                  = (TH1D*)histoMCPtInNLMBins[i]->Clone("ClusterMergedNLMALLMC");
                    histoMCPtInNLMBins[12]->Sumw2();
                } else {
                    histoMCPtInNLMBins[12]->Add(histoMCPtInNLMBins[i]);            
                }    
            }
        }
        for (Int_t i = 0; i< 12; i++ ){
            if (histoMCPtInNLMBins[i] && histoMCPtInNLMBins[12]){
                histoMCRatioPtInNLMBins[i]            = (TH1D*)histoMCPtInNLMBins[i]->Clone(Form("RatioClusterMergedNLM%dMC",i+1));
                histoMCRatioPtInNLMBins[i]->Divide(histoMCRatioPtInNLMBins[i],histoMCPtInNLMBins[12],1.,1.,"B");
            }
        }
    }
    
    TH1D* histoTrueClustersBGPt[10]             = { NULL, NULL, NULL, NULL, NULL,
                                                    NULL, NULL, NULL, NULL, NULL};                                        
    TH1D* histoRatioTrueClustersBGPt[10]        = { NULL, NULL, NULL, NULL, NULL,
                                                    NULL, NULL, NULL, NULL, NULL};                                        
    TH1D* histoTrueYieldPi0M02                  = NULL;
    TH1D* histoTrueYieldEtaM02                  = NULL;
    
    TH1D* histoTrueYieldPi0DCM02                = NULL;
    TH1D* histoTrueYieldPi0GGM02                = NULL;
    TH1D* histoTrueYieldPi0DalitzM02            = NULL;
    TH1D* histoTrueYieldEtaDCM02                = NULL;
    TH1D* histoTrueYieldEtaGGM02                = NULL;
    TH1D* histoTrueYieldEtaDalitzM02            = NULL;
    TH1D* histoTrueYieldGammaM02                = NULL;
    TH1D* histoTrueYieldDirGammaM02             = NULL;
    TH1D* histoTrueYieldElectronM02             = NULL;
    TH1D* histoRatioTrueYieldEtaM02             = NULL;
    TH1D* histoRatioTrueYieldPi0M02             = NULL;
    TH1D* histoRatioTrueYieldGammaM02           = NULL;
    TH1D* histoRatioTrueYieldElectronM02        = NULL;
    TH1D* histoRatioPi0DCFrac                   = NULL;
    TH1D* histoRatioPi0GGFrac                   = NULL;
    TH1D* histoRatioPi0DalitzFrac               = NULL;
    TH1D* histoRatioEtaDCFrac                   = NULL;
    TH1D* histoRatioEtaGGFrac                   = NULL;
    TH1D* histoRatioEtaDalitzFrac               = NULL;
    TH1D* histoMergedAll                        = NULL;
    TH1D* histoMergedPure                       = NULL;
    TH1D* histoMergedPartConv                   = NULL;
    TH1D* histoMergedOneGamma                   = NULL;
    TH1D* histoMergedOneElectron                = NULL;
    TH1D* histoRatioMergedPureFrac              = NULL;
    TH1D* histoRatioMergedPartConvFrac          = NULL;
    TH1D* histoRatioMergedOneGammaFrac          = NULL;
    TH1D* histoRatioMergedOneElectronFrac       = NULL;
    Bool_t  isUpdatedOutputFormat               = kFALSE;
    if (kIsMC){
        for (Int_t i = 0; i < 10; i++){
            histoTrueClustersBGPt[i]                = (TH1D*)fileCorrections->Get(Form("TrueClusBG_%s_Pt",labelsBG[i].Data()));
            histoTrueClustersBGPt[i]->Sumw2();            
            histoRatioTrueClustersBGPt[i]           = (TH1D*)fileCorrections->Get(Form("RatioTrueClusBG_%s_Pt",labelsBG[i].Data()));
        }
        
        histoTrueYieldEtaM02                        = (TH1D*)fileCorrections->Get("histoTrueYieldEtaM02");
        histoRatioTrueYieldEtaM02                   = (TH1D*)fileCorrections->Get("RatioTrueYieldEtaM02");
        histoTrueYieldPi0M02                        = (TH1D*)fileCorrections->Get("histoTrueYieldPi0M02");
        histoRatioTrueYieldPi0M02                   = (TH1D*)fileCorrections->Get("RatioTrueYieldPi0M02");
        
        histoTrueYieldGammaM02                      = (TH1D*)fileCorrections->Get("histoTrueYieldGammaM02");
        histoRatioTrueYieldGammaM02                 = (TH1D*)fileCorrections->Get("RatioTrueYieldGammaM02");
                
        histoTrueYieldElectronM02                   = (TH1D*)fileCorrections->Get("histoTrueYieldElectronM02");
        histoRatioTrueYieldElectronM02              = (TH1D*)fileCorrections->Get("RatioTrueYieldElectronM02");
        
        histoTrueYieldPi0DCM02                      = (TH1D*)fileCorrections->Get("histoTrueYieldPi0DCM02");
        histoRatioPi0DCFrac                         = (TH1D*)fileCorrections->Get("RatioPi0DCFrac");

        histoTrueYieldPi0GGM02                      = (TH1D*)fileCorrections->Get("histoTrueYieldPi0GGM02");
        histoRatioPi0GGFrac                         = (TH1D*)fileCorrections->Get("RatioPi0GGFrac");

        histoTrueYieldPi0DalitzM02                  = (TH1D*)fileCorrections->Get("histoTrueYieldPi0DalitzM02");
        histoRatioPi0DalitzFrac                     = (TH1D*)fileCorrections->Get("RatioPi0DalitzFrac");
        
        histoTrueYieldEtaDCM02                      = (TH1D*)fileCorrections->Get("histoTrueYieldEtaDCM02");
        histoRatioEtaDCFrac                         = (TH1D*)fileCorrections->Get("RatioEtaDCFrac");
        
        histoTrueYieldEtaGGM02                      = (TH1D*)fileCorrections->Get("histoTrueYieldEtaGGM02");
        histoRatioEtaGGFrac                         = (TH1D*)fileCorrections->Get("RatioEtaGGFrac");

        histoTrueYieldEtaDalitzM02                  = (TH1D*)fileCorrections->Get("histoTrueYieldEtaDalitzM02");
        histoRatioEtaDalitzFrac                     = (TH1D*)fileCorrections->Get("RatioEtaDalitzFrac");

        histoMergedAll                              = (TH1D*)fileCorrections->Get("histoTrueYieldMergedM02");
        histoMergedPure                             = (TH1D*)fileCorrections->Get(Form("histoTrueYieldMergedPureFrom%sM02",nameMeson.Data()));
        cout << "maximum of pure merged: "<< histoMergedPure->GetMaximum() << endl;
        if (!(histoMergedPure->GetMaximum() > 0)){
            histoMergedPure                         = (TH1D*)fileCorrections->Get("histoTrueYieldMergedPureM02");
            histoRatioMergedPureFrac                = (TH1D*)fileCorrections->Get(Form("Ratio%sMergedPure",nameMeson.Data()));
        } else { 
            isUpdatedOutputFormat                   = kTRUE;
            histoRatioMergedPureFrac                = (TH1D*)histoMergedPure->Clone("RatioMergedPure");
            histoRatioMergedPureFrac->Divide(histoRatioMergedPureFrac,histoTrueYieldPi0M02,1.,1.,"B");
            histoRatioMergedPureFrac->Scale(100.);
        }
        if (!isUpdatedOutputFormat){
            histoMergedPartConv                     = (TH1D*)fileCorrections->Get("histoTrueYieldMergedPartConvM02");
            histoRatioMergedPartConvFrac            = (TH1D*)histoMergedPartConv->Clone("RatioMergedPartConv");
            histoRatioMergedPartConvFrac->Divide(histoRatioMergedPartConvFrac,histoMergedAll,1.,1.,"B");
            histoRatioMergedPartConvFrac->Scale(100.);
            histoMergedOneGamma                     = (TH1D*)fileCorrections->Get("histoTrueYieldMergedOneGammaFromMesonM02");
            histoRatioMergedOneGammaFrac            = (TH1D*)histoMergedOneGamma->Clone("RatioMergedOneGamma");
            histoRatioMergedOneGammaFrac->Divide(histoRatioMergedOneGammaFrac,histoMergedAll,1.,1.,"B");
            histoRatioMergedOneGammaFrac->Scale(100.);
            histoMergedOneElectron                  = (TH1D*)fileCorrections->Get("histoTrueYieldMergedOneElectronFromMesonM02");
            histoRatioMergedOneElectronFrac         = (TH1D*)histoMergedOneElectron->Clone("RatioMergedOneElectron");
            histoRatioMergedOneElectronFrac->Divide(histoRatioMergedOneElectronFrac,histoMergedAll,1.,1.,"B");
            histoRatioMergedOneElectronFrac->Scale(100.);
        } else {
            cout << "went into identified meson routine" << endl;
            histoMergedPartConv                     = (TH1D*)fileCorrections->Get(Form("histoTrueYieldMergedPartConvFrom%sM02",nameMeson.Data()));
            histoRatioMergedPartConvFrac            = (TH1D*)fileCorrections->Get(Form("Ratio%sMergedPartConv",nameMeson.Data()));
            histoMergedOneGamma                     = (TH1D*)fileCorrections->Get(Form("histoTrueYieldMergedOneGammaFrom%sM02",nameMeson.Data()));
            histoRatioMergedOneGammaFrac            = (TH1D*)fileCorrections->Get(Form("Ratio%sMergedOneGamma",nameMeson.Data()));
            histoMergedOneElectron                  = (TH1D*)fileCorrections->Get(Form("histoTrueYieldMergedOneElectronFrom%sM02",nameMeson.Data()));
            histoRatioMergedOneElectronFrac         = (TH1D*)fileCorrections->Get(Form("Ratio%sMergedOneElectron",nameMeson.Data()));
        }    
    }
    
    // scale additional photon component according to theory calculations, as only FF photons are in our MC's
    TH1D* histoRatioAdditionalGammaCorrM02          = NULL;
    TH1D* histoMesonPurityUnmodPt                   = NULL;
    if (!kIsMC){
        histoRatioAdditionalGammaCorrM02            = (TH1D*)fileCorrections->Get("RatioTrueYieldGammaM02");
        if (fitPromptdivFragTheo){
            cout << "found theo scaling fac" <<  endl;
            cout << "adjusting gamma contribution according theory predictions" <<  endl;
            histoRatioAdditionalGammaCorrM02->Multiply(fitPromptdivFragTheo);
            histoMesonPurityUnmodPt                 = (TH1D*)histoMesonPurityPt->Clone("histoMesonPurityUnmodPt");
            histoMesonPurityPt->Add(histoRatioAdditionalGammaCorrM02,-1);
        }
        if (optionEnergy.CompareTo("8TeV") == 0 && nameMeson.CompareTo("Pi0") == 0 ){
            cout << "adjusting eta contribution according data/MC comparison for eta/pi0" <<  endl;
            TH1D* histoRatioAdditionalEtaCorrM02    = (TH1D*)fileCorrections->Get("RatioTrueYieldEtaM02");
            for (Int_t i = histoRatioAdditionalEtaCorrM02->GetXaxis()->FindBin(11); i< histoRatioAdditionalEtaCorrM02->GetNbinsX(); i++){
                cout << histoRatioAdditionalEtaCorrM02->GetBinCenter(i) << "\t" << histoRatioAdditionalEtaCorrM02->GetBinContent(i) << endl;
            }    
            histoRatioAdditionalEtaCorrM02->Scale(0.43/0.407-1.); // this factor is first the measured eta/pi0 in data and the eta/pi0 from the Pythia8, JJ at high pt
            for (Int_t i = histoRatioAdditionalEtaCorrM02->GetXaxis()->FindBin(11); i< histoRatioAdditionalEtaCorrM02->GetNbinsX(); i++){
                cout << histoRatioAdditionalEtaCorrM02->GetBinCenter(i) << "\t" << histoRatioAdditionalEtaCorrM02->GetBinContent(i) << endl;
            }    
            histoMesonPurityPt->Add(histoRatioAdditionalEtaCorrM02,-1);
        }
    }
    
    // loading efficiency
    TH1D *histoTrueEffiPt                       = NULL;
    TH1D *histoTrueEffiPrimMesonPt              = NULL;
    histoTrueEffiPrimMesonPt                    = (TH1D*)fileCorrections->Get("TrueEfficiencyPrimMesonPt"); 
    
    // loading MC input spectra
    TH1D* histoInputMesonPt                     = (TH1D*)fileCorrections->Get("MC_Meson_Rebin");
    TH1D* histoInputMesonOldBinPt               = (TH1D*)fileCorrections->Get("MC_Meson");
        

    Float_t nEvtMC                              = 0;
    if (kCollisionSystem == 1){
        nEvtMC                                  = histoEventQualityMC->GetBinContent(1);
    } else {
        nEvtMC                                  = GetNEvents(histoEventQualityMC);
        // BinContent 5 - Zvertex-position, BinContent 4 - no Trigger Bit, BinContent 7 - PileUp 
    }
  
    TH1D *histoTrueSecFracMeson[4]              = { NULL, NULL, NULL, NULL };
    TH1D *histoTrueSecFracMeson_Or[4]           = { NULL, NULL, NULL, NULL };
    TF1* fitSecFrac[4]                          = { NULL, NULL, NULL, NULL };
    for (Int_t j = 0; j < 4; j++){
        fitSecFrac[j]                           = new TF1(Form("fitSecFracFrom%s",nameSecMeson[j].Data()),"[0]/pow(x,[1])+[2]+[3]*x");
        fitSecFrac[j]->SetParLimits(2,0,10);
        fitSecFrac[j]->SetParLimits(3,0,1e-2);
    }    
    TH1D* histoYieldSecMeson[4]                 = { NULL, NULL, NULL, NULL };
    TH1D* histoYieldSecMesonFromToy[4]          = { NULL, NULL, NULL, NULL };
    TH1D* histoRatioYieldSecMeson[4]            = { NULL, NULL, NULL, NULL };
    TH1D* histoRatioYieldSecMesonFromToy[4]     = { NULL, NULL, NULL, NULL };
    TH1D* histoSecEffiPt[4]                     = { NULL, NULL, NULL, NULL };
    Int_t nEffHistSec                               = 0;
    Bool_t modifiedSecEff[4]                    = {kFALSE, kFALSE, kFALSE, kFALSE};
    TH1D* histoRatioSecEffDivTrueEff[4]         = { NULL, NULL, NULL, NULL };
    TH1D* histoSecAccPt[4]                      = { NULL, NULL, NULL, NULL };
    Bool_t modifiedSecAcc[4]                    = {kFALSE, kFALSE, kFALSE, kFALSE};
    Int_t nAccHistSec                           = 0;
    Double_t scalingFacSec[4]                   = { 1+doubleAddFactorK0s, 1, 1, 1};
    Bool_t isNewOutput                          = kFALSE;
    Bool_t doK0SecCorrection                    = kFALSE;
    Bool_t haveSec[4]                           = { kTRUE, kTRUE, kTRUE, kTRUE };
    Bool_t haveSecUsed[4]                       = { kTRUE, kTRUE, kTRUE, kTRUE };

    Int_t doK0SecCorrectionWithDefaultHisto     = 0;
    
    if ( !kIsEta ) {
        doK0SecCorrection                       = kTRUE;
    }
    
    if (doK0SecCorrection){
        for (Int_t j = 0; j< 4; j++){
            cout << "Fit standard range - secondaries: " << nameSecMeson[j].Data() << endl;  
            histoTrueSecFracMeson[j]            = (TH1D*)fileCorrections->Get(Form("TrueSecFracFrom%s",nameSecMeson[j].Data()));
            if (histoTrueSecFracMeson[j]){
                histoTrueSecFracMeson_Or[j]         = (TH1D*)histoTrueSecFracMeson[j]->Clone(Form("TrueSecFracFrom%s_Or",nameSecMeson[j].Data()));
                Double_t maxPtSecondaries           = histoTrueSecFracMeson[j]->GetXaxis()->GetBinUpEdge(histoTrueSecFracMeson[j]->GetNbinsX());
                fitSecFrac[j]->SetRange(minPtMesonSec, maxPtSecondaries);
                TFitResultPtr resultSecFrac         = histoTrueSecFracMeson[j]->Fit(fitSecFrac[j], "SNRME+", "", minPtMesonSec, maxPtSecondaries);
                
                for (Int_t i = 1; i < histoTrueSecFracMeson[j]->GetNbinsX()+1; i++){
                    Double_t ptStart                    = histoTrueSecFracMeson[j]->GetXaxis()->GetBinLowEdge(i);
                    Double_t ptEnd                      = histoTrueSecFracMeson[j]->GetXaxis()->GetBinUpEdge(i);
                    Double_t binWidth                   = ptEnd-ptStart;
                    Double_t secFrac                    = fitSecFrac[j]->Integral(ptStart, ptEnd, resultSecFrac->GetParams()) / binWidth;
                    Double_t errorSecFrac               = fitSecFrac[j]->IntegralError(ptStart, ptEnd, resultSecFrac->GetParams(), resultSecFrac->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                    histoTrueSecFracMeson[j]->SetBinContent(i, secFrac);
                    histoTrueSecFracMeson[j]->SetBinError(i, errorSecFrac);
                }
            } else {
                haveSec[j]                          = kFALSE;
                haveSecUsed[j]                      = kFALSE;
                histoTrueSecFracMeson[j]            = (TH1D*)histoTrueSecFracMeson[0]->Clone(Form("TrueSecFracFrom%s",nameSecMeson[j].Data()));
                histoTrueSecFracMeson_Or[j]         = (TH1D*)histoTrueSecFracMeson_Or[0]->Clone(Form("TrueSecFracFrom%s_Or",nameSecMeson[j].Data()));
                for (Int_t i = 1; i < histoTrueSecFracMeson[j]->GetNbinsX()+1; i++){
                    histoTrueSecFracMeson[j]->SetBinContent(i, 0);
                    histoTrueSecFracMeson[j]->SetBinError(i, 0);
                    histoTrueSecFracMeson_Or[j]->SetBinContent(i, 0);
                    histoTrueSecFracMeson_Or[j]->SetBinError(i, 0);
                }
            }    
            
            histoYieldSecMeson[j]                   = (TH1D*)histoUnCorrectedYield->Clone(Form("SecYieldMesonFrom%s", nameSecMeson[j].Data()));
            histoYieldSecMeson[j]->Sumw2();
            histoYieldSecMeson[j]->Multiply(histoMesonPurityPt);
            histoYieldSecMeson[j]->Multiply(histoTrueSecFracMeson[j]);
            histoYieldSecMeson[j]->Scale(scalingFacSec[j]*1./nEvt);       
        }    
        
        for (Int_t j = 0; j< 4; j++){
            histoSecAccPt[j]                        = (TH1D*)fileCorrections->Get(Form("fMCSecPi0From%sAccepPt",nameSecMeson[j].Data()));
            histoSecEffiPt[j]                       = (TH1D*)fileCorrections->Get(Form("TrueMesonEffiSecFrom%sPt",nameSecMeson[j].Data()));
            // correcting sec acceptance due to statistical fluctuations
            if (histoSecAccPt[j]){
                nAccHistSec++;
                if ( j == 0 ){
                    histoSecAccPt[j]                    = (TH1D*)histoAcceptance->Clone(Form("AcceptanceSecPi0From%s",nameSecMeson[j].Data()));
                    modifiedSecAcc[j]                   = kTRUE;
                } else if ( j == 2 || j == 1 ){
                    histoSecAccPt[j]->SetName(Form("AcceptanceSecPi0From%s",nameSecMeson[j].Data()));
                    Double_t    accSec                  = ReturnDeltaEtaCalo(fClusterCutSelection)/deltaRapid*ReturnDeltaPhiCalo(fClusterCutSelection)/(2*TMath::Pi());
                    for (Int_t iPt = histoSecAccPt[j]->FindBin(minPtMeson); iPt< histoSecAccPt[j]->GetNbinsX()+1; iPt++ ){
                        histoSecAccPt[j]->SetBinContent(iPt, accSec);
                        histoSecAccPt[j]->SetBinError(iPt, accSec*0.1);
                    }
                    modifiedSecAcc[j]                   = kTRUE; 
                }
            }
            // adjusting sec eff if needed 
            if (histoSecEffiPt[j]){
                nEffHistSec++;
                histoRatioSecEffDivTrueEff[j]       = (TH1D*)histoSecEffiPt[j]->Clone(Form("ratioSecEffDivTrueEff%s",nameSecMeson[j].Data()));
                histoRatioSecEffDivTrueEff[j]->Divide(histoRatioSecEffDivTrueEff[j],histoTrueEffiPrimMesonPt);

                TF1*  fitConst                      = new TF1("fitConst","[0]");
                fitConst->SetLineColor(colorSec[j]);
                histoRatioSecEffDivTrueEff[j]->Fit(fitConst);
                cout << fitConst->GetParameter(0) << "\t +-" << fitConst->GetParError(0) << endl;
                if (j == 0){
                    for (Int_t iPt = histoSecEffiPt[j]->FindBin(minPtMeson); iPt< histoSecEffiPt[j]->GetNbinsX()+1; iPt++ ){
                        if (histoSecEffiPt[j]->GetBinContent(iPt) == 0){
                            Double_t ratioPrevBin           = 2;
                            if (iPt-1 > histoSecEffiPt[j]->FindBin(minPtMeson) && histoSecEffiPt[j]->GetBinContent(iPt-1) > 0){
                                ratioPrevBin                = histoSecEffiPt[j]->GetBinContent(iPt-1)/histoTrueEffiPrimMesonPt->GetBinContent(iPt-1);
                            }    
                            histoSecEffiPt[j]->SetBinContent(iPt, ratioPrevBin*histoTrueEffiPrimMesonPt->GetBinContent(iPt));
                            histoSecEffiPt[j]->SetBinError(iPt, ratioPrevBin*histoTrueEffiPrimMesonPt->GetBinError(iPt));
                            modifiedSecEff[j]               = kTRUE;
                        }
                    }    
                } else if (j == 1){
                    histoSecEffiPt[j]               = (TH1D*)histoSecEffiPt[0]->Clone(Form("TrueMesonEffiSecFrom%sPt",nameSecMeson[j].Data()));
                    modifiedSecEff[j]               = kTRUE;
                } else if ( j == 2 ){
                    for (Int_t iPt = histoSecEffiPt[j]->FindBin(minPtMeson); iPt< histoSecEffiPt[j]->GetNbinsX()+1; iPt++ ){
                        histoSecEffiPt[j]->SetBinContent(iPt, 1.2);
                        histoSecEffiPt[j]->SetBinError(iPt, 0.12);
                    }
                    modifiedSecEff[j]               = kTRUE;
                }    
                
                if (modifiedSecEff[j])
                    cout << "adjusted sec effi, due to to little stat: " << j << endl;
            }
        }

        if (histoSecAccPt[0] && histoSecEffiPt[0] )
            isNewOutput                             = kTRUE;
        
        if (isNewOutput){
            for (Int_t j = 0; j < 3; j++){
                if (histoSecAccPt[j] && histoSecEffiPt[j] && histoToyMCInputSecPi0[j]){
                    histoYieldSecMesonFromToy[j]    = (TH1D*)histoToyMCInputSecPi0[j]->Clone(Form("SecYieldFrom%sMesonFromToy", nameSecMeson[j].Data()));
                    histoYieldSecMesonFromToy[j]->Sumw2();
                    histoYieldSecMesonFromToy[j]->Multiply(histoSecAccPt[j]);
                    histoYieldSecMesonFromToy[j]->Multiply(histoSecEffiPt[j]);                                        
                }    
            }    
        }
    }
    
    Int_t nSecComp                  = 0;
    Int_t nSecCompUsed              = 0;
    for (Int_t j = 0; j < 4; j++){
        if (haveSec[j]) nSecComp++;
        if (haveSecUsed[j]) nSecCompUsed++;
    }    

    
    //**********************************************************************************
    //******************** Acceptance Plot *********************************************
    //**********************************************************************************
    
    if (kIsMC){
        cout << "Plotting acceptance" << endl;
        TCanvas* canvasAcceptance               = new TCanvas("canvasAcceptance2","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasAcceptance, 0.07, 0.01, 0.02, 0.08);

        DrawAutoGammaMesonHistos( histoAcceptance, 
                                    "", "#it{p}_{T} (GeV/#it{c})", Form("A_{%s} in |#it{y}| < %s",textMeson.Data(),rapidityRange.Data()), 
                                    kTRUE, 1.3, 3e-6, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
        histoAcceptance->GetYaxis()->SetTitleOffset(0.85);        
        DrawGammaSetMarker(histoAcceptance, 20, 1.5, kAzure-6, kAzure-6);
        histoAcceptance->DrawCopy("e1"); 

        PutProcessLabelAndEnergyOnPlot(0.7, 0.25, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);

        canvasAcceptance->Update();
        canvasAcceptance->SaveAs(Form("%s/%s_Acceptance_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
        
        if (isNewOutput){
            canvasAcceptance->cd();

                TLegend* legendSecAcc = GetAndSetLegend2(0.1, 0.12, 0.3, 0.12+(nAccHistSec+1)*0.03*1.05, 32,1); 
                histoAcceptance->GetYaxis()->SetRangeUser(0,histoAcceptance->GetMaximum()*1.8);
                DrawGammaSetMarker(histoAcceptance, 20, 1.5, kAzure-6, kAzure-6);
                histoAcceptance->DrawCopy("e1"); 
                legendSecAcc->AddEntry(histoAcceptance,"prim.");
                for (Int_t j = 0; j < 4; j++){
                    if (histoSecAccPt[j]){
                        DrawGammaSetMarker(histoSecAccPt[j], markerStyleSec[j], markerSizeSec[j], colorSec[j], colorSec[j]);
                        histoSecAccPt[j]->DrawCopy("e1,same"); 
                        if (!modifiedSecAcc[j]){
                            legendSecAcc->AddEntry(histoSecAccPt[j],Form("sec. from %s", nameSecMesonPlot[j].Data()));
                        } else {
                            legendSecAcc->AddEntry(histoSecAccPt[j],Form("sec. from %s, adj.", nameSecMesonPlot[j].Data()));
                        }    
                    }    
                }
                legendSecAcc->Draw();

                PutProcessLabelAndEnergyOnPlot(0.7, 0.25, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);

            canvasAcceptance->Update();
            canvasAcceptance->SaveAs(Form("%s/%s_AcceptanceWithSec_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
            
            
        }
        
        delete canvasAcceptance;
        
        //**********************************************************************************
        //******************** Secondary Fraction     **************************************
        //**********************************************************************************        
        if ( nameMeson.Contains("Pi0")){
            cout << "Plotting secondary fractions" << endl;
            TCanvas* canvasSecFrac = new TCanvas("canvasSecFrac","",200,10,1350,900);  // gives the page size
            DrawGammaCanvasSettings( canvasSecFrac, 0.09, 0.02, 0.04, 0.08);
            canvasSecFrac->SetLogy(1);
            
            Double_t rangeSecRatio[2]       = {0, 1};    
            Double_t rangeSecRatioLin[2]    = {0, 0.15};    
            
            TH2F* histo2DDummySecHad;
            histo2DDummySecHad         = new TH2F("histo2DDummySecHad","histo2DDummySecHad",1000,0, maxPtMeson,
                                                                                    100000, rangeSecRatio[0], rangeSecRatio[1]*100);
            SetStyleHistoTH2ForGraphs(histo2DDummySecHad, "#it{p}_{T} (GeV/#it{c})", "#it{r}_{X} = #frac{X->#pi^{0}}{#pi^{0}}", 0.035   ,0.04, 0.035,0.04, 0.9,1.,510,505);
            histo2DDummySecHad->GetYaxis()->SetRangeUser(0.001,rangeSecRatio[1]);
            histo2DDummySecHad->DrawCopy();         

            TLegend* legendSecFrac  = GetAndSetLegend2(0.65, 0.93-nSecComp*0.035*1.15, 0.94, 0.93, 0.035, 2, "", 42, 0.125);            
            for (Int_t j = 0; j < 4; j++){
                if (haveSec[j]){
                    if (histoTrueSecFracMeson_Or[j]->GetMaximum() > 1e-4){
                        DrawGammaSetMarker(histoTrueSecFracMeson_Or[j], markerStyleSec[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                        histoTrueSecFracMeson_Or[j]->DrawCopy("p,e1,same");  
                        legendSecFrac->AddEntry(histoTrueSecFracMeson_Or[j],Form("#it{r}_{%s}",nameSecMesonPlot[j].Data()));
                        if (fitSecFrac[j]){
                            fitSecFrac[j]->SetLineColor(colorSec[j]);  
                            fitSecFrac[j]->Draw("same");
                            legendSecFrac->AddEntry(fitSecFrac[j],Form("fit to #it{r}_{%s}", nameSecMesonPlot[j].Data()),"l");
                        } else {
                            legendSecFrac->AddEntry((TObject*)0, "","");
                        }
                    }    
                }    
            }    
            legendSecFrac->Draw();

            PutProcessLabelAndEnergyOnPlot(0.14, 0.94, 0.03, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data());
            canvasSecFrac->Update();
            canvasSecFrac->SaveAs(Form("%s/%s_FracSecondaries_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
            
            canvasSecFrac->SetLogy(0);
            histo2DDummySecHad->GetYaxis()->SetRangeUser(rangeSecRatioLin[0],rangeSecRatioLin[1]);
            histo2DDummySecHad->DrawCopy();
            for (Int_t j = 0; j < 4; j++){
                if (haveSec[j]){
                    if (histoTrueSecFracMeson_Or[j]->GetMaximum() > 1e-4){
                        histoTrueSecFracMeson_Or[j]->DrawCopy("p,e1,same");  
                        if (fitSecFrac[j]){
                            fitSecFrac[j]->Draw("same");
                        } 
                    }    
                }    
            }    
            legendSecFrac->Draw();
            PutProcessLabelAndEnergyOnPlot(0.14, 0.94, 0.03, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data());
            
            canvasSecFrac->Update();
            canvasSecFrac->SaveAs(Form("%s/%s_FracSecondariesLin_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
            delete canvasSecFrac;
        }
        
        //**********************************************************************************
        //******************** Efficiency Simple Plot **************************************
        //**********************************************************************************
        cout << "Plotting efficiency" << endl;
        TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasEffSimple, 0.065, 0.01, 0.01, 0.08);
        canvasEffSimple->SetLogy(1);
        
        DrawAutoGammaMesonHistos(   histoTrueEffiPrimMesonPt, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#epsilon_{eff}", 
                                    kTRUE, 2., 3e-5, kFALSE,
                                    kFALSE, 0., 0.3, 
                                    kFALSE, 0., 10.);
        histoTrueEffiPrimMesonPt->GetYaxis()->SetTitleOffset(0.8);        
        DrawGammaSetMarker(histoTrueEffiPrimMesonPt, 20, 1., kBlack, kBlack);
        histoTrueEffiPrimMesonPt->DrawCopy("e1");
        
        TLegend* legendEff = GetAndSetLegend2(0.4,0.125,0.6,0.205, 28);
        legendEff->SetMargin(0.2);
        legendEff->AddEntry(histoTrueEffiPrimMesonPt,Form("prim. %s efficiency",textMeson.Data()));            
        legendEff->Draw();
        
        PutProcessLabelAndEnergyOnPlot(0.7, 0.25, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasEffSimple->Update();
        canvasEffSimple->SaveAs(Form("%s/%s_TrueEffSimple_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data())); 

        // plotting efficiency linearly
        canvasEffSimple->SetLogy(0);

        DrawAutoGammaMesonHistos(   histoTrueEffiPrimMesonPt, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#epsilon_{eff}", 
                                    kTRUE, 0.55, 0, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
        histoTrueEffiPrimMesonPt->GetYaxis()->SetTitleOffset(0.8);        
        DrawGammaSetMarker(histoTrueEffiPrimMesonPt, 20, 1., kBlack, kBlack);
        
        histoTrueEffiPrimMesonPt->DrawCopy("e1");

        TLegend* legendEffLinY = GetAndSetLegend2(0.7,0.125,0.9,0.205, 28);
        legendEffLinY->SetMargin(0.2);
        legendEffLinY->AddEntry(histoTrueEffiPrimMesonPt,Form("prim. %s efficiency",textMeson.Data()));                    
        legendEffLinY->Draw();
        PutProcessLabelAndEnergyOnPlot(0.13, 0.95, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasEffSimple->Update();
        canvasEffSimple->SaveAs(Form("%s/%s_TrueEffSimpleLinY_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data())); 
        
        if (isNewOutput){
            canvasEffSimple->SetLogy(0);

            TLegend* legendEffLinYSec = GetAndSetLegend2(0.13,0.65,0.35,0.8, 28);
            legendEffLinYSec->SetMargin(0.15);
            
            histoTrueEffiPrimMesonPt->GetYaxis()->SetRangeUser(0,3.);
            histoTrueEffiPrimMesonPt->GetYaxis()->SetTitleOffset(0.8);        
            DrawGammaSetMarker(histoTrueEffiPrimMesonPt, 20, 1., kBlack, kBlack);
            histoTrueEffiPrimMesonPt->DrawCopy("e1");
            legendEffLinYSec->AddEntry(histoTrueEffiPrimMesonPt,Form("prim. %s",textMeson.Data()));                    
            
            for (Int_t j = 0; j < 4; j++){
                if (histoSecEffiPt[j]){
                    DrawGammaSetMarker(histoSecEffiPt[j], markerStyleSec[j], markerSizeSec[j], colorSec[j], colorSec[j]);
                    histoSecEffiPt[j]->DrawCopy("e1,same"); 
                    if (!modifiedSecEff[j])
                        legendEffLinYSec->AddEntry(histoSecEffiPt[j],Form("sec #pi^{0} from %s",nameSecMesonPlot[j].Data()),"pe");
                    else 
                        legendEffLinYSec->AddEntry(histoSecEffiPt[j],Form("sec #pi^{0} from %s, adj",nameSecMesonPlot[j].Data()),"pe");
                }    
            }
        
            legendEffLinYSec->Draw();
            PutProcessLabelAndEnergyOnPlot(0.13, 0.95, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
            canvasEffSimple->Update();
            canvasEffSimple->SaveAs(Form("%s/%s_TrueEffSimpleLinYWithSec_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data())); 
        
            canvasEffSimple->cd();
            
            TH2F* histo2DDummyEffiRatio;
            histo2DDummyEffiRatio       = new TH2F("histo2DDummyEffiRatio","histo2DDummyEffiRatio",1000,0, histoTrueEffiPrimMesonPt->GetXaxis()->GetBinUpEdge(histoTrueEffiPrimMesonPt->GetNbinsX()),
                                                                                    1000, 1, 10);
            SetStyleHistoTH2ForGraphs(histo2DDummyEffiRatio, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{sec,eff}/#epsilon_{eff}", 0.035,0.04, 0.035,0.04, 0.9,0.8, 510, 505);
            histo2DDummyEffiRatio->DrawCopy();         

            
            TLegend* legendEffWithSecRatio = GetAndSetLegend2(0.12,0.8-(nEffHistSec)*0.035,0.4,0.8 , 0.035, 1, "", 42, 0.1);
            Bool_t  plotted             = kFALSE;
            for (Int_t j = 0; j < 4; j++){
                if (histoSecEffiPt[j]){
                    plotted                             = kTRUE;
                    DrawGammaSetMarker(histoRatioSecEffDivTrueEff[j],  markerStyleSec[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                    histoRatioSecEffDivTrueEff[j]->DrawCopy("same,e1");  
                    legendEffWithSecRatio->AddEntry(histoRatioSecEffDivTrueEff[j],Form("val. #pi^{0} from %s",nameSecMesonPlot[j].Data()),"p");
                }    
            }    
            legendEffWithSecRatio->Draw();
            
            canvasEffSimple->Update();
            PutProcessLabelAndEnergyOnPlot(0.12, 0.95, 0.035, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());

            if (plotted)
                canvasEffSimple->SaveAs(Form("%s/%s_RatioSecEffiToTrueEff_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));                 
        }    
        delete canvasEffSimple;
    }    
        
    //*********************************************************************************
    //************************** Purity Plot ******************************************
    //*********************************************************************************
    cout << "Plotting purity" << endl;
    TCanvas* canvasPurity = new TCanvas("canvasPurity","",0,0,1000,900);// gives the page size
    DrawGammaCanvasSettings( canvasPurity, 0.09, 0.017, 0.015, 0.08);
    canvasPurity->SetLogy(0);

    DrawAutoGammaMesonHistos( histoMesonPurityPt, 
                                "", "#it{p}_{T} (GeV/#it{c})", "#it{P}_{#pi^{0}}", 
                                kFALSE, 0.75, 3e-6, kFALSE,
                                kTRUE, 0.6, 1.02, 
                                kFALSE, 0., 10.);
    histoMesonPurityPt->GetYaxis()->SetTitleOffset(1.05);        
    DrawGammaSetMarker(histoMesonPurityPt,  20 , 1.5, kAzure+2, kAzure+2);
    histoMesonPurityPt->Draw("e1");
    if (histoMesonPurityUnmodPt){
        DrawGammaSetMarker(histoMesonPurityUnmodPt,  24 , 1.5, kRed-6, kRed-6);  
        histoMesonPurityUnmodPt->Draw("same,e1");
    }    
    TLegend* legendPurity = GetAndSetLegend2(0.2, 0.125, 0.65, 0.205, 28);
    legendPurity->SetMargin(0.12);
    legendPurity->AddEntry(histoMesonPurityPt,Form("%s Purity",textMeson.Data()));
    if (histoMesonPurityUnmodPt)
        legendPurity->AddEntry(histoMesonPurityUnmodPt,Form("%s Purity, uncorr add. #gamma,#eta",textMeson.Data()));
    legendPurity->Draw();   
    
    PutProcessLabelAndEnergyOnPlot(0.68, 0.25, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
    canvasPurity->Update();

    canvasPurity->SaveAs(Form("%s/%s_%s_TruePurity_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
    delete legendPurity;
    delete canvasPurity;
    
     
    if (kIsMC){ 
        //**********************************************************************************
        //********************** BG distribution *******************************************
        //**********************************************************************************
        cout << "Plotting BG yield" << endl;
        TCanvas* canvasBGYield = new TCanvas("canvasBGYield","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasBGYield, 0.07, 0.01, 0.02, 0.08); 
        canvasBGYield->SetLogy(1);

        DrawAutoGammaMesonHistos( histoUnCorrectedYield, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "Yield", 
                                    kTRUE, 4., 4e-10, kTRUE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
        histoUnCorrectedYield->SetLineWidth(0.5); 
        histoUnCorrectedYield->GetYaxis()->SetTitleOffset(0.8);
        DrawGammaSetMarker(histoUnCorrectedYield, 20, 1, kBlack, kBlack);  
        histoUnCorrectedYield->DrawCopy("e1");

        TLegend* legendBG = GetAndSetLegend2(0.55, 0.78, 0.95, 0.95, 28);
        legendBG->SetMargin(0.2);
        legendBG->SetNColumns(4);
        Int_t nLegendBG = 0;
        legendBG->AddEntry(histoUnCorrectedYield,"Signal+BG","p");
        nLegendBG++;
        for (Int_t i = 0; i < 9; i++){
            DrawGammaSetMarker(histoTrueClustersBGPt[i], markerStyleBGPlot[i], 1, colorBGPlot[i], colorBGPlot[i]);  
            histoTrueClustersBGPt[i]->DrawCopy("e1,same");
            legendBG->AddEntry(histoTrueClustersBGPt[i],labelsBGPlot[i], "p");
            nLegendBG++;
            if (nLegendBG % 4 == 0){
                legendBG->AddEntry((TObject*)0,"", "");
                nLegendBG++;
            }    
        }    
        if (kIsEta){
            DrawGammaSetMarker(histoTrueYieldPi0M02, 20, 1, kMagenta+2, kMagenta+2);  
            histoTrueYieldPi0M02->DrawCopy("e1,same");
            legendBG->AddEntry(histoTrueYieldPi0M02,"#pi^{0}","p");
        } else {
            DrawGammaSetMarker(histoTrueYieldEtaM02, 20, 1, kMagenta+2, kMagenta+2);  
            histoTrueYieldEtaM02->DrawCopy("e1,same");            
            legendBG->AddEntry(histoTrueYieldEtaM02,"#eta","p");
        }
        DrawGammaSetMarker(histoTrueYieldGammaM02, 21, 1, kAzure, kAzure);  
        histoTrueYieldGammaM02->DrawCopy("e1,same");            
        legendBG->AddEntry(histoTrueYieldGammaM02,"#gamma","p");
        DrawGammaSetMarker(histoTrueYieldElectronM02, 21, 1, kGreen+3, kGreen+3);  
        histoTrueYieldElectronM02->DrawCopy("e1,same");            
        legendBG->AddEntry(histoTrueYieldElectronM02,"e^{#pm}","p");

        legendBG->Draw();   
        
        PutProcessLabelAndEnergyOnPlot(0.12, 0.95, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasBGYield->Update();

        canvasBGYield->SaveAs(Form("%s/%s_%s_BGYieldPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        
        DrawGammaSetMarker(histoUnCorrectedYield, 20, 1, kBlack, kBlack);  
        histoUnCorrectedYield->DrawCopy("e1");

        TLegend* legendBG2 = GetAndSetLegend2(0.55, 0.85, 0.95, 0.95, 28);
        legendBG2->SetMargin(0.2);
        legendBG2->SetNColumns(2);
        legendBG2->AddEntry(histoUnCorrectedYield,"Signal+BG","p");
        DrawGammaSetMarker(histoTrueClustersBGPt[9], markerStyleBGPlot[9], 1, colorBGPlot[9], colorBGPlot[9]);  
        histoTrueClustersBGPt[9]->DrawCopy("e1,same");
        legendBG2->AddEntry(histoTrueClustersBGPt[9],labelsBGPlot[9], "p");
        if (kIsEta){
            DrawGammaSetMarker(histoTrueYieldPi0M02, 20, 1, kMagenta+2, kMagenta+2);  
            histoTrueYieldPi0M02->DrawCopy("e1,same");
            legendBG2->AddEntry(histoTrueYieldPi0M02,"#pi^{0}","p");
        } else {
            DrawGammaSetMarker(histoTrueYieldEtaM02, 20, 1, kMagenta+2, kMagenta+2);  
            histoTrueYieldEtaM02->DrawCopy("e1,same");            
            legendBG2->AddEntry(histoTrueYieldEtaM02,"#eta","p");
        }
        DrawGammaSetMarker(histoTrueYieldGammaM02, 24, 1, kAzure, kAzure);  
        histoTrueYieldGammaM02->DrawCopy("e1,same");            
        legendBG2->AddEntry(histoTrueYieldGammaM02,"#gamma","p");
        DrawGammaSetMarker(histoTrueYieldElectronM02, 25, 1, kGreen+3, kGreen+3);  
        histoTrueYieldElectronM02->DrawCopy("e1,same");            
        legendBG2->AddEntry(histoTrueYieldElectronM02,"e^{#pm}","p");

        legendBG2->Draw();   
        
        PutProcessLabelAndEnergyOnPlot(0.12, 0.25, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasBGYield->Update();

        canvasBGYield->SaveAs(Form("%s/%s_%s_BGYieldPtCleaner_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        delete canvasBGYield;
    
        cout << "Plotting BG ratio" << endl;
        TCanvas* canvasBGRatio = new TCanvas("canvasBGRatio","",200,10,900,900);  // gives the page size
        DrawGammaCanvasSettings( canvasBGRatio, 0.1, 0.02, 0.02, 0.08); 
        canvasBGRatio->SetLogy(1);

            DrawAutoGammaMesonHistos( histoRatioTrueYieldGammaM02,
                                        "", "#it{p}_{T} (GeV/#it{c})", "#it{c}_{i}", 
                                        kFALSE, 4., 4e-10, kTRUE,
                                        kTRUE, 1e-5, 10, 
                                        kFALSE, 0., 10.);
            DrawGammaSetMarker(histoRatioTrueYieldGammaM02, 21, 1, kAzure, kAzure);  
            histoRatioTrueYieldGammaM02->DrawCopy("e1");

            TLegend* legendBGRatio = GetAndSetLegend2(0.15, 0.78, 0.55, 0.95, 28);
            legendBGRatio->SetMargin(0.2);
            legendBGRatio->SetNColumns(3);
            legendBGRatio->AddEntry(histoRatioTrueYieldGammaM02,"#gamma","p");

            if (kIsEta){
                DrawGammaSetMarker(histoRatioTrueYieldPi0M02, 20, 1, kMagenta+2, kMagenta+2);  
                histoRatioTrueYieldPi0M02->DrawCopy("e1,same");
                legendBGRatio->AddEntry(histoRatioTrueYieldPi0M02,"#pi^{0}","p");
            } else {
                DrawGammaSetMarker(histoRatioTrueYieldEtaM02, 20, 1, kMagenta+2, kMagenta+2);  
                histoRatioTrueYieldEtaM02->DrawCopy("e1,same");            
                legendBGRatio->AddEntry(histoRatioTrueYieldEtaM02,"#eta","p");
            }
            DrawGammaSetMarker(histoRatioTrueYieldElectronM02, 21, 1, kGreen+3, kGreen+3);  
            histoRatioTrueYieldElectronM02->DrawCopy("e1,same");            
            legendBGRatio->AddEntry(histoRatioTrueYieldElectronM02,"e^{#pm}","p");
            
            for (Int_t i = 0; i < 9; i++){
                DrawGammaSetMarker(histoRatioTrueClustersBGPt[i], markerStyleBGPlot[i], 1, colorBGPlot[i], colorBGPlot[i]);  
                histoRatioTrueClustersBGPt[i]->DrawCopy("e1,same");
                legendBGRatio->AddEntry(histoRatioTrueClustersBGPt[i],labelsBGPlot[i], "p");
            }    
            legendBGRatio->Draw();   
            
            PutProcessLabelAndEnergyOnPlot(0.14, 0.25, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasBGRatio->Update();
        canvasBGRatio->SaveAs(Form("%s/%s_%s_BGRatioPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));

            DrawAutoGammaMesonHistos( histoRatioTrueClustersBGPt[0],
                                        "", "#it{p}_{T} (GeV/#it{c})", "#it{c}_{i}", 
                                        kFALSE, 4., 4e-10, kTRUE,
                                        kTRUE, 1e-5, 10, 
                                        kFALSE, 0., 10.);
            DrawGammaSetMarker(histoRatioTrueClustersBGPt[0], markerStyleBGPlot[0], 1, colorBGPlot[0], colorBGPlot[0]);  
            histoRatioTrueClustersBGPt[0]->DrawCopy("e1");

            TLegend* legendBGRatio3 = GetAndSetLegend2(0.15, 0.83, 0.55, 0.95, 28);
            legendBGRatio3->SetMargin(0.2);
            legendBGRatio3->SetNColumns(3);
            legendBGRatio3->AddEntry(histoRatioTrueClustersBGPt[0],labelsBGPlot[0],"p");

            for (Int_t i = 1; i < 9; i++){
                DrawGammaSetMarker(histoRatioTrueClustersBGPt[i], markerStyleBGPlot[i], 1, colorBGPlot[i], colorBGPlot[i]);  
                histoRatioTrueClustersBGPt[i]->DrawCopy("e1,same");
                legendBGRatio3->AddEntry(histoRatioTrueClustersBGPt[i],labelsBGPlot[i], "p");
            }
            legendBGRatio3->Draw();   
            
            PutProcessLabelAndEnergyOnPlot(0.14, 0.25, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasBGRatio->Update();
        canvasBGRatio->SaveAs(Form("%s/%s_%s_BGRatioPtPureBG_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        
            DrawAutoGammaMesonHistos(   histoRatioTrueYieldGammaM02,
                                        "", "#it{p}_{T} (GeV/#it{c})", "#it{c}_{i}", 
                                        kFALSE, 4., 4e-10, kTRUE,
                                        kTRUE, 1e-5, 10, 
                                        kFALSE, 0., 10.);
            DrawGammaSetMarker(histoRatioTrueYieldGammaM02, 24, 1, kAzure, kAzure);  
            histoRatioTrueYieldGammaM02->DrawCopy("e1");

            TLegend* legendBGRatio2 = GetAndSetLegend2(0.15, 0.87, 0.55, 0.95, 28);
            legendBGRatio2->SetMargin(0.2);
            legendBGRatio2->SetNColumns(2);
            legendBGRatio2->AddEntry(histoRatioTrueYieldGammaM02,"#gamma","p");

            if (kIsEta){
                DrawGammaSetMarker(histoRatioTrueYieldPi0M02, 20, 1, kMagenta+2, kMagenta+2);  
                histoRatioTrueYieldPi0M02->DrawCopy("e1,same");
                legendBGRatio2->AddEntry(histoRatioTrueYieldPi0M02,"#pi^{0}","p");
            } else {
                DrawGammaSetMarker(histoRatioTrueYieldEtaM02, 20, 1, kMagenta+2, kMagenta+2);  
                histoRatioTrueYieldEtaM02->DrawCopy("e1,same");            
                legendBGRatio2->AddEntry(histoRatioTrueYieldEtaM02,"#eta","p");
            }
            DrawGammaSetMarker(histoRatioTrueYieldElectronM02, 25, 1, kGreen+3, kGreen+3);
            histoRatioTrueYieldElectronM02->DrawCopy("e1,same");            
            legendBGRatio2->AddEntry(histoRatioTrueYieldElectronM02,"e^{#pm}","p");
            
        
            DrawGammaSetMarker(histoRatioTrueClustersBGPt[9], markerStyleBGPlot[9], 1, colorBGPlot[9], colorBGPlot[9]);  
            histoRatioTrueClustersBGPt[9]->DrawCopy("e1,same");
            legendBGRatio2->AddEntry(histoRatioTrueClustersBGPt[9],labelsBGPlot[9], "p");
            legendBGRatio2->Draw();   
            
            PutProcessLabelAndEnergyOnPlot(0.14, 0.25, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasBGRatio->Update();

        canvasBGRatio->SaveAs(Form("%s/%s_%s_BGRatioPtCleaner_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));

        canvasBGRatio->SetLogy(0);

            DrawAutoGammaMesonHistos(   histoRatioTrueYieldGammaM02,
                                "", "#it{p}_{T} (GeV/#it{c})", "#it{c}_{i}", 
                                kTRUE, 0.02, 0, kTRUE,
                                kFALSE, 0, 0.2, 
                                kFALSE, 0., 10.);
            DrawGammaSetMarker(histoRatioTrueYieldGammaM02, 24, 1, kAzure, kAzure);  
            histoRatioTrueYieldGammaM02->DrawCopy("e1");

            if (kIsEta){
                histoRatioTrueYieldPi0M02->DrawCopy("e1,same");
            } else {
                histoRatioTrueYieldEtaM02->DrawCopy("e1,same");            
            }
            histoRatioTrueYieldElectronM02->DrawCopy("e1,same");                    
            histoRatioTrueClustersBGPt[9]->DrawCopy("e1,same");
            legendBGRatio2->Draw();   
            
            PutProcessLabelAndEnergyOnPlot(0.14, 0.25, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);

        canvasBGRatio->Update();
        canvasBGRatio->SaveAs(Form("%s/%s_%s_BGRatioPtCleanerLinY_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));

            DrawGammaSetMarker(histoRatioTrueYieldGammaM02, 21, 1, kAzure, kAzure);  
            histoRatioTrueYieldGammaM02->DrawCopy("e1");

            if (kIsEta){
                DrawGammaSetMarker(histoRatioTrueYieldPi0M02, 20, 1, kMagenta+2, kMagenta+2);  
                histoRatioTrueYieldPi0M02->DrawCopy("e1,same");
            } else {
                DrawGammaSetMarker(histoRatioTrueYieldEtaM02, 20, 1, kMagenta+2, kMagenta+2);  
                histoRatioTrueYieldEtaM02->DrawCopy("e1,same");            
            }
            DrawGammaSetMarker(histoRatioTrueYieldElectronM02, 21, 1, kGreen+3, kGreen+3);  
            histoRatioTrueYieldElectronM02->DrawCopy("e1,same");            
            
            for (Int_t i = 0; i < 9; i++){
                DrawGammaSetMarker(histoRatioTrueClustersBGPt[i], markerStyleBGPlot[i], 1, colorBGPlot[i], colorBGPlot[i]);  
                histoRatioTrueClustersBGPt[i]->DrawCopy("e1,same");
            }    
            legendBGRatio->Draw();   
            
            PutProcessLabelAndEnergyOnPlot(0.70, 0.94, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasBGRatio->Update();
        canvasBGRatio->SaveAs(Form("%s/%s_%s_BGRatioPtLinY_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));

        //**********************************************************************************
        //********************** Plot double counting fraction     ************************
        //**********************************************************************************
        TCanvas* canvasDCFrac = new TCanvas("canvasDCFrac","",200,10,900,900);  // gives the page size
        DrawGammaCanvasSettings( canvasDCFrac, 0.1, 0.02, 0.02, 0.08); 
        canvasDCFrac->SetLogy(0);

            histoRatioPi0DCFrac->Scale(100);
            DrawAutoGammaMesonHistos( histoRatioPi0DCFrac,
                                        "", "#it{p}_{T} (GeV/#it{c})", "D = double counted #pi^{0}/ all reconstructed #pi^{0} (%)", 
                                        kFALSE, 4., 4e-10, kTRUE,
                                        kFALSE, 1e-5, 10, 
                                        kFALSE, 0., 10.);
            DrawGammaSetMarker(histoRatioPi0DCFrac, 21, 1, kAzure, kAzure);  
            histoRatioPi0DCFrac->DrawCopy("e1");
            
            PutProcessLabelAndEnergyOnPlot(0.14, 0.95, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasDCFrac->Update();
        canvasDCFrac->SaveAs(Form("%s/%s_%s_Pi0DCFrac_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));

            histoRatioEtaDCFrac->Scale(100);
            DrawAutoGammaMesonHistos( histoRatioEtaDCFrac,
                                        "", "#it{p}_{T} (GeV/#it{c})", "D = double counted #eta/ all reconstructed #eta (%)", 
                                        kFALSE, 4., 4e-10, kTRUE,
                                        kFALSE, 1e-5, 10, 
                                        kFALSE, 0., 10.);
            DrawGammaSetMarker(histoRatioEtaDCFrac, 21, 1, kAzure, kAzure);  
            histoRatioEtaDCFrac->DrawCopy("e1");
            
            PutProcessLabelAndEnergyOnPlot(0.14, 0.95, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasDCFrac->Update();
        canvasDCFrac->SaveAs(Form("%s/%s_%s_EtaDCFrac_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));

        //**********************************************************************************
        //********************** Plot decay channel decompositions     ************************
        //**********************************************************************************
        TCanvas* canvasDecayChannelFrac = new TCanvas("canvasDecayChannelFrac","",200,10,900,900);  // gives the page size
        DrawGammaCanvasSettings( canvasDecayChannelFrac, 0.1, 0.02, 0.02, 0.08); 
        canvasDecayChannelFrac->SetLogy(1);

            histoRatioPi0GGFrac->Scale(100);
            DrawAutoGammaMesonHistos( histoRatioPi0GGFrac,
                                        "", "#it{p}_{T} (GeV/#it{c})", "K_{X} = rec. #pi^{0} #rightarrow X / all rec. #pi^{0} (%)", 
                                        kFALSE, 4., 4e-10, kTRUE,
                                        kTRUE, 0.01, 205, 
                                        kFALSE, 0., 10.);
            DrawGammaSetMarker(histoRatioPi0GGFrac, 21, 1, kAzure, kAzure);  
            histoRatioPi0GGFrac->DrawCopy("e1");      
            histoRatioPi0DalitzFrac->Scale(100);
            DrawGammaSetMarker(histoRatioPi0DalitzFrac, 20, 1, kRed+2, kRed+2);  
            histoRatioPi0DalitzFrac->DrawCopy("e1,same");

            TLegend* legendDecayChannelDecomp = GetAndSetLegend2(0.65, 0.12, 0.95, 0.2, 28);
            legendDecayChannelDecomp->SetMargin(0.2);
            legendDecayChannelDecomp->AddEntry(histoRatioPi0GGFrac,"K_{#gamma#gamma}", "p");
            legendDecayChannelDecomp->AddEntry(histoRatioPi0DalitzFrac,"K_{#gamma e^{+}e^{-}}", "p");
            legendDecayChannelDecomp->Draw();
            
            PutProcessLabelAndEnergyOnPlot(0.14, 0.23, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasDecayChannelFrac->Update();
        canvasDecayChannelFrac->SaveAs(Form("%s/%s_%s_Pi0DecayChannelDecomposition_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));

            histoRatioEtaGGFrac->Scale(100);
            DrawAutoGammaMesonHistos( histoRatioEtaGGFrac,
                                        "", "#it{p}_{T} (GeV/#it{c})", "K_{X} = rec. #eta #rightarrow X / all rec. #eta (%)", 
                                        kFALSE, 4., 4e-10, kTRUE,
                                        kTRUE, 0.01, 205, 
                                        kFALSE, 0., 10.);
            DrawGammaSetMarker(histoRatioEtaGGFrac, 21, 1, kAzure, kAzure);  
            histoRatioEtaGGFrac->DrawCopy("e1");
            histoRatioEtaDalitzFrac->Scale(100);
            DrawGammaSetMarker(histoRatioEtaDalitzFrac, 20, 1, kRed+2, kRed+2);  
            histoRatioEtaDalitzFrac->DrawCopy("e1,same");
            legendDecayChannelDecomp->Draw();
            
            PutProcessLabelAndEnergyOnPlot(0.14, 0.23, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasDecayChannelFrac->Update();
        canvasDecayChannelFrac->SaveAs(Form("%s/%s_%s_EtaDecayChannelDecomposition_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));

        //**********************************************************************************
        //********************** Plot merged cluster definition decomposition **************
        //**********************************************************************************
        TCanvas* canvasClusterClassificationFrac = new TCanvas("canvasClusterClassificationFrac","",200,10,900,900);  // gives the page size
        DrawGammaCanvasSettings( canvasClusterClassificationFrac, 0.1, 0.02, 0.02, 0.08); 
        canvasClusterClassificationFrac->SetLogy(1);

            
            DrawAutoGammaMesonHistos( histoRatioMergedPureFrac,
                                        "", "#it{p}_{T} (GeV/#it{c})", "L_{X} = clus. rec from X / all clus. (%)", 
                                        kFALSE, 4., 4e-10, kTRUE,
                                        kTRUE, 3, 205, 
                                        kFALSE, 0., 10.);
            DrawGammaSetMarker(histoRatioMergedPureFrac, 21, 1, kAzure, kAzure);  
            histoRatioMergedPureFrac->DrawCopy("e1");      
            DrawGammaSetMarker(histoRatioMergedPartConvFrac, 20, 1, kRed+2, kRed+2);  
            histoRatioMergedPartConvFrac->DrawCopy("e1,same");
            DrawGammaSetMarker(histoRatioMergedOneGammaFrac, 25, 1, kGreen+2, kGreen+2);  
            histoRatioMergedOneGammaFrac->DrawCopy("e1,same");
            DrawGammaSetMarker(histoRatioMergedOneElectronFrac, 24, 1, kViolet+2, kViolet+2);  
            histoRatioMergedOneElectronFrac->DrawCopy("e1,same");

            TLegend* legendClusterClassificationDecomp = GetAndSetLegend2(0.65, 0.82, 0.95, 0.95, 28);
            legendClusterClassificationDecomp->SetMargin(0.2);
            legendClusterClassificationDecomp->AddEntry(histoRatioMergedPureFrac,"L_{merged}", "p");
            legendClusterClassificationDecomp->AddEntry(histoRatioMergedPartConvFrac,"L_{merged part. conv}", "p");
            legendClusterClassificationDecomp->AddEntry(histoRatioMergedOneGammaFrac,"L_{1#gamma from decay}", "p");
            legendClusterClassificationDecomp->AddEntry(histoRatioMergedOneElectronFrac,"L_{1e^{#pm} from decay}", "p");
            legendClusterClassificationDecomp->Draw();
            
            for (Int_t n=1; n< histoRatioMergedPureFrac->GetNbinsX()+1; n++){
                Double_t all     = histoRatioMergedPureFrac->GetBinContent(n) + histoRatioMergedPartConvFrac->GetBinContent(n) + histoRatioMergedOneGammaFrac->GetBinContent(n) + 
                                   histoRatioMergedOneElectronFrac->GetBinContent(n);
                if (all != 0){
                    cout << histoRatioMergedPureFrac->GetBinCenter(n) << ": \t"<< histoRatioMergedPureFrac->GetBinContent(n) << "\t" << histoRatioMergedPartConvFrac->GetBinContent(n) << "\t" 
                        << histoRatioMergedOneGammaFrac->GetBinContent(n) << "\t" << histoRatioMergedOneElectronFrac->GetBinContent(n) << "\t" << all << endl;
                }        
            }
            PutProcessLabelAndEnergyOnPlot(0.14, 0.95, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasClusterClassificationFrac->Update();
        canvasClusterClassificationFrac->SaveAs(Form("%s/%s_%s_ClusterClassificationDecomposition_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        
    }    
    
    
    //**********************************************************************************
    //*************************** MC Yield *********************************************
    //***** need to do it for MC and data in order to have it in the output file  ******
    //**********************************************************************************

    cout << "correct MC yield" << endl;
    TH1D *histoMCYieldMesonOldBin = (TH1D*)histoInputMesonOldBinPt->Clone();
    histoMCYieldMesonOldBin->SetName("MCYield_Meson_oldBin");
    ScaleMCYield(histoMCYieldMesonOldBin,  deltaRapid,  scaling,  nEvtMC,  nameMeson );
    Float_t integralMB = 0;
    integralMB = histoMCYieldMesonOldBin->Integral(histoMCYieldMesonOldBin->FindBin(8),histoMCYieldMesonOldBin->FindBin(10));

    TH1D *histoMCYieldMeson = (TH1D*)histoInputMesonPt->Clone();
    histoMCYieldMeson->SetName("MCYield_Meson");
    ScaleMCYield(histoMCYieldMeson,  deltaRapid,  scaling,  nEvtMC,  nameMeson );
        
    //**********************************************************************************
    //******************** RAW Yield spectrum ******************************************
    //**********************************************************************************
    cout << "plotting raw yield" << endl;    
    TCanvas* canvasRAWYield = new TCanvas("canvasRAWYield","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRAWYield, 0.07, 0.01, 0.02, 0.08); 
    canvasRAWYield->SetLogy(1);

    TH1D* histoUnCorrectedYieldDrawing = (TH1D*)histoUnCorrectedYield->Clone();
    histoUnCorrectedYieldDrawing->Multiply(histoMesonPurityPt);
    histoUnCorrectedYieldDrawing->Scale(1./nEvt);
    DrawAutoGammaMesonHistos( histoUnCorrectedYieldDrawing, 
                                "", "#it{p}_{T} (GeV/#it{c})", "RAW Yield/ #it{N}_{Evt}", 
                                kTRUE, 2.5, 0.1/nEvt, kTRUE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    histoUnCorrectedYieldDrawing->SetLineWidth(0.5); 
    histoUnCorrectedYieldDrawing->GetYaxis()->SetTitleOffset(0.8);
    DrawGammaSetMarker(histoUnCorrectedYieldDrawing, 20, 1, kBlack, kBlack);  
    histoUnCorrectedYieldDrawing->DrawCopy("e1");

    PutProcessLabelAndEnergyOnPlot(0.68, 0.95, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
    canvasRAWYield->Update();

    canvasRAWYield->SaveAs(Form("%s/%s_%s_RAWYieldPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
    delete canvasRAWYield;

    //***********************************************************************************************
    //***************************  Secondary RAW Yield  *********************************************
    //***********************************************************************************************    
    if ( nameMeson.Contains("Pi0") ){ //&& !kIsMC
        TCanvas* canvasRAWYieldSec              = new TCanvas("canvasRAWYieldSec","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRAWYieldSec, 0.07, 0.01, 0.02, 0.08); 
        canvasRAWYieldSec->SetLogy(1);
       
        Int_t nColumnsSec           = 1;
        Int_t nSecPlot              = nSecCompUsed;
        if (!kIsMC){
            nColumnsSec             = 2;
            nSecPlot                = 4;
        }

        // calculate ratio to Uncorrected yi
        for (Int_t j = 0; j < 4; j++){
            if (histoYieldSecMeson[j]){
                histoRatioYieldSecMeson[j]    = (TH1D*)histoYieldSecMeson[j]->Clone(Form("RatioSecYieldFrom%sMesonToRaw", nameSecMeson[j].Data()));
                histoRatioYieldSecMeson[j]->Sumw2();
                histoRatioYieldSecMeson[j]->Divide(histoRatioYieldSecMeson[j],histoUnCorrectedYieldDrawing);
            }    
            if (j < 3){
                if (histoYieldSecMesonFromToy[j]){
                    histoRatioYieldSecMesonFromToy[j]    = (TH1D*)histoYieldSecMesonFromToy[j]->Clone(Form("RatioSecYieldFrom%sMesonFromToyToRaw", nameSecMeson[j].Data()));
                    histoRatioYieldSecMesonFromToy[j]->Sumw2();
                    histoRatioYieldSecMesonFromToy[j]->Divide(histoRatioYieldSecMesonFromToy[j],histoUnCorrectedYieldDrawing);
                }    
            }    
        }
        
        
        DrawAutoGammaMesonHistos( histoUnCorrectedYieldDrawing, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "RAW #pi^{0} Yield/ #it{N}_{Evt}", 
                                    kTRUE, 1, 0.0001/nEvt, kTRUE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
        histoUnCorrectedYieldDrawing->SetLineWidth(0.5); 
        histoUnCorrectedYieldDrawing->GetYaxis()->SetTitleOffset(0.8);
        DrawGammaSetMarker(histoUnCorrectedYieldDrawing, 20, 1., kBlack, kBlack);
        histoUnCorrectedYieldDrawing->Draw("e1");

        
        TLegend* legendSecRAWYield  = GetAndSetLegend2(0.65,0.93-(nSecPlot+1)*0.035,0.93,0.93, 0.035, nColumnsSec, "", 42, 0.12);
        for (Int_t j = 0; j < 4; j++){
            if (histoYieldSecMeson[j] && haveSecUsed[j]){
                DrawGammaSetMarker(histoYieldSecMeson[j],  markerStyleSecWithToy[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                histoYieldSecMeson[j]->DrawCopy("same,e1");  
                legendSecRAWYield->AddEntry(histoYieldSecMeson[j],Form("#pi^{0} from %s",nameSecMesonPlot[j].Data()),"p");
            } else if (!kIsMC && j < 3){
                if (histoYieldSecMesonFromToy[j] && histoYieldSecMesonFromToy[j]->GetEntries()){
                    if (!kIsMC) legendSecRAWYield->AddEntry((TObject*)0,Form("#pi^{0} from %s",nameSecMesonPlot[j].Data()),"");
                }    
            }    
            if (j < 3){
                if (histoYieldSecMesonFromToy[j] && !kIsMC ){
                    if (histoYieldSecMesonFromToy[j]->GetEntries() > 0){
                        DrawGammaSetMarker(histoYieldSecMesonFromToy[j],  markerStyleSecFromToy[j] , markerSizeSecFromToy[j], colorSecFromToy[j], colorSecFromToy[j]);  
                        histoYieldSecMesonFromToy[j]->DrawCopy("same,e1");  
                        legendSecRAWYield->AddEntry(histoYieldSecMesonFromToy[j],"Toy appr.","p");
                    } else if (histoYieldSecMeson[j] && !kIsMC && haveSecUsed[j]){
                        legendSecRAWYield->AddEntry((TObject*)0,"","");
                    } 
                } else if (histoYieldSecMeson[j] && !kIsMC && haveSecUsed[j]){
                    legendSecRAWYield->AddEntry((TObject*)0,"","");
                } 
            }             
        }    
        legendSecRAWYield->Draw();
        PutProcessLabelAndEnergyOnPlot(0.13, 0.94, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        
        canvasRAWYieldSec->Update();
        canvasRAWYieldSec->SaveAs(Form("%s/%s_%s_RAWYieldSecPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        delete canvasRAWYieldSec;
        
        if ( histoYieldSecMesonFromToy[0] && !kIsMC){
            TCanvas* canvasSecFrac2 = new TCanvas("canvasSecFrac","",200,10,1350,900);  // gives the page size
            DrawGammaCanvasSettings( canvasSecFrac2, 0.09, 0.018, 0.04, 0.08);

            Double_t rangeSecRatio[2]   = {0, 0.2};    
            
            TH2F* histo2DDummySecHad2;
            histo2DDummySecHad2         = new TH2F("histo2DDummySecHad2","histo2DDummySecHad2",1000,0, maxPtMeson,
                                                                                    100000, rangeSecRatio[0], rangeSecRatio[1]*100);
            SetStyleHistoTH2ForGraphs(histo2DDummySecHad2, "#it{p}_{T} (GeV/#it{c})", "#it{r}_{X} = #frac{X->#pi^{0}}{#pi^{0}}", 0.035   ,0.04, 0.035,0.04, 0.9,1.,510,505);
            histo2DDummySecHad2->GetYaxis()->SetRangeUser(0,rangeSecRatio[1]);
            histo2DDummySecHad2->DrawCopy();         
            
            TLegend* legendSecRAWRatio  = GetAndSetLegend2(0.65,0.93-(nSecPlot)*0.035,0.93,0.93, 0.035, nColumnsSec, "", 42, 0.12);
            for (Int_t j = 0; j < 4; j++){
                if (histoRatioYieldSecMeson[j] && haveSecUsed[j]){
                    DrawGammaSetMarker(histoRatioYieldSecMeson[j],  markerStyleSecWithToy[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                    histoRatioYieldSecMeson[j]->DrawCopy("same,e1");  
                    legendSecRAWRatio->AddEntry(histoRatioYieldSecMeson[j],Form("#pi^{0} from %s",nameSecMesonPlot[j].Data()),"p");
                } else if (!kIsMC && j < 3){
                    if (histoRatioYieldSecMesonFromToy[j] && histoRatioYieldSecMesonFromToy[j]->GetEntries()){
                        if (!kIsMC) legendSecRAWRatio->AddEntry((TObject*)0,Form("#pi^{0} from %s",nameSecMesonPlot[j].Data()),"");
                    }    
                }    
                if (j < 3){
                    if (histoRatioYieldSecMesonFromToy[j] && !kIsMC ){
                        if (histoRatioYieldSecMesonFromToy[j]->GetEntries() > 0){
                            DrawGammaSetMarker(histoRatioYieldSecMesonFromToy[j],  markerStyleSecFromToy[j] , markerSizeSecFromToy[j], colorSecFromToy[j], colorSecFromToy[j]);  
                            histoRatioYieldSecMesonFromToy[j]->DrawCopy("same,e1");  
                            legendSecRAWRatio->AddEntry(histoRatioYieldSecMesonFromToy[j],"Toy appr.","p");
                        } else if (histoRatioYieldSecMeson[j] && !kIsMC && haveSecUsed[j]){
                            legendSecRAWRatio->AddEntry((TObject*)0,"","");
                        } 
                    } else if (histoRatioYieldSecMeson[j] && !kIsMC && haveSecUsed[j]){
                        legendSecRAWRatio->AddEntry((TObject*)0,"","");
                    } 
                }     
            }
            legendSecRAWRatio->Draw();
            PutProcessLabelAndEnergyOnPlot(0.15, 0.93, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            
            canvasSecFrac2->Update();
            canvasSecFrac2->SaveAs(Form("%s/%s_%s_EffectiveSecCorrPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        }    

        
    }
    
    //***********************************************************************************************
    //*********************************** correction for yield **************************************
    //***********************************************************************************************
    cout << "calculating corrected yield" << endl;    
    
    TH1D* histoCorrectedYieldTrue       = (TH1D*)histoUnCorrectedYield->Clone();
    histoCorrectedYieldTrue->SetName("CorrectedYieldTrueEff");
    
    CorrectYield(histoCorrectedYieldTrue, histoMesonPurityPt, histoYieldSecMeson, histoYieldSecMesonFromToy, histoTrueEffiPrimMesonPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson, kIsMC);
    
    // **************************************************************************************
    // ************** Plot corrected yield with differnt yield extraction methods ***********
    // **************************************************************************************
    cout << "plotting corrected yield" << endl;    

    TCanvas* canvasCorrecftedYield = new TCanvas("canvasCorrecftedYield","",1350,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasCorrecftedYield, 0.13, 0.02, 0.02, 0.09);
    canvasCorrecftedYield->SetLogy();

    TPad* padCorrectedYieldHistos = new TPad("padCorrectedYieldHistos", "", 0., 0.3, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padCorrectedYieldHistos, 0.12, 0.02, 0.02, 0.);
    padCorrectedYieldHistos->Draw();

    TPad* padCorrectedYieldRatios = new TPad("padCorrectedYieldRatios", "", 0., 0., 1., 0.3,-1, -1, -2);
    DrawGammaPadSettings( padCorrectedYieldRatios, 0.12, 0.02, 0., 0.21);
    padCorrectedYieldRatios->Draw();

    padCorrectedYieldHistos->cd();
    padCorrectedYieldHistos->SetLogy(); 

    DrawAutoGammaMesonHistos( histoCorrectedYieldTrue, 
                                "", "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 
                                kTRUE, 3., 4e-10, kTRUE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    DrawGammaSetMarker(histoCorrectedYieldTrue, 20, 1., kBlack, kBlack);  
    histoCorrectedYieldTrue->DrawCopy("e1");  
//     DrawGammaSetMarker(histoMCYieldMeson, 24, 1., kRed+2, kRed+2);  
//     histoMCYieldMeson->DrawCopy("same,e1");  

    TLegend* legendYield3 = GetAndSetLegend2(0.15, 0.03, 0.66, 0.19, 28);
    legendYield3->SetMargin(0.15);
    legendYield3->AddEntry(histoCorrectedYieldTrue,"corrected yield");
//     legendYield3->AddEntry(histoMCYieldMeson,"MC input yield");
    legendYield3->Draw();

    PutProcessLabelAndEnergyOnPlot(0.62, 0.95, 0.03, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data());
    
    padCorrectedYieldRatios->cd();
    padCorrectedYieldRatios->SetTickx();
    padCorrectedYieldRatios->SetTicky();
    padCorrectedYieldRatios->SetLogy(0);
    TH1D *ratioTrue         = (TH1D*) histoCorrectedYieldTrue->Clone(); 
    ratioTrue->Divide(ratioTrue,histoCorrectedYieldTrue,1.,1.,"");
    TH1D *ratioTrueMCInput  = (TH1D*) histoMCYieldMeson->Clone();
    ratioTrueMCInput->Divide(histoMCYieldMeson, histoCorrectedYieldTrue,1.,1.,"");

    ratioTrue->SetYTitle("#frac{modified}{standard}"); 
    SetStyleHistoTH1ForGraphs(ratioTrue, "#it{p}_{T} (GeV/#it{c})","#frac{modified}{standard}", 0.07,0.1, 0.07,0.11, 0.9,0.4, 510, 505);
    ratioTrue->GetYaxis()->SetRangeUser(0.8,1.23);
    DrawGammaSetMarker(ratioTrue, 20, 1., kBlack, kBlack); 
    ratioTrue->DrawCopy("p,e1");  
//     DrawGammaSetMarker(ratioTrueMCInput, 24, 1., kRed+2, kRed+2); 
//     ratioTrueMCInput->DrawCopy("same,p,e1");  
    
    canvasCorrecftedYield->Update();
    canvasCorrecftedYield->SaveAs(Form("%s/%s_%s_CorrectedYield_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(),  fCutSelection.Data(), suffix.Data()));

    // **************************************************************************************
    // ************** Plot corrected yield with differnt efficiencies & MC yield ************
    // **************************** Sanity check for MC *************************************
    // **************************************************************************************
    if (kIsMC){
        cout << "plotting sanity check" << endl;    
    
        canvasCorrecftedYield->cd();

        padCorrectedYieldHistos->cd();
        padCorrectedYieldHistos->SetLogy(); 

        DrawGammaSetMarker(histoCorrectedYieldTrue, 20, 1., kBlack, kBlack);  
        histoCorrectedYieldTrue->DrawCopy("e1");  
        DrawGammaSetMarker(histoMCYieldMeson, 24, 1., kRed+2, kRed+2);  
        histoMCYieldMeson->DrawCopy("e1,same"); 
        
        cout << "here" << endl; 
        TLegend* legendYield4 = GetAndSetLegend2(0.15, 0.03, 0.66, 0.19, 28);
        legendYield4->SetMargin(0.15);
        legendYield4->AddEntry(histoCorrectedYieldTrue,"corr true eff");
        legendYield4->AddEntry(histoMCYieldMeson,"MC input (possibly weighted)");
        legendYield4->Draw();

        PutProcessLabelAndEnergyOnPlot(0.62, 0.95, 0.03, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data());
        padCorrectedYieldRatios->cd();
        padCorrectedYieldRatios->SetTickx();
        padCorrectedYieldRatios->SetTicky();
        padCorrectedYieldRatios->SetLogy(0);
        
        DrawGammaSetMarker(ratioTrue, 20, 1., kBlack, kBlack);
        for(Int_t b = 0; b< ratioTrue->GetNbinsX(); b++){
            ratioTrue->SetBinError(b+1,histoCorrectedYieldTrue->GetBinError(b+1)/histoCorrectedYieldTrue->GetBinContent(b+1));
        }
        ratioTrue->SetFillColor(kGray+2);
        ratioTrue->SetFillStyle(1);
        ratioTrue->DrawCopy("p,e2");
        
        DrawGammaSetMarker(ratioTrueMCInput, 24, 1., kRed+2, kRed+2);
        ratioTrueMCInput->DrawCopy("e1,same"); 

        canvasCorrecftedYield->Update();
        canvasCorrecftedYield->SaveAs(Form("%s/%s_%s_CorrectedYield_SanityCheck_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(),  fCutSelection.Data(), suffix.Data()));
    }

    delete canvasCorrecftedYield;
    delete legendYield3;
    
    // ********************************************************************************************************************************
    // **************************************** plotting NLM fractions ****************************************************************
    // ********************************************************************************************************************************
    if (enableNLMBinsPlotting && !kIsMC){
        TCanvas* canvasNLMFractions = new TCanvas("canvasNLMFractions","",200,10,900,900);  // gives the page size
        DrawGammaCanvasSettings( canvasNLMFractions, 0.1, 0.02, 0.02, 0.08); 
//         canvasNLMFractions->SetLogy(1);

        Color_t colorNLM[11]        = { kBlack, kRed+1, kBlue+1, kGreen+2, kViolet+1, 
                                        kOrange, kCyan+2, kMagenta+1, kOrange+7, kGray+1,
                                        kPink-8 };
        Style_t markerStyleNLM[11]  = { 20, 21, 29, 34, 33, 
                                        20, 21, 29, 34, 33,
                                        20 };
        Style_t markerStyleNLMMC[11]= { 24, 25, 28, 28, 27, 
                                        24, 25, 28, 28, 27,
                                        24 };                                   
        TH2F* histo2DDummyNLMFrac;
        histo2DDummyNLMFrac         = new TH2F("histo2DDummyNLMFrac","histo2DDummyNLMFrac",1000,0, maxPtMeson,
                                                                                100000, 0, 1.2);
            SetStyleHistoTH2ForGraphs(histo2DDummyNLMFrac, "#it{p}_{T} (GeV/#it{c})", "#it{f}_{X} = cluster_{NLM=X}/cluster_{all NLM}", 0.035   ,0.04, 0.035,0.04, 0.9,1.,510,505);
//             histo2DDummyNLMFrac->GetYaxis()->SetRangeUser(0,1);
            histo2DDummyNLMFrac->DrawCopy();         
        
            TLegend* legendNLMs = GetAndSetLegend2(0.13, 0.96-4*0.035, 0.45, 0.96, 28,3);
            for (Int_t i = 0; i < 3; i++){
                if (histoRatioPtInNLMBins[i]){
                    legendNLMs->AddEntry((TObject*)0,Form("X=%d",i+1) ,"");
                    DrawGammaSetMarker(histoRatioPtInNLMBins[i], markerStyleNLM[i], 1, colorNLM[i], colorNLM[i]);  
                    histoRatioPtInNLMBins[i]->DrawCopy("same, e1");      
                    legendNLMs->AddEntry(histoRatioPtInNLMBins[i],"Data" ,"p");
                    if (histoMCRatioPtInNLMBins[i]){
                        DrawGammaSetMarker(histoMCRatioPtInNLMBins[i], markerStyleNLMMC[i], 1, colorNLM[i], colorNLM[i]);  
                        histoMCRatioPtInNLMBins[i]->DrawCopy("e1,same");
                        legendNLMs->AddEntry(histoMCRatioPtInNLMBins[i],"MC","p");
                    } else {
                        legendNLMs->AddEntry((TObject*)0,"","");
                    }    
                }
                
            }    
            if (histoRatioPtInNLMBins[11]){
                legendNLMs->AddEntry((TObject*)0,"X>3" ,"");
                DrawGammaSetMarker(histoRatioPtInNLMBins[11], markerStyleNLM[3], 1, colorNLM[3], colorNLM[3]);  
                histoRatioPtInNLMBins[11]->DrawCopy("same, e1");      
                legendNLMs->AddEntry(histoRatioPtInNLMBins[11],"Data","p");
                if (histoMCRatioPtInNLMBins[11]){
                    DrawGammaSetMarker(histoMCRatioPtInNLMBins[11], markerStyleNLMMC[3], 1, colorNLM[3], colorNLM[3]);  
                    histoMCRatioPtInNLMBins[11]->DrawCopy("e1,same");
                    legendNLMs->AddEntry(histoMCRatioPtInNLMBins[11],"MC","p");
                } else {
                    legendNLMs->AddEntry((TObject*)0,"","");
                }    
            }
            legendNLMs->Draw();
            PutProcessLabelAndEnergyOnPlot(0.65, 0.95, 28, collisionSystem.Data(),fDetectionProcess.Data(), "" , 63, 0.03);
        canvasNLMFractions->Update();
        canvasNLMFractions->SaveAs(Form("%s/%s_%s_NLMFractions_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
    }
    
    // ********************************************************************************************************************************
    // ************************** find example bin in original file and put it into output ********************************************
    // ********************************************************************************************************************************
    Double_t scaleFactorSingleBin           = 1.0;
    Int_t fExampleBin                       = ReturnSingleInvariantMassBinPlotting (nameMeson, optionEnergy, mode, trigger.Atoi(), scaleFactorSingleBin, -1);
    TH1D* histoDataM02ExBin                 = (TH1D*)fileUncorrected.Get(Form("M02_PtBin%d",fExampleBin));    
    if(histoDataM02ExBin){
        histoDataM02ExBin->SetName(Form("Data_M02_in_Pt_Bin%d",fExampleBin));\
        histoDataM02ExBin->Scale(1./histoDataM02ExBin->Integral());
    }
    TH1D* histoMCrecM02ExBin                = (TH1D*)fileCorrections->Get(Form("M02_in_Pt_Bin%d",fExampleBin));
    TH1D* histoTruePi0M02ExBin              = (TH1D*)fileCorrections->Get(Form("TrueClusPi0_M02_in_Pt_Bin%d",fExampleBin));
    TH1D* histoTrueEtaM02ExBin              = (TH1D*)fileCorrections->Get(Form("TrueClusEta_M02_in_Pt_Bin%d",fExampleBin));
    TH1D* histoTrueGammaM02ExBin            = (TH1D*)fileCorrections->Get(Form("TrueClusGamma_M02_in_Pt_Bin%d",fExampleBin));
    TH1D* histoTrueElectronM02ExBin         = (TH1D*)fileCorrections->Get(Form("TrueClusElectron_M02_in_Pt_Bin%d",fExampleBin));
    TH1D* histoTrueBGM02ExBin               = (TH1D*)fileCorrections->Get(Form("TrueClusBG_M02_in_Pt_Bin%d",fExampleBin));
    TH1D* histoTruePi0PureMergedM02ExBin    = (TH1D*)fileCorrections->Get(Form("TrueClusMergedPureFromPi0_M02_in_Pt_Bin%d",fExampleBin));
    TH1D* histoTruePi0PConvMergedM02ExBin   = (TH1D*)fileCorrections->Get(Form("TrueClusPartConvMergedFromPi0_M02_in_Pt_Bin%d",fExampleBin));
    TH1D* histoTruePi0OneGammaM02ExBin      = (TH1D*)fileCorrections->Get(Form("TrueClusOneGammaFromPi0_M02_in_Pt_Bin%d",fExampleBin));
    TH1D* histoTruePi0OneElectronM02ExBin   = (TH1D*)fileCorrections->Get(Form("TrueClusOneElectronFromPi0_M02_in_Pt_Bin%d",fExampleBin));
    if(histoMCrecM02ExBin){
        histoMCrecM02ExBin->SetName(Form("MCrec_M02_in_Pt_Bin%d",fExampleBin));
        Int_t integMC = histoMCrecM02ExBin->Integral();
        histoMCrecM02ExBin->Scale(1./integMC);
        histoTruePi0M02ExBin->Scale(1./integMC);
        histoTrueEtaM02ExBin->Scale(1./integMC);
        histoTrueGammaM02ExBin->Scale(1./integMC);
        histoTrueBGM02ExBin->Scale(1./integMC);
        histoTrueElectronM02ExBin->Scale(1./integMC);
        histoTruePi0PureMergedM02ExBin->Scale(1./integMC);
        histoTruePi0PConvMergedM02ExBin->Scale(1./integMC);
        histoTruePi0OneGammaM02ExBin->Scale(1./integMC);
        histoTruePi0OneElectronM02ExBin->Scale(1./integMC);
    }
    
    if (histoDataM02ExBin && histoMCrecM02ExBin && !kIsMC){
        Double_t textSizeLabelsPixel                 = 50*3/5;
        TCanvas* canvasM02SamplePlot    = new TCanvas("canvasM02SamplePlot","",0,0,750,750);  // gives the page size
        DrawGammaCanvasSettings( canvasM02SamplePlot,  0.08, 0.01, 0.015, 0.085);
        canvasM02SamplePlot->SetLogy(1);
        
        Style_t markerStyleM02MCrec     = 1;
        Size_t markerSizeM02MCrec       = 0;
        Color_t markerColorM02MCrec     = kRed+2;
        Style_t markerStyleM02Data      = 20;
        Size_t markerSizeM02Data        = 1.2;
        Color_t markerColorM02Data      = kBlack;
        Style_t markerStyleM02MC        = 24;
        Size_t markerSizeM02MC          = 1.3;
        Color_t markerColorM02MC        = kRed+2;
        Style_t markerStyleM02Pi0       = 20;
        Size_t markerSizeM02Pi0         = 1.2;
        Color_t markerColorM02Pi0       = kBlack;
        Style_t markerStyleM02Eta       = 24;
        Size_t markerSizeM02Eta         = 1.2;
        Color_t markerColorM02Eta       = kAzure+2;
        Style_t markerStyleM02Gamma     = 24;
        Size_t markerSizeM02Gamma       = 0;
        Color_t markerColorM02Gamma     = 807;
        Style_t markerStyleM02Elec      = 24;
        Size_t markerSizeM02Elec        = 0;
        Color_t markerColorM02Elec      = kGreen-5;
        Style_t markerStyleM02BG        = 24;
        Size_t markerSizeM02BG          = 0;
        Color_t markerColorM02BG        = kBlue+2;
                
        Double_t marginM02          = 0.07*750;
        Double_t textsizeLabelsM02  = 0;
        Double_t textsizeFacM02     = 0;
        Double_t minYAxisM02        = 0.00001;
        if (canvasM02SamplePlot->XtoPixel(canvasM02SamplePlot->GetX2()) < canvasM02SamplePlot->YtoPixel(canvasM02SamplePlot->GetY1())){
            textsizeLabelsM02       = (Double_t)textSizeLabelsPixel/canvasM02SamplePlot->XtoPixel(canvasM02SamplePlot->GetX2()) ;
            textsizeFacM02          = (Double_t)1./canvasM02SamplePlot->XtoPixel(canvasM02SamplePlot->GetX2()) ;
        } else {
            textsizeLabelsM02       = (Double_t)textSizeLabelsPixel/canvasM02SamplePlot->YtoPixel(canvasM02SamplePlot->GetY1());
            textsizeFacM02          = (Double_t)1./canvasM02SamplePlot->YtoPixel(canvasM02SamplePlot->GetY1());
        }
        cout << textsizeLabelsM02 << endl;

        //****************************************************************************
        // plotting everyting together MC
        //****************************************************************************
        TH2F * histo2DPi0M02Dummy;
        histo2DPi0M02Dummy             = new TH2F("histo2DPi0M02Dummy","histo2DPi0M02Dummy",11000,0.0,2.,10000,0,1.2);
        SetStyleHistoTH2ForGraphs(histo2DPi0M02Dummy, "#it{#sigma}_{long}^{2}","#it{P}",0.85*textsizeLabelsM02, textsizeLabelsM02,
                                0.85*textsizeLabelsM02, textsizeLabelsM02,0.88, 0.115/(textsizeFacM02*marginM02));
        
        canvasM02SamplePlot->cd();
        histo2DPi0M02Dummy->GetYaxis()->SetRangeUser(minYAxisM02,1.2);
        histo2DPi0M02Dummy->GetYaxis()->SetTickLength(0.025);
        histo2DPi0M02Dummy->GetXaxis()->SetTickLength(0.025);
        histo2DPi0M02Dummy->DrawCopy();

        TString range = histoMCrecM02ExBin->GetTitle();
//         range.ReplaceAll("#it{p}_{T}","#it{E}_{T}");
        
        TLatex *labelM02PtRange = new TLatex(0.965,0.925,Form("%s", range.Data()));
//         TLatex *labelM02PtRange = new TLatex(0.965,0.93,Form("#pi^{0}: %s", histoMCrecM02ExBin->GetTitle()));
        SetStyleTLatex( labelM02PtRange, 0.85*textSizeLabelsPixel,4);
        labelM02PtRange->SetTextAlign(31);
        labelM02PtRange->SetTextFont(43);

        TLatex *labelM02Energy      = new TLatex(0.11,0.925-1*0.8*textsizeLabelsM02,collisionSystem.Data());
        SetStyleTLatex( labelM02Energy, 0.85*textSizeLabelsPixel,4);
        labelM02Energy->SetTextFont(43);
        
        TLatex *labelM02Trigger      = new TLatex(0.11,0.925-2*0.8*textsizeLabelsM02,Form("%s triggered",nameTrigger.Data()));
        SetStyleTLatex( labelM02Trigger, 0.85*textSizeLabelsPixel,4);
        labelM02Trigger->SetTextFont(43);

        TLatex *labelM02Reco  = new TLatex(0.11,0.925-3*0.8*textsizeLabelsM02,"mEMC");
        SetStyleTLatex( labelM02Reco, 0.85*textSizeLabelsPixel,4);
        labelM02Reco->SetTextFont(43);
        TLatex *labelM02Simulation  = new TLatex(0.11,0.925,"ALICE simulation");
        SetStyleTLatex( labelM02Simulation, 0.85*textSizeLabelsPixel,4);
        labelM02Simulation->SetTextFont(43);
        TLatex *labelM02Performance  = new TLatex(0.11,0.925,"ALICE performance");
        SetStyleTLatex( labelM02Performance, 0.85*textSizeLabelsPixel,4);
        labelM02Performance->SetTextFont(43);
        
        histoMCrecM02ExBin->SetMinimum(minYAxisM02);
        histoTruePi0M02ExBin->SetMinimum(minYAxisM02);
        histoTrueGammaM02ExBin->SetMinimum(minYAxisM02);
        histoTrueEtaM02ExBin->SetMinimum(minYAxisM02);
        histoTrueElectronM02ExBin->SetMinimum(minYAxisM02);
        histoTrueBGM02ExBin->SetMinimum(minYAxisM02);
        DrawGammaSetMarker(histoMCrecM02ExBin, markerStyleM02MCrec, markerSizeM02MCrec, markerColorM02MCrec, markerColorM02MCrec);
        histoMCrecM02ExBin->SetLineWidth(3.5);
        histoMCrecM02ExBin->Draw("hist,e,same");
        DrawGammaSetMarker(histoTrueGammaM02ExBin, markerStyleM02Gamma, markerSizeM02Gamma, markerColorM02Gamma, markerColorM02Gamma);
        histoTrueGammaM02ExBin->SetLineWidth(3);
        histoTrueGammaM02ExBin->SetFillColor(markerColorM02Gamma);
        histoTrueGammaM02ExBin->SetFillStyle(3154);
        histoTrueGammaM02ExBin->Draw("hist,B,same");
        histoTrueGammaM02ExBin->Draw("hist,same");
        DrawGammaSetMarker(histoTrueElectronM02ExBin, markerStyleM02Elec, markerSizeM02Elec, markerColorM02Elec, markerColorM02Elec);
        histoTrueElectronM02ExBin->SetLineWidth(3);
        histoTrueElectronM02ExBin->SetFillColor(markerColorM02Elec);
        histoTrueElectronM02ExBin->SetFillStyle(3145);
        histoTrueElectronM02ExBin->Draw("b,same,hist");
        histoTrueElectronM02ExBin->Draw("same,hist");
        DrawGammaSetMarker(histoTrueBGM02ExBin, markerStyleM02BG, markerSizeM02BG, markerColorM02BG, markerColorM02BG);
        histoTrueBGM02ExBin->SetLineWidth(3.5);
        histoTrueBGM02ExBin->Draw("hist,same");
        DrawGammaSetMarker(histoTruePi0M02ExBin, markerStyleM02Pi0, markerSizeM02Pi0, markerColorM02Pi0, markerColorM02Pi0);
        histoTruePi0M02ExBin->Draw("p,same");
        DrawGammaSetMarker(histoTrueEtaM02ExBin, markerStyleM02Eta, markerSizeM02Eta, markerColorM02Eta, markerColorM02Eta);
        histoTrueEtaM02ExBin->Draw("p,same");
        labelM02PtRange->Draw();
        labelM02Energy->Draw();
        labelM02Trigger->Draw();
        labelM02Reco->Draw();
        labelM02Simulation->Draw();
        
        histo2DPi0M02Dummy->Draw("AXIS,same");
        
        TLegend* legendM02Val  = GetAndSetLegend2(0.75, 0.89-6*0.75*textsizeLabelsM02, 0.9, 0.89, 0.85*textSizeLabelsPixel);
        legendM02Val->SetMargin(0.05/(0.9-0.75));
        legendM02Val->AddEntry(histoMCrecM02ExBin,"Clusters","l");
        legendM02Val->AddEntry(histoTruePi0M02ExBin,"#pi^{0}","p");
        legendM02Val->AddEntry(histoTrueEtaM02ExBin,"#eta BG","p");
        legendM02Val->AddEntry(histoTrueGammaM02ExBin,"#gamma BG","f");
        legendM02Val->AddEntry(histoTrueElectronM02ExBin,"e^{#pm} BG","f");
        legendM02Val->AddEntry(histoTrueBGM02ExBin,"Had. BG","l");
        legendM02Val->Draw();
        canvasM02SamplePlot->SaveAs(Form("%s/Pi0_ValidatedM02BinMEMC_%s.%s",outputDir.Data(), nameTrigger.Data(), suffix.Data()));

        //****************************************************************************
        //******* Set ranges for sub components***************************************
        //****************************************************************************
        histoTruePi0PConvMergedM02ExBin->SetMinimum(minYAxisM02);
        histoTruePi0OneGammaM02ExBin->SetMinimum(minYAxisM02);
        histoTruePi0PureMergedM02ExBin->SetMinimum(minYAxisM02);
        histoTruePi0OneElectronM02ExBin->SetMinimum(minYAxisM02);
        TH1D* histoTruePi0MergedM02ExBin = (TH1D*)histoTruePi0PureMergedM02ExBin->Clone(Form("TrueClusMergedPi0_M02_in_Pt_Bin%d",fExampleBin));
        histoTruePi0MergedM02ExBin->Add(histoTruePi0PConvMergedM02ExBin);
        TH1D* histoTruePi0SinglePartM02ExBin = (TH1D*)histoTruePi0OneGammaM02ExBin->Clone(Form("TrueClusSinglePartPi0_M02_in_Pt_Bin%d",fExampleBin));
        histoTruePi0SinglePartM02ExBin->Add(histoTruePi0OneElectronM02ExBin);
        
        //****************************************************************************
        // plotting everything together MC + decomposed pi0 signal *******************
        //****************************************************************************
        canvasM02SamplePlot->cd();
        histo2DPi0M02Dummy->DrawCopy();
        
        DrawGammaSetMarker(histoMCrecM02ExBin, markerStyleM02MCrec, markerSizeM02MCrec, markerColorM02MCrec, markerColorM02MCrec);
        histoMCrecM02ExBin->SetLineWidth(3.5);
        histoMCrecM02ExBin->Draw("hist,e,same");
        DrawGammaSetMarker(histoTrueGammaM02ExBin, markerStyleM02Gamma, markerSizeM02Gamma, markerColorM02Gamma, markerColorM02Gamma);
        histoTrueGammaM02ExBin->SetLineWidth(3);
        histoTrueGammaM02ExBin->SetFillColor(markerColorM02Gamma);
        histoTrueGammaM02ExBin->SetFillStyle(3154);
        histoTrueGammaM02ExBin->Draw("hist,B,same");
        histoTrueGammaM02ExBin->Draw("hist,same");
        DrawGammaSetMarker(histoTrueElectronM02ExBin, markerStyleM02Elec, markerSizeM02Elec, markerColorM02Elec, markerColorM02Elec);
        histoTrueElectronM02ExBin->SetLineWidth(3);
        histoTrueElectronM02ExBin->SetFillColor(markerColorM02Elec);
        histoTrueElectronM02ExBin->SetFillStyle(3145);
        histoTrueElectronM02ExBin->Draw("b,same,hist");
        histoTrueElectronM02ExBin->Draw("same,hist");
        DrawGammaSetMarker(histoTrueBGM02ExBin, markerStyleM02BG, markerSizeM02BG, markerColorM02BG, markerColorM02BG);
        histoTrueBGM02ExBin->SetLineWidth(3.5);
        histoTrueBGM02ExBin->Draw("hist,same");
        DrawGammaSetMarker(histoTruePi0MergedM02ExBin, markerStyleM02Pi0, markerSizeM02Pi0, markerColorM02Pi0, markerColorM02Pi0);
        histoTruePi0MergedM02ExBin->Draw("p,same");
        DrawGammaSetMarker(histoTruePi0SinglePartM02ExBin, markerStyleM02Pi0+4, markerSizeM02Pi0, markerColorM02Pi0, markerColorM02Pi0);
        histoTruePi0SinglePartM02ExBin->Draw("p,same");

        DrawGammaSetMarker(histoTrueEtaM02ExBin, markerStyleM02Eta, markerSizeM02Eta, markerColorM02Eta, markerColorM02Eta);
        histoTrueEtaM02ExBin->Draw("p,same");
        labelM02PtRange->Draw();
        labelM02Energy->Draw();
        labelM02Trigger->Draw();
        labelM02Reco->Draw();
        labelM02Simulation->Draw();
        
        histo2DPi0M02Dummy->Draw("AXIS,same");
        
        TLegend* legendM02Val2  = GetAndSetLegend2(0.6, 0.89-7*0.75*textsizeLabelsM02, 0.9, 0.89, 0.85*textSizeLabelsPixel);
        legendM02Val2->SetMargin(0.05/(0.9-0.7));
        legendM02Val2->AddEntry(histoMCrecM02ExBin,"Clusters","l");
        legendM02Val2->AddEntry(histoTruePi0MergedM02ExBin,"#pi^{0} (merged showers)","p");
        legendM02Val2->AddEntry(histoTruePi0SinglePartM02ExBin,"#pi^{0} (single showers)","p");
        legendM02Val2->AddEntry(histoTrueEtaM02ExBin,"#eta BG","p");
        legendM02Val2->AddEntry(histoTrueGammaM02ExBin,"#gamma BG","f");
        legendM02Val2->AddEntry(histoTrueElectronM02ExBin,"e^{#pm} BG","f");
        legendM02Val2->AddEntry(histoTrueBGM02ExBin,"Had. BG","l");
        legendM02Val2->Draw();
        cout << "M02 plotting sample bin" << endl;
        canvasM02SamplePlot->SaveAs(Form("%s/Pi0_ValidatedM02BinPlusPi0DecompMEMC_%s.%s",outputDir.Data(), nameTrigger.Data(), suffix.Data()));
        
        
        //****************************************************************************
        // plotting decomposed pi0 signal ********************************************
        //****************************************************************************
        histo2DPi0M02Dummy->DrawCopy();
        
        DrawGammaSetMarker(histoTruePi0MergedM02ExBin, markerStyleM02Gamma, markerSizeM02Gamma, markerColorM02Gamma, markerColorM02Gamma);
        histoTruePi0MergedM02ExBin->SetLineWidth(3);
        histoTruePi0MergedM02ExBin->SetFillColor(markerColorM02Gamma);
        histoTruePi0MergedM02ExBin->SetFillStyle(3154);
        histoTruePi0MergedM02ExBin->Draw("hist,B,same");
        histoTruePi0MergedM02ExBin->Draw("hist,same");
        DrawGammaSetMarker(histoTruePi0SinglePartM02ExBin, markerStyleM02Elec, markerSizeM02Elec, markerColorM02Elec, markerColorM02Elec);
        histoTruePi0SinglePartM02ExBin->SetLineWidth(3);
        histoTruePi0SinglePartM02ExBin->SetFillColor(markerColorM02Elec);
        histoTruePi0SinglePartM02ExBin->SetFillStyle(3145);
        histoTruePi0SinglePartM02ExBin->Draw("b,same,hist");
        histoTruePi0SinglePartM02ExBin->Draw("same,hist");
        DrawGammaSetMarker(histoTruePi0M02ExBin, markerStyleM02Pi0, markerSizeM02Pi0, markerColorM02Pi0, markerColorM02Pi0);
        histoTruePi0M02ExBin->Draw("p,same");

        histo2DPi0M02Dummy->Draw("AXIS,same");
        labelM02PtRange->Draw();
        labelM02Energy->Draw();
        labelM02Trigger->Draw();
        labelM02Reco->Draw();
        labelM02Simulation->Draw();
        
        TLegend* legendM02ValPi0  = GetAndSetLegend2(0.67, 0.89-3*0.75*textsizeLabelsM02, 0.9, 0.89, 0.85*textSizeLabelsPixel);
        legendM02ValPi0->SetMargin(0.05/(0.9-0.67));
        legendM02ValPi0->AddEntry(histoTruePi0M02ExBin,"All #pi^{0}","p");
        legendM02ValPi0->AddEntry(histoTruePi0MergedM02ExBin,"Both #gamma in clus.","f");
        legendM02ValPi0->AddEntry(histoTruePi0SinglePartM02ExBin,"1 #gamma in clus.","f");
        legendM02ValPi0->Draw();
        canvasM02SamplePlot->SaveAs(Form("%s/Pi0_ValidatedM02BinOnlyPi0MEMC_%s.%s",outputDir.Data(), nameTrigger.Data(), suffix.Data()));

        
        //****************************************************************************
        // plotting data vs MC comp **************************************************
        //****************************************************************************
        histo2DPi0M02Dummy->DrawCopy();
        histoDataM02ExBin->SetMinimum(minYAxisM02);
        DrawGammaSetMarker(histoDataM02ExBin, markerStyleM02Data, markerSizeM02Data, markerColorM02Data, markerColorM02Data);
        histoDataM02ExBin->Draw("p,same");
        DrawGammaSetMarker(histoMCrecM02ExBin, markerStyleM02MC, markerSizeM02MC, markerColorM02MC, markerColorM02MC);
        histoMCrecM02ExBin->Draw("p,same");

        histo2DPi0M02Dummy->Draw("AXIS,same");
        labelM02PtRange->Draw();
        labelM02Energy->Draw();
        labelM02Trigger->Draw();
        labelM02Reco->Draw();
        labelM02Performance->Draw();
        
        TLegend* legendM02Data  = GetAndSetLegend2(0.8, 0.89-2*0.75*textsizeLabelsM02, 0.9, 0.89, 0.85*textSizeLabelsPixel);
        legendM02Data->SetMargin(0.05/(0.9-0.8));
        legendM02Data->AddEntry(histoDataM02ExBin,"Data","p");
        legendM02Data->AddEntry(histoMCrecM02ExBin,"MC","p");
        legendM02Data->Draw();
        canvasM02SamplePlot->SaveAs(Form("%s/Pi0_DataVsMCM02BinOnlyPi0MEMC_%s.%s",outputDir.Data(), nameTrigger.Data(), suffix.Data()));
        
    }
    
    // ********************************************************************************************************************************
    // ****************************** Write file with all further needed histograms ***************************************************
    // ********************************************************************************************************************************    
    cout << "writing corrected file" << endl;    
    
    const char* nameOutput = Form("%s/%s/%s_%s_GammaMergedCorrection%s_%s.root",fCutSelection.Data(),optionEnergy.Data(),nameMeson.Data(),prefix2.Data(),optionPeriod.Data(),fCutSelection.Data());
    TFile* correctedOutput = new TFile(nameOutput,"RECREATE");  

    if (histoCorrectedYieldTrue)            histoCorrectedYieldTrue->Write("CorrectedYieldTrueEff");
    
    for (Int_t j = 0; j<4; j++){
        if (histoYieldSecMeson[j])          histoYieldSecMeson[j]->Write();
        if (histoRatioYieldSecMeson[j])     histoRatioYieldSecMeson[j]->Write();
        if (histoSecEffiPt[j])              histoSecEffiPt[j]->Write();
        if (histoSecAccPt[j])               histoSecAccPt[j]->Write();
        if (j < 3){
            if (histoYieldSecMesonFromToy[j])        histoYieldSecMesonFromToy[j]->Write();
            if (histoRatioYieldSecMesonFromToy[j])   histoRatioYieldSecMesonFromToy[j]->Write();
        }
    }
    
    if (histoUnCorrectedYield)              histoUnCorrectedYield->Write();

    if (histoAcceptance)                    histoAcceptance->Write();
    if (histoTrueEffiPrimMesonPt)           histoTrueEffiPrimMesonPt->Write("PrimaryMesonEfficiency");
    if (histoMesonPurityPt)                 histoMesonPurityPt->Write("MesonPurity");
    if (histoEventQuality)                  histoEventQuality->Write();
    if (histoInputMesonPt)                  histoInputMesonPt->Write();
    if (histoMCYieldMeson)                  histoMCYieldMeson->Write();
    if (histoMCYieldMesonOldBin)            histoMCYieldMesonOldBin->Write();
    if (histoUnCorrectedYieldDrawing){
        histoUnCorrectedYieldDrawing->SetName("histoYieldMesonPerEvent");
        histoUnCorrectedYieldDrawing->Write();
    }
    
    if (deltaPt)                            deltaPt->Write("deltaPt");
    if (kIsMC){
        if (histoRatioMergedPureFrac)           histoRatioMergedPureFrac->Write();
        if (histoRatioMergedPartConvFrac)       histoRatioMergedPartConvFrac->Write();
        if (histoRatioMergedOneGammaFrac)       histoRatioMergedOneGammaFrac->Write();
        if (histoRatioMergedOneElectronFrac)    histoRatioMergedOneElectronFrac->Write();
    }
    
    //write example bins
    if (histoDataM02ExBin)                      histoDataM02ExBin->Write();
    if (histoMCrecM02ExBin)                     histoMCrecM02ExBin->Write();
    if (histoTruePi0M02ExBin)                   histoTruePi0M02ExBin->Write();
    if (histoTrueEtaM02ExBin)                   histoTrueEtaM02ExBin->Write();
    if (histoTrueGammaM02ExBin)                 histoTrueGammaM02ExBin->Write();
    if (histoTrueElectronM02ExBin)              histoTrueElectronM02ExBin->Write();
    if (histoTrueBGM02ExBin)                    histoTrueBGM02ExBin->Write();
    if (histoTruePi0PureMergedM02ExBin)         histoTruePi0PureMergedM02ExBin->Write();
    if (histoTruePi0PConvMergedM02ExBin)        histoTruePi0PConvMergedM02ExBin->Write();
    if (histoTruePi0OneGammaM02ExBin)           histoTruePi0OneGammaM02ExBin->Write();
    if (histoTruePi0OneElectronM02ExBin)        histoTruePi0OneElectronM02ExBin->Write();
    
    correctedOutput->Write();
    correctedOutput->Close();
    
}
