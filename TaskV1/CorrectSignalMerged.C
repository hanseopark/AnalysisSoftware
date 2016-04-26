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
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

struct SysErrorConversion {
    Double_t value;
    Double_t error;
};

void CorrectYield( TH1D* histoCorrectedYield, 
                   TH1D* histoPurity,
                   TH1D* histoRawSecYield, 
                   TH1D* histoRawAddSecYieldFromK0s, 
                   TH1D* histoEffiPt, 
                   TH1D* histoAcceptance, 
                   Double_t deltaRapid, 
                   Double_t scaling, 
                   Double_t nEvt, 
                   TString nameMeson
                 ){
    histoCorrectedYield->Sumw2();
    histoCorrectedYield->Multiply(histoPurity);
    histoCorrectedYield->Add(histoRawSecYield,-1.);
    histoCorrectedYield->Add(histoRawAddSecYieldFromK0s,-1.);
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
        histoCorrectedYield->Scale(1./0.98798);
    }else{
        histoCorrectedYield->Scale(1./0.3931);
    }
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
    if (nameMeson.CompareTo("Pi0") == 0 ||nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
        histoCorrectedToBeScaled->Scale(1./0.98798);
    }else{
        histoCorrectedToBeScaled->Scale(1./0.3931);
        
    }
}


void  CorrectSignalMerged(  TString fileNameUnCorrectedFile = "myOutput", 
                            TString fileNameCorrectionFile  = "", 
                            TString fCutSelection           = "", 
                            TString suffix                  = "gif", 
                            TString nameMeson               = "", 
                            Bool_t  kIsMC                   = kFALSE, 
                            TString optionEnergy            = "", 
                            TString optionPeriod            = "",
                            Int_t mode                      = 10
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
    
    //Variable defintion
    Double_t scaling                            = 1./(2.*TMath::Pi());
    
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
    TH1D *histoPurityPt                         = (TH1D*)fileCorrections->Get("TruePurityMergedPt");
    TString namePurity                          = "TruePurityPi0Pt";
    if (kIsEta) namePurity                      = "TruePurityEtaPt";
    TH1D* histoMesonPurityPt                    = (TH1D*)fileCorrections->Get(namePurity.Data());
    TH1D *histoAcceptance                       = (TH1D*)fileCorrections->Get("AcceptancePt");
    
    TString labelsBG[10]                        = { "PiCh", "P", "KCh", "N", "K0s", 
                                                    "Lambda", "Mu", "K0l", "Rest", "AllNonEM"};
    TString labelsBGPlot[10]                    = { "#pi^{#pm}", "p/#bar{p}", "K^{#pm}", "n", "K^{0}_{s}", 
                                                    "#Lambda", "#mu^{#pm}", "K^{0}_{L}", "rest", "all BG non em."};
    Color_t colorBGPlot[10]                     = { kRed+1, kBlue+1, kCyan-2, kTeal+4, 800,
                                                    kOrange+7, kGreen+2, kViolet+2, kAzure+5, kGray+2 };
    Style_t markerStyleBGPlot[10]               = { 20, 21, 24, 25, 33, 
                                                    27, 34, 24, 28, 20 };
    TH1D* histoTrueClustersBGPt[10]             = { NULL, NULL, NULL, NULL, NULL,
                                                    NULL, NULL, NULL, NULL, NULL};                                        
    TH1D* histoRatioTrueClustersBGPt[10]        = { NULL, NULL, NULL, NULL, NULL,
                                                    NULL, NULL, NULL, NULL, NULL};                                        
    TH1D* histoTrueYieldEtaFullM02              = NULL;
    TH1D* histoTrueYieldPi0FullM02              = NULL;
    TH1D* histoTrueYieldGammaM02                = NULL;
    TH1D* histoTrueYieldElectronM02             = NULL;
    TH1D* histoRatioTrueYieldEtaFullM02         = NULL;
    TH1D* histoRatioTrueYieldPi0FullM02         = NULL;
    TH1D* histoRatioTrueYieldGammaM02           = NULL;
    TH1D* histoRatioTrueYieldElectronM02        = NULL;

    if (kIsMC){
        for (Int_t i = 0; i < 9; i++){
            histoTrueClustersBGPt[i]                = (TH1D*)fileCorrections->Get(Form("TrueClusBG_%s_Pt",labelsBG[i].Data()));
            histoTrueClustersBGPt[i]->Sumw2();
            
            histoRatioTrueClustersBGPt[i]           = (TH1D*)histoTrueClustersBGPt[i]->Clone(Form("RatioTrueClusBG_%s_Pt",labelsBG[i].Data()));
            histoRatioTrueClustersBGPt[i]->Divide(histoRatioTrueClustersBGPt[i],histoUnCorrectedYield,1.,1.,"");
            if (i == 0) histoTrueClustersBGPt[9]    = (TH1D*)histoTrueClustersBGPt[i]->Clone(Form("TrueClusBG_%s_Pt",labelsBG[9].Data()));
                else histoTrueClustersBGPt[9]->Add(histoTrueClustersBGPt[i]); 
        }
        histoRatioTrueClustersBGPt[9]               = (TH1D*)histoTrueClustersBGPt[9]->Clone(Form("RatioTrueClusBG_%s_Pt",labelsBG[9].Data()));
        histoRatioTrueClustersBGPt[9]->Divide(histoRatioTrueClustersBGPt[9],histoUnCorrectedYield,1.,1.,"");
        
        histoTrueYieldEtaFullM02                    = (TH1D*)fileCorrections->Get("histoTrueYieldEtaFullM02");
        histoRatioTrueYieldEtaFullM02               = (TH1D*)histoTrueYieldEtaFullM02->Clone("RatioTrueYieldEtaFullM02");
        histoRatioTrueYieldEtaFullM02->Divide(histoRatioTrueYieldEtaFullM02,histoUnCorrectedYield,1.,1.,"");
        
        histoTrueYieldPi0FullM02                    = (TH1D*)fileCorrections->Get("histoTrueYieldPi0FullM02");
        histoRatioTrueYieldPi0FullM02               = (TH1D*)histoTrueYieldPi0FullM02->Clone("RatioTrueYieldPi0FullM02");
        histoRatioTrueYieldPi0FullM02->Divide(histoRatioTrueYieldPi0FullM02,histoUnCorrectedYield,1.,1.,"");
        
        histoTrueYieldGammaM02                      = (TH1D*)fileCorrections->Get("histoTrueYieldGammaM02");
        histoRatioTrueYieldGammaM02                 = (TH1D*)histoTrueYieldGammaM02->Clone("RatioTrueYieldGammaM02");
        histoRatioTrueYieldGammaM02->Divide(histoRatioTrueYieldGammaM02,histoUnCorrectedYield,1.,1.,"");
        
        histoTrueYieldElectronM02                   = (TH1D*)fileCorrections->Get("histoTrueYieldElectronM02");
        histoRatioTrueYieldElectronM02              = (TH1D*)histoTrueYieldElectronM02->Clone("RatioTrueYieldElectronM02");
        histoRatioTrueYieldElectronM02->Divide(histoRatioTrueYieldElectronM02,histoUnCorrectedYield,1.,1.,"");
    }
    
    // loading efficiency
    TH1D *histoTrueEffiPt                       = NULL;
    TH1D *histoTrueEffiPrimMesonPt              = NULL;
    histoTrueEffiPt                             = (TH1D*)fileCorrections->Get("TrueEfficiencyMergedPt"); 
    histoTrueEffiPrimMesonPt                    = (TH1D*)fileCorrections->Get("TrueEfficiencyPrimMesonPt"); 
    TH1D* histoRatioEffi                        = (TH1D*)histoTrueEffiPrimMesonPt->Clone("histoRatioEffi");
    histoRatioEffi->Divide(histoRatioEffi,histoTrueEffiPt,1.,1.,"B");

    
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
  
    TH1D *histoYieldTrueSecFracMeson            = NULL;
    TH1D *histoYieldTrueSecFracFK0SMeson        = NULL;
    TH1D *histoYieldTrueSecFracMeson_Or         = NULL;
    TH1D *histoYieldTrueSecFracFK0SMeson_Or     = NULL;
    TH1D *histoYieldTrueSecFracFLambdaMeson     = NULL;
    TH1D *histoYieldTrueSecFracFLambdaMeson_Or  = NULL;
    TF1* fitSecFrac                             = new TF1("fitSecFrac","[0]/pow(x,[1])+[2]+[3]*x");
    fitSecFrac->SetParLimits(2,0,10);
    fitSecFrac->SetParLimits(3,0,1e-2);
    TF1* fitSecFracFromK0                       = new TF1("fitSecFrac","[0]/pow(x,[1])+[2]+[3]*x");
    fitSecFracFromK0->SetParLimits(2,0,10);
    fitSecFracFromK0->SetParLimits(3,0,1e-2);
    TH1D* histoYieldSecMeson                    = NULL;
    TH1D* histoYieldSecFromK0SMeson             = NULL;
    

    Bool_t doK0SecCorrection                    = kFALSE;
    Int_t doK0SecCorrectionWithDefaultHisto     = 0;
    
    if ( !kIsEta ) {
        doK0SecCorrection                       = kTRUE;
    }
    
    if (doK0SecCorrection){
        histoYieldTrueSecFracMeson              = (TH1D*)fileCorrections->Get("TrueSecFrac");
        histoYieldTrueSecFracMeson_Or           = (TH1D*)histoYieldTrueSecFracMeson->Clone("TrueSecFrac_Or");
        histoYieldTrueSecFracFK0SMeson          = (TH1D*)fileCorrections->Get("TrueSecFracFromK0S");
        histoYieldTrueSecFracFLambdaMeson       = (TH1D*)fileCorrections->Get("TrueSecFracFromLambda");
        histoYieldTrueSecFracFK0SMeson_Or       = (TH1D*)histoYieldTrueSecFracFK0SMeson->Clone("TrueSecFracFromK0S_Or");
        histoYieldTrueSecFracFLambdaMeson_Or    = (TH1D*)histoYieldTrueSecFracFLambdaMeson->Clone("TrueSecFracFromLambda_Or");
        Double_t maxPtSecondaries               = histoYieldTrueSecFracMeson->GetXaxis()->GetBinUpEdge(histoYieldTrueSecFracMeson->GetNbinsX());
        
        cout << "Fit standard range - secondaries"<< endl;  
        fitSecFrac->SetRange(minPtMesonSec, maxPtSecondaries);
        TFitResultPtr resultSecFrac             = histoYieldTrueSecFracMeson->Fit(fitSecFrac, "SINRME+", "", minPtMesonSec, maxPtSecondaries);
                    
        cout << "Fit standard range - from K0s"<< endl;  
        fitSecFracFromK0->SetRange(minPtMesonSec, maxPtSecondaries);
        TFitResultPtr resultSecFracFromK0       = histoYieldTrueSecFracFK0SMeson->Fit(fitSecFracFromK0, "SINRME+", "", minPtMesonSec, maxPtSecondaries);

        for (Int_t i = 2; i < histoYieldTrueSecFracMeson->GetNbinsX()+1; i++){
            Double_t ptStart                    = histoYieldTrueSecFracMeson->GetXaxis()->GetBinLowEdge(i);
            Double_t ptEnd                      = histoYieldTrueSecFracMeson->GetXaxis()->GetBinUpEdge(i);
            Double_t binWidth                   = ptEnd-ptStart;
            Double_t secFrac                    = fitSecFrac->Integral(ptStart, ptEnd, resultSecFrac->GetParams()) / binWidth;
            Double_t errorSecFrac               = fitSecFrac->IntegralError(ptStart, ptEnd, resultSecFrac->GetParams(), resultSecFrac->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
            histoYieldTrueSecFracMeson->SetBinContent(i, secFrac);
            histoYieldTrueSecFracMeson->SetBinError(i, errorSecFrac);
                                    
            secFrac                             = fitSecFracFromK0->Integral(ptStart, ptEnd, resultSecFracFromK0->GetParams()) / binWidth;
            errorSecFrac                        = fitSecFracFromK0->IntegralError(ptStart, ptEnd, resultSecFracFromK0->GetParams(), resultSecFracFromK0->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
            histoYieldTrueSecFracFK0SMeson->SetBinContent(i, secFrac);
            histoYieldTrueSecFracFK0SMeson->SetBinError(i, errorSecFrac);
        }
        
        histoYieldSecMeson                      = (TH1D*)histoUnCorrectedYield->Clone("SecFracMeson");
        histoYieldSecMeson->Sumw2();
        histoYieldSecMeson->Multiply(histoMesonPurityPt);
        histoYieldSecMeson->Multiply(histoYieldTrueSecFracMeson);

        histoYieldSecFromK0SMeson               = (TH1D*)histoUnCorrectedYield->Clone("SecFracFromK0SMeson");
        histoYieldSecFromK0SMeson->Sumw2();
        histoYieldSecFromK0SMeson->Multiply(histoMesonPurityPt);
        histoYieldSecFromK0SMeson->Multiply(histoYieldTrueSecFracFK0SMeson);
        histoYieldSecFromK0SMeson->Scale(doubleAddFactorK0s);                
        
        cout << "I am changing the sec yield" << endl;
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
        delete canvasAcceptance;
        
        //**********************************************************************************
        //******************** Secondary Fraction     **************************************
        //**********************************************************************************        
        if ( nameMeson.Contains("Pi0")){
            cout << "Plotting secondary fractions" << endl;
            TCanvas* canvasSecFrac = new TCanvas("canvasSecFrac","",200,10,1350,900);  // gives the page size
            DrawGammaCanvasSettings( canvasSecFrac, 0.083, 0.02, 0.04, 0.09);
            canvasSecFrac->SetLogy(1);
            
            DrawAutoGammaMesonHistos( histoYieldTrueSecFracMeson_Or, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#it{r}_{X} = #frac{X->#pi^{0}}{#pi^{0}}", // (%)", 
                                    kTRUE, 1000, 1e-5, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
            histoYieldTrueSecFracMeson_Or->GetYaxis()->SetTitleOffset(0.75);
            DrawGammaSetMarker(histoYieldTrueSecFracMeson_Or, 20, 1., kBlack, kBlack);  
            DrawGammaSetMarker(histoYieldTrueSecFracFK0SMeson_Or, 24, 1., kBlue, kBlue);
            DrawGammaSetMarker(histoYieldTrueSecFracFLambdaMeson_Or, 24, 1., kViolet+2, kViolet+2);
            histoYieldTrueSecFracMeson_Or->DrawCopy("e1");  
            histoYieldTrueSecFracFK0SMeson_Or->DrawCopy("e1,same");  
            histoYieldTrueSecFracFLambdaMeson_Or->DrawCopy("e1,same"); 
            fitSecFrac->SetLineColor(kBlack);  
            fitSecFracFromK0->SetLineColor(kBlue); 
            fitSecFrac->Draw("same");
            fitSecFracFromK0->Draw("same");
            
            TLegend* legendSecFrac = GetAndSetLegend2(0.7, 0.79, 0.97, 0.93, 32,2); 
            legendSecFrac->AddEntry(histoYieldTrueSecFracMeson_Or,"#it{r}_{All}","pe");
            legendSecFrac->AddEntry(fitSecFrac,"fit to #it{r}_{All}","l");
            legendSecFrac->AddEntry(histoYieldTrueSecFracFK0SMeson_Or,"#it{r}_{K_{s}^{0}}","pe");
            legendSecFrac->AddEntry(fitSecFracFromK0,"fit to #it{r}_{K_{s}^{0}}","l");
            legendSecFrac->AddEntry(histoYieldTrueSecFracFLambdaMeson_Or,"#it{r}_{#Lambda}", "pe"); 
            legendSecFrac->Draw();

            PutProcessLabelAndEnergyOnPlot(0.62, 0.75, 0.03, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data());
            canvasSecFrac->Update();
            canvasSecFrac->SaveAs(Form("%s/%s_FracSecondaries_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
            
            canvasSecFrac->SetLogy(0);

            DrawAutoGammaMesonHistos( histoYieldTrueSecFracMeson_Or, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#it{r}_{X} = #frac{X->#pi^{0}}{#pi^{0}}", // (%)", 
                                    kFALSE, 1e-2, 0, kFALSE,
                                    kTRUE, 0., 0.3, 
                                    kFALSE, 0., 10.);
            histoYieldTrueSecFracMeson_Or->GetYaxis()->SetTitleOffset(0.75);
            DrawGammaSetMarker(histoYieldTrueSecFracMeson_Or, 20, 1., kBlack, kBlack);  
            histoYieldTrueSecFracMeson_Or->DrawCopy("e1");  
            histoYieldTrueSecFracFK0SMeson_Or->DrawCopy("e1,same");  
            histoYieldTrueSecFracFLambdaMeson_Or->DrawCopy("e1,same"); 
            fitSecFrac->Draw("same");
            fitSecFracFromK0->Draw("same");
            legendSecFrac->Draw();
            PutProcessLabelAndEnergyOnPlot(0.62, 0.75, 0.03, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data());
            
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
        
        DrawAutoGammaMesonHistos(   histoTrueEffiPt, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#epsilon_{eff}", 
                                    kTRUE, 2., 3e-5, kFALSE,
                                    kFALSE, 0., 0.3, 
                                    kFALSE, 0., 10.);
        histoTrueEffiPt->GetYaxis()->SetTitleOffset(0.8);        
        DrawGammaSetMarker(histoTrueEffiPt, 20, 1., kBlack, kBlack);
        histoTrueEffiPt->DrawCopy("e1");
        DrawGammaSetMarker(histoTrueEffiPrimMesonPt, 24, 1., kBlue+1, kBlue+1);
        histoTrueEffiPrimMesonPt->DrawCopy("e1,same");
        
        TLegend* legendEff = GetAndSetLegend2(0.4,0.125,0.6,0.205, 28);
        legendEff->SetMargin(0.2);
        legendEff->AddEntry(histoTrueEffiPt,"merged cluster efficiency");            
        legendEff->AddEntry(histoTrueEffiPrimMesonPt,Form("prim. %s efficiency",textMeson.Data()));            
        legendEff->Draw();
        
        PutProcessLabelAndEnergyOnPlot(0.7, 0.25, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasEffSimple->Update();
        canvasEffSimple->SaveAs(Form("%s/%s_TrueEffSimple_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data())); 

        // plotting efficiency linearly
        canvasEffSimple->SetLogy(0);

        DrawAutoGammaMesonHistos(   histoTrueEffiPt, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#epsilon_{eff}", 
                                    kTRUE, 0.55, 0, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
        histoTrueEffiPt->GetYaxis()->SetTitleOffset(0.8);        
        DrawGammaSetMarker(histoTrueEffiPt, 20, 1., kBlack, kBlack);
        
        histoTrueEffiPt->DrawCopy("e1");
        DrawGammaSetMarker(histoTrueEffiPrimMesonPt, 24, 1., kBlue+1, kBlue+1);
        histoTrueEffiPrimMesonPt->DrawCopy("e1,same");

        legendEff->Draw();
        PutProcessLabelAndEnergyOnPlot(0.7, 0.25, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasEffSimple->Update();
        canvasEffSimple->SaveAs(Form("%s/%s_TrueEffSimpleLinY_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data())); 
        delete canvasEffSimple;

        //*********************************************************************************
        //************************** Efficiency Ratio merged vs pi0 merged ********************
        //*********************************************************************************                
        TCanvas* canvasEfficiencyRatio = new TCanvas("canvasEfficiencyRatio","",200,10,1350,900);
        DrawGammaCanvasSettings( canvasEfficiencyRatio, 0.06, 0.015, 0.015, 0.08);

        DrawAutoGammaMesonHistos( histoRatioEffi, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#epsilon_{eff,A}/#epsilon_{eff,B}", 
                                    kFALSE, 0.75, 3e-6, kFALSE,
                                    kTRUE, 0.8, 1.05, 
                                    kFALSE, 0., 10.);
        histoRatioEffi->GetYaxis()->SetTitleOffset(0.7);        
        DrawGammaSetMarker(histoRatioEffi, 20, 1., kAzure+2, kAzure+2); 
        histoRatioEffi->DrawCopy("e1"); 
        DrawGammaLines(0., maxPtMeson,1., 1.,2,kGray+2,7);
        histoRatioEffi->DrawCopy("same,e1"); 
        
        TLegend* legendEfficiencyRatio = GetAndSetLegend2(0.55, 0.125, 0.85, 0.205, 28);
        legendEfficiencyRatio->SetMargin(0.2);
        legendEfficiencyRatio->AddEntry(histoRatioEffi,Form("prim %s Efficiency/Merged Efficiency",textMeson.Data()));
        legendEfficiencyRatio->Draw();   
        
        PutProcessLabelAndEnergyOnPlot(0.125, 0.95, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasEfficiencyRatio->Update();

        canvasEfficiencyRatio->SaveAs(Form("%s/%s_%s_TrueEfficiencyRatio_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        delete legendEfficiencyRatio;
        delete canvasEfficiencyRatio;
        
        
        //*********************************************************************************
        //************************** Purity Plot ******************************************
        //*********************************************************************************
        cout << "Plotting purity" << endl;
        TCanvas* canvasPurity = new TCanvas("canvasPurity","",200,10,1350,900);
        DrawGammaCanvasSettings( canvasPurity, 0.06, 0.015, 0.015, 0.08);

        DrawAutoGammaMesonHistos( histoPurityPt, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#epsilon_{pur}", 
                                    kFALSE, 0.75, 3e-6, kFALSE,
                                    kTRUE, 0., 1, 
                                    kFALSE, 0., 10.);
        histoPurityPt->GetYaxis()->SetTitleOffset(0.7);        
        DrawGammaSetMarker(histoPurityPt, 20, 1., kAzure+2, kAzure+2); 
        histoPurityPt->DrawCopy("e1"); 
        DrawGammaSetMarker(histoMesonPurityPt, 24, 1., kRed+2, kRed+2); 
        histoMesonPurityPt->DrawCopy("e1,same"); 

        TLegend* legendPurity = GetAndSetLegend2(0.5, 0.125, 0.65, 0.205, 28);
        legendPurity->SetMargin(0.2);
        legendPurity->AddEntry(histoPurityPt,"Merged Purity");
        legendPurity->AddEntry(histoMesonPurityPt,Form("%s Purity",textMeson.Data()));
        legendPurity->Draw();   
        
        PutProcessLabelAndEnergyOnPlot(0.68, 0.25, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasPurity->Update();

        canvasPurity->SaveAs(Form("%s/%s_%s_TruePurity_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        delete legendPurity;
        delete canvasPurity;

        
        //*********************************************************************************
        //************************** Purity Ratio merged vs pi0 merged ********************
        //*********************************************************************************        
        TH1D* histoRatioPurities = (TH1D*)histoMesonPurityPt->Clone("histoRatioPurities");
        histoRatioPurities->Divide(histoRatioPurities,histoPurityPt,1.,1.,"B");
        
        TCanvas* canvasPurityRatio = new TCanvas("canvasPurityRatio","",200,10,1350,900);
        DrawGammaCanvasSettings( canvasPurityRatio, 0.06, 0.015, 0.015, 0.08);

        DrawAutoGammaMesonHistos( histoRatioPurities, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#epsilon_{pur,A}/#epsilon_{pur,B}", 
                                    kFALSE, 0.75, 3e-6, kFALSE,
                                    kTRUE, 0.8, 1.05, 
                                    kFALSE, 0., 10.);
        histoRatioPurities->GetYaxis()->SetTitleOffset(0.7);        
        DrawGammaSetMarker(histoRatioPurities, 20, 1., kAzure+2, kAzure+2); 
        histoRatioPurities->DrawCopy("e1"); 
        DrawGammaLines(0., maxPtMeson,1., 1.,2,kGray+2,7);
        histoRatioPurities->DrawCopy("same,e1"); 
        
        TLegend* legendPurityRatio = GetAndSetLegend2(0.65, 0.125, 0.85, 0.205, 28);
        legendPurityRatio->SetMargin(0.2);
        legendPurityRatio->AddEntry(histoRatioPurities,Form("%s Purity/Merged Purity",textMeson.Data()));
        legendPurityRatio->Draw();   
        
        PutProcessLabelAndEnergyOnPlot(0.125, 0.95, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasPurityRatio->Update();

        canvasPurityRatio->SaveAs(Form("%s/%s_%s_TruePurityRatio_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        delete legendPurityRatio;
        delete canvasPurityRatio;
        
        
    }
     
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

        TLegend* legendBG = GetAndSetLegend2(0.55, 0.8, 0.95, 0.95, 28);
        legendBG->SetMargin(0.2);
        legendBG->SetNColumns(3);
        legendBG->AddEntry(histoUnCorrectedYield,"Signal+BG","p");
        for (Int_t i = 0; i < 9; i++){
           DrawGammaSetMarker(histoTrueClustersBGPt[i], markerStyleBGPlot[i], 1, colorBGPlot[i], colorBGPlot[i]);  
           histoTrueClustersBGPt[i]->DrawCopy("e1,same");
           legendBG->AddEntry(histoTrueClustersBGPt[i],labelsBGPlot[i], "p");
        }    
        if (kIsEta){
            DrawGammaSetMarker(histoTrueYieldPi0FullM02, 20, 1, kMagenta+2, kMagenta+2);  
            histoTrueYieldPi0FullM02->DrawCopy("e1,same");
            legendBG->AddEntry(histoTrueYieldPi0FullM02,"#pi^{0}","p");
        } else {
            DrawGammaSetMarker(histoTrueYieldEtaFullM02, 20, 1, kMagenta+2, kMagenta+2);  
            histoTrueYieldEtaFullM02->DrawCopy("e1,same");            
            legendBG->AddEntry(histoTrueYieldEtaFullM02,"#eta","p");
        }
        DrawGammaSetMarker(histoTrueYieldGammaM02, 21, 1, kAzure, kAzure);  
        histoTrueYieldGammaM02->DrawCopy("e1,same");            
        legendBG->AddEntry(histoTrueYieldGammaM02,"#gamma","p");
        DrawGammaSetMarker(histoTrueYieldElectronM02, 21, 1, kGreen+3, kGreen+3);  
        histoTrueYieldElectronM02->DrawCopy("e1,same");            
        legendBG->AddEntry(histoTrueYieldElectronM02,"e^{#pm}","p");

        legendBG->Draw();   
        
        PutProcessLabelAndEnergyOnPlot(0.12, 0.25, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        canvasBGYield->Update();

        canvasBGYield->SaveAs(Form("%s/%s_%s_BGYieldPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        
        DrawGammaSetMarker(histoUnCorrectedYield, 20, 1, kBlack, kBlack);  
        histoUnCorrectedYield->DrawCopy("e1");

        TLegend* legendBG2 = GetAndSetLegend2(0.65, 0.85, 0.95, 0.95, 28);
        legendBG2->SetMargin(0.2);
        legendBG2->SetNColumns(2);
        legendBG2->AddEntry(histoUnCorrectedYield,"Signal+BG","p");
        DrawGammaSetMarker(histoTrueClustersBGPt[9], markerStyleBGPlot[9], 1, colorBGPlot[9], colorBGPlot[9]);  
        histoTrueClustersBGPt[9]->DrawCopy("e1,same");
        legendBG2->AddEntry(histoTrueClustersBGPt[9],labelsBGPlot[9], "p");
        if (kIsEta){
            DrawGammaSetMarker(histoTrueYieldPi0FullM02, 20, 1, kMagenta+2, kMagenta+2);  
            histoTrueYieldPi0FullM02->DrawCopy("e1,same");
            legendBG2->AddEntry(histoTrueYieldPi0FullM02,"#pi^{0}","p");
        } else {
            DrawGammaSetMarker(histoTrueYieldEtaFullM02, 20, 1, kMagenta+2, kMagenta+2);  
            histoTrueYieldEtaFullM02->DrawCopy("e1,same");            
            legendBG2->AddEntry(histoTrueYieldEtaFullM02,"#eta","p");
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
                                        "", "#it{p}_{T} (GeV/#it{c})", "BG Source/ S+B", 
                                        kFALSE, 4., 4e-10, kTRUE,
                                        kTRUE, 1e-5, 10, 
                                        kFALSE, 0., 10.);
            DrawGammaSetMarker(histoRatioTrueYieldGammaM02, 21, 1, kAzure, kAzure);  
            histoRatioTrueYieldGammaM02->DrawCopy("e1");

            TLegend* legendBGRatio = GetAndSetLegend2(0.55, 0.8, 0.95, 0.95, 28);
            legendBGRatio->SetMargin(0.2);
            legendBGRatio->SetNColumns(3);
            legendBGRatio->AddEntry(histoRatioTrueYieldGammaM02,"#gamma","p");

            if (kIsEta){
                DrawGammaSetMarker(histoRatioTrueYieldPi0FullM02, 20, 1, kMagenta+2, kMagenta+2);  
                histoRatioTrueYieldPi0FullM02->DrawCopy("e1,same");
                legendBGRatio->AddEntry(histoRatioTrueYieldPi0FullM02,"#pi^{0}","p");
            } else {
                DrawGammaSetMarker(histoRatioTrueYieldEtaFullM02, 20, 1, kMagenta+2, kMagenta+2);  
                histoRatioTrueYieldEtaFullM02->DrawCopy("e1,same");            
                legendBGRatio->AddEntry(histoRatioTrueYieldEtaFullM02,"#eta","p");
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
                                        "", "#it{p}_{T} (GeV/#it{c})", "BG Source/ S+B", 
                                        kFALSE, 4., 4e-10, kTRUE,
                                        kTRUE, 1e-5, 10, 
                                        kFALSE, 0., 10.);
            DrawGammaSetMarker(histoRatioTrueClustersBGPt[0], markerStyleBGPlot[0], 1, colorBGPlot[0], colorBGPlot[0]);  
            histoRatioTrueClustersBGPt[0]->DrawCopy("e1");

            TLegend* legendBGRatio3 = GetAndSetLegend2(0.55, 0.8, 0.95, 0.95, 28);
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
                                        "", "#it{p}_{T} (GeV/#it{c})", "BG Source/ S+B", 
                                        kFALSE, 4., 4e-10, kTRUE,
                                        kTRUE, 1e-5, 10, 
                                        kFALSE, 0., 10.);
            DrawGammaSetMarker(histoRatioTrueYieldGammaM02, 24, 1, kAzure, kAzure);  
            histoRatioTrueYieldGammaM02->DrawCopy("e1");

            TLegend* legendBGRatio2 = GetAndSetLegend2(0.65, 0.85, 0.95, 0.95, 28);
            legendBGRatio2->SetMargin(0.2);
            legendBGRatio2->SetNColumns(2);
            legendBGRatio2->AddEntry(histoRatioTrueYieldGammaM02,"#gamma","p");

            if (kIsEta){
                DrawGammaSetMarker(histoRatioTrueYieldPi0FullM02, 20, 1, kMagenta+2, kMagenta+2);  
                histoRatioTrueYieldPi0FullM02->DrawCopy("e1,same");
                legendBGRatio2->AddEntry(histoRatioTrueYieldPi0FullM02,"#pi^{0}","p");
            } else {
                DrawGammaSetMarker(histoRatioTrueYieldEtaFullM02, 20, 1, kMagenta+2, kMagenta+2);  
                histoRatioTrueYieldEtaFullM02->DrawCopy("e1,same");            
                legendBGRatio2->AddEntry(histoRatioTrueYieldEtaFullM02,"#eta","p");
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
    if ( doubleAddFactorK0s >= 0 && nameMeson.Contains("Pi0") ){ //&& !kIsMC
        TCanvas* canvasRAWYieldSec = new TCanvas("canvasRAWYieldSec","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRAWYieldSec, 0.07, 0.01, 0.02, 0.08); 
        canvasRAWYieldSec->SetLogy(1);
        DrawAutoGammaMesonHistos( histoUnCorrectedYieldDrawing, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "RAW Yield/ #it{N}_{Evt}", 
                                    kTRUE, 1, 0.01/nEvt, kTRUE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
        histoUnCorrectedYieldDrawing->SetLineWidth(0.5); 
        histoUnCorrectedYieldDrawing->GetYaxis()->SetTitleOffset(0.8);
        DrawGammaSetMarker(histoUnCorrectedYieldDrawing, 20, 1., kBlack, kBlack);
        histoUnCorrectedYieldDrawing->Draw("e1");
        histoYieldSecMeson->Scale(1./nEvt);
        DrawGammaSetMarker(histoYieldSecMeson, 20, 1., kBlue, kBlue);
        histoYieldSecMeson->DrawCopy("same,e1");  
        histoYieldSecFromK0SMeson->Scale(1./nEvt);
        DrawGammaSetMarker(histoYieldSecFromK0SMeson, 24, 1., kCyan, kCyan);  
        histoYieldSecFromK0SMeson->DrawCopy("same,e1");  

        TLegend* legendSecRAWYield = GetAndSetLegend2(0.6, 0.8, 0.93, 0.93, 28);
        legendSecRAWYield->SetMargin(0.15);
        legendSecRAWYield->AddEntry(histoUnCorrectedYieldDrawing,"RAW yield");
        legendSecRAWYield->AddEntry(histoYieldSecMeson,"total secondaries");
        legendSecRAWYield->AddEntry(histoYieldSecFromK0SMeson,"additional secondaries from K^{0}_{s}");
        legendSecRAWYield->Draw();

        PutProcessLabelAndEnergyOnPlot(0.65, 0.78, 28, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 63, 0.03);
        
        canvasRAWYieldSec->Update();
        canvasRAWYieldSec->SaveAs(Form("%s/%s_%s_RAWYieldSecPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        delete canvasRAWYieldSec;
    } 
    
    //***********************************************************************************************
    //*********************************** correction for yield **************************************
    //***********************************************************************************************
    cout << "calculating corrected yield" << endl;    
    
    TH1D* histoCorrectedYieldTrue = (TH1D*)histoUnCorrectedYield->Clone();
    histoCorrectedYieldTrue->SetName("CorrectedYieldTrueEff");
    TH1D* histoCorrectedYieldTrueAlter = (TH1D*)histoUnCorrectedYield->Clone();
    histoCorrectedYieldTrueAlter->SetName("CorrectedYieldTrueEff");
    
    CorrectYield(histoCorrectedYieldTrue, histoMesonPurityPt, histoYieldSecMeson, histoYieldSecFromK0SMeson, histoTrueEffiPrimMesonPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
    CorrectYield(histoCorrectedYieldTrueAlter, histoPurityPt, histoYieldSecMeson, histoYieldSecFromK0SMeson, histoTrueEffiPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
    
    histoCorrectedYieldTrue->Multiply(histoRatioEffi);
    
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
    DrawGammaSetMarker(histoCorrectedYieldTrueAlter, 24, 1., kAzure+2, kAzure+2);  
    histoCorrectedYieldTrueAlter->DrawCopy("same,e1");  
//     DrawGammaSetMarker(histoMCYieldMeson, 24, 1., kRed+2, kRed+2);  
//     histoMCYieldMeson->DrawCopy("same,e1");  

    TLegend* legendYield3 = GetAndSetLegend2(0.15, 0.03, 0.66, 0.19, 28);
    legendYield3->SetMargin(0.15);
    legendYield3->AddEntry(histoCorrectedYieldTrue,"corrected yield");
    legendYield3->AddEntry(histoCorrectedYieldTrueAlter,"corrected yield, altern. corr facs");
//     legendYield3->AddEntry(histoMCYieldMeson,"corrected yield, altern. corr facs");
    legendYield3->Draw();

    PutProcessLabelAndEnergyOnPlot(0.62, 0.95, 0.03, collisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data());
    
    padCorrectedYieldRatios->cd();
    padCorrectedYieldRatios->SetTickx();
    padCorrectedYieldRatios->SetTicky();
    padCorrectedYieldRatios->SetLogy(0);
    TH1D *ratioTrue = (TH1D*) histoCorrectedYieldTrue->Clone(); 
    ratioTrue->Divide(ratioTrue,histoCorrectedYieldTrue,1.,1.,"");
    TH1D *ratioTrueAlter = (TH1D*) histoCorrectedYieldTrueAlter->Clone(); 
    ratioTrueAlter->Divide(ratioTrueAlter,histoCorrectedYieldTrue,1.,1.,"");
    TH1D *ratioTrueMCInput = (TH1D*) histoMCYieldMeson->Clone();
    ratioTrueMCInput->Divide(histoMCYieldMeson, histoCorrectedYieldTrue,1.,1.,"");

    ratioTrue->SetYTitle("#frac{modified}{standard}"); 
    SetStyleHistoTH1ForGraphs(ratioTrue, "#it{p}_{T} (GeV/#it{c})","#frac{modified}{standard}", 0.07,0.1, 0.07,0.11, 0.9,0.4, 510, 505);
    ratioTrue->GetYaxis()->SetRangeUser(0.8,1.23);
    DrawGammaSetMarker(ratioTrue, 20, 1., kBlack, kBlack); 
    ratioTrue->DrawCopy("p,e1");  
//     DrawGammaSetMarker(ratioTrueMCInput, 24, 1., kRed+2, kRed+2); 
//     ratioTrueMCInput->DrawCopy("same,p,e1");  
    DrawGammaSetMarker(ratioTrueAlter, 24, 1., kAzure+2, kAzure+2); 
    ratioTrueAlter->DrawCopy("same,p,e1");  
    
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
        DrawGammaSetMarker(histoCorrectedYieldTrueAlter, 24, 1., kAzure+2, kAzure+2);  
        histoCorrectedYieldTrueAlter->DrawCopy("same,e1");  
        
        cout << "here" << endl; 
        TLegend* legendYield4 = GetAndSetLegend2(0.15, 0.03, 0.66, 0.19, 28);
        legendYield4->SetMargin(0.15);
        legendYield4->AddEntry(histoCorrectedYieldTrue,"corr true eff");
        legendYield4->AddEntry(histoCorrectedYieldTrueAlter,"corrected yield, altern. corr facs");
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
        DrawGammaSetMarker(ratioTrueAlter, 24, 1., kAzure+2, kAzure+2); 
        ratioTrueAlter->DrawCopy("same,p,e1");  
                
        DrawGammaSetMarker(ratioTrueMCInput, 24, 1., kRed+2, kRed+2);
        ratioTrueMCInput->DrawCopy("e1,same"); 

        canvasCorrecftedYield->Update();
        canvasCorrecftedYield->SaveAs(Form("%s/%s_%s_CorrectedYield_SanityCheck_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(),  fCutSelection.Data(), suffix.Data()));
    }

    delete canvasCorrecftedYield;
    delete legendYield3;
    
    // ********************************************************************************************************************************
    // ****************************** Write file with all further needed histograms ***************************************************
    // ********************************************************************************************************************************    
    cout << "writing corrected file" << endl;    
    
    const char* nameOutput = Form("%s/%s/%s_%s_GammaMergedCorrection%s_%s.root",fCutSelection.Data(),optionEnergy.Data(),nameMeson.Data(),prefix2.Data(),optionPeriod.Data(),fCutSelection.Data());
    TFile* correctedOutput = new TFile(nameOutput,"RECREATE");  

    if (histoCorrectedYieldTrue)            histoCorrectedYieldTrue->Write();
    
    if (histoYieldSecMeson)                 histoYieldSecMeson->Write();
    if (histoYieldSecFromK0SMeson)          histoYieldSecFromK0SMeson->Write();
    
    if (histoUnCorrectedYield)              histoUnCorrectedYield->Write();

    if (histoAcceptance)                    histoAcceptance->Write();
    if (histoTrueEffiPt)                    histoTrueEffiPt->Write("TrueMesonEffiPt");
    if (histoTrueEffiPrimMesonPt)           histoTrueEffiPrimMesonPt->Write();
    if (histoPurityPt)                      histoPurityPt->Write();
    if (histoMesonPurityPt)                 histoMesonPurityPt->Write();
    if (histoEventQuality)                  histoEventQuality->Write();
    if (histoInputMesonPt)                  histoInputMesonPt->Write();
    if (histoMCYieldMeson)                  histoMCYieldMeson->Write();
    if (histoMCYieldMesonOldBin)            histoMCYieldMesonOldBin->Write();
    if (histoUnCorrectedYieldDrawing){
        histoUnCorrectedYieldDrawing->SetName("histoYieldMesonPerEvent");
        histoUnCorrectedYieldDrawing->Write();
    }
    if (deltaPt)                            deltaPt->Write("deltaPt");
    
    correctedOutput->Write();
    correctedOutput->Close();
    
}
