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
    // TString name;
};

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


void CorrectYield(TH1D* histoCorrectedYield, TH1D* histoRawSecYield, TH1D* histoRawAddSecYieldFromK0s, TH1D* histoEffiPt, TH1D* histoAcceptance, Double_t deltaRapid, Double_t scaling, Double_t nEvt, TString nameMeson){
    histoCorrectedYield->Sumw2();
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
                        Int_t mode                      = 9
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
    Bool_t kIsMC = kFALSE;
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
    Double_t deltaRapid     =  ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
    
    Int_t fExampleBin       = 2;
    if(optionEnergy.CompareTo("7TeV") == 0){
        if (nameMeson.Contains("Eta") )fExampleBin=6;
            else fExampleBin=7;
    } else if( optionEnergy.CompareTo("8TeV") == 0) {
        if (nameMeson.Contains("Eta") )fExampleBin=6;
            else fExampleBin=7;
    } else if( optionEnergy.CompareTo("2.76TeV") == 0) {
        if (nameMeson.Contains("Eta") )fExampleBin=4;
            else fExampleBin=7;  
    } else if( optionEnergy.CompareTo("900GeV") == 0) {
        if (nameMeson.Contains("Eta") )fExampleBin=2;
            else fExampleBin=7;
    } else if( optionEnergy.CompareTo("pPb_5.023TeV") == 0 ) {
        if (nameMeson.Contains("Eta") )fExampleBin=4;
            else fExampleBin=7;
    } else if( optionEnergy.CompareTo("PbPb_2.76TeV") == 0 ) {
        if (nameMeson.Contains("Eta") )fExampleBin=6;
            else fExampleBin=7;  
    }
    
    //Variable defintion
    Double_t scaling        = 1./(2.*TMath::Pi());
    
    //*******************************************************************************************************
    //***********************************Reading data file **************************************************
    //*******************************************************************************************************
    // File definitions
    TFile fileUncorrected(fileNameUnCorrectedFile.Data());  
    if (fileUncorrected.IsZombie()) return;
    TH1F *histoNumberOfGoodESDTracksVtx             = (TH1F*)fileUncorrected.Get("GoodESDTracks");
    TH1D *histoEventQuality                         = (TH1D*)fileUncorrected.Get("NEvents");
    TH1D *histoUnCorrectedYield                     = (TH1D*)fileUncorrected.Get("histoYieldMeson");
    TH1D *histoUnCorrectedYieldWide                 = (TH1D*)fileUncorrected.Get("histoYieldMesonWide");
    TH1D *histoUnCorrectedYieldNarrow               = (TH1D*)fileUncorrected.Get("histoYieldMesonNarrow");
    TH1D *histoUnCorrectedYieldLeft                 = (TH1D*)fileUncorrected.Get("histoYieldMesonLeft");
    TH1D *histoUnCorrectedYieldLeftWide             = (TH1D*)fileUncorrected.Get("histoYieldMesonLeftWide");
    TH1D *histoUnCorrectedYieldLeftNarrow           = (TH1D*)fileUncorrected.Get("histoYieldMesonLeftNarrow");
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
    
    // set min and max values for pT
    Double_t maxPtMeson     = histoUnCorrectedYield->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield->GetNbinsX());
    Double_t minPtMeson     = 0;
    Int_t ptBin             = 1;
    while (histoUnCorrectedYield->GetBinContent(ptBin) == 0.){
        ptBin++;
        minPtMeson          = histoUnCorrectedYield->GetXaxis()->GetBinLowEdge(ptBin);
    }
    
    Double_t minPtMesonSec  = 0.3;
    if (mode == 2)
        minPtMesonSec       = minPtMeson;
    else if (mode == 4)
        minPtMesonSec       = minPtMeson;
    cout << "minimum pT: " << minPtMeson << ", maximum pT: " << maxPtMeson << endl;
    
    //*******************************************************************************************************
    //***********************************Reading additional MC file file ************************************
    //*******************************************************************************************************
    
    TFile* fileCorrectionsSecPi0 = new TFile("ExternalInput/PCM/SecondaryFractionHistogramms7TeV.root");
    if (fileCorrectionsSecPi0->IsZombie()) return;
    TH1D* histoDefaultTrueSecFracMeson              = (TH1D*)fileCorrectionsSecPi0->Get("TrueSecFrac");
    TH1D* histoDefaultTrueSecFracMesonWide          = (TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracWide");
    TH1D* histoDefaultTrueSecFracMesonNarrow        = (TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracNarrow");
    TH1D* histoDefaultTrueSecFracFromK0SMeson       = (TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracFromK0S");
    TH1D* histoDefaultTrueSecFracFromK0SMesonWide   = (TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracFromK0SWide");
    TH1D* histoDefaultTrueSecFracFromK0SMesonNarrow = (TH1D*)fileCorrectionsSecPi0->Get("TrueSecFracFromK0SNarrow");
    
    TH1D* histoDefaultTrueSecFracNotFromK0sMeson    = (TH1D*)histoDefaultTrueSecFracMeson->Clone("histoDefaultTrueSecFracNotFromK0sMeson");
    histoDefaultTrueSecFracNotFromK0sMeson->Sumw2();
    histoDefaultTrueSecFracNotFromK0sMeson->Add(histoDefaultTrueSecFracFromK0SMeson,-1);
    TH1D* histoDefaultTrueSecFracModeledToData      = (TH1D*)histoDefaultTrueSecFracFromK0SMeson->Clone("histoDefaultTrueSecFracModeledToData");
    histoDefaultTrueSecFracModeledToData->Sumw2();
    histoDefaultTrueSecFracModeledToData->Scale(1+doubleAddFactorK0s);
    histoDefaultTrueSecFracModeledToData->Add(histoDefaultTrueSecFracNotFromK0sMeson,1);
    
    TF1* fitDefaultSecFrac                          = new TF1("fitDefaultSecFrac","[0]/pow(x,[1])");
    fitDefaultSecFrac->SetRange(minPtMesonSec, maxPtMeson);
    TFitResultPtr resultSecFrac                     = histoDefaultTrueSecFracMeson->Fit(fitDefaultSecFrac,"SINRME+","",minPtMesonSec, maxPtMeson);
    TF1* fitDefaultSecFracWide                      = new TF1("fitDefaultSecFrac","[0]/pow(x,[1])");
    fitDefaultSecFracWide->SetRange(minPtMesonSec, maxPtMeson);
    TFitResultPtr resultSecFracWide                 = histoDefaultTrueSecFracMesonWide->Fit(fitDefaultSecFracWide,"SINRME+","",minPtMesonSec, maxPtMeson);
    TF1* fitDefaultSecFracNarrow                    = new TF1("fitDefaultSecFrac","[0]/pow(x,[1])");
    fitDefaultSecFracNarrow->SetRange(minPtMesonSec, maxPtMeson);
    TFitResultPtr resultSecFracNarrow               = histoDefaultTrueSecFracMesonNarrow->Fit(fitDefaultSecFracNarrow,"SINRME+","",minPtMesonSec, maxPtMeson);
    TF1* fitDefaultSecFracFromK0                    = new TF1("fitDefaultSecFrac","[0]/pow(x,[1])");
    fitDefaultSecFracFromK0->SetRange(minPtMesonSec, maxPtMeson);
    TFitResultPtr resultSecFracFromK0               = histoDefaultTrueSecFracFromK0SMeson->Fit(fitDefaultSecFracFromK0,"SINRME+","",minPtMesonSec, maxPtMeson);
    TF1* fitDefaultSecFracFromK0Wide                = new TF1("fitDefaultSecFrac","[0]/pow(x,[1])");
    fitDefaultSecFracFromK0Wide->SetRange(minPtMesonSec, maxPtMeson);
    TFitResultPtr resultSecFracFromK0Wide           = histoDefaultTrueSecFracFromK0SMesonWide->Fit(fitDefaultSecFracFromK0Wide,"SINRME+","",minPtMesonSec, maxPtMeson);
    TF1* fitDefaultSecFracFromK0Narrow              = new TF1("fitDefaultSecFrac","[0]/pow(x,[1])");
    fitDefaultSecFracFromK0Narrow->SetRange(minPtMesonSec, maxPtMeson);
    TFitResultPtr resultSecFracFromK0Narrow         = histoDefaultTrueSecFracFromK0SMesonNarrow->Fit(fitDefaultSecFracFromK0Narrow,"SINRME+","",minPtMesonSec, maxPtMeson);


    //*******************************************************************************************************
    //***********************************Reading MC correction file *****************************************
    //*******************************************************************************************************    
    TFile* fileCorrections =         new TFile(fileNameCorrectionFile.Data());
    if (fileCorrections->IsZombie()) return;
    TH1F *histoEventQualityMC                       = (TH1F*)fileCorrections->Get("NEvents");
    TH1D *histoEffiPt                               = (TH1D*)fileCorrections->Get("MesonEffiPt"); //not yet correct MesonEffiPt
    TH1D *histoEffiNarrowPt                         = (TH1D*)fileCorrections->Get("MesonNarrowEffiPt");
    TH1D *histoEffiWidePt                           = (TH1D*)fileCorrections->Get("MesonWideEffiPt");
    TH1D *histoEffiLeftPt                           = (TH1D*)fileCorrections->Get("MesonLeftEffiPt");
    TH1D *histoEffiLeftNarrowPt                     = (TH1D*)fileCorrections->Get("MesonLeftNarrowEffiPt");
    TH1D *histoEffiLeftWidePt                       = (TH1D*)fileCorrections->Get("MesonLeftWideEffiPt");
    TH1D *histoAcceptance                           = (TH1D*)fileCorrections->Get("fMCMesonAccepPt");
    
    // load histograms without weighting for comparison of efficiencies & enable the comparison if those are present in the input file
    Bool_t containsWOWeights                        = kFALSE;
    TH1D* histoAcceptanceWOWeights                  = NULL;
    TH1D *histoTrueEffiPtWOWeights                  = NULL;
    TH1D *histoTrueEffiNarrowPtWOWeights            = NULL;
    TH1D *histoTrueEffiWidePtWOWeights              = NULL;
    histoAcceptanceWOWeights                        = (TH1D*)fileCorrections->Get("fMCMesonAccepPtWOWeights"); 
    histoTrueEffiPtWOWeights                        = (TH1D*)fileCorrections->Get("TrueMesonEffiPtUnweighted"); 
    histoTrueEffiNarrowPtWOWeights                  = (TH1D*)fileCorrections->Get("TrueMesonNarrowEffiPtUnweighted");
    histoTrueEffiWidePtWOWeights                    = (TH1D*)fileCorrections->Get("TrueMesonWideEffiPtUnweighted");
    if (histoAcceptanceWOWeights && histoTrueEffiPtWOWeights && histoTrueEffiNarrowPtWOWeights && histoTrueEffiWidePtWOWeights) 
        containsWOWeights                           = kTRUE;
    
    TH1D *histoTrueEffiPt                           = NULL;
    TH1D *histoTrueEffiNarrowPt                     = NULL;
    TH1D *histoTrueEffiWidePt                       = NULL;
    histoTrueEffiPt                                 = (TH1D*)fileCorrections->Get("TrueMesonEffiPt"); 
    histoTrueEffiNarrowPt                           = (TH1D*)fileCorrections->Get("TrueMesonNarrowEffiPt");
    histoTrueEffiWidePt                             = (TH1D*)fileCorrections->Get("TrueMesonWideEffiPt");

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
    
    TH1D *histoTrueEffiPtFixed                      = (TH1D*)histoTrueEffiPt->Clone("histoTrueEffiPtFixed");
    TH1D *histoTrueEffiNarrowPtFixed                = (TH1D*)histoTrueEffiNarrowPt->Clone("TrueMesonNarrowEffiPtFixed");
    TH1D *histoTrueEffiWidePtFixed                  = (TH1D*)histoTrueEffiWidePt->Clone("TrueMesonWideEffiPt");
    TH1D *histoEffiPtFixed                          = (TH1D*)histoEffiPt->Clone("histoEffiPtFixed");
    TH1D *histoEffiNarrowPtFixed                    = (TH1D*)histoEffiNarrowPt->Clone("MesonNarrowEffiPtFixed");
    TH1D *histoEffiWidePtFixed                      = (TH1D*)histoEffiWidePt->Clone("MesonWideEffiPt");

    histoTrueEffiPtFixed                            = FixEfficiency(histoTrueEffiPtFixed,histoTrueEffiPt,optionEnergy,centralityString);
    histoTrueEffiNarrowPtFixed                      = FixEfficiency(histoTrueEffiNarrowPtFixed,histoTrueEffiNarrowPt,optionEnergy,centralityString);
    histoTrueEffiWidePtFixed                        = FixEfficiency(histoTrueEffiWidePtFixed,histoTrueEffiWidePt,optionEnergy,centralityString);
    histoEffiPtFixed                                = FixEfficiency(histoEffiPtFixed,histoEffiPt,optionEnergy,centralityString);
    histoEffiNarrowPtFixed                          = FixEfficiency(histoEffiNarrowPtFixed,histoEffiNarrowPt,optionEnergy,centralityString);
    histoEffiWidePtFixed                            = FixEfficiency(histoEffiWidePtFixed,histoEffiWidePt,optionEnergy,centralityString);

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

    TH1D *histoYieldTrueSecFracMeson                = NULL;
    TH1D *histoYieldTrueSecFracMesonWide            = NULL;
    TH1D *histoYieldTrueSecFracMesonNarrow          = NULL;
    TH1D *histoYieldTrueSecFracFromK0SMesonNarrow   = NULL;
    TH1D *histoYieldTrueSecFracFromK0SMeson         = NULL;
    TH1D *histoYieldTrueSecFracFromK0SMesonWide     = NULL;
    TH1D *histoYieldTrueSecFracMeson_orig           = NULL;
    TH1D *histoYieldTrueSecFracFromK0SMeson_orig    = NULL;
    TH1D *histoYieldTrueSecFracFromLambdaMeson      = NULL;
    TH1D *histoYieldTrueSecFracFromLambdaMeson_orig = NULL;
    TF1* fitSecFracPurePowerlaw                     = new TF1("fitSecFracPurePowerlaw","[0]/pow(x,[1])");
    TF1* fitSecFracPurePowerlawWide                 = new TF1("fitSecFracPurePowerlawWide","[0]/pow(x,[1])");
    TF1* fitSecFracPurePowerlawNarrow               = new TF1("fitSecFracPurePowerlawNarrow","[0]/pow(x,[1])");
    TF1* fitSecFracPurePowerlawFromK0               = new TF1("fitSecFracPurePowerlawFromK0s","[0]/pow(x,[1])");
    TF1* fitSecFracPurePowerlawFromK0Wide           = new TF1("fitSecFracPurePowerlawFromK0Wide","[0]/pow(x,[1])");
    TF1* fitSecFracPurePowerlawFromK0Narrow         = new TF1("fitSecFracPurePowerlawFromK0Narrow","[0]/pow(x,[1])");
    TF1* fitSecFracPLWithConst                      = new TF1("fitSecFracPLWithConst","[0]/pow(x,[1])+[2]");
    fitSecFracPLWithConst->SetParLimits(1,0,1000);
    TF1* fitSecFracPLWithConstWide                  = new TF1("fitSecFracPLWithConstWide","[0]/pow(x,[1])+[2]");
    fitSecFracPLWithConstWide->SetParLimits(1,0,1000);
    TF1* fitSecFracPLWithConstNarrow                = new TF1("fitSecFracPLWithConstNarrow","[0]/pow(x,[1])+[2]");
    fitSecFracPLWithConstNarrow->SetParLimits(1,0,1000);
    TF1* fitSecFracPLWithConstFromK0                = new TF1("fitSecFracPLWithConstFromK0","[0]/pow(x,[1])+[2]");
    fitSecFracPLWithConstFromK0->SetParLimits(1,0,1000);
    TF1* fitSecFracPLWithConstFromK0Wide            = new TF1("fitSecFracPLWithConstFromK0Wide","[0]/pow(x,[1])+[2]");
    fitSecFracPLWithConstFromK0Wide->SetParLimits(1,0,1000);
    TF1* fitSecFracPLWithConstFromK0Narrow          = new TF1("fitSecFracPLWithConstFromK0Narrow","[0]/pow(x,[1])+[2]");
    fitSecFracPLWithConstFromK0Narrow->SetParLimits(1,0,1000);
    TH1D* histoYieldSecMesonLeft                    = NULL;
    TH1D* histoYieldSecMeson                        = NULL;
    TH1D* histoYieldSecMesonLeftNarrow              = NULL;
    TH1D* histoYieldSecMesonNarrow                  = NULL;
    TH1D* histoYieldSecMesonLeftWide                = NULL;
    TH1D* histoYieldSecMesonWide                    = NULL;
    TH1D* histoYieldSecFromK0SMeson                 = NULL;
    TH1D* histoYieldSecFromK0SMesonLeft             = NULL;
    TH1D* histoYieldSecFromK0SMesonLeftNarrow       = NULL;
    TH1D* histoYieldSecFromK0SMesonLeftWide         = NULL;
    TH1D* histoYieldSecFromK0SMesonNarrow           = NULL;
    TH1D* histoYieldSecFromK0SMesonWide             = NULL;
    TH1D *histoYieldTrueGGFracMeson                 = NULL;
    TH1D *histoYieldTrueGGFracMesonWide             = NULL;
    TH1D *histoYieldTrueGGFracMesonNarrow           = NULL;
    TH1D* histoYieldGGMesonLeft                     = NULL;
    TH1D* histoYieldGGMeson                         = NULL;
    TH1D* histoYieldGGMesonLeftNarrow               = NULL;
    TH1D* histoYieldGGMesonNarrow                   = NULL;
    TH1D* histoYieldGGMesonLeftWide                 = NULL;
    TH1D* histoYieldGGMesonWide                     = NULL;
    

    Bool_t doK0SecCorrection                        = kFALSE;
    Int_t doK0SecCorrectionWithDefaultHisto         = 0;
    Bool_t doGGCorrection                           = kFALSE;
    
    if ( !kIsEta ) {
        doK0SecCorrection = kTRUE;
        if ( (mode == 0 || mode == 9) && !kIsMC && kCollisionSystem != 2 )                  doK0SecCorrectionWithDefaultHisto = 1;
        if ( (mode == 0 || mode == 9 || mode == 1) && !kIsMC && kCollisionSystem == 2 )     doK0SecCorrectionWithDefaultHisto = 2;
    }
    if (optDalitz) {
        doGGCorrection                        = kTRUE;
        doK0SecCorrection                     = kFALSE;
        doK0SecCorrectionWithDefaultHisto     = kFALSE;
    }    
    
    if (doK0SecCorrection){
        histoYieldTrueSecFracMeson                  = (TH1D*)fileCorrections->Get("TrueSecFrac");
        histoYieldTrueSecFracMeson_orig             = (TH1D*)histoYieldTrueSecFracMeson->Clone("TrueSecFrac_orig");
        histoYieldTrueSecFracMesonWide              = (TH1D*)fileCorrections->Get("TrueSecFracWide");
        histoYieldTrueSecFracMesonNarrow            = (TH1D*)fileCorrections->Get("TrueSecFracNarrow");
        histoYieldTrueSecFracFromK0SMeson           = (TH1D*)fileCorrections->Get("TrueSecFracFromK0S");
        histoYieldTrueSecFracFromLambdaMeson        = (TH1D*)fileCorrections->Get("TrueSecFracFromLambda");
        histoYieldTrueSecFracFromK0SMeson_orig      = (TH1D*)histoYieldTrueSecFracFromK0SMeson->Clone("TrueSecFracFromK0S_orig");
        histoYieldTrueSecFracFromLambdaMeson_orig   = (TH1D*)histoYieldTrueSecFracFromLambdaMeson->Clone("TrueSecFracFromLambda_orig");
        histoYieldTrueSecFracFromK0SMesonWide       = (TH1D*)fileCorrections->Get("TrueSecFracFromK0SWide");
        histoYieldTrueSecFracFromK0SMesonNarrow     = (TH1D*)fileCorrections->Get("TrueSecFracFromK0SNarrow");
        if (doK0SecCorrectionWithDefaultHisto == 1){
            for (Int_t i = 2; i < histoYieldTrueSecFracMeson->GetNbinsX()+1; i++){
                Double_t ptStart        = histoYieldTrueSecFracMeson->GetXaxis()->GetBinLowEdge(i);
                Double_t ptEnd          = histoYieldTrueSecFracMeson->GetXaxis()->GetBinUpEdge(i);
                Double_t binWidth       = ptEnd-ptStart;
                Double_t secFrac        = fitDefaultSecFrac->Integral(ptStart, ptEnd, resultSecFrac->GetParams()) / binWidth;
                Double_t errorSecFrac   = fitDefaultSecFrac->IntegralError(ptStart, ptEnd, resultSecFrac->GetParams(), resultSecFrac->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracMeson->SetBinContent(i, secFrac);
                histoYieldTrueSecFracMeson->SetBinError(i, errorSecFrac);
                
                secFrac                 = fitDefaultSecFracWide->Integral(ptStart, ptEnd, resultSecFracWide->GetParams()) / binWidth;
                errorSecFrac            = fitDefaultSecFracWide->IntegralError(ptStart, ptEnd, resultSecFracWide->GetParams(), resultSecFracWide->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracMesonWide->SetBinContent(i, secFrac);
                histoYieldTrueSecFracMesonWide->SetBinError(i, errorSecFrac);
                
                secFrac                 = fitDefaultSecFracNarrow->Integral(ptStart, ptEnd, resultSecFracNarrow->GetParams()) / binWidth;
                errorSecFrac            = fitDefaultSecFracNarrow->IntegralError(ptStart, ptEnd, resultSecFracNarrow->GetParams(), resultSecFracNarrow->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracMesonNarrow->SetBinContent(i, secFrac);
                histoYieldTrueSecFracMesonNarrow->SetBinError(i, errorSecFrac);
                
                secFrac                 = fitDefaultSecFracFromK0->Integral(ptStart, ptEnd, resultSecFracFromK0->GetParams()) / binWidth;
                errorSecFrac            = fitDefaultSecFracFromK0->IntegralError(ptStart, ptEnd, resultSecFracFromK0->GetParams(), resultSecFracFromK0->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracFromK0SMeson->SetBinContent(i, secFrac);
                histoYieldTrueSecFracFromK0SMeson->SetBinError(i, errorSecFrac);
                
                secFrac                 = fitDefaultSecFracFromK0Wide->Integral(ptStart, ptEnd, resultSecFracFromK0Wide->GetParams()) / binWidth;
                errorSecFrac            = fitDefaultSecFracFromK0Wide->IntegralError(ptStart, ptEnd, resultSecFracFromK0Wide->GetParams(), resultSecFracFromK0Wide->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracFromK0SMesonWide->SetBinContent(i, secFrac);
                histoYieldTrueSecFracFromK0SMesonWide->SetBinError(i, errorSecFrac);
                
                secFrac                 = fitDefaultSecFracFromK0Narrow->Integral(ptStart, ptEnd, resultSecFracFromK0Narrow->GetParams()) / binWidth;
                errorSecFrac            = fitDefaultSecFracFromK0Narrow->IntegralError(ptStart, ptEnd, resultSecFracFromK0Narrow->GetParams(), resultSecFracFromK0Narrow->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracFromK0SMesonNarrow->SetBinContent(i, secFrac);
                histoYieldTrueSecFracFromK0SMesonNarrow->SetBinError(i, errorSecFrac);
            }    
        } else if ( doK0SecCorrectionWithDefaultHisto == 2 ){
            cout << "Fit1"<< endl;  
            fitSecFracPurePowerlaw->SetRange(minPtMesonSec, maxPtMeson);
            TFitResultPtr resultpPbSecFrac = histoYieldTrueSecFracMeson->Fit(fitSecFracPurePowerlaw,"SINRME+","",minPtMesonSec, maxPtMeson);
            cout << "Fit2"<< endl;  
            fitSecFracPurePowerlawWide->SetRange(minPtMesonSec, maxPtMeson);
            TFitResultPtr resultpPbSecFracWide = histoYieldTrueSecFracMesonWide->Fit(fitSecFracPurePowerlawWide,"SINRME+","",minPtMesonSec, maxPtMeson);

            fitSecFracPurePowerlawNarrow->SetRange(minPtMesonSec, maxPtMeson);
            TFitResultPtr resultpPbSecFracNarrow = histoYieldTrueSecFracMesonNarrow->Fit(fitSecFracPurePowerlawNarrow,"SINRME+","",minPtMesonSec, maxPtMeson);
            cout << "Fit3"<< endl;
            fitSecFracPurePowerlawFromK0->SetRange(minPtMesonSec, maxPtMeson);
            TFitResultPtr resultpPbSecFracFromK0 = histoYieldTrueSecFracFromK0SMeson->Fit(fitSecFracPurePowerlawFromK0,"SINRME+","",minPtMesonSec, maxPtMeson);

            fitSecFracPurePowerlawFromK0Wide->SetRange(minPtMesonSec, maxPtMeson);
            TFitResultPtr resultpPbSecFracFromK0Wide =histoYieldTrueSecFracFromK0SMesonWide->Fit(fitSecFracPurePowerlawFromK0Wide,"SINRME+","",minPtMesonSec, maxPtMeson);

            fitSecFracPurePowerlawFromK0Narrow->SetRange(minPtMesonSec, maxPtMeson);
            TFitResultPtr resultpPbSecFracFromK0Narrow = histoYieldTrueSecFracFromK0SMesonNarrow->Fit(fitSecFracPurePowerlawFromK0Narrow,"SINRME+","",minPtMesonSec, maxPtMeson);
            
            for (Int_t i = 2; i < histoYieldTrueSecFracMeson->GetNbinsX()+1; i++){
                Double_t ptStart        = histoYieldTrueSecFracMeson->GetXaxis()->GetBinLowEdge(i);
                Double_t ptEnd          = histoYieldTrueSecFracMeson->GetXaxis()->GetBinUpEdge(i);
                Double_t binWidth       = ptEnd-ptStart;
                Double_t secFrac        = fitSecFracPurePowerlaw->Integral(ptStart, ptEnd, resultpPbSecFrac->GetParams()) / binWidth;
                Double_t errorSecFrac   = fitSecFracPurePowerlaw->IntegralError(ptStart, ptEnd, resultpPbSecFrac->GetParams(), resultpPbSecFrac->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracMeson->SetBinContent(i, secFrac);
                histoYieldTrueSecFracMeson->SetBinError(i, errorSecFrac);
                    
                secFrac                 = fitSecFracPurePowerlawWide->Integral(ptStart, ptEnd, resultpPbSecFracWide->GetParams()) / binWidth;
                errorSecFrac            = fitSecFracPurePowerlawWide->IntegralError(ptStart, ptEnd, resultpPbSecFracWide->GetParams(), resultpPbSecFracWide->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracMesonWide->SetBinContent(i, secFrac);
                histoYieldTrueSecFracMesonWide->SetBinError(i, errorSecFrac);
                    
                secFrac                 = fitSecFracPurePowerlawNarrow->Integral(ptStart, ptEnd, resultpPbSecFracNarrow->GetParams()) / binWidth;
                errorSecFrac            = fitSecFracPurePowerlawNarrow->IntegralError(ptStart, ptEnd, resultpPbSecFracNarrow->GetParams(), resultpPbSecFracNarrow->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracMesonNarrow->SetBinContent(i, secFrac);
                histoYieldTrueSecFracMesonNarrow->SetBinError(i, errorSecFrac);
                    
                secFrac                 = fitSecFracPurePowerlawFromK0->Integral(ptStart, ptEnd, resultpPbSecFracFromK0->GetParams()) / binWidth;
                errorSecFrac            = fitSecFracPurePowerlawFromK0->IntegralError(ptStart, ptEnd, resultpPbSecFracFromK0->GetParams(), resultpPbSecFracFromK0->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracFromK0SMeson->SetBinContent(i, secFrac);
                histoYieldTrueSecFracFromK0SMeson->SetBinError(i, errorSecFrac);
                    
                secFrac                 = fitSecFracPurePowerlawFromK0Wide->Integral(ptStart, ptEnd, resultpPbSecFracFromK0Wide->GetParams()) / binWidth;
                errorSecFrac            = fitSecFracPurePowerlawFromK0Wide->IntegralError(ptStart, ptEnd, resultpPbSecFracFromK0Wide->GetParams(), resultpPbSecFracFromK0Wide->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracFromK0SMesonWide->SetBinContent(i, secFrac);
                histoYieldTrueSecFracFromK0SMesonWide->SetBinError(i, errorSecFrac);
                    
                secFrac                 = fitSecFracPurePowerlawFromK0Narrow->Integral(ptStart, ptEnd, resultpPbSecFracFromK0Narrow->GetParams()) / binWidth;
                errorSecFrac            = fitSecFracPurePowerlawFromK0Narrow->IntegralError(ptStart, ptEnd, resultpPbSecFracFromK0Narrow->GetParams(), resultpPbSecFracFromK0Narrow->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracFromK0SMesonNarrow->SetBinContent(i, secFrac);
                histoYieldTrueSecFracFromK0SMesonNarrow->SetBinError(i, errorSecFrac);
            }
        } else if ( doK0SecCorrectionWithDefaultHisto == 0 ){            
            cout << "Fit standard range - secondaries"<< endl;  
            fitSecFracPLWithConst->SetRange(minPtMesonSec, maxPtMeson);
            TFitResultPtr resultPbPbLHC11hSecFrac               = histoYieldTrueSecFracMeson->Fit(fitSecFracPLWithConst,"SINRME+","",minPtMesonSec, maxPtMeson);
            
            cout << "Fit wide range - secondaries"<< endl;  
            fitSecFracPLWithConstWide->SetRange(minPtMesonSec, maxPtMeson);
            TFitResultPtr resultPbPbLHC11hSecFracWide           = histoYieldTrueSecFracMesonWide->Fit(fitSecFracPLWithConstWide,"SINRME+","",minPtMesonSec, maxPtMeson);

            cout << "Fit narrow range - secondaries"<< endl;  
            fitSecFracPLWithConstNarrow->SetRange(minPtMesonSec, maxPtMeson);
            TFitResultPtr resultPbPbLHC11hSecFracNarrow         = histoYieldTrueSecFracMesonNarrow->Fit(fitSecFracPLWithConstNarrow,"SINRME+","",minPtMesonSec, maxPtMeson);
            
            cout << "Fit standard range - from K0s"<< endl;  
            fitSecFracPLWithConstFromK0->SetRange(minPtMesonSec, maxPtMeson);
            TFitResultPtr resultPbPbLHC11hSecFracFromK0         = histoYieldTrueSecFracFromK0SMeson->Fit(fitSecFracPLWithConstFromK0,"SINRME+","",minPtMesonSec, maxPtMeson);

            cout << "Fit wide range - from K0s"<< endl;  
            fitSecFracPLWithConstFromK0Wide->SetRange(minPtMesonSec, maxPtMeson);
            TFitResultPtr resultPbPbLHC11hSecFracFromK0Wide     = histoYieldTrueSecFracFromK0SMesonWide->Fit(fitSecFracPLWithConstFromK0Wide,"SINRME+","",minPtMesonSec, maxPtMeson);

            cout << "Fit narrow range - from K0s"<< endl;  
            fitSecFracPLWithConstFromK0Narrow->SetRange(minPtMesonSec, maxPtMeson);
            TFitResultPtr resultPbPbLHC11hSecFracFromK0Narrow = histoYieldTrueSecFracFromK0SMesonNarrow->Fit(fitSecFracPLWithConstFromK0Narrow,"SINRME+","",minPtMesonSec, maxPtMeson);
            
            for (Int_t i = 2; i < histoYieldTrueSecFracMeson->GetNbinsX()+1; i++){
                Double_t ptStart        = histoYieldTrueSecFracMeson->GetXaxis()->GetBinLowEdge(i);
                Double_t ptEnd          = histoYieldTrueSecFracMeson->GetXaxis()->GetBinUpEdge(i);
                Double_t binWidth       = ptEnd-ptStart;
                Double_t secFrac        = fitSecFracPLWithConst->Integral(ptStart, ptEnd, resultPbPbLHC11hSecFrac->GetParams()) / binWidth;
                Double_t errorSecFrac   = fitSecFracPLWithConst->IntegralError(ptStart, ptEnd, resultPbPbLHC11hSecFrac->GetParams(), resultPbPbLHC11hSecFrac->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracMeson->SetBinContent(i, secFrac);
                histoYieldTrueSecFracMeson->SetBinError(i, errorSecFrac);
                    
                secFrac                 = fitSecFracPLWithConstWide->Integral(ptStart, ptEnd, resultPbPbLHC11hSecFracWide->GetParams()) / binWidth;
                errorSecFrac            = fitSecFracPLWithConstWide->IntegralError(ptStart, ptEnd, resultPbPbLHC11hSecFracWide->GetParams(), resultPbPbLHC11hSecFracWide->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracMesonWide->SetBinContent(i, secFrac);
                histoYieldTrueSecFracMesonWide->SetBinError(i, errorSecFrac);
                    
                secFrac                 = fitSecFracPLWithConstNarrow->Integral(ptStart, ptEnd, resultPbPbLHC11hSecFracNarrow->GetParams()) / binWidth;
                errorSecFrac            = fitSecFracPLWithConstNarrow->IntegralError(ptStart, ptEnd, resultPbPbLHC11hSecFracNarrow->GetParams(), resultPbPbLHC11hSecFracNarrow->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracMesonNarrow->SetBinContent(i, secFrac);
                histoYieldTrueSecFracMesonNarrow->SetBinError(i, errorSecFrac);
                    
                secFrac                 = fitSecFracPLWithConstFromK0->Integral(ptStart, ptEnd, resultPbPbLHC11hSecFracFromK0->GetParams()) / binWidth;
                errorSecFrac            = fitSecFracPLWithConstFromK0->IntegralError(ptStart, ptEnd, resultPbPbLHC11hSecFracFromK0->GetParams(), resultPbPbLHC11hSecFracFromK0->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracFromK0SMeson->SetBinContent(i, secFrac);
                histoYieldTrueSecFracFromK0SMeson->SetBinError(i, errorSecFrac);
                    
                secFrac                 = fitSecFracPLWithConstFromK0Wide->Integral(ptStart, ptEnd, resultPbPbLHC11hSecFracFromK0Wide->GetParams()) / binWidth;
                errorSecFrac            = fitSecFracPLWithConstFromK0Wide->IntegralError(ptStart, ptEnd, resultPbPbLHC11hSecFracFromK0Wide->GetParams(), resultPbPbLHC11hSecFracFromK0Wide->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracFromK0SMesonWide->SetBinContent(i, secFrac);
                histoYieldTrueSecFracFromK0SMesonWide->SetBinError(i, errorSecFrac);
                    
                secFrac                 = fitSecFracPLWithConstFromK0Narrow->Integral(ptStart, ptEnd, resultPbPbLHC11hSecFracFromK0Narrow->GetParams()) / binWidth;
                errorSecFrac            = fitSecFracPLWithConstFromK0Narrow->IntegralError(ptStart, ptEnd, resultPbPbLHC11hSecFracFromK0Narrow->GetParams(), resultPbPbLHC11hSecFracFromK0Narrow->GetCovarianceMatrix().GetMatrixArray() ) / binWidth;
                histoYieldTrueSecFracFromK0SMesonNarrow->SetBinContent(i, secFrac);
                histoYieldTrueSecFracFromK0SMesonNarrow->SetBinError(i, errorSecFrac);
            }
        }
        
        histoYieldSecMeson                          = (TH1D*)histoUnCorrectedYield->Clone("SecFracMeson");
        histoYieldSecMeson->Sumw2();
        histoYieldSecMeson->Multiply(histoYieldTrueSecFracMeson);
        histoYieldSecMesonLeft                      = (TH1D*)histoUnCorrectedYieldLeft->Clone("SecFracMesonLeft");
        histoYieldSecMesonLeft->Sumw2();
        histoYieldSecMesonLeft->Multiply(histoYieldTrueSecFracMeson);
        histoYieldSecMesonNarrow                    = (TH1D*)histoUnCorrectedYieldNarrow->Clone("SecFracMesonNarrow");
        histoYieldSecMesonNarrow->Sumw2();
        histoYieldSecMesonNarrow->Multiply(histoYieldTrueSecFracMesonNarrow);
        histoYieldSecMesonLeftNarrow                = (TH1D*)histoUnCorrectedYieldLeftNarrow->Clone("SecFracMesonLeftNarrow");
        histoYieldSecMesonLeftNarrow->Sumw2();
        histoYieldSecMesonLeftNarrow->Multiply(histoYieldTrueSecFracMesonNarrow);
        histoYieldSecMesonWide                      = (TH1D*)histoUnCorrectedYieldWide->Clone("SecFracMesonWide");
        histoYieldSecMesonWide->Sumw2();
        histoYieldSecMesonWide->Multiply(histoYieldTrueSecFracMesonWide);
        histoYieldSecMesonLeftWide                  = (TH1D*)histoUnCorrectedYieldLeftWide->Clone("SecFracMesonLeftWide");
        histoYieldSecMesonLeftWide->Sumw2();
        histoYieldSecMesonLeftWide->Multiply(histoYieldTrueSecFracMesonWide);

        histoYieldSecFromK0SMeson                   = (TH1D*)histoUnCorrectedYield->Clone("SecFracFromK0SMeson");
        histoYieldSecFromK0SMeson->Sumw2();
        histoYieldSecFromK0SMeson->Multiply(histoYieldTrueSecFracFromK0SMeson);
        histoYieldSecFromK0SMeson->Scale(doubleAddFactorK0s);
        histoYieldSecFromK0SMesonLeft               = (TH1D*)histoUnCorrectedYieldLeft->Clone("SecFracFromK0SMesonLeft");
        histoYieldSecFromK0SMesonLeft->Sumw2();
        histoYieldSecFromK0SMesonLeft->Multiply(histoYieldTrueSecFracFromK0SMeson);
        histoYieldSecFromK0SMesonLeft->Scale(doubleAddFactorK0s);
        histoYieldSecFromK0SMesonNarrow             = (TH1D*)histoUnCorrectedYieldNarrow->Clone("SecFracFromK0SMesonNarrow");
        histoYieldSecFromK0SMesonNarrow->Sumw2();
        histoYieldSecFromK0SMesonNarrow->Multiply(histoYieldTrueSecFracFromK0SMesonNarrow);
        histoYieldSecFromK0SMesonNarrow->Scale(doubleAddFactorK0s);
        histoYieldSecFromK0SMesonLeftNarrow         = (TH1D*)histoUnCorrectedYieldLeftNarrow->Clone("SecFracFromK0SMesonLeftNarrow");
        histoYieldSecFromK0SMesonLeftNarrow->Sumw2();
        histoYieldSecFromK0SMesonLeftNarrow->Multiply(histoYieldTrueSecFracFromK0SMesonNarrow);
        histoYieldSecFromK0SMesonLeftNarrow->Scale(doubleAddFactorK0s);
        histoYieldSecFromK0SMesonWide               = (TH1D*)histoUnCorrectedYieldWide->Clone("SecFracFromK0SMesonWide");
        histoYieldSecFromK0SMesonWide->Sumw2();
        histoYieldSecFromK0SMesonWide->Multiply(histoYieldTrueSecFracFromK0SMesonWide);
        histoYieldSecFromK0SMesonWide->Scale(doubleAddFactorK0s);
        histoYieldSecFromK0SMesonLeftWide           = (TH1D*)histoUnCorrectedYieldLeftWide->Clone("SecFracFromK0SMesonLeftWide");
        histoYieldSecFromK0SMesonLeftWide->Sumw2();
        histoYieldSecFromK0SMesonLeftWide->Multiply(histoYieldTrueSecFracFromK0SMesonWide);
        histoYieldSecFromK0SMesonLeftWide->Scale(doubleAddFactorK0s);
                
    } else if (doGGCorrection){
        histoYieldTrueGGFracMeson                   = (TH1D*)fileCorrections->Get("TrueGGFrac");
        histoYieldTrueGGFracMesonWide               = (TH1D*)fileCorrections->Get("TrueGGFracWide");
        histoYieldTrueGGFracMesonNarrow             = (TH1D*)fileCorrections->Get("TrueGGFracNarrow");
        histoYieldGGMeson                           = (TH1D*)histoUnCorrectedYield->Clone("GGFracMeson");
        histoYieldGGMeson->Sumw2();
        histoYieldGGMeson->Multiply(histoYieldTrueGGFracMeson);
        histoYieldGGMesonLeft                       = (TH1D*)histoUnCorrectedYieldLeft->Clone("GGFracMesonLeft");
        histoYieldGGMesonLeft->Sumw2();
        histoYieldGGMesonLeft->Multiply(histoYieldTrueGGFracMeson);
        histoYieldGGMesonNarrow                     = (TH1D*)histoUnCorrectedYieldNarrow->Clone("GGFracMesonNarrow");
        histoYieldGGMesonNarrow->Sumw2();
        histoYieldGGMesonNarrow->Multiply(histoYieldTrueGGFracMesonNarrow);
        histoYieldGGMesonLeftNarrow                 = (TH1D*)histoUnCorrectedYieldLeftNarrow->Clone("GGFracMesonLeftNarrow");
        histoYieldGGMesonLeftNarrow->Sumw2();
        histoYieldGGMesonLeftNarrow->Multiply(histoYieldTrueGGFracMesonNarrow);
        histoYieldGGMesonWide                       = (TH1D*)histoUnCorrectedYieldWide->Clone("GGFracMesonWide");
        histoYieldGGMesonWide->Sumw2();
        histoYieldGGMesonWide->Multiply(histoYieldTrueGGFracMesonWide);
        histoYieldGGMesonLeftWide                   = (TH1D*)histoUnCorrectedYieldLeftWide->Clone("GGFracMesonLeftWide");
        histoYieldGGMesonLeftWide->Sumw2();
        histoYieldGGMesonLeftWide->Multiply(histoYieldTrueGGFracMesonWide);
    }
    
    Double_t mesonMassExpect = 0;
    if( !kIsEta )     mesonMassExpect = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    if( kIsEta )     mesonMassExpect = TDatabasePDG::Instance()->GetParticle(221)->Mass();
    
    cout << "made it!!" << endl;
    
    
    //*******************************************************************************************************
    //***********************************Reading pileup correction file data ********************************
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
    //***********************************Reading pileup correction file MC **********************************
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
    TH1D* histoMCDCAZUnderMesonAllCat_AllPt                         = NULL;
    TH1D* histoMCDCAZGarbageAllCat_AllPt                            = NULL;
    TH1D* histoMCDCAZTrueBackgroundAllCat_AllPt                     = NULL;
    TH1D* histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt    = NULL;
    TH1D* histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt          = NULL;
    TH1D* histoMCDCAZTruePrimaryMesonDalitzAllCat_AllPt             = NULL;
    TH1D* histoMCDCAZTruePrimaryMesonGammaGammaAllCat_AllPt         = NULL;

    if (kDCAFileMCExists){
        histoMCDCAZUnderMesonAllCat_AllPt                           = (TH1D*)fileDCAAnalysisMonteCarlo->Get("HistDCAZUnderMesonAllCat_AllPt");
        histoMCDCAZGarbageAllCat_AllPt                              = (TH1D*)fileDCAAnalysisMonteCarlo->Get("HistDCAZGarbageAllCat_AllPt");
        histoMCDCAZTrueBackgroundAllCat_AllPt                       = (TH1D*)fileDCAAnalysisMonteCarlo->Get("HistDCAZTrueBackgroundAllCat_AllPt");
        histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt      = (TH1D*)fileDCAAnalysisMonteCarlo->Get("HistDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt");
        histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt            = (TH1D*)fileDCAAnalysisMonteCarlo->Get("HistDCAZTrueSecondaryMesonFromK0sAllCat_AllPt");
        histoMCDCAZTruePrimaryMesonDalitzAllCat_AllPt               = (TH1D*)fileDCAAnalysisMonteCarlo->Get("HistDCAZTruePrimaryMesonDalitzAllCat_AllPt");
        histoMCDCAZTruePrimaryMesonGammaGammaAllCat_AllPt           = (TH1D*)fileDCAAnalysisMonteCarlo->Get("HistDCAZTruePrimaryMesonGammaGammaAllCat_AllPt");
        for (Int_t i = 0; i < 6 ; i++){
            histoMCFracCatvsPt[i]                                   = (TH1D*)fileDCAAnalysisMonteCarlo->Get(Form("fHistFracCat_%i_vsPt",i+1));
            histoMCFracIntHistBGvsPt[i]                             = (TH1D*)fileDCAAnalysisMonteCarlo->Get(Form("fHistFracIntHistBGvsPt_Cat_%i_Variant_1",i+1));
        }
    }
    
    Color_t  colorCat[6]        = { kRed+1, 807, 800, kGreen+2, kCyan+2, kBlue+1};
    Color_t  colorCatMC[6]      = { kRed+3, 807+2, 800+2, kGreen+4, kCyan+4, kBlue+3};
    Style_t  styleCat[6]        = { 20, 21, 29, 33, 20, 21};
    Style_t  styleCatMC[6]      = { 24, 25, 30, 27, 24, 25};

    TH1D *histoBGEstimateA                          = (TH1D*)histoUnCorrectedYield->Clone("histoBGEstimateA");
    TH1D *histoBGEstimateB                          = (TH1D*)histoUnCorrectedYield->Clone("histoBGEstimateB");
    TH1D *histoBGEstimateCatA                       = (TH1D*)histoUnCorrectedYield->Clone("histoBGEstimateCatA");
    TH1D *histoBGEstimateCatC                       = (TH1D*)histoUnCorrectedYield->Clone("histoBGEstimateCatC");
    TH1D *histoBGEstimateCatD                       = (TH1D*)histoUnCorrectedYield->Clone("histoBGEstimateCatD");
    
    // ************************************************************************************************
    // ********************** Plot dca distribution with MC component for all pT **********************
    // ************************************************************************************************
    if (kDCAFileDataExists && kDCAFileMCExists){
        cout << "Plotting dca distribution with MC component for all pT" << endl;
        TCanvas* canvasDCAMCComponents = new TCanvas("canvasDCAMCComponents","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasDCAMCComponents, 0.08, 0.02, 0.02, 0.09);
        canvasDCAMCComponents->SetLogy();
            if (histoDCAZUnderMesonAllCat_AllPt){
                DrawAutoGammaMesonHistos( histoDCAZUnderMesonAllCat_AllPt, 
                                "","dca_{z} #gamma (cm)", "d(dca_{z})/#it{N}_{evt}", 
                                kFALSE, 2.,1, kFALSE,
                                kTRUE,1e-8,10*histoMCDCAZUnderMesonAllCat_AllPt->GetMaximum(), 
                                kTRUE, -6., 6.);
                DrawGammaSetMarker(histoDCAZUnderMesonAllCat_AllPt, 20, 1., kBlack, kBlack);
                histoDCAZUnderMesonAllCat_AllPt->GetYaxis()->SetTitleOffset(0.8);
                histoDCAZUnderMesonAllCat_AllPt->DrawCopy("p,e1"); 

            }
            if (histoMCDCAZUnderMesonAllCat_AllPt){
                DrawGammaSetMarker(histoMCDCAZUnderMesonAllCat_AllPt, 24, 1., kGray, kGray);
                histoMCDCAZUnderMesonAllCat_AllPt->DrawCopy("same,p,e1"); 
            }   
            if (histoMCDCAZTruePrimaryMesonGammaGammaAllCat_AllPt){
                DrawGammaSetMarker(histoMCDCAZTruePrimaryMesonGammaGammaAllCat_AllPt, 20, 1., kRed+2, kRed+2);
                histoMCDCAZTruePrimaryMesonGammaGammaAllCat_AllPt->DrawCopy("same,p,e1"); 
            } 
            if (histoMCDCAZTruePrimaryMesonDalitzAllCat_AllPt){
                DrawGammaSetMarker(histoMCDCAZTruePrimaryMesonDalitzAllCat_AllPt, 20, 1., kGreen+2, kGreen+2);
                histoMCDCAZTruePrimaryMesonDalitzAllCat_AllPt->DrawCopy("same,p,e1"); 
            }
            if (histoMCDCAZTrueBackgroundAllCat_AllPt){
                DrawGammaSetMarker(histoMCDCAZTrueBackgroundAllCat_AllPt, 20, 1., kPink+2, kPink+2);
                histoMCDCAZTrueBackgroundAllCat_AllPt->DrawCopy("same,p,e1"); 
            }
            if (histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt){
                DrawGammaSetMarker(histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt, 20, 1., 807, 807);
                histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt->DrawCopy("same,p,e1"); 
            }
            if (histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt){
                DrawGammaSetMarker(histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt, 20, 1., kViolet+2, kViolet+2);
                histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt->DrawCopy("same,p,e1"); 
            }
            if (histoMCDCAZGarbageAllCat_AllPt){
                DrawGammaSetMarker(histoMCDCAZGarbageAllCat_AllPt, 20, 1., kCyan+2, kCyan+2);
                histoMCDCAZGarbageAllCat_AllPt->DrawCopy("same,p,e1"); 
            }
                
            TLatex *labelEnergy = new TLatex(0.11,0.9,Form("%s",collisionSystem.Data()));
            SetStyleTLatex( labelEnergy, 0.04,4);
            labelEnergy->Draw();

            TLegend* legendDCAMCComponents0 = new TLegend(0.7,0.7,0.85,0.95);
            legendDCAMCComponents0->SetTextSize(0.04);
            legendDCAMCComponents0->SetFillColor(0);
            legendDCAMCComponents0->SetLineColor(0);
            legendDCAMCComponents0->SetNColumns(1);
            if (histoDCAZUnderMesonAllCat_AllPt)                    legendDCAMCComponents0->AddEntry(histoDCAZUnderMesonAllCat_AllPt,"Data","p");
            if (histoMCDCAZUnderMesonAllCat_AllPt)                    legendDCAMCComponents0->AddEntry(histoMCDCAZUnderMesonAllCat_AllPt,"Total MC","p");
            if (histoMCDCAZTruePrimaryMesonGammaGammaAllCat_AllPt)    legendDCAMCComponents0->AddEntry(histoMCDCAZTruePrimaryMesonGammaGammaAllCat_AllPt,Form("Prim"),"p");
            if (histoMCDCAZTruePrimaryMesonDalitzAllCat_AllPt)        legendDCAMCComponents0->AddEntry(histoMCDCAZTruePrimaryMesonDalitzAllCat_AllPt,Form("Dalitz"),"l");
            if (histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt && histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt->GetEntries() > 0)
                                                                    legendDCAMCComponents0->AddEntry(histoMCDCAZTrueSecondaryMesonFromK0sAllCat_AllPt,"Sec. #pi^{0} from K^{0}_{s}","p");
            if (histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt && histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt->GetEntries() > 0)
                                                                    legendDCAMCComponents0->AddEntry(histoMCDCAZTrueSecondaryMesonFromSomethingAllCat_AllPt,"Sec. #pi^{0} from X","p");
            if (histoMCDCAZTrueBackgroundAllCat_AllPt)                legendDCAMCComponents0->AddEntry(histoMCDCAZTrueBackgroundAllCat_AllPt,"#gamma#gamma BG","p");
            if (histoMCDCAZGarbageAllCat_AllPt)                        legendDCAMCComponents0->AddEntry(histoMCDCAZGarbageAllCat_AllPt,"garbage","p");
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
        if(nameMeson.Contains("Eta")){
            DrawAutoGammaMesonHistos( histoFracCatvsPt[0], 
                        "", "#it{p}_{T,#eta} (GeV/#it{c})", "#it{N}_{#eta per cat}/(#it{N}_{#eta}) (%)", 
                        kFALSE, 2.,1e-8, kFALSE,
                        kTRUE, 0, 100., 
                        kFALSE, 0., 10.);  
            
        } else {
            DrawAutoGammaMesonHistos( histoFracCatvsPt[0], 
                                    "", "#it{p}_{T,#pi^{0}} (GeV/#it{c})", "#it{N}_{#pi^{0} per cat}/(#it{N}_{#pi^{0}}) (%)", 
                                    kFALSE, 2.,1e-8, kFALSE,
                                    kTRUE, 0, 100., 
                                    kFALSE, 0., 10.);
        }
        DrawGammaSetMarker(histoFracCatvsPt[0], styleCat[0], 1., colorCat[0], colorCat[0]);
        histoFracCatvsPt[0]->GetYaxis()->SetTitleOffset(0.8);
        histoFracCatvsPt[0]->DrawCopy("p,e1"); 
    
        if(histoMCFracCatvsPt[0]){
            DrawGammaSetMarker(histoMCFracCatvsPt[0], styleCatMC[0], 1., colorCatMC[0], colorCatMC[0]);
            histoMCFracCatvsPt[0]->GetYaxis()->SetTitleOffset(0.8);
            histoMCFracCatvsPt[0]->DrawCopy("same,p,e1"); 
        }
        TLegend* legendFractionCat = new TLegend(0.65,0.7,0.95,0.95);
        legendFractionCat->SetTextSize(0.04);
        legendFractionCat->SetFillColor(0);
        legendFractionCat->SetLineColor(0);
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
        if(nameMeson.Contains("Eta")){
            DrawAutoGammaMesonHistos( histoFracIntHistBGvsPt[0], 
                        "", "#it{p}_{T,#eta} (GeV/#it{c})", "BG/Total (%)", 
                        kFALSE, 2.,1e-8, kFALSE,
                        kTRUE, 0, 20, 
                        kFALSE, 0., 10.);
        } else {
            DrawAutoGammaMesonHistos( histoFracIntHistBGvsPt[0], 
                                    "", "#it{p}_{T,#pi^{0}} (GeV/#it{c})", "BG/Total (%)", 
                                    kFALSE, 2.,1e-8, kFALSE,
                                    kTRUE, 0, 20, 
                                    kFALSE, 0., 10.);
        }
        DrawGammaSetMarker(histoFracIntHistBGvsPt[0], styleCat[0], 1., colorCat[0], colorCat[0]);
        histoFracIntHistBGvsPt[0]->GetYaxis()->SetTitleOffset(0.8);
        histoFracIntHistBGvsPt[0]->DrawCopy("p,e1"); 
    
        if(histoMCFracIntHistBGvsPt[0]){
            DrawGammaSetMarker(histoMCFracIntHistBGvsPt[0], styleCatMC[0], 1., colorCatMC[0], colorCatMC[0]);
            histoMCFracIntHistBGvsPt[0]->DrawCopy("same,p,e1"); 
        }
        TLegend* legendFractionCat = new TLegend(0.65,0.7,0.95,0.95);
        legendFractionCat->SetTextSize(0.04);
        legendFractionCat->SetFillColor(0);
        legendFractionCat->SetLineColor(0);
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
    if (kDCAFileDataExists){
        cout << "Plotting total contamination from  out of bunch pileup with different methods" << endl;
        TCanvas* canvasCorrFrac = new TCanvas("canvasCorrFrac","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasCorrFrac, 0.08, 0.02, 0.02, 0.09);

        canvasCorrFrac->cd();
        if(nameMeson.Contains("Eta")){
            DrawAutoGammaMesonHistos( histoCorrectionFactorsHistvsPt, 
                        "", "#it{p}_{T,#eta} (GeV/#it{c})", "Contamination from Pileup (%)", 
                        kFALSE, 2.,1e-8, kFALSE,
                        kTRUE, 0, 8., 
                        kFALSE, 0., 10.);
        } else {
            DrawAutoGammaMesonHistos( histoCorrectionFactorsHistvsPt, 
                                    "", "#it{p}_{T,#pi^{0}} (GeV/#it{c})", "Contamination from Pileup (%)", 
                                    kFALSE, 2.,1e-8, kFALSE,
                                    kTRUE, 0, 8., 
                                    kFALSE, 0., 10.);
        }
        DrawGammaSetMarker(histoCorrectionFactorsHistvsPt, styleMethod[0], 1.2, colorMethod[0], colorMethod[0]);
        histoCorrectionFactorsHistvsPt->GetYaxis()->SetTitleOffset(0.8);
        histoCorrectionFactorsHistvsPt->DrawCopy("p,e1"); 
        TF1* fitCorrectionFactorsHistvsPt = new TF1("fitCorrectionFactorsHistvsPt","[0]/pow(x,[1])+[2]");
        fitCorrectionFactorsHistvsPt->SetRange(0.4, maxPtMeson);
        TFitResultPtr resultCorrectionFactorsHistvsPt = histoCorrectionFactorsHistvsPt->Fit(fitCorrectionFactorsHistvsPt,"SINRME+","",0.4, maxPtMeson);
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
            resultCorrectionFactorsFitvsPt = histoCorrectionFactorsFitvsPt->Fit(fitCorrectionFactorsFitvsPt,"SINRME+","",0.4, maxPtMeson);
            fitCorrectionFactorsFitvsPt->SetLineColor(colorMethod[1]);
            fitCorrectionFactorsFitvsPt->Draw("same");
        }
        DrawGammaSetMarker(histoCorrectionFactorsHistvsPtCatA, styleMethod[2], 1.2, colorMethod[2], colorMethod[2]);
        histoCorrectionFactorsHistvsPtCatA->DrawCopy("same,p,e1"); 
        TF1* fitCorrectionFactorsHistvsPtCatA = new TF1("fitCorrectionFactorsHistvsPtCatA","[0]/pow(x,[1])+[2]");
        fitCorrectionFactorsHistvsPtCatA->SetRange(0.4, maxPtMeson);
        
        TFitResultPtr resultCorrectionFactorsHistvsPtCatA = histoCorrectionFactorsHistvsPtCatA->Fit(fitCorrectionFactorsHistvsPtCatA,"SINRME+","",0.4, maxPtMeson);
        fitCorrectionFactorsHistvsPtCatA->SetLineColor(colorMethod[2]);
        fitCorrectionFactorsHistvsPtCatA->Draw("same");

        DrawGammaSetMarker(histoCorrectionFactorsHistvsPtCatC, styleMethod[3], 1.2, colorMethod[3], colorMethod[3]);
        histoCorrectionFactorsHistvsPtCatC->DrawCopy("same,p,e1"); 
        TF1* fitCorrectionFactorsHistvsPtCatC = new TF1("fitCorrectionFactorsHistvsPtCatC","[0]/pow(x,[1])+[2]");
        fitCorrectionFactorsHistvsPtCatC->SetRange(0.4, maxPtMeson);
        TFitResultPtr resultCorrectionFactorsHistvsPtCatC = histoCorrectionFactorsHistvsPtCatC->Fit(fitCorrectionFactorsHistvsPtCatC,"SINRME+","",0.4, maxPtMeson);
        fitCorrectionFactorsHistvsPtCatC->SetLineColor(colorMethod[3]);
        fitCorrectionFactorsHistvsPtCatC->Draw("same");

        DrawGammaSetMarker(histoCorrectionFactorsHistvsPtCatD, styleMethod[4], 1.2, colorMethod[4], colorMethod[4]);
        histoCorrectionFactorsHistvsPtCatD->DrawCopy("same,p,e1"); 
        TF1* fitCorrectionFactorsHistvsPtCatD = new TF1("fitCorrectionFactorsHistvsPtCatD","[0]/pow(x,[1])+[2]");
        fitCorrectionFactorsHistvsPtCatD->SetRange(0.4, maxPtMeson);
        TFitResultPtr resultCorrectionFactorsHistvsPtCatD = histoCorrectionFactorsHistvsPtCatD->Fit(fitCorrectionFactorsHistvsPtCatD,"SINRME+","",0.4, maxPtMeson);
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
            Double_t errorBGEstimate = (fitCorrectionFactorsHistvsPt->IntegralError(ptStart, ptEnd, resultCorrectionFactorsHistvsPt->GetParams(), resultCorrectionFactorsHistvsPt->GetCovarianceMatrix().GetMatrixArray() ) / binWidth )/100.;
            histoBGEstimateA->SetBinContent(i, bgEstimate);
            histoBGEstimateA->SetBinError(i, errorBGEstimate);
            
            if (fitCorrectionFactorsFitvsPt){
                bgEstimate = (100-fitCorrectionFactorsFitvsPt->Integral(ptStart, ptEnd, resultCorrectionFactorsFitvsPt->GetParams()) / binWidth )/100.;
                errorBGEstimate = (fitCorrectionFactorsFitvsPt->IntegralError(ptStart, ptEnd, resultCorrectionFactorsFitvsPt->GetParams(), resultCorrectionFactorsFitvsPt->GetCovarianceMatrix().GetMatrixArray() ) / binWidth )/100.;
                histoBGEstimateB->SetBinContent(i, bgEstimate);
                histoBGEstimateB->SetBinError(i, errorBGEstimate);
            }
            
            bgEstimate = (100-fitCorrectionFactorsHistvsPtCatA->Integral(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatA->GetParams()) / binWidth )/100.;
            errorBGEstimate = (fitCorrectionFactorsHistvsPtCatA->IntegralError(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatA->GetParams(), resultCorrectionFactorsHistvsPtCatA->GetCovarianceMatrix().GetMatrixArray() ) / binWidth )/100.;
            histoBGEstimateCatA->SetBinContent(i, bgEstimate);
            histoBGEstimateCatA->SetBinError(i, errorBGEstimate);
            
            bgEstimate = (100-fitCorrectionFactorsHistvsPtCatC->Integral(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatC->GetParams()) / binWidth )/100.;
            errorBGEstimate = (fitCorrectionFactorsHistvsPtCatC->IntegralError(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatC->GetParams(), resultCorrectionFactorsHistvsPtCatC->GetCovarianceMatrix().GetMatrixArray() ) / binWidth )/100.;
            histoBGEstimateCatC->SetBinContent(i, bgEstimate);
            histoBGEstimateCatC->SetBinError(i, errorBGEstimate);
            
            bgEstimate = (100-fitCorrectionFactorsHistvsPtCatD->Integral(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatD->GetParams()) / binWidth )/100.;
            errorBGEstimate = (fitCorrectionFactorsHistvsPtCatD->IntegralError(ptStart, ptEnd, resultCorrectionFactorsHistvsPtCatD->GetParams(), resultCorrectionFactorsHistvsPtCatD->GetCovarianceMatrix().GetMatrixArray() ) / binWidth )/100.;
            histoBGEstimateCatD->SetBinContent(i, bgEstimate);
            histoBGEstimateCatD->SetBinError(i, errorBGEstimate);
        }

        canvasCorrFrac->cd();
        if(nameMeson.Contains("Eta")){
            DrawAutoGammaMesonHistos( histoBGEstimateA, 
                                "", "#it{p}_{T,#eta} (GeV/#it{c})", "Correction factor", 
                                kFALSE, 2.,1e-8, kFALSE,
                                kTRUE, 0.9, 1, 
                                kFALSE, 0., 10.);
        } else {
            DrawAutoGammaMesonHistos( histoBGEstimateA, 
                                "", "#it{p}_{T,#pi^{0}} (GeV/#it{c})", "Correction factor", 
                                kFALSE, 2.,1e-8, kFALSE,
                                kTRUE, 0.9, 1, 
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
    
    //**********************************************************************************
    //******************** Mass Plot *********************************************
    //**********************************************************************************
    
    TH1D* histoRatioRecMass         = NULL;
    TH1D* histoRatioValRecMass      = NULL;
    TH1D* histoRatioRecMassGauss    = NULL; 
    TH1D* histoRatioValRecMassGauss = NULL;
    if (!kIsMC){ 
        TCanvas* canvasMass = new TCanvas("canvasMass","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasMass, 0.10, 0.01, 0.02, 0.10);
        
        if ( !kIsEta ){
            histoMassMeson->GetYaxis()->SetRangeUser(0.130,0.140);
            if (mode == 2 || mode == 4 ) histoMassMeson->GetYaxis()->SetRangeUser(0.120,0.140);
            if (kCollisionSystem == 1 && mode > 1) histoMassMeson->GetYaxis()->SetRangeUser(0.130,0.155);
        } else {
            histoMassMeson->GetYaxis()->SetRangeUser(0.54,0.56);
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
        
        TLegend* legendMass = new TLegend(0.15,0.12,0.5,0.25);
        legendMass->SetTextSize(0.02);
        legendMass->SetFillColor(0);
        legendMass->SetFillStyle(0);
        legendMass->SetLineColor(0);
        legendMass->AddEntry(histoMassMeson,"reconstructed Data");

        if (histoMCrecMassMeson){
            DrawGammaSetMarker(histoMCrecMassMeson, 21, 0.8, kRed-4, kRed-4);
            histoMCrecMassMeson->DrawCopy("same,e1,p"); 
            legendMass->AddEntry(histoMCrecMassMeson,"reconstructed MC");
        }
        if( !kIsEta ) legendMass->AddEntry(histoTrueMassMeson,"True reconstructed #pi^{0}");
        if( kIsEta ) legendMass->AddEntry(histoTrueMassMeson,"True reconstructed #eta");

        
        legendMass->Draw();
        PutProcessLabelAndEnergyOnPlot(0.18, 0.97, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
        
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
            
            TLegend* legendMassRatio = new TLegend(0.15,0.12,0.5,0.25);
            legendMassRatio->SetTextSize(0.02);
            legendMassRatio->SetFillColor(0);
            legendMassRatio->SetFillStyle(0);
            legendMassRatio->SetLineColor(0);
            legendMassRatio->SetNColumns(2);
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
            
            PutProcessLabelAndEnergyOnPlot(0.18, 0.97, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            
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
                if (mode == 2 || mode == 4 ) histoMassMeson->GetYaxis()->SetRangeUser(0.120,0.140);
                if (kCollisionSystem == 1 && mode > 1) histoMassMeson->GetYaxis()->SetRangeUser(0.130,0.155);
            } else {
                histoMassMeson->GetYaxis()->SetRangeUser(0.52,0.58);
            }               

            histoMassMeson->DrawCopy("e1,p"); 
            histoTrueMassMeson->DrawCopy("same,e1,p"); 
                
            TLegend* legendMass4 = new TLegend(0.55,0.12,0.95,0.25);
            legendMass4->SetTextSize(0.02);
            legendMass4->SetFillColor(0);
            legendMass4->SetFillStyle(0);
            legendMass4->SetLineColor(0);
            legendMass4->AddEntry(histoMassMeson,"reconstructed Data");

            if (histoMCrecMassMeson){
                histoMCrecMassMeson->DrawCopy("same,e1,p"); 
                legendMass4->AddEntry(histoMCrecMassMeson,"reconstructed MC");
            }

            if(!kIsEta ) legendMass4->AddEntry(histoTrueMassMeson,"True reconstructed #pi^{0}");
            if(kIsEta ) legendMass4->AddEntry(histoTrueMassMeson,"True reconstructed #eta");
            
            if (histoMassGaussianMeson){
                DrawGammaSetMarker(histoMassGaussianMeson, 20, 0.8, kGray+2, kGray+2);
                histoMassGaussianMeson->DrawCopy("same,e1,p"); 
                legendMass4->AddEntry(histoMassGaussianMeson,"reconstructed Data, pure Gauss"); 
            }
            
            if (histoMCrecMassGaussianMeson){
                DrawGammaSetMarker(histoMCrecMassGaussianMeson, 20, 0.8, kGreen-2, kGreen-2);
                histoMCrecMassGaussianMeson->DrawCopy("same,e1,p"); 
                legendMass4->AddEntry(histoMCrecMassGaussianMeson,"reconstructed MC, pure Gauss");
            }
            
            if (histoTrueMassGaussianMeson){
                DrawGammaSetMarker(histoTrueMassGaussianMeson, 24, 0.8, kGreen+4, kGreen+4);
                histoTrueMassGaussianMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendMass4->AddEntry(histoTrueMassGaussianMeson,"True reconstructed #pi^{0}, pure Gauss");
                if(kIsEta ) legendMass4->AddEntry(histoTrueMassGaussianMeson,"True reconstructed #eta, pure Gauss");
            }
            
            DrawGammaLines(0., maxPtMeson,mesonMassExpect, mesonMassExpect,0.1);
            PutProcessLabelAndEnergyOnPlot(0.55, 0.4, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            legendMass4->Draw();
            canvasMass->Update(); 
            canvasMass->SaveAs(Form("%s/%s_MassComparisonPureGaussian_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
        }
        
        //**********************************************************************************
        //******************** Mass Plot further decomposed for PCM + Calo *****************
        //**********************************************************************************        
        canvasMass->cd();
        if (mode==2 || mode == 3){
            if ( !kIsEta ){
                histoMassMeson->GetYaxis()->SetRangeUser(0.120,0.140);
                if (kCollisionSystem == 1 && mode > 1) histoMassMeson->GetYaxis()->SetRangeUser(0.130,0.155);
            } else {
                histoMassMeson->GetYaxis()->SetRangeUser(0.52,0.58);
            }               

            histoMassMeson->DrawCopy("e1,p"); 
            histoTrueMassMeson->DrawCopy("same,e1,p"); 
                
            TLegend* legendMass2 = new TLegend(0.55,0.12,0.95,0.25);
            legendMass2->SetTextSize(0.02);
            legendMass2->SetFillColor(0);
            legendMass2->SetFillStyle(0);
            legendMass2->SetLineColor(0);
            legendMass2->AddEntry(histoMassMeson,"reconstructed Data");

            if (histoMCrecMassMeson){
                histoMCrecMassMeson->DrawCopy("same,e1,p"); 
                legendMass2->AddEntry(histoMCrecMassMeson,"reconstructed MC");
            }

            if(!kIsEta ) legendMass2->AddEntry(histoTrueMassMeson,"True reconstructed #pi^{0}");
            if(kIsEta ) legendMass2->AddEntry(histoTrueMassMeson,"True reconstructed #eta");
            
            TH1D* histoTrueMassCaloPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueMassMesonCaloPhoton");
            if (histoTrueMassCaloPhotonMeson){
                DrawGammaSetMarker(histoTrueMassCaloPhotonMeson, 25, 0.8, kGreen+2, kGreen+2);
                histoTrueMassCaloPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendMass2->AddEntry(histoTrueMassCaloPhotonMeson,"True reconstructed #pi^{0}, cluster real #gamma");
                if(kIsEta ) legendMass2->AddEntry(histoTrueMassCaloPhotonMeson,"True reconstructed #eta, cluster real #gamma");
            }
            TH1D* histoTrueMassCaloConvPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueMassMesonCaloConvPhoton");
            if (histoTrueMassCaloConvPhotonMeson){
                DrawGammaSetMarker(histoTrueMassCaloConvPhotonMeson, 25, 0.8, kCyan+2, kCyan+2);
                histoTrueMassCaloConvPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendMass2->AddEntry(histoTrueMassCaloConvPhotonMeson,"True reconstructed #pi^{0}, cluster conv #gamma");
                if(kIsEta ) legendMass2->AddEntry(histoTrueMassCaloConvPhotonMeson,"True reconstructed #eta, cluster conv #gamma");
            }
            TH1D* histoTrueMassCaloMergedClusterMeson =          (TH1D*)fileCorrections->Get("histoTrueMassMesonCaloMergedCluster");
            if (histoTrueMassCaloMergedClusterMeson){
                DrawGammaSetMarker(histoTrueMassCaloMergedClusterMeson, 25, 0.8, kViolet+2, kViolet+2);
                histoTrueMassCaloMergedClusterMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendMass2->AddEntry(histoTrueMassCaloMergedClusterMeson,"True reconstructed #pi^{0}, merged cluster #gamma");
                if(kIsEta ) legendMass2->AddEntry(histoTrueMassCaloMergedClusterMeson,"True reconstructed #eta, merged cluster #gamma");
            }
            
            DrawGammaLines(0., maxPtMeson,mesonMassExpect, mesonMassExpect,0.1);
            PutProcessLabelAndEnergyOnPlot(0.55, 0.4, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            legendMass2->Draw();
            canvasMass->Update(); 
            canvasMass->SaveAs(Form("%s/%s_MassAddedInfos_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
        }

        //**********************************************************************************
        //******************** Mass Plot further decomposed for Calo + Calo *****************
        //**********************************************************************************        
        if (mode==4 || mode == 5){
            if (!kIsEta){
                histoMassMeson->GetYaxis()->SetRangeUser(0.120,0.140);
                if (kCollisionSystem == 1 && mode > 1) histoMassMeson->GetYaxis()->SetRangeUser(0.130,0.155);
            } else {
                histoMassMeson->GetYaxis()->SetRangeUser(0.52,0.58);
            }               

            histoMassMeson->DrawCopy("e1,p"); 
            histoTrueMassMeson->DrawCopy("same,e1,p"); 
                
            TLegend* legendMass3 = new TLegend(0.55,0.12,0.95,0.25);
            legendMass3->SetTextSize(0.02);
            legendMass3->SetFillColor(0);
            legendMass3->SetFillStyle(0);
            legendMass3->SetLineColor(0);
            legendMass3->AddEntry(histoMassMeson,"reconstructed Data");

            if (histoMCrecMassMeson){
                histoMCrecMassMeson->DrawCopy("same,e1,p"); 
                legendMass3->AddEntry(histoMCrecMassMeson,"reconstructed MC");
            }

            if(!kIsEta ) legendMass3->AddEntry(histoTrueMassMeson,"True reconstructed #pi^{0}");
            if(kIsEta ) legendMass3->AddEntry(histoTrueMassMeson,"True reconstructed #eta");
            
            TH1D* histoTrueMassCaloPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueMassMesonCaloPhoton");
            if (histoTrueMassCaloPhotonMeson){
                DrawGammaSetMarker(histoTrueMassCaloPhotonMeson, 25, 0.8, kGreen+2, kGreen+2);
                histoTrueMassCaloPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendMass3->AddEntry(histoTrueMassCaloPhotonMeson,"True reconstructed #pi^{0}, #gamma#gamma");
                if(kIsEta ) legendMass3->AddEntry(histoTrueMassCaloPhotonMeson,"True reconstructed #eta, #gamma#gamma");
            }
            TH1D* histoTrueMassCaloConvPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueMassMesonCaloConvPhoton");
            if (histoTrueMassCaloConvPhotonMeson){
                DrawGammaSetMarker(histoTrueMassCaloConvPhotonMeson, 25, 0.8, kCyan+2, kCyan+2);
                histoTrueMassCaloConvPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendMass3->AddEntry(histoTrueMassCaloConvPhotonMeson,"True reconstructed #pi^{0}, #gamma_{conv}#gamma_{conv}");
                if(kIsEta ) legendMass3->AddEntry(histoTrueMassCaloConvPhotonMeson,"True reconstructed #eta, #gamma_{conv}#gamma_{conv}");
            }
            TH1D* histoTrueMassMixedCaloConvPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueMassMesonMixedCaloConvPhoton");
            if (histoTrueMassMixedCaloConvPhotonMeson){
                DrawGammaSetMarker(histoTrueMassMixedCaloConvPhotonMeson, 25, 0.8, kBlue+2, kBlue+2);
                histoTrueMassMixedCaloConvPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendMass3->AddEntry(histoTrueMassMixedCaloConvPhotonMeson,"True reconstructed #pi^{0}, #gamma#gamma_{conv}");
                if(kIsEta ) legendMass3->AddEntry(histoTrueMassMixedCaloConvPhotonMeson,"True reconstructed #eta, #gamma#gamma_{conv}");
            }
            
            DrawGammaLines(0., maxPtMeson,mesonMassExpect, mesonMassExpect,0.1);
            
            PutProcessLabelAndEnergyOnPlot(0.55, 0.4, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
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
        DrawGammaCanvasSettings( canvasFWHM, 0.10, 0.01, 0.04, 0.10);
        
        histoFWHMMeson->Sumw2();
        histoFWHMMeson->Scale(1./2.35);
        DrawAutoGammaMesonHistos( histoFWHMMeson, 
                                    "", "#it{p}_{T} (GeV/#it{c})", Form("FWHM/2.35 for %s in |#it{y}| < %s (GeV/#it{c}^{2})",textMeson.Data(), rapidityRange.Data()), 
                                    kFALSE, 1.5,-20., kFALSE,
                                    kTRUE, -0.004, 0.020, 
                                    kFALSE, 0., 10.);  
        histoFWHMMeson->GetYaxis()->SetNdivisions(510); 
        
        TLegend* legendFWHM = new TLegend(0.15,0.1,0.5,0.2);
        legendFWHM->SetTextSize(0.02);
        legendFWHM->SetFillColor(0);
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
        PutProcessLabelAndEnergyOnPlot(0.18, 0.94, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
        canvasFWHM->Update();
        canvasFWHM->SaveAs(Form("%s/%s_FWHM_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));

        //**********************************************************************************
        //******************** FWHM Plot comparison to pure Gaussian width *****************
        //**********************************************************************************        
        
        if (histoWidthGaussianMeson){
            histoFWHMMeson->GetYaxis()->SetRangeUser(-0.004, 0.050);
            histoFWHMMeson->DrawCopy("e1,p"); 
            histoTrueFWHMMeson->DrawCopy("same,e1,p"); 

            TLegend* legendFWHM4 = new TLegend(0.55,0.12,0.95,0.25);
            legendFWHM4->SetTextSize(0.02);
            legendFWHM4->SetFillColor(0);
            legendFWHM4->SetFillStyle(0);
            legendFWHM4->SetLineColor(0);
            legendFWHM4->AddEntry(histoFWHMMeson,"reconstructed Data");

            if (histoMCrecFWHMMeson){
                histoMCrecFWHMMeson->DrawCopy("same,e1,p"); 
                legendFWHM4->AddEntry(histoMCrecFWHMMeson,"reconstructed MC");
            }

            if(!kIsEta ) legendFWHM4->AddEntry(histoTrueFWHMMeson,"True reconstructed #pi^{0}");
            if(kIsEta ) legendFWHM4->AddEntry(histoTrueFWHMMeson,"True reconstructed #eta");
            
            
            if (histoWidthGaussianMeson){
                DrawGammaSetMarker(histoWidthGaussianMeson, 20, 1.0, kGray+2, kGray+2);
                histoWidthGaussianMeson->DrawCopy("same,e1,p"); 
                legendFWHM4->AddEntry(histoWidthGaussianMeson,"reconstructed Data, #sigma pure Gauss");
            }
            if (histoMCrecWidthGaussMeson){
                DrawGammaSetMarker(histoMCrecWidthGaussMeson, 20, 1.0, kGreen-2, kGreen-2);
                histoMCrecWidthGaussMeson->DrawCopy("same,e1,p"); 
                legendFWHM4->AddEntry(histoMCrecWidthGaussMeson,"reconstructed MC, #sigma pure Gauss");
            }

            if (histoTrueWidthGaussianMeson){
                DrawGammaSetMarker(histoTrueWidthGaussianMeson, 24, 1.0, kGreen+4, kGreen+4);
                histoTrueWidthGaussianMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendFWHM4->AddEntry(histoTrueWidthGaussianMeson,"True reconstructed #pi^{0}, #sigma pure Gauss");
                if(kIsEta ) legendFWHM4->AddEntry(histoTrueWidthGaussianMeson,"True reconstructed #eta, #sigma pure Gauss");
            }
            
            PutProcessLabelAndEnergyOnPlot(0.18, 0.94, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            legendFWHM4->Draw();
            canvasFWHM->Update(); 
            canvasFWHM->SaveAs(Form("%s/%s_FWHMComparisonPureGauss_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
        }
        //**********************************************************************************
        //******************** FWHM Plot further decomposed for PCM + Calo *****************
        //**********************************************************************************        
        if (mode==2 || mode == 3){
        
            histoFWHMMeson->GetYaxis()->SetRangeUser(-0.004, 0.030); 
            histoFWHMMeson->DrawCopy("e1,p"); 
            histoTrueFWHMMeson->DrawCopy("same,e1,p"); 

            TLegend* legendFWHM2 = new TLegend(0.55,0.12,0.95,0.25);
            legendFWHM2->SetTextSize(0.02);
            legendFWHM2->SetFillColor(0);
            legendFWHM2->SetFillStyle(0);
            legendFWHM2->SetLineColor(0);
            legendFWHM2->AddEntry(histoFWHMMeson,"reconstructed Data");

            if (histoMCrecFWHMMeson){
                histoMCrecFWHMMeson->DrawCopy("same,e1,p"); 
                legendFWHM2->AddEntry(histoMCrecFWHMMeson,"reconstructed MC");
            }

            if(!kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMMeson,"True reconstructed #pi^{0}");
            if(kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMMeson,"True reconstructed #eta");
            
            TH1D* histoTrueFWHMCaloPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueFWHMMesonCaloPhoton");
            if (histoTrueFWHMCaloPhotonMeson){
                histoTrueFWHMCaloPhotonMeson->Scale(1./2.35);
                DrawGammaSetMarker(histoTrueFWHMCaloPhotonMeson, 25, 0.8, kGreen+2, kGreen+2);
                histoTrueFWHMCaloPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMCaloPhotonMeson,"True reconstructed #pi^{0}, cluster real #gamma");
                if(kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMCaloPhotonMeson,"True reconstructed #eta, cluster real #gamma");
            }
            TH1D* histoTrueFWHMCaloConvPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueFWHMMesonCaloConvPhoton");
            if (histoTrueFWHMCaloConvPhotonMeson){
                histoTrueFWHMCaloConvPhotonMeson->Scale(1./2.35);
                DrawGammaSetMarker(histoTrueFWHMCaloConvPhotonMeson, 25, 0.8, kCyan+2, kCyan+2);
                histoTrueFWHMCaloConvPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMCaloConvPhotonMeson,"True reconstructed #pi^{0}, cluster conv #gamma");
                if(kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMCaloConvPhotonMeson,"True reconstructed #eta, cluster conv #gamma");
            }
            TH1D* histoTrueFWHMCaloMergedClusterMeson =          (TH1D*)fileCorrections->Get("histoTrueFWHMMesonCaloMergedCluster");
            if (histoTrueFWHMCaloMergedClusterMeson){
                histoTrueFWHMCaloMergedClusterMeson->Scale(1./2.35);
                DrawGammaSetMarker(histoTrueFWHMCaloMergedClusterMeson, 25, 0.8, kViolet+2, kViolet+2);
                histoTrueFWHMCaloMergedClusterMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMCaloMergedClusterMeson,"True reconstructed #pi^{0}, merged cluster #gamma");
                if(kIsEta ) legendFWHM2->AddEntry(histoTrueFWHMCaloMergedClusterMeson,"True reconstructed #eta, merged cluster #gamma");
            }
            
            PutProcessLabelAndEnergyOnPlot(0.18, 0.94, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            legendFWHM2->Draw();
            canvasFWHM->Update(); 
            canvasFWHM->SaveAs(Form("%s/%s_FWHMAddedInfos_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
        }

        //**********************************************************************************
        //******************** FWHM Plot further decomposed for Calo + Calo *****************
        //**********************************************************************************        

        if (mode==4 || mode == 5){
            histoFWHMMeson->GetYaxis()->SetRangeUser(-0.004, 0.030); 
            histoFWHMMeson->DrawCopy("e1,p"); 
            histoTrueFWHMMeson->DrawCopy("same,e1,p"); 
                
            TLegend* legendFWHM3 = new TLegend(0.55,0.12,0.95,0.25);
            legendFWHM3->SetTextSize(0.02);
            legendFWHM3->SetFillColor(0);
            legendFWHM3->SetFillStyle(0);
            legendFWHM3->SetLineColor(0);
            legendFWHM3->AddEntry(histoFWHMMeson,"reconstructed Data");

            if (histoMCrecFWHMMeson){
                histoMCrecFWHMMeson->DrawCopy("same,e1,p"); 
                legendFWHM3->AddEntry(histoMCrecFWHMMeson,"reconstructed MC");
            }

            if(!kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMMeson,"True reconstructed #pi^{0}");
            if(kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMMeson,"True reconstructed #eta");

            TH1D* histoTrueFWHMCaloPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueFWHMMesonCaloPhoton");
            if (histoTrueFWHMCaloPhotonMeson){
                histoTrueFWHMCaloPhotonMeson->Scale(1./2.35);
                DrawGammaSetMarker(histoTrueFWHMCaloPhotonMeson, 25, 0.8, kGreen+2, kGreen+2);
                histoTrueFWHMCaloPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMCaloPhotonMeson,"True reconstructed #pi^{0}, #gamma#gamma");
                if(kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMCaloPhotonMeson,"True reconstructed #eta, #gamma#gamma");
            }
            TH1D* histoTrueFWHMCaloConvPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueFWHMMesonCaloConvPhoton");
            if (histoTrueFWHMCaloConvPhotonMeson){
                histoTrueFWHMCaloConvPhotonMeson->Scale(1./2.35);
                DrawGammaSetMarker(histoTrueFWHMCaloConvPhotonMeson, 25, 0.8, kCyan+2, kCyan+2);
                histoTrueFWHMCaloConvPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMCaloConvPhotonMeson,"True reconstructed #pi^{0}, #gamma_{conv}#gamma_{conv}");
                if(kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMCaloConvPhotonMeson,"True reconstructed #eta, #gamma_{conv}#gamma_{conv}");
            }
            TH1D* histoTrueFWHMMixedCaloConvPhotonMeson =          (TH1D*)fileCorrections->Get("histoTrueFWHMMesonMixedCaloConvPhoton");
            if (histoTrueFWHMMixedCaloConvPhotonMeson){
                histoTrueFWHMMixedCaloConvPhotonMeson->Scale(1./2.35);
                DrawGammaSetMarker(histoTrueFWHMMixedCaloConvPhotonMeson, 25, 0.8, kViolet+2, kViolet+2);
                histoTrueFWHMMixedCaloConvPhotonMeson->DrawCopy("same,e1,p"); 
                if(!kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMMixedCaloConvPhotonMeson,"True reconstructed #pi^{0}, #gamma#gamma_{conv}");
                if(kIsEta ) legendFWHM3->AddEntry(histoTrueFWHMMixedCaloConvPhotonMeson,"True reconstructed #eta, #gamma#gamma_{conv}");
            }
            
            legendFWHM3->Draw();
            PutProcessLabelAndEnergyOnPlot(0.18, 0.94, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());

            canvasFWHM->Update(); 
            canvasFWHM->SaveAs(Form("%s/%s_FWHMAddedInfos_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
            delete legendFWHM3;
        }

        
        delete canvasFWHM;
        delete legendFWHM;
    }    

    TH1D* histoTrueEffiPtUnmod                         = (TH1D*) histoTrueEffiPt->Clone("histoTrueEffiPtUnmod"); 
    TH1D* histoTrueEffiNarrowPtUnmod                 = (TH1D*) histoTrueEffiNarrowPt->Clone("histoTrueEffiNarrowPtUnmod"); 
    TH1D* histoTrueEffiWidePtUnmod                 = (TH1D*) histoTrueEffiWidePt->Clone("histoTrueEffiWidePtUnmod"); 
    if (containsWOWeights && (mode!=0 || mode!=9) && nameMeson.Contains("Pi0")){

        TCanvas* canvasCompEffSimple = new TCanvas("canvasCompEffSimple","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasCompEffSimple, 0.10, 0.01, 0.035, 0.09);

        TH1D* histoRatioEffWOWeightingNormalEff        = (TH1D*) histoEffiPt->Clone(); 
        histoRatioEffWOWeightingNormalEff->Divide(histoRatioEffWOWeightingNormalEff, histoTrueEffiPtWOWeights, 1., 1., "B");
        TH1D* histoRatioEffWOWeightingNormalEffNarrow        = (TH1D*) histoEffiNarrowPt->Clone(); 
        histoRatioEffWOWeightingNormalEffNarrow->Divide(histoRatioEffWOWeightingNormalEffNarrow, histoTrueEffiNarrowPtWOWeights, 1., 1., "B");
        TH1D* histoRatioEffWOWeightingNormalEffWide        = (TH1D*) histoEffiWidePt->Clone(); 
        histoRatioEffWOWeightingNormalEffWide->Divide(histoRatioEffWOWeightingNormalEffWide, histoTrueEffiWidePtWOWeights, 1., 1., "B");

        // Calculation & Plotting of correction factor for Normal integration window
        DrawAutoGammaMesonHistos( histoRatioEffWOWeightingNormalEff, 
                                    "", "#it{p}_{T} (GeV/#it{c})", Form("#epsilon_{eff,%s, rec}/#epsilon_{eff,%s, true wo weights} ", textMeson.Data(), textMeson.Data()), 
                                    kFALSE, 1.3, 3e-6, kFALSE,
                                    kTRUE, 0.8, 1.5, 
                                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoRatioEffWOWeightingNormalEff, 24, 1., 807, 807);
        histoRatioEffWOWeightingNormalEff->Draw("e1");

        
        TF1* fitEffiBiasWOWeightsNormalPol0         = new TF1("fitEffiBiasWOWeightsNormalPol0","[0]",0.4,maxPtMeson);
        TF1* fitEffiBiasWOWeightsNormalPol1         = new TF1("fitEffiBiasWOWeightsNormalPol1","[0]/pow(x,[1])+[2]",0.4,maxPtMeson);
        fitEffiBiasWOWeightsNormalPol1->SetParLimits(2,0.5,1.5);
        
        histoRatioEffWOWeightingNormalEff->Fit(fitEffiBiasWOWeightsNormalPol0,"NRME+","",0.4,maxPtMeson);
        cout << WriteParameterToFile(fitEffiBiasWOWeightsNormalPol0) << endl;
        TH1D* histoRatioEffWOWeightingNormalEffCFPol0 = (TH1D*)histoRatioEffWOWeightingNormalEff->Clone("histoRatioEffWOWeightingNormalEffCFPol0");
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(histoRatioEffWOWeightingNormalEffCFPol0);
        histoRatioEffWOWeightingNormalEffCFPol0->SetStats(kFALSE);
        histoRatioEffWOWeightingNormalEffCFPol0->SetFillColor(806);
        histoRatioEffWOWeightingNormalEffCFPol0->SetMarkerSize(0);
        histoRatioEffWOWeightingNormalEffCFPol0->Draw("e3,same");
        
        histoRatioEffWOWeightingNormalEff->Fit(fitEffiBiasWOWeightsNormalPol1,"NRME+","",0.4,maxPtMeson    );
        cout << WriteParameterToFile(fitEffiBiasWOWeightsNormalPol1) << endl;
        TH1D* histoRatioEffWOWeightingNormalEffCFPol1 = (TH1D*)histoRatioEffWOWeightingNormalEff->Clone("histoRatioEffWOWeightingNormalEffCFPol1");
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(histoRatioEffWOWeightingNormalEffCFPol1);
        histoRatioEffWOWeightingNormalEffCFPol1->SetStats(kFALSE);
        histoRatioEffWOWeightingNormalEffCFPol1->SetFillColor(791);
        histoRatioEffWOWeightingNormalEffCFPol1->SetFillStyle(3003);
        histoRatioEffWOWeightingNormalEffCFPol1->SetMarkerSize(0);
        for (Int_t i=1; i< histoTrueEffiPt->GetNbinsX(); i++){
            if (histoTrueEffiPt->GetBinContent(i) == 0){
                histoRatioEffWOWeightingNormalEffCFPol1->SetBinContent(i,1);
                histoRatioEffWOWeightingNormalEffCFPol1->SetBinContent(i,0);
                histoRatioEffWOWeightingNormalEffCFPol0->SetBinContent(i,1);
                histoRatioEffWOWeightingNormalEffCFPol0->SetBinContent(i,0); 
            }    
        }    
        histoRatioEffWOWeightingNormalEffCFPol1->Draw("e3,same");  
        fitEffiBiasWOWeightsNormalPol0->SetLineColor(807);
        fitEffiBiasWOWeightsNormalPol0->SetLineStyle(1);
        fitEffiBiasWOWeightsNormalPol0->Draw("same");
        fitEffiBiasWOWeightsNormalPol1->SetLineColor(797);
        fitEffiBiasWOWeightsNormalPol1->SetLineStyle(7);
        fitEffiBiasWOWeightsNormalPol1->Draw("same");
        histoRatioEffWOWeightingNormalEff->Draw("e1,same");
        
        PutProcessLabelAndEnergyOnPlot(0.72, 0.25, 28, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data(), 63, 0.03);
        
        canvasCompEffSimple->Update();
        canvasCompEffSimple->SaveAs(Form("%s/%s_EffiCompW0WeightingNormalRatio_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));  
        
        // Calculation & Plotting of correction factor for narrow integration window
        DrawGammaSetMarker(histoRatioEffWOWeightingNormalEffNarrow, 24, 1., kGreen+2, kGreen+2);
        histoRatioEffWOWeightingNormalEffNarrow->Draw("e1,same");

        DrawAutoGammaMesonHistos( histoRatioEffWOWeightingNormalEffNarrow, 
                                    "", "#it{p}_{T} (GeV/#it{c})", Form("#epsilon_{eff,%s, rec}/#epsilon_{eff,%s, true wo weights} ", textMeson.Data(), textMeson.Data()), 
                                    kFALSE, 1.3, 3e-6, kFALSE,
                                    kTRUE, 0.8, 1.5, 
                                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoRatioEffWOWeightingNormalEffNarrow, 24, 1., kGreen+2, kGreen+2);
        histoRatioEffWOWeightingNormalEffNarrow->Draw("e1");

        
        TF1* fitEffiBiasWOWeightsNormalPol0Nar         = new TF1("fitEffiBiasWOWeightsNormalPol0Nar","[0]",0.4,maxPtMeson);
        TF1* fitEffiBiasWOWeightsNormalPol1Nar         = new TF1("fitEffiBiasWOWeightsNormalPol1Nar","[0]/pow(x,[1])+[2]",0.4,maxPtMeson);
        fitEffiBiasWOWeightsNormalPol1Nar->SetParLimits(2,0.5,1.5);
        
        histoRatioEffWOWeightingNormalEffNarrow->Fit(fitEffiBiasWOWeightsNormalPol0Nar,"NRME+","",0.4,maxPtMeson);
        cout << WriteParameterToFile(fitEffiBiasWOWeightsNormalPol0Nar) << endl;
        TH1D* histoRatioEffWOWeightingNormalEffCFPol0Nar = (TH1D*)histoRatioEffWOWeightingNormalEffNarrow->Clone("histoRatioEffWOWeightingNormalEffCFPol0Nar");
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(histoRatioEffWOWeightingNormalEffCFPol0Nar);
        histoRatioEffWOWeightingNormalEffCFPol0Nar->SetStats(kFALSE);
        histoRatioEffWOWeightingNormalEffCFPol0Nar->SetFillColor(kGreen-7);
        histoRatioEffWOWeightingNormalEffCFPol0Nar->SetMarkerSize(0);
        histoRatioEffWOWeightingNormalEffCFPol0Nar->Draw("e3,same");
        
        histoRatioEffWOWeightingNormalEffNarrow->Fit(fitEffiBiasWOWeightsNormalPol1Nar,"NRME+","",0.4,maxPtMeson    );
        cout << WriteParameterToFile(fitEffiBiasWOWeightsNormalPol1Nar) << endl;
        TH1D* histoRatioEffWOWeightingNormalEffCFPol1Nar = (TH1D*)histoRatioEffWOWeightingNormalEffNarrow->Clone("histoRatioEffWOWeightingNormalEffCFPol1Nar");
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(histoRatioEffWOWeightingNormalEffCFPol1Nar);
        histoRatioEffWOWeightingNormalEffCFPol1Nar->SetStats(kFALSE);
        histoRatioEffWOWeightingNormalEffCFPol1Nar->SetFillColor(kGreen-6);
        histoRatioEffWOWeightingNormalEffCFPol1Nar->SetFillStyle(3003);
        histoRatioEffWOWeightingNormalEffCFPol1Nar->SetMarkerSize(0);
        for (Int_t i=1; i< histoTrueEffiPt->GetNbinsX(); i++){
            if (histoTrueEffiPt->GetBinContent(i) == 0){
                histoRatioEffWOWeightingNormalEffCFPol1Nar->SetBinContent(i,1);
                histoRatioEffWOWeightingNormalEffCFPol1Nar->SetBinContent(i,0);
                histoRatioEffWOWeightingNormalEffCFPol0Nar->SetBinContent(i,1);
                histoRatioEffWOWeightingNormalEffCFPol0Nar->SetBinContent(i,0); 
            }    
        }    
        histoRatioEffWOWeightingNormalEffCFPol1Nar->Draw("e3,same");  
        fitEffiBiasWOWeightsNormalPol0Nar->SetLineColor(kGreen+1);
        fitEffiBiasWOWeightsNormalPol0Nar->SetLineStyle(1);
        fitEffiBiasWOWeightsNormalPol0Nar->Draw("same");
        fitEffiBiasWOWeightsNormalPol1Nar->SetLineColor(kGreen+3);
        fitEffiBiasWOWeightsNormalPol1Nar->SetLineStyle(7);
        fitEffiBiasWOWeightsNormalPol1Nar->Draw("same");
        histoRatioEffWOWeightingNormalEffNarrow->Draw("e1,same");
        
        PutProcessLabelAndEnergyOnPlot(0.72, 0.25, 28, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data(), 63, 0.03);
        
        canvasCompEffSimple->Update();
        canvasCompEffSimple->SaveAs(Form("%s/%s_EffiCompW0WeightingNormalRatioNarrow_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));  

        // Calculation & Plotting of correction factor for wide integration window
        DrawAutoGammaMesonHistos( histoRatioEffWOWeightingNormalEffWide, 
                                    "", "#it{p}_{T} (GeV/#it{c})", Form("#epsilon_{eff,%s, rec}/#epsilon_{eff,%s, true wo weights} ", textMeson.Data(), textMeson.Data()), 
                                    kFALSE, 1.3, 3e-6, kFALSE,
                                    kTRUE, 0.5, 1.5, 
                                    kFALSE, 0., 10.);
        DrawGammaSetMarker(histoRatioEffWOWeightingNormalEffWide, 24, 1., kCyan+2, kCyan+2);
        histoRatioEffWOWeightingNormalEffWide->Draw("e1");

        
        TF1* fitEffiBiasWOWeightsNormalPol0Wi         = new TF1("fitEffiBiasWOWeightsNormalPol0Wi","[0]",0.4,maxPtMeson);
        TF1* fitEffiBiasWOWeightsNormalPol1Wi         = new TF1("fitEffiBiasWOWeightsNormalPol1Wi","[0]/pow(x,[1])+[2]",0.4,maxPtMeson);
        fitEffiBiasWOWeightsNormalPol1Wi->SetParLimits(2,0.5,1.5);
        
        histoRatioEffWOWeightingNormalEff->Fit(fitEffiBiasWOWeightsNormalPol0Wi,"NRME+","",0.4,maxPtMeson);
        cout << WriteParameterToFile(fitEffiBiasWOWeightsNormalPol0Wi) << endl;
        TH1D* histoRatioEffWOWeightingNormalEffCFPol0Wi = (TH1D*)histoRatioEffWOWeightingNormalEffWide->Clone("histoRatioEffWOWeightingNormalEffCFPol0Wi");
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(histoRatioEffWOWeightingNormalEffCFPol0Wi);
        histoRatioEffWOWeightingNormalEffCFPol0Wi->SetStats(kFALSE);
        histoRatioEffWOWeightingNormalEffCFPol0Wi->SetFillColor(kCyan-7);
        histoRatioEffWOWeightingNormalEffCFPol0Wi->SetMarkerSize(0);
        histoRatioEffWOWeightingNormalEffCFPol0Wi->Draw("e3,same");
        
        histoRatioEffWOWeightingNormalEff->Fit(fitEffiBiasWOWeightsNormalPol1Wi,"NRME+","",0.4,maxPtMeson    );
        cout << WriteParameterToFile(fitEffiBiasWOWeightsNormalPol1Wi) << endl;
        TH1D* histoRatioEffWOWeightingNormalEffCFPol1Wi = (TH1D*)histoRatioEffWOWeightingNormalEffWide->Clone("histoRatioEffWOWeightingNormalEffCFPol1Wi");
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(histoRatioEffWOWeightingNormalEffCFPol1Wi);
        histoRatioEffWOWeightingNormalEffCFPol1Wi->SetStats(kFALSE);
        histoRatioEffWOWeightingNormalEffCFPol1Wi->SetFillColor(kCyan-6);
        histoRatioEffWOWeightingNormalEffCFPol1Wi->SetFillStyle(3003);
        histoRatioEffWOWeightingNormalEffCFPol1Wi->SetMarkerSize(0);
        for (Int_t i=1; i< histoTrueEffiPt->GetNbinsX(); i++){
            if (histoTrueEffiPt->GetBinContent(i) == 0){
                histoRatioEffWOWeightingNormalEffCFPol1Wi->SetBinContent(i,1);
                histoRatioEffWOWeightingNormalEffCFPol1Wi->SetBinContent(i,0);
                histoRatioEffWOWeightingNormalEffCFPol0Wi->SetBinContent(i,1);
                histoRatioEffWOWeightingNormalEffCFPol0Wi->SetBinContent(i,0); 
            }    
        }    
        histoRatioEffWOWeightingNormalEffCFPol1Wi->Draw("e3,same");  
        fitEffiBiasWOWeightsNormalPol0Wi->SetLineColor(kCyan+1);
        fitEffiBiasWOWeightsNormalPol0Wi->SetLineStyle(1);
        fitEffiBiasWOWeightsNormalPol0Wi->Draw("same");
        fitEffiBiasWOWeightsNormalPol1Wi->SetLineColor(kCyan+3);
        fitEffiBiasWOWeightsNormalPol1Wi->SetLineStyle(7);
        fitEffiBiasWOWeightsNormalPol1Wi->Draw("same");
        histoRatioEffWOWeightingNormalEffWide->Draw("e1,same");
        
        PutProcessLabelAndEnergyOnPlot(0.72, 0.25, 28, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data(), 63, 0.03);
        
        canvasCompEffSimple->Update();
        canvasCompEffSimple->SaveAs(Form("%s/%s_EffiCompW0WeightingNormalRatioWide_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));  
        
        
        
        histoTrueEffiPt->Multiply(histoTrueEffiPt,histoRatioEffWOWeightingNormalEffCFPol1 );
        histoTrueEffiNarrowPt->Multiply(histoTrueEffiNarrowPt,histoRatioEffWOWeightingNormalEffCFPol1Nar );
        histoTrueEffiWidePt->Multiply(histoTrueEffiWidePt,histoRatioEffWOWeightingNormalEffCFPol1Wi );
        
        // plotting of final comparison
        TH1D* histoRatioEffWOWeightingTrueEffCorr        = (TH1D*) histoEffiPt->Clone(); 
        histoRatioEffWOWeightingTrueEffCorr->Divide(histoRatioEffWOWeightingTrueEffCorr, histoTrueEffiPt, 1., 1., "B");

        histoRatioEffWOWeightingNormalEff->Draw("e1");
        DrawGammaSetMarker(histoRatioEffWOWeightingTrueEffCorr, 20, 1.5, kAzure-6, kAzure-6);
        histoRatioEffWOWeightingTrueEffCorr->Draw("same,e1");
        
        canvasCompEffSimple->Update();
        canvasCompEffSimple->SaveAs(Form("%s/%s_EffiCompW0WeightingNormalRatioAfterFix_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));  
        
    }    
    
    //**********************************************************************************
    //******************** Acceptance Plot *********************************************
    //**********************************************************************************
    if (kIsMC){
        TCanvas* canvasAcceptance = new TCanvas("canvasAcceptance2","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasAcceptance, 0.1, 0.01, 0.02, 0.10);

        DrawAutoGammaMesonHistos( histoAcceptance, 
                                    "", "#it{p}_{T} (GeV/#it{c})", Form("A_{%s} in |#it{y}| < %s",textMeson.Data(),rapidityRange.Data()), 
                                    kTRUE, 1.3, 3e-6, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
        if (!kIsEta){
            if (optionEnergy.CompareTo("pPb_5.023TeV")==0 && mode == 0) histoAcceptance->GetYaxis()->SetRangeUser(0.2,1.02);
            else if (mode == 0) histoAcceptance->GetYaxis()->SetRangeUser(0.7,1.02);
            else if (mode == 4 || mode == 2) histoAcceptance->GetYaxis()->SetRangeUser(0.,0.3);
            else histoAcceptance->GetYaxis()->SetRangeUser(0.7,1.02);
        } else {
            if (optionEnergy.CompareTo("pPb_5.023TeV")==0 && mode == 0) histoAcceptance->GetYaxis()->SetRangeUser(0.1,1.);
            else if (mode == 0) histoAcceptance->GetYaxis()->SetRangeUser(0.5,1.);
            else if (mode == 4 || mode == 2) histoAcceptance->GetYaxis()->SetRangeUser(0.,0.3);
            else histoAcceptance->GetYaxis()->SetRangeUser(0.5,1.02);
        }  
                
        DrawGammaSetMarker(histoAcceptance, 20, 1.5, kAzure-6, kAzure-6);
        histoAcceptance->DrawCopy("e1"); 

        PutProcessLabelAndEnergyOnPlot(0.72, 0.25, 28, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data(), 63, 0.03);

        canvasAcceptance->Update();
        canvasAcceptance->SaveAs(Form("%s/%s_Acceptance_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
        
        if (containsWOWeights){
            DrawGammaSetMarker(histoAcceptanceWOWeights, 24, 1., kBlack, kBlack);  
            histoAcceptanceWOWeights->DrawCopy("e1, same"); 

            TLegend* legendAccComp = GetAndSetLegend2(0.45, 0.12, 0.65, 0.22, 28);
            legendAccComp->AddEntry(histoAcceptance,"with weights");
            legendAccComp->AddEntry(histoAcceptanceWOWeights,"without weights");
            legendAccComp->Draw();
        
            
            canvasAcceptance->Update();
            canvasAcceptance->SaveAs(Form("%s/%s_AcceptanceCompWAndW0Weighting_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));

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
        
        delete canvasAcceptance;
        
        //**********************************************************************************
        //******************** Secondary Fraction     **************************************
        //**********************************************************************************        
        if (!kDalitz && nameMeson.Contains("Pi0")){
            TCanvas* canvasSecFrac = new TCanvas("canvasSecFrac","",200,10,1350,900);  // gives the page size
            DrawGammaCanvasSettings( canvasSecFrac, 0.09, 0.02, 0.04, 0.09);
//             canvasSecFrac->SetLogy(1); 
            
//             histoYieldTrueSecFracMeson_orig->Scale(100.);
//             histoYieldTrueSecFracFromK0SMeson_orig->Scale(100.);
//             histoYieldTrueSecFracFromLambdaMeson_orig->Scale(100.);
            DrawAutoGammaMesonHistos( histoYieldTrueSecFracMeson_orig, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#it{r}_{X} = #frac{X->#pi^{0}}{#pi^{0}}", // (%)", 
                                    kTRUE, 1.5, 0, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
            histoYieldTrueSecFracMeson_orig->GetYaxis()->SetTitleOffset(0.9);
            DrawGammaSetMarker(histoYieldTrueSecFracMeson_orig, 20, 1., kBlack, kBlack);  
            DrawGammaSetMarker(histoYieldTrueSecFracFromK0SMeson_orig, 24, 1., kBlue, kBlue);
            DrawGammaSetMarker(histoYieldTrueSecFracFromLambdaMeson_orig, 24, 1., kViolet+2, kViolet+2);
            histoYieldTrueSecFracMeson_orig->DrawCopy("e1");  
            histoYieldTrueSecFracFromK0SMeson_orig->DrawCopy("e1,same");  
            histoYieldTrueSecFracFromLambdaMeson_orig->DrawCopy("e1,same"); 
            if (doK0SecCorrectionWithDefaultHisto == 0) {
                fitSecFracPLWithConst->SetLineColor(kBlack);  
//                 fitSecFracPLWithConst->Scale(100.);
                fitSecFracPLWithConstFromK0->SetLineColor(kBlue); 
//                 fitSecFracPLWithConstFromK0->Scale(100.);
                fitSecFracPLWithConst->Draw("same");
                fitSecFracPLWithConstFromK0->Draw("same");
            } else if (doK0SecCorrectionWithDefaultHisto == 2) {
                fitSecFracPurePowerlaw->SetLineColor(kBlack);  
                fitSecFracPurePowerlawFromK0->SetLineColor(kBlue); 
                fitSecFracPurePowerlaw->Draw("same");
                fitSecFracPurePowerlawFromK0->Draw("same");
            } else if (doK0SecCorrectionWithDefaultHisto == 1) {    
                fitDefaultSecFrac->SetLineColor(kBlack);  
                fitDefaultSecFracFromK0->SetLineColor(kBlue);
                fitDefaultSecFrac->Draw("same");
                fitDefaultSecFracFromK0->Draw("same");
            }
            
            TLegend* legendSecFrac = new TLegend(0.6,0.76,0.94,0.93);
            legendSecFrac->SetTextSize(0.03);
            legendSecFrac->SetLineColor(0);
            legendSecFrac->SetLineWidth(0);
            legendSecFrac->SetFillColor(0);
            legendSecFrac->SetFillStyle(0);
            legendSecFrac->AddEntry(histoYieldTrueSecFracMeson_orig,"#it{r}_{All}");
            if (doK0SecCorrectionWithDefaultHisto == 0) legendSecFrac->AddEntry(fitSecFracPLWithConst,"fit to #it{r}_{All}");
                else if (doK0SecCorrectionWithDefaultHisto == 2) legendSecFrac->AddEntry(fitSecFracPurePowerlaw,"fit to #it{r}_{All}");
                else if (doK0SecCorrectionWithDefaultHisto == 1) legendSecFrac->AddEntry(fitDefaultSecFrac,"fit to #it{r}_{All}");
            legendSecFrac->AddEntry(histoYieldTrueSecFracFromK0SMeson_orig,"#it{r}_{K_{s}^{0}}");
            if (doK0SecCorrectionWithDefaultHisto == 0) legendSecFrac->AddEntry(fitSecFracPLWithConstFromK0,"fit to #it{r}_{K_{s}^{0}}");
                else if (doK0SecCorrectionWithDefaultHisto == 2) legendSecFrac->AddEntry(fitSecFracPurePowerlawFromK0,"fit to #it{r}_{K_{s}^{0}}");
                else if (doK0SecCorrectionWithDefaultHisto == 1) legendSecFrac->AddEntry(fitDefaultSecFracFromK0,"fit to #it{r}_{K_{s}^{0}}");
            legendSecFrac->AddEntry(histoYieldTrueSecFracFromLambdaMeson_orig,"#it{r}_{#Lambda}"); 
            legendSecFrac->Draw();

            PutProcessLabelAndEnergyOnPlot(0.62, 0.75, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
            canvasSecFrac->Update();
            canvasSecFrac->SaveAs(Form("%s/%s_FracSecondaries_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
            delete canvasSecFrac;
        }
        
        //**********************************************************************************
        //******************** Efficiency Simple Plot **************************************
        //**********************************************************************************
        TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasEffSimple, 0.10, 0.01, 0.035, 0.09);
                    
        DrawAutoGammaMesonHistos( histoTrueEffiPtUnmod, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#epsilon_{eff}", 
                                    kTRUE, 1.3, 3e-6, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);
                
        DrawGammaSetMarker(histoTrueEffiPtUnmod, 20, 1., kBlack, kBlack);
        histoTrueEffiPtUnmod->DrawCopy("e1");
        if (containsWOWeights){
            DrawGammaSetMarker(histoEffiPt, 25, 1., kGreen+2, kGreen+2);
            histoEffiPt->DrawCopy("same,e1,p");
            DrawGammaSetMarker(histoTrueEffiPtWOWeights, 24, 1., 807, 807);
            histoTrueEffiPtWOWeights->DrawCopy("same,e1,p");
            DrawGammaSetMarker(histoTrueEffiPt, 21, 1., kBlue+1, kBlue+1);
            histoTrueEffiPt->DrawCopy("same,e1,p");
        } else if (mode == 4){
            DrawGammaSetMarker(histoEffiPt, 25, 1., kGreen+2, kGreen+2);
            histoEffiPt->DrawCopy("same,e1,p");            
        }    
        TLegend* legendEff = GetAndSetLegend2(0.25,0.13,0.45,0.24, 28);
        legendEff->SetMargin(0.15);
        
        if (containsWOWeights){
            legendEff->AddEntry(histoTrueEffiPtUnmod,"validated efficiency, w weights");
            legendEff->AddEntry(histoTrueEffiPtWOWeights,"validated efficiency, w/o weights"); 
            legendEff->AddEntry(histoEffiPt,"reconstructed efficiency, as in Data"); 
            legendEff->AddEntry(histoTrueEffiPt,"corr validated efficiency"); 
        } else if (mode == 4) {
            legendEff->AddEntry(histoEffiPt,"reconstructed efficiency, as in Data"); 
            legendEff->AddEntry(histoTrueEffiPtUnmod,"validated efficiency");            
        } else {
            legendEff->AddEntry(histoTrueEffiPtUnmod,"validated efficiency");
        }    
        legendEff->Draw();
        
        canvasEffSimple->Update();
        PutProcessLabelAndEnergyOnPlot(0.72, 0.25, 28, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data(), 63, 0.03);

        canvasEffSimple->SaveAs(Form("%s/%s_TrueEffSimple_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data())); 
        delete canvasEffSimple;

        //*********************************************************************************
        //********************** True Efficiency Plot ******************************************
        //*********************************************************************************
        
        TCanvas* canvasEffi = new TCanvas("canvasEffi","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasEffi, 0.10, 0.015, 0.035, 0.10);

        DrawAutoGammaMesonHistos( histoTrueEffiPt, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#epsilon_{eff}", 
                                    kFALSE, 0.75, 3e-6, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 10.);

        DrawGammaSetMarker(histoTrueEffiPt, 20, 1., kBlack, kBlack); 
        histoTrueEffiPt->DrawCopy("e1"); 
        
        //right Side Normalization narrow
        DrawGammaSetMarker(histoTrueEffiNarrowPt, 24, 1., kGray+1, kGray+1); 
        histoTrueEffiNarrowPt->DrawCopy("e1,same"); 
        
        //       //right Side Normalization wide
        DrawGammaSetMarker(histoTrueEffiWidePt, 24, 1., kGray+3, kGray+3); 
        histoTrueEffiWidePt->DrawCopy("e1,same"); 
                
        TLegend* legendTrueEff = new TLegend(0.6,0.13,0.98,0.24);
        legendTrueEff->SetTextSize(0.02);
        legendTrueEff->SetFillColor(0);
        legendTrueEff->SetFillStyle(0);
        legendTrueEff->SetLineColor(0);
        legendTrueEff->AddEntry(histoTrueEffiPt,"true normal");
        legendTrueEff->AddEntry(histoTrueEffiWidePt,"true wide int");
        legendTrueEff->AddEntry(histoTrueEffiNarrowPt,"true narrow int");
        legendTrueEff->Draw();
        
        PutProcessLabelAndEnergyOnPlot(0.62, 0.35, 0.02, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
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
                            kFALSE, 0., histoUnCorrectedYield->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield->GetNbinsX()));

    TF1* fitTsallisMC;
    if (nameMeson.CompareTo("Pi0")==0 || nameMeson.CompareTo("Pi0EtaBinning")==0 ){
        fitTsallisMC= FitObject("l","fitTsallisMC","Pi0",histoMCYieldMesonOldBin,0.3,histoUnCorrectedYield->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield->GetNbinsX()),NULL,"QNRME+");
    } else { 
        fitTsallisMC= FitObject("l","fitTsallisMC","Eta",histoMCYieldMesonOldBin,0.3,histoUnCorrectedYield->GetXaxis()->GetBinUpEdge(histoUnCorrectedYield->GetNbinsX()),NULL,"QNRME+");
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

    TH1D* histoUnCorrectedYieldDrawing = (TH1D*)histoUnCorrectedYield->Clone();
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

    //***********************************************************************************************
    //*********************************** correction for yield **************************************
    //***********************************************************************************************
    
    TH1D* histoCorrectedYieldNorm = (TH1D*)histoUnCorrectedYield->Clone();
    histoCorrectedYieldNorm->SetName("CorrectedYieldNormEff");
    TH1D* histoCorrectedYieldNormNarrow = (TH1D*)histoUnCorrectedYieldNarrow->Clone();
    histoCorrectedYieldNormNarrow->SetName("CorrectedYieldNormEffNarrow");
    TH1D* histoCorrectedYieldNormWide = (TH1D*)histoUnCorrectedYieldWide->Clone();
    histoCorrectedYieldNormWide->SetName("CorrectedYieldNormEffWide");
    TH1D* histoCorrectedYieldNormLeft = (TH1D*)histoUnCorrectedYieldLeft->Clone();
    histoCorrectedYieldNormLeft->SetName("CorrectedYieldNormEffLeft");
    TH1D* histoCorrectedYieldNormLeftNarrow = (TH1D*)histoUnCorrectedYieldLeftNarrow->Clone();
    histoCorrectedYieldNormLeftNarrow->SetName("CorrectedYieldNormEffLeftNarrow");
    TH1D* histoCorrectedYieldNormLeftWide = (TH1D*)histoUnCorrectedYieldLeftWide->Clone();
    histoCorrectedYieldNormLeftWide->SetName("CorrectedYieldNormEffLeftWide");
    
    TH1D* histoCorrectedYieldTrue = (TH1D*)histoUnCorrectedYield->Clone();
    histoCorrectedYieldTrue->SetName("CorrectedYieldTrueEff");
    TH1D* histoCorrectedYieldTrueNarrow = (TH1D*)histoUnCorrectedYieldNarrow->Clone();
    histoCorrectedYieldTrueNarrow->SetName("CorrectedYieldTrueEffNarrow");
    TH1D* histoCorrectedYieldTrueWide = (TH1D*)histoUnCorrectedYieldWide->Clone();
    histoCorrectedYieldTrueWide->SetName("CorrectedYieldTrueEffWide");
    TH1D* histoCorrectedYieldFixed = (TH1D*)histoUnCorrectedYield->Clone();
    histoCorrectedYieldFixed->SetName("CorrectedYieldEffFixed");
    TH1D* histoCorrectedYieldNarrowFixed = (TH1D*)histoUnCorrectedYieldNarrow->Clone();
    histoCorrectedYieldNarrowFixed->SetName("CorrectedYieldEffNarrowFixed");
    TH1D* histoCorrectedYieldWideFixed = (TH1D*)histoUnCorrectedYieldWide->Clone();
    histoCorrectedYieldWideFixed->SetName("CorrectedYieldEffWideFixed");
    TH1D* histoCorrectedYieldTrueFixed = (TH1D*)histoUnCorrectedYield->Clone();
    histoCorrectedYieldTrueFixed->SetName("CorrectedYieldTrueEffFixed");
    TH1D* histoCorrectedYieldTrueNarrowFixed = (TH1D*)histoUnCorrectedYieldNarrow->Clone();
    histoCorrectedYieldTrueNarrowFixed->SetName("CorrectedYieldTrueEffNarrowFixed");
    TH1D* histoCorrectedYieldTrueWideFixed = (TH1D*)histoUnCorrectedYieldWide->Clone();
    histoCorrectedYieldTrueWideFixed->SetName("CorrectedYieldTrueEffWideFixed");
    TH1D* histoCorrectedYieldTrueLeft = (TH1D*)histoUnCorrectedYieldLeft->Clone();
    histoCorrectedYieldTrueLeft->SetName("CorrectedYieldTrueEffLeft");
    TH1D* histoCorrectedYieldTrueLeftNarrow = (TH1D*)histoUnCorrectedYieldLeftNarrow->Clone();
    histoCorrectedYieldTrueLeftNarrow->SetName("CorrectedYieldTrueEffLeftNarrow");
    TH1D* histoCorrectedYieldTrueLeftWide = (TH1D*)histoUnCorrectedYieldLeftWide->Clone();
    histoCorrectedYieldTrueLeftWide->SetName("CorrectedYieldTrueEffLeftWide");
    
    TH1D* histoCompleteCorr = (TH1D*)histoTrueEffiPt->Clone();

    if (!optDalitz){
        CorrectYield(histoCorrectedYieldNorm, histoYieldSecMeson, histoYieldSecFromK0SMeson, histoEffiPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
        CorrectYield(histoCorrectedYieldNormNarrow, histoYieldSecMesonNarrow, histoYieldSecFromK0SMesonNarrow, histoEffiNarrowPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
        CorrectYield(histoCorrectedYieldNormWide, histoYieldSecMesonWide, histoYieldSecFromK0SMesonWide, histoEffiWidePt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
        CorrectYield(histoCorrectedYieldNormLeft, histoYieldSecMeson, histoYieldSecFromK0SMeson, histoEffiLeftPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
        CorrectYield(histoCorrectedYieldNormLeftNarrow, histoYieldSecMesonNarrow, histoYieldSecFromK0SMesonNarrow, histoEffiLeftNarrowPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
        CorrectYield(histoCorrectedYieldNormLeftWide, histoYieldSecMesonWide, histoYieldSecFromK0SMesonWide, histoEffiLeftWidePt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
        
        CorrectYield(histoCorrectedYieldTrue, histoYieldSecMeson, histoYieldSecFromK0SMeson, histoTrueEffiPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
        CompileFullCorrectionFactor( histoCompleteCorr, histoAcceptance, deltaRapid);
        CorrectYield(histoCorrectedYieldTrueNarrow, histoYieldSecMesonNarrow, histoYieldSecFromK0SMesonNarrow , histoTrueEffiNarrowPt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
        CorrectYield(histoCorrectedYieldTrueWide, histoYieldSecMesonWide, histoYieldSecFromK0SMesonWide, histoTrueEffiWidePt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
        
        CorrectYield(histoCorrectedYieldFixed, histoYieldSecMeson, histoYieldSecFromK0SMeson, histoEffiPtFixed, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
        CorrectYield(histoCorrectedYieldNarrowFixed, histoYieldSecMesonNarrow, histoYieldSecFromK0SMesonNarrow , histoEffiNarrowPtFixed, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
        CorrectYield(histoCorrectedYieldWideFixed, histoYieldSecMesonWide, histoYieldSecFromK0SMesonWide, histoEffiWidePtFixed, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);

        CorrectYield(histoCorrectedYieldTrueFixed, histoYieldSecMeson, histoYieldSecFromK0SMeson, histoTrueEffiPtFixed, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
        CorrectYield(histoCorrectedYieldTrueNarrowFixed, histoYieldSecMesonNarrow, histoYieldSecFromK0SMesonNarrow , histoTrueEffiNarrowPtFixed, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
        CorrectYield(histoCorrectedYieldTrueWideFixed, histoYieldSecMesonWide, histoYieldSecFromK0SMesonWide, histoTrueEffiWidePtFixed, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
    
        CorrectYield(histoCorrectedYieldTrueLeft, histoYieldSecMesonLeft, histoYieldSecFromK0SMesonLeft, histoTrueEffiPt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
        CorrectYield(histoCorrectedYieldTrueLeftNarrow, histoYieldSecMesonLeftNarrow, histoYieldSecFromK0SMesonLeftNarrow, histoTrueEffiNarrowPt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
        CorrectYield(histoCorrectedYieldTrueLeftWide, histoYieldSecMesonLeftWide, histoYieldSecFromK0SMesonLeftWide, histoTrueEffiWidePt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);

        if (!kIsMC && kDCAFileDataExists){
                histoCorrectedYieldNorm->Multiply(histoBGEstimateCatA);
                histoCorrectedYieldTrue->Multiply(histoBGEstimateCatA);
                histoCorrectedYieldTrueNarrow->Multiply(histoBGEstimateCatA);
                histoCorrectedYieldTrueWide->Multiply(histoBGEstimateCatA);
                histoCorrectedYieldFixed->Multiply(histoBGEstimateCatA);
                histoCorrectedYieldNarrowFixed->Multiply(histoBGEstimateCatA);
                histoCorrectedYieldWideFixed->Multiply(histoBGEstimateCatA);
                histoCorrectedYieldTrueFixed->Multiply(histoBGEstimateCatA);
                histoCorrectedYieldTrueNarrowFixed->Multiply(histoBGEstimateCatA);
                histoCorrectedYieldTrueWideFixed->Multiply(histoBGEstimateCatA);
                histoCorrectedYieldTrueLeft->Multiply(histoBGEstimateCatA);
                histoCorrectedYieldTrueLeftNarrow->Multiply(histoBGEstimateCatA);
                histoCorrectedYieldTrueLeftWide->Multiply(histoBGEstimateCatA);
        }   
    } else {
        CorrectYieldDalitz(histoCorrectedYieldNorm, histoYieldGGMeson, histoEffiPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
        CorrectYieldDalitz(histoCorrectedYieldTrue, histoYieldGGMeson,histoTrueEffiPt, histoAcceptance, deltaRapid, scaling, nEvt, nameMeson);
        CorrectYieldDalitz(histoCorrectedYieldTrueNarrow, histoYieldGGMesonNarrow, histoTrueEffiNarrowPt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
        CorrectYieldDalitz(histoCorrectedYieldTrueWide, histoYieldGGMesonWide, histoTrueEffiWidePt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);

        CorrectYieldDalitz(histoCorrectedYieldTrueLeft, histoYieldGGMesonLeft, histoTrueEffiPt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
        CorrectYieldDalitz(histoCorrectedYieldTrueLeftNarrow, histoYieldGGMesonLeftNarrow,histoTrueEffiNarrowPt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
        CorrectYieldDalitz(histoCorrectedYieldTrueLeftWide, histoYieldGGMesonLeftWide, histoTrueEffiWidePt, histoAcceptance,  deltaRapid, scaling, nEvt, nameMeson);
    } 
    
    // **************************************************************************************
    // ************** Plot corrected yield with differnt yield extraction methods ***********
    // **************************************************************************************
    TCanvas* canvasCorrecftedYield = new TCanvas("canvasCorrecftedYield","",1350,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasCorrecftedYield, 0.13, 0.02, 0.02, 0.09);
    canvasCorrecftedYield->SetLogy();

    TPad* padCorrectedYieldHistos = new TPad("padCorrectedYieldHistos", "", 0., 0.3, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padCorrectedYieldHistos, 0.12, 0.02, 0.02, 0.);
    padCorrectedYieldHistos->Draw();

    TPad* padCorrectedYieldRatios = new TPad("padCorrectedYieldRatios", "", 0., 0., 1., 0.3,-1, -1, -2);
    DrawGammaPadSettings( padCorrectedYieldRatios, 0.12, 0.02, 0., 0.26);
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
    //right Side Normalization narrow
    DrawGammaSetMarker(histoCorrectedYieldTrueNarrow, 24, 1., kGray+1, kGray+1);  
    histoCorrectedYieldTrueNarrow->DrawCopy("e1,same"); 
    //right Side Normalization wide
    DrawGammaSetMarker(histoCorrectedYieldTrueWide, 24, 1., kGray+3, kGray+3);  
    histoCorrectedYieldTrueWide->DrawCopy("e1,same"); 
    //left Side Normalization 
    DrawGammaSetMarker(histoCorrectedYieldTrueLeft, 20, 1., kBlue, kBlue);
    histoCorrectedYieldTrueLeft->DrawCopy("e1,same"); 
    //left Side Normalization narrow
    DrawGammaSetMarker(histoCorrectedYieldTrueLeftNarrow, 24, 1., kBlue-5, kBlue-5); 
    histoCorrectedYieldTrueLeftNarrow->DrawCopy("e1,same"); 
    //left Side Normalization wide
    DrawGammaSetMarker(histoCorrectedYieldTrueLeftWide, 24, 1., kBlue+2, kBlue+2); 
    histoCorrectedYieldTrueLeftWide->DrawCopy("e1,same"); 
    
    TLegend* legendYield3 = new TLegend(0.15,0.03,0.66,0.19);
    legendYield3->SetTextSize(0.02); 
    legendYield3->SetFillColor(0);
    legendYield3->AddEntry(histoCorrectedYieldTrue,"corr true eff/right norm");
    legendYield3->AddEntry(histoCorrectedYieldTrueWide,"corr true eff wide int /right norm");
    legendYield3->AddEntry(histoCorrectedYieldTrueNarrow,"corr true eff narrow int /right norm");
    legendYield3->AddEntry(histoCorrectedYieldTrueLeft,Form("corr true eff /left norm"));
    legendYield3->AddEntry(histoCorrectedYieldTrueLeftWide,"corr true eff wide int /left norm");
    legendYield3->AddEntry(histoCorrectedYieldTrueLeftNarrow,"corr true eff narrow int /left norm");
    legendYield3->Draw();

    PutProcessLabelAndEnergyOnPlot(0.62, 0.95, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
    
    padCorrectedYieldRatios->cd();
    padCorrectedYieldRatios->SetTickx();
    padCorrectedYieldRatios->SetTicky();
    padCorrectedYieldRatios->SetLogy(0);
    TH1D *RatioTrue = (TH1D*) histoCorrectedYieldTrue->Clone(); 
    RatioTrue->Divide(RatioTrue,histoCorrectedYieldTrue,1.,1.,"");
    TH1D *RatioTrueWide = (TH1D*) histoCorrectedYieldTrue->Clone();
    RatioTrueWide->Divide(RatioTrueWide,histoCorrectedYieldTrueWide,1.,1.,"");
    TH1D *RatioTrueNarrow = (TH1D*) histoCorrectedYieldTrue->Clone(); 
    RatioTrueNarrow->Divide(RatioTrueNarrow,histoCorrectedYieldTrueNarrow,1.,1.,"");
    TH1D *RatioTrueLeft = (TH1D*) histoCorrectedYieldTrue->Clone();
    RatioTrueLeft->Divide(RatioTrueLeft,histoCorrectedYieldTrueLeft,1.,1.,"");
    TH1D *RatioTrueLeftWide = (TH1D*) histoCorrectedYieldTrue->Clone();  
    RatioTrueLeftWide->Divide(RatioTrueLeftWide,histoCorrectedYieldTrueLeftWide,1.,1.,"");
    TH1D *RatioTrueLeftNarrow = (TH1D*) histoCorrectedYieldTrue->Clone();
    RatioTrueLeftNarrow->Divide(RatioTrueLeftNarrow,histoCorrectedYieldTrueLeftNarrow,1.,1.,"");
    TH1D *RatioTrueMCInput = (TH1D*) histoCorrectedYieldTrue->Clone();
    RatioTrueMCInput->Divide(RatioTrueMCInput,histoMCYieldMeson,1.,1.,"");
    TH1D *RatioNormal = (TH1D*) histoCorrectedYieldTrue->Clone();  
    RatioNormal->Divide(RatioNormal,histoCorrectedYieldNorm,1.,1.,"");

    RatioTrue->SetYTitle("#frac{standard}{modified}"); 
    RatioTrue->GetYaxis()->SetRangeUser(0.8,1.23);
    RatioTrue->GetYaxis()->SetLabelSize(0.07);
    RatioTrue->GetYaxis()->SetNdivisions(505);
    RatioTrue->GetYaxis()->SetTitleSize(0.1); 
    RatioTrue->GetYaxis()->SetDecimals();
    RatioTrue->GetYaxis()->SetTitleOffset(0.55);
    RatioTrue->GetXaxis()->SetTitleSize(0.11);
    RatioTrue->GetXaxis()->SetLabelSize(0.07);
    RatioTrue->SetMarkerStyle(20);
    RatioTrue->SetMarkerSize(1.);
    RatioTrue->SetMarkerColor(kBlack);
    RatioTrue->SetLineColor(kBlack);
    RatioTrue->DrawCopy("p,e1"); 

    DrawGammaSetMarker(RatioTrue, 20, 1., kBlack, kBlack); 
    RatioTrue->DrawCopy("p,e1");  
    //right Side Normalization narrow
    DrawGammaSetMarker(RatioTrueNarrow, 24, 1., kGray+1, kGray+1); 
    RatioTrueNarrow->DrawCopy("e1,same"); 
    //right Side Normalization wide
    DrawGammaSetMarker(RatioTrueWide, 24, 1., kGray+3, kGray+3); 
    RatioTrueWide->DrawCopy("e1,same"); 
    //left Side Normalization 
    DrawGammaSetMarker(RatioTrueLeft, 20, 1., kBlue, kBlue);  
    RatioTrueLeft->DrawCopy("e1,same"); 
    //left Side Normalization narrow
    DrawGammaSetMarker(RatioTrueLeftNarrow, 24, 1., kBlue-5, kBlue-5);
    RatioTrueLeftNarrow->DrawCopy("e1,same"); 
    //left Side Normalization wide
    DrawGammaSetMarker(RatioTrueLeftWide, 24, 1., kBlue+2, kBlue+2);
    RatioTrueLeftWide->DrawCopy("e1,same"); 

    canvasCorrecftedYield->Update();
    canvasCorrecftedYield->SaveAs(Form("%s/%s_%s_CorrectedYieldTrueEff_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(),  fCutSelection.Data(), suffix.Data()));

    // **************************************************************************************
    // ************** Plot corrected yield with differnt efficiencies & MC yield ************
    // **************************** Sanity check for MC *************************************
    // **************************************************************************************
    
    if (kIsMC){
        canvasCorrecftedYield->cd();

        padCorrectedYieldHistos->cd();
        padCorrectedYieldHistos->SetLogy(); 

        DrawGammaSetMarker(histoCorrectedYieldTrue, 20, 1., kBlack, kBlack);  
        histoCorrectedYieldTrue->DrawCopy("e1");  
        DrawGammaSetMarker(histoMCYieldMeson, 24, 1., kRed+2, kRed+2);  
        histoMCYieldMeson->DrawCopy("e1,same"); 
        DrawGammaSetMarker(histoCorrectedYieldNorm, 24, 1., kGreen+2, kGreen+2); 
        histoCorrectedYieldNorm->DrawCopy("e1,same"); 
        
        cout << "here" << endl; 
        TLegend* legendYield4 = new TLegend(0.15,0.03,0.66,0.19);
        legendYield4->SetTextSize(0.02); 
        legendYield4->SetFillColor(0);
        legendYield4->SetFillStyle(0);
        legendYield4->AddEntry(histoCorrectedYieldTrue,"corr true eff");
        legendYield4->AddEntry(histoMCYieldMeson,"MC input (possibly weighted)");
        legendYield4->AddEntry(histoCorrectedYieldNorm,"normal eff");
        legendYield4->Draw();

        PutProcessLabelAndEnergyOnPlot(0.62, 0.95, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
        padCorrectedYieldRatios->cd();
        padCorrectedYieldRatios->SetTickx();
        padCorrectedYieldRatios->SetTicky();
        padCorrectedYieldRatios->SetLogy(0);
        
        DrawGammaSetMarker(RatioTrue, 20, 1., kBlack, kBlack);
        for(Int_t b = 0; b< RatioTrue->GetNbinsX(); b++){
            RatioTrue->SetBinError(b+1,histoCorrectedYieldTrue->GetBinError(b+1)/histoCorrectedYieldTrue->GetBinContent(b+1));
        }
        RatioTrue->SetFillColor(kGray+2);
        RatioTrue->SetFillStyle(1);
        
        RatioTrue->DrawCopy("p,e2");  
        DrawGammaSetMarker(RatioTrueMCInput, 24, 1., kRed+2, kRed+2);
        RatioTrueMCInput->DrawCopy("e1,same"); 
        DrawGammaSetMarker(RatioNormal, 25, 1., kGreen+2, kGreen+2); 
        RatioNormal->DrawCopy("e1,same"); 

        canvasCorrecftedYield->Update();
        canvasCorrecftedYield->SaveAs(Form("%s/%s_%s_CorrectedYield_SanityCheck_%s.%s",outputDir.Data(), nameMeson.Data(), prefix2.Data(),  fCutSelection.Data(), suffix.Data()));
    }

    delete canvasCorrecftedYield;
    delete legendYield3;

    //***********************************************************************************************
    //***************************  Secondary RAW Yield  *********************************************
    //***********************************************************************************************
//     cout << histoYieldSecMeson << endl;
//     cout << histoYieldSecFromK0SMeson << endl;
    
    if (!optDalitz && doubleAddFactorK0s >= 0 && nameMeson.Contains("Pi0") ){ //&& !kIsMC
        TCanvas* canvasRAWYieldSec = new TCanvas("canvasRAWYieldSec","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRAWYieldSec, 0.10, 0.01, 0.02, 0.08); 
        canvasRAWYieldSec->SetLogy(1);
        histoUnCorrectedYieldDrawing->GetYaxis()->SetRangeUser(1e-9,8e-1);
        DrawGammaSetMarker(histoUnCorrectedYieldDrawing, 20, 1., kBlack, kBlack);
        histoUnCorrectedYieldDrawing->Draw("e1");
        histoYieldSecMeson->Scale(1./nEvt);
        DrawGammaSetMarker(histoYieldSecMeson, 20, 1., kBlue, kBlue);
        histoYieldSecMeson->DrawCopy("same,e1");  
        histoYieldSecFromK0SMeson->Scale(1./nEvt);
        DrawGammaSetMarker(histoYieldSecFromK0SMeson, 24, 1., kCyan, kCyan);  
        histoYieldSecFromK0SMeson->DrawCopy("same,e1");  

        TLegend* legendSecRAWYield = new TLegend(0.6,0.8,0.93,0.93);
        legendSecRAWYield->SetTextSize(0.03);  
        legendSecRAWYield->SetFillColor(0);
        legendSecRAWYield->SetFillStyle(0);
        legendSecRAWYield->SetLineColor(0);
        legendSecRAWYield->SetBorderSize(0);
        legendSecRAWYield->AddEntry(histoUnCorrectedYieldDrawing,"RAW yield");
        legendSecRAWYield->AddEntry(histoYieldSecMeson,"total secondaries");
        legendSecRAWYield->AddEntry(histoYieldSecFromK0SMeson,"additional secondaries from K^{0}_{s}");
        legendSecRAWYield->Draw();

        PutProcessLabelAndEnergyOnPlot(0.62, 0.78, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
        
        canvasRAWYieldSec->Update();
        canvasRAWYieldSec->SaveAs(Form("%s/%s_%s_RAWYieldSecPt_%s.%s",outputDir.Data(),nameMeson.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
        delete canvasRAWYieldSec;
    } else if (optDalitz){
        TCanvas* canvasRAWYieldSec = new TCanvas("canvasRAWYieldSec","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRAWYieldSec, 0.10, 0.01, 0.02, 0.08); 
        canvasRAWYieldSec->SetLogy(1);

        DrawGammaSetMarker(histoUnCorrectedYieldDrawing, 20, 1., kBlack, kBlack);
        histoUnCorrectedYieldDrawing->Draw("e1");
        histoYieldGGMeson->Scale(1./nEvt);
        DrawGammaSetMarker(histoYieldGGMeson, 20, 1., kBlue, kBlue); 
        histoYieldGGMeson->DrawCopy("same,e1");
            
        TLegend* legendSecRAWYield = new TLegend(0.6,0.8,0.93,0.93);
        legendSecRAWYield->SetTextSize(0.03);  
        legendSecRAWYield->SetFillColor(0);
        legendSecRAWYield->SetBorderSize(0);
        legendSecRAWYield->AddEntry(histoUnCorrectedYieldDrawing,"RAW yield");
        legendSecRAWYield->AddEntry(histoYieldGGMeson,"Contamination from #gamma#gamma");
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
        TH1D* histoTrueTotalRecYield =          (TH1D*)fileCorrections->Get("histoYieldTrueMesonFixedWindow");
        if (histoTrueTotalRecYield){
            TH1D* histoTrueTotalRecYieldGamma =          (TH1D*)fileCorrections->Get("histoYieldTrueMesonGammaFixedWindow");
            TH1D* histoTrueTotalRecYieldConvGamma =          (TH1D*)fileCorrections->Get("histoYieldTrueMesonGammaConvGammaFixedWindow");
            
            TH1D* ratioTrueGammaDivTotal = NULL;
            TH1D* ratioTrueConvGammaDivTotal = NULL;
            if (histoTrueTotalRecYieldGamma){
                ratioTrueGammaDivTotal = (TH1D*)histoTrueTotalRecYieldGamma->Clone("ratioTrueGammaDivTotal");
                ratioTrueGammaDivTotal->Divide(ratioTrueGammaDivTotal,histoTrueTotalRecYield,1.,1.,"B");
            }
            if (histoTrueTotalRecYieldConvGamma){
                ratioTrueConvGammaDivTotal = (TH1D*)histoTrueTotalRecYieldConvGamma->Clone("ratioTrueConvGammaDivTotal");
                ratioTrueConvGammaDivTotal->Divide(ratioTrueConvGammaDivTotal,histoTrueTotalRecYield,1.,1.,"B");
            }
            
            if (histoTrueTotalRecYieldGamma){
                TCanvas* canvasFracDifferentContrib = new TCanvas("canvasFracDifferentContrib","",200,10,1350,900);  // gives the page size
                DrawGammaCanvasSettings( canvasFracDifferentContrib, 0.09, 0.02, 0.04, 0.09);
                //       canvasSecFrac->SetLogy(1); 
                            
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
                
                TLegend* legendFracDifferentContrib = new TLegend(0.6,0.8,0.94,0.93);
                legendFracDifferentContrib->SetTextSize(0.03);
                legendFracDifferentContrib->SetLineColor(0);
                legendFracDifferentContrib->SetLineWidth(0);
                legendFracDifferentContrib->SetFillColor(0);
                legendFracDifferentContrib->SetFillStyle(0);
                legendFracDifferentContrib->AddEntry(ratioTrueGammaDivTotal,"XX= #gamma_{PCM}, #gamma_{cluster}");
                legendFracDifferentContrib->AddEntry(ratioTrueConvGammaDivTotal,"XX= #gamma_{PCM}, #gamma_{conv,cluster}");
                legendFracDifferentContrib->Draw();

                PutProcessLabelAndEnergyOnPlot(0.15, 0.92, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
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
        TH1D* histoTrueTotalRecYield =          (TH1D*)fileCorrections->Get("histoYieldTrueMesonFixedWindow");
        if (histoTrueTotalRecYield){
            TH1D* histoTrueTotalRecYieldGamma =          (TH1D*)fileCorrections->Get("histoYieldTrueMesonGammaFixedWindow");
            TH1D* histoTrueTotalRecYieldConvGamma =          (TH1D*)fileCorrections->Get("histoYieldTrueMesonGammaConvGammaFixedWindow");
            TH1D* histoTrueTotalRecYieldConvGamma2 =          (TH1D*)fileCorrections->Get("histoYieldTrueMesonConvGammaConvGammaFixedWindow");
            TH1D* ratioTrueGammaDivTotal = NULL;
            TH1D* ratioTrueConvGammaDivTotal = NULL;
            TH1D* ratioTrueConvGamma2DivTotal = NULL;
            if (histoTrueTotalRecYieldGamma){
                ratioTrueGammaDivTotal = (TH1D*)histoTrueTotalRecYieldGamma->Clone("ratioTrueGammaDivTotal");
                ratioTrueGammaDivTotal->Divide(ratioTrueGammaDivTotal,histoTrueTotalRecYield,1.,1.,"B");
            }
            if (histoTrueTotalRecYieldConvGamma){
                ratioTrueConvGammaDivTotal = (TH1D*)histoTrueTotalRecYieldConvGamma->Clone("ratioTrueConvGammaDivTotal");
                ratioTrueConvGammaDivTotal->Divide(ratioTrueConvGammaDivTotal,histoTrueTotalRecYield,1.,1.,"B");
            }
            if (histoTrueTotalRecYieldConvGamma2){
                ratioTrueConvGamma2DivTotal = (TH1D*)histoTrueTotalRecYieldConvGamma2->Clone("ratioTrueConvGamma2DivTotal");
                ratioTrueConvGamma2DivTotal->Divide(ratioTrueConvGamma2DivTotal,histoTrueTotalRecYield,1.,1.,"B");
            }
            
            if (histoTrueTotalRecYieldGamma){
                TCanvas* canvasFracDifferentContrib = new TCanvas("canvasFracDifferentContrib","",200,10,1350,900);  // gives the page size
                DrawGammaCanvasSettings( canvasFracDifferentContrib, 0.09, 0.02, 0.04, 0.09);
                //       canvasSecFrac->SetLogy(1); 
                            
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
                
                TLegend* legendFracDifferentContrib = new TLegend(0.6,0.8,0.94,0.93);
                legendFracDifferentContrib->SetTextSize(0.03);
                legendFracDifferentContrib->SetLineColor(0);
                legendFracDifferentContrib->SetLineWidth(0);
                legendFracDifferentContrib->SetFillColor(0);
                legendFracDifferentContrib->SetFillStyle(0);
                legendFracDifferentContrib->AddEntry(ratioTrueGammaDivTotal,"XX= #gamma_{cluster}, #gamma_{cluster}");
                legendFracDifferentContrib->AddEntry(ratioTrueConvGammaDivTotal,"XX= #gamma_{cluster}, #gamma_{conv,cluster}");
                legendFracDifferentContrib->AddEntry(ratioTrueConvGamma2DivTotal,"XX= #gamma_{conv,cluster}, #gamma_{conv,cluster}");
                legendFracDifferentContrib->Draw();

                PutProcessLabelAndEnergyOnPlot(0.15, 0.92, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
                canvasFracDifferentContrib->Update();
                canvasFracDifferentContrib->SaveAs(Form("%s/%s_RelativeContributionsToTruePeak_%s.%s",outputDir.Data(),nameMeson.Data(),fCutSelection.Data(),suffix.Data()));
                delete canvasFracDifferentContrib;
            }
            
        }        
    }    

    
    Int_t nBinsPt =   histoCorrectedYieldTrue->GetNbinsX();
    TH1D *histoCorrectedYieldTrueMt = NULL;
    TH1D *histoCorrectedYieldTrueNarrowMt = NULL;
    TH1D *histoCorrectedYieldTrueWideMt = NULL;
    TH1D *histoCorrectedYieldTrueLeftMt = NULL;
    TH1D *histoCorrectedYieldTrueLeftNarrowMt = NULL;
    TH1D *histoCorrectedYieldTrueLeftWideMt = NULL;
    TH1D *histoRatioTrueMt = NULL;
    TH1D *histoRatioTrueWideMt = NULL;
    TH1D *histoRatioTrueNarrowMt = NULL;
    TH1D *histoRatioTrueLeftMt = NULL;
    TH1D *histoRatioTrueLeftWideMt = NULL;
    TH1D *histoRatioTrueLeftNarrowMt = NULL;
    

    if (!kDalitz){
        //***********************************************************************************************
        //***************************  correction for yield in mt bins **********************************
        //***********************************************************************************************
        

        Double_t * binsMt = NULL;
        binsMt=        new Double_t[nBinsPt+1];
        binsMt[0] =       mesonMassExpect;


        for(Int_t iPt=1;iPt<nBinsPt+1;iPt++){
            binsMt[iPt]=pow((histoCorrectedYieldTrue->GetXaxis()->GetBinUpEdge(iPt)*histoCorrectedYieldTrue->GetXaxis()->GetBinUpEdge(iPt)+mesonMassExpect*mesonMassExpect),0.5);
            cout << "recalculation pt to mt:    " << iPt<<"     pt      "<< histoCorrectedYieldTrue->GetXaxis()->GetBinUpEdge(iPt)<< "      mt    " << binsMt[iPt]<< endl;
        }

        TH1F *deltaMt =   new TH1F("deltaMt","",nBinsPt,binsMt);

        histoCorrectedYieldTrueMt = new TH1D("CorrectedYieldTrueEff_Mt","",nBinsPt,binsMt);
        histoCorrectedYieldTrueMt->SetName("CorrectedYieldTrueEff_Mt");
        histoCorrectedYieldTrueNarrowMt = new TH1D("CorrectedYieldTrueEffNarrow_Mt","",nBinsPt,binsMt);
        histoCorrectedYieldTrueNarrowMt->SetName("CorrectedYieldTrueEffNarrow_Mt");
        histoCorrectedYieldTrueWideMt = new TH1D("CorrectedYieldTrueEffWide_Mt","",nBinsPt,binsMt);
        histoCorrectedYieldTrueWideMt->SetName("CorrectedYieldTrueEffWide_Mt");

        histoCorrectedYieldTrueLeftMt = new TH1D("CorrectedYieldTrueEffLeft_Mt","",nBinsPt,binsMt);
        histoCorrectedYieldTrueLeftMt->SetName("CorrectedYieldTrueEffLeft_Mt");
        histoCorrectedYieldTrueLeftNarrowMt = new TH1D("CorrectedYieldTrueEffLeftNarrow_Mt","",nBinsPt,binsMt);
        histoCorrectedYieldTrueLeftNarrowMt->SetName("CorrectedYieldTrueEffLeftNarrow_Mt");
        histoCorrectedYieldTrueLeftWideMt = new TH1D("CorrectedYieldTrueEffLeftWide_Mt","",nBinsPt,binsMt);
        histoCorrectedYieldTrueLeftWideMt->SetName("CorrectedYieldTrueEffLeftWide_Mt");

        for(Int_t iPt=1;iPt<nBinsPt+1;iPt++){
            deltaMt->SetBinContent(iPt,binsMt[iPt]-binsMt[iPt-1]);
            deltaMt->SetBinError(iPt,0);
        }

        for(Int_t iPt=1;iPt<nBinsPt+1;iPt++){
            histoCorrectedYieldTrueMt->SetBinContent(iPt,histoCorrectedYieldTrue->GetBinContent(iPt));
            histoCorrectedYieldTrueMt->SetBinError(iPt,histoCorrectedYieldTrue->GetBinError(iPt));
            histoCorrectedYieldTrueNarrowMt->SetBinContent(iPt,histoCorrectedYieldTrueNarrow->GetBinContent(iPt));
            histoCorrectedYieldTrueNarrowMt->SetBinError(iPt,histoCorrectedYieldTrueNarrow->GetBinError(iPt));
            histoCorrectedYieldTrueWideMt->SetBinContent(iPt,histoCorrectedYieldTrueWide->GetBinContent(iPt));
            histoCorrectedYieldTrueWideMt->SetBinError(iPt,histoCorrectedYieldTrueWide->GetBinError(iPt));
            histoCorrectedYieldTrueLeftMt->SetBinContent(iPt,histoCorrectedYieldTrueLeft->GetBinContent(iPt));
            histoCorrectedYieldTrueLeftMt->SetBinError(iPt,histoCorrectedYieldTrueLeft->GetBinError(iPt));
            histoCorrectedYieldTrueLeftNarrowMt->SetBinContent(iPt,histoCorrectedYieldTrueLeftNarrow->GetBinContent(iPt));
            histoCorrectedYieldTrueLeftNarrowMt->SetBinError(iPt,histoCorrectedYieldTrueLeftNarrow->GetBinError(iPt));
            histoCorrectedYieldTrueLeftWideMt->SetBinContent(iPt,histoCorrectedYieldTrueLeftWide->GetBinContent(iPt));
            histoCorrectedYieldTrueLeftWideMt->SetBinError(iPt,histoCorrectedYieldTrueLeftWide->GetBinError(iPt));
        }

        histoRatioTrueMt = (TH1D*) histoCorrectedYieldTrueMt->Clone(); 
        histoRatioTrueMt->Divide(histoRatioTrueMt,histoCorrectedYieldTrueMt,1.,1.,"");
        histoRatioTrueWideMt = (TH1D*) histoCorrectedYieldTrueMt->Clone();
        histoRatioTrueWideMt->Divide(histoRatioTrueWideMt,histoCorrectedYieldTrueWideMt,1.,1.,"");
        histoRatioTrueNarrowMt = (TH1D*) histoCorrectedYieldTrueMt->Clone(); 
        histoRatioTrueNarrowMt->Divide(histoRatioTrueNarrowMt,histoCorrectedYieldTrueNarrowMt,1.,1.,"");
        histoRatioTrueLeftMt = (TH1D*) histoCorrectedYieldTrueMt->Clone();
        histoRatioTrueLeftMt->Divide(histoRatioTrueMt,histoCorrectedYieldTrueLeftMt,1.,1.,"");
        histoRatioTrueLeftWideMt = (TH1D*) histoCorrectedYieldTrueMt->Clone();  
        histoRatioTrueLeftWideMt->Divide(histoRatioTrueLeftWideMt,histoCorrectedYieldTrueLeftWideMt,1.,1.,"");
        histoRatioTrueLeftNarrowMt = (TH1D*) histoCorrectedYieldTrueMt->Clone();
        histoRatioTrueLeftNarrowMt->Divide(histoRatioTrueLeftNarrowMt,histoCorrectedYieldTrueLeftNarrowMt,1.,1.,"");

        //*************************************************************************************************
        //*********************** Plotting Corrected Yield in mt - bins normal eff ************************
        //*************************************************************************************************

        TCanvas* canvasCorrecftedYieldMt = new TCanvas("canvasCorrecftedYieldMt","",1350,1500);  // gives the page size
        DrawGammaCanvasSettings( canvasCorrecftedYieldMt, 0.13, 0.02, 0.02, 0.09); 
        canvasCorrecftedYieldMt->SetLogy(); 

        TPad* padCorrectedYieldHistosMt = new TPad("padCorrectedYieldHistosMt", "", 0., 0.3, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padCorrectedYieldHistosMt, 0.12, 0.02, 0.02, 0.);
        padCorrectedYieldHistosMt->Draw();

        TPad* padCorrectedYieldRatiosMt = new TPad("padCorrectedYieldRatiosMt", "", 0., 0., 1., 0.3,-1, -1, -2);
        DrawGammaPadSettings( padCorrectedYieldRatiosMt, 0.12, 0.02, 0., 0.21);
        padCorrectedYieldRatiosMt->Draw();

        padCorrectedYieldHistosMt->cd();
        padCorrectedYieldHistosMt->SetLogy();  

        
        DrawAutoGammaMesonHistos( histoCorrectedYieldTrueMt, 
                                    "", "#it{m}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{m}_{T}d#it{m}_{T}d#it{y}} (#it{c}/GeV)^{2}", 
                                    kTRUE, 3., 4e-10, kTRUE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., 15.);
        DrawGammaSetMarker(histoCorrectedYieldTrueMt, 20, 1., kBlack, kBlack);
        histoCorrectedYieldTrueMt->DrawCopy("e1,same"); 
        //right Side Normalization narrow
        DrawGammaSetMarker(histoCorrectedYieldTrueNarrowMt, 24, 1., kGray+1, kGray+1);
        histoCorrectedYieldTrueNarrowMt->DrawCopy("e1,same"); 
        //right Side Normalization wide
        DrawGammaSetMarker(histoCorrectedYieldTrueWideMt, 24, 1., kGray+3, kGray+3);
        histoCorrectedYieldTrueWideMt->DrawCopy("e1,same"); 
        //left Side Normalization
        DrawGammaSetMarker(histoCorrectedYieldTrueLeftMt, 20, 1., kBlue, kBlue); 
        histoCorrectedYieldTrueLeftMt->DrawCopy("e1,same");
        //left Side Normalization narrow
        DrawGammaSetMarker(histoCorrectedYieldTrueLeftNarrowMt, 24, 1., kBlue-5, kBlue-5);  
        histoCorrectedYieldTrueLeftNarrowMt->DrawCopy("e1,same"); 
        //left Side Normalization wide
        DrawGammaSetMarker(histoCorrectedYieldTrueLeftWideMt, 24, 1., kBlue+2, kBlue+2);  
        histoCorrectedYieldTrueLeftWideMt->DrawCopy("e1,same"); 


        TLegend* legendYield_Mt = new TLegend(0.15,0.03,0.66,0.19);
        legendYield_Mt->SetTextSize(0.02);  
        legendYield_Mt->SetFillColor(0);
        legendYield_Mt->AddEntry(histoCorrectedYieldTrueMt,"corr true /right norm");
        legendYield_Mt->AddEntry(histoCorrectedYieldTrueWideMt,"corr true wide int /right norm");
        legendYield_Mt->AddEntry(histoCorrectedYieldTrueNarrowMt,"corr true narrow int /right norm");
        legendYield_Mt->AddEntry(histoCorrectedYieldTrueLeftMt,Form("corr true  /left norm"));
        legendYield_Mt->AddEntry(histoCorrectedYieldTrueLeftWideMt,"corr true wide int /left norm");
        legendYield_Mt->AddEntry(histoCorrectedYieldTrueLeftNarrowMt,"corr true narrow int /left norm");

        legendYield_Mt->Draw();

        PutProcessLabelAndEnergyOnPlot(0.62, 0.95, 0.03, collisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
        padCorrectedYieldRatiosMt->cd();
        histoRatioTrueMt->SetYTitle("#frac{standard}{modified}");
        histoRatioTrueMt->SetXTitle("#it{m}_{T} (GeV/#it{c})");
        if(!kIsEta ) histoRatioTrueMt->GetYaxis()->SetRangeUser(0.8,1.23);  
        if(kIsEta ) histoRatioTrueMt->GetYaxis()->SetRangeUser(0.4,1.63);
//         histoRatioTrueMt->GetXaxis()->SetRangeUser(0.,15.);
        histoRatioTrueMt->GetYaxis()->SetNdivisions(505);
        //histoRatioTrueMt->GetXaxis()->SetRangeUser(0.,1.5);
        histoRatioTrueMt->GetYaxis()->SetLabelSize(0.08);
        histoRatioTrueMt->GetYaxis()->SetTitleSize(0.1);
        histoRatioTrueMt->GetYaxis()->SetDecimals();
        histoRatioTrueMt->GetYaxis()->SetTitleOffset(0.55);
        histoRatioTrueMt->GetXaxis()->SetTitleSize(0.11);  
        histoRatioTrueMt->GetXaxis()->SetLabelSize(0.08);

        DrawGammaSetMarker(histoRatioTrueMt, 20, 1., kBlack, kBlack);
        histoRatioTrueMt->DrawCopy("e1");
        //right Side Normalization narrow
        DrawGammaSetMarker(histoRatioTrueNarrowMt, 24, 1., kGray+1, kGray+1);
        histoRatioTrueNarrowMt->DrawCopy("e1,same"); 
        //right Side Normalization wide
        DrawGammaSetMarker(histoRatioTrueWideMt, 24, 1., kGray+3, kGray+3);
        histoRatioTrueWideMt->DrawCopy("e1,same"); 
        //left Side Normalization
        DrawGammaSetMarker(histoRatioTrueLeftMt, 20, 1., kBlue, kBlue); 
        histoRatioTrueLeftMt->DrawCopy("e1,same");
        //left Side Normalization narrow
        DrawGammaSetMarker(histoRatioTrueLeftNarrowMt, 24, 1., kBlue-5, kBlue-5);  
        histoRatioTrueLeftNarrowMt->DrawCopy("e1,same"); 
        //left Side Normalization wide
        DrawGammaSetMarker(histoRatioTrueLeftWideMt, 24, 1., kBlue+2, kBlue+2);  
        histoRatioTrueLeftWideMt->DrawCopy("e1,same"); 
        DrawGammaLines(0., maxPtMeson,1., 1.,0.1);

        canvasCorrecftedYieldMt->Update();

        canvasCorrecftedYieldMt->SaveAs(Form("%s/%s_%s_CorrectedYield_Mtspectra_TrueEff_%s.%s", outputDir.Data(), nameMeson.Data() ,prefix2.Data() , fCutSelection.Data(), suffix.Data()));
        delete canvasCorrecftedYieldMt;
        delete legendYield_Mt;
            
    }
    
    //*************************************************************************************************
    //******************** Output of the systematic Error due to Signal extraction ********************
    //*************************************************************************************************
    Double_t  binsXCenter[50];
    Double_t  binsXWidth[50];
    //    Double_t binYValue[50];
    //    binsXCenter[0] =     0;
    binsXWidth[0]=       0.;
    //    binYValue[0]=        0.;
    for (Int_t i = 1; i < nBinsPt +1; i++){
        binsXCenter[i] =  histoCorrectedYieldTrue->GetBinCenter(i);
        binsXWidth[i]=    histoCorrectedYieldTrue->GetBinWidth(i)/2.;
    }


    SysErrorConversion sysErrNormal[50];
    SysErrorConversion sysErrLeft[50];
    SysErrorConversion sysErrWide[50];
    SysErrorConversion sysErrNarrow[50];
    SysErrorConversion sysErrLeftWide[50];
    SysErrorConversion sysErrLeftNarrow[50];

    for (Int_t i = 1; i < nBinsPt +1; i++){
        if (mode == 9 || mode == 0 || mode == 1){
//          binYValue[i] = histoCorrectedYieldTrue->GetBinContent(i);
            sysErrNormal[i].value = histoCorrectedYieldTrue->GetBinContent(i);
            sysErrNormal[i].error = histoCorrectedYieldTrue->GetBinError(i);
            sysErrLeft[i].value = histoCorrectedYieldTrueLeft->GetBinContent(i);
            sysErrLeft[i].error = histoCorrectedYieldTrueLeft->GetBinError(i);
            sysErrNarrow[i].value = histoCorrectedYieldTrueNarrow->GetBinContent(i);
            sysErrNarrow[i].error = histoCorrectedYieldTrueNarrow->GetBinError(i);
            sysErrLeftNarrow[i].value = histoCorrectedYieldTrueLeftNarrow->GetBinContent(i);
            sysErrLeftNarrow[i].error = histoCorrectedYieldTrueLeftNarrow->GetBinError(i);
            sysErrWide[i].value = histoCorrectedYieldTrueWide->GetBinContent(i);
            sysErrWide[i].error = histoCorrectedYieldTrueWide->GetBinError(i);
            sysErrLeftWide[i].value = histoCorrectedYieldTrueLeftWide->GetBinContent(i);
            sysErrLeftWide[i].error = histoCorrectedYieldTrueLeftWide->GetBinError(i); 
        } else {
            sysErrNormal[i].value = histoCorrectedYieldNorm->GetBinContent(i);
            sysErrNormal[i].error = histoCorrectedYieldNorm->GetBinError(i);
            sysErrLeft[i].value = histoCorrectedYieldNormLeft->GetBinContent(i);
            sysErrLeft[i].error = histoCorrectedYieldNormLeft->GetBinError(i);
            sysErrNarrow[i].value = histoCorrectedYieldNormNarrow->GetBinContent(i);
            sysErrNarrow[i].error = histoCorrectedYieldNormNarrow->GetBinError(i);
            sysErrLeftNarrow[i].value = histoCorrectedYieldNormLeftNarrow->GetBinContent(i);
            sysErrLeftNarrow[i].error = histoCorrectedYieldNormLeftNarrow->GetBinError(i);
            sysErrWide[i].value = histoCorrectedYieldNormWide->GetBinContent(i);
            sysErrWide[i].error = histoCorrectedYieldNormWide->GetBinError(i);
            sysErrLeftWide[i].value = histoCorrectedYieldNormLeftWide->GetBinContent(i);
            sysErrLeftWide[i].error = histoCorrectedYieldNormLeftWide->GetBinError(i); 
        }    
        
    }
    Double_t differenceLeft[50];
    Double_t differenceLeftWide[50];
    Double_t differenceLeftNarrow[50];
    Double_t differenceWide[50];
    Double_t differenceNarrow[50];
    Double_t largestDifferenceNeg[50];
    Double_t largestDifferencePos[50];
    Double_t differenceLeftError[50];
    Double_t differenceLeftWideError[50];
    Double_t differenceLeftNarrowError[50];
    Double_t differenceWideError[50];
    Double_t differenceNarrowError[50];
    Double_t largestDifferenceNegError[50];
    Double_t largestDifferencePosError[50];
    Double_t relDifferenceLeft[50];
    Double_t relDifferenceLeftWide[50];
    Double_t relDifferenceLeftNarrow[50];
    Double_t relDifferenceWide[50];
    Double_t relDifferenceNarrow[50];
    Double_t relLargestDifferenceNeg[50];
    Double_t relLargestDifferencePos[50];
    Double_t relDifferenceLeftError[50];
    Double_t relDifferenceLeftWideError[50];
    Double_t relDifferenceLeftNarrowError[50];
    Double_t relDifferenceWideError[50];
    Double_t relDifferenceNarrowError[50];
    Double_t relLargestDifferenceNegError[50];
    Double_t relLargestDifferencePosError[50];
    for (Int_t i = 1; i < nBinsPt +1; i++){
        largestDifferenceNeg[i] =     0;
        largestDifferencePos[i] =     0;
        largestDifferenceNegError[i] =   0;
        largestDifferencePosError[i] =   0;
        relLargestDifferenceNeg[i] =     0;
        relLargestDifferencePos[i] =     0;
        relLargestDifferenceNegError[i] =   0;
        relLargestDifferencePosError[i] =   0;
        //Calculate differences
        differenceLeft[i] =        sysErrLeft[i].value - sysErrNormal[i].value;
        differenceLeftError[i] =      TMath::Sqrt(TMath::Abs(TMath::Power(sysErrLeft[i].error,2)-TMath::Power(sysErrNormal[i].error,2)));
        differenceLeftNarrow[i] =     sysErrLeftNarrow[i].value - sysErrNormal[i].value;
        differenceLeftNarrowError[i] =   TMath::Sqrt(TMath::Abs(TMath::Power(sysErrLeftNarrow[i].error,2)-TMath::Power(sysErrNormal[i].error,2)));
        differenceLeftWide[i] =          sysErrLeftWide[i].value - sysErrNormal[i].value;
        differenceLeftWideError[i] =     TMath::Sqrt(TMath::Abs(TMath::Power(sysErrLeftWide[i].error,2)-TMath::Power(sysErrNormal[i].error,2)));
        differenceNarrow[i] =         sysErrNarrow[i].value - sysErrNormal[i].value;
        differenceNarrowError[i] =       TMath::Sqrt(TMath::Abs(TMath::Power(sysErrNarrow[i].error,2)-TMath::Power(sysErrNormal[i].error,2)));
        differenceWide[i] =        sysErrWide[i].value - sysErrNormal[i].value;
        differenceWideError[i] =      TMath::Sqrt(TMath::Abs(TMath::Power(sysErrWide[i].error,2)-TMath::Power(sysErrNormal[i].error,2)));

        if (sysErrNormal[i].value != 0){
            relDifferenceLeft[i] =     (differenceLeft[i]/sysErrNormal[i].value)*100.;
            relDifferenceLeftError[i] =   (differenceLeftError[i]/sysErrNormal[i].value)*100.;
            relDifferenceLeftNarrow[i] =  (differenceLeftNarrow[i]/sysErrNormal[i].value)*100.;
            relDifferenceLeftNarrowError[i] = (differenceLeftNarrowError[i]/sysErrNormal[i].value)*100.;
            relDifferenceLeftWide[i] =    (differenceLeftWide[i]/sysErrNormal[i].value)*100.;
            relDifferenceLeftWideError[i] = (differenceLeftWideError[i]/sysErrNormal[i].value)*100.;
            relDifferenceNarrow[i] =   (differenceNarrow[i]/sysErrNormal[i].value)*100.;
            relDifferenceNarrowError[i] = (differenceNarrowError[i]/sysErrNormal[i].value)*100.;
            relDifferenceWide[i] =     (differenceWide[i]/sysErrNormal[i].value)*100.;
            relDifferenceWideError[i] =   (differenceWideError[i]/sysErrNormal[i].value)*100.;
        } else {
            relDifferenceLeft[i] =     0.;
            relDifferenceLeftError[i] =   0.;
            relDifferenceLeftNarrow[i] =  0.;
            relDifferenceLeftNarrowError[i] = 0.;
            relDifferenceLeftWide[i] =    0.;
            relDifferenceLeftWideError[i] = 0.;
            relDifferenceNarrow[i] =      0.;
            relDifferenceNarrowError[i] = 0.;
            relDifferenceWide[i] =     0.;
            relDifferenceWideError[i] =   0.;
        }

        //Find biggest Deviation
        if (TMath::Abs(relDifferenceLeft[i]) < 75. ){
            if(differenceLeft[i] < 0){
                largestDifferenceNeg[i] =     differenceLeft[i];
                largestDifferenceNegError[i] =   differenceLeftError[i];
                relLargestDifferenceNeg[i] =     relDifferenceLeft[i];
                relLargestDifferenceNegError[i] =   relDifferenceLeftError[i];
            }else{   
                largestDifferencePos[i] =     differenceLeft[i]; 
                largestDifferencePosError[i] =   differenceLeftError[i]; 
                relLargestDifferencePos[i] =     relDifferenceLeft[i]; 
                relLargestDifferencePosError[i] =   relDifferenceLeftError[i]; 
            }  
        }
        if (TMath::Abs(relDifferenceNarrow[i]) < 75.){
            if(differenceNarrow[i] < 0){
                if(differenceNarrow[i] < largestDifferenceNeg[i]){
                largestDifferenceNeg[i] =     differenceNarrow[i]; 
                largestDifferenceNegError[i] =   differenceNarrowError[i];  
                relLargestDifferenceNeg[i] =     relDifferenceNarrow[i]; 
                relLargestDifferenceNegError[i] =   relDifferenceNarrowError[i];  
                }
            } else {
                if(differenceNarrow[i] > largestDifferencePos[i]){
                largestDifferencePos[i] =     differenceNarrow[i];
                largestDifferencePosError[i] =   differenceNarrowError[i];
                relLargestDifferencePos[i] =     relDifferenceNarrow[i];
                relLargestDifferencePosError[i] =   relDifferenceNarrowError[i];
                }
            }  
        }
        if (TMath::Abs(relDifferenceWide[i]) < 75.){
            if(differenceWide[i] < 0){
                if(differenceWide[i] < largestDifferenceNeg[i]){
                largestDifferenceNeg[i] =     differenceWide[i];
                largestDifferenceNegError[i] =   differenceWideError[i];
                relLargestDifferenceNeg[i] =     relDifferenceWide[i];
                relLargestDifferenceNegError[i] =   relDifferenceWideError[i];
                }
            } else {
                if(differenceWide[i] > largestDifferencePos[i]) {
                largestDifferencePos[i] =     differenceWide[i];
                largestDifferencePosError[i] =   differenceWideError[i]; 
                relLargestDifferencePos[i] =     relDifferenceWide[i];
                relLargestDifferencePosError[i] =   relDifferenceWideError[i]; 
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
        fileSysErrDat << i << "\t" << sysErrNormal[i].value << "\t" << sysErrNormal[i].error << endl;
    }
    fileSysErrDat << endl;
    fileSysErrDat << "Bin" << "\t" << "Narrow value" << "\t" << "Narrow error" << "\t" << "Diff to Norm" << endl;
    for(Int_t i = 1; i < (nBinsPt +1); i++){
        fileSysErrDat << i << "\t" << sysErrNarrow[i].value << "\t" << sysErrNarrow[i].error << "\t" << differenceNarrow[i] << "\t" << differenceNarrowError[i] << "\t" << relDifferenceNarrow[i] << "\t" << relDifferenceNarrowError[i] <<endl;
    }
    fileSysErrDat << endl;
    fileSysErrDat << "Bin" << "\t" << "Wide value" << "\t" << "Wide error" << "\t" << "Diff to Norm" << endl;
    for(Int_t i = 1; i < (nBinsPt +1); i++){
        fileSysErrDat << i << "\t" << sysErrWide[i].value << "\t" << sysErrWide[i].error << "\t" << differenceWide[i] << "\t" << differenceWideError[i] << "\t" << relDifferenceWide[i] << "\t" << relDifferenceWideError[i] <<endl;
    }
    fileSysErrDat << endl;
    fileSysErrDat << "Bin" << "\t" << "Left value" << "\t" << "Left error" << "\t" << "Diff to Norm" << endl;
    for(Int_t i = 1; i < (nBinsPt +1); i++){
        fileSysErrDat << i << "\t" << sysErrLeft[i].value << "\t" << sysErrLeft[i].error << "\t" << differenceLeft[i] << "\t" << differenceLeftError[i] << "\t" << relDifferenceLeft[i] << "\t" << relDifferenceLeftError[i] <<endl;
    }
    fileSysErrDat << endl;
    fileSysErrDat << "Bin" << "\t" << "Left Narrow value" << "\t" << "Left Narrow error" << "\t" << "Diff to Norm" << endl;
    for(Int_t i = 1; i < (nBinsPt +1); i++){
        fileSysErrDat << i << "\t" << sysErrLeftNarrow[i].value << "\t" << sysErrLeftNarrow[i].error << "\t" << differenceLeftNarrow[i] <<  "\t" << differenceLeftNarrowError[i] << "\t" << relDifferenceLeftNarrow[i] << "\t" << relDifferenceLeftNarrowError[i] <<endl;
    }
    fileSysErrDat << endl;
    fileSysErrDat << "Bin" << "\t" << "Left Wide value" << "\t" << "Left Wide error" << "\t" << "Diff to Norm" << endl;
    for(Int_t i = 1; i < (nBinsPt +1); i++){
        fileSysErrDat << i << "\t" << sysErrLeftWide[i].value << "\t" << sysErrLeftWide[i].error << "\t" << differenceLeftWide[i] << "\t" << differenceLeftWideError[i] << "\t" << relDifferenceLeftWide[i] << "\t" << relDifferenceLeftWideError[i] <<endl;
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

    Double_t relBGEstimate[50];
    Double_t relBGEstimateError[50];
    for (Int_t i = 1; i < nBinsPt +1; i++){
            relBGEstimateError[i] = 0.;
            if ( TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateCatC->GetBinContent(i)) > TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateCatD->GetBinContent(i)) && TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateA->GetBinContent(i)) > TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateA->GetBinContent(i))){
                relBGEstimate[i] = TMath::Abs(histoBGEstimateCatA->GetBinContent(i)- histoBGEstimateCatC->GetBinContent(i)) *100;
            } else if ( TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateCatD->GetBinContent(i)) > TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateCatC->GetBinContent(i)) && TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateA->GetBinContent(i)) > TMath::Abs(histoBGEstimateCatA->GetBinContent(i)-histoBGEstimateA->GetBinContent(i))) {
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
        if (histoTrueEffiPt)                histoTrueEffiPt->Write("TrueMesonEffiPt");
        if (histoCompleteCorr)              histoCompleteCorr->Write("EffiTimesAcceptanceTimesDeltaY");
    correctedOutput2->Write();
    correctedOutput2->Close();
    
    // ********************************************************************************************************************************
    // ****************************** Write file with all further needed histograms ***************************************************
    // ********************************************************************************************************************************    
    const char* nameOutput = Form("%s/%s/%s_%s_GammaConvV1%sCorrection%s_%s.root",fCutSelection.Data(),optionEnergy.Data(),nameMeson.Data(),prefix2.Data(),fDalitz.Data(),optionPeriod.Data(),fCutSelection.Data());
    TFile* correctedOutput = new TFile(nameOutput,"RECREATE");  

    if (histoCorrectedYieldNorm)            histoCorrectedYieldNorm->Write();
    if (SystErrGraphPos)                    SystErrGraphPos->Write(Form("%s_SystErrorRelPos_YieldExtraction_%s",nameMeson.Data(),centralityString2.Data()),TObject::kOverwrite);
    if (SystErrGraphNeg)                    SystErrGraphNeg->Write(Form("%s_SystErrorRelNeg_YieldExtraction_%s",nameMeson.Data(),centralityString2.Data()),TObject::kOverwrite);
    if (SystErrGraphBGEstimate)             SystErrGraphBGEstimate->Write(Form("%s_SystErrorRel_BGEstimate_%s",nameMeson.Data(),centralityString2.Data()),TObject::kOverwrite);
    if (histoBGEstimateCatA)                histoBGEstimateCatA->Write("BGEstimateFromPileup");
    if (histoCorrectionFactorsHistvsPtCatA) histoCorrectionFactorsHistvsPtCatA->Write("PileupContamination");
    if (histoCorrectedYieldTrue)            histoCorrectedYieldTrue->Write();
    if (histoCorrectedYieldTrueNarrow)      histoCorrectedYieldTrueNarrow->Write();
    if (histoCorrectedYieldTrueWide)        histoCorrectedYieldTrueWide->Write();
    if (histoCorrectedYieldFixed)           histoCorrectedYieldFixed->Write();
    if (histoCorrectedYieldNarrowFixed)     histoCorrectedYieldNarrowFixed->Write();
    if (histoCorrectedYieldWideFixed)       histoCorrectedYieldWideFixed->Write();
    if (histoCorrectedYieldTrueFixed)       histoCorrectedYieldTrueFixed->Write();
    if (histoCorrectedYieldTrueNarrowFixed) histoCorrectedYieldTrueNarrowFixed->Write();
    if (histoCorrectedYieldTrueWideFixed)   histoCorrectedYieldTrueWideFixed->Write();
    if (histoCorrectedYieldTrueLeft)        histoCorrectedYieldTrueLeft->Write();
    if (histoCorrectedYieldTrueLeftWide)    histoCorrectedYieldTrueLeftWide->Write();
    if (histoCorrectedYieldTrueLeftNarrow)  histoCorrectedYieldTrueLeftNarrow->Write();
    
    if (histoYieldSecMeson)                 histoYieldSecMeson->Write();
    if (histoYieldSecFromK0SMeson)          histoYieldSecFromK0SMeson->Write();
    
    
    if (!kDalitz){
        if (histoCorrectedYieldTrueMt)              histoCorrectedYieldTrueMt->Write();
        if (histoCorrectedYieldTrueNarrowMt)        histoCorrectedYieldTrueNarrowMt->Write();
        if (histoCorrectedYieldTrueWideMt)          histoCorrectedYieldTrueWideMt->Write();
        if (histoCorrectedYieldTrueLeftMt)          histoCorrectedYieldTrueLeftMt->Write();
        if (histoCorrectedYieldTrueLeftWideMt)      histoCorrectedYieldTrueLeftWideMt->Write();
        if (histoCorrectedYieldTrueLeftNarrowMt)    histoCorrectedYieldTrueLeftNarrowMt->Write();
    }    
    if (histoUnCorrectedYield)              histoUnCorrectedYield->Write();
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
    if (histoTrueEffiPt)                    histoTrueEffiPt->Write("TrueMesonEffiPt");
    if (histoEffiPt)                        histoEffiPt->Write("MesonEffiPt");
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
    
    if (histoEffiNarrowPt || histoEffiWidePt || histoEffiLeftPt || histoEffiLeftNarrowPt || histoEffiLeftWidePt){}
}
