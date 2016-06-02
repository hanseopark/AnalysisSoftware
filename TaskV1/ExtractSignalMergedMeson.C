// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

#include <stdlib.h>
#include <iostream>
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
#include "TProfile2D.h"
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
#include "TMath.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "../CommonHeaders/PlottingMeson.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "ExtractSignalMergedMeson.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "../CommonHeaders/ExtractSignalPlotting.h"
#include "THnSparse.h"

                                                                  
//****************************************************************************
//********* Main function for extraction of signal for merged mesons *********
//****************************************************************************
void ExtractSignalMergedMeson(  TString meson                   = "", 
                                TString file                    = "", 
                                TString cutSelection            = "", 
                                TString suffix                  = "", 
                                Bool_t  optionMC                = kFALSE, 
                                TString optionEnergy            = "", 
                                TString optionCrystalBall       = "", 
                                TString optionPeriod            = "", 
                                TString optionAdvancedMesonQA   = "",
                                Int_t   numberOfBins            = 30,
                                Int_t   mode                    = 10,
                                Bool_t  kJJGammaTriggMC         = 0
                                
                             ) {
    gROOT->Reset();

    //*****************************************************************************************************************
    //******************************* Set global variables ************************************************************
    //*****************************************************************************************************************
    fMode                           = mode;
    // mode:   10 // new output merged EMCal
    //         11 // new output merged PHOS

    if(optionAdvancedMesonQA.Contains("AdvancedMesonQA")){
        fAdvancedMesonQA            = kTRUE;
    }
    
    if (kJJGammaTriggMC){
        fAdditionalName             = "JJGammaTrigg";
    }
        
    fAdditionalLabels               = kTRUE;    
        
    fCutSelection                   = cutSelection;
    TString fCutSelectionRead       = cutSelection;
    TString dummyString             = "";
    ReturnSeparatedCutNumberAdvanced( cutSelection, fEventCutSelection, fClusterCutSelection, fClusterMergedCutSelection, dummyString, fMesonCutSelection, mode);

    StyleSettingsThesis(suffix);
    SetPlotStyle();
    
    fEnergyFlag                     = optionEnergy;
    fPrefix                         = meson;
    fPeriodFlag                     = optionPeriod;

    TString outputDir               = Form("%s/%s/%s/ExtractSignalMergedMeson",cutSelection.Data(),optionEnergy.Data(),suffix.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    
    cout<<"Pictures are saved as "<< suffix.Data()<< endl;
    fdate                           = ReturnDateString();
    
    //****************************** Specification of collision system ************************************************
    TString textProcess             = ReturnMesonString (fPrefix);
    if(textProcess.CompareTo("") == 0 ){
        cout << "Meson unknown" << endl;
        return ;
    }
    
    fTextMeasurement                = Form("%s #rightarrow #gamma#gamma", textProcess.Data());
    fCollisionSystem                = ReturnFullCollisionsSystem(fEnergyFlag);
    fDecayChannel                   = "#gamma#gamma";
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    fDetectionProcess               = ReturnFullTextReconstructionProcess(fMode,0,textProcess.Data());
    fDetectionProcessPtBins         = ReturnFullTextReconstructionProcess(fMode,0,textProcess.Data(), fClusterMergedCutSelection);

    fNLMString                      = "";
    fNLMmin                         = ReturnClusterNLM(fClusterMergedCutSelection);
    if (ReturnClusterNLM(fClusterMergedCutSelection) == 1)
        fNLMString                  = Form("%i local maximum", ReturnClusterNLM(fClusterMergedCutSelection));        
    else 
        fNLMString                  = Form("%i local maxima", ReturnClusterNLM(fClusterMergedCutSelection));
    
    Double_t startXLabelInvMass     = 0.12;
    if (ReturnClusterNLM(fClusterMergedCutSelection) == 2)
        startXLabelInvMass          = 0.55;
    
    //****************************** Choice of Fitting procedure ******************************************************
    if(optionCrystalBall.CompareTo("CrystalBall") == 0){// means we want to plot values for the pi0
        fCrysFitting                = 1;
        cout << "CrystalBall fit chosen ..." << endl;
    } else   {
        fCrysFitting                = 0;
        cout << "Gaussian fit chosen ..." << endl;
    }
    
    if(cutSelection.Length() == 0){
        cout<<"ERROR: Cut selection is not set, please do!"<<endl;
        return;
    }

    //***************************** Specification Data/MC ************************************************************
    if( optionMC ){
        fIsMC                       = 1;
        fIsMCGammaTrig              = 0;
        fPrefix2                    = "MC";
        if ( (file.Contains("LHC15a3b") || file.Contains("LHC15g1b")) && mode == 10){
            fIsMCGammaTrig          = 1;
            cout << "WARNING: running on JetJet MC trigger on EMC photons, need to correct for acceptance factor" << endl;
        }    
    } else {
        fIsMC                       = 0;
        fPrefix2                    = "data";
    }
    
    //**************************** Determine Centrality *************************************************************
    fCentralityString                = GetCentralityString(fEventCutSelection);
    if (fCentralityString.CompareTo("pp")==0){
        fTextCent                   = "MinBias";
    } else {
        fTextCent                   = Form("%s central", fCentralityString.Data());
    }
    if (fCentralityString.CompareTo("pp")!=0 && !fCentralityString.Contains("0-100%") ){
        fCollisionSystem            = Form("%s %s", fCentralityString.Data(), fCollisionSystem.Data());
    }
    
    //***************************** Initialization of variables according to meson type ******************************
    if(meson.CompareTo("Pi0") == 0){
        Initialize("Pi0",numberOfBins);
    } else if (meson.CompareTo("Eta") == 0) {
        Initialize("Eta",numberOfBins);
    } else if(meson.CompareTo("Pi0EtaBinning") == 0) {
        Initialize("Pi0EtaBinning",numberOfBins);
    } else   {
        cout<<"ERROR: First argument in the ExtractSignalMergedMeson(....) has to be either Pi0 or Eta or Pi0EtaBinning  or EtaPrim"<<endl;
        return;
    }
            
    //************************* Start of Main routine ***************************************************************
    const char* fFileErrLogDatname  = Form("%s/%s/%s_%s_FileErrLog%s_%s.dat", cutSelection.Data(), fEnergyFlag.Data(), fPrefix.Data(), fPrefix2.Data(), 
                                           fPeriodFlag.Data(), fCutSelectionRead.Data());
    fFileErrLog.open(fFileErrLogDatname, ios::out);

    TFile f(file.Data());
    
    TString nameMainDir             = "";
    if (mode == 10 || mode == 11){
        nameMainDir                 = "GammaCaloMerged";
    } else {
        cout << "ERROR: Wrong mode aborting here!" << endl;
        return;
    }
    
    TList *TopDir                   = (TList*)f.Get(nameMainDir.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }
    
    TList *HistosGammaConversion    = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if(HistosGammaConversion == NULL){
        cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in File"<<endl;
        return;
    }

    TList *ESDContainer                 = (TList*) HistosGammaConversion->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));
    fHistoClustersPt                    = (TH1D*)ESDContainer->FindObject("ClusGamma_Pt");
    if (fIsMC){
        fHistoClustersOverlapHeadersPt  = (TH1D*)ESDContainer->FindObject("ClusGammaOverlapHeaders_Pt");
    }
    fHistoClustersMergedPtM02           = (TH2F*)ESDContainer->FindObject("ClusMerged_Pt_M02");
    fHistoClustersMergedPtM02AccMeson   = (TH2F*)ESDContainer->FindObject("ClusMerged_Pt_M02_AcceptedMeson");
    fHistoClusterCandidates             = (TH1D*)ESDContainer->FindObject("GammaCandidates");
    fHistoClusterMergedCandidates       = (TH1D*)ESDContainer->FindObject("MergedCandidates");
    fHistoTrackVsClusterCandidates      = (TH2F*)ESDContainer->FindObject("GoodESDTracksVsGammaCandidates");
    fHistoSPDtrackletvsSPDclusters      = (TH2F*)ESDContainer->FindObject("SPD tracklets vs SPD clusters");
    
    fHistoClusterNCellsPt               = (TH2F*)ESDContainer->FindObject("Clus_NCells_Pt");
    fHistoClusterMergedNCellsPt         = (TH2F*)ESDContainer->FindObject("ClusMerged_NCells_Pt");
    fHistoClusterMergedNCellsArClPt     = (TH2F*)ESDContainer->FindObject("ClusMerged_NCellsAroundClus_Pt");
    fHistoClusterMergedNCellsArAInclPt  = (TH2F*)ESDContainer->FindObject("ClusMerged_NCellsAroundAndInClus_Pt");
    fHistoClusterMergedEAroundClE       = (TH2F*)ESDContainer->FindObject("ClusMerged_EAroundClus_E");
    if (fIsMC){
        fHistoClusterMergedNPartPt      = (TH2F*)ESDContainer->FindObject("ClusMerged_NPart_Pt");
    }
    if (fHistoClusterNCellsPt && fHistoClusterMergedNCellsPt && fHistoClusterMergedNCellsArClPt && fHistoClusterMergedNCellsArAInclPt && fHistoClusterMergedEAroundClE ){
        fAdvancedClusterQA              = kTRUE;
        cout << "**************************************************************" << endl;
        cout << "INFO: Enabeling QA for clusters" << endl;
        cout << "**************************************************************" << endl;
    }
        
    fNumberOfGoodESDTracks              = (TH1D*)ESDContainer->FindObject("GoodESDTracks");
    fEventQuality                       = (TH1D*)ESDContainer->FindObject("NEvents");
        
    TString rapidityRange;
    fYMaxMeson                          = ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
        
    TString ObjectNameESD               = "ESD_Mother_InvMass_Pt";
    
    cout << "line " << __LINE__ << endl;
    fHistoInvMassVsPt                   = (TH2F*)ESDContainer->FindObject(ObjectNameESD.Data());
    fHistoInvMassVsPt->Sumw2();
    
    const char* FileDataLogname         = Form("%s/%s/%s_%s_EffiCheck_RAWDATA%s_%s.dat", cutSelection.Data(), fEnergyFlag.Data(), fPrefix.Data(), fPrefix2.Data(), fPeriodFlag.Data(),
                                               fCutSelectionRead.Data());
    fFileDataLog.open(FileDataLogname, ios::out);
    
    
    fMesonMassExpect                    = TDatabasePDG::Instance()->GetParticle(fMesonId)->Mass();
    if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 || fEnergyFlag.CompareTo("pPb_5.023TeV") == 0){
        fNEvents                        = fEventQuality->GetBinContent(1);
    } else {
        fNEvents                        =  GetNEvents(fEventQuality);
    }

    FillDataHistosArray( fHistoInvMassVsPt, fHistoClustersMergedPtM02AccMeson );

    if ( fIsMC ){
        SetCorrectMCHistogrammNames(meson);
        
        TList *TrueContainer                    = (TList*) HistosGammaConversion->FindObject(Form("%s True histograms",fCutSelectionRead.Data()));
        fHistoTrueClustersMergedPtM02           = (TH2F*) TrueContainer->FindObject(fObjectNameTrueMergedM02.Data());
        if (fHistoTrueClustersMergedPtM02) fHistoTrueClustersMergedPtM02->Sumw2();
        fHistoTrueClustersMergedInvMassPt       = (TH2F*) TrueContainer->FindObject(fObjectNameTrueMergedInvMass.Data());
        if (fHistoTrueClustersMergedInvMassPt) fHistoTrueClustersMergedInvMassPt->Sumw2();
        fHistoTrueClustersPi0PtM02              = (TH2F*) TrueContainer->FindObject(fObjectNameTrueFromPi0M02.Data());
        if (fHistoTrueClustersPi0PtM02) fHistoTrueClustersPi0PtM02->Sumw2();
        fHistoTrueClustersPi0InvMassPt          = (TH2F*) TrueContainer->FindObject(fObjectNameTrueFromPi0InvMass.Data());
        if (fHistoTrueClustersPi0InvMassPt) fHistoTrueClustersPi0InvMassPt->Sumw2();
        fHistoTrueClustersEtaPtM02              = (TH2F*) TrueContainer->FindObject(fObjectNameTrueFromEtaM02.Data());
        if (fHistoTrueClustersEtaPtM02) fHistoTrueClustersEtaPtM02->Sumw2();
        fHistoTrueClustersEtaInvMassPt          = (TH2F*) TrueContainer->FindObject(fObjectNameTrueFromEtaInvMass.Data());
        if (fHistoTrueClustersEtaInvMassPt) fHistoTrueClustersEtaInvMassPt->Sumw2();
        fHistoTrueClustersBGPtM02               = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusBGM02.Data());
        if (fHistoTrueClustersBGPtM02) fHistoTrueClustersBGPtM02->Sumw2();
        fHistoTrueClustersBGInvMassPt           = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusBGInvMass.Data());
        if (fHistoTrueClustersBGInvMassPt) fHistoTrueClustersBGInvMassPt->Sumw2();
        fHistoTrueClustersGammaPtM02            = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusGammaM02.Data());
        if (fHistoTrueClustersGammaPtM02) fHistoTrueClustersGammaPtM02->Sumw2();
        fHistoTrueClustersGammaInvMassPt        = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusGammaInvMass.Data());
        if (fHistoTrueClustersGammaInvMassPt) fHistoTrueClustersGammaInvMassPt->Sumw2();
        fHistoTrueClustersElectronPtM02         = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusElectronM02.Data());
        if (fHistoTrueClustersElectronPtM02) fHistoTrueClustersElectronPtM02->Sumw2();
        fHistoTrueClustersElectronInvMassPt     = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusElectronInvMass.Data());
        if (fHistoTrueClustersElectronInvMassPt) fHistoTrueClustersElectronInvMassPt->Sumw2();
        
        
        fHistoTrueClusPartConvMergedPtM02       = (TH2F*) TrueContainer->FindObject(fObjectNameTrueMergedPartConvM02.Data());
        if (fHistoTrueClusPartConvMergedPtM02) fHistoTrueClusPartConvMergedPtM02->Sumw2();
        fHistoTrueClusPartConvMergedInvMassPt   = (TH2F*) TrueContainer->FindObject(fObjectNameTrueMergedPartConvInvMass.Data());
        if (fHistoTrueClusPartConvMergedInvMassPt) fHistoTrueClusPartConvMergedInvMassPt->Sumw2();
        fHistoTrueClusPartConvPi0PtM02          = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusPartConvFromPi0M02.Data());
        if (fHistoTrueClusPartConvPi0PtM02) fHistoTrueClusPartConvPi0PtM02->Sumw2();
        fHistoTrueClusPartConvPi0InvMassPt      = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusPartConvFromPi0InvMass.Data());
        if (fHistoTrueClusPartConvPi0InvMassPt) fHistoTrueClusPartConvPi0InvMassPt->Sumw2();
        fHistoTrueClusPartConvEtaPtM02          = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusPartConvFromEtaM02.Data());
        if (fHistoTrueClusPartConvEtaPtM02) fHistoTrueClusPartConvEtaPtM02->Sumw2();
        fHistoTrueClusPartConvEtaInvMassPt      = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusPartConvFromEtaInvMass.Data());
        if (fHistoTrueClusPartConvEtaInvMassPt) fHistoTrueClusPartConvEtaInvMassPt->Sumw2();
        FillMCInvMassHistosArray(   fHistoTrueClustersMergedInvMassPt, fHistoTrueClustersPi0InvMassPt, fHistoTrueClustersEtaInvMassPt, fHistoTrueClustersGammaInvMassPt, 
                                    fHistoTrueClustersElectronInvMassPt, fHistoTrueClustersBGInvMassPt, fHistoTrueClusPartConvMergedInvMassPt, fHistoTrueClusPartConvPi0InvMassPt,
                                    fHistoTrueClusPartConvEtaInvMassPt
                                );
        FillMCM02HistosArray(       fHistoTrueClustersMergedPtM02, fHistoTrueClustersPi0PtM02, fHistoTrueClustersEtaPtM02, fHistoTrueClustersGammaPtM02, 
                                    fHistoTrueClustersElectronPtM02, fHistoTrueClustersBGPtM02, fHistoTrueClusPartConvMergedPtM02, fHistoTrueClusPartConvPi0PtM02,
                                    fHistoTrueClusPartConvEtaPtM02
                            );

        
        fHistoTrueClustersBGPtSource            = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusBG_Source.Data());
        if (fHistoTrueClustersBGPtSource){
            fHistoTrueClustersBGPtSource->Sumw2();
            FillMCBGSeparated(fHistoTrueClustersBGPtSource);
        }
        
        if (meson.Contains("Pi0")){
            fHistoTrueClustersPrimPi0PtM02                  = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusPrimMesonM02.Data());
            if (fHistoTrueClustersPrimPi0PtM02) fHistoTrueClustersPrimPi0PtM02->Sumw2();
            fHistoTrueClustersPrimPi0InvMassPt              = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusPrimMesonInvMass.Data());
            if (fHistoTrueClustersPrimPi0InvMassPt) fHistoTrueClustersPrimPi0InvMassPt->Sumw2();
            fHistoTrueClustersSecPi0PtM02                   = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusSecMesonM02.Data());
            if (fHistoTrueClustersSecPi0PtM02) fHistoTrueClustersSecPi0PtM02->Sumw2();
            fHistoTrueClustersSecPi0InvMassPt               = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusSecMesonInvMass.Data());
            if (fHistoTrueClustersSecPi0InvMassPt) fHistoTrueClustersSecPi0InvMassPt->Sumw2();
            fHistoTrueClustersSecPi0FK0sPtM02               = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusSecMesonFromK0sM02.Data());
            if (fHistoTrueClustersSecPi0FK0sPtM02) fHistoTrueClustersSecPi0FK0sPtM02->Sumw2();
            fHistoTrueClustersSecPi0FK0sInvMassPt           = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusSecMesonFromK0sInvMass.Data());
            if (fHistoTrueClustersSecPi0FK0sInvMassPt) fHistoTrueClustersSecPi0FK0sInvMassPt->Sumw2();
            fHistoTrueClustersSecPi0FLambdaPtM02            = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusSecMesonFromLambdaM02.Data());
            if (fHistoTrueClustersSecPi0FLambdaPtM02) fHistoTrueClustersSecPi0FLambdaPtM02->Sumw2();
            fHistoTrueClustersSecPi0FLambdaInvMassPt        = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusSecMesonFromLambdaInvMass.Data());
            if (fHistoTrueClustersSecPi0FLambdaInvMassPt) fHistoTrueClustersSecPi0FLambdaInvMassPt->Sumw2();

            fHistoTrueClusPartConvPrimPi0PtM02              = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusPartConvPrimMesonM02.Data());
            if (fHistoTrueClusPartConvPrimPi0PtM02) fHistoTrueClusPartConvPrimPi0PtM02->Sumw2();
            fHistoTrueClusPartConvPrimPi0InvMassPt          = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusPartConvPrimMesonInvMass.Data());
            if (fHistoTrueClusPartConvPrimPi0InvMassPt) fHistoTrueClusPartConvPrimPi0InvMassPt->Sumw2();
            fHistoTrueClusPartConvSecPi0PtM02               = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusPartConvSecMesonM02.Data());
            if (fHistoTrueClusPartConvSecPi0PtM02) fHistoTrueClusPartConvSecPi0PtM02->Sumw2();
            fHistoTrueClusPartConvSecPi0InvMassPt           = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusPartConvSecMesonInvMass.Data());
            if (fHistoTrueClusPartConvSecPi0InvMassPt) fHistoTrueClusPartConvSecPi0InvMassPt->Sumw2();
            fHistoTrueClusPartConvSecPi0FK0sPtM02           = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusPartConvSecMesonFromK0sM02.Data());
            if (fHistoTrueClusPartConvSecPi0FK0sPtM02) fHistoTrueClusPartConvSecPi0FK0sPtM02->Sumw2();
            fHistoTrueClusPartConvSecPi0FK0sInvMassPt       = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusPartConvSecMesonFromK0sInvMass.Data());
            if (fHistoTrueClusPartConvSecPi0FK0sInvMassPt) fHistoTrueClusPartConvSecPi0FK0sInvMassPt->Sumw2();
            fHistoTrueClusPartConvSecPi0FLambdaPtM02        = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusPartConvSecMesonFromLambdaM02.Data());
            if (fHistoTrueClusPartConvSecPi0FLambdaPtM02) fHistoTrueClusPartConvSecPi0FLambdaPtM02->Sumw2();
            fHistoTrueClusPartConvSecPi0FLambdaInvMassPt    = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusPartConvSecMesonFromLambdaInvMass.Data());
            if (fHistoTrueClusPartConvSecPi0FLambdaInvMassPt) fHistoTrueClusPartConvSecPi0FLambdaInvMassPt->Sumw2();            

            FillMCPrimSecInvMassHistosArray(   fHistoTrueClustersPrimPi0InvMassPt, fHistoTrueClustersSecPi0InvMassPt, fHistoTrueClustersSecPi0FK0sInvMassPt, fHistoTrueClustersSecPi0FLambdaInvMassPt, 
                                               fHistoTrueClusPartConvPrimPi0InvMassPt, fHistoTrueClusPartConvSecPi0InvMassPt, fHistoTrueClusPartConvSecPi0FK0sInvMassPt, 
                                               fHistoTrueClusPartConvSecPi0FLambdaInvMassPt
                                           );
            FillMCPrimSecM02HistosArray(    fHistoTrueClustersPrimPi0PtM02, fHistoTrueClustersSecPi0PtM02, fHistoTrueClustersSecPi0FK0sPtM02, fHistoTrueClustersSecPi0FLambdaPtM02, 
                                            fHistoTrueClusPartConvPrimPi0PtM02, fHistoTrueClusPartConvSecPi0PtM02, fHistoTrueClusPartConvSecPi0FK0sPtM02, fHistoTrueClusPartConvSecPi0FLambdaPtM02
                                       );
        }
        
        
        TList *MCContainer                      = (TList*)HistosGammaConversion->FindObject(Form("%s MC histograms",fCutSelectionRead.Data()));
        
        fHistoMCMesonWithinAccepPt              = (TH1D*)MCContainer->FindObject(fObjectNameMCMesonAcc.Data());
        fHistoMCMesonPt                         = (TH1D*)MCContainer->FindObject(fObjectNameMCMeson.Data());  
        fHistoMCMesonPtWOWeights                = (TH1D*)MCContainer->FindObject(fObjectNameMCMesonWOWeights.Data());
        if (fIsMCGammaTrig){
            Double_t acceptanceCorrFacEMC       = 1/0.277; //Jet-Jet Gamma trigg on EMC phi (80-180), \eta \pm 0.7
            fHistoMCMesonPt->Sumw2();
            fHistoMCMesonPtWOWeights->Sumw2();
            fHistoMCMesonPt->Scale(acceptanceCorrFacEMC);
            fHistoMCMesonPtWOWeights->Scale(acceptanceCorrFacEMC);
        }    
    }

    //******************************************************************************************
    //************************* Analysing single pt bin ****************************************
    //******************************************************************************************    
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){ // BEGIN ANALYSIS for each Pt bin
        cout << "BIN: " << iPt << endl; 
        IntegrateHistoM02(fHistoM02PtBin[iPt], fMesonM02IntRange );
        fMesonM02Yields[iPt]        = fYields;
        fMesonM02YieldsError[iPt]   = fYieldsError;
        
        if (fIsMC){
            IntegrateHistoM02(fHistoTrueClusMergedM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueMergedYields[iPt]               = fYields;
            fMesonM02TrueMergedYieldsError[iPt]          = fYieldsError;
            
            IntegrateHistoM02(fHistoTrueClusPartConvMergedM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueMergedPartConvYields[iPt]       = fYields;
            fMesonM02TrueMergedPartConvYieldsError[iPt]  = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusFullMergedM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueMergedFullYields[iPt]           = fYields;
            fMesonM02TrueMergedFullYieldsError[iPt]      = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusPi0M02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TruePi0Yields[iPt]                  = fYields;
            fMesonM02TruePi0YieldsError[iPt]             = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusPartConvPi0M02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TruePi0PartConvYields[iPt]          = fYields;
            fMesonM02TruePi0PartConvYieldsError[iPt]     = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusFullPi0M02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TruePi0FullYields[iPt]              = fYields;
            fMesonM02TruePi0FullYieldsError[iPt]         = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusEtaM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueEtaYields[iPt]                  = fYields;
            fMesonM02TrueEtaYieldsError[iPt]             = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusPartConvEtaM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueEtaPartConvYields[iPt]          = fYields;
            fMesonM02TrueEtaPartConvYieldsError[iPt]     = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusFullEtaM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueEtaFullYields[iPt]              = fYields;
            fMesonM02TrueEtaFullYieldsError[iPt]         = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusGammaM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueGammaYields[iPt]                = fYields;
            fMesonM02TrueGammaYieldsError[iPt]           = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusElectronM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueElectronYields[iPt]             = fYields;
            fMesonM02TrueElectronYieldsError[iPt]        = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusBGM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueBGYields[iPt]                   = fYields;
            fMesonM02TrueBGYieldsError[iPt]              = fYieldsError;
            
            if (meson.Contains("Pi0")){
                IntegrateHistoM02(fHistoTrueClusPrimPi0M02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TruePrimPi0Yields[iPt]                     = fYields;
                fMesonM02TruePrimPi0YieldsError[iPt]                = fYieldsError;

                IntegrateHistoM02(fHistoTrueClusPartConvPrimPi0M02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TruePrimPi0PartConvYields[iPt]             = fYields;
                fMesonM02TruePrimPi0PartConvYieldsError[iPt]        = fYieldsError;

                IntegrateHistoM02(fHistoTrueClusFullPrimPi0M02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TruePrimPi0FullYields[iPt]                 = fYields;
                fMesonM02TruePrimPi0FullYieldsError[iPt]            = fYieldsError;

                IntegrateHistoM02(fHistoTrueClusSecPi0M02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TrueSecPi0Yields[iPt]                      = fYields;
                fMesonM02TrueSecPi0YieldsError[iPt]                 = fYieldsError;

                IntegrateHistoM02(fHistoTrueClusPartConvSecPi0M02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TrueSecPi0PartConvYields[iPt]              = fYields;
                fMesonM02TrueSecPi0PartConvYieldsError[iPt]         = fYieldsError;

                IntegrateHistoM02(fHistoTrueClusFullSecPi0M02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TrueSecPi0FullYields[iPt]                  = fYields;
                fMesonM02TrueSecPi0FullYieldsError[iPt]             = fYieldsError;

                IntegrateHistoM02(fHistoTrueClusSecPi0FK0sM02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TrueSecPi0FK0sYields[iPt]                  = fYields;
                fMesonM02TrueSecPi0FK0sYieldsError[iPt]             = fYieldsError;

                IntegrateHistoM02(fHistoTrueClusPartConvSecPi0FK0sM02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TrueSecPi0FK0sPartConvYields[iPt]          = fYields;
                fMesonM02TrueSecPi0FK0sPartConvYieldsError[iPt]     = fYieldsError;

                IntegrateHistoM02(fHistoTrueClusFullSecPi0FK0sM02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TrueSecPi0FK0sFullYields[iPt]              = fYields;
                fMesonM02TrueSecPi0FK0sFullYieldsError[iPt]         = fYieldsError;

                IntegrateHistoM02(fHistoTrueClusSecPi0FLambdaM02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TrueSecPi0FLambdaYields[iPt]               = fYields;
                fMesonM02TrueSecPi0FLambdaYieldsError[iPt]          = fYieldsError;

                IntegrateHistoM02(fHistoTrueClusPartConvSecPi0FLambdaM02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TrueSecPi0FLambdaPartConvYields[iPt]       = fYields;
                fMesonM02TrueSecPi0FLambdaPartConvYieldsError[iPt]  = fYieldsError;

                IntegrateHistoM02(fHistoTrueClusFullSecPi0FLambdaM02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TrueSecPi0FLambdaFullYields[iPt]           = fYields;
                fMesonM02TrueSecPi0FLambdaFullYieldsError[iPt]      = fYieldsError;                
            }
        }
    }
    
    CreatePtHistos();
    FillPtHistos();
    
    if (fIsMC){
        FillHistosArrayMC(fHistoMCMesonWithinAccepPt, fHistoMCMesonPt, fDeltaPt);
    } 
    
    //******************************************************************************************
    //************************* Plotting gamma/merged candidates per event *********************
    //******************************************************************************************
    fHistoClusterCandidates->Sumw2();
    fHistoClusterCandidates->Scale(1./fNEvents);
    fHistoClusterMergedCandidates->Sumw2();
    fHistoClusterMergedCandidates->Scale(1./fNEvents);
    TCanvas* canvasNClusters = new TCanvas("canvasNClusters","",2250,1500);
    DrawGammaCanvasSettings(canvasNClusters, 0.07, 0.01, 0.01, 0.075);
    canvasNClusters->SetLogy(1);
        Double_t maxYNClus = 30;
        if (fEnergyFlag.Contains("Pb")) maxYNClus = 100;
        DrawAutoGammaMesonHistos(   fHistoClusterCandidates, 
                                    "", "cluster candidates", "counts/event", 
                                    kTRUE, 2.,1e-6, kFALSE,
                                    kFALSE, 0, 20, 
                                    kTRUE, 0., maxYNClus);
        fHistoClusterCandidates->GetYaxis()->SetTitleOffset(0.9);
        DrawGammaSetMarker(fHistoClusterCandidates, 20, 1.5, kBlack, kBlack); 
        fHistoClusterCandidates->Draw("e,p");
        DrawGammaSetMarker(fHistoClusterMergedCandidates, 24, 1.5, kRed+2, kRed+2); 
        fHistoClusterMergedCandidates->Draw("same,e,p");

        PutProcessLabelAndEnergyOnPlot(0.55, 0.97, 0.045, fCollisionSystem.Data(), fDetectionProcess.Data(), "", 42, 0.03, "", 1, 1.1);
        
        TLegend* legendNClusters = GetAndSetLegend(0.5,0.75,2);
        legendNClusters->AddEntry(fHistoClusterCandidates,"passed general clusters cuts");
        legendNClusters->AddEntry(fHistoClusterMergedCandidates,"passed merged cluster cuts");
        legendNClusters->Draw();

        
    canvasNClusters->Update();
    canvasNClusters->SaveAs(Form("%s/%s_NClusterPerEvent%s.%s", outputDir.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    
    
    //******************************************************************************************
    //************************* Plotting 2D M02 vs Pt ******************************************
    //******************************************************************************************
    Double_t minZM02 = 0;
    Double_t maxZM02 = 0;
    
    TCanvas* canvasPtM02 = new TCanvas("canvasPtM02","",2250,1500);
    DrawGammaCanvasSettings(canvasPtM02, 0.085, 0.12, 0.02, 0.1);
    canvasPtM02->SetLogx(1); 
    canvasPtM02->SetLogy(0); 
    canvasPtM02->SetLogz(1); 
            
    canvasPtM02->cd();
    
    if(fHistoClustersMergedPtM02){
        fHistoClustersMergedPtM02->Sumw2();
        fHistoClustersMergedPtM02->Scale(1./fNEvents);
        maxZM02     = fHistoClustersMergedPtM02->GetMaximum();
        minZM02     = FindSmallestEntryIn2D(fHistoClustersMergedPtM02);
        DrawAutoGammaHistoPaper2D(fHistoClustersMergedPtM02,
                                " ",
                                "#it{p}_{T} (GeV/#it{c})",
                                "#lambda_{0}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,2.95,50.05,0.8,0.8);
        fHistoClustersMergedPtM02->GetXaxis()->SetMoreLogLabels();
        fHistoClustersMergedPtM02->GetXaxis()->SetLabelOffset(-0.02);
        fHistoClustersMergedPtM02->GetZaxis()->SetLabelOffset(-0.008);
        fHistoClustersMergedPtM02->GetZaxis()->SetLabelSize(0.051);
        fHistoClustersMergedPtM02->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        fHistoClustersMergedPtM02->GetXaxis()->SetTickLength(0.05);
        fHistoClustersMergedPtM02->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.55, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin);
        
        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_PtVsM02_AllAcceptedClusters%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //**************** Plotting 2D M02 vs Pt merged clusters passing meson cuts ****************
    //******************************************************************************************    
    if(fHistoClustersMergedPtM02AccMeson){
        fHistoClustersMergedPtM02AccMeson->Sumw2();
        fHistoClustersMergedPtM02AccMeson->Scale(1./fNEvents);
        if (minZM02 == 0) minZM02     = FindSmallestEntryIn2D(fHistoClustersMergedPtM02AccMeson);
        if (maxZM02 == 0) maxZM02     = fHistoClustersMergedPtM02AccMeson->GetMaximum();
        DrawAutoGammaHistoPaper2D(fHistoClustersMergedPtM02AccMeson,
                                " ",
                                "#it{p}_{T} (GeV/#it{c})",
                                "#lambda_{0}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,2.95,50.05,0.8,0.8);
        fHistoClustersMergedPtM02AccMeson->GetXaxis()->SetMoreLogLabels();
        fHistoClustersMergedPtM02AccMeson->GetXaxis()->SetLabelOffset(-0.02);
        fHistoClustersMergedPtM02AccMeson->GetZaxis()->SetLabelSize(0.051);
        fHistoClustersMergedPtM02AccMeson->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        fHistoClustersMergedPtM02AccMeson->GetXaxis()->SetTickLength(0.05);
        fHistoClustersMergedPtM02AccMeson->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.12, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);
        
        if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin);
        
        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_PtVsM02_AllAcceptedMesons%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //******************** Plotting 2D M02 vs Pt validated merged clusters *********************
    //******************************************************************************************        
    if (fHistoTrueClustersMergedPtM02 && fHistoTrueClusPartConvMergedPtM02){
        fHistoTrueClustersMergedPtM02->Sumw2();
        fHistoTrueClusPartConvMergedPtM02->Sumw2();
        TH2F* dumm2D = (TH2F*)fHistoTrueClustersMergedPtM02->Clone("forplotting");
        dumm2D->Add(fHistoTrueClusPartConvMergedPtM02);
        
        dumm2D->Scale(1./fNEvents);
        DrawAutoGammaHistoPaper2D(dumm2D,
                                " ",
                                "#it{p}_{T} (GeV/#it{c})",
                                "#lambda_{0}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,2.95,50.05,0.8,0.8);
        dumm2D->GetXaxis()->SetMoreLogLabels();
        dumm2D->GetXaxis()->SetLabelOffset(-0.02);
        dumm2D->GetZaxis()->SetLabelOffset(-0.008);
        dumm2D->GetZaxis()->SetLabelSize(0.051);
        dumm2D->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        dumm2D->GetXaxis()->SetTickLength(0.05);
        dumm2D->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.12, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);
        
        TLatex *labelM02 = new TLatex(0.11, 0.15, "val. merged clusters");
        SetStyleTLatex( labelM02, 0.05,4);
        labelM02->Draw();
        if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin);
        
        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_PtVsM02_TrueMerged%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }
    
    //******************************************************************************************
    //******************** Plotting 2D M02 vs Pt validated merged pi0s *************************
    //******************************************************************************************        
    if (fHistoTrueClustersPi0PtM02 && fHistoTrueClusPartConvPi0PtM02){
        fHistoTrueClustersPi0PtM02->Sumw2();
        fHistoTrueClusPartConvPi0PtM02->Sumw2();
        TH2F* dumm2D = (TH2F*)fHistoTrueClustersPi0PtM02->Clone("forplotting");
        dumm2D->Add(fHistoTrueClusPartConvPi0PtM02);
        
        dumm2D->Scale(1./fNEvents);
        DrawAutoGammaHistoPaper2D(dumm2D,
                                " ",
                                "#it{p}_{T} (GeV/#it{c})",
                                "#lambda_{0}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,2.95,50.05,0.8,0.8);
        dumm2D->GetXaxis()->SetMoreLogLabels();
        dumm2D->GetXaxis()->SetLabelOffset(-0.02);
        dumm2D->GetZaxis()->SetLabelOffset(-0.008);
        dumm2D->GetZaxis()->SetLabelSize(0.051);
        dumm2D->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        dumm2D->GetXaxis()->SetTickLength(0.05);
        dumm2D->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.12, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin);
        TLatex *labelM02 = new TLatex(0.11, 0.15, "val. merged #pi^{0}");
        SetStyleTLatex( labelM02, 0.05,4);
        labelM02->Draw();
        
        
        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_PtVsM02_TruePi0%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //******************** Plotting 2D M02 vs Pt validated merged etas *************************
    //******************************************************************************************            
    if (fHistoTrueClustersEtaPtM02 && fHistoTrueClusPartConvEtaPtM02){
        fHistoTrueClustersEtaPtM02->Sumw2();
        fHistoTrueClusPartConvEtaPtM02->Sumw2();
        TH2F* dumm2D = (TH2F*)fHistoTrueClustersEtaPtM02->Clone("forplotting");
        dumm2D->Add(fHistoTrueClusPartConvEtaPtM02);
        dumm2D->Scale(1./fNEvents);
        DrawAutoGammaHistoPaper2D(dumm2D,
                                " ",
                                "#it{p}_{T} (GeV/#it{c})",
                                "#lambda_{0}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,2.95,50.05,0.8,0.8);
        dumm2D->GetXaxis()->SetMoreLogLabels();
        dumm2D->GetXaxis()->SetLabelOffset(-0.02);
        dumm2D->GetZaxis()->SetLabelOffset(-0.008);
        dumm2D->GetZaxis()->SetLabelSize(0.051);
        dumm2D->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        dumm2D->GetXaxis()->SetTickLength(0.05);
        dumm2D->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.12, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        TLatex *labelM02 = new TLatex(0.11, 0.15, "val. merged #eta");
        SetStyleTLatex( labelM02, 0.05,4);
        labelM02->Draw();
        
        if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin);
        
        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_PtVsM02_TrueEta%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //******************** Plotting 2D M02 vs Pt validated gammas ******************************
    //******************************************************************************************            
    if (fHistoTrueClustersGammaPtM02){
        fHistoTrueClustersGammaPtM02->Scale(1./fNEvents);
        DrawAutoGammaHistoPaper2D(fHistoTrueClustersGammaPtM02,
                                " ",
                                "#it{p}_{T} (GeV/#it{c})",
                                "#lambda_{0}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,2.95,50.05,0.8,0.8);
        fHistoTrueClustersGammaPtM02->GetXaxis()->SetMoreLogLabels();
        fHistoTrueClustersGammaPtM02->GetXaxis()->SetLabelOffset(-0.02);
        fHistoTrueClustersGammaPtM02->GetZaxis()->SetLabelOffset(-0.008);
        fHistoTrueClustersGammaPtM02->GetZaxis()->SetLabelSize(0.051);
        fHistoTrueClustersGammaPtM02->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        fHistoTrueClustersGammaPtM02->GetXaxis()->SetTickLength(0.05);
        fHistoTrueClustersGammaPtM02->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.12, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        TLatex *labelM02 = new TLatex(0.11, 0.15, "val. #gamma");
        SetStyleTLatex( labelM02, 0.05,4);
        labelM02->Draw();
        if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin);
        
        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_PtVsM02_TrueGamma%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //******************** Plotting 2D M02 vs Pt validated electrons ***************************
    //******************************************************************************************            
    if (fHistoTrueClustersElectronPtM02){
        fHistoTrueClustersElectronPtM02->Sumw2();
        fHistoTrueClustersElectronPtM02->Scale(1./fNEvents);
        DrawAutoGammaHistoPaper2D(fHistoTrueClustersElectronPtM02,
                                " ",
                                "#it{p}_{T} (GeV/#it{c})",
                                "#lambda_{0}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,2.95,50.05,0.8,0.8);
        fHistoTrueClustersElectronPtM02->GetXaxis()->SetMoreLogLabels();
        fHistoTrueClustersElectronPtM02->GetXaxis()->SetLabelOffset(-0.02);
        fHistoTrueClustersElectronPtM02->GetZaxis()->SetLabelOffset(-0.008);
        fHistoTrueClustersElectronPtM02->GetZaxis()->SetLabelSize(0.051);
        fHistoTrueClustersElectronPtM02->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        fHistoTrueClustersElectronPtM02->GetXaxis()->SetTickLength(0.05);
        fHistoTrueClustersElectronPtM02->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.12, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        TLatex *labelM02 = new TLatex(0.11, 0.15, "val. e^{#pm}");
        SetStyleTLatex( labelM02, 0.05,4);
        labelM02->Draw();
     
        if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin);
        
        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_PtVsM02_TrueElectron%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }
     
    //******************************************************************************************
    //******************** Plotting 2D M02 vs Pt validated merged BG ***************************
    //******************************************************************************************            
    if (fHistoTrueClustersBGPtM02){
        fHistoTrueClustersBGPtM02->Sumw2();
        fHistoTrueClustersBGPtM02->Scale(1./fNEvents);
        DrawAutoGammaHistoPaper2D(fHistoTrueClustersBGPtM02,
                                " ",
                                "#it{p}_{T} (GeV/#it{c})",
                                "#lambda_{0}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,2.95,50.05,0.8,0.8);
        fHistoTrueClustersBGPtM02->GetXaxis()->SetMoreLogLabels();
        fHistoTrueClustersBGPtM02->GetXaxis()->SetLabelOffset(-0.02);
        fHistoTrueClustersBGPtM02->GetZaxis()->SetLabelOffset(-0.008);
        fHistoTrueClustersBGPtM02->GetZaxis()->SetLabelSize(0.051);
        fHistoTrueClustersBGPtM02->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        fHistoTrueClustersBGPtM02->GetXaxis()->SetTickLength(0.05);
        fHistoTrueClustersBGPtM02->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.12, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        TLatex *labelM02 = new TLatex(0.11, 0.15, "background");
        SetStyleTLatex( labelM02, 0.05,4);
        labelM02->Draw();
        
        if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin);
        
        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_PtVsM02_TrueBG%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    
    //******************************************************************************************
    //***************** Plotting gamma candidates per event vs track candidates ****************
    //******************************************************************************************
    if(fHistoTrackVsClusterCandidates){
        TCanvas* canvasTracksVsGamma = new TCanvas("canvasTracksVsGamma","",2250,1500);
        DrawGammaCanvasSettings(canvasTracksVsGamma, 0.08, 0.12, 0.015, 0.1);
        canvasTracksVsGamma->SetLogz(1);
        
            Double_t maxXNTracks = 100;
            fHistoTrackVsClusterCandidates->Sumw2();
            fHistoTrackVsClusterCandidates->Scale(1./fNEvents);
            Double_t minZClus     = FindSmallestEntryIn2D(fHistoTrackVsClusterCandidates);
            Double_t maxZClus     = fHistoTrackVsClusterCandidates->GetMaximum();
            DrawAutoGammaHistoPaper2D(fHistoTrackVsClusterCandidates,
                                    "",
                                    "ESD track candidates",
                                    "cluster candidates",
                                    0,0,0,
                                    1,0,maxYNClus,
                                    1,0,maxXNTracks,0.8,0.75);
            fHistoTrackVsClusterCandidates->GetZaxis()->SetLabelSize(0.051);
            fHistoTrackVsClusterCandidates->GetZaxis()->SetLabelOffset(-0.008);
            fHistoTrackVsClusterCandidates->GetZaxis()->SetRangeUser(minZClus,maxZClus);
            fHistoTrackVsClusterCandidates->GetXaxis()->SetTickLength(0.05);
            fHistoTrackVsClusterCandidates->DrawCopy("COLZ");
            PutProcessLabelAndEnergyOnPlot(0.12, 0.93, 0.045, fCollisionSystem.Data(), fDetectionProcess.Data(), "", 42, 0.03, "", 1, 1.1);

        canvasTracksVsGamma->Update();
        canvasTracksVsGamma->SaveAs(Form("%s/%s_TrackVsGamma%s.%s", outputDir.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }
    
    //******************************************************************************************
    //**************************** Plotting SPD tracklet vs SPD clusters ***********************
    //******************************************************************************************    
    if (fHistoSPDtrackletvsSPDclusters){
        TCanvas* canvasSPDTracksVsCluster = new TCanvas("canvasSPDTracksVsCluster","",2250,1500);
        DrawGammaCanvasSettings(canvasSPDTracksVsCluster, 0.08, 0.12, 0.045, 0.1);
        canvasSPDTracksVsCluster->SetLogz(1);
        
            fHistoSPDtrackletvsSPDclusters->Sumw2();
            fHistoSPDtrackletvsSPDclusters->Scale(1./fNEvents);
            Double_t minZSPD        = FindSmallestEntryIn2D(fHistoSPDtrackletvsSPDclusters);
            Double_t maxZSPD        = fHistoSPDtrackletvsSPDclusters->GetMaximum();
            DrawAutoGammaHistoPaper2D(fHistoSPDtrackletvsSPDclusters,
                                    "",
                                    "SPD tracklets",
                                    "SPD clusters",
                                    0,0,0,
                                    0,0,0,
                                    0,0,0,0.8,0.75);
            fHistoSPDtrackletvsSPDclusters->GetZaxis()->SetLabelSize(0.051);
            fHistoSPDtrackletvsSPDclusters->GetZaxis()->SetLabelOffset(-0.008);
            fHistoSPDtrackletvsSPDclusters->GetZaxis()->SetRangeUser(minZSPD,maxZSPD);
            fHistoSPDtrackletvsSPDclusters->GetXaxis()->SetTickLength(0.05);
            fHistoSPDtrackletvsSPDclusters->DrawCopy("COLZ");
            PutProcessLabelAndEnergyOnPlot(0.6, 0.93, 0.045, fCollisionSystem.Data(), "", "", 42, 0.03, "", 1, 1.1);

        canvasSPDTracksVsCluster->Update();
        canvasSPDTracksVsCluster->SaveAs(Form("%s/%s_SPDtrackletVsSPDclusters%s.%s", outputDir.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //******************************** Plotting InvMass versus Pt ******************************
    //******************************************************************************************        
    TCanvas* canvasInvMassVsPt = new TCanvas("canvasInvMassVsPt","",2250,1500);
    DrawGammaCanvasSettings(canvasInvMassVsPt, 0.08, 0.12, 0.045, 0.1);
    canvasInvMassVsPt->SetLogz(1);
    Double_t minZInv        = 0;
    Double_t maxZInv        = 0;
 
    if (fHistoInvMassVsPt){        
        fHistoInvMassVsPt->Sumw2();
        minZInv             = FindSmallestEntryIn2D(fHistoInvMassVsPt);
        maxZInv             = fHistoInvMassVsPt->GetMaximum();
        DrawAutoGammaHistoPaper2D(fHistoInvMassVsPt,
                                "",
                                "#it{M}_{inv} (GeV/#it{c})",
                                "#it{p}_{T} (GeV/#it{c})",
                                0,0,0,
                                0,0,0,
                                0,0,0,0.85,0.75);
        fHistoInvMassVsPt->GetZaxis()->SetLabelSize(0.051);
        fHistoInvMassVsPt->GetZaxis()->SetLabelOffset(-0.008);
        fHistoInvMassVsPt->GetZaxis()->SetRangeUser(minZInv,maxZInv);
        fHistoInvMassVsPt->GetXaxis()->SetTickLength(0.05);
        fHistoInvMassVsPt->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(startXLabelInvMass, 0.93, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        canvasInvMassVsPt->Update();
        canvasInvMassVsPt->SaveAs(Form("%s/%s_%s_InvMassVsPt%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //********************* Plotting true merged clusters InvMass versus Pt ********************
    //******************************************************************************************            
    if (fHistoTrueClustersMergedInvMassPt){
        fHistoTrueClustersMergedInvMassPt->Sumw2();
        DrawAutoGammaHistoPaper2D(fHistoTrueClustersMergedInvMassPt,
                                "",
                                "#it{M}_{inv} (GeV/#it{c})",
                                "#it{p}_{T} (GeV/#it{c})",
                                0,0,0,
                                0,0,0,
                                0,0,0,0.85,0.75);
        fHistoTrueClustersMergedInvMassPt->GetZaxis()->SetLabelSize(0.051);
        fHistoTrueClustersMergedInvMassPt->GetZaxis()->SetLabelOffset(-0.008);
        fHistoTrueClustersMergedInvMassPt->GetZaxis()->SetRangeUser(minZInv,maxZInv);
        fHistoTrueClustersMergedInvMassPt->GetXaxis()->SetTickLength(0.05);
        fHistoTrueClustersMergedInvMassPt->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(startXLabelInvMass, 0.93, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        TLatex *labelInvMass = new TLatex(0.55, 0.18, "val. merged clusters");
        SetStyleTLatex( labelInvMass, 0.05,4);
        labelInvMass->Draw();

        
        canvasInvMassVsPt->Update();
        canvasInvMassVsPt->SaveAs(Form("%s/%s_%s_TrueMerged_InvMassVsPt%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //************************ Plotting true pi0 clusters InvMass versus Pt ********************
    //******************************************************************************************                
    if (fHistoTrueClustersPi0InvMassPt){
        fHistoTrueClustersPi0InvMassPt->Sumw2();
        DrawAutoGammaHistoPaper2D(fHistoTrueClustersPi0InvMassPt,
                                "",
                                "#it{M}_{inv} (GeV/#it{c})",
                                "#it{p}_{T} (GeV/#it{c})",
                                0,0,0,
                                0,0,0,
                                0,0,0,0.85,0.75);
        fHistoTrueClustersPi0InvMassPt->GetZaxis()->SetLabelSize(0.051);
        fHistoTrueClustersPi0InvMassPt->GetZaxis()->SetLabelOffset(-0.008);
        fHistoTrueClustersPi0InvMassPt->GetZaxis()->SetRangeUser(minZInv,maxZInv);
        fHistoTrueClustersPi0InvMassPt->GetXaxis()->SetTickLength(0.05);
        fHistoTrueClustersPi0InvMassPt->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(startXLabelInvMass, 0.93, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        TLatex *labelInvMass = new TLatex(0.65, 0.18, "val. merged #pi^{0}");
        SetStyleTLatex( labelInvMass, 0.05,4);
        labelInvMass->Draw();
        
        canvasInvMassVsPt->Update();
        canvasInvMassVsPt->SaveAs(Form("%s/%s_%s_TruePi0_InvMassVsPt%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //************************ Plotting true eta clusters InvMass versus Pt ********************
    //******************************************************************************************                
    if (fHistoTrueClustersEtaInvMassPt){
        fHistoTrueClustersEtaInvMassPt->Sumw2();
        DrawAutoGammaHistoPaper2D(fHistoTrueClustersEtaInvMassPt,
                                "",
                                "#it{M}_{inv} (GeV/#it{c})",
                                "#it{p}_{T} (GeV/#it{c})",
                                0,0,0,
                                0,0,0,
                                0,0,0,0.85,0.75);
        fHistoTrueClustersEtaInvMassPt->GetZaxis()->SetLabelSize(0.051);
        fHistoTrueClustersEtaInvMassPt->GetZaxis()->SetLabelOffset(-0.008);
        fHistoTrueClustersEtaInvMassPt->GetZaxis()->SetRangeUser(minZInv,maxZInv);
        fHistoTrueClustersEtaInvMassPt->GetXaxis()->SetTickLength(0.05);
        fHistoTrueClustersEtaInvMassPt->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(startXLabelInvMass, 0.93, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        TLatex *labelInvMass = new TLatex(0.65, 0.18, "val. merged #eta");
        SetStyleTLatex( labelInvMass, 0.05,4);
        labelInvMass->Draw();
        
        canvasInvMassVsPt->Update();
        canvasInvMassVsPt->SaveAs(Form("%s/%s_%s_TrueEta_InvMassVsPt%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //************************ Plotting true gamma clusters InvMass versus Pt ********************
    //******************************************************************************************                
    if (fHistoTrueClustersGammaInvMassPt){
        fHistoTrueClustersGammaInvMassPt->Sumw2();
        DrawAutoGammaHistoPaper2D(fHistoTrueClustersGammaInvMassPt,
                                "",
                                "#it{M}_{inv} (GeV/#it{c})",
                                "#it{p}_{T} (GeV/#it{c})",
                                0,0,0,
                                0,0,0,
                                0,0,0,0.85,0.75);
        fHistoTrueClustersGammaInvMassPt->GetZaxis()->SetLabelSize(0.051);
        fHistoTrueClustersGammaInvMassPt->GetZaxis()->SetLabelOffset(-0.008);
        fHistoTrueClustersGammaInvMassPt->GetZaxis()->SetRangeUser(minZInv,maxZInv);
        fHistoTrueClustersGammaInvMassPt->GetXaxis()->SetTickLength(0.05);
        fHistoTrueClustersGammaInvMassPt->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(startXLabelInvMass, 0.93, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        TLatex *labelInvMass = new TLatex(0.75, 0.18, "val. #gamma");
        SetStyleTLatex( labelInvMass, 0.05,4);
        labelInvMass->Draw();

        
        canvasInvMassVsPt->Update();
        canvasInvMassVsPt->SaveAs(Form("%s/%s_%s_TrueGamma_InvMassVsPt%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //************************ Plotting true eta clusters InvMass versus Pt ********************
    //******************************************************************************************                
    if (fHistoTrueClustersElectronInvMassPt){
        fHistoTrueClustersElectronInvMassPt->Sumw2();
        DrawAutoGammaHistoPaper2D(fHistoTrueClustersElectronInvMassPt,
                                "",
                                "#it{M}_{inv} (GeV/#it{c})",
                                "#it{p}_{T} (GeV/#it{c})",
                                0,0,0,
                                0,0,0,
                                0,0,0,0.85,0.75);
        fHistoTrueClustersElectronInvMassPt->GetZaxis()->SetLabelSize(0.051);
        fHistoTrueClustersElectronInvMassPt->GetZaxis()->SetLabelOffset(-0.008);
        fHistoTrueClustersElectronInvMassPt->GetZaxis()->SetRangeUser(minZInv,maxZInv);
        fHistoTrueClustersElectronInvMassPt->GetXaxis()->SetTickLength(0.05);
        fHistoTrueClustersElectronInvMassPt->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(startXLabelInvMass, 0.93, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        TLatex *labelInvMass = new TLatex(0.75, 0.18, "val. e^{#pm}");
        SetStyleTLatex( labelInvMass, 0.05,4);
        labelInvMass->Draw();
        
        canvasInvMassVsPt->Update();
        canvasInvMassVsPt->SaveAs(Form("%s/%s_%s_TrueElectron_InvMassVsPt%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }
    
    //******************************************************************************************
    //************************* Plotting true BG clusters InvMass versus Pt ********************
    //******************************************************************************************                    
    if (fHistoTrueClustersBGInvMassPt){
        fHistoTrueClustersBGInvMassPt->Sumw2();
        DrawAutoGammaHistoPaper2D(fHistoTrueClustersBGInvMassPt,
                                "",
                                "#it{M}_{inv} (GeV/#it{c})",
                                "#it{p}_{T} (GeV/#it{c})",
                                0,0,0,
                                0,0,0,
                                0,0,0,0.85,0.75);
        fHistoTrueClustersBGInvMassPt->GetZaxis()->SetLabelSize(0.051);
        fHistoTrueClustersBGInvMassPt->GetZaxis()->SetLabelOffset(-0.008);
        fHistoTrueClustersBGInvMassPt->GetZaxis()->SetRangeUser(minZInv,maxZInv);
        fHistoTrueClustersBGInvMassPt->GetXaxis()->SetTickLength(0.05);
        fHistoTrueClustersBGInvMassPt->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(startXLabelInvMass, 0.93, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        TLatex *labelInvMass = new TLatex(0.65, 0.18, "background");
        SetStyleTLatex( labelInvMass, 0.05,4);
        labelInvMass->Draw();
        
        canvasInvMassVsPt->Update();
        canvasInvMassVsPt->SaveAs(Form("%s/%s_%s_TrueBG_InvMassVsPt%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }
    
    if (fAdvancedClusterQA){
        TCanvas* canvasNCellsVsPt = new TCanvas("canvasNCellsVsPt","",1900,1500);
        DrawGammaCanvasSettings(canvasNCellsVsPt, 0.09, 0.13, 0.02, 0.1);
        canvasNCellsVsPt->SetLogx(0); 
        canvasNCellsVsPt->SetLogy(1); 
        canvasNCellsVsPt->SetLogz(1);         
        canvasNCellsVsPt->cd();
        Double_t minZNcells     = 0;
        Double_t maxZNcells     = 0;
        
        if(fHistoClusterNCellsPt){
            fHistoClusterNCellsPt->Sumw2();
            fHistoClusterNCellsPt->Scale(1./fNEvents);
            maxZNcells     = fHistoClusterNCellsPt->GetMaximum();
            minZNcells     = FindSmallestEntryIn2D(fHistoClusterNCellsPt);

            DrawAutoGammaHistoPaper2D(fHistoClusterNCellsPt,
                                    " ",
                                    "#it{N}_{cells in cluster}",
                                    "#it{p}_{T} (GeV/#it{c})",
                                    0,0,0,
                                    1,0.7,50.0,
                                    1,-0.5,99.5,
                                    0.9,0.8);
            fHistoClusterNCellsPt->GetYaxis()->SetMoreLogLabels();
            fHistoClusterNCellsPt->GetYaxis()->SetLabelOffset(-0.001);
            fHistoClusterNCellsPt->GetZaxis()->SetLabelOffset(-0.008);
            fHistoClusterNCellsPt->GetZaxis()->SetLabelSize(0.051);
            fHistoClusterNCellsPt->GetXaxis()->SetTickLength(0.03);
            fHistoClusterNCellsPt->GetZaxis()->SetRangeUser(minZNcells,maxZNcells);
            fHistoClusterNCellsPt->DrawCopy("COLZ");
            PutProcessLabelAndEnergyOnPlot(0.52, 0.96, 0.045, fCollisionSystem.Data(), fDetectionProcess.Data(), "", 42, 0.03, "", 1, 1.1);

            canvasNCellsVsPt->Update();
            canvasNCellsVsPt->SaveAs(Form("%s/%s_%s_NCellsPt_AllClusters%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
        }

        Double_t startPtMergedCl    = 2.95;
        if (fNLMmin == 1){
            startPtMergedCl         = 8;
        }    
            
        if(fHistoClusterMergedNCellsPt){
            fHistoClusterMergedNCellsPt->Sumw2();
            fHistoClusterMergedNCellsPt->Scale(1./fNEvents);
            maxZNcells     = fHistoClusterMergedNCellsPt->GetMaximum();
            minZNcells     = FindSmallestEntryIn2D(fHistoClusterMergedNCellsPt);

            DrawAutoGammaHistoPaper2D(fHistoClusterMergedNCellsPt,
                                    " ",
                                    "#it{N}_{cells in merged cluster}",
                                    "#it{p}_{T} (GeV/#it{c})",
                                    0,0,0,
                                    1,startPtMergedCl,50.0,
                                    1,-0.5,99.5,
                                    0.9,0.8);
            fHistoClusterMergedNCellsPt->GetYaxis()->SetMoreLogLabels();
            fHistoClusterMergedNCellsPt->GetYaxis()->SetLabelOffset(-0.001);
            fHistoClusterMergedNCellsPt->GetZaxis()->SetLabelOffset(-0.008);
            fHistoClusterMergedNCellsPt->GetZaxis()->SetLabelSize(0.051);
            fHistoClusterMergedNCellsPt->GetXaxis()->SetTickLength(0.03);
            fHistoClusterMergedNCellsPt->GetZaxis()->SetRangeUser(minZNcells,maxZNcells);
            fHistoClusterMergedNCellsPt->DrawCopy("COLZ");
            PutProcessLabelAndEnergyOnPlot(0.52, 0.96, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

            canvasNCellsVsPt->Update();
            canvasNCellsVsPt->SaveAs(Form("%s/%s_%s_NCellsPt_MergedClusters%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
        }

        if(fHistoClusterMergedNCellsArClPt){
            fHistoClusterMergedNCellsArClPt->Sumw2();
            fHistoClusterMergedNCellsArClPt->Scale(1./fNEvents);
            maxZNcells     = fHistoClusterMergedNCellsArClPt->GetMaximum();
            minZNcells     = FindSmallestEntryIn2D(fHistoClusterMergedNCellsArClPt);

            DrawAutoGammaHistoPaper2D(fHistoClusterMergedNCellsArClPt,
                                    " ",
                                    "#it{N}_{cells around merged cluster (#Delta#it{R} < 0.15)}",
                                    "#it{p}_{T} (GeV/#it{c})",
                                    0,0,0,
                                    1,startPtMergedCl,50.0,
                                    1,-0.5,99.5,
                                    0.9,0.8);
            fHistoClusterMergedNCellsArClPt->GetYaxis()->SetMoreLogLabels();
            fHistoClusterMergedNCellsArClPt->GetYaxis()->SetLabelOffset(-0.001);
            fHistoClusterMergedNCellsArClPt->GetZaxis()->SetLabelOffset(-0.008);
            fHistoClusterMergedNCellsArClPt->GetZaxis()->SetLabelSize(0.051);
            fHistoClusterMergedNCellsArClPt->GetXaxis()->SetTickLength(0.03);
            fHistoClusterMergedNCellsArClPt->GetZaxis()->SetRangeUser(minZNcells,maxZNcells);
            fHistoClusterMergedNCellsArClPt->DrawCopy("COLZ");
            PutProcessLabelAndEnergyOnPlot(0.52, 0.96, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

            canvasNCellsVsPt->Update();
            canvasNCellsVsPt->SaveAs(Form("%s/%s_%s_NCellsAroundClPt_MergedClusters%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
        }

         if(fHistoClusterMergedNCellsArAInclPt){
            fHistoClusterMergedNCellsArAInclPt->Sumw2();
            fHistoClusterMergedNCellsArAInclPt->Scale(1./fNEvents);
            maxZNcells     = fHistoClusterMergedNCellsArAInclPt->GetMaximum();
            minZNcells     = FindSmallestEntryIn2D(fHistoClusterMergedNCellsArAInclPt);

            DrawAutoGammaHistoPaper2D(fHistoClusterMergedNCellsArAInclPt,
                                    " ",
                                    "#it{N}_{cells around merged cluster (#Delta#it{R} < 0.15)} + #it{N}_{cells in merged cluster}",
                                    "#it{p}_{T} (GeV/#it{c})",
                                    0,0,0,
                                    1,startPtMergedCl,50.0,
                                    1,-0.5,99.5,
                                    0.9,0.8);
            fHistoClusterMergedNCellsArAInclPt->GetYaxis()->SetMoreLogLabels();
            fHistoClusterMergedNCellsArAInclPt->GetYaxis()->SetLabelOffset(-0.001);
            fHistoClusterMergedNCellsArAInclPt->GetZaxis()->SetLabelOffset(-0.008);
            fHistoClusterMergedNCellsArAInclPt->GetZaxis()->SetLabelSize(0.051);
            fHistoClusterMergedNCellsArAInclPt->GetXaxis()->SetTickLength(0.03);
            fHistoClusterMergedNCellsArAInclPt->GetZaxis()->SetRangeUser(minZNcells,maxZNcells);
            fHistoClusterMergedNCellsArAInclPt->DrawCopy("COLZ");
            PutProcessLabelAndEnergyOnPlot(0.52, 0.96, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

            canvasNCellsVsPt->Update();
            canvasNCellsVsPt->SaveAs(Form("%s/%s_%s_NCellsAroundAndInClPt_MergedClusters%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
        }

        if(fHistoClusterMergedEAroundClE){
            fHistoClusterMergedEAroundClE->Sumw2();
            fHistoClusterMergedEAroundClE->Scale(1./fNEvents);
            maxZNcells     = fHistoClusterMergedEAroundClE->GetMaximum();
            minZNcells     = FindSmallestEntryIn2D(fHistoClusterMergedEAroundClE);

            DrawAutoGammaHistoPaper2D(fHistoClusterMergedEAroundClE,
                                    " ",
                                    "#it{E}_{around merged cluster (#Delta#it{R} < 0.15)} (GeV/#it{c})",
                                    "#it{E}_{merged cluster} (GeV/#it{c})",
                                    0,0,0,
                                    1,startPtMergedCl,50.0,
                                    1,0,100,
                                    0.9,0.8);
            fHistoClusterMergedEAroundClE->GetYaxis()->SetMoreLogLabels();
            fHistoClusterMergedEAroundClE->GetYaxis()->SetLabelOffset(-0.001);
            fHistoClusterMergedEAroundClE->GetZaxis()->SetLabelOffset(-0.008);
            fHistoClusterMergedEAroundClE->GetZaxis()->SetLabelSize(0.051);
            fHistoClusterMergedEAroundClE->GetXaxis()->SetTickLength(0.03);
            fHistoClusterMergedEAroundClE->GetZaxis()->SetRangeUser(minZNcells,maxZNcells);
            fHistoClusterMergedEAroundClE->DrawCopy("COLZ");
            PutProcessLabelAndEnergyOnPlot(0.52, 0.96, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

            canvasNCellsVsPt->Update();
            canvasNCellsVsPt->SaveAs(Form("%s/%s_%s_EAroundClE_MergedClusters%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
        }
        
        if(fHistoClusterMergedNPartPt){
            fHistoClusterMergedNPartPt->Sumw2();
            fHistoClusterMergedNPartPt->Scale(1./fNEvents);
            maxZNcells     = fHistoClusterMergedNPartPt->GetMaximum();
            minZNcells     = FindSmallestEntryIn2D(fHistoClusterMergedNPartPt);

            DrawAutoGammaHistoPaper2D(fHistoClusterMergedNPartPt,
                                    " ",
                                    "#it{N}_{MC part. point. to cl.}",
                                    "#it{p}_{T} (GeV/#it{c})",
                                    0,0,0,
                                    1,startPtMergedCl,50.0,
                                    1,0,50,
                                    0.9,0.8);
            fHistoClusterMergedNPartPt->GetYaxis()->SetMoreLogLabels();
            fHistoClusterMergedNPartPt->GetYaxis()->SetLabelOffset(-0.001);
            fHistoClusterMergedNPartPt->GetZaxis()->SetLabelOffset(-0.008);
            fHistoClusterMergedNPartPt->GetZaxis()->SetLabelSize(0.051);
            fHistoClusterMergedNPartPt->GetXaxis()->SetTickLength(0.03);
            fHistoClusterMergedNPartPt->GetZaxis()->SetRangeUser(minZNcells,maxZNcells);
            fHistoClusterMergedNPartPt->DrawCopy("COLZ");
            PutProcessLabelAndEnergyOnPlot(0.52, 0.96, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

            canvasNCellsVsPt->Update();
            canvasNCellsVsPt->SaveAs(Form("%s/%s_%s_NParticlePointingToClPt_MergedClusters%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
        }
        
        
        
    }    
    
    fFileErrLog.close();
    fFileDataLog.close();
    
    
    //******************************************************************************************
    //************************ Plotting transverse momentum bins separately ********************
    //******************************************************************************************                
    TString plotPrefix  = Form("%s/%s_%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data());
    TString plotSuffix  = Form("%s%s_%s.%s", fPeriodFlag.Data(), fAdditionalName.Data(), fCutSelection.Data(), suffix.Data());
    
    TString nameMeson   = Form("%s_MesonInvMass%s", plotPrefix.Data(), plotSuffix.Data());
    TString nameCanvas  = "MesonCanvas";
    TString namePad     = "MesonPad";
    cout << nameMeson.Data() << endl;
    PlotInvMassMergedInPtBins(  fHistoInvMassPtBin, nameMeson, nameCanvas, namePad, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, 
                                fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem);

    if (fIsMC){
        nameMeson       = Form("%s_MesonInvMassSplit%s", plotPrefix.Data(), plotSuffix.Data());
        cout << nameMeson.Data() << endl;
        PlotInvMassMergedMCInPtBins(    fHistoInvMassPtBin, fHistoTrueClusMergedInvMassPtBin, fHistoTrueClusPartConvMergedInvMassPtBin, 
                                        fHistoTrueClusBGInvMassPtBin, fHistoTrueClusGammaInvMassPtBin, fHistoTrueClusElectronInvMassPtBin,
                                        nameMeson, nameCanvas, namePad, fMesonMassPlotRange, 
                                        fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem);

        nameMeson       = Form("%s_MesonFullInvMassSplit%s", plotPrefix.Data(), plotSuffix.Data());
        cout << nameMeson.Data() << endl;
        PlotInvMassMergedMCInPtBins(    fHistoInvMassPtBin, fHistoTrueClusFullMergedInvMassPtBin, NULL, 
                                        fHistoTrueClusBGInvMassPtBin, fHistoTrueClusGammaInvMassPtBin, fHistoTrueClusElectronInvMassPtBin,
                                        nameMeson, nameCanvas, namePad, fMesonMassPlotRange, 
                                        fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem);

        nameMeson       = Form("%s_MesonInvMassValidatedDisentangled%s", plotPrefix.Data(), plotSuffix.Data());
        cout << nameMeson.Data() << endl;
        PlotInvMassMergedTrueInPtBins(  fHistoTrueClusMergedInvMassPtBin, fHistoTrueClusPartConvMergedInvMassPtBin, fHistoTrueClusPi0InvMassPtBin, 
                                        fHistoTrueClusPartConvPi0InvMassPtBin, fHistoTrueClusEtaInvMassPtBin, fHistoTrueClusPartConvEtaInvMassPtBin,
                                        nameMeson, nameCanvas, namePad, fMesonMassPlotRange, 
                                        fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem); 

        nameMeson       = Form("%s_MesonFullInvMassValidatedDisentangled%s", plotPrefix.Data(), plotSuffix.Data());
        cout << nameMeson.Data() << endl;
        PlotInvMassMergedTrueInPtBins(  fHistoTrueClusFullMergedInvMassPtBin, NULL, fHistoTrueClusFullPi0InvMassPtBin, 
                                        NULL, fHistoTrueClusFullEtaInvMassPtBin, NULL,
                                        nameMeson, nameCanvas, namePad, fMesonMassPlotRange, 
                                        fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem);
    }
    
    nameMeson   = Form("%s_MesonM02%s", plotPrefix.Data(), plotSuffix.Data());
    cout << nameMeson.Data() << endl;
    PlotM02MergedInPtBins(  fHistoM02PtBin, nameMeson, nameCanvas, namePad, fMesonM02PlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, 
                            fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem);

    if (fIsMC){
        nameMeson       = Form("%s_MesonM02Split%s", plotPrefix.Data(), plotSuffix.Data());
        cout << nameMeson.Data() << endl;
        PlotM02MergedMCInPtBins(    fHistoM02PtBin, fHistoTrueClusMergedM02PtBin, fHistoTrueClusPartConvMergedM02PtBin, 
                                        fHistoTrueClusBGM02PtBin, fHistoTrueClusGammaM02PtBin, fHistoTrueClusElectronM02PtBin,
                                        nameMeson, nameCanvas, namePad, fMesonM02PlotRange, 
                                        fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem);

        nameMeson       = Form("%s_MesonFullM02Split%s", plotPrefix.Data(), plotSuffix.Data());
        cout << nameMeson.Data() << endl;
        PlotM02MergedMCInPtBins(    fHistoM02PtBin, fHistoTrueClusFullMergedM02PtBin, NULL, 
                                        fHistoTrueClusBGM02PtBin, fHistoTrueClusGammaM02PtBin, fHistoTrueClusElectronM02PtBin,
                                        nameMeson, nameCanvas, namePad, fMesonM02PlotRange, 
                                        fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem);

        nameMeson       = Form("%s_MesonM02ValidatedDisentangled%s", plotPrefix.Data(), plotSuffix.Data());
        cout << nameMeson.Data() << endl;
        PlotM02MergedTrueInPtBins(  fHistoTrueClusMergedM02PtBin, fHistoTrueClusPartConvMergedM02PtBin, fHistoTrueClusPi0M02PtBin, 
                                        fHistoTrueClusPartConvPi0M02PtBin, fHistoTrueClusEtaM02PtBin, fHistoTrueClusPartConvEtaM02PtBin,
                                        nameMeson, nameCanvas, namePad, fMesonM02PlotRange, 
                                        fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem); 

        nameMeson       = Form("%s_MesonFullM02ValidatedDisentangled%s", plotPrefix.Data(), plotSuffix.Data());
//         cout << nameMeson.Data() << endl;
        PlotM02MergedTrueInPtBins(  fHistoTrueClusFullMergedM02PtBin, NULL, fHistoTrueClusFullPi0M02PtBin, 
                                        NULL, fHistoTrueClusFullEtaM02PtBin, NULL,
                                        nameMeson, nameCanvas, namePad, fMesonM02PlotRange, 
                                        fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem);
        if (meson.Contains("Pi0")){
            nameMeson       = Form("%s_MesonPi0ValidatedDisentangled%s", plotPrefix.Data(), plotSuffix.Data());
            cout << nameMeson.Data() << endl;            
            PlotM02MergedTruePrimSecInPtBins(   fHistoTrueClusPi0M02PtBin, fHistoTrueClusPrimPi0M02PtBin, fHistoTrueClusSecPi0M02PtBin, 
                                                fHistoTrueClusSecPi0FK0sM02PtBin, fHistoTrueClusSecPi0FLambdaM02PtBin,
                                                nameMeson, nameCanvas, namePad, fMesonM02PlotRange, 
                                                fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem);

            nameMeson       = Form("%s_MesonPi0PartConvValidatedDisentangled%s", plotPrefix.Data(), plotSuffix.Data());
            cout << nameMeson.Data() << endl;            
            PlotM02MergedTruePrimSecInPtBins(   fHistoTrueClusPartConvPi0M02PtBin, fHistoTrueClusPartConvPrimPi0M02PtBin, fHistoTrueClusPartConvSecPi0M02PtBin, 
                                                fHistoTrueClusPartConvSecPi0FK0sM02PtBin, fHistoTrueClusPartConvSecPi0FLambdaM02PtBin,
                                                nameMeson, nameCanvas, namePad, fMesonM02PlotRange, 
                                                fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem);

            nameMeson       = Form("%s_MesonPi0FullValidatedDisentangled%s", plotPrefix.Data(), plotSuffix.Data());
            cout << nameMeson.Data() << endl;            
            PlotM02MergedTruePrimSecInPtBins(   fHistoTrueClusFullPi0M02PtBin, fHistoTrueClusFullPrimPi0M02PtBin, fHistoTrueClusFullSecPi0M02PtBin, 
                                                fHistoTrueClusFullSecPi0FK0sM02PtBin, fHistoTrueClusFullSecPi0FLambdaM02PtBin,
                                                nameMeson, nameCanvas, namePad, fMesonM02PlotRange, 
                                                fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem);            
        }
    }

    //******************************************************************************************
    //************************ Saving histograms for further processing ************************
    //******************************************************************************************                    
    SaveHistos(fIsMC, fCutSelection, fPrefix2);

    if (fIsMC){
        CalculateMesonAcceptance();
        
        TString fNameHistoEffi      = "TrueMesonEffiMergedPt";
        cout << fNameHistoEffi.Data() << endl;
        fHistoTrueEffiMerged        = CalculateMesonEfficiency(fHistoTrueYieldMergedFullM02, NULL, fNameHistoEffi);

        fNameHistoEffi              = "TrueMesonEffiPrimPt";
        cout << fNameHistoEffi.Data() << endl;
        if (meson.Contains("Pi0")){
            fHistoTrueEffiPrimMeson = CalculateMesonEfficiency(fHistoTrueYieldPrimPi0FullM02, NULL, fNameHistoEffi);
        } else {
            fHistoTrueEffiPrimMeson = CalculateMesonEfficiency(fHistoTrueYieldEtaFullM02, NULL, fNameHistoEffi);
        }    
        
        TString fNameHistoPur       = "TrueMesonPurityMergedPt";
        cout << fNameHistoPur.Data() << endl;
        fHistoTruePurityMerged      = CalculatePurity(fHistoYieldMesonM02, fHistoTrueYieldMergedFullM02, fNameHistoPur);

        fNameHistoPur               = "TruePi0PurityMergedPt";
        cout << fNameHistoPur.Data() << endl;
        fHistoTruePi0PurityMerged   = CalculatePurity(fHistoYieldMesonM02, fHistoTrueYieldPi0FullM02, fNameHistoPur);

        fNameHistoPur               = "TrueEtaPurityMergedPt";
        cout << fNameHistoPur.Data() << endl;
        fHistoTrueEtaPurityMerged   = CalculatePurity(fHistoYieldMesonM02, fHistoTrueYieldEtaFullM02, fNameHistoPur);
        
        // Calculation of secondary fractions
        if (meson.Contains("Pi0")){
            TString fNameHistoFrac      ="TrueSecFrac";
            fHistoTruePi0SecFrac        = CalculateSecondaryFractions(fHistoTrueYieldPi0FullM02, fHistoTrueYieldSecPi0FullM02, fNameHistoFrac);
            fNameHistoFrac              ="TrueSecFracFromK0S";
            fHistoTruePi0SecFracFK0S    = CalculateSecondaryFractions(fHistoTrueYieldPi0FullM02, fHistoTrueYieldSecPi0FK0sFullM02, fNameHistoFrac);
            fNameHistoFrac              ="TrueSecFracFromLambda";
            fHistoTruePi0SecFracFLambda = CalculateSecondaryFractions(fHistoTrueYieldPi0FullM02, fHistoTrueYieldSecPi0FLambdaFullM02, fNameHistoFrac);
        }
        SaveCorrectionHistos(fCutSelection, fPrefix2);   
    }
    //******************************************************************************************
    //************************ Deleting histograms for further processing **********************
    //******************************************************************************************                    
    Delete();
    
    return;
}

//****************************************************************************
//*************** Draw 2D histogram for paper ********************************
//****************************************************************************
void DrawAutoGammaHistoPaper2D( TH2* histo1,
                    TString Title, TString XTitle, TString YTitle,
                    Bool_t YRangeMax, Double_t YMaxFactor, Double_t YMinimum,
                    Bool_t YRange, Double_t YMin ,Double_t YMax,
                    Bool_t XRange, Double_t XMin, Double_t XMax, Double_t xOffset, Double_t yOffset) {
    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR = histo1->GetMinimum();
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR = histo1->GetMinimum();
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }

    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }

    histo1->SetTitle(Title.Data());

    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }

    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);
    histo1->GetYaxis()->SetLabelSize(0.045);
    histo1->GetYaxis()->SetTitleSize(0.05);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetXaxis()->SetTitleOffset(xOffset);
    histo1->GetYaxis()->SetTitleOffset(yOffset);
    histo1->GetXaxis()->SetTitleSize(0.05);
    histo1->GetXaxis()->SetLabelSize(0.045);
}


//**********************************************************************************
//******************* return minimum for 2 D histo  ********************************
//**********************************************************************************
Double_t FindSmallestEntryIn2D(TH2F* histo){
    Double_t minimum = 1;
    for (Int_t i = 1; i<histo->GetNbinsX(); i++){
        for (Int_t j = 1; j<histo->GetNbinsY(); j++){
            if (histo->GetBinContent(i,j) < minimum && histo->GetBinContent(i,j) > 0){
                minimum = histo->GetBinContent(i,j);
            }
        }
    }
    return minimum;
}



//****************************************************************************
//****** Initialization of arrays and variables, histograms for analysis *****
//****** depending on mesonType, number of bins, mode, energy, centrality ****
//****************************************************************************
void Initialize(TString setPi0,Int_t numberOfBins){
    
    InitializeBinning(setPi0, numberOfBins, fEnergyFlag, "", fMode, fEventCutSelection, fClusterMergedCutSelection);
    if (setPi0.CompareTo("Pi0") == 0 || setPi0.CompareTo("Pi0EtaBinning") == 0){
        fMesonId                        = 111;
        fMesonMassPlotRange             = new Double_t[2]; 
        fMesonMassPlotRange[0]          = 0.; 
        fMesonMassPlotRange[1]          = 0.3;
        fMesonMassIntRange              = new Double_t[2]; 
        fMesonMassIntRange[0]           = 0.; 
        fMesonMassIntRange[1]           = 0.3;

        fMesonM02PlotRange              = new Double_t[2];
        fMesonM02PlotRange[0]           = 0;
        fMesonM02PlotRange[1]           = 4.8;
        fMesonM02IntRange               = new Double_t[2];
        fMesonM02IntRange[0]            = 0;
        fMesonM02IntRange[1]            = 4.8;

        if (fNLMmin == 1) {
            fMesonM02PlotRange[1]       = 2.;
        }
    } else if (setPi0.CompareTo("Eta") == 0){
        fMesonId                        = 221;
        fMesonMassPlotRange             = new Double_t[2]; 
        fMesonMassPlotRange[0]          = 0.35; 
        fMesonMassPlotRange[1]          = 0.79;
        fMesonMassIntRange              = new Double_t[2]; 
        fMesonMassIntRange[0]           = 0.35; 
        fMesonMassIntRange[1]           = 0.79;
        fMesonM02PlotRange              = new Double_t[2];
        fMesonM02PlotRange[0]           = 0;
        fMesonM02PlotRange[1]           = 4.5;
        fMesonM02IntRange               = new Double_t[2];
        fMesonM02IntRange[0]            = 0;
        fMesonM02IntRange[1]            = 4.5;
    }
 
    fHistoInvMassPtBin                  = new TH1D*[fNBinsPt];
    fHistoM02PtBin                      = new TH1D*[fNBinsPt];
    fMesonM02Yields                     = new Double_t[fNBinsPt];
    fMesonM02YieldsError                = new Double_t[fNBinsPt];
    
    if (fIsMC){
        fHistoTrueClusMergedInvMassPtBin            = new TH1D*[fNBinsPt];
        fHistoTrueClusMergedM02PtBin                = new TH1D*[fNBinsPt];
        fHistoTrueClusPi0InvMassPtBin               = new TH1D*[fNBinsPt];
        fHistoTrueClusPi0M02PtBin                   = new TH1D*[fNBinsPt];
        fHistoTrueClusEtaInvMassPtBin               = new TH1D*[fNBinsPt];
        fHistoTrueClusEtaM02PtBin                   = new TH1D*[fNBinsPt];
        fHistoTrueClusGammaInvMassPtBin             = new TH1D*[fNBinsPt];
        fHistoTrueClusGammaM02PtBin                 = new TH1D*[fNBinsPt];
        fHistoTrueClusElectronInvMassPtBin          = new TH1D*[fNBinsPt];
        fHistoTrueClusElectronM02PtBin              = new TH1D*[fNBinsPt];
        fHistoTrueClusBGInvMassPtBin                = new TH1D*[fNBinsPt];
        fHistoTrueClusBGM02PtBin                    = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvMergedInvMassPtBin    = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvMergedM02PtBin        = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvPi0InvMassPtBin       = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvPi0M02PtBin           = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvEtaInvMassPtBin       = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvEtaM02PtBin           = new TH1D*[fNBinsPt];
        fHistoTrueClusFullMergedInvMassPtBin        = new TH1D*[fNBinsPt];
        fHistoTrueClusFullMergedM02PtBin            = new TH1D*[fNBinsPt];
        fHistoTrueClusFullPi0InvMassPtBin           = new TH1D*[fNBinsPt];
        fHistoTrueClusFullPi0M02PtBin               = new TH1D*[fNBinsPt];
        fHistoTrueClusFullEtaInvMassPtBin           = new TH1D*[fNBinsPt];
        fHistoTrueClusFullEtaM02PtBin               = new TH1D*[fNBinsPt];
        fHistoTrueClusPrimPi0InvMassPtBin           = new TH1D*[fNBinsPt];
        fHistoTrueClusSecPi0InvMassPtBin            = new TH1D*[fNBinsPt];
        fHistoTrueClusSecPi0FK0sInvMassPtBin        = new TH1D*[fNBinsPt];
        fHistoTrueClusSecPi0FLambdaInvMassPtBin     = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvPrimPi0InvMassPtBin   = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvSecPi0InvMassPtBin    = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvSecPi0FK0sInvMassPtBin    = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvSecPi0FLambdaInvMassPtBin = new TH1D*[fNBinsPt];
        fHistoTrueClusFullPrimPi0InvMassPtBin       = new TH1D*[fNBinsPt];
        fHistoTrueClusFullSecPi0InvMassPtBin        = new TH1D*[fNBinsPt];
        fHistoTrueClusFullSecPi0FK0sInvMassPtBin    = new TH1D*[fNBinsPt];
        fHistoTrueClusFullSecPi0FLambdaInvMassPtBin = new TH1D*[fNBinsPt];
        fHistoTrueClusPrimPi0M02PtBin               = new TH1D*[fNBinsPt];
        fHistoTrueClusSecPi0M02PtBin                = new TH1D*[fNBinsPt];
        fHistoTrueClusSecPi0FK0sM02PtBin            = new TH1D*[fNBinsPt];
        fHistoTrueClusSecPi0FLambdaM02PtBin         = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvPrimPi0M02PtBin       = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvSecPi0M02PtBin        = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvSecPi0FK0sM02PtBin    = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvSecPi0FLambdaM02PtBin = new TH1D*[fNBinsPt];
        fHistoTrueClusFullPrimPi0M02PtBin           = new TH1D*[fNBinsPt];
        fHistoTrueClusFullSecPi0M02PtBin            = new TH1D*[fNBinsPt];
        fHistoTrueClusFullSecPi0FK0sM02PtBin        = new TH1D*[fNBinsPt];
        fHistoTrueClusFullSecPi0FLambdaM02PtBin     = new TH1D*[fNBinsPt];
        
        fMesonM02TrueMergedYields                   = new Double_t[fNBinsPt];
        fMesonM02TruePi0Yields                      = new Double_t[fNBinsPt];
        fMesonM02TrueEtaYields                      = new Double_t[fNBinsPt];
        fMesonM02TrueGammaYields                    = new Double_t[fNBinsPt];
        fMesonM02TrueElectronYields                 = new Double_t[fNBinsPt];
        fMesonM02TrueBGYields                       = new Double_t[fNBinsPt];
        fMesonM02TrueMergedYieldsError              = new Double_t[fNBinsPt];
        fMesonM02TruePi0YieldsError                 = new Double_t[fNBinsPt];
        fMesonM02TrueEtaYieldsError                 = new Double_t[fNBinsPt];
        fMesonM02TrueGammaYieldsError               = new Double_t[fNBinsPt];
        fMesonM02TrueElectronYieldsError            = new Double_t[fNBinsPt];
        fMesonM02TrueBGYieldsError                  = new Double_t[fNBinsPt];
        fMesonM02TrueMergedPartConvYields           = new Double_t[fNBinsPt];
        fMesonM02TruePi0PartConvYields              = new Double_t[fNBinsPt];
        fMesonM02TrueEtaPartConvYields              = new Double_t[fNBinsPt];
        fMesonM02TrueMergedPartConvYieldsError      = new Double_t[fNBinsPt];
        fMesonM02TruePi0PartConvYieldsError         = new Double_t[fNBinsPt];
        fMesonM02TrueEtaPartConvYieldsError         = new Double_t[fNBinsPt];
        fMesonM02TrueMergedFullYields               = new Double_t[fNBinsPt];
        fMesonM02TruePi0FullYields                  = new Double_t[fNBinsPt];
        fMesonM02TrueEtaFullYields                  = new Double_t[fNBinsPt];
        fMesonM02TrueMergedFullYieldsError          = new Double_t[fNBinsPt];
        fMesonM02TruePi0FullYieldsError             = new Double_t[fNBinsPt];
        fMesonM02TrueEtaFullYieldsError             = new Double_t[fNBinsPt];
        fMesonM02TruePrimPi0Yields                  = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0Yields                   = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0FK0sYields               = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0FLambdaYields            = new Double_t[fNBinsPt];
        fMesonM02TruePrimPi0YieldsError             = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0YieldsError              = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0FK0sYieldsError          = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0FLambdaYieldsError       = new Double_t[fNBinsPt];
        fMesonM02TruePrimPi0PartConvYields          = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0PartConvYields           = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0FK0sPartConvYields       = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0FLambdaPartConvYields    = new Double_t[fNBinsPt];
        fMesonM02TruePrimPi0PartConvYieldsError     = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0PartConvYieldsError      = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0FK0sPartConvYieldsError  = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0FLambdaPartConvYieldsError   = new Double_t[fNBinsPt];
        fMesonM02TruePrimPi0FullYields              = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0FullYields               = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0FK0sFullYields           = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0FLambdaFullYields        = new Double_t[fNBinsPt];
        fMesonM02TruePrimPi0FullYieldsError         = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0FullYieldsError          = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0FK0sFullYieldsError      = new Double_t[fNBinsPt];
        fMesonM02TrueSecPi0FLambdaFullYieldsError   = new Double_t[fNBinsPt];
    }
    
    for(Int_t i = 0;i<fNBinsPt; i++){
        fHistoInvMassPtBin[i]               = NULL;
        fHistoM02PtBin[i]                   = NULL;
        fMesonM02Yields[i]                  = 0.;
        fMesonM02YieldsError[i]             = 0.;
        if (fIsMC){
            fHistoTrueClusMergedInvMassPtBin[i]         = NULL;
            fHistoTrueClusMergedM02PtBin[i]             = NULL;
            fHistoTrueClusPi0InvMassPtBin[i]            = NULL;
            fHistoTrueClusPi0M02PtBin[i]                = NULL;
            fHistoTrueClusEtaInvMassPtBin[i]            = NULL;
            fHistoTrueClusEtaM02PtBin[i]                = NULL;
            fHistoTrueClusGammaInvMassPtBin[i]          = NULL;
            fHistoTrueClusGammaM02PtBin[i]              = NULL;
            fHistoTrueClusElectronInvMassPtBin[i]       = NULL;
            fHistoTrueClusElectronM02PtBin[i]           = NULL;
            fHistoTrueClusBGInvMassPtBin[i]             = NULL;
            fHistoTrueClusBGM02PtBin[i]                 = NULL;
            fHistoTrueClusPartConvMergedInvMassPtBin[i] = NULL;
            fHistoTrueClusPartConvMergedM02PtBin[i]     = NULL;
            fHistoTrueClusPartConvPi0InvMassPtBin[i]    = NULL;
            fHistoTrueClusPartConvPi0M02PtBin[i]        = NULL;
            fHistoTrueClusPartConvEtaInvMassPtBin[i]    = NULL;
            fHistoTrueClusPartConvEtaM02PtBin[i]        = NULL;
            fHistoTrueClusFullMergedInvMassPtBin[i]     = NULL;
            fHistoTrueClusFullMergedM02PtBin[i]         = NULL;
            fHistoTrueClusFullPi0InvMassPtBin[i]        = NULL;
            fHistoTrueClusFullPi0M02PtBin[i]            = NULL;
            fHistoTrueClusFullEtaInvMassPtBin[i]        = NULL;
            fHistoTrueClusFullEtaM02PtBin[i]            = NULL;
            fHistoTrueClusPrimPi0InvMassPtBin[i]        = NULL;
            fHistoTrueClusSecPi0InvMassPtBin[i]         = NULL;
            fHistoTrueClusSecPi0FK0sInvMassPtBin[i]     = NULL;
            fHistoTrueClusSecPi0FLambdaInvMassPtBin[i]  = NULL;
            fHistoTrueClusPartConvPrimPi0InvMassPtBin[i]= NULL;
            fHistoTrueClusPartConvSecPi0InvMassPtBin[i] = NULL;
            fHistoTrueClusPartConvSecPi0FK0sInvMassPtBin[i]     = NULL;
            fHistoTrueClusPartConvSecPi0FLambdaInvMassPtBin[i]  = NULL;
            fHistoTrueClusFullPrimPi0InvMassPtBin[i]    = NULL;
            fHistoTrueClusFullSecPi0InvMassPtBin[i]     = NULL;
            fHistoTrueClusFullSecPi0FK0sInvMassPtBin[i] = NULL;
            fHistoTrueClusFullSecPi0FLambdaInvMassPtBin[i]      = NULL;
            fHistoTrueClusPrimPi0M02PtBin[i]            = NULL;
            fHistoTrueClusSecPi0M02PtBin[i]             = NULL;
            fHistoTrueClusSecPi0FK0sM02PtBin[i]         = NULL;
            fHistoTrueClusSecPi0FLambdaM02PtBin[i]      = NULL;
            fHistoTrueClusPartConvPrimPi0M02PtBin[i]    = NULL;
            fHistoTrueClusPartConvSecPi0M02PtBin[i]     = NULL;
            fHistoTrueClusPartConvSecPi0FK0sM02PtBin[i] = NULL;
            fHistoTrueClusPartConvSecPi0FLambdaM02PtBin[i]      = NULL;
            fHistoTrueClusFullPrimPi0M02PtBin[i]        = NULL;
            fHistoTrueClusFullSecPi0M02PtBin[i]         = NULL;
            fHistoTrueClusFullSecPi0FK0sM02PtBin[i]     = NULL;
            fHistoTrueClusFullSecPi0FLambdaM02PtBin[i]  = NULL;

            fMesonM02TrueMergedYields[i]                = 0.;
            fMesonM02TruePi0Yields[i]                   = 0.;
            fMesonM02TrueEtaYields[i]                   = 0.;
            fMesonM02TrueGammaYields[i]                 = 0.;
            fMesonM02TrueElectronYields[i]              = 0.;
            fMesonM02TrueBGYields[i]                    = 0.;
            fMesonM02TrueMergedYieldsError[i]           = 0.;
            fMesonM02TruePi0YieldsError[i]              = 0.;
            fMesonM02TrueEtaYieldsError[i]              = 0.;
            fMesonM02TrueGammaYieldsError[i]            = 0.;
            fMesonM02TrueElectronYieldsError[i]         = 0.;
            fMesonM02TrueBGYieldsError[i]               = 0.;
            fMesonM02TrueMergedPartConvYields[i]        = 0.;
            fMesonM02TruePi0PartConvYields[i]           = 0.;
            fMesonM02TrueEtaPartConvYields[i]           = 0.;
            fMesonM02TrueMergedPartConvYieldsError[i]   = 0.;
            fMesonM02TruePi0PartConvYieldsError[i]      = 0.;
            fMesonM02TrueEtaPartConvYieldsError[i]      = 0.;
            fMesonM02TrueMergedFullYields[i]            = 0.;
            fMesonM02TruePi0FullYields[i]               = 0.;
            fMesonM02TrueEtaFullYields[i]               = 0.;
            fMesonM02TrueMergedFullYieldsError[i]       = 0.;
            fMesonM02TruePi0FullYieldsError[i]          = 0.;
            fMesonM02TrueEtaFullYieldsError[i]          = 0.;
            fMesonM02TruePrimPi0Yields[i]               = 0.;
            fMesonM02TrueSecPi0Yields[i]                = 0.;
            fMesonM02TrueSecPi0FK0sYields[i]            = 0.;
            fMesonM02TrueSecPi0FLambdaYields[i]         = 0.;
            fMesonM02TruePrimPi0YieldsError[i]          = 0.;
            fMesonM02TrueSecPi0YieldsError[i]           = 0.;
            fMesonM02TrueSecPi0FK0sYieldsError[i]       = 0.;
            fMesonM02TrueSecPi0FLambdaYieldsError[i]    = 0.;
            fMesonM02TruePrimPi0PartConvYields[i]       = 0.;
            fMesonM02TrueSecPi0PartConvYields[i]        = 0.;
            fMesonM02TrueSecPi0FK0sPartConvYields[i]    = 0.;
            fMesonM02TrueSecPi0FLambdaPartConvYields[i] = 0.;
            fMesonM02TruePrimPi0PartConvYieldsError[i]  = 0.;
            fMesonM02TrueSecPi0PartConvYieldsError[i]   = 0.;
            fMesonM02TrueSecPi0FK0sPartConvYieldsError[i]       = 0.;
            fMesonM02TrueSecPi0FLambdaPartConvYieldsError[i]    = 0.;
            fMesonM02TruePrimPi0FullYields[i]           = 0.;
            fMesonM02TrueSecPi0FullYields[i]            = 0.;
            fMesonM02TrueSecPi0FK0sFullYields[i]        = 0.;
            fMesonM02TrueSecPi0FLambdaFullYields[i]     = 0.;
            fMesonM02TruePrimPi0FullYieldsError[i]      = 0.;
            fMesonM02TrueSecPi0FullYieldsError[i]       = 0.;
            fMesonM02TrueSecPi0FK0sFullYieldsError[i]   = 0.;
            fMesonM02TrueSecPi0FLambdaFullYieldsError[i]        = 0.;
            
        }
    }
}

//****************************************************************************
//************** Initializiation of MC histogram names according *************
//****************************************************************************
void SetCorrectMCHistogrammNames(TString mesonType){

    fObjectNameTrueMergedM02                    = "ESD_TrueClusMerged_Pt_M02";
    fObjectNameTrueMergedInvMass                = "ESD_TrueClusMerged_InvMass_Pt";
    fObjectNameTrueFromPi0M02                   = "ESD_TrueClusFromPi0_Pt_M02";
    fObjectNameTrueFromPi0InvMass               = "ESD_TrueClusFromPi0_InvMass_Pt";
    fObjectNameTrueFromEtaM02                   = "ESD_TrueClusFromEta_Pt_M02";
    fObjectNameTrueFromEtaInvMass               = "ESD_TrueClusFromEta_InvMass_Pt";
    fObjectNameTrueMergedPartConvM02            = "ESD_TrueClusMergedPartConv_Pt_M02";
    fObjectNameTrueMergedPartConvLeadEM02       = "ESD_TrueClusMergedPartConvLeadE_Pt_M02";
    fObjectNameTrueMergedPartConvLeadEInvMass   = "ESD_TrueClusMergedPartConvLeadE_InvMass_Pt";
    fObjectNameTrueMergedPartConvInvMass        = "ESD_TrueClusMergedPartConv_InvMass_Pt";
    fObjectNameTrueClusPartConvFromPi0M02       = "ESD_TrueClusPartConvFromPi0_Pt_M02";
    fObjectNameTrueClusPartConvFromPi0InvMass   = "ESD_TrueClusPartConvFromPi0_InvMass_Pt";
    fObjectNameTrueClusPartConvFromEtaM02       = "ESD_TrueClusPartConvFromEta_Pt_M02";
    fObjectNameTrueClusPartConvFromEtaInvMass   = "ESD_TrueClusPartConvFromEta_InvMass_Pt";
    fObjectNameTrueClusBGM02                    = "ESD_TrueClusBG_Pt_M02";
    fObjectNameTrueClusBGInvMass                = "ESD_TrueClusBG_InvMass_Pt";
    fObjectNameTrueClusGammaM02                 = "ESD_TrueClusGamma_Pt_M02";
    fObjectNameTrueClusGammaInvMass             = "ESD_TrueClusGamma_InvMass_Pt";
    fObjectNameTrueClusElectronM02              = "ESD_TrueClusElectron_Pt_M02";
    fObjectNameTrueClusElectronInvMass          = "ESD_TrueClusElectron_InvMass_Pt";
    fObjectNameTrueClusBG_Source                = "ESD_TrueClusBG_Pt_Source";
    if (mesonType.Contains("Pi0")){
        fObjectNameMCMesonAcc                                   = "MC_Pi0InAcc_Pt";
        fObjectNameMCMeson                                      = "MC_Pi0_Pt";
        fObjectNameMCMesonWOWeights                             = "MC_Pi0_WOWeights_Pt";
        fObjectNameTrueClusPrimMesonM02                         = "ESD_TrueClusFromPrimPi0_Pt_M02";
        fObjectNameTrueClusPrimMesonInvMass                     = "ESD_TrueClusFromPrimPi0_InvMass_Pt";
        fObjectNameTrueClusSecMesonM02                          = "ESD_TrueClusFromSecPi0_Pt_M02";
        fObjectNameTrueClusSecMesonInvMass                      = "ESD_TrueClusFromSecPi0_InvMass_Pt";
        fObjectNameTrueClusSecMesonFromK0sM02                   = "ESD_TrueClusFromSecPi0FromK0s_Pt_M02";
        fObjectNameTrueClusSecMesonFromK0sInvMass               = "ESD_TrueClusFromSecPi0FromK0s_InvMass_Pt";
        fObjectNameTrueClusSecMesonFromLambdaM02                = "ESD_TrueClusFromSecPi0FromLambda_Pt_M02";
        fObjectNameTrueClusSecMesonFromLambdaInvMass            = "ESD_TrueClusFromSecPi0FromLambda_InvMass_Pt";
        fObjectNameTrueClusPartConvPrimMesonM02                 = "ESD_TrueClusPartConvFromPrimPi0_Pt_M02";
        fObjectNameTrueClusPartConvPrimMesonInvMass             = "ESD_TrueClusPartConvFromPrimPi0_InvMass_Pt";
        fObjectNameTrueClusPartConvSecMesonM02                  = "ESD_TrueClusPartConvFromSecPi0_Pt_M02";
        fObjectNameTrueClusPartConvSecMesonInvMass              = "ESD_TrueClusPartConvFromSecPi0_InvMass_Pt";
        fObjectNameTrueClusPartConvSecMesonFromK0sM02           = "ESD_TrueClusPartConvFromSecPi0FromK0s_Pt_M02";
        fObjectNameTrueClusPartConvSecMesonFromK0sInvMass       = "ESD_TrueClusPartConvFromSecPi0FromK0s_InvMass_Pt";
        fObjectNameTrueClusPartConvSecMesonFromLambdaM02        = "ESD_TrueClusPartConvFromSecPi0FromLamdba_Pt_M02";
        fObjectNameTrueClusPartConvSecMesonFromLambdaInvMass    = "ESD_TrueClusPartConvFromSecPi0FromLamdba_InvMass_Pt";
    } else {
        fObjectNameMCMesonAcc                       = "MC_EtaInAcc_Pt";
        fObjectNameMCMeson                          = "MC_Eta_Pt";
        fObjectNameMCMesonWOWeights                 = "MC_Eta_WOWeights_Pt";
    }
    return;
}

//****************************************************************************
//*********************** creation of momentum dependent histos **************
//****************************************************************************
void CreatePtHistos(){

    fDeltaPt                            = new TH1D("deltaPt", "", fNBinsPt, fBinsPt);
    fDeltaPt->Sumw2();

    fHistoYieldMesonM02                 = new TH1D("histoYieldMesonM02", "", fNBinsPt, fBinsPt);
    fHistoYieldMesonM02->Sumw2();
    
    if (fIsMC){
        fHistoTrueYieldMergedM02                  = new TH1D("histoTrueYieldMergedM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedM02->Sumw2();
        fHistoTrueYieldMergedPartConvM02          = new TH1D("histoTrueYieldMergedPartConvM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedPartConvM02->Sumw2();
        fHistoTrueYieldMergedFullM02              = new TH1D("histoTrueYieldMergedFullM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedFullM02->Sumw2();
        fHistoTrueYieldPi0M02                     = new TH1D("histoTrueYieldPi0M02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldPi0M02->Sumw2();
        fHistoTrueYieldPi0PartConvM02             = new TH1D("histoTrueYieldPi0PartConvM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldPi0PartConvM02->Sumw2();
        fHistoTrueYieldPi0FullM02                 = new TH1D("histoTrueYieldPi0FullM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldPi0FullM02->Sumw2();
        fHistoTrueYieldEtaM02                     = new TH1D("histoTrueYieldEtaM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldEtaM02->Sumw2();
        fHistoTrueYieldEtaPartConvM02             = new TH1D("histoTrueYieldEtaPartConvM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldEtaPartConvM02->Sumw2();
        fHistoTrueYieldEtaFullM02                 = new TH1D("histoTrueYieldEtaFullM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldEtaFullM02->Sumw2();
        fHistoTrueYieldGammaM02                   = new TH1D("histoTrueYieldGammaM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldGammaM02->Sumw2();
        fHistoTrueYieldElectronM02                = new TH1D("histoTrueYieldElectronM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldElectronM02->Sumw2();
        fHistoTrueYieldBGM02                      = new TH1D("histoTrueYieldBGM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldBGM02->Sumw2();
        if (fPrefix.Contains("Pi0")){
            fHistoTrueYieldPrimPi0M02                   = new TH1D("histoTrueYieldPrimPi0M02", "", fNBinsPt, fBinsPt);
            fHistoTrueYieldPrimPi0M02->Sumw2();
            fHistoTrueYieldPrimPi0PartConvM02           = new TH1D("histoTrueYieldPrimPi0PartConvM02", "", fNBinsPt, fBinsPt);
            fHistoTrueYieldPrimPi0PartConvM02->Sumw2();
            fHistoTrueYieldPrimPi0FullM02               = new TH1D("histoTrueYieldPrimPi0FullM02", "", fNBinsPt, fBinsPt);
            fHistoTrueYieldPrimPi0FullM02->Sumw2();
            fHistoTrueYieldSecPi0M02                    = new TH1D("histoTrueYieldSecPi0M02", "", fNBinsPt, fBinsPt);
            fHistoTrueYieldSecPi0M02->Sumw2();
            fHistoTrueYieldSecPi0PartConvM02            = new TH1D("histoTrueYieldSecPi0PartConvM02", "", fNBinsPt, fBinsPt);
            fHistoTrueYieldSecPi0PartConvM02->Sumw2();
            fHistoTrueYieldSecPi0FullM02                = new TH1D("histoTrueYieldSecPi0FullM02", "", fNBinsPt, fBinsPt);
            fHistoTrueYieldSecPi0FullM02->Sumw2();
            fHistoTrueYieldSecPi0FK0sM02                = new TH1D("histoTrueYieldSecPi0FK0sM02", "", fNBinsPt, fBinsPt);
            fHistoTrueYieldSecPi0FK0sM02->Sumw2();
            fHistoTrueYieldSecPi0FK0sPartConvM02        = new TH1D("histoTrueYieldSecPi0FK0sPartConvM02", "", fNBinsPt, fBinsPt);
            fHistoTrueYieldSecPi0FK0sPartConvM02->Sumw2();
            fHistoTrueYieldSecPi0FK0sFullM02            = new TH1D("histoTrueYieldSecPi0FK0sFullM02", "", fNBinsPt, fBinsPt);
            fHistoTrueYieldSecPi0FK0sFullM02->Sumw2();
            fHistoTrueYieldSecPi0FLambdaM02             = new TH1D("histoTrueYieldSecPi0FLambdaM02", "", fNBinsPt, fBinsPt);
            fHistoTrueYieldSecPi0FLambdaM02->Sumw2();
            fHistoTrueYieldSecPi0FLambdaPartConvM02     = new TH1D("histoTrueYieldSecPi0FLambdaPartConvM02", "", fNBinsPt, fBinsPt);
            fHistoTrueYieldSecPi0FLambdaPartConvM02->Sumw2();
            fHistoTrueYieldSecPi0FLambdaFullM02         = new TH1D("histoTrueYieldSecPi0FLambdaFullM02", "", fNBinsPt, fBinsPt);
            fHistoTrueYieldSecPi0FLambdaFullM02->Sumw2();
        }
    }
}


//****************************************************************************
//*************** Fill momentum dependent histograms from arrays *************
//****************************************************************************
void FillPtHistos(){
    for(Int_t iPt=fStartPtBin+1;iPt<fNBinsPt+1;iPt++){
        fDeltaPt->SetBinContent(iPt,fBinsPt[iPt]-fBinsPt[iPt-1]);
        fDeltaPt->SetBinError(iPt,0);
        
        fHistoYieldMesonM02->SetBinContent(iPt, fMesonM02Yields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonM02->SetBinError(iPt, fMesonM02YieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        
        if (fIsMC){
            fHistoTrueYieldMergedM02->SetBinContent(iPt, fMesonM02TrueMergedYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedM02->SetBinError(iPt, fMesonM02TrueMergedYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldMergedPartConvM02->SetBinContent(iPt, fMesonM02TrueMergedPartConvYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedPartConvM02->SetBinError(iPt, fMesonM02TrueMergedPartConvYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            
            fHistoTrueYieldMergedFullM02->SetBinContent(iPt, fMesonM02TrueMergedFullYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedFullM02->SetBinError(iPt, fMesonM02TrueMergedFullYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldPi0M02->SetBinContent(iPt, fMesonM02TruePi0Yields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldPi0M02->SetBinError(iPt, fMesonM02TruePi0YieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldPi0PartConvM02->SetBinContent(iPt, fMesonM02TruePi0PartConvYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldPi0PartConvM02->SetBinError(iPt, fMesonM02TruePi0PartConvYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            
            fHistoTrueYieldPi0FullM02->SetBinContent(iPt, fMesonM02TruePi0FullYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldPi0FullM02->SetBinError(iPt, fMesonM02TruePi0FullYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldEtaM02->SetBinContent(iPt, fMesonM02TrueEtaYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldEtaM02->SetBinError(iPt, fMesonM02TrueEtaYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldEtaPartConvM02->SetBinContent(iPt, fMesonM02TrueEtaPartConvYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldEtaPartConvM02->SetBinError(iPt, fMesonM02TrueEtaPartConvYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            
            fHistoTrueYieldEtaFullM02->SetBinContent(iPt, fMesonM02TrueEtaFullYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldEtaFullM02->SetBinError(iPt, fMesonM02TrueEtaFullYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldGammaM02->SetBinContent(iPt, fMesonM02TrueGammaYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldGammaM02->SetBinError(iPt, fMesonM02TrueGammaYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldElectronM02->SetBinContent(iPt, fMesonM02TrueElectronYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldElectronM02->SetBinError(iPt, fMesonM02TrueElectronYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldBGM02->SetBinContent(iPt, fMesonM02TrueBGYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldBGM02->SetBinError(iPt, fMesonM02TrueBGYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            
            if (fPrefix.Contains("Pi0")){
                fHistoTrueYieldPrimPi0M02->SetBinContent(iPt, fMesonM02TruePrimPi0Yields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoTrueYieldPrimPi0M02->SetBinError(iPt, fMesonM02TruePrimPi0YieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

                fHistoTrueYieldPrimPi0PartConvM02->SetBinContent(iPt, fMesonM02TruePrimPi0PartConvYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoTrueYieldPrimPi0PartConvM02->SetBinError(iPt, fMesonM02TruePrimPi0PartConvYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                
                fHistoTrueYieldPrimPi0FullM02->SetBinContent(iPt, fMesonM02TruePrimPi0FullYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoTrueYieldPrimPi0FullM02->SetBinError(iPt, fMesonM02TruePrimPi0FullYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

                fHistoTrueYieldSecPi0M02->SetBinContent(iPt, fMesonM02TrueSecPi0Yields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoTrueYieldSecPi0M02->SetBinError(iPt, fMesonM02TrueSecPi0YieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

                fHistoTrueYieldSecPi0PartConvM02->SetBinContent(iPt, fMesonM02TrueSecPi0PartConvYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoTrueYieldSecPi0PartConvM02->SetBinError(iPt, fMesonM02TrueSecPi0PartConvYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                
                fHistoTrueYieldSecPi0FullM02->SetBinContent(iPt, fMesonM02TrueSecPi0FullYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoTrueYieldSecPi0FullM02->SetBinError(iPt, fMesonM02TrueSecPi0FullYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

                fHistoTrueYieldSecPi0FK0sM02->SetBinContent(iPt, fMesonM02TrueSecPi0FK0sYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoTrueYieldSecPi0FK0sM02->SetBinError(iPt, fMesonM02TrueSecPi0FK0sYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

                fHistoTrueYieldSecPi0FK0sPartConvM02->SetBinContent(iPt, fMesonM02TrueSecPi0FK0sPartConvYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoTrueYieldSecPi0FK0sPartConvM02->SetBinError(iPt, fMesonM02TrueSecPi0FK0sPartConvYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                
                fHistoTrueYieldSecPi0FK0sFullM02->SetBinContent(iPt, fMesonM02TrueSecPi0FK0sFullYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoTrueYieldSecPi0FK0sFullM02->SetBinError(iPt, fMesonM02TrueSecPi0FK0sFullYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

                fHistoTrueYieldSecPi0FLambdaM02->SetBinContent(iPt, fMesonM02TrueSecPi0FLambdaYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoTrueYieldSecPi0FLambdaM02->SetBinError(iPt, fMesonM02TrueSecPi0FLambdaYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

                fHistoTrueYieldSecPi0FLambdaPartConvM02->SetBinContent(iPt, fMesonM02TrueSecPi0FLambdaPartConvYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoTrueYieldSecPi0FLambdaPartConvM02->SetBinError(iPt, fMesonM02TrueSecPi0FLambdaPartConvYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                
                fHistoTrueYieldSecPi0FLambdaFullM02->SetBinContent(iPt, fMesonM02TrueSecPi0FLambdaFullYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoTrueYieldSecPi0FLambdaFullM02->SetBinError(iPt, fMesonM02TrueSecPi0FLambdaFullYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                
            }    
        }
    }
}

//****************************************************************************
//******************** Projection out of 2D in X *****************************
//****************************************************************************
TH1D* FillProjectionX (TH2F* fDummy2D, TString name, Double_t minY, Double_t maxY, Int_t rebin){
    TH1D* dummy1D           = new TH1D(name.Data(), name.Data(), fDummy2D->GetNbinsX(), 0., fDummy2D->GetXaxis()->GetBinUpEdge(fDummy2D->GetNbinsX()));
    dummy1D->Sumw2();
    Int_t startBin          = fDummy2D->GetYaxis()->FindBin(minY+0.001);
    Int_t endBin            = fDummy2D->GetYaxis()->FindBin(maxY-0.001);
    fDummy2D->ProjectionX(name.Data(),startBin,endBin,"e");
    dummy1D                 = (TH1D*)gDirectory->Get(name.Data());
    if(rebin>1){
        dummy1D->Rebin(rebin);
    }
    return dummy1D; 
}

//****************************************************************************
//******************** Projection out of 2D in X *****************************
//****************************************************************************
TH1D* FillProjectionY (TH2F* fDummy2D, TString name, Double_t minX, Double_t maxX, Int_t rebin){
    TH1D* dummy1D           = new TH1D(name.Data(), name.Data(), fDummy2D->GetNbinsY(), 0., fDummy2D->GetYaxis()->GetBinUpEdge(fDummy2D->GetNbinsY()));
    dummy1D->Sumw2();
    Int_t startBin          = fDummy2D->GetXaxis()->FindBin(minX+0.001);
    Int_t endBin            = fDummy2D->GetXaxis()->FindBin(maxX-0.001);
    fDummy2D->ProjectionY(name.Data(),startBin,endBin,"e");
    dummy1D                 = (TH1D*)gDirectory->Get(name.Data());
    if(rebin>1){
        dummy1D->Rebin(rebin);
    }
    return dummy1D; 
}

//****************************************************************************
//******************** Add possibly not existing histos **********************
//****************************************************************************
TH1D* AddPossiblyNotExistingHists(TH1D* fDummy1, TH1D* fDummy2, TString name){
    TH1D* fDummy3 = NULL;
    if ( fDummy1!= NULL && fDummy2!= NULL ){
        fDummy3       = (TH1D*)fDummy1->Clone(name.Data());
        fDummy3->Sumw2();
        fDummy3->Add(fDummy2);
    } else if (fDummy1!= NULL) {
        fDummy3       = (TH1D*)fDummy1->Clone(name.Data());
        fDummy3->Sumw2();
    } else if ( fDummy2!= NULL ) {
        fDummy3       = (TH1D*)fDummy2->Clone(name.Data());
        fDummy3->Sumw2();
    }
    return fDummy3; 
}

//****************************************************************************
//*************** Check if histo already exists if not clear *****************
//****************************************************************************
void CheckForNULLForPointer(TH1D* fDummy1){
    if(fDummy1!= NULL){
        delete fDummy1;
        fDummy1           = NULL;
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//****************************************************************************
void FillDataHistosArray(TH2F* fInvMassVSPtDummy, TH2F* fM02VsPtDummy) {
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        TString fNameHistoInvMass       = Form("InvMass_in_Pt_Bin%02d", iPt);
        TString fNameHistoM02           = Form("M02_in_Pt_Bin%02d", iPt);

        CheckForNULLForPointer(fHistoInvMassPtBin[iPt]);
        if (fInvMassVSPtDummy){
            fHistoInvMassPtBin[iPt]     = FillProjectionX(fInvMassVSPtDummy, fNameHistoInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }

        CheckForNULLForPointer(fHistoM02PtBin[iPt]);
        if (fM02VsPtDummy){
            fHistoM02PtBin[iPt]         = FillProjectionY(fM02VsPtDummy, fNameHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
    }
}

//****************************************************************************
//**** Fill array of invariant mass histograms in pT slices for MC ***********
//****************************************************************************
void FillMCInvMassHistosArray(  TH2F* fTrueMergedInvMassVsPtDummy, 
                                TH2F* fTruePi0InvMassVsPt, 
                                TH2F* fTrueEtaInvMassVsPt,
                                TH2F* fTrueGammaInvMassVsPt,
                                TH2F* fTrueElectronInvMassVsPt,
                                TH2F* fTrueBGInvMassVsPt,
                                TH2F* fTruePartConvMergedInvMassVsPtDummy, 
                                TH2F* fTruePartConvPi0InvMassVsPt, 
                                TH2F* fTruePartConvEtaInvMassVsPt
                              ) {

    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        TString fNameHistoInvMass           = Form("TrueClusMerged_InvMass_in_Pt_Bin%02d", iPt);
        TString fNamePi0HistoInvMass        = Form("TrueClusPi0_InvMass_in_Pt_Bin%02d", iPt);
        TString fNameEtaHistoInvMass        = Form("TrueClusEta_InvMass_in_Pt_Bin%02d", iPt);
        TString fNameGammaHistoInvMass      = Form("TrueClusGamma_InvMass_in_Pt_Bin%02d", iPt);
        TString fNameElectronHistoInvMass   = Form("TrueClusElectron_InvMass_in_Pt_Bin%02d", iPt);
        TString fNameBGHistoInvMass         = Form("TrueClusBG_InvMass_in_Pt_Bin%02d", iPt);
        TString fNamePartConvHistoInvMass   = Form("TrueClusPartConvMerged_InvMass_in_Pt_Bin%02d", iPt);
        TString fNamePartConvPi0HistoInvMass= Form("TrueClusPartConvPi0_InvMass_in_Pt_Bin%02d", iPt);
        TString fNamePartConvEtaHistoInvMass= Form("TrueClusPartConvEta_InvMass_in_Pt_Bin%02d", iPt);
        
        // true merged clusters
        CheckForNULLForPointer(fHistoTrueClusMergedInvMassPtBin[iPt]);
        if (fTrueMergedInvMassVsPtDummy){
            fHistoTrueClusMergedInvMassPtBin[iPt]               = FillProjectionX(fTrueMergedInvMassVsPtDummy, fNameHistoInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // true merged clusters part conv
        CheckForNULLForPointer(fHistoTrueClusPartConvMergedInvMassPtBin[iPt]);
        if (fTruePartConvMergedInvMassVsPtDummy){
            fHistoTrueClusPartConvMergedInvMassPtBin[iPt]       = FillProjectionX(fTruePartConvMergedInvMassVsPtDummy, fNamePartConvHistoInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // adding merged cluster histograms
        CheckForNULLForPointer(fHistoTrueClusFullMergedInvMassPtBin[iPt]);
        fHistoTrueClusFullMergedInvMassPtBin[iPt]               = AddPossiblyNotExistingHists( fHistoTrueClusMergedInvMassPtBin[iPt], fHistoTrueClusPartConvMergedInvMassPtBin[iPt],
                                                                                                Form("TrueClusMergedFull_InvMass_in_Pt_Bin%02d", iPt));
        
        // true pi0 clusters
        CheckForNULLForPointer(fHistoTrueClusPi0InvMassPtBin[iPt]);
        if (fTruePi0InvMassVsPt){
            fHistoTrueClusPi0InvMassPtBin[iPt]                  = FillProjectionX(fTruePi0InvMassVsPt, fNamePi0HistoInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // true pi0 clusters part conv
        CheckForNULLForPointer(fHistoTrueClusPartConvPi0InvMassPtBin[iPt]);
        if (fTruePartConvPi0InvMassVsPt){
            fHistoTrueClusPartConvPi0InvMassPtBin[iPt]          = FillProjectionX(fTruePartConvPi0InvMassVsPt, fNamePartConvPi0HistoInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // adding pi0 cluster histograms
        CheckForNULLForPointer(fHistoTrueClusFullPi0InvMassPtBin[iPt]);
        fHistoTrueClusFullPi0InvMassPtBin[iPt]                  = AddPossiblyNotExistingHists( fHistoTrueClusPi0InvMassPtBin[iPt], fHistoTrueClusPartConvPi0InvMassPtBin[iPt],
                                                                                                Form("TrueClusFullPi0_InvMass_in_Pt_Bin%02d", iPt));
        
        // true eta clusters
        CheckForNULLForPointer(fHistoTrueClusEtaInvMassPtBin[iPt]);
        if (fTrueEtaInvMassVsPt){
            fHistoTrueClusEtaInvMassPtBin[iPt]                  = FillProjectionX(fTrueEtaInvMassVsPt, fNameEtaHistoInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // true eta clusters part conv
        CheckForNULLForPointer(fHistoTrueClusPartConvEtaInvMassPtBin[iPt]);
        if (fTruePartConvEtaInvMassVsPt){
            fHistoTrueClusPartConvEtaInvMassPtBin[iPt]          = FillProjectionX(fTruePartConvEtaInvMassVsPt, fNamePartConvEtaHistoInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // adding eta cluster histograms
        CheckForNULLForPointer(fHistoTrueClusFullEtaInvMassPtBin[iPt]);
        fHistoTrueClusFullEtaInvMassPtBin[iPt]                  = AddPossiblyNotExistingHists( fHistoTrueClusEtaInvMassPtBin[iPt], fHistoTrueClusPartConvEtaInvMassPtBin[iPt],
                                                                                                Form("TrueClusFullEta_InvMass_in_Pt_Bin%02d", iPt));
        
        // true gamma clusters
        CheckForNULLForPointer(fHistoTrueClusGammaInvMassPtBin[iPt]);
        if (fTrueGammaInvMassVsPt){
            fHistoTrueClusGammaInvMassPtBin[iPt]                = FillProjectionX(fTrueGammaInvMassVsPt, fNameGammaHistoInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // true electron clusters
        CheckForNULLForPointer(fHistoTrueClusElectronInvMassPtBin[iPt]);
        if (fTrueElectronInvMassVsPt){
            fHistoTrueClusElectronInvMassPtBin[iPt]             = FillProjectionX(fTrueElectronInvMassVsPt, fNameElectronHistoInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // true BG clusters
        CheckForNULLForPointer(fHistoTrueClusBGInvMassPtBin[iPt]);
        if (fTrueBGInvMassVsPt){
            fHistoTrueClusBGInvMassPtBin[iPt]                   = FillProjectionX(fTrueBGInvMassVsPt, fNameBGHistoInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
    }
}


//****************************************************************************
//** Fill array of invariant mass histograms in pT slices for Pi0 Prim/Sec ***
//****************************************************************************
void FillMCPrimSecInvMassHistosArray(   TH2F* fTruePi0PrimInvMassVsPt, 
                                        TH2F* fTruePi0SecInvMassVsPt, 
                                        TH2F* fTruePi0SecFromK0sInvMassVsPt,
                                        TH2F* fTruePi0SecFromLambdaInvMassVsPt,
                                        TH2F* fTruePi0PrimPartConvInvMassVsPt, 
                                        TH2F* fTruePi0SecPartConvInvMassVsPt, 
                                        TH2F* fTruePi0SecFromK0sPartConvInvMassVsPt,
                                        TH2F* fTruePi0SecFromLambdaPartConvInvMassVsPt
                              ) {

    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        TString fNamePi0PrimInvMass                 = Form("TrueClusPrimPi0_InvMass_in_Pt_Bin%02d", iPt);
        TString fNamePi0SecInvMass                  = Form("TrueClusSecPi0_InvMass_in_Pt_Bin%02d", iPt);
        TString fNameSecPi0FK0sInvMass              = Form("TrueClusSecPi0FromK0s_InvMass_in_Pt_Bin%02d", iPt);
        TString fNameSecPi0FLambdaInvMass           = Form("TrueClusSecPi0FromLambda_InvMass_in_Pt_Bin%02d", iPt);
        TString fNamePi0PrimPartConvInvMass         = Form("TrueClusPartConvPrimPi0_InvMass_in_Pt_Bin%02d", iPt);
        TString fNamePi0SecPartConvInvMass          = Form("TrueClusPartConvSecPi0_InvMass_in_Pt_Bin%02d", iPt);
        TString fNameSecPi0FK0sPartConvInvMass      = Form("TrueClusPartConvSecPi0FromK0s_InvMass_in_Pt_Bin%02d", iPt);
        TString fNameSecPi0FLambdaPartConvInvMass   = Form("TrueClusPartConvSecPi0FromLambda_InvMass_in_Pt_Bin%02d", iPt);
        
        // true prim pi0 clusters
        CheckForNULLForPointer(fHistoTrueClusPrimPi0InvMassPtBin[iPt]);
        if (fTruePi0PrimInvMassVsPt){
            fHistoTrueClusPrimPi0InvMassPtBin[iPt]              = FillProjectionX(fTruePi0PrimInvMassVsPt, fNamePi0PrimInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // true pi0 prim clusters part conv
        CheckForNULLForPointer(fHistoTrueClusPartConvPrimPi0InvMassPtBin[iPt]);
        if (fTruePi0PrimPartConvInvMassVsPt){
            fHistoTrueClusPartConvPrimPi0InvMassPtBin[iPt]      = FillProjectionX(fTruePi0PrimPartConvInvMassVsPt, fNamePi0PrimPartConvInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // adding prim pi0 cluster histograms
        CheckForNULLForPointer(fHistoTrueClusFullPrimPi0InvMassPtBin[iPt]);
        fHistoTrueClusFullPrimPi0InvMassPtBin[iPt]              = AddPossiblyNotExistingHists( fHistoTrueClusPrimPi0InvMassPtBin[iPt], fHistoTrueClusPartConvPrimPi0InvMassPtBin[iPt],
                                                                                               Form("TrueClusFullPi0Prim_InvMass_in_Pt_Bin%02d", iPt));

        // true pi0 clusters from sec
        CheckForNULLForPointer(fHistoTrueClusSecPi0InvMassPtBin[iPt]);
        if (fTruePi0SecInvMassVsPt){
            fHistoTrueClusSecPi0InvMassPtBin[iPt]               = FillProjectionX(fTruePi0SecInvMassVsPt, fNamePi0SecInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // true pi0 from sec clusters part conv
        CheckForNULLForPointer(fHistoTrueClusPartConvSecPi0InvMassPtBin[iPt]);
        if (fTruePi0SecPartConvInvMassVsPt){
            fHistoTrueClusPartConvSecPi0InvMassPtBin[iPt]       = FillProjectionX(fTruePi0SecPartConvInvMassVsPt, fNamePi0SecPartConvInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // adding pi0 from sec cluster histograms
        CheckForNULLForPointer(fHistoTrueClusFullSecPi0InvMassPtBin[iPt]);
        fHistoTrueClusFullSecPi0InvMassPtBin[iPt]               = AddPossiblyNotExistingHists( fHistoTrueClusSecPi0InvMassPtBin[iPt], fHistoTrueClusPartConvSecPi0InvMassPtBin[iPt],
                                                                                               Form("TrueClusFullPi0Sec_InvMass_in_Pt_Bin%02d", iPt));

        // true pi0 from sec from K0s clusters
        CheckForNULLForPointer(fHistoTrueClusSecPi0FK0sInvMassPtBin[iPt]);
        if (fTruePi0SecFromK0sInvMassVsPt){
            fHistoTrueClusSecPi0FK0sInvMassPtBin[iPt]           = FillProjectionX(fTruePi0SecFromK0sInvMassVsPt, fNameSecPi0FK0sInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // true pi0 from sec from K0s clusters part conv
        CheckForNULLForPointer(fHistoTrueClusPartConvSecPi0FK0sInvMassPtBin[iPt]);
        if (fTruePi0SecFromK0sPartConvInvMassVsPt){
            fHistoTrueClusPartConvSecPi0FK0sInvMassPtBin[iPt]   = FillProjectionX(fTruePi0SecFromK0sPartConvInvMassVsPt, fNameSecPi0FK0sPartConvInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // adding pi0 from sec from K0s cluster histograms
        CheckForNULLForPointer(fHistoTrueClusFullSecPi0FK0sInvMassPtBin[iPt]);
        fHistoTrueClusFullSecPi0FK0sInvMassPtBin[iPt]           = AddPossiblyNotExistingHists( fHistoTrueClusSecPi0FK0sInvMassPtBin[iPt], fHistoTrueClusPartConvSecPi0FK0sInvMassPtBin[iPt],
                                                                                               Form("TrueClusFullPi0SecFK0s_InvMass_in_Pt_Bin%02d", iPt));

        // true pi0 from sec from Lambda clusters
        CheckForNULLForPointer(fHistoTrueClusSecPi0FLambdaInvMassPtBin[iPt]);
        if (fTruePi0SecFromLambdaInvMassVsPt){
            fHistoTrueClusSecPi0FLambdaInvMassPtBin[iPt]        = FillProjectionX(fTruePi0SecFromLambdaInvMassVsPt, fNameSecPi0FLambdaInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // true pi0 from sec from Lambda clusters part conv
        CheckForNULLForPointer(fHistoTrueClusPartConvSecPi0FLambdaInvMassPtBin[iPt]);
        if (fTruePi0SecFromLambdaPartConvInvMassVsPt){
            fHistoTrueClusPartConvSecPi0FLambdaInvMassPtBin[iPt]= FillProjectionX(fTruePi0SecFromLambdaPartConvInvMassVsPt, fNameSecPi0FLambdaPartConvInvMass, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        }
        // adding pi0 from sec from Lambda cluster histograms
        CheckForNULLForPointer(fHistoTrueClusFullSecPi0FLambdaInvMassPtBin[iPt]);
        fHistoTrueClusFullSecPi0FLambdaInvMassPtBin[iPt]        = AddPossiblyNotExistingHists( fHistoTrueClusSecPi0FLambdaInvMassPtBin[iPt], fHistoTrueClusPartConvSecPi0FLambdaInvMassPtBin[iPt],
                                                                                               Form("TrueClusFullPi0SecFLambda_InvMass_in_Pt_Bin%02d", iPt));
    }
}    

//****************************************************************************
//*************** Fill array of M02 histograms in pT slices for MC ***********
//****************************************************************************
void FillMCM02HistosArray(  TH2F* fTrueMergedPtVsM02Dummy, 
                            TH2F* fTruePi0PtVsM02, 
                            TH2F* fTrueEtaPtVsM02,
                            TH2F* fTrueGammaPtVsM02,
                            TH2F* fTrueElectronPtVsM02,
                            TH2F* fTrueBGPtVsM02,
                            TH2F* fTruePartConvMergedVsM02Dummy,
                            TH2F* fTruePartConvPi0PtVsM02, 
                            TH2F* fTruePartConvEtaPtVsM02
                         ) {
    cout << __LINE__ << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        TString fNameHistoM02                           = Form("TrueClusMerged_M02_in_Pt_Bin%02d", iPt);
        TString fNamePi0HistoM02                        = Form("TrueClusPi0_M02_in_Pt_Bin%02d", iPt);
        TString fNameEtaHistoM02                        = Form("TrueClusEta_M02_in_Pt_Bin%02d", iPt);
        TString fNameGammaHistoM02                      = Form("TrueClusGamma_M02_in_Pt_Bin%02d", iPt);
        TString fNameElectronHistoM02                   = Form("TrueClusElectron_M02_in_Pt_Bin%02d", iPt);
        TString fNameBGHistoM02                         = Form("TrueClusBG_M02_in_Pt_Bin%02d", iPt);
        TString fNamePartConvHistoM02                   = Form("TrueClusPartConvMerged_M02_in_Pt_Bin%02d", iPt);
        TString fNamePartConvPi0HistoM02                = Form("TrueClusPartConvPi0_M02_in_Pt_Bin%02d", iPt);
        TString fNamePartConvEtaHistoM02                = Form("TrueClusPartConvEta_M02_in_Pt_Bin%02d", iPt);

        // true merged clusters        
        CheckForNULLForPointer(fHistoTrueClusMergedM02PtBin[iPt]);
        if (fTrueMergedPtVsM02Dummy){
            fHistoTrueClusMergedM02PtBin[iPt]                   = FillProjectionY(fTrueMergedPtVsM02Dummy, fNameHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
        // true merged clusters part conv       
        CheckForNULLForPointer(fHistoTrueClusPartConvMergedM02PtBin[iPt]);
        if (fTruePartConvMergedVsM02Dummy){
            fHistoTrueClusPartConvMergedM02PtBin[iPt]           = FillProjectionY(fTruePartConvMergedVsM02Dummy, fNamePartConvHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
        // adding merged cluster histograms
        CheckForNULLForPointer(fHistoTrueClusFullMergedM02PtBin[iPt]);
        fHistoTrueClusFullMergedM02PtBin[iPt]                   = AddPossiblyNotExistingHists( fHistoTrueClusMergedM02PtBin[iPt], fHistoTrueClusPartConvMergedM02PtBin[iPt],
                                                                                               Form("TrueClusMergedFull_M02_in_Pt_Bin%02d", iPt));
        
        // true merged pi0        
        CheckForNULLForPointer(fHistoTrueClusPi0M02PtBin[iPt]);
        if (fTruePi0PtVsM02){
            fHistoTrueClusPi0M02PtBin[iPt]                      = FillProjectionY(fTruePi0PtVsM02, fNamePi0HistoM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
        // true merged pi0 part conv
        CheckForNULLForPointer(fHistoTrueClusPartConvPi0M02PtBin[iPt]);
        if (fTruePartConvPi0PtVsM02){
            fHistoTrueClusPartConvPi0M02PtBin[iPt]              = FillProjectionY(fTruePartConvPi0PtVsM02, fNamePartConvPi0HistoM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
        // adding pi0 cluster histograms
        CheckForNULLForPointer(fHistoTrueClusFullPi0M02PtBin[iPt]);
        fHistoTrueClusFullPi0M02PtBin[iPt]                      = AddPossiblyNotExistingHists( fHistoTrueClusPi0M02PtBin[iPt], fHistoTrueClusPartConvPi0M02PtBin[iPt],
                                                                                               Form("TrueClusFullPi0_M02_in_Pt_Bin%02d", iPt));

        // true merged eta        
        CheckForNULLForPointer(fHistoTrueClusEtaM02PtBin[iPt]);
        if (fTrueEtaPtVsM02){
            fHistoTrueClusEtaM02PtBin[iPt]                      = FillProjectionY(fTrueEtaPtVsM02, fNameEtaHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
        // true merged eta part conv
        CheckForNULLForPointer(fHistoTrueClusPartConvEtaM02PtBin[iPt]);
        if (fTruePartConvEtaPtVsM02){
            fHistoTrueClusPartConvEtaM02PtBin[iPt]              = FillProjectionY(fTruePartConvEtaPtVsM02, fNamePartConvEtaHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
        // adding eta cluster histograms
        CheckForNULLForPointer(fHistoTrueClusFullEtaM02PtBin[iPt]);
        fHistoTrueClusFullEtaM02PtBin[iPt]                      = AddPossiblyNotExistingHists( fHistoTrueClusEtaM02PtBin[iPt], fHistoTrueClusPartConvEtaM02PtBin[iPt],
                                                                                               Form("TrueClusFullEta_M02_in_Pt_Bin%02d", iPt));

        // true gamma
        CheckForNULLForPointer(fHistoTrueClusGammaM02PtBin[iPt]);
        if (fTrueGammaPtVsM02){
            fHistoTrueClusGammaM02PtBin[iPt]                    = FillProjectionY(fTrueGammaPtVsM02, fNameGammaHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }

        // true electron
        CheckForNULLForPointer(fHistoTrueClusElectronM02PtBin[iPt]);
        if (fTrueElectronPtVsM02){
            fHistoTrueClusElectronM02PtBin[iPt]                 = FillProjectionY(fTrueElectronPtVsM02, fNameElectronHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }

        // true BG
        CheckForNULLForPointer(fHistoTrueClusBGM02PtBin[iPt]);
        if (fTrueBGPtVsM02){
            fHistoTrueClusBGM02PtBin[iPt]                       = FillProjectionY(fTrueBGPtVsM02, fNameBGHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
    }
}

//****************************************************************************
//******** Fill array of M02 histograms in pT slices for Pi0 Prim/Sec ********
//****************************************************************************
void FillMCPrimSecM02HistosArray(   TH2F* fTruePi0PrimM02VsPt, 
                                        TH2F* fTruePi0SecM02VsPt, 
                                        TH2F* fTruePi0SecFromK0sM02VsPt,
                                        TH2F* fTruePi0SecFromLambdaM02VsPt,
                                        TH2F* fTruePi0PrimPartConvM02VsPt, 
                                        TH2F* fTruePi0SecPartConvM02VsPt, 
                                        TH2F* fTruePi0SecFromK0sPartConvM02VsPt,
                                        TH2F* fTruePi0SecFromLambdaPartConvM02VsPt
                              ) {

    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        TString fNamePi0PrimM02                 = Form("TrueClusPrimPi0_M02_in_Pt_Bin%02d", iPt);
        TString fNamePi0SecM02                  = Form("TrueClusSecPi0_M02_in_Pt_Bin%02d", iPt);
        TString fNameSecPi0FK0sM02              = Form("TrueClusSecPi0FromK0s_M02_in_Pt_Bin%02d", iPt);
        TString fNameSecPi0FLambdaM02           = Form("TrueClusSecPi0FromLambda_M02_in_Pt_Bin%02d", iPt);
        TString fNamePi0PrimPartConvM02         = Form("TrueClusPartConvPrimPi0_M02_in_Pt_Bin%02d", iPt);
        TString fNamePi0SecPartConvM02          = Form("TrueClusPartConvSecPi0_M02_in_Pt_Bin%02d", iPt);
        TString fNameSecPi0FK0sPartConvM02      = Form("TrueClusPartConvSecPi0FromK0s_M02_in_Pt_Bin%02d", iPt);
        TString fNameSecPi0FLambdaPartConvM02   = Form("TrueClusPartConvSecPi0FromLambda_M02_in_Pt_Bin%02d", iPt);
        
        // true prim pi0 clusters
        CheckForNULLForPointer(fHistoTrueClusPrimPi0M02PtBin[iPt]);
        if (fTruePi0PrimM02VsPt){
            fHistoTrueClusPrimPi0M02PtBin[iPt]                  = FillProjectionY(fTruePi0PrimM02VsPt, fNamePi0PrimM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
        // true pi0 prim clusters part conv
        CheckForNULLForPointer(fHistoTrueClusPartConvPrimPi0M02PtBin[iPt]);
        if (fTruePi0PrimPartConvM02VsPt){
            fHistoTrueClusPartConvPrimPi0M02PtBin[iPt]          = FillProjectionY(fTruePi0PrimPartConvM02VsPt, fNamePi0PrimPartConvM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
        // adding prim pi0 cluster histograms
        CheckForNULLForPointer(fHistoTrueClusFullPrimPi0M02PtBin[iPt]);
        fHistoTrueClusFullPrimPi0M02PtBin[iPt]                  = AddPossiblyNotExistingHists( fHistoTrueClusPrimPi0M02PtBin[iPt], fHistoTrueClusPartConvPrimPi0M02PtBin[iPt],
                                                                                               Form("TrueClusFullPi0Prim_M02_in_Pt_Bin%02d", iPt));

        // true pi0 clusters from sec
        CheckForNULLForPointer(fHistoTrueClusSecPi0M02PtBin[iPt]);
        if (fTruePi0SecM02VsPt){
            fHistoTrueClusSecPi0M02PtBin[iPt]                   = FillProjectionY(fTruePi0SecM02VsPt, fNamePi0SecM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
        // true pi0 from sec clusters part conv
        CheckForNULLForPointer(fHistoTrueClusPartConvSecPi0M02PtBin[iPt]);
        if (fTruePi0SecPartConvM02VsPt){
            fHistoTrueClusPartConvSecPi0M02PtBin[iPt]           = FillProjectionY(fTruePi0SecPartConvM02VsPt, fNamePi0SecPartConvM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
        // adding pi0 from sec cluster histograms
        CheckForNULLForPointer(fHistoTrueClusFullSecPi0M02PtBin[iPt]);
        fHistoTrueClusFullSecPi0M02PtBin[iPt]                  = AddPossiblyNotExistingHists( fHistoTrueClusSecPi0M02PtBin[iPt], fHistoTrueClusPartConvSecPi0M02PtBin[iPt],
                                                                                               Form("TrueClusFullPi0Sec_M02_in_Pt_Bin%02d", iPt));

        // true pi0 from sec from K0s clusters
        CheckForNULLForPointer(fHistoTrueClusSecPi0FK0sM02PtBin[iPt]);
        if (fTruePi0SecFromK0sM02VsPt){
            fHistoTrueClusSecPi0FK0sM02PtBin[iPt]               = FillProjectionY(fTruePi0SecFromK0sM02VsPt, fNameSecPi0FK0sM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
        // true pi0 from sec from K0s clusters part conv
        CheckForNULLForPointer(fHistoTrueClusPartConvSecPi0FK0sM02PtBin[iPt]);
        if (fTruePi0SecFromK0sPartConvM02VsPt){
            fHistoTrueClusPartConvSecPi0FK0sM02PtBin[iPt]       = FillProjectionY(fTruePi0SecFromK0sPartConvM02VsPt, fNameSecPi0FK0sPartConvM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
        // adding pi0 from sec from K0s cluster histograms
        CheckForNULLForPointer(fHistoTrueClusFullSecPi0FK0sM02PtBin[iPt]);
        fHistoTrueClusFullSecPi0FK0sM02PtBin[iPt]               = AddPossiblyNotExistingHists( fHistoTrueClusSecPi0FK0sM02PtBin[iPt], fHistoTrueClusPartConvSecPi0FK0sM02PtBin[iPt],
                                                                                               Form("TrueClusFullPi0SecFK0s_M02_in_Pt_Bin%02d", iPt));

        // true pi0 from sec from Lambda clusters
        CheckForNULLForPointer(fHistoTrueClusSecPi0FLambdaM02PtBin[iPt]);
        if (fTruePi0SecFromLambdaM02VsPt){
            fHistoTrueClusSecPi0FLambdaM02PtBin[iPt]            = FillProjectionY(fTruePi0SecFromLambdaM02VsPt, fNameSecPi0FLambdaM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
        // true pi0 from sec from Lambda clusters part conv
        CheckForNULLForPointer(fHistoTrueClusPartConvSecPi0FLambdaM02PtBin[iPt]);
        if (fTruePi0SecFromLambdaPartConvM02VsPt){
            fHistoTrueClusPartConvSecPi0FLambdaM02PtBin[iPt]    = FillProjectionY(fTruePi0SecFromLambdaPartConvM02VsPt, fNameSecPi0FLambdaPartConvM02, fBinsPt[iPt], fBinsPt[iPt+1], 4);
        }
        // adding pi0 from sec from Lambda cluster histograms
        CheckForNULLForPointer(fHistoTrueClusFullSecPi0FLambdaM02PtBin[iPt]);
        fHistoTrueClusFullSecPi0FLambdaM02PtBin[iPt]            = AddPossiblyNotExistingHists( fHistoTrueClusSecPi0FLambdaM02PtBin[iPt], fHistoTrueClusPartConvSecPi0FLambdaM02PtBin[iPt],
                                                                                               Form("TrueClusFullPi0SecFLambda_M02_in_Pt_Bin%02d", iPt));
    }
}

//****************************************************************************
//***************** Filling MC BG separated **********************************
//****************************************************************************
void FillMCBGSeparated (TH2F* dummy2D){
    if (dummy2D){
        for (Int_t i = 0; i< 9; i++){
            TH1D* dummy1D               = (TH1D*)dummy2D->ProjectionX(Form("TrueClusBG_%s_Pt",labelsBG[i].Data()),i+1,i+1,"e");
            fHistoTrueClustersBGPt[i]   = (TH1D*)dummy1D->Rebin(fNBinsPt,Form("TrueClusBG_%s_Pt",labelsBG[i].Data()),fBinsPt);
        }
    }
}    


//****************************************************************************
//***************** Filling of MC histograms in proper binning ***************
//****************************************************************************
void FillHistosArrayMC(TH1D* fHistoMCMesonWithinAccepPtFill, TH1D * fHistoMCMesonPtFill, TH1D * fDeltaPtFill) {
    fHistoMCMesonWithinAccepPtFill->Sumw2();
    fHistoMCMesonWithinAccepPtRebin = (TH1D*)fHistoMCMesonWithinAccepPtFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
    fHistoMCMesonWithinAccepPtRebin->Divide(fDeltaPtFill);
    fHistoMCMesonPtFill->Sumw2();
    fHistoMCMesonPtRebin = (TH1D*)fHistoMCMesonPtFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
    fHistoMCMesonPtRebin->Divide(fDeltaPtFill);
}


//****************************************************************************
//*** Integration of Invariant Mass Histogram in given integration window ****
//****************************************************************************
void IntegrateHistoInvMass(TH1D * fHistoDummy, Double_t * fMesonIntRangeInt){
    if (fHistoDummy){
        Int_t binLowMassMeson   = fHistoDummy->GetXaxis()->FindBin(fMesonIntRangeInt[0]);
        Int_t binHighMassMeson  = fHistoDummy->GetXaxis()->FindBin(fMesonIntRangeInt[1]);
        fYields                 = fHistoDummy->IntegralAndError(binLowMassMeson,binHighMassMeson,fYieldsError);
    } else {
        fYields                 = 0.;
        fYieldsError            = 0.;
    }
}

//****************************************************************************
//********* Integration of M02 Histogram in given integration window *********
//****************************************************************************
void IntegrateHistoM02(TH1D * fHistoDummy, Double_t * fMesonIntRangeInt){
    if (fHistoDummy){
        Int_t binLowM02Meson    = fHistoDummy->GetXaxis()->FindBin(fMesonIntRangeInt[0]);
        Int_t binHighM02Meson   = fHistoDummy->GetXaxis()->FindBin(fMesonIntRangeInt[1]);
        fYields                 = fHistoDummy->IntegralAndError(binLowM02Meson,binHighM02Meson,fYieldsError);
    } else {
        fYieldsError            = 0.;
        fYields                 = 0.;
    }    
}

//****************************************************************************
//***************** Calculation of Meson Acceptance **************************
//****************************************************************************
void CalculateMesonAcceptance() {
    fHistoMCAcceptancePt        = new TH1D("fHistoMCAcceptancePt","",fNBinsPt,fBinsPt);
    fHistoMCAcceptancePt->Sumw2();
    fHistoMCAcceptancePt->Divide(fHistoMCMesonWithinAccepPtRebin,fHistoMCMesonPtRebin,1.,1.,"B");
}

//****************************************************************************
//***************** Calculation of Meson Efficiency **************************
//****************************************************************************
TH1D* CalculateMesonEfficiency(TH1D* fMC_fMesonYieldsPt, TH1D* fMC_SecondaryYieldPt, TString nameEfi ) {
    TH1D* fHistoMCMesonEffiPt   = (TH1D*)fMC_fMesonYieldsPt->Clone(nameEfi.Data());
    fHistoMCMesonEffiPt->Sumw2();
    if (fMC_SecondaryYieldPt) fHistoMCMesonEffiPt->Add(fMC_SecondaryYieldPt,-1.);
    fHistoMCMesonEffiPt->Divide(fHistoMCMesonEffiPt,fHistoMCMesonWithinAccepPtRebin,1.,1.,"B");
    return fHistoMCMesonEffiPt;
}

//****************************************************************************
//***************** Calculation of Meson Purity **************************
//****************************************************************************
TH1D* CalculatePurity(TH1D* histoYieldsRec, TH1D* histoYieldVal, TString namePur ){
    TH1D* fHistoPur             = (TH1D*) histoYieldVal->Clone(namePur.Data());
    fHistoPur->Sumw2();
    fHistoPur->Divide(fHistoPur,histoYieldsRec,1.,1.,"B");
    return fHistoPur;
}

//****************************************************************************
//****** Definition Crystal ball function for signal +linear background  *****
//*******parameters are:                                                   *****
//*******               - 0 normalization                                  *****
//*******               - 1 mean                                           *****
//*******               - 2 sigma                                           *****
//*******               - 3 n                                               *****
//*******               - 4 alpha                                           *****
//****************************************************************************
Double_t CrystalBall(Double_t *x,Double_t *par) {
    // The Crystal Ball shape is a Gaussian that is 'connected' to an exponential tail at
    // 'alpha' sigma of the Gau   ssian. The sign determines if it happens on the left or
    // right side. The 'n' parameter controls the slope of the exponential part.
    // typical par limits:
    //    1.0 < alpha < 5.0
    // and 0.5 < n < 100.0
    Double_t alpha      = par[4];
    Double_t n          = par[3];
    Double_t meanx      = par[1];
    Double_t sigma      = par[2];
    Double_t nn         = par[0];
    Double_t a          = TMath::Power((n/TMath::Abs(alpha)), n) * TMath::Exp(-0.5*alpha*alpha);
    Double_t b          = n/TMath::Abs(alpha) - TMath::Abs(alpha);
    Double_t arg        = (x[0] - meanx)/sigma;
    Double_t fitval     = 0;
    if (arg > -1.0*alpha) {
        fitval          = nn * TMath::Exp(-0.5*arg*arg);
    } else {
        fitval          = nn * a * TMath::Power((b-arg), (-1*n));
    }
    return fitval;
}

//****************************************************************************
//****** Definition Crystal ball function for signal +linear background  *****
//*******parameters are:                                                   *****
//*******               - 0 normalization                                  *****
//*******               - 1 mean                                         *****
//*******               - 2 sigma                                         *****
//*******               - 3 n                                            *****
//*******               - 4 alpha                                         *****
//*******               - 5 constant BG                                   *****
//*******               - 6 linear BG                                      *****
//****************************************************************************
Double_t CrystalBallBck(Double_t *x,Double_t *par) {
    return CrystalBall(x,par) + LinearBackground(x,&par[5]);
}


//****************************************************************************
//******************** Definition of linear BG function **********************
//*******parameters are:                                               *****
//*******               - 0 constant BG                                *****
//*******               - 1 linear BG                                  *****
//****************************************************************************
Double_t LinearBackground(Double_t *x,Double_t *par) {
    return par[0] + par[1]*x[0];
}

//****************************************************************************
//******** Calculation of FWHM for Gaussian + left side exponential  *********
//****************************************************************************
void CalculateFWHM(TF1 * fFunc){
// Default function
    if (fCrysFitting == 0){
        TF1* fFunc_def;
        fFunc_def               = new TF1("fFunc_def","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))", 
                                          fMesonFitRange[0], fMesonFitRange[1]);
        fFunc_def->SetParameter(0,fFunc->GetParameter(0));
        fFunc_def->SetParameter(1,fFunc->GetParameter(1));
        fFunc_def->SetParameter(2,fFunc->GetParameter(2));
        fFunc_def->SetParameter(3,fFunc->GetParameter(3));
    
        //FWHM
        fFWHMFunc               = fFunc_def->GetX( fFunc_def->GetParameter(0)*0.5, fFunc_def->GetParameter(1), fMesonFitRange[1]) - fFunc_def->GetX( fFunc_def->GetParameter(0)*0.5, 
                                  fMesonFitRange[0],fFunc_def->GetParameter(1));

        //FWHM error +
        TF1* fFunc_plus;
        fFunc_plus              = new TF1("fFunc_plus","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",
                                          fMesonFitRange[0], fMesonFitRange[1]);
        fFunc_plus->SetParameter(0,fFunc->GetParameter(0) + fFunc->GetParError(0));
        fFunc_plus->SetParameter(1,fFunc->GetParameter(1) + fFunc->GetParError(1));
        fFunc_plus->SetParameter(2,fFunc->GetParameter(2) + fFunc->GetParError(2));
        fFunc_plus->SetParameter(3,fFunc->GetParameter(3) + fFunc->GetParError(3));
        Double_t FWHM_plus      = fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,fFunc_plus->GetParameter(1), fMesonFitRange[1]) - fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5, 
                                  fMesonFitRange[0],fFunc_plus->GetParameter(1));

        //FWHM error -
        TF1* fFunc_minus;
        fFunc_minus             = new TF1("fFunc_minus","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))", 
                                          fMesonFitRange[0], fMesonFitRange[1]);
        fFunc_minus->SetParameter(0,fFunc->GetParameter(0) - fFunc->GetParError(0));
        fFunc_minus->SetParameter(1,fFunc->GetParameter(1) - fFunc->GetParError(1));
        fFunc_minus->SetParameter(2,fFunc->GetParameter(2) - fFunc->GetParError(2));
        fFunc_minus->SetParameter(3,fFunc->GetParameter(3) - fFunc->GetParError(3));
        
        Double_t FWHM_minus     = fFunc_minus->GetX( fFunc_minus->GetParameter(0)*0.5, fFunc_minus->GetParameter(1), fMesonFitRange[1] ) - fFunc_minus->GetX( fFunc_minus->GetParameter(0)*0.5, 
                                  fMesonFitRange[0],fFunc_minus->GetParameter(1));
        Double_t Error1         = TMath::Abs(fFWHMFunc-FWHM_plus);
        Double_t Error2         = TMath::Abs(fFWHMFunc-FWHM_minus);
        if(Error1>=Error2) 
            fFWHMFuncError      = Error1;
        if(Error1<Error2) 
            fFWHMFuncError      = Error2;
    } else {
        fFWHMFunc               = fFunc->GetParameter(2)*2.35;
        fFWHMFuncError          = fFunc->GetParError(2)*2.35;
    }
}

//****************************************************************************
//********************** Calculate Fraction of Secondaries *******************
//****************************************************************************
TH1D* CalculateSecondaryFractions(TH1D* histoRawYield, TH1D* histoRawYieldSec, TString nameHistoFrac){
    histoRawYield->Sumw2();
    histoRawYieldSec->Sumw2();
    TH1D* histoFracSec = (TH1D*)histoRawYieldSec->Clone(nameHistoFrac.Data());
    histoFracSec->Divide(histoFracSec,histoRawYield,1.,1.,"B");
    return histoFracSec;
}


//****************************************************************************
//****** Saving of general histograms, fits needed in the analysis ***********
//****** RAW output file, no correction histograms ***************************
//****************************************************************************
void SaveHistos(Int_t optionMC, TString fCutID, TString fPrefix3) {

    Int_t fNBinsClusterPt           = 64;
    Double_t fBinsClusterPt[65]     = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                       1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                       2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8,
                                       4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8,
                                       6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.4, 7.8, 8.2, 8.6,
                                       9.0, 9.5, 10,  11,  12,  14,  16,  18,  20,  25,
                                       30,  35,  40,  45,  50 
                                      };

    TH1D*   fDeltaPtCluster         = new TH1D("fDeltaPtCluster", "", fNBinsClusterPt, fBinsClusterPt);
    for(Int_t iPt=1;iPt<fNBinsClusterPt+1;iPt++){
        fDeltaPtCluster->SetBinContent(iPt,fBinsClusterPt[iPt]-fBinsClusterPt[iPt-1]);
        fDeltaPtCluster->SetBinError(iPt,0);
    }
    
    TH1D*   fHistoClustersPtPerEvent   = NULL;
    if (fHistoClustersPt){
        fHistoClustersPtPerEvent    = (TH1D*)fHistoClustersPt->Rebin(fNBinsClusterPt,"ClusterPtPerEvent",fBinsClusterPt);
        fHistoClustersPtPerEvent->Divide(fDeltaPtCluster);
        fHistoClustersPtPerEvent->Scale(1./fNEvents);
    }
    TH1D*   fHistoClustersOverlapHeadersPtPerEvent   = NULL;
    if (fHistoClustersOverlapHeadersPt){
        fHistoClustersOverlapHeadersPtPerEvent   = (TH1D*)fHistoClustersOverlapHeadersPt->Rebin(fNBinsClusterPt,"ClusterOverlapHeadersPtPerEvent",fBinsClusterPt);
        fHistoClustersOverlapHeadersPtPerEvent->Divide(fDeltaPtCluster);
        fHistoClustersOverlapHeadersPtPerEvent->Scale(1./fNEvents);
    }
    const char* nameOutput          = Form("%s/%s/%s_%s_GammaMergedWithoutCorrection%s%s_%s.root", fCutSelection.Data(), fEnergyFlag.Data(), fPrefix.Data(), fPrefix3.Data(), 
                                           fPeriodFlag.Data(), fAdditionalName.Data(), fCutID.Data());
    fOutput1                        = new TFile(nameOutput, "RECREATE");

        cout << "Begin writing Uncorrected File" << endl;
        if (fEventQuality)                              fEventQuality->Write("NEvents");

        if (fHistoClustersPt)                           fHistoClustersPt->Write("ClusterPt");
        if (fHistoClustersPtPerEvent)                   fHistoClustersPtPerEvent->Write("ClusterPtPerEvent");
        if (fDeltaPtCluster)                            fDeltaPtCluster->Write();
        
        if (fHistoClustersOverlapHeadersPt)             fHistoClustersOverlapHeadersPt->Write("ClusterOverlapHeadersPt");
        if (fHistoClustersOverlapHeadersPtPerEvent)     fHistoClustersOverlapHeadersPtPerEvent->Write("ClusterOverlapHeadersPtPerEvent");
        
        if(fHistoClustersMergedPtM02)                   fHistoClustersMergedPtM02->Write("ClusterMergedPtM02");
        if(fHistoClustersMergedPtM02AccMeson)           fHistoClustersMergedPtM02AccMeson->Write("ClusterMergedPtM02AccMeson");
        if(fHistoTrackVsClusterCandidates)              fHistoTrackVsClusterCandidates->Write("TrackCandidatesVsClusters");

        if (fAdvancedClusterQA){
            if (fHistoClusterNCellsPt)                  fHistoClusterNCellsPt->Write("ClusterNCellsPt");
            if (fHistoClusterMergedNCellsPt)            fHistoClusterMergedNCellsPt->Write("ClusterMergedNCellsPt");        
            if (fHistoClusterMergedNCellsArClPt)        fHistoClusterMergedNCellsArClPt->Write("ClusterMergedNCellsAroundClusterPt");        
            if (fHistoClusterMergedNCellsArAInclPt)     fHistoClusterMergedNCellsArAInclPt->Write("ClusterMergedNCellsAroundAndInClusterPt");        
            if (fHistoClusterMergedEAroundClE)          fHistoClusterMergedEAroundClE->Write("ClusterMergedEAroundClusterE");        
            if (fHistoClusterMergedNPartPt)             fHistoClusterMergedNPartPt->Write("ClusterMergedNPartPointingToClusPt");        
        }    
        
        for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
            if (fHistoInvMassPtBin[ii])                     fHistoInvMassPtBin[ii]->Write(Form("InvMass_PtBin%i",ii));
            if (fHistoM02PtBin[ii])                         fHistoM02PtBin[ii]->Write(Form("M02_PtBin%i",ii));
            if (fIsMC){
                if (fHistoTrueClusFullMergedInvMassPtBin[ii])   fHistoTrueClusFullMergedInvMassPtBin[ii]->Write(Form("ValidatedMerged_InvMass_PtBin%i",ii));
                if (fHistoTrueClusFullMergedM02PtBin[ii])       fHistoTrueClusFullMergedM02PtBin[ii]->Write(Form("ValidatedMerged_M02_PtBin%i",ii));
            }
        }
        cout << "End writing Uncorrected File" << endl;
    
        if (fDeltaPt)                                   fDeltaPt->Write();
        if (fHistoYieldMesonM02)                        fHistoYieldMesonM02->Write();
        
        
    fOutput1->Write();
    fOutput1->Close();
}


//****************************************************************************
//****** Saving of MC histograms needed for corrections  *********************
//****** or comparisons to data at a later stage *****************************
//****************************************************************************
void SaveCorrectionHistos(TString fCutID, TString fPrefix3){
    const char* nameOutput = Form("%s/%s/%s_%s_GammaMergedCorrectionHistos%s%s_%s.root", fCutSelection.Data(), fEnergyFlag.Data(), fPrefix.Data(), fPrefix3.Data(),
                                  fPeriodFlag.Data(), fAdditionalName.Data(), fCutID.Data());
    fOutput2 = new TFile(nameOutput,"RECREATE");
    cout << "======================================================" << endl;
    cout << nameOutput << endl;
    cout << "======================================================" << endl;

    if (fIsMC){
        if (fEventQuality)                          fEventQuality->Write("NEvents");
        if (fHistoTrueYieldMergedM02)               fHistoTrueYieldMergedM02->Write();
        if (fHistoTrueYieldMergedPartConvM02)       fHistoTrueYieldMergedPartConvM02->Write();
        if (fHistoTrueYieldMergedFullM02)           fHistoTrueYieldMergedFullM02->Write();
        if (fHistoTrueYieldPi0M02)                  fHistoTrueYieldPi0M02->Write();
        if (fHistoTrueYieldPi0PartConvM02)          fHistoTrueYieldPi0PartConvM02->Write();
        if (fHistoTrueYieldPi0FullM02)              fHistoTrueYieldPi0FullM02->Write();
        if (fHistoTrueYieldEtaM02)                  fHistoTrueYieldEtaM02->Write();
        if (fHistoTrueYieldEtaPartConvM02)          fHistoTrueYieldEtaPartConvM02->Write();
        if (fHistoTrueYieldEtaFullM02)              fHistoTrueYieldEtaFullM02->Write();
        if (fHistoTrueYieldGammaM02)                fHistoTrueYieldGammaM02->Write();
        if (fHistoTrueYieldElectronM02)             fHistoTrueYieldElectronM02->Write();
        if (fHistoTrueYieldBGM02)                   fHistoTrueYieldBGM02->Write();
        if (fHistoTrueYieldPrimPi0M02)              fHistoTrueYieldPrimPi0M02->Write();
        if (fHistoTrueYieldPrimPi0PartConvM02)      fHistoTrueYieldPrimPi0PartConvM02->Write();
        if (fHistoTrueYieldPrimPi0FullM02)          fHistoTrueYieldPrimPi0FullM02->Write();
        if (fHistoTrueYieldSecPi0M02)               fHistoTrueYieldSecPi0M02->Write();
        if (fHistoTrueYieldSecPi0PartConvM02)       fHistoTrueYieldSecPi0PartConvM02->Write();
        if (fHistoTrueYieldSecPi0FullM02)           fHistoTrueYieldSecPi0FullM02->Write();
        if (fHistoTrueYieldSecPi0FK0sM02)           fHistoTrueYieldSecPi0FK0sM02->Write();
        if (fHistoTrueYieldSecPi0FK0sPartConvM02)   fHistoTrueYieldSecPi0FK0sPartConvM02->Write();
        if (fHistoTrueYieldSecPi0FK0sFullM02)       fHistoTrueYieldSecPi0FK0sFullM02->Write();
        if (fHistoTrueYieldSecPi0FLambdaM02)        fHistoTrueYieldSecPi0FLambdaM02->Write();
        if (fHistoTrueYieldSecPi0FLambdaPartConvM02)fHistoTrueYieldSecPi0FLambdaPartConvM02->Write();
        if (fHistoTrueYieldSecPi0FLambdaFullM02)    fHistoTrueYieldSecPi0FLambdaFullM02->Write();
        if (fHistoTrueClustersBGPtSource)           fHistoTrueClustersBGPtSource->Write();

        if (fHistoMCMesonPt)                        fHistoMCMesonPt->Write("MC_Meson");
        if (fHistoMCMesonPtRebin)                   fHistoMCMesonPtRebin->Write("MC_Meson_Rebin");
        if (fHistoMCMesonWithinAccepPt)             fHistoMCMesonWithinAccepPt->Write("MC_Meson_WithinAccept");
        if (fHistoMCMesonWithinAccepPtRebin)        fHistoMCMesonWithinAccepPtRebin->Write("MC_Meson_WithinAccept_Rebin");
        if (fHistoMCMesonPtWOWeights)               fHistoMCMesonPtWOWeights->Write("MC_Meson_WOWeights");
        if (fHistoMCAcceptancePt)                   fHistoMCAcceptancePt->Write("AcceptancePt");
        if (fHistoTrueEffiMerged)                   fHistoTrueEffiMerged->Write("TrueEfficiencyMergedPt");
        if (fHistoTrueEffiPrimMeson)                fHistoTrueEffiPrimMeson->Write("TrueEfficiencyPrimMesonPt");
        if (fHistoTruePurityMerged)                 fHistoTruePurityMerged->Write("TruePurityMergedPt");
        if (fHistoTruePi0PurityMerged)              fHistoTruePi0PurityMerged->Write("TruePurityPi0Pt");
        if (fHistoTrueEtaPurityMerged)              fHistoTrueEtaPurityMerged->Write("TruePurityEtaPt");
        if (fHistoTruePi0SecFrac)                   fHistoTruePi0SecFrac->Write();
        if (fHistoTruePi0SecFracFK0S)               fHistoTruePi0SecFracFK0S->Write();
        if (fHistoTruePi0SecFracFLambda)            fHistoTruePi0SecFracFLambda->Write();
        for (Int_t i = 0; i< 9; i++){
            if (fHistoTrueClustersBGPt[i]) fHistoTrueClustersBGPt[i]->Write();
        }
    }
    
    fOutput2->Write();
    fOutput2->Close();
}


//****************************************************************************
//************************* Deleting pointers ********************************
//****************************************************************************
void Delete(){
    if (fMesonMassPlotRange)                                    delete fMesonMassPlotRange;
    if (fMesonM02PlotRange)                                     delete fMesonM02PlotRange;
    
    if (fHistoClustersPt)                                       delete fHistoClustersPt;
    if (fHistoClustersOverlapHeadersPt)                         delete fHistoClustersOverlapHeadersPt;
    if (fHistoClustersMergedPtM02)                              delete fHistoClustersMergedPtM02;
    if (fHistoClustersMergedPtM02AccMeson)                      delete fHistoClustersMergedPtM02AccMeson;
    if (fHistoTrackVsClusterCandidates)                         delete fHistoTrackVsClusterCandidates;
    if (fHistoInvMassVsPt)                                      delete fHistoInvMassVsPt;
    if (fHistoSPDtrackletvsSPDclusters)                         delete fHistoSPDtrackletvsSPDclusters;
    if (fMesonM02Yields)                                        delete fMesonM02Yields;
    if (fMesonM02YieldsError)                                   delete fMesonM02YieldsError;

    if (fIsMC){
        if (fHistoTrueClustersMergedInvMassPt)                  delete fHistoTrueClustersMergedInvMassPt;
        if (fHistoTrueClustersMergedPtM02)                      delete fHistoTrueClustersMergedPtM02;
        if (fHistoTrueClustersPi0InvMassPt)                     delete fHistoTrueClustersPi0InvMassPt;
        if (fHistoTrueClustersPi0PtM02)                         delete fHistoTrueClustersPi0PtM02;
        if (fHistoTrueClustersEtaInvMassPt)                     delete fHistoTrueClustersEtaInvMassPt;
        if (fHistoTrueClustersEtaPtM02)                         delete fHistoTrueClustersEtaPtM02;
        if (fHistoTrueClustersGammaInvMassPt)                   delete fHistoTrueClustersGammaInvMassPt;
        if (fHistoTrueClustersGammaPtM02)                       delete fHistoTrueClustersGammaPtM02;
        if (fHistoTrueClustersElectronInvMassPt)                delete fHistoTrueClustersElectronInvMassPt;
        if (fHistoTrueClustersElectronPtM02)                    delete fHistoTrueClustersElectronPtM02;
        if (fHistoTrueClustersBGInvMassPt)                      delete fHistoTrueClustersBGInvMassPt;
        if (fHistoTrueClustersBGPtM02)                          delete fHistoTrueClustersBGPtM02;
        if (fHistoTrueClusPartConvMergedInvMassPt)              delete fHistoTrueClusPartConvMergedInvMassPt;
        if (fHistoTrueClusPartConvMergedPtM02)                  delete fHistoTrueClusPartConvMergedPtM02;
        if (fHistoTrueClusPartConvPi0InvMassPt)                 delete fHistoTrueClusPartConvPi0InvMassPt;
        if (fHistoTrueClusPartConvPi0PtM02)                     delete fHistoTrueClusPartConvPi0PtM02;
        if (fHistoTrueClusPartConvEtaInvMassPt)                 delete fHistoTrueClusPartConvEtaInvMassPt;
        if (fHistoTrueClusPartConvEtaPtM02)                     delete fHistoTrueClusPartConvEtaPtM02;
        if (fHistoTrueClustersBGPtSource)                       delete fHistoTrueClustersBGPtSource;
        if (fMesonM02TrueMergedYields)                          delete fMesonM02TrueMergedYields;
        if (fMesonM02TruePi0Yields)                             delete fMesonM02TruePi0Yields;
        if (fMesonM02TrueEtaYields)                             delete fMesonM02TrueEtaYields;
        if (fMesonM02TrueGammaYields)                           delete fMesonM02TrueGammaYields;
        if (fMesonM02TrueElectronYields)                        delete fMesonM02TrueElectronYields;
        if (fMesonM02TrueBGYields)                              delete fMesonM02TrueBGYields;
        if (fMesonM02TrueMergedYieldsError)                     delete fMesonM02TrueMergedYieldsError;
        if (fMesonM02TruePi0YieldsError)                        delete fMesonM02TruePi0YieldsError;
        if (fMesonM02TrueEtaYieldsError)                        delete fMesonM02TrueEtaYieldsError;
        if (fMesonM02TrueGammaYieldsError)                      delete fMesonM02TrueGammaYieldsError;
        if (fMesonM02TrueElectronYieldsError)                   delete fMesonM02TrueElectronYieldsError;
        if (fMesonM02TrueBGYieldsError)                         delete fMesonM02TrueBGYieldsError;
        if (fMesonM02TrueMergedPartConvYields)                  delete fMesonM02TrueMergedPartConvYields;
        if (fMesonM02TruePi0PartConvYields)                     delete fMesonM02TruePi0PartConvYields;
        if (fMesonM02TrueEtaPartConvYields)                     delete fMesonM02TrueEtaPartConvYields;
        if (fMesonM02TrueMergedPartConvYieldsError)             delete fMesonM02TrueMergedPartConvYieldsError;
        if (fMesonM02TruePi0PartConvYieldsError)                delete fMesonM02TruePi0PartConvYieldsError;
        if (fMesonM02TrueEtaPartConvYieldsError)                delete fMesonM02TrueEtaPartConvYieldsError;
        if (fMesonM02TrueMergedFullYields)                      delete fMesonM02TrueMergedFullYields;
        if (fMesonM02TruePi0FullYields)                         delete fMesonM02TruePi0FullYields;
        if (fMesonM02TrueEtaFullYields)                         delete fMesonM02TrueEtaFullYields;
        if (fMesonM02TrueMergedFullYieldsError)                 delete fMesonM02TrueMergedFullYieldsError;
        if (fMesonM02TruePi0FullYieldsError)                    delete fMesonM02TruePi0FullYieldsError;
        if (fMesonM02TrueEtaFullYieldsError)                    delete fMesonM02TrueEtaFullYieldsError;
        if (fMesonM02TruePrimPi0Yields)                         delete fMesonM02TruePrimPi0Yields;
        if (fMesonM02TrueSecPi0Yields)                          delete fMesonM02TrueSecPi0Yields;
        if (fMesonM02TrueSecPi0FK0sYields)                      delete fMesonM02TrueSecPi0FK0sYields;
        if (fMesonM02TrueSecPi0FLambdaYields)                   delete fMesonM02TrueSecPi0FLambdaYields;
        if (fMesonM02TruePrimPi0YieldsError)                    delete fMesonM02TruePrimPi0YieldsError;
        if (fMesonM02TrueSecPi0YieldsError)                     delete fMesonM02TrueSecPi0YieldsError;
        if (fMesonM02TrueSecPi0FK0sYieldsError)                 delete fMesonM02TrueSecPi0FK0sYieldsError;
        if (fMesonM02TrueSecPi0FLambdaYieldsError)              delete fMesonM02TrueSecPi0FLambdaYieldsError;
        if (fMesonM02TruePrimPi0PartConvYields)                 delete fMesonM02TruePrimPi0PartConvYields;
        if (fMesonM02TrueSecPi0PartConvYields)                  delete fMesonM02TrueSecPi0PartConvYields;
        if (fMesonM02TrueSecPi0FK0sPartConvYields)              delete fMesonM02TrueSecPi0FK0sPartConvYields;
        if (fMesonM02TrueSecPi0FLambdaPartConvYields)           delete fMesonM02TrueSecPi0FLambdaPartConvYields;
        if (fMesonM02TruePrimPi0PartConvYieldsError)            delete fMesonM02TruePrimPi0PartConvYieldsError;
        if (fMesonM02TrueSecPi0PartConvYieldsError)             delete fMesonM02TrueSecPi0PartConvYieldsError;
        if (fMesonM02TrueSecPi0FK0sPartConvYieldsError)         delete fMesonM02TrueSecPi0FK0sPartConvYieldsError;
        if (fMesonM02TrueSecPi0FLambdaPartConvYieldsError)      delete fMesonM02TrueSecPi0FLambdaPartConvYieldsError;
        if (fMesonM02TruePrimPi0FullYields)                     delete fMesonM02TruePrimPi0FullYields;
        if (fMesonM02TrueSecPi0FullYields)                      delete fMesonM02TrueSecPi0FullYields;
        if (fMesonM02TrueSecPi0FK0sFullYields)                  delete fMesonM02TrueSecPi0FK0sFullYields;
        if (fMesonM02TrueSecPi0FLambdaFullYields)               delete fMesonM02TrueSecPi0FLambdaFullYields;
        if (fMesonM02TruePrimPi0FullYieldsError)                delete fMesonM02TruePrimPi0FullYieldsError;
        if (fMesonM02TrueSecPi0FullYieldsError)                 delete fMesonM02TrueSecPi0FullYieldsError;
        if (fMesonM02TrueSecPi0FK0sFullYieldsError)             delete fMesonM02TrueSecPi0FK0sFullYieldsError;
        if (fMesonM02TrueSecPi0FLambdaFullYieldsError)          delete fMesonM02TrueSecPi0FLambdaFullYieldsError;
    }
    
    for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
        if (fHistoInvMassPtBin[ii])                             delete fHistoInvMassPtBin[ii];
        if (fHistoM02PtBin[ii])                                 delete fHistoM02PtBin[ii];
        if (fIsMC){
            if (fHistoTrueClusMergedInvMassPtBin[ii])               delete fHistoTrueClusMergedInvMassPtBin[ii];
            if (fHistoTrueClusMergedM02PtBin[ii])                   delete fHistoTrueClusMergedM02PtBin[ii];
            if (fHistoTrueClusPi0InvMassPtBin[ii])                  delete fHistoTrueClusPi0InvMassPtBin[ii];
            if (fHistoTrueClusPi0M02PtBin[ii])                      delete fHistoTrueClusPi0M02PtBin[ii];
            if (fHistoTrueClusEtaInvMassPtBin[ii])                  delete fHistoTrueClusEtaInvMassPtBin[ii];
            if (fHistoTrueClusEtaM02PtBin[ii])                      delete fHistoTrueClusEtaM02PtBin[ii];
            if (fHistoTrueClusGammaInvMassPtBin[ii])                delete fHistoTrueClusGammaInvMassPtBin[ii];
            if (fHistoTrueClusGammaM02PtBin[ii])                    delete fHistoTrueClusGammaM02PtBin[ii];
            if (fHistoTrueClusElectronInvMassPtBin[ii])             delete fHistoTrueClusElectronInvMassPtBin[ii];
            if (fHistoTrueClusElectronM02PtBin[ii])                 delete fHistoTrueClusElectronM02PtBin[ii];
            if (fHistoTrueClusBGInvMassPtBin[ii])                   delete fHistoTrueClusBGInvMassPtBin[ii];
            if (fHistoTrueClusBGM02PtBin[ii])                       delete fHistoTrueClusBGM02PtBin[ii];
            if (fHistoTrueClusPartConvMergedInvMassPtBin[ii])       delete fHistoTrueClusPartConvMergedInvMassPtBin[ii];
            if (fHistoTrueClusPartConvMergedM02PtBin[ii])           delete fHistoTrueClusPartConvMergedM02PtBin[ii];
            if (fHistoTrueClusPartConvPi0InvMassPtBin[ii])          delete fHistoTrueClusPartConvPi0InvMassPtBin[ii];
            if (fHistoTrueClusPartConvPi0M02PtBin[ii])              delete fHistoTrueClusPartConvPi0M02PtBin[ii];
            if (fHistoTrueClusPartConvEtaInvMassPtBin[ii])          delete fHistoTrueClusPartConvEtaInvMassPtBin[ii];
            if (fHistoTrueClusPartConvEtaM02PtBin[ii])              delete fHistoTrueClusPartConvEtaM02PtBin[ii];
            if (fHistoTrueClusFullMergedInvMassPtBin[ii])           delete fHistoTrueClusFullMergedInvMassPtBin[ii];
            if (fHistoTrueClusFullMergedM02PtBin[ii])               delete fHistoTrueClusFullMergedM02PtBin[ii];
            if (fHistoTrueClusFullPi0InvMassPtBin[ii])              delete fHistoTrueClusFullPi0InvMassPtBin[ii];
            if (fHistoTrueClusFullPi0M02PtBin[ii])                  delete fHistoTrueClusFullPi0M02PtBin[ii];
            if (fHistoTrueClusFullEtaInvMassPtBin[ii])              delete fHistoTrueClusFullEtaInvMassPtBin[ii];
            if (fHistoTrueClusFullEtaM02PtBin[ii])                  delete fHistoTrueClusFullEtaM02PtBin[ii];
            if (fHistoTrueClusPrimPi0InvMassPtBin[ii])              delete fHistoTrueClusPrimPi0InvMassPtBin[ii];
            if (fHistoTrueClusSecPi0InvMassPtBin[ii])               delete fHistoTrueClusSecPi0InvMassPtBin[ii];
            if (fHistoTrueClusSecPi0FK0sInvMassPtBin[ii])           delete fHistoTrueClusSecPi0FK0sInvMassPtBin[ii];
            if (fHistoTrueClusSecPi0FLambdaInvMassPtBin[ii])        delete fHistoTrueClusSecPi0FLambdaInvMassPtBin[ii];
            if (fHistoTrueClusPartConvPrimPi0InvMassPtBin[ii])      delete fHistoTrueClusPartConvPrimPi0InvMassPtBin[ii];
            if (fHistoTrueClusPartConvSecPi0InvMassPtBin[ii])       delete fHistoTrueClusPartConvSecPi0InvMassPtBin[ii];
            if (fHistoTrueClusPartConvSecPi0FK0sInvMassPtBin[ii])   delete fHistoTrueClusPartConvSecPi0FK0sInvMassPtBin[ii];
            if (fHistoTrueClusPartConvSecPi0FLambdaInvMassPtBin[ii])delete fHistoTrueClusPartConvSecPi0FLambdaInvMassPtBin[ii];
            if (fHistoTrueClusFullPrimPi0InvMassPtBin[ii])          delete fHistoTrueClusFullPrimPi0InvMassPtBin[ii];
            if (fHistoTrueClusFullSecPi0InvMassPtBin[ii])           delete fHistoTrueClusFullSecPi0InvMassPtBin[ii];
            if (fHistoTrueClusFullSecPi0FK0sInvMassPtBin[ii])       delete fHistoTrueClusFullSecPi0FK0sInvMassPtBin[ii];
            if (fHistoTrueClusFullSecPi0FLambdaInvMassPtBin[ii])    delete fHistoTrueClusFullSecPi0FLambdaInvMassPtBin[ii];
            if (fHistoTrueClusPrimPi0M02PtBin[ii])                  delete fHistoTrueClusPrimPi0M02PtBin[ii];
            if (fHistoTrueClusSecPi0M02PtBin[ii])                   delete fHistoTrueClusSecPi0M02PtBin[ii];
            if (fHistoTrueClusSecPi0FK0sM02PtBin[ii])               delete fHistoTrueClusSecPi0FK0sM02PtBin[ii];
            if (fHistoTrueClusSecPi0FLambdaM02PtBin[ii])            delete fHistoTrueClusSecPi0FLambdaM02PtBin[ii];
            if (fHistoTrueClusPartConvPrimPi0M02PtBin[ii])          delete fHistoTrueClusPartConvPrimPi0M02PtBin[ii];
            if (fHistoTrueClusPartConvSecPi0M02PtBin[ii])           delete fHistoTrueClusPartConvSecPi0M02PtBin[ii];
            if (fHistoTrueClusPartConvSecPi0FK0sM02PtBin[ii])       delete fHistoTrueClusPartConvSecPi0FK0sM02PtBin[ii];
            if (fHistoTrueClusPartConvSecPi0FLambdaM02PtBin[ii])    delete fHistoTrueClusPartConvSecPi0FLambdaM02PtBin[ii];
            if (fHistoTrueClusFullPrimPi0M02PtBin[ii])              delete fHistoTrueClusFullPrimPi0M02PtBin[ii];
            if (fHistoTrueClusFullSecPi0M02PtBin[ii])               delete fHistoTrueClusFullSecPi0M02PtBin[ii];
            if (fHistoTrueClusFullSecPi0FK0sM02PtBin[ii])           delete fHistoTrueClusFullSecPi0FK0sM02PtBin[ii];
            if (fHistoTrueClusFullSecPi0FLambdaM02PtBin[ii])        delete fHistoTrueClusFullSecPi0FLambdaM02PtBin[ii];
        }
    }
    
    if (fDeltaPt)                                               delete fDeltaPt;
    if (fHistoYieldMesonM02)                                    delete fHistoYieldMesonM02;
    if (fIsMC){
        if (fHistoTrueYieldMergedM02)                           delete fHistoTrueYieldMergedM02;
        if (fHistoTrueYieldMergedPartConvM02)                   delete fHistoTrueYieldMergedPartConvM02;
        if (fHistoTrueYieldMergedFullM02)                       delete fHistoTrueYieldMergedFullM02;
        if (fHistoTrueYieldPi0M02)                              delete fHistoTrueYieldPi0M02;
        if (fHistoTrueYieldPi0PartConvM02)                      delete fHistoTrueYieldPi0PartConvM02;
        if (fHistoTrueYieldPi0FullM02)                          delete fHistoTrueYieldPi0FullM02;
        if (fHistoTrueYieldEtaM02)                              delete fHistoTrueYieldEtaM02;
        if (fHistoTrueYieldEtaPartConvM02)                      delete fHistoTrueYieldEtaPartConvM02;
        if (fHistoTrueYieldEtaFullM02)                          delete fHistoTrueYieldEtaFullM02;
        if (fHistoTrueYieldGammaM02)                            delete fHistoTrueYieldGammaM02;
        if (fHistoTrueYieldElectronM02)                         delete fHistoTrueYieldElectronM02;
        if (fHistoTrueYieldBGM02)                               delete fHistoTrueYieldBGM02;
        if (fHistoTrueYieldPrimPi0M02)                          delete fHistoTrueYieldPrimPi0M02;
        if (fHistoTrueYieldPrimPi0PartConvM02)                  delete fHistoTrueYieldPrimPi0PartConvM02;
        if (fHistoTrueYieldPrimPi0FullM02)                      delete fHistoTrueYieldPrimPi0FullM02;
        if (fHistoTrueYieldSecPi0M02)                           delete fHistoTrueYieldSecPi0M02;
        if (fHistoTrueYieldSecPi0PartConvM02)                   delete fHistoTrueYieldSecPi0PartConvM02;
        if (fHistoTrueYieldSecPi0FullM02)                       delete fHistoTrueYieldSecPi0FullM02;
        if (fHistoTrueYieldSecPi0FK0sM02)                       delete fHistoTrueYieldSecPi0FK0sM02;
        if (fHistoTrueYieldSecPi0FK0sPartConvM02)               delete fHistoTrueYieldSecPi0FK0sPartConvM02;
        if (fHistoTrueYieldSecPi0FK0sFullM02)                   delete fHistoTrueYieldSecPi0FK0sFullM02;
        if (fHistoTrueYieldSecPi0FLambdaM02)                    delete fHistoTrueYieldSecPi0FLambdaM02;
        if (fHistoTrueYieldSecPi0FLambdaPartConvM02)            delete fHistoTrueYieldSecPi0FLambdaPartConvM02;
        if (fHistoTrueYieldSecPi0FLambdaFullM02)                delete fHistoTrueYieldSecPi0FLambdaFullM02;
        if (fHistoMCAcceptancePt)                               delete fHistoMCAcceptancePt;
        if (fHistoTruePi0SecFrac)                               delete fHistoTruePi0SecFrac;
        if (fHistoTruePi0SecFracFK0S)                           delete fHistoTruePi0SecFracFK0S;
        if (fHistoTruePi0SecFracFLambda)                        delete fHistoTruePi0SecFracFLambda;
        if (fHistoTruePurityMerged)                             delete fHistoTruePurityMerged;
        if (fHistoTruePi0PurityMerged)                          delete fHistoTruePi0PurityMerged;
        if (fHistoTrueEtaPurityMerged)                          delete fHistoTrueEtaPurityMerged;
        if (fHistoTrueEffiPrimMeson)                            delete fHistoTrueEffiPrimMeson;
        
    }
}