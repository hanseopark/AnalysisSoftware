// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

//This file is not supposed to be run on outputfiles of the GammaConv-Software before the 30th Sept 2010.

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
#include "../CommonHeaders/PlottingMeson.h"
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
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
//#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "ExtractSignalV2.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "../CommonHeaders/ExtractSignalPlotting.h"
#include "THnSparse.h"

                                                                  
//****************************************************************************
//************** Main function for extraction of signal **********************
//****************************************************************************
void ExtractSignalV2(   TString meson                   = "", 
                        TString file                    = "", 
                        TString cutSelection            = "", 
                        TString Suffix                  = "", 
                        TString optionMC                = "", 
                        TString optionEnergy            = "", 
                        TString optionCrystalBall       = "", 
                        TString directphotonPlots       = "", 
                        TString optionUseMinBiasEff     = "", 
                        TString optionPeriod            = "", 
                        TString optionAdvancedMesonQA   = "",
                        Int_t numberOfBins              = 30,
                        Bool_t addSig                   = kFALSE, 
                        Int_t mode                      = 9, 
                        Bool_t UseTHnSparse             = kTRUE,
                        Int_t triggerSet                = -1
                    ) {
    gROOT->Reset();

    fMode = mode;
    // mode:   0 // new output PCM-PCM
    //         1 // new output PCM dalitz
    //         2 // new output PCM-Calo
    //         3 // new output Calo-Calo
    //         4 // new output EMCAL-EMCAL
    //         5 // new output PHOS-PHOS
    //         9 // old output PCM-PCM
    
    if(directphotonPlots){}

    if(optionAdvancedMesonQA.Contains("AdvancedMesonQA")){fAdvancedMesonQA = kTRUE;}
    
    fCutSelection = cutSelection;
    TString fCutSelectionRead = cutSelection;
    if (mode == 9){
        ReturnSeparatedCutNumber(cutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection);
        fEventCutSelection = fGammaCutSelection(0,7);
        fGammaCutSelection = fGammaCutSelection(7,fGammaCutSelection.Length()-7);
        cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
    } else {
        ReturnSeparatedCutNumberAdvanced(cutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
    }
    //if(fMesonCutSelection.Length() < 16 ) fMesonCutSelection.Append("0");
    TString fEventCutSelectionRead = fEventCutSelection.Data();
    TString fGammaCutSelectionRead = fGammaCutSelection.Data();
    TString fMesonCutSelectionRead = fMesonCutSelection.Data();
    if (addSig) {
        cout << "running added Signal" << endl;
        cout << fEventCutSelection.Data() << endl;
        fEventCutSelection.Replace(GetEventRejectExtraSignalsCutPosition(),1,"2");
        cout << fEventCutSelection.Data() << endl;
        fEventCutSelectionRead = fEventCutSelection;
        fGammaCutSelectionRead = fGammaCutSelection;
        fMesonCutSelectionRead = fMesonCutSelection;
        if (mode==9)fCutSelectionRead = Form("%s%s_%s", fEventCutSelection.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
        else if (mode==0) fCutSelectionRead = Form("%s_%s_%s",fEventCutSelection.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
        if (mode==2 || mode==3)fCutSelectionRead = Form("%s_%s_%s_%s",fEventCutSelection.Data(), fGammaCutSelection.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
        cout << fCutSelectionRead.Data() << endl;
    }
    if(optionUseMinBiasEff.CompareTo("MinBiasEffOnly")==0 && optionMC.CompareTo("kTRUE") == 0){
        cout << "calculating MinBias Eff" << endl;
        cout << fEventCutSelection.Data() << endl;
        fEventCutSelection.Replace(GetEventCentralityMinCutPosition(),2,"00");
        fEventCutSelectionRead = fEventCutSelection;
        fGammaCutSelectionRead = fGammaCutSelection;
        fMesonCutSelectionRead = fMesonCutSelection;      
        cout << fGammaCutSelection.Data() << endl;
        if (mode==9)fCutSelectionRead = Form("%s%s_%s", fEventCutSelection.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
        else if (mode==0)fCutSelectionRead = Form("%s_%s_%s", fEventCutSelection.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
        cout << fCutSelectionRead.Data() << endl;
    }
    
    if(mode == 4 && optionMC.CompareTo("kTRUE") == 0){
        cout << "changing trigger for MC (no EMC trigger in MC)" << endl;
        TString fEventCutSelectionTriggerRead = fEventCutSelection(3,fEventCutSelection.Length()-5);
        cout << fEventCutSelectionTriggerRead.Data() << endl;
        if(fEventCutSelectionTriggerRead.CompareTo("83")==0 || fEventCutSelectionTriggerRead.CompareTo("85")==0 || fEventCutSelectionTriggerRead.CompareTo("93")==0 || fEventCutSelectionTriggerRead.CompareTo("95")==0 || fEventCutSelectionTriggerRead.CompareTo("52")==0 || fEventCutSelectionTriggerRead.CompareTo("51")==0 || fEventCutSelectionTriggerRead.CompareTo("00")==0){
            fEventCutSelection.Replace(GetEventSelectSpecialTriggerCutPosition (),2,"03");
            fEventCutSelectionRead = fEventCutSelection;
            fGammaCutSelectionRead = fGammaCutSelection;
            fMesonCutSelectionRead = fMesonCutSelection;
            cout << fEventCutSelection.Data() << endl;
            fCutSelectionRead = Form("%s_%s_%s",fEventCutSelection.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
        }
        if(file.Contains("LHC12i3") && !file.Contains("LHC12f1a")){
            fEventCutSelection.Replace(GetEventRejectExtraSignalsCutPosition(),1,"2");
            fEventCutSelectionRead = fEventCutSelection;
            fGammaCutSelectionRead = fGammaCutSelection;
            fMesonCutSelectionRead = fMesonCutSelection;
            cout << fEventCutSelection.Data() << endl;
            fCutSelectionRead = Form("%s_%s_%s",fEventCutSelection.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
        }
    }

    
    StyleSettingsThesis(Suffix);
    SetPlotStyle();
        
    fEnergyFlag = optionEnergy;
    fPrefix=meson;

    fPeriodFlag = optionPeriod;
    fdirectphoton = directphotonPlots;

    TString outputDir = Form("%s/%s/%s/ExtractSignal",cutSelection.Data(),optionEnergy.Data(),Suffix.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    
    cout<<"Pictures are saved as "<< Suffix.Data()<< endl;
    fdate = ReturnDateString();
    
    //****************************** Specification of collision system ************************************************
    TString textProcess = ReturnMesonString (fPrefix);
    if(textProcess.CompareTo("") == 0 ){
        cout << "Meson unknown" << endl;
        return ;
    }
    
    fTextMeasurement = Form("%s #rightarrow #gamma#gamma", textProcess.Data());
    fCollisionSystem = ReturnFullCollisionsSystem(fEnergyFlag);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    fDetectionProcess = ReturnFullTextReconstructionProcess(fMode);
    
    //****************************** Choice of Fitting procedure ******************************************************
    if(optionCrystalBall.CompareTo("CrystalBall") == 0){// means we want to plot values for the pi0
        fCrysFitting=1;
        cout << "CrystalBall fit chosen ..." << endl;
    } else   {
        fCrysFitting=0;
        cout << "Gaussian fit chosen ..." << endl;
    }
    
    if(cutSelection.Length() == 0){
        cout<<"ERROR: Cut selection is not set, please do!"<<endl;
        return;
    }

    //***************************** Specification Data/MC ************************************************************
    if(optionMC.CompareTo("kTRUE") == 0){
        fIsMC = 1;
        fPrefix2 = "MC";
    } else {
        fIsMC = 0;
        fPrefix2 = "data";
    }
    
    //**************************** Determine Centrality *************************************************************
    centralityString = GetCentralityString(fEventCutSelection);
    if (centralityString.CompareTo("pp")==0){
        fTextCent = "MinBias";  
    } else {
        fTextCent = Form("%s central", centralityString.Data());
    }
    if (centralityString.CompareTo("pp")!=0 && !centralityString.Contains("0-100%") ){
        fCollisionSystem = Form("%s %s", centralityString.Data(), fCollisionSystem.Data());
    }
    
    cout << "line " << __LINE__ << endl;
    
    //***************************** Initialization of variables according to meson type ******************************
    if(meson.CompareTo("Pi0") == 0){
        Initialize("Pi0",numberOfBins, triggerSet);
    } else if (meson.CompareTo("Eta") == 0) {
        Initialize("Eta",numberOfBins, triggerSet);
    } else if (meson.CompareTo("EtaPrim") == 0) {
        Initialize("EtaPrim",numberOfBins, triggerSet);
    } else if(meson.CompareTo("Pi0EtaBinning") == 0) {
        Initialize("Pi0EtaBinning",numberOfBins, triggerSet);
    } else   {
        cout<<"ERROR: First argument in the ExtractSignal(....) has to be either Pi0 or Eta or Pi0EtaBinning  or EtaPrim"<<endl;
        return;
    }
    
    
    cout << "Integration window normal: "<< fMesonIntDeltaRange[0] << "\t" << fMesonIntDeltaRange[1] << endl;
    cout << "Integration window narrow: "<< fMesonIntDeltaRangeNarrow[0] << "\t" << fMesonIntDeltaRangeNarrow[1] << endl;
    cout << "Integration window wide: "<< fMesonIntDeltaRangeWide[0] << "\t" << fMesonIntDeltaRangeWide[1] << endl;
    
    cout << "line " << __LINE__ << endl;
        
    //************************* Start of Main routine ***************************************************************
    const char* fFileErrLogDatname = Form("%s/%s/%s_%s_FileErrLog%s_%s.dat",cutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelectionRead.Data());
    fFileErrLog.open(fFileErrLogDatname, ios::out);

    TFile f(file.Data());
    
    TString nameMainDir = "";
    if (mode == 9 || mode == 0) nameMainDir = "GammaConvV1";
    else if (mode == 2 || mode == 3) nameMainDir = "GammaConvCalo";
    else if (mode == 4 || mode == 5) nameMainDir = "GammaCalo";
    
    TList *TopDir =(TList*)f.Get(nameMainDir.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }
    
    TList *HistosGammaConversion       = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if(HistosGammaConversion == NULL){
        cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in File"<<endl;
        return;
    }
    if (meson.CompareTo("Pi0") == 0 || meson.CompareTo("Pi0EtaBinning") == 0 ){
        SetCorrectMCHistogrammNames("Pi0");
    } else if (meson.CompareTo("Eta") == 0 ){
        SetCorrectMCHistogrammNames("Eta");
    }

    TList *ESDContainer                 = (TList*) HistosGammaConversion->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));
    TList *BackgroundContainer          = (TList*) HistosGammaConversion->FindObject(Form("%s Back histograms",fCutSelectionRead.Data()));
    TList *MotherContainer              = (TList*) HistosGammaConversion->FindObject(Form("%s Mother histograms",fCutSelectionRead.Data()));
    if (fMode == 2 || fMode == 3 ){
        TList *ClusterContainer             = (TList*) HistosGammaConversion->FindObject(Form("%s Cluster Output",fCutSelectionRead.Data()));
        if (ClusterContainer){
            fHistoClustersPt                = (TH1D*)ClusterContainer->FindObject("ClusGamma_Pt");
            cout << "line " << __LINE__ << endl;
            fHistoClustersOverlapHeadersPt  = (TH1D*)ClusterContainer->FindObject("ClusGammaOverlapHeaders_Pt");
            cout << "line " << __LINE__ << endl;
            TH2F* fHistoTrue2DGammaDCClusPt = (TH2F*)ClusterContainer->FindObject(ObjectNameDCGammaClusPt.Data());
            cout << "line " << __LINE__ << endl;
            if (fHistoTrue2DGammaDCClusPt!=NULL) fEnableDCCluster= kTRUE;
            if (fEnableDCCluster){
                fHistoTrue2DGammaDCClusPt->Sumw2();
                fHistoTrueGammaDCClusPt             = (TH1D*)fHistoTrue2DGammaDCClusPt->ProjectionX("TrueClusGamma_Pt",0,-1,"e");
                cout << "line " << __LINE__ << endl;
                cout << "Cluster DC found " << endl;
                fHistoTrueGammaClusMultipleCount    = (TH1F*)ClusterContainer->FindObject(ObjectNameGammaClusMultipleCount.Data());
                cout << "line " << __LINE__ << endl;
                fHistoTrueGammaClusPt               = (TH1F*)ClusterContainer->FindObject("TrueClusGamma_Pt");
                cout << "line " << __LINE__ << endl;
            }
        }
    }
    if ( fMode == 4  || fMode == 5 ){
            fHistoClustersPt                = (TH1D*)ESDContainer->FindObject("ClusGamma_Pt");
            fHistoClustersOverlapHeadersPt  = (TH1D*)ESDContainer->FindObject("ClusGammaOverlapHeaders_Pt");
        
    }
    
    cout << fMesonCutSelectionRead.Data() << endl;
    cout << fGammaCutSelectionRead.Data() << endl;   
    fNumberOfGoodESDTracks              = (TH1D*)ESDContainer->FindObject("GoodESDTracks");
    fEventQuality                       = (TH1D*)ESDContainer->FindObject("NEvents");
        
    TString rapidityRange;
    fYMaxMeson                          = ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
    fBackgroundMultNumber               = ReturnBackgroundMult(fMesonCutSelection);
        
    TString ObjectNameESD               = "ESD_Mother_InvMass_Pt";
    TString ObjectNameBck               = "ESD_Background_InvMass_Pt";
    
    cout << "line " << __LINE__ << endl;
    fGammaGammaInvMassVSPt              = (TH2D*)ESDContainer->FindObject(ObjectNameESD.Data());
    fGammaGammaInvMassVSPt->Sumw2();
    fBckInvMassVSPt                     = (TH2D*)ESDContainer->FindObject(ObjectNameBck.Data());
    fBckInvMassVSPt->Sumw2();
    cout << "line " << __LINE__ << endl;
    
    const char* FileDataLogname         = Form("%s/%s/%s_%s_EffiCheck_RAWDATA%s_%s.dat", cutSelection.Data(), fEnergyFlag.Data(), fPrefix.Data(), fPrefix2.Data(), fPeriodFlag.Data(),
                                        fCutSelectionRead.Data());
    fFileDataLog.open(FileDataLogname, ios::out);

    if(UseTHnSparse) ProduceBckProperWeighting(BackgroundContainer,MotherContainer,UseTHnSparse);
    else ProduceBckProperWeighting(ESDContainer,ESDContainer,UseTHnSparse);
    cout << "line " << __LINE__ << endl;
    
    if(fIsMC){
        TList *MCContainer              = (TList*)HistosGammaConversion->FindObject(Form("%s MC histograms",fCutSelectionRead.Data()));
        TList *TrueConversionContainer  = (TList*)HistosGammaConversion->FindObject(Form("%s True histograms",fCutSelectionRead.Data()));
        cout << "line " << __LINE__ << endl;
        
        if( fMesonId == 111){
            fHistoMCMesonPtWithinAcceptance     = (TH1D*)MCContainer->FindObject(ObjectNameMCPi0Acc.Data());
            cout << "line " << __LINE__ << endl;
            if ( fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5 ){
                fHistoMCMesonPtWithinAcceptanceWOWeights    = (TH1D*)MCContainer->FindObject(ObjectNameMCPi0AccWOWeights.Data());
            }
            fHistoMCMesonPt                             = (TH1D*)MCContainer->FindObject(ObjectNameMCPi0.Data());   // Not the best; better having a 2D Pt_vs_Rapid in case we change limits
            fHistoMCMesonPtWOWeights                    = (TH1D*)MCContainer->FindObject(ObjectNameMCPi0WOWeights.Data());
        }
        if( fMesonId == 221){
            fHistoMCMesonPtWithinAcceptance     = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaAcc.Data());
//          if ( fMode == 2 || fMode == 3  ){
                fHistoMCMesonPtWithinAcceptanceWOWeights    = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaAccWOWeights.Data());
//          }
        cout << "line " << __LINE__ << endl;
            fHistoMCMesonPt                     = (TH1D*)MCContainer->FindObject(ObjectNameMCEta.Data());   // Not the best; better having a 2D Pt_vs_Rapid in case we change limits
            fHistoMCMesonPtWOWeights            = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaWOWeights.Data());
        }
        fHistoMCMesonPt->Sumw2();
        fHistoMCMesonPtWithinAcceptance->Sumw2();

        cout << "line " << __LINE__ << endl;
    
        if (fHistoMCMesonPtWithinAcceptanceWOWeights){ 
            cout << "found: " << ObjectNameMCPi0AccWOWeights.Data() << endl;
            fHistoMCMesonPtWithinAcceptanceWOWeights->Sumw2();
        }
        if (fHistoMCMesonPtWOWeights){
            fHistoMCMesonPtWeights              = (TH1D*)fHistoMCMesonPtWOWeights->Clone("WeightsMeson");
            fHistoMCMesonPtWeights->Divide(fHistoMCMesonPt,fHistoMCMesonPtWOWeights, 1.,1.,"B");
        }   

        if (fMode == 4 || fMode == 5){
            TH2F* fHistoTrue2DGammaDCClusPt     = (TH2F*)TrueConversionContainer->FindObject(ObjectNameDCGammaClusPt.Data());
            if (fHistoTrue2DGammaDCClusPt!=NULL) fEnableDCCluster= kTRUE;
            if (fEnableDCCluster){
                fHistoTrue2DGammaDCClusPt->Sumw2();
                fHistoTrueGammaDCClusPt         = (TH1D*)fHistoTrue2DGammaDCClusPt->ProjectionX("TrueClusGamma_Pt",0,-1,"e");
                cout << "line " << __LINE__ << endl;
                cout << "Cluster DC found " << endl;
                fHistoTrueGammaClusMultipleCount    = (TH1F*)TrueConversionContainer->FindObject(ObjectNameGammaClusMultipleCount.Data());
                fHistoTrueGammaClusPt               = (TH1F*)TrueConversionContainer->FindObject("TrueClusGamma_Pt");
            }
        }
        
        fHistoTrueMesonInvMassVSPt                  = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrue.Data());
        fHistoTrueMesonInvMassVSPt->Sumw2();
        fHistoTrueFullMesonInvMassVSPt              = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueFull.Data());
        fHistoTrueFullMesonInvMassVSPt->Sumw2();
        fHistoTrueMesonInvMassVSPtWOWeights         = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueWOWeights.Data());
        fHistoTrueMesonInvMassVSPtWOWeights->Sumw2();
        fProfileTrueMesonInvMassVSPtWeights         = (TProfile2D*)TrueConversionContainer->FindObject(ObjectNameProfileWeights.Data());
        fProfileTrueMesonInvMassVSPtWeights->Sumw2();
        fHistoTrueMesonInvMassVSPtReweighted        = (TH2D*)fHistoTrueMesonInvMassVSPtWOWeights->Clone("Reweighted");
        fHistoTrueMesonInvMassVSPtReweighted->Sumw2();
        fHistoTrueMesonInvMassVSPtReweighted->Multiply(fProfileTrueMesonInvMassVSPtWeights);
            
        cout << ObjectNameTrue.Data() << endl;
        FillMassMCTrueMesonHistosArray(fHistoTrueMesonInvMassVSPt);
        FillMassMCTrueFullMesonHistosArray(fHistoTrueFullMesonInvMassVSPt);
        FillMassMCTrueReweightedMesonHistosArray(fHistoTrueMesonInvMassVSPtReweighted);
        FillMassMCTrueUnweightedMesonHistosArray(fHistoTrueMesonInvMassVSPtWOWeights);

        cout << "line " << __LINE__ << endl;
        cout << fMode << endl;
//         return;
        
        if (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5){
            fHistoTrueMesonCaloPhotonInvMassVSPt                = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloPhoton.Data());
            if (fHistoTrueMesonCaloPhotonInvMassVSPt==NULL) fAdvancedMesonQA = kFALSE;
            else fAdvancedMesonQA = kTRUE;
            cout << fAdvancedMesonQA << endl;
            // temp fix 
//             fAdvancedMesonQA = kFALSE;
        }
        
        if (fMode == 0 || fMode == 1 || fMode == 9){
            fHistoTrueContBckInvMassVSPt                        = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueContBck.Data());
            if (fHistoTrueContBckInvMassVSPt==NULL) fAdvancedMesonQA = kFALSE;
        }
                
        if (fAdvancedMesonQA) {
            if (fMode == 2 || fMode == 3){
                fHistoTrueMesonCaloPhotonInvMassVSPt                = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloPhoton.Data());
                FillMassMCTrueMesonCaloPhotonHistosArray(fHistoTrueMesonCaloPhotonInvMassVSPt);
                fHistoTrueMesonCaloConvPhotonInvMassVSPt            = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloConvPhoton.Data());
                FillMassMCTrueMesonCaloConvPhotonHistosArray(fHistoTrueMesonCaloConvPhotonInvMassVSPt);
                fHistoTrueMesonMergedClusterInvMassVSPt               = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloMerged.Data());
                FillMassMCTrueMesonCaloMergedClusterHistosArray(fHistoTrueMesonMergedClusterInvMassVSPt);
                fHistoTrueMesonMergedClusterPartConvInvMassVSPt     = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloMergedPartConv.Data());
                FillMassMCTrueMesonCaloMergedClusterPartConvHistosArray(fHistoTrueMesonMergedClusterPartConvInvMassVSPt);
            } else if (fMode == 4 || fMode == 5){
                fHistoTrueMesonCaloPhotonInvMassVSPt                = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloPhoton.Data());
                FillMassMCTrueMesonCaloPhotonHistosArray(fHistoTrueMesonCaloPhotonInvMassVSPt);
                fHistoTrueMesonCaloConvPhotonInvMassVSPt            = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloConvPhoton.Data());
                FillMassMCTrueMesonCaloConvPhotonHistosArray(fHistoTrueMesonCaloConvPhotonInvMassVSPt);
                fHistoTrueMesonMixedCaloConvPhotonInvMassVSPt       = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueMixedCaloConvPhoton.Data());
                FillMassMCTrueMesonMixedCaloConvPhotonHistosArray(fHistoTrueMesonMixedCaloConvPhotonInvMassVSPt);
                fHistoTrueMesonMergedClusterInvMassVSPt             = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloMerged.Data());
                FillMassMCTrueMesonCaloMergedClusterHistosArray(fHistoTrueMesonMergedClusterInvMassVSPt);
                fHistoTrueMesonMergedClusterPartConvInvMassVSPt     = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloMergedPartConv.Data());
                FillMassMCTrueMesonCaloMergedClusterPartConvHistosArray(fHistoTrueMesonMergedClusterPartConvInvMassVSPt);
            } else {
                fHistoTrueContBckInvMassVSPt                        = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueContBck.Data());
                FillMassMCTrueContBckHistosArray(fHistoTrueContBckInvMassVSPt);
                fHistoTrueGGBckInvMassVSPt                          = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueGGBck.Data());
                FillMassMCTrueGGBckHistosArray(fHistoTrueGGBckInvMassVSPt);
                fHistoTrueAllBckInvMassVSPt                         = (TH2D*)fHistoTrueGGBckInvMassVSPt->Clone(ObjectNameTrueAllBck.Data());
                fHistoTrueAllBckInvMassVSPt->Sumw2();
                fHistoTrueAllBckInvMassVSPt->Add(fHistoTrueContBckInvMassVSPt);
                FillMassMCTrueAllBckHistosArray(fHistoTrueAllBckInvMassVSPt);
            }
            fHistoYieldK0sWithPi0DaughterRec                        = (TH1D*)TrueConversionContainer->FindObject(ObjectNameK0sRecPi0.Data());
            if(fHistoYieldK0sWithPi0DaughterRec) fHistoYieldK0sWithPi0DaughterRec->Sumw2();
            fHistoYieldLambdaWithPi0DaughterRec                     = (TH1D*)TrueConversionContainer->FindObject(ObjectNameLambdaRecPi0.Data());
            if(fHistoYieldLambdaWithPi0DaughterRec) fHistoYieldLambdaWithPi0DaughterRec->Sumw2();
        }
    
        if (meson.Contains("Pi0")){
            fHistoTrueSecMesonInvMassVSPt                           = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueSec.Data());
            FillMassMCTrueSecMesonHistosArray(fHistoTrueSecMesonInvMassVSPt);
            fHistoTrueSecFromK0SMesonInvMassVSPt                    = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueSecFromK0S.Data());
            FillMassMCTrueSecFromK0SMesonHistosArray(fHistoTrueSecFromK0SMesonInvMassVSPt);   
            fHistoTrueSecFromLambdaMesonInvMassVSPt                 = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueSecFromLambda.Data());
            FillMassMCTrueSecFromLambdaMesonHistosArray(fHistoTrueSecFromLambdaMesonInvMassVSPt);   
        }

        cout << "line " << __LINE__ << endl;
        fHistoTrueMesonDCInvMassVSPt                                = (TH2D*)TrueConversionContainer->FindObject(ObjectNameDCMesonInvMassPt.Data());
        if (fHistoTrueMesonDCInvMassVSPt!= NULL) fEnableDCMeson = kTRUE;
        cout << "line " << __LINE__ << endl;
        if (fEnableDCMeson){
            FillMassMCTrueMesonDCHistosArray(fHistoTrueMesonDCInvMassVSPt);
            fHistoTrueMesonMultipleCount = (TH1F*) TrueConversionContainer->FindObject(ObjectNameMesonMultipleCount.Data());
        }
        
        
        cout << "line " << __LINE__ << endl;
    }
    
    
    fMesonMassExpect                            = TDatabasePDG::Instance()->GetParticle(fMesonId)->Mass();
    if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 || fEnergyFlag.CompareTo("pPb_5.023TeV") == 0){
        fNEvents        = fEventQuality->GetBinContent(1);
    } else {
        fNEvents        =  GetNEvents(fEventQuality);
    }
    
    TH1D *fBck                                  = (TH1D*)fBckInvMassVSPt->ProjectionX("ESD_Background_InvMass");
    fBck->Sumw2();
    TH1D *fGammaGamma                           = (TH1D*)fGammaGammaInvMassVSPt->ProjectionX("ESD_Mother_InvMass");
    fGammaGamma->Sumw2();
    fPeakPosAlpha01                             = (TH2D*)ESDContainer->FindObject("ESD_Mother_InvMass_vs_E_alpha");
    if (fPeakPosAlpha01==NULL) fPeakPosAlpha01  = (TH2D*)ESDContainer->FindObject("ESD_Mother_InvMass_vs_Pt_Alpha");
    
    cout<< "The mass of the meson is: "<< fMesonMassExpect<< " Events analysed: "<< fNEvents<< endl;
    cout << "line " << __LINE__ << endl;
    // Process the 1D invariant mass histos
    fGammaGamma->SetTitle(Form("%s %s",fGammaGamma->GetTitle(),fCutSelection.Data()));
    cout << "line " << __LINE__ << endl;
    fGammaGamma->Rebin(fNRebinGlobal);
    //fGammaGamma->Scale(1./fNRebinGlobal);
    fBck->Rebin(fNRebinGlobal);
    //fBck->Scale(1./fNRebinGlobal);
    ProcessEM( fGammaGamma , fBck, fBGFitRange);
    fHistoMappingBackNormInvMass    = fBckNorm;
    fHistoMappingBackNormInvMass->Sumw2();
    fHistoMappingSignalInvMass      = fSignal;
    fHistoMappingSignalInvMass->Sumw2();
    
    fGammaGamma->DrawCopy();
    fHistoMappingBackNormInvMass->DrawCopy("same");
    fHistoMappingSignalInvMass->DrawCopy("same");
    
    // Function to Project the 2D histos InvariantMass VS Pt into Invariant Mass spectrum
    FillMassHistosArray(fGammaGammaInvMassVSPt,fPeakPosAlpha01);
    cout << "line " << __LINE__ << endl;
    
    ProcessEM( fMesonFullPtSignal, fMesonFullPtBackground, fBGFitRange);
    fMesonFullPtBackNorm            = fBckNorm;
    
    ProcessEM( fFittingHistMidPtSignal, fFittingHistMidPtBackground, fBGFitRange);
    fFittingHistMidPtSignalSub      = fSignal;
    if(fCrysFitting==0){
        fFileErrLog << "Using exp fit"<<endl;
        FitSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub, fMesonIntDeltaRange,200,kTRUE);
        fFitSignalInvMassMidPt      = fFitReco;
//         if (fMode == 4){
//             GausFitSubtractedInvMassInPtBinsNew(fFittingHistMidPtSignalSub, fMesonIntDeltaRange,200,kTRUE,"SinglefitfunctionMidPt",kFALSE);
//             fFitSignalInvMassMidPt2         = fFitReco;
//         }
    } else {
        fFileErrLog << "Using Crystal Ball function"<<endl;
        FitCBSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub, fMesonIntDeltaRange,200,kTRUE,"SinglefitfunctionMidPt",kFALSE);
        fFitSignalInvMassMidPt      = fFitReco;
    }

    if (fIsMC){
        TH1D* fHistoMappingTrueMesonInvMassPtMidPt= NULL;
        fHistoMappingTrueMesonInvMassPtMidPt= new TH1D("TrueMassMidPt", "TrueMassMidPt", fHistoTrueMesonInvMassVSPtReweighted->GetNbinsX(), 0., 
                                                    fHistoTrueMesonInvMassVSPtReweighted->GetXaxis()->GetBinUpEdge(fHistoTrueMesonInvMassVSPtReweighted->GetNbinsX()));
        fHistoMappingTrueMesonInvMassPtMidPt->Sumw2();
        Int_t startBin  = fHistoTrueMesonInvMassVSPtReweighted->GetYaxis()->FindBin(fMidPt[0]+0.001);
        Int_t endBin    = fHistoTrueMesonInvMassVSPtReweighted->GetYaxis()->FindBin(fMidPt[1]-0.001);
        fHistoTrueMesonInvMassVSPtReweighted->ProjectionX("TrueMassMidPt",startBin,endBin,"e");
            
        fHistoMappingTrueMesonInvMassPtMidPt=(TH1D*)gDirectory->Get("TrueMassMidPt");
        fHistoMappingTrueMesonInvMassPtMidPt->Rebin(fNRebin[5]);
        fHistoMappingTrueMesonInvMassPtMidPt->SetLineWidth(1);
        fHistoMappingTrueMesonInvMassPtMidPt->SetLineColor(2);
        FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtMidPt, fMesonIntDeltaRange,50,kTRUE);
        
        if(fCrysFitting==0){
            FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtMidPt, fMesonIntDeltaRange,50,kTRUE);
        } else {
            FitCBSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtMidPt, fMesonIntDeltaRange,50,kFALSE,"CBFitFuncTrueMidPt",kTRUE);
        }
        
    }
    
    
    TString nameMesonFittingMidPt= Form("%s/%s_%s_MesonSubtractedFittingMidPt%s_%s_%02d.%s",outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(), 200,
                                Suffix.Data());
    TString nameCanvasFittingMidPt= "MesonCanvasSubtractedFittingMidPt";
    TString namePadFittingMidPt= "MesonPadSubtractedFittingMidPt";
    //    
    TString fDecayChannel = "#gamma#gamma";
    delete fMidPt;

    if(fEnergyFlag.CompareTo("7TeV") == 0 ){
        if(fPrefix.CompareTo("Pi0") == 0){
            for (Int_t iPt=fStartPtBin;iPt<fNBinsPeakPt;iPt++){
                fFitPeakPosPtBin[iPt] =0x00;
                if(fCrysFitting==0){
                    fFileErrLog << "Using exp fit"<<endl;
                    FitPeakPosInvMassInPtBins(fHistoMappingPeakPosInvMassPtBin[iPt],iPt,kTRUE);
                    fFitPeakPosPtBin[iPt]       = fFitReco;
//                     if (fMode == 4){
//                         GausFitSubtractedInvMassInPtBinsNew(fHistoMappingPeakPosInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncPeakPos%02d",iPt),kFALSE);
//                         fFitPeakPosPtBin2[iPt]      = fFitReco;
//                     }
                } else {
                    fFileErrLog << "Using Crystal Ball function"<<endl;
                    FitCBSubtractedInvMassInPtBins(fHistoMappingPeakPosInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncPeakPos%02d",iPt),kFALSE);
                    fFitPeakPosPtBin[iPt]       = fFitReco;
                }
                if (fFitPeakPosPtBin[iPt] !=0x00){
                    fMesonMassPeakPos[iPt]      = fFitPeakPosPtBin[iPt]->GetParameter(1);
                    fMesonMassPeakPosError[iPt] = fFitPeakPosPtBin[iPt]->GetParError(1);
                    CalculateFWHM( fFitPeakPosPtBin[iPt]);
                    fMesonFWHMAlpha01[iPt]      = fFWHMFunc;
                    fMesonFWHMAlpha01Error[iPt] = fFWHMFuncError;
                } else {
                    fMesonMassPeakPos[iPt]      = 0.;
                    fMesonMassPeakPosError[iPt] = 0.;
                    fMesonFWHMAlpha01[iPt]      = 0.;
                    fMesonFWHMAlpha01Error[iPt] = 0.;
                }
            }
        }
    }
    
    
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){ // BEGIN ANALYSIS for each Pt bin


        cout << "---------------------------------------------------------------------------------" << endl;
        cout << "Begin Analysis Pt Bin " << iPt <<endl;
        cout << "---------------------------------------------------------------------------------" << endl;
        // Function to subtract GG minus Bck
        fFileDataLog << "---------------------------------------------------------------------------------" << endl;
        fFileDataLog << "----------------------------------new pT bin ------------------------------------" << endl;
        fFileDataLog << "---------------------------------------------------------------------------------" << endl;
        
        ProcessEM( fHistoMappingGGInvMassPtBin[iPt], fHistoMappingBackInvMassPtBin[iPt], fBGFitRange);
        fHistoMappingSignalInvMassPtBin[iPt] = fSignal;
        fHistoMappingBackNormInvMassPtBin[iPt] = fBckNorm;

        fHistoMappingSignalInvMassPtBin[iPt]->SetName(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",iPt));
        
        // Fitting the subtracted spectra
        fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "normal range/right normalization" << endl;

        fFitSignalInvMassPtBin[iPt]=0x00;
        if(fCrysFitting==0){
            fFileErrLog << "Using exp fit"<<endl;
            fFileDataLog << "Subtracted mixed event" << endl;
            FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
            fFitSignalInvMassPtBin[iPt]             = fFitReco;
            fFitSignalPeakPosInvMassPtBin[iPt]      = fFitGausExp;
            fFitBckInvMassPtBin[iPt]                = fFitLinearBck;
            fMesonYieldsResidualBckFunc[iPt]        = fIntLinearBck;
            fMesonYieldsResidualBckFuncError[iPt]   = fIntLinearBckError;
//             if (fMode == 4){
//                 fFitSignalInvMassPtBin2[iPt]                        = fFitReco;
//                 fFitSignalPeakPosInvMassPtBin2[iPt]                 = fFitGausExp;
//                 fFitBckInvMassPtBin2[iPt]                           = fFitLinearBck;
//                 GausFitSubtractedInvMassInPtBinsNew(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("GausFitFuncNormalBin%02d",iPt),kFALSE);
//                 fHistoMappingSignalRemainingBGSubInvMassPtBin[iPt]  = (TH1D*)fCopySignal->Clone(Form("histoSignalRemainingBGSubtractedBin%02d",iPt));
//                 fHistoMappingRemainingBGInvMassPtBin[iPt]           = (TH1D*)fCopyOnlyBG->Clone(Form("histoRemainingBGBin%02d",iPt));
//                 fFitSignalInvMassPtBin[iPt]                         = fFitReco;
//                 fFitRemainingBGInvMassPtBin[iPt]                    = fFitLinearBck;
//                 fFitSignalPeakPosInvMassPtBin[iPt]                  = fFitGausExp;
//                 fFitBckInvMassPtBin[iPt]                            = fFitLinearBck;
//             }
        } else {
            fFileErrLog << "Using Crystal Ball function"<<endl;
            FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncNormalBin%02d",iPt),kFALSE);
            fHistoMappingSignalRemainingBGSubInvMassPtBin[iPt]  = (TH1D*)fCopySignal->Clone(Form("histoSignalRemainingBGSubtractedBin%02d",iPt));
            fHistoMappingRemainingBGInvMassPtBin[iPt]           = (TH1D*)fCopyOnlyBG->Clone(Form("histoRemainingBGBin%02d",iPt));
            fFitSignalInvMassPtBin[iPt]                         = fFitReco;
            fFitRemainingBGInvMassPtBin[iPt]                    = fFitLinearBck;
            fFitSignalPeakPosInvMassPtBin[iPt]                  = fFitGausExp;
            fFitBckInvMassPtBin[iPt]                            = fFitLinearBck;
            fMesonYieldsResidualBckFunc[iPt]                    = 0;
            fMesonYieldsResidualBckFuncError[iPt]               = 0;
            
        }

        //Get FWHM
        CalculateFWHM( fFitSignalInvMassPtBin[iPt]);
        fMesonFWHM[iPt] = fFWHMFunc;
        fMesonFWHMError[iPt] = fFWHMFuncError;

        if (fFitSignalInvMassPtBin[iPt] !=0x00){
            fMesonMass[iPt]                 = fFitSignalInvMassPtBin[iPt]->GetParameter(1);
            fMesonMassError[iPt]            = fFitSignalInvMassPtBin[iPt]->GetParError(1);
            
            fMesonLambdaTailpar[iPt]        = fFitSignalInvMassPtBin[iPt]->GetParameter(3);
            fMesonLambdaTailparError[iPt]   = fFitSignalInvMassPtBin[iPt]->GetParError(3);

            fMesonCurIntRange[0]            = fMesonMass[iPt] + fMesonIntDeltaRange[0];
            fMesonCurIntRangeWide[0]        = fMesonMass[iPt] + fMesonIntDeltaRangeWide[0];
            fMesonCurIntRangeNarrow[0]      = fMesonMass[iPt] + fMesonIntDeltaRangeNarrow[0];
            fMesonCurIntRange[1]            = fMesonMass[iPt] + fMesonIntDeltaRange[1];
            fMesonCurIntRangeWide[1]        = fMesonMass[iPt] + fMesonIntDeltaRangeWide[1];
            fMesonCurIntRangeNarrow[1]      = fMesonMass[iPt] + fMesonIntDeltaRangeNarrow[1];
        } else {
            fMesonMass[iPt]                 = fMesonMassExpect;
            fMesonMassError[iPt]            = 0.;
            fMesonCurIntRange[0]            = fMesonMassExpect + fMesonIntDeltaRange[0];
            fMesonCurIntRangeWide[0]        = fMesonMassExpect + fMesonIntDeltaRangeWide[0];
            fMesonCurIntRangeNarrow[0]      = fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
            fMesonCurIntRange[1]            = fMesonMassExpect + fMesonIntDeltaRange[1];
            fMesonCurIntRangeWide[1]        = fMesonMassExpect + fMesonIntDeltaRangeWide[1];
            fMesonCurIntRangeNarrow[1]      = fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];
        }


        fMassWindowHigh[iPt]        = fMesonCurIntRange[1];
        fMassWindowLow[iPt]         = fMesonCurIntRange[0];
        fMassWindowWideHigh[iPt]    = fMesonCurIntRangeWide[1];
        fMassWindowWideLow[iPt]     = fMesonCurIntRangeWide[0];
        fMassWindowNarrowHigh[iPt]  = fMesonCurIntRangeNarrow[1];
        fMassWindowNarrowLow[iPt]   = fMesonCurIntRangeNarrow[0];
        
//         fMassWindowHigh[iPt]         = fMesonMass[iPt] + fMesonIntDeltaRange[1];
//         fMassWindowLow[iPt]          = fMesonMass[iPt] + fMesonIntDeltaRange[0]; 
//         fMassWindowWideHigh[iPt]     = fMesonMass[iPt] + fMesonIntDeltaRangeWide[1];
//         fMassWindowWideLow[iPt]      = fMesonMass[iPt] + fMesonIntDeltaRangeWide[0];
//         fMassWindowNarrowHigh[iPt]   = fMesonMass[iPt] + fMesonIntDeltaRangeNarrow[1];
//         fMassWindowNarrowLow[iPt]    = fMesonMass[iPt] + fMesonIntDeltaRangeNarrow[0];
        
        FitSubtractedPureGaussianInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt],iPt);
        fFitSignalGaussianInvMassPtBin[iPt] = fFitReco;
        if (fFitSignalGaussianInvMassPtBin[iPt] !=0x00){
            fMesonMassGaussian[iPt]         = fFitSignalGaussianInvMassPtBin[iPt]->GetParameter(1);
            fMesonMassGaussianError[iPt]    = fFitSignalGaussianInvMassPtBin[iPt]->GetParError(1);
            fMesonWidthGaussian[iPt]        = fFitSignalGaussianInvMassPtBin[iPt]->GetParameter(2);
            fMesonWidthGaussianError[iPt]   = fFitSignalGaussianInvMassPtBin[iPt]->GetParError(2);
        } else {
            fMesonMassGaussian[iPt]         = 0.;
            fMesonMassGaussianError[iPt]    = 0.;
            fMesonWidthGaussian[iPt]        = 0.;
            fMesonWidthGaussianError[iPt]   = 0.;
        }

        
        if (fCrysFitting == 0){
            IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRange);
            fGGYields[iPt]                  = fYields;
            fGGYieldsError[iPt]             = fYieldsError;
            IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRangeWide);
            fGGYieldsWide[iPt]              = fYields;
            fGGYieldsWideError[iPt]         = fYieldsError;
            IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRangeNarrow);
            fGGYieldsNarrow[iPt]            = fYields;
            fGGYieldsNarrowError[iPt]       = fYieldsError;
    
            // Integrate the bck histo
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRange);
            fBckYields[iPt]                 = fYields;
            fBckYieldsError[iPt]            = fYieldsError;
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRangeWide);
            fBckYieldsWide[iPt]             = fYields;
            fBckYieldsWideError[iPt]        = fYieldsError;
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRangeNarrow);
            fBckYieldsNarrow[iPt]           = fYields;
            fBckYieldsNarrowError[iPt]      = fYieldsError;

            // Integrate the signal histo
            fFileDataLog<< endl <<"Signal histo normal range, right norm "<< fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRange);
            fMesonYields[iPt]               = fYields;
            fMesonYieldsError[iPt]          = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
            fFileDataLog<< endl <<"Signal histo wide range, right norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRangeWide);
            fMesonYieldsWide[iPt]           = fYields;
            fMesonYieldsWideError[iPt]      = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
            fFileDataLog<< endl <<"Signal histo narrow range, right norm" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRangeNarrow);
            fMesonYieldsNarrow[iPt]         = fYields;
            fMesonYieldsNarrowError[iPt]    = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
        } else {
            IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRange);
            fGGYields[iPt]                  = fYields;
            fGGYieldsError[iPt]             = fYieldsError;
            IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRangeWide);
            fGGYieldsWide[iPt]              = fYields;
            fGGYieldsWideError[iPt]         = fYieldsError;
            IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRangeNarrow);
            fGGYieldsNarrow[iPt]            = fYields;
            fGGYieldsNarrowError[iPt]       = fYieldsError;
    
            // Integrate the bck histo
//             if (fMode != 4){
                fHistoMappingBackNormInvMassPtBin[iPt]->Sumw2();
                fHistoMappingBackNormInvMassPtBin[iPt]->Add(fHistoMappingRemainingBGInvMassPtBin[iPt],1);
//             }
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRange);
            
            fBckYields[iPt]                 = fYields;
            fBckYieldsError[iPt]            = fYieldsError;
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRangeWide);
            fBckYieldsWide[iPt]             = fYields;
            fBckYieldsWideError[iPt]        = fYieldsError;
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRangeNarrow);
            fBckYieldsNarrow[iPt]           = fYields;
            fBckYieldsNarrowError[iPt]      = fYieldsError;

            // Integrate the signal histo
            fFileDataLog<< endl <<"Signal histo normal range, right norm "<< fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalRemainingBGSubInvMassPtBin[iPt], fMesonCurIntRange);
            fMesonYields[iPt]               = fYields;
            fMesonYieldsError[iPt]          = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
            fFileDataLog<< endl <<"Signal histo wide range, right norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalRemainingBGSubInvMassPtBin[iPt], fMesonCurIntRangeWide);
            fMesonYieldsWide[iPt]           = fYields;
            fMesonYieldsWideError[iPt]      = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
            fFileDataLog<< endl <<"Signal histo narrow range, right norm" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalRemainingBGSubInvMassPtBin[iPt], fMesonCurIntRangeNarrow);
            fMesonYieldsNarrow[iPt]         = fYields;
            fMesonYieldsNarrowError[iPt]    = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
        }


        if(fIsMC){
            fFileDataLog<< endl <<"True histo normal range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            fFitTrueSignalInvMassPtBin[iPt]=0x00;
            if(fCrysFitting==0){
                fFileErrLog << "Using exp fit"<<endl;
                FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);
//                 if(fMode == 4){
//                     GausFitSubtractedInvMassInPtBinsNew(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("GausFitFuncMCTrueBin%02d",iPt),kTRUE);
//                 }
            } else {
                fFileErrLog << "Using Crystal Ball function"<<endl;
                FitCBSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncMCTrueBin%02d",iPt),kTRUE);
            }

            //   FitSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);
            if (fHistoMappingTrueMesonInvMassPtBins[iPt]->GetEntries() !=0){
                fFitTrueSignalInvMassPtBin[iPt]     = fFitReco;
                if (fFitTrueSignalInvMassPtBin[iPt] != 0x00){
                    fMesonTrueMass[iPt]             = fFitTrueSignalInvMassPtBin[iPt]->GetParameter(1);
                    fMesonTrueMassError[iPt]        = fFitTrueSignalInvMassPtBin[iPt]->GetParError(1);
                    
                    fMesonLambdaTailMCpar[iPt]      = fFitTrueSignalInvMassPtBin[iPt]->GetParameter(3);
                    fMesonLambdaTailMCparError[iPt] = fFitTrueSignalInvMassPtBin[iPt]->GetParError(3);

                    CalculateFWHM(fFitTrueSignalInvMassPtBin[iPt]);
                    fMesonTrueFWHM[iPt]             = fFWHMFunc;
                    fMesonTrueFWHMError[iPt]        = fFWHMFuncError;
                    fFileDataLog << "TrueFWHM \t" << fMesonTrueFWHM[iPt] << "\t +-" << fMesonTrueFWHMError[iPt] << endl;
                    fMesonTrueIntRange[0]           = fMesonTrueMass[iPt] + fMesonIntDeltaRange[0];
                    fMesonTrueIntRangeWide[0]       = fMesonTrueMass[iPt] + fMesonIntDeltaRangeWide[0];
                    fMesonTrueIntRangeNarrow[0]     = fMesonTrueMass[iPt] + fMesonIntDeltaRangeNarrow[0];
                    fMesonTrueIntRange[1]           = fMesonTrueMass[iPt] + fMesonIntDeltaRange[1] ;
                    fMesonTrueIntRangeWide[1]       = fMesonTrueMass[iPt] + fMesonIntDeltaRangeWide[1];
                    fMesonTrueIntRangeNarrow[1]     = fMesonTrueMass[iPt] + fMesonIntDeltaRangeNarrow[1];
                } else {
                    fMesonTrueMass[iPt]             = 0.;
                    fMesonTrueMassError[iPt]        = 1.;
                    fMesonTrueFWHM[iPt]             = 0.;
                    fMesonTrueFWHMError[iPt]        = 0.;
                    fMesonTrueIntRange[0]           = fMesonMassExpect + fMesonIntDeltaRange[0];
                    fMesonTrueIntRangeWide[0]       = fMesonMassExpect + fMesonIntDeltaRangeWide[0];
                    fMesonTrueIntRangeNarrow[0]     = fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
                    fMesonTrueIntRange[1]           = fMesonMassExpect + fMesonIntDeltaRange[1];
                    fMesonTrueIntRangeWide[1]       = fMesonMassExpect + fMesonIntDeltaRangeWide[1];
                    fMesonTrueIntRangeNarrow[1]     = fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];
                }
            }

            IntegrateHistoInvMassStream( fHistoMappingTrueFullMesonInvMassPtBins[iPt], fIntFixedRange);
            fMesonTrueYieldFixedWindow[iPt]         = fYields;
            fMesonTrueYieldErrorFixedWindow[iPt]    = fYieldsError;

            FitTrueInvMassPureGaussianInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt],iPt);
            if (fHistoMappingTrueMesonInvMassPtBins[iPt]->GetEntries() !=0){
                fFitTrueSignalGaussianInvMassPtBin[iPt] = fFitReco;
                if (fFitTrueSignalGaussianInvMassPtBin[iPt] !=0x00){
                    fMesonTrueMassGaussian[iPt]         = fFitTrueSignalGaussianInvMassPtBin[iPt]->GetParameter(1);
                    fMesonTrueMassGaussianError[iPt]    = fFitTrueSignalGaussianInvMassPtBin[iPt]->GetParError(1);
                    fMesonTrueWidthGaussian[iPt]        = fFitTrueSignalGaussianInvMassPtBin[iPt]->GetParameter(2);
                    fMesonTrueWidthGaussianError[iPt]   = fFitTrueSignalGaussianInvMassPtBin[iPt]->GetParError(2);
                } else {
                    fMesonTrueMassGaussian[iPt]         = 0.;
                    fMesonTrueMassGaussianError[iPt]    = 0.;
                    fMesonTrueWidthGaussian[iPt]        = 0.;
                    fMesonTrueWidthGaussianError[iPt]   = 0.;
                }
            }
            
            fFitTrueSignalInvMassPtReweightedBin[iPt]   = 0x00;
            cout << "line " << __LINE__ << endl;
            if(fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]){
                cout << "Using exp fit"<<endl;
                fFileErrLog << "Using exp fit"<<endl;
                if(fCrysFitting==0){
                    fFileErrLog << "Using exp fit"<<endl;
                    FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtReweightedBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);
//                     if( fMode == 4 ){
//                     GausFitSubtractedInvMassInPtBinsNew(fHistoMappingTrueMesonInvMassPtReweightedBins[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("GausFitFuncMCTrueBinReweighted%02d",iPt),kTRUE);
//                     }
                } else {
                    fFileErrLog << "Using Crystal Ball function"<<endl;
                    FitCBSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtReweightedBins[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncMCTrueBinReweighted%02d",iPt),kTRUE);
                }
                
                if (fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->GetEntries() !=0){
                    fFitTrueSignalInvMassPtReweightedBin[iPt]   = fFitReco;
                    if (fFitTrueSignalInvMassPtReweightedBin[iPt] != 0x00){
                        fMesonTrueMassReweighted[iPt]           = fFitTrueSignalInvMassPtReweightedBin[iPt]->GetParameter(1);
                        fMesonTrueMassReweightedError[iPt]      = fFitTrueSignalInvMassPtReweightedBin[iPt]->GetParError(1);
                        CalculateFWHM(fFitTrueSignalInvMassPtReweightedBin[iPt]);
                        fMesonTrueFWHMReweighted[iPt]           = fFWHMFunc;
                        fMesonTrueFWHMReweightedError[iPt]      = fFWHMFuncError;
                        fMesonTrueIntReweightedRange[0]         = fMesonTrueMassReweighted[iPt] + fMesonIntDeltaRange[0];
                        fMesonTrueIntReweightedRangeWide[0]     = fMesonTrueMassReweighted[iPt] + fMesonIntDeltaRangeWide[0];
                        fMesonTrueIntReweightedRangeNarrow[0]   = fMesonTrueMassReweighted[iPt] + fMesonIntDeltaRangeNarrow[0];
                        fMesonTrueIntReweightedRange[1]         = fMesonTrueMassReweighted[iPt] + fMesonIntDeltaRange[1];
                        fMesonTrueIntReweightedRangeWide[1]     = fMesonTrueMassReweighted[iPt] + fMesonIntDeltaRangeWide[1];
                        fMesonTrueIntReweightedRangeNarrow[1]   = fMesonTrueMassReweighted[iPt] + fMesonIntDeltaRangeNarrow[1];
                    } else {
                        fMesonTrueMassReweighted[iPt]           = 0.;
                        fMesonTrueMassReweightedError[iPt]      = 1.;
                        fMesonTrueFWHMReweighted[iPt]           = 0.;
                        fMesonTrueFWHMReweightedError[iPt]      = 0.;
                        fMesonTrueIntReweightedRange[0]         = fMesonMassExpect + fMesonIntDeltaRange[0];
                        fMesonTrueIntReweightedRangeWide[0]     = fMesonMassExpect + fMesonIntDeltaRangeWide[0];
                        fMesonTrueIntReweightedRangeNarrow[0]   = fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
                        fMesonTrueIntReweightedRange[1]         = fMesonMassExpect + fMesonIntDeltaRange[1];
                        fMesonTrueIntReweightedRangeWide[1]     = fMesonMassExpect + fMesonIntDeltaRangeWide[1];
                        fMesonTrueIntReweightedRangeNarrow[1]   = fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];

                    }
                }
            }

            fFitTrueSignalInvMassPtUnweightedBin[iPt]=0x00;
            cout << "line " << __LINE__ << endl;
            if(fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt]){
                cout << "Using exp fit"<<endl;
                fFileErrLog << "Using exp fit"<<endl;
                if(fCrysFitting==0){
                    fFileErrLog << "Using exp fit"<<endl;
                    FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);
//                     if( fMode == 4 ){
//                     GausFitSubtractedInvMassInPtBinsNew(fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("GausFitFuncMCTrueBinUnweighted%02d",iPt),kTRUE);
//                     }
                } else {
                    fFileErrLog << "Using Crystal Ball function"<<endl;
                    FitCBSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncMCTrueBinUnweighted%02d",iPt),kTRUE);
                }
                
                if (fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt]->GetEntries() !=0){
                    fFitTrueSignalInvMassPtUnweightedBin[iPt] = fFitReco;
                    if (fFitTrueSignalInvMassPtUnweightedBin[iPt] != 0x00){
                    fMesonTrueMassUnweighted[iPt]               = fFitTrueSignalInvMassPtUnweightedBin[iPt]->GetParameter(1);
                    fMesonTrueMassUnweightedError[iPt]          = fFitTrueSignalInvMassPtUnweightedBin[iPt]->GetParError(1);
                    CalculateFWHM(fFitTrueSignalInvMassPtUnweightedBin[iPt]);
                    fMesonTrueFWHMUnweighted[iPt]               = fFWHMFunc;
                    fMesonTrueFWHMUnweightedError[iPt]          = fFWHMFuncError;
                        fMesonTrueIntUnweightedRange[0]         = fMesonTrueMassUnweighted[iPt] + fMesonIntDeltaRange[0];
                        fMesonTrueIntUnweightedRangeWide[0]     = fMesonTrueMassUnweighted[iPt] + fMesonIntDeltaRangeWide[0];
                        fMesonTrueIntUnweightedRangeNarrow[0]   = fMesonTrueMassUnweighted[iPt] + fMesonIntDeltaRangeNarrow[0];
                        fMesonTrueIntUnweightedRange[1]         = fMesonTrueMassUnweighted[iPt] + fMesonIntDeltaRange[1];
                        fMesonTrueIntUnweightedRangeWide[1]     = fMesonTrueMassUnweighted[iPt] + fMesonIntDeltaRangeWide[1];
                        fMesonTrueIntUnweightedRangeNarrow[1]   = fMesonTrueMassUnweighted[iPt] + fMesonIntDeltaRangeNarrow[1];
                    } else {
                        fMesonTrueMassUnweighted[iPt]           = 0.;
                        fMesonTrueMassUnweightedError[iPt]      = 1.;
                        fMesonTrueFWHMUnweighted[iPt]           = 0.;
                        fMesonTrueFWHMUnweightedError[iPt]      = 0.;
                        fMesonTrueIntUnweightedRange[0]         = fMesonMassExpect + fMesonIntDeltaRange[0];
                        fMesonTrueIntUnweightedRangeWide[0]     = fMesonMassExpect + fMesonIntDeltaRangeWide[0];
                        fMesonTrueIntUnweightedRangeNarrow[0]   = fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
                        fMesonTrueIntUnweightedRange[1]         = fMesonMassExpect + fMesonIntDeltaRange[1];
                        fMesonTrueIntUnweightedRangeWide[1]     = fMesonMassExpect + fMesonIntDeltaRangeWide[1];
                        fMesonTrueIntUnweightedRangeNarrow[1]   = fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];

                    }
                }
            }
            
            if (fAdvancedMesonQA && (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5)){
                if (fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]->GetEntries() !=0){
                    if(fCrysFitting==0){
                        fFileErrLog << "Using exp fit"<<endl;
                        FitTrueInvMassInPtBins(fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);
//                         if( fMode == 4 ) {
//                         GausFitSubtractedInvMassInPtBinsNew(fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("GausFitFuncMCTrueBinCaloPhoton%02d",iPt),kTRUE);
//                         }
                    } else {
                        fFileErrLog << "Using Crystal Ball function"<<endl;
                        FitCBSubtractedInvMassInPtBins(fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncMCTrueBinCaloPhoton%02d",iPt),kTRUE);
                    }
                    
                    fFitTrueSignalCaloPhotonInvMassPtBin[iPt]   = fFitReco;
                    if (fFitTrueSignalCaloPhotonInvMassPtBin[iPt] != 0x00){
                        fMesonTrueMassCaloPhoton[iPt]           = fFitTrueSignalCaloPhotonInvMassPtBin[iPt]->GetParameter(1);
                        fMesonTrueMassErrorCaloPhoton[iPt]      = fFitTrueSignalCaloPhotonInvMassPtBin[iPt]->GetParError(1);
                        CalculateFWHM(fFitTrueSignalCaloPhotonInvMassPtBin[iPt]);
                        fMesonTrueFWHMCaloPhoton[iPt]           = fFWHMFunc;
                        fMesonTrueFWHMErrorCaloPhoton[iPt]      = fFWHMFuncError;
                    } else {
                        fMesonTrueMassCaloPhoton[iPt]           = 0.;
                        fMesonTrueMassErrorCaloPhoton[iPt]      = 1.;
                        fMesonTrueFWHMCaloPhoton[iPt]           = 0.;
                        fMesonTrueFWHMErrorCaloPhoton[iPt]      = 0.;
                    }
                    
                    IntegrateHistoInvMassStream( fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt], fIntFixedRange);
                    fMesonTrueYieldGammaFixedWindow[iPt]        = fYields;
                    fMesonTrueYieldGammaErrorFixedWindow[iPt]   = fYieldsError;
                        
                }
                if (fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]->GetEntries() !=0){
                    if(fCrysFitting==0){
                        fFileErrLog << "Using exp fit"<<endl;
                        FitTrueInvMassInPtBins(fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);
//                         if( fMode == 4 ){
//                         GausFitSubtractedInvMassInPtBinsNew(fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("GausFitFuncMCTrueBinConvCaloPhoton%02d",iPt),kTRUE);
//                         }
                    } else {
                        fFileErrLog << "Using Crystal Ball function"<<endl;
                        FitCBSubtractedInvMassInPtBins(fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt],
                                                    fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncMCTrueBinConvCaloPhoton%02d",iPt),kTRUE);
                    }
                    fFitTrueSignalCaloConvPhotonInvMassPtBin[iPt]   = fFitReco;
                    if (fFitTrueSignalCaloConvPhotonInvMassPtBin[iPt] != 0x00){
                        fMesonTrueMassCaloConvPhoton[iPt]           = fFitTrueSignalCaloConvPhotonInvMassPtBin[iPt]->GetParameter(1);
                        fMesonTrueMassErrorCaloConvPhoton[iPt]      = fFitTrueSignalCaloConvPhotonInvMassPtBin[iPt]->GetParError(1);
                        CalculateFWHM(fFitTrueSignalCaloConvPhotonInvMassPtBin[iPt]);
                        fMesonTrueFWHMCaloConvPhoton[iPt]           = fFWHMFunc;
                        fMesonTrueFWHMErrorCaloConvPhoton[iPt]      = fFWHMFuncError;
                    } else {
                        fMesonTrueMassCaloConvPhoton[iPt]           = 0.;
                        fMesonTrueMassErrorCaloConvPhoton[iPt]      = 1.;
                        fMesonTrueFWHMCaloConvPhoton[iPt]           = 0.;
                        fMesonTrueFWHMErrorCaloConvPhoton[iPt]      = 0.;
                    }
                    
                    IntegrateHistoInvMassStream( fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt], fIntFixedRange);
                    if (fMode == 4 || fMode ==5 ){
                        fMesonTrueYieldConvGammaConvGammaFixedWindow[iPt]       = fYields;
                        fMesonTrueYieldConvGammaConvGammaErrorFixedWindow[iPt]  = fYieldsError;
                    } else {
                        fMesonTrueYieldGammaConvGammaFixedWindow[iPt]           = fYields;
                        fMesonTrueYieldGammaConvGammaErrorFixedWindow[iPt]      = fYieldsError;
                        fMesonTrueYieldConvGammaConvGammaFixedWindow[iPt]       = 0;
                        fMesonTrueYieldConvGammaConvGammaErrorFixedWindow[iPt]  = 0;
                    }
                    
                }
                if (fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt]->GetEntries() !=0){
                    if(fCrysFitting==0){
                        fFileErrLog << "Using exp fit"<<endl;
                        FitTrueInvMassInPtBins(fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);
//                         if( fMode == 4 ) {
//                         GausFitSubtractedInvMassInPtBinsNew(fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("GausFitFuncMCTrueBinCaloMergedCluster%02d",iPt),kTRUE);
//                         }
                    } else {
                        fFileErrLog << "Using Crystal Ball function"<<endl;
                        FitCBSubtractedInvMassInPtBins(fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt],
                                                    fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncMCTrueBinCaloMergedCluster%02d",iPt),kTRUE);
                    }

                    fFitTrueSignalCaloMergedClusterInvMassPtBin[iPt]    = fFitReco;
                    if (fFitTrueSignalCaloMergedClusterInvMassPtBin[iPt] != 0x00){
                        fMesonTrueMassCaloMergedCluster[iPt]            = fFitTrueSignalCaloMergedClusterInvMassPtBin[iPt]->GetParameter(1);
                        fMesonTrueMassErrorCaloMergedCluster[iPt]       = fFitTrueSignalCaloMergedClusterInvMassPtBin[iPt]->GetParError(1);
                        CalculateFWHM(fFitTrueSignalCaloMergedClusterInvMassPtBin[iPt]);
                        fMesonTrueFWHMCaloMergedCluster[iPt]            = fFWHMFunc;
                        fMesonTrueFWHMErrorCaloMergedCluster[iPt]       = fFWHMFuncError;
                    } else {
                        fMesonTrueMassCaloMergedCluster[iPt]            = 0.;
                        fMesonTrueMassErrorCaloMergedCluster[iPt]       = 1.;
                        fMesonTrueFWHMCaloMergedCluster[iPt]            = 0.;
                        fMesonTrueFWHMErrorCaloMergedCluster[iPt]       = 0.;
                    }
                }
                if (fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt]->GetEntries() !=0){
                    if(fCrysFitting==0){
                        fFileErrLog << "Using exp fit"<<endl;
                        FitTrueInvMassInPtBins(fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);
//                         if( fMode == 4 ) {
//                         GausFitSubtractedInvMassInPtBinsNew(fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt],
//                                                     fMesonIntDeltaRange,iPt,kFALSE,Form("GausFitFuncMCTrueBinCaloMergedClusterPartConv%02d",iPt),kTRUE);
//                         }
                    } else {
                        fFileErrLog << "Using Crystal Ball function"<<endl;
                        FitCBSubtractedInvMassInPtBins(fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt],
                                                    fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncMCTrueBinCaloMergedClusterPartConv%02d",iPt),kTRUE);
                    }

                    fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[iPt]    = fFitReco;
                    if (fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[iPt] != 0x00){
                        fMesonTrueMassCaloMergedClusterPartConv[iPt]            = fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[iPt]->GetParameter(1);
                        fMesonTrueMassErrorCaloMergedClusterPartConv[iPt]       = fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[iPt]->GetParError(1);
                        CalculateFWHM(fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[iPt]);
                        fMesonTrueFWHMCaloMergedClusterPartConv[iPt]            = fFWHMFunc;
                        fMesonTrueFWHMErrorCaloMergedClusterPartConv[iPt]       = fFWHMFuncError;
                    } else {
                        fMesonTrueMassCaloMergedClusterPartConv[iPt]            = 0.;
                        fMesonTrueMassErrorCaloMergedClusterPartConv[iPt]       = 1.;
                        fMesonTrueFWHMCaloMergedClusterPartConv[iPt]            = 0.;
                        fMesonTrueFWHMErrorCaloMergedClusterPartConv[iPt]       = 0.;
                    }
                }
            }
            if (fAdvancedMesonQA && (fMode == 4 || fMode == 5)){
                if (fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[iPt]->GetEntries() !=0){
                    if(fCrysFitting==0){
                        fFileErrLog << "Using exp fit"<<endl;
                        FitTrueInvMassInPtBins(fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);
//                         if( fMode == 4) {
//                         GausFitSubtractedInvMassInPtBinsNew(fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[iPt],
//                                                     fMesonIntDeltaRange,iPt,kFALSE,Form("GausFitFuncMCTrueBinCaloMixedCaloConvPhoton%02d",iPt),kTRUE);
//                         }
                    } else {
                        fFileErrLog << "Using Crystal Ball function"<<endl;
                        FitCBSubtractedInvMassInPtBins(fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[iPt],
                                                    fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncMCTrueBinCaloMixedCaloConvPhoton%02d",iPt),kTRUE);
                    }

                    IntegrateHistoInvMassStream( fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[iPt], fIntFixedRange);
                    fMesonTrueYieldGammaConvGammaFixedWindow[iPt]           = fYields;
                    fMesonTrueYieldGammaConvGammaErrorFixedWindow[iPt]      = fYieldsError;
                    
                    fFitTrueSignalMixedCaloConvPhotonInvMassPtBin[iPt]      = fFitReco;
                    if (fFitTrueSignalMixedCaloConvPhotonInvMassPtBin[iPt] != 0x00){
                        fMesonTrueMassMixedCaloConvPhoton[iPt]              = fFitTrueSignalMixedCaloConvPhotonInvMassPtBin[iPt]->GetParameter(1);
                        fMesonTrueMassErrorMixedCaloConvPhoton[iPt]         = fFitTrueSignalMixedCaloConvPhotonInvMassPtBin[iPt]->GetParError(1);
                        CalculateFWHM(fFitTrueSignalMixedCaloConvPhotonInvMassPtBin[iPt]);
                        fMesonTrueFWHMMixedCaloConvPhoton[iPt]              = fFWHMFunc;
                        fMesonTrueFWHMErrorMixedCaloConvPhoton[iPt]         = fFWHMFuncError;
                    } else {
                        fMesonTrueMassMixedCaloConvPhoton[iPt]              = 0.;
                        fMesonTrueMassErrorMixedCaloConvPhoton[iPt]         = 1.;
                        fMesonTrueFWHMMixedCaloConvPhoton[iPt]              = 0.;
                        fMesonTrueFWHMErrorMixedCaloConvPhoton[iPt]         = 0.;
                    }
                }
            }
    
            
            IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonTrueIntRange);
            fMesonTrueYields[iPt]                       = fYields;
            fMesonTrueYieldsError[iPt]                  = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

            fFileDataLog<< endl <<"True histo wide range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonTrueIntRangeWide);
            fMesonTrueYieldsWide[iPt]                   = fYields;
            fMesonTrueYieldsWideError[iPt]              = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

            fFileDataLog<< endl <<"True histo narrow range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonTrueIntRangeNarrow);
            fMesonTrueYieldsNarrow[iPt]                 = fYields;
            fMesonTrueYieldsNarrowError[iPt]            = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

            cout<< "Analyse reweighted histos" << endl;
            
            fFileDataLog<< endl <<"True histo normal range reweighted" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtReweightedBins[iPt], fMesonTrueIntReweightedRange);
            fMesonTrueYieldsReweighted[iPt]             = fYields;
            fMesonTrueYieldsReweightedError[iPt]        = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

            fFileDataLog<< endl <<"True histo wide range reweighted" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtReweightedBins[iPt], fMesonTrueIntReweightedRangeWide);
            fMesonTrueYieldsReweightedWide[iPt]         = fYields;
            fMesonTrueYieldsReweightedWideError[iPt]    = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

            fFileDataLog<< endl <<"True histo narrow range reweighted" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtReweightedBins[iPt], fMesonTrueIntReweightedRangeNarrow);
            fMesonTrueYieldsReweightedNarrow[iPt]       = fYields;
            fMesonTrueYieldsReweightedNarrowError[iPt]  = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;


            fFileDataLog<< endl <<"True histo normal range unweighted" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt], fMesonTrueIntUnweightedRange);
            fMesonTrueYieldsUnweighted[iPt]             = fYields;
            fMesonTrueYieldsUnweightedError[iPt]        = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

            fFileDataLog<< endl <<"True histo wide range unweighted" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt], fMesonTrueIntUnweightedRangeWide);
            fMesonTrueYieldsUnweightedWide[iPt]         = fYields;
            fMesonTrueYieldsUnweightedWideError[iPt]    = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

            fFileDataLog<< endl <<"True histo narrow range unweighted" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt], fMesonTrueIntUnweightedRangeNarrow);
            fMesonTrueYieldsUnweightedNarrow[iPt]       = fYields;
            fMesonTrueYieldsUnweightedNarrowError[iPt]  = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
            
            if (fEnableDCMeson){
        
                IntegrateHistoInvMassStream( fHistoMappingTrueMesonDCInvMassPtBins[iPt], fMesonTrueIntRange);
                fMesonTrueYieldsDC[iPt]                         = fYields;
                fMesonTrueYieldsDCError[iPt]                    = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
                cout << "***********************************************" << endl;
                cout << "***********************************************" << endl;
                cout << "***********************************************" << endl;
        
                cout << "***********************************************" << endl;
            }
            
            if (meson.Contains("Pi0")){
                fFileDataLog<< endl <<"TrueSecTotal histo" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingTrueSecMesonInvMassPtBins[iPt], fMesonTrueIntRange);
                fMesonTrueSecYields[iPt]                    = fYields;
                fMesonTrueSecYieldsError[iPt]               = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
                    
                fFileDataLog<< endl <<"TrueSecTotal histo wide range  " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingTrueSecMesonInvMassPtBins[iPt], fMesonTrueIntRangeWide);
                fMesonTrueSecYieldsWide[iPt]                = fYields;
                fMesonTrueSecYieldsWideError[iPt]           = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
                    
                fFileDataLog<< endl <<"TrueSecTotal histo narrow range " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingTrueSecMesonInvMassPtBins[iPt], fMesonTrueIntRangeNarrow);
                fMesonTrueSecYieldsNarrow[iPt]              = fYields;
                fMesonTrueSecYieldsNarrowError[iPt]         = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

                fFileDataLog<< endl <<"TrueSecK0s histo " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt], fMesonTrueIntRange);
                fMesonTrueSecFromK0SYields[iPt]             = fYields;
                fMesonTrueSecFromK0SYieldsError[iPt]        = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
            
                fFileDataLog<< endl <<"TrueSecK0s histo wide range " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt], fMesonTrueIntRangeWide);
                fMesonTrueSecFromK0SYieldsWide[iPt]         = fYields;
                fMesonTrueSecFromK0SYieldsWideError[iPt]    = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
            
                fFileDataLog<< endl <<"TrueSecK0s histo narrow range " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt], fMesonTrueIntRangeNarrow);
                fMesonTrueSecFromK0SYieldsNarrow[iPt]       = fYields;
                fMesonTrueSecFromK0SYieldsNarrowError[iPt]  = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
                
                fFileDataLog<< endl <<"TrueSecLambda histo " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingTrueSecFromLambdaMesonInvMassPtBins[iPt], fMesonTrueIntRange);
                fMesonTrueSecFromLambdaYields[iPt]          = fYields;
                fMesonTrueSecFromLambdaYieldsError[iPt]     = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
            
                fFileDataLog<< endl <<"TrueSecLambda histo wide range " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingTrueSecFromLambdaMesonInvMassPtBins[iPt], fMesonTrueIntRangeWide);
                fMesonTrueSecFromLambdaYieldsWide[iPt]      = fYields;
                fMesonTrueSecFromLambdaYieldsWideError[iPt] = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
            
                fFileDataLog<< endl <<"TrueSecLambda histo narrow range " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingTrueSecFromLambdaMesonInvMassPtBins[iPt], fMesonTrueIntRangeNarrow);
                fMesonTrueSecFromLambdaYieldsNarrow[iPt]        = fYields;
                fMesonTrueSecFromLambdaYieldsNarrowError[iPt]   = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
            } else {
                fMesonTrueSecFromK0SYields[iPt]                     = 0;
                fMesonTrueSecFromK0SYieldsError[iPt]                = 0;
                fMesonTrueSecFromK0SYieldsWide[iPt]                 = 0;
                fMesonTrueSecFromK0SYieldsWideError[iPt]            = 0;
                fMesonTrueSecFromK0SYieldsNarrow[iPt]               = 0;
                fMesonTrueSecFromK0SYieldsNarrowError[iPt]          = 0;
                fMesonTrueSecFromLambdaYields[iPt]                  = 0;
                fMesonTrueSecFromLambdaYieldsError[iPt]             = 0;
                fMesonTrueSecFromLambdaYieldsWide[iPt]              = 0;
                fMesonTrueSecFromLambdaYieldsWideError[iPt]         = 0;
                fMesonTrueSecFromLambdaYieldsNarrow[iPt]            = 0;
                fMesonTrueSecFromLambdaYieldsNarrowError[iPt]       = 0;
            }


            if( (fGGYields[iPt] - fMesonTrueYields[iPt]) > 0) {
                fMesonTrueSB[iPt]               = fMesonTrueYields[iPt] / ( fGGYields[iPt] - fMesonTrueYields[iPt] );
                fMesonTrueSign[iPt]             = fMesonTrueYields[iPt] / pow( ( fGGYields[iPt] - fMesonTrueYields[iPt] ) , 0.5);
                fMesonTrueSBError[iPt]          = 0;
                fMesonTrueSignError[iPt]        = 0;
            }
            else {
                fMesonTrueSB[iPt]               = 0.;
                fMesonTrueSign[iPt]             = 0.;
                fMesonTrueSBError[iPt]          = 0.;
                fMesonTrueSignError[iPt]        = 0.;
            }
        }
        fFileDataLog<< "Residual Background leftover norm integration/right norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFunc[iPt] << "\t +- \t" << fMesonYieldsResidualBckFuncError[iPt] << endl<< endl;
        
        /////////////////////// added to check yields //////////////////////////////////////////////////////////   
        fTotalBckYields[iPt] = fBckYields[iPt] + fMesonYieldsResidualBckFunc[iPt];
        fTotalBckYieldsError[iPt] = pow(fBckYieldsError[iPt]*fBckYieldsError[iPt] + fMesonYieldsResidualBckFuncError[iPt]*fMesonYieldsResidualBckFuncError[iPt],0.5);
        fFileDataLog<< "Total Background norm integration/right norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fTotalBckYields[iPt] << "\t +- \t" << fTotalBckYieldsError[iPt] << endl<< endl;
        fFileDataLog<< "Background norm integration/right norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fBckYields[iPt] << "\t +- \t" << fBckYieldsError[iPt] << endl<< endl;
        ///////////////////////////////////////////////////////////////////////////////////////////////////////
                
        fMesonYieldsCorResidualBckFunc[iPt]         = fMesonYields[iPt]- fMesonYieldsResidualBckFunc[iPt];
        fMesonYieldsCorResidualBckFuncError[iPt]    = pow(( fMesonYieldsError[iPt]*fMesonYieldsError[iPt] +
                                                            fMesonYieldsResidualBckFuncError[iPt]*fMesonYieldsResidualBckFuncError[iPt]),0.5);
        fMesonYieldsPerEvent[iPt]                   = fMesonYieldsCorResidualBckFunc[iPt]/fNEvents;
        fMesonYieldsPerEventError[iPt]              = fMesonYieldsCorResidualBckFuncError[iPt]/fNEvents;
        
        //Integrate Fit Function
        IntegrateFitFunc( fFitSignalPeakPosInvMassPtBin[iPt], fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRange);
        fMesonYieldsFunc[iPt]                       = fYieldsFunc;

        //SB default
        fMesonSBdefault[iPt]                        = fMesonYieldsCorResidualBckFunc[iPt]/fTotalBckYields[iPt];
        fMesonSBdefaultError[iPt]                   = pow( pow(fMesonYieldsCorResidualBckFuncError[iPt]/fTotalBckYields[iPt], 2.) +
                                                           pow((fTotalBckYieldsError[iPt]*fMesonYieldsCorResidualBckFunc[iPt])/(fTotalBckYields[iPt] *fTotalBckYields[iPt]), 2.), 0.5);
        
        //Significance default
        fMesonSigndefault[iPt]                      = fMesonYieldsCorResidualBckFunc[iPt]/pow(fMesonYieldsCorResidualBckFunc[iPt] + fTotalBckYields[iPt],0.5);
        Double_t a                                  = ( pow(fMesonYieldsCorResidualBckFunc[iPt] + fTotalBckYields[iPt], -0.5) -
                                                        0.5*fMesonYieldsCorResidualBckFunc[iPt]*pow(fMesonYieldsCorResidualBckFunc[iPt] +
                                                                                                    fTotalBckYields[iPt], -1.5) * fMesonYieldsCorResidualBckFuncError[iPt]);
        Double_t b                                  = 0.5*fMesonYieldsCorResidualBckFunc[iPt]*pow(fMesonYieldsCorResidualBckFunc[iPt]
                                                                                                  + fTotalBckYields[iPt],-1.5) * fTotalBckYieldsError[iPt];
        fMesonSigndefaultError[iPt]                 = pow( a*a + b*b, 0.5);

        
        if( fFitBckInvMassPtBin[iPt]->Integral(fMesonMass[iPt]-fMesonFWHM[iPt], fMesonMass[iPt]+fFitSignalInvMassPtBin[iPt]->GetParameter(2))!=0){
            Double_t background         = fFitBckInvMassPtBin[iPt]->Integral(fMesonMass[iPt]-fMesonFWHM[iPt], 
                                                                    fMesonMass[iPt]+fFitSignalInvMassPtBin[iPt]->GetParameter(2));
            Double_t backgroundErr      = fFitBckInvMassPtBin[iPt]->IntegralError(fMesonMass[iPt]-fMesonFWHM[iPt],
                                                                    fMesonMass[iPt]+fFitSignalInvMassPtBin[iPt]->GetParameter(2),
                                                                    fFitBckInvMassPtBin[iPt]->GetParameters(), NULL /* cov matrix */);
            Double_t signal             = fFitSignalInvMassPtBin[iPt]->Integral(fMesonMass[iPt]-fMesonFWHM[iPt], 
                                                                       fMesonMass[iPt]+fFitSignalInvMassPtBin[iPt]->GetParameter(2)) - background;
            Double_t signalErr          = pow( pow(fFitSignalInvMassPtBin[iPt]->IntegralError(fMesonMass[iPt]-fMesonFWHM[iPt],
                                                                                fMesonMass[iPt]+fFitSignalInvMassPtBin[iPt]->GetParameter(2),
                                                                                fFitSignalInvMassPtBin[iPt]->GetParameters(),
                                                                                NULL /* cov matrix */), 2) + pow(backgroundErr, 2), 0.5);
            fMesonSB[iPt]               = signal/ background;
            fMesonSBError[iPt]          = pow( pow(signalErr/background, 2.) + pow(signal*backgroundErr/(background*background), 2.) ,0.5);
            fMesonSign[iPt]             = signal/ pow(background + signal,0.5);
            fMesonSignError[iPt]        = pow(pow( (pow(background + signal, -0.5) - 0.5*signal*pow(background+signal, -1.5)) * signalErr, 2) + pow( 0.5*signal*backgroundErr*pow(signal+background, -1.5), 2), 0.5);
        }else{
            fMesonSB[iPt]               = 0.;
            fMesonSBError[iPt]          = 0.;
            fMesonSBdefault[iPt]        = 0.;
            fMesonSBdefaultError[iPt]   = 0.;
            fMesonSign[iPt]             = 0.;
            fMesonSignError[iPt]        = 0.;
            fMesonSigndefault[iPt]      = 0.;
            fMesonSigndefaultError[iPt] = 0.;
        }

        //       cout<< "iPt"<< iPt<< " "<< "FWHM done"<<endl;

        // Wide integration mass window
        //       cout<< "iPt"<< iPt<< " "<< "wide range"<<endl;
        fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "wide range/right normalization" << endl;

        if(fCrysFitting==0){
            fFileErrLog << "Using exp fit"<<endl;
            FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRangeWide,iPt,kFALSE);
            fMesonYieldsResidualBckFuncWide[iPt]        = fIntLinearBck;
            fMesonYieldsResidualBckFuncWideError[iPt]   = fIntLinearBckError;
        } else {
//           fFileErrLog << "Using Crystal Ball function"<<endl;
//           FitCBSubtractedInvMassInPtBins(fHistoMappingSignalRemainingBGSubInvMassPtBin[iPt], fMesonIntDeltaRangeWide,iPt,kFALSE,Form("CBFitFuncNormalWideBin%02d",iPt),kFALSE);
            fMesonYieldsResidualBckFuncWide[iPt]        = 0;
            fMesonYieldsResidualBckFuncWideError[iPt]   = 0;

        }

        //FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt],fMesonIntDeltaRangeWide,iPt,kFALSE);
        fFileDataLog<< "Residual Background leftover wide integration/right norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncWide[iPt] <<"\t +- \t" << fMesonYieldsResidualBckFuncWideError[iPt] <<endl<< endl;
        fMesonYieldsCorResidualBckFuncWide[iPt]         = fMesonYieldsWide[iPt]- fMesonYieldsResidualBckFuncWide[iPt];

        fTotalBckYieldsWide[iPt]                        = fBckYieldsWide[iPt] + fMesonYieldsResidualBckFuncWide[iPt];
        fTotalBckYieldsWideError[iPt]                   = pow( fBckYieldsWideError[iPt]*fBckYieldsWideError[iPt] + 
                                                               fMesonYieldsResidualBckFuncWideError[iPt]*fMesonYieldsResidualBckFuncWideError[iPt],0.5);

        fMesonYieldsCorResidualBckFuncWideError[iPt]    = pow((fMesonYieldsWideError[iPt]*fMesonYieldsWideError[iPt]+
                                                               fMesonYieldsResidualBckFuncWideError[iPt]*fMesonYieldsResidualBckFuncWideError[iPt]),0.5);
        fMesonYieldsPerEventWide[iPt]                   = fMesonYieldsCorResidualBckFuncWide[iPt]/fNEvents;
        fMesonYieldsPerEventWideError[iPt]              = fMesonYieldsCorResidualBckFuncWideError[iPt]/fNEvents;

        if( fTotalBckYieldsWide[iPt]!=0){
            fMesonSBWide[iPt]           = fMesonYieldsCorResidualBckFuncWide[iPt]/fTotalBckYieldsWide[iPt];
            fMesonSBWideError[iPt]      = pow(pow(fMesonYieldsCorResidualBckFuncWideError[iPt]/fTotalBckYieldsWide[iPt],2.)+pow(fMesonYieldsCorResidualBckFuncWide[iPt]/(fTotalBckYieldsWide[iPt]*fTotalBckYieldsWide[iPt])*fTotalBckYieldsWideError[iPt],2.) ,0.5);
            fMesonSignWide[iPt]         = fMesonYieldsCorResidualBckFuncWide[iPt]/pow(fTotalBckYieldsWide[iPt],0.5);
            fMesonSignWideError[iPt]    = pow(pow(fMesonYieldsCorResidualBckFuncWideError[iPt]/pow(fTotalBckYieldsWide[iPt],0.5),2.)+pow(0.5*fMesonYieldsCorResidualBckFuncWide[iPt]/pow(fTotalBckYieldsWide[iPt],1.5)*fTotalBckYieldsWideError[iPt],2.) ,0.5);

        }else{
            fMesonSBWide[iPt]           = 0.;
            fMesonSBWideError[iPt]      = 0.;
            fMesonSignWide[iPt]         = 0.;
            fMesonSignWideError[iPt]    = 0.;
        }


        // Narrow integration mass window
        fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "narrow range/right normalization" << endl; ;
        //       cout<< "iPt"<< iPt<< " "<< "narrow range"<<endl;

        if(fCrysFitting==0){
            fFileErrLog << "Using exp fit"<<endl;
            FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRangeNarrow,iPt,kFALSE);
            fMesonYieldsResidualBckFuncNarrow[iPt]          = fIntLinearBck;
            fMesonYieldsResidualBckFuncNarrowError[iPt]     = fIntLinearBckError;

        } else {
//          fFileErrLog << "Using Crystal Ball function"<<endl;
//          FitCBSubtractedInvMassInPtBins(fHistoMappingSignalRemainingBGSubInvMassPtBin[iPt], fMesonIntDeltaRangeNarrow,iPt,kFALSE, Form("CBFitFuncNormalNarrowBin%02d",iPt),kFALSE);
            fMesonYieldsResidualBckFuncNarrow[iPt]          = 0;
            fMesonYieldsResidualBckFuncNarrowError[iPt]     = 0;

        }

        //FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt],fMesonIntDeltaRangeNarrow,iPt,kFALSE);
        fFileDataLog<< "Residual Background leftover narrow integration/right norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncNarrow[iPt] <<"\t +- \t" << fMesonYieldsResidualBckFuncNarrowError[iPt]<<endl<< endl;
        fMesonYieldsCorResidualBckFuncNarrow[iPt]           = fMesonYieldsNarrow[iPt]- fMesonYieldsResidualBckFuncNarrow[iPt];
        fMesonYieldsCorResidualBckFuncNarrowError[iPt]      = pow(( fMesonYieldsNarrowError[iPt]*fMesonYieldsNarrowError[iPt]+
                                                                    fMesonYieldsResidualBckFuncNarrowError[iPt]*fMesonYieldsResidualBckFuncNarrowError[iPt]),0.5);
        fMesonYieldsPerEventNarrow[iPt]                     = fMesonYieldsCorResidualBckFuncNarrow[iPt]/fNEvents;
        fMesonYieldsPerEventNarrowError[iPt]                = fMesonYieldsCorResidualBckFuncNarrowError[iPt]/fNEvents;

        fTotalBckYieldsNarrow[iPt]                          = fBckYieldsNarrow[iPt] + fMesonYieldsResidualBckFuncNarrow[iPt];
        fTotalBckYieldsNarrowError[iPt]                     = pow(  fBckYieldsNarrowError[iPt]*fBckYieldsNarrowError[iPt] + 
                                                                    fMesonYieldsResidualBckFuncNarrowError[iPt]*fMesonYieldsResidualBckFuncNarrowError[iPt],0.5);

        //SB default Narrow
        fMesonSBdefaultNarrow[iPt]                          = fMesonYieldsCorResidualBckFuncNarrow[iPt]/fTotalBckYieldsNarrow[iPt];
        fMesonSBdefaultNarrowError[iPt]                     = pow( pow(fMesonYieldsCorResidualBckFuncNarrowError[iPt]/fTotalBckYieldsNarrow[iPt],2.) + 
                                                                   pow((fTotalBckYieldsNarrowError[iPt]*fMesonYieldsCorResidualBckFuncNarrow[iPt])/
                                                                   (fTotalBckYieldsNarrow[iPt] *fTotalBckYieldsNarrow[iPt] ),2.) ,0.5); 
        
        //Significance default Narrow
        fMesonSigndefaultNarrow[iPt]                        = fMesonYieldsCorResidualBckFuncNarrow[iPt]/pow(fMesonYieldsCorResidualBckFuncNarrow[iPt] + fTotalBckYieldsNarrow[iPt],0.5);
        Double_t aNarrow                                    = pow( ( pow(fMesonYieldsCorResidualBckFuncNarrow[iPt] + fTotalBckYieldsNarrow[iPt],-0.5) - 
                                                                     0.5*fMesonYieldsCorResidualBckFuncNarrow[iPt]*pow(fMesonYieldsCorResidualBckFuncNarrow[iPt] + 
                                                                     fTotalBckYieldsNarrow[iPt],-1.5)/(fMesonYieldsCorResidualBckFuncNarrow[iPt] + fTotalBckYieldsNarrow[iPt]) )
                                                                   *fMesonYieldsCorResidualBckFuncNarrowError[iPt] ,2);
        Double_t bNarrow                                    = pow( (fMesonSigndefaultNarrow[iPt])*(0.5*fTotalBckYieldsNarrowError[iPt]* 
                                                                    pow(fMesonYieldsCorResidualBckFuncNarrow[iPt] + fTotalBckYieldsNarrow[iPt],-1.5))  ,2.);
        fMesonSigndefaultNarrowError[iPt]                   = pow( aNarrow + bNarrow ,0.5); 

        
        if( fTotalBckYieldsNarrow[iPt]!=0){
            fMesonSBNarrow[iPt]                 = fMesonYieldsCorResidualBckFuncNarrow[iPt]/fTotalBckYieldsNarrow[iPt];
            fMesonSBNarrowError[iPt]            = pow(pow(fMesonYieldsCorResidualBckFuncNarrowError[iPt]/fTotalBckYieldsNarrow[iPt],2.) + 
                                                      pow(fMesonYieldsCorResidualBckFuncNarrow[iPt]/(fTotalBckYieldsNarrow[iPt]*fTotalBckYieldsNarrow[iPt])*fTotalBckYieldsNarrowError[iPt],2.) ,0.5);
            fMesonSignNarrow[iPt]               = fMesonYieldsCorResidualBckFuncNarrow[iPt]/pow(fTotalBckYieldsNarrow[iPt],0.5);
            fMesonSignNarrowError[iPt]          = pow(pow(fMesonYieldsCorResidualBckFuncNarrowError[iPt]/pow(fTotalBckYieldsNarrow[iPt],0.5),2.) + 
                                                      pow(0.5*fMesonYieldsCorResidualBckFuncNarrow[iPt]/pow(fTotalBckYieldsNarrow[iPt],1.5)*fTotalBckYieldsNarrowError[iPt],2.) ,0.5);
        } else {
            fMesonSBNarrow[iPt]                 = 0.;
            fMesonSBNarrowError[iPt]            = 0.;
            fMesonSignNarrow[iPt]               = 0.;
            fMesonSignNarrowError[iPt]          = 0.;
            fMesonSigndefaultNarrow[iPt]        = 0.;
            fMesonSigndefaultNarrowError[iPt]   = 0.;
            fMesonSBdefaultNarrow[iPt]          = 0.;
            fMesonSBdefaultNarrowError[iPt]     = 0.;
        }

        //////////////////////////////// Start Analysis with  Normalization at the left of the Meson Peak
        // Function to subtract GG minus Bck
        cout << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
        ProcessEM( fHistoMappingGGInvMassPtBin[iPt], fHistoMappingBackInvMassPtBin[iPt], fBGFitRangeLeft);
        fHistoMappingSignalInvMassLeftPtBin[iPt]        = fSignal;
        fHistoMappingBackNormInvMassLeftPtBin[iPt]      = fBckNorm;


        fHistoMappingSignalInvMassLeftPtBin[iPt]->SetName(Form("fHistoMappingSignalLeftInvMass_in_Pt_Bin%02d",iPt));
        //       cout<< "iPt"<< iPt<< " "<< "standard range"<<endl;
        // Fitting the subtracted spectra
        fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "normal range/left normalization" << endl;

        fFitInvMassLeftPtBin[iPt] =0x00;
        if(fCrysFitting==0){
//         if( fMode == 4 ) {
//             fFileErrLog << "Using exp fit"<<endl;
//             FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
//             fFitInvMassLeftPtBin2[iPt]                              = fFitReco;
//             fFitSignalPeakPosInvMassLeftPtBin2[iPt]                 = fFitGausExp;
//             fFitBckInvMassLeftPtBin2[iPt]                           = fFitLinearBck;
//             fMesonYieldsResidualBckFuncLeft[iPt]                    = fIntLinearBck;
//             fMesonYieldsResidualBckFuncLeftError[iPt]               = fIntLinearBckError;
//             GausFitSubtractedInvMassInPtBinsNew(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE, Form("GausFitFuncLeftBin%02d",iPt),kFALSE);
//             fHistoMappingSignalRemainingBGSubInvMassLeftPtBin[iPt]  = (TH1D*)fCopySignal->Clone(Form("histoSignalRemainingBGSubtractedLeftBin%02d",iPt));
//             fHistoMappingRemainingBGInvMassLeftPtBin[iPt]           = (TH1D*)fCopyOnlyBG->Clone(Form("histoRemainingBGLeftBin%02d",iPt));
//             fFitRemainingBGInvMassLeftPtBin[iPt]                    = fFitLinearBck;
//             fFitSignalPeakPosInvMassLeftPtBin[iPt]                  = fFitGausExp;
//             fFitInvMassLeftPtBin[iPt]                               = fFitReco;
//             fFitBckInvMassLeftPtBin[iPt]                            = fFitLinearBck;
//         } else {
            fFileErrLog << "Using exp fit"<<endl;
            FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
            fFitInvMassLeftPtBin[iPt]                               = fFitReco;
            fFitSignalPeakPosInvMassLeftPtBin[iPt]                  = fFitGausExp;
            fFitBckInvMassLeftPtBin[iPt]                            = fFitLinearBck;
            fMesonYieldsResidualBckFuncLeft[iPt]                    = fIntLinearBck;
            fMesonYieldsResidualBckFuncLeftError[iPt]               = fIntLinearBckError;
//         }
        } else {
            fFileErrLog << "Using Crystal Ball function"<<endl;
            FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE, Form("CBFitFuncLeftBin%02d",iPt),kFALSE);
            fHistoMappingSignalRemainingBGSubInvMassLeftPtBin[iPt]  = (TH1D*)fCopySignal->Clone(Form("histoSignalRemainingBGSubtractedLeftBin%02d",iPt));
            fHistoMappingRemainingBGInvMassLeftPtBin[iPt]           = (TH1D*)fCopyOnlyBG->Clone(Form("histoRemainingBGLeftBin%02d",iPt));
            fFitRemainingBGInvMassLeftPtBin[iPt]                    = fFitLinearBck;
            fFitSignalPeakPosInvMassLeftPtBin[iPt]                  = fFitGausExp;
            fFitInvMassLeftPtBin[iPt]                               = fFitReco;
            fFitSignalPeakPosInvMassLeftPtBin[iPt]                  = fFitGausExp;
            fFitBckInvMassLeftPtBin[iPt]                            = fFitLinearBck;
            fMesonYieldsResidualBckFuncLeft[iPt]                    = 0;
            fMesonYieldsResidualBckFuncLeftError[iPt]               = 0;

        }
        //FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
        CalculateFWHM(fFitInvMassLeftPtBin[iPt]);
        fMesonFWHMLeft[iPt]         = fFWHMFunc;
        fMesonFWHMLeftError[iPt]    = fFWHMFuncError;
        
        if (fFitInvMassLeftPtBin[iPt] !=0x00){
            fMesonMassLeft[iPt]             = fFitInvMassLeftPtBin[iPt]->GetParameter(1);
            fMesonMassLeftError[iPt]        = fFitInvMassLeftPtBin[iPt]->GetParError(1);
            fMesonCurLeftIntRange[0]        = fMesonMassLeft[iPt] + fMesonIntDeltaRange[0];
            fMesonCurLeftIntRangeWide[0]    = fMesonMassLeft[iPt] + fMesonIntDeltaRangeWide[0];
            fMesonCurLeftIntRangeNarrow[0]  = fMesonMassLeft[iPt] + fMesonIntDeltaRangeNarrow[0];
            fMesonCurLeftIntRange[1]        = fMesonMassLeft[iPt] + fMesonIntDeltaRange[1];
            fMesonCurLeftIntRangeWide[1]    = fMesonMassLeft[iPt] + fMesonIntDeltaRangeWide[1];
            fMesonCurLeftIntRangeNarrow[1]  = fMesonMassLeft[iPt] + fMesonIntDeltaRangeNarrow[1];
        } else {
            fMesonMassLeft[iPt]             = 0.;
            fMesonMassLeftError[iPt]        = 0.;
            fMesonCurLeftIntRange[0]        = fMesonMassExpect + fMesonIntDeltaRange[0];
            fMesonCurLeftIntRangeWide[0]    = fMesonMassExpect + fMesonIntDeltaRangeWide[0];
            fMesonCurLeftIntRangeNarrow[0]  = fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
            fMesonCurLeftIntRange[1]        = fMesonMassExpect + fMesonIntDeltaRange[1];
            fMesonCurLeftIntRangeWide[1]    = fMesonMassExpect + fMesonIntDeltaRangeWide[1];
            fMesonCurLeftIntRangeNarrow[1]  = fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];
        }

        // Integrate the bck histo
        if (fCrysFitting ==0){
//         if( fMode == 4 ){
//             fHistoMappingBackNormInvMassLeftPtBin[iPt]->Sumw2();//added
//             fHistoMappingBackNormInvMassLeftPtBin[iPt]->Add(fHistoMappingRemainingBGInvMassLeftPtBin[iPt],1);//added
//         }
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurLeftIntRange);
            fBckYieldsLeft[iPt]                 = fYields;
            fBckYieldsLeftError[iPt]            = fYieldsError;
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeWide);
            fBckYieldsLeftWide[iPt]             = fYields;
            fBckYieldsLeftWideError[iPt]        = fYieldsError;
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeNarrow);
            fBckYieldsLeftNarrow[iPt]           = fYields;
            fBckYieldsLeftNarrowError[iPt]      = fYieldsError;

            // Integrate the signal histo
            fFileDataLog<< endl <<"Signal histo normal range, left norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRange);
            fMesonYieldsLeft[iPt]               = fYields;
            fMesonYieldsLeftError[iPt]          = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
            fFileDataLog<< endl <<"Signal histo wide range, left norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeWide);
            fMesonYieldsLeftWide[iPt]           = fYields;
            fMesonYieldsLeftWideError[iPt]      = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
            fFileDataLog<< endl <<"Signal histo narrow range, left norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeNarrow);
            fMesonYieldsLeftNarrow[iPt]         = fYields;
            fMesonYieldsLeftNarrowError[iPt]    = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
        } else {
            fHistoMappingBackNormInvMassLeftPtBin[iPt]->Sumw2();
            fHistoMappingBackNormInvMassLeftPtBin[iPt]->Add(fHistoMappingRemainingBGInvMassLeftPtBin[iPt],1);
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurLeftIntRange);
            fBckYieldsLeft[iPt]                 = fYields;
            fBckYieldsLeftError[iPt]            = fYieldsError;
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeWide);
            fBckYieldsLeftWide[iPt]             = fYields;
            fBckYieldsLeftWideError[iPt]        = fYieldsError;
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeNarrow);
            fBckYieldsLeftNarrow[iPt]           = fYields;
            fBckYieldsLeftNarrowError[iPt]      = fYieldsError;

            // Integrate the signal histo
            fFileDataLog<< endl <<"Signal histo normal range, left norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalRemainingBGSubInvMassLeftPtBin[iPt], fMesonCurLeftIntRange);
            fMesonYieldsLeft[iPt]               = fYields;
            fMesonYieldsLeftError[iPt]          = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
            fFileDataLog<< endl <<"Signal histo wide range, left norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalRemainingBGSubInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeWide);
            fMesonYieldsLeftWide[iPt]           = fYields;
            fMesonYieldsLeftWideError[iPt]      = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
            fFileDataLog<< endl <<"Signal histo narrow range, left norm " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalRemainingBGSubInvMassLeftPtBin[iPt], fMesonCurLeftIntRangeNarrow);
            fMesonYieldsLeftNarrow[iPt]         = fYields;
            fMesonYieldsLeftNarrowError[iPt]    = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
        }
        fFileDataLog<< "Residual Background leftover norm integration/left norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncLeft[iPt] << "\t +- \t" <<
        fMesonYieldsResidualBckFuncLeftError[iPt]<<endl<< endl;
        fMesonYieldsCorResidualBckFuncLeft[iPt]         = fMesonYieldsLeft[iPt]- fMesonYieldsResidualBckFuncLeft[iPt];
        fMesonYieldsCorResidualBckFuncLeftError[iPt]    = pow( ( fMesonYieldsLeftError[iPt] * fMesonYieldsLeftError[iPt] + fMesonYieldsResidualBckFuncLeftError[iPt] *
                                                                 fMesonYieldsResidualBckFuncLeftError[iPt]), 0.5);
        fMesonYieldsLeftPerEvent[iPt]                   = fMesonYieldsCorResidualBckFuncLeft[iPt]/fNEvents;
        fMesonYieldsLeftPerEventError[iPt]              = fMesonYieldsCorResidualBckFuncLeftError[iPt]/fNEvents;

        fTotalBckYieldsLeft[iPt]                        = fBckYieldsLeft[iPt] + fMesonYieldsResidualBckFuncLeft[iPt];
        fTotalBckYieldsLeftError[iPt]                   = pow( fBckYieldsLeftError[iPt] * fBckYieldsLeftError[iPt] + fMesonYieldsResidualBckFuncLeftError[iPt] * 
                                                               fMesonYieldsResidualBckFuncLeftError[iPt], 0.5);

        //Integrate Fit Function
        IntegrateFitFunc( fFitSignalPeakPosInvMassLeftPtBin[iPt], fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurLeftIntRange);
        fMesonYieldsFuncLeft[iPt]                       = fYieldsFunc;

        //GetFWHM
        
        if( fFitBckInvMassLeftPtBin[iPt]->Integral(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2))!=0){
            Double_t background         = fFitBckInvMassLeftPtBin[iPt]->Integral( fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2));
            Double_t backgroundErr      = fFitBckInvMassLeftPtBin[iPt]->IntegralError(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2));
            Double_t signal             = fFitInvMassLeftPtBin[iPt]->Integral(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2)) - background;
            Double_t signalErr          = pow( pow(fFitInvMassLeftPtBin[iPt]->IntegralError(fMesonMassLeft[iPt]-fMesonFWHMLeft[iPt], 
                                                                                            fMesonMassLeft[iPt]+fFitInvMassLeftPtBin[iPt]->GetParameter(2)),2 )+ pow(backgroundErr,2),0.5);
            fMesonSBLeft[iPt]           = signal/ background;
            fMesonSBLeftError[iPt]      = pow( pow(signalErr/background,2.)+pow(signal/(background *background )*backgroundErr ,2.) ,0.5); 
            fMesonSignLeft[iPt]         = signal/ pow(background + signal,0.5);
            fMesonSignLeftError[iPt]    = pow(pow( (pow(background + signal ,0.5) - 0.5*signal*pow(background+signal,-0.5))/(background+signal) * signalErr ,2) 
                                                    + pow( 0.5*pow(signal+background,-1.5),2) ,0.5); 
        }else{
            fMesonSBLeft[iPt]           = 0.;
            fMesonSBLeftError[iPt]      = 0.;
            fMesonSignLeft[iPt]         = 0.;
            fMesonSignLeftError[iPt]    = 0.;
        }

        // Wide integration mass window
        //       cout<< "iPt"<< iPt<< " "<< "wide range"<<endl;
        fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "wide range/left normalization" << endl;

        if(fCrysFitting==0){
            fFileErrLog << "Using exp fit"<<endl;
            FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRangeWide,iPt,kFALSE);
            fMesonYieldsResidualBckFuncLeftWide[iPt]        = fIntLinearBck;
            fMesonYieldsResidualBckFuncLeftWideError[iPt]   = fIntLinearBckError;

        } else {
//           fFileErrLog << "Using Crystal Ball function"<<endl;
//           FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRangeWide,iPt,kFALSE, Form("CBFitFuncLeftWideBin%02d",iPt),kFALSE);
            fMesonYieldsResidualBckFuncLeftWide[iPt]       = 0;
            fMesonYieldsResidualBckFuncLeftWideError[iPt]    = 0;

        }

        //FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt],fMesonIntDeltaRangeWide,iPt,kFALSE);
        fFileDataLog<< "Residual Background leftover wide integration/left norm in iPt " << fBinsPt[iPt] << "-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncLeftWide[iPt] << "\t +- \t" 
        << fMesonYieldsResidualBckFuncLeftWideError[iPt] <<endl<< endl;
        fMesonYieldsCorResidualBckFuncLeftWide[iPt]         = fMesonYieldsLeftWide[iPt]- fMesonYieldsResidualBckFuncLeftWide[iPt];
        fMesonYieldsCorResidualBckFuncLeftWideError[iPt]    = pow( ( fMesonYieldsLeftWideError[iPt] * fMesonYieldsLeftWideError[iPt] + fMesonYieldsResidualBckFuncLeftWideError[iPt] *
                                                                     fMesonYieldsResidualBckFuncLeftWideError[iPt]), 0.5);
        fMesonYieldsLeftPerEventWide[iPt]                   = fMesonYieldsCorResidualBckFuncLeftWide[iPt]/fNEvents;
        fMesonYieldsLeftPerEventWideError[iPt]              = fMesonYieldsCorResidualBckFuncLeftWideError[iPt]/fNEvents;

        fTotalBckYieldsLeftWide[iPt]                        = fBckYieldsLeftWide[iPt] + fMesonYieldsResidualBckFuncLeftWide[iPt];
        fTotalBckYieldsLeftWideError[iPt]                   = pow( fBckYieldsLeftWideError[iPt]*fBckYieldsLeftWideError[iPt] + fMesonYieldsResidualBckFuncLeftWideError[iPt] *
                                                                   fMesonYieldsResidualBckFuncLeftWideError[iPt], 0.5);


        if( fTotalBckYieldsLeftWide[iPt]!=0){
            fMesonSBLeftWide[iPt]               = fMesonYieldsCorResidualBckFuncLeftWide[iPt]/fTotalBckYieldsLeftWide[iPt];
            fMesonSBLeftWideError[iPt]          = pow( pow( fMesonYieldsCorResidualBckFuncLeftWideError[iPt] / fTotalBckYieldsLeftWide[iPt], 2.) + pow( fMesonYieldsCorResidualBckFuncLeftWide[iPt] /
                                                    ( fTotalBckYieldsLeftWide[iPt] * fTotalBckYieldsLeftWide[iPt]) * fTotalBckYieldsLeftWideError[iPt], 2.) ,0.5);
            fMesonSignLeftWide[iPt]             = fMesonYieldsCorResidualBckFuncLeftWide[iPt]/pow(fTotalBckYieldsLeftWide[iPt],0.5);
            fMesonSignLeftWideError[iPt]        = pow( pow( fMesonYieldsCorResidualBckFuncLeftWideError[iPt] / pow( fTotalBckYieldsLeftWide[iPt], 0.5), 2.) + pow( 0.5 *
                                                    fMesonYieldsCorResidualBckFuncLeftWide[iPt] / pow( fTotalBckYieldsLeftWide[iPt], 1.5) * fTotalBckYieldsLeftWideError[iPt], 2.) ,0.5);
        }else{
            fMesonSBLeftWide[iPt]               = 0.;
            fMesonSBLeftWideError[iPt]          = 0.;
            fMesonSignLeftWide[iPt]             = 0.;
            fMesonSignLeftWideError[iPt]        = 0.;
        }


        // Narrow integration mass window
        //       cout<< "iPt"<< iPt<< " "<< "narrow range"<<endl;
        fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "narrow range/left normalization" << endl;

        if(fCrysFitting==0){
            fFileErrLog << "Using exp fit"<<endl;
            FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRangeNarrow,iPt,kFALSE);
            fMesonYieldsResidualBckFuncLeftNarrow[iPt]          = fIntLinearBck;
            fMesonYieldsResidualBckFuncLeftNarrowError[iPt]     = fIntLinearBckError;

        } else {
            fMesonYieldsResidualBckFuncLeftNarrow[iPt]          = 0;
            fMesonYieldsResidualBckFuncLeftNarrowError[iPt]     = 0;
        }
        
        //FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt],fMesonIntDeltaRangeNarrow,iPt,kFALSE);
        fFileDataLog<< "Residual Background leftover narrow integration/left norm in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t" << fMesonYieldsResidualBckFuncLeftNarrow[iPt] 
        << "\t +- \t" << fMesonYieldsResidualBckFuncLeftNarrowError[iPt] << endl << endl;
        fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]       = fMesonYieldsLeftNarrow[iPt]- fMesonYieldsResidualBckFuncLeftNarrow[iPt];
        fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt]  = pow( ( fMesonYieldsLeftNarrowError[iPt] * fMesonYieldsLeftNarrowError[iPt] + fMesonYieldsResidualBckFuncLeftNarrowError[iPt] * 
                                                            fMesonYieldsResidualBckFuncLeftNarrowError[iPt]), 0.5);
        fMesonYieldsLeftPerEventNarrow[iPt]                 = fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/fNEvents;
        fMesonYieldsLeftPerEventNarrowError[iPt]            = fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt]/fNEvents;

        fTotalBckYieldsLeftNarrow[iPt]                      = fBckYieldsLeftNarrow[iPt] + fMesonYieldsResidualBckFuncLeftNarrow[iPt];
        fTotalBckYieldsLeftNarrowError[iPt]                 = pow( fBckYieldsLeftNarrowError[iPt]*fBckYieldsLeftNarrowError[iPt] + fMesonYieldsResidualBckFuncLeftNarrowError[iPt] * 
                                                            fMesonYieldsResidualBckFuncLeftNarrowError[iPt], 0.5);

        if( fTotalBckYieldsLeftNarrow[iPt]!=0){
            fMesonSBLeftNarrow[iPt]         = fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/fTotalBckYieldsLeftNarrow[iPt];
            fMesonSBLeftNarrowError[iPt]    = pow( pow( fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt] / fTotalBckYieldsLeftNarrow[iPt], 2.) + pow(
                                            fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/ (fTotalBckYieldsLeftNarrow[iPt] * fTotalBckYieldsLeftNarrow[iPt]) * fTotalBckYieldsLeftNarrowError[iPt],
                                            2.), 0.5);
            fMesonSignLeftNarrow[iPt]       = fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]/pow(fTotalBckYieldsLeftNarrow[iPt],0.5);
            fMesonSignLeftNarrowError[iPt]  = pow( pow( fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt] / pow( fTotalBckYieldsLeftNarrow[iPt] , 0.5), 2.) + pow( 0.5 *
                                            fMesonYieldsCorResidualBckFuncLeftNarrow[iPt] / pow( fTotalBckYieldsLeftNarrow[iPt], 1.5) * fTotalBckYieldsLeftNarrowError[iPt], 2.) ,0.5);
        }else{
            fMesonSBLeftNarrow[iPt]         = 0.;
            fMesonSBLeftNarrowError[iPt]    = 0.;
            fMesonSignLeftNarrow[iPt]       = 0.;
            fMesonSignLeftNarrowError[iPt]  = 0.;
        }

    }

    //******************** Data OUTPUTFILE ***************************************************
    const char* fileNameSysErrDat = Form("%s/%s/%s_%s_SystematicErrorYieldExtraction_RAWDATA_%s.dat",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix2.Data(), cutSelection.Data());
    fstream fileSysErrDat;
    fileSysErrDat.open(fileNameSysErrDat, ios::out);
    fileSysErrDat << "Calculation of the systematic error due to the yield extraction RAWDATA " << endl;
    fileSysErrDat <<  endl;
    fileSysErrDat << "fGGYields" << endl;
    fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr" << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << fGGYields[iPt] << "+-" << fGGYieldsError[iPt] << "\t" <<
            fGGYieldsWide[iPt] << "+-" << fGGYieldsWideError[iPt] << "\t" <<
            fGGYieldsNarrow[iPt] << "+-" << fGGYieldsNarrowError[iPt] << endl;

    }
    fileSysErrDat <<  endl;
    fileSysErrDat << "fTotalBckYields" << endl;
    fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
            fTotalBckYields[iPt] << "+-" << fTotalBckYieldsError[iPt] << "\t" <<
            fTotalBckYieldsWide[iPt] << "+-" << fTotalBckYieldsWideError[iPt] << "\t" <<
            fTotalBckYieldsNarrow[iPt] << "+-" << fTotalBckYieldsNarrowError[iPt] << "\t" <<
            fTotalBckYieldsLeft[iPt] << "+-" << fTotalBckYieldsLeftError[iPt]<< "\t" <<
            fTotalBckYieldsLeftWide[iPt]<< "+-" << fTotalBckYieldsLeftWideError[iPt]<< "\t" <<
            fTotalBckYieldsLeftNarrow[iPt]<< "+-" << fTotalBckYieldsLeftNarrowError[iPt] << endl;
    }
    fileSysErrDat <<  endl;
    fileSysErrDat << "fMesonYields" << endl;
    fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
            fMesonYieldsCorResidualBckFunc[iPt] << "+-" << fMesonYieldsCorResidualBckFuncError[iPt] << "\t" <<
            fMesonYieldsCorResidualBckFuncWide[iPt] << "+-" << fMesonYieldsCorResidualBckFuncWideError[iPt] << "\t" <<
            fMesonYieldsCorResidualBckFuncNarrow[iPt] << "+-" << fMesonYieldsCorResidualBckFuncNarrowError[iPt] << "\t" <<
            fMesonYieldsCorResidualBckFuncLeft[iPt]<< "+-" << fMesonYieldsCorResidualBckFuncLeftError[iPt]<< "\t" <<
            fMesonYieldsCorResidualBckFuncLeftWide[iPt]<< "+-" << fMesonYieldsCorResidualBckFuncLeftWideError[iPt]<< "\t" <<
            fMesonYieldsCorResidualBckFuncLeftNarrow[iPt]<< "+-" << fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt] << endl;
    }
    if(fIsMC){
        fileSysErrDat <<  endl;
        fileSysErrDat << "TrueYields" << endl;
        fileSysErrDat << "Bin \t True \t True Wide \t True Narr " << endl;
        for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
            fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
                fMesonTrueYields[iPt] << "\t" <<
                fMesonTrueYieldsWide[iPt] << "\t" <<
                fMesonTrueYieldsNarrow[iPt] << endl;
        }
    }
    fileSysErrDat.close();
    //******************************** OUTPUT END ******************************************************
    TString plotPrefix  = Form("%s/%s_%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data());
    TString plotSuffix  = Form("%s_%s.%s",fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    
    TString nameMeson   = Form("%s_MesonWithBck%s", plotPrefix.Data(), plotSuffix.Data());
    TString nameCanvas  = "MesonWithBckCanvas";
    TString namePad     = "MesonWithBckPad";
    cout << nameMeson.Data() << endl;
    PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingBackNormInvMassPtBin, nameMeson, nameCanvas, namePad, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, 
                        fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcess, fCollisionSystem);
    
    TString nameMesonSub    = "";
    TString nameCanvasSub   = "";
    TString namePadSub      = "";
    if (fCrysFitting == 0){
        nameMesonSub    = Form("%s_MesonSubtracted%s", plotPrefix.Data(), plotSuffix.Data());
        nameCanvasSub   = "MesonCanvasSubtracted";
        namePadSub      = "MesonPadSubtracted";
        cout << nameMesonSub.Data() << endl;
//         if( fMode == 4 ) {
//             PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassPtBin2, nameMesonSub, nameCanvasSub, namePadSub,
//                                                 fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess,
//                                                 fCollisionSystem, "MC validated");
//         } else {
            PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassPtBin, nameMesonSub, nameCanvasSub, namePadSub,
                                                fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess,
                                                fCollisionSystem, "MC validated");
//         }           
        nameMesonSub    = Form("%s_MesonSubtractedWithFits%s", plotPrefix.Data(), plotSuffix.Data());
        nameCanvasSub   = "MesonCanvasSubtractedWithFits";
        namePadSub      = "MesonPadSubtractedWithFits";
        PlotWith2FitsSubtractedInvMassInPtBins( fHistoMappingSignalInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassPtBin, fFitBckInvMassPtBin, nameMesonSub, nameCanvasSub, namePadSub,
                                            fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess,
                                            fCollisionSystem, "MC validated");
            
        nameMesonSub    = Form("%s_MesonSubtractedPureGaussianFit%s", plotPrefix.Data(), plotSuffix.Data());
        nameCanvasSub   = "MesonCanvasSubtractedPureGaussianFit";
        namePadSub      = "MesonPadSubtractedPureGaussianFit";
        cout << nameMesonSub.Data() << endl;
        PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalGaussianInvMassPtBin, nameMesonSub, nameCanvasSub, namePadSub,
                                            fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess,
                                            fCollisionSystem, "MC validated");
        
        nameMesonSub    = Form("%s_MesonSubtractedLeft%s", plotPrefix.Data(), plotSuffix.Data());
        nameCanvasSub   = "MesonCanvasSubtractedLeft";
        namePadSub      = "MesonPadSubtractedLeft";
        cout << nameMesonSub.Data() << endl;
//         if( fMode == 4) {
//             PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassLeftPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitInvMassLeftPtBin2, nameMesonSub, nameCanvasSub, namePadSub,
//                                                 fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess,
//                                                 fCollisionSystem, "MC validated");
// 
//             PlotExampleInvMassBins(fHistoMappingGGInvMassPtBin[fExampleBin], fHistoMappingSignalInvMassPtBin[fExampleBin], fHistoMappingBackNormInvMassPtBin[fExampleBin],
//                                 fFitSignalInvMassPtBin2[fExampleBin], fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassPlotRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
//                                 fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess);
//                                 
//                             ///////Added///////
//             nameMesonSub    = Form("%s_MesonSubtractednew%s", plotPrefix.Data(), plotSuffix.Data());
//             nameCanvasSub   = "MesonCanvasSubtractednew";
//             namePadSub      = "MesonPadSubtractednew";
//             cout << nameMesonSub.Data() << endl;
//             PlotWithBGFitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin, fHistoMappingRemainingBGInvMassPtBin, fHistoMappingSignalRemainingBGSubInvMassPtBin, fFitRemainingBGInvMassPtBin,
//                                                 nameMesonSub, nameCanvasSub, namePadSub, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,
//                                                 fDecayChannel, fDetectionProcess, fCollisionSystem);
// 
//             nameMesonSub    = Form("%s_MesonSubtractedLeftnew%s", plotPrefix.Data(), plotSuffix.Data());
//             nameCanvasSub   = "MesonCanvasSubtractednew";
//             namePadSub      = "MesonPadSubtractednew";
//             cout << nameMesonSub.Data() << endl;
//             PlotWithBGFitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin, fHistoMappingRemainingBGInvMassLeftPtBin, fHistoMappingSignalRemainingBGSubInvMassLeftPtBin,
//                                                 fFitRemainingBGInvMassLeftPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt,
//                                                 fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem);
// 
//             nameMesonSub    = Form("%s_MesonSubtractedRemaingBGSubtracted%s", plotPrefix.Data(), plotSuffix.Data());
//             nameCanvasSub   = "MesonCanvasSubtracted";
//             namePadSub      = "MesonPadSubtracted";
//             cout << nameMesonSub.Data() << endl;
//             PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalRemainingBGSubInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassPtBin, nameMesonSub, nameCanvasSub, namePadSub,
//                                                 fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess,
//                                                 fCollisionSystem,"MC validated");
//             
//             nameMesonSub    = Form("%s_MesonSubtractedRemaingBGSubtractedLeft%s", plotPrefix.Data(), plotSuffix.Data());
//             nameCanvasSub   = "MesonCanvasSubtractedLeft";
//             namePadSub      = "MesonPadSubtractedLeft";
//             cout << nameMesonSub.Data() << endl;
//             PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalRemainingBGSubInvMassLeftPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitInvMassLeftPtBin, nameMesonSub, nameCanvasSub, namePadSub,
//                                                 fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess,
//                                                 fCollisionSystem, "MC validated");
//             
//         } else {
            PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassLeftPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitInvMassLeftPtBin, nameMesonSub, nameCanvasSub, namePadSub,
                                                fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess,
                                                fCollisionSystem, "MC validated");

            PlotExampleInvMassBins(fHistoMappingGGInvMassPtBin[fExampleBin], fHistoMappingSignalInvMassPtBin[fExampleBin], fHistoMappingBackNormInvMassPtBin[fExampleBin],
                                fFitSignalInvMassPtBin[fExampleBin], fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassPlotRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                                fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess);
//         }
        
    } else {
        nameMesonSub    = Form("%s_MesonSubtracted%s", plotPrefix.Data(), plotSuffix.Data());
        nameCanvasSub   = "MesonCanvasSubtracted";
        namePadSub      = "MesonPadSubtracted";
        cout << nameMesonSub.Data() << endl;
        PlotWithBGFitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin, fHistoMappingRemainingBGInvMassPtBin, fHistoMappingSignalRemainingBGSubInvMassPtBin, fFitRemainingBGInvMassPtBin,
                                            nameMesonSub, nameCanvasSub, namePadSub, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,
                                            fDecayChannel, fDetectionProcess, fCollisionSystem);

        nameMesonSub    = Form("%s_MesonSubtractedLeft%s", plotPrefix.Data(), plotSuffix.Data());
        nameCanvasSub   = "MesonCanvasSubtracted";
        namePadSub      = "MesonPadSubtracted";
        cout << nameMesonSub.Data() << endl;
        PlotWithBGFitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin, fHistoMappingRemainingBGInvMassLeftPtBin, fHistoMappingSignalRemainingBGSubInvMassLeftPtBin,
                                            fFitRemainingBGInvMassLeftPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt,
                                            fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem);

        nameMesonSub    = Form("%s_MesonSubtractedRemaingBGSubtracted%s", plotPrefix.Data(), plotSuffix.Data());
        nameCanvasSub   = "MesonCanvasSubtracted";
        namePadSub      = "MesonPadSubtracted";
        cout << nameMesonSub.Data() << endl;
        PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalRemainingBGSubInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassPtBin, nameMesonSub, nameCanvasSub, namePadSub,
                                            fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess,
                                            fCollisionSystem,"MC validated");
        
        nameMesonSub    = Form("%s_MesonSubtractedRemaingBGSubtractedLeft%s", plotPrefix.Data(), plotSuffix.Data());
        nameCanvasSub   = "MesonCanvasSubtractedLeft";
        namePadSub      = "MesonPadSubtractedLeft";
        cout << nameMesonSub.Data() << endl;
        PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalRemainingBGSubInvMassLeftPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitInvMassLeftPtBin, nameMesonSub, nameCanvasSub, namePadSub,
                                            fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess,
                                            fCollisionSystem, "MC validated");
        
        PlotExampleInvMassBins(fHistoMappingGGInvMassPtBin[fExampleBin], fHistoMappingSignalRemainingBGSubInvMassPtBin[fExampleBin], fHistoMappingBackNormInvMassPtBin[fExampleBin], 
                            fFitSignalInvMassPtBin[fExampleBin], fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassPlotRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                            fThesis, fCollisionSystem, fBinsPt, fDecayChannel);

    }

    nameMeson       = Form("%s_MesonWithBckLeft%s", plotPrefix.Data(), plotSuffix.Data());
    nameCanvas      = "MesonWithBckCanvasLeft";
    namePad         = "MesonWithBckPadLeft";
    cout << nameMeson.Data() << endl;
    PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingBackNormInvMassLeftPtBin, nameMeson, nameCanvas, namePad,  fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin,
                        fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem);
        
    if(fIsMC){
        TString nameMesonTrue   = Form("%s_TrueMesonFitted%s", plotPrefix.Data(), plotSuffix.Data());
        TString nameCanvasTrue  = "TrueMesonCanvasFitted";
        TString namePadTrue     = "TrueMesonPadFitted";
        cout << nameMesonTrue.Data() << endl;
        PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins, fHistoMappingTrueMesonInvMassPtBins, fFitTrueSignalInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue,
                                            fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess,
                                            fCollisionSystem, "MC validated weighted",kFALSE);

        nameMesonTrue           = Form("%s_TrueMesonReweightedFitted%s", plotPrefix.Data(), plotSuffix.Data());
        nameCanvasTrue          = "TrueMesonCanvasReweightedFitted";
        namePadTrue             = "TrueMesonPadReweightedFitted";
        cout << nameMesonTrue.Data() << endl;
        PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtReweightedBins, fHistoMappingTrueMesonInvMassPtReweightedBins, fFitTrueSignalInvMassPtReweightedBin, nameMesonTrue,
                                            nameCanvasTrue, namePadTrue, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel,
                                            fDetectionProcess, fCollisionSystem, "MC validated reweighted",kFALSE);

        nameMesonTrue           = Form("%s_TrueMesonUnweightedFitted%s", plotPrefix.Data(), plotSuffix.Data());
        nameCanvasTrue          = "TrueMesonCanvasUnweightedFitted";
        namePadTrue             = "TrueMesonPadUnweightedFitted";
        cout << nameMesonTrue.Data() << endl;
        PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtUnweightedBins, fHistoMappingTrueMesonInvMassPtUnweightedBins, fFitTrueSignalInvMassPtUnweightedBin, nameMesonTrue,
                                            nameCanvasTrue, namePadTrue, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel,
                                            fDetectionProcess, fCollisionSystem, "MC validated unweighted",kFALSE);
        
        nameMesonTrue           = Form("%s_TrueMesonFittedPureGaussian%s", plotPrefix.Data(), plotSuffix.Data());
        nameCanvasTrue          = "TrueMesonCanvasFittedPureGaussian";
        namePadTrue             = "TrueMesonPadFittedPureGaussian";
        cout << nameMesonTrue.Data() << endl;
        PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins, fHistoMappingTrueMesonInvMassPtBins, fFitTrueSignalGaussianInvMassPtBin, nameMesonTrue, nameCanvasTrue,
                                            namePadTrue, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel,
                                            fDetectionProcess, fCollisionSystem, "MC validated weighted",kFALSE);

        if (meson.Contains("Pi0")){
            nameMesonTrue       = Form("%s_TrueMesonSecondary%s", plotPrefix.Data(), plotSuffix.Data());
            nameCanvasTrue      = "TrueMesonCanvasSec";
            namePadTrue         = "TrueMesonPadSec";
            cout << nameMesonTrue.Data() << endl;
            PlotInvMassSecondaryInPtBins( fHistoMappingTrueMesonInvMassPtBins, fHistoMappingTrueSecMesonInvMassPtBins, fHistoMappingTrueSecFromK0SMesonInvMassPtBins, nameMesonTrue, nameCanvasTrue,
                                        namePadTrue, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess,
                                        fCollisionSystem);
        }
        
        if (fEnableDCMeson){
            nameMesonTrue       = Form("%s_TrueMesonDoubleCounting%s", plotPrefix.Data(), plotSuffix.Data());
            nameCanvasTrue      = "TrueMesonCanvasDC";
            namePadTrue         = "TrueMesonPadDC";
            cout << nameMesonTrue.Data() << endl;
            PlotInvMassDoubleCountingInPtBins( fHistoMappingTrueMesonInvMassPtBins, fHistoMappingTrueMesonDCInvMassPtBins, nameMesonTrue, nameCanvasTrue,
                                        namePadTrue, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess,
                                        fCollisionSystem);
        }
        
        if (fAdvancedMesonQA && (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5)){
            nameMesonTrue       = Form("%s_TrueMesonCaloPhoton%s", plotPrefix.Data(), plotSuffix.Data());
            nameCanvasTrue      = "TrueMesonCaloPhotonCanvasFitted";
            namePadTrue         = "TrueMesonCaloPhotonPadFitted";
            cout << nameMesonTrue.Data() << endl;
            PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonCaloPhotonInvMassPtBins, fHistoMappingTrueMesonCaloPhotonInvMassPtBins, fFitTrueSignalCaloPhotonInvMassPtBin, nameMesonTrue,
                                                nameCanvasTrue, namePadTrue, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,    
                                                fDecayChannel, fDetectionProcess, fCollisionSystem, "validated #gamma#gamma",kFALSE);

            nameMesonTrue       = Form("%s_TrueMesonCaloConvPhoton%s", plotPrefix.Data(), plotSuffix.Data());
            nameCanvasTrue      = "TrueMesonCaloConvPhotonCanvasFitted";
            namePadTrue         = "TrueMesonCaloConvPhotonPadFitted";
            cout << nameMesonTrue.Data() << endl;
            PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins, fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins, fFitTrueSignalCaloConvPhotonInvMassPtBin,
                                                nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement,
                                                fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem, "validated #gamma_{conv}#gamma_{conv}",kFALSE);

            nameMesonTrue       = Form("%s_TrueMesonCaloMergedCluster%s", plotPrefix.Data(), plotSuffix.Data());
            nameCanvasTrue      = "TrueMesonCaloMergedClusterCanvasFitted";
            namePadTrue         = "TrueMesonCaloMergedClusterPadFitted";
            cout << nameMesonTrue.Data() << endl;
            PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins, fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins,
                                                fFitTrueSignalCaloMergedClusterInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn,
                                                fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem, "validated #gamma's merged",kFALSE);

            nameMesonTrue       = Form("%s_TrueMesonCaloMergedClusterPartConv%s", plotPrefix.Data(), plotSuffix.Data());
            nameCanvasTrue      = "TrueMesonCaloMergedClusterPartConvCanvasFitted";
            namePadTrue         = "TrueMesonCaloMergedClusterPartConvPadFitted";
            cout << nameMesonTrue.Data() << endl;
            PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins, fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins,
                                                fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn,
                                                fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem, "val. #gamma's mer., part. conv",
                                                kFALSE);


            nameMesonTrue       = Form("%s_TrueMesonDecomposedMerged%s", plotPrefix.Data(), plotSuffix.Data());
            cout << nameMesonTrue.Data() << endl;
            PlotTrueInvMassSplittedInMergedInPtBins(fHistoMappingTrueFullMesonInvMassPtBins, fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins,
                                                    fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassPlotRange, fdate, fPrefix,
                                                    fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem, fMode);
            
        }
        cout << "line" << __LINE__ << endl;
        if (fAdvancedMesonQA && (fMode == 4 || fMode == 5)){
            nameMesonTrue       = Form("%s_TrueMesonMixedCaloConvPhoton%s", plotPrefix.Data(), plotSuffix.Data());
            nameCanvasTrue      = "TrueMesonMixedCaloConvPhotonCanvasFitted";
            namePadTrue         = "TrueMesonMixedCaloConvPhotonPadFitted";
            cout << nameMesonTrue.Data() << endl;
            PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins, fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins,
                                                fFitTrueSignalMixedCaloConvPhotonInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn,
                                                fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem, "validated #gamma#gamma_{conv}",
                                                kFALSE);
            
            nameMesonTrue       = Form("%s_TrueMesonDecomposedPhotonsAndElectron%s", plotPrefix.Data(), plotSuffix.Data());
            cout << nameMesonTrue.Data() << endl;
            PlotTrueInvMassSplittedInPhotonAndElectronInPtBins(fHistoMappingTrueFullMesonInvMassPtBins, fHistoMappingTrueMesonCaloPhotonInvMassPtBins,  NULL,
                                                            fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins, fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins, nameMesonTrue,
                                                            nameCanvasTrue, namePadTrue, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement,
                                                            fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem, fMode);

            PlotExampleInvMassBinsMC(fHistoMappingTrueFullMesonInvMassPtBins[fExampleBin], fHistoMappingTrueMesonCaloPhotonInvMassPtBins[fExampleBin],  NULL,
                                    fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[fExampleBin], fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[fExampleBin], fExampleBin,
                                    outputDir.Data(),Suffix.Data(), fMesonMassPlotRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                                    fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, fMode);

            
        } else if (fAdvancedMesonQA && (fMode == 2 || fMode == 3)){
            nameMesonTrue       = Form("%s_TrueMesonDecomposedPhotonsAndElectron%s", plotPrefix.Data(), plotSuffix.Data());
            cout << nameMesonTrue.Data() << endl;
            PlotTrueInvMassSplittedInPhotonAndElectronInPtBins(fHistoMappingTrueFullMesonInvMassPtBins, fHistoMappingTrueMesonCaloPhotonInvMassPtBins,  NULL,
                                                            fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins, NULL, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassPlotRange, fdate,
                                                            fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem,
                                                            fMode);
            PlotExampleInvMassBinsMC(fHistoMappingTrueFullMesonInvMassPtBins[fExampleBin], fHistoMappingTrueMesonCaloPhotonInvMassPtBins[fExampleBin],  NULL,
                                    fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[fExampleBin], NULL, fExampleBin,
                                    outputDir.Data(),Suffix.Data(), fMesonMassPlotRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                                    fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, fMode);
            
        }
        
    }

    CreatePtHistos();
    FillPtHistos();

    ///************************ Plotting fit control plot ************************///
    ///*********************** Mass integration windows
    TCanvas* canvasMassWindows = new TCanvas("canvasMassWindows","",1550,1200);  // gives the page size
    canvasMassWindows->SetTickx();
    canvasMassWindows->SetTicky();
//   canvasMassWindows->SetLeftMargin(0.13);
//   canvasMassWindows->SetRightMargin(0.02);
//   canvasMassWindows->SetTopMargin(0.02);
//   canvasMassWindows->SetFillColor(0);


    DrawGammaSetMarker(fHistoMassWindowLow, 20, 2., kBlack, kBlack);
    if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0){
        DrawAutoGammaMesonHistos( fHistoMassWindowLow, 
                                "", "p_{T} (GeV/c)", Form("Mass (GeV/c^{2})"), 
                                kFALSE, 3.,0.,  kFALSE,
                                kTRUE, 0.,0.3, 
                                kFALSE, 0., 10.);
    } else {
        DrawAutoGammaMesonHistos( fHistoMassWindowLow, 
                                "", "p_{T} (GeV/c)", Form("Mass (GeV/c^{2})"), 
                                kFALSE, 3.,0.,  kFALSE,
                                kTRUE, 0.3,0.7, 
                                kFALSE, 0., 10.);
    }
    fHistoMassWindowLow->DrawCopy("hist");
    DrawGammaSetMarker(fHistoMassWindowHigh, 20, 2., kBlack, kBlack);
    fHistoMassWindowHigh->SetLineColor(kBlack);
    fHistoMassWindowHigh->DrawCopy("same,hist");
    DrawGammaSetMarker(fHistoMassWindowWideLow, 20, 2., kGreen+2, kGreen+2);
    fHistoMassWindowWideLow->DrawCopy("same,hist");
    DrawGammaSetMarker(fHistoMassWindowWideHigh, 20, 2., kGreen+2, kGreen+2);
    fHistoMassWindowWideHigh->DrawCopy("same,hist");
    DrawGammaSetMarker(fHistoMassWindowNarrowLow, 20, 2., kBlue+2, kBlue+2);
    fHistoMassWindowNarrowLow->DrawCopy("same,hist");
    DrawGammaSetMarker(fHistoMassWindowNarrowHigh, 20, 2., kBlue+2, kBlue+2);
    fHistoMassWindowNarrowHigh->DrawCopy("same,hist");
    
    canvasMassWindows->Update();
    
    TLegend* legendMassWindows = new TLegend(0.15,0.8,0.4,0.95);
    legendMassWindows->SetFillColor(0);
    legendMassWindows->SetLineColor(0);
    legendMassWindows->SetTextSize(0.04);
    legendMassWindows->AddEntry(fHistoMassWindowLow,Form("Standard mass window for %s",fPrefix.Data()),"l");
    legendMassWindows->AddEntry(fHistoMassWindowWideLow,"Wide mass window","l");
    legendMassWindows->AddEntry(fHistoMassWindowNarrowLow, "Narrow mass window","l");
    legendMassWindows->Draw();

    if (fIsMC) canvasMassWindows->SaveAs(Form("%s/%s_MC_MassWindows_%s.%s",outputDir.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));
    else canvasMassWindows->SaveAs(Form("%s/%s_data_MassWindows_%s.%s",outputDir.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));
    
    
    
    ///*********************** Lambda tail
    TCanvas* canvasLambdaTail = new TCanvas("canvasLambdaTail","",1550,1200);  // gives the page size
    canvasLambdaTail->SetTickx();
    canvasLambdaTail->SetTicky();
//   canvasLambdaTail->SetLeftMargin(0.13);
//   canvasLambdaTail->SetRightMargin(0.02);
//   canvasLambdaTail->SetTopMargin(0.02);
//   canvasLambdaTail->SetFillColor(0);

    DrawGammaSetMarker(fHistoLambdaTail, 20, 1., kBlack, kBlack);
    if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0){
        DrawAutoGammaMesonHistos( fHistoLambdaTail, 
                                "", "p_{T} (GeV/c)", "#lambda", 
                                kFALSE, 3.,0.,  kFALSE,
                                kTRUE, 0.,fMesonLambdaTailRange[1]*1.2, 
                                kFALSE, 0., 10.);
    } else {
        DrawAutoGammaMesonHistos( fHistoLambdaTail, 
                                "", "p_{T} (GeV/c)", "#lambda", 
                                kFALSE, 3.,0.,  kFALSE,
                                kTRUE, 5e-3,fMesonLambdaTailRange[1]*1.2, 
                                kFALSE, 0., 10.);
    }
    canvasLambdaTail->Update();
    
    TLegend* legendLambdaTail = new TLegend(0.15,0.8,0.4,0.95);
    legendLambdaTail->SetFillColor(0);
    legendLambdaTail->SetLineColor(0);
    legendLambdaTail->SetTextSize(0.04);
    legendLambdaTail->AddEntry(fHistoLambdaTail,Form("Lambda tail parameter for %s",fPrefix.Data()),"p");
    legendLambdaTail->Draw();

    if (fIsMC) canvasLambdaTail->SaveAs(Form("%s/%s_MC_LambdaTail_%s.%s",outputDir.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));
    else canvasLambdaTail->SaveAs(Form("%s/%s_data_LambdaTail_%s.%s",outputDir.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));


    ///*********************** Mass   
    TCanvas* canvasMesonMass = new TCanvas("canvasMesonMass","",1550,1200);  // gives the page size
    canvasMesonMass->SetTickx();
    canvasMesonMass->SetTicky();
//    canvasMesonMass->SetLeftMargin(0.13);
//    canvasMesonMass->SetRightMargin(0.02);
//    canvasMesonMass->SetTopMargin(0.02);
//    canvasMesonMass->SetFillColor(0);

    DrawGammaSetMarker(fHistoMassMeson, 20, 1., kBlack, kBlack);
    if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0){
        DrawAutoGammaMesonHistos( fHistoMassMeson, 
                                "", "p_{T} (GeV/c)", Form("Mass (GeV/c^{2})"), 
                                kFALSE, 3.,0.,  kFALSE,
                                kTRUE, 0.132,0.140, 
                                kFALSE, 0., 10.);
    } else {
        DrawAutoGammaMesonHistos( fHistoMassMeson, 
                                "", "p_{T} (GeV/c)", Form("Mass (GeV/c^{2})"), 
                                kFALSE, 3.,0.,  kFALSE,
                                kTRUE, 0.52,0.58, 
                                kFALSE, 0., 10.);
    }
    canvasMesonMass->Update();
    
    TLegend* legendMesonMass = new TLegend(0.15,0.8,0.4,0.95);
    legendMesonMass->SetFillColor(0);
    legendMesonMass->SetLineColor(0);
    legendMesonMass->SetTextSize(0.04);
    legendMesonMass->AddEntry(fHistoMassMeson,Form("%s mass",fPrefix.Data()),"p");
    legendMesonMass->Draw();

    if (fIsMC) canvasMesonMass->SaveAs(Form("%s/%s_MC_MesonMass_%s.%s",outputDir.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));
    else canvasMesonMass->SaveAs(Form("%s/%s_data_MesonMass_%s.%s",outputDir.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));

    
    ///*********************** Width
    TCanvas* canvasMesonFWHM = new TCanvas("canvasMesonFWHM","",1550,1200);  // gives the page size
    canvasMesonFWHM->SetTickx();
    canvasMesonFWHM->SetTicky();
//   canvasMesonFWHM->SetLeftMargin(0.13);
//   canvasMesonFWHM->SetRightMargin(0.02);
//   canvasMesonFWHM->SetTopMargin(0.02);
//   canvasMesonFWHM->SetFillColor(0);

    DrawGammaSetMarker(fHistoFWHMMeson, 20, 1., kBlack, kBlack);
    if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0){
        DrawAutoGammaMesonHistos( fHistoFWHMMeson, 
                                    "", "p_{T} (GeV/c)","FWHM (GeV/c^{2})", 
                                    kFALSE, 3.,0., kFALSE,
                                    kTRUE, -0.004, fMesonWidthRange[1]*1.5, 
                                    kFALSE, 0., 10.);
    } else {
        DrawAutoGammaMesonHistos( fHistoFWHMMeson, 
                                    "", "p_{T} (GeV/c)","FWHM (GeV/c^{2})", 
                                    kFALSE, 3.,0., kFALSE,
                                    kTRUE, 0., fMesonWidthRange[1]*1.5, 
                                    kFALSE, 0., 10.);
    }
    canvasMesonFWHM->Update();
    
    TLegend* legendMesonFWHM = new TLegend(0.2,0.12,0.45,0.26);
    legendMesonFWHM->SetFillColor(0);
    legendMesonFWHM->SetLineColor(0);
    legendMesonFWHM->SetTextSize(0.04);
    legendMesonFWHM->AddEntry(fHistoFWHMMeson,Form("%s FWHM",fPrefix.Data()),"p");
    legendMesonFWHM->Draw();

    if (fIsMC) canvasMesonFWHM->SaveAs(Form("%s/%s_MC_MesonFWHM_%s.%s",outputDir.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));
    else canvasMesonFWHM->SaveAs(Form("%s/%s_data_MesonFWHM_%s.%s",outputDir.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));
    
    
    if(fIsMC){
        // Rebin MC histograms for acceptance and input with possible weights
        FillHistosArrayMC(fHistoMCMesonPtWithinAcceptance, fHistoMCMesonPt, fDeltaPt);
        // Rebin MC histograms for acceptance and input without weights
        cout << "came till here" << endl;
        if (fHistoMCMesonPtWithinAcceptanceWOWeights) FillHistosArrayMCWOWeights(fHistoMCMesonPtWithinAcceptanceWOWeights, fHistoMCMesonPtWOWeights, fDeltaPt);
        
        // Calculation of meson acceptance with possible weighted input
        CalculateMesonAcceptance();
        // Calculation of meson acceptance without weights for input
        if (fHistoMCMesonPtWithinAcceptanceWOWeights) CalculateMesonAcceptanceWOWeights();

        // calculate pure rec efficiency as in data with fully unweighted histograms if possible
        // ATTENTION: if unweighted histograms are not available this efficiency should not be used for anything!!!!
        fNameHistoEffi="MesonEffiPt";
        cout << fNameHistoEffi.Data() << endl;
        if (fHistoMCMesonPtWithinAcceptanceWOWeights) CalculateMesonEfficiencyWOWeights(fHistoYieldMeson, fHistoYieldTrueSecMeson, fNameHistoEffi);
            else CalculateMesonEfficiency(fHistoYieldMeson, fHistoYieldTrueSecMeson, fNameHistoEffi);
        fHistoMonteMesonEffiPt = fHistoMCMesonEffiPt;
        fNameHistoEffi="MesonWideEffiPt";
        cout << fNameHistoEffi.Data() << endl;
        if (fHistoMCMesonPtWithinAcceptanceWOWeights) CalculateMesonEfficiencyWOWeights(fHistoYieldMesonWide, fHistoYieldTrueSecMesonWide, fNameHistoEffi);
            else CalculateMesonEfficiency(fHistoYieldMesonWide, fHistoYieldTrueSecMesonWide, fNameHistoEffi);
        fHistoMonteMesonWideEffiPt = fHistoMCMesonEffiPt;
        fNameHistoEffi="MesonNarrowEffiPt";
        cout << fNameHistoEffi.Data() << endl;
        if (fHistoMCMesonPtWithinAcceptanceWOWeights) CalculateMesonEfficiencyWOWeights(fHistoYieldMesonNarrow, fHistoYieldTrueSecMesonNarrow, fNameHistoEffi);
            else CalculateMesonEfficiency(fHistoYieldMesonNarrow, fHistoYieldTrueSecMesonNarrow, fNameHistoEffi);
        fHistoMonteMesonNarrowEffiPt = fHistoMCMesonEffiPt;
        fNameHistoEffi="MesonLeftEffiPt";
        cout << fNameHistoEffi.Data() << endl;
        if (fHistoMCMesonPtWithinAcceptanceWOWeights) CalculateMesonEfficiencyWOWeights(fHistoYieldMesonLeft,fHistoYieldTrueSecMeson,fNameHistoEffi);
            else CalculateMesonEfficiency(fHistoYieldMesonLeft,fHistoYieldTrueSecMeson,fNameHistoEffi);
        fHistoMonteMesonLeftEffiPt = fHistoMCMesonEffiPt;
        fNameHistoEffi="MesonLeftNarrowEffiPt";
        cout << fNameHistoEffi.Data() << endl;
        if (fHistoMCMesonPtWithinAcceptanceWOWeights) CalculateMesonEfficiencyWOWeights(fHistoYieldMesonLeftNarrow,fHistoYieldTrueSecMesonNarrow,fNameHistoEffi);
            else CalculateMesonEfficiency(fHistoYieldMesonLeftNarrow,fHistoYieldTrueSecMesonNarrow,fNameHistoEffi);
        fHistoMonteMesonLeftNarrowEffiPt = fHistoMCMesonEffiPt;
        fNameHistoEffi="MesonLeftWideEffiPt";
        cout << fNameHistoEffi.Data() << endl;
        if (fHistoMCMesonPtWithinAcceptanceWOWeights) CalculateMesonEfficiencyWOWeights(fHistoYieldMesonLeftWide,fHistoYieldTrueSecMesonWide,fNameHistoEffi);
            else CalculateMesonEfficiency(fHistoYieldMesonLeftWide,fHistoYieldTrueSecMesonWide,fNameHistoEffi);
        fHistoMonteMesonLeftWideEffiPt = fHistoMCMesonEffiPt;
        
        // True meson efficiencies for fully unweighted MC, to be compared to MesonEffiPt if those have been created with unweighted histograms
        fNameHistoEffi="TrueMesonEffiPtUnweighted";
        if (fHistoMCMesonPtWithinAcceptanceWOWeights) CalculateMesonEfficiencyWOWeights(fHistoYieldTrueMesonUnweighted,NULL,fNameHistoEffi);
        fHistoMCTrueMesonEffiPtUnweighted = fHistoMCMesonEffiPt; 
        fNameHistoEffi="TrueMesonWideEffiPtUnweighted";
        if (fHistoMCMesonPtWithinAcceptanceWOWeights) CalculateMesonEfficiencyWOWeights(fHistoYieldTrueMesonUnweightedWide,NULL,fNameHistoEffi);
        fHistoMCTrueMesonWideEffiPtUnweighted = fHistoMCMesonEffiPt; 
        fNameHistoEffi="TrueMesonNarrowEffiPtUnweighted";
        if (fHistoMCMesonPtWithinAcceptanceWOWeights) CalculateMesonEfficiencyWOWeights(fHistoYieldTrueMesonUnweightedNarrow,NULL,fNameHistoEffi);
        fHistoMCTrueMesonNarrowEffiPtUnweighted = fHistoMCMesonEffiPt; 
        
        // True meson efficiencies with possibly fully weighted inputs on a meson by meson basis in the aliphysics task, should always be used if you start weighting the MC
        // True Meson (only once case, because no normalization)
        fNameHistoEffi="TrueMesonEffiPt";
        cout << fNameHistoEffi.Data() << endl;
        CalculateMesonEfficiency(fHistoYieldTrueMeson,NULL,fNameHistoEffi);
        fHistoMCTrueMesonEffiPt = fHistoMCMesonEffiPt; 
        fNameHistoEffi="TrueMesonWideEffiPt";
        cout << fNameHistoEffi.Data() << endl;
        CalculateMesonEfficiency(fHistoYieldTrueMesonWide,NULL,fNameHistoEffi);
        fHistoMCTrueMesonWideEffiPt = fHistoMCMesonEffiPt;
        fNameHistoEffi="TrueMesonNarrowEffiPt";
        cout << fNameHistoEffi.Data() << endl;
        CalculateMesonEfficiency(fHistoYieldTrueMesonNarrow,NULL,fNameHistoEffi);
        fHistoMCTrueMesonNarrowEffiPt = fHistoMCMesonEffiPt;

        // True meson efficiencies with possibly fully weighted inputs taking the average weight per inv mass bin in the original binning of the TrueMesonInvMass vs pT plot
        // should give on average the same as TrueMesonEffiPt
        fNameHistoEffi="TrueMesonEffiPtReweighted";
        cout << fNameHistoEffi.Data() << endl;
        CalculateMesonEfficiency(fHistoYieldTrueMesonReweighted,NULL,fNameHistoEffi);
        fHistoMCTrueMesonEffiPtReweighted = fHistoMCMesonEffiPt; 
        fNameHistoEffi="TrueMesonWideEffiPtReweighted";
        cout << fNameHistoEffi.Data() << endl;
        CalculateMesonEfficiency(fHistoYieldTrueMesonReweightedWide,NULL,fNameHistoEffi);
        fHistoMCTrueMesonWideEffiPtReweighted = fHistoMCMesonEffiPt;
        fNameHistoEffi="TrueMesonNarrowEffiPtReweighted";
        cout << fNameHistoEffi.Data() << endl;
        CalculateMesonEfficiency(fHistoYieldTrueMesonReweightedNarrow,NULL,fNameHistoEffi);
        fHistoMCTrueMesonNarrowEffiPtReweighted = fHistoMCMesonEffiPt;

        // Calculation of secondary fractions using unweighted histograms, as secondaries are never weighted
        fNameHistoFrac="TrueSecFrac";
        TH1D* fHistoYieldTrueMesonSecPlusPrim = (TH1D*)fHistoYieldTrueMesonUnweighted->Clone("fHistoYieldTrueMesonSecPlusPrim");
        fHistoYieldTrueMesonSecPlusPrim->Add(fHistoYieldTrueSecMeson);
        fHistoYieldTrueSecFracMeson= CalculateSecondaryFractions(fHistoYieldTrueMesonSecPlusPrim, fHistoYieldTrueSecMeson, fNameHistoFrac);
        fNameHistoFrac="TrueSecFracFromK0S";
        fHistoYieldTrueSecFracFromK0SMeson= CalculateSecondaryFractions(fHistoYieldTrueMesonSecPlusPrim, fHistoYieldTrueSecFromK0SMeson, fNameHistoFrac);
        fNameHistoFrac="TrueSecFracFromLambda";
        fHistoYieldTrueSecFracFromLambdaMeson= CalculateSecondaryFractions(fHistoYieldTrueMesonSecPlusPrim, fHistoYieldTrueSecFromLambdaMeson, fNameHistoFrac);
        fNameHistoFrac="TrueSecFracNarrow";
        TH1D* fHistoYieldTrueMesonSecPlusPrimNarrow = (TH1D*)fHistoYieldTrueMesonUnweightedNarrow->Clone("fHistoYieldTrueMesonSecPlusPrimNarrow");
        fHistoYieldTrueMesonSecPlusPrimNarrow->Add(fHistoYieldTrueSecMesonNarrow);
        fHistoYieldTrueSecFracMesonNarrow= CalculateSecondaryFractions(fHistoYieldTrueMesonSecPlusPrimNarrow, fHistoYieldTrueSecMesonNarrow, fNameHistoFrac);
        fNameHistoFrac="TrueSecFracFromK0SNarrow";
        fHistoYieldTrueSecFracFromK0SMesonNarrow= CalculateSecondaryFractions(fHistoYieldTrueMesonSecPlusPrimNarrow, fHistoYieldTrueSecFromK0SMesonNarrow, fNameHistoFrac);
        fNameHistoFrac="TrueSecFracFromLambdaNarrow";
        fHistoYieldTrueSecFracFromLambdaMesonNarrow= CalculateSecondaryFractions(fHistoYieldTrueMesonSecPlusPrimNarrow, fHistoYieldTrueSecFromLambdaMesonNarrow, fNameHistoFrac);
        fNameHistoFrac="TrueSecFracWide";
        TH1D* fHistoYieldTrueMesonSecPlusPrimWide = (TH1D*)fHistoYieldTrueMesonUnweightedWide->Clone("fHistoYieldTrueMesonSecPlusPrimWide");
        fHistoYieldTrueMesonSecPlusPrimWide->Add(fHistoYieldTrueSecMesonWide);
        fHistoYieldTrueSecFracMesonWide= CalculateSecondaryFractions(fHistoYieldTrueMesonSecPlusPrimWide, fHistoYieldTrueSecMesonWide, fNameHistoFrac);
        fNameHistoFrac="TrueSecFracFromK0SWide";
        fHistoYieldTrueSecFracFromK0SMesonWide= CalculateSecondaryFractions(fHistoYieldTrueMesonSecPlusPrimWide, fHistoYieldTrueSecFromK0SMesonWide, fNameHistoFrac);
        fNameHistoFrac="TrueSecFracFromLambdaWide";
        fHistoYieldTrueSecFracFromLambdaMesonWide= CalculateSecondaryFractions(fHistoYieldTrueMesonSecPlusPrimWide, fHistoYieldTrueSecFromLambdaMesonWide, fNameHistoFrac);

        SaveCorrectionHistos(fCutSelection, fPrefix2);
    }
    SaveHistos(fIsMC, fCutSelection, fPrefix2, UseTHnSparse);
    cout << "line " << __LINE__ << endl;
    fFileErrLog.close();
    cout << "line " << __LINE__ << endl;
    fFileDataLog.close();
    cout << "line " << __LINE__ << endl;
    Delete();
    cout << "line " << __LINE__ << endl;
}


//****************************************************************************
//************** Produce background with proper weighting ********************
//****************************************************************************
void ProduceBckProperWeighting(TList* fBackgroundContainer,TList* fMotherContainer, Bool_t UseTHnSparse){

    if(UseTHnSparse){
        cout << "Using THnSparse for the background" << endl;
        THnSparseF* fSparseMotherZM;
        THnSparseF* fSparseBckZM;
        THnSparseF* fSparseMotherZPsi;
        THnSparseF* fSparseBckZPsi;
        
        fSparseMotherZM = (THnSparseF*)fMotherContainer->FindObject("Back_Mother_InvMass_Pt_z_m");
        fSparseBckZM = (THnSparseF*)fBackgroundContainer->FindObject("Back_Back_InvMass_Pt_z_m");
        if(fSparseMotherZM && fSparseBckZM){
           fUseRPBackground = kFALSE;
           cout << "with ZM bins estimation" << endl;
        }
        fSparseMotherZPsi = (THnSparseF*)fMotherContainer->FindObject("Back_Mother_InvMass_Pt_z_psi");
        fSparseBckZPsi = (THnSparseF*)fBackgroundContainer->FindObject("Back_Back_InvMass_Pt_z_psi");
        if(fSparseMotherZPsi && fSparseBckZPsi){
           fUseRPBackground = kTRUE;
           cout << "with ZPsi bins estimation" << endl;
        }
        
        for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
          if(!fUseRPBackground){
            //with ZM bins estimation
            fHistoWeightsBGZbinVsMbin[iPt] = new  TH2F("BGWeights", "", fSparseMotherZM->GetAxis(2)->GetNbins(),  0, fSparseMotherZM->GetAxis(2)->GetNbins(), 
                                            fSparseMotherZM->GetAxis(3)->GetNbins(),  0, fSparseMotherZM->GetAxis(3)->GetNbins());
            fHistoWeightsBGZbinVsMbin[iPt]->GetYaxis()->SetTitle("M-bins");
            fHistoWeightsBGZbinVsMbin[iPt]->GetXaxis()->SetTitle("Z-bins");
            fHistoWeightsBGZbinVsMbin[iPt]->Sumw2();
            fHistoFillPerEventBGZbinVsMbin[iPt] = new  TH2F("BGPoolsFillstatus", "", fSparseMotherZM->GetAxis(2)->GetNbins(),  0, fSparseMotherZM->GetAxis(2)->GetNbins(), 
                                                fSparseMotherZM->GetAxis(3)->GetNbins(),  0, fSparseMotherZM->GetAxis(3)->GetNbins());
            fHistoFillPerEventBGZbinVsMbin[iPt]->GetYaxis()->SetTitle("M-bins");
            fHistoFillPerEventBGZbinVsMbin[iPt]->GetXaxis()->SetTitle("Z-bins");
            fHistoFillPerEventBGZbinVsMbin[iPt]->Sumw2();
        
            for (Int_t z=0;z < fSparseMotherZM->GetAxis(2)->GetNbins();z++){
                for (Int_t m = 0; m < fSparseMotherZM->GetAxis(3)->GetNbins(); m++){ 
                    // pt
                    fSparseMotherZM->GetAxis(1)->SetRange((fSparseMotherZM->GetAxis(1))->FindBin(fBinsPt[iPt]+0.001),(fSparseMotherZM->GetAxis(1))->FindBin(fBinsPt[iPt+1]-0.001));
                    fSparseBckZM->GetAxis(1)->SetRange((fSparseBckZM->GetAxis(1))->FindBin(fBinsPt[iPt]+0.001),(fSparseBckZM->GetAxis(1))->FindBin(fBinsPt[iPt+1]-0.001));
                    // z
                    fSparseMotherZM->GetAxis(2)->SetRange(z, z);
                    fSparseBckZM->GetAxis(2)->SetRange(z, z);
                    // m
                    fSparseMotherZM->GetAxis(3)->SetRange(m,m);
                    fSparseBckZM->GetAxis(3)->SetRange(m,m);
                    
                    fHistoMotherZMProj = (TH1D*)fSparseMotherZM->Projection(0);
                    fHistoMotherZMProj->Sumw2();
                    fHistoBckZMProj = (TH1D*)fSparseBckZM->Projection(0);
                    fHistoBckZMProj->Sumw2();
                
                    fScalingFactorBck[z][m]= 1./fBackgroundMultNumber;
                    if (m==0 && z ==0){
                        if(fHistoMappingBackInvMassPtBin[iPt]!= NULL){
                            delete fHistoMappingBackInvMassPtBin[iPt];
                            fHistoMappingBackInvMassPtBin[iPt]=NULL;
                        }
                        fNameHistoBack = Form("Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
                        fHistoMappingBackInvMassPtBin[iPt]= (TH1D*)fHistoBckZMProj->Clone(fNameHistoBack);
                        fHistoMappingBackInvMassPtBin[iPt]->Sumw2();
                        for (Int_t ii = 0; ii < fHistoBckZMProj->GetNbinsX()+1; ii++){
                            fHistoMappingBackInvMassPtBin[iPt]->SetBinContent(ii,0.);
                            fHistoMappingBackInvMassPtBin[iPt]->SetBinError(ii,0.);
                        }
                    }
                    Int_t startBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[0]);
                    Int_t endBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[1]);
                    if (fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral) != 0) {
                        fScalingFactorBck[z][m] = fHistoMotherZMProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral);
                        if ( fScalingFactorBck[z][m]> (20./fBackgroundMultNumber) ){
                            fScalingFactorBck[z][m]=1./fBackgroundMultNumber;
                        }
                    }
                    fHistoMappingBackInvMassPtBin[iPt]->Add(fHistoBckZMProj,fScalingFactorBck[z][m]);
                    fHistoWeightsBGZbinVsMbin[iPt]->Fill(z+0.5,m+0.5,fScalingFactorBck[z][m]);
                    fHistoFillPerEventBGZbinVsMbin[iPt]->Fill(z+0.5,m+0.5,fHistoBckZMProj->GetEntries());
                    fHistoMotherZMProj->Clear();
                    fHistoBckZMProj->Clear();
                }
            }
            fHistoMappingBackInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
            for (Int_t ii = 0; ii < fHistoMappingBackInvMassPtBin[iPt]->GetNbinsX()+1; ii++){
                if(fHistoMappingBackInvMassPtBin[iPt]->GetBinContent(ii) == 0){
                    fHistoMappingBackInvMassPtBin[iPt]->SetBinContent(ii,0.);
                    fHistoMappingBackInvMassPtBin[iPt]->SetBinError(ii,1.);
                }
            }
            fFileDataLog << "Scaling Background factors for Pt bin " << iPt << " z m " << endl;
            for (Int_t z=0; z < fSparseMotherZM->GetAxis(2)->GetNbins(); z++){
                fFileDataLog << fScalingFactorBck[z][0] << "\t" << fScalingFactorBck[z][1] << "\t" << fScalingFactorBck[z][2] << "\t" << fScalingFactorBck[z][3] << endl;
            }
            
          } else {
            //with ZPsi bins estimation
            fHistoWeightsBGZbinVsPsibin[iPt] = new  TH2F("BGWeights", "", fSparseMotherZPsi->GetAxis(2)->GetNbins(),  0, fSparseMotherZPsi->GetAxis(2)->GetNbins(), 
                                            fSparseMotherZPsi->GetAxis(3)->GetNbins(),  0, fSparseMotherZPsi->GetAxis(3)->GetNbins());
            fHistoWeightsBGZbinVsPsibin[iPt]->GetYaxis()->SetTitle("Psi-bins");
            fHistoWeightsBGZbinVsPsibin[iPt]->GetXaxis()->SetTitle("Z-bins");
            fHistoWeightsBGZbinVsPsibin[iPt]->Sumw2();
            fHistoFillPerEventBGZbinVsPsibin[iPt] = new  TH2F("BGPoolsFillstatus", "", fSparseMotherZPsi->GetAxis(2)->GetNbins(),  0, fSparseMotherZPsi->GetAxis(2)->GetNbins(), 
                                                fSparseMotherZPsi->GetAxis(3)->GetNbins(),  0, fSparseMotherZPsi->GetAxis(3)->GetNbins());
            fHistoFillPerEventBGZbinVsPsibin[iPt]->GetYaxis()->SetTitle("Psi-bins");
            fHistoFillPerEventBGZbinVsPsibin[iPt]->GetXaxis()->SetTitle("Z-bins");
            fHistoFillPerEventBGZbinVsPsibin[iPt]->Sumw2();

            for (Int_t z=0;z < fSparseMotherZPsi->GetAxis(2)->GetNbins();z++){
                for (Int_t psi = 0; psi < fSparseMotherZPsi->GetAxis(3)->GetNbins(); psi++){ 
                    // pt
                    fSparseMotherZPsi->GetAxis(1)->SetRange((fSparseMotherZPsi->GetAxis(1))->FindBin(fBinsPt[iPt]+0.001),(fSparseMotherZPsi->GetAxis(1))->FindBin(fBinsPt[iPt+1]-0.001));
                    fSparseBckZPsi->GetAxis(1)->SetRange((fSparseBckZPsi->GetAxis(1))->FindBin(fBinsPt[iPt]+0.001),(fSparseBckZPsi->GetAxis(1))->FindBin(fBinsPt[iPt+1]-0.001));
                    // z
                    fSparseMotherZPsi->GetAxis(2)->SetRange(z, z);
                    fSparseBckZPsi->GetAxis(2)->SetRange(z, z);
                    // psi
                    fSparseMotherZPsi->GetAxis(3)->SetRange(psi,psi);
                    fSparseBckZPsi->GetAxis(3)->SetRange(psi,psi);

                    fHistoMotherZPsiProj = (TH1D*)fSparseMotherZPsi->Projection(0);
                    fHistoMotherZPsiProj->Sumw2();
                    fHistoBckZPsiProj = (TH1D*)fSparseBckZPsi->Projection(0);
                    fHistoBckZPsiProj->Sumw2();
                
                    fScalingFactorBck[z][psi]= 1./fBackgroundMultNumber;
                    if (psi==0 && z ==0){
                        if(fHistoMappingBackInvMassPtBin[iPt]!= NULL){
                            delete fHistoMappingBackInvMassPtBin[iPt];
                            fHistoMappingBackInvMassPtBin[iPt]=NULL;
                        }
                        fNameHistoBack = Form("Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
                        fHistoMappingBackInvMassPtBin[iPt]= (TH1D*)fHistoBckZPsiProj->Clone(fNameHistoBack);
                        fHistoMappingBackInvMassPtBin[iPt]->Sumw2();
                        for (Int_t ii = 0; ii < fHistoBckZPsiProj->GetNbinsX()+1; ii++){
                            fHistoMappingBackInvMassPtBin[iPt]->SetBinContent(ii,0.);
                            fHistoMappingBackInvMassPtBin[iPt]->SetBinError(ii,0.);
                        }
                    }
                    Int_t startBinIntegral = fHistoMotherZPsiProj->GetXaxis()->FindBin(fBGFitRange[0]);
                    Int_t endBinIntegral = fHistoMotherZPsiProj->GetXaxis()->FindBin(fBGFitRange[1]);
                    if (fHistoBckZPsiProj->Integral(startBinIntegral,endBinIntegral) != 0) {
                        fScalingFactorBck[z][psi] = fHistoMotherZPsiProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZPsiProj->Integral(startBinIntegral,endBinIntegral);
                        if ( fScalingFactorBck[z][psi]> (20./fBackgroundMultNumber) ){
                            fScalingFactorBck[z][psi]=1./fBackgroundMultNumber;
                        }
                    }
                    fHistoMappingBackInvMassPtBin[iPt]->Add(fHistoBckZPsiProj,fScalingFactorBck[z][psi]);
                    fHistoWeightsBGZbinVsPsibin[iPt]->Fill(z+0.5,psi+0.5,fScalingFactorBck[z][psi]);
                    fHistoFillPerEventBGZbinVsPsibin[iPt]->Fill(z+0.5,psi+0.5,fHistoBckZPsiProj->GetEntries());
                    fHistoMotherZPsiProj->Clear();
                    fHistoBckZPsiProj->Clear();
                }
            }
            fHistoMappingBackInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
            for (Int_t ii = 0; ii < fHistoMappingBackInvMassPtBin[iPt]->GetNbinsX()+1; ii++){
                if(fHistoMappingBackInvMassPtBin[iPt]->GetBinContent(ii) == 0){
                    fHistoMappingBackInvMassPtBin[iPt]->SetBinContent(ii,0.);
                    fHistoMappingBackInvMassPtBin[iPt]->SetBinError(ii,1.);
                }
            }
            fFileDataLog << "Scaling Background factors for Pt bin " << iPt << " z psi " << endl;
            for (Int_t z=0; z < fSparseMotherZPsi->GetAxis(2)->GetNbins(); z++){
                fFileDataLog << fScalingFactorBck[z][0] << "\t" << fScalingFactorBck[z][1] << "\t" << fScalingFactorBck[z][2] << "\t" << fScalingFactorBck[z][3] << endl;
            }         
           
          }
        }
        
        if(!fUseRPBackground){
          for (Int_t z=0;z < fSparseMotherZM->GetAxis(2)->GetNbins();z++){
              for (Int_t m = 0; m < fSparseMotherZM->GetAxis(3)->GetNbins(); m++) {
                  // pt
                  fSparseMotherZM->GetAxis(1)->SetRange((fSparseMotherZM->GetAxis(1))->FindBin(fFullPt[0]+0.001),(fSparseMotherZM->GetAxis(1))->FindBin(fFullPt[1]-0.001));
                  fSparseBckZM->GetAxis(1)->SetRange((fSparseBckZM->GetAxis(1))->FindBin(fFullPt[0]+0.001),(fSparseBckZM->GetAxis(1))->FindBin(fFullPt[1]-0.001));
                  // z
                  fSparseMotherZM->GetAxis(2)->SetRange(z+1, z+1);
                  fSparseBckZM->GetAxis(2)->SetRange(z+1, z+1);
                  // m
                  fSparseMotherZM->GetAxis(3)->SetRange(m+1,m+1);
                  fSparseBckZM->GetAxis(3)->SetRange(m+1,m+1);
                  
                  fHistoMotherZMProj = (TH1D*)fSparseMotherZM->Projection(0);
                  fHistoMotherZMProj->Sumw2();
                  fHistoBckZMProj = (TH1D*)fSparseBckZM->Projection(0);
                  fHistoBckZMProj->Sumw2();
                  
                  fScalingFactorBck[z][m]= 1./fBackgroundMultNumber;
                  if (m==0 && z ==0){
                      fNameHistoBack = "Mapping_Back_InvMass_FullPt";
                      fMesonFullPtBackground = (TH1D*)fHistoBckZMProj->Clone(fNameHistoBack);
                      fMesonFullPtBackground->Sumw2();
                      for (Int_t ii = 0; ii < fHistoBckZMProj->GetNbinsX()+1; ii++){
                      fMesonFullPtBackground->SetBinContent(ii,0.);
                      fMesonFullPtBackground->SetBinError(ii,0.);
                      }
                  }
                  Int_t startBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[0]);
                  Int_t endBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[1]);
                  if (fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral) != 0) {
                      fScalingFactorBck[z][m] = fHistoMotherZMProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral);
                      if ( fScalingFactorBck[z][m]>20./fBackgroundMultNumber ){
                      fScalingFactorBck[z][m]=1./fBackgroundMultNumber;
                      }
                  }
                  fMesonFullPtBackground->Add(fHistoBckZMProj,fScalingFactorBck[z][m]);
                  fHistoMotherZMProj->Clear();
                  fHistoBckZMProj->Clear();
              }
            }
        } else {
          for (Int_t z=0;z < fSparseMotherZPsi->GetAxis(2)->GetNbins();z++){
            for (Int_t psi = 0; psi < fSparseMotherZPsi->GetAxis(3)->GetNbins(); psi++) {
              // pt
              fSparseMotherZPsi->GetAxis(1)->SetRange((fSparseMotherZPsi->GetAxis(1))->FindBin(fFullPt[0]+0.001),(fSparseMotherZPsi->GetAxis(1))->FindBin(fFullPt[1]-0.001));
              fSparseBckZPsi->GetAxis(1)->SetRange((fSparseBckZPsi->GetAxis(1))->FindBin(fFullPt[0]+0.001),(fSparseBckZPsi->GetAxis(1))->FindBin(fFullPt[1]-0.001));
              // z
              fSparseMotherZPsi->GetAxis(2)->SetRange(z+1, z+1);
              fSparseBckZPsi->GetAxis(2)->SetRange(z+1, z+1);
              // psi
              fSparseMotherZPsi->GetAxis(3)->SetRange(psi+1,psi+1);
              fSparseBckZPsi->GetAxis(3)->SetRange(psi+1,psi+1);
              
              fHistoMotherZPsiProj = (TH1D*)fSparseMotherZPsi->Projection(0);
              fHistoMotherZPsiProj->Sumw2();
              fHistoBckZPsiProj = (TH1D*)fSparseBckZPsi->Projection(0);
              fHistoBckZPsiProj->Sumw2();
              
              fScalingFactorBck[z][psi]= 1./fBackgroundMultNumber;
              if (psi==0 && z ==0){
                  fNameHistoBack = "Mapping_Back_InvMass_FullPt";
                  fMesonFullPtBackground = (TH1D*)fHistoBckZPsiProj->Clone(fNameHistoBack);
                  fMesonFullPtBackground->Sumw2();
                  for (Int_t ii = 0; ii < fHistoBckZPsiProj->GetNbinsX()+1; ii++){
                  fMesonFullPtBackground->SetBinContent(ii,0.);
                  fMesonFullPtBackground->SetBinError(ii,0.);
                  }
              }
              Int_t startBinIntegral = fHistoMotherZPsiProj->GetXaxis()->FindBin(fBGFitRange[0]);
              Int_t endBinIntegral = fHistoMotherZPsiProj->GetXaxis()->FindBin(fBGFitRange[1]);
              if (fHistoBckZPsiProj->Integral(startBinIntegral,endBinIntegral) != 0) {
                  fScalingFactorBck[z][psi] = fHistoMotherZPsiProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZPsiProj->Integral(startBinIntegral,endBinIntegral);
                  if ( fScalingFactorBck[z][psi]>20./fBackgroundMultNumber ){
                  fScalingFactorBck[z][psi]=1./fBackgroundMultNumber;
                  }
              }
              fMesonFullPtBackground->Add(fHistoBckZPsiProj,fScalingFactorBck[z][psi]);
              fHistoMotherZPsiProj->Clear();
              fHistoBckZPsiProj->Clear();
            }
          }
        }
        fMesonFullPtBackground->Rebin(fNRebin[4]);
        
        if(!fUseRPBackground){
          for (Int_t z=0;z < fSparseMotherZM->GetAxis(2)->GetNbins();z++){
            for (Int_t m = 0; m < fSparseMotherZM->GetAxis(3)->GetNbins(); m++) {
                
                // pt
                fSparseMotherZM->GetAxis(1)->SetRange((fSparseMotherZM->GetAxis(1))->FindBin(fMidPt[0]+0.001),(fSparseMotherZM->GetAxis(1))->FindBin(fMidPt[1]-0.001));
                fSparseBckZM->GetAxis(1)->SetRange((fSparseBckZM->GetAxis(1))->FindBin(fMidPt[0]+0.001),(fSparseBckZM->GetAxis(1))->FindBin(fMidPt[1]-0.001));
                // z
                fSparseMotherZM->GetAxis(2)->SetRange(z+1, z+1);
                fSparseBckZM->GetAxis(2)->SetRange(z+1, z+1);
                // m
                fSparseMotherZM->GetAxis(3)->SetRange(m+1,m+1);
                fSparseBckZM->GetAxis(3)->SetRange(m+1,m+1);
                
                fHistoMotherZMProj = (TH1D*)fSparseMotherZM->Projection(0);
                fHistoMotherZMProj->Sumw2();
                fHistoBckZMProj = (TH1D*)fSparseBckZM->Projection(0);
                fHistoBckZMProj->Sumw2();
                
                fScalingFactorBck[z][m]= 1./fBackgroundMultNumber;
                if (m==0 && z ==0){
                    fNameHistoBack = "Mapping_Back_InvMass_MidPt";
                    fFittingHistMidPtBackground = (TH1D*)fHistoBckZMProj->Clone(fNameHistoBack);
                    fFittingHistMidPtBackground->Sumw2();
                    for (Int_t ii = 0; ii < fHistoBckZMProj->GetNbinsX()+1; ii++){
                    fFittingHistMidPtBackground->SetBinContent(ii,0.);
                    fFittingHistMidPtBackground->SetBinError(ii,0.);
                    }
                }
                Int_t startBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[0]);
                Int_t endBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[1]);
                if (fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral) != 0) {
                    fScalingFactorBck[z][m] = fHistoMotherZMProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral);
                    if ( fScalingFactorBck[z][m]>20./fBackgroundMultNumber ){
                    fScalingFactorBck[z][m]=1./fBackgroundMultNumber;
                    }
                }
                fFittingHistMidPtBackground->Add(fHistoBckZMProj,fScalingFactorBck[z][m]);
                
                fHistoMotherZMProj->Clear();
                fHistoBckZMProj->Clear();
            }
          }
        } else {
          for (Int_t z=0;z < fSparseMotherZPsi->GetAxis(2)->GetNbins();z++){
            for (Int_t psi = 0; psi < fSparseMotherZPsi->GetAxis(3)->GetNbins(); psi++) {
                
                // pt
                fSparseMotherZPsi->GetAxis(1)->SetRange((fSparseMotherZPsi->GetAxis(1))->FindBin(fMidPt[0]+0.001),(fSparseMotherZPsi->GetAxis(1))->FindBin(fMidPt[1]-0.001));
                fSparseBckZPsi->GetAxis(1)->SetRange((fSparseBckZPsi->GetAxis(1))->FindBin(fMidPt[0]+0.001),(fSparseBckZPsi->GetAxis(1))->FindBin(fMidPt[1]-0.001));
                // z
                fSparseMotherZPsi->GetAxis(2)->SetRange(z+1, z+1);
                fSparseBckZPsi->GetAxis(2)->SetRange(z+1, z+1);
                // psi
                fSparseMotherZPsi->GetAxis(3)->SetRange(psi+1,psi+1);
                fSparseBckZPsi->GetAxis(3)->SetRange(psi+1,psi+1);
                
                fHistoMotherZPsiProj = (TH1D*)fSparseMotherZPsi->Projection(0);
                fHistoMotherZPsiProj->Sumw2();
                fHistoBckZPsiProj = (TH1D*)fSparseBckZPsi->Projection(0);
                fHistoBckZPsiProj->Sumw2();
                
                fScalingFactorBck[z][psi]= 1./fBackgroundMultNumber;
                if (psi==0 && z ==0){
                    fNameHistoBack = "Mapping_Back_InvMass_MidPt";
                    fFittingHistMidPtBackground = (TH1D*)fHistoBckZPsiProj->Clone(fNameHistoBack);
                    fFittingHistMidPtBackground->Sumw2();
                    for (Int_t ii = 0; ii < fHistoBckZPsiProj->GetNbinsX()+1; ii++){
                    fFittingHistMidPtBackground->SetBinContent(ii,0.);
                    fFittingHistMidPtBackground->SetBinError(ii,0.);
                    }
                }
                Int_t startBinIntegral = fHistoMotherZPsiProj->GetXaxis()->FindBin(fBGFitRange[0]);
                Int_t endBinIntegral = fHistoMotherZPsiProj->GetXaxis()->FindBin(fBGFitRange[1]);
                if (fHistoBckZPsiProj->Integral(startBinIntegral,endBinIntegral) != 0) {
                    fScalingFactorBck[z][psi] = fHistoMotherZPsiProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZPsiProj->Integral(startBinIntegral,endBinIntegral);
                    if ( fScalingFactorBck[z][psi]>20./fBackgroundMultNumber ){
                    fScalingFactorBck[z][psi]=1./fBackgroundMultNumber;
                    }
                }
                fFittingHistMidPtBackground->Add(fHistoBckZPsiProj,fScalingFactorBck[z][psi]);
                
                fHistoMotherZPsiProj->Clear();
                fHistoBckZPsiProj->Clear();
            }
          }
        }
        fFittingHistMidPtBackground->Rebin(fNRebin[4]);
        
    } else {
        cout << "Using TH2 for the background" << endl;
        
        fHistoMotherZM = (TH2D*)fMotherContainer->FindObject("ESD_Mother_InvMass_Pt");
        fHistoMotherZM->Sumw2();
        fHistoBckZM = (TH2D*)fBackgroundContainer->FindObject("ESD_Background_InvMass_Pt");
        fHistoBckZM->Sumw2();
        
        for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){

            Int_t startBin = fHistoMotherZM->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
            Int_t endBin = fHistoMotherZM->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

            fHistoMotherZMProj = fHistoMotherZM->ProjectionX("ProjectMother",startBin,endBin);
            fHistoMotherZMProj->Sumw2();
            fHistoBckZMProj = fHistoBckZM->ProjectionX("ProjectBck",startBin,endBin);
            fHistoBckZMProj->Sumw2();

            fScalingFactorBck[0][0]= 1./fBackgroundMultNumber;
                if(fHistoMappingBackInvMassPtBin[iPt]!= NULL){
                    delete fHistoMappingBackInvMassPtBin[iPt];
                    fHistoMappingBackInvMassPtBin[iPt]=NULL;
                }
            fNameHistoBack = Form("Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
            fHistoMappingBackInvMassPtBin[iPt]= (TH1D*)fHistoBckZMProj->Clone(fNameHistoBack);
            fHistoMappingBackInvMassPtBin[iPt]->Sumw2();
            
            Int_t startBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[0]);
            Int_t endBinIntegral = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[1]);
            if (fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral) != 0) {
                fScalingFactorBck[0][0] = fHistoMotherZMProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral);
                if ( fScalingFactorBck[0][0]>20./fBackgroundMultNumber ){
                    fScalingFactorBck[0][0]=1./fBackgroundMultNumber;
                }
            }
            fHistoMappingBackInvMassPtBin[iPt]->Add(fHistoBckZMProj,fScalingFactorBck[0][0]);
            fHistoMappingBackInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
            //fHistoMappingBackInvMassPtBin[iPt]->Scale(1./fNRebin[iPt]);
            for (Int_t ii = 0; ii < fHistoMappingBackInvMassPtBin[iPt]->GetNbinsX()+1; ii++){
                if(fHistoMappingBackInvMassPtBin[iPt]->GetBinContent(ii) == 0){
                    fHistoMappingBackInvMassPtBin[iPt]->SetBinContent(ii,0.);
                    fHistoMappingBackInvMassPtBin[iPt]->SetBinError(ii,1.);
                }
            }

            fFileDataLog << "Scaling Background factors for Pt bin " << iPt << endl;
            fFileDataLog << fScalingFactorBck[0][0] << endl;
            
        }

        
        Int_t startBinFullPt = fHistoMotherZM->GetYaxis()->FindBin(fFullPt[0]+0.001);
        Int_t endBinFullPt = fHistoMotherZM->GetYaxis()->FindBin(fFullPt[1]-0.001);

        fHistoMotherZMProjFullPt = fHistoMotherZM->ProjectionX("ProjectMother",startBinFullPt,endBinFullPt);
        fHistoMotherZMProjFullPt->Sumw2();
        fHistoBckZMProjFullPt = fHistoBckZM->ProjectionX("ProjectBck",startBinFullPt,endBinFullPt);
        fHistoBckZMProjFullPt->Sumw2();

        fScalingFactorBckFullPt= 1./fBackgroundMultNumber;
        fNameHistoBack = "Mapping_Back_InvMass_FullPt";
        fMesonFullPtBackground = (TH1D*)fHistoBckZMProjFullPt->Clone(fNameHistoBack);
        fMesonFullPtBackground->Sumw2();
        for (Int_t ii = 0; ii < fHistoBckZMProjFullPt->GetNbinsX()+1; ii++){
                fMesonFullPtBackground->SetBinContent(ii,0.);
                fMesonFullPtBackground->SetBinError(ii,0.);
        }
        Int_t startBinIntegralFullPt = fHistoMotherZMProjFullPt->GetXaxis()->FindBin(fBGFitRange[0]);
        Int_t endBinIntegralFullPt = fHistoMotherZMProjFullPt->GetXaxis()->FindBin(fBGFitRange[1]);
        if (fHistoBckZMProjFullPt->Integral(startBinIntegralFullPt,endBinIntegralFullPt) != 0) {
            fScalingFactorBckFullPt = fHistoMotherZMProjFullPt->Integral(startBinIntegralFullPt,endBinIntegralFullPt)/fHistoBckZMProjFullPt->Integral(startBinIntegralFullPt,endBinIntegralFullPt);
            if ( fScalingFactorBckFullPt>20./fBackgroundMultNumber ){
                fScalingFactorBckFullPt=1./fBackgroundMultNumber;
            }
        }
        fMesonFullPtBackground->Add(fHistoBckZMProjFullPt,fScalingFactorBckFullPt);

        fMesonFullPtBackground->Rebin(fNRebin[4]);
        //fMesonFullPtBackground->Scale(1./fNRebin[4]);

        Int_t startBinMidPt = fHistoMotherZM->GetYaxis()->FindBin(fMidPt[0]+0.001);
        Int_t endBinMidPt = fHistoMotherZM->GetYaxis()->FindBin(fMidPt[1]-0.001);

        fHistoMotherZMProjMidPt = fHistoMotherZM->ProjectionX("ProjectMother",startBinMidPt,endBinMidPt);
        fHistoMotherZMProjMidPt->Sumw2();
        fHistoBckZMProjMidPt = fHistoBckZM->ProjectionX("ProjectBck",startBinMidPt,endBinMidPt);
        fHistoBckZMProjMidPt->Sumw2();

        fScalingFactorBckMidPt= 1./fBackgroundMultNumber;
                
        fNameHistoBack = "Mapping_Back_InvMass_MidPt";
        fFittingHistMidPtBackground = (TH1D*)fHistoBckZMProjMidPt->Clone(fNameHistoBack);
        fFittingHistMidPtBackground->Sumw2();
        for (Int_t ii = 0; ii < fHistoBckZMProjMidPt->GetNbinsX()+1; ii++){
            fFittingHistMidPtBackground->SetBinContent(ii,0.);
            fFittingHistMidPtBackground->SetBinError(ii,0.);
        }
                
        Int_t startBinIntegralMidPt = fHistoMotherZMProjMidPt->GetXaxis()->FindBin(fBGFitRange[0]);
        Int_t endBinIntegralMidPt = fHistoMotherZMProjMidPt->GetXaxis()->FindBin(fBGFitRange[1]);
        if (fHistoBckZMProjMidPt->Integral(startBinIntegralMidPt,endBinIntegralMidPt) != 0) {
            fScalingFactorBckMidPt = fHistoMotherZMProjMidPt->Integral(startBinIntegralMidPt,endBinIntegralMidPt)/fHistoBckZMProjMidPt->Integral(startBinIntegralMidPt,endBinIntegralMidPt);
            if ( fScalingFactorBckMidPt>20./fBackgroundMultNumber ){
                fScalingFactorBckMidPt=1./fBackgroundMultNumber;
            }
        }
        fFittingHistMidPtBackground->Add(fHistoBckZMProjMidPt,fScalingFactorBckMidPt);

        fFittingHistMidPtBackground->Rebin(fNRebin[4]);
        //fFittingHistMidPtBackground->Scale(1./fNRebin[4]);
    }
    
}

//****************************************************************************
//****** Initialization of arrays and variables, histograms for analysis *****
//****** depending on mesonType, number of bins, mode, energy, centrality ****
//****************************************************************************
void Initialize(TString setPi0, Int_t numberOfBins, Int_t triggerSet){
    
    InitializeBinning(setPi0, numberOfBins, fEnergyFlag, fdirectphoton, fMode, fEventCutSelection, fClusterCutSelection, triggerSet);
    
    if (setPi0.CompareTo("Pi0") == 0 || setPi0.CompareTo("Pi0EtaBinning") == 0){
        fPeakRange                  = new Double_t[2];
        fPeakRange[0]               = 0.1;
        fPeakRange[1]               = 0.145; 
        fIntFixedRange              = new Double_t[2];
        fIntFixedRange[0]           = 0.08; 
        fIntFixedRange[1]           = 0.2; 
        fFitRange                   = new Double_t[2];
        fFitRange[0]                = 0.05;
        fFitRange[1]                = 0.25; 
        fBGFitRange                 = new Double_t[2];
        fBGFitRange[0]              = 0.17;
        fBGFitRange[1]              = 0.3; 
        fBGFitRangeLeft             = new Double_t[2];
        fBGFitRangeLeft[0]          = 0.05;
        fBGFitRangeLeft[1]          = 0.08;  
        fMesonPlotRange             = new Double_t[2]; 
        fMesonPlotRange[0]          = 0.13;
        fMesonPlotRange[1]          = 0.138;
        fMesonIntDeltaRange         = new Double_t[2];
        fMesonIntDeltaRange[0]      = -0.035;
        fMesonIntDeltaRange[1]      = 0.010;
        fMesonIntDeltaRangeWide     = new Double_t[2];
        fMesonIntDeltaRangeWide[0]  = -0.055; 
        fMesonIntDeltaRangeWide[1]  = 0.025;
        fMesonIntDeltaRangeNarrow   = new Double_t[2]; 
        fMesonIntDeltaRangeNarrow[0]=-0.015; 
        fMesonIntDeltaRangeNarrow[1]=0.005;
        fMesonMassRange             = new Double_t[2]; 
        fMesonMassRange[0]          = 0.;
        fMesonMassRange[1]          = 0.3;
        if (fMode == 0 && fEnergyFlag.CompareTo("pPb_5.023TeV") == 0)
            fMesonMassRange[1]=   0.14;
        fMesonMassPlotRange         = new Double_t[2]; 
        fMesonMassPlotRange[0]      = 0.; 
        fMesonMassPlotRange[1]      = 0.3;
        fMesonFitRange              = new Double_t[2]; 
        if( fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 || fEnergyFlag.CompareTo("pPb_5.023TeV") == 0){
            fMesonFitRange[0]       = 0.07;
            fMesonFitRange[1]       = 0.22;
        } else { 
            fMesonFitRange[0]       = 0.05; 
            fMesonFitRange[1]       = 0.25; 
        }
        //      fMesonFitRangeWithoutPeak = new Double_t[2]; fMesonFitRangeWithoutPeak[0] = 0.05; fMesonFitRangeWithoutPeak[0] = 0.3;
        fMesonId                    = 111;
        fMesonWidthExpect           = 0.003;
        fMesonLambdaTail            = 0.012;
        fMesonWidthRange            = new Double_t[2];
        fMesonWidthRange[0]         = 0.001; 
        fMesonWidthRange[1]         = 0.009;
        fMesonLambdaTailRange       = new Double_t[2];
        fMesonLambdaTailRange[0]    = 0.001; 
        fMesonLambdaTailRange[1]    = 0.02;
        fMidPt                      = new Double_t[2]; 
        fMidPt[0]                   = 0.8;
        fMidPt[1]                   = 2.5;
        if (fMode == 2){
            fPeakRange[0]                   = 0.05;
            fPeakRange[1]                   = 0.145; 
            fFitRange[0]                    = 0.02;
            fFitRange[1]                    = 0.25; 
            fBGFitRange[0]                  = 0.19;
            fBGFitRange[1]                  = 0.3; 
            fBGFitRangeLeft[0]              = 0.03; 
            fBGFitRangeLeft[1]              = 0.05;  
            fMesonFitRange[0]               = 0.01; 
            fMesonFitRange[1]               = 0.25; 
            fMesonWidthExpect               = 0.005;
            fMesonLambdaTail                = 0.012;
            fMesonWidthRange[0]             = 0.001;
            fMesonWidthRange[1]             = 0.025;
            fMesonLambdaTailRange[0]        = 0.001;
            fMesonLambdaTailRange[1]        = 0.09;
            fMesonMassRange[0]              = 0.;
            fMesonMassRange[1]              = 0.3;
            fMesonMassPlotRange[0]          = 0.;
            fMesonMassPlotRange[1]          = 0.3;
            fMesonIntDeltaRange[0]          = -0.032;
            if( fEnergyFlag.CompareTo("8TeV") == 0 ){ fMesonIntDeltaRange[0] = -0.034;}
            fMesonIntDeltaRange[1]          = 0.022;
            fMesonIntDeltaRangeWide[0]      = -0.048; 
            fMesonIntDeltaRangeWide[1]      = 0.028;
            fMesonIntDeltaRangeNarrow[0]    = -0.016; 
            if( fEnergyFlag.CompareTo("8TeV") == 0 ){ fMesonIntDeltaRangeNarrow[0] = -0.020;}
            fMesonIntDeltaRangeNarrow[1]    = 0.016;

            if ( fEnergyFlag.CompareTo("8TeV") == 0 ){
              TString trigger         = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
              if( trigger.CompareTo("81") == 0 ){
                fFitRange[0]=0.04;
                fMesonFitRange[0] = 0.04;
              }
            }

        }
        if (fMode == 4 || fMode == 5 ){
            fPeakRange[0]                   = 0.05;
            fPeakRange[1]                   = 0.145; 
            fFitRange[0]                    = 0.07; 
            fFitRange[1]                    = 0.25; 
            fBGFitRange[0]                  = 0.19;
            fBGFitRange[1]                  = 0.3; 
            fBGFitRangeLeft[0]              = 0.04;
            fBGFitRangeLeft[1]              = 0.08;  
            fMesonFitRange[0]               = 0.08;
            fMesonFitRange[1]               = 0.25; 
            fMesonWidthExpect               = 0.005;
            fMesonLambdaTail                = 0.012;
            fMesonWidthRange[0]             = 0.001; 
            fMesonWidthRange[1]             = 0.040;
            fMesonLambdaTailRange[0]        = 0.001; 
            fMesonLambdaTailRange[1]        = 0.03;
            fMesonMassRange[0]              = 0.;
            fMesonMassRange[1]              = 0.3;
            fMesonMassPlotRange[0]          = 0.;
            fMesonMassPlotRange[1]          = 0.3;
/*          if ( fMode == 4){
                fMesonIntDeltaRange[0]      = -0.042;
                fMesonIntDeltaRange[1]      = 0.030;
                fMesonIntDeltaRangeWide[0]  = -0.058; 
                fMesonIntDeltaRangeWide[1]  = 0.035;
                fMesonIntDeltaRangeNarrow[0]= -0.026; 
                fMesonIntDeltaRangeNarrow[1]= 0.025;
            }*/
            if ( fMode == 4){
                fPeakRange[0]                   = 0.05;
                fPeakRange[1]                   = 0.17; //0.155; //0.145; 
                fFitRange[0]                    = 0.07; 
                fFitRange[1]                    = 0.25;  
                fBGFitRange[0]                  = 0.19;
                fBGFitRange[1]                  = 0.3; 
                fBGFitRangeLeft[0]              = 0.05; //0.04;
                fBGFitRangeLeft[1]              = 0.08;  
                fMesonFitRange[0]               = 0.06;
                fMesonFitRange[1]               = 0.25; 
                fMesonWidthExpect               = 0.01; //0.005;
                fMesonLambdaTail                = 0.02;
                fMesonWidthRange[0]             = 0.006; //0.001; 
                fMesonWidthRange[1]             = 0.028;
                fMesonLambdaTailRange[0]        = 0.01; // 0.005; //0.001;
                if( fEnergyFlag.CompareTo("8TeV") == 0 ){ fMesonLambdaTailRange[0] = 0.005;}
                fMesonLambdaTailRange[1]        = 0.03;
                fMesonMassRange[0]              = 0.;
                fMesonMassRange[1]              = 0.3;
                fMesonMassPlotRange[0]          = 0.;
                fMesonMassPlotRange[1]          = 0.3;
                fBGFitRangeLeft[0]              = 0.06;
                fBGFitRangeLeft[1]              = 0.08; 
                fMesonIntDeltaRange[0]          = -0.05; 
                fMesonIntDeltaRange[1]          = 0.04; // mod.
                fMesonIntDeltaRangeWide[0]      = fMesonIntDeltaRange[0]*1.4; // mod.
                fMesonIntDeltaRangeWide[1]      = fMesonIntDeltaRange[1]*1.4;// mod.
                fMesonIntDeltaRangeNarrow[0]    = fMesonIntDeltaRange[0]*0.6; // mod.
                fMesonIntDeltaRangeNarrow[1]    = fMesonIntDeltaRange[1]*0.6; // mod.

                if ( fEnergyFlag.CompareTo("8TeV") == 0 ){
                  TString trigger         = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
                  if( trigger.CompareTo("52") == 0 || trigger.CompareTo("81") == 0 ){
                    fBGFitRange[0]=0.21;
                  }
                }
            }
        }
    } else if (setPi0.CompareTo("Eta") == 0){
        fPeakRange                  = new Double_t[2];
        fPeakRange[0]               = 0.48; 
        fPeakRange[1]               = 0.58; 
        fIntFixedRange              = new Double_t[2];
        fIntFixedRange[0]           = 0.48;
        fIntFixedRange[1]           = 0.58; 
        fFitRange                   = new Double_t[2];
        fFitRange[0]                = 0.4;
        fFitRange[1]                = 0.65; 
        fBGFitRange                 = new Double_t[2];
        fBGFitRange[0]              = 0.58; 
        fBGFitRange[1]              = 0.79; 
        fBGFitRangeLeft             = new Double_t[2];
        fBGFitRangeLeft[0]          = 0.35;
        fBGFitRangeLeft[1]          = 0.48;  
        fMesonPlotRange             = new Double_t[2];
        fMesonPlotRange[0]          = 0.53; 
        fMesonPlotRange[1]          = 0.560;
        fMesonIntDeltaRange         = new Double_t[2]; 
        fMesonIntDeltaRange[0]      = -0.048;
        fMesonIntDeltaRange[1]      = 0.022;
        fMesonIntDeltaRangeWide     = new Double_t[2]; 
        fMesonIntDeltaRangeWide[0]  = -0.068;
        fMesonIntDeltaRangeWide[1]  = 0.032;
        fMesonIntDeltaRangeNarrow   = new Double_t[2]; 
        fMesonIntDeltaRangeNarrow[0]= -0.033;
        fMesonIntDeltaRangeNarrow[1]= 0.012;
        fMesonMassRange             = new Double_t[2];
        fMesonMassRange[0]          = 0.35; 
        fMesonMassRange[1]          = 0.79;
        fMesonMassPlotRange         = new Double_t[2]; 
        fMesonMassPlotRange[0]      = 0.35; 
        fMesonMassPlotRange[1]      = 0.79;
        fMesonFitRange              = new Double_t[2]; 
        fMesonFitRange[0]           = 0.4;
        fMesonFitRange[1]           = 0.7;
        fMesonId                    = 221;
        if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0){
            fMesonWidthExpect = 0.010;
        } else {
            fMesonWidthExpect = 0.005;
        }   
        fMesonLambdaTail            = 0.007;
        fMesonWidthRange            = new Double_t[2];
        fMesonWidthRange[0]         = 0.002;   
        fMesonWidthRange[1]         = 0.020;
        fMesonLambdaTailRange       = new Double_t[2];
        fMesonLambdaTailRange[0]    = 0.004;
        fMesonLambdaTailRange[1]    = 0.03;
        fMidPt                      = new Double_t[2]; 
        fMidPt[0]                   = 1.5; 
        fMidPt[1]                   = 2.5;
        if (fMode == 2){
            fPeakRange[0]                   = 0.48; 
            fPeakRange[1]                   = 0.58;
            fFitRange[0]                    = 0.4;
            fFitRange[1]                    = 0.7;
            fBGFitRange[0]                  = 0.65; 
            fBGFitRange[1]                  = 0.75;
            fBGFitRangeLeft[0]              = 0.35;
            fBGFitRangeLeft[1]              = 0.42; 
            if( fEnergyFlag.CompareTo("8TeV") == 0 ){ fBGFitRangeLeft[1] = 0.46;}
            fMesonFitRange[0]               = 0.35; 
            fMesonFitRange[1]               = 0.7;
            fMesonWidthExpect               = 0.020;
            fMesonLambdaTail                = 0.020;
            fMesonWidthRange[0]             = 0.010; 
            fMesonWidthRange[1]             = 0.040;
            fMesonLambdaTailRange[0]        = 0.018;
            if( fEnergyFlag.CompareTo("8TeV") == 0 ){ fMesonLambdaTailRange[0] = 0.005;}
            fMesonLambdaTailRange[1]        = 0.020;
            if( fEnergyFlag.CompareTo("8TeV") == 0 ){ fMesonLambdaTailRange[1] = 0.030;}
            fMesonIntDeltaRange[0]          = -0.060;
            fMesonIntDeltaRange[1]          = 0.055;
            fMesonIntDeltaRangeWide[0]      = -0.080; 
            fMesonIntDeltaRangeWide[1]      = 0.065;
            fMesonIntDeltaRangeNarrow[0]    = -0.040;
            fMesonIntDeltaRangeNarrow[1]    = 0.045;
        }
        if (fMode == 4){
            fPeakRange[0]                   = 0.51;
            fPeakRange[1]                   = 0.59;
            fFitRange[0]                    = 0.41;
            if( fEnergyFlag.CompareTo("8TeV") == 0 ){ fFitRange[0] = 0.39;}
            fFitRange[1]                    = 0.73;
            if( fEnergyFlag.CompareTo("8TeV") == 0 ){ fFitRange[1] = 0.71;}
            fBGFitRange[0]                  = 0.67;
            fBGFitRange[1]                  = 0.799;
            if( fEnergyFlag.CompareTo("8TeV") == 0 ){ fBGFitRange[1] = 0.795;}
            fBGFitRangeLeft[0]              = 0.34;
            fBGFitRangeLeft[1]              = 0.44;
            fMesonFitRange[0]               = 0.41;
            if( fEnergyFlag.CompareTo("8TeV") == 0 ){ fMesonFitRange[0] = 0.39;}
            fMesonFitRange[1]               = 0.73;
            if( fEnergyFlag.CompareTo("8TeV") == 0 ){ fMesonFitRange[1] = 0.71;}
            fMesonWidthExpect               = 0.025;
            fMesonLambdaTail                = 0.012;
            fMesonWidthRange[0]             = 0.018;
            fMesonWidthRange[1]             = 0.070;
            fMesonLambdaTailRange[0]        = 0.001;
            fMesonLambdaTailRange[1]        = 0.025;
            fMesonIntDeltaRange[0]          = -0.080;
            fMesonIntDeltaRange[1]          = 0.080;
            fMesonIntDeltaRangeWide[0]      = -0.10;
            fMesonIntDeltaRangeWide[1]      = 0.10;
            fMesonIntDeltaRangeNarrow[0]    = -0.060;
            fMesonIntDeltaRangeNarrow[1]    = 0.060;
            if( fEnergyFlag.CompareTo("8TeV") == 0 ) {
              fMesonIntDeltaRange[0]          = -0.070;
              fMesonIntDeltaRange[1]          = 0.070;
              fMesonIntDeltaRangeWide[0]      = -0.09;
              fMesonIntDeltaRangeWide[1]      = 0.09;
              fMesonIntDeltaRangeNarrow[0]    = -0.050;
              fMesonIntDeltaRangeNarrow[1]    = 0.050;
            }
        }
    } else if (setPi0.CompareTo("EtaPrim") == 0){
        fPeakRange                      = new Double_t[2]; 
        fPeakRange[0]                   = 0.1;
        fPeakRange[1]                   = 0.145; 
        fIntFixedRange                  = new Double_t[2];
        fIntFixedRange[0]               = 0.08;
        fIntFixedRange[1]               = 0.2; 
        fFitRange                       = new Double_t[2];
        fFitRange[0]                    = 0.05; 
        fFitRange[1]                    = 0.25; 
        fBGFitRange                     = new Double_t[2];
        fBGFitRange[0]                  = 0.978; 
        fBGFitRange[1]                  = 1.; 
        fBGFitRangeLeft                 = new Double_t[2];
        fBGFitRangeLeft[0]              = 0.9; 
        fBGFitRangeLeft[1]              = 0.94;  
        fMesonPlotRange                 = new Double_t[2];
        fMesonPlotRange[0]              = 0.9;
        fMesonPlotRange[1]              = 1.;
        fMesonIntDeltaRange             = new Double_t[2]; 
        fMesonIntDeltaRange[0]          = 0.92;
        fMesonIntDeltaRange[1]          = 0.96;
        fMesonIntDeltaRangeWide         = new Double_t[2]; 
        fMesonIntDeltaRangeWide[0]      = 0.91;
        fMesonIntDeltaRangeWide[1]      = 0.965;
        fMesonIntDeltaRangeNarrow       = new Double_t[2];
        fMesonIntDeltaRangeNarrow[0]    = 0.92; 
        fMesonIntDeltaRangeNarrow[1]    = 0.959;
        fMesonMassRange                 = new Double_t[2];
        fMesonMassRange[0]              = 0.9;
        fMesonMassRange[1]              = 1.;
        fMesonMassPlotRange             = new Double_t[2];
        fMesonMassPlotRange[0]          = 0.9;
        fMesonMassPlotRange[1]          = 1.;;
        fMesonFitRange                  = new Double_t[2];
        fMesonFitRange[0]               = 0.9;
        fMesonFitRange[1]               = 0.98;
        fMesonId                        = 331;
        fMesonWidthExpect               = 0.01;
        fMesonLambdaTail                = 0.007;
        fMesonWidthRange                = new Double_t[2];
        fMesonWidthRange[0]             = 0.002;
        fMesonWidthRange[1]             = 0.02;
        fMesonLambdaTailRange           = new Double_t[2];
        fMesonLambdaTailRange[0]        = 0.0005; 
        fMesonLambdaTailRange[1]        = 0.03;
        fMidPt                          = new Double_t[2];
        fMidPt[0]                       = 1.2; 
        fMidPt[1]                       = 2.5;
    }
    fMesonLambdaTailMC      = fMesonLambdaTail;
    fMesonWidthExpectMC     = fMesonWidthExpect;
    fMesonWidthRangeMC      = new Double_t[2]; 
    fMesonWidthRangeMC[0]   = fMesonWidthRange[0];
    fMesonWidthRangeMC[1]   = fMesonWidthRange[1];
    fFullPt                 = new Double_t[2]; 
    fFullPt[0]              = 0.4;
    fFullPt[1]              = 15.;

    fMesonCurIntRange                                               = new Double_t[2];
    fMesonCurIntRangeWide                                           = new Double_t[2];
    fMesonCurIntRangeNarrow                                         = new Double_t[2];
    fMesonCurLeftIntRange                                           = new Double_t[2];
    fMesonCurLeftIntRangeWide                                       = new Double_t[2];
    fMesonCurLeftIntRangeNarrow                                     = new Double_t[2];
    fMesonTrueIntRange                                              = new Double_t[2];
    fMesonTrueIntRangeWide                                          = new Double_t[2];
    fMesonTrueIntRangeNarrow                                        = new Double_t[2];
    fMesonTrueIntReweightedRange                                    = new Double_t[2];
    fMesonTrueIntReweightedRangeWide                                = new Double_t[2];
    fMesonTrueIntReweightedRangeNarrow                              = new Double_t[2];
    fMesonTrueIntUnweightedRange                                    = new Double_t[2];
    fMesonTrueIntUnweightedRangeWide                                = new Double_t[2];
    fMesonTrueIntUnweightedRangeNarrow                              = new Double_t[2];

    fGGYields                                                       = new Double_t[fNBinsPt];
    fBckYields                                                      = new Double_t[fNBinsPt];
    fTotalBckYields                                                 = new Double_t[fNBinsPt];
    fMesonYields                                                    = new Double_t[fNBinsPt];
    fMesonTrueYields                                                = new Double_t[fNBinsPt];
    fMesonTrueYieldsDC                                              = new Double_t[fNBinsPt];
    fMesonTrueYieldsReweighted                                      = new Double_t[fNBinsPt];
    fMesonTrueYieldsUnweighted                                      = new Double_t[fNBinsPt];
    
    fMesonTrueYieldFixedWindow                                      = new Double_t[fNBinsPt];
    fMesonTrueYieldGammaFixedWindow                                 = new Double_t[fNBinsPt];
    fMesonTrueYieldGammaConvGammaFixedWindow                        = new Double_t[fNBinsPt];
    fMesonTrueYieldConvGammaConvGammaFixedWindow                    = new Double_t[fNBinsPt];
    fMesonTrueYieldErrorFixedWindow                                 = new Double_t[fNBinsPt];
    fMesonTrueYieldGammaErrorFixedWindow                            = new Double_t[fNBinsPt];
    fMesonTrueYieldGammaConvGammaErrorFixedWindow                   = new Double_t[fNBinsPt];
    fMesonTrueYieldConvGammaConvGammaErrorFixedWindow               = new Double_t[fNBinsPt];

    fMesonTrueSecYields                                             = new Double_t[fNBinsPt];
    fMesonTrueSecFromK0SYields                                      = new Double_t[fNBinsPt];
    fMesonTrueSecFromLambdaYields                                   = new Double_t[fNBinsPt];
    fMesonYieldsFunc                                                = new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFunc                                     = new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFunc                                  = new Double_t[fNBinsPt];
    fMesonYieldsPerEvent                                            = new Double_t[fNBinsPt];
    fMesonMass                                                      = new Double_t[fNBinsPt];
    fMesonLambdaTailpar                                             = new Double_t[fNBinsPt];
    fMesonLambdaTailparError                                        = new Double_t[fNBinsPt];
    fMesonLambdaTailMCpar                                           = new Double_t[fNBinsPt];
    fMesonLambdaTailMCparError                                      = new Double_t[fNBinsPt];
    fMesonSB                                                        = new Double_t[fNBinsPt];
    fMesonSBdefault                                                 = new Double_t[fNBinsPt];
    fMesonSign                                                      = new Double_t[fNBinsPt];
    fMesonSigndefault                                               = new Double_t[fNBinsPt];
    fMassWindowHigh                                                 = new Double_t[fNBinsPt];
    fMassWindowLow                                                   = new Double_t[fNBinsPt];
    fMassWindowWideHigh                                             = new Double_t[fNBinsPt];
    fMassWindowWideLow                                              = new Double_t[fNBinsPt];
    fMassWindowNarrowHigh                                           = new Double_t[fNBinsPt];
    fMassWindowNarrowLow                                            = new Double_t[fNBinsPt];
    fMesonFWHM                                                      = new Double_t[fNBinsPt];
    fMesonFWHMAlpha01                                               = new Double_t[fNBinsPt];
    fMesonTrueMass                                                  = new Double_t[fNBinsPt];
    fMesonTrueMassReweighted                                        = new Double_t[fNBinsPt];
    fMesonTrueMassUnweighted                                        = new Double_t[fNBinsPt];
    fMesonTrueFWHM                                                  = new Double_t[fNBinsPt];
    fMesonTrueFWHMReweighted                                        = new Double_t[fNBinsPt];
    fMesonTrueFWHMUnweighted                                        = new Double_t[fNBinsPt];
    fMesonTrueSB                                                    = new Double_t[fNBinsPt];
    fMesonTrueSign                                                  = new Double_t[fNBinsPt];
    
    fMesonTrueMassCaloPhoton                                        = new Double_t[fNBinsPt];
    fMesonTrueMassCaloElectron                                      = new Double_t[fNBinsPt];
    fMesonTrueMassCaloConvPhoton                                    = new Double_t[fNBinsPt];
    fMesonTrueMassCaloMergedCluster                                 = new Double_t[fNBinsPt];
    fMesonTrueMassCaloMergedClusterPartConv                         = new Double_t[fNBinsPt];
    fMesonTrueMassMixedCaloConvPhoton                               = new Double_t[fNBinsPt];
    fMesonTrueFWHMCaloPhoton                                        = new Double_t[fNBinsPt];
    fMesonTrueFWHMCaloElectron                                      = new Double_t[fNBinsPt];
    fMesonTrueFWHMCaloConvPhoton                                    = new Double_t[fNBinsPt];
    fMesonTrueFWHMCaloMergedCluster                                 = new Double_t[fNBinsPt];
    fMesonTrueFWHMCaloMergedClusterPartConv                         = new Double_t[fNBinsPt];
    fMesonTrueFWHMMixedCaloConvPhoton                               = new Double_t[fNBinsPt];
    
    // Normalization at the left of the peak
    fGGYieldsLeft                                                   = new Double_t[fNBinsPt];
    fBckYieldsLeft                                                  = new Double_t[fNBinsPt];
    fTotalBckYieldsLeft                                             = new Double_t[fNBinsPt];
    fMesonYieldsLeft                                                = new Double_t[fNBinsPt];
    fMesonYieldsFuncLeft                                            = new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncLeft                                 = new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncLeft                              = new Double_t[fNBinsPt];
    fMesonYieldsLeftPerEvent                                        = new Double_t[fNBinsPt];
    fMesonMassLeft                                                  = new Double_t[fNBinsPt];
    fMesonSBLeft                                                    = new Double_t[fNBinsPt];
    fMesonSignLeft                                                  = new Double_t[fNBinsPt];
    fMesonFWHMLeft                                                  = new Double_t[fNBinsPt];

    // Narrow Integration Window
    fGGYieldsNarrow                                                 = new Double_t[fNBinsPt];
    fBckYieldsNarrow                                                = new Double_t[fNBinsPt];
    fTotalBckYieldsNarrow                                           = new Double_t[fNBinsPt];
    fMesonYieldsNarrow                                              = new Double_t[fNBinsPt];
    fMesonTrueYieldsNarrow                                          = new Double_t[fNBinsPt];
    fMesonTrueYieldsReweightedNarrow                                = new Double_t[fNBinsPt];
    fMesonTrueYieldsUnweightedNarrow                                = new Double_t[fNBinsPt];
    fMesonTrueSecYieldsNarrow                                       = new Double_t[fNBinsPt];
    fMesonTrueSecFromK0SYieldsNarrow                                = new Double_t[fNBinsPt];
    fMesonTrueSecFromLambdaYieldsNarrow                             = new Double_t[fNBinsPt];
    fMesonYieldsFuncNarrow                                          = new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncNarrow                               = new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncNarrow                            = new Double_t[fNBinsPt];
    fMesonYieldsPerEventNarrow                                      = new Double_t[fNBinsPt];
    fMesonSBNarrow                                                  = new Double_t[fNBinsPt];
    fMesonSBdefaultNarrow                                           = new Double_t[fNBinsPt];
    fMesonSignNarrow                                                = new Double_t[fNBinsPt];
    fMesonSigndefaultNarrow                                         = new Double_t[fNBinsPt];
    
    fGGYieldsLeftNarrow                                             = new Double_t[fNBinsPt];
    fBckYieldsLeftNarrow                                            = new Double_t[fNBinsPt];
    fTotalBckYieldsLeftNarrow                                       = new Double_t[fNBinsPt];
    fMesonYieldsLeftNarrow                                          = new Double_t[fNBinsPt];
    fMesonYieldsFuncLeftNarrow                                      = new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncLeftNarrow                           = new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncLeftNarrow                        = new Double_t[fNBinsPt];
    fMesonYieldsLeftPerEventNarrow                                  = new Double_t[fNBinsPt];
    fMesonSBLeftNarrow                                              = new Double_t[fNBinsPt];
    fMesonSignLeftNarrow                                            = new Double_t[fNBinsPt];

    // Wide Integration Window
    fGGYieldsWide                                                   = new Double_t[fNBinsPt];
    fBckYieldsWide                                                  = new Double_t[fNBinsPt];
    fTotalBckYieldsWide                                             = new Double_t[fNBinsPt];
    fMesonYieldsWide                                                = new Double_t[fNBinsPt];
    fMesonTrueYieldsWide                                            = new Double_t[fNBinsPt];
    fMesonTrueYieldsReweightedWide                                  = new Double_t[fNBinsPt];
    fMesonTrueYieldsUnweightedWide                                  = new Double_t[fNBinsPt];
    fMesonTrueSecYieldsWide                                         = new Double_t[fNBinsPt];
    fMesonTrueSecFromK0SYieldsWide                                  = new Double_t[fNBinsPt];
    fMesonTrueSecFromLambdaYieldsWide                               = new Double_t[fNBinsPt];
    fMesonYieldsFuncWide                                            = new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncWide                                 = new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncWide                              = new Double_t[fNBinsPt];
    fMesonYieldsPerEventWide                                        = new Double_t[fNBinsPt];
    fMesonSBWide                                                    = new Double_t[fNBinsPt];
    fMesonSignWide                                                  = new Double_t[fNBinsPt];

    fGGYieldsLeftWide                                               = new Double_t[fNBinsPt];
    fBckYieldsLeftWide                                              = new Double_t[fNBinsPt];
    fTotalBckYieldsLeftWide                                         = new Double_t[fNBinsPt];
    fMesonYieldsLeftWide                                            = new Double_t[fNBinsPt];
    fMesonYieldsFuncLeftWide                                        = new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncLeftWide                             = new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncLeftWide                          = new Double_t[fNBinsPt];
    fMesonYieldsLeftPerEventWide                                    = new Double_t[fNBinsPt];
    fMesonSBLeftWide                                                = new Double_t[fNBinsPt];
    fMesonSignLeftWide                                              = new Double_t[fNBinsPt];

    fGGYieldsError                                                  = new Double_t[fNBinsPt];
    fBckYieldsError                                                 = new Double_t[fNBinsPt];
    fTotalBckYieldsError                                            = new Double_t[fNBinsPt];
    fMesonYieldsError                                               = new Double_t[fNBinsPt];
    fMesonTrueYieldsError                                           = new Double_t[fNBinsPt];
    fMesonTrueYieldsDCError                                         = new Double_t[fNBinsPt];
    fMesonTrueYieldsReweightedError                                 = new Double_t[fNBinsPt];
    fMesonTrueYieldsUnweightedError                                 = new Double_t[fNBinsPt];
    fMesonTrueSecYieldsError                                        = new Double_t[fNBinsPt];
    fMesonTrueSecFromK0SYieldsError                                 = new Double_t[fNBinsPt];
    fMesonTrueSecFromLambdaYieldsError                              = new Double_t[fNBinsPt];
    fMesonYieldsFuncError                                           = new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncError                                = new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncError                             = new Double_t[fNBinsPt];
    fMesonYieldsPerEventError                                       = new Double_t[fNBinsPt];
    fMesonMassError                                                 = new Double_t[fNBinsPt];
    fMesonSBError                                                   = new Double_t[fNBinsPt];
    fMesonSBdefaultError                                            = new Double_t[fNBinsPt];
    fMesonSignError                                                 = new Double_t[fNBinsPt];
    fMesonSigndefaultError                                          = new Double_t[fNBinsPt];
    fMesonFWHMError                                                 = new Double_t[fNBinsPt];
    fMesonFWHMAlpha01Error                                          = new Double_t[fNBinsPt];
    fMesonTrueMassError                                             = new Double_t[fNBinsPt];
    fMesonTrueMassReweightedError                                   = new Double_t[fNBinsPt];
    fMesonTrueMassUnweightedError                                   = new Double_t[fNBinsPt];
    fMesonTrueFWHMError                                             = new Double_t[fNBinsPt];
    fMesonTrueFWHMReweightedError                                   = new Double_t[fNBinsPt];
    fMesonTrueFWHMUnweightedError                                   = new Double_t[fNBinsPt];
    fMesonTrueSBError                                               = new Double_t[fNBinsPt];
    fMesonTrueSignError                                             = new Double_t[fNBinsPt];
    fMesonTrueMassErrorCaloPhoton                                   = new Double_t[fNBinsPt];
    fMesonTrueMassErrorCaloElectron                                 = new Double_t[fNBinsPt];
    fMesonTrueMassErrorCaloConvPhoton                               = new Double_t[fNBinsPt];
    fMesonTrueMassErrorCaloMergedCluster                            = new Double_t[fNBinsPt];
    fMesonTrueMassErrorCaloMergedClusterPartConv                    = new Double_t[fNBinsPt];
    fMesonTrueMassErrorMixedCaloConvPhoton                          = new Double_t[fNBinsPt];
    fMesonTrueFWHMErrorCaloPhoton                                   = new Double_t[fNBinsPt];
    fMesonTrueFWHMErrorCaloElectron                                 = new Double_t[fNBinsPt];
    fMesonTrueFWHMErrorCaloConvPhoton                               = new Double_t[fNBinsPt];
    fMesonTrueFWHMErrorCaloMergedCluster                            = new Double_t[fNBinsPt];
    fMesonTrueFWHMErrorCaloMergedClusterPartConv                    = new Double_t[fNBinsPt];
    fMesonTrueFWHMErrorMixedCaloConvPhoton                          = new Double_t[fNBinsPt];
    
    fGGYieldsLeftError                                              = new Double_t[fNBinsPt];
    fBckYieldsLeftError                                             = new Double_t[fNBinsPt];
    fTotalBckYieldsLeftError                                        = new Double_t[fNBinsPt];
    fMesonYieldsLeftError                                           = new Double_t[fNBinsPt];
    fMesonYieldsFuncLeftError                                       = new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncLeftError                            = new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncLeftError                         = new Double_t[fNBinsPt];
    fMesonYieldsLeftPerEventError                                   = new Double_t[fNBinsPt];
    fMesonMassLeftError                                             = new Double_t[fNBinsPt];
    fMesonSBLeftError                                               = new Double_t[fNBinsPt];
    fMesonSignLeftError                                             = new Double_t[fNBinsPt];
    fMesonFWHMLeftError                                             = new Double_t[fNBinsPt];

    // Narrow integration Window
    fGGYieldsNarrowError                                            = new Double_t[fNBinsPt];
    fBckYieldsNarrowError                                           = new Double_t[fNBinsPt];
    fTotalBckYieldsNarrowError                                      = new Double_t[fNBinsPt];
    fMesonYieldsNarrowError                                         = new Double_t[fNBinsPt];
    fMesonTrueYieldsNarrowError                                     = new Double_t[fNBinsPt];
    fMesonTrueYieldsReweightedNarrowError                           = new Double_t[fNBinsPt];
    fMesonTrueYieldsUnweightedNarrowError                           = new Double_t[fNBinsPt];
    fMesonTrueSecYieldsNarrowError                                  = new Double_t[fNBinsPt];
    fMesonTrueSecFromK0SYieldsNarrowError                           = new Double_t[fNBinsPt];
    fMesonTrueSecFromLambdaYieldsNarrowError                        = new Double_t[fNBinsPt];
    fMesonYieldsFuncNarrowError                                     = new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncNarrowError                          = new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncNarrowError                       = new Double_t[fNBinsPt];
    fMesonYieldsPerEventNarrowError                                 = new Double_t[fNBinsPt];
    fMesonSBNarrowError                                             = new Double_t[fNBinsPt];
    fMesonSBdefaultNarrowError                                      = new Double_t[fNBinsPt];
    fMesonSigndefaultNarrowError                                    = new Double_t[fNBinsPt];
    fMesonSignNarrowError                                           = new Double_t[fNBinsPt];

    fGGYieldsLeftNarrowError                                        = new Double_t[fNBinsPt];
    fBckYieldsLeftNarrowError                                       = new Double_t[fNBinsPt];
    fTotalBckYieldsLeftNarrowError                                  = new Double_t[fNBinsPt];
    fMesonYieldsLeftNarrowError                                     = new Double_t[fNBinsPt];
    fMesonYieldsFuncLeftNarrowError                                 = new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncLeftNarrowError                      = new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncLeftNarrowError                   = new Double_t[fNBinsPt];
    fMesonYieldsLeftPerEventNarrowError                             = new Double_t[fNBinsPt];
    fMesonSBLeftNarrowError                                         = new Double_t[fNBinsPt];
    fMesonSignLeftNarrowError                                       = new Double_t[fNBinsPt];

    // Wide integration Window
    fGGYieldsWideError                                              = new Double_t[fNBinsPt];
    fBckYieldsWideError                                             = new Double_t[fNBinsPt];
    fTotalBckYieldsWideError                                        = new Double_t[fNBinsPt];
    fMesonYieldsWideError                                           = new Double_t[fNBinsPt];
    fMesonTrueYieldsWideError                                       = new Double_t[fNBinsPt];
    fMesonTrueYieldsReweightedWideError                             = new Double_t[fNBinsPt];
    fMesonTrueYieldsUnweightedWideError                             = new Double_t[fNBinsPt];
    fMesonTrueSecYieldsWideError                                    = new Double_t[fNBinsPt];
    fMesonTrueSecFromK0SYieldsWideError                             = new Double_t[fNBinsPt];
    fMesonTrueSecFromLambdaYieldsWideError                          = new Double_t[fNBinsPt];
    fMesonYieldsFuncWideError                                       = new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncWideError                            = new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncWideError                         = new Double_t[fNBinsPt];
    fMesonYieldsPerEventWideError                                   = new Double_t[fNBinsPt];
    fMesonSBWideError                                               = new Double_t[fNBinsPt];
    fMesonSignWideError                                             = new Double_t[fNBinsPt];

    fGGYieldsLeftWideError                                          = new Double_t[fNBinsPt];
    fBckYieldsLeftWideError                                         = new Double_t[fNBinsPt];
    fTotalBckYieldsLeftWideError                                    = new Double_t[fNBinsPt];
    fMesonYieldsLeftWideError                                       = new Double_t[fNBinsPt];
    fMesonYieldsFuncLeftWideError                                   = new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncLeftWideError                        = new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncLeftWideError                     = new Double_t[fNBinsPt];
    fMesonYieldsLeftPerEventWideError                               = new Double_t[fNBinsPt];
    fMesonSBLeftWideError                                           = new Double_t[fNBinsPt];
    fMesonSignLeftWideError                                         = new Double_t[fNBinsPt];

    fHistoMappingTrueMesonInvMassPtBins                             = new TH1D*[fNBinsPt];    
    fHistoMappingTrueMesonDCInvMassPtBins                           = new TH1D*[fNBinsPt];    
    fHistoMappingTrueFullMesonInvMassPtBins                         = new TH1D*[fNBinsPt];    
    fHistoMappingTrueMesonInvMassPtReweightedBins                   = new TH1D*[fNBinsPt];    
    fHistoMappingTrueMesonInvMassPtUnweightedBins                   = new TH1D*[fNBinsPt];    
    fHistoMappingTrueGGBckInvMassPtBins                             = new TH1D*[fNBinsPt];    
    fHistoMappingTrueContBckInvMassPtBins                           = new TH1D*[fNBinsPt];    
    fHistoMappingTrueAllBckInvMassPtBins                            = new TH1D*[fNBinsPt];    
    fHistoMappingTrueSecMesonInvMassPtBins                          = new TH1D*[fNBinsPt];    
    fHistoMappingTrueSecFromK0SMesonInvMassPtBins                   = new TH1D*[fNBinsPt];    
    fHistoMappingTrueSecFromLambdaMesonInvMassPtBins                = new TH1D*[fNBinsPt];    
    fHistoMappingTrueMesonCaloPhotonInvMassPtBins                   = new TH1D*[fNBinsPt];    
    fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins               = new TH1D*[fNBinsPt];    
    fHistoMappingTrueMesonCaloElectronInvMassPtBins                 = new TH1D*[fNBinsPt];    
    fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins            = new TH1D*[fNBinsPt];    
    fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins    = new TH1D*[fNBinsPt];    
    fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins          = new TH1D*[fNBinsPt];    

    fHistoWeightsBGZbinVsMbin                                       = new TH2F*[fNBinsPt];    
    fHistoFillPerEventBGZbinVsMbin                                  = new TH2F*[fNBinsPt];    
    fHistoWeightsBGZbinVsPsibin                                     = new TH2F*[fNBinsPt];    
    fHistoFillPerEventBGZbinVsPsibin                                = new TH2F*[fNBinsPt];    
    
    fHistoMappingGGInvMassPtBin                                     = new TH1D*[fNBinsPt];    
    fHistoMappingBackInvMassPtBin                                   = new TH1D*[fNBinsPt];
    fHistoMappingBackNormInvMassPtBin                               = new TH1D*[fNBinsPt];
    fHistoMappingSignalInvMassPtBin                                 = new TH1D*[fNBinsPt];
    fHistoMappingSignalRemainingBGSubInvMassPtBin                   = new TH1D*[fNBinsPt];
    fHistoMappingRemainingBGInvMassPtBin                            = new TH1D*[fNBinsPt];
    fHistoMappingRatioSBInvMassPtBin                                = new TH1D*[fNBinsPt];

    fFitSignalInvMassPtBin                                          = new TF1*[fNBinsPt];
    fFitRemainingBGInvMassPtBin                                     = new TF1*[fNBinsPt];
    fFitTrueSignalInvMassPtBin                                      = new TF1*[fNBinsPt];
    fFitTrueSignalInvMassPtReweightedBin                            = new TF1*[fNBinsPt];
    fFitTrueSignalInvMassPtUnweightedBin                            = new TF1*[fNBinsPt];
    fFitTrueSignalCaloConvPhotonInvMassPtBin                        = new TF1*[fNBinsPt];
    fFitTrueSignalCaloElectronInvMassPtBin                          = new TF1*[fNBinsPt];
    fFitTrueSignalMixedCaloConvPhotonInvMassPtBin                   = new TF1*[fNBinsPt];
    fFitTrueSignalCaloMergedClusterInvMassPtBin                     = new TF1*[fNBinsPt];
    fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin             = new TF1*[fNBinsPt];
    fFitTrueSignalCaloPhotonInvMassPtBin                            = new TF1*[fNBinsPt];
    
    fFitSignalPeakPosInvMassPtBin                                   = new TF1*[fNBinsPt];
    fFitBckInvMassPtBin                                             = new TF1*[fNBinsPt];
    fFitRatioInvMassPtBin                                           = new TF1*[fNBinsPt];
    // Histograms for normalization on the left of the peak
    fHistoMappingBackNormInvMassLeftPtBin                           = new TH1D*[fNBinsPt];
    fHistoMappingSignalInvMassLeftPtBin                             = new TH1D*[fNBinsPt];
    fHistoMappingSignalRemainingBGSubInvMassLeftPtBin               = new TH1D*[fNBinsPt];
    fHistoMappingRemainingBGInvMassLeftPtBin                        = new TH1D*[fNBinsPt];
    fHistoMappingPeakPosInvMassPtBin                                = new TH1D*[fNBinsPt];

    fFitInvMassLeftPtBin                                            = new TF1*[fNBinsPt];
    fFitRemainingBGInvMassLeftPtBin                                 = new TF1*[fNBinsPt];
    fFitSignalPeakPosInvMassLeftPtBin                               = new TF1*[fNBinsPt];
    fFitBckInvMassLeftPtBin                                         = new TF1*[fNBinsPt];
    fFitPeakPosPtBin                                                = new TF1*[fNBinsPt];
    fMesonMassPeakPos                                               = new Double_t[fNBinsPt];
    fMesonMassPeakPosError                                          = new Double_t[fNBinsPt];

    fFitSignalGaussianInvMassPtBin                                  = new TF1*[fNBinsPt];
    fFitTrueSignalGaussianInvMassPtBin                              = new TF1*[fNBinsPt];
    fMesonMassGaussian                                              = new Double_t[fNBinsPt];
    fMesonMassGaussianError                                         = new Double_t[fNBinsPt];
    fMesonWidthGaussian                                             = new Double_t[fNBinsPt];
    fMesonWidthGaussianError                                        = new Double_t[fNBinsPt];
    fMesonTrueMassGaussian                                          = new Double_t[fNBinsPt];
    fMesonTrueMassGaussianError                                     = new Double_t[fNBinsPt];
    fMesonTrueWidthGaussian                                         = new Double_t[fNBinsPt];
    fMesonTrueWidthGaussianError                                    = new Double_t[fNBinsPt];
//     if( fMode == 4 ) {
//     fFitSignalInvMassPtBin2                                         = new TF1*[fNBinsPt];
//     fFitSignalPeakPosInvMassPtBin2                                  = new TF1*[fNBinsPt];
//     fFitBckInvMassPtBin2                                            = new TF1*[fNBinsPt];
//     fFitInvMassLeftPtBin2                                           = new TF1*[fNBinsPt];
//     fFitSignalPeakPosInvMassLeftPtBin2                              = new TF1*[fNBinsPt];
//     fFitBckInvMassLeftPtBin2                                        = new TF1*[fNBinsPt];
//     fFitPeakPosPtBin2                                               = new TF1*[fNBinsPt];
//     }
    for(Int_t i = 0;i<fNBinsPt; i++){
//         if( fMode == 4 ) {
//             fFitSignalInvMassPtBin2[i]                                  = NULL;
//             fFitSignalPeakPosInvMassPtBin2[i]                           = NULL;
//             fFitBckInvMassPtBin2[i]                                     = NULL;
//             fFitInvMassLeftPtBin2[i]                                    = NULL;
//             fFitSignalPeakPosInvMassLeftPtBin2[i]                       = NULL;
//             fFitBckInvMassLeftPtBin2[i]                                 = NULL;
//             fFitPeakPosPtBin2[i]                                        = NULL;
//             
//         }
        fHistoMappingTrueMesonInvMassPtBins[i]                              = NULL;
        fHistoMappingTrueMesonDCInvMassPtBins[i]                            = NULL;
        fHistoMappingTrueFullMesonInvMassPtBins[i]                          = NULL;
        fHistoMappingTrueMesonInvMassPtReweightedBins[i]                    = NULL;
        fHistoMappingTrueMesonInvMassPtUnweightedBins[i]                    = NULL;
        fHistoMappingTrueGGBckInvMassPtBins[i]                              = NULL;
        fHistoMappingTrueContBckInvMassPtBins[i]                            = NULL;
        fHistoMappingTrueAllBckInvMassPtBins[i]                             = NULL;
        fHistoMappingTrueSecMesonInvMassPtBins[i]                           = NULL;
        fHistoMappingTrueSecFromK0SMesonInvMassPtBins[i]                    = NULL;
        fHistoMappingTrueSecFromLambdaMesonInvMassPtBins[i]                 = NULL;
        fHistoMappingTrueMesonCaloPhotonInvMassPtBins[i]                    = NULL;    
        fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[i]                = NULL;    
        fHistoMappingTrueMesonCaloElectronInvMassPtBins[i]                  = NULL;    
        fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[i]             = NULL;    
        fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[i]     = NULL;    
        fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[i]           = NULL;    
        
        fHistoMappingGGInvMassPtBin[i]                                      = NULL;
        fHistoMappingBackInvMassPtBin[i]                                    = NULL;
        fHistoMappingBackNormInvMassPtBin[i]                               = NULL;
        fHistoMappingSignalInvMassPtBin[i]                                  = NULL;
        fHistoMappingSignalRemainingBGSubInvMassPtBin[i]                    = NULL;
        fHistoMappingRemainingBGInvMassPtBin[i]                             = NULL;
        fHistoMappingRatioSBInvMassPtBin[i]                                 = NULL;

        fFitSignalInvMassPtBin[i]                                           = NULL;
        fFitRemainingBGInvMassPtBin[i]                                      = NULL;
        fFitTrueSignalInvMassPtBin[i]                                       = NULL;
        fFitTrueSignalInvMassPtReweightedBin[i]                             = NULL;
        fFitTrueSignalInvMassPtUnweightedBin[i]                             = NULL;
        fFitTrueSignalCaloConvPhotonInvMassPtBin[i]                         = NULL;
        fFitTrueSignalCaloElectronInvMassPtBin[i]                           = NULL;
        fFitTrueSignalMixedCaloConvPhotonInvMassPtBin[i]                    = NULL;
        fFitTrueSignalCaloMergedClusterInvMassPtBin[i]                      = NULL;
        fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[i]              = NULL;
        fFitTrueSignalCaloPhotonInvMassPtBin[i]                             = NULL;

        fFitSignalPeakPosInvMassPtBin[i]                                    = NULL;
        fFitBckInvMassPtBin[i]                                              = NULL;
        fFitRatioInvMassPtBin[i]                                            = NULL;
        // Histograms for normalization on the left of the peak
        fHistoMappingBackNormInvMassLeftPtBin[i]                            = NULL;
        fHistoMappingSignalInvMassLeftPtBin[i]                              = NULL;
        fHistoMappingSignalRemainingBGSubInvMassLeftPtBin[i]                = NULL;
        fHistoMappingRemainingBGInvMassLeftPtBin[i]                         = NULL;
        fHistoMappingPeakPosInvMassPtBin[i]                                 = NULL;
        
        fFitInvMassLeftPtBin[i]                                             = NULL;
        fFitRemainingBGInvMassLeftPtBin[i]                                  = NULL;
        fFitSignalPeakPosInvMassLeftPtBin[i]                                = NULL;
        fFitBckInvMassLeftPtBin[i]                                          = NULL;
        fFitPeakPosPtBin[i]                                                 = NULL;
        
        fFitSignalGaussianInvMassPtBin[i]                                   = NULL;
        fFitTrueSignalGaussianInvMassPtBin[i]                               = NULL;
    }
}


//****************************************************************************
//****** Initializiation of MC histogram names according  *****
//****************************************************************************
void SetCorrectMCHistogrammNames(TString mesonType){
    cout << "standard MC chosen" << endl;

    ObjectNameTrue                      = "ESD_TruePrimaryMother_InvMass_Pt";
    ObjectNameTrueFull                  = "ESD_TrueMother_InvMass_Pt";
    ObjectNameTrueWOWeights             = "ESD_TruePrimaryMotherW0Weights_InvMass_Pt";
    ObjectNameProfileWeights            = "ESD_TruePrimaryMotherWeights_InvMass_Pt";
    ObjectNameTrueSec                   = "ESD_TrueSecondaryMother_InvMass_Pt";;
    ObjectNameTrueSecFromK0S            = "ESD_TrueSecondaryMotherFromK0s_InvMass_Pt";
    ObjectNameTrueSecFromLambda         = "ESD_TrueSecondaryMotherFromLambda_InvMass_Pt";
    ObjectNameMCPi0Acc                  = "MC_Pi0InAcc_Pt";
    ObjectNameMCPi0AccWOWeights         = "MC_Pi0WOWeightInAcc_Pt";
    ObjectNameMCEtaAcc                  = "MC_EtaInAcc_Pt";
    ObjectNameMCEtaAccWOWeights         = "MC_EtaWOWeightInAcc_Pt";
    ObjectNameMCPi0                     = "MC_Pi0_Pt";
    ObjectNameMCPi0WOWeights            = "MC_Pi0_WOWeights_Pt";
    ObjectNameMCEta                     = "MC_Eta_Pt";
    ObjectNameMCEtaWOWeights            = "MC_Eta_WOWeights_Pt";
    ObjectNameTrueGGBck                 = "ESD_TrueBckGG_InvMass_Pt";
    ObjectNameTrueContBck               = "ESD_TrueBckCont_InvMass_Pt";
    ObjectNameTrueAllBck                = "ESD_TrueAllCont_InvMass_Pt";
    ObjectNameTrueCaloPhoton            = "ESD_TrueMotherCaloPhoton_InvMass_Pt";
    ObjectNameTrueCaloConvPhoton        = "ESD_TrueMotherCaloConvertedPhoton_InvMass_Pt";
    ObjectNameTrueMixedCaloConvPhoton   = "ESD_TrueMotherCaloMixedPhotonConvertedPhoton_InvMass_Pt";
    ObjectNameTrueCaloElectron          = "ESD_TrueMotherCaloElectron_InvMass_Pt";
    ObjectNameTrueCaloMerged            = "ESD_TrueMotherCaloMergedCluster_InvMass_Pt";
    ObjectNameTrueCaloMergedPartConv    = "ESD_TrueMotherCaloMergedClusterPartConv_InvMass_Pt";
    ObjectNameK0sRecPi0                 = "TrueK0sWithPi0Daughter_MCPt";
    ObjectNameLambdaRecPi0              = "TrueLambdaWithPi0Daughter_MCPt";
    if (fMode > 1 && fMode !=9){
        ObjectNameTrueCaloPhoton            = Form("ESD_True%sCaloPhoton_InvMass_Pt", mesonType.Data());
        ObjectNameTrueCaloConvPhoton        = Form("ESD_True%sCaloConvertedPhoton_InvMass_Pt", mesonType.Data());
        ObjectNameTrueMixedCaloConvPhoton   = Form("ESD_True%sCaloMixedPhotonConvertedPhoton_InvMass_Pt", mesonType.Data());
        if (mesonType.CompareTo("Pi0") == 0){
            if (fMode == 4 || fMode == 5){
                ObjectNameTrueCaloElectron          = Form("ESD_True%srCaloElectron_InvMass_Pt", mesonType.Data());
            } else {
                ObjectNameTrueCaloElectron          = Form("ESD_True%sCaloElectron_InvMass_Pt", mesonType.Data());
            }
        } else {
            ObjectNameTrueCaloElectron      = Form("ESD_True%sCaloElectron_InvMass_Pt", mesonType.Data());
        }
        ObjectNameTrueCaloMerged            = Form("ESD_True%sCaloMergedCluster_InvMass_Pt", mesonType.Data());
        ObjectNameTrueCaloMergedPartConv    = Form("ESD_True%sCaloMergedClusterPartConv_InvMass_Pt", mesonType.Data());
        ObjectNameTrue                      = Form("ESD_TruePrimary%s_InvMass_Pt", mesonType.Data());
        ObjectNameTrueFull                  = Form("ESD_True%s_InvMass_Pt", mesonType.Data());
        ObjectNameTrueWOWeights             = Form("ESD_TruePrimary%sW0Weights_InvMass_Pt", mesonType.Data());
        ObjectNameProfileWeights            = Form("ESD_TruePrimary%sWeights_InvMass_Pt", mesonType.Data());
        ObjectNameTrueSec                   = Form("ESD_TrueSecondary%s_InvMass_Pt", mesonType.Data());
        ObjectNameTrueSecFromK0S            = Form("ESD_TrueSecondary%sFromK0s_InvMass_Pt", mesonType.Data());
        ObjectNameTrueSecFromLambda         = Form("ESD_TrueSecondary%sFromLambda_InvMass_Pt", mesonType.Data());
        
        
    }
    cout <<    ObjectNameTrueCaloPhoton.Data()          << "\t" << 
                ObjectNameTrueCaloConvPhoton.Data()       << "\t" << 
                ObjectNameTrueMixedCaloConvPhoton.Data()    << "\t" << 
                ObjectNameTrueCaloElectron.Data()          << "\t" << 
                ObjectNameTrueCaloMerged.Data()          << "\t" << 
                ObjectNameTrueCaloMergedPartConv.Data() << endl;
    ObjectNameDCMesonInvMassPt          = Form("ESD_TrueDoubleCount%s_InvMass_Pt", mesonType.Data());
    ObjectNameDCGammaClusPt             = "TrueDoubleCountClusterGamma_Pt";
    ObjectNameMesonMultipleCount        = Form("ESD_TrueMultipleCount%s", mesonType.Data());
    ObjectNameGammaClusMultipleCount    = "ESD_TrueMultipleCountClusterGamma";
    
    return;
}


//****************************************************************************
//******************** Projection out of 2D in X *****************************
//****************************************************************************
TH1D* FillProjectionX (TH2* fDummy2D, TString name, Double_t minY, Double_t maxY, Int_t rebin){
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
//*************** Check if histo already exists if not clear *****************
//****************************************************************************
void CheckForNULLForPointer(TH1D* fDummy1){
    if(fDummy1!= NULL){
        delete fDummy1;
        fDummy1           = NULL;
    }
}

//****************************************************************************
//************** Produce background without proper weighting *****************
//****************************************************************************
void ProduceBckWithoutWeighting(TH2D *fBckInvMassVSPtDummy){
    //calculation background for midPt without weighting
    fBckInvMassVSPtDummy->Sumw2();  
    fFittingHistMidPtBackground     = FillProjectionX(fBckInvMassVSPtDummy, "Mapping_Back_InvMass_MidPt", fMidPt[0], fMidPt[1], fNRebin[4]);
    //calulation background for fullPt without weighting
    fMesonFullPtBackground          = FillProjectionX(fBckInvMassVSPtDummy, "Mapping_Back_InvMass_FullPt", fFullPt[0], fFullPt[1], fNRebin[4]);
  
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoBack              = Form("Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingBackInvMassPtBin[iPt]);
        fHistoMappingBackInvMassPtBin[iPt] = FillProjectionX(fBckInvMassVSPtDummy, fNameHistoBack, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
    }
}

//****************************************************************************
//************** Remove BG from Signal ***************************************
//****************************************************************************
void ProcessEM(TH1D* fGammaGamma, TH1D* fBck, Double_t * fBGFitRangeEM) {
    for (Int_t binx= 0; binx < fGammaGamma->GetNbinsX()+1; binx++){
        if(fGammaGamma->GetBinContent(binx) == 0){
            fGammaGamma->SetBinError(binx,1.);
            fGammaGamma->SetBinContent(binx,0.);
        }
    }
    fBckNorm = (TH1D*)fBck->Clone("fBckNorm");
    fGammaGamma->Sumw2();
    fBck->Sumw2();
    fBckNorm->Sumw2();
    
    Double_t    r       = fGammaGamma->Integral(fGammaGamma->GetXaxis()->FindBin(fBGFitRangeEM[0]),fGammaGamma->GetXaxis()->FindBin(fBGFitRangeEM[1]));
    Double_t    b       = fBck->Integral(fBck->GetXaxis()->FindBin(fBGFitRangeEM[0]),fBck->GetXaxis()->FindBin(fBGFitRangeEM[1]));
    Double_t    norm    = 1;
    
    if(b != 0) norm     = r/b;
    fBckNorm->Scale(norm);
    
    Int_t numberOfZeros = 0;
    for (Int_t i = 1; i < fBckNorm->GetNbinsX()+1; i++){
        if (fBckNorm->GetBinContent(i) == 0){
            numberOfZeros++;
            if (norm > 1.){
                fBckNorm->SetBinError(i,1.);
                fBckNorm->SetBinContent(i,0.);
            }
        }
    }
    fSignal             = (TH1D*)fGammaGamma->Clone("fSignal");
    fSignal->Sumw2();
    if ((Double_t)numberOfZeros/fBck->GetNbinsX()< 0.25) fSignal->Add(fBckNorm,-1.);
}


//****************************************************************************
//************** Calculate ratio of signal/ background ***********************
//****************************************************************************
void ProcessRatioSignalBackground(TH1D* fGammaGamma, TH1D* fBck) {
    for (Int_t binx= 0; binx < fGammaGamma->GetNbinsX()+1; binx++){
        if(fGammaGamma->GetBinContent(binx) == 0){
            fGammaGamma->SetBinError(binx,1.);
            fGammaGamma->SetBinContent(binx,0.);
        }
    }
    fGammaGamma->Sumw2();
    fBck->Sumw2();
    fRatioSB = (TH1D*)fGammaGamma->Clone("RatioSB");
    fRatioSB->Divide(fGammaGamma, fBck, 1.,1.,"B");
}


//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//****************************************************************************
void FillMassHistosArray(TH2D* fGammaGammaInvMassVSPtDummy, TH2D *fAlphaPeakPosDummy) {
    fGammaGammaInvMassVSPtDummy->Sumw2();
    fAlphaPeakPosDummy->Sumw2();
    fFittingHistMidPtSignal     = FillProjectionX(fGammaGammaInvMassVSPtDummy, "Mapping_GG_InvMass_MidPt", fMidPt[0], fMidPt[1], fNRebin[4]);
    fMesonFullPtSignal          = FillProjectionX(fGammaGammaInvMassVSPtDummy, "Mapping_GG_InvMass_FullPt", fFullPt[0], fFullPt[1], fNRebin[4]);

    if(fEnergyFlag.CompareTo("PbPb_2.76TeV") != 0 && fEnergyFlag.CompareTo("900GeV") != 0 && fEnergyFlag.CompareTo("2.76TeV") != 0 && fEnergyFlag.CompareTo("pPb_5.023TeV") != 0){
        if(fPrefix.CompareTo("Pi0") == 0){
            for(Int_t iPt=fStartPtBin;iPt<fNBinsPeakPt;iPt++){
                fNameHistoPP    = Form("Mapping_PP_InvMass_in_Pt_Bin%02d", iPt);
                CheckForNULLForPointer(fHistoMappingPeakPosInvMassPtBin[iPt]);
                fHistoMappingPeakPosInvMassPtBin[iPt]=  FillProjectionX(fAlphaPeakPosDummy, fNameHistoPP, fBinsPeakPt[iPt], fBinsPeakPt[iPt+1], fBinsPeakPtRebin[iPt]);
            }
        }
    }
    cout << "nach Peak Pos" << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoGG    = Form("Mapping_GG_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingGGInvMassPtBin[iPt]);
        fHistoMappingGGInvMassPtBin[iPt]=  FillProjectionX(fGammaGammaInvMassVSPtDummy, fNameHistoGG, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated true mesons **********************************************
//****************************************************************************
void FillMassMCTrueMesonHistosArray(TH2D* fHistoTrueMesonInvMassVSPtFill) {
    fHistoTrueMesonInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrue  = Form("Mapping_TrueMeson_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueMesonInvMassPtBins[iPt]);
        fHistoMappingTrueMesonInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueMesonInvMassVSPtFill, fNameHistoTrue, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonInvMassPtBins[iPt]->GetEntries() << endl;
        fHistoMappingTrueMesonInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueMesonInvMassPtBins[iPt]->SetLineColor(2);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated true mesons **********************************************
//****************************************************************************
void FillMassMCTrueFullMesonHistosArray(TH2D* fHistoTrueMesonInvMassVSPtFill) {
    fHistoTrueMesonInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrue  = Form("Mapping_TrueFullMeson_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueFullMesonInvMassPtBins[iPt]);
        fHistoMappingTrueFullMesonInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueMesonInvMassVSPtFill, fNameHistoTrue, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueFullMesonInvMassPtBins[iPt]->GetEntries() << endl;
        fHistoMappingTrueFullMesonInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueFullMesonInvMassPtBins[iPt]->SetLineColor(2);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated true mesons **********************************************
//****************************************************************************
void FillMassMCTrueMesonDCHistosArray(TH2D* fHistoTrueMesonInvMassVSPtFill) {
    fHistoTrueMesonInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrue  = Form("Mapping_TrueMesonDC_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueMesonDCInvMassPtBins[iPt]);
        fHistoMappingTrueMesonDCInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueMesonInvMassVSPtFill, fNameHistoTrue, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonDCInvMassPtBins[iPt]->GetEntries() << endl;
        fHistoMappingTrueMesonDCInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueMesonDCInvMassPtBins[iPt]->SetLineColor(4);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated true mesons, clusters real gammas ************************
//****************************************************************************
void FillMassMCTrueMesonCaloPhotonHistosArray(TH2D* fHistoTrueMesonCaloPhotonInvMassVSPtFill) {
    fHistoTrueMesonCaloPhotonInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrue  = Form("Mapping_TrueMesonCaloPhoton_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]);
        fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueMesonCaloPhotonInvMassVSPtFill, fNameHistoTrue, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]->GetEntries() << endl;
        fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]->SetLineColor(2);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated true mesons, clusters real converted gammas **************
//****************************************************************************
void FillMassMCTrueMesonCaloConvPhotonHistosArray(TH2D* fHistoTrueMesonCaloConvPhotonInvMassVSPtFill) {
    fHistoTrueMesonCaloConvPhotonInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrue  = Form("Mapping_TrueMesonCaloConvPhoton_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]);
        fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueMesonCaloConvPhotonInvMassVSPtFill, fNameHistoTrue, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]->GetEntries() << endl;
        fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]->SetLineColor(2);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated true mesons, clusters real electron **********************
//****************************************************************************
void FillMassMCTrueMesonCaloElectronHistosArray(TH2D* fHistoTrueMesonCaloElectronInvMassVSPtFill) {
    fHistoTrueMesonCaloElectronInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrue  = Form("Mapping_TrueMesonCaloElectron_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt]);
        fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueMesonCaloElectronInvMassVSPtFill, fNameHistoTrue, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt]->GetEntries() << endl;
        fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt]->SetLineColor(2);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//** validated true mesons, mixed clusters real photon, cluster conv photon **
//****************************************************************************
void FillMassMCTrueMesonMixedCaloConvPhotonHistosArray(TH2D* fHistoTrueMesonMixedCaloConvPhotonInvMassVSPtFill) {
    
    fHistoTrueMesonMixedCaloConvPhotonInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrue  = Form("Mapping_TrueMesonMixedCaloConvPhoton_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[iPt]);
        fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueMesonMixedCaloConvPhotonInvMassVSPtFill, fNameHistoTrue, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);

        fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[iPt]->SetLineColor(2);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated true mesons, merged clusters ***************************** 
//****************************************************************************
void FillMassMCTrueMesonCaloMergedClusterHistosArray(TH2D* fHistoTrueMesonCaloMergedClusterInvMassVSPtFill) {
    fHistoTrueMesonCaloMergedClusterInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrue  = Form("Mapping_TrueMesonCaloMergedCluster_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt]);
        fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueMesonCaloMergedClusterInvMassVSPtFill, fNameHistoTrue, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[iPt]->SetLineColor(2);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated true mesons, merged clusters part conv ******************* 
//****************************************************************************
void FillMassMCTrueMesonCaloMergedClusterPartConvHistosArray(TH2D* fHistoTrueMesonCaloMergedClusterPartConvInvMassVSPtFill) {
    fHistoTrueMesonCaloMergedClusterPartConvInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrue  = Form("Mapping_TrueMesonCaloMergedClusterPartConv_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt]);
        fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueMesonCaloMergedClusterPartConvInvMassVSPtFill, fNameHistoTrue, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[iPt]->SetLineColor(2);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated true mesons reweighted ***********************************
//****************************************************************************
void FillMassMCTrueReweightedMesonHistosArray(TH2D* fHistoTrueMesonInvMassVSPtFill) {
    fHistoTrueMesonInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrue  = Form("Mapping_TrueMeson_InvMassReweighted_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]);
        fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]=  FillProjectionX(fHistoTrueMesonInvMassVSPtFill, fNameHistoTrue, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->GetEntries() << endl;
        fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->SetLineColor(2);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated true mesons unweighted ***********************************
//****************************************************************************
void FillMassMCTrueUnweightedMesonHistosArray(TH2D* fHistoTrueMesonInvMassVSPtFill) {
    fHistoTrueMesonInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrue  = Form("Mapping_TrueMeson_InvMassUnweighted_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt]);
        fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt]=  FillProjectionX(fHistoTrueMesonInvMassVSPtFill, fNameHistoTrue, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt]->GetEntries() << endl;
        fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt]->SetLineColor(2);
    }
}


//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated gamma gamma BG *******************************************
//****************************************************************************
void FillMassMCTrueGGBckHistosArray(TH2D* fHistoTrueGGBckInvMassVSPtFill) {
    fHistoTrueGGBckInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrueGGBck     = Form("Mapping_TrueGGBck_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueGGBckInvMassPtBins[iPt]);
        fHistoMappingTrueGGBckInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueGGBckInvMassVSPtFill, fNameHistoTrueGGBck, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        fHistoMappingTrueGGBckInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueGGBckInvMassPtBins[iPt]->SetLineColor(3);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated contamination BG *****************************************
//****************************************************************************
void FillMassMCTrueContBckHistosArray(TH2D* fHistoTrueContBckInvMassVSPtFill) {
    fHistoTrueContBckInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrueContBck   = Form("Mapping_TrueContBck_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueContBckInvMassPtBins[iPt]);
        fHistoMappingTrueContBckInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueContBckInvMassVSPtFill, fNameHistoTrueContBck, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        fHistoMappingTrueContBckInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueContBckInvMassPtBins[iPt]->SetLineColor(5);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated BG *******************************************************
//****************************************************************************
void FillMassMCTrueAllBckHistosArray(TH2D* fHistoTrueAllBckInvMassVSPtFill) {
    fHistoTrueAllBckInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrueAllBck = Form("Mapping_TrueAllBck_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueAllBckInvMassPtBins[iPt]);
        fHistoMappingTrueAllBckInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueAllBckInvMassVSPtFill, fNameHistoTrueAllBck, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        fHistoMappingTrueAllBckInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueAllBckInvMassPtBins[iPt]->SetLineColor(2);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated secondary for mesons from any source *********************
//****************************************************************************
void FillMassMCTrueSecMesonHistosArray(TH2D* fHistoTrueSecMesonInvMassVSPtFill) {
    fHistoTrueSecMesonInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrueSec = Form("Mapping_TrueSecMeson_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueSecMesonInvMassPtBins[iPt]);
        fHistoMappingTrueSecMesonInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueSecMesonInvMassVSPtFill, fNameHistoTrueSec, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        fHistoMappingTrueSecMesonInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueSecMesonInvMassPtBins[iPt]->SetLineColor(2);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated secondary for mesons from K0s ****************************
//****************************************************************************
void FillMassMCTrueSecFromK0SMesonHistosArray(TH2D* fHistoTrueSecFromK0SMesonInvMassVSPtFill) {
    fHistoTrueSecFromK0SMesonInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrueSecFromK0S = Form("Mapping_TrueSecFromK0SMeson_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt]);
        fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueSecFromK0SMesonInvMassVSPtFill, fNameHistoTrueSecFromK0S, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueSecFromK0SMesonInvMassPtBins[iPt]->SetLineColor(2);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated secondary for mesons from Lambdas ************************
//****************************************************************************
void FillMassMCTrueSecFromLambdaMesonHistosArray(TH2D* fHistoTrueSecFromLambdaMesonInvMassVSPtFill) {
    fHistoTrueSecFromLambdaMesonInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrueSecFromLambda = Form("Mapping_TrueSecFromLambdaMeson_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueSecFromLambdaMesonInvMassPtBins[iPt]);
        fHistoMappingTrueSecFromLambdaMesonInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueSecFromLambdaMesonInvMassVSPtFill, fNameHistoTrueSecFromLambda, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        fHistoMappingTrueSecFromLambdaMesonInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueSecFromLambdaMesonInvMassPtBins[iPt]->SetLineColor(2);
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
//*** Create momentum dependent histograms with variable momentum binning ****
//****************************************************************************
void CreatePtHistos(){

    fDeltaPt                            = new TH1D("deltaPt","",fNBinsPt,fBinsPt);
    fDeltaPt->Sumw2();

    fHistoYieldMeson                    = new TH1D("histoYieldMeson","",fNBinsPt,fBinsPt);
    fHistoYieldMeson->Sumw2();
    fHistoYieldTrueMeson                = new TH1D("histoYieldTrueMeson","",fNBinsPt,fBinsPt);
    fHistoYieldTrueMeson->Sumw2();
    fHistoYieldTrueMesonDC              = new TH1D("histoYieldTrueMesonDC","",fNBinsPt,fBinsPt);
    fHistoYieldTrueMesonDC->Sumw2();

    fHistoYieldTrueMesonReweighted      = new TH1D("histoYieldTrueMesonReweighted","",fNBinsPt,fBinsPt);
    fHistoYieldTrueMesonReweighted->Sumw2();
    fHistoYieldTrueMesonUnweighted      = new TH1D("histoYieldTrueMesonUnweighted","",fNBinsPt,fBinsPt);
    fHistoYieldTrueMesonUnweighted->Sumw2();
    fHistoYieldTrueSecMeson             = new TH1D("histoYieldTrueSecMeson","",fNBinsPt,fBinsPt);
    fHistoYieldTrueSecMeson->Sumw2();
    fHistoYieldTrueSecFromK0SMeson      = new TH1D("histoYieldTrueSecFromK0SMeson","",fNBinsPt,fBinsPt);
    fHistoYieldTrueSecFromK0SMeson->Sumw2();
    fHistoYieldTrueSecFromLambdaMeson   = new TH1D("histoYieldTrueSecFromLambdaMeson","",fNBinsPt,fBinsPt);
    fHistoYieldTrueSecFromLambdaMeson->Sumw2();

    fHistoYieldMesonPerEvent            = new TH1D("histoYieldMesonPerEvent","",fNBinsPt,fBinsPt);
    fHistoYieldMesonPerEvent->Sumw2();
    fHistoSignMeson                     = new TH1D("histoSignMeson","",fNBinsPt,fBinsPt);
    fHistoSignMeson->Sumw2();
    fHistoSigndefaultMeson              = new TH1D("histoSigndefaultMeson","",fNBinsPt,fBinsPt);
    fHistoSigndefaultMeson->Sumw2();
    fHistoSBMeson                       = new TH1D("histoSBMeson","",fNBinsPt,fBinsPt);
    fHistoSBMeson->Sumw2();
    fHistoSBdefaultMeson                = new TH1D("histoSBdefaultMeson","",fNBinsPt,fBinsPt);
    fHistoSBdefaultMeson->Sumw2();
    if(fPrefix.CompareTo("Pi0") == 0 && fEnergyFlag.CompareTo("7TeV") == 0){
        fHistoMassPosition              = new TH1D("histoMassPosition","",fNBinsPeakPt,fBinsPeakPtHalf);
        fHistoMassPosition->Sumw2();
    }
    if(fPrefix.CompareTo("Pi0") == 0 && fEnergyFlag.CompareTo("7TeV") == 0){
        fHistoFWHMMesonAlpha01          = new TH1D("histoFWHMMesonAlpha01","",fNBinsPeakPt,fBinsPeakPtHalf); 
        fHistoFWHMMesonAlpha01->Sumw2();
    }
    fHistoMassMeson                     = new TH1D("histoMassMeson","",fNBinsPt,fBinsPt);
    fHistoMassMeson->Sumw2();
    fHistoMassGaussianMeson             = new TH1D("histoMassGaussianMeson","",fNBinsPt,fBinsPt);
    fHistoMassGaussianMeson->Sumw2();
    fHistoWidthGaussianMeson            = new TH1D("histoWidthGaussianMeson","",fNBinsPt,fBinsPt);
    fHistoWidthGaussianMeson->Sumw2();
    fHistoFWHMMeson                     = new TH1D("histoFWHMMeson","",fNBinsPt,fBinsPt);
    fHistoFWHMMeson->Sumw2();
    fHistoTrueMassMeson                 = new TH1D("histoTrueMassMeson","",fNBinsPt,fBinsPt);
    fHistoTrueMassMeson->Sumw2();
    fHistoTrueMassGaussianMeson         = new TH1D("histoTrueMassGaussianMeson","",fNBinsPt,fBinsPt);
    fHistoTrueMassGaussianMeson->Sumw2();
    fHistoTrueMassMesonReweighted       = new TH1D("histoTrueMassMesonReweighted","",fNBinsPt,fBinsPt);
    fHistoTrueMassMesonReweighted->Sumw2();
    fHistoTrueMassMesonUnweighted       = new TH1D("histoTrueMassMesonUnweighted","",fNBinsPt,fBinsPt);
    fHistoTrueMassMesonUnweighted->Sumw2();
    fHistoTrueFWHMMeson                 = new TH1D("histoTrueFWHMMeson","",fNBinsPt,fBinsPt);
    fHistoTrueFWHMMeson->Sumw2();
    fHistoTrueWidthGaussianMeson        = new TH1D("histoTrueWidthGaussianMeson","",fNBinsPt,fBinsPt);
    fHistoTrueWidthGaussianMeson->Sumw2();
    fHistoTrueFWHMMesonReweighted       = new TH1D("histoTrueFWHMMesonReweighted","",fNBinsPt,fBinsPt);
    fHistoTrueFWHMMesonReweighted->Sumw2();
    fHistoTrueFWHMMesonUnweighted       = new TH1D("histoTrueFWHMMesonUnweighted","",fNBinsPt,fBinsPt);
    fHistoTrueFWHMMesonUnweighted->Sumw2();
    fHistoTrueSignMeson                 = new TH1D("histoTrueSignMeson","",fNBinsPt,fBinsPt);
    fHistoTrueSignMeson->Sumw2();
    fHistoTrueSBMeson                   = new TH1D("histoTrueSBMeson","",fNBinsPt,fBinsPt);
    fHistoTrueSBMeson->Sumw2();
    if (fAdvancedMesonQA && (fMode ==2 || fMode == 3 || fMode ==4 || fMode == 5)){
        fHistoTrueMassMesonCaloConvPhoton       = new TH1D("histoTrueMassMesonCaloConvPhoton","",fNBinsPt,fBinsPt);
        fHistoTrueMassMesonCaloConvPhoton->Sumw2();
        fHistoTrueMassMesonCaloElectron         = new TH1D("histoTrueMassMesonCaloConvElectron","",fNBinsPt,fBinsPt);
        fHistoTrueMassMesonCaloElectron->Sumw2();
        fHistoTrueMassMesonCaloMergedCluster    = new TH1D("histoTrueMassMesonCaloMergedCluster","",fNBinsPt,fBinsPt);
        fHistoTrueMassMesonCaloMergedCluster->Sumw2();
        fHistoTrueMassMesonCaloMergedPartConvCluster    = new TH1D("histoTrueMassMesonCaloMergedPartConvCluster","",fNBinsPt,fBinsPt);
        fHistoTrueMassMesonCaloMergedPartConvCluster->Sumw2();
        fHistoTrueMassMesonCaloPhoton                   = new TH1D("histoTrueMassMesonCaloPhoton","",fNBinsPt,fBinsPt);
        fHistoTrueMassMesonCaloPhoton->Sumw2();
        fHistoTrueFWHMMesonCaloConvPhoton               = new TH1D("histoTrueFWHMMesonCaloConvPhoton","",fNBinsPt,fBinsPt);
        fHistoTrueFWHMMesonCaloConvPhoton->Sumw2();
        fHistoTrueFWHMMesonCaloElectron                 = new TH1D("histoTrueFWHMMesonCaloConvElectron","",fNBinsPt,fBinsPt);
        fHistoTrueFWHMMesonCaloElectron->Sumw2();
        fHistoTrueFWHMMesonCaloMergedCluster            = new TH1D("histoTrueFWHMMesonCaloMergedCluster","",fNBinsPt,fBinsPt);
        fHistoTrueFWHMMesonCaloMergedCluster->Sumw2();
        fHistoTrueFWHMMesonCaloMergedPartConvCluster    = new TH1D("histoTrueFWHMMesonCaloMergedPartConvCluster","",fNBinsPt,fBinsPt);
        fHistoTrueFWHMMesonCaloMergedPartConvCluster->Sumw2();
        fHistoTrueFWHMMesonCaloPhoton                   = new TH1D("histoTrueFWHMMesonCaloPhoton","",fNBinsPt,fBinsPt);
        fHistoTrueFWHMMesonCaloPhoton->Sumw2();
        fHistoYieldTrueMesonFixedWindow                 = new TH1D("histoYieldTrueMesonFixedWindow","",fNBinsPt,fBinsPt);
        fHistoYieldTrueMesonFixedWindow->Sumw2();
        fHistoYieldTrueMesonGammaFixedWindow            = new TH1D("histoYieldTrueMesonGammaFixedWindow","",fNBinsPt,fBinsPt);
        fHistoYieldTrueMesonGammaFixedWindow->Sumw2();
        fHistoYieldTrueMesonGammaConvGammaFixedWindow   = new TH1D("histoYieldTrueMesonGammaConvGammaFixedWindow","",fNBinsPt,fBinsPt);
        fHistoYieldTrueMesonGammaConvGammaFixedWindow->Sumw2();
        fHistoYieldTrueMesonConvGammaConvGammaFixedWindow=   new TH1D("histoYieldTrueMesonConvGammaConvGammaFixedWindow","",fNBinsPt,fBinsPt);
        fHistoYieldTrueMesonConvGammaConvGammaFixedWindow->Sumw2();
    }
    if (fAdvancedMesonQA && (fMode ==4 || fMode == 5)){
        fHistoTrueMassMesonMixedCaloConvPhoton =       new TH1D("histoTrueMassMesonMixedCaloConvPhoton","",fNBinsPt,fBinsPt);
        fHistoTrueMassMesonMixedCaloConvPhoton->Sumw2();
        fHistoTrueFWHMMesonMixedCaloConvPhoton =       new TH1D("histoTrueFWHMMesonMixedCaloConvPhoton","",fNBinsPt,fBinsPt);
        fHistoTrueFWHMMesonMixedCaloConvPhoton->Sumw2();
    }

    fHistoLambdaTail                    = new TH1D("histoLambdaTail","",fNBinsPt,fBinsPt);
    fHistoLambdaTail->Sumw2();

    fHistoMassWindowHigh                = new TH1D("histoMassWindowHigh","",fNBinsPt,fBinsPt);
    fHistoMassWindowHigh->Sumw2();
    fHistoMassWindowLow                 = new TH1D("histoMassWindowLow","",fNBinsPt,fBinsPt);
    fHistoMassWindowLow->Sumw2();
    fHistoMassWindowWideHigh            = new TH1D("histoMassWindowWideHigh","",fNBinsPt,fBinsPt);
    fHistoMassWindowWideHigh->Sumw2();
    fHistoMassWindowWideLow             = new TH1D("histoMassWindowWideLow","",fNBinsPt,fBinsPt);
    fHistoMassWindowWideLow->Sumw2();
    fHistoMassWindowNarrowHigh          = new TH1D("histoMassWindowNarrowHigh","",fNBinsPt,fBinsPt);
    fHistoMassWindowNarrowHigh->Sumw2();
    fHistoMassWindowNarrowLow           = new TH1D("histoMassWindowNarrowLow","",fNBinsPt,fBinsPt);
    fHistoMassWindowNarrowLow->Sumw2();
    
    fHistoYieldMesonNarrow              = new TH1D("histoYieldMesonNarrow","",fNBinsPt,fBinsPt);
    fHistoYieldMesonNarrow->Sumw2();
    fHistoYieldTrueMesonNarrow          = new TH1D("histoYieldTrueMesonNarrow","",fNBinsPt,fBinsPt);
    fHistoYieldTrueMesonNarrow->Sumw2();
    fHistoYieldTrueMesonReweightedNarrow    = new TH1D("histoYieldTrueMesonNarrowReweighted","",fNBinsPt,fBinsPt);
    fHistoYieldTrueMesonReweightedNarrow->Sumw2();
    fHistoYieldTrueMesonUnweightedNarrow    = new TH1D("histoYieldTrueMesonNarrowUnweighted","",fNBinsPt,fBinsPt);
    fHistoYieldTrueMesonUnweightedNarrow->Sumw2();
    fHistoYieldTrueSecMesonNarrow           = new TH1D("histoYieldTrueSecMesonNarrow","",fNBinsPt,fBinsPt);
    fHistoYieldTrueSecMesonNarrow->Sumw2();
    fHistoYieldTrueSecFromK0SMesonNarrow    = new TH1D("histoYieldTrueSecFromK0SMesonNarrow","",fNBinsPt,fBinsPt);
    fHistoYieldTrueSecFromK0SMesonNarrow->Sumw2();
    fHistoYieldTrueSecFromLambdaMesonNarrow = new TH1D("histoYieldTrueSecFromLambdaMesonNarrow","",fNBinsPt,fBinsPt);
    fHistoYieldTrueSecFromLambdaMesonNarrow->Sumw2();

    fHistoYieldMesonPerEventNarrow          = new TH1D("histoYieldMesonPerEventNarrow","",fNBinsPt,fBinsPt);
    fHistoYieldMesonPerEventNarrow->Sumw2();
    fHistoSignMesonNarrow                   = new TH1D("histoSignMesonNarrow","",fNBinsPt,fBinsPt);
    fHistoSignMesonNarrow->Sumw2();
    fHistoSBMesonNarrow                     = new TH1D("histoSBMesonNarrow","",fNBinsPt,fBinsPt);
    fHistoSBMesonNarrow->Sumw2();
    fHistoSBdefaultNarrowMeson              = new TH1D("histoSBdefaultNarrowMeson","",fNBinsPt,fBinsPt);
    fHistoSBdefaultNarrowMeson->Sumw2();
    fHistoSigndefaultNarrowMeson            = new TH1D("histoSigndefaultNarrowMeson","",fNBinsPt,fBinsPt);
    fHistoSigndefaultNarrowMeson->Sumw2();

    fHistoYieldMesonWide                    = new TH1D("histoYieldMesonWide","",fNBinsPt,fBinsPt);
    fHistoYieldMesonWide->Sumw2();
    fHistoYieldTrueMesonWide                = new TH1D("histoYieldTrueMesonWide","",fNBinsPt,fBinsPt);
    fHistoYieldTrueMesonWide->Sumw2();
    fHistoYieldTrueMesonReweightedWide      = new TH1D("histoYieldTrueMesonWideReweighted","",fNBinsPt,fBinsPt);
    fHistoYieldTrueMesonReweightedWide->Sumw2();
    fHistoYieldTrueMesonUnweightedWide      = new TH1D("histoYieldTrueMesonWideUnweighted","",fNBinsPt,fBinsPt);
    fHistoYieldTrueMesonUnweightedWide->Sumw2();
    fHistoYieldTrueSecMesonWide             = new TH1D("histoYieldTrueSecMesonWide","",fNBinsPt,fBinsPt);
    fHistoYieldTrueSecMesonWide->Sumw2();
    fHistoYieldTrueSecFromK0SMesonWide      = new TH1D("histoYieldTrueSecFromK0SMesonWide","",fNBinsPt,fBinsPt);
    fHistoYieldTrueSecFromK0SMesonWide->Sumw2();
    fHistoYieldTrueSecFromLambdaMesonWide   = new TH1D("histoYieldTrueSecFromLambdaMesonWide","",fNBinsPt,fBinsPt);
    fHistoYieldTrueSecFromLambdaMesonWide->Sumw2();
    fHistoYieldMesonPerEventWide            = new TH1D("histoYieldMesonPerEventWide","",fNBinsPt,fBinsPt);
    fHistoYieldMesonPerEventWide->Sumw2();
    fHistoSignMesonWide                     = new TH1D("histoSignMesonWide","",fNBinsPt,fBinsPt);
    fHistoSignMesonWide->Sumw2();
    fHistoSBMesonWide                       = new TH1D("histoSBMesonWide","",fNBinsPt,fBinsPt);
    fHistoSBMesonWide->Sumw2();

    // Histos for normalization at the left of the peak

    fHistoYieldMesonLeft                    = new TH1D("histoYieldMesonLeft","",fNBinsPt,fBinsPt);
    fHistoYieldMesonLeft->Sumw2();
    fHistoYieldMesonLeftPerEvent            = new TH1D("histoYieldMesonLeftPerEvent","",fNBinsPt,fBinsPt);
    fHistoYieldMesonLeftPerEvent->Sumw2();
    fHistoSignMesonLeft                     = new TH1D("histoSignMesonLeft","",fNBinsPt,fBinsPt);
    fHistoSignMesonLeft->Sumw2();
    fHistoSBMesonLeft                       = new TH1D("histoSBMesonLeft","",fNBinsPt,fBinsPt);
    fHistoSBMesonLeft->Sumw2();
    fHistoMassMesonLeft                     = new TH1D("histoMassMesonLeft","",fNBinsPt,fBinsPt);
    fHistoMassMesonLeft->Sumw2();
    fHistoFWHMMesonLeft                     = new TH1D("histoFWHMMesonLeft","",fNBinsPt,fBinsPt);
    fHistoFWHMMesonLeft->Sumw2();


    fHistoYieldMesonLeftNarrow              = new TH1D("histoYieldMesonLeftNarrow","",fNBinsPt,fBinsPt);
    fHistoYieldMesonLeftNarrow->Sumw2();
    fHistoYieldMesonLeftPerEventNarrow      = new TH1D("histoYieldMesonLeftPerEventNarrow","",fNBinsPt,fBinsPt);
    fHistoYieldMesonLeftPerEventNarrow->Sumw2();
    fHistoSignMesonLeftNarrow               = new TH1D("histoSignMesonLeftNarrow","",fNBinsPt,fBinsPt);
    fHistoSignMesonLeftNarrow->Sumw2();
    fHistoSBMesonLeftNarrow                 = new TH1D("histoSBMesonLeftNarrow","",fNBinsPt,fBinsPt);
    fHistoSBMesonLeftNarrow->Sumw2();


    fHistoYieldMesonLeftWide                = new TH1D("histoYieldMesonLeftWide","",fNBinsPt,fBinsPt);
    fHistoYieldMesonLeftWide->Sumw2();
    fHistoYieldMesonLeftPerEventWide        = new TH1D("histoYieldMesonLeftPerEventWide","",fNBinsPt,fBinsPt);
    fHistoYieldMesonLeftPerEventWide->Sumw2();
    fHistoSignMesonLeftWide                 = new TH1D("histoSignMesonLeftWide","",fNBinsPt,fBinsPt);
    fHistoSignMesonLeftWide->Sumw2();
    fHistoSBMesonLeftWide                   = new TH1D("histoSBMesonLeftWide","",fNBinsPt,fBinsPt);
    fHistoSBMesonLeftWide->Sumw2();

}


//****************************************************************************
//*************** Fill momentum dependent histograms from arrays *************
//****************************************************************************
void FillPtHistos()
{
    if(fPrefix.CompareTo("Pi0") == 0 && fEnergyFlag.CompareTo("7TeV") == 0){
        for (Int_t iPt=fStartPtBin+1;iPt<fNBinsPeakPt+1;iPt++){
            fHistoMassPosition->SetBinContent(iPt,fMesonMassPeakPos[iPt-1]);
            fHistoMassPosition->SetBinError(iPt,fMesonMassPeakPosError[iPt-1]);
            fHistoFWHMMesonAlpha01->SetBinContent(iPt,fMesonFWHMAlpha01[iPt-1]);
            fHistoFWHMMesonAlpha01->SetBinError(iPt,fMesonFWHMAlpha01Error[iPt-1]);
        }
    }
    for(Int_t iPt=fStartPtBin+1;iPt<fNBinsPt+1;iPt++){

        fDeltaPt->SetBinContent(iPt,fBinsPt[iPt]-fBinsPt[iPt-1]);
        fDeltaPt->SetBinError(iPt,0);


        fHistoMassMeson->SetBinContent(iPt,fMesonMass[iPt-1]);
        fHistoMassMeson->SetBinError(iPt,fMesonMassError[iPt-1]);
        fHistoMassGaussianMeson->SetBinContent(iPt,fMesonMassGaussian[iPt-1]);
        fHistoMassGaussianMeson->SetBinError(iPt,fMesonMassGaussianError[iPt-1]);
        fHistoWidthGaussianMeson->SetBinContent(iPt,fMesonWidthGaussian[iPt-1]);
        fHistoWidthGaussianMeson->SetBinError(iPt,fMesonWidthGaussianError[iPt-1]);
        fHistoFWHMMeson->SetBinContent(iPt,fMesonFWHM[iPt-1]);
        fHistoFWHMMeson->SetBinError(iPt,fMesonFWHMError[iPt-1]);

        if (fIsMC) {
            fHistoTrueMassMeson->SetBinContent(iPt,fMesonTrueMass[iPt-1]);
            fHistoTrueMassMeson->SetBinError(iPt,fMesonTrueMassError[iPt-1]);
            fHistoTrueMassGaussianMeson->SetBinContent(iPt,fMesonTrueMassGaussian[iPt-1]);
            fHistoTrueMassGaussianMeson->SetBinError(iPt,fMesonTrueMassGaussianError[iPt-1]);

            fHistoTrueMassMesonReweighted->SetBinContent(iPt,fMesonTrueMassReweighted[iPt-1]);
            fHistoTrueMassMesonReweighted->SetBinError(iPt,fMesonTrueMassReweightedError[iPt-1]);
            fHistoTrueMassMesonUnweighted->SetBinContent(iPt,fMesonTrueMassUnweighted[iPt-1]);
            fHistoTrueMassMesonUnweighted->SetBinError(iPt,fMesonTrueMassUnweightedError[iPt-1]);

            fHistoTrueFWHMMeson->SetBinContent(iPt,fMesonTrueFWHM[iPt-1]);
            fHistoTrueFWHMMeson->SetBinError(iPt,fMesonTrueFWHMError[iPt-1]);
            fHistoTrueWidthGaussianMeson->SetBinContent(iPt,fMesonTrueWidthGaussian[iPt-1]);
            fHistoTrueWidthGaussianMeson->SetBinError(iPt,fMesonTrueWidthGaussianError[iPt-1]);
            fHistoTrueFWHMMesonReweighted->SetBinContent(iPt,fMesonTrueFWHMReweighted[iPt-1]);
            fHistoTrueFWHMMesonReweighted->SetBinError(iPt,fMesonTrueFWHMReweightedError[iPt-1]);
            fHistoTrueFWHMMesonUnweighted->SetBinContent(iPt,fMesonTrueFWHMUnweighted[iPt-1]);
            fHistoTrueFWHMMesonUnweighted->SetBinError(iPt,fMesonTrueFWHMUnweightedError[iPt-1]);
            
            fHistoTrueSignMeson->SetBinContent(iPt,fMesonTrueSign[iPt-1]);
            fHistoTrueSignMeson->SetBinError(iPt,fMesonTrueSignError[iPt-1]);
            fHistoTrueSBMeson->SetBinContent(iPt,fMesonTrueSB[iPt-1]);
            fHistoTrueSBMeson->SetBinError(iPt,fMesonTrueSBError[iPt-1]);
            
            if (fAdvancedMesonQA && (fMode == 2 || fMode == 3 || fMode ==4 || fMode == 5)){
                fHistoTrueMassMesonCaloConvPhoton->SetBinContent(iPt,fMesonTrueMassCaloConvPhoton[iPt-1]);
                fHistoTrueMassMesonCaloConvPhoton->SetBinError(iPt,fMesonTrueMassErrorCaloConvPhoton[iPt-1]);
                fHistoTrueMassMesonCaloElectron->SetBinContent(iPt,fMesonTrueMassCaloElectron[iPt-1]);
                fHistoTrueMassMesonCaloElectron->SetBinError(iPt,fMesonTrueMassErrorCaloElectron[iPt-1]);
                fHistoTrueMassMesonCaloMergedCluster->SetBinContent(iPt,fMesonTrueMassCaloMergedCluster[iPt-1]);
                fHistoTrueMassMesonCaloMergedCluster->SetBinError(iPt,fMesonTrueMassErrorCaloMergedCluster[iPt-1]);
                fHistoTrueMassMesonCaloMergedPartConvCluster->SetBinContent(iPt,fMesonTrueMassCaloMergedClusterPartConv[iPt-1]);
                fHistoTrueMassMesonCaloMergedPartConvCluster->SetBinError(iPt,fMesonTrueMassErrorCaloMergedClusterPartConv[iPt-1]);
                fHistoTrueMassMesonCaloPhoton->SetBinContent(iPt,fMesonTrueMassCaloPhoton[iPt-1]);
                fHistoTrueMassMesonCaloPhoton->SetBinError(iPt,fMesonTrueMassErrorCaloPhoton[iPt-1]);
                fHistoTrueFWHMMesonCaloConvPhoton->SetBinContent(iPt,fMesonTrueFWHMCaloConvPhoton[iPt-1]);
                fHistoTrueFWHMMesonCaloConvPhoton->SetBinError(iPt,fMesonTrueFWHMErrorCaloConvPhoton[iPt-1]);
                fHistoTrueFWHMMesonCaloElectron->SetBinContent(iPt,fMesonTrueFWHMCaloElectron[iPt-1]);
                fHistoTrueFWHMMesonCaloElectron->SetBinError(iPt,fMesonTrueFWHMErrorCaloElectron[iPt-1]);
                fHistoTrueFWHMMesonCaloMergedCluster->SetBinContent(iPt,fMesonTrueFWHMCaloMergedCluster[iPt-1]);
                fHistoTrueFWHMMesonCaloMergedCluster->SetBinError(iPt,fMesonTrueFWHMErrorCaloMergedCluster[iPt-1]);
                fHistoTrueFWHMMesonCaloMergedPartConvCluster->SetBinContent(iPt,fMesonTrueFWHMCaloMergedClusterPartConv[iPt-1]);
                fHistoTrueFWHMMesonCaloMergedPartConvCluster->SetBinError(iPt,fMesonTrueFWHMErrorCaloMergedClusterPartConv[iPt-1]);
                fHistoTrueFWHMMesonCaloPhoton->SetBinContent(iPt,fMesonTrueFWHMCaloPhoton[iPt-1]);
                fHistoTrueFWHMMesonCaloPhoton->SetBinError(iPt,fMesonTrueFWHMErrorCaloPhoton[iPt-1]);

                fHistoYieldTrueMesonFixedWindow->SetBinContent(iPt,fMesonTrueYieldFixedWindow[iPt-1]);
                fHistoYieldTrueMesonFixedWindow->SetBinError(iPt,fMesonTrueYieldErrorFixedWindow[iPt-1]);
                fHistoYieldTrueMesonGammaFixedWindow->SetBinContent(iPt,fMesonTrueYieldGammaFixedWindow[iPt-1]);
                fHistoYieldTrueMesonGammaFixedWindow->SetBinError(iPt,fMesonTrueYieldGammaErrorFixedWindow[iPt-1]);
                fHistoYieldTrueMesonGammaConvGammaFixedWindow->SetBinContent(iPt,fMesonTrueYieldGammaConvGammaFixedWindow[iPt-1]);
                fHistoYieldTrueMesonGammaConvGammaFixedWindow->SetBinError(iPt,fMesonTrueYieldGammaConvGammaErrorFixedWindow[iPt-1]);
                fHistoYieldTrueMesonConvGammaConvGammaFixedWindow->SetBinContent(iPt,fMesonTrueYieldConvGammaConvGammaFixedWindow[iPt-1]);
                fHistoYieldTrueMesonConvGammaConvGammaFixedWindow->SetBinError(iPt,fMesonTrueYieldConvGammaConvGammaErrorFixedWindow[iPt-1]);
            } 
            if (fAdvancedMesonQA && ( fMode ==4 || fMode == 5)){
                fHistoTrueFWHMMesonMixedCaloConvPhoton->SetBinContent(iPt,fMesonTrueFWHMMixedCaloConvPhoton[iPt-1]);
                fHistoTrueFWHMMesonMixedCaloConvPhoton->SetBinError(iPt,fMesonTrueFWHMErrorMixedCaloConvPhoton[iPt-1]);
                fHistoTrueMassMesonMixedCaloConvPhoton->SetBinContent(iPt,fMesonTrueMassMixedCaloConvPhoton[iPt-1]);
                fHistoTrueMassMesonMixedCaloConvPhoton->SetBinError(iPt,fMesonTrueMassErrorMixedCaloConvPhoton[iPt-1]);
            }
        }
        
        if (fIsMC) {
            fHistoLambdaTail->SetBinContent(iPt,fMesonLambdaTailMCpar[iPt-1]);
            fHistoLambdaTail->SetBinError(iPt,fMesonLambdaTailMCparError[iPt-1]);
        }
        fHistoLambdaTail->SetBinContent(iPt,fMesonLambdaTailpar[iPt-1]);
        fHistoLambdaTail->SetBinError(iPt,fMesonLambdaTailparError[iPt-1]);

        fHistoSignMeson->SetBinContent(iPt,fMesonSign[iPt-1]);
        fHistoSignMeson->SetBinError(iPt,fMesonSignError[iPt-1]);
        fHistoSBMeson->SetBinContent(iPt,fMesonSB[iPt-1]);
        fHistoSBMeson->SetBinError(iPt,fMesonSBError[iPt-1]);
        fHistoSBdefaultMeson->SetBinContent(iPt,fMesonSBdefault[iPt-1]);
        fHistoSBdefaultMeson->SetBinError(iPt,fMesonSBdefaultError[iPt-1]);
        fHistoSigndefaultMeson->SetBinContent(iPt,fMesonSigndefault[iPt-1]);
        fHistoSigndefaultMeson->SetBinError(iPt,fMesonSigndefaultError[iPt-1]);
        
        fHistoMassWindowHigh->SetBinContent(iPt,fMassWindowHigh[iPt-1]);
        fHistoMassWindowLow->SetBinContent(iPt,fMassWindowLow[iPt-1]);
        fHistoMassWindowWideHigh->SetBinContent(iPt,fMassWindowWideHigh[iPt-1]);
        fHistoMassWindowWideLow->SetBinContent(iPt,fMassWindowWideLow[iPt-1]);
        fHistoMassWindowNarrowHigh->SetBinContent(iPt,fMassWindowNarrowHigh[iPt-1]);
        fHistoMassWindowNarrowLow->SetBinContent(iPt,fMassWindowNarrowLow[iPt-1]);

        
        fHistoYieldMeson->SetBinContent(iPt,fMesonYieldsCorResidualBckFunc[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMeson->SetBinError(iPt,fMesonYieldsCorResidualBckFuncError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        
        if (fIsMC) {
            fHistoYieldTrueMeson->SetBinContent(iPt,fMesonTrueYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueMeson->SetBinError(iPt,fMesonTrueYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            if (fEnableDCMeson) fHistoYieldTrueMesonDC->SetBinContent(iPt,fMesonTrueYieldsDC[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            if (fEnableDCMeson) fHistoYieldTrueMesonDC->SetBinError(iPt,fMesonTrueYieldsDCError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueMesonReweighted->SetBinContent(iPt,fMesonTrueYieldsReweighted[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueMesonReweighted->SetBinError(iPt,fMesonTrueYieldsReweightedError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueMesonUnweighted->SetBinContent(iPt,fMesonTrueYieldsUnweighted[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueMesonUnweighted->SetBinError(iPt,fMesonTrueYieldsUnweightedError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecMeson->SetBinContent(iPt,fMesonTrueSecYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecMeson->SetBinError(iPt,fMesonTrueSecYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecFromK0SMeson->SetBinContent(iPt,fMesonTrueSecFromK0SYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecFromK0SMeson->SetBinError(iPt,fMesonTrueSecFromK0SYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecFromLambdaMeson->SetBinContent(iPt,fMesonTrueSecFromLambdaYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecFromLambdaMeson->SetBinError(iPt,fMesonTrueSecFromLambdaYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        
        }
        fHistoYieldMesonPerEvent->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFunc[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonPerEvent->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

        // Narrow integration window
        fHistoSignMesonNarrow->SetBinContent(iPt,fMesonSignNarrow[iPt-1]);
        fHistoSignMesonNarrow->SetBinError(iPt,fMesonSignNarrowError[iPt-1]);
        fHistoSBMesonNarrow->SetBinContent(iPt,fMesonSBNarrow[iPt-1]);
        fHistoSBMesonNarrow->SetBinError(iPt,fMesonSBNarrowError[iPt-1]);
        
        fHistoSBdefaultNarrowMeson->SetBinContent(iPt,fMesonSBdefaultNarrow[iPt-1]);
        fHistoSBdefaultNarrowMeson->SetBinError(iPt,fMesonSBdefaultNarrowError[iPt-1]);
        fHistoSigndefaultNarrowMeson->SetBinContent(iPt,fMesonSigndefaultNarrow[iPt-1]);
        fHistoSigndefaultNarrowMeson->SetBinError(iPt,fMesonSigndefaultNarrowError[iPt-1]);


        fHistoYieldMesonNarrow->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonNarrow->SetBinError(iPt,fMesonYieldsCorResidualBckFuncNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

        if (fIsMC) {
            fHistoYieldTrueMesonNarrow->SetBinContent(iPt,fMesonTrueYieldsNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueMesonNarrow->SetBinError(iPt,fMesonTrueYieldsNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueMesonReweightedNarrow->SetBinContent(iPt,fMesonTrueYieldsReweightedNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueMesonReweightedNarrow->SetBinError(iPt,fMesonTrueYieldsReweightedNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueMesonUnweightedNarrow->SetBinContent(iPt,fMesonTrueYieldsUnweightedNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueMesonUnweightedNarrow->SetBinError(iPt,fMesonTrueYieldsUnweightedNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecMesonNarrow->SetBinContent(iPt,fMesonTrueSecYieldsNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecMesonNarrow->SetBinError(iPt,fMesonTrueSecYieldsNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecFromK0SMesonNarrow->SetBinContent(iPt,fMesonTrueSecFromK0SYieldsNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecFromK0SMesonNarrow->SetBinError(iPt,fMesonTrueSecFromK0SYieldsNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecFromLambdaMesonNarrow->SetBinContent(iPt,fMesonTrueSecFromLambdaYieldsNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecFromLambdaMesonNarrow->SetBinError(iPt,fMesonTrueSecFromLambdaYieldsNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        
        }
        fHistoYieldMesonPerEventNarrow->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonPerEventNarrow->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

        // Wide integration window
        fHistoSignMesonWide->SetBinContent(iPt,fMesonSignWide[iPt-1]);
        fHistoSignMesonWide->SetBinError(iPt,fMesonSignWideError[iPt-1]);
        fHistoSBMesonWide->SetBinContent(iPt,fMesonSBWide[iPt-1]);
        fHistoSBMesonWide->SetBinError(iPt,fMesonSBWideError[iPt-1]);

        fHistoYieldMesonWide->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonWide->SetBinError(iPt,fMesonYieldsCorResidualBckFuncWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

        if (fIsMC) {
            fHistoYieldTrueMesonWide->SetBinContent(iPt,fMesonTrueYieldsWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueMesonWide->SetBinError(iPt,fMesonTrueYieldsWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueMesonReweightedWide->SetBinContent(iPt,fMesonTrueYieldsReweightedWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueMesonReweightedWide->SetBinError(iPt,fMesonTrueYieldsReweightedWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueMesonUnweightedWide->SetBinContent(iPt,fMesonTrueYieldsUnweightedWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueMesonUnweightedWide->SetBinError(iPt,fMesonTrueYieldsUnweightedWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecMesonWide->SetBinContent(iPt,fMesonTrueSecYieldsWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecMesonWide->SetBinError(iPt,fMesonTrueSecYieldsWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecFromK0SMesonWide->SetBinContent(iPt,fMesonTrueSecFromK0SYieldsWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecFromK0SMesonWide->SetBinError(iPt,fMesonTrueSecFromK0SYieldsWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecFromLambdaMesonWide->SetBinContent(iPt,fMesonTrueSecFromLambdaYieldsWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldTrueSecFromLambdaMesonWide->SetBinError(iPt,fMesonTrueSecFromLambdaYieldsWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        }
        fHistoYieldMesonPerEventWide->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonPerEventWide->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

        // Histos for integration at the left of the peak
        fHistoMassMesonLeft->SetBinContent(iPt,fMesonMassLeft[iPt-1]);
        fHistoMassMesonLeft->SetBinError(iPt,fMesonMassLeftError[iPt-1]);
        fHistoFWHMMesonLeft->SetBinContent(iPt,fMesonFWHMLeft[iPt-1]);
        fHistoFWHMMesonLeft->SetBinError(iPt,fMesonFWHMLeftError[iPt-1]);

        fHistoSignMesonLeft->SetBinContent(iPt,fMesonSignLeft[iPt-1]);
        fHistoSignMesonLeft->SetBinError(iPt,fMesonSignLeftError[iPt-1]);
        fHistoSBMesonLeft->SetBinContent(iPt,fMesonSBLeft[iPt-1]);
        fHistoSBMesonLeft->SetBinError(iPt,fMesonSBLeftError[iPt-1]);

        fHistoYieldMesonLeft->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncLeft[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonLeft->SetBinError(iPt,fMesonYieldsCorResidualBckFuncLeftError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonLeftPerEvent->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeft[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonLeftPerEvent->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

        // Narrow integration window
        fHistoSignMesonLeftNarrow->SetBinContent(iPt,fMesonSignLeftNarrow[iPt-1]);
        fHistoSignMesonLeftNarrow->SetBinError(iPt,fMesonSignLeftNarrowError[iPt-1]);
        fHistoSBMesonLeftNarrow->SetBinContent(iPt,fMesonSBLeftNarrow[iPt-1]);
        fHistoSBMesonLeftNarrow->SetBinError(iPt,fMesonSBLeftNarrowError[iPt-1]);

        fHistoYieldMesonLeftNarrow->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncLeftNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonLeftNarrow->SetBinError(iPt,fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonLeftPerEventNarrow->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftNarrow[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonLeftPerEventNarrow->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftNarrowError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

        // Wide integration window
        fHistoSignMesonLeftWide->SetBinContent(iPt,fMesonSignLeftWide[iPt-1]);
        fHistoSignMesonLeftWide->SetBinError(iPt,fMesonSignLeftWideError[iPt-1]);
        fHistoSBMesonLeftWide->SetBinContent(iPt,fMesonSBLeftWide[iPt-1]);
        fHistoSBMesonLeftWide->SetBinError(iPt,fMesonSBLeftWideError[iPt-1]);

        fHistoYieldMesonLeftWide->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncLeftWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonLeftWide->SetBinError(iPt,fMesonYieldsCorResidualBckFuncLeftWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonLeftPerEventWide->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftWide[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonLeftPerEventWide->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncLeftWideError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
    }
}

//****************************************************************************
//******** Fit of Signal+ BG with Gaussian + Exponential + Linear BG *********
//****************************************************************************
void FitSubtractedInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Double_t* fMesonIntDeltaRangeFit, Int_t ptBin, Bool_t vary){

    //    cout<<"Start Fitting spectra"<<endl;
    fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude = fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin;
    Double_t mesonAmplitudeMax;
    
    if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0){
        if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
            if(ptBin == 1) mesonAmplitudeMin = mesonAmplitude*70./100.;
            else if(ptBin > 1 && ptBin < 4) mesonAmplitudeMin = mesonAmplitude*90./100.;
            else if(ptBin > 17) mesonAmplitudeMin = mesonAmplitude*80./100.;
            else mesonAmplitudeMin = mesonAmplitude*98./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            if (fMode == 2 || fMode == 3) {
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*600./100.;
            }
            if (fMode == 4 || fMode == 5) {
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                mesonAmplitudeMax = mesonAmplitude*400./100.;
            }

        } else {
            if(ptBin == 2) mesonAmplitudeMin = mesonAmplitude*40./100.;
            else if(ptBin == 3) mesonAmplitudeMin = mesonAmplitude*60./100.;
            else mesonAmplitudeMin = mesonAmplitude*30./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            if (fMode == 2 || fMode == 3){
                mesonAmplitudeMin = mesonAmplitude*10./100.;
            }
            if (fMode == 4 || fMode == 5){
                mesonAmplitudeMin = mesonAmplitude*5./100.;
            }
        }
    } else {
        if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
            mesonAmplitudeMin = mesonAmplitude*98./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            if (fEnergyFlag.CompareTo("pPb_5.023TeV") == 0) mesonAmplitudeMin = mesonAmplitude*92./100.;
            if (fMode == 2 || fMode == 3) {
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*600./100.;
            }
            if (fMode == 4 || fMode == 5) {
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                mesonAmplitudeMax = mesonAmplitude*400./100.;
                if( fEnergyFlag.CompareTo("8TeV") == 0 ){
                  mesonAmplitudeMin = mesonAmplitude*90./100.;
                  mesonAmplitudeMax = mesonAmplitude*400./100.;
                }
            }
            
        } else {
            mesonAmplitudeMin = mesonAmplitude*50./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            if (fMode == 2 || fMode == 3){
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                if( fEnergyFlag.CompareTo("8TeV") == 0 ){
                  mesonAmplitudeMin = mesonAmplitude*20./100.;
                  mesonAmplitudeMax = mesonAmplitude*115./100.;
                }
            }
            if (fMode == 4 || fMode == 5){
                mesonAmplitudeMin = mesonAmplitude*5./100.;
            }
        }
    }
    
    fFitReco= NULL;
    fFitReco = new TF1("GaussExpLinear","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",fMesonFitRange[0],fMesonFitRange[1]);

    fFitGausExp =NULL;
    fFitGausExp = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

    fFitLinearBck = NULL;
    fFitLinearBck = new TF1("Linear","[0]+[1]*x",fMesonFitRange[0],fMesonFitRange[1]);


    fFitReco->SetParameter(0,mesonAmplitude);
    fFitReco->SetParameter(1,fMesonMassExpect);
    fFitReco->SetParameter(2,fMesonWidthExpect);
    //    if(vary){
        fFitReco->SetParameter(3,fMesonLambdaTail);
    //    } else {
    //       fFitReco->FixParameter(3,fMesonLambdaTail);
    //    }
    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
//    fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
    fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.15);
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
    //    if(vary){
        fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);
    // }

    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");

    
    fFitReco->SetLineColor(3);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);

    if (vary && !fIsMC && (fMode == 0 || fMode == 9)){
        if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 && fPrefix.CompareTo("Pi0") ==0 && ptBin >=17){
            cout << "Skipping the vary option for this case, pt: " << ptBin << endl;
        }
        else if (fEnergyFlag.CompareTo("pPb_5.023TeV") == 0 && (ptBin >= 20) ){//
        cout << "Skipping the vary option for this case" << endl;
        } else {// ...do what you are supposed to....
        
        fMesonLambdaTail = fFitReco->GetParameter(3);
        fMesonLambdaTailRange[0] = 0.9*fFitReco->GetParameter(3);
        fMesonLambdaTailRange[1] = 1.1*fFitReco->GetParameter(3);
        fMesonWidthExpect = fFitReco->GetParameter(2);
        fMesonWidthRange[0] = 0.5*fFitReco->GetParameter(2);
        fMesonWidthRange[1] = 1.5*fFitReco->GetParameter(2);
    }
    }
    fFitGausExp->SetParameter(0,fFitReco->GetParameter(0));
    fFitGausExp->SetParameter(1,fFitReco->GetParameter(1));
    fFitGausExp->SetParameter(2,fFitReco->GetParameter(2));
    fFitGausExp->SetParameter(3,fFitReco->GetParameter(3));

    fFitGausExp->SetParError(0,fFitReco->GetParError(0));
    fFitGausExp->SetParError(1,fFitReco->GetParError(1));
    fFitGausExp->SetParError(2,fFitReco->GetParError(2));
    fFitGausExp->SetParError(3,fFitReco->GetParError(3));

    fFitLinearBck->SetParameter(0,fFitReco->GetParameter(4));
    fFitLinearBck->SetParameter(1,fFitReco->GetParameter(5));

    fFitLinearBck->SetParError(0,fFitReco->GetParError(4));
    fFitLinearBck->SetParError(1,fFitReco->GetParError(5));

    Int_t binCenterStart;
    Double_t startBinEdge;
    Int_t binCenterEnd;
    Double_t endBinEdge;

    TVirtualFitter * fitter = TVirtualFitter::GetFitter();

    fIntLinearBck = 0;
    fIntLinearBckError = 0;
    if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
        binCenterStart = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+fMesonIntDeltaRangeFit[0]);
        startBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
        binCenterEnd = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+fMesonIntDeltaRangeFit[1]);
        endBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

        Int_t nFreePar = fFitReco->GetNumberFreeParameters();
        double * covMatrix = fitter->GetCovarianceMatrix();

        Float_t intLinearBack = fFitLinearBck->GetParameter(0)*(endBinEdge-startBinEdge)+
            0.5*fFitLinearBck->GetParameter(1)*(endBinEdge*endBinEdge-startBinEdge*startBinEdge);

        Float_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fFitReco->GetParError(4),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fFitReco->GetParError(5),2)+2*covMatrix[nFreePar*nFreePar-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5);

        fFileDataLog << "Parameter for bin " << ptBin << endl;
        fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
        fFileDataLog << "Linear: \t"<<fFitReco->GetParameter(4)<<"+-" << fFitReco->GetParError(4) << "\t "<<fFitReco->GetParameter(5) <<"+-" << fFitReco->GetParError(5)<< endl;

        fIntLinearBck = intLinearBack/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
        fIntLinearBckError = errorLinearBck/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
    } else {
        fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }
    fFitReco->DrawCopy("same");
}

//****************************************************************************
//************* Fit of Signal+ BG with Gaussian + Linear BG ******************
//****************************************************************************
void FitSubtractedPureGaussianInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Int_t ptBin ){

    fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin;
    Double_t mesonAmplitudeMax;
    if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
        mesonAmplitudeMin = mesonAmplitude*98./100.;
        mesonAmplitudeMax = mesonAmplitude*115./100.;
        if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 || fEnergyFlag.CompareTo("pPb_5.023TeV") == 0) mesonAmplitudeMin = mesonAmplitude*92./100.;
        if (fMode == 2 || fMode == 3) {
            mesonAmplitudeMin = mesonAmplitude*98./100.;
            mesonAmplitudeMax = mesonAmplitude*400./100.;
        }
        if (fMode == 4 || fMode == 5) {
            mesonAmplitudeMin = mesonAmplitude*10./100.;
            mesonAmplitudeMax = mesonAmplitude*400./100.;
            if( fEnergyFlag.CompareTo("8TeV") == 0  ){
              mesonAmplitudeMin = mesonAmplitude*90./100.;
              mesonAmplitudeMax = mesonAmplitude*400./100.;
            }
        }
        
    } else {
        mesonAmplitudeMin = mesonAmplitude*50./100.;
        mesonAmplitudeMax = mesonAmplitude*120./100.;
        if (fMode == 2 || fMode == 3){
            mesonAmplitudeMin = mesonAmplitude*10./100.;
            if( fEnergyFlag.CompareTo("8TeV") == 0  ){
              mesonAmplitudeMin = mesonAmplitude*10./100.;
              mesonAmplitudeMax = mesonAmplitude*200./100.;
            }
        }
        if(fMode == 4){
            mesonAmplitudeMin = mesonAmplitude*30./100.;
            mesonAmplitudeMax = mesonAmplitude*220./100.;
        }
    }

    Double_t linBckg = 0.05;
    if(fEnergyFlag.CompareTo("8TeV") == 0 && fPrefix.CompareTo("Eta") == 0) linBckg = 0.1;

    fFitReco= NULL;
    fFitReco = new TF1("GaussLinearBG","gaus(0)+[3]+[4]*x",fMesonFitRange[0],fMesonFitRange[1]);

    fFitLinearBck = NULL;
    fFitLinearBck = new TF1("Linear","[0]+[1]*x",fMesonFitRange[1]-linBckg,fMesonFitRange[1]);

    fFitReco->SetParameter(0,mesonAmplitude);
    fFitReco->SetParameter(1,fMesonMassExpect);
    fFitReco->SetParameter(2,fMesonWidthExpect);
    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    fFitReco->SetParLimits(1,fMesonMassExpect*0.8,fMesonMassExpect*1.2);
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]*2);

//    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitLinearBck,"QRME0","",fMesonFitRange[1]-linBckg,fMesonFitRange[1]);
    fFitReco->SetParameter(3,fFitLinearBck->GetParameter(0));
    fFitReco->SetParLimits(3,fFitLinearBck->GetParameter(0)-2*fFitLinearBck->GetParError(0),fFitLinearBck->GetParameter(0)+2*fFitLinearBck->GetParError(0));
    fFitReco->SetParameter(4,fFitLinearBck->GetParameter(1));
    fFitReco->SetParLimits(4,fFitLinearBck->GetParameter(1)-2*fFitLinearBck->GetParError(1),fFitLinearBck->GetParameter(1)+2*fFitLinearBck->GetParError(1));
    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    //exclude second iteration of fitting, otherwise fits go completely wrong in 8 TeV
    //fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0","",fFitReco->GetParameter(1)-2*fFitReco->GetParameter(2),fFitReco->GetParameter(1)+2*fFitReco->GetParameter(2));
    
    fFitReco->SetLineColor(5);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);

    fFitLinearBck->SetParameter(0,fFitReco->GetParameter(3));
    fFitLinearBck->SetParameter(1,fFitReco->GetParameter(4));
    fFitLinearBck->SetParError(0,fFitReco->GetParError(3));
    fFitLinearBck->SetParError(1,fFitReco->GetParError(4));

    if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
        fFileDataLog << "Parameter for pure Gaussian bin " << ptBin << endl;
        fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<<endl;
        fFileDataLog << "Linear: \t"<<fFitReco->GetParameter(3)<<"+-" << fFitReco->GetParError(3) << "\t "<<fFitReco->GetParameter(4) <<"+-" << fFitReco->GetParError(4)<< endl;
    } else {
        fFileErrLog << "Pure Gaussian fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }
}

// Analog to the Funktion used for Crystalball
//****************************************************************************
//*** Fit of subtracted Signal+ BG with Gaus + tail + Lin BG        ******
//*** linear BG subtracted in this function after initial fit without   ******
//*** peak region, final fit only with Gaus + tail              ******
//*** additional outputs: fCopySignal - only Signal         ******
//***                     fCopyOnlyBG - only remaining BG       ******
//****************************************************************************
void GausFitSubtractedInvMassInPtBinsNew(TH1D* fHistoMappingSignalInvMassPtBinSingle,Double_t * fMesonIntDeltaRangeFit, Int_t ptBin,Bool_t vary,TString functionname ,Bool_t kMC){
    if(vary){}; 
    cout <<"Start Fitting spectra"<<endl;

    fCopySignal = (TH1D*)fHistoMappingSignalInvMassPtBinSingle->Clone("fCopySignal");
    fCopySignal->Sumw2();
    fCopyOnlyBG = (TH1D*)fHistoMappingSignalInvMassPtBinSingle->Clone("fCopyOnlyBG");
    fCopyOnlyBG->Sumw2();
    
    fFileErrLog<<"Start Fitting spectra with Gaus fit"<<endl;
    fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin = mesonAmplitude*10./100.;
    Double_t mesonAmplitudeMax = mesonAmplitude*400./100.;
    
//  fFitReco = NULL;
//  fFitReco = new TF1(functionname,"(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",fMesonFitRange[0],fMesonFitRange[1]);
    fFitReco = NULL;
    fFitReco = new TF1(functionname,"(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

    
    fFitGausExp = NULL;
    fFitGausExp = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

    fFitLinearBck = NULL;
    fFitLinearBck = new TF1("Linear","[0]+[1]*x",fMesonFitRange[0],fMesonFitRange[1]);

    fFitReco->SetParameter(0,mesonAmplitude);
    fFitReco->SetParameter(1,fMesonMassExpect);
    fFitReco->SetParameter(2,fMesonWidthExpect);
    fFitReco->SetParameter(3,fMesonLambdaTail);
    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.15);
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
    fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);
    
    if (!kMC){
        fFitLinearBckExcl = NULL;
        fFitLinearBckExcl = new TF1("LinearEx",LinearBGExclusionnew,fMesonMassPlotRange[0],fMesonMassPlotRange[1],2);
        fCopyOnlyBG->Fit(fFitLinearBckExcl,"QRME0","",fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
        
        fFitLinearBck->SetParameter(0,fFitLinearBckExcl->GetParameter(0));
        fFitLinearBck->SetParameter(1,fFitLinearBckExcl->GetParameter(1));
        fFitLinearBck->SetParError(0,fFitLinearBckExcl->GetParError(0));
        fFitLinearBck->SetParError(1,fFitLinearBckExcl->GetParError(1));

        
        TVirtualFitter * fitter2 = TVirtualFitter::GetFitter();
        Int_t nFreePar2 = fFitLinearBckExcl->GetNumberFreeParameters();
        double * covMatrix2 = fitter2->GetCovarianceMatrix();
        for (Int_t i = 1; i < fCopySignal->GetXaxis()->FindBin(fMesonMassRange[1])+1; i++){
            Double_t startBinEdge = fCopySignal->GetXaxis()->GetBinLowEdge(i);
            Double_t endBinEdge = fCopySignal->GetXaxis()->GetBinUpEdge(i);
            Double_t intLinearBack = fFitLinearBck->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
            Double_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fFitLinearBckExcl->GetParError(0),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fFitLinearBckExcl->GetParError(1),2)+2*covMatrix2[nFreePar2*nFreePar2-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
            fCopyOnlyBG->SetBinContent(i,intLinearBack);
            fCopyOnlyBG->SetBinError(i,errorLinearBck);
            fCopySignal->SetBinContent(i,fCopySignal->GetBinContent(i)-intLinearBack);
            fCopySignal->SetBinError(i,TMath::Sqrt(errorLinearBck*errorLinearBck+ fCopySignal->GetBinError(i)*fCopySignal->GetBinError(i)));
            //          fCopySignal->Add(fCopyOnlyBG,-1.);
//          cout << fFitLinearBck->Eval(startBinEdge) << "\t" <<fFitLinearBck->Eval(endBinEdge) << "\t" <<fCopySignal->GetBinContent(i) << "\t" <<fCopySignal->GetBinContent(i) << endl;
        }
        fCopySignal->Fit(fFitReco,"QRME0");
    } else {
        fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    }
    
//  fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");

    fFitReco->SetLineColor(3);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);

    fFitGausExp->SetParameter(0,fFitReco->GetParameter(0));
    fFitGausExp->SetParameter(1,fFitReco->GetParameter(1));
    fFitGausExp->SetParameter(2,fFitReco->GetParameter(2));
    fFitGausExp->SetParameter(3,fFitReco->GetParameter(3));
//  fFitGausExp->SetParameter(4,fFitReco->GetParameter(4));

    fFitGausExp->SetParError(0,fFitReco->GetParError(0));
    fFitGausExp->SetParError(1,fFitReco->GetParError(1));
    fFitGausExp->SetParError(2,fFitReco->GetParError(2));
    fFitGausExp->SetParError(3,fFitReco->GetParError(3));
//  fFitGausExp->SetParError(4,fFitReco->GetParError(4));

    fFitLinearBck->SetParameter(0,fFitLinearBckExcl->GetParameter(0));
    fFitLinearBck->SetParameter(1,fFitLinearBckExcl->GetParameter(1));
    fFitLinearBck->SetParError(0,fFitLinearBckExcl->GetParError(0));
    fFitLinearBck->SetParError(1,fFitLinearBckExcl->GetParError(1));

    Int_t binCenterStart;
    Double_t startBinEdge;
    Int_t binCenterEnd;
    Double_t endBinEdge;

    TVirtualFitter * fitter = TVirtualFitter::GetFitter();

    if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ){
        binCenterStart = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+fMesonIntDeltaRangeFit[0]);
        startBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
        binCenterEnd = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+fMesonIntDeltaRangeFit[1]);
        endBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

        Int_t nFreePar = fFitReco->GetNumberFreeParameters();
        double * covMatrix = fitter->GetCovarianceMatrix();

        if (!kMC){
            Float_t intLinearBack = fFitLinearBck->GetParameter(0)*(endBinEdge-startBinEdge)+
                0.5*fFitLinearBck->GetParameter(1)*(endBinEdge*endBinEdge-startBinEdge*startBinEdge);   
            
            Double_t errorConst = fFitReco->GetParError(5);
            Double_t errorLin = fFitReco->GetParError(6);
            if (errorConst == 0) errorConst = fFitLinearBck->GetParError(0);
            if (errorLin == 0) errorLin = fFitLinearBck->GetParError(1);
            if (errorConst == 0) errorConst = abs(fFitLinearBck->GetParameter(0)*0.005);
            if (errorLin == 0) errorLin = abs(fFitLinearBck->GetParameter(1)*0.005);
            Float_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*errorConst,2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*errorLin,2)+2*covMatrix[nFreePar*nFreePar-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5);
    
            fIntLinearBck = intLinearBack/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
            fIntLinearBckError = errorLinearBck/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
        
        } else {
            fIntLinearBck = 0;
            fIntLinearBckError = 0;         
        }   
    } 
    fFitReco->DrawCopy("same");
}




//****************************************************************************
//*** Fit of Pure MC Signal with Gaussian + Exponential **********************
//****************************************************************************
void FitTrueInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Double_t* fMesonIntDeltaRangeFit, Int_t ptBin, Bool_t vary)
{
    //    cout<<"Start Fitting spectra"<<endl;
    fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin;
    Double_t mesonAmplitudeMax;
    if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0){
        mesonAmplitudeMin = mesonAmplitude*99./100.;
        mesonAmplitudeMax = mesonAmplitude*110./100.;
        if (fMode == 2 || fMode == 3){
            mesonAmplitudeMin = mesonAmplitude*20./100.;
            mesonAmplitudeMax = mesonAmplitude*1000./100.;
        }
        if (fMode == 4 || fMode == 5) {
            mesonAmplitudeMin = mesonAmplitude*10./100.;
            mesonAmplitudeMax = mesonAmplitude*400./100.;
        }
    } else {
        mesonAmplitudeMin = mesonAmplitude*95./100.;
        mesonAmplitudeMax = mesonAmplitude*130./100.;
        if (fMode == 2 || fMode == 3){
            mesonAmplitudeMin = mesonAmplitude*95./100.;
            mesonAmplitudeMax = mesonAmplitude*1000./100.;
        }
        if (fMode == 4 || fMode == 5) {
            mesonAmplitudeMin = mesonAmplitude*10./100.;
            mesonAmplitudeMax = mesonAmplitude*400./100.;
            if( fEnergyFlag.CompareTo("8TeV") == 0  ){
              mesonAmplitudeMin = mesonAmplitude*90./100.;
              mesonAmplitudeMax = mesonAmplitude*400./100.;
            }
        }
    }
    
    fFitReco = NULL;
    TF1* fFitRecoPre = new TF1("fGauss","([0]*exp(-0.5*((x-[1])/[2])^2))", fMesonFitRange[0], fMesonFitRange[1]);
    if (fMode == 2 || fMode == 4){
        fFitReco = new TF1("GaussExpLinear","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",
                        fMesonFitRange[0], fMesonFitRange[1]);
    } else {
        fFitReco = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))", fMesonFitRange[0],
                        fMesonFitRange[1]);
    }

    fFitRecoPre->SetParameter(0,mesonAmplitude);
    fFitRecoPre->SetParameter(1,fMesonMassExpect);
    fFitRecoPre->SetParameter(2,fMesonWidthExpect);
    fFitRecoPre->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    fFitRecoPre->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
    if (fMode == 2 || fMode == 4) fFitRecoPre->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.15);
    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitRecoPre,"QRME0");
    
    fFitReco->SetParameter(3,fMesonLambdaTailMC);
    Double_t mass = fMesonMassExpect;
    if (fMode == 4){
        mass = fFitRecoPre->GetParameter(1);
        fFitReco->SetParameter(0,fFitRecoPre->GetParameter(0));
        fFitReco->SetParameter(1,mass);
        fFitReco->SetParameter(2,fFitRecoPre->GetParameter(2));
    } else {
        fFitReco->SetParameter(0,mesonAmplitude);
        fFitReco->SetParameter(1,fMesonMassExpect);
        fFitReco->SetParameter(2,fMesonWidthExpect);
    }
    
    if (fMode == 2){
        if (ptBin > fStartPtBin) fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    } else fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    
    if ( !(fMode == 2 || fMode == 4)) fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
    if (fMode == 2 ) fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.15);
    if (fMode == 4 ) fFitReco->SetParLimits(1,mass*0.95,mass*1.08);
//    if (fMode == 4) fFitReco->SetParLimits(1,fMesonMassExpect*0.97,fMesonMassExpect*1.05);
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
    fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);
    // //    fFitReco->SetParLimits(1,fMesonMassExpect*0.8,fMesonMassExpect*1.3);
    // //    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
    // //    if(vary){
    //       fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);
    
//    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");

    //    cout << TString(gMinuit->fCstatu.Data()).Data() << endl;

    if (vary && (fMode==9 || fMode ==0)){
        fMesonLambdaTailMC = fFitReco->GetParameter(3);
    //       fMesonWidthExpectMC = fMesonMassExpect*0.03;
        }
    fFitReco->SetLineColor(3);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);

    //    Int_t binCenterStart = 0;
    //    Double_t startBinEdge = 0;;
    //    Int_t binCenterEnd = 0;
    //    Double_t endBinEdge = 0;

    if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ){
    //       binCenterStart = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+fMesonIntDeltaRangeFit[0]);
    //       startBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
    //       binCenterEnd = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+fMesonIntDeltaRangeFit[1]);
    //       endBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

        fFileDataLog << "Parameter for bin " << ptBin << endl;
        fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
    } else {
        fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }
    fFitReco->DrawCopy("same");
    if (fMesonIntDeltaRangeFit){}
}

//****************************************************************************
//*** Fit of Pure MC Signal with Gaussian ************************************
//****************************************************************************
void FitTrueInvMassPureGaussianInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Int_t ptBin ){

    fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin;
    Double_t mesonAmplitudeMax;
    if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
        mesonAmplitudeMin = mesonAmplitude*98./100.;
        mesonAmplitudeMax = mesonAmplitude*115./100.;
        if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 || fEnergyFlag.CompareTo("pPb_5.023TeV") == 0) mesonAmplitudeMin = mesonAmplitude*92./100.;
        if (fMode == 2 || fMode == 3) {
            mesonAmplitudeMin = mesonAmplitude*98./100.;
            mesonAmplitudeMax = mesonAmplitude*400./100.;
        }
        if (fMode == 4 || fMode == 5) {
            mesonAmplitudeMin = mesonAmplitude*10./100.;
            mesonAmplitudeMax = mesonAmplitude*400./100.;
            if( fEnergyFlag.CompareTo("8TeV") == 0  ){
              mesonAmplitudeMin = mesonAmplitude*90./100.;
              mesonAmplitudeMax = mesonAmplitude*400./100.;
            }
        }
        
    } else {
        mesonAmplitudeMin = mesonAmplitude*50./100.;
        mesonAmplitudeMax = mesonAmplitude*120./100.;
        if (fMode == 2 || fMode == 3){
            mesonAmplitudeMin = mesonAmplitude*10./100.;
            if( fEnergyFlag.CompareTo("8TeV") == 0  ){
              mesonAmplitudeMin = mesonAmplitude*10./100.;
              mesonAmplitudeMax = mesonAmplitude*200./100.;
            }
        }
    }
    fFitReco= NULL;
    fFitReco = new TF1("GaussLinearBG","gaus(0)",fMesonFitRange[0]-0.05,fMesonFitRange[1]);

    fFitReco->SetParameter(0,mesonAmplitude);
    fFitReco->SetParameter(1,fMesonMassExpect);
    fFitReco->SetParameter(2,fMesonWidthExpect);
    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    fFitReco->SetParLimits(1,fMesonMassExpect*0.8,fMesonMassExpect*1.2);
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]*2);

    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    //exclude second iteration of fitting, otherwise fits go completely wrong in 8 TeV
    //fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0","",fFitReco->GetParameter(1)-2*fFitReco->GetParameter(2),fFitReco->GetParameter(1)+2*fFitReco->GetParameter(2));
    
    fFitReco->SetLineColor(5);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);

    if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
        fFileDataLog << "Parameter for pure Gaussian bin " << ptBin << endl;
        fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<<endl;
    } else {
        fFileErrLog << "Pure Gaussian fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }
}

//****************************************************************************
//*** Fit of Signal+ BG for symmetric decays (alpha < 0.1) with **************
//*** Gaussian + Exponential + quadratic BG **********************************
//****************************************************************************
void FitPeakPosInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Int_t ptBin, Bool_t vary)
{

    //    cout<<"Start Fitting spectra"<<endl;
    fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(0.05,0.3);
    Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin;
    Double_t mesonAmplitudeMax;
    if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
        mesonAmplitudeMin = mesonAmplitude*80./100.;
        mesonAmplitudeMax = mesonAmplitude*115./100.;
        if (fMode == 2 || fMode == 3) {
            mesonAmplitudeMin = mesonAmplitude*98./100.;
            mesonAmplitudeMax = mesonAmplitude*400./100.;
        }
        if (fMode == 4 || fMode == 5) {
            mesonAmplitudeMin = mesonAmplitude*10./100.;
            mesonAmplitudeMax = mesonAmplitude*400./100.;
            if( fEnergyFlag.CompareTo("8TeV") == 0  ){
              mesonAmplitudeMin = mesonAmplitude*90./100.;
              mesonAmplitudeMax = mesonAmplitude*400./100.;
            }
        }

    } else {
        mesonAmplitudeMin = mesonAmplitude*80./110.;
        mesonAmplitudeMax = mesonAmplitude*115./100.;
    }
    fFitReco = NULL;
    fFitReco = new TF1("GaussExpLinear","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x+[6]*x*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x+[6]*x*x)",0.05,0.3);

    fFitGausExp = NULL;
    fFitGausExp = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",0.05,0.3);

    fFitLinearBck = NULL;
    fFitLinearBck = new TF1("Linear","[0]+[1]*x+[2]*x*x",0.05,0.3);


    fFitReco->SetParameter(0,mesonAmplitude);
    fFitReco->SetParameter(1,fMesonMassExpect);
    fFitReco->SetParameter(2,fMesonWidthExpect);
    if(vary){
        fFitReco->SetParameter(3,fMesonLambdaTail);
    } else {
        fFitReco->FixParameter(3,fMesonLambdaTail);
    }
    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
    if(vary){fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);}

    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");

    fFitReco->SetLineColor(3);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);

    //if (vary) fMesonLambdaTail = fFitReco->GetParameter(3);

    fFitGausExp->SetParameter(0,fFitReco->GetParameter(0));
    fFitGausExp->SetParameter(1,fFitReco->GetParameter(1));
    fFitGausExp->SetParameter(2,fFitReco->GetParameter(2));
    fFitGausExp->SetParameter(3,fFitReco->GetParameter(3));

    fFitGausExp->SetParError(0,fFitReco->GetParError(0));
    fFitGausExp->SetParError(1,fFitReco->GetParError(1));
    fFitGausExp->SetParError(2,fFitReco->GetParError(2));
    fFitGausExp->SetParError(3,fFitReco->GetParError(3));

    fFitLinearBck->SetParameter(0,fFitReco->GetParameter(4));
    fFitLinearBck->SetParameter(1,fFitReco->GetParameter(5));
    fFitLinearBck->SetParameter(1,fFitReco->GetParameter(6));

    fFitLinearBck->SetParError(0,fFitReco->GetParError(4));
    fFitLinearBck->SetParError(1,fFitReco->GetParError(5));
    fFitLinearBck->SetParError(1,fFitReco->GetParError(6));

    if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ){
        fFileDataLog << "Parameter for bin " << ptBin << endl;
        fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
        fFileDataLog << "Quadratic: \t"<<fFitReco->GetParameter(4)<<"+-" << fFitReco->GetParError(4) << "\t "<<fFitReco->GetParameter(5) <<"+-" << fFitReco->GetParError(5) << "\t "<<fFitReco->GetParameter(6) <<"+-" << fFitReco->GetParError(6)<< endl;
    } else {
        fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }
    fFitReco->DrawCopy("same");
}

//****************************************************************************
//*** Fit of subtracted Signal+ BG with CrystalBall tail + Lin BG ************
//*** linear BG subtracted in this function after initial fit without ********
//*** peak region, final fit only with CrystalBall ***************************
//*** additional outputs: fCopySignal - only Signal                      ******
//***                     fCopyOnlyBG - only remaining BG                 ******
//****************************************************************************
void FitCBSubtractedInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle,Double_t * fMesonIntDeltaRangeFit, Int_t ptBin,Bool_t vary ,TString functionname, Bool_t kMC)
{
    if(vary){}; //dummy case to remove warning

    fCopySignal = (TH1D*)fHistoMappingSignalInvMassPtBinSingle->Clone("fCopySignal");
    fCopySignal->Sumw2();
    fCopyOnlyBG = (TH1D*)fHistoMappingSignalInvMassPtBinSingle->Clone("fCopyOnlyBG");
    fCopyOnlyBG->Sumw2();
    
    fFileErrLog<<"Start Fitting spectra with CB fit"<<endl;
    fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin = mesonAmplitude*50./100.;
    Double_t mesonAmplitudeMax = mesonAmplitude*200./100.;
    if (fMode == 4) mesonAmplitudeMin = mesonAmplitude*5./100.;
    
    fFitReco = NULL;
//    if (!kMC) {
//       fFitReco = new TF1(functionname,CrystalBallBck,fMesonFitRange[0],fMesonFitRange[1],7);
//    } else {
        fFitReco = new TF1(functionname,CrystalBall,fMesonFitRange[0],fMesonFitRange[1],5);
//    }
    
    fFitGausExp = NULL;
    fFitGausExp = new TF1("CrystalBall",CrystalBall,fMesonFitRange[0],fMesonFitRange[1],5);

    fFitLinearBck = NULL;
    fFitLinearBck = new TF1("Linear","[0]+[1]*x",fMesonFitRange[0],fMesonFitRange[1]);

    fFitReco->SetParameter(0,mesonAmplitude);
    fFitReco->SetParameter(1,fMesonMassExpect);
    fFitReco->SetParameter(2,fMesonWidthExpect);
    fFitReco->SetParameter(3,2.);  // n
    fFitReco->SetParameter(4,2. ); // alpha

    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
    
    if (!kMC){
        fFitLinearBckExcl = NULL;
        fFitLinearBckExcl = new TF1("LinearEx",LinearBGExclusion,fMesonFitRange[0],fMesonFitRange[1],2);
        fCopyOnlyBG->Fit(fFitLinearBckExcl,"QRME0","",fMesonFitRange[0],fMesonFitRange[1]);
        
        fFitLinearBck->SetParameter(0,fFitLinearBckExcl->GetParameter(0));
        fFitLinearBck->SetParameter(1,fFitLinearBckExcl->GetParameter(1));
        fFitLinearBck->SetParError(0,fFitLinearBckExcl->GetParError(0));
        fFitLinearBck->SetParError(1,fFitLinearBckExcl->GetParError(1));

        
        TVirtualFitter * fitter2 = TVirtualFitter::GetFitter();
        Int_t nFreePar2 = fFitLinearBckExcl->GetNumberFreeParameters();
        double * covMatrix2 = fitter2->GetCovarianceMatrix();
        for (Int_t i = 1; i < fCopySignal->GetXaxis()->FindBin(fMesonMassRange[1])+1; i++){
            Double_t startBinEdge = fCopySignal->GetXaxis()->GetBinLowEdge(i);
            Double_t endBinEdge = fCopySignal->GetXaxis()->GetBinUpEdge(i);
            Double_t intLinearBack = fFitLinearBck->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
            Double_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fFitLinearBckExcl->GetParError(0),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fFitLinearBckExcl->GetParError(1),2)+2*covMatrix2[nFreePar2*nFreePar2-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
            fCopyOnlyBG->SetBinContent(i,intLinearBack);
            fCopyOnlyBG->SetBinError(i,errorLinearBck);
            fCopySignal->SetBinContent(i,fCopySignal->GetBinContent(i)-intLinearBack);
            fCopySignal->SetBinError(i,TMath::Sqrt(errorLinearBck*errorLinearBck+ fCopySignal->GetBinError(i)*fCopySignal->GetBinError(i)));
            //          fCopySignal->Add(fCopyOnlyBG,-1.);
//          cout << fFitLinearBck->Eval(startBinEdge) << "\t" <<fFitLinearBck->Eval(endBinEdge) << "\t" <<fCopySignal->GetBinContent(i) << "\t" <<fCopySignal->GetBinContent(i) << endl;
        }
        fCopySignal->Fit(fFitReco,"QRME0");
    } else {
        fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    }
    
//    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");

    fFitReco->SetLineColor(3);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);

    fFitGausExp->SetParameter(0,fFitReco->GetParameter(0));
    fFitGausExp->SetParameter(1,fFitReco->GetParameter(1));
    fFitGausExp->SetParameter(2,fFitReco->GetParameter(2));
    fFitGausExp->SetParameter(3,fFitReco->GetParameter(3));
    fFitGausExp->SetParameter(4,fFitReco->GetParameter(4));

    fFitGausExp->SetParError(0,fFitReco->GetParError(0));
    fFitGausExp->SetParError(1,fFitReco->GetParError(1));
    fFitGausExp->SetParError(2,fFitReco->GetParError(2));
    fFitGausExp->SetParError(3,fFitReco->GetParError(3));
    fFitGausExp->SetParError(4,fFitReco->GetParError(4));

    fFitLinearBck->SetParameter(0,fFitLinearBckExcl->GetParameter(0));
    fFitLinearBck->SetParameter(1,fFitLinearBckExcl->GetParameter(1));
    fFitLinearBck->SetParError(0,fFitLinearBckExcl->GetParError(0));
    fFitLinearBck->SetParError(1,fFitLinearBckExcl->GetParError(1));

    Int_t binCenterStart;
    Double_t startBinEdge;
    Int_t binCenterEnd;
    Double_t endBinEdge;

    TVirtualFitter * fitter = TVirtualFitter::GetFitter();

    if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ){
        binCenterStart = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+fMesonIntDeltaRangeFit[0]);
        startBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
        binCenterEnd = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+fMesonIntDeltaRangeFit[1]);
        endBinEdge = fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

        Int_t nFreePar = fFitReco->GetNumberFreeParameters();
        double * covMatrix = fitter->GetCovarianceMatrix();

        if (!kMC){
            Float_t intLinearBack = fFitLinearBck->GetParameter(0)*(endBinEdge-startBinEdge)+
                0.5*fFitLinearBck->GetParameter(1)*(endBinEdge*endBinEdge-startBinEdge*startBinEdge);
            
            Double_t errorConst = fFitReco->GetParError(5);
            Double_t errorLin = fFitReco->GetParError(6);
            if (errorConst == 0) errorConst = fFitLinearBck->GetParError(0);
            if (errorLin == 0) errorLin = fFitLinearBck->GetParError(1);
            if (errorConst == 0) errorConst = abs(fFitLinearBck->GetParameter(0)*0.005);
            if (errorLin == 0) errorLin = abs(fFitLinearBck->GetParameter(1)*0.005);
            Float_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*errorConst,2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*errorLin,2)+2*covMatrix[nFreePar*nFreePar-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5);

            
            
            fFileDataLog << "Parameter for bin " << ptBin << endl;
            fFileDataLog << "CrystalBall: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<< "\t "<< fFitReco->GetParameter(4) <<"+-" << fFitReco->GetParError(4)<<endl;
            fFileDataLog << "Linear: \t"<<fFitReco->GetParameter(5)<<"+-" << errorConst << "\t "<<fFitReco->GetParameter(6) <<"+-" << errorLin<< endl;
        
            fIntLinearBck = intLinearBack/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
            fIntLinearBckError = errorLinearBck/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
            fFileDataLog << "Integrated BG: \t" << intLinearBack << "+-" <<  errorLinearBck << "\t bin width" <<fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10) <<endl;
            
        } else {
            fIntLinearBck = 0;
            fIntLinearBckError = 0;
        }
    } else {
        fFileErrLog << "Fitting failed in " << ptBin << " with status::" << gMinuit->fCstatu.Data() <<"why failed?"<<endl << endl;
    }
    fFitReco->DrawCopy("same");
}


//****************************************************************************
//*** Integration of Invariant Mass Histogram in given integration window ****
//****************************************************************************
void IntegrateHistoInvMass(TH1D * fHistoMappingSignalInvMassPtBinSingle, Double_t * fMesonIntRangeInt)
{
    Int_t binLowMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[0]);
    Int_t binHighMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[1]);
    fYields = fHistoMappingSignalInvMassPtBinSingle->IntegralAndError(binLowMassMeson,binHighMassMeson,fYieldsError);
}


//****************************************************************************
//*** Integration of Invariant Mass Histogram in given integration window ****
//*** with detailed output to log file ***************************************
//****************************************************************************
void IntegrateHistoInvMassStream(TH1D * fHistoMappingSignalInvMassPtBinSingle, Double_t * fMesonIntRangeInt) {
    Int_t binLowMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[0]);
    Int_t binHighMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[1]);
    fYields = fHistoMappingSignalInvMassPtBinSingle->IntegralAndError(binLowMassMeson,binHighMassMeson,fYieldsError);
    for ( Int_t M = binLowMassMeson; M < binHighMassMeson+1; M++){
        fFileDataLog << M << "\t" << fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(M) <<"\t" <<fHistoMappingSignalInvMassPtBinSingle->GetBinContent(M)<< "+-"<< fHistoMappingSignalInvMassPtBinSingle->GetBinError(M)<< endl;
    }
}

//****************************************************************************
//********* Integration of Fit function in given integration window **********
//****************************************************************************
void IntegrateFitFunc(TF1 * fFunc, TH1D *  fHistoMappingSignalInvMassPtBinSingle,Double_t * fMesonIntRangeInt) {
    fYieldsFunc = fFunc->Integral(fMesonIntRangeInt[0],fMesonIntRangeInt[1])/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
}


//****************************************************************************
//***************** Filling of MC histograms in proper binning ***************
//****************************************************************************
void FillHistosArrayMC(TH1D* fHistoMCMesonPtWithinAcceptanceFill, TH1D * fHistoMCMesonPtFill, TH1D * fDeltaPtFill) {
//    Char_t nameHisto[100] = "fHistoMCMesonPtEtaWithinAcceptance";
    fHistoMCMesonPtWithinAcceptanceFill->Sumw2();
    fHistoMCMesonWithinAccepPt = (TH1D*)fHistoMCMesonPtWithinAcceptanceFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
    fHistoMCMesonWithinAccepPt->Divide(fDeltaPtFill);
    fHistoMCMesonPtFill->Sumw2();
    fHistoMCMesonPt1 = (TH1D*)fHistoMCMesonPtFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
    fHistoMCMesonPt1->Divide(fDeltaPtFill);

}

//****************************************************************************
//***************** Filling of MC histograms in proper binning ***************
//****************************************************************************
void FillHistosArrayMCWOWeights(TH1D* fHistoMCMesonPtWithinAcceptanceFill, TH1D * fHistoMCMesonPtFill, TH1D * fDeltaPtFill) {
//    Char_t nameHisto[100] = "fHistoMCMesonPtEtaWithinAcceptance";
    fHistoMCMesonPtWithinAcceptanceFill->Sumw2();
    fHistoMCMesonWithinAccepPtWOWeights = (TH1D*)fHistoMCMesonPtWithinAcceptanceFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
    fHistoMCMesonWithinAccepPtWOWeights->Divide(fDeltaPtFill);
    fHistoMCMesonPtFill->Sumw2();
    fHistoMCMesonPt1WOWeights = (TH1D*)fHistoMCMesonPtFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
    fHistoMCMesonPt1WOWeights->Divide(fDeltaPtFill);

}

//****************************************************************************
//***************** Calculation of Meson Acceptance **************************
//****************************************************************************
void CalculateMesonAcceptance() {
    fHistoMCMesonAcceptPt = new TH1D("fMCMesonAccepPt","",fNBinsPt,fBinsPt);
    fHistoMCMesonAcceptPt->Sumw2();

    fHistoMCMesonAcceptPt->Divide(fHistoMCMesonWithinAccepPt,fHistoMCMesonPt1,1.,1.,"B");
    fHistoMCMesonAcceptPt->DrawCopy();
    fFileDataLog << endl << "Calculation of the Acceptance" << endl;
    for ( Int_t i = 1; i < fHistoMCMesonAcceptPt->GetNbinsX()+1 ; i++){
        fFileDataLog << "Bin " << i << "\t"<< fHistoMCMesonAcceptPt->GetBinCenter(i)<< "\t" << fHistoMCMesonAcceptPt->GetBinContent(i) << "\t" << fHistoMCMesonAcceptPt->GetBinError(i) <<endl;
    }
}

//****************************************************************************
//***************** Calculation of Meson Acceptance **************************
//****************************************************************************
void CalculateMesonAcceptanceWOWeights() {
    fHistoMCMesonAcceptPtWOWeights = new TH1D("fMCMesonAccepPtWOWeights","",fNBinsPt,fBinsPt);
    fHistoMCMesonAcceptPtWOWeights->Sumw2();

    fHistoMCMesonAcceptPtWOWeights->Divide(fHistoMCMesonWithinAccepPtWOWeights,fHistoMCMesonPt1WOWeights,1.,1.,"B");
    fHistoMCMesonAcceptPtWOWeights->DrawCopy();
    fFileDataLog << endl << "Calculation of the Acceptance wo weights" << endl;
    for ( Int_t i = 1; i < fHistoMCMesonAcceptPtWOWeights->GetNbinsX()+1 ; i++){
        fFileDataLog << "Bin " << i << "\t"<< fHistoMCMesonAcceptPtWOWeights->GetBinCenter(i)<< "\t" << fHistoMCMesonAcceptPtWOWeights->GetBinContent(i) << "\t" << fHistoMCMesonAcceptPtWOWeights->GetBinError(i) <<endl;
    }
}

//****************************************************************************
//***************** Calculation of Meson Efficiency **************************
//****************************************************************************
void CalculateMesonEfficiency(TH1D* fMC_fMesonYieldsPt, TH1D* fMC_SecondaryYieldPt, TString nameEfi ) {
    fHistoMCMesonEffiPt = new TH1D(nameEfi.Data(),"",fNBinsPt,fBinsPt);
        
    fHistoMCMesonEffiPt->Sumw2();
    fHistoMCMesonEffiPt->Add(fMC_fMesonYieldsPt,1.);
    if (fMC_SecondaryYieldPt) fHistoMCMesonEffiPt->Add(fMC_SecondaryYieldPt,-1.);
        
    fHistoMCMesonEffiPt->Divide(fHistoMCMesonEffiPt,fHistoMCMesonWithinAccepPt,1.,1.,"B");
    fFileDataLog << endl << "Calculation of the Efficiency" << nameEfi.Data()<< endl;
    for ( Int_t i = 1; i < fHistoMCMesonEffiPt->GetNbinsX()+1 ; i++){
        fFileDataLog << "Bin " << i << "\t" << fHistoMCMesonEffiPt->GetBinCenter(i)<< "\t"<< fHistoMCMesonEffiPt->GetBinContent(i) << "\t" << fHistoMCMesonEffiPt->GetBinError(i) <<endl;
    }
}

//****************************************************************************
//***************** Calculation of Meson Efficiency **************************
//****************************************************************************
void CalculateMesonEfficiencyWOWeights(TH1D* fMC_fMesonYieldsPt, TH1D* fMC_SecondaryYieldPt, TString nameEfi ) {
    fHistoMCMesonEffiPt = new TH1D(nameEfi.Data(),"",fNBinsPt,fBinsPt);
        
    fHistoMCMesonEffiPt->Sumw2();
    fHistoMCMesonEffiPt->Add(fMC_fMesonYieldsPt,1.);
    if (fMC_SecondaryYieldPt) fHistoMCMesonEffiPt->Add(fMC_SecondaryYieldPt,-1.);
        
    fHistoMCMesonEffiPt->Divide(fHistoMCMesonEffiPt,fHistoMCMesonWithinAccepPtWOWeights,1.,1.,"B");
    fFileDataLog << endl << "Calculation of the Efficiency" << nameEfi.Data()<< endl;
    for ( Int_t i = 1; i < fHistoMCMesonEffiPt->GetNbinsX()+1 ; i++){
        fFileDataLog << "Bin " << i << "\t" << fHistoMCMesonEffiPt->GetBinCenter(i)<< "\t"<< fHistoMCMesonEffiPt->GetBinContent(i) << "\t" << fHistoMCMesonEffiPt->GetBinError(i) <<endl;
    }
}

//****************************************************************************
//****** Saving of general histograms, fits needed in the analysis ***********
//****** RAW output file, no correction histograms ***************************
//****************************************************************************
void SaveHistos(Int_t optionMC, TString fCutID, TString fPrefix3, Bool_t UseTHnSparse ) {
    const char* nameOutput = Form("%s/%s/%s_%s_GammaConvV1WithoutCorrection%s_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix3.Data(),fPeriodFlag.Data(),fCutID.Data());
    fOutput1 = new TFile(nameOutput,"RECREATE");

    cout << "Begin writing Uncorrected File" << endl;
    
    Int_t fNBinsClusterPt           =  64;
    Double_t fBinsClusterPt[65]     =  {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                        1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                        2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8,
                                        4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8,
                                        6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.4, 7.8, 8.2, 8.6,
                                        9.0, 9.5, 10,  11,  12,  14,  16,  18,  20,  25, 
                                        30 , 35, 40, 45, 50 };

    TH1D*   fDeltaPtCluster       = new TH1D("fDeltaPtCluster","",fNBinsClusterPt,fBinsClusterPt);
    for(Int_t iPt=1;iPt<fNBinsClusterPt+1;iPt++){
        fDeltaPtCluster->SetBinContent(iPt,fBinsClusterPt[iPt]-fBinsClusterPt[iPt-1]);
        fDeltaPtCluster->SetBinError(iPt,0);
    }
    
    if (fHistoClustersPt){
        fHistoClustersPt->Write("ClusterPt");
        TH1D*   fHistoClustersPtPerEvent   = (TH1D*)fHistoClustersPt->Rebin(fNBinsClusterPt,"fHistoClustersPtPerEvent",fBinsClusterPt);
        fHistoClustersPtPerEvent->Divide(fDeltaPtCluster);
        fHistoClustersPtPerEvent->Scale(1./fNEvents);
        fHistoClustersPtPerEvent->Write("ClusterPtPerEvent");
    }
    if (fEnableDCCluster){
        TH1D*   fHistoTrueGammaClusPtRebinned   = NULL;
        TH1D*   fHistoTrueGammaDCClusPtRebinned   = NULL;
        if (fHistoTrueGammaClusPt){
            fHistoTrueGammaClusPt->Write();
            fHistoTrueGammaClusPtRebinned   = (TH1D*)fHistoTrueGammaClusPt->Rebin(fNBinsClusterPt,"fHistoTrueGammaClusPtRebinned",fBinsClusterPt);
        }
        if (fHistoTrueGammaDCClusPt){
            cout << "line" << __LINE__ << endl;
            fHistoTrueGammaDCClusPt->Write();
            cout << "line" << __LINE__ << endl;
            fHistoTrueGammaDCClusPtRebinned   = (TH1D*)fHistoTrueGammaDCClusPt->Rebin(fNBinsClusterPt,"fHistoTrueGammaDCClusPtRebinned",fBinsClusterPt);
        }
        cout << "line" << __LINE__ << endl;
        if (fHistoTrueGammaClusMultipleCount)fHistoTrueGammaClusMultipleCount->Write();
        if (fHistoTrueGammaClusPt && fHistoTrueGammaDCClusPt){
            cout << "line" << __LINE__ << endl;
            TH1F* fHistoRatioDCTrueGammaClus = (TH1F*)fHistoTrueGammaDCClusPt->Clone("fHistoRatioDCTrueGammaClus");
            fHistoRatioDCTrueGammaClus->Divide(fHistoRatioDCTrueGammaClus,fHistoTrueGammaClusPt,1,1,"B");
            fHistoRatioDCTrueGammaClus->Write("FractionDoubleCountedClusters");
            cout << "line" << __LINE__ << endl;
            TH1F* fHistoRatioDCTrueGammaClusRebinned = (TH1F*)fHistoTrueGammaDCClusPtRebinned->Clone("fHistoRatioDCTrueGammaClusRebinned");
            fHistoRatioDCTrueGammaClusRebinned->Divide(fHistoRatioDCTrueGammaClusRebinned,fHistoTrueGammaClusPtRebinned,1,1,"B");
            fHistoRatioDCTrueGammaClusRebinned->Write("FractionDoubleCountedClustersRebinned");
            cout << "line" << __LINE__ << endl;
        }
    }
    
    if (fHistoClustersOverlapHeadersPt){
        fHistoClustersOverlapHeadersPt->Write("ClusterOverlapHeadersPt");
        TH1D*   fHistoClustersOverlapHeadersPtPerEvent   = (TH1D*)fHistoClustersOverlapHeadersPt->Rebin(fNBinsClusterPt,"fHistoClustersOverlapHeadersPtPerEvent",fBinsClusterPt);
        fHistoClustersOverlapHeadersPtPerEvent->Divide(fDeltaPtCluster);
        fHistoClustersOverlapHeadersPtPerEvent->Scale(1./fNEvents);
        fHistoClustersOverlapHeadersPtPerEvent->Write("ClusterOverlapHeadersPtPerEvent");
    }
    fHistoYieldMeson->Write();
    fHistoYieldMesonPerEvent->Write();
    fHistoSignMeson->Write();
    fHistoSigndefaultMeson->Write();
    fHistoSBMeson->Write();
    fHistoSBdefaultMeson->Write();
    
    fHistoMassWindowHigh->Write();
    fHistoMassWindowLow->Write();
    fHistoMassWindowWideHigh->Write();
    fHistoMassWindowWideLow->Write();
    fHistoMassWindowNarrowHigh->Write();
    fHistoMassWindowNarrowLow->Write();

    fHistoYieldMesonNarrow->Write();
    fHistoYieldMesonPerEventNarrow->Write();
    fHistoSignMesonNarrow->Write();
    fHistoSBMesonNarrow->Write();
    fHistoSBdefaultNarrowMeson->Write();
    fHistoSigndefaultNarrowMeson->Write();

    fHistoLambdaTail->Write();

    fHistoYieldMesonWide->Write();
    fHistoYieldMesonPerEventWide->Write();
    fHistoSignMesonWide->Write();
    fHistoSBMesonWide->Write();

    if(fPrefix.CompareTo("Pi0") == 0 && fEnergyFlag.CompareTo("7TeV") == 0 ){fHistoMassPosition->Write();}
    if(fPrefix.CompareTo("Pi0") == 0 && fEnergyFlag.CompareTo("7TeV") == 0 ){fHistoFWHMMesonAlpha01->Write();}
    fHistoMassMeson->Write();
    fHistoMassGaussianMeson->Write();
    fHistoWidthGaussianMeson->Write();
    fHistoFWHMMeson->Write();
    fDeltaPt->Write();
    
    fHistoYieldMesonLeft->Write();
    fHistoYieldMesonLeftPerEvent->Write();
    fHistoSignMesonLeft->Write();
    fHistoSBMesonLeft->Write();

    fHistoYieldMesonLeftNarrow->Write();
    fHistoYieldMesonLeftPerEventNarrow->Write();
    fHistoSignMesonLeftNarrow->Write();
    fHistoSBMesonLeftNarrow->Write();

    fHistoYieldMesonLeftWide->Write();
    fHistoYieldMesonLeftPerEventWide->Write();
    fHistoSignMesonLeftWide->Write();
    fHistoSBMesonLeftWide->Write();

    fHistoMassMesonLeft->Write();
    fHistoFWHMMesonLeft->Write();
    fMesonFullPtSignal->Write();
    fMesonFullPtBackground->Write();
    fMesonFullPtBackNorm->SetName("Mapping_BackNorm_InvMass_FullPt");
    fMesonFullPtBackNorm->Write();
    fNumberOfGoodESDTracks->Write();
    fEventQuality->Write();

    TString nameHistoSignal;
    TString nameHistoSignalLeft;
    TString nameHistoPeakPos;
    TString nameHistoBckNorm;
    TString fitnameSignal;
//    if( fMode == 4) {
    TString nameHistoBckNormLeft;
    TString fitnameSignalLeft;
//     TString fitnameSignalLeft2;
//     TString fitnameSignal2;
//    }
    TString nameHistoSignalPos;
    if (fPrefix.CompareTo("Pi0") == 0 && fEnergyFlag.CompareTo("7TeV") == 0) {
        for (Int_t ii=fStartPtBin;ii<fNBinsPeakPt;ii++){
            nameHistoPeakPos = Form("InvMassAlpha01_in_Pt_Bin%02d", ii);
            fHistoMappingPeakPosInvMassPtBin[ii]->Write(nameHistoPeakPos.Data());
            fNameFitSignalPos = Form("Signal_InvMassFitPos_in_Pt_Bin%02d", ii);
            if(fFitPeakPosPtBin[ii]!=0x00) fFitPeakPosPtBin[ii]->Write(fNameFitSignalPos.Data());
        }
    }

    for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
        if(UseTHnSparse){
          if(!fUseRPBackground){
            fHistoWeightsBGZbinVsMbin[ii]->Write(Form("BGWeights_%02d", ii));
            fHistoFillPerEventBGZbinVsMbin[ii]->Scale(1./fNEvents);
            fHistoFillPerEventBGZbinVsMbin[ii]->Write(Form("BGPoolsFillstatus_%02d", ii));
          } else {
            fHistoWeightsBGZbinVsPsibin[ii]->Write(Form("BGWeights_%02d", ii));
            fHistoFillPerEventBGZbinVsPsibin[ii]->Scale(1./fNEvents);
            fHistoFillPerEventBGZbinVsPsibin[ii]->Write(Form("BGPoolsFillstatus_%02d", ii));
          }
        }

        fHistoMappingGGInvMassPtBin[ii]->Write();
        nameHistoBckNorm = Form("Mapping_BckNorm_InvMass_in_Pt_Bin%02d", ii);
        fHistoMappingBackNormInvMassPtBin[ii]->Write(nameHistoBckNorm.Data());
        nameHistoSignal = Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d", ii);
        fHistoMappingSignalInvMassPtBin[ii]->Write(nameHistoSignal.Data());
        nameHistoSignalLeft = Form("fHistoMappingSignalInvMassLeft_in_Pt_Bin%02d", ii); //Added 12.12.2014
        fHistoMappingSignalInvMassLeftPtBin[ii]->Write(nameHistoSignalLeft.Data());//Added 12.12.2014
        fitnameSignal = Form("Signal_InvMassFit_in_Pt_Bin%02d", ii);
        if( fMode == 4 ) {
        nameHistoBckNormLeft = Form("Mapping_BckNormLeft_InvMass_in_Pt_Bin%02d", ii);
        fHistoMappingBackNormInvMassLeftPtBin[ii]->Write(nameHistoBckNormLeft.Data());
        fitnameSignal = Form("Signalnew_InvMassFit_in_Pt_Bin%02d", ii);
        
        fitnameSignalLeft = Form("SignalnewLeft_InvMassFit_in_Pt_Bin%02d", ii);
        if(fFitInvMassLeftPtBin[ii]!=0x00) fFitInvMassLeftPtBin[ii]->Write(fitnameSignalLeft.Data());
//         fitnameSignalLeft2 = Form("SignalLeft_InvMassFit_in_Pt_Bin%02d", ii);
//         if(fFitInvMassLeftPtBin2[ii]!=0x00) fFitInvMassLeftPtBin2[ii]->Write(fitnameSignalLeft2.Data());
//         fitnameSignal2 = Form("Signal_InvMassFit_in_Pt_Bin%02d", ii);
//         if(fFitSignalInvMassPtBin2[ii]!=0x00) fFitSignalInvMassPtBin2[ii]->Write(fitnameSignal2.Data());
        
        }
        if(fFitSignalInvMassPtBin[ii]!=0x00) fFitSignalInvMassPtBin[ii]->Write(fitnameSignal.Data());
        
        if (fCrysFitting==1 ){ //|| fMode == 4
            if (fHistoMappingSignalRemainingBGSubInvMassPtBin[ii]) fHistoMappingSignalRemainingBGSubInvMassPtBin[ii]->Write();
            if (fHistoMappingSignalRemainingBGSubInvMassLeftPtBin[ii]) fHistoMappingSignalRemainingBGSubInvMassLeftPtBin[ii]->Write();
            if (fHistoMappingRemainingBGInvMassPtBin[ii]) fHistoMappingRemainingBGInvMassPtBin[ii]->Write();
            if (fHistoMappingRemainingBGInvMassLeftPtBin[ii]) fHistoMappingRemainingBGInvMassLeftPtBin[ii]->Write();
            if (fFitRemainingBGInvMassLeftPtBin[ii]) fFitRemainingBGInvMassLeftPtBin[ii]->Write();
            if (fFitRemainingBGInvMassPtBin[ii]) fFitRemainingBGInvMassPtBin[ii]->Write();
        }
    }

    if(optionMC){
        fHistoTrueSignMeson->Write();
        fHistoTrueSBMeson->Write();
        fHistoMCMesonPtWithinAcceptance->Write();
        fHistoMCMesonWithinAccepPt->Write(); // Proper bins in Pt
        if (fHistoMCMesonWithinAccepPtWOWeights) fHistoMCMesonWithinAccepPtWOWeights->Write();
        fHistoMCMesonPt1->Write(); // Proper bins in Pt
        if (fHistoMCMesonPt1WOWeights) fHistoMCMesonPt1WOWeights->Write(); // Proper bins in Pt
        fHistoYieldTrueMeson->Write();
        fHistoYieldTrueMesonWide->Write();
        fHistoYieldTrueMesonNarrow->Write();
        fHistoYieldTrueMesonReweighted->Write();
        fHistoYieldTrueMesonReweightedWide->Write();
        fHistoYieldTrueMesonReweightedNarrow->Write();
        fHistoYieldTrueMesonUnweighted->Write();
        fHistoYieldTrueMesonUnweightedWide->Write();
        fHistoYieldTrueMesonUnweightedNarrow->Write();
        
        fHistoTrueMesonInvMassVSPt->Write();
        fHistoYieldTrueSecMeson->Write();
        fHistoYieldTrueSecFromK0SMeson->Write();
        fHistoYieldTrueSecFromLambdaMeson->Write();
        fHistoYieldTrueSecMesonWide->Write();
        fHistoYieldTrueSecFromK0SMesonWide->Write();
        fHistoYieldTrueSecFromLambdaMesonWide->Write();
        fHistoYieldTrueSecMesonNarrow->Write();
        fHistoYieldTrueSecFromK0SMesonNarrow->Write();
        fHistoYieldTrueSecFromLambdaMesonNarrow->Write();
        for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
            fHistoMappingTrueMesonInvMassPtBins[ii]->Write();
            if (fEnableDCMeson) fHistoMappingTrueMesonDCInvMassPtBins[ii]->Write();
            fHistoMappingTrueFullMesonInvMassPtBins[ii]->Write();
            fHistoMappingTrueMesonInvMassPtReweightedBins[ii]->Write();
            fHistoMappingTrueMesonInvMassPtUnweightedBins[ii]->Write();
            if (fAdvancedMesonQA){
                if (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5){
                    if (fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[ii])fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[ii]->Write();
//                    if (fHistoMappingTrueMesonCaloElectronInvMassPtBins[ii])fHistoMappingTrueMesonCaloElectronInvMassPtBins[ii]->Write();
                    if (fHistoMappingTrueMesonCaloPhotonInvMassPtBins[ii])fHistoMappingTrueMesonCaloPhotonInvMassPtBins[ii]->Write();
                    if (fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[ii])fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[ii]->Write();
                    if (fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[ii])fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[ii]->Write();
                    if (fFitTrueSignalCaloConvPhotonInvMassPtBin[ii]!=0x00) fFitTrueSignalCaloConvPhotonInvMassPtBin[ii]->Write();
                    if (fFitTrueSignalCaloPhotonInvMassPtBin[ii]!=0x00) fFitTrueSignalCaloPhotonInvMassPtBin[ii]->Write();
                    if (fFitTrueSignalCaloElectronInvMassPtBin[ii]!=0x00) fFitTrueSignalCaloElectronInvMassPtBin[ii]->Write();
                    if (fFitTrueSignalCaloMergedClusterInvMassPtBin[ii]!=0x00) fFitTrueSignalCaloMergedClusterInvMassPtBin[ii]->Write();
                    if (fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[ii]!=0x00) fFitTrueSignalCaloMergedClusterPartConvInvMassPtBin[ii]->Write();
                } else {
                    fHistoMappingTrueGGBckInvMassPtBins[ii]->Write();
                    fHistoMappingTrueContBckInvMassPtBins[ii]->Write();
                    fHistoMappingTrueAllBckInvMassPtBins[ii]->Write();
                }
                if (fMode == 4 || fMode == 5){
                    if (fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[ii])fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[ii]->Write();
                    if (fFitTrueSignalMixedCaloConvPhotonInvMassPtBin[ii]!=0x00) fFitTrueSignalMixedCaloConvPhotonInvMassPtBin[ii]->Write();
                }
                
            }
            if (fHistoMappingTrueSecMesonInvMassPtBins[ii] != 0x00) fHistoMappingTrueSecMesonInvMassPtBins[ii]->Write();
            if (fHistoMappingTrueSecFromK0SMesonInvMassPtBins[ii] != 0x00) fHistoMappingTrueSecFromK0SMesonInvMassPtBins[ii]->Write();
            if (fHistoMappingTrueSecFromLambdaMesonInvMassPtBins[ii] != 0x00) fHistoMappingTrueSecFromLambdaMesonInvMassPtBins[ii]->Write();
            if (fFitTrueSignalInvMassPtBin[ii]!=0x00) fFitTrueSignalInvMassPtBin[ii]->Write();
            if (fFitTrueSignalInvMassPtReweightedBin[ii]!=0x00) fFitTrueSignalInvMassPtReweightedBin[ii]->Write();
            if (fFitTrueSignalInvMassPtUnweightedBin[ii]!=0x00) fFitTrueSignalInvMassPtUnweightedBin[ii]->Write();
        }
    }

    cout << "End writing Uncorrected File" << endl;
    
    fOutput1->Write();
    fOutput1->Close();
}


//****************************************************************************
//****** Saving of MC histograms needed for corrections  *********************
//****** or comparisons to data at a later stage *****************************
//****************************************************************************
void SaveCorrectionHistos(TString fCutID, TString fPrefix3){
    const char* nameOutput = Form("%s/%s/%s_%s_GammaConvV1CorrectionHistos%s_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix3.Data(),fPeriodFlag.Data(),fCutID.Data());
    fOutput2 = new TFile(nameOutput,"RECREATE");
    cout<<"======================================================"<<endl;
    cout<<"======================================================"<<endl;
    cout<<"======================================================"<<endl;
    cout<<nameOutput<<endl;
    cout<<"======================================================"<<endl;
    cout<<"======================================================"<<endl;
    cout<<"======================================================"<<endl;

    cout << "Begin writing Correction File" << endl;
    fHistoMCMesonAcceptPt->Write();
    fHistoTrueMesonInvMassVSPt->Write();
    fHistoMonteMesonEffiPt->Write();
    fHistoMonteMesonNarrowEffiPt->Write();
    fHistoMonteMesonWideEffiPt->Write();
    fHistoMonteMesonLeftEffiPt->Write();
    fHistoMonteMesonLeftNarrowEffiPt->Write();
    fHistoMonteMesonLeftWideEffiPt->Write();
    fHistoMCTrueMesonEffiPt->Write();
    fHistoMCTrueMesonEffiPtReweighted->Write();
    fHistoMCTrueMesonNarrowEffiPt->Write();
    fHistoMCTrueMesonNarrowEffiPtReweighted->Write();
    fHistoMCTrueMesonWideEffiPt->Write();
    fHistoMCTrueMesonWideEffiPtReweighted->Write();
    cout << "line" << __LINE__ << endl;
    if (fHistoMCMesonPtWithinAcceptanceWOWeights){
        fHistoMCMesonAcceptPtWOWeights->Write();
        fHistoMCTrueMesonEffiPtUnweighted->Write();
        fHistoMCTrueMesonNarrowEffiPtUnweighted->Write();
        fHistoMCTrueMesonWideEffiPtUnweighted->Write();
        fHistoMCMesonPt1WOWeights->SetName("MC_Meson_genPt_properBinning_WOWeights");
        fHistoMCMesonPt1WOWeights->Write(); // Proper bins in Pt
    }
    cout << "line" << __LINE__ << endl;
    fHistoMassMeson->Write();
    fHistoFWHMMeson->Write();
    fHistoMassGaussianMeson->Write();
    fHistoWidthGaussianMeson->Write();
    fHistoTrueMassMeson->Write();
    fHistoTrueMassGaussianMeson->Write();
    fHistoTrueWidthGaussianMeson->Write();
    fHistoTrueMassMesonReweighted->Write();
    fHistoTrueMassMesonUnweighted->Write();
    fHistoTrueFWHMMeson->Write();
    fHistoTrueFWHMMesonReweighted->Write();
    fHistoTrueFWHMMesonUnweighted->Write();
    cout << "line" << __LINE__ << endl;
    if (fAdvancedMesonQA && (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5)){
        if (fHistoTrueMassMesonCaloPhoton) fHistoTrueMassMesonCaloPhoton->Write();
        if (fHistoTrueMassMesonCaloElectron) fHistoTrueMassMesonCaloElectron->Write();
        if (fHistoTrueMassMesonCaloConvPhoton) fHistoTrueMassMesonCaloConvPhoton->Write();
        if (fHistoTrueMassMesonCaloMergedCluster) fHistoTrueMassMesonCaloMergedCluster->Write();
        if (fHistoTrueMassMesonCaloMergedPartConvCluster) fHistoTrueMassMesonCaloMergedPartConvCluster->Write();
        if (fHistoTrueFWHMMesonCaloPhoton) fHistoTrueFWHMMesonCaloPhoton->Write();
        if (fHistoTrueFWHMMesonCaloElectron) fHistoTrueFWHMMesonCaloElectron->Write();
        if (fHistoTrueFWHMMesonCaloConvPhoton) fHistoTrueFWHMMesonCaloConvPhoton->Write();
        if (fHistoTrueFWHMMesonCaloMergedCluster) fHistoTrueFWHMMesonCaloMergedCluster->Write();
        if (fHistoTrueFWHMMesonCaloMergedPartConvCluster) fHistoTrueFWHMMesonCaloMergedPartConvCluster->Write();
        if (fHistoYieldTrueMesonFixedWindow)fHistoYieldTrueMesonFixedWindow->Write();
        if (fHistoYieldTrueMesonGammaFixedWindow)fHistoYieldTrueMesonGammaFixedWindow->Write();
        if (fHistoYieldTrueMesonGammaConvGammaFixedWindow)fHistoYieldTrueMesonGammaConvGammaFixedWindow->Write();
        if (fHistoYieldTrueMesonConvGammaConvGammaFixedWindow)fHistoYieldTrueMesonConvGammaConvGammaFixedWindow->Write();
    }
    cout << "line" << __LINE__ << endl;
    if (fAdvancedMesonQA && (fMode == 4 || fMode == 5)){
        if (fHistoTrueMassMesonMixedCaloConvPhoton) fHistoTrueMassMesonMixedCaloConvPhoton->Write();
        if (fHistoTrueFWHMMesonMixedCaloConvPhoton) fHistoTrueFWHMMesonMixedCaloConvPhoton->Write();
    }
    cout << "line" << __LINE__ << endl;
    fHistoMCMesonPt1->SetName("MC_Meson_genPt");
    fHistoMCMesonPt1->Write(); // Proper bins in Pt
    fHistoMCMesonPt->SetName("MC_Meson_genPt_oldBin");
    fHistoMCMesonPt->Scale(1./fHistoMCMesonPt->GetBinWidth(5));
    //    fHistoMCMesonPt->GetXaxis()->SetRangeUser(0.,25.);
    fHistoMCMesonPt->Write(); 
    if (fHistoMCMesonPtWOWeights){
        fHistoMCMesonPtWOWeights->SetName("MC_Meson_genPt_WOWeights");
        fHistoMCMesonPtWOWeights->Scale(1./fHistoMCMesonPtWOWeights->GetBinWidth(5));
    //       fHistoMCMesonPtWOWeights->GetXaxis()->SetRangeUser(0.,25.);
        fHistoMCMesonPtWOWeights->Write(); 
    }
    if (fHistoMCMesonPtWeights){
        fHistoMCMesonPtWeights->Write("MC_Meson_genPt_Weights"); 
    }   
    cout << "line" << __LINE__ << endl;
    fEventQuality->Write();
    fHistoYieldTrueSecMeson->Write();
    fHistoYieldTrueSecFromK0SMeson->Write();
    fHistoYieldTrueSecFromLambdaMeson->Write();
    fHistoYieldTrueSecMesonWide->Write();
    fHistoYieldTrueSecFromK0SMesonWide->Write();
    fHistoYieldTrueSecFromLambdaMesonWide->Write();
    fHistoYieldTrueSecMesonNarrow->Write();
    fHistoYieldTrueSecFromK0SMesonNarrow->Write();
    fHistoYieldTrueSecFromLambdaMesonNarrow->Write();
    fHistoYieldTrueSecFracMeson->Write();
    fHistoYieldTrueSecFracFromK0SMeson->Write();
    fHistoYieldTrueSecFracFromLambdaMeson->Write();
    fHistoYieldTrueSecFracMesonWide->Write();
    fHistoYieldTrueSecFracFromK0SMesonWide->Write();
    fHistoYieldTrueSecFracFromLambdaMesonWide->Write();
    fHistoYieldTrueSecFracMesonNarrow->Write();
    fHistoYieldTrueSecFracFromK0SMesonNarrow->Write();
    fHistoYieldTrueSecFracFromLambdaMesonNarrow->Write();
    cout << "line" << __LINE__ << endl;
    fHistoYieldTrueMeson->Write();
    if (fEnableDCMeson){
        if (fHistoYieldTrueMesonDC) fHistoYieldTrueMesonDC->Write();
        if (fHistoYieldTrueMesonDC){
            TH1D* fHistoRatioDCTrueMeson = (TH1D*)fHistoYieldTrueMesonDC->Clone("fHistoRatioDCTrueMeson");
            fHistoRatioDCTrueMeson->Divide(fHistoRatioDCTrueMeson, fHistoYieldTrueMeson, 1, 1, "B");
            fHistoRatioDCTrueMeson->Write("FractionDoubleCountedMesons");
        }
        if (fHistoTrueMesonMultipleCount) fHistoTrueMesonMultipleCount->Write();
    }
    cout << "line" << __LINE__ << endl;
    fHistoYieldTrueMesonWide->Write();
    fHistoYieldTrueMesonNarrow->Write();
    fHistoYieldTrueMesonReweighted->Write();
    fHistoYieldTrueMesonReweightedWide->Write();
    fHistoYieldTrueMesonReweightedNarrow->Write();
    fHistoYieldTrueMesonUnweighted->Write();
    fHistoYieldTrueMesonUnweightedWide->Write();
    fHistoYieldTrueMesonUnweightedNarrow->Write();
    cout << "line" << __LINE__ << endl;
    if (fHistoYieldK0sWithPi0DaughterRec)fHistoYieldK0sWithPi0DaughterRec->Write("K0sWithPi0DaughterRec");
    if (fHistoYieldLambdaWithPi0DaughterRec)fHistoYieldLambdaWithPi0DaughterRec->Write("LambdaWithPi0DaughterRec");
    cout << "end writing Correction File" << endl;
    
    fOutput2->Write();
    fOutput2->Close();
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
    Double_t alpha = par[4];
    Double_t n = par[3];
    Double_t meanx = par[1];
    Double_t sigma = par[2];
    Double_t nn = par[0];
    Double_t a = TMath::Power((n/TMath::Abs(alpha)), n) * TMath::Exp(-0.5*alpha*alpha);
    Double_t b = n/TMath::Abs(alpha) - TMath::Abs(alpha);
    Double_t arg = (x[0] - meanx)/sigma;
    Double_t fitval = 0;
    if (arg > -1.0*alpha) {
        fitval = nn * TMath::Exp(-0.5*arg*arg);
    } else {
        fitval = nn * a * TMath::Power((b-arg), (-1*n));
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
//******************** Definition of linear BG fit with expluded  ************
//******************** region fBGFitRangeLeft[1]-fBGFitRange[0] **************
//*******parameters are:                                                *****
//*******               - 0 constant BG                                 *****
//*******               - 1 linear BG                                   *****
//****************************************************************************
Double_t LinearBGExclusion(Double_t *x, Double_t *par) {
    if (x[0] > fBGFitRangeLeft[1] && x[0] < fBGFitRange[0]) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0];
}

// Using an different Fitrange for Mode 4
Double_t LinearBGExclusionnew(Double_t *x, Double_t *par) {
    if (x[0] > fBGFitRangeLeft[1] && x[0] < fBGFitRange[0]) {
        TF1::RejectPoint();
        return 0;
    } 
    if (x[0] < fBGFitRangeLeft[0]) {
        TF1::RejectPoint();
        return 0;
    } 
    if (x[0] > fBGFitRange[1]) {
        TF1::RejectPoint();
        return 0;
    } 
    return par[0] + par[1]*x[0];
}


//****************************************************************************
//******** Calculation of FWHM for Gaussian + left side exponential  *********
//****************************************************************************
void CalculateFWHM(TF1 * fFunc){
// Default function
    if (fCrysFitting == 0){
        TF1* fFunc_def;
        fFunc_def = new TF1("fFunc_def","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);
        fFunc_def->SetParameter(0,fFunc->GetParameter(0));
        fFunc_def->SetParameter(1,fFunc->GetParameter(1));
        fFunc_def->SetParameter(2,fFunc->GetParameter(2));
        fFunc_def->SetParameter(3,fFunc->GetParameter(3));
    
        //FWHM
        fFWHMFunc = fFunc_def->GetX(fFunc_def->GetParameter(0)*0.5,fFunc_def->GetParameter(1), fMesonFitRange[1]) - fFunc_def->GetX(fFunc_def->GetParameter(0)*0.5,fMesonFitRange[0],fFunc_def->GetParameter(1));

        //FWHM error +
        TF1* fFunc_plus;
        fFunc_plus = new TF1("fFunc_plus","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);
        fFunc_plus->SetParameter(0,fFunc->GetParameter(0) + fFunc->GetParError(0));
        fFunc_plus->SetParameter(1,fFunc->GetParameter(1) + fFunc->GetParError(1));
        fFunc_plus->SetParameter(2,fFunc->GetParameter(2) + fFunc->GetParError(2));
        fFunc_plus->SetParameter(3,fFunc->GetParameter(3) + fFunc->GetParError(3));
        Double_t FWHM_plus = fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,fFunc_plus->GetParameter(1), fMesonFitRange[1]) - fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,fMesonFitRange[0],fFunc_plus->GetParameter(1));

        //FWHM error -
        TF1* fFunc_minus;
        //   fFunc_minus = fFunc;
        fFunc_minus = new TF1("fFunc_minus","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);
        fFunc_minus->SetParameter(0,fFunc->GetParameter(0) - fFunc->GetParError(0));
        fFunc_minus->SetParameter(1,fFunc->GetParameter(1) - fFunc->GetParError(1));
        fFunc_minus->SetParameter(2,fFunc->GetParameter(2) - fFunc->GetParError(2));
        fFunc_minus->SetParameter(3,fFunc->GetParameter(3) - fFunc->GetParError(3));
        
        Double_t FWHM_minus =  fFunc_minus->GetX(fFunc_minus->GetParameter(0)*0.5,fFunc_minus->GetParameter(1), fMesonFitRange[1]) -fFunc_minus->GetX(fFunc_minus->GetParameter(0)*0.5,fMesonFitRange[0],fFunc_minus->GetParameter(1));
        Double_t Error1 = TMath::Abs(fFWHMFunc-FWHM_plus);
        Double_t Error2 = TMath::Abs(fFWHMFunc-FWHM_minus);
        if(Error1>=Error2) fFWHMFuncError = Error1;
        if(Error1<Error2) fFWHMFuncError = Error2;
    } else {
        fFWHMFunc = fFunc->GetParameter(2)*2.35;
        fFWHMFuncError = fFunc->GetParError(2)*2.35;
    }
}

//****************************************************************************
//****** Deleting all pointers generated during this analysis ****************
//****************************************************************************
void Delete(){
    if (fBinsPt)                                                delete fBinsPt;
    if (fPeakRange)                                             delete fPeakRange;
    if (fIntFixedRange)                                         delete fIntFixedRange;
    if (fFitRange)                                              delete fFitRange;
    if (fBGFitRange)                                            delete fBGFitRange;
    if (fBGFitRangeLeft)                                        delete fBGFitRangeLeft;
    if (fMesonPlotRange)                                        delete fMesonPlotRange;
    if (fMesonIntDeltaRange)                                    delete fMesonIntDeltaRange;
    if (fMesonIntDeltaRangeWide)                                delete fMesonIntDeltaRangeWide;
    if (fMesonIntDeltaRangeNarrow)                              delete fMesonIntDeltaRangeNarrow;
    if (fMesonMassRange)                                        delete fMesonMassRange;
    if (fMesonMassPlotRange)                                    delete fMesonMassPlotRange;
    if (fMesonFitRange)                                         delete fMesonFitRange;
    if (fMesonWidthRange)                                       delete fMesonWidthRange;
    if (fMesonLambdaTailRange)                                  delete fMesonLambdaTailRange;
    if (fNRebin)                                                delete fNRebin;
    if (fGGYields)                                              delete fGGYields;
    if (fBckYields)                                             delete fBckYields;
    if (fMesonYields)                                           delete fMesonYields;
    if (fMesonTrueYields)                                       delete fMesonTrueYields;
    if (fMesonTrueYieldsDC)                                     delete fMesonTrueYieldsDC;
    if (fMesonTrueYieldFixedWindow)                             delete fMesonTrueYieldFixedWindow;
    if (fMesonTrueYieldGammaConvGammaFixedWindow)               delete fMesonTrueYieldGammaConvGammaFixedWindow;
    if (fMesonTrueYieldConvGammaConvGammaFixedWindow)           delete fMesonTrueYieldConvGammaConvGammaFixedWindow;
    if (fMesonTrueYieldGammaFixedWindow)                        delete fMesonTrueYieldGammaFixedWindow;
    if (fMesonTrueYieldErrorFixedWindow)                        delete fMesonTrueYieldErrorFixedWindow;
    if (fMesonTrueYieldGammaErrorFixedWindow)                   delete fMesonTrueYieldGammaErrorFixedWindow;
    if (fMesonTrueYieldGammaConvGammaErrorFixedWindow)          delete fMesonTrueYieldGammaConvGammaErrorFixedWindow;
    if (fMesonTrueYieldConvGammaConvGammaErrorFixedWindow)      delete fMesonTrueYieldConvGammaConvGammaErrorFixedWindow;
    if (fMesonTrueYieldsReweighted)                             delete fMesonTrueYieldsReweighted;
    if (fMesonTrueYieldsUnweighted)                             delete fMesonTrueYieldsUnweighted;
    if (fMesonYieldsFunc)                                       delete fMesonYieldsFunc;
    if (fMesonYieldsResidualBckFunc)                           delete fMesonYieldsResidualBckFunc;
    if (fMesonYieldsCorResidualBckFunc)                         delete fMesonYieldsCorResidualBckFunc;
    if (fMesonYieldsPerEvent)                                   delete fMesonYieldsPerEvent;
    if (fMesonMass)                                             delete fMesonMass;
    if (fMesonSB)                                               delete fMesonSB;
    if (fMesonSBdefault)                                        delete fMesonSBdefault;
    if (fMesonSign)                                             delete fMesonSign;
    if (fMesonSigndefault)                                      delete fMesonSigndefault;
    if (fMassWindowHigh)                                        delete fMassWindowHigh;
    if (fMassWindowLow)                                         delete fMassWindowLow;
    if (fMassWindowWideHigh)                                    delete fMassWindowWideHigh;
    if (fMassWindowWideLow)                                     delete fMassWindowWideLow;
    if (fMassWindowNarrowHigh)                                  delete fMassWindowNarrowHigh;
    if (fMassWindowNarrowLow)                                   delete fMassWindowNarrowLow;
    if (fMesonLambdaTailpar)                                    delete fMesonLambdaTailpar;
    if (fMesonLambdaTailparError)                               delete fMesonLambdaTailparError;
    if (fMesonLambdaTailMCpar)                                  delete fMesonLambdaTailMCpar;
    if (fMesonLambdaTailMCparError)                             delete fMesonLambdaTailMCparError;

    if (fMesonTrueSB)                                           delete fMesonTrueSB;
    if (fMesonTrueSign)                                         delete fMesonTrueSign;
    if (fMesonFWHM)                                             delete fMesonFWHM;
    if (fGGYieldsLeft)                                          delete fGGYieldsLeft;
    if (fBckYieldsLeft)                                         delete fBckYieldsLeft;
    if (fMesonYieldsLeft)                                       delete fMesonYieldsLeft;
    if (fMesonYieldsFuncLeft)                                   delete fMesonYieldsFuncLeft;
    if (fMesonYieldsResidualBckFuncLeft)                        delete fMesonYieldsResidualBckFuncLeft;
    if (fMesonYieldsCorResidualBckFuncLeft)                     delete fMesonYieldsCorResidualBckFuncLeft;
    if (fMesonYieldsLeftPerEvent)                               delete fMesonYieldsLeftPerEvent;
    if (fMesonMassLeft)                                         delete fMesonMassLeft;
    if (fMesonSBLeft)                                           delete fMesonSBLeft;
    if (fMesonSignLeft)                                         delete fMesonSignLeft;
    if (fMesonFWHMLeft)                                         delete fMesonFWHMLeft;
    if (fGGYieldsNarrow)                                        delete fGGYieldsNarrow;
    if (fBckYieldsNarrow)                                       delete fBckYieldsNarrow;
    if (fMesonYieldsNarrow)                                     delete fMesonYieldsNarrow;
    if (fMesonTrueYieldsNarrow)                                 delete fMesonTrueYieldsNarrow;
    if (fMesonTrueYieldsReweightedNarrow)                       delete fMesonTrueYieldsReweightedNarrow;
    if (fMesonTrueYieldsUnweightedNarrow)                       delete fMesonTrueYieldsUnweightedNarrow;
    if (fMesonYieldsFuncNarrow)                                 delete fMesonYieldsFuncNarrow;
    if (fMesonYieldsResidualBckFuncNarrow)                      delete fMesonYieldsResidualBckFuncNarrow;
    if (fMesonYieldsCorResidualBckFuncNarrow)                   delete fMesonYieldsCorResidualBckFuncNarrow;
    if (fMesonYieldsPerEventNarrow)                               delete fMesonYieldsPerEventNarrow;
    if (fMesonSBNarrow)                                         delete fMesonSBNarrow;
    if (fMesonSBdefaultNarrow)                                  delete fMesonSBdefaultNarrow;
    if (fMesonSignNarrow)                                       delete fMesonSignNarrow;
    if (fMesonSigndefaultNarrow)                                delete fMesonSigndefaultNarrow;
    if (fGGYieldsLeftNarrow)                                    delete fGGYieldsLeftNarrow;
    if (fBckYieldsLeftNarrow)                                   delete fBckYieldsLeftNarrow;
    if (fMesonYieldsLeftNarrow)                                 delete fMesonYieldsLeftNarrow;
    if (fMesonYieldsFuncLeftNarrow)                             delete fMesonYieldsFuncLeftNarrow;
    if (fMesonYieldsResidualBckFuncLeftNarrow)                  delete fMesonYieldsResidualBckFuncLeftNarrow;
    if (fMesonYieldsCorResidualBckFuncLeftNarrow)               delete fMesonYieldsCorResidualBckFuncLeftNarrow;
    if (fMesonYieldsLeftPerEventNarrow)                         delete fMesonYieldsLeftPerEventNarrow;
    if (fMesonSBLeftNarrow)                                     delete fMesonSBLeftNarrow;
    if (fMesonSignLeftNarrow)                                   delete fMesonSignLeftNarrow;
    if (fGGYieldsWide)                                          delete fGGYieldsWide;
    if (fBckYieldsWide)                                         delete fBckYieldsWide;
    if (fMesonYieldsWide)                                       delete fMesonYieldsWide;
    if (fMesonTrueYieldsWide)                                   delete fMesonTrueYieldsWide;
    if (fMesonTrueYieldsReweightedWide)                         delete fMesonTrueYieldsReweightedWide;
    if (fMesonTrueYieldsUnweightedWide)                         delete fMesonTrueYieldsUnweightedWide;
    if (fMesonYieldsFuncWide)                                   delete fMesonYieldsFuncWide;
    if (fMesonYieldsResidualBckFuncWide)                        delete fMesonYieldsResidualBckFuncWide;
    if (fMesonYieldsCorResidualBckFuncWide)                     delete fMesonYieldsCorResidualBckFuncWide;
    if (fMesonYieldsPerEventWide)                               delete fMesonYieldsPerEventWide;
    if (fMesonSBWide)                                           delete fMesonSBWide;
    if (fMesonSignWide)                                         delete fMesonSignWide;
    if (fGGYieldsLeftWide)                                      delete fGGYieldsLeftWide;
    if (fBckYieldsLeftWide)                                     delete fBckYieldsLeftWide;
    if (fMesonYieldsLeftWide)                                   delete fMesonYieldsLeftWide;
    if (fMesonYieldsFuncLeftWide)                               delete fMesonYieldsFuncLeftWide;
    if (fMesonYieldsResidualBckFuncLeftWide)                    delete fMesonYieldsResidualBckFuncLeftWide;
    if (fMesonYieldsCorResidualBckFuncLeftWide)                 delete fMesonYieldsCorResidualBckFuncLeftWide;
    if (fMesonYieldsLeftPerEventWide)                           delete fMesonYieldsLeftPerEventWide;
    if (fMesonSBLeftWide)                                       delete fMesonSBLeftWide;
    if (fMesonSignLeftWide)                                     delete fMesonSignLeftWide;
    if (fBckYieldsError)                                        delete fBckYieldsError;
    if (fMesonYieldsError)                                      delete fMesonYieldsError;
    if (fMesonYieldsFuncError)                                  delete fMesonYieldsFuncError;
    if (fMesonYieldsResidualBckFuncError)                       delete fMesonYieldsResidualBckFuncError;
    if (fMesonYieldsCorResidualBckFuncError)                    delete fMesonYieldsCorResidualBckFuncError;
    if (fMesonYieldsPerEventError)                              delete fMesonYieldsPerEventError;
    if (fMesonMassError)                                        delete fMesonMassError;
    if (fMesonSBError)                                          delete fMesonSBError;
    if (fMesonSBdefaultError)                                   delete fMesonSBdefaultError;
    if (fMesonSBdefaultNarrowError)                             delete fMesonSBdefaultNarrowError;
    if (fMesonSigndefaultNarrowError)                           delete fMesonSigndefaultNarrowError;
    if (fMesonSignError)                                       delete fMesonSignError;
    if (fMesonSigndefaultError)                                   delete fMesonSigndefaultError;
    if (fMesonTrueSBError)                                      delete fMesonTrueSBError;
    if (fMesonTrueSignError)                                    delete fMesonTrueSignError;
    if (fMesonFWHMError)                                        delete fMesonFWHMError;
    if (fGGYieldsLeftError)                                     delete fGGYieldsLeftError;
    if (fBckYieldsLeftError)                                    delete fBckYieldsLeftError;
    if (fMesonYieldsLeftError)                                  delete fMesonYieldsLeftError;
    if (fMesonYieldsFuncLeftError)                              delete fMesonYieldsFuncLeftError;
    if (fMesonYieldsResidualBckFuncLeftError)                   delete fMesonYieldsResidualBckFuncLeftError;
    if (fMesonYieldsCorResidualBckFuncLeftError)                delete fMesonYieldsCorResidualBckFuncLeftError;
    if (fMesonYieldsLeftPerEventError)                          delete fMesonYieldsLeftPerEventError;
    if (fMesonMassLeftError)                                    delete fMesonMassLeftError;
    if (fMesonSBLeftError)                                      delete fMesonSBLeftError;
    if (fMesonSignLeftError)                                    delete fMesonSignLeftError;
    if (fMesonFWHMLeftError)                                    delete fMesonFWHMLeftError;
    if (fGGYieldsNarrowError)                                   delete fGGYieldsNarrowError;
    if (fBckYieldsNarrowError)                                  delete fBckYieldsNarrowError;
    if (fMesonYieldsNarrowError)                                delete fMesonYieldsNarrowError;
    if (fMesonYieldsFuncNarrowError)                            delete fMesonYieldsFuncNarrowError;
    if (fMesonYieldsResidualBckFuncNarrowError)                 delete fMesonYieldsResidualBckFuncNarrowError;
    if (fMesonYieldsCorResidualBckFuncNarrowError)              delete fMesonYieldsCorResidualBckFuncNarrowError;
    if (fMesonYieldsPerEventNarrowError)                       delete fMesonYieldsPerEventNarrowError;
    if (fGGYieldsLeftNarrowError)                               delete fGGYieldsLeftNarrowError;
    if (fBckYieldsLeftNarrowError)                              delete fBckYieldsLeftNarrowError;
    if (fMesonYieldsLeftNarrowError)                            delete fMesonYieldsLeftNarrowError;
    if (fMesonYieldsFuncLeftNarrowError)                        delete fMesonYieldsFuncLeftNarrowError;
    if (fMesonYieldsResidualBckFuncLeftNarrowError)             delete fMesonYieldsResidualBckFuncLeftNarrowError;
    if (fMesonYieldsCorResidualBckFuncLeftNarrowError)          delete fMesonYieldsCorResidualBckFuncLeftNarrowError;
    if (fMesonYieldsLeftPerEventNarrowError)                    delete fMesonYieldsLeftPerEventNarrowError;
    if (fMesonSBLeftNarrowError)                                delete fMesonSBLeftNarrowError;
    if (fMesonSignLeftNarrowError)                              delete fMesonSignLeftNarrowError;
    if (fGGYieldsWideError)                                     delete fGGYieldsWideError;
    if (fBckYieldsWideError)                                    delete fBckYieldsWideError;
    if (fMesonYieldsWideError)                                  delete fMesonYieldsWideError;
    if (fMesonYieldsFuncWideError)                              delete fMesonYieldsFuncWideError;
    if (fMesonYieldsResidualBckFuncWideError)                   delete fMesonYieldsResidualBckFuncWideError;
    if (fMesonYieldsCorResidualBckFuncWideError)                delete fMesonYieldsCorResidualBckFuncWideError;
    if (fMesonYieldsPerEventWideError)                          delete fMesonYieldsPerEventWideError;
    if (fMesonSBWideError)                                      delete fMesonSBWideError;
    if (fMesonSignWideError)                                    delete fMesonSignWideError;
    if (fGGYieldsLeftWideError)                                 delete fGGYieldsLeftWideError;
    if (fBckYieldsLeftWideError)                                delete fBckYieldsLeftWideError;
    if (fMesonYieldsLeftWideError)                              delete fMesonYieldsLeftWideError;
    if (fMesonYieldsFuncLeftWideError)                          delete fMesonYieldsFuncLeftWideError;
    if (fMesonYieldsResidualBckFuncLeftWideError)               delete fMesonYieldsResidualBckFuncLeftWideError;
    if (fMesonYieldsCorResidualBckFuncLeftWideError)            delete fMesonYieldsCorResidualBckFuncLeftWideError;
    if (fMesonYieldsLeftPerEventWideError)                      delete fMesonYieldsLeftPerEventWideError;
    if (fMesonSBLeftWideError)                                  delete fMesonSBLeftWideError;
    if (fMesonSignLeftWideError)                                delete fMesonSignLeftWideError;
    if (fHistoMappingTrueMesonInvMassPtBins)                    delete fHistoMappingTrueMesonInvMassPtBins;
    if (fHistoMappingTrueMesonDCInvMassPtBins)                  delete fHistoMappingTrueMesonDCInvMassPtBins;
    if (fHistoMappingTrueFullMesonInvMassPtBins)                delete fHistoMappingTrueFullMesonInvMassPtBins;
    if (fHistoMappingTrueMesonInvMassPtReweightedBins)          delete fHistoMappingTrueMesonInvMassPtReweightedBins;
    if (fHistoMappingTrueMesonInvMassPtUnweightedBins)          delete fHistoMappingTrueMesonInvMassPtUnweightedBins;
    if (fHistoMappingTrueGGBckInvMassPtBins)                    delete fHistoMappingTrueGGBckInvMassPtBins;
    if (fHistoMappingTrueContBckInvMassPtBins)                  delete fHistoMappingTrueContBckInvMassPtBins;
    if (fHistoMappingTrueAllBckInvMassPtBins)                   delete fHistoMappingTrueAllBckInvMassPtBins;
    if (fHistoMappingGGInvMassPtBin)                            delete fHistoMappingGGInvMassPtBin;
    if (fHistoMappingBackInvMassPtBin)                          delete fHistoMappingBackInvMassPtBin;
    if (fHistoMappingBackNormInvMassPtBin)                      delete fHistoMappingBackNormInvMassPtBin;
    if (fHistoMappingSignalInvMassPtBin)                        delete fHistoMappingSignalInvMassPtBin;
    if (fHistoMappingSignalRemainingBGSubInvMassPtBin)          delete fHistoMappingSignalRemainingBGSubInvMassPtBin;
    if (fHistoMappingSignalRemainingBGSubInvMassLeftPtBin)      delete fHistoMappingSignalRemainingBGSubInvMassLeftPtBin;
    if (fHistoMappingRemainingBGInvMassPtBin)                   delete fHistoMappingRemainingBGInvMassPtBin;
    if (fHistoMappingRemainingBGInvMassLeftPtBin)               delete fHistoMappingRemainingBGInvMassLeftPtBin;
    if (fHistoMappingRatioSBInvMassPtBin)                       delete fHistoMappingRatioSBInvMassPtBin;
    if (fFitSignalInvMassPtBin)                                 delete fFitSignalInvMassPtBin;
    if (fFitRemainingBGInvMassPtBin)                            delete fFitRemainingBGInvMassPtBin;
    if (fFitRemainingBGInvMassLeftPtBin)                        delete fFitRemainingBGInvMassLeftPtBin;
    if (fFitSignalPeakPosInvMassPtBin)                          delete fFitSignalPeakPosInvMassPtBin;
    if (fFitBckInvMassPtBin)                                    delete fFitBckInvMassPtBin;
    if (fHistoMappingBackNormInvMassLeftPtBin)                  delete fHistoMappingBackNormInvMassLeftPtBin;
    if (fHistoMappingSignalInvMassLeftPtBin)                    delete fHistoMappingSignalInvMassLeftPtBin;
    if (fFitInvMassLeftPtBin)                                   delete fFitInvMassLeftPtBin;
    if (fFitSignalPeakPosInvMassLeftPtBin)                      delete fFitSignalPeakPosInvMassLeftPtBin;
    if (fFitBckInvMassLeftPtBin)                                delete fFitBckInvMassLeftPtBin;
    if (fFitRatioInvMassPtBin)                                  delete fFitRatioInvMassPtBin;
    if (fMesonFWHMAlpha01)                                      delete fMesonFWHMAlpha01;
    if (fMesonFWHMAlpha01Error)                                 delete fMesonFWHMAlpha01Error;
    if (fHistoWeightsBGZbinVsMbin)                              delete fHistoWeightsBGZbinVsMbin;
    if (fHistoFillPerEventBGZbinVsMbin)                         delete fHistoFillPerEventBGZbinVsMbin;
    if (fHistoWeightsBGZbinVsPsibin)                            delete fHistoWeightsBGZbinVsPsibin;
    if (fHistoFillPerEventBGZbinVsPsibin)                       delete fHistoFillPerEventBGZbinVsPsibin;
    // delete Gaussian fit histograms
    if (fMesonMassGaussian)                                     delete fMesonMassGaussian;
    if (fMesonMassGaussianError)                                delete fMesonMassGaussianError;
    if (fMesonWidthGaussian)                                    delete fMesonWidthGaussian;
    if (fMesonWidthGaussianError)                               delete fMesonWidthGaussianError;
    if (fMesonTrueMassGaussian)                                 delete fMesonTrueMassGaussian;
    if (fMesonTrueMassGaussianError)                            delete fMesonTrueMassGaussianError;
    if (fMesonTrueWidthGaussian)                                delete fMesonTrueWidthGaussian;
    if (fMesonTrueWidthGaussianError)                           delete fMesonTrueWidthGaussianError;
    if (fHistoMassGaussianMeson)                                delete fHistoMassGaussianMeson;
    if (fHistoTrueMassGaussianMeson)                            delete fHistoTrueMassGaussianMeson;
    if (fHistoWidthGaussianMeson)                               delete fHistoWidthGaussianMeson;
    if (fHistoTrueWidthGaussianMeson)                           delete fHistoTrueWidthGaussianMeson;
    if (fFitSignalGaussianInvMassPtBin)                         delete fFitSignalGaussianInvMassPtBin;
    if (fFitTrueSignalGaussianInvMassPtBin)                     delete fFitTrueSignalGaussianInvMassPtBin;
//     if( fMode == 4 ) {
//     if (fFitSignalInvMassPtBin2)                                delete fFitSignalInvMassPtBin2;
//     if (fFitSignalPeakPosInvMassPtBin2)                         delete fFitSignalPeakPosInvMassPtBin2;
//     if (fFitBckInvMassPtBin2)                                   delete fFitBckInvMassPtBin2;
//     if (fFitInvMassLeftPtBin2)                                  delete fFitInvMassLeftPtBin2;
//     if (fFitSignalPeakPosInvMassLeftPtBin2)                     delete fFitSignalPeakPosInvMassLeftPtBin2;
//     if (fFitBckInvMassLeftPtBin2)                               delete fFitBckInvMassLeftPtBin2;
//     }
}

