// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion
//This file is not supposed to be run on outputfiles of the GammaConv-Software before the 30th Sept 2010.

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TStyle.h>
#include <TString.h>
#include <TLatex.h>
// #include "TGaxis.h"
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
#include "../CommonHeaders/ConversionFunctions.h"
#include "ExtractSignalV2.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "../CommonHeaders/ExtractSignalPlotting.h"
#include "THnSparse.h"
using namespace std;

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

    fMode      = mode;
    fModeHeavy = mode; // added to store mode >100 info
    // modes:
    // 0    new output PCM-PCM
    // 1    new output PCM dalitz
    // 2    new output PCM-Calo
    // 3    new output Calo-Calo
    // 4    new output EMCAL-EMCAL
    // 5    new output PHOS-PHOS
    // 9    old output PCM-PCM
    // 12   new output DCal-DCal
    // 13   new output PCM-DCal

    // Heavy meson fix
    if(mode>=100) fMode -= 100;


    //******************************************************************************************************
    //***************************** Get selected file and main dir *****************************************
    //******************************************************************************************************
    TFile* f = new TFile(file.Data());
    TString autoDetectedMainDir;
    if(mode<100) autoDetectedMainDir = AutoDetectMainTList(mode,f);
    else         autoDetectedMainDir = AutoDetectMainTList(mode,f,meson); // heavy meson analysis
    if (autoDetectedMainDir.CompareTo("") == 0){
        cout << "ERROR: trying to read file, which is incompatible with mode selected (mode " << mode << ")" << endl;;
        return;
    }
    TList *TopDir                   = (TList*)f->Get(autoDetectedMainDir.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }
    cout << "Reading from TopDir \"" << autoDetectedMainDir << "\"" << endl;


    if(optionAdvancedMesonQA.Contains("AdvancedMesonQA")){fAdvancedMesonQA = kTRUE;}

    // adjusting Cutnumber
    fCutSelection = cutSelection;
    TString fCutSelectionRead = cutSelection;
    if (mode == 9){
        ReturnSeparatedCutNumber(cutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection);
        fEventCutSelection  = fGammaCutSelection(0,7);
        fGammaCutSelection  = fGammaCutSelection(7,fGammaCutSelection.Length()-7);
        cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
    } else {
        ReturnSeparatedCutNumberAdvanced(cutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
    }

    //****************************** Specification of collision system ************************************************
    fEnergyFlag         = optionEnergy;
    fPrefix             = meson;

    fPeriodFlag         = optionPeriod;
    fdirectphoton       = directphotonPlots;

    TString textProcess = ReturnMesonString (fPrefix);
    if(textProcess.CompareTo("") == 0 ){
        cout << "Meson unknown" << endl;
        return ;
    }

    fTextMeasurement    = Form("%s #rightarrow #gamma#gamma", textProcess.Data());
    fCollisionSystem    = ReturnFullCollisionsSystem(fEnergyFlag);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    fDetectionProcess   = ReturnFullTextReconstructionProcess(mode);
    //**************************** Determine Centrality *************************************************************
    centralityString        = GetCentralityString(fEventCutSelection);
    cout << "Centrality: " << centralityString.Data() << endl;
    if (centralityString.CompareTo("pp")!=0 && centralityString.CompareTo("0-100%") != 0){
        fCollisionSystem    = Form("%s %s", centralityString.Data(), fCollisionSystem.Data());
    }
    cout << "Collisions system: " << fCollisionSystem.Data() << endl;

    // ******************************* Adjust cutstrings if needed **************************************************
    TString fEventCutSelectionRead  = fEventCutSelection.Data();
    TString fGammaCutSelectionRead  = fGammaCutSelection.Data();
    TString fMesonCutSelectionRead  = fMesonCutSelection.Data();
    TString fClusterCutSelectionRead= fClusterCutSelection.Data();
    if (addSig) {
        if(directphotonPlots.CompareTo("directPhoton")==0){
            cout << "running added Signal for photons, be careful" << endl;
            cout << fEventCutSelectionRead.Data() << endl;
            fEventCutSelectionRead.Replace(GetEventRejectExtraSignalsCutPosition(),1,"2");
            cout << fEventCutSelectionRead.Data() << endl;
            if (mode==9 || mode==109)
                fCutSelectionRead       = Form("%s%s_%s", fEventCutSelectionRead.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
            else if (mode==0 || mode==100)
                fCutSelectionRead       = Form("%s_%s_%s",fEventCutSelectionRead.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
            if (mode==2 || mode==3 || mode==102 || mode==103)
                fCutSelectionRead       = Form("%s_%s_%s_%s",fEventCutSelectionRead.Data(), fGammaCutSelection.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
            cout << fCutSelectionRead.Data() << endl;
        } else {
            cout << "running added Signal" << endl;
            cout << fEventCutSelection.Data() << endl;
            fEventCutSelectionRead.Replace(GetEventRejectExtraSignalsCutPosition(),1,"2");
            cout << fEventCutSelectionRead.Data() << endl;
            if (mode==9 || mode==109)
                fCutSelectionRead       = Form("%s%s_%s", fEventCutSelectionRead.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
            else if (mode==0 || mode==100)
                fCutSelectionRead       = Form("%s_%s_%s",fEventCutSelectionRead.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
            if (mode==2 || mode==3 || mode==102 || mode==103)
                fCutSelectionRead       = Form("%s_%s_%s_%s",fEventCutSelectionRead.Data(), fGammaCutSelection.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
            cout << fCutSelectionRead.Data() << endl;
        }
    }

    if(optionUseMinBiasEff.CompareTo("MinBiasEffOnly")==0 && optionMC.CompareTo("kTRUE") == 0){
        cout << "calculating MinBias Eff" << endl;
        // take out the cent numbers
        fEventCutSelectionRead.Replace(GetEventCentralityMinCutPosition(),2,"00");
        // set to full MB cut in pPb
        if (fEventCutSelectionRead.BeginsWith("a")) fEventCutSelectionRead.Replace(0,1,"8");
        if (fEventCutSelectionRead.BeginsWith("c")) fEventCutSelectionRead.Replace(0,1,"8");
        if (fEventCutSelectionRead.BeginsWith("9")) fEventCutSelectionRead.Replace(0,1,"8");
        if (fEventCutSelectionRead.BeginsWith("e")) fEventCutSelectionRead.Replace(0,1,"8");
        cout << fEventCutSelectionRead.Data() << endl;
        if (mode==0 || mode==100)
            fCutSelectionRead       = Form("%s_%s_%s",fEventCutSelectionRead.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
        if (mode==2 || mode==3 || mode==102 || mode==103)
            fCutSelectionRead       = Form("%s_%s_%s_%s",fEventCutSelectionRead.Data(), fGammaCutSelection.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
        if (mode==4 || mode==5 || mode==104 || mode==105)
            fCutSelectionRead       = Form("%s_%s_%s",fEventCutSelectionRead.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
        cout << fCutSelectionRead.Data() << endl;
    }

    //******************************************************************************************************
    //***************************** Get Main folder for cut ************************************************
    //******************************************************************************************************
    TList *HistosGammaConversion       = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    // check if cutnumber available, otherwise adjust for pileup cut
    if (HistosGammaConversion == NULL){
        TString fEventCutSelectionPileUpRejection   = fEventCutSelection(5,1);
        cout << "cutnumber for PileUpRejection is: " << fEventCutSelectionPileUpRejection << endl;
        if( CutNumberToInteger(fEventCutSelectionPileUpRejection) > 1  && optionMC.CompareTo("kTRUE") == 0){
        cout << "changing PileUpCut for MC" << endl;
        cout << fEventCutSelectionRead.Data() << endl;
        fEventCutSelectionRead.Replace(GetEventRemovePileUpCutPosition(),1,"1");
        cout << fEventCutSelectionRead.Data() << endl;
        fClusterCutSelectionRead  = fClusterCutSelection;
        if (mode==0 || mode==100)
            fCutSelectionRead       = Form("%s_%s_%s",fEventCutSelectionRead.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
        if (mode==2 || mode==3 || mode==102 || mode==103 || mode==100)
            fCutSelectionRead       = Form("%s_%s_%s_%s",fEventCutSelectionRead.Data(), fGammaCutSelection.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
        if (mode==4 || mode==5 || mode==104 || mode==105)
            fCutSelectionRead       = Form("%s_%s_%s",fEventCutSelectionRead.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
        cout << fCutSelectionRead.Data() << endl;
        }
    }
    HistosGammaConversion       = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if (HistosGammaConversion == NULL){
        TString fEventCutSelectionBaseEvent     = fEventCutSelectionRead(GetEventSystemCutPosition(),1);
        cout << "cutnumber for basic event selection is: " << fEventCutSelectionBaseEvent << endl;
        if ( CutNumberToInteger(fEventCutSelectionBaseEvent) == 1 && optionMC.CompareTo("kTRUE") == 0){
            cout << "changing basic event cut for MC" << endl;
            cout << fEventCutSelectionRead.Data() << endl;
            fEventCutSelectionRead.Replace(GetEventSystemCutPosition(),1,"5");
            cout << fEventCutSelectionRead.Data() << endl;
            fClusterCutSelectionRead  = fClusterCutSelection;
            if (mode==0 || mode==100)
                fCutSelectionRead       = Form("%s_%s_%s",fEventCutSelectionRead.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
            if (mode==2 || mode==3 || mode==102 || mode==103)
                fCutSelectionRead       = Form("%s_%s_%s_%s",fEventCutSelectionRead.Data(), fGammaCutSelection.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
            if (mode==4 || mode==5 || mode==104 || mode==105)
                fCutSelectionRead       = Form("%s_%s_%s",fEventCutSelectionRead.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
            cout << fCutSelectionRead.Data() << endl;
        } else if ( CutNumberToInteger(fEventCutSelectionBaseEvent) == 5 && optionMC.CompareTo("kTRUE") == 0){
            cout << "changing basic event cut for MC" << endl;
            cout << fEventCutSelectionRead.Data() << endl;
            fEventCutSelectionRead.Replace(GetEventSystemCutPosition(),1,"1");
            cout << fEventCutSelectionRead.Data() << endl;
            fClusterCutSelectionRead  = fClusterCutSelection;
            if (mode==0 || mode==100)
                fCutSelectionRead       = Form("%s_%s_%s",fEventCutSelectionRead.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
            if (mode==2 || mode==3 || mode==102 || mode==103)
                fCutSelectionRead       = Form("%s_%s_%s_%s",fEventCutSelectionRead.Data(), fGammaCutSelection.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
            if (mode==4 || mode==5 || mode==104 || mode==105)
                fCutSelectionRead       = Form("%s_%s_%s",fEventCutSelectionRead.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
            cout << fCutSelectionRead.Data() << endl;
        }
    }

    // Prepend digit in case of heavy meson analysis
    if( mode >= 100) fEventCutSelectionRead.Prepend( Form("%d_",GetHeavyMesonDigit(meson)) );

    StyleSettingsThesis(Suffix);
    SetPlotStyle();

    TString outputDir   = Form("%s/%s/%s/ExtractSignal",cutSelection.Data(),optionEnergy.Data(),Suffix.Data());
    TString outputDirMon= Form("%s/%s/%s/ExtractSignal/Monitoring/",cutSelection.Data(),optionEnergy.Data(),Suffix.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec("mkdir -p "+outputDirMon);

    cout<<"Pictures are saved as "<< Suffix.Data()<< endl;
    fdate = ReturnDateString();


    //****************************** Choice of Fitting procedure ******************************************************
    if(optionCrystalBall.CompareTo("CrystalBall") == 0){// means we want to plot values for the pi0
        fCrysFitting        = 1;
        cout << "CrystalBall fit chosen ..." << endl;
    } else   {
        fCrysFitting        = 0;
        cout << "Gaussian fit chosen ..." << endl;
    }

    if(cutSelection.Length() == 0){
        cout<<"ERROR: Cut selection is not set, please do!"<<endl;
        return;
    }

    //***************************** Specification Data/MC ************************************************************
    if(optionMC.CompareTo("kTRUE") == 0){
        fIsMC               = 1;
        fPrefix2            = "MC";
    } else {
        fIsMC               = 0;
        fPrefix2            = "data";
    }
    cout << "line " << __LINE__ << endl;

    //***************************** Initialization of variables according to meson type ******************************
    if(meson.CompareTo("Pi0") == 0){
        cout << "entering Pi0" << endl;
        Initialize("Pi0",numberOfBins, triggerSet);
    } else if (meson.CompareTo("Eta") == 0) {
        Initialize("Eta",numberOfBins, triggerSet);
    } else if (meson.CompareTo("EtaPrime") == 0) {
        Initialize("EtaPrime",numberOfBins, triggerSet);
    } else if(meson.CompareTo("Pi0EtaBinning") == 0) {
        Initialize("Pi0EtaBinning",numberOfBins, triggerSet);
    } else   {
        cout<<"ERROR: First argument in the ExtractSignal(....) has to be either Pi0 or Eta or Pi0EtaBinning  or EtaPrime"<<endl;
        return;
    }

    // set global variables for rap and BG number
    TString rapidityRange;
    fYMaxMeson                          = ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
    fBackgroundMultNumber               = ReturnBackgroundMult(fMesonCutSelection);

    cout << "Integration window normal: "<< fMesonIntDeltaRange[0] << "\t" << fMesonIntDeltaRange[1] << endl;
    cout << "Integration window narrow: "<< fMesonIntDeltaRangeNarrow[0] << "\t" << fMesonIntDeltaRangeNarrow[1] << endl;
    cout << "Integration window wide: "<< fMesonIntDeltaRangeWide[0] << "\t" << fMesonIntDeltaRangeWide[1] << endl;

    cout << "line " << __LINE__ << endl;


    //************************* Start of Main routine ***************************************************************
    const char* fFileErrLogDatname = Form("%s/%s/%s_%s_FileErrLog%s_%s.dat",cutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelectionRead.Data());
    fFileErrLog.open(fFileErrLogDatname, ios::out);

    //******************************************************************************************************
    //***************************** Get Main folder for cut ************************************************
    //******************************************************************************************************
    HistosGammaConversion       = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    // couldn't find main folder for cut
    if(HistosGammaConversion == NULL){
        //******************************************************************************************************
        //Check whether MC would contain a different time cut as timing cuts aren't implemented for EMC in MC **
        //******************************************************************************************************
        if ( (mode == 2 || mode == 4 || mode == 12 || mode == 13 ) && optionMC.CompareTo("kTRUE") == 0 ){
            cout << "testing whether output with diff time cut exists in file for MC" << endl;
            fEventCutSelectionRead          = fEventCutSelection;
            fGammaCutSelectionRead          = fGammaCutSelection;
            fMesonCutSelectionRead          = fMesonCutSelection;
            fClusterCutSelectionRead        = fClusterCutSelection;
            // Alternative time cuts for EMC
            TString mostProbableTimeCuts[10]= {"5", "6", "0", "1", "2", "3", "4", "7", "8", "9"};
            Int_t i = 0;
            while (HistosGammaConversion == NULL && i < 10){
                // replace time cut
                fClusterCutSelection.Replace(GetClusterTimingCutPosition(fClusterCutSelectionRead),1,mostProbableTimeCuts[i].Data());
                fClusterCutSelectionRead= fClusterCutSelection;
                if (mode==2 || mode==13 || mode==102 || mode==113)
                    fCutSelectionRead       = Form("%s_%s_%s_%s",fEventCutSelection.Data(), fGammaCutSelection.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
                else
                    fCutSelectionRead       = Form("%s_%s_%s",fEventCutSelection.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
                cout << "testing cutnumber:" << fCutSelectionRead.Data() << endl;
                HistosGammaConversion   = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
                i++;
            }
            i = 0;
            while (HistosGammaConversion == NULL && i < 10){
                // replace time cut
                fEventCutSelection.Replace(GetEventSystemCutPosition(),1,"1");
                cout << fEventCutSelection.Data() << endl;
                fEventCutSelectionRead    = fEventCutSelection;
                fClusterCutSelection.Replace(GetClusterTimingCutPosition(fClusterCutSelectionRead),1,mostProbableTimeCuts[i].Data());
                fClusterCutSelectionRead= fClusterCutSelection;
                if (mode==2 || mode==13 || mode==102 || mode==113)
                    fCutSelectionRead       = Form("%s_%s_%s_%s",fEventCutSelection.Data(), fGammaCutSelection.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
                else
                    fCutSelectionRead       = Form("%s_%s_%s",fEventCutSelection.Data(), fClusterCutSelection.Data(), fMesonCutSelection.Data());
                cout << "testing cutnumber:" << fCutSelectionRead.Data() << endl;
                HistosGammaConversion   = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
                i++;
            }
            if (HistosGammaConversion == NULL){
                cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in File"<<endl;
                return;
            }
        //******************************************************************************************************
        //***** Abort processing due to missing cutnumber in files**********************************************
        //******************************************************************************************************
        } else {
            cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in File"<<endl;
            return;
        }
    }

    //*********************************************************************************************************
    //******************* Set MC histo names ******************************************************************
    //*********************************************************************************************************
    if (meson.CompareTo("Pi0") == 0 || meson.CompareTo("Pi0EtaBinning") == 0 ){
        SetCorrectMCHistogrammNames("Pi0");
    } else if (meson.CompareTo("Eta") == 0 ){
        SetCorrectMCHistogrammNames("Eta");
    } else if (meson.CompareTo("EtaPrime") == 0 ){
        SetCorrectMCHistogrammNames("EtaPrime");
    }

    TList *ESDContainer                 = (TList*) HistosGammaConversion->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));
    TList *BackgroundContainer          = (TList*) HistosGammaConversion->FindObject(Form("%s Back histograms",fCutSelectionRead.Data()));
    TList *MotherContainer              = (TList*) HistosGammaConversion->FindObject(Form("%s Mother histograms",fCutSelectionRead.Data()));
    if (fMode == 2 || fMode == 13 || fMode == 3 ){
        TList *ClusterContainer             = (TList*) HistosGammaConversion->FindObject(Form("%s Cluster Output",fCutSelectionRead.Data()));
        if (ClusterContainer){
            fHistoClustersPt                = (TH1D*)ClusterContainer->FindObject("ClusGamma_Pt");
            fHistoClustersE                 = (TH1D*)ClusterContainer->FindObject("ClusGamma_E");
            fHistoClustersOverlapHeadersPt  = (TH1D*)ClusterContainer->FindObject("ClusGammaOverlapHeaders_Pt");
            TH2F* fHistoTrue2DGammaDCClusPt = (TH2F*)ClusterContainer->FindObject(ObjectNameDCGammaClusPt.Data());
            if (fHistoTrue2DGammaDCClusPt!=NULL) fEnableDCCluster= kTRUE;
            if (fEnableDCCluster){
                fHistoTrue2DGammaDCClusPt->Sumw2();
                fHistoTrueGammaDCClusPt             = (TH1D*)fHistoTrue2DGammaDCClusPt->ProjectionX("TrueClusGamma_Pt",0,-1,"e");
                cout << "Cluster DC found " << endl;
                fHistoTrueGammaClusMultipleCount    = (TH1F*)ClusterContainer->FindObject(ObjectNameGammaClusMultipleCount.Data());
                fHistoTrueGammaClusPt               = (TH1F*)ClusterContainer->FindObject("TrueClusGamma_Pt");
            }
        }
    }
    if ( fMode == 4 || fMode == 12  || fMode == 5 ){
        fHistoClustersPt                = (TH1D*)ESDContainer->FindObject("ClusGamma_Pt");
        fHistoClustersE                 = (TH1D*)ESDContainer->FindObject("ClusGamma_E");
        fHistoClustersOverlapHeadersPt  = (TH1D*)ESDContainer->FindObject("ClusGammaOverlapHeaders_Pt");
    }
    if ( fMode == 0 ){
        fHistoClustersPt                = (TH1D*)ESDContainer->FindObject("ClusGamma_Pt");
        fHistoClustersE                 = (TH1D*)ESDContainer->FindObject("ClusGamma_E");
        if (fHistoClustersPt){
            cout << "INFO: found cluster output in PCM stream, adding it to the raw data file." << endl;
        }

    }

    TList* EventCuts                    = (TList*)HistosGammaConversion->FindObject(Form("ConvEventCuts_%s",fEventCutSelectionRead.Data()));
    if (EventCuts){
      fHistoPileUpVertexDistance               = (TH1F*) EventCuts->FindObject(Form("PileupVertexDistance %s",fEventCutSelectionRead.Data()));
      fHistoPileUpVertexDistance_SPDPileup     = (TH1F*) EventCuts->FindObject(Form("PileupVertexDistance_SPDPileup %s",fEventCutSelectionRead.Data()));
      fHistoPileUpVertexDistance_TrackletHits  = (TH1F*) EventCuts->FindObject(Form("PileupVertexDistance_TrackletvsHits %s",fEventCutSelectionRead.Data()));
    }

    cout << fMesonCutSelectionRead.Data() << endl;
    cout << fGammaCutSelectionRead.Data() << endl;
    fNumberOfGoodESDTracks              = (TH1D*)ESDContainer->FindObject("GoodESDTracks");
    fEventQuality                       = (TH1D*)ESDContainer->FindObject("NEvents");

    TString ObjectNameESD               = "ESD_Mother_InvMass_Pt";
    TString ObjectNameBck               = "ESD_Background_InvMass_Pt";

    fGammaGammaInvMassVSPt              = (TH2D*)ESDContainer->FindObject(ObjectNameESD.Data());
    fGammaGammaInvMassVSPt->Sumw2();
    fBckInvMassVSPt                     = (TH2D*)ESDContainer->FindObject(ObjectNameBck.Data());
    fBckInvMassVSPt->Sumw2();

    const char* FileDataLogname         = Form("%s/%s/%s_%s_EffiCheck_RAWDATA%s_%s.dat", cutSelection.Data(), fEnergyFlag.Data(), fPrefix.Data(), fPrefix2.Data(), fPeriodFlag.Data(),
                                        fCutSelectionRead.Data());
    fFileDataLog.open(FileDataLogname, ios::out);

    if(UseTHnSparse) ProduceBckProperWeighting(BackgroundContainer,MotherContainer,UseTHnSparse);
    else ProduceBckProperWeighting(ESDContainer,ESDContainer,UseTHnSparse);
    cout << "line " << __LINE__ << endl;

    // enter pure simulation routines
    if(fIsMC){
        // load containers for simulation
        TList *MCContainer              = (TList*)HistosGammaConversion->FindObject(Form("%s MC histograms",fCutSelectionRead.Data()));
        TList *TrueConversionContainer  = (TList*)HistosGammaConversion->FindObject(Form("%s True histograms",fCutSelectionRead.Data()));

        // loading histograms for pi0
        if( fMesonId == 111){
            // histos without acceptance requirement
            fHistoMCMesonPt                     = (TH1D*)MCContainer->FindObject(ObjectNameMCPi0.Data());   // Not the best; better having a 2D Pt_vs_Rapid in case we change limits
            fHistoMCMesonPtWOWeights            = (TH1D*)MCContainer->FindObject(ObjectNameMCPi0WOWeights.Data());
            fHistoMCMesonPtWOEvtWeights         = (TH1D*)MCContainer->FindObject(ObjectNameMCPi0WOEvtWeights.Data());

            // histos with gamma's in acceptance
            fHistoMCMesonPtWithinAcceptance     = (TH1D*)MCContainer->FindObject(ObjectNameMCPi0Acc.Data());
            // if ( fMode == 2 || fMode == 13 || fMode == 3 || fMode == 4 || fMode == 12 || fMode == 5 ){
                fHistoMCMesonPtWithinAcceptanceWOWeights    = (TH1D*)MCContainer->FindObject(ObjectNameMCPi0AccWOWeights.Data());
                fHistoMCMesonPtWithinAcceptanceWOEvtWeights = (TH1D*)MCContainer->FindObject(ObjectNameMCPi0AccWOEvtWeights.Data());
            // }

            // secondary neutral pions histograms
            fHistoMCSecPi0SourcePt              = (TH2D*)MCContainer->FindObject(ObjectNameMCSecPi0.Data());
            fHistoMCSecPi0WAccSourcePt          = (TH2D*)MCContainer->FindObject(ObjectNameMCSecPi0Acc.Data());
            if (fHistoMCSecPi0SourcePt && fHistoMCSecPi0WAccSourcePt){
                fHistoMCSecPi0SourcePt->Sumw2();
                fHistoMCSecPi0WAccSourcePt->Sumw2();
                fNewMCOutput                     =  kTRUE;
            }
            cout << "line " << __LINE__ << endl;
        }

        // Loading histograms for Eta
        if( fMesonId == 221){
            // Histograms without acceptance requirement
            fHistoMCMesonPt                     = (TH1D*)MCContainer->FindObject(ObjectNameMCEta.Data()); // (not the best; better having a 2D Pt_vs_Rapid in case we change limits)
            fHistoMCMesonPtWOWeights            = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaWOWeights.Data());
            fHistoMCMesonPtWOEvtWeights         = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaWOEvtWeights.Data());
            // Histograms with gammas in acceptance
            fHistoMCMesonPtWithinAcceptance     = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaAcc.Data());
            fHistoMCMesonPtWithinAcceptanceWOWeights    = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaAccWOWeights.Data());
            fHistoMCMesonPtWithinAcceptanceWOEvtWeights = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaAccWOEvtWeights.Data());
        }

        // Loading histograms for EtaPrime
        else if( fMesonId == 331 ){
            // Histograms without acceptance requirement
            fHistoMCMesonPt                     = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaPrime.Data()); // (not the best; better having a 2D Pt_vs_Rapid in case we change limits)
            fHistoMCMesonPtWOWeights            = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaPrimeWOWeights.Data());
            fHistoMCMesonPtWOEvtWeights         = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaPrimeWOEvtWeights.Data());
            // Histograms with gammas in acceptance
            fHistoMCMesonPtWithinAcceptance     = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaPrimeAcc.Data());
            fHistoMCMesonPtWithinAcceptanceWOWeights    = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaPrimeAccWOWeights.Data());
            fHistoMCMesonPtWithinAcceptanceWOEvtWeights = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaPrimeAccWOEvtWeights.Data());
        }

        // prepare histos for correct error calculation
        fHistoMCMesonPt->Sumw2();
        fHistoMCMesonPtWithinAcceptance->Sumw2();
        if (fHistoMCMesonPtWithinAcceptanceWOWeights) fHistoMCMesonPtWithinAcceptanceWOWeights->Sumw2();
        if (fHistoMCMesonPtWithinAcceptanceWOEvtWeights) fHistoMCMesonPtWithinAcceptanceWOEvtWeights->Sumw2();

        // calculate applied weights
        if (fHistoMCMesonPtWOWeights){
            fHistoMCMesonPtWeights              = (TH1D*)fHistoMCMesonPtWOWeights->Clone("WeightsMeson");
            fHistoMCMesonPtWeights->Divide(fHistoMCMesonPt,fHistoMCMesonPtWOWeights, 1.,1.,"B");
        }

        // load double counting histograms for Calo mode
        if (fMode == 4 || fMode == 12 || fMode == 5){
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

        // load reconstructed meson histograms
        fHistoTrueMesonInvMassVSPt                  = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrue.Data());
        fHistoTrueFullMesonInvMassVSPt              = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueFull.Data());
        fHistoTrueMesonInvMassVSPtWOWeights         = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueWOWeights.Data());
        fProfileTrueMesonInvMassVSPtWeights         = (TProfile2D*)TrueConversionContainer->FindObject(ObjectNameProfileWeights.Data());
        fHistoTrueMesonInvMassVSPtReweighted        = (TH2D*)fHistoTrueMesonInvMassVSPtWOWeights->Clone("Reweighted");
        fHistoTrueMesonInvMassVSPt->Sumw2();
        fHistoTrueFullMesonInvMassVSPt->Sumw2();
        fHistoTrueMesonInvMassVSPtWOWeights->Sumw2();
        fProfileTrueMesonInvMassVSPtWeights->Sumw2();
        fHistoTrueMesonInvMassVSPtReweighted->Sumw2();
        fHistoTrueMesonInvMassVSPtReweighted->Multiply(fProfileTrueMesonInvMassVSPtWeights);

        cout << ObjectNameTrue.Data() << endl;
        FillMassMCTrueMesonHistosArray(fHistoTrueMesonInvMassVSPt);
        FillMassMCTrueFullMesonHistosArray(fHistoTrueFullMesonInvMassVSPt);
        FillMassMCTrueReweightedMesonHistosArray(fHistoTrueMesonInvMassVSPtReweighted);
        FillMassMCTrueUnweightedMesonHistosArray(fHistoTrueMesonInvMassVSPtWOWeights);

        cout << "line " << __LINE__ << endl;
        cout << "Mode: " << fMode << endl;

        if (fMode == 2 || fMode == 13 || fMode == 3 || fMode == 4 || fMode == 12 || fMode == 5){
            fHistoTrueMesonCaloPhotonInvMassVSPt                = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloPhoton.Data());
            if (fHistoTrueMesonCaloPhotonInvMassVSPt==NULL) fAdvancedMesonQA = kFALSE;
            else fAdvancedMesonQA = kTRUE;
            cout << fAdvancedMesonQA << endl;
        }

        if (fMode == 0 || fMode == 1 || fMode == 9){
            fHistoTrueContBckInvMassVSPt                        = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueContBck.Data());
            if (fHistoTrueContBckInvMassVSPt==NULL) fAdvancedMesonQA = kFALSE;
        }

        if (fAdvancedMesonQA) {
            if (fMode == 2 || fMode == 13 || fMode == 3){
                fHistoTrueMesonCaloPhotonInvMassVSPt                = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloPhoton.Data());
                FillMassMCTrueMesonCaloPhotonHistosArray(fHistoTrueMesonCaloPhotonInvMassVSPt);
                fHistoTrueMesonCaloConvPhotonInvMassVSPt            = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloConvPhoton.Data());
                FillMassMCTrueMesonCaloConvPhotonHistosArray(fHistoTrueMesonCaloConvPhotonInvMassVSPt);
                fHistoTrueMesonMergedClusterInvMassVSPt               = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloMerged.Data());
                FillMassMCTrueMesonCaloMergedClusterHistosArray(fHistoTrueMesonMergedClusterInvMassVSPt);
                fHistoTrueMesonMergedClusterPartConvInvMassVSPt     = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueCaloMergedPartConv.Data());
                FillMassMCTrueMesonCaloMergedClusterPartConvHistosArray(fHistoTrueMesonMergedClusterPartConvInvMassVSPt);
            } else if (fMode == 4 || fMode == 12 || fMode == 5){
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
            }
        //    else {
        //        fHistoTrueContBckInvMassVSPt                        = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueContBck.Data());
        //        FillMassMCTrueContBckHistosArray(fHistoTrueContBckInvMassVSPt);
        //        fHistoTrueGGBckInvMassVSPt                          = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueGGBck.Data());
        //        FillMassMCTrueGGBckHistosArray(fHistoTrueGGBckInvMassVSPt);
        //        fHistoTrueAllBckInvMassVSPt                         = (TH2D*)fHistoTrueGGBckInvMassVSPt->Clone(ObjectNameTrueAllBck.Data());
        //        fHistoTrueAllBckInvMassVSPt->Sumw2();
        //        fHistoTrueAllBckInvMassVSPt->Add(fHistoTrueContBckInvMassVSPt);
        //        FillMassMCTrueAllBckHistosArray(fHistoTrueAllBckInvMassVSPt);
        //    }
            fHistoYieldK0sWithPi0DaughterRec                        = (TH1D*)TrueConversionContainer->FindObject(ObjectNameK0sRecPi0.Data());
            if(fHistoYieldK0sWithPi0DaughterRec) fHistoYieldK0sWithPi0DaughterRec->Sumw2();
            fHistoYieldLambdaWithPi0DaughterRec                     = (TH1D*)TrueConversionContainer->FindObject(ObjectNameLambdaRecPi0.Data());
            if(fHistoYieldLambdaWithPi0DaughterRec) fHistoYieldLambdaWithPi0DaughterRec->Sumw2();
        }

        fHistoTrueContBckInvMassVSPt                        = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueContBck.Data());
        fHistoTrueGGBckInvMassVSPt                          = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueGGBck.Data());
        if(fHistoTrueContBckInvMassVSPt && fHistoTrueGGBckInvMassVSPt){
            fEnableNormBckHistoComparisonToTrueBck = kTRUE;
            FillMassMCTrueContBckHistosArray(fHistoTrueContBckInvMassVSPt);
            FillMassMCTrueGGBckHistosArray(fHistoTrueGGBckInvMassVSPt);
            fHistoTrueAllBckInvMassVSPt                         = (TH2D*)fHistoTrueGGBckInvMassVSPt->Clone(ObjectNameTrueAllBck.Data());
            fHistoTrueAllBckInvMassVSPt->Sumw2();
            fHistoTrueAllBckInvMassVSPt->Add(fHistoTrueContBckInvMassVSPt);
            FillMassMCTrueAllBckHistosArray(fHistoTrueAllBckInvMassVSPt);
            fHistoTrueFullMesonContainedInvMassVSPt      = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueBckFullMesonContained.Data());
            fHistoTrueAsymEClusInvMassVSPt               = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueBckAsymEClus.Data());
            if(fHistoTrueFullMesonContainedInvMassVSPt) FillMassMCTrueFullMesonContainedHistosArray(fHistoTrueFullMesonContainedInvMassVSPt);
            if(fHistoTrueAsymEClusInvMassVSPt) FillMassMCTrueAsymEClusHistosArray(fHistoTrueAsymEClusInvMassVSPt);
        }


        if (meson.Contains("Pi0")){
            fHistoTrueSecMesonInvMassVSPt[0]                        = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueSecFromK0S.Data());
            fHistoTrueSecMesonInvMassVSPt[1]                        = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueSecFromLambda.Data());
            fHistoTrueSecMesonInvMassVSPt[2]                        = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueSecFromK0L.Data());
            fHistoTrueSecMesonInvMassVSPt[3]                        = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueSec.Data());

            for (Int_t j = 0; j<3; j++){
                if (fHistoTrueSecMesonInvMassVSPt[j]) fHistoTrueSecMesonInvMassVSPt[3]->Add(fHistoTrueSecMesonInvMassVSPt[j],-1);
            }
            FillMassMCTrueSecMesonHistosArray(fHistoTrueSecMesonInvMassVSPt);

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

    // calculate meson mass from pdg code
    cout << (TDatabasePDG::Instance())->GetParticle(fMesonId) << endl;
    fMesonMassExpect                            = (TDatabasePDG::Instance())->GetParticle(fMesonId)->Mass();
    // calculate number of events for normalization
    if (fEnergyFlag.Contains("PbPb") || fEnergyFlag.Contains("pPb")){
        fNEvents        = fEventQuality->GetBinContent(1);
    } else {
        fNEvents        = GetNEvents(fEventQuality);
    }
    cout<< "The mass of the meson is: "<< fMesonMassExpect<< " Events analysed: "<< fNEvents<< endl;

    // Function to Project the 2D histos InvariantMass VS Pt into Invariant Mass spectrum
    FillMassHistosArray(fGammaGammaInvMassVSPt);
    cout << "line " << __LINE__ << endl;

    ProcessEM( fMesonFullPtSignal, fMesonFullPtBackground, fBGFitRange);
    fMesonFullPtBackNorm                        = fBckNorm;

    ProcessEM( fFittingHistMidPtSignal, fFittingHistMidPtBackground, fBGFitRange);
    fFittingHistMidPtSignalSub                  = fSignal;
    if(fCrysFitting==0){
        fFileErrLog << "Using exp fit"<<endl;
        FitSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub, fMesonIntDeltaRange,200,kTRUE);
        fFitSignalInvMassMidPt                  = fFitReco;
    } else {
        fFileErrLog << "Using Crystal Ball function"<<endl;
        FitCBSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub, fMesonIntDeltaRange,200,kTRUE,"SinglefitfunctionMidPt",kFALSE);
        fFitSignalInvMassMidPt                  = fFitReco;
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
        if(fEnableNormBckHistoComparisonToTrueBck) fHistoMappingBackNormAndRemainingBGInvMassPtBin[iPt] = (TH1D*) fHistoMappingBackNormInvMassPtBin[iPt]->Clone(Form("histoBackNormAndRemainingBGBin%02d",iPt));

        fHistoMappingSignalInvMassPtBin[iPt]->SetName(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",iPt));

        // Fitting the subtracted spectra
        fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "normal range/right normalization" << endl;

        fFitSignalInvMassPtBin[iPt]             = 0x00;
        fMesonResidualBGlin[iPt]                = 0;
        fMesonResidualBGlinError[iPt]           = 0;
        fMesonResidualBGcon[iPt]                = 0;
        fMesonResidualBGconError[iPt]           = 0;
        if(fCrysFitting==0){
            fFileErrLog << "Using exp fit"<<endl;
            fFileDataLog << "Subtracted mixed event" << endl;
            FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
            fFitSignalInvMassPtBin[iPt]                 = fFitReco;
            fFitSignalPeakPosInvMassPtBin[iPt]          = fFitGausExp;
            fFitBckInvMassPtBin[iPt]                    = fFitLinearBck;
            fMesonYieldsResidualBckFunc[0][iPt]         = fIntLinearBck;
            fMesonYieldsResidualBckFuncError[0][iPt]    = fIntLinearBckError;
            fMesonResidualBGlin[iPt]                    = fFitBckInvMassPtBin[iPt]->GetParameter(1);
            fMesonResidualBGlinError[iPt]               = fFitBckInvMassPtBin[iPt]->GetParError(1);
            fMesonResidualBGcon[iPt]                    = fFitBckInvMassPtBin[iPt]->GetParameter(0);
            fMesonResidualBGconError[iPt]               = fFitBckInvMassPtBin[iPt]->GetParError(0);
            if (fFitReco)
                fMesonChi2[0][iPt]                      = fFitReco->GetChisquare()/fFitReco->GetNDF();
            else
                fMesonChi2[0][iPt]                      = -1;

            if(fEnableNormBckHistoComparisonToTrueBck){
                for (Int_t j = fHistoMappingBackNormAndRemainingBGInvMassPtBin[iPt]->GetXaxis()->FindBin(fMesonMassPlotRange[0]); j < fHistoMappingBackNormAndRemainingBGInvMassPtBin[iPt]->GetXaxis()->FindBin(fMesonMassPlotRange[1])+1; j++){
                    Double_t startBinEdge                                   = fHistoMappingBackNormAndRemainingBGInvMassPtBin[iPt]->GetXaxis()->GetBinLowEdge(j);
                    Double_t endBinEdge                                     = fHistoMappingBackNormAndRemainingBGInvMassPtBin[iPt]->GetXaxis()->GetBinUpEdge(j);
                    Double_t intLinearBack                                  = fFitBckInvMassPtBin[iPt]->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
                    fHistoMappingBackNormAndRemainingBGInvMassPtBin[iPt]->SetBinContent(j,intLinearBack+fHistoMappingBackNormAndRemainingBGInvMassPtBin[iPt]->GetBinContent(j));
                }
            }
            FitSubtractedPol2InvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
            fFitSignalWithOtherBGInvMassPtBin[0][iPt]   = fFitReco;
            fFitBckOtherInvMassPtBin[0][iPt]            = fFitLinearBck;
            fMesonYieldsResBckOtherFunc[0][iPt]         = fIntLinearBck;
            fMesonYieldsResBckOtherFuncError[0][iPt]    = fIntLinearBckError;
            if (fFitReco)
                fMesonChi2[1][iPt]                      = fFitReco->GetChisquare()/fFitReco->GetNDF();
            else
                fMesonChi2[1][iPt]                      = -1;
            FitSubtractedExp1InvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
            fFitSignalWithOtherBGInvMassPtBin[1][iPt]   = fFitReco;
            fFitBckOtherInvMassPtBin[1][iPt]            = fFitLinearBck;
            fMesonYieldsResBckOtherFunc[1][iPt]         = fIntLinearBck;
            fMesonYieldsResBckOtherFuncError[1][iPt]    = fIntLinearBckError;
            if (fFitReco)
                fMesonChi2[2][iPt]                      = fFitReco->GetChisquare()/fFitReco->GetNDF();
            else
                fMesonChi2[2][iPt]                      = -1;
            FitSubtractedExp2InvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
            fFitSignalWithOtherBGInvMassPtBin[2][iPt]   = fFitReco;
            fFitBckOtherInvMassPtBin[2][iPt]            = fFitLinearBck;
            fMesonYieldsResBckOtherFunc[2][iPt]         = fIntLinearBck;
            fMesonYieldsResBckOtherFuncError[2][iPt]    = fIntLinearBckError;
            if (fFitReco)
                fMesonChi2[3][iPt]                      = fFitReco->GetChisquare()/fFitReco->GetNDF();
            else
                fMesonChi2[3][iPt]                      = -1;
        } else {
            fFileErrLog << "Using Crystal Ball function"<<endl;
            FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncNormalBin%02d",iPt),kFALSE);
            fHistoMappingSignalRemainingBGSubInvMassPtBin[iPt]  = (TH1D*)fCopySignal->Clone(Form("histoSignalRemainingBGSubtractedBin%02d",iPt));
            fHistoMappingRemainingBGInvMassPtBin[iPt]           = (TH1D*)fCopyOnlyBG->Clone(Form("histoRemainingBGBin%02d",iPt));
            fFitSignalInvMassPtBin[iPt]                         = fFitReco;
            fFitRemainingBGInvMassPtBin[iPt]                    = fFitLinearBck;
            fFitSignalPeakPosInvMassPtBin[iPt]                  = fFitGausExp;
            fFitBckInvMassPtBin[iPt]                            = fFitLinearBck;
            fMesonYieldsResidualBckFunc[0][iPt]                 = 0;
            fMesonYieldsResidualBckFuncError[0][iPt]            = 0;

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

            fMesonAmplitudepar[iPt]         = fFitSignalInvMassPtBin[iPt]->GetParameter(0);
            fMesonAmplitudeparError[iPt]    = fFitSignalInvMassPtBin[iPt]->GetParError(0);
            fMesonSigmapar[iPt]             = fFitSignalInvMassPtBin[iPt]->GetParameter(2);
            fMesonSigmaparError[iPt]        = fFitSignalInvMassPtBin[iPt]->GetParError(2);

            fMesonCurIntRange[0][0]         = fMesonMass[iPt] + fMesonIntDeltaRange[0];
            fMesonCurIntRange[1][0]         = fMesonMass[iPt] + fMesonIntDeltaRangeWide[0];
            fMesonCurIntRange[2][0]         = fMesonMass[iPt] + fMesonIntDeltaRangeNarrow[0];
            fMesonCurIntRange[0][1]         = fMesonMass[iPt] + fMesonIntDeltaRange[1];
            fMesonCurIntRange[1][1]         = fMesonMass[iPt] + fMesonIntDeltaRangeWide[1];
            fMesonCurIntRange[2][1]         = fMesonMass[iPt] + fMesonIntDeltaRangeNarrow[1];

        } else {
            fMesonMass[iPt]                 = fMesonMassExpect;
            fMesonMassError[iPt]            = 0.;
            fMesonCurIntRange[0][0]         = fMesonMassExpect + fMesonIntDeltaRange[0];
            fMesonCurIntRange[1][0]         = fMesonMassExpect + fMesonIntDeltaRangeWide[0];
            fMesonCurIntRange[2][0]         = fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
            fMesonCurIntRange[0][1]         = fMesonMassExpect + fMesonIntDeltaRange[1];
            fMesonCurIntRange[1][1]         = fMesonMassExpect + fMesonIntDeltaRangeWide[1];
            fMesonCurIntRange[2][1]         = fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];
        }

        for (Int_t k = 0; k < 3; k++){
            fMassWindowHigh[k][iPt]        = fMesonCurIntRange[k][1];
            fMassWindowLow[k][iPt]         = fMesonCurIntRange[k][0];
        }

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
            for (Int_t k = 0; k < 3; k++){
                IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRange[k]);
                fGGYields[k][iPt]               = fYields;
                fGGYieldsError[k][iPt]          = fYieldsError;

                // Integrate the bck histo
                IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRange[k]);
                fBckYields[k][iPt]              = fYields;
                fBckYieldsError[k][iPt]         = fYieldsError;

                // Integrate the signal histo
                fFileDataLog<< endl <<"Signal histo "<< nameIntRange[k].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
                IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRange[k]);
                fMesonYields[k][iPt]            = fYields;
                fMesonYieldsError[k][iPt]       = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
            }
        } else {
            for (Int_t k = 0; k < 3; k++){
                IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRange[k]);
                fGGYields[k][iPt]               = fYields;
                fGGYieldsError[k][iPt]          = fYieldsError;

                if (k == 0){
                    fHistoMappingBackNormInvMassPtBin[iPt]->Sumw2();
                    fHistoMappingBackNormInvMassPtBin[iPt]->Add(fHistoMappingRemainingBGInvMassPtBin[iPt],1);
                }
                IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRange[k]);
                fBckYields[k][iPt]              = fYields;
                fBckYieldsError[k][iPt]         = fYieldsError;

                // Integrate the signal histo
                fFileDataLog<< endl <<"Signal histo "<< nameIntRange[k].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
                IntegrateHistoInvMassStream( fHistoMappingSignalRemainingBGSubInvMassPtBin[iPt], fMesonCurIntRange[k]);
                fMesonYields[k][iPt]            = fYields;
                fMesonYieldsError[k][iPt]       = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
            }
        }


        if(fIsMC){
            fFileDataLog<< endl <<"True histo normal range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            fFitTrueSignalInvMassPtBin[iPt]=0x00;
            if(fCrysFitting==0){
                fFileErrLog << "Using exp fit"<<endl;
                cout << "default" << endl;
                FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);
            } else {
                fFileErrLog << "Using Crystal Ball function"<<endl;
                FitCBSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncMCTrueBin%02d",iPt),kTRUE);
            }

            if (fHistoMappingTrueMesonInvMassPtBins[iPt]->GetEntries() !=0){
                fFitTrueSignalInvMassPtBin[iPt]     = fFitReco;
                if (fFitTrueSignalInvMassPtBin[iPt] != 0x00){
                    fMesonTrueMass[iPt]             = fFitTrueSignalInvMassPtBin[iPt]->GetParameter(1);
                    fMesonTrueMassError[iPt]        = fFitTrueSignalInvMassPtBin[iPt]->GetParError(1);
                    fMesonTrueSigmapar[iPt]         = fFitTrueSignalInvMassPtBin[iPt]->GetParameter(2);
                    fMesonTrueSigmaparError[iPt]    = fFitTrueSignalInvMassPtBin[iPt]->GetParError(2);

                    fMesonLambdaTailMCpar[iPt]      = fFitTrueSignalInvMassPtBin[iPt]->GetParameter(3);
                    fMesonLambdaTailMCparError[iPt] = fFitTrueSignalInvMassPtBin[iPt]->GetParError(3);

                    CalculateFWHM(fFitTrueSignalInvMassPtBin[iPt]);
                    fMesonTrueFWHM[iPt]             = fFWHMFunc;
                    fMesonTrueFWHMError[iPt]        = fFWHMFuncError;
                    fFileDataLog << "TrueFWHM \t" << fMesonTrueFWHM[iPt] << "\t +-" << fMesonTrueFWHMError[iPt] << endl;
                    fMesonTrueIntRange[0][0]        = fMesonTrueMass[iPt] + fMesonIntDeltaRange[0];
                    fMesonTrueIntRange[1][0]        = fMesonTrueMass[iPt] + fMesonIntDeltaRangeWide[0];
                    fMesonTrueIntRange[2][0]        = fMesonTrueMass[iPt] + fMesonIntDeltaRangeNarrow[0];
                    fMesonTrueIntRange[0][1]        = fMesonTrueMass[iPt] + fMesonIntDeltaRange[1] ;
                    fMesonTrueIntRange[1][1]        = fMesonTrueMass[iPt] + fMesonIntDeltaRangeWide[1];
                    fMesonTrueIntRange[2][1]        = fMesonTrueMass[iPt] + fMesonIntDeltaRangeNarrow[1];
                } else {
                    fMesonTrueMass[iPt]             = 0.;
                    fMesonTrueMassError[iPt]        = 1.;
                    fMesonTrueFWHM[iPt]             = 0.;
                    fMesonTrueFWHMError[iPt]        = 0.;
                    fMesonTrueIntRange[0][0]        = fMesonMassExpect + fMesonIntDeltaRange[0];
                    fMesonTrueIntRange[1][0]        = fMesonMassExpect + fMesonIntDeltaRangeWide[0];
                    fMesonTrueIntRange[2][0]        = fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
                    fMesonTrueIntRange[0][1]        = fMesonMassExpect + fMesonIntDeltaRange[1];
                    fMesonTrueIntRange[1][1]        = fMesonMassExpect + fMesonIntDeltaRangeWide[1];
                    fMesonTrueIntRange[2][1]        = fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];
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
                cout << "reweighted" << endl;
                if(fCrysFitting==0){
                    fFileErrLog << "Using exp fit"<<endl;
                    FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtReweightedBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);
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
                        fMesonTrueIntReweightedRange[0][0]      = fMesonTrueMassReweighted[iPt] + fMesonIntDeltaRange[0];
                        fMesonTrueIntReweightedRange[1][0]      = fMesonTrueMassReweighted[iPt] + fMesonIntDeltaRangeWide[0];
                        fMesonTrueIntReweightedRange[2][0]      = fMesonTrueMassReweighted[iPt] + fMesonIntDeltaRangeNarrow[0];
                        fMesonTrueIntReweightedRange[0][1]      = fMesonTrueMassReweighted[iPt] + fMesonIntDeltaRange[1];
                        fMesonTrueIntReweightedRange[1][1]      = fMesonTrueMassReweighted[iPt] + fMesonIntDeltaRangeWide[1];
                        fMesonTrueIntReweightedRange[2][1]      = fMesonTrueMassReweighted[iPt] + fMesonIntDeltaRangeNarrow[1];
                    } else {
                        fMesonTrueMassReweighted[iPt]           = 0.;
                        fMesonTrueMassReweightedError[iPt]      = 1.;
                        fMesonTrueFWHMReweighted[iPt]           = 0.;
                        fMesonTrueFWHMReweightedError[iPt]      = 0.;
                        fMesonTrueIntReweightedRange[0][0]      = fMesonMassExpect + fMesonIntDeltaRange[0];
                        fMesonTrueIntReweightedRange[1][0]      = fMesonMassExpect + fMesonIntDeltaRangeWide[0];
                        fMesonTrueIntReweightedRange[2][0]      = fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
                        fMesonTrueIntReweightedRange[0][1]      = fMesonMassExpect + fMesonIntDeltaRange[1];
                        fMesonTrueIntReweightedRange[1][1]      = fMesonMassExpect + fMesonIntDeltaRangeWide[1];
                        fMesonTrueIntReweightedRange[2][1]      = fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];

                    }
                }
            }

            fFitTrueSignalInvMassPtUnweightedBin[iPt]=0x00;
            cout << "line " << __LINE__ << endl;
            if(fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt]){
                cout << "Using exp fit"<<endl;
                cout << "unweighted" << endl;
                fFileErrLog << "Using exp fit"<<endl;
                if(fCrysFitting==0){
                    fFileErrLog << "Using exp fit"<<endl;
                    FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);
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
                        fMesonTrueIntUnweightedRange[0][0]      = fMesonTrueMassUnweighted[iPt] + fMesonIntDeltaRange[0];
                        fMesonTrueIntUnweightedRange[1][0]      = fMesonTrueMassUnweighted[iPt] + fMesonIntDeltaRangeWide[0];
                        fMesonTrueIntUnweightedRange[2][0]      = fMesonTrueMassUnweighted[iPt] + fMesonIntDeltaRangeNarrow[0];
                        fMesonTrueIntUnweightedRange[0][1]      = fMesonTrueMassUnweighted[iPt] + fMesonIntDeltaRange[1];
                        fMesonTrueIntUnweightedRange[1][1]      = fMesonTrueMassUnweighted[iPt] + fMesonIntDeltaRangeWide[1];
                        fMesonTrueIntUnweightedRange[2][1]      = fMesonTrueMassUnweighted[iPt] + fMesonIntDeltaRangeNarrow[1];
                    } else {
                        fMesonTrueMassUnweighted[iPt]           = 0.;
                        fMesonTrueMassUnweightedError[iPt]      = 1.;
                        fMesonTrueFWHMUnweighted[iPt]           = 0.;
                        fMesonTrueFWHMUnweightedError[iPt]      = 0.;
                        fMesonTrueIntUnweightedRange[0][0]      = fMesonMassExpect + fMesonIntDeltaRange[0];
                        fMesonTrueIntUnweightedRange[1][0]      = fMesonMassExpect + fMesonIntDeltaRangeWide[0];
                        fMesonTrueIntUnweightedRange[2][0]      = fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
                        fMesonTrueIntUnweightedRange[0][1]      = fMesonMassExpect + fMesonIntDeltaRange[1];
                        fMesonTrueIntUnweightedRange[1][1]      = fMesonMassExpect + fMesonIntDeltaRangeWide[1];
                        fMesonTrueIntUnweightedRange[2][1]      = fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];

                    }
                }
            }

            if (fAdvancedMesonQA && (fMode == 2 || fMode == 13 || fMode == 3 || fMode == 4 || fMode == 12 || fMode == 5)){
                if (fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]->GetEntries() !=0){
                    if(fCrysFitting==0){
                        fFileErrLog << "Using exp fit"<<endl;
                        FitTrueInvMassInPtBins(fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);
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
                    if (fMode == 4 || fMode == 12 || fMode ==5 ){
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
            if (fAdvancedMesonQA && (fMode == 4 || fMode == 12 || fMode == 5)){
                if (fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[iPt]->GetEntries() !=0){
                    if(fCrysFitting==0){
                        fFileErrLog << "Using exp fit"<<endl;
                        FitTrueInvMassInPtBins(fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);
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

            // fill pt array with integrated validated MC yield for different integration windows: normal, wide, narrow
            for (Int_t k = 0; k< 3;k++){
                if (k > 0 )
                    fFileDataLog<< endl <<"True histo " << nameIntRange[k].Data() << " range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonTrueIntRange[k]);
                fMesonTrueYields[k][iPt]                        = fYields;
                fMesonTrueYieldsError[k][iPt]                   = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

                fFileDataLog<< endl <<"True histo " << nameIntRange[k].Data() << " range reweighted" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtReweightedBins[iPt], fMesonTrueIntReweightedRange[k]);
                fMesonTrueYieldsReweighted[k][iPt]              = fYields;
                fMesonTrueYieldsReweightedError[k][iPt]         = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

                fFileDataLog<< endl <<"True histo " << nameIntRange[k].Data() << " range unweighted" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt], fMesonTrueIntUnweightedRange[k]);
                fMesonTrueYieldsUnweighted[k][iPt]              = fYields;
                fMesonTrueYieldsUnweightedError[k][iPt]         = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

            }

            // fill double counting pt-array  if available
            if (fEnableDCMeson){
                IntegrateHistoInvMassStream( fHistoMappingTrueMesonDCInvMassPtBins[iPt], fMesonTrueIntRange[0]);
                fMesonTrueYieldsDC[iPt]                         = fYields;
                fMesonTrueYieldsDCError[iPt]                    = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
                cout << "***********************************************" << endl;
                cout << "***********************************************" << endl;
                cout << "***********************************************" << endl;
                cout << "***********************************************" << endl;
            }

            // fill secondary yield pt-arrays for neutral pion for different sources and integration windows
            if (meson.Contains("Pi0")){
                for (Int_t j = 0; j < 4; j++){
                    if (fHistoMappingTrueSecMesonInvMassPtBins[j][iPt]){
                        for (Int_t k = 0; k < 3; k++){
                            fFileDataLog<< endl <<"TrueSec " << nameSecondaries[j].Data() << " histo " << nameIntRange[k].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                            IntegrateHistoInvMassStream( fHistoMappingTrueSecMesonInvMassPtBins[j][iPt], fMesonTrueIntRange[k]);
                            fMesonTrueSecYields[k][j][iPt]              = fYields;
                            fMesonTrueSecYieldsError[k][j][iPt]         = fYieldsError;
                            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;
                        }
                    } else {
                        for (Int_t k= 0; k< 3; k++){
                            fMesonTrueSecYields[k][j][iPt]          = 0;
                            fMesonTrueSecYieldsError[k][j][iPt]     = 0;
                        }
                    }
                }
            } else {
                for (Int_t j = 0; j < 4; j++){
                    for (Int_t k= 0; k< 3; k++){
                        fMesonTrueSecYields[k][j][iPt]              = 0;
                        fMesonTrueSecYieldsError[k][j][iPt]         = 0;
                    }
                }
            }


            if( (fGGYields[0][iPt] - fMesonTrueYields[0][iPt]) > 0) {
                fMesonTrueSB[iPt]               = fMesonTrueYields[0][iPt] / ( fGGYields[0][iPt] - fMesonTrueYields[0][iPt] );
                fMesonTrueSign[iPt]             = fMesonTrueYields[0][iPt] / pow( ( fGGYields[0][iPt] - fMesonTrueYields[0][iPt] ) , 0.5);
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

        for (Int_t k = 0; k < 3; k++){
            if (k > 0){
                fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << nameIntRange[k].Data()<<  endl;
                Double_t intRange[2]    = {0,0};
                if (k == 1){
                    intRange[0]         = fMesonIntDeltaRangeWide[0];
                    intRange[1]         = fMesonIntDeltaRangeWide[1];
                } else if ( k == 2){
                    intRange[0]         = fMesonIntDeltaRangeNarrow[0];
                    intRange[1]         = fMesonIntDeltaRangeNarrow[1];
                }
                if(fCrysFitting==0){
                    fFileErrLog << "Using exp fit"<<endl;
                    FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], intRange, iPt, kFALSE);
                    fMesonYieldsResidualBckFunc[k][iPt]         = fIntLinearBck;
                    fMesonYieldsResidualBckFuncError[k][iPt]    = fIntLinearBckError;
                } else {
                    fMesonYieldsResidualBckFunc[k][iPt]         = 0;
                    fMesonYieldsResidualBckFuncError[k][iPt]    = 0;
                }
            }

            fFileDataLog << "Residual Background leftover " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                         << fMesonYieldsResidualBckFunc[k][iPt] << "\t +- \t" << fMesonYieldsResidualBckFuncError[k][iPt] << endl<< endl;

            /////////////////////// added to check yields //////////////////////////////////////////////////////////
            fTotalBckYields[k][iPt]                         = fBckYields[k][iPt] + fMesonYieldsResidualBckFunc[k][iPt];
            fTotalBckYieldsError[k][iPt]                    = pow(fBckYieldsError[k][iPt]*fBckYieldsError[k][iPt] + fMesonYieldsResidualBckFuncError[k][iPt]*fMesonYieldsResidualBckFuncError[k][iPt],0.5);
            fFileDataLog << "Total Background " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                         << fTotalBckYields[k][iPt] << "\t +- \t" << fTotalBckYieldsError[k][iPt] << endl<< endl;
            fFileDataLog << "Background " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                         << fBckYields[k][iPt] << "\t +- \t" << fBckYieldsError[k][iPt] << endl<< endl;
            ///////////////////////////////////////////////////////////////////////////////////////////////////////

            fMesonYieldsCorResidualBckFunc[k][iPt]          = fMesonYields[k][iPt]- fMesonYieldsResidualBckFunc[k][iPt];
            fMesonYieldsCorResidualBckFuncError[k][iPt]     = pow(( fMesonYieldsError[k][iPt]*fMesonYieldsError[k][iPt]
                                                                    + fMesonYieldsResidualBckFuncError[k][iPt]*fMesonYieldsResidualBckFuncError[k][iPt]),0.5);
            fMesonYieldsPerEvent[k][iPt]                    = fMesonYieldsCorResidualBckFunc[k][iPt]/fNEvents;
            fMesonYieldsPerEventError[k][iPt]               = fMesonYieldsCorResidualBckFuncError[k][iPt]/fNEvents;

            //Integrate Fit Function
            IntegrateFitFunc( fFitSignalPeakPosInvMassPtBin[iPt], fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRange[k]);
            fMesonYieldsFunc[k][iPt]                        = fYieldsFunc;

            //SB default
            if (fTotalBckYields[k][iPt] != 0){
                fMesonSBdefault[k][iPt]                     = fMesonYieldsCorResidualBckFunc[k][iPt]/fTotalBckYields[k][iPt];
                fMesonSBdefaultError[k][iPt]                = pow( pow(fMesonYieldsCorResidualBckFuncError[k][iPt]/fTotalBckYields[k][iPt], 2.) +
                                                            pow((fTotalBckYieldsError[k][iPt]*fMesonYieldsCorResidualBckFunc[k][iPt])/(fTotalBckYields[k][iPt] *fTotalBckYields[k][iPt]), 2.), 0.5);
            } else {
                fMesonSBdefault[k][iPt]         = 0.;
                fMesonSBdefaultError[k][iPt]    = 0.;
            }
            //Significance default
            if ( pow(fMesonYieldsCorResidualBckFunc[k][iPt] + fTotalBckYields[k][iPt],0.5) != 0){
                fMesonSigndefault[k][iPt]                   = fMesonYieldsCorResidualBckFunc[k][iPt]/pow(fMesonYieldsCorResidualBckFunc[k][iPt] + fTotalBckYields[k][iPt],0.5);
                Double_t a                                  = ( pow(fMesonYieldsCorResidualBckFunc[k][iPt] + fTotalBckYields[k][iPt], -0.5) -
                                                                0.5*fMesonYieldsCorResidualBckFunc[k][iPt]*pow(fMesonYieldsCorResidualBckFunc[k][iPt] +
                                                                fTotalBckYields[k][iPt], -1.5) * fMesonYieldsCorResidualBckFuncError[k][iPt]);
                Double_t b                                  = 0.5*fMesonYieldsCorResidualBckFunc[k][iPt]*pow(fMesonYieldsCorResidualBckFunc[k][iPt]
                                                              + fTotalBckYields[k][iPt],-1.5) * fTotalBckYieldsError[k][iPt];
                fMesonSigndefaultError[k][iPt]              = pow( a*a + b*b, 0.5);
            } else {
                fMesonSigndefault[k][iPt]       = 0.;
                fMesonSigndefaultError[k][iPt]  = 0.;
            }
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
            fFileErrLog << "Using exp fit"<<endl;
            FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
            fFitInvMassLeftPtBin[iPt]                               = fFitReco;
            fFitSignalPeakPosInvMassLeftPtBin[iPt]                  = fFitGausExp;
            fFitBckInvMassLeftPtBin[iPt]                            = fFitLinearBck;
            fMesonYieldsResidualBckFunc[3][iPt]                     = fIntLinearBck;
            fMesonYieldsResidualBckFuncError[3][iPt]                = fIntLinearBckError;
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
            fMesonYieldsResidualBckFunc[3][iPt]                     = 0;
            fMesonYieldsResidualBckFuncError[3][iPt]                = 0;

        }
        CalculateFWHM(fFitInvMassLeftPtBin[iPt]);
        fMesonFWHMLeft[iPt]         = fFWHMFunc;
        fMesonFWHMLeftError[iPt]    = fFWHMFuncError;

        if (fFitInvMassLeftPtBin[iPt] !=0x00){
            fMesonMassLeft[iPt]             = fFitInvMassLeftPtBin[iPt]->GetParameter(1);
            fMesonMassLeftError[iPt]        = fFitInvMassLeftPtBin[iPt]->GetParError(1);
            fMesonCurIntRange[3][0]         = fMesonMassLeft[iPt] + fMesonIntDeltaRange[0];
            fMesonCurIntRange[4][0]         = fMesonMassLeft[iPt] + fMesonIntDeltaRangeWide[0];
            fMesonCurIntRange[5][0]         = fMesonMassLeft[iPt] + fMesonIntDeltaRangeNarrow[0];
            fMesonCurIntRange[3][1]         = fMesonMassLeft[iPt] + fMesonIntDeltaRange[1];
            fMesonCurIntRange[4][1]         = fMesonMassLeft[iPt] + fMesonIntDeltaRangeWide[1];
            fMesonCurIntRange[5][1]         = fMesonMassLeft[iPt] + fMesonIntDeltaRangeNarrow[1];
        } else {
            fMesonMassLeft[iPt]             = 0.;
            fMesonMassLeftError[iPt]        = 0.;
            fMesonCurIntRange[3][0]         = fMesonMassExpect + fMesonIntDeltaRange[0];
            fMesonCurIntRange[4][0]         = fMesonMassExpect + fMesonIntDeltaRangeWide[0];
            fMesonCurIntRange[5][0]         = fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
            fMesonCurIntRange[3][1]         = fMesonMassExpect + fMesonIntDeltaRange[1];
            fMesonCurIntRange[4][1]         = fMesonMassExpect + fMesonIntDeltaRangeWide[1];
            fMesonCurIntRange[5][1]         = fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];
        }

        // Integrate the bck histo
        if (fCrysFitting == 0){
            for (Int_t k = 0; k < 3; k++){
                IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurIntRange[k+3]);
                fBckYields[k+3][iPt]              = fYields;
                fBckYieldsError[k+3][iPt]         = fYieldsError;

                fFileDataLog<< endl <<"Signal histo " << nameIntRange[k+3].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurIntRange[k+3]);
                fMesonYields[k+3][iPt]            = fYields;
                fMesonYieldsError[k+3][iPt]       = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;

            }
        } else {
            for (Int_t k = 0; k < 3; k++){
                if (k == 0){
                    fHistoMappingBackNormInvMassLeftPtBin[iPt]->Sumw2();
                    fHistoMappingBackNormInvMassLeftPtBin[iPt]->Add(fHistoMappingRemainingBGInvMassLeftPtBin[iPt],1);
                }
                IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurIntRange[k+3]);
                fBckYields[k+3][iPt]              = fYields;
                fBckYieldsError[k+3][iPt]         = fYieldsError;

                fFileDataLog<< endl <<"Signal histo " << nameIntRange[k+3].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingSignalRemainingBGSubInvMassLeftPtBin[iPt], fMesonCurIntRange[k+3]);
                fMesonYields[k+3][iPt]            = fYields;
                fMesonYieldsError[k+3][iPt]       = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;

            }
        }

        for (Int_t k = 0; k < 3; k++){
            if (k > 0){
                fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << nameIntRange[k+3].Data()<<  endl;
                Double_t intRange[2]    = {0,0};
                if (k == 1){
                    intRange[0]         = fMesonIntDeltaRangeWide[0];
                    intRange[1]         = fMesonIntDeltaRangeWide[1];
                } else if ( k == 2){
                    intRange[0]         = fMesonIntDeltaRangeNarrow[0];
                    intRange[1]         = fMesonIntDeltaRangeNarrow[1];
                }
                if(fCrysFitting==0){
                    fFileErrLog << "Using exp fit"<<endl;
                    FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], intRange, iPt, kFALSE);
                    fMesonYieldsResidualBckFunc[k+3][iPt]         = fIntLinearBck;
                    fMesonYieldsResidualBckFuncError[k+3][iPt]    = fIntLinearBckError;
                } else {
                    fMesonYieldsResidualBckFunc[k+3][iPt]         = 0;
                    fMesonYieldsResidualBckFuncError[k+3][iPt]    = 0;
                }
            }

            fFileDataLog << "Residual Background leftover " << nameIntRange[k+3].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                         << fMesonYieldsResidualBckFunc[k+3][iPt] << "\t +- \t" << fMesonYieldsResidualBckFuncError[k+3][iPt] << endl<< endl;

            fTotalBckYields[k+3][iPt]                       = fBckYields[k+3][iPt] + fMesonYieldsResidualBckFunc[k+3][iPt];
            fTotalBckYieldsError[k+3][iPt]                  = pow(fBckYieldsError[k+3][iPt]*fBckYieldsError[k+3][iPt] + fMesonYieldsResidualBckFuncError[k+3][iPt]*fMesonYieldsResidualBckFuncError[k+3][iPt],0.5);
            fFileDataLog << "Total Background " << nameIntRange[k+3].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                         << fTotalBckYields[k+3][iPt] << "\t +- \t" << fTotalBckYieldsError[k+3][iPt] << endl<< endl;
            fFileDataLog << "Background " << nameIntRange[k+3].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                         << fBckYields[k+3][iPt] << "\t +- \t" << fBckYieldsError[k+3][iPt] << endl<< endl;



            fMesonYieldsCorResidualBckFunc[k+3][iPt]        = fMesonYields[k+3][iPt]- fMesonYieldsResidualBckFunc[k+3][iPt];
            fMesonYieldsCorResidualBckFuncError[k+3][iPt]   = pow(( fMesonYieldsError[k+3][iPt]*fMesonYieldsError[k+3][iPt]
                                                                    + fMesonYieldsResidualBckFuncError[k+3][iPt]*fMesonYieldsResidualBckFuncError[k+3][iPt]),0.5);
            fMesonYieldsPerEvent[k+3][iPt]                  = fMesonYieldsCorResidualBckFunc[k+3][iPt]/fNEvents;
            fMesonYieldsPerEventError[k+3][iPt]             = fMesonYieldsCorResidualBckFuncError[k+3][iPt]/fNEvents;

            //Integrate Fit Function
            IntegrateFitFunc( fFitSignalPeakPosInvMassLeftPtBin[iPt], fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurIntRange[k+3]);
            fMesonYieldsFunc[k+3][iPt]                      = fYieldsFunc;
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
        fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
            fGGYields[0][iPt] << "+-" << fGGYieldsError[0][iPt] << "\t" <<
            fGGYields[1][iPt] << "+-" << fGGYieldsError[1][iPt] << "\t" <<
            fGGYields[2][iPt] << "+-" << fGGYieldsError[2][iPt] << endl;

    }
    fileSysErrDat <<  endl;
    fileSysErrDat << "fTotalBckYields" << endl;
    fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
            fTotalBckYields[0][iPt] << "+-" << fTotalBckYieldsError[0][iPt] << "\t" <<
            fTotalBckYields[1][iPt] << "+-" << fTotalBckYieldsError[1][iPt] << "\t" <<
            fTotalBckYields[2][iPt] << "+-" << fTotalBckYieldsError[2][iPt] << "\t" <<
            fTotalBckYields[3][iPt] << "+-" << fTotalBckYieldsError[3][iPt]<< "\t" <<
            fTotalBckYields[4][iPt] << "+-" << fTotalBckYieldsError[4][iPt]<< "\t" <<
            fTotalBckYields[5][iPt] << "+-" << fTotalBckYieldsError[5][iPt] << endl;
    }
    fileSysErrDat <<  endl;
    fileSysErrDat << "fMesonYields" << endl;
    fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
            fMesonYieldsCorResidualBckFunc[0][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError[0][iPt] << "\t" <<
            fMesonYieldsCorResidualBckFunc[1][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError[1][iPt] << "\t" <<
            fMesonYieldsCorResidualBckFunc[2][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError[2][iPt] << "\t" <<
            fMesonYieldsCorResidualBckFunc[3][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError[3][iPt]<< "\t" <<
            fMesonYieldsCorResidualBckFunc[4][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError[4][iPt]<< "\t" <<
            fMesonYieldsCorResidualBckFunc[5][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError[5][iPt] << endl;
    }
    if(fIsMC){
        fileSysErrDat <<  endl;
        fileSysErrDat << "TrueYields" << endl;
        fileSysErrDat << "Bin \t True \t True Wide \t True Narr " << endl;
        for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
            fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
                fMesonTrueYields[0][iPt] << "\t" <<
                fMesonTrueYields[1][iPt] << "\t" <<
                fMesonTrueYields[2][iPt] << endl;
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
    TString labelsOtherFits[3]  = {"pol2 BG", "a exp(bx) BG", "a + b exp(cx) BG"};
    if (fCrysFitting == 0){
        nameMesonSub    = Form("%s_MesonSubtracted%s", plotPrefix.Data(), plotSuffix.Data());
        nameCanvasSub   = "MesonCanvasSubtracted";
        namePadSub      = "MesonPadSubtracted";
        cout << nameMesonSub.Data() << endl;
        PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassPtBin, nameMesonSub, nameCanvasSub, namePadSub,
                                            fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess,
                                            fCollisionSystem, "MC validated");
        nameMesonSub    = Form("%s_MesonSubtractedWithOther%s", plotPrefix.Data(), plotSuffix.Data());
        nameCanvasSub   = "MesonCanvasSubtractedWithOther";
        namePadSub      = "MesonPadSubtractedWithOther";
        cout << nameMesonSub.Data() << endl;

        PlotWithManyFitSubtractedInvMassInPtBins(   fHistoMappingSignalInvMassPtBin, fFitSignalInvMassPtBin, fFitSignalWithOtherBGInvMassPtBin, 3, labelsOtherFits, nameMesonSub,
                                                    nameCanvasSub, namePadSub, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt,
                                                    fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem, "MC validated", kTRUE, "pol1 BG");

        nameMesonSub    = Form("%s_MesonSubtractedWithOtherOnlyBGFits%s", plotPrefix.Data(), plotSuffix.Data());
        nameCanvasSub   = "MesonCanvasSubtractedWithOtherOnlyBGFits";
        namePadSub      = "MesonPadSubtractedWithOtherOnlyBGFits";
        cout << nameMesonSub.Data() << endl;

        PlotWithManyFitSubtractedInvMassInPtBins(   fHistoMappingSignalInvMassPtBin, fFitBckInvMassPtBin, fFitBckOtherInvMassPtBin, 3, labelsOtherFits, nameMesonSub,
                                                    nameCanvasSub, namePadSub, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt,
                                                    fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem, "MC validated", kTRUE, "pol1 BG");

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
        PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassLeftPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitInvMassLeftPtBin, nameMesonSub, nameCanvasSub, namePadSub,
                                            fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess,
                                            fCollisionSystem, "MC validated");

        cout << "Example bin: " << fExampleBin << endl;
        TString triggerInt         = fEventCutSelectionRead(GetEventSelectSpecialTriggerCutPosition(),2);
        PlotExampleInvMassBinsV2(fHistoMappingGGInvMassPtBin[fExampleBin], fHistoMappingSignalInvMassPtBin[fExampleBin], fHistoMappingBackNormInvMassPtBin[fExampleBin],
                            fFitSignalInvMassPtBin[fExampleBin], fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassPlotRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                            fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, triggerInt.Atoi(), fExampleBinScaleFac, fMode, addSig );

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

        PlotExampleInvMassBinsV2(fHistoMappingGGInvMassPtBin[fExampleBin], fHistoMappingSignalRemainingBGSubInvMassPtBin[fExampleBin], fHistoMappingBackNormInvMassPtBin[fExampleBin],
                            fFitSignalInvMassPtBin[fExampleBin], fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassPlotRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                            fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, triggerSet, fExampleBinScaleFac, fMode);
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
            PlotInvMassSecondaryInPtBins( fHistoMappingTrueMesonInvMassPtBins, fHistoMappingTrueSecMesonInvMassPtBins[3], fHistoMappingTrueSecMesonInvMassPtBins[1], nameMesonTrue, nameCanvasTrue,
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

        if (fAdvancedMesonQA && (fMode == 2 || fMode == 13 || fMode == 3 || fMode == 4 || fMode == 12 || fMode == 5)){

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
        if (fAdvancedMesonQA && (fMode == 4 || fMode == 12 || fMode == 5)){
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
            TString triggerInt         = fEventCutSelectionRead(GetEventSelectSpecialTriggerCutPosition(),2);
            PlotExampleInvMassBinsMC(fHistoMappingTrueFullMesonInvMassPtBins[fExampleBin], fHistoMappingTrueMesonCaloPhotonInvMassPtBins[fExampleBin],  NULL,
                                    fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[fExampleBin], fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[fExampleBin], fExampleBin,
                                    outputDir.Data(),Suffix.Data(), fMesonMassPlotRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                                    fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, triggerInt.Atoi(), fMode, addSig);


        } else if (fAdvancedMesonQA && (fMode == 2 || fMode == 13 || fMode == 3)){
            nameMesonTrue       = Form("%s_TrueMesonDecomposedPhotonsAndElectron%s", plotPrefix.Data(), plotSuffix.Data());
            cout << nameMesonTrue.Data() << endl;
            PlotTrueInvMassSplittedInPhotonAndElectronInPtBins(fHistoMappingTrueFullMesonInvMassPtBins, fHistoMappingTrueMesonCaloPhotonInvMassPtBins,  NULL,
                                                            fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins, NULL, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassPlotRange, fdate,
                                                            fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem,
                                                            fMode);
            TString triggerInt  = fEventCutSelectionRead(GetEventSelectSpecialTriggerCutPosition(),2);
            PlotExampleInvMassBinsMC(fHistoMappingTrueFullMesonInvMassPtBins[fExampleBin], fHistoMappingTrueMesonCaloPhotonInvMassPtBins[fExampleBin],  NULL,
                                    fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[fExampleBin], NULL, fExampleBin,
                                    outputDir.Data(),Suffix.Data(), fMesonMassPlotRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                                    fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, triggerInt.Atoi(), fMode, addSig);

        }

        if( fEnableNormBckHistoComparisonToTrueBck ){
            nameMeson       = Form("%s_MesonWithBckAndTrueBck%s", plotPrefix.Data(), plotSuffix.Data());
            nameCanvas      = "MesonWithBckAndTrueBckCanvas";
            namePad         = "MesonWithBckAndTrueBckPad";
            cout << nameMeson.Data() << endl;

            PlotInvMassInPtBins( fHistoMappingBackNormAndRemainingBGInvMassPtBin, fHistoMappingBackNormInvMassPtBin, fHistoMappingTrueAllBckInvMassPtBins, fHistoMappingTrueGGBckInvMassPtBins, fHistoMappingTrueContBckInvMassPtBins, fHistoMappingTrueMesonContainedInvMassPtBins, fHistoMappingTrueAsymEClusInvMassPtBins, nameMeson, nameCanvas, namePad,  fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem);
        }
    }

    CreatePtHistos();
    FillPtHistos();


    //Calculate the raw yield assuming different background fit functions and taking the ratio to the default pol1
    //normal
    fHistoYieldDiffBckRatios[0]->Add(fHistoYieldMeson[0]);
    fHistoYieldDiffBckRatios[0]->Divide(fHistoYieldDiffBckRatios[0],fHistoYieldDiffBck[0],1,1,"B");
    fHistoYieldDiffBckRatios[1]->Add(fHistoYieldMeson[0]);
    fHistoYieldDiffBckRatios[1]->Divide(fHistoYieldDiffBckRatios[1],fHistoYieldDiffBck[1],1,1,"B");
    //normal exp
    fHistoYieldDiffBckRatios[2]->Add(fHistoYieldMeson[0]);
    fHistoYieldDiffBckRatios[2]->Divide(fHistoYieldDiffBckRatios[2],fHistoYieldDiffBck[2],1,1,"B");
    fHistoYieldDiffBckRatios[3]->Add(fHistoYieldMeson[0]);
    fHistoYieldDiffBckRatios[3]->Divide(fHistoYieldDiffBckRatios[3],fHistoYieldDiffBck[3],1,1,"B");
    //normal exp2
    fHistoYieldDiffBckRatios[4]->Add(fHistoYieldMeson[0]);
    fHistoYieldDiffBckRatios[4]->Divide(fHistoYieldDiffBckRatios[4],fHistoYieldDiffBck[4],1,1,"B");
    fHistoYieldDiffBckRatios[5]->Add(fHistoYieldMeson[0]);
    fHistoYieldDiffBckRatios[5]->Divide(fHistoYieldDiffBckRatios[5],fHistoYieldDiffBck[5],1,1,"B");

    //normal
    for(Int_t m=1;m<fHistoYieldMeson[0]->GetNbinsX()+1;m++){
        fHistoYieldDiffBckResult[0]->SetBinContent(m,TMath::Sqrt(pow(1-fHistoYieldDiffBckRatios[0]->GetBinContent(m),2)+pow(1-fHistoYieldDiffBckRatios[1]->GetBinContent(m),2))/TMath::Sqrt(2));
        fHistoYieldDiffBckResult[0]->SetBinError(m,fHistoYieldDiffBckRatios[0]->GetBinError(m));
        if(fHistoYieldDiffBckResult[0]->GetBinContent(m)==1){ fHistoYieldDiffBckResult[0]->SetBinContent(m,0); fHistoYieldDiffBckResult[0]->SetBinError(m,0); }
    }
    //normal exp
    for(Int_t m=1;m<fHistoYieldMeson[0]->GetNbinsX()+1;m++){
        fHistoYieldDiffBckResult[1]->SetBinContent(m,TMath::Sqrt(pow(1-fHistoYieldDiffBckRatios[2]->GetBinContent(m),2)+pow(1-fHistoYieldDiffBckRatios[3]->GetBinContent(m),2))/TMath::Sqrt(2));
        fHistoYieldDiffBckResult[1]->SetBinError(m,fHistoYieldDiffBckRatios[2]->GetBinError(m));
        if(fHistoYieldDiffBckResult[1]->GetBinContent(m)==1){ fHistoYieldDiffBckResult[1]->SetBinContent(m,0); fHistoYieldDiffBckResult[1]->SetBinError(m,0); }
    }
    //normal exp2
    for(Int_t m=1;m<fHistoYieldMeson[0]->GetNbinsX()+1;m++){
        fHistoYieldDiffBckResult[2]->SetBinContent(m,TMath::Sqrt(pow(1-fHistoYieldDiffBckRatios[4]->GetBinContent(m),2)+pow(1-fHistoYieldDiffBckRatios[5]->GetBinContent(m),2))/TMath::Sqrt(2));
        fHistoYieldDiffBckResult[2]->SetBinError(m,fHistoYieldDiffBckRatios[4]->GetBinError(m));
        if(fHistoYieldDiffBckResult[2]->GetBinContent(m)==1){ fHistoYieldDiffBckResult[2]->SetBinContent(m,0); fHistoYieldDiffBckResult[2]->SetBinError(m,0); }
    }

    //normal
    TCanvas* canvasDiffBck = new TCanvas("canvasDiffBck","",1550,1200);
    canvasDiffBck->SetTickx();
    canvasDiffBck->SetTicky();

    DrawGammaSetMarker(fHistoYieldDiffBckResult[0], 20, 1.3, kBlack, kBlack);
    DrawAutoGammaMesonHistos( fHistoYieldDiffBckResult[0],
                                "", "p_{T} (GeV/c)", "dev. raw yield",
                                kFALSE, 3.,0.,  kFALSE,
                                kTRUE, -0.2,0.5,
                                kFALSE, 0., 16.);

    DrawGammaSetMarker(fHistoYieldDiffBckResult[1], 33, 1.3, kRed+2, kRed+2);
    fHistoYieldDiffBckResult[1]->DrawCopy("e1,p,SAME");

    DrawGammaSetMarker(fHistoYieldDiffBckResult[2], 34, 1.3, kGreen+2, kGreen+2);
    fHistoYieldDiffBckResult[2]->DrawCopy("e1,p,SAME");

    canvasDiffBck->Update();

    TLegend* legendDiffBck = new TLegend(0.15,0.8,0.4,0.95);
    legendDiffBck->SetFillColor(0);
    legendDiffBck->SetLineColor(0);
    legendDiffBck->SetTextSize(0.04);
    legendDiffBck->AddEntry(fHistoYieldDiffBckResult[0],Form("dev. raw yield wrt pol2 Bck fit for %s",fPrefix.Data()),"p");
    legendDiffBck->AddEntry(fHistoYieldDiffBckResult[1],Form("dev. raw yield wrt exp Bck fit for %s",fPrefix.Data()),"p");
    legendDiffBck->AddEntry(fHistoYieldDiffBckResult[2],Form("dev. raw yield wrt exp2 Bck fit for %s",fPrefix.Data()),"p");
    legendDiffBck->Draw();

    if (fIsMC) canvasDiffBck->SaveAs(Form("%s/%s_MC_DiffBck_%s.%s",outputDir.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));
    else canvasDiffBck->SaveAs(Form("%s/%s_data_DiffBck_%s.%s",outputDir.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));



    if (!fIsMC && meson.Contains("Pi0") ){
        fHaveCocktailInputForSec     = LoadSecondaryPionsFromCocktailFile(cutSelection,optionEnergy);
        if(!fHaveCocktailInputForSec)
            fHaveToyMCInputForSec     = LoadSecondaryPionsFromExternalFile();
        if (fHaveCocktailInputForSec){
            cout << "SECONDARIES: I am gonna add the cocktail output to the uncorrected file" << endl;
        } else if (fHaveToyMCInputForSec){
            cout << "SECONDARIES: I am gonna add the toy MC output to the uncorrected file" << endl;
        } else {
            cout << "SECONDARIES: no ToyMC or cocktail input has been found for the secondaries" << endl;
        }
    }

    ///*********************** Lambda tail
    TCanvas* canvasLambdaTail = new TCanvas("canvasLambdaTail","",1550,1200);  // gives the page size
    canvasLambdaTail->SetTickx();
    canvasLambdaTail->SetTicky();

    DrawGammaSetMarker(fHistoLambdaTail, 20, 1., kBlack, kBlack);
    Float_t maxPlotLambda               = fMesonLambdaTailRangeNominal[1]*1.2;
    if (fMesonLambdaTailRangeNominal[1] == fMesonLambdaTailRangeNominal[0])
        maxPlotLambda                   = fMesonLambdaTailRangeNominal[1]*2;
    if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0){
        DrawAutoGammaMesonHistos( fHistoLambdaTail,
                                "", "p_{T} (GeV/c)", "#lambda",
                                kFALSE, 3.,0.,  kFALSE,
                                kTRUE, 0.,maxPlotLambda,
                                kFALSE, 0., 10.);
    } else {
        DrawAutoGammaMesonHistos( fHistoLambdaTail,
                                "", "p_{T} (GeV/c)", "#lambda",
                                kFALSE, 3.,0.,  kFALSE,
                                kTRUE, 5e-3,maxPlotLambda,
                                kFALSE, 0., 10.);
    }
    if (fIsMC){
        DrawGammaSetMarker(fHistoTrueLambdaTail, 24, 1., kRed+2, kRed+2);
        fHistoTrueLambdaTail->Draw("same,pe");
    }
    canvasLambdaTail->Update();

    TLegend* legendLambdaTail = GetAndSetLegend2(0.15,0.90,0.4,0.94, 0.04*1200,1);
    legendLambdaTail->AddEntry(fHistoLambdaTail,Form("Lambda tail parameter for %s",fPrefix.Data()),"p");
    legendLambdaTail->Draw();

    DrawGammaLines(0., fBinsPt[fNBinsPt], fMesonLambdaTailRangeNominal[0], fMesonLambdaTailRangeNominal[0], 1, kRed+1, 2);
    DrawGammaLines(0., fBinsPt[fNBinsPt], fMesonLambdaTail, fMesonLambdaTail, 1, kGray+2, 2);
    DrawGammaLines(0., fBinsPt[fNBinsPt], fMesonLambdaTailRangeNominal[1], fMesonLambdaTailRangeNominal[1], 1, kRed+1, 2);

    if (fIsMC) canvasLambdaTail->SaveAs(Form("%s/%s_MC_LambdaTail_%s.%s",outputDirMon.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));
    else canvasLambdaTail->SaveAs(Form("%s/%s_data_LambdaTail_%s.%s",outputDirMon.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));


    ///*********************** Mass
    TCanvas* canvasMesonMass = new TCanvas("canvasMesonMass","",1550,1200);  // gives the page size
    canvasMesonMass->SetTickx();
    canvasMesonMass->SetTicky();

    Double_t maxMesonMassRange = 0.140;
    Double_t minMesonMassRange = 0.132;
    if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0){
        if (fEnergyFlag.Contains("PbPb")){
            maxMesonMassRange = 0.160;
            minMesonMassRange = 0.130;
        } else if (fEnergyFlag.CompareTo("13TeVLowB") == 0 && fMode == 0 ) {
            maxMesonMassRange = 0.150;
            minMesonMassRange = 0.132;
        } else {
            maxMesonMassRange = 0.140;
            minMesonMassRange = 0.132;
        }
    } else if (fPrefix.CompareTo("Eta") ==0) {
        maxMesonMassRange = 0.64;
        minMesonMassRange = 0.46;
    } else if (fPrefix.CompareTo("EtaPrime") ==0) {
        maxMesonMassRange = 0.8;
        minMesonMassRange = 1.1;
    }


    DrawGammaSetMarker(fHistoMassMeson, 20, 1., kBlack, kBlack);
    DrawAutoGammaMesonHistos(   fHistoMassMeson,
                                "", "p_{T} (GeV/c)", Form("Mass (GeV/c^{2})"),
                                kFALSE, 3.,0., kFALSE,
                                kTRUE, minMesonMassRange, maxMesonMassRange,
                                kFALSE, 0., 10.);

    if (fIsMC > 0){
        DrawGammaSetMarker(fHistoTrueMassMeson, 24, 1., kRed+2, kRed+2);
        fHistoTrueMassMeson->Draw("same,pe");
    }
    canvasMesonMass->Update();

    TLegend* legendMesonMass = GetAndSetLegend2(0.15,0.90,0.4,0.94, 0.04*1200,1);
    legendMesonMass->AddEntry(fHistoMassMeson,Form("%s mass",fPrefix.Data()),"p");
    legendMesonMass->Draw();

    DrawGammaLines(0., fBinsPt[fNBinsPt], fMesonMassRange[0], fMesonMassRange[0], 1, kRed+1, 2);
    DrawGammaLines(0., fBinsPt[fNBinsPt], fMesonMassExpect, fMesonMassExpect, 1, kGray+2, 2);
    DrawGammaLines(0., fBinsPt[fNBinsPt], fMesonMassRange[1], fMesonMassRange[1], 1, kRed+1, 2);

    if (fIsMC) canvasMesonMass->SaveAs(Form("%s/%s_MC_MesonMass_%s.%s",outputDirMon.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));
    else canvasMesonMass->SaveAs(Form("%s/%s_data_MesonMass_%s.%s",outputDirMon.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));


    ///*********************** Width
    TCanvas* canvasMesonFWHM = new TCanvas("canvasMesonFWHM","",1550,1200);  // gives the page size
    canvasMesonFWHM->SetTickx();
    canvasMesonFWHM->SetTicky();

    DrawGammaSetMarker(fHistoFWHMMeson, 20, 1., kBlack, kBlack);
    if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0){
        DrawAutoGammaMesonHistos( fHistoFWHMMeson,
                                    "", "p_{T} (GeV/c)","FWHM (GeV/c^{2})",
                                    kFALSE, 3.,0., kFALSE,
                                    kTRUE, -0.004, fMesonWidthRange[1]*4, // *4 because fMesonWidthRange[1] corresponds to fit parameter sigma, but histogram is FWHM
                                    kFALSE, 0., 10.);
    } else {
        DrawAutoGammaMesonHistos( fHistoFWHMMeson,
                                    "", "p_{T} (GeV/c)","FWHM (GeV/c^{2})",
                                    kFALSE, 3.,0., kFALSE,
                                    kTRUE, 0., fMesonWidthRange[1]*4,
                                    kFALSE, 0., 10.);
    }
    canvasMesonFWHM->Update();

    TLegend* legendMesonFWHM = GetAndSetLegend2(0.2,0.12,0.45,0.16, 0.04*1200,1);
    legendMesonFWHM->AddEntry(fHistoFWHMMeson,Form("%s FWHM",fPrefix.Data()),"p");
    legendMesonFWHM->Draw();

    if (fIsMC) canvasMesonFWHM->SaveAs(Form("%s/%s_MC_MesonFWHM_%s.%s",outputDirMon.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));
    else canvasMesonFWHM->SaveAs(Form("%s/%s_data_MesonFWHM_%s.%s",outputDirMon.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));

    ///*********************** Amplitude fit parameter monitoring
    TCanvas* canvasAmplitude = new TCanvas("canvasAmplitude","",1550,1200);  // gives the page size
    canvasAmplitude->SetTickx();
    canvasAmplitude->SetTicky();
    canvasAmplitude->SetLogy();
    DrawGammaSetMarker(fHistoAmplitude, 20, 1., kBlack, kBlack);
    DrawAutoGammaMesonHistos( fHistoAmplitude,
                            "", "p_{T} (GeV/c)", "A",
                            kFALSE, 3.,0.,  kTRUE,
                            kFALSE, 0., 0.,
                            kFALSE, 0., 10.);

    canvasAmplitude->Update();

    TLegend* legendAmplitude = new TLegend(0.45,0.8,0.7,0.95);
    legendAmplitude->SetFillColor(0);
    legendAmplitude->SetLineColor(0);
    legendAmplitude->SetTextSize(0.04);
    legendAmplitude->AddEntry(fHistoAmplitude,Form("Amplitude parameter for %s",fPrefix.Data()),"p");
    legendAmplitude->Draw();

    if (fIsMC) canvasAmplitude->SaveAs(Form("%s/%s_MC_Amplitude_%s.%s",outputDirMon.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));
    else canvasAmplitude->SaveAs(Form("%s/%s_data_Amplitude_%s.%s",outputDirMon.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));

    ///*********************** sigma fit parameter monitoring
    TCanvas* canvasSigma = new TCanvas("canvasSigma","",1550,1200);  // gives the page size
    canvasSigma->SetTickx();
    canvasSigma->SetTicky();

    DrawGammaSetMarker(fHistoSigma, 20, 1., kBlack, kBlack);
if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0){
        DrawAutoGammaMesonHistos(   fHistoSigma,
                                    "", "p_{T} (GeV/c)","#sigma (GeV/c^{2})",
                                    kFALSE, 3.,0., kFALSE,
                                    kTRUE, -0.004, fMesonWidthRange[1]*1.7,
                                    kFALSE, 0., 10.);
    } else {
        DrawAutoGammaMesonHistos(   fHistoSigma,
                                    "", "p_{T} (GeV/c)","sigma (GeV/c^{2})",
                                    kFALSE, 3.,0., kFALSE,
                                    kTRUE, 0., fMesonWidthRange[1]*1.7,
                                    kFALSE, 0., 10.);
    }
    if (fIsMC > 0){
        DrawGammaSetMarker(fHistoTrueSigma, 24, 1., kRed+2, kRed+2);
        fHistoTrueSigma->Draw("same,pe");
    }

    canvasSigma->Update();

    TLegend* legendSigma = new TLegend(0.15,0.8,0.4,0.95);
    legendSigma->SetFillColor(0);
    legendSigma->SetLineColor(0);
    legendSigma->SetTextSize(0.04);
    legendSigma->AddEntry(fHistoSigma,Form("Sigma parameter for %s",fPrefix.Data()),"p");
    legendSigma->Draw();

    DrawGammaLines(0., fBinsPt[fNBinsPt], fMesonWidthRange[0], fMesonWidthRange[0], 1, kRed+1, 2);
    DrawGammaLines(0., fBinsPt[fNBinsPt], fMesonWidthExpect, fMesonWidthExpect, 1, kGray+2, 2);
    DrawGammaLines(0., fBinsPt[fNBinsPt], fMesonWidthRange[1], fMesonWidthRange[1], 1, kRed+1, 2);

    if (fIsMC) canvasSigma->SaveAs(Form("%s/%s_MC_Sigma_%s.%s",outputDirMon.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));
    else canvasSigma->SaveAs(Form("%s/%s_data_Sigma_%s.%s",outputDirMon.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));


    ///*********************** Mid Pt fit (for monitoring)
    TCanvas* canvasMesonMidPt = new TCanvas("canvasMesonMidPt","",1550,1200);  // gives the page size
    canvasMesonMidPt->SetTickx();
    canvasMesonMidPt->SetTicky();

    SetStyleHisto(fFittingHistMidPtSignalSub,1,1,kBlack);
    fFitSignalInvMassMidPt->SetLineColor(kRed+1);
    fFittingHistMidPtSignalSub->Draw("same");
    fFitSignalInvMassMidPt->Draw("same");
    canvasMesonMidPt->Update();

    TLegend* legendMesonMidPt = new TLegend(0.15,0.62,0.45,0.76);
    legendMesonMidPt->SetFillColor(0);
    legendMesonMidPt->SetLineColor(0);
    legendMesonMidPt->SetTextSize(0.04);
    legendMesonMidPt->AddEntry(fFittingHistMidPtSignalSub,Form("%s Raw Mid-Pt",fPrefix.Data()),"l");
    legendMesonMidPt->AddEntry(fFitSignalInvMassMidPt,Form("%s Fit Mid-Pt",fPrefix.Data()),"l");
    legendMesonMidPt->Draw();

    if (fIsMC) canvasMesonMidPt->SaveAs(Form("%s/%s_MC_MesonSubtractedFittingMidPt_%s.%s",outputDirMon.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));
    else canvasMesonMidPt->SaveAs(Form("%s/%s_data_MesonSubtractedFittingMidPt_%s.%s",outputDirMon.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));


    ///*********************** SPD pileup monitoring
    if(fHistoPileUpVertexDistance && fHistoPileUpVertexDistance_SPDPileup && fHistoPileUpVertexDistance_TrackletHits){
        TCanvas* canvasSPDPileUp = new TCanvas("canvasSPDPileUp","",1550,1200);  // gives the page size
        canvasSPDPileUp->SetTickx();
        canvasSPDPileUp->SetTicky();
        canvasSPDPileUp->SetLogy();

        DrawGammaSetMarker(fHistoPileUpVertexDistance, 20, 1., kBlack, kBlack);
        DrawGammaSetMarker(fHistoPileUpVertexDistance_SPDPileup, 20, 1, kBlue, kBlue);
        DrawGammaSetMarker(fHistoPileUpVertexDistance_TrackletHits, 20, 1, kRed, kRed);
        DrawAutoGammaMesonHistos(   fHistoPileUpVertexDistance,
                                    "", "distance between vertices", "numer of entries / bin",
                                    kFALSE, 3.,0., kTRUE,
                                    kFALSE, -0.004, 1.,
                                    kFALSE, 0., 10.);
        fHistoPileUpVertexDistance_SPDPileup->Draw("same");
        fHistoPileUpVertexDistance_TrackletHits->Draw("same");
        canvasSPDPileUp->Update();

        TLegend* legendSPDPileUp = new TLegend(0.15,0.8,0.4,0.95);
        legendSPDPileUp->SetFillColor(0);
        legendSPDPileUp->SetLineColor(0);
        legendSPDPileUp->SetTextSize(0.04);
        legendSPDPileUp->AddEntry(fHistoPileUpVertexDistance,"PileUpVertexDistance","p");
        legendSPDPileUp->AddEntry(fHistoPileUpVertexDistance_SPDPileup,"removed by SPDPileUp");
        legendSPDPileUp->AddEntry(fHistoPileUpVertexDistance_TrackletHits,"removed by SPDTrackletHits");
        legendSPDPileUp->Draw();

        if (fIsMC) canvasSPDPileUp->SaveAs(Form("%s/MC_SPDPileUpTotal_%s.%s",outputDirMon.Data(),fCutSelection.Data(),Suffix.Data()));
        else canvasSPDPileUp->SaveAs(Form("%s/data_SPDPileUpTOtal_%s.%s",outputDirMon.Data(),fCutSelection.Data(),Suffix.Data()));

        TF1* fFitGausPileUp = new TF1("gausFitPileUp",fitGaussianPileUp,-15,15,3);
        fFitGausPileUp->SetParameters(1,0,30);
        fHistoPileUpVertexDistance_SPDPileup->Fit("gausFitPileUp","QRMNE0");

        TF1* fFitGausTrackletHits = new TF1("gausFitPileUp2",fitGaussianPileUp2,-0.8,0.8,3);
        fFitGausTrackletHits->SetParameters(1,0,1000);
        fHistoPileUpVertexDistance_TrackletHits->Fit("gausFitPileUp2","QRMNE0");

        DrawGammaSetMarker(fHistoPileUpVertexDistance_SPDPileup, 20, 1, kBlue, kBlue);
        DrawGammaSetMarker(fHistoPileUpVertexDistance_TrackletHits, 20, 1, kGreen+3, kGreen+3);
        DrawAutoGammaMesonHistos(   fHistoPileUpVertexDistance_TrackletHits,
                                    "", "distance between vertices", "numer of entries / bin",
                                    kFALSE, 3.,0., kTRUE,
                                    kFALSE, -0.004, 1.,
                                    kFALSE, 0., 10.);
        fHistoPileUpVertexDistance_SPDPileup->SetLineWidth(2);
        fHistoPileUpVertexDistance_SPDPileup->Draw("same");
        // fFitGausPileUp->SetLineColor(kGreen+2);
        // fFitGausPileUp->Draw("same");
        canvasSPDPileUp->Update();

        TF1* fFitGausPileUpFull = new TF1("gausFitPileUpFull","gaus",-15,15);
        fFitGausPileUpFull->SetParameters(fFitGausPileUp->GetParameters());
        fFitGausPileUpFull->SetLineColor(kBlack);
        fFitGausPileUpFull->SetLineStyle(2);
        fFitGausPileUpFull->SetLineWidth(5);
        fFitGausPileUpFull->Draw("same");

        TF1* fFitGausTrackletHitsFull = new TF1("gausFitPileUpFull2","gaus",-0.8,0.8);
        fFitGausTrackletHitsFull->SetParameters(fFitGausTrackletHits->GetParameters());
        fFitGausTrackletHitsFull->SetLineColor(kRed);
        fFitGausTrackletHitsFull->SetLineStyle(2);
        fFitGausTrackletHitsFull->SetLineWidth(5);
        fFitGausTrackletHitsFull->Draw("same");

        Double_t pileUpIntegralTotal = fFitGausPileUpFull->Integral(-15,15);
        Double_t pileUpIntegral = fFitGausPileUpFull->Integral(-15,-0.9);
        pileUpIntegral += fFitGausPileUpFull->Integral(0.9,15);
        pileUpIntegral += fFitGausTrackletHitsFull->Integral(-0.8,0.8);

        Double_t ratioPileUpIntegrals = pileUpIntegral / pileUpIntegralTotal;

        // Double_t integralPileUpVertexDistance_TrackletHits = fHistoPileUpVertexDistance_TrackletHits->Integral(1,fHistoPileUpVertexDistance_TrackletHits->GetNbinsX());
        // Double_t ratioIntegralPileUpIntegrals = (integralPileUpVertexDistance_TrackletHits+pileUpIntegral) / pileUpIntegralTotal;

        TLegend* legendSPDPileUp2 = new TLegend(0.15,0.7,0.5,0.95);
        legendSPDPileUp2->SetFillColor(0);
        legendSPDPileUp2->SetLineColor(0);
        legendSPDPileUp2->SetTextSize(0.04);
        legendSPDPileUp2->AddEntry(fHistoPileUpVertexDistance_SPDPileup,"removed by SPDPileUp");
        legendSPDPileUp2->AddEntry(fHistoPileUpVertexDistance_TrackletHits,"removed by SPDTrackletHits");
        legendSPDPileUp2->AddEntry(fFitGausPileUpFull,"fit of SPDPileUp");
        legendSPDPileUp2->AddEntry(fFitGausTrackletHitsFull,"fit of SPDTrackletHits");
        legendSPDPileUp2->AddEntry((TObject*)0,"","");
        legendSPDPileUp2->AddEntry((TObject*)0,"SPD pile-up cuts","");
        legendSPDPileUp2->AddEntry((TObject*)0,Form("efficiency: %.2f",ratioPileUpIntegrals),"");
        // legendSPDPileUp2->AddEntry((TObject*)0,Form("%.2f",ratioIntegralPileUpIntegrals),"");
        legendSPDPileUp2->Draw();

        if (fIsMC) canvasSPDPileUp->SaveAs(Form("%s/MC_SPDPileUp_%s.%s",outputDirMon.Data(),fCutSelection.Data(),Suffix.Data()));
        else canvasSPDPileUp->SaveAs(Form("%s/data_SPDPileUp_%s.%s",outputDirMon.Data(),fCutSelection.Data(),Suffix.Data()));
    }
    // **************************************************************************************************************
    // ************************ Chi2/ndf compared MC vs Data ********************************************************
    // **************************************************************************************************************
    TCanvas* canvasChi2 = new TCanvas("canvasChi2","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasChi2, 0.092, 0.01, 0.02, 0.082);

    Double_t maxChi2    = fHistoChi2[0]->GetMaximum();
    for (Int_t m = 1; m < 4; m++){
        if (maxChi2 < fHistoChi2[m]->GetMaximum())
            maxChi2     = fHistoChi2[m]->GetMaximum();
    }
    maxChi2             = maxChi2*1.2;

    TLegend* legendChi2 = GetAndSetLegend2(0.75, 0.95-(0.035*4), 0.95, 0.95, 0.035, 1, "", 42, 0.25);

    fHistoChi2[0]->GetYaxis()->SetRangeUser(0, maxChi2);
    DrawAutoGammaMesonHistos(   fHistoChi2[0],
                                "", "#it{p}_{T} (GeV/#it{c})", "#it{#chi}^{2}/ndf",
                                kFALSE, 0., 0.7, kFALSE,
                                kFALSE, 0., 0.7,
                                kFALSE, 0., 10.);
    DrawGammaSetMarker(fHistoChi2[0], 20, 2, kCyan+2, kCyan+2);
    fHistoChi2[0]->DrawCopy("same,e1,p");
    legendChi2->AddEntry(fHistoChi2[0],"pol1 BG","p");

    Color_t colorFit[3] = {kRed+1, kAzure+2, 807};
    Style_t styleFit[3] = {34, 21, 33};

    for (Int_t m = 1; m < 4; m++){
        DrawGammaSetMarker(fHistoChi2[m], styleFit[m-1], 2, colorFit[m-1], colorFit[m-1]);
        fHistoChi2[m]->DrawCopy("same,e1,p");
        legendChi2->AddEntry(fHistoChi2[m],labelsOtherFits[m-1],"p");
    }
    fHistoChi2[0]->DrawCopy("same,e1,p");
    legendChi2->Draw();


    PutProcessLabelAndEnergyOnPlot(0.15, 0.25, 0.035, fCollisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
    canvasChi2->Update();
    if (fIsMC) canvasChi2->SaveAs(Form("%s/%s_MC_Chi2FitComp_%s.%s",outputDir.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));
    else canvasChi2->SaveAs(Form("%s/%s_data_Chi2FitComp_%s.%s",outputDir.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));

    // **************************************************************************************************************
    // ************************ ResBG compared MC vs Data ********************************************************
    // **************************************************************************************************************
    TCanvas* canvasResBG = new TCanvas("canvasResBG","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasResBG, 0.092, 0.01, 0.02, 0.082);

    Double_t maxResBG    = fHistoResBGYield[0]->GetMaximum();
    for (Int_t m = 1; m < 4; m++){
        if (maxResBG < fHistoResBGYield[m]->GetMaximum())
            maxResBG     = fHistoResBGYield[m]->GetMaximum();
    }
    maxResBG             = maxResBG*1.2;

    Double_t minResBG    = fHistoResBGYield[0]->GetMinimum();
    for (Int_t m = 1; m < 4; m++){
        if (minResBG > fHistoResBGYield[m]->GetMinimum())
            minResBG     = fHistoResBGYield[m]->GetMinimum();
    }
    if (minResBG < 0)
        minResBG             = minResBG*1.4;
    else
        minResBG             = minResBG*0.8;

    TLegend* legendResBG = GetAndSetLegend2(0.75, 0.80-(0.035*4), 0.95, 0.80, 0.035, 1, "", 42, 0.25);

    fHistoResBGYield[0]->GetYaxis()->SetRangeUser(minResBG, maxResBG);
    DrawAutoGammaMesonHistos( fHistoResBGYield[0],
                                "", "#it{p}_{T} (GeV/#it{c})", "Res BG yield",
                                kFALSE, 0., 0.7, kFALSE,
                                kFALSE, 0., 0.7,
                                kFALSE, 0., 10.);
    DrawGammaSetMarker(fHistoResBGYield[0], 20, 2, kCyan+2, kCyan+2);
    fHistoResBGYield[0]->DrawCopy("same,e1,p");
    legendResBG->AddEntry(fHistoResBGYield[0],"pol1 BG","p");

    for (Int_t m = 1; m < 4; m++){
        DrawGammaSetMarker(fHistoResBGYield[m], styleFit[m-1], 2, colorFit[m-1], colorFit[m-1]);
        fHistoResBGYield[m]->DrawCopy("same,e1,p");
        legendResBG->AddEntry(fHistoResBGYield[m],labelsOtherFits[m-1],"p");
    }
    fHistoResBGYield[0]->DrawCopy("same,e1,p");
    legendResBG->Draw();

    PutProcessLabelAndEnergyOnPlot(0.70, 0.95, 0.035, fCollisionSystem.Data(), fTextMeasurement.Data(), fDetectionProcess.Data());
    canvasResBG->Update();
    if (fIsMC) canvasResBG->SaveAs(Form("%s/%s_MC_ResBGYieldFitComp_%s.%s",outputDir.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));
    else canvasResBG->SaveAs(Form("%s/%s_data_ResBGYieldFitComp_%s.%s",outputDir.Data(),fPrefix.Data(),fCutSelection.Data(),Suffix.Data()));

    // **************************************************************************************************************
    // Filling MC hists
    // **************************************************************************************************************

    if(fIsMC){
        // Rebin MC histograms for acceptance and input with possible weights
        FillHistosArrayMC(fHistoMCMesonPtWithinAcceptance, fHistoMCMesonPt, fDeltaPt);
        // Rebin MC histograms for acceptance and input without weights
        if (fHistoMCMesonPtWithinAcceptanceWOWeights) FillHistosArrayMCWOWeights(fHistoMCMesonPtWithinAcceptanceWOWeights, fHistoMCMesonPtWOWeights, fDeltaPt);
        if (fHistoMCMesonPtWithinAcceptanceWOEvtWeights) FillHistosArrayMCWOEvtWeights(fHistoMCMesonPtWithinAcceptanceWOEvtWeights, fHistoMCMesonPtWOEvtWeights, fDeltaPt);

        // Calculation of meson acceptance with possible weighted input
        CalculateMesonAcceptance();
        // Calculation of meson acceptance without weights for input
        if (fHistoMCMesonPtWithinAcceptanceWOWeights) CalculateMesonAcceptanceWOWeights();
        if (fHistoMCMesonPtWithinAcceptanceWOEvtWeights) CalculateMesonAcceptanceWOEvtWeights();

        if (fNewMCOutput){
            FillMCSecondaryHistAndCalculateAcceptance(fHistoMCSecPi0SourcePt, fHistoMCSecPi0WAccSourcePt);
        }

        // calculate pure rec efficiency as in data with fully unweighted histograms if possible
        // ATTENTION: if unweighted histograms are not available this efficiency should not be used for anything!!!!

        for (Int_t k = 0; k < 3; k++){
            fNameHistoEffi                      = Form("Meson%sEffiPt",nameIntRange[k].Data());
            cout << fNameHistoEffi.Data() << endl;
            if (fHistoMCMesonPtWithinAcceptanceWOWeights)
                fHistoMonteMesonEffiPt[k]       = CalculateMesonEfficiency(fHistoYieldMeson[k], fHistoYieldTrueSecMeson[k], fHistoMCMesonWithinAccepPtWOWeights, fNameHistoEffi);
            else
                fHistoMonteMesonEffiPt[k]       = CalculateMesonEfficiency(fHistoYieldMeson[k], fHistoYieldTrueSecMeson[k], fHistoMCMesonWithinAccepPt, fNameHistoEffi);

            fNameHistoEffi                      = Form("Meson%sEffiPt",nameIntRange[k+3].Data());
            cout << fNameHistoEffi.Data() << endl;
            if (fHistoMCMesonPtWithinAcceptanceWOWeights)
                fHistoMonteMesonEffiPt[k+3]       = CalculateMesonEfficiency(fHistoYieldMeson[k+3], fHistoYieldTrueSecMeson[k], fHistoMCMesonWithinAccepPtWOWeights, fNameHistoEffi);
            else
                fHistoMonteMesonEffiPt[k+3]       = CalculateMesonEfficiency(fHistoYieldMeson[k+3], fHistoYieldTrueSecMeson[k], fHistoMCMesonWithinAccepPt, fNameHistoEffi);
        }

        for (Int_t k = 0; k < 3; k++){
            // True meson efficiencies for fully unweighted MC, to be compared to MesonEffiPt if those have been created with unweighted histograms
            fNameHistoEffi                          = Form("TrueMeson%sEffiPtUnweighted",nameIntRange[k].Data());
            cout << fNameHistoEffi.Data() << endl;
            if (fHistoMCMesonPtWithinAcceptanceWOWeights)
                fHistoMCTrueMesonEffiPtUnweighted[k]= CalculateMesonEfficiency(fHistoYieldTrueMesonUnweighted[k], NULL, fHistoMCMesonWithinAccepPtWOWeights, fNameHistoEffi);

            // True meson efficiencies with possibly fully weighted inputs on a meson by meson basis in the aliphysics task, should always be used if you start weighting the MC
            // True Meson (only once case, because no normalization)
            fNameHistoEffi                          = Form("TrueMeson%sEffiPt",nameIntRange[k].Data());
            cout << fNameHistoEffi.Data() << endl;
            fHistoMCTrueMesonEffiPt[k]              = CalculateMesonEfficiency(fHistoYieldTrueMeson[k], NULL, fHistoMCMesonWithinAccepPt, fNameHistoEffi);


            // True meson efficiencies with possibly fully weighted inputs taking the average weight per inv mass bin in the original binning of the TrueMesonInvMass vs pT plot
            // should give on average the same as TrueMesonEffiPt
            fNameHistoEffi                          = Form("TrueMeson%sEffiPtReweighted",nameIntRange[k].Data());
            cout << fNameHistoEffi.Data() << endl;
            fHistoMCTrueMesonEffiPtReweighted[k]    = CalculateMesonEfficiency(fHistoYieldTrueMesonReweighted[k], NULL, fHistoMCMesonWithinAccepPt, fNameHistoEffi);

            if (fNewMCOutput){
                for (Int_t j = 0; j < 4; j++){
                    fNameHistoEffi                          = Form("TrueSecFrom%s%sEffiPt",nameSecondaries[j].Data(), nameIntRange[k].Data());
                    cout << "trying to create: "<< fNameHistoEffi.Data() << endl;
                    if (fHistoYieldTrueSecMeson[k][j] && fHistoMCSecPi0PtWAccReb[j]){
                        fHistoMCTrueSecMesonEffiPt[k][j]        = CalculateMesonEfficiency(fHistoYieldTrueSecMeson[k][j], NULL, fHistoMCSecPi0PtWAccReb[j], fNameHistoEffi);
                    } else {
                        cout << Form("TrueSecFrom%s%sEffiPtReweighted",nameSecondaries[j].Data(), nameIntRange[k].Data()) << " could not be created " << endl;
                        if (!fHistoYieldTrueSecMeson[k][j])
                            cout << "true rec yield for " << nameSecondaries[j].Data() << " in int range: " << nameIntRange[k].Data() << " was missing" << endl;
                        if (!fHistoMCSecPi0PtWAccReb[j])
                            cout << "MC yield for " << nameSecondaries[j].Data() << " in acceptance was missing" << endl;
                    }
                }
            }
        }

        // Calculation of secondary fractions using unweighted histograms, as secondaries are never weighted
        TH1D* fHistoYieldTrueMesonSecPlusPrim[3];
        for (Int_t k = 0; k < 3; k++){
            fHistoYieldTrueMesonSecPlusPrim[k]          = (TH1D*)fHistoYieldTrueMesonUnweighted[k]->Clone(Form("fHistoYieldTrueMesonSecPlusPrim%s",nameIntRange[k].Data()));
            for (Int_t j = 0; j< 4; j++){
                fHistoYieldTrueMesonSecPlusPrim[k]->Add(fHistoYieldTrueSecMeson[k][j]);
            }
            for (Int_t j = 0; j< 4; j++){
                fHistoYieldTrueSecFracMeson[k][j] = CalculateSecondaryFractions(    fHistoYieldTrueMesonSecPlusPrim[k], fHistoYieldTrueSecMeson[k][j],
                                                                                    Form("TrueSecFrom%sFrac%s",nameSecondaries[j].Data(),nameIntRange[k].Data()));
            }
        }

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
void ProduceBckProperWeighting(TList* backgroundContainer,TList* motherContainer, Bool_t UseTHnSparse){

    if(UseTHnSparse){
        cout << "Using THnSparse for the background" << endl;
        THnSparseF* fSparseMotherZM;
        THnSparseF* fSparseBckZM;
        THnSparseF* fSparseMotherZPsi;
        THnSparseF* fSparseBckZPsi;

        fSparseMotherZM = (THnSparseF*)motherContainer->FindObject("Back_Mother_InvMass_Pt_z_m");
        fSparseBckZM = (THnSparseF*)backgroundContainer->FindObject("Back_Back_InvMass_Pt_z_m");
        if(fSparseMotherZM && fSparseBckZM){
            fUseRPBackground = kFALSE;
            cout << "with ZM bins estimation" << endl;
        }
        fSparseMotherZPsi = (THnSparseF*)motherContainer->FindObject("Back_Mother_InvMass_Pt_z_psi");
        fSparseBckZPsi = (THnSparseF*)backgroundContainer->FindObject("Back_Back_InvMass_Pt_z_psi");
        if(fSparseMotherZPsi && fSparseBckZPsi){
            fUseRPBackground = kTRUE;
            cout << "with ZPsi bins estimation" << endl;
        }

        for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
            //cout << "Pt:"<< iPt << endl;
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
                        //cout << m << "\t" << z << endl;
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
                        //cout << "Z:"<<  z << "\t psi: " << psi << endl;
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

        fHistoMotherZM = (TH2D*)motherContainer->FindObject("ESD_Mother_InvMass_Pt");
        fHistoMotherZM->Sumw2();
        fHistoBckZM = (TH2D*)backgroundContainer->FindObject("ESD_Background_InvMass_Pt");
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
    }
}

//****************************************************************************
//****** Initialization of arrays and variables, histograms for analysis *****
//****** depending on mesonType, number of bins, mode, energy, centrality ****
//****************************************************************************
void Initialize(TString setPi0, Int_t numberOfBins, Int_t triggerSet){

    cout << "meson in intialize function" <<  setPi0.Data() << endl;
    InitializeBinning(setPi0, numberOfBins, fEnergyFlag, fdirectphoton, fModeHeavy, fEventCutSelection, fClusterCutSelection, triggerSet, kFALSE, "", "", fGammaCutSelection);

    TString trigger         = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
    InitializeWindows(setPi0, fMode, trigger, triggerSet);

    // initialize integration array for integration windows: normal, wide, narrow, left, left wide, left narrow
    for (Int_t k = 0; k < 6; k++){
        fMesonCurIntRange[k]                                        = new Double_t[2];
    }

    // initialize integration pt-arrays for integration windows: normal, wide, narrow, left, left wide, left narrow
    for (Int_t k = 0; k < 6; k++){
        // Initialize yield arrays
        fGGYields[k]                                                = new Double_t[fNBinsPt];
        fBckYields[k]                                               = new Double_t[fNBinsPt];
        fTotalBckYields[k]                                          = new Double_t[fNBinsPt];
        fMesonYields[k]                                             = new Double_t[fNBinsPt];
        fMesonYieldsFunc[k]                                         = new Double_t[fNBinsPt];
        fMesonYieldsResidualBckFunc[k]                              = new Double_t[fNBinsPt];
        fMesonYieldsCorResidualBckFunc[k]                           = new Double_t[fNBinsPt];
        fMesonYieldsPerEvent[k]                                     = new Double_t[fNBinsPt];

        // Initialize error arrays
        fGGYieldsError[k]                                           = new Double_t[fNBinsPt];
        fBckYieldsError[k]                                          = new Double_t[fNBinsPt];
        fTotalBckYieldsError[k]                                     = new Double_t[fNBinsPt];
        fMesonYieldsError[k]                                        = new Double_t[fNBinsPt];
        fMesonYieldsFuncError[k]                                    = new Double_t[fNBinsPt];
        fMesonYieldsResidualBckFuncError[k]                         = new Double_t[fNBinsPt];
        fMesonYieldsCorResidualBckFuncError[k]                      = new Double_t[fNBinsPt];
        fMesonYieldsPerEventError[k]                                = new Double_t[fNBinsPt];
    }

    for (Int_t m = 0; m < 3; m++){
        fMesonYieldsResBckOtherFunc[m]                              = new Double_t[fNBinsPt];
        fMesonYieldsResBckOtherFuncError[m]                         = new Double_t[fNBinsPt];
    }

    // initialize variable for different integration windows: normal, wide, narrow
    for (Int_t k = 0; k < 3; k++ ){
        // initialize pt arrays for integration ranges
        fMesonTrueIntRange[k]                                        = new Double_t[2];
        fMesonTrueIntReweightedRange[k]                              = new Double_t[2];
        fMesonTrueIntUnweightedRange[k]                              = new Double_t[2];
        // initialize pt arrays for mass windows
        fMassWindowHigh[k]                                           = new Double_t[fNBinsPt];
        fMassWindowLow[k]                                            = new Double_t[fNBinsPt];
        // initialize pt arrays for S/B & Significance
        fMesonSBdefault[k]                                           = new Double_t[fNBinsPt];
        fMesonSigndefault[k]                                         = new Double_t[fNBinsPt];
        fMesonSBdefaultError[k]                                      = new Double_t[fNBinsPt];
        fMesonSigndefaultError[k]                                    = new Double_t[fNBinsPt];

        // initialize pt arrays for reconstructed validated yields
        fMesonTrueYields[k]                                          = new Double_t[fNBinsPt];
        fMesonTrueYieldsReweighted[k]                                = new Double_t[fNBinsPt];
        fMesonTrueYieldsUnweighted[k]                                = new Double_t[fNBinsPt];
        fMesonTrueYieldsError[k]                                     = new Double_t[fNBinsPt];
        fMesonTrueYieldsReweightedError[k]                           = new Double_t[fNBinsPt];
        fMesonTrueYieldsUnweightedError[k]                           = new Double_t[fNBinsPt];

        // initialize pt arrays for reconstructed validated yields from secondaries
        for (Int_t j = 0; j < 4; j++){
            fMesonTrueSecYields[k][j]                               = new Double_t[fNBinsPt];
            fMesonTrueSecYieldsError[k][j]                          = new Double_t[fNBinsPt];
        }
    }

    // initialize pt arrays for additional yields in different fixed windows
    fMesonTrueYieldsDC                                              = new Double_t[fNBinsPt];
    fMesonTrueYieldsDCError                                         = new Double_t[fNBinsPt];
    fMesonTrueYieldFixedWindow                                      = new Double_t[fNBinsPt];
    fMesonTrueYieldGammaFixedWindow                                 = new Double_t[fNBinsPt];
    fMesonTrueYieldGammaConvGammaFixedWindow                        = new Double_t[fNBinsPt];
    fMesonTrueYieldConvGammaConvGammaFixedWindow                    = new Double_t[fNBinsPt];
    fMesonTrueYieldErrorFixedWindow                                 = new Double_t[fNBinsPt];
    fMesonTrueYieldGammaErrorFixedWindow                            = new Double_t[fNBinsPt];
    fMesonTrueYieldGammaConvGammaErrorFixedWindow                   = new Double_t[fNBinsPt];
    fMesonTrueYieldConvGammaConvGammaErrorFixedWindow               = new Double_t[fNBinsPt];

    fMesonLambdaTailpar                                             = new Double_t[fNBinsPt];
    fMesonLambdaTailparError                                        = new Double_t[fNBinsPt];
    fMesonLambdaTailMCpar                                           = new Double_t[fNBinsPt];
    fMesonLambdaTailMCparError                                      = new Double_t[fNBinsPt];
    fMesonAmplitudepar                                              = new Double_t[fNBinsPt];
    fMesonAmplitudeparError                                         = new Double_t[fNBinsPt];
    fMesonSigmapar                                                  = new Double_t[fNBinsPt];
    fMesonSigmaparError                                             = new Double_t[fNBinsPt];
    fMesonTrueSigmapar                                              = new Double_t[fNBinsPt];
    fMesonTrueSigmaparError                                         = new Double_t[fNBinsPt];
    fMesonResidualBGlin                                             = new Double_t[fNBinsPt];
    fMesonResidualBGlinError                                        = new Double_t[fNBinsPt];
    fMesonResidualBGcon                                             = new Double_t[fNBinsPt];
    fMesonResidualBGconError                                        = new Double_t[fNBinsPt];
    for (Int_t m = 0; m < 4; m++){
        fMesonChi2[m]                                               = new Double_t[fNBinsPt];
    }

    // initialize pt-arrays for different mass & width fitting procedures
    fMesonMass                                                      = new Double_t[fNBinsPt];
    fMesonMassError                                                 = new Double_t[fNBinsPt];
    fMesonFWHM                                                      = new Double_t[fNBinsPt];
    fMesonFWHMError                                                 = new Double_t[fNBinsPt];
    fMesonTrueMass                                                  = new Double_t[fNBinsPt];
    fMesonTrueMassReweighted                                        = new Double_t[fNBinsPt];
    fMesonTrueMassUnweighted                                        = new Double_t[fNBinsPt];
    fMesonTrueFWHM                                                  = new Double_t[fNBinsPt];
    fMesonTrueFWHMReweighted                                        = new Double_t[fNBinsPt];
    fMesonTrueFWHMUnweighted                                        = new Double_t[fNBinsPt];
    fMesonTrueMassError                                             = new Double_t[fNBinsPt];
    fMesonTrueMassReweightedError                                   = new Double_t[fNBinsPt];
    fMesonTrueMassUnweightedError                                   = new Double_t[fNBinsPt];
    fMesonTrueFWHMError                                             = new Double_t[fNBinsPt];
    fMesonTrueFWHMReweightedError                                   = new Double_t[fNBinsPt];
    fMesonTrueFWHMUnweightedError                                   = new Double_t[fNBinsPt];
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

    // nitialize pt-arrays for different mass & width fitting procedures, normalization at the left of the peak
    fMesonMassLeft                                                  = new Double_t[fNBinsPt];
    fMesonFWHMLeft                                                  = new Double_t[fNBinsPt];
    fMesonMassLeftError                                             = new Double_t[fNBinsPt];
    fMesonFWHMLeftError                                             = new Double_t[fNBinsPt];

    // initialize pt-arrays for validated Significance and S/B
    fMesonTrueSB                                                    = new Double_t[fNBinsPt];
    fMesonTrueSign                                                  = new Double_t[fNBinsPt];
    fMesonTrueSBError                                               = new Double_t[fNBinsPt];
    fMesonTrueSignError                                             = new Double_t[fNBinsPt];


    // initialize different histo-pt-arrays for validated meson quantities
    fHistoMappingTrueMesonInvMassPtBins                             = new TH1D*[fNBinsPt];
    fHistoMappingTrueMesonDCInvMassPtBins                           = new TH1D*[fNBinsPt];
    fHistoMappingTrueFullMesonInvMassPtBins                         = new TH1D*[fNBinsPt];
    fHistoMappingTrueMesonInvMassPtReweightedBins                   = new TH1D*[fNBinsPt];
    fHistoMappingTrueMesonInvMassPtUnweightedBins                   = new TH1D*[fNBinsPt];
    fHistoMappingTrueGGBckInvMassPtBins                             = new TH1D*[fNBinsPt];
    fHistoMappingTrueContBckInvMassPtBins                           = new TH1D*[fNBinsPt];
    fHistoMappingTrueMesonContainedInvMassPtBins                    = new TH1D*[fNBinsPt];
    fHistoMappingTrueAsymEClusInvMassPtBins                         = new TH1D*[fNBinsPt];
    fHistoMappingTrueAllBckInvMassPtBins                            = new TH1D*[fNBinsPt];
    for (Int_t i = 0; i < 4; i++){
        fHistoMappingTrueSecMesonInvMassPtBins[i]                   = new TH1D*[fNBinsPt];
    }
    fHistoMappingTrueMesonCaloPhotonInvMassPtBins                   = new TH1D*[fNBinsPt];
    fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins               = new TH1D*[fNBinsPt];
    fHistoMappingTrueMesonCaloElectronInvMassPtBins                 = new TH1D*[fNBinsPt];
    fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins            = new TH1D*[fNBinsPt];
    fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins    = new TH1D*[fNBinsPt];
    fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins          = new TH1D*[fNBinsPt];

    // initialize arrays for BG weighting
    fHistoWeightsBGZbinVsMbin                                       = new TH2F*[fNBinsPt];
    fHistoFillPerEventBGZbinVsMbin                                  = new TH2F*[fNBinsPt];
    fHistoWeightsBGZbinVsPsibin                                     = new TH2F*[fNBinsPt];
    fHistoFillPerEventBGZbinVsPsibin                                = new TH2F*[fNBinsPt];

    // initial  histo-pt-arrays for  data inv mass
    fHistoMappingGGInvMassPtBin                                     = new TH1D*[fNBinsPt];
    fHistoMappingBackInvMassPtBin                                   = new TH1D*[fNBinsPt];
    fHistoMappingBackNormInvMassPtBin                               = new TH1D*[fNBinsPt];
    fHistoMappingSignalInvMassPtBin                                 = new TH1D*[fNBinsPt];
    fHistoMappingSignalRemainingBGSubInvMassPtBin                   = new TH1D*[fNBinsPt];
    fHistoMappingRemainingBGInvMassPtBin                            = new TH1D*[fNBinsPt];
    fHistoMappingRatioSBInvMassPtBin                                = new TH1D*[fNBinsPt];
    fHistoMappingBackNormAndRemainingBGInvMassPtBin                 = new TH1D*[fNBinsPt];

    // initial  fit-pt-arrays for  data inv mass
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

    for (Int_t m = 0; m < 3; m++){
        fFitSignalWithOtherBGInvMassPtBin[m]                        = new TF1*[fNBinsPt];
        fFitBckOtherInvMassPtBin[m]                                 = new TF1*[fNBinsPt];
    }

    // Histograms for normalization on the left of the peak
    fHistoMappingBackNormInvMassLeftPtBin                           = new TH1D*[fNBinsPt];
    fHistoMappingSignalInvMassLeftPtBin                             = new TH1D*[fNBinsPt];
    fHistoMappingSignalRemainingBGSubInvMassLeftPtBin               = new TH1D*[fNBinsPt];
    fHistoMappingRemainingBGInvMassLeftPtBin                        = new TH1D*[fNBinsPt];

    fFitInvMassLeftPtBin                                            = new TF1*[fNBinsPt];
    fFitRemainingBGInvMassLeftPtBin                                 = new TF1*[fNBinsPt];
    fFitSignalPeakPosInvMassLeftPtBin                               = new TF1*[fNBinsPt];
    fFitBckInvMassLeftPtBin                                         = new TF1*[fNBinsPt];

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
    for(Int_t i = 0;i<fNBinsPt; i++){
        fHistoMappingTrueMesonInvMassPtBins[i]                              = NULL;
        fHistoMappingTrueMesonDCInvMassPtBins[i]                            = NULL;
        fHistoMappingTrueFullMesonInvMassPtBins[i]                          = NULL;
        fHistoMappingTrueMesonInvMassPtReweightedBins[i]                    = NULL;
        fHistoMappingTrueMesonInvMassPtUnweightedBins[i]                    = NULL;
        fHistoMappingTrueGGBckInvMassPtBins[i]                              = NULL;
        fHistoMappingTrueMesonContainedInvMassPtBins[i]                     = NULL;
        fHistoMappingTrueAsymEClusInvMassPtBins[i]                          = NULL;
        fHistoMappingTrueContBckInvMassPtBins[i]                            = NULL;
        fHistoMappingTrueAllBckInvMassPtBins[i]                             = NULL;
        for (Int_t j = 0; j < 4; j++ ){
            fHistoMappingTrueSecMesonInvMassPtBins[j][i]                    = NULL;
        }
        fHistoMappingTrueMesonCaloPhotonInvMassPtBins[i]                    = NULL;
        fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[i]                = NULL;
        fHistoMappingTrueMesonCaloElectronInvMassPtBins[i]                  = NULL;
        fHistoMappingTrueMesonCaloMergedClusterInvMassPtBins[i]             = NULL;
        fHistoMappingTrueMesonCaloMergedClusterPartConvInvMassPtBins[i]     = NULL;
        fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[i]           = NULL;

        fHistoMappingGGInvMassPtBin[i]                                      = NULL;
        fHistoMappingBackInvMassPtBin[i]                                    = NULL;
        fHistoMappingBackNormInvMassPtBin[i]                                = NULL;
        fHistoMappingSignalInvMassPtBin[i]                                  = NULL;
        fHistoMappingSignalRemainingBGSubInvMassPtBin[i]                    = NULL;
        fHistoMappingRemainingBGInvMassPtBin[i]                             = NULL;
        fHistoMappingRatioSBInvMassPtBin[i]                                 = NULL;
        fHistoMappingBackNormAndRemainingBGInvMassPtBin[i]                  = NULL;

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

        for (Int_t m = 0; m < 3; m++){
            fFitSignalWithOtherBGInvMassPtBin[m][i]                         = NULL;
            fFitBckOtherInvMassPtBin[m][i]                                  = NULL;
        }

        // Histograms for normalization on the left of the peak
        fHistoMappingBackNormInvMassLeftPtBin[i]                            = NULL;
        fHistoMappingSignalInvMassLeftPtBin[i]                              = NULL;
        fHistoMappingSignalRemainingBGSubInvMassLeftPtBin[i]                = NULL;
        fHistoMappingRemainingBGInvMassLeftPtBin[i]                         = NULL;

        fFitInvMassLeftPtBin[i]                                             = NULL;
        fFitRemainingBGInvMassLeftPtBin[i]                                  = NULL;
        fFitSignalPeakPosInvMassLeftPtBin[i]                                = NULL;
        fFitBckInvMassLeftPtBin[i]                                          = NULL;

        fFitSignalGaussianInvMassPtBin[i]                                   = NULL;
        fFitTrueSignalGaussianInvMassPtBin[i]                               = NULL;
    }
}


//****************************************************************************
//****** Initializiation of MC histogram names according  *****
//****************************************************************************
void SetCorrectMCHistogrammNames(TString mesonType){
    cout << "standard MC chosen" << endl;

    // MC histograms primaries
    // Pi0
        ObjectNameMCPi0Acc                  = "MC_Pi0InAcc_Pt";
        ObjectNameMCPi0AccWOWeights         = "MC_Pi0WOWeightInAcc_Pt";
        if(fMode == 4 && (fEnergyFlag.CompareTo("pPb_5.023TeV") == 0|| fEnergyFlag.CompareTo("pPb_5.023TeVRun2") == 0) ) ObjectNameMCPi0AccWOWeights = "MC_Pi0InAcc_Pt";
        if(fMode == 4 || fMode == 12 || fMode == 5)
            ObjectNameMCPi0AccWOEvtWeights    = "MC_Pi0WOEvtWeightInAcc_Pt";
        else
            ObjectNameMCPi0AccWOEvtWeights    = "MC_Pi0_WOEventWeightsInAcc_Pt";
        ObjectNameMCPi0                     = "MC_Pi0_Pt";
        ObjectNameMCPi0WOWeights            = "MC_Pi0_WOWeights_Pt";
        ObjectNameMCPi0WOEvtWeights         = "MC_Pi0_WOEventWeights_Pt";
    // Eta
        if( fModeHeavy<100 ) {
            ObjectNameMCEtaAcc                  = "MC_EtaInAcc_Pt";
            ObjectNameMCEtaAccWOWeights         = "MC_EtaWOWeightInAcc_Pt";
            if(fMode == 4 && (fEnergyFlag.CompareTo("pPb_5.023TeV") == 0|| fEnergyFlag.CompareTo("pPb_5.023TeVRun2") == 0) )
                ObjectNameMCEtaAccWOWeights = "MC_EtaInAcc_Pt";
            if(fMode == 4 || fMode == 12 || fMode == 5)
                ObjectNameMCEtaAccWOEvtWeights    = "MC_EtaWOEvtWeightInAcc_Pt";
            else
                ObjectNameMCEtaAccWOEvtWeights    = "MC_Eta_WOEventWeightsInAcc_Pt";
            ObjectNameMCEta                     = "MC_Eta_Pt";
            ObjectNameMCEtaWOWeights            = "MC_Eta_WOWeights_Pt";
            ObjectNameMCEtaWOEvtWeights         = "MC_Eta_WOEventWeights_Pt";
        } else {
            ObjectNameMCEtaAcc                  = "MC_MesonInAcc_Pt";
            ObjectNameMCEtaAccWOWeights         = "MC_MesonWOWeightInAcc_Pt";
            if(fMode == 4 && (fEnergyFlag.CompareTo("pPb_5.023TeV") == 0|| fEnergyFlag.CompareTo("pPb_5.023TeVRun2") == 0) )
                ObjectNameMCEtaAccWOWeights = "MC_MesonInAcc_Pt";
            if(fMode == 4 || fMode == 12 || fMode == 5)
                ObjectNameMCEtaAccWOEvtWeights    = "MC_MesonWOEvtWeightInAcc_Pt";
            else
                ObjectNameMCEtaAccWOEvtWeights    = "MC_Meson_WOEventWeightsInAcc_Pt";
            ObjectNameMCEta                = "MC_Meson_Pt";
            ObjectNameMCEtaWOWeights       = "MC_Meson_WOWeights_Pt";
            ObjectNameMCEtaWOEvtWeights    = "MC_Meson_WOEventWeights_Pt";
        }
    // EtaPrime / Heavy meson
        ObjectNameMCEtaPrimeAcc                  = "MC_MesonInAcc_Pt";
        ObjectNameMCEtaPrimeAccWOWeights         = "MC_MesonWOWeightInAcc_Pt";
        if(fMode == 4 && (fEnergyFlag.CompareTo("pPb_5.023TeV") == 0|| fEnergyFlag.CompareTo("pPb_5.023TeVRun2") == 0) )
            ObjectNameMCEtaPrimeAccWOWeights = "MC_MesonInAcc_Pt";
        if(fMode == 4 || fMode == 12 || fMode == 5)
            ObjectNameMCEtaPrimeAccWOEvtWeights    = "MC_MesonWOEvtWeightInAcc_Pt";
        else
            ObjectNameMCEtaPrimeAccWOEvtWeights    = "MC_Meson_WOEventWeightsInAcc_Pt";
        ObjectNameMCEtaPrime                = "MC_Meson_Pt";
        ObjectNameMCEtaPrimeWOWeights       = "MC_Meson_WOWeights_Pt";
        ObjectNameMCEtaPrimeWOEvtWeights    = "MC_Meson_WOEventWeights_Pt";

    // MC histograms secondaries
    if (mesonType.Contains("Pi0")){
        ObjectNameMCSecPi0Acc           = "MC_SecPi0InAcc_Pt_Source";
        ObjectNameMCSecPi0              = "MC_SecPi0_Pt_Source";
    }

    // reconstructed validated histograms main
    ObjectNameTrue                      = "ESD_TruePrimaryMother_InvMass_Pt";
    ObjectNameTrueFull                  = "ESD_TrueMother_InvMass_Pt";
    ObjectNameTrueWOWeights             = "ESD_TruePrimaryMotherW0Weights_InvMass_Pt";
    ObjectNameProfileWeights            = "ESD_TruePrimaryMotherWeights_InvMass_Pt";
    ObjectNameTrueSec                   = "ESD_TrueSecondaryMother_InvMass_Pt";;
    ObjectNameTrueSecFromK0S            = "ESD_TrueSecondaryMotherFromK0s_InvMass_Pt";
    ObjectNameTrueSecFromK0L            = "ESD_TrueSecondaryMotherFromK0l_InvMass_Pt";
    ObjectNameTrueSecFromLambda         = "ESD_TrueSecondaryMotherFromLambda_InvMass_Pt";

    // reconstructed validated histograms additional conversions
    ObjectNameTrueGGBck                 = "ESD_TrueBckGG_InvMass_Pt";
    ObjectNameTrueContBck               = "ESD_TrueBckCont_InvMass_Pt";
    ObjectNameTrueBckFullMesonContained = "ESD_TrueBckFullMesonContained_InvMass_Pt";
    ObjectNameTrueBckAsymEClus          = "ESD_TrueBckAsymEClus_InvMass_Pt";
    ObjectNameTrueAllBck                = "ESD_TrueAllCont_InvMass_Pt";

    // reconstructed validated histograms additional calorimeters
    ObjectNameTrueCaloPhoton            = "ESD_TrueMotherCaloPhoton_InvMass_Pt";
    ObjectNameTrueCaloConvPhoton        = "ESD_TrueMotherCaloConvertedPhoton_InvMass_Pt";
    ObjectNameTrueMixedCaloConvPhoton   = "ESD_TrueMotherCaloMixedPhotonConvertedPhoton_InvMass_Pt";
    ObjectNameTrueCaloMerged            = "ESD_TrueMotherCaloMergedCluster_InvMass_Pt";
    ObjectNameTrueCaloMergedPartConv    = "ESD_TrueMotherCaloMergedClusterPartConv_InvMass_Pt";
    if (fMode > 1 && fMode !=9){
        ObjectNameTrueCaloPhoton            = Form("ESD_True%sCaloPhoton_InvMass_Pt", mesonType.Data());
        ObjectNameTrueCaloConvPhoton        = Form("ESD_True%sCaloConvertedPhoton_InvMass_Pt", mesonType.Data());
        ObjectNameTrueMixedCaloConvPhoton   = Form("ESD_True%sCaloMixedPhotonConvertedPhoton_InvMass_Pt", mesonType.Data());
        ObjectNameTrueCaloMerged            = Form("ESD_True%sCaloMergedCluster_InvMass_Pt", mesonType.Data());
        ObjectNameTrueCaloMergedPartConv    = Form("ESD_True%sCaloMergedClusterPartConv_InvMass_Pt", mesonType.Data());

        cout <<    ObjectNameTrueCaloPhoton.Data()          << "\t" <<
                    ObjectNameTrueCaloConvPhoton.Data()       << "\t" <<
                    ObjectNameTrueMixedCaloConvPhoton.Data()    << "\t" <<
                    ObjectNameTrueCaloMerged.Data()          << "\t" <<
                    ObjectNameTrueCaloMergedPartConv.Data() << endl;

        // set correct names for Calo modes
        ObjectNameTrue                      = Form("ESD_TruePrimary%s_InvMass_Pt", mesonType.Data());
        ObjectNameTrueFull                  = Form("ESD_True%s_InvMass_Pt", mesonType.Data());
        ObjectNameTrueWOWeights             = Form("ESD_TruePrimary%sW0Weights_InvMass_Pt", mesonType.Data());
        ObjectNameProfileWeights            = Form("ESD_TruePrimary%sWeights_InvMass_Pt", mesonType.Data());
        ObjectNameTrueSec                   = Form("ESD_TrueSecondary%s_InvMass_Pt", mesonType.Data());
        ObjectNameTrueSecFromK0S            = Form("ESD_TrueSecondary%sFromK0s_InvMass_Pt", mesonType.Data());
        ObjectNameTrueSecFromK0L            = Form("ESD_TrueSecondary%sFromK0l_InvMass_Pt", mesonType.Data());

        ObjectNameTrueSecFromLambda         = Form("ESD_TrueSecondary%sFromLambda_InvMass_Pt", mesonType.Data());
    }

    // Correction for MC histograms in EtaPrime
    if( mesonType.EqualTo("EtaPrime") || fModeHeavy>=100 ){
        ObjectNameTrue                      = "ESD_TruePrimaryMeson_InvMass_Pt";
        ObjectNameTrueFull                  = "ESD_TrueMeson_InvMass_Pt";
        ObjectNameTrueWOWeights             = "ESD_TruePrimaryMesonW0Weights_InvMass_Pt";
        ObjectNameProfileWeights            = "ESD_TruePrimaryMesonWeights_InvMass_Pt";
    }

    // input spectra for secondary particles
    ObjectNameK0sRecPi0                 = "TrueK0sWithPi0Daughter_MCPt";
    ObjectNameLambdaRecPi0              = "TrueLambdaWithPi0Daughter_MCPt";


    // double counting histograms
    ObjectNameDCMesonInvMassPt          = Form("ESD_TrueDoubleCount%s_InvMass_Pt", mesonType.Data());
    ObjectNameDCGammaClusPt             = "TrueDoubleCountClusterGamma_Pt";
    ObjectNameMesonMultipleCount        = Form("ESD_TrueMultipleCount%s", mesonType.Data());
    ObjectNameGammaClusMultipleCount    = "ESD_TrueMultipleCountClusterGamma";

    return;
}

//****************************************************************************
//Load secondary neutral pions from toy MC file as generated from data spectra
// - put them in proper scaling
// - rebin them according to current pi0 binning
//****************************************************************************
Bool_t LoadSecondaryPionsFromExternalFile(){
    ifstream in("ToyMCOutputs.txt");

    // number of ToyMC inputs
    Int_t nrOfToyMCInput        = 0;
    // number of triggers which are really used for the respective analysis
    TString nameToyMCInputs[10];
    Bool_t usedInput[10]        = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE };
    while(!in.eof() && nrOfToyMCInput<10 ){
        in >> nameToyMCInputs[nrOfToyMCInput];
        cout<< nameToyMCInputs[nrOfToyMCInput] << endl;
        nrOfToyMCInput++;
    }
    if (nrOfToyMCInput==0 || (nrOfToyMCInput==1 && nameToyMCInputs[nrOfToyMCInput].CompareTo("") == 0) )
        return kFALSE;

    Int_t nSecInputHistsFound       = 0;
    for (Int_t j = 0; j < 3; j++){
        cout << "searching for input for " << nameSecondaries[j].Data() << endl;
        Bool_t foundSourceFile      = kFALSE;
        TString nameSourceFile      = "";
        // find correct source file only first one for respective particle will be used from the list
        for (Int_t f= 0; (f< nrOfToyMCInput && !foundSourceFile); f++){
            if (!usedInput[f]){
                if ( nameToyMCInputs[f].Contains(nameSecondaries[j].Data()) ){
                    foundSourceFile = kTRUE;
                    nameSourceFile  = nameToyMCInputs[f];
                }
            }
        }
        if (foundSourceFile){
            cout << "found correct input: " << nameSourceFile.Data() << endl;
            cout << "trying to find " << Form("histoPiZeroDaughtersPt_InRap_%1.2f", fYMaxMeson/2) << endl;
            fFileToyMCInput[j]                  = new TFile(nameSourceFile.Data());
            fHistoYieldExternSecInput[j]         = (TH1D*)fFileToyMCInput[j]->Get(Form("histoPiZeroDaughtersPt_InRap_%1.2f", fYMaxMeson/2));
            if (fHistoYieldExternSecInput[j]){
                nSecInputHistsFound++;
                fHistoYieldExternSecInput[j]->Sumw2();
                fHistoYieldExternSecInput[j]->SetName(Form("histoSecPi0YieldFrom%s_FromToy_orgBinning",nameSecondaries[j].Data()));

                cout << "found it, rebinning" << endl;
                fHistoYieldExternSecInputReb[j]  =  (TH1D*)fHistoYieldExternSecInput[j]->Rebin(fNBinsPt,Form("histoSecPi0YieldFrom%s_FromToy",nameSecondaries[j].Data()),fBinsPt); // Proper bins in Pt
                if (fHistoYieldExternSecInputReb[j]){
                    fHistoYieldExternSecInputReb[j]->Divide(fDeltaPt);
                    fHistoYieldExternSecInput[j]->Scale(1./fHistoYieldExternSecInput[j]->GetBinWidth(1));
                    cout << "that worked" << endl;
                }
            } else {
                cout << "file didn't contain proper histo" << endl;
            }
        } else {
            cout << "could not find correct input file" << endl;
        }
    }
    // cout << nSecInputHistsFound << endl;
    if (nSecInputHistsFound == 0)
        return kFALSE;

    return kTRUE;
}


//******************************************************************************
// Load secondary neutral pions from cocktail file
// - put them in proper scaling
// - rebin them according to current pi0 binning
//******************************************************************************
Bool_t LoadSecondaryPionsFromCocktailFile(TString cutSelection, TString optionEnergy){
    TString nameCocktailFile                     = Form("%s/%s/SecondaryPi0%s_%.2f_%s.root",cutSelection.Data(),optionEnergy.Data(),fPeriodFlag.Data(),fYMaxMeson/2,cutSelection.Data());
    fFileCocktailInput                           = new TFile(nameCocktailFile.Data());
    if (!fFileCocktailInput->IsZombie()){
        cout << "found correct input: " << nameCocktailFile.Data() << endl;

        // secondary neutral pions
        for (Int_t j = 0; j < 3; j++){
            cout << "trying to find " << Form("Pi0_From_%s_Pt_OrBin", nameSecondariesCocktail[j].Data()) << endl;

            fHistoYieldExternSecInput[j]           = (TH1D*)fFileCocktailInput->Get(Form("Pi0_From_%s_Pt_OrBin", nameSecondariesCocktail[j].Data()));
            if (fHistoYieldExternSecInput[j]){
                fHistoYieldExternSecInput[j]->Sumw2();
                fHistoYieldExternSecInput[j]->SetName(Form("histoSecPi0YieldFrom%s_FromCocktail_orgBinning",nameSecondaries[j].Data())); // Proper bins in Pt

                fHistoYieldExternSecInputReb[j]     = (TH1D*)fHistoYieldExternSecInput[j]->Rebin(fNBinsPt,Form("histoSecPi0YieldFrom%s_FromCocktail",nameSecondaries[j].Data()),fBinsPt); // Proper bins in Pt
                if (fHistoYieldExternSecInputReb[j]){
                    fHistoYieldExternSecInputReb[j]->Divide(fDeltaPt);
                    fHistoYieldExternSecInput[j]->Scale(1./fHistoYieldExternSecInput[j]->GetBinWidth(1));
                }
            } else {
                cout << "file didn't contain " << Form("Pi0_From_%s_Pt_OrBin", nameSecondariesCocktail[j].Data()) << endl;
            }
        }

        // neutral pions from resonance decays
        for (Int_t j = 0; j < 15; j++) {
            cout << "trying to find " << Form("Pi0_From_%s_Pt_OrBin", nameResonanceFeedDownContributionsCocktail[j].Data()) << endl;

            fHistoYieldExternResonanceFeedDownInput[j]           = (TH1D*)fFileCocktailInput->Get(Form("Pi0_From_%s_Pt_OrBin", nameResonanceFeedDownContributionsCocktail[j].Data()));
            if (fHistoYieldExternResonanceFeedDownInput[j]) {
                fHistoYieldExternResonanceFeedDownInput[j]->Sumw2();
                fHistoYieldExternResonanceFeedDownInput[j]->SetName(Form("histoResonanceFeedDownPi0YieldFrom%s_FromCocktail_orgBinning",nameResonanceFeedDownContributions[j].Data()));

                fHistoYieldExternResonanceFeedDownInputReb[j]    = (TH1D*)fHistoYieldExternResonanceFeedDownInput[j]->Rebin(fNBinsPt,Form("histoResonanceFeedDownPi0YieldFrom%s_FromCocktail",nameResonanceFeedDownContributions[j].Data()),fBinsPt);
                if (fHistoYieldExternResonanceFeedDownInputReb[j]) {
                    fHistoYieldExternResonanceFeedDownInputReb[j]->Divide(fDeltaPt);
                    fHistoYieldExternResonanceFeedDownInput[j]->Scale(1./fHistoYieldExternResonanceFeedDownInput[j]->GetBinWidth(1));
                }
            } else {
                cout << "file didn't contain " << Form("Pi0_From_%s_Pt_OrBin", nameResonanceFeedDownContributionsCocktail[j].Data()) << endl;
            }
        }

        return kTRUE;
    }
    return kFALSE;
}


//****************************************************************************
//******************** Projection out of 2D in X *****************************
//****************************************************************************
TH1D* FillProjectionX (TH2* dummy2D, TString name, Double_t minY, Double_t maxY, Int_t rebin){
    TH1D* dummy1D           = new TH1D(name.Data(), name.Data(), dummy2D->GetNbinsX(), 0., dummy2D->GetXaxis()->GetBinUpEdge(dummy2D->GetNbinsX()));
    dummy1D->Sumw2();
    Int_t startBin          = dummy2D->GetYaxis()->FindBin(minY+0.001);
    Int_t endBin            = dummy2D->GetYaxis()->FindBin(maxY-0.001);
    dummy2D->ProjectionX(name.Data(),startBin,endBin,"e");
    dummy1D                 = (TH1D*)gDirectory->Get(name.Data());
    if(rebin>1){
        dummy1D->Rebin(rebin);
    }
    return dummy1D;
}

//****************************************************************************
//*************** Check if histo already exists if not clear *****************
//****************************************************************************
void CheckForNULLForPointer(TH1D* dummy1){
    if(dummy1!= NULL){
        delete dummy1;
        dummy1           = NULL;
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
void FillMassHistosArray(TH2D* fGammaGammaInvMassVSPtDummy) {
    fGammaGammaInvMassVSPtDummy->Sumw2();
    fFittingHistMidPtSignal     = FillProjectionX(fGammaGammaInvMassVSPtDummy, "Mapping_GG_InvMass_MidPt", fMidPt[0], fMidPt[1], fNRebin[4]);
    fMesonFullPtSignal          = FillProjectionX(fGammaGammaInvMassVSPtDummy, "Mapping_GG_InvMass_FullPt", fFullPt[0], fFullPt[1], fNRebin[4]);

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
//         cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonInvMassPtBins[iPt]->GetEntries() << endl;
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
//         cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueFullMesonInvMassPtBins[iPt]->GetEntries() << endl;
        fHistoMappingTrueFullMesonInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueFullMesonInvMassPtBins[iPt]->SetLineColor(2);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated true mesons **********************************************
//****************************************************************************
void FillMassMCTrueMesonDCHistosArray(TH2D* fHistoTrueMesonInvMassVSPtFill) {
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrue  = Form("Mapping_TrueMesonDC_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueMesonDCInvMassPtBins[iPt]);
        fHistoMappingTrueMesonDCInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueMesonInvMassVSPtFill, fNameHistoTrue, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
//         cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonDCInvMassPtBins[iPt]->GetEntries() << endl;
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
//         cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonCaloPhotonInvMassPtBins[iPt]->GetEntries() << endl;
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
//         cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonCaloConvPhotonInvMassPtBins[iPt]->GetEntries() << endl;
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
//         cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonCaloElectronInvMassPtBins[iPt]->GetEntries() << endl;
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
//         cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->GetEntries() << endl;
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
//         cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonInvMassPtUnweightedBins[iPt]->GetEntries() << endl;
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
//******* validated BG *******************************************************
//****************************************************************************
void FillMassMCTrueFullMesonContainedHistosArray(TH2D* fHistoTrueFullMesonContainedInvMassVSPtFill) {
    fHistoTrueFullMesonContainedInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrueMesonContained = Form("Mapping_TrueFullMesonContained_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueMesonContainedInvMassPtBins[iPt]);
        fHistoMappingTrueMesonContainedInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueFullMesonContainedInvMassVSPtFill, fNameHistoTrueMesonContained, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        fHistoMappingTrueMesonContainedInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueMesonContainedInvMassPtBins[iPt]->SetLineColor(2);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated BG *******************************************************
//****************************************************************************
void FillMassMCTrueAsymEClusHistosArray(TH2D* fHistoTrueAsymEClusInvMassVSPtFill) {
    fHistoTrueAsymEClusInvMassVSPtFill->Sumw2();
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoTrueAsymEClus = Form("Mapping_TrueAsymEClus_InvMass_in_Pt_Bin%02d", iPt);
        CheckForNULLForPointer(fHistoMappingTrueAsymEClusInvMassPtBins[iPt]);
        fHistoMappingTrueAsymEClusInvMassPtBins[iPt]=  FillProjectionX(fHistoTrueAsymEClusInvMassVSPtFill, fNameHistoTrueAsymEClus, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
        fHistoMappingTrueAsymEClusInvMassPtBins[iPt]->SetLineWidth(1);
        fHistoMappingTrueAsymEClusInvMassPtBins[iPt]->SetLineColor(2);
    }
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//******* validated secondary for mesons from any source *********************
//****************************************************************************
void FillMassMCTrueSecMesonHistosArray(TH2D** fHistoTrueSecMesonInvMassVSPtFill) {

    for (Int_t j = 0; j < 4; j++){
        if (fHistoTrueSecMesonInvMassVSPtFill[j]){
            fHistoTrueSecMesonInvMassVSPtFill[j]->Sumw2();
            for(Int_t iPt=fStartPtBin; iPt<fNBinsPt; iPt++){
                fNameHistoTrueSec                               = Form("Mapping_TrueSecPi0From%s_InvMass_in_Pt_Bin%02d",nameSecondaries[j].Data(), iPt);
                CheckForNULLForPointer(fHistoMappingTrueSecMesonInvMassPtBins[j][iPt]);
                fHistoMappingTrueSecMesonInvMassPtBins[j][iPt]  =  FillProjectionX(fHistoTrueSecMesonInvMassVSPtFill[j], fNameHistoTrueSec, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt]);
                fHistoMappingTrueSecMesonInvMassPtBins[j][iPt]->SetLineWidth(1);
                fHistoMappingTrueSecMesonInvMassPtBins[j][iPt]->SetLineColor(2);
            }
        } else {
            cout << "reconstructed histo for: " << nameSecondaries[j].Data() << "was not available" << endl;
        }
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

    // create histos for different integration windows: normal, wide, narrow, left, left wide, left narrow
    for (Int_t k = 0; k < 6; k++){
        // reconstructed yields in different integration windows
        fHistoYieldMeson[k]                    = new TH1D(Form("histoYieldMeson%s",nameIntRange[k].Data()),"",fNBinsPt,fBinsPt);
        fHistoYieldMeson[k]->Sumw2();
        fHistoYieldMesonPerEvent[k]            = new TH1D(Form("histoYieldMeson%sPerEvent",nameIntRange[k].Data()),"",fNBinsPt,fBinsPt);
        fHistoYieldMesonPerEvent[k]->Sumw2();
    }

    // yield assumptions with different backgrounds
    for (Int_t k = 0; k < 6; k++){
        fHistoYieldDiffBck[k]                    = new TH1D(Form("histoYieldMesonDiffBck%s",nameIntBck[k].Data()),"",fNBinsPt,fBinsPt);
        fHistoYieldDiffBck[k]->Sumw2();
        fHistoYieldDiffBckRatios[k]                    = new TH1D(Form("histoYieldMesonDiffBckRatios%s",nameIntBckRatios[k].Data()),"",fNBinsPt,fBinsPt);
        fHistoYieldDiffBckRatios[k]->Sumw2();
    }
    for (Int_t k = 0; k < 3; k++){
        fHistoYieldDiffBckResult[k]                    = new TH1D(Form("histoYieldMesonDiffBckResult%s",nameIntBckResult[k].Data()),"",fNBinsPt,fBinsPt);
        fHistoYieldDiffBckResult[k]->Sumw2();
    }


    // create histos for different integration windows: normal, wide, narrow
    for (Int_t k = 0; k < 3; k++){
        // validated yields
        fHistoYieldTrueMeson[k]                = new TH1D(Form("histoYieldTrueMeson%s",nameIntRange[k].Data()),"",fNBinsPt,fBinsPt);
        fHistoYieldTrueMeson[k]->Sumw2();
        fHistoYieldTrueMesonReweighted[k]      = new TH1D(Form("histoYieldTrueMeson%sReweighted",nameIntRange[k].Data()), "",fNBinsPt,fBinsPt);
        fHistoYieldTrueMesonReweighted[k]->Sumw2();
        fHistoYieldTrueMesonUnweighted[k]      = new TH1D(Form("histoYieldTrueMeson%sUnweighted",nameIntRange[k].Data()), "",fNBinsPt,fBinsPt);
        fHistoYieldTrueMesonUnweighted[k]->Sumw2();

        // Significance and S/B for reconstructed signal
        fHistoSigndefaultMeson[k]              = new TH1D(Form("histoSigndefault%sMeson",nameIntRange[k].Data()),"",fNBinsPt,fBinsPt);
        fHistoSigndefaultMeson[k]->Sumw2();
        fHistoSBdefaultMeson[k]                = new TH1D(Form("histoSBdefault%sMeson",nameIntRange[k].Data()),"",fNBinsPt,fBinsPt);
        fHistoSBdefaultMeson[k]->Sumw2();

        // mass window monitoring histos
        fHistoMassWindowHigh[k]                = new TH1D(Form("histoMassWindow%sHigh",nameIntRange[k].Data()),"",fNBinsPt,fBinsPt);
        fHistoMassWindowHigh[k]->Sumw2();
        fHistoMassWindowLow[k]                 = new TH1D(Form("histoMassWindow%sLow",nameIntRange[k].Data()),"",fNBinsPt,fBinsPt);
        fHistoMassWindowLow[k]->Sumw2();
    }

    fHistoYieldTrueMesonDC              = new TH1D("histoYieldTrueMesonDC","",fNBinsPt,fBinsPt);
    fHistoYieldTrueMesonDC->Sumw2();

    // Mass & Width histos for normalization at right side
    fHistoMassMeson                     = new TH1D("histoMassMeson","",fNBinsPt,fBinsPt);
    fHistoMassMeson->Sumw2();
    fHistoMassGaussianMeson             = new TH1D("histoMassGaussianMeson","",fNBinsPt,fBinsPt);
    fHistoMassGaussianMeson->Sumw2();
    fHistoWidthGaussianMeson            = new TH1D("histoWidthGaussianMeson","",fNBinsPt,fBinsPt);
    fHistoWidthGaussianMeson->Sumw2();
    fHistoFWHMMeson                     = new TH1D("histoFWHMMeson","",fNBinsPt,fBinsPt);
    fHistoFWHMMeson->Sumw2();
    // Mass & Width histos for normalization at left side
    fHistoMassMesonLeft                     = new TH1D("histoMassMesonLeft","",fNBinsPt,fBinsPt);
    fHistoMassMesonLeft->Sumw2();
    fHistoFWHMMesonLeft                     = new TH1D("histoFWHMMesonLeft","",fNBinsPt,fBinsPt);
    fHistoFWHMMesonLeft->Sumw2();


    // Mass & Width histos for real signal normalization at right side
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
    // Significance and S/B for MC
    fHistoTrueSignMeson                 = new TH1D("histoTrueSignMeson","",fNBinsPt,fBinsPt);
    fHistoTrueSignMeson->Sumw2();
    fHistoTrueSBMeson                   = new TH1D("histoTrueSBMeson","",fNBinsPt,fBinsPt);
    fHistoTrueSBMeson->Sumw2();

    // Mass & Width histos for real signal in different reco methods
    if (fAdvancedMesonQA && (fMode == 2 || fMode == 13 || fMode == 3 || fMode == 4 || fMode == 12 || fMode == 5)){
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
    // Mass & Width histos for real signal in different reco methods
    if (fAdvancedMesonQA && (fMode == 4 || fMode == 12 || fMode == 5)){
        fHistoTrueMassMesonMixedCaloConvPhoton =       new TH1D("histoTrueMassMesonMixedCaloConvPhoton","",fNBinsPt,fBinsPt);
        fHistoTrueMassMesonMixedCaloConvPhoton->Sumw2();
        fHistoTrueFWHMMesonMixedCaloConvPhoton =       new TH1D("histoTrueFWHMMesonMixedCaloConvPhoton","",fNBinsPt,fBinsPt);
        fHistoTrueFWHMMesonMixedCaloConvPhoton->Sumw2();
    }

    // monitoring histos
    fHistoLambdaTail                    = new TH1D("histoLambdaTail","",fNBinsPt,fBinsPt);
    fHistoLambdaTail->Sumw2();
    fHistoTrueLambdaTail                = new TH1D("histoTrueLambdaTail","",fNBinsPt,fBinsPt);
    fHistoTrueLambdaTail->Sumw2();
    fHistoAmplitude                     = new TH1D("histoAmplitude","",fNBinsPt,fBinsPt);
    fHistoAmplitude->Sumw2();
    fHistoSigma                         = new TH1D("histoSigma","",fNBinsPt,fBinsPt);
    fHistoSigma->Sumw2();
    fHistoTrueSigma                     = new TH1D("histoTrueSigma","",fNBinsPt,fBinsPt);
    fHistoTrueSigma->Sumw2();

    fHistoResidualBGlin                 = new TH1D("histoResidualBGlin","",fNBinsPt,fBinsPt);
    fHistoResidualBGlin->Sumw2();
    fHistoResidualBGcon                 = new TH1D("histoResidualBGcon","",fNBinsPt,fBinsPt);
    fHistoResidualBGcon->Sumw2();
    fHistoRatioResBGYield               = new TH1D("histoRatioResBGYield","",fNBinsPt,fBinsPt);
    fHistoRatioResBGYield->Sumw2();
    fHistoRatioResBGYieldToSPlusResBG   = new TH1D("histoRatioResBGYieldToSPlusResBG","",fNBinsPt,fBinsPt);
    fHistoRatioResBGYieldToSPlusResBG->Sumw2();
    for (Int_t m = 0; m < 4; m++){
        fHistoChi2[m]                   = new TH1D(Form("histoChi2_%d",m),"",fNBinsPt,fBinsPt);
        fHistoChi2[m]->Sumw2();
        fHistoResBGYield[m]             = new TH1D(Form("histoResBGYield_%d",m),"",fNBinsPt,fBinsPt);
        fHistoResBGYield[m]->Sumw2();
    }

    // create secondary histos
    for (Int_t k = 0; k < 3; k++){
        for (Int_t j = 0; j < 4; j++){
            fHistoYieldTrueSecMeson[k][j]   = new TH1D(Form("histoYieldTrueFrom%sSecMeson%s",nameSecondaries[j].Data(),nameIntRange[k].Data()),"",fNBinsPt,fBinsPt);
            fHistoYieldTrueSecMeson[k][j]->Sumw2();
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

            if (fAdvancedMesonQA && (fMode == 2 || fMode == 13 || fMode == 3 || fMode == 4 || fMode == 12 || fMode == 5)){
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
            if (fAdvancedMesonQA && ( fMode == 4 || fMode == 12 || fMode == 5)){
                fHistoTrueFWHMMesonMixedCaloConvPhoton->SetBinContent(iPt,fMesonTrueFWHMMixedCaloConvPhoton[iPt-1]);
                fHistoTrueFWHMMesonMixedCaloConvPhoton->SetBinError(iPt,fMesonTrueFWHMErrorMixedCaloConvPhoton[iPt-1]);
                fHistoTrueMassMesonMixedCaloConvPhoton->SetBinContent(iPt,fMesonTrueMassMixedCaloConvPhoton[iPt-1]);
                fHistoTrueMassMesonMixedCaloConvPhoton->SetBinError(iPt,fMesonTrueMassErrorMixedCaloConvPhoton[iPt-1]);
            }
        }

        if (fIsMC) {
            fHistoTrueLambdaTail->SetBinContent(iPt,fMesonLambdaTailMCpar[iPt-1]);
            fHistoTrueLambdaTail->SetBinError(iPt,fMesonLambdaTailMCparError[iPt-1]);
        }
        fHistoLambdaTail->SetBinContent(iPt,fMesonLambdaTailpar[iPt-1]);
        fHistoLambdaTail->SetBinError(iPt,fMesonLambdaTailparError[iPt-1]);

        fHistoAmplitude->SetBinContent(iPt,fMesonAmplitudepar[iPt-1]);
        fHistoAmplitude->SetBinError(iPt,fMesonAmplitudeparError[iPt-1]);
        fHistoSigma->SetBinContent(iPt,fMesonSigmapar[iPt-1]);
        fHistoSigma->SetBinError(iPt,fMesonSigmaparError[iPt-1]);
        if (fIsMC) {
            fHistoTrueSigma->SetBinContent(iPt,fMesonTrueSigmapar[iPt-1]);
            fHistoTrueSigma->SetBinError(iPt,fMesonTrueSigmaparError[iPt-1]);
        }
        fHistoResidualBGlin->SetBinContent(iPt,fMesonResidualBGlin[iPt-1]);
        fHistoResidualBGlin->SetBinError(iPt,fMesonResidualBGlinError[iPt-1]);
        fHistoResidualBGcon->SetBinContent(iPt,fMesonResidualBGcon[iPt-1]);
        fHistoResidualBGcon->SetBinError(iPt,fMesonResidualBGconError[iPt-1]);
        if (fTotalBckYields[0][iPt-1] != 0){
            Double_t ratio      = fMesonYieldsResidualBckFunc[0][iPt-1]/fTotalBckYields[0][iPt-1];
            fHistoRatioResBGYield->SetBinContent(iPt,ratio);

            Double_t relErrorA  = fMesonYieldsResidualBckFuncError[0][iPt-1]/fMesonYieldsResidualBckFunc[0][iPt-1];
            Double_t relErrorB  = fTotalBckYieldsError[0][iPt-1]/fTotalBckYields[0][iPt-1];
            Double_t error      = ratio * TMath::Sqrt(relErrorA*relErrorA+relErrorB*relErrorB);
            fHistoRatioResBGYield->SetBinError(iPt,error);
        }
        if ((fMesonYieldsResidualBckFunc[0][iPt-1] + fMesonYieldsCorResidualBckFunc[0][iPt-1]) > 0){
            Double_t ratio      = fMesonYieldsResidualBckFunc[0][iPt-1]/(fMesonYieldsResidualBckFunc[0][iPt-1] + fMesonYieldsCorResidualBckFunc[0][iPt-1]);
            fHistoRatioResBGYieldToSPlusResBG->SetBinContent(iPt,ratio);

            Double_t relErrorA  = fMesonYieldsResidualBckFuncError[0][iPt-1]/fMesonYieldsResidualBckFunc[0][iPt-1];
            Double_t relErrorB  = TMath::Sqrt(fMesonYieldsCorResidualBckFuncError[0][iPt-1]*fMesonYieldsCorResidualBckFuncError[0][iPt-1]
                                               + fMesonYieldsResidualBckFuncError[0][iPt-1]*fMesonYieldsResidualBckFuncError[0][iPt-1]) /
                                    (fMesonYieldsResidualBckFunc[0][iPt-1]+fMesonYieldsCorResidualBckFunc[0][iPt-1]);
            Double_t error      = ratio * TMath::Sqrt(relErrorA*relErrorA+relErrorB*relErrorB);
            fHistoRatioResBGYieldToSPlusResBG->SetBinError(iPt,error);
        }
        for (Int_t m = 0; m < 4; m++){
            fHistoChi2[m]->SetBinContent(iPt,fMesonChi2[m][iPt-1]);
            fHistoChi2[m]->SetBinError(iPt,0);
        }

        fHistoResBGYield[0]->SetBinContent(iPt,fMesonYieldsResidualBckFunc[0][iPt-1]);
        fHistoResBGYield[0]->SetBinError(iPt,fMesonYieldsResidualBckFuncError[0][iPt-1]);
        for (Int_t m = 0; m < 3; m++){
            fHistoResBGYield[m+1]->SetBinContent(iPt,fMesonYieldsResBckOtherFunc[m][iPt-1]);
            fHistoResBGYield[m+1]->SetBinError(iPt,fMesonYieldsResBckOtherFuncError[m][iPt-1]);
        }


        // filling histogram arrays for normal, wide, narrow
        for (Int_t k = 0; k < 3; k++){
            fHistoMassWindowHigh[k]->SetBinContent(iPt,fMassWindowHigh[k][iPt-1]);
            fHistoMassWindowLow[k]->SetBinContent(iPt,fMassWindowLow[k][iPt-1]);
            fHistoSBdefaultMeson[k]->SetBinContent(iPt,fMesonSBdefault[k][iPt-1]);
            fHistoSBdefaultMeson[k]->SetBinError(iPt,fMesonSBdefaultError[k][iPt-1]);
            fHistoSigndefaultMeson[k]->SetBinContent(iPt,fMesonSigndefault[k][iPt-1]);
            fHistoSigndefaultMeson[k]->SetBinError(iPt,fMesonSigndefaultError[k][iPt-1]);
        }

        // filling histogram arrays for normal, wide, narrow, left, left wide, left narrow
        for (Int_t k = 0; k < 6; k++){
            fHistoYieldMeson[k]->SetBinContent(iPt,fMesonYieldsCorResidualBckFunc[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldMeson[k]->SetBinError(iPt,fMesonYieldsCorResidualBckFuncError[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldMesonPerEvent[k]->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFunc[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldMesonPerEvent[k]->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncError[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        }

        //normal
        fHistoYieldDiffBck[0]->SetBinContent(iPt,(fMesonYieldsCorResidualBckFunc[0][iPt-1]-(fMesonYieldsResidualBckFunc[0][iPt-1]-fMesonYieldsResBckOtherFunc[0][iPt-1]))/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldDiffBck[1]->SetBinContent(iPt,(fMesonYieldsCorResidualBckFunc[0][iPt-1]+(fMesonYieldsResidualBckFunc[0][iPt-1]-fMesonYieldsResBckOtherFunc[0][iPt-1]))/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        //normal exp
        fHistoYieldDiffBck[2]->SetBinContent(iPt,(fMesonYieldsCorResidualBckFunc[0][iPt-1]-(fMesonYieldsResidualBckFunc[0][iPt-1]-fMesonYieldsResBckOtherFunc[1][iPt-1]))/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldDiffBck[3]->SetBinContent(iPt,(fMesonYieldsCorResidualBckFunc[0][iPt-1]+(fMesonYieldsResidualBckFunc[0][iPt-1]-fMesonYieldsResBckOtherFunc[1][iPt-1]))/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        //normal exp2
        fHistoYieldDiffBck[4]->SetBinContent(iPt,(fMesonYieldsCorResidualBckFunc[0][iPt-1]-(fMesonYieldsResidualBckFunc[0][iPt-1]-fMesonYieldsResBckOtherFunc[2][iPt-1]))/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldDiffBck[5]->SetBinContent(iPt,(fMesonYieldsCorResidualBckFunc[0][iPt-1]+(fMesonYieldsResidualBckFunc[0][iPt-1]-fMesonYieldsResBckOtherFunc[2][iPt-1]))/(fBinsPt[iPt]-fBinsPt[iPt-1]));


        if (fIsMC) {
            for (Int_t k = 0; k < 3; k++){
                fHistoYieldTrueMeson[k]->SetBinContent(iPt,fMesonTrueYields[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoYieldTrueMeson[k]->SetBinError(iPt,fMesonTrueYieldsError[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoYieldTrueMesonReweighted[k]->SetBinContent(iPt,fMesonTrueYieldsReweighted[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoYieldTrueMesonReweighted[k]->SetBinError(iPt,fMesonTrueYieldsReweightedError[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoYieldTrueMesonUnweighted[k]->SetBinContent(iPt,fMesonTrueYieldsUnweighted[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoYieldTrueMesonUnweighted[k]->SetBinError(iPt,fMesonTrueYieldsUnweightedError[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                if (fEnableDCMeson && k == 0){
                    fHistoYieldTrueMesonDC->SetBinContent(iPt,fMesonTrueYieldsDC[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                    fHistoYieldTrueMesonDC->SetBinError(iPt,fMesonTrueYieldsDCError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                }
            }
        }

        // Histos for integration at the left of the peak
        fHistoMassMesonLeft->SetBinContent(iPt,fMesonMassLeft[iPt-1]);
        fHistoMassMesonLeft->SetBinError(iPt,fMesonMassLeftError[iPt-1]);
        fHistoFWHMMesonLeft->SetBinContent(iPt,fMesonFWHMLeft[iPt-1]);
        fHistoFWHMMesonLeft->SetBinError(iPt,fMesonFWHMLeftError[iPt-1]);

        // fill secondary histograms
        if (fIsMC){
            for (Int_t k = 0; k< 3; k++){
                for (Int_t j = 0; j< 4; j++){
                    fHistoYieldTrueSecMeson[k][j]->SetBinContent(iPt,fMesonTrueSecYields[k][j][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                    fHistoYieldTrueSecMeson[k][j]->SetBinError(iPt,fMesonTrueSecYieldsError[k][j][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                }
            }
        }
    }
}


//****************************************************************************
//******** Fit of Signal+ BG with Gaussian + Exponential + Linear BG *********
//****************************************************************************
void FitSubtractedInvMassInPtBins(TH1D* histoMappingSignalInvMassPtBinSingle, Double_t* mesonIntDeltaRangeFit, Int_t ptBin, Bool_t vary){

    //--------------------------------------------------------------------------------------
    // determine fit start values for amplitude and special settings for different collision
    // systems and energies
    //--------------------------------------------------------------------------------------
    histoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude     = histoMappingSignalInvMassPtBinSingle->GetMaximum();
    if( fMode == 4 || fMode == 5 ){
        mesonAmplitude = 0;
        for(Int_t i=histoMappingSignalInvMassPtBinSingle->FindBin(fMesonFitRange[0]); i<histoMappingSignalInvMassPtBinSingle->FindBin(fMesonFitRange[1]) ; i++){
            if(histoMappingSignalInvMassPtBinSingle->GetBinContent(i)>mesonAmplitude) mesonAmplitude = histoMappingSignalInvMassPtBinSingle->GetBinContent(i);
        }
    }
    Double_t mesonAmplitudeMin  = mesonAmplitude*10./100.;
    Double_t mesonAmplitudeMax  = mesonAmplitude*400./100.;

    // Settings specific for PbPb collisions
    if (fEnergyFlag.Contains("PbPb")){
        if (fPrefix.CompareTo("Pi0") == 0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
            //pi0 meson amplitude settings
            if(fMode == 0){
                if(((TString)histoMappingSignalInvMassPtBinSingle->GetName()).Contains("Left") && (GetCentralityString(fEventCutSelection)).CompareTo("20-50%")==0){
                    if(ptBin == 3) mesonAmplitudeMin = mesonAmplitude*80./100.;
                    else  mesonAmplitudeMin = mesonAmplitude*98./100.;
                } else {
                    if(ptBin > 1) mesonAmplitudeMin = mesonAmplitude*98./100.;
                    else  mesonAmplitudeMin = mesonAmplitude*80./100.;
                }
                mesonAmplitudeMax = mesonAmplitude*115./100.;
            }
            if (fMode == 2 || fMode == 13 || fMode == 3) {
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*600./100.;
            }
            if (fMode == 4 || fMode == 12 || fMode == 5) {
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                mesonAmplitudeMax = mesonAmplitude*800./100.;
            }

        } else {
            //eta meson amplitude settings
            if(fMode == 0){
                if(!(GetCentralityString(fEventCutSelection)).BeginsWith("2")){
                    if(ptBin>2) mesonAmplitudeMin = mesonAmplitude*60./100.;
                    else mesonAmplitudeMin = mesonAmplitude*80./100.;
                } else {
                    if(ptBin>2) mesonAmplitudeMin = mesonAmplitude*80./100.;
                    else mesonAmplitudeMin = mesonAmplitude*40./100.;
                }
                mesonAmplitudeMax = mesonAmplitude*115./100.;
            }
            if (fMode == 2 || fMode == 13 || fMode == 3){
                mesonAmplitudeMin = mesonAmplitude*10./100.;
            }
            if (fMode == 4 || fMode == 12 || fMode == 5){
                mesonAmplitudeMin = mesonAmplitude*5./100.;
            }
        }
    // Settings specific for pPb collisions
    } else if (fEnergyFlag.Contains("pPb") ){
        // pi0 reconstruction
        if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
            if ( fMode == 0) {                                      // PCM
                mesonAmplitudeMin = mesonAmplitude*92./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
            } else if ( fMode == 2 || fMode == 13) {  // PCM-EMC, PCM-DMC
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*600./100.;
                if(fBinsPt[ptBin] >= 3.0 ) {
                    if (fIsMC == 0){
                        fMesonLambdaTail            = 0.011;
                        fMesonLambdaTailRange[0]    = 0.011;
                        fMesonLambdaTailRange[1]    = 0.011;
                    } else {
                        fMesonLambdaTail            = 0.0095;
                        fMesonLambdaTailRange[0]    = 0.0095;
                        fMesonLambdaTailRange[1]    = 0.0095;
                    }
                } else if(fBinsPt[ptBin] >= 20.0 ) {
                    fMesonLambdaTail                = 0.0;
                    fMesonLambdaTailRange[0]        = 0.0;
                    fMesonLambdaTailRange[1]        = 0.0;
                }
            } else if (fMode == 3) {  // PCM-PHOS
                mesonAmplitudeMin = mesonAmplitude*90./100.;
                mesonAmplitudeMax = mesonAmplitude*120./100.;
                if(fBinsPt[ptBin] >= 4.0 ) {
                    if (fIsMC == 0){
                        fMesonLambdaTail            = 0.007;
                        fMesonLambdaTailRange[0]    = 0.007;
                        fMesonLambdaTailRange[1]    = 0.007;
                    } else {
                        fMesonLambdaTail            = 0.007;
                        fMesonLambdaTailRange[0]    = 0.007;
                        fMesonLambdaTailRange[1]    = 0.007;
                    }
                }
            } else if (fMode == 4 || fMode == 12) {   // EMC, DMC
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                mesonAmplitudeMax = mesonAmplitude*1000./100.;
                if(fBinsPt[ptBin] >= 10) {
                    fMesonLambdaTail            = 0.015;
                    fMesonLambdaTailRange[0]    = 0.015;
                    fMesonLambdaTailRange[1]    = 0.015;
                }
            } else if (fMode == 5) {                                // PHOS
                if (ptBin < 21) mesonAmplitudeMin = 1./100.;
                else mesonAmplitudeMin = mesonAmplitude*10./100.;
                mesonAmplitudeMax = mesonAmplitude*1000./100.;
                if(fBinsPt[ptBin] >= 10) {
                    fMesonLambdaTail            = 0.005;
                    fMesonLambdaTailRange[0]    = 0.005;
                    fMesonLambdaTailRange[1]    = 0.005;
                }
            } else {                                                // defaults
                mesonAmplitudeMin = mesonAmplitude*92./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
            }
        // eta reconstruction
        } else {
            if ( fMode == 0) {                                      // PCM
                mesonAmplitudeMin = mesonAmplitude*50./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
            } else if ( fMode == 2 || fMode == 13 ) {               // PCM-EMC, PCM-DMC
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                mesonAmplitudeMax = mesonAmplitude*600./100.;
            } else if ( fMode == 3) {                               // PCM-PHOS
                mesonAmplitudeMin = mesonAmplitude*0.1/100.;
                if (fBinsPt[ptBin] < 2.)
                    mesonAmplitudeMax = mesonAmplitude*105./100.;
                else
                    mesonAmplitudeMax = mesonAmplitude*110./100.;
            } else if (fMode == 4 || fMode == 12 ) {                //EMC, DMC
                mesonAmplitudeMin = mesonAmplitude*0.1/100.;
                mesonAmplitudeMax = mesonAmplitude*300./100.;
                if(fBinsPt[ptBin] >= 4) {
                    fMesonLambdaTail            = 0.033;
                    fMesonLambdaTailRange[0]    = 0.033;
                    fMesonLambdaTailRange[1]    = 0.033;
                }
            } else if (fMode == 5) {                                // PHOS
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                mesonAmplitudeMax = mesonAmplitude*1000./100.;
                fMesonLambdaTail            = 0.01;
                fMesonLambdaTailRange[0]    = 0.01;
                fMesonLambdaTailRange[1]    = 0.01;
            } else {                                                // defaults
                mesonAmplitudeMin = mesonAmplitude*50./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
            }
        }
    // settings specific for pp collisions
    } else {
        // pi0 reconstruction
        if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
            if (fMode == 0){                                        // PCM
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
                if ( !fEnergyFlag.CompareTo("8TeV")){
                    if ((ptBin > 2)&&ptBin<100 ){
                        fMesonWidthRange[0]         = 0.001;
                        fMesonWidthRange[1]         = 0.009;
                    }
                } else  if ( !fEnergyFlag.CompareTo("7TeV")){
                    if ((ptBin > 2)&&ptBin<100 ){
                        fMesonWidthRange[0]         = 0.001;
                        fMesonWidthRange[1]         = 0.009;
                    }
                } else  if ( !fEnergyFlag.CompareTo("5TeV2017")){
                    if ((ptBin < 3)){
                        mesonAmplitudeMin = mesonAmplitude*60./100.;
                        mesonAmplitudeMax = mesonAmplitude*105./100.;
                    }
                    fMesonWidthRange[0]         = 0.001;
                    fMesonWidthRange[1]         = 0.009;
                }
            } else if (fMode == 2 || fMode == 13 ) {                // PCM-EMC, PCM-DMC
                if( fEnergyFlag.CompareTo("8TeV") == 0 ){
                  mesonAmplitudeMin = mesonAmplitude*90./100.;
                  mesonAmplitudeMax = mesonAmplitude*600./100.;
                } else {
                    mesonAmplitudeMin = mesonAmplitude*98./100.;
                    mesonAmplitudeMax = mesonAmplitude*600./100.;
                }
            } else if (fMode == 3) {                                // PCM-PHOS
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*150./100.;
            } else if (fMode == 4 || fMode == 12 || fMode == 5) {   // EMC, DMC, PHOS
                if( fEnergyFlag.CompareTo("8TeV") == 0 ){
                    mesonAmplitudeMin = mesonAmplitude*90./100.;
                    mesonAmplitudeMax = mesonAmplitude*400./100.;
                    TString trigger = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
                    if(trigger.CompareTo("52") == 0 && fBinsPt[ptBin]==13. && fPrefix.CompareTo("Pi0EtaBinning")){
                        fMesonFitRange[0] = 0.08;
                        fMesonFitRange[1] = 0.25;
                    }else if(trigger.CompareTo("52") == 0 && fBinsPt[ptBin]>=14. && fPrefix.CompareTo("Pi0EtaBinning")){
                        fMesonFitRange[0] = 0.08;
                        fMesonFitRange[1] = 0.29;
                    }
                } else {
                    mesonAmplitudeMin = mesonAmplitude*10./100.;
                    mesonAmplitudeMax = mesonAmplitude*400./100.;
                }
            } else {                                                // defaults
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
            }
        // eta reconstruction
        } else {
            if (fMode == 0){                                        // PCM
                if ( !fEnergyFlag.CompareTo("7TeV") ){
                    mesonAmplitudeMin = mesonAmplitude*80./100.;
                    mesonAmplitudeMax = mesonAmplitude*115./100.;
                    if(ptBin < 3) mesonAmplitudeMax = mesonAmplitude*100/100.;
                } else if ( !fEnergyFlag.CompareTo("8TeV") ){
                    mesonAmplitudeMin = mesonAmplitude*65./100.;
                    if(ptBin > 2)  mesonAmplitudeMin = mesonAmplitude*85./100.;
                    mesonAmplitudeMax = mesonAmplitude*115./100.;
                    if(ptBin < 3) mesonAmplitudeMax = mesonAmplitude*100/100.;
                } else if ( !fEnergyFlag.CompareTo("5TeV2017") ){
                    mesonAmplitudeMin = mesonAmplitude*80./100.;
                    mesonAmplitudeMax = mesonAmplitude*120./100.;
                } else {
                    mesonAmplitudeMin = mesonAmplitude*50./100.;
                    mesonAmplitudeMax = mesonAmplitude*115./100.;
                }
            } else if (fMode == 2 || fMode == 13 ){                 // PCM-EMC, PCM-DMC
                if( fEnergyFlag.CompareTo("8TeV") == 0 ){
                    mesonAmplitudeMin = mesonAmplitude*20./100.;
                    mesonAmplitudeMax = mesonAmplitude*115./100.;
                } else {
                    mesonAmplitudeMin = mesonAmplitude*10./100.;
                    mesonAmplitudeMax = mesonAmplitude*115./100.;
                }
            } else if (fMode == 3 ){                                // PCM-PHOS
                mesonAmplitudeMin = mesonAmplitude*0.1/100.;
                if (fBinsPt[ptBin] < 2.)
                    mesonAmplitudeMax = mesonAmplitude*150./100.;
                else
                    mesonAmplitudeMax = mesonAmplitude*400./100.;
            } else if (fMode == 4 || fMode == 12 || fMode == 5){    // EMC, DMC, PHOS
                mesonAmplitudeMin = mesonAmplitude*5./100.;
                mesonAmplitudeMax = mesonAmplitude*300./100.;
            } else {                                                // defaults
                mesonAmplitudeMin = mesonAmplitude*50./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
            }
        }
    }

    //--------------------------------------------------------------------------------------
    // define complete fitting functions
    //--------------------------------------------------------------------------------------
    fFitReco= NULL;
    fFitReco = new TF1("GaussExpLinear","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",fMesonFitRange[0],fMesonFitRange[1]);

    fFitGausExp =NULL;
    fFitGausExp = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

    fFitLinearBck = NULL;
    fFitLinearBck = new TF1("Linear","[0]+[1]*x",fMesonFitRange[0],fMesonFitRange[1]);
    //--------------------------------------------------------------------------------------

    //--------------------------------------------------------------------------------------
    //------------------------- Set parameter boundaries -----------------------------------
    //--------------------------------------------------------------------------------------
    // set amplitude start value
    fFitReco->SetParameter(0,mesonAmplitude);
    // set ranges for amplitude fitting
    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);

    // set mass start value
    fFitReco->SetParameter(1,fMesonMassExpect);
    if (fEnergyFlag.CompareTo("13TeVLowB") == 0 && fMode == 0 && ptBin < 2 )
        fFitReco->SetParameter(1,fMesonMassExpect*0.7);
    // set ranges for mass fitting
    fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.15);
    if( (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0)){
        if(fEnergyFlag.Contains("pPb") ){
            if(fMode == 4 ) {
                fFitReco->SetParLimits(1,fMesonMassExpect*0.5,fMesonMassExpect*2);
                fFitReco->SetParLimits(2,fMesonWidthRange[0]*0.5,fMesonWidthRange[1]*2.);
                if(ptBin >= 31) {fFitReco->FixParameter(3,0.015); fFitReco->SetParameter(1,0.156); }
            }else if(fMode == 3) {
                if (fBinsPt[ptBin] > 7. ) fFitReco->SetParLimits(1,fMesonMassExpect*0.995,fMesonMassExpect*1.15);
                if (fBinsPt[ptBin] > 9. ) fFitReco->SetParLimits(1,fMesonMassExpect*0.8,fMesonMassExpect*1.2);
            }
        } else if (fEnergyFlag.CompareTo("8TeV") == 0 ){
            if (fMode == 4 || fMode == 12 )
                fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.3);
        } else if (fEnergyFlag.CompareTo("13TeVLowB") == 0 && fMode == 0 && ptBin < 2 ){
                 fFitReco->SetParLimits(1,fMesonMassExpect*0.5,fMesonMassExpect*1);
        } else if ( fEnergyFlag.Contains("PbPb") || fEnergyFlag.Contains("XeXe") ){
            if (fMode == 4 || fMode == 12 )
                fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.5);
        }
    } else {
        if( fEnergyFlag.Contains("pPb") ){
            if ( fMode == 2 )
                fFitReco->SetParLimits(1,fMesonMassExpect*0.8,fMesonMassExpect*1.3);
            else if  ( fMode == 3 )
                fFitReco->SetParLimits(1,fMesonMassExpect*0.985,fMesonMassExpect*1.15);
            else if  ( fMode == 4 || fMode == 12 )
                fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.3);
        } else if (fEnergyFlag.CompareTo("8TeV") == 0 || fEnergyFlag.Contains("PbPb") || fEnergyFlag.Contains("XeXe") ) {
            if ( fMode == 4 || fMode == 12 )
                fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.3);
        }
    }

    // set width start value
    fFitReco->SetParameter(2,fMesonWidthExpect);
    // set ranges for width fitting
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);

    // set start value and ranges for Lambda fitting
    if (fMesonLambdaTail == fMesonLambdaTailRange[0] && fMesonLambdaTail == fMesonLambdaTailRange[1] ){
        fFitReco->FixParameter(3,fMesonLambdaTail);
    } else {
        fFitReco->SetParameter(3,fMesonLambdaTail);
        fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);
    }
    //--------------------------------------------------------------------------------------
    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);

    //--------------------------------------------------------------------------------------
    //------------------------- start fitting & set style fit ------------------------------
    //--------------------------------------------------------------------------------------
    histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");


    fFitReco->SetLineColor(3);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);
    fFitReco->SetNpx(10000);

    //--------------------------------------------------------------------------------------
    // change settign for next pt bin if selected
    //--------------------------------------------------------------------------------------
    if (vary && !fIsMC && (fMode == 0 || fMode == 9 )){
        if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 && fPrefix.CompareTo("Pi0") ==0 && ptBin >=17){
            cout << "Skipping the vary option for this case, pt: " << ptBin << endl;
        } else if (fEnergyFlag.Contains("pPb") && (ptBin >= 20) ){
            cout << "Skipping the vary option for this case" << endl;
        } else if (fEnergyFlag.Contains("13TeV")&&(fMode == 0)){
            cout << "Skipping the vary option for this case" << endl;
        } else {// ...do what you are supposed to....
            if (!(fMesonLambdaTail == fMesonLambdaTailRange[0] && fMesonLambdaTail == fMesonLambdaTailRange[1]) ){
                fMesonLambdaTail = fFitReco->GetParameter(3);
                if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 && fPrefix.CompareTo("Eta") ==0){
                    fMesonLambdaTailRange[0] = 0.7*fFitReco->GetParameter(3);
                    fMesonLambdaTailRange[1] = 1.5*fFitReco->GetParameter(3);
                } else {
                    fMesonLambdaTailRange[0] = 0.9*fFitReco->GetParameter(3);
                    fMesonLambdaTailRange[1] = 1.1*fFitReco->GetParameter(3);
                }
            }
            fMesonWidthExpect = fFitReco->GetParameter(2);
            fMesonWidthRange[0] = 0.5*fFitReco->GetParameter(2);
            fMesonWidthRange[1] = 1.5*fFitReco->GetParameter(2);
        }
    }

    //--------------------------------------------------------------------------------------
    //------------------- set parameters of subfunctions -----------------------------------
    //--------------------------------------------------------------------------------------
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
    TVirtualFitter* fitter  = TVirtualFitter::GetFitter();

    //--------------------------------------------------------------------------------------
    //--------------- calculate integrals with proper errors if fit converged---------------
    //--------------------------------------------------------------------------------------
    fIntLinearBck           = 0;
    fIntLinearBckError      = 0;
    if(TString(gMinuit->fCstatu.Data()).Contains("CONVERGED") == 1 || TString(gMinuit->fCstatu.Data()).Contains("SUCCESSFUL") == 1 || TString(gMinuit->fCstatu.Data()).Contains("PROBLEMS") == 1){
        binCenterStart          = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+mesonIntDeltaRangeFit[0]);
        startBinEdge            = histoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
        binCenterEnd            = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+mesonIntDeltaRangeFit[1]);
        endBinEdge              = histoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

        Int_t nFreePar          = fFitReco->GetNumberFreeParameters();
        double *covMatrix       = fitter->GetCovarianceMatrix();
        Float_t intLinearBack   = fFitLinearBck->GetParameter(0)*(endBinEdge-startBinEdge)+
            0.5*fFitLinearBck->GetParameter(1)*(endBinEdge*endBinEdge-startBinEdge*startBinEdge);
        Float_t errorLinearBck  = pow((pow( (endBinEdge-startBinEdge)*fFitReco->GetParError(4),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fFitReco->GetParError(5),2)
                                            +2*covMatrix[nFreePar*nFreePar-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5);

        fFileDataLog << "Parameter for bin " << ptBin << endl;
        fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
        fFileDataLog << "Linear: \t"<<fFitReco->GetParameter(4)<<"+-" << fFitReco->GetParError(4) << "\t "<<fFitReco->GetParameter(5) <<"+-" << fFitReco->GetParError(5)<< endl;

        fIntLinearBck           = intLinearBack/histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
        fIntLinearBckError      = errorLinearBck/histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
    } else {
        fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }
    fFitReco->DrawCopy("same");
}

//****************************************************************************
//******** Fit of Signal+ BG with Gaussian + Exponential + Pol2 BG *********
//****************************************************************************
void FitSubtractedPol2InvMassInPtBins(TH1D* histoMappingSignalInvMassPtBinSingle, Double_t* mesonIntDeltaRangeFit, Int_t ptBin, Bool_t vary){

    //    cout<<"Start Fitting spectra"<<endl;
    histoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude         = histoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin      = 0;
    Double_t mesonAmplitudeMax      = 0;

    if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0){
        if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
            if(ptBin == 1) mesonAmplitudeMin = mesonAmplitude*70./100.;
            else if(ptBin > 1 && ptBin < 4) mesonAmplitudeMin = mesonAmplitude*90./100.;
            else if(ptBin > 17) mesonAmplitudeMin = mesonAmplitude*80./100.;
            else mesonAmplitudeMin = mesonAmplitude*98./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            if (fMode == 2 || fMode == 13 || fMode == 3) {
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*600./100.;
            }
            if (fMode == 4 || fMode == 12 || fMode == 5) {
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                mesonAmplitudeMax = mesonAmplitude*400./100.;
            }
        } else {
            mesonAmplitudeMin = mesonAmplitude*30./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            if (fMode == 2 || fMode == 13 || fMode == 3){
                mesonAmplitudeMin = mesonAmplitude*10./100.;
            }
            if (fMode == 4 || fMode == 12 || fMode == 5){
                mesonAmplitudeMin = mesonAmplitude*5./100.;
            }
        }
    } else {
        if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
            mesonAmplitudeMin = mesonAmplitude*98./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            if (fEnergyFlag.Contains("pPb_5.023TeV") ) mesonAmplitudeMin = mesonAmplitude*92./100.;
            if ( (fMode == 0) && !fEnergyFlag.CompareTo("8TeV")){
                if ((ptBin > 2)&&ptBin<100 ){
                    fMesonWidthRange[0]         = 0.001;
                    fMesonWidthRange[1]         = 0.009;
                }
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
            }
            if ((fMode == 0) && !fEnergyFlag.CompareTo("7TeV")){
                if ((ptBin > 2)&&ptBin<100 ){
                    fMesonWidthRange[0]         = 0.001;
                    fMesonWidthRange[1]         = 0.009;
                }
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
            }
            if (fMode == 2 || fMode == 13 || fMode == 3 ) {
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*600./100.;
                if( fEnergyFlag.CompareTo("8TeV") == 0 ){
                  mesonAmplitudeMin = mesonAmplitude*90./100.;
                  mesonAmplitudeMax = mesonAmplitude*600./100.;
                }
            }
            if (fMode == 4 || fMode == 12 || fMode == 5) {
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                mesonAmplitudeMax = mesonAmplitude*400./100.;
                if( fEnergyFlag.CompareTo("8TeV") == 0 ){
                    mesonAmplitudeMin = mesonAmplitude*90./100.;
                    mesonAmplitudeMax = mesonAmplitude*400./100.;
                    TString trigger = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
                    if(trigger.CompareTo("52") == 0 && fBinsPt[ptBin]>=13. && fPrefix.CompareTo("Pi0EtaBinning")){
                        fMesonFitRange[0] = 0.08;
                        fMesonFitRange[1] = 0.29;
                    }
                } else if (fEnergyFlag.Contains("pPb_5.023TeVRun2") ){
                    if (ptBin < 8)  mesonAmplitudeMin = mesonAmplitude*1./100.;
                    else            mesonAmplitudeMin = mesonAmplitude*10./100.;
                }
            }

        } else {
            mesonAmplitudeMin = mesonAmplitude*50./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            if (fMode == 0 && !fEnergyFlag.CompareTo("8TeV") ){
                mesonAmplitudeMin = mesonAmplitude*65./100.;
                if(ptBin > 2)  mesonAmplitudeMin = mesonAmplitude*85./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
                if(ptBin < 3) mesonAmplitudeMax = mesonAmplitude*100/100.;
            }
            if (fMode == 2 || fMode == 13 || fMode == 3 ){
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                if( fEnergyFlag.CompareTo("8TeV") == 0 ){
                  mesonAmplitudeMin = mesonAmplitude*20./100.;
                  mesonAmplitudeMax = mesonAmplitude*115./100.;
                }
            }
            if (fMode == 4 || fMode == 12 || fMode == 5 ){
                mesonAmplitudeMin = mesonAmplitude*5./100.;
            }
        }
    }

    fFitReco        = NULL;
    fFitReco        = new TF1("GaussExpPol2","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x+[6]*x*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x+[6]*x*x)",fMesonFitRange[0],fMesonFitRange[1]);

    fFitGausExp     = NULL;
    fFitGausExp     = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

    fFitLinearBck   = NULL;
    fFitLinearBck   = new TF1("BGfitPol2","[0]+[1]*x+[2]*x*x",fMesonFitRange[0],fMesonFitRange[1]);

    fFitReco->SetParameter(0,mesonAmplitude);
    fFitReco->SetParameter(1,fMesonMassExpect);
    fFitReco->SetParameter(2,fMesonWidthExpect);
    if (fMesonLambdaTail == fMesonLambdaTailRange[0] && fMesonLambdaTail == fMesonLambdaTailRange[1] ){
        fFitReco->FixParameter(3,fMesonLambdaTail);
    } else {
        fFitReco->SetParameter(3,fMesonLambdaTail);
        fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);
    }
    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.15);
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);

    histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    fFitReco->SetLineColor(3);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);
    fFitReco->SetNpx(10000);

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
    fFitLinearBck->SetParameter(2,fFitReco->GetParameter(6));
    fFitLinearBck->SetParError(0,fFitReco->GetParError(4));
    fFitLinearBck->SetParError(1,fFitReco->GetParError(5));
    fFitLinearBck->SetParError(2,fFitReco->GetParError(6));
    Int_t binCenterStart;
    Double_t startBinEdge;
    Int_t binCenterEnd;
    Double_t endBinEdge;

    TVirtualFitter * fitter = TVirtualFitter::GetFitter();

    fIntLinearBck = 0;
    fIntLinearBckError = 0;
    if(TString(gMinuit->fCstatu.Data()).Contains("CONVERGED") == 1 || TString(gMinuit->fCstatu.Data()).Contains("SUCCESSFUL") == 1 || TString(gMinuit->fCstatu.Data()).Contains("PROBLEMS") == 1){
        binCenterStart  = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+mesonIntDeltaRangeFit[0]);
        startBinEdge    = histoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
        binCenterEnd    = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+mesonIntDeltaRangeFit[1]);
        endBinEdge      = histoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

        Int_t nFreePar      = fFitReco->GetNumberFreeParameters();
        double * covMatrix  = fitter->GetCovarianceMatrix();

        Float_t intBack     = fFitLinearBck->Integral(startBinEdge, endBinEdge);
        Float_t errorBck    = fFitLinearBck->IntegralError(startBinEdge, endBinEdge);

        fFileDataLog << "Parameter for bin " << ptBin << endl;
        fFileDataLog << "Pol2: \t"<<fFitReco->GetParameter(4)<<"+-" << fFitReco->GetParError(4) << "\t "<<fFitReco->GetParameter(5) <<"+-" << fFitReco->GetParError(5)<< "\t "<<fFitReco->GetParameter(6) <<"+-" << fFitReco->GetParError(6)<<endl;

        fIntLinearBck       = intBack/histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
        fIntLinearBckError  = errorBck/histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
    } else {
        fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }
}


//****************************************************************************
//******** Fit of Signal+ BG with Gaussian + Exponential + Exp1 *********
//****************************************************************************
void FitSubtractedExp1InvMassInPtBins(TH1D* histoMappingSignalInvMassPtBinSingle, Double_t* mesonIntDeltaRangeFit, Int_t ptBin, Bool_t vary){

    //    cout<<"Start Fitting spectra"<<endl;
    histoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude         = histoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin      = 0;
    Double_t mesonAmplitudeMax      = 0;

    if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0){
        if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
            if(ptBin == 1) mesonAmplitudeMin = mesonAmplitude*70./100.;
            else if(ptBin > 1 && ptBin < 4) mesonAmplitudeMin = mesonAmplitude*90./100.;
            else if(ptBin > 17) mesonAmplitudeMin = mesonAmplitude*80./100.;
            else mesonAmplitudeMin = mesonAmplitude*98./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            if (fMode == 2 || fMode == 13 || fMode == 3) {
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*600./100.;
            }
            if (fMode == 4 || fMode == 12 || fMode == 5) {
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                mesonAmplitudeMax = mesonAmplitude*400./100.;
            }
        } else {
            mesonAmplitudeMin = mesonAmplitude*30./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            if (fMode == 2 || fMode == 13 || fMode == 3){
                mesonAmplitudeMin = mesonAmplitude*10./100.;
            }
            if (fMode == 4 || fMode == 12 || fMode == 5){
                mesonAmplitudeMin = mesonAmplitude*5./100.;
            }
        }
    } else {
        if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
            mesonAmplitudeMin = mesonAmplitude*98./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            if (fEnergyFlag.Contains("pPb_5.023TeV") ) mesonAmplitudeMin = mesonAmplitude*92./100.;
            if (fMode == 0 && !fEnergyFlag.CompareTo("8TeV")){
                if ((ptBin > 2)&&ptBin<100 ){
                    fMesonWidthRange[0]         = 0.001;
                    fMesonWidthRange[1]         = 0.009;
                }
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
            }
            if (fMode == 0 && !fEnergyFlag.CompareTo("7TeV")){
                if ((ptBin > 2)&&ptBin<100 ){
                    fMesonWidthRange[0]         = 0.001;
                    fMesonWidthRange[1]         = 0.009;
                }
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
            }
            if (fMode == 2 || fMode == 13 || fMode == 3) {
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*600./100.;
                if( fEnergyFlag.CompareTo("8TeV") == 0 ){
                  mesonAmplitudeMin = mesonAmplitude*90./100.;
                  mesonAmplitudeMax = mesonAmplitude*600./100.;
                }
            }
            if (fMode == 4 || fMode == 12 || fMode == 5) {
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                mesonAmplitudeMax = mesonAmplitude*400./100.;
                if( fEnergyFlag.CompareTo("8TeV") == 0 ){
                  mesonAmplitudeMin = mesonAmplitude*90./100.;
                  mesonAmplitudeMax = mesonAmplitude*400./100.;
                    TString trigger = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
                    if(trigger.CompareTo("52") == 0 && fBinsPt[ptBin]>=13. && fPrefix.CompareTo("Pi0EtaBinning")){
                        fMesonFitRange[0] = 0.08;
                        fMesonFitRange[1] = 0.29;
                    }
                }
            }

        } else {
            mesonAmplitudeMin = mesonAmplitude*50./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            if (fMode == 0 && !fEnergyFlag.CompareTo("8TeV") ){
                mesonAmplitudeMin = mesonAmplitude*65./100.;
                if(ptBin > 2)  mesonAmplitudeMin = mesonAmplitude*85./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
                if(ptBin < 3) mesonAmplitudeMax = mesonAmplitude*100/100.;
            }
            if (fMode == 2 || fMode == 13 || fMode == 3){
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                if( fEnergyFlag.CompareTo("8TeV") == 0 ){
                  mesonAmplitudeMin = mesonAmplitude*20./100.;
                  mesonAmplitudeMax = mesonAmplitude*115./100.;
                }
            }
            if (fMode == 4 || fMode == 12 || fMode == 5){
                mesonAmplitudeMin = mesonAmplitude*5./100.;
            }
        }
    }

    fFitReco        = NULL;
    fFitReco        = new TF1("GaussExpPol2","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]*exp([5]*x))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]*exp([5]*x))",fMesonFitRange[0],fMesonFitRange[1]);

    fFitGausExp     = NULL;
    fFitGausExp     = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

    fFitLinearBck   = NULL;
    fFitLinearBck   = new TF1("BGfitExp1","[0]*exp([1]*x)",fMesonFitRange[0],fMesonFitRange[1]);


    fFitReco->SetParameter(0,mesonAmplitude);
    fFitReco->SetParameter(1,fMesonMassExpect);
    fFitReco->SetParameter(2,fMesonWidthExpect);
    if (fMesonLambdaTail == fMesonLambdaTailRange[0] && fMesonLambdaTail == fMesonLambdaTailRange[1] ){
        fFitReco->FixParameter(3,fMesonLambdaTail);
    } else {
        fFitReco->SetParameter(3,fMesonLambdaTail);
        fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);
    }
    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.15);
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);

    TVirtualFitter * fitter2= NULL;
    histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    fitter2 = TVirtualFitter::GetFitter();

    fFitReco->SetLineColor(3);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);
    fFitReco->SetNpx(10000);
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

    fIntLinearBck = 0;
    fIntLinearBckError = 0;
    if( TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ||
        TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
        binCenterStart  = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+mesonIntDeltaRangeFit[0]);
        startBinEdge    = histoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
        binCenterEnd    = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+mesonIntDeltaRangeFit[1]);
        endBinEdge      = histoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

        Int_t nFreePar      = fFitReco->GetNumberFreeParameters();
        double * covMatrix  = fitter2->GetCovarianceMatrix();

        Float_t intBack     = fFitLinearBck->Integral(startBinEdge, endBinEdge);
        Float_t errorBck    = fFitLinearBck->IntegralError(startBinEdge, endBinEdge);

        fFileDataLog << "Parameter for bin " << ptBin << endl;
        fFileDataLog << "Exp1: \t"<<fFitReco->GetParameter(4)<<"+-" << fFitReco->GetParError(4) << "\t "<<fFitReco->GetParameter(5) <<"+-" << fFitReco->GetParError(5)<<endl;

        fIntLinearBck       = intBack/histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
        fIntLinearBckError  = errorBck/histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
    } else {
        fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }
}

//****************************************************************************
//******** Fit of Signal+ BG with Gaussian + Exponential + Exp2 *********
//****************************************************************************
void FitSubtractedExp2InvMassInPtBins(TH1D* histoMappingSignalInvMassPtBinSingle, Double_t* mesonIntDeltaRangeFit, Int_t ptBin, Bool_t vary){

    //    cout<<"Start Fitting spectra"<<endl;
    histoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude         = histoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin      = 0;
    Double_t mesonAmplitudeMax      = 0;

    if (fEnergyFlag.Contains("PbPb")){
        if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
            if(ptBin == 1) mesonAmplitudeMin = mesonAmplitude*70./100.;
            else if(ptBin > 1 && ptBin < 4) mesonAmplitudeMin = mesonAmplitude*90./100.;
            else if(ptBin > 17) mesonAmplitudeMin = mesonAmplitude*80./100.;
            else mesonAmplitudeMin = mesonAmplitude*98./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            if (fMode == 2 || fMode == 13 || fMode == 3) {
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*600./100.;
            }
            if (fMode == 4 || fMode == 12 || fMode == 5) {
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                mesonAmplitudeMax = mesonAmplitude*400./100.;
            }
        } else {
            mesonAmplitudeMin = mesonAmplitude*30./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            if (fMode == 2 || fMode == 13 || fMode == 3){
                mesonAmplitudeMin = mesonAmplitude*10./100.;
            }
            if (fMode == 4 || fMode == 12 || fMode == 5){
                mesonAmplitudeMin = mesonAmplitude*5./100.;
            }
        }
    } else {
        if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
            mesonAmplitudeMin = mesonAmplitude*98./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            // set global max for pPb
            if (fEnergyFlag.Contains("pPb_5.023TeV") ) mesonAmplitudeMin = mesonAmplitude*92./100.;
            if (fMode == 0 && !fEnergyFlag.CompareTo("8TeV")){
                if ((ptBin > 2)&&ptBin<100 ){
                    fMesonWidthRange[0]         = 0.001;
                    fMesonWidthRange[1]         = 0.009;
                }
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
            }
            if (fMode == 0 && !fEnergyFlag.CompareTo("7TeV")){
                if ((ptBin > 2)&&ptBin<100 ){
                    fMesonWidthRange[0]         = 0.001;
                    fMesonWidthRange[1]         = 0.009;
                }
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
            }
            if (fMode == 2 || fMode == 13 || fMode == 3) {
                mesonAmplitudeMin = mesonAmplitude*98./100.;
                mesonAmplitudeMax = mesonAmplitude*600./100.;
                if( fEnergyFlag.CompareTo("8TeV") == 0 ){
                  mesonAmplitudeMin = mesonAmplitude*90./100.;
                  mesonAmplitudeMax = mesonAmplitude*600./100.;
                }
            }
            if (fMode == 4 || fMode == 12 || fMode == 5) {
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                mesonAmplitudeMax = mesonAmplitude*400./100.;
                if( fEnergyFlag.CompareTo("8TeV") == 0 ){
                  mesonAmplitudeMin = mesonAmplitude*90./100.;
                  mesonAmplitudeMax = mesonAmplitude*400./100.;
                    TString trigger = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
                    if(trigger.CompareTo("52") == 0 && fBinsPt[ptBin]>=13. && fPrefix.CompareTo("Pi0EtaBinning")){
                        fMesonFitRange[0] = 0.08;
                        fMesonFitRange[1] = 0.29;
                    }
                }
            }

        } else {
            mesonAmplitudeMin = mesonAmplitude*50./100.;
            mesonAmplitudeMax = mesonAmplitude*115./100.;
            if (fMode == 0 && !fEnergyFlag.CompareTo("8TeV") ){
                mesonAmplitudeMin = mesonAmplitude*65./100.;
                if(ptBin > 2)  mesonAmplitudeMin = mesonAmplitude*85./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
                if(ptBin < 3) mesonAmplitudeMax = mesonAmplitude*100/100.;
            }
            if (fMode == 2 || fMode == 13 || fMode == 3){
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                if( fEnergyFlag.CompareTo("8TeV") == 0 ){
                  mesonAmplitudeMin = mesonAmplitude*20./100.;
                  mesonAmplitudeMax = mesonAmplitude*115./100.;
                }
            }
            if (fMode == 4 || fMode == 12 || fMode == 5){
                mesonAmplitudeMin = mesonAmplitude*5./100.;
            }
        }
    }

    fFitReco        = NULL;
    fFitReco        = new TF1("GaussExpPol2","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*exp([6]*x))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*exp([6]*x))",fMesonFitRange[0],fMesonFitRange[1]);

    fFitGausExp     = NULL;
    fFitGausExp     = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

    fFitLinearBck   = NULL;
    fFitLinearBck   = new TF1("BGfit","[0]+[1]*exp([2]*x)",fMesonFitRange[0],fMesonFitRange[1]);


    fFitReco->SetParameter(0,mesonAmplitude);
    fFitReco->SetParameter(1,fMesonMassExpect);
    fFitReco->SetParameter(2,fMesonWidthExpect);
    if (fMesonLambdaTail == fMesonLambdaTailRange[0] && fMesonLambdaTail == fMesonLambdaTailRange[1] ){
        fFitReco->FixParameter(3,fMesonLambdaTail);
    } else {
        fFitReco->SetParameter(3,fMesonLambdaTail);
        fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);
    }
    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.15);
    if (fEnergyFlag.Contains("PbPb")){
        fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.5);
    }
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);


    histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");

    fFitReco->SetLineColor(3);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);
    fFitReco->SetNpx(10000);

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
    fFitLinearBck->SetParameter(2,fFitReco->GetParameter(6));

    fFitLinearBck->SetParError(0,fFitReco->GetParError(4));
    fFitLinearBck->SetParError(1,fFitReco->GetParError(5));
    fFitLinearBck->SetParError(2,fFitReco->GetParError(6));

    Int_t binCenterStart;
    Double_t startBinEdge;
    Int_t binCenterEnd;
    Double_t endBinEdge;

    TVirtualFitter * fitter = TVirtualFitter::GetFitter();

    fIntLinearBck = 0;
    fIntLinearBckError = 0;
    if( TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ||
        TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
        binCenterStart  = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+mesonIntDeltaRangeFit[0]);
        startBinEdge    = histoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
        binCenterEnd    = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+mesonIntDeltaRangeFit[1]);
        endBinEdge      = histoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

        Int_t nFreePar      = fFitReco->GetNumberFreeParameters();
        double * covMatrix  = fitter->GetCovarianceMatrix();

        Float_t intBack     = fFitLinearBck->Integral(startBinEdge, endBinEdge);
        Float_t errorBck    = fFitLinearBck->IntegralError(startBinEdge, endBinEdge);

        fFileDataLog << "Parameter for bin " << ptBin << endl;
        fFileDataLog << "Exp1: \t"<<fFitReco->GetParameter(4)<<"+-" << fFitReco->GetParError(4) << "\t "<<fFitReco->GetParameter(5) <<"+-" << fFitReco->GetParError(5)<< "\t "
                     << fFitReco->GetParameter(6) <<"+-" << fFitReco->GetParError(6)<<endl;

        fIntLinearBck       = intBack/histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
        fIntLinearBckError  = errorBck/histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
    } else {
        fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }
}


//****************************************************************************
//************* Fit of Signal+ BG with Gaussian + Linear BG ******************
//****************************************************************************
void FitSubtractedPureGaussianInvMassInPtBins(TH1D* histoMappingSignalInvMassPtBinSingle, Int_t ptBin ){

    histoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude         = histoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin      = 0;
    Double_t mesonAmplitudeMax      = 0;
    if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
        mesonAmplitudeMin = mesonAmplitude*98./100.;
        mesonAmplitudeMax = mesonAmplitude*115./100.;
        if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 || fEnergyFlag.Contains("pPb_5.023TeV") ) mesonAmplitudeMin = mesonAmplitude*92./100.;
        if (fMode == 2 || fMode == 13 || fMode == 3) {
            mesonAmplitudeMin = mesonAmplitude*98./100.;
            mesonAmplitudeMax = mesonAmplitude*400./100.;
        }
        if (fMode == 4 || fMode == 12 || fMode == 5) {
            mesonAmplitudeMin = mesonAmplitude*10./100.;
            mesonAmplitudeMax = mesonAmplitude*400./100.;
            if( fEnergyFlag.CompareTo("8TeV") == 0  ){
                mesonAmplitudeMin = mesonAmplitude*90./100.;
                mesonAmplitudeMax = mesonAmplitude*400./100.;
                TString trigger = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
                if(trigger.CompareTo("52") == 0 && fBinsPt[ptBin]>=13. && fPrefix.CompareTo("Pi0EtaBinning")){
                    fMesonFitRange[0] = 0.08;
                    fMesonFitRange[1] = 0.29;
                }
            } else if (fEnergyFlag.Contains("PbPb")){
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                mesonAmplitudeMax = mesonAmplitude*400./100.;
            }
        }

    } else {
        mesonAmplitudeMin = mesonAmplitude*50./100.;
        mesonAmplitudeMax = mesonAmplitude*120./100.;
        if (fMode == 2 || fMode == 13 || fMode == 3){
            mesonAmplitudeMin = mesonAmplitude*10./100.;
            if( fEnergyFlag.CompareTo("8TeV") == 0  ){
              mesonAmplitudeMin = mesonAmplitude*10./100.;
              mesonAmplitudeMax = mesonAmplitude*200./100.;
            }
        }
        if(fMode == 4 || fMode == 12){
            mesonAmplitudeMin = mesonAmplitude*30./100.;
            mesonAmplitudeMax = mesonAmplitude*220./100.;
        }
    }

    Double_t linBckg = 0.05;
    if(fEnergyFlag.CompareTo("8TeV") == 0 && fPrefix.CompareTo("Eta") == 0) linBckg = 0.1;
    if ( fEnergyFlag.CompareTo("8TeV") == 0 && fPrefix.CompareTo("Pi0") == 0){
        TString trigger         = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
        if( trigger.CompareTo("52") == 0 ){
            linBckg = 0.07;
        }
    }

    fFitReco= NULL;
    fFitReco = new TF1("GaussLinearBG","gaus(0)+[3]+[4]*x",fMesonFitRange[0],fMesonFitRange[1]);

    fFitLinearBck = NULL;
    fFitLinearBck = new TF1("Linear","[0]+[1]*x",fMesonFitRange[1]-linBckg,fMesonFitRange[1]);

    fFitReco->SetParameter(0,mesonAmplitude);
    fFitReco->SetParameter(1,fMesonMassExpect);
    fFitReco->SetParameter(2,fMesonWidthExpect);
    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    fFitReco->SetParLimits(1,fMesonMassExpect*0.8,fMesonMassExpect*1.2);
    if( fEnergyFlag.CompareTo("8TeV") == 0 && fMode == 4 ){
        fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.3);
    } else if( fEnergyFlag.Contains("PbPb") && fMode == 4 ){
        fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.5);
    }
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]*2);

    histoMappingSignalInvMassPtBinSingle->Fit(fFitLinearBck,"QRME0","",fMesonFitRange[1]-linBckg,fMesonFitRange[1]);
    fFitReco->SetParameter(3,fFitLinearBck->GetParameter(0));
    fFitReco->SetParLimits(3,fFitLinearBck->GetParameter(0)-2*fFitLinearBck->GetParError(0),fFitLinearBck->GetParameter(0)+2*fFitLinearBck->GetParError(0));
    fFitReco->SetParameter(4,fFitLinearBck->GetParameter(1));
    fFitReco->SetParLimits(4,fFitLinearBck->GetParameter(1)-2*fFitLinearBck->GetParError(1),fFitLinearBck->GetParameter(1)+2*fFitLinearBck->GetParError(1));
    histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    //exclude second iteration of fitting, otherwise fits go completely wrong in 8 TeV

    fFitReco->SetLineColor(5);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);
    fFitReco->SetNpx(10000);

    fFitLinearBck->SetParameter(0,fFitReco->GetParameter(3));
    fFitLinearBck->SetParameter(1,fFitReco->GetParameter(4));
    fFitLinearBck->SetParError(0,fFitReco->GetParError(3));
    fFitLinearBck->SetParError(1,fFitReco->GetParError(4));

    if(TString(gMinuit->fCstatu.Data()).Contains("CONVERGED") == 1 || TString(gMinuit->fCstatu.Data()).Contains("SUCCESSFUL") == 1 || TString(gMinuit->fCstatu.Data()).Contains("PROBLEMS") == 1){
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
void GausFitSubtractedInvMassInPtBinsNew(TH1D* histoMappingSignalInvMassPtBinSingle,Double_t * mesonIntDeltaRangeFit, Int_t ptBin,Bool_t vary,TString functionname ,Bool_t kMC){
    if(vary){};
    cout <<"Start Fitting spectra"<<endl;

    fCopySignal = (TH1D*)histoMappingSignalInvMassPtBinSingle->Clone("fCopySignal");
    fCopySignal->Sumw2();
    fCopyOnlyBG = (TH1D*)histoMappingSignalInvMassPtBinSingle->Clone("fCopyOnlyBG");
    fCopyOnlyBG->Sumw2();

    fFileErrLog<<"Start Fitting spectra with Gaus fit"<<endl;
    histoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude     = histoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin  = mesonAmplitude*10./100.;
    Double_t mesonAmplitudeMax  = mesonAmplitude*400./100.;

    fFitReco = NULL;
    fFitReco = new TF1(functionname,"(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

    fFitGausExp = NULL;
    fFitGausExp = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

    fFitLinearBck = NULL;
    fFitLinearBck = new TF1("Linear","[0]+[1]*x",fMesonFitRange[0],fMesonFitRange[1]);

    fFitReco->SetParameter(0,mesonAmplitude);
    fFitReco->SetParameter(1,fMesonMassExpect);
    fFitReco->SetParameter(2,fMesonWidthExpect);
    if (fMesonLambdaTail == fMesonLambdaTailRange[0] && fMesonLambdaTail == fMesonLambdaTailRange[1] ){
        fFitReco->FixParameter(3,fMesonLambdaTail);
    } else {
        fFitReco->SetParameter(3,fMesonLambdaTail);
        fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);
    }

    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.15);
    if( fEnergyFlag.CompareTo("8TeV") == 0 && fMode == 4 ){
      fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.3);
    }
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);

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
        }
        fCopySignal->Fit(fFitReco,"QRME0");
    } else {
        histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    }

    fFitReco->SetLineColor(3);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);
    fFitReco->SetNpx(10000);

    fFitGausExp->SetParameter(0,fFitReco->GetParameter(0));
    fFitGausExp->SetParameter(1,fFitReco->GetParameter(1));
    fFitGausExp->SetParameter(2,fFitReco->GetParameter(2));
    fFitGausExp->SetParameter(3,fFitReco->GetParameter(3));

    fFitGausExp->SetParError(0,fFitReco->GetParError(0));
    fFitGausExp->SetParError(1,fFitReco->GetParError(1));
    fFitGausExp->SetParError(2,fFitReco->GetParError(2));
    fFitGausExp->SetParError(3,fFitReco->GetParError(3));

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
        binCenterStart = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+mesonIntDeltaRangeFit[0]);
        startBinEdge = histoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
        binCenterEnd = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+mesonIntDeltaRangeFit[1]);
        endBinEdge = histoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

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

            fIntLinearBck = intLinearBack/histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
            fIntLinearBckError = errorLinearBck/histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

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
void FitTrueInvMassInPtBins(TH1D* histoMappingSignalInvMassPtBinSingle, Double_t* mesonIntDeltaRangeFit, Int_t ptBin, Bool_t vary)
{
    //    cout<<"Start Fitting spectra"<<endl;
    histoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude         = histoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin      = 0;
    Double_t mesonAmplitudeMax      = 0;
    if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 || fEnergyFlag.Contains("XeXe")){
        mesonAmplitudeMin = mesonAmplitude*99./100.;
        mesonAmplitudeMax = mesonAmplitude*110./100.;
        if (fMode == 2 || fMode == 13 || fMode == 3){
            mesonAmplitudeMin = mesonAmplitude*20./100.;
            mesonAmplitudeMax = mesonAmplitude*1000./100.;
        }
        if (fMode == 4 || fMode == 12 || fMode == 5) {
            mesonAmplitudeMin = mesonAmplitude*10./100.;
            mesonAmplitudeMax = mesonAmplitude*800./100.;
        }
    } else {
        mesonAmplitudeMin       = mesonAmplitude*95./100.;
        mesonAmplitudeMax       = mesonAmplitude*130./100.;
        if (fMode == 2 || fMode == 13 || fMode == 3){
            mesonAmplitudeMin   = mesonAmplitude*8./100.;
            mesonAmplitudeMax   = mesonAmplitude*1000./100.;
            if( fEnergyFlag.CompareTo("8TeV") == 0 ){
              mesonAmplitudeMin                     = mesonAmplitude*80./100.;
              mesonAmplitudeMax                     = mesonAmplitude*1000./100.;
            } else if( fEnergyFlag.CompareTo("2.76TeV") == 0 ){
                  mesonAmplitudeMin                 = mesonAmplitude*80./100.;
                  mesonAmplitudeMax                 = mesonAmplitude*1000./100.;
            } else if (fEnergyFlag.Contains("pPb_5.023TeV")  ){
                mesonAmplitudeMin                   = mesonAmplitude*90./100.;
                if(fBinsPt[ptBin] >= 12.0 && fPrefix.Contains("Pi0")) {
                    fMesonLambdaTailMC              = 0.0095;
                    fMesonLambdaTailRangeMC[0]      = 0.0095;
                    fMesonLambdaTailRangeMC[1]      = 0.0095;
                } else if( fPrefix.CompareTo("Eta") == 0) {
                    mesonAmplitudeMin               = mesonAmplitude*70./100.;
                }
            }
        } else if (fMode == 4 || fMode == 12 || fMode == 5) {
            if( fEnergyFlag.CompareTo("8TeV") == 0  ){
                mesonAmplitudeMin = mesonAmplitude*90./100.;
                mesonAmplitudeMax = mesonAmplitude*400./100.;
            } else if (fEnergyFlag.Contains("pPb_5.023TeV")  ){
                mesonAmplitudeMin = mesonAmplitude*90./100.;
                mesonAmplitudeMax = mesonAmplitude*700./100.;
                if( !fPrefix.CompareTo("Eta")) {
                    mesonAmplitudeMin           = mesonAmplitude*80./100.;
                    if (fBinsPt[ptBin] >= 2 ){
                        fMesonLambdaTailMC          = 0.033;
                        fMesonLambdaTailRangeMC[0]  = 0.033;
                        fMesonLambdaTailRangeMC[1]  = 0.033;
                    }
                }
            } else {
                mesonAmplitudeMin = mesonAmplitude*10./100.;
                mesonAmplitudeMax = mesonAmplitude*400./100.;
            }
        }
    }

    fFitReco = NULL;
    TF1* fFitRecoPre = new TF1("fGauss","([0]*exp(-0.5*((x-[1])/[2])^2))", fMesonFitRange[0], fMesonFitRange[1]);
    if (fMode == 2 || fMode == 13 || fMode == 4 || fMode == 12 || fMode == 5 || fMode == 3){
        fFitReco = new TF1("GaussExpLinear","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",
                        fMesonFitRange[0], fMesonFitRange[1]);
    } else {
        fFitReco = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))", fMesonFitRange[0],
                        fMesonFitRange[1]);
    }

    // prefitting the true distributions
    fFitRecoPre->SetParameter(0,mesonAmplitude);
    fFitRecoPre->SetParameter(1,fMesonMassExpect);
    fFitRecoPre->SetParameter(2,fMesonWidthExpect);
    fFitRecoPre->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    fFitRecoPre->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
    if (fMode == 2 || fMode == 13 ){
        fFitRecoPre->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.2);
    } else if(  fMode == 4 || fMode == 12){
        if ( fEnergyFlag.CompareTo("8TeV") == 0 )
            fFitRecoPre->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.3);
        else if ( fEnergyFlag.Contains("PbPb") )
            fFitRecoPre->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.5);
        else if ( fEnergyFlag.Contains("XeXe") )
            fFitRecoPre->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.5);
        else
            fFitRecoPre->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.2);
    }
    histoMappingSignalInvMassPtBinSingle->Fit(fFitRecoPre,"QRME0");



    cout << "MC tail params: " << fMesonLambdaTailMC  << "\t" << fMesonLambdaTailRangeMC[0] << "\t" << fMesonLambdaTailRangeMC[1] << endl;
    if (fMesonLambdaTailMC == fMesonLambdaTailRangeMC[0] && fMesonLambdaTailMC == fMesonLambdaTailRangeMC[1] ){
//         cout << "fixing tail parameter for MC: "<< fBinsPt[ptBin] << endl;
        fFitReco->FixParameter(3,fMesonLambdaTailMC);
    } else {
        if (ptBin > fStartPtBin+1 && fMode == 4 && fEnergyFlag.Contains("pPb_5.023TeV") ){
            fFitReco->SetParameter(3,fMesonLambdaTailMCpar[ptBin-1]);
            fFitReco->SetParLimits(3,fMesonLambdaTailMCpar[ptBin-1]*0.5,fMesonLambdaTailMCpar[ptBin-1]*1.5);
        } else {
            fFitReco->SetParameter(3,fMesonLambdaTailMC);
            fFitReco->SetParLimits(3,fMesonLambdaTailRangeMC[0],fMesonLambdaTailRangeMC[1]);
        }
    }

    // take parameters from prefit
    Double_t mass   = fMesonMassExpect;
    if (fMode == 4 || fMode == 12){
        mass        = fFitRecoPre->GetParameter(1);
        fFitReco->SetParameter(0,fFitRecoPre->GetParameter(0));
        fFitReco->SetParameter(1,mass);
        fFitReco->SetParameter(2,fFitRecoPre->GetParameter(2));
    } else {
        fFitReco->SetParameter(0,mesonAmplitude);
        fFitReco->SetParameter(1,fMesonMassExpect);
        fFitReco->SetParameter(2,fMesonWidthExpect);
    }
    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    if (fMode == 2 || fMode == 13){
        if (fBinsPt[ptBin] > 12 && fEnergyFlag.Contains("pPb_5.023TeV")  && fPrefix.Contains("Pi0") ){
            fFitReco->FixParameter(4,0);
            fFitReco->FixParameter(5,0);
        } else if (fBinsPt[ptBin] > 9 && fEnergyFlag.Contains("pPb_5.023TeV")   ){
            fFitReco->FixParameter(4,0);
            fFitReco->FixParameter(5,0);
        }
    } else if (fMode == 4 || fMode == 12){
        if (fBinsPt[ptBin] > 10 && fEnergyFlag.Contains("pPb_5.023TeV")  ){
            fFitReco->FixParameter(4,0);
            fFitReco->FixParameter(5,0);
        }
    }

    fFitReco->SetNpx(10000);

    if ( !(fMode == 2 || fMode == 13 || fMode == 4 || fMode == 12))
        fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
    else if (fMode == 2 || fMode == 13 )
        fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.2);
    else if (fMode == 4 || fMode == 12 )
        fFitReco->SetParLimits(1,mass*0.95,mass*1.08);

    fFitReco->SetParLimits(2,fMesonWidthRangeMC[0],fMesonWidthRangeMC[1]);

    cout << "amp range: " << mesonAmplitude  << "\t" << mesonAmplitudeMin << "\t" << mesonAmplitudeMax << endl;
    cout << "width range params: " << fMesonWidthExpect  << "\t" << fMesonWidthRangeMC[0] << "\t" << fMesonWidthRangeMC[1] << endl;
    if ( (fMode == 2 || fMode == 4 || fMode == 12 || fMode == 13) && fEnergyFlag.Contains("pPb_5.023TeV")  )
        histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"WLRME0");
    else if ( (fMode == 5 || fMode == 3) && fEnergyFlag.CompareTo("pPb_5.023TeVRun2") == 0  )
        histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"WLRME0");
    else if ( fMode == 2  && fEnergyFlag.CompareTo("2.76TeV") == 0 )
        histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"WLRME0");
    else if ( (fMode == 2 || fMode == 4 || fMode == 3 || fMode == 5 || fMode == 12 || fMode == 13) && fEnergyFlag.Contains("XeXe")  )
        histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"WLRME0");
    else
        histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    if (vary && (fMode==9 || fMode ==0)){
        fMesonLambdaTailMC = fFitReco->GetParameter(3);
    }
    fFitReco->SetLineColor(3);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);

    if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ){
        fFileDataLog << "Parameter for bin " << ptBin << endl;
        fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
    } else {
        fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }
    fFitReco->DrawCopy("same");
    if (mesonIntDeltaRangeFit){}
}

//****************************************************************************
//*** Fit of Pure MC Signal with Gaussian ************************************
//****************************************************************************
void FitTrueInvMassPureGaussianInPtBins(TH1D* histoMappingSignalInvMassPtBinSingle, Int_t ptBin ){

    histoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude         = histoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin      = 0;
    Double_t mesonAmplitudeMax      = 0;
    if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
        mesonAmplitudeMin = mesonAmplitude*98./100.;
        mesonAmplitudeMax = mesonAmplitude*115./100.;
        if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 || fEnergyFlag.Contains("pPb_5.023TeV") ) mesonAmplitudeMin = mesonAmplitude*92./100.;
        if (fMode == 2 || fMode == 13 || fMode == 3) {
            mesonAmplitudeMin = mesonAmplitude*98./100.;
            mesonAmplitudeMax = mesonAmplitude*400./100.;
        }
        if (fMode == 4 || fMode == 12 || fMode == 5) {
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
        if (fMode == 2 || fMode == 13 || fMode == 3){
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
    if( fEnergyFlag.CompareTo("8TeV") == 0 && fMode == 4 ){
      fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.3);
    }
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]*2);

    histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    //exclude second iteration of fitting, otherwise fits go completely wrong in 8 TeV

    fFitReco->SetLineColor(5);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);

    if(TString(gMinuit->fCstatu.Data()).Contains("CONVERGED") == 1 || TString(gMinuit->fCstatu.Data()).Contains("SUCCESSFUL") == 1 || TString(gMinuit->fCstatu.Data()).Contains("PROBLEMS") == 1){
        fFileDataLog << "Parameter for pure Gaussian bin " << ptBin << endl;
        fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<<endl;
    } else {
        fFileErrLog << "Pure Gaussian fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }
}


//****************************************************************************
//*** Fit of subtracted Signal+ BG with CrystalBall tail + Lin BG ************
//*** linear BG subtracted in this function after initial fit without ********
//*** peak region, final fit only with CrystalBall ***************************
//*** additional outputs: fCopySignal - only Signal                      ******
//***                     fCopyOnlyBG - only remaining BG                 ******
//****************************************************************************
void FitCBSubtractedInvMassInPtBins(TH1D* histoMappingSignalInvMassPtBinSingle,Double_t * mesonIntDeltaRangeFit, Int_t ptBin,Bool_t vary ,TString functionname, Bool_t kMC)
{
    if(vary){}; //dummy case to remove warning

    fCopySignal = (TH1D*)histoMappingSignalInvMassPtBinSingle->Clone("fCopySignal");
    fCopySignal->Sumw2();
    fCopyOnlyBG = (TH1D*)histoMappingSignalInvMassPtBinSingle->Clone("fCopyOnlyBG");
    fCopyOnlyBG->Sumw2();

    fFileErrLog<<"Start Fitting spectra with CB fit"<<endl;
    histoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude         = histoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin      = mesonAmplitude*50./100.;
    Double_t mesonAmplitudeMax      = mesonAmplitude*200./100.;
    if (fMode == 4 || fMode == 12)
        mesonAmplitudeMin           = mesonAmplitude*5./100.;

    fFitReco = NULL;
//    if (!kMC) {
    //   fFitReco = new TF1(functionname,CrystalBallBck,fMesonFitRange[0],fMesonFitRange[1],7);
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
        //  cout << fFitLinearBck->Eval(startBinEdge) << "\t" <<fFitLinearBck->Eval(endBinEdge) << "\t" <<fCopySignal->GetBinContent(i) << "\t" <<fCopySignal->GetBinContent(i) << endl;
        }
        fCopySignal->Fit(fFitReco,"QRME0");
    } else {
        histoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    }

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
        binCenterStart = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+mesonIntDeltaRangeFit[0]);
        startBinEdge = histoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterStart)- 0.5*histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
        binCenterEnd = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+mesonIntDeltaRangeFit[1]);
        endBinEdge = histoMappingSignalInvMassPtBinSingle->GetBinCenter(binCenterEnd)+ 0.5*histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

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

            fIntLinearBck = intLinearBack/histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
            fIntLinearBckError = errorLinearBck/histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
            fFileDataLog << "Integrated BG: \t" << intLinearBack << "+-" <<  errorLinearBck << "\t bin width" <<histoMappingSignalInvMassPtBinSingle->GetBinWidth(10) <<endl;

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
void IntegrateHistoInvMass(TH1D * histoMappingSignalInvMassPtBinSingle, Double_t * fMesonIntRangeInt)
{
    Int_t binLowMassMeson = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[0]);
    Int_t binHighMassMeson = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[1]);
    fYields = histoMappingSignalInvMassPtBinSingle->IntegralAndError(binLowMassMeson,binHighMassMeson,fYieldsError);
}


//****************************************************************************
//*** Integration of Invariant Mass Histogram in given integration window ****
//*** with detailed output to log file ***************************************
//****************************************************************************
void IntegrateHistoInvMassStream(TH1D * histoMappingSignalInvMassPtBinSingle, Double_t * fMesonIntRangeInt) {
    Int_t binLowMassMeson = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[0]);
    Int_t binHighMassMeson = histoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[1]);
    fYields = histoMappingSignalInvMassPtBinSingle->IntegralAndError(binLowMassMeson,binHighMassMeson,fYieldsError);
    for ( Int_t M = binLowMassMeson; M < binHighMassMeson+1; M++){
        fFileDataLog << M << "\t" << histoMappingSignalInvMassPtBinSingle->GetBinCenter(M) <<"\t" <<histoMappingSignalInvMassPtBinSingle->GetBinContent(M)<< "+-"<< histoMappingSignalInvMassPtBinSingle->GetBinError(M)<< endl;
    }
}

//****************************************************************************
//********* Integration of Fit function in given integration window **********
//****************************************************************************
void IntegrateFitFunc(TF1 * fFunc, TH1D *  histoMappingSignalInvMassPtBinSingle,Double_t * fMesonIntRangeInt) {
    fYieldsFunc = fFunc->Integral(fMesonIntRangeInt[0],fMesonIntRangeInt[1])/histoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
}


//****************************************************************************
//***************** Filling of MC histograms in proper binning ***************
//****************************************************************************
void FillHistosArrayMC(TH1D* fHistoMCMesonPtWithinAcceptanceFill, TH1D * fHistoMCMesonPtFill, TH1D * fDeltaPtFill) {
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
    fHistoMCMesonPtWithinAcceptanceFill->Sumw2();
    fHistoMCMesonWithinAccepPtWOWeights = (TH1D*)fHistoMCMesonPtWithinAcceptanceFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
    fHistoMCMesonWithinAccepPtWOWeights->Divide(fDeltaPtFill);
    fHistoMCMesonPtFill->Sumw2();
    fHistoMCMesonPt1WOWeights = (TH1D*)fHistoMCMesonPtFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
    fHistoMCMesonPt1WOWeights->Divide(fDeltaPtFill);

}

//****************************************************************************
//***************** Filling of MC histograms in proper binning ***************
//****************************************************************************
void FillHistosArrayMCWOEvtWeights(TH1D* fHistoMCMesonPtWithinAcceptanceFill, TH1D * fHistoMCMesonPtFill, TH1D * fDeltaPtFill) {
    fHistoMCMesonPtWithinAcceptanceFill->Sumw2();
    fHistoMCMesonWithinAccepPtWOEvtWeights = (TH1D*)fHistoMCMesonPtWithinAcceptanceFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
    fHistoMCMesonWithinAccepPtWOEvtWeights->Divide(fDeltaPtFill);
    fHistoMCMesonPtFill->Sumw2();
    fHistoMCMesonPt1WOEvtWeights = (TH1D*)fHistoMCMesonPtFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
    fHistoMCMesonPt1WOEvtWeights->Divide(fDeltaPtFill);

}

//****************************************************************************
//***************** Filling of MC secondary histograms ***********************
//***************** rebin to proper binning & calculate acc ******************
//****************************************************************************
void FillMCSecondaryHistAndCalculateAcceptance(TH2D* mcSecInputSourcePt, TH2D* mcSecInputInAccSourcePt){
    // MC histo pt vs sourcs
    // 1-> K0s
    // 2-> Lambda
    // 3-> K0L
    // 5-> p
    // 6-> n
    // 7-> pi
    // 8-> K
    // 9-> rho
    // 10-> Delta
    // 11-> K*
    // 15-> Rest
    Int_t startBin[4]   = {1,2,3,4};
    Int_t endBin[4]     = {1,2,3,15};

    for (Int_t j = 0; j < 4; j++){
        // project correct bin
        fHistoMCSecPi0Pt[j]           = (TH1D*)mcSecInputSourcePt->ProjectionX( Form("MCSecPi0From%s",nameSecondaries[j].Data()), mcSecInputSourcePt->GetYaxis()->FindBin(startBin[j]),
                                                                                mcSecInputSourcePt->GetYaxis()->FindBin(endBin[j]),"e");
        fHistoMCSecPi0Pt[j]->SetTitle(Form("MCSecPi0From%s",nameSecondaries[j].Data()));
        fHistoMCSecPi0Pt[j]->Sumw2();
        fHistoMCSecPi0PtWAcc[j]       = (TH1D*)mcSecInputInAccSourcePt->ProjectionX(Form("MCSecPi0From%s_InAcc",nameSecondaries[j].Data()), mcSecInputInAccSourcePt->GetYaxis()->FindBin(startBin[j]),
                                                                                    mcSecInputInAccSourcePt->GetYaxis()->FindBin(endBin[j]),"e");
        fHistoMCSecPi0PtWAcc[j]->SetTitle(Form("MCSecPi0From%s_InAcc",nameSecondaries[j].Data()));
        fHistoMCSecPi0PtWAcc[j]->Sumw2();

        // rebin in to current analysis binning
        fHistoMCSecPi0PtReb[j]          = (TH1D*)fHistoMCSecPi0Pt[j]->Rebin(fNBinsPt,Form("MCSecPi0From%s_Rebinned",nameSecondaries[j].Data()),fBinsPt); // Proper bins in Pt
        fHistoMCSecPi0PtReb[j]->SetTitle(Form("MCSecPi0From%s_Rebinned",nameSecondaries[j].Data()));
        fHistoMCSecPi0PtReb[j]->Divide(fDeltaPt);
        fHistoMCSecPi0PtWAccReb[j]      = (TH1D*)fHistoMCSecPi0PtWAcc[j]->Rebin(fNBinsPt,Form("MCSecPi0From%s_InAcc_Rebinned",nameSecondaries[j].Data()),fBinsPt); // Proper bins in Pt
        fHistoMCSecPi0PtWAccReb[j]->SetTitle(Form("MCSecPi0From%s_InAcc_Rebinned",nameSecondaries[j].Data()));
        fHistoMCSecPi0PtWAccReb[j]->Divide(fDeltaPt);

        // calculate acceptance
        fHistoMCSecPi0AcceptPt[j]       = (TH1D*)fHistoMCSecPi0PtWAccReb[j]->Clone(Form("fMCSecPi0From%sAccepPt",nameSecondaries[j].Data()));
        fHistoMCSecPi0AcceptPt[j]->SetTitle(Form("fMCSecPi0From%sAccepPt",nameSecondaries[j].Data()));
        fHistoMCSecPi0AcceptPt[j]->Divide(fHistoMCSecPi0PtWAccReb[j],fHistoMCSecPi0PtReb[j],1.,1.,"B");
    }
}


//****************************************************************************
//***************** Calculation of Meson Acceptance **************************
//****************************************************************************
void CalculateMesonAcceptance() {
    cout << "calculating acceptance" << endl;
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
    cout << "calculating acceptance w/o weights" << endl;
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
//***************** Calculation of Meson Acceptance **************************
//****************************************************************************
void CalculateMesonAcceptanceWOEvtWeights() {
    cout << "calculating acceptance w/o event weights" << endl;
    fHistoMCMesonAcceptPtWOEvtWeights = new TH1D("fMCMesonAccepPtWOEvtWeights","",fNBinsPt,fBinsPt);
    fHistoMCMesonAcceptPtWOEvtWeights->Sumw2();

    fHistoMCMesonAcceptPtWOEvtWeights->Divide(fHistoMCMesonWithinAccepPtWOEvtWeights,fHistoMCMesonPt1WOEvtWeights,1.,1.,"B");
    fHistoMCMesonAcceptPtWOEvtWeights->DrawCopy();
    fFileDataLog << endl << "Calculation of the Acceptance wo weights" << endl;
    for ( Int_t i = 1; i < fHistoMCMesonAcceptPtWOEvtWeights->GetNbinsX()+1 ; i++){
        fFileDataLog << "Bin " << i << "\t"<< fHistoMCMesonAcceptPtWOEvtWeights->GetBinCenter(i)<< "\t" << fHistoMCMesonAcceptPtWOEvtWeights->GetBinContent(i) << "\t" << fHistoMCMesonAcceptPtWOEvtWeights->GetBinError(i) <<endl;
    }
}

//****************************************************************************
//***************** Calculation of Meson Efficiency **************************
//****************************************************************************
TH1D* CalculateMesonEfficiency(TH1D* MC_fMesonYieldsPt, TH1D** MC_SecondaryYieldPt, TH1D* histAcc, TString nameEfi ) {

    // create histo with proper binning
    TH1D* histoMCMesonEffiPt = new TH1D(nameEfi.Data(),"",fNBinsPt,fBinsPt);
    histoMCMesonEffiPt->Sumw2();
    // add original reconstructed yields
    histoMCMesonEffiPt->Add(MC_fMesonYieldsPt,1.);
    // subtract secondary yield to get only primary efficiency if wanted
    if (MC_SecondaryYieldPt){
        for (Int_t j = 0; j<4; j++){
            if(MC_SecondaryYieldPt[j]) histoMCMesonEffiPt->Add(MC_SecondaryYieldPt[j],-1.);
        }
    }
    // devide by MC input yield in acceptance
    histoMCMesonEffiPt->Divide(histoMCMesonEffiPt, histAcc, 1.,1.,"B");
    // write efficiency to output file as text
    fFileDataLog << endl << "Calculation of the Efficiency: " << nameEfi.Data()<< endl;
    for ( Int_t i = 1; i < histoMCMesonEffiPt->GetNbinsX()+1 ; i++){
        fFileDataLog << "Bin " << i << "\t" << histoMCMesonEffiPt->GetBinCenter(i)<< "\t"<< histoMCMesonEffiPt->GetBinContent(i) << "\t" << histoMCMesonEffiPt->GetBinError(i) <<endl;
    }
    // return pointer to hist
    return histoMCMesonEffiPt;
}

//****************************************************************************
//****** Saving of general histograms, fits needed in the analysis ***********
//****** RAW output file, no correction histograms ***************************
//****************************************************************************
void SaveHistos(Int_t optionMC, TString cutID, TString prefix3, Bool_t UseTHnSparse ) {
    TString nameOutput = Form("%s/%s/%s_%s_GammaConvV1WithoutCorrection%s_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),prefix3.Data(),fPeriodFlag.Data(),cutID.Data());
    fOutput1 = new TFile(nameOutput.Data(),"RECREATE");

    cout << "----------------------------------------------------------------------------------" << endl;
    cout << "Created output file for rec data: " << nameOutput.Data() << endl;
    cout << "----------------------------------------------------------------------------------" << endl;
    cout << "Begin writing Uncorrected File" << endl;

    TH1D*   deltaPtCluster       = new TH1D("fDeltaPtCluster","",fNBinsClusterPt,fBinsClusterPt);
    for(Int_t iPt=1;iPt<fNBinsClusterPt+1;iPt++){
        deltaPtCluster->SetBinContent(iPt,fBinsClusterPt[iPt]-fBinsClusterPt[iPt-1]);
        deltaPtCluster->SetBinError(iPt,0);
    }

    if (fHistoClustersPt){
        cout << "writing ClusterPt" << endl;
        fHistoClustersPt->Write("ClusterPt");
        TGraphErrors* graphClusterPt        = new TGraphErrors(fHistoClustersPt);
//         graphClusterPt->Print();
        TH1D*   fHistoClustersPtPerEvent   = (TH1D*)fHistoClustersPt->Rebin(fNBinsClusterPt,"fHistoClustersPtPerEvent",fBinsClusterPt);
        fHistoClustersPtPerEvent->Divide(deltaPtCluster);
        fHistoClustersPtPerEvent->Scale(1./fNEvents);
        fHistoClustersPtPerEvent->Write("ClusterPtPerEvent");
        TGraphErrors* graphClusterPtReb        = new TGraphErrors(fHistoClustersPtPerEvent);
//         graphClusterPtReb->Print();
    }
    if (fHistoClustersE){
        cout << "writing ClusterE" << endl;
        fHistoClustersE->Write("ClusterE");

        TH1D*   fHistoClustersEPerEvent   = (TH1D*)fHistoClustersE->Rebin(fNBinsClusterPt,"fHistoClustersEPerEvent",fBinsClusterPt);
        fHistoClustersEPerEvent->Divide(deltaPtCluster);
        fHistoClustersEPerEvent->Scale(1./fNEvents);
        fHistoClustersEPerEvent->Write("ClusterEPerEvent");
    }
    if (fEnableDCCluster){
        cout << "writing ClusterPtReb" << endl;
        TH1D*   fHistoTrueGammaClusPtRebinned   = NULL;
        TH1D*   fHistoTrueGammaDCClusPtRebinned   = NULL;
        if (fHistoTrueGammaClusPt){
            fHistoTrueGammaClusPt->Write();
            fHistoTrueGammaClusPtRebinned   = (TH1D*)fHistoTrueGammaClusPt->Rebin(fNBinsClusterPt,"fHistoTrueGammaClusPtRebinned",fBinsClusterPt);
        }
        if (fHistoTrueGammaDCClusPt){
            fHistoTrueGammaDCClusPt->Write();
            fHistoTrueGammaDCClusPtRebinned   = (TH1D*)fHistoTrueGammaDCClusPt->Rebin(fNBinsClusterPt,"fHistoTrueGammaDCClusPtRebinned",fBinsClusterPt);
        }
        if (fHistoTrueGammaClusMultipleCount)fHistoTrueGammaClusMultipleCount->Write();
        if (fHistoTrueGammaClusPt && fHistoTrueGammaDCClusPt){
            TH1F* fHistoRatioDCTrueGammaClus = (TH1F*)fHistoTrueGammaDCClusPt->Clone("fHistoRatioDCTrueGammaClus");
            fHistoRatioDCTrueGammaClus->Divide(fHistoRatioDCTrueGammaClus,fHistoTrueGammaClusPt,1,1,"B");
            fHistoRatioDCTrueGammaClus->Write("FractionDoubleCountedClusters");
            TH1F* fHistoRatioDCTrueGammaClusRebinned = (TH1F*)fHistoTrueGammaDCClusPtRebinned->Clone("fHistoRatioDCTrueGammaClusRebinned");
            fHistoRatioDCTrueGammaClusRebinned->Divide(fHistoRatioDCTrueGammaClusRebinned,fHistoTrueGammaClusPtRebinned,1,1,"B");
            fHistoRatioDCTrueGammaClusRebinned->Write("FractionDoubleCountedClustersRebinned");
        }
    }

    if (fHistoClustersOverlapHeadersPt){
        fHistoClustersOverlapHeadersPt->Write("ClusterOverlapHeadersPt");
        TH1D*   fHistoClustersOverlapHeadersPtPerEvent   = (TH1D*)fHistoClustersOverlapHeadersPt->Rebin(fNBinsClusterPt,"fHistoClustersOverlapHeadersPtPerEvent",fBinsClusterPt);
        fHistoClustersOverlapHeadersPtPerEvent->Divide(deltaPtCluster);
        fHistoClustersOverlapHeadersPtPerEvent->Scale(1./fNEvents);
        fHistoClustersOverlapHeadersPtPerEvent->Write("ClusterOverlapHeadersPtPerEvent");
    }

    // write histograms for all integration windows: normal, wide, narrow, left, left wide, left narrow
    for (Int_t k = 0; k < 6; k++){
        if (fHistoYieldMeson[k])            fHistoYieldMeson[k]->Write();
        if (fHistoYieldMesonPerEvent[k])    fHistoYieldMesonPerEvent[k]->Write();
    }

    // write histograms for assumption of different backgrounds
    for (Int_t k = 0; k < 6; k++){
        if (fHistoYieldDiffBck[k]) fHistoYieldDiffBck[k]->Write();
        if (fHistoYieldDiffBckRatios[k]) fHistoYieldDiffBckRatios[k]->Write();
    }
    for (Int_t k = 0; k < 3; k++){
        if (fHistoYieldDiffBckResult[k]) fHistoYieldDiffBckResult[k]->Write();
    }


    // write histograms for integration windows: normal, wide, narrow
    for (Int_t k = 0; k < 3; k++){
        if (fHistoMassWindowHigh[k])        fHistoMassWindowHigh[k]->Write();
        if (fHistoMassWindowLow[k])         fHistoMassWindowLow[k]->Write();
        if (fHistoSigndefaultMeson[k])      fHistoSigndefaultMeson[k]->Write();
        if (fHistoSBdefaultMeson[k])        fHistoSBdefaultMeson[k]->Write();
    }
    fHistoLambdaTail->Write();
    fHistoAmplitude->Write();
    fHistoSigma->Write();
    fHistoResidualBGcon->Write();
    fHistoResidualBGlin->Write();
    fHistoRatioResBGYield->Write();
    fHistoRatioResBGYieldToSPlusResBG->Write();
    for (Int_t m = 0; m < 4; m++){
        fHistoChi2[m]->Write();
        fHistoResBGYield[m]->Write();
    }

    fHistoMassMeson->Write();
    fHistoMassGaussianMeson->Write();
    fHistoWidthGaussianMeson->Write();
    fHistoFWHMMeson->Write();
    fDeltaPt->Write();


    fHistoMassMesonLeft->Write();
    fHistoFWHMMesonLeft->Write();
    fMesonFullPtSignal->Write();
    fMesonFullPtBackground->Write();
    fMesonFullPtBackNorm->SetName("Mapping_BackNorm_InvMass_FullPt");
    fMesonFullPtBackNorm->Write();
    fNumberOfGoodESDTracks->Write();
    fEventQuality->Write();

    TString nameHistoSignal;
    TString titleHistoSignal;
    TString nameHistoSignalLeft;
    TString nameHistoBckNorm;
    TString fitnameSignal;
    TString fitnameSignalOther;
    TString nameHistoBckNormLeft;
    TString fitnameSignalLeft;

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
        nameHistoBckNorm    = Form("Mapping_BckNorm_InvMass_in_Pt_Bin%02d", ii);
        nameHistoSignal     = Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d", ii);
        nameHistoSignalLeft = Form("fHistoMappingSignalInvMassLeft_in_Pt_Bin%02d", ii); //Added 12.12.2014
        fitnameSignal       = Form("Signal_InvMassFit_in_Pt_Bin%02d", ii);
        fitnameSignalOther   = Form("Signal_InvMassFitWithOther_in_Pt_Bin%02d", ii);
        titleHistoSignal    = Form("%3.1f GeV/#it{c} < #it{p}_{T} < %3.1f GeV/#it{c}",fBinsPt[ii],fBinsPt[ii+1]),
        fHistoMappingBackNormInvMassPtBin[ii]->Write(nameHistoBckNorm.Data());
        fHistoMappingSignalInvMassPtBin[ii]->SetTitle(titleHistoSignal.Data());
        fHistoMappingSignalInvMassPtBin[ii]->Write(nameHistoSignal.Data());
        fHistoMappingSignalInvMassLeftPtBin[ii]->Write(nameHistoSignalLeft.Data());//Added 12.12.2014
        if( fMode == 4 || fMode == 12 ) {
            nameHistoBckNormLeft = Form("Mapping_BckNormLeft_InvMass_in_Pt_Bin%02d", ii);
            fHistoMappingBackNormInvMassLeftPtBin[ii]->Write(nameHistoBckNormLeft.Data());
            fitnameSignalLeft = Form("SignalnewLeft_InvMassFit_in_Pt_Bin%02d", ii);
            if(fFitInvMassLeftPtBin[ii]!=0x00) fFitInvMassLeftPtBin[ii]->Write(fitnameSignalLeft.Data());
        }
        if (fFitSignalInvMassPtBin[ii]!=0x00) fFitSignalInvMassPtBin[ii]->Write(fitnameSignal.Data());
        fitnameSignalOther   = Form("Signal_InvMassFitWithPol2_in_Pt_Bin%02d", ii);
        if (fFitSignalWithOtherBGInvMassPtBin[0][ii] != 0x00) fFitSignalWithOtherBGInvMassPtBin[0][ii]->Write(fitnameSignalOther.Data());
        fitnameSignalOther   = Form("Signal_InvMassFitWithExp1_in_Pt_Bin%02d", ii);
        if (fFitSignalWithOtherBGInvMassPtBin[1][ii] != 0x00) fFitSignalWithOtherBGInvMassPtBin[1][ii]->Write(fitnameSignalOther.Data());
        fitnameSignalOther   = Form("Signal_InvMassFitWithExp2_in_Pt_Bin%02d", ii);
        if (fFitSignalWithOtherBGInvMassPtBin[2][ii] != 0x00) fFitSignalWithOtherBGInvMassPtBin[2][ii]->Write(fitnameSignalOther.Data());

        if (fCrysFitting==1 ){
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
        fHistoTrueMesonInvMassVSPt->Write();

        for (Int_t k = 0; k < 3; k++){ // different integration windows: normal, wide, narrow
            if (fHistoYieldTrueMeson[k])            fHistoYieldTrueMeson[k]->Write();
            if (fHistoYieldTrueMesonReweighted[k])  fHistoYieldTrueMesonReweighted[k]->Write();
            if (fHistoYieldTrueMesonUnweighted[k])  fHistoYieldTrueMesonUnweighted[k]->Write();
            for (Int_t j = 0; j<4; j++){ // different secondary types: K0s, Lambda, K0L, rest
                if(fHistoYieldTrueSecMeson[k][j])    fHistoYieldTrueSecMeson[k][j]->Write();
            }
        }
        for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
            fHistoMappingTrueMesonInvMassPtBins[ii]->Write();
            if (fEnableDCMeson) fHistoMappingTrueMesonDCInvMassPtBins[ii]->Write();
            fHistoMappingTrueFullMesonInvMassPtBins[ii]->Write();
            fHistoMappingTrueMesonInvMassPtReweightedBins[ii]->Write();
            fHistoMappingTrueMesonInvMassPtUnweightedBins[ii]->Write();
            if (fAdvancedMesonQA){
                if (fMode == 2 || fMode == 13 || fMode == 3 || fMode == 4 || fMode == 12 || fMode == 5){
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
                if (fMode == 4 || fMode == 12 || fMode == 5){
                    if (fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[ii])fHistoMappingTrueMesonMixedCaloConvPhotonInvMassPtBins[ii]->Write();
                    if (fFitTrueSignalMixedCaloConvPhotonInvMassPtBin[ii]!=0x00) fFitTrueSignalMixedCaloConvPhotonInvMassPtBin[ii]->Write();
                }

            }
            if(fMode == 4 || fMode == 12 || fMode == 5){
                if(fHistoMappingTrueGGBckInvMassPtBins[ii]!=0x00 && fHistoMappingTrueContBckInvMassPtBins[ii]!=0x00 && fHistoMappingTrueAllBckInvMassPtBins[ii]!=0x00){
                    fHistoMappingTrueGGBckInvMassPtBins[ii]->Write();
                    if(fHistoMappingTrueMesonContainedInvMassPtBins[ii]!=0x00)fHistoMappingTrueMesonContainedInvMassPtBins[ii]->Write();
                    if(fHistoMappingTrueAsymEClusInvMassPtBins[ii]!=0x00)fHistoMappingTrueAsymEClusInvMassPtBins[ii]->Write();
                    fHistoMappingTrueContBckInvMassPtBins[ii]->Write();
                    fHistoMappingTrueAllBckInvMassPtBins[ii]->Write();
                }
            }
            for (Int_t j = 0; j < 4; j++){
                if (fHistoMappingTrueSecMesonInvMassPtBins[j][ii] != 0x00) fHistoMappingTrueSecMesonInvMassPtBins[j][ii]->Write();
            }
            if (fFitTrueSignalInvMassPtBin[ii]!=0x00) fFitTrueSignalInvMassPtBin[ii]->Write();
            if (fFitTrueSignalInvMassPtReweightedBin[ii]!=0x00) fFitTrueSignalInvMassPtReweightedBin[ii]->Write();
            if (fFitTrueSignalInvMassPtUnweightedBin[ii]!=0x00) fFitTrueSignalInvMassPtUnweightedBin[ii]->Write();
        }
    }


    if (fHaveToyMCInputForSec||fHaveCocktailInputForSec){
        cout << "Writing input for secondary pi0 correction" << endl;
        for (Int_t j = 0; j < 3; j++){
            if (fHistoYieldExternSecInput[j]){
                cout << "writing external sec input: " << nameSecondaries[j].Data() << endl;
                fHistoYieldExternSecInput[j]->Write();
            }
            if (fHistoYieldExternSecInputReb[j])     fHistoYieldExternSecInputReb[j]->Write();
        }

        cout << "Writing input for pi0 resonance feed down correction" << endl;
        for (Int_t j = 0; j < 15; j++){
            if (fHistoYieldExternResonanceFeedDownInput[j])        fHistoYieldExternResonanceFeedDownInput[j]->Write();
            if (fHistoYieldExternResonanceFeedDownInputReb[j])     fHistoYieldExternResonanceFeedDownInputReb[j]->Write();
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
void SaveCorrectionHistos(TString cutID, TString prefix3){
    const char* nameOutput = Form("%s/%s/%s_%s_GammaConvV1CorrectionHistos%s_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),prefix3.Data(),fPeriodFlag.Data(),cutID.Data());
    fOutput2 = new TFile(nameOutput,"RECREATE");
    cout<<"======================================================"<<endl;
    cout<<"======================================================"<<endl;
    cout<<"======================================================"<<endl;
    cout<<nameOutput<<endl;
    cout<<"======================================================"<<endl;
    cout<<"======================================================"<<endl;
    cout<<"======================================================"<<endl;

    cout << "Begin writing Correction File" << endl;

    // write event counting histo
    if (fEventQuality)          fEventQuality->Write();

    // write acceptance
    if (fHistoMCMesonAcceptPt)  fHistoMCMesonAcceptPt->Write();
    // write input spectra with weights
    if (fHistoMCMesonPt1){
        fHistoMCMesonPt1->SetName("MC_Meson_genPt");
        fHistoMCMesonPt1->Write(); // Proper bins in Pt
    }
    if (fHistoMCMesonPt){
        fHistoMCMesonPt->SetName("MC_Meson_genPt_oldBin");
	DivideTH1ByBinWidth(fHistoMCMesonPt);
        fHistoMCMesonPt->Write();
    }

    // write input spectra w/o particle weights
    if (fHistoMCMesonPtWOWeights){
        fHistoMCMesonPtWOWeights->SetName("MC_Meson_genPt_WOWeights");
	DivideTH1ByBinWidth(fHistoMCMesonPtWOWeights);
        fHistoMCMesonPtWOWeights->Write();
    }
    // write particle weights
    if (fHistoMCMesonPtWeights){
        fHistoMCMesonPtWeights->Write("MC_Meson_genPt_Weights");
    }

    // write acceptances & input spectra w/o particle weights
    if (fHistoMCMesonPtWithinAcceptanceWOWeights){
        fHistoMCMesonAcceptPtWOWeights->Write();
        fHistoMCMesonPt1WOWeights->SetName("MC_Meson_genPt_properBinning_WOWeights");
        fHistoMCMesonPt1WOWeights->Write(); // Proper bins in Pt
    }
    // write acceptances & input spectra w/o event weights
    if (fHistoMCMesonPtWithinAcceptanceWOEvtWeights){
        fHistoMCMesonAcceptPtWOEvtWeights->Write();
        fHistoMCMesonPt1WOEvtWeights->SetName("MC_Meson_genPt_properBinning_WOEvtWeights");
        fHistoMCMesonPt1WOEvtWeights->Write(); // Proper bins in Pt
    }

    // write efficiencies depending on integration windows: normal, wide, narrow, left, left wide, left narrow
    for (Int_t k = 0; k < 6; k++){
        if (fHistoMonteMesonEffiPt[k])          fHistoMonteMesonEffiPt[k]->Write();
        if (k < 3){
            if (fHistoMCTrueMesonEffiPt[k])                 fHistoMCTrueMesonEffiPt[k]->Write();
            if (fHistoMCTrueMesonEffiPtReweighted[k])       fHistoMCTrueMesonEffiPtReweighted[k]->Write();
            if (fHistoMCMesonPtWithinAcceptanceWOWeights){
                if (fHistoMCTrueMesonEffiPtUnweighted[k])   fHistoMCTrueMesonEffiPtUnweighted[k]->Write();
            }
        }
    }

    // Write new secondary pi0 MC quantities
    if( fNewMCOutput ){
        if (fHistoMCSecPi0SourcePt)         fHistoMCSecPi0SourcePt->Write();
        if (fHistoMCSecPi0WAccSourcePt)     fHistoMCSecPi0WAccSourcePt->Write();
        for (Int_t j = 0; j < 4; j++){
            if(fHistoMCSecPi0Pt[j])         fHistoMCSecPi0Pt[j]->Write();
            if(fHistoMCSecPi0PtWAcc[j])     fHistoMCSecPi0PtWAcc[j]->Write();
            if(fHistoMCSecPi0PtReb[j])      fHistoMCSecPi0PtReb[j]->Write();
            if(fHistoMCSecPi0PtWAccReb[j])  fHistoMCSecPi0PtWAccReb[j]->Write();
            if(fHistoMCSecPi0AcceptPt[j])   fHistoMCSecPi0AcceptPt[j]->Write();
        }
        for (Int_t k = 0; k < 3; k++){
            for (Int_t j = 0; j < 4; j++){
                if (fHistoMCTrueSecMesonEffiPt[k][j])   fHistoMCTrueSecMesonEffiPt[k][j]->Write();
            }
        }
    }

    // write different mass and width histos reco
    if (fHistoMassMeson)            fHistoMassMeson->Write();
    if (fHistoFWHMMeson)            fHistoFWHMMeson->Write();
    if (fHistoMassGaussianMeson)    fHistoMassGaussianMeson->Write();
    if (fHistoWidthGaussianMeson)   fHistoWidthGaussianMeson->Write();

    // write different mass and width histos validated reco
    if (fHistoTrueMassMeson) fHistoTrueMassMeson->Write();
    if (fHistoTrueMassGaussianMeson)    fHistoTrueMassGaussianMeson->Write();
    if (fHistoTrueWidthGaussianMeson)   fHistoTrueWidthGaussianMeson->Write();
    if (fHistoTrueMassMesonReweighted)  fHistoTrueMassMesonReweighted->Write();
    if (fHistoTrueMassMesonUnweighted)  fHistoTrueMassMesonUnweighted->Write();
    if (fHistoTrueFWHMMeson)            fHistoTrueFWHMMeson->Write();
    if (fHistoTrueFWHMMesonReweighted)  fHistoTrueFWHMMesonReweighted->Write();
    if (fHistoTrueFWHMMesonUnweighted)  fHistoTrueFWHMMesonUnweighted->Write();
    if (fAdvancedMesonQA && (fMode == 2 || fMode == 13 || fMode == 3 || fMode == 4 || fMode == 12 || fMode == 5)){
        if (fHistoTrueMassMesonCaloPhoton)                  fHistoTrueMassMesonCaloPhoton->Write();
        if (fHistoTrueMassMesonCaloElectron)                fHistoTrueMassMesonCaloElectron->Write();
        if (fHistoTrueMassMesonCaloConvPhoton)              fHistoTrueMassMesonCaloConvPhoton->Write();
        if (fHistoTrueMassMesonCaloMergedCluster)           fHistoTrueMassMesonCaloMergedCluster->Write();
        if (fHistoTrueMassMesonCaloMergedPartConvCluster)   fHistoTrueMassMesonCaloMergedPartConvCluster->Write();
        if (fHistoTrueFWHMMesonCaloPhoton)                  fHistoTrueFWHMMesonCaloPhoton->Write();
        if (fHistoTrueFWHMMesonCaloElectron)                fHistoTrueFWHMMesonCaloElectron->Write();
        if (fHistoTrueFWHMMesonCaloConvPhoton)              fHistoTrueFWHMMesonCaloConvPhoton->Write();
        if (fHistoTrueFWHMMesonCaloMergedCluster)           fHistoTrueFWHMMesonCaloMergedCluster->Write();
        if (fHistoTrueFWHMMesonCaloMergedPartConvCluster)   fHistoTrueFWHMMesonCaloMergedPartConvCluster->Write();
        if (fHistoYieldTrueMesonFixedWindow)                fHistoYieldTrueMesonFixedWindow->Write();
        if (fHistoYieldTrueMesonGammaFixedWindow)           fHistoYieldTrueMesonGammaFixedWindow->Write();
        if (fHistoYieldTrueMesonGammaConvGammaFixedWindow)  fHistoYieldTrueMesonGammaConvGammaFixedWindow->Write();
        if (fHistoYieldTrueMesonConvGammaConvGammaFixedWindow)fHistoYieldTrueMesonConvGammaConvGammaFixedWindow->Write();
    }
    if (fAdvancedMesonQA && (fMode == 4 || fMode == 12 || fMode == 5)){
        if (fHistoTrueMassMesonMixedCaloConvPhoton)         fHistoTrueMassMesonMixedCaloConvPhoton->Write();
        if (fHistoTrueFWHMMesonMixedCaloConvPhoton)         fHistoTrueFWHMMesonMixedCaloConvPhoton->Write();
    }

    // write reconstructed yield for reconstructed mesons for different integration ranges: normal, wide, narrow
    for (Int_t k = 0; k < 3; k++){
        // primary yield
        if (fHistoYieldTrueMeson[k])            fHistoYieldTrueMeson[k]->Write();
        if (fHistoYieldTrueMesonReweighted[k])  fHistoYieldTrueMesonReweighted[k]->Write();
        if (fHistoYieldTrueMesonUnweighted[k])  fHistoYieldTrueMesonUnweighted[k]->Write();
        // secondary yield
        for (Int_t j = 0; j<4; j++){
            if (fHistoYieldTrueSecMeson[k][j])fHistoYieldTrueSecMeson[k][j]->Write();
            if (fHistoYieldTrueSecFracMeson[k][j])fHistoYieldTrueSecFracMeson[k][j]->Write();
        }
    }

    // write double counting histos
    if (fEnableDCMeson){
        if (fHistoYieldTrueMesonDC)             fHistoYieldTrueMesonDC->Write();
        if (fHistoYieldTrueMesonDC){
            TH1D* fHistoRatioDCTrueMeson    = (TH1D*)fHistoYieldTrueMesonDC->Clone("fHistoRatioDCTrueMeson");
            fHistoRatioDCTrueMeson->Divide(fHistoRatioDCTrueMeson, fHistoYieldTrueMeson[0], 1, 1, "B");
            fHistoRatioDCTrueMeson->Write("FractionDoubleCountedMesons");
        }
        if (fHistoTrueMesonMultipleCount) fHistoTrueMesonMultipleCount->Write();
    }

    // write mother distributions for reconstructed mesons
    if (fHistoYieldK0sWithPi0DaughterRec)       fHistoYieldK0sWithPi0DaughterRec->Write("K0sWithPi0DaughterRec");
    if (fHistoYieldLambdaWithPi0DaughterRec)    fHistoYieldLambdaWithPi0DaughterRec->Write("LambdaWithPi0DaughterRec");

    // write 2D invariant mass histo validate meson
    if (fHistoTrueMesonInvMassVSPt)             fHistoTrueMesonInvMassVSPt->Write();

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
//******** definition of Gaussian to fit SPD pileup distribution  ************
//****************************************************************************
Double_t fitGaussianPileUp(Double_t *x, Double_t *par)
{
    if (x[0] > -0.9 && x[0] < 0.9) {
        TF1::RejectPoint();
        return 0;
    }
   return par[0]*(exp(-0.5*pow(((x[0]-par[1])/par[2]),2)));
}

Double_t fitGaussianPileUp2(Double_t *x, Double_t *par)
{
    if (x[0] > -0.2 && x[0] < 0.2) {
        TF1::RejectPoint();
        return 0;
    }
    if (x[0] < -0.8 || x[0] > 0.8) {
        TF1::RejectPoint();
        return 0;
    }
   return par[0]*(exp(-0.5*pow(((x[0]-par[1])/par[2]),2)));
}

//****************************************************************************
//****** Convert heavy meson string to identification digit ******************
//****************************************************************************
Int_t GetHeavyMesonDigit(TString mesonString) {
    if(mesonString.CompareTo("Pi0"))      return 0;
    if(mesonString.CompareTo("Eta"))      return 1;
    if(mesonString.CompareTo("EtaPrime")) return 2;
    // If invalid mode was chosen
    std::cout << "Not chosen a valid particle name (\"" << mesonString << "\")" << std::endl;
    return -1;
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

        if (fMode == 3){
            fFWHMFunc = fFunc->GetParameter(2)*2.35;
            fFWHMFuncError = fFunc->GetParError(2)*2.35;
        }
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

    for (Int_t m = 0; m< 3; m++){
        if (fMesonYieldsResBckOtherFunc[m])                     delete fMesonYieldsResBckOtherFunc[m];
        if (fMesonYieldsResBckOtherFuncError[m])                delete fMesonYieldsResBckOtherFuncError[m];
    }

    for (Int_t k = 0; k < 6; k++){
        // delete arrays for yields
        if (fGGYields[k])                                       delete fGGYields[k];
        if (fBckYields[k])                                      delete fBckYields[k];
        if (fMesonYields[k])                                    delete fMesonYields[k];
        if (fMesonYieldsFunc[k])                                delete fMesonYieldsFunc[k];
        if (fMesonYieldsResidualBckFunc[k])                     delete fMesonYieldsResidualBckFunc[k];
        if (fMesonYieldsCorResidualBckFunc[k])                  delete fMesonYieldsCorResidualBckFunc[k];
        if (fMesonYieldsPerEvent[k])                            delete fMesonYieldsPerEvent[k];

        // delete arrays for errors
        if (fGGYieldsError[k])                                  delete fGGYieldsError[k];
        if (fBckYieldsError[k])                                 delete fBckYieldsError[k];
        if (fMesonYieldsError[k])                               delete fMesonYieldsError[k];
        if (fMesonYieldsFuncError[k])                           delete fMesonYieldsFuncError[k];
        if (fMesonYieldsResidualBckFuncError[k])                delete fMesonYieldsResidualBckFuncError[k];
        if (fMesonYieldsCorResidualBckFuncError[k])             delete fMesonYieldsCorResidualBckFuncError[k];
        if (fMesonYieldsPerEventError[k])                       delete fMesonYieldsPerEventError[k];

        // delete mass range array
        if (fMesonCurIntRange[k])                               delete fMesonCurIntRange[k];

    }
    for (Int_t k = 0; k < 3; k++){
        // delete mass window arrays
        if (fMassWindowHigh[k])                                 delete fMassWindowHigh[k];
        if (fMassWindowLow[k])                                  delete fMassWindowLow[k];

        // delete integration window arrays
        if (fMesonTrueIntRange[k])                              delete fMesonTrueIntRange[k];
        if (fMesonTrueIntReweightedRange[k])                    delete fMesonTrueIntReweightedRange[k];
        if (fMesonTrueIntUnweightedRange[k])                    delete fMesonTrueIntUnweightedRange[k];

        // delete true meson yield arrays
        if (fMesonTrueYields[k])                                delete fMesonTrueYields[k];
        if (fMesonTrueYieldsReweighted[k])                      delete fMesonTrueYieldsReweighted[k];
        if (fMesonTrueYieldsUnweighted[k])                      delete fMesonTrueYieldsUnweighted[k];

        // delete array for S/B and Significance
        if (fMesonSBdefault[k])                                 delete fMesonSBdefault[k];
        if (fMesonSigndefault[k])                               delete fMesonSigndefault[k];
        if (fMesonSBdefaultError[k])                            delete fMesonSBdefaultError[k];
        if (fMesonSigndefaultError[k])                          delete fMesonSigndefaultError[k];

    }
    if (fMesonTrueYieldsDC)                                     delete fMesonTrueYieldsDC;
    if (fMesonTrueYieldFixedWindow)                             delete fMesonTrueYieldFixedWindow;
    if (fMesonTrueYieldGammaConvGammaFixedWindow)               delete fMesonTrueYieldGammaConvGammaFixedWindow;
    if (fMesonTrueYieldConvGammaConvGammaFixedWindow)           delete fMesonTrueYieldConvGammaConvGammaFixedWindow;
    if (fMesonTrueYieldGammaFixedWindow)                        delete fMesonTrueYieldGammaFixedWindow;
    if (fMesonTrueYieldErrorFixedWindow)                        delete fMesonTrueYieldErrorFixedWindow;
    if (fMesonTrueYieldGammaErrorFixedWindow)                   delete fMesonTrueYieldGammaErrorFixedWindow;
    if (fMesonTrueYieldGammaConvGammaErrorFixedWindow)          delete fMesonTrueYieldGammaConvGammaErrorFixedWindow;
    if (fMesonTrueYieldConvGammaConvGammaErrorFixedWindow)      delete fMesonTrueYieldConvGammaConvGammaErrorFixedWindow;
    if (fMesonMass)                                             delete fMesonMass;
    if (fMesonLambdaTailpar)                                    delete fMesonLambdaTailpar;
    if (fMesonLambdaTailparError)                               delete fMesonLambdaTailparError;
    if (fMesonLambdaTailMCpar)                                  delete fMesonLambdaTailMCpar;
    if (fMesonLambdaTailMCparError)                             delete fMesonLambdaTailMCparError;
    if (fMesonSigmapar)                                         delete fMesonSigmapar;
    if (fMesonSigmaparError)                                    delete fMesonSigmaparError;
    if (fMesonTrueSigmapar)                                     delete fMesonTrueSigmapar;
    if (fMesonTrueSigmaparError)                                delete fMesonTrueSigmaparError;
    if (fMesonAmplitudepar)                                     delete fMesonAmplitudepar;
    if (fMesonAmplitudeparError)                                delete fMesonAmplitudeparError;
    if (fMesonResidualBGlin)                                    delete fMesonResidualBGlin;
    if (fMesonResidualBGlinError)                               delete fMesonResidualBGlinError;
    if (fMesonResidualBGcon)                                    delete fMesonResidualBGcon;
    if (fMesonResidualBGconError)                               delete fMesonResidualBGconError;
    for (Int_t m = 0; m < 4; m++){
        if (fMesonChi2[m])                                      delete fMesonChi2[m];
    }

    if (fMesonTrueSB)                                           delete fMesonTrueSB;
    if (fMesonTrueSign)                                         delete fMesonTrueSign;
    if (fMesonFWHM)                                             delete fMesonFWHM;
    if (fMesonMassLeft)                                         delete fMesonMassLeft;
    if (fMesonFWHMLeft)                                         delete fMesonFWHMLeft;

    if (fMesonMassError)                                        delete fMesonMassError;
    if (fMesonTrueSBError)                                      delete fMesonTrueSBError;
    if (fMesonTrueSignError)                                    delete fMesonTrueSignError;
    if (fMesonFWHMError)                                        delete fMesonFWHMError;
    if (fMesonMassLeftError)                                    delete fMesonMassLeftError;
    if (fMesonFWHMLeftError)                                    delete fMesonFWHMLeftError;
    if (fHistoMappingTrueMesonInvMassPtBins)                    delete fHistoMappingTrueMesonInvMassPtBins;
    if (fHistoMappingTrueMesonDCInvMassPtBins)                  delete fHistoMappingTrueMesonDCInvMassPtBins;
    if (fHistoMappingTrueFullMesonInvMassPtBins)                delete fHistoMappingTrueFullMesonInvMassPtBins;
    if (fHistoMappingTrueMesonInvMassPtReweightedBins)          delete fHistoMappingTrueMesonInvMassPtReweightedBins;
    if (fHistoMappingTrueMesonInvMassPtUnweightedBins)          delete fHistoMappingTrueMesonInvMassPtUnweightedBins;
    if (fHistoMappingTrueGGBckInvMassPtBins)                    delete fHistoMappingTrueGGBckInvMassPtBins;
    if (fHistoMappingTrueContBckInvMassPtBins)                  delete fHistoMappingTrueContBckInvMassPtBins;
    if (fHistoMappingTrueAllBckInvMassPtBins)                   delete fHistoMappingTrueAllBckInvMassPtBins;
    if (fHistoMappingTrueMesonContainedInvMassPtBins)           delete fHistoMappingTrueMesonContainedInvMassPtBins;
    if (fHistoMappingTrueAsymEClusInvMassPtBins)                delete fHistoMappingTrueAsymEClusInvMassPtBins;
    if (fHistoMappingGGInvMassPtBin)                            delete fHistoMappingGGInvMassPtBin;
    if (fHistoMappingBackNormAndRemainingBGInvMassPtBin)        delete fHistoMappingBackNormAndRemainingBGInvMassPtBin;
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
    if (fHistoWeightsBGZbinVsMbin)                              delete fHistoWeightsBGZbinVsMbin;
    if (fHistoFillPerEventBGZbinVsMbin)                         delete fHistoFillPerEventBGZbinVsMbin;
    if (fHistoWeightsBGZbinVsPsibin)                            delete fHistoWeightsBGZbinVsPsibin;
    if (fHistoFillPerEventBGZbinVsPsibin)                       delete fHistoFillPerEventBGZbinVsPsibin;
    for (Int_t m = 0; m < 3; m++){
        if (fFitSignalWithOtherBGInvMassPtBin[m])               delete fFitSignalWithOtherBGInvMassPtBin[m];

        if (fFitBckOtherInvMassPtBin[m])                        delete fFitBckOtherInvMassPtBin[m];
    }
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

    for (Int_t j = 0; j < 3; j++){
        if (fFileToyMCInput[j] )                                   delete fFileToyMCInput[j];
    }
    if (fFileCocktailInput)                                     delete fFileCocktailInput;

}

// MAIN FUNCTION for non-ROOT compilation (prototype)
int main( int argc, char* argv[] )
{
    // Default arguments for ExtractSignal
        TString meson                   = "";
        TString file                    = "";
        TString cutSelection            = "";
        TString Suffix                  = "";
        TString optionMC                = "";
        TString optionEnergy            = "";
        TString optionCrystalBall       = "";
        TString directphotonPlots       = "";
        TString optionUseMinBiasEff     = "";
        TString optionPeriod            = "";
        TString optionAdvancedMesonQA   = "";
        Int_t numberOfBins              = 30;
        Bool_t addSig                   = kFALSE;
        Int_t mode                      = 9;
        Bool_t UseTHnSparse             = kTRUE;
        Int_t triggerSet                = -1;

    // Import main call arguments
        TString import;
        if( argc >  1 ) meson                 = argv[1];
        if( argc >  2 ) file                  = argv[2];
        if( argc >  3 ) cutSelection          = argv[3];
        if( argc >  4 ) Suffix                = argv[4];
        if( argc >  5 ) optionMC              = argv[5];
        if( argc >  6 ) optionEnergy          = argv[6];
        if( argc >  7 ) optionCrystalBall     = argv[7];
        if( argc >  8 ) directphotonPlots     = argv[8];
        if( argc >  9 ) optionUseMinBiasEff   = argv[9];
        if( argc > 10 ) optionPeriod          = argv[10];
        if( argc > 11 ) optionAdvancedMesonQA = argv[11];
        if( argc > 12 ) { // numberOfBins
            istringstream sstr(argv[12]);
            sstr >> numberOfBins;
        } if( argc > 13 ) { // addSig
            import = argv[13];
            if( import.EqualTo("kTRUE") )  addSig = kTRUE;
            if( import.EqualTo("kFALSE") ) addSig = kFALSE;
        } if( argc > 14 ) { // mode
            istringstream sstr(argv[14]);
            sstr >> mode;
        } if( argc > 15 ) { // UseTHnSparse
            import = argv[15];
            if( import.EqualTo("kTRUE") )  UseTHnSparse = kTRUE;
            if( import.EqualTo("kFALSE") ) UseTHnSparse = kFALSE;
        } if( argc > 16 ) { // triggerSet
            istringstream sstr(argv[16]);
            sstr >> triggerSet;
        }
    return 0;

    // Function call ExtractSignalV2
        ExtractSignalV2(
            meson,
            file,
            cutSelection,
            Suffix,
            optionMC,
            optionEnergy,
            optionCrystalBall,
            directphotonPlots,
            optionUseMinBiasEff,
            optionPeriod,
            optionAdvancedMesonQA,
            numberOfBins,
            addSig,
            mode,
            UseTHnSparse,
            triggerSet
        );
}