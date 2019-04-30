// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion
// Author: Friederike Bock, fbock@cern.ch
// This version is supposed to be used for files produced with AliPhysics vAN-20160527-1 and newer

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
#include "ExtractSignalMergedMesonV2.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "../CommonHeaders/ExtractSignalPlotting.h"
#include "THnSparse.h"


//****************************************************************************
//********* Main function for extraction of signal for merged mesons *********
//****************************************************************************
void ExtractSignalMergedMesonV2(    TString meson                   = "",
                                    TString file                    = "",
                                    TString cutSelection            = "",
                                    TString suffix                  = "",
                                    Bool_t  optionMC                = kFALSE,
                                    TString optionEnergy            = "",
                                    TString optionPeriod            = "",
                                    TString optionAdvancedMesonQA   = "",
                                    Int_t   numberOfBins            = 30,
                                    Int_t   mode                    = 10,
                                    Bool_t  kJJGammaTriggMC         = 0,
                                    Int_t   triggerSet              = -1,
                                    TString optionCorrFrameworkDir  = ""
                                ) {
    gROOT->Reset();

    //*****************************************************************************************************************
    //******************************* Set global variables ************************************************************
    //*****************************************************************************************************************
    fMode                           = mode;
    // mode:   10 // new output merged EMCal
    //         11 // new output merged PHOS
    if (!(fMode == 10 || fMode == 11)){
        cout << "ERROR: You are running the wrong macro, this macro is only designed for merged cluster analysis." << endl;
        return;
    }

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

    fDoExtendedExamplePlots = kTRUE;

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
        Initialize("Pi0", numberOfBins, triggerSet);
    } else if (meson.CompareTo("Eta") == 0) {
        Initialize("Eta", numberOfBins, triggerSet);
    } else if(meson.CompareTo("Pi0EtaBinning") == 0) {
        Initialize("Pi0EtaBinning", numberOfBins, triggerSet);
    } else   {
        cout<<"ERROR: First argument in the ExtractSignalMergedMeson(....) has to be either Pi0 or Eta or Pi0EtaBinning  or EtaPrime"<<endl;
        return;
    }
    TString fEventCutSelectionRead          = fEventCutSelection.Data();
    TString fClusterCutSelectionRead        = fClusterCutSelection.Data();
    TString fClusterMergedCutSelectionRead  = fClusterMergedCutSelection.Data();
    TString fMesonCutSelectionRead          = fMesonCutSelection.Data();

    //************************* Start of Main routine ***************************************************************
    const char* fFileErrLogDatname  = Form("%s/%s/%s_%s_FileErrLog%s_%s.dat", cutSelection.Data(), fEnergyFlag.Data(), fPrefix.Data(), fPrefix2.Data(),
                                           fPeriodFlag.Data(), fCutSelectionRead.Data());
    fFileErrLog.open(fFileErrLogDatname, ios::out);

    TFile* f                        = new TFile(file.Data());
    TString autoDetectedMainDir     = AutoDetectMainTList(mode , f,"",optionCorrFrameworkDir);
    if (autoDetectedMainDir.CompareTo("") == 0){
        cout << "ERROR: trying to read file, which is incompatible with mode selected" << endl;;
        return;
    }

    TList *TopDir                   = (TList*)f->Get(autoDetectedMainDir.Data());
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
    fHistoClustersE                     = (TH1D*)ESDContainer->FindObject("ClusGamma_E");
    if (fIsMC){
        fHistoClustersOverlapHeadersPt  = (TH1D*)ESDContainer->FindObject("ClusGammaOverlapHeaders_Pt");
    }
    fHistoClustersMergedPtM02           = (TH2F*)ESDContainer->FindObject("ClusMerged_Pt_M02");
    fHistoClustersMergedPtM02AccMeson   = (TH2F*)ESDContainer->FindObject("ClusMerged_Pt_M02_AcceptedMeson");
    fHistoClustersMergedEM02AccMeson    = (TH2F*)ESDContainer->FindObject("ClusMerged_E_M02_AcceptedMeson");
    if (fHistoClustersMergedEM02AccMeson)
        fHistoClustersMergedE           = (TH1D*)fHistoClustersMergedEM02AccMeson->ProjectionX("ClusMerged_E");

    fHistoClusterCandidates             = (TH1D*)ESDContainer->FindObject("GammaCandidates");
    fHistoClusterMergedCandidates       = (TH1D*)ESDContainer->FindObject("MergedCandidates");
    fHistoTrackVsClusterCandidates      = (TH2F*)ESDContainer->FindObject("GoodESDTracksVsGammaCandidates");
    fHistoSPDtrackletvsSPDclusters      = (TH2F*)ESDContainer->FindObject("SPD tracklets vs SPD clusters");

    fHistoClusterNCellsPt               = (TH2F*)ESDContainer->FindObject("Clus_NCells_Pt");
    fHistoClusterMergedNCellsPt         = (TH2F*)ESDContainer->FindObject("ClusMerged_NCells_Pt");
    fHistoClusterMergedNCellsArClPt     = (TH2F*)ESDContainer->FindObject("ClusMerged_NCellsAroundClus_Pt");
    fHistoClusterMergedNCellsArAInclPt  = (TH2F*)ESDContainer->FindObject("ClusMerged_NCellsAroundAndInClus_Pt");
    fHistoClusterMergedEAroundClE       = (TH2F*)ESDContainer->FindObject("ClusMerged_EAroundClus_E");
    fHistoClustersMergedNLMPt           = (TH2F*)ESDContainer->FindObject("ClusMerged_NLM_Pt_AcceptedMeson");
    if (fHistoClustersMergedNLMPt){
        for (Int_t i = 0; i< 11; i++ ){
            fHistoClustersMergedPtNLMBins[i] = FillProjectionY(fHistoClustersMergedNLMPt, Form("ClusterMergedNLM%d",i+1), i+1, i+1, fNBinsClusterPt, fBinsClusterPt);
        }
    }
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
    if (fEnergyFlag.Contains("PbPb") || fEnergyFlag.Contains("pPb") ){
        fNEvents                        = fEventQuality->GetBinContent(1);
    } else {
        fNEvents                        =  GetNEvents(fEventQuality);
    }

    cout << Form("CaloCuts_%s",fClusterCutSelection.Data()) << endl;
    TList *CaloPreCuts                  = (TList*) HistosGammaConversion->FindObject(Form("CaloCuts_%s",fClusterCutSelection.Data()));
    if(CaloPreCuts == NULL){
        cout << "couldn't find ClusterCuts output, disabling Exotics QA" << endl;
        fEnableExoticsQA                = kFALSE;
    } else {
        fHistoExoticsEtaPhi             = (TH2F*)CaloPreCuts->FindObject(Form("EtaPhi_Exotics %s",fClusterCutSelection.Data()));
        if (fHistoExoticsEtaPhi) cout << "found EtaPhi for exotics" << endl;
        fHistoExoticsEvsM02             = (TH2F*)CaloPreCuts->FindObject(Form("EVsM02_Exotics %s",fClusterCutSelection.Data()));
        if (fHistoExoticsEvsM02) cout << "found EvsM02 for exotics" << endl;
        fHistoExoticsEvsNCells          = (TH2F*)CaloPreCuts->FindObject(Form("ClusterEnergyVsNCells_Exotics %s",fClusterCutSelection.Data()));
        if (fHistoExoticsEvsNCells) cout << "found EvsNCells for exotics" << endl;
        fHistoExoticsEvsEstar           = (TH2F*)CaloPreCuts->FindObject(Form("ClusterEnergyVsEnergystar_Exotics %s",fClusterCutSelection.Data()));
        if (fHistoExoticsEvsEstar) cout << "found EvEstar for exotics" << endl;
        if (fHistoExoticsEtaPhi && fHistoExoticsEvsM02 && fHistoExoticsEvsNCells && fHistoExoticsEvsEstar)
            fEnableExoticsQA            = kTRUE;

        if (fEnableExoticsQA){
            FillExoticsM02HistosArray(fHistoExoticsEvsM02);
        }


    }

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

    FillDataHistosArray( fHistoInvMassVsPt, fHistoClustersMergedPtM02AccMeson );

    if ( fIsMC ){
        SetCorrectMCHistogrammNames(meson);

        TList *TrueContainer                    = (TList*) HistosGammaConversion->FindObject(Form("%s True histograms",fCutSelectionRead.Data()));
        fHistoTrueClustersMergedPtM02           = (TH2F*) TrueContainer->FindObject(fObjectNameTrueMergedM02.Data());
        if (fHistoTrueClustersMergedPtM02) fHistoTrueClustersMergedPtM02->Sumw2();
        fHistoTrueClustersPi0PtM02              = (TH2F*) TrueContainer->FindObject(fObjectNameTrueFromPi0M02.Data());
        if (fHistoTrueClustersPi0PtM02) fHistoTrueClustersPi0PtM02->Sumw2();
        fHistoTrueClustersPi0DCPtM02            = (TH2F*) TrueContainer->FindObject(fObjectNameTrueFromPi0DCM02.Data());
        if (fHistoTrueClustersPi0DCPtM02) fHistoTrueClustersPi0DCPtM02->Sumw2();
        fHistoTrueClustersPi0DalitzPtM02        = (TH2F*) TrueContainer->FindObject(fObjectNameTrueFromPi0DalitzM02.Data());
        if (fHistoTrueClustersPi0DalitzPtM02) fHistoTrueClustersPi0DalitzPtM02->Sumw2();
        if (fHistoTrueClustersPi0PtM02 && fHistoTrueClustersPi0DalitzPtM02){
            fHistoTrueClustersPi0GGPtM02        = (TH2F*) fHistoTrueClustersPi0PtM02->Clone(fObjectNameTrueFromPi0GGM02.Data());
            fHistoTrueClustersPi0GGPtM02->Sumw2();
            fHistoTrueClustersPi0GGPtM02->Add(fHistoTrueClustersPi0DalitzPtM02,-1);
        }

        fHistoTrueClustersEtaPtM02              = (TH2F*) TrueContainer->FindObject(fObjectNameTrueFromEtaM02.Data());
        if (fHistoTrueClustersEtaPtM02) fHistoTrueClustersEtaPtM02->Sumw2();
        fHistoTrueClustersEtaDCPtM02            = (TH2F*) TrueContainer->FindObject(fObjectNameTrueFromEtaDCM02.Data());
        if (fHistoTrueClustersEtaDCPtM02) fHistoTrueClustersEtaDCPtM02->Sumw2();
        fHistoTrueClustersEtaDalitzPtM02        = (TH2F*) TrueContainer->FindObject(fObjectNameTrueFromEtaDalitzM02.Data());
        if (fHistoTrueClustersEtaDalitzPtM02) fHistoTrueClustersEtaDalitzPtM02->Sumw2();
        if (fHistoTrueClustersEtaPtM02 && fHistoTrueClustersEtaDalitzPtM02){
            fHistoTrueClustersEtaGGPtM02        = (TH2F*) fHistoTrueClustersEtaPtM02->Clone(fObjectNameTrueFromEtaGGM02.Data());
            fHistoTrueClustersEtaGGPtM02->Sumw2();
            fHistoTrueClustersEtaGGPtM02->Add(fHistoTrueClustersEtaDalitzPtM02,-1);
        }

        fHistoTrueClustersBGPtM02               = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusBGM02.Data());
        if (fHistoTrueClustersBGPtM02) fHistoTrueClustersBGPtM02->Sumw2();
        fHistoTrueClustersGammaPtM02            = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusGammaM02.Data());
        if (fHistoTrueClustersGammaPtM02) fHistoTrueClustersGammaPtM02->Sumw2();
        fHistoTrueClustersElectronPtM02         = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusElectronM02.Data());
        if (fHistoTrueClustersElectronPtM02) fHistoTrueClustersElectronPtM02->Sumw2();
        fHistoTrueClusPartConvMergedPtM02       = (TH2F*) TrueContainer->FindObject(fObjectNameTrueMergedPartConvM02.Data());
        if (fHistoTrueClusPartConvMergedPtM02){
            fHistoTrueClusPartConvMergedPtM02->Sumw2();
        } else {
            fHistoTrueClusPartConvMergedFromPi0PtM02    = (TH2F*) TrueContainer->FindObject(fObjectNameTrueMergedPartConvFromPi0M02.Data());
            if (fHistoTrueClusPartConvMergedFromPi0PtM02){
                fHistoTrueClusPartConvMergedFromPi0PtM02->Sumw2();
                fHistoTrueClusPartConvMergedPtM02       = (TH2F*) fHistoTrueClusPartConvMergedFromPi0PtM02->Clone(fObjectNameTrueMergedPartConvM02.Data());
                fHistoTrueClusPartConvMergedPtM02->Sumw2();
            }
            fHistoTrueClusPartConvMergedFromEtaPtM02    = (TH2F*) TrueContainer->FindObject(fObjectNameTrueMergedPartConvFromEtaM02.Data());
            if (fHistoTrueClusPartConvMergedFromEtaPtM02){
                fHistoTrueClusPartConvMergedFromEtaPtM02->Sumw2();
                if (fHistoTrueClusPartConvMergedPtM02){
                    fHistoTrueClusPartConvMergedPtM02->Add(fHistoTrueClusPartConvMergedFromEtaPtM02);
                }

            }
        }
        fHistoTrueClusPureMergedPtM02           = (TH2F*) TrueContainer->FindObject(fObjectNameTrueMergedPureM02.Data());
        if (fHistoTrueClusPureMergedPtM02){
            fHistoTrueClusPureMergedPtM02->Sumw2();
        } else {
            fHistoTrueClusPureMergedFromPi0PtM02        = (TH2F*) TrueContainer->FindObject(fObjectNameTrueMergedPureFromPi0M02.Data());
            if (fHistoTrueClusPureMergedFromPi0PtM02){
                fHistoTrueClusPureMergedFromPi0PtM02->Sumw2();
                fHistoTrueClusPureMergedPtM02           = (TH2F*) fHistoTrueClusPureMergedFromPi0PtM02->Clone(fObjectNameTrueMergedPureM02.Data());
                fHistoTrueClusPureMergedPtM02->Sumw2();
            }
            fHistoTrueClusPureMergedFromEtaPtM02        = (TH2F*) TrueContainer->FindObject(fObjectNameTrueMergedPureFromEtaM02.Data());
            if (fHistoTrueClusPureMergedFromEtaPtM02){
                fHistoTrueClusPureMergedFromEtaPtM02->Sumw2();
                if (fHistoTrueClusPureMergedPtM02){
                    fHistoTrueClusPureMergedPtM02->Add(fHistoTrueClusPureMergedFromEtaPtM02);
                }

            }
        }
        fHistoTrueClusOneGammaFromPi0PtM02      = (TH2F*) TrueContainer->FindObject(fObjectNameTrueMergedOneGammaFromPi0M02.Data());
        if (fHistoTrueClusOneGammaFromPi0PtM02) fHistoTrueClusOneGammaFromPi0PtM02->Sumw2();
        fHistoTrueClusOneGammaFromEtaPtM02      = (TH2F*) TrueContainer->FindObject(fObjectNameTrueMergedOneGammaFromEtaM02.Data());
        if (fHistoTrueClusOneGammaFromEtaPtM02) fHistoTrueClusOneGammaFromEtaPtM02->Sumw2();
        fHistoTrueClusOneElectronFromPi0PtM02   = (TH2F*) TrueContainer->FindObject(fObjectNameTrueMergedOneElectronFromPi0M02.Data());
        if (fHistoTrueClusOneElectronFromPi0PtM02) fHistoTrueClusOneElectronFromPi0PtM02->Sumw2();
        fHistoTrueClusOneElectronFromEtaPtM02   = (TH2F*) TrueContainer->FindObject(fObjectNameTrueMergedOneElectronFromEtaM02.Data());
        if (fHistoTrueClusOneElectronFromEtaPtM02) fHistoTrueClusOneElectronFromEtaPtM02->Sumw2();

        if (fHistoTrueClusOneGammaFromPi0PtM02){
            fHistoTrueClusOneGammaPtM02             = (TH2F*)fHistoTrueClusOneGammaFromPi0PtM02->Clone(fObjectNameTrueMergedOneGammaM02.Data());
            fHistoTrueClusOneGammaPtM02->Sumw2();
            if (fHistoTrueClusOneGammaFromEtaPtM02)
                fHistoTrueClusOneGammaPtM02->Add(fHistoTrueClusOneGammaFromEtaPtM02);
        }
        if (fHistoTrueClusOneElectronFromPi0PtM02){
            fHistoTrueClusOneElectronPtM02          = (TH2F*)fHistoTrueClusOneElectronFromPi0PtM02->Clone(fObjectNameTrueMergedOneElectronM02.Data());
            fHistoTrueClusOneElectronPtM02->Sumw2();
            if (fHistoTrueClusOneElectronFromEtaPtM02)
                fHistoTrueClusOneElectronPtM02->Add(fHistoTrueClusOneElectronFromEtaPtM02);
        }


        FillMCM02HistosArray(       fHistoTrueClustersMergedPtM02, fHistoTrueClustersPi0PtM02, fHistoTrueClustersEtaPtM02, fHistoTrueClustersGammaPtM02,
                                    fHistoTrueClustersElectronPtM02, fHistoTrueClustersBGPtM02
                            );
        FillMCM02AdditionHistosArray(   fHistoTrueClusPureMergedPtM02, fHistoTrueClusPartConvMergedPtM02, fHistoTrueClusOneGammaPtM02, fHistoTrueClusOneGammaFromPi0PtM02,
                                        fHistoTrueClusOneGammaFromEtaPtM02, fHistoTrueClusOneElectronPtM02, fHistoTrueClusOneElectronFromPi0PtM02, fHistoTrueClusOneElectronFromEtaPtM02,
                                        fHistoTrueClustersPi0GGPtM02, fHistoTrueClustersPi0DalitzPtM02, fHistoTrueClustersEtaGGPtM02, fHistoTrueClustersEtaDalitzPtM02,
                                        fHistoTrueClustersPi0DCPtM02, fHistoTrueClustersEtaDCPtM02, fHistoTrueClusPureMergedFromPi0PtM02, fHistoTrueClusPureMergedFromEtaPtM02,
                                        fHistoTrueClusPartConvMergedFromPi0PtM02, fHistoTrueClusPartConvMergedFromEtaPtM02
                                    );


        fHistoTrueClustersBGPtSource            = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusBG_Source.Data());
        fHistoTrueClustersElectronPtSource       = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusElectron_Source.Data());
        fHistoTrueClustersGammaPtSource          = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusGamma_Source.Data());

        if (meson.Contains("Pi0")){
            fHistoTrueClustersPrimPi0PtM02                  = (TH2F*) TrueContainer->FindObject(fObjectNameTrueClusPrimMesonM02.Data());
            if (fHistoTrueClustersPrimPi0PtM02) fHistoTrueClustersPrimPi0PtM02->Sumw2();
            fHistoTrueClustersSecPi0PtM02[0]                        = (TH2F*)TrueContainer->FindObject(fObjectNameTrueClusSecMesonFromK0sM02.Data());
            fHistoTrueClustersSecPi0PtM02[1]                        = (TH2F*)TrueContainer->FindObject(fObjectNameTrueClusSecMesonFromLambdaM02.Data());
            fHistoTrueClustersSecPi0PtM02[2]                        = (TH2F*)TrueContainer->FindObject(fObjectNameTrueClusSecMesonFromK0lM02.Data());
            fHistoTrueClustersSecPi0PtM02[3]                        = (TH2F*)TrueContainer->FindObject(fObjectNameTrueClusSecMesonM02.Data());
            for (Int_t j = 0; j < 4; j++){
                if (fHistoTrueClustersSecPi0PtM02[j]) fHistoTrueClustersSecPi0PtM02[j]->Sumw2();
            }
            // remove pi0 from other source (K0s, K0l, Lambda) from secondary yield
            for (Int_t j = 0; j<3; j++){
                if (fHistoTrueClustersSecPi0PtM02[j]) fHistoTrueClustersSecPi0PtM02[3]->Add(fHistoTrueClustersSecPi0PtM02[j],-1);
            }


            FillMCPrimSecM02HistosArray(    fHistoTrueClustersPrimPi0PtM02, fHistoTrueClustersSecPi0PtM02 );
        }


        TList *MCContainer                      = (TList*)HistosGammaConversion->FindObject(Form("%s MC histograms",fCutSelectionRead.Data()));

        fHistoMCMesonGGWithinAccepPt            = (TH1D*)MCContainer->FindObject(fObjectNameMCMesonAcc.Data());
        fHistoMCMesonGGPt                       = (TH1D*)MCContainer->FindObject(fObjectNameMCMeson.Data());
        fHistoMCMultipleMesonGGWithinAccepPt    = (TH1D*)MCContainer->FindObject(fObjectNameMCMultipleMesonAcc.Data());
        fHistoMCMultipleMesonGGPt               = (TH1D*)MCContainer->FindObject(fObjectNameMCMultipleMeson.Data());
        fHistoMCMesonGGPtWOWeights              = (TH1D*)MCContainer->FindObject(fObjectNameMCMesonWOWeights.Data());
        fHistoMCMesonDalitzWithinAccepPt        = (TH1D*)MCContainer->FindObject(fObjectNameMCMesonDalitzAcc.Data());
        fHistoMCMesonDalitzPt                   = (TH1D*)MCContainer->FindObject(fObjectNameMCMesonDalitz.Data());
        fHistoMCMesonDalitzPtWOWeights          = (TH1D*)MCContainer->FindObject(fObjectNameMCMesonDalitzWOWeights.Data());

        // sum all mesons including Dalitz
        fHistoMCMesonWithinAccepPt              = (TH1D*)fHistoMCMesonGGWithinAccepPt->Clone("Meson_inAcc");
        fHistoMCMesonWithinAccepPt->Sumw2();
        fHistoMCMesonWithinAccepPt->Add(fHistoMCMesonDalitzWithinAccepPt);
        if(fHistoMCMultipleMesonGGWithinAccepPt)fHistoMCMesonWithinAccepPt->Add(fHistoMCMultipleMesonGGWithinAccepPt,-1);
        fHistoMCMesonPt                         = (TH1D*)fHistoMCMesonGGPt->Clone("Meson_All");
        fHistoMCMesonPt->Sumw2();
        fHistoMCMesonPt->Add(fHistoMCMesonDalitzPt);
        if(fHistoMCMultipleMesonGGPt)fHistoMCMesonPt->Add(fHistoMCMultipleMesonGGPt,-1);
        fHistoMCMesonPtWOWeights                = (TH1D*)fHistoMCMesonGGPtWOWeights->Clone("Meson_All_WOWeights");
        fHistoMCMesonPtWOWeights->Sumw2();
        fHistoMCMesonPtWOWeights->Add(fHistoMCMesonDalitzPtWOWeights);

        if (meson.Contains("Pi0")){
            cout << "trying to load sec output" << endl;
            fHistoMCSecPi0WithinAccepPtSource   = (TH2D*)MCContainer->FindObject(fObjectNameMCMesonAccSecPi0.Data());
            fHistoMCSecPi0PtSource              = (TH2D*)MCContainer->FindObject(fObjectNameMCMesonSecPi0.Data());
            if (fHistoMCSecPi0WithinAccepPtSource && fHistoMCSecPi0PtSource){
                fHistoMCSecPi0WithinAccepPtSource->Sumw2();
                fHistoMCSecPi0PtSource->Sumw2();
                fNewMCOutput                    = kTRUE;
                cout << "new MC output is being processed" << endl;
            }
        }

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

            if (fHistoTrueClusPartConvMergedFromPi0M02PtBin[iPt]){
                IntegrateHistoM02(fHistoTrueClusPartConvMergedFromPi0M02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TrueMergedPartConvFromPi0Yields[iPt]       = fYields;
                fMesonM02TrueMergedPartConvFromPi0YieldsError[iPt]  = fYieldsError;
            }
            if (fHistoTrueClusPartConvMergedFromEtaM02PtBin[iPt]){
                IntegrateHistoM02(fHistoTrueClusPartConvMergedFromEtaM02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TrueMergedPartConvFromEtaYields[iPt]       = fYields;
                fMesonM02TrueMergedPartConvFromEtaYieldsError[iPt]  = fYieldsError;
            }

            IntegrateHistoM02(fHistoTrueClusPureMergedM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueMergedPureYields[iPt]       = fYields;
            fMesonM02TrueMergedPureYieldsError[iPt]  = fYieldsError;

            if (fHistoTrueClusPureMergedFromPi0M02PtBin[iPt]){
                IntegrateHistoM02(fHistoTrueClusPureMergedFromPi0M02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TrueMergedPureFromPi0Yields[iPt]           = fYields;
                fMesonM02TrueMergedPureFromPi0YieldsError[iPt]      = fYieldsError;
            }
            if (fHistoTrueClusPureMergedFromEtaM02PtBin[iPt]){
                IntegrateHistoM02(fHistoTrueClusPureMergedFromEtaM02PtBin[iPt], fMesonM02IntRange );
                fMesonM02TrueMergedPureFromEtaYields[iPt]           = fYields;
                fMesonM02TrueMergedPureFromEtaYieldsError[iPt]      = fYieldsError;
            }


            IntegrateHistoM02(fHistoTrueClusOneGammaM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueMergedOneGammaYields[iPt]       = fYields;
            fMesonM02TrueMergedOneGammaYieldsError[iPt]  = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusOneGammaFromPi0M02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueMergedOneGammaFromPi0Yields[iPt]       = fYields;
            fMesonM02TrueMergedOneGammaFromPi0YieldsError[iPt]  = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusOneGammaFromEtaM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueMergedOneGammaFromEtaYields[iPt]       = fYields;
            fMesonM02TrueMergedOneGammaFromEtaYieldsError[iPt]  = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusOneElectronM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueMergedOneElectronYields[iPt]       = fYields;
            fMesonM02TrueMergedOneElectronYieldsError[iPt]  = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusOneElectronFromPi0M02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueMergedOneElectronFromPi0Yields[iPt]       = fYields;
            fMesonM02TrueMergedOneElectronFromPi0YieldsError[iPt]  = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusOneElectronFromEtaM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueMergedOneElectronFromEtaYields[iPt]       = fYields;
            fMesonM02TrueMergedOneElectronFromEtaYieldsError[iPt]  = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusPi0M02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TruePi0Yields[iPt]                  = fYields;
            fMesonM02TruePi0YieldsError[iPt]             = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusPi0DCM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TruePi0DCYields[iPt]                = fYields;
            fMesonM02TruePi0DCYieldsError[iPt]           = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusPi0GGM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TruePi0GGYields[iPt]                = fYields;
            fMesonM02TruePi0GGYieldsError[iPt]           = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusPi0DalitzM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TruePi0DalitzYields[iPt]            = fYields;
            fMesonM02TruePi0DalitzYieldsError[iPt]       = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusEtaM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueEtaYields[iPt]                  = fYields;
            fMesonM02TrueEtaYieldsError[iPt]             = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusEtaDCM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueEtaDCYields[iPt]                = fYields;
            fMesonM02TrueEtaDCYieldsError[iPt]           = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusEtaGGM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueEtaGGYields[iPt]                = fYields;
            fMesonM02TrueEtaGGYieldsError[iPt]           = fYieldsError;

            IntegrateHistoM02(fHistoTrueClusEtaDalitzM02PtBin[iPt], fMesonM02IntRange );
            fMesonM02TrueEtaDalitzYields[iPt]            = fYields;
            fMesonM02TrueEtaDalitzYieldsError[iPt]       = fYieldsError;

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

                for (Int_t j = 0; j< 4; j++){
                    IntegrateHistoM02(fHistoTrueClusSecPi0M02PtBin[j][iPt], fMesonM02IntRange );
                    fMesonM02TrueSecPi0Yields[j][iPt]               = fYields;
                    fMesonM02TrueSecPi0YieldsError[j][iPt]          = fYieldsError;
                }
            }
        }
    }

    CreatePtHistos();
    FillPtHistos();

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
    if (fIsMC){
        if (fHistoTrueClustersBGPtSource){
            fHistoTrueClustersBGPtSource->Sumw2();
            FillMCBGSeparated(fHistoTrueClustersBGPtSource);
        }
        if (fHistoTrueClustersElectronPtSource){
            fHistoTrueClustersElectronPtSource->Sumw2();
            FillMCElectronSeparated(fHistoTrueClustersElectronPtSource);
        }
        if (fHistoTrueClustersGammaPtSource){
            fHistoTrueClustersGammaPtSource->Sumw2();
            FillMCGammaSeparated(fHistoTrueClustersGammaPtSource);
        }
        FillHistosArrayMC(fHistoMCMesonWithinAccepPt, fHistoMCMesonPt, fDeltaPt);
        if (meson.Contains("Pi0") && fNewMCOutput){
            FillHistosArrayMCSecAndCalcAcceptance(fHistoMCSecPi0WithinAccepPtSource, fHistoMCSecPi0PtSource);
        }
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
    Double_t minZM02        = 0;
    Double_t maxZM02        = 0;
    Double_t minPtPlotting  = 2.95;
    Double_t maxPtPlotting  = 50;
    if (fEnergyFlag.CompareTo("8TeV") == 0){
        minPtPlotting       = 10;
        maxPtPlotting       = 70;
    }
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
                                "#sigma_{long}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,minPtPlotting, maxPtPlotting,0.8,0.8);
        fHistoClustersMergedPtM02->GetXaxis()->SetMoreLogLabels();
        fHistoClustersMergedPtM02->GetXaxis()->SetLabelOffset(-0.02);
        fHistoClustersMergedPtM02->GetZaxis()->SetLabelOffset(-0.008);
        fHistoClustersMergedPtM02->GetZaxis()->SetLabelSize(0.051);
        fHistoClustersMergedPtM02->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        fHistoClustersMergedPtM02->GetXaxis()->SetTickLength(0.05);
        fHistoClustersMergedPtM02->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.55, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

//         if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin, maxPtPlotting);

        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_PtVsM02_AllAcceptedClusters%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //**************** Plotting 2D M02 vs Pt merged clusters passing meson cuts ****************
    //******************************************************************************************
    if(fHistoClustersMergedPtM02AccMeson){
        fHistoClustersMergedPtM02AccMeson->Sumw2();
        fHistoClustersMergedPtM02AccMeson->Scale(1./fNEvents);
        if (minZM02 == 0) minZM02       = FindSmallestEntryIn2D(fHistoClustersMergedPtM02AccMeson);
        if (maxZM02 == 0) maxZM02       = fHistoClustersMergedPtM02AccMeson->GetMaximum();

        DrawAutoGammaHistoPaper2D(fHistoClustersMergedPtM02AccMeson,
                                " ",
                                "#it{p}_{T} (GeV/#it{c})",
                                "#sigma_{long}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,minPtPlotting, maxPtPlotting,0.8,0.8);
        fHistoClustersMergedPtM02AccMeson->GetXaxis()->SetMoreLogLabels();
        fHistoClustersMergedPtM02AccMeson->GetXaxis()->SetLabelOffset(-0.02);
        fHistoClustersMergedPtM02AccMeson->GetZaxis()->SetLabelSize(0.051);
        fHistoClustersMergedPtM02AccMeson->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        fHistoClustersMergedPtM02AccMeson->GetXaxis()->SetTickLength(0.05);
        fHistoClustersMergedPtM02AccMeson->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.12, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

//         if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin, maxPtPlotting);

        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_PtVsM02_AllAcceptedMesons%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //**************** Plotting 2D M02 vs E merged clusters passing meson cuts ****************
    //******************************************************************************************
    if(fHistoClustersMergedEM02AccMeson){
        fHistoClustersMergedEM02AccMeson->Sumw2();
        fHistoClustersMergedEM02AccMeson->Scale(1./fNEvents);
        if (minZM02 == 0) minZM02     = FindSmallestEntryIn2D(fHistoClustersMergedEM02AccMeson);
        if (maxZM02 == 0) maxZM02     = fHistoClustersMergedEM02AccMeson->GetMaximum();
        DrawAutoGammaHistoPaper2D(fHistoClustersMergedEM02AccMeson,
                                " ",
                                "#it{E} (GeV)",
                                "#sigma_{long}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,2.95, maxPtPlotting,0.8,0.8);
        fHistoClustersMergedEM02AccMeson->GetXaxis()->SetMoreLogLabels();
        fHistoClustersMergedEM02AccMeson->GetXaxis()->SetLabelOffset(-0.02);
        fHistoClustersMergedEM02AccMeson->GetZaxis()->SetLabelSize(0.051);
        fHistoClustersMergedEM02AccMeson->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        fHistoClustersMergedEM02AccMeson->GetXaxis()->SetTickLength(0.05);
        fHistoClustersMergedEM02AccMeson->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.12, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin, maxPtPlotting);

        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_EVsM02_AllAcceptedMesons%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //************** Plotting 2D M02 vs Pt validated merged clusters part conv *****************
    //******************************************************************************************
    if (fHistoTrueClusPartConvMergedPtM02){
        fHistoTrueClusPartConvMergedPtM02->Scale(1./fNEvents);
        DrawAutoGammaHistoPaper2D(fHistoTrueClusPartConvMergedPtM02,
                                " ",
                                "#it{p}_{T} (GeV/#it{c})",
                                "#sigma_{long}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,minPtPlotting, maxPtPlotting,0.8,0.8);
        fHistoTrueClusPartConvMergedPtM02->GetXaxis()->SetMoreLogLabels();
        fHistoTrueClusPartConvMergedPtM02->GetXaxis()->SetLabelOffset(-0.02);
        fHistoTrueClusPartConvMergedPtM02->GetZaxis()->SetLabelOffset(-0.008);
        fHistoTrueClusPartConvMergedPtM02->GetZaxis()->SetLabelSize(0.051);
        fHistoTrueClusPartConvMergedPtM02->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        fHistoTrueClusPartConvMergedPtM02->GetXaxis()->SetTickLength(0.05);
        fHistoTrueClusPartConvMergedPtM02->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.12, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        TLatex *labelM02 = new TLatex(0.11, 0.15, "val. merged clusters, part conv");
        SetStyleTLatex( labelM02, 0.05,4);
        labelM02->Draw();
        if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin, maxPtPlotting);

        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_PtVsM02_TrueMergedPartConv%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //************** Plotting 2D M02 vs Pt validated merged clusters part conv *****************
    //******************************************************************************************
    if (fHistoTrueClusPureMergedPtM02){
        fHistoTrueClusPureMergedPtM02->Scale(1./fNEvents);
        DrawAutoGammaHistoPaper2D(fHistoTrueClusPureMergedPtM02,
                                " ",
                                "#it{p}_{T} (GeV/#it{c})",
                                "#sigma_{long}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,minPtPlotting, maxPtPlotting,0.8,0.8);
        fHistoTrueClusPureMergedPtM02->GetXaxis()->SetMoreLogLabels();
        fHistoTrueClusPureMergedPtM02->GetXaxis()->SetLabelOffset(-0.02);
        fHistoTrueClusPureMergedPtM02->GetZaxis()->SetLabelOffset(-0.008);
        fHistoTrueClusPureMergedPtM02->GetZaxis()->SetLabelSize(0.051);
        fHistoTrueClusPureMergedPtM02->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        fHistoTrueClusPureMergedPtM02->GetXaxis()->SetTickLength(0.05);
        fHistoTrueClusPureMergedPtM02->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.12, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        TLatex *labelM02 = new TLatex(0.11, 0.15, "val. merged clusters, pure");
        SetStyleTLatex( labelM02, 0.05,4);
        labelM02->Draw();
        if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin, maxPtPlotting);

        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_PtVsM02_TrueMergedPure%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //********** Plotting 2D M02 vs Pt validated merged clusters one decay gamma ***************
    //******************************************************************************************
    if (fHistoTrueClusOneGammaPtM02){
        fHistoTrueClusOneGammaPtM02->Scale(1./fNEvents);
        DrawAutoGammaHistoPaper2D(fHistoTrueClusOneGammaPtM02,
                                " ",
                                "#it{p}_{T} (GeV/#it{c})",
                                "#sigma_{long}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,minPtPlotting, maxPtPlotting,0.8,0.8);
        fHistoTrueClusOneGammaPtM02->GetXaxis()->SetMoreLogLabels();
        fHistoTrueClusOneGammaPtM02->GetXaxis()->SetLabelOffset(-0.02);
        fHistoTrueClusOneGammaPtM02->GetZaxis()->SetLabelOffset(-0.008);
        fHistoTrueClusOneGammaPtM02->GetZaxis()->SetLabelSize(0.051);
        fHistoTrueClusOneGammaPtM02->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        fHistoTrueClusOneGammaPtM02->GetXaxis()->SetTickLength(0.05);
        fHistoTrueClusOneGammaPtM02->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.12, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        TLatex *labelM02 = new TLatex(0.11, 0.15, "val. merged clusters, 1 #gamma from decay");
        SetStyleTLatex( labelM02, 0.05,4);
        labelM02->Draw();
        if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin, maxPtPlotting);

        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_PtVsM02_TrueMergedOneDecayGamma%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //******** Plotting 2D M02 vs Pt validated merged clusters one electron from decay *********
    //******************************************************************************************
    if (fHistoTrueClusOneElectronPtM02){
        fHistoTrueClusOneElectronPtM02->Scale(1./fNEvents);
        DrawAutoGammaHistoPaper2D(fHistoTrueClusOneElectronPtM02,
                                " ",
                                "#it{p}_{T} (GeV/#it{c})",
                                "#sigma_{long}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,minPtPlotting, maxPtPlotting,0.8,0.8);
        fHistoTrueClusOneElectronPtM02->GetXaxis()->SetMoreLogLabels();
        fHistoTrueClusOneElectronPtM02->GetXaxis()->SetLabelOffset(-0.02);
        fHistoTrueClusOneElectronPtM02->GetZaxis()->SetLabelOffset(-0.008);
        fHistoTrueClusOneElectronPtM02->GetZaxis()->SetLabelSize(0.051);
        fHistoTrueClusOneElectronPtM02->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        fHistoTrueClusOneElectronPtM02->GetXaxis()->SetTickLength(0.05);
        fHistoTrueClusOneElectronPtM02->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.12, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

        TLatex *labelM02 = new TLatex(0.11, 0.15, "val. merged clusters, 1 e^{+/-} from decay");
        SetStyleTLatex( labelM02, 0.05,4);
        labelM02->Draw();
        if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin, maxPtPlotting);

        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_PtVsM02_TrueMergedOneDecayElectron%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //******************** Plotting 2D M02 vs Pt validated merged pi0s *************************
    //******************************************************************************************
    if (fHistoTrueClustersPi0PtM02 ){
        fHistoTrueClustersPi0PtM02->Sumw2();
        TH2F* dumm2D = (TH2F*)fHistoTrueClustersPi0PtM02->Clone("forplotting");

        dumm2D->Scale(1./fNEvents);
        DrawAutoGammaHistoPaper2D(dumm2D,
                                " ",
                                "#it{p}_{T} (GeV/#it{c})",
                                "#sigma_{long}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,minPtPlotting, maxPtPlotting,0.8,0.8);
        dumm2D->GetXaxis()->SetMoreLogLabels();
        dumm2D->GetXaxis()->SetLabelOffset(-0.02);
        dumm2D->GetZaxis()->SetLabelOffset(-0.008);
        dumm2D->GetZaxis()->SetLabelSize(0.051);
        dumm2D->GetZaxis()->SetRangeUser(minZM02,maxZM02);
        dumm2D->GetXaxis()->SetTickLength(0.05);
        dumm2D->DrawCopy("COLZ");
        PutProcessLabelAndEnergyOnPlot(0.12, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

//         if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin, maxPtPlotting);
        TLatex *labelM02 = new TLatex(0.11, 0.15, "val. merged #pi^{0}");
        SetStyleTLatex( labelM02, 0.05,4);
        labelM02->Draw();


        canvasPtM02->Update();
        canvasPtM02->SaveAs(Form("%s/%s_%s_PtVsM02_TruePi0%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
    }

    //******************************************************************************************
    //******************** Plotting 2D M02 vs Pt validated merged etas *************************
    //******************************************************************************************
    if (fHistoTrueClustersEtaPtM02 ){
        fHistoTrueClustersEtaPtM02->Sumw2();
        TH2F* dumm2D = (TH2F*)fHistoTrueClustersEtaPtM02->Clone("forplotting");
        dumm2D->Scale(1./fNEvents);
        DrawAutoGammaHistoPaper2D(dumm2D,
                                " ",
                                "#it{p}_{T} (GeV/#it{c})",
                                "#sigma_{long}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,minPtPlotting, maxPtPlotting,0.8,0.8);
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

//         if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin, maxPtPlotting);

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
                                "#sigma_{long}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,minPtPlotting, maxPtPlotting,0.8,0.8);
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
//         if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin, maxPtPlotting);

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
                                "#sigma_{long}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,minPtPlotting, maxPtPlotting,0.8,0.8);
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

//         if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin, maxPtPlotting);

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
                                "#sigma_{long}^{2}",
                                0,0,0,
                                1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                1,minPtPlotting, maxPtPlotting,0.8,0.8);
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

//         if (fAdditionalLabels) DrawMergedClusterLambdaCuts(fNLMmin, maxPtPlotting);

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
    //********************** Plotting extended exotics QA if available *************************
    //******************************************************************************************
    if (fEnableExoticsQA){
        if(fHistoExoticsEvsM02){
            canvasPtM02->cd();

            fHistoExoticsEvsM02->Sumw2();
            fHistoExoticsEvsM02->Scale(1./fNEvents);
            if (minZM02 == 0) minZM02     = FindSmallestEntryIn2D(fHistoExoticsEvsM02);
            if (maxZM02 == 0) maxZM02     = fHistoExoticsEvsM02->GetMaximum();
            DrawAutoGammaHistoPaper2D(fHistoExoticsEvsM02,
                                    " ",
                                    "#it{E} (GeV)",
                                    "#sigma_{long}^{2}",
                                    0,0,0,
                                    1,fMesonM02PlotRange[0],fMesonM02PlotRange[1],
                                    1,2.95, maxPtPlotting,0.8,0.8);
            fHistoExoticsEvsM02->GetXaxis()->SetMoreLogLabels();
            fHistoExoticsEvsM02->GetXaxis()->SetLabelOffset(-0.02);
            fHistoExoticsEvsM02->GetZaxis()->SetLabelSize(0.051);
            fHistoExoticsEvsM02->GetZaxis()->SetRangeUser(minZM02,maxZM02);
            fHistoExoticsEvsM02->GetXaxis()->SetTickLength(0.05);
            fHistoExoticsEvsM02->DrawCopy("COLZ");
            PutProcessLabelAndEnergyOnPlot(0.12, 0.97, 0.045, fCollisionSystem.Data(), fNLMString.Data(), fDetectionProcess.Data(), 42, 0.03, "", 1, 1.1);

            canvasPtM02->Update();
            canvasPtM02->SaveAs(Form("%s/%s_%s_EVsM02_RejectedExotics%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
        }

        TCanvas* canvasYieldExo = new TCanvas("canvasYieldExo","",1900,2100);
        DrawGammaCanvasSettings(canvasYieldExo, 0.08, 0.015, 0.02, 0.08);
        canvasYieldExo->SetLogx(0);
        canvasYieldExo->SetLogy(1);
        canvasYieldExo->cd();

        if (fHistoExoticsAboveM02MinCut && fHistoExoticsBelowM02MinBaseCut && fHistoExoticsBelowM02MinCut){
            fHistoExoticsAboveM02MinCut->Sumw2();
            fHistoExoticsAboveM02MinCut->Scale(1./fNEvents);
            fHistoExoticsAboveM02MinCut->Divide(fDeltaPt);
            fHistoExoticsBelowM02MinBaseCut->Sumw2();
            fHistoExoticsBelowM02MinBaseCut->Scale(1./fNEvents);
            fHistoExoticsBelowM02MinBaseCut->Divide(fDeltaPt);
            fHistoExoticsBelowM02MinCut->Sumw2();
            fHistoExoticsBelowM02MinCut->Scale(1./fNEvents);
            fHistoExoticsBelowM02MinCut->Divide(fDeltaPt);

            TH1D* clusterEdummy     = NULL;
            if (fHistoClustersE){
                clusterEdummy   = (TH1D*)fHistoClustersE->Clone("fHistoClusterEDummy");
                clusterEdummy   = (TH1D*)clusterEdummy->Rebin(fNBinsPt,"fHistoClusterEDummy2",fBinsPt);
                clusterEdummy->Sumw2();
                clusterEdummy->Scale(1./fNEvents);
                clusterEdummy->Divide(fDeltaPt);
            }
            TH1D* clusterEmdummy  = NULL;
            if (fHistoClustersMergedE){
                clusterEmdummy   = (TH1D*)fHistoClustersMergedE->Clone("fHistoClusterEmergedDummy");
                clusterEmdummy   = (TH1D*)clusterEmdummy->Rebin(fNBinsPt,"fHistoClusterEmergedDummy2",fBinsPt);
                clusterEmdummy->Sumw2();
                clusterEmdummy->Scale(1./fNEvents);
                clusterEmdummy->Divide(fDeltaPt);
            }
            Double_t maxYExo   = fHistoExoticsAboveM02MinCut->GetMaximum();
            Double_t minYExo   = FindSmallestEntryIn1D(fHistoExoticsAboveM02MinCut);
            if (maxYExo < fHistoExoticsBelowM02MinBaseCut->GetMaximum())
                maxYExo = fHistoExoticsBelowM02MinBaseCut->GetMaximum();
            if (maxYExo < fHistoExoticsBelowM02MinCut->GetMaximum())
                maxYExo = fHistoExoticsBelowM02MinCut->GetMaximum();
            maxYExo = 100*maxYExo;
            if (minYExo > FindSmallestEntryIn1D(fHistoExoticsBelowM02MinBaseCut))
                minYExo = FindSmallestEntryIn1D(fHistoExoticsBelowM02MinBaseCut);
            if (minYExo > FindSmallestEntryIn1D(fHistoExoticsBelowM02MinCut))
                minYExo = FindSmallestEntryIn1D(fHistoExoticsBelowM02MinCut);
            minYExo = 0.1*minYExo;

            DrawAutoGammaMesonHistos(   fHistoExoticsBelowM02MinBaseCut,
                                        "", "E (GeV)", "counts/event",
                                        kFALSE, 2.,1e-6, kFALSE,
                                        kTRUE, minYExo, maxYExo,
                                        kFALSE, 0, 50);
            fHistoExoticsBelowM02MinBaseCut->GetYaxis()->SetTitleOffset(1);
            DrawGammaSetMarker(fHistoExoticsBelowM02MinBaseCut, 21, 1.5, kGreen+2, kGreen+2);
            fHistoExoticsBelowM02MinBaseCut->Draw("e,p");
            DrawGammaSetMarker(fHistoExoticsAboveM02MinCut, 25, 1.5, kBlue+1, kBlue+1);
            fHistoExoticsAboveM02MinCut->Draw("same,e,p");
            DrawGammaSetMarker(fHistoExoticsBelowM02MinCut, 24, 1.5, kRed+1, kRed+1);
            fHistoExoticsBelowM02MinCut->Draw("same,e,p");

            TLegend* legendExotics = GetAndSetLegend2(0.65,0.8-0.035*3*1.05,0.95,0.8,0.035,1,"",42,0.1);
            legendExotics->AddEntry(fHistoExoticsBelowM02MinBaseCut,"0 < #sigma_{long}^{2} < 0.1","p");
            legendExotics->AddEntry(fHistoExoticsBelowM02MinCut,"0.1 < #sigma_{long}^{2} < 0.27","p");
            legendExotics->AddEntry(fHistoExoticsAboveM02MinCut,"0.27 < #sigma_{long}^{2} < #infty","p");
            legendExotics->Draw();


            PutProcessLabelAndEnergyOnPlot(0.65, 0.95, 0.035, fCollisionSystem.Data(), fDetectionProcess.Data(), "rejected exotics", 42, 0.03, "", 1, 1.1);

            canvasYieldExo->Update();
            canvasYieldExo->SaveAs(Form("%s/%s_%s_ExoticsRej_YieldPerEvent%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));


            maxYExo = maxYExo/100;
            if (clusterEdummy){
                if (maxYExo < clusterEdummy->GetMaximum())
                    maxYExo = clusterEdummy->GetMaximum();
            }
            if (clusterEmdummy){
                if (maxYExo < clusterEmdummy->GetMaximum())
                    maxYExo = clusterEmdummy->GetMaximum();
            }
            maxYExo = maxYExo*100;
            minYExo = minYExo/0.1;
            if (clusterEdummy){
                if (minYExo > FindSmallestEntryIn1D(clusterEdummy))
                    minYExo = FindSmallestEntryIn1D(clusterEdummy);
            }
            if (clusterEmdummy){
                if (minYExo > FindSmallestEntryIn1D(clusterEmdummy))
                    minYExo = FindSmallestEntryIn1D(clusterEmdummy);
            }
            minYExo = 0.1*minYExo;

            DrawAutoGammaMesonHistos(   fHistoExoticsBelowM02MinBaseCut,
                                        "", "E (GeV)", "counts/event",
                                        kFALSE, 2.,1e-6, kFALSE,
                                        kTRUE, minYExo, maxYExo,
                                        kFALSE, 0, 50);
            fHistoExoticsBelowM02MinBaseCut->GetYaxis()->SetTitleOffset(1);
            DrawGammaSetMarker(fHistoExoticsBelowM02MinBaseCut, 21, 1.5, kGreen+2, kGreen+2);
            fHistoExoticsBelowM02MinBaseCut->Draw("e,p");
            DrawGammaSetMarker(fHistoExoticsAboveM02MinCut, 25, 1.5, kBlue+1, kBlue+1);
            fHistoExoticsAboveM02MinCut->Draw("same,e,p");
            DrawGammaSetMarker(fHistoExoticsBelowM02MinCut, 24, 1.5, kRed+1, kRed+1);
            fHistoExoticsBelowM02MinCut->Draw("same,e,p");
            if (clusterEdummy){
                DrawGammaSetMarker(clusterEdummy, 20, 2, kBlack, kBlack);
                clusterEdummy->Draw("same,e,p");
            }
            if (clusterEmdummy){
                DrawGammaSetMarker(clusterEmdummy, 30, 3.5, kBlack, kBlack);
                clusterEmdummy->Draw("same,e,p");
            }

            legendExotics->Draw();

            TLegend* legendExotics2 = GetAndSetLegend2(0.13,0.94-0.035*3*1.05,0.43,0.94,0.035,1,"Good clusters",42,0.1);
            if (clusterEdummy)legendExotics2->AddEntry(clusterEdummy,"#sigma_{long}^{2} > 0.1","p");
            if (clusterEmdummy)legendExotics2->AddEntry(clusterEmdummy,"merged cl. #sigma_{long}^{2} > 0.27","p");
            legendExotics2->Draw();

            PutProcessLabelAndEnergyOnPlot(0.65, 0.95, 0.035, fCollisionSystem.Data(), fDetectionProcess.Data(), "rejected exotics", 42, 0.03, "", 1, 1.1);

            canvasYieldExo->Update();
            canvasYieldExo->SaveAs(Form("%s/%s_%s_ExoticsRejWithGood_YieldPerEvent%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));

        }
        delete canvasYieldExo;



        TCanvas* canvasNCellsVsPt = new TCanvas("canvasNCellsVsPt","",1900,1500);
        DrawGammaCanvasSettings(canvasNCellsVsPt, 0.09, 0.13, 0.02, 0.1);
        canvasNCellsVsPt->SetLogx(1);
        canvasNCellsVsPt->SetLogy(0);
        canvasNCellsVsPt->SetLogz(1);
        canvasNCellsVsPt->cd();
        Double_t minZNcells     = 0;
        Double_t maxZNcells     = 0;

        if(fHistoExoticsEvsNCells){
            fHistoExoticsEvsNCells->Sumw2();
            fHistoExoticsEvsNCells->Scale(1./fNEvents);
            maxZNcells     = fHistoExoticsEvsNCells->GetMaximum();
            minZNcells     = FindSmallestEntryIn2D(fHistoExoticsEvsNCells);

            DrawAutoGammaHistoPaper2D(fHistoExoticsEvsNCells,
                                    " ",
                                    "#it{E} (GeV)",
                                    "#it{N}_{cells in cluster}",
                                    0,0,0,
                                    1,0.7,50.0,
                                    1,-0.5,99.5,
                                    0.9,0.8);
            fHistoExoticsEvsNCells->GetYaxis()->SetMoreLogLabels();
            fHistoExoticsEvsNCells->GetYaxis()->SetLabelOffset(-0.001);
            fHistoExoticsEvsNCells->GetZaxis()->SetLabelOffset(-0.008);
            fHistoExoticsEvsNCells->GetZaxis()->SetLabelSize(0.051);
            fHistoExoticsEvsNCells->GetXaxis()->SetTickLength(0.03);
            fHistoExoticsEvsNCells->GetZaxis()->SetRangeUser(minZNcells,maxZNcells);
            fHistoExoticsEvsNCells->DrawCopy("COLZ");
            PutProcessLabelAndEnergyOnPlot(0.52, 0.96, 0.045, fCollisionSystem.Data(), fDetectionProcess.Data(), "", 42, 0.03, "", 1, 1.1);

            canvasNCellsVsPt->Update();
            canvasNCellsVsPt->SaveAs(Form("%s/%s_%s_ExoticsRej_ENCells%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
        }
        delete canvasNCellsVsPt;

        TCanvas* canvasEvsEstar = new TCanvas("canvasEvsEstar","",1900,1500);
        DrawGammaCanvasSettings(canvasEvsEstar, 0.09, 0.13, 0.02, 0.1);
        canvasEvsEstar->SetLogx(1);
        canvasEvsEstar->SetLogy(1);
        canvasEvsEstar->SetLogz(1);
        canvasEvsEstar->cd();
        Double_t minZEstar     = 0;
        Double_t maxZEstar     = 0;

        if(fHistoExoticsEvsEstar){
            fHistoExoticsEvsEstar->Sumw2();
            fHistoExoticsEvsEstar->Scale(1./fNEvents);
            maxZEstar     = fHistoExoticsEvsEstar->GetMaximum();
            minZEstar     = FindSmallestEntryIn2D(fHistoExoticsEvsEstar);

            DrawAutoGammaHistoPaper2D(fHistoExoticsEvsEstar,
                                    " ",
                                    "#it{E} (GeV)",
                                    "#it{E}_{star} (GeV)",
                                    0,0,0,
                                    1,0.7,50.0,
                                    1,0.7,50.0,
                                    0.9,0.8);
            fHistoExoticsEvsEstar->GetYaxis()->SetMoreLogLabels();
            fHistoExoticsEvsEstar->GetYaxis()->SetLabelOffset(-0.001);
            fHistoExoticsEvsEstar->GetZaxis()->SetLabelOffset(-0.008);
            fHistoExoticsEvsEstar->GetZaxis()->SetLabelSize(0.051);
            fHistoExoticsEvsEstar->GetXaxis()->SetTickLength(0.03);
            fHistoExoticsEvsEstar->GetZaxis()->SetRangeUser(minZEstar,maxZEstar);
            fHistoExoticsEvsEstar->DrawCopy("COLZ");
            PutProcessLabelAndEnergyOnPlot(0.52, 0.96, 0.045, fCollisionSystem.Data(), fDetectionProcess.Data(), "", 42, 0.03, "", 1, 1.1);

            canvasEvsEstar->Update();
            canvasEvsEstar->SaveAs(Form("%s/%s_%s_ExoticsRej_EvsEstar%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
        }
        delete canvasEvsEstar;

        TCanvas* canvasExoticsEtaPhi = new TCanvas("canvasExoticsEtaPhi","",1900,1500);
        DrawGammaCanvasSettings(canvasExoticsEtaPhi, 0.09, 0.13, 0.02, 0.1);
        canvasExoticsEtaPhi->SetLogx(0);
        canvasExoticsEtaPhi->SetLogy(0);
        canvasExoticsEtaPhi->SetLogz(1);
        canvasExoticsEtaPhi->cd();
        Double_t minZExoEtaPhi     = 0;
        Double_t maxZExoEtaPhi     = 0;

        if(fHistoExoticsEtaPhi){
            fHistoExoticsEtaPhi->Sumw2();
            fHistoExoticsEtaPhi->Scale(1./fNEvents);
            maxZExoEtaPhi     = fHistoExoticsEtaPhi->GetMaximum();
            minZExoEtaPhi     = FindSmallestEntryIn2D(fHistoExoticsEtaPhi);

            DrawAutoGammaHistoPaper2D(fHistoExoticsEtaPhi,
                                    " ",
                                    "#it{E} (GeV)",
                                    "#it{E}_{star} (GeV)",
                                    0,0,0,
                                    1,0.7,50.0,
                                    1,0.7,50.0,
                                    0.9,0.8);
            fHistoExoticsEtaPhi->GetYaxis()->SetMoreLogLabels();
            fHistoExoticsEtaPhi->GetYaxis()->SetLabelOffset(-0.001);
            fHistoExoticsEtaPhi->GetZaxis()->SetLabelOffset(-0.008);
            fHistoExoticsEtaPhi->GetZaxis()->SetLabelSize(0.051);
            fHistoExoticsEtaPhi->GetXaxis()->SetTickLength(0.03);
            fHistoExoticsEtaPhi->GetZaxis()->SetRangeUser(minZExoEtaPhi,maxZExoEtaPhi);
            fHistoExoticsEtaPhi->DrawCopy("COLZ");
            PutProcessLabelAndEnergyOnPlot(0.52, 0.96, 0.045, fCollisionSystem.Data(), fDetectionProcess.Data(), "", 42, 0.03, "", 1, 1.1);

            canvasExoticsEtaPhi->Update();
            canvasExoticsEtaPhi->SaveAs(Form("%s/%s_%s_ExoticsRej_EtaPhi%s.%s", outputDir.Data(), fPrefix.Data(), fPrefix2.Data(), fAdditionalName.Data(), suffix.Data()));
        }
        delete canvasExoticsEtaPhi;
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


    nameMeson   = Form("%s_MesonM02%s", plotPrefix.Data(), plotSuffix.Data());
    cout << nameMeson.Data() << endl;
    PlotM02MergedInPtBins(  fHistoM02PtBin, nameMeson, nameCanvas, namePad, fMesonM02PlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt,
                            fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem);

    if (fEnableExoticsQA){
        nameMeson   = Form("%s_RejectedExoticsM02%s", plotPrefix.Data(), plotSuffix.Data());
        cout << nameMeson.Data() << endl;
        PlotM02ExoticsInEBins(  fHistoExoticsM02EBin, nameMeson, nameCanvas, namePad, fMesonM02PlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt,
                                fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem);
    }

    if (fIsMC){
        nameMeson       = Form("%s_MesonM02Split%s", plotPrefix.Data(), plotSuffix.Data());
        cout << nameMeson.Data() << endl;
        PlotM02MergedMCInPtBins(    fHistoM02PtBin, fHistoTrueClusMergedM02PtBin, fHistoTrueClusPartConvMergedM02PtBin,
                                        fHistoTrueClusBGM02PtBin, fHistoTrueClusGammaM02PtBin, fHistoTrueClusElectronM02PtBin,
                                        nameMeson, nameCanvas, namePad, fMesonM02PlotRange,
                                        fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem);

        nameMeson       = Form("%s_MesonM02ValidatedDisentangled%s", plotPrefix.Data(), plotSuffix.Data());

        PlotM02MergedTrueInPtBins(  fHistoTrueClusMergedM02PtBin, NULL, fHistoTrueClusPi0M02PtBin,
                                        NULL, fHistoTrueClusEtaM02PtBin, NULL,
                                        nameMeson, nameCanvas, namePad, fMesonM02PlotRange,
                                        fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcessPtBins, fCollisionSystem);
        if (meson.Contains("Pi0")){
            nameMeson       = Form("%s_MesonPi0ValidatedDisentangled%s", plotPrefix.Data(), plotSuffix.Data());
            cout << nameMeson.Data() << endl;
            PlotM02MergedTruePrimSecInPtBins(   fHistoTrueClusPi0M02PtBin, fHistoTrueClusPrimPi0M02PtBin, fHistoTrueClusSecPi0M02PtBin[3],
                                                fHistoTrueClusSecPi0M02PtBin[0], fHistoTrueClusSecPi0M02PtBin[1],
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
        fHistoTrueEffiMerged        = CalculateMesonEfficiency(fHistoTrueYieldMergedM02, NULL, fNameHistoEffi);

        fNameHistoEffi              = "TrueMesonEffiPrimPt";
        cout << fNameHistoEffi.Data() << endl;
        if (meson.Contains("Pi0")){
            fHistoTrueEffiPrimMeson = CalculateMesonEfficiency(fHistoTrueYieldPrimPi0M02, NULL, fNameHistoEffi);
        } else {
            fHistoTrueEffiPrimMeson = CalculateMesonEfficiency(fHistoTrueYieldEtaM02, NULL, fNameHistoEffi);
        }

        TString fNameHistoPur       = "TrueMesonPurityMergedPt";
        cout << fNameHistoPur.Data() << endl;
        fHistoTruePurityMerged      = CalculatePurity(fHistoYieldMesonM02, fHistoTrueYieldMergedM02, fNameHistoPur);

        fNameHistoPur               = "TruePi0PurityMergedPt";
        cout << fNameHistoPur.Data() << endl;
        fHistoTruePi0PurityMerged   = CalculatePurity(fHistoYieldMesonM02, fHistoTrueYieldPi0M02, fNameHistoPur);

        fNameHistoPur               = "TrueEtaPurityMergedPt";
        cout << fNameHistoPur.Data() << endl;
        fHistoTrueEtaPurityMerged   = CalculatePurity(fHistoYieldMesonM02, fHistoTrueYieldEtaM02, fNameHistoPur);

        // Calculation of secondary fractions
        if (meson.Contains("Pi0")){
            for (Int_t j = 0; j < 4; j++){
                TString fNameHistoFrac      = Form("TrueSecFracFrom%s",nameSecondaries[j].Data());
                fHistoTruePi0SecFrac[j]     = CalculateSecondaryFractions(fHistoTrueYieldPi0M02, fHistoTrueYieldSecPi0M02[j], fNameHistoFrac);
            }
            if (fNewMCOutput){
                for (Int_t j = 0; j < 4; j++){
                    fNameHistoEffi              = Form("TrueMesonEffiSecFrom%sPt",nameSecondaries[j].Data());
                    fHistoTrueEffiSecPi0[j]     = CalculateMesonEfficiencySec(fHistoTrueYieldSecPi0M02[j], fHistoMCSecPi0PtWAccReb[j], fNameHistoEffi);
                }
            }
        }

        CreateRatioHistos();

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


//**********************************************************************************
//******************* return minimum for 2 D histo  ********************************
//**********************************************************************************
Double_t FindSmallestEntryIn1D(TH1* histo){
    Double_t minimum = 1;
    for (Int_t i = 1; i<histo->GetNbinsX(); i++){
        if (histo->GetBinContent(i) < minimum && histo->GetBinContent(i) > 0){
            minimum = histo->GetBinContent(i);
        }
    }
    return minimum;
}


//****************************************************************************
//****** Initialization of arrays and variables, histograms for analysis *****
//****** depending on mesonType, number of bins, mode, energy, centrality ****
//****************************************************************************
void Initialize(TString setPi0, Int_t numberOfBins, Int_t triggerSet){

    InitializeBinning(setPi0, numberOfBins, fEnergyFlag, "", fMode, fEventCutSelection, fClusterMergedCutSelection, triggerSet);
    if (setPi0.CompareTo("Pi0") == 0 || setPi0.CompareTo("Pi0EtaBinning") == 0){
        fMesonId                        = 111;
        fMesonMassPlotRange             = new Double_t[2];
        fMesonMassPlotRange[0]          = 0.;
        fMesonMassPlotRange[1]          = 0.5;
        fMesonMassIntRange              = new Double_t[2];
        fMesonMassIntRange[0]           = 0.;
        fMesonMassIntRange[1]           = 0.5;

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
    fHistoExoticsM02EBin                = new TH1D*[fNBinsPt];
    fMesonM02Yields                     = new Double_t[fNBinsPt];
    fMesonM02YieldsError                = new Double_t[fNBinsPt];

    if (fIsMC){
        fHistoTrueClusMergedM02PtBin                = new TH1D*[fNBinsPt];
        fHistoTrueClusPi0M02PtBin                   = new TH1D*[fNBinsPt];
        fHistoTrueClusPi0DCM02PtBin                 = new TH1D*[fNBinsPt];
        fHistoTrueClusPi0GGM02PtBin                 = new TH1D*[fNBinsPt];
        fHistoTrueClusPi0DalitzM02PtBin             = new TH1D*[fNBinsPt];
        fHistoTrueClusEtaM02PtBin                   = new TH1D*[fNBinsPt];
        fHistoTrueClusEtaDCM02PtBin                 = new TH1D*[fNBinsPt];
        fHistoTrueClusEtaGGM02PtBin                 = new TH1D*[fNBinsPt];
        fHistoTrueClusEtaDalitzM02PtBin             = new TH1D*[fNBinsPt];
        fHistoTrueClusGammaM02PtBin                 = new TH1D*[fNBinsPt];
        fHistoTrueClusElectronM02PtBin              = new TH1D*[fNBinsPt];
        fHistoTrueClusBGM02PtBin                    = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvMergedM02PtBin        = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvMergedFromPi0M02PtBin = new TH1D*[fNBinsPt];
        fHistoTrueClusPartConvMergedFromEtaM02PtBin = new TH1D*[fNBinsPt];
        fHistoTrueClusPureMergedM02PtBin            = new TH1D*[fNBinsPt];
        fHistoTrueClusPureMergedFromPi0M02PtBin     = new TH1D*[fNBinsPt];
        fHistoTrueClusPureMergedFromEtaM02PtBin     = new TH1D*[fNBinsPt];
        fHistoTrueClusOneGammaM02PtBin              = new TH1D*[fNBinsPt];
        fHistoTrueClusOneGammaFromPi0M02PtBin       = new TH1D*[fNBinsPt];
        fHistoTrueClusOneGammaFromEtaM02PtBin       = new TH1D*[fNBinsPt];
        fHistoTrueClusOneElectronM02PtBin           = new TH1D*[fNBinsPt];
        fHistoTrueClusOneElectronFromPi0M02PtBin    = new TH1D*[fNBinsPt];
        fHistoTrueClusOneElectronFromEtaM02PtBin    = new TH1D*[fNBinsPt];
        fHistoTrueClusPrimPi0M02PtBin               = new TH1D*[fNBinsPt];
        for (Int_t j = 0; j < 4; j++){
            fHistoTrueClusSecPi0M02PtBin[j]         = new TH1D*[fNBinsPt];
        }

        fMesonM02TrueMergedYields                   = new Double_t[fNBinsPt];
        fMesonM02TruePi0Yields                      = new Double_t[fNBinsPt];
        fMesonM02TruePi0DCYields                    = new Double_t[fNBinsPt];
        fMesonM02TruePi0GGYields                    = new Double_t[fNBinsPt];
        fMesonM02TruePi0DalitzYields                = new Double_t[fNBinsPt];
        fMesonM02TrueEtaYields                      = new Double_t[fNBinsPt];
        fMesonM02TrueEtaDCYields                    = new Double_t[fNBinsPt];
        fMesonM02TrueEtaGGYields                    = new Double_t[fNBinsPt];
        fMesonM02TrueEtaDalitzYields                = new Double_t[fNBinsPt];
        fMesonM02TrueGammaYields                    = new Double_t[fNBinsPt];
        fMesonM02TrueElectronYields                 = new Double_t[fNBinsPt];
        fMesonM02TrueBGYields                       = new Double_t[fNBinsPt];
        fMesonM02TrueMergedYieldsError              = new Double_t[fNBinsPt];
        fMesonM02TruePi0YieldsError                 = new Double_t[fNBinsPt];
        fMesonM02TruePi0DCYieldsError               = new Double_t[fNBinsPt];
        fMesonM02TruePi0GGYieldsError               = new Double_t[fNBinsPt];
        fMesonM02TruePi0DalitzYieldsError           = new Double_t[fNBinsPt];
        fMesonM02TrueEtaYieldsError                 = new Double_t[fNBinsPt];
        fMesonM02TrueEtaDCYieldsError               = new Double_t[fNBinsPt];
        fMesonM02TrueEtaGGYieldsError               = new Double_t[fNBinsPt];
        fMesonM02TrueEtaDalitzYieldsError           = new Double_t[fNBinsPt];
        fMesonM02TrueGammaYieldsError               = new Double_t[fNBinsPt];
        fMesonM02TrueElectronYieldsError            = new Double_t[fNBinsPt];
        fMesonM02TrueBGYieldsError                  = new Double_t[fNBinsPt];
        fMesonM02TrueMergedPartConvYields           = new Double_t[fNBinsPt];
        fMesonM02TrueMergedPartConvYieldsError      = new Double_t[fNBinsPt];
        fMesonM02TrueMergedPartConvFromPi0Yields        = new Double_t[fNBinsPt];
        fMesonM02TrueMergedPartConvFromPi0YieldsError   = new Double_t[fNBinsPt];
        fMesonM02TrueMergedPartConvFromEtaYields        = new Double_t[fNBinsPt];
        fMesonM02TrueMergedPartConvFromEtaYieldsError   = new Double_t[fNBinsPt];
        fMesonM02TrueMergedPureYields               = new Double_t[fNBinsPt];
        fMesonM02TrueMergedPureYieldsError          = new Double_t[fNBinsPt];
        fMesonM02TrueMergedPureFromPi0Yields        = new Double_t[fNBinsPt];
        fMesonM02TrueMergedPureFromPi0YieldsError   = new Double_t[fNBinsPt];
        fMesonM02TrueMergedPureFromEtaYields        = new Double_t[fNBinsPt];
        fMesonM02TrueMergedPureFromEtaYieldsError   = new Double_t[fNBinsPt];
        fMesonM02TrueMergedOneGammaYields           = new Double_t[fNBinsPt];
        fMesonM02TrueMergedOneGammaYieldsError      = new Double_t[fNBinsPt];
        fMesonM02TrueMergedOneGammaFromPi0Yields            = new Double_t[fNBinsPt];
        fMesonM02TrueMergedOneGammaFromPi0YieldsError       = new Double_t[fNBinsPt];
        fMesonM02TrueMergedOneGammaFromEtaYields            = new Double_t[fNBinsPt];
        fMesonM02TrueMergedOneGammaFromEtaYieldsError       = new Double_t[fNBinsPt];
        fMesonM02TrueMergedOneElectronYields        = new Double_t[fNBinsPt];
        fMesonM02TrueMergedOneElectronYieldsError   = new Double_t[fNBinsPt];
        fMesonM02TrueMergedOneElectronFromPi0Yields         = new Double_t[fNBinsPt];
        fMesonM02TrueMergedOneElectronFromPi0YieldsError    = new Double_t[fNBinsPt];
        fMesonM02TrueMergedOneElectronFromEtaYields         = new Double_t[fNBinsPt];
        fMesonM02TrueMergedOneElectronFromEtaYieldsError    = new Double_t[fNBinsPt];

        fMesonM02TruePrimPi0Yields                  = new Double_t[fNBinsPt];
        fMesonM02TruePrimPi0YieldsError             = new Double_t[fNBinsPt];
        for (Int_t j = 0; j< 4; j++){
            fMesonM02TrueSecPi0Yields[j]            = new Double_t[fNBinsPt];
            fMesonM02TrueSecPi0YieldsError[j]       = new Double_t[fNBinsPt];
        }

    }

    for(Int_t i = 0;i<fNBinsPt; i++){
        fHistoInvMassPtBin[i]               = NULL;
        fHistoM02PtBin[i]                   = NULL;
        fHistoExoticsM02EBin[i]             = NULL;
        fMesonM02Yields[i]                  = 0.;
        fMesonM02YieldsError[i]             = 0.;
        if (fIsMC){
            fHistoTrueClusMergedM02PtBin[i]             = NULL;
            fHistoTrueClusPi0M02PtBin[i]                = NULL;
            fHistoTrueClusPi0DCM02PtBin[i]              = NULL;
            fHistoTrueClusPi0GGM02PtBin[i]              = NULL;
            fHistoTrueClusPi0DalitzM02PtBin[i]          = NULL;
            fHistoTrueClusEtaM02PtBin[i]                = NULL;
            fHistoTrueClusEtaDCM02PtBin[i]              = NULL;
            fHistoTrueClusEtaGGM02PtBin[i]              = NULL;
            fHistoTrueClusEtaDalitzM02PtBin[i]          = NULL;
            fHistoTrueClusGammaM02PtBin[i]              = NULL;
            fHistoTrueClusElectronM02PtBin[i]           = NULL;
            fHistoTrueClusBGM02PtBin[i]                 = NULL;
            fHistoTrueClusPartConvMergedM02PtBin[i]     = NULL;
            fHistoTrueClusPartConvMergedFromPi0M02PtBin[i]  = NULL;
            fHistoTrueClusPartConvMergedFromEtaM02PtBin[i]  = NULL;
            fHistoTrueClusPureMergedM02PtBin[i]         = NULL;
            fHistoTrueClusPureMergedFromPi0M02PtBin[i]  = NULL;
            fHistoTrueClusPureMergedFromEtaM02PtBin[i]  = NULL;
            fHistoTrueClusOneGammaM02PtBin[i]           = NULL;
            fHistoTrueClusOneGammaFromPi0M02PtBin[i]    = NULL;
            fHistoTrueClusOneGammaFromEtaM02PtBin[i]    = NULL;
            fHistoTrueClusOneElectronM02PtBin[i]        = NULL;
            fHistoTrueClusOneElectronFromPi0M02PtBin[i] = NULL;
            fHistoTrueClusOneElectronFromEtaM02PtBin[i] = NULL;
            fHistoTrueClusPrimPi0M02PtBin[i]            = NULL;
            for (Int_t j = 0; j < 4; j++){
                fHistoTrueClusSecPi0M02PtBin[j][i]      = NULL;
            }

            fMesonM02TrueMergedYields[i]                = 0.;
            fMesonM02TruePi0Yields[i]                   = 0.;
            fMesonM02TruePi0DCYields[i]                 = 0.;
            fMesonM02TruePi0GGYields[i]                 = 0.;
            fMesonM02TruePi0DalitzYields[i]             = 0.;
            fMesonM02TrueEtaYields[i]                   = 0.;
            fMesonM02TrueEtaDCYields[i]                 = 0.;
            fMesonM02TrueEtaGGYields[i]                 = 0.;
            fMesonM02TrueEtaDalitzYields[i]             = 0.;
            fMesonM02TrueGammaYields[i]                 = 0.;
            fMesonM02TrueElectronYields[i]              = 0.;
            fMesonM02TrueBGYields[i]                    = 0.;
            fMesonM02TrueMergedYieldsError[i]           = 0.;
            fMesonM02TruePi0YieldsError[i]              = 0.;
            fMesonM02TruePi0DCYieldsError[i]            = 0.;
            fMesonM02TruePi0GGYieldsError[i]            = 0.;
            fMesonM02TruePi0DalitzYieldsError[i]        = 0.;
            fMesonM02TrueEtaYieldsError[i]              = 0.;
            fMesonM02TrueEtaDCYieldsError[i]            = 0.;
            fMesonM02TrueEtaGGYieldsError[i]            = 0.;
            fMesonM02TrueEtaDalitzYieldsError[i]        = 0.;
            fMesonM02TrueGammaYieldsError[i]            = 0.;
            fMesonM02TrueElectronYieldsError[i]         = 0.;
            fMesonM02TrueBGYieldsError[i]               = 0.;
            fMesonM02TrueMergedPartConvYields[i]        = 0.;
            fMesonM02TrueMergedPartConvYieldsError[i]   = 0.;
            fMesonM02TrueMergedPartConvFromPi0Yields[i]         = 0.;
            fMesonM02TrueMergedPartConvFromPi0YieldsError[i]    = 0.;
            fMesonM02TrueMergedPartConvFromEtaYields[i]         = 0.;
            fMesonM02TrueMergedPartConvFromEtaYieldsError[i]    = 0.;
            fMesonM02TrueMergedPureYields[i]            = 0.;
            fMesonM02TrueMergedPureYieldsError[i]       = 0.;
            fMesonM02TrueMergedPureFromPi0Yields[i]             = 0.;
            fMesonM02TrueMergedPureFromPi0YieldsError[i]        = 0.;
            fMesonM02TrueMergedPureFromEtaYields[i]             = 0.;
            fMesonM02TrueMergedPureFromEtaYieldsError[i]        = 0.;
            fMesonM02TrueMergedOneGammaYields[i]        = 0.;
            fMesonM02TrueMergedOneGammaYieldsError[i]   = 0.;
            fMesonM02TrueMergedOneGammaFromPi0Yields[i]         = 0.;
            fMesonM02TrueMergedOneGammaFromPi0YieldsError[i]    = 0.;
            fMesonM02TrueMergedOneGammaFromEtaYields[i]         = 0.;
            fMesonM02TrueMergedOneGammaFromEtaYieldsError[i]    = 0.;
            fMesonM02TrueMergedOneElectronYields[i]             = 0.;
            fMesonM02TrueMergedOneElectronYieldsError[i]        = 0.;
            fMesonM02TrueMergedOneElectronFromPi0Yields[i]      = 0.;
            fMesonM02TrueMergedOneElectronFromPi0YieldsError[i] = 0.;
            fMesonM02TrueMergedOneElectronFromEtaYields[i]      = 0.;
            fMesonM02TrueMergedOneElectronFromEtaYieldsError[i] = 0.;
            fMesonM02TruePrimPi0Yields[i]               = 0.;
            fMesonM02TruePrimPi0YieldsError[i]          = 0.;
            for (Int_t j = 0; j<4; j++){
                fMesonM02TrueSecPi0Yields[j][i]         = 0.;
                fMesonM02TrueSecPi0YieldsError[j][i]    = 0.;
            }
        }
    }
}

//****************************************************************************
//************** Initializiation of MC histogram names according *************
//****************************************************************************
void SetCorrectMCHistogrammNames(TString mesonType){

    fObjectNameTrueMergedM02                    = "ESD_TrueClusMerged_Pt_M02";
    fObjectNameTrueFromPi0M02                   = "ESD_TrueClusFromPi0_Pt_M02";
    fObjectNameTrueFromPi0DCM02                 = "ESD_TrueDoubleCountPi0_Pt_M02";
    fObjectNameTrueFromPi0GGM02                 = "ESD_TrueClusFromPi0GG_Pt_M02";
    fObjectNameTrueFromPi0DalitzM02             = "ESD_TrueClusFromPi0Dalitz_Pt_M02";
    fObjectNameTrueFromEtaM02                   = "ESD_TrueClusFromEta_Pt_M02";
    fObjectNameTrueFromEtaDCM02                 = "ESD_TrueDoubleCountEta_Pt_M02";
    fObjectNameTrueFromEtaGGM02                 = "ESD_TrueClusFromEtaGG_Pt_M02";
    fObjectNameTrueFromEtaDalitzM02             = "ESD_TrueClusFromEtaDalitz_Pt_M02";
    fObjectNameTrueMergedPureM02                = "ESD_TrueClusMergedPure_Pt_M02";
    fObjectNameTrueMergedPureFromPi0M02         = "ESD_TrueClusMergedPureFromPi0_Pt_M02";
    fObjectNameTrueMergedPureFromEtaM02         = "ESD_TrueClusMergedPureFromEta_Pt_M02";
    fObjectNameTrueMergedPartConvM02            = "ESD_TrueClusMergedPartConv_Pt_M02";
    fObjectNameTrueMergedPartConvFromPi0M02     = "ESD_TrueClusMergedPartConvFromPi0_Pt_M02";
    fObjectNameTrueMergedPartConvFromEtaM02     = "ESD_TrueClusMergedPartConvFromEta_Pt_M02";
    fObjectNameTrueMergedOneGammaM02            = "ESD_TrueClusGamma_FromMeson_Pt_M02";
    fObjectNameTrueMergedOneGammaFromPi0M02     = "ESD_TrueClusGamma_FromPi0_Pt_M02";
    fObjectNameTrueMergedOneGammaFromEtaM02     = "ESD_TrueClusGamma_FromEta_Pt_M02";
    fObjectNameTrueMergedOneElectronM02         = "ESD_TrueClusElectron_FromMeson_Pt_M02";
    fObjectNameTrueMergedOneElectronFromPi0M02  = "ESD_TrueClusElectron_FromPi0_Pt_M02";
    fObjectNameTrueMergedOneElectronFromEtaM02  = "ESD_TrueClusElectron_FromEta_Pt_M02";
    fObjectNameTrueClusBGM02                    = "ESD_TrueClusBG_Pt_M02";
    fObjectNameTrueClusGammaM02                 = "ESD_TrueClusGamma_Pt_M02";
    fObjectNameTrueClusElectronM02              = "ESD_TrueClusElectron_Pt_M02";
    fObjectNameTrueClusBG_Source                = "ESD_TrueClusBG_Pt_Source";
    fObjectNameTrueClusGamma_Source             = "ESD_TrueClusGamma_Pt_Source";
    fObjectNameTrueClusElectron_Source          = "ESD_TrueClusElectron_Pt_Source";
    if (mesonType.Contains("Pi0")){
        fObjectNameMCMesonAcc                       = "MC_Pi0InAcc_Pt";
        fObjectNameMCMeson                          = "MC_Pi0_Pt";
        // fObjectNameMCMultipleMesonAcc               = "MC_Pi0InAccMulipleSubtraction_Pt";
        // fObjectNameMCMultipleMeson                  = "MC_Pi0MultipleSubtraction_Pt";
        fObjectNameMCMesonAccSecPi0                 = "MC_SecPi0InAcc_Pt_Source";
        fObjectNameMCMesonSecPi0                    = "MC_SecPi0_Pt_Source";
        fObjectNameMCMesonWOWeights                 = "MC_Pi0_WOWeights_Pt";
        fObjectNameMCMesonDalitzAcc                 = "MC_Pi0DalitzInAcc_Pt";
        fObjectNameMCMesonDalitz                    = "MC_Pi0Dalitz_Pt";
        fObjectNameMCMesonDalitzWOWeights           = "MC_Pi0Dalitz_WOWeights_Pt";
        fObjectNameTrueClusPrimMesonM02             = "ESD_TrueClusFromPrimPi0_Pt_M02";
        fObjectNameTrueClusSecMesonM02              = "ESD_TrueClusFromSecPi0_Pt_M02";
        fObjectNameTrueClusSecMesonFromK0sM02       = "ESD_TrueClusFromSecPi0FromK0s_Pt_M02";
        fObjectNameTrueClusSecMesonFromK0lM02       = "ESD_TrueClusFromSecPi0FromK0l_Pt_M02";
        fObjectNameTrueClusSecMesonFromLambdaM02    = "ESD_TrueClusFromSecPi0FromLambda_Pt_M02";
    } else {
        fObjectNameMCMesonAcc                       = "MC_EtaInAcc_Pt";
        fObjectNameMCMeson                          = "MC_Eta_Pt";
        fObjectNameMCMesonWOWeights                 = "MC_Eta_WOWeights_Pt";
        fObjectNameMCMesonDalitzAcc                 = "MC_EtaDalitzInAcc_Pt";
        fObjectNameMCMesonDalitz                    = "MC_EtaDalitz_Pt";
        fObjectNameMCMesonDalitzWOWeights           = "MC_EtaDalitz_WOWeights_Pt";

    }
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
            fHistoYieldToyMCSecInput[j]         = (TH1D*)fFileToyMCInput[j]->Get(Form("histoPiZeroDaughtersPt_InRap_%1.2f", fYMaxMeson/2));
            if (fHistoYieldToyMCSecInput[j]){
                nSecInputHistsFound++;
                fHistoYieldToyMCSecInput[j]->Sumw2();
                fHistoYieldToyMCSecInput[j]->SetName(Form("histoSecPi0YieldFrom%s_FromToy_orgBinning",nameSecondaries[j].Data()));

                cout << "found it, rebinning" << endl;
                fHistoYieldToyMCSecInputReb[j]  =  (TH1D*)fHistoYieldToyMCSecInput[j]->Rebin(fNBinsPt,Form("histoSecPi0YieldFrom%s_FromToy",nameSecondaries[j].Data()),fBinsPt); // Proper bins in Pt
                if (fHistoYieldToyMCSecInputReb[j]){
                    fHistoYieldToyMCSecInputReb[j]->Divide(fDeltaPt);
                    fHistoYieldToyMCSecInput[j]->Scale(1./fHistoYieldToyMCSecInput[j]->GetBinWidth(1));
                    cout << "that worked" << endl;
                }
            } else {
                cout << "file didn't contain proper histo" << endl;
            }
        } else {
           cout << "could not find correct input file" << endl;
        }
    }
//     cout << nSecInputHistsFound << endl;
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

            // fHistoYieldExternSecInput[j]           = (TH1D*)fFileCocktailInput->Get(Form("Pi0_From_%s_Pt_OrBin", nameSecondariesCocktail[j].Data()));
            // if (fHistoYieldExternSecInput[j]){
            //     fHistoYieldExternSecInput[j]->Sumw2();
            //     fHistoYieldExternSecInput[j]->SetName(Form("histoSecPi0YieldFrom%s_FromCocktail_orgBinning",nameSecondaries[j].Data())); // Proper bins in Pt

            //     fHistoYieldExternSecInputReb[j]     = (TH1D*)fHistoYieldExternSecInput[j]->Rebin(fNBinsPt,Form("histoSecPi0YieldFrom%s_FromCocktail",nameSecondaries[j].Data()),fBinsPt); // Proper bins in Pt
            //     if (fHistoYieldExternSecInputReb[j]){
            //         fHistoYieldExternSecInputReb[j]->Divide(fDeltaPt);
            //         fHistoYieldExternSecInput[j]->Scale(1./fHistoYieldExternSecInput[j]->GetBinWidth(1));
            //     }

            TH1D *nCocktailEvents = (TH1D*)fFileCocktailInput->Get("NEvents");

            // fHistoYieldExternSecInput[j]           = (TH1D*)fFileCocktailInput->Get(Form("Pi0_From_%s_Pt_Y_OrBin", nameSecondariesCocktail[j].Data()));
            TH2F *h2temp = 0x0;

            h2temp           = (TH2F*)fFileCocktailInput->Get(Form("Pi0_From_%s_Pt_Y_OrBin", nameSecondariesCocktail[j].Data()));
            if (h2temp){
                h2temp->Sumw2();
                h2temp->Scale(1/(h2temp->GetXaxis()->GetBinWidth(1)));
                fHistoYieldExternSecInput[j] = (TH1D*)h2temp->ProjectionX(Form("Pi0_from_%s_cocktail", nameSecondaries[j].Data()));
                fHistoYieldExternSecInput[j]->SetName(Form("histoSecPi0YieldFrom%s_FromCocktail_orgBinning",nameSecondaries[j].Data())); // Proper bins in Pt

                fHistoYieldExternSecInputReb[j]     = (TH1D*)fHistoYieldExternSecInput[j]->Rebin(fNBinsPt,Form("histoSecPi0YieldFrom%s_FromCocktail",nameSecondaries[j].Data()),fBinsPt); // Proper bins in Pt
                if (fHistoYieldExternSecInputReb[j]){
                    fHistoYieldExternSecInputReb[j]->Divide(fDeltaPt);
                    fHistoYieldExternSecInputReb[j]->Scale(1./nCocktailEvents->GetBinContent(1));
                    fHistoYieldExternSecInputReb[j]->SetDirectory(0);
                    fHistoYieldExternSecInput[j]->Scale(1./fHistoYieldExternSecInput[j]->GetBinWidth(1));

                printf("secondary histogram found for pi0 feed-down from %s\n", nameSecondaries[j].Data());
                }else{
                    printf("secondary histogram NOT FOUND for pi0 feed-down from %s\n", nameSecondaries[j].Data());
                }

            } else {
                cout << "file didn't contain " << Form("Pi0_From_%s_Pt_OrBin", nameSecondariesCocktail[j].Data()) << endl;
            }
        }
        return kTRUE;
    }
    return kFALSE;
}

//****************************************************************************
//*********************** creation of momentum dependent histos **************
//****************************************************************************
void CreatePtHistos(){

    fDeltaPt                            = new TH1D("deltaPt", "", fNBinsPt, fBinsPt);
    fDeltaPt->Sumw2();

    fHistoYieldMesonM02                 = new TH1D("histoYieldMesonM02", "", fNBinsPt, fBinsPt);
    fHistoYieldMesonM02->Sumw2();

    fHistoExoticsBelowM02MinBaseCut     = new TH1D("histoYieldExoticsBelowMinBaseCut", "", fNBinsPt, fBinsPt);
    fHistoExoticsBelowM02MinBaseCut->Sumw2();

    fHistoExoticsBelowM02MinCut         = new TH1D("histoYieldExoticsBelowMinCut", "", fNBinsPt, fBinsPt);
    fHistoExoticsBelowM02MinCut->Sumw2();
    fHistoExoticsAboveM02MinCut         = new TH1D("histoYieldExoticsAboveMinCut", "", fNBinsPt, fBinsPt);
    fHistoExoticsAboveM02MinCut->Sumw2();

    if (fIsMC){
        fHistoTrueYieldMergedM02                  = new TH1D("histoTrueYieldMergedM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedM02->Sumw2();
        fHistoTrueYieldMergedPartConvM02          = new TH1D("histoTrueYieldMergedPartConvM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedPartConvM02->Sumw2();
        fHistoTrueYieldMergedPartConvFromPi0M02   = new TH1D("histoTrueYieldMergedPartConvFromPi0M02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedPartConvFromPi0M02->Sumw2();
        fHistoTrueYieldMergedPartConvFromEtaM02   = new TH1D("histoTrueYieldMergedPartConvFromEtaM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedPartConvFromEtaM02->Sumw2();
        fHistoTrueYieldMergedPureM02              = new TH1D("histoTrueYieldMergedPureM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedPureM02->Sumw2();
        fHistoTrueYieldMergedPureFromPi0M02       = new TH1D("histoTrueYieldMergedPureFromPi0M02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedPureFromPi0M02->Sumw2();
        fHistoTrueYieldMergedPureFromEtaM02       = new TH1D("histoTrueYieldMergedPureFromEtaM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedPureFromEtaM02->Sumw2();
        fHistoTrueYieldMergedOneGammaM02          = new TH1D("histoTrueYieldMergedOneGammaFromMesonM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedOneGammaM02->Sumw2();
        fHistoTrueYieldMergedOneGammaFromPi0M02   = new TH1D("histoTrueYieldMergedOneGammaFromPi0M02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedOneGammaFromPi0M02->Sumw2();
        fHistoTrueYieldMergedOneGammaFromEtaM02   = new TH1D("histoTrueYieldMergedOneGammaFromEtaM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedOneGammaFromEtaM02->Sumw2();
        fHistoTrueYieldMergedOneElectronM02       = new TH1D("histoTrueYieldMergedOneElectronFromMesonM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedOneElectronM02->Sumw2();
        fHistoTrueYieldMergedOneElectronFromPi0M02= new TH1D("histoTrueYieldMergedOneElectronFromPi0M02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedOneElectronFromPi0M02->Sumw2();
        fHistoTrueYieldMergedOneElectronFromEtaM02= new TH1D("histoTrueYieldMergedOneElectronFromEtaM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldMergedOneElectronFromEtaM02->Sumw2();
        fHistoTrueYieldPi0M02                     = new TH1D("histoTrueYieldPi0M02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldPi0M02->Sumw2();
        fHistoTrueYieldPi0DCM02                   = new TH1D("histoTrueYieldPi0DCM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldPi0DCM02->Sumw2();
        fHistoTrueYieldPi0GGM02                   = new TH1D("histoTrueYieldPi0GGM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldPi0GGM02->Sumw2();
        fHistoTrueYieldPi0DalitzM02               = new TH1D("histoTrueYieldPi0DalitzM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldPi0DalitzM02->Sumw2();
        fHistoTrueYieldEtaM02                     = new TH1D("histoTrueYieldEtaM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldEtaM02->Sumw2();
        fHistoTrueYieldEtaDCM02                   = new TH1D("histoTrueYieldEtaDCM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldEtaDCM02->Sumw2();
        fHistoTrueYieldEtaGGM02                   = new TH1D("histoTrueYieldEtaGGM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldEtaGGM02->Sumw2();
        fHistoTrueYieldEtaDalitzM02               = new TH1D("histoTrueYieldEtaDalitzM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldEtaDalitzM02->Sumw2();
        fHistoTrueYieldGammaM02                   = new TH1D("histoTrueYieldGammaM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldGammaM02->Sumw2();
        fHistoTrueYieldElectronM02                = new TH1D("histoTrueYieldElectronM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldElectronM02->Sumw2();
        fHistoTrueYieldBGM02                      = new TH1D("histoTrueYieldBGM02", "", fNBinsPt, fBinsPt);
        fHistoTrueYieldBGM02->Sumw2();
        if (fPrefix.Contains("Pi0")){
            fHistoTrueYieldPrimPi0M02                   = new TH1D("histoTrueYieldPrimPi0M02", "", fNBinsPt, fBinsPt);
            fHistoTrueYieldPrimPi0M02->Sumw2();
            for (Int_t j = 0; j< 4; j++){
                fHistoTrueYieldSecPi0M02[j]             = new TH1D(Form("histoTrueYieldSecPi0F%sM02", nameSecondaries[j].Data()), "", fNBinsPt, fBinsPt);
                fHistoTrueYieldSecPi0M02[j]->Sumw2();
            }
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


        if (fEnableExoticsQA){
            if (fHistoExoticsM02EBin[iPt-1]){
                Double_t integErr   = 0;
                Double_t integ      = fHistoExoticsM02EBin[iPt-1]->IntegralAndError(fHistoExoticsM02EBin[iPt-1]->FindBin(0), fHistoExoticsM02EBin[iPt-1]->FindBin(0.1),integErr);

                fHistoExoticsBelowM02MinBaseCut->SetBinContent(iPt, integ);
                fHistoExoticsBelowM02MinBaseCut->SetBinError(iPt, integErr);

                integErr            = 0;
                integ               = fHistoExoticsM02EBin[iPt-1]->IntegralAndError(fHistoExoticsM02EBin[iPt-1]->FindBin(0.1+0.0001), fHistoExoticsM02EBin[iPt-1]->FindBin(0.27),integErr);

                fHistoExoticsBelowM02MinCut->SetBinContent(iPt, integ);
                fHistoExoticsBelowM02MinCut->SetBinError(iPt, integErr);

                integErr            = 0;
                integ               = fHistoExoticsM02EBin[iPt-1]->IntegralAndError(fHistoExoticsM02EBin[iPt-1]->FindBin(0.27+0.0001), fHistoExoticsM02EBin[iPt-1]->GetNbinsX(),integErr);
                fHistoExoticsAboveM02MinCut->SetBinContent(iPt, integ);
                fHistoExoticsAboveM02MinCut->SetBinError(iPt, integErr);
            }

        }
        if (fIsMC){
            fHistoTrueYieldMergedM02->SetBinContent(iPt, fMesonM02TrueMergedYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedM02->SetBinError(iPt, fMesonM02TrueMergedYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldMergedPartConvM02->SetBinContent(iPt, fMesonM02TrueMergedPartConvYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedPartConvM02->SetBinError(iPt, fMesonM02TrueMergedPartConvYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedPartConvFromPi0M02->SetBinContent(iPt, fMesonM02TrueMergedPartConvFromPi0Yields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedPartConvFromPi0M02->SetBinError(iPt, fMesonM02TrueMergedPartConvFromPi0YieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedPartConvFromEtaM02->SetBinContent(iPt, fMesonM02TrueMergedPartConvFromEtaYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedPartConvFromEtaM02->SetBinError(iPt, fMesonM02TrueMergedPartConvFromEtaYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldMergedPureM02->SetBinContent(iPt, fMesonM02TrueMergedPureYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedPureM02->SetBinError(iPt, fMesonM02TrueMergedPureYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedPureFromPi0M02->SetBinContent(iPt, fMesonM02TrueMergedPureFromPi0Yields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedPureFromPi0M02->SetBinError(iPt, fMesonM02TrueMergedPureFromPi0YieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedPureFromEtaM02->SetBinContent(iPt, fMesonM02TrueMergedPureFromEtaYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedPureFromEtaM02->SetBinError(iPt, fMesonM02TrueMergedPureFromEtaYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldMergedOneGammaM02->SetBinContent(iPt, fMesonM02TrueMergedOneGammaYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedOneGammaM02->SetBinError(iPt, fMesonM02TrueMergedOneGammaYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldMergedOneGammaFromPi0M02->SetBinContent(iPt, fMesonM02TrueMergedOneGammaFromPi0Yields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedOneGammaFromPi0M02->SetBinError(iPt, fMesonM02TrueMergedOneGammaFromPi0YieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldMergedOneGammaFromEtaM02->SetBinContent(iPt, fMesonM02TrueMergedOneGammaFromEtaYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedOneGammaFromEtaM02->SetBinError(iPt, fMesonM02TrueMergedOneGammaFromEtaYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldMergedOneElectronM02->SetBinContent(iPt, fMesonM02TrueMergedOneElectronYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedOneElectronM02->SetBinError(iPt, fMesonM02TrueMergedOneElectronYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldMergedOneElectronFromPi0M02->SetBinContent(iPt, fMesonM02TrueMergedOneElectronFromPi0Yields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedOneElectronFromPi0M02->SetBinError(iPt, fMesonM02TrueMergedOneElectronFromPi0YieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldMergedOneElectronFromEtaM02->SetBinContent(iPt, fMesonM02TrueMergedOneElectronFromEtaYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldMergedOneElectronFromEtaM02->SetBinError(iPt, fMesonM02TrueMergedOneElectronFromEtaYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldPi0M02->SetBinContent(iPt, fMesonM02TruePi0Yields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldPi0M02->SetBinError(iPt, fMesonM02TruePi0YieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldPi0DCM02->SetBinContent(iPt, fMesonM02TruePi0DCYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldPi0DCM02->SetBinError(iPt, fMesonM02TruePi0DCYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldPi0GGM02->SetBinContent(iPt, fMesonM02TruePi0GGYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldPi0GGM02->SetBinError(iPt, fMesonM02TruePi0GGYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldPi0DalitzM02->SetBinContent(iPt, fMesonM02TruePi0DalitzYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldPi0DalitzM02->SetBinError(iPt, fMesonM02TruePi0DalitzYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldEtaM02->SetBinContent(iPt, fMesonM02TrueEtaYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldEtaM02->SetBinError(iPt, fMesonM02TrueEtaYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldEtaDCM02->SetBinContent(iPt, fMesonM02TrueEtaDCYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldEtaDCM02->SetBinError(iPt, fMesonM02TrueEtaDCYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldEtaGGM02->SetBinContent(iPt, fMesonM02TrueEtaGGYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldEtaGGM02->SetBinError(iPt, fMesonM02TrueEtaGGYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldEtaDalitzM02->SetBinContent(iPt, fMesonM02TrueEtaDalitzYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldEtaDalitzM02->SetBinError(iPt, fMesonM02TrueEtaDalitzYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldGammaM02->SetBinContent(iPt, fMesonM02TrueGammaYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldGammaM02->SetBinError(iPt, fMesonM02TrueGammaYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldElectronM02->SetBinContent(iPt, fMesonM02TrueElectronYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldElectronM02->SetBinError(iPt, fMesonM02TrueElectronYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            fHistoTrueYieldBGM02->SetBinContent(iPt, fMesonM02TrueBGYields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoTrueYieldBGM02->SetBinError(iPt, fMesonM02TrueBGYieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

            if (fPrefix.Contains("Pi0")){
                fHistoTrueYieldPrimPi0M02->SetBinContent(iPt, fMesonM02TruePrimPi0Yields[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoTrueYieldPrimPi0M02->SetBinError(iPt, fMesonM02TruePrimPi0YieldsError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

                for (Int_t j = 0; j < 4; j++){
                    fHistoTrueYieldSecPi0M02[j]->SetBinContent(iPt, fMesonM02TrueSecPi0Yields[j][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                    fHistoTrueYieldSecPi0M02[j]->SetBinError(iPt, fMesonM02TrueSecPi0YieldsError[j][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                }
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
TH1D* FillProjectionY (TH2F* fDummy2D, TString name, Double_t minX, Double_t maxX, Int_t rebin, Double_t* array){
    TH1D* dummy1D           = new TH1D(name.Data(), name.Data(), fDummy2D->GetNbinsY(), 0., fDummy2D->GetYaxis()->GetBinUpEdge(fDummy2D->GetNbinsY()));
    dummy1D->Sumw2();
    Int_t startBin          = fDummy2D->GetXaxis()->FindBin(minX+0.001);
    Int_t endBin            = fDummy2D->GetXaxis()->FindBin(maxX-0.001);
    fDummy2D->ProjectionY(name.Data(),startBin,endBin,"e");
    dummy1D                 = (TH1D*)gDirectory->Get(name.Data());
    if (array != NULL){
        TH1D* dummy2        = (TH1D*)dummy1D->Rebin(rebin,name.Data(),array);
        return dummy2;
    } else {
        if(rebin>1){
            dummy1D->Rebin(rebin);
        }
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
            fHistoM02PtBin[iPt]         = FillProjectionY(fM02VsPtDummy, fNameHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }
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
                            TH2F* fTrueBGPtVsM02
                         ) {
    cout << __LINE__ << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        TString fNameHistoM02                           = Form("TrueClusMerged_M02_in_Pt_Bin%02d", iPt);
        TString fNamePi0HistoM02                        = Form("TrueClusPi0_M02_in_Pt_Bin%02d", iPt);
        TString fNameEtaHistoM02                        = Form("TrueClusEta_M02_in_Pt_Bin%02d", iPt);
        TString fNameGammaHistoM02                      = Form("TrueClusGamma_M02_in_Pt_Bin%02d", iPt);
        TString fNameElectronHistoM02                   = Form("TrueClusElectron_M02_in_Pt_Bin%02d", iPt);
        TString fNameBGHistoM02                         = Form("TrueClusBG_M02_in_Pt_Bin%02d", iPt);

        // true merged clusters
        CheckForNULLForPointer(fHistoTrueClusMergedM02PtBin[iPt]);
        if (fTrueMergedPtVsM02Dummy){
            fHistoTrueClusMergedM02PtBin[iPt]                   = FillProjectionY(fTrueMergedPtVsM02Dummy, fNameHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }
        // true merged pi0
        CheckForNULLForPointer(fHistoTrueClusPi0M02PtBin[iPt]);
        if (fTruePi0PtVsM02){
            fHistoTrueClusPi0M02PtBin[iPt]                      = FillProjectionY(fTruePi0PtVsM02, fNamePi0HistoM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // true merged eta
        CheckForNULLForPointer(fHistoTrueClusEtaM02PtBin[iPt]);
        if (fTrueEtaPtVsM02){
            fHistoTrueClusEtaM02PtBin[iPt]                      = FillProjectionY(fTrueEtaPtVsM02, fNameEtaHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // true gamma
        CheckForNULLForPointer(fHistoTrueClusGammaM02PtBin[iPt]);
        if (fTrueGammaPtVsM02){
            fHistoTrueClusGammaM02PtBin[iPt]                    = FillProjectionY(fTrueGammaPtVsM02, fNameGammaHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // true electron
        CheckForNULLForPointer(fHistoTrueClusElectronM02PtBin[iPt]);
        if (fTrueElectronPtVsM02){
            fHistoTrueClusElectronM02PtBin[iPt]                 = FillProjectionY(fTrueElectronPtVsM02, fNameElectronHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // true BG
        CheckForNULLForPointer(fHistoTrueClusBGM02PtBin[iPt]);
        if (fTrueBGPtVsM02){
            fHistoTrueClusBGM02PtBin[iPt]                       = FillProjectionY(fTrueBGPtVsM02, fNameBGHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }
    }
}

//****************************************************************************
//*************** Fill array of M02 histograms in pT slices for MC ***********
//****************************************************************************
void FillMCM02AdditionHistosArray(      TH2F* fTrueMergedPureDummy,
                                        TH2F* fTruePartConvMergedVsM02Dummy,
                                        TH2F* fTrueOneGammaDummy,
                                        TH2F* fTrueOneGammaFromPi0Dummy,
                                        TH2F* fTrueOneGammaFromEtaDummy,
                                        TH2F* fTrueOneElectronDummy,
                                        TH2F* fTrueOneElectronFromPi0Dummy,
                                        TH2F* fTrueOneElectronFromEtaDummy,
                                        TH2F* fTruePi0GGDummy,
                                        TH2F* fTruePi0DalitzDummy,
                                        TH2F* fTrueEtaGGDummy,
                                        TH2F* fTrueEtaDalitzDummy,
                                        TH2F* fTruePi0DCDummy,
                                        TH2F* fTrueEtaDCDummy,
                                        TH2F* fTrueMergedPurePi0Dummy,
                                        TH2F* fTrueMergedPureEtaDummy,
                                        TH2F* fTruePartConvMergedPi0VsM02Dummy,
                                        TH2F* fTruePartConvMergedEtaVsM02Dummy

                                  ) {
    cout << __LINE__ << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        TString fNameHistoM02                           = Form("TrueClusMergedPure_M02_in_Pt_Bin%02d", iPt);
        TString fNamePartConvHistoM02                   = Form("TrueClusPartConvMerged_M02_in_Pt_Bin%02d", iPt);
        TString fNameOneGammaHistoM02                   = Form("TrueClusOneGamma_M02_in_Pt_Bin%02d", iPt);
        TString fNameOneGammaFromPi0HistoM02            = Form("TrueClusOneGammaFromPi0_M02_in_Pt_Bin%02d", iPt);
        TString fNameOneGammaFromEtaHistoM02            = Form("TrueClusOneGammaFromEta_M02_in_Pt_Bin%02d", iPt);
        TString fNameOneElectronHistoM02                = Form("TrueClusOneElectron_M02_in_Pt_Bin%02d", iPt);
        TString fNameOneElectronFromPi0HistoM02         = Form("TrueClusOneElectronFromPi0_M02_in_Pt_Bin%02d", iPt);
        TString fNameOneElectronFromEtaHistoM02         = Form("TrueClusOneElectronFromEta_M02_in_Pt_Bin%02d", iPt);
        TString fNameTruePi0GGM02                       = Form("TrueClusPi0GG_M02_in_Pt_Bin%02d", iPt);
        TString fNameTruePi0DalitzM02                   = Form("TrueClusPi0Dalitz_M02_in_Pt_Bin%02d", iPt);
        TString fNameTrueEtaGGM02                       = Form("TrueClusEtaGG_M02_in_Pt_Bin%02d", iPt);
        TString fNameTrueEtaDalitzM02                   = Form("TrueClusEtaDalitz_M02_in_Pt_Bin%02d", iPt);
        TString fNameTruePi0DCM02                       = Form("TrueClusPi0DC_M02_in_Pt_Bin%02d", iPt);
        TString fNameTrueEtaDCM02                       = Form("TrueClusEtaDC_M02_in_Pt_Bin%02d", iPt);
        TString fNameHistoM02FromPi0                    = Form("TrueClusMergedPureFromPi0_M02_in_Pt_Bin%02d", iPt);
        TString fNamePartConvHistoM02FromPi0            = Form("TrueClusPartConvMergedFromPi0_M02_in_Pt_Bin%02d", iPt);
        TString fNameHistoM02FromEta                    = Form("TrueClusMergedPureFromEta_M02_in_Pt_Bin%02d", iPt);
        TString fNamePartConvHistoM02FromEta            = Form("TrueClusPartConvMergedFromEta_M02_in_Pt_Bin%02d", iPt);

        // true merged clusters
        CheckForNULLForPointer(fHistoTrueClusPureMergedM02PtBin[iPt]);
        if (fTrueMergedPureDummy){
            fHistoTrueClusPureMergedM02PtBin[iPt]               = FillProjectionY(fTrueMergedPureDummy, fNameHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }
        CheckForNULLForPointer(fHistoTrueClusPureMergedFromPi0M02PtBin[iPt]);
        if (fTrueMergedPurePi0Dummy){
            fHistoTrueClusPureMergedFromPi0M02PtBin[iPt]        = FillProjectionY(fTrueMergedPurePi0Dummy, fNameHistoM02FromPi0, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }
        CheckForNULLForPointer(fHistoTrueClusPureMergedFromEtaM02PtBin[iPt]);
        if (fTrueMergedPureEtaDummy){
            fHistoTrueClusPureMergedFromEtaM02PtBin[iPt]        = FillProjectionY(fTrueMergedPureEtaDummy, fNameHistoM02FromEta, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }
        // true merged clusters part conv
        CheckForNULLForPointer(fHistoTrueClusPartConvMergedM02PtBin[iPt]);
        if (fTruePartConvMergedVsM02Dummy){
            fHistoTrueClusPartConvMergedM02PtBin[iPt]           = FillProjectionY(fTruePartConvMergedVsM02Dummy, fNamePartConvHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }
        CheckForNULLForPointer(fHistoTrueClusPartConvMergedFromPi0M02PtBin[iPt]);
        if (fTruePartConvMergedPi0VsM02Dummy){
            fHistoTrueClusPartConvMergedFromPi0M02PtBin[iPt]    = FillProjectionY(fTruePartConvMergedPi0VsM02Dummy, fNamePartConvHistoM02FromPi0, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }
        CheckForNULLForPointer(fHistoTrueClusPartConvMergedFromEtaM02PtBin[iPt]);
        if (fTruePartConvMergedEtaVsM02Dummy){
            fHistoTrueClusPartConvMergedFromEtaM02PtBin[iPt]    = FillProjectionY(fTruePartConvMergedEtaVsM02Dummy, fNamePartConvHistoM02FromEta, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }
        // true cluster with 1 gamma from mesons
        CheckForNULLForPointer(fHistoTrueClusOneGammaM02PtBin[iPt]);
        if (fTrueOneGammaDummy){
            fHistoTrueClusOneGammaM02PtBin[iPt]                 = FillProjectionY(fTrueOneGammaDummy, fNameOneGammaHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // true cluster with 1 gamma from pi0
        CheckForNULLForPointer(fHistoTrueClusOneGammaFromPi0M02PtBin[iPt]);
        if (fTrueOneGammaFromPi0Dummy){
            fHistoTrueClusOneGammaFromPi0M02PtBin[iPt]          = FillProjectionY(fTrueOneGammaFromPi0Dummy, fNameOneGammaFromPi0HistoM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // true cluster with 1 gamma from eta
        CheckForNULLForPointer(fHistoTrueClusOneGammaFromEtaM02PtBin[iPt]);
        if (fTrueOneGammaFromEtaDummy){
            fHistoTrueClusOneGammaFromEtaM02PtBin[iPt]          = FillProjectionY(fTrueOneGammaFromEtaDummy, fNameOneGammaFromEtaHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // true cluster with 1 electron from meson
        CheckForNULLForPointer(fHistoTrueClusOneElectronM02PtBin[iPt]);
        if (fTrueOneElectronDummy){
            fHistoTrueClusOneElectronM02PtBin[iPt]              = FillProjectionY(fTrueOneElectronDummy, fNameOneElectronHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // true cluster with 1 electron from pi0
        CheckForNULLForPointer(fHistoTrueClusOneElectronFromPi0M02PtBin[iPt]);
        if (fTrueOneElectronFromPi0Dummy){
            fHistoTrueClusOneElectronFromPi0M02PtBin[iPt]       = FillProjectionY(fTrueOneElectronFromPi0Dummy, fNameOneElectronFromPi0HistoM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // true cluster with 1 electron from eta
        CheckForNULLForPointer(fHistoTrueClusOneElectronFromEtaM02PtBin[iPt]);
        if (fTrueOneElectronFromEtaDummy){
            fHistoTrueClusOneElectronFromEtaM02PtBin[iPt]       = FillProjectionY(fTrueOneElectronFromEtaDummy, fNameOneElectronFromEtaHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // true cluster from pi0 in GG channel
        CheckForNULLForPointer(fHistoTrueClusPi0GGM02PtBin[iPt]);
        if (fTruePi0GGDummy){
            fHistoTrueClusPi0GGM02PtBin[iPt]                    = FillProjectionY(fTruePi0GGDummy, fNameTruePi0GGM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // true cluster from pi0 in Dalitz channel
        CheckForNULLForPointer(fHistoTrueClusPi0DalitzM02PtBin[iPt]);
        if (fTruePi0DalitzDummy){
            fHistoTrueClusPi0DalitzM02PtBin[iPt]                = FillProjectionY(fTruePi0DalitzDummy, fNameTruePi0DalitzM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // true cluster from eta in GG channel
        CheckForNULLForPointer(fHistoTrueClusEtaGGM02PtBin[iPt]);
        if (fTrueEtaGGDummy){
            fHistoTrueClusEtaGGM02PtBin[iPt]                    = FillProjectionY(fTrueEtaGGDummy, fNameTrueEtaGGM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // true cluster from eta in Dalitz channel
        CheckForNULLForPointer(fHistoTrueClusEtaDalitzM02PtBin[iPt]);
        if (fTrueEtaDalitzDummy){
            fHistoTrueClusEtaDalitzM02PtBin[iPt]                = FillProjectionY(fTrueEtaDalitzDummy, fNameTrueEtaDalitzM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // true cluster from pi0 double counted
        CheckForNULLForPointer(fHistoTrueClusPi0DCM02PtBin[iPt]);
        if (fTruePi0DCDummy){
            fHistoTrueClusPi0DCM02PtBin[iPt]                    = FillProjectionY(fTruePi0DCDummy, fNameTruePi0DCM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // true cluster from eta double counted
        CheckForNULLForPointer(fHistoTrueClusEtaDCM02PtBin[iPt]);
        if (fTrueEtaDCDummy){
            fHistoTrueClusEtaDCM02PtBin[iPt]                    = FillProjectionY(fTrueEtaDCDummy, fNameTrueEtaDCM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

    }
}


//****************************************************************************
//******** Fill array of M02 histograms in pT slices for Pi0 Prim/Sec ********
//****************************************************************************
void FillMCPrimSecM02HistosArray(   TH2F* fTruePi0PrimM02VsPt,
                                    TH2F** fTruePi0SecM02VsPt
                                ) {

    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        TString fNamePi0PrimM02                 = Form("TrueClusPrimPi0_M02_in_Pt_Bin%02d", iPt);

        // true prim pi0 clusters
        CheckForNULLForPointer(fHistoTrueClusPrimPi0M02PtBin[iPt]);
        if (fTruePi0PrimM02VsPt){
            fHistoTrueClusPrimPi0M02PtBin[iPt]                  = FillProjectionY(fTruePi0PrimM02VsPt, fNamePi0PrimM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
        }

        // fill secondary pt differential histo for pi0 from: K0s, Lambda, K0l, Rest
        if (fTruePi0SecM02VsPt){
            for (Int_t j = 0; j< 4; j++){
                TString fNamePi0SecM02                  = Form("TrueClusSecPi0From%s_M02_in_Pt_Bin%02d", nameSecondaries[j].Data(), iPt);
                CheckForNULLForPointer(fHistoTrueClusSecPi0M02PtBin[j][iPt]);
                if (fTruePi0SecM02VsPt[j]){
                    fHistoTrueClusSecPi0M02PtBin[j][iPt]            = FillProjectionY(fTruePi0SecM02VsPt[j], fNamePi0SecM02, fBinsPt[iPt], fBinsPt[iPt+1], fNRebin[iPt], NULL);
                }
            }
        }
    }
}

//****************************************************************************
//*************** Fill array of M02 histograms in pT slices for MC ***********
//****************************************************************************
void FillExoticsM02HistosArray(  TH2F* fExoticsM02 ) {
    cout << __LINE__ << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        TString fNameHistoM02                           = Form("ExoticClusters_M02_in_E_Bin%02d", iPt);
        // true merged clusters
        CheckForNULLForPointer(fHistoExoticsM02EBin[iPt]);
        if (fExoticsM02){
            fHistoExoticsM02EBin[iPt]                   = FillProjectionY(fExoticsM02, fNameHistoM02, fBinsPt[iPt], fBinsPt[iPt+1], 1, NULL);
        }
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
            fHistoTrueClustersBGPt[i]->Divide(fDeltaPt);
            if (i == 0) fHistoTrueClustersBGPt[9]    = (TH1D*)fHistoTrueClustersBGPt[i]->Clone(Form("TrueClusBG_%s_Pt",labelsBG[9].Data()));
                else fHistoTrueClustersBGPt[9]->Add(fHistoTrueClustersBGPt[i]);
        }
    }
}

//****************************************************************************
//***************** Filling MC Electron separated **********************************
//****************************************************************************
void FillMCElectronSeparated (TH2F* dummy2D){
    if (dummy2D){
        for (Int_t i = 0; i< 9; i++){
            TH1D* dummy1D                   = (TH1D*)dummy2D->ProjectionX(Form("TrueClusElectron_%s_Pt",labelsElectron[i].Data()),i+1,i+1,"e");
            fHistoTrueClustersElectronPt[i] = (TH1D*)dummy1D->Rebin(fNBinsPt,Form("TrueClusElectron_%s_Pt",labelsElectron[i].Data()),fBinsPt);
            fHistoTrueClustersElectronPt[i]->Divide(fDeltaPt);
        }
    }
}

//****************************************************************************
//***************** Filling MC Gamma separated **********************************
//****************************************************************************
void FillMCGammaSeparated (TH2F* dummy2D){
    if (dummy2D){
        for (Int_t i = 0; i< 8; i++){
            TH1D* dummy1D                   = (TH1D*)dummy2D->ProjectionX(Form("TrueClusGamma_%s_Pt",labelsGamma[i].Data()),i+1,i+1,"e");
            fHistoTrueClustersGammaPt[i]    = (TH1D*)dummy1D->Rebin(fNBinsPt,Form("TrueClusGamma_%s_Pt",labelsGamma[i].Data()),fBinsPt);
            fHistoTrueClustersGammaPt[i]->Divide(fDeltaPt);
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
//***************** Filling of MC histograms in proper binning ***************
//****************************************************************************
void FillHistosArrayMCSecAndCalcAcceptance( TH2D* mcSecInputInAccSourcePt, TH2D * mcSecInputSourcePt ) {

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

    // Filling of secondary histos with different sources
    for (Int_t j = 0; j < 4; j++){
        // project correct bin
        fHistoMCSecPi0Pt[j]             = (TH1D*)mcSecInputSourcePt->ProjectionX(Form("MCSecPi0From%s",nameSecondaries[j].Data()), mcSecInputSourcePt->GetYaxis()->FindBin(startBin[j]),
                                                                                 mcSecInputSourcePt->GetYaxis()->FindBin(endBin[j]),"e");
        fHistoMCSecPi0Pt[j]->SetTitle(Form("MCSecPi0From%s",nameSecondaries[j].Data()));
        fHistoMCSecPi0Pt[j]->Sumw2();
        fHistoMCSecPi0PtWAcc[j]         = (TH1D*)mcSecInputInAccSourcePt->ProjectionX(Form("MCSecPi0From%s_InAcc",nameSecondaries[j].Data()), mcSecInputInAccSourcePt->GetYaxis()->FindBin(startBin[j]),
                                                                                      mcSecInputInAccSourcePt->GetYaxis()->FindBin(endBin[j]),"e");
        fHistoMCSecPi0PtWAcc[j]->SetTitle(Form("MCSecPi0From%s_InAcc",nameSecondaries[j].Data()));
        fHistoMCSecPi0PtWAcc[j]->Sumw2();

        // rebin in to current analysis binning
        fHistoMCSecPi0PtReb[j]          = (TH1D*)fHistoMCSecPi0Pt[j]->Rebin(fNBinsPt,Form("MCSecPi0From%s_Rebinned",nameSecondaries[j].Data()),fBinsPt); // Proper bins in Pt
        fHistoMCSecPi0PtReb[j]->SetTitle(Form("MCSecPi0From%s_Rebinned",nameSecondaries[j].Data()));
        fHistoMCSecPi0PtReb[j] ->Divide(fDeltaPt);
        fHistoMCSecPi0PtWAccReb[j]      = (TH1D*)fHistoMCSecPi0PtWAcc[j]->Rebin(fNBinsPt,Form("MCSecPi0From%s_InAcc_Rebinned",nameSecondaries[j].Data()),fBinsPt); // Proper bins in Pt
        fHistoMCSecPi0PtWAccReb[j]->SetTitle(Form("MCSecPi0From%s_InAcc_Rebinned",nameSecondaries[j].Data()));
        fHistoMCSecPi0PtWAccReb[j] ->Divide(fDeltaPt);

        // calculate acceptance
        fHistoMCSecPi0AcceptPt[j]       = (TH1D*)fHistoMCSecPi0PtWAccReb[j]->Clone(Form("fMCSecPi0From%sAccepPt",nameSecondaries[j].Data()));
        fHistoMCSecPi0AcceptPt[j]->SetTitle(Form("fMCSecPi0From%sAccepPt",nameSecondaries[j].Data()));
        fHistoMCSecPi0AcceptPt[j]->Divide(fHistoMCSecPi0PtWAccReb[j],fHistoMCSecPi0PtReb[j],1.,1.,"B");
    }
}

//****************************************************************************
//***** Create ratios of reconstructed yield and different backgrounds *******
//****************************************************************************
void CreateRatioHistos(){
    for (Int_t i = 0; i < 10; i++){
        fHistoTrueClustersBGPt[i]->Sumw2();
        fHistoRatioTrueClustersBGPt[i]           = (TH1D*)fHistoTrueClustersBGPt[i]->Clone(Form("RatioTrueClusBG_%s_Pt",labelsBG[i].Data()));
        fHistoRatioTrueClustersBGPt[i]->Divide(fHistoRatioTrueClustersBGPt[i],fHistoYieldMesonM02,1.,1.,"e");
    }

    fHistoRatioTrueYieldEtaM02                   = (TH1D*)fHistoTrueYieldEtaM02->Clone("RatioTrueYieldEtaM02");
    fHistoRatioTrueYieldEtaM02->Divide(fHistoRatioTrueYieldEtaM02,fHistoYieldMesonM02,1.,1.,"");

    fHistoRatioTrueYieldPi0M02                   = (TH1D*)fHistoTrueYieldPi0M02->Clone("RatioTrueYieldPi0M02");
    fHistoRatioTrueYieldPi0M02->Divide(fHistoRatioTrueYieldPi0M02,fHistoYieldMesonM02,1.,1.,"");
    fHistoRatioTrueYieldGammaM02                 = (TH1D*)fHistoTrueYieldGammaM02->Clone("RatioTrueYieldGammaM02");
    fHistoRatioTrueYieldGammaM02->Divide(fHistoRatioTrueYieldGammaM02,fHistoYieldMesonM02,1.,1.,"");
    fHistoRatioTrueYieldElectronM02              = (TH1D*)fHistoTrueYieldElectronM02->Clone("RatioTrueYieldElectronM02");
    fHistoRatioTrueYieldElectronM02->Divide(fHistoRatioTrueYieldElectronM02,fHistoYieldMesonM02,1.,1.,"");

    fHistoRatioPi0DCFrac                         = (TH1D*)fHistoTrueYieldPi0DCM02->Clone("RatioPi0DCFrac");
    fHistoRatioPi0DCFrac->Divide(fHistoRatioPi0DCFrac,fHistoTrueYieldPi0M02,1.,1.,"");

    fHistoRatioPi0GGFrac                         = (TH1D*)fHistoTrueYieldPi0GGM02->Clone("RatioPi0GGFrac");
    fHistoRatioPi0GGFrac->Divide(fHistoRatioPi0GGFrac,fHistoTrueYieldPi0M02,1.,1.,"B");
    fHistoRatioPi0DalitzFrac                     = (TH1D*)fHistoTrueYieldPi0DalitzM02->Clone("RatioPi0DalitzFrac");
    fHistoRatioPi0DalitzFrac->Divide(fHistoRatioPi0DalitzFrac,fHistoTrueYieldPi0M02,1.,1.,"B");

    fHistoRatioEtaDCFrac                         = (TH1D*)fHistoTrueYieldEtaDCM02->Clone("RatioEtaDCFrac");
    fHistoRatioEtaDCFrac->Divide(fHistoRatioEtaDCFrac,fHistoTrueYieldEtaM02,1.,1.,"B");

    fHistoRatioEtaGGFrac                         = (TH1D*)fHistoTrueYieldEtaGGM02->Clone("RatioEtaGGFrac");
    fHistoRatioEtaGGFrac->Divide(fHistoRatioEtaGGFrac,fHistoTrueYieldEtaM02,1.,1.,"B");

    fHistoRatioEtaDalitzFrac                     = (TH1D*)fHistoTrueYieldEtaDalitzM02->Clone("RatioEtaDalitzFrac");
    fHistoRatioEtaDalitzFrac->Divide(fHistoRatioEtaDalitzFrac,fHistoTrueYieldEtaM02,1.,1.,"B");

    fHistoRatioMergedPureFracPi0                    = (TH1D*)fHistoTrueYieldMergedPureFromPi0M02->Clone("RatioPi0MergedPure");
    fHistoRatioMergedPureFracPi0->Divide(fHistoRatioMergedPureFracPi0,fHistoTrueYieldPi0M02,1.,1.,"B");
    fHistoRatioMergedPureFracPi0->Scale(100.);
    fHistoRatioMergedPartConvFracPi0                = (TH1D*)fHistoTrueYieldMergedPartConvFromPi0M02->Clone("RatioPi0MergedPartConv");
    fHistoRatioMergedPartConvFracPi0->Divide(fHistoRatioMergedPartConvFracPi0,fHistoTrueYieldPi0M02,1.,1.,"B");
    fHistoRatioMergedPartConvFracPi0->Scale(100.);
    fHistoRatioMergedOneGammaFracPi0                = (TH1D*)fHistoTrueYieldMergedOneGammaFromPi0M02->Clone("RatioPi0MergedOneGamma");
    fHistoRatioMergedOneGammaFracPi0->Divide(fHistoRatioMergedOneGammaFracPi0,fHistoTrueYieldPi0M02,1.,1.,"B");
    fHistoRatioMergedOneGammaFracPi0->Scale(100.);
    fHistoRatioMergedOneElectronFracPi0             = (TH1D*)fHistoTrueYieldMergedOneElectronFromPi0M02->Clone("RatioPi0MergedOneElectron");
    fHistoRatioMergedOneElectronFracPi0->Divide(fHistoRatioMergedOneElectronFracPi0,fHistoTrueYieldPi0M02,1.,1.,"B");
    fHistoRatioMergedOneElectronFracPi0->Scale(100.);

    fHistoRatioMergedPureFracEta                    = (TH1D*)fHistoTrueYieldMergedPureFromEtaM02->Clone("RatioEtaMergedPure");
    fHistoRatioMergedPureFracEta->Divide(fHistoRatioMergedPureFracEta,fHistoTrueYieldEtaM02,1.,1.,"B");
    fHistoRatioMergedPureFracEta->Scale(100.);
    fHistoRatioMergedPartConvFracEta                = (TH1D*)fHistoTrueYieldMergedPartConvFromEtaM02->Clone("RatioEtaMergedPartConv");
    fHistoRatioMergedPartConvFracEta->Divide(fHistoRatioMergedPartConvFracEta,fHistoTrueYieldEtaM02,1.,1.,"B");
    fHistoRatioMergedPartConvFracEta->Scale(100.);
    fHistoRatioMergedOneGammaFracEta                = (TH1D*)fHistoTrueYieldMergedOneGammaFromEtaM02->Clone("RatioEtaMergedOneGamma");
    fHistoRatioMergedOneGammaFracEta->Divide(fHistoRatioMergedOneGammaFracEta,fHistoTrueYieldEtaM02,1.,1.,"B");
    fHistoRatioMergedOneGammaFracEta->Scale(100.);
    fHistoRatioMergedOneElectronFracEta             = (TH1D*)fHistoTrueYieldMergedOneElectronFromEtaM02->Clone("RatioEtaMergedOneElectron");
    fHistoRatioMergedOneElectronFracEta->Divide(fHistoRatioMergedOneElectronFracEta,fHistoTrueYieldEtaM02,1.,1.,"B");
    fHistoRatioMergedOneElectronFracEta->Scale(100.);

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
//***************** Calculation of Meson Efficiency **************************
//****************************************************************************
TH1D* CalculateMesonEfficiencySec(TH1D* fMC_fMesonYieldsPt, TH1D* fMCYieldWithinAcc ,TString nameEfi ) {
    TH1D* fHistoMCMesonEffiPt   = (TH1D*)fMC_fMesonYieldsPt->Clone(nameEfi.Data());
    fHistoMCMesonEffiPt->Sumw2();
    fHistoMCMesonEffiPt->Divide(fHistoMCMesonEffiPt,fMCYieldWithinAcc,1.,1.,"B");
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

    TH1D*   fHistoClustersEPerEvent   = NULL;
    if (fHistoClustersE){
        fHistoClustersEPerEvent    = (TH1D*)fHistoClustersE->Rebin(fNBinsClusterPt,"ClusterEPerEvent",fBinsClusterPt);
        fHistoClustersEPerEvent->Divide(fDeltaPtCluster);
        fHistoClustersEPerEvent->Scale(1./fNEvents);
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
        if (fHistoClustersE)                            fHistoClustersE->Write("ClusterE");
        if (fHistoClustersEPerEvent)                    fHistoClustersEPerEvent->Write("ClusterEPerEvent");
        if (fDeltaPtCluster)                            fDeltaPtCluster->Write();

        if (fHistoClustersOverlapHeadersPt)             fHistoClustersOverlapHeadersPt->Write("ClusterOverlapHeadersPt");
        if (fHistoClustersOverlapHeadersPtPerEvent)     fHistoClustersOverlapHeadersPtPerEvent->Write("ClusterOverlapHeadersPtPerEvent");

        if(fHistoClustersMergedPtM02)                   fHistoClustersMergedPtM02->Write("ClusterMergedPtM02");
        if(fHistoClustersMergedPtM02AccMeson)           fHistoClustersMergedPtM02AccMeson->Write("ClusterMergedPtM02AccMeson");
        if(fHistoClustersMergedEM02AccMeson)            fHistoClustersMergedEM02AccMeson->Write("ClusterMergedEM02AccMeson");
        if(fHistoTrackVsClusterCandidates)              fHistoTrackVsClusterCandidates->Write("TrackCandidatesVsClusters");

        if (fAdvancedClusterQA){
            if (fHistoClusterNCellsPt)                  fHistoClusterNCellsPt->Write("ClusterNCellsPt");
            if (fHistoClusterMergedNCellsPt)            fHistoClusterMergedNCellsPt->Write("ClusterMergedNCellsPt");
            if (fHistoClusterMergedNCellsArClPt)        fHistoClusterMergedNCellsArClPt->Write("ClusterMergedNCellsAroundClusterPt");
            if (fHistoClusterMergedNCellsArAInclPt)     fHistoClusterMergedNCellsArAInclPt->Write("ClusterMergedNCellsAroundAndInClusterPt");
            if (fHistoClusterMergedEAroundClE)          fHistoClusterMergedEAroundClE->Write("ClusterMergedEAroundClusterE");
            if (fHistoClusterMergedNPartPt)             fHistoClusterMergedNPartPt->Write("ClusterMergedNPartPointingToClusPt");
            if (fHistoClustersMergedNLMPt)              fHistoClustersMergedNLMPt->Write("ClusterMergedNLMPt");
            for (Int_t i = 0; i< 11; i++ ){
                if (fHistoClustersMergedPtNLMBins[i]){
                     if (fHistoClustersMergedPtNLMBins[i]->GetEntries() > 0){
                         fHistoClustersMergedPtNLMBins[i]->Write();
                     }
                }
            }
        }

        for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
            if (fHistoInvMassPtBin[ii])                     fHistoInvMassPtBin[ii]->Write(Form("InvMass_PtBin%i",ii));
            if (fHistoM02PtBin[ii])                         fHistoM02PtBin[ii]->Write(Form("M02_PtBin%i",ii));
            if (fHistoExoticsM02EBin[ii])                   fHistoExoticsM02EBin[ii]->Write(Form("M02_ExoticsRejected_EBin%i",ii));
            if (fIsMC){
                if (fHistoTrueClusMergedM02PtBin[ii])       fHistoTrueClusMergedM02PtBin[ii]->Write(Form("ValidatedMerged_M02_PtBin%i",ii));
            }
        }
        cout << "End writing Uncorrected File" << endl;

        if (fDeltaPt)                                   fDeltaPt->Write();
        if (fHistoYieldMesonM02)                        fHistoYieldMesonM02->Write();

        if (fHaveToyMCInputForSec){
            cout << "Writing ToyMC input for secondary pi0 correction" << endl;
            for (Int_t j = 0; j < 3; j++){
                if (fHistoYieldToyMCSecInput[j])        fHistoYieldToyMCSecInput[j]->Write();
                if (fHistoYieldToyMCSecInputReb[j])     fHistoYieldToyMCSecInputReb[j]->Write();
            }
        }
        if (fHaveCocktailInputForSec){
            cout << "Writing input for secondary pi0 correction" << endl;
            for (Int_t j = 0; j < 3; j++){
                if (fHistoYieldExternSecInput[j]){
                    cout << "writing external sec input: " << nameSecondaries[j].Data() << endl;
                    fHistoYieldExternSecInput[j]->Write();
                }
                if (fHistoYieldExternSecInputReb[j])     fHistoYieldExternSecInputReb[j]->Write();
            }
        }
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
        if (fHistoTrueYieldMergedPartConvFromPi0M02)    fHistoTrueYieldMergedPartConvFromPi0M02->Write();
        if (fHistoTrueYieldMergedPartConvFromEtaM02)    fHistoTrueYieldMergedPartConvFromEtaM02->Write();
        if (fHistoTrueYieldMergedPureM02)           fHistoTrueYieldMergedPureM02->Write();
        if (fHistoTrueYieldMergedPureFromPi0M02)    fHistoTrueYieldMergedPureFromPi0M02->Write();
        if (fHistoTrueYieldMergedPureFromEtaM02)    fHistoTrueYieldMergedPureFromEtaM02->Write();
        if (fHistoTrueYieldMergedOneGammaM02)       fHistoTrueYieldMergedOneGammaM02->Write();
        if (fHistoTrueYieldMergedOneGammaFromPi0M02)    fHistoTrueYieldMergedOneGammaFromPi0M02->Write();
        if (fHistoTrueYieldMergedOneGammaFromEtaM02)    fHistoTrueYieldMergedOneGammaFromEtaM02->Write();
        if (fHistoTrueYieldMergedOneElectronM02)    fHistoTrueYieldMergedOneElectronM02->Write();
        if (fHistoTrueYieldMergedOneElectronFromPi0M02) fHistoTrueYieldMergedOneElectronFromPi0M02->Write();
        if (fHistoTrueYieldMergedOneElectronFromEtaM02) fHistoTrueYieldMergedOneElectronFromEtaM02->Write();
        if (fHistoTrueYieldPi0M02)                  fHistoTrueYieldPi0M02->Write();
        if (fHistoTrueYieldPi0DCM02)                fHistoTrueYieldPi0DCM02->Write();
        if (fHistoTrueYieldPi0GGM02)                fHistoTrueYieldPi0GGM02->Write();
        if (fHistoTrueYieldPi0DalitzM02)            fHistoTrueYieldPi0DalitzM02->Write();
        if (fHistoTrueYieldEtaM02)                  fHistoTrueYieldEtaM02->Write();
        if (fHistoTrueYieldEtaDCM02)                fHistoTrueYieldEtaDCM02->Write();
        if (fHistoTrueYieldEtaGGM02)                fHistoTrueYieldEtaGGM02->Write();
        if (fHistoTrueYieldEtaDalitzM02)            fHistoTrueYieldEtaDalitzM02->Write();
        if (fHistoTrueYieldGammaM02)                fHistoTrueYieldGammaM02->Write();
        if (fHistoTrueYieldElectronM02)             fHistoTrueYieldElectronM02->Write();
        if (fHistoTrueYieldBGM02)                   fHistoTrueYieldBGM02->Write();
        if (fHistoTrueYieldPrimPi0M02)              fHistoTrueYieldPrimPi0M02->Write();
        if (fHistoTrueClustersBGPtSource)           fHistoTrueClustersBGPtSource->Write();
        if (fHistoTrueClustersGammaPtSource)        fHistoTrueClustersGammaPtSource->Write();
        if (fHistoTrueClustersElectronPtSource)     fHistoTrueClustersElectronPtSource->Write();

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

        if (fHistoRatioTrueYieldEtaM02)             fHistoRatioTrueYieldEtaM02->Write();
        if (fHistoRatioTrueYieldPi0M02)             fHistoRatioTrueYieldPi0M02->Write();
        if (fHistoRatioTrueYieldGammaM02)           fHistoRatioTrueYieldGammaM02->Write();
        if (fHistoRatioTrueYieldElectronM02)        fHistoRatioTrueYieldElectronM02->Write();
        if (fHistoRatioPi0DCFrac)                   fHistoRatioPi0DCFrac->Write();
        if (fHistoRatioPi0GGFrac)                   fHistoRatioPi0GGFrac->Write();
        if (fHistoRatioPi0DalitzFrac)               fHistoRatioPi0DalitzFrac->Write();
        if (fHistoRatioEtaDCFrac)                   fHistoRatioEtaDCFrac->Write();
        if (fHistoRatioEtaGGFrac)                   fHistoRatioEtaGGFrac->Write();
        if (fHistoRatioEtaDalitzFrac)               fHistoRatioEtaDalitzFrac->Write();
        if (fHistoRatioMergedPureFracPi0)           fHistoRatioMergedPureFracPi0->Write();
        if (fHistoRatioMergedPartConvFracPi0)       fHistoRatioMergedPartConvFracPi0->Write();
        if (fHistoRatioMergedOneGammaFracPi0)       fHistoRatioMergedOneGammaFracPi0->Write();
        if (fHistoRatioMergedOneElectronFracPi0)    fHistoRatioMergedOneElectronFracPi0->Write();
        if (fHistoRatioMergedPureFracEta)           fHistoRatioMergedPureFracEta->Write();
        if (fHistoRatioMergedPartConvFracEta)       fHistoRatioMergedPartConvFracEta->Write();
        if (fHistoRatioMergedOneGammaFracEta)       fHistoRatioMergedOneGammaFracEta->Write();
        if (fHistoRatioMergedOneElectronFracEta)    fHistoRatioMergedOneElectronFracEta->Write();

        for (Int_t i = 0; i< 11; i++ ){
            if (fHistoClustersMergedPtNLMBins[i]){
                if (fHistoClustersMergedPtNLMBins[i]->GetEntries() > 0){
                    fHistoClustersMergedPtNLMBins[i]->Write();
                }
            }
        }


        for (Int_t i = 0; i< 10; i++){
            if (fHistoTrueClustersBGPt[i])          fHistoTrueClustersBGPt[i]->Write();
            if (fHistoRatioTrueClustersBGPt[i])     fHistoRatioTrueClustersBGPt[i]->Write();
        }
        for (Int_t i = 0; i< 8; i++){
            if (fHistoTrueClustersGammaPt[i])       fHistoTrueClustersGammaPt[i]->Write();
        }
        for (Int_t i = 0; i< 9; i++){
            if (fHistoTrueClustersElectronPt[i])    fHistoTrueClustersElectronPt[i]->Write();
        }

        for (Int_t j = 0; j< 4; j++){
            if (fHistoTrueYieldSecPi0M02[j])        fHistoTrueYieldSecPi0M02[j]->Write();
            if (fHistoTruePi0SecFrac[j])            fHistoTruePi0SecFrac[j]->Write();
            if (fHistoTrueEffiSecPi0[j])            fHistoTrueEffiSecPi0[j]->Write();
        }

        if (fNewMCOutput){
            if (fHistoMCSecPi0PtSource)                 fHistoMCSecPi0PtSource->Write("MC_SecPi0_Source");
            if (fHistoMCSecPi0WithinAccepPtSource)      fHistoMCSecPi0WithinAccepPtSource->Write("MC_SecPi0_WithinAccept_Source");
            for (Int_t j = 0; j < 4; j++){
                if(fHistoMCSecPi0Pt[j])         fHistoMCSecPi0Pt[j]->Write();
                if(fHistoMCSecPi0PtWAcc[j])     fHistoMCSecPi0PtWAcc[j]->Write();
                if(fHistoMCSecPi0PtReb[j])      fHistoMCSecPi0PtReb[j]->Write();
                if(fHistoMCSecPi0PtWAccReb[j])  fHistoMCSecPi0PtWAccReb[j]->Write();
                if(fHistoMCSecPi0AcceptPt[j])   fHistoMCSecPi0AcceptPt[j]->Write();
            }
        }
        for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
            if(!fDoExtendedExamplePlots && ii!=fExampleBin)
                continue;
            TString titleHistoSignal    = Form("%3.1f GeV/#it{c} < #it{p}_{T} < %3.1f GeV/#it{c}",fBinsPt[ii],fBinsPt[ii+1]);
            if (fHistoM02PtBin[ii]){
                fHistoM02PtBin[ii]->SetTitle(titleHistoSignal.Data());
                fHistoM02PtBin[ii]->Write();
            }
            if (fHistoTrueClusMergedM02PtBin[ii])                  fHistoTrueClusMergedM02PtBin[ii]->Write();
            if (fHistoTrueClusPi0M02PtBin[ii])                     fHistoTrueClusPi0M02PtBin[ii]->Write();
            if (fHistoTrueClusPi0DCM02PtBin[ii])                   fHistoTrueClusPi0DCM02PtBin[ii]->Write();
            if (fHistoTrueClusPi0GGM02PtBin[ii])                   fHistoTrueClusPi0GGM02PtBin[ii]->Write();
            if (fHistoTrueClusPi0DalitzM02PtBin[ii])               fHistoTrueClusPi0DalitzM02PtBin[ii]->Write();
            if (fHistoTrueClusEtaM02PtBin[ii])                     fHistoTrueClusEtaM02PtBin[ii]->Write();
            if (fHistoTrueClusEtaDCM02PtBin[ii])                   fHistoTrueClusEtaDCM02PtBin[ii]->Write();
            if (fHistoTrueClusEtaGGM02PtBin[ii])                   fHistoTrueClusEtaGGM02PtBin[ii]->Write();
            if (fHistoTrueClusEtaDalitzM02PtBin[ii])               fHistoTrueClusEtaDalitzM02PtBin[ii]->Write();
            if (fHistoTrueClusGammaM02PtBin[ii])                   fHistoTrueClusGammaM02PtBin[ii]->Write();
            if (fHistoTrueClusElectronM02PtBin[ii])                fHistoTrueClusElectronM02PtBin[ii]->Write();
            if (fHistoTrueClusBGM02PtBin[ii])                      fHistoTrueClusBGM02PtBin[ii]->Write();
            if (fHistoTrueClusPartConvMergedM02PtBin[ii])          fHistoTrueClusPartConvMergedM02PtBin[ii]->Write();
            if (fHistoTrueClusPartConvMergedFromPi0M02PtBin[ii])   fHistoTrueClusPartConvMergedFromPi0M02PtBin[ii]->Write();
            if (fHistoTrueClusPartConvMergedFromEtaM02PtBin[ii])   fHistoTrueClusPartConvMergedFromEtaM02PtBin[ii]->Write();
            if (fHistoTrueClusPureMergedM02PtBin[ii])              fHistoTrueClusPureMergedM02PtBin[ii]->Write();
            if (fHistoTrueClusPureMergedFromPi0M02PtBin[ii])       fHistoTrueClusPureMergedFromPi0M02PtBin[ii]->Write();
            if (fHistoTrueClusPureMergedFromEtaM02PtBin[ii])       fHistoTrueClusPureMergedFromEtaM02PtBin[ii]->Write();
            if (fHistoTrueClusOneGammaM02PtBin[ii])                fHistoTrueClusOneGammaM02PtBin[ii]->Write();
            if (fHistoTrueClusOneGammaFromPi0M02PtBin[ii])         fHistoTrueClusOneGammaFromPi0M02PtBin[ii]->Write();
            if (fHistoTrueClusOneGammaFromEtaM02PtBin[ii])         fHistoTrueClusOneGammaFromEtaM02PtBin[ii]->Write();
            if (fHistoTrueClusOneElectronM02PtBin[ii])             fHistoTrueClusOneElectronM02PtBin[ii]->Write();
            if (fHistoTrueClusOneElectronFromPi0M02PtBin[ii])      fHistoTrueClusOneElectronFromPi0M02PtBin[ii]->Write();
            if (fHistoTrueClusOneElectronFromEtaM02PtBin[ii])      fHistoTrueClusOneElectronFromEtaM02PtBin[ii]->Write();
            if (fHistoTrueClusPrimPi0M02PtBin[ii])                 fHistoTrueClusPrimPi0M02PtBin[ii]->Write();
            for (Int_t j = 0; j<4; j++){
                if (fHistoTrueClusSecPi0M02PtBin[j][ii])           fHistoTrueClusSecPi0M02PtBin[j][ii]->Write();
            }
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
    if (fHistoClustersE)                                        delete fHistoClustersE;
    if (fHistoClustersOverlapHeadersPt)                         delete fHistoClustersOverlapHeadersPt;
    if (fHistoClustersMergedPtM02)                              delete fHistoClustersMergedPtM02;
    if (fHistoClustersMergedPtM02AccMeson)                      delete fHistoClustersMergedPtM02AccMeson;
    if (fHistoClustersMergedEM02AccMeson)                       delete fHistoClustersMergedEM02AccMeson;
    if (fHistoTrackVsClusterCandidates)                         delete fHistoTrackVsClusterCandidates;
    if (fHistoInvMassVsPt)                                      delete fHistoInvMassVsPt;
    if (fHistoSPDtrackletvsSPDclusters)                         delete fHistoSPDtrackletvsSPDclusters;
    if (fMesonM02Yields)                                        delete fMesonM02Yields;
    if (fMesonM02YieldsError)                                   delete fMesonM02YieldsError;

    if (fIsMC){
        if (fHistoTrueClustersMergedPtM02)                      delete fHistoTrueClustersMergedPtM02;
        if (fHistoTrueClustersPi0PtM02)                         delete fHistoTrueClustersPi0PtM02;
        if (fHistoTrueClustersPi0DCPtM02)                       delete fHistoTrueClustersPi0DCPtM02;
        if (fHistoTrueClustersPi0GGPtM02)                       delete fHistoTrueClustersPi0GGPtM02;
        if (fHistoTrueClustersPi0DalitzPtM02)                   delete fHistoTrueClustersPi0DalitzPtM02;
        if (fHistoTrueClustersEtaPtM02)                         delete fHistoTrueClustersEtaPtM02;
        if (fHistoTrueClustersEtaDCPtM02)                       delete fHistoTrueClustersEtaDCPtM02;
        if (fHistoTrueClustersEtaGGPtM02)                       delete fHistoTrueClustersEtaGGPtM02;
        if (fHistoTrueClustersEtaDalitzPtM02)                   delete fHistoTrueClustersEtaDalitzPtM02;
        if (fHistoTrueClustersGammaPtM02)                       delete fHistoTrueClustersGammaPtM02;
        if (fHistoTrueClustersElectronPtM02)                    delete fHistoTrueClustersElectronPtM02;
        if (fHistoTrueClustersBGPtM02)                          delete fHistoTrueClustersBGPtM02;
        if (fHistoTrueClusPartConvMergedPtM02)                  delete fHistoTrueClusPartConvMergedPtM02;
        if (fHistoTrueClustersBGPtSource)                       delete fHistoTrueClustersBGPtSource;
        if (fHistoTrueClustersGammaPtSource)                    delete fHistoTrueClustersGammaPtSource;
        if (fHistoTrueClustersElectronPtSource)                 delete fHistoTrueClustersElectronPtSource;
        if (fMesonM02TrueMergedYields)                          delete fMesonM02TrueMergedYields;
        if (fMesonM02TruePi0Yields)                             delete fMesonM02TruePi0Yields;
        if (fMesonM02TruePi0DCYields)                           delete fMesonM02TruePi0DCYields;
        if (fMesonM02TruePi0GGYields)                           delete fMesonM02TruePi0GGYields;
        if (fMesonM02TruePi0DalitzYields)                       delete fMesonM02TruePi0DalitzYields;
        if (fMesonM02TrueEtaYields)                             delete fMesonM02TrueEtaYields;
        if (fMesonM02TrueEtaDCYields)                           delete fMesonM02TrueEtaDCYields;
        if (fMesonM02TrueEtaGGYields)                           delete fMesonM02TrueEtaGGYields;
        if (fMesonM02TrueEtaDalitzYields)                       delete fMesonM02TrueEtaDalitzYields;
        if (fMesonM02TrueGammaYields)                           delete fMesonM02TrueGammaYields;
        if (fMesonM02TrueElectronYields)                        delete fMesonM02TrueElectronYields;
        if (fMesonM02TrueBGYields)                              delete fMesonM02TrueBGYields;
        if (fMesonM02TrueMergedYieldsError)                     delete fMesonM02TrueMergedYieldsError;
        if (fMesonM02TruePi0YieldsError)                        delete fMesonM02TruePi0YieldsError;
        if (fMesonM02TruePi0DCYieldsError)                      delete fMesonM02TruePi0DCYieldsError;
        if (fMesonM02TruePi0GGYieldsError)                      delete fMesonM02TruePi0GGYieldsError;
        if (fMesonM02TruePi0DalitzYieldsError)                  delete fMesonM02TruePi0DalitzYieldsError;
        if (fMesonM02TrueEtaYieldsError)                        delete fMesonM02TrueEtaYieldsError;
        if (fMesonM02TrueEtaDCYieldsError)                      delete fMesonM02TrueEtaDCYieldsError;
        if (fMesonM02TrueEtaGGYieldsError)                      delete fMesonM02TrueEtaGGYieldsError;
        if (fMesonM02TrueEtaDalitzYieldsError)                  delete fMesonM02TrueEtaDalitzYieldsError;
        if (fMesonM02TrueGammaYieldsError)                      delete fMesonM02TrueGammaYieldsError;
        if (fMesonM02TrueElectronYieldsError)                   delete fMesonM02TrueElectronYieldsError;
        if (fMesonM02TrueBGYieldsError)                         delete fMesonM02TrueBGYieldsError;
        if (fMesonM02TrueMergedPartConvYields)                  delete fMesonM02TrueMergedPartConvYields;
        if (fMesonM02TrueMergedPartConvYieldsError)             delete fMesonM02TrueMergedPartConvYieldsError;
        if (fMesonM02TrueMergedPartConvFromPi0Yields)           delete fMesonM02TrueMergedPartConvFromPi0Yields;
        if (fMesonM02TrueMergedPartConvFromPi0YieldsError)      delete fMesonM02TrueMergedPartConvFromPi0YieldsError;
        if (fMesonM02TrueMergedPartConvFromEtaYields)           delete fMesonM02TrueMergedPartConvFromEtaYields;
        if (fMesonM02TrueMergedPartConvFromEtaYieldsError)      delete fMesonM02TrueMergedPartConvFromEtaYieldsError;
        if (fMesonM02TrueMergedPureYields)                      delete fMesonM02TrueMergedPureYields;
        if (fMesonM02TrueMergedPureYieldsError)                 delete fMesonM02TrueMergedPureYieldsError;
        if (fMesonM02TrueMergedPureFromPi0Yields)               delete fMesonM02TrueMergedPureFromPi0Yields;
        if (fMesonM02TrueMergedPureFromPi0YieldsError)          delete fMesonM02TrueMergedPureFromPi0YieldsError;
        if (fMesonM02TrueMergedPureFromEtaYields)               delete fMesonM02TrueMergedPureFromEtaYields;
        if (fMesonM02TrueMergedPureFromEtaYieldsError)          delete fMesonM02TrueMergedPureFromEtaYieldsError;
        if (fMesonM02TrueMergedOneGammaYields)                  delete fMesonM02TrueMergedOneGammaYields;
        if (fMesonM02TrueMergedOneGammaYieldsError)             delete fMesonM02TrueMergedOneGammaYieldsError;
        if (fMesonM02TrueMergedOneGammaFromPi0Yields)           delete fMesonM02TrueMergedOneGammaFromPi0Yields;
        if (fMesonM02TrueMergedOneGammaFromPi0YieldsError)      delete fMesonM02TrueMergedOneGammaFromPi0YieldsError;
        if (fMesonM02TrueMergedOneGammaFromEtaYields)           delete fMesonM02TrueMergedOneGammaFromEtaYields;
        if (fMesonM02TrueMergedOneGammaFromEtaYieldsError)      delete fMesonM02TrueMergedOneGammaFromEtaYieldsError;
        if (fMesonM02TrueMergedOneElectronYields)               delete fMesonM02TrueMergedOneElectronYields;
        if (fMesonM02TrueMergedOneElectronYieldsError)          delete fMesonM02TrueMergedOneElectronYieldsError;
        if (fMesonM02TrueMergedOneElectronFromPi0Yields)        delete fMesonM02TrueMergedOneElectronFromPi0Yields;
        if (fMesonM02TrueMergedOneElectronFromPi0YieldsError)   delete fMesonM02TrueMergedOneElectronFromPi0YieldsError;
        if (fMesonM02TrueMergedOneElectronFromEtaYields)        delete fMesonM02TrueMergedOneElectronFromEtaYields;
        if (fMesonM02TrueMergedOneElectronFromEtaYieldsError)   delete fMesonM02TrueMergedOneElectronFromEtaYieldsError;
        if (fMesonM02TruePrimPi0Yields)                         delete fMesonM02TruePrimPi0Yields;
        if (fMesonM02TruePrimPi0YieldsError)                    delete fMesonM02TruePrimPi0YieldsError;
        for (Int_t j = 0; j<4; j++){
            if (fMesonM02TrueSecPi0Yields[j])                   delete fMesonM02TrueSecPi0Yields[j];
            if (fMesonM02TrueSecPi0YieldsError[j])              delete fMesonM02TrueSecPi0YieldsError[j];
        }
    }

    for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
        if (fHistoInvMassPtBin[ii])                             delete fHistoInvMassPtBin[ii];
        if (fHistoM02PtBin[ii])                                 delete fHistoM02PtBin[ii];
        if (fHistoExoticsM02EBin[ii])                           delete fHistoExoticsM02EBin[ii];
        if (fIsMC){
            if (fHistoTrueClusMergedM02PtBin[ii])                   delete fHistoTrueClusMergedM02PtBin[ii];
            if (fHistoTrueClusPi0M02PtBin[ii])                      delete fHistoTrueClusPi0M02PtBin[ii];
            if (fHistoTrueClusPi0DCM02PtBin[ii])                    delete fHistoTrueClusPi0DCM02PtBin[ii];
            if (fHistoTrueClusPi0GGM02PtBin[ii])                    delete fHistoTrueClusPi0GGM02PtBin[ii];
            if (fHistoTrueClusPi0DalitzM02PtBin[ii])                delete fHistoTrueClusPi0DalitzM02PtBin[ii];
            if (fHistoTrueClusEtaM02PtBin[ii])                      delete fHistoTrueClusEtaM02PtBin[ii];
            if (fHistoTrueClusEtaDCM02PtBin[ii])                    delete fHistoTrueClusEtaDCM02PtBin[ii];
            if (fHistoTrueClusEtaGGM02PtBin[ii])                    delete fHistoTrueClusEtaGGM02PtBin[ii];
            if (fHistoTrueClusEtaDalitzM02PtBin[ii])                delete fHistoTrueClusEtaDalitzM02PtBin[ii];
            if (fHistoTrueClusGammaM02PtBin[ii])                    delete fHistoTrueClusGammaM02PtBin[ii];
            if (fHistoTrueClusElectronM02PtBin[ii])                 delete fHistoTrueClusElectronM02PtBin[ii];
            if (fHistoTrueClusBGM02PtBin[ii])                       delete fHistoTrueClusBGM02PtBin[ii];
            if (fHistoTrueClusPartConvMergedM02PtBin[ii])           delete fHistoTrueClusPartConvMergedM02PtBin[ii];
            if (fHistoTrueClusPartConvMergedFromPi0M02PtBin[ii])    delete fHistoTrueClusPartConvMergedFromPi0M02PtBin[ii];
            if (fHistoTrueClusPartConvMergedFromEtaM02PtBin[ii])    delete fHistoTrueClusPartConvMergedFromEtaM02PtBin[ii];
            if (fHistoTrueClusPureMergedM02PtBin[ii])               delete fHistoTrueClusPureMergedM02PtBin[ii];
            if (fHistoTrueClusPureMergedFromPi0M02PtBin[ii])        delete fHistoTrueClusPureMergedFromPi0M02PtBin[ii];
            if (fHistoTrueClusPureMergedFromEtaM02PtBin[ii])        delete fHistoTrueClusPureMergedFromEtaM02PtBin[ii];
            if (fHistoTrueClusOneGammaM02PtBin[ii])                 delete fHistoTrueClusOneGammaM02PtBin[ii];
            if (fHistoTrueClusOneGammaFromPi0M02PtBin[ii])          delete fHistoTrueClusOneGammaFromPi0M02PtBin[ii];
            if (fHistoTrueClusOneGammaFromEtaM02PtBin[ii])          delete fHistoTrueClusOneGammaFromEtaM02PtBin[ii];
            if (fHistoTrueClusOneElectronM02PtBin[ii])              delete fHistoTrueClusOneElectronM02PtBin[ii];
            if (fHistoTrueClusOneElectronFromPi0M02PtBin[ii])       delete fHistoTrueClusOneElectronFromPi0M02PtBin[ii];
            if (fHistoTrueClusOneElectronFromEtaM02PtBin[ii])       delete fHistoTrueClusOneElectronFromEtaM02PtBin[ii];
            if (fHistoTrueClusPrimPi0M02PtBin[ii])                  delete fHistoTrueClusPrimPi0M02PtBin[ii];
            for (Int_t j = 0; j<4; j++){
                if (fHistoTrueClusSecPi0M02PtBin[j][ii])            delete fHistoTrueClusSecPi0M02PtBin[j][ii];
            }
        }
    }

    if (fDeltaPt)                                               delete fDeltaPt;
    if (fHistoYieldMesonM02)                                    delete fHistoYieldMesonM02;
    if (fIsMC){
        if (fHistoTrueYieldMergedM02)                           delete fHistoTrueYieldMergedM02;
        if (fHistoTrueYieldMergedPartConvM02)                   delete fHistoTrueYieldMergedPartConvM02;
        if (fHistoTrueYieldMergedPartConvFromPi0M02)            delete fHistoTrueYieldMergedPartConvFromPi0M02;
        if (fHistoTrueYieldMergedPartConvFromEtaM02)            delete fHistoTrueYieldMergedPartConvFromEtaM02;
        if (fHistoTrueYieldMergedPureM02)                       delete fHistoTrueYieldMergedPureM02;
        if (fHistoTrueYieldMergedPureFromPi0M02)                delete fHistoTrueYieldMergedPureFromPi0M02;
        if (fHistoTrueYieldMergedPureFromEtaM02)                delete fHistoTrueYieldMergedPureFromEtaM02;
        if (fHistoTrueYieldMergedOneGammaM02)                   delete fHistoTrueYieldMergedOneGammaM02;
        if (fHistoTrueYieldMergedOneGammaFromPi0M02)            delete fHistoTrueYieldMergedOneGammaFromPi0M02;
        if (fHistoTrueYieldMergedOneGammaFromEtaM02)            delete fHistoTrueYieldMergedOneGammaFromEtaM02;
        if (fHistoTrueYieldMergedOneElectronM02)                delete fHistoTrueYieldMergedOneElectronM02;
        if (fHistoTrueYieldMergedOneElectronFromPi0M02)         delete fHistoTrueYieldMergedOneElectronFromPi0M02;
        if (fHistoTrueYieldMergedOneElectronFromEtaM02)         delete fHistoTrueYieldMergedOneElectronFromEtaM02;
        if (fHistoTrueYieldPi0M02)                              delete fHistoTrueYieldPi0M02;
        if (fHistoTrueYieldPi0DCM02)                            delete fHistoTrueYieldPi0DCM02;
        if (fHistoTrueYieldPi0GGM02)                            delete fHistoTrueYieldPi0GGM02;
        if (fHistoTrueYieldPi0DalitzM02)                        delete fHistoTrueYieldPi0DalitzM02;
        if (fHistoTrueYieldEtaM02)                              delete fHistoTrueYieldEtaM02;
        if (fHistoTrueYieldEtaDCM02)                            delete fHistoTrueYieldEtaDCM02;
        if (fHistoTrueYieldEtaGGM02)                            delete fHistoTrueYieldEtaGGM02;
        if (fHistoTrueYieldEtaDalitzM02)                        delete fHistoTrueYieldEtaDalitzM02;
        if (fHistoTrueYieldGammaM02)                            delete fHistoTrueYieldGammaM02;
        if (fHistoTrueYieldElectronM02)                         delete fHistoTrueYieldElectronM02;
        if (fHistoTrueYieldBGM02)                               delete fHistoTrueYieldBGM02;
        if (fHistoTrueYieldPrimPi0M02)                          delete fHistoTrueYieldPrimPi0M02;

        if (fHistoMCAcceptancePt)                               delete fHistoMCAcceptancePt;
        if (fHistoTruePurityMerged)                             delete fHistoTruePurityMerged;
        if (fHistoTruePi0PurityMerged)                          delete fHistoTruePi0PurityMerged;
        if (fHistoTrueEtaPurityMerged)                          delete fHistoTrueEtaPurityMerged;
        if (fHistoTrueEffiPrimMeson)                            delete fHistoTrueEffiPrimMeson;
        for (Int_t j = 0; j<4; j++){
            if (fHistoTrueYieldSecPi0M02[j])                    delete fHistoTrueYieldSecPi0M02[j];
            if (fHistoTrueEffiSecPi0[j])                        delete fHistoTrueEffiSecPi0[j];
            if (fHistoTruePi0SecFrac[j])                        delete fHistoTruePi0SecFrac[j];
        }

    }

    // Delete Toy MC files
    for (Int_t j = 0; j < 3; j++){
        if (fFileToyMCInput[j] )                                   delete fFileToyMCInput[j];
    }
    if (fFileCocktailInput)                                     delete fFileCocktailInput;

}
