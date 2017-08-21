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
#include "ExtractSignalPiPlPiMiPiZero.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctions.h" // changed to standard CommonHeaders
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ExtractSignalPlotting.h"
#include "THnSparse.h"

Double_t FunctionBGExclusion(Double_t *x, Double_t *par){
    if (x[0] > fBGFitRangeLeft[1] && x[0] < fBGFitRange[0]) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
}


// Main Function
void ExtractSignalPiPlPiMiPiZero(   TString meson                  = "",
                                    TString file                   = "",
                                    TString cutSelection           = "",
                                    TString Suffix                 = "",
                                    TString optionMC               = "",
                                    TString optionEnergy           = "",
                                    TString optionCrystalBall      = "",
                                    TString optionUseMinBiasEff    = "",
                                    TString optionPeriod           = "",
                                    TString optionAdvancedMesonQA  = "",
                                    Int_t numberOfBins             = 30,
                                    Bool_t addSig                  = kFALSE,
                                    Int_t UsrMode                  = 40         // Mode given by user
                                ){
    gROOT->Reset();

    //------------NEW LABELING ---------------
    // mode:	40 // new output PCM-PCM
    //			41 // new output PCM EMCAL
    //			42 // new output PCM-PHOS
    //			43 // new output PCM-DCAL
    //			44 // new output EMCAL-EMCAL
    //			45 // new output PHOS-PHOS
    //			46 // old output DCAL-DCAL
    //      	47 // new output PCM-DALITZ
    //			48 // new output EMCAL-DALITZ
    //			49 // new output PHOS-DALITZ
    //			50 // new output DCAL-DALITZ

    fCutSelection               = cutSelection;
    TString fCutSelectionRead   = cutSelection;

    //Int_t fMode = -1;
    fMode                       = ReturnSeparatedCutNumberPiPlPiMiPiZero( cutSelection, fTypeCutSelection, fEventCutSelection, fGammaCutSelection, fClusterCutSelection,
                                                                          fPionCutSelection, fNeutralPionCutSelection, fMesonCutSelection);
    if(fMode!=UsrMode){
        cout << "ERROR: Chosen mode is not identical to mode extracted from cutstring! " << endl;
        return;
    } // check if extracted mode from cutnumber is the same mode given by user
    Int_t mode                  = fMode;


    cout << "\t MODE = " << mode << endl;

    StyleSettingsThesis();
    SetPlotStyle();

    fEnergyFlag             = optionEnergy;
    fPrefix                 = meson;

    fPeriodFlag             = optionPeriod;

    TString outputDir       = Form("%s/%s/%s/ExtractSignal",cutSelection.Data(),optionEnergy.Data(),Suffix.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    //gSystem->Exec("mkdir -p "+outputDir+"/Plots_Data_QA");
    //gSystem->Exec("mkdir -p "+outputDir+"/Plots_MC_QA");

    cout<<"Pictures are saved as "<< Suffix.Data()<< endl;
    fdate = ReturnDateString();

    //****************************** Specification of collision system ************************************************
    TString textProcess = ReturnMesonString (fPrefix);
    if(textProcess.CompareTo("") == 0 ){
        cout << "Meson unknown" << endl;
        return ;
    }

    fTextMeasurement = Form("%s #rightarrow #gamma#gamma", textProcess.Data());
    if (meson.CompareTo("Omega")==0) fTextMeasurement = "#omega #rightarrow #pi^{+} #pi^{-} #pi^{0}";
    if (meson.CompareTo("Eta")==0) fTextMeasurement = "#eta #rightarrow #pi^{+} #pi^{-} #pi^{0}";

    fCollisionSystem = ReturnFullCollisionsSystem(fEnergyFlag);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    fDetectionProcess = ReturnFullTextReconstructionProcess(mode);

    cout << "Detection process is: " << fDetectionProcess.Data() << endl;

    //****************************** Choice of Fitting procedure ******************************************************
    if(optionCrystalBall.CompareTo("CrystalBall") == 0){// means we want to plot values for the pi0
        fCrysFitting=1;
        cout << "CrystalBall fit chosen ..." << endl;
    } else	{
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
    if (centralityString.CompareTo("pp")!=0 && !centralityString.Contains("0-100%") ){ // if collision system if not pp, mention it in string
        fCollisionSystem    = Form("%s %s", centralityString.Data(), fCollisionSystem.Data());
    }

    //***************************** Initialization of variables according to meson type ******************************
    if (meson.CompareTo("Eta") == 0) {
        Initialize("Eta", numberOfBins);
    } else if (meson.CompareTo("Omega") == 0) {
        Initialize("Omega", numberOfBins);
    } else {
        cout<<"ERROR: First argument in the ExtractSignal(....) has to be either Eta or Omega"<<endl;
        return;
    }
    //************************* Start of Main routine ***************************************************************
    const char* fFileErrLogDatname = Form("%s/%s/%s_%s_FileErrLog%s_%s.dat",cutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelectionRead.Data());
    fFileErrLog.open(fFileErrLogDatname, ios::out);

    TFile* f = new TFile(file.Data());

    TString nameMainDir = "";
    //if (mode == 9 || mode == 0) {nameMainDir = "GammaConvNeutralMesonPiPlPiMiPiZero_0_9";}
    //else if (mode == 2 || mode == 3) nameMainDir = "GammaConvCalo";

    //detect subdir
    nameMainDir     = AutoDetectMainTList(40 , f); // TODO: Change 40 to mode if modes are changed in ConversionFunctions correctly
    if(nameMainDir.CompareTo("")!=0){
        cout << "Succesfully detected " << nameMainDir <<" as MainDir!" << endl;
    } else{
        printf("ERROR: MainDir not found!");
        return;
    }
    //if (mode == 6) nameMainDir = "GammaConvNeutralMesonPiPlPiMiPiZero_1";

    TList *TopDir =(TList*)f->Get(nameMainDir.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }

    TList *HistosGammaConversion = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if(HistosGammaConversion == NULL){
        cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in File"<<endl;
        return;
    }
    TList *ESDContainer             = (TList*) HistosGammaConversion->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));
    TList *BackgroundContainer      = (TList*) HistosGammaConversion->FindObject(Form("%s Back histograms",fCutSelectionRead.Data()));
    TList *MotherContainer          = (TList*) HistosGammaConversion->FindObject(Form("%s Mother histograms",fCutSelectionRead.Data()));
    fNumberOfGoodESDTracks          = (TH1D*)ESDContainer->FindObject("GoodESDTracks");
    fEventQuality                   = (TH1D*)ESDContainer->FindObject("NEvents");

    TString rapidityRange;
    fYMaxMeson                      =  ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
    fBackgroundMultNumber           = ReturnBackgroundMult(fMesonCutSelection);

    TString ObjectNameESD           = "ESD_Mother_InvMass_Pt";
    TString ObjectNameBck[4]        = { "ESD_Background_1_InvMass_Pt",
                                        "ESD_Background_2_InvMass_Pt",
                                        "ESD_Background_3_InvMass_Pt",
                                        "ESD_Background_4_InvMass_Pt"
                                      };
    // Get Names of alternative InvMass histos

    // InvMass minus Pi0 inv Mass
    TString ObjectNameESDSubPiZero           =   "ESD_InvMass_Mother_Sub_InvMass(NeutralPion)_Pt";
    TString ObjectNameBckSubPiZero[4]        = { "ESD_Background_1_InvMass_Sub_InvMass(NeutralPion)_Pt",
                                                 "ESD_Background_2_InvMass_Sub_InvMass(NeutralPion)_Pt",
                                                 "ESD_Background_3_InvMass_Sub_InvMass(NeutralPion)_Pt",
                                                 "ESD_Background_4_InvMass_Sub_InvMass(NeutralPion)_Pt"
                                      };

    // InvMass with fixed pi0 p_z
    TString ObjectNameESDFixedPzPiZero       =   "ESD_InvMass_Mother_FixedPz(NeutralPion)_Pt";
    TString ObjectNameBckFixedPzPiZero[4]    = { "ESD_Background_1_InvMass_FixedPz(NeutralPion)_Pt",
                                                 "ESD_Background_2_InvMass_FixedPz(NeutralPion)_Pt",
                                                 "ESD_Background_3_InvMass_FixedPz(NeutralPion)_Pt",
                                                 "ESD_Background_4_InvMass_FixedPz(NeutralPion)_Pt"
                                      };

    for(Int_t k=0;k<4;k++){
       hist_bck[k]                       = (TH2D*)ESDContainer->FindObject(ObjectNameBck[k].Data());
       hist_bck_SubPiZero[k]             = (TH2D*)ESDContainer->FindObject(ObjectNameBckSubPiZero[k].Data());
       hist_bck_FixedPzPiZero[k]         = (TH2D*)ESDContainer->FindObject(ObjectNameBckFixedPzPiZero[k].Data());
    }

    for(Int_t k=1;k<5;k++){
        fBckInvMassVSPt[k]                = (TH2D*)hist_bck[k-1]->Clone(Form("hist_bck%i",k));
        fBckInvMassVSPt_SubPiZero[k]      = (TH2D*)hist_bck_SubPiZero[k-1]->Clone(Form("hist_bck%i_SubPiZero",k));
        fBckInvMassVSPt_FixedPzPiZero[k]  = (TH2D*)hist_bck_FixedPzPiZero[k-1]->Clone(Form("hist_bck%i_FixedPzPiZero",k));
    }

    SetCorrectMCHistogrammNames();

    fGammaGammaInvMassVSPt                = (TH2D*)ESDContainer->FindObject(ObjectNameESD.Data());
    fGammaGammaInvMassVSPt_SubPiZero      = (TH2D*)ESDContainer->FindObject(ObjectNameESDSubPiZero.Data());
    fGammaGammaInvMassVSPt_FixedPzPiZero  = (TH2D*)ESDContainer->FindObject(ObjectNameESDFixedPzPiZero.Data());


    fBckInvMassVSPt[0]  =(TH2D*) fBckInvMassVSPt[1]->Clone("hist_bck0");
    fBckInvMassVSPt[0]->Add(fBckInvMassVSPt[2]);
    fBckInvMassVSPt[0]->Add(fBckInvMassVSPt[3]);
    fBckInvMassVSPt[0]->Add(fBckInvMassVSPt[4]);

    fBckInvMassVSPt_SubPiZero[0]  =(TH2D*) fBckInvMassVSPt_SubPiZero[1]->Clone("hist_bck0_SubPiZero");
    fBckInvMassVSPt_SubPiZero[0]->Add(fBckInvMassVSPt_SubPiZero[2]);
    fBckInvMassVSPt_SubPiZero[0]->Add(fBckInvMassVSPt_SubPiZero[3]);
    fBckInvMassVSPt_SubPiZero[0]->Add(fBckInvMassVSPt_SubPiZero[4]);

    fBckInvMassVSPt_FixedPzPiZero[0]  =(TH2D*) fBckInvMassVSPt_FixedPzPiZero[1]->Clone("hist_bck0_FixedPzPiZero");
    fBckInvMassVSPt_FixedPzPiZero[0]->Add(fBckInvMassVSPt_FixedPzPiZero[2]);
    fBckInvMassVSPt_FixedPzPiZero[0]->Add(fBckInvMassVSPt_FixedPzPiZero[3]);
    fBckInvMassVSPt_FixedPzPiZero[0]->Add(fBckInvMassVSPt_FixedPzPiZero[4]);


    const char* FileDataLogname     = Form("%s/%s/%s_%s_EffiCheck_RAWDATA%s_%s.dat", cutSelection.Data(), fEnergyFlag.Data(), fPrefix.Data(), fPrefix2.Data(), fPeriodFlag.Data(), fCutSelectionRead.Data());
    fFileDataLog.open(FileDataLogname, ios::out);

    ProduceBckWithoutWeighting(fBckInvMassVSPt,fBckInvMassVSPt_SubPiZero,fBckInvMassVSPt_FixedPzPiZero); //background without weighting because THSparse wasn't used


    if(fIsMC){
        cout<<"STARTED BLOCK fIsMC"<<endl;
        TList *MCContainer                  = (TList*)HistosGammaConversion->FindObject(Form("%s MC histograms",fCutSelectionRead.Data()));
        TList *TrueConversionContainer      = (TList*)HistosGammaConversion->FindObject(Form("%s True histograms",fCutSelectionRead.Data()));

        if(MCContainer == NULL){
            cout<<"ERROR: " << Form("%s MC histograms",fCutSelectionRead.Data()) << " not Found in File"<<endl;
            return;
        } else
            cout<<"MC analysis: MCContainer successfully initialized..."<<endl;


        if(TrueConversionContainer == NULL){
            cout<<"ERROR: " << Form("%s True histograms",fCutSelectionRead.Data()) << " not Found in File"<<endl;
            return;
        } else
            cout<<"MC analysis: TrueConversionContainer successfully initialized..."<<endl;

        if( fMesonId == 221){
            fHistoMCMesonPtWithinAcceptance = (TH1D*)MCContainer->FindObject(ObjectNameMCEtaAcc.Data());
            fHistoMCMesonPt                 = (TH1D*)MCContainer->FindObject(ObjectNameMCEta.Data());   // Not the best; better having a 2D Pt_vs_Rapid in case we change limits
            fHistoMCMesonPtWOWeights =(TH1D*)MCContainer->FindObject(ObjectNameMCEtaWOWeights.Data());
        } else if( fMesonId == 223){
            fHistoMCMesonPtWithinAcceptance = (TH1D*)MCContainer->FindObject(ObjectNameMCOmegaAcc.Data());
            fHistoMCMesonPt                 = (TH1D*)MCContainer->FindObject(ObjectNameMCOmega.Data());   // Not the best; better having a 2D Pt_vs_Rapid in case we change limits
            fHistoMCMesonPtWOWeights        = (TH1D*)MCContainer->FindObject(ObjectNameMCOmegaWOWeights.Data());
        }
        fHistoMCMesonPt->Sumw2();
        fHistoMCMesonPtWithinAcceptance->Sumw2();
        if (fHistoMCMesonPtWOWeights){
            fHistoMCMesonPtWeights = (TH1D*)fHistoMCMesonPtWOWeights->Clone("WeightsMeson");
            fHistoMCMesonPtWeights->Divide(fHistoMCMesonPt,fHistoMCMesonPtWOWeights, 1.,1.,"B");
        }

        fHistoTrueMesonInvMassVSPt          = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrue.Data());
        fHistoTrueMesonInvMassVSPtWOWeights = (TH2D*)TrueConversionContainer->FindObject(ObjectNameTrueWOWeights.Data());
        if (fHistoTrueMesonInvMassVSPtWOWeights == NULL) fHistoTrueMesonInvMassVSPtWOWeights = (TH2D*)fHistoTrueMesonInvMassVSPt->Clone("WOWeights");

        fProfileTrueMesonInvMassVSPtWeights = (TProfile2D*)TrueConversionContainer->FindObject(ObjectNameProfileWeights.Data());
        if (fProfileTrueMesonInvMassVSPtWeights == NULL){
            fProfileTrueMesonInvMassVSPtWeights = new TProfile2D("name_TP2D","name_TP2D", fHistoTrueMesonInvMassVSPt->GetNbinsX(), 0., 2., fHistoTrueMesonInvMassVSPt->GetNbinsY(), 0., 25., "");
            for (Int_t i = 1; i<fProfileTrueMesonInvMassVSPtWeights->GetNbinsX(); i++)
                for (Int_t j = 1; j<fProfileTrueMesonInvMassVSPtWeights->GetNbinsY(); j++)
                    fProfileTrueMesonInvMassVSPtWeights->SetBinContent(i,j,1);
        }

        fHistoTrueMesonInvMassVSPtReweighted = (TH2D*)fHistoTrueMesonInvMassVSPtWOWeights->Clone("Reweighted");
        fHistoTrueMesonInvMassVSPtReweighted->Multiply(fProfileTrueMesonInvMassVSPtWeights);
        FillMassMCTrueMesonHistosArray(fHistoTrueMesonInvMassVSPt);
        FillMassMCTrueReweightedMesonHistosArray(fHistoTrueMesonInvMassVSPtReweighted);
    }

    fMesonMassExpect = TDatabasePDG::Instance()->GetParticle(fMesonId)->Mass();
    if (fEnergyFlag.Contains("PbPb") || fEnergyFlag.Contains("pPb")){
        fNEvents = fEventQuality->GetBinContent(1);
    } else {
        fNEvents = GetNEvents(fEventQuality);
    }

    TH1D *fBck[5] = {NULL,NULL,NULL,NULL,NULL};
    TH1D *fBck_SubPiZero[5] = {NULL,NULL,NULL,NULL,NULL};
    TH1D *fBck_FixedPzPiZero[5] = {NULL,NULL,NULL,NULL,NULL};
    for(Int_t k = 0; k< 5; k++){
        fBck[k]                     =  (TH1D*)fBckInvMassVSPt[k]->ProjectionX("");
        fBck[k]->Rebin(fNRebinGlobal);
        fBck_SubPiZero[k]           =  (TH1D*)fBckInvMassVSPt[k]->ProjectionX("");
        fBck_SubPiZero[k]->Rebin(fNRebinGlobal);
        fBck_FixedPzPiZero[k]       =  (TH1D*)fBckInvMassVSPt[k]->ProjectionX("");
        fBck_FixedPzPiZero[k]->Rebin(fNRebinGlobal);
    }
    TH1D *fGammaGamma                 = (TH1D*)fGammaGammaInvMassVSPt->ProjectionX("ESD_Mother_InvMass");
    TH1D *fGammaGamma_SubPiZero       = (TH1D*)fGammaGammaInvMassVSPt->ProjectionX("ESD_Mother_InvMass_SubPiZero");
    TH1D *fGammaGamma_FixedPzPiZero   = (TH1D*)fGammaGammaInvMassVSPt->ProjectionX("ESD_Mother_InvMass_FixedPzPiZero");

    cout<< "The mass of the meson is: "<< fMesonMassExpect<< " Events analysed: "<< fNEvents<< endl;
    cout << "here" << endl;
    // Process the 1D invariant mass histos
    fGammaGamma->SetTitle(Form("%s %s",fGammaGamma->GetTitle(),fCutSelection.Data()));
    fGammaGamma_SubPiZero->SetTitle(Form("%s %s",fGammaGamma->GetTitle(),fCutSelection.Data()));
    fGammaGamma_FixedPzPiZero->SetTitle(Form("%s %s",fGammaGamma->GetTitle(),fCutSelection.Data()));
    fGammaGamma->Rebin(fNRebinGlobal);
    fGammaGamma_SubPiZero->Rebin(fNRebinGlobal);
    fGammaGamma_FixedPzPiZero->Rebin(fNRebinGlobal);
    //fBck->Scale(1./fNRebinGlobal);
    ProcessEM( fGammaGamma , fBck[0], fBGFitRange); // scale the added background
    fHistoMappingBackNormInvMass = fBckNorm;
    fHistoMappingSignalInvMass = fSignal;

    ProcessEM( fGammaGamma_SubPiZero , fBck_SubPiZero[0], fBGFitRange_SubPiZero);
    fHistoMappingBackNormInvMass_SubPiZero = fBckNorm;
    fHistoMappingSignalInvMass_SubPiZero = fSignal;

    ProcessEM( fGammaGamma_FixedPzPiZero , fBck_FixedPzPiZero[0], fBGFitRange_FixedPzPiZero);
    fHistoMappingBackNormInvMass_FixedPzPiZero = fBckNorm;
    fHistoMappingSignalInvMass_FixedPzPiZero = fSignal;

    fGammaGamma->DrawCopy();
    fHistoMappingBackNormInvMass->DrawCopy("same");
    fHistoMappingSignalInvMass->DrawCopy("same");

    fGammaGamma_SubPiZero->DrawCopy();
    fHistoMappingBackNormInvMass_SubPiZero->DrawCopy("same");
    fHistoMappingSignalInvMass_SubPiZero->DrawCopy("same");

    fGammaGamma_FixedPzPiZero->DrawCopy();
    fHistoMappingBackNormInvMass_FixedPzPiZero->DrawCopy("same");
    fHistoMappingSignalInvMass_FixedPzPiZero->DrawCopy("same");

    // Function to Project the 2D histos InvariantMass VS Pt into Invariant Mass spectrum
    // Only do it for added background for now
    // TODO: Check if I have done this correctly
    FillMassHistosArray(fGammaGammaInvMassVSPt,fGammaGammaInvMassVSPt_SubPiZero,fGammaGammaInvMassVSPt_FixedPzPiZero);
    ProcessEM( fMesonFullPtSignal, fMesonFullPtBackground[0], fBGFitRange);
    fMesonFullPtBackNorm            = fBckNorm;
    ProcessEM( fFittingHistMidPtSignal, fFittingHistMidPtBackground[0], fBGFitRange);
    fFittingHistMidPtSignalSub = fSignal;

    ProcessEM( fMesonFullPtSignal_SubPiZero, fMesonFullPtBackground_SubPiZero[0], fBGFitRange_SubPiZero);
    fMesonFullPtBackNorm_SubPiZero            = fBckNorm;
    ProcessEM( fFittingHistMidPtSignal_SubPiZero, fFittingHistMidPtBackground_SubPiZero[0], fBGFitRange_SubPiZero);
    fFittingHistMidPtSignalSub_SubPiZero = fSignal;

    ProcessEM( fMesonFullPtSignal_FixedPzPiZero, fMesonFullPtBackground_FixedPzPiZero[0], fBGFitRange_FixedPzPiZero);
    fMesonFullPtBackNorm_FixedPzPiZero            = fBckNorm;
    ProcessEM( fFittingHistMidPtSignal_FixedPzPiZero, fFittingHistMidPtBackground_FixedPzPiZero[0], fBGFitRange_FixedPzPiZero);
    fFittingHistMidPtSignalSub_FixedPzPiZero = fSignal;
    if(fCrysFitting==0){
        fFileErrLog << "Using exp fit"<<endl;
        FitSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub, fMesonIntDeltaRange,200,kTRUE,0);
        fFitSignalInvMassMidPt = fFitReco;

        FitSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub_SubPiZero, fMesonIntDeltaRange,200,kTRUE,1);
        fFitSignalInvMassMidPt_SubPiZero = fFitReco;

        FitSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub_FixedPzPiZero, fMesonIntDeltaRange,200,kTRUE,2);
        fFitSignalInvMassMidPt_FixedPzPiZero = fFitReco;
    } else {
        fFileErrLog << "Using Crystal Ball function"<<endl;
        FitCBSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub, fMesonIntDeltaRange,200,kTRUE,"SinglefitfunctionMidPt",0);
        fFitSignalInvMassMidPt = fFitReco;

        FitCBSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub_SubPiZero, fMesonIntDeltaRange,200,kTRUE,"SinglefitfunctionMidPt",1);
        fFitSignalInvMassMidPt_SubPiZero = fFitReco;

        FitCBSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub_FixedPzPiZero, fMesonIntDeltaRange,200,kTRUE,"SinglefitfunctionMidPt",2);
        fFitSignalInvMassMidPt_FixedPzPiZero = fFitReco;
    }

    if (fIsMC){
        TH1D* fHistoMappingTrueMesonInvMassPtMidPt= NULL;
        fHistoMappingTrueMesonInvMassPtMidPt= new TH1D("TrueMassMidPt","TrueMassMidPt",fHistoTrueMesonInvMassVSPtReweighted->GetNbinsX(),0.,fHistoTrueMesonInvMassVSPtReweighted->GetXaxis()->GetBinUpEdge(fHistoTrueMesonInvMassVSPtReweighted->GetNbinsX()));
        fHistoMappingTrueMesonInvMassPtMidPt->Sumw2();
        Int_t startBin = fHistoTrueMesonInvMassVSPtReweighted->GetYaxis()->FindBin(fMidPt[0]+0.001);
        Int_t endBin = fHistoTrueMesonInvMassVSPtReweighted->GetYaxis()->FindBin(fMidPt[1]-0.001);

        fHistoTrueMesonInvMassVSPtReweighted->ProjectionX("TrueMassMidPt",startBin,endBin,"e");

        fHistoMappingTrueMesonInvMassPtMidPt=(TH1D*)gDirectory->Get(fNameHistoTrue.Data());
        fHistoMappingTrueMesonInvMassPtMidPt->Rebin(fNRebin[5]);
        fHistoMappingTrueMesonInvMassPtMidPt->SetLineWidth(1);
        fHistoMappingTrueMesonInvMassPtMidPt->SetLineColor(2);
        FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtMidPt, fMesonIntDeltaRange,50,kTRUE);


    }

    TString fDecayChannel = "#gamma#gamma";
    if (meson.CompareTo("Omega")==0) fDecayChannel = "#pi^{+} #pi^{-} #pi^{0}";
    if (meson.CompareTo("Eta")==0) fDecayChannel = "#pi^{+} #pi^{-} #pi^{0}";
    delete fMidPt;

    TFile *filefull = TFile::Open("histbackfull.root", "recreate");
    //writes background  and signal for full pt in a root file - carolina's modification
    filefull->cd();
    fOmegaFullPtSignal               = (TH1D*)fMesonFullPtSignal->Clone("fMesonFullPtSignal");
    fOmegaFullPtSignal_SubPiZero     = (TH1D*)fMesonFullPtSignal_SubPiZero->Clone("fMesonFullPtSignal_SubPiZero");
    fOmegaFullPtSignal_FixedPzPiZero = (TH1D*)fMesonFullPtSignal_FixedPzPiZero->Clone("fMesonFullPtSignal_FixedPzPiZero");
    fOmegaFullPtSignal->Write("fSignal");
    fOmegaFullPtSignal_SubPiZero->Write("fSignal_SubPiZero");
    fOmegaFullPtSignal_FixedPzPiZero->Write("fSignal_FixedPzPiZero");
    // Different event mixing background groups
    for(Int_t k = 0; k<4;k++){
        fOmegaFullPtBack[k]                 = (TH1D*)hist_bck[k]->Clone(Form("hist_bck%i",k));
        fOmegaFullPtBack_SubPiZero[k]       = (TH1D*)hist_bck_SubPiZero[k]->Clone(Form("hist_bck%i_SubPiZero",k));
        fOmegaFullPtBack_FixedPzPiZero[k]   = (TH1D*)hist_bck_FixedPzPiZero[k]->Clone(Form("hist_bck%i_FixedPzPiZero",k));
        fOmegaFullPtBack[k]->Write(Form("fBack%i",k));
        fOmegaFullPtBack_SubPiZero[k]->Write(Form("fBack%i_SubPiZero",k));
        fOmegaFullPtBack_FixedPzPiZero[k]->Write(Form("fBack%i_FixedPzPiZero",k));
    }

    fOmegaFullPtBackNorm               = (TH1D*)fMesonFullPtBackNorm->Clone("fMesonFullPtBackNorm");
    fOmegaFullPtBackNorm_SubPiZero     = (TH1D*)fMesonFullPtBackNorm_SubPiZero->Clone("fMesonFullPtBackNorm_SubPiZero");
    fOmegaFullPtBackNorm_FixedPzPiZero = (TH1D*)fMesonFullPtBackNorm_FixedPzPiZero->Clone("fMesonFullPtBackNorm_FixedPzPiZero");
    fOmegaFullPtBackNorm->Write("fBackNorm");
    fOmegaFullPtBackNorm_SubPiZero->Write("fBackNorm_SubPiZero");
    fOmegaFullPtBackNorm_FixedPzPiZero->Write("fBackNorm_FixedPzPiZero");
    filefull->Close();


    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){ // BEGIN ANALYSIS for each Pt bin

        cout << "---------------------------------------------------------------------------------" << endl;
        cout << "Begin Analysis Pt Bin " << iPt <<endl;
        cout << "---------------------------------------------------------------------------------" << endl;
        // Function to subtract GG minus Bck
        fFileDataLog << "---------------------------------------------------------------------------------" << endl;
        fFileDataLog << "----------------------------------new pT bin ------------------------------------" << endl;
        fFileDataLog << "---------------------------------------------------------------------------------" << endl;

        ProcessBckFitSubtraction(fHistoMappingGGInvMassPtBin[iPt],iPt,fPeakRange,fFitRange,optionEnergy,Suffix,cutSelection,meson,0);
        ProcessBckFitSubtraction(fHistoMappingGGInvMassPtBin_SubPiZero[iPt],iPt,fPeakRange_SubPiZero,fFitRange_SubPiZero,optionEnergy,Suffix,cutSelection,meson,1);
        ProcessBckFitSubtraction(fHistoMappingGGInvMassPtBin_FixedPzPiZero[iPt],iPt,fPeakRange_FixedPzPiZero,fFitRange_FixedPzPiZero,optionEnergy,Suffix,cutSelection,meson,2);
        Double_t fNormTot = 1.; // variable containing the factor used to scale the total background
        for(Int_t k=0;k<5;k++){
          if((fHistoMappingGGInvMassBackFitPtBin[iPt]!=NULL) && (fHistoMappingBackInvMassPtBin[k][iPt]!=NULL )){
            ProcessEM( fHistoMappingGGInvMassPtBin[iPt], fHistoMappingBackInvMassPtBin[k][iPt], fBGFitRange);
            if(k==0) fNormTot = fNorm; // k==0 means tot back was just scaled -> store value
            fHistoMappingSignalInvMassPtBin[iPt] = fSignal;
            fHistoMappingBackNormInvMassPtBin[k][iPt] = fBckNorm;

            fHistoMappingBackSameNormInvMassPtBin[k][iPt] = (TH1D*) fHistoMappingBackInvMassPtBin[k][iPt]->Clone(Form("BckSameNorm%i",k)); // clone all the unscaled bckgroups
            fHistoMappingBackSameNormInvMassPtBin[k][iPt]->Sumw2();
            fHistoMappingBackSameNormInvMassPtBin[k][iPt]->Scale(fNormTot); // scale verything with the SAME factor -> factor obtained from scaling of total background
          }

          if((fHistoMappingGGInvMassBackFitPtBin_SubPiZero[iPt]!=NULL) && (fHistoMappingBackInvMassPtBin_SubPiZero[k][iPt]!=NULL )){
            ProcessEM( fHistoMappingGGInvMassPtBin_SubPiZero[iPt], fHistoMappingBackInvMassPtBin_SubPiZero[k][iPt], fBGFitRange_SubPiZero);
            if(k==0) fNormTot = fNorm; // k==0 means tot back was just scaled -> store value
            fHistoMappingSignalInvMassPtBin_SubPiZero[iPt] = fSignal;
            fHistoMappingBackNormInvMassPtBin_SubPiZero[k][iPt] = fBckNorm;

            fHistoMappingBackSameNormInvMassPtBin_SubPiZero[k][iPt] = (TH1D*) fHistoMappingBackInvMassPtBin_SubPiZero[k][iPt]->Clone(Form("BckSameNorm%i_SubPiZero",k)); // clone all the unscaled bckgroups
            fHistoMappingBackSameNormInvMassPtBin_SubPiZero[k][iPt]->Sumw2();
            fHistoMappingBackSameNormInvMassPtBin_SubPiZero[k][iPt]->Scale(fNormTot); // scale verything with the SAME factor -> factor obtained from scaling of total background
          }

          if((fHistoMappingGGInvMassBackFitPtBin_FixedPzPiZero[iPt]!=NULL) && (fHistoMappingBackInvMassPtBin_FixedPzPiZero[k][iPt]!=NULL )){
            ProcessEM( fHistoMappingGGInvMassPtBin_FixedPzPiZero[iPt], fHistoMappingBackInvMassPtBin_FixedPzPiZero[k][iPt], fBGFitRange_FixedPzPiZero);
            if(k==0) fNormTot = fNorm; // k==0 means tot back was just scaled -> store value
            fHistoMappingSignalInvMassPtBin_FixedPzPiZero[iPt] = fSignal;
            fHistoMappingBackNormInvMassPtBin_FixedPzPiZero[k][iPt] = fBckNorm;

            fHistoMappingBackSameNormInvMassPtBin_FixedPzPiZero[k][iPt] = (TH1D*) fHistoMappingBackInvMassPtBin_FixedPzPiZero[k][iPt]->Clone(Form("BckSameNorm%i_FixedPzPiZero",k)); // clone all the unscaled bckgroups
            fHistoMappingBackSameNormInvMassPtBin_FixedPzPiZero[k][iPt]->Sumw2();
            fHistoMappingBackSameNormInvMassPtBin_FixedPzPiZero[k][iPt]->Scale(fNormTot); // scale verything with the SAME factor -> factor obtained from scaling of total background
          }
        }

        TString namesecHistoBckNorm;
        TString namesecHistoSignal;
        TString namesecRatio;

        TString tempStr = (iPt==fStartPtBin)?"recreate":"update";

        TFile *file1 = TFile::Open("histback.root",tempStr.Data());
        //writes background  and signal for each pt bin in a root file - carolina's modification
        file1->cd();
        for(Int_t k=0;k<5;k++){
            namesecHistoBckNorm                         = Form("BckNorm%i_InvMass_in_Pt_Bin%02i",k, iPt);
            fHistoMappingBackNormInvMassPtBin[k][iPt]->Write(namesecHistoBckNorm.Data());
            namesecHistoBckNorm                         = Form("BckNorm%i_InvMass_SubPiZero_in_Pt_Bin%02i",k, iPt);
            if (fHistoMappingBackNormInvMassPtBin_SubPiZero[k][iPt]!=NULL) fHistoMappingBackNormInvMassPtBin_SubPiZero[k][iPt]->Write(namesecHistoBckNorm.Data());
            namesecHistoBckNorm                         = Form("BckNorm%i_InvMass_FixedPzPiZero_in_Pt_Bin%02i",k, iPt);
            if (fHistoMappingBackNormInvMassPtBin_FixedPzPiZero[k][iPt]!=NULL )fHistoMappingBackNormInvMassPtBin_FixedPzPiZero[k][iPt]->Write(namesecHistoBckNorm.Data());

        }
        namesecHistoSignal                          = Form("Signal_InvMass_in_Pt_Bin%02i", iPt);
        fHistoMappingGGInvMassPtBin[iPt]->Write(namesecHistoSignal.Data());
        namesecHistoSignal                          = Form("Signal_InvMass_SubPiZero_in_Pt_Bin%02i", iPt);
        fHistoMappingGGInvMassPtBin_SubPiZero[iPt]->Write(namesecHistoSignal.Data());
        namesecHistoSignal                          = Form("Signal_InvMass_FixedPzPiZero_in_Pt_Bin%02i", iPt);
        fHistoMappingGGInvMassPtBin_FixedPzPiZero[iPt]->Write(namesecHistoSignal.Data());

        file1->Close();

        // TODO how do you use this function without using fMesonFitRange
        FitWithPol2ForBG(fHistoMappingGGInvMassPtBin[iPt], fMesonFitRange,iPt,kFALSE,0);
        fFitWithPol2ForBG[iPt]                       = fFitReco;

        FitWithPol2ForBG(fHistoMappingGGInvMassPtBin_SubPiZero[iPt], fMesonFitRange_SubPiZero,iPt,kFALSE,1);
        fFitWithPol2ForBG_SubPiZero[iPt]                       = fFitReco;

        FitWithPol2ForBG(fHistoMappingGGInvMassPtBin_FixedPzPiZero[iPt], fMesonFitRange_FixedPzPiZero,iPt,kFALSE,2);
        fFitWithPol2ForBG_FixedPzPiZero[iPt]                       = fFitReco;

        ProcessRatioSignalBackground(fHistoMappingGGInvMassPtBin[iPt], fHistoMappingBackNormInvMassPtBin[0][iPt]);
        fHistoMappingRatioSBInvMassPtBin[0][iPt]        = fRatioSB; // only implemented for added background so far

        ProcessRatioSignalBackground(fHistoMappingGGInvMassPtBin_SubPiZero[iPt], fHistoMappingBackNormInvMassPtBin_SubPiZero[0][iPt]);
        fHistoMappingRatioSBInvMassPtBin_SubPiZero[0][iPt]        = fRatioSB; // only implemented for added background so far

        ProcessRatioSignalBackground(fHistoMappingGGInvMassPtBin_FixedPzPiZero[iPt], fHistoMappingBackNormInvMassPtBin_FixedPzPiZero[0][iPt]);
        fHistoMappingRatioSBInvMassPtBin_FixedPzPiZero[0][iPt]        = fRatioSB; // only implemented for added background so far

        fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "normal range/right normalization" << endl;

        fFitSignalInvMassPtBin[iPt]                       = 0x00;
        fFitSignalInvMassPtBin_SubPiZero[iPt]             = 0x00;
        fFitSignalInvMassPtBin_FixedPzPiZero[iPt]         = 0x00;
        fFitSignalInvMassBackFitPtBin[iPt]                = 0x00;
        fFitSignalInvMassBackFitPtBin_SubPiZero[iPt]      = 0x00;
        fFitSignalInvMassBackFitPtBin_FixedPzPiZero[iPt]  = 0x00;
        if(fCrysFitting==0){
            fFileErrLog << "Using exp fit"<<endl;
            fFileDataLog << "Subtracted mixed event" << endl;
            if(fHistoMappingSignalInvMassPtBin[iPt]!=NULL){
              FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE,0);
              fFitSignalInvMassPtBin[iPt]                     = fFitReco;
              fFitBckInvMassPtBin[iPt]                        = fFitLinearBck;
              fMesonYieldsResidualBckFunc[0][iPt]             = fIntLinearBck;
              fMesonYieldsResidualBckFuncError[0][iPt]        = fIntLinearBckError;
              fFileDataLog << "Subtracted polinomial fit" << endl;
              FitSubtractedInvMassInPtBins(fHistoMappingGGInvMassBackFitPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE,0);
              fFitSignalInvMassBackFitPtBin[iPt]              = fFitReco;
              fFitBckInvMassBackFitPtBin[iPt]                 = fFitLinearBck;
              fMesonYieldsResidualBckFuncBackFit[iPt]         = fIntLinearBck;
              fMesonYieldsResidualBckFuncBackFitError[iPt]    = fIntLinearBckError;
            }
            // Calculation for pi0 mass subtracted
            if(fHistoMappingSignalInvMassPtBin_SubPiZero[iPt]!=NULL){
              FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin_SubPiZero[iPt], fMesonIntDeltaRange,iPt,kFALSE,1);
              fFitSignalInvMassPtBin_SubPiZero[iPt]                     = fFitReco;
              fFitBckInvMassPtBin_SubPiZero[iPt]                        = fFitLinearBck;
              fMesonYieldsResidualBckFunc_SubPiZero[0][iPt]             = fIntLinearBck;
              fMesonYieldsResidualBckFuncError_SubPiZero[0][iPt]        = fIntLinearBckError;
              fFileDataLog << "Subtracted polinomial fit" << endl;
              FitSubtractedInvMassInPtBins(fHistoMappingGGInvMassBackFitPtBin_SubPiZero[iPt], fMesonIntDeltaRange,iPt,kFALSE,1);
              fFitSignalInvMassBackFitPtBin_SubPiZero[iPt]              = fFitReco;
              fFitBckInvMassBackFitPtBin_SubPiZero[iPt]                 = fFitLinearBck;
              fMesonYieldsResidualBckFuncBackFit_SubPiZero[iPt]         = fIntLinearBck;
              fMesonYieldsResidualBckFuncBackFitError_SubPiZero[iPt]    = fIntLinearBckError;
            }

            // Calculation for pz of Pi0 fixed
            if(fHistoMappingSignalInvMassPtBin_FixedPzPiZero[iPt]!=NULL){
              FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin_FixedPzPiZero[iPt], fMesonIntDeltaRange,iPt,kFALSE,2);
              fFitSignalInvMassPtBin_FixedPzPiZero[iPt]                     = fFitReco;
              fFitBckInvMassPtBin_FixedPzPiZero[iPt]                        = fFitLinearBck;
              fMesonYieldsResidualBckFunc_FixedPzPiZero[0][iPt]             = fIntLinearBck;
              fMesonYieldsResidualBckFuncError_FixedPzPiZero[0][iPt]        = fIntLinearBckError;
              fFileDataLog << "Subtracted polinomial fit" << endl;
              FitSubtractedInvMassInPtBins(fHistoMappingGGInvMassBackFitPtBin_FixedPzPiZero[iPt], fMesonIntDeltaRange,iPt,kFALSE,2);
              fFitSignalInvMassBackFitPtBin_FixedPzPiZero[iPt]              = fFitReco;
              fFitBckInvMassBackFitPtBin_FixedPzPiZero[iPt]                 = fFitLinearBck;
              fMesonYieldsResidualBckFuncBackFit_FixedPzPiZero[iPt]         = fIntLinearBck;
              fMesonYieldsResidualBckFuncBackFitError_FixedPzPiZero[iPt]    = fIntLinearBckError;
            }

        } else {
          if (fHistoMappingSignalInvMassPtBin[iPt] != NULL){
            fFileErrLog << "Using Crystal Ball function"<<endl;
            FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncNormalBin%02d",iPt),0);
            fFitSignalInvMassPtBin[iPt]                     = fFitReco;
            fFitBckInvMassPtBin[iPt]                        = fFitLinearBck;
            fMesonYieldsResidualBckFunc[0][iPt]             = fIntLinearBck;
            fMesonYieldsResidualBckFuncError[0][iPt]        = fIntLinearBckError;
          }

          // Calculation for Inv Mass Pi0 subtracted
          if (fHistoMappingSignalInvMassPtBin_SubPiZero[iPt] != NULL){
            fFileErrLog << "Using Crystal Ball function"<<endl;
            FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin_SubPiZero[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncNormal_SubPiZero_Bin%02d",iPt),1);
            fFitSignalInvMassPtBin_SubPiZero[iPt]                     = fFitReco;
            fFitBckInvMassPtBin_SubPiZero[iPt]                        = fFitLinearBck;
            fMesonYieldsResidualBckFunc_SubPiZero[0][iPt]             = fIntLinearBck;
            fMesonYieldsResidualBckFuncError_SubPiZero[0][iPt]        = fIntLinearBckError;
          }

          //Calculation for Pz of Pi0 fixe
          if (fHistoMappingSignalInvMassPtBin_FixedPzPiZero[iPt] != NULL){
            fFileErrLog << "Using Crystal Ball function"<<endl;
            FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin_FixedPzPiZero[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncNormal_FixedPzPiZero_Bin%02d",iPt),2);
            fFitSignalInvMassPtBin_FixedPzPiZero[iPt]                     = fFitReco;
            fFitBckInvMassPtBin_FixedPzPiZero[iPt]                        = fFitLinearBck;
            fMesonYieldsResidualBckFunc_FixedPzPiZero[0][iPt]             = fIntLinearBck;
            fMesonYieldsResidualBckFuncError_FixedPzPiZero[0][iPt]        = fIntLinearBckError;
          }
        }

        //GetFWHM
        CalculateFWHM( fFitSignalInvMassPtBin[iPt],fMesonFitRange);
        fMesonFWHM[iPt]                                     = fFWHMFunc;
        fMesonFWHMError[iPt]                                = fFWHMFuncError;

        // Calculation for Pi0 InvMass subtracted
        CalculateFWHM( fFitSignalInvMassPtBin_SubPiZero[iPt],fMesonFitRange_SubPiZero);
        fMesonFWHM_SubPiZero[iPt]                           = fFWHMFunc;
        fMesonFWHMError_SubPiZero[iPt]                      = fFWHMFuncError;

        // Calculation for fixed Pi0 Pz
        CalculateFWHM( fFitSignalInvMassPtBin_FixedPzPiZero[iPt],fMesonFitRange_FixedPzPiZero);
        fMesonFWHM_FixedPzPiZero[iPt]                       = fFWHMFunc;
        fMesonFWHMError_FixedPzPiZero[iPt]                  = fFWHMFuncError;

        if (fFitSignalInvMassPtBin[iPt] !=0x00){
            fMesonMass[iPt]                 = fFitSignalInvMassPtBin[iPt]->GetParameter(1);
            fMesonMassError[iPt]            = fFitSignalInvMassPtBin[iPt]->GetParError(1);
            fMesonWidth[iPt]                = fFitSignalInvMassPtBin[iPt]->GetParameter(2);
            fMesonWidthError[iPt]           = fFitSignalInvMassPtBin[iPt]->GetParError(2);

            fMesonCurIntRange[0][0]         = fMesonMass[iPt] + fMesonIntDeltaRange[0];
            fMesonCurIntRange[1][0]         = fMesonMass[iPt] + fMesonIntDeltaRangeWide[0];
            fMesonCurIntRange[2][0]         = fMesonMass[iPt] + fMesonIntDeltaRangeNarrow[0];
            fMesonCurIntRange[0][1]         = fMesonMass[iPt] + fMesonIntDeltaRange[1];
            fMesonCurIntRange[1][1]         = fMesonMass[iPt] + fMesonIntDeltaRangeWide[1];
            fMesonCurIntRange[2][1]         = fMesonMass[iPt] + fMesonIntDeltaRangeNarrow[1];

            if (fFitSignalInvMassPtBin_SubPiZero !=0x00){
              fMesonMass_SubPiZero[iPt]                 = fFitSignalInvMassPtBin_SubPiZero[iPt]->GetParameter(1);
              fMesonMassError_SubPiZero[iPt]            = fFitSignalInvMassPtBin_SubPiZero[iPt]->GetParError(1);
              fMesonWidth_SubPiZero[iPt]                = fFitSignalInvMassPtBin_SubPiZero[iPt]->GetParameter(2);
              fMesonWidthError_SubPiZero[iPt]           = fFitSignalInvMassPtBin_SubPiZero[iPt]->GetParError(2);

              fMesonCurIntRange_SubPiZero[0][0]         = fMesonMass_SubPiZero[iPt] + fMesonIntDeltaRange[0];
              fMesonCurIntRange_SubPiZero[1][0]         = fMesonMass_SubPiZero[iPt] + fMesonIntDeltaRangeWide[0];
              fMesonCurIntRange_SubPiZero[2][0]         = fMesonMass_SubPiZero[iPt] + fMesonIntDeltaRangeNarrow[0];
              fMesonCurIntRange_SubPiZero[0][1]         = fMesonMass_SubPiZero[iPt] + fMesonIntDeltaRange[1];
              fMesonCurIntRange_SubPiZero[1][1]         = fMesonMass_SubPiZero[iPt] + fMesonIntDeltaRangeWide[1];
              fMesonCurIntRange_SubPiZero[2][1]         = fMesonMass_SubPiZero[iPt] + fMesonIntDeltaRangeNarrow[1];
            }

            if (fFitSignalInvMassPtBin_FixedPzPiZero !=0x00){
              fMesonMass_FixedPzPiZero[iPt]                 = fFitSignalInvMassPtBin_FixedPzPiZero[iPt]->GetParameter(1);
              fMesonMassError_FixedPzPiZero[iPt]            = fFitSignalInvMassPtBin_FixedPzPiZero[iPt]->GetParError(1);
              fMesonWidth_FixedPzPiZero[iPt]                = fFitSignalInvMassPtBin_FixedPzPiZero[iPt]->GetParameter(2);
              fMesonWidthError_FixedPzPiZero[iPt]           = fFitSignalInvMassPtBin_FixedPzPiZero[iPt]->GetParError(2);

              fMesonCurIntRange_FixedPzPiZero[0][0]         = fMesonMass_FixedPzPiZero[iPt] + fMesonIntDeltaRange[0];
              fMesonCurIntRange_FixedPzPiZero[1][0]         = fMesonMass_FixedPzPiZero[iPt] + fMesonIntDeltaRangeWide[0];
              fMesonCurIntRange_FixedPzPiZero[2][0]         = fMesonMass_FixedPzPiZero[iPt] + fMesonIntDeltaRangeNarrow[0];
              fMesonCurIntRange_FixedPzPiZero[0][1]         = fMesonMass_FixedPzPiZero[iPt] + fMesonIntDeltaRange[1];
              fMesonCurIntRange_FixedPzPiZero[1][1]         = fMesonMass_FixedPzPiZero[iPt] + fMesonIntDeltaRangeWide[1];
              fMesonCurIntRange_FixedPzPiZero[2][1]         = fMesonMass_FixedPzPiZero[iPt] + fMesonIntDeltaRangeNarrow[1];
            }
        } else {
            fMesonMass[iPt]                 = fMesonMassExpect;
            fMesonMassError[iPt]            = 0.;
            fMesonCurIntRange[0][0]         = fMesonMassExpect + fMesonIntDeltaRange[0];
            fMesonCurIntRange[1][0]         = fMesonMassExpect + fMesonIntDeltaRangeWide[0];
            fMesonCurIntRange[2][0]         = fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
            fMesonCurIntRange[0][1]         = fMesonMassExpect + fMesonIntDeltaRange[1];
            fMesonCurIntRange[1][1]         = fMesonMassExpect + fMesonIntDeltaRangeWide[1];
            fMesonCurIntRange[2][1]         = fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];

            fMesonMass_SubPiZero[iPt]                 = fMesonMassExpect-0.134;
            fMesonMassError_SubPiZero[iPt]            = 0.;
            fMesonCurIntRange_SubPiZero[0][0]         = fMesonMass_SubPiZero[iPt]   + fMesonIntDeltaRange[0];
            fMesonCurIntRange_SubPiZero[1][0]         = fMesonMass_SubPiZero[iPt]   + fMesonIntDeltaRangeWide[0];
            fMesonCurIntRange_SubPiZero[2][0]         = fMesonMass_SubPiZero[iPt]   + fMesonIntDeltaRangeNarrow[0];
            fMesonCurIntRange_SubPiZero[0][1]         = fMesonMass_SubPiZero[iPt]   + fMesonIntDeltaRange[1];
            fMesonCurIntRange_SubPiZero[1][1]         = fMesonMass_SubPiZero[iPt]   + fMesonIntDeltaRangeWide[1];
            fMesonCurIntRange_SubPiZero[2][1]         = fMesonMass_SubPiZero[iPt]   + fMesonIntDeltaRangeNarrow[1];

            fMesonMass_FixedPzPiZero[iPt]                 = fMesonMassExpect;
            fMesonMassError_FixedPzPiZero[iPt]            = 0.;
            fMesonCurIntRange_FixedPzPiZero[0][0]         = fMesonMass_FixedPzPiZero[iPt]   + fMesonIntDeltaRange[0];
            fMesonCurIntRange_FixedPzPiZero[1][0]         = fMesonMass_FixedPzPiZero[iPt]   + fMesonIntDeltaRangeWide[0];
            fMesonCurIntRange_FixedPzPiZero[2][0]         = fMesonMass_FixedPzPiZero[iPt]   + fMesonIntDeltaRangeNarrow[0];
            fMesonCurIntRange_FixedPzPiZero[0][1]         = fMesonMass_FixedPzPiZero[iPt]   + fMesonIntDeltaRange[1];
            fMesonCurIntRange_FixedPzPiZero[1][1]         = fMesonMass_FixedPzPiZero[iPt]   + fMesonIntDeltaRangeWide[1];
            fMesonCurIntRange_FixedPzPiZero[2][1]         = fMesonMass_FixedPzPiZero[iPt]   + fMesonIntDeltaRangeNarrow[1];

        }

        if (fFitSignalInvMassBackFitPtBin[iPt] !=0x00){ //TODO  look into the whole BackFit part
            cout<<fFitSignalInvMassBackFitPtBin[iPt]->GetParameter(1)<<endl;
            fMesonMassBackFit[iPt]          = fFitSignalInvMassBackFitPtBin[iPt]->GetParameter(1);
            fMesonMassBackFitError[iPt]     = fFitSignalInvMassBackFitPtBin[iPt]->GetParError(1);
            fMesonWidthBackFit[iPt]         = fFitSignalInvMassBackFitPtBin[iPt]->GetParameter(2);
            fMesonWidthBackFitError[iPt]    = fFitSignalInvMassBackFitPtBin[iPt]->GetParError(2);
            //fMesonCurIntRangeBackFit[0]     = fMesonIntRange[0] - (fMesonMassExpect-fMesonMassBackFit[iPt]);
            //fMesonCurIntRangeBackFit[1]     = fMesonIntRange[1] - (fMesonMassExpect-fMesonMassBackFit[iPt]);
            fMesonCurIntRangeBackFit[0]     = fMesonMass[iPt] + fMesonIntDeltaRange[0]; // TODO look into this (for now same as fMesonCurIntRange)
            fMesonCurIntRangeBackFit[1]     = fMesonMass[iPt] + fMesonIntDeltaRange[1];

            if(fFitSignalInvMassBackFitPtBin_SubPiZero !=0x00){
              fMesonMassBackFit_SubPiZero[iPt]          = fFitSignalInvMassBackFitPtBin_SubPiZero[iPt]->GetParameter(1);
              fMesonMassBackFitError_SubPiZero[iPt]     = fFitSignalInvMassBackFitPtBin_SubPiZero[iPt]->GetParError(1);
              fMesonWidthBackFit_SubPiZero[iPt]         = fFitSignalInvMassBackFitPtBin_SubPiZero[iPt]->GetParameter(2);
              fMesonWidthBackFitError_SubPiZero[iPt]    = fFitSignalInvMassBackFitPtBin_SubPiZero[iPt]->GetParError(2);
              //fMesonCurIntRangeBackFit[0]     = fMesonIntRange[0] - (fMesonMassExpect-fMesonMassBackFit[iPt]);
              //fMesonCurIntRangeBackFit[1]     = fMesonIntRange[1] - (fMesonMassExpect-fMesonMassBackFit[iPt]);
              fMesonCurIntRangeBackFit_SubPiZero[0]     = fMesonMass_SubPiZero[iPt] + fMesonIntDeltaRange[0]; // TODO look into this (for now same as fMesonCurIntRange)
              fMesonCurIntRangeBackFit_SubPiZero[1]     = fMesonMass_SubPiZero[iPt] + fMesonIntDeltaRange[1];
            }

            if(fFitSignalInvMassBackFitPtBin_FixedPzPiZero !=0x00){
              fMesonMassBackFit_FixedPzPiZero[iPt]          = fFitSignalInvMassBackFitPtBin_FixedPzPiZero[iPt]->GetParameter(1);
              fMesonMassBackFitError_FixedPzPiZero[iPt]     = fFitSignalInvMassBackFitPtBin_FixedPzPiZero[iPt]->GetParError(1);
              fMesonWidthBackFit_FixedPzPiZero[iPt]         = fFitSignalInvMassBackFitPtBin_FixedPzPiZero[iPt]->GetParameter(2);
              fMesonWidthBackFitError_FixedPzPiZero[iPt]    = fFitSignalInvMassBackFitPtBin_FixedPzPiZero[iPt]->GetParError(2);
              //fMesonCurIntRangeBackFit[0]     = fMesonIntRange[0] - (fMesonMassExpect-fMesonMassBackFit[iPt]);
              //fMesonCurIntRangeBackFit[1]     = fMesonIntRange[1] - (fMesonMassExpect-fMesonMassBackFit[iPt]);
              fMesonCurIntRangeBackFit_FixedPzPiZero[0]     = fMesonMass_FixedPzPiZero[iPt] + fMesonIntDeltaRange[0]; // TODO look into this (for now same as fMesonCurIntRange)
              fMesonCurIntRangeBackFit_FixedPzPiZero[1]     = fMesonMass_FixedPzPiZero[iPt] + fMesonIntDeltaRange[1];
            }
        } else {
            fMesonMassBackFit[iPt]          = 0;
            fMesonMassBackFitError[iPt]     = 0;
            fMesonWidthBackFit[iPt]         = 0;
            fMesonWidthBackFitError[iPt]    = 0;
            //MesonCurIntRangeBackFit[0]     = fMesonIntRange[0];
            //MesonCurIntRangeBackFit[1]     = fMesonIntRange[1];
            fMesonCurIntRangeBackFit[0]     = fMesonMass[iPt] + fMesonIntDeltaRange[0]; // TODO look into this (for now same as fMesonCurIntRange)
            fMesonCurIntRangeBackFit[1]     = fMesonMass[iPt] + fMesonIntDeltaRange[1];

            fMesonMassBackFit_SubPiZero[iPt]          = 0;
            fMesonMassBackFitError_SubPiZero[iPt]     = 0;
            fMesonWidthBackFit_SubPiZero[iPt]         = 0;
            fMesonWidthBackFitError_SubPiZero[iPt]    = 0;
            //MesonCurIntRangeBackFit[0]     = fMesonIntRange[0];
            //MesonCurIntRangeBackFit[1]     = fMesonIntRange[1];
            fMesonCurIntRangeBackFit_SubPiZero[0]     = fMesonMass[iPt] + fMesonIntDeltaRange[0]; // TODO look into this (for now same as fMesonCurIntRange)
            fMesonCurIntRangeBackFit_SubPiZero[1]     = fMesonMass[iPt] + fMesonIntDeltaRange[1];

            fMesonMassBackFit_FixedPzPiZero[iPt]          = 0;
            fMesonMassBackFitError_FixedPzPiZero[iPt]     = 0;
            fMesonWidthBackFit_FixedPzPiZero[iPt]         = 0;
            fMesonWidthBackFitError_FixedPzPiZero[iPt]    = 0;
            //MesonCurIntRangeBackFit[0]     = fMesonIntRange[0];
            //MesonCurIntRangeBackFit[1]     = fMesonIntRange[1];
            fMesonCurIntRangeBackFit_FixedPzPiZero[0]     = fMesonMass[iPt] + fMesonIntDeltaRange[0]; // TODO look into this (for now same as fMesonCurIntRange)
            fMesonCurIntRangeBackFit_FixedPzPiZero[1]     = fMesonMass[iPt] + fMesonIntDeltaRange[1];

        }

        for (Int_t k = 0; k < 3; k++){
            IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRange[k]);
            fGGYields[k][iPt]               = fYields;
            fGGYieldsError[k][iPt]          = fYieldsError;

            // Integrate histo with InvMass pi0 subtracted
            if( fHistoMappingGGInvMassPtBin_SubPiZero[iPt] != NULL){
              IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin_SubPiZero[iPt], fMesonCurIntRange_SubPiZero[k]);
              fGGYields_SubPiZero[k][iPt]               = fYields;
              fGGYieldsError_SubPiZero[k][iPt]          = fYieldsError;
            }

            // Integrate histo with pz of pi0 fixed
            if(fHistoMappingGGInvMassPtBin_FixedPzPiZero[iPt] != NULL){
              IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin_FixedPzPiZero[iPt], fMesonCurIntRange_FixedPzPiZero[k]);
              fGGYields_FixedPzPiZero[k][iPt]               = fYields;
              fGGYieldsError_FixedPzPiZero[k][iPt]          = fYieldsError;
            }

            // Integrate the bck histo
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[0][iPt], fMesonCurIntRange[k]);
            fBckYields[k][iPt]              = fYields;
            fBckYieldsError[k][iPt]         = fYieldsError;

            // Integrate the bck histo with InvMass pi0 subtracted
            if(  fHistoMappingBackNormInvMassPtBin_SubPiZero[0][iPt] != NULL){
              IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin_SubPiZero[0][iPt], fMesonCurIntRange_SubPiZero[k]);
              fBckYields_SubPiZero[k][iPt]              = fYields;
              fBckYieldsError_SubPiZero[k][iPt]         = fYieldsError;
            }

            // Integrate the bck histo with Pz of pi0 fixed
            if(fHistoMappingBackNormInvMassPtBin_FixedPzPiZero[0][iPt] != NULL){
              IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin_FixedPzPiZero[0][iPt], fMesonCurIntRange_FixedPzPiZero[k]);
              fBckYields_FixedPzPiZero[k][iPt]              = fYields;
              fBckYieldsError_FixedPzPiZero[k][iPt]         = fYieldsError;
            }

            if (k == 0){
                // Integrate the Backfit histo
                fFileDataLog<< endl <<"BackFit histo "<< nameIntRange[k].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
                IntegrateHistoInvMassStream( fHistoMappingGGInvMassBackFitPtBin[iPt], fMesonCurIntRangeBackFit); // fMesonCurIntRangeBackFit no array ?
                fMesonYieldsBackFit[iPt]      = fYields; // change this part because it is written each time of the loop -> TODO Change to array or move outside of loop
                fMesonYieldsBackFitError[iPt] = fYieldsError;

                // Integrate the Backfit histo with InvMass pi0 subtracted
                fFileDataLog<< endl <<"BackFit histo "<< nameIntRange[k].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
                IntegrateHistoInvMassStream( fHistoMappingGGInvMassBackFitPtBin_SubPiZero[iPt], fMesonCurIntRangeBackFit_SubPiZero); // fMesonCurIntRangeBackFit no array ?
                fMesonYieldsBackFit_SubPiZero[iPt]      = fYields; // change this part because it is written each time of the loop -> TODO Change to array or move outside of loop
                fMesonYieldsBackFitError_SubPiZero[iPt] = fYieldsError;

                // Integrate the Backfit histo with Pz of pi0 fixed
                fFileDataLog<< endl <<"BackFit histo "<< nameIntRange[k].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
                IntegrateHistoInvMassStream( fHistoMappingGGInvMassBackFitPtBin_FixedPzPiZero[iPt], fMesonCurIntRangeBackFit_FixedPzPiZero); // fMesonCurIntRangeBackFit no array ?
                fMesonYieldsBackFit_FixedPzPiZero[iPt]      = fYields; // change this part because it is written each time of the loop -> TODO Change to array or move outside of loop
                fMesonYieldsBackFitError_FixedPzPiZero[iPt] = fYieldsError;
            }
            // Integrate the signal histo
            fFileDataLog<< endl <<"Signal histo "<< nameIntRange[k].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRange[k]); //cout<<"2. call of IntegrateHistoInvMassStream"<<endl;

            fMesonYields[k][iPt]            = fYields;
            fMesonYieldsError[k][iPt]       = fYieldsError;
            fFileDataLog << "Integrated value signal histo: \t" << fYields <<"+-" <<fYieldsError <<endl;

            // Integrate the signal histo invMass pi0 subtracted
            fFileDataLog<< endl <<"Signal histo "<< nameIntRange[k].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin_SubPiZero[iPt], fMesonCurIntRange_SubPiZero[k]); //cout<<"2. call of IntegrateHistoInvMassStream"<<endl;

            fMesonYields_SubPiZero[k][iPt]            = fYields;
            fMesonYieldsError_SubPiZero[k][iPt]       = fYieldsError;
            fFileDataLog << "Integrated value signal (pi0 subtracted) histo: \t" << fYields <<"+-" <<fYieldsError <<endl;

            // Integrate the signal histo pz of pi0 fixed
            fFileDataLog<< endl <<"Signal histo "<< nameIntRange[k].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin_FixedPzPiZero[iPt], fMesonCurIntRange_FixedPzPiZero[k]); //cout<<"2. call of IntegrateHistoInvMassStream"<<endl;

            fMesonYields_FixedPzPiZero[k][iPt]            = fYields;
            fMesonYieldsError_FixedPzPiZero[k][iPt]       = fYieldsError;
            fFileDataLog << "Integrated value signal (pi0 pz fixed) histo: \t" << fYields <<"+-" <<fYieldsError <<endl;

            if(fIsMC){
                fFileDataLog<< endl <<"True histo normal range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                fFitTrueSignalInvMassPtBin[iPt]=0x00;
                if(fCrysFitting==0){
                    fFileErrLog << "Using exp fit"<<endl;
                    FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntRange,iPt,kFALSE);
                } else {
                    fFileErrLog << "Using Crystal Ball function"<<endl;
                    FitCBSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntRange,iPt,kFALSE,Form("CBFitFuncMCTrueBin%02d",iPt),0);
                }

                if (fHistoMappingTrueMesonInvMassPtBins[iPt]->GetEntries() !=0){
                    fFitTrueSignalInvMassPtBin[iPt]                  = fFitReco;
                    if (fFitTrueSignalInvMassPtBin[iPt] != 0x00){
                        fMesonTrueMass[iPt]                          = fFitTrueSignalInvMassPtBin[iPt]->GetParameter(1);
                        fMesonTrueMassError[iPt]                     = fFitTrueSignalInvMassPtBin[iPt]->GetParError(1);
                        CalculateFWHM(fFitTrueSignalInvMassPtBin[iPt],fMesonFitRange);
                        fMesonTrueFWHM[iPt]                          = fFWHMFunc;
                        fMesonTrueFWHMError[iPt]                     = fFWHMFuncError;
                        fFileDataLog << "TrueFWHM \t" << fMesonTrueFWHM[iPt] << "\t +-" << fMesonTrueFWHMError[iPt] << endl;

                        fMesonTrueIntRange[1][0]                     = fMesonTrueMass[iPt] + fMesonIntDeltaRangeWide[0];
                        fMesonTrueIntRange[2][0]                     = fMesonTrueMass[iPt] + fMesonIntDeltaRangeNarrow[0];
                        fMesonTrueIntRange[0][1]                     = fMesonTrueMass[iPt] + fMesonIntDeltaRange[1] ;
                        fMesonTrueIntRange[1][1]                     = fMesonTrueMass[iPt] + fMesonIntDeltaRangeWide[1];
                        fMesonTrueIntRange[2][1]                     = fMesonTrueMass[iPt] + fMesonIntDeltaRangeNarrow[1];
                    } else {
                        fMesonTrueMass[iPt]                          = 0.;
                        fMesonTrueMassError[iPt]                     = 1.;
                        fMesonTrueFWHM[iPt]                          = 0.;
                        fMesonTrueFWHMError[iPt]                     = 0.;
                        fMesonTrueIntRange[0][0]                     = fMesonMassExpect + fMesonIntDeltaRange[0];
                        fMesonTrueIntRange[1][0]                     = fMesonMassExpect + fMesonIntDeltaRangeWide[0];
                        fMesonTrueIntRange[2][0]                     = fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
                        fMesonTrueIntRange[0][1]                     = fMesonMassExpect + fMesonIntDeltaRange[1];
                        fMesonTrueIntRange[1][1]                     = fMesonMassExpect + fMesonIntDeltaRangeWide[1];
                        fMesonTrueIntRange[2][1]                     = fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];
                    }
                }

                fFitTrueSignalInvMassPtReweightedBin[iPt]=0x00;

                if(fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]){
                    cout << "Using exp fit"<<endl;
                    fFileErrLog << "Using exp fit"<<endl;
                    FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtReweightedBins[iPt], fMesonIntDeltaRange,iPt,kFALSE);

                    if (fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->GetEntries() !=0){
                        fFitTrueSignalInvMassPtReweightedBin[iPt]   = fFitReco;

                        if (fFitTrueSignalInvMassPtReweightedBin[iPt] != 0x00){
                            fMesonTrueMassReweighted[iPt]           = fFitTrueSignalInvMassPtReweightedBin[iPt]->GetParameter(1);
                            fMesonTrueMassReweightedError[iPt]      = fFitTrueSignalInvMassPtReweightedBin[iPt]->GetParError(1);
                            CalculateFWHM(fFitTrueSignalInvMassPtReweightedBin[iPt],fMesonFitRange);
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

                fFileDataLog<< endl <<"True histo " << nameIntRange[k].Data() << " range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonTrueIntRange[k]);
                fMesonTrueYields[k][iPt]                        = fYields;
                fMesonTrueYieldsError[k][iPt]                   = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

                // Yields reweighted in different ranges
                fFileDataLog<< endl <<"True histo " << nameIntRange[k].Data() << " range reweighted" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                IntegrateHistoInvMassStream( fHistoMappingTrueMesonInvMassPtReweightedBins[iPt], fMesonTrueIntReweightedRange[k]);
                fMesonTrueYieldsReweighted[k][iPt]              = fYields;
                fMesonTrueYieldsReweightedError[k][iPt]             = fYieldsError;
                fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError<<endl;

                if (k == 0 ){
                    if( (fGGYields[0][iPt] - fMesonTrueYields[0][iPt]) > 0) {
                        fMesonTrueSB[iPt]           = fMesonTrueYields[0][iPt] / ( fGGYields[0][iPt] - fMesonTrueYields[0][iPt] );
                        fMesonTrueSign[iPt]         = fMesonTrueYields[0][iPt] / pow( ( fGGYields[0][iPt] - fMesonTrueYields[0][iPt] ) , 0.5);
                        fMesonTrueSBError[iPt]      = 0;
                        fMesonTrueSignError[iPt]    = 0;
                    }
                }
            }

            fFileDataLog << "Residual Background leftover " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                                << fMesonYieldsResidualBckFunc[k][iPt] << "\t +- \t" << fMesonYieldsResidualBckFuncError[k][iPt] << endl<< endl;

            fTotalBckYields[k][iPt]                         = fBckYields[k][iPt] + fMesonYieldsResidualBckFunc[k][iPt];
            fTotalBckYields_SubPiZero[k][iPt]               = fBckYields_SubPiZero[k][iPt] + fMesonYieldsResidualBckFunc_SubPiZero[k][iPt];
            fTotalBckYields_FixedPzPiZero[k][iPt]           = fBckYields_FixedPzPiZero[k][iPt] + fMesonYieldsResidualBckFunc_FixedPzPiZero[k][iPt];
            fTotalBckYieldsError[k][iPt]                    = pow(fBckYieldsError[k][iPt]*fBckYieldsError[k][iPt] + fMesonYieldsResidualBckFuncError[k][iPt]*fMesonYieldsResidualBckFuncError[k][iPt],0.5);
            fTotalBckYieldsError_SubPiZero[k][iPt]          = pow(fBckYieldsError_SubPiZero[k][iPt]*fBckYieldsError_SubPiZero[k][iPt] + fMesonYieldsResidualBckFuncError_SubPiZero[k][iPt]*fMesonYieldsResidualBckFuncError_SubPiZero[k][iPt],0.5);
            fTotalBckYieldsError_FixedPzPiZero[k][iPt]     = pow(fBckYieldsError_FixedPzPiZero[k][iPt]*fBckYieldsError_FixedPzPiZero[k][iPt] + fMesonYieldsResidualBckFuncError_FixedPzPiZero[k][iPt]*fMesonYieldsResidualBckFuncError_FixedPzPiZero[k][iPt],0.5);

            fFileDataLog << "Total Background " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fTotalBckYields[k][iPt] << "\t +- \t" << fTotalBckYieldsError[k][iPt] << endl<< endl;
            fFileDataLog << "Background " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fBckYields[k][iPt] << "\t +- \t" << fBckYieldsError[k][iPt] << endl<< endl;

            fFileDataLog << "Total Background (Pi0mass subtracted)" << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fTotalBckYields_SubPiZero[k][iPt] << "\t +- \t" << fTotalBckYieldsError_SubPiZero[k][iPt] << endl<< endl;
            fFileDataLog << "Background (pi0mass subtracted)" << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fBckYields_SubPiZero[k][iPt] << "\t +- \t" << fBckYieldsError_SubPiZero[k][iPt] << endl<< endl;

            fFileDataLog << "Total Background (Pi0 pz fixed)" << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fTotalBckYields_FixedPzPiZero[k][iPt] << "\t +- \t" << fTotalBckYieldsError_FixedPzPiZero[k][iPt] << endl<< endl;
            fFileDataLog << "Background (Pi0 pz fixed)" << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fBckYields_FixedPzPiZero[k][iPt] << "\t +- \t" << fBckYieldsError_FixedPzPiZero[k][iPt] << endl<< endl;

            fMesonYieldsCorResidualBckFunc[k][iPt]          = fMesonYields[k][iPt]- fMesonYieldsResidualBckFunc[k][iPt];
            fMesonYieldsCorResidualBckFuncError[k][iPt]     = pow(( fMesonYieldsError[k][iPt]*fMesonYieldsError[k][iPt] +
                                                                fMesonYieldsResidualBckFuncError[k][iPt]*fMesonYieldsResidualBckFuncError[k][iPt]),0.5);
            fMesonYieldsPerEvent[k][iPt]                    = fMesonYieldsCorResidualBckFunc[k][iPt]/fNEvents;
            fMesonYieldsPerEventError[k][iPt]               = fMesonYieldsCorResidualBckFuncError[k][iPt]/fNEvents;

            // Calculation for InvMass pi0 subtracted
            fMesonYieldsCorResidualBckFunc_SubPiZero[k][iPt]          = fMesonYields_SubPiZero[k][iPt]- fMesonYieldsResidualBckFunc_SubPiZero[k][iPt];
            fMesonYieldsCorResidualBckFuncError_SubPiZero[k][iPt]     = pow(( fMesonYieldsError_SubPiZero[k][iPt]*fMesonYieldsError_SubPiZero[k][iPt] +
                                                                fMesonYieldsResidualBckFuncError[k][iPt]*fMesonYieldsResidualBckFuncError[k][iPt]),0.5);
            fMesonYieldsPerEvent_SubPiZero[k][iPt]                    = fMesonYieldsCorResidualBckFunc_SubPiZero[k][iPt]/fNEvents;
            fMesonYieldsPerEventError_SubPiZero[k][iPt]               = fMesonYieldsCorResidualBckFuncError_SubPiZero[k][iPt]/fNEvents;

            // Calculation for InvMass pz of pi0 fixes
            fMesonYieldsCorResidualBckFunc_FixedPzPiZero[k][iPt]          = fMesonYields_FixedPzPiZero[k][iPt]- fMesonYieldsResidualBckFunc_FixedPzPiZero[k][iPt];
            fMesonYieldsCorResidualBckFuncError_FixedPzPiZero[k][iPt]     = pow(( fMesonYieldsError_FixedPzPiZero[k][iPt]*fMesonYieldsError_FixedPzPiZero[k][iPt] +
                                                                fMesonYieldsResidualBckFuncError[k][iPt]*fMesonYieldsResidualBckFuncError[k][iPt]),0.5);
            fMesonYieldsPerEvent_FixedPzPiZero[k][iPt]                    = fMesonYieldsCorResidualBckFunc_FixedPzPiZero[k][iPt]/fNEvents;
            fMesonYieldsPerEventError_FixedPzPiZero[k][iPt]               = fMesonYieldsCorResidualBckFuncError_SubPiZero[k][iPt]/fNEvents;

            // check this part, maybe move to an array later on. Dirty fix -> only do it once for k == 0)
            if (k == 0){
                fMesonYieldsCorResidualBckFuncBackFit[iPt]          = fMesonYieldsBackFit[iPt]- fMesonYieldsResidualBckFuncBackFit[iPt];
                fMesonYieldsCorResidualBckFuncBackFitError[iPt]     =
                pow((fMesonYieldsBackFitError[iPt]*fMesonYieldsBackFitError[iPt]+
                    fMesonYieldsResidualBckFuncBackFitError[iPt]*fMesonYieldsResidualBckFuncBackFitError[iPt]),0.5);
                fMesonYieldsPerEventBackFit[iPt]                    = fMesonYieldsCorResidualBckFuncBackFit[iPt]/fNEvents;
                fMesonYieldsPerEventBackFitError[iPt]               = fMesonYieldsCorResidualBckFuncBackFitError[iPt]/fNEvents;

                // Calculation for InvMass pi0 subtracted
                fMesonYieldsCorResidualBckFuncBackFit_SubPiZero[iPt]          = fMesonYieldsBackFit_SubPiZero[iPt]- fMesonYieldsResidualBckFuncBackFit_SubPiZero[iPt];
                fMesonYieldsCorResidualBckFuncBackFitError_SubPiZero[iPt]     =
                pow((fMesonYieldsBackFitError_SubPiZero[iPt]*fMesonYieldsBackFitError_SubPiZero[iPt]+
                    fMesonYieldsResidualBckFuncBackFitError_SubPiZero[iPt]*fMesonYieldsResidualBckFuncBackFitError_SubPiZero[iPt]),0.5);
                fMesonYieldsPerEventBackFit_SubPiZero[iPt]                    = fMesonYieldsCorResidualBckFuncBackFit_SubPiZero[iPt]/fNEvents;
                fMesonYieldsPerEventBackFitError_SubPiZero[iPt]               = fMesonYieldsCorResidualBckFuncBackFitError_SubPiZero[iPt]/fNEvents;


                // Calculation for pz of pi0 fixed
                fMesonYieldsCorResidualBckFuncBackFit_FixedPzPiZero[iPt]          = fMesonYieldsBackFit_FixedPzPiZero[iPt]- fMesonYieldsResidualBckFuncBackFit_FixedPzPiZero[iPt];
                fMesonYieldsCorResidualBckFuncBackFitError_FixedPzPiZero[iPt]     =
                pow((fMesonYieldsBackFitError_FixedPzPiZero[iPt]*fMesonYieldsBackFitError_FixedPzPiZero[iPt]+
                    fMesonYieldsResidualBckFuncBackFitError_FixedPzPiZero[iPt]*fMesonYieldsResidualBckFuncBackFitError_FixedPzPiZero[iPt]),0.5);
                fMesonYieldsPerEventBackFit_FixedPzPiZero[iPt]                    = fMesonYieldsCorResidualBckFuncBackFit_FixedPzPiZero[iPt]/fNEvents;
                fMesonYieldsPerEventBackFitError_FixedPzPiZero[iPt]               = fMesonYieldsCorResidualBckFuncBackFitError_SubPiZero[iPt]/fNEvents;
            }


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
                    cout << "k =" << k << endl;
                    FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], intRange, iPt, kFALSE,0);
                    fMesonYieldsResidualBckFunc[k][iPt]         = fIntLinearBck;
                    fMesonYieldsResidualBckFuncError[k+3][iPt]    = fIntLinearBckError;

                    // Fitting for InvMass pi0 subtracted
                    FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin_SubPiZero[iPt], intRange, iPt, kFALSE,1);
                    fMesonYieldsResidualBckFunc_SubPiZero[k][iPt]           = fIntLinearBck;
                    fMesonYieldsResidualBckFuncError_SubPiZero[k+3][iPt]    = fIntLinearBckError;

                    // Fitting for InvMass pz of pi0 fixed
                    FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin_FixedPzPiZero[iPt], intRange, iPt, kFALSE,2);
                    fMesonYieldsResidualBckFunc_FixedPzPiZero[k][iPt]           = fIntLinearBck;
                    fMesonYieldsResidualBckFuncError_FixedPzPiZero[k+3][iPt]    = fIntLinearBckError;
                } else {
                    fMesonYieldsResidualBckFunc[k][iPt]         = 0;
                    fMesonYieldsResidualBckFuncError[k][iPt]    = 0;
                    fMesonYieldsResidualBckFunc_SubPiZero[k][iPt]         = 0;
                    fMesonYieldsResidualBckFuncError_SubPiZero[k][iPt]    = 0;
                    fMesonYieldsResidualBckFunc_FixedPzPiZero[k][iPt]         = 0;
                    fMesonYieldsResidualBckFuncError_FixedPzPiZero[k][iPt]    = 0;
                }
            }
            fFileDataLog << "Residual Background leftover " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                            << fMesonYieldsResidualBckFunc[k][iPt] << "\t +- \t" << fMesonYieldsResidualBckFuncError[k][iPt] << endl<< endl;
            fFileDataLog << "Residual Background leftover (pi0 mass subtracted) " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                            << fMesonYieldsResidualBckFunc_SubPiZero[k][iPt] << "\t +- \t" << fMesonYieldsResidualBckFuncError_SubPiZero[k][iPt] << endl<< endl;
            fFileDataLog << "Residual Background leftover (pz of pi0 fixed)" << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                            << fMesonYieldsResidualBckFunc_SubPiZero[k][iPt] << "\t +- \t" << fMesonYieldsResidualBckFuncError_SubPiZero[k][iPt] << endl<< endl;


            fTotalBckYields[k][iPt]                         = fBckYields[k][iPt] + fMesonYieldsResidualBckFunc[k][iPt];
            fTotalBckYieldsError[k][iPt]                    = pow(fBckYieldsError[k][iPt]*fBckYieldsError[k][iPt] + fMesonYieldsResidualBckFuncError[k][iPt]*fMesonYieldsResidualBckFuncError[k][iPt],0.5);

            fTotalBckYields_SubPiZero[k][iPt]               = fBckYields_SubPiZero[k][iPt] + fMesonYieldsResidualBckFunc_SubPiZero[k][iPt];
            fTotalBckYieldsError_SubPiZero[k][iPt]                    = pow(fBckYieldsError_SubPiZero[k][iPt]*fBckYieldsError_SubPiZero[k][iPt] + fMesonYieldsResidualBckFuncError_SubPiZero[k][iPt]*fMesonYieldsResidualBckFuncError_SubPiZero[k][iPt],0.5);

            fTotalBckYields_FixedPzPiZero[k][iPt]                         = fBckYields_FixedPzPiZero[k][iPt] + fMesonYieldsResidualBckFunc_FixedPzPiZero[k][iPt];
            fTotalBckYieldsError_FixedPzPiZero[k][iPt]                    = pow(fBckYieldsError_FixedPzPiZero[k][iPt]*fBckYieldsError_FixedPzPiZero[k][iPt] + fMesonYieldsResidualBckFuncError_FixedPzPiZero[k][iPt]*fMesonYieldsResidualBckFuncError_FixedPzPiZero[k][iPt],0.5);

            fFileDataLog << "Total Background " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fTotalBckYields[k][iPt] << "\t +- \t" << fTotalBckYieldsError[k][iPt] << endl<< endl;
            fFileDataLog << "Background " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fBckYields[k][iPt] << "\t +- \t" << fBckYieldsError[k][iPt] << endl<< endl;

            fFileDataLog << "Total Background (pi0 mass subtracted) " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fTotalBckYields_SubPiZero[k][iPt] << "\t +- \t" << fTotalBckYieldsError_SubPiZero[k][iPt] << endl<< endl;
            fFileDataLog << "Background (pi0 mass subtracted) " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fBckYields_SubPiZero[k][iPt] << "\t +- \t" << fBckYieldsError_SubPiZero[k][iPt] << endl<< endl;

            fFileDataLog << "Total Background (pz of pi0 fixed) " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fTotalBckYields_FixedPzPiZero[k][iPt] << "\t +- \t" << fTotalBckYieldsError_FixedPzPiZero[k][iPt] << endl<< endl;
            fFileDataLog << "Background (piz of pizero fixed) " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fBckYields_FixedPzPiZero[k][iPt] << "\t +- \t" << fBckYieldsError_FixedPzPiZero[k][iPt] << endl<< endl;

            fMesonYieldsCorResidualBckFunc[k][iPt]          = fMesonYields[k][iPt]- fMesonYieldsResidualBckFunc[k][iPt];
            fMesonYieldsCorResidualBckFuncError[k][iPt]     = pow(( fMesonYieldsError[k][iPt]*fMesonYieldsError[k][iPt] +
                                                                fMesonYieldsResidualBckFuncError[k][iPt]*fMesonYieldsResidualBckFuncError[k][iPt]),0.5);
            fMesonYieldsPerEvent[k][iPt]                    = fMesonYieldsCorResidualBckFunc[k][iPt]/fNEvents;
            fMesonYieldsPerEventError[k][iPt]               = fMesonYieldsCorResidualBckFuncError[k][iPt]/fNEvents;

            // pi0 mass subtracted
            fMesonYieldsCorResidualBckFunc_SubPiZero[k][iPt]          = fMesonYields_SubPiZero[k][iPt]- fMesonYieldsResidualBckFunc_SubPiZero[k][iPt];
            fMesonYieldsCorResidualBckFuncError_SubPiZero[k][iPt]     = pow(( fMesonYieldsError_SubPiZero[k][iPt]*fMesonYieldsError_SubPiZero[k][iPt] +
                                                                fMesonYieldsResidualBckFuncError_SubPiZero[k][iPt]*fMesonYieldsResidualBckFuncError_SubPiZero[k][iPt]),0.5);
            fMesonYieldsPerEvent_SubPiZero[k][iPt]                    = fMesonYieldsCorResidualBckFunc_SubPiZero[k][iPt]/fNEvents;
            fMesonYieldsPerEventError_SubPiZero[k][iPt]               = fMesonYieldsCorResidualBckFuncError_SubPiZero[k][iPt]/fNEvents;

            // pz of pi0 fixed
            fMesonYieldsCorResidualBckFunc_FixedPzPiZero[k][iPt]          = fMesonYields_FixedPzPiZero[k][iPt]- fMesonYieldsResidualBckFunc_FixedPzPiZero[k][iPt];
            fMesonYieldsCorResidualBckFuncError_FixedPzPiZero[k][iPt]     = pow(( fMesonYieldsError_FixedPzPiZero[k][iPt]*fMesonYieldsError_FixedPzPiZero[k][iPt] +
                                                                fMesonYieldsResidualBckFuncError_FixedPzPiZero[k][iPt]*fMesonYieldsResidualBckFuncError_FixedPzPiZero[k][iPt]),0.5);
            fMesonYieldsPerEvent_FixedPzPiZero[k][iPt]                    = fMesonYieldsCorResidualBckFunc_FixedPzPiZero[k][iPt]/fNEvents;
            fMesonYieldsPerEventError_FixedPzPiZero[k][iPt]               = fMesonYieldsCorResidualBckFuncError_FixedPzPiZero[k][iPt]/fNEvents;
        }


        //////////////////////////////// Start Analysis with  Normalization at the left of the Meson Peak

        // Function to subtract GG minus Bck
        cout << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
        ProcessEMLeftRight( fHistoMappingGGInvMassPtBin[iPt], fHistoMappingBackInvMassPtBin[0][iPt], fBGFitRangeLeft, fBGFitRange); // only done for added background so far
        fHistoMappingSignalInvMassLeftPtBin[iPt] = fSignal;
        fHistoMappingBackNormInvMassLeftPtBin[iPt] = fBckNorm;

        // Function to subtract GG minus Bck (pi0 mass subtracted)
        cout << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
        ProcessEMLeftRight( fHistoMappingGGInvMassPtBin_SubPiZero[iPt], fHistoMappingBackInvMassPtBin_SubPiZero[0][iPt], fBGFitRangeLeft_SubPiZero, fBGFitRange_SubPiZero); // only done for added background so far
        fHistoMappingSignalInvMassLeftPtBin_SubPiZero[iPt] = fSignal;
        fHistoMappingBackNormInvMassLeftPtBin_SubPiZero[iPt] = fBckNorm;

        // Function to subtract GG minus Bck (pi0 pz fixed)
        cout << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
        ProcessEMLeftRight( fHistoMappingGGInvMassPtBin_FixedPzPiZero[iPt], fHistoMappingBackInvMassPtBin_FixedPzPiZero[0][iPt], fBGFitRangeLeft_FixedPzPiZero, fBGFitRange_FixedPzPiZero); // only done for added background so far
        fHistoMappingSignalInvMassLeftPtBin_FixedPzPiZero[iPt] = fSignal;
        fHistoMappingBackNormInvMassLeftPtBin_FixedPzPiZero[iPt] = fBckNorm;

        // Fitting the subtracted spectra
        fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "normal range/left normalization" << endl;

        fFitInvMassLeftPtBin[iPt]               = 0x00;
        fFitInvMassLeftPtBin_SubPiZero[iPt]     = 0x00;
        fFitInvMassLeftPtBin_FixedPzPiZero[iPt] = 0x00;
        if(fCrysFitting==0){
            fFileErrLog << "Using exp fit"<<endl;
            FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE,0);
            fFitInvMassLeftPtBin[iPt]                   = fFitReco;
            fFitSignalPeakPosInvMassLeftPtBin[iPt]      = fFitGausExp;
            fFitBckInvMassLeftPtBin[iPt]                = fFitLinearBck;
            fMesonYieldsResidualBckFunc[3][iPt]         = fIntLinearBck;
            fMesonYieldsResidualBckFuncError[3][iPt]    = fIntLinearBckError;

            // Fitting pi0 inv mass subtracted
            if(fHistoMappingSignalInvMassLeftPtBin_SubPiZero[iPt] != NULL){
              fFileErrLog << "Using exp fit (pi0 inv mass subtracted"<<endl;
              FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin_SubPiZero[iPt], fMesonIntDeltaRange,iPt,kFALSE,1);
              fFitInvMassLeftPtBin_SubPiZero[iPt]                   = fFitReco;
              fFitSignalPeakPosInvMassLeftPtBin_SubPiZero[iPt]      = fFitGausExp;
              fFitBckInvMassLeftPtBin_SubPiZero[iPt]                = fFitLinearBck;
              fMesonYieldsResidualBckFunc_SubPiZero[3][iPt]         = fIntLinearBck;
              fMesonYieldsResidualBckFuncError_SubPiZero[3][iPt]    = fIntLinearBckError;
            }

            // Fitting with pz of pi0 fixed
            if(fHistoMappingSignalInvMassLeftPtBin_FixedPzPiZero[iPt] != NULL){
              fFileErrLog << "Using exp fit (pi0 inv mass subtracted"<<endl;
              FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin_FixedPzPiZero[iPt], fMesonIntDeltaRange,iPt,kFALSE,2);
              fFitInvMassLeftPtBin_FixedPzPiZero[iPt]                   = fFitReco;
              fFitSignalPeakPosInvMassLeftPtBin_FixedPzPiZero[iPt]      = fFitGausExp;
              fFitBckInvMassLeftPtBin_FixedPzPiZero[iPt]                = fFitLinearBck;
              fMesonYieldsResidualBckFunc_FixedPzPiZero[3][iPt]         = fIntLinearBck;
              fMesonYieldsResidualBckFuncError_FixedPzPiZero[3][iPt]    = fIntLinearBckError;
            }

        } else {
            fFileErrLog << "Using Crystal Ball function"<<endl;
            FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE, Form("CBFitFuncLeftBin%02d",iPt),0);
            fFitInvMassLeftPtBin[iPt]                   = fFitReco;
            fFitSignalPeakPosInvMassLeftPtBin[iPt]      = fFitGausExp;
            fFitBckInvMassLeftPtBin[iPt]                = fFitLinearBck;
            fMesonYieldsResidualBckFunc[3][iPt]         = fIntLinearBck;
            fMesonYieldsResidualBckFuncError[3][iPt]    = fIntLinearBckError;

            // Fitting with crystal ball and pi0 inv mass subtracted
            if(fHistoMappingSignalInvMassLeftPtBin_SubPiZero[iPt]!=NULL){
              FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin_SubPiZero[iPt], fMesonIntDeltaRange,iPt,kFALSE, Form("CBFitFuncLeft_SubPiZero_Bin%02d",iPt),1);
              fFitInvMassLeftPtBin_SubPiZero[iPt]                   = fFitReco;
              fFitSignalPeakPosInvMassLeftPtBin_SubPiZero[iPt]      = fFitGausExp;
              fFitBckInvMassLeftPtBin_SubPiZero[iPt]                = fFitLinearBck;
              fMesonYieldsResidualBckFunc_SubPiZero[3][iPt]         = fIntLinearBck;
              fMesonYieldsResidualBckFuncError_SubPiZero[3][iPt]    = fIntLinearBckError;
            }

            // Fitting with crystal ball and pi0 inv mass subtracted
            if(fHistoMappingSignalInvMassLeftPtBin_FixedPzPiZero[iPt]!=NULL){
              FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin_FixedPzPiZero[iPt], fMesonIntDeltaRange,iPt,kFALSE, Form("CBFitFuncLeft_FixedPzPiZero_Bin%02d",iPt),1);
              fFitInvMassLeftPtBin_FixedPzPiZero[iPt]                   = fFitReco;
              fFitSignalPeakPosInvMassLeftPtBin_FixedPzPiZero[iPt]      = fFitGausExp;
              fFitBckInvMassLeftPtBin_FixedPzPiZero[iPt]                = fFitLinearBck;
              fMesonYieldsResidualBckFunc_FixedPzPiZero[3][iPt]         = fIntLinearBck;
              fMesonYieldsResidualBckFuncError_FixedPzPiZero[3][iPt]    = fIntLinearBckError;
            }

        }

        CalculateFWHM(fFitInvMassLeftPtBin[iPt],fMesonFitRange);
        fMesonFWHMLeft[iPt]      = fFWHMFunc;
        fMesonFWHMLeftError[iPt] = fFWHMFuncError;

        CalculateFWHM(fFitInvMassLeftPtBin_SubPiZero[iPt],fMesonFitRange_SubPiZero);
        fMesonFWHMLeft_SubPiZero[iPt]      = fFWHMFunc;
        fMesonFWHMLeftError_SubPiZero[iPt] = fFWHMFuncError;

        CalculateFWHM(fFitInvMassLeftPtBin_FixedPzPiZero[iPt],fMesonFitRange_FixedPzPiZero);
        fMesonFWHMLeft_FixedPzPiZero[iPt]      = fFWHMFunc;
        fMesonFWHMLeftError_FixedPzPiZero[iPt] = fFWHMFuncError;

        if (fFitInvMassLeftPtBin[iPt] !=0x00){
            fMesonMassLeft[iPt] = fFitInvMassLeftPtBin[iPt]->GetParameter(1);
            fMesonMassLeftError[iPt] = fFitInvMassLeftPtBin[iPt]->GetParError(1);
            fMesonWidthLeft[iPt] = fFitInvMassLeftPtBin[iPt]->GetParameter(2);
            fMesonWidthLeftError[iPt] = fFitInvMassLeftPtBin[iPt]->GetParError(2);

            fMesonCurIntRange[3][0]         = fMesonMassLeft[iPt] + fMesonIntDeltaRange[0];
            fMesonCurIntRange[4][0]         = fMesonMassLeft[iPt] + fMesonIntDeltaRangeWide[0];
            fMesonCurIntRange[5][0]         = fMesonMassLeft[iPt] + fMesonIntDeltaRangeNarrow[0];
            fMesonCurIntRange[3][1]         = fMesonMassLeft[iPt] + fMesonIntDeltaRange[1];
            fMesonCurIntRange[4][1]         = fMesonMassLeft[iPt] + fMesonIntDeltaRangeWide[1];
            fMesonCurIntRange[5][1]         = fMesonMassLeft[iPt] + fMesonIntDeltaRangeNarrow[1];

            if(fFitInvMassLeftPtBin_SubPiZero !=0x00){
              fMesonMassLeft_SubPiZero[iPt]       = fFitInvMassLeftPtBin_SubPiZero[iPt]->GetParameter(1);
              fMesonMassLeftError_SubPiZero[iPt]  = fFitInvMassLeftPtBin_SubPiZero[iPt]->GetParError(1);
              fMesonWidthLeft_SubPiZero[iPt]      = fFitInvMassLeftPtBin_SubPiZero[iPt]->GetParameter(2);
              fMesonWidthLeftError_SubPiZero[iPt] = fFitInvMassLeftPtBin_SubPiZero[iPt]->GetParError(2);

              fMesonCurIntRange_SubPiZero[3][0]         = fMesonMassLeft_SubPiZero[iPt] + fMesonIntDeltaRange[0];
              fMesonCurIntRange_SubPiZero[4][0]         = fMesonMassLeft_SubPiZero[iPt] + fMesonIntDeltaRangeWide[0];
              fMesonCurIntRange_SubPiZero[5][0]         = fMesonMassLeft_SubPiZero[iPt] + fMesonIntDeltaRangeNarrow[0];
              fMesonCurIntRange_SubPiZero[3][1]         = fMesonMassLeft_SubPiZero[iPt] + fMesonIntDeltaRange[1];
              fMesonCurIntRange_SubPiZero[4][1]         = fMesonMassLeft_SubPiZero[iPt] + fMesonIntDeltaRangeWide[1];
              fMesonCurIntRange_SubPiZero[5][1]         = fMesonMassLeft_SubPiZero[iPt] + fMesonIntDeltaRangeNarrow[1];
            }

            if(fFitInvMassLeftPtBin_FixedPzPiZero !=0x00){
              fMesonMassLeft_FixedPzPiZero[iPt]       = fFitInvMassLeftPtBin_FixedPzPiZero[iPt]->GetParameter(1);
              fMesonMassLeftError_FixedPzPiZero[iPt]  = fFitInvMassLeftPtBin_FixedPzPiZero[iPt]->GetParError(1);
              fMesonWidthLeft_FixedPzPiZero[iPt]      = fFitInvMassLeftPtBin_FixedPzPiZero[iPt]->GetParameter(2);
              fMesonWidthLeftError_FixedPzPiZero[iPt] = fFitInvMassLeftPtBin_FixedPzPiZero[iPt]->GetParError(2);

              fMesonCurIntRange_FixedPzPiZero[3][0]         = fMesonMassLeft_FixedPzPiZero[iPt] + fMesonIntDeltaRange[0];
              fMesonCurIntRange_FixedPzPiZero[4][0]         = fMesonMassLeft_FixedPzPiZero[iPt] + fMesonIntDeltaRangeWide[0];
              fMesonCurIntRange_FixedPzPiZero[5][0]         = fMesonMassLeft_FixedPzPiZero[iPt] + fMesonIntDeltaRangeNarrow[0];
              fMesonCurIntRange_FixedPzPiZero[3][1]         = fMesonMassLeft_FixedPzPiZero[iPt] + fMesonIntDeltaRange[1];
              fMesonCurIntRange_FixedPzPiZero[4][1]         = fMesonMassLeft_FixedPzPiZero[iPt] + fMesonIntDeltaRangeWide[1];
              fMesonCurIntRange_FixedPzPiZero[5][1]         = fMesonMassLeft_FixedPzPiZero[iPt] + fMesonIntDeltaRangeNarrow[1];
            }
        } else {
            fMesonMassLeft[iPt]             = 0.;
            fMesonMassLeftError[iPt]        = 0.;
            fMesonWidthLeft[iPt]            = 0.;
            fMesonWidthLeftError[iPt]       = 0.;
            fMesonCurIntRange[3][0]         = fMesonMassExpect + fMesonIntDeltaRange[0];
            fMesonCurIntRange[4][0]         = fMesonMassExpect + fMesonIntDeltaRangeWide[0];
            fMesonCurIntRange[5][0]         = fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
            fMesonCurIntRange[3][1]         = fMesonMassExpect + fMesonIntDeltaRange[1];
            fMesonCurIntRange[4][1]         = fMesonMassExpect + fMesonIntDeltaRangeWide[1];
            fMesonCurIntRange[5][1]         = fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];

            fMesonMassLeft_SubPiZero[iPt]             = 0.;
            fMesonMassLeftError_SubPiZero[iPt]        = 0.;
            fMesonWidthLeft_SubPiZero[iPt]            = 0.;
            fMesonWidthLeftError_SubPiZero[iPt]       = 0.;
            fMesonCurIntRange_SubPiZero[3][0]         = fMesonMassExpect - 0.134 + fMesonIntDeltaRange[0];
            fMesonCurIntRange_SubPiZero[4][0]         = fMesonMassExpect - 0.134 + fMesonIntDeltaRangeWide[0];
            fMesonCurIntRange_SubPiZero[5][0]         = fMesonMassExpect - 0.134 + fMesonIntDeltaRangeNarrow[0];
            fMesonCurIntRange_SubPiZero[3][1]         = fMesonMassExpect - 0.134 + fMesonIntDeltaRange[1];
            fMesonCurIntRange_SubPiZero[4][1]         = fMesonMassExpect - 0.134 + fMesonIntDeltaRangeWide[1];
            fMesonCurIntRange_SubPiZero[5][1]         = fMesonMassExpect - 0.134 + fMesonIntDeltaRangeNarrow[1];

            fMesonMassLeft_FixedPzPiZero[iPt]             = 0.;
            fMesonMassLeftError_FixedPzPiZero[iPt]        = 0.;
            fMesonWidthLeft_FixedPzPiZero[iPt]            = 0.;
            fMesonWidthLeftError_FixedPzPiZero[iPt]       = 0.;
            fMesonCurIntRange_FixedPzPiZero[3][0]         = fMesonMassExpect + fMesonIntDeltaRange[0];
            fMesonCurIntRange_FixedPzPiZero[4][0]         = fMesonMassExpect + fMesonIntDeltaRangeWide[0];
            fMesonCurIntRange_FixedPzPiZero[5][0]         = fMesonMassExpect + fMesonIntDeltaRangeNarrow[0];
            fMesonCurIntRange_FixedPzPiZero[3][1]         = fMesonMassExpect + fMesonIntDeltaRange[1];
            fMesonCurIntRange_FixedPzPiZero[4][1]         = fMesonMassExpect + fMesonIntDeltaRangeWide[1];
            fMesonCurIntRange_FixedPzPiZero[5][1]         = fMesonMassExpect + fMesonIntDeltaRangeNarrow[1];
        }

        // Integrate the bck histo
        for (Int_t k = 0; k < 3; k++){
            // Normal InvMass
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurIntRange[k+3]);
            fBckYields[k+3][iPt]              = fYields;
            fBckYieldsError[k+3][iPt]         = fYieldsError;

            fFileDataLog<< endl <<"Signal histo " << nameIntRange[k+3].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurIntRange[k+3]);
            fMesonYields[k+3][iPt]            = fYields;
            fMesonYieldsError[k+3][iPt]       = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;

            // InvMass pi0 inv mass subtracted
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin_SubPiZero[iPt], fMesonCurIntRange_SubPiZero[k+3]);
            fBckYields_SubPiZero[k+3][iPt]              = fYields;
            fBckYieldsError_SubPiZero[k+3][iPt]         = fYieldsError;

            fFileDataLog<< endl <<"Signal histo (pi0 inv mass subtracted)" << nameIntRange[k+3].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin_SubPiZero[iPt], fMesonCurIntRange_SubPiZero[k+3]);
            fMesonYields_SubPiZero[k+3][iPt]            = fYields;
            fMesonYieldsError_SubPiZero[k+3][iPt]       = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;

            // InvMass fixed pz of pi0
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin_FixedPzPiZero[iPt], fMesonCurIntRange_FixedPzPiZero[k+3]);
            fBckYields_FixedPzPiZero[k+3][iPt]              = fYields;
            fBckYieldsError_FixedPzPiZero[k+3][iPt]         = fYieldsError;

            fFileDataLog<< endl <<"Signal histo (pi0 with fixed pz)" << nameIntRange[k+3].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin_FixedPzPiZero[iPt], fMesonCurIntRange_FixedPzPiZero[k+3]);
            fMesonYields_FixedPzPiZero[k+3][iPt]            = fYields;
            fMesonYieldsError_FixedPzPiZero[k+3][iPt]       = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
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
                    FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], intRange, iPt, kFALSE,0);
                    fMesonYieldsResidualBckFunc[k+3][iPt]         = fIntLinearBck;
                    fMesonYieldsResidualBckFuncError[k+3][iPt]    = fIntLinearBckError;

                    FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin_SubPiZero[iPt], intRange, iPt, kFALSE,1);
                    fMesonYieldsResidualBckFunc_SubPiZero[k+3][iPt]         = fIntLinearBck;
                    fMesonYieldsResidualBckFuncError_SubPiZero[k+3][iPt]    = fIntLinearBckError;

                    FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin_FixedPzPiZero[iPt], intRange, iPt, kFALSE,2);
                    fMesonYieldsResidualBckFunc_FixedPzPiZero[k+3][iPt]         = fIntLinearBck;
                    fMesonYieldsResidualBckFuncError_FixedPzPiZero[k+3][iPt]    = fIntLinearBckError;
                } else {
                    fMesonYieldsResidualBckFunc[k+3][iPt]         = 0;
                    fMesonYieldsResidualBckFuncError[k+3][iPt]    = 0;

                    fMesonYieldsResidualBckFunc_SubPiZero[k+3][iPt]         = 0;
                    fMesonYieldsResidualBckFuncError_SubPiZero[k+3][iPt]    = 0;

                    fMesonYieldsResidualBckFunc_FixedPzPiZero[k+3][iPt]         = 0;
                    fMesonYieldsResidualBckFuncError_FixedPzPiZero[k+3][iPt]    = 0;
                }
            }

            // Normal Inv mass
            fFileDataLog << "Residual Background leftover " << nameIntRange[k+3].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fMesonYieldsResidualBckFunc[k+3][iPt] << "\t +- \t" << fMesonYieldsResidualBckFuncError[k+3][iPt] << endl<< endl;

            fTotalBckYields[k+3][iPt]                       = fBckYields[k+3][iPt] + fMesonYieldsResidualBckFunc[k+3][iPt];
            fTotalBckYieldsError[k+3][iPt]                  = pow(fBckYieldsError[k+3][iPt]*fBckYieldsError[k+3][iPt] + fMesonYieldsResidualBckFuncError[k+3][iPt]*fMesonYieldsResidualBckFuncError[k+3][iPt],0.5);
            fFileDataLog << "Total Background " << nameIntRange[k+3].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fTotalBckYields[k+3][iPt] << "\t +- \t" << fTotalBckYieldsError[k+3][iPt] << endl<< endl;
            fFileDataLog << "Background " << nameIntRange[k+3].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fBckYields[k+3][iPt] << "\t +- \t" << fBckYieldsError[k+3][iPt] << endl<< endl;

            fMesonYieldsCorResidualBckFunc[k+3][iPt]        = fMesonYields[k+3][iPt]- fMesonYieldsResidualBckFunc[k+3][iPt];
            fMesonYieldsCorResidualBckFuncError[k+3][iPt]   = pow(( fMesonYieldsError[k+3][iPt]*fMesonYieldsError[k+3][iPt] +
                                                                    fMesonYieldsResidualBckFuncError[k+3][iPt]*fMesonYieldsResidualBckFuncError[k+3][iPt]),0.5);
            fMesonYieldsPerEvent[k+3][iPt]                  = fMesonYieldsCorResidualBckFunc[k+3][iPt]/fNEvents;
            fMesonYieldsPerEventError[k+3][iPt]             = fMesonYieldsCorResidualBckFuncError[k+3][iPt]/fNEvents;

            // Inv mass pi0 mass subtracted
            fFileDataLog << "Residual Background leftover (pi0 subtracted" << nameIntRange[k+3].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fMesonYieldsResidualBckFunc_SubPiZero[k+3][iPt] << "\t +- \t" << fMesonYieldsResidualBckFuncError_SubPiZero[k+3][iPt] << endl<< endl;

            fTotalBckYields_SubPiZero[k+3][iPt]                       = fBckYields_SubPiZero[k+3][iPt] + fMesonYieldsResidualBckFunc_SubPiZero[k+3][iPt];
            fTotalBckYieldsError_SubPiZero[k+3][iPt]                  = pow(fBckYieldsError_SubPiZero[k+3][iPt]*fBckYieldsError_SubPiZero[k+3][iPt] + fMesonYieldsResidualBckFuncError_SubPiZero[k+3][iPt]*fMesonYieldsResidualBckFuncError_SubPiZero[k+3][iPt],0.5);
            fFileDataLog << "Total Background " << nameIntRange[k+3].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fTotalBckYields_SubPiZero[k+3][iPt] << "\t +- \t" << fTotalBckYieldsError_SubPiZero[k+3][iPt] << endl<< endl;
            fFileDataLog << "Background " << nameIntRange[k+3].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fBckYields_SubPiZero[k+3][iPt] << "\t +- \t" << fBckYieldsError_SubPiZero[k+3][iPt] << endl<< endl;

            fMesonYieldsCorResidualBckFunc_SubPiZero[k+3][iPt]        = fMesonYields_SubPiZero[k+3][iPt]- fMesonYieldsResidualBckFunc_SubPiZero[k+3][iPt];
            fMesonYieldsCorResidualBckFuncError_SubPiZero[k+3][iPt]   = pow(( fMesonYieldsError_SubPiZero[k+3][iPt]*fMesonYieldsError_SubPiZero[k+3][iPt] +
                                                                    fMesonYieldsResidualBckFuncError_SubPiZero[k+3][iPt]*fMesonYieldsResidualBckFuncError_SubPiZero[k+3][iPt]),0.5);
            fMesonYieldsPerEvent_SubPiZero[k+3][iPt]                  = fMesonYieldsCorResidualBckFunc_SubPiZero[k+3][iPt]/fNEvents;
            fMesonYieldsPerEventError_SubPiZero[k+3][iPt]             = fMesonYieldsCorResidualBckFuncError_SubPiZero[k+3][iPt]/fNEvents;

            // Inv mass with pz of pi0 fixed
            fFileDataLog << "Residual Background leftover (pi0 subtracted" << nameIntRange[k+3].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fMesonYieldsResidualBckFunc_FixedPzPiZero[k+3][iPt] << "\t +- \t" << fMesonYieldsResidualBckFuncError_FixedPzPiZero[k+3][iPt] << endl<< endl;

            fTotalBckYields_FixedPzPiZero[k+3][iPt]                       = fBckYields_FixedPzPiZero[k+3][iPt] + fMesonYieldsResidualBckFunc_FixedPzPiZero[k+3][iPt];
            fTotalBckYieldsError_FixedPzPiZero[k+3][iPt]                  = pow(fBckYieldsError_FixedPzPiZero[k+3][iPt]*fBckYieldsError_FixedPzPiZero[k+3][iPt] + fMesonYieldsResidualBckFuncError_FixedPzPiZero[k+3][iPt]*fMesonYieldsResidualBckFuncError_FixedPzPiZero[k+3][iPt],0.5);
            fFileDataLog << "Total Background " << nameIntRange[k+3].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fTotalBckYields_FixedPzPiZero[k+3][iPt] << "\t +- \t" << fTotalBckYieldsError_FixedPzPiZero[k+3][iPt] << endl<< endl;
            fFileDataLog << "Background " << nameIntRange[k+3].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fBckYields_FixedPzPiZero[k+3][iPt] << "\t +- \t" << fBckYieldsError_FixedPzPiZero[k+3][iPt] << endl<< endl;

            fMesonYieldsCorResidualBckFunc_FixedPzPiZero[k+3][iPt]        = fMesonYields_FixedPzPiZero[k+3][iPt]- fMesonYieldsResidualBckFunc_FixedPzPiZero[k+3][iPt];
            fMesonYieldsCorResidualBckFuncError_FixedPzPiZero[k+3][iPt]   = pow(( fMesonYieldsError_FixedPzPiZero[k+3][iPt]*fMesonYieldsError_FixedPzPiZero[k+3][iPt] +
                                                                    fMesonYieldsResidualBckFuncError_FixedPzPiZero[k+3][iPt]*fMesonYieldsResidualBckFuncError_FixedPzPiZero[k+3][iPt]),0.5);
            fMesonYieldsPerEvent_FixedPzPiZero[k+3][iPt]                  = fMesonYieldsCorResidualBckFunc_FixedPzPiZero[k+3][iPt]/fNEvents;
            fMesonYieldsPerEventError_FixedPzPiZero[k+3][iPt]             = fMesonYieldsCorResidualBckFuncError_FixedPzPiZero[k+3][iPt]/fNEvents;
        }
    }

    //******************** Data OUTPUTFILE ***************************************************
    const char* fileNameSysErrDat = Form("%s/%s/%s_%s_SystematicErrorYieldExtraction_RAWDATA_%s.dat",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix2.Data(), cutSelection.Data());
    fstream fileSysErrDat;
    fileSysErrDat.open(fileNameSysErrDat, ios::out);
    fileSysErrDat << "Calculation of the systematic error due to the yield extraction RAWDATA " << endl;
    fileSysErrDat <<  endl;
    fileSysErrDat << "fPiPlPiMiPiZeroYields" << endl;
    fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr" << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
            fGGYields[0][iPt] << "+-" << fGGYieldsError[0][iPt] << "\t" <<
            fGGYields[1][iPt] << "+-" << fGGYieldsError[1][iPt] << "\t" <<
            fGGYields[2][iPt] << "+-" << fGGYieldsError[2][iPt] << endl;
    }

    fileSysErrDat << "fPiPlPiMiPiZeroYields (pi0 mass subtracted)" << endl;
    fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr" << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
            fGGYields_SubPiZero[0][iPt] << "+-" << fGGYieldsError_SubPiZero[0][iPt] << "\t" <<
            fGGYields_SubPiZero[1][iPt] << "+-" << fGGYieldsError_SubPiZero[1][iPt] << "\t" <<
            fGGYields_SubPiZero[2][iPt] << "+-" << fGGYieldsError_SubPiZero[2][iPt] << endl;
    }

    fileSysErrDat << "fPiPlPiMiPiZeroYields (pz of pi0 fixed)" << endl;
    fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr" << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
            fGGYields_FixedPzPiZero[0][iPt] << "+-" << fGGYieldsError_FixedPzPiZero[0][iPt] << "\t" <<
            fGGYields_FixedPzPiZero[1][iPt] << "+-" << fGGYieldsError_FixedPzPiZero[1][iPt] << "\t" <<
            fGGYields_FixedPzPiZero[2][iPt] << "+-" << fGGYieldsError_FixedPzPiZero[2][iPt] << endl;
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
    fileSysErrDat << "fTotalBckYields (pi0 mass subtracted)" << endl;
    fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
            fTotalBckYields_SubPiZero[0][iPt] << "+-" << fTotalBckYieldsError_SubPiZero[0][iPt] << "\t" <<
            fTotalBckYields_SubPiZero[1][iPt] << "+-" << fTotalBckYieldsError_SubPiZero[1][iPt] << "\t" <<
            fTotalBckYields_SubPiZero[2][iPt] << "+-" << fTotalBckYieldsError_SubPiZero[2][iPt] << "\t" <<
            fTotalBckYields_SubPiZero[3][iPt] << "+-" << fTotalBckYieldsError_SubPiZero[3][iPt]<< "\t" <<
            fTotalBckYields_SubPiZero[4][iPt] << "+-" << fTotalBckYieldsError_SubPiZero[4][iPt]<< "\t" <<
            fTotalBckYields_SubPiZero[5][iPt] << "+-" << fTotalBckYieldsError_SubPiZero[5][iPt] << endl;
    }

    fileSysErrDat <<  endl;
    fileSysErrDat << "fTotalBckYields (pz of pi0 fixed)" << endl;
    fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
            fTotalBckYields_FixedPzPiZero[0][iPt] << "+-" << fTotalBckYieldsError_FixedPzPiZero[0][iPt] << "\t" <<
            fTotalBckYields_FixedPzPiZero[1][iPt] << "+-" << fTotalBckYieldsError_FixedPzPiZero[1][iPt] << "\t" <<
            fTotalBckYields_FixedPzPiZero[2][iPt] << "+-" << fTotalBckYieldsError_FixedPzPiZero[2][iPt] << "\t" <<
            fTotalBckYields_FixedPzPiZero[3][iPt] << "+-" << fTotalBckYieldsError_FixedPzPiZero[3][iPt]<< "\t" <<
            fTotalBckYields_FixedPzPiZero[4][iPt] << "+-" << fTotalBckYieldsError_FixedPzPiZero[4][iPt]<< "\t" <<
            fTotalBckYields_FixedPzPiZero[5][iPt] << "+-" << fTotalBckYieldsError_FixedPzPiZero[5][iPt] << endl;
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

    fileSysErrDat <<  endl;
    fileSysErrDat << "fMesonYields (pi0 mass subtracted)" << endl;
    fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
            fMesonYieldsCorResidualBckFunc_SubPiZero[0][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError_SubPiZero[0][iPt] << "\t" <<
            fMesonYieldsCorResidualBckFunc_SubPiZero[1][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError_SubPiZero[1][iPt] << "\t" <<
            fMesonYieldsCorResidualBckFunc_SubPiZero[2][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError_SubPiZero[2][iPt] << "\t" <<
            fMesonYieldsCorResidualBckFunc_SubPiZero[3][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError_SubPiZero[3][iPt]<< "\t" <<
            fMesonYieldsCorResidualBckFunc_SubPiZero[4][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError_SubPiZero[4][iPt]<< "\t" <<
            fMesonYieldsCorResidualBckFunc_SubPiZero[5][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError_SubPiZero[5][iPt] << endl;
    }


    fileSysErrDat <<  endl;
    fileSysErrDat << "fMesonYields (pz of pi0 fixed) "<< endl;
    fileSysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fileSysErrDat << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" <<
            fMesonYieldsCorResidualBckFunc_FixedPzPiZero[0][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError_FixedPzPiZero[0][iPt] << "\t" <<
            fMesonYieldsCorResidualBckFunc_FixedPzPiZero[1][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError_FixedPzPiZero[1][iPt] << "\t" <<
            fMesonYieldsCorResidualBckFunc_FixedPzPiZero[2][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError_FixedPzPiZero[2][iPt] << "\t" <<
            fMesonYieldsCorResidualBckFunc_FixedPzPiZero[3][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError_FixedPzPiZero[3][iPt]<< "\t" <<
            fMesonYieldsCorResidualBckFunc_FixedPzPiZero[4][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError_FixedPzPiZero[4][iPt]<< "\t" <<
            fMesonYieldsCorResidualBckFunc_FixedPzPiZero[5][iPt] << "+-" << fMesonYieldsCorResidualBckFuncError_FixedPzPiZero[5][iPt] << endl;
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

    // Plot one for each background group
    TString nameMeson;
    TString nameCanvas;
    TString namePad;
    for(Int_t k=0;k<5;k++){
        nameMeson           = Form("%s/%s_%s_MesonWith_Group%i_Bck%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),k,fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
        nameCanvas          = "MesonWithBckCanvas";
        namePad             = "MesonWithBckPad";
        PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingBackNormInvMassPtBin[k], nameMeson, nameCanvas, namePad, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcess, fCollisionSystem,kFALSE,k);

        // Plotting for InvMass Pi0 subtracted
        nameMeson           = Form("%s/%s_%s_MesonWith_Group%i_Bck%s_SubPiZero_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),k,fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
        nameCanvas          = "MesonWithBckCanvas";
        namePad             = "MesonWithBckPad";
        PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin_SubPiZero, fHistoMappingBackNormInvMassPtBin_SubPiZero[k], nameMeson, nameCanvas, namePad, fMesonMassRange_SubPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcess, fCollisionSystem,kFALSE,k);

        // Plotting for InvMass with pz of pi0 fixed
        nameMeson           = Form("%s/%s_%s_MesonWith_Group%i_Bck%s_FixedPzPiZero_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),k,fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
        nameCanvas          = "MesonWithBckCanvas";
        namePad             = "MesonWithBckPad";
        PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin_FixedPzPiZero, fHistoMappingBackNormInvMassPtBin_FixedPzPiZero[k], nameMeson, nameCanvas, namePad, fMesonMassRange_FixedPzPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcess, fCollisionSystem,kFALSE,k);

    }
    // Draw all backgroups in one plot (individual scaling)
    nameMeson           = Form("%s/%s_%s_MesonWith_AllGroups_Bck%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    nameCanvas          = "MesonWithBckCanvas";
    namePad             = "MesonWithBckPad";
    PlotInvMassInPtBinsBckGroups(fHistoMappingGGInvMassPtBin, fHistoMappingBackNormInvMassPtBin, nameMeson, nameCanvas, namePad, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcess, fCollisionSystem);

    // Draw all backgroups in one plot (individual scaling and pi0 mass subtracted)
    nameMeson           = Form("%s/%s_%s_MesonWith_AllGroups_Bck%s_SubPiZero_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    nameCanvas          = "MesonWithBckCanvas";
    namePad             = "MesonWithBckPad";
    PlotInvMassInPtBinsBckGroups(fHistoMappingGGInvMassPtBin_SubPiZero, fHistoMappingBackNormInvMassPtBin_SubPiZero, nameMeson, nameCanvas, namePad, fMesonMassRange_SubPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcess, fCollisionSystem);

    // Draw all backgroups in one plot (individual scaling and pz of pi0 fixed)
    nameMeson           = Form("%s/%s_%s_MesonWith_AllGroups_Bck%s_FixedPzPiZero_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    nameCanvas          = "MesonWithBckCanvas";
    namePad             = "MesonWithBckPad";
    PlotInvMassInPtBinsBckGroups(fHistoMappingGGInvMassPtBin_FixedPzPiZero, fHistoMappingBackNormInvMassPtBin_FixedPzPiZero, nameMeson, nameCanvas, namePad, fMesonMassRange_FixedPzPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcess, fCollisionSystem);

    // Meson Subtracted
    TString nameMesonSub        = Form("%s/%s_%s_MesonSubtracted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
    TString nameCanvasSub       = "MesonCanvasSubtracted";
    TString namePadSub          = "MesonPadSubtracted";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");

    // Meson Subtracted minus InvMass pi0
    nameMesonSub        = Form("%s/%s_%s_MesonSubtracted_SubPiZero_%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
    nameCanvasSub       = "MesonCanvasSubtracted";
    namePadSub          = "MesonPadSubtracted";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassPtBin_SubPiZero, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassPtBin_SubPiZero, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange_SubPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");

    // Meson Subtracted with pz of pi0 fixed
    nameMesonSub        = Form("%s/%s_%s_MesonSubtracted_FixedPzPiZero_%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
    nameCanvasSub       = "MesonCanvasSubtracted";
    namePadSub          = "MesonPadSubtracted";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassPtBin_FixedPzPiZero, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassPtBin_FixedPzPiZero, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange_FixedPzPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, kFALSE,fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");

    // Meson Subtracted backfit
    nameMesonSub                = Form("%s/%s_%s_MesonSubtractedBackFit%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasSubtractedBackFit";
    namePadSub                  = "MesonPadSubtractedBackFit";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassBackFitPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassBackFitPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");

    // Meson Subtracted backfit (minus pi0 inv mass) (NOTE: isMC is always disabled because no true histos exist for minus pi0 inv mass)
    nameMesonSub                = Form("%s/%s_%s_MesonSubtractedBackFit_SubPiZero_%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasSubtractedBackFit";
    namePadSub                  = "MesonPadSubtractedBackFit";
    cout << " fHistoMappingGGInvMassBachFitPtBin_SubPiZero" << fHistoMappingGGInvMassBackFitPtBin_SubPiZero << endl;
    cout << " fFitSignalInvMassBackFitPtBin_SubPiZero" << fFitSignalInvMassBackFitPtBin_SubPiZero << endl;
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassBackFitPtBin_SubPiZero, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassBackFitPtBin_SubPiZero, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange_SubPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, kFALSE,fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");

    // Meson Subtracted backfit (pz of pi0 is fixed) (NOTE: isMC is always disabled because no true histos exist)
    nameMesonSub                = Form("%s/%s_%s_MesonSubtractedBackFit_FixedPzPiZero_%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasSubtractedBackFit";
    namePadSub                  = "MesonPadSubtractedBackFit";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassBackFitPtBin_FixedPzPiZero, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassBackFitPtBin_FixedPzPiZero, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange_FixedPzPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, kFALSE,fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");

    //loop through Pt slices and plot inv mass around peak (only in DATA, not MC)

    // Meson Back Fit
    nameMesonSub                = Form("%s/%s_%s_MesonBackFit%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasBackFit";
    namePadSub                  = "MesonPadBackFit";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fBackgroundFitPol, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");

   // Meson Back Fit (minus pi0 inv mass)
    nameMesonSub                = Form("%s/%s_%s_MesonBackFit_SubPiZero_%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasBackFit";
    namePadSub                  = "MesonPadBackFit";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassPtBin_SubPiZero, fHistoMappingTrueMesonInvMassPtBins, fBackgroundFitPol_SubPiZero, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange_SubPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, kFALSE, fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");

    // Meson Back Fit (pz of pi0 fixed)
    nameMesonSub                = Form("%s/%s_%s_MesonBackFit_FixedPzPiZero_%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasBackFit";
    namePadSub                  = "MesonPadBackFit";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassPtBin_FixedPzPiZero, fHistoMappingTrueMesonInvMassPtBins, fBackgroundFitPol_FixedPzPiZero, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange_FixedPzPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, kFALSE, fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");

    // Meson with back left
    nameMeson                   = Form("%s/%s_%s_MesonWithBckLeft%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    nameCanvas                  = "MesonWithBckCanvasLeft";
    namePad                     = "MesonWithBckPadLeft";
    PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingBackNormInvMassLeftPtBin, nameMeson, nameCanvas, namePad,  fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcess, fCollisionSystem);

    // Meson with back left (minus pi0 inv mass)
    nameMeson                   = Form("%s/%s_%s_MesonWithBckLeft_SubPiZero_%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    nameCanvas                  = "MesonWithBckCanvasLeft";
    namePad                     = "MesonWithBckPadLeft";
    PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin_SubPiZero, fHistoMappingBackNormInvMassLeftPtBin_SubPiZero, nameMeson, nameCanvas, namePad,  fMesonMassRange_SubPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, kFALSE ,fDecayChannel, fDetectionProcess, fCollisionSystem);

    // Meson with back left (pz of pi0 fixed)
    nameMeson                   = Form("%s/%s_%s_MesonWithBckLeft_FixedPzPiZero_%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    nameCanvas                  = "MesonWithBckCanvasLeft";
    namePad                     = "MesonWithBckPadLeft";
    PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin_FixedPzPiZero, fHistoMappingBackNormInvMassLeftPtBin_FixedPzPiZero, nameMeson, nameCanvas, namePad,  fMesonMassRange_FixedPzPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, kFALSE ,fDecayChannel, fDetectionProcess, fCollisionSystem);

    // Meson with background and fit
    nameMesonSub                = Form("%s/%s_%s_MesonWithBGAndFit%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasWithBGAndFit";
    namePadSub                  = "MesonPadWithBGAndFit";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins,fFitWithPol2ForBG, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC , fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");

    // Meson with background and fit (minus pi0 mass)
    nameMesonSub                = Form("%s/%s_%s_MesonWithBGAndFit_SubPiZero%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasWithBGAndFit";
    namePadSub                  = "MesonPadWithBGAndFit";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassPtBin_SubPiZero, fHistoMappingTrueMesonInvMassPtBins,fFitWithPol2ForBG_SubPiZero, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange_SubPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, kFALSE , fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");

    // Meson with background and fit (pz of pi0 fixed)
    nameMesonSub                = Form("%s/%s_%s_MesonWithBGAndFit_FixedPzPiZero%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasWithBGAndFit";
    namePadSub                  = "MesonPadWithBGAndFit";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassPtBin_FixedPzPiZero, fHistoMappingTrueMesonInvMassPtBins,fFitWithPol2ForBG_FixedPzPiZero, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange_FixedPzPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, kFALSE , fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");

    // Meson Subtracted left side
    nameMesonSub                = Form("%s/%s_%s_MesonSubtractedLeft%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasSubtractedLeft";
    namePadSub                  = "MesonPadSubtractedLeft";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassLeftPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitInvMassLeftPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");

    // Meson Subtracted left side (pi0 inv mass subtracted)
    nameMesonSub                = Form("%s/%s_%s_MesonSubtractedLeft_SubPiZero_%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasSubtractedLeft";
    namePadSub                  = "MesonPadSubtractedLeft";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassLeftPtBin_SubPiZero, fHistoMappingTrueMesonInvMassPtBins, fFitInvMassLeftPtBin_SubPiZero, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange_SubPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, kFALSE, fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");

    // Meson Subtracted left side (pz of pi0 fixed)
    nameMesonSub                = Form("%s/%s_%s_MesonSubtractedLeft_FixedPzPiZero_%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasSubtractedLeft";
    namePadSub                  = "MesonPadSubtractedLeft";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassLeftPtBin_FixedPzPiZero, fHistoMappingTrueMesonInvMassPtBins, fFitInvMassLeftPtBin_FixedPzPiZero, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange_FixedPzPiZero, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, kFALSE, fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");

    PlotExampleInvMassBinsV2(fHistoMappingGGInvMassPtBin[fExampleBin], fHistoMappingSignalInvMassPtBin[fExampleBin], fHistoMappingBackNormInvMassPtBin[0][fExampleBin],
                        fFitSignalInvMassPtBin[fExampleBin], fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                        fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, 0, fScaleFac, fMode, addSig );
    // In example bin with pi0 inv mass subtracted

    PlotExampleInvMassBinsV2(fHistoMappingGGInvMassPtBin_SubPiZero[fExampleBin], fHistoMappingSignalInvMassPtBin_SubPiZero[fExampleBin], fHistoMappingBackNormInvMassPtBin_SubPiZero[0][fExampleBin],
                        fFitSignalInvMassPtBin_SubPiZero[fExampleBin], fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassRange_SubPiZero, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                        fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, 0, fScaleFac, fMode, addSig );
    // In example bin with pz of pi0 fixed

    PlotExampleInvMassBinsV2(fHistoMappingGGInvMassPtBin_FixedPzPiZero[fExampleBin], fHistoMappingSignalInvMassPtBin_FixedPzPiZero[fExampleBin], fHistoMappingBackNormInvMassPtBin_FixedPzPiZero[0][fExampleBin],
                        fFitSignalInvMassPtBin_FixedPzPiZero[fExampleBin], fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMass_FixedPzPiZero, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                        fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, 0, fScaleFac, fMode, addSig );

    PlotExampleInvMassBinsBckFit(fHistoMappingGGInvMassPtBin[fExampleBin], fHistoMappingGGInvMassBackFitPtBin[fExampleBin], fBackgroundFitPol[fExampleBin],fHistoBckFitConfidence[fExampleBin],
                         fFitSignalInvMassBackFitPtBin[fExampleBin], fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                        fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, 0, fScaleFac, fMode, addSig );
    // In example bin with pi0 inv mass subtracted
    PlotExampleInvMassBinsBckFit(fHistoMappingGGInvMassPtBin_SubPiZero[fExampleBin], fHistoMappingGGInvMassBackFitPtBin_SubPiZero[fExampleBin], fBackgroundFitPol_SubPiZero[fExampleBin],fHistoBckFitConfidence_SubPiZero[fExampleBin],
                         fFitSignalInvMassBackFitPtBin_SubPiZero[fExampleBin], fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassRange_SubPiZero, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                        fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, 0, fScaleFac, fMode, addSig );

    // In example bin with pz of pi0 fixed
    PlotExampleInvMassBinsBckFit(fHistoMappingGGInvMassPtBin_FixedPzPiZero[fExampleBin], fHistoMappingGGInvMassBackFitPtBin_FixedPzPiZero[fExampleBin], fBackgroundFitPol_FixedPzPiZero[fExampleBin],fHistoBckFitConfidence_FixedPzPiZero[fExampleBin],
                         fFitSignalInvMassBackFitPtBin_FixedPzPiZero[fExampleBin], fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMass_FixedPzPiZero, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                        fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, 0, fScaleFac, fMode, addSig );

    // Background groups individually scaled in example bin
    PlotExampleInvMassBinsBckGroups(fHistoMappingGGInvMassPtBin[fExampleBin], fHistoMappingBackNormInvMassPtBin,fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                                    fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, 0, fScaleFac, fMode, addSig ,"InvMassBinBckGroups");

    // Background groups individually scaled in example bin (pi0 subtracted)
    PlotExampleInvMassBinsBckGroups(fHistoMappingGGInvMassPtBin_SubPiZero[fExampleBin], fHistoMappingBackNormInvMassPtBin_SubPiZero,fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassRange_SubPiZero, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                                    fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, 0, fScaleFac, fMode, addSig ,"InvMassBinBckGroups_SubPiZero");

    // Background groups individually scaled in example bin (pz of pi0 fixed)
    PlotExampleInvMassBinsBckGroups(fHistoMappingGGInvMassPtBin_FixedPzPiZero[fExampleBin], fHistoMappingBackNormInvMassPtBin_FixedPzPiZero,fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassRange_FixedPzPiZero, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                                    fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, 0, fScaleFac, fMode, addSig ,"InvMassBinBckGroups_FixedPzPiZero");

    // Background groups scaled with same factor in example bin
    PlotExampleInvMassBinsBckGroups(fHistoMappingGGInvMassPtBin[fExampleBin], fHistoMappingBackSameNormInvMassPtBin,fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                                    fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, 0, fScaleFac, fMode, addSig ,"InvMassBinBckGroupsSameScale");
    // Background groups scaled with same factor in example bin (pi0 subtracted)
    PlotExampleInvMassBinsBckGroups(fHistoMappingGGInvMassPtBin_SubPiZero[fExampleBin], fHistoMappingBackSameNormInvMassPtBin_SubPiZero,fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassRange_SubPiZero, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                                    fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, 0, fScaleFac, fMode, addSig ,"InvMassBinBckGroupsSameScale_SubPiZero");

    // Background groups scaled with same factor in example bin (pz of pi0 fixed)
    PlotExampleInvMassBinsBckGroups(fHistoMappingGGInvMassPtBin_FixedPzPiZero[fExampleBin], fHistoMappingBackSameNormInvMassPtBin_FixedPzPiZero,fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassRange_FixedPzPiZero, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                                    fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, 0, fScaleFac, fMode, addSig ,"InvMassBinBckGroupsSameScale_FixedPzPiZero");

    if(fIsMC){
        TString nameMesonTrue   = Form("%s/%s_%s_TrueMesonFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
        TString nameCanvasTrue  = "TrueMesonCanvasFitted";
        TString namePadTrue     = "TrueMesonPadFitted";
        PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins, fHistoMappingTrueMesonInvMassPtBins, fFitTrueSignalInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem,"MC validated signal",kTRUE,"Fit","mixed evt. subtr. #it{M}_{#pi^{+}#pi^{-}#pi^{0}}");
    }

    CreatePtHistos();
    FillPtHistos();

    if(fIsMC){

        FillHistosArrayMC(fHistoMCMesonPtWithinAcceptance,fHistoMCMesonPt,fDeltaPt,Form("%s/%s_%s_MCYieldMesonFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data()));

        CalculateMesonAcceptance();
        //       cout << "Calculated MesonAcceptance" << endl;
        for (Int_t k = 0; k < 6; k++){ //TODO is 6 correct ?
            fNameHistoEffi                      = Form("Meson%sEffiPt",nameIntRange[k].Data());
            cout << fNameHistoEffi.Data() << endl;
            fHistoMonteMesonEffiPt[k]           = CalculateMesonEfficiency(fHistoYieldMeson[k], fHistoMCMesonWithinAccepPt, fNameHistoEffi );
        }
        fNameHistoEffi="MesonEffiBackFitPt";
        fHistoMonteMesonEffiBackFitPt           = CalculateMesonEfficiency(fHistoYieldMesonBackFit, fHistoMCMesonWithinAccepPt, fNameHistoEffi );

        // True Meson (only once case, because no normalization
        for (Int_t k = 0; k < 3; k++){
            fNameHistoEffi                          = Form("TrueMeson%sEffiPt",nameIntRange[k].Data());
            fHistoMCTrueMesonEffiPt[k]              = CalculateMesonEfficiency(fHistoYieldTrueMeson[k], fHistoMCMesonWithinAccepPt, fNameHistoEffi);
            cout << fNameHistoEffi.Data() << endl;

            // True meson efficiencies with possibly fully weighted inputs taking the average weight per inv mass bin in the original binning of the TrueMesonInvMass vs pT plot
            // should give on average the same as TrueMesonEffiPt
            fNameHistoEffi                          = Form("TrueMeson%sEffiPtReweighted",nameIntRange[k].Data());
            cout << fNameHistoEffi.Data() << endl;
            fHistoMCTrueMesonEffiPtReweighted[k]    = CalculateMesonEfficiency(fHistoYieldTrueMesonReweighted[k], fHistoMCMesonWithinAccepPt, fNameHistoEffi);
        }

        SaveCorrectionHistos(fCutSelection, fPrefix2);
    }
    SaveHistos(fIsMC, fCutSelection, fPrefix2);
    fFileErrLog.close();
    fFileDataLog.close();
    Delete();
}

void ProcessBckFitSubtraction(TH1D *fGammaGamma, Int_t i, Double_t * fPeakRangeDummy, Double_t *fFitRangeDummy, TString energy, TString suffix, TString cutSelection, TString meson,Int_t InvMassType){

    fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i] = (TH1D*)fGammaGamma->Clone(Form("GG_WithoutSigal_%i",i));
    fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->Sumw2();

    for (Int_t binx= 0; binx < fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetNbinsX()+1; binx++){
        if(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinContent(binx) == 0){
            fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinError(binx,1.);
            fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinContent(binx,0.);
        }
        if(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) > fPeakRangeDummy[0] && fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) < fPeakRangeDummy[1]){
            fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinContent(binx,0.0);
            fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinError(binx,0.0);
        }
    }

    TF1 *FitBackFunc;//, *FitFuncCenter, *FitFuncHigh;
    FitBackFunc = new TF1("BGfit",FunctionBGExclusion,fFitRangeDummy[0],fFitRangeDummy[1],5);
    TFitResultPtr resultBckFit =fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->Fit(FitBackFunc,"QMRES0","",fFitRangeDummy[0],fFitRangeDummy[1]);//QMRES0
    Double_t FitParams[5];

    if(InvMassType==0){
      fBackgroundFitPol[i] = NULL;
      fBackgroundFitPol[i] = new TF1("BGfit","pol4",fFitRangeDummy[0],fFitRangeDummy[1]);

      FitBackFunc->GetParameters(&FitParams[0]);
      fBackgroundFitPol[i]->SetParameters(FitParams);

      fBackgroundFitPol[i]->SetRange(fMesonFitRange[0],fMesonFitRange[1]);

      for (Int_t binx= 0; binx < fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetNbinsX()+1; binx++){
          if(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) > fFitRangeDummy[0] && fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) < fFitRangeDummy[1]){

              Double_t area = fBackgroundFitPol[i]->Integral(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinLowEdge(binx),
                                                             (fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinLowEdge(binx))+fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinWidth(binx));
              Double_t area_err = fBackgroundFitPol[i]->IntegralError(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinLowEdge(binx),
                                                                     (fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinLowEdge(binx))+fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinWidth(binx),
                                                                      resultBckFit->GetParams(), resultBckFit->GetCovarianceMatrix().GetMatrixArray() );
              fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinContent(binx,area/(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinWidth(binx)));
              fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinError(binx,area_err/(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinWidth(binx)));

          }
          // Clone histo for error band
          fHistoBckFitConfidence[i] =  (TH1D*)fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->Clone(Form(" fHistoBckFitConfidence_%i",i));
      }

      fHistoMappingGGInvMassBackFitPtBin[i] = (TH1D*) fGammaGamma->Clone(Form("GG_SubtractedSignal_%i",i));
      fHistoMappingGGInvMassBackFitPtBin[i]->Sumw2();
      fHistoMappingGGInvMassBackFitPtBin[i]->Add(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i],-1);
    }else if(InvMassType==1){
      fBackgroundFitPol_SubPiZero[i] = NULL;
      fBackgroundFitPol_SubPiZero[i] = new TF1("BGfit","pol4",fFitRangeDummy[0],fFitRangeDummy[1]);

      FitBackFunc->GetParameters(&FitParams[0]);
      fBackgroundFitPol_SubPiZero[i]->SetParameters(FitParams);

      fBackgroundFitPol_SubPiZero[i]->SetRange(fMesonFitRange_SubPiZero[0],fMesonFitRange_SubPiZero[1]);

      for (Int_t binx= 0; binx < fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetNbinsX()+1; binx++){
          if(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) > fFitRangeDummy[0] && fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) < fFitRangeDummy[1]){

              Double_t area = fBackgroundFitPol_SubPiZero[i]->Integral(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinLowEdge(binx),
                                                             (fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinLowEdge(binx))+fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinWidth(binx));
              Double_t area_err = fBackgroundFitPol_SubPiZero[i]->IntegralError(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinLowEdge(binx),
                                                                     (fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinLowEdge(binx))+fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinWidth(binx),
                                                                      resultBckFit->GetParams(), resultBckFit->GetCovarianceMatrix().GetMatrixArray() );
              fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinContent(binx,area/(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinWidth(binx)));
              fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinError(binx,area_err/(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinWidth(binx)));

          }
          // Clone histo for error band
          fHistoBckFitConfidence_SubPiZero[i] =  (TH1D*)fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->Clone(Form(" fHistoBckFitConfidence_SubPiZero_%i",i));
      }
      fHistoMappingGGInvMassBackFitPtBin_SubPiZero[i] = (TH1D*) fGammaGamma->Clone(Form("GG_SubtractedSignal_SubPiZero_%i",i));
      fHistoMappingGGInvMassBackFitPtBin_SubPiZero[i]->Sumw2();
      fHistoMappingGGInvMassBackFitPtBin_SubPiZero[i]->Add(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i],-1);
    }else if(InvMassType==2){
      fBackgroundFitPol_FixedPzPiZero[i] = NULL;
      fBackgroundFitPol_FixedPzPiZero[i] = new TF1("BGfit","pol4",fFitRangeDummy[0],fFitRangeDummy[1]);

      FitBackFunc->GetParameters(&FitParams[0]);
      fBackgroundFitPol_FixedPzPiZero[i]->SetParameters(FitParams);

      fBackgroundFitPol_FixedPzPiZero[i]->SetRange(fMesonFitRange_FixedPzPiZero[0],fMesonFitRange_FixedPzPiZero[1]);

      for (Int_t binx= 0; binx < fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetNbinsX()+1; binx++){
          if(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) > fFitRangeDummy[0] && fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) < fFitRangeDummy[1]){

              Double_t area = fBackgroundFitPol_FixedPzPiZero[i]->Integral(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinLowEdge(binx),
                                                             (fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinLowEdge(binx))+fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinWidth(binx));
              Double_t area_err = fBackgroundFitPol_FixedPzPiZero[i]->IntegralError(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinLowEdge(binx),
                                                                     (fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinLowEdge(binx))+fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinWidth(binx),
                                                                      resultBckFit->GetParams(), resultBckFit->GetCovarianceMatrix().GetMatrixArray() );
              fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinContent(binx,area/(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinWidth(binx)));
              fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinError(binx,area_err/(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinWidth(binx)));

          }
          // Clone histo for error band
          fHistoBckFitConfidence_FixedPzPiZero[i] =  (TH1D*)fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->Clone(Form(" fHistoBckFitConfidence_FixedPzPiZero_%i",i));
      }
      fHistoMappingGGInvMassBackFitPtBin_FixedPzPiZero[i] = (TH1D*) fGammaGamma->Clone(Form("GG_SubtractedSignal_FixedPzPiZero_%i",i));
      fHistoMappingGGInvMassBackFitPtBin_FixedPzPiZero[i]->Sumw2();
      fHistoMappingGGInvMassBackFitPtBin_FixedPzPiZero[i]->Add(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i],-1);
    }

}

void ProcessEM( TH1D* fGammaGamma,
                TH1D* fBck,
                Double_t * fBGFitRangeEM
              ){
    //cout<<"{"<<endl<<"ProcessEM started"<<endl;
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

    Double_t 	r= fGammaGamma->Integral(fGammaGamma->GetXaxis()->FindBin(fBGFitRangeEM[0]),fGammaGamma->GetXaxis()->FindBin(fBGFitRangeEM[1]));
    Double_t 	b= fBck->Integral(fBck->GetXaxis()->FindBin(fBGFitRangeEM[0]),fBck->GetXaxis()->FindBin(fBGFitRangeEM[1]));

    if(b != 0) fNorm = r/b;

    fBckNorm->Scale(fNorm);
    //    cout<<"r="<<r<<" b="<<b<<" r/b="<<r/b<< " " << endl;

    Int_t numberOfZeros = 0;
    for (Int_t i = 1; i < fBckNorm->GetNbinsX()+1; i++){
        if (fBckNorm->GetBinContent(i) == 0){
            numberOfZeros++;
            if (fNorm > 1.){
                fBckNorm->SetBinError(i,1.);
                fBckNorm->SetBinContent(i,0.);
            }
        }
    }
    //    cout << "counted " << numberOfZeros << " in the normalized BG : "<< (Double_t)numberOfZeros/fBck->GetNbinsX() << endl;

    fSignal = (TH1D*)fGammaGamma->Clone("fSignal");
    fSignal->Sumw2();
    if ((Double_t)numberOfZeros/fBck->GetNbinsX()< 0.25) fSignal->Add(fBckNorm,-1.);
    //cout<<"ProcessEM finished"<<endl<<"}"<<endl;
}


void ProcessEMLeftRight(    TH1D* fFgr,
                            TH1D* fBck,
                            Double_t* fBGFitRangeEMLeft,
                            Double_t* fBGFitRangeEMRight){
    //cout<<"{"<<endl<<"ProcessEM started"<<endl;
    for (Int_t binx= 0; binx < fFgr->GetNbinsX()+1; binx++){
        if(fFgr->GetBinContent(binx) == 0){
            fFgr->SetBinError(binx,1.);
            fFgr->SetBinContent(binx,0.);
        }
    }

    fBckNorm = (TH1D*)fBck->Clone("fBckNorm");
    fFgr->Sumw2();
    fBck->Sumw2();
    fBckNorm->Sumw2();

    Double_t 	norm = 1;

    Double_t 	r= fFgr->Integral(fFgr->GetXaxis()->FindBin(fBGFitRangeEMLeft[0]),fFgr->GetXaxis()->FindBin(fBGFitRangeEMLeft[1]));
    Double_t 	b= fBck->Integral(fBck->GetXaxis()->FindBin(fBGFitRangeEMLeft[0]),fBck->GetXaxis()->FindBin(fBGFitRangeEMLeft[1]));

    r+= fFgr->Integral(fFgr->GetXaxis()->FindBin(fBGFitRangeEMRight[0]),fFgr->GetXaxis()->FindBin(fBGFitRangeEMRight[1]));
    b+= fBck->Integral(fBck->GetXaxis()->FindBin(fBGFitRangeEMRight[0]),fBck->GetXaxis()->FindBin(fBGFitRangeEMRight[1]));

    if(b != 0) norm = r/b;

    fBckNorm->Scale(norm);
    //    cout<<"r="<<r<<" b="<<b<<" r/b="<<r/b<< " " << endl;

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
    //    cout << "counted " << numberOfZeros << " in the normalized BG : "<< (Double_t)numberOfZeros/fBck->GetNbinsX() << endl;

    fSignal = (TH1D*)fFgr->Clone("fSignal");
    fSignal->Sumw2();

    if ((Double_t)numberOfZeros/fBck->GetNbinsX()< 0.25) fSignal->Add(fBckNorm,-1.);
    //cout<<"ProcessEM finished"<<endl<<"}"<<endl;
                            }


void ProcessRatioSignalBackground(TH1D* fGammaGamma, TH1D* fBck)
{
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
/*
void ProduceBckProperWeighting(TList* fBackgroundContainer,TList* fMotherContainer){

    THnSparseF* fSparseMotherZM;
    THnSparseF* fSparseBckZM;
    fSparseMotherZM = (THnSparseF*)fMotherContainer->FindObject("Back_Mother_InvMass_Pt_z_m");
    fSparseBckZM = (THnSparseF*)fBackgroundContainer->FindObject("Back_Back_InvMass_Pt_z_m");



    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fHistoWeightsBGZbinVsMbin[iPt] = new  TH2F("BGWeights", "", fSparseMotherZM->GetAxis(2)->GetNbins(),  0, fSparseMotherZM->GetAxis(2)->GetNbins(),fSparseMotherZM->GetAxis(3)->GetNbins(),  0, fSparseMotherZM->GetAxis(3)->GetNbins());
        fHistoWeightsBGZbinVsMbin[iPt]->GetYaxis()->SetTitle("M-bins");
        fHistoWeightsBGZbinVsMbin[iPt]->GetXaxis()->SetTitle("Z-bins");
        fHistoWeightsBGZbinVsMbin[iPt]->Sumw2();
        fHistoFillPerEventBGZbinVsMbin[iPt] = new  TH2F("BGPoolsFillstatus", "", fSparseMotherZM->GetAxis(2)->GetNbins(),  0, fSparseMotherZM->GetAxis(2)->GetNbins(),fSparseMotherZM->GetAxis(3)->GetNbins(),  0, fSparseMotherZM->GetAxis(3)->GetNbins());
        fHistoFillPerEventBGZbinVsMbin[iPt]->GetYaxis()->SetTitle("M-bins");
        fHistoFillPerEventBGZbinVsMbin[iPt]->GetXaxis()->SetTitle("Z-bins");
        fHistoFillPerEventBGZbinVsMbin[iPt]->Sumw2();

        for (Int_t z=0;z < fSparseMotherZM->GetAxis(2)->GetNbins()-1;z++){
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
    //                cout << z << "\t" << m << "\t"  << fScalingFactorBck[z][m] << "\t" << 20./fBackgroundMultNumber << endl;
                if ( fScalingFactorBck[z][m]> (20./fBackgroundMultNumber) ){
    //                   cout << "fail safe entered" << endl;
                    fScalingFactorBck[z][m]=1./fBackgroundMultNumber;
    //                   cout << "\t" <<  z << "\t" << m << "\t"  << fScalingFactorBck[z][m] << "\t" << 20./fBackgroundMultNumber << endl;
                }
                }
                fHistoMappingBackInvMassPtBin[iPt]->Add(fHistoBckZMProj,fScalingFactorBck[z][m]);
                //cout<<"ADDING 7 (iPt factor [z][m])"<<endl;
                fHistoWeightsBGZbinVsMbin[iPt]->Fill(z+0.5,m+0.5,fScalingFactorBck[z][m]);
                fHistoFillPerEventBGZbinVsMbin[iPt]->Fill(z+0.5,m+0.5,fHistoBckZMProj->GetEntries());
                fHistoMotherZMProj->Clear();
                fHistoBckZMProj->Clear();

            }
        }
        fHistoMappingBackInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
        //fHistoMappingBackInvMassPtBin[iPt]->Scale(1./fNRebin[iPt]);
        for (Int_t ii = 0; ii < fHistoMappingBackInvMassPtBin[iPt]->GetNbinsX()+1; ii++){
            if(fHistoMappingBackInvMassPtBin[iPt]->GetBinContent(ii) == 0){
                fHistoMappingBackInvMassPtBin[iPt]->SetBinContent(ii,0.);
                fHistoMappingBackInvMassPtBin[iPt]->SetBinError(ii,1.);
            }
        }
        fFileDataLog << "Scaling Background factors for Pt bin " << iPt << " z m " << endl;
        // 		cout << "Scaling Background factors for Pt bin " << iPt << " z m " << endl;
        for (Int_t z=0; z < fSparseMotherZM->GetAxis(2)->GetNbins()-1; z++){
            fFileDataLog << fScalingFactorBck[z][0] << "\t" << fScalingFactorBck[z][1] << "\t" << fScalingFactorBck[z][2] << "\t" << fScalingFactorBck[z][3] << endl;
            // 			cout << fScalingFactorBck[z][0] << "\t" << fScalingFactorBck[z][1] << "\t" << fScalingFactorBck[z][2] << "\t" << fScalingFactorBck[z][3] << endl;
        }
    }

    for (Int_t z=0;z < fSparseMotherZM->GetAxis(2)->GetNbins()-1;z++){
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
            //cout<<"ADDING 8"<<endl;

            fHistoMotherZMProj->Clear();
            fHistoBckZMProj->Clear();
        }
    }
    fMesonFullPtBackground->Rebin(fNRebin[4]);
    //fMesonFullPtBackground->Scale(1./fNRebin[4]);


    for (Int_t z=0;z < fSparseMotherZM->GetAxis(2)->GetNbins()-1;z++){
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
            //cout<<"ADDING 9"<<endl;

            fHistoMotherZMProj->Clear();
            fHistoBckZMProj->Clear();
        }
    }
    fFittingHistMidPtBackground->Rebin(fNRebin[4]);
    //fFittingHistMidPtBackground->Scale(1./fNRebin[4]);
}
*/
void ProduceBckWithoutWeighting(TH2D **fBckInvMassVSPtDummy,TH2D **fBckInvMassVSPtDummy_SubPiZero, TH2D **fBckInvMassVSPtDummy_FixedPzPiZero){
    //calculation background for midPt without weighting
    for(Int_t k = 0; k<5; k++){
        if(fBckInvMassVSPtDummy[k]==NULL){
            continue;
        }
        //calculation background for midPt without weighting
        Int_t startBinMidPt = fBckInvMassVSPtDummy[k]->GetYaxis()->FindBin(fMidPt[0]+0.001);
        Int_t endBinMidPt = fBckInvMassVSPtDummy[k]->GetYaxis()->FindBin(fMidPt[1]-0.001);
        fFittingHistMidPtBackground[k] = new TH1D(Form("Mapping_Back%i_InvMass_MidPt",k),Form("Mapping_Back%i_InvMass_MidPt",k),fBckInvMassVSPtDummy[k]->GetNbinsX(),0.,1.);
        fFittingHistMidPtBackground[k]->Sumw2();
        fBckInvMassVSPtDummy[k]->ProjectionX(Form("Mapping_Back%i_InvMass_MidPt",k),startBinMidPt,endBinMidPt);
        fFittingHistMidPtBackground[k]=(TH1D*)gDirectory->Get(Form("Mapping_Back%i_InvMass_MidPt",k));
        fFittingHistMidPtBackground[k]->Rebin(fNRebin[4]);

        //calulation background for fullPt without weighting
        fMesonFullPtBackground[k] = new TH1D(Form("Mapping_Back%i_InvMass_FullPt",k),Form("Mapping_Back%i_InvMass_FullPt",k),fBckInvMassVSPtDummy[k]->GetNbinsX(),0.,1.);
        fMesonFullPtBackground[k]->Sumw2();
        Int_t startBinFullPt = fBckInvMassVSPtDummy[k]->GetYaxis()->FindBin(fFullPt[0]+0.001);
        Int_t endBinFullPt = fBckInvMassVSPtDummy[k]->GetYaxis()->FindBin(fFullPt[1]-0.001);
        fBckInvMassVSPtDummy[k]->ProjectionX(Form("Mapping_Back%i_InvMass_FullPt",k),startBinFullPt,endBinFullPt);
        fMesonFullPtBackground[k]=(TH1D*)gDirectory->Get(Form("Mapping_Back%i_InvMass_FullPt",k));
        fMesonFullPtBackground[k]->Rebin(fNRebin[4]);

        for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
            fNameHistoBack = Form("Mapping_Back%i_InvMass_in_Pt_Bin%02d",k ,iPt);
            if(fHistoMappingBackInvMassPtBin[k][iPt]!= NULL){
                delete fHistoMappingBackInvMassPtBin[k][iPt];
                fHistoMappingBackInvMassPtBin[k][iPt]=NULL;
            }
            fHistoMappingBackInvMassPtBin[k][iPt]=new TH1D(fNameHistoBack.Data(),fNameHistoBack.Data(),fBckInvMassVSPtDummy[k]->GetNbinsX(),0.,1.);
            fHistoMappingBackInvMassPtBin[k][iPt]->Sumw2();
            Int_t startBin = fBckInvMassVSPtDummy[k]->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
            Int_t endBin = fBckInvMassVSPtDummy[k]->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

            fBckInvMassVSPtDummy[k]->ProjectionX(fNameHistoBack.Data(),startBin,endBin);
            fHistoMappingBackInvMassPtBin[k][iPt]=(TH1D*)gDirectory->Get(fNameHistoBack.Data());
            if(fNRebin[iPt]>1){
                fHistoMappingBackInvMassPtBin[k][iPt]->Rebin(fNRebin[iPt]);
            }

        }

        // Begin Calculation for histos with PiZero mass subtracted
        if(fBckInvMassVSPtDummy_SubPiZero[k]==NULL){
            continue;
        }
        //calculation background for midPt without weighting
        startBinMidPt = fBckInvMassVSPtDummy_SubPiZero[k]->GetYaxis()->FindBin(fMidPt[0]+0.001);
        endBinMidPt = fBckInvMassVSPtDummy_SubPiZero[k]->GetYaxis()->FindBin(fMidPt[1]-0.001);
        fFittingHistMidPtBackground_SubPiZero[k] = new TH1D(Form("Mapping_Back%i_InvMass_SubPiZero_MidPt",k),Form("Mapping_Back%i_InvMass_SubPiZero_MidPt",k),fBckInvMassVSPtDummy_SubPiZero[k]->GetNbinsX(),0.,1.);
        fFittingHistMidPtBackground_SubPiZero[k]->Sumw2();
        fBckInvMassVSPtDummy_SubPiZero[k]->ProjectionX(Form("Mapping_Back%i_InvMass_SubPiZero_MidPt",k),startBinMidPt,endBinMidPt);
        fFittingHistMidPtBackground_SubPiZero[k]=(TH1D*)gDirectory->Get(Form("Mapping_Back%i_InvMass_SubPiZero_MidPt",k));
        fFittingHistMidPtBackground_SubPiZero[k]->Rebin(fNRebin[4]);

        //calulation background for fullPt without weighting
        fMesonFullPtBackground_SubPiZero[k] = new TH1D(Form("Mapping_Back%i_InvMass_SubPiZero_FullPt",k),Form("Mapping_Back%i_InvMass_SubPiZero_FullPt",k),fBckInvMassVSPtDummy_SubPiZero[k]->GetNbinsX(),0.,1.);
        fMesonFullPtBackground_SubPiZero[k]->Sumw2();
        startBinFullPt = fBckInvMassVSPtDummy_SubPiZero[k]->GetYaxis()->FindBin(fFullPt[0]+0.001);
        endBinFullPt = fBckInvMassVSPtDummy_SubPiZero[k]->GetYaxis()->FindBin(fFullPt[1]-0.001);
        fBckInvMassVSPtDummy_SubPiZero[k]->ProjectionX(Form("Mapping_Back%i_InvMass_SubPiZero_FullPt",k),startBinFullPt,endBinFullPt);
        fMesonFullPtBackground_SubPiZero[k]=(TH1D*)gDirectory->Get(Form("Mapping_Back%i_InvMass_SubPiZero_FullPt",k));
        fMesonFullPtBackground_SubPiZero[k]->Rebin(fNRebin[4]);

        for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
            fNameHistoBack = Form("Mapping_Back%i_InvMass_SubPiZero_in_Pt_Bin%02d",k ,iPt);
            if(fHistoMappingBackInvMassPtBin_SubPiZero[k][iPt]!= NULL){
                delete fHistoMappingBackInvMassPtBin_SubPiZero[k][iPt];
                fHistoMappingBackInvMassPtBin_SubPiZero[k][iPt]=NULL;
            }
            fHistoMappingBackInvMassPtBin_SubPiZero[k][iPt]=new TH1D(fNameHistoBack.Data(),fNameHistoBack.Data(),fBckInvMassVSPtDummy_SubPiZero[k]->GetNbinsX(),0.,1.);
            fHistoMappingBackInvMassPtBin_SubPiZero[k][iPt]->Sumw2();
            Int_t startBin = fBckInvMassVSPtDummy_SubPiZero[k]->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
            Int_t endBin = fBckInvMassVSPtDummy_SubPiZero[k]->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

            fBckInvMassVSPtDummy_SubPiZero[k]->ProjectionX(fNameHistoBack.Data(),startBin,endBin);
            fHistoMappingBackInvMassPtBin_SubPiZero[k][iPt]=(TH1D*)gDirectory->Get(fNameHistoBack.Data());
            if(fNRebin[iPt]>1){
                fHistoMappingBackInvMassPtBin_SubPiZero[k][iPt]->Rebin(fNRebin[iPt]);
            }

        }

        // Begin Calculation for histos with fixed Pz of Pi0
        if(fBckInvMassVSPtDummy_FixedPzPiZero[k]==NULL){
            continue;
        }
        //calculation background for midPt without weighting
        startBinMidPt = fBckInvMassVSPtDummy_FixedPzPiZero[k]->GetYaxis()->FindBin(fMidPt[0]+0.001);
        endBinMidPt = fBckInvMassVSPtDummy_FixedPzPiZero[k]->GetYaxis()->FindBin(fMidPt[1]-0.001);
        fFittingHistMidPtBackground_FixedPzPiZero[k] = new TH1D(Form("Mapping_Back%i_InvMass_FixedPzPiZero_MidPt",k),Form("Mapping_Back%i_InvMass_FixedPzPiZero_MidPt",k),fBckInvMassVSPtDummy_FixedPzPiZero[k]->GetNbinsX(),0.,1.);
        fFittingHistMidPtBackground_FixedPzPiZero[k]->Sumw2();
        fBckInvMassVSPtDummy_FixedPzPiZero[k]->ProjectionX(Form("Mapping_Back%i_InvMass_FixedPzPiZero_MidPt",k),startBinMidPt,endBinMidPt);
        fFittingHistMidPtBackground_FixedPzPiZero[k]=(TH1D*)gDirectory->Get(Form("Mapping_Back%i_InvMass_FixedPzPiZero_MidPt",k));
        fFittingHistMidPtBackground_FixedPzPiZero[k]->Rebin(fNRebin[4]);

        //calulation background for fullPt without weighting
        fMesonFullPtBackground_FixedPzPiZero[k] = new TH1D(Form("Mapping_Back%i_InvMass_FixedPzPiZero_FullPt",k),Form("Mapping_Back%i_InvMass_FixedPzPiZero_FullPt",k),fBckInvMassVSPtDummy_FixedPzPiZero[k]->GetNbinsX(),0.,1.);
        fMesonFullPtBackground_FixedPzPiZero[k]->Sumw2();
        startBinFullPt = fBckInvMassVSPtDummy_FixedPzPiZero[k]->GetYaxis()->FindBin(fFullPt[0]+0.001);
        endBinFullPt = fBckInvMassVSPtDummy_FixedPzPiZero[k]->GetYaxis()->FindBin(fFullPt[1]-0.001);
        fBckInvMassVSPtDummy_FixedPzPiZero[k]->ProjectionX(Form("Mapping_Back%i_InvMass_FixedPzPiZero_FullPt",k),startBinFullPt,endBinFullPt);
        fMesonFullPtBackground_FixedPzPiZero[k]=(TH1D*)gDirectory->Get(Form("Mapping_Back%i_InvMass_FixedPzPiZero_FullPt",k));
        fMesonFullPtBackground_FixedPzPiZero[k]->Rebin(fNRebin[4]);

        for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
            fNameHistoBack = Form("Mapping_Back%i_InvMass_FixedPzPiZero_in_Pt_Bin%02d",k ,iPt);
            if(fHistoMappingBackInvMassPtBin_FixedPzPiZero[k][iPt]!= NULL){
                delete fHistoMappingBackInvMassPtBin_FixedPzPiZero[k][iPt];
                fHistoMappingBackInvMassPtBin_FixedPzPiZero[k][iPt]=NULL;
            }
            fHistoMappingBackInvMassPtBin_FixedPzPiZero[k][iPt]=new TH1D(fNameHistoBack.Data(),fNameHistoBack.Data(),fBckInvMassVSPtDummy_FixedPzPiZero[k]->GetNbinsX(),0.,1.);
            fHistoMappingBackInvMassPtBin_FixedPzPiZero[k][iPt]->Sumw2();
            Int_t startBin = fBckInvMassVSPtDummy_FixedPzPiZero[k]->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
            Int_t endBin = fBckInvMassVSPtDummy_FixedPzPiZero[k]->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

            fBckInvMassVSPtDummy_FixedPzPiZero[k]->ProjectionX(fNameHistoBack.Data(),startBin,endBin);
            fHistoMappingBackInvMassPtBin_FixedPzPiZero[k][iPt]=(TH1D*)gDirectory->Get(fNameHistoBack.Data());
            if(fNRebin[iPt]>1){
                fHistoMappingBackInvMassPtBin_FixedPzPiZero[k][iPt]->Rebin(fNRebin[iPt]);
            }

        }
    }
}

void ProduceBckWithoutWeightingMinimal(TH2D *fBckInvMassVSPtDummy,TH1D **fHistoMappingBackGroupInvMassPtBin){ // Produce Background in bins without midPt and without fullPt without weighting
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoBack = Form("Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
        if(fHistoMappingBackGroupInvMassPtBin[iPt]!= NULL){
            delete fHistoMappingBackGroupInvMassPtBin[iPt];
            fHistoMappingBackGroupInvMassPtBin[iPt]=NULL;
        }
        fHistoMappingBackGroupInvMassPtBin[iPt]=new TH1D(fNameHistoBack.Data(),fNameHistoBack.Data(),fBckInvMassVSPtDummy->GetNbinsX(),0.,1.);
            fHistoMappingBackGroupInvMassPtBin[iPt]->Sumw2();
        Int_t startBin = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
        Int_t endBin = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

        fBckInvMassVSPtDummy->ProjectionX(fNameHistoBack.Data(),startBin,endBin);
        fHistoMappingBackGroupInvMassPtBin[iPt]=(TH1D*)gDirectory->Get(fNameHistoBack.Data());
        if(fNRebin[iPt]>1){
            fHistoMappingBackGroupInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
        }
    }
}

void FillMassHistosArray(TH2D* fGammaGammaInvMassVSPtDummy,TH2D* fGammaGammaInvMassVSPtDummy_SubPiZero,TH2D* fGammaGammaInvMassVSPtDummy_FixedPzPiZero)
{
    fFittingHistMidPtSignal = new TH1D("Mapping_GG_InvMass_MidPt","Mapping_GG_InvMass_MidPt",fGammaGammaInvMassVSPtDummy->GetNbinsX(),0.,(double)fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()));
    if(fGammaGammaInvMassVSPtDummy_SubPiZero!=NULL) fFittingHistMidPtSignal_SubPiZero = new TH1D("Mapping_GG_InvMass_SubPiZero_MidPt","Mapping_GG_InvMass_SubPiZero_MidPt",fGammaGammaInvMassVSPtDummy_SubPiZero->GetNbinsX(),0.
                                                                                                 ,(double)fGammaGammaInvMassVSPtDummy_SubPiZero->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy_SubPiZero->GetNbinsX()));
    if(fGammaGammaInvMassVSPtDummy_FixedPzPiZero!=NULL) fFittingHistMidPtSignal_FixedPzPiZero = new TH1D("Mapping_GG_InvMass_FixedPzPiZero_MidPt","Mapping_GG_InvMass_FixedPzPiZero_MidPt",fGammaGammaInvMassVSPtDummy_FixedPzPiZero->GetNbinsX(),0.
                                                                                                         ,(double)fGammaGammaInvMassVSPtDummy_FixedPzPiZero->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy_FixedPzPiZero->GetNbinsX()));

    // Calculation for normal InvMass
    Int_t startBinMidPt = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fMidPt[0]+0.001);
    Int_t endBinMidPt = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fMidPt[1]-0.001);
        fFittingHistMidPtSignal->Sumw2();
    fGammaGammaInvMassVSPtDummy->ProjectionX("Mapping_GG_InvMass_MidPt",startBinMidPt,endBinMidPt);

    fFittingHistMidPtSignal=(TH1D*)gDirectory->Get("Mapping_GG_InvMass_MidPt");
    fFittingHistMidPtSignal->Rebin(fNRebin[4]);
    //fFittingHistMidPtSignal->Scale(1./fNRebin[4]);
        cout << "Mid pt geschrieben" << endl;

    fMesonFullPtSignal = new TH1D("Mapping_GG_InvMass_FullPt","Mapping_GG_InvMass_FullPt",fGammaGammaInvMassVSPtDummy->GetNbinsX(),0.,fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()));
        fMesonFullPtSignal->Sumw2();
    Int_t startBinFullPt = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fFullPt[0]+0.001);
    Int_t endBinFullPt = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fFullPt[1]-0.001);

    fGammaGammaInvMassVSPtDummy->ProjectionX("Mapping_GG_InvMass_FullPt",startBinFullPt,endBinFullPt);

    fMesonFullPtSignal=(TH1D*)gDirectory->Get("Mapping_GG_InvMass_FullPt");
    fMesonFullPtSignal->Rebin(fNRebin[4]);

    // Calculation for InvMass - InvMass Pi0
    if(fGammaGammaInvMassVSPtDummy_SubPiZero!=NULL){
      Int_t startBinMidPt = fGammaGammaInvMassVSPtDummy_SubPiZero->GetYaxis()->FindBin(fMidPt[0]+0.001);
      Int_t endBinMidPt = fGammaGammaInvMassVSPtDummy_SubPiZero->GetYaxis()->FindBin(fMidPt[1]-0.001);
      fFittingHistMidPtSignal_SubPiZero->Sumw2();
      fGammaGammaInvMassVSPtDummy_SubPiZero->ProjectionX("Mapping_GG_InvMass_SubPiZero_MidPt",startBinMidPt,endBinMidPt);

      fFittingHistMidPtSignal_SubPiZero=(TH1D*)gDirectory->Get("Mapping_GG_InvMass_SubPiZero_MidPt");
      fFittingHistMidPtSignal_SubPiZero->Rebin(fNRebin[4]);
      //fFittingHistMidPtSignal->Scale(1./fNRebin[4]);
      cout << "Mid pt geschrieben" << endl;

      fMesonFullPtSignal_SubPiZero = new TH1D("Mapping_GG_InvMass_SubPiZero_FullPt","Mapping_GG_InvMass_SubPiZero_FullPt",fGammaGammaInvMassVSPtDummy_SubPiZero->GetNbinsX(),0.,fGammaGammaInvMassVSPtDummy_SubPiZero->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy_SubPiZero->GetNbinsX()));
      fMesonFullPtSignal_SubPiZero->Sumw2();
      Int_t startBinFullPt = fGammaGammaInvMassVSPtDummy_SubPiZero->GetYaxis()->FindBin(fFullPt[0]+0.001);
      Int_t endBinFullPt = fGammaGammaInvMassVSPtDummy_SubPiZero->GetYaxis()->FindBin(fFullPt[1]-0.001);

      fGammaGammaInvMassVSPtDummy_SubPiZero->ProjectionX("Mapping_GG_InvMass_SubPiZero_FullPt",startBinFullPt,endBinFullPt);

      fMesonFullPtSignal_SubPiZero=(TH1D*)gDirectory->Get("Mapping_GG_InvMass_SubPiZero_FullPt");
      fMesonFullPtSignal_SubPiZero->Rebin(fNRebin[4]);
    }

    // Calculation for InvMass with fixed pz of pi0
    if(fGammaGammaInvMassVSPtDummy_FixedPzPiZero!=NULL){
      Int_t startBinMidPt = fGammaGammaInvMassVSPtDummy_FixedPzPiZero->GetYaxis()->FindBin(fMidPt[0]+0.001);
      Int_t endBinMidPt = fGammaGammaInvMassVSPtDummy_FixedPzPiZero->GetYaxis()->FindBin(fMidPt[1]-0.001);
      fFittingHistMidPtSignal_FixedPzPiZero->Sumw2();
      fGammaGammaInvMassVSPtDummy_FixedPzPiZero->ProjectionX("Mapping_GG_InvMass_FixedPzPiZero_MidPt",startBinMidPt,endBinMidPt);

      fFittingHistMidPtSignal_FixedPzPiZero=(TH1D*)gDirectory->Get("Mapping_GG_InvMass_FixedPzPiZero_MidPt");
      fFittingHistMidPtSignal_FixedPzPiZero->Rebin(fNRebin[4]);
      //fFittingHistMidPtSignal->Scale(1./fNRebin[4]);
      cout << "Mid pt geschrieben" << endl;

      fMesonFullPtSignal_FixedPzPiZero = new TH1D("Mapping_GG_InvMass_FixedPzPiZero_FullPt","Mapping_GG_InvMass_FixedPzPiZero_FullPt",fGammaGammaInvMassVSPtDummy_FixedPzPiZero->GetNbinsX(),0.,fGammaGammaInvMassVSPtDummy_FixedPzPiZero->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy_FixedPzPiZero->GetNbinsX()));
      fMesonFullPtSignal_FixedPzPiZero->Sumw2();
      Int_t startBinFullPt = fGammaGammaInvMassVSPtDummy_FixedPzPiZero->GetYaxis()->FindBin(fFullPt[0]+0.001);
      Int_t endBinFullPt = fGammaGammaInvMassVSPtDummy_FixedPzPiZero->GetYaxis()->FindBin(fFullPt[1]-0.001);

      fGammaGammaInvMassVSPtDummy_FixedPzPiZero->ProjectionX("Mapping_GG_InvMass_FixedPzPiZero_FullPt",startBinFullPt,endBinFullPt);

      fMesonFullPtSignal_FixedPzPiZero=(TH1D*)gDirectory->Get("Mapping_GG_InvMass_FixedPzPiZero_FullPt");
      fMesonFullPtSignal_FixedPzPiZero->Rebin(fNRebin[4]);
    }

    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){

        // Calculation for norm InvMass
        fNameHistoGG = Form("Mapping_GG_InvMass_in_Pt_Bin%02d", iPt);

        if(fHistoMappingGGInvMassPtBin[iPt]!= NULL){
            delete fHistoMappingGGInvMassPtBin[iPt];
            fHistoMappingGGInvMassPtBin[iPt]=NULL;
        }
        fHistoMappingGGInvMassPtBin[iPt]=new TH1D(fNameHistoGG.Data(),fNameHistoGG.Data(),fGammaGammaInvMassVSPtDummy->GetNbinsX(),0.,fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()));
            fHistoMappingGGInvMassPtBin[iPt]->Sumw2();

        Int_t startBin = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
        Int_t endBin = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

    //             cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
    //             cout<< "bin values::"<< fGammaGammaInvMassVSPtDummy->GetYaxis()->GetBinCenter(startBin)<< " "
    //                 << fGammaGammaInvMassVSPtDummy->GetYaxis()->GetBinCenter(endBin)<< endl;

        fGammaGammaInvMassVSPtDummy->ProjectionX(fNameHistoGG.Data(),startBin,endBin);


        fHistoMappingGGInvMassPtBin[iPt]=(TH1D*)gDirectory->Get(fNameHistoGG.Data());

        if(fNRebin[iPt]>1){
            fHistoMappingGGInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
            //fHistoMappingGGInvMassPtBin[iPt]->Scale(1./fNRebin[iPt]);
        }

        // Calculation for InvMass - InvMass Pi0
        if(fGammaGammaInvMassVSPtDummy_SubPiZero!=NULL){
          fNameHistoGG = Form("Mapping_GG_InvMass_SubPiZero_in_Pt_Bin%02d", iPt);

          if(fHistoMappingGGInvMassPtBin_SubPiZero[iPt]!= NULL){
            delete fHistoMappingGGInvMassPtBin_SubPiZero[iPt];
            fHistoMappingGGInvMassPtBin_SubPiZero[iPt]=NULL;
          }
          fHistoMappingGGInvMassPtBin_SubPiZero[iPt]=new TH1D(fNameHistoGG.Data(),fNameHistoGG.Data(),fGammaGammaInvMassVSPtDummy_SubPiZero->GetNbinsX(),0.,fGammaGammaInvMassVSPtDummy_SubPiZero->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()));
          fHistoMappingGGInvMassPtBin_SubPiZero[iPt]->Sumw2();

          Int_t startBin = fGammaGammaInvMassVSPtDummy_SubPiZero->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
          Int_t endBin = fGammaGammaInvMassVSPtDummy_SubPiZero->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

          //             cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
          //             cout<< "bin values::"<< fGammaGammaInvMassVSPtDummy->GetYaxis()->GetBinCenter(startBin)<< " "
          //                 << fGammaGammaInvMassVSPtDummy->GetYaxis()->GetBinCenter(endBin)<< endl;

          fGammaGammaInvMassVSPtDummy_SubPiZero->ProjectionX(fNameHistoGG.Data(),startBin,endBin);


          fHistoMappingGGInvMassPtBin_SubPiZero[iPt]=(TH1D*)gDirectory->Get(fNameHistoGG.Data());

          if(fNRebin[iPt]>1){
            fHistoMappingGGInvMassPtBin_SubPiZero[iPt]->Rebin(fNRebin[iPt]);
            //fHistoMappingGGInvMassPtBin[iPt]->Scale(1./fNRebin[iPt]);
          }
        }

        // Calculation for Fixed Pz of Pi0
        if(fGammaGammaInvMassVSPtDummy_FixedPzPiZero!=NULL){
          fNameHistoGG = Form("Mapping_GG_InvMass_FixedPzPiZero_in_Pt_Bin%02d", iPt);

          if(fHistoMappingGGInvMassPtBin_FixedPzPiZero[iPt]!= NULL){
            delete fHistoMappingGGInvMassPtBin_FixedPzPiZero[iPt];
            fHistoMappingGGInvMassPtBin_FixedPzPiZero[iPt]=NULL;
          }
          fHistoMappingGGInvMassPtBin_FixedPzPiZero[iPt]=new TH1D(fNameHistoGG.Data(),fNameHistoGG.Data(),fGammaGammaInvMassVSPtDummy_FixedPzPiZero->GetNbinsX(),0.,fGammaGammaInvMassVSPtDummy_FixedPzPiZero->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()));
          fHistoMappingGGInvMassPtBin_FixedPzPiZero[iPt]->Sumw2();

          Int_t startBin = fGammaGammaInvMassVSPtDummy_FixedPzPiZero->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
          Int_t endBin = fGammaGammaInvMassVSPtDummy_FixedPzPiZero->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

          //             cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
          //             cout<< "bin values::"<< fGammaGammaInvMassVSPtDummy->GetYaxis()->GetBinCenter(startBin)<< " "
          //                 << fGammaGammaInvMassVSPtDummy->GetYaxis()->GetBinCenter(endBin)<< endl;

          fGammaGammaInvMassVSPtDummy_FixedPzPiZero->ProjectionX(fNameHistoGG.Data(),startBin,endBin);


          fHistoMappingGGInvMassPtBin_FixedPzPiZero[iPt]=(TH1D*)gDirectory->Get(fNameHistoGG.Data());

          if(fNRebin[iPt]>1){
            fHistoMappingGGInvMassPtBin_FixedPzPiZero[iPt]->Rebin(fNRebin[iPt]);
            //fHistoMappingGGInvMassPtBin[iPt]->Scale(1./fNRebin[iPt]);
          }
        }
    }
    //    cout << "each pt written" << endl;

}

void FillMassMCTrueMesonHistosArray(TH2D* fHistoTrueMesonInvMassVSPtFill)
{
    fHistoTrueMesonInvMassVSPtFill->Sumw2();
for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
    fNameHistoTrue = Form("Mapping_TrueMeson_InvMass_in_Pt_Bin%02d", iPt);
    if(fHistoMappingTrueMesonInvMassPtBins[iPt]!= NULL){
        delete fHistoMappingTrueMesonInvMassPtBins[iPt];
        fHistoMappingTrueMesonInvMassPtBins[iPt]=NULL;
    }
    fHistoMappingTrueMesonInvMassPtBins[iPt] = new TH1D(fNameHistoTrue.Data(),fNameHistoTrue.Data(),fHistoTrueMesonInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueMesonInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueMesonInvMassVSPtFill->GetNbinsX()));
        fHistoMappingTrueMesonInvMassPtBins[iPt]->Sumw2();
    Int_t startBin = fHistoTrueMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
    Int_t endBin = fHistoTrueMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

    //       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
    //       cout<< "bin values::"<< fHistoTrueMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
    //           << fHistoTrueMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

    fHistoTrueMesonInvMassVSPtFill->ProjectionX(fNameHistoTrue.Data(),startBin,endBin,"e");

    fHistoMappingTrueMesonInvMassPtBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrue.Data());
        cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonInvMassPtBins[iPt]->GetEntries() << endl;
    if(fNRebin[iPt]>1){
        fHistoMappingTrueMesonInvMassPtBins[iPt]->Rebin(fNRebin[iPt]);
    }
    fHistoMappingTrueMesonInvMassPtBins[iPt]->SetLineWidth(1);
    fHistoMappingTrueMesonInvMassPtBins[iPt]->SetLineColor(2);
}

}

void FillMassMCTrueReweightedMesonHistosArray(TH2D* fHistoTrueMesonInvMassVSPtFill)
{
cout << "\t \t \tfHistoTrueMesonInvMassVSPtFill="<< fHistoTrueMesonInvMassVSPtFill->GetEntries() << endl;
fHistoTrueMesonInvMassVSPtFill->Sumw2();
for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
    fNameHistoTrue = Form("Mapping_TrueMeson_InvMassReweighted_in_Pt_Bin%02d", iPt);
    if(fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]!= NULL){
        delete fHistoMappingTrueMesonInvMassPtReweightedBins[iPt];
        fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]=NULL;
    }
    fHistoMappingTrueMesonInvMassPtReweightedBins[iPt] = new TH1D(fNameHistoTrue.Data(),fNameHistoTrue.Data(),fHistoTrueMesonInvMassVSPtFill->GetNbinsX(),0.,fHistoTrueMesonInvMassVSPtFill->GetXaxis()->GetBinUpEdge(fHistoTrueMesonInvMassVSPtFill->GetNbinsX()));
    fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->Sumw2();
    Int_t startBin = fHistoTrueMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
    Int_t endBin = fHistoTrueMesonInvMassVSPtFill->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

    //       cout<< "bins::"<< startBin<< " " << endBin<<" "<< fBinsPt[iPt]<< " "<<fBinsPt[iPt+1]<<  endl;
    //       cout<< "bin values::"<< fHistoTrueMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(startBin)<< " "
    //           << fHistoTrueMesonInvMassVSPtFill->GetYaxis()->GetBinCenter(endBin)<< endl;

    fHistoTrueMesonInvMassVSPtFill->ProjectionX(fNameHistoTrue.Data(),startBin,endBin,"e"); 
    fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]=(TH1D*)gDirectory->Get(fNameHistoTrue.Data());
    cout << "bin: " << iPt << "\t Entries in projection: " << fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->GetEntries() << endl;
    if(fNRebin[iPt]>1){
        fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->Rebin(fNRebin[iPt]);
    }
    fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->SetLineWidth(1);
    fHistoMappingTrueMesonInvMassPtReweightedBins[iPt]->SetLineColor(2);
}

}


void CreatePtHistos(){

    fDeltaPt =			 	new TH1D("deltaPt","",fNBinsPt,fBinsPt);
    fDeltaPt->Sumw2();

    // create histos for different integration windows: normal, wide, narrow, left, left wide, left narrow
    for (Int_t k = 0; k < 6; k++){
        // reconstructed yields in different integration windows
        fHistoYieldMeson[k]                    = new TH1D(Form("histoYieldMeson%s",nameIntRange[k].Data()),"",fNBinsPt,fBinsPt);
        fHistoYieldMeson[k]->Sumw2();
        fHistoYieldMesonPerEvent[k]            = new TH1D(Form("histoYieldMeson%sPerEvent",nameIntRange[k].Data()),"",fNBinsPt,fBinsPt);
        fHistoYieldMesonPerEvent[k]->Sumw2();
    }

    fHistoYieldMesonBackFit = 			new TH1D("histoYieldMesonBackFit","",fNBinsPt,fBinsPt);

    // create histos for different integration windows: normal, wide, narrow
    for (Int_t k = 0; k < 3; k++){
        // validated yields
        fHistoYieldTrueMeson[k]                = new TH1D(Form("histoYieldTrueMeson%s",nameIntRange[k].Data()),"",fNBinsPt,fBinsPt);
        fHistoYieldTrueMeson[k]->Sumw2();
        fHistoYieldTrueMesonReweighted[k]      = new TH1D(Form("histoYieldTrueMeson%sReweighted",nameIntRange[k].Data()), "",fNBinsPt,fBinsPt);
        fHistoYieldTrueMesonReweighted[k]->Sumw2();

        // Significance and S/B for reconstructed signal
        fHistoSigndefaultMeson[k]              = new TH1D(Form("histoSigndefault%sMeson",nameIntRange[k].Data()),"",fNBinsPt,fBinsPt);
        fHistoSigndefaultMeson[k]->Sumw2();
        fHistoSBdefaultMeson[k]                = new TH1D(Form("histoSBdefault%sMeson",nameIntRange[k].Data()),"",fNBinsPt,fBinsPt);
        fHistoSBdefaultMeson[k]->Sumw2();
    }

    fHistoYieldMesonPerEventBackFit = 	new TH1D("histoYieldMesonPerEventBackFit","",fNBinsPt,fBinsPt);
    fHistoMassMeson = 			new TH1D("histoMassMeson","",fNBinsPt,fBinsPt);
    fHistoWidthMeson = 			new TH1D("histoWidthMeson","",fNBinsPt,fBinsPt);
    fHistoFWHMMeson = 			new TH1D("histoFWHMMeson","",fNBinsPt,fBinsPt);
    fHistoTrueMassMeson = 		new TH1D("histoTrueMassMeson","",fNBinsPt,fBinsPt);
    fHistoTrueMassMesonReweighted =      new TH1D("histoTrueMassMesonReweighted","",fNBinsPt,fBinsPt);
    fHistoTrueFWHMMeson = 		new TH1D("histoTrueFWHMMeson","",fNBinsPt,fBinsPt);
    fHistoTrueFWHMMesonReweighted =      new TH1D("histoTrueFWHMMesonReweighted","",fNBinsPt,fBinsPt);
    fHistoTrueSignMeson = 		new TH1D("histoTrueSignMeson","",fNBinsPt,fBinsPt);
    fHistoTrueSBMeson = 		new TH1D("histoTrueSBMeson","",fNBinsPt,fBinsPt);
//
    // Histos for normalization at the left of the peak

    fHistoMassMesonLeft = 		new TH1D("histoMassMesonLeft","",fNBinsPt,fBinsPt);
    fHistoWidthMesonLeft = 		new TH1D("histoWidthMesonLeft","",fNBinsPt,fBinsPt);
    fHistoFWHMMesonLeft = 		new TH1D("histoFWHMMesonLeft","",fNBinsPt,fBinsPt);
}

void FillPtHistos()
{
    for(Int_t iPt=fStartPtBin+1;iPt<fNBinsPt+1;iPt++){

        fDeltaPt->SetBinContent(iPt,fBinsPt[iPt]-fBinsPt[iPt-1]);
        fDeltaPt->SetBinError(iPt,0);


        fHistoMassMeson->SetBinContent(iPt,fMesonMass[iPt-1]);
        fHistoMassMeson->SetBinError(iPt,fMesonMassError[iPt-1]);
        // fHistoWidthMeson->SetBinContent(iPt);
        fHistoFWHMMeson->SetBinContent(iPt,fMesonFWHM[iPt-1]);
        fHistoFWHMMeson->SetBinError(iPt,fMesonFWHMError[iPt-1]);

        if (fIsMC) {
            fHistoTrueMassMeson->SetBinContent(iPt,fMesonTrueMass[iPt-1]);
            fHistoTrueMassMeson->SetBinError(iPt,fMesonTrueMassError[iPt-1]);
            fHistoTrueMassMesonReweighted->SetBinContent(iPt,fMesonTrueMassReweighted[iPt-1]);
            fHistoTrueMassMesonReweighted->SetBinError(iPt,fMesonTrueMassReweightedError[iPt-1]);

            fHistoTrueFWHMMeson->SetBinContent(iPt,fMesonTrueFWHM[iPt-1]);
            fHistoTrueFWHMMeson->SetBinError(iPt,fMesonTrueFWHMError[iPt-1]);
            fHistoTrueFWHMMesonReweighted->SetBinContent(iPt,fMesonTrueFWHMReweighted[iPt-1]);
            fHistoTrueFWHMMesonReweighted->SetBinError(iPt,fMesonTrueFWHMReweightedError[iPt-1]);

            fHistoTrueSignMeson->SetBinContent(iPt,fMesonTrueSign[iPt-1]);
            fHistoTrueSignMeson->SetBinError(iPt,fMesonTrueSignError[iPt-1]);
            fHistoTrueSBMeson->SetBinContent(iPt,fMesonTrueSB[iPt-1]);
            fHistoTrueSBMeson->SetBinError(iPt,fMesonTrueSBError[iPt-1]);

            for (Int_t k = 0; k < 3; k++){
                fHistoYieldTrueMeson[k]->SetBinContent(iPt,fMesonTrueYields[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                cout << "-----> SetBinContent (should be filled) = " << fMesonTrueYields[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]) <<endl;
                fHistoYieldTrueMeson[k]->SetBinError(iPt,fMesonTrueYieldsError[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                cout << "-----> SetBinContent = " << fMesonTrueYieldsReweighted[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]) <<endl;
                fHistoYieldTrueMesonReweighted[k]->SetBinContent(iPt,fMesonTrueYieldsReweighted[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
                fHistoYieldTrueMesonReweighted[k]->SetBinError(iPt,fMesonTrueYieldsReweightedError[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            }
        }

        // filling histogram arrays for normal, wide, narrow, left, left wide, left narrow
        for (Int_t k = 0; k < 6; k++){
            fHistoYieldMeson[k]->SetBinContent(iPt,fMesonYieldsCorResidualBckFunc[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldMeson[k]->SetBinError(iPt,fMesonYieldsCorResidualBckFuncError[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldMesonPerEvent[k]->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFunc[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
            fHistoYieldMesonPerEvent[k]->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncError[k][iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        }
        for (Int_t k = 0; k < 3; k++){
            fHistoSigndefaultMeson[k]->SetBinContent(iPt,fMesonSigndefault[k][iPt]-1);
            fHistoSigndefaultMeson[k]->SetBinError(iPt,fMesonSigndefaultError[k][iPt]-1);
            fHistoSBdefaultMeson[k]->SetBinContent(iPt,fMesonSBdefault[k][iPt]-1);
            fHistoSBdefaultMeson[k]->SetBinError(iPt,fMesonSBdefaultError[k][iPt]-1);

        }
        fHistoYieldMesonBackFit->SetBinContent(iPt,fMesonYieldsCorResidualBckFuncBackFit[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonBackFit->SetBinError(iPt,fMesonYieldsCorResidualBckFuncBackFitError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

        fHistoYieldMesonPerEventBackFit->SetBinContent(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncBackFit[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));
        fHistoYieldMesonPerEventBackFit->SetBinError(iPt,(1./fNEvents)*fMesonYieldsCorResidualBckFuncBackFitError[iPt-1]/(fBinsPt[iPt]-fBinsPt[iPt-1]));

        // Histos for integration at the left of the peak
        fHistoMassMesonLeft->SetBinContent(iPt,fMesonMassLeft[iPt-1]); // TODO change to array
        fHistoMassMesonLeft->SetBinError(iPt,fMesonMassLeftError[iPt-1]);
        fHistoFWHMMesonLeft->SetBinContent(iPt,fMesonFWHMLeft[iPt-1]);
        fHistoFWHMMesonLeft->SetBinError(iPt,fMesonFWHMLeftError[iPt-1]);
    }
}

void FitSubtractedInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Double_t* fMesonIntDeltaRangeFit, Int_t ptBin, Bool_t vary, Int_t InvMassType)
{
    /*
     * InvMassType=0 : fHistoMappingSignalInvMassPtBinSingle is normal InvMass distribution
     * InvMassType=1 : fHistoMappingSignalInvMassPtBinSingle is InvMass minus pi0 inv Mass
     * InvMassType=2 : fHistoMappingSignalInvMassPtBinSingle is InvMass with fixed pi0 pz
     */

    // Set Correct Meson Plot Range
    if(InvMassType==0) {
      fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    } else if(InvMassType==1){
      fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange_SubPiZero[0],fMesonMassPlotRange_SubPiZero[1]);
    } else if(InvMassType==2){
      fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange_FixedPzPiZero[0],fMesonMassPlotRange_FixedPzPiZero[1]);
    } else {
      cout << " ERROR! Invalid InvMassType was given. FitSubtractedInvMassInPtBins won't be run!" << endl;
      return;
    }

    Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin;
    Double_t mesonAmplitudeMax;
    if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
        mesonAmplitudeMin = mesonAmplitude*98./100.;
        mesonAmplitudeMax = mesonAmplitude*115./100.;
        if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 || fEnergyFlag.CompareTo("pPb_5.023TeV") == 0) mesonAmplitudeMin = mesonAmplitude*92./100.;
    } else {
        mesonAmplitudeMin = mesonAmplitude*50./100.;
        mesonAmplitudeMax = mesonAmplitude*120./100.;
        if (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5 ||  fMode == 41 || fMode == 42 || fMode == 44 || fMode == 45){
            mesonAmplitudeMin = mesonAmplitude*10./100.;
        }
    }
    Double_t FitRangeTmp[2];
    if(InvMassType==0) {
      FitRangeTmp[0] = fMesonFitRange[0];
      FitRangeTmp[1] = fMesonFitRange[1];
    } else if(InvMassType==1){
      FitRangeTmp[0] = fMesonFitRange_SubPiZero[0];
      FitRangeTmp[1] = fMesonFitRange_SubPiZero[1];
    } else if(InvMassType==2){
      FitRangeTmp[0] = fMesonFitRange_FixedPzPiZero[0];
      FitRangeTmp[1] = fMesonFitRange_FixedPzPiZero[1];
    } else {
      cout << " ERROR! Invalid InvMassType was given. FitSubtractedInvMassInPtBins won't be run!" << endl;
      return;
    }
    fFitReco= NULL;
    fFitReco = new TF1("GaussExpLinear","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",FitRangeTmp[0],FitRangeTmp[1]);

    fFitGausExp =NULL;
    fFitGausExp = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",FitRangeTmp[0],FitRangeTmp[1]);

    fFitLinearBck = NULL;
    fFitLinearBck = new TF1("Linear","[0]+[1]*x",FitRangeTmp[0],FitRangeTmp[1]);


    fFitReco->SetParameter(0,mesonAmplitude);

    if(InvMassType == 1){
      fFitReco->SetParameter(1,fMesonMassExpect-0.134);
    } else{
      fFitReco->SetParameter(1,fMesonMassExpect);
    }

    fFitReco->SetParameter(2,fMesonWidthExpect);
    //    if(vary){
        fFitReco->SetParameter(3,fMesonLambdaTail);
    //    } else {
    //       fFitReco->FixParameter(3,fMesonLambdaTail);
    //    }
    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    //    fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
     if(InvMassType == 1){
       fFitReco->SetParLimits(1,(fMesonMassExpect-0.134)*0.9,(fMesonMassExpect-0.134)*1.15);
     } else{
       fFitReco->SetParLimits(1,fMesonMassExpect*0.9,fMesonMassExpect*1.15);
     }
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
    //    if(vary){
        fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);
    // }

    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");


    fFitReco->SetLineColor(3);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);

    if (vary && !fIsMC){
        fMesonLambdaTail = fFitReco->GetParameter(3);
        fMesonLambdaTailRange[0] = 0.9*fFitReco->GetParameter(3);
        fMesonLambdaTailRange[1] = 1.1*fFitReco->GetParameter(3);
        fMesonWidthExpect = fFitReco->GetParameter(2);
        fMesonWidthRange[0] = 0.5*fFitReco->GetParameter(2);
        fMesonWidthRange[1] = 1.5*fFitReco->GetParameter(2);
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

void FitTrueInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Double_t* fMesonIntDeltaRangeFit, Int_t ptBin, Bool_t vary)
{
    //    cout<<"Start Fitting spectra"<<endl;
    fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange[0],fMesonMassRange[1]);
    Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
    Double_t mesonAmplitudeMin;
    Double_t mesonAmplitudeMax;
    mesonAmplitudeMin = mesonAmplitude*95./100.;
    mesonAmplitudeMax = mesonAmplitude*130./100.;
    if (fMode == 2 || fMode == 3 || fMode == 4 || fMode == 5 || fMode == 41 || fMode == 42 || fMode == 44 || fMode == 45){
        mesonAmplitudeMin = mesonAmplitude*20./100.;
    }
    fFitReco = NULL;
    fFitReco = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

    fFitReco->SetParameter(0,mesonAmplitude);
    fFitReco->SetParameter(1,fMesonMassExpect);
    fFitReco->SetParameter(2,fMesonWidthExpect);
    fFitReco->SetParameter(3,fMesonLambdaTailMC);

    fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
    fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
    fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);
    // //    fFitReco->SetParLimits(1,fMesonMassExpect*0.8,fMesonMassExpect*1.3);
    // //    fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);
    // //    if(vary){
    // 	   fFitReco->SetParLimits(3,fMesonLambdaTailRange[0],fMesonLambdaTailRange[1]);

    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRME0");
    fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"RME0");

    //    cout << TString(gMinuit->fCstatu.Data()).Data() << endl;

    if (vary){
        fMesonLambdaTailMC = fFitReco->GetParameter(3);
    // 	   fMesonWidthExpectMC = fMesonMassExpect*0.03;
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
}

void FitCBSubtractedInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle,Double_t * fMesonIntDeltaRangeFit, Int_t ptBin,Bool_t vary ,TString functionname, Int_t InvMassType)
{

fFileErrLog<<"Start Fitting spectra with CB fit"<<endl;

if(InvMassType==0) {
  fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange[0],fMesonMassRange[1]);
} else if(InvMassType==1){
  fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange_SubPiZero[0],fMesonMassRange_SubPiZero[1]);
} else if(InvMassType==2){
  fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange_FixedPzPiZero[0],fMesonMassRange_FixedPzPiZero[1]);
} else {
  cout << " ERROR! Invalid InvMassType was given. FitCBSubtractedInvMassInPtBins won't be run!" << endl;
  return;
}

Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
Double_t mesonAmplitudeMin = mesonAmplitude*50./100.;
Double_t mesonAmplitudeMax = mesonAmplitude*115./100.;

Double_t FitRangeTmp[2];
if(InvMassType==0) {
  FitRangeTmp[0] = fMesonFitRange[0];
  FitRangeTmp[1] = fMesonFitRange[1];
} else if(InvMassType==1){
  FitRangeTmp[0] = fMesonFitRange_SubPiZero[0];
  FitRangeTmp[1] = fMesonFitRange_SubPiZero[1];
} else if(InvMassType==2){
  FitRangeTmp[0] = fMesonFitRange_FixedPzPiZero[0];
  FitRangeTmp[1] = fMesonFitRange_FixedPzPiZero[1];
} else {
  cout << " ERROR! Invalid InvMassType was given. FitSubtractedInvMassInPtBins won't be run!" << endl;
  return;
}

fFitReco = NULL;
fFitReco = new TF1(functionname,CrystalBallBck,FitRangeTmp[0],FitRangeTmp[1],7);

fFitGausExp = NULL;
fFitGausExp = new TF1("optionCrystalBall",CrystalBall,FitRangeTmp[0],FitRangeTmp[1],5);

fFitLinearBck = NULL;
fFitLinearBck = new TF1("Linear","[0]+[1]*x",FitRangeTmp[0],FitRangeTmp[1]);

fFitReco->SetParameter(0,mesonAmplitude);
if(InvMassType == 1){
  fFitReco->SetParameter(1,fMesonMassExpect-0.134);
} else{
  fFitReco->SetParameter(1,fMesonMassExpect);
}
fFitReco->SetParameter(2,fMesonWidthExpect);
fFitReco->SetParameter(3,2.);
fFitReco->SetParameter(4,0.7);
fFitReco->SetParameter(5,0.);
fFitReco->SetParameter(6,1.);

if (!vary) {
    fFitReco->FixParameter(3,fCBn);
    fFitReco->FixParameter(4,fCBAlpha);
}

fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
if(InvMassType == 0){
  fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
}else if(InvMassType == 1){
  fFitReco->SetParLimits(1,fMesonMassRange_SubPiZero[0],fMesonMassRange_SubPiZero[1]);
}else if(InvMassType == 2){
  fFitReco->SetParLimits(1,fMesonMassRange_FixedPzPiZero[0],fMesonMassRange_FixedPzPiZero[1]);
}else{
  fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
}

fFitReco->SetParLimits(2,fMesonWidthRange[0],fMesonWidthRange[1]);

fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRE0");
fHistoMappingSignalInvMassPtBinSingle->Fit(fFitReco,"QRE0");

fFitReco->SetLineColor(3);
fFitReco->SetLineWidth(1);
fFitReco->SetLineStyle(1);

if (vary) {
    fCBAlpha = fFitReco->GetParameter(4);
    fCBn = fFitReco->GetParError(3);
}

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

fFitLinearBck->SetParameter(0,fFitReco->GetParameter(5));
fFitLinearBck->SetParameter(1,fFitReco->GetParameter(6));

fFitLinearBck->SetParError(0,fFitReco->GetParError(5));
fFitLinearBck->SetParError(1,fFitReco->GetParError(6));

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

    Float_t intLinearBack = fFitLinearBck->GetParameter(0)*(endBinEdge-startBinEdge)+
        0.5*fFitLinearBck->GetParameter(1)*(endBinEdge*endBinEdge-startBinEdge*startBinEdge);

    Float_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fFitReco->GetParError(5),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fFitReco->GetParError(6),2)+2*covMatrix[nFreePar*nFreePar-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5);

    fFileDataLog << "Parameter for bin " << ptBin << endl;
    fFileDataLog << "CrystalBall: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<< "\t "<< fFitReco->GetParameter(4) <<"+-" << fFitReco->GetParError(4)<<endl;
    fFileDataLog << "Linear: \t"<<fFitReco->GetParameter(5)<<"+-" << fFitReco->GetParError(5) << "\t "<<fFitReco->GetParameter(6) <<"+-" << fFitReco->GetParError(6)<< endl;

    fIntLinearBck = intLinearBack/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
    fIntLinearBckError = errorLinearBck/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);
} else {
    fFileErrLog << "Fitting failed in " << ptBin << " with status::" << gMinuit->fCstatu.Data() <<"why failed?"<<endl << endl;
}
fFitReco->DrawCopy("same");

}

void FitWithPol2ForBG(TH1D* fHistoMappingSignalInvMassPtBinSingle,Double_t * fMesonFitRangeCur, Int_t ptBin, Bool_t vary, Int_t InvMassType)
{
//    cout<<"Start Fitting spectra"<<endl;
Int_t startBinSearch = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonFitRangeCur[0]);
Int_t endBinSearch = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonFitRangeCur[1]);
Double_t mesonAmplitude = 0;
for (Int_t i = startBinSearch; i < endBinSearch+1; i++){
    if (mesonAmplitude < fHistoMappingSignalInvMassPtBinSingle->GetBinContent(i)){
        mesonAmplitude = fHistoMappingSignalInvMassPtBinSingle->GetBinContent(i);
    }
}
mesonAmplitude = mesonAmplitude-fHistoMappingSignalInvMassPtBinSingle->GetBinContent(endBinSearch);

fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange[0] ,fMesonMassRange[1]);
Double_t mesonAmplitudeMin;
Double_t mesonAmplitudeMax;
if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
    mesonAmplitudeMin = mesonAmplitude*80./100.;
    mesonAmplitudeMax = mesonAmplitude*115./100.;
} else {
    mesonAmplitudeMin = mesonAmplitude*70./100.;
    mesonAmplitudeMax = mesonAmplitude*110./100.;
}
Double_t fitRange[2];
if (fPrefix.CompareTo("Pi0") ==0 || fPrefix.CompareTo("Pi0EtaBinning")==0 ){
    fitRange[0] = 	fMesonMassRange[0];
    fitRange[1] = 	fMesonMassRange[1];
} else {
    fitRange[0] = 	0.3;
    fitRange[1] = 	0.8;
}

if (fPrefix.CompareTo("Omega") ==0) {fitRange[0]=0.4; fitRange[1]=0.9;}

fFitReco = NULL;
fFitReco = new TF1("GaussExpLinear","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x+[6]*x*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x+[6]*x*x)",fitRange[0] ,fitRange[1]);

fFitGausExp = NULL;
fFitGausExp = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",0.05,0.3);

fFitLinearBck = NULL;
fFitLinearBck = new TF1("Linear","[0]+[1]*x+[2]*x*x",fitRange[0] ,fitRange[1]);

fFitReco->SetParameter(0,mesonAmplitude);
if(InvMassType==0){
fFitReco->SetParameter(1,fMesonMassExpect);
}else if(InvMassType==1){
  fFitReco->SetParameter(1,fMesonMassExpect-0.134);
}else if(InvMassType==2){
  fFitReco->SetParameter(1,fMesonMassExpect);
}else{
  cout << "ERROR: InvMassType " << InvMassType << "is not a valid option!" << endl;
  return;
}
fFitReco->SetParameter(2,fMesonWidthExpect);
if(vary){
    fFitReco->SetParameter(3,fMesonLambdaTail);
} else {
    fFitReco->FixParameter(3,fMesonLambdaTail);
}
fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
if(InvMassType==0){
  fFitReco->SetParLimits(1,fMesonMassExpect*95./100.,fMesonMassExpect*105./100.);
}else if(InvMassType==1){
  fFitReco->SetParLimits(1,(fMesonMassExpect-0.134)*95./100.,(fMesonMassExpect-0.134)*105./100.);
}else if(InvMassType==2){
  fFitReco->SetParLimits(1,fMesonMassExpect*95./100.,fMesonMassExpect*105./100.);
}else{
  cout << "ERROR: InvMassType " << InvMassType << "is not a valid option!" << endl;
  return;
}
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
fFitLinearBck->SetParameter(2,fFitReco->GetParameter(6));

fFitLinearBck->SetParError(0,fFitReco->GetParError(4));
fFitLinearBck->SetParError(1,fFitReco->GetParError(5));
fFitLinearBck->SetParError(2,fFitReco->GetParError(6));

if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 ){
    fFileDataLog << "Parameter for bin " << ptBin << endl;
    fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
    fFileDataLog << "Quadratic: \t"<<fFitReco->GetParameter(4)<<"+-" << fFitReco->GetParError(4) << "\t "<<fFitReco->GetParameter(5) <<"+-" << fFitReco->GetParError(5) << "\t "<<fFitReco->GetParameter(6) <<"+-" << fFitReco->GetParError(6)<< endl;
} else {
    fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
}
fFitReco->DrawCopy("same");
}

void IntegrateHistoInvMass(TH1D * fHistoMappingSignalInvMassPtBinSingle, Double_t * fMesonIntRangeInt)
{
    Int_t binLowMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[0]);
    Int_t binHighMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[1]);
    fYields = fHistoMappingSignalInvMassPtBinSingle->IntegralAndError(binLowMassMeson,binHighMassMeson,fYieldsError);
}

void IntegrateHistoInvMassStream(TH1D * fHistoMappingSignalInvMassPtBinSingle, Double_t * fMesonIntRangeInt)
{
    Int_t binLowMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[0]);
    Int_t binHighMassMeson = fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[1]);
    fYields = fHistoMappingSignalInvMassPtBinSingle->IntegralAndError(binLowMassMeson,binHighMassMeson,fYieldsError);
    for ( Int_t M = binLowMassMeson; M < binHighMassMeson+1; M++){
        fFileDataLog << M << "\t" << fHistoMappingSignalInvMassPtBinSingle->GetBinCenter(M) <<"\t" <<fHistoMappingSignalInvMassPtBinSingle->GetBinContent(M)<< "+-"<< fHistoMappingSignalInvMassPtBinSingle->GetBinError(M)<< endl;
    }
}

void IntegrateFitFunc(TF1 * fFunc, TH1D *  fHistoMappingSignalInvMassPtBinSingle,Double_t * fMesonIntRangeInt)
{

    fYieldsFunc = fFunc->Integral(fMesonIntRangeInt[0],fMesonIntRangeInt[1])/fHistoMappingSignalInvMassPtBinSingle->GetBinWidth(10);

}

void FillHistosArrayMC(TH1D* fHistoMCMesonPtWithinAcceptanceFill, TH1D * fHistoMCMesonPtFill, TH1D * fDeltaPtFill, TString nameCanvas)
{
//    Char_t nameHisto[100] = "fHistoMCMesonPtEtaWithinAcceptance";
    fHistoMCMesonPtWithinAcceptanceFill->Sumw2();
    fHistoMCMesonWithinAccepPt = (TH1D*)fHistoMCMesonPtWithinAcceptanceFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
    fHistoMCMesonWithinAccepPt->Divide(fDeltaPtFill);
    fHistoMCMesonPtFill->Sumw2();
    fHistoMCMesonPt1 = (TH1D*)fHistoMCMesonPtFill->Rebin(fNBinsPt,"",fBinsPt); // Proper bins in Pt
    fHistoMCMesonPt1->Divide(fDeltaPtFill);

}

void CalculateMesonAcceptance(){
    fHistoMCMesonAcceptPt = new TH1D("fMCMesonAccepPt","",fNBinsPt,fBinsPt);
    fHistoMCMesonAcceptPt->Sumw2();

    fHistoMCMesonAcceptPt->Divide(fHistoMCMesonWithinAccepPt,fHistoMCMesonPt1,1.,1.,"B");
    fHistoMCMesonAcceptPt->DrawCopy();
    fFileDataLog << endl << "Calculation of the Acceptance" << endl;
    for ( Int_t i = 1; i < fHistoMCMesonAcceptPt->GetNbinsX()+1 ; i++){
        fFileDataLog << "Bin " << i << "\t"<< fHistoMCMesonAcceptPt->GetBinCenter(i)<< "\t" << fHistoMCMesonAcceptPt->GetBinContent(i) << "\t" << fHistoMCMesonAcceptPt->GetBinError(i) <<endl;
    }
}

TH1D* CalculateMesonEfficiency( TH1D* fMC_fMesonYieldsPt,
                                TH1D* fHistoMCMesonWithinAccepPt,
                                TString nameEfi
                              ){
    cout << "---> Begin of CalculateMesonEfficiency " << endl;
    fHistoMCMesonEffiPt = new TH1D(nameEfi.Data(),"",fNBinsPt,fBinsPt);

    fHistoMCMesonEffiPt->Sumw2();
    fHistoMCMesonEffiPt->Add(fMC_fMesonYieldsPt,1.);
    //cout<<"ADDING 10"<<endl;


    TH1D* fHistoMCMesonEffiScaledByPt = (TH1D*)fHistoMCMesonEffiPt->Clone("ScaledByPt");
    for (Int_t i = 1; i < fHistoMCMesonEffiScaledByPt->GetNbinsX()+1; i++){
        fHistoMCMesonEffiScaledByPt->SetBinContent(i, fHistoMCMesonEffiScaledByPt->GetBinContent(i)/fNEvents);//fHistoMCMesonEffiScaledByPt->GetBinCenter(i)/2/TMath::Pi()
        fHistoMCMesonEffiScaledByPt->SetBinError(i, fHistoMCMesonEffiScaledByPt->GetBinError(i)/fNEvents	);//fHistoMCMesonEffiScaledByPt->GetBinCenter(i)/2/TMath::Pi()
    }

    fHistoMCMesonEffiPt->Divide(fHistoMCMesonEffiPt,fHistoMCMesonWithinAccepPt,1.,1.,"B");
    fFileDataLog << endl << "Calculation of the Efficiency" << nameEfi.Data()<< endl;
    for ( Int_t i = 1; i < fHistoMCMesonEffiPt->GetNbinsX()+1 ; i++){
        fFileDataLog << "Bin " << i << "\t" << fHistoMCMesonEffiPt->GetBinCenter(i)<< "\t"<< fHistoMCMesonEffiPt->GetBinContent(i) << "\t" << fHistoMCMesonEffiPt->GetBinError(i) <<endl;
    }
    fHistoMCMesonEffiFitPt = (TH1D*)fHistoMCMesonEffiPt->Clone(Form("%s_Fit",nameEfi.Data()));

    return fHistoMCMesonEffiPt;
}

void SaveHistos(Int_t optionMC, TString fCutID, TString fPrefix3)
{
    const char* nameOutput = Form("%s/%s/%s_%s_GammaConvV1WithoutCorrection%s_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPrefix3.Data(),fPeriodFlag.Data(),fCutID.Data());
    fOutput1 = new TFile(nameOutput,"RECREATE");

    cout << "Begin writing Uncorrected File" << endl;

    // write histograms for all integration windows: normal, wide, narrow, left, left wide, left narrow
    for (Int_t k = 0; k < 6; k++){
        if (fHistoYieldMeson[k])            fHistoYieldMeson[k]->Write();
        if (fHistoYieldMesonPerEvent[k])    fHistoYieldMesonPerEvent[k]->Write();
    }

    fHistoYieldMesonBackFit->Write();
    fHistoYieldMesonPerEventBackFit->Write();
    // write histograms for integration windows: normal, wide, narrow
    for (Int_t k = 0; k < 3; k++){
        if (fHistoSigndefaultMeson[k])      fHistoSigndefaultMeson[k]->Write();
        if (fHistoSBdefaultMeson[k])        fHistoSBdefaultMeson[k]->Write();
    }

    fHistoMassMeson->Write();
    fHistoWidthMeson->Write();
    fHistoFWHMMeson->Write();
    fDeltaPt->Write();


    fHistoMassMesonLeft->Write();
    fHistoWidthMesonLeft->Write();
    fHistoFWHMMesonLeft->Write();
    fMesonFullPtSignal->Write();
    fMesonFullPtSignal_SubPiZero->Write();
    fMesonFullPtSignal_FixedPzPiZero->Write();

    fMesonFullPtBackNorm->SetName("Mapping_BackNorm_InvMass_FullPt");
    fMesonFullPtBackNorm_SubPiZero->SetName("Mapping_BackNorm_InvMass_SubPiZero_FullPt");
    fMesonFullPtBackNorm_FixedPzPiZero->SetName("Mapping_BackNorm_InvMass_FixedPzPiZero_FullPt");
    fMesonFullPtBackNorm->Write();
    fMesonFullPtBackNorm_SubPiZero->Write();
    fMesonFullPtBackNorm_FixedPzPiZero->Write();

    fNumberOfGoodESDTracks->Write();
    fEventQuality->Write();

    TString nameHistoSignal;
    TString nameHistoBckNorm;
    TString fitnameSignal;
    TString nameHistoSignalPos;
    for(Int_t k=0;k<5;k++){
        fMesonFullPtBackground[k]->Write();
        for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
            nameHistoBckNorm = Form("Mapping_BckNorm_Group%i_InvMass_in_Pt_Bin%02d",k, ii);
            fHistoMappingBackNormInvMassPtBin[k][ii]->Write(nameHistoBckNorm.Data());
        }
    }
    for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
        fHistoMappingGGInvMassPtBin[ii]->Write();
        fHistoMappingGGInvMassPtBin_SubPiZero[ii]->Write();
        fHistoMappingGGInvMassPtBin_FixedPzPiZero[ii]->Write();
        nameHistoSignal = Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d", ii);
        fHistoMappingSignalInvMassPtBin[ii]->Write(nameHistoSignal.Data());
        nameHistoSignal = Form("fHistoMappingSignalInvMass_SubPiZero_in_Pt_Bin%02d", ii);
        fHistoMappingSignalInvMassPtBin_SubPiZero[ii]->Write(nameHistoSignal.Data());
        nameHistoSignal = Form("fHistoMappingSignalInvMass_FixedPzPiZero_in_Pt_Bin%02d", ii);
        fHistoMappingSignalInvMassPtBin_FixedPzPiZero[ii]->Write(nameHistoSignal.Data());
        fitnameSignal = Form("Signal_InvMassFit_in_Pt_Bin%02d", ii);
        if(fFitSignalInvMassPtBin[ii]!=0x00) fFitSignalInvMassPtBin[ii]->Write(fitnameSignal.Data());

        fitnameSignal = Form("Signal_InvMassFit_SubPiZero_in_Pt_Bin%02d", ii);
        if(fFitSignalInvMassPtBin_SubPiZero[ii]!=0x00) fFitSignalInvMassPtBin_SubPiZero[ii]->Write(fitnameSignal.Data());

        fitnameSignal = Form("Signal_InvMassFit_FixedPzPiZero_in_Pt_Bin%02d", ii);
        if(fFitSignalInvMassPtBin_FixedPzPiZero[ii]!=0x00) fFitSignalInvMassPtBin_FixedPzPiZero[ii]->Write(fitnameSignal.Data());
    }

    if(optionMC){
        fHistoTrueSignMeson->Write();
        fHistoTrueSBMeson->Write();
        fHistoMCMesonPtWithinAcceptance->Write();
        fHistoMCMesonWithinAccepPt->Write(); // Proper bins in Pt
        fHistoMCMesonPt1->Write(); // Proper bins in Pt

        for (Int_t k = 0; k < 3; k++){ // different integration windows: normal, wide, narrow
            if (fHistoYieldTrueMeson[k])            fHistoYieldTrueMeson[k]->Write();
            if (fHistoYieldTrueMesonReweighted[k])  fHistoYieldTrueMesonReweighted[k]->Write();
        }

        fHistoTrueMesonInvMassVSPt->Write();
        for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
            fHistoMappingTrueMesonInvMassPtBins[ii]->Write();
            fHistoMappingTrueMesonInvMassPtReweightedBins[ii]->Write();
            if (fFitTrueSignalInvMassPtBin[ii]!=0x00) fFitTrueSignalInvMassPtBin[ii]->Write();
            if (fFitTrueSignalInvMassPtReweightedBin[ii]!=0x00) fFitTrueSignalInvMassPtReweightedBin[ii]->Write();
        }
    }

    cout << "End writing Uncorrected File" << endl;

    fOutput1->Write();
    fOutput1->Close();
}

void SaveCorrectionHistos(TString fCutID, TString fPrefix3)
{
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
    fHistoMonteMesonEffiBackFitPt->Write();

    // write efficiencies depending on integration windows: normal, wide, narrow, left, left wide, left narrow
    for (Int_t k = 0; k < 6; k++){
        if (fHistoMonteMesonEffiPt[k])          fHistoMonteMesonEffiPt[k]->Write();
        if (k < 3){
            if (fHistoMCTrueMesonEffiPt[k])                 fHistoMCTrueMesonEffiPt[k]->Write();
            if (fHistoMCTrueMesonEffiPtReweighted[k])       fHistoMCTrueMesonEffiPtReweighted[k]->Write();
        }
    }

    fHistoTrueMassMeson->Write();
    fHistoTrueMassMesonReweighted->Write();
    fHistoTrueFWHMMeson->Write();
    fHistoTrueFWHMMesonReweighted->Write();
    fHistoMCMesonPt1->SetName("MC_Meson_genPt");
    fHistoMCMesonPt1->Write(); // Proper bins in Pt
    fHistoMCMesonPt->SetName("MC_Meson_genPt_oldBin");
    fHistoMCMesonPt->Scale(1./fHistoMCMesonPt->GetBinWidth(5));
    //    fHistoMCMesonPt->GetXaxis()->SetRangeUser(0.,25.);
    fHistoMCMesonPt->Write(); // Proper bins in Pt
    if (fHistoMCMesonPtWOWeights){
        fHistoMCMesonPtWOWeights->SetName("MC_Meson_genPt_WOWeights");
        fHistoMCMesonPtWOWeights->Scale(1./fHistoMCMesonPtWOWeights->GetBinWidth(5));
    //       fHistoMCMesonPtWOWeights->GetXaxis()->SetRangeUser(0.,25.);
        fHistoMCMesonPtWOWeights->Write(); // Proper bins in Pt
    }
    if (fHistoMCMesonPtWeights){
        fHistoMCMesonPtWeights->Write("MC_Meson_genPt_Weights"); // Proper bins in Pt
    }
    fEventQuality->Write();

    // write reconstructed yield for reconstructed mesons for different integration ranges: normal, wide, narrow
    for (Int_t k = 0; k < 3; k++){
        // primary yield
        if (fHistoYieldTrueMeson[k])            fHistoYieldTrueMeson[k]->Write();
        if (fHistoYieldTrueMesonReweighted[k])  fHistoYieldTrueMesonReweighted[k]->Write();
    }

    cout << "end writing Correction File" << endl;

    fOutput2->Write();
    fOutput2->Close();
}

void Initialize(TString setPi0, Int_t numberOfBins){

    //cout << "MODE in INITIALIZE = " << fMode << endl;
    // New way to initialize
    InitializeBinning(setPi0, numberOfBins, fEnergyFlag, kFALSE, fMode, fEventCutSelection, fClusterCutSelection);
    InitializeWindows(setPi0, fMode, ""); // trigger left empty

    // initialize integration array for integration windows: normal, wide, narrow, left, left wide, left narrow
    for (Int_t k = 0; k < 6; k++){
        fMesonCurIntRange[k]                                        = new Double_t[2];
        fMesonCurIntRange_SubPiZero[k]                              = new Double_t[2];
        fMesonCurIntRange_FixedPzPiZero[k]                          = new Double_t[2];
    }

    // initialize integration pt-arrays for integration windows: normal, wide, narrow, left, left wide, left narrow
    for (Int_t k = 0; k < 6; k++){
        // Initialize yield arrays
        fGGYields[k] = 											new Double_t[fNBinsPt];
        fGGYields_SubPiZero[k] =                				new Double_t[fNBinsPt];
        fGGYields_FixedPzPiZero[k] =        					new Double_t[fNBinsPt];
        fBckYields[k] =                                         new Double_t[fNBinsPt];
        fBckYields_SubPiZero[k] =                               new Double_t[fNBinsPt];
        fBckYields_FixedPzPiZero[k] =                           new Double_t[fNBinsPt];
        fTotalBckYields[k] = 								    new Double_t[fNBinsPt];
        fTotalBckYields_SubPiZero[k] =						    new Double_t[fNBinsPt];
        fTotalBckYields_FixedPzPiZero[k] =					    new Double_t[fNBinsPt];
        fMesonYields[k] = 									    new Double_t[fNBinsPt];
        fMesonYields_SubPiZero[k] =							    new Double_t[fNBinsPt];
        fMesonYields_FixedPzPiZero[k] =						    new Double_t[fNBinsPt];
        fMesonYieldsResidualBckFunc[k] = 						new Double_t[fNBinsPt];
        fMesonYieldsResidualBckFunc_SubPiZero[k] =				new Double_t[fNBinsPt];
        fMesonYieldsResidualBckFunc_FixedPzPiZero[k] =			new Double_t[fNBinsPt];
        fMesonYieldsCorResidualBckFunc[k] = 					new Double_t[fNBinsPt];
        fMesonYieldsCorResidualBckFunc_SubPiZero[k] = 			new Double_t[fNBinsPt];
        fMesonYieldsCorResidualBckFunc_FixedPzPiZero[k] =		new Double_t[fNBinsPt];
        fMesonYieldsPerEvent[k] = 								new Double_t[fNBinsPt];
        fMesonYieldsPerEvent_SubPiZero[k] =						new Double_t[fNBinsPt];
        fMesonYieldsPerEvent_FixedPzPiZero[k] =     			new Double_t[fNBinsPt];

        // Initialize error arrays
        fGGYieldsError[k] = 								    new Double_t[fNBinsPt];
        fGGYieldsError_SubPiZero[k] = 						    new Double_t[fNBinsPt];
        fGGYieldsError_FixedPzPiZero[k] = 					    new Double_t[fNBinsPt];
        fBckYieldsError[k] = 									new Double_t[fNBinsPt];
        fBckYieldsError_SubPiZero[k] =  						new Double_t[fNBinsPt];
        fBckYieldsError_FixedPzPiZero[k] =						new Double_t[fNBinsPt];
        fTotalBckYieldsError[k] = 								new Double_t[fNBinsPt];
        fTotalBckYieldsError_SubPiZero[k] = 					new Double_t[fNBinsPt];
        fTotalBckYieldsError_FixedPzPiZero[k] =					new Double_t[fNBinsPt];
        fMesonYieldsError[k] = 									new Double_t[fNBinsPt];
        fMesonYieldsError_SubPiZero[k] =						new Double_t[fNBinsPt];
        fMesonYieldsError_FixedPzPiZero[k] =					new Double_t[fNBinsPt];
        fMesonYieldsResidualBckFuncError[k] = 					new Double_t[fNBinsPt]; 
        fMesonYieldsResidualBckFuncError_SubPiZero[k] = 		new Double_t[fNBinsPt];
        fMesonYieldsResidualBckFuncError_FixedPzPiZero[k] = 	new Double_t[fNBinsPt];
        fMesonYieldsCorResidualBckFuncError[k] =				new Double_t[fNBinsPt];
        fMesonYieldsCorResidualBckFuncError_SubPiZero[k] =    	new Double_t[fNBinsPt];
        fMesonYieldsCorResidualBckFuncError_FixedPzPiZero[k] =	new Double_t[fNBinsPt];
        fMesonYieldsPerEventError[k] = 							new Double_t[fNBinsPt];
        fMesonYieldsPerEventError_SubPiZero[k] =				new Double_t[fNBinsPt];
        fMesonYieldsPerEventError_FixedPzPiZero[k] =			new Double_t[fNBinsPt];
    }

    // initialize variable for different integration windows: normal, wide, narrow
    for (Int_t k = 0; k < 3; k++ ){
        // initialize pt arrays for integration ranges
        fMesonTrueIntRange[k]                                        = new Double_t[2];
        fMesonTrueIntReweightedRange[k]                              = new Double_t[2];

        // initialize pt arrays for reconstructed validated yields
        fMesonTrueYields[k]                                          = new Double_t[fNBinsPt];
        fMesonTrueYieldsReweighted[k]                                = new Double_t[fNBinsPt];
        fMesonTrueYieldsError[k]                                     = new Double_t[fNBinsPt];
        fMesonTrueYieldsReweightedError[k]                           = new Double_t[fNBinsPt];

        fMesonSBdefault[k]                                           = new Double_t[fNBinsPt];
        fMesonSigndefault[k]                                         = new Double_t[fNBinsPt];
        fMesonSBdefaultError[k]                                      = new Double_t[fNBinsPt];
        fMesonSigndefaultError[k]                                    = new Double_t[fNBinsPt];

    }
    fMesonCurIntRangeBackFit = 								new Double_t[2];
    fMesonCurIntRangeBackFit_SubPiZero = 					new Double_t[2];
    fMesonCurIntRangeBackFit_FixedPzPiZero = 				new Double_t[2];
    fMesonYieldsBackFit =									new Double_t[fNBinsPt];
    fMesonYieldsBackFit_SubPiZero =           				new Double_t[fNBinsPt];
    fMesonYieldsBackFit_FixedPzPiZero =						new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncBackFit = 					new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncBackFit_SubPiZero = 			new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncBackFit_FixedPzPiZero =		new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncBackFit = 				new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncBackFit_SubPiZero = 		new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncBackFit_FixedPzPiZero =	new Double_t[fNBinsPt];
    fMesonYieldsPerEventBackFit = 							new Double_t[fNBinsPt];
    fMesonYieldsPerEventBackFit_SubPiZero = 				new Double_t[fNBinsPt];
    fMesonYieldsPerEventBackFit_FixedPzPiZero =				new Double_t[fNBinsPt];
    fMesonMass = 											new Double_t[fNBinsPt];
    fMesonMass_SubPiZero =									new Double_t[fNBinsPt];
    fMesonMass_FixedPzPiZero =								new Double_t[fNBinsPt];
    fMesonMassBackFit = 									new Double_t[fNBinsPt];
    fMesonMassBackFit_SubPiZero =							new Double_t[fNBinsPt];
    fMesonMassBackFit_FixedPzPiZero =						new Double_t[fNBinsPt];
    fMesonWidth = 											new Double_t[fNBinsPt];
    fMesonWidth_SubPiZero =									new Double_t[fNBinsPt];
    fMesonWidth_FixedPzPiZero =								new Double_t[fNBinsPt];
    fMesonWidthBackFit = 									new Double_t[fNBinsPt];
    fMesonWidthBackFit_SubPiZero =							new Double_t[fNBinsPt];
    fMesonWidthBackFit_FixedPzPiZero =						new Double_t[fNBinsPt];
    fMesonFWHM = 											new Double_t[fNBinsPt];
    fMesonFWHM_SubPiZero =									new Double_t[fNBinsPt];
    fMesonFWHM_FixedPzPiZero =								new Double_t[fNBinsPt];
    fMesonTrueMass = 										new Double_t[fNBinsPt];
    fMesonTrueMassReweighted = 								new Double_t[fNBinsPt];
    fMesonTrueFWHM = 										new Double_t[fNBinsPt];
    fMesonTrueFWHMReweighted =								new Double_t[fNBinsPt];
    fMesonTrueSB = 											new Double_t[fNBinsPt];
    fMesonTrueSBError =                                     new Double_t[fNBinsPt];
    fMesonTrueSign = 										new Double_t[fNBinsPt];
    fMesonTrueSignError = 									new Double_t[fNBinsPt];

    // Normalization at the left of the peak
    fMesonMassLeft = 										new Double_t[fNBinsPt];
    fMesonMassLeft_SubPiZero =                  			new Double_t[fNBinsPt];
    fMesonMassLeft_FixedPzPiZero = 							new Double_t[fNBinsPt];
    fMesonWidthLeft = 										new Double_t[fNBinsPt];
    fMesonWidthLeft_SubPiZero = 							new Double_t[fNBinsPt];
    fMesonWidthLeft_FixedPzPiZero = 						new Double_t[fNBinsPt];
    fMesonFWHMLeft = 										new Double_t[fNBinsPt];
    fMesonFWHMLeft_SubPiZero = 								new Double_t[fNBinsPt];
    fMesonFWHMLeft_FixedPzPiZero = 							new Double_t[fNBinsPt];

    fMesonYieldsBackFitError =								new Double_t[fNBinsPt];
    fMesonYieldsBackFitError_SubPiZero =					new Double_t[fNBinsPt];
    fMesonYieldsBackFitError_FixedPzPiZero =				new Double_t[fNBinsPt];

    fMesonYieldsResidualBckFuncBackFitError = 				new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncBackFitError_SubPiZero = 	new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncBackFitError_FixedPzPiZero =	new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncBackFitError =			new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncBackFitError_SubPiZero =		new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncBackFitError_FixedPzPiZero =	new Double_t[fNBinsPt];
    fMesonYieldsPerEventBackFitError = 						new Double_t[fNBinsPt];
    fMesonYieldsPerEventBackFitError_SubPiZero =			new Double_t[fNBinsPt];
    fMesonYieldsPerEventBackFitError_FixedPzPiZero = 		new Double_t[fNBinsPt];
    fMesonMassError = 										new Double_t[fNBinsPt];
    fMesonMassError_SubPiZero =								new Double_t[fNBinsPt];
    fMesonMassError_FixedPzPiZero =							new Double_t[fNBinsPt];
    fMesonMassBackFitError = 								new Double_t[fNBinsPt];
    fMesonMassBackFitError_SubPiZero =  					new Double_t[fNBinsPt];
    fMesonMassBackFitError_FixedPzPiZero =					new Double_t[fNBinsPt];
    fMesonWidthError = 										new Double_t[fNBinsPt];
    fMesonWidthError_SubPiZero = 							new Double_t[fNBinsPt];
    fMesonWidthError_FixedPzPiZero =      					new Double_t[fNBinsPt];
    fMesonWidthBackFitError = 								new Double_t[fNBinsPt];
    fMesonWidthBackFitError_SubPiZero = 					new Double_t[fNBinsPt];
    fMesonWidthBackFitError_FixedPzPiZero = 				new Double_t[fNBinsPt];
    fMesonFWHMError = 										new Double_t[fNBinsPt];
    fMesonFWHMError_SubPiZero = 							new Double_t[fNBinsPt];
    fMesonFWHMError_FixedPzPiZero =							new Double_t[fNBinsPt];
    fMesonTrueMassError = 									new Double_t[fNBinsPt];
    fMesonTrueMassReweightedError =		 					new Double_t[fNBinsPt];
    fMesonTrueFWHMError = 									new Double_t[fNBinsPt];
    fMesonTrueFWHMReweightedError = 						new Double_t[fNBinsPt];

    fMesonMassLeftError = 									new Double_t[fNBinsPt];
    fMesonMassLeftError_SubPiZero =							new Double_t[fNBinsPt];
    fMesonMassLeftError_FixedPzPiZero =						new Double_t[fNBinsPt];
    fMesonWidthLeftError = 									new Double_t[fNBinsPt];
    fMesonWidthLeftError_SubPiZero = 						new Double_t[fNBinsPt];
    fMesonWidthLeftError_FixedPzPiZero =					new Double_t[fNBinsPt];
    fMesonFWHMLeftError = 									new Double_t[fNBinsPt];
    fMesonFWHMLeftError_SubPiZero = 						new Double_t[fNBinsPt];
    fMesonFWHMLeftError_FixedPzPiZero =						new Double_t[fNBinsPt];

    // nitialize pt-arrays for different mass & width fitting procedures, normalization at the left of the peak
    fMesonMassLeft                                                  = new Double_t[fNBinsPt];
    fMesonMassLeft_SubPiZero                                        = new Double_t[fNBinsPt];
    fMesonMassLeft_FixedPzPiZero                                    = new Double_t[fNBinsPt];
    fMesonFWHMLeft                                                  = new Double_t[fNBinsPt];
    fMesonFWHMLeft_SubPiZero                                        = new Double_t[fNBinsPt];
    fMesonFWHMLeft_FixedPzPiZero                                    = new Double_t[fNBinsPt];
    fMesonMassLeftError                                             = new Double_t[fNBinsPt];
    fMesonMassLeftError_SubPiZero                                   = new Double_t[fNBinsPt];
    fMesonMassLeftError_FixedPzPiZero                                 = new Double_t[fNBinsPt];
    fMesonFWHMLeftError                                             = new Double_t[fNBinsPt];
    fMesonFWHMLeftError_SubPiZero                                   = new Double_t[fNBinsPt];
    fMesonFWHMLeftError_FixedPzPiZero                               = new Double_t[fNBinsPt];

    // initialize pt-arrays for validated Significance and S/B
    fMesonTrueSB                                                    = new Double_t[fNBinsPt];
    fMesonTrueSign                                                  = new Double_t[fNBinsPt];

    fHistoMappingTrueMesonInvMassPtBins = 					new TH1D*[fNBinsPt];
    fHistoMappingTrueMesonInvMassPtReweightedBins =  		new TH1D*[fNBinsPt];

    fHistoWeightsBGZbinVsMbin = 							new TH2F*[fNBinsPt];
    fHistoFillPerEventBGZbinVsMbin = 						new TH2F*[fNBinsPt];

    fHistoMappingGGInvMassPtBin = 							new TH1D*[fNBinsPt];
    fHistoMappingGGInvMassPtBin_SubPiZero = 				new TH1D*[fNBinsPt];
    fHistoMappingGGInvMassPtBin_FixedPzPiZero =				new TH1D*[fNBinsPt];

    fHistoMappingGGInvMassBackFitPtBin =					new TH1D*[fNBinsPt];
    fHistoMappingGGInvMassBackFitPtBin_SubPiZero =   		new TH1D*[fNBinsPt];
    fHistoMappingGGInvMassBackFitPtBin_FixedPzPiZero =		new TH1D*[fNBinsPt];
    fHistoMappingGGInvMassBackFitWithoutSignalPtBin =		new TH1D*[fNBinsPt];
    fHistoMappingSignalInvMassPtBin =				    new TH1D*[fNBinsPt];
    fHistoMappingSignalInvMassPtBin_SubPiZero =		    new TH1D*[fNBinsPt];
    fHistoMappingSignalInvMassPtBin_FixedPzPiZero =	    new TH1D*[fNBinsPt];
    for(Int_t k = 0; k<5; k++){
        fHistoMappingBackInvMassPtBin[k] = 						new TH1D*[fNBinsPt];
        fHistoMappingBackInvMassPtBin_SubPiZero[k] =			new TH1D*[fNBinsPt];
        fHistoMappingBackInvMassPtBin_FixedPzPiZero[k] =		new TH1D*[fNBinsPt];
        fHistoMappingBackNormInvMassPtBin[k] = 					new TH1D*[fNBinsPt];
        fHistoMappingBackNormInvMassPtBin_SubPiZero[k] =		new TH1D*[fNBinsPt];
        fHistoMappingBackNormInvMassPtBin_FixedPzPiZero[k] =	new TH1D*[fNBinsPt];
        fHistoMappingBackSameNormInvMassPtBin[k] =              new TH1D*[fNBinsPt];
        fHistoMappingBackSameNormInvMassPtBin_SubPiZero[k]=     new TH1D*[fNBinsPt];
        fHistoMappingBackSameNormInvMassPtBin_FixedPzPiZero[k]= new TH1D*[fNBinsPt];
        fHistoMappingRatioSBInvMassPtBin[k]= 				    new TH1D*[fNBinsPt];
        fHistoMappingRatioSBInvMassPtBin_SubPiZero[k]= 		    new TH1D*[fNBinsPt];
        fHistoMappingRatioSBInvMassPtBin_FixedPzPiZero[k]=	    new TH1D*[fNBinsPt];
        for(Int_t i = 0;i<fNBinsPt; i++){
            fHistoMappingBackInvMassPtBin[k][i] = 							NULL;
            fHistoMappingBackInvMassPtBin_SubPiZero[k][i] = 				NULL;
            fHistoMappingBackInvMassPtBin_FixedPzPiZero[k][i] = 			NULL;
            fHistoMappingBackNormInvMassPtBin[k][i] = 						NULL;
    //      fRatio[i] =                                                 NULL;
            fHistoMappingRatioSBInvMassPtBin[k][i] = 						NULL;
            fHistoMappingRatioSBInvMassPtBin_SubPiZero[k][i] = 						NULL;
            fHistoMappingRatioSBInvMassPtBin_FixedPzPiZero[k][i] = 						NULL;
        }
    }
//  fRatio =                                                new TH1D*[fNBinsPt];

    fBackgroundFitPol = 									new TF1*[fNBinsPt];
    fBackgroundFitPol_SubPiZero = 							new TF1*[fNBinsPt];
    fBackgroundFitPol_FixedPzPiZero =						new TF1*[fNBinsPt];
    fHistoBckFitConfidence =                                new TH1D*[fNBinsPt];
    fHistoBckFitConfidence_SubPiZero =                      new TH1D*[fNBinsPt];
    fHistoBckFitConfidence_FixedPzPiZero =                  new TH1D*[fNBinsPt];
    fFitSignalInvMassPtBin = 								new TF1*[fNBinsPt];
    fFitSignalInvMassPtBin_SubPiZero =						new TF1*[fNBinsPt];
    fFitSignalInvMassPtBin_FixedPzPiZero =   				new TF1*[fNBinsPt];
    fFitSignalInvMassBackFitPtBin = 						new TF1*[fNBinsPt];
    fFitSignalInvMassBackFitPtBin_SubPiZero =				new TF1*[fNBinsPt];
    fFitSignalInvMassBackFitPtBin_FixedPzPiZero =			new TF1*[fNBinsPt];
    fFitSignalPeakPosInvMassLeftPtBin                     = new TF1*[fNBinsPt];
    fFitSignalPeakPosInvMassLeftPtBin_SubPiZero           = new TF1*[fNBinsPt];
    fFitSignalPeakPosInvMassLeftPtBin_FixedPzPiZero       = new TF1*[fNBinsPt];
    fFitTrueSignalInvMassPtBin = 							new TF1*[fNBinsPt];
    fFitTrueSignalInvMassPtReweightedBin =					new TF1*[fNBinsPt];

    fFitBckInvMassPtBin = 									new TF1*[fNBinsPt];
    fFitBckInvMassPtBin_SubPiZero = 						new TF1*[fNBinsPt];
    fFitBckInvMassPtBin_FixedPzPiZero =						new TF1*[fNBinsPt];
    fFitBckInvMassBackFitPtBin =							new TF1*[fNBinsPt];
    fFitBckInvMassBackFitPtBin_SubPiZero =					new TF1*[fNBinsPt];
    fFitBckInvMassBackFitPtBin_FixedPzPiZero =				new TF1*[fNBinsPt];
    fFitRatioInvMassPtBin = 								new TF1*[fNBinsPt];
    // Histograms for normalization on the left of the peak
    fHistoMappingBackNormInvMassLeftPtBin = 				new TH1D*[fNBinsPt];
    fHistoMappingBackNormInvMassLeftPtBin_SubPiZero = 		new TH1D*[fNBinsPt];
    fHistoMappingBackNormInvMassLeftPtBin_FixedPzPiZero = 	new TH1D*[fNBinsPt];
    fHistoMappingSignalInvMassLeftPtBin = 					new TH1D*[fNBinsPt];
    fHistoMappingSignalInvMassLeftPtBin_SubPiZero = 		new TH1D*[fNBinsPt];
    fHistoMappingSignalInvMassLeftPtBin_FixedPzPiZero =		new TH1D*[fNBinsPt];

    fFitInvMassLeftPtBin = 									new TF1*[fNBinsPt];
    fFitInvMassLeftPtBin_SubPiZero =						new TF1*[fNBinsPt];
    fFitInvMassLeftPtBin_FixedPzPiZero =					new TF1*[fNBinsPt];
    fFitBckInvMassLeftPtBin = 								new TF1*[fNBinsPt];
    fFitBckInvMassLeftPtBin_SubPiZero =						new TF1*[fNBinsPt];
    fFitBckInvMassLeftPtBin_FixedPzPiZero =					new TF1*[fNBinsPt];
    fFitWithPol2ForBG = 									new TF1*[fNBinsPt];
    fFitWithPol2ForBG_SubPiZero =							new TF1*[fNBinsPt];
    fFitWithPol2ForBG_FixedPzPiZero =						new TF1*[fNBinsPt];

    for(Int_t i = 0;i<fNBinsPt; i++){
        fHistoMappingTrueMesonInvMassPtBins[i] = 					NULL;
        fHistoMappingTrueMesonInvMassPtReweightedBins[i] =			NULL;// array of histos for pt slices

        fHistoMappingGGInvMassPtBin[i] = 							NULL;
        fHistoMappingGGInvMassPtBin_SubPiZero[i] = 					NULL;
        fHistoMappingGGInvMassPtBin_FixedPzPiZero[i] =				NULL;
        fHistoMappingGGInvMassBackFitPtBin[i] = 					NULL;
        fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i] = 		NULL;
        fBackgroundFitPol[i] = 										NULL;
        fBackgroundFitPol_SubPiZero[i] = 										NULL;
        fBackgroundFitPol_FixedPzPiZero[i] = 										NULL;

        fFitSignalInvMassPtBin[i] = 								NULL;
        fFitSignalInvMassPtBin_SubPiZero[i] =						NULL;
        fFitSignalInvMassPtBin_FixedPzPiZero[i] =					NULL;
        fFitSignalInvMassBackFitPtBin[i] = 							NULL;
        fFitSignalInvMassBackFitPtBin_SubPiZero[i] =				NULL;
        fFitSignalInvMassBackFitPtBin_FixedPzPiZero[i] =			NULL;
        fFitTrueSignalInvMassPtBin[i] = 							NULL;
        fFitTrueSignalInvMassPtReweightedBin[i] = 					NULL;

        fFitBckInvMassPtBin[i] = 									NULL;
        fFitBckInvMassPtBin_SubPiZero[i] =							NULL;
        fFitBckInvMassPtBin_FixedPzPiZero[i] =						NULL;
        fFitBckInvMassBackFitPtBin[i] = 							NULL;
        fFitBckInvMassBackFitPtBin_SubPiZero[i] =					NULL;
        fFitBckInvMassBackFitPtBin_FixedPzPiZero[i] =				NULL;
        fFitRatioInvMassPtBin[i] = 									NULL;
        // Histograms for normalization on the left of the peak
        fHistoMappingBackNormInvMassLeftPtBin[i] = 					NULL;
        fHistoMappingBackNormInvMassLeftPtBin_SubPiZero[i] =		NULL;
        fHistoMappingBackNormInvMassLeftPtBin_FixedPzPiZero[i] =	NULL;
        fHistoMappingSignalInvMassLeftPtBin[i] = 					NULL;
        fHistoMappingSignalInvMassLeftPtBin_SubPiZero[i] =			NULL;
        fHistoMappingSignalInvMassLeftPtBin_FixedPzPiZero[i] =		NULL;

        fFitInvMassLeftPtBin[i] = 									NULL;
        fFitInvMassLeftPtBin_SubPiZero[i] = 									NULL;
        fFitInvMassLeftPtBin_FixedPzPiZero[i] = 									NULL;
        fFitBckInvMassLeftPtBin[i] = 								NULL;
        fFitBckInvMassLeftPtBin_SubPiZero[i] = 								NULL;
        fFitBckInvMassLeftPtBin_FixedPzPiZero[i] = 								NULL;
        fFitWithPol2ForBG[i] = 										NULL;
    }
}

void CalculateFWHM(TF1 * fFunc, Double_t* InputMesonFitRange){
// Default function
    TF1* fFunc_def;
    fFunc_def = new TF1("fFunc_def","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",InputMesonFitRange[0],InputMesonFitRange[1]);
    fFunc_def->SetParameter(0,fFunc->GetParameter(0));
    fFunc_def->SetParameter(1,fFunc->GetParameter(1));
    fFunc_def->SetParameter(2,fFunc->GetParameter(2));
    fFunc_def->SetParameter(3,fFunc->GetParameter(3));



    //FWHM
    fFWHMFunc = fFunc_def->GetX(fFunc_def->GetParameter(0)*0.5,fFunc_def->GetParameter(1), InputMesonFitRange[1]) - fFunc_def->GetX(fFunc_def->GetParameter(0)*0.5,InputMesonFitRange[0],fFunc_def->GetParameter(1));

    //FWHM error +
    TF1* fFunc_plus;
    //	fFunc_plus = fFunc;
    fFunc_plus = new TF1("fFunc_plus","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",InputMesonFitRange[0],InputMesonFitRange[1]);

    fFunc_plus->SetParameter(0,fFunc->GetParameter(0) + fFunc->GetParError(0));
    fFunc_plus->SetParameter(1,fFunc->GetParameter(1) + fFunc->GetParError(1));
    fFunc_plus->SetParameter(2,fFunc->GetParameter(2) + fFunc->GetParError(2));
    fFunc_plus->SetParameter(3,fFunc->GetParameter(3) + fFunc->GetParError(3));
    Double_t FWHM_plus = fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,fFunc_plus->GetParameter(1), InputMesonFitRange[1]) - fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,InputMesonFitRange[0],fFunc_plus->GetParameter(1));

    //FWHM error -
    TF1* fFunc_minus;
    //	fFunc_minus = fFunc;
    fFunc_minus = new TF1("fFunc_minus","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",InputMesonFitRange[0],InputMesonFitRange[1]);
    fFunc_minus->SetParameter(0,fFunc->GetParameter(0) - fFunc->GetParError(0));
    fFunc_minus->SetParameter(1,fFunc->GetParameter(1) - fFunc->GetParError(1));
    fFunc_minus->SetParameter(2,fFunc->GetParameter(2) - fFunc->GetParError(2));
    fFunc_minus->SetParameter(3,fFunc->GetParameter(3) - fFunc->GetParError(3));
    Double_t FWHM_minus =  fFunc_minus->GetX(fFunc_minus->GetParameter(0)*0.5,fFunc_minus->GetParameter(1), InputMesonFitRange[1]) -fFunc_minus->GetX(fFunc_minus->GetParameter(0)*0.5,InputMesonFitRange[0],fFunc_minus->GetParameter(1));

    Double_t Error1 = TMath::Abs(fFWHMFunc-FWHM_plus);
    Double_t Error2 = TMath::Abs(fFWHMFunc-FWHM_minus);

    if(Error1>=Error2) fFWHMFuncError = Error1;
    if(Error1<Error2) fFWHMFuncError = Error2;

}

//Crystal ball function for signal +linear background, parameters are 0:normalization,1:mean,2:sigma,3:n,4:alpha;
Double_t CrystalBallBck(Double_t *x,Double_t *par) {

    Double_t t = (x[0]-par[1])/par[2];
    if (par[4] < 0) t = -t;

    Double_t absAlpha = fabs((Double_t)par[4]);

    if (t >= -absAlpha) {
        return par[0]*exp(-0.5*t*t)+par[5]+par[6]*x[0];
    } else {
        Double_t a =  TMath::Power(par[3]/absAlpha,par[3])*exp(-0.5*absAlpha*absAlpha);
        Double_t b= par[3]/absAlpha - absAlpha;

        return par[0]*(a/TMath::Power(b - t, par[3]))+par[5]+par[6]*x[0];
    }
}

//Crystal ball function for signal, parameters are 0:normalization,1:mean,2:sigma,3:n,4:alpha;
Double_t CrystalBall(Double_t *x,Double_t *par) {

    Double_t t = (x[0]-par[1])/par[2];
    if (par[4] < 0) t = -t;

    Double_t absAlpha = fabs((Double_t)par[4]);

    if (t >= -absAlpha) {
        return par[0]*exp(-0.5*t*t);
    }
    else {
        Double_t a =  TMath::Power(par[3]/absAlpha,par[3])*exp(-0.5*absAlpha*absAlpha);
        Double_t b= par[3]/absAlpha - absAlpha;

        return par[0]*(a/TMath::Power(b - t, par[3]));
    }
}


void Delete(){
    if (fBinsPt) delete fBinsPt;
    if (fPeakRange) delete fPeakRange;
    if (fPeakRange_SubPiZero) delete fPeakRange_SubPiZero;
    if (fPeakRange_FixedPzPiZero) delete fPeakRange_FixedPzPiZero;
    if (fFitRange) delete fFitRange;
    if (fFitRange_SubPiZero) delete fFitRange_SubPiZero;
    if (fFitRange_FixedPzPiZero) delete fFitRange_FixedPzPiZero;
    if (fBGFitRange) delete fBGFitRange;
    if (fBGFitRangeLeft) delete fBGFitRangeLeft;
    if (fBGFitRangeLeft_SubPiZero)     delete fBGFitRangeLeft_SubPiZero;
    if (fBGFitRangeLeft_FixedPzPiZero) delete fBGFitRangeLeft_FixedPzPiZero;
    if (fMesonPlotRange) delete fMesonPlotRange;
    if (fMesonIntRange) delete fMesonIntRange;
    if (fMesonMassRange) delete fMesonMassRange;
    if (fMesonMassRange_SubPiZero) delete fMesonMassRange_SubPiZero;
    if (fMesonMassRange_FixedPzPiZero) delete fMesonMassRange_FixedPzPiZero;
    if (fMesonFitRange) delete fMesonFitRange;
    if (fMesonWidthRange) delete fMesonWidthRange;
    if (fMesonLambdaTailRange) delete fMesonLambdaTailRange;
    if (fNRebin) delete fNRebin;
    for (Int_t k = 0; k< 6; k++){
        if (fGGYields[k])                     delete fGGYields[k];
        if (fGGYields_SubPiZero[k])           delete fGGYields_SubPiZero[k];
        if (fGGYields_FixedPzPiZero[k])       delete fGGYields_FixedPzPiZero[k];
        if (fGGYieldsError[k])                delete fGGYieldsError[k];
        if (fGGYieldsError_SubPiZero[k])                    delete fGGYieldsError_SubPiZero[k];
        if (fGGYieldsError_FixedPzPiZero[k])                delete fGGYieldsError_FixedPzPiZero[k];
        if (fBckYields[k])                        delete fBckYields[k];
        if (fBckYields_SubPiZero[k])              delete fBckYields_SubPiZero[k];
        if (fBckYields_FixedPzPiZero[k])          delete fBckYields_FixedPzPiZero[k];
        if (fBckYieldsError[k])                   delete fBckYieldsError[k];
        if (fBckYieldsError_SubPiZero[k])         delete fBckYieldsError_SubPiZero[k];
        if (fBckYieldsError_FixedPzPiZero[k])     delete fBckYieldsError_FixedPzPiZero[k];
        if (fMesonYields[k])                      delete fMesonYields[k];
        if (fMesonYields_SubPiZero[k])            delete fMesonYields_SubPiZero[k];
        if (fMesonYields_FixedPzPiZero[k])        delete fMesonYields_FixedPzPiZero[k];
        if (fMesonYieldsResidualBckFunc[k])                   delete fMesonYieldsResidualBckFunc[k];
        if (fMesonYieldsResidualBckFunc_SubPiZero[k])         delete fMesonYieldsResidualBckFunc_SubPiZero[k];
        if (fMesonYieldsResidualBckFunc_FixedPzPiZero[k])     delete fMesonYieldsResidualBckFunc_FixedPzPiZero[k];
        if (fMesonYieldsResidualBckFuncError[k])    delete fMesonYieldsResidualBckFuncError[k];
        if (fMesonYieldsResidualBckFuncError_SubPiZero[k])        delete fMesonYieldsResidualBckFuncError_SubPiZero[k];
        if (fMesonYieldsResidualBckFuncError_FixedPzPiZero[k])    delete fMesonYieldsResidualBckFuncError_FixedPzPiZero[k];
        if (fMesonYieldsPerEventError[k])                         delete fMesonYieldsPerEventError[k];
        if (fMesonYieldsPerEventError_SubPiZero[k])               delete fMesonYieldsPerEventError_SubPiZero[k];
        if (fMesonYieldsPerEventError_FixedPzPiZero[k])           delete fMesonYieldsPerEventError_FixedPzPiZero[k];
        if (fMesonYieldsPerEvent[k])                delete fMesonYieldsPerEvent[k];
        if (fMesonYieldsError[k])                   delete fMesonYieldsError[k];
        if (fMesonYieldsError_SubPiZero[k])         delete fMesonYieldsError_SubPiZero[k];
        if (fMesonYieldsError_FixedPzPiZero[k])     delete fMesonYieldsError_FixedPzPiZero[k];
        if (fMesonYieldsCorResidualBckFuncError[k]) delete fMesonYieldsCorResidualBckFuncError[k];
        if (fMesonYieldsCorResidualBckFuncError_SubPiZero[k])     delete fMesonYieldsCorResidualBckFuncError_SubPiZero[k];
        if (fMesonYieldsCorResidualBckFuncError_FixedPzPiZero[k]) delete fMesonYieldsCorResidualBckFuncError_FixedPzPiZero[k];
    }
    for (Int_t k = 0; k< 3; k++){
        if (fMesonTrueYields[k])            delete fMesonTrueYields[k];
        if (fMesonTrueYieldsReweighted[k])  delete fMesonTrueYieldsReweighted[k];
    }


    if (fMesonYieldsBackFit)                              delete fMesonYieldsBackFit;
    if (fMesonYieldsBackFit_SubPiZero)                    delete fMesonYieldsBackFit_SubPiZero;
    if (fMesonYieldsBackFit_FixedPzPiZero)                delete fMesonYieldsBackFit_FixedPzPiZero;
    if (fMesonYieldsResidualBckFuncBackFit)               delete fMesonYieldsResidualBckFuncBackFit;
    if (fMesonYieldsResidualBckFuncBackFit_SubPiZero)     delete fMesonYieldsResidualBckFuncBackFit_SubPiZero;
    if (fMesonYieldsResidualBckFuncBackFit_FixedPzPiZero) delete fMesonYieldsResidualBckFuncBackFit_FixedPzPiZero;

    if (fMesonYieldsCorResidualBckFuncBackFit)                delete fMesonYieldsCorResidualBckFuncBackFit;
    if (fMesonYieldsCorResidualBckFuncBackFit_SubPiZero)      delete fMesonYieldsCorResidualBckFuncBackFit_SubPiZero;
    if (fMesonYieldsCorResidualBckFuncBackFit_FixedPzPiZero)  delete fMesonYieldsCorResidualBckFuncBackFit_FixedPzPiZero;
    if (fMesonYieldsPerEventBackFit)                          delete fMesonYieldsPerEventBackFit;
    if (fMesonYieldsPerEventBackFit_SubPiZero)                delete fMesonYieldsPerEventBackFit_SubPiZero;
    if (fMesonYieldsPerEventBackFit_FixedPzPiZero)            delete fMesonYieldsPerEventBackFit_FixedPzPiZero;
    if (fMesonMass)                             delete fMesonMass;
    if (fMesonMass_SubPiZero)                   delete fMesonMass_SubPiZero;
    if (fMesonMass_FixedPzPiZero)               delete fMesonMass_FixedPzPiZero;
    if (fMesonWidth)                            delete fMesonWidth;
    if (fMesonWidth_SubPiZero)                  delete fMesonWidth_SubPiZero;
    if (fMesonWidth_FixedPzPiZero)              delete fMesonWidth_FixedPzPiZero;
    if (fMesonTrueSB)                           delete fMesonTrueSB;
    if (fMesonTrueSign)                         delete fMesonTrueSign;
    if (fMesonFWHM)                             delete fMesonFWHM;
    if (fMesonFWHM_SubPiZero)                   delete fMesonFWHM_SubPiZero;
    if (fMesonFWHM_FixedPzPiZero)               delete fMesonFWHM_FixedPzPiZero;
    if (fMesonMassLeft)                         delete fMesonMassLeft;
    if (fMesonMassLeft_SubPiZero)               delete fMesonMassLeft_SubPiZero;
    if (fMesonMassLeft_FixedPzPiZero)           delete fMesonMassLeft_FixedPzPiZero;
    if (fMesonWidthLeft)                        delete fMesonWidthLeft;
    if (fMesonWidthLeft_SubPiZero)              delete fMesonWidthLeft_SubPiZero;
    if (fMesonWidthLeft_FixedPzPiZero)          delete fMesonWidthLeft_FixedPzPiZero;
    if (fMesonFWHMLeft)                         delete fMesonFWHMLeft;
    if (fMesonFWHMLeft_SubPiZero)               delete fMesonFWHMLeft_SubPiZero;
    if (fMesonFWHMLeft_FixedPzPiZero)           delete fMesonFWHMLeft_FixedPzPiZero;
    if (fMesonYieldsBackFitError)               delete fMesonYieldsBackFitError;
    if (fMesonYieldsBackFitError_SubPiZero)     delete fMesonYieldsBackFitError_SubPiZero;
    if (fMesonYieldsBackFitError_FixedPzPiZero) delete fMesonYieldsBackFitError_FixedPzPiZero;


    if (fMesonYieldsResidualBckFuncBackFitError)    delete fMesonYieldsResidualBckFuncBackFitError;
    if (fMesonYieldsCorResidualBckFuncBackFitError) delete fMesonYieldsCorResidualBckFuncBackFitError;
    if (fMesonYieldsCorResidualBckFuncBackFitError_SubPiZero)     delete fMesonYieldsCorResidualBckFuncBackFitError_SubPiZero;
    if (fMesonYieldsCorResidualBckFuncBackFitError_FixedPzPiZero) delete fMesonYieldsCorResidualBckFuncBackFitError_FixedPzPiZero;
    if (fMesonYieldsPerEventBackFitError)       delete fMesonYieldsPerEventBackFitError;
    if (fMesonYieldsPerEventBackFitError_SubPiZero)       delete fMesonYieldsPerEventBackFitError_SubPiZero;
    if (fMesonYieldsPerEventBackFitError_FixedPzPiZero)   delete fMesonYieldsPerEventBackFitError_FixedPzPiZero;
    if (fMesonMassError)                        delete fMesonMassError;
    if (fMesonMassError_SubPiZero)              delete fMesonMassError_SubPiZero;
    if (fMesonMassError_FixedPzPiZero)          delete fMesonMassError_FixedPzPiZero;
    if (fMesonWidthError)                       delete fMesonWidthError;
    if (fMesonWidthError_SubPiZero)             delete fMesonWidthError_SubPiZero;
    if (fMesonWidthError_FixedPzPiZero)         delete fMesonWidthError_FixedPzPiZero;
    if (fMesonTrueSBError)                      delete fMesonTrueSBError;
    if (fMesonTrueSignError)                    delete fMesonTrueSignError;
    if (fMesonFWHMError)                        delete fMesonFWHMError;
    if (fMesonFWHMError_SubPiZero)              delete fMesonFWHMError_SubPiZero;
    if (fMesonFWHMError_FixedPzPiZero)          delete fMesonFWHMError_FixedPzPiZero;
    if (fMesonFWHMLeftError)                    delete fMesonFWHMLeftError;
    if (fMesonFWHMLeftError_SubPiZero)          delete fMesonFWHMLeftError_SubPiZero;
    if (fMesonFWHMLeftError_FixedPzPiZero)      delete fMesonFWHMLeftError_FixedPzPiZero;
    if (fHistoMappingTrueMesonInvMassPtBins)    delete fHistoMappingTrueMesonInvMassPtBins;
    if (fHistoMappingTrueMesonInvMassPtReweightedBins) delete fHistoMappingTrueMesonInvMassPtReweightedBins;
    if (fHistoMappingGGInvMassPtBin)                   delete fHistoMappingGGInvMassPtBin;
    if (fHistoMappingGGInvMassPtBin_SubPiZero)         delete fHistoMappingGGInvMassPtBin_SubPiZero;
    if (fHistoMappingGGInvMassPtBin_FixedPzPiZero)     delete fHistoMappingGGInvMassPtBin_FixedPzPiZero;
    for (Int_t k = 0; k< 5; k++){
        if (fHistoMappingBackInvMassPtBin[k])               delete fHistoMappingBackInvMassPtBin[k];
        if (fHistoMappingBackInvMassPtBin_SubPiZero[k])     delete fHistoMappingBackInvMassPtBin_SubPiZero[k];
        if (fHistoMappingBackInvMassPtBin_FixedPzPiZero[k]) delete fHistoMappingBackInvMassPtBin_FixedPzPiZero[k];
        if (fHistoMappingBackSameNormInvMassPtBin[k]) delete fHistoMappingBackSameNormInvMassPtBin[k];
        if (fHistoMappingBackSameNormInvMassPtBin_SubPiZero[k]) delete fHistoMappingBackSameNormInvMassPtBin_SubPiZero[k];
        if (fHistoMappingBackSameNormInvMassPtBin_FixedPzPiZero[k]) delete fHistoMappingBackSameNormInvMassPtBin_FixedPzPiZero[k];
        if (fHistoMappingBackNormInvMassPtBin[k]) delete fHistoMappingBackNormInvMassPtBin[k];
        if (fHistoMappingBackNormInvMassPtBin_SubPiZero[k]) delete fHistoMappingBackNormInvMassPtBin_SubPiZero[k];
        if (fHistoMappingBackNormInvMassPtBin_FixedPzPiZero[k]) delete fHistoMappingBackNormInvMassPtBin_FixedPzPiZero[k];
        if (fHistoMappingRatioSBInvMassPtBin[k]) delete fHistoMappingRatioSBInvMassPtBin[k];
        if (fHistoMappingRatioSBInvMassPtBin_SubPiZero[k]) delete fHistoMappingRatioSBInvMassPtBin_SubPiZero[k];
        if (fHistoMappingRatioSBInvMassPtBin_FixedPzPiZero[k]) delete fHistoMappingRatioSBInvMassPtBin_FixedPzPiZero[k];
    }
    if (fFitSignalInvMassPtBin)               delete fFitSignalInvMassPtBin;
    if (fFitSignalInvMassPtBin_SubPiZero)     delete fFitSignalInvMassPtBin_SubPiZero;
    if (fFitSignalInvMassPtBin_FixedPzPiZero) delete fFitSignalInvMassPtBin_FixedPzPiZero;
    if (fFitBckInvMassPtBin)                  delete fFitBckInvMassPtBin;
    if (fFitBckInvMassPtBin_SubPiZero)                       delete fFitBckInvMassPtBin_SubPiZero;
    if (fFitBckInvMassPtBin_FixedPzPiZero)                   delete fFitBckInvMassPtBin_FixedPzPiZero;
    if (fHistoMappingBackNormInvMassLeftPtBin)               delete fHistoMappingBackNormInvMassLeftPtBin;
    if (fHistoMappingBackNormInvMassLeftPtBin_SubPiZero)     delete fHistoMappingBackNormInvMassLeftPtBin_SubPiZero;
    if (fHistoMappingBackNormInvMassLeftPtBin_FixedPzPiZero) delete fHistoMappingBackNormInvMassLeftPtBin_FixedPzPiZero;
    if (fHistoMappingSignalInvMassLeftPtBin)               delete fHistoMappingSignalInvMassLeftPtBin;
    if (fHistoMappingSignalInvMassLeftPtBin_SubPiZero)     delete fHistoMappingSignalInvMassLeftPtBin_SubPiZero;
    if (fHistoMappingSignalInvMassLeftPtBin_FixedPzPiZero) delete fHistoMappingSignalInvMassLeftPtBin_FixedPzPiZero;
    if (fFitInvMassLeftPtBin)               delete fFitInvMassLeftPtBin;
    if (fFitInvMassLeftPtBin_SubPiZero)     delete fFitInvMassLeftPtBin_SubPiZero;
    if (fFitInvMassLeftPtBin_FixedPzPiZero) delete fFitInvMassLeftPtBin_FixedPzPiZero;
    if (fFitSignalPeakPosInvMassLeftPtBin)                delete fFitSignalPeakPosInvMassLeftPtBin;
    if (fFitSignalPeakPosInvMassLeftPtBin_SubPiZero)      delete fFitSignalPeakPosInvMassLeftPtBin_SubPiZero;
    if (fFitSignalPeakPosInvMassLeftPtBin_FixedPzPiZero)  delete fFitSignalPeakPosInvMassLeftPtBin_FixedPzPiZero;
    if (fFitBckInvMassLeftPtBin)                          delete fFitBckInvMassLeftPtBin;
    if (fFitBckInvMassLeftPtBin_SubPiZero)                delete fFitBckInvMassLeftPtBin_SubPiZero;
    if (fFitBckInvMassLeftPtBin_FixedPzPiZero)            delete fFitBckInvMassLeftPtBin_FixedPzPiZero;
    if (fFitRatioInvMassPtBin) delete fFitRatioInvMassPtBin;
    if (fFitWithPol2ForBG)                  delete fFitWithPol2ForBG;
    if (fFitWithPol2ForBG_SubPiZero)        delete fFitWithPol2ForBG_SubPiZero;
    if (fFitWithPol2ForBG_FixedPzPiZero)    delete fFitWithPol2ForBG_FixedPzPiZero;
    if (fHistoWeightsBGZbinVsMbin)          delete fHistoWeightsBGZbinVsMbin;
    if (fHistoFillPerEventBGZbinVsMbin)     delete fHistoFillPerEventBGZbinVsMbin;
}


void SetCorrectMCHistogrammNames(){
    cout << "standard MC chosen" << endl;
    ObjectNameTrue 						= "ESD_TrueMotherPiPlPiMiPiZero_InvMass_Pt";
    ObjectNameTrueWOWeights 			= "ESD_TrueMotherPiPlPiMiPiZeroWOWeights_InvMass_Pt";	//  n/a
    ObjectNameProfileWeights 			= "ESD_TruePrimaryMotherWeights_InvMass_Pt";			//  n/a
    ObjectNameMCEtaAcc 					= "MC_EtaInAcc_Pt";
    ObjectNameMCOmegaAcc 				= "MC_OmegaInAcc_Pt";
    ObjectNameMCEta 					= "MC_Eta_Pt";
    ObjectNameMCEtaWOWeights 			= "MC_Eta_WOWeights_Pt";								//  n/a
    ObjectNameMCOmega 					= "MC_Omega_Pt";
    ObjectNameMCOmegaWOWeights 			= "MC_Omega_WOWeights_Pt";								//  n/a
}

