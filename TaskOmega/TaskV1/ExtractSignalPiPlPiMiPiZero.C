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
#include "../../CommonHeaders/PlottingMeson.h"
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
#include "../../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "ExtractSignalPiPlPiMiPiZero.h"
#include "../../CommonHeaders/ExtractSignalBinning.h"
#include "../../CommonHeaders/FittingGammaConversion.h"
#include "../../CommonHeaders/ConversionFunctions.h" // changed to standard CommonHeaders
#include "../../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../../CommonHeaders/ExtractSignalPlotting.h"
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
    gSystem->Exec("mkdir -p "+outputDir+"/Plots_Data_QA");
    gSystem->Exec("mkdir -p "+outputDir+"/Plots_MC_QA");

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
    TString ObjectNameBck[3]        = { "ESD_Background_2_InvMass_Pt",
                                        "ESD_Background_3_InvMass_Pt",
                                        "ESD_Background_4_InvMass_Pt"
                                      };

    hist_bck2                       = (TH2D*)ESDContainer->FindObject(ObjectNameBck[0].Data());
    hist_bck3                       = (TH2D*)ESDContainer->FindObject(ObjectNameBck[1].Data());
    hist_bck4                       = (TH2D*)ESDContainer->FindObject(ObjectNameBck[2].Data());

    fBckInvMassVSPt2                = (TH2D*)hist_bck2->Clone("hist_bck2");
    fBckInvMassVSPt3                = (TH2D*)hist_bck3->Clone("hist_bck3");
    fBckInvMassVSPt4                = (TH2D*)hist_bck3->Clone("hist_bck4");

    SetCorrectMCHistogrammNames();

    fGammaGammaInvMassVSPt= (TH2D*)ESDContainer->FindObject(ObjectNameESD.Data());

    fBckInvMassVSPt                 = fBckInvMassVSPt2;
    fBckInvMassVSPt->Add(fBckInvMassVSPt3);
    fBckInvMassVSPt->Add(fBckInvMassVSPt4);

    const char* FileDataLogname     = Form("%s/%s/%s_%s_EffiCheck_RAWDATA%s_%s.dat", cutSelection.Data(), fEnergyFlag.Data(), fPrefix.Data(), fPrefix2.Data(), fPeriodFlag.Data(), fCutSelectionRead.Data());
    fFileDataLog.open(FileDataLogname, ios::out);

    ProduceBckWithoutWeighting(fBckInvMassVSPt); //background without weighting because THSparse wasn't used


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

    TH1D *fBck              = (TH1D*)fBckInvMassVSPt->ProjectionX("");
    TH1D *fGammaGamma       = (TH1D*)fGammaGammaInvMassVSPt->ProjectionX("ESD_Mother_InvMass");

    cout<< "The mass of the meson is: "<< fMesonMassExpect<< " Events analysed: "<< fNEvents<< endl;
    cout << "here" << endl;
    // Process the 1D invariant mass histos
    fGammaGamma->SetTitle(Form("%s %s",fGammaGamma->GetTitle(),fCutSelection.Data()));
    cout << "here" << endl;
    fGammaGamma->Rebin(fNRebinGlobal);
    //fGammaGamma->Scale(1./fNRebinGlobal);
    fBck->Rebin(fNRebinGlobal);
    //fBck->Scale(1./fNRebinGlobal);
    ProcessEM( fGammaGamma , fBck, fBGFitRange);
    fHistoMappingBackNormInvMass = fBckNorm;
    fHistoMappingSignalInvMass = fSignal;

    fGammaGamma->DrawCopy();
    fHistoMappingBackNormInvMass->DrawCopy("same");
    fHistoMappingSignalInvMass->DrawCopy("same");

    // Function to Project the 2D histos InvariantMass VS Pt into Invariant Mass spectrum
    FillMassHistosArray(fGammaGammaInvMassVSPt);
    ProcessEM( fMesonFullPtSignal, fMesonFullPtBackground, fBGFitRange);
    fMesonFullPtBackNorm            = fBckNorm;

    ProcessEM( fFittingHistMidPtSignal, fFittingHistMidPtBackground, fBGFitRange);
    fFittingHistMidPtSignalSub = fSignal;
    if(fCrysFitting==0){
        fFileErrLog << "Using exp fit"<<endl;
        FitSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub, fMesonIntDeltaRange,200,kTRUE);
        fFitSignalInvMassMidPt = fFitReco;
    } else {
        fFileErrLog << "Using Crystal Ball function"<<endl;
        FitCBSubtractedInvMassInPtBins(fFittingHistMidPtSignalSub, fMesonIntDeltaRange,200,kTRUE,"SinglefitfunctionMidPt");
        fFitSignalInvMassMidPt = fFitReco;
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
    fOmegaFullPtSignal = (TH1D*)fMesonFullPtSignal->Clone("fMesonFullPtSignal");
    fOmegaFullPtSignal->Write("fSignal");
    // Different event mixing background groups
    fOmegaFullPtBack2 = (TH1D*)hist_bck2->Clone("hist_bck2");
    fOmegaFullPtBack2->Write("fBack2");
    fOmegaFullPtBack3 = (TH1D*)hist_bck3->Clone("hist_bck3");
    fOmegaFullPtBack3->Write("fBack3");
    fOmegaFullPtBack4 = (TH1D*)hist_bck4->Clone("hist_bck4");
    fOmegaFullPtBack4->Write("fBack4");
    fOmegaFullPtBackNorm = (TH1D*)fMesonFullPtBackNorm->Clone("fMesonFullPtBackNorm");
    fOmegaFullPtBackNorm->Write("fBackNorm");
    filefull->Close();


    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){ // BEGIN ANALYSIS for each Pt bin

        cout << "---------------------------------------------------------------------------------" << endl;
        cout << "Begin Analysis Pt Bin " << iPt <<endl;
        cout << "---------------------------------------------------------------------------------" << endl;
        // Function to subtract GG minus Bck
        fFileDataLog << "---------------------------------------------------------------------------------" << endl;
        fFileDataLog << "----------------------------------new pT bin ------------------------------------" << endl;
        fFileDataLog << "---------------------------------------------------------------------------------" << endl;

        ProcessBckFitSubtraction(fHistoMappingGGInvMassPtBin[iPt],iPt,fPeakRange,fFitRange,optionEnergy,Suffix,cutSelection,meson);

        ProcessEM( fHistoMappingGGInvMassPtBin[iPt], fHistoMappingBackInvMassPtBin[iPt], fBGFitRange);
        fHistoMappingSignalInvMassPtBin[iPt] = fSignal;
        fHistoMappingBackNormInvMassPtBin[iPt] = fBckNorm;

        TString namesecHistoBckNorm;
        TString namesecHistoSignal;
        TString namesecRatio;

        TString tempStr = (iPt==fStartPtBin)?"recreate":"update";

        TFile *file1 = TFile::Open("histback.root",tempStr.Data());
        //writes background  and signal for each pt bin in a root file - carolina's modification
        file1->cd();
        namesecHistoBckNorm                         = Form("BckNorm_InvMass_in_Pt_Bin%02i", iPt);
        fHistoMappingBackNormInvMassPtBin[iPt]->Write(namesecHistoBckNorm.Data());
        namesecHistoSignal                          = Form("Signal_InvMass_in_Pt_Bin%02i", iPt);
        fHistoMappingGGInvMassPtBin[iPt]->Write(namesecHistoSignal.Data());

        file1->Close();

        // TODO how do you use this function without using fMesonFitRange
        FitWithPol2ForBG(fHistoMappingGGInvMassPtBin[iPt], fMesonFitRange,iPt,kFALSE);
        fFitWithPol2ForBG[iPt]                       = fFitReco;

        ProcessRatioSignalBackground(fHistoMappingGGInvMassPtBin[iPt], fHistoMappingBackNormInvMassPtBin[iPt]);
        fHistoMappingRatioSBInvMassPtBin[iPt]        = fRatioSB;

        fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "normal range/right normalization" << endl;

        fFitSignalInvMassPtBin[iPt]                 = 0x00;
        fFitSignalInvMassBackFitPtBin[iPt]          = 0x00;
        if(fCrysFitting==0){
            fFileErrLog << "Using exp fit"<<endl;
            fFileDataLog << "Subtracted mixed event" << endl;
            FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
            fFitSignalInvMassPtBin[iPt]                     = fFitReco;
            fFitBckInvMassPtBin[iPt]                        = fFitLinearBck;
            fMesonYieldsResidualBckFunc[0][iPt]             = fIntLinearBck;
            fMesonYieldsResidualBckFuncError[0][iPt]        = fIntLinearBckError;
            fFileDataLog << "Subtracted polinomial fit" << endl;
            FitSubtractedInvMassInPtBins(fHistoMappingGGInvMassBackFitPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
            fFitSignalInvMassBackFitPtBin[iPt]              = fFitReco;
            fFitBckInvMassBackFitPtBin[iPt]                 = fFitLinearBck;
            fMesonYieldsResidualBckFuncBackFit[iPt]         = fIntLinearBck;
            fMesonYieldsResidualBckFuncBackFitError[iPt]    = fIntLinearBckError;

        } else {
            fFileErrLog << "Using Crystal Ball function"<<endl;
            FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE,Form("CBFitFuncNormalBin%02d",iPt));
            fFitSignalInvMassPtBin[iPt]                     = fFitReco;
            fFitBckInvMassPtBin[iPt]                        = fFitLinearBck;
            fMesonYieldsResidualBckFunc[0][iPt]             = fIntLinearBck;
            fMesonYieldsResidualBckFuncError[0][iPt]        = fIntLinearBckError;
        }

        //GetFWHM
        CalculateFWHM( fFitSignalInvMassPtBin[iPt]);
        fMesonFWHM[iPt]                                     = fFWHMFunc;
        fMesonFWHMError[iPt]                                = fFWHMFuncError;

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
        } else {
            fMesonMassBackFit[iPt]          = 0;
            fMesonMassBackFitError[iPt]     = 0;
            fMesonWidthBackFit[iPt]         = 0;
            fMesonWidthBackFitError[iPt]    = 0;
            //MesonCurIntRangeBackFit[0]     = fMesonIntRange[0];
            //MesonCurIntRangeBackFit[1]     = fMesonIntRange[1];
            fMesonCurIntRangeBackFit[0]     = fMesonMass[iPt] + fMesonIntDeltaRange[0]; // TODO look into this (for now same as fMesonCurIntRange)
            fMesonCurIntRangeBackFit[1]     = fMesonMass[iPt] + fMesonIntDeltaRange[1];
        }

        for (Int_t k = 0; k < 3; k++){
            IntegrateHistoInvMass( fHistoMappingGGInvMassPtBin[iPt], fMesonCurIntRange[k]);
            fGGYields[k][iPt]               = fYields;
            fGGYieldsError[k][iPt]          = fYieldsError;

            // Integrate the bck histo
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassPtBin[iPt], fMesonCurIntRange[k]);
            fBckYields[k][iPt]              = fYields;
            fBckYieldsError[k][iPt]         = fYieldsError;

            if (k == 0){
                // Integrate the Backfit histo
                fFileDataLog<< endl <<"BackFit histo "<< nameIntRange[k].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
                IntegrateHistoInvMassStream( fHistoMappingGGInvMassBackFitPtBin[iPt], fMesonCurIntRangeBackFit); // fMesonCurIntRangeBackFit no array ?
                fMesonYieldsBackFit[iPt]      = fYields; // change this part because it is written each time of the loop -> TODO Change to array or move outside of loop
                fMesonYieldsBackFitError[iPt] = fYieldsError;
            }
            // Integrate the signal histo
            fFileDataLog<< endl <<"Signal histo "<< nameIntRange[k].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalInvMassPtBin[iPt], fMesonCurIntRange[k]); //cout<<"2. call of IntegrateHistoInvMassStream"<<endl;
            fMesonYields[k][iPt]            = fYields;
            fMesonYieldsError[k][iPt]       = fYieldsError;
            fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;


            if(fIsMC){
                fFileDataLog<< endl <<"True histo normal range" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                fFitTrueSignalInvMassPtBin[iPt]=0x00;
                if(fCrysFitting==0){
                    fFileErrLog << "Using exp fit"<<endl;
                    FitTrueInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntRange,iPt,kFALSE);
                } else {
                    fFileErrLog << "Using Crystal Ball function"<<endl;
                    FitCBSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins[iPt], fMesonIntRange,iPt,kFALSE,Form("CBFitFuncMCTrueBin%02d",iPt));
                }

                if (fHistoMappingTrueMesonInvMassPtBins[iPt]->GetEntries() !=0){
                    fFitTrueSignalInvMassPtBin[iPt]                  = fFitReco;
                    if (fFitTrueSignalInvMassPtBin[iPt] != 0x00){
                        fMesonTrueMass[iPt]                          = fFitTrueSignalInvMassPtBin[iPt]->GetParameter(1);
                        fMesonTrueMassError[iPt]                     = fFitTrueSignalInvMassPtBin[iPt]->GetParError(1);
                        CalculateFWHM(fFitTrueSignalInvMassPtBin[iPt]);
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
            fTotalBckYieldsError[k][iPt]                    = pow(fBckYieldsError[k][iPt]*fBckYieldsError[k][iPt] + fMesonYieldsResidualBckFuncError[k][iPt]*fMesonYieldsResidualBckFuncError[k][iPt],0.5);
            fFileDataLog << "Total Background " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fTotalBckYields[k][iPt] << "\t +- \t" << fTotalBckYieldsError[k][iPt] << endl<< endl;
            fFileDataLog << "Background " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fBckYields[k][iPt] << "\t +- \t" << fBckYieldsError[k][iPt] << endl<< endl;

            fMesonYieldsCorResidualBckFunc[k][iPt]          = fMesonYields[k][iPt]- fMesonYieldsResidualBckFunc[k][iPt];
            fMesonYieldsCorResidualBckFuncError[k][iPt]     = pow(( fMesonYieldsError[k][iPt]*fMesonYieldsError[k][iPt] +
                                                                fMesonYieldsResidualBckFuncError[k][iPt]*fMesonYieldsResidualBckFuncError[k][iPt]),0.5);
            fMesonYieldsPerEvent[k][iPt]                    = fMesonYieldsCorResidualBckFunc[k][iPt]/fNEvents;
            fMesonYieldsPerEventError[k][iPt]               = fMesonYieldsCorResidualBckFuncError[k][iPt]/fNEvents;

            // check this part, maybe move to an array later on. Dirty fix -> only do it once for k == 0)
            if (k == 0){
                fMesonYieldsCorResidualBckFuncBackFit[iPt]          = fMesonYieldsBackFit[iPt]- fMesonYieldsResidualBckFuncBackFit[iPt];
                fMesonYieldsCorResidualBckFuncBackFitError[iPt]     =
                pow((fMesonYieldsBackFitError[iPt]*fMesonYieldsBackFitError[iPt]+
                    fMesonYieldsResidualBckFuncBackFitError[iPt]*fMesonYieldsResidualBckFuncBackFitError[iPt]),0.5);
                fMesonYieldsPerEventBackFit[iPt]                    = fMesonYieldsCorResidualBckFuncBackFit[iPt]/fNEvents;
                fMesonYieldsPerEventBackFitError[iPt]               = fMesonYieldsCorResidualBckFuncBackFitError[iPt]/fNEvents;
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
                    FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassPtBin[iPt], intRange, iPt, kFALSE);
                    fMesonYieldsResidualBckFunc[k][iPt]         = fIntLinearBck;
                    fMesonYieldsResidualBckFuncError[k+3][iPt]    = fIntLinearBckError;
                } else {
                    fMesonYieldsResidualBckFunc[k][iPt]         = 0;
                    fMesonYieldsResidualBckFuncError[k][iPt]    = 0;
                }
            }
            fFileDataLog << "Residual Background leftover " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                            << fMesonYieldsResidualBckFunc[k][iPt] << "\t +- \t" << fMesonYieldsResidualBckFuncError[k][iPt] << endl<< endl;


            fTotalBckYields[k][iPt]                         = fBckYields[k][iPt] + fMesonYieldsResidualBckFunc[k][iPt];
            fTotalBckYieldsError[k][iPt]                    = pow(fBckYieldsError[k][iPt]*fBckYieldsError[k][iPt] + fMesonYieldsResidualBckFuncError[k][iPt]*fMesonYieldsResidualBckFuncError[k][iPt],0.5);
            fFileDataLog << "Total Background " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fTotalBckYields[k][iPt] << "\t +- \t" << fTotalBckYieldsError[k][iPt] << endl<< endl;
            fFileDataLog << "Background " << nameIntRange[k].Data() << " in iPt " << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << ":\t"
                        << fBckYields[k][iPt] << "\t +- \t" << fBckYieldsError[k][iPt] << endl<< endl;

            fMesonYieldsCorResidualBckFunc[k][iPt]          = fMesonYields[k][iPt]- fMesonYieldsResidualBckFunc[k][iPt];
            fMesonYieldsCorResidualBckFuncError[k][iPt]     = pow(( fMesonYieldsError[k][iPt]*fMesonYieldsError[k][iPt] +
                                                                fMesonYieldsResidualBckFuncError[k][iPt]*fMesonYieldsResidualBckFuncError[k][iPt]),0.5);
            fMesonYieldsPerEvent[k][iPt]                    = fMesonYieldsCorResidualBckFunc[k][iPt]/fNEvents;
            fMesonYieldsPerEventError[k][iPt]               = fMesonYieldsCorResidualBckFuncError[k][iPt]/fNEvents;
        }


        //////////////////////////////// Start Analysis with  Normalization at the left of the Meson Peak

        // Function to subtract GG minus Bck
        cout << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
        ProcessEMLeftRight( fHistoMappingGGInvMassPtBin[iPt], fHistoMappingBackInvMassPtBin[iPt], fBGFitRangeLeft, fBGFitRange);
        fHistoMappingSignalInvMassLeftPtBin[iPt] = fSignal;
        fHistoMappingBackNormInvMassLeftPtBin[iPt] = fBckNorm;

        // Fitting the subtracted spectra
        fFileErrLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "normal range/left normalization" << endl;

        fFitInvMassLeftPtBin[iPt] =0x00;
        if(fCrysFitting==0){
            fFileErrLog << "Using exp fit"<<endl;
            FitSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
            fFitInvMassLeftPtBin[iPt]                   = fFitReco;
            fFitSignalPeakPosInvMassLeftPtBin[iPt]      = fFitGausExp;
            fFitBckInvMassLeftPtBin[iPt]                = fFitLinearBck;
            fMesonYieldsResidualBckFunc[3][iPt]         = fIntLinearBck;
            fMesonYieldsResidualBckFuncError[3][iPt]    = fIntLinearBckError;

        } else {
            fFileErrLog << "Using Crystal Ball function"<<endl;
            FitCBSubtractedInvMassInPtBins(fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonIntDeltaRange,iPt,kFALSE, Form("CBFitFuncLeftBin%02d",iPt));
            fFitInvMassLeftPtBin[iPt]                   = fFitReco;
            fFitSignalPeakPosInvMassLeftPtBin[iPt]      = fFitGausExp;
            fFitBckInvMassLeftPtBin[iPt]                = fFitLinearBck;
            fMesonYieldsResidualBckFunc[3][iPt]         = fIntLinearBck;
            fMesonYieldsResidualBckFuncError[3][iPt]    = fIntLinearBckError;
        }

        CalculateFWHM(fFitInvMassLeftPtBin[iPt]);
        fMesonFWHMLeft[iPt] = fFWHMFunc;
        fMesonFWHMLeftError[iPt] = fFWHMFuncError;

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
        }

        // Integrate the bck histo
        for (Int_t k = 0; k < 3; k++){
            IntegrateHistoInvMass( fHistoMappingBackNormInvMassLeftPtBin[iPt], fMesonCurIntRange[k+3]);
            fBckYields[k+3][iPt]              = fYields;
            fBckYieldsError[k+3][iPt]         = fYieldsError;

            fFileDataLog<< endl <<"Signal histo " << nameIntRange[k+3].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
            IntegrateHistoInvMassStream( fHistoMappingSignalInvMassLeftPtBin[iPt], fMesonCurIntRange[k+3]);
            fMesonYields[k+3][iPt]            = fYields;
            fMesonYieldsError[k+3][iPt]       = fYieldsError;
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
            fMesonYieldsCorResidualBckFuncError[k+3][iPt]   = pow(( fMesonYieldsError[k+3][iPt]*fMesonYieldsError[k+3][iPt] +
                                                                    fMesonYieldsResidualBckFuncError[k+3][iPt]*fMesonYieldsResidualBckFuncError[k+3][iPt]),0.5);
            fMesonYieldsPerEvent[k+3][iPt]                  = fMesonYieldsCorResidualBckFunc[k+3][iPt]/fNEvents;
            fMesonYieldsPerEventError[k+3][iPt]             = fMesonYieldsCorResidualBckFuncError[k+3][iPt]/fNEvents;

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

    TString nameMeson           = Form("%s/%s_%s_MesonWithBck%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    TString nameCanvas          = "MesonWithBckCanvas";
    TString namePad             = "MesonWithBckPad";
    PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingBackNormInvMassPtBin, nameMeson, nameCanvas, namePad, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcess, fCollisionSystem);

    TString nameMesonSub        = Form("%s/%s_%s_MesonSubtracted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
    TString nameCanvasSub       = "MesonCanvasSubtracted";
    TString namePadSub          = "MesonPadSubtracted";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem);

    nameMesonSub                = Form("%s/%s_%s_MesonSubtractedBackFit%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasSubtractedBackFit";
    namePadSub                  = "MesonPadSubtractedBackFit";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassBackFitPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitSignalInvMassBackFitPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem);
    //loop through Pt slices and plot inv mass around peak (only in DATA, not MC)

    nameMesonSub                = Form("%s/%s_%s_MesonBackFit%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(), fPeriodFlag.Data(), fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasBackFit";
    namePadSub                  = "MesonPadBackFit";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins, fBackgroundFitPol, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem);

    nameMeson                   = Form("%s/%s_%s_MesonWithBckLeft%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    nameCanvas                  = "MesonWithBckCanvasLeft";
    namePad                     = "MesonWithBckPadLeft";
    PlotInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingBackNormInvMassLeftPtBin, nameMeson, nameCanvas, namePad,  fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcess, fCollisionSystem);

    nameMesonSub                = Form("%s/%s_%s_MesonWithBGAndFit%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasWithBGAndFit";
    namePadSub                  = "MesonPadWithBGAndFit";
    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingGGInvMassPtBin, fHistoMappingTrueMesonInvMassPtBins,fFitWithPol2ForBG, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC , fDecayChannel, fDetectionProcess, fCollisionSystem);

    nameMesonSub                = Form("%s/%s_%s_MesonSubtractedLeft%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
    nameCanvasSub               = "MesonCanvasSubtractedLeft";
    namePadSub                  = "MesonPadSubtractedLeft";

    PlotWithFitSubtractedInvMassInPtBins( fHistoMappingSignalInvMassLeftPtBin, fHistoMappingTrueMesonInvMassPtBins, fFitInvMassLeftPtBin, nameMesonSub, nameCanvasSub, namePadSub, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem);

    PlotExampleInvMassBinsV2(fHistoMappingGGInvMassPtBin[fExampleBin], fHistoMappingSignalInvMassPtBin[fExampleBin], fHistoMappingBackNormInvMassPtBin[fExampleBin],
                        fFitSignalInvMassPtBin[fExampleBin], fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                        fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, 0, fScaleFac, fMode, addSig );
    PlotExampleInvMassBinsBckFit(fHistoMappingGGInvMassPtBin[fExampleBin], fHistoMappingGGInvMassBackFitPtBin[fExampleBin], fBackgroundFitPol[fExampleBin],
                         fFitSignalInvMassBackFitPtBin[fExampleBin], fExampleBin, outputDir.Data(),Suffix.Data(), fMesonMassRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                        fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, 0, fScaleFac, fMode, addSig );
    if(fIsMC){
        TString nameMesonTrue   = Form("%s/%s_%s_TrueMesonFitted%s_%s.%s",outputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fPeriodFlag.Data(),fCutSelection.Data(),Suffix.Data());
        TString nameCanvasTrue  = "TrueMesonCanvasFitted";
        TString namePadTrue     = "TrueMesonPadFitted";
        PlotWithFitSubtractedInvMassInPtBins(fHistoMappingTrueMesonInvMassPtBins, fHistoMappingTrueMesonInvMassPtBins, fFitTrueSignalInvMassPtBin, nameMesonTrue, nameCanvasTrue, namePadTrue, fMesonMassRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess, fCollisionSystem);
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

void ProcessBckFitSubtraction(TH1D *fGammaGamma, Int_t i, Double_t * fPeakRangeDummy, Double_t *fFitRangeDummy, TString energy, TString suffix, TString cutSelection, TString meson)
{

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
    fBackgroundFitPol[i] = NULL;

    //2ND OPTION
    FitBackFunc = new TF1("BGfit",FunctionBGExclusion,fFitRangeDummy[0],fFitRangeDummy[1],5);
    fBackgroundFitPol[i] = new TF1("BGfit","pol4",fFitRangeDummy[0],fFitRangeDummy[1]);
    fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->Fit(FitBackFunc,"QMRE0","",fFitRangeDummy[0],fFitRangeDummy[1]);//QMRE0

    Double_t FitParams[5];
    FitBackFunc->GetParameters(&FitParams[0]);
    fBackgroundFitPol[i]->SetParameters(FitParams);

    fBackgroundFitPol[i]->SetRange(fMesonFitRange[0],fMesonFitRange[1]);
    /*delete FitFuncLow;
    delete FitFuncHigh;
        */
    for (Int_t binx= 0; binx < fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetNbinsX()+1; binx++){
        if(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) > fFitRangeDummy[0] && fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx) < fFitRangeDummy[1]){
            fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinContent(binx,fBackgroundFitPol[i]->Eval(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->GetBinCenter(binx)));
            fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i]->SetBinError(binx,fGammaGamma->GetBinError(binx));
        }
    }
    fHistoMappingGGInvMassBackFitPtBin[i] = (TH1D*) fGammaGamma->Clone(Form("GG_SubtractedSignal_%i",i));
    fHistoMappingGGInvMassBackFitPtBin[i]->Sumw2();

    fHistoMappingGGInvMassBackFitPtBin[i]->Add(fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i],-1);

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
    Double_t 	norm = 1;

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

void ProduceBckWithoutWeighting(TH2D *fBckInvMassVSPtDummy){
    //calculation background for midPt without weighting
    Int_t startBinMidPt = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fMidPt[0]+0.001);
    Int_t endBinMidPt = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fMidPt[1]-0.001);
    fFittingHistMidPtBackground = new TH1D("Mapping_Back_InvMass_MidPt","Mapping_Back_InvMass_MidPt",fBckInvMassVSPtDummy->GetNbinsX(),0.,1.);
    fFittingHistMidPtBackground->Sumw2();
    fBckInvMassVSPtDummy->ProjectionX("Mapping_Back_InvMass_MidPt",startBinMidPt,endBinMidPt);
    fFittingHistMidPtBackground=(TH1D*)gDirectory->Get("Mapping_Back_InvMass_MidPt");
    fFittingHistMidPtBackground->Rebin(fNRebin[4]);

    //calulation background for fullPt without weighting
    fMesonFullPtBackground = new TH1D("Mapping_Back_InvMass_FullPt","Mapping_Back_InvMass_FullPt",fBckInvMassVSPtDummy->GetNbinsX(),0.,1.);
    fMesonFullPtBackground->Sumw2();
    Int_t startBinFullPt = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fFullPt[0]+0.001);
    Int_t endBinFullPt = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fFullPt[1]-0.001);
    fBckInvMassVSPtDummy->ProjectionX("Mapping_Back_InvMass_FullPt",startBinFullPt,endBinFullPt);
    fMesonFullPtBackground=(TH1D*)gDirectory->Get("Mapping_Back_InvMass_FullPt");
    fMesonFullPtBackground->Rebin(fNRebin[4]);

    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoBack = Form("Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
        if(fHistoMappingBackInvMassPtBin[iPt]!= NULL){
            delete fHistoMappingBackInvMassPtBin[iPt];
            fHistoMappingBackInvMassPtBin[iPt]=NULL;
        }
        fHistoMappingBackInvMassPtBin[iPt]=new TH1D(fNameHistoBack.Data(),fNameHistoBack.Data(),fBckInvMassVSPtDummy->GetNbinsX(),0.,1.);
            fHistoMappingBackInvMassPtBin[iPt]->Sumw2();
        Int_t startBin = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
        Int_t endBin = fBckInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

        fBckInvMassVSPtDummy->ProjectionX(fNameHistoBack.Data(),startBin,endBin);
        fHistoMappingBackInvMassPtBin[iPt]=(TH1D*)gDirectory->Get(fNameHistoBack.Data());
        if(fNRebin[iPt]>1){
            fHistoMappingBackInvMassPtBin[iPt]->Rebin(fNRebin[iPt]);
        }

    }
}


void FillMassHistosArray(TH2D* fGammaGammaInvMassVSPtDummy)
{
    fFittingHistMidPtSignal = new TH1D("Mapping_GG_InvMass_MidPt","Mapping_GG_InvMass_MidPt",fGammaGammaInvMassVSPtDummy->GetNbinsX(),0.,(double)fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()));
    cout << fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()) << endl;
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
    //fMesonFullPtSignal->Scale(1./fNRebin[4]);
        cout << "Full pt geschrieben" << endl;

    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
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

void FitSubtractedInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle, Double_t* fMesonIntDeltaRangeFit, Int_t ptBin, Bool_t vary)
{

//    cout<<"Start Fitting spectra"<<endl;
    fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
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

void FitCBSubtractedInvMassInPtBins(TH1D* fHistoMappingSignalInvMassPtBinSingle,Double_t * fMesonIntDeltaRangeFit, Int_t ptBin,Bool_t vary ,TString functionname)
{


fFileErrLog<<"Start Fitting spectra with CB fit"<<endl;
fHistoMappingSignalInvMassPtBinSingle->GetXaxis()->SetRangeUser(fMesonMassRange[0],fMesonMassRange[1]);
Double_t mesonAmplitude =fHistoMappingSignalInvMassPtBinSingle->GetMaximum();
Double_t mesonAmplitudeMin = mesonAmplitude*50./100.;
Double_t mesonAmplitudeMax = mesonAmplitude*115./100.;

fFitReco = NULL;
fFitReco = new TF1(functionname,CrystalBallBck,fMesonFitRange[0],fMesonFitRange[1],7);

fFitGausExp = NULL;
fFitGausExp = new TF1("optionCrystalBall",CrystalBall,fMesonFitRange[0],fMesonFitRange[1],5);

fFitLinearBck = NULL;
fFitLinearBck = new TF1("Linear","[0]+[1]*x",fMesonFitRange[0],fMesonFitRange[1]);

fFitReco->SetParameter(0,mesonAmplitude);
fFitReco->SetParameter(1,fMesonMassExpect);
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
fFitReco->SetParLimits(1,fMesonMassRange[0],fMesonMassRange[1]);
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

void FitWithPol2ForBG(TH1D* fHistoMappingSignalInvMassPtBinSingle,Double_t * fMesonFitRangeCur, Int_t ptBin, Bool_t vary)
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
fFitReco->SetParameter(1,fMesonMassExpect);
fFitReco->SetParameter(2,fMesonWidthExpect);
if(vary){
    fFitReco->SetParameter(3,fMesonLambdaTail);
} else {
    fFitReco->FixParameter(3,fMesonLambdaTail);
}
fFitReco->SetParLimits(0,mesonAmplitudeMin,mesonAmplitudeMax);
fFitReco->SetParLimits(1,fMesonMassExpect*95./100.,fMesonMassExpect*105./100.);
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
    fMesonFullPtBackground->Write();
    fMesonFullPtBackNorm->SetName("Mapping_BackNorm_InvMass_FullPt");
    fMesonFullPtBackNorm->Write();
    fNumberOfGoodESDTracks->Write();
    fEventQuality->Write();

    TString nameHistoSignal;
    TString nameHistoBckNorm;
    TString fitnameSignal;
    TString nameHistoSignalPos;
    for(Int_t ii =fStartPtBin;ii<fNBinsPt;ii++){
        fHistoMappingGGInvMassPtBin[ii]->Write();
        nameHistoBckNorm = Form("Mapping_BckNorm_InvMass_in_Pt_Bin%02d", ii);
        fHistoMappingBackNormInvMassPtBin[ii]->Write(nameHistoBckNorm.Data());
        nameHistoSignal = Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d", ii);
        fHistoMappingSignalInvMassPtBin[ii]->Write(nameHistoSignal.Data());
        fitnameSignal = Form("Signal_InvMassFit_in_Pt_Bin%02d", ii);
        if(fFitSignalInvMassPtBin[ii]!=0x00) fFitSignalInvMassPtBin[ii]->Write(fitnameSignal.Data());
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
    }

    // initialize integration pt-arrays for integration windows: normal, wide, narrow, left, left wide, left narrow
    for (Int_t k = 0; k < 6; k++){
        // Initialize yield arrays
        fGGYields[k] = 											new Double_t[fNBinsPt];
        fBckYields[k] =                                         new Double_t[fNBinsPt];
        fTotalBckYields[k] = 								    new Double_t[fNBinsPt];
        fMesonYields[k] = 									    new Double_t[fNBinsPt];
        fMesonYieldsResidualBckFunc[k] = 						new Double_t[fNBinsPt];
        fMesonYieldsCorResidualBckFunc[k] = 					new Double_t[fNBinsPt];
        fMesonYieldsPerEvent[k] = 								new Double_t[fNBinsPt];

        // Initialize error arrays
        fGGYieldsError[k] = 								    new Double_t[fNBinsPt];
        fBckYieldsError[k] = 									new Double_t[fNBinsPt];
        fTotalBckYieldsError[k] = 								new Double_t[fNBinsPt];
        fMesonYieldsError[k] = 									new Double_t[fNBinsPt];
        fMesonYieldsResidualBckFuncError[k] = 					new Double_t[fNBinsPt];
        fMesonYieldsCorResidualBckFuncError[k] =				new Double_t[fNBinsPt];
        fMesonYieldsPerEventError[k] = 							new Double_t[fNBinsPt];
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
    fMesonYieldsBackFit =									new Double_t[fNBinsPt];
    fMesonYieldsResidualBckFuncBackFit = 					new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncBackFit = 				new Double_t[fNBinsPt];
    fMesonYieldsPerEventBackFit = 							new Double_t[fNBinsPt];
    fMesonMass = 											new Double_t[fNBinsPt];
    fMesonMassBackFit = 									new Double_t[fNBinsPt];
    fMesonWidth = 											new Double_t[fNBinsPt];
    fMesonWidthBackFit = 									new Double_t[fNBinsPt];
    fMesonFWHM = 											new Double_t[fNBinsPt];
    fMesonTrueMass = 										new Double_t[fNBinsPt];
    fMesonTrueMassReweighted = 								new Double_t[fNBinsPt];
    fMesonTrueFWHM = 										new Double_t[fNBinsPt];
    fMesonTrueFWHMReweighted =								new Double_t[fNBinsPt];
    fMesonTrueSB = 											new Double_t[fNBinsPt];
    fMesonTrueSign = 										new Double_t[fNBinsPt];

    // Normalization at the left of the peak
    fMesonMassLeft = 										new Double_t[fNBinsPt];
    fMesonWidthLeft = 										new Double_t[fNBinsPt];
    fMesonFWHMLeft = 										new Double_t[fNBinsPt];

    fMesonYieldsBackFitError =								new Double_t[fNBinsPt];

    fMesonYieldsResidualBckFuncBackFitError = 				new Double_t[fNBinsPt];
    fMesonYieldsCorResidualBckFuncBackFitError =			new Double_t[fNBinsPt];
    fMesonYieldsPerEventBackFitError = 						new Double_t[fNBinsPt];
    fMesonMassError = 										new Double_t[fNBinsPt];
    fMesonMassBackFitError = 								new Double_t[fNBinsPt];
    fMesonWidthError = 										new Double_t[fNBinsPt];
    fMesonWidthBackFitError = 								new Double_t[fNBinsPt];
    fMesonFWHMError = 										new Double_t[fNBinsPt];
    fMesonTrueMassError = 									new Double_t[fNBinsPt];
    fMesonTrueMassReweightedError =		 					new Double_t[fNBinsPt];
    fMesonTrueFWHMError = 									new Double_t[fNBinsPt];
    fMesonTrueFWHMReweightedError = 						new Double_t[fNBinsPt];

    fMesonMassLeftError = 									new Double_t[fNBinsPt];
    fMesonWidthLeftError = 									new Double_t[fNBinsPt];
    fMesonFWHMLeftError = 									new Double_t[fNBinsPt];

    // nitialize pt-arrays for different mass & width fitting procedures, normalization at the left of the peak
    fMesonMassLeft                                                  = new Double_t[fNBinsPt];
    fMesonFWHMLeft                                                  = new Double_t[fNBinsPt];
    fMesonMassLeftError                                             = new Double_t[fNBinsPt];
    fMesonFWHMLeftError                                             = new Double_t[fNBinsPt];

    // initialize pt-arrays for validated Significance and S/B
    fMesonTrueSB                                                    = new Double_t[fNBinsPt];
    fMesonTrueSign                                                  = new Double_t[fNBinsPt];

    fHistoMappingTrueMesonInvMassPtBins = 					new TH1D*[fNBinsPt];
    fHistoMappingTrueMesonInvMassPtReweightedBins =  		new TH1D*[fNBinsPt];

    fHistoWeightsBGZbinVsMbin = 							new TH2F*[fNBinsPt];
    fHistoFillPerEventBGZbinVsMbin = 						new TH2F*[fNBinsPt];

    fHistoMappingGGInvMassPtBin = 							new TH1D*[fNBinsPt];
    fHistoMappingGGInvMassBackFitPtBin =					new TH1D*[fNBinsPt];
    fHistoMappingGGInvMassBackFitWithoutSignalPtBin =		new TH1D*[fNBinsPt];
    fHistoMappingBackInvMassPtBin = 						new TH1D*[fNBinsPt];
    fHistoMappingBackNormInvMassPtBin = 					new TH1D*[fNBinsPt];
    fHistoMappingSignalInvMassPtBin =						new TH1D*[fNBinsPt];
//  fRatio =                                                new TH1D*[fNBinsPt];
    fHistoMappingRatioSBInvMassPtBin= 						new TH1D*[fNBinsPt];

    fBackgroundFitPol = 									new TF1*[fNBinsPt];
    fFitSignalInvMassPtBin = 								new TF1*[fNBinsPt];
    fFitSignalInvMassBackFitPtBin = 						new TF1*[fNBinsPt];
    fFitSignalPeakPosInvMassLeftPtBin                     = new TF1*[fNBinsPt];
    fFitTrueSignalInvMassPtBin = 							new TF1*[fNBinsPt];
    fFitTrueSignalInvMassPtReweightedBin =					new TF1*[fNBinsPt];

    fFitBckInvMassPtBin = 									new TF1*[fNBinsPt];
    fFitBckInvMassBackFitPtBin =							new TF1*[fNBinsPt];
    fFitRatioInvMassPtBin = 								new TF1*[fNBinsPt];
    // Histograms for normalization on the left of the peak
    fHistoMappingBackNormInvMassLeftPtBin = 				new TH1D*[fNBinsPt];
    fHistoMappingSignalInvMassLeftPtBin = 					new TH1D*[fNBinsPt];

    fFitInvMassLeftPtBin = 									new TF1*[fNBinsPt];
    fFitBckInvMassLeftPtBin = 								new TF1*[fNBinsPt];
    fFitWithPol2ForBG = 									new TF1*[fNBinsPt];

    for(Int_t i = 0;i<fNBinsPt; i++){
        fHistoMappingTrueMesonInvMassPtBins[i] = 					NULL;
        fHistoMappingTrueMesonInvMassPtReweightedBins[i] =			NULL;// array of histos for pt slices

        fHistoMappingGGInvMassPtBin[i] = 							NULL;
        fHistoMappingGGInvMassBackFitPtBin[i] = 					NULL;
        fHistoMappingGGInvMassBackFitWithoutSignalPtBin[i] = 		NULL;
        fBackgroundFitPol[i] = 										NULL;
        fHistoMappingBackInvMassPtBin[i] = 							NULL;
        fHistoMappingBackNormInvMassPtBin[i] = 						NULL;
//      fRatio[i] =                                                 NULL;
        fHistoMappingSignalInvMassPtBin[i] = 						NULL;
        fHistoMappingRatioSBInvMassPtBin[i] = 						NULL;

        fFitSignalInvMassPtBin[i] = 								NULL;
        fFitSignalInvMassBackFitPtBin[i] = 							NULL;
        fFitTrueSignalInvMassPtBin[i] = 							NULL;
        fFitTrueSignalInvMassPtReweightedBin[i] = 					NULL;

        fFitBckInvMassPtBin[i] = 									NULL;
        fFitBckInvMassBackFitPtBin[i] = 							NULL;
        fFitRatioInvMassPtBin[i] = 									NULL;
        // Histograms for normalization on the left of the peak
        fHistoMappingBackNormInvMassLeftPtBin[i] = 					NULL;
        fHistoMappingSignalInvMassLeftPtBin[i] = 					NULL;

        fFitInvMassLeftPtBin[i] = 									NULL;
        fFitBckInvMassLeftPtBin[i] = 								NULL;
        fFitWithPol2ForBG[i] = 										NULL;
    }
}

void CalculateFWHM(TF1 * fFunc){
// Default function
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
    //	fFunc_plus = fFunc;
    fFunc_plus = new TF1("fFunc_plus","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",fMesonFitRange[0],fMesonFitRange[1]);

    fFunc_plus->SetParameter(0,fFunc->GetParameter(0) + fFunc->GetParError(0));
    fFunc_plus->SetParameter(1,fFunc->GetParameter(1) + fFunc->GetParError(1));
    fFunc_plus->SetParameter(2,fFunc->GetParameter(2) + fFunc->GetParError(2));
    fFunc_plus->SetParameter(3,fFunc->GetParameter(3) + fFunc->GetParError(3));
    Double_t FWHM_plus = fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,fFunc_plus->GetParameter(1), fMesonFitRange[1]) - fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,fMesonFitRange[0],fFunc_plus->GetParameter(1));

    //FWHM error -
    TF1* fFunc_minus;
    //	fFunc_minus = fFunc;
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
    if (fFitRange) delete fFitRange;
    if (fBGFitRange) delete fBGFitRange;
    if (fBGFitRangeLeft) delete fBGFitRangeLeft;
    if (fMesonPlotRange) delete fMesonPlotRange;
    if (fMesonIntRange) delete fMesonIntRange;
    if (fMesonMassRange) delete fMesonMassRange;
    if (fMesonFitRange) delete fMesonFitRange;
    if (fMesonWidthRange) delete fMesonWidthRange;
    if (fMesonLambdaTailRange) delete fMesonLambdaTailRange;
    if (fNRebin) delete fNRebin;
    for (Int_t k = 0; k< 6; k++){
        if (fGGYields[k])           delete fGGYields[k];
        if (fBckYields[k])          delete fBckYields[k];
        if (fBckYieldsError[k])     delete fBckYieldsError[k];
        if (fMesonYields[k])        delete fMesonYields[k];
        if (fMesonYieldsResidualBckFunc[k])         delete fMesonYieldsResidualBckFunc[k];
        if (fMesonYieldsResidualBckFuncError[k])    delete fMesonYieldsResidualBckFuncError[k];
        if (fMesonYieldsPerEventError[k])           delete fMesonYieldsPerEventError[k];
        if (fMesonYieldsPerEvent[k])                delete fMesonYieldsPerEvent[k];
        if (fMesonYieldsError[k])                   delete fMesonYieldsError[k];
        if (fMesonYieldsCorResidualBckFuncError[k]) delete fMesonYieldsCorResidualBckFuncError[k];
    }
    for (Int_t k = 0; k< 3; k++){
        if (fMesonTrueYields[k])            delete fMesonTrueYields[k];
        if (fMesonTrueYieldsReweighted[k])  delete fMesonTrueYieldsReweighted[k];
    }


    if (fMesonYieldsBackFit)                    delete fMesonYieldsBackFit;
    if (fMesonYieldsResidualBckFuncBackFit)     delete fMesonYieldsResidualBckFuncBackFit;

    if (fMesonYieldsCorResidualBckFuncBackFit)  delete fMesonYieldsCorResidualBckFuncBackFit;
    if (fMesonYieldsPerEventBackFit)            delete fMesonYieldsPerEventBackFit;
    if (fMesonMass)                             delete fMesonMass;
    if (fMesonWidth)                            delete fMesonWidth;
    if (fMesonTrueSB)                           delete fMesonTrueSB;
    if (fMesonTrueSign)                         delete fMesonTrueSign;
    if (fMesonFWHM)                             delete fMesonFWHM;
    if (fMesonMassLeft)                         delete fMesonMassLeft;
    if (fMesonWidthLeft)                        delete fMesonWidthLeft;
    if (fMesonFWHMLeft)                         delete fMesonFWHMLeft;
    if (fMesonYieldsBackFitError)               delete fMesonYieldsBackFitError;


    if (fMesonYieldsResidualBckFuncBackFitError)    delete fMesonYieldsResidualBckFuncBackFitError;
    if (fMesonYieldsCorResidualBckFuncBackFitError) delete fMesonYieldsCorResidualBckFuncBackFitError;
    if (fMesonYieldsPerEventBackFitError)       delete fMesonYieldsPerEventBackFitError;
    if (fMesonMassError)                        delete fMesonMassError;
    if (fMesonWidthError)                       delete fMesonWidthError;
    if (fMesonTrueSBError)                      delete fMesonTrueSBError;
    if (fMesonTrueSignError)                    delete fMesonTrueSignError;
    if (fMesonFWHMError) delete fMesonFWHMError;
    if (fMesonFWHMLeftError) delete fMesonFWHMLeftError;
    if (fHistoMappingTrueMesonInvMassPtBins) delete fHistoMappingTrueMesonInvMassPtBins;
    if (fHistoMappingTrueMesonInvMassPtReweightedBins) delete fHistoMappingTrueMesonInvMassPtReweightedBins;
    if (fHistoMappingGGInvMassPtBin) delete fHistoMappingGGInvMassPtBin;
    if (fHistoMappingBackInvMassPtBin) delete fHistoMappingBackInvMassPtBin;
    if (fHistoMappingBackNormInvMassPtBin) delete fHistoMappingBackNormInvMassPtBin;
    if (fHistoMappingSignalInvMassPtBin) delete fHistoMappingSignalInvMassPtBin;
    if (fHistoMappingRatioSBInvMassPtBin) delete fHistoMappingRatioSBInvMassPtBin;
    if (fFitSignalInvMassPtBin) delete fFitSignalInvMassPtBin;
    if (fFitBckInvMassPtBin) delete fFitBckInvMassPtBin;
    if (fHistoMappingBackNormInvMassLeftPtBin) delete fHistoMappingBackNormInvMassLeftPtBin;
    if (fHistoMappingSignalInvMassLeftPtBin) delete fHistoMappingSignalInvMassLeftPtBin;
    if (fFitInvMassLeftPtBin) delete fFitInvMassLeftPtBin;
    if (fFitSignalPeakPosInvMassLeftPtBin) delete fFitSignalPeakPosInvMassLeftPtBin;
    if (fFitBckInvMassLeftPtBin) delete fFitBckInvMassLeftPtBin;
    if (fFitRatioInvMassPtBin) delete fFitRatioInvMassPtBin;
    if (fFitWithPol2ForBG)                  delete fFitWithPol2ForBG;
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

