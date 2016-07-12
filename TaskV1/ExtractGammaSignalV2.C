// provided by Gamma Conversion Group, $ALICE_PHYSCIS/PWGGA/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion


#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
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
#include "TFitResultPtr.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TDatabasePDG.h"
#include "TVirtualFitter.h"
#include "TMinuit.h"
#include "TTree.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/PlottingMeson.h"
#include "../CommonHeaders/FittingGammaConversion.h"
//#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "ExtractSignalV2.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "../CommonHeaders/ExtractSignalPlotting.h"
#include "ExtractGammaSignalV2.h"

//**************************************************************************************************
//**********************************  Main Function ************************************************
//**************************************************************************************************
void ExtractGammaSignalV2(      TString meson               = "",
                                TString file                = "",
                                TString cutSelection        = "",
                                TString suffix              = "",
                                TString isMC                = "",
                                TString option              = "",
                                TString directphotonPlots   = "",
                                TString period              = "",
                                Int_t numberOfBins          = 30,
                                Bool_t addSig               = 0,
                                Int_t mode                  = 0
                            ) {

    //********************************* Catch modes which are not supported ****************************
    if (mode == 9) {
        cout << "ERROR: this mode is not supported anymore" << endl;
        return;
    } else if ( mode == 1 ){
        cout << "ERROR: you can't run the photon extraction in the Dalitz mode" << endl;
        return;
    } else if ( mode == 2 ||  mode == 3 ){
        cout << "WARNING: you are running in hybrid mode the software is still under construction for this one" << endl;
        fEnablePCM  = 1;
        fEnableCalo = 1;
    } else if ( mode == 4 ||  mode == 5 ){
        cout << "WARNING: you are running in calo mode the software is still under construction for this one" << endl;
        fEnableCalo = 1;
    } else if ( mode == 0){
        fEnablePCM  = 1;
    }    
    
    //************************************* Set general style settings *********************************
    StyleSettingsThesis();
    SetPlotStyle();

    //************************************ Define Output directory *************************************
    fOutputDir = Form("%s/%s/%s/ExtractGammaSignal",cutSelection.Data(),option.Data(),suffix.Data());
    gSystem->Exec("mkdir -p "+fOutputDir);

    //************************************ Set global variables ****************************************
    fDate                                       = ReturnDateString();
    fDirectPhoton                               = directphotonPlots;
    fEnergyFlag                                 = option;
    fPrefix                                     = meson;
    fPeriodFlag                                 = period;
    fSuffix                                     = suffix;
    fMeson                                      = meson;
    fMode                                       = mode;
    cout << "Pictures are saved as " << suffix.Data() << endl;

    //************************ Detect correct folder name ****************************************
    TString nameMainDir                         = "";
    if (fMode == 0) 
        nameMainDir                             = "GammaConvV1";
    else if (fMode == 2 || fMode == 3) 
        nameMainDir                             = "GammaConvCalo";
    else if (fMode == 4 || fMode == 5) 
        nameMainDir                             = "GammaCalo";

    
    //************************************ Separate cutstrings ***********************************
    fCutSelection                               = cutSelection;
    fCutSelectionRead                           = cutSelection;
    ReturnSeparatedCutNumberAdvanced( fCutSelection,fEventCutNumber, fGammaCutNumber, fClusterCutNumber, fElectronCutNumber, fMesonCutNumber, fMode);
    
    fEventCutSelectionRead                      = fEventCutNumber.Data();
    fGammaCutSelectionRead                      = fGammaCutNumber.Data();
    fMesonCutSelectionRead                      = fMesonCutNumber.Data();
    if (addSig) {
        cout << "running added Signal" << endl;
        cout << fEventCutNumber.Data() << endl;
        fEventCutNumber.Replace(GetEventRejectExtraSignalsCutPosition(),1,"2");
        cout << fEventCutNumber.Data() << endl;
        fEventCutSelectionRead                  = fEventCutNumber;
        fGammaCutSelectionRead                  = fGammaCutNumber;
        fMesonCutSelectionRead                  = fMesonCutNumber;
        if (fMode==9) fCutSelectionRead         = Form("%s%s_%s", fEventCutNumber.Data(), fGammaCutNumber.Data(), fMesonCutNumber.Data());
        else if (fMode==0) fCutSelectionRead    = Form("%s_%s_%s",fEventCutNumber.Data(), fGammaCutNumber.Data(), fMesonCutNumber.Data());
        cout << fCutSelectionRead.Data() << endl;
    }

    //***************************** Load binning for spectrum *******************************************
    Initialize(fMeson, fEnergyFlag, numberOfBins, fMode, addSig);

    //****************************** Set specific histogram names****************************************
    ObjectNameDCGammaConvRPt                = "ESD_TrueDoubleCountConvGamma_R_Pt";
    ObjectNameGammaConvMultipleCount        = "ESD_TrueMultipleCountConvGamma";

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
    
    //***************************** Specification Data/MC ***********************************************
    if(isMC.CompareTo("kTRUE") == 0){
        fIsMC                               = 1;
        fPrefix2                            = "MC";
    } else {
        fIsMC                               = 0;
        fPrefix2                            = "data";
    }

    //************************************** Read file ***************************************************
    TFile f(file.Data());

    TList *TopDir =(TList*)f.Get(nameMainDir.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }
    TList* HistosGammaConversion            = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if (!HistosGammaConversion){
      cout << "ERROR: folder with Cutnumber - " <<   fCutSelectionRead.Data() << "not contained in file " << file.Data() << endl;
      return;
    }    
    TList* ESDContainer                     = (TList*)HistosGammaConversion->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));
    if (fEnablePCM){
        TList* ConvCutsContainer            = (TList*) HistosGammaConversion->FindObject(Form("ConvCuts_%s",fGammaCutSelectionRead.Data()));

        fHistoPhotonIsSelected              = (TH1D*)ConvCutsContainer->FindObject(Form("IsPhotonSelected %s",fGammaCutSelectionRead.Data()));
        fHistoPhotonIsSelected->Scale(1./fHistoPhotonIsSelected->GetEntries());
    }
    
    // general histograms
    fNumberOfGoodESDTracks                  = (TH1D*)ESDContainer->FindObject("GoodESDTracks");
    fEventQuality                           = (TH1D*)ESDContainer->FindObject("NEvents");
    
    if (option.CompareTo("PbPb_2.76TeV") == 0){
        fNEvnt                              = fEventQuality->GetBinContent(1);
    } else {
        fNEvnt                              = GetNEvents(fEventQuality);
    }
    
    // *******************************************************************************************************
    // ******************************* Read histograms from data for PCM *************************************
    // *******************************************************************************************************
    if (fEnablePCM){
        // read reconstructed gamma histograms
        fHistoGammaConvPt                   = (TH1D*)ESDContainer->FindObject("ESD_ConvGamma_Pt");
        fHistoGammaConvPt->Sumw2();
        fHistoGammaConvPtOrBin              = (TH1D*)fHistoGammaConvPt->Clone("ESD_ConvGamma_Pt_OriginalBinning");
        RebinSpectrum(fHistoGammaConvPt,"");
        
        // read dca tree for conversions
        if(!addSig && !(mode == 4 || mode == 5)){
            DCAContainer                    = (TList*) HistosGammaConversion->FindObject(Form("%s Photon DCA tree",fCutSelectionRead.Data()));
            if(DCAContainer){
                dcaTree                     = (TTree*)DCAContainer->FindObject("ESD_ConvGamma_Pt_Dcaz_R_Eta");
                FillDCAHistogramsFromTree(dcaTree,kFALSE);
                CalculatePileUpBackground(kFALSE);
                pileUpCorrection            = kTRUE;
            }         
        }
    }

    // *******************************************************************************************************
    // ******************************* Read histograms from data for Calo ************************************
    // *******************************************************************************************************
    if (fEnableCalo){
        if (mode == 2 || mode == 3){
            CaloContainer                           = (TList*) HistosGammaConversion->FindObject(Form("%s Cluster Output",fCutSelectionRead.Data()));
            fHistoGammaCaloPt                       = (TH1D*) CaloContainer->FindObject("ClusGamma_Pt");
            fHistoGammaCaloPt->Sumw2();
            fHistoGammaCaloPtOrBin                  = (TH1D*)fHistoGammaCaloPt->Clone("ClusGamma_Pt_OriginalBinning");
            RebinSpectrum(fHistoGammaCaloPt,"");
        } else if ( mode == 4 || mode == 5){
            fHistoGammaCaloPt                       = (TH1D*) ESDContainer->FindObject("ClusGamma_Pt");
            fHistoGammaCaloPt->Sumw2();
            fHistoGammaCaloPtOrBin                  = (TH1D*)fHistoGammaCaloPt->Clone("ClusGamma_Pt_OriginalBinning");
            RebinSpectrum(fHistoGammaCaloPt,"");
        }    
    }
    
    // read  MC quantities from file
    if(fIsMC){
        // copy reconstructed photon histo for MC
        if (fEnablePCM){
            fHistoGammaMCrecConvPt                              = (TH1D*) fHistoGammaConvPt->Clone("MCrec_ConvGamma_Pt");
            fHistoGammaMCrecConvPtOrBin                         = (TH1D*) fHistoGammaConvPtOrBin->Clone("MCrec_ConvGamma_Pt_OriginalBinning");
        }
        if (fEnableCalo){
            fHistoGammaMCrecCaloPt                              = (TH1D*) fHistoGammaCaloPt->Clone("MCrec_CaloGamma_Pt");
        }
        
        // MC container (contains all input MC histograms)    
        TList *MCContainer                                      = (TList*)HistosGammaConversion->FindObject(Form("%s MC histograms",fCutSelectionRead.Data()));
        // reading input distributions
        if (fEnablePCM){
            fHistoGammaMCConvPt                                 = (TH1D*)MCContainer->FindObject("MC_ConvGamma_Pt");
            fHistoGammaMCConvPt->Sumw2();
            RebinSpectrum(fHistoGammaMCConvPt);
        }    
        if (fEnableCalo && mode == 2 ){
            fHistoGammaMCAllInEMCAccPt                          = (TH1D*)MCContainer->FindObject("MC_AllGammaEMCALAcc_Pt");
            fHistoGammaMCAllInEMCAccPt->Sumw2();
            fHistoGammaMCAllInEMCAccPtOrBin                     = (TH1D*)fHistoGammaMCAllInEMCAccPt->Clone("MC_AllGammaEMCALAcc_OriginalBinning_MCPt");
            fHistoGammaMCAllInEMCAccPtOrBin->Scale(1./fHistoGammaMCAllInEMCAccPtOrBin->GetBinWidth(5));
            fHistoGammaMCAllInEMCAccPtOrBin->GetXaxis()->SetRangeUser(0.,25.);
            RebinSpectrum(fHistoGammaMCAllInEMCAccPt);
        }    
            
        fHistoGammaMCAllPt                                      = (TH1D*)MCContainer->FindObject("MC_AllGamma_Pt");
        fHistoGammaMCAllPt->Sumw2();
        fHistoGammaMCAllPtOrBin                                 = (TH1D*)fHistoGammaMCAllPt->Clone("MC_AllGamma_OriginalBinning_MCPt");
        fHistoGammaMCAllPtOrBin->Scale(1./fHistoGammaMCAllPtOrBin->GetBinWidth(5));
//         fHistoGammaMCAllPtOrBin->GetXaxis()->SetRangeUser(0.,25.);
        RebinSpectrum(fHistoGammaMCAllPt);
        
        // Gamma from Decay
        fHistoGammaMCDecayPt                                    = new TH1D*[7];
        for(Int_t i = 0; i<7; i++){
            fHistoGammaMCDecayPt[i]                             = (TH1D*)MCContainer->FindObject(Form("MC_DecayGamma%s_Pt",fDecays[i].Data()));
            fHistoGammaMCDecayPt[i]->Sumw2();
            fHistoGammaMCDecayPt[i]->Scale(1./fHistoGammaMCDecayPt[i]->GetBinWidth(5));
        }
        
        // container with histos for validated reconstructed photons
        TList *TrueConversionContainer                          = (TList*)HistosGammaConversion->FindObject(Form("%s True histograms",fCutSelectionRead.Data()));
        // reading reconstructed validated distributions
        if (fEnablePCM){
            fHistoGammaTrueConvPt                                   = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueConvGamma_Pt");
            fHistoGammaTrueConvPt->Sumw2();
            fHistoGammaTrueConvPtOrBin                              = (TH1D*)fHistoGammaTrueConvPt->Clone("TrueConvGamma_Pt_OriginalBinning");
            fHistoGammaTrueConvPtOrBin->Sumw2();
            RebinSpectrum(fHistoGammaTrueConvPt);
    
            fHistoGammaTruePrimaryConvPt                            = (TH1D*)TrueConversionContainer->FindObject("ESD_TruePrimaryConvGamma_Pt");
            fHistoGammaTruePrimaryConvPt->Sumw2();
            fHistoGammaTruePrimaryConvPtOrBin                       = (TH1D*)fHistoGammaTruePrimaryConvPt->Clone("TruePrimaryConvGamma_Pt_OriginalBinning");
            RebinSpectrum(fHistoGammaTruePrimaryConvPt);
    
            
            fHistoTrueGammaConvDCRVSPt                              = (TH2D*)TrueConversionContainer->FindObject(ObjectNameDCGammaConvRPt.Data());
            if (fHistoTrueGammaConvDCRVSPt != NULL) 
                fEnableDCConv                                       = kTRUE;
            if (fEnableDCConv){
                fHistoTrueGammaConvDCRVSPt->Sumw2();
                fHistoTrueGammaConvDCPt                             = (TH1D*)fHistoTrueGammaConvDCRVSPt->ProjectionY("fHistoTrueGammaConvDCPt");
                fHistoTrueGammaConvDCR                              = (TH1D*)fHistoTrueGammaConvDCRVSPt->ProjectionX("fHistoTrueGammaConvDCR");
                fHistoTrueGammaConvMultipleCount                    = (TH1F*)TrueConversionContainer->FindObject(ObjectNameGammaConvMultipleCount.Data());
            }
            
            fHistoGammaTrueSecondaryConvPt                          = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryConvGamma_Pt");
            fHistoGammaTrueSecondaryConvPt->Sumw2();
            fHistoGammaTrueSecondaryConvPtOrBin                     = (TH1D*)fHistoGammaTrueSecondaryConvPt->Clone("ESD_TrueSecondaryConvGamma_Pt_OriginalBinning");
            RebinSpectrum(fHistoGammaTrueSecondaryConvPt);    
    
            fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt         = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryConvGammaFromXFromK0s_Pt");
            fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt->Sumw2();
            fHistoGammaTrueSecondaryConvGammaFromXFromK0sPtOrBin    = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt->Clone("ESD_TrueSecondaryConvGammaFromXFromK0s_Pt_OriginalBinning");
            RebinSpectrum(fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt);
    

            fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt      = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryConvGammaFromXFromLambda_Pt");
            fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt->Sumw2();
            fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPtOrBin = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt->Clone("ESD_TrueSecondaryConvGammaFromXFromLambda_Pt_OriginalBinning");
            RebinSpectrum(fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt);
            
            // combinatorial BG distributions
            fHistoCombinatorialBackground                           = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueCombinatorial_Pt");
            fHistoCombinatorialSpecies                              = new TH1D*[17];
            fHistoCombinatorialSpecies[nCombinatorics]              = (TH1D*)fHistoCombinatorialBackground->ProjectionX(Form("ESD_TrueComb%s_Pt",fCombinatorics[nCombinatorics].Data()));
            fHistoCombinatorialSpecies[nCombinatorics]->Sumw2();
            RebinSpectrum(fHistoCombinatorialSpecies[nCombinatorics]);
    
            for(Int_t i = 0; i<nCombinatorics; i++){
                fHistoCombinatorialSpecies[i]                       = (TH1D*)fHistoCombinatorialBackground->ProjectionX(Form("ESD_TrueComb%s_Pt",fCombinatorics[i].Data()),i+1,i+1);
                fHistoCombinatorialSpecies[i]->Sumw2();
                RebinSpectrum(fHistoCombinatorialSpecies[i]);
        
            }

            // Response matrix PCM
            fHistoGammaTruePrimaryConv_recPt_MCPt_MC                = (TH2D*)TrueConversionContainer->FindObject("ESD_TruePrimaryConvGammaESD_PtMCPt");
            fHistoGammaTruePrimaryConv_recPt_MCPt_MC->Sumw2();
            fHistoGammaTruePrimaryConvMCPt                          = (TH1D*)fHistoGammaTruePrimaryConv_recPt_MCPt_MC->ProjectionY("ESD_TruePrimaryConvGamma_MCPt");
            RebinSpectrum(fHistoGammaTruePrimaryConvMCPt);
        }
        
        if (fEnableCalo){
            if (mode == 2 || mode == 3){
                fHistoGammaTrueCaloPt                           = (TH1D*)CaloContainer->FindObject("TrueClusGamma_Pt");
                fHistoGammaTrueCaloPt->Sumw2();
                fHistoGammaTrueCaloPtOrBin                      = (TH1D*)fHistoGammaTrueCaloPt->Clone("TrueClusGamma_Pt_OriginalBinning");
                fHistoGammaTrueCaloPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTrueCaloPt);

                fHistoGammaTruePrimaryCaloPt                    = (TH1D*)CaloContainer->FindObject("TruePrimaryClusGamma_Pt");
                fHistoGammaTruePrimaryCaloPt->Sumw2();
                fHistoGammaTruePrimaryCaloPtOrBin               = (TH1D*)fHistoGammaTruePrimaryCaloPt->Clone("TruePrimaryClusGamma_Pt_OriginalBinning");
                RebinSpectrum(fHistoGammaTruePrimaryCaloPt);

                fHistoGammaTrueSecondaryCaloPt                  = (TH1D*)CaloContainer->FindObject("TrueSecondaryClusGamma_Pt");
                fHistoGammaTrueSecondaryCaloPt->Sumw2();
                fHistoGammaTrueSecondaryCaloPtOrBin             = (TH1D*)fHistoGammaTrueSecondaryCaloPt->Clone("TrueSecondaryClusGamma_Pt_OriginalBinning");
                RebinSpectrum(fHistoGammaTrueSecondaryCaloPt);

                fHistoGammaTrueSecondaryCaloFromK0sPt           = (TH1D*)CaloContainer->FindObject("TrueSecondaryClusGammaFromK0s_Pt");
                fHistoGammaTrueSecondaryCaloFromK0sPt->Sumw2();
                fHistoGammaTrueSecondaryCaloFromK0sPtOrBin      = (TH1D*)fHistoGammaTrueSecondaryCaloFromK0sPt->Clone("TrueSecondaryClusGammaFromK0s_Pt_OriginalBinning");
                RebinSpectrum(fHistoGammaTrueSecondaryCaloFromK0sPt);

                fHistoGammaTrueSecondaryCaloFromLambdaPt        = (TH1D*)CaloContainer->FindObject("TrueSecondaryClusGammaFromLambda_Pt");
                fHistoGammaTrueSecondaryCaloFromLambdaPt->Sumw2();
                fHistoGammaTrueSecondaryCaloFromLambdaPtOrBin   = (TH1D*)fHistoGammaTrueSecondaryCaloFromLambdaPt->Clone("TrueSecondaryClusGammaFromLambda_Pt_OriginalBinning");
                RebinSpectrum(fHistoGammaTrueSecondaryCaloFromLambdaPt);
                
                // Response matrix calo
                fHistoGammaTruePrimaryCalo_recPt_MCPt_MC        = (TH2D*)CaloContainer->FindObject("TruePrimaryClusGamma_Pt_MCPt");
                fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->Sumw2();
                fHistoGammaTruePrimaryCaloMCPt                  = (TH1D*)fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->ProjectionY("TruePrimaryCaloGamma_MCPt");
                RebinSpectrum(fHistoGammaTruePrimaryCaloMCPt);

                // true conversion histograms
                fHistoGammaTrueCaloConvPt                       = (TH1D*)CaloContainer->FindObject("TrueClusConvGamma_Pt");
                fHistoGammaTrueCaloConvPt->Sumw2();
                fHistoGammaTrueCaloConvPtOrBin                  = (TH1D*)fHistoGammaTrueCaloConvPt->Clone("TrueClusConvGamma_Pt_OriginalBinning");
                fHistoGammaTrueCaloConvPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTrueCaloConvPt);

                fHistoGammaTruePrimaryCaloConvPt                = (TH1D*)CaloContainer->FindObject("TruePrimaryClusConvGamma_Pt");
                fHistoGammaTruePrimaryCaloConvPt->Sumw2();
                fHistoGammaTruePrimaryCaloConvPtOrBin           = (TH1D*)fHistoGammaTruePrimaryCaloConvPt->Clone("TruePrimaryClusConvGamma_Pt_OriginalBinning");
                fHistoGammaTruePrimaryCaloConvPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTruePrimaryCaloConvPt);

                fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC    = (TH2D*)CaloContainer->FindObject("TruePrimaryClusConvGamma_Pt_MCPt");
                fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->Sumw2();
                fHistoGammaTruePrimaryCaloConvMCPt              = (TH1D*)fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->ProjectionY("TruePrimaryCaloConvGamma_MCPt");
                RebinSpectrum(fHistoGammaTruePrimaryCaloConvMCPt);
                
                // true gamma histograms
                fHistoGammaTrueCaloUnConvPt                     = (TH1D*)fHistoGammaTrueCaloPt->Clone("TrueClusUnConvGamma_Pt");
                fHistoGammaTrueCaloUnConvPt->Sumw2();
                fHistoGammaTrueCaloUnConvPt->Add(fHistoGammaTrueCaloConvPt,-1);
                fHistoGammaTrueCaloUnConvPtOrBin                = (TH1D*)fHistoGammaTrueCaloPtOrBin->Clone("TrueClusUnConvGamma_Pt_OriginalBinning");
                fHistoGammaTrueCaloUnConvPtOrBin->Sumw2();
                fHistoGammaTrueCaloUnConvPtOrBin->Add(fHistoGammaTrueCaloConvPtOrBin,-1);

                fHistoGammaTruePrimaryCaloUnConvPt              = (TH1D*)fHistoGammaTruePrimaryCaloPt->Clone("TruePrimaryClusUnConvGamma_Pt");
                fHistoGammaTruePrimaryCaloUnConvPt->Sumw2();
                fHistoGammaTruePrimaryCaloUnConvPt->Add(fHistoGammaTruePrimaryCaloConvPt,-1);
                fHistoGammaTruePrimaryCaloUnConvPtOrBin         = (TH1D*)fHistoGammaTruePrimaryCaloPtOrBin->Clone("TruePrimaryClusUnConvGamma_Pt_OriginalBinning");
                fHistoGammaTruePrimaryCaloUnConvPtOrBin->Sumw2();
                fHistoGammaTruePrimaryCaloUnConvPtOrBin->Add(fHistoGammaTruePrimaryCaloConvPtOrBin,-1);

                fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC  = (TH2D*)fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->Clone("TruePrimaryClusUnConvGamma_Pt_MCPt");
                fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->Sumw2();
                fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->Add(fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC);
                fHistoGammaTruePrimaryCaloUnConvMCPt            = (TH1D*)fHistoGammaTruePrimaryCaloMCPt->Clone("TruePrimaryCaloUnConvGamma_MCPt");
                fHistoGammaTruePrimaryCaloUnConvMCPt->Sumw2();
                fHistoGammaTruePrimaryCaloUnConvMCPt->Add(fHistoGammaTruePrimaryCaloConvMCPt);
                
            } else {
                fHistoGammaTrueCaloUnConvPt                     = (TH1D*)TrueConversionContainer->FindObject("TrueClusGamma_Pt");
                fHistoGammaTrueCaloUnConvPt->SetName("TrueClusUnConvGamma_Pt");
                fHistoGammaTrueCaloUnConvPt->Sumw2();
                fHistoGammaTrueCaloUnConvPtOrBin                = (TH1D*)fHistoGammaTrueCaloUnConvPt->Clone("TrueClusUnConvGamma_Pt_OriginalBinning");
                fHistoGammaTrueCaloUnConvPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTrueCaloUnConvPt);

                fHistoGammaTrueCaloConvPt                       = (TH1D*)TrueConversionContainer->FindObject("TrueClusConvGamma_Pt");
                fHistoGammaTrueCaloConvPt->Sumw2();
                fHistoGammaTrueCaloConvPtOrBin                  = (TH1D*)fHistoGammaTrueCaloConvPt->Clone("TrueClusConvGamma_Pt_OriginalBinning");
                fHistoGammaTrueCaloConvPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTrueCaloConvPt);
                
                fHistoGammaTrueCaloPt                           = (TH1D*)fHistoGammaTrueCaloUnConvPt->Clone("TrueClusGamma_Pt");
                fHistoGammaTrueCaloPt->Add(fHistoGammaTrueCaloConvPt);
                fHistoGammaTrueCaloPtOrBin                      = (TH1D*)fHistoGammaTrueCaloUnConvPtOrBin->Clone("TrueClusGamma_Pt_OriginalBinning");
                fHistoGammaTrueCaloPtOrBin->Add(fHistoGammaTrueCaloConvPtOrBin);
                
                fHistoGammaTruePrimaryCaloUnConvPt              = (TH1D*)TrueConversionContainer->FindObject("TruePrimaryClusGamma_Pt");
                fHistoGammaTruePrimaryCaloUnConvPt->SetName("TruePrimaryClusUnConvGamma_Pt");
                fHistoGammaTruePrimaryCaloUnConvPt->Sumw2();
                fHistoGammaTruePrimaryCaloUnConvPtOrBin         = (TH1D*)fHistoGammaTruePrimaryCaloUnConvPt->Clone("TruePrimaryClusUnConvGamma_Pt_OriginalBinning");
                fHistoGammaTruePrimaryCaloUnConvPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTruePrimaryCaloUnConvPt);

                fHistoGammaTruePrimaryCaloConvPt                = (TH1D*)TrueConversionContainer->FindObject("TruePrimaryClusConvGamma_Pt");
                fHistoGammaTruePrimaryCaloConvPt->Sumw2();
                fHistoGammaTruePrimaryCaloConvPtOrBin           = (TH1D*)fHistoGammaTruePrimaryCaloConvPt->Clone("TruePrimaryClusConvGamma_Pt_OriginalBinning");
                fHistoGammaTruePrimaryCaloConvPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTruePrimaryCaloConvPt);
                
                fHistoGammaTruePrimaryCaloPt                    = (TH1D*)fHistoGammaTruePrimaryCaloUnConvPt->Clone("TruePrimaryClusGamma_Pt");
                fHistoGammaTruePrimaryCaloPt->Add(fHistoGammaTruePrimaryCaloConvPt);
                fHistoGammaTruePrimaryCaloPtOrBin               = (TH1D*)fHistoGammaTruePrimaryCaloUnConvPtOrBin->Clone("TruePrimaryClusGamma_Pt_OriginalBinning");
                fHistoGammaTruePrimaryCaloPtOrBin->Add(fHistoGammaTruePrimaryCaloConvPtOrBin);

                fHistoGammaTrueSecondaryCaloUnConvPt            = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusGamma_Pt");
                fHistoGammaTrueSecondaryCaloUnConvPt->SetName("TrueSecondaryClusUnConvGamma_Pt");
                fHistoGammaTrueSecondaryCaloUnConvPt->Sumw2();
                fHistoGammaTrueSecondaryCaloUnConvPtOrBin       = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvPt->Clone("TrueSecondaryClusUnConvGamma_Pt_OriginalBinning");
                fHistoGammaTrueSecondaryCaloUnConvPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTrueSecondaryCaloUnConvPt);

                fHistoGammaTrueSecondaryCaloConvPt              = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusConvGamma_Pt");
                fHistoGammaTrueSecondaryCaloConvPt->Sumw2();
                fHistoGammaTrueSecondaryCaloConvPtOrBin         = (TH1D*)fHistoGammaTrueSecondaryCaloConvPt->Clone("TrueSecondaryClusConvGamma_Pt_OriginalBinning");
                fHistoGammaTrueSecondaryCaloConvPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTrueSecondaryCaloConvPt);
                
                fHistoGammaTrueSecondaryCaloPt                  = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvPt->Clone("TrueSecondaryClusGamma_Pt");
                fHistoGammaTrueSecondaryCaloPt->Add(fHistoGammaTrueSecondaryCaloConvPt);
                fHistoGammaTrueSecondaryCaloPtOrBin             = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvPtOrBin->Clone("TrueSecondaryClusGamma_Pt_OriginalBinning");
                fHistoGammaTrueSecondaryCaloPtOrBin->Add(fHistoGammaTrueSecondaryCaloConvPtOrBin);
                
                fHistoGammaTrueSecondaryCaloUnConvFromK0sPt     = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusGammaFromXFromK0s_Pt");
                fHistoGammaTrueSecondaryCaloUnConvFromK0sPt->SetName("TrueSecondaryClusUnConvGammaFromK0s_Pt");
                fHistoGammaTrueSecondaryCaloUnConvFromK0sPt->Sumw2();
                fHistoGammaTrueSecondaryCaloUnConvFromK0sPtOrBin= (TH1D*)fHistoGammaTrueSecondaryCaloUnConvFromK0sPt->Clone("TrueSecondaryClusUnConvGammaFromK0s_Pt_OriginalBinning");
                fHistoGammaTrueSecondaryCaloUnConvFromK0sPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTrueSecondaryCaloUnConvFromK0sPt);

                fHistoGammaTrueSecondaryCaloConvFromK0sPt       = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusConvGammaFromXFromK0s_Pt");
                fHistoGammaTrueSecondaryCaloConvFromK0sPt->Sumw2();
                fHistoGammaTrueSecondaryCaloConvFromK0sPtOrBin  = (TH1D*)fHistoGammaTrueSecondaryCaloConvFromK0sPt->Clone("TrueSecondaryClusConvGammaFromK0s_Pt_OriginalBinning");
                fHistoGammaTrueSecondaryCaloConvFromK0sPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTrueSecondaryCaloConvFromK0sPt);
                
                fHistoGammaTrueSecondaryCaloFromK0sPt           = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvFromK0sPt->Clone("TrueSecondaryClusGammaFromK0s_Pt");
                fHistoGammaTrueSecondaryCaloFromK0sPt->Add(fHistoGammaTrueSecondaryCaloConvFromK0sPt);
                fHistoGammaTrueSecondaryCaloFromK0sPtOrBin      = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvFromK0sPtOrBin->Clone("TrueSecondaryClusGammaFromK0s_Pt_OriginalBinning");
                fHistoGammaTrueSecondaryCaloFromK0sPtOrBin->Add(fHistoGammaTrueSecondaryCaloConvFromK0sPtOrBin);

                fHistoGammaTrueSecondaryCaloUnConvFromLambdaPt        = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusGammaFromXFromLambda_Pt");
                fHistoGammaTrueSecondaryCaloUnConvFromLambdaPt->SetName("TrueSecondaryClusUnConvGammaFromLambda_Pt");
                fHistoGammaTrueSecondaryCaloUnConvFromLambdaPt->Sumw2();
                fHistoGammaTrueSecondaryCaloUnConvFromLambdaPtOrBin   = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvFromLambdaPt->Clone("TrueSecondaryClusUnConvGammaFromLambda_Pt_OriginalBinning");
                fHistoGammaTrueSecondaryCaloUnConvFromLambdaPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTrueSecondaryCaloUnConvFromLambdaPt);
                
                fHistoGammaTrueSecondaryCaloConvFromLambdaPt          = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusConvGammaFromXFromLambda_Pt");
                fHistoGammaTrueSecondaryCaloConvFromLambdaPt->Sumw2();
                fHistoGammaTrueSecondaryCaloConvFromLambdaPtOrBin     = (TH1D*)fHistoGammaTrueSecondaryCaloConvFromLambdaPt->Clone("TrueSecondaryClusConvGammaFromLambda_Pt_OriginalBinning");
                fHistoGammaTrueSecondaryCaloConvFromLambdaPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTrueSecondaryCaloConvFromLambdaPt);
                
                fHistoGammaTrueSecondaryCaloFromLambdaPt        = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvFromLambdaPt->Clone("TrueSecondaryClusGammaFromXFromLambda_Pt");
                fHistoGammaTrueSecondaryCaloFromLambdaPt->Add(fHistoGammaTrueSecondaryCaloConvFromLambdaPt);
                fHistoGammaTrueSecondaryCaloFromLambdaPtOrBin   = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvFromLambdaPtOrBin->Clone("TrueSecondaryClusGammaFromLambda_Pt_OriginalBinning");
                fHistoGammaTrueSecondaryCaloFromLambdaPtOrBin->Add(fHistoGammaTrueSecondaryCaloConvFromLambdaPtOrBin);

                // Response matrix calo
                fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC  = (TH2D*)TrueConversionContainer->FindObject("TruePrimaryClusGamma_Pt_MCPt");
                fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->SetName("TruePrimaryClusUnConvGamma_Pt_MCPt");
                fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->Sumw2();
                fHistoGammaTruePrimaryCaloUnConvMCPt            = (TH1D*)fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->ProjectionY("ESD_TruePrimaryCaloUnConvGamma_MCPt");
                RebinSpectrum(fHistoGammaTruePrimaryCaloUnConvMCPt);

                fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC    = (TH2D*)TrueConversionContainer->FindObject("TruePrimaryClusConvGamma_Pt_MCPt");
                fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->Sumw2();
                fHistoGammaTruePrimaryCaloConvMCPt              = (TH1D*)fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->ProjectionY("ESD_TruePrimaryCaloConvGamma_MCPt");
                RebinSpectrum(fHistoGammaTruePrimaryCaloConvMCPt);

                fHistoGammaTruePrimaryCalo_recPt_MCPt_MC        = (TH2D*)fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->Clone("TruePrimaryClusGamma_Pt_MCPt");
                fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->Add(fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC);
                fHistoGammaTruePrimaryCaloMCPt                  = (TH1D*)fHistoGammaTruePrimaryCaloUnConvMCPt->Clone("ESD_TruePrimaryCaloGamma_MCPt");
                fHistoGammaTruePrimaryCaloMCPt->Add(fHistoGammaTruePrimaryCaloConvMCPt);
                
                // combinatorial BG distributions
                fHistoCombinatorialBackground                   = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueClusPhotonBG_Pt");
                fHistoCombinatorialSpecies                      = new TH1D*[10];
                fHistoCombinatorialSpecies[nContamination]      = (TH1D*)fHistoCombinatorialBackground->ProjectionX(Form("ESD_TrueComb%s_Pt",fContamination[nContamination].Data()));
                fHistoCombinatorialSpecies[nContamination]->Sumw2();
                RebinSpectrum(fHistoCombinatorialSpecies[9]);
                for(Int_t i = 0; i<nContamination; i++){
                    fHistoCombinatorialSpecies[i]               = (TH1D*)fHistoCombinatorialBackground->ProjectionX(Form("ESD_TrueComb%s_Pt",fContamination[i].Data()),i+1,i+1);
                    fHistoCombinatorialSpecies[i]->Sumw2();
                    RebinSpectrum(fHistoCombinatorialSpecies[i]);
                }
            }    
        }
        
        //************************************** Calculate correction factors **************************************
        CalculateGammaCorrection();
        
        //**************************** Calculate pilepup correction factors for MC *********************************
        if(pileUpCorrection && !addSig && !(mode == 4 || mode == 5)){
            FillDCAHistogramsFromTree(dcaTree,kTRUE);
            CalculatePileUpBackground(kTRUE);
            CalculatePileUpGammaCorrection();
        }

        //**************************** Save correction histograms for gammas ***************************************
        SaveCorrectionHistos(cutSelection,fPrefix2,pileUpCorrection);
    }

    //******************************** Save raw histograms for gammas ***************************************
    SaveHistos(fIsMC,cutSelection,fPrefix2,pileUpCorrection);

    //******************************** Save pileup related histograms for raw gammas ************************
    if(pileUpCorrection && !addSig){
        SaveDCAHistos(fIsMC,cutSelection,fPrefix2);
        PlotAdditionalDCAz(fIsMC,cutSelection);
    }
}

//**************************************************************************************************
//******************************** Function to rebin spectra histos ********************************
//**************************************************************************************************
void RebinSpectrum(TH1D *Spectrum, TString NewName){
    if(NewName.CompareTo(""))
        NewName = Spectrum->GetName();

    *Spectrum = *((TH1D*)Spectrum->Rebin(fNBinsPt,NewName,fBinsPt));
    Spectrum->Divide(fDeltaPt);
}

void RebinSpectrumToDCAzDistBinning(TH1D *Spectrum, TString NewName){
    if(NewName.CompareTo(""))
        NewName = Spectrum->GetName();
    
    *Spectrum = *((TH1D*)Spectrum->Rebin(fNBinsPtDummy,NewName,fBinsPtDummy));
    Spectrum->Divide(fDeltaPtDummy);
}


//**************************************************************************************************
//******** Function to calculate BG from pileup based on DCAz distribution of converted photons ****
//******** including distribution for MC (which doesn't contain pileup) to estimate fake pileup ****
//**************************************************************************************************
void CalculatePileUpBackground(Bool_t doMC){
    
    if(fIsMC && doMC){

        fMCrecGammaPtDCAzBins                                                               = new TH1D**[4];
        fMCrecGammaPtDCAzBinsBack                                                           = new TH1D**[4];
        fMCrecSubGammaPtDCAzBins                                                            = new TH1D**[4];
        
        fTruePrimaryGammaPtDCAzBins                                                         = new TH1D**[4];
        fTruePrimarySubGammaPtDCAzBins                                                      = new TH1D**[4];
        fTrueSecondaryGammaPtDCAzBins                                                       = new TH1D**[4];
        fTrueSecondarySubGammaPtDCAzBins                                                    = new TH1D**[4];
        fTrueSecondaryGammaFromXFromK0sPtDCAzBins                                           = new TH1D**[4];
        fTrueSecondarySubGammaFromXFromK0sPtDCAzBins                                        = new TH1D**[4];
        
        fTrueBackgroundPtDCAzBins                                                           = new TH1D**[4];
        fTrueGammaPtDCAzBins                                                                = new TH1D**[4];
        fTrueSubGammaPtDCAzBins                                                             = new TH1D**[4];
        
        for (Int_t i = 0; i < 4; i++) {
            fMCrecGammaPtDCAzBins[i]                                                        = new TH1D*[fNBinsPtDummy+1];
            fMCrecGammaPtDCAzBinsBack[i]                                                    = new TH1D*[fNBinsPtDummy+1];
            fMCrecSubGammaPtDCAzBins[i]                                                     = new TH1D*[fNBinsPtDummy+1];
            
            fTruePrimaryGammaPtDCAzBins[i]                                                  = new TH1D*[fNBinsPtDummy+1];
            fTruePrimarySubGammaPtDCAzBins[i]                                               = new TH1D*[fNBinsPtDummy+1];
            fTrueSecondaryGammaPtDCAzBins[i]                                                = new TH1D*[fNBinsPtDummy+1];
            fTrueSecondarySubGammaPtDCAzBins[i]                                             = new TH1D*[fNBinsPtDummy+1];
            fTrueSecondaryGammaFromXFromK0sPtDCAzBins[i]                                    = new TH1D*[fNBinsPtDummy+1];
            fTrueSecondarySubGammaFromXFromK0sPtDCAzBins[i]                                 = new TH1D*[fNBinsPtDummy+1];
            
            fTrueBackgroundPtDCAzBins[i]                                                    = new TH1D*[fNBinsPtDummy+1];
            fTrueGammaPtDCAzBins[i]                                                         = new TH1D*[fNBinsPtDummy+1];
            fTrueSubGammaPtDCAzBins[i]                                                      = new TH1D*[fNBinsPtDummy+1];
        }
        
        fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat                            = new TH1D("MCrec_ConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_AllCatComb", "", fNBinsPtDummy, fBinsPtDummy);
        fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat->Sumw2();
        fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinning                                  = new TH1D("MCrec_ConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning", "", fNBinsPtDummy, fBinsPtDummy);
        fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinning->Sumw2();
        
        fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat                  = new TH1D("ESD_TruePrimaryConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_AllCatComb", "", fNBinsPtDummy, fBinsPtDummy);
        fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat->Sumw2();
        fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning                        = new TH1D("ESD_TruePrimaryConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning", "", fNBinsPtDummy, fBinsPtDummy);
        fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning->Sumw2();

        fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat                = new TH1D("ESD_TrueSecondaryConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_AllCatComb", "", fNBinsPtDummy, fBinsPtDummy);
        fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat->Sumw2();
        fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning                      = new TH1D("ESD_TrueSecondaryConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning", "", fNBinsPtDummy, fBinsPtDummy);
        fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning->Sumw2();

        fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat    = new TH1D("ESD_TrueSecondaryFromXFromK0sConvGamma_Pt__Ratio_WithWithoutPileUp_DCAzDistBinning_AllCatComb", "", fNBinsPtDummy, fBinsPtDummy);
        fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat->Sumw2();
        fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinning          = new TH1D("ESD_TrueSecondaryFromXFromK0sConvGamma_Pt__Ratio_WithWithoutPileUp_DCAzDistBinning", "", fNBinsPtDummy, fBinsPtDummy);
        fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinning->Sumw2();
        
        // loop over photon categories
        Int_t category;
        for (Int_t catIter = 0; catIter < 4; catIter++) {
            if (catIter == 0) {
                category            = 5;
            } else {
                category            = catIter;
            }
            
            // MCrec gamma
            fMCrecGammaPtDCAzBins[catIter][0]                                               = (TH1D*)fESDGammaPtDCAz[category]->ProjectionY(Form("MCrec_GammaPtDCAzBin_Full_%s", categoryName[catIter].Data()));
            fMCrecGammaPtDCAzBins[catIter][0]->Sumw2();
            
            // fake pileup background estimate from MCrec gamma
            if (catIter < 3){
                fMCrecGammaPtDCAzBinsBack[catIter][0]                                       = (TH1D*)fMCrecGammaPtDCAzBins[catIter][0]->ShowBackground(nIterationsShowBackground[catIter],optionShowBackground[0].Data());
                fMCrecGammaPtDCAzBinsBack[catIter][0]->Sumw2();
                fMCrecGammaPtDCAzBinsBack[catIter][0]->SetName(Form("MCrec_GammaPtDCAzBackBin_Full_%s", categoryName[catIter].Data()));
                if (fMCrecGammaPtDCAzBinsBack[catIter][0]->Integral() < 1 || fMCrecGammaPtDCAzBinsBack[catIter][0]->GetEntries() > fMCrecGammaPtDCAzBins[catIter][0]->GetEntries()) {
                    fMCrecGammaPtDCAzBinsBack[catIter][0]->Reset("ICES");
                }
            } else {
                fMCrecGammaPtDCAzBinsBack[catIter][0] = (TH1D*)fMCrecGammaPtDCAzBins[catIter][0]->Clone(Form("MCrec_GammaPtDCAzBackBin_Full_%s", categoryName[catIter].Data()));
                fMCrecGammaPtDCAzBinsBack[catIter][0]->Reset();
            }
            
            // true primary
            fTruePrimaryGammaPtDCAzBins[catIter][0]                                         = (TH1D*)fTruePrimaryPhotonPtDCAz[category]->ProjectionY(Form("ESD_TruePrimaryGammaPtDCAzBin_Full_%s", categoryName[catIter].Data()));
            fTruePrimaryGammaPtDCAzBins[catIter][0]->Sumw2();
            fTruePrimarySubGammaPtDCAzBins[catIter][0]                                      = (TH1D*)fTruePrimaryGammaPtDCAzBins[catIter][0]->Clone(Form("ESD_TruePrimarySubGammaPtDCAzBin_Full_%s", categoryName[catIter].Data()));
            fTruePrimarySubGammaPtDCAzBins[catIter][0]->Sumw2();
            
            // true secondary
            fTrueSecondaryGammaPtDCAzBins[catIter][0]                                       = (TH1D*)fTrueSecondaryPhotonPtDCAz[category]->ProjectionY(Form("ESD_TrueSecondaryGammaPtDCAzBin_Full_%s", categoryName[catIter].Data()));
            fTrueSecondaryGammaPtDCAzBins[catIter][0]->Sumw2();
            fTrueSecondarySubGammaPtDCAzBins[catIter][0]                                    = (TH1D*)fTrueSecondaryGammaPtDCAzBins[catIter][0]->Clone(Form("ESD_TrueSecondarySubGammaPtDCAzBin_Full_%s", categoryName[catIter].Data()));
            fTrueSecondarySubGammaPtDCAzBins[catIter][0]->Sumw2();
            
            // true secondary from X from K0s
            fTrueSecondaryGammaFromXFromK0sPtDCAzBins[catIter][0]                           = (TH1D*)fTrueSecondaryPhotonFromXFromK0sPtDCAz[category]->ProjectionY(Form("ESD_TrueSecondaryGammaFromXFromK0sPtDCAzBin_Full_%s", categoryName[catIter].Data()));
            fTrueSecondaryGammaFromXFromK0sPtDCAzBins[catIter][0]->Sumw2();
            fTrueSecondarySubGammaFromXFromK0sPtDCAzBins[catIter][0]                        = (TH1D*)fTrueSecondaryGammaPtDCAzBins[catIter][0]->Clone(Form("ESD_TrueSecondarySubGammaFromXFromK0sPtDCAzBin_Full_%s", categoryName[catIter].Data()));
            fTrueSecondarySubGammaFromXFromK0sPtDCAzBins[catIter][0]->Sumw2();
            
            // true background  = MCrec - true primary -  true secondary
            fTrueBackgroundPtDCAzBins[catIter][0]                                           = (TH1D*)fMCrecGammaPtDCAzBins[catIter][0]->Clone(Form("ESD_TrueGammaBackgroundPtDCAzBin_Full_%s", categoryName[catIter].Data()));
            fTrueBackgroundPtDCAzBins[catIter][0]->Sumw2();
            fTrueBackgroundPtDCAzBins[catIter][0]->Add(fTruePrimaryGammaPtDCAzBins[catIter][0],-1);
            fTrueBackgroundPtDCAzBins[catIter][0]->Add(fTrueSecondaryGammaPtDCAzBins[catIter][0],-1);
            
            // true gamma = true primary + true secondary
            fTrueGammaPtDCAzBins[catIter][0]                                                = (TH1D*)fTruePrimaryGammaPtDCAzBins[catIter][0]->Clone(Form("ESD_TrueGammaPtDCAz_Full_%s", categoryName[catIter].Data()));
            fTrueGammaPtDCAzBins[catIter][0]->Sumw2();
            fTrueGammaPtDCAzBins[catIter][0]->Add(fTrueSecondaryGammaPtDCAzBins[catIter][0],1);
            
            // MCrec fake pileup subtracted = MCrec - fake pileup
            fMCrecSubGammaPtDCAzBins[catIter][0]                                            = (TH1D*)fMCrecGammaPtDCAzBins[catIter][0]->Clone(Form("MCrec_SubGammaPtDCAzBin_Full_%s", categoryName[catIter].Data()));
            fMCrecSubGammaPtDCAzBins[catIter][0]->Sumw2();
            fMCrecSubGammaPtDCAzBins[catIter][0]->Add(fMCrecGammaPtDCAzBinsBack[catIter][0],-1);
            
            // true gamma fake pileup subtracted = true gamma - fake pileup
            fTrueSubGammaPtDCAzBins[catIter][0]                                             = (TH1D*)fTrueGammaPtDCAzBins[catIter][0]->Clone(Form("ESD_TrueSubGammaPtDCAz_Full_%s", categoryName[catIter].Data()));
            fTrueSubGammaPtDCAzBins[catIter][0]->Sumw2();
            fTrueSubGammaPtDCAzBins[catIter][0]->Add(fMCrecGammaPtDCAzBinsBack[catIter][0],-1);
            
            // fake pileup subtracted true primary, true secondary, true secondary from X from K0s
            for(Int_t i = 0; i<fTrueGammaPtDCAzBins[catIter][0]->GetNbinsX();i++){
                if(fTrueSubGammaPtDCAzBins[catIter][0]->GetBinContent(i+1)<0){
                    fTrueSubGammaPtDCAzBins[catIter][0]->SetBinContent(i+1,0);
                    fTrueSubGammaPtDCAzBins[catIter][0]->SetBinError(i+1,0);
                }
            }
                
            CalculatePileUpSubtractedDCAz(fTrueGammaPtDCAzBins[catIter][0], fTrueSubGammaPtDCAzBins[catIter][0], fTruePrimaryGammaPtDCAzBins[catIter][0], fTruePrimarySubGammaPtDCAzBins[catIter][0]);
            CalculatePileUpSubtractedDCAz(fTrueGammaPtDCAzBins[catIter][0], fTrueSubGammaPtDCAzBins[catIter][0], fTrueSecondaryGammaPtDCAzBins[catIter][0], fTrueSecondarySubGammaPtDCAzBins[catIter][0]);
            CalculatePileUpSubtractedDCAz(fTrueGammaPtDCAzBins[catIter][0], fTrueSubGammaPtDCAzBins[catIter][0], fTrueSecondaryGammaFromXFromK0sPtDCAzBins[catIter][0], fTrueSecondarySubGammaFromXFromK0sPtDCAzBins[catIter][0]);
            
            // loop over pt bins
            for(Int_t bin = 1; bin<fNBinsPtDummy+1; bin++){
                Int_t startBin                                                              = fHistoGammaConvPtOrBin->FindBin(fBinsPtDummy[bin-1]+0.001);
                Int_t endBin                                                                = fHistoGammaConvPtOrBin->FindBin(fBinsPtDummy[bin]-0.001);
                
                // MCrec gamma
                fMCrecGammaPtDCAzBins[catIter][bin]                                         = (TH1D*)fESDGammaPtDCAz[category]->ProjectionY(Form("MCrec_GammaPtDCAzBin_%.1f_%.1f_%s",fBinsPtDummy[bin-1],fBinsPtDummy[bin], categoryName[catIter].Data()),startBin,endBin);
                fMCrecGammaPtDCAzBins[catIter][bin]->Sumw2();
                
                // fake pileup background estimate from MCrec gamma
                if (catIter < 3){
                    fMCrecGammaPtDCAzBinsBack[catIter][bin]                                 = (TH1D*)fMCrecGammaPtDCAzBins[catIter][bin]->ShowBackground(nIterationsShowBackground[catIter],optionShowBackground[0].Data());
                    fMCrecGammaPtDCAzBinsBack[catIter][bin]->Sumw2();
                    fMCrecGammaPtDCAzBinsBack[catIter][bin]->SetName(Form("MCrec_GammaPtDCAzBackBin_%.1f_%.1f_%s", fBinsPtDummy[bin-1], fBinsPtDummy[bin], categoryName[catIter].Data()));
                    if (fMCrecGammaPtDCAzBinsBack[catIter][bin]->Integral() < 1 || fMCrecGammaPtDCAzBinsBack[catIter][bin]->GetEntries() > fMCrecGammaPtDCAzBins[catIter][bin]->GetEntries()) {
                        fMCrecGammaPtDCAzBinsBack[catIter][bin]->Reset("ICES");
                    }
                } else {
                    fMCrecGammaPtDCAzBinsBack[catIter][bin] = (TH1D*)fMCrecGammaPtDCAzBins[catIter][bin]->Clone(Form("MCrec_GammaPtDCAzBackBin_%.1f_%.1f_%s", fBinsPtDummy[bin-1], fBinsPtDummy[bin], categoryName[catIter].Data()));
                    fMCrecGammaPtDCAzBinsBack[catIter][bin]->Reset();
                }
                
                // true primary
                fTruePrimaryGammaPtDCAzBins[catIter][bin]                                   = (TH1D*)fTruePrimaryPhotonPtDCAz[category]->ProjectionY(Form("ESD_TruePrimaryGammaPtDCAzBin_%.1f_%.1f_%s",fBinsPtDummy[bin-1],fBinsPtDummy[bin], categoryName[catIter].Data()),startBin,endBin);
                fTruePrimaryGammaPtDCAzBins[catIter][bin]->Sumw2();
                fTruePrimarySubGammaPtDCAzBins[catIter][bin]                                = (TH1D*)fTruePrimaryGammaPtDCAzBins[catIter][bin]->Clone(Form("ESD_TruePrimarySubGammaPtDCAzBin_%.1f_%.1f_%s",fBinsPtDummy[bin-1],fBinsPtDummy[bin], categoryName[catIter].Data()));
                fTruePrimarySubGammaPtDCAzBins[catIter][bin]->Sumw2();
                
                // true secondary
                fTrueSecondaryGammaPtDCAzBins[catIter][bin]                                 = (TH1D*)fTrueSecondaryPhotonPtDCAz[category]->ProjectionY(Form("ESD_TrueSecondaryGammaPtDCAzBin_%.1f_%.1f_%s",fBinsPtDummy[bin-1],fBinsPtDummy[bin], categoryName[catIter].Data()),startBin,endBin);
                fTrueSecondaryGammaPtDCAzBins[catIter][bin]->Sumw2();
                fTrueSecondarySubGammaPtDCAzBins[catIter][bin]                              = (TH1D*)fTrueSecondaryGammaPtDCAzBins[catIter][bin]->Clone(Form("ESD_TrueSecondarySubGammaPtDCAzBin_%.1f_%.1f_%s",fBinsPtDummy[bin-1],fBinsPtDummy[bin], categoryName[catIter].Data()));
                fTrueSecondarySubGammaPtDCAzBins[catIter][bin]->Sumw2();
                
                // true secondary from X from K0s
                fTrueSecondaryGammaFromXFromK0sPtDCAzBins[catIter][bin]                     = (TH1D*)fTrueSecondaryPhotonFromXFromK0sPtDCAz[category]->ProjectionY(Form("ESD_TrueSecondaryGammaFromXFromK0sPtDCAzBin_%.1f_%.1f_%s",fBinsPtDummy[bin-1],fBinsPtDummy[bin], categoryName[catIter].Data()),startBin,endBin);
                fTrueSecondaryGammaFromXFromK0sPtDCAzBins[catIter][bin]->Sumw2();
                fTrueSecondarySubGammaFromXFromK0sPtDCAzBins[catIter][bin]                  = (TH1D*)fTrueSecondaryGammaFromXFromK0sPtDCAzBins[catIter][bin]->Clone(Form("ESD_TrueSecondarySubGammaFromXFromK0sPtDCAzBin_%.1f_%.1f_%s",fBinsPtDummy[bin-1],fBinsPtDummy[bin], categoryName[catIter].Data()));
                fTrueSecondarySubGammaFromXFromK0sPtDCAzBins[catIter][bin]->Sumw2();
                
                // true background  = MCrec - true primary -  true secondary
                fTrueBackgroundPtDCAzBins[catIter][bin]                                     = (TH1D*)fMCrecGammaPtDCAzBins[catIter][bin]->Clone(Form("ESD_TrueGammaBackgroundPtDCAzBin_%.1f_%.1f_%s",fBinsPtDummy[bin-1],fBinsPtDummy[bin], categoryName[catIter].Data()));
                fTrueBackgroundPtDCAzBins[catIter][bin]->Sumw2();
                fTrueBackgroundPtDCAzBins[catIter][bin]->Add(fTruePrimaryGammaPtDCAzBins[catIter][bin],-1);
                fTrueBackgroundPtDCAzBins[catIter][bin]->Add(fTrueSecondaryGammaPtDCAzBins[catIter][bin],-1);
                
                // true gamma = true primary + true secondary
                fTrueGammaPtDCAzBins[catIter][bin]                                      = (TH1D*)fTruePrimaryGammaPtDCAzBins[catIter][bin]->Clone(Form("ESD_TrueGammaPtDCAz_%.1f_%.1f_%s",fBinsPtDummy[bin-1],fBinsPtDummy[bin], categoryName[catIter].Data()));
                fTrueGammaPtDCAzBins[catIter][bin]->Sumw2();
                fTrueGammaPtDCAzBins[catIter][bin]->Add(fTrueSecondaryGammaPtDCAzBins[catIter][bin],1);

                // MCrec fake pileup subtracted = MCrec - fake pileup
                fMCrecSubGammaPtDCAzBins[catIter][bin]                                  = (TH1D*)fMCrecGammaPtDCAzBins[catIter][bin]->Clone(Form("MCrec_SubGammaPtDCAzBin_%.1f_%.1f_%s",fBinsPtDummy[bin-1],fBinsPtDummy[bin], categoryName[catIter].Data()));
                fMCrecSubGammaPtDCAzBins[catIter][bin]->Sumw2();
                fMCrecSubGammaPtDCAzBins[catIter][bin]->Add(fMCrecGammaPtDCAzBinsBack[catIter][bin],-1);
                
                // true gamma fake pileup subtracted = true gamma - fake pileup
                fTrueSubGammaPtDCAzBins[catIter][bin]                                   = (TH1D*)fTrueGammaPtDCAzBins[catIter][bin]->Clone(Form("ESD_TrueSubGammaPtDCAz_%.1f_%.1f_%s",fBinsPtDummy[bin-1],fBinsPtDummy[bin], categoryName[catIter].Data()));
                fTrueSubGammaPtDCAzBins[catIter][bin]->Sumw2();
                fTrueSubGammaPtDCAzBins[catIter][bin]->Add(fMCrecGammaPtDCAzBinsBack[catIter][bin],-1);
                
                // fake pileup subtracted true primary, true secondary, true secondary from X from K0s
                for(Int_t i = 0; i<fTrueGammaPtDCAzBins[catIter][bin]->GetNbinsX();i++){
                    if(fTrueSubGammaPtDCAzBins[catIter][bin]->GetBinContent(i+1)<0){
                        fTrueSubGammaPtDCAzBins[catIter][bin]->SetBinContent(i+1,0);
                        fTrueSubGammaPtDCAzBins[catIter][bin]->SetBinError(i+1,0);
                    }
                }
                
                CalculatePileUpSubtractedDCAz(fTrueGammaPtDCAzBins[catIter][bin], fTrueSubGammaPtDCAzBins[catIter][bin], fTruePrimaryGammaPtDCAzBins[catIter][bin], fTruePrimarySubGammaPtDCAzBins[catIter][bin]);
                CalculatePileUpSubtractedDCAz(fTrueGammaPtDCAzBins[catIter][bin], fTrueSubGammaPtDCAzBins[catIter][bin], fTrueSecondaryGammaPtDCAzBins[catIter][bin], fTrueSecondarySubGammaPtDCAzBins[catIter][bin]);
                CalculatePileUpSubtractedDCAz(fTrueGammaPtDCAzBins[catIter][bin], fTrueSubGammaPtDCAzBins[catIter][bin], fTrueSecondaryGammaFromXFromK0sPtDCAzBins[catIter][bin], fTrueSecondarySubGammaFromXFromK0sPtDCAzBins[catIter][bin]);
            }
            
            //plotting DCAz distributions for MC rec and identified particle in pt slices
            TString nameFile    = Form("%s/%s_%s_MCrec_DCAz_vs_Pt_%s_%s.%s", fOutputDir.Data(), fPrefix.Data(), fPrefix2.Data(), categoryName[catIter].Data(), fCutSelection.Data(), fSuffix.Data());
            PlotDCAzInPtBinsWithBack( fMCrecGammaPtDCAzBins[catIter], fMCrecGammaPtDCAzBinsBack[catIter], NULL, nameFile, "CanvasESDDCAz", "PadESDDCAz",
            fDate, fMeson, 1, fNBinsPtDummy, fBinsPtDummy, "#gamma --> e^{+}e^{-}", fIsMC, "MinBias");
            
            nameFile            = Form("%s/%s_%s_SignalAfterSubtraction_DCAz_vs_Pt_%s_%s.%s", fOutputDir.Data(), fPrefix.Data(), fPrefix2.Data(), categoryName[catIter].Data(), fCutSelection.Data(), fSuffix.Data());
            PlotDCAzInPtBinsWithBack( fMCrecSubGammaPtDCAzBins[catIter], fTrueSubGammaPtDCAzBins[catIter], NULL, nameFile, "CanvasESDDCAz", "PadESDDCAz",
            fDate, fMeson, 1, fNBinsPtDummy, fBinsPtDummy, "#gamma --> e^{+}e^{-}", fIsMC, "MinBias");
            
            nameFile            = Form("%s/%s_%s_TrueBackDCAz_vs_Pt_%s_%s.%s", fOutputDir.Data(), fPrefix.Data(), fPrefix2.Data(), categoryName[catIter].Data(), fCutSelection.Data(), fSuffix.Data());
            PlotDCAzInPtBinsWithBack( fTrueBackgroundPtDCAzBins[catIter], fMCrecGammaPtDCAzBinsBack[catIter], NULL, nameFile, "CanvasESDDCAz", "PadESDDCAz",
            fDate, fMeson, 1, fNBinsPtDummy, fBinsPtDummy, "#gamma --> e^{+}e^{-}", fIsMC, "MinBias");
             
            nameFile            = Form("%s/%s_%s_TrueSignalDCAz_vs_Pt_%s_%s.%s", fOutputDir.Data(), fPrefix.Data(), fPrefix2.Data(), categoryName[catIter].Data(), fCutSelection.Data(), fSuffix.Data());
            PlotDCAzInPtBinsWithBack( fMCrecGammaPtDCAzBins[catIter], fTrueGammaPtDCAzBins[catIter], NULL, nameFile, "CanvasESDDCAz", "PadESDDCAz",
            fDate, fMeson, 1, fNBinsPtDummy, fBinsPtDummy, "#gamma --> e^{+}e^{-}", fIsMC,  "MinBias");
        }
        
        // building ratios with / without fake pileup
        CalculateDCAzDistributionRatio(fMCrecGammaPtDCAzBins, fMCrecSubGammaPtDCAzBins, 0, 0, fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat);
        CalculateDCAzDistributionRatio(fMCrecGammaPtDCAzBins, fMCrecSubGammaPtDCAzBins, 1, 3, fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinning);
        
        CalculateDCAzDistributionRatio(fTruePrimaryGammaPtDCAzBins, fTruePrimarySubGammaPtDCAzBins, 0, 0, fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat);
        CalculateDCAzDistributionRatio(fTruePrimaryGammaPtDCAzBins, fTruePrimarySubGammaPtDCAzBins, 1, 3, fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning);

        CalculateDCAzDistributionRatio(fTrueSecondaryGammaPtDCAzBins, fTrueSecondarySubGammaPtDCAzBins, 0, 0, fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat);
        CalculateDCAzDistributionRatio(fTrueSecondaryGammaPtDCAzBins, fTrueSecondarySubGammaPtDCAzBins, 1, 3, fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning);

        CalculateDCAzDistributionRatio(fTrueSecondaryGammaFromXFromK0sPtDCAzBins, fTrueSecondarySubGammaFromXFromK0sPtDCAzBins, 0, 0, fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat);
        CalculateDCAzDistributionRatio(fTrueSecondaryGammaFromXFromK0sPtDCAzBins, fTrueSecondarySubGammaFromXFromK0sPtDCAzBins, 1, 3, fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinning);
        
        // define pileup correction factors
        fMCrecGammaPileUpCorrFactorAllCat                                           = new TH1D("fMCrecGammaPileUpCorrFactorAllCatComb", "fMCrecGammaPileUpCorrFactorAllCatComb", fNBinsPt, fBinsPt);
        fMCrecGammaPileUpCorrFactorAllCat->Sumw2();
        fMCrecGammaPileUpCorrFactor                                                 = new TH1D("fMCrecGammaPileUpCorrFactor", "fMCrecGammaPileUpCorrFactor", fNBinsPt, fBinsPt);
        fMCrecGammaPileUpCorrFactor->Sumw2();
        fTruePrimaryConvGammaPileUpCorrFactorAllCat                                 = new TH1D("fTruePrimaryConvGammaPileUpCorrFactorAllCatComb", "fTruePrimaryConvGammaPileUpCorrFactorAllCatComb", fNBinsPt, fBinsPt);
        fTruePrimaryConvGammaPileUpCorrFactorAllCat->Sumw2();
        fTruePrimaryConvGammaPileUpCorrFactor                                       = new TH1D("fTruePrimaryConvGammaPileUpCorrFactor", "fTruePrimaryConvGammaPileUpCorrFactor", fNBinsPt, fBinsPt);
        fTruePrimaryConvGammaPileUpCorrFactor->Sumw2();
        fTrueSecondaryConvGammaPileUpCorrFactorAllCat                               = new TH1D("fTrueSecondaryConvGammaPileUpCorrFactorAllCatComb", "fTrueSecondaryConvGammaPileUpCorrFactorAllCatComb", fNBinsPt, fBinsPt);
        fTrueSecondaryConvGammaPileUpCorrFactorAllCat->Sumw2();
        fTrueSecondaryConvGammaPileUpCorrFactor                                     = new TH1D("fTrueSecondaryConvGammaPileUpCorrFactor", "fTrueSecondaryConvGammaPileUpCorrFactor", fNBinsPt, fBinsPt);
        fTrueSecondaryConvGammaPileUpCorrFactor->Sumw2();
        fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactorAllCat                   = new TH1D("fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactorAllCatComb", "fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactorAllCatComb", fNBinsPt, fBinsPt);
        fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactorAllCat->Sumw2();
        fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactor                         = new TH1D("fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactor", "fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactor", fNBinsPt, fBinsPt);
        fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactor->Sumw2();
        
        // calculating pileup correction factors
        CalculatePileUpCorrectionFactor(fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat, fMCrecGammaPileUpCorrFactorAllCat, fMCrecGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat);
        CalculatePileUpCorrectionFactor(fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinning, fMCrecGammaPileUpCorrFactor, fMCrecGammaPtRatioWithWithoutPileUpFitDCAzDistBinning);
        
        CalculatePileUpCorrectionFactor(fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat, fTruePrimaryConvGammaPileUpCorrFactorAllCat, fTruePrimaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat);
        CalculatePileUpCorrectionFactor(fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning, fTruePrimaryConvGammaPileUpCorrFactor, fTruePrimaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning);
        
        CalculatePileUpCorrectionFactor(fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat, fTrueSecondaryConvGammaPileUpCorrFactorAllCat, fTrueSecondaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat);
        CalculatePileUpCorrectionFactor(fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning, fTrueSecondaryConvGammaPileUpCorrFactor, fTrueSecondaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning);
        
        CalculatePileUpCorrectionFactor(fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat, fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactorAllCat, fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat);
        CalculatePileUpCorrectionFactor(fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinning, fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactor, fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning);
        
        // calculate spectra w/o fake pileup
        fMCrecGammaPtPileUpAllCat                                                   = (TH1D*)fHistoGammaMCrecConvPt->Clone("MCrec_ConvGamma_Pt_PileUp_AllCatComb");
        fMCrecGammaPtPileUpAllCat->Sumw2();
        fMCrecGammaPtPileUpAllCat->Multiply(fMCrecGammaPileUpCorrFactorAllCat);
        
        fMCrecGammaPtPileUp                                                         = (TH1D*)fHistoGammaMCrecConvPt->Clone("MCrec_ConvGamma_Pt_PileUp");
        fMCrecGammaPtPileUp->Sumw2();
        fMCrecGammaPtPileUp->Multiply(fMCrecGammaPileUpCorrFactor);
        
        fTruePrimaryConvGammaPtPileUpAllCat                                         = (TH1D*)fHistoGammaTruePrimaryConvPt->Clone("ESD_TruePrimaryConvGamma_Pt_PileUp_AllCatComb");
        fTruePrimaryConvGammaPtPileUpAllCat->Sumw2();
        fTruePrimaryConvGammaPtPileUpAllCat->Multiply(fTruePrimaryConvGammaPileUpCorrFactorAllCat);
        
        fTruePrimaryConvGammaPtPileUp                                               = (TH1D*)fHistoGammaTruePrimaryConvPt->Clone("ESD_TruePrimaryConvGamma_Pt_PileUp");
        fTruePrimaryConvGammaPtPileUp->Sumw2();
        fTruePrimaryConvGammaPtPileUp->Multiply(fTruePrimaryConvGammaPileUpCorrFactor);

        fTrueSecondaryConvGammaPtPileUpAllCat                                       = (TH1D*)fHistoGammaTrueSecondaryConvPt->Clone("ESD_TrueSecondaryConvGamma_Pt_PileUp_AllCatComb");
        fTrueSecondaryConvGammaPtPileUpAllCat->Sumw2();
        fTrueSecondaryConvGammaPtPileUpAllCat->Multiply(fTrueSecondaryConvGammaPileUpCorrFactorAllCat);
        
        fTrueSecondaryConvGammaPtPileUp                                             = (TH1D*)fHistoGammaTrueSecondaryConvPt->Clone("ESD_TrueSecondaryConvGamma_Pt_PileUp");
        fTrueSecondaryConvGammaPtPileUp->Sumw2();
        fTrueSecondaryConvGammaPtPileUp->Multiply(fTrueSecondaryConvGammaPileUpCorrFactor);

        fTrueSecondaryFromXFromK0sConvGammaPtPileUpAllCat                           = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt->Clone("ESD_TrueSecondaryConvGammaFromXFromK0s_Pt_PileUp_AllCatComb");
        fTrueSecondaryFromXFromK0sConvGammaPtPileUpAllCat->Sumw2();
        fTrueSecondaryFromXFromK0sConvGammaPtPileUpAllCat->Multiply(fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactorAllCat);
        
        fTrueSecondaryFromXFromK0sConvGammaPtPileUp                                 = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt->Clone("ESD_TrueSecondaryConvGammaFromXFromK0s_Pt_PileUp");
        fTrueSecondaryFromXFromK0sConvGammaPtPileUp->Sumw2();
        fTrueSecondaryFromXFromK0sConvGammaPtPileUp->Multiply(fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactor);
        
        // plotting ratios + fits
        TCanvas *RatioWithWithoutPileUpCanvasMC                                     = GetAndSetCanvas("canvasRatioWithWithoutPileUpMC");
        
        SetHistogramm(fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinning,"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.25);
        SetHistogramm(fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning,"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.25);
        SetHistogramm(fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning,"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.25);
        SetHistogramm(fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinning,"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.25);
        
        DrawGammaSetMarker(fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinning, 24, 1.0, kBlack, kBlack);
        DrawGammaSetMarker(fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning, 24, 1.0, kRed, kRed);
        DrawGammaSetMarker(fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning, 24, 1.0, kBlue-9, kBlue-9);
        DrawGammaSetMarker(fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinning, 24, 1.0, kBlue+3, kBlue+3);

        if (fMCrecGammaPtRatioWithWithoutPileUpFitDCAzDistBinning) fMCrecGammaPtRatioWithWithoutPileUpFitDCAzDistBinning->SetLineColor(kBlack);
        if (fTruePrimaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning) fTruePrimaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning->SetLineColor(kRed);
        if (fTrueSecondaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning) fTrueSecondaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning->SetLineColor(kBlue-9);
        if (fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning) fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning->SetLineColor(kBlue+3);
        
        fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinning->DrawCopy();
        if (fMCrecGammaPtRatioWithWithoutPileUpFitDCAzDistBinning) fMCrecGammaPtRatioWithWithoutPileUpFitDCAzDistBinning->Draw("same");
        
        fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning->DrawCopy("same");
        if (fTruePrimaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning) fTruePrimaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning->Draw("same");
        
        fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning->DrawCopy("same");
        if (fTrueSecondaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning) fTrueSecondaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning->Draw("same");
        
        fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinning->DrawCopy("same");
        if (fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning) fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning->Draw("same");
        
        TLegend* legend     = GetAndSetLegend(0.6,0.75,4,1);
        legend->AddEntry(fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinning,"rec. #gamma","lp");
        legend->AddEntry(fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning,"true prim. #gamma","lp");
        legend->AddEntry(fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning,"true sec. #gamma","lp");
        legend->AddEntry(fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinning,"true sec. #gamma from X from K^{0}_{s}","lp");
        legend->Draw("same");
        
        RatioWithWithoutPileUpCanvasMC->Print(Form("%s/%s_%s_With_vs_Without_Pileup_pT_%s.%s",fOutputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fCutSelection.Data(),fSuffix.Data()));
        delete RatioWithWithoutPileUpCanvasMC;
    } else {
        // *****************************************************************************************************
        // ************************* Processing pileup estimation based on DCAz for Data ***********************
        // *****************************************************************************************************
        // create histos with dca z distribution for each bin
        fESDGammaPtDCAzBins                                                 = new TH1D**[4];
        fESDGammaPtDCAzBinsBack                                             = new TH1D***[4];
        fESDSubGammaPtDCAzBins                                              = new TH1D***[4];
        
        for (Int_t i = 0; i < 4; i++) {
            fESDGammaPtDCAzBins[i]                                          = new TH1D*[fNBinsPtDummy+1];
            fESDGammaPtDCAzBinsBack[i]                                      = new TH1D**[fNBinsPtDummy+1];
            fESDSubGammaPtDCAzBins[i]                                       = new TH1D**[fNBinsPtDummy+1];
            
            for (Int_t j = 0; j < fNBinsPtDummy+1; j++) {
                fESDGammaPtDCAzBinsBack[i][j]                               = new TH1D*[3];
                fESDSubGammaPtDCAzBins[i][j]                                = new TH1D*[3];
            }
        }
        
        fESDGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat              = new TH1D*[3];
        fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning                    = new TH1D*[3];
        
        for (Int_t i = 0; i < 3; i++) {
            fESDGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat[i]       = new TH1D(Form("ESD_ConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_AllCatComb_%s", backgroundExtractionMethod[i].Data()), "", fNBinsPtDummy, fBinsPtDummy);
            fESDGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat[i]->Sumw2();
            fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[i]             = new TH1D(Form("ESD_ConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_%s", backgroundExtractionMethod[i].Data()), "", fNBinsPtDummy, fBinsPtDummy);
            fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[i]->Sumw2();
        }
        
        fESDGammaPerCatPtDCAzBins                                           = new TH1D*[4];
        fESDGammaRatioCatToCombinedPtDCAzBins                               = new TH1D*[3];

        // loop over photon categories
        Int_t category;
        for (Int_t catIter = 0; catIter < 4; catIter++) {
            if (catIter == 0) {
                category = 5;
            } else {
                category = catIter;
            }
            
            fESDGammaPtDCAzBins[catIter][0]                                 = (TH1D*)fESDGammaPtDCAz[category]->ProjectionY(Form("ESD_GammaPtDCAzBin_Full_%s", categoryName[catIter].Data()));
            fESDGammaPtDCAzBins[catIter][0]->Sumw2();
            
            fESDGammaPerCatPtDCAzBins[catIter]                              = new TH1D(Form("ESD_GammaPtDCAzBin_%s", categoryName[catIter].Data()), "", fNBinsPtDummy, fBinsPtDummy);
            fESDGammaPerCatPtDCAzBins[catIter]->Sumw2();
            
            // loop over background extraction methods
            for (Int_t i = 0; i < 3; i++) {

                // estimate pileup BG
                if (catIter < 3) {
                    fESDGammaPtDCAzBinsBack[catIter][0][i]                  = (TH1D*)fESDGammaPtDCAzBins[catIter][0]->ShowBackground(nIterationsShowBackground[catIter],optionShowBackground[i].Data());
                    fESDGammaPtDCAzBinsBack[catIter][0][i]->SetName(Form("ESD_GammaPtDCAzBackBin_Full_%s_%s", categoryName[catIter].Data(), backgroundExtractionMethod[i].Data()));
                    fESDGammaPtDCAzBinsBack[catIter][0][i]->Sumw2();
                    if (fESDGammaPtDCAzBinsBack[catIter][0][i]->Integral() < 1 || fESDGammaPtDCAzBinsBack[catIter][0][i]->GetEntries() > fESDGammaPtDCAzBins[catIter][0]->GetEntries()) {
                        fESDGammaPtDCAzBinsBack[catIter][0][i]->Reset("ICES");
                    }
                } else {
                    fESDGammaPtDCAzBinsBack[catIter][0][i]                  = (TH1D*)fESDGammaPtDCAzBins[catIter][0]->Clone(Form("ESD_GammaPtDCAzBackBin_Full_%s_%s", categoryName[catIter].Data(), backgroundExtractionMethod[i].Data()));
                    fESDGammaPtDCAzBinsBack[catIter][0][i]->Reset();
                }
                
                // subtract estimated pileup BG
                fESDSubGammaPtDCAzBins[catIter][0][i]                       = (TH1D*)fESDGammaPtDCAzBins[catIter][0]->Clone(Form("ESD_SubGammaPtDCAzBin_Full_%s_%s", categoryName[catIter].Data(), backgroundExtractionMethod[i].Data()));
                fESDSubGammaPtDCAzBins[catIter][0][i]->Sumw2();
                fESDSubGammaPtDCAzBins[catIter][0][i]->Add(fESDGammaPtDCAzBinsBack[catIter][0][i],-1);
            }

            // loop over pt bins
            for(Int_t bin = 1; bin<fNBinsPtDummy+1; bin++){
                Int_t startBin                                              = fHistoGammaConvPtOrBin->FindBin(fBinsPtDummy[bin-1]+0.001);
                Int_t endBin                                                = fHistoGammaConvPtOrBin->FindBin(fBinsPtDummy[bin]-0.001);
                
                fESDGammaPtDCAzBins[catIter][bin]                           = (TH1D*)fESDGammaPtDCAz[category]->ProjectionY(Form("ESD_GammaPtDCAzBin_%.1f_%.1f_%s",fBinsPtDummy[bin-1],fBinsPtDummy[bin], categoryName[catIter].Data()),startBin,endBin);
                fESDGammaPtDCAzBins[catIter][bin]->Sumw2();
                
                // raw yields per category
                Double_t tempBinError                                       = 0;
                Double_t tempBinContent                                     = fESDGammaPtDCAzBins[catIter][bin]->IntegralAndError(-1000,1000,tempBinError);
                fESDGammaPerCatPtDCAzBins[catIter]->SetBinContent(bin,      tempBinContent);
                fESDGammaPerCatPtDCAzBins[catIter]->SetBinError(bin,        tempBinError);

                // loop over background extraction methods
                for (Int_t i = 0; i < 3; i++) {
                    
                    // estimate pileup BG
                    if (catIter < 3) {
                        fESDGammaPtDCAzBinsBack[catIter][bin][i]            = (TH1D*)fESDGammaPtDCAzBins[catIter][bin]->ShowBackground(nIterationsShowBackground[catIter],optionShowBackground[i].Data());
                        fESDGammaPtDCAzBinsBack[catIter][bin][i]->SetName(Form("ESD_GammaPtDCAzBackBin_%.1f_%.1f_%s_%s", fBinsPtDummy[bin-1], fBinsPtDummy[bin], categoryName[catIter].Data(), backgroundExtractionMethod[i].Data()));
                        fESDGammaPtDCAzBinsBack[catIter][bin][i]->Sumw2();
                        if (fESDGammaPtDCAzBinsBack[catIter][bin][i]->Integral() < 1 || fESDGammaPtDCAzBinsBack[catIter][bin][i]->GetEntries() > fESDGammaPtDCAzBins[catIter][bin]->GetEntries()) {
                            fESDGammaPtDCAzBinsBack[catIter][bin][i]->Reset("ICES");
                        }
                    } else {
                        fESDGammaPtDCAzBinsBack[catIter][bin][i]            = (TH1D*)fESDGammaPtDCAzBins[catIter][bin]->Clone(Form("ESD_GammaPtDCAzBackBin_%.1f_%0.1f_%s_%s", fBinsPtDummy[bin-1],fBinsPtDummy[bin], categoryName[catIter].Data(), backgroundExtractionMethod[i].Data()));
                        fESDGammaPtDCAzBinsBack[catIter][bin][i]->Reset();
                    }
                    
                    // subtract estimated pileup BG
                    fESDSubGammaPtDCAzBins[catIter][bin][i]                 = (TH1D*)fESDGammaPtDCAzBins[catIter][bin]->Clone(Form("ESD_SubGammaPtDCAzBin_%.1f_%.1f_%s_%s",fBinsPtDummy[bin-1],fBinsPtDummy[bin], categoryName[catIter].Data(), backgroundExtractionMethod[i].Data()));
                    fESDSubGammaPtDCAzBins[catIter][bin][i]->Sumw2();
                    fESDSubGammaPtDCAzBins[catIter][bin][i]->Add(fESDGammaPtDCAzBinsBack[catIter][bin][i],-1);
                }
            }
            
            // plotting DCAz distributions for rec gamma with estimated BG in pt slices
            TString nameFile                                                = Form("%s/%s_%s_ESD_DCAz_vs_Pt_%s_%s.%s", fOutputDir.Data(), fPrefix.Data(), fPrefix2.Data(), categoryName[catIter].Data(), fCutSelection.Data(), fSuffix.Data());
            PlotDCAzInPtBinsWithBack( fESDGammaPtDCAzBins[catIter], fESDGammaPtDCAzBinsBack[catIter], NULL, nameFile, "CanvasESDDCAz", "PadESDDCAz",
                                     fDate, fMeson, 1, fNBinsPtDummy, fBinsPtDummy, "#gamma --> e^{+}e^{-}", fIsMC, "MinBias");
        }
        
        // calculate fractions per category
        for (Int_t i=0; i<4; i++) fESDGammaPerCatPtDCAzBins[i]->Divide(fDeltaPtDummy);
        for (Int_t i=0; i<3; i++) {
            fESDGammaRatioCatToCombinedPtDCAzBins[i]                        = (TH1D*)fESDGammaPerCatPtDCAzBins[i+1]->Clone(Form("ESD_GammaPtDCAzBin_Ratio_%s_to_%s", categoryName[i+1].Data(), categoryName[0].Data()));
            fESDGammaRatioCatToCombinedPtDCAzBins[i]->Sumw2();
            fESDGammaRatioCatToCombinedPtDCAzBins[i]->Divide(fESDGammaRatioCatToCombinedPtDCAzBins[i],fESDGammaPerCatPtDCAzBins[0],1,1,"B");
        }
        
        // draw fractions per category
        DrawFractionPerCat(fESDGammaRatioCatToCombinedPtDCAzBins, fOutputDir, fPrefix, fPrefix2, fCutSelection, fSuffix);
        
        // calculate ratio with/without pileup
        Double_t binContent, binError;
        for (Int_t i = 0; i < 3; i++) {
            CalculateDCAzDistributionRatio(fESDGammaPtDCAzBins, fESDSubGammaPtDCAzBins, i, 0, 0, fESDGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat[i]);
            CalculateDCAzDistributionRatio(fESDGammaPtDCAzBins, fESDSubGammaPtDCAzBins, i, 1, 3, fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[i]);
        }
        
        // define pileup correction factors
        fESDGammaPileUpCorrFactorAllCat                                     = new TH1D*[3];
        fESDGammaPileUpCorrFactor                                           = new TH1D*[3];
        for (Int_t i = 0; i < 3; i++) {
            fESDGammaPileUpCorrFactorAllCat[i]                              = new TH1D(Form("fESDGammaPileUpCorrFactorAllCatComb%i", i), Form("fESDGammaPileUpCorrFactorAllCatComb%i", i), fNBinsPt, fBinsPt);
            fESDGammaPileUpCorrFactorAllCat[i]->Sumw2();
            fESDGammaPileUpCorrFactor[i]                                    = new TH1D(Form("fESDGammaPileUpCorrFactor%i", i), Form("fESDGammaPileUpCorrFactor%i", i), fNBinsPt, fBinsPt);
            fESDGammaPileUpCorrFactor[i]->Sumw2();
        }

        // calculate pileup correction factors
        fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat           = new TF1*[3];
        fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinning                 = new TF1*[3];
        for (Int_t i = 0; i < 3; i++) {
            CalculatePileUpCorrectionFactor(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat[i], fESDGammaPileUpCorrFactorAllCat[i], fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat[i]);
            CalculatePileUpCorrectionFactor(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[i], fESDGammaPileUpCorrFactor[i], fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinning[i]);
        }

        // calculate spectra w/o pileup (using standard background extraction method)
        fESDGammaPtPileUpAllCat                                             = (TH1D*)fHistoGammaConvPt->Clone("ESD_ConvGamma_Pt_PileUp_AllCatComb");
        fESDGammaPtPileUpAllCat->Sumw2();
        fESDGammaPtPileUpAllCat->Multiply(fESDGammaPileUpCorrFactorAllCat[0]);

        fESDGammaPtPileUp                                                   = (TH1D*)fHistoGammaConvPt->Clone("ESD_ConvGamma_Pt_PileUp");
        fESDGammaPtPileUp->Sumw2();
        fESDGammaPtPileUp->Multiply(fESDGammaPileUpCorrFactor[0]);
        
        // plotting ratio + fit
        TCanvas *RatioWithWithoutPileUpCanvas                               = GetAndSetCanvas("canvasRatioWithWithoutPileUp");

        SetHistogramm(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat[0],"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.25);
        SetHistogramm(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[0],"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.25);
        SetHistogramm(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[1],"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.25);
        SetHistogramm(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[2],"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.25);

        DrawGammaSetMarker(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat[0], 25, 1.0, kGray+2, kGray+2);
        DrawGammaSetMarker(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[0], 24, 1.0, kBlack, kBlack);
        DrawGammaSetMarker(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[1], 24, 1.0, kBlue-2, kBlue-2);
        DrawGammaSetMarker(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[2], 24, 1.0, kGreen+2, kGreen+2);

        fESDGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat[0]->DrawCopy("");
        fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[0]->DrawCopy("same");
        fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[1]->DrawCopy("same");
        fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[2]->DrawCopy("same");

        if (fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat[0]) {
            fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat[0]->SetLineColor(kGray+2);
            fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat[0]->Draw("same");
        }
        if (fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinning[0]) {
            fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinning[0]->SetLineColor(kBlack);
            fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinning[0]->Draw("same");
        }
        if (fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinning[1]) {
            fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinning[1]->SetLineColor(kBlue-2);
            fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinning[1]->Draw("same");
        }
        if (fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinning[2]) {
            fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinning[2]->SetLineColor(kGreen+2);
            fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinning[2]->Draw("same");
        }

        TLegend* legendDCAZData                                             = GetAndSetLegend(0.6,0.75,4,1);
        legendDCAZData->AddEntry(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[0],"standard, sep. cat.","lp");
        legendDCAZData->AddEntry(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[1],"variation 1, sep. cat.","lp");
        legendDCAZData->AddEntry(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[2],"variation 2, sep. cat.","lp");
        legendDCAZData->AddEntry(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat[0],"standard, combined cat.","lp");
        legendDCAZData->Draw();
        
        RatioWithWithoutPileUpCanvas->Print(Form("%s/%s_%s_ESD_With_vs_Without_Pileup_pT_%s.%s",fOutputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fCutSelection.Data(),fSuffix.Data()));
        delete RatioWithWithoutPileUpCanvas;
    }
}

// *****************************************************************************************************
// *********************** Rescaling of Correction factors due to misidentified pielup *****************
// *****************************************************************************************************
void CalculatePileUpGammaCorrection(){
    
    fHistoFracAllGammaToSecPileUp                       = (TH1D*) fESDGammaPtPileUp->Clone("FracAllGammaToSecPileUp");
    fHistoFracAllGammaToSecPileUp->Divide(fTrueSecondaryConvGammaPtPileUp,fHistoFracAllGammaToSecPileUp,1,1,"B");

    fHistoFracAllGammaToSecFromXFromK0sPileUp           = (TH1D*) fESDGammaPtPileUp->Clone("FracAllGammaToSecFromXFromK0sPileUp");
    fHistoFracAllGammaToSecFromXFromK0sPileUp->Divide(fTrueSecondaryFromXFromK0sConvGammaPtPileUp,fHistoFracAllGammaToSecFromXFromK0sPileUp,1,1,"B");
    
    // ================= PURITY =================
    fHistoGammaMCPurityPileUp                           = new TH1D("GammaPurity_PileUp_Pt","",fNBinsPt,fBinsPt);
    fHistoGammaMCPurityPileUp->Sumw2();
    fHistoGammaMCPurityPileUp->Add(fTruePrimaryConvGammaPtPileUp);
    fHistoGammaMCPurityPileUp->Add(fTrueSecondaryConvGammaPtPileUp);
    fHistoGammaMCPurityPileUp->Divide(fHistoGammaMCPurityPileUp,fMCrecGammaPtPileUp,1,1,"B");

    fHistoGammaMCrecPrimaryConvPtPileUp                 = (TH1D*) fMCrecGammaPtPileUp->Clone("MCrec_PrimaryConvGamma_PtPileUp");
    fHistoGammaMCrecPrimaryConvPtPileUp->Add(fTrueSecondaryConvGammaPtPileUp,-1);
    fHistoGammaMCTruePurityPileUp                       = new TH1D("GammaTruePurity_PileUp_Pt","",fNBinsPt,fBinsPt);
    fHistoGammaMCTruePurityPileUp->Sumw2();
    fHistoGammaMCTruePurityPileUp->Divide(fTruePrimaryConvGammaPtPileUp,fHistoGammaMCrecPrimaryConvPtPileUp,1,1,"B");
    // ==========================================

    // ================ Reco Eff ================
    fHistoGammaMCRecoEffPileUp                          = new TH1D("GammaRecoEff_PileUp_Pt","",fNBinsPt,fBinsPt);
    fHistoGammaMCRecoEffPileUp->Sumw2();
    fHistoGammaMCRecoEffPileUp->Add(fTruePrimaryConvGammaPtPileUp);
    fHistoGammaMCRecoEffPileUp->Add(fTrueSecondaryConvGammaPtPileUp);
    fHistoGammaMCRecoEffPileUp->Divide(fHistoGammaMCRecoEffPileUp,fHistoGammaMCConvPt,1,1,"B");

    fHistoGammaMCPrimaryRecoEffPileUp                   = new TH1D("GammaPrimaryRecoEff_PileUp_Pt","",fNBinsPt,fBinsPt);
    fHistoGammaMCPrimaryRecoEffPileUp->Sumw2();
    fHistoGammaMCPrimaryRecoEffPileUp->Divide(fTruePrimaryConvGammaPtPileUp,fHistoGammaMCConvPt,1,1,"B");
    // ==========================================
}

// *****************************************************************************************************
// ********************** Calculation of correction factors for Gamma **********************************
// *****************************************************************************************************
void CalculateGammaCorrection(){

    if(fEnablePCM){
        TAxis *xAxis                                        = fHistoGammaTruePrimaryConv_recPt_MCPt_MC->GetXaxis();
        TAxis *yAxis                                        = fHistoGammaTruePrimaryConv_recPt_MCPt_MC->GetYaxis();

        // =========== Response matrix ==============
        fHistoGammaTruePrimaryConv_recPt_MCPt_MC_Rebin          = new TH2D("TruePrimaryConvGamma_recPt_MCPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
        for(Int_t x = 1; x<fHistoGammaTruePrimaryConv_recPt_MCPt_MC->GetNbinsX()+1; x++){
            for(Int_t y = 1; y<fHistoGammaTruePrimaryConv_recPt_MCPt_MC->GetNbinsY()+1; y++){
                Double_t binContent                         = fHistoGammaTruePrimaryConv_recPt_MCPt_MC->GetBinContent(x,y);
                Double_t xcenter                            = xAxis->GetBinCenter(x);
                Double_t ycenter                            = yAxis->GetBinCenter(y);
                fHistoGammaTruePrimaryConv_recPt_MCPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
            }
        }
        // ==========================================
        
        // ======== Secondary fractions =============
        fHistoFracAllGammaToSecOrBin                        = (TH1D*) fHistoGammaConvPtOrBin->Clone("FracAllGammaToSecOriginalBinning");
        fHistoFracAllGammaToSecOrBin->Divide(fHistoGammaTrueSecondaryConvPtOrBin,fHistoFracAllGammaToSecOrBin,1,1,"B");

        fHistoFracAllGammaToSec                             = (TH1D*) fHistoGammaConvPt->Clone("FracAllGammaToSec");
        fHistoFracAllGammaToSec->Divide(fHistoGammaTrueSecondaryConvPt,fHistoFracAllGammaToSec,1,1,"B");


        fHistoFracAllGammaToSecFromXFromK0s                 = (TH1D*) fHistoGammaConvPt->Clone("FracAllGammaToSecFromXFromK0s");
        fHistoFracAllGammaToSecFromXFromK0s->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt,fHistoFracAllGammaToSecFromXFromK0s,1,1,"B");


        fHistoFracAllGammaToSecFromXFromK0sOrBin            = (TH1D*) fHistoGammaConvPtOrBin->Clone("FracAllGammaToSecFromXFromK0sOriginalBinning");
        fHistoFracAllGammaToSecFromXFromK0sOrBin->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromK0sPtOrBin,fHistoFracAllGammaToSecFromXFromK0sOrBin,1,1,"B");


        fHistoFracAllGammaToSecFromXFromLambda              = (TH1D*) fHistoGammaConvPt->Clone("FracAllGammaToSecFromXFromLambda");
        fHistoFracAllGammaToSecFromXFromLambda->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt,fHistoFracAllGammaToSecFromXFromLambda,1,1,"B");


        fHistoFracAllGammaToSecFromXFromLambdaOrBin         = (TH1D*) fHistoGammaConvPtOrBin->Clone("FracAllGammaToSecFromXFromLambdaOriginalBinning");
        fHistoFracAllGammaToSecFromXFromLambdaOrBin->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPtOrBin,fHistoFracAllGammaToSecFromXFromLambdaOrBin,1,1,"B");
        // ==========================================
        
        // =============== Conv Prob ================
        fHistoGammaMCConvProb                               = new TH1D("MCGammaConvProb_MCPt","",fNBinsPt,fBinsPt);
        fHistoGammaMCConvProb->Sumw2();
        fHistoGammaMCConvProb->Divide(fHistoGammaMCConvPt,fHistoGammaMCAllPt,1,1,"B");
        // ==========================================

        // ================= PURITY =================
        fHistoGammaMCPurity                                 = new TH1D("GammaPurity_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCPurity->Sumw2();
        fHistoGammaMCPurity->Divide(fHistoGammaTrueConvPt,fHistoGammaConvPt,1,1,"B");


        fHistoGammaMCrecPrimaryConvPt                       = (TH1D*) fHistoGammaConvPt->Clone("MCrec_PrimaryConvGamma_Pt");
        fHistoGammaMCrecPrimaryConvPt->Add(fHistoGammaTrueSecondaryConvPt,-1);
        fHistoGammaMCTruePurity                             = new TH1D("GammaTruePurity_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCTruePurity->Sumw2();
        fHistoGammaMCTruePurity->Divide(fHistoGammaTruePrimaryConvPt,fHistoGammaMCrecPrimaryConvPt,1,1,"B");


        fHistoGammaMCrecPrimaryConvPtOrBin                  = (TH1D*) fHistoGammaConvPtOrBin->Clone("MC_ESDPrimaryConvGammaPt");
        fHistoGammaMCrecPrimaryConvPtOrBin->Add(fHistoGammaTrueSecondaryConvPtOrBin,-1);
        fHistoGammaMCTruePurityOrBin                        = (TH1D*)fHistoGammaMCrecPrimaryConvPtOrBin->Clone("GammaTruePurity_OriginalBinning_Pt");
        fHistoGammaMCTruePurityOrBin->Sumw2();
        fHistoGammaMCTruePurityOrBin->Divide(fHistoGammaTruePrimaryConvPtOrBin,fHistoGammaMCrecPrimaryConvPtOrBin,1,1,"B");
        // ==========================================

        // ================ Reco Eff ================
        fHistoGammaMCRecoEff                                = new TH1D("GammaRecoEff_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCRecoEff->Sumw2();
        fHistoGammaMCRecoEff->Divide(fHistoGammaTrueConvPt,fHistoGammaMCConvPt,1,1,"B");


        fHistoGammaMCPrimaryRecoEff                         = new TH1D("GammaPrimaryRecoEff_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCPrimaryRecoEff->Sumw2();
        fHistoGammaMCPrimaryRecoEff->Divide(fHistoGammaTruePrimaryConvPt,fHistoGammaMCConvPt,1,1,"B");


        fHistoGammaMCPrimaryRecoEffMCPt                     = new TH1D("GammaPrimaryRecoEff_MCPt","",fNBinsPt,fBinsPt);
        fHistoGammaMCPrimaryRecoEffMCPt->Sumw2();
        fHistoGammaMCPrimaryRecoEffMCPt->Divide(fHistoGammaTruePrimaryConvMCPt,fHistoGammaMCConvPt,1,1,"B");
        // ==========================================

        // ========== identified MC BG ==============
        fHistoGammaMCBackground                             = new TH1D("MCrec_Background","",fNBinsPt,fBinsPt);
        fHistoGammaMCBackground->Sumw2();
        fHistoGammaMCBackground = (TH1D*)fHistoGammaMCrecConvPt->Clone("MCrec_Background");
        fHistoGammaMCBackground->Add(fHistoGammaTrueConvPt,-1);
        // ==========================================
    } 
    
    if (fEnableCalo && fEnablePCM){
        TAxis *xAxis                                        = fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetXaxis();
        TAxis *yAxis                                        = fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetYaxis();

        // =========== Response matrix ==============
        fHistoGammaTruePrimaryCalo_recPt_MCPt_MC_Rebin      = new TH2D("TruePrimaryCaloGamma_recPt_MCPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
        for(Int_t x = 1; x<fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetNbinsX()+1; x++){
            for(Int_t y = 1; y<fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetNbinsY()+1; y++){
                Double_t binContent                         = fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetBinContent(x,y);
                Double_t xcenter                            = xAxis->GetBinCenter(x);
                Double_t ycenter                            = yAxis->GetBinCenter(y);
                fHistoGammaTruePrimaryCalo_recPt_MCPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
            }
        }
        // ==========================================
        
        // ======== Secondary fractions =============
        fHistoFracAllGammaCaloToSecOrBin                    = (TH1D*) fHistoGammaCaloPtOrBin->Clone("FracAllGammaCaloToSecOriginalBinning");
        fHistoFracAllGammaCaloToSecOrBin->Divide(fHistoGammaTrueSecondaryCaloPtOrBin,fHistoFracAllGammaCaloToSecOrBin,1,1,"B");


        fHistoFracAllGammaCaloToSec                         = (TH1D*) fHistoGammaCaloPt->Clone("FracAllGammaCaloToSec");
        fHistoFracAllGammaCaloToSec->Divide(fHistoGammaTrueSecondaryCaloPt,fHistoFracAllGammaCaloToSec,1,1,"B");


        fHistoFracAllGammaCaloToSecFromK0sOrBin             = (TH1D*) fHistoGammaCaloPtOrBin->Clone("FracAllGammaCaloToSecFromK0sOriginalBinning");
        fHistoFracAllGammaCaloToSecFromK0sOrBin->Divide(fHistoGammaTrueSecondaryCaloFromK0sPtOrBin,fHistoFracAllGammaCaloToSecFromK0sOrBin,1,1,"B");


        fHistoFracAllGammaCaloToSecFromK0s                  = (TH1D*) fHistoGammaCaloPt->Clone("fHistoFracAllGammaCaloToSecFromK0s");
        fHistoFracAllGammaCaloToSecFromK0s->Divide(fHistoGammaTrueSecondaryCaloFromK0sPt,fHistoFracAllGammaCaloToSecFromK0s,1,1,"B");


        fHistoFracAllGammaCaloToSecFromLambdaOrBin          = (TH1D*) fHistoGammaCaloPtOrBin->Clone("FracAllGammaCaloToSecFromLambdaOriginalBinning");
        fHistoFracAllGammaCaloToSecFromLambdaOrBin->Divide(fHistoGammaTrueSecondaryCaloFromLambdaPtOrBin,fHistoFracAllGammaCaloToSecFromLambdaOrBin,1,1,"B");


        fHistoFracAllGammaCaloToSecFromLambda               = (TH1D*) fHistoGammaCaloPt->Clone("fHistoFracAllGammaCaloToSecFromLambda");
        fHistoFracAllGammaCaloToSecFromLambda->Divide(fHistoGammaTrueSecondaryCaloFromLambdaPt,fHistoFracAllGammaCaloToSecFromLambda,1,1,"B");

        
         // ================= PURITY =================
        fHistoGammaCaloMCPurity                             = new TH1D("GammaCaloPurity_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaCaloMCPurity->Sumw2();
        fHistoGammaCaloMCPurity->Divide(fHistoGammaTrueCaloPt,fHistoGammaCaloPt,1,1,"B");


        fHistoGammaMCrecPrimaryCaloPt                       = (TH1D*) fHistoGammaCaloPt->Clone("MCrec_PrimaryCaloGamma_Pt");
        fHistoGammaMCrecPrimaryCaloPt->Add(fHistoGammaTrueSecondaryCaloPt,-1);
        fHistoGammaCaloMCTruePurity                         = new TH1D("GammaCaloTruePurity_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaCaloMCTruePurity->Sumw2();
        fHistoGammaCaloMCTruePurity->Divide(fHistoGammaTruePrimaryCaloPt,fHistoGammaMCrecPrimaryCaloPt,1,1,"B");


        fHistoGammaMCrecPrimaryCaloPtOrBin                  = (TH1D*) fHistoGammaCaloPtOrBin->Clone("MC_ESDPrimaryCaloGammaPt");
        fHistoGammaMCrecPrimaryCaloPtOrBin->Add(fHistoGammaTrueSecondaryCaloPtOrBin,-1);
        fHistoGammaCaloMCTruePurityOrBin                    = (TH1D*)fHistoGammaMCrecPrimaryCaloPtOrBin->Clone("GammaCaloTruePurity_OriginalBinning_Pt");
        fHistoGammaCaloMCTruePurityOrBin->Sumw2();
        fHistoGammaCaloMCTruePurityOrBin->Divide(fHistoGammaTruePrimaryCaloPtOrBin,fHistoGammaMCrecPrimaryCaloPtOrBin,1,1,"B");
        // ==========================================
       
        // ================ Reco Eff ================
        fHistoGammaCaloMCRecoEff                            = new TH1D("GammaCaloRecoEff_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaCaloMCRecoEff->Sumw2();
        fHistoGammaCaloMCRecoEff->Divide(fHistoGammaTrueCaloPt,fHistoGammaMCAllInEMCAccPt,1,1,"B");


        fHistoGammaCaloMCPrimaryRecoEff                     = new TH1D("GammaCaloPrimaryRecoEff_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaCaloMCPrimaryRecoEff->Sumw2();
        fHistoGammaCaloMCPrimaryRecoEff->Divide(fHistoGammaTruePrimaryCaloPt,fHistoGammaMCAllInEMCAccPt,1,1,"B");


        fHistoGammaCaloMCPrimaryRecoEffMCPt                 = new TH1D("GammaCaloPrimaryRecoEff_MCPt","",fNBinsPt,fBinsPt);
        fHistoGammaCaloMCPrimaryRecoEffMCPt->Sumw2();
        fHistoGammaCaloMCPrimaryRecoEffMCPt->Divide(fHistoGammaTruePrimaryCaloMCPt,fHistoGammaMCAllInEMCAccPt,1,1,"B");
        // ==========================================

        // ========== identified MC BG ==============
        fHistoGammaCaloMCBackground                         = new TH1D("MCrec_Calo_Background","",fNBinsPt,fBinsPt);
        fHistoGammaCaloMCBackground->Sumw2();
        fHistoGammaCaloMCBackground = (TH1D*)fHistoGammaMCrecCaloPt->Clone("MCrec_Calo_Background");
        fHistoGammaCaloMCBackground->Add(fHistoGammaTrueCaloPt,-1);
        // ==========================================

        
    }    
    
    if(fEnableCalo && !fEnablePCM){
        TAxis *xAxis                                        = fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetXaxis();
        TAxis *yAxis                                        = fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetYaxis();

        // =========== Response matrix ==============
        fHistoGammaTruePrimaryCalo_recPt_MCPt_MC_Rebin      = new TH2D("TruePrimaryCaloGamma_recPt_MCPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
        for(Int_t x = 1; x<fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetNbinsX()+1; x++){
            for(Int_t y = 1; y<fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetNbinsY()+1; y++){
                Double_t binContent                         = fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetBinContent(x,y);
                Double_t xcenter                            = xAxis->GetBinCenter(x);
                Double_t ycenter                            = yAxis->GetBinCenter(y);
                fHistoGammaTruePrimaryCalo_recPt_MCPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
            }
        }
        
        xAxis                                               = fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->GetXaxis();
        yAxis                                               = fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->GetYaxis();
        // =========== Response matrix ==============
        fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC_Rebin  = new TH2D("TruePrimaryCaloConvGamma_recPt_MCPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
        for(Int_t x = 1; x<fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->GetNbinsX()+1; x++){
            for(Int_t y = 1; y<fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->GetNbinsY()+1; y++){
                Double_t binContent                         = fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->GetBinContent(x,y);
                Double_t xcenter                            = xAxis->GetBinCenter(x);
                Double_t ycenter                            = yAxis->GetBinCenter(y);
                fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
            }
        }

        xAxis                                               = fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->GetXaxis();
        yAxis                                               = fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->GetYaxis();
        // =========== Response matrix ==============
        fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC_Rebin  = new TH2D("TruePrimaryCaloUnConvGamma_recPt_MCPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
        for(Int_t x = 1; x<fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->GetNbinsX()+1; x++){
            for(Int_t y = 1; y<fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->GetNbinsY()+1; y++){
                Double_t binContent                         = fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->GetBinContent(x,y);
                Double_t xcenter                            = xAxis->GetBinCenter(x);
                Double_t ycenter                            = yAxis->GetBinCenter(y);
                fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
            }
        }
        // ==========================================
        
        // ======== Secondary fractions =============
        fHistoFracAllGammaToSecOrBin                        = (TH1D*) fHistoGammaCaloPtOrBin->Clone("FracAllGammaToSecOriginalBinning");
        fHistoFracAllGammaToSecOrBin->Divide(fHistoGammaTrueSecondaryCaloPtOrBin,fHistoFracAllGammaToSecOrBin,1,1,"B");


        fHistoFracAllGammaToSec                             = (TH1D*) fHistoGammaCaloPt->Clone("FracAllGammaToSec");
        fHistoFracAllGammaToSec->Divide(fHistoGammaTrueSecondaryCaloPt,fHistoFracAllGammaToSec,1,1,"B");


        fHistoFracAllGammaToSecFromXFromK0s                 = (TH1D*) fHistoGammaCaloPt->Clone("FracAllGammaToSecFromXFromK0s");
        fHistoFracAllGammaToSecFromXFromK0s->Divide(fHistoGammaTrueSecondaryCaloFromK0sPt,fHistoFracAllGammaToSecFromXFromK0s,1,1,"B");


        fHistoFracAllGammaToSecFromXFromK0sOrBin            = (TH1D*) fHistoGammaCaloPtOrBin->Clone("FracAllGammaToSecFromXFromK0sOriginalBinning");
        fHistoFracAllGammaToSecFromXFromK0sOrBin->Divide(fHistoGammaTrueSecondaryCaloFromK0sPtOrBin,fHistoFracAllGammaToSecFromXFromK0sOrBin,1,1,"B");


        fHistoFracAllGammaToSecFromXFromLambda              = (TH1D*) fHistoGammaCaloPt->Clone("FracAllGammaToSecFromXFromLambda");
        fHistoFracAllGammaToSecFromXFromLambda->Divide(fHistoGammaTrueSecondaryCaloFromLambdaPt,fHistoFracAllGammaToSecFromXFromLambda,1,1,"B");


        fHistoFracAllGammaToSecFromXFromLambdaOrBin         = (TH1D*) fHistoGammaCaloPtOrBin->Clone("FracAllGammaToSecFromXFromLambdaOriginalBinning");
        fHistoFracAllGammaToSecFromXFromLambdaOrBin->Divide(fHistoGammaTrueSecondaryCaloFromLambdaPtOrBin,fHistoFracAllGammaToSecFromXFromLambdaOrBin,1,1,"B");
        // ==========================================
        
        // ================= PURITY =================
        fHistoGammaMCPurity                                 = new TH1D("GammaPurity_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCPurity->Sumw2();
        fHistoGammaMCPurity->Divide(fHistoGammaTrueCaloPt,fHistoGammaCaloPt,1,1,"B");


        fHistoGammaMCrecPrimaryCaloPt                       = (TH1D*) fHistoGammaCaloPt->Clone("MCrec_PrimaryCaloGamma_Pt");
        fHistoGammaMCrecPrimaryCaloPt->Add(fHistoGammaTrueSecondaryCaloPt,-1);
        fHistoGammaMCTruePurity                             = new TH1D("GammaTruePurity_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCTruePurity->Sumw2();
        fHistoGammaMCTruePurity->Divide(fHistoGammaTruePrimaryCaloPt,fHistoGammaMCrecPrimaryCaloPt,1,1,"B");


        fHistoGammaMCrecPrimaryCaloPtOrBin                  = (TH1D*) fHistoGammaCaloPtOrBin->Clone("MC_ESDPrimaryCaloGammaPt");
        fHistoGammaMCrecPrimaryCaloPtOrBin->Add(fHistoGammaTrueSecondaryCaloPtOrBin,-1);
        fHistoGammaMCTruePurityOrBin                        = (TH1D*)fHistoGammaMCrecPrimaryCaloPtOrBin->Clone("GammaTruePurity_OriginalBinning_Pt");
        fHistoGammaMCTruePurityOrBin->Sumw2();
        fHistoGammaMCTruePurityOrBin->Divide(fHistoGammaTruePrimaryCaloPtOrBin,fHistoGammaMCrecPrimaryCaloPtOrBin,1,1,"B");
        // ==========================================

        // ================ Reco Eff ================
        fHistoGammaMCRecoEff                                = new TH1D("GammaRecoEff_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCRecoEff->Sumw2();
        fHistoGammaMCRecoEff->Divide(fHistoGammaTrueCaloPt,fHistoGammaMCAllPt,1,1,"B");


        fHistoGammaMCPrimaryRecoEff                         = new TH1D("GammaPrimaryRecoEff_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCPrimaryRecoEff->Sumw2();
        fHistoGammaMCPrimaryRecoEff->Divide(fHistoGammaTruePrimaryCaloPt,fHistoGammaMCAllPt,1,1,"B");


        fHistoGammaMCPrimaryRecoEffMCPt                     = new TH1D("GammaPrimaryRecoEff_MCPt","",fNBinsPt,fBinsPt);
        fHistoGammaMCPrimaryRecoEffMCPt->Sumw2();
        fHistoGammaMCPrimaryRecoEffMCPt->Divide(fHistoGammaTruePrimaryCaloMCPt,fHistoGammaMCAllPt,1,1,"B");
        // ==========================================

        // ========== identified MC BG ==============
        fHistoGammaMCBackground                             = new TH1D("MCrec_Background","",fNBinsPt,fBinsPt);
        fHistoGammaMCBackground->Sumw2();
        fHistoGammaMCBackground = (TH1D*)fHistoGammaMCrecCaloPt->Clone("MCrec_Background");
        fHistoGammaMCBackground->Add(fHistoGammaTrueCaloPt,-1);
        // ==========================================
    } 
}

// *****************************************************************************************************
// *********************** Initialize histograms and binning *******************************************
// *****************************************************************************************************
void Initialize(TString setPi0, TString energy , Int_t numberOfBins, Int_t mode, Bool_t addSig){

    InitializeBinning(setPi0, numberOfBins, energy, fDirectPhoton, fMode, fEventCutNumber, fClusterCutNumber);

    fDeltaPt                                            = new TH1D("deltaPt","",fNBinsPt,fBinsPt);
    for(Int_t iPt=fStartPtBin+1;iPt<fNBinsPt+1;iPt++){
        fDeltaPt->SetBinContent(iPt,fBinsPt[iPt]-fBinsPt[iPt-1]);
        fDeltaPt->SetBinError(iPt,0);
    }

    // initializing binning used for the DCAz distributions
    if(!addSig && !(mode == 4 || mode == 5)){
        if (fBinsPtDCAzDist && fNBinsPtDCAzDist ) {
            if (fBinsPtDCAzDist[0] == fBinsPt[0] && fBinsPtDCAzDist[fNBinsPtDCAzDist] == fBinsPt[fNBinsPt]) {
                cout << "A different binning will be used for the DCAz distributions." << endl;

                fNBinsPtDummy                           = fNBinsPtDCAzDist;
                fBinsPtDummy                            = fBinsPtDCAzDist;
            } else {
                cout << "WARNING: The bin range chosen for the DCAz distributions doesn't coincide with the one used for the spectra, using the usual binning." << endl;
                
                fNBinsPtDummy                           = fNBinsPt;
                fBinsPtDummy                            = fBinsPt;
            }
        } else {
            cout << "There is no binning for the DCAz distributions defined, using the usual binning." << endl;
            fNBinsPtDummy                               = fNBinsPt;
            fBinsPtDummy                                = fBinsPt;
        }
    
        fDeltaPtDummy                                   = new TH1D("deltaPtDummy","",fNBinsPtDummy,fBinsPtDummy);
        for(Int_t iPt=fStartPtBin+1;iPt<fNBinsPtDummy+1;iPt++){
            fDeltaPtDummy->SetBinContent(iPt,fBinsPtDummy[iPt]-fBinsPtDummy[iPt-1]);
            fDeltaPtDummy->SetBinError(iPt,0);
        }
    }
    
    // initialize ShowBackground for DCAz distributions
    if ((fEnergyFlag.CompareTo("13TeV") == 0) && (fMeson.CompareTo("Pi0") == 0) && (fDirectPhoton.CompareTo("directPhoton") == 0)) {
        nIterationsShowBackground[0]                    = 13;
        nIterationsShowBackground[1]                    = 13;
        nIterationsShowBackground[2]                    = 18;
        nIterationsShowBackground[3]                    = 20;
        optionShowBackground[0]                         = "BackDecreasingWindow";                   // standard
        optionShowBackground[1]                         = "nosmoothing";
        optionShowBackground[2]                         = "BackDecreasingWindow, BackSmoothing5";
    } else {
        cout << "WARNING: No ShowBackground-options defined, using the default ones." << endl;
        nIterationsShowBackground[0]                    = 12;
        nIterationsShowBackground[1]                    = 12;
        nIterationsShowBackground[2]                    = 18;
        nIterationsShowBackground[3]                    = 20;
        optionShowBackground[0]                         = "BackDecreasingWindow, BackSmoothing5";   // standard
        optionShowBackground[1]                         = "BackDecreasingWindow, BackSmoothing3";
        optionShowBackground[2]                         = "BackDecreasingWindow, BackSmoothing7";
    }
    
    // initialize mass histo array
    fHistoGconvGInvMassPtGConvBin                       = new TH1D*[fNBinsPt];
    for(Int_t i = 0;i<fNBinsPt; i++){
        fHistoGconvGInvMassPtGConvBin[i]                = NULL;
    }
}

// *****************************************************************************************************
// ************************ Create histos from DCA-tree for pileup estimate ****************************
// *****************************************************************************************************
void FillDCAHistogramsFromTree(TTree *dcaTree,Bool_t isMC){

    Float_t dcaZPhoton, pt;
    UChar_t cat, photonMCInfo;
    dcaTree->SetBranchAddress("Pt",&pt);
    dcaTree->SetBranchAddress("DcaZPhoton",&dcaZPhoton);
    dcaTree->SetBranchAddress("cat",&cat);
    if (isMC) dcaTree->SetBranchAddress("photonMCInfo",&photonMCInfo);

    if(!isMC){
        fESDGammaPtDCAz     = new TH2F*[6];
        fESDGammaPtDCAz[0]  = new TH2F("ESD_GammaPtDCAz_cat0","ESD_GammaPtDCAz_cat0",250,0,25,201,-10,10);
        fESDGammaPtDCAz[0]->Sumw2();
        fESDGammaPtDCAz[1]  = new TH2F("ESD_GammaPtDCAz_cat1","ESD_GammaPtDCAz_cat1",250,0,25,201,-10,10);
        fESDGammaPtDCAz[1]->Sumw2();
        fESDGammaPtDCAz[2]  = new TH2F("ESD_GammaPtDCAz_cat2","ESD_GammaPtDCAz_cat2",250,0,25,201,-10,10);
        fESDGammaPtDCAz[2]->Sumw2();
        fESDGammaPtDCAz[3]  = new TH2F("ESD_GammaPtDCAz_cat3","ESD_GammaPtDCAz_cat3",250,0,25,201,-10,10);
        fESDGammaPtDCAz[3]->Sumw2();
        fESDGammaPtDCAz[4]  = new TH2F("ESD_GammaPtDCAz_cat23","ESD_GammaPtDCAz_cat23",250,0,25,201,-10,10);
        fESDGammaPtDCAz[4]->Sumw2();
        fESDGammaPtDCAz[5]  = new TH2F("ESD_GammaPtDCAz_all","ESD_GammaPtDCAz_all",250,0,25,201,-10,10);
        fESDGammaPtDCAz[5]->Sumw2();

        Long64_t nentries = dcaTree->GetEntries();
        for (Long64_t l=0;l<nentries;l++) {
            dcaTree->GetEntry(l);
            if(cat == 0) fESDGammaPtDCAz[0]->Fill(pt,dcaZPhoton);
            if(cat == 1) fESDGammaPtDCAz[1]->Fill(pt,dcaZPhoton);
            if(cat == 2) fESDGammaPtDCAz[2]->Fill(pt,dcaZPhoton);
            if(cat == 3) fESDGammaPtDCAz[3]->Fill(pt,dcaZPhoton);
            if(cat == 2 || cat == 3) fESDGammaPtDCAz[4]->Fill(pt,dcaZPhoton);
            if(cat != 0)fESDGammaPtDCAz[5]->Fill(pt,dcaZPhoton);
        }
    }
    if(isMC){
        fMCrecGammaPtDCAz                           = new TH2F("MCrec_GammaPtDCAz_all","MCrec_GammaPtDCAz_all",250,0,25,201,-10,10);
        fMCrecGammaPtDCAz->Sumw2();
        fTruePrimaryPhotonPtDCAz                    = new TH2F*[6];
        fTrueSecondaryPhotonPtDCAz                  = new TH2F*[6];
        fTrueSecondaryPhotonFromXFromK0sPtDCAz      = new TH2F*[6];
        fTruePrimaryPhotonPtDCAz[0]                 = new TH2F("ESD_TruePrimaryGammaPtDCAz_cat0","ESD_TruePrimaryGammaPtDCAz_cat0",250,0,25,201,-10,10);
        fTrueSecondaryPhotonPtDCAz[0]               = new TH2F("ESD_TrueSecondaryGammaPtDCAz_cat0","ESD_TrueSecondaryGammaPtDCAz_cat0",250,0,25,201,-10,10);
        fTrueSecondaryPhotonFromXFromK0sPtDCAz[0]   = new TH2F("ESD_TrueSecondaryGammaFromXFromK0sDCAz_cat0","ESD_TrueSecondaryGammaFromXFromK0sDCAz_cat0",250,0,25,201,-10,10);
        fTruePrimaryPhotonPtDCAz[0]->Sumw2();
        fTrueSecondaryPhotonPtDCAz[0]->Sumw2();
        fTrueSecondaryPhotonFromXFromK0sPtDCAz[0]->Sumw2();
        fTruePrimaryPhotonPtDCAz[1]                 = new TH2F("ESD_TruePrimaryGammaPtDCAz_cat1","ESD_TruePrimaryGammaPtDCAz_cat1",250,0,25,201,-10,10);
        fTrueSecondaryPhotonPtDCAz[1]               = new TH2F("ESD_TrueSecondaryGammaPtDCAz_cat1","ESD_TrueSecondaryGammaPtDCAz_cat1",250,0,25,201,-10,10);
        fTrueSecondaryPhotonFromXFromK0sPtDCAz[1]   = new TH2F("ESD_TrueSecondaryGammaFromXFromK0sDCAz_cat1","ESD_TrueSecondaryGammaFromXFromK0sDCAz_cat1",250,0,25,201,-10,10);
        fTruePrimaryPhotonPtDCAz[1]->Sumw2();
        fTrueSecondaryPhotonPtDCAz[1]->Sumw2();
        fTrueSecondaryPhotonFromXFromK0sPtDCAz[1]->Sumw2();
        fTruePrimaryPhotonPtDCAz[2]                 = new TH2F("ESD_TruePrimaryGammaPtDCAz_cat2","ESD_TruePrimaryGammaPtDCAz_cat2",250,0,25,201,-10,10);
        fTrueSecondaryPhotonPtDCAz[2]               = new TH2F("ESD_TrueSecondaryGammaPtDCAz_cat2","ESD_TrueSecondaryGammaPtDCAz_cat2",250,0,25,201,-10,10);
        fTrueSecondaryPhotonFromXFromK0sPtDCAz[2]   = new TH2F("ESD_TrueSecondaryGammaFromXFromK0sDCAz_cat2","ESD_TrueSecondaryGammaFromXFromK0sDCAz_cat2",250,0,25,201,-10,10);
        fTruePrimaryPhotonPtDCAz[2]->Sumw2();
        fTrueSecondaryPhotonPtDCAz[2]->Sumw2();
        fTrueSecondaryPhotonFromXFromK0sPtDCAz[2]->Sumw2();
        fTruePrimaryPhotonPtDCAz[3]                 = new TH2F("ESD_TruePrimaryGammaPtDCAz_cat3","ESD_TruePrimaryGammaPtDCAz_cat3",250,0,25,201,-10,10);
        fTrueSecondaryPhotonPtDCAz[3]               = new TH2F("ESD_TrueSecondaryGammaPtDCAz_cat3","ESD_TrueSecondaryGammaPtDCAz_cat3",250,0,25,201,-10,10);
        fTrueSecondaryPhotonFromXFromK0sPtDCAz[3]   = new TH2F("ESD_TrueSecondaryGammaFromXFromK0sDCAz_cat3","ESD_TrueSecondaryGammaFromXFromK0sDCAz_cat3",250,0,25,201,-10,10);
        fTruePrimaryPhotonPtDCAz[3]->Sumw2();
        fTrueSecondaryPhotonPtDCAz[3]->Sumw2();
        fTrueSecondaryPhotonFromXFromK0sPtDCAz[3]->Sumw2();
        fTruePrimaryPhotonPtDCAz[4]                 = new TH2F("ESD_TruePrimaryGammaPtDCAz_cat23","ESD_TruePrimaryGammaPtDCAz_cat23",250,0,25,201,-10,10);
        fTrueSecondaryPhotonPtDCAz[4]               = new TH2F("ESD_TrueSecondaryGammaPtDCAz_cat23","ESD_TrueSecondaryGammaPtDCAz_cat23",250,0,25,201,-10,10);
        fTrueSecondaryPhotonFromXFromK0sPtDCAz[4]   = new TH2F("ESD_TrueSecondaryGammaFromXFromK0sDCAz_cat23","ESD_TrueSecondaryGammaFromXFromK0sDCAz_cat23",250,0,25,101,-10,10);
        fTruePrimaryPhotonPtDCAz[4]->Sumw2();
        fTrueSecondaryPhotonPtDCAz[4]->Sumw2();
        fTrueSecondaryPhotonFromXFromK0sPtDCAz[4]->Sumw2();
        fTruePrimaryPhotonPtDCAz[5]                 = new TH2F("ESD_TruePrimaryGammaPtDCAz_all","ESD_TruePrimaryGammaPtDCAz_all",250,0,25,201,-10,10);
        fTrueSecondaryPhotonPtDCAz[5]               = new TH2F("ESD_TrueSecondaryGammaPtDCAz_all","ESD_TrueSecondaryGammaPtDCAz_all",250,0,25,201,-10,10);
        fTrueSecondaryPhotonFromXFromK0sPtDCAz[5]   = new TH2F("ESD_TrueSecondaryGammaFromXFromK0sDCAz_all","ESD_TrueSecondaryGammaFromXFromK0sDCAz_all",250,0,25,201,-10,10);
        fTruePrimaryPhotonPtDCAz[5]->Sumw2();
        fTrueSecondaryPhotonPtDCAz[5]->Sumw2();
        fTrueSecondaryPhotonFromXFromK0sPtDCAz[5]->Sumw2();

        Long64_t nentries = dcaTree->GetEntries();
        for (Long64_t l=0;l<nentries;l++) {
            dcaTree->GetEntry(l);
            fMCrecGammaPtDCAz->Fill(pt,dcaZPhoton);
            if(cat == 0){
                if(photonMCInfo == 6) fTruePrimaryPhotonPtDCAz[0]->Fill(pt,dcaZPhoton);
                if(photonMCInfo == 2 || photonMCInfo == 3 || photonMCInfo == 4 || photonMCInfo == 5)
                fTrueSecondaryPhotonPtDCAz[0]->Fill(pt,dcaZPhoton);
                if(photonMCInfo == 4) fTrueSecondaryPhotonFromXFromK0sPtDCAz[0]->Fill(pt,dcaZPhoton);
            }
            if(cat == 1){
                if(photonMCInfo == 6) fTruePrimaryPhotonPtDCAz[1]->Fill(pt,dcaZPhoton);
                if(photonMCInfo == 2 || photonMCInfo == 3 || photonMCInfo == 4 || photonMCInfo == 5)
                fTrueSecondaryPhotonPtDCAz[1]->Fill(pt,dcaZPhoton);
                if(photonMCInfo == 4) fTrueSecondaryPhotonFromXFromK0sPtDCAz[1]->Fill(pt,dcaZPhoton);
            }
            if(cat == 2){
                if(photonMCInfo == 6) fTruePrimaryPhotonPtDCAz[2]->Fill(pt,dcaZPhoton);
                if(photonMCInfo == 2 || photonMCInfo == 3 || photonMCInfo == 4 || photonMCInfo == 5)
                fTrueSecondaryPhotonPtDCAz[2]->Fill(pt,dcaZPhoton);
                if(photonMCInfo == 4) fTrueSecondaryPhotonFromXFromK0sPtDCAz[2]->Fill(pt,dcaZPhoton);
            }
            if(cat == 3){
                if(photonMCInfo == 6) fTruePrimaryPhotonPtDCAz[3]->Fill(pt,dcaZPhoton);
                if(photonMCInfo == 2 || photonMCInfo == 3 || photonMCInfo == 4 || photonMCInfo == 5)
                fTrueSecondaryPhotonPtDCAz[3]->Fill(pt,dcaZPhoton);
                if(photonMCInfo == 4) fTrueSecondaryPhotonFromXFromK0sPtDCAz[3]->Fill(pt,dcaZPhoton);
            }
            if(cat == 2 || cat == 3){
                if(photonMCInfo == 6) fTruePrimaryPhotonPtDCAz[4]->Fill(pt,dcaZPhoton);
                if(photonMCInfo == 2 || photonMCInfo == 3 || photonMCInfo == 4 || photonMCInfo == 5)
                fTrueSecondaryPhotonPtDCAz[4]->Fill(pt,dcaZPhoton);
                if(photonMCInfo == 4) fTrueSecondaryPhotonFromXFromK0sPtDCAz[4]->Fill(pt,dcaZPhoton);
            }
            if(cat != 0){
                if(photonMCInfo == 6) fTruePrimaryPhotonPtDCAz[5]->Fill(pt,dcaZPhoton);
                if(photonMCInfo == 2 || photonMCInfo == 3 || photonMCInfo == 4 || photonMCInfo == 5)
                fTrueSecondaryPhotonPtDCAz[5]->Fill(pt,dcaZPhoton);
                if(photonMCInfo == 4) fTrueSecondaryPhotonFromXFromK0sPtDCAz[5]->Fill(pt,dcaZPhoton);
            }
        }
    }
}

//**************************************************************************************************
//****************** Routine for saving histograms concerning photons mainly for data **************
//**************************************************************************************************
void SaveHistos(Int_t isMC, TString fCutID, TString fPrefix3,Bool_t PileUpCorrection){
    const char* nameOutput  = Form("%s/%s/%s_%s_GammaConvV1WithoutCorrection%s_%s.root", fCutSelection.Data(), fEnergyFlag.Data(), fPrefix.Data(),
                                   fPrefix3.Data(), fPeriodFlag.Data(), fCutID.Data());
    cout << "INFO: writing into: " << nameOutput << endl;
    TFile *Output1          = new TFile(nameOutput,"UPDATE");
        
        // write event histogram for data
        fEventQuality->Write("NEvents",TObject::kOverwrite);
        // write binning histogram
        fDeltaPt->Write("deltaPt",TObject::kOverwrite);
        if (fDeltaPtDummy) fDeltaPtDummy->Write("deltaPtDCAzDistBinning",TObject::kOverwrite);

        // write histo with reasons for photon rejection
        if (fHistoPhotonIsSelected) fHistoPhotonIsSelected->Write(Form("IsPhotonSelected %s",fGammaCutSelectionRead.Data()),TObject::kOverwrite);
        
        // write raw photons distribution
        if (fHistoGammaConvPt) fHistoGammaConvPt->Write("ESD_ConvGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaConvPtOrBin) fHistoGammaConvPtOrBin->Write("ESD_ConvGamma_Pt_OriginalBinning",TObject::kOverwrite);
        
        // write raw photon distributions after pileup subtraction
        if(PileUpCorrection && fEnablePCM){
            fESDGammaPtPileUpAllCat->Write("ESD_ConvGamma_Pt_PileUp_AllCatComb",TObject::kOverwrite);
            fESDGammaPtPileUp->Write("ESD_ConvGamma_Pt_PileUp",TObject::kOverwrite);
            for (Int_t i = 0; i < 4; i++) {
                fESDGammaPtDCAzBins[i][0]->Write(Form("ESD_GammaPtDCAzBin_Full_%s", categoryName[i].Data()),TObject::kOverwrite);
                fESDGammaPerCatPtDCAzBins[i]->Write(Form("ESD_GammaPtDCAzBin_%s", categoryName[i].Data()),TObject::kOverwrite);
                if (i < 3) fESDGammaRatioCatToCombinedPtDCAzBins[i]->Write(Form("ESD_GammaPtDCAzBin_Ratio_%s_to_%s", categoryName[i+1].Data(), categoryName[0].Data()),TObject::kOverwrite);
            }
        }

        // write basics MC quantities if possible (converted input spectrum & reconstructed validated (&primary) photons)
        if(isMC){
            if (fHistoGammaMCConvPt) fHistoGammaMCConvPt->Write("MC_ConvGamma_MCPt",TObject::kOverwrite);
            if (fHistoGammaTrueConvPt) fHistoGammaTrueConvPt->Write("TrueConvGamma_Pt",TObject::kOverwrite);
            if (fHistoGammaTruePrimaryConvPt) fHistoGammaTruePrimaryConvPt->Write("TruePrimaryConvGamma_Pt",TObject::kOverwrite);
            if (fHistoGammaTrueCaloPt) fHistoGammaTrueCaloPt->Write("TrueCaloGamma_Pt",TObject::kOverwrite);
            if (fHistoGammaTruePrimaryCaloPt) fHistoGammaTruePrimaryCaloPt->Write("TruePrimaryCaloGamma_Pt",TObject::kOverwrite);
        }
        
        // write raw photons distribution
        if (fHistoGammaCaloPt) fHistoGammaCaloPt->Write("ESD_CaloGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaCaloPtOrBin) fHistoGammaCaloPtOrBin->Write("ESD_CaloGamma_Pt_OriginalBinning",TObject::kOverwrite);
    
    Output1->Write();
    Output1->Close();
}


//**************************************************************************************************
//************* Routine for saving correction histograms concerning photons, only run for MC *******
//**************************************************************************************************
void SaveCorrectionHistos(TString fCutID, TString fPrefix3,Bool_t PileUpCorrection){
    
    const char* nameOutput  = Form("%s/%s/%s_%s_GammaConvV1CorrectionHistos%s_%s.root",fCutSelection.Data(), fEnergyFlag.Data(), fPrefix.Data(), fPrefix3.Data(), fPeriodFlag.Data(), fCutID.Data());
    cout << "INFO: writing into: " << nameOutput << endl;
    TFile *Output2          = new TFile(nameOutput,"UPDATE");
    
        // write input gamma distributions
        if (fHistoGammaMCAllPt) fHistoGammaMCAllPt->Write("MC_AllGamma_MCPt",TObject::kOverwrite);
        if (fHistoGammaMCAllPtOrBin) fHistoGammaMCAllPtOrBin->Write("MC_AllGamma_OriginalBinning_MCPt",TObject::kOverwrite);
        if (fHistoGammaMCAllInEMCAccPt) fHistoGammaMCAllInEMCAccPt->Write("MC_AllGammaEMCAcc_MCPt",TObject::kOverwrite);
        if (fHistoGammaMCAllInEMCAccPtOrBin) fHistoGammaMCAllInEMCAccPtOrBin->Write("MC_AllGammaEMCAcc_OriginalBinning_MCPt",TObject::kOverwrite);

        // write input converted photon distributions
        if (fHistoGammaMCConvPt) fHistoGammaMCConvPt->Write("MC_ConvGamma_MCPt",TObject::kOverwrite);

        // write photon candidates in MC
        if (fHistoGammaMCrecConvPt) fHistoGammaMCrecConvPt->Write("MCrec_ConvGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaMCrecConvPtOrBin) fHistoGammaMCrecConvPtOrBin->Write("MCrec_ConvGamma_OriginalBinning_Pt",TObject::kOverwrite);
        if (fHistoGammaMCrecCaloPt) fHistoGammaMCrecCaloPt->Write("MCrec_CaloGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaMCrecPrimaryConvPt) fHistoGammaMCrecPrimaryConvPt->Write("MCrec_PrimaryConvGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaMCrecPrimaryCaloPt) fHistoGammaMCrecPrimaryCaloPt->Write("MCrec_PrimaryCaloGamma_Pt",TObject::kOverwrite);

        // write reconstructed real photons
        if (fHistoGammaTrueConvPt) fHistoGammaTrueConvPt->Write("TrueConvGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueConvPtOrBin) fHistoGammaTrueConvPtOrBin->Write("TrueConvGamma_Pt_OriginalBinning",TObject::kOverwrite);
        if (fHistoGammaTrueCaloPt) fHistoGammaTrueCaloPt->Write("TrueCaloGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueCaloPtOrBin) fHistoGammaTrueCaloPtOrBin->Write("TrueCaloGamma_Pt_OriginalBinning",TObject::kOverwrite);

        // write primary and secondary reconstructed photons
        if (fHistoGammaTruePrimaryConvPt) fHistoGammaTruePrimaryConvPt->Write("TruePrimaryConvGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryConvPt) fHistoGammaTrueSecondaryConvPt->Write("TrueSecondaryConvGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt) fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt->Write("TrueSecondaryConvGammaFromXFromK0s_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt) fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt->Write("TrueSecondaryConvGammaFromXFromLambda_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryConvPtOrBin) fHistoGammaTrueSecondaryConvPtOrBin->Write("TrueSecondaryConvGamma_Pt_OriginalBinning",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryConvGammaFromXFromK0sPtOrBin) fHistoGammaTrueSecondaryConvGammaFromXFromK0sPtOrBin->Write("TrueSecondaryConvGammaFromXFromK0s_Pt_OriginalBinning",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPtOrBin) fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPtOrBin->Write("TrueSecondaryConvGammaFromXFromLambda_Pt_OriginalBinning",TObject::kOverwrite);
        
        if (fHistoGammaTruePrimaryCaloPt) fHistoGammaTruePrimaryCaloPt->Write("TruePrimaryCaloGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaTruePrimaryCaloPtOrBin) fHistoGammaTruePrimaryCaloPtOrBin->Write("TruePrimaryCaloGamma_Pt_OriginalBinning",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryCaloPt) fHistoGammaTrueSecondaryCaloPt->Write("TrueSecondaryCaloGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryCaloPtOrBin) fHistoGammaTrueSecondaryCaloPtOrBin->Write("TrueSecondaryCaloGamma_Pt_OriginalBinning",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryCaloFromK0sPt) fHistoGammaTrueSecondaryCaloFromK0sPt->Write("TrueSecondaryCaloGammaFromK0s_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryCaloFromK0sPtOrBin) fHistoGammaTrueSecondaryCaloFromK0sPtOrBin->Write("TrueSecondaryCaloGammaFromK0s_Pt_OriginalBinning",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryCaloFromLambdaPt) fHistoGammaTrueSecondaryCaloFromLambdaPt->Write("TrueSecondaryCaloGammaFromLambda_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryCaloFromLambdaPtOrBin) fHistoGammaTrueSecondaryCaloFromLambdaPtOrBin->Write("TrueSecondaryCaloGammaFromLambda_Pt_OriginalBinning",TObject::kOverwrite);
        
        // write fractions of secondary photons
        if (fHistoFracAllGammaToSec) fHistoFracAllGammaToSec->Write("FracAllGammaToSec",TObject::kOverwrite);
        if (fHistoFracAllGammaCaloToSec) fHistoFracAllGammaCaloToSec->Write("FracAllGammaCaloToSec",TObject::kOverwrite);
        if (fHistoFracAllGammaToSecFromXFromK0s) fHistoFracAllGammaToSecFromXFromK0s->Write("FracAllGammaToSecFromXFromK0s",TObject::kOverwrite);
        if (fHistoFracAllGammaCaloToSecFromK0s) fHistoFracAllGammaCaloToSecFromK0s->Write("FracAllGammaCaloToSecFromXFromK0s",TObject::kOverwrite);
        if (fHistoFracAllGammaToSecFromXFromLambda) fHistoFracAllGammaToSecFromXFromLambda->Write("FracAllGammaToSecFromXFromLambda",TObject::kOverwrite);
        if (fHistoFracAllGammaCaloToSecFromLambda) fHistoFracAllGammaCaloToSecFromLambda->Write("FracAllGammaCaloToSecFromXFromLambda",TObject::kOverwrite);
        if (fHistoFracAllGammaToSecOrBin) fHistoFracAllGammaToSecOrBin->Write("FracAllGammaToSecOriginalBinning",TObject::kOverwrite);
        if (fHistoFracAllGammaCaloToSecOrBin) fHistoFracAllGammaCaloToSecOrBin->Write("FracAllGammaCaloToSecOriginalBinning",TObject::kOverwrite);
        if (fHistoFracAllGammaToSecFromXFromK0sOrBin) fHistoFracAllGammaToSecFromXFromK0sOrBin->Write("FracAllGammaToSecFromXFromK0sOriginalBinning",TObject::kOverwrite);
        if (fHistoFracAllGammaCaloToSecFromK0sOrBin) fHistoFracAllGammaCaloToSecFromK0sOrBin->Write("FracAllGammaCaloToSecFromXFromK0sOriginalBinning",TObject::kOverwrite);
        if (fHistoFracAllGammaToSecFromXFromLambdaOrBin) fHistoFracAllGammaToSecFromXFromLambdaOrBin->Write("FracAllGammaToSecFromXFromLambdaOriginalBinning",TObject::kOverwrite);
        if (fHistoFracAllGammaCaloToSecFromLambdaOrBin) fHistoFracAllGammaCaloToSecFromLambdaOrBin->Write("FracAllGammaCaloToSecFromXFromLambdaOriginalBinning",TObject::kOverwrite);
        
        // write double counting histograms
        if (fEnableDCConv){
            TH1D*    fHistoTrueGammaConvDCPtRebinned    = NULL;
            if (fHistoTrueGammaConvDCPt){
                fHistoTrueGammaConvDCPt->Write("fHistoTrueGammaConvDCPt",TObject::kOverwrite);
                fHistoTrueGammaConvDCPtRebinned = (TH1D*)fHistoTrueGammaConvDCPt->Clone("fHistoTrueGammaConvDCPtRebinned");
                RebinSpectrum(fHistoTrueGammaConvDCPtRebinned);
            }
            if (fHistoTrueGammaConvDCR) fHistoTrueGammaConvDCR->Write("fHistoTrueGammaConvDCR",TObject::kOverwrite);
            if (fHistoTrueGammaConvDCRVSPt) fHistoTrueGammaConvDCRVSPt->Write("fHistoTrueGammaConvDCRVSPt",TObject::kOverwrite);
            if (fHistoTrueGammaConvMultipleCount)fHistoTrueGammaConvMultipleCount->Write("fHistoTrueGammaConvMultipleCount",TObject::kOverwrite);
            if (fHistoGammaTrueConvPtOrBin && fHistoTrueGammaConvDCPt ){
                TH1D* fHistoRatioDCGammaConv = (TH1D*)fHistoTrueGammaConvDCPt->Clone("fHistoRatioDCGammaConv");
                fHistoRatioDCGammaConv->Divide(fHistoRatioDCGammaConv,fHistoGammaTrueConvPtOrBin,1,1,"B" );
                fHistoRatioDCGammaConv->Write("FractionDCGammaConv",TObject::kOverwrite);
                
                if (fHistoTrueGammaConvDCPtRebinned){
                    TH1D* fHistoRatioDCGammaConvRebinned = (TH1D*)fHistoTrueGammaConvDCPtRebinned->Clone("fHistoRatioDCGammaConvRebinned");
                    fHistoRatioDCGammaConvRebinned->Divide(fHistoRatioDCGammaConvRebinned,fHistoGammaTrueConvPt,1,1,"B" );
                    fHistoRatioDCGammaConvRebinned->Write("FractionDCGammaConvRebinned",TObject::kOverwrite);
                }    
            }    
        }    
        
        // write source splitting of input photons
        for(Int_t i = 0; i<7; i++)
            if (fHistoGammaMCDecayPt[i]) fHistoGammaMCDecayPt[i]->Write(Form("MC_DecayGamma%s_Pt",fDecays[i].Data()),TObject::kOverwrite);
        
        // write correction histograms
        // ---> Purity
        if (fHistoGammaMCPurity) fHistoGammaMCPurity->Write("GammaPurity_Pt",TObject::kOverwrite);
        if (fHistoGammaMCTruePurity) fHistoGammaMCTruePurity->Write("GammaTruePurity_Pt",TObject::kOverwrite);
        if (fHistoGammaMCTruePurityOrBin) fHistoGammaMCTruePurityOrBin->Write("GammaTruePurity_OriginalBinning_Pt",TObject::kOverwrite);
        if (fHistoGammaCaloMCPurity) fHistoGammaCaloMCPurity->Write("GammaCaloPurity_Pt",TObject::kOverwrite);
        if (fHistoGammaCaloMCTruePurity) fHistoGammaCaloMCTruePurity->Write("GammaCaloTruePurity_Pt",TObject::kOverwrite);
        if (fHistoGammaCaloMCTruePurityOrBin) fHistoGammaCaloMCTruePurityOrBin->Write("GammaCaloTruePurity_OriginalBinning_Pt",TObject::kOverwrite);
        
        // ---> Conversion probability
        if (fHistoGammaMCConvProb) fHistoGammaMCConvProb->Write("GammaConvProb_MCPt",TObject::kOverwrite);
        
        // ---> Reco efficiency
        if (fHistoGammaMCRecoEff) fHistoGammaMCRecoEff->Write("GammaRecoEff_Pt",TObject::kOverwrite);
        if (fHistoGammaMCPrimaryRecoEff) fHistoGammaMCPrimaryRecoEff->Write("GammaPrimaryRecoEff_Pt",TObject::kOverwrite);
        if (fHistoGammaMCPrimaryRecoEffMCPt) fHistoGammaMCPrimaryRecoEffMCPt->Write("GammaPrimaryRecoEff_MCPt",TObject::kOverwrite);
        if (fHistoGammaCaloMCRecoEff) fHistoGammaCaloMCRecoEff->Write("GammaCaloRecoEff_Pt",TObject::kOverwrite);
        if (fHistoGammaCaloMCPrimaryRecoEff) fHistoGammaCaloMCPrimaryRecoEff->Write("GammaCaloPrimaryRecoEff_Pt",TObject::kOverwrite);
        if (fHistoGammaCaloMCPrimaryRecoEffMCPt) fHistoGammaCaloMCPrimaryRecoEffMCPt->Write("GammaCaloPrimaryRecoEff_MCPt",TObject::kOverwrite);
        
        // write response matices
        if (fHistoGammaTruePrimaryConv_recPt_MCPt_MC) fHistoGammaTruePrimaryConv_recPt_MCPt_MC->Write("TruePrimaryConvGamma_recPt_MCPt",TObject::kOverwrite);
        if (fHistoGammaTruePrimaryConv_recPt_MCPt_MC_Rebin) fHistoGammaTruePrimaryConv_recPt_MCPt_MC_Rebin->Write("TruePrimaryConvGamma_recPt_MCPt_Rebin",TObject::kOverwrite);
        if (fHistoGammaTruePrimaryCalo_recPt_MCPt_MC) fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->Write("TruePrimaryCaloGamma_recPt_MCPt",TObject::kOverwrite);
        if (fHistoGammaTruePrimaryCalo_recPt_MCPt_MC_Rebin) fHistoGammaTruePrimaryCalo_recPt_MCPt_MC_Rebin->Write("TruePrimaryCaloGamma_recPt_MCPt_Rebin",TObject::kOverwrite);
        
        // write total MC rec BG
        if (fHistoGammaMCBackground) fHistoGammaMCBackground->Write("MCrec_Background",TObject::kOverwrite);
        if (fHistoGammaCaloMCBackground) fHistoGammaCaloMCBackground->Write("MCrec_Calo_Background",TObject::kOverwrite);
        
        // write event histogram for MC
        fEventQuality->Write("NEvents",TObject::kOverwrite);

        // write pileup corrected histograms and DCA distributions for MC
        if(PileUpCorrection){
            for (Int_t catIter = 0; catIter < 4; catIter++) fMCrecGammaPtDCAzBins[catIter][0]->Write(Form("MCrec_GammaPtDCAzBin_Full_%s", categoryName[catIter].Data()),TObject::kOverwrite);
            fMCrecGammaPtPileUpAllCat->Write("MCrec_ConvGamma_Pt_PileUp_AllCatComb",TObject::kOverwrite);
            fMCrecGammaPtPileUp->Write("MCrec_ConvGamma_Pt_PileUp",TObject::kOverwrite);
            fHistoFracAllGammaToSecPileUp->Write("FracAllGammaToSecPileUp",TObject::kOverwrite);
            fHistoFracAllGammaToSecFromXFromK0sPileUp->Write("FracAllGammaToSecFromXFromK0sPileUp",TObject::kOverwrite);
            fHistoGammaMCPurityPileUp->Write("GammaPurity_PileUp_Pt",TObject::kOverwrite);
            fHistoGammaMCrecPrimaryConvPtPileUp->Write("MCrec_PrimaryConvGamma_PtPileUp",TObject::kOverwrite);
            fHistoGammaMCTruePurityPileUp->Write("GammaTruePurity_PileUp_Pt",TObject::kOverwrite);
            fHistoGammaMCRecoEffPileUp->Write("GammaRecoEff_PileUp_Pt",TObject::kOverwrite);
            fHistoGammaMCPrimaryRecoEffPileUp->Write("GammaPrimaryRecoEff_PileUp_Pt",TObject::kOverwrite);
        }
    
        // write combintorial histos
        if (fHistoCombinatorialBackground) fHistoCombinatorialBackground->Write("ESD_TrueCombinatorial_Pt",TObject::kOverwrite);
        if (fEnablePCM){
            for(Int_t i = 0; i<nCombinatorics+1; i++)
                if (fHistoCombinatorialSpecies[i]) fHistoCombinatorialSpecies[i]->Write(Form("ESD_TrueComb%s_Pt",fCombinatorics[i].Data()),TObject::kOverwrite);
        }
        if (fEnableCalo && !fEnablePCM){
            for(Int_t i = 0; i<nContamination+1; i++)
                if (fHistoCombinatorialSpecies[i]) fHistoCombinatorialSpecies[i]->Write(Form("ESD_TrueComb%s_Pt",fContamination[i].Data()),TObject::kOverwrite);
        }
        
    Output2->Write();
    Output2->Close();
}

//**************************************************************************************************
//************ Routine for saving histograms concerning pileup determination for conversion ********
//************ photons mainly for data                                                  ************
//**************************************************************************************************
void SaveDCAHistos(Int_t isMC, TString fCutID, TString fPrefix3){
    const char* nameOutput                              = Form("%s/%s/%s_%s_GammaConvV1DCAHistogramms%s_%s.root",fCutSelection.Data(), fEnergyFlag.Data(), fPrefix.Data(), fPrefix3.Data(),
                                                       fPeriodFlag.Data(), fCutID.Data());
    cout << "INFO: writing into: " << nameOutput << endl;
    TFile *Output1                                      = new TFile(nameOutput,"RECREATE");

        // write ratio of spectra with and without pileup
        if (!isMC) {
            
            for (Int_t i = 0; i < 3; i++) {
                fESDGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat[i]->Write(Form("ESD_ConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_AllCatComb_%s", backgroundExtractionMethod[i].Data()), TObject::kOverwrite);
                if (fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat[i])
                    fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat[i]->Write(Form("ESD_ConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_Fit_AllCatComb_%s", backgroundExtractionMethod[i].Data()), TObject::kOverwrite);
                fESDGammaPileUpCorrFactorAllCat[i]->Write(Form("ESD_ConvGamma_PileUpCorFactor_AllCatComb_%s", backgroundExtractionMethod[i].Data()), TObject::kOverwrite);
                fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[i]->Write(Form("ESD_ConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_%s", backgroundExtractionMethod[i].Data()), TObject::kOverwrite);
                if (fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinning[i])
                    fESDGammaPtRatioWithWithoutPileUpFitDCAzDistBinning[i]->Write(Form("ESD_ConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_Fit_%s", backgroundExtractionMethod[i].Data()), TObject::kOverwrite);
                if (i == 0)
                    fESDGammaPileUpCorrFactor[i]->Write("ESD_ConvGamma_PileUpCorFactor", TObject::kOverwrite);      // standard background estimation method
                else
                    fESDGammaPileUpCorrFactor[i]->Write(Form("ESD_ConvGamma_PileUpCorFactor_%s", backgroundExtractionMethod[i].Data()), TObject::kOverwrite);
            }
            
            for (Int_t catIter = 0; catIter < 4; catIter++) {
                if (catIter == 0) {
                    fESDGammaPtDCAz[5]->Write(fESDGammaPtDCAz[5]->GetName(), TObject::kOverwrite);
                    TH2D *ESDGammaPtDCAzPerEvent                                = (TH2D*)fESDGammaPtDCAz[5]->Clone(Form("%s_perEvent",fESDGammaPtDCAz[5]->GetName()));
                    ESDGammaPtDCAzPerEvent->Scale(1./fNEvnt);
                    ESDGammaPtDCAzPerEvent->Write(ESDGammaPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                } else {
                    fESDGammaPtDCAz[catIter]->Write(fESDGammaPtDCAz[catIter]->GetName(), TObject::kOverwrite);
                    TH2D *ESDGammaPtDCAzPerEvent                                = (TH2D*)fESDGammaPtDCAz[catIter]->Clone(Form("%s_perEvent",fESDGammaPtDCAz[catIter]->GetName()));
                    ESDGammaPtDCAzPerEvent->Scale(1./fNEvnt);
                    ESDGammaPtDCAzPerEvent->Write(ESDGammaPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                }
                
                Int_t ptBin;
                for (Int_t i = 0; i < 3; i++) {
                    if (i == 0)         ptBin                                   = 0;
                    else if (i == 1)    ptBin                                   = exemplaryLowPtBin;
                    else                ptBin                                   = exemplaryHighPtBin;
                    
                    fESDGammaPtDCAzBins[catIter][ptBin]->Write(fESDGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *ESDGammaPtDCAzBinsPerEvent                            = (TH1D*)fESDGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fESDGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    ESDGammaPtDCAzBinsPerEvent->Scale(1./fNEvnt);
                    ESDGammaPtDCAzBinsPerEvent->Write(ESDGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    fESDGammaPtDCAzBinsBack[catIter][ptBin][0]->Write(fESDGammaPtDCAzBinsBack[catIter][ptBin][0]->GetName(), TObject::kOverwrite);
                    TH1D *ESDGammaPtDCAzBinsBackPerEvent                        = (TH1D*)fESDGammaPtDCAzBinsBack[catIter][ptBin][0]->Clone(Form("%s_perEvent",fESDGammaPtDCAzBinsBack[catIter][ptBin][0]->GetName()));
                    ESDGammaPtDCAzBinsBackPerEvent->Scale(1./fNEvnt);
                    ESDGammaPtDCAzBinsBackPerEvent->Write(ESDGammaPtDCAzBinsBackPerEvent->GetName(), TObject::kOverwrite);
                    
                    fESDSubGammaPtDCAzBins[catIter][ptBin][0]->Write(fESDSubGammaPtDCAzBins[catIter][ptBin][0]->GetName(), TObject::kOverwrite);
                    TH1D *ESDGammaDCAzAllSubperEvent                            = (TH1D*)fESDSubGammaPtDCAzBins[catIter][ptBin][0]->Clone(Form("%s_perEvent",fESDSubGammaPtDCAzBins[catIter][ptBin][0]->GetName()));
                    ESDGammaDCAzAllSubperEvent->Scale(1./fNEvnt);
                    ESDGammaDCAzAllSubperEvent->Write(ESDGammaDCAzAllSubperEvent->GetName(), TObject::kOverwrite);
                }
            }
        } else {
            
            fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat->Write("MCrec_ConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_AllCatComb", TObject::kOverwrite);
            if(fMCrecGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat)
                fMCrecGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat->Write("MCrec_ConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_Fit_AllCatComb", TObject::kOverwrite);
            fMCrecGammaPileUpCorrFactorAllCat->Write("MCrec_ConvGamma_PileUpCorrFactor_AllCatComb", TObject::kOverwrite);
            fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinning->Write("MCrec_ConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning", TObject::kOverwrite);
            if(fMCrecGammaPtRatioWithWithoutPileUpFitDCAzDistBinning)
                fMCrecGammaPtRatioWithWithoutPileUpFitDCAzDistBinning->Write("MCrec_ConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_Fit", TObject::kOverwrite);
            fMCrecGammaPileUpCorrFactor->Write("MCrec_ConvGamma_PileUpCorrFactor", TObject::kOverwrite);

            fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat->Write("ESD_TruePrimaryConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_AllCatComb", TObject::kOverwrite);
            if(fTruePrimaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat)
                fTruePrimaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat->Write("ESD_TruePrimaryConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_Fit_AllCatComb", TObject::kOverwrite);
            fTruePrimaryConvGammaPileUpCorrFactorAllCat->Write("ESD_TruePrimaryConvGamma_PileUpCorrFactor_AllCatComb", TObject::kOverwrite);
            fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning->Write("ESD_TruePrimaryConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning", TObject::kOverwrite);
            if(fTruePrimaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning)
                fTruePrimaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning->Write("ESD_TruePrimaryConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_Fit", TObject::kOverwrite);
            fTruePrimaryConvGammaPileUpCorrFactor->Write("ESD_TruePrimaryConvGamma_PileUpCorrFactor", TObject::kOverwrite);

            fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat->Write("ESD_TrueSecondaryConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_AllCatComb", TObject::kOverwrite);
            if(fTrueSecondaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat)
                fTrueSecondaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat->Write("ESD_TrueSecondaryConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_Fit_AllCatComb", TObject::kOverwrite);
            fTrueSecondaryConvGammaPileUpCorrFactorAllCat->Write("ESD_TrueSecondaryConvGamma_PileUpCorrFactor_AllCatComb", TObject::kOverwrite);
            fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning->Write("ESD_TrueSecondaryConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning", TObject::kOverwrite);
            if(fTrueSecondaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning)
                fTrueSecondaryConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning->Write("ESD_TrueSecondaryConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_Fit", TObject::kOverwrite);
            fTrueSecondaryConvGammaPileUpCorrFactor->Write("ESD_TrueSecondaryConvGamma_PileUpCorrFactor", TObject::kOverwrite);

            fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat->Write("ESD_TrueSecondaryFromXFromK0sConvGamma_Pt__Ratio_WithWithoutPileUp_DCAzDistBinning_AllCatComb", TObject::kOverwrite);
            if(fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat)
                fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinningAllCat->Write("ESD_TrueSecondaryFromXFromK0sConvGamma_Pt__Ratio_WithWithoutPileUp_DCAzDistBinning_Fit_AllCatComb", TObject::kOverwrite);
            fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactorAllCat->Write("ESD_TrueSecondaryFromXFromK0sConvGamma_PileUpCorrFactor_AllCatComb", TObject::kOverwrite);
            fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinning->Write("ESD_TrueSecondaryFromXFromK0sConvGamma_Pt__Ratio_WithWithoutPileUp_DCAzDistBinning", TObject::kOverwrite);
            if(fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning)
                fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpFitDCAzDistBinning->Write("ESD_TrueSecondaryFromXFromK0sConvGamma_Pt__Ratio_WithWithoutPileUp_DCAzDistBinning_Fit", TObject::kOverwrite);
            fTrueSecondaryFromXFromK0sConvGammaPileUpCorrFactor->Write("ESD_TrueSecondaryFromXFromK0sConvGamma_PileUpCorrFactor", TObject::kOverwrite);

            for (Int_t catIter = 0; catIter < 4; catIter++) {
                if (catIter == 0) {
                    fESDGammaPtDCAz[5]->Write("MCrec_GammaPtDCAz_all", TObject::kOverwrite);
                    TH2D *ESDGammaPtDCAzPerEvent                                = (TH2D*)fESDGammaPtDCAz[5]->Clone("MCrec_GammaPtDCAz_all_perEvent");
                    ESDGammaPtDCAzPerEvent->Scale(1./fNEvnt);
                    ESDGammaPtDCAzPerEvent->Write(ESDGammaPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTruePrimaryPhotonPtDCAz[5]->Write(fTruePrimaryPhotonPtDCAz[5]->GetName(), TObject::kOverwrite);
                    TH2D *TruePrimaryPhotonPtDCAzPerEvent                       = (TH2D*)fTruePrimaryPhotonPtDCAz[5]->Clone(Form("%s_perEvent",fTruePrimaryPhotonPtDCAz[5]->GetName()));
                    TruePrimaryPhotonPtDCAzPerEvent->Scale(1./fNEvnt);
                    TruePrimaryPhotonPtDCAzPerEvent->Write(TruePrimaryPhotonPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueSecondaryPhotonPtDCAz[5]->Write(fTrueSecondaryPhotonPtDCAz[5]->GetName(), TObject::kOverwrite);
                    TH2D *TrueSecondaryPhotonPtDCAzPerEvent                     = (TH2D*)fTrueSecondaryPhotonPtDCAz[5]->Clone(Form("%s_perEvent",fTrueSecondaryPhotonPtDCAz[5]->GetName()));
                    TrueSecondaryPhotonPtDCAzPerEvent->Scale(1./fNEvnt);
                    TrueSecondaryPhotonPtDCAzPerEvent->Write(TrueSecondaryPhotonPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueSecondaryPhotonFromXFromK0sPtDCAz[5]->Write(fTrueSecondaryPhotonFromXFromK0sPtDCAz[5]->GetName(), TObject::kOverwrite);
                    TH2D *TrueSecondaryPhotonFromXFromK0sPtDCAzPerEvent         = (TH2D*)fTrueSecondaryPhotonFromXFromK0sPtDCAz[5]->Clone(Form("%s_perEvent",fTrueSecondaryPhotonFromXFromK0sPtDCAz[5]->GetName()));
                    TrueSecondaryPhotonFromXFromK0sPtDCAzPerEvent->Scale(1./fNEvnt);
                    TrueSecondaryPhotonFromXFromK0sPtDCAzPerEvent->Write(TrueSecondaryPhotonFromXFromK0sPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                } else {
                    fESDGammaPtDCAz[catIter]->Write(Form("MCrec_GammaPtDCAz_%s", categoryName[catIter].Data()), TObject::kOverwrite);
                    TH2D *ESDGammaPtDCAzPerEvent                                = (TH2D*)fESDGammaPtDCAz[catIter]->Clone(Form("MCrec_GammaPtDCAz_%s_perEvent",categoryName[catIter].Data()));
                    ESDGammaPtDCAzPerEvent->Scale(1./fNEvnt);
                    ESDGammaPtDCAzPerEvent->Write(ESDGammaPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTruePrimaryPhotonPtDCAz[catIter]->Write(fTruePrimaryPhotonPtDCAz[catIter]->GetName(), TObject::kOverwrite);
                    TH2D *TruePrimaryPhotonPtDCAzPerEvent                       = (TH2D*)fTruePrimaryPhotonPtDCAz[catIter]->Clone(Form("%s_perEvent",fTruePrimaryPhotonPtDCAz[catIter]->GetName()));
                    TruePrimaryPhotonPtDCAzPerEvent->Scale(1./fNEvnt);
                    TruePrimaryPhotonPtDCAzPerEvent->Write(TruePrimaryPhotonPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueSecondaryPhotonPtDCAz[catIter]->Write(fTrueSecondaryPhotonPtDCAz[catIter]->GetName(), TObject::kOverwrite);
                    TH2D *TrueSecondaryPhotonPtDCAzPerEvent                     = (TH2D*)fTrueSecondaryPhotonPtDCAz[catIter]->Clone(Form("%s_perEvent",fTrueSecondaryPhotonPtDCAz[catIter]->GetName()));
                    TrueSecondaryPhotonPtDCAzPerEvent->Scale(1./fNEvnt);
                    TrueSecondaryPhotonPtDCAzPerEvent->Write(TrueSecondaryPhotonPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueSecondaryPhotonFromXFromK0sPtDCAz[catIter]->Write(fTrueSecondaryPhotonFromXFromK0sPtDCAz[catIter]->GetName(), TObject::kOverwrite);
                    TH2D *TrueSecondaryPhotonFromXFromK0sPtDCAzPerEvent         = (TH2D*)fTrueSecondaryPhotonFromXFromK0sPtDCAz[catIter]->Clone(Form("%s_perEvent",fTrueSecondaryPhotonFromXFromK0sPtDCAz[catIter]->GetName()));
                    TrueSecondaryPhotonFromXFromK0sPtDCAzPerEvent->Scale(1./fNEvnt);
                    TrueSecondaryPhotonFromXFromK0sPtDCAzPerEvent->Write(TrueSecondaryPhotonFromXFromK0sPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                }
                
                Int_t ptBin;
                for (Int_t i = 0; i < 3; i++) {
                    if (i == 0)         ptBin                                   = 0;
                    else if (i == 1)    ptBin                                   = exemplaryLowPtBin;
                    else                ptBin                                   = exemplaryHighPtBin;

                    // MCrec
                    fMCrecGammaPtDCAzBins[catIter][ptBin]->Write(fMCrecGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *MCrecGammaPtDCAzBinsPerEvent                          = (TH1D*)fMCrecGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fMCrecGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    MCrecGammaPtDCAzBinsPerEvent->Scale(1./fNEvnt);
                    MCrecGammaPtDCAzBinsPerEvent->Write(MCrecGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);

                    fMCrecGammaPtDCAzBinsBack[catIter][ptBin]->Write(fMCrecGammaPtDCAzBinsBack[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *MCrecGammaPtDCAzBinsBackPerEvent                      = (TH1D*)fMCrecGammaPtDCAzBinsBack[catIter][ptBin]->Clone(Form("%s_perEvent",fMCrecGammaPtDCAzBinsBack[catIter][ptBin]->GetName()));
                    MCrecGammaPtDCAzBinsBackPerEvent->Scale(1./fNEvnt);
                    MCrecGammaPtDCAzBinsBackPerEvent->Write(MCrecGammaPtDCAzBinsBackPerEvent->GetName(), TObject::kOverwrite);
                    
                    fMCrecSubGammaPtDCAzBins[catIter][ptBin]->Write(fMCrecSubGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *MCrecSubGammaPtDCAzBinsPerEvent                       = (TH1D*)fMCrecSubGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fMCrecSubGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    MCrecSubGammaPtDCAzBinsPerEvent->Scale(1./fNEvnt);
                    MCrecSubGammaPtDCAzBinsPerEvent->Write(MCrecSubGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    // true gamma
                    fTrueGammaPtDCAzBins[catIter][ptBin]->Write(fTrueGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TrueGammaPtDCAzBinsPerEvent                       = (TH1D*)fTrueGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTrueGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    TrueGammaPtDCAzBinsPerEvent->Scale(1./fNEvnt);
                    TrueGammaPtDCAzBinsPerEvent->Write(TrueGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueSubGammaPtDCAzBins[catIter][ptBin]->Write(fTrueSubGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TrueSubGammaPtDCAzBinsPerEvent                       = (TH1D*)fTrueSubGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTrueSubGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    TrueSubGammaPtDCAzBinsPerEvent->Scale(1./fNEvnt);
                    TrueSubGammaPtDCAzBinsPerEvent->Write(TrueSubGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueBackgroundPtDCAzBins[catIter][ptBin]->Write(fTrueBackgroundPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TrueBackgroundPtDCAzBinsPerEvent                       = (TH1D*)fTrueBackgroundPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTrueBackgroundPtDCAzBins[catIter][ptBin]->GetName()));
                    TrueBackgroundPtDCAzBinsPerEvent->Scale(1./fNEvnt);
                    TrueBackgroundPtDCAzBinsPerEvent->Write(TrueBackgroundPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);                    
                    
                    // true primary
                    fTruePrimaryGammaPtDCAzBins[catIter][ptBin]->Write(fTruePrimaryGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TruePrimaryGammaPtDCAzBinsPerEvent                    = (TH1D*)fTruePrimaryGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTruePrimaryGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    TruePrimaryGammaPtDCAzBinsPerEvent->Scale(1./fNEvnt);
                    TruePrimaryGammaPtDCAzBinsPerEvent->Write(TruePrimaryGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTruePrimarySubGammaPtDCAzBins[catIter][ptBin]->Write(fTruePrimarySubGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TruePrimarySubGammaPtDCAzBinsPerEvent                 = (TH1D*)fTruePrimarySubGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTruePrimarySubGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    TruePrimarySubGammaPtDCAzBinsPerEvent->Scale(1./fNEvnt);
                    TruePrimarySubGammaPtDCAzBinsPerEvent->Write(TruePrimarySubGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    // true secondary
                    fTrueSecondaryGammaPtDCAzBins[catIter][ptBin]->Write(fTrueSecondaryGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TrueSecondaryGammaPtDCAzBinsPerEvent                  = (TH1D*)fTrueSecondaryGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTrueSecondaryGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    TrueSecondaryGammaPtDCAzBinsPerEvent->Scale(1./fNEvnt);
                    TrueSecondaryGammaPtDCAzBinsPerEvent->Write(TrueSecondaryGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueSecondarySubGammaPtDCAzBins[catIter][ptBin]->Write(fTrueSecondarySubGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TrueSecondarySubGammaPtDCAzBinsPerEvent               = (TH1D*)fTrueSecondarySubGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTrueSecondarySubGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    TrueSecondarySubGammaPtDCAzBinsPerEvent->Scale(1./fNEvnt);
                    TrueSecondarySubGammaPtDCAzBinsPerEvent->Write(TrueSecondarySubGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    // true secondary from X from K0s
                    fTrueSecondaryGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->Write(fTrueSecondaryGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TrueSecondaryGammaFromXFromK0sPtDCAzBinsPerEvent      = (TH1D*)fTrueSecondaryGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTrueSecondaryGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->GetName()));
                    TrueSecondaryGammaFromXFromK0sPtDCAzBinsPerEvent->Scale(1./fNEvnt);
                    TrueSecondaryGammaFromXFromK0sPtDCAzBinsPerEvent->Write(TrueSecondaryGammaFromXFromK0sPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueSecondarySubGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->Write(fTrueSecondarySubGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TrueSecondarySubGammaFromXFromK0sPtDCAzBinsPerEvent   = (TH1D*)fTrueSecondarySubGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTrueSecondarySubGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->GetName()));
                    TrueSecondarySubGammaFromXFromK0sPtDCAzBinsPerEvent->Scale(1./fNEvnt);
                    TrueSecondarySubGammaFromXFromK0sPtDCAzBinsPerEvent->Write(TrueSecondarySubGammaFromXFromK0sPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);

                }
            }
        }

    Output1->Write();
    Output1->Close();
}

//**************************************************************************************************
//************* Routine to produce additional DCA plots ********************************************
//**************************************************************************************************
void PlotAdditionalDCAz(Int_t isMC, TString fCutID){

    gStyle->SetOptTitle(0);

    TString nameInputMC                                 = Form("%s/%s/%s_MC_GammaConvV1DCAHistogramms%s_%s.root",fCutSelection.Data(), fEnergyFlag.Data(), fPrefix.Data(), fPeriodFlag.Data(), fCutID.Data());
    TFile *InputMC                                      = new TFile(nameInputMC);
    if (InputMC->IsZombie())
        InputMC                                         = NULL;

    TString nameInputData                               = Form("%s/%s/%s_data_GammaConvV1DCAHistogramms%s_%s.root",fCutSelection.Data(),fEnergyFlag.Data(),fPrefix.Data(),fPeriodFlag.Data(),fCutID.Data());
    TFile *InputData                                    = new TFile(nameInputData);
    if (InputData->IsZombie()) 
        InputData                                       = NULL;
    
    TH1D *ESDGammaDCAzBack;
    TH1D *ESDGammaDCAzBackperEvent;
    TH1D *ESDGammaDCAzAll;
    TH1D *ESDGammaDCAzAllperEvent;
    TH1D *ESDGammaDCAzAllSub;
    TH1D *ESDGammaDCAzAllSubperEvent;
    
    // Plot Data
    if(!isMC && InputData){
        for (Int_t catIter = 0; catIter < 4; catIter++) {
            
            ESDGammaDCAzAll                             = (TH1D*)InputData->Get(Form("ESD_GammaPtDCAzBin_Full_%s", categoryName[catIter].Data()));
            ESDGammaDCAzAllperEvent                     = (TH1D*)InputData->Get(Form("ESD_GammaPtDCAzBin_Full_%s_perEvent", categoryName[catIter].Data()));
            ESDGammaDCAzBack                            = (TH1D*)InputData->Get(Form("ESD_GammaPtDCAzBackBin_Full_%s_%s", categoryName[catIter].Data(), backgroundExtractionMethod[0].Data()));
            ESDGammaDCAzBackperEvent                    = (TH1D*)InputData->Get(Form("ESD_GammaPtDCAzBackBin_Full_%s_%s_perEvent", categoryName[catIter].Data(), backgroundExtractionMethod[0].Data()));
            ESDGammaDCAzAllSub                          = (TH1D*)InputData->Get(Form("ESD_SubGammaPtDCAzBin_Full_%s_%s", categoryName[catIter].Data(), backgroundExtractionMethod[0].Data()));
            ESDGammaDCAzAllSubperEvent                  = (TH1D*)InputData->Get(Form("ESD_SubGammaPtDCAzBin_Full_%s_%s_perEvent", categoryName[catIter].Data(), backgroundExtractionMethod[0].Data()));
            
            DrawGammaSetMarker(ESDGammaDCAzAll, 23, 1.0, kBlack, kBlack);
            DrawGammaSetMarker(ESDGammaDCAzAllperEvent, 23, 1.0, kBlack, kBlack);
            DrawGammaSetMarker(ESDGammaDCAzBack, 23, 1.0, kRed, kRed);
            DrawGammaSetMarker(ESDGammaDCAzBackperEvent, 23, 1.0, kRed, kRed);
            DrawGammaSetMarker(ESDGammaDCAzAllSub, 23, 1.0, kOrange+5, kOrange+5);
            DrawGammaSetMarker(ESDGammaDCAzAllSubperEvent, 23, 1.0,  kOrange+5, kOrange+5);
            
            TCanvas *canvasDCAzData                     = GetAndSetCanvas("canvasDCAzData");
            canvasDCAzData->SetLogy();
            
            SetHistogramm(ESDGammaDCAzAll,"DCA z (cm)","counts");
            SetHistogramm(ESDGammaDCAzAllperEvent,"DCA z (cm)","counts per event");
            SetHistogramm(ESDGammaDCAzBack,"DCA z (cm)","counts");
            SetHistogramm(ESDGammaDCAzBackperEvent,"DCA z (cm)","counts per event");
            SetHistogramm(ESDGammaDCAzAllSub,"DCA z (cm)","counts");
            SetHistogramm(ESDGammaDCAzAllSubperEvent,"DCA z (cm)","counts per event");
            
            ESDGammaDCAzAll->DrawCopy();
            ESDGammaDCAzBack->DrawCopy("same");
            ESDGammaDCAzAllSub->DrawCopy("same");

            TLegend* legendDCAZData                     = GetAndSetLegend(0.6,0.75,3,1);
            legendDCAZData->AddEntry(ESDGammaDCAzAll,"DCA z","lp");
            legendDCAZData->AddEntry(ESDGammaDCAzBack,"Estimated Pile-Up Bckg","lp");
            legendDCAZData->AddEntry(ESDGammaDCAzAllSub,"Subtracted Pile-Up Bckg","lp");
            legendDCAZData->Draw();
            
            canvasDCAzData->Print(Form("%s/%s_data_ESD_DCAz_%s_%s.%s",fOutputDir.Data(),fPrefix.Data(),categoryName[catIter].Data(),fCutSelection.Data(),fSuffix.Data()));
            delete canvasDCAzData;
            
            TCanvas *canvasDCAzDataPerEvent             = GetAndSetCanvas("canvasDCAzDataPerEvent");
            canvasDCAzDataPerEvent->SetLogy();
            
            ESDGammaDCAzAllperEvent->DrawCopy();
            ESDGammaDCAzBackperEvent->DrawCopy("same");
            ESDGammaDCAzAllSubperEvent->DrawCopy("same");
            
            legendDCAZData->Draw();
            
            canvasDCAzDataPerEvent->Print(Form("%s/%s_data_ESD_DCAz_PerEvent_%s_%s.%s",fOutputDir.Data(),fPrefix.Data(),categoryName[catIter].Data(),fCutSelection.Data(),fSuffix.Data()));
            delete legendDCAZData;
        }
    }
    
    TH1D *MCrecGammaDCAzAll;
    TH1D *MCrecGammaDCAzAllperEvent;
    TH1D *MCrecGammaDCAzBack;
    TH1D *MCrecGammaDCAzBackperEvent;
    TH1D *TruePrimaryGammaDCAzAll;
    TH1D *TruePrimaryGammaDCAzAllperEvent;
    TH1D *TrueSecondaryGammaDCAzAll;
    TH1D *TrueSecondaryGammaDCAzAllperEvent;
    TH1D *TrueSecondaryFromXFromK0sGammaDCAzAll;
    TH1D *TrueSecondaryFromXFromK0sGammaDCAzAllperEvent;
    TH1D *MCrecGammaBackgroundDCAzAll;
    TH1D *MCrecGammaBackgroundDCAzAllperEvent;

    if(isMC && InputMC){
        for (Int_t catIter = 0; catIter < 4; catIter++) {

            MCrecGammaDCAzAll                               = (TH1D*)InputMC->Get(Form("MCrec_GammaPtDCAzBin_Full_%s", categoryName[catIter].Data()));
            MCrecGammaDCAzAllperEvent                       = (TH1D*)InputMC->Get(Form("MCrec_GammaPtDCAzBin_Full_%s_perEvent", categoryName[catIter].Data()));
            MCrecGammaDCAzBack                              = (TH1D*)InputMC->Get(Form("MCrec_GammaPtDCAzBackBin_Full_%s", categoryName[catIter].Data()));
            MCrecGammaDCAzBackperEvent                      = (TH1D*)InputMC->Get(Form("MCrec_GammaPtDCAzBackBin_Full_%s_perEvent", categoryName[catIter].Data()));
            TruePrimaryGammaDCAzAll                         = (TH1D*)InputMC->Get(Form("ESD_TruePrimaryGammaPtDCAzBin_Full_%s", categoryName[catIter].Data()));
            TruePrimaryGammaDCAzAllperEvent                 = (TH1D*)InputMC->Get(Form("ESD_TruePrimaryGammaPtDCAzBin_Full_%s_perEvent", categoryName[catIter].Data()));
            TrueSecondaryGammaDCAzAll                       = (TH1D*)InputMC->Get(Form("ESD_TrueSecondaryGammaPtDCAzBin_Full_%s", categoryName[catIter].Data()));
            TrueSecondaryGammaDCAzAllperEvent               = (TH1D*)InputMC->Get(Form("ESD_TrueSecondaryGammaPtDCAzBin_Full_%s_perEvent", categoryName[catIter].Data()));
            TrueSecondaryFromXFromK0sGammaDCAzAll           = (TH1D*)InputMC->Get(Form("ESD_TrueSecondaryGammaFromXFromK0sPtDCAzBin_Full_%s", categoryName[catIter].Data()));
            TrueSecondaryFromXFromK0sGammaDCAzAllperEvent   = (TH1D*)InputMC->Get(Form("ESD_TrueSecondaryGammaFromXFromK0sPtDCAzBin_Full_%s_perEvent", categoryName[catIter].Data()));
            
            MCrecGammaBackgroundDCAzAll                     = (TH1D*)MCrecGammaDCAzAll->Clone("ESD_TrueGammaBackgroundDCAz_all");
            MCrecGammaBackgroundDCAzAll->Add(TruePrimaryGammaDCAzAll,-1);
            
            MCrecGammaBackgroundDCAzAllperEvent             = (TH1D*)MCrecGammaDCAzAllperEvent->Clone("ESD_TrueGammaBackgroundDCAz_all");
            MCrecGammaBackgroundDCAzAllperEvent->Add(TruePrimaryGammaDCAzAllperEvent,-1);

            TrueSecondaryGammaDCAzAll->Add(TrueSecondaryFromXFromK0sGammaDCAzAll,-1);
            TrueSecondaryGammaDCAzAllperEvent->Add(TrueSecondaryFromXFromK0sGammaDCAzAllperEvent,-1);

            DrawGammaSetMarker(MCrecGammaDCAzAll, 23, 1.0, kYellow+1, kYellow+1);
            DrawGammaSetMarker(MCrecGammaDCAzAllperEvent, 23, 1.0, kYellow+1, kYellow+1);
            DrawGammaSetMarker(MCrecGammaDCAzBack, 23, 1.0, kRed, kRed);
            DrawGammaSetMarker(MCrecGammaDCAzBackperEvent, 23, 1.0, kRed, kRed);
            DrawGammaSetMarker(TruePrimaryGammaDCAzAll, 23, 1.0, kBlue+1, kBlue+1);
            DrawGammaSetMarker(TruePrimaryGammaDCAzAllperEvent, 23, 1.0, kBlue+1, kBlue+1);
            DrawGammaSetMarker(TrueSecondaryGammaDCAzAll, 23, 1.0, kCyan+2, kCyan+2);
            DrawGammaSetMarker(TrueSecondaryGammaDCAzAllperEvent, 23, 1.0, kCyan+2, kCyan+2);
            DrawGammaSetMarker(TrueSecondaryFromXFromK0sGammaDCAzAll, 23, 1.0, kViolet+2, kViolet+2);
            DrawGammaSetMarker(TrueSecondaryFromXFromK0sGammaDCAzAllperEvent, 23, 1.0, kViolet+2, kViolet+2);
            DrawGammaSetMarker(MCrecGammaBackgroundDCAzAll, 23, 1.0, kGreen+2, kGreen+2);
            DrawGammaSetMarker(MCrecGammaBackgroundDCAzAllperEvent, 23, 1.0, kGreen+2, kGreen+2);
            
            TCanvas *canvasDCAzMC = GetAndSetCanvas("canvasDCAzMC");
            canvasDCAzMC->SetLogy();
            
            SetHistogramm(MCrecGammaDCAzAll,"DCA z (cm)","counts");
            SetHistogramm(MCrecGammaDCAzAllperEvent,"DCA z (cm)","counts per event");
            SetHistogramm(MCrecGammaDCAzBack,"DCA z (cm)","counts");
            SetHistogramm(MCrecGammaDCAzBackperEvent,"DCA z (cm)","counts per event");
            SetHistogramm(TruePrimaryGammaDCAzAll,"DCA z (cm)","counts");
            SetHistogramm(TruePrimaryGammaDCAzAllperEvent,"DCA z (cm)","counts per event");
            SetHistogramm(TrueSecondaryGammaDCAzAll,"DCA z (cm)","counts");
            SetHistogramm(TrueSecondaryGammaDCAzAllperEvent,"DCA z (cm)","counts per event");
            SetHistogramm(TrueSecondaryFromXFromK0sGammaDCAzAll,"DCA z (cm)","counts");
            SetHistogramm(TrueSecondaryFromXFromK0sGammaDCAzAllperEvent,"DCA z (cm)","counts per event");
            SetHistogramm(MCrecGammaBackgroundDCAzAll,"DCA z (cm)","counts");
            SetHistogramm(MCrecGammaBackgroundDCAzAllperEvent,"DCA z (cm)","counts per event");
            
            MCrecGammaDCAzAll->DrawCopy();
            MCrecGammaBackgroundDCAzAll->DrawCopy("same");
            MCrecGammaDCAzBack->DrawCopy("same");
            TruePrimaryGammaDCAzAll->DrawCopy("same");
            TrueSecondaryGammaDCAzAll->DrawCopy("same");
            TrueSecondaryFromXFromK0sGammaDCAzAll->DrawCopy("same");
            
            TLegend* legendDCAZMC = GetAndSetLegend(0.6,0.65,6,1);
            legendDCAZMC->AddEntry(MCrecGammaDCAzAll,"DCA z","lp");
            legendDCAZMC->AddEntry(MCrecGammaDCAzBack,"Estimated Pile-Up Bckg","lp");
            legendDCAZMC->AddEntry(MCrecGammaBackgroundDCAzAll,"Background + Secondary #gamma","lp");
            legendDCAZMC->AddEntry(TruePrimaryGammaDCAzAll,"True Primary #gamma","lp");
            legendDCAZMC->AddEntry(TrueSecondaryFromXFromK0sGammaDCAzAll,"True Secondary #gamma from K_{s}^{0}","lp");
            legendDCAZMC->AddEntry(TrueSecondaryGammaDCAzAll,"Additional Secondary #gamma","lp");
            legendDCAZMC->Draw();
            
            canvasDCAzMC->Print(Form("%s/%s_MC_MCrec_DCAz_%s_%s.%s",fOutputDir.Data(),fPrefix.Data(),categoryName[catIter].Data(),fCutSelection.Data(),fSuffix.Data()));
            delete canvasDCAzMC;
            
            TCanvas *canvasDCAzMCPerEvent                   = GetAndSetCanvas("canvasDCAzMCPerEvent");
            canvasDCAzMCPerEvent->SetLogy();
            
            TLegend* legendDCAZMCperEvent                   = GetAndSetLegend(0.6,0.65,5,1,"MC");
            TLegend* legendDCAZperEvent                     = GetAndSetLegend(0.15,0.85,1,1,"Data");
            
            if(InputData){
                ESDGammaDCAzAllSubperEvent                  = (TH1D*)InputData->Get(Form("ESD_SubGammaPtDCAzBin_Full_%s_%s_perEvent", categoryName[catIter].Data(), backgroundExtractionMethod[0].Data()));
                DrawGammaSetMarker(ESDGammaDCAzAllSubperEvent, 23, 1.0,  kBlack, kBlack);
                SetHistogramm(ESDGammaDCAzAllSubperEvent,"DCA z (cm)","counts per event",1e-8,5e-2);
                ESDGammaDCAzAllSubperEvent->DrawCopy("");
                MCrecGammaDCAzAllperEvent->DrawCopy("same");
                legendDCAZperEvent->AddEntry(ESDGammaDCAzAllSubperEvent,"DCA z (Pile-Up Subtracted)","lp");
            }
            else MCrecGammaDCAzAllperEvent->DrawCopy("");
            
            MCrecGammaBackgroundDCAzAllperEvent->DrawCopy("same");
            MCrecGammaDCAzBackperEvent->DrawCopy("same");
            TruePrimaryGammaDCAzAllperEvent->DrawCopy("same");
            TrueSecondaryGammaDCAzAllperEvent->DrawCopy("same");
            TrueSecondaryFromXFromK0sGammaDCAzAllperEvent->DrawCopy("same");
            
            legendDCAZMCperEvent->AddEntry(MCrecGammaDCAzAllperEvent,"DCA z","lp");
            legendDCAZMCperEvent->AddEntry(MCrecGammaDCAzBackperEvent,"Estimated Pile-Up Bckg","lp");
            legendDCAZMCperEvent->AddEntry(MCrecGammaBackgroundDCAzAllperEvent,"Background + Secondary #gamma","lp");
            legendDCAZMCperEvent->AddEntry(TruePrimaryGammaDCAzAllperEvent,"True Primary #gamma","lp");
            legendDCAZMCperEvent->AddEntry(TrueSecondaryFromXFromK0sGammaDCAzAllperEvent,"True Secondary #gamma from K_{s}^{0}","lp");
            legendDCAZMCperEvent->AddEntry(TrueSecondaryGammaDCAzAllperEvent,"Additional Secondary #gamma","lp");
            legendDCAZMCperEvent->Draw();
            legendDCAZperEvent->Draw();
            canvasDCAzMCPerEvent->Print(Form("%s/%s_MC_MCrec_DCAz_PerEvent_%s_%s.%s",fOutputDir.Data(),fPrefix.Data(),categoryName[catIter].Data(),fCutSelection.Data(),fSuffix.Data()));
            delete canvasDCAzMCPerEvent;
        }
    }
    
    TH1D *RatioWithWithoutPileUpData;
    TH1D *RatioWithWithoutPileUpMC;

    if(InputMC && InputData){
        RatioWithWithoutPileUpData                          = (TH1D*)InputData->Get(Form("ESD_ConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning_%s", backgroundExtractionMethod[0].Data()));
        RatioWithWithoutPileUpMC                            = (TH1D*)InputMC->Get(Form("MCrec_ConvGamma_Pt_Ratio_WithWithoutPileUp_DCAzDistBinning"));

        DrawGammaSetMarker(RatioWithWithoutPileUpData, 20, 1.0, kBlack, kBlack);
        DrawGammaSetMarker(RatioWithWithoutPileUpMC, 20, 1.0, kRed, kRed);

        TCanvas *canvasComparisonWithWithoutPileUp      = GetAndSetCanvas("canvasComparisonWithWithoutPileUp");
        TLegend* legendComparisonWithWithoutPileUp      = GetAndSetLegend(0.3,0.75,2.4,1);
        legendComparisonWithWithoutPileUp->AddEntry(RatioWithWithoutPileUpData,"#gamma / #gamma Pile-Up Cor. data","lp");
        legendComparisonWithWithoutPileUp->AddEntry(RatioWithWithoutPileUpMC,"#gamma / #gamma Pile-Up Cor. MC","lp");

        RatioWithWithoutPileUpData->DrawCopy();
        RatioWithWithoutPileUpMC->DrawCopy("same");

        legendComparisonWithWithoutPileUp->Draw();

        canvasComparisonWithWithoutPileUp->Print(Form("%s/%s_PileUpComparisonMCDate_%s.%s",fOutputDir.Data(),fPrefix.Data(),fCutSelection.Data(),fSuffix.Data()));
        delete canvasComparisonWithWithoutPileUp;
    }
    
    gStyle->SetOptTitle(1);
}

//**************************************************************************************************
//************* Routine to produce DCAz plots in pt bins *******************************************
//**************************************************************************************************
void PlotDCAzInPtBinsWithBack(TH1D** ESDGammaPtDCAzBins, TH1D** ESDGammaPtDCAzBinsBack,TH1D** ESDGammaPtDCAzBinsBackB, TString namePlot, TString nameCanvas, TString namePad, TString dateDummy, TString fMesonType,  Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, Double_t* fRangeBinsPt, TString fDecayChannel, Bool_t fMonteCarloInfo, TString textCent){
    
    cout << textCent.Data() << endl;
    TGaxis::SetMaxDigits(3);
    
    TCanvas * canvasDataSpectra                 = new TCanvas(nameCanvas.Data(),"",2800,1800);  // gives the page size
    canvasDataSpectra->SetTopMargin(0.02);
    canvasDataSpectra->SetBottomMargin(0.02);
    canvasDataSpectra->SetRightMargin(0.02);
    canvasDataSpectra->SetLeftMargin(0.02);
    
    TPad * padDataSpectra                       = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(1);
    padDataSpectra->Divide(fColumnPlot,fRowPlot);
    padDataSpectra->Draw();
    
    Double_t relWidthLogo;
    if (fMesonType.CompareTo("Pi0") == 0){
        relWidthLogo                            = 0.3;
    } else {
        relWidthLogo                            = 0.3;
    }
    Double_t padXWidth                          = 2800/fColumnPlot;
    Double_t padYWidth                          = 1800/fRowPlot;
    
    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    
    Int_t place                                 = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins+1;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt                        = fRangeBinsPt[iPt];
        Double_t endPt                          = fRangeBinsPt[iPt+1];
        
        place                                   = place + 1;                                    //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);
            
            TString textAlice                   = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo) {
                textEvents                      = "MC";
            } else {
                textEvents                      = "Data";
            }
            
            Double_t nPixels                    = 13;
            Double_t textHeight                 = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight                      = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight                      = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }
            
            Double_t startTextX                 = 0.1;
            Double_t startTextY                 = 0.9;
            Double_t differenceText             = textHeight*1.25;
            
            TLatex *alice                       = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate                   = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy                      = new TLatex(startTextX, (startTextY-2.25*differenceText), fCollisionSystem);
            TLatex *process                     = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess                  = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionProcess);
            TLatex *events                      = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvnt));
            
            alice->SetNDC();
            alice->SetTextColor(1);
            alice->SetTextSize(textHeight*1.3);
            alice->Draw();
            
            latexDate->SetNDC();
            latexDate->SetTextColor(1);
            latexDate->SetTextSize(textHeight);
            latexDate->Draw();
            
            energy->SetNDC();
            energy->SetTextColor(1);
            energy->SetTextSize(textHeight);
            energy->Draw();
            
            process->SetNDC();
            process->SetTextColor(1);
            process->SetTextSize(textHeight);
            process->Draw();
            
            detprocess->SetNDC();
            detprocess->SetTextColor(1);
            detprocess->SetTextSize(textHeight);
            detprocess->Draw();
            
            events->SetNDC();
            events->SetTextColor(1);
            events->SetTextSize(textHeight);
            events->Draw();
            
            TLegend* legendData                 = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+2.)*differenceText);
            legendData->SetTextSize(textHeight);
            legendData->SetTextFont(62);
            legendData->SetFillColor(0);
            legendData->SetFillStyle(0);
            legendData->SetLineWidth(0);
            legendData->SetLineColor(0);
            legendData->SetMargin(0.15);
            legendData->AddEntry(ESDGammaPtDCAzBins[iPt],ESDGammaPtDCAzBins[iPt]->GetName(),"l");
            legendData->AddEntry(ESDGammaPtDCAzBinsBack[iPt],ESDGammaPtDCAzBinsBack[iPt]->GetName(),"l");
            if(ESDGammaPtDCAzBinsBackB)legendData->AddEntry(ESDGammaPtDCAzBinsBackB[iPt],ESDGammaPtDCAzBinsBackB[iPt]->GetName(),"l");
            legendData->Draw();
        } else {
            padDataSpectra->cd(place);
            padDataSpectra->cd(place)->SetLogy();
            padDataSpectra->cd(place)->SetTopMargin(0.12);
            padDataSpectra->cd(place)->SetBottomMargin(0.15);
            padDataSpectra->cd(place)->SetRightMargin(0.02);
            int remaining                       = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);
            
            DrawDCAzHisto(  ESDGammaPtDCAzBins[iPt],
                          Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                          "DCA z (cm)", "dN/dDCA z",
                          -10,10,0);
            
            DrawDCAzHisto(  ESDGammaPtDCAzBinsBack[iPt],
                            Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                            "DCA z (cm)", "dN/dDCA z",
                            -10,10,1);
            
            if (ESDGammaPtDCAzBinsBackB) {
                DrawDCAzHisto(  ESDGammaPtDCAzBinsBackB[iPt],
                              Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                              "DCA z (cm)", "dN/dDCA z",
                              -10,10,2);
            }
        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;
}

void PlotDCAzInPtBinsWithBack(TH1D** ESDGammaPtDCAzBins, TH1D*** ESDGammaPtDCAzBinsBack, TH1D** ESDGammaPtDCAzBinsBackB, TString namePlot, TString nameCanvas, TString namePad, TString dateDummy, TString fMesonType, Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, Double_t* fRangeBinsPt, TString fDecayChannel, Bool_t fMonteCarloInfo, TString textCent){
    
    cout << textCent.Data() << endl;
    TGaxis::SetMaxDigits(3);
    
    TCanvas * canvasDataSpectra                 = new TCanvas(nameCanvas.Data(),"",2800,1800);  // gives the page size
    canvasDataSpectra->SetTopMargin(0.02);
    canvasDataSpectra->SetBottomMargin(0.02);
    canvasDataSpectra->SetRightMargin(0.02);
    canvasDataSpectra->SetLeftMargin(0.02);
    
    TPad * padDataSpectra                       = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(1);
    padDataSpectra->Divide(fColumnPlot,fRowPlot);
    padDataSpectra->Draw();
    
    Double_t relWidthLogo;
    if (fMesonType.CompareTo("Pi0") == 0){
        relWidthLogo                            = 0.3;
    } else {
        relWidthLogo                            = 0.3;
    }
    Double_t padXWidth                          = 2800/fColumnPlot;
    Double_t padYWidth                          = 1800/fRowPlot;
    
    cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
    
    Int_t place                                 = 0;
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins+1;iPt++){
        cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
        Double_t startPt                        = fRangeBinsPt[iPt];
        Double_t endPt                          = fRangeBinsPt[iPt+1];
        
        place                                   = place + 1;                                    //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataSpectra->cd(place);
            
            TString textAlice                   = "ALICE performance";
            TString textEvents;
            if(fMonteCarloInfo) {
                textEvents                      = "MC";
            } else {
                textEvents                      = "Data";
            }
            
            Double_t nPixels                    = 13;
            Double_t textHeight                 = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight                      = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight                      = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }
            
            Double_t startTextX                 = 0.1;
            Double_t startTextY                 = 0.9;
            Double_t differenceText             = textHeight*1.25;
            
            TLatex *alice                       = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
            TLatex *latexDate                   = new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
            TLatex *energy                      = new TLatex(startTextX, (startTextY-2.25*differenceText), fCollisionSystem);
            TLatex *process                     = new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
            TLatex *detprocess                  = new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionProcess);
            TLatex *events                      = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvnt));
            
            alice->SetNDC();
            alice->SetTextColor(1);
            alice->SetTextSize(textHeight*1.3);
            alice->Draw();
            
            latexDate->SetNDC();
            latexDate->SetTextColor(1);
            latexDate->SetTextSize(textHeight);
            latexDate->Draw();
            
            energy->SetNDC();
            energy->SetTextColor(1);
            energy->SetTextSize(textHeight);
            energy->Draw();
            
            process->SetNDC();
            process->SetTextColor(1);
            process->SetTextSize(textHeight);
            process->Draw();
            
            detprocess->SetNDC();
            detprocess->SetTextColor(1);
            detprocess->SetTextSize(textHeight);
            detprocess->Draw();
            
            events->SetNDC();
            events->SetTextColor(1);
            events->SetTextSize(textHeight);
            events->Draw();
            
            TLegend* legendData                 = GetAndSetLegend(startTextX, startTextY-12*differenceText, 4);
            legendData->SetTextSize(textHeight);
            legendData->SetTextFont(62);
            legendData->SetFillColor(0);
            legendData->SetFillStyle(0);
            legendData->SetLineWidth(0);
            legendData->SetLineColor(0);
            legendData->SetMargin(0.15);
            legendData->AddEntry(ESDGammaPtDCAzBins[iPt],ESDGammaPtDCAzBins[iPt]->GetName(),"l");
            for (Int_t i = 0; i < 3; i++)
                legendData->AddEntry(ESDGammaPtDCAzBinsBack[iPt][i],Form("background extraction %s", backgroundExtractionMethod[i].Data()),"l");
            if(ESDGammaPtDCAzBinsBackB)legendData->AddEntry(ESDGammaPtDCAzBinsBackB[iPt],ESDGammaPtDCAzBinsBackB[iPt]->GetName(),"l");
            legendData->Draw();
        } else {
            padDataSpectra->cd(place);
            padDataSpectra->cd(place)->SetLogy();
            padDataSpectra->cd(place)->SetTopMargin(0.12);
            padDataSpectra->cd(place)->SetBottomMargin(0.15);
            padDataSpectra->cd(place)->SetRightMargin(0.02);
            int remaining                       = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);
            
            DrawDCAzHisto(  ESDGammaPtDCAzBins[iPt],
                          Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                          "DCA z (cm)", "dN/dDCA z",
                          -10,10,0);
            
            for (Int_t i = 0; i < 3; i++) {
                DrawDCAzHisto(  ESDGammaPtDCAzBinsBack[iPt][i],
                              Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                              "DCA z (cm)", "dN/dDCA z",
                              -10,10,1,backgroundColor[i]);
            }
            
            if (ESDGammaPtDCAzBinsBackB) {
                DrawDCAzHisto(  ESDGammaPtDCAzBinsBackB[iPt],
                              Form("%3.2f GeV/#it{c} < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt),
                              "DCA z (cm)", "dN/dDCA z",
                              -10,10,2);
            }
        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;
}

// overloading PlotDCAzInPtBinsWithBack()
void PlotDCAzInPtBinsWithBack(TH1D** ESDGammaPtDCAzBins, TH1D** ESDGammaPtDCAzBinsBack,TH1D** ESDGammaPtDCAzBinsBackB, TString namePlot, TString nameCanvas, TString namePad, TString dateDummy, TString fMesonType, Int_t fStartBinPtRange, Int_t fNumberPtBins, Double_t* fRangeBinsPt, TString fDecayChannel, Bool_t fMonteCarloInfo, TString textCent) {
    
    Int_t nPads = fNumberPtBins + 2;
    
    Int_t nColumns  = 2;
    
    for (Int_t i = 0; i < nPads; i++) {
        if (((nColumns+1) * CalculateNumberOfRowsForDCAzPlots(nPads, nColumns+1) - nPads <= (nColumns) * CalculateNumberOfRowsForDCAzPlots(nPads, nColumns) - nPads) || (TMath::Abs(nColumns+1 - CalculateNumberOfRowsForDCAzPlots(nPads, nColumns+1)) < TMath::Abs(nColumns - CalculateNumberOfRowsForDCAzPlots(nPads, nColumns)))) {
            nColumns++;
        } else {
            break;
        }
    }
    
    Int_t nRows = CalculateNumberOfRowsForDCAzPlots(nPads, nColumns);
    
    PlotDCAzInPtBinsWithBack(ESDGammaPtDCAzBins, ESDGammaPtDCAzBinsBack,ESDGammaPtDCAzBinsBackB, namePlot, nameCanvas, namePad, dateDummy, fMesonType, nRows, nColumns, fStartBinPtRange, fNumberPtBins, fRangeBinsPt, fDecayChannel, fMonteCarloInfo, textCent);
}

void PlotDCAzInPtBinsWithBack(TH1D** ESDGammaPtDCAzBins, TH1D*** ESDGammaPtDCAzBinsBack,TH1D** ESDGammaPtDCAzBinsBackB, TString namePlot, TString nameCanvas, TString namePad, TString dateDummy, TString fMesonType, Int_t fStartBinPtRange, Int_t fNumberPtBins, Double_t* fRangeBinsPt, TString fDecayChannel, Bool_t fMonteCarloInfo, TString textCent) {
    
    Int_t nPads = fNumberPtBins + 2;
    
    Int_t nColumns  = 2;
    
    for (Int_t i = 0; i < nPads; i++) {
        if (((nColumns+1) * CalculateNumberOfRowsForDCAzPlots(nPads, nColumns+1) - nPads <= (nColumns) * CalculateNumberOfRowsForDCAzPlots(nPads, nColumns) - nPads) || (TMath::Abs(nColumns+1 - CalculateNumberOfRowsForDCAzPlots(nPads, nColumns+1)) < TMath::Abs(nColumns - CalculateNumberOfRowsForDCAzPlots(nPads, nColumns)))) {
            nColumns++;
        } else {
            break;
        }
    }
    
    Int_t nRows = CalculateNumberOfRowsForDCAzPlots(nPads, nColumns);
    
    PlotDCAzInPtBinsWithBack(ESDGammaPtDCAzBins, ESDGammaPtDCAzBinsBack,ESDGammaPtDCAzBinsBackB, namePlot, nameCanvas, namePad, dateDummy, fMesonType, nRows, nColumns, fStartBinPtRange, fNumberPtBins, fRangeBinsPt, fDecayChannel, fMonteCarloInfo, textCent);
}


//**************************************************************************************************
//******* Function to calculate number of rows for given number of bins and columns ****************
//**************************************************************************************************
Int_t CalculateNumberOfRowsForDCAzPlots(Int_t numberOfPads, Int_t numberOfColumns) {
    // this function returns the number of rows
    // for a given number of pads and columns,
    
    Int_t over = 0;
    Int_t rows = 0;
    
    for (Int_t i = 0; i < numberOfPads; i++) {
        if ((numberOfPads + over)%numberOfColumns != 0) {
            over++;
        } else if ((numberOfPads + over)%numberOfColumns == 0) {
            break;
        }
    }
    
    return rows = (numberOfPads + over)/numberOfColumns;
}


//****************************************************************************
//******* Function to draw DCAz histograms in pT-bins ************************
//****************************************************************************
void DrawDCAzHisto( TH1* histo1,
                    TString Title,
                    TString XTitle,
                    TString YTitle,
                    Float_t xMin,
                    Float_t xMax,
                    Int_t bck,
                    Color_t color) {
    
    histo1->GetXaxis()->SetRangeUser(xMin, xMax);
    
    if(XTitle.Length() > 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.Length() > 0){
        histo1->SetYTitle(YTitle.Data());
    }
    histo1->GetYaxis()->SetLabelSize(0.02);
    histo1->GetYaxis()->SetTitleSize(0.025);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(0.5);
    histo1->GetXaxis()->SetTitleSize(0.025);
    histo1->GetXaxis()->SetLabelSize(0.02);
    histo1->SetMarkerStyle(20);
    histo1->SetMarkerColor(1);
    histo1->SetLineColor(1);
    histo1->SetLineWidth(0.5);
    histo1->SetMarkerSize(0.5);
    histo1->SetTitleOffset(1.2,"xy");
    histo1->SetTitleSize(0.05,"xy");
    histo1->GetYaxis()->SetLabelSize(0.05);
    histo1->GetXaxis()->SetLabelSize(0.05);
    histo1->GetXaxis()->SetNdivisions(507,kTRUE);
    if( bck == 1 ){
        histo1->SetLineStyle(1);
        histo1->SetLineColor(color);
        histo1->SetMarkerColor(color);
        histo1->SetMarkerStyle(24);
        histo1->SetLineWidth(0.9);
        histo1->DrawCopy("hist,same");
    } else {
        if( bck == 2 ){
            histo1->DrawCopy("same");
        } else {
            if(Title.Length() > 0){
                histo1->SetTitle("");
            }
            histo1->DrawCopy("e1,p");
            if(Title.Length() > 0){
                TLatex *alice = new TLatex(0.1,0.95,Form("%s",Title.Data())); // Bo: this was
                alice->SetNDC();
                alice->SetTextColor(1);
                alice->SetTextSize(0.062);
                alice->Draw();
            }
        }
    }
}

//**************************************************************************************************
//************* Routine to produce fraction per category vs pt plots  ******************************
//**************************************************************************************************
void DrawFractionPerCat(TH1D** frac, TString fOutputDir, TString fPrefix, TString fPrefix2, TString fCutSelection, TString fSuffix) {
    
    TCanvas *canvas             = GetAndSetCanvas("canvas");
    Color_t markerColor[3]      = {kBlack, kRed, kBlue};
    
    TLegend* legend             = GetAndSetLegend(0.7,0.75,3,1);
    for (Int_t i=0; i<3; i++) {
        SetHistogramm(frac[i],"p_{T} (GeV/c)","#gamma_{cat i} / #gamma_{all cat}",0.0,1.0);
        DrawGammaSetMarker(frac[i], 20, 1.0, markerColor[i], markerColor[i]);
        frac[i]->Draw("same");
        legend->AddEntry(frac[i],Form("#gamma_{cat %i} / #gamma_{all cat}", i+1),"lp");
    }
    DrawGammaLines(0., frac[0]->GetXaxis()->GetBinUpEdge(frac[0]->GetNbinsX()), 1., 1., 0.5, kBlack);
    legend->Draw("same");
    
    canvas->Print(Form("%s/%s_%s_ESD_FractionPerCategory_%s.%s",fOutputDir.Data(),fPrefix.Data(),fPrefix2.Data(),fCutSelection.Data(),fSuffix.Data()));
    delete legend;
    delete canvas;
}

//****************************************************************************
//******* Fill array of invariant mass histograms in pT slices ***************
//****************************************************************************
void FillMassHistosArray(TH2D* fGammaGammaInvMassVSPtDummy) {
    TString fNameHistoGG                                        = "";
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        fNameHistoGG                                            = Form("Mapping_GconvG_InvMass_in_Pt_Bin%02d", iPt);

        if(fHistoGconvGInvMassPtGConvBin[iPt]!= NULL){
            delete fHistoGconvGInvMassPtGConvBin[iPt];
            fHistoGconvGInvMassPtGConvBin[iPt]                  = NULL;
        }
        fHistoGconvGInvMassPtGConvBin[iPt]                      = new TH1D(fNameHistoGG.Data(),fNameHistoGG.Data(),fGammaGammaInvMassVSPtDummy->GetNbinsX(),
                                                                           0.,fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()));
        fHistoGconvGInvMassPtGConvBin[iPt]->Sumw2();

        Int_t startBin                                          = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
        Int_t endBin                                            = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

        fGammaGammaInvMassVSPtDummy->ProjectionX(fNameHistoGG.Data(),startBin,endBin);
        fHistoGconvGInvMassPtGConvBin[iPt]                      = (TH1D*)gDirectory->Get(fNameHistoGG.Data());
        if(fNRebin[iPt]>1){
            fHistoGconvGInvMassPtGConvBin[iPt]->Rebin(fNRebin[iPt]);
        }
    }
}

//****************************************************************************
//******* Function to compare two arrays *************************************
//****************************************************************************
Bool_t CompareArrays(Int_t nEntriesA, Double_t* arrayA, Int_t nEntriesB, Double_t* arrayB) {
    Bool_t returnValue = kTRUE;
    
    if (nEntriesA == nEntriesB) {
        for (Int_t i = 0; i < nEntriesA; i++) {
            if (arrayA[i] == arrayB[i]) {
                continue;
            } else {
                returnValue = kFALSE;
                break;
            }
        }
    } else {
        returnValue = kFALSE;
    }
    
    return returnValue;
}

//****************************************************************************
//******* Function to calculate pileup subtracted DCAz ***********************
//****************************************************************************
Bool_t CalculatePileUpSubtractedDCAz(TH1D* trueGamma, TH1D* trueSubGamma, TH1D* trueGammaX, TH1D* &trueSubGammaX) {
    Bool_t returnValue = kTRUE;
    
    TH1D* trueGammaClone = (TH1D*)trueGamma->Clone("trueGammaClone");
    trueGammaClone->Sumw2();
    
    TH1D* trueGammaXClone = (TH1D*)trueGammaX->Clone("trueGammaXClone");
    trueGammaXClone->Sumw2();
    
    for (Int_t i = 1; i <= trueGammaClone->GetNbinsX(); i++) {
        if (!trueGammaClone->GetBinContent(i)) {
            returnValue = kFALSE;
            
            trueGammaClone->SetBinContent(i, 1);
            trueGammaXClone->SetBinContent(i, 0);
        } else
            continue;
    }
    
    TH1D* ratioTrueXToTrue = (TH1D*)trueGammaXClone->Clone("ratioTrueXToTrue");
    ratioTrueXToTrue->Sumw2();
    ratioTrueXToTrue->Divide(ratioTrueXToTrue,trueGammaClone,1,1,"B");
    
    trueSubGammaX = (TH1D*)trueSubGamma->Clone(trueSubGammaX->GetName());
    trueSubGammaX->Sumw2();
    trueSubGammaX->Multiply(trueSubGammaX,ratioTrueXToTrue,1,1,"B");
    
    return returnValue;
}

//****************************************************************************
//******* Function to calculate ratios DCAz distributions ********************
//****************************************************************************
Bool_t CalculateDCAzDistributionRatio(TH1D*** inputNum, TH1D*** inputDenom, Int_t categoryFirst, Int_t categoryLast, TH1D* &ratio) {
    Bool_t returnValue              = kTRUE;
    
    TH1D* numerator                  = new TH1D("numerator", "numerator", fNBinsPtDummy, fBinsPtDummy);
    TH1D* denominator                = new TH1D("denominator", "denominator", fNBinsPtDummy, fBinsPtDummy);

    Double_t binContentNum, binContentDenom;
    Double_t binErrorTemp, binErrorNum, binErrorDenom;
    
    for (Int_t ptBin = 1; ptBin <= fNBinsPtDummy; ptBin++) {
        
        binContentNum               = 0;
        binContentDenom             = 0;
        
        binErrorNum                 = 0;
        binErrorDenom               = 0;
        
        for (Int_t cat = categoryFirst; cat <= categoryLast; cat++) {
            
            binContentNum           += inputNum[cat][ptBin]->IntegralAndError(-1000,1000,binErrorTemp);
            binErrorNum             += binErrorTemp*binErrorTemp;
            
            binContentDenom         += inputDenom[cat][ptBin]->IntegralAndError(-1000,1000,binErrorTemp);
            binErrorDenom           += binErrorTemp*binErrorTemp;
        }
        
        binErrorNum                 = TMath::Sqrt(binErrorNum);
        binErrorDenom               = TMath::Sqrt(binErrorDenom);
        
        numerator->SetBinContent(ptBin,         binContentNum);
        numerator->SetBinError(ptBin,           binErrorNum);
        
        if (binContentDenom) {
            denominator->SetBinContent(ptBin,   binContentDenom);
            denominator->SetBinError(ptBin,     binErrorDenom);
        } else {
            returnValue             = kFALSE;

            numerator->SetBinContent(ptBin,     0);
            denominator->SetBinContent(ptBin,   1);
            denominator->SetBinError(ptBin,     binErrorDenom);
        }
    }
    
    // if the histograms (numerator and denominator) are identical, the error is set to 0
    // this is a feature of the binomial-option of divide
    ratio->Divide(numerator,denominator,1,1,"B");

    return returnValue;
}

// overloading of function for different background estimation methods
Bool_t CalculateDCAzDistributionRatio(TH1D*** inputNum, TH1D**** inputDenom, Int_t backgroundExtractionMethod, Int_t categoryFirst, Int_t categoryLast, TH1D* &ratio) {
    Bool_t returnValue              = kTRUE;
    
    TH1D* numerator                  = new TH1D("numerator", "numerator", fNBinsPtDummy, fBinsPtDummy);
    TH1D* denominator                = new TH1D("denominator", "denominator", fNBinsPtDummy, fBinsPtDummy);
    
    Double_t binContentNum, binContentDenom;
    Double_t binErrorTemp, binErrorNum, binErrorDenom;
    
    for (Int_t ptBin = 1; ptBin <= fNBinsPtDummy; ptBin++) {
        
        binContentNum               = 0;
        binContentDenom             = 0;
        
        binErrorNum                 = 0;
        binErrorDenom               = 0;
        
        for (Int_t cat = categoryFirst; cat <= categoryLast; cat++) {
            
            binContentNum           += inputNum[cat][ptBin]->IntegralAndError(-1000,1000,binErrorTemp);
            binErrorNum             += binErrorTemp*binErrorTemp;
            
            binContentDenom         += inputDenom[cat][ptBin][backgroundExtractionMethod]->IntegralAndError(-1000,1000,binErrorTemp);
            binErrorDenom           += binErrorTemp*binErrorTemp;
        }
        
        binErrorNum                 = TMath::Sqrt(binErrorNum);
        binErrorDenom               = TMath::Sqrt(binErrorDenom);
        
        numerator->SetBinContent(ptBin,         binContentNum);
        numerator->SetBinError(ptBin,           binErrorNum);
        
        if (binContentDenom) {
            denominator->SetBinContent(ptBin,   binContentDenom);
            denominator->SetBinError(ptBin,     binErrorDenom);
        } else {
            returnValue             = kFALSE;
            
            numerator->SetBinContent(ptBin,     0);
            denominator->SetBinContent(ptBin,   1);
            denominator->SetBinError(ptBin,     binErrorDenom);
        }
    }
    
    // if the histograms (numerator and denominator) are identical, the error is set to 0
    // this is a feature of the binomial-option of divide
    ratio->Divide(numerator,denominator,1,1,"B");
    
    return returnValue;
}

//****************************************************************************
//******* Function to calculate pileup correction factors ********************
//****************************************************************************
Bool_t CalculatePileUpCorrectionFactor(TH1D* ratioWithWithoutPileUp, TH1D* &pileupCorrectionFactor, TF1* &fitToRatio) {
    
    TH1D* unityHisto                        = (TH1D*)pileupCorrectionFactor->Clone("unityHisto");
    unityHisto->Reset("ICES");
    unityHisto->Sumw2();
    for (Int_t i = 1; i < unityHisto->GetNbinsX()+1; i++) unityHisto->SetBinContent(i, 1);

    Bool_t returnValue                      = kTRUE;
    
    Int_t fFitStartBin, fitStatus;
    Double_t binContent, binError;
    Double_t stopX = 0;
    
    if (!CompareArrays(fNBinsPt, fBinsPt, fNBinsPtDummy, fBinsPtDummy)) {
        
        // binning between spectra and DCAz distributions differs, extract pileup correction factor from (partial) fit to the ratio
        Int_t iMax = (fNBinsPt >= fNBinsPtDummy) ? fNBinsPt : fNBinsPtDummy;
        for (Int_t i = 1; i < iMax+1; i++) {
            if ( (fBinsPt[i-1] == fBinsPtDummy[i-1]) && (fBinsPt[i] == fBinsPtDummy[i]) ) {
                
                pileupCorrectionFactor->SetBinContent(i,            1/ratioWithWithoutPileUp->GetBinContent(i));
                pileupCorrectionFactor->SetBinError(i,              ratioWithWithoutPileUp->GetBinError(i));
            } else {
                
                fFitStartBin = i;
                break;
            }
        }
        
        for (Int_t i = 1; i < ratioWithWithoutPileUp->GetNbinsX()+1; i++) {
            if (ratioWithWithoutPileUp->GetBinContent(i)) {
                continue;
            } else {
                stopX = ratioWithWithoutPileUp->GetXaxis()->GetBinLowEdge(i);
                break;
            }
        }
        
        if (stopX == 0) stopX               = pileupCorrectionFactor->GetXaxis()->GetBinUpEdge(pileupCorrectionFactor->GetNbinsX());
        
        // fit over whole range, otherwise sharp onset possible
        fitToRatio                          = new TF1("fitRatio", "1+[0]/TMath::Power((x-[1]), [2])", fBinsPtDummy[1], fBinsPtDummy[fNBinsPtDummy]);
        fitToRatio->SetParameters(1, 0, 1);
        fitToRatio->SetName(Form("%s_fit", ratioWithWithoutPileUp->GetName()));
        TFitResultPtr fitToRatioResult      = ratioWithWithoutPileUp->Fit(fitToRatio, "SIMNRE");
        fitStatus                           = fitToRatioResult;
        
        // fit status: https://root.cern.ch/doc/master/classTH1.html and https://root.cern.ch/doc/master/classTMinuit.html
        // fitStatus = migradResult + 10*minosResult + 100*hesseResult + 1000*improveResult.
        //  0: command executed normally
        //  1: command is blank, ignored
        //  2: command line unreadable, ignored
        //  3: unknown command, ignored
        //  4: abnormal termination (e.g., MIGRAD not converged)
        //  9: reserved
        //  10: END command
        //  11: EXIT or STOP command
        //  12: RETURN command
        
        if (fitStatus == 0 || fitStatus >= 1000) {
            // accepting fits, if everything went fine (i.e. fitstatus = 0) of if only improve had problems (i.e. fitstatus >= 1000)
            
            for (Int_t i = fFitStartBin; i < fNBinsPt+1; i++) {

                binContent                  = fitToRatio->Integral(         pileupCorrectionFactor->GetXaxis()->GetBinLowEdge(i),
                                                                            pileupCorrectionFactor->GetXaxis()->GetBinUpEdge(i),
                                                                            fitToRatioResult->GetParams()) / pileupCorrectionFactor->GetBinWidth(i);
                
                binError                    = fitToRatio->IntegralError(    pileupCorrectionFactor->GetXaxis()->GetBinLowEdge(i),
                                                                            pileupCorrectionFactor->GetXaxis()->GetBinUpEdge(i),
                                                                            fitToRatioResult->GetParams(),
                                                                            fitToRatioResult->GetCovarianceMatrix().GetMatrixArray()) / pileupCorrectionFactor->GetBinWidth(i);
                
                if (pileupCorrectionFactor->GetXaxis()->GetBinUpEdge(i) <= stopX) {
                    pileupCorrectionFactor->SetBinContent(i,        1/binContent);
                    pileupCorrectionFactor->SetBinError(i,          binError/binContent/binContent);
                } else {
                    pileupCorrectionFactor->SetBinContent(i,        0);
                    pileupCorrectionFactor->SetBinError(i,          0);
                }
            }
        } else {
            // if fit failed, set correction factors to one
            
            returnValue                     = kFALSE;
            cout << "WARNING: fit to " << ratioWithWithoutPileUp->GetName() << " failed, correction factor set to one!" << endl;

            for (Int_t i = fFitStartBin; i < fNBinsPt+1; i++) {
                pileupCorrectionFactor->SetBinContent(i,            1);
                pileupCorrectionFactor->SetBinError(i,              0);     // not good, how to treat the error?
            }
        }
    } else {
        // same binning between spectra and DCAz distributions, pileup correction factor is inverse of ratio
        
        for (Int_t i = 1; i < fNBinsPt+1; i++) {
            if (!ratioWithWithoutPileUp->GetBinContent(i))
                ratioWithWithoutPileUp->SetBinContent(i,            1);
        }

        pileupCorrectionFactor->Divide(unityHisto,ratioWithWithoutPileUp,1,1,"B");
        
        // set fit to NULL
        fitToRatio                          = NULL;
    }
    
    delete unityHisto;
    
    return returnValue;
}