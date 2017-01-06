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
                                Int_t mode                  = 0,
                                Int_t nPileupMethod         = 0,
                                TString purityFileName      = ""
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
    fDate                                                                       = ReturnDateString();
    fDirectPhoton                                                               = directphotonPlots;
    fEnergyFlag                                                                 = option;
    fPrefix                                                                     = meson;
    fPeriodFlag                                                                 = period;
    fSuffix                                                                     = suffix;
    fMeson                                                                      = meson;
    fMode                                                                       = mode;
    nPileupMethodUsed                                                           = nPileupMethod;
    cout << "Pictures are saved as " << suffix.Data() << endl;
    
    
    //************************************ Separate cutstrings ***********************************
    fCutSelection                                                               = cutSelection;
    fCutSelectionRead                                                           = cutSelection;
    TString dummy                                                               = "" ;
    ReturnSeparatedCutNumberAdvanced( fCutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, dummy, fMesonCutSelection, fMode);
    
    fEventCutSelectionRead                                                      = fEventCutSelection.Data();
    fGammaCutSelectionRead                                                      = fGammaCutSelection.Data();
    fMesonCutSelectionRead                                                      = fMesonCutSelection.Data();
    if (addSig) {
        cout << "running added Signal" << endl;
        cout << fEventCutSelection.Data() << endl;
        fEventCutSelection.Replace(GetEventRejectExtraSignalsCutPosition(),1,"2");
        cout << fEventCutSelection.Data() << endl;
        fEventCutSelectionRead                                                  = fEventCutSelection;
        fGammaCutSelectionRead                                                  = fGammaCutSelection;
        fMesonCutSelectionRead                                                  = fMesonCutSelection;
        if (fMode==9)       fCutSelectionRead                                   = Form("%s%s_%s", fEventCutSelection.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
        else if (fMode==0)  fCutSelectionRead                                   = Form("%s_%s_%s",fEventCutSelection.Data(), fGammaCutSelection.Data(), fMesonCutSelection.Data());
        cout << fCutSelectionRead.Data() << endl;
    }
    // set global variables for rap and BG number
    TString rapidityRange;
    fYMaxMeson                                                                  = ReturnRapidityStringAndDouble(fMesonCutSelectionRead, rapidityRange);
    fBackgroundMultNumber                                                       = ReturnBackgroundMult(fMesonCutSelection);
    //***************************** Check for data driven purity ****************************************
    TString namePurityHistogram                                                 = "";
    TFile*  purityFile                                                          = NULL;
    if (purityFileName.CompareTo("")) {
        
        // intercept dEdx cut
        TString dEdxCut( fGammaCutSelection(9,1) );
        if (dEdxCut.Atoi() != 0) {
            cout << "Data driven purity with kappa cut requested but usual dEdx cut given, returning." << endl;
            return;
        }

        // load purity file
        cout << "loading data driven purity file: " << purityFileName << endl;
        purityFile                                                              = new TFile(purityFileName.Data());
        
        // check for kappa cut and load corresponding purity histo
        TString kappaCut( fGammaCutSelection(8,1) );
        if (kappaCut.Atoi() == 3) {
            namePurityHistogram                                                 = "hSignalPurity1";
        } else if (kappaCut.Atoi() == 5) {
            namePurityHistogram                                                 = "hSignalPurity2";
        } else if (kappaCut.Atoi() == 6) {
            namePurityHistogram                                                 = "hSignalPurity3";
        } else {
            cout << "Kappa cut not recognized, not implemented yet." << endl;
            return;
        }
        fHistoPurityKappaTemplates                                              = (TH1D*)purityFile->Get(namePurityHistogram.Data());   // replaces "GammaTruePurity_Pt" in MC output, will be used
        if (fHistoPurityKappaTemplates) fUseDataDrivenPurity                     = kTRUE;                                               // for bin-by-bin correction in CorrectGammaV2
    }
    
    //***************************** Load binning for spectrum *******************************************
    Initialize(fMeson, fEnergyFlag, numberOfBins, fMode, addSig);

    //****************************** Set specific histogram names****************************************
    ObjectNameDCGammaConvRPt                                                    = "ESD_TrueDoubleCountConvGamma_R_Pt";
    ObjectNameGammaConvMultipleCount                                            = "ESD_TrueMultipleCountConvGamma";

    //****************************** Specification of collision system **********************************
    TString textProcess = ReturnMesonString (fPrefix);
    if(textProcess.CompareTo("") == 0 ){
        cout << "Meson unknown" << endl;
        return ;
    }
    
    fTextMeasurement                                                            = Form("%s #rightarrow #gamma#gamma", textProcess.Data());
    fCollisionSystem                                                            = ReturnFullCollisionsSystem(fEnergyFlag);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    fDetectionProcess                                                           = ReturnFullTextReconstructionProcess(fMode);
    TString fDecayChannel                                                       = "#gamma#gamma";
    //***************************** Specification Data/MC ***********************************************
    if(isMC.CompareTo("kTRUE") == 0){
        fIsMC                                                                   = 1;
        fPrefix2                                                                = "MC";
    } else {
        fIsMC                                                                   = 0;
        fPrefix2                                                                = "data";
    }
    const char* FileDataLogname         = Form("%s/%s/%s_%s_GammaExtractionEffiCheck_RAWDATA%s_%s.dat", cutSelection.Data(), fEnergyFlag.Data(), fPrefix.Data(), fPrefix2.Data(), fPeriodFlag.Data(),
                                        fCutSelectionRead.Data());
    fFileDataLog.open(FileDataLogname, ios::out);


    //************************************** Read file ***************************************************
    TFile* f                                                                    = new TFile(file.Data());
    TString autoDetectedMainDir                                                 = AutoDetectMainTList(mode , f);
    if (autoDetectedMainDir.CompareTo("") == 0){
        cout << "ERROR: trying to read file, which is incompatible with mode selected" << endl;;
        return;
    }

    TList *TopDir =(TList*)f->Get(autoDetectedMainDir.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }
    TList* HistosGammaConversion                                                = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if (!HistosGammaConversion){
      cout << "ERROR: folder with Cutnumber - " <<   fCutSelectionRead.Data() << "not contained in file " << file.Data() << endl;
      return;
    }    
    TList* ESDContainer                                                         = (TList*)HistosGammaConversion->FindObject(Form("%s ESD histograms",fCutSelectionRead.Data()));
    if (fEnablePCM){
        TList* ConvCutsContainer                                                = (TList*)HistosGammaConversion->FindObject(Form("ConvCuts_%s",fGammaCutSelectionRead.Data()));
        fHistoPhotonIsSelected                                                  = (TH1D*)ConvCutsContainer->FindObject(Form("IsPhotonSelected %s",fGammaCutSelectionRead.Data()));
        fHistoPhotonIsSelected->Scale(1./fHistoPhotonIsSelected->GetEntries());
    }
    
    // general histograms
    fNumberOfGoodESDTracks                                                      = (TH1D*)ESDContainer->FindObject("GoodESDTracks");
    fEventQuality                                                               = (TH1D*)ESDContainer->FindObject("NEvents");
    if (option.CompareTo("PbPb_2.76TeV") == 0)  fNEvents                        = fEventQuality->GetBinContent(1);
    else                                        fNEvents                        = GetNEvents(fEventQuality);
    
    fMesonMassExpect                            = TDatabasePDG::Instance()->GetParticle(fMesonId)->Mass();
    // calculate number of events for normalization
    if (fEnergyFlag.CompareTo("PbPb_2.76TeV") == 0 || fEnergyFlag.CompareTo("pPb_5.023TeV") == 0){
        fNEvents        = fEventQuality->GetBinContent(1);
    } else {
        fNEvents        =  GetNEvents(fEventQuality);
    }
    cout<< "The mass of the meson is: "<< fMesonMassExpect<< " Events analysed: "<< fNEvents<< endl;
    
    // *******************************************************************************************************
    // ******************************* Read histograms from data for PCM *************************************
    // *******************************************************************************************************
    if (fEnablePCM){
        // read reconstructed gamma histograms
        fHistoGammaConvPt                                                       = (TH1D*)ESDContainer->FindObject("ESD_ConvGamma_Pt");
        fHistoGammaConvPt->Sumw2();
        fHistoGammaConvPtOrBin                                                  = (TH1D*)fHistoGammaConvPt->Clone("ESD_ConvGamma_Pt_OriginalBinning");
        RebinSpectrum(fHistoGammaConvPt,"");
        
        // read dca tree for conversions
        if(!addSig && !(mode == 4 || mode == 5)){
            DCAContainer                                                        = (TList*)HistosGammaConversion->FindObject(Form("%s Photon DCA tree",fCutSelectionRead.Data()));
            if(DCAContainer){
                dcaTree                                                         = (TTree*)DCAContainer->FindObject("ESD_ConvGamma_Pt_Dcaz_R_Eta");
                FillDCAHistogramsFromTree(dcaTree,kFALSE);
                CalculatePileUpBackground(kFALSE);
                pileUpCorrection                                                = kTRUE;
            }         
        }
    }
    
    // *******************************************************************************************************
    // ******************************* Read histograms from data for Calo ************************************
    // *******************************************************************************************************
    if (fEnableCalo){
        if (mode == 2 || mode == 3){
            CaloContainer                                                       = (TList*)HistosGammaConversion->FindObject(Form("%s Cluster Output",fCutSelectionRead.Data()));
            fHistoGammaCaloPt                                                   = (TH1D*)CaloContainer->FindObject("ClusGamma_Pt");
            fHistoGammaCaloPt->Sumw2();
            fHistoGammaCaloPtOrBin                                              = (TH1D*)fHistoGammaCaloPt->Clone("ClusGamma_Pt_OriginalBinning");
            RebinSpectrum(fHistoGammaCaloPt,"");
        } else if ( mode == 4 || mode == 5){
            fHistoGammaCaloPt                                                   = (TH1D*)ESDContainer->FindObject("ClusGamma_Pt");
            fHistoGammaCaloPt->Sumw2();
            fHistoGammaCaloPtOrBin                                              = (TH1D*)fHistoGammaCaloPt->Clone("ClusGamma_Pt_OriginalBinning");
            RebinSpectrum(fHistoGammaCaloPt,"");
        }
    }
    
    // *******************************************************************************************************
    // ********************* Read Meson hists data for pi0 extraction for PCM-Calo mode***********************
    // *******************************************************************************************************
    if ( mode == 2 || mode == 3){
    
        TString ObjectNameESD                                                   = "ESD_Mother_InvMass_PtConv";
        TString ObjectNameBck                                                   = "ESD_Background_InvMass_PtConv";
        
        fGammaGammaInvMassVSPt                                                  = (TH2D*)ESDContainer->FindObject(ObjectNameESD.Data());
        fGammaGammaInvMassVSPt->Sumw2();
        fBckInvMassVSPt                                                         = (TH2D*)ESDContainer->FindObject(ObjectNameBck.Data());
        fBckInvMassVSPt->Sumw2();

        FillMassHistosArray(fGammaGammaInvMassVSPt);
        ProduceBckProperWeighting(fGammaGammaInvMassVSPt, fBckInvMassVSPt);
    }
    
    // *******************************************************************************************************
    // ************************************* Load cocktail input (if available) ******************************
    // *******************************************************************************************************
    fUseCocktail                                                                = LoadSecondariesFromCocktailFile(cutSelection, option);
    
    // *******************************************************************************************************
    // ******************************* Read MC quantities ****************************************************
    // *******************************************************************************************************
    if(fIsMC){
        // copy reconstructed photon histo for MC
        if (fEnablePCM){
            fHistoGammaMCrecConvPt                                              = (TH1D*)fHistoGammaConvPt->Clone("MCrec_ConvGamma_Pt");
            fHistoGammaMCrecConvPtOrBin                                         = (TH1D*)fHistoGammaConvPtOrBin->Clone("MCrec_ConvGamma_Pt_OriginalBinning");
        }
        if (fEnableCalo){
            fHistoGammaMCrecCaloPt                                              = (TH1D*)fHistoGammaCaloPt->Clone("MCrec_CaloGamma_Pt");
            fHistoGammaMCrecCaloPtOrBin                                         = (TH1D*)fHistoGammaCaloPtOrBin->Clone("MCrec_CaloGamma_OriginalBinning_Pt");
        }
        
        // MC container (contains all input MC histograms)    
        TList *MCContainer                                                      = (TList*)HistosGammaConversion->FindObject(Form("%s MC histograms",fCutSelectionRead.Data()));
        
        // reading input distributions
        if (fEnablePCM){
            fHistoGammaMCConvPt                                                 = (TH1D*)MCContainer->FindObject("MC_ConvGamma_Pt");
            fHistoGammaMCConvPt->Sumw2();
            fHistoGammaMCConvPtOrBin                                            = (TH1D*)fHistoGammaMCConvPt->Clone("MC_ConvGamma_Pt_OriginalBinning");
            fHistoGammaMCConvPtOrBin->Scale(1./fHistoGammaMCConvPtOrBin->GetBinWidth(5));
            RebinSpectrum(fHistoGammaMCConvPt);
            
            // secondary conv gammas
            if(fUseCocktail){
                f2DHistoSecondaryGammaMCConvPt                                  = (TH2D*)MCContainer->FindObject("MC_SecondaryConvGamma_Pt");
                if(f2DHistoSecondaryGammaMCConvPt){
                    fHistoSecondaryGammaConvFromXFromK0sPt                      = (TH1D*)f2DHistoSecondaryGammaMCConvPt->ProjectionX("ESD_SecondaryConvGammaFromXFromK0s_Pt",1,1,"e");
                    fHistoSecondaryGammaConvFromXFromK0sPt->Sumw2();
                    fHistoSecondaryGammaConvFromXFromK0sPtOrBin                 = (TH1D*)fHistoSecondaryGammaConvFromXFromK0sPt->Clone("fHistoSecondaryGammaConvFromXFromK0sPtOrBin");
                    RebinSpectrum(fHistoSecondaryGammaConvFromXFromK0sPt);
                    fHistoSecondaryGammaConvFromXFromK0lPt                      = (TH1D*)f2DHistoSecondaryGammaMCConvPt->ProjectionX("ESD_SecondaryConvGammaFromXFromK0l_Pt",2,2,"e");
                    fHistoSecondaryGammaConvFromXFromK0lPt->Sumw2();
                    fHistoSecondaryGammaConvFromXFromK0lPtOrBin                 = (TH1D*)fHistoSecondaryGammaConvFromXFromK0lPt->Clone("fHistoSecondaryGammaConvFromXFromK0lPtOrBin");
                    RebinSpectrum(fHistoSecondaryGammaConvFromXFromK0lPt);
                    fHistoSecondaryGammaConvFromXFromLambdaPt                   = (TH1D*)f2DHistoSecondaryGammaMCConvPt->ProjectionX("ESD_SecondaryConvGammaFromXFromLambda_Pt",3,3,"e");
                    fHistoSecondaryGammaConvFromXFromLambdaPt->Sumw2();
                    fHistoSecondaryGammaConvFromXFromLambdaPtOrBin              = (TH1D*)fHistoSecondaryGammaConvFromXFromLambdaPt->Clone("fHistoSecondaryGammaConvFromXFromLambdaPtOrBin");
                    RebinSpectrum(fHistoSecondaryGammaConvFromXFromLambdaPt);
                }
            }
        }
        if (fEnableCalo && mode == 2 ){
            fHistoGammaMCAllInEMCAccPt                                          = (TH1D*)MCContainer->FindObject("MC_AllGammaEMCALAcc_Pt");
            fHistoGammaMCAllInEMCAccPt->Sumw2();
            fHistoGammaMCAllInEMCAccPtOrBin                                     = (TH1D*)fHistoGammaMCAllInEMCAccPt->Clone("MC_AllGammaEMCALAcc_OriginalBinning_MCPt");
            fHistoGammaMCAllInEMCAccPtOrBin->Scale(1./fHistoGammaMCAllInEMCAccPtOrBin->GetBinWidth(5));
            fHistoGammaMCAllInEMCAccPtOrBin->GetXaxis()->SetRangeUser(0.,25.);
            RebinSpectrum(fHistoGammaMCAllInEMCAccPt);
        }    
        
        // all MC gammas
        fHistoGammaMCAllPt                                                      = (TH1D*)MCContainer->FindObject("MC_AllGamma_Pt");
        fHistoGammaMCAllPt->Sumw2();
        fHistoGammaMCAllPtOrBin                                                 = (TH1D*)fHistoGammaMCAllPt->Clone("MC_AllGamma_OriginalBinning_MCPt");
        fHistoGammaMCAllPtOrBin->Scale(1./fHistoGammaMCAllPtOrBin->GetBinWidth(5));
//         fHistoGammaMCAllPtOrBin->GetXaxis()->SetRangeUser(0.,25.);
        RebinSpectrum(fHistoGammaMCAllPt);
        
        // all secondary MC gammas
        if(fUseCocktail){
            f2DHistoAllSecondaryGammaMCPt                                       = (TH2D*)MCContainer->FindObject("MC_AllSecondaryGamma_Pt");
            if(f2DHistoAllSecondaryGammaMCPt){
                fHistoAllSecondaryGammaFromXFromK0sPt                           = (TH1D*)f2DHistoAllSecondaryGammaMCPt->ProjectionX("ESD_AllSecondaryGammaFromXFromK0s_Pt",1,1,"e");
                fHistoAllSecondaryGammaFromXFromK0sPt->Sumw2();
                fHistoAllSecondaryGammaFromXFromK0sPtOrBin                      = (TH1D*)fHistoAllSecondaryGammaFromXFromK0sPt->Clone("ESD_AllSecondaryGammaFromXFromK0s_PtOrBin");
                RebinSpectrum(fHistoAllSecondaryGammaFromXFromK0sPt);
                fHistoAllSecondaryGammaFromXFromK0lPt                           = (TH1D*)f2DHistoAllSecondaryGammaMCPt->ProjectionX("ESD_AllSecondaryGammaFromXFromK0l_Pt",2,2,"e");
                fHistoAllSecondaryGammaFromXFromK0lPt->Sumw2();
                fHistoAllSecondaryGammaFromXFromK0lPtOrBin                      = (TH1D*)fHistoAllSecondaryGammaFromXFromK0lPt->Clone("ESD_AllSecondaryGammaFromXFromK0l_PtOrBin");
                RebinSpectrum(fHistoAllSecondaryGammaFromXFromK0lPt);
                fHistoAllSecondaryGammaFromXFromLambdaPt                        = (TH1D*)f2DHistoAllSecondaryGammaMCPt->ProjectionX("ESD_AllSecondaryGammaFromXFromLambda_Pt",3,3,"e");
                fHistoAllSecondaryGammaFromXFromLambdaPt->Sumw2();
                fHistoAllSecondaryGammaFromXFromLambdaPtOrBin                   = (TH1D*)fHistoAllSecondaryGammaFromXFromLambdaPt->Clone("ESD_AllSecondaryGammaFromXFromLambda_PtOrBin");
                RebinSpectrum(fHistoAllSecondaryGammaFromXFromLambdaPt);
            }
        }
        
        // Gamma from Decay
        fHistoGammaMCDecayPt                                                    = new TH1D*[7];
        for(Int_t i = 0; i<7; i++){
            fHistoGammaMCDecayPt[i]                                             = (TH1D*)MCContainer->FindObject(Form("MC_DecayGamma%s_Pt",fDecays[i].Data()));
            fHistoGammaMCDecayPt[i]->Sumw2();
            fHistoGammaMCDecayPt[i]->Scale(1./fHistoGammaMCDecayPt[i]->GetBinWidth(5));
        }
        
        // container with histos for validated reconstructed photons
        TList *TrueConversionContainer                                          = (TList*)HistosGammaConversion->FindObject(Form("%s True histograms",fCutSelectionRead.Data()));
        
        // reading reconstructed validated distributions
        if (fEnablePCM){
            fHistoGammaTrueConvPt                                               = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueConvGamma_Pt");
            fHistoGammaTrueConvPt->Sumw2();
            fHistoGammaTrueConvPtOrBin                                          = (TH1D*)fHistoGammaTrueConvPt->Clone("TrueConvGamma_Pt_OriginalBinning");
            fHistoGammaTrueConvPtOrBin->Sumw2();
            RebinSpectrum(fHistoGammaTrueConvPt);
    
            fHistoGammaTruePrimaryConvPt                                        = (TH1D*)TrueConversionContainer->FindObject("ESD_TruePrimaryConvGamma_Pt");
            fHistoGammaTruePrimaryConvPt->Sumw2();
            fHistoGammaTruePrimaryConvPtOrBin                                   = (TH1D*)fHistoGammaTruePrimaryConvPt->Clone("TruePrimaryConvGamma_Pt_OriginalBinning");
            RebinSpectrum(fHistoGammaTruePrimaryConvPt);
            
            fHistoTrueGammaConvDCRVSPt                                          = (TH2D*)TrueConversionContainer->FindObject(ObjectNameDCGammaConvRPt.Data());
            if (fHistoTrueGammaConvDCRVSPt != NULL) fEnableDCConv               = kTRUE;
            if (fEnableDCConv){
                fHistoTrueGammaConvDCRVSPt->Sumw2();
                fHistoTrueGammaConvDCPt                                         = (TH1D*)fHistoTrueGammaConvDCRVSPt->ProjectionY("fHistoTrueGammaConvDCPt");
                fHistoTrueGammaConvDCR                                          = (TH1D*)fHistoTrueGammaConvDCRVSPt->ProjectionX("fHistoTrueGammaConvDCR");
                fHistoTrueGammaConvMultipleCount                                = (TH1F*)TrueConversionContainer->FindObject(ObjectNameGammaConvMultipleCount.Data());
            }
            
            // check for new trainoutput, secondary historgrams in 2D
            fHistogramDimension                                                 = TrueConversionContainer->FindObject("ESD_TrueSecondaryConvGamma_Pt")->ClassName();
            if (fHistogramDimension.Contains("2")){
               nHistogramDimension                                              = 2;
               cout << "Found 2D histgrams -> using new trainoutput" << endl;
            }
            
            if (nHistogramDimension==2) {
                
                // true secondary conv gamma rec Pt
                f2DHistoGammaTrueSecondaryConvPt                                = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryConvGamma_Pt");
                fHistoGammaTrueSecondaryConvPt                                  = (TH1D*)f2DHistoGammaTrueSecondaryConvPt->ProjectionX("ESD_TrueSecondaryConvGamma_Pt",1,4,"e");
                fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt                 = (TH1D*)f2DHistoGammaTrueSecondaryConvPt->ProjectionX("ESD_TrueSecondaryConvGammaFromXFromK0s_Pt",1,1,"e");
                fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt->Sumw2();
                fHistoGammaTrueSecondaryConvGammaFromXFromK0lPt                 = (TH1D*)f2DHistoGammaTrueSecondaryConvPt->ProjectionX("ESD_TrueSecondaryConvGammaFromXFromK0l_Pt",2,2,"e");
                fHistoGammaTrueSecondaryConvGammaFromXFromK0lPt->Sumw2();
                fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt              = (TH1D*)f2DHistoGammaTrueSecondaryConvPt->ProjectionX("ESD_TrueSecondaryConvGammaFromXFromLambda_Pt",3,3,"e");
                fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt->Sumw2();
                fHistoGammaTrueSecondaryConvGammaRestPt                         = (TH1D*)f2DHistoGammaTrueSecondaryConvPt->ProjectionX("ESD_TrueSecondaryConvGammaRest_Pt",4,4,"e");
                fHistoGammaTrueSecondaryConvGammaRestPt->Sumw2();
                
                // MC histograms for calculation of raw secondary spectra from cocktail
                if(fUseCocktail){
                    
                    // true secondary conv gamma in MC Pt
                    f2DHistoGammaTrueSecondaryConvMCPt                          = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryConvGamma_MCPt");
                    fHistoGammaTrueSecondaryConvGammaFromXFromK0sMCPt           = (TH1D*)f2DHistoGammaTrueSecondaryConvMCPt->ProjectionX("ESD_TrueSecondaryConvGammaFromXFromK0s_MCPt",1,1,"e");
                    fHistoGammaTrueSecondaryConvGammaFromXFromK0sMCPt->Sumw2();
                    fHistoGammaTrueSecondaryConvGammaFromXFromK0sMCPtOrBin      = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromK0sMCPt->Clone("ESD_TrueSecondaryConvGammaFromXFromK0s_MCPtOrBin");
                    RebinSpectrum(fHistoGammaTrueSecondaryConvGammaFromXFromK0sMCPt);
                    fHistoGammaTrueSecondaryConvGammaFromXFromK0lMCPt           = (TH1D*)f2DHistoGammaTrueSecondaryConvMCPt->ProjectionX("ESD_TrueSecondaryConvGammaFromXFromK0l_MCPt",2,2,"e");
                    fHistoGammaTrueSecondaryConvGammaFromXFromK0lMCPt->Sumw2();
                    fHistoGammaTrueSecondaryConvGammaFromXFromK0lMCPtOrBin      = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromK0lMCPt->Clone("ESD_TrueSecondaryConvGammaFromXFromK0l_MCPtOrBin");
                    RebinSpectrum(fHistoGammaTrueSecondaryConvGammaFromXFromK0lMCPt);
                    fHistoGammaTrueSecondaryConvGammaFromXFromLambdaMCPt        = (TH1D*)f2DHistoGammaTrueSecondaryConvMCPt->ProjectionX("ESD_TrueSecondaryConvGammaFromXFromLambda_MCPt",3,3,"e");
                    fHistoGammaTrueSecondaryConvGammaFromXFromLambdaMCPt->Sumw2();
                    fHistoGammaTrueSecondaryConvGammaFromXFromLambdaMCPtOrBin   = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromLambdaMCPt->Clone("ESD_TrueSecondaryConvGammaFromXFromLambda_MCPtOrBin");
                    RebinSpectrum(fHistoGammaTrueSecondaryConvGammaFromXFromLambdaMCPt);
                
                    // secondary response matrices
                    fHistoGammaTrueSecondaryFromXFromK0sConv_MCPt_recPt_MC      = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryConvGammaFromXFromK0sESD_MCPtPt");
                    fHistoGammaTrueSecondaryFromXFromK0sConv_MCPt_recPt_MC->Sumw2();
                    fHistoGammaTrueSecondaryFromXFromK0lConv_MCPt_recPt_MC      = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryConvGammaFromXFromK0lESD_MCPtPt");
                    fHistoGammaTrueSecondaryFromXFromK0lConv_MCPt_recPt_MC->Sumw2();
                    fHistoGammaTrueSecondaryFromXFromLambdaConv_MCPt_recPt_MC   = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryConvGammaFromXFromLambdaESD_MCPtPt");
                    fHistoGammaTrueSecondaryFromXFromLambdaConv_MCPt_recPt_MC->Sumw2();
                }
            } else {
                // old train output, secondary histos not 2D
                fHistoGammaTrueSecondaryConvPt                                  = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryConvGamma_Pt");
                fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt                 = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryConvGammaFromXFromK0s_Pt");
                fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt              = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryConvGammaFromXFromLambda_Pt");
            }
            
            // rebin true secondary histos
            fHistoGammaTrueSecondaryConvPt->Sumw2();
            fHistoGammaTrueSecondaryConvPtOrBin                                 = (TH1D*)fHistoGammaTrueSecondaryConvPt->Clone("ESD_TrueSecondaryConvGamma_Pt_OriginalBinning");
            RebinSpectrum(fHistoGammaTrueSecondaryConvPt);    
    
            fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt->Sumw2();
            fHistoGammaTrueSecondaryConvGammaFromXFromK0sPtOrBin                = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt->Clone("ESD_TrueSecondaryConvGammaFromXFromK0s_Pt_OriginalBinning");
            RebinSpectrum(fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt);
            
            if(fHistoGammaTrueSecondaryConvGammaFromXFromK0lPt){
                fHistoGammaTrueSecondaryConvGammaFromXFromK0lPt->Sumw2();
                fHistoGammaTrueSecondaryConvGammaFromXFromK0lPtOrBin            = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromK0lPt->Clone("ESD_TrueSecondaryConvGammaFromXFromK0s_Pt_OriginalBinning");
                RebinSpectrum(fHistoGammaTrueSecondaryConvGammaFromXFromK0lPt);
            }

            fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt->Sumw2();
            fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPtOrBin             = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt->Clone("ESD_TrueSecondaryConvGammaFromXFromLambda_Pt_OriginalBinning");
            RebinSpectrum(fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt);

            if(fHistoGammaTrueSecondaryConvGammaRestPt){
                fHistoGammaTrueSecondaryConvGammaRestPt->Sumw2();
                fHistoGammaTrueSecondaryConvGammaRestPtOrBin                    = (TH1D*)fHistoGammaTrueSecondaryConvGammaRestPt->Clone("ESD_TrueSecondaryConvGammaFromRest_Pt_OriginalBinning");
                RebinSpectrum(fHistoGammaTrueSecondaryConvGammaRestPt);
            }

            // combinatorial BG distributions
            fHistoCombinatorialBackground                                       = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueCombinatorial_Pt");
            fHistoCombinatorialSpecies                                          = new TH1D*[17];
            fHistoCombinatorialSpecies[nCombinatorics]                          = (TH1D*)fHistoCombinatorialBackground->ProjectionX(Form("ESD_TrueComb%s_Pt",fCombinatorics[nCombinatorics].Data()));
            fHistoCombinatorialSpecies[nCombinatorics]->Sumw2();
            RebinSpectrum(fHistoCombinatorialSpecies[nCombinatorics]);
            for(Int_t i = 0; i<nCombinatorics; i++){
                fHistoCombinatorialSpecies[i]                                   = (TH1D*)fHistoCombinatorialBackground->ProjectionX(Form("ESD_TrueComb%s_Pt",fCombinatorics[i].Data()),i+1,i+1);
                fHistoCombinatorialSpecies[i]->Sumw2();
                RebinSpectrum(fHistoCombinatorialSpecies[i]);
        
            }

            // primary response matrix PCM
            fHistoGammaTruePrimaryConv_recPt_MCPt_MC                            = (TH2D*)TrueConversionContainer->FindObject("ESD_TruePrimaryConvGammaESD_PtMCPt");
            fHistoGammaTruePrimaryConv_recPt_MCPt_MC->Sumw2();
            fHistoGammaTruePrimaryConvMCPt                                      = (TH1D*)fHistoGammaTruePrimaryConv_recPt_MCPt_MC->ProjectionY("ESD_TruePrimaryConvGamma_MCPt");
            fHistoGammaTruePrimaryConvMCPt->Sumw2();
            fHistoGammaTruePrimaryConvMCPtOrBin                                 = (TH1D*)fHistoGammaTruePrimaryConvMCPt->Clone("ESD_TruePrimaryConvGamma_MCPt_OriginalBinning");
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

                fHistoGammaTrueSecondaryCaloFromXFromK0sPt           = (TH1D*)CaloContainer->FindObject("TrueSecondaryClusGammaFromK0s_Pt");
                fHistoGammaTrueSecondaryCaloFromXFromK0sPt->Sumw2();
                fHistoGammaTrueSecondaryCaloFromXFromK0sPtOrBin      = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromK0sPt->Clone("TrueSecondaryClusGammaFromK0s_Pt_OriginalBinning");
                RebinSpectrum(fHistoGammaTrueSecondaryCaloFromXFromK0sPt);

                fHistoGammaTrueSecondaryCaloFromXFromLambdaPt        = (TH1D*)CaloContainer->FindObject("TrueSecondaryClusGammaFromLambda_Pt");
                fHistoGammaTrueSecondaryCaloFromXFromLambdaPt->Sumw2();
                fHistoGammaTrueSecondaryCaloFromXFromLambdaPtOrBin   = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromLambdaPt->Clone("TrueSecondaryClusGammaFromLambda_Pt_OriginalBinning");
                RebinSpectrum(fHistoGammaTrueSecondaryCaloFromXFromLambdaPt);
                
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
                
                // MC true cluster gammas
                fHistoGammaTrueCaloUnConvPt                                                 = (TH1D*)TrueConversionContainer->FindObject("TrueClusUnConvGamma_Pt");
                fHistoGammaTrueCaloUnConvPt->SetName("TrueClusUnConvGamma_Pt");
                fHistoGammaTrueCaloUnConvPt->Sumw2();
                fHistoGammaTrueCaloUnConvPtOrBin                                            = (TH1D*)fHistoGammaTrueCaloUnConvPt->Clone("TrueClusUnConvGamma_Pt_OriginalBinning");
                fHistoGammaTrueCaloUnConvPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTrueCaloUnConvPt);

                fHistoGammaTrueCaloConvPt                                                   = (TH1D*)TrueConversionContainer->FindObject("TrueClusConvGamma_Pt");
                fHistoGammaTrueCaloConvPt->SetName("TrueClusConvGamma_Pt");
                fHistoGammaTrueCaloConvPt->Sumw2();
                fHistoGammaTrueCaloConvPtOrBin                                              = (TH1D*)fHistoGammaTrueCaloConvPt->Clone("TrueClusConvGamma_Pt_OriginalBinning");
                fHistoGammaTrueCaloConvPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTrueCaloConvPt);
                
                fHistoGammaTrueCaloPt                                                       = (TH1D*)TrueConversionContainer->FindObject("TrueClusGamma_Pt");
                fHistoGammaTrueCaloPt->Sumw2();
                fHistoGammaTrueCaloPtOrBin                                                  = (TH1D*)fHistoGammaTrueCaloPt->Clone("TrueClusGamma_Pt_OriginalBinning");
                fHistoGammaTrueCaloPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTrueCaloPt);
                
                // MC true primary cluster gammas
                fHistoGammaTruePrimaryCaloUnConvPt                                          = (TH1D*)TrueConversionContainer->FindObject("TruePrimaryClusGamma_Pt");
                fHistoGammaTruePrimaryCaloUnConvPt->SetName("TruePrimaryClusUnConvGamma_Pt");
                fHistoGammaTruePrimaryCaloUnConvPt->Sumw2();
                fHistoGammaTruePrimaryCaloUnConvPtOrBin                                     = (TH1D*)fHistoGammaTruePrimaryCaloUnConvPt->Clone("TruePrimaryClusUnConvGamma_Pt_OriginalBinning");
                fHistoGammaTruePrimaryCaloUnConvPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTruePrimaryCaloUnConvPt);

                fHistoGammaTruePrimaryCaloConvPt                                            = (TH1D*)TrueConversionContainer->FindObject("TruePrimaryClusConvGamma_Pt");
                fHistoGammaTruePrimaryCaloConvPt->Sumw2();
                fHistoGammaTruePrimaryCaloConvPtOrBin                                       = (TH1D*)fHistoGammaTruePrimaryCaloConvPt->Clone("TruePrimaryClusConvGamma_Pt_OriginalBinning");
                fHistoGammaTruePrimaryCaloConvPtOrBin->Sumw2();
                RebinSpectrum(fHistoGammaTruePrimaryCaloConvPt);

                fHistoGammaTruePrimaryCaloPt                                                = (TH1D*)fHistoGammaTruePrimaryCaloUnConvPt->Clone("TruePrimaryClusGamma_Pt");
                fHistoGammaTruePrimaryCaloPt->Add(fHistoGammaTruePrimaryCaloConvPt);
                fHistoGammaTruePrimaryCaloPtOrBin                                           = (TH1D*)fHistoGammaTruePrimaryCaloUnConvPtOrBin->Clone("TruePrimaryClusGamma_Pt_OriginalBinning");
                fHistoGammaTruePrimaryCaloPtOrBin->Add(fHistoGammaTruePrimaryCaloConvPtOrBin);

                // check for new trainoutput, secondary historgrams in 2D
                fHistogramDimension                                                         = TrueConversionContainer->FindObject("ESD_TrueSecondaryClusGamma_Pt")->ClassName();
                if (fHistogramDimension.Contains("2")){
                    nHistogramDimension                                                     = 2;
                    cout << "Found 2D histgrams -> using new trainoutput" << endl;
                }

                if (nHistogramDimension==2) {
 
                    // true secondary gammas in rec Pt
                    f2DHistoGammaTrueSecondaryCaloUnConvPt                                  = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusGamma_Pt");
                    f2DHistoGammaTrueSecondaryCaloUnConvPt->SetName("ESD_2D_TrueSecondaryClusUnConvGamma_Pt");
                    f2DHistoGammaTrueSecondaryCaloConvPt                                    = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusConvGamma_Pt");
                    f2DHistoGammaTrueSecondaryCaloConvPt->SetName("ESD_2D_TrueSecondaryClusConvGamma_Pt");
                    f2DHistoGammaTrueSecondaryCaloPt                                        = (TH2D*)f2DHistoGammaTrueSecondaryCaloUnConvPt->Clone("ESD_2D_TrueSecondaryClusGamma_Pt");
                    f2DHistoGammaTrueSecondaryCaloPt->Add(f2DHistoGammaTrueSecondaryCaloConvPt);
                    
                    // all
                    fHistoGammaTrueSecondaryCaloPt                                          = (TH1D*)f2DHistoGammaTrueSecondaryCaloPt->ProjectionX("ESD_TrueSecondaryClusGamma_Pt",1,5,"e");
                    fHistoGammaTrueSecondaryCaloPt->Sumw2();
                    fHistoGammaTrueSecondaryCaloPtOrBin                                     = (TH1D*)fHistoGammaTrueSecondaryCaloPt->Clone("ESD_TrueSecondaryClusGamma_PtOrBin");
                    RebinSpectrum(fHistoGammaTrueSecondaryCaloPt);
                    
                    // K0s
                    fHistoGammaTrueSecondaryCaloFromXFromK0sPt                              = (TH1D*)f2DHistoGammaTrueSecondaryCaloPt->ProjectionX("ESD_TrueSecondaryClusGammaFromXFromK0s_Pt",1,1,"e");
                    fHistoGammaTrueSecondaryCaloFromXFromK0sPt->Sumw2();
                    fHistoGammaTrueSecondaryCaloFromXFromK0sPtOrBin                         = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromK0sPt->Clone("ESD_TrueSecondaryClusGammaFromXFromK0s_PtOrBin");
                    RebinSpectrum(fHistoGammaTrueSecondaryCaloFromXFromK0sPt);
                    
                    // K0l
                    fHistoGammaTrueSecondaryCaloFromXFromK0lPt                              = (TH1D*)f2DHistoGammaTrueSecondaryCaloPt->ProjectionX("ESD_TrueSecondaryClusGammaFromXFromK0l_Pt",2,2,"e");
                    fHistoGammaTrueSecondaryCaloFromXFromK0lPt->Sumw2();
                    fHistoGammaTrueSecondaryCaloFromXFromK0lPtOrBin                         = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromK0lPt->Clone("ESD_TrueSecondaryClusGammaFromXFromK0l_PtOrBin");
                    RebinSpectrum(fHistoGammaTrueSecondaryCaloFromXFromK0lPt);

                    // Lambda
                    fHistoGammaTrueSecondaryCaloFromXFromLambdaPt                           = (TH1D*)f2DHistoGammaTrueSecondaryCaloPt->ProjectionX("ESD_TrueSecondaryClusGammaFromXFromLambda_Pt",3,3,"e");
                    fHistoGammaTrueSecondaryCaloFromXFromLambdaPt->Sumw2();
                    fHistoGammaTrueSecondaryCaloFromXFromLambdaPtOrBin                      = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromLambdaPt->Clone("ESD_TrueSecondaryClusGammaFromXFromLambda_PtOrBin");
                    RebinSpectrum(fHistoGammaTrueSecondaryCaloFromXFromLambdaPt);
                    
                    // Eta
                    fHistoGammaTrueSecondaryCaloFromXFromEtaPt                              = (TH1D*)f2DHistoGammaTrueSecondaryCaloPt->ProjectionX("ESD_TrueSecondaryClusGammaFromXFromEta_Pt",4,4,"e");
                    fHistoGammaTrueSecondaryCaloFromXFromEtaPt->Sumw2();
                    fHistoGammaTrueSecondaryCaloFromXFromEtaPtOrBin                         = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromEtaPt->Clone("ESD_TrueSecondaryClusGammaFromXFromEta_PtOrBin");
                    RebinSpectrum(fHistoGammaTrueSecondaryCaloFromXFromEtaPt);
                    
                    // Rest (including Eta)
                    fHistoGammaTrueSecondaryCaloRestPt                                      = (TH1D*)f2DHistoGammaTrueSecondaryCaloPt->ProjectionX("ESD_TrueSecondaryClusGammaRest_Pt",4,5,"e");
                    fHistoGammaTrueSecondaryCaloRestPt->Sumw2();
                    fHistoGammaTrueSecondaryCaloRestPtOrBin                                 = (TH1D*)fHistoGammaTrueSecondaryCaloRestPt->Clone("ESD_TrueSecondaryClusGammaRest_PtOrBin");
                    RebinSpectrum(fHistoGammaTrueSecondaryCaloRestPt);

                    // Rest (without Eta)
                    fHistoGammaTrueSecondaryCaloRestNoEtaPt                                 = (TH1D*)f2DHistoGammaTrueSecondaryCaloPt->ProjectionX("ESD_TrueSecondaryClusGammaRestNoEta_Pt",5,5,"e");
                    fHistoGammaTrueSecondaryCaloRestNoEtaPt->Sumw2();
                    fHistoGammaTrueSecondaryCaloRestNoEtaPtOrBin                            = (TH1D*)fHistoGammaTrueSecondaryCaloRestNoEtaPt->Clone("ESD_TrueSecondaryClusGammaRestNoEta_PtOrBin");
                    RebinSpectrum(fHistoGammaTrueSecondaryCaloRestNoEtaPt);

                    // used to calculate raw spectra from cocktail
                    if (fUseCocktail) {
                        
                        // true secondary gammas in MC Pt
                        f2DHistoGammaTrueSecondaryCaloUnConvMCPt                            = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusGamma_MCPt");
                        f2DHistoGammaTrueSecondaryCaloUnConvMCPt->SetName("ESD_2D_TrueSecondaryClusUnConvGamma_MCPt");
                        f2DHistoGammaTrueSecondaryCaloConvMCPt                              = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusConvGamma_MCPt");
                        f2DHistoGammaTrueSecondaryCaloConvMCPt->SetName("ESD_2D_TrueSecondaryClusUnConvGamma_MCPt");
                        f2DHistoGammaTrueSecondaryCaloMCPt                                  = (TH2D*)f2DHistoGammaTrueSecondaryCaloUnConvMCPt->Clone("ESD_2D_TrueSecondaryClusGamma_MCPt");
                        f2DHistoGammaTrueSecondaryCaloMCPt->Add(f2DHistoGammaTrueSecondaryCaloConvMCPt);

                        // All
                        fHistoGammaTrueSecondaryCaloMCPt                                    = (TH1D*)f2DHistoGammaTrueSecondaryCaloMCPt->ProjectionX("ESD_TrueSecondaryClusGamma_MCPt",1,5,"e");
                        fHistoGammaTrueSecondaryCaloMCPt->Sumw2();
                        fHistoGammaTrueSecondaryCaloMCPtOrBin                               = (TH1D*)fHistoGammaTrueSecondaryCaloMCPt->Clone("ESD_TrueSecondaryClusGamma_MCPtOrBin");
                        RebinSpectrum(fHistoGammaTrueSecondaryCaloMCPt);
                        
                        // K0s
                        fHistoGammaTrueSecondaryCaloFromXFromK0sMCPt                        = (TH1D*)f2DHistoGammaTrueSecondaryCaloMCPt->ProjectionX("ESD_TrueSecondaryClusGammaFromXFromK0s_MCPt",1,1,"e");
                        fHistoGammaTrueSecondaryCaloFromXFromK0sMCPt->Sumw2();
                        fHistoGammaTrueSecondaryCaloFromXFromK0sMCPtOrBin                   = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromK0sMCPt->Clone("ESD_TrueSecondaryClusGammaFromXFromK0s_MCPtOrBin");
                        RebinSpectrum(fHistoGammaTrueSecondaryCaloFromXFromK0sMCPt);
                        
                        // K0l
                        fHistoGammaTrueSecondaryCaloFromXFromK0lMCPt                        = (TH1D*)f2DHistoGammaTrueSecondaryCaloMCPt->ProjectionX("ESD_TrueSecondaryClusGammaFromXFromK0l_MCPt",2,2,"e");
                        fHistoGammaTrueSecondaryCaloFromXFromK0lMCPt->Sumw2();
                        fHistoGammaTrueSecondaryCaloFromXFromK0lMCPtOrBin                   = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromK0lMCPt->Clone("ESD_TrueSecondaryClusGammaFromXFromK0l_MCPtOrBin");
                        RebinSpectrum(fHistoGammaTrueSecondaryCaloFromXFromK0lMCPt);

                        // Lambda
                        fHistoGammaTrueSecondaryCaloFromXFromLambdaMCPt                     = (TH1D*)f2DHistoGammaTrueSecondaryCaloMCPt->ProjectionX("ESD_TrueSecondaryClusGammaFromXFromLambda_MCPt",3,3,"e");
                        fHistoGammaTrueSecondaryCaloFromXFromLambdaMCPt->Sumw2();
                        fHistoGammaTrueSecondaryCaloFromXFromLambdaMCPtOrBin                = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromLambdaMCPt->Clone("ESD_TrueSecondaryClusGammaFromXFromLambda_MCPtOrBin");
                        RebinSpectrum(fHistoGammaTrueSecondaryCaloFromXFromLambdaMCPt);
                        
                        // Eta
                        fHistoGammaTrueSecondaryCaloFromXFromEtaMCPt                        = (TH1D*)f2DHistoGammaTrueSecondaryCaloMCPt->ProjectionX("ESD_TrueSecondaryClusGammaFromXFromEta_MCPt",4,4,"e");
                        fHistoGammaTrueSecondaryCaloFromXFromEtaMCPt->Sumw2();
                        fHistoGammaTrueSecondaryCaloFromXFromEtaMCPtOrBin                   = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromEtaMCPt->Clone("ESD_TrueSecondaryClusGammaFromXFromEta_MCPtOrBin");
                        RebinSpectrum(fHistoGammaTrueSecondaryCaloFromXFromEtaMCPt);
                        
                        // Rest (including Eta)
                        fHistoGammaTrueSecondaryCaloRestMCPt                                = (TH1D*)f2DHistoGammaTrueSecondaryCaloMCPt->ProjectionX("ESD_TrueSecondaryClusGammaRest_MCPt",4,5,"e");
                        fHistoGammaTrueSecondaryCaloRestMCPt->Sumw2();
                        fHistoGammaTrueSecondaryCaloRestMCPtOrBin                           = (TH1D*)fHistoGammaTrueSecondaryCaloRestMCPt->Clone("ESD_TrueSecondaryClusGammaRest_MCPtOrBin");
                        RebinSpectrum(fHistoGammaTrueSecondaryCaloRestMCPt);

                        // Rest (without Eta)
                        fHistoGammaTrueSecondaryCaloRestNoEtaMCPt                           = (TH1D*)f2DHistoGammaTrueSecondaryCaloMCPt->ProjectionX("ESD_TrueSecondaryClusGammaRest_MCPt",5,5,"e");
                        fHistoGammaTrueSecondaryCaloRestNoEtaMCPt->Sumw2();
                        fHistoGammaTrueSecondaryCaloRestNoEtaMCPtOrBin                      = (TH1D*)fHistoGammaTrueSecondaryCaloRestNoEtaMCPt->Clone("ESD_TrueSecondaryClusGammaRestNoEta_MCPtOrBin");
                        RebinSpectrum(fHistoGammaTrueSecondaryCaloRestNoEtaMCPt);

                        // secondary response matrices
                        fHistoGammaTrueSecondaryFromXFromK0sCaloUnConv_MCPt_recPt_MC        = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusGammaFromXFromK0s_MCPt_Pt");
                        fHistoGammaTrueSecondaryFromXFromK0sCaloUnConv_MCPt_recPt_MC->SetName("ESD_TrueSecondaryClusUnConvGammaFromXFromK0s_MCPt_Pt");
                        fHistoGammaTrueSecondaryFromXFromK0sCaloConv_MCPt_recPt_MC          = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusConvGammaFromXFromK0s_MCPt_Pt");
                        fHistoGammaTrueSecondaryFromXFromK0sCalo_MCPt_recPt_MC              = (TH2D*)fHistoGammaTrueSecondaryFromXFromK0sCaloUnConv_MCPt_recPt_MC->Clone("ESD_TrueSecondaryClusGammaFromXFromK0s_MCPt_recPt");
                        fHistoGammaTrueSecondaryFromXFromK0sCalo_MCPt_recPt_MC->Add(fHistoGammaTrueSecondaryFromXFromK0sCaloConv_MCPt_recPt_MC);
                        
                        fHistoGammaTrueSecondaryFromXFromK0lCaloUnConv_MCPt_recPt_MC        = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusGammaFromXFromK0l_MCPt_Pt");
                        fHistoGammaTrueSecondaryFromXFromK0lCaloUnConv_MCPt_recPt_MC->SetName("ESD_TrueSecondaryClusUnConvGammaFromXFromK0l_MCPt_Pt");
                        fHistoGammaTrueSecondaryFromXFromK0lCaloConv_MCPt_recPt_MC          = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusConvGammaFromXFromK0l_MCPt_Pt");
                        fHistoGammaTrueSecondaryFromXFromK0lCalo_MCPt_recPt_MC              = (TH2D*)fHistoGammaTrueSecondaryFromXFromK0lCaloUnConv_MCPt_recPt_MC->Clone("ESD_TrueSecondaryClusGammaFromXFromK0l_MCPt_recPt");
                        fHistoGammaTrueSecondaryFromXFromK0lCalo_MCPt_recPt_MC->Add(fHistoGammaTrueSecondaryFromXFromK0lCaloConv_MCPt_recPt_MC);

                        fHistoGammaTrueSecondaryFromXFromLambdaCaloUnConv_MCPt_recPt_MC     = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusGammaFromXFromLambda_MCPt_Pt");
                        fHistoGammaTrueSecondaryFromXFromLambdaCaloUnConv_MCPt_recPt_MC->SetName("ESD_TrueSecondaryClusUnConvGammaFromXFromLambda_MCPt_Pt");
                        fHistoGammaTrueSecondaryFromXFromLambdaCaloConv_MCPt_recPt_MC       = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusConvGammaFromXFromLambda_MCPt_Pt");
                        fHistoGammaTrueSecondaryFromXFromLambdaCalo_MCPt_recPt_MC           = (TH2D*)fHistoGammaTrueSecondaryFromXFromLambdaCaloUnConv_MCPt_recPt_MC->Clone("ESD_TrueSecondaryClusGammaFromXFromLambda_MCPt_recPt");
                        fHistoGammaTrueSecondaryFromXFromLambdaCalo_MCPt_recPt_MC->Add(fHistoGammaTrueSecondaryFromXFromLambdaCaloConv_MCPt_recPt_MC);
                    }
                } else {
                    // MC true secondary cluster gammas
                    fHistoGammaTrueSecondaryCaloUnConvPt                                    = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusGamma_Pt");
                    fHistoGammaTrueSecondaryCaloUnConvPt->SetName("TrueSecondaryClusUnConvGamma_Pt");
                    fHistoGammaTrueSecondaryCaloUnConvPt->Sumw2();
                    fHistoGammaTrueSecondaryCaloUnConvPtOrBin                               = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvPt->Clone("TrueSecondaryClusUnConvGamma_Pt_OriginalBinning");
                    fHistoGammaTrueSecondaryCaloUnConvPtOrBin->Sumw2();
                    RebinSpectrum(fHistoGammaTrueSecondaryCaloUnConvPt);
                    
                    fHistoGammaTrueSecondaryCaloConvPt                                      = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusConvGamma_Pt");
                    fHistoGammaTrueSecondaryCaloConvPt->Sumw2();
                    fHistoGammaTrueSecondaryCaloConvPtOrBin                                 = (TH1D*)fHistoGammaTrueSecondaryCaloConvPt->Clone("TrueSecondaryClusConvGamma_Pt_OriginalBinning");
                    fHistoGammaTrueSecondaryCaloConvPtOrBin->Sumw2();
                    RebinSpectrum(fHistoGammaTrueSecondaryCaloConvPt);
                    
                    fHistoGammaTrueSecondaryCaloPt                                          = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvPt->Clone("TrueSecondaryClusGamma_Pt");
                    fHistoGammaTrueSecondaryCaloPt->Add(fHistoGammaTrueSecondaryCaloConvPt);
                    fHistoGammaTrueSecondaryCaloPtOrBin                                     = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvPtOrBin->Clone("TrueSecondaryClusGamma_Pt_OriginalBinning");
                    fHistoGammaTrueSecondaryCaloPtOrBin->Add(fHistoGammaTrueSecondaryCaloConvPtOrBin);
                    
                    // MC true secondary cluster gammas from X from K0s
                    fHistoGammaTrueSecondaryCaloUnConvFromXFromK0sPt                        = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusGammaFromXFromK0s_Pt");
                    fHistoGammaTrueSecondaryCaloUnConvFromXFromK0sPt->SetName("TrueSecondaryClusUnConvGammaFromK0s_Pt");
                    fHistoGammaTrueSecondaryCaloUnConvFromXFromK0sPt->Sumw2();
                    fHistoGammaTrueSecondaryCaloUnConvFromXFromK0sPtOrBin                   = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvFromXFromK0sPt->Clone("TrueSecondaryClusUnConvGammaFromK0s_Pt_OriginalBinning");
                    fHistoGammaTrueSecondaryCaloUnConvFromXFromK0sPtOrBin->Sumw2();
                    RebinSpectrum(fHistoGammaTrueSecondaryCaloUnConvFromXFromK0sPt);
                    
                    fHistoGammaTrueSecondaryCaloConvFromXFromK0sPt                          = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusConvGammaFromXFromK0s_Pt");
                    fHistoGammaTrueSecondaryCaloConvFromXFromK0sPt->Sumw2();
                    fHistoGammaTrueSecondaryCaloConvFromXFromK0sPtOrBin                     = (TH1D*)fHistoGammaTrueSecondaryCaloConvFromXFromK0sPt->Clone("TrueSecondaryClusConvGammaFromK0s_Pt_OriginalBinning");
                    fHistoGammaTrueSecondaryCaloConvFromXFromK0sPtOrBin->Sumw2();
                    RebinSpectrum(fHistoGammaTrueSecondaryCaloConvFromXFromK0sPt);
                    
                    fHistoGammaTrueSecondaryCaloFromXFromK0sPt                              = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvFromXFromK0sPt->Clone("TrueSecondaryClusGammaFromK0s_Pt");
                    fHistoGammaTrueSecondaryCaloFromXFromK0sPt->Add(fHistoGammaTrueSecondaryCaloConvFromXFromK0sPt);
                    fHistoGammaTrueSecondaryCaloFromXFromK0sPtOrBin                         = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvFromXFromK0sPtOrBin->Clone("TrueSecondaryClusGammaFromK0s_Pt_OriginalBinning");
                    fHistoGammaTrueSecondaryCaloFromXFromK0sPtOrBin->Add(fHistoGammaTrueSecondaryCaloConvFromXFromK0sPtOrBin);

                    // MC true secondary cluster gammas from X from Lambda
                    fHistoGammaTrueSecondaryCaloUnConvFromXFromLambdaPt                     = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusGammaFromXFromLambda_Pt");
                    fHistoGammaTrueSecondaryCaloUnConvFromXFromLambdaPt->SetName("TrueSecondaryClusUnConvGammaFromLambda_Pt");
                    fHistoGammaTrueSecondaryCaloUnConvFromXFromLambdaPt->Sumw2();
                    fHistoGammaTrueSecondaryCaloUnConvFromXFromLambdaPtOrBin                = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvFromXFromLambdaPt->Clone("TrueSecondaryClusUnConvGammaFromLambda_Pt_OriginalBinning");
                    fHistoGammaTrueSecondaryCaloUnConvFromXFromLambdaPtOrBin->Sumw2();
                    RebinSpectrum(fHistoGammaTrueSecondaryCaloUnConvFromXFromLambdaPt);
                    
                    fHistoGammaTrueSecondaryCaloConvFromXFromLambdaPt                       = (TH1D*)TrueConversionContainer->FindObject("ESD_TrueSecondaryClusConvGammaFromXFromLambda_Pt");
                    fHistoGammaTrueSecondaryCaloConvFromXFromLambdaPt->Sumw2();
                    fHistoGammaTrueSecondaryCaloConvFromXFromLambdaPtOrBin                  = (TH1D*)fHistoGammaTrueSecondaryCaloConvFromXFromLambdaPt->Clone("TrueSecondaryClusConvGammaFromLambda_Pt_OriginalBinning");
                    fHistoGammaTrueSecondaryCaloConvFromXFromLambdaPtOrBin->Sumw2();
                    RebinSpectrum(fHistoGammaTrueSecondaryCaloConvFromXFromLambdaPt);
                    
                    fHistoGammaTrueSecondaryCaloFromXFromLambdaPt                           = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvFromXFromLambdaPt->Clone("TrueSecondaryClusGammaFromXFromLambda_Pt");
                    fHistoGammaTrueSecondaryCaloFromXFromLambdaPt->Add(fHistoGammaTrueSecondaryCaloConvFromXFromLambdaPt);
                    fHistoGammaTrueSecondaryCaloFromXFromLambdaPtOrBin                      = (TH1D*)fHistoGammaTrueSecondaryCaloUnConvFromXFromLambdaPtOrBin->Clone("TrueSecondaryClusGammaFromLambda_Pt_OriginalBinning");
                    fHistoGammaTrueSecondaryCaloFromXFromLambdaPtOrBin->Add(fHistoGammaTrueSecondaryCaloConvFromXFromLambdaPtOrBin);
                }

                // response matrix calo
                fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC                              = (TH2D*)TrueConversionContainer->FindObject("TruePrimaryClusGamma_Pt_MCPt");
                fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->SetName("TruePrimaryClusUnConvGamma_Pt_MCPt");
                fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->Sumw2();

                fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC                                = (TH2D*)TrueConversionContainer->FindObject("TruePrimaryClusConvGamma_Pt_MCPt");
                fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->Sumw2();

                fHistoGammaTruePrimaryCalo_recPt_MCPt_MC                                    = (TH2D*)fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->Clone("TruePrimaryClusGamma_Pt_MCPt");
                fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->Add(fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC);

                // true primary calo gamma
                fHistoGammaTruePrimaryCaloUnConvMCPt                                        = (TH1D*)fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->ProjectionY("ESD_TruePrimaryCaloUnConvGamma_MCPt");
                fHistoGammaTruePrimaryCaloUnConvMCPt->Sumw2();
                fHistoGammaTruePrimaryCaloUnConvMCPtOrBin                                   = (TH1D*)fHistoGammaTruePrimaryCaloUnConvMCPt->Clone("ESD_TruePrimaryCaloUnConvGamma_MCPt_OriginalBinning");
                RebinSpectrum(fHistoGammaTruePrimaryCaloUnConvMCPt);

                fHistoGammaTruePrimaryCaloConvMCPt                                          = (TH1D*)fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->ProjectionY("ESD_TruePrimaryCaloConvGamma_MCPt");
                fHistoGammaTruePrimaryCaloConvMCPt->Sumw2();
                fHistoGammaTruePrimaryCaloConvMCPtOrBin                                     = (TH1D*)fHistoGammaTruePrimaryCaloConvMCPt->Clone("ESD_TruePrimaryCaloConvGamma_MCPt_OriginalBinning");
                RebinSpectrum(fHistoGammaTruePrimaryCaloConvMCPt);
                
                fHistoGammaTruePrimaryCaloMCPt                                              = (TH1D*)fHistoGammaTruePrimaryCaloUnConvMCPt->Clone("ESD_TruePrimaryCaloGamma_MCPt");
                fHistoGammaTruePrimaryCaloMCPt->Add(fHistoGammaTruePrimaryCaloConvMCPt);

                fHistoGammaTruePrimaryCaloMCPtOrBin                                         = (TH1D*)fHistoGammaTruePrimaryCaloUnConvMCPtOrBin->Clone("ESD_TruePrimaryCaloGamma_MCPt_OriginalBinning");
                fHistoGammaTruePrimaryCaloMCPtOrBin->Add(fHistoGammaTruePrimaryCaloConvMCPtOrBin);

                // combinatorial BG distributions
                fHistoCombinatorialBackground                                               = (TH2D*)TrueConversionContainer->FindObject("ESD_TrueClusPhotonPlusConvBG_Pt");
                fHistoCombinatorialSpecies                                                  = new TH1D*[11];
                fHistoCombinatorialSpecies[nContamination]                                  = (TH1D*)fHistoCombinatorialBackground->ProjectionX(Form("ESD_TrueComb%s_Pt",fContamination[nContamination].Data()));
                fHistoCombinatorialSpecies[nContamination]->Sumw2();
                RebinSpectrum(fHistoCombinatorialSpecies[nContamination]);
                for(Int_t i = 0; i<nContamination; i++){
                    fHistoCombinatorialSpecies[i]                                           = (TH1D*)fHistoCombinatorialBackground->ProjectionX(Form("ESD_TrueComb%s_Pt",fContamination[i].Data()),i+1,i+1);
                    fHistoCombinatorialSpecies[i]->Sumw2();
                    RebinSpectrum(fHistoCombinatorialSpecies[i]);
                }
            }    
        }
        
        if (mode == 2 || mode == 3){
            for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){ // BEGIN ANALYSIS for each Pt bin
                cout << "---------------------------------------------------------------------------------" << endl;
                cout << "Begin Analysis Pt Bin " << iPt <<endl;
                cout << "---------------------------------------------------------------------------------" << endl;
                // Function to subtract GG minus Bck
                fFileDataLog << "---------------------------------------------------------------------------------" << endl;
                fFileDataLog << "----------------------------------new pT bin ------------------------------------" << endl;
                fFileDataLog << "---------------------------------------------------------------------------------" << endl;
                
                ProcessEM( fHistoGGInvMassPtGConvBin[iPt], fHistoBackInvMassPtGconvBin[iPt], fBGFitRange);
                fHistoSignalInvMassPtGConvBin[iPt]          = fSignal;
                fHistoBackNormInvMassPtGconvBin[iPt]        = fBckNorm;

                fHistoSignalInvMassPtGConvBin[iPt]->SetName(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",iPt));
                
                // Fitting the subtracted spectra
                fFileDataLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "normal range/right normalization" << endl;
                fFitSignalInvMassPtBin[iPt]                 = 0x00;
                fFileDataLog << "Using exp fit"<<endl;
                fFileDataLog << "Subtracted mixed event" << endl;
                FitSubtractedInvMassInPtBins(fHistoSignalInvMassPtGConvBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
                fFitSignalInvMassPtBin[iPt]                 = fFitReco;
                fFitBckInvMassPtBin[iPt]                    = fFitLinearBck;
                fMesonYieldsResidualBckFunc[0][iPt]         = fIntLinearBck;
                fMesonYieldsResidualBckFuncError[0][iPt]    = fIntLinearBckError;
                if (fFitReco)
                    fMesonChi2[0][iPt]                      = fFitReco->GetChisquare()/fFitReco->GetNDF();
                else 
                    fMesonChi2[0][iPt]                      = -1;                
                
                //Get FWHM
                CalculateFWHM( fFitSignalInvMassPtBin[iPt]);
                fMesonFWHM[iPt]                     = fFWHMFunc;
                fMesonFWHMError[iPt]                = fFWHMFuncError;

                if (fFitSignalInvMassPtBin[iPt] !=0x00){
                    fMesonMass[iPt]                 = fFitSignalInvMassPtBin[iPt]->GetParameter(1);
                    fMesonMassError[iPt]            = fFitSignalInvMassPtBin[iPt]->GetParError(1);
                    
                    fMesonLambdaTailpar[iPt]        = fFitSignalInvMassPtBin[iPt]->GetParameter(3);
                    fMesonLambdaTailparError[iPt]   = fFitSignalInvMassPtBin[iPt]->GetParError(3);

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
                
                for (Int_t k = 0; k < 3; k++){
                    IntegrateHistoInvMass( fHistoGGInvMassPtGConvBin[iPt], fMesonCurIntRange[k]);
                    fGGYields[k][iPt]               = fYields;
                    fGGYieldsError[k][iPt]          = fYieldsError;
                    
                    // Integrate the bck histo
                    IntegrateHistoInvMass( fHistoBackNormInvMassPtGconvBin[iPt], fMesonCurIntRange[k]);
                    fBckYields[k][iPt]              = fYields;
                    fBckYieldsError[k][iPt]         = fYieldsError;
                    
                    // Integrate the signal histo
                    fFileDataLog<< endl <<"Signal histo "<< nameIntRange[k].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
                    IntegrateHistoInvMassStream( fHistoSignalInvMassPtGConvBin[iPt], fMesonCurIntRange[k]);
                    fMesonYields[k][iPt]            = fYields;
                    fMesonYieldsError[k][iPt]       = fYieldsError;
                    fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
                }    
            
                for (Int_t k = 0; k < 3; k++){
                    if (k > 0){
                        fFileDataLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << nameIntRange[k].Data()<<  endl;
                        Double_t intRange[2]    = {0,0};
                        if (k == 1){
                            intRange[0]         = fMesonIntDeltaRangeWide[0];
                            intRange[1]         = fMesonIntDeltaRangeWide[1];
                        } else if ( k == 2){
                            intRange[0]         = fMesonIntDeltaRangeNarrow[0];
                            intRange[1]         = fMesonIntDeltaRangeNarrow[1];                    
                        }    
                        fFileDataLog << "Using exp fit"<<endl;
                        FitSubtractedInvMassInPtBins(fHistoSignalInvMassPtGConvBin[iPt], intRange, iPt, kFALSE);
                        fMesonYieldsResidualBckFunc[k][iPt]         = fIntLinearBck;
                        fMesonYieldsResidualBckFuncError[k][iPt]    = fIntLinearBckError;
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
                    fMesonYieldsCorResidualBckFuncError[k][iPt]     = pow(( fMesonYieldsError[k][iPt]*fMesonYieldsError[k][iPt] +
                                                                        fMesonYieldsResidualBckFuncError[k][iPt]*fMesonYieldsResidualBckFuncError[k][iPt]),0.5);
                    fMesonYieldsPerEvent[k][iPt]                    = fMesonYieldsCorResidualBckFunc[k][iPt]/fNEvents;
                    fMesonYieldsPerEventError[k][iPt]               = fMesonYieldsCorResidualBckFuncError[k][iPt]/fNEvents;
                    
                }    


                //////////////////////////////// Start Analysis with  Normalization at the left of the Meson Peak
                // Function to subtract GG minus Bck
                fFileDataLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << endl;
                ProcessEM( fHistoGGInvMassPtGConvBin[iPt], fHistoBackInvMassPtGconvBin[iPt], fBGFitRangeLeft);
                fHistoSignalInvMassLeftPtGConvBin[iPt]          = fSignal;
                fHistoBackNormInvMassLeftPtGconvBin[iPt]        = fBckNorm;


                fHistoSignalInvMassLeftPtGConvBin[iPt]->SetName(Form("fHistoMappingSignalLeftInvMass_in_Pt_Bin%02d",iPt));
                //       fFileDataLog<< "iPt"<< iPt<< " "<< "standard range"<<endl;
                // Fitting the subtracted spectra
                fFileDataLog  << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << "normal range/left normalization" << endl;

                fFitInvMassLeftPtBin[iPt] =0x00;
                fFileDataLog << "Using exp fit"<<endl;
                FitSubtractedInvMassInPtBins(fHistoSignalInvMassLeftPtGConvBin[iPt], fMesonIntDeltaRange,iPt,kFALSE);
                fFitInvMassLeftPtBin[iPt]                               = fFitReco;
                fFitBckInvMassLeftPtBin[iPt]                            = fFitLinearBck;
                fMesonYieldsResidualBckFunc[3][iPt]                     = fIntLinearBck;
                fMesonYieldsResidualBckFuncError[3][iPt]                = fIntLinearBckError;
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
                for (Int_t k = 0; k < 3; k++){
                    IntegrateHistoInvMass( fHistoBackNormInvMassLeftPtGconvBin[iPt], fMesonCurIntRange[k+3]);
                    fBckYields[k+3][iPt]              = fYields;
                    fBckYieldsError[k+3][iPt]         = fYieldsError;

                    fFileDataLog<< endl <<"Signal histo " << nameIntRange[k+3].Data() << ":\t" << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1]<< endl;
                    IntegrateHistoInvMassStream( fHistoSignalInvMassLeftPtGConvBin[iPt], fMesonCurIntRange[k+3]);
                    fMesonYields[k+3][iPt]            = fYields;
                    fMesonYieldsError[k+3][iPt]       = fYieldsError;
                    fFileDataLog << "Integrated value: \t" << fYields <<"+-" <<fYieldsError <<endl;
                    
                }
            

                for (Int_t k = 0; k < 3; k++){
                    if (k > 0){
                        fFileDataLog << fBinsPt[iPt] <<"-" << fBinsPt[iPt+1] << "\t" << nameIntRange[k+3].Data()<<  endl;
                        Double_t intRange[2]    = {0,0};
                        if (k == 1){
                            intRange[0]         = fMesonIntDeltaRangeWide[0];
                            intRange[1]         = fMesonIntDeltaRangeWide[1];
                        } else if ( k == 2){
                            intRange[0]         = fMesonIntDeltaRangeNarrow[0];
                            intRange[1]         = fMesonIntDeltaRangeNarrow[1];                    
                        }    
                        fFileDataLog << "Using exp fit"<<endl;
                        FitSubtractedInvMassInPtBins(fHistoSignalInvMassLeftPtGConvBin[iPt], intRange, iPt, kFALSE);
                        fMesonYieldsResidualBckFunc[k+3][iPt]         = fIntLinearBck;
                        fMesonYieldsResidualBckFuncError[k+3][iPt]    = fIntLinearBckError;
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
            TString plotPrefix  = Form("%s/%s_%s",fOutputDir.Data(),fPrefix.Data(),fPrefix2.Data());
            TString plotSuffix  = Form("%s_%s.%s",fPeriodFlag.Data(),fCutSelection.Data(),fSuffix.Data());
            
            TString nameMeson   = Form("%s_GammaConvPt_MesonWithBck%s", plotPrefix.Data(), plotSuffix.Data());
            TString nameCanvas  = "MesonWithBckCanvas";
            TString namePad     = "MesonWithBckPad";
            cout << nameMeson.Data() << endl;
            PlotInvMassInPtBins( fHistoGGInvMassPtGConvBin, fHistoBackNormInvMassPtGconvBin, nameMeson, nameCanvas, namePad, fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, 
                                fBinsPt, fTextMeasurement, fIsMC ,fDecayChannel, fDetectionProcess, fCollisionSystem);

            nameMeson       = Form("%s_GammaConvPt_MesonWithBckLeft%s", plotPrefix.Data(), plotSuffix.Data());
            nameCanvas      = "MesonWithBckCanvasLeft";
            namePad         = "MesonWithBckPadLeft";
            cout << nameMeson.Data() << endl;
            PlotInvMassInPtBins( fHistoGGInvMassPtGConvBin, fHistoBackNormInvMassLeftPtGconvBin, nameMeson, nameCanvas, namePad,  fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin,
                                fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess, fCollisionSystem);

            TString nameMesonSub    = Form("%s_GammaConvPt_MesonSubtracted%s", plotPrefix.Data(), plotSuffix.Data());
            TString nameCanvasSub   = "MesonCanvasSubtracted";
            TString namePadSub      = "MesonPadSubtracted";
            cout << nameMesonSub.Data() << endl;
            PlotWithFitSubtractedInvMassInPtBins(   fHistoSignalInvMassPtGConvBin, NULL, fFitSignalInvMassPtBin, nameMesonSub, nameCanvasSub, namePadSub,
                                                    fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC,fDecayChannel, fDetectionProcess,
                                                    fCollisionSystem, "MC validated");
                        
            nameMesonSub    = Form("%s_GammaConvPt_MesonSubtractedLeft%s", plotPrefix.Data(), plotSuffix.Data());
            nameCanvasSub   = "MesonCanvasSubtractedLeft";
            namePadSub      = "MesonPadSubtractedLeft";
            cout << nameMesonSub.Data() << endl;
            PlotWithFitSubtractedInvMassInPtBins(   fHistoSignalInvMassLeftPtGConvBin, NULL, fFitInvMassLeftPtBin, nameMesonSub, nameCanvasSub, namePadSub,
                                                    fMesonMassPlotRange, fdate, fPrefix, fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fTextMeasurement, fIsMC, fDecayChannel, fDetectionProcess,
                                                    fCollisionSystem, "MC validated");

            cout << "Example bin: "<< fExampleBin << endl;
            TString triggerInt         = fEventCutSelectionRead(GetEventSelectSpecialTriggerCutPosition(),2);
            PlotExampleInvMassBinsV2(   fHistoGGInvMassPtGConvBin[fExampleBin], fHistoSignalInvMassPtGConvBin[fExampleBin], fHistoBackNormInvMassPtGconvBin[fExampleBin],
                                        fFitSignalInvMassPtBin[fExampleBin], fExampleBin, fOutputDir.Data(),fSuffix.Data(), fMesonMassPlotRange, pictDrawingCoordinatesFWHM, fNEvents, fdate, fPrefix, fPrefix2,
                                        fThesis, fCollisionSystem, fBinsPt, fDecayChannel, fDetectionProcess, triggerInt.Atoi(), fExampleBinScaleFac, fMode, addSig );

            
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
            if (catIter == 0)   category                                                    = 5;
            else                category                                                    = catIter;
            
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
        
        SetHistogramm(fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinning,"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.1);
        SetHistogramm(fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning,"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.1);
        SetHistogramm(fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning,"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.1);
        SetHistogramm(fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinning,"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.1);
        
        DrawGammaSetMarker(fMCrecGammaPtRatioWithWithoutPileUpDCAzDistBinning, 20, 1.0, kBlack, kBlack);
        DrawGammaSetMarker(fTruePrimaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning, 24, 1.0, kRed, kRed);
        DrawGammaSetMarker(fTrueSecondaryConvGammaPtRatioWithWithoutPileUpDCAzDistBinning, 25, 1.0, kBlue-9, kBlue-9);
        DrawGammaSetMarker(fTrueSecondaryFromXFromK0sConvGammaPtRatioWithWithoutPileUpDCAzDistBinning, 28, 1.0, kBlue+3, kBlue+3);

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
        fESDGammaPtPileUp->Multiply(fESDGammaPileUpCorrFactor[nPileupMethodUsed]);
        
        // plotting ratio + fit
        TCanvas *RatioWithWithoutPileUpCanvas                               = GetAndSetCanvas("canvasRatioWithWithoutPileUp");

        SetHistogramm(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat[0],"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.1);
        SetHistogramm(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[0],"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.1);
        SetHistogramm(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[1],"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.1);
        SetHistogramm(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[2],"p_{T} (GeV/c)","#gamma / #gamma Pile-Up correted (1/#it{C}_{pileup})",0.95,1.1);

        DrawGammaSetMarker(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinningAllCat[0], 28, 1.0, kGray+2, kGray+2);
        DrawGammaSetMarker(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[0], 20, 1.0, kBlack, kBlack);
        DrawGammaSetMarker(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[1], 24, 1.0, kBlue-2, kBlue-2);
        DrawGammaSetMarker(fESDGammaPtRatioWithWithoutPileUpDCAzDistBinning[2], 25, 1.0, kGreen+2, kGreen+2);

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
        TAxis *xAxis                                                = fHistoGammaTruePrimaryConv_recPt_MCPt_MC->GetXaxis();
        TAxis *yAxis                                                = fHistoGammaTruePrimaryConv_recPt_MCPt_MC->GetYaxis();

        // =========== Response matrix ==============
        fHistoGammaTruePrimaryConv_recPt_MCPt_MC_Rebin              = new TH2D("TruePrimaryConvGamma_recPt_MCPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
        for(Int_t x = 1; x<fHistoGammaTruePrimaryConv_recPt_MCPt_MC->GetNbinsX()+1; x++){
            for(Int_t y = 1; y<fHistoGammaTruePrimaryConv_recPt_MCPt_MC->GetNbinsY()+1; y++){
                Double_t binContent                                 = fHistoGammaTruePrimaryConv_recPt_MCPt_MC->GetBinContent(x,y);
                Double_t xcenter                                    = xAxis->GetBinCenter(x);
                Double_t ycenter                                    = yAxis->GetBinCenter(y);
                fHistoGammaTruePrimaryConv_recPt_MCPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
            }
        }
        
        // secondary response matrices
        if(fUseCocktail && nHistogramDimension==2){
            fHistoGammaTrueSecondaryFromXFromK0sConv_MCPt_recPt_MC_Rebin = new TH2D("TrueSecondaryFromXFromK0sConvGamma_MCPt_rectPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
            for(Int_t x = 1; x<fHistoGammaTrueSecondaryFromXFromK0sConv_MCPt_recPt_MC->GetNbinsX()+1; x++){
                for(Int_t y = 1; y<fHistoGammaTrueSecondaryFromXFromK0sConv_MCPt_recPt_MC->GetNbinsY()+1; y++){
                    Double_t binContent                             = fHistoGammaTrueSecondaryFromXFromK0sConv_MCPt_recPt_MC->GetBinContent(x,y);
                    Double_t xcenter                                = xAxis->GetBinCenter(x);
                    Double_t ycenter                                = yAxis->GetBinCenter(y);
                    fHistoGammaTrueSecondaryFromXFromK0sConv_MCPt_recPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
                }
            }
            fHistoGammaTrueSecondaryFromXFromK0lConv_MCPt_recPt_MC_Rebin = new TH2D("TrueSecondaryFromXFromK0sConvGamma_MCPt_rectPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
            for(Int_t x = 1; x<fHistoGammaTrueSecondaryFromXFromK0lConv_MCPt_recPt_MC->GetNbinsX()+1; x++){
                for(Int_t y = 1; y<fHistoGammaTrueSecondaryFromXFromK0lConv_MCPt_recPt_MC->GetNbinsY()+1; y++){
                    Double_t binContent                             = fHistoGammaTrueSecondaryFromXFromK0lConv_MCPt_recPt_MC->GetBinContent(x,y);
                    Double_t xcenter                                = xAxis->GetBinCenter(x);
                    Double_t ycenter                                = yAxis->GetBinCenter(y);
                    fHistoGammaTrueSecondaryFromXFromK0lConv_MCPt_recPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
                }
            }
            fHistoGammaTrueSecondaryFromXFromLambdaConv_MCPt_recPt_MC_Rebin = new TH2D("TrueSecondaryFromXFromK0sConvGamma_MCPt_rectPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
            for(Int_t x = 1; x<fHistoGammaTrueSecondaryFromXFromLambdaConv_MCPt_recPt_MC->GetNbinsX()+1; x++){
                for(Int_t y = 1; y<fHistoGammaTrueSecondaryFromXFromLambdaConv_MCPt_recPt_MC->GetNbinsY()+1; y++){
                    Double_t binContent                             = fHistoGammaTrueSecondaryFromXFromLambdaConv_MCPt_recPt_MC->GetBinContent(x,y);
                    Double_t xcenter                                = xAxis->GetBinCenter(x);
                    Double_t ycenter                                = yAxis->GetBinCenter(y);
                    fHistoGammaTrueSecondaryFromXFromLambdaConv_MCPt_recPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
                }
            }
        }
        // ==========================================
        
        // ======== Secondary fractions =============
        fHistoFracAllGammaToSecOrBin                                = (TH1D*)fHistoGammaConvPtOrBin->Clone("FracAllGammaToSecOriginalBinning");
        fHistoFracAllGammaToSecOrBin->Divide(fHistoGammaTrueSecondaryConvPtOrBin,fHistoFracAllGammaToSecOrBin,1,1,"B");
        fHistoFracAllGammaToSec                                     = (TH1D*) fHistoGammaConvPt->Clone("FracAllGammaToSec");
        fHistoFracAllGammaToSec->Divide(fHistoGammaTrueSecondaryConvPt,fHistoFracAllGammaToSec,1,1,"B");

        fHistoFracAllGammaToSecFromXFromK0s                         = (TH1D*)fHistoGammaConvPt->Clone("FracAllGammaToSecFromXFromK0s");
        fHistoFracAllGammaToSecFromXFromK0s->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt,fHistoFracAllGammaToSecFromXFromK0s,1,1,"B");

        fHistoFracAllGammaToSecFromXFromK0sOrBin                    = (TH1D*)fHistoGammaConvPtOrBin->Clone("FracAllGammaToSecFromXFromK0sOriginalBinning");
        fHistoFracAllGammaToSecFromXFromK0sOrBin->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromK0sPtOrBin,fHistoFracAllGammaToSecFromXFromK0sOrBin,1,1,"B");
        
        if(nHistogramDimension==2){
            fHistoFracAllGammaToSecFromXFromK0l                     = (TH1D*)fHistoGammaConvPt->Clone("FracAllGammaToSecFromXFromK0l");
            fHistoFracAllGammaToSecFromXFromK0l->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromK0lPt,fHistoFracAllGammaToSecFromXFromK0l,1,1,"B");

            fHistoFracAllGammaToSecFromXFromK0lOrBin                = (TH1D*)fHistoGammaConvPtOrBin->Clone("FracAllGammaToSecFromXFromK0lOriginalBinning");
            fHistoFracAllGammaToSecFromXFromK0lOrBin->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromK0lPtOrBin,fHistoFracAllGammaToSecFromXFromK0lOrBin,1,1,"B");
        }
        
        fHistoFracAllGammaToSecFromXFromLambda                      = (TH1D*)fHistoGammaConvPt->Clone("FracAllGammaToSecFromXFromLambda");
        fHistoFracAllGammaToSecFromXFromLambda->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt,fHistoFracAllGammaToSecFromXFromLambda,1,1,"B");

        fHistoFracAllGammaToSecFromXFromLambdaOrBin                 = (TH1D*)fHistoGammaConvPtOrBin->Clone("FracAllGammaToSecFromXFromLambdaOriginalBinning");
        fHistoFracAllGammaToSecFromXFromLambdaOrBin->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPtOrBin,fHistoFracAllGammaToSecFromXFromLambdaOrBin,1,1,"B");
        
        if(fUseCocktail){
            fHistoFracAllGammaToSecRest                             = (TH1D*)fHistoGammaConvPt->Clone("FracAllGammaToSecRest");
            fHistoFracAllGammaToSecRest->Divide(fHistoGammaTrueSecondaryConvGammaRestPt,fHistoFracAllGammaToSecRest,1,1,"B");
            
            fHistoFracAllGammaToSecRestOrBin                        = (TH1D*)fHistoGammaConvPtOrBin->Clone("FracAllGammaToSecRestOriginalBinning");
            fHistoFracAllGammaToSecRestOrBin->Divide(fHistoGammaTrueSecondaryConvGammaRestPtOrBin,fHistoFracAllGammaToSecRestOrBin,1,1,"B");
        }
        // ==========================================
        
        // =============== Conv Prob ================
        fHistoGammaMCConvProb                                       = new TH1D("MCGammaConvProb_MCPt","",fNBinsPt,fBinsPt);
        fHistoGammaMCConvProb->Sumw2();
        fHistoGammaMCConvProb->Divide(fHistoGammaMCConvPt,fHistoGammaMCAllPt,1,1,"B");

        fHistoGammaMCConvProbOrBin                                  = (TH1D*)fHistoGammaMCConvPtOrBin->Clone("MCGammaConvProb_MCPt_OriginalBinning");
        fHistoGammaMCConvProbOrBin->Sumw2();
        fHistoGammaMCConvProbOrBin->Divide(fHistoGammaMCConvProbOrBin,fHistoGammaMCAllPtOrBin,1,1,"B");

        // secondary conversion probabilities
        if(fUseCocktail && nHistogramDimension==2){
            // K0s
            fHistoSecondaryGammaFromXFromK0sMCConvProb              = new TH1D("SecondaryGammaFromXFromK0sMCGammaConvProb_MCPt","",fNBinsPt,fBinsPt);
            fHistoSecondaryGammaFromXFromK0sMCConvProb->Sumw2();
            fHistoSecondaryGammaFromXFromK0sMCConvProb->Divide(fHistoSecondaryGammaConvFromXFromK0sPt,fHistoAllSecondaryGammaFromXFromK0sPt,1,1,"B");

            fHistoSecondaryGammaFromXFromK0sMCConvProbOrBin         = (TH1D*)fHistoSecondaryGammaConvFromXFromK0sPtOrBin->Clone("SecondaryGammaFromXFromK0sMCGammaConvProb_MCPtOrBin");
            fHistoSecondaryGammaFromXFromK0sMCConvProbOrBin->Sumw2();
            fHistoSecondaryGammaFromXFromK0sMCConvProbOrBin->Divide(fHistoSecondaryGammaConvFromXFromK0sPtOrBin,fHistoAllSecondaryGammaFromXFromK0sPtOrBin,1,1,"B");

            // K0l
            fHistoSecondaryGammaFromXFromK0lMCConvProb              = new TH1D("SecondaryGammaFromXFromK0lMCGammaConvProb_MCPt","",fNBinsPt,fBinsPt);
            fHistoSecondaryGammaFromXFromK0lMCConvProb->Sumw2();
            fHistoSecondaryGammaFromXFromK0lMCConvProb->Divide(fHistoSecondaryGammaConvFromXFromK0lPt,fHistoAllSecondaryGammaFromXFromK0lPt,1,1,"B");
            
            fHistoSecondaryGammaFromXFromK0lMCConvProbOrBin         = (TH1D*)fHistoSecondaryGammaConvFromXFromK0lPtOrBin->Clone("SecondaryGammaFromXFromK0lMCGammaConvProb_MCPtOrBin");
            fHistoSecondaryGammaFromXFromK0lMCConvProbOrBin->Sumw2();
            fHistoSecondaryGammaFromXFromK0lMCConvProbOrBin->Divide(fHistoSecondaryGammaConvFromXFromK0lPtOrBin,fHistoAllSecondaryGammaFromXFromK0lPtOrBin,1,1,"B");

            // Lambda
            fHistoSecondaryGammaFromXFromLambdaMCConvProb           = new TH1D("SecondaryGammaFromXFromLambdaMCGammaConvProb_MCPt","",fNBinsPt,fBinsPt);
            fHistoSecondaryGammaFromXFromLambdaMCConvProb->Sumw2();
            fHistoSecondaryGammaFromXFromLambdaMCConvProb->Divide(fHistoSecondaryGammaConvFromXFromLambdaPt,fHistoAllSecondaryGammaFromXFromLambdaPt,1,1,"B");

            fHistoSecondaryGammaFromXFromLambdaMCConvProbOrBin      = (TH1D*)fHistoSecondaryGammaConvFromXFromLambdaPtOrBin->Clone("SecondaryGammaFromXFromLambdaMCGammaConvProb_MCPtOrBin");
            fHistoSecondaryGammaFromXFromLambdaMCConvProbOrBin->Sumw2();
            fHistoSecondaryGammaFromXFromLambdaMCConvProbOrBin->Divide(fHistoSecondaryGammaConvFromXFromLambdaPtOrBin,fHistoAllSecondaryGammaFromXFromLambdaPtOrBin,1,1,"B");
        }
        // ==========================================

        // ================= PURITY =================
        fHistoGammaMCPurity                                         = new TH1D("GammaPurity_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCPurity->Sumw2();
        fHistoGammaMCPurity->Divide(fHistoGammaTrueConvPt,fHistoGammaConvPt,1,1,"B");

        fHistoGammaMCrecPrimaryConvPt                               = (TH1D*) fHistoGammaConvPt->Clone("MCrec_PrimaryConvGamma_Pt");
        fHistoGammaMCrecPrimaryConvPt->Add(fHistoGammaTrueSecondaryConvPt,-1);
        fHistoGammaMCTruePurity                                     = new TH1D("GammaTruePurity_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCTruePurity->Sumw2();
        fHistoGammaMCTruePurity->Divide(fHistoGammaTruePrimaryConvPt,fHistoGammaMCrecPrimaryConvPt,1,1,"B");

        fHistoGammaMCrecPrimaryConvPtOrBin                          = (TH1D*) fHistoGammaConvPtOrBin->Clone("MC_ESDPrimaryConvGammaPt");
        fHistoGammaMCrecPrimaryConvPtOrBin->Add(fHistoGammaTrueSecondaryConvPtOrBin,-1);
        fHistoGammaMCTruePurityOrBin                                = (TH1D*)fHistoGammaMCrecPrimaryConvPtOrBin->Clone("GammaTruePurity_OriginalBinning_Pt");
        fHistoGammaMCTruePurityOrBin->Sumw2();
        fHistoGammaMCTruePurityOrBin->Divide(fHistoGammaTruePrimaryConvPtOrBin,fHistoGammaMCrecPrimaryConvPtOrBin,1,1,"B");
        // ==========================================

        // ================ Reco Eff ================
        fHistoGammaMCRecoEff                                        = new TH1D("GammaRecoEff_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCRecoEff->Sumw2();
        fHistoGammaMCRecoEff->Divide(fHistoGammaTrueConvPt,fHistoGammaMCConvPt,1,1,"B");

        fHistoGammaMCPrimaryRecoEff                                 = new TH1D("GammaPrimaryRecoEff_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCPrimaryRecoEff->Sumw2();
        fHistoGammaMCPrimaryRecoEff->Divide(fHistoGammaTruePrimaryConvPt,fHistoGammaMCConvPt,1,1,"B");

        TH1D* fHistoGammaTruePrimaryConvPtOrBinTemp                 = (TH1D*)fHistoGammaTruePrimaryConvPtOrBin->Clone("fHistoGammaTruePrimaryConvPtOrBinTemp");
        fHistoGammaTruePrimaryConvPtOrBinTemp->Scale(1./fHistoGammaTruePrimaryConvPtOrBinTemp->GetBinWidth(5));
        fHistoGammaMCPrimaryRecoEffOrBin                            = (TH1D*)fHistoGammaTruePrimaryConvPtOrBinTemp->Clone("GammaPrimaryRecoEff_Pt_OriginalBinning");
        fHistoGammaMCPrimaryRecoEffOrBin->Sumw2();
        fHistoGammaMCPrimaryRecoEffOrBin->Divide(fHistoGammaMCPrimaryRecoEffOrBin,fHistoGammaMCConvPtOrBin,1,1,"B");

        fHistoGammaMCPrimaryRecoEffMCPt                             = new TH1D("GammaPrimaryRecoEff_MCPt","",fNBinsPt,fBinsPt);
        fHistoGammaMCPrimaryRecoEffMCPt->Sumw2();
        fHistoGammaMCPrimaryRecoEffMCPt->Divide(fHistoGammaTruePrimaryConvMCPt,fHistoGammaMCConvPt,1,1,"B");

        TH1D* fHistoGammaTruePrimaryConvMCPtOrBinTemp               = (TH1D*)fHistoGammaTruePrimaryConvMCPtOrBin->Clone("fHistoGammaTruePrimaryConvMCPtOrBinTemp");
        fHistoGammaTruePrimaryConvMCPtOrBinTemp->Scale(1./fHistoGammaTruePrimaryConvMCPtOrBinTemp->GetBinWidth(5));
        fHistoGammaMCPrimaryRecoEffMCPtOrBin                        = (TH1D*)fHistoGammaTruePrimaryConvMCPtOrBinTemp->Clone("GammaPrimaryRecoEff_Pt_OriginalBinning");
        fHistoGammaMCPrimaryRecoEffMCPtOrBin->Sumw2();
        fHistoGammaMCPrimaryRecoEffMCPtOrBin->Divide(fHistoGammaMCPrimaryRecoEffMCPtOrBin,fHistoGammaMCConvPtOrBin,1,1,"B");

        // secondary reconstruction efficiencies
        if(fUseCocktail && nHistogramDimension==2){
            // K0s
            fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPt           = new TH1D("SecondaryGammaFromXFromK0sRecoEff_MCPt","",fNBinsPt,fBinsPt);
            fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPt->Sumw2();
            fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPt->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromK0sMCPt,fHistoSecondaryGammaConvFromXFromK0sPt,1,1,"B");

            fHistoSecondaryGammaFromXFromK0sMCRecoEffPt             = new TH1D("SecondaryGammaFromXFromK0sRecoEff_Pt","",fNBinsPt,fBinsPt);
            fHistoSecondaryGammaFromXFromK0sMCRecoEffPt->Sumw2();
            fHistoSecondaryGammaFromXFromK0sMCRecoEffPt->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt,fHistoSecondaryGammaConvFromXFromK0sPt,1,1,"B");
            
            fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPtOrBin      = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromK0sMCPtOrBin->Clone("SecondaryGammaFromXFromK0sRecoEff_MCPtOrBin");
            fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPtOrBin->Sumw2();
            fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPtOrBin->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromK0sMCPtOrBin,fHistoSecondaryGammaConvFromXFromK0sPtOrBin,1,1,"B");
            
            fHistoSecondaryGammaFromXFromK0sMCRecoEffPtOrBin        = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromK0sPtOrBin->Clone("SecondaryGammaFromXFromK0sRecoEff_PtOrBin");
            fHistoSecondaryGammaFromXFromK0sMCRecoEffPtOrBin->Sumw2();
            fHistoSecondaryGammaFromXFromK0sMCRecoEffPtOrBin->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromK0sPtOrBin,fHistoSecondaryGammaConvFromXFromK0sPtOrBin,1,1,"B");

            // K0l
            fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPt           = new TH1D("SecondaryGammaFromXFromK0lRecoEff_MCPt","",fNBinsPt,fBinsPt);
            fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPt->Sumw2();
            fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPt->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromK0lMCPt,fHistoSecondaryGammaConvFromXFromK0lPt,1,1,"B");

            fHistoSecondaryGammaFromXFromK0lMCRecoEffPt             = new TH1D("SecondaryGammaFromXFromK0lRecoEff_Pt","",fNBinsPt,fBinsPt);
            fHistoSecondaryGammaFromXFromK0lMCRecoEffPt->Sumw2();
            fHistoSecondaryGammaFromXFromK0lMCRecoEffPt->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromK0lPt,fHistoSecondaryGammaConvFromXFromK0lPt,1,1,"B");

            fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPtOrBin     = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromK0lMCPtOrBin->Clone("SecondaryGammaFromXFromK0lRecoEff_MCPtOrBin");
            fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPtOrBin->Sumw2();
            fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPtOrBin->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromK0lMCPtOrBin,fHistoSecondaryGammaConvFromXFromK0lPtOrBin,1,1,"B");

            fHistoSecondaryGammaFromXFromK0lMCRecoEffPtOrBin        = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromK0lPtOrBin->Clone("SecondaryGammaFromXFromK0lRecoEff_PtOrBin");
            fHistoSecondaryGammaFromXFromK0lMCRecoEffPtOrBin->Sumw2();
            fHistoSecondaryGammaFromXFromK0lMCRecoEffPtOrBin->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromK0lPtOrBin,fHistoSecondaryGammaConvFromXFromK0lPtOrBin,1,1,"B");

            // Lambda
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPt        = new TH1D("SecondaryGammaFromXFromLambdaRecoEff_MCPt","",fNBinsPt,fBinsPt);
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPt->Sumw2();
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPt->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromLambdaMCPt,fHistoSecondaryGammaConvFromXFromLambdaPt,1,1,"B");
            
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffPt          = new TH1D("SecondaryGammaFromXFromLambdaRecoEff_Pt","",fNBinsPt,fBinsPt);
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffPt->Sumw2();
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffPt->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt,fHistoSecondaryGammaConvFromXFromLambdaPt,1,1,"B");
            
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPtOrBin   = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromLambdaMCPtOrBin->Clone("SecondaryGammaFromXFromLambdaRecoEff_MCPtOrBin");
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPtOrBin->Sumw2();
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPtOrBin->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromLambdaMCPtOrBin,fHistoSecondaryGammaConvFromXFromLambdaPtOrBin,1,1,"B");
            
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffPtOrBin     = (TH1D*)fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPtOrBin->Clone("SecondaryGammaFromXFromLambdaRecoEff_PtOrBin");
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffPtOrBin->Sumw2();
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffPtOrBin->Divide(fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPtOrBin,fHistoSecondaryGammaConvFromXFromLambdaPtOrBin,1,1,"B");
        }
        // ==========================================

        // ========== identified MC BG ==============
        fHistoGammaMCBackground                                     = new TH1D("MCrec_Background","",fNBinsPt,fBinsPt);
        fHistoGammaMCBackground->Sumw2();
        fHistoGammaMCBackground                                     = (TH1D*)fHistoGammaMCrecConvPt->Clone("MCrec_Background");
        fHistoGammaMCBackground->Add(fHistoGammaTrueConvPt,-1);
        // ==========================================
    } 
    
    if (fEnableCalo && fEnablePCM) {
        TAxis *xAxis                                                = fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetXaxis();
        TAxis *yAxis                                                = fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetYaxis();

        // =========== Response matrix ==============
        fHistoGammaTruePrimaryCalo_recPt_MCPt_MC_Rebin              = new TH2D("TruePrimaryCaloGamma_recPt_MCPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
        for(Int_t x = 1; x<fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetNbinsX()+1; x++){
            for(Int_t y = 1; y<fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetNbinsY()+1; y++){
                Double_t binContent                                 = fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetBinContent(x,y);
                Double_t xcenter                                    = xAxis->GetBinCenter(x);
                Double_t ycenter                                    = yAxis->GetBinCenter(y);
                fHistoGammaTruePrimaryCalo_recPt_MCPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
            }
        }
        // ==========================================
        
        // ======== Secondary fractions =============
        fHistoFracAllGammaCaloToSecOrBin                            = (TH1D*) fHistoGammaCaloPtOrBin->Clone("FracAllGammaCaloToSecOriginalBinning");
        fHistoFracAllGammaCaloToSecOrBin->Divide(fHistoGammaTrueSecondaryCaloPtOrBin,fHistoFracAllGammaCaloToSecOrBin,1,1,"B");

        fHistoFracAllGammaCaloToSec                                 = (TH1D*) fHistoGammaCaloPt->Clone("FracAllGammaCaloToSec");
        fHistoFracAllGammaCaloToSec->Divide(fHistoGammaTrueSecondaryCaloPt,fHistoFracAllGammaCaloToSec,1,1,"B");

        fHistoFracAllGammaCaloToSecFromXFromK0sOrBin                     = (TH1D*) fHistoGammaCaloPtOrBin->Clone("FracAllGammaCaloToSecFromK0sOriginalBinning");
        fHistoFracAllGammaCaloToSecFromXFromK0sOrBin->Divide(fHistoGammaTrueSecondaryCaloFromXFromK0sPtOrBin,fHistoFracAllGammaCaloToSecFromXFromK0sOrBin,1,1,"B");

        fHistoFracAllGammaCaloToSecFromXFromK0s                          = (TH1D*) fHistoGammaCaloPt->Clone("fHistoFracAllGammaCaloToSecFromXFromK0s");
        fHistoFracAllGammaCaloToSecFromXFromK0s->Divide(fHistoGammaTrueSecondaryCaloFromXFromK0sPt,fHistoFracAllGammaCaloToSecFromXFromK0s,1,1,"B");

        fHistoFracAllGammaCaloToSecFromXFromLambdaOrBin                  = (TH1D*) fHistoGammaCaloPtOrBin->Clone("FracAllGammaCaloToSecFromLambdaOriginalBinning");
        fHistoFracAllGammaCaloToSecFromXFromLambdaOrBin->Divide(fHistoGammaTrueSecondaryCaloFromXFromLambdaPtOrBin,fHistoFracAllGammaCaloToSecFromXFromLambdaOrBin,1,1,"B");

        fHistoFracAllGammaCaloToSecFromXFromLambda                       = (TH1D*) fHistoGammaCaloPt->Clone("fHistoFracAllGammaCaloToSecFromXFromLambda");
        fHistoFracAllGammaCaloToSecFromXFromLambda->Divide(fHistoGammaTrueSecondaryCaloFromXFromLambdaPt,fHistoFracAllGammaCaloToSecFromXFromLambda,1,1,"B");
        // ==========================================
        
         // ================= PURITY =================
        fHistoGammaCaloMCPurity                                     = new TH1D("GammaCaloPurity_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaCaloMCPurity->Sumw2();
        fHistoGammaCaloMCPurity->Divide(fHistoGammaTrueCaloPt,fHistoGammaCaloPt,1,1,"B");

        fHistoGammaMCrecPrimaryCaloPt                               = (TH1D*) fHistoGammaCaloPt->Clone("MCrec_PrimaryCaloGamma_Pt");
        fHistoGammaMCrecPrimaryCaloPt->Add(fHistoGammaTrueSecondaryCaloPt,-1);
        fHistoGammaCaloMCTruePurity                                 = new TH1D("GammaCaloTruePurity_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaCaloMCTruePurity->Sumw2();
        fHistoGammaCaloMCTruePurity->Divide(fHistoGammaTruePrimaryCaloPt,fHistoGammaMCrecPrimaryCaloPt,1,1,"B");

        fHistoGammaMCrecPrimaryCaloPtOrBin                          = (TH1D*) fHistoGammaCaloPtOrBin->Clone("MC_ESDPrimaryCaloGammaPt");
        fHistoGammaMCrecPrimaryCaloPtOrBin->Add(fHistoGammaTrueSecondaryCaloPtOrBin,-1);
        fHistoGammaCaloMCTruePurityOrBin                            = (TH1D*)fHistoGammaMCrecPrimaryCaloPtOrBin->Clone("GammaCaloTruePurity_OriginalBinning_Pt");
        fHistoGammaCaloMCTruePurityOrBin->Sumw2();
        fHistoGammaCaloMCTruePurityOrBin->Divide(fHistoGammaTruePrimaryCaloPtOrBin,fHistoGammaMCrecPrimaryCaloPtOrBin,1,1,"B");
        // ==========================================
       
        // ================ Reco Eff ================
        fHistoGammaCaloMCRecoEff                                    = new TH1D("GammaCaloRecoEff_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaCaloMCRecoEff->Sumw2();
        fHistoGammaCaloMCRecoEff->Divide(fHistoGammaTrueCaloPt,fHistoGammaMCAllInEMCAccPt,1,1,"B");

        fHistoGammaCaloMCPrimaryRecoEff                             = new TH1D("GammaCaloPrimaryRecoEff_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaCaloMCPrimaryRecoEff->Sumw2();
        fHistoGammaCaloMCPrimaryRecoEff->Divide(fHistoGammaTruePrimaryCaloPt,fHistoGammaMCAllInEMCAccPt,1,1,"B");

        fHistoGammaCaloMCPrimaryRecoEffMCPt                         = new TH1D("GammaCaloPrimaryRecoEff_MCPt","",fNBinsPt,fBinsPt);
        fHistoGammaCaloMCPrimaryRecoEffMCPt->Sumw2();
        fHistoGammaCaloMCPrimaryRecoEffMCPt->Divide(fHistoGammaTruePrimaryCaloMCPt,fHistoGammaMCAllInEMCAccPt,1,1,"B");
        // ==========================================

        // ========== identified MC BG ==============
        fHistoGammaCaloMCBackground                                 = new TH1D("MCrec_Calo_Background","",fNBinsPt,fBinsPt);
        fHistoGammaCaloMCBackground->Sumw2();
        fHistoGammaCaloMCBackground = (TH1D*)fHistoGammaMCrecCaloPt->Clone("MCrec_Calo_Background");
        fHistoGammaCaloMCBackground->Add(fHistoGammaTrueCaloPt,-1);
        // ==========================================
    }
    
    if(fEnableCalo && !fEnablePCM){

        // =========== Response matrix ==============
        TAxis *xAxis                                                = fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetXaxis();
        TAxis *yAxis                                                = fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetYaxis();
        fHistoGammaTruePrimaryCalo_recPt_MCPt_MC_Rebin              = new TH2D("TruePrimaryCaloGamma_recPt_MCPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
        for(Int_t x = 1; x<fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetNbinsX()+1; x++){
            for(Int_t y = 1; y<fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetNbinsY()+1; y++){
                Double_t binContent                                 = fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->GetBinContent(x,y);
                Double_t xcenter                                    = xAxis->GetBinCenter(x);
                Double_t ycenter                                    = yAxis->GetBinCenter(y);
                fHistoGammaTruePrimaryCalo_recPt_MCPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
            }
        }
        
        xAxis                                                       = fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->GetXaxis();
        yAxis                                                       = fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->GetYaxis();
        fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC_Rebin          = new TH2D("TruePrimaryCaloConvGamma_recPt_MCPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
        for(Int_t x = 1; x<fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->GetNbinsX()+1; x++){
            for(Int_t y = 1; y<fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->GetNbinsY()+1; y++){
                Double_t binContent                                 = fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC->GetBinContent(x,y);
                Double_t xcenter                                    = xAxis->GetBinCenter(x);
                Double_t ycenter                                    = yAxis->GetBinCenter(y);
                fHistoGammaTruePrimaryCaloConv_recPt_MCPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
            }
        }

        xAxis                                                       = fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->GetXaxis();
        yAxis                                                       = fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->GetYaxis();
        fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC_Rebin        = new TH2D("TruePrimaryCaloUnConvGamma_recPt_MCPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
        for(Int_t x = 1; x<fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->GetNbinsX()+1; x++){
            for(Int_t y = 1; y<fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->GetNbinsY()+1; y++){
                Double_t binContent                                 = fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC->GetBinContent(x,y);
                Double_t xcenter                                    = xAxis->GetBinCenter(x);
                Double_t ycenter                                    = yAxis->GetBinCenter(y);
                fHistoGammaTruePrimaryCaloUnConv_recPt_MCPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
            }
        }
        
        // secondary response matrices
        if(fUseCocktail && nHistogramDimension==2){
            xAxis                                                   = fHistoGammaTrueSecondaryFromXFromK0sCalo_MCPt_recPt_MC->GetXaxis();
            yAxis                                                   = fHistoGammaTrueSecondaryFromXFromK0sCalo_MCPt_recPt_MC->GetYaxis();

            fHistoGammaTrueSecondaryFromXFromK0sCalo_MCPt_recPt_MC_Rebin = new TH2D("ESD_TrueSecondaryClusGammaFromXFromK0s_MCPt_recPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
            for(Int_t x = 1; x<fHistoGammaTrueSecondaryFromXFromK0sCalo_MCPt_recPt_MC->GetNbinsX()+1; x++){
                for(Int_t y = 1; y<fHistoGammaTrueSecondaryFromXFromK0sCalo_MCPt_recPt_MC->GetNbinsY()+1; y++){
                    Double_t binContent                             = fHistoGammaTrueSecondaryFromXFromK0sCalo_MCPt_recPt_MC->GetBinContent(x,y);
                    Double_t xcenter                                = xAxis->GetBinCenter(x);
                    Double_t ycenter                                = yAxis->GetBinCenter(y);
                    fHistoGammaTrueSecondaryFromXFromK0sCalo_MCPt_recPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
                }
            }
            
            fHistoGammaTrueSecondaryFromXFromK0lCalo_MCPt_recPt_MC_Rebin = new TH2D("ESD_TrueSecondaryClusGammaFromXFromK0l_MCPt_recPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
            for(Int_t x = 1; x<fHistoGammaTrueSecondaryFromXFromK0lCalo_MCPt_recPt_MC->GetNbinsX()+1; x++){
                for(Int_t y = 1; y<fHistoGammaTrueSecondaryFromXFromK0lCalo_MCPt_recPt_MC->GetNbinsY()+1; y++){
                    Double_t binContent                             = fHistoGammaTrueSecondaryFromXFromK0lCalo_MCPt_recPt_MC->GetBinContent(x,y);
                    Double_t xcenter                                = xAxis->GetBinCenter(x);
                    Double_t ycenter                                = yAxis->GetBinCenter(y);
                    fHistoGammaTrueSecondaryFromXFromK0lCalo_MCPt_recPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
                }
            }

            fHistoGammaTrueSecondaryFromXFromLambdaCalo_MCPt_recPt_MC_Rebin = new TH2D("ESD_TrueSecondaryClusGammaFromXFromLambda_MCPt_recPt_Rebin","",fNBinsPt,fBinsPt,fNBinsPt,fBinsPt);
            for(Int_t x = 1; x<fHistoGammaTrueSecondaryFromXFromLambdaCalo_MCPt_recPt_MC->GetNbinsX()+1; x++){
                for(Int_t y = 1; y<fHistoGammaTrueSecondaryFromXFromLambdaCalo_MCPt_recPt_MC->GetNbinsY()+1; y++){
                    Double_t binContent                             = fHistoGammaTrueSecondaryFromXFromLambdaCalo_MCPt_recPt_MC->GetBinContent(x,y);
                    Double_t xcenter                                = xAxis->GetBinCenter(x);
                    Double_t ycenter                                = yAxis->GetBinCenter(y);
                    fHistoGammaTrueSecondaryFromXFromLambdaCalo_MCPt_recPt_MC_Rebin->Fill(xcenter,ycenter,binContent);
                }
            }
        }
        // ==========================================

        // ======== Secondary fractions =============
        fHistoFracAllGammaToSec                                     = (TH1D*)fHistoGammaCaloPt->Clone("FracAllGammaToSec");
        fHistoFracAllGammaToSec->Divide(fHistoGammaTrueSecondaryCaloPt,fHistoFracAllGammaToSec,1,1,"B");

        fHistoFracAllGammaToSecOrBin                                = (TH1D*)fHistoGammaCaloPtOrBin->Clone("FracAllGammaToSecOriginalBinning");
        fHistoFracAllGammaToSecOrBin->Divide(fHistoGammaTrueSecondaryCaloPtOrBin,fHistoFracAllGammaToSecOrBin,1,1,"B");

        fHistoFracAllGammaToSecFromXFromK0s                         = (TH1D*)fHistoGammaCaloPt->Clone("FracAllGammaToSecFromXFromK0s");
        fHistoFracAllGammaToSecFromXFromK0s->Divide(fHistoGammaTrueSecondaryCaloFromXFromK0sPt,fHistoFracAllGammaToSecFromXFromK0s,1,1,"B");

        fHistoFracAllGammaToSecFromXFromK0sOrBin                    = (TH1D*)fHistoGammaCaloPtOrBin->Clone("FracAllGammaToSecFromXFromK0sOriginalBinning");
        fHistoFracAllGammaToSecFromXFromK0sOrBin->Divide(fHistoGammaTrueSecondaryCaloFromXFromK0sPtOrBin,fHistoFracAllGammaToSecFromXFromK0sOrBin,1,1,"B");

        if (nHistogramDimension==2) {
            fHistoFracAllGammaToSecFromXFromK0l                     = (TH1D*)fHistoGammaCaloPt->Clone("FracAllGammaToSecFromXFromK0l");
            fHistoFracAllGammaToSecFromXFromK0l->Divide(fHistoGammaTrueSecondaryCaloFromXFromK0lPt,fHistoFracAllGammaToSecFromXFromK0l,1,1,"B");

            fHistoFracAllGammaToSecFromXFromK0lOrBin                = (TH1D*)fHistoGammaCaloPtOrBin->Clone("FracAllGammaToSecFromXFromK0lOriginalBinning");
            fHistoFracAllGammaToSecFromXFromK0lOrBin->Divide(fHistoGammaTrueSecondaryCaloFromXFromK0lPtOrBin,fHistoFracAllGammaToSecFromXFromK0lOrBin,1,1,"B");
        }

        fHistoFracAllGammaToSecFromXFromLambda                      = (TH1D*)fHistoGammaCaloPt->Clone("FracAllGammaToSecFromXFromLambda");
        fHistoFracAllGammaToSecFromXFromLambda->Divide(fHistoGammaTrueSecondaryCaloFromXFromLambdaPt,fHistoFracAllGammaToSecFromXFromLambda,1,1,"B");

        fHistoFracAllGammaToSecFromXFromLambdaOrBin                 = (TH1D*)fHistoGammaCaloPtOrBin->Clone("FracAllGammaToSecFromXFromLambdaOriginalBinning");
        fHistoFracAllGammaToSecFromXFromLambdaOrBin->Divide(fHistoGammaTrueSecondaryCaloFromXFromLambdaPtOrBin,fHistoFracAllGammaToSecFromXFromLambdaOrBin,1,1,"B");

        if(nHistogramDimension==2 && fUseCocktail){
            fHistoFracAllGammaToSecRest                             = (TH1D*)fHistoGammaCaloPt->Clone("FracAllGammaToSecRest");
            fHistoFracAllGammaToSecRest->Divide(fHistoGammaTrueSecondaryCaloRestPt,fHistoFracAllGammaToSecRest,1,1,"B");
            
            fHistoFracAllGammaToSecRestOrBin                        = (TH1D*)fHistoGammaCaloPtOrBin->Clone("FracAllGammaToSecRestOriginalBinning");
            fHistoFracAllGammaToSecRestOrBin->Divide(fHistoGammaTrueSecondaryCaloRestPtOrBin,fHistoFracAllGammaToSecRestOrBin,1,1,"B");
        }
        // ==========================================

        // ================= PURITY =================
        fHistoGammaMCPurity                                         = new TH1D("GammaPurity_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCPurity->Sumw2();
        fHistoGammaMCPurity->Divide(fHistoGammaTrueCaloPt,fHistoGammaCaloPt,1,1,"B");

        fHistoGammaMCrecPrimaryCaloPt                               = (TH1D*)fHistoGammaCaloPt->Clone("MCrec_PrimaryCaloGamma_Pt");
        fHistoGammaMCrecPrimaryCaloPt->Add(fHistoGammaTrueSecondaryCaloPt,-1);
        fHistoGammaMCTruePurity                                     = new TH1D("GammaTruePurity_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCTruePurity->Sumw2();
        fHistoGammaMCTruePurity->Divide(fHistoGammaTruePrimaryCaloPt,fHistoGammaMCrecPrimaryCaloPt,1,1,"B");

        fHistoGammaMCrecPrimaryCaloPtOrBin                          = (TH1D*)fHistoGammaCaloPtOrBin->Clone("MC_ESDPrimaryCaloGammaPt");
        fHistoGammaMCrecPrimaryCaloPtOrBin->Add(fHistoGammaTrueSecondaryCaloPtOrBin,-1);
        fHistoGammaMCTruePurityOrBin                                = (TH1D*)fHistoGammaMCrecPrimaryCaloPtOrBin->Clone("GammaTruePurity_OriginalBinning_Pt");
        fHistoGammaMCTruePurityOrBin->Sumw2();
        fHistoGammaMCTruePurityOrBin->Divide(fHistoGammaTruePrimaryCaloPtOrBin,fHistoGammaMCrecPrimaryCaloPtOrBin,1,1,"B");
        // ==========================================

        // ================ Reco Eff ================
        fHistoGammaMCRecoEff                                        = new TH1D("GammaRecoEff_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCRecoEff->Sumw2();
        fHistoGammaMCRecoEff->Divide(fHistoGammaTrueCaloPt,fHistoGammaMCAllPt,1,1,"B");

        fHistoGammaMCPrimaryRecoEff                                 = new TH1D("GammaPrimaryRecoEff_Pt","",fNBinsPt,fBinsPt);
        fHistoGammaMCPrimaryRecoEff->Sumw2();
        fHistoGammaMCPrimaryRecoEff->Divide(fHistoGammaTruePrimaryCaloPt,fHistoGammaMCAllPt,1,1,"B");

        TH1D* fHistoGammaTruePrimaryCaloPtOrBinTemp                 = (TH1D*)fHistoGammaTruePrimaryCaloPtOrBin->Clone("fHistoGammaTruePrimaryCaloPtOrBinTemp");
        fHistoGammaTruePrimaryCaloPtOrBinTemp->Scale(1./fHistoGammaTruePrimaryCaloPtOrBinTemp->GetBinWidth(5));
        fHistoGammaMCPrimaryRecoEffOrBin                            = (TH1D*)fHistoGammaTruePrimaryCaloPtOrBinTemp->Clone("GammaPrimaryRecoEff_Pt_OriginalBinning");
        fHistoGammaMCPrimaryRecoEffOrBin->Sumw2();
        fHistoGammaMCPrimaryRecoEffOrBin->Divide(fHistoGammaMCPrimaryRecoEffOrBin,fHistoGammaMCAllPtOrBin,1,1,"B");

        fHistoGammaMCPrimaryRecoEffMCPt                             = new TH1D("GammaPrimaryRecoEff_MCPt","",fNBinsPt,fBinsPt);
        fHistoGammaMCPrimaryRecoEffMCPt->Sumw2();
        fHistoGammaMCPrimaryRecoEffMCPt->Divide(fHistoGammaTruePrimaryCaloMCPt,fHistoGammaMCAllPt,1,1,"B");

        TH1D* fHistoGammaTruePrimaryCaloMCPtOrBinTemp               = (TH1D*)fHistoGammaTruePrimaryCaloMCPtOrBin->Clone("fHistoGammaTruePrimaryCaloMCPtOrBinTemp");
        fHistoGammaTruePrimaryCaloMCPtOrBinTemp->Scale(1./fHistoGammaTruePrimaryCaloMCPtOrBinTemp->GetBinWidth(5));
        fHistoGammaMCPrimaryRecoEffMCPtOrBin                        = (TH1D*)fHistoGammaTruePrimaryCaloMCPtOrBinTemp->Clone("GammaPrimaryRecoEff_Pt_OriginalBinning");
        fHistoGammaMCPrimaryRecoEffMCPtOrBin->Sumw2();
        fHistoGammaMCPrimaryRecoEffMCPtOrBin->Divide(fHistoGammaMCPrimaryRecoEffMCPtOrBin,fHistoGammaMCAllPtOrBin,1,1,"B");

        // secondary reco eff
        if(fUseCocktail && nHistogramDimension==2){
            
            // K0s
            fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPt           = new TH1D("SecondaryGammaFromXFromK0sRecoEff_MCPt","",fNBinsPt,fBinsPt);
            fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPt->Sumw2();
            fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPt->Divide(fHistoGammaTrueSecondaryCaloFromXFromK0sMCPt,fHistoAllSecondaryGammaFromXFromK0sPt,1,1,"B");

            fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPtOrBin      = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromK0sMCPtOrBin->Clone("SecondaryGammaFromXFromK0sRecoEff_MCPtOrBin");
            fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPtOrBin->Sumw2();
            fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPtOrBin->Divide(fHistoGammaTrueSecondaryCaloFromXFromK0sMCPtOrBin,fHistoAllSecondaryGammaFromXFromK0sPtOrBin,1,1,"B");
            
            fHistoSecondaryGammaFromXFromK0sMCRecoEffPt             = new TH1D("SecondaryGammaFromXFromK0sRecoEff_Pt","",fNBinsPt,fBinsPt);
            fHistoSecondaryGammaFromXFromK0sMCRecoEffPt->Sumw2();
            fHistoSecondaryGammaFromXFromK0sMCRecoEffPt->Divide(fHistoGammaTrueSecondaryCaloFromXFromK0sPt,fHistoAllSecondaryGammaFromXFromK0sPt,1,1,"B");
            
            fHistoSecondaryGammaFromXFromK0sMCRecoEffPtOrBin        = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromK0sPtOrBin->Clone("SecondaryGammaFromXFromK0sRecoEff_PtOrBin");
            fHistoSecondaryGammaFromXFromK0sMCRecoEffPtOrBin->Sumw2();
            fHistoSecondaryGammaFromXFromK0sMCRecoEffPtOrBin->Divide(fHistoGammaTrueSecondaryCaloFromXFromK0sPtOrBin,fHistoAllSecondaryGammaFromXFromK0sPtOrBin,1,1,"B");

            // K0l
            fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPt           = new TH1D("SecondaryGammaFromXFromK0lRecoEff_MCPt","",fNBinsPt,fBinsPt);
            fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPt->Sumw2();
            fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPt->Divide(fHistoGammaTrueSecondaryCaloFromXFromK0lMCPt,fHistoAllSecondaryGammaFromXFromK0lPt,1,1,"B");
            
            fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPtOrBin      = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromK0lMCPtOrBin->Clone("SecondaryGammaFromXFromK0lRecoEff_MCPtOrBin");
            fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPtOrBin->Sumw2();
            fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPtOrBin->Divide(fHistoGammaTrueSecondaryCaloFromXFromK0lMCPtOrBin,fHistoAllSecondaryGammaFromXFromK0lPtOrBin,1,1,"B");

            fHistoSecondaryGammaFromXFromK0lMCRecoEffPt             = new TH1D("SecondaryGammaFromXFromK0lRecoEff_Pt","",fNBinsPt,fBinsPt);
            fHistoSecondaryGammaFromXFromK0lMCRecoEffPt->Sumw2();
            fHistoSecondaryGammaFromXFromK0lMCRecoEffPt->Divide(fHistoGammaTrueSecondaryCaloFromXFromK0lPt,fHistoAllSecondaryGammaFromXFromK0lPt,1,1,"B");
            
            fHistoSecondaryGammaFromXFromK0lMCRecoEffPtOrBin        = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromK0lPtOrBin->Clone("SecondaryGammaFromXFromK0lRecoEff_PtOrBin");
            fHistoSecondaryGammaFromXFromK0lMCRecoEffPtOrBin->Sumw2();
            fHistoSecondaryGammaFromXFromK0lMCRecoEffPtOrBin->Divide(fHistoGammaTrueSecondaryCaloFromXFromK0lPtOrBin,fHistoAllSecondaryGammaFromXFromK0lPtOrBin,1,1,"B");

            // Lambda
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPt        = new TH1D("SecondaryGammaFromXFromLambdaRecoEff_MCPt","",fNBinsPt,fBinsPt);
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPt->Sumw2();
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPt->Divide(fHistoGammaTrueSecondaryCaloFromXFromLambdaMCPt,fHistoAllSecondaryGammaFromXFromLambdaPt,1,1,"B");
            
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPtOrBin   = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromLambdaMCPtOrBin->Clone("SecondaryGammaFromXFromLambdaRecoEff_MCPtOrBin");
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPtOrBin->Sumw2();
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPtOrBin->Divide(fHistoGammaTrueSecondaryCaloFromXFromLambdaMCPtOrBin,fHistoAllSecondaryGammaFromXFromLambdaPtOrBin,1,1,"B");

            fHistoSecondaryGammaFromXFromLambdaMCRecoEffPt          = new TH1D("SecondaryGammaFromXFromLambdaRecoEff_Pt","",fNBinsPt,fBinsPt);
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffPt->Sumw2();
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffPt->Divide(fHistoGammaTrueSecondaryCaloFromXFromLambdaPt,fHistoAllSecondaryGammaFromXFromLambdaPt,1,1,"B");
            
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffPtOrBin     = (TH1D*)fHistoGammaTrueSecondaryCaloFromXFromLambdaPtOrBin->Clone("SecondaryGammaFromXFromLambdaRecoEff_PtOrBin");
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffPtOrBin->Sumw2();
            fHistoSecondaryGammaFromXFromLambdaMCRecoEffPtOrBin->Divide(fHistoGammaTrueSecondaryCaloFromXFromLambdaPtOrBin,fHistoAllSecondaryGammaFromXFromLambdaPtOrBin,1,1,"B");
        }
        // ==========================================

        // ========== identified MC BG ==============
        fHistoGammaMCBackground                                     = new TH1D("MCrec_Background","",fNBinsPt,fBinsPt);
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

    InitializeBinning(setPi0, numberOfBins, energy, fDirectPhoton, fMode, fEventCutSelection, fClusterCutSelection);

    
    fDeltaPt                                            = new TH1D("deltaPt","",fNBinsPt,fBinsPt);
    for(Int_t iPt=fStartPtBin+1;iPt<fNBinsPt+1;iPt++){
        fDeltaPt->SetBinContent(iPt,fBinsPt[iPt]-fBinsPt[iPt-1]);
        fDeltaPt->SetBinError(iPt,0);
    }

    // initializing binning used for the DCAz distributions
    if(!addSig && !(mode == 4 || mode == 5)){
        if (fBinsPtDCAzDist && fNBinsPtDCAzDist) {
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
    } else if ((fEnergyFlag.CompareTo("7TeV") == 0) && (fMeson.CompareTo("Pi0") == 0) && (fDirectPhoton.CompareTo("directPhoton") == 0)) {
        nIterationsShowBackground[0]                    = 12;
        nIterationsShowBackground[1]                    = 12;
        nIterationsShowBackground[2]                    = 19;
        nIterationsShowBackground[3]                    = 20;
        optionShowBackground[0]                         = "BackDecreasingWindow, BackSmoothing3";   // standard
        optionShowBackground[1]                         = "nosmoothing";
        optionShowBackground[2]                         = "BackDecreasingWindow, BackSmoothing7";
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

    if (mode == 2 || mode == 3){
        InitializeWindows(setPi0, fMode, fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2) );
        
        // initialize mass histo array
        fHistoGGInvMassPtGConvBin                           = new TH1D*[fNBinsPt];
        fHistoSignalInvMassPtGConvBin                       = new TH1D*[fNBinsPt];
        fHistoSignalInvMassLeftPtGConvBin                   = new TH1D*[fNBinsPt];
        fHistoBackInvMassPtGconvBin                         = new TH1D*[fNBinsPt];
        fHistoBackNormInvMassPtGconvBin                     = new TH1D*[fNBinsPt];
        fHistoBackNormInvMassLeftPtGconvBin                 = new TH1D*[fNBinsPt];
        fFitSignalInvMassPtBin                              = new TF1*[fNBinsPt];
        fFitInvMassLeftPtBin                                = new TF1*[fNBinsPt];
        fFitBckInvMassPtBin                                 = new TF1*[fNBinsPt];
        fFitBckInvMassLeftPtBin                             = new TF1*[fNBinsPt];
        for(Int_t i = 0;i<fNBinsPt; i++){
            fHistoGGInvMassPtGConvBin[i]                    = NULL;
            fHistoSignalInvMassPtGConvBin[i]                = NULL;
            fHistoSignalInvMassLeftPtGConvBin[i]            = NULL;
            fHistoBackInvMassPtGconvBin[i]                  = NULL;
            fHistoBackNormInvMassPtGconvBin[i]              = NULL;
            fHistoBackNormInvMassLeftPtGconvBin[i]          = NULL;
            fFitSignalInvMassPtBin[i]                       = NULL;
            fFitInvMassLeftPtBin[i]                         = NULL;
            fFitBckInvMassPtBin[i]                          = NULL;
            fFitBckInvMassLeftPtBin[i]                      = NULL;
        }
        
        fMesonMass                                                      = new Double_t[fNBinsPt];
        fMesonMassError                                                 = new Double_t[fNBinsPt];
        fMesonFWHM                                                      = new Double_t[fNBinsPt];
        fMesonFWHMError                                                 = new Double_t[fNBinsPt];
        fMesonMassLeft                                                  = new Double_t[fNBinsPt];
        fMesonFWHMLeft                                                  = new Double_t[fNBinsPt];
        fMesonMassLeftError                                             = new Double_t[fNBinsPt];
        fMesonFWHMLeftError                                             = new Double_t[fNBinsPt];
        fMesonLambdaTailpar                                             = new Double_t[fNBinsPt];
        fMesonLambdaTailparError                                        = new Double_t[fNBinsPt];
        fMesonLambdaTailMCpar                                           = new Double_t[fNBinsPt];
        fMesonLambdaTailMCparError                                      = new Double_t[fNBinsPt];

        // chi2 test for different function 
        for (Int_t m = 0; m < 4; m++){
            fMesonChi2[m]                                               = new Double_t[fNBinsPt];
        }
        
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
            fMesonYieldsResidualBckFunc[k]                              = new Double_t[fNBinsPt];
            fMesonYieldsCorResidualBckFunc[k]                           = new Double_t[fNBinsPt];
            fMesonYieldsPerEvent[k]                                     = new Double_t[fNBinsPt];
            
            // Initialize error arrays
            fGGYieldsError[k]                                           = new Double_t[fNBinsPt];
            fBckYieldsError[k]                                          = new Double_t[fNBinsPt];
            fTotalBckYieldsError[k]                                     = new Double_t[fNBinsPt];
            fMesonYieldsError[k]                                        = new Double_t[fNBinsPt];
            fMesonYieldsResidualBckFuncError[k]                         = new Double_t[fNBinsPt];
            fMesonYieldsCorResidualBckFuncError[k]                      = new Double_t[fNBinsPt];
            fMesonYieldsPerEventError[k]                                = new Double_t[fNBinsPt];
        }
        
        for (Int_t k = 0; k < 3; k++){
            fMassWindowHigh[k]                                          = new Double_t[fNBinsPt];
            fMassWindowLow[k]                                           = new Double_t[fNBinsPt];
        }
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
    
        // write secondary gamma cocktail spectra
        if(fUseCocktail){
            if (fHistoSecondaryGammaCocktailFromXFromK0sPt)             fHistoSecondaryGammaCocktailFromXFromK0sPt->Write("CocktailSecondaryGammaFromXFromK0s_Pt",TObject::kOverwrite);
            if (fHistoSecondaryGammaCocktailFromXFromK0lPt)             fHistoSecondaryGammaCocktailFromXFromK0lPt->Write("CocktailSecondaryGammaFromXFromK0l_Pt",TObject::kOverwrite);
            if (fHistoSecondaryGammaCocktailFromXFromLambdaPt)          fHistoSecondaryGammaCocktailFromXFromLambdaPt->Write("CocktailSecondaryGammaFromXFromLambda_Pt",TObject::kOverwrite);
            if (fHistoSecondaryGammaCocktailFromXFromK0sPtOrBin)        fHistoSecondaryGammaCocktailFromXFromK0sPtOrBin->Write("CocktailSecondaryGammaFromXFromK0s_PtOrBin",TObject::kOverwrite);
            if (fHistoSecondaryGammaCocktailFromXFromK0lPtOrBin)        fHistoSecondaryGammaCocktailFromXFromK0lPtOrBin->Write("CocktailSecondaryGammaFromXFromK0l_PtOrBin",TObject::kOverwrite);
            if (fHistoSecondaryGammaCocktailFromXFromLambdaPtOrBin)     fHistoSecondaryGammaCocktailFromXFromLambdaPtOrBin->Write("CocktailSecondaryGammaFromXFromLambda_PtOrBin",TObject::kOverwrite);
        }
    
        // write basics MC quantities if possible (converted input spectrum & reconstructed validated (&primary) photons)
        if(isMC){
            if (fHistoGammaMCConvPt)                                    fHistoGammaMCConvPt->Write("MC_ConvGamma_MCPt",TObject::kOverwrite);
            if (fHistoGammaMCConvPtOrBin)                               fHistoGammaMCConvPtOrBin->Write("MC_ConvGamma_MCPt_OriginalBinning",TObject::kOverwrite);
            if (fHistoGammaTrueConvPt)                                  fHistoGammaTrueConvPt->Write("TrueConvGamma_Pt",TObject::kOverwrite);
            if (fHistoGammaTruePrimaryConvPt)                           fHistoGammaTruePrimaryConvPt->Write("TruePrimaryConvGamma_Pt",TObject::kOverwrite);
            if (fHistoGammaTrueCaloPt)                                  fHistoGammaTrueCaloPt->Write("TrueCaloGamma_Pt",TObject::kOverwrite);
            if (fHistoGammaTruePrimaryCaloPt)                           fHistoGammaTruePrimaryCaloPt->Write("TruePrimaryCaloGamma_Pt",TObject::kOverwrite);
            
            // write secondary spectr
            if(fUseCocktail && nHistogramDimension==2){
                
                // secondary conv gamma
                if (fHistoSecondaryGammaConvFromXFromK0sPt)             fHistoSecondaryGammaConvFromXFromK0sPt->Write("fHistoSecondaryGammaConvFromXFromK0sPt",TObject::kOverwrite);
                if (fHistoSecondaryGammaConvFromXFromK0lPt)             fHistoSecondaryGammaConvFromXFromK0lPt->Write("fHistoSecondaryGammaConvFromXFromK0lPt",TObject::kOverwrite);
                if (fHistoSecondaryGammaConvFromXFromLambdaPt)          fHistoSecondaryGammaConvFromXFromLambdaPt->Write("fHistoSecondaryGammaConvFromXFromLambdaPt",TObject::kOverwrite);
                if (fHistoSecondaryGammaConvFromXFromK0sPtOrBin)        fHistoSecondaryGammaConvFromXFromK0sPtOrBin->Write("fHistoSecondaryGammaConvFromXFromK0sPtOrBin",TObject::kOverwrite);
                if (fHistoSecondaryGammaConvFromXFromK0lPtOrBin)        fHistoSecondaryGammaConvFromXFromK0lPtOrBin->Write("fHistoSecondaryGammaConvFromXFromK0lPtOrBin",TObject::kOverwrite);
                if (fHistoSecondaryGammaConvFromXFromLambdaPtOrBin)     fHistoSecondaryGammaConvFromXFromLambdaPtOrBin->Write("fHistoSecondaryGammaConvFromXFromLambdaPtOrBin",TObject::kOverwrite);
                
                // all secondary gamma
                if (fHistoAllSecondaryGammaFromXFromK0sPt)              fHistoAllSecondaryGammaFromXFromK0sPt->Write("fHistoAllSecondaryGammaFromXFromK0sPt",TObject::kOverwrite);
                if (fHistoAllSecondaryGammaFromXFromK0lPt)              fHistoAllSecondaryGammaFromXFromK0lPt->Write("fHistoAllSecondaryGammaFromXFromK0lPt",TObject::kOverwrite);
                if (fHistoAllSecondaryGammaFromXFromLambdaPt)           fHistoAllSecondaryGammaFromXFromLambdaPt->Write("fHistoAllSecondaryGammaFromXFromLambdaPt",TObject::kOverwrite);
                if (fHistoAllSecondaryGammaFromXFromK0sPtOrBin)         fHistoAllSecondaryGammaFromXFromK0sPtOrBin->Write("fHistoAllSecondaryGammaFromXFromK0sPtOrBin",TObject::kOverwrite);
                if (fHistoAllSecondaryGammaFromXFromK0lPtOrBin)         fHistoAllSecondaryGammaFromXFromK0lPtOrBin->Write("fHistoAllSecondaryGammaFromXFromK0lPtOrBin",TObject::kOverwrite);
                if (fHistoAllSecondaryGammaFromXFromLambdaPtOrBin)      fHistoAllSecondaryGammaFromXFromLambdaPtOrBin->Write("fHistoAllSecondaryGammaFromXFromLambdaPtOrBin",TObject::kOverwrite);
                
                // MC validated secondary conv gamma rec Pt
                if (fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt)    fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt->Write("fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt",TObject::kOverwrite);
                if (fHistoGammaTrueSecondaryConvGammaFromXFromK0lPt)    fHistoGammaTrueSecondaryConvGammaFromXFromK0lPt->Write("fHistoGammaTrueSecondaryConvGammaFromXFromK0lPt",TObject::kOverwrite);
                if (fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt) fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt->Write("fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt",TObject::kOverwrite);
                if (fHistoGammaTrueSecondaryConvGammaFromXFromK0sPtOrBin)
                    fHistoGammaTrueSecondaryConvGammaFromXFromK0sPtOrBin->Write("fHistoGammaTrueSecondaryConvGammaFromXFromK0sPtOrBin",TObject::kOverwrite);
                if (fHistoGammaTrueSecondaryConvGammaFromXFromK0lPtOrBin)
                    fHistoGammaTrueSecondaryConvGammaFromXFromK0lPtOrBin->Write("fHistoGammaTrueSecondaryConvGammaFromXFromK0lPtOrBin",TObject::kOverwrite);
                if (fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPtOrBin)
                    fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPtOrBin->Write("fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPtOrBin",TObject::kOverwrite);
                
                // MC validated secondary conv gamma MC Pt
                if (fHistoGammaTrueSecondaryConvGammaFromXFromK0sMCPt)
                    fHistoGammaTrueSecondaryConvGammaFromXFromK0sMCPt->Write("fHistoGammaTrueSecondaryConvGammaFromXFromK0sMCPt",TObject::kOverwrite);
                if (fHistoGammaTrueSecondaryConvGammaFromXFromK0lMCPt)
                    fHistoGammaTrueSecondaryConvGammaFromXFromK0lMCPt->Write("fHistoGammaTrueSecondaryConvGammaFromXFromK0lMCPt",TObject::kOverwrite);
                if (fHistoGammaTrueSecondaryConvGammaFromXFromLambdaMCPt)
                    fHistoGammaTrueSecondaryConvGammaFromXFromLambdaMCPt->Write("fHistoGammaTrueSecondaryConvGammaFromXFromLambdaMCPt",TObject::kOverwrite);
                if (fHistoGammaTrueSecondaryConvGammaFromXFromK0sMCPtOrBin)
                    fHistoGammaTrueSecondaryConvGammaFromXFromK0sMCPtOrBin->Write("fHistoGammaTrueSecondaryConvGammaFromXFromK0sMCPtOrBin",TObject::kOverwrite);
                if (fHistoGammaTrueSecondaryConvGammaFromXFromK0lMCPtOrBin)
                    fHistoGammaTrueSecondaryConvGammaFromXFromK0lMCPtOrBin->Write("fHistoGammaTrueSecondaryConvGammaFromXFromK0lMCPtOrBin",TObject::kOverwrite);
                if (fHistoGammaTrueSecondaryConvGammaFromXFromLambdaMCPtOrBin)
                    fHistoGammaTrueSecondaryConvGammaFromXFromLambdaMCPtOrBin->Write("fHistoGammaTrueSecondaryConvGammaFromXFromLambdaMCPtOrBin",TObject::kOverwrite);
            }
        }

        // write raw photons distribution
        if (fHistoGammaCaloPt)          fHistoGammaCaloPt->Write("ESD_CaloGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaCaloPtOrBin)     fHistoGammaCaloPtOrBin->Write("ESD_CaloGamma_Pt_OriginalBinning",TObject::kOverwrite);
    
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
        if (fHistoGammaMCAllPt)                                         fHistoGammaMCAllPt->Write("MC_AllGamma_MCPt",TObject::kOverwrite);
        if (fHistoGammaMCAllPtOrBin)                                    fHistoGammaMCAllPtOrBin->Write("MC_AllGamma_OriginalBinning_MCPt",TObject::kOverwrite);
        if (fHistoGammaMCAllInEMCAccPt)                                 fHistoGammaMCAllInEMCAccPt->Write("MC_AllGammaEMCAcc_MCPt",TObject::kOverwrite);
        if (fHistoGammaMCAllInEMCAccPtOrBin)                            fHistoGammaMCAllInEMCAccPtOrBin->Write("MC_AllGammaEMCAcc_OriginalBinning_MCPt",TObject::kOverwrite);

        // write input converted photon distributions
        if (fHistoGammaMCConvPt)                                        fHistoGammaMCConvPt->Write("MC_ConvGamma_MCPt",TObject::kOverwrite);
        if (fHistoGammaMCConvPtOrBin)                                   fHistoGammaMCConvPtOrBin->Write("MC_ConvGamma_MCPt_OriginalBinning",TObject::kOverwrite);

        // write photon candidates in MC
        if (fHistoGammaMCrecConvPt)                                     fHistoGammaMCrecConvPt->Write("MCrec_ConvGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaMCrecConvPtOrBin)                                fHistoGammaMCrecConvPtOrBin->Write("MCrec_ConvGamma_OriginalBinning_Pt",TObject::kOverwrite);
        if (fHistoGammaMCrecCaloPt)                                     fHistoGammaMCrecCaloPt->Write("MCrec_CaloGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaMCrecCaloPtOrBin)                                fHistoGammaMCrecCaloPtOrBin->Write("MCrec_CaloGamma_OriginalBinning_Pt",TObject::kOverwrite);
        if (fHistoGammaMCrecPrimaryConvPt)                              fHistoGammaMCrecPrimaryConvPt->Write("MCrec_PrimaryConvGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaMCrecPrimaryCaloPt)                              fHistoGammaMCrecPrimaryCaloPt->Write("MCrec_PrimaryCaloGamma_Pt",TObject::kOverwrite);

        // write reconstructed real photons
        if (fHistoGammaTrueConvPt)                                      fHistoGammaTrueConvPt->Write("TrueConvGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueConvPtOrBin)                                 fHistoGammaTrueConvPtOrBin->Write("TrueConvGamma_Pt_OriginalBinning",TObject::kOverwrite);
        if (fHistoGammaTrueCaloPt)                                      fHistoGammaTrueCaloPt->Write("TrueCaloGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueCaloPtOrBin)                                 fHistoGammaTrueCaloPtOrBin->Write("TrueCaloGamma_Pt_OriginalBinning",TObject::kOverwrite);

        // write primary and secondary reconstructed photons
        if (fHistoGammaTruePrimaryConvPt)                               fHistoGammaTruePrimaryConvPt->Write("TruePrimaryConvGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryConvPt)                             fHistoGammaTrueSecondaryConvPt->Write("TrueSecondaryConvGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt)            fHistoGammaTrueSecondaryConvGammaFromXFromK0sPt->Write("TrueSecondaryConvGammaFromXFromK0s_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryConvGammaFromXFromK0lPt)            fHistoGammaTrueSecondaryConvGammaFromXFromK0lPt->Write("TrueSecondaryConvGammaFromXFromK0l_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt)         fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPt->Write("TrueSecondaryConvGammaFromXFromLambda_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryConvPtOrBin)                        fHistoGammaTrueSecondaryConvPtOrBin->Write("TrueSecondaryConvGamma_Pt_OriginalBinning",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryConvGammaFromXFromK0sPtOrBin)       fHistoGammaTrueSecondaryConvGammaFromXFromK0sPtOrBin->Write("TrueSecondaryConvGammaFromXFromK0s_Pt_OriginalBinning",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPtOrBin)    fHistoGammaTrueSecondaryConvGammaFromXFromLambdaPtOrBin->Write("TrueSecondaryConvGammaFromXFromLambda_Pt_OriginalBinning",TObject::kOverwrite);
        
        if (fHistoGammaTruePrimaryCaloPt)                               fHistoGammaTruePrimaryCaloPt->Write("TruePrimaryCaloGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaTruePrimaryCaloPtOrBin)                          fHistoGammaTruePrimaryCaloPtOrBin->Write("TruePrimaryCaloGamma_Pt_OriginalBinning",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryCaloPt)                             fHistoGammaTrueSecondaryCaloPt->Write("TrueSecondaryCaloGamma_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryCaloPtOrBin)                        fHistoGammaTrueSecondaryCaloPtOrBin->Write("TrueSecondaryCaloGamma_Pt_OriginalBinning",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryCaloFromXFromK0sPt)                 fHistoGammaTrueSecondaryCaloFromXFromK0sPt->Write("TrueSecondaryCaloGammaFromXFromK0s_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryCaloFromXFromK0sPtOrBin)            fHistoGammaTrueSecondaryCaloFromXFromK0sPtOrBin->Write("TrueSecondaryCaloGammaFromXFromK0s_Pt_OriginalBinning",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryCaloFromXFromLambdaPt)              fHistoGammaTrueSecondaryCaloFromXFromLambdaPt->Write("TrueSecondaryCaloGammaFromXFromLambda_Pt",TObject::kOverwrite);
        if (fHistoGammaTrueSecondaryCaloFromXFromLambdaPtOrBin)         fHistoGammaTrueSecondaryCaloFromXFromLambdaPtOrBin->Write("TrueSecondaryCaloGammaFromXFromLambda_Pt_OriginalBinning",TObject::kOverwrite);
        
        // write fractions of secondary photons
        if (fHistoFracAllGammaToSec)                                    fHistoFracAllGammaToSec->Write("FracAllGammaToSec",TObject::kOverwrite);
        if (fHistoFracAllGammaCaloToSec)                                fHistoFracAllGammaCaloToSec->Write("FracAllGammaCaloToSec",TObject::kOverwrite);
        if (fHistoFracAllGammaToSecFromXFromK0s)                        fHistoFracAllGammaToSecFromXFromK0s->Write("FracAllGammaToSecFromXFromK0s",TObject::kOverwrite);
        if (fHistoFracAllGammaToSecFromXFromK0l)                        fHistoFracAllGammaToSecFromXFromK0l->Write("FracAllGammaToSecFromXFromK0l",TObject::kOverwrite);
        if (fHistoFracAllGammaCaloToSecFromXFromK0s)                    fHistoFracAllGammaCaloToSecFromXFromK0s->Write("FracAllGammaCaloToSecFromXFromK0s",TObject::kOverwrite);
        if (fHistoFracAllGammaToSecFromXFromLambda)                     fHistoFracAllGammaToSecFromXFromLambda->Write("FracAllGammaToSecFromXFromLambda",TObject::kOverwrite);
        if (fHistoFracAllGammaCaloToSecFromXFromLambda)                 fHistoFracAllGammaCaloToSecFromXFromLambda->Write("FracAllGammaCaloToSecFromXFromLambda",TObject::kOverwrite);
        if (fHistoFracAllGammaToSecOrBin)                               fHistoFracAllGammaToSecOrBin->Write("FracAllGammaToSecOriginalBinning",TObject::kOverwrite);
        if (fHistoFracAllGammaCaloToSecOrBin)                           fHistoFracAllGammaCaloToSecOrBin->Write("FracAllGammaCaloToSecOriginalBinning",TObject::kOverwrite);
        if (fHistoFracAllGammaToSecFromXFromK0sOrBin)                   fHistoFracAllGammaToSecFromXFromK0sOrBin->Write("FracAllGammaToSecFromXFromK0sOriginalBinning",TObject::kOverwrite);
        if (fHistoFracAllGammaToSecFromXFromK0lOrBin)                   fHistoFracAllGammaToSecFromXFromK0lOrBin->Write("FracAllGammaToSecFromXFromK0lOriginalBinning",TObject::kOverwrite);
        if (fHistoFracAllGammaCaloToSecFromXFromK0sOrBin)               fHistoFracAllGammaCaloToSecFromXFromK0sOrBin->Write("FracAllGammaCaloToSecFromXFromK0sOriginalBinning",TObject::kOverwrite);
        if (fHistoFracAllGammaToSecFromXFromLambdaOrBin)                fHistoFracAllGammaToSecFromXFromLambdaOrBin->Write("FracAllGammaToSecFromXFromLambdaOriginalBinning",TObject::kOverwrite);
        if (fHistoFracAllGammaCaloToSecFromXFromLambdaOrBin)            fHistoFracAllGammaCaloToSecFromXFromLambdaOrBin->Write("FracAllGammaCaloToSecFromXFromLambdaOriginalBinning",TObject::kOverwrite);
    
        // write correction factors to calculate raw from cocktail secondary spectra
        if(fUseCocktail && nHistogramDimension==2){
            
            // secondary response matrices
            if (fHistoGammaTrueSecondaryFromXFromK0sConv_MCPt_recPt_MC_Rebin)
                fHistoGammaTrueSecondaryFromXFromK0sConv_MCPt_recPt_MC_Rebin->Write("TrueSecondaryConvGammaFromXFromK0s_MCPt_recPt",TObject::kOverwrite);
            if (fHistoGammaTrueSecondaryFromXFromK0lConv_MCPt_recPt_MC_Rebin)
                fHistoGammaTrueSecondaryFromXFromK0lConv_MCPt_recPt_MC_Rebin->Write("TrueSecondaryConvGammaFromXFromK0l_MCPt_recPt",TObject::kOverwrite);
            if (fHistoGammaTrueSecondaryFromXFromLambdaConv_MCPt_recPt_MC_Rebin)
                fHistoGammaTrueSecondaryFromXFromLambdaConv_MCPt_recPt_MC_Rebin->Write("TrueSecondaryConvGammaFromXFromLambda_MCPt_recPt",TObject::kOverwrite);
            if (fHistoGammaTrueSecondaryFromXFromK0sConv_MCPt_recPt_MC)
                fHistoGammaTrueSecondaryFromXFromK0sConv_MCPt_recPt_MC->Write("TrueSecondaryConvGammaFromXFromK0s_MCPt_recPt_orBin",TObject::kOverwrite);
            if (fHistoGammaTrueSecondaryFromXFromK0lConv_MCPt_recPt_MC)
                fHistoGammaTrueSecondaryFromXFromK0lConv_MCPt_recPt_MC->Write("TrueSecondaryConvGammaFromXFromK0l_MCPt_recPt_orBin",TObject::kOverwrite);
            if (fHistoGammaTrueSecondaryFromXFromLambdaConv_MCPt_recPt_MC)
                fHistoGammaTrueSecondaryFromXFromLambdaConv_MCPt_recPt_MC->Write("TrueSecondaryConvGammaFromXFromLambda_MCPt_recPt_orBin",TObject::kOverwrite);

            if (fHistoGammaTrueSecondaryFromXFromK0sCalo_MCPt_recPt_MC_Rebin)
                fHistoGammaTrueSecondaryFromXFromK0sCalo_MCPt_recPt_MC_Rebin->Write("TrueSecondaryCaloGammaFromXFromK0s_MCPt_recPt",TObject::kOverwrite);
            if (fHistoGammaTrueSecondaryFromXFromK0lCalo_MCPt_recPt_MC_Rebin)
                fHistoGammaTrueSecondaryFromXFromK0lCalo_MCPt_recPt_MC_Rebin->Write("TrueSecondaryCaloGammaFromXFromK0l_MCPt_recPt",TObject::kOverwrite);
            if (fHistoGammaTrueSecondaryFromXFromLambdaCalo_MCPt_recPt_MC_Rebin)
                fHistoGammaTrueSecondaryFromXFromLambdaCalo_MCPt_recPt_MC_Rebin->Write("TrueSecondaryCaloGammaFromXFromLambda_MCPt_recPt",TObject::kOverwrite);
            if (fHistoGammaTrueSecondaryFromXFromK0sCalo_MCPt_recPt_MC)
                fHistoGammaTrueSecondaryFromXFromK0sCalo_MCPt_recPt_MC->Write("TrueSecondaryCaloGammaFromXFromK0s_MCPt_recPt_orBin",TObject::kOverwrite);
            if (fHistoGammaTrueSecondaryFromXFromK0lCalo_MCPt_recPt_MC)
                fHistoGammaTrueSecondaryFromXFromK0lCalo_MCPt_recPt_MC->Write("TrueSecondaryCaloGammaFromXFromK0l_MCPt_recPt_orBin",TObject::kOverwrite);
            if (fHistoGammaTrueSecondaryFromXFromLambdaCalo_MCPt_recPt_MC)
                fHistoGammaTrueSecondaryFromXFromLambdaCalo_MCPt_recPt_MC->Write("TrueSecondaryCaloGammaFromXFromLambda_MCPt_recPt_orBin",TObject::kOverwrite);
            
            // secondary conv prob in MC Pt
            if (fHistoSecondaryGammaFromXFromK0sMCConvProb)             fHistoSecondaryGammaFromXFromK0sMCConvProb->Write("SecondaryGammaFromXFromK0sConvProb_MCPt",TObject::kOverwrite);
            if (fHistoSecondaryGammaFromXFromK0lMCConvProb)             fHistoSecondaryGammaFromXFromK0lMCConvProb->Write("SecondaryGammaFromXFromK0lConvProb_MCPt",TObject::kOverwrite);
            if (fHistoSecondaryGammaFromXFromLambdaMCConvProb)          fHistoSecondaryGammaFromXFromLambdaMCConvProb->Write("SecondaryGammaFromXFromLambdaConvProb_MCPt",TObject::kOverwrite);
            if (fHistoSecondaryGammaFromXFromK0sMCConvProbOrBin)        fHistoSecondaryGammaFromXFromK0sMCConvProbOrBin->Write("SecondaryGammaFromXFromK0sConvProb_MCPtOrBin",TObject::kOverwrite);
            if (fHistoSecondaryGammaFromXFromK0lMCConvProbOrBin)        fHistoSecondaryGammaFromXFromK0lMCConvProbOrBin->Write("SecondaryGammaFromXFromK0lConvProb_MCPtOrBin",TObject::kOverwrite);
            if (fHistoSecondaryGammaFromXFromLambdaMCConvProbOrBin)     fHistoSecondaryGammaFromXFromLambdaMCConvProbOrBin->Write("SecondaryGammaFromXFromLambdaConvProb_MCPtOrBin",TObject::kOverwrite);
            
            // secondary reco eff in MC Pt
            if (fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPt)          fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPt->Write("SecondaryGammaFromXFromK0sRecoEff_MCPt",TObject::kOverwrite);
            if (fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPt)          fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPt->Write("SecondaryGammaFromXFromK0lRecoEff_MCPt",TObject::kOverwrite);
            if (fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPt)       fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPt->Write("SecondaryGammaFromXFromLambdaRecoEff_MCPt",TObject::kOverwrite);
            if (fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPtOrBin)     fHistoSecondaryGammaFromXFromK0sMCRecoEffMCPtOrBin->Write("SecondaryGammaFromXFromK0sRecoEff_MCPtOrBin",TObject::kOverwrite);
            if (fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPtOrBin)     fHistoSecondaryGammaFromXFromK0lMCRecoEffMCPtOrBin->Write("SecondaryGammaFromXFromK0lRecoEff_MCPtOrBin",TObject::kOverwrite);
            if (fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPtOrBin)  fHistoSecondaryGammaFromXFromLambdaMCRecoEffMCPtOrBin->Write("SecondaryGammaFromXFromLambdaRecoEff_MCPtOrBin",TObject::kOverwrite);

            // secondary reco eff in rec Pt
            if (fHistoSecondaryGammaFromXFromK0sMCRecoEffPt)            fHistoSecondaryGammaFromXFromK0sMCRecoEffPt->Write("SecondaryGammaFromXFromK0sRecoEff_Pt",TObject::kOverwrite);
            if (fHistoSecondaryGammaFromXFromK0lMCRecoEffPt)            fHistoSecondaryGammaFromXFromK0lMCRecoEffPt->Write("SecondaryGammaFromXFromK0lRecoEff_Pt",TObject::kOverwrite);
            if (fHistoSecondaryGammaFromXFromLambdaMCRecoEffPt)         fHistoSecondaryGammaFromXFromLambdaMCRecoEffPt->Write("SecondaryGammaFromXFromLambdaRecoEff_Pt",TObject::kOverwrite);
            if (fHistoSecondaryGammaFromXFromK0sMCRecoEffPtOrBin)       fHistoSecondaryGammaFromXFromK0sMCRecoEffPtOrBin->Write("SecondaryGammaFromXFromK0sRecoEff_PtOrBin",TObject::kOverwrite);
            if (fHistoSecondaryGammaFromXFromK0lMCRecoEffPtOrBin)       fHistoSecondaryGammaFromXFromK0lMCRecoEffPtOrBin->Write("SecondaryGammaFromXFromK0lRecoEff_PtOrBin",TObject::kOverwrite);
            if (fHistoSecondaryGammaFromXFromLambdaMCRecoEffPtOrBin)    fHistoSecondaryGammaFromXFromLambdaMCRecoEffPtOrBin->Write("SecondaryGammaFromXFromLambdaRecoEff_PtOrBin",TObject::kOverwrite);

            // true secondary conv rest (= All - K0s - K0l - Lambda)
            if (fHistoGammaTrueSecondaryConvGammaRestPt)                fHistoGammaTrueSecondaryConvGammaRestPt->Write("fHistoGammaTrueSecondaryConvGammaRestPt",TObject::kOverwrite);
            if (fHistoGammaTrueSecondaryConvGammaRestPtOrBin)           fHistoGammaTrueSecondaryConvGammaRestPtOrBin->Write("fHistoGammaTrueSecondaryConvGammaRestPtOrBin",TObject::kOverwrite);

            // true secondary calo rest (= All - K0s - K0l - Lambda)
            if (fHistoGammaTrueSecondaryCaloRestPt)                     fHistoGammaTrueSecondaryCaloRestPt->Write("fHistoGammaTrueSecondaryCaloRestPt",TObject::kOverwrite);
            if (fHistoGammaTrueSecondaryCaloRestPtOrBin)                fHistoGammaTrueSecondaryCaloRestPtOrBin->Write("fHistoGammaTrueSecondaryCaloRestPtOrBin",TObject::kOverwrite);
        }
        
        
        // write double counting histograms
        if (fEnableDCConv){
            TH1D* fHistoTrueGammaConvDCPtRebinned               = NULL;
            if (fHistoTrueGammaConvDCPt){
                fHistoTrueGammaConvDCPt->Write("fHistoTrueGammaConvDCPt",TObject::kOverwrite);
                fHistoTrueGammaConvDCPtRebinned                 = (TH1D*)fHistoTrueGammaConvDCPt->Clone("fHistoTrueGammaConvDCPtRebinned");
                RebinSpectrum(fHistoTrueGammaConvDCPtRebinned);
            }
            if (fHistoTrueGammaConvDCR)             fHistoTrueGammaConvDCR->Write("fHistoTrueGammaConvDCR",TObject::kOverwrite);
            if (fHistoTrueGammaConvDCRVSPt)         fHistoTrueGammaConvDCRVSPt->Write("fHistoTrueGammaConvDCRVSPt",TObject::kOverwrite);
            if (fHistoTrueGammaConvMultipleCount)   fHistoTrueGammaConvMultipleCount->Write("fHistoTrueGammaConvMultipleCount",TObject::kOverwrite);
            if (fHistoGammaTrueConvPtOrBin && fHistoTrueGammaConvDCPt ){
                TH1D* fHistoRatioDCGammaConv                    = (TH1D*)fHistoTrueGammaConvDCPt->Clone("fHistoRatioDCGammaConv");
                fHistoRatioDCGammaConv->Divide(fHistoRatioDCGammaConv,fHistoGammaTrueConvPtOrBin,1,1,"B" );
                fHistoRatioDCGammaConv->Write("FractionDCGammaConv",TObject::kOverwrite);
                
                if (fHistoTrueGammaConvDCPtRebinned){
                    TH1D* fHistoRatioDCGammaConvRebinned        = (TH1D*)fHistoTrueGammaConvDCPtRebinned->Clone("fHistoRatioDCGammaConvRebinned");
                    fHistoRatioDCGammaConvRebinned->Divide(fHistoRatioDCGammaConvRebinned,fHistoGammaTrueConvPt,1,1,"B" );
                    fHistoRatioDCGammaConvRebinned->Write("FractionDCGammaConvRebinned",TObject::kOverwrite);
                }    
            }    
        }    
        
        // write source splitting of input photons
        for(Int_t i = 0; i<7; i++)
            if (fHistoGammaMCDecayPt[i])                    fHistoGammaMCDecayPt[i]->Write(Form("MC_DecayGamma%s_Pt",fDecays[i].Data()),TObject::kOverwrite);
        
        // write correction histograms
        // ---> Purity
        if (fHistoGammaMCPurity)                            fHistoGammaMCPurity->Write("GammaPurity_Pt",TObject::kOverwrite);
        if (!fUseDataDrivenPurity) {
            if (fHistoGammaMCTruePurity)                    fHistoGammaMCTruePurity->Write("GammaTruePurity_Pt",TObject::kOverwrite);
        } else {
            if (fHistoPurityKappaTemplates)                 fHistoPurityKappaTemplates->Write("GammaTruePurity_Pt",TObject::kOverwrite);
        }
        if (fHistoGammaMCTruePurityOrBin)                   fHistoGammaMCTruePurityOrBin->Write("GammaTruePurity_OriginalBinning_Pt",TObject::kOverwrite);
        if (fHistoGammaCaloMCPurity)                        fHistoGammaCaloMCPurity->Write("GammaCaloPurity_Pt",TObject::kOverwrite);
        if (fHistoGammaCaloMCTruePurity)                    fHistoGammaCaloMCTruePurity->Write("GammaCaloTruePurity_Pt",TObject::kOverwrite);
        if (fHistoGammaCaloMCTruePurityOrBin)               fHistoGammaCaloMCTruePurityOrBin->Write("GammaCaloTruePurity_OriginalBinning_Pt",TObject::kOverwrite);
        
        // ---> Conversion probability
        if (fHistoGammaMCConvProb)                          fHistoGammaMCConvProb->Write("GammaConvProb_MCPt",TObject::kOverwrite);
        if (fHistoGammaMCConvProbOrBin)                     fHistoGammaMCConvProbOrBin->Write("GammaConvProb_MCPt_OriginalBinning",TObject::kOverwrite);
    
        // ---> Reco efficiency
        if (fHistoGammaMCRecoEff)                           fHistoGammaMCRecoEff->Write("GammaRecoEff_Pt",TObject::kOverwrite);
        if (fHistoGammaMCPrimaryRecoEff)                    fHistoGammaMCPrimaryRecoEff->Write("GammaPrimaryRecoEff_Pt",TObject::kOverwrite);
        if (fHistoGammaMCPrimaryRecoEffOrBin)               fHistoGammaMCPrimaryRecoEffOrBin->Write("GammaPrimaryRecoEff_Pt_OriginalBinning",TObject::kOverwrite);
        if (fHistoGammaMCPrimaryRecoEffMCPt)                fHistoGammaMCPrimaryRecoEffMCPt->Write("GammaPrimaryRecoEff_MCPt",TObject::kOverwrite);
        if (fHistoGammaMCPrimaryRecoEffMCPtOrBin)           fHistoGammaMCPrimaryRecoEffMCPtOrBin->Write("GammaPrimaryRecoEff_MCPt_OriginalBinning",TObject::kOverwrite);
        if (fHistoGammaCaloMCRecoEff)                       fHistoGammaCaloMCRecoEff->Write("GammaCaloRecoEff_Pt",TObject::kOverwrite);
        if (fHistoGammaCaloMCPrimaryRecoEff)                fHistoGammaCaloMCPrimaryRecoEff->Write("GammaCaloPrimaryRecoEff_Pt",TObject::kOverwrite);
        if (fHistoGammaCaloMCPrimaryRecoEffMCPt)            fHistoGammaCaloMCPrimaryRecoEffMCPt->Write("GammaCaloPrimaryRecoEff_MCPt",TObject::kOverwrite);
        
        // write response matices
        if (fHistoGammaTruePrimaryConv_recPt_MCPt_MC)       fHistoGammaTruePrimaryConv_recPt_MCPt_MC->Write("TruePrimaryConvGamma_recPt_MCPt",TObject::kOverwrite);
        if (fHistoGammaTruePrimaryConv_recPt_MCPt_MC_Rebin) fHistoGammaTruePrimaryConv_recPt_MCPt_MC_Rebin->Write("TruePrimaryConvGamma_recPt_MCPt_Rebin",TObject::kOverwrite);
        if (fHistoGammaTruePrimaryCalo_recPt_MCPt_MC)       fHistoGammaTruePrimaryCalo_recPt_MCPt_MC->Write("TruePrimaryCaloGamma_recPt_MCPt",TObject::kOverwrite);
        if (fHistoGammaTruePrimaryCalo_recPt_MCPt_MC_Rebin) fHistoGammaTruePrimaryCalo_recPt_MCPt_MC_Rebin->Write("TruePrimaryCaloGamma_recPt_MCPt_Rebin",TObject::kOverwrite);
        
        // write total MC rec BG
        if (fHistoGammaMCBackground)                        fHistoGammaMCBackground->Write("MCrec_Background",TObject::kOverwrite);
        if (fHistoGammaCaloMCBackground)                    fHistoGammaCaloMCBackground->Write("MCrec_Calo_Background",TObject::kOverwrite);
        
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
                    ESDGammaPtDCAzPerEvent->Scale(1./fNEvents);
                    ESDGammaPtDCAzPerEvent->Write(ESDGammaPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                } else {
                    fESDGammaPtDCAz[catIter]->Write(fESDGammaPtDCAz[catIter]->GetName(), TObject::kOverwrite);
                    TH2D *ESDGammaPtDCAzPerEvent                                = (TH2D*)fESDGammaPtDCAz[catIter]->Clone(Form("%s_perEvent",fESDGammaPtDCAz[catIter]->GetName()));
                    ESDGammaPtDCAzPerEvent->Scale(1./fNEvents);
                    ESDGammaPtDCAzPerEvent->Write(ESDGammaPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                }
                
                Int_t ptBin;
                for (Int_t i = 0; i < 3; i++) {
                    if (i == 0)         ptBin                                   = 0;
                    else if (i == 1)    ptBin                                   = exemplaryLowPtBin;
                    else                ptBin                                   = exemplaryHighPtBin;
                    
                    fESDGammaPtDCAzBins[catIter][ptBin]->Write(fESDGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *ESDGammaPtDCAzBinsPerEvent                            = (TH1D*)fESDGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fESDGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    ESDGammaPtDCAzBinsPerEvent->Scale(1./fNEvents);
                    ESDGammaPtDCAzBinsPerEvent->Write(ESDGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    fESDGammaPtDCAzBinsBack[catIter][ptBin][0]->Write(fESDGammaPtDCAzBinsBack[catIter][ptBin][0]->GetName(), TObject::kOverwrite);
                    TH1D *ESDGammaPtDCAzBinsBackPerEvent                        = (TH1D*)fESDGammaPtDCAzBinsBack[catIter][ptBin][0]->Clone(Form("%s_perEvent",fESDGammaPtDCAzBinsBack[catIter][ptBin][0]->GetName()));
                    ESDGammaPtDCAzBinsBackPerEvent->Scale(1./fNEvents);
                    ESDGammaPtDCAzBinsBackPerEvent->Write(ESDGammaPtDCAzBinsBackPerEvent->GetName(), TObject::kOverwrite);
                    
                    fESDSubGammaPtDCAzBins[catIter][ptBin][0]->Write(fESDSubGammaPtDCAzBins[catIter][ptBin][0]->GetName(), TObject::kOverwrite);
                    TH1D *ESDGammaDCAzAllSubperEvent                            = (TH1D*)fESDSubGammaPtDCAzBins[catIter][ptBin][0]->Clone(Form("%s_perEvent",fESDSubGammaPtDCAzBins[catIter][ptBin][0]->GetName()));
                    ESDGammaDCAzAllSubperEvent->Scale(1./fNEvents);
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
                    ESDGammaPtDCAzPerEvent->Scale(1./fNEvents);
                    ESDGammaPtDCAzPerEvent->Write(ESDGammaPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTruePrimaryPhotonPtDCAz[5]->Write(fTruePrimaryPhotonPtDCAz[5]->GetName(), TObject::kOverwrite);
                    TH2D *TruePrimaryPhotonPtDCAzPerEvent                       = (TH2D*)fTruePrimaryPhotonPtDCAz[5]->Clone(Form("%s_perEvent",fTruePrimaryPhotonPtDCAz[5]->GetName()));
                    TruePrimaryPhotonPtDCAzPerEvent->Scale(1./fNEvents);
                    TruePrimaryPhotonPtDCAzPerEvent->Write(TruePrimaryPhotonPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueSecondaryPhotonPtDCAz[5]->Write(fTrueSecondaryPhotonPtDCAz[5]->GetName(), TObject::kOverwrite);
                    TH2D *TrueSecondaryPhotonPtDCAzPerEvent                     = (TH2D*)fTrueSecondaryPhotonPtDCAz[5]->Clone(Form("%s_perEvent",fTrueSecondaryPhotonPtDCAz[5]->GetName()));
                    TrueSecondaryPhotonPtDCAzPerEvent->Scale(1./fNEvents);
                    TrueSecondaryPhotonPtDCAzPerEvent->Write(TrueSecondaryPhotonPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueSecondaryPhotonFromXFromK0sPtDCAz[5]->Write(fTrueSecondaryPhotonFromXFromK0sPtDCAz[5]->GetName(), TObject::kOverwrite);
                    TH2D *TrueSecondaryPhotonFromXFromK0sPtDCAzPerEvent         = (TH2D*)fTrueSecondaryPhotonFromXFromK0sPtDCAz[5]->Clone(Form("%s_perEvent",fTrueSecondaryPhotonFromXFromK0sPtDCAz[5]->GetName()));
                    TrueSecondaryPhotonFromXFromK0sPtDCAzPerEvent->Scale(1./fNEvents);
                    TrueSecondaryPhotonFromXFromK0sPtDCAzPerEvent->Write(TrueSecondaryPhotonFromXFromK0sPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                } else {
                    fESDGammaPtDCAz[catIter]->Write(Form("MCrec_GammaPtDCAz_%s", categoryName[catIter].Data()), TObject::kOverwrite);
                    TH2D *ESDGammaPtDCAzPerEvent                                = (TH2D*)fESDGammaPtDCAz[catIter]->Clone(Form("MCrec_GammaPtDCAz_%s_perEvent",categoryName[catIter].Data()));
                    ESDGammaPtDCAzPerEvent->Scale(1./fNEvents);
                    ESDGammaPtDCAzPerEvent->Write(ESDGammaPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTruePrimaryPhotonPtDCAz[catIter]->Write(fTruePrimaryPhotonPtDCAz[catIter]->GetName(), TObject::kOverwrite);
                    TH2D *TruePrimaryPhotonPtDCAzPerEvent                       = (TH2D*)fTruePrimaryPhotonPtDCAz[catIter]->Clone(Form("%s_perEvent",fTruePrimaryPhotonPtDCAz[catIter]->GetName()));
                    TruePrimaryPhotonPtDCAzPerEvent->Scale(1./fNEvents);
                    TruePrimaryPhotonPtDCAzPerEvent->Write(TruePrimaryPhotonPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueSecondaryPhotonPtDCAz[catIter]->Write(fTrueSecondaryPhotonPtDCAz[catIter]->GetName(), TObject::kOverwrite);
                    TH2D *TrueSecondaryPhotonPtDCAzPerEvent                     = (TH2D*)fTrueSecondaryPhotonPtDCAz[catIter]->Clone(Form("%s_perEvent",fTrueSecondaryPhotonPtDCAz[catIter]->GetName()));
                    TrueSecondaryPhotonPtDCAzPerEvent->Scale(1./fNEvents);
                    TrueSecondaryPhotonPtDCAzPerEvent->Write(TrueSecondaryPhotonPtDCAzPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueSecondaryPhotonFromXFromK0sPtDCAz[catIter]->Write(fTrueSecondaryPhotonFromXFromK0sPtDCAz[catIter]->GetName(), TObject::kOverwrite);
                    TH2D *TrueSecondaryPhotonFromXFromK0sPtDCAzPerEvent         = (TH2D*)fTrueSecondaryPhotonFromXFromK0sPtDCAz[catIter]->Clone(Form("%s_perEvent",fTrueSecondaryPhotonFromXFromK0sPtDCAz[catIter]->GetName()));
                    TrueSecondaryPhotonFromXFromK0sPtDCAzPerEvent->Scale(1./fNEvents);
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
                    MCrecGammaPtDCAzBinsPerEvent->Scale(1./fNEvents);
                    MCrecGammaPtDCAzBinsPerEvent->Write(MCrecGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);

                    fMCrecGammaPtDCAzBinsBack[catIter][ptBin]->Write(fMCrecGammaPtDCAzBinsBack[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *MCrecGammaPtDCAzBinsBackPerEvent                      = (TH1D*)fMCrecGammaPtDCAzBinsBack[catIter][ptBin]->Clone(Form("%s_perEvent",fMCrecGammaPtDCAzBinsBack[catIter][ptBin]->GetName()));
                    MCrecGammaPtDCAzBinsBackPerEvent->Scale(1./fNEvents);
                    MCrecGammaPtDCAzBinsBackPerEvent->Write(MCrecGammaPtDCAzBinsBackPerEvent->GetName(), TObject::kOverwrite);
                    
                    fMCrecSubGammaPtDCAzBins[catIter][ptBin]->Write(fMCrecSubGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *MCrecSubGammaPtDCAzBinsPerEvent                       = (TH1D*)fMCrecSubGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fMCrecSubGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    MCrecSubGammaPtDCAzBinsPerEvent->Scale(1./fNEvents);
                    MCrecSubGammaPtDCAzBinsPerEvent->Write(MCrecSubGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    // true gamma
                    fTrueGammaPtDCAzBins[catIter][ptBin]->Write(fTrueGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TrueGammaPtDCAzBinsPerEvent                       = (TH1D*)fTrueGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTrueGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    TrueGammaPtDCAzBinsPerEvent->Scale(1./fNEvents);
                    TrueGammaPtDCAzBinsPerEvent->Write(TrueGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueSubGammaPtDCAzBins[catIter][ptBin]->Write(fTrueSubGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TrueSubGammaPtDCAzBinsPerEvent                       = (TH1D*)fTrueSubGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTrueSubGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    TrueSubGammaPtDCAzBinsPerEvent->Scale(1./fNEvents);
                    TrueSubGammaPtDCAzBinsPerEvent->Write(TrueSubGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueBackgroundPtDCAzBins[catIter][ptBin]->Write(fTrueBackgroundPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TrueBackgroundPtDCAzBinsPerEvent                       = (TH1D*)fTrueBackgroundPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTrueBackgroundPtDCAzBins[catIter][ptBin]->GetName()));
                    TrueBackgroundPtDCAzBinsPerEvent->Scale(1./fNEvents);
                    TrueBackgroundPtDCAzBinsPerEvent->Write(TrueBackgroundPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);                    
                    
                    // true primary
                    fTruePrimaryGammaPtDCAzBins[catIter][ptBin]->Write(fTruePrimaryGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TruePrimaryGammaPtDCAzBinsPerEvent                    = (TH1D*)fTruePrimaryGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTruePrimaryGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    TruePrimaryGammaPtDCAzBinsPerEvent->Scale(1./fNEvents);
                    TruePrimaryGammaPtDCAzBinsPerEvent->Write(TruePrimaryGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTruePrimarySubGammaPtDCAzBins[catIter][ptBin]->Write(fTruePrimarySubGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TruePrimarySubGammaPtDCAzBinsPerEvent                 = (TH1D*)fTruePrimarySubGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTruePrimarySubGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    TruePrimarySubGammaPtDCAzBinsPerEvent->Scale(1./fNEvents);
                    TruePrimarySubGammaPtDCAzBinsPerEvent->Write(TruePrimarySubGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    // true secondary
                    fTrueSecondaryGammaPtDCAzBins[catIter][ptBin]->Write(fTrueSecondaryGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TrueSecondaryGammaPtDCAzBinsPerEvent                  = (TH1D*)fTrueSecondaryGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTrueSecondaryGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    TrueSecondaryGammaPtDCAzBinsPerEvent->Scale(1./fNEvents);
                    TrueSecondaryGammaPtDCAzBinsPerEvent->Write(TrueSecondaryGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueSecondarySubGammaPtDCAzBins[catIter][ptBin]->Write(fTrueSecondarySubGammaPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TrueSecondarySubGammaPtDCAzBinsPerEvent               = (TH1D*)fTrueSecondarySubGammaPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTrueSecondarySubGammaPtDCAzBins[catIter][ptBin]->GetName()));
                    TrueSecondarySubGammaPtDCAzBinsPerEvent->Scale(1./fNEvents);
                    TrueSecondarySubGammaPtDCAzBinsPerEvent->Write(TrueSecondarySubGammaPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    // true secondary from X from K0s
                    fTrueSecondaryGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->Write(fTrueSecondaryGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TrueSecondaryGammaFromXFromK0sPtDCAzBinsPerEvent      = (TH1D*)fTrueSecondaryGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTrueSecondaryGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->GetName()));
                    TrueSecondaryGammaFromXFromK0sPtDCAzBinsPerEvent->Scale(1./fNEvents);
                    TrueSecondaryGammaFromXFromK0sPtDCAzBinsPerEvent->Write(TrueSecondaryGammaFromXFromK0sPtDCAzBinsPerEvent->GetName(), TObject::kOverwrite);
                    
                    fTrueSecondarySubGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->Write(fTrueSecondarySubGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->GetName(), TObject::kOverwrite);
                    TH1D *TrueSecondarySubGammaFromXFromK0sPtDCAzBinsPerEvent   = (TH1D*)fTrueSecondarySubGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->Clone(Form("%s_perEvent",fTrueSecondarySubGammaFromXFromK0sPtDCAzBins[catIter][ptBin]->GetName()));
                    TrueSecondarySubGammaFromXFromK0sPtDCAzBinsPerEvent->Scale(1./fNEvents);
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
        DrawGammaSetMarker(RatioWithWithoutPileUpMC, 24, 1.0, kRed, kRed);

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
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
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
            TLatex *events                      = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));
            
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
            legendData->AddEntry(ESDGammaPtDCAzBins[iPt],textEvents,"l");
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
            if (remaining > 0)  padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else                padDataSpectra->cd(place)->SetLeftMargin(0.25);
            
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
    for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
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
            TLatex *events                      = new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));
            
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
            //legendData->AddEntry(ESDGammaPtDCAzBins[iPt],ESDGammaPtDCAzBins[iPt]->GetName(),"l");
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

        if(fHistoGGInvMassPtGConvBin[iPt]!= NULL){
            delete fHistoGGInvMassPtGConvBin[iPt];
            fHistoGGInvMassPtGConvBin[iPt]                      = NULL;
        }
        fHistoGGInvMassPtGConvBin[iPt]                          = new TH1D(fNameHistoGG.Data(),fNameHistoGG.Data(),fGammaGammaInvMassVSPtDummy->GetNbinsX(),
                                                                           0.,fGammaGammaInvMassVSPtDummy->GetXaxis()->GetBinUpEdge(fGammaGammaInvMassVSPtDummy->GetNbinsX()));
        fHistoGGInvMassPtGConvBin[iPt]->Sumw2();

        Int_t startBin                                          = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
        Int_t endBin                                            = fGammaGammaInvMassVSPtDummy->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

        fGammaGammaInvMassVSPtDummy->ProjectionX(fNameHistoGG.Data(),startBin,endBin);
        fHistoGGInvMassPtGConvBin[iPt]                          = (TH1D*)gDirectory->Get(fNameHistoGG.Data());
        if(fNRebin[iPt]>1){
            fHistoGGInvMassPtGConvBin[iPt]->Rebin(fNRebin[iPt]);
        }
    }
}

//****************************************************************************
//************** Produce background with proper weighting ********************
//****************************************************************************
void ProduceBckProperWeighting(TH2D* fHistoMotherZM, TH2D* fHistoBckZM){

    cout << "Using TH2 for the background" << endl;
    for(Int_t iPt=fStartPtBin;iPt<fNBinsPt;iPt++){
        Int_t startBin                                          = fHistoMotherZM->GetYaxis()->FindBin(fBinsPt[iPt]+0.001);
        Int_t endBin                                            = fHistoMotherZM->GetYaxis()->FindBin(fBinsPt[iPt+1]-0.001);

        fHistoMotherZMProj                                      = fHistoMotherZM->ProjectionX("ProjectMother",startBin,endBin);
        fHistoMotherZMProj->Sumw2();
        fHistoBckZMProj                                         = fHistoBckZM->ProjectionX("ProjectBck",startBin,endBin);
        fHistoBckZMProj->Sumw2();

        fScalingFactorBck[0][0]= 1./fBackgroundMultNumber;
        if(fHistoBackInvMassPtGconvBin[iPt]!= NULL){
            delete fHistoBackInvMassPtGconvBin[iPt];
            fHistoBackInvMassPtGconvBin[iPt]                    = NULL;
        }
        fNameHistoBack = Form("Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
        fHistoBackInvMassPtGconvBin[iPt]= (TH1D*)fHistoBckZMProj->Clone(fNameHistoBack);
        fHistoBackInvMassPtGconvBin[iPt]->Sumw2();
        
        Int_t startBinIntegral                                  = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[0]);
        Int_t endBinIntegral                                    = fHistoMotherZMProj->GetXaxis()->FindBin(fBGFitRange[1]);
        if (fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral) != 0) {
            fScalingFactorBck[0][0]                             = fHistoMotherZMProj->Integral(startBinIntegral,endBinIntegral)/fHistoBckZMProj->Integral(startBinIntegral,endBinIntegral);
            if ( fScalingFactorBck[0][0]>20./fBackgroundMultNumber ){
                fScalingFactorBck[0][0]=1./fBackgroundMultNumber;
            }
        }
        
        fHistoBackInvMassPtGconvBin[iPt]->Add(fHistoBckZMProj,fScalingFactorBck[0][0]);
        fHistoBackInvMassPtGconvBin[iPt]->Rebin(fNRebin[iPt]);
        for (Int_t ii = 0; ii < fHistoBackInvMassPtGconvBin[iPt]->GetNbinsX()+1; ii++){
            if(fHistoBackInvMassPtGconvBin[iPt]->GetBinContent(ii) == 0){
                fHistoBackInvMassPtGconvBin[iPt]->SetBinContent(ii,0.);
                fHistoBackInvMassPtGconvBin[iPt]->SetBinError(ii,1.);
            }
        }
        cout << "Scaling Background factors for Pt bin " << iPt << endl;
        cout << fScalingFactorBck[0][0] << endl;
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
    
    Int_t fFitStartBin = -1;
    Int_t fitStatus;
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


//******************************************************************************
// Load secondary gamma histos from cocktail file
// - put them in proper scaling
// - rebin them according to current gamma binning
//******************************************************************************
Bool_t LoadSecondariesFromCocktailFile(TString cutSelection, TString optionEnergy){
    
    Double_t xMax                                               = 20;
    if (fEnablePCM) xMax                                        = fHistoGammaConvPtOrBin->GetXaxis()->GetXmax();
    else            xMax                                        = fHistoGammaCaloPtOrBin->GetXaxis()->GetXmax();

    // search for cocktail file
    TString nameCocktailFile                                    = Form("%s/%s/SecondaryGamma%s_%.2f_%s.root",cutSelection.Data(),optionEnergy.Data(),fPeriodFlag.Data(),fYMaxMeson/2,cutSelection.Data());
    fFileCocktailInput                                          = new TFile(nameCocktailFile.Data());
    if (fFileCocktailInput->IsZombie()) fFileCocktailInput      = NULL;

    // get secondary spectra from cocktail file
    if (fFileCocktailInput){
        cout << "Found cocktail file: " << nameCocktailFile.Data() << " -> will add cocktail histos to output" << endl;
        
        // K0s
        fHistoSecondaryGammaCocktailFromXFromK0sPt              = (TH1D*)fFileCocktailInput->Get("Gamma_From_X_From_K0s_Pt_OrBin");
        if(!fHistoSecondaryGammaCocktailFromXFromK0sPt){
            cout << "Gamma_From_X_From_K0s_Pt_OrBin not found in cocktail file! Cocktail will not be used." << endl;
            return kFALSE;
        }
        fHistoSecondaryGammaCocktailFromXFromK0sPt->Sumw2();
        fHistoSecondaryGammaCocktailFromXFromK0sPt->Scale(fHistoSecondaryGammaCocktailFromXFromK0sPt->GetBinWidth(1));
        fHistoSecondaryGammaCocktailFromXFromK0sPtOrBin         = (TH1D*)fHistoSecondaryGammaCocktailFromXFromK0sPt->Clone("CocktailSecondaryGammaFromXFromK0s_PtOrBin");
        RebinSpectrum(fHistoSecondaryGammaCocktailFromXFromK0sPt);
        
        // K0l
        fHistoSecondaryGammaCocktailFromXFromK0lPt              = (TH1D*)fFileCocktailInput->Get("Gamma_From_X_From_K0l_Pt_OrBin");
        if(!fHistoSecondaryGammaCocktailFromXFromK0lPt){
            cout << "Gamma_From_X_From_K0l_Pt_OrBin not found in cocktail file! Cocktail will not be used." << endl;
            return kFALSE;
        }
        fHistoSecondaryGammaCocktailFromXFromK0lPt->Sumw2();
        fHistoSecondaryGammaCocktailFromXFromK0lPt->Scale(fHistoSecondaryGammaCocktailFromXFromK0lPt->GetBinWidth(1));
        fHistoSecondaryGammaCocktailFromXFromK0lPtOrBin         = (TH1D*)fHistoSecondaryGammaCocktailFromXFromK0lPt->Clone("CocktailSecondaryGammaFromXFromK0l_PtOrBin");
        RebinSpectrum(fHistoSecondaryGammaCocktailFromXFromK0lPt);
        
        // Lambda
        fHistoSecondaryGammaCocktailFromXFromLambdaPt           = (TH1D*)fFileCocktailInput->Get("Gamma_From_X_From_Lambda_Pt_OrBin");
        if(!fHistoSecondaryGammaCocktailFromXFromLambdaPt){
            cout << "Gamma_From_X_From_Lambda_Pt_OrBin not found in cocktail file! Cocktail will not be used." << endl;
            return kFALSE;
        }
        fHistoSecondaryGammaCocktailFromXFromLambdaPt->Sumw2();
        fHistoSecondaryGammaCocktailFromXFromLambdaPt->Scale(fHistoSecondaryGammaCocktailFromXFromLambdaPt->GetBinWidth(1));
        fHistoSecondaryGammaCocktailFromXFromLambdaPtOrBin      = (TH1D*)fHistoSecondaryGammaCocktailFromXFromLambdaPt->Clone("CocktailSecondaryGammaFromXFromLambda_PtOrBin");
        RebinSpectrum(fHistoSecondaryGammaCocktailFromXFromLambdaPt);
        
        // set bins
        Int_t       nBins   = 0;
        Double_t    xMin    = 0.;
        Double_t    xMax    = 0.;
        if (fEnablePCM && !fEnableCalo) {
            nBins           = fHistoGammaConvPtOrBin->GetNbinsX();
            xMin            = fHistoGammaConvPtOrBin->GetXaxis()->GetXmin();
            xMax            = fHistoGammaConvPtOrBin->GetXaxis()->GetXmax();
        } else if (fEnableCalo && !fEnablePCM) {
            nBins           = fHistoGammaCaloPtOrBin->GetNbinsX();
            xMin            = fHistoGammaCaloPtOrBin->GetXaxis()->GetXmin();
            xMax            = fHistoGammaCaloPtOrBin->GetXaxis()->GetXmax();
        } else {
            cout << "Will use " << fHistoGammaConvPtOrBin->GetName() << " for original binning of secondary gamma spectra from cocktail." << endl;
            nBins           = fHistoGammaConvPtOrBin->GetNbinsX();
            xMin            = fHistoGammaConvPtOrBin->GetXaxis()->GetXmin();
            xMax            = fHistoGammaConvPtOrBin->GetXaxis()->GetXmax();
        }
        fHistoSecondaryGammaCocktailFromXFromK0sPtOrBin->SetBins(   nBins,xMin,xMax);
        fHistoSecondaryGammaCocktailFromXFromK0lPtOrBin->SetBins(   nBins,xMin,xMax);
        fHistoSecondaryGammaCocktailFromXFromLambdaPtOrBin->SetBins(nBins,xMin,xMax);
        
        // all spectra found
        return kTRUE;
    }
}

//****************************************************************************
//***************** Delete intialized arrays *********************************
//****************************************************************************
void Delete(){
    if (fBGFitRange)                                            delete fBGFitRange;
    if (fBGFitRangeLeft)                                        delete fBGFitRangeLeft;

    for (Int_t i = 0; i< fNBinsPt; i++){
        if (fHistoGGInvMassPtGConvBin[i])                       delete fHistoGGInvMassPtGConvBin[i];
        if (fHistoSignalInvMassPtGConvBin[i])                   delete fHistoSignalInvMassPtGConvBin[i];
        if (fHistoSignalInvMassLeftPtGConvBin[i])               delete fHistoSignalInvMassLeftPtGConvBin[i];
        if (fHistoBackInvMassPtGconvBin[i])                     delete fHistoBackInvMassPtGconvBin[i];
        if (fHistoBackNormInvMassLeftPtGconvBin[i])             delete fHistoBackNormInvMassLeftPtGconvBin[i];
        if (fHistoBackNormInvMassPtGconvBin[i])                 delete fHistoBackNormInvMassPtGconvBin[i];
        if (fFitSignalInvMassPtBin[i])                          delete fFitSignalInvMassPtBin[i];
        if (fFitInvMassLeftPtBin[i])                            delete fFitInvMassLeftPtBin[i];
        if (fFitBckInvMassPtBin[i])                             delete fFitBckInvMassPtBin[i];
        if (fFitBckInvMassLeftPtBin[i])                         delete fFitBckInvMassLeftPtBin[i];
    } 
    if (fHistoGGInvMassPtGConvBin)                              delete fHistoGGInvMassPtGConvBin;
    if (fHistoSignalInvMassPtGConvBin)                          delete fHistoSignalInvMassPtGConvBin;
    if (fHistoSignalInvMassLeftPtGConvBin)                      delete fHistoSignalInvMassLeftPtGConvBin;
    if (fHistoBackInvMassPtGconvBin)                            delete fHistoBackInvMassPtGconvBin;
    if (fHistoBackNormInvMassLeftPtGconvBin)                    delete fHistoBackNormInvMassLeftPtGconvBin;
    if (fHistoBackNormInvMassPtGconvBin)                        delete fHistoBackNormInvMassPtGconvBin;
    if (fFitSignalInvMassPtBin)                                 delete fFitSignalInvMassPtBin;
    if (fFitInvMassLeftPtBin)                                   delete fFitInvMassLeftPtBin;
    if (fFitBckInvMassPtBin)                                    delete fFitBckInvMassPtBin;
    if (fFitBckInvMassLeftPtBin)                                delete fFitBckInvMassLeftPtBin;
    if (fMesonMass)                                             delete fMesonMass;
    if (fMesonMassError)                                        delete fMesonMassError;
    if (fMesonMassLeft)                                         delete fMesonMassLeft;
    if (fMesonMassLeftError)                                    delete fMesonMassLeftError;
    if (fMesonFWHM)                                             delete fMesonFWHM;
    if (fMesonFWHMError)                                        delete fMesonFWHMError;
    if (fMesonFWHMLeft)                                         delete fMesonFWHMLeft;
    if (fMesonFWHMLeftError)                                    delete fMesonFWHMLeftError;
    if (fMesonLambdaTailpar)                                    delete fMesonLambdaTailpar;
    if (fMesonLambdaTailparError)                               delete fMesonLambdaTailparError;
    if (fMesonLambdaTailMCpar)                                  delete fMesonLambdaTailMCpar;
    if (fMesonLambdaTailMCparError)                             delete fMesonLambdaTailMCparError;
    
    for (Int_t k = 0; k < 6; k++){   
        // delete arrays for yields
        if (fGGYields[k])                                       delete fGGYields[k];
        if (fBckYields[k])                                      delete fBckYields[k];
        if (fMesonYields[k])                                    delete fMesonYields[k];
        if (fMesonYieldsResidualBckFunc[k])                     delete fMesonYieldsResidualBckFunc[k];
        if (fMesonYieldsCorResidualBckFunc[k])                  delete fMesonYieldsCorResidualBckFunc[k];
        if (fMesonYieldsPerEvent[k])                            delete fMesonYieldsPerEvent[k];
        
        // delete arrays for errors
        if (fGGYieldsError[k])                                  delete fGGYieldsError[k];
        if (fBckYieldsError[k])                                 delete fBckYieldsError[k];
        if (fMesonYieldsError[k])                               delete fMesonYieldsError[k];
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
    }    

    for (Int_t m = 0; m < 4; m++){
        if (fMesonChi2[m])                                      delete fMesonChi2[m];
    }    
}

//****************************************************************************
//******** Fit of Signal+ BG with Gaussian + Exponential + Linear BG *********
//****************************************************************************
void FitSubtractedInvMassInPtBins(TH1D* fHistoSignalInvMassPtGConvBinSingle, Double_t* fMesonIntDeltaRangeFit, Int_t ptBin, Bool_t vary){

    //    cout<<"Start Fitting spectra"<<endl;
    fHistoSignalInvMassPtGConvBinSingle->GetXaxis()->SetRangeUser(fMesonMassPlotRange[0],fMesonMassPlotRange[1]);
    Double_t mesonAmplitude = fHistoSignalInvMassPtGConvBinSingle->GetMaximum();
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
            mesonAmplitudeMin = mesonAmplitude*30./100.;
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
            if (fMode == 0 && !fEnergyFlag.CompareTo("8TeV") ){
                mesonAmplitudeMin = mesonAmplitude*65./100.;
                if(ptBin > 2)  mesonAmplitudeMin = mesonAmplitude*85./100.;
                mesonAmplitudeMax = mesonAmplitude*115./100.;
                if(ptBin < 3) mesonAmplitudeMax = mesonAmplitude*100/100.;
            }
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


    cout << mesonAmplitude << "\t" << mesonAmplitudeMin << "-"<< mesonAmplitudeMax << endl; 
    cout << fMesonMassExpect << "\t" << mesonAmplitudeMin << "-"<< mesonAmplitudeMax << endl; 
    
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
    

    fHistoSignalInvMassPtGConvBinSingle->Fit(fFitReco,"QRME0");
    fHistoSignalInvMassPtGConvBinSingle->Fit(fFitReco,"QRME0");

    
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
        
        if (!(fMesonLambdaTail == fMesonLambdaTailRange[0] && fMesonLambdaTail == fMesonLambdaTailRange[1]) ){
            fMesonLambdaTail = fFitReco->GetParameter(3);
            fMesonLambdaTailRange[0] = 0.9*fFitReco->GetParameter(3);
            fMesonLambdaTailRange[1] = 1.1*fFitReco->GetParameter(3);
        }    
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
        binCenterStart = fHistoSignalInvMassPtGConvBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+fMesonIntDeltaRangeFit[0]);
        startBinEdge = fHistoSignalInvMassPtGConvBinSingle->GetBinCenter(binCenterStart)- 0.5*fHistoSignalInvMassPtGConvBinSingle->GetBinWidth(10);
        binCenterEnd = fHistoSignalInvMassPtGConvBinSingle->GetXaxis()->FindBin(fFitReco->GetParameter(1)+fMesonIntDeltaRangeFit[1]);
        endBinEdge = fHistoSignalInvMassPtGConvBinSingle->GetBinCenter(binCenterEnd)+ 0.5*fHistoSignalInvMassPtGConvBinSingle->GetBinWidth(10);

        Int_t nFreePar = fFitReco->GetNumberFreeParameters();
        double * covMatrix = fitter->GetCovarianceMatrix();

        Float_t intLinearBack = fFitLinearBck->GetParameter(0)*(endBinEdge-startBinEdge)+
            0.5*fFitLinearBck->GetParameter(1)*(endBinEdge*endBinEdge-startBinEdge*startBinEdge);

        Float_t errorLinearBck = pow((pow( (endBinEdge-startBinEdge)*fFitReco->GetParError(4),2)+pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fFitReco->GetParError(5),2)+2*covMatrix[nFreePar*nFreePar-2]*(endBinEdge-startBinEdge)*0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5);

        fFileDataLog << "Parameter for bin " << ptBin << endl;
        fFileDataLog << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<< "\t "<< fFitReco->GetParameter(3) <<"+-" << fFitReco->GetParError(3)<<endl;
        fFileDataLog << "Linear: \t"<<fFitReco->GetParameter(4)<<"+-" << fFitReco->GetParError(4) << "\t "<<fFitReco->GetParameter(5) <<"+-" << fFitReco->GetParError(5)<< endl;

        fIntLinearBck = intLinearBack/fHistoSignalInvMassPtGConvBinSingle->GetBinWidth(10);
        fIntLinearBckError = errorLinearBck/fHistoSignalInvMassPtGConvBinSingle->GetBinWidth(10);
    } else {
        fFileErrLog << "Fitting failed in " << ptBin << " with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }
    fFitReco->DrawCopy("same");
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
//*** Integration of Invariant Mass Histogram in given integration window ****
//****************************************************************************
void IntegrateHistoInvMass(TH1D * fHistoSignalInvMassPtGConvBinSingle, Double_t * fMesonIntRangeInt)
{
    Int_t binLowMassMeson = fHistoSignalInvMassPtGConvBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[0]);
    Int_t binHighMassMeson = fHistoSignalInvMassPtGConvBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[1]);
    fYields = fHistoSignalInvMassPtGConvBinSingle->IntegralAndError(binLowMassMeson,binHighMassMeson,fYieldsError);
}


//****************************************************************************
//*** Integration of Invariant Mass Histogram in given integration window ****
//*** with detailed output to log file ***************************************
//****************************************************************************
void IntegrateHistoInvMassStream(TH1D * fHistoSignalInvMassPtGConvBinSingle, Double_t * fMesonIntRangeInt) {
    Int_t binLowMassMeson = fHistoSignalInvMassPtGConvBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[0]);
    Int_t binHighMassMeson = fHistoSignalInvMassPtGConvBinSingle->GetXaxis()->FindBin(fMesonIntRangeInt[1]);
    fYields = fHistoSignalInvMassPtGConvBinSingle->IntegralAndError(binLowMassMeson,binHighMassMeson,fYieldsError);
    for ( Int_t M = binLowMassMeson; M < binHighMassMeson+1; M++){
        fFileDataLog << M << "\t" << fHistoSignalInvMassPtGConvBinSingle->GetBinCenter(M) <<"\t" <<fHistoSignalInvMassPtGConvBinSingle->GetBinContent(M)<< "+-"<< fHistoSignalInvMassPtGConvBinSingle->GetBinError(M)<< endl;
    }
}
