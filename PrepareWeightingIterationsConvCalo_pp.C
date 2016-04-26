/****************************************************************************************************************************
******      Friederike Bock, friederike.bock@cern.ch                                      *****
*****************************************************************************************************************************/

#include <Riostream.h>
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
#include "TBenchmark.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"

extern TRandom*      gRandom;
extern TBenchmark*   gBenchmark;
extern TSystem*      gSystem;
extern TMinuit*      gMinuit;

void ScaleMCYield(TH1D* histoCorrectedToBeScaled, Double_t deltaRapid, Double_t scaling, Double_t nEvtMC, TString nameMeson, TString optionDalitz ){
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
        if (optionDalitz.CompareTo("kFALSE")==0){
            histoCorrectedToBeScaled->Scale(1./0.98798);
        } else {
            histoCorrectedToBeScaled->Scale(1./0.01198);
        }
    }else{
        if (optionDalitz.CompareTo("kFALSE")==0){
            histoCorrectedToBeScaled->Scale(1./0.3931);
        } else {
            histoCorrectedToBeScaled->Scale(1./6.8e-5);
        }
    }
}


void PrepareWeightingIterationsConvCalo_pp( TString suffix = "pdf", 
                                            TString nameFilepp2760GeV = "data_PCMResults_pp.root", 
                                            TString nameFilepp8TeV = "data_PCMResults_pp.root", 
                                            Bool_t runDrawReweighted = kTRUE, 
                                            TString stringIterationNumber = ""){

    gROOT->Reset();   
    gROOT->SetStyle("Plain");

    TString dateForOutput                                           = ReturnDateStringForOutput();
    TString outputDir                                               = Form("%s/%s/PrepareWeightingIterationsConvCalo_pp",suffix.Data(),dateForOutput.Data());
    
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputFilePCMEMCALPionpp2760GeV.root ",nameFilepp2760GeV.Data(),outputDir.Data() ));
//     gSystem->Exec(Form("cp %s %s/InputFilePCMPionpp8TeV.root ",nameFilepp8TeV.Data(),outputDir.Data() ));

    StyleSettingsThesis();  
    SetPlotStyle();
    
    Color_t     color2760GeV                                        = kMagenta+2;
//     Color_t     color2760GeVBox                                  = color2760GeV-10;

//     Color_t     colorMCPythiaPP900GeV                            = color900GeV-4;
    Color_t     colorMCPythiaPP2760GeV                              = color2760GeV+2;
    Color_t     colorMCPythia2PP2760GeV                             = color2760GeV+1;
    Color_t     colorMCPythiaAddSigPP2760GeV                        = color2760GeV-6;
    Color_t     colorMCPhojetPP2760GeV                              = color2760GeV-4;
        
    Style_t     markerStyleSpectrum2760GeV                          = 29 ;
    Size_t      markerSizePP2760GeV                                 = 2.2;

    TString     collisionSystemPP2760GeV                            = ReturnFullCollisionsSystem("2.76TeV");
    TString     detectionSystemEMCALPCM                             = ReturnFullTextReconstructionProcess(2);
    Color_t     colorPCMEMCAL                                       = GetDefaultColorDiffDetectors("PCM-EMCal", kFALSE, kFALSE, kTRUE);
    Style_t     markerStylePCMEMCAL                                 = GetDefaultMarkerStyleDiffDetectors("PCM-EMCal", kFALSE);
    Size_t      markerSizePCMEMCAL                                  = GetDefaultMarkerSizeDiffDetectors("PCM-EMCal", kFALSE);
    
    
    TString     nameHistoPCMEMCAL                                   = "CorrectedYieldPi0";
    TString     nameHistoPCMEMCALEta                                = "CorrectedYieldEta";

    Style_t     lineStyleMCA                                       = 1;
    Style_t     lineStyleMCB                                       = 2;
    Style_t     lineStyleMCC                                       = 3;
    Style_t     lineStyleMCAddSig                                  = 4 ;
    
    Style_t     markerStyleMCA                                     = 24;
    Style_t     markerStyleMCB                                     = 25;
    Style_t     markerStyleMCC                                     = 30;
    Style_t     markerStyleMCAddSig                                = 29;
    
    // read reconstructed data
    TFile*      filePCMEMCALpp2760GeV                               = new TFile(nameFilepp2760GeV);
    TDirectory* directoryPCMEMCALPi02760GeV                         = (TDirectory*)filePCMEMCALpp2760GeV->Get("Pi02.76TeV"); 
    TH1D*       histoPCMEMCALYieldPi02760GeV                        = (TH1D*)directoryPCMEMCALPi02760GeV->Get(nameHistoPCMEMCAL.Data());
    TH1D*       histoPCMEMCALRAWYieldPi02760GeV                     = (TH1D*)directoryPCMEMCALPi02760GeV->Get("RAWYieldPerEventsPi0");
    TH1D*       histoPi0InputFullReweighted2760GeV                  = (TH1D*)directoryPCMEMCALPi02760GeV->Get("Pi0_Input_Reweighted");
    TH1D*       histoPi0InputFull2760GeV                            = (TH1D*)directoryPCMEMCALPi02760GeV->Get("Pi0_Input");
    TH1D*       histoPi0InputFullAddSigReweighted2760GeV            = (TH1D*)directoryPCMEMCALPi02760GeV->Get("Pi0_Input_Reweighted_AddedSig");
    TH1D*       histoPi0InputFullAddSig2760GeV                      = (TH1D*)directoryPCMEMCALPi02760GeV->Get("Pi0_Input_AddedSig");
    TH1D*       histoPi0Mass2760GeV                                 = (TH1D*)directoryPCMEMCALPi02760GeV->Get("MassPi0");
    histoPi0Mass2760GeV->Scale(1000.);
    TH1D*       histoPi0FWHM2760GeV                                 = (TH1D*)directoryPCMEMCALPi02760GeV->Get("FWHMPi0MeV");
    TH1D*       histoPi0TrueMass2760GeV                             = (TH1D*)directoryPCMEMCALPi02760GeV->Get("TrueMassPi0");
    histoPi0TrueMass2760GeV->Scale(1000);
    TH1D*       histoPi0TrueFWHM2760GeV                             = (TH1D*)directoryPCMEMCALPi02760GeV->Get("TrueFWHMPi0MeV");
    TH1D*       histoPi0Efficiency2760GeV                           = (TH1D*)directoryPCMEMCALPi02760GeV->Get("EfficiencyPi0");
    TH1D*       histoPi0Acceptance2760GeV                           = (TH1D*)directoryPCMEMCALPi02760GeV->Get("AcceptancePi0");
    
    TH1D*       histoPi0AcceptanceTimesEff2760GeV                   = (TH1D*)histoPi0Efficiency2760GeV->Clone("histoPi0AcceptanceTimesEff2760GeV");
    histoPi0AcceptanceTimesEff2760GeV->Multiply(histoPi0Acceptance2760GeV);
    
    TDirectory* directoryPCMEMCALEta2760GeV                         = (TDirectory*)filePCMEMCALpp2760GeV->Get("Eta2.76TeV"); 
    TH1D*       histoPCMEMCALYieldEta2760GeV                        = (TH1D*)directoryPCMEMCALEta2760GeV->Get(nameHistoPCMEMCALEta.Data());
    TH1D*       histoPCMEMCALRAWYieldEta2760GeV                     = (TH1D*)directoryPCMEMCALEta2760GeV->Get("RAWYieldPerEventsEta");
    TH1D*       histoEtaInputFullReweighted2760GeV                  = (TH1D*)directoryPCMEMCALEta2760GeV->Get("Eta_Input_Reweighted");
    TH1D*       histoEtaInputFull2760GeV                            = (TH1D*)directoryPCMEMCALEta2760GeV->Get("Eta_Input");
    TH1D*       histoEtaInputFullAddSigReweighted2760GeV            = (TH1D*)directoryPCMEMCALEta2760GeV->Get("Eta_Input_Reweighted_AddedSig");
    TH1D*       histoEtaInputFullAddSig2760GeV                      = (TH1D*)directoryPCMEMCALEta2760GeV->Get("Eta_Input_AddedSig");
//     TH1D*        histoEtaMass2760GeV                                 = (TH1D*)directoryPCMEMCALEta2760GeV->Get("MassEta");
//     TH1D*        histoEtaFWHM2760GeV                                 = (TH1D*)directoryPCMEMCALEta2760GeV->Get("FWHMEtaMeV");
//     TH1D*        histoEtaTrueMass2760GeV                             = (TH1D*)directoryPCMEMCALEta2760GeV->Get("TrueMassEta");
//     TH1D*        histoEtaTrueFWHM2760GeV                             = (TH1D*)directoryPCMEMCALEta2760GeV->Get("TrueFWHMEtaMeV");
    TH1D*       histoEtaEfficiency2760GeV                           = (TH1D*)directoryPCMEMCALEta2760GeV->Get("EfficiencyEta");
    TH1D*       histoEtaAcceptance2760GeV                           = (TH1D*)directoryPCMEMCALEta2760GeV->Get("AcceptanceEta");

    // Pythia 8
    TFile*      filePi0PCMEMCAL2760GeV_LHC12f1a_PYT8_wSDD                   = new TFile("0thIteration/0000311_00200009327000008250400000_10000053032230000_01631031000000_LHC12f1a/2.76TeV/Pi0_data_GammaConvV1Correction_0000311_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD         = (TH1D*)filePi0PCMEMCAL2760GeV_LHC12f1a_PYT8_wSDD->Get("MCYield_Meson_oldBinWOWeights");
    TH1D*       histoPi0Efficiency_LHC12f1a_PYT8_2760GeV_wSDD               = (TH1D*)filePi0PCMEMCAL2760GeV_LHC12f1a_PYT8_wSDD->Get("TrueMesonEffiPt");
    TFile*      fileEtaPCMEMCAL2760GeV_LHC12f1a_PYT8_wSDD                   = new TFile("0thIteration/0000311_00200009327000008250400000_10000053032230000_01631031000000_LHC12f1a/2.76TeV/Eta_data_GammaConvV1Correction_0000311_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD         = (TH1D*)fileEtaPCMEMCAL2760GeV_LHC12f1a_PYT8_wSDD->Get("MCYield_Meson_oldBinWOWeights");
    TH1D*       histoEtaEfficiency_LHC12f1a_PYT8_2760GeV_wSDD               = (TH1D*)fileEtaPCMEMCAL2760GeV_LHC12f1a_PYT8_wSDD->Get("TrueMesonEffiPt");

    // Phojet
    TFile*      filePi0PCMEMCAL2760GeV_LHC12f1b_PHO_wSDD                    = new TFile("0thIteration/0000311_00200009327000008250400000_10000053032230000_01631031000000_LHC12f1b/2.76TeV/Pi0_data_GammaConvV1Correction_0000311_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD          = (TH1D*)filePi0PCMEMCAL2760GeV_LHC12f1b_PHO_wSDD->Get("MCYield_Meson_oldBinWOWeights");
    TH1D*       histoPi0Efficiency_LHC12f1b_PHO_2760GeV_wSDD                = (TH1D*)filePi0PCMEMCAL2760GeV_LHC12f1b_PHO_wSDD->Get("TrueMesonEffiPt");
    TFile*      fileEtaPCMEMCAL2760GeV_LHC12f1b_PHO_wSDD                    = new TFile("0thIteration/0000311_00200009327000008250400000_10000053032230000_01631031000000_LHC12f1b/2.76TeV/Eta_data_GammaConvV1Correction_0000311_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD          = (TH1D*)fileEtaPCMEMCAL2760GeV_LHC12f1b_PHO_wSDD->Get("MCYield_Meson_oldBinWOWeights");
    TH1D*       histoEtaEfficiency_LHC12f1b_PHO_2760GeV_wSDD                = (TH1D*)fileEtaPCMEMCAL2760GeV_LHC12f1b_PHO_wSDD->Get("TrueMesonEffiPt");

    // Pythia 8, LHC12i3
    TFile*      filePi0PCMEMCAL2760GeV_LHC12i3_PYT8_wSDD                    = new TFile("0thIteration/0000311_00200009327000008250400000_10000053032230000_01631031000000_LHC12i3/2.76TeV/Pi0_data_GammaConvV1Correction_0000311_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD          = (TH1D*)filePi0PCMEMCAL2760GeV_LHC12i3_PYT8_wSDD->Get("MCYield_Meson_oldBinWOWeights");
    TH1D*       histoPi0Efficiency_LHC12i3_PYT8_2760GeV_wSDD                = (TH1D*)filePi0PCMEMCAL2760GeV_LHC12i3_PYT8_wSDD->Get("TrueMesonEffiPt");
    TFile*      fileEtaPCMEMCAL2760GeV_LHC12i3_PYT8_wSDD                    = new TFile("0thIteration/0000311_00200009327000008250400000_10000053032230000_01631031000000_LHC12i3/2.76TeV/Eta_data_GammaConvV1Correction_0000311_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD          = (TH1D*)fileEtaPCMEMCAL2760GeV_LHC12i3_PYT8_wSDD->Get("MCYield_Meson_oldBinWOWeights");
    TH1D*       histoEtaEfficiency_LHC12i3_PYT8_2760GeV_wSDD                = (TH1D*)fileEtaPCMEMCAL2760GeV_LHC12i3_PYT8_wSDD->Get("TrueMesonEffiPt");

    // Pythia 8, LHC12i3 added signals
    TFile*      filePi0PCMEMCAL2760GeV_LHC12i3_PYT8Addsig_wSDD              = new TFile("0thIteration/0000312_00200009327000008250400000_10000053032230000_01631031000000_LHC12i3/2.76TeV/Pi0_data_GammaConvV1Correction_0000312_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoPi0InputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD    = (TH1D*)filePi0PCMEMCAL2760GeV_LHC12i3_PYT8Addsig_wSDD->Get("MCYield_Meson_oldBinWOWeights");
    TH1D*       histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV_wSDD          = (TH1D*)filePi0PCMEMCAL2760GeV_LHC12i3_PYT8Addsig_wSDD->Get("TrueMesonEffiPt");
    TFile*      fileEtaPCMEMCAL2760GeV_LHC12i3_PYT8Addsig_wSDD              = new TFile("0thIteration/0000312_00200009327000008250400000_10000053032230000_01631031000000_LHC12i3/2.76TeV/Eta_data_GammaConvV1Correction_0000312_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoEtaInputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD    = (TH1D*)fileEtaPCMEMCAL2760GeV_LHC12i3_PYT8Addsig_wSDD->Get("MCYield_Meson_oldBinWOWeights");
    TH1D*       histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV_wSDD          = (TH1D*)fileEtaPCMEMCAL2760GeV_LHC12i3_PYT8Addsig_wSDD->Get("TrueMesonEffiPt");

    // Pythia 8
    TFile*      filePi0PCMEMCAL2760GeV_LHC12f1a_PYT8_wSDDRW                 = new TFile("0000311_00200009327000008250400000_10000053032230000_01631031000000_LHC12f1a/2.76TeV/Pi0_data_GammaConvV1Correction_0000311_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoPi0Efficiency_LHC12f1a_PYT8_2760GeV_wSDDRW             = (TH1D*)filePi0PCMEMCAL2760GeV_LHC12f1a_PYT8_wSDDRW->Get("TrueMesonEffiPt");
    TFile*      fileEtaPCMEMCAL2760GeV_LHC12f1a_PYT8_wSDDRW                 = new TFile("0000311_00200009327000008250400000_10000053032230000_01631031000000_LHC12f1a/2.76TeV/Eta_data_GammaConvV1Correction_0000311_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoEtaEfficiency_LHC12f1a_PYT8_2760GeV_wSDDRW             = (TH1D*)fileEtaPCMEMCAL2760GeV_LHC12f1a_PYT8_wSDDRW->Get("TrueMesonEffiPt");

    // Phojet
    TFile*      filePi0PCMEMCAL2760GeV_LHC12f1b_PHO_wSDDRW                  = new TFile("0000311_00200009327000008250400000_10000053032230000_01631031000000_LHC12f1b/2.76TeV/Pi0_data_GammaConvV1Correction_0000311_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoPi0Efficiency_LHC12f1b_PHO_2760GeV_wSDDRW              = (TH1D*)filePi0PCMEMCAL2760GeV_LHC12f1b_PHO_wSDDRW->Get("TrueMesonEffiPt");
    TFile*      fileEtaPCMEMCAL2760GeV_LHC12f1b_PHO_wSDDRW                  = new TFile("0000311_00200009327000008250400000_10000053032230000_01631031000000_LHC12f1b/2.76TeV/Eta_data_GammaConvV1Correction_0000311_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoEtaEfficiency_LHC12f1b_PHO_2760GeV_wSDDRW              = (TH1D*)fileEtaPCMEMCAL2760GeV_LHC12f1b_PHO_wSDDRW->Get("TrueMesonEffiPt");

    // Pythia 8, LHC12i3
    TFile*      filePi0PCMEMCAL2760GeV_LHC12i3_PYT8_wSDDRW                  = new TFile("0000311_00200009327000008250400000_10000053032230000_01631031000000_LHC12i3/2.76TeV/Pi0_data_GammaConvV1Correction_0000311_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoPi0Efficiency_LHC12i3_PYT8_2760GeV_wSDDRW              = (TH1D*)filePi0PCMEMCAL2760GeV_LHC12i3_PYT8_wSDDRW->Get("TrueMesonEffiPt");
    TFile*      fileEtaPCMEMCAL2760GeV_LHC12i3_PYT8_wSDDRW                  = new TFile("0000311_00200009327000008250400000_10000053032230000_01631031000000_LHC12i3/2.76TeV/Eta_data_GammaConvV1Correction_0000311_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoEtaEfficiency_LHC12i3_PYT8_2760GeV_wSDDRW              = (TH1D*)fileEtaPCMEMCAL2760GeV_LHC12i3_PYT8_wSDDRW->Get("TrueMesonEffiPt");

    // Pythia 8, LHC12i3 added signals
    TFile*      filePi0PCMEMCAL2760GeV_LHC12i3_PYT8Addsig_wSDDRW            = new TFile("0000312_00200009327000008250400000_10000053032230000_01631031000000_LHC12i3/2.76TeV/Pi0_data_GammaConvV1Correction_0000312_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV_wSDDRW        = (TH1D*)filePi0PCMEMCAL2760GeV_LHC12i3_PYT8Addsig_wSDDRW->Get("TrueMesonEffiPt");
    TFile*      fileEtaPCMEMCAL2760GeV_LHC12i3_PYT8Addsig_wSDDRW            = new TFile("0000312_00200009327000008250400000_10000053032230000_01631031000000_LHC12i3/2.76TeV/Eta_data_GammaConvV1Correction_0000312_00200009327000008250400000_10000053032230000_01631031000000.root");
    TH1D*       histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV_wSDDRW        = (TH1D*)fileEtaPCMEMCAL2760GeV_LHC12i3_PYT8Addsig_wSDDRW->Get("TrueMesonEffiPt");
    
    //    **********************************************************************************************************************
    //    ****************************************Pi0 Spectra 2.76TeV compared to MC********************************************
    //    **********************************************************************************************************************
    TCanvas* canvasPi0Spectra2760GeV                                 = new TCanvas("canvasPi0Spectra2760GeV", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasPi0Spectra2760GeV,  0.13, 0.01, 0.015, 0.08);
    canvasPi0Spectra2760GeV->SetLogy();
    canvasPi0Spectra2760GeV->SetLogx();
    
    TH2F * histo2DPi0Spectra2760GeV;
    histo2DPi0Spectra2760GeV                                             = new TH2F("histo2DPi0Spectra2760GeV", "histo2DPi0Spectra2760GeV",1000, 0.23, 20., 1000, 1e-8, 1e1 );
    SetStyleHistoTH2ForGraphs( histo2DPi0Spectra2760GeV, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 
                            0.03, 0.04, 0.03, 0.04, 0.83, 1.4);
    histo2DPi0Spectra2760GeV->GetXaxis()->SetLabelOffset(-0.01);
    histo2DPi0Spectra2760GeV->GetYaxis()->SetLabelOffset(0.01);
    histo2DPi0Spectra2760GeV->DrawCopy(); 

    DrawGammaSetMarker(histoPCMEMCALYieldPi02760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, color2760GeV , color2760GeV);
    histoPCMEMCALYieldPi02760GeV->Draw("p,same,e1");
    
    SetStyleHisto(histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD, 1., lineStyleMCA, colorMCPythiaPP2760GeV);  
    histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD->Draw("same,hist,c");    
    SetStyleHisto(histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD, 1., lineStyleMCB, colorMCPhojetPP2760GeV);  
    histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD->Draw("same,hist,c");    
//     SetStyleHisto(histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD, 1., lineStyleMCC, colorMCPythia2PP2760GeV);  
//     histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD->Draw("same,hist,c");    
//     SetStyleHisto(histoPi0InputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD, 1., lineStyleMCAddSig, colorMCPythiaAddSigPP2760GeV);  
//     histoPi0InputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD->Draw("same,hist,c");    

    TLatex *labelSpectraPi0Label                                     = new TLatex(0.75,0.92,"#pi^{0} #rightarrow #gamma #gamma(e^{+}e^{-})");
    SetStyleTLatex( labelSpectraPi0Label, 0.035  ,4);
    labelSpectraPi0Label->Draw();
    
    TLegend* legendSpectraPi02760GeV                                 = new TLegend(0.16,0.09,0.73,0.3);
    legendSpectraPi02760GeV->SetFillColor(0);
    legendSpectraPi02760GeV->SetFillStyle(0);
    legendSpectraPi02760GeV->SetLineColor(0);
    legendSpectraPi02760GeV->SetTextSize(0.025);
    legendSpectraPi02760GeV->SetMargin(0.2);
    legendSpectraPi02760GeV->AddEntry(histoPCMEMCALYieldPi02760GeV,"PCM-EMCal, stat","pe");
    legendSpectraPi02760GeV->AddEntry(histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD,"Phojet, LHC12f1b","l");
    legendSpectraPi02760GeV->AddEntry(histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD,"Pythia 8, LHC12f1a","l");
//     legendSpectraPi02760GeV->AddEntry(histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD,"Pythia 8, LHC12i3","l");
//     legendSpectraPi02760GeV->AddEntry(histoPi0InputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD,"Pythia 8, LHC12i3, added Signals","l");
    legendSpectraPi02760GeV->Draw();
    
    canvasPi0Spectra2760GeV->Update();
    canvasPi0Spectra2760GeV->Print(Form("%s/Pi0_Spectra_2760GeV_wSDD.%s",outputDir.Data(),suffix.Data()));
// 
//     //    **********************************************************************************************************************
//     //    ****************************************Fit Pi0 Spectra 2.76TeV and plot *********************************************
//     //    **********************************************************************************************************************
//     
//     // fit spectrum with Tsallis function
//     TF1* fitInvYieldDataPi0Comb2760GeV                             = FitObject("l","fitInvYieldDataPi0Comb2760GeV","Pi0");
//     histoPCMEMCALYieldPi02760GeV->Fit(fitInvYieldDataPi0Comb2760GeV,"QNRMEI+","",0.4,10.);
//     fitInvYieldDataPi0Comb2760GeV->SetRange(0.1,20.);
//     // print fit result to shell
//     cout << WriteParameterToFile(fitInvYieldDataPi0Comb2760GeV)<< endl;
// 
//     // plot result with input
//     canvasPi0Spectra2760GeV->cd();
//     histo2DPi0Spectra2760GeV->DrawCopy(); 
// 
//     DrawGammaSetMarker(histoPCMEMCALYieldPi02760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, color2760GeV , color2760GeV);
//     histoPCMEMCALYieldPi02760GeV->Draw("p,same,e1");
// 
//     
//     fitInvYieldDataPi0Comb2760GeV->SetLineColor(color2760GeV);
//     fitInvYieldDataPi0Comb2760GeV->Draw("same");
//     
//     SetStyleHisto(histoPi0InputFullReweighted2760GeV, 1., lineStyleMCAddSig, kBlue+1);  
//     histoPi0InputFullReweighted2760GeV->Draw("same,hist,c");    
//     SetStyleHisto(histoPi0InputFullAddSigReweighted2760GeV, 1., lineStyleMCAddSig, kBlue+2);  
//     histoPi0InputFullAddSigReweighted2760GeV->Draw("same,hist,c");    
// 
//     labelSpectraPi0Label->Draw();
//     
//     canvasPi0Spectra2760GeV->Print(Form("%s/Pi0_Spectra_WithFit_2760GeV.%s",outputDir.Data(),suffix.Data()));
//     
//     // **********************************************************************************************************************
//     // ******************************* Ratio of data to fit and MC input to fit for pi0 2.76 TeV ****************************
//     // **********************************************************************************************************************
//     // definition of canvas for ratios
//     TCanvas* canvasRatioToFit = new TCanvas("canvasRatioToFit","",1550,1200);  // gives the page size
//     DrawGammaCanvasSettings( canvasRatioToFit,  0.08, 0.015, 0.015, 0.08);
//     canvasRatioToFit->SetGridx(0);
//     canvasRatioToFit->SetGridy(0);
// 
//     // Calculation of ratio histograms
//     TH1D* histoRatioPi0DatatoFit2760GeV         = CalculateHistoRatioToFit (histoPCMEMCALYieldPi02760GeV, fitInvYieldDataPi0Comb2760GeV,kTRUE);
//     TH1D* histoRatioPi0MCtoDataFit2760GeV       = CalculateHistoRatioToFit (histoPi0InputFullReweighted2760GeV, fitInvYieldDataPi0Comb2760GeV,kTRUE);
//     TH1D* histoRatioPi0MCAddSigtoDataFit2760GeV = CalculateHistoRatioToFit (histoPi0InputFullAddSigReweighted2760GeV, fitInvYieldDataPi0Comb2760GeV,kTRUE);
//     TH1D* histoRatioPi0MCUnweightedtoDataFit2760GeV = NULL;
//     if (histoPi0InputFull2760GeV) histoRatioPi0MCUnweightedtoDataFit2760GeV = CalculateHistoRatioToFit (histoPi0InputFull2760GeV, fitInvYieldDataPi0Comb2760GeV,kTRUE);
//     if (histoRatioPi0MCUnweightedtoDataFit2760GeV) SetStyleHisto(histoRatioPi0MCUnweightedtoDataFit2760GeV, 2, lineStyleMCB, 807 );
//     canvasRatioToFit->cd();
//     SetStyleHisto(histoRatioPi0MCtoDataFit2760GeV, 2, lineStyleMCA, kRed+2 );
//     SetStyleHisto(histoRatioPi0MCAddSigtoDataFit2760GeV, 2, lineStyleMCAddSig, kBlue+2 );
//     DrawGammaSetMarker(histoRatioPi0DatatoFit2760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, kBlack , kBlack);
//     DrawAutoGammaMesonHistos( histoRatioPi0DatatoFit2760GeV,
//                 "", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
//                 kFALSE, 1.5, 0, kTRUE,
//                 kTRUE, 0, 2.1,
//                 kTRUE, 0.,histoRatioPi0DatatoFit2760GeV->GetXaxis()->GetBinUpEdge(histoRatioPi0DatatoFit2760GeV->GetNbinsX()));
//     histoRatioPi0DatatoFit2760GeV->GetYaxis()->SetTitleOffset(0.9);
//     histoRatioPi0DatatoFit2760GeV->Draw("e,p");  
//     if (runDrawReweighted) histoRatioPi0MCtoDataFit2760GeV->Draw("same,hist,l");  
//     if (runDrawReweighted) histoRatioPi0MCAddSigtoDataFit2760GeV->Draw("same,hist,l");  
//     if (histoRatioPi0MCUnweightedtoDataFit2760GeV) histoRatioPi0MCUnweightedtoDataFit2760GeV->Draw("same,hist,l");  
// 
//     TLatex *labelSpectraPi0LabelRatio                                    = new TLatex(0.65,0.9,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
//     SetStyleTLatex( labelSpectraPi0LabelRatio, 0.035  ,4);
//     labelSpectraPi0LabelRatio->Draw();
//     TLatex *labelEnergy2760GeVRatio                                      = new TLatex(0.65,0.86,collisionSystemPP2760GeV.Data());
//     SetStyleTLatex( labelEnergy2760GeVRatio, 0.035  ,4);
//     labelEnergy2760GeVRatio->Draw();
//     
//     TLegend* legendRatioPi02760GeV = new TLegend(0.11,0.12,0.4,0.30);
//     legendRatioPi02760GeV->SetFillColor(0);
//     legendRatioPi02760GeV->SetLineColor(0);
//     legendRatioPi02760GeV->SetTextSize(0.035);
//     legendRatioPi02760GeV->SetMargin(0.2);
//     legendRatioPi02760GeV->AddEntry(histoRatioPi0DatatoFit2760GeV,"Data/Tsallis fit to Data (0.4 <#it{p}_{T}<10)","p");
//     if (runDrawReweighted) legendRatioPi02760GeV->AddEntry(histoRatioPi0MCtoDataFit2760GeV,Form("MC weighted %s/Tsallis  fit to Data (0.4 <#it{p}_{T}<10)",stringIterationNumber.Data()),"l");
//     if (runDrawReweighted) legendRatioPi02760GeV->AddEntry(histoRatioPi0MCAddSigtoDataFit2760GeV,Form("MC add Sig weighted %s/Tsallis  fit to Data (0.4 <#it{p}_{T}<10)", 
//                                                            stringIterationNumber.Data()), "l");
//     if (histoRatioPi0MCUnweightedtoDataFit2760GeV) legendRatioPi02760GeV->AddEntry(histoRatioPi0MCUnweightedtoDataFit2760GeV,"MC/Tsallis  fit to Data (0.4 <#it{p}_{T}<10)","l");
//     legendRatioPi02760GeV->Draw();
//     
//     DrawGammaLines(0., histoRatioPi0DatatoFit2760GeV->GetXaxis()->GetBinUpEdge(histoRatioPi0DatatoFit2760GeV->GetNbinsX()) ,1., 1.,0.1);
//     canvasRatioToFit->Update();
//     canvasRatioToFit->SaveAs(Form("%s/Pi0_RatioToDataFit_2760GeV.%s",outputDir.Data(),suffix.Data()));
    
    //    **********************************************************************************************************************
    //    ***********************Compare Efficiencies for pi0 at 2.76TeV for different MC without weighting ********************
    //    **********************************************************************************************************************
    TCanvas* canvasPi0Efficiencies2760GeVDiffMC                     = new TCanvas("canvasPi0Efficiencies2760GeVDiffMC", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasPi0Efficiencies2760GeVDiffMC,  0.09, 0.01, 0.015, 0.08);
    canvasPi0Efficiencies2760GeVDiffMC->SetLogy();
//     canvasPi0Efficiencies2760GeVDiffMC->SetLogx();
    
    TH2F * histo2DPi0Effi2760GeV;
    histo2DPi0Effi2760GeV                                           = new TH2F("histo2DPi0Effi2760GeV", "histo2DPi0Effi2760GeV",1000, 0., 9.5, 1000, 1e-5, 1e-1 );
    SetStyleHistoTH2ForGraphs( histo2DPi0Effi2760GeV, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{#pi^{0}}", 
                            0.03, 0.04, 0.03, 0.04, 0.83, 1.05);
    histo2DPi0Effi2760GeV->GetYaxis()->SetLabelOffset(0.01);
    histo2DPi0Effi2760GeV->DrawCopy(); 

    canvasPi0Efficiencies2760GeVDiffMC->cd();
    histo2DPi0Effi2760GeV->DrawCopy(); 

    DrawGammaSetMarker(histoPi0Efficiency_LHC12f1a_PYT8_2760GeV_wSDD, markerStyleMCA, markerSizePP2760GeV, colorMCPythiaPP2760GeV , colorMCPythiaPP2760GeV);
    histoPi0Efficiency_LHC12f1a_PYT8_2760GeV_wSDD->Draw("p,same,e1");

//     DrawGammaSetMarker(histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV_wSDD, markerStyleMCAddSig, markerSizePP2760GeV, colorMCPythiaAddSigPP2760GeV , colorMCPythiaAddSigPP2760GeV);
//     histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV_wSDD->Draw("p,same,e1");

//     DrawGammaSetMarker(histoPi0Efficiency_LHC12i3_PYT8_2760GeV_wSDD, markerStyleMCC, markerSizePP2760GeV, colorMCPythia2PP2760GeV , colorMCPythia2PP2760GeV);
//     histoPi0Efficiency_LHC12i3_PYT8_2760GeV_wSDD->Draw("p,same,e1");
    
    DrawGammaSetMarker(histoPi0Efficiency_LHC12f1b_PHO_2760GeV_wSDD, markerStyleMCB, markerSizePP2760GeV, colorMCPhojetPP2760GeV , colorMCPhojetPP2760GeV);
    histoPi0Efficiency_LHC12f1b_PHO_2760GeV_wSDD->Draw("p,same,e1");

    TLegend* legendEfficiencyPi02760GeV                             = new TLegend(0.56,0.09,0.93,0.3);
    legendEfficiencyPi02760GeV->SetFillColor(0);
    legendEfficiencyPi02760GeV->SetFillStyle(0);
    legendEfficiencyPi02760GeV->SetLineColor(0);
    legendEfficiencyPi02760GeV->SetTextSize(0.025);
    legendEfficiencyPi02760GeV->SetMargin(0.2);
    legendEfficiencyPi02760GeV->AddEntry(histoPi0Efficiency_LHC12f1b_PHO_2760GeV_wSDD,"Phojet, LHC12f1b","p");
    legendEfficiencyPi02760GeV->AddEntry(histoPi0Efficiency_LHC12f1a_PYT8_2760GeV_wSDD,"Pythia 8, LHC12f1a","p");
//     legendEfficiencyPi02760GeV->AddEntry(histoPi0Efficiency_LHC12i3_PYT8_2760GeV_wSDD,"Pythia 8, LHC12i3","p");
//     legendEfficiencyPi02760GeV->AddEntry(histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV_wSDD,"Pythia 8, LHC12i3, added Signals","p");
    legendEfficiencyPi02760GeV->Draw();
    
    canvasPi0Efficiencies2760GeVDiffMC->Update();
    canvasPi0Efficiencies2760GeVDiffMC->Print(Form("%s/Pi0_Efficiency_WSDD_2760GeV.%s",outputDir.Data(),suffix.Data()));
    
    canvasPi0Efficiencies2760GeVDiffMC->SetLogy(0);
    canvasPi0Efficiencies2760GeVDiffMC->SetTopMargin(0.035);
    
    histo2DPi0Effi2760GeV->GetYaxis()->SetRangeUser(1e-5,2E-2);
    histo2DPi0Effi2760GeV->DrawCopy(); 
    
    histoPi0Efficiency_LHC12f1a_PYT8_2760GeV_wSDD->Draw("p,same,e1");
//     histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV_wSDD->Draw("p,same,e1");
//     histoPi0Efficiency_LHC12i3_PYT8_2760GeV_wSDD->Draw("p,same,e1");
    histoPi0Efficiency_LHC12f1b_PHO_2760GeV_wSDD->Draw("p,same,e1");
    
    legendEfficiencyPi02760GeV->SetX1NDC(0.1);
    legendEfficiencyPi02760GeV->SetX2NDC(0.1+0.37);
    legendEfficiencyPi02760GeV->SetY1NDC(0.94-0.21);
    legendEfficiencyPi02760GeV->SetY2NDC(0.94);
    legendEfficiencyPi02760GeV->Draw();
    
    canvasPi0Efficiencies2760GeVDiffMC->Update();
    canvasPi0Efficiencies2760GeVDiffMC->Print(Form("%s/Pi0_Efficiency_WSDD_2760GeV_LinY.%s",outputDir.Data(),suffix.Data()));


    //    **********************************************************************************************************************
    //    ***********************Compare Efficiencies for pi0 at 2.76TeV for different MC with weighting ***********************
    //    **********************************************************************************************************************

    canvasPi0Efficiencies2760GeVDiffMC->cd();
    canvasPi0Efficiencies2760GeVDiffMC->SetLogy(1);
    histo2DPi0Effi2760GeV->DrawCopy(); 

    DrawGammaSetMarker(histoPi0Efficiency_LHC12f1a_PYT8_2760GeV_wSDDRW, markerStyleMCA, markerSizePP2760GeV, colorMCPythiaPP2760GeV , colorMCPythiaPP2760GeV);
    histoPi0Efficiency_LHC12f1a_PYT8_2760GeV_wSDDRW->Draw("p,same,e1");

//     DrawGammaSetMarker(histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV_wSDDRW, markerStyleMCAddSig, markerSizePP2760GeV, colorMCPythiaAddSigPP2760GeV , colorMCPythiaAddSigPP2760GeV);
//     histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV_wSDDRW->Draw("p,same,e1");

//     DrawGammaSetMarker(histoPi0Efficiency_LHC12i3_PYT8_2760GeV_wSDDRW, markerStyleMCC, markerSizePP2760GeV, colorMCPythia2PP2760GeV , colorMCPythia2PP2760GeV);
//     histoPi0Efficiency_LHC12i3_PYT8_2760GeV_wSDDRW->Draw("p,same,e1");
    
    DrawGammaSetMarker(histoPi0Efficiency_LHC12f1b_PHO_2760GeV_wSDDRW, markerStyleMCB, markerSizePP2760GeV, colorMCPhojetPP2760GeV , colorMCPhojetPP2760GeV);
    histoPi0Efficiency_LHC12f1b_PHO_2760GeV_wSDDRW->Draw("p,same,e1");

    DrawGammaSetMarker(histoPi0Efficiency2760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, color2760GeV , color2760GeV);
    histoPi0Efficiency2760GeV->Draw("p,same,e1");
    
    legendEfficiencyPi02760GeV->AddEntry(histoPi0Efficiency2760GeV,"final merged efficiency","p");
    legendEfficiencyPi02760GeV->SetX1NDC(0.56);
    legendEfficiencyPi02760GeV->SetX2NDC(0.56+0.37);
    legendEfficiencyPi02760GeV->SetY1NDC(0.09-0.21);
    legendEfficiencyPi02760GeV->SetY2NDC(0.09);
    legendEfficiencyPi02760GeV->Draw();
    
    canvasPi0Efficiencies2760GeVDiffMC->Update();
    canvasPi0Efficiencies2760GeVDiffMC->Print(Form("%s/Pi0_Efficiency_WSDD_2760GeV_reweighted.%s",outputDir.Data(),suffix.Data()));
    
    canvasPi0Efficiencies2760GeVDiffMC->SetLogy(0);
    canvasPi0Efficiencies2760GeVDiffMC->SetTopMargin(0.035);
    
    histo2DPi0Effi2760GeV->GetYaxis()->SetRangeUser(1e-5,2E-2);
    histo2DPi0Effi2760GeV->DrawCopy(); 
    
    histoPi0Efficiency_LHC12f1a_PYT8_2760GeV_wSDDRW->Draw("p,same,e1");
//     histoPi0Efficiency_LHC12i3_PYT8addSig_2760GeV_wSDDRW->Draw("p,same,e1");
//     histoPi0Efficiency_LHC12i3_PYT8_2760GeV_wSDDRW->Draw("p,same,e1");
    histoPi0Efficiency_LHC12f1b_PHO_2760GeV_wSDDRW->Draw("p,same,e1");
    histoPi0Efficiency2760GeV->Draw("p,same,e1");
    
    legendEfficiencyPi02760GeV->SetX1NDC(0.1);
    legendEfficiencyPi02760GeV->SetX2NDC(0.1+0.37);
    legendEfficiencyPi02760GeV->SetY1NDC(0.94-0.21);
    legendEfficiencyPi02760GeV->SetY2NDC(0.94);
    legendEfficiencyPi02760GeV->Draw();
    
    canvasPi0Efficiencies2760GeVDiffMC->Update();
    canvasPi0Efficiencies2760GeVDiffMC->Print(Form("%s/Pi0_Efficiency_WSDD_2760GeV_reweighted_LinY.%s",outputDir.Data(),suffix.Data()));
    
    //    **********************************************************************************************************************
    //    ****************************************Eta Spectra 2.76TeV compared to MC********************************************
    //    **********************************************************************************************************************    
    TCanvas* canvasEtaSpectra2760GeV                                 = new TCanvas("canvasEtaSpectra2760GeV", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasEtaSpectra2760GeV,  0.13, 0.01, 0.015, 0.08);
    canvasEtaSpectra2760GeV->SetLogy();
    canvasEtaSpectra2760GeV->SetLogx();
    
    TH2F * histo2DEtaSpectraAll;
    histo2DEtaSpectraAll                                             = new TH2F("histo2DEtaSpectraAll", "histo2DEtaSpectraAll",1000, 0.23, 20., 1000, 1e-8, 1e1 );
    SetStyleHistoTH2ForGraphs( histo2DEtaSpectraAll, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 
                            0.03, 0.04, 0.03, 0.04, 0.83, 1.4);
    histo2DEtaSpectraAll->GetXaxis()->SetLabelOffset(-0.01);
    histo2DEtaSpectraAll->GetYaxis()->SetLabelOffset(0.01);
    histo2DEtaSpectraAll->DrawCopy(); 
    
    DrawGammaSetMarker(histoPCMEMCALYieldEta2760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, color2760GeV , color2760GeV);
    histoPCMEMCALYieldEta2760GeV->Draw("p,same,e1");

    SetStyleHisto(histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD, 1., lineStyleMCA, colorMCPythiaPP2760GeV);  
    histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD->Draw("same,hist,c");    
    SetStyleHisto(histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD, 1., lineStyleMCB, colorMCPhojetPP2760GeV);  
    histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD->Draw("same,hist,c");    
//     SetStyleHisto(histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD, 1., lineStyleMCC, colorMCPythia2PP2760GeV);  
//     histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD->Draw("same,hist,c");    
//     SetStyleHisto(histoEtaInputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD, 1., lineStyleMCAddSig, colorMCPythiaAddSigPP2760GeV);  
//     histoEtaInputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD->Draw("same,hist,c");    
    
    TLatex *labelSpectraEtaLabel                             = new TLatex(0.55,0.92,"#eta #rightarrow #gamma #gamma(e^{+}e^{-})");
    SetStyleTLatex( labelSpectraEtaLabel, 0.035  ,4);
    labelSpectraEtaLabel->Draw();
    
    TLegend* legendSpectraEta2760GeV                         = new TLegend(0.16,0.09,0.73,0.3);
    legendSpectraEta2760GeV->SetFillColor(0);
    legendSpectraEta2760GeV->SetLineColor(0);
    legendSpectraEta2760GeV->SetTextSize(0.025);
    legendSpectraEta2760GeV->SetMargin(0.2);
    legendSpectraEta2760GeV->AddEntry(histoPCMEMCALYieldEta2760GeV,"PCM-EMCal, stat","pe");
    legendSpectraEta2760GeV->AddEntry(histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD,"Phojet, LHC12f1b","l");
    legendSpectraEta2760GeV->AddEntry(histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD,"Pythia 8, LHC12f1a","l");
//     legendSpectraEta2760GeV->AddEntry(histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD,"Pythia 8, LHC12i3","l");
//     legendSpectraEta2760GeV->AddEntry(histoEtaInputFullReweighted2760GeV,"full MC reweighted","l");
//     legendSpectraEta2760GeV->AddEntry(histoEtaInputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD,"Pythia 8, LHC12i3, added Signals","l");
//     legendSpectraEta2760GeV->AddEntry(histoEtaInputFullAddSigReweighted2760GeV,"added Signals, reweighted","l");
    legendSpectraEta2760GeV->Draw();
    
    canvasEtaSpectra2760GeV->Update();
    canvasEtaSpectra2760GeV->Print(Form("%s/Eta_Spectra_2760GeV_wSDD.%s",outputDir.Data(),suffix.Data()));
// 
//     //    **********************************************************************************************************************
//     //    ****************************************Fit Pi0 Spectra 2.76TeV and plot *********************************************
//     //    **********************************************************************************************************************
//     
//     // fit spectrum with Tsallis function - "n" fixed by result from pi0 at same energy    
//     TF1* fitInvYieldDataEtaComb2760GeV                                 = FitObject("l","fitInvYieldDataEtaComb2760GeV","Eta");
// //     fitInvYieldDataEtaComb2760GeV->FixParameter(1,fitInvYieldDataPi0Comb2760GeV->GetParameter(1));
// //     fitInvYieldDataEtaComb2760GeV->SetParameter(1,fitInvYieldDataPi0Comb2760GeV->GetParameter(1));
// //     fitInvYieldDataEtaComb2760GeV->SetParLimits(1,fitInvYieldDataPi0Comb2760GeV->GetParameter(1)*0.9,fitInvYieldDataPi0Comb2760GeV->GetParameter(1)*1.1);
//     histoPCMEMCALYieldEta2760GeV->Fit(fitInvYieldDataEtaComb2760GeV,"QNRMEI+","",0.6,6.);
//     fitInvYieldDataEtaComb2760GeV->SetRange(0.1,20.);
//     // print result to shell
//     cout << WriteParameterToFile(fitInvYieldDataEtaComb2760GeV)<< endl;    
//     
//     // plot result with input
//     canvasEtaSpectra2760GeV->cd();
//     histo2DEtaSpectraAll->DrawCopy(); 
// 
//     DrawGammaSetMarker(histoPCMEMCALYieldEta2760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, color2760GeV , color2760GeV);
//     histoPCMEMCALYieldEta2760GeV->Draw("p,same,e1");
// 
//     fitInvYieldDataEtaComb2760GeV->SetLineColor(color2760GeV);
//     fitInvYieldDataEtaComb2760GeV->Draw("same");
//     SetStyleHisto(histoEtaInputFullReweighted2760GeV, 1., lineStyleMCAddSig, kBlue+1);  
//     histoEtaInputFullReweighted2760GeV->Draw("same,hist,c");    
//     SetStyleHisto(histoEtaInputFullAddSigReweighted2760GeV, 1., lineStyleMCAddSig, kBlue+2);  
//     histoEtaInputFullAddSigReweighted2760GeV->Draw("same,hist,c");    
// 
//     
//     labelSpectraEtaLabel->Draw();
//     
//     canvasEtaSpectra2760GeV->Print(Form("%s/Eta_Spectra_WithFit_2760GeV.%s",outputDir.Data(),suffix.Data()));
// 
//     // **********************************************************************************************************************
//     // ******************************* Ratio of data to fit and MC input to fit for Eta 2.76 TeV ****************************
//     // **********************************************************************************************************************
//     
//     // Calculation of ratio histograms
//     TH1D* histoRatioEtaDatatoFit2760GeV         = CalculateHistoRatioToFit (histoPCMEMCALYieldEta2760GeV, fitInvYieldDataEtaComb2760GeV,kTRUE);
//     TH1D* histoRatioEtaMCtoDataFit2760GeV         = CalculateHistoRatioToFit (histoEtaInputFullReweighted2760GeV, fitInvYieldDataEtaComb2760GeV,kTRUE);
//     TH1D* histoRatioEtaMCAddSigtoDataFit2760GeV = CalculateHistoRatioToFit (histoEtaInputFullAddSigReweighted2760GeV, fitInvYieldDataEtaComb2760GeV,kTRUE);
//     TH1D* histoRatioEtaMCUnweightedtoDataFit2760GeV = NULL;
//     if (histoEtaInputFull2760GeV) histoRatioEtaMCUnweightedtoDataFit2760GeV = CalculateHistoRatioToFit (histoEtaInputFull2760GeV, fitInvYieldDataEtaComb2760GeV,kTRUE);
//     if (histoRatioEtaMCUnweightedtoDataFit2760GeV) SetStyleHisto(histoRatioEtaMCUnweightedtoDataFit2760GeV, 2, lineStyleMCB, 807 );
//     // plotting ratios
//     canvasRatioToFit->cd();
//     SetStyleHisto(histoRatioEtaMCtoDataFit2760GeV, 2, lineStyleMCA, kRed+2 );
//     SetStyleHisto(histoRatioEtaMCAddSigtoDataFit2760GeV, 2, lineStyleMCAddSig, kBlue+2 );
//     DrawGammaSetMarker(histoRatioEtaDatatoFit2760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, kBlack , kBlack);
//     DrawAutoGammaMesonHistos( histoRatioEtaDatatoFit2760GeV,
//                 "", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
//                 kFALSE, 1.5, 0, kTRUE,
//                 kTRUE, 0, 2.1,
//                 kTRUE, 0.,histoRatioEtaDatatoFit2760GeV->GetXaxis()->GetBinUpEdge(histoRatioEtaDatatoFit2760GeV->GetNbinsX()));
//     histoRatioEtaDatatoFit2760GeV->GetYaxis()->SetTitleOffset(0.9);
//     histoRatioEtaDatatoFit2760GeV->Draw("e,p");  
//     if (runDrawReweighted) histoRatioEtaMCtoDataFit2760GeV->Draw("same,hist,l");  
//     if (runDrawReweighted) histoRatioEtaMCAddSigtoDataFit2760GeV->Draw("same,hist,l");  
//     if (histoRatioEtaMCUnweightedtoDataFit2760GeV) histoRatioEtaMCUnweightedtoDataFit2760GeV->Draw("same,hist,l");  
// 
//     // labeling
//     TLatex *labelSpectraEtaLabelRatio                                     = new TLatex(0.65,0.9,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
//     SetStyleTLatex( labelSpectraEtaLabelRatio, 0.035  ,4);
//     labelSpectraEtaLabelRatio->Draw();
//     labelEnergy2760GeVRatio->Draw();
//     
//     TLegend* legendRatioEta2760GeV = new TLegend(0.11,0.12,0.4,0.30);
//     legendRatioEta2760GeV->SetFillStyle(0);
//     legendRatioEta2760GeV->SetFillColor(0);
//     legendRatioEta2760GeV->SetLineColor(0);
//     legendRatioEta2760GeV->SetTextSize(0.035);
//     legendRatioEta2760GeV->SetMargin(0.2);
//     legendRatioEta2760GeV->AddEntry(histoRatioEtaDatatoFit2760GeV,"Data/Tsallis fit to Data (0.6 <#it{p}_{T}<6)","p");
//     if (runDrawReweighted) legendRatioEta2760GeV->AddEntry(histoRatioEtaMCtoDataFit2760GeV,Form("MC weighted %s/Tsallis  fit to Data (0.6 <#it{p}_{T}<6)",stringIterationNumber.Data()),"l");
//     if (runDrawReweighted) legendRatioEta2760GeV->AddEntry(histoRatioEtaMCAddSigtoDataFit2760GeV,Form("MC add Sig weighted %s/Tsallis  fit to Data (0.6 <#it{p}_{T}<6)",
//                                                            stringIterationNumber.Data()), "l");
//     if (histoRatioEtaMCUnweightedtoDataFit2760GeV) legendRatioEta2760GeV->AddEntry(histoRatioEtaMCUnweightedtoDataFit2760GeV,"MC/Tsallis  fit to Data (0.6 <#it{p}_{T}<6)","l");
//     legendRatioEta2760GeV->Draw();
//     
//     DrawGammaLines(0., histoRatioEtaDatatoFit2760GeV->GetXaxis()->GetBinUpEdge(histoRatioEtaDatatoFit2760GeV->GetNbinsX()) ,1., 1.,0.1);
//     canvasRatioToFit->Update();
//     canvasRatioToFit->SaveAs(Form("%s/Eta_RatioToDataFit_2760GeV.%s",outputDir.Data(),suffix.Data()));
//     
    
    //    **********************************************************************************************************************
    //    ******************************Compare Efficiencies for Eta at 2.76TeV for different MC *******************************
    //    **********************************************************************************************************************
    TCanvas* canvasEtaEfficiencies2760GeVDiffMC                     = new TCanvas("canvasEtaEfficiencies2760GeVDiffMC", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasEtaEfficiencies2760GeVDiffMC,  0.09, 0.01, 0.015, 0.08);
    canvasEtaEfficiencies2760GeVDiffMC->SetLogy();
    
    TH2F * histo2DEtaEffi2760GeV;
    histo2DEtaEffi2760GeV                                                     = new TH2F("histo2DEtaEffi2760GeV", "histo2DEtaEffi2760GeV",1000, 0., 6.5, 1000, 1e-5, 1e-1 );
    SetStyleHistoTH2ForGraphs( histo2DEtaEffi2760GeV, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{#eta}", 
                            0.03, 0.04, 0.03, 0.04, 0.83, 1.05);
    histo2DEtaEffi2760GeV->GetYaxis()->SetLabelOffset(0.01);
    histo2DEtaEffi2760GeV->DrawCopy(); 

    DrawGammaSetMarker(histoEtaEfficiency_LHC12f1a_PYT8_2760GeV_wSDD, markerStyleMCA, markerSizePP2760GeV, colorMCPythiaPP2760GeV , colorMCPythiaPP2760GeV);
    histoEtaEfficiency_LHC12f1a_PYT8_2760GeV_wSDD->Draw("p,same,e1");

//     DrawGammaSetMarker(histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV_wSDD, markerStyleMCAddSig, markerSizePP2760GeV, colorMCPythiaAddSigPP2760GeV , colorMCPythiaAddSigPP2760GeV);
//     histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV_wSDD->Draw("p,same,e1");

    DrawGammaSetMarker(histoEtaEfficiency_LHC12i3_PYT8_2760GeV_wSDD, markerStyleMCC, markerSizePP2760GeV, colorMCPythia2PP2760GeV , colorMCPythia2PP2760GeV);
    histoEtaEfficiency_LHC12i3_PYT8_2760GeV_wSDD->Draw("p,same,e1");
    
    DrawGammaSetMarker(histoEtaEfficiency_LHC12f1b_PHO_2760GeV_wSDD, markerStyleMCB, markerSizePP2760GeV, colorMCPhojetPP2760GeV , colorMCPhojetPP2760GeV);
    histoEtaEfficiency_LHC12f1b_PHO_2760GeV_wSDD->Draw("p,same,e1");

    TLegend* legendEfficiencyEta2760GeV                                         = new TLegend(0.56,0.09,0.93,0.3);
    legendEfficiencyEta2760GeV->SetFillColor(0);
    legendEfficiencyEta2760GeV->SetLineColor(0);
    legendEfficiencyEta2760GeV->SetTextSize(0.025);
    legendEfficiencyEta2760GeV->SetMargin(0.2);
    legendEfficiencyEta2760GeV->AddEntry(histoEtaEfficiency_LHC12f1b_PHO_2760GeV_wSDD,"Phojet, LHC12f1b","p");
    legendEfficiencyEta2760GeV->AddEntry(histoEtaEfficiency_LHC12f1a_PYT8_2760GeV_wSDD,"Pythia 8, LHC12f1a","p");
//     legendEfficiencyEta2760GeV->AddEntry(histoEtaEfficiency_LHC12i3_PYT8_2760GeV_wSDD,"Pythia 8, LHC12i3","p");
//     legendEfficiencyEta2760GeV->AddEntry(histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV_wSDD,"Pythia 8, LHC12i3, added Signals","p");
    legendEfficiencyEta2760GeV->Draw();
    
    canvasEtaEfficiencies2760GeVDiffMC->Update();
    canvasEtaEfficiencies2760GeVDiffMC->Print(Form("%s/Eta_Efficiency_WSDD_2760GeV.%s",outputDir.Data(),suffix.Data()));
    
    canvasEtaEfficiencies2760GeVDiffMC->SetLogy(0);
    canvasEtaEfficiencies2760GeVDiffMC->SetTopMargin(0.035);
    
    histo2DEtaEffi2760GeV->GetYaxis()->SetRangeUser(1e-5,4E-2);
    histo2DEtaEffi2760GeV->DrawCopy(); 
    histoEtaEfficiency_LHC12f1a_PYT8_2760GeV_wSDD->Draw("p,same,e1");
//     histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV_wSDD->Draw("p,same,e1");
    histoEtaEfficiency_LHC12i3_PYT8_2760GeV_wSDD->Draw("p,same,e1");
    histoEtaEfficiency_LHC12f1b_PHO_2760GeV_wSDD->Draw("p,same,e1");
    legendEfficiencyEta2760GeV->SetX1NDC(0.1);
    legendEfficiencyEta2760GeV->SetX2NDC(0.1+0.37);
    legendEfficiencyEta2760GeV->SetY1NDC(0.94-0.21);
    legendEfficiencyEta2760GeV->SetY2NDC(0.94);
    legendEfficiencyEta2760GeV->Draw();
    
    canvasEtaEfficiencies2760GeVDiffMC->Update();
    canvasEtaEfficiencies2760GeVDiffMC->Print(Form("%s/Eta_Efficiency_WSDD_2760GeV_LinY.%s",outputDir.Data(),suffix.Data()));

    //    **********************************************************************************************************************
    //    ***********************Compare Efficiencies for eta at 2.76TeV for different MC with weighting ***********************
    //    **********************************************************************************************************************

    canvasEtaEfficiencies2760GeVDiffMC->cd();
    canvasEtaEfficiencies2760GeVDiffMC->SetLogy(1);
    histo2DEtaEffi2760GeV->DrawCopy(); 

    DrawGammaSetMarker(histoEtaEfficiency_LHC12f1a_PYT8_2760GeV_wSDDRW, markerStyleMCA, markerSizePP2760GeV, colorMCPythiaPP2760GeV , colorMCPythiaPP2760GeV);
    histoEtaEfficiency_LHC12f1a_PYT8_2760GeV_wSDDRW->Draw("p,same,e1");

//     DrawGammaSetMarker(histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV_wSDDRW, markerStyleMCAddSig, markerSizePP2760GeV, colorMCPythiaAddSigPP2760GeV , colorMCPythiaAddSigPP2760GeV);
//     histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV_wSDDRW->Draw("p,same,e1");

//     DrawGammaSetMarker(histoEtaEfficiency_LHC12i3_PYT8_2760GeV_wSDDRW, markerStyleMCC, markerSizePP2760GeV, colorMCPythia2PP2760GeV , colorMCPythia2PP2760GeV);
//     histoEtaEfficiency_LHC12i3_PYT8_2760GeV_wSDDRW->Draw("p,same,e1");
    
    DrawGammaSetMarker(histoEtaEfficiency_LHC12f1b_PHO_2760GeV_wSDDRW, markerStyleMCB, markerSizePP2760GeV, colorMCPhojetPP2760GeV , colorMCPhojetPP2760GeV);
    histoEtaEfficiency_LHC12f1b_PHO_2760GeV_wSDDRW->Draw("p,same,e1");

    DrawGammaSetMarker(histoEtaEfficiency2760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, color2760GeV , color2760GeV);
    histoEtaEfficiency2760GeV->Draw("p,same,e1");
    
    legendEfficiencyEta2760GeV->AddEntry(histoEtaEfficiency2760GeV,"final merged efficiency","p");
    legendEfficiencyEta2760GeV->SetX1NDC(0.56);
    legendEfficiencyEta2760GeV->SetX2NDC(0.56+0.37);
    legendEfficiencyEta2760GeV->SetY1NDC(0.09-0.21);
    legendEfficiencyEta2760GeV->SetY2NDC(0.09);
    legendEfficiencyEta2760GeV->Draw();
    
    canvasEtaEfficiencies2760GeVDiffMC->Update();
    canvasEtaEfficiencies2760GeVDiffMC->Print(Form("%s/Eta_Efficiency_WSDD_2760GeV_reweighted.%s",outputDir.Data(),suffix.Data()));
    
    canvasEtaEfficiencies2760GeVDiffMC->SetLogy(0);
    canvasEtaEfficiencies2760GeVDiffMC->SetTopMargin(0.035);
    
    histo2DEtaEffi2760GeV->GetYaxis()->SetRangeUser(1e-5,2E-2);
    histo2DEtaEffi2760GeV->DrawCopy(); 
    
    histoEtaEfficiency_LHC12f1a_PYT8_2760GeV_wSDDRW->Draw("p,same,e1");
//     histoEtaEfficiency_LHC12i3_PYT8addSig_2760GeV_wSDDRW->Draw("p,same,e1");
//     histoEtaEfficiency_LHC12i3_PYT8_2760GeV_wSDDRW->Draw("p,same,e1");
    histoEtaEfficiency_LHC12f1b_PHO_2760GeV_wSDDRW->Draw("p,same,e1");
    histoEtaEfficiency2760GeV->Draw("p,same,e1");
    
    legendEfficiencyEta2760GeV->SetX1NDC(0.1);
    legendEfficiencyEta2760GeV->SetX2NDC(0.1+0.37);
    legendEfficiencyEta2760GeV->SetY1NDC(0.94-0.21);
    legendEfficiencyEta2760GeV->SetY2NDC(0.94);
    legendEfficiencyEta2760GeV->Draw();
    
    canvasEtaEfficiencies2760GeVDiffMC->Update();
    canvasEtaEfficiencies2760GeVDiffMC->Print(Form("%s/Eta_Efficiency_WSDD_2760GeV_reweighted_LinY.%s",outputDir.Data(),suffix.Data()));
    
    //    **********************************************************************************************************************
    //    ******************************** Efficiencies for eta and pi0 at 2.76TeV after weighting *****************************
    //    **********************************************************************************************************************
    
    TCanvas* canvasEfficiencies2760GeVDiffMC                         = new TCanvas("canvasEfficiencies2760GeVDiffMC", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasEfficiencies2760GeVDiffMC,  0.09, 0.01, 0.035, 0.08);
    canvasEfficiencies2760GeVDiffMC->SetLogy(0);
    
    TH2F * histo2DEffi2760GeV;
    histo2DEffi2760GeV                                                     = new TH2F("histo2DEffi2760GeV", "histo2DEffi2760GeV",1000, 0., 8.5, 1000, 1e-6, 1e-1 );
    SetStyleHistoTH2ForGraphs( histo2DEffi2760GeV, "#it{p}_{T} (GeV/#it{c})", "#epsilon", 
                            0.03, 0.04, 0.03, 0.04, 0.83, 1.05);
    histo2DEffi2760GeV->GetYaxis()->SetLabelOffset(0.01);
    histo2DEffi2760GeV->GetYaxis()->SetRangeUser(1e-5,1.6E-2);
    histo2DEffi2760GeV->DrawCopy(); 

    DrawGammaSetMarker(histoPi0Efficiency2760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, color2760GeV , color2760GeV);
    histoPi0Efficiency2760GeV->Draw("p,same,e1");
    
    DrawGammaSetMarker(histoEtaEfficiency2760GeV, markerStyleMCC, markerSizePP2760GeV, color2760GeV , color2760GeV);
    histoEtaEfficiency2760GeV->Draw("p,same,e1");

    TLatex *labelEnergyEff = new TLatex(0.64,0.30,collisionSystemPP2760GeV.Data());
    SetStyleTLatex( labelEnergyEff, 0.03,4);
    labelEnergyEff->Draw();
    TLatex *labelDetSysEff = new TLatex(0.64,0.26,detectionSystemEMCALPCM.Data());
    SetStyleTLatex( labelDetSysEff, 0.03,4);
    labelDetSysEff->Draw();

    TLegend* legendEfficiency2760GeV                                         = new TLegend(0.72,0.15,0.9,0.25);
    legendEfficiency2760GeV->SetFillColor(0);
    legendEfficiency2760GeV->SetLineColor(0);
    legendEfficiency2760GeV->SetTextSize(0.03);
    legendEfficiency2760GeV->SetMargin(0.3);
    legendEfficiency2760GeV->AddEntry(histoPi0Efficiency2760GeV,"#pi^{0}","p");
    legendEfficiency2760GeV->AddEntry(histoEtaEfficiency2760GeV,"#eta","p");
    legendEfficiency2760GeV->Draw();

    canvasEfficiencies2760GeVDiffMC->Update();
    canvasEfficiencies2760GeVDiffMC->Print(Form("%s/Efficiency_WSDD_2760GeV_reweighted_LinY.%s",outputDir.Data(),suffix.Data()));

    //    **********************************************************************************************************************
    //    ******************************** Acceptances for eta and pi0 at 2.76TeV after weighting ******************************
    //    **********************************************************************************************************************
    
    TCanvas* canvasAcceptance2760GeVDiffMC                         = new TCanvas("canvasAcceptance2760GeVDiffMC", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasAcceptance2760GeVDiffMC,  0.09, 0.01, 0.015, 0.08);
    canvasAcceptance2760GeVDiffMC->SetLogy(0);
    
    TH2F * histo2DAcc2760GeV;
    histo2DAcc2760GeV                                             = new TH2F("histo2DAcc2760GeV", "histo2DAcc2760GeV",1000, 0., 8.5, 1000, 0.6, 1.05 );
    SetStyleHistoTH2ForGraphs( histo2DAcc2760GeV, "#it{p}_{T} (GeV/#it{c})", "A", 
                            0.03, 0.04, 0.03, 0.04, 0.83, 1.2);
    histo2DAcc2760GeV->GetYaxis()->SetLabelOffset(0.01);
//     histo2DAcc2760GeV->GetYaxis()->SetRangeUser(1e-5,1.6E-2);
    histo2DAcc2760GeV->DrawCopy(); 

    DrawGammaSetMarker(histoPi0Acceptance2760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, color2760GeV , color2760GeV);
    histoPi0Acceptance2760GeV->Draw("p,same,e1");
    
    DrawGammaSetMarker(histoEtaAcceptance2760GeV, markerStyleMCC, markerSizePP2760GeV, color2760GeV , color2760GeV);
    histoEtaAcceptance2760GeV->Draw("p,same,e1");

    labelEnergyEff->Draw();
    labelDetSysEff->Draw();
    
    TLegend* legendAcceptance2760GeV                                         = new TLegend(0.72,0.15,0.9,0.25);
    legendAcceptance2760GeV->SetFillColor(0);
    legendAcceptance2760GeV->SetLineColor(0);
    legendAcceptance2760GeV->SetTextSize(0.03);
    legendAcceptance2760GeV->SetMargin(0.3);
    legendAcceptance2760GeV->AddEntry(histoPi0Acceptance2760GeV,"#pi^{0}","p");
    legendAcceptance2760GeV->AddEntry(histoEtaAcceptance2760GeV,"#eta","p");
    legendAcceptance2760GeV->Draw();

    canvasAcceptance2760GeVDiffMC->Update();
    canvasAcceptance2760GeVDiffMC->Print(Form("%s/Acceptance_WSDD_2760GeV_reweighted_LinY.%s",outputDir.Data(),suffix.Data()));
    
    TCanvas* canvasRawYield = new TCanvas("canvasRawYield","",200,10,1350*1.4,1350);  // gives the page size
    DrawGammaCanvasSettings( canvasRawYield, 0.12, 0.02, 0.035, 0.09);
    TH2F * histo2DRaw;
//       canvasRawYield->SetLogx();
    canvasRawYield->SetLogy();
    histo2DRaw = new TH2F("histo2DRaw","histo2DRawEta",1000,0.,9,2000,1.e-7,3e-3  );
    SetStyleHistoTH2ForGraphs(histo2DRaw, "#it{p}_{T} (GeV/#it{c})","#frac{d#it{N}_{raw}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
    histo2DRaw->Draw("copy");


        DrawGammaSetMarker(histoPCMEMCALRAWYieldPi02760GeV, markerStyleSpectrum2760GeV, markerSizePP2760GeV, color2760GeV , color2760GeV);
        histoPCMEMCALRAWYieldPi02760GeV->Draw("p,same,e1");
    
        DrawGammaSetMarker(histoPCMEMCALRAWYieldEta2760GeV, markerStyleMCC, markerSizePP2760GeV, color2760GeV , color2760GeV);
        histoPCMEMCALRAWYieldEta2760GeV->Draw("p,same,e1");

        TLatex *labelEnergyRawYield = new TLatex(0.64,0.88,collisionSystemPP2760GeV.Data());
        SetStyleTLatex( labelEnergyRawYield, 0.04,4);
        labelEnergyRawYield->Draw();
        TLatex *labelDetSysRawYield = new TLatex(0.64,0.84,detectionSystemEMCALPCM.Data());
        SetStyleTLatex( labelDetSysRawYield, 0.04,4);
        labelDetSysRawYield->Draw();

        
        TLegend* legendRawYield = new TLegend(0.72,0.72,0.9,0.82);
        legendRawYield->SetFillColor(0);
        legendRawYield->SetLineColor(0);
        legendRawYield->SetTextSize(0.04);
        legendRawYield->AddEntry(histoPCMEMCALRAWYieldPi02760GeV,"#pi^{0}","p");
        legendRawYield->AddEntry(histoPCMEMCALRAWYieldEta2760GeV,"#eta","p");
        legendRawYield->Draw();

canvasRawYield->SaveAs(Form("%s/RawYieldCompEta.%s",outputDir.Data(),suffix.Data()));

    TCanvas* canvasAcceptanceTimesEff                         = new TCanvas("canvasAcceptanceTimesEff", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasAcceptanceTimesEff,  0.09, 0.01, 0.015, 0.08);
    canvasAcceptanceTimesEff->SetLogy(1);
    
    TH2F * histo2DAccEff2760GeV;
    histo2DAccEff2760GeV                                             = new TH2F("histo2DAccEff2760GeV", "histo2DAccEff2760GeV",1000, 0.,  8.5, 1000, 8e-5, 3e-2 );
    SetStyleHistoTH2ForGraphs( histo2DAccEff2760GeV, "#it{p}_{T} (GeV/#it{c})", "A #times #epsilon_{eff}", 
                            0.03, 0.04, 0.03, 0.04, 0.83, 1.1);
    histo2DAccEff2760GeV->GetYaxis()->SetLabelOffset(0.001);
    histo2DAccEff2760GeV->DrawCopy(); 

    
    DrawGammaSetMarker(histoPi0AcceptanceTimesEff2760GeV, markerStylePCMEMCAL, markerSizePCMEMCAL, colorPCMEMCAL , colorPCMEMCAL);
    histoPi0AcceptanceTimesEff2760GeV->Draw("p,same,e1");
    
    labelEnergyEff->Draw();
    labelDetSysEff->Draw();

    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/AcceptanceTimesEff_WSDD_2760GeV_reweighted.%s",outputDir.Data(),suffix.Data()));
    
    //    **********************************************************************************************************************
    //    ****************************************Write fits & input MC spectra to file ****************************************
    //    **********************************************************************************************************************
    TFile fMCSpectraInput("MCSpectraInputPCMEMCALpp.root","UPDATE");
//      if (fitInvYieldDataPi0Comb2760GeV){
//          fitInvYieldDataPi0Comb2760GeV->SetRange(0,30);
//          fitInvYieldDataPi0Comb2760GeV->Write("Pi0_Fit_Data_2760GeV",TObject::kOverwrite);
//      }
//      if (fitInvYieldDataEtaComb2760GeV){
//          fitInvYieldDataEtaComb2760GeV->SetRange(0,30);
//          fitInvYieldDataEtaComb2760GeV->Write("Eta_Fit_Data_2760GeV",TObject::kOverwrite);
//      }
        if (histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD){
            histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD->SetTitle("Pi0_Phojet_LHC12f1b_WSDD_2760GeV");
            histoPi0InputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD->Write("Pi0_Phojet_LHC12f1b_WSDD_2760GeV",TObject::kOverwrite);
        }
        if (histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD){
            histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD->SetTitle("Pi0_Pythia8_LHC12f1a_WSDD_2760GeV");
            histoPi0InputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD->Write("Pi0_Pythia8_LHC12f1a_WSDD_2760GeV",TObject::kOverwrite);
        }
        if (histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD){
            histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD->SetTitle("Pi0_Pythia8_LHC12i3_WSDD_2760GeV");
            histoPi0InputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD->Write("Pi0_Pythia8_LHC12i3_WSDD_2760GeV",TObject::kOverwrite);
        }
        if (histoPi0InputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD){
            histoPi0InputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD->SetTitle("Pi0_Pythia8_LHC12i3_WSDD_addSig_2760GeV");
            histoPi0InputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD->Write("Pi0_Pythia8_LHC12i3_WSDD_addSig_2760GeV",TObject::kOverwrite);
        }
        if (histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD){
            histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD->SetTitle("Eta_Phojet_LHC12f1b_WSDD_2760GeV");
            histoEtaInputMCWOWeights_LHC12f1b_PHO_2760GeV_wSDD->Write("Eta_Phojet_LHC12f1b_WSDD_2760GeV",TObject::kOverwrite);
        }
        if (histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD){
            histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD->SetTitle("Eta_Pythia8_LHC12f1a_WSDD_2760GeV");
            histoEtaInputMCWOWeights_LHC12f1a_PYT8_2760GeV_wSDD->Write("Eta_Pythia8_LHC12f1a_WSDD_2760GeV",TObject::kOverwrite);
        }
        if (histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD){
            histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD->SetTitle("Eta_Pythia8_LHC12i3_WSDD_2760GeV");
            histoEtaInputMCWOWeights_LHC12i3_PYT8_2760GeV_wSDD->Write("Eta_Pythia8_LHC12i3_WSDD_2760GeV",TObject::kOverwrite);
        }
        if (histoEtaInputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD){
            histoEtaInputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD->SetTitle("Eta_Pythia8_LHC12i3_WSDD_addSig_2760GeV");
            histoEtaInputMCWOWeights_LHC12i3_PYT8addSig_2760GeV_wSDD->Write("Eta_Pythia8_LHC12i3_WSDD_addSig_2760GeV",TObject::kOverwrite);
        }
    fMCSpectraInput.Close();
        
}

