// provided by Gamma Conversion Group, $ALICE_ROOT/PWGGA/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion
// ***************************************************************************************************************
// **************** this macro has been written by                                                          ******
// ***************** Friederike Bock, fbock@cern.ch  on the basis of work by Martin Wilde                   ******
// ***************** as of 29th Oct 2015 it superceeds CorrectGamma.C ********************************************
// ***************************************************************************************************************

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
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#ifdef __CLING__
// include something for ROOT6 here
R__ADD_INCLUDE_PATH($ALICE_ROOT)
R__ADD_INCLUDE_PATH($ROOUNFOLD_ROOT)
#else
#include "../RooUnfold/src/RooUnfold.h"
#include "../RooUnfold/src/RooUnfoldResponse.h"
#include "../RooUnfold/src/RooUnfoldBayes.h"
#include "../RooUnfold/src/RooUnfoldBinByBin.h"
#endif
#include "CorrectGammaV2.h"

//*********************************************************************************************************
//***************************** Main routine to correct inclusive photons *********************************
//*********************************************************************************************************
void  CorrectGammaV2(   const char *nameUnCorrectedFile     = "myOutput",
                        const char *nameCorrectionFile      = "",
                        TString cutSelection                = "",
                        TString suffix                      = "eps",
                        TString nameMeson                   = "Pi0",
                        TString isMC                        = "kFALSE",
                        TString energy                      = "",
                        TString optionPeriod                = "",
                        TString fEstimatePileup             = "",
                        Int_t mode                          = 9,
                        TString directPhoton                = ""
                    ){

    gROOT->Reset();

    Int_t   nIterationsUnfolding                = 4;
    Bool_t  useUnfoldingForCocktailSecondaries  = kFALSE;

    //****************************************************************************************
    //*************************** Catch old versions *****************************************
    //****************************************************************************************
    if (mode == 1 || mode > 5){
        cout << "This macro isn't designed to run for Dalitz or the old software version" << endl;
        return;
    } else if ( mode == 4 || mode == 5){
        cout << "WARNING: this macro is still under construction for this mode" << endl;
    } else if ( mode == 2 || mode == 3){
        cout << "WARNING: running hybrid mode, this macro is still under construction for this mode" << endl;
    }

    //*****************************************************************************************
    //************************** Set general style settings ***********************************
    //*****************************************************************************************
    StyleSettingsThesis();
    SetPlotStyle();

    TString nameRec;
    Bool_t isRunMC                  = kFALSE;
    if (isMC.CompareTo("kTRUE")==0){
        nameRec                     = "MC";
        isRunMC                     = kTRUE;
    } else {
        nameRec                     = "data";
    }

    //*****************************************************************************************
    //*************************** Determine K0s scaling factor ********************************
    //*****************************************************************************************
    Double_t doubleAddFactorK0s     = ReturnCorrectK0ScalingFactor( energy,  cutSelection);
    if(isRunMC){
        doubleAddFactorK0s = 0.;
        cout<<" MC --> doubleAddFactorK0s = 0"<<endl;
    } else {
        cout << "The additional K0 correction factor is: "  << doubleAddFactorK0s<<endl;
    }
    Double_t scaleFactorsSec[4]     = {1.+doubleAddFactorK0s, 1., 1., 1.};

    //*****************************************************************************************
    //*************************** Find out cutSelections for different things *****************
    //*****************************************************************************************
    TString fEventCutSelection      = "";
    TString fGammaCutSelection      = "";
    TString fClusterCutSelection    = "";
    TString fElectronCutSelection   = "";
    TString fMesonCutSelection      = "";
    TString fCutSelection           = cutSelection;
    ReturnSeparatedCutNumberAdvanced(cutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);

    Bool_t isPCM                    = 0;
    Bool_t isCalo                   = 0;
    if ( mode == 0 || mode == 2 || mode == 3)
        isPCM                       = 1;
    if ( mode == 2 || mode == 3 || mode == 4 || mode == 5)
        isCalo                      = 1;

    //*****************************************************************************************
    //*************************** Determine centrality from CutString *************************
    //*****************************************************************************************
    TString centrality              = GetCentralityString(fEventCutSelection);
    TString collisionSystem         = ReturnFullCollisionsSystem(energy);
    TString cent                    = "";
    TString textMeasurement         = "#gamma";
    TString detectionProcess        = ReturnFullTextReconstructionProcess(mode);
    TString detectionProcess1       = "";
    TString detectionProcess2       = "";
    if (isPCM && isCalo){
        detectionProcess1           = ReturnFullTextReconstructionProcess(mode,1);
        detectionProcess            = detectionProcess1;
        detectionProcess2           = ReturnFullTextReconstructionProcess(mode,2);
    }
    if(energy.Contains("PbPb") || (energy.Contains("pPb") && centrality.CompareTo("0-100%") != 0)){
        cent                        = Form("%s %s", centrality.Data(), collisionSystem.Data());
    } else {
        cent                        = collisionSystem;
    }

    //******************************************************************************************
    //***************************** Set output directory ***************************************
    //******************************************************************************************
    TString outputDir               = Form("%s/%s/%s/CorrectGammaV2",cutSelection.Data(),energy.Data(),suffix.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    //******************************************************************************************
    //************************ Define additional global variables ******************************
    //******************************************************************************************
    TString textPi0New              = Form("Gamma_%s",nameMeson.Data());

    Double_t deltaEta               = ReturnDeltaEta(fGammaCutSelection);
    Double_t eta                    = deltaEta*0.5;
    Double_t deltaEtaCalo           = 0;
    Double_t deltaPhiCalo           = 0;
    Double_t etaCalo                = 0;
    if (isCalo){
        deltaEtaCalo                = ReturnDeltaEtaCalo(fClusterCutSelection);
        etaCalo                     = deltaEtaCalo*0.5;
        deltaPhiCalo                = ReturnDeltaPhiCalo(fClusterCutSelection);
    }

    Double_t scaling                = 1./(2.*TMath::Pi());
    Double_t scalingCalo            = 0.;
    if (isCalo) scalingCalo         = 1./deltaPhiCalo;

    Bool_t kDoPileup                = kFALSE;
    if (fEstimatePileup.CompareTo("EstimateTrainPileUp") == 0) {
        kDoPileup                   = kTRUE;
        if (mode == 4 || mode == 5) {
            cout << "requested out-of-bunch pileup correction for calo mode, skipping" << endl;
            kDoPileup               = kFALSE;
        }
    } else if ( mode == 2 || mode == 3 ) {
        kDoPileup                   = kTRUE;
    }

    //******************************************************************************************
    //********************************** secondary labels **************************************
    //******************************************************************************************
    TString nameSecondaries[4]                              = { "K0s", "K0l", "Lambda", "Rest" };
    TString nameLabelSecondaries[4]                         = { "K^{0}_{s}", "K^{0}_{l}", "#Lambda", "rest" };
    Style_t markerStyleSec[4]                               = { 21, 33, 29, 34};
    Style_t markerStyleSecWithToy[4]                        = { 25, 27, 30, 28};
    Size_t  markerSizeSec[4]                                = { 1.5, 1.75, 2., 1.5};
    Color_t colorSecFromToy[4]                              = { kRed-2, kOrange+1, kCyan-2, kBlue+2};

    //******************************************************************************************
    //************************ True Combinatorial Background labels ****************************
    //******************************************************************************************
    TString combinatorics[17]                               = { "Elec+Elec","Elec+Pion","Elec+Kaon","Elec+Proton","Elec+Muon","Pion+Pion","Pion+Kaon","Pion+Proton",
                                                                "Pion+Muon","Kaon+Kaon","Kaon+Proton","Kaon+Muon","Proton+Proton","Proton+Muon","Muon+Muon","Rest","All"};
    TString combinatoricsLabels[17]                         = { "e^{#pm}e^{#mp}","e^{#pm}#pi^{#mp}","e^{#pm}K^{#mp}","e^{#pm}#bar{p}(p)","e^{#pm}#mu^{#mp}","#pi^{#pm}#pi^{#mp}","#pi^{#pm}K^{#mp}",
                                                                "#pi^{#pm}#bar{p}(p)","#pi^{#pm}#mu^{#mp}","K^{#pm}K^{#mp}","K^{#pm}#bar{p}(p)","K^{#pm}#mu^{#mp}","p(#bar{p})#bar{p}(p)","p(#bar{p})#mu^{#mp}","#mu^{#pm}#mu^{#mp}","rest","all"};
    Color_t colorsCombinatorics[17]                         = { kAzure-1, kRed+1, kOrange+7, kMagenta+2, kRed-9,
                                                                kBlue-6, kAzure+5, kPink, kCyan-3, kGreen-3,
                                                                kSpring+9, kGreen+2, kBlue+2, kMagenta-6, kSpring+4,
                                                                kCyan+2, 809};
    Style_t markersCombinatorics[17]                        = { 20, 21, 24, 25, 27,
                                                                28, 29, 30, 33, 34,
                                                                20, 21, 24, 25, 27,
                                                                28, 29};
    TString combinatoricsCalo[11]                           = { "Electron","Pion","Proton","Kaon","Neutron","K0s","Lambda","Muon","K0l","Rest","All"};
    TString combinatoricsLabelsCalo[11]                     = { "e^{#pm}","#pi^{#pm}","p(#bar{p})","K^{#pm}","n","K^{0}_{s}","#Lambda","#mu^{#pm}","K^{0}_{l}","rest","all"};
    Color_t colorsCombinatoricsCalo[11]                     = { kAzure-1, kRed+1, kOrange+7, kMagenta+2, kRed-9, kBlue,
                                                                kBlue-6, kAzure+5, kPink, kCyan-3, kGreen-3};
    Style_t markersCombinatoricsCalo[11]                    = { 20, 21, 24, 25, 27,
                                                                28, 29, 30, 33, 34, 20};

    //******************************************************************************************
    //********************************** decay labels ******************************************
    //******************************************************************************************
    TString decays[9]                                       = { "Pi0","Eta","Etap","Omega","Rho",
                                                                "Phi","Sigma", "All decays","direct #gamma"};
    TString decaysLabels[9]                                 = { "#gamma from #pi^{0}","#gamma from #eta","#gamma from #eta'","#gamma from #omega","#gamma from #rho^{0}",
                                                                "#gamma from #phi", "#gamma from #Sigma^{0}", "All decays","direct #gamma"};
    Color_t colorsDecay[9]                                  = { kRed+1, kAzure-1, kOrange+7, kGreen+2, kCyan-3, 809,
                                                                kBlue+2, kBlack, kRed-1};

    //******************************************************************************************
    //*******************File definitions and reading histograms from files ********************
    //******************************************************************************************
    TFile*  fileUnCorrected                                         = new TFile(nameUnCorrectedFile);
    TFile*  fileCorrections                                         = new TFile(nameCorrectionFile);
    if (!fileUnCorrected) {
        cout << "file " << nameUnCorrectedFile << " not found!" << endl;
        return;
    }
    if (!nameCorrectionFile) {
        cout << "file " << nameCorrectionFile << " not found!" << endl;
        return;
    }
    TFile*  doPileUpCorr                                            = NULL;
    if (kDoPileup){
        doPileUpCorr                                                = new  TFile(Form("%s/%s/%s_%s_GammaConvV1DCAHistogramms%s_%s.root", cutSelection.Data(), energy.Data(), nameMeson.Data(),
                                                                                      nameRec.Data(), optionPeriod.Data(), cutSelection.Data()));
        if (doPileUpCorr->IsZombie())   doPileUpCorr                = 0;
    } else                              doPileUpCorr                = 0;


    //******************* Calculate number of events for normalization *************************
    TH1D*   histoEventQuality                                       = (TH1D*)fileUnCorrected->Get("NEvents");
    Float_t nEvt                                                    = 0.;
    if (energy.Contains("PbPb") || energy.Contains("pPb")  )  nEvt  = histoEventQuality->GetBinContent(1);
    else                                        nEvt                = GetNEvents(histoEventQuality);

    TH1F*   histoEventQualityMC                                     = (TH1F*)fileCorrections->Get("NEvents");
    Float_t nEvtMC                                                  = 0.;
    if (energy.Contains("PbPb") || energy.Contains("pPb") )  nEvtMC = histoEventQualityMC->GetBinContent(1);
    else                                        nEvtMC              = GetNEvents(histoEventQualityMC);

    //******************* Binning histogram ****************************************************
    TH1D*   deltaPt                                                 = (TH1D*)fileUnCorrected->Get("deltaPt");
    if(!deltaPt) deltaPt                                            = (TH1D*)fileUnCorrected->Get("deltaPtGamma");
    for (Int_t i = 0; i < deltaPt->GetNbinsX() +1; i++) deltaPt->SetBinError(i, 0);

    //******************* Conversion gamma spectra *********************************************
    TH1D*   histoESDConvGammaPt                                     = NULL;
    TH1D*   histoESDConvGammaPt_OrBin                               = NULL;
    if (isPCM) {
        histoESDConvGammaPt                                         = (TH1D*)fileUnCorrected->Get("ESD_ConvGamma_Pt");
        histoESDConvGammaPt_OrBin                                   = (TH1D*)fileUnCorrected->Get("ESD_ConvGamma_Pt_OriginalBinning");
    }

    //******************* Calo gamma spectra ***************************************************
    TH1D*   histoESDCaloGammaPt                                     = NULL;
    TH1D*   histoESDCaloGammaPt_OrBin                               = NULL;
    if (isCalo) {
        histoESDCaloGammaPt                                         = (TH1D*)fileUnCorrected->Get("ESD_CaloGamma_Pt");
        histoESDCaloGammaPt_OrBin                                   = (TH1D*)fileUnCorrected->Get("ESD_CaloGamma_Pt_OriginalBinning");
    }

    //******************* Determine max Pt *****************************************************
    Double_t                maxPtGamma                              = 0.;
    Double_t                minPtGamma                              = 0.;
    if (isPCM && !isCalo){
        maxPtGamma                              = histoESDConvGammaPt->GetXaxis()->GetBinUpEdge(histoESDConvGammaPt->GetNbinsX());
        if(energy.Contains("PbPb"))
            minPtGamma                          = histoESDConvGammaPt->GetXaxis()->GetBinLowEdge(1);
    } else if (isPCM && isCalo){
        maxPtGamma                              = histoESDConvGammaPt->GetXaxis()->GetBinUpEdge(histoESDConvGammaPt->GetNbinsX());
    } else if (isCalo && !isPCM){
        maxPtGamma                              = histoESDCaloGammaPt->GetXaxis()->GetBinUpEdge(histoESDCaloGammaPt->GetNbinsX());
    }
    //******************* Primary gamma correction factors (PCM and Calo) **********************
    TH1D*   histoGammaPurity_Pt                                     = NULL;
    TH1D*   histoGammaTruePurity_Pt                                 = NULL;
    TH1D*   histoGammaTruePurity_Pt_OrBin                           = NULL;
    TH1D*   histoGammaPrimaryRecoEff_Pt                             = NULL;
    TH1D*   histoGammaPrimaryRecoEff_Pt_OrBin                       = NULL;
    TH1D*   histoGammaPrimaryRecoEff_MCPt                           = NULL;
    TH1D*   histoGammaPrimaryRecoEff_MCPt_OrBin                     = NULL;
    if ( isPCM || isCalo ) {
        // reconstructed validated photons / all reconstructed photon candidates in MC
        histoGammaPurity_Pt                                         = (TH1D*)fileCorrections->Get("GammaPurity_Pt");
        // reconstructed validated primary photons / (all reconstructed photon candidates in MC - validated reconstructed secondary photons)
        // rescaling of histoGammaPurity_Pt for primary photons
        histoGammaTruePurity_Pt                                     = (TH1D*)fileCorrections->Get("GammaTruePurity_Pt");
        histoGammaTruePurity_Pt_OrBin                               = (TH1D*)fileCorrections->Get("GammaTruePurity_OriginalBinning_Pt");
        // reconstructed validated primary photons vs rec pt / converted MC photons vs MCpt
        // this reconstruction efficiency includes the transverse momentum resolution correction and the acceptance correction
        histoGammaPrimaryRecoEff_Pt                                 = (TH1D*)fileCorrections->Get("GammaPrimaryRecoEff_Pt");
        histoGammaPrimaryRecoEff_Pt_OrBin                           = (TH1D*)fileCorrections->Get("GammaPrimaryRecoEff_Pt_OriginalBinning");
        // reconstructed validated primary photons vs MC pt / converted MC photons vs MCpt
        // this reconstruction efficiency includes the acceptance correction
        histoGammaPrimaryRecoEff_MCPt                               = (TH1D*)fileCorrections->Get("GammaPrimaryRecoEff_MCPt");
        histoGammaPrimaryRecoEff_MCPt_OrBin                         = (TH1D*)fileCorrections->Get("GammaPrimaryRecoEff_MCPt_OriginalBinning");
    }
    TH1D*   histoGammaConvProb_MCPt                                 = NULL;
    TH1D*   histoGammaConvProb_MCPt_OrBin                           = NULL;
    TH2D*   histoGammaTruePrimaryConv_recPt_MCPt                    = NULL;
    if ( isPCM ) {
        // converted photons MC/ all generated photons (always based on primaries only)
        histoGammaConvProb_MCPt                                     = (TH1D*)fileCorrections->Get("GammaConvProb_MCPt");
        histoGammaConvProb_MCPt_OrBin                               = (TH1D*)fileCorrections->Get("GammaConvProb_MCPt_OriginalBinning");
        // response matrix for photons recPt vs MC pt (conversions)
        histoGammaTruePrimaryConv_recPt_MCPt                        = (TH2D*)fileCorrections->Get("TruePrimaryConvGamma_recPt_MCPt");
    }
    TH2D*   histoGammaTruePrimaryCalo_recPt_MCPt                    = NULL;
    if (isCalo) {
        // response matrix for photons recPt vs MC pt (calorimeter)
        histoGammaTruePrimaryCalo_recPt_MCPt                        = (TH2D*)fileCorrections->Get("TruePrimaryCaloGamma_recPt_MCPt");
    }

    //******************* Primary gamma correction factors (PCM-Calo hybrid) *******************
    TH1D*   histoGammaCaloPurity_Pt                                 = NULL;
    TH1D*   histoGammaCaloTruePurity_Pt                             = NULL;
    TH1D*   histoGammaCaloTruePurity_Pt_OrBin                       = NULL;
    TH1D*   histoGammaCaloPrimaryRecoEff_Pt                         = NULL;
    TH1D*   histoGammaCaloPrimaryRecoEff_MCPt                       = NULL;
    // Calorimeter photons
    if ( isPCM && isCalo ) {
        // reconstructed validated photons / all reconstructed photon candidates in MC
        histoGammaCaloPurity_Pt                                     = (TH1D*)fileCorrections->Get("GammaCaloPurity_Pt");
        // reconstructed validated primary photons / (all reconstructed photon candidates in MC - validated reconstructed secondary photons)
        // rescaling of histoGammaPurity_Pt for primary photons
        histoGammaCaloTruePurity_Pt                                 = (TH1D*)fileCorrections->Get("GammaCaloTruePurity_Pt");
        histoGammaCaloTruePurity_Pt_OrBin                           = (TH1D*)fileCorrections->Get("GammaCaloTruePurity_OriginalBinning_Pt");
        // reconstructed validated primary photons vs rec pt / converted MC photons vs MCpt
        // this reconstruction efficiency includes the transverse momentum resolution correction and the acceptance correction
        histoGammaCaloPrimaryRecoEff_Pt                             = (TH1D*)fileCorrections->Get("GammaCaloPrimaryRecoEff_Pt");
        // reconstructed validated primary photons vs MC pt / converted MC photons vs MCpt
        // this reconstruction efficiency includes the acceptance correction
        histoGammaCaloPrimaryRecoEff_MCPt                           = (TH1D*)fileCorrections->Get("GammaCaloPrimaryRecoEff_MCPt");
    }

    //******************* MC histograms (PCM and Calo) *****************************************
    TH1D*   histoMCrecGamma_Pt                                      = NULL;
    TH1D*   histoMCrecGamma_Pt_OrBin                                = NULL;
    TH1D*   histoMCrecPrimaryGamma_Pt                               = NULL;
    TH1D*   histoMCGammaConv_MCPt                                   = NULL;
    TH1D*   histoGammaTrueConv_Pt                                   = NULL;
    TH1D*   histoGammaTruePrimaryConv_Pt                            = NULL;
    if (isPCM) {
        // reconstructed photon candidates in MC
        histoMCrecGamma_Pt                                          = (TH1D*)fileCorrections->Get("MCrec_ConvGamma_Pt");
        histoMCrecGamma_Pt_OrBin                                    = (TH1D*)fileCorrections->Get("MCrec_ConvGamma_OriginalBinning_Pt");
        // reconstructed photon candidates in MC - validated secondary photons
        histoMCrecPrimaryGamma_Pt                                   = (TH1D*)fileCorrections->Get("MCrec_PrimaryConvGamma_Pt");
        // MC converted photons regardless of source
        histoMCGammaConv_MCPt                                       = (TH1D*)fileCorrections->Get("MC_ConvGamma_MCPt");
        // validated reconstructed photons
        histoGammaTrueConv_Pt                                       = (TH1D*)fileCorrections->Get("TrueConvGamma_Pt");
        // validated reconstructed primary photons
        histoGammaTruePrimaryConv_Pt                                = (TH1D*)fileCorrections->Get("TruePrimaryConvGamma_Pt");
    }
    TH1D*   histoMCrecBackground_Pt                                 = NULL;
    TH1D*   histoMCAllGamma_MCPt                                    = NULL;
    TH1D*   histoMCAllGamma_OriginalBin_MCPt                        = NULL;
    if (isPCM ||isCalo) {
        // reconstructed background in MC
        histoMCrecBackground_Pt                                     = (TH1D*)fileCorrections->Get("MCrec_Background");
        // MC input for all gamma, regardless of source
        histoMCAllGamma_MCPt                                        = (TH1D*)fileCorrections->Get("MC_AllGamma_MCPt");
        histoMCAllGamma_OriginalBin_MCPt                            = (TH1D*)fileCorrections->Get("MC_AllGamma_OriginalBinning_MCPt");
    }

    //******************* MC histograms (Calo) *************************************************
    TH1D*   histoMCrecGammaCalo_Pt                                  = NULL;
    TH1D*   histoMCrecGammaCalo_Pt_OrBin                            = NULL;
    TH1D*   histoMCrecPrimaryGammaCalo_Pt                           = NULL;
    TH1D*   histoGammaTrueCalo_Pt                                   = NULL;
    TH1D*   histoGammaTruePrimaryCalo_Pt                            = NULL;
    if (isCalo) {
        // reconstructed photon candidates in MC
        histoMCrecGammaCalo_Pt                                      = (TH1D*)fileCorrections->Get("MCrec_CaloGamma_Pt");
        histoMCrecGammaCalo_Pt_OrBin                                = (TH1D*)fileCorrections->Get("MCrec_CaloGamma_OriginalBinning_Pt");
        // reconstructed photon candidates in MC - validated secondary photons
        histoMCrecPrimaryGammaCalo_Pt                               = (TH1D*)fileCorrections->Get("MCrec_PrimaryCaloGamma_Pt");
        // validated reconstructed photons
        histoGammaTrueCalo_Pt                                       = (TH1D*)fileCorrections->Get("TrueCaloGamma_Pt");
        // validated reconstructed primary photons
        histoGammaTruePrimaryCalo_Pt                                = (TH1D*)fileCorrections->Get("TruePrimaryCaloGamma_Pt");
    }
    TH1D*   histoMCrecCaloBackground_Pt                             = NULL;
    TH1D*   histoMCAllGammaCalo_MCPt                                = NULL;
    TH1D*   histoMCAllGammaCalo_OriginalBin_MCPt                    = NULL;
    if (isPCM && isCalo) {
        // reconstructed background in MC
        histoMCrecCaloBackground_Pt                                 = (TH1D*)fileCorrections->Get("MCrec_Calo_Background");
        // MC input for all gamma, regardless of source
        histoMCAllGammaCalo_MCPt                                    = (TH1D*)fileCorrections->Get("MC_AllGammaEMCAcc_MCPt");
        histoMCAllGammaCalo_OriginalBin_MCPt                        = (TH1D*)fileCorrections->Get("MC_AllGammaEMCAcc_OriginalBinning_MCPt");
    }

    //******************* MC decay gammas ******************************************************
    TH1D**  histoPhotonSource_MCPt                                  = new TH1D*[9];
    for (Int_t i = 0;i<7;i++) {
        histoPhotonSource_MCPt[i]                                   = (TH1D*)fileCorrections->Get(Form("MC_DecayGamma%s_Pt",decays[i].Data()));
        histoPhotonSource_MCPt[i]->Sumw2();
        if (i == 0) histoPhotonSource_MCPt[7]                       = (TH1D*)histoPhotonSource_MCPt[0]->Clone("MC_DecayGammaAll_Pt");
        else        histoPhotonSource_MCPt[7]->Add(histoPhotonSource_MCPt[i]);
    }
    histoPhotonSource_MCPt[8]                                       = (TH1D*)histoMCAllGamma_OriginalBin_MCPt->Clone("MC_DirectPhotons");
    histoPhotonSource_MCPt[8]->Sumw2();
    histoPhotonSource_MCPt[8]->Add(histoPhotonSource_MCPt[7],-1);

    //******************* MC combinatorial gammas (PCM) ****************************************
    TH1D**  histoCombinatorialSpecies_Pt                            = NULL;
    if (isPCM) {
        histoCombinatorialSpecies_Pt                                = new TH1D*[17];
        for(Int_t i = 0;i<17;i++){
            histoCombinatorialSpecies_Pt[i]                         = (TH1D*)fileCorrections->Get(Form("ESD_TrueComb%s_Pt",combinatorics[i].Data()));
            histoCombinatorialSpecies_Pt[i]->SetMinimum(1e-10);
        }
    }

    //******************* MC combinatorial gammas (Calo) ***************************************
    TH1D**  histoCombinatorialSpeciesCalo_Pt                        = NULL;
    if (isCalo && !isPCM) {
        histoCombinatorialSpeciesCalo_Pt                            = new TH1D*[11];
        for(Int_t i = 0;i<11;i++){
            histoCombinatorialSpeciesCalo_Pt[i]                     = (TH1D*)fileCorrections->Get(Form("ESD_TrueComb%s_Pt",combinatoricsCalo[i].Data()));
            histoCombinatorialSpeciesCalo_Pt[i]->SetMinimum(1e-10);
        }
    }

    //******************* MC true secondary gammas (PCM) ***************************************
    TH1D*   histoGammaTrueSecConvGammaFromX_Pt[4]                   = {NULL, NULL, NULL, NULL};
    //TH1D*   histoGammaTrueSecConvGammaFromX_Pt_OrBin[4]             = {NULL, NULL, NULL, NULL};
    if (isPCM) {
        for (Int_t k = 0; k < 4; k++){
            histoGammaTrueSecConvGammaFromX_Pt[k]                   = (TH1D*)fileCorrections->Get(Form("TrueSecondaryConvGammaFromXFrom%s_Pt", nameSecondaries[k].Data()));
            //histoGammaTrueSecConvGammaFromX_Pt_OrBin[k]             = (TH1D*)fileCorrections->Get(Form("TrueSecondaryConvGammaFromXFrom%s_Pt_OriginalBinning", nameSecondaries[k].Data()));
        }
    }

    //******************* MC true secondary gammas (Calo) **************************************
    TH1D*   histoGammaTrueSecCaloGammaFromX_Pt[4]                   = {NULL, NULL, NULL, NULL};
    //TH1D*   histoGammaTrueSecCaloGammaFromX_PtOrgBin[4]             = {NULL, NULL, NULL, NULL};
    if (isCalo && !isPCM) {
        for (Int_t k = 0; k < 4; k++){
            histoGammaTrueSecCaloGammaFromX_Pt[k]                   = (TH1D*)fileCorrections->Get(Form("TrueSecondaryCaloGammaFromXFrom%s_Pt", nameSecondaries[k].Data()));
            //histoGammaTrueSecCaloGammaFromX_PtOrgBin[k]             = (TH1D*)fileCorrections->Get(Form("TrueSecondaryCaloGammaFromXFrom%s_Pt_OriginalBinning", nameSecondaries[k].Data()));
        }
    }

    //******************* MC true secondary gamma fractions (PCM or Calo) **********************
    TH1D* histoFracAllGammaToSecFromX_Pt[4]                         = {NULL, NULL, NULL, NULL};
    TH1D* histoFracAllGammaToSecFromX_Pt_OrBin[4]                   = {NULL, NULL, NULL, NULL};
    if ( isPCM || isCalo ) {
        for (Int_t k = 0; k<4; k++){
            histoFracAllGammaToSecFromX_Pt[k]                       = (TH1D*)fileCorrections->Get(Form("FracAllGammaToSecFromXFrom%s",nameSecondaries[k].Data()));
            histoFracAllGammaToSecFromX_Pt_OrBin[k]                 = (TH1D*)fileCorrections->Get(Form("FracAllGammaToSecFromXFrom%sOriginalBinning",nameSecondaries[k].Data()));
        }
    }
    TH1D* histoFracAllGammaToSecFromX_PileUp_Pt[4]                  = {NULL, NULL, NULL, NULL};
    if (isPCM && doPileUpCorr) {
        for (Int_t k = 0; k<4; k++){
            histoFracAllGammaToSecFromX_PileUp_Pt[k]                = (TH1D*)fileCorrections->Get(Form("FracAllGammaToSecFromXFrom%sPileUp",nameSecondaries[k].Data()));
        }
    }

    //******************************************************************************************
    //******************* Load Secondary gamma spectra from cocktail ***************************
    //******************************************************************************************
    Bool_t  hasCocktailInput                                       = kTRUE;
    TH1D*   histoGammaTrueSecCocktailGammaFromX_Pt[4]              = { NULL, NULL, NULL, NULL};
    TH1D*   histoGammaTrueSecCocktailGammaFromX_PtOrBin[4]         = { NULL, NULL, NULL, NULL};
    if( (isPCM || isCalo) && !isRunMC ){
        for (Int_t k = 0; k < 4; k++){
            if (k < 3){
                histoGammaTrueSecCocktailGammaFromX_Pt[k]          = (TH1D*)fileUnCorrected->Get(Form("CocktailSecondaryGammaFromX%s_Pt",nameSecondaries[k].Data()));
                histoGammaTrueSecCocktailGammaFromX_PtOrBin[k]     = (TH1D*)fileUnCorrected->Get(Form("CocktailSecondaryGammaFromX%s_PtOrBin",nameSecondaries[k].Data()));
            } else if ( isPCM ) {
                histoGammaTrueSecCocktailGammaFromX_Pt[k]          = (TH1D*)histoESDConvGammaPt->Clone("TrueSecondaryConvGammaCocktailFromXFromRest_Pt");
                histoGammaTrueSecCocktailGammaFromX_PtOrBin[k]     = (TH1D*)histoESDConvGammaPt_OrBin->Clone("TrueSecondaryConvGammaCocktailFromXFromRest_Pt_OriginalBinning");

                if(histoGammaTrueSecCocktailGammaFromX_Pt[k])
                   histoGammaTrueSecCocktailGammaFromX_Pt[k]->Multiply(histoFracAllGammaToSecFromX_Pt[k]);
               if(histoGammaTrueSecCocktailGammaFromX_PtOrBin[k])
                   histoGammaTrueSecCocktailGammaFromX_PtOrBin[k]->Multiply(histoFracAllGammaToSecFromX_Pt_OrBin[k]);
            } else if ( isCalo && !isPCM ) {
                histoGammaTrueSecCocktailGammaFromX_Pt[k]          = (TH1D*)histoESDCaloGammaPt->Clone("TrueSecondaryCaloGammaCocktailFromXFromRest_Pt");
                histoGammaTrueSecCocktailGammaFromX_PtOrBin[k]     = (TH1D*)histoESDCaloGammaPt_OrBin->Clone("TrueSecondaryCaloGammaCocktailFromXFromRest_Pt_OriginalBinning");

                if(histoGammaTrueSecCocktailGammaFromX_Pt[k])
                   histoGammaTrueSecCocktailGammaFromX_Pt[k]->Multiply(histoFracAllGammaToSecFromX_Pt[k]);
                if(histoGammaTrueSecCocktailGammaFromX_PtOrBin[k])
                   histoGammaTrueSecCocktailGammaFromX_PtOrBin[k]->Multiply(histoFracAllGammaToSecFromX_Pt_OrBin[k]);
            }
            if (!histoGammaTrueSecCocktailGammaFromX_Pt[k] || !histoGammaTrueSecCocktailGammaFromX_PtOrBin[k])
                hasCocktailInput                                   = kFALSE;
        }
        if (!hasCocktailInput) cout << "secondary spectra from cocktail not found, will not use" << endl;
    } else {
        hasCocktailInput                                           = kFALSE;
    }

    if (hasCocktailInput)
        cout << "Will use cocktail-based secondary correction!" << endl;
    else
        cout << "Will use MC-based secondary correction!" << endl;

    //****************************************************************************************** // would be nice if this would be put inside a fct. in the header at some point
    //******************* Secondary gamma reco eff *********************************************
    //******************************************************************************************
    TH1D* histoGammaSecFromXRecoEff_MCPt[3]                         = { NULL, NULL, NULL };
    TH1D* histoGammaSecFromXRecoEff_MCPt_Unscaled[3]                = { NULL, NULL, NULL };
    TH1D* histoGammaSecFromXRecoEff_MCPt_OrBin[3]                   = { NULL, NULL, NULL };
    TH1D* histoGammaSecFromXRecoEff_RecPt[3]                        = { NULL, NULL, NULL };
    TH1D* histoGammaSecFromXRecoEff_RecPt_Unscaled[3]               = { NULL, NULL, NULL };
    TH1D* histoGammaSecFromXRecoEff_RecPt_OrBin[3]                  = { NULL, NULL, NULL };
    TH1D* histoGammaSecFromXRecoEff_RecPt_OrBin_Unscaled[3]         = { NULL, NULL, NULL };
    TH1D* ratioGammaSecEffMCPt[3]                                   = { NULL, NULL, NULL };
    TH1D* ratioGammaSecEffRecPt[3]                                  = { NULL, NULL, NULL };
    Double_t constOffsetEffMCPt[3]                                  = { 1, 1, 1};
    Double_t constOffsetEffRecPt[3]                                 = { 1, 1, 1};

    // fit settings, also used in conv. prob. fit
    enum        funcFitSec {kFitSecConst, kFitSecLinear, kFitSecExp, kFitSecPower};
    Double_t    minPtFitSec[3]                                      = {  1.5,   1.0,  0.0}; // K0s, K0l, Lambda
    Double_t    maxPtFitSec[3]                                      = { 50.0,  50.0, 50.0};
    funcFitSec  FitSecFunction[3]                                   = {kFitSecConst, kFitSecConst, kFitSecConst};

    // adjust settings for different energies and modes
    if (energy.CompareTo("pPb_5.023TeVRun2") == 0 && mode == 0){
        maxPtFitSec[0]                                              = maxPtGamma;
        maxPtFitSec[1]                                              = maxPtGamma;
        maxPtFitSec[2]                                              = 2.5;
    } else if (energy.Contains("pPb")){
        maxPtFitSec[0]                                              = maxPtGamma;
        maxPtFitSec[1]                                              = maxPtGamma;
        maxPtFitSec[2]                                              = maxPtGamma;
    } else if (energy.Contains("2.76TeV") && mode == 4) {
        minPtFitSec[1]                                              = 1.8;
    } else if (!energy.CompareTo("900GeV")) {
        minPtFitSec[1]                                              = 0.4;

        maxPtFitSec[0]                                              = 1.1;
        maxPtFitSec[1]                                              = 1.1;
        maxPtFitSec[2]                                              = 1.1;
    } else if (!energy.CompareTo("7TeV") && mode == 0) {
        minPtFitSec[2]                                              = 0.4;
        maxPtFitSec[2]                                              = 1.6;
    } else if (!energy.CompareTo("7TeV") && mode == 2) {
        FitSecFunction[0]                                           = kFitSecConst;
        minPtFitSec[0]                                              = 0.8;
        maxPtFitSec[0]                                              = 16.0;
        FitSecFunction[1]                                           = kFitSecConst;
        minPtFitSec[1]                                              = 0.8;
        maxPtFitSec[1]                                              = 5.0;
        minPtFitSec[2]                                              = 0.8;
        maxPtFitSec[2]                                              = 2.0;
    } else if (!energy.CompareTo("7TeV") && mode == 4) {
        minPtFitSec[0]                                              = 1.2;
        maxPtFitSec[0]                                              = 16.0;
        FitSecFunction[1]                                           = kFitSecExp;
        minPtFitSec[1]                                              = 1.2;
        maxPtFitSec[1]                                              = 16.0;
        minPtFitSec[2]                                              = 1.2;
        maxPtFitSec[2]                                              = 4.0;
    } else if (energy.BeginsWith("8TeV") && mode == 0) {
        minPtFitSec[1]                                              = 1.0;
        maxPtFitSec[1]                                              = 2.0;
    } else if (energy.BeginsWith("8TeV") && mode == 2) {
        FitSecFunction[0]                                           = kFitSecPower;
        minPtFitSec[0]                                              = 0.8;
        maxPtFitSec[0]                                              = 12.0;
        FitSecFunction[1]                                           = kFitSecLinear;
        minPtFitSec[1]                                              = 0.8;
        maxPtFitSec[1]                                              = 10.0;
        minPtFitSec[2]                                              = 0.8;
        maxPtFitSec[2]                                              = 2.0;
    } else if (energy.BeginsWith("8TeV") && mode == 4) {
        minPtFitSec[0]                                              = 1.2;
        maxPtFitSec[0]                                              = 16.0;
        FitSecFunction[1]                                           = kFitSecLinear;
        minPtFitSec[1]                                              = 1.2;
        maxPtFitSec[1]                                              = 4.5;
        FitSecFunction[2]                                           = kFitSecLinear;
        minPtFitSec[2]                                              = 1.0;
        maxPtFitSec[2]                                              = 6.0;
    } else if (energy.Contains("PbPb_2.76TeV")) {
        FitSecFunction[0]                                           = kFitSecExp;
        maxPtFitSec[0]                                              = 5.;
        maxPtFitSec[1]                                              = 6.;
        maxPtFitSec[2]                                              = 2.4;
    }

    if ( hasCocktailInput && (isPCM || isCalo) ) {
        TF1* constantRecPt                                          = new TF1("constantRecPt", "[0]",                          0, 50);
        TF1* constantMCPt                                           = new TF1("constantMCPt",  "[0]",                          0, 50);
        TF1* linearRecPt                                            = new TF1("linearRecPt",   "[0]+[1]*x",                    0, 50);
        TF1* linearMCPt                                             = new TF1("linearMCPt",    "[0]+[1]*x",                    0, 50);
        TF1* exponentRecPt                                          = new TF1("exponentRecPt", "[0]-TMath::Exp(-[1]*x+[2])",   0, 50);
        TF1* exponentMCPt                                           = new TF1("exponentMCPt",  "[0]-TMath::Exp(-[1]*x+[2])",   0, 50);
        TF1* powerRecPt                                             = new TF1("powerRecPt",    "[0]/TMath::Power(x,[1])+[2]",           0, 50);
        TF1* powerMCPt                                              = new TF1("powerMCPt",     "[0]/TMath::Power(x,[1])+[2]",           0, 50);

        // taken directly from MC (if statistics sufficient)
        for (Int_t k = 0; k < 3; k++ ){
            histoGammaSecFromXRecoEff_MCPt[k]                       = (TH1D*)fileCorrections->Get(Form(                     "SecondaryGammaFromXFrom%sRecoEff_MCPt",
                                                                                                                            nameSecondaries[k].Data()));
            histoGammaSecFromXRecoEff_MCPt_Unscaled[k]              = (TH1D*)histoGammaSecFromXRecoEff_MCPt[k]->Clone(Form( "SecondaryGammaFromXFrom%sRecoEff_MCPt_Unscaled",
                                                                                                                            nameSecondaries[k].Data()));
            histoGammaSecFromXRecoEff_MCPt_OrBin[k]                 = (TH1D*)fileCorrections->Get(Form(                     "SecondaryGammaFromXFrom%sRecoEff_MCPtOrBin",
                                                                                                                            nameSecondaries[k].Data()));

            histoGammaSecFromXRecoEff_RecPt[k]                      = (TH1D*)fileCorrections->Get(Form(                             "SecondaryGammaFromXFrom%sRecoEff_Pt",
                                                                                                                                    nameSecondaries[k].Data()));
            histoGammaSecFromXRecoEff_RecPt_Unscaled[k]             = (TH1D*)histoGammaSecFromXRecoEff_RecPt[k]->Clone(Form(        "SecondaryGammaFromXFrom%sRecoEff_Pt_Unscaled",
                                                                                                                                    nameSecondaries[k].Data()));
            histoGammaSecFromXRecoEff_RecPt_OrBin[k]                = (TH1D*)fileCorrections->Get(Form(                             "SecondaryGammaFromXFrom%sRecoEff_PtOrBin",
                                                                                                                                    nameSecondaries[k].Data()));
            histoGammaSecFromXRecoEff_RecPt_OrBin_Unscaled[k]       = (TH1D*)histoGammaSecFromXRecoEff_RecPt_OrBin[k]->Clone(Form(  "SecondaryGammaFromXFrom%sRecoEff_PtOrBin_Unscaled",
                                                                                                                                    nameSecondaries[k].Data()));

            ratioGammaSecEffMCPt[k]                                 = (TH1D*)histoGammaSecFromXRecoEff_MCPt[k]->Clone(Form( "RatioSecEffFrom%sToPrim_MCPt", nameSecondaries[k].Data()));
            ratioGammaSecEffMCPt[k]->Divide(histoGammaPrimaryRecoEff_MCPt);
            ratioGammaSecEffRecPt[k]                                = (TH1D*)histoGammaSecFromXRecoEff_RecPt[k]->Clone(Form("RatioSecEffFrom%sToPrim_RecPt", nameSecondaries[k].Data()));
            ratioGammaSecEffRecPt[k]->Divide(histoGammaPrimaryRecoEff_Pt);

            // constant fit to scale primary reco eff (statistics insufficient) MC pt
            ratioGammaSecEffMCPt[k]->Fit(constantMCPt,"SMNR0E+","", minPtFitSec[k], maxPtFitSec[k]);
            constOffsetEffMCPt[k]                                   = constantMCPt->GetParameter(0);

            // exponential fit to scale primary reco eff (statistics insufficient) MC pt
            if (FitSecFunction[k] == kFitSecExp) {
                exponentMCPt->SetParameter(0, constOffsetEffMCPt[k]);
                exponentMCPt->SetParameter(1, 1.);
                exponentMCPt->SetParameter(2, minPtFitSec[k]);
                ratioGammaSecEffMCPt[k]->Fit(exponentMCPt,"SMNR0E+","", minPtFitSec[k], maxPtFitSec[k]);
            }else if(FitSecFunction[k] == kFitSecLinear){
                linearMCPt->SetParameter(1,constOffsetEffMCPt[k]);
                linearMCPt->SetParLimits(1,0.,2.);
                ratioGammaSecEffRecPt[k]->Fit(linearMCPt,"SMNR0E+","", minPtFitSec[k], maxPtFitSec[k]);
            }else if(FitSecFunction[k] == kFitSecPower){
                powerMCPt->SetParameter(0,1.);
                powerMCPt->SetParameter(1,2.);
                powerMCPt->SetParameter(2,constOffsetEffRecPt[k]);
                powerMCPt->SetParLimits(2,0.,2.);
                ratioGammaSecEffRecPt[k]->Fit(powerMCPt,"SMNR0E+","", minPtFitSec[k], maxPtFitSec[k]);
            }

            // constant fit to scale primary reco eff (statistics insufficient) MC pt
            ratioGammaSecEffRecPt[k]->Fit(constantRecPt,"SMNR0E+","", minPtFitSec[k], maxPtFitSec[k]);
            constOffsetEffRecPt[k]                                  = constantRecPt->GetParameter(0);

            // exponential fit to scale primary reco eff (statistics insufficient) MC pt
            if (FitSecFunction[k] == kFitSecExp) {
                exponentRecPt->SetParameter(0, constOffsetEffRecPt[k]);
                exponentRecPt->SetParameter(1, 1.);
                exponentRecPt->SetParameter(2, minPtFitSec[k]);
                ratioGammaSecEffRecPt[k]->Fit(exponentRecPt,"SMNR0E+","", minPtFitSec[k], maxPtFitSec[k]);
            }else if(FitSecFunction[k] == kFitSecLinear){
                linearRecPt->SetParameter(1,constOffsetEffRecPt[k]);
                linearRecPt->SetParLimits(1,0.,2.);
                ratioGammaSecEffRecPt[k]->Fit(linearRecPt,"SMNR0E+","", minPtFitSec[k], maxPtFitSec[k]);
            }else if(FitSecFunction[k] == kFitSecPower){
                powerRecPt->SetParameter(0,1.);
                powerRecPt->SetParameter(1,2.);
                powerRecPt->SetParameter(2,constOffsetEffRecPt[k]);
                powerRecPt->SetParLimits(2,0.,2.);
                ratioGammaSecEffRecPt[k]->Fit(powerRecPt,"SMNR0E+","", minPtFitSec[k], maxPtFitSec[k]);
            }

            // fixed ratios
            if (energy.Contains("2.76TeV") && mode == 4) {
                if (k == 1) {
                    constOffsetEffMCPt[k]                           = 1.35;
                    constOffsetEffRecPt[k]                          = 1.35;
                } else if (k == 2) {
                    constOffsetEffMCPt[k]                           = 1.75;
                    constOffsetEffRecPt[k]                          = 1.75;
                }
            } else if (energy.Contains("PbPb_2.76TeV")) {
                if (k == 2 && constOffsetEffRecPt[k]>1.2) {
                    ratioGammaSecEffRecPt[k]->Fit(constantRecPt,"SMNR0E+","", minPtFitSec[k], 1.8);
                    constOffsetEffRecPt[k]                          = constantRecPt->GetParameter(0);
                    constOffsetEffMCPt[k]                           = constantRecPt->GetParameter(0);
                }
            }

            TF1* tempFunctionRec                            = 0x0;
            TF1* tempFunctionMC                             = 0x0;
            if(FitSecFunction[k] == kFitSecConst){
              tempFunctionRec = constantRecPt;
              tempFunctionMC  = constantMCPt;
            } else if(FitSecFunction[k] == kFitSecLinear){
              tempFunctionRec = linearRecPt;
              tempFunctionMC  = linearMCPt;
            } else if(FitSecFunction[k] == kFitSecExp){
              tempFunctionRec = exponentRecPt;
              tempFunctionMC  = exponentMCPt;
            } else if(FitSecFunction[k] == kFitSecPower){
              tempFunctionRec = powerRecPt;
              tempFunctionMC  = powerMCPt;
            }

            if(!tempFunctionRec || !tempFunctionMC){
              cout << "\n\tERROR: tempFunctionRec or tempFunctionMC not set in L" << __LINE__ << " in CorrectGammaV2.C. Returning..." << endl;
              return;
            }

            // rescale secondary efficiencies from primary efficiency
            Double_t tempScaleRecPt                                 = 1.;
            Double_t tempScaleMCPt                                  = 1.;
            for (Int_t i = histoGammaSecFromXRecoEff_RecPt[k]->FindBin(minPtFitSec[k]); i < histoGammaSecFromXRecoEff_RecPt[k]->GetNbinsX()+1; i++){
                if (FitSecFunction[k] == kFitSecConst) {
                    histoGammaSecFromXRecoEff_RecPt[k]->SetBinContent(  i, histoGammaPrimaryRecoEff_Pt->GetBinContent(i)    * constOffsetEffRecPt[k]);
                    histoGammaSecFromXRecoEff_RecPt[k]->SetBinError(    i, histoGammaPrimaryRecoEff_Pt->GetBinError(i)      * constOffsetEffRecPt[k]);
                    histoGammaSecFromXRecoEff_MCPt[k]->SetBinContent(   i, histoGammaPrimaryRecoEff_MCPt->GetBinContent(i)  * constOffsetEffMCPt[k]);
                    histoGammaSecFromXRecoEff_MCPt[k]->SetBinError(     i, histoGammaPrimaryRecoEff_MCPt->GetBinError(i)    * constOffsetEffMCPt[k]);
                } else {
                    // evaluate fits
                    tempScaleRecPt                                  = tempFunctionRec->Integral(histoGammaSecFromXRecoEff_RecPt[k]->GetXaxis()->GetBinLowEdge(i),
                                                                                                histoGammaSecFromXRecoEff_RecPt[k]->GetXaxis()->GetBinUpEdge(i));
                    tempScaleRecPt                                  = tempScaleRecPt / histoGammaSecFromXRecoEff_RecPt[k]->GetXaxis()->GetBinWidth(i);

                    tempScaleMCPt                                   = tempFunctionMC->Integral(histoGammaSecFromXRecoEff_MCPt[k]->GetXaxis()->GetBinLowEdge(i),
                                                                                               histoGammaSecFromXRecoEff_MCPt[k]->GetXaxis()->GetBinUpEdge(i));
                    tempScaleMCPt                                   = tempScaleMCPt / histoGammaSecFromXRecoEff_MCPt[k]->GetXaxis()->GetBinWidth(i);

                    histoGammaSecFromXRecoEff_RecPt[k]->SetBinContent(  i, histoGammaPrimaryRecoEff_Pt->GetBinContent(i)    * tempScaleRecPt);
                    histoGammaSecFromXRecoEff_RecPt[k]->SetBinError(    i, histoGammaPrimaryRecoEff_Pt->GetBinError(i)      * tempScaleRecPt);
                    histoGammaSecFromXRecoEff_MCPt[k]->SetBinContent(   i, histoGammaPrimaryRecoEff_MCPt->GetBinContent(i)  * tempScaleMCPt);
                    histoGammaSecFromXRecoEff_MCPt[k]->SetBinError(     i, histoGammaPrimaryRecoEff_MCPt->GetBinError(i)    * tempScaleMCPt);
                }
            }
            for (Int_t i = histoGammaSecFromXRecoEff_RecPt_OrBin[k]->FindBin(minPtFitSec[k]); i < histoGammaSecFromXRecoEff_RecPt_OrBin[k]->GetNbinsX()+1; i++){
                if (FitSecFunction[k]==kFitSecConst) {
                    histoGammaSecFromXRecoEff_RecPt_OrBin[k]->SetBinContent( i, histoGammaPrimaryRecoEff_Pt_OrBin->GetBinContent(i)  * constOffsetEffRecPt[k]);
                    histoGammaSecFromXRecoEff_RecPt_OrBin[k]->SetBinError(   i, histoGammaPrimaryRecoEff_Pt_OrBin->GetBinError(i)    * constOffsetEffRecPt[k]);
                    histoGammaSecFromXRecoEff_MCPt_OrBin[k]->SetBinContent(  i, histoGammaPrimaryRecoEff_MCPt_OrBin->GetBinContent(i)* constOffsetEffMCPt[k]);
                    histoGammaSecFromXRecoEff_MCPt_OrBin[k]->SetBinError(    i, histoGammaPrimaryRecoEff_MCPt_OrBin->GetBinError(i)  * constOffsetEffMCPt[k]);
                } else {
                    // evaluate fits
                    tempScaleRecPt                                  = tempFunctionRec->Integral(histoGammaSecFromXRecoEff_RecPt_OrBin[k]->GetXaxis()->GetBinLowEdge(i),
                                                                                                histoGammaSecFromXRecoEff_RecPt_OrBin[k]->GetXaxis()->GetBinUpEdge(i));
                    tempScaleRecPt                                  = tempScaleRecPt / histoGammaSecFromXRecoEff_RecPt_OrBin[k]->GetXaxis()->GetBinWidth(i);

                    tempScaleMCPt                                   = tempFunctionMC->Integral(histoGammaSecFromXRecoEff_MCPt_OrBin[k]->GetXaxis()->GetBinLowEdge(i),
                                                                                               histoGammaSecFromXRecoEff_MCPt_OrBin[k]->GetXaxis()->GetBinUpEdge(i));
                    tempScaleMCPt                                   = tempScaleMCPt / histoGammaSecFromXRecoEff_MCPt_OrBin[k]->GetXaxis()->GetBinWidth(i);

                    histoGammaSecFromXRecoEff_RecPt_OrBin[k]->SetBinContent( i, histoGammaPrimaryRecoEff_Pt_OrBin->GetBinContent(i)  * tempScaleRecPt);
                    histoGammaSecFromXRecoEff_RecPt_OrBin[k]->SetBinError(   i, histoGammaPrimaryRecoEff_Pt_OrBin->GetBinError(i)    * tempScaleRecPt);
                    histoGammaSecFromXRecoEff_MCPt_OrBin[k]->SetBinContent(  i, histoGammaPrimaryRecoEff_MCPt_OrBin->GetBinContent(i)* tempScaleMCPt);
                    histoGammaSecFromXRecoEff_MCPt_OrBin[k]->SetBinError(    i, histoGammaPrimaryRecoEff_MCPt_OrBin->GetBinError(i)  * tempScaleMCPt);
                }
            }

            //******************************************************************************************
            //********************************* plot ratio efficiencies ********************************
            //******************************************************************************************
            TCanvas *canvasSecEffiRatio                             = GetAndSetCanvas("canvasSecEffiRatio");

                Double_t ratioGammaSecEffiToPrimMaxY                = (ratioGammaSecEffRecPt[k]->GetMaximum() + ratioGammaSecEffRecPt[k]->GetBinError(ratioGammaSecEffRecPt[k]->GetMaximumBin()))*1.2;

                SetHistogramm(ratioGammaSecEffRecPt[k],"#it{p}_{T} (GeV/#it{c})","#epsilon_{eff,#gamma, sec}/#epsilon_{eff,#gamma, prim}",0.0,ratioGammaSecEffiToPrimMaxY);
                DrawGammaSetMarker(ratioGammaSecEffRecPt[k], 20, 1, kRed+2, kRed+2);
                ratioGammaSecEffRecPt[k]->Draw();
                DrawGammaSetMarkerTF1( tempFunctionRec, 9, 2, kRed-6);
                tempFunctionRec->SetRange(minPtFitSec[k], maxPtFitSec[k]);
                tempFunctionRec->Draw("same");

                TLegend* legendSecEffRatio = GetAndSetLegend2(0.15,0.93-3*1.1*0.035, 0.4,0.93, 0.035, 1, cent, 42, 0.1);
                legendSecEffRatio->AddEntry(ratioGammaSecEffRecPt[k], Form("sec %s reco. eff/prim",nameLabelSecondaries[k].Data()),"lp");
                legendSecEffRatio->AddEntry(tempFunctionRec, Form("fit: %s",((TString)tempFunctionRec->GetExpFormula("P")).Data()),"l");
                legendSecEffRatio->Draw();

            canvasSecEffiRatio->SaveAs(Form("%s/%s_RatioSecEffiToPrim%sPt_%s_%s.%s",outputDir.Data(),textPi0New.Data(),nameSecondaries[k].Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
            delete canvasSecEffiRatio;

            //******************************************************************************************
            //************************* plot reco effs and fits for secondaries ************************
            //******************************************************************************************
            TCanvas *canvasSecEffiFit                               = GetAndSetCanvas("canvasSecEffiFit");

                if ( isPCM )                    SetHistogramm(histoGammaPrimaryRecoEff_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.02);
                else if ( isCalo && !isPCM )    SetHistogramm(histoGammaPrimaryRecoEff_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 2.0);
                else                            SetHistogramm(histoGammaPrimaryRecoEff_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.02);

                DrawGammaSetMarker(histoGammaPrimaryRecoEff_Pt,         20, 1, kGray+2, kGray+2);
                DrawGammaSetMarker(histoGammaPrimaryRecoEff_Pt_OrBin,   24, 1, kBlack, kBlack);
                histoGammaPrimaryRecoEff_Pt->Draw();
                histoGammaPrimaryRecoEff_Pt_OrBin->Draw("same");

                DrawGammaSetMarker(histoGammaSecFromXRecoEff_RecPt_Unscaled[k],         20, 1, kBlue-8, kBlue-8);
                DrawGammaSetMarker(histoGammaSecFromXRecoEff_RecPt_OrBin_Unscaled[k],   24, 1, kBlue+2, kBlue+2);
                histoGammaSecFromXRecoEff_RecPt_Unscaled[k]->Draw("same");
                histoGammaSecFromXRecoEff_RecPt_OrBin_Unscaled[k]->Draw("same");

                DrawGammaSetMarker(histoGammaSecFromXRecoEff_RecPt[k],       20, 1, kRed-8, kRed-8);
                DrawGammaSetMarker(histoGammaSecFromXRecoEff_RecPt_OrBin[k], 20, 1, kRed+2, kRed+2);
                histoGammaSecFromXRecoEff_RecPt[k]->Draw("same");
                histoGammaSecFromXRecoEff_RecPt_OrBin[k]->Draw("same");

                TLegend* legendSecEffFits = GetAndSetLegend2(0.15,0.93-4*1.1*0.035, 0.7,0.93, 0.035, 2, Form("%s, reco. eff.", cent.Data()), 42, 0.1);
                legendSecEffFits->SetBorderSize(0);
                legendSecEffFits->AddEntry(histoGammaPrimaryRecoEff_Pt_OrBin,                   "primary" ,                                                     "lp");
                legendSecEffFits->AddEntry(histoGammaPrimaryRecoEff_Pt,                         "rebin primary",                                                "lp");
                legendSecEffFits->AddEntry(histoGammaSecFromXRecoEff_RecPt_OrBin_Unscaled[k],   Form("sec. from %s",nameLabelSecondaries[k].Data()),            "lp");
                legendSecEffFits->AddEntry(histoGammaSecFromXRecoEff_RecPt_Unscaled[k],         Form("rebin sec. from %s",nameLabelSecondaries[k].Data()),      "lp");
                legendSecEffFits->AddEntry(histoGammaSecFromXRecoEff_RecPt_OrBin[k],            Form("new sec. from %s",nameLabelSecondaries[k].Data()),        "lp");
                legendSecEffFits->AddEntry(histoGammaSecFromXRecoEff_RecPt[k],                  Form("rebin new sec. from %s",nameLabelSecondaries[k].Data()),  "lp");
                legendSecEffFits->Draw();

            canvasSecEffiFit->SaveAs(Form("%s/%s_SecEffiFits%sPt_%s_%s.%s",outputDir.Data(),textPi0New.Data(),nameSecondaries[k].Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
            delete canvasSecEffiFit;
        }

        delete constantRecPt;
        delete constantMCPt;
        delete linearRecPt;
        delete linearMCPt;
        delete exponentRecPt;
        delete exponentMCPt;
        delete powerRecPt;
        delete powerMCPt;
    }

    //****************************************************************************************** // would be nice if this would be put inside a fct. in the header at some point
    //******************* Secondary gamma conv prob ********************************************
    //******************************************************************************************
    TH1D* histoGammaSecondaryFromXConvProb_MCPt[3]                  = { NULL, NULL, NULL };
    TH1D* histoGammaSecondaryFromXConvProb_MCPt_Unscaled[3]         = { NULL, NULL, NULL };
    TH1D* histoGammaSecondaryFromXConvProb_MCPt_OrBin[3]            = { NULL, NULL, NULL };
    TH1D* histoGammaSecondaryFromXConvProb_MCPt_OrBin_Unscaled[3]   = { NULL, NULL, NULL };
    TH1D* ratioGammaConvProbMCPt[3]                                 = { NULL, NULL, NULL };

    enum        funcFitConvProb {kFitConvProbConst, kFitConvProbLinear, kFitConvProbPower, kFitConvProbPowerOffset};
    funcFitConvProb  FitConvProbFunction[3]                                   = {kFitConvProbPower, kFitConvProbPower, kFitConvProbPower};

    // defaults
    FitConvProbFunction[0]                                          = kFitConvProbPower;
    FitConvProbFunction[1]                                          = kFitConvProbPower;
    FitConvProbFunction[2]                                          = kFitConvProbConst;

    if (energy.Contains("PbPb_2.76TeV")) {
        if(centrality.CompareTo("20-40%")==0 || centrality.CompareTo("20-50%")==0) {
            minPtFitSec[0]                                              = 1.8;
            maxPtFitSec[0]                                              = 50.;
            minPtFitSec[1]                                              = 1.5;
            minPtFitSec[2]                                              = 1.2;
            maxPtFitSec[2]                                              = 3.;
        } else {
            maxPtFitSec[0]                                              = 50.;
            minPtFitSec[1]                                              = 1.5;
            maxPtFitSec[1]                                              = 50.;
            maxPtFitSec[2]                                              = 50.;
        }
    }else if(energy.Contains("7TeV") && mode == 2){
        FitConvProbFunction[0]                                          = kFitConvProbPowerOffset;
        minPtFitSec[0]                                                  = 1.0;
        maxPtFitSec[0]                                                  = 16.;
        FitConvProbFunction[1]                                          = kFitConvProbPowerOffset;
        minPtFitSec[1]                                                  = 0.8;
        maxPtFitSec[1]                                                  = 6.;
        minPtFitSec[2]                                                  = 0.8;
        maxPtFitSec[2]                                                  = 4.;
    }else if(energy.Contains("8TeV") && mode == 2){
        FitConvProbFunction[0]                                          = kFitConvProbPowerOffset;
        minPtFitSec[0]                                                  = 1.0;
        maxPtFitSec[0]                                                  = 16.;
        FitConvProbFunction[1]                                          = kFitConvProbPowerOffset;
    } else if(energy.CompareTo("pPb_5.023TeVRun2") == 0 && mode == 0){
        minPtFitSec[2]                                                  = 0.9;
        maxPtFitSec[2]                                                  = 3.0;
    }

    if ( hasCocktailInput && isPCM ) {
        TF1* powerCP                                                = new TF1("powerCP",       "[0]/TMath::Power((x-[1]), [2])", 0, 50);
        TF1* constantCP                                             = new TF1("constantCP",    "[0]",                            0, 50);
        TF1* linearCP                                               = new TF1("linearCP",      "[0]+[1]*x",                      0, 50);
        TF1* powerOffsetCP                                          = new TF1("powerOffsetCP", "[0]/TMath::Power(x,[1])+[2]",    0, 50);

        for (Int_t k = 0; k < 3; k++) {
            histoGammaSecondaryFromXConvProb_MCPt[k]                = (TH1D*)fileCorrections->Get(Form(                                 "SecondaryGammaFromXFrom%sConvProb_MCPt",
                                                                                                                                        nameSecondaries[k].Data()));
            histoGammaSecondaryFromXConvProb_MCPt_OrBin[k]          = (TH1D*)fileCorrections->Get(Form(                                 "SecondaryGammaFromXFrom%sConvProb_MCPtOrBin",
                                                                                                                                        nameSecondaries[k].Data()));
            histoGammaSecondaryFromXConvProb_MCPt_Unscaled[k]       = (TH1D*)histoGammaSecondaryFromXConvProb_MCPt[k]->Clone(Form(      "SecondaryGammaFromXFrom%sConvProb_MCPt_Unscaled",
                                                                                                                                        nameSecondaries[k].Data()));
            histoGammaSecondaryFromXConvProb_MCPt_OrBin_Unscaled[k] = (TH1D*)histoGammaSecondaryFromXConvProb_MCPt_OrBin[k]->Clone(Form("SecondaryGammaFromXFrom%sConvProb_MCPtOrBin_Unscaled",
                                                                                                                                        nameSecondaries[k].Data()));

            ratioGammaConvProbMCPt[k]                               = (TH1D*)histoGammaSecondaryFromXConvProb_MCPt[k]->Clone(Form(      "RatioConvProbFrom%sToPrim_MCPt",
                                                                                                                                        nameSecondaries[k].Data()));

            ratioGammaConvProbMCPt[k]->Divide(histoGammaConvProb_MCPt);

            ratioGammaConvProbMCPt[k]->Fit(constantCP,"SMNR0E+","", minPtFitSec[k], maxPtFitSec[k]);
            if(FitConvProbFunction[k] == kFitConvProbPower){
              if (energy.Contains("PbPb_2.76TeV") && k==1) powerCP->FixParameter(2,1);
              if (energy.CompareTo("2.76TeV") == 0 && mode == 2 && k ==2 ){
                  powerCP->FixParameter(0,0.3);
              }
              ratioGammaConvProbMCPt[k]->Fit(powerCP,"SMNR0E+", "", minPtFitSec[k], maxPtFitSec[k]);
            }else if(FitConvProbFunction[k] == kFitConvProbConst){
              ratioGammaConvProbMCPt[k]->Fit(constantCP,"SMNR0E+", "", minPtFitSec[k], maxPtFitSec[k]);
            }else if(FitConvProbFunction[k] == kFitConvProbLinear){
              ratioGammaConvProbMCPt[k]->Fit(linearCP,"SMNR0E+", "", minPtFitSec[k], maxPtFitSec[k]);
            }else if(FitConvProbFunction[k] == kFitConvProbPowerOffset){
              powerOffsetCP->SetParameter(0,1.);
              powerOffsetCP->SetParameter(1,2.);
              powerOffsetCP->SetParameter(2,constantCP->GetParameter(0));
              powerOffsetCP->SetParLimits(2,-1.,1.);
              ratioGammaConvProbMCPt[k]->Fit(powerOffsetCP,"SMNR0E+", "", minPtFitSec[k], maxPtFitSec[k]);
            }

            TF1* tempFunction                            = 0x0;
            if(FitConvProbFunction[k] == kFitConvProbPower){
              tempFunction = powerCP;
            } else if(FitConvProbFunction[k] == kFitConvProbConst){
              tempFunction = constantCP;
            } else if(FitConvProbFunction[k] == kFitConvProbLinear){
              tempFunction = linearCP;
            } else if(FitConvProbFunction[k] == kFitConvProbPowerOffset){
              tempFunction = powerOffsetCP;
            }

            if(!tempFunction){
              cout << "\n\tERROR: tempFunction not set in L" << __LINE__ << " in CorrectGammaV2.C. Returning..." << endl;
              return;
            }

            Double_t tempEval                                       = 1.;
            for (Int_t i=histoGammaSecondaryFromXConvProb_MCPt[k]->FindBin(minPtFitSec[k]); i<histoGammaSecondaryFromXConvProb_MCPt[k]->GetNbinsX()+1; i++) {
                // was: tempEval = powerCP->Eval(histoGammaSecondaryFromXConvProb_MCPt[k]->GetBinCenter(i))
                tempEval                                            = tempFunction->Integral(histoGammaSecondaryFromXConvProb_MCPt[k]->GetXaxis()->GetBinLowEdge(i),
                                                                                      histoGammaSecondaryFromXConvProb_MCPt[k]->GetXaxis()->GetBinUpEdge(i));
                tempEval                                            = tempEval / histoGammaSecondaryFromXConvProb_MCPt[k]->GetBinWidth(i);

                histoGammaSecondaryFromXConvProb_MCPt[k]->SetBinContent(i, histoGammaConvProb_MCPt->GetBinContent(i)    * tempEval);
                histoGammaSecondaryFromXConvProb_MCPt[k]->SetBinError(  i, histoGammaConvProb_MCPt->GetBinError(i)      * tempEval);
            }
            for (Int_t i=histoGammaSecondaryFromXConvProb_MCPt_OrBin[k]->FindBin(minPtFitSec[k]); i<histoGammaSecondaryFromXConvProb_MCPt_OrBin[k]->GetNbinsX()+1; i++) {
                // was: tempEval = powerCP->Eval(histoGammaSecondaryFromXConvProb_MCPt_OrBin[k]->GetBinCenter(i))
                tempEval                                            = tempFunction->Integral(histoGammaSecondaryFromXConvProb_MCPt_OrBin[k]->GetXaxis()->GetBinLowEdge(i),
                                                                                      histoGammaSecondaryFromXConvProb_MCPt_OrBin[k]->GetXaxis()->GetBinUpEdge(i));
                tempEval                                            = tempEval / histoGammaSecondaryFromXConvProb_MCPt_OrBin[k]->GetBinWidth(i);

                histoGammaSecondaryFromXConvProb_MCPt_OrBin[k]->SetBinContent(i, histoGammaConvProb_MCPt_OrBin->GetBinContent(i) * tempEval);
                histoGammaSecondaryFromXConvProb_MCPt_OrBin[k]->SetBinError(  i, histoGammaConvProb_MCPt_OrBin->GetBinError(i)   * tempEval);
            }
            if (energy.Contains("PbPb_2.76TeV") && k==1) tempFunction->ReleaseParameter(2);

            //******************************************************************************************
            //************************* plot ratio conversion probabilies ******************************
            //******************************************************************************************
            TCanvas *canvasSecConvProbRatio                         = GetAndSetCanvas("canvasSecConvProbRatio");

                SetHistogramm(ratioGammaConvProbMCPt[k],"#it{p}_{T} (GeV/#it{c})","P_{conv, sec}/P_{conv, prim}",0.0,3.0);
                DrawGammaSetMarker(ratioGammaConvProbMCPt[k], 20, 1, kRed+2, kRed+2);
                ratioGammaConvProbMCPt[k]->Draw();
                DrawGammaSetMarkerTF1( tempFunction, 9, 2, kRed-6);
                tempFunction->SetRange(minPtFitSec[k], maxPtFitSec[k]);
                tempFunction->Draw("same");

                TLegend* legendSecConvProbRatio                     = GetAndSetLegend2(0.15,0.93-3*1.1*0.035, 0.4,0.93, 0.035, 1, cent, 42, 0.1);
                legendSecConvProbRatio->AddEntry(ratioGammaConvProbMCPt[k], Form("sec %s P_{conv}/prim",nameLabelSecondaries[k].Data()),"lp");
                legendSecConvProbRatio->AddEntry(tempFunction, Form("fit: %s",((TString)tempFunction->GetExpFormula("P")).Data()),"l");
                legendSecConvProbRatio->Draw();

            canvasSecConvProbRatio->SaveAs(Form("%s/%s_RatioSecConvProbToPrim%sPt_%s_%s.%s",outputDir.Data(),textPi0New.Data(),nameSecondaries[k].Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
            delete canvasSecConvProbRatio;


            //******************************************************************************************
            //*************** plot conversion probabilities for secondaries and fits *******************
            //******************************************************************************************
            TCanvas *canvasSecConvProbFit                           = GetAndSetCanvas("canvasSecConvProbFit");
            canvasSecConvProbFit->SetTopMargin(0.035);

                SetHistogramm(histoGammaConvProb_MCPt,"#it{p}_{T} (GeV/#it{c})","P_{conv} in |#eta| < 0.9",0.0,1.5e-1);
                DrawGammaSetMarker(histoGammaConvProb_MCPt, 20, 1, kGray+2, kGray+2);
                histoGammaConvProb_MCPt->Draw();

                DrawGammaSetMarker(histoGammaConvProb_MCPt_OrBin, 24, 1, kBlack, kBlack);
                histoGammaConvProb_MCPt_OrBin->Draw("same");

                DrawGammaSetMarker(histoGammaSecondaryFromXConvProb_MCPt_Unscaled[k],       20, 1, kBlue-8, kBlue-8);
                DrawGammaSetMarker(histoGammaSecondaryFromXConvProb_MCPt_OrBin_Unscaled[k], 24, 1, kBlue+2, kBlue+2);
                DrawGammaSetMarker(histoGammaSecondaryFromXConvProb_MCPt[k],                20, 1, kRed-8, kRed-8);
                DrawGammaSetMarker(histoGammaSecondaryFromXConvProb_MCPt_OrBin[k],          24, 1, kRed+2, kRed+2);

                histoGammaSecondaryFromXConvProb_MCPt_Unscaled[k]->Draw("same");
                histoGammaSecondaryFromXConvProb_MCPt_OrBin_Unscaled[k]->Draw("same");
                histoGammaSecondaryFromXConvProb_MCPt[k]->Draw("same");
                histoGammaSecondaryFromXConvProb_MCPt_OrBin[k]->Draw("same");

                TLegend* legendSecConvProbFits = GetAndSetLegend2(0.15,0.93-4*1.1*0.035, 0.7,0.93, 0.035, 2, Form("%s, conv. prob.", cent.Data()), 42, 0.1);
                legendSecConvProbFits->SetBorderSize(0);
                legendSecConvProbFits->AddEntry(histoGammaConvProb_MCPt_OrBin,                              "primary",                                                      "lp");
                legendSecConvProbFits->AddEntry(histoGammaConvProb_MCPt,                                    "rebin primary",                                                "lp");
                legendSecConvProbFits->AddEntry(histoGammaSecondaryFromXConvProb_MCPt_OrBin_Unscaled[k],    Form("sec. from %s",nameLabelSecondaries[k].Data()),            "lp");
                legendSecConvProbFits->AddEntry(histoGammaSecondaryFromXConvProb_MCPt_Unscaled[k],          Form("rebin sec. from %s",nameLabelSecondaries[k].Data()),      "lp");
                legendSecConvProbFits->AddEntry(histoGammaSecondaryFromXConvProb_MCPt_OrBin[k],             Form("new sec. from %s",nameLabelSecondaries[k].Data()),        "lp");
                legendSecConvProbFits->AddEntry(histoGammaSecondaryFromXConvProb_MCPt[k],                   Form("rebin new sec. from %s",nameLabelSecondaries[k].Data()),  "lp");
                legendSecConvProbFits->Draw();

            canvasSecConvProbFit->SaveAs(Form("%s/%s_ConversionProbFits%sMCPt_%s_%s.%s",outputDir.Data(),textPi0New.Data(),nameSecondaries[k].Data(), nameRec.Data(),cutSelection.Data(),suffix.Data()));
            delete canvasSecConvProbFit;
        }

        delete powerCP;
        delete constantCP;
        delete powerOffsetCP;
        delete linearCP;
    }

    //******************************************************************************************
    //******************* Secondary gamma response matrices ************************************
    //******************************************************************************************
    TH2D* histoGammaTrueSecondaryFromX_MCPt_recPt[3]                = { NULL, NULL, NULL };
    TH2D* histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[3]          = { NULL, NULL, NULL };
    // PCM
    if ( hasCocktailInput && isPCM ) {
        for (Int_t k = 0; k <3; k++){
            histoGammaTrueSecondaryFromX_MCPt_recPt[k]              = (TH2D*)fileCorrections->Get(Form("TrueSecondaryConvGammaFromXFrom%s_MCPt_recPt",nameSecondaries[k].Data()));
            histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[k]        = (TH2D*)fileCorrections->Get(Form("TrueSecondaryConvGammaFromXFrom%s_MCPt_recPt_orBin",nameSecondaries[k].Data()));
        }
    }
    // Calo
    if ( hasCocktailInput && isCalo && !isPCM ) {
        for (Int_t k = 0; k <3; k++){
            histoGammaTrueSecondaryFromX_MCPt_recPt[k]              = (TH2D*)fileCorrections->Get(Form("TrueSecondaryCaloGammaFromXFrom%s_MCPt_recPt",nameSecondaries[k].Data()));
            histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[k]        = (TH2D*)fileCorrections->Get(Form("TrueSecondaryCaloGammaFromXFrom%s_MCPt_recPt_orBin",nameSecondaries[k].Data()));
        }
    }

    //******************************************************************************************
    //******************* Calculate raw secondary spectra from cocktail input ******************
    //******************************************************************************************
    TH1D* histoGammaSecGammaFromX_Cocktail_Raw_Pt[4]                = { NULL, NULL, NULL, NULL };
    TH1D* histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin[4]          = { NULL, NULL, NULL, NULL };
    if ( hasCocktailInput && (isPCM || isCalo) ) {
        cout << "calculating raw secondary spectra from cocktail" << endl;
        // K0s: clone cocktail spectra
        for (Int_t k = 0; k < 3; k++){
            histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]              = (TH1D*) histoGammaTrueSecCocktailGammaFromX_Pt[k]->Clone(Form("histoGammaTrueSecGammaFromXFrom%s_Cocktail_Raw_Pt",
                                                                                                                                    nameSecondaries[k].Data()));
            histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]->Sumw2();
            histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin[k]        = (TH1D*) histoGammaTrueSecCocktailGammaFromX_PtOrBin[k]->Clone(Form("histoGammaTrueSecGammaFromXFrom%s_Cocktail_Raw_PtOrBin",
                                                                                                                                         nameSecondaries[k].Data()));
            histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin[k]->Sumw2();
        }

        if ( isPCM ) {
            // K0s: calculate raw yield
            if (useUnfoldingForCocktailSecondaries) {
                hasCocktailInput                                    = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt[0], histoGammaSecondaryFromXConvProb_MCPt[0],
                                                                                                    histoGammaSecFromXRecoEff_MCPt[0], histoGammaTrueSecondaryFromX_MCPt_recPt[0], nEvt,
                                                                                                    kTRUE, nIterationsUnfolding);

                hasCocktailInput                                    = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin[0], histoGammaSecondaryFromXConvProb_MCPt_OrBin[0],
                                                                                                    histoGammaSecFromXRecoEff_MCPt_OrBin[0], histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[0], nEvt,
                                                                                                    kTRUE, nIterationsUnfolding);
            } else {
                hasCocktailInput                                    = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt[0], histoGammaSecondaryFromXConvProb_MCPt[0],
                                                                                                    histoGammaSecFromXRecoEff_RecPt[0], histoGammaTrueSecondaryFromX_MCPt_recPt[0], nEvt,
                                                                                                    kFALSE);

                hasCocktailInput                                    = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin[0], histoGammaSecondaryFromXConvProb_MCPt_OrBin[0],
                                                                                                    histoGammaSecFromXRecoEff_RecPt_OrBin[0], histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[0], nEvt,
                                                                                                    kFALSE);
            }
            // K0l and Lambda
            for (Int_t k = 1; k < 3; k++){
                hasCocktailInput                                    = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k], histoGammaSecondaryFromXConvProb_MCPt[k],
                                                                                                    histoGammaSecFromXRecoEff_RecPt[k], histoGammaTrueSecondaryFromX_MCPt_recPt[k], nEvt,
                                                                                                    kFALSE);

                hasCocktailInput                                    = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin[k], histoGammaSecondaryFromXConvProb_MCPt_OrBin[k],
                                                                                                    histoGammaSecFromXRecoEff_RecPt_OrBin[k], histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[k], nEvt,
                                                                                                    kFALSE);
            }
        }

        if ( isCalo && !isPCM ) {
            // K0s: calculate raw yield
            if (useUnfoldingForCocktailSecondaries) {
                hasCocktailInput                                    = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt[0], histoGammaSecFromXRecoEff_MCPt[0],
                                                                                                    histoGammaTrueSecondaryFromX_MCPt_recPt[0], nEvt, kTRUE, nIterationsUnfolding);
                hasCocktailInput                                    = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin[0], histoGammaSecFromXRecoEff_MCPt_OrBin[0],
                                                                                                    histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[0], nEvt, kTRUE, nIterationsUnfolding);
            } else {
                hasCocktailInput                                    = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt[0], histoGammaSecFromXRecoEff_RecPt[0],
                                                                                                    histoGammaTrueSecondaryFromX_MCPt_recPt[0], nEvt, kFALSE);
                hasCocktailInput                                    = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin[0], histoGammaSecFromXRecoEff_RecPt_OrBin[0],
                                                                                                    histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[0], nEvt, kFALSE);
            }
            // K0l and Lambda
            for (Int_t k = 1; k < 3; k++){
                hasCocktailInput                                    = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k], histoGammaSecFromXRecoEff_RecPt[k],
                                                                                                    histoGammaTrueSecondaryFromX_MCPt_recPt[k], nEvt, kFALSE);
                hasCocktailInput                                    = ConvertCocktailSecondaryToRaw(histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin[k], histoGammaSecFromXRecoEff_RecPt_OrBin[k],
                                                                                                    histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[k], nEvt, kFALSE);
            }
        }
    }

    //******************************************************************************************
    //******************* Calculate secondary fractions from cocktail input ********************
    //******************************************************************************************
    TH1D* histoFracAllGammaToSecFromX_Cocktail_Pt[4]                = { NULL, NULL, NULL, NULL };
    TH1D* histoFracAllGammaToSecFromX_Cocktail_PtOrBin[3]           = { NULL, NULL, NULL };
    if ( hasCocktailInput && isPCM ) {
        cout << "calculating secondary fractions from cocktail" << endl;
        for (Int_t k = 0; k < 3; k++){
            histoFracAllGammaToSecFromX_Cocktail_Pt[k]              = (TH1D*)histoESDConvGammaPt->Clone(Form("FracAllGammaToSecFromXFrom%s", nameSecondaries[k].Data()));
            histoFracAllGammaToSecFromX_Cocktail_Pt[k]->Divide(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k],histoFracAllGammaToSecFromX_Cocktail_Pt[k],1,1,"B");
            histoFracAllGammaToSecFromX_Cocktail_PtOrBin[k]         = (TH1D*)histoESDConvGammaPt_OrBin->Clone(Form("FracAllGammaToSecFromXFrom%sOriginalBinning", nameSecondaries[k].Data()));
            histoFracAllGammaToSecFromX_Cocktail_PtOrBin[k]->Divide(histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin[k],histoFracAllGammaToSecFromX_Cocktail_PtOrBin[k],1,1,"B");
        }

        histoFracAllGammaToSecFromX_Cocktail_Pt[3]                  = (TH1D*)histoESDConvGammaPt->Clone("FracAllGammaToSecRest");
        histoFracAllGammaToSecFromX_Cocktail_Pt[3]->Divide(histoGammaTrueSecCocktailGammaFromX_Pt[3],histoFracAllGammaToSecFromX_Cocktail_Pt[3],1,1,"B");
    }
    if ( hasCocktailInput && isCalo && !isPCM ) {
        cout << "calculating secondary fractions from cocktail" << endl;
        for (Int_t k = 0; k < 3; k++){
            histoFracAllGammaToSecFromX_Cocktail_Pt[k]              = (TH1D*)histoESDCaloGammaPt->Clone(Form("FracAllGammaToSecFromXFrom%s", nameSecondaries[k].Data()));
            histoFracAllGammaToSecFromX_Cocktail_Pt[k]->Divide(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k],histoFracAllGammaToSecFromX_Cocktail_Pt[k],1,1,"B");
            histoFracAllGammaToSecFromX_Cocktail_PtOrBin[k]         = (TH1D*)histoESDCaloGammaPt_OrBin->Clone(Form("FracAllGammaToSecFromXFrom%sOriginalBinning", nameSecondaries[k].Data()));
            histoFracAllGammaToSecFromX_Cocktail_PtOrBin[k]->Divide(histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin[k],histoFracAllGammaToSecFromX_Cocktail_PtOrBin[k],1,1,"B");
        }

        histoFracAllGammaToSecFromX_Cocktail_Pt[3]                  = (TH1D*)histoESDCaloGammaPt->Clone("FracAllGammaToSecRest");
        histoFracAllGammaToSecFromX_Cocktail_Pt[3]->Divide(histoGammaTrueSecCocktailGammaFromX_Pt[3],histoFracAllGammaToSecFromX_Cocktail_Pt[3],1,1,"B");
    }

    //******************************************************************************************
    //******************* Pileup correction factors ********************************************
    //******************************************************************************************
    TH1D*   histoESDConvGammaPtPileUp                               = NULL;
    TH1D*   histoRatioWithWithoutPileUp                             = NULL;
    TF1*    histoRatioWithWithoutPileUpFit                          = NULL;
    TH1D*   histoPileUpCorrectionFactor_Pt                          = NULL;
    TH1D*   histoPileUpCorrectionFactor_Pt_OrBin                    = NULL;
    TH1D*   histoPileUpCorrectionFactorNoFit_Pt                     = NULL;
    TGraphAsymmErrors* graphGammaSysErrOOBPileupDown                = NULL;
    TGraphAsymmErrors* graphGammaSysErrOOBPileupUp                  = NULL;
    if(doPileUpCorr && isPCM){
        graphGammaSysErrOOBPileupDown                               = (TGraphAsymmErrors*)fileUnCorrected->Get("Gamma_OOBPileupSysDown");
        graphGammaSysErrOOBPileupUp                                 = (TGraphAsymmErrors*)fileUnCorrected->Get("Gamma_OOBPileupSysUp");
        histoESDConvGammaPtPileUp                                   = (TH1D*)fileUnCorrected->Get("ESD_ConvGamma_Pt_PileUp");

        // ratio raw yield to raw yield after pileup BG subtraction
        histoRatioWithWithoutPileUp                                 = (TH1D*)histoESDConvGammaPt->Clone("histoRatioWithWithoutPileUp");
        histoRatioWithWithoutPileUp->Divide(histoRatioWithWithoutPileUp,histoESDConvGammaPtPileUp,1,1,"B");

        // pileup correction factor directly from ratio to cross check correction factors from fit
        histoPileUpCorrectionFactorNoFit_Pt                         = (TH1D*)histoESDConvGammaPtPileUp->Clone("PileUpCorrectionFactorNoFit");
        histoPileUpCorrectionFactorNoFit_Pt->Divide(histoPileUpCorrectionFactorNoFit_Pt,histoESDConvGammaPt,1,1,"B");

        // fit correction factor to get back to original binning
        cout << "fitting ratio gamma raw yield to raw yield after pileup subtraction to extract the pileup correction factor" << endl;
        //otherwise function is null
        Double_t rangeShift = 0.;
        if( energy.Contains("PbPb") ||
            (energy.CompareTo("pPb_5.023TeVRun2")==0 ) ||

            (energy.BeginsWith("8TeV") && mode == 2))
            rangeShift = 0.5;
        else if ((energy.CompareTo("pPb_5.023TeV")==0 && centrality.CompareTo("0-100%") != 0 && mode == 2) )
            rangeShift = 0.8;

        Int_t   fitStatus                                           = 0;
                histoRatioWithWithoutPileUpFit                      = new TF1("histoRatioWithWithoutPileUpFit", "1+[0]/TMath::Power((x-[1]), [2])",
                                                                              histoESDConvGammaPt_OrBin->GetXaxis()->GetXmin()+rangeShift,
                                                                              histoESDConvGammaPt_OrBin->GetXaxis()->GetXmax());
        histoRatioWithWithoutPileUpFit->SetParameters(1, 0, 1);
        histoRatioWithWithoutPileUpFit->SetName(Form("%s_fit", histoRatioWithWithoutPileUp->GetName()));
        TFitResultPtr fitResultPileup                               = histoRatioWithWithoutPileUp->Fit(histoRatioWithWithoutPileUpFit, "SMNRE+","",
                                                                                                       histoRatioWithWithoutPileUp->GetXaxis()->GetXmin()+rangeShift,
                                                                                                       histoRatioWithWithoutPileUp->GetXaxis()->GetXmax());
        fitStatus                                                   = fitResultPileup;
        // accepting fits, if everything went fine (i.e. fitstatus = 0) of if only improve had problems (i.e. fitstatus >= 1000)
        cout << "fit status: " << fitStatus << endl;

        if (fitStatus == 0 || fitStatus >= 1000 || fitStatus == 4) {

            // pileup correction factor in analysis binning
            histoPileUpCorrectionFactor_Pt                          = (TH1D*)histoESDConvGammaPt->Clone("PileUpCorrectionFactorOrBin");
            histoPileUpCorrectionFactor_Pt->Reset("ICES");
            histoPileUpCorrectionFactor_Pt->Sumw2();

            Double_t binContent                                     = 0.;
            Double_t binError                                       = 0.;
            for (Int_t i=1; i<histoPileUpCorrectionFactor_Pt->GetNbinsX()+1; i++) {
                histoRatioWithWithoutPileUpFit->SetParameters(fitResultPileup->GetParams());
                binContent                                          = histoRatioWithWithoutPileUpFit->Integral(histoPileUpCorrectionFactor_Pt->GetXaxis()->GetBinLowEdge(i),
                                                                                                               histoPileUpCorrectionFactor_Pt->GetXaxis()->GetBinUpEdge(i)) /
                                                                                                               histoPileUpCorrectionFactor_Pt->GetBinWidth(i);
                binError                                            = histoRatioWithWithoutPileUpFit->IntegralError(histoPileUpCorrectionFactor_Pt->GetXaxis()->GetBinLowEdge(i),
                                                                                                                    histoPileUpCorrectionFactor_Pt->GetXaxis()->GetBinUpEdge(i),
                                                                                                                    fitResultPileup->GetParams(),
                                                                                                                    fitResultPileup->GetCovarianceMatrix().GetMatrixArray()) /
                                                                                                                    histoPileUpCorrectionFactor_Pt->GetBinWidth(i);
                if (binContent > 0) {
                    histoPileUpCorrectionFactor_Pt->SetBinContent(i, 1./binContent);
                    histoPileUpCorrectionFactor_Pt->SetBinError(  i, binError/binContent/binContent);
                } else {
                    histoPileUpCorrectionFactor_Pt->SetBinContent(i, 1.);
                    histoPileUpCorrectionFactor_Pt->SetBinError(  i, 0.);
                }
            }

            // pileup correction factor in original binning
            histoPileUpCorrectionFactor_Pt_OrBin                    = (TH1D*)histoESDConvGammaPt_OrBin->Clone("PileUpCorrectionFactorOrBin");
            histoPileUpCorrectionFactor_Pt_OrBin->Reset("ICES");
            histoPileUpCorrectionFactor_Pt_OrBin->Sumw2();

            binContent                                              = 0.;
            binError                                                = 0.;
            for (Int_t i=1; i<histoPileUpCorrectionFactor_Pt_OrBin->GetNbinsX()+1; i++) {
                histoRatioWithWithoutPileUpFit->SetParameters(fitResultPileup->GetParams());
                binContent                                          = histoRatioWithWithoutPileUpFit->Integral(histoPileUpCorrectionFactor_Pt_OrBin->GetXaxis()->GetBinLowEdge(i),
                                                                                                               histoPileUpCorrectionFactor_Pt_OrBin->GetXaxis()->GetBinUpEdge(i)) /
                                                                                                               histoPileUpCorrectionFactor_Pt_OrBin->GetBinWidth(i);
                binError                                            = histoRatioWithWithoutPileUpFit->IntegralError(histoPileUpCorrectionFactor_Pt_OrBin->GetXaxis()->GetBinLowEdge(i),
                                                                                                                    histoPileUpCorrectionFactor_Pt_OrBin->GetXaxis()->GetBinUpEdge(i),
                                                                                                                    fitResultPileup->GetParams(),
                                                                                                                    fitResultPileup->GetCovarianceMatrix().GetMatrixArray()) /
                                                                                                                    histoPileUpCorrectionFactor_Pt_OrBin->GetBinWidth(i);
                if (binContent > 0) {
                    histoPileUpCorrectionFactor_Pt_OrBin->SetBinContent(i, 1./binContent);
                    histoPileUpCorrectionFactor_Pt_OrBin->SetBinError(  i, binError/binContent/binContent);
                } else {
                    histoPileUpCorrectionFactor_Pt_OrBin->SetBinContent(i, 1.);
                    histoPileUpCorrectionFactor_Pt_OrBin->SetBinError(  i, 0.);
                }
            }
        } else {
            cout << "ERROR! Pileup correction factor fit failed, no pileup correction will be applied!" << endl;
            doPileUpCorr                                            = 0;
        }
    }

    //******************************************************************************************
    //******************* MC pileup histograms *************************************************
    //******************************************************************************************
    TH1D*   histoGammaPurity_PileUp_Pt                              = NULL;
    TH1D*   histoGammaTruePurity_PileUp_Pt                          = NULL;
    TH1D*   histoGammaRecoEff_PileUp_Pt                             = NULL;
    TH1D*   histoGammaPrimaryRecoEff_PileUp_Pt                      = NULL;
    if( doPileUpCorr && isPCM ){
        // (true primary + true secondary conv. gamma)_(after MC pileup corr.) / (rec. MC conv. gamma)_(after MC pileup corr.)
        histoGammaPurity_PileUp_Pt                                  = (TH1D*)fileCorrections->Get("GammaPurity_PileUp_Pt");
        // (true primary conv. gamma)_(after MC pileup corr.) / (rec. MC conv. gamma - true secondary conv. gamma)_(after MC pileup corr.)
        histoGammaTruePurity_PileUp_Pt                              = (TH1D*)fileCorrections->Get("GammaTruePurity_PileUp_Pt");
        // (true primary + true secondary conv. gamma)_(after MC pileup corr.) / (MC conv. gamma)
        histoGammaRecoEff_PileUp_Pt                                 = (TH1D*)fileCorrections->Get("GammaRecoEff_PileUp_Pt");
        // (true primary conv. gamma)_(after MC pileup corr.) / (MC conv. gamma)
        histoGammaPrimaryRecoEff_PileUp_Pt                          = (TH1D*)fileCorrections->Get("GammaPrimaryRecoEff_PileUp_Pt");
    }
    TH1D*   histoMCrecGamma_PileUp_Pt                               = NULL;
    TH1D*   histoPileUpCorrectionFactorMC_Pt                        = NULL;
    if(doPileUpCorr && isPCM){
        histoMCrecGamma_PileUp_Pt                                   = (TH1D*)fileCorrections->Get("MCrec_ConvGamma_Pt_PileUp");
        if (histoMCrecGamma_PileUp_Pt){
            // MC pileup correction factor (there is no pileup in MC, used as cross check/to disentangle effect from real pileup and prim. vertex resolution)
            histoPileUpCorrectionFactorMC_Pt                        = (TH1D*)histoMCrecGamma_PileUp_Pt->Clone("PileUpCorrectionFactorMC");
            histoPileUpCorrectionFactorMC_Pt->Divide(histoPileUpCorrectionFactorMC_Pt,histoMCrecGamma_Pt,1,1,"B");
        }
    }

    //******************************************************************************************
    //******************* Proper scaling of background *****************************************
    //******************************************************************************************
    TH1D *ScalingGammaBackground_Pt                                 = NULL;
    if (isPCM ) {
        ScalingGammaBackground_Pt                                   = (TH1D*) histoESDConvGammaPt->Clone("ScalingGammaBackground_Pt");
        ScalingGammaBackground_Pt->Divide(ScalingGammaBackground_Pt, histoMCrecGamma_Pt, 1., 1, "");
    }
    if (isCalo && !isPCM) {
        ScalingGammaBackground_Pt                                   = (TH1D*) histoESDCaloGammaPt->Clone("ScalingGammaBackground_Pt");
        ScalingGammaBackground_Pt->Divide(ScalingGammaBackground_Pt, histoMCrecGammaCalo_Pt, 1., 1, "");
    }
    ScalingGammaBackground_Pt->Scale(nEvtMC/nEvt);
    histoMCrecBackground_Pt->Scale(1./nEvtMC);
    TH1D *histoGammaMCBackground_Pt                                 = (TH1D*)histoMCrecBackground_Pt->Clone("histoGammaMCBackground_Pt");
    histoMCrecBackground_Pt->Multiply(ScalingGammaBackground_Pt);

    //******************************************************************************************
    //******************* Calculate secondary spectra from data ********************************
    //******************************************************************************************
    TH1D* histoSecondaryGammaFromXSpecPt[4]                         = { NULL, NULL, NULL, NULL };
    TH1D* histoSecondaryGammaFromXSpecPtOrBin[4]                    = { NULL, NULL, NULL, NULL };
    TH1D* histoSecondaryGammaFromXSpecPileUpPt[4]                   = { NULL, NULL, NULL, NULL };
    if (isPCM ) {
        for (Int_t k = 0; k < 4; k++){
            histoSecondaryGammaFromXSpecPt[k]                       = (TH1D*)histoESDConvGammaPt->Clone(        Form("SecondaryGammaSpecFromXFrom%sPt",nameSecondaries[k].Data()));
            histoSecondaryGammaFromXSpecPtOrBin[k]                  = (TH1D*)histoESDConvGammaPt_OrBin->Clone(  Form("SecondaryGammaSpecFromXFrom%sPt_OrBin",nameSecondaries[k].Data()));
            if (doPileUpCorr)
                histoSecondaryGammaFromXSpecPileUpPt[k]             = (TH1D*)histoESDConvGammaPtPileUp->Clone(  Form("SecondaryGammaSpecFromXFrom%sPileUpPt",nameSecondaries[k].Data()));

            // overwrite spectra loaded from pure MC if running data
            if (!isRunMC) histoGammaTrueSecConvGammaFromX_Pt[k]     = (TH1D*)histoESDConvGammaPt->Clone(        Form("SecondaryMCGammaSpecFromXFrom%sPt",nameSecondaries[k].Data()));
        }
    }
    if (isCalo && !isPCM) {
        for (Int_t k = 0; k < 4; k++){
            histoSecondaryGammaFromXSpecPt[k]                       = (TH1D*)histoESDCaloGammaPt->Clone(        Form("SecondaryGammaSpecFromXFrom%sPt",nameSecondaries[k].Data()));
            histoSecondaryGammaFromXSpecPtOrBin[k]                  = (TH1D*)histoESDCaloGammaPt_OrBin->Clone(  Form("SecondaryGammaSpecFromXFrom%sPt",nameSecondaries[k].Data()));

            // overwrite spectra loaded from pure MC if running data
            if (!isRunMC) histoGammaTrueSecCaloGammaFromX_Pt[k]     = (TH1D*)histoESDCaloGammaPt->Clone(        Form("SecondaryMCGammaSpecFromXFrom%sPt",nameSecondaries[k].Data()));
        }
    }
    for (Int_t k = 0; k < 4; k++){
        // secondary spectra
        if (histoFracAllGammaToSecFromX_Pt[k]) {
            histoSecondaryGammaFromXSpecPt[k]->Multiply(histoFracAllGammaToSecFromX_Pt[k]);
            histoSecondaryGammaFromXSpecPt[k]->Scale(1./nEvt);
            histoSecondaryGammaFromXSpecPt[k]->Scale(scaleFactorsSec[k]);
        } else {
            histoSecondaryGammaFromXSpecPt[k]                       = NULL;
        }

        // scale reconstructed gamma with secondaries
        if (!isRunMC){
            if(isPCM){
                if (histoFracAllGammaToSecFromX_Pt[k]) {
                    histoGammaTrueSecConvGammaFromX_Pt[k]->Multiply(histoFracAllGammaToSecFromX_Pt[k]);
                    histoGammaTrueSecConvGammaFromX_Pt[k]->Scale(1./nEvt);
                } else histoGammaTrueSecConvGammaFromX_Pt[k] = NULL;
            }else if (isCalo && !isPCM){
                if (histoFracAllGammaToSecFromX_Pt[k]) {
                    histoGammaTrueSecCaloGammaFromX_Pt[k]->Multiply(histoFracAllGammaToSecFromX_Pt[k]);
                    histoGammaTrueSecCaloGammaFromX_Pt[k]->Scale(1./nEvt);
                } else histoGammaTrueSecCaloGammaFromX_Pt[k] = NULL;
            }
        }else{
            if (isPCM ) {
                histoGammaTrueSecConvGammaFromX_Pt[k]->Scale(1./nEvtMC);
            } else if (isCalo && !isPCM) {
                histoGammaTrueSecCaloGammaFromX_Pt[k]->Scale(1./nEvtMC);
            }
        }

        // secondary spectra original binning
        if (histoFracAllGammaToSecFromX_Pt_OrBin[k]) {
            histoSecondaryGammaFromXSpecPtOrBin[k]->Multiply(histoFracAllGammaToSecFromX_Pt_OrBin[k]);
            histoSecondaryGammaFromXSpecPtOrBin[k]->Scale(scaleFactorsSec[k]);
        } else {
            histoSecondaryGammaFromXSpecPtOrBin[k]                  = NULL;
        }

        // secondary spectra pileup corrected
        if (histoFracAllGammaToSecFromX_PileUp_Pt[k] && doPileUpCorr) {
            histoSecondaryGammaFromXSpecPileUpPt[k]->Multiply(histoFracAllGammaToSecFromX_PileUp_Pt[k]);
	    histoSecondaryGammaFromXSpecPileUpPt[k]->Scale(1./nEvt);
            histoSecondaryGammaFromXSpecPileUpPt[k]->Scale(scaleFactorsSec[k]);
        } else {
            histoSecondaryGammaFromXSpecPileUpPt[k]                 = NULL;
        }
    }

    //******************************************************************************************
    //******** cocktail: use pileup corrected "rest" contribution for PCM in doPileUp case *****
    //******************************************************************************************
    if (hasCocktailInput && isPCM && doPileUpCorr) {
        if(histoGammaTrueSecCocktailGammaFromX_Pt[3])
            delete histoGammaTrueSecCocktailGammaFromX_Pt[3];
        if(histoGammaTrueSecCocktailGammaFromX_PtOrBin[3])
            delete histoGammaTrueSecCocktailGammaFromX_PtOrBin[3];

        histoGammaTrueSecCocktailGammaFromX_Pt[3]          = (TH1D*)histoESDConvGammaPtPileUp->Clone("TrueSecondaryConvGammaCocktailFromXFromRest_Pt");
        if(histoGammaTrueSecCocktailGammaFromX_Pt[3])
            histoGammaTrueSecCocktailGammaFromX_Pt[3]->Multiply(histoFracAllGammaToSecFromX_PileUp_Pt[3]);

        histoGammaTrueSecCocktailGammaFromX_PtOrBin[3]     = (TH1D*)histoESDConvGammaPt_OrBin->Clone("TrueSecondaryConvGammaCocktailFromXFromRest_Pt_OriginalBinning");
        if(histoGammaTrueSecCocktailGammaFromX_PtOrBin[3]) {
            histoGammaTrueSecCocktailGammaFromX_PtOrBin[3]->Multiply(histoPileUpCorrectionFactor_Pt_OrBin);
            histoGammaTrueSecCocktailGammaFromX_PtOrBin[3]->Multiply(histoFracAllGammaToSecFromX_Pt_OrBin[3]); // histoFracAllGammaToSecFromX_PileUp_Pt_OrBin should be used here, could be constructed in the future (small difference)
	}
    }

    //******************************************************************************************
    //******************* Scale cocktail secondary spectra *************************************
    //******************************************************************************************
    if (hasCocktailInput) {
        for (Int_t k = 0; k<3; k++){
            histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]->Scale(1./nEvt);
        }
        histoGammaTrueSecCocktailGammaFromX_Pt[3]->Scale(1./nEvt);
    }

    //******************************************************************************************
    //******************** Pileup correction factor plot ***************************************
    //******************************************************************************************
    if( doPileUpCorr && isPCM) {
        // rebin pileup correction factor used for unfolding correction strand
        TH1D* histoPileUpCorrectionFactor_PtTemp                    = (TH1D*)histoPileUpCorrectionFactor_Pt_OrBin->Clone("histoPileUpCorrectionFactor_PtTemp");
        histoPileUpCorrectionFactor_PtTemp                          = RebinTH1D(histoPileUpCorrectionFactor_PtTemp,histoPileUpCorrectionFactor_Pt,kFALSE);
        histoPileUpCorrectionFactor_PtTemp->SetBinContent( 1, 0);
        histoPileUpCorrectionFactor_PtTemp->SetBinError(   1, 0);

        // pileup correction factor plot (rebinned version of factor used in unfolding correction strand)
        Bool_t includeHistoPileUpCorrectionFactorMC_Pt              = kTRUE;
        TCanvas* canvasPileUpCorrFactor                             = GetAndSetCanvas("canvasPileUpCorrFactor");


	if( energy.Contains("13TeV") )  {
	  SetHistogramm(histoPileUpCorrectionFactor_PtTemp,"#it{p}_{T} (GeV/#it{c})","Correction Factor (%)",0.64,1.02);
	} else {
	  SetHistogramm(histoPileUpCorrectionFactor_PtTemp,"#it{p}_{T} (GeV/#it{c})","Correction Factor (%)",0.84,1.02);
	}

        DrawGammaSetMarker(histoPileUpCorrectionFactor_PtTemp, 24, 2.0, kBlack, kBlack);

        histoPileUpCorrectionFactor_PtTemp->Draw("e1");
        DrawGammaLines(minPtGamma, maxPtGamma,1.0, 1.0, 1, kGray+2, 2);
        histoPileUpCorrectionFactor_PtTemp->Draw("e1,same");

        if (includeHistoPileUpCorrectionFactorMC_Pt && !isRunMC) {
            DrawGammaSetMarker(histoPileUpCorrectionFactorMC_Pt, 24, 2.0, kBlue, kBlue);
            histoPileUpCorrectionFactorMC_Pt->Draw("e1,same");

            TLegend* legendPileUp1                                  = GetAndSetLegend(0.15,0.15,2,1);
            legendPileUp1->AddEntry(histoPileUpCorrectionFactor_PtTemp,"data","pl");
            legendPileUp1->AddEntry(histoPileUpCorrectionFactorMC_Pt,"MC","pl");
            legendPileUp1->Draw("same");
        }

        PutProcessLabelAndEnergyOnPlot( 0.935, 0.25, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,31);

        canvasPileUpCorrFactor->SaveAs(Form("%s/%s_%s_PileUpCorrFactor_%s.%s",outputDir.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasPileUpCorrFactor;

        // inverse of pileup correction factor plot with fit
        TCanvas* canvasInversePileUpCorrFactor                      = GetAndSetCanvas("canvasInversePileUpCorrFactor");

        SetHistogramm(histoRatioWithWithoutPileUp,"#it{p}_{T} (GeV/#it{c})","1 / Correction Factor (%)",0.9,histoRatioWithWithoutPileUp->GetMaximum()*1.1);
        DrawGammaSetMarker(histoRatioWithWithoutPileUp, 24, 2.0, kBlack, kBlack);

        histoRatioWithWithoutPileUpFit->SetLineColor(kBlue-2);
        histoRatioWithWithoutPileUpFit->SetLineWidth(2);

        histoRatioWithWithoutPileUp->Draw("e1");
        DrawGammaLines(minPtGamma, maxPtGamma,1.0, 1.0, 1, kGray+2, 2);
        histoRatioWithWithoutPileUp->Draw("e1,same");
        histoRatioWithWithoutPileUpFit->Draw("same");

        TLegend* legendPileUp2                                      = GetAndSetLegend(0.15,0.15,3,1);
        legendPileUp2->AddEntry(histoRatioWithWithoutPileUp,    "#gamma_{raw} / #gamma_{raw, pileup sub.}","pl");
        legendPileUp2->AddEntry(histoRatioWithWithoutPileUpFit, Form("fit, #chi^{2}/ndf = %.2f", histoRatioWithWithoutPileUpFit->GetChisquare() / histoRatioWithWithoutPileUpFit->GetNDF()),"pl");
        legendPileUp2->Draw("same");

        PutProcessLabelAndEnergyOnPlot( 0.935, 0.25, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,31);

        canvasInversePileUpCorrFactor->SaveAs(Form("%s/%s_%s_PileUpCorrFactorFit_%s.%s",outputDir.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasInversePileUpCorrFactor;

        // pileup correction factor comparison plots
        TCanvas* canvasPileUpCorrFactorComp                         = GetAndSetCanvas("canvasPileUpCorrFactorComp");

        SetHistogramm(histoPileUpCorrectionFactor_PtTemp,"#it{p}_{T} (GeV/#it{c})","Correction Factor (%)",0.82,1.02);
        DrawGammaSetMarker(histoPileUpCorrectionFactor_PtTemp,      24, 2.0, kBlack,    kBlack);
        DrawGammaSetMarker(histoPileUpCorrectionFactor_Pt,          25, 2.0, kBlue-2,   kBlue-2);
        DrawGammaSetMarker(histoPileUpCorrectionFactorNoFit_Pt,     25, 2.5, kAzure+2,  kAzure+2);
        DrawGammaSetMarker(histoPileUpCorrectionFactor_Pt_OrBin,    27, 2.5, kRed+2,    kRed+2);
        DrawGammaSetMarker(histoPileUpCorrectionFactorMC_Pt,        28, 2.5, kOrange+2, kOrange+2);

        histoPileUpCorrectionFactor_PtTemp->Draw("e1");
        DrawGammaLines(minPtGamma, maxPtGamma,1.0, 1.0, 1, kGray+2, 2);
        histoPileUpCorrectionFactor_PtTemp->Draw("e1,same");
        histoPileUpCorrectionFactor_Pt->Draw("e1,same");
        histoPileUpCorrectionFactorNoFit_Pt->Draw("e1,same");
        histoPileUpCorrectionFactor_Pt_OrBin->Draw("e1,same");
        histoPileUpCorrectionFactorMC_Pt->Draw("e1,same");

        TLegend* legendPileUp3                                      = GetAndSetLegend(0.2,0.15,5,1);
        legendPileUp3->AddEntry(histoPileUpCorrectionFactorNoFit_Pt,    "from raw gamma ratio",         "pl");
        legendPileUp3->AddEntry(histoPileUpCorrectionFactor_Pt_OrBin,   "from fit in or. bin.",         "pl");
        legendPileUp3->AddEntry(histoPileUpCorrectionFactor_PtTemp,     "rebin from fit in or. bin.",   "pl");
        legendPileUp3->AddEntry(histoPileUpCorrectionFactor_Pt,         "from fit in ana. bin.",        "pl");
        legendPileUp3->AddEntry(histoPileUpCorrectionFactorMC_Pt,       "MC from raw gamma ratio",      "pl");
        legendPileUp3->Draw("same");

        PutProcessLabelAndEnergyOnPlot( 0.935, 0.25, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,31);

        canvasPileUpCorrFactorComp->SaveAs(Form("%s/%s_%s_PileUpCorrFactorComp_%s.%s",outputDir.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasPileUpCorrFactorComp;
    }

    //**********************************************************************************
    //******************** Background Plot *********************************************
    //**********************************************************************************
    TCanvas* canvasBackground                                       = new TCanvas("canvasBackground","",200,10,1000*1.25,1100*1.25);  // gives the page size
    DrawGammaCanvasSettings( canvasBackground, 0.1, 0.02, 0.02, 0.09);
    canvasBackground->SetLogy();

        TH1D* histoGammaRawSpectrum_PtMC                            = NULL;
        TH1D* histoGammaRawSpectrum_Pt                              = NULL;
        TLegend* legendBackground                                   = GetAndSetLegend(0.6,0.75,2,1);
        if (isPCM ) {
            histoGammaRawSpectrum_PtMC                              = (TH1D*)histoMCrecGamma_Pt->Clone("histoGammaRawSpectrum_PtMC");
            histoGammaRawSpectrum_PtMC->Scale(1./nEvtMC);
            histoGammaRawSpectrum_Pt                                = (TH1D*)histoESDConvGammaPt->Clone("histoGammaRawSpectrum_Pt");
            histoGammaRawSpectrum_Pt->Scale(1./nEvt);

            DrawGammaSetMarker(histoGammaRawSpectrum_Pt, 20, 1.0, 1, 1);
            DrawGammaSetMarker(histoMCrecBackground_Pt, 20, 1.0, 2, 2);

            SetHistogramm(histoMCrecBackground_Pt,"#it{p}_{T} (GeV/#it{c})",Form("Background for #gamma in |#eta| < %g",eta),-99,-99,1.0,1.7);
            SetHistogramm(histoGammaRawSpectrum_Pt,"#it{p}_{T} (GeV/#it{c})","Raw yield",histoGammaRawSpectrum_Pt->GetMinimum(0)/1000., histoGammaRawSpectrum_Pt->GetMaximum()*10.,1.0,1.2);

            histoGammaRawSpectrum_Pt->DrawCopy("");
            histoMCrecBackground_Pt->Draw("same");
        }
        if (isCalo && !isPCM) {
            histoGammaRawSpectrum_PtMC                              = (TH1D*) histoMCrecGammaCalo_Pt->Clone("histoGammaRawSpectrum_PtMC");
            histoGammaRawSpectrum_PtMC->Scale(1./nEvtMC);
            histoGammaRawSpectrum_Pt                                = (TH1D*) histoESDCaloGammaPt->Clone("histoGammaRawSpectrum_Pt");
            histoGammaRawSpectrum_Pt->Scale(1./nEvt);

            DrawGammaSetMarker(histoGammaRawSpectrum_Pt, 20, 1.0, 1, 1);
            DrawGammaSetMarker(histoMCrecBackground_Pt, 20, 1.0, 2, 2);

            SetHistogramm(histoMCrecBackground_Pt,"#it{p}_{T} (GeV/#it{c})",Form("Background for #gamma in |#eta| < %g",etaCalo),-99,-99,1.0,1.7);
            SetHistogramm(histoGammaRawSpectrum_Pt,"#it{p}_{T} (GeV/#it{c})","Raw yield",histoGammaRawSpectrum_Pt->GetMinimum(0)/300., histoGammaRawSpectrum_Pt->GetMaximum()*10., 1.0, 1.2);

            histoGammaRawSpectrum_Pt->DrawCopy("");
            histoMCrecBackground_Pt->Draw("same");
        }

        legendBackground->AddEntry(histoGammaRawSpectrum_Pt,"raw #gamma spectrum","pl");
        legendBackground->AddEntry(histoMCrecBackground_Pt,"#gamma background","pl");
        legendBackground->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.935, 0.95, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,31);

    canvasBackground->SaveAs(Form("%s/%s_Background_%s_%s.%s",outputDir.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasBackground;

    //**********************************************************************************
    //******************** Secondary Spectra Plot **************************************
    //**********************************************************************************
    TCanvas* canvasSecSpec                                          = new TCanvas("canvasSecSpec","",200,10,1000*1.25,1100*1.25);  // gives the page size
    DrawGammaCanvasSettings( canvasSecSpec, 0.15, 0.02, 0.02, 0.09);
    canvasSecSpec->SetLogy();

        TLegend* legendSecSpec                                      = GetAndSetLegend2(0.35, 0.935-0.035*1.1*4, 0.95, 0.935, 0.035, 3, "", 42, 0.1);

        // secondary gamma spectra plotting as would they would be subtracted using MC fractions
        for (Int_t k = 0; k < 4; k++){
            legendSecSpec->AddEntry((TObject*)0, Form("sec #gamma from %s:", nameLabelSecondaries[k].Data()),"");
            // plotting MC unscaled secondaries
            if (isPCM ) {
                if (histoGammaTrueSecConvGammaFromX_Pt[k]){
                    SetHistogramm(histoGammaTrueSecConvGammaFromX_Pt[k],"#it{p}_{T} (GeV/#it{c})","#frac{1}{#it{N}_{ev.}} #frac{d#it{N}}{d#it{p}_{T}} (#it{c}/GeV)",histoGammaTrueSecConvGammaFromX_Pt[k]->GetMinimum(0)/1000.,histoGammaTrueSecConvGammaFromX_Pt[k]->GetMaximum()*10.,1.0,1.6);
                    DrawGammaSetMarker(histoGammaTrueSecConvGammaFromX_Pt[k], markerStyleSecWithToy[k], markerSizeSec[k]*1.5, colorSecFromToy[k], colorSecFromToy[k]);
                    histoGammaTrueSecConvGammaFromX_Pt[k]->Draw("same");
                    legendSecSpec->AddEntry(histoGammaTrueSecConvGammaFromX_Pt[k], "MC     ","pl");
                } else {
                    legendSecSpec->AddEntry((TObject*)0, "","");
                }
            } else if (isCalo && !isPCM) {
                if (histoGammaTrueSecCaloGammaFromX_Pt[k]){
                    SetHistogramm(histoGammaTrueSecCaloGammaFromX_Pt[k],"#it{p}_{T} (GeV/#it{c})","#frac{1}{#it{N}_{ev.}} #frac{d#it{N}}{d#it{p}_{T}} (#it{c}/GeV)",histoGammaTrueSecCaloGammaFromX_Pt[k]->GetMinimum(0)/1000.,histoGammaTrueSecCaloGammaFromX_Pt[k]->GetMaximum()*100.,1.0,1.6);
                    DrawGammaSetMarker(histoGammaTrueSecCaloGammaFromX_Pt[k], markerStyleSecWithToy[k], markerSizeSec[k]*1.5, colorSecFromToy[k], colorSecFromToy[k]);
                    histoGammaTrueSecCaloGammaFromX_Pt[k]->Draw("same");
                    legendSecSpec->AddEntry(histoGammaTrueSecCaloGammaFromX_Pt[k], "MC    ","pl");
                } else {
                    legendSecSpec->AddEntry((TObject*)0, "","");
                }
            }

            // plotting cocktail on top of the MC unscaled secondaries
            if (hasCocktailInput) {
                if (histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]){
                    if (isPCM ) {
                        SetHistogramm(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k],"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
                    } else if (isCalo && !isPCM) {
                        SetHistogramm(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k],"#it{p}_{T} (GeV/#it{c})","Secondary #gamma");
                    }

                    DrawGammaSetMarker(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k], markerStyleSec[k] , markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                    histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]->DrawCopy("same");
                    legendSecSpec->AddEntry(histoGammaSecGammaFromX_Cocktail_Raw_Pt[k],"cocktail","pl");
                } else {
                    legendSecSpec->AddEntry((TObject*)0, "","");
                }
            } else {
                if (histoSecondaryGammaFromXSpecPt[k]){
                    if (isPCM ) {
                        SetHistogramm(histoSecondaryGammaFromXSpecPt[k], "#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
                    } else if (isCalo && !isPCM) {
                        SetHistogramm(histoSecondaryGammaFromXSpecPt[k], "#it{p}_{T} (GeV/#it{c})","Secondary #gamma");
                    }
                    DrawGammaSetMarker(histoSecondaryGammaFromXSpecPt[k], markerStyleSec[k], markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                    histoSecondaryGammaFromXSpecPt[k]->Draw("same");
                    legendSecSpec->AddEntry(histoSecondaryGammaFromXSpecPt[k],"data scaled","pl");
                } else {
                    legendSecSpec->AddEntry((TObject*)0, "","");
                }
            }
        }
        legendSecSpec->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.95, 0.935-0.035*1.1*4, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,31);

    canvasSecSpec->SaveAs(Form("%s/%s_%s_SecondarySpectra_%s.%s", outputDir.Data(), textPi0New.Data(), nameRec.Data(), cutSelection.Data(), suffix.Data()));
    delete canvasSecSpec;

    //**********************************************************************************
    //******************** Secondary Fractions Plot ************************************
    //**********************************************************************************
    TCanvas* canvasSecFrac = new TCanvas("canvasSecFrac","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasSecFrac, 0.072, 0.02, 0.035, 0.09);

        Double_t plotYrange[2]              = {-99.,-99.};
        if(energy.Contains("PbPb_2.76TeV")){
            plotYrange[0] = 0.;
            plotYrange[1] = 35.e-3;
        }
        if(energy.BeginsWith("8TeV")){
            if(mode==2){
                plotYrange[0] = 0.;
                plotYrange[1] = 35.e-3;
            }else if(mode==4){
                plotYrange[0] = 0.;
                plotYrange[1] = 80.e-3;
            }
        }
        Int_t nColumnsSec                   = 1;
        Double_t startXSecLegend            = 0.65;
        Int_t nRowsSec                      = 4;
        Double_t marginSecLegend            = 0.15;
        if (hasCocktailInput){
            nColumnsSec                     = 3;
            startXSecLegend                 = 0.55;
            marginSecLegend                 = 0.15 ;
        }
        TLegend* legendSecFrac              = GetAndSetLegend2(startXSecLegend, 0.935-0.035*1.1*nRowsSec, 0.95, 0.935, 0.035, nColumnsSec, "", 42, marginSecLegend);
        if (!hasCocktailInput) {
            for (Int_t k = 0; k < 4; k++){
                if (histoFracAllGammaToSecFromX_Pt[k]){
                    SetHistogramm(histoFracAllGammaToSecFromX_Pt[k],"#it{p}_{T} (GeV/#it{c})","C_{sec, #gamma}", plotYrange[0] ,plotYrange[1], 1.0, 0.9);
                    DrawGammaSetMarker(histoFracAllGammaToSecFromX_Pt[k], markerStyleSec[k], markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                    histoFracAllGammaToSecFromX_Pt[k]->Draw("same");
                    legendSecFrac->AddEntry(histoFracAllGammaToSecFromX_Pt[k],Form("#gamma from %s", nameLabelSecondaries[k].Data()),"pl");
                }
            }
        } else {
            for (Int_t k = 0; k < 4; k++){
                if (histoFracAllGammaToSecFromX_Cocktail_Pt[k]){
                    SetHistogramm(histoFracAllGammaToSecFromX_Cocktail_Pt[k],"#it{p}_{T} (GeV/#it{c})","C_{sec, #gamma}", plotYrange[0] ,plotYrange[1], 1.0, 0.9);
                    SetHistogramm(histoFracAllGammaToSecFromX_Pt[k],"#it{p}_{T} (GeV/#it{c})","C_{sec, #gamma}", plotYrange[0] ,plotYrange[1], 1.0, 0.9);
                    DrawGammaSetMarker(histoFracAllGammaToSecFromX_Cocktail_Pt[k], markerStyleSec[k], markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                    DrawGammaSetMarker(histoFracAllGammaToSecFromX_Pt[k], markerStyleSecWithToy[k], markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                    if (k == 0){
                        histoFracAllGammaToSecFromX_Pt[k]->Draw();
                        histoFracAllGammaToSecFromX_Cocktail_Pt[k]->Draw("same");
                    } else {
                        histoFracAllGammaToSecFromX_Pt[k]->Draw("same");
                        histoFracAllGammaToSecFromX_Cocktail_Pt[k]->Draw("same");
                    }
                    legendSecFrac->AddEntry((TObject*)0,Form("#gamma from %s", nameLabelSecondaries[k].Data()),"");
                    legendSecFrac->AddEntry(histoFracAllGammaToSecFromX_Cocktail_Pt[k],"cocktail","pl");
                    legendSecFrac->AddEntry(histoFracAllGammaToSecFromX_Pt[k],"MC","pl");
                }
            }
        }
        legendSecFrac->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.94, 0.94-0.035*1.1*nRowsSec, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,31);

    canvasSecFrac->SaveAs(Form("%s/%s_%s_SecondaryGammaFraction_%s.%s", outputDir.Data(), textPi0New.Data(), nameRec.Data(), cutSelection.Data(), suffix.Data()));
    delete canvasSecFrac;

    //**********************************************************************************
    //******************** Purity Plot *************************************************
    //**********************************************************************************
    TCanvas *canvasPurity                       = GetAndSetCanvas("canvasPurity");

        TLegend* legendPurity                   = NULL;

        if ( (isPCM ) || (isCalo && !isPCM) ){
            DrawGammaSetMarker(histoGammaPurity_Pt, 24, 1.0, kRed+2, kRed+2);
            DrawGammaSetMarker(histoGammaTruePurity_Pt, 20, 1.0, 1, 1);
            if(doPileUpCorr)DrawGammaSetMarker(histoGammaTruePurity_PileUp_Pt, 20, 1.0, kBlue+2, kBlue+2);

            if (isPCM )             SetHistogramm(histoGammaTruePurity_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{pur,#gamma} in |#eta| < %g",eta),0.6, 1.);
            if (isCalo && !isPCM)   SetHistogramm(histoGammaTruePurity_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{pur,#gamma} in |#eta| < %g",etaCalo),0.6, 1.);
            histoGammaTruePurity_Pt->Draw();
            histoGammaPurity_Pt->Draw("same");
            if (doPileUpCorr) histoGammaTruePurity_PileUp_Pt->Draw("same");

            if (doPileUpCorr)   legendPurity    = GetAndSetLegend(0.16,0.15,3);
            else                legendPurity    = GetAndSetLegend(0.16,0.15,2);
            legendPurity->AddEntry(histoGammaPurity_Pt,"purity");
            legendPurity->AddEntry(histoGammaTruePurity_Pt,"rescaled purity, secondaries removed");
            if (doPileUpCorr) legendPurity->AddEntry(histoGammaTruePurity_PileUp_Pt,"rescaled purity, secondaries removed, pileup");
            legendPurity->Draw();

            PutProcessLabelAndEnergyOnPlot( 0.18, 0.32, 0.035, cent, detectionProcess, "", 42, 0.03);

            canvasPurity->SaveAs(Form("%s/%s_Purity_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        }
        if ( isPCM && isCalo ){
            DrawGammaSetMarker(histoGammaCaloPurity_Pt, 24, 1.0, kRed+2, kRed+2);
            DrawGammaSetMarker(histoGammaCaloTruePurity_Pt, 20, 1.0, 1, 1);

            SetHistogramm(histoGammaCaloTruePurity_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{pur,#gamma} in |#eta| < %g",etaCalo),0.5, 1.);
            histoGammaCaloTruePurity_Pt->Draw();
            histoGammaCaloPurity_Pt->Draw("same");

            legendPurity                        = GetAndSetLegend(0.45,0.15,2);
            legendPurity->AddEntry(histoGammaPurity_Pt,"purity");
            legendPurity->AddEntry(histoGammaTruePurity_Pt,"rescaled purity, secondaries removed");
            legendPurity->Draw();

            PutProcessLabelAndEnergyOnPlot( 0.95, 0.32, 0.035, cent, detectionProcess2, "", 42, 0.03,"", 1, 1.25,31);

            canvasPurity->SaveAs(Form("%s/%s_PurityCalo_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        }

    delete canvasPurity;

    TCanvas*    canvasPurity2                   = GetAndSetCanvas("canvasPurity2");
    Bool_t      doPurityPlotSimple              = kTRUE;
    if (doPurityPlotSimple && ((isPCM && !isCalo) || (isCalo && !isPCM))) {

        if (isPCM )             SetHistogramm(histoGammaTruePurity_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{pur,#gamma} in |#eta| < %g",eta),0.8, 1.1);
        if (isCalo && !isPCM)   SetHistogramm(histoGammaTruePurity_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{pur,#gamma} in |#eta| < %g",etaCalo),0.8, 1.1);

        histoGammaTruePurity_Pt->Draw();
        DrawGammaSetMarker(histoGammaTruePurity_Pt, 24, 1.5, 1, 1);
        DrawGammaLines(minPtGamma, maxPtGamma,1.0, 1.0, 1, kGray+2, 2);
        histoGammaTruePurity_Pt->Draw("same");

        PutProcessLabelAndEnergyOnPlot( 0.18, 0.85, 0.035, cent, detectionProcess, "", 42, 0.03);

        canvasPurity2->SaveAs(Form("%s/%s_TruePurity_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    }
    delete canvasPurity2;

    //**********************************************************************************
    //******************** Conversion Prob Plot ****************************************
    //**********************************************************************************
    if ( isPCM ){
        TCanvas *canvasConvProb         = GetAndSetCanvas("canvasConvProb");
            DrawGammaSetMarker(histoGammaConvProb_MCPt, 24, 2.0, 1, 1);
            SetHistogramm(histoGammaConvProb_MCPt, "#it{p}_{T} (GeV/#it{c})",Form("#it{P}_{conv} in |#eta| < %g",eta), 0.04, 0.10);
            histoGammaConvProb_MCPt->Draw();

            Double_t minFit_fConv = 2.5;
            if(energy.BeginsWith("8TeV") && mode == 2) minFit_fConv = 3.;
            TF1 *fConv                  = new TF1("line","[0]",minFit_fConv,25.);
            histoGammaConvProb_MCPt->Fit(fConv,"QRME0");
            Double_t parameterProb[1];
            fConv->GetParameters(parameterProb);

            DrawGammaLines(minPtGamma, maxPtGamma,parameterProb[0], parameterProb[0], 1, kGray+2, 2);
            histoGammaConvProb_MCPt->Draw("same");

            PutProcessLabelAndEnergyOnPlot( 0.95, 0.23, 0.035, cent, detectionProcess, "", 42, 0.03,"", 1, 1.25,31);

        canvasConvProb->SaveAs(Form("%s/%s_ConversionProb_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasConvProb;

        // secondary reco eff plot (only used for cocktail sec corr)
        if (hasCocktailInput) {
            TCanvas *canvasConvProbSec  = GetAndSetCanvas("canvasConvProbSec");
            canvasConvProbSec->SetTopMargin(0.035);
            TLegend* legendSecConvProb       = GetAndSetLegend2(0.25, 0.935-0.035*1.4*2, 0.65, 0.935, 0.035, 2, "", 42, 0.12);

            DrawGammaSetMarker(histoGammaConvProb_MCPt, 20, 2.0, 1, 1);
            SetHistogramm(histoGammaConvProb_MCPt, "#it{p}_{T} (GeV/#it{c})",Form("#it{P}_{conv} in |#eta| < %g",eta), 0.0, 0.15);
            histoGammaConvProb_MCPt->Draw();
            legendSecConvProb->AddEntry(histoGammaConvProb_MCPt,"Prim. #gamma","p");

            for (Int_t k = 0; k < 3; k++){
                SetHistogramm(histoGammaSecondaryFromXConvProb_MCPt[k],"#it{p}_{T} (GeV/#it{c})",Form("#it{P}_{conv} in |#eta| < %g",eta), 0., 1.0);
                DrawGammaSetMarker(histoGammaSecondaryFromXConvProb_MCPt[k], markerStyleSecWithToy[k], markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                histoGammaSecondaryFromXConvProb_MCPt[k]->Draw("same");
                legendSecConvProb->AddEntry(histoGammaSecondaryFromXConvProb_MCPt[k],Form("Sec. #gamma from %s",nameLabelSecondaries[k].Data()),"p");
            }

            PutProcessLabelAndEnergyOnPlot( 0.95, 0.95, 0.035, cent, detectionProcess,"", 42, 0.03,"", 1, 1.25,31);
            legendSecConvProb->Draw();

            canvasConvProbSec->SaveAs(Form("%s/%s_ConversionProbSecondaries_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
            delete canvasConvProbSec;
        }
    }

    //**********************************************************************************
    //******************** Reconstruction Eff Plot *************************************
    //**********************************************************************************
    TCanvas *canvasRecoEff = GetAndSetCanvas("canvasRecoEff");
        if ( isPCM ){
            DrawGammaSetMarker(histoGammaPrimaryRecoEff_MCPt, 20, 1.0, 1, 1);
            SetHistogramm(histoGammaPrimaryRecoEff_MCPt,"#it{p}_{T,MC} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);
            histoGammaPrimaryRecoEff_MCPt->Draw();

            PutProcessLabelAndEnergyOnPlot( 0.95, 0.95, 0.035, cent, detectionProcess,"", 42, 0.03, "", 1, 1.25,31);

            canvasRecoEff->SaveAs(Form("%s/%s_ReconstructionEff_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        }
        if ( isCalo && !isPCM ){
            DrawGammaSetMarker(histoGammaPrimaryRecoEff_MCPt, 20, 1.0, 1, 1);
            SetHistogramm(histoGammaPrimaryRecoEff_MCPt,"#it{p}_{T,MC} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
            histoGammaPrimaryRecoEff_MCPt->Draw();

            PutProcessLabelAndEnergyOnPlot( 0.95, 0.95, 0.035, cent, detectionProcess,"", 42, 0.03, "", 1, 1.25,31);

            canvasRecoEff->SaveAs(Form("%s/%s_ReconstructionEff_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        }
        if ( isPCM && isCalo ){
            DrawGammaSetMarker(histoGammaCaloPrimaryRecoEff_MCPt, 20, 1.0, 1, 1);
            SetHistogramm(histoGammaCaloPrimaryRecoEff_MCPt,"#it{p}_{T,MC} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
            histoGammaCaloPrimaryRecoEff_MCPt->Draw();

            PutProcessLabelAndEnergyOnPlot( 0.95, 0.95, 0.035, cent, detectionProcess,"", 42, 0.03, "", 1, 1.25,31);

            canvasRecoEff->SaveAs(Form("%s/%s_ReconstructionEffCalo_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));

        }
    delete canvasRecoEff;

    // secondary reco eff plot (only used for cocktail sec corr)
    if (hasCocktailInput) {
        TCanvas *canvasRecoEffSec       = GetAndSetCanvas("canvasRecoEffSec");
        TLegend* legendSecRecoEff       = GetAndSetLegend2(0.25, 0.935-0.035*1.4*2, 0.65, 0.935, 0.035, 2, "", 42, 0.12);

        DrawGammaSetMarker(histoGammaPrimaryRecoEff_Pt, 20, 1.0, 1, 1);
        if ( isPCM ){
            SetHistogramm(histoGammaPrimaryRecoEff_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);
        } else if ( isCalo && !isPCM ){
            SetHistogramm(histoGammaPrimaryRecoEff_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 2.0);
        } else {
            SetHistogramm(histoGammaPrimaryRecoEff_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
        }
        legendSecRecoEff->AddEntry(histoGammaPrimaryRecoEff_Pt,"Prim. #gamma","p");

        histoGammaPrimaryRecoEff_Pt->Draw();

        if ( isPCM || isCalo ){
            for (Int_t k = 0; k < 3; k++){
                if (isPCM ) SetHistogramm(histoGammaSecFromXRecoEff_RecPt[k],"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);
                if (isCalo && !isPCM) SetHistogramm(histoGammaSecFromXRecoEff_RecPt[k],"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 2.0);
                DrawGammaSetMarker(histoGammaSecFromXRecoEff_RecPt[k], markerStyleSecWithToy[k], markerSizeSec[k], colorSecFromToy[k], colorSecFromToy[k]);
                histoGammaSecFromXRecoEff_RecPt[k]->Draw("same");
                legendSecRecoEff->AddEntry(histoGammaSecFromXRecoEff_RecPt[k],Form("Sec. #gamma from %s",nameLabelSecondaries[k].Data()),"p");
            }

        }
        PutProcessLabelAndEnergyOnPlot( 0.95, 0.95, 0.035, cent, detectionProcess, "", 42, 0.03,"", 1, 1.25,31);
        legendSecRecoEff->Draw();

        canvasRecoEffSec->SaveAs(Form("%s/%s_ReconstructionEffSec_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasRecoEffSec;
    }

    //**********************************************************************************
    //******************** Response Matrix for detector resolution *********************
    //**********************************************************************************
    if (isPCM){
        SetStyleHistoTH2ForGraphs(  histoGammaTruePrimaryConv_recPt_MCPt, "Reconstructed #it{p}_{T} (GeV/#it{c})","MC #it{p}_{T} (GeV/#it{c})", 0.035, 0.04,
                                    0.035, 0.04, 0.9, 1.0, 510, 510);
        histoGammaTruePrimaryConv_recPt_MCPt->GetYaxis()->SetRangeUser(0,25);
        histoGammaTruePrimaryConv_recPt_MCPt->GetXaxis()->SetRangeUser(0,25);
    }
    if (isCalo){
        SetStyleHistoTH2ForGraphs(  histoGammaTruePrimaryCalo_recPt_MCPt, "Reconstructed #it{p}_{T} (GeV/#it{c})","MC #it{p}_{T} (GeV/#it{c})", 0.035, 0.04,
                                    0.035, 0.04, 0.9, 1.0, 510, 510);
        histoGammaTruePrimaryCalo_recPt_MCPt->GetYaxis()->SetRangeUser(0,25);
        histoGammaTruePrimaryCalo_recPt_MCPt->GetXaxis()->SetRangeUser(0,25);
    }

    //**********************************************************************************
    //******************** Response Matrix Plot ****************************************
    //**********************************************************************************
    TCanvas * canvasResponseMatrix = new TCanvas("canvasResponseMatrix","",480,440);  // gives the page size
    DrawGammaCanvasSettings( canvasResponseMatrix, 0.09, 0.11, 0.02, 0.085);
    canvasResponseMatrix->SetLogz(1);
    canvasResponseMatrix->cd();

        if (isPCM){
            cout << "Responsematrix conv entries: "<< histoGammaTruePrimaryConv_recPt_MCPt->GetMinimum(0) << "\t"<< histoGammaTruePrimaryConv_recPt_MCPt->GetMaximum() << endl;
            histoGammaTruePrimaryConv_recPt_MCPt->GetZaxis()->SetRangeUser(histoGammaTruePrimaryConv_recPt_MCPt->GetMinimum(0), 2*histoGammaTruePrimaryConv_recPt_MCPt->GetMaximum());
            histoGammaTruePrimaryConv_recPt_MCPt->Draw("colz");
            PutProcessLabelAndEnergyOnPlot( 0.15, 0.95, 0.035, cent, detectionProcess, "", 42, 0.03);
            canvasResponseMatrix->SaveAs(Form("%s/%s_ResponseMatrix_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        }

        if (isCalo){
            cout << "Responsematrix calo entries: "<< histoGammaTruePrimaryCalo_recPt_MCPt->GetMinimum(0) << "\t"<< histoGammaTruePrimaryCalo_recPt_MCPt->GetMaximum() << endl;
            histoGammaTruePrimaryCalo_recPt_MCPt->GetZaxis()->SetRangeUser(histoGammaTruePrimaryCalo_recPt_MCPt->GetMinimum(0), 2*histoGammaTruePrimaryCalo_recPt_MCPt->GetMaximum());
            histoGammaTruePrimaryCalo_recPt_MCPt->Draw("colz");
            if (isPCM)  PutProcessLabelAndEnergyOnPlot( 0.15, 0.95, 0.035, cent, detectionProcess2, "", 42, 0.03);
            else        PutProcessLabelAndEnergyOnPlot( 0.15, 0.95, 0.035, cent, detectionProcess,  "", 42, 0.03);
            canvasResponseMatrix->SaveAs(Form("%s/%s_ResponseMatrixCalo_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));

        }
    //delete canvasResponseMatrix;

    // response matrix for cocktail secondary corr
    if (hasCocktailInput) {
        TCanvas* canvasResponseMatrixSec = new TCanvas("canvasResponseMatrixSec","",480,440);  // gives the page size
        DrawGammaCanvasSettings( canvasResponseMatrixSec, 0.09, 0.105, 0.02, 0.085);
        canvasResponseMatrixSec->SetLogz(1);
        canvasResponseMatrixSec->cd();

        if (isPCM || isCalo){
            cout << "Secondary responsematrix entries: "<< histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[0]->GetMinimum(0) << "\t"<< histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[0]->GetMaximum() << endl;
            histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[0]->GetZaxis()->SetRangeUser(histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[0]->GetMinimum(0), 2*histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[0]->GetMaximum());
            SetStyleHistoTH2ForGraphs(  histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[0], "Reconstructed #it{p}_{T} (GeV/#it{c})","MC #it{p}_{T} (GeV/#it{c})", 0.035, 0.04,
                                        0.035, 0.04, 0.9, 1.0, 510, 510);
            histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[0]->GetYaxis()->SetRangeUser(0,25);
            histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[0]->GetXaxis()->SetRangeUser(0,25);
            histoGammaTrueSecondaryFromX_MCPt_recPt_OrBin[0]->Draw("colz");
            PutProcessLabelAndEnergyOnPlot( 0.15, 0.95, 0.035, cent, detectionProcess,"", 42, 0.03);
        }

        canvasResponseMatrixSec->SaveAs(Form("%s/%s_ResponseMatrixSecK0s_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasResponseMatrixSec;
    }

    //**********************************************************************************
    //************************ Unfolding of inclusive gamma spectrum *******************
    //**********************************************************************************
    TH1D* histoGammaCorrUnfoldReso_PtCopy                       = NULL;
    TH1D* histoGammaCorrUnfoldReso_Pt                           = NULL;
    TH1D* histoGammaCorrUnfoldReso_BinByBin_Pt                  = NULL;
    TH1D* histoGammaResolCorrUnfold_Pt                          = NULL;
    TH1D* histoGammaResolCorrUnfold_BinByBin_Pt                 = NULL;
    TH1D* histoGammaRaw_PileUpPuritySecondaryCorrected_Pt       = NULL;
    // do the same with MC rec gammas as a sanity check
    TH1D* histoMCrecGammaCorr_Pt                                = NULL;

    if (isPCM ){
        // correct inclusive photon spectrum for pileup contribution
        if (doPileUpCorr) histoESDConvGammaPt_OrBin->Multiply(histoPileUpCorrectionFactor_Pt_OrBin);

        // correct measured gamma spectra for secondaries and purity
        if(!hasCocktailInput){
            //subtract secondary contribution from inclusive spectrum
            CorrectGammaSecAndPurity(   histoESDConvGammaPt_OrBin,
                                        histoSecondaryGammaFromXSpecPtOrBin,
                                        histoGammaTruePurity_Pt_OrBin
                                    );
            CorrectGammaSecAndPurity(   histoMCrecGamma_Pt_OrBin,
                                        histoSecondaryGammaFromXSpecPtOrBin,
                                        histoGammaTruePurity_Pt_OrBin
                                    );
        } else {
            //subtract secondary contribution from cocktail from inclusive spectrum
            cout << "will use cocktail for secondary correction" << endl;
            CorrectGammaSecAndPurityCocktail(   histoESDConvGammaPt_OrBin,
                                                histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin,
                                                histoGammaTrueSecCocktailGammaFromX_PtOrBin[3],
                                                histoGammaTruePurity_Pt_OrBin
                                            );
            CorrectGammaSecAndPurityCocktail(   histoMCrecGamma_Pt_OrBin,
                                                histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin,
                                                histoGammaTrueSecCocktailGammaFromX_PtOrBin[3],
                                                histoGammaTruePurity_Pt_OrBin
                                            );
        }
        cout << "done with secondary corrections" << endl;
        // create histograms for unfolding for different techniques
        histoGammaCorrUnfoldReso_PtCopy                         = (TH1D*)histoESDConvGammaPt_OrBin->Clone("histoGammaCorrUnfoldReso_PtCopy");
        histoGammaCorrUnfoldReso_Pt                             = (TH1D*)histoESDConvGammaPt_OrBin->Clone("histoGammaCorrUnfoldReso_Pt");
        histoGammaCorrUnfoldReso_BinByBin_Pt                    = (TH1D*)histoESDConvGammaPt_OrBin->Clone("histoGammaCorrUnfoldReso_BinByBin_Pt");
        histoMCrecGammaCorr_Pt                                  = (TH1D*)histoMCrecGamma_Pt_OrBin->Clone("GammaSpecCorrESDMC");

        // gamma spectrum after secondary and purity correction
        histoGammaRaw_PileUpPuritySecondaryCorrected_Pt = (TH1D*)histoESDConvGammaPt_OrBin->Clone("GammaRaw_PileUpPuritySecondaryCorrected");
        histoGammaRaw_PileUpPuritySecondaryCorrected_Pt = RebinTH1D(histoGammaRaw_PileUpPuritySecondaryCorrected_Pt,histoESDConvGammaPt,kTRUE);
        histoGammaRaw_PileUpPuritySecondaryCorrected_Pt->Scale(1./nEvt);
        for(Int_t iBin=0; iBin<histoGammaRawSpectrum_Pt->GetNbinsX(); iBin++){
          if(histoGammaRawSpectrum_Pt->GetBinContent(iBin)==0){
            histoGammaRaw_PileUpPuritySecondaryCorrected_Pt->SetBinContent(iBin,0);
            histoGammaRaw_PileUpPuritySecondaryCorrected_Pt->SetBinError(iBin,0);
          }
        }
        // TH1D *histoGammaCorrUnfoldReso_SvD_Pt                = (TH1D*)histoESDConvGammaPt_OrBin->Clone("histoGammaCorrUnfoldReso_SvD_Pt");
        // TH1D *histoGammaCorrUnfoldReso_TUnfold_Pt            = (TH1D*)histoESDConvGammaPt_OrBin->Clone("histoGammaCorrUnfoldReso_TUnfold_Pt");

        // Unfolding setup:
        // RooUnfoldResponse constructor - create from already-filled histograms
        // "response" gives the response matrix, measured X truth.
        // "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
        // but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
        // in "truth" for unmeasured events (inefficiency).
        // "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
        // to indicate, respectively, no fakes and/or no inefficiency.
        RooUnfoldResponse   response(0,0,histoGammaTruePrimaryConv_recPt_MCPt);

        // unfold with different techniques, but same response matrix
        RooUnfoldBayes      unfold_Spectrum (           &response,histoGammaCorrUnfoldReso_Pt, nIterationsUnfolding);
        RooUnfoldBinByBin   unfold_SpectrumBinByBin (   &response,histoGammaCorrUnfoldReso_BinByBin_Pt);
        RooUnfoldBayes      unfold_SpectrumMCrec (      &response,histoMCrecGammaCorr_Pt, nIterationsUnfolding);
        //RooUnfoldSvd      unfold_SpectrumSvD (&response, histoGammaCorrUnfoldReso_SvD_Pt, 20);
        //RooUnfoldTUnfold  unfold_SpectrumTUnfold (&response,histoGammaCorrUnfoldReso_TUnfold_Pt);

        // get histograms from RooUnfold and rebin them
        histoGammaCorrUnfoldReso_Pt                             = (TH1D*)unfold_Spectrum.Hreco();
        histoGammaResolCorrUnfold_Pt                            = (TH1D*)histoESDConvGammaPt_OrBin->Clone("histoGammaResolCorrUnfold_Pt");
        histoGammaResolCorrUnfold_Pt->Divide(histoGammaCorrUnfoldReso_Pt);
        histoGammaCorrUnfoldReso_Pt                             = RebinTH1D(histoGammaCorrUnfoldReso_Pt,histoESDConvGammaPt,kTRUE);

        histoGammaCorrUnfoldReso_BinByBin_Pt                    = (TH1D*)unfold_SpectrumBinByBin.Hreco();
        histoGammaResolCorrUnfold_BinByBin_Pt                   = (TH1D*)histoESDConvGammaPt_OrBin->Clone("histoGammaResolCorrUnfold_BinByBin_Pt");
        histoGammaResolCorrUnfold_BinByBin_Pt->Divide(histoGammaCorrUnfoldReso_BinByBin_Pt);
        histoGammaCorrUnfoldReso_BinByBin_Pt                    = RebinTH1D(histoGammaCorrUnfoldReso_BinByBin_Pt,histoESDConvGammaPt,kTRUE);

        histoMCrecGammaCorr_Pt                                  = (TH1D*)unfold_SpectrumMCrec.Hreco();
        histoMCrecGammaCorr_Pt                                  = RebinTH1D(histoMCrecGammaCorr_Pt,histoESDConvGammaPt,kTRUE);

        // histoGammaCorrUnfoldReso_SvD_Pt                      = (TH1D*)unfold_SpectrumSvD.Hreco();
        // histoGammaCorrUnfoldReso_SvD_Pt                      = RebinTH1D(histoGammaCorrUnfoldReso_SvD_Pt,histoESDConvGammaPt,kTRUE);
        // histoGammaCorrUnfoldReso_TUnfold_Pt                  = (TH1D*)unfold_SpectrumTUnfold.Hreco();
        // histoGammaCorrUnfoldReso_TUnfold_Pt                  = RebinTH1D(histoGammaCorrUnfoldReso_TUnfold_Pt,histoESDConvGammaPt,kTRUE);

        // Correct inclusive photon spectrum with conversion probability & reco efficiency (both versus MC pt)
        CorrectGammaUnfoldResol( histoGammaCorrUnfoldReso_Pt,
                                 histoGammaConvProb_MCPt,
                                 histoGammaPrimaryRecoEff_MCPt,
                                 deltaEta, scaling, nEvt);
        CorrectGammaUnfoldResol( histoGammaCorrUnfoldReso_BinByBin_Pt,
                                 histoGammaConvProb_MCPt,
                                 histoGammaPrimaryRecoEff_MCPt,
                                 deltaEta, scaling, nEvt);
        CorrectGammaUnfoldResol( histoMCrecGammaCorr_Pt,
                                 histoGammaConvProb_MCPt,
                                 histoGammaPrimaryRecoEff_MCPt,
                                 deltaEta, scaling, nEvt);
        //CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_SvD_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
        //CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_TUnfold_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);

        SetHistogramm(histoGammaCorrUnfoldReso_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                      histoGammaCorrUnfoldReso_Pt->GetMinimum(0)/100., histoGammaCorrUnfoldReso_Pt->GetMaximum()*10, 1.0, 1.7);
        SetHistogramm(histoGammaCorrUnfoldReso_BinByBin_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", -99, -99, 1.0, 1.7);
        SetHistogramm(histoMCrecGammaCorr_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", -99, -99, 1.0, 1.7);

        DrawGammaSetMarker(histoGammaCorrUnfoldReso_Pt, 20, 1.0, 1, 1);
        DrawGammaSetMarker(histoGammaCorrUnfoldReso_BinByBin_Pt, 24, 1.0, kBlue, kBlue);
        DrawGammaSetMarker(histoMCrecGammaCorr_Pt, 20, 1.0, kGreen-1, kGreen-1);
        //SetHistogramm(histoGammaCorrUnfoldReso_SvD_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        //SetHistogramm(histoGammaCorrUnfoldReso_TUnfold_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
    }

    if (isCalo && !isPCM) {
        // correct measured gamma spectra for secondaries and purity
        if(!hasCocktailInput){
            //subtract secondary contribution from inclusive spectrum
            CorrectGammaSecAndPurity(   histoESDCaloGammaPt_OrBin,
                                        histoSecondaryGammaFromXSpecPtOrBin,
                                        histoGammaTruePurity_Pt_OrBin
                                    );
            CorrectGammaSecAndPurity(   histoMCrecGammaCalo_Pt_OrBin,
                                        histoSecondaryGammaFromXSpecPtOrBin,
                                        histoGammaTruePurity_Pt_OrBin
                                    );
        } else {
            //subtract secondary contribution from cocktail from inclusive spectrum
            cout << "will use cocktail for secondary correction" << endl;
            CorrectGammaSecAndPurityCocktail(histoESDCaloGammaPt_OrBin,
                                             histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin,
                                             histoGammaTrueSecCocktailGammaFromX_PtOrBin[3],
                                             histoGammaTruePurity_Pt_OrBin );
            CorrectGammaSecAndPurityCocktail(histoMCrecGammaCalo_Pt_OrBin,
                                             histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin,
                                             histoGammaTrueSecCocktailGammaFromX_PtOrBin[3],
                                             histoGammaTruePurity_Pt_OrBin );
        }

        // create histograms for unfolding for different techniques
        histoGammaCorrUnfoldReso_Pt                                 = (TH1D*)histoESDCaloGammaPt_OrBin->Clone("histoGammaCorrUnfoldReso_Pt");
        histoGammaCorrUnfoldReso_BinByBin_Pt                        = (TH1D*)histoESDCaloGammaPt_OrBin->Clone("histoGammaCorrUnfoldReso_BinByBin_Pt");
        histoMCrecGammaCorr_Pt                                      = (TH1D*)histoMCrecGammaCalo_Pt_OrBin->Clone("GammaSpecCorrESDMC");

        // gamma spectrum after secondary and purity correction
        histoGammaRaw_PileUpPuritySecondaryCorrected_Pt             = (TH1D*)histoESDCaloGammaPt_OrBin->Clone("GammaRaw_PileUpPuritySecondaryCorrected");
        histoGammaRaw_PileUpPuritySecondaryCorrected_Pt             = RebinTH1D(histoGammaRaw_PileUpPuritySecondaryCorrected_Pt,histoESDCaloGammaPt,kTRUE);
        histoGammaRaw_PileUpPuritySecondaryCorrected_Pt->Scale(1./nEvt);
        for(Int_t iBin=0; iBin<histoGammaRawSpectrum_Pt->GetNbinsX(); iBin++){
            if(histoGammaRawSpectrum_Pt->GetBinContent(iBin)==0){
                histoGammaRaw_PileUpPuritySecondaryCorrected_Pt->SetBinContent(iBin,0);
                histoGammaRaw_PileUpPuritySecondaryCorrected_Pt->SetBinError(iBin,0);
            }
        }

        // Unfolding setup:
        // RooUnfoldResponse constructor - create from already-filled histograms
        // "response" gives the response matrix, measured X truth.
        // "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
        // but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
        // in "truth" for unmeasured events (inefficiency).
        // "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
        // to indicate, respectively, no fakes and/or no inefficiency.
        RooUnfoldResponse   response(0,0,histoGammaTruePrimaryCalo_recPt_MCPt);

        // unfold with different techniques, but same response matrix
        RooUnfoldBayes      unfold_Spectrum (&response,histoGammaCorrUnfoldReso_Pt, nIterationsUnfolding);
        RooUnfoldBayes      unfold_SpectrumMCrec (&response,histoMCrecGammaCorr_Pt, nIterationsUnfolding);
        RooUnfoldBinByBin   unfold_SpectrumBinByBin (&response,histoGammaCorrUnfoldReso_BinByBin_Pt);

        // get histograms from RooUnfold and rebin them
        histoGammaCorrUnfoldReso_Pt                                 = (TH1D*)unfold_Spectrum.Hreco();
        histoGammaResolCorrUnfold_Pt                                = (TH1D*)histoESDCaloGammaPt_OrBin->Clone("histoGammaResolCorrUnfold_Pt");
        histoGammaResolCorrUnfold_Pt->Divide(histoGammaCorrUnfoldReso_Pt);
        histoGammaCorrUnfoldReso_Pt                                 = RebinTH1D(histoGammaCorrUnfoldReso_Pt,histoESDCaloGammaPt,kTRUE);

        histoGammaCorrUnfoldReso_BinByBin_Pt                        = (TH1D*)unfold_SpectrumBinByBin.Hreco();
        histoGammaResolCorrUnfold_BinByBin_Pt                       = (TH1D*)histoESDCaloGammaPt_OrBin->Clone("histoGammaResolCorrUnfold_BinByBin_Pt");
        histoGammaResolCorrUnfold_BinByBin_Pt->Divide(histoGammaCorrUnfoldReso_BinByBin_Pt);
        histoGammaCorrUnfoldReso_BinByBin_Pt                        = RebinTH1D(histoGammaCorrUnfoldReso_BinByBin_Pt,histoESDCaloGammaPt,kTRUE);

        histoMCrecGammaCorr_Pt                                      = (TH1D*)unfold_SpectrumMCrec.Hreco();
        histoMCrecGammaCorr_Pt                                      = RebinTH1D(histoMCrecGammaCorr_Pt,histoESDCaloGammaPt,kTRUE);

        // Correct inclusive photon spectrum with reco efficiency (versus MC pt)
        CorrectGammaUnfoldResol(    histoGammaCorrUnfoldReso_Pt,
                                    histoGammaPrimaryRecoEff_MCPt,
                                    deltaEtaCalo,
                                    scalingCalo,
                                    nEvt
                                );
        CorrectGammaUnfoldResol(    histoMCrecGammaCorr_Pt,
                                    histoGammaPrimaryRecoEff_MCPt,
                                    deltaEtaCalo,
                                    scalingCalo,
                                    nEvt
                                );
        CorrectGammaUnfoldResol(    histoGammaCorrUnfoldReso_BinByBin_Pt,
                                    histoGammaPrimaryRecoEff_MCPt,
                                    deltaEtaCalo,
                                    scalingCalo,
                                    nEvt
                                );

        SetHistogramm(histoGammaCorrUnfoldReso_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                      histoGammaCorrUnfoldReso_Pt->GetMinimum(0)/100., histoGammaCorrUnfoldReso_Pt->GetMaximum()*10, 1.0, 1.7);
        SetHistogramm(histoGammaCorrUnfoldReso_BinByBin_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", -99, -99, 1.0, 1.7);
        SetHistogramm(histoMCrecGammaCorr_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", -99, -99, 1.0, 1.7);

        DrawGammaSetMarker(histoGammaCorrUnfoldReso_Pt, 20, 1.0, 1, 1);
        DrawGammaSetMarker(histoGammaCorrUnfoldReso_BinByBin_Pt, 24, 1.0, kBlue, kBlue);
        DrawGammaSetMarker(histoMCrecGammaCorr_Pt, 20, 1.0, kGreen-1, kGreen-1);
    }

    //**********************************************************************************
    //******************** Rebin resolution correction from unfolding ******************
    //**********************************************************************************
    TH1D* histoGammaResolCorrUnfold_Pt_Rebin                        = (TH1D*)histoGammaResolCorrUnfold_Pt->Clone(Form("%s_Rebin", histoGammaResolCorrUnfold_Pt->GetName()));
    histoGammaResolCorrUnfold_Pt_Rebin->Sumw2();
    if (isPCM)
        histoGammaResolCorrUnfold_Pt_Rebin                          = RebinTH1D(histoGammaResolCorrUnfold_Pt_Rebin,histoESDConvGammaPt,kFALSE);
    else if (isCalo && !isPCM)
        histoGammaResolCorrUnfold_Pt_Rebin                          = RebinTH1D(histoGammaResolCorrUnfold_Pt_Rebin,histoESDCaloGammaPt,kFALSE);

    TH1D* histoGammaResolCorrUnfold_BinByBin_Pt_Rebin               = (TH1D*)histoGammaResolCorrUnfold_BinByBin_Pt->Clone(Form("%s_Rebin", histoGammaResolCorrUnfold_BinByBin_Pt->GetName()));
    histoGammaResolCorrUnfold_BinByBin_Pt_Rebin->Sumw2();
    if (isPCM)
        histoGammaResolCorrUnfold_BinByBin_Pt_Rebin                 = RebinTH1D(histoGammaResolCorrUnfold_BinByBin_Pt_Rebin,histoESDConvGammaPt,kFALSE);
    else if (isCalo && !isPCM)
        histoGammaResolCorrUnfold_BinByBin_Pt_Rebin                 = RebinTH1D(histoGammaResolCorrUnfold_BinByBin_Pt_Rebin,histoESDCaloGammaPt,kFALSE);

    // resolution from binwise corrections (no unfolding)
    TH1D* histoGammaResolCorrEff_Pt                                 = (TH1D*)histoGammaPrimaryRecoEff_Pt->Clone("histoGammaResolCorrEff_Pt");
    histoGammaResolCorrEff_Pt->Divide(histoGammaPrimaryRecoEff_MCPt);

    //**********************************************************************************
    //******************** Reconstruction Eff Comparison Plot **************************
    //**********************************************************************************
    TCanvas*    canvasCompRecoEff                                   = GetAndSetCanvas("canvasCompRecoEff");
    TPad*       padCompRecoEff                                      = new TPad("padCompRecoEff",        "", 0., 0.25, 1., 1.00, -1, -1, -2);
    TPad*       padCompRecoEffRatio                                 = new TPad("padCompRecoEffRatio",   "", 0., 0.00, 1., 0.25, -1, -1, -2);

    DrawGammaCanvasSettings(    canvasCompRecoEff,  0.10, 0.015, 0.02, 0.09);
    DrawGammaPadSettings(       padCompRecoEff,     0.10, 0.015, 0.02, 0.00);
    DrawGammaPadSettings(       padCompRecoEffRatio,0.10, 0.015, 0.00, 0.30);

    padCompRecoEff->Draw();
    padCompRecoEffRatio->Draw();

    // convert primary reco. eff. in MC pT to rec. pT from unfolding resolution for a more reasonable comparison
    TH1D* histoGammaPrimaryRecoEff_MCPt_ConvertToRecPt              = (TH1D*)histoGammaPrimaryRecoEff_MCPt->Clone(Form("%s_ConvertToRecPt", histoGammaPrimaryRecoEff_MCPt->GetName()));
    histoGammaPrimaryRecoEff_MCPt_ConvertToRecPt->Multiply(histoGammaResolCorrUnfold_Pt_Rebin);

    // in case of PCM or Calo method, or PCM part of hybrid PCM-Calo
    if ( (isPCM ) || (isCalo && !isPCM) ){

        // reconstruction efficiencies
        padCompRecoEff->cd();

        if (isPCM ) {
            SetHistogramm(histoGammaPrimaryRecoEff_MCPt,                        "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);
            SetHistogramm(histoGammaPrimaryRecoEff_MCPt_ConvertToRecPt,         "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);
            SetHistogramm(histoGammaPrimaryRecoEff_Pt,                          "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0.05, 1.0);
            if(doPileUpCorr)SetHistogramm(histoGammaPrimaryRecoEff_PileUp_Pt,   "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);
        }
        if (isCalo && !isPCM) {
            SetHistogramm(histoGammaPrimaryRecoEff_MCPt,                        "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
            SetHistogramm(histoGammaPrimaryRecoEff_MCPt_ConvertToRecPt,         "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
            SetHistogramm(histoGammaPrimaryRecoEff_Pt,                          "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0.05, 1.0);
        }

        DrawGammaSetMarker(histoGammaPrimaryRecoEff_MCPt,                       20, 1.0, kBlack,    kBlack);
        DrawGammaSetMarker(histoGammaPrimaryRecoEff_MCPt_ConvertToRecPt,        34, 1.5, kBlue-2,   kBlue-2);
        DrawGammaSetMarker(histoGammaPrimaryRecoEff_Pt,                         24, 1.0, kOrange-3, kOrange-3);
        if(doPileUpCorr)DrawGammaSetMarker(histoGammaPrimaryRecoEff_PileUp_Pt,  25, 1.0, kRed+2,    kRed+2);

        histoGammaPrimaryRecoEff_Pt->Draw("e1");
        histoGammaPrimaryRecoEff_MCPt_ConvertToRecPt->Draw("e1,same");
        histoGammaPrimaryRecoEff_Pt->Draw("e1,same");
        histoGammaPrimaryRecoEff_MCPt->Draw("e1,same");
        if(doPileUpCorr)histoGammaPrimaryRecoEff_PileUp_Pt->Draw("e1,same");

        TLegend* legendCompRecoEff;
        if(doPileUpCorr)    legendCompRecoEff                       = GetAndSetLegend(0.15 ,0.05,4,1);
        else                legendCompRecoEff                       = GetAndSetLegend(0.15 ,0.05,3,1);
        legendCompRecoEff->SetTextSize(0.04);
        legendCompRecoEff->AddEntry(histoGammaPrimaryRecoEff_MCPt,                      "unfolding prim. reco. eff. in #it{p}_{T, MC}", "p");
        legendCompRecoEff->AddEntry(histoGammaPrimaryRecoEff_MCPt_ConvertToRecPt,       "unfolding prim. reco. eff. in #it{p}_{T}", "p");
        legendCompRecoEff->AddEntry(histoGammaPrimaryRecoEff_Pt,                        "binwise prim. reco. eff. in #it{p}_{T}", "p");
        if(doPileUpCorr)legendCompRecoEff->AddEntry(histoGammaPrimaryRecoEff_PileUp_Pt, "binwise prim. reco. eff. in #it{p}_{T}, pile-up", "p");
        legendCompRecoEff->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.935, 0.95, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,31);

        // ratios of reconstruction efficiencies to reco. eff. from unfolding in rec. pT
        padCompRecoEffRatio->cd();

        TH1D* histoRatioTemp[3]                                     = {NULL, NULL, NULL};
        histoRatioTemp[0]                                           = (TH1D*)histoGammaPrimaryRecoEff_MCPt->Clone("tempRatio0");
        histoRatioTemp[0]->Divide(histoGammaPrimaryRecoEff_MCPt_ConvertToRecPt);
        histoRatioTemp[1]                                           = (TH1D*)histoGammaPrimaryRecoEff_Pt->Clone("tempRatio1");
        histoRatioTemp[1]->Divide(histoGammaPrimaryRecoEff_MCPt_ConvertToRecPt);
        if (doPileUpCorr) {
            histoRatioTemp[2]                                       = (TH1D*)histoGammaPrimaryRecoEff_PileUp_Pt->Clone("tempRatio2");
            histoRatioTemp[2]->Divide(histoGammaPrimaryRecoEff_MCPt_ConvertToRecPt);
        }

        SetStyleHistoTH1ForGraphs(histoRatioTemp[0], "#it{p}_{T} (GeV/#it{c})","#frac{reco eff x}{prim. reco. eff. unfold.}" , 0.14, 0.15,
                                  0.12, 0.10,  0.85, 0.4, 510, 505);
        histoRatioTemp[0]->GetYaxis()->SetRangeUser(0.65,1.35);

        DrawGammaSetMarker(histoRatioTemp[0],                   20, 1.0, kBlack,    kBlack);
        DrawGammaSetMarker(histoRatioTemp[1],                   24, 1.0, kOrange-3, kOrange-3);
        if(doPileUpCorr)DrawGammaSetMarker(histoRatioTemp[2],   25, 1.0, kRed+2,    kRed+2);

        histoRatioTemp[0]->Draw("e1");

        DrawGammaLines(minPtGamma, maxPtGamma, 0.8, 0.8, 1,kGray+2,2);
        DrawGammaLines(minPtGamma, maxPtGamma, 1.0, 1.0, 1,kGray+2,2);
        DrawGammaLines(minPtGamma, maxPtGamma, 1.2, 1.2, 1,kGray+2,2);

        histoRatioTemp[0]->Draw("e1,same");
        histoRatioTemp[1]->Draw("e1,same");
        if(doPileUpCorr) histoRatioTemp[2]->Draw("e1,same");

        canvasCompRecoEff ->SaveAs(Form("%s/%s_CompRecEff_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    }

    // in case of hybrid PCM-Calo method
    if (isPCM && isCalo){
        padCompRecoEff->cd();
        SetHistogramm(histoGammaCaloPrimaryRecoEff_MCPt, "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
        SetHistogramm(histoGammaCaloPrimaryRecoEff_Pt, "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0.05, 1.0);

        DrawGammaSetMarker(histoGammaCaloPrimaryRecoEff_MCPt, 20, 1.0, kBlack, kBlack);
        DrawGammaSetMarker(histoGammaCaloPrimaryRecoEff_Pt, 24, 1.0, kBlue+1, kBlue+1);

        histoGammaCaloPrimaryRecoEff_Pt->Draw("e1");
        histoGammaCaloPrimaryRecoEff_MCPt->Draw("same,e1");

        TLegend* legendCompRecoEff = GetAndSetLegend(0.15 ,0.05,3,1);
        legendCompRecoEff->SetTextSize(0.04);
        legendCompRecoEff->AddEntry(histoGammaCaloPrimaryRecoEff_MCPt,"MC pT Reconstruction Efficiency (unfolding)");
        legendCompRecoEff->AddEntry(histoGammaCaloPrimaryRecoEff_Pt,"Reconstruction Efficiency for primary  #gamma");
        legendCompRecoEff->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.935, 0.95, 0.035, cent, detectionProcess2, "", 42, 0.03,"",1,1.25,31);
        canvasCompRecoEff ->SaveAs(Form("%s/%s_CompCaloRecEff_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    }

    TH1D* histoGammaCaloResolCorrEff_Pt = NULL;
    if (isPCM && isCalo){
        histoGammaCaloResolCorrEff_Pt = (TH1D*)histoGammaCaloPrimaryRecoEff_Pt->Clone("histoGammaCaloResolCorrEff_Pt");
        histoGammaCaloResolCorrEff_Pt->Divide(histoGammaCaloPrimaryRecoEff_MCPt);

        padCompRecoEffRatio->cd();
        SetStyleHistoTH1ForGraphs(  histoGammaCaloResolCorrEff_Pt, "#it{p}_{T} (GeV/#it{c})","#frac{Reco Eff x}{Primary MC Reco Eff}" , 0.14, 0.15,
                                    0.12, 0.10,  0.85, 0.4, 510, 505);
        histoGammaCaloResolCorrEff_Pt->GetYaxis()->SetRangeUser(0.7,1.25);

        histoGammaCaloResolCorrEff_Pt->DrawCopy("e1");

        DrawGammaLines(minPtGamma, maxPtGamma,1, 1,0.5,kGray+2,2);

        canvasCompRecoEff ->SaveAs(Form("%s/%s_CompCaloResolCorr_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    }

    delete padCompRecoEff;
    delete padCompRecoEffRatio;
    delete canvasCompRecoEff;

    //**********************************************************************************
    //******************** Corrected Photon Spectrum Plot ******************************
    //**********************************************************************************
    TCanvas *canvasCorrGammaSpecUnfold     = new TCanvas("canvasDecayGammaSpecMC","",200,10,1000*1.25,1100*1.25);  // gives the page size
    DrawGammaCanvasSettings( canvasCorrGammaSpecUnfold, 0.16, 0.02, 0.02, 0.09);
    canvasCorrGammaSpecUnfold->SetLogy();

        TH1D* Dummy                             = NULL;
        if (isPCM ) Dummy                       = (TH1D*)histoESDConvGammaPt->Clone("Dummy");
        if (isCalo && !isPCM) Dummy             = (TH1D*)histoESDCaloGammaPt->Clone("Dummy");

        if (isPCM )             SetHistogramm(Dummy,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", histoGammaCorrUnfoldReso_Pt->GetMinimum(0), histoGammaCorrUnfoldReso_Pt->GetMaximum(), 1.0, 1.7);
        if (isCalo && !isPCM)   SetHistogramm(Dummy,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", histoGammaCorrUnfoldReso_Pt->GetMinimum(0), histoGammaCorrUnfoldReso_Pt->GetMaximum(), 1.0, 1.7);
        Dummy->Draw();
        histoGammaCorrUnfoldReso_Pt->DrawCopy("same");
        histoGammaCorrUnfoldReso_BinByBin_Pt->DrawCopy("same");
        //histoGammaCorrUnfoldReso_SvD_Pt->DrawCopy("same");
        //histoGammaCorrUnfoldReso_TUnfold_Pt->DrawCopy("same");

        TLegend* legendGammaSpecUnfoldComp      = GetAndSetLegend2(0.6,0.945-2*1.1*0.035, 0.95,0.945, 0.035, 1, "", 42, 0.1);
        legendGammaSpecUnfoldComp->AddEntry(histoGammaCorrUnfoldReso_Pt,"Bayesian unfolding");
        legendGammaSpecUnfoldComp->AddEntry(histoGammaCorrUnfoldReso_BinByBin_Pt,"Bin-by-Bin unfolding");
        legendGammaSpecUnfoldComp->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.935, 0.945-2*1.1*0.035, 0.035, cent, detectionProcess,"", 42, 0.03, "", 1, 1.25,31);

    canvasCorrGammaSpecUnfold->SaveAs(Form("%s/%s_%s_CorrGammaUnfoldSpectrumPurityMinusSec_%s.%s",outputDir.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasCorrGammaSpecUnfold;

    //**********************************************************************************
    //******************** Resolution Correction Plot **********************************
    //**********************************************************************************
    TCanvas *canvasResolutionCorr               = GetAndSetCanvas("canvasResolutionCorr");
        TH1D*                   Dummy2          = NULL;
        if (isPCM )             Dummy2          = (TH1D*)histoESDConvGammaPt->Clone("Dummy2");
        if (isCalo && !isPCM)   Dummy2          = (TH1D*)histoESDCaloGammaPt->Clone("Dummy2");
        SetHistogramm(Dummy2,"#it{p}_{T} (GeV/#it{c})", "resolution correction",0,2);
        Dummy2->Draw();

        DrawGammaSetMarker(histoGammaResolCorrUnfold_Pt,            20, 1.0, kBlue+1, kBlue+1);
        DrawGammaSetMarker(histoGammaResolCorrUnfold_BinByBin_Pt,   21, 1.0, kGray+2, kGray+2);
        DrawGammaSetMarker(histoGammaResolCorrEff_Pt,               24, 1.5, kRed+2, kRed+2);

        DrawGammaLines(minPtGamma, maxPtGamma,1, 1,1, kGray);
        DrawGammaLines(minPtGamma, maxPtGamma,0.8, 0.8,1, kGray, 7);
        DrawGammaLines(minPtGamma, maxPtGamma,1.2, 1.2,1, kGray, 7);

        histoGammaResolCorrUnfold_BinByBin_Pt->DrawCopy("same");
        histoGammaResolCorrUnfold_Pt->DrawCopy("same");
        histoGammaResolCorrEff_Pt->DrawCopy("same");
        //histoGammaCorrUnfoldReso_SvD_Pt->DrawCopy("same");
        //histoGammaCorrUnfoldReso_TUnfold_Pt->DrawCopy("same");

        TLegend* legendResolutionCorr                = GetAndSetLegend(0.15,0.80,3);
        legendResolutionCorr->AddEntry(histoGammaResolCorrEff_Pt,"from effi");
        legendResolutionCorr->AddEntry(histoGammaResolCorrUnfold_Pt,"from Bayesian unfolding");
        legendResolutionCorr->AddEntry(histoGammaResolCorrUnfold_BinByBin_Pt,"from bin-by-bin unfolding");
        legendResolutionCorr->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.18, 0.25, 0.035, cent, detectionProcess,"", 42, 0.03);

    canvasResolutionCorr->SaveAs(Form("%s/%s_%s_ResolutionCorrection_%s.%s",outputDir.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasResolutionCorr;

    //******************************************************************************************
    // copy raw inc gamma spectrum and correct for;
    // - secondary contamination,
    // - purity,
    // - reco effi (vs rec pt with sec correction)
    // - conversion probability &
    //******************************************************************************************
    TH1D* histoGammaCorrEffiReso_Pt                 = NULL;
    // correct alternate way for conversion reco
    if (isPCM ) {
        histoGammaCorrEffiReso_Pt                   = (TH1D*)histoESDConvGammaPt->Clone("CorrGammaSpecPurityMinusSec");
        if (!hasCocktailInput)
            CorrectGammaEffiResol(  histoGammaCorrEffiReso_Pt,
                                    histoSecondaryGammaFromXSpecPt,
                                    histoGammaTruePurity_Pt,
                                    histoGammaConvProb_MCPt,
                                    histoGammaPrimaryRecoEff_Pt,
                                    deltaEta, scaling, nEvt
                                    );
        else
            CorrectGammaEffiResolCocktail(  histoGammaCorrEffiReso_Pt,
                                            histoGammaSecGammaFromX_Cocktail_Raw_Pt,    // already scaled per event
                                            histoGammaTrueSecCocktailGammaFromX_Pt[3],  // already scaled per event
                                            histoGammaTruePurity_Pt,
                                            histoGammaConvProb_MCPt,
                                            histoGammaPrimaryRecoEff_Pt,
                                            deltaEta, scaling, nEvt
                                            );
    }
    // correct alternate way for calo reco
    if (isCalo && !isPCM) {
        histoGammaCorrEffiReso_Pt                   = (TH1D*)histoESDCaloGammaPt->Clone("CorrGammaSpecPurityMinusSec");

        if (!hasCocktailInput)
            CorrectGammaEffiResol(  histoGammaCorrEffiReso_Pt,
                                    histoSecondaryGammaFromXSpecPt,
                                    histoGammaTruePurity_Pt,
                                    histoGammaPrimaryRecoEff_Pt,
                                    deltaEtaCalo, scalingCalo, nEvt
                                    );
        else
            CorrectGammaEffiResolCocktail(  histoGammaCorrEffiReso_Pt,
                                            histoGammaSecGammaFromX_Cocktail_Raw_Pt,      // already scaled per event
                                            histoGammaTrueSecCocktailGammaFromX_Pt[3],    // already scaled per event
                                            histoGammaTruePurity_Pt,
                                            histoGammaPrimaryRecoEff_Pt,
                                            deltaEtaCalo, scalingCalo, nEvt
                                            );
    }

    DrawGammaSetMarker(histoGammaCorrEffiReso_Pt, 20, 1.0, kGreen+2, kGreen+2);
    SetHistogramm(histoGammaCorrEffiReso_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                  histoGammaCorrEffiReso_Pt->GetMinimum(0)/100., histoGammaCorrEffiReso_Pt->GetMaximum()*10, 1.0,1.7);

    //******************************************************************************************
    TH1D* histoGammaCorrEffiReso_PileUp_Pt              = NULL;
    TH1D* histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt    = NULL;
    if(doPileUpCorr){
        //******************************************************************************************
        // copy raw inc gamma spectrum after pure data driven pileup correction,
        // - secondary contamination,
        // - pileup corrected purity, (pileup corrected, vs rec pT with sec correction)
        // - reco effi (vs rec pt with sec correction)
        // - conversion probability
        //******************************************************************************************
        histoGammaCorrEffiReso_PileUp_Pt                = (TH1D*)histoESDConvGammaPt->Clone("GammaCorrEffiResolPileup_Pt");
        histoGammaCorrEffiReso_PileUp_Pt->Multiply(histoPileUpCorrectionFactor_Pt);

        // instead now using correction factor to match unfolding strand
        //histoGammaCorrEffiReso_PileUp_Pt                = (TH1D*)histoESDConvGammaPtPileUp->Clone("CorrGammaSpecPurityMinusSecPileUp");

        if (!hasCocktailInput) {
            CorrectGammaEffiResol(  histoGammaCorrEffiReso_PileUp_Pt,
                                    histoSecondaryGammaFromXSpecPileUpPt,
                                    histoGammaTruePurity_PileUp_Pt,
                                    histoGammaConvProb_MCPt,
                                    histoGammaPrimaryRecoEff_PileUp_Pt,
                                    deltaEta, scaling, nEvt
                                    );
        } else {
            CorrectGammaEffiResolCocktail(  histoGammaCorrEffiReso_PileUp_Pt,
                                            histoGammaSecGammaFromX_Cocktail_Raw_Pt,      // cocktail spectra are from fully corrected, i.e. don't contain pileup contribution
                                            histoSecondaryGammaFromXSpecPileUpPt[3],
                                            histoGammaTruePurity_PileUp_Pt,
                                            histoGammaConvProb_MCPt,
                                            histoGammaPrimaryRecoEff_PileUp_Pt,
                                            deltaEta, scaling, nEvt
                                        );
        }

        DrawGammaSetMarker(histoGammaCorrEffiReso_PileUp_Pt, 20, 1.0, kMagenta, kMagenta);
        SetHistogramm(histoGammaCorrEffiReso_PileUp_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", -99, -99, 1.0, 1.7);

        //******************************************************************************************
        // copy raw inc gamma spectrum after pure data driven pileup correction,
        // - secondary contamination,
        // - purity,
        // - reco effi (vs rec pt with sec correction)
        // - conversion probability
        //******************************************************************************************
        histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt      = (TH1D*)histoESDConvGammaPt->Clone("GammaCorrEffiResolPileup_NoMCUpdate_Pt");
        histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt->Multiply(histoPileUpCorrectionFactor_Pt);

        // instead now using correction factor to match unfolding strand
        //histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt      = (TH1D*)histoESDConvGammaPtPileUp->Clone("CorrGammaSpecPurityMinusSecPileUpNoMCUpdate");

        if (!hasCocktailInput) {
            CorrectGammaEffiResol(  histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt,
                                    histoSecondaryGammaFromXSpecPt,
                                    histoGammaTruePurity_Pt,
                                    histoGammaConvProb_MCPt,
                                    histoGammaPrimaryRecoEff_Pt,
                                    deltaEta, scaling, nEvt
                                );
        } else {
            CorrectGammaEffiResolCocktail(  histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt,
                                            histoGammaSecGammaFromX_Cocktail_Raw_Pt,     // cocktail spectra are from fully corrected, i.e. don't contain pileup contribution
                                            histoGammaTrueSecCocktailGammaFromX_Pt[3],
                                            histoGammaTruePurity_Pt,
                                            histoGammaConvProb_MCPt,
                                            histoGammaPrimaryRecoEff_Pt,
                                            deltaEta, scaling, nEvt
                                        );
        }

        DrawGammaSetMarker(histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt, 20, 1.0, kBlue+2, kBlue+2);
        SetHistogramm(histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
    }

    //*************************************************************************************************
    //********** Sanity check for photon spectrum (do we get back what we put in) *********************
    //*************************************************************************************************
    // correct input inc gamma spectrum for delta Eta, delta phi & number of events
    TH1D* histoMCGammaSpec_MCPt                      = (TH1D*)histoMCAllGamma_MCPt->Clone("GammaSpecMC");
    if (isPCM )             CorrectGammaMC(histoMCGammaSpec_MCPt, deltaEta,       scaling,        nEvtMC);
    if (isCalo && !isPCM)   CorrectGammaMC(histoMCGammaSpec_MCPt, deltaEtaCalo,   scalingCalo,    nEvtMC);

    //*************************************************************************************************
    //********************************* Comparison Gamma Spec *****************************************
    //*************************************************************************************************
    TH1D* histoRatioGammaCorrMinusSecGammaUnfold            = (TH1D*)histoGammaCorrEffiReso_Pt->Clone("histoRatioGammaCorrMinusSecGammaUnfold");
    histoRatioGammaCorrMinusSecGammaUnfold->Divide(histoGammaCorrUnfoldReso_Pt);
    TH1D* histoRatioGammaCorrMinusSecGammaPileUp            = NULL;
    TH1D* histoRatioGammaCorrMinusSecGammaPileUpNoMCUpdate  = NULL;
    if(doPileUpCorr){
        histoRatioGammaCorrMinusSecGammaPileUp              = (TH1D*)histoGammaCorrEffiReso_Pt->Clone("histoRatioGammaCorrMinusSecGammaPileUp");
        histoRatioGammaCorrMinusSecGammaPileUp->Divide(histoGammaCorrEffiReso_PileUp_Pt);

        histoRatioGammaCorrMinusSecGammaPileUpNoMCUpdate    = (TH1D*)histoGammaCorrEffiReso_Pt->Clone("histoRatioGammaCorrMinusSecGammaPileUp");
        histoRatioGammaCorrMinusSecGammaPileUpNoMCUpdate->Divide(histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt);
    }

    //*************************************************************************************************
    //********************************* Decay gamma to all gamma **************************************
    //*************************************************************************************************
    TCanvas *canvasRatioAllDiffDecay    = GetAndSetCanvas("canvasRatioAllDiffDecay");
    canvasRatioAllDiffDecay->SetLogy(0);
    canvasRatioAllDiffDecay->cd();

        TH1D* histoAllDiffDecayGamma    = (TH1D*)histoMCAllGamma_OriginalBin_MCPt->Clone("RatioAllToDecay");
        histoAllDiffDecayGamma->Scale(histoAllDiffDecayGamma->GetBinWidth(5));
        histoAllDiffDecayGamma->Rebin(5);
        histoAllDiffDecayGamma->Scale(1./histoAllDiffDecayGamma->GetBinWidth(5));

        TH1D* histoAllDecayGamma        = (TH1D*)histoPhotonSource_MCPt[7]->Clone("AllDecay");
        histoAllDecayGamma->Scale(histoAllDecayGamma->GetBinWidth(5));
        histoAllDecayGamma->Rebin(5);
        histoAllDecayGamma->Scale(1./histoAllDecayGamma->GetBinWidth(5));

        histoAllDiffDecayGamma->Divide(histoAllDiffDecayGamma,histoAllDecayGamma,1,1,"B");
        histoAllDiffDecayGamma->GetYaxis()->SetRangeUser(0.8,1.5);

        SetHistogramm(histoAllDiffDecayGamma,"#it{p}_{T} (GeV/#it{c})", "#frac{N_{#gamma_{incl}}}{N_{#gamma_{decay}}}");
        histoAllDiffDecayGamma->Draw("p,e1");

        DrawGammaLines(0., 50,1, 1,0.5, kGray+2, 7);

        PutProcessLabelAndEnergyOnPlot( 0.95, 0.95, 0.035, cent, "#gamma", "", 42, 0.03,"", 1, 1.25,31);

    canvasRatioAllDiffDecay->SaveAs(Form("%s/%s_AllGammaDivDecayGammaSpectrumMC_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasRatioAllDiffDecay;

    //*************************************************************************************************
    //********************************* Gamma from Decay **********************************************
    //*************************************************************************************************
    // correct input spectra of different decay channels with delta eta, delta phi & events number
    for (Int_t i = 0; i< 9; i++){
        CorrectGammaMC(histoPhotonSource_MCPt[i], deltaEta, scaling, nEvtMC);
        SetHistogramm(histoPhotonSource_MCPt[i],"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",-99., -99., 1.0, 1.65);
        DrawGammaSetMarker(histoPhotonSource_MCPt[i], 22, 1.0, colorsDecay[i] , colorsDecay[i]);
    }

    TCanvas *canvasDecayGammaSpecMC     = new TCanvas("canvasDecayGammaSpecMC","",200,10,1000*1.25,1100*1.25);  // gives the page size
    DrawGammaCanvasSettings( canvasDecayGammaSpecMC, 0.16, 0.02, 0.02, 0.09);
    canvasDecayGammaSpecMC->SetLogy();

        TLegend* legendDecaySpectra     = GetAndSetLegend2(0.5, 0.935-0.035*1.1*5, 0.95, 0.935, 0.035, 2, "", 42, 0.2);
        histoPhotonSource_MCPt[7]->DrawCopy("chist");
        for (Int_t i = 0; i< 9; i++){
            histoPhotonSource_MCPt[i]->DrawCopy("chist,same");
            legendDecaySpectra->AddEntry(histoPhotonSource_MCPt[i],decaysLabels[i].Data(),"l");
        }
        legendDecaySpectra->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.94, 0.935-0.035*1.1*5, 0.035, cent, detectionProcess,"", 42, 0.03,"", 1, 1.25,31);

    canvasDecayGammaSpecMC->SaveAs(Form("%s/%s_DecayGammaSpectrumMC_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasDecayGammaSpecMC;

    //*************************************************************************************************
    //**************** Compare cocktail secondaries to MC secondaries *********************************
    //*************************************************************************************************
    if(hasCocktailInput){
        // divide cocktail secondary spectra by MC ones
        TH1D* histoRatioSecondariesCocktailMC_Raw_Pt[3]     = {NULL, NULL, NULL};
        Double_t ratioSecondariesCocktailMCYMax             = 1.0;
        Double_t ratioSecondariesCocktailMCYMin             = 0.0;
        for (Int_t k = 0; k < 3; k++) {
            histoRatioSecondariesCocktailMC_Raw_Pt[k]       = (TH1D*)histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]->Clone(Form("histoRatioSecondariesCocktailMC%s_Raw_Pt", nameSecondaries[k].Data()));
            if (isPCM )             histoRatioSecondariesCocktailMC_Raw_Pt[k]->Divide(histoGammaTrueSecConvGammaFromX_Pt[k]);
            if (isCalo && !isPCM)   histoRatioSecondariesCocktailMC_Raw_Pt[k]->Divide(histoGammaTrueSecCaloGammaFromX_Pt[k]);

            // get y-range min
            if (k==0 || histoRatioSecondariesCocktailMC_Raw_Pt[k]->GetMinimum(0) < ratioSecondariesCocktailMCYMin)
                ratioSecondariesCocktailMCYMin              = histoRatioSecondariesCocktailMC_Raw_Pt[k]->GetMinimum(0);

            // get y-range max
            if (k==0 || histoRatioSecondariesCocktailMC_Raw_Pt[k]->GetMaximum() > ratioSecondariesCocktailMCYMax)
                ratioSecondariesCocktailMCYMax              = histoRatioSecondariesCocktailMC_Raw_Pt[k]->GetMaximum();
        }
        ratioSecondariesCocktailMCYMin                      = ratioSecondariesCocktailMCYMin;
        ratioSecondariesCocktailMCYMax                      = ratioSecondariesCocktailMCYMax;

        TCanvas *canvasSecondaryComparison                  = GetAndSetCanvas("canvasSecondaryComparison");
        DrawGammaCanvasSettings( canvasSecondaryComparison, 0.07, 0.02, 0.02, 0.085);
        TLegend* legendCompareSecCocktailMC                 = GetAndSetLegend2(0.75, 0.935-0.035*1.1*3, 0.95, 0.935, 0.035, 1, "", 42, 0.15);
        legendCompareSecCocktailMC->SetBorderSize(0);

        for (Int_t k = 0; k < 3; k++) {
            SetHistogramm(histoRatioSecondariesCocktailMC_Raw_Pt[k],"#it{p}_{T} (GeV/#it{c})", "Cocktail/MC",-99,-99,1.0,0.9);

            // set y-range depending on max. ratio value
            if (ratioSecondariesCocktailMCYMax > 2.0) {
                canvasSecondaryComparison->SetLogy();
                histoRatioSecondariesCocktailMC_Raw_Pt[k]->GetYaxis()->SetRangeUser(ratioSecondariesCocktailMCYMin*0.1, ratioSecondariesCocktailMCYMax*5.0);
            } else {
                histoRatioSecondariesCocktailMC_Raw_Pt[k]->GetYaxis()->SetRangeUser(0.0, ratioSecondariesCocktailMCYMax*1.5);
            }

            DrawGammaSetMarker(histoRatioSecondariesCocktailMC_Raw_Pt[k], markerStyleSec[k], markerSizeSec[k]*2, colorSecFromToy[k] , colorSecFromToy[k]);
            if (k == 0) histoRatioSecondariesCocktailMC_Raw_Pt[k]->Draw();
            else        histoRatioSecondariesCocktailMC_Raw_Pt[k]->Draw("same");
            legendCompareSecCocktailMC->AddEntry(histoRatioSecondariesCocktailMC_Raw_Pt[k],Form("Sec. #gamma from %s", nameLabelSecondaries[k].Data()),"p");
        }
        legendCompareSecCocktailMC->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.935, 0.24, 0.035, cent, detectionProcess,"", 42, 0.03,"", 1, 1.25,31);

        canvasSecondaryComparison->SaveAs(Form("%s/%s_SecondaryComparisonCocktailMC_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasSecondaryComparison;
    }

    //*************************************************************************************************
    //******************************* Gamma Combinatorial Background **********************************
    //*************************************************************************************************
    // pure BG plotting
    TCanvas *canvasCombBackSpecMC       = new TCanvas("canvasCombBackSpecMC","",200,10,1000*1.25,1100*1.25);  // gives the page size
    DrawGammaCanvasSettings( canvasCombBackSpecMC, 0.16, 0.02, 0.02, 0.09);
    canvasCombBackSpecMC->SetLogy();

        TLegend* legendCombSpectra      = NULL;
        if (isPCM ) {
            legendCombSpectra           = GetAndSetLegend2(0.35,0.945-4*1.1*0.035, 0.95,0.945, 0.035, 5, "", 42, 0.1);

            SetHistogramm(histoCombinatorialSpecies_Pt[16],"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",-99., -99., 1.0, 1.7);
            DrawGammaSetMarker(histoCombinatorialSpecies_Pt[16], markersCombinatorics[16], 1., colorsCombinatorics[16], colorsCombinatorics[16]);
            histoCombinatorialSpecies_Pt[16]->Scale(1./nEvtMC);
            histoCombinatorialSpecies_Pt[16]->GetYaxis()->SetRangeUser(histoCombinatorialSpecies_Pt[16]->GetMinimum(0)*1e-3, histoCombinatorialSpecies_Pt[16]->GetMaximum()*10.);
            histoCombinatorialSpecies_Pt[16]->DrawCopy("");
            legendCombSpectra->AddEntry(histoCombinatorialSpecies_Pt[16],combinatoricsLabels[16]);

            for(Int_t i = 0;i<16;i++){
                SetHistogramm(histoCombinatorialSpecies_Pt[i],"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",-99., -99., 1.0, 1.7);
                DrawGammaSetMarker(histoCombinatorialSpecies_Pt[i], markersCombinatorics[i], 1., colorsCombinatorics[i], colorsCombinatorics[i]);
                histoCombinatorialSpecies_Pt[i]->Scale(1./nEvtMC);
                histoCombinatorialSpecies_Pt[i]->DrawCopy("same");
                legendCombSpectra->AddEntry(histoCombinatorialSpecies_Pt[i],combinatoricsLabels[i]);
            }
        }
        if (isCalo && !isPCM) {
            legendCombSpectra           = GetAndSetLegend2(0.45,0.945-3*1.1*0.035, 0.95,0.945, 0.035, 5, "", 42, 0.13);
            SetHistogramm(histoCombinatorialSpeciesCalo_Pt[10],"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",-99., -99., 1.0, 1.7);
            DrawGammaSetMarker(histoCombinatorialSpeciesCalo_Pt[10], markersCombinatoricsCalo[10], 1., colorsCombinatoricsCalo[10], colorsCombinatoricsCalo[10]);
            histoCombinatorialSpeciesCalo_Pt[10]->Scale(1./nEvtMC);
            histoCombinatorialSpeciesCalo_Pt[10]->GetYaxis()->SetRangeUser(histoCombinatorialSpeciesCalo_Pt[10]->GetMinimum(0)*1e-3, histoCombinatorialSpeciesCalo_Pt[10]->GetMaximum()*10.);
            histoCombinatorialSpeciesCalo_Pt[10]->DrawCopy("");
            legendCombSpectra->AddEntry(histoCombinatorialSpeciesCalo_Pt[10],combinatoricsLabelsCalo[10]);

            for(Int_t i = 0;i<10;i++){
                SetHistogramm(histoCombinatorialSpeciesCalo_Pt[i],"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",-99., -99., 1.0, 1.7);
                DrawGammaSetMarker(histoCombinatorialSpeciesCalo_Pt[i], markersCombinatoricsCalo[i], 1., colorsCombinatoricsCalo[i], colorsCombinatoricsCalo[i]);
                histoCombinatorialSpeciesCalo_Pt[i]->Scale(1./nEvtMC);
                histoCombinatorialSpeciesCalo_Pt[i]->DrawCopy("same");
                legendCombSpectra->AddEntry(histoCombinatorialSpeciesCalo_Pt[i],combinatoricsLabelsCalo[i]);
            }
        }

        legendCombSpectra->Draw();
        PutProcessLabelAndEnergyOnPlot( 0.935, 0.93-4*1.1*0.035, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,31);

    canvasCombBackSpecMC->SaveAs(Form("%s/%s_CombBackSpectrumMC_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasCombBackSpecMC;

    // plot ratio to real primary photons
    TCanvas *canvasSignalToCombBackgroundRatio              = new TCanvas("canvasSignalToCombBackgroundRatio", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasSignalToCombBackgroundRatio,  0.095, 0.015, 0.015, 0.095);
    canvasSignalToCombBackgroundRatio->SetLogy();

        TLegend* legendSignalToCombBackgroundRatio  = NULL;
        TH1D **histoSignalToCombBackgroundRatio     = NULL;
        TH1D *SummedSmallContributionsCombBack      = NULL; //10,11,12,13,14,15
        if (isPCM ) {
            histoSignalToCombBackgroundRatio        = new TH1D*[17];

            legendSignalToCombBackgroundRatio       = GetAndSetLegend2(0.15,0.945-2*1.1*0.035, 0.95,0.945, 0.035, 6, "", 42, 0.1);
            for(Int_t i = 0;i<16;i++){
                histoSignalToCombBackgroundRatio[i] = (TH1D*) histoCombinatorialSpecies_Pt[i]->Clone(Form("ESD_TrueCombRatioSignal_%s_Pt",combinatorics[i].Data()));
                histoSignalToCombBackgroundRatio[i]->Scale(nEvtMC);

                if(i==9){
                    SummedSmallContributionsCombBack = (TH1D*)histoSignalToCombBackgroundRatio[i]->Clone("SummedSmallContributions");
                    SetHistogramm(SummedSmallContributionsCombBack,"#it{p}_{T} (GeV/#it{c})","SummedSmallContributions",10,5e7);
                    SummedSmallContributionsCombBack->SetMinimum(1e-5);
                } else if(i>9){
                    SummedSmallContributionsCombBack->Add(histoSignalToCombBackgroundRatio[i]);
                }
                histoSignalToCombBackgroundRatio[i]->Divide(histoSignalToCombBackgroundRatio[i],histoGammaTruePrimaryConv_Pt,1,1,"");
                SetHistogramm(histoSignalToCombBackgroundRatio[i],"#it{p}_{T} (GeV/#it{c})","#it{C}_{i}",10,5e7);
                histoSignalToCombBackgroundRatio[i]->SetMinimum(1e-5);

                if(i==0){
                    histoSignalToCombBackgroundRatio[i]->GetYaxis()->SetRangeUser(1e-5,100);
                    if(energy.BeginsWith("8TeV") && (mode==2 || mode==4)) histoSignalToCombBackgroundRatio[i]->GetYaxis()->SetRangeUser(1e-5,1);
                    histoSignalToCombBackgroundRatio[i]->DrawCopy("e1");
                } else if(i<9){
                    DrawGammaSetMarker(histoSignalToCombBackgroundRatio[i], markersCombinatorics[i], 1., colorsCombinatorics[i], colorsCombinatorics[i]);
                    histoSignalToCombBackgroundRatio[i]->DrawCopy("e1same");
                } else continue;

                legendSignalToCombBackgroundRatio->AddEntry(histoSignalToCombBackgroundRatio[i],combinatoricsLabels[i]);
            }

            SummedSmallContributionsCombBack->Divide(SummedSmallContributionsCombBack,histoGammaTruePrimaryConv_Pt,1,1,"");
            legendSignalToCombBackgroundRatio->AddEntry(SummedSmallContributionsCombBack,"p(#bar{p})K^{#pm}#mu^{#pm}");
        }
        if (isCalo && !isPCM) {
            histoSignalToCombBackgroundRatio        = new TH1D*[10];
            legendSignalToCombBackgroundRatio       = GetAndSetLegend2(0.15,0.93-1*1.1*0.035, 0.95,0.93, 0.035, 9, "", 42, 0.15);

            for(Int_t i = 0;i<10;i++){
                histoSignalToCombBackgroundRatio[i] = (TH1D*)histoCombinatorialSpeciesCalo_Pt[i]->Clone(Form("ESD_TrueCombRatioSignal_%s_Pt",combinatoricsCalo[i].Data()));
                histoSignalToCombBackgroundRatio[i]->Scale(nEvtMC);

                if(i==7){
                    SummedSmallContributionsCombBack = (TH1D*)histoSignalToCombBackgroundRatio[i]->Clone("SummedSmallContributions");
                    SetHistogramm(SummedSmallContributionsCombBack,"#it{p}_{T} (GeV/#it{c})","SummedSmallContributions",10,5e7);
                    SummedSmallContributionsCombBack->SetMinimum(1e-5);
                } else if(  i==9){
                    SummedSmallContributionsCombBack->Add(histoSignalToCombBackgroundRatio[i]);
                }

                histoSignalToCombBackgroundRatio[i]->Divide(histoSignalToCombBackgroundRatio[i],histoGammaTruePrimaryCalo_Pt,1,1,"");
                SetHistogramm(histoSignalToCombBackgroundRatio[i],"#it{p}_{T} (GeV/#it{c})","#it{C}_{i}",10,5e7);
                histoSignalToCombBackgroundRatio[i]->SetMinimum(1e-5);

                if(i==0){
                    histoSignalToCombBackgroundRatio[i]->GetYaxis()->SetRangeUser(1e-5,100);
                    histoSignalToCombBackgroundRatio[i]->DrawCopy("e1");
                } else if( i<7 || i==8){
                    DrawGammaSetMarker(histoSignalToCombBackgroundRatio[i], markersCombinatoricsCalo[i], 1., colorsCombinatoricsCalo[i], colorsCombinatoricsCalo[i]);
                    histoSignalToCombBackgroundRatio[i]->DrawCopy("e1same");
                } else continue;

                legendSignalToCombBackgroundRatio->AddEntry(histoSignalToCombBackgroundRatio[i],combinatoricsLabelsCalo[i]);
            }

            SummedSmallContributionsCombBack->Divide(SummedSmallContributionsCombBack,histoGammaTruePrimaryCalo_Pt,1,1,"");
            legendSignalToCombBackgroundRatio->AddEntry(SummedSmallContributionsCombBack,"#mu+rest");
        }

        DrawGammaSetMarker(SummedSmallContributionsCombBack, markersCombinatorics[10], 1., colorsCombinatorics[10], colorsCombinatorics[10]);
        SummedSmallContributionsCombBack->DrawCopy("e1same");
        legendSignalToCombBackgroundRatio->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.15, 0.945-0.035*1.05*2, 0.035, cent, detectionProcess, "" , 42, 0.03,"",1,1.25,11);

    canvasSignalToCombBackgroundRatio->SaveAs(Form("%s/%s_CombBackgroundRatioToSignal_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasSignalToCombBackgroundRatio;

    // plot ratio to summed total MC BG
    TCanvas *canvasRatioCombBackToBack              = new TCanvas("canvasRatioCombBackToBack", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioCombBackToBack,  0.095, 0.015, 0.015, 0.095);
    canvasRatioCombBackToBack->SetLogy();

        TLegend* legendRatioCombBackToBack          = NULL;

        TH1D **histoRatioCombBackToBack                 = NULL;
        TH1D *SummedSmallContributionsCombBackToBack    = NULL; //10,11,12,13,14,15
        if (isPCM ) {
            histoRatioCombBackToBack                    = new TH1D*[17];
            legendRatioCombBackToBack                   = GetAndSetLegend2(0.15,0.945-2*1.1*0.035, 0.95,0.945, 0.035, 6, "", 42, 0.1);
            for(Int_t i = 0;i<16;i++){

                histoRatioCombBackToBack[i] = (TH1D*) histoCombinatorialSpecies_Pt[i]->Clone(Form("ESD_TrueCombRatioSignal_%s_Pt",combinatorics[i].Data()));

                if(i==9){
                    SummedSmallContributionsCombBackToBack = (TH1D*)histoRatioCombBackToBack[i]->Clone("SummedSmallContributions");
                    SetHistogramm(SummedSmallContributionsCombBackToBack,"#it{p}_{T} (GeV/#it{c})","SummedSmallContributions",10,5e7);
                    SummedSmallContributionsCombBackToBack->SetMinimum(1e-5);
                } else if(i>9){
                    SummedSmallContributionsCombBackToBack->Add(histoRatioCombBackToBack[i]);
                }
                histoRatioCombBackToBack[i]->Divide(histoRatioCombBackToBack[i],histoGammaMCBackground_Pt,1,1,"B");
                SetHistogramm(histoRatioCombBackToBack[i],"#it{p}_{T} (GeV/#it{c})","#it{K}_{i}",10,5e7);
                histoRatioCombBackToBack[i]->SetMinimum(1e-5);

                if(i==0){
                    histoRatioCombBackToBack[i]->GetYaxis()->SetRangeUser(1e-4,40);
                    histoRatioCombBackToBack[i]->DrawCopy("e1");
                } else if(i<9){
                    DrawGammaSetMarker(histoRatioCombBackToBack[i], markersCombinatorics[i], 1., colorsCombinatorics[i], colorsCombinatorics[i]);
                    histoRatioCombBackToBack[i]->DrawCopy("e1same");
                } else continue;

                legendRatioCombBackToBack->AddEntry(histoRatioCombBackToBack[i],combinatoricsLabels[i]);
            }

            legendRatioCombBackToBack->AddEntry(SummedSmallContributionsCombBackToBack,"p(#bar{p})K^{#pm}#mu^{#pm}");
        }
        if (isCalo && !isPCM) {
            histoRatioCombBackToBack                    = new TH1D*[10];
            legendRatioCombBackToBack                   = GetAndSetLegend2(0.15,0.93-1*1.1*0.035, 0.95,0.93, 0.035, 9, "", 42, 0.15);
            for(Int_t i = 0;i<10;i++){
                histoRatioCombBackToBack[i]             = (TH1D*)histoCombinatorialSpeciesCalo_Pt[i]->Clone(Form("ESD_TrueCombRatioSignal_%s_Pt",combinatoricsCalo[i].Data()));

                if(i==7){
                    SummedSmallContributionsCombBackToBack = (TH1D*)histoRatioCombBackToBack[i]->Clone("SummedSmallContributions");
                    SetHistogramm(SummedSmallContributionsCombBackToBack,"#it{p}_{T} (GeV/#it{c})","SummedSmallContributions",10,5e7);
                    SummedSmallContributionsCombBackToBack->SetMinimum(1e-5);
                } else if( i==9){
                    SummedSmallContributionsCombBackToBack->Add(histoRatioCombBackToBack[i]);
                }

                histoRatioCombBackToBack[i]->Divide(histoRatioCombBackToBack[i],histoGammaMCBackground_Pt,1,1,"B");
                SetHistogramm(histoRatioCombBackToBack[i],"#it{p}_{T} (GeV/#it{c})","#it{K}_{i}",10,5e7);
                histoRatioCombBackToBack[i]->SetMinimum(1e-5);

                if(i==0){
                    histoRatioCombBackToBack[i]->GetYaxis()->SetRangeUser(1e-4,40);
                    histoRatioCombBackToBack[i]->DrawCopy("e1");
                } else if( i<7 || i==8){
                    DrawGammaSetMarker(histoRatioCombBackToBack[i], markersCombinatoricsCalo[i], 1., colorsCombinatoricsCalo[i], colorsCombinatoricsCalo[i]);
                    histoRatioCombBackToBack[i]->DrawCopy("e1same");
                } else continue;

                legendRatioCombBackToBack->AddEntry(histoRatioCombBackToBack[i],combinatoricsLabelsCalo[i]);
            }

            legendRatioCombBackToBack->AddEntry(SummedSmallContributionsCombBackToBack,"#mu+rest");
        }

        SummedSmallContributionsCombBackToBack->Divide(SummedSmallContributionsCombBackToBack,histoGammaMCBackground_Pt,1,1,"B");
        DrawGammaSetMarker(SummedSmallContributionsCombBackToBack, markersCombinatorics[10], 1., colorsCombinatorics[10], colorsCombinatorics[10]);
        SummedSmallContributionsCombBackToBack->DrawCopy("e1same");
        legendRatioCombBackToBack->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.15, 0.945-0.035*1.05*2, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,11);

    canvasRatioCombBackToBack->SaveAs(Form("%s/%s_RatioCombBackToBack_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));


    canvasRatioCombBackToBack->cd();
    if (isPCM ) {
        for(Int_t i = 0;i<16;i++){
            if(i==0){
                histoRatioCombBackToBack[i]->GetXaxis()->SetRangeUser(1.,3.);
                histoRatioCombBackToBack[i]->GetYaxis()->SetRangeUser(1e-4,40);
                histoRatioCombBackToBack[i]->DrawCopy("e1");
            }
            if(i<9){
                DrawGammaSetMarker(histoRatioCombBackToBack[i], markersCombinatorics[i], 1., colorsCombinatorics[i], colorsCombinatorics[i]);
                histoRatioCombBackToBack[i]->DrawCopy("e1same");
            }
            else continue;
        }
        SummedSmallContributionsCombBackToBack->DrawCopy("e1same");
        legendRatioCombBackToBack->Draw();
        PutProcessLabelAndEnergyOnPlot( 0.15, 0.945-0.035*1.05*2, 0.035, cent, detectionProcess, "", 42, 0.03,"",1,1.25,11);
        canvasRatioCombBackToBack->SaveAs(Form("%s/%s_RatioCombBackToBackZoomed_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    }
    delete canvasRatioCombBackToBack;

    //*************************************************************************************************
    //***************** Compare Gamma Yields alternate corr method vs unfolded ************************
    //*************************************************************************************************
    TCanvas *canvaskGammaSpecAlternateCorrMethods   = new TCanvas("canvaskGammaSpecAlternateCorrMethods","",200,10,1000*1.25,1300*1.25);  // gives the page size
    DrawGammaCanvasSettings(canvaskGammaSpecAlternateCorrMethods,0.16, 0.015, 0.02, 0.09);
    TPad* padAlternateCorrMethods       = new TPad("padAlternateCorrMethods", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings(padAlternateCorrMethods, 0.16, 0.015, 0.02, 0.);
    padAlternateCorrMethods->Draw();
    TPad* padAlternateCorrMethodsRatio  = new TPad("padAlternateCorrMethodsRatio", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings(padAlternateCorrMethodsRatio,  0.16, 0.015, 0.0, 0.27);
    padAlternateCorrMethodsRatio->Draw();

    padAlternateCorrMethodsRatio->Draw();
    padAlternateCorrMethods->cd();
    padAlternateCorrMethods->SetLogy();

        TLegend* legendAlternateCorrMethods = GetAndSetLegend(0.50,0.8,3);
        // plot unfolded spectrum
        histoGammaCorrUnfoldReso_Pt->DrawCopy("e1");
        legendAlternateCorrMethods->AddEntry(histoGammaCorrUnfoldReso_Pt,"#gamma with Bayes unfolding","p");
        // copy pileup correct hist if pileup correction enabled
        if (doPileUpCorr){
            DrawGammaSetMarker(histoGammaCorrEffiReso_PileUp_Pt,  24, 1.0, kAzure+7,  kAzure+7);
            histoGammaCorrEffiReso_PileUp_Pt->DrawCopy("e1,same");
            legendAlternateCorrMethods->AddEntry(histoGammaCorrEffiReso_PileUp_Pt,"#gamma with bin-by-bin MC","p");
            legendAlternateCorrMethods->AddEntry((TObject*)0,"based corrections","");
        // take hist without pileup corr if disabled
        } else {
            DrawGammaSetMarker(histoGammaCorrEffiReso_Pt,  24, 1.0, kAzure+7,  kAzure+7);
            histoGammaCorrEffiReso_Pt->DrawCopy("e1,same");
            legendAlternateCorrMethods->AddEntry(histoGammaCorrEffiReso_Pt,"#gamma with bin-by-bin MC","p");
            legendAlternateCorrMethods->AddEntry((TObject*)0,"based corrections","");
        }
        legendAlternateCorrMethods->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.20, 0.15, 0.035, cent, detectionProcess,"", 42, 0.03);

    padAlternateCorrMethodsRatio->cd();

        TH1D*               histoRatioGammaLegacyCorrVsUnfold   = NULL;
        if (doPileUpCorr)   histoRatioGammaLegacyCorrVsUnfold   = (TH1D*)histoGammaCorrEffiReso_PileUp_Pt->Clone("histoRatioGammaLegacyCorrVsUnfold");
        else                histoRatioGammaLegacyCorrVsUnfold   = (TH1D*)histoGammaCorrEffiReso_Pt->Clone("histoRatioGammaLegacyCorrVsUnfold");
        histoRatioGammaLegacyCorrVsUnfold->Divide(histoRatioGammaLegacyCorrVsUnfold,histoGammaCorrUnfoldReso_Pt,1,1,"");
        histoRatioGammaLegacyCorrVsUnfold->GetYaxis()->SetRangeUser(0.88,1.12);

        SetStyleHistoTH1ForGraphs(  histoRatioGammaLegacyCorrVsUnfold, "#it{p}_{T} (GeV/#it{c})","#frac{#gamma corr legacy}{#gamma corr unfold}" , 0.10, 0.12,
                                    0.10, 0.12,  0.85, 0.6, 510, 505);

        DrawGammaSetMarker(histoRatioGammaLegacyCorrVsUnfold, 24, 1.0, kAzure+7,  kAzure+7);
        histoRatioGammaLegacyCorrVsUnfold->DrawCopy("");

        DrawGammaLines(minPtGamma, maxPtGamma,1, 1,0.5);

    canvaskGammaSpecAlternateCorrMethods->SaveAs(Form("%s/%s_%s_GammaSpectraComparison_%s.%s",outputDir.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
    delete canvaskGammaSpecAlternateCorrMethods;


    //*************************************************************************************************
    //***************** Compare reconstructed gamma yields to MC input ********************************
    //*************************************************************************************************
    if (isPCM ) {
        CorrectGammaMC(histoMCGammaConv_MCPt, deltaEta, scaling, nEvtMC);
        CorrectGammaMC(histoGammaTrueConv_Pt, deltaEta, scaling, nEvtMC);
    }
    if (isCalo && !isPCM) {
        // this might not be true since deltaEtaCalo and scalingCalo are read from the cutstring, doesn't take into account missing modules for e.g. 2010 data!
        CorrectGammaMC(histoGammaTrueCalo_Pt, deltaEtaCalo, scalingCalo, nEvtMC);
    }

    TCanvas *canvasGammaPurityConvBin   = new TCanvas(  "canvasGammaPurityConvBin",     "",200,10,1000*1.25,1300*1.25);  // gives the page size
    TPad* padGammaPurityConvBin         = new TPad(     "padGammaPurityConvBin",        "", 0., 0.25, 1., 1.,-1, -1, -2);
    TPad* padGammaPurityConvBinRatio    = new TPad(     "padGammaPurityConvBinRatio",   "", 0., 0., 1., 0.25,-1, -1, -2);

    DrawGammaCanvasSettings(canvasGammaPurityConvBin,   0.16, 0.015, 0.02,  0.09);
    DrawGammaPadSettings(   padGammaPurityConvBin,      0.16, 0.015, 0.02,  0.);
    DrawGammaPadSettings(   padGammaPurityConvBinRatio, 0.16, 0.015, 0.0,   0.27);

    padGammaPurityConvBin->Draw();
    padGammaPurityConvBinRatio->Draw();

    padGammaPurityConvBin->cd();
    padGammaPurityConvBin->SetLogy();

        DrawGammaSetMarker(histoMCGammaSpec_MCPt,       20, 1.0, kBlack, kBlack);
        DrawGammaSetMarker(histoGammaCorrUnfoldReso_Pt, 25, 1.0, kBlue+2, kBlue+2);
        if (doPileUpCorr)   DrawGammaSetMarker(histoGammaCorrEffiReso_PileUp_Pt,    28, 1.5, kRed-2, kRed-2);
        else                DrawGammaSetMarker(histoGammaCorrEffiReso_Pt,           28, 1.5, kRed-2, kRed-2);

        SetHistogramm(histoMCGammaSpec_MCPt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", -99, -99, 1.0, 1.7);
        histoMCGammaSpec_MCPt->GetYaxis()->SetRangeUser(histoMCGammaSpec_MCPt->GetMinimum(0)/100.,histoMCGammaSpec_MCPt->GetMaximum()*100);

        histoMCGammaSpec_MCPt->DrawCopy("e1");
        histoGammaCorrUnfoldReso_Pt->DrawCopy("e1,same");
        if (doPileUpCorr)   histoGammaCorrEffiReso_PileUp_Pt->DrawCopy("e1,same");
        else                histoGammaCorrEffiReso_Pt->DrawCopy("e1,same");

        TLegend *legendGammaSpectraConvBin          = GetAndSetLegend(0.325,0.82,3);
        if (isRunMC) {
            legendGammaSpectraConvBin->AddEntry(histoGammaCorrUnfoldReso_Pt,"corr. rec. MC #gamma spectrum (unfolding)");
            if (doPileUpCorr)   legendGammaSpectraConvBin->AddEntry(histoGammaCorrEffiReso_PileUp_Pt,   "corr. rec. MC #gamma spectrum (binwise)");
            else                legendGammaSpectraConvBin->AddEntry(histoGammaCorrEffiReso_Pt,          "corr. rec. MC #gamma spectrum (binwise)");
        } else {
            legendGammaSpectraConvBin->AddEntry(histoGammaCorrUnfoldReso_Pt,"corr. data #gamma spectrum (unfolding)");
            if (doPileUpCorr)   legendGammaSpectraConvBin->AddEntry(histoGammaCorrEffiReso_PileUp_Pt,   "corr. data #gamma spectrum (binwise)");
            else                legendGammaSpectraConvBin->AddEntry(histoGammaCorrEffiReso_Pt,          "corr. data #gamma spectrum (binwise)");
        }
        legendGammaSpectraConvBin->AddEntry(histoMCGammaSpec_MCPt,"input MC #gamma spectrum");
        legendGammaSpectraConvBin->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.20, 0.15, 0.035, cent, detectionProcess, "", 42, 0.03);

    padGammaPurityConvBinRatio->cd();

        TH1D* tempRatioDataMCUnfold                 = NULL;
        tempRatioDataMCUnfold                       = (TH1D*)histoGammaCorrUnfoldReso_Pt->Clone("tempRatioDataMCUnfold");
        tempRatioDataMCUnfold->Divide(tempRatioDataMCUnfold,histoMCGammaSpec_MCPt,1.,1.,"");

        TH1D*               tempRatioDataMCBinwise  = NULL;
        if (doPileUpCorr)   tempRatioDataMCBinwise  = (TH1D*)histoGammaCorrEffiReso_PileUp_Pt->Clone("tempRatioDataMCBinwise");
        else                tempRatioDataMCBinwise  = (TH1D*)histoGammaCorrEffiReso_Pt->Clone("tempRatioDataMCBinwise");
        tempRatioDataMCBinwise->Divide(tempRatioDataMCBinwise,histoMCGammaSpec_MCPt,1.,1.,"");

        DrawGammaSetMarker(tempRatioDataMCUnfold,   25, 1.0, kBlue+2, kBlue+2);
        DrawGammaSetMarker(tempRatioDataMCBinwise,  28, 1.5, kRed-2, kRed-2);

        if (isRunMC) {
            SetStyleHistoTH1ForGraphs(tempRatioDataMCUnfold,    "#it{p}_{T} (GeV/#it{c})", "#frac{rec. MC}{MC}", 0.10, 0.14, 0.10, 0.14, 0.85, 0.5, 510, 505);
            SetStyleHistoTH1ForGraphs(tempRatioDataMCBinwise,   "#it{p}_{T} (GeV/#it{c})", "#frac{rec. MC}{MC}", 0.10, 0.14, 0.10, 0.14, 0.85, 0.5, 510, 505);
        } else {
            SetStyleHistoTH1ForGraphs(tempRatioDataMCUnfold,    "#it{p}_{T} (GeV/#it{c})", "#frac{data}{MC}",    0.10, 0.14, 0.10, 0.14, 0.85, 0.5, 510, 505);
            SetStyleHistoTH1ForGraphs(tempRatioDataMCBinwise,   "#it{p}_{T} (GeV/#it{c})", "#frac{data}{MC}",    0.10, 0.14, 0.10, 0.14, 0.85, 0.5, 510, 505);
        }

        Double_t minYLines, maxYLines;
        if (isRunMC){
            tempRatioDataMCUnfold->GetYaxis()->SetRangeUser(0.85, 1.15);
            minYLines                               = 0.9;
            maxYLines                               = 1.1;
        } else if (energy.Contains("pPb")) {
            tempRatioDataMCUnfold->GetYaxis()->SetRangeUser(0.5, 2.7);
            minYLines                               = 0.8;
            maxYLines                               = 1.2;
        } else {
            tempRatioDataMCUnfold->GetYaxis()->SetRangeUser(0.5, 1.5);
            if(energy.BeginsWith("8TeV") && mode==2) tempRatioDataMCUnfold->GetYaxis()->SetRangeUser(0.8, 1.8);
            minYLines                               = 0.8;
            maxYLines                               = 1.2;
        }

        tempRatioDataMCUnfold->DrawCopy("e1");

        DrawGammaLines(0., tempRatioDataMCUnfold->GetXaxis()->GetBinUpEdge(tempRatioDataMCUnfold->GetNbinsX()), minYLines,  minYLines,  1.0, kGray+2, 8);
        DrawGammaLines(0., tempRatioDataMCUnfold->GetXaxis()->GetBinUpEdge(tempRatioDataMCUnfold->GetNbinsX()), 1.0,        1.0,        1.0, kGray+2, 2);
        DrawGammaLines(0., tempRatioDataMCUnfold->GetXaxis()->GetBinUpEdge(tempRatioDataMCUnfold->GetNbinsX()), maxYLines,  maxYLines,  1.0, kGray+2, 8);

        tempRatioDataMCUnfold->DrawCopy("e1,same");
        tempRatioDataMCBinwise->DrawCopy("e1,same");

    canvasGammaPurityConvBin->SaveAs(Form("%s/%s_%s_GammaSpectrum_SanityCheck_%s.%s",outputDir.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasGammaPurityConvBin;

    // is this the full correction factor that is applied to the raw spectrum?
    TH1D*       histoGammaCorrFac_Pt                = (TH1D*) histoGammaPrimaryRecoEff_Pt->Clone("GammaCorrFac_Pt");
    if (isPCM ) histoGammaCorrFac_Pt->Multiply(histoGammaConvProb_MCPt);

    //***********************************************************************************************************
    //******************************* Write output file for photons *********************************************
    //***********************************************************************************************************
    TString nameOutput = Form("%s/%s/%s_%s_GammaConvV1Correction_%s.root",cutSelection.Data(),energy.Data(),textPi0New.Data(),nameRec.Data(),cutSelection.Data());
    TFile* fileCorrectedOutput = new TFile(nameOutput,"RECREATE");

        //________________________ writing MC quantities to file _____________________________
        // input spectrum corrected for deta, dPhi, Nevt
        if (histoMCGammaSpec_MCPt)                              histoMCGammaSpec_MCPt->Write(                       "GammaSpecMC",                                  TObject::kOverwrite);
        if (histoSecondaryGammaFromXSpecPt[3])                  histoSecondaryGammaFromXSpecPt[3]->Write(           "histoSecondaryGammaFromXSpecPtRest",           TObject::kOverwrite);
        if (histoGammaTrueSecCocktailGammaFromX_Pt[3])          histoGammaTrueSecCocktailGammaFromX_Pt[3]->Write(   "histoGammaTrueSecCocktailGammaRest_Pt",        TObject::kOverwrite);
        // reconstructed MC gamma spectrum corrected deta, dPhi, Nevt
        if (histoMCrecGammaCorr_Pt)                             histoMCrecGammaCorr_Pt->Write(                      "GammaSpecCorrESDMC",                           TObject::kOverwrite);
        // split in different source (pi0,eta,...)
        for (Int_t i= 0; i< 8; i++)
            if (histoPhotonSource_MCPt[i])                      histoPhotonSource_MCPt[i]->Write(                   histoPhotonSource_MCPt[i]->GetName(),           TObject::kOverwrite);
        // pure direct photons
        if (histoPhotonSource_MCPt[8])                          histoPhotonSource_MCPt[8]->Write(                   "MC_DirectPhoton_Pt",                           TObject::kOverwrite);
        // input spectrum of converted photons corrected for deta, dphi, Nevt
        if (histoMCGammaConv_MCPt)                              histoMCGammaConv_MCPt->Write(                       "MC_ConvGamma_MCPt",                            TObject::kOverwrite);
        // reconstructed MC true photons corrected for deta, dphi, Nevt
        if (histoGammaTrueConv_Pt)                              histoGammaTrueConv_Pt->Write(                       "TrueConvGamma_Pt",                             TObject::kOverwrite);
        if (histoGammaTrueCalo_Pt)                              histoGammaTrueCalo_Pt->Write(                       "TrueCaloGamma_Pt",                             TObject::kOverwrite);
        // recontructed MC photon candidates, uncorrected
        if (histoMCrecGamma_Pt)                                 histoMCrecGamma_Pt->Write(                          "MCrec_ConvGamma_Pt",                           TObject::kOverwrite);
        if (histoMCrecGammaCalo_Pt)                             histoMCrecGammaCalo_Pt->Write(                      "MCrec_CaloGamma_Pt",                           TObject::kOverwrite);
        // summed reconstructed BG in MC
        if (histoMCrecBackground_Pt)                            histoMCrecBackground_Pt->Write(                     "MCrec_Background",                             TObject::kOverwrite);
        // combinatorial histos split
        if (histoCombinatorialSpecies_Pt)
            for (Int_t i = 0;i<17;i++)                          histoCombinatorialSpecies_Pt[i]->Write(             histoCombinatorialSpecies_Pt[i]->GetName(),     TObject::kOverwrite);
        if (histoCombinatorialSpeciesCalo_Pt)
            for (Int_t i=0; i<11; i++)                          histoCombinatorialSpeciesCalo_Pt[i]->Write(         histoCombinatorialSpeciesCalo_Pt[i]->GetName(), TObject::kOverwrite);


        //_________________________ writing secondary correction factors to file _________________________
        for (Int_t k = 0; k < 3; k++){
            if (histoGammaSecFromXRecoEff_MCPt[k])                  histoGammaSecFromXRecoEff_MCPt[k]->Write(               histoGammaSecFromXRecoEff_MCPt[k]->GetName(),
                                                                                                                            TObject::kOverwrite);
            if (histoGammaSecFromXRecoEff_MCPt_OrBin[k])            histoGammaSecFromXRecoEff_MCPt_OrBin[k]->Write(         histoGammaSecFromXRecoEff_MCPt_OrBin[k]->GetName(),
                                                                                                                            TObject::kOverwrite);
            if (histoGammaSecFromXRecoEff_RecPt[k])                 histoGammaSecFromXRecoEff_RecPt[k]->Write(              histoGammaSecFromXRecoEff_RecPt[k]->GetName(),
                                                                                                                            TObject::kOverwrite);
            if (histoGammaSecFromXRecoEff_RecPt_OrBin[k])           histoGammaSecFromXRecoEff_RecPt_OrBin[k]->Write(        histoGammaSecFromXRecoEff_RecPt_OrBin[k]->GetName(),
                                                                                                                            TObject::kOverwrite);
            if (histoGammaSecondaryFromXConvProb_MCPt[k])           histoGammaSecondaryFromXConvProb_MCPt[k]->Write(        histoGammaSecondaryFromXConvProb_MCPt[k]->GetName(),
                                                                                                                            TObject::kOverwrite);
            if (histoGammaSecondaryFromXConvProb_MCPt_OrBin[k])     histoGammaSecondaryFromXConvProb_MCPt_OrBin[k]->Write(  histoGammaSecondaryFromXConvProb_MCPt_OrBin[k]->GetName(),
                                                                                                                            TObject::kOverwrite);
            if (histoGammaSecGammaFromX_Cocktail_Raw_Pt[k])         histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]->Write(      histoGammaSecGammaFromX_Cocktail_Raw_Pt[k]->GetName(),
                                                                                                                            TObject::kOverwrite);
            if (histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin[k])   histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin[k]->Write(histoGammaSecGammaFromX_Cocktail_Raw_Pt_OrBin[k]->GetName(),
                                                                                                                            TObject::kOverwrite);
            if (histoGammaTrueSecCocktailGammaFromX_Pt[k])          histoGammaTrueSecCocktailGammaFromX_Pt[k]->Write(       Form("histoGammaTrueSecCocktailGammaFromXFrom%s_Pt",
                                                                                                                            nameSecondaries[k].Data()), TObject::kOverwrite);
            if (histoSecondaryGammaFromXSpecPt[k])                  histoSecondaryGammaFromXSpecPt[k]->Write(               Form("histoSecondaryGammaFromXFrom%sSpecPt", nameSecondaries[k].Data()),
                                                                                                                            TObject::kOverwrite);
            if (histoFracAllGammaToSecFromX_Cocktail_Pt[k])         histoFracAllGammaToSecFromX_Cocktail_Pt[k]->Write(      Form("histoSecCocktailGammaFromXFrom%sEffCorr",
                                                                                                                            nameSecondaries[k].Data()),TObject::kOverwrite);
            if (histoFracAllGammaToSecFromX_Pt[k])                  histoFracAllGammaToSecFromX_Pt[k]->Write(               Form("histoSecondaryGammaFromXFrom%sEffCorr", nameSecondaries[k].Data()),
                                                                                                                            TObject::kOverwrite);


        }
        if (histoFracAllGammaToSecFromX_Cocktail_Pt[3])         histoFracAllGammaToSecFromX_Cocktail_Pt[3]->Write( Form("histoSecCocktailGammaFromXFrom%sEffCorr", nameSecondaries[3].Data()),
                                                                                                                   TObject::kOverwrite);
        if (histoFracAllGammaToSecFromX_Pt[3])                  histoFracAllGammaToSecFromX_Pt[3]->Write(          Form("histoSecondaryGammaFromXFrom%sEffCorr", nameSecondaries[3].Data()),
                                                                                                                   TObject::kOverwrite);

        if (histoAllDiffDecayGamma)                             histoAllDiffDecayGamma->Write("histoRatioAllGammaDivDecayGammaSpectrumMC", TObject::kOverwrite);

        //_________________________ writing correction factors to file _________________________
        // photon purity without secondary subtraction
        if (histoGammaPurity_Pt)                                histoGammaPurity_Pt->Write(                 "GammaPurityWSec_Pt",                   TObject::kOverwrite);
        if (histoGammaCaloPurity_Pt)                            histoGammaCaloPurity_Pt->Write(             "histoGammaCaloPurityWSec_Pt",          TObject::kOverwrite);
        // photon purity with secondary subtraction
        if (histoGammaTruePurity_Pt)                            histoGammaTruePurity_Pt->Write(             "GammaPurityWOSec_Pt",                  TObject::kOverwrite);
        if (histoGammaCaloTruePurity_Pt)                        histoGammaCaloTruePurity_Pt->Write(         "GammaCaloPurityWOSec_Pt",              TObject::kOverwrite);
        // photon reconstruction efficiency including resolution correction
        if (histoGammaPrimaryRecoEff_Pt)                        histoGammaPrimaryRecoEff_Pt->Write(         "GammaRecoEff_WithResolCorr_Pt",        TObject::kOverwrite);
        if (histoGammaCaloPrimaryRecoEff_Pt)                    histoGammaCaloPrimaryRecoEff_Pt->Write(     "GammaCaloRecoEff_WithResolCorr_Pt",    TObject::kOverwrite);
        // photon reconstruction efficiency without resolution correction
        if (histoGammaPrimaryRecoEff_MCPt)                      histoGammaPrimaryRecoEff_MCPt->Write(       "GammaRecoEff_MCPt",                    TObject::kOverwrite);
        if (histoGammaCaloPrimaryRecoEff_MCPt)                  histoGammaCaloPrimaryRecoEff_MCPt->Write(   "GammaCaloRecoEff_MCPt",                TObject::kOverwrite);
        // photon conversion probability
        if (histoGammaConvProb_MCPt)                            histoGammaConvProb_MCPt->Write(             "GammaConvProb_Pt",                     TObject::kOverwrite);
        // resolution correction in case of unfolding
        if (histoGammaResolCorrUnfold_Pt){
            for (Int_t ipt = 1; ipt < histoGammaResolCorrUnfold_Pt->GetNbinsX()+1; ipt++){
                if (TMath::IsNaN(histoGammaResolCorrUnfold_Pt->GetBinContent(ipt)) || TMath::Finite(histoGammaResolCorrUnfold_Pt->GetBinContent(ipt))){
//                     cout << "needed correction" << endl;
                    histoGammaResolCorrUnfold_Pt->SetBinContent(ipt,-10000);
                }
            }
            histoGammaResolCorrUnfold_Pt->Write(                                                            "GammaResolCorrUnfold_Pt",              TObject::kOverwrite);
        }
        // photon correction factors (conv Prob, efficiency incl. resolution correction)
        if (histoGammaCorrFac_Pt)                               histoGammaCorrFac_Pt->Write(                "GammaCorrFac_Pt",                      TObject::kOverwrite);

        // ________________________ writing data quantities to file
        // uncorrected spectrum (scaled by 1/Nevt)
        if (histoGammaRawSpectrum_Pt)                           histoGammaRawSpectrum_Pt->Write(            "GammaRaw_Pt",                          TObject::kOverwrite);
        // corrected spectrum (legacy corrections: purity, effi incl resolution correction, secondaries, conv prob, plus trivial factors  )
        if (histoGammaCorrEffiReso_Pt)                          histoGammaCorrEffiReso_Pt->Write(           "GammaCorrEffiResol_Pt",                TObject::kOverwrite);
        // corrected spectrum (unfolding corrections: purity, secondaries, unfolding resolution correction, effi without unfolding corr, conv prob, plus trivial factors )
        if (histoGammaCorrUnfoldReso_PtCopy)                    histoGammaCorrUnfoldReso_PtCopy->Write(     "GammaCorrUnfold_PtControl",            TObject::kOverwrite);
        if (histoGammaCorrUnfoldReso_Pt)                        histoGammaCorrUnfoldReso_Pt->Write(         "GammaCorrUnfold_Pt",                   TObject::kOverwrite);
        if (histoGammaCorrUnfoldReso_BinByBin_Pt)               histoGammaCorrUnfoldReso_BinByBin_Pt->Write("GammaCorrUnfold_BinByBin_Pt",          TObject::kOverwrite);
        if(doPileUpCorr){
            // same as histoGammaCorrEffiReso_Pt with additional pileup correction
            if (histoGammaCorrEffiReso_PileUp_Pt)               histoGammaCorrEffiReso_PileUp_Pt->Write(            "GammaCorrEffiResolPileup_Pt",              TObject::kOverwrite);
            if (histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt)     histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt->Write(  "GammaCorrEffiResolPileup_NoMCUpdate_Pt",   TObject::kOverwrite);
            // pileup correction factor
            if (histoPileUpCorrectionFactor_Pt)                 histoPileUpCorrectionFactor_Pt->Write(              "PileUpCorrectionFactor",                   TObject::kOverwrite);
            if (histoPileUpCorrectionFactor_Pt_OrBin)           histoPileUpCorrectionFactor_Pt_OrBin->Write(        "PileUpCorrectionFactorOrBin",              TObject::kOverwrite);
            // ratio raw yield to raw yield after pileup subtraction
            if (histoRatioWithWithoutPileUp)                    histoRatioWithWithoutPileUp->Write(                 "RatioWithWithoutPileUp",                   TObject::kOverwrite);
            if (histoRatioWithWithoutPileUpFit)                 histoRatioWithWithoutPileUpFit->Write(              "RatioWithWithoutPileUpFit",                TObject::kOverwrite);
            if (graphGammaSysErrOOBPileupDown)                  graphGammaSysErrOOBPileupDown->Write(               "Gamma_OOBPileupSysDown",                   TObject::kOverwrite);
            if (graphGammaSysErrOOBPileupUp)                    graphGammaSysErrOOBPileupUp->Write(                 "Gamma_OOBPileupSysUp",                     TObject::kOverwrite);
        }
        // spectrum for pi0tagging (corrections: purity, secondaries)
        if (histoGammaRaw_PileUpPuritySecondaryCorrected_Pt)    histoGammaRaw_PileUpPuritySecondaryCorrected_Pt->Write("GammaPileUpPuritySecondaryCorr",        TObject::kOverwrite);

    fileCorrectedOutput->Close();
}
