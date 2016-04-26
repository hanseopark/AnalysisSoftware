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
//#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"

void CorrectGammaEffiResol(TH1D* histoGammaCorr, TH1D *histoAllSec, TH1D *histoK0sSec, TH1D* histoPurity, TH1D* histoConvProb, TH1D* histoRecoEff, Double_t deltaEta, Double_t scaling, Double_t nEvt){
    histoGammaCorr->Scale(1./nEvt);
    histoGammaCorr->Add(histoAllSec,-1);
    histoGammaCorr->Add(histoK0sSec,-1);
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
    histoGammaCorr->Divide(histoGammaCorr,histoConvProb,1.,1.,"");
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    histoGammaCorr->Scale(1./deltaEta);
    histoGammaCorr->Scale(scaling);
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoGammaCorr->GetBinContent(i)/histoGammaCorr->GetBinCenter(i);
        Double_t newBinError = histoGammaCorr->GetBinError(i)/histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,newBinError);
    }
}

void CorrectGammaSecAndPurity(TH1D* histoGammaCorr, TH1D* histoSecGamma, TH1D* histoSecGammaAddK0s, TH1D* histoPurity){
    histoGammaCorr->Sumw2();
    histoGammaCorr->Add(histoSecGamma,-1);
    histoGammaCorr->Add(histoSecGammaAddK0s,-1);
    histoGammaCorr->Multiply(histoGammaCorr,histoPurity,1.,1.,"");
}


void CorrectGammaUnfoldResol(TH1D* histoGammaCorr, TH1D* histoConvProb, TH1D* histoRecoEff, Double_t deltaEta, Double_t scaling, Double_t nEvt){
    histoGammaCorr->Divide(histoGammaCorr,histoConvProb,1.,1.,"");
    histoGammaCorr->Divide(histoGammaCorr,histoRecoEff,1.,1.,"");
    histoGammaCorr->Scale(1./deltaEta);
    histoGammaCorr->Scale(scaling);
    histoGammaCorr->Scale(1./nEvt);
    for (Int_t i = 1; i < histoGammaCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoGammaCorr->GetBinContent(i)/histoGammaCorr->GetBinCenter(i);
        Double_t newBinError = histoGammaCorr->GetBinError(i)/histoGammaCorr->GetBinCenter(i);
        histoGammaCorr->SetBinContent(i,newBinContent);
        histoGammaCorr->SetBinError(i,newBinError);
    }
}


void CorrectGammaMC(TH1D* histoMCGammaSpec_MCPtCorr,  Double_t deltaEta, Double_t scaling, Double_t nEvtMC){
    histoMCGammaSpec_MCPtCorr->Scale(1./deltaEta);
    histoMCGammaSpec_MCPtCorr->Scale(scaling);
    histoMCGammaSpec_MCPtCorr->Scale(1./nEvtMC);
    for (Int_t i = 1; i < histoMCGammaSpec_MCPtCorr->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoMCGammaSpec_MCPtCorr->GetBinContent(i)/histoMCGammaSpec_MCPtCorr->GetBinCenter(i);
        Double_t newBinError = histoMCGammaSpec_MCPtCorr->GetBinError(i)/histoMCGammaSpec_MCPtCorr->GetBinCenter(i);
        histoMCGammaSpec_MCPtCorr->SetBinContent(i,newBinContent);
        histoMCGammaSpec_MCPtCorr->SetBinError(i,newBinError);
    }
}


//*********************************************************************************************************
//***************************** Main routine to correct inclusive photons *********************************
//*********************************************************************************************************
void  CorrectGammaV2(   const char *nameUnCorrectedFile     = "myOutput", 
                        const char *nameCorrectionFile      = "", 
                        TString cutSelection                ="", 
                        TString suffix                      = "eps", 
                        TString nameMeson                   = "Pi0", 
                        TString isMC                        = "kFALSE",
                        TString option                      = "", 
                        TString optionPeriod                = "", 
                        TString fEstimatePileup             = "",//Bool_t kDoPileup = kFALSE, 
                        Int_t mode                          = 9
                    ){
    
    gROOT->Reset();
    
    //****************************************************************************************
    //*************************** Catch old versions *****************************************
    //****************************************************************************************
    if (mode == 1 || mode > 5){
        cout << "This macro isn't designed to run for Dalitz or the old software version" << endl;
        return;
    } else if ( mode == 4 || mode == 5){
        cout << "This macro can't yet deal with these modi" << endl;
        return;
    } else if ( mode == 2 || mode == 3){
        cout << "WARNING: running hybrid mode, this macro is still under construction for this mode" << endl;    
    }
    //*****************************************************************************************
    //************************** Set general style settings ***********************************
    //*****************************************************************************************
    StyleSettingsThesis();  
    SetPlotStyle();

    //*****************************************************************************************
    //*************************** Determine K0s scaling factor ********************************
    //*****************************************************************************************
    Double_t doubleAddFactorK0s = ReturnCorrectK0ScalingFactor( option,  cutSelection);
    if(isMC.CompareTo("kTRUE") == 0){
        doubleAddFactorK0s = 0.;
        cout<<" MC --> doubleAddFactorK0s = 0"<<endl;
    } else {
        cout << "The additional K0 correction factor is: "  << doubleAddFactorK0s<<endl;
    }

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
    TString centrality          = GetCentralityString(fEventCutSelection);
    TString collisionSystem     = ReturnFullCollisionsSystem(option);
    TString cent                = "";
    TString textMeasurement     = "#gamma";
    TString detectionProcess    = ReturnFullTextReconstructionProcess(mode);
    TString detectionProcess1   = "";
    TString detectionProcess2   = "";
    if (isPCM && isCalo){
        detectionProcess1       = ReturnFullTextReconstructionProcess(mode,1);
        detectionProcess        = detectionProcess1;
        detectionProcess2       = ReturnFullTextReconstructionProcess(mode,2);
    }
    if(option.Contains("PbPb")){
        cent        = Form("%s %s", centrality.Data(), collisionSystem.Data());
    } else {
        cent        = collisionSystem; 
    }    
    
    //******************************************************************************************
    //***************************** Set output directory ***************************************
    //******************************************************************************************
    TString outputDir           = Form("%s/%s/%s/CorrectGammaV2",cutSelection.Data(),option.Data(),suffix.Data());
    gSystem->Exec("mkdir "+outputDir);

    //******************************************************************************************
    //************************ Define additional global variables ******************************
    //******************************************************************************************
    TString textPi0New      = Form("Gamma_%s",nameMeson.Data());
    
    TString textPrefix2;
    if (isMC.CompareTo("kTRUE")==0){
        textPrefix2         = "MC";
    } else {
        textPrefix2         = "data";
    }
    
    Double_t deltaEta       = ReturnDeltaEta(fGammaCutSelection);
    Double_t eta            = deltaEta*0.5;
    Double_t deltaEtaCalo   = 0;
    Double_t deltaPhiCalo   = 0;
    Double_t minPhiCalo     = 0;
    Double_t maxPhiCalo     = 0;
    Double_t etaCalo        = 0;
    if (isPCM && isCalo){
        deltaEtaCalo   = ReturnDeltaEtaCalo(fClusterCutSelection);
        etaCalo        = deltaEtaCalo*0.5;
        deltaPhiCalo   = ReturnDeltaPhiCalo(fClusterCutSelection);
        TString phiMinCut(fClusterCutSelection(GetClusterPhiMinCutPosition(fClusterCutSelection),1));
        TString phiMaxCut(fClusterCutSelection(GetClusterPhiMaxCutPosition(fClusterCutSelection),1));
        minPhiCalo     = AnalyseClusterMinPhiCut(phiMinCut.Atoi());
        maxPhiCalo     = AnalyseClusterMaxPhiCut(phiMaxCut.Atoi());
    }
    
    Double_t scaling        = 1./(2.*TMath::Pi());
    Bool_t kDoPileup = kFALSE;
    if (fEstimatePileup.CompareTo("EstimateTrainPileUp") == 0)
        kDoPileup = kTRUE;
    
    // Read True Combinatorial Background
    TString combinatorics[17]                               = { "Elec+Elec","Elec+Pion","Elec+Kaon","Elec+Proton","Elec+Muon","Pion+Pion","Pion+Kaon","Pion+Proton",
                                                                "Pion+Muon","Kaon+Kaon","Kaon+Proton","Kaon+Muon","Proton+Proton","Proton+Muon","Muon+Muon","Rest","All"};

    Color_t colorsCombinatorics[17]                         = { kAzure-1, kRed+1, kOrange+7, kMagenta+2, kRed-9,
                                                                kBlue-6, kAzure+5, kPink, kCyan-3, kGreen-3,
                                                                kSpring+9, kGreen+2, kBlue+2, kMagenta-6, kSpring+4,
                                                                kCyan+2, 809
                                                              };
    Style_t markersCombinatorics[17]                        = { 20, 21, 24, 25, 27,
                                                                28, 29, 30, 33, 34, 
                                                                20, 21, 24, 25, 27,
                                                                28, 29
                                                              };
    TString decays[9]                                       = { "Pi0","Eta","Etap","Omega","Rho",
                                                                "Phi","Sigma", "All decays","direct #gamma"
                                                              };
    TString decaysLabels[9]                                  = { "#gamma from #pi^{0}","#gamma from #eta","#gamma from #eta'","#gamma from #omega","#gamma from #rho^{0}",
                                                                "#gamma from #phi", "#gamma from #Sigma^{0}", "All decays","direct #gamma"
                                                              };

    Color_t colorsDecay[9]                                  = { kRed+1, kAzure-1, kOrange+7, kGreen+2, kCyan-3, 809, 
                                                                kBlue+2, kBlack, kRed-1
                                                              };
    
    
    //****************************************************************************************** 
    //*******************File definitions and reading histograms from files ********************
    //******************************************************************************************
    // reading uncorrected quantities
    TString fileNameUncorrected             = nameUnCorrectedFile;
    TFile*  fileUnCorrected                 = new TFile(fileNameUncorrected.Data());
    TH1D*   histoESDConvGammaPt             = (TH1D*)fileUnCorrected->Get("ESD_ConvGamma_Pt");
    TH1D*   histoGamma_OriginalBin_Pt       = (TH1D*)fileUnCorrected->Get("ESD_ConvGamma_Pt_OriginalBinning");
    TH1D*   histoEventQuality               = (TH1D*)fileUnCorrected->Get("NEvents");

    // reading pileup correction factors
    TFile *doPileUpCorr;
    if (kDoPileup){ 
        doPileUpCorr                        = new  TFile(Form("%s/%s/%s_%s_GammaConvV1DCAHistogramms%s_%s.root", cutSelection.Data(), option.Data(), nameMeson.Data(),
                                                              textPrefix2.Data(), optionPeriod.Data(), cutSelection.Data()));
        if(doPileUpCorr->IsZombie()) 
            doPileUpCorr                    = 0;
    } else {
        doPileUpCorr                        = 0;
    }    
    
    TH1D*   histoESDConvGammaPtPileUp;
    TH1D*   histoPhotonDCAzFullPt;
    TH1D*   histoPileUpCorrectionFactor;
    if(doPileUpCorr){
        histoESDConvGammaPtPileUp           = (TH1D*)fileUnCorrected->Get("ESD_ConvGamma_Pt_PileUp");
        histoPhotonDCAzFullPt               = (TH1D*)fileUnCorrected->Get("ESD_GammaPtDCAzBin_Full");
        histoPileUpCorrectionFactor         = (TH1D*)histoESDConvGammaPtPileUp->Clone("PileUpCorrectionFactor");
        histoPileUpCorrectionFactor->Divide(histoPileUpCorrectionFactor,histoESDConvGammaPt,1,1,"B");
    }

    // Calculate number of events used for normalization
    Float_t nEvt;
    if (option.CompareTo("PbPb_2.76TeV") == 0){
        nEvt                                = histoEventQuality->GetBinContent(1);
    } else {
        nEvt                                = GetNEvents(histoEventQuality);
    }
    
    // Read binning histogram
    TH1D*   deltaPt                         = (TH1D*)fileUnCorrected->Get("deltaPt");
    for (Int_t i = 0; i < deltaPt->GetNbinsX() +1; i++){
        deltaPt->SetBinError(i, 0);
    }

    // Read correction histograms
    TFile*  fileCorrections                     = new TFile(nameCorrectionFile);
    TH1D*   histoGammaPurity_Pt                 = NULL;
    TH1D*   histoGammaTruePurity_Pt             = NULL;
    TH1D*   histoGammaTruePurity_OriginalBin_Pt = NULL;
    TH1D*   histoGammaConvProb_MCPt             = NULL;
    TH1D*   histoGammaPrimaryRecoEff_Pt         = NULL;
    TH1D*   histoGammaPrimaryRecoEff_MCPt       = NULL;
    
    if ( isPCM ) {
        // reconstructed validated photons / all reconstructed photon candidates in MC
        histoGammaPurity_Pt                 = (TH1D*)fileCorrections->Get("GammaPurity_Pt"); 
        // reconstructed validated primary photons / (all reconstructed photon candidates in MC - validated reconstructed secondary photons)
        // rescaling of histoGammaPurity_Pt for primary photons
        histoGammaTruePurity_Pt             = (TH1D*)fileCorrections->Get("GammaTruePurity_Pt"); 
        // reconstructed validated primary photons / (all reconstructed photon candidates in MC - validated reconstructed secondary photons)
        // rescaling of histoGammaPurity_Pt for primary photons
        // rebinned
        histoGammaTruePurity_OriginalBin_Pt = (TH1D*)fileCorrections->Get("GammaTruePurity_OriginalBinning_Pt"); 
        // converted photons MC/ all generated photons (always based on primaries only)
        histoGammaConvProb_MCPt             = (TH1D*)fileCorrections->Get("GammaConvProb_MCPt");
        // reconstructed validated primary photons vs rec pt / converted MC photons vs MCpt 
        // this reconstruction efficiency includes the transverse momentum resolution correction and the acceptance correction    
        histoGammaPrimaryRecoEff_Pt         = (TH1D*)fileCorrections->Get("GammaPrimaryRecoEff_Pt");
        // reconstructed validated primary photons vs MC pt / converted MC photons vs MCpt 
        // this reconstruction efficiency includes the acceptance correction        
        histoGammaPrimaryRecoEff_MCPt       = (TH1D*)fileCorrections->Get("GammaPrimaryRecoEff_MCPt");
    }
    // Pileup histograms from MC, same definition as above just corrected for MC pileup contribution
    TH1D*   histoGammaPurity_PileUp_Pt;
    TH1D*   histoGammaTruePurity_PileUp_Pt;
    TH1D*   histoGammaRecoEff_PileUp_Pt;
    TH1D*   histoGammaPrimaryRecoEff_PileUp_Pt;
    if( doPileUpCorr && isPCM ){
        histoGammaPurity_PileUp_Pt              = (TH1D*)fileCorrections->Get("GammaPurity_PileUp_Pt"); 
        histoGammaTruePurity_PileUp_Pt          = (TH1D*)fileCorrections->Get("GammaTruePurity_PileUp_Pt"); 
        histoGammaRecoEff_PileUp_Pt             = (TH1D*)fileCorrections->Get("GammaRecoEff_PileUp_Pt");
        histoGammaPrimaryRecoEff_PileUp_Pt      = (TH1D*)fileCorrections->Get("GammaPrimaryRecoEff_PileUp_Pt");
    }
        
    TH1D*   histoGammaCaloPurity_Pt             = NULL;
    TH1D*   histoGammaCaloTruePurity_Pt         = NULL;
    TH1D*   histoGammaCaloTruePurity_OrigBin_Pt = NULL;
    TH1D*   histoGammaCaloPrimaryRecoEff_Pt     = NULL;
    TH1D*   histoGammaCaloPrimaryRecoEff_MCPt   = NULL;
    // Calorimeter photons 
    if ( isPCM && isCalo ){
        // reconstructed validated photons / all reconstructed photon candidates in MC
        histoGammaCaloPurity_Pt             = (TH1D*)fileCorrections->Get("GammaCaloPurity_Pt"); 
        // reconstructed validated primary photons / (all reconstructed photon candidates in MC - validated reconstructed secondary photons)
        // rescaling of histoGammaPurity_Pt for primary photons
        histoGammaCaloTruePurity_Pt         = (TH1D*)fileCorrections->Get("GammaCaloTruePurity_Pt"); 
        // reconstructed validated primary photons / (all reconstructed photon candidates in MC - validated reconstructed secondary photons)
        // rescaling of histoGammaPurity_Pt for primary photons
        // rebinned
        histoGammaCaloTruePurity_OrigBin_Pt = (TH1D*)fileCorrections->Get("GammaCaloTruePurity_OriginalBinning_Pt"); 
        // reconstructed validated primary photons vs rec pt / converted MC photons vs MCpt 
        // this reconstruction efficiency includes the transverse momentum resolution correction and the acceptance correction    
        histoGammaCaloPrimaryRecoEff_Pt     = (TH1D*)fileCorrections->Get("GammaCaloPrimaryRecoEff_Pt");
        // reconstructed validated primary photons vs MC pt / converted MC photons vs MCpt 
        // this reconstruction efficiency includes the acceptance correction        
        histoGammaCaloPrimaryRecoEff_MCPt   = (TH1D*)fileCorrections->Get("GammaCaloPrimaryRecoEff_MCPt");
    
    }
    

    TH1D*   histoMCrecBackground_Pt                         = NULL;
    TH1D*   histoMCrecGamma_Pt                              = NULL;
    TH1D*   histoMCrecPrimaryGamma_Pt                       = NULL;
    TH1D*   histoMCAllGamma_MCPt                            = NULL;
    TH1D*   histoMCAllGamma_OriginalBin_MCPt                = NULL;
    TH1D*   histoMCGammaConv_MCPt                           = NULL;
    TH1D*   histoGammaTrueConv_Pt                           = NULL;
    TH1D*   histoGammaTruePrimaryConv_Pt                    = NULL;
    TH2D*   histoGammaTruePrimary_recPt_MCPt                = NULL;
    if (isPCM){
        // reconstructed background in MC
        histoMCrecBackground_Pt                             = (TH1D*)fileCorrections->Get("MCrec_Background");
        // reconstructed photon candidates in MC
        histoMCrecGamma_Pt                                  = (TH1D*)fileCorrections->Get("MCrec_ConvGamma_Pt");
        histoMCrecGamma_OriginalBin_Pt                      = (TH1D*)fileCorrections->Get("MCrec_ConvGamma_OriginalBinning_Pt");
        // reconstructed photon candidates in MC - validated secondary photons
        histoMCrecPrimaryGamma_Pt                           = (TH1D*)fileCorrections->Get("MCrec_PrimaryConvGamma_Pt");
        // MC input for all gamma, regardless of source
        histoMCAllGamma_MCPt                                = (TH1D*)fileCorrections->Get("MC_AllGamma_MCPt");
        histoMCAllGamma_OriginalBin_MCPt                    = (TH1D*)fileCorrections->Get("MC_AllGamma_OriginalBinning_MCPt");
        // MC converted photons regardless of source
        histoMCGammaConv_MCPt                               = (TH1D*)fileCorrections->Get("MC_ConvGamma_MCPt");
        // validated reconstructed photons
        histoGammaTrueConv_Pt                               = (TH1D*)fileCorrections->Get("TrueConvGamma_Pt");
        // validated reconstructed primary photons
        histoGammaTruePrimaryConv_Pt                        = (TH1D*)fileCorrections->Get("TruePrimaryConvGamma_Pt");
        // response matrix for photons recPt vs MC pt
        histoGammaTruePrimary_recPt_MCPt                    = (TH2D*)fileCorrections->Get("TruePrimaryConvGamma_recPt_MCPt");
        // response matrix for photons recPt vs MC pt rebinned
    }
    
    TH1D*   histoMCrecCaloBackground_Pt                     = NULL;
    TH1D*   histoMCrecGammaCalo_Pt                          = NULL;
    TH1D*   histoMCrecPrimaryGammaCalo_Pt                   = NULL;
    TH1D*   histoMCAllGammaCalo_MCPt                        = NULL;
    TH1D*   histoMCAllGammaCalo_OriginalBin_MCPt            = NULL;
    TH1D*   histoGammaCaloTrue_Pt                           = NULL;
    TH1D*   histoGammaCaloTruePrimary_Pt                    = NULL;
    TH2D*   histoGammaCaloTruePrimary_recPt_MCPt            = NULL;
    if (isPCM && isCalo){
        // reconstructed background in MC
        histoMCrecCaloBackground_Pt                         = (TH1D*)fileCorrections->Get("MCrec_Calo_Background");
        // reconstructed photon candidates in MC
        histoMCrecGammaCalo_Pt                              = (TH1D*)fileCorrections->Get("MCrec_CaloGamma_Pt");
        // reconstructed photon candidates in MC - validated secondary photons
        histoMCrecPrimaryGammaCalo_Pt                       = (TH1D*)fileCorrections->Get("MCrec_PrimaryCaloGamma_Pt");
        // MC input for all gamma, regardless of source
        histoMCAllGammaCalo_MCPt                            = (TH1D*)fileCorrections->Get("MC_AllGammaEMCAcc_MCPt");
        histoMCAllGammaCalo_OriginalBin_MCPt                = (TH1D*)fileCorrections->Get("MC_AllGammaEMCAcc_OriginalBinning_MCPt");
        // validated reconstructed photons
        histoGammaCaloTrue_Pt                               = (TH1D*)fileCorrections->Get("TrueCaloGamma_Pt");
        // validated reconstructed primary photons
        histoGammaCaloTruePrimary_Pt                        = (TH1D*)fileCorrections->Get("TruePrimaryCaloGamma_Pt");
        // response matrix for photons recPt vs MC pt
        histoGammaCaloTruePrimary_recPt_MCPt                = (TH2D*)fileCorrections->Get("TruePrimaryCaloGamma_recPt_MCPt");
    }
    
    
    // Pileup histograms, same definition as above just corrected for MC pileup contributions
    TH1D*   histoMCrecGamma_PileUp_Pt;
    TH1D*   histoPileUpCorrectionFactorMC;
    if(doPileUpCorr){
        histoMCrecGamma_PileUp_Pt                           = (TH1D*)fileCorrections->Get("MCrec_ConvGamma_Pt_PileUp");
        if (histoMCrecGamma_PileUp_Pt){
            histoPileUpCorrectionFactorMC                   = (TH1D*)histoMCrecGamma_PileUp_Pt->Clone("PileUpCorrectionFactorMC");
            histoPileUpCorrectionFactorMC->Divide(histoPileUpCorrectionFactorMC,histoMCrecGamma_Pt,1,1,"B");
        }
    }
    
    TH1F*   histoEventQualityMC                                 = (TH1F*)fileCorrections->Get("NEvents");

    TH1D**  histoPhotonSource_MCPt                              = new TH1D*[9];
    
    for(Int_t i = 0;i<7;i++){
        cout << __LINE__ << endl;
        histoPhotonSource_MCPt[i]                               = (TH1D*)fileCorrections->Get(Form("MC_DecayGamma%s_Pt",decays[i].Data()));
        histoPhotonSource_MCPt[i]->Sumw2();
        if (i == 0) histoPhotonSource_MCPt[7]                   = (TH1D*)histoPhotonSource_MCPt[0]->Clone("MC_DecayGammaAll_Pt");
        else histoPhotonSource_MCPt[7]->Add(histoPhotonSource_MCPt[i]);
    }
    cout << __LINE__ << endl;
    histoPhotonSource_MCPt[8]                                   = (TH1D*)histoMCAllGamma_OriginalBin_MCPt->Clone("MC_DirectPhotons");
    cout << __LINE__ << endl;
    histoPhotonSource_MCPt[8]->Sumw2();
    histoPhotonSource_MCPt[8]->Add(histoPhotonSource_MCPt[7],-1);
    
    TH1D**  histoCombinatorialSpecies_Pt                        = new TH1D*[17];
    for(Int_t i = 0;i<17;i++){
        histoCombinatorialSpecies_Pt[i]                         = (TH1D*)fileCorrections->Get(Form("ESD_TrueComb%s_Pt",combinatorics[i].Data()));
        histoCombinatorialSpecies_Pt[i]->SetMinimum(1e-10);
    }
    // Read secondary contamination histos
    TH1D*   histoGammaTrueSecConv_Pt                            = (TH1D*)fileCorrections->Get("TrueSecondaryConvGamma_Pt");
    TH1D*   histoGammaTrueSecConvGammaFromXFromK0s_Pt           = (TH1D*)fileCorrections->Get("TrueSecondaryConvGammaFromXFromK0s_Pt");
    TH1D*   histoGammaTrueSecConvGammaFromXFromLambda_Pt        = (TH1D*)fileCorrections->Get("TrueSecondaryConvGammaFromXFromLambda_Pt");
    TH1D*   histoFracAllGammaToSec_Pt                           = (TH1D*)fileCorrections->Get("FracAllGammaToSec");
    TH1D*   histoFracAllGammaToSecFromXFromK0s_Pt               = (TH1D*)fileCorrections->Get("FracAllGammaToSecFromXFromK0s");
    TH1D*   histoFracAllGammaToSecFromXFromLambda_Pt            = (TH1D*)fileCorrections->Get("FracAllGammaToSecFromXFromLambda");
    TH1D*   histoFracAllGammaToSec_OriginalBin_Pt               = (TH1D*)fileCorrections->Get("FracAllGammaToSecOriginalBinning");
    TH1D*   histoFracAllGammaToSecFromXFromK0s_OriginalBin_Pt   = (TH1D*)fileCorrections->Get("FracAllGammaToSecFromXFromK0sOriginalBinning");
    TH1D*   histoFracAllGammaToSecFromXFromLambda_OriginalBin_Pt= (TH1D*)fileCorrections->Get("FracAllGammaToSecFromXFromLambdaOriginalBinning");

    TH1D*   histoFracAllGammaToSec_PileUp_Pt;
    TH1D*   histoFracAllGammaToSecFromXFromK0s_PileUp_Pt;
    TH1D*   histoMCrecPhotonDCAzFullPt;
    if(doPileUpCorr){
        histoFracAllGammaToSec_PileUp_Pt                        = (TH1D*)fileCorrections->Get("FracAllGammaToSecPileUp");
        histoFracAllGammaToSecFromXFromK0s_PileUp_Pt            = (TH1D*)fileCorrections->Get("FracAllGammaToSecFromXFromK0sPileUp");
        histoMCrecPhotonDCAzFullPt                              = (TH1D*)fileCorrections->Get("MCrec_GammaPtDCAzBin_Full");     // category must be adapted
    }
    // Determine number of events in MC
    Float_t nEvtMC;
    if (option.CompareTo("PbPb_2.76TeV") == 0){
        nEvtMC = histoEventQualityMC->GetBinContent(1);
    } else {
        nEvtMC = GetNEvents(histoEventQualityMC);
    }

    // Determine maximum pT
    Double_t maxPtGamma                                     = histoESDConvGammaPt->GetXaxis()->GetBinUpEdge(histoESDConvGammaPt->GetNbinsX());
    
    
    // Proper Scaling Background
    TH1D *ScalingGammaBackground_Pt                         = (TH1D*) histoESDConvGammaPt->Clone("ScalingGammaBackground_Pt");
    ScalingGammaBackground_Pt->Divide(ScalingGammaBackground_Pt, histoMCrecGamma_Pt, 1., 1, "");
    ScalingGammaBackground_Pt->Scale(nEvtMC/nEvt);
    histoMCrecBackground_Pt->Scale(1./nEvtMC);
    TH1D *histoGammaMCBackground_Pt                         = (TH1D*)histoMCrecBackground_Pt->Clone("histoGammaMCBackground_Pt");
    histoMCrecBackground_Pt->Multiply(ScalingGammaBackground_Pt);

    //**********************************************************************************
    //******************** PrimVtx DCA Plot ********************************************
    //**********************************************************************************

    if( doPileUpCorr && isPCM ){      
        TCanvas *canvasPileUpCorrFactor = GetAndSetCanvas("canvasPileUpCorrFactor");

            DrawGammaSetMarker(histoPileUpCorrectionFactor, 20, 3, 1, 1);
            if (histoMCrecGamma_PileUp_Pt)DrawGammaSetMarker(histoPileUpCorrectionFactorMC, 24, 3, 2, 2);

            SetHistogramm(histoPileUpCorrectionFactor,"#it{p}_{T} (GeV/#it{c})","Correction Factor (%)",0.84,1.02);
            if (histoMCrecGamma_PileUp_Pt)SetHistogramm(histoPileUpCorrectionFactorMC,"#it{p}_{T} (GeV/#it{c})","Correction Factor (%)",0.84,1.02);

            histoPileUpCorrectionFactor->DrawCopy("");
            if (histoMCrecGamma_PileUp_Pt)histoPileUpCorrectionFactorMC->Draw("same");

            TLegend* legendPileUpCorrFactor = GetAndSetLegend(0.3,0.2,2.2,1,cent);
            legendPileUpCorrFactor->AddEntry(histoPileUpCorrectionFactor,"Correction Factor Data","lp");
            if (histoMCrecGamma_PileUp_Pt)legendPileUpCorrFactor->AddEntry(histoPileUpCorrectionFactorMC,"Correction Factor MC","lp");
            legendPileUpCorrFactor->Draw();
        
        canvasPileUpCorrFactor->SaveAs(Form("%s/%s_PileUpCorrFactor_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasPileUpCorrFactor;
    }
    
    //**********************************************************************************
    //******************** Background Plot *********************************************
    //**********************************************************************************
    TCanvas *canvasBackground = GetAndSetCanvas("canvasBackground");
    canvasBackground->SetLogy();

        TH1D* histoGammaRawSpectrum_PtMC = (TH1D*) histoMCrecGamma_Pt->Clone("histoGammaRawSpectrum_PtMC");
        histoGammaRawSpectrum_PtMC->Scale(1./nEvtMC);
        TH1D* histoGammaRawSpectrum_Pt = (TH1D*) histoESDConvGammaPt->Clone("histoGammaRawSpectrum_Pt");
        histoGammaRawSpectrum_Pt->Scale(1./nEvt);

        DrawGammaSetMarker(histoGammaRawSpectrum_Pt, 20, 1.0, 1, 1);
        DrawGammaSetMarker(histoMCrecBackground_Pt, 20, 1.0, 2, 2);

        SetHistogramm(histoMCrecBackground_Pt,"#it{p}_{T} (GeV/#it{c})",Form("Background for #gamma in |#eta| < %g",eta));
        SetHistogramm(histoGammaRawSpectrum_Pt,"#it{p}_{T} (GeV/#it{c})","Raw yield",1e-8,5);

        histoGammaRawSpectrum_Pt->DrawCopy("");
        histoMCrecBackground_Pt->Draw("same");
        
        TLegend* legendBackground = GetAndSetLegend(0.6,0.7,2,1);
        legendBackground->AddEntry(histoGammaRawSpectrum_Pt,"raw #gamma spectrum","pl");
        legendBackground->AddEntry(histoMCrecBackground_Pt,"#gamma background","pl");
        legendBackground->Draw();

        // labeling
        PutProcessLabelAndEnergyOnPlot( 0.6, 0.95, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
    canvasBackground->SaveAs(Form("%s/%s_Background_%s_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasBackground;
    
    //****************************************************************************************** 
    //************************ Calculating background from Secondaries *************************
    //****************************************************************************************** 
    TH1D *histoSecondaryGammaSpecPt                 = (TH1D*) histoESDConvGammaPt->Clone("SecondaryGammaSpecPt");
    TH1D *histoSecondaryGammaFromXFromK0sSpecPt     = (TH1D*) histoESDConvGammaPt->Clone("SecondaryGammaSpecFromXFromK0sPt");
    histoSecondaryGammaSpecPt->Multiply(histoFracAllGammaToSec_Pt);
    histoSecondaryGammaFromXFromK0sSpecPt->Multiply(histoFracAllGammaToSecFromXFromK0s_Pt);
    
    histoGammaTrueSecConv_Pt->Scale(1./nEvtMC);
    histoGammaTrueSecConvGammaFromXFromK0s_Pt->Scale(1./nEvtMC);
    histoSecondaryGammaSpecPt->Scale(1./nEvt);
    histoSecondaryGammaFromXFromK0sSpecPt->Scale(1./nEvt);
    histoSecondaryGammaFromXFromK0sSpecPt->Scale(doubleAddFactorK0s);

    // Correct the secondary fractions for pileup
    TH1D *histoSecondaryGammaSpecPtPileUp;
    TH1D *histoSecondaryGammaFromXFromK0sSpecPtPileUp;
    if(doPileUpCorr){
        histoSecondaryGammaSpecPtPileUp             = (TH1D*) histoESDConvGammaPtPileUp->Clone("SecondaryGammaSpecPtPileUp");
        histoSecondaryGammaFromXFromK0sSpecPtPileUp = (TH1D*) histoESDConvGammaPtPileUp->Clone("SecondaryGammaSpecFromXFromK0sPtPileUp");
        histoSecondaryGammaSpecPtPileUp->Multiply(histoFracAllGammaToSec_PileUp_Pt);
        histoSecondaryGammaFromXFromK0sSpecPtPileUp->Multiply(histoFracAllGammaToSecFromXFromK0s_PileUp_Pt);
        histoSecondaryGammaSpecPtPileUp->Scale(1./nEvt);
        histoSecondaryGammaFromXFromK0sSpecPtPileUp->Scale(1./nEvt);
        histoSecondaryGammaFromXFromK0sSpecPtPileUp->Scale(doubleAddFactorK0s);
    }
    
    // plotting secondary spectra
    TCanvas *canvasSecSpec = GetAndSetCanvas("canvasSecSpec");
    canvasSecSpec->SetLogy();

        SetHistogramm(histoSecondaryGammaSpecPt, "#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
        SetHistogramm(histoSecondaryGammaFromXFromK0sSpecPt,"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
        SetHistogramm(histoGammaTrueSecConvGammaFromXFromK0s_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
        SetHistogramm(histoGammaTrueSecConvGammaFromXFromLambda_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma");
        SetHistogramm(histoGammaTrueSecConv_Pt,"#it{p}_{T} (GeV/#it{c})","Secondary Converted #gamma",1e-10,1e-1);

        DrawGammaSetMarker(histoGammaTrueSecConv_Pt, 20, 1.0, kRed+1, kRed+1);
        DrawGammaSetMarker(histoGammaTrueSecConvGammaFromXFromK0s_Pt, 24, 1.0, kRed-6, kRed-6);
        DrawGammaSetMarker(histoGammaTrueSecConvGammaFromXFromLambda_Pt, 21, 1.0, kRed-6, kRed-6);
        DrawGammaSetMarker(histoSecondaryGammaSpecPt, 21, 1.0, kBlue+1, kBlue+1);
        DrawGammaSetMarker(histoSecondaryGammaFromXFromK0sSpecPt, 25, 1.0, kBlue-6, kBlue-6);

        histoGammaTrueSecConv_Pt->Draw();
        histoGammaTrueSecConvGammaFromXFromK0s_Pt->Draw("same");
        histoGammaTrueSecConvGammaFromXFromLambda_Pt->Draw("same");
        histoSecondaryGammaSpecPt->Draw("same");
        histoSecondaryGammaFromXFromK0sSpecPt->DrawCopy("same");

        TLegend* legendSecSpec = GetAndSetLegend(0.45,0.7,5);
        legendSecSpec->AddEntry(histoGammaTrueSecConv_Pt,"Raw MC Sec #gamma Spectrum","pl");
        legendSecSpec->AddEntry(histoGammaTrueSecConvGammaFromXFromK0s_Pt,"Raw MC Sec #gamma from X from K^{0}_{s}","pl");
        legendSecSpec->AddEntry(histoGammaTrueSecConvGammaFromXFromLambda_Pt,"Raw MC Sec #gamma from X from #Lambda","pl");
        legendSecSpec->AddEntry(histoSecondaryGammaSpecPt,"Raw data Sec #gamma Spectrum ","pl");
        legendSecSpec->AddEntry(histoSecondaryGammaFromXFromK0sSpecPt,"Raw data Sec #gamma from X from K^{0}_{s}","pl");
        legendSecSpec->AddEntry((TObject*)0, "(scaled with K^{0}_{s} factor)","");
        legendSecSpec->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.18, 0.3, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
    canvasSecSpec->SaveAs(Form("%s/%s_SecondarySpectra_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasSecSpec;
    
    // plotting secondary fractions
    TCanvas *canvasSecFrac = GetAndSetCanvas("canvasSecFrac");
    canvasSecFrac->SetTopMargin(0.035);

        SetHistogramm(histoFracAllGammaToSecFromXFromK0s_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary Converted #gamma");
        SetHistogramm(histoFracAllGammaToSec_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary Converted #gamma");
        SetHistogramm(histoFracAllGammaToSecFromXFromLambda_Pt,"#it{p}_{T} (GeV/#it{c})","Fraction of Secondary Converted #gamma");
    
        DrawGammaSetMarker(histoFracAllGammaToSec_Pt, 20, 1.0, 1, 1);
        DrawGammaSetMarker(histoFracAllGammaToSecFromXFromK0s_Pt, 24, 1.0, kBlue-6, kBlue-6);
        DrawGammaSetMarker(histoFracAllGammaToSecFromXFromLambda_Pt, 21, 1.0, kRed-6, kRed-6);

        histoFracAllGammaToSec_Pt->Draw();
        histoFracAllGammaToSecFromXFromK0s_Pt->Draw("same");
        histoFracAllGammaToSecFromXFromLambda_Pt->Draw("same");

        TLegend* legendSecFrac = GetAndSetLegend(0.45,0.8,3,1);
        legendSecFrac->AddEntry(histoFracAllGammaToSec_Pt,"Fraction Sec #gamma","pl");
        legendSecFrac->AddEntry(histoFracAllGammaToSecFromXFromK0s_Pt,Form("Fraction Sec #gamma from X from K^{0}_{s}"),"pl");
        legendSecFrac->AddEntry(histoFracAllGammaToSecFromXFromLambda_Pt,Form("Fraction Sec #gamma from X from #Lambda"),"pl");
        legendSecFrac->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.6, 0.75, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
    canvasSecFrac->SaveAs(Form("%s/%s_SecondaryGammaFraction_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasSecFrac;

    //**********************************************************************************
    //******************** Purity Plot *************************************************
    //**********************************************************************************
    // Plain Purity plotted
    TCanvas *canvasPurity = GetAndSetCanvas("canvasPurity");

        if ( isPCM ){
            DrawGammaSetMarker(histoGammaPurity_Pt, 24, 1.0, kRed+2, kRed+2);
            DrawGammaSetMarker(histoGammaTruePurity_Pt, 20, 1.0, 1, 1);
            if(doPileUpCorr)DrawGammaSetMarker(histoGammaTruePurity_PileUp_Pt, 20, 1.0, kBlue+2, kBlue+2);
            
            SetHistogramm(histoGammaTruePurity_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{pur,#gamma} in |#eta| < %g",eta),0.6, 1.);
            histoGammaTruePurity_Pt->Draw();
            histoGammaPurity_Pt->Draw("same");
            if(doPileUpCorr)histoGammaTruePurity_PileUp_Pt->Draw("same");
            
            TLegend* legendPurity;
            if (doPileUpCorr)
                legendPurity = GetAndSetLegend(0.18,0.15,3);
            else
                legendPurity = GetAndSetLegend(0.18,0.15,2);
            legendPurity->AddEntry(histoGammaPurity_Pt,"purity");
            legendPurity->AddEntry(histoGammaTruePurity_Pt,"rescaled purity, secondaries removed");
            if(doPileUpCorr)legendPurity->AddEntry(histoGammaTruePurity_PileUp_Pt,"rescaled purity, secondaries removed, pileup");
            legendPurity->Draw();

            PutProcessLabelAndEnergyOnPlot( 0.18, 0.5, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
            canvasPurity->SaveAs(Form("%s/%s_Purity_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        } 
        
        if ( isPCM && isCalo ){
            DrawGammaSetMarker(histoGammaCaloPurity_Pt, 24, 1.0, kRed+2, kRed+2);
            DrawGammaSetMarker(histoGammaCaloTruePurity_Pt, 20, 1.0, 1, 1);
            
            SetHistogramm(histoGammaCaloTruePurity_Pt,"#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{pur,#gamma} in |#eta| < %g",etaCalo),0.5, 1.);
            histoGammaCaloTruePurity_Pt->Draw();
            histoGammaCaloPurity_Pt->Draw("same");
            
            TLegend* legendPurity = GetAndSetLegend(0.45,0.15,2);
            legendPurity->AddEntry(histoGammaPurity_Pt,"purity");
            legendPurity->AddEntry(histoGammaTruePurity_Pt,"rescaled purity, secondaries removed");
            legendPurity->Draw();

            PutProcessLabelAndEnergyOnPlot( 0.7, 0.4, 0.035, cent, textMeasurement, detectionProcess2, 42, 0.03);
        
            canvasPurity->SaveAs(Form("%s/%s_PurityCalo_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        }
    delete canvasPurity;

    if ( isPCM ){
        //**********************************************************************************
        //******************** Conversion Prob Plot ****************************************
        //**********************************************************************************
        TCanvas *canvasConvProb = GetAndSetCanvas("canvasConvProb");
            DrawGammaSetMarker(histoGammaConvProb_MCPt, 20, 1.0, 1, 1);
            SetHistogramm(histoGammaConvProb_MCPt, "#it{p}_{T} (GeV/#it{c})",Form("#it{P}_{conv} in |#eta| < %g",eta), 0.04, 0.10);
            histoGammaConvProb_MCPt->Draw();

            TF1 *fConv  = new TF1("line","[0]",2.5,25.);
            histoGammaConvProb_MCPt->Fit(fConv,"QRME0");
            Double_t parameterProb[1];
            fConv->GetParameters(parameterProb);

            DrawGammaLines(0., maxPtGamma,parameterProb[0], parameterProb[0],1);

            PutProcessLabelAndEnergyOnPlot( 0.75, 0.3, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);

        canvasConvProb->SaveAs(Form("%s/%s_ConversionProb_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        delete canvasConvProb;
    }

    //**********************************************************************************
    //******************** Reconstructin Eff Plot **************************************
    //**********************************************************************************
    // drawing reco Efficiency
    TCanvas *canvasRecoEff = GetAndSetCanvas("canvasRecoEff");
        if ( isPCM ){
            DrawGammaSetMarker(histoGammaPrimaryRecoEff_MCPt, 20, 1.0, 1, 1);
            SetHistogramm(histoGammaPrimaryRecoEff_MCPt,"#it{p}_{T,MC} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);
            histoGammaPrimaryRecoEff_MCPt->Draw();
        
            PutProcessLabelAndEnergyOnPlot( 0.75, 0.95, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
            canvasRecoEff->SaveAs(Form("%s/%s_ReconstructionEff_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        }
        if ( isPCM && isCalo ){
            DrawGammaSetMarker(histoGammaCaloPrimaryRecoEff_MCPt, 20, 1.0, 1, 1);
            SetHistogramm(histoGammaCaloPrimaryRecoEff_MCPt,"#it{p}_{T,MC} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",etaCalo), 0., 1.0);
            histoGammaCaloPrimaryRecoEff_MCPt->Draw();
        
            PutProcessLabelAndEnergyOnPlot( 0.75, 0.95, 0.035, cent, textMeasurement, detectionProcess2, 42, 0.03);
        
            canvasRecoEff->SaveAs(Form("%s/%s_ReconstructionEffCalo_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        
        }    
    delete canvasRecoEff;
    
    // draw comparison of reco efficiencies
    TCanvas* canvasCompRecoEff = GetAndSetCanvas("canvasCompRecoEff");
    DrawGammaCanvasSettings( canvasCompRecoEff,0.1, 0.015, 0.02, 0.09);
    TPad* padCompRecoEff = new TPad("padCompRecoEff", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padCompRecoEff, 0.10, 0.015, 0.02, 0.);
    padCompRecoEff->Draw();
    TPad* padBinCompRecoEffRatio = new TPad("padBinCompRecoEffRatio", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings( padBinCompRecoEffRatio,  0.1, 0.015, 0.0, 0.3);
    padBinCompRecoEffRatio->Draw();

    if (isPCM){
        padCompRecoEff->cd();
        SetHistogramm(histoGammaPrimaryRecoEff_MCPt, "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);
        SetHistogramm(histoGammaPrimaryRecoEff_Pt, "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0.05, 1.0);
        if(doPileUpCorr)SetHistogramm(histoGammaPrimaryRecoEff_PileUp_Pt, "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);

        DrawGammaSetMarker(histoGammaPrimaryRecoEff_MCPt, 20, 1.0, kBlack, kBlack);
        DrawGammaSetMarker(histoGammaPrimaryRecoEff_Pt, 24, 1.0, kBlue+1, kBlue+1);
        if(doPileUpCorr)DrawGammaSetMarker(histoGammaPrimaryRecoEff_PileUp_Pt, 20, 1.0, 5, 5);

        histoGammaPrimaryRecoEff_Pt->Draw("e1");
        if(doPileUpCorr)histoGammaPrimaryRecoEff_PileUp_Pt->Draw("e1,same");
        histoGammaPrimaryRecoEff_MCPt->Draw("same,e1");

        TLegend* legendCompRecoEff;
        if(doPileUpCorr)
            legendCompRecoEff = GetAndSetLegend(0.15 ,0.05,3,1);
        else
            legendCompRecoEff = GetAndSetLegend(0.15 ,0.05,2,1);
        legendCompRecoEff->SetTextSize(0.04);
        legendCompRecoEff->AddEntry(histoGammaPrimaryRecoEff_MCPt,"MC pT Reconstruction Efficiency (unfolding)");
        legendCompRecoEff->AddEntry(histoGammaPrimaryRecoEff_Pt,"Reconstruction Efficiency for primary  #gamma");
        if(doPileUpCorr)legendCompRecoEff->AddEntry(histoGammaPrimaryRecoEff_PileUp_Pt,"Reconstruction Efficiency for primary  #gamma Pile Up");
        legendCompRecoEff->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.75, 0.95, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
    }

    TH1D* histoGammaResolCorrEff_Pt = NULL;
    if (isPCM){
        histoGammaResolCorrEff_Pt = (TH1D*)histoGammaPrimaryRecoEff_Pt->Clone("histoGammaResolCorrEff_Pt");
        histoGammaResolCorrEff_Pt->Divide(histoGammaPrimaryRecoEff_MCPt);

        padBinCompRecoEffRatio->cd();
        TH1D* histoGammaPrimaryRecoEff_PileUp_PtRatio;
        if(doPileUpCorr){
            histoGammaPrimaryRecoEff_PileUp_PtRatio = (TH1D*)histoGammaPrimaryRecoEff_PileUp_Pt->Clone("GammaPrimaryRecoEffRatioPileUp");
            histoGammaPrimaryRecoEff_PileUp_PtRatio->Divide(histoGammaPrimaryRecoEff_MCPt);
        }
        SetStyleHistoTH1ForGraphs(  histoGammaResolCorrEff_Pt, "#it{p}_{T} (GeV/#it{c})","#frac{Reco Eff x}{Primary MC Reco Eff}" , 0.14, 0.15,
                                    0.12, 0.10,  0.85, 0.4, 510, 505);
        histoGammaResolCorrEff_Pt->GetYaxis()->SetRangeUser(0.7,1.25);

        histoGammaResolCorrEff_Pt->DrawCopy("e1");
        if(doPileUpCorr)histoGammaPrimaryRecoEff_PileUp_PtRatio->DrawCopy("e1,same");
        
        DrawGammaLines(0., maxPtGamma,1, 1,0.5);
    
        canvasCompRecoEff ->SaveAs(Form("%s/%s_CompRecEff_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    }
    if (isPCM && isCalo){
        padCompRecoEff->cd();
        SetHistogramm(histoGammaCaloPrimaryRecoEff_MCPt, "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0., 1.0);
        SetHistogramm(histoGammaCaloPrimaryRecoEff_Pt, "#it{p}_{T} (GeV/#it{c})",Form("#epsilon_{eff,#gamma} in |#eta| < %g",eta), 0.05, 1.0);
        
        DrawGammaSetMarker(histoGammaCaloPrimaryRecoEff_MCPt, 20, 1.0, kBlack, kBlack);
        DrawGammaSetMarker(histoGammaCaloPrimaryRecoEff_Pt, 24, 1.0, kBlue+1, kBlue+1);
        
        histoGammaCaloPrimaryRecoEff_Pt->Draw("e1");
        histoGammaCaloPrimaryRecoEff_MCPt->Draw("same,e1");

        TLegend* legendCompRecoEff = GetAndSetLegend(0.15 ,0.05,3,1);
        legendCompRecoEff->SetTextSize(0.04);
        legendCompRecoEff->AddEntry(histoGammaCaloPrimaryRecoEff_MCPt,"MC pT Reconstruction Efficiency (unfolding)");
        legendCompRecoEff->AddEntry(histoGammaCaloPrimaryRecoEff_Pt,"Reconstruction Efficiency for primary  #gamma");
        legendCompRecoEff->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.75, 0.95, 0.035, cent, textMeasurement, detectionProcess2, 42, 0.03);
    }

    TH1D* histoGammaCaloResolCorrEff_Pt = NULL;
    if (isPCM && isCalo){
        histoGammaCaloResolCorrEff_Pt = (TH1D*)histoGammaCaloPrimaryRecoEff_Pt->Clone("histoGammaCaloResolCorrEff_Pt");
        histoGammaCaloResolCorrEff_Pt->Divide(histoGammaCaloPrimaryRecoEff_MCPt);

        padBinCompRecoEffRatio->cd();
        SetStyleHistoTH1ForGraphs(  histoGammaCaloResolCorrEff_Pt, "#it{p}_{T} (GeV/#it{c})","#frac{Reco Eff x}{Primary MC Reco Eff}" , 0.14, 0.15,  
                                    0.12, 0.10,  0.85, 0.4, 510, 505);
        histoGammaCaloResolCorrEff_Pt->GetYaxis()->SetRangeUser(0.7,1.25);

        histoGammaCaloResolCorrEff_Pt->DrawCopy("e1");
        
        DrawGammaLines(0., maxPtGamma,1, 1,0.5);
    
        canvasCompRecoEff ->SaveAs(Form("%s/%s_CompCaloRecEff_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    }
        
    delete padCompRecoEff;
    delete padBinCompRecoEffRatio;
    delete canvasCompRecoEff;
    
    //**********************************************************************************
    //******************** Response Matrix for detector resolution *********************
    //**********************************************************************************
    if (isPCM){
        SetStyleHistoTH2ForGraphs(  histoGammaTruePrimary_recPt_MCPt, "Reconstructed #it{p}_{T} (GeV/#it{c})","MC #it{p}_{T} (GeV/#it{c})", 0.035, 0.04,  
                                    0.035, 0.04, 0.9, 1.0, 510, 510);
        histoGammaTruePrimary_recPt_MCPt->GetYaxis()->SetRangeUser(0,25);
        histoGammaTruePrimary_recPt_MCPt->GetXaxis()->SetRangeUser(0,25);
    }
    if (isPCM && isCalo){
        SetStyleHistoTH2ForGraphs(  histoGammaCaloTruePrimary_recPt_MCPt, "Reconstructed #it{p}_{T} (GeV/#it{c})","MC #it{p}_{T} (GeV/#it{c})", 0.035, 0.04,  
                                    0.035, 0.04, 0.9, 1.0, 510, 510);
        histoGammaCaloTruePrimary_recPt_MCPt->GetYaxis()->SetRangeUser(0,25);
        histoGammaCaloTruePrimary_recPt_MCPt->GetXaxis()->SetRangeUser(0,25);
        
    }
    // Plotting response matrix
    TCanvas * canvasResponseMatrix = new TCanvas("canvasResponseMatrix","",480,440);  // gives the page size
    DrawGammaCanvasSettings( canvasResponseMatrix, 0.09, 0.105, 0.02, 0.085);    
    canvasResponseMatrix->SetLogz(1);
    canvasResponseMatrix->cd();

        if (isPCM){
            histoGammaTruePrimary_recPt_MCPt->Draw("colz");
            PutProcessLabelAndEnergyOnPlot( 0.18, 0.95, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
            canvasResponseMatrix->SaveAs(Form("%s/%s_ResponseMatrix_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
        } 
        
        if (isPCM && isCalo){
            histoGammaCaloTruePrimary_recPt_MCPt->Draw("colz");
            PutProcessLabelAndEnergyOnPlot( 0.18, 0.95, 0.035, cent, textMeasurement, detectionProcess2, 42, 0.03);
            canvasResponseMatrix->SaveAs(Form("%s/%s_ResponseMatrixCalo_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
            
        }    
    
//     delete canvasResponseMatrix;
    
    //**********************************************************************************
    //************************ Unfolding of inclusive gamma spectrum *******************
    //**********************************************************************************
    
    TH1D* histoGammaCorrUnfoldReso_Pt                   = NULL;
    TH1D* histoGammaCorrUnfoldReso_BinByBin_Pt          = NULL;
    TH1D* histoGammaCorrUnfoldResoPileUp_Pt             = NULL;
    TH1D* histoGammaResolCorrUnfold_Pt                  = NULL; 
    TH1D* histoGammaCorrUnfoldReso_PtNotCorrected       = NULL;
    TH1D* histoGammaResolCorrUnfold_BinByBin_Pt         = NULL;
    
    // do the same with MC rec gammas as a sanity check
    TH1D* histoMCrecGammaCorr_Pt                        = NULL;

    if (isPCM){
        // determine secondary contribution in orginal binning
        TH1D *histoSecondaryGammaSpecPtOriginalBin             = (TH1D*) histoGamma_OriginalBin_Pt->Clone("SecondaryGammaSpecPt");
        TH1D *histoSecondaryGammaFromXFromK0sSpecPtOriginalBin = (TH1D*) histoGamma_OriginalBin_Pt->Clone("SecondaryGammaSpecFromXFromK0sPt");
        histoSecondaryGammaSpecPtOriginalBin->Multiply(histoFracAllGammaToSec_OriginalBin_Pt);
        histoSecondaryGammaFromXFromK0sSpecPtOriginalBin->Multiply(histoFracAllGammaToSecFromXFromK0s_OriginalBin_Pt);
        histoSecondaryGammaFromXFromK0sSpecPtOriginalBin->Scale(doubleAddFactorK0s);

        //subtract secondary contribution from inclusive spectrum
        CorrectGammaSecAndPurity(histoGamma_OriginalBin_Pt, histoSecondaryGammaSpecPtOriginalBin, histoSecondaryGammaFromXFromK0sSpecPtOriginalBin, histoGammaTruePurity_OriginalBin_Pt );
        CorrectGammaSecAndPurity(histoMCrecGamma_OriginalBin_Pt, histoSecondaryGammaSpecPtOriginalBin, histoSecondaryGammaFromXFromK0sSpecPtOriginalBin, histoGammaTruePurity_OriginalBin_Pt );
        
        // create histograms for unfolding for different techniques
        histoGammaCorrUnfoldReso_Pt             = (TH1D*) histoGamma_OriginalBin_Pt->Clone("histoGammaCorrUnfoldReso_Pt");
        histoMCrecGammaCorr_Pt                  = (TH1D*) histoMCrecGamma_OriginalBin_Pt->Clone("GammaSpecCorrESDMC");
        histoGammaCorrUnfoldReso_BinByBin_Pt    = (TH1D*) histoGamma_OriginalBin_Pt->Clone("histoGammaCorrUnfoldReso_BinByBin_Pt");
        // TH1D *histoGammaCorrUnfoldReso_SvD_Pt       = (TH1D*) histoGamma_OriginalBin_Pt->Clone("histoGammaCorrUnfoldReso_SvD_Pt");
        // TH1D *histoGammaCorrUnfoldReso_TUnfold_Pt   = (TH1D*) histoGamma_OriginalBin_Pt->Clone("histoGammaCorrUnfoldReso_TUnfold_Pt");

        // Unfolding setup:
        // RooUnfoldResponse constructor - create from already-filled histograms
        // "response" gives the response matrix, measured X truth.
        // "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
        // but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
        // in "truth" for unmeasured events (inefficiency).
        // "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
        // to indicate, respectively, no fakes and/or no inefficiency.
        RooUnfoldResponse   response(0,0,histoGammaTruePrimary_recPt_MCPt);
        
        // unfold with different techniques, but same response matrix
        RooUnfoldBayes      unfold_Spectrum (&response,histoGammaCorrUnfoldReso_Pt, 4);
        RooUnfoldBayes      unfold_Spectrum (&response,histoMCrecGammaCorr_Pt, 4);
        RooUnfoldBinByBin   unfold_SpectrumBinByBin (&response,histoGammaCorrUnfoldReso_BinByBin_Pt);
        //RooUnfoldSvd      unfold_SpectrumSvD (&response, histoGammaCorrUnfoldReso_SvD_Pt, 20);   
        //RooUnfoldTUnfold  unfold_SpectrumTUnfold (&response,histoGammaCorrUnfoldReso_TUnfold_Pt);
        
        // get histograms from RooUnfold and rebin them 
        histoGammaCorrUnfoldReso_Pt             = (TH1D*) unfold_Spectrum.Hreco();
        histoMCrecGammaCorr_Pt                  = (TH1D*) unfold_Spectrum.Hreco();
        histoMCrecGammaCorr_Pt                  = RebinTH1D(histoMCrecGammaCorr_Pt,histoESDConvGammaPt,kTRUE);
        histoGammaResolCorrUnfold_Pt            = (TH1D*) histoGamma_OriginalBin_Pt->Clone("histoGammaResolCorrUnfold_Pt");
        histoGammaResolCorrUnfold_Pt->Divide(histoGammaCorrUnfoldReso_Pt);    
        histoGammaCorrUnfoldReso_Pt             = RebinTH1D(histoGammaCorrUnfoldReso_Pt,histoESDConvGammaPt,kTRUE);

        histoGammaCorrUnfoldReso_BinByBin_Pt    = (TH1D*) unfold_SpectrumBinByBin.Hreco();
        histoGammaResolCorrUnfold_BinByBin_Pt   = (TH1D*) histoGamma_OriginalBin_Pt->Clone("histoGammaResolCorrUnfold_BinByBin_Pt");
        histoGammaResolCorrUnfold_BinByBin_Pt->Divide(histoGammaCorrUnfoldReso_BinByBin_Pt);

        histoGammaCorrUnfoldReso_BinByBin_Pt    = RebinTH1D(histoGammaCorrUnfoldReso_BinByBin_Pt,histoESDConvGammaPt,kTRUE);
        // histoGammaCorrUnfoldReso_SvD_Pt      = (TH1D*) unfold_SpectrumSvD.Hreco();
        // histoGammaCorrUnfoldReso_SvD_Pt      = RebinTH1D(histoGammaCorrUnfoldReso_SvD_Pt,histoESDConvGammaPt,kTRUE);
        // histoGammaCorrUnfoldReso_TUnfold_Pt  = (TH1D*) unfold_SpectrumTUnfold.Hreco();
        // histoGammaCorrUnfoldReso_TUnfold_Pt  = RebinTH1D(histoGammaCorrUnfoldReso_TUnfold_Pt,histoESDConvGammaPt,kTRUE);

        // Correct inclusive photon spectrum with conversion probability & reco efficiency (both versus MC pt)
        histoGammaCorrUnfoldReso_PtNotCorrected = (TH1D*) histoGammaCorrUnfoldReso_Pt->Clone("GammaUnfoldNotCorrected");
        CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
        CorrectGammaUnfoldResol(histoMCrecGammaCorr_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
        CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_BinByBin_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
        //CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_SvD_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
        //CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_TUnfold_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
        
        // Correct inclusive photon spectrum for pileup contribution
        if(doPileUpCorr){
            histoGammaCorrUnfoldResoPileUp_Pt         = (TH1D*) histoGammaCorrUnfoldReso_Pt->Clone("GammaUnfoldPileUp");
            histoGammaCorrUnfoldResoPileUp_Pt->Multiply(histoPileUpCorrectionFactor);
        }
        
        SetHistogramm(histoGammaCorrUnfoldReso_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        DrawGammaSetMarker(histoGammaCorrUnfoldReso_Pt, 20, 1.0, 1, 1);
        SetHistogramm(histoGammaCorrUnfoldReso_BinByBin_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        DrawGammaSetMarker(histoGammaCorrUnfoldReso_BinByBin_Pt, 24, 1.0, kBlue, kBlue);
        SetHistogramm(histoMCrecGammaCorr_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        DrawGammaSetMarker(histoMCrecGammaCorr_Pt, 20, 1.0, kGreen-1, kGreen-1);
        //SetHistogramm(histoGammaCorrUnfoldReso_SvD_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        //SetHistogramm(histoGammaCorrUnfoldReso_TUnfold_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
    }    

    TH1D *histoGammaCaloCorrUnfoldReso_Pt                   = NULL;
    TH1D *histoGammaCaloCorrUnfoldReso_BinByBin_Pt          = NULL;
    TH1D *histoGammaCaloCorrUnfoldResoPileUp_Pt             = NULL;
    TH1D* histoGammaCaloResolCorrUnfold_Pt                  = NULL; 
    TH1D *histoGammaCaloCorrUnfoldReso_PtNotCorrected       = NULL;
    TH1D* histoGammaCaloResolCorrUnfold_BinByBin_Pt         = NULL; 
    
//     if (isPCM && isCalo){
//         // determine secondary contribution in orginal binning
//         TH1D *histoSecondaryGammaSpecPtOriginalBin             = (TH1D*) histoGamma_OriginalBin_Pt->Clone("SecondaryGammaSpecPt");
//         TH1D *histoSecondaryGammaFromXFromK0sSpecPtOriginalBin = (TH1D*) histoGamma_OriginalBin_Pt->Clone("SecondaryGammaSpecFromXFromK0sPt");
//         histoSecondaryGammaSpecPtOriginalBin->Multiply(histoFracAllGammaToSec_OriginalBin_Pt);
//         histoSecondaryGammaFromXFromK0sSpecPtOriginalBin->Multiply(histoFracAllGammaToSecFromXFromK0s_OriginalBin_Pt);
//         histoSecondaryGammaFromXFromK0sSpecPtOriginalBin->Scale(doubleAddFactorK0s);
// 
//         //subtract secondary contribution from inclusive spectrum
//         CorrectGammaSecAndPurity(histoGamma_OriginalBin_Pt, histoSecondaryGammaSpecPtOriginalBin, histoSecondaryGammaFromXFromK0sSpecPtOriginalBin, histoGammaTruePurity_OriginalBin_Pt );
//         
//         
//         // create histograms for unfolding for different techniques
//         histoGammaCorrUnfoldReso_Pt             = (TH1D*) histoGamma_OriginalBin_Pt->Clone("histoGammaCorrUnfoldReso_Pt");
//         histoGammaCorrUnfoldReso_BinByBin_Pt    = (TH1D*) histoGamma_OriginalBin_Pt->Clone("histoGammaCorrUnfoldReso_BinByBin_Pt");
//         // TH1D *histoGammaCorrUnfoldReso_SvD_Pt       = (TH1D*) histoGamma_OriginalBin_Pt->Clone("histoGammaCorrUnfoldReso_SvD_Pt");
//         // TH1D *histoGammaCorrUnfoldReso_TUnfold_Pt   = (TH1D*) histoGamma_OriginalBin_Pt->Clone("histoGammaCorrUnfoldReso_TUnfold_Pt");
// 
//         // Unfolding setup:
//         // RooUnfoldResponse constructor - create from already-filled histograms
//         // "response" gives the response matrix, measured X truth.
//         // "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
//         // but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
//         // in "truth" for unmeasured events (inefficiency).
//         // "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
//         // to indicate, respectively, no fakes and/or no inefficiency.
//         RooUnfoldResponse   response(0,0,histoGammaTruePrimary_recPt_MCPt);
//         
//         // unfold with different techniques, but same response matrix
//         RooUnfoldBayes      unfold_Spectrum (&response,histoGammaCorrUnfoldReso_Pt, 4);
//         RooUnfoldBinByBin   unfold_SpectrumBinByBin (&response,histoGammaCorrUnfoldReso_BinByBin_Pt);
//         //RooUnfoldSvd      unfold_SpectrumSvD (&response, histoGammaCorrUnfoldReso_SvD_Pt, 20);   
//         //RooUnfoldTUnfold  unfold_SpectrumTUnfold (&response,histoGammaCorrUnfoldReso_TUnfold_Pt);
//         
//         // get histograms from RooUnfold and rebin them 
//         histoGammaCorrUnfoldReso_Pt             = (TH1D*) unfold_Spectrum.Hreco();
//         histoGammaResolCorrUnfold_Pt            = (TH1D*) histoGamma_OriginalBin_Pt->Clone("histoGammaResolCorrUnfold_Pt");
//         histoGammaResolCorrUnfold_Pt->Divide(histoGammaCorrUnfoldReso_Pt);    
//         histoGammaCorrUnfoldReso_Pt             = RebinTH1D(histoGammaCorrUnfoldReso_Pt,histoESDConvGammaPt,kTRUE);
// 
//         histoGammaCorrUnfoldReso_BinByBin_Pt    = (TH1D*) unfold_SpectrumBinByBin.Hreco();
//         histoGammaResolCorrUnfold_BinByBin_Pt   = (TH1D*) histoGamma_OriginalBin_Pt->Clone("histoGammaResolCorrUnfold_BinByBin_Pt");
//         histoGammaResolCorrUnfold_BinByBin_Pt->Divide(histoGammaCorrUnfoldReso_BinByBin_Pt);
// 
//         histoGammaCorrUnfoldReso_BinByBin_Pt    = RebinTH1D(histoGammaCorrUnfoldReso_BinByBin_Pt,histoESDConvGammaPt,kTRUE);
//         // histoGammaCorrUnfoldReso_SvD_Pt      = (TH1D*) unfold_SpectrumSvD.Hreco();
//         // histoGammaCorrUnfoldReso_SvD_Pt      = RebinTH1D(histoGammaCorrUnfoldReso_SvD_Pt,histoESDConvGammaPt,kTRUE);
//         // histoGammaCorrUnfoldReso_TUnfold_Pt  = (TH1D*) unfold_SpectrumTUnfold.Hreco();
//         // histoGammaCorrUnfoldReso_TUnfold_Pt  = RebinTH1D(histoGammaCorrUnfoldReso_TUnfold_Pt,histoESDConvGammaPt,kTRUE);
// 
//         // Correct inclusive photon spectrum with conversion probability & reco efficiency (both versus MC pt)
//         histoGammaCorrUnfoldReso_PtNotCorrected = (TH1D*) histoGammaCorrUnfoldReso_Pt->Clone("GammaUnfoldNotCorrected");
//         CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
//         CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_BinByBin_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
//         //CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_SvD_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
//         //CorrectGammaUnfoldResol(histoGammaCorrUnfoldReso_TUnfold_Pt,histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_MCPt, deltaEta, scaling, nEvt);
//         
//         // Correct inclusive photon spectrum for pileup contribution
//         
//         if(doPileUpCorr){
//             histoGammaCorrUnfoldResoPileUp_Pt         = (TH1D*) histoGammaCorrUnfoldReso_Pt->Clone("GammaUnfoldPileUp");
//             histoGammaCorrUnfoldResoPileUp_Pt->Multiply(histoPileUpCorrectionFactor);
//         }
//         SetHistogramm(histoGammaCorrUnfoldReso_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
//         DrawGammaSetMarker(histoGammaCorrUnfoldReso_Pt, 20, 1.0, 1, 1);
//         SetHistogramm(histoGammaCorrUnfoldReso_BinByBin_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
//         DrawGammaSetMarker(histoGammaCorrUnfoldReso_BinByBin_Pt, 24, 1.0, kBlue, kBlue);
//         //SetHistogramm(histoGammaCorrUnfoldReso_SvD_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
//         //SetHistogramm(histoGammaCorrUnfoldReso_TUnfold_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
//     }    
    
    
    // Plotting corrected photon spectrum
    TCanvas *canvasCorrGammaSpecUnfold          = GetAndSetCanvas("canvasCorrGammaSpecUnfold");
    canvasCorrGammaSpecUnfold->SetLogy();
        TH1D *Dummy                             = (TH1D*)histoESDConvGammaPt->Clone("Dummy");
        SetHistogramm(Dummy,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",1e-8,10);
        Dummy->Reset();
        Dummy->Draw();
        histoGammaCorrUnfoldReso_Pt->DrawCopy("same");
        histoGammaCorrUnfoldReso_BinByBin_Pt->DrawCopy("same");
        //histoGammaCorrUnfoldReso_SvD_Pt->DrawCopy("same");
        //histoGammaCorrUnfoldReso_TUnfold_Pt->DrawCopy("same");

        TLegend* legendGammaSpecUnfoldComp                = GetAndSetLegend(0.6,0.8,2);
        legendGammaSpecUnfoldComp->AddEntry(histoGammaCorrUnfoldReso_Pt,"Bayesian unfolding");
        legendGammaSpecUnfoldComp->AddEntry(histoGammaCorrUnfoldReso_BinByBin_Pt,"Bin-by-Bin unfolding");
        legendGammaSpecUnfoldComp->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.68, 0.78, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
    canvasCorrGammaSpecUnfold->SaveAs(Form("%s/%s_%s_CorrGammaUnfoldSpectrumPurityMinusSec_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasCorrGammaSpecUnfold;

    // plotting resolution correction
    TCanvas *canvasResolutionCorr          = GetAndSetCanvas("canvasResolutionCorr");
        TH1D *Dummy2                            = (TH1D*)histoESDConvGammaPt->Clone("Dummy2");
        SetHistogramm(Dummy2,"#it{p}_{T} (GeV/#it{c})", "resolution correction",0,2);
        Dummy2->Reset();
        Dummy2->Draw();

        DrawGammaSetMarker(histoGammaResolCorrUnfold_Pt, 20, 1.0, kBlue+1, kBlue+1);
        DrawGammaSetMarker(histoGammaResolCorrUnfold_BinByBin_Pt, 21, 1.0, kGray+2, kGray+2);
        DrawGammaSetMarker(histoGammaResolCorrEff_Pt, 24, 1.5, kBlack, kBlack);

        DrawGammaLines(0., maxPtGamma,1, 1,1, kGray);
        DrawGammaLines(0., maxPtGamma,0.8, 0.8,1, kGray, 7);
        DrawGammaLines(0., maxPtGamma,1.2, 1.2,1, kGray, 7);
        
        histoGammaResolCorrUnfold_BinByBin_Pt->DrawCopy("same");
        histoGammaResolCorrUnfold_Pt->DrawCopy("same");
        histoGammaResolCorrEff_Pt->DrawCopy("same");
        //histoGammaCorrUnfoldReso_SvD_Pt->DrawCopy("same");
        //histoGammaCorrUnfoldReso_TUnfold_Pt->DrawCopy("same");

        TLegend* legendResolutionCorr                = GetAndSetLegend(0.15,0.80,3);
        legendResolutionCorr->AddEntry(histoGammaResolCorrEff_Pt,"from effi ");
        legendResolutionCorr->AddEntry(histoGammaResolCorrUnfold_Pt,"from Bayesian unfolding");
        legendResolutionCorr->AddEntry(histoGammaResolCorrUnfold_BinByBin_Pt,"from bin-by-bin unfolding");
        legendResolutionCorr->Draw();

        PutProcessLabelAndEnergyOnPlot( 0.18, 0.3, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);

 
    canvasResolutionCorr->SaveAs(Form("%s/%s_%s_ResolutionCorrection_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasResolutionCorr;

    
    //****************************************************************************************** 
    //****************************************************************************************** 
    //****************************************************************************************** 
    cout << __LINE__ << endl;

    
    // copy raw inc gamma spectrum and correct for secondary contamination, purity, conversion probability & reco effi (vs rec pt with sec correction)
    TH1D* histoGammaCorrEffiReso_Pt      = (TH1D*)histoESDConvGammaPt->Clone("CorrGammaSpecPurityMinusSec");
    CorrectGammaEffiResol( histoGammaCorrEffiReso_Pt, histoSecondaryGammaSpecPt, histoSecondaryGammaFromXFromK0sSpecPt, 
                          histoGammaTruePurity_Pt, histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_Pt, deltaEta, scaling, nEvt);
    DrawGammaSetMarker(histoGammaCorrEffiReso_Pt, 20, 1.0, kGreen+2, kGreen+2);
    SetHistogramm(histoGammaCorrEffiReso_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
    
    // Do same as before with secondary corrections applied
    TH1D* histoGammaCorrEffiReso_PileUp_Pt;
    TH1D* histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt;
    if(doPileUpCorr){
        // copy raw inc gamma spectrum after pure data driven pileup correction, correct with pileup corrected secondary histograms, pileup corrected purity, conversion prob, reco effi 
        // (pileup corrected, vs rec pT with sec correction) all contain MC update of pileup correction
        histoGammaCorrEffiReso_PileUp_Pt              = (TH1D*)histoESDConvGammaPtPileUp->Clone("CorrGammaSpecPurityMinusSecPileUp");
        CorrectGammaEffiResol( histoGammaCorrEffiReso_PileUp_Pt, histoSecondaryGammaSpecPtPileUp, histoSecondaryGammaFromXFromK0sSpecPtPileUp, 
                              histoGammaTruePurity_PileUp_Pt, histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_PileUp_Pt, deltaEta, scaling, nEvt);
        
        DrawGammaSetMarker(histoGammaCorrEffiReso_PileUp_Pt, 20, 1.0, kMagenta, kMagenta);
        SetHistogramm(histoGammaCorrEffiReso_PileUp_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");

        // copy raw inc gamma spectrum after pure data driven pileup correction, correct with pileup corrected secondary histograms, pileup corrected purity, conversion prob, reco effi 
        // (pileup corrected, vs rec pT with sec correction)
        histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt    = (TH1D*)histoESDConvGammaPtPileUp->Clone("CorrGammaSpecPurityMinusSecPileUpNoMCUpdate");
        CorrectGammaEffiResol( histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt, histoSecondaryGammaSpecPt, histoSecondaryGammaFromXFromK0sSpecPt, 
                              histoGammaTruePurity_Pt, histoGammaConvProb_MCPt, histoGammaPrimaryRecoEff_Pt, deltaEta, scaling, nEvt);
        
        DrawGammaSetMarker(histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt, 20, 1.0, kBlue+2, kBlue+2);
        SetHistogramm(histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt->Draw("same");
    }

   
    //*************************************************************************************************    
    //********** Sanity check for photon spectrum (do we get back what we put in) *********************
    //*************************************************************************************************
    // correct input inc gamma spectrum for delta Eta, delta phi & number of events
    TH1D* histoMCGammaSpec_MCPt                      = (TH1D*)histoMCAllGamma_MCPt->Clone("GammaSpecMC");
    CorrectGammaMC(histoMCGammaSpec_MCPt, deltaEta, scaling, nEvtMC);
    DrawGammaSetMarker(histoMCGammaSpec_MCPt, 20, 1.0, 2, 2);
    SetHistogramm(histoMCGammaSpec_MCPt,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
            
    //*************************************************************************************************
    //********************************* Comparison Gamma Spec *****************************************
    //*************************************************************************************************

    TH1D* histoRatioGammaCorrMinusSecGammaUnfold    = (TH1D*)histoGammaCorrEffiReso_Pt->Clone("histoRatioGammaCorrMinusSecGammaUnfold");
    histoRatioGammaCorrMinusSecGammaUnfold->Divide(histoGammaCorrUnfoldReso_Pt);
    TH1D* histoRatioGammaCorrMinusSecGammaPileUp;
    TH1D* histoRatioGammaCorrMinusSecGammaPileUpNoMCUpdate;
    if(doPileUpCorr){
        histoRatioGammaCorrMinusSecGammaPileUp           = (TH1D*)histoGammaCorrEffiReso_Pt->Clone("histoRatioGammaCorrMinusSecGammaPileUp");
        histoRatioGammaCorrMinusSecGammaPileUp->Divide(histoGammaCorrEffiReso_PileUp_Pt);

        histoRatioGammaCorrMinusSecGammaPileUpNoMCUpdate = (TH1D*)histoGammaCorrEffiReso_Pt->Clone("histoRatioGammaCorrMinusSecGammaPileUp");
        histoRatioGammaCorrMinusSecGammaPileUpNoMCUpdate->Divide(histoGammaCorrEffiReso_PileUpNoMCUpdate_Pt);
    }
    
//     TCanvas *canvasGammaSpecComp                    = GetAndSetCanvas("canvasGammaSpecComp");
// //     canvasGammaSpecComp->SetGridx();
// //     canvasGammaSpecComp->SetGridy();
//     
//         SetHistogramm(histoRatioGammaCorrMinusSecGammaUnfold,"#it{p}_{T} (GeV/#it{c})", "#gamma / gamma_{Mod}",0.8,1.3);
//         DrawGammaSetMarker(histoRatioGammaCorrMinusSecGammaUnfold, 24, 1.0, 1, 1);
//         histoRatioGammaCorrMinusSecGammaUnfold->DrawCopy();
// 
//         if(doPileUpCorr){
//             SetHistogramm(histoRatioGammaCorrMinusSecGammaPileUp,"#it{p}_{T} (GeV/#it{c})", "#gamma / gamma_{Mod}");
//             DrawGammaSetMarker(histoRatioGammaCorrMinusSecGammaPileUp, 26, 1.0, 3, 3);
//             SetHistogramm(histoRatioGammaCorrMinusSecGammaPileUpNoMCUpdate,"#it{p}_{T} (GeV/#it{c})", "#gamma / gamma_{Mod}");
//             DrawGammaSetMarker(histoRatioGammaCorrMinusSecGammaPileUpNoMCUpdate, 26, 1.0, kOrange, kOrange);
//             histoRatioGammaCorrMinusSecGammaPileUp->DrawCopy("same");
//             histoRatioGammaCorrMinusSecGammaPileUpNoMCUpdate->DrawCopy("same");
//         }
// 
//         TLegend* legendGammaSpecComp                = GetAndSetLegend(0.4,0.25,5);
//         legendGammaSpecComp->AddEntry(histoRatioGammaCorrMinusSecGammaUnfold,"Gamma Minus Sec - Unfold");
//         if(doPileUpCorr)legendGammaSpecComp->AddEntry(histoRatioGammaCorrMinusSecGammaPileUp,"Gamma Minus Sec - Pile Up");
//         legendGammaSpecComp->Draw();
// 
//     canvasGammaSpecComp->SaveAs(Form("%s/%s_%s_GammaSpectraRatios_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data(),suffix.Data()));
//     delete canvasGammaSpecComp;
//     
    //*************************************************************************************************
    //********************************* Gamma from Decay **********************************************
    //*************************************************************************************************
    // correct input spectra of different decay channels with delta eta, delta phi & events number
    for (Int_t i = 0; i< 9; i++){
        CorrectGammaMC(histoPhotonSource_MCPt[i], deltaEta, scaling, nEvtMC);
        SetHistogramm(histoPhotonSource_MCPt[i],"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
        DrawGammaSetMarker(histoPhotonSource_MCPt[i], 22, 1.0, colorsDecay[i] , colorsDecay[i]);
    }
        
    TCanvas *canvasDecayGammaSpecMC = GetAndSetCanvas("canvasDecayGammaSpecMC");
    canvasDecayGammaSpecMC->SetLogy();
    
        TLegend* legendDecaySpectra = GetAndSetLegend(0.5,0.7,5,2);
        histoPhotonSource_MCPt[7]->DrawCopy("chist");
        for (Int_t i = 0; i< 9; i++){   
            histoPhotonSource_MCPt[i]->DrawCopy("chist,same");
            legendDecaySpectra->AddEntry(histoPhotonSource_MCPt[i],decaysLabels[i].Data());
        }
        legendDecaySpectra->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.7, 0.68, 0.035, cent, "#gamma", detectionProcess, 42, 0.03);
        
    canvasDecayGammaSpecMC->SaveAs(Form("%s/%s_DecayGammaSpectrumMC_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasDecayGammaSpecMC;

    //*************************************************************************************************
    //******************************* Gamma Combinatorial Background **********************************
    //*************************************************************************************************
    // pure BG plotting
    TCanvas *canvasCombBackSpecMC = GetAndSetCanvas("canvasCombBackSpecMC");
    canvasCombBackSpecMC->SetLogy();

        TLegend* legendCombSpectra = GetAndSetLegend(0.25,0.68,6,3);
        for(Int_t i = 0;i<17;i++){
            SetHistogramm(histoCombinatorialSpecies_Pt[i],"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}");
            DrawGammaSetMarker(histoCombinatorialSpecies_Pt[i], markersCombinatorics[i], 1., colorsCombinatorics[i], colorsCombinatorics[i]);
            histoCombinatorialSpecies_Pt[i]->Scale(1./nEvtMC);

            if(i==0) histoCombinatorialSpecies_Pt[i]->DrawCopy("");
            else histoCombinatorialSpecies_Pt[i]->DrawCopy("same");

            legendCombSpectra->AddEntry(histoCombinatorialSpecies_Pt[i],combinatorics[i]);
        }

        legendCombSpectra->Draw();
        PutProcessLabelAndEnergyOnPlot( 0.75, 0.68, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);

        
    canvasCombBackSpecMC->SaveAs(Form("%s/%s_CombBackSpectrumMC_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasCombBackSpecMC;

    // plot ratio to real primary photons
    TCanvas *canvasSignalToCombBackgroundRatio = GetAndSetCanvas("canvasSignalToCombBackgroundRatio");
    canvasSignalToCombBackgroundRatio->SetLogy();

        TLegend* legendSignalToCombBackgroundRatio = GetAndSetLegend(0.15,0.76,4,3);
        TH1D **histoSignalToCombBackgroundRatio = new TH1D*[17];
        TH1D *SummedSmallContributionsCombBack; //10,11,12,13,14,15
        
        for(Int_t i = 0;i<16;i++){
            histoSignalToCombBackgroundRatio[i] = (TH1D*) histoCombinatorialSpecies_Pt[i]->Clone(Form("ESD_TrueCombRatioSignal_%s_Pt",combinatorics[i].Data()));
            histoSignalToCombBackgroundRatio[i]->Scale(nEvtMC);
            
            if(i==9){
                SummedSmallContributionsCombBack = (TH1D*)histoSignalToCombBackgroundRatio[i]->Clone("SummedSmallContributions");
                SetHistogramm(SummedSmallContributionsCombBack,"#it{p}_{T} (GeV/#it{c})","SummedSmallContributions",10,5e7);
                SummedSmallContributionsCombBack->SetMinimum(1e-5);
            }
            
            if(i>9){
                SummedSmallContributionsCombBack->Add(histoSignalToCombBackgroundRatio[i]);
            }
            histoSignalToCombBackgroundRatio[i]->Divide(histoSignalToCombBackgroundRatio[i],histoGammaTruePrimaryConv_Pt,1,1,"");
            SetHistogramm(histoSignalToCombBackgroundRatio[i],"#it{p}_{T} (GeV/#it{c})","identified BG/ real primary photons",10,5e7);
            histoSignalToCombBackgroundRatio[i]->SetMinimum(1e-5);
            
            if(i==0){
                histoSignalToCombBackgroundRatio[i]->GetYaxis()->SetRangeUser(1e-5,100);
                histoSignalToCombBackgroundRatio[i]->DrawCopy("e1");
            }
            if(i<9){
                DrawGammaSetMarker(histoSignalToCombBackgroundRatio[i], markersCombinatorics[i], 1., colorsCombinatorics[i], colorsCombinatorics[i]);
                histoSignalToCombBackgroundRatio[i]->DrawCopy("e1same");
            }
            else continue;
            
            legendSignalToCombBackgroundRatio->AddEntry(histoSignalToCombBackgroundRatio[i],combinatorics[i]);
        }

        SummedSmallContributionsCombBack->Divide(SummedSmallContributionsCombBack,histoGammaTruePrimaryConv_Pt,1,1,"");
        DrawGammaSetMarker(SummedSmallContributionsCombBack, markersCombinatorics[10], 1., colorsCombinatorics[10], colorsCombinatorics[10]);
       
        SummedSmallContributionsCombBack->DrawCopy("e1same");
        legendSignalToCombBackgroundRatio->AddEntry(SummedSmallContributionsCombBack,"Protos+Kaons+Muons");
        legendSignalToCombBackgroundRatio->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.2, 0.76, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
    canvasSignalToCombBackgroundRatio->SaveAs(Form("%s/%s_CombBackgroundRatioToSignal_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasSignalToCombBackgroundRatio;

    // plot ratio to summed total MC BG
    TCanvas *canvasRatioCombBackToBack = GetAndSetCanvas("canvasRatioCombBackToBack");
    canvasRatioCombBackToBack->SetLogy();

        TLegend* legendRatioCombBackToBack = GetAndSetLegend(0.15,0.76,4,3);
        
        TH1D **histoRatioCombBackToBack = new TH1D*[17];
        TH1D *SummedSmallContributionsCombBackToBack; //10,11,12,13,14,15
        for(Int_t i = 0;i<16;i++){

            histoRatioCombBackToBack[i] = (TH1D*) histoCombinatorialSpecies_Pt[i]->Clone(Form("ESD_TrueCombRatioSignal_%s_Pt",combinatorics[i].Data()));
            
            if(i==9){
                SummedSmallContributionsCombBackToBack = (TH1D*)histoRatioCombBackToBack[i]->Clone("SummedSmallContributions");
                SetHistogramm(SummedSmallContributionsCombBackToBack,"#it{p}_{T} (GeV/#it{c})","SummedSmallContributions",10,5e7);
                SummedSmallContributionsCombBackToBack->SetMinimum(1e-5);
            }
            
            if(i>9){
                SummedSmallContributionsCombBackToBack->Add(histoRatioCombBackToBack[i]);
            }
            histoRatioCombBackToBack[i]->Divide(histoRatioCombBackToBack[i],histoGammaMCBackground_Pt,1,1,"B");
            SetHistogramm(histoRatioCombBackToBack[i],"#it{p}_{T} (GeV/#it{c})","identified BG/Total BG",10,5e7);
            histoRatioCombBackToBack[i]->SetMinimum(1e-5);
            
            if(i==0){
                histoRatioCombBackToBack[i]->GetYaxis()->SetRangeUser(1e-4,400);
                histoRatioCombBackToBack[i]->DrawCopy("e1");
            }
            if(i<9){
                DrawGammaSetMarker(histoRatioCombBackToBack[i], markersCombinatorics[i], 1., colorsCombinatorics[i], colorsCombinatorics[i]);
                histoRatioCombBackToBack[i]->DrawCopy("e1same");
            }
            else continue;
            
            legendRatioCombBackToBack->AddEntry(histoRatioCombBackToBack[i],combinatorics[i]);
        }
        
        SummedSmallContributionsCombBackToBack->Divide(SummedSmallContributionsCombBackToBack,histoGammaMCBackground_Pt,1,1,"B");
        DrawGammaSetMarker(SummedSmallContributionsCombBackToBack, markersCombinatorics[10], 1., colorsCombinatorics[10], colorsCombinatorics[10]);
        SummedSmallContributionsCombBackToBack->DrawCopy("e1same");
        legendRatioCombBackToBack->AddEntry(SummedSmallContributionsCombBackToBack,"Protos+Kaons+Muons");
        legendRatioCombBackToBack->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.75, 0.3, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
    
    canvasRatioCombBackToBack->SaveAs(Form("%s/%s_RatioCombBackToBack_%s.%s",outputDir.Data(),textPi0New.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasRatioCombBackToBack;

    
    //*************************************************************************************************
    //*********************** Compare Gamma Yields MC/Reconstructed *************************************
    //*************************************************************************************************
    TCanvas *canvasPurBackGammaSpec = GetAndSetCanvas("canvasPurBackGammaSpec");
    DrawGammaCanvasSettings( canvasPurBackGammaSpec, 0.1, 0.015, 0.02, 0.09);
        TPad *padPurBackSpec = new TPad("padNLOHistos", "", 0., 0.25, 1., 1.,-1, -1, -2);
        DrawGammaPadSettings( padPurBackSpec, 0.10, 0.015, 0.02, 0.);
        padPurBackSpec->Draw();
        TPad *padPurBackRatio = new TPad("padNLORatios", "", 0., 0., 1., 0.25,-1, -1, -2);
        DrawGammaPadSettings( padPurBackRatio, 0.1, 0.015, 0.0, 0.35);
        padPurBackRatio->Draw();

        padPurBackSpec->cd();
        padPurBackSpec->SetLogy();

        DrawGammaSetMarker(histoGammaCorrEffiReso_Pt, 20, 1.0, kBlue+2,  kBlue+2);
        DrawGammaSetMarker(histoGammaCorrUnfoldReso_Pt, 24, 1.0, kAzure+7,  kAzure+7);

        histoGammaCorrEffiReso_Pt->DrawCopy("e1");
        histoGammaCorrUnfoldReso_Pt->DrawCopy("e1,same");
        
        TLegend* legendPurBackSpectra = GetAndSetLegend(0.5,0.85,2);
        legendPurBackSpectra->AddEntry(histoGammaCorrEffiReso_Pt,"#gamma with Purity and corrected for Secondaries","p");
        legendPurBackSpectra->AddEntry(histoGammaCorrUnfoldReso_Pt,"#gamma unfolded with Secondary correction","p");
        legendPurBackSpectra->Draw();
        
        PutProcessLabelAndEnergyOnPlot( 0.18, 0.2, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
        padPurBackRatio->cd();
        TH1D *histGammaRatioUnfoldMinusSec = (TH1D*) histoGammaCorrEffiReso_Pt->Clone("histGammaRatioMinusSec");
        histGammaRatioUnfoldMinusSec->Divide(histGammaRatioUnfoldMinusSec,histoGammaCorrUnfoldReso_Pt,1,1,"");
        histGammaRatioUnfoldMinusSec->GetYaxis()->SetRangeUser(0.88,1.12);
        SetStyleHistoTH1ForGraphs(  histGammaRatioUnfoldMinusSec, "#it{p}_{T} (GeV/#it{c})","#frac{#gamma corr legacy}{#gamma corr unfold}" , 0.14, 0.15,  
                                    0.12, 0.10,  0.85, 0.4, 510, 505);
        
        DrawGammaSetMarker(histGammaRatioUnfoldMinusSec, 24, 1.0, kAzure+7,  kAzure+7);
        histGammaRatioUnfoldMinusSec->DrawCopy("");
        
        DrawGammaLines(0., maxPtGamma,1, 1,0.5);
    
    canvasPurBackGammaSpec->SaveAs(Form("%s/%s_%s_GammaSpectraComparison_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasPurBackGammaSpec;
        
    
    CorrectGammaMC(histoMCGammaConv_MCPt, deltaEta, scaling, nEvtMC);
    CorrectGammaMC(histoGammaTrueConv_Pt, deltaEta, scaling, nEvtMC);

    TCanvas *canvasGammaPurityConvBin = GetAndSetCanvas("canvasGammaPurityConvBin");
    canvasGammaPurityConvBin->SetLogy();

        DrawGammaSetMarker(histoGammaCorrEffiReso_Pt, 24, 1.0, 1, 1);
        
        histoMCGammaSpec_MCPt->GetYaxis()->SetRangeUser(5e-10,100);
        histoMCGammaSpec_MCPt->DrawCopy("e1");
        histoGammaCorrUnfoldReso_Pt->SetMarkerColor(1);
        histoGammaCorrUnfoldReso_Pt->SetLineColor(1);
        histoGammaCorrUnfoldReso_Pt->DrawCopy("e1,same");
        TLegend *legendGammaSpectraConvBin = GetAndSetLegend(0.55,0.8,2);
        legendGammaSpectraConvBin->AddEntry(histoGammaCorrUnfoldReso_Pt,"corrected data #gamma spectrum");
        legendGammaSpectraConvBin->AddEntry(histoMCGammaSpec_MCPt,"input MC #gamma spectrum");
        legendGammaSpectraConvBin->Draw();
    
        PutProcessLabelAndEnergyOnPlot( 0.18, 0.3, 0.035, cent, textMeasurement, detectionProcess, 42, 0.03);
        
    canvasGammaPurityConvBin->SaveAs(Form("%s/%s_%s_GammaSpec_%s.%s",outputDir.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data(),suffix.Data()));
    delete canvasGammaPurityConvBin;
    
    TH1D* histoGammaCorrFac_Pt           = (TH1D*) histoGammaPrimaryRecoEff_Pt->Clone("GammaCorrFac_Pt");
    histoGammaCorrFac_Pt->Multiply(histoGammaConvProb_MCPt);
    
    
    //***********************************************************************************************************
    //******************************* Write output file for photons *********************************************
    //***********************************************************************************************************
    TString nameOutput = Form("%s/%s/%s_%s_GammaConvV1Correction_%s.root",cutSelection.Data(),option.Data(),textPi0New.Data(),textPrefix2.Data(),cutSelection.Data());
    TFile* fileCorrectedOutput = new TFile(nameOutput,"RECREATE");

        //________________________ writing MC quantities to file _____________________________
        // input spectrum corrected for deta, dPhi, Nevt
        histoMCGammaSpec_MCPt->Write("GammaSpecMC");
    
        // reconstructed MC gamma spectrum corrected deta, dPhi, Nevt
        histoMCrecGammaCorr_Pt->Write("GammaSpecCorrESDMC");
    
        // split in different source (pi0,eta,...)
        for (Int_t i= 0; i< 8; i++){
            histoPhotonSource_MCPt[i]->Write();
        }    
        // pure direct photons
        histoPhotonSource_MCPt[8]->Write("MC_DirectPhoton_Pt");
        // input spectrum of converted photons corrected for deta, dphi, Nevt
        histoMCGammaConv_MCPt->Write("MC_ConvGamma_MCPt");
        // reconstructed MC true photons corrected for deta, dphi, Nevt
        histoGammaTrueConv_Pt->Write("TrueConvGamma_Pt");
        // recontructed MC photon candidates, uncorrected
        histoMCrecGamma_Pt->Write("MCrec_ConvGamma_Pt");
        
        // summed reconstructed BG in MC
        histoMCrecBackground_Pt->Write("MCrec_Background");
        // combinatorial histos split
        for(Int_t i = 0;i<17;i++)
            histoCombinatorialSpecies_Pt[i]->Write();
        
        //_________________________ writing correction factors to file _________________________        
        // photon purity without secondary subtraction
        histoGammaPurity_Pt->Write("GammaPurityWSec_Pt");
        // photon purity with secondary subtraction
        histoGammaTruePurity_Pt->Write("GammaPurityWOSec_Pt");
        // photon reconstruction efficiency including resolution correction
        histoGammaPrimaryRecoEff_Pt->Write("GammaRecoEff_WithResolCorr_Pt");
        // photon reconstruction efficiency without resolution correction
        histoGammaPrimaryRecoEff_MCPt->Write("GammaRecoEff_MCPt");
        // photon conversion probability
        histoGammaConvProb_MCPt->Write("GammaConvProb_Pt");
        // resolution correction in case of unfolding
        histoGammaResolCorrUnfold_Pt->Write("GammaResolCorrUnfold_Pt");
        // photon correction factors (conv Prob, efficiency incl. resolution correction)
        histoGammaCorrFac_Pt->Write("GammaCorrFac_Pt");
        
        // ________________________ writing data quantities to file
        // uncorrected spectrum (scaled by 1/Nevt)
        histoGammaRawSpectrum_Pt->Write("GammaRaw_Pt");
        
        // corrected spectrum (legacy corrections: purity, effi incl resolution correction, secondaries, conv prob, plus trivial factors  )
        histoGammaCorrEffiReso_Pt->Write("GammaCorrEffiResol_Pt");
        // corrected spectrum (unfolding corrections: purity, secondaries, unfolding resolution correction, effi without unfolding corr, conv prob, plus trivial factors )
        histoGammaCorrUnfoldReso_Pt->Write("GammaCorrUnfold_Pt");
        if(doPileUpCorr){
            // same as histoGammaCorrEffiReso_Pt with additional pileup correction
            histoGammaCorrEffiReso_PileUp_Pt->Write("GammaCorrEffiResolPileup_Pt");
            // same as histoGammaCorrUnfoldReso_Pt with additional pileup correction
            histoGammaCorrUnfoldResoPileUp_Pt->Write("GammaCorrUnfoldPileUp_Pt");
        }
        
    fileCorrectedOutput->Close();

}
