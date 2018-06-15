/****************************************************************************************************************************
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

extern TRandom* gRandom;
extern TBenchmark* gBenchmark;
extern TSystem* gSystem;
extern TMinuit* gMinuit;

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
    } else {
        if (optionDalitz.CompareTo("kFALSE")==0){
            histoCorrectedToBeScaled->Scale(1./0.3931);
        } else {
            histoCorrectedToBeScaled->Scale(1./6.8e-5);
        }
    }
}

void PrepareWeightingIterations_pp5TeV2017(TString suffix = "pdf",
                                           TString nameFilepp5TeV2017 = "00010113_00200009227300008250404000_0152103500000000-LHC17pq-fastwoSDD/data_PCMResultsFullCorrection_PP.root",  Bool_t runDrawReweighted = kTRUE,
                                           TString stringIterationNumber = "Iter. 0",
                                           TString nameRecCls="fastwoSDD"){

    gROOT->Reset();
    gROOT->SetStyle("Plain");

    TString dateForOutput = ReturnDateStringForOutput();
    TString outputDir = Form("%s/%s/PrepareWeightingIterations_pp",suffix.Data(),dateForOutput.Data());
    TString cutNumber = "00010113_00200009227300008250404000_0152103500000000";
    TString cutNumbAdd= "";
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputFilePCMPP5TeV2017.root ",nameFilepp5TeV2017.Data(),outputDir.Data() ));

    Bool_t Iteration0 = kFALSE;

    if(  stringIterationNumber.CompareTo("0") == 0 ){
        Iteration0 = kTRUE;
    }


    StyleSettingsThesis();
    SetPlotStyle();

    Color_t color5TeV2017                   = kMagenta+2;
    Color_t colorMCPythiaPP5TeV2017         = color5TeV2017+2;
    Color_t colorMCPythiaAddSigPP5TeV2017   = color5TeV2017-6;
    Color_t colorMCPythiaJetJetPP5TeV2017   = color5TeV2017+6;
    Color_t colorMCPhojetPP5TeV2017         = color5TeV2017-4;

    Style_t markerStyleSpectrum5TeV2017     = 29 ;
    Size_t  markerSizePP5TeV2017            = 2.2;

    TString collisionSystemPP5TeV2017       = "pp #sqrt{#it{s}} = 5 TeV (2017)";

    Style_t lineStyleMCA                    = 1;
    Style_t lineStyleMCB                    = 2;
    Style_t lineStyleMCC                    = 3;
    Style_t lineStyleMCAddSig               = 4;
    Style_t lineStyleMCJetJet               = 5;

    Style_t markerStyleMCA                  = 20;
    Style_t markerStyleMCB                  = 25;
    Style_t markerStyleMCC                  = 30;
    Style_t markerStyleMCAddSig             = 29;
    Style_t markerStyleMCJetJet             = 32;

    // read reconstructed data
    TFile*      filePCMpp5TeV2017                           = new TFile(nameFilepp5TeV2017);
    TDirectory* directoryPCMPi05TeV2017                     = (TDirectory*)filePCMpp5TeV2017->Get("Pi05TeV2017");
    TH1D*       histoPCMYieldPi05TeV2017                    = (TH1D*)directoryPCMPi05TeV2017->Get("CorrectedYieldPi0");
    TH1D*       histoPi0InputFullReweighted5TeV2017         = (TH1D*)directoryPCMPi05TeV2017->Get("Pi0_Input_Reweighted");
    TH1D*       histoPi0InputFull5TeV2017                   = (TH1D*)directoryPCMPi05TeV2017->Get("Pi0_Input");
    TH1D*       histoPi0InputFullAddSigReweighted5TeV2017   = (TH1D*)directoryPCMPi05TeV2017->Get("Pi0_Input_Reweighted_AddedSig");
    TH1D*       histoPi0InputFullAddSig5TeV2017             = (TH1D*)directoryPCMPi05TeV2017->Get("Pi0_Input_AddedSig");

    TDirectory* directoryPCMEta5TeV2017                     = (TDirectory*)filePCMpp5TeV2017->Get("Eta5TeV2017");
    TH1D*       histoPCMYieldEta5TeV2017                    = (TH1D*)directoryPCMEta5TeV2017->Get("CorrectedYieldEta");
    TH1D*       histoEtaInputFullReweighted5TeV2017         = (TH1D*)directoryPCMEta5TeV2017->Get("Eta_Input_Reweighted");
    TH1D*       histoEtaInputFull5TeV2017                   = (TH1D*)directoryPCMEta5TeV2017->Get("Eta_Input");
    TH1D*       histoEtaInputFullAddSigReweighted5TeV2017   = (TH1D*)directoryPCMEta5TeV2017->Get("Eta_Input_Reweighted_AddedSig");
    TH1D*       histoEtaInputFullAddSig5TeV2017             = (TH1D*)directoryPCMEta5TeV2017->Get("Eta_Input_AddedSig");

    // Pythia 8
    TString     fileNamePi0Pythia8 = Form("%s-LHC17pq-%s/5TeV2017/Pi0_data_GammaConvV1Correction_%s.root",cutNumber.Data(),nameRecCls.Data(),cutNumber.Data());
    TFile*      filePi0PCM5TeV2017_LHC17pq_PYT8                 = new TFile(fileNamePi0Pythia8.Data());
    TH1D*       histoPi0InputMCWOWeights_LHC17pq_PYT8_5TeV2017  = (TH1D*)filePi0PCM5TeV2017_LHC17pq_PYT8->Get("MCYield_Meson_oldBinWOWeights");
    TH1D*       histoPi0Efficiency_LHC17pq_PYT8_5TeV2017        = (TH1D*)filePi0PCM5TeV2017_LHC17pq_PYT8->Get("TrueMesonEffiPt");

    TString     fileNameEtaPythia8 = Form("%s-LHC17pq-%s/5TeV2017/Eta_data_GammaConvV1Correction_%s.root",cutNumber.Data(),nameRecCls.Data(),cutNumber.Data());
    TFile*      fileEtaPCM5TeV2017_LHC17pq_PYT8                 = new TFile(fileNameEtaPythia8.Data());
    TH1D*       histoEtaInputMCWOWeights_LHC17pq_PYT8_5TeV2017  = (TH1D*)fileEtaPCM5TeV2017_LHC17pq_PYT8->Get("MCYield_Meson_oldBinWOWeights");
    TH1D*       histoEtaEfficiency_LHC17pq_PYT8_5TeV2017        = (TH1D*)fileEtaPCM5TeV2017_LHC17pq_PYT8->Get("TrueMesonEffiPt");

    // Phojet
    TString     fileNamePi0Phojet = Form("%s-LHC17pLowInt-%s/5TeV2017/Pi0_data_GammaConvV1Correction_%s.root",cutNumber.Data(),nameRecCls.Data(),cutNumber.Data());
    TFile*      filePi0PCM5TeV2017_LHC17pLowInt_PHO = new TFile(fileNamePi0Phojet.Data());
    TH1D*       histoPi0InputMCWOWeights_LHC17pLowInt_PHO_5TeV2017 = NULL;
    TH1D*       histoPi0Efficiency_LHC17pLowInt_PHO_5TeV2017       = NULL;
    if( !filePi0PCM5TeV2017_LHC17pLowInt_PHO->IsZombie() ){
        histoPi0InputMCWOWeights_LHC17pLowInt_PHO_5TeV2017      = (TH1D*)filePi0PCM5TeV2017_LHC17pLowInt_PHO->Get("MCYield_Meson_oldBinWOWeights");
        histoPi0Efficiency_LHC17pLowInt_PHO_5TeV2017            = (TH1D*)filePi0PCM5TeV2017_LHC17pLowInt_PHO->Get("TrueMesonEffiPt");
    }

    TString     fileNameEtaPhojet = Form("%s-LHC17pLowInt-%s/5TeV2017/Eta_data_GammaConvV1Correction_%s.root",cutNumber.Data(),nameRecCls.Data(),cutNumber.Data());
	TFile*      fileEtaPCM5TeV2017_LHC17pLowInt_PHO = new TFile(fileNameEtaPhojet.Data());
    TH1D*       histoEtaInputMCWOWeights_LHC17pLowInt_PHO_5TeV2017 = NULL;
    TH1D*       histoEtaEfficiency_LHC17pLowInt_PHO_5TeV2017       = NULL;
    if( !fileEtaPCM5TeV2017_LHC17pLowInt_PHO->IsZombie() ){
        histoEtaInputMCWOWeights_LHC17pLowInt_PHO_5TeV2017      = (TH1D*)fileEtaPCM5TeV2017_LHC17pLowInt_PHO->Get("MCYield_Meson_oldBinWOWeights");
        histoEtaEfficiency_LHC17pLowInt_PHO_5TeV2017            = (TH1D*)fileEtaPCM5TeV2017_LHC17pLowInt_PHO->Get("TrueMesonEffiPt");
    }

    // Jet-Jet
    TString     fileNamePi0JetJet = Form("%s-LHC18b8-%s/5TeV2017/Pi0_data_GammaConvV1Correction_%s.root",cutNumber.Data(),nameRecCls.Data(),cutNumber.Data());
    TFile*      filePi0PCM5TeV2017_LHC18b8_PYT8                = new TFile(fileNamePi0JetJet.Data());
    TH1D*       histoPi0InputMCWOWeights_LHC18b8_PYT8_5TeV2017 = NULL;
    TH1D*       histoPi0Efficiency_LHC18b8_PYT8_5TeV2017       = NULL;
    if( !filePi0PCM5TeV2017_LHC18b8_PYT8->IsZombie() ){
        histoPi0InputMCWOWeights_LHC18b8_PYT8_5TeV2017         = (TH1D*)filePi0PCM5TeV2017_LHC18b8_PYT8->Get("MCYield_Meson_oldBinWOWeights");
        histoPi0Efficiency_LHC18b8_PYT8_5TeV2017               = (TH1D*)filePi0PCM5TeV2017_LHC18b8_PYT8->Get("TrueMesonEffiPt");
    }

    TString fileNameEtaJetJet = Form("%s-LHC18b8-%s/5TeV2017/Eta_data_GammaConvV1Correction_%s.root",cutNumber.Data(),nameRecCls.Data(),cutNumber.Data());
    TFile*      fileEtaPCM5TeV2017_LHC18b8_PYT8                = new TFile(fileNameEtaJetJet.Data());
    TH1D*       histoEtaInputMCWOWeights_LHC18b8_PYT8_5TeV2017 = NULL;
    TH1D*       histoEtaEfficiency_LHC18b8_PYT8_5TeV2017       = NULL;
    if( !fileEtaPCM5TeV2017_LHC18b8_PYT8->IsZombie() ){
        histoEtaInputMCWOWeights_LHC18b8_PYT8_5TeV2017         = (TH1D*)fileEtaPCM5TeV2017_LHC18b8_PYT8->Get("MCYield_Meson_oldBinWOWeights");
        histoEtaEfficiency_LHC18b8_PYT8_5TeV2017               = (TH1D*)fileEtaPCM5TeV2017_LHC18b8_PYT8->Get("TrueMesonEffiPt");
    }


// 	TString nameFilePi0Fullmergedpp5TeV2017 = Form("%s_FullMerged-%s/5TeV2017/Pi0_MC_GammaConvV1Correction_%s.root",cutNumber.Data(),nameRecCls.Data(),cutNumber.Data());
// 	TString nameFileEtaFullmergedpp5TeV2017 = Form("%s_FullMerged-%s/5TeV2017/Eta_MC_GammaConvV1Correction_%s.root",cutNumber.Data(),nameRecCls.Data(),cutNumber.Data());
//
//
// 	if( Iteration0  ){
// 	  nameFilePi0Fullmergedpp5TeV2017 = Form("%s_FullAdded-%s/5TeV2017/Pi0_MC_GammaConvV1Correction_%s.root",cutNumber.Data(),nameRecCls.Data(),cutNumber.Data());
// 	  nameFileEtaFullmergedpp5TeV2017 = Form("%s_FullAdded-%s/5TeV2017/Eta_MC_GammaConvV1Correction_%s.root",cutNumber.Data(),nameRecCls.Data(),cutNumber.Data());
// 	}



    cout<< "Plotting Pi0 spectrum..."<< endl;
    //*******************************************************************************************************************
    //****************************************Pi0 Spectra 5TeV compared to MC********************************************
    //*******************************************************************************************************************
    TCanvas* canvasPi0Spectra5TeV2017  = new TCanvas("canvasPi0Spectra5TeV2017", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasPi0Spectra5TeV2017,  0.13, 0.01, 0.015, 0.08);
    canvasPi0Spectra5TeV2017->SetLogy();
    canvasPi0Spectra5TeV2017->SetLogx();

    TH2F * histo2DPi0Spectra5TeV2017;
    histo2DPi0Spectra5TeV2017    = new TH2F("histo2DPi0Spectra5TeV2017", "histo2DPi0Spectra5TeV2017",1000, 0.23, 20., 1000, 1e-8, 1e1 );
    SetStyleHistoTH2ForGraphs( histo2DPi0Spectra5TeV2017, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",
                                0.03, 0.04, 0.03, 0.04, 0.83, 1.4);
    histo2DPi0Spectra5TeV2017->GetXaxis()->SetLabelOffset(-0.01);
    histo2DPi0Spectra5TeV2017->GetYaxis()->SetLabelOffset(0.01);
    histo2DPi0Spectra5TeV2017->DrawCopy();

        DrawGammaSetMarker(histoPCMYieldPi05TeV2017, markerStyleSpectrum5TeV2017, markerSizePP5TeV2017, color5TeV2017 , color5TeV2017);
        histoPCMYieldPi05TeV2017->Draw("p,same,e1");

        SetStyleHisto(histoPi0InputFullReweighted5TeV2017, 1., lineStyleMCA, kBlue+1);
        histoPi0InputFullReweighted5TeV2017->Draw("same,hist,c");

        SetStyleHisto(histoPi0InputMCWOWeights_LHC17pq_PYT8_5TeV2017, 1., lineStyleMCA, colorMCPythiaPP5TeV2017);
        histoPi0InputMCWOWeights_LHC17pq_PYT8_5TeV2017->Draw("same,hist,c");

        if(histoPi0InputMCWOWeights_LHC17pLowInt_PHO_5TeV2017){
            SetStyleHisto(histoPi0InputMCWOWeights_LHC17pLowInt_PHO_5TeV2017, 1., lineStyleMCB, colorMCPhojetPP5TeV2017);
            histoPi0InputMCWOWeights_LHC17pLowInt_PHO_5TeV2017->Draw("same,hist,c");
        }

        if(histoPi0InputMCWOWeights_LHC18b8_PYT8_5TeV2017){
            SetStyleHisto(histoPi0InputMCWOWeights_LHC18b8_PYT8_5TeV2017,1.,lineStyleMCJetJet, colorMCPythiaJetJetPP5TeV2017);
            histoPi0InputMCWOWeights_LHC18b8_PYT8_5TeV2017->Draw("same,hist,c");
        }


        TLatex *labelSpectraPi0Label 	= new TLatex(0.55,0.92,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
        SetStyleTLatex( labelSpectraPi0Label, 0.035  ,4);
        labelSpectraPi0Label->Draw();

        TLegend* legendSpectraPi05TeV2017 = new TLegend(0.16,0.09,0.73,0.3);
        legendSpectraPi05TeV2017->SetFillColor(0);
        legendSpectraPi05TeV2017->SetLineColor(0);
        legendSpectraPi05TeV2017->SetTextSize(0.025);
        // 	legendSpectraPi05TeV2017->SetNColumns(2);
        legendSpectraPi05TeV2017->SetMargin(0.2);
        legendSpectraPi05TeV2017->AddEntry(histoPCMYieldPi05TeV2017,collisionSystemPP5TeV2017.Data(),"pf");
        legendSpectraPi05TeV2017->AddEntry(histoPi0InputFullReweighted5TeV2017,"full MC reweighted","l");
        legendSpectraPi05TeV2017->AddEntry(histoPi0InputMCWOWeights_LHC17pq_PYT8_5TeV2017,"Pythia 8, LHC17pq","l");
        if(histoPi0InputMCWOWeights_LHC17pLowInt_PHO_5TeV2017) legendSpectraPi05TeV2017->AddEntry(histoPi0InputMCWOWeights_LHC17pLowInt_PHO_5TeV2017,"Phojet, LHC17pLowInt","l");
        if( histoPi0InputMCWOWeights_LHC18b8_PYT8_5TeV2017 ) legendSpectraPi05TeV2017->AddEntry(histoPi0InputMCWOWeights_LHC18b8_PYT8_5TeV2017,"Pythia 8, LHC18b8, Jet-Jet","l");
        legendSpectraPi05TeV2017->Draw();

    canvasPi0Spectra5TeV2017->Update();
    canvasPi0Spectra5TeV2017->Print(Form("%s/Pi0_Spectra_5TeV2017.%s",outputDir.Data(),suffix.Data()));

    canvasPi0Spectra5TeV2017->cd();
    histo2DPi0Spectra5TeV2017->DrawCopy();

        histoPCMYieldPi05TeV2017->Draw("p,same,e1");
        histoPi0InputFullReweighted5TeV2017->Draw("same,hist,c");
        if(histoPi0InputMCWOWeights_LHC18b8_PYT8_5TeV2017) histoPi0InputMCWOWeights_LHC18b8_PYT8_5TeV2017->Draw("same,hist,c");

        TLegend* legendSpectraPi0WithJetJet5TeV2017 = new TLegend(0.16,0.09,0.73,0.3);
        legendSpectraPi0WithJetJet5TeV2017->SetFillColor(0);
        legendSpectraPi0WithJetJet5TeV2017->SetLineColor(0);
        legendSpectraPi0WithJetJet5TeV2017->SetTextSize(0.025);
        legendSpectraPi0WithJetJet5TeV2017->SetMargin(0.2);
        legendSpectraPi0WithJetJet5TeV2017->AddEntry(histoPCMYieldPi05TeV2017,collisionSystemPP5TeV2017.Data(),"pf");
        legendSpectraPi0WithJetJet5TeV2017->AddEntry(histoPi0InputFullReweighted5TeV2017,"full MC reweighted","l");
        if( histoPi0InputMCWOWeights_LHC18b8_PYT8_5TeV2017) legendSpectraPi0WithJetJet5TeV2017->AddEntry(histoPi0InputMCWOWeights_LHC18b8_PYT8_5TeV2017,"Pythia 8, LHC18b8, Jet-Jet","l");
        legendSpectraPi0WithJetJet5TeV2017->Draw();

        labelSpectraPi0Label->Draw();

    canvasPi0Spectra5TeV2017->Update();
    canvasPi0Spectra5TeV2017->Print(Form("%s/Pi0_Spectra_With_Jet_Jet_5TeV2017.%s",outputDir.Data(),suffix.Data()));


    cout << "Fitting the Pi0 spectrum..." << endl;
    //	**********************************************************************************************************************
    //	****************************************Fit Pi0 Spectra 5TeV and plot *********************************************
    //	**********************************************************************************************************************
    // fit spectrum with Tsallis function
    TF1* fitInvYieldDataPi0Comb5TeV2017 = FitObject("l","fitInvYieldDataPi0Comb5TeV2017","Pi0");
    histoPCMYieldPi05TeV2017->Fit(fitInvYieldDataPi0Comb5TeV2017,"QNRMEI+","",0.3,12.);
    fitInvYieldDataPi0Comb5TeV2017->SetRange(0.1,20.);

    cout << WriteParameterToFile(fitInvYieldDataPi0Comb5TeV2017)<< endl;

    canvasPi0Spectra5TeV2017->cd();
    histo2DPi0Spectra5TeV2017->DrawCopy();

        DrawGammaSetMarker(histoPCMYieldPi05TeV2017, markerStyleSpectrum5TeV2017, markerSizePP5TeV2017, color5TeV2017 , color5TeV2017);
        histoPCMYieldPi05TeV2017->Draw("p,same,e1");

        SetStyleHisto(histoPi0InputFullReweighted5TeV2017, 1., lineStyleMCAddSig, kBlue+1);
        histoPi0InputFullReweighted5TeV2017->Draw("same,hist,c");

        fitInvYieldDataPi0Comb5TeV2017->SetLineColor(color5TeV2017);
        fitInvYieldDataPi0Comb5TeV2017->SetLineStyle(2);
        fitInvYieldDataPi0Comb5TeV2017->Draw("same");

        TLegend* legendSpectraPi0WithFit5TeV2017 = new TLegend(0.16,0.09,0.73,0.3);
        legendSpectraPi0WithFit5TeV2017->SetFillColor(0);
        legendSpectraPi0WithFit5TeV2017->SetLineColor(0);
        legendSpectraPi0WithFit5TeV2017->SetTextSize(0.025);
        legendSpectraPi0WithFit5TeV2017->SetMargin(0.2);
        legendSpectraPi0WithFit5TeV2017->AddEntry(histoPCMYieldPi05TeV2017,collisionSystemPP5TeV2017.Data(),"pf");
        legendSpectraPi0WithFit5TeV2017->AddEntry(histoPi0InputFullReweighted5TeV2017,"full MC reweighted","l");
        legendSpectraPi0WithFit5TeV2017->AddEntry(fitInvYieldDataPi0Comb5TeV2017,"Tsallis fit","l");
        legendSpectraPi0WithFit5TeV2017->Draw();

        labelSpectraPi0Label->Draw();

    canvasPi0Spectra5TeV2017->Print(Form("%s/Pi0_Spectra_WithFit_5TeV2017.%s",outputDir.Data(),suffix.Data()));


    cout << "Calculating ratio fit to Pi0 spectrum..." << endl;
    // **********************************************************************************************************************
    // ******************************* Ratio of data to fit and MC input to fit for pi0 2.76 TeV ****************************
    // **********************************************************************************************************************
    TCanvas* canvasRatioToFit = new TCanvas("canvasRatioToFit","",1550,1200);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToFit,  0.08, 0.015, 0.015, 0.08);
    canvasRatioToFit->SetGridx(0);
    canvasRatioToFit->SetGridy(0);
    canvasRatioToFit->cd();

        TLatex *labelSpectraPi0LabelRatio = new TLatex(0.65,0.9,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
        SetStyleTLatex( labelSpectraPi0LabelRatio, 0.035  ,4);
        labelSpectraPi0LabelRatio->Draw();
        TLatex *labelEnergy5TeV2017Ratio = new TLatex(0.65,0.86,collisionSystemPP5TeV2017.Data());
        SetStyleTLatex( labelEnergy5TeV2017Ratio, 0.035  ,4);
        labelEnergy5TeV2017Ratio->Draw();

        TH1D* histoRatioPi0DatatoFit5TeV2017             = CalculateHistoRatioToFit (histoPCMYieldPi05TeV2017, fitInvYieldDataPi0Comb5TeV2017,kTRUE);
        DrawGammaSetMarker(histoRatioPi0DatatoFit5TeV2017, markerStyleSpectrum5TeV2017, markerSizePP5TeV2017, kBlack , kBlack);

        TH1D* histoRatioPi0MCtoDataFit5TeV2017           = CalculateHistoRatioToFit (histoPi0InputFullReweighted5TeV2017, fitInvYieldDataPi0Comb5TeV2017,kTRUE);
        SetStyleHisto(histoRatioPi0MCtoDataFit5TeV2017, 2, lineStyleMCA, kRed+2 );

        TH1D* histoRatioPi0MCUnweightedtoDataFit5TeV2017 = NULL;
        if (histoPi0InputFull5TeV2017) histoRatioPi0MCUnweightedtoDataFit5TeV2017 = CalculateHistoRatioToFit (histoPi0InputFull5TeV2017, fitInvYieldDataPi0Comb5TeV2017,kTRUE);
        if (histoRatioPi0MCUnweightedtoDataFit5TeV2017) SetStyleHisto(histoRatioPi0MCUnweightedtoDataFit5TeV2017, 2, lineStyleMCB, 807 );

        DrawAutoGammaMesonHistos( histoRatioPi0DatatoFit5TeV2017,
                    "", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
                    kFALSE, 1.5, 0, kTRUE,
                    kTRUE, 0, 2.1,
                    kTRUE, 0.,histoRatioPi0DatatoFit5TeV2017->GetXaxis()->GetBinUpEdge(histoRatioPi0DatatoFit5TeV2017->GetNbinsX()));
        histoRatioPi0DatatoFit5TeV2017->GetYaxis()->SetTitleOffset(0.9);
        histoRatioPi0DatatoFit5TeV2017->Draw("e,p");
        if (runDrawReweighted) histoRatioPi0MCtoDataFit5TeV2017->Draw("same,hist,l");
        if (histoRatioPi0MCUnweightedtoDataFit5TeV2017) histoRatioPi0MCUnweightedtoDataFit5TeV2017->Draw("same,hist,l");

        TLegend* legendRatioPi05TeV2017 = new TLegend(0.11,0.12,0.4,0.30);
        legendRatioPi05TeV2017->SetFillColor(0);
        legendRatioPi05TeV2017->SetLineColor(0);
        legendRatioPi05TeV2017->SetTextSize(0.035);
        legendRatioPi05TeV2017->SetMargin(0.2);
        legendRatioPi05TeV2017->AddEntry(histoRatioPi0DatatoFit5TeV2017,"Data/Tsallis fit to Data (0.3 < #it{p}_{T} < 12 GeV/#it{c})","p");
        if (runDrawReweighted) legendRatioPi05TeV2017->AddEntry(histoRatioPi0MCtoDataFit5TeV2017,Form("MC weighted %s/Tsallis  fit to Data (0.3 < #it{p}_{T} < 12  GeV/#it{c})",stringIterationNumber.Data()),"l");
        if (histoRatioPi0MCUnweightedtoDataFit5TeV2017) legendRatioPi05TeV2017->AddEntry(histoRatioPi0MCUnweightedtoDataFit5TeV2017,"MC/Tsallis  fit to Data (0.3 < #it{p}_{T} < 12  GeV/#it{c})","l");
        legendRatioPi05TeV2017->Draw();
        DrawGammaLines(0., 12. ,1., 1.,0.1);

    canvasRatioToFit->Update();
    canvasRatioToFit->SaveAs(Form("%s/Pi0_RatioToDataFit_5TeV2017.%s",outputDir.Data(),suffix.Data()));


    cout << "Compare efficiency for Pi0..." << endl;
	//	**********************************************************************************************************************
	//	******************************Compare Efficiencies for pi0 at 5TeV for different MC *******************************
	//	**********************************************************************************************************************
	TCanvas* canvasPi0Efficiencies5TeV2017DiffMC 					= new TCanvas("canvasPi0Efficiencies5TeV2017DiffMC", "", 200, 10, 1200, 1100);  // gives the page size
	DrawGammaCanvasSettings( canvasPi0Efficiencies5TeV2017DiffMC,  0.09, 0.01, 0.015, 0.08);
	canvasPi0Efficiencies5TeV2017DiffMC->SetLogy();
// 	canvasPi0Efficiencies5TeV2017DiffMC->SetLogx();

	TH2F * histo2DPi0Effi5TeV2017 = new TH2F("histo2DPi0Effi5TeV2017", "histo2DPi0Effi5TeV2017",1000, 0., 10.5, 1000, 1e-5, 1e-2 );
	SetStyleHistoTH2ForGraphs( histo2DPi0Effi5TeV2017, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{#pi^{0}}",
							   0.03, 0.04, 0.03, 0.04, 0.83, 1.05);
	histo2DPi0Effi5TeV2017->GetYaxis()->SetLabelOffset(0.01);
	histo2DPi0Effi5TeV2017->DrawCopy();

        DrawGammaSetMarker(histoPi0Efficiency_LHC17pq_PYT8_5TeV2017, markerStyleMCA, markerSizePP5TeV2017, colorMCPythiaPP5TeV2017 , colorMCPythiaPP5TeV2017);
        histoPi0Efficiency_LHC17pq_PYT8_5TeV2017->Draw("p,same,e1");

        if(histoPi0Efficiency_LHC17pLowInt_PHO_5TeV2017){
            DrawGammaSetMarker(histoPi0Efficiency_LHC17pLowInt_PHO_5TeV2017, markerStyleMCB, markerSizePP5TeV2017, colorMCPhojetPP5TeV2017 , colorMCPhojetPP5TeV2017);
            histoPi0Efficiency_LHC17pLowInt_PHO_5TeV2017->Draw("p,same,e1");
        }
        if( histoPi0Efficiency_LHC18b8_PYT8_5TeV2017 ){
            DrawGammaSetMarker(histoPi0Efficiency_LHC18b8_PYT8_5TeV2017, markerStyleMCJetJet, markerSizePP5TeV2017, colorMCPythiaJetJetPP5TeV2017 , colorMCPythiaJetJetPP5TeV2017);
            histoPi0Efficiency_LHC18b8_PYT8_5TeV2017->Draw("p,same,e1");
        }

        TLegend* legendEfficiencyPi05TeV2017 = new TLegend(0.56,0.09,0.93,0.3);
        legendEfficiencyPi05TeV2017->SetFillColor(0);
        legendEfficiencyPi05TeV2017->SetLineColor(0);
        legendEfficiencyPi05TeV2017->SetTextSize(0.025);
        legendEfficiencyPi05TeV2017->SetMargin(0.2);
        legendEfficiencyPi05TeV2017->AddEntry(histoPi0Efficiency_LHC17pq_PYT8_5TeV2017,"Pythia 8, LHC17pq","p");
        if(histoPi0Efficiency_LHC17pLowInt_PHO_5TeV2017)legendEfficiencyPi05TeV2017->AddEntry(histoPi0Efficiency_LHC17pLowInt_PHO_5TeV2017,"Phojet, LHC17pLowInt","p");
        if( histoPi0Efficiency_LHC18b8_PYT8_5TeV2017 )legendEfficiencyPi05TeV2017->AddEntry(histoPi0Efficiency_LHC18b8_PYT8_5TeV2017,"Pythia 8, LHC18b8, Jet-Jet","p");
        legendEfficiencyPi05TeV2017->Draw();

	canvasPi0Efficiencies5TeV2017DiffMC->Update();
	canvasPi0Efficiencies5TeV2017DiffMC->Print(Form("%s/Pi0_Efficiency_%s_5TeV2017.%s",outputDir.Data(),nameRecCls.Data(),suffix.Data()));


    cout<< "Plotting Eta spectrum..."<< endl;
    //	**********************************************************************************************************************
    //	****************************************Eta Spectra 5TeV compared to MC********************************************
    //	**********************************************************************************************************************
    TCanvas* canvasEtaSpectra5TeV2017 = new TCanvas("canvasEtaSpectra5TeV2017", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasEtaSpectra5TeV2017,  0.13, 0.01, 0.015, 0.08);
    canvasEtaSpectra5TeV2017->SetLogy();
    canvasEtaSpectra5TeV2017->SetLogx();

    TH2F * histo2DEtaSpectraAll = new TH2F("histo2DEtaSpectraAll", "histo2DEtaSpectraAll",1000, 0.23, 20., 1000, 1e-8, 1e1 );
    SetStyleHistoTH2ForGraphs( histo2DEtaSpectraAll, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",
                                0.03, 0.04, 0.03, 0.04, 0.83, 1.4);
    histo2DEtaSpectraAll->GetXaxis()->SetLabelOffset(-0.01);
    histo2DEtaSpectraAll->GetYaxis()->SetLabelOffset(0.01);
    histo2DEtaSpectraAll->DrawCopy();

        DrawGammaSetMarker(histoPCMYieldEta5TeV2017, markerStyleSpectrum5TeV2017, markerSizePP5TeV2017, color5TeV2017 , color5TeV2017);
        histoPCMYieldEta5TeV2017->Draw("p,same,e1");

        SetStyleHisto(histoEtaInputFullReweighted5TeV2017, 1., lineStyleMCA, kBlue+1);
        histoEtaInputFullReweighted5TeV2017->Draw("same,hist,c");

        SetStyleHisto(histoEtaInputMCWOWeights_LHC17pq_PYT8_5TeV2017, 1., lineStyleMCA, colorMCPythiaPP5TeV2017);
        histoEtaInputMCWOWeights_LHC17pq_PYT8_5TeV2017->Draw("same,hist,c");

        if( histoEtaInputMCWOWeights_LHC17pLowInt_PHO_5TeV2017) {
            SetStyleHisto(histoEtaInputMCWOWeights_LHC17pLowInt_PHO_5TeV2017, 1., lineStyleMCB, colorMCPhojetPP5TeV2017);
            histoEtaInputMCWOWeights_LHC17pLowInt_PHO_5TeV2017->Draw("same,hist,c");
        }
        if( histoEtaInputMCWOWeights_LHC18b8_PYT8_5TeV2017) {
            SetStyleHisto(histoEtaInputMCWOWeights_LHC18b8_PYT8_5TeV2017, 1., lineStyleMCJetJet, colorMCPythiaJetJetPP5TeV2017);
            histoEtaInputMCWOWeights_LHC18b8_PYT8_5TeV2017->Draw("same,hist,c");
        }


        TLatex *labelSpectraEtaLabel 		       	= new TLatex(0.55,0.92,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
        SetStyleTLatex( labelSpectraEtaLabel, 0.035  ,4);
        labelSpectraEtaLabel->Draw();

        TLegend* legendSpectraEta5TeV2017 = new TLegend(0.16,0.09,0.73,0.3);
        legendSpectraEta5TeV2017->SetFillColor(0);
        legendSpectraEta5TeV2017->SetLineColor(0);
        legendSpectraEta5TeV2017->SetTextSize(0.025);
        legendSpectraEta5TeV2017->SetMargin(0.2);
        legendSpectraEta5TeV2017->AddEntry(histoPCMYieldEta5TeV2017,collisionSystemPP5TeV2017.Data(),"pf");
        legendSpectraEta5TeV2017->AddEntry(histoEtaInputFullReweighted5TeV2017,"full MC reweighted","l");
        legendSpectraEta5TeV2017->AddEntry(histoEtaInputMCWOWeights_LHC17pq_PYT8_5TeV2017,"Pythia 8, LHC17pq","l");
        if( histoEtaInputMCWOWeights_LHC17pLowInt_PHO_5TeV2017) legendSpectraEta5TeV2017->AddEntry(histoEtaInputMCWOWeights_LHC17pLowInt_PHO_5TeV2017,"Phojet, LHC17pLowInt","l");
        if( histoEtaInputMCWOWeights_LHC18b8_PYT8_5TeV2017) legendSpectraEta5TeV2017->AddEntry(histoEtaInputMCWOWeights_LHC18b8_PYT8_5TeV2017,"Pythia 8, LHC18b8","l");
        legendSpectraEta5TeV2017->Draw();

    canvasEtaSpectra5TeV2017->Update();
    canvasEtaSpectra5TeV2017->Print(Form("%s/Eta_Spectra_5TeV2017.%s",outputDir.Data(),suffix.Data()));

    canvasEtaSpectra5TeV2017->cd();
    histo2DEtaSpectraAll->DrawCopy();

        histoPCMYieldEta5TeV2017->Draw("p,same,e1");
        histoEtaInputFullReweighted5TeV2017->Draw("same,hist,c");
        if( histoEtaInputMCWOWeights_LHC18b8_PYT8_5TeV2017 ) histoEtaInputMCWOWeights_LHC18b8_PYT8_5TeV2017->Draw("same,hist,c");

        TLegend* legendSpectraEtaWithJetJet5TeV2017 = new TLegend(0.16,0.09,0.73,0.3);
        legendSpectraEtaWithJetJet5TeV2017->SetFillColor(0);
        legendSpectraEtaWithJetJet5TeV2017->SetLineColor(0);
        legendSpectraEtaWithJetJet5TeV2017->SetTextSize(0.025);
        legendSpectraEtaWithJetJet5TeV2017->SetMargin(0.2);
        legendSpectraEtaWithJetJet5TeV2017->AddEntry(histoPCMYieldEta5TeV2017,collisionSystemPP5TeV2017.Data(),"pf");
        legendSpectraEtaWithJetJet5TeV2017->AddEntry(histoEtaInputFullReweighted5TeV2017,"full MC reweighted","l");
        if( histoEtaInputMCWOWeights_LHC18b8_PYT8_5TeV2017) legendSpectraEtaWithJetJet5TeV2017->AddEntry(histoEtaInputMCWOWeights_LHC18b8_PYT8_5TeV2017,"Pythia 8, LHC18b8, Jet-Jet","l");
        legendSpectraEtaWithJetJet5TeV2017->Draw();

        labelSpectraEtaLabel->Draw();

    canvasEtaSpectra5TeV2017->Update();
    canvasEtaSpectra5TeV2017->Print(Form("%s/Eta_Spectra_With_Jet_Jet_5TeV2017.%s",outputDir.Data(),suffix.Data()));


    cout << "Fitting the Eta spectrum..." << endl;
    //	**********************************************************************************************************************
    //	****************************************Fit Pi0 Spectra 5TeV and plot *********************************************
    //	**********************************************************************************************************************
    // fit spectrum with Tsallis function - "n" fixed by result from pi0 at same energy
    TF1* fitInvYieldDataEtaComb5TeV2017 = FitObject("l","fitInvYieldDataEtaComb5TeV2017","Eta");
    // 	fitInvYieldDataEtaComb5TeV2017->FixParameter(1,fitInvYieldDataPi0Comb5TeV2017->GetParameter(1));
    // 	fitInvYieldDataEtaComb5TeV2017->SetParameter(1,fitInvYieldDataPi0Comb5TeV2017->GetParameter(1));
    // 	fitInvYieldDataEtaComb5TeV2017->SetParLimits(1,fitInvYieldDataPi0Comb5TeV2017->GetParameter(1)*0.9,fitInvYieldDataPi0Comb5TeV2017->GetParameter(1)*1.1);
    histoPCMYieldEta5TeV2017->Fit(fitInvYieldDataEtaComb5TeV2017,"QNRMEI+","",0.4,12.);
    fitInvYieldDataEtaComb5TeV2017->SetRange(0.1,20.);
    cout << WriteParameterToFile(fitInvYieldDataEtaComb5TeV2017)<< endl;

    canvasEtaSpectra5TeV2017->cd();
    histo2DEtaSpectraAll->DrawCopy();

        DrawGammaSetMarker(histoPCMYieldEta5TeV2017, markerStyleSpectrum5TeV2017, markerSizePP5TeV2017, color5TeV2017 , color5TeV2017);
        histoPCMYieldEta5TeV2017->Draw("p,same,e1");

        SetStyleHisto(histoEtaInputFullReweighted5TeV2017, 1., lineStyleMCAddSig, kBlue+1);
        histoEtaInputFullReweighted5TeV2017->Draw("same,hist,c");

        fitInvYieldDataEtaComb5TeV2017->SetLineColor(color5TeV2017);
        fitInvYieldDataEtaComb5TeV2017->SetLineStyle(2);
        fitInvYieldDataEtaComb5TeV2017->Draw("same");

        TLegend* legendSpectraEtaWithFit5TeV2017 = new TLegend(0.16,0.09,0.73,0.3);
        legendSpectraEtaWithFit5TeV2017->SetFillColor(0);
        legendSpectraEtaWithFit5TeV2017->SetLineColor(0);
        legendSpectraEtaWithFit5TeV2017->SetTextSize(0.025);
        legendSpectraEtaWithFit5TeV2017->SetMargin(0.2);
        legendSpectraEtaWithFit5TeV2017->AddEntry(histoPCMYieldEta5TeV2017,collisionSystemPP5TeV2017.Data(),"pf");
        legendSpectraEtaWithFit5TeV2017->AddEntry(histoEtaInputFullReweighted5TeV2017,"full MC reweighted","l");
        legendSpectraEtaWithFit5TeV2017->AddEntry(fitInvYieldDataEtaComb5TeV2017,"Tsallis fit","l");
        legendSpectraEtaWithFit5TeV2017->Draw();
        labelSpectraEtaLabel->Draw();

    canvasEtaSpectra5TeV2017->Print(Form("%s/Eta_Spectra_WithFit_5TeV2017.%s",outputDir.Data(),suffix.Data()));


    cout << "Calculating ratio fit to Eta spectrum..." << endl;
	// **********************************************************************************************************************
    // ******************************* Ratio of data to fit and MC input to fit for Eta 2.76 TeV ****************************
	// **********************************************************************************************************************

	TH1D* histoRatioEtaDatatoFit5TeV2017 		= CalculateHistoRatioToFit (histoPCMYieldEta5TeV2017, fitInvYieldDataEtaComb5TeV2017,kTRUE);
	DrawGammaSetMarker(histoRatioEtaDatatoFit5TeV2017, markerStyleSpectrum5TeV2017, markerSizePP5TeV2017, kBlack , kBlack);

	TH1D* histoRatioEtaMCtoDataFit5TeV2017 		= CalculateHistoRatioToFit (histoEtaInputFullReweighted5TeV2017, fitInvYieldDataEtaComb5TeV2017,kTRUE);
	SetStyleHisto(histoRatioEtaMCtoDataFit5TeV2017, 2, lineStyleMCA, kRed+2 );

	TH1D* histoRatioEtaMCUnweightedtoDataFit5TeV2017 = NULL;
	if (histoEtaInputFull5TeV2017) histoRatioEtaMCUnweightedtoDataFit5TeV2017 = CalculateHistoRatioToFit (histoEtaInputFull5TeV2017, fitInvYieldDataEtaComb5TeV2017,kTRUE);
	if (histoRatioEtaMCUnweightedtoDataFit5TeV2017) SetStyleHisto(histoRatioEtaMCUnweightedtoDataFit5TeV2017, 2, lineStyleMCB, 807 );

	canvasRatioToFit->cd();

        TLatex *labelSpectraEtaLabelRatio 									= new TLatex(0.65,0.9,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
        SetStyleTLatex( labelSpectraEtaLabelRatio, 0.035  ,4);
        labelSpectraEtaLabelRatio->Draw();
        labelEnergy5TeV2017Ratio->Draw();

        DrawAutoGammaMesonHistos( histoRatioEtaDatatoFit5TeV2017,
                    "", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
                    kFALSE, 1.5, 0, kTRUE,
                    kTRUE, 0, 2.1,
                    kTRUE, 0.,histoRatioEtaDatatoFit5TeV2017->GetXaxis()->GetBinUpEdge(histoRatioEtaDatatoFit5TeV2017->GetNbinsX()));
        histoRatioEtaDatatoFit5TeV2017->GetYaxis()->SetTitleOffset(0.9);
        histoRatioEtaDatatoFit5TeV2017->Draw("e,p");
        if (runDrawReweighted) histoRatioEtaMCtoDataFit5TeV2017->Draw("same,hist,l");
        if (histoRatioEtaMCUnweightedtoDataFit5TeV2017) histoRatioEtaMCUnweightedtoDataFit5TeV2017->Draw("same,hist,l");

        TLegend* legendRatioEta5TeV2017 = new TLegend(0.11,0.12,0.4,0.30);
        legendRatioEta5TeV2017->SetFillStyle(0);
        legendRatioEta5TeV2017->SetFillColor(0);
        legendRatioEta5TeV2017->SetLineColor(0);
        legendRatioEta5TeV2017->SetTextSize(0.035);
        legendRatioEta5TeV2017->SetMargin(0.2);
        legendRatioEta5TeV2017->AddEntry(histoRatioEtaDatatoFit5TeV2017,"Data/Tsallis fit to Data (0.4 < #it{p}_{T} < 12 GeV/#it{c})","p");
        if (runDrawReweighted) legendRatioEta5TeV2017->AddEntry(histoRatioEtaMCtoDataFit5TeV2017,Form("MC weighted %s/Tsallis  fit to Data (0.4 < #it{p}_{T} < 12 GeV/#it{c})",stringIterationNumber.Data()),"l");
        if (histoRatioEtaMCUnweightedtoDataFit5TeV2017) legendRatioEta5TeV2017->AddEntry(histoRatioEtaMCUnweightedtoDataFit5TeV2017,"MC/Tsallis  fit to Data (0.4 < #it{p}_{T} < 12 GeV/#it{c})","l");
        legendRatioEta5TeV2017->Draw();
        DrawGammaLines(0.,12.,1., 1.,0.1);

	canvasRatioToFit->Update();
	canvasRatioToFit->SaveAs(Form("%s/Eta_RatioToDataFit_5TeV2017.%s",outputDir.Data(),suffix.Data()));


    cout << "Compare efficiency for Eta..." << endl;
	//	**********************************************************************************************************************
	//	******************************Compare Efficiencies for Eta at 5TeV for different MC *******************************
	//	**********************************************************************************************************************
	TCanvas* canvasEtaEfficiencies5TeV2017DiffMC 					= new TCanvas("canvasEtaEfficiencies5TeV2017DiffMC", "", 200, 10, 1200, 1100);  // gives the page size
	DrawGammaCanvasSettings( canvasEtaEfficiencies5TeV2017DiffMC,  0.09, 0.01, 0.015, 0.08);
	canvasEtaEfficiencies5TeV2017DiffMC->SetLogy();
// 	canvasEtaEfficiencies5TeV2017DiffMC->SetLogx();

	TH2F * histo2DEtaEffi5TeV2017 = new TH2F("histo2DEtaEffi5TeV2017", "histo2DEtaEffi5TeV2017",1000, 0., 6.5, 1000, 1e-5, 1e-2 );
	SetStyleHistoTH2ForGraphs( histo2DEtaEffi5TeV2017, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{#eta}",
							   0.03, 0.04, 0.03, 0.04, 0.83, 1.05);
	histo2DEtaEffi5TeV2017->GetYaxis()->SetLabelOffset(0.01);
	histo2DEtaEffi5TeV2017->DrawCopy();

        DrawGammaSetMarker(histoEtaEfficiency_LHC17pq_PYT8_5TeV2017, markerStyleMCA, markerSizePP5TeV2017, colorMCPythiaPP5TeV2017 , colorMCPythiaPP5TeV2017);
        histoEtaEfficiency_LHC17pq_PYT8_5TeV2017->Draw("p,same,e1");

        if( histoEtaEfficiency_LHC17pLowInt_PHO_5TeV2017 ){
            DrawGammaSetMarker(histoEtaEfficiency_LHC17pLowInt_PHO_5TeV2017, markerStyleMCB, markerSizePP5TeV2017, colorMCPhojetPP5TeV2017 , colorMCPhojetPP5TeV2017);
            histoEtaEfficiency_LHC17pLowInt_PHO_5TeV2017->Draw("p,same,e1");
        }

        if( histoEtaEfficiency_LHC18b8_PYT8_5TeV2017 ){
            DrawGammaSetMarker(histoEtaEfficiency_LHC18b8_PYT8_5TeV2017, markerStyleMCJetJet, markerSizePP5TeV2017, colorMCPythiaJetJetPP5TeV2017 , colorMCPythiaJetJetPP5TeV2017);
            histoEtaEfficiency_LHC18b8_PYT8_5TeV2017->Draw("p,same,e1");
        }

        TLegend* legendEfficiencyEta5TeV2017 										= new TLegend(0.56,0.09,0.93,0.3);
        legendEfficiencyEta5TeV2017->SetFillColor(0);
        legendEfficiencyEta5TeV2017->SetLineColor(0);
        legendEfficiencyEta5TeV2017->SetTextSize(0.025);
        legendEfficiencyEta5TeV2017->SetMargin(0.2);
        legendEfficiencyEta5TeV2017->AddEntry(histoEtaEfficiency_LHC17pq_PYT8_5TeV2017,"Pythia 8, LHC17pq","p");
        if(histoEtaEfficiency_LHC17pLowInt_PHO_5TeV2017) legendEfficiencyEta5TeV2017->AddEntry(histoEtaEfficiency_LHC17pLowInt_PHO_5TeV2017,"Phojet, LHC17pLowInt","p");
        if(histoEtaEfficiency_LHC18b8_PYT8_5TeV2017) legendEfficiencyEta5TeV2017->AddEntry(histoEtaEfficiency_LHC18b8_PYT8_5TeV2017,"Pythia 8, LHC18b8, Jet-Jet","p");
        legendEfficiencyEta5TeV2017->Draw();

	canvasEtaEfficiencies5TeV2017DiffMC->Update();
	canvasEtaEfficiencies5TeV2017DiffMC->Print(Form("%s/Eta_Efficiency_%s_5TeV2017.%s",outputDir.Data(),nameRecCls.Data(),suffix.Data()));


	//	**********************************************************************************************************************
	//	****************************************Write fits & input MC spectra to file ****************************************
	//	**********************************************************************************************************************
    TFile *fMCSpectraInput = new TFile(Form("%s/MCSpectraInputpp.root",outputDir.Data()),"UPDATE");
    fMCSpectraInput->cd();

		if (fitInvYieldDataPi0Comb5TeV2017){
			fitInvYieldDataPi0Comb5TeV2017->SetRange(0,30);
			fitInvYieldDataPi0Comb5TeV2017->Write("Pi0_Fit_Data_5TeV2017",TObject::kOverwrite);
		}
		if (fitInvYieldDataEtaComb5TeV2017){
			fitInvYieldDataEtaComb5TeV2017->SetRange(0,30);
			fitInvYieldDataEtaComb5TeV2017->Write("Eta_Fit_Data_5TeV2017",TObject::kOverwrite);
		}
		if (histoPi0InputMCWOWeights_LHC17pq_PYT8_5TeV2017){
			histoPi0InputMCWOWeights_LHC17pq_PYT8_5TeV2017->SetTitle(Form("Pi0_Pythia8_LHC17pq_%s_5TeV2017",nameRecCls.Data()));
			histoPi0InputMCWOWeights_LHC17pq_PYT8_5TeV2017->Write(Form("Pi0_Pythia8_LHC17pq_%s_5TeV2017",nameRecCls.Data()),TObject::kOverwrite);
		}
		if (histoEtaInputMCWOWeights_LHC17pq_PYT8_5TeV2017){
			histoEtaInputMCWOWeights_LHC17pq_PYT8_5TeV2017->SetTitle(Form("Eta_Pythia8_LHC17pq_%s_5TeV2017",nameRecCls.Data()));
			histoEtaInputMCWOWeights_LHC17pq_PYT8_5TeV2017->Write(Form("Eta_Pythia8_LHC17pq_%s_5TeV2017",nameRecCls.Data()),TObject::kOverwrite);
		}
		if (histoPi0InputMCWOWeights_LHC17pLowInt_PHO_5TeV2017){
			histoPi0InputMCWOWeights_LHC17pLowInt_PHO_5TeV2017->SetTitle(Form("Pi0_Phojet_LHC17pLowInt_%s_5TeV2017",nameRecCls.Data()));
			histoPi0InputMCWOWeights_LHC17pLowInt_PHO_5TeV2017->Write(Form("Pi0_Phojet_LHC17pLowInt_%s_5TeV2017",nameRecCls.Data()),TObject::kOverwrite);
		}
		if (histoEtaInputMCWOWeights_LHC17pLowInt_PHO_5TeV2017){
			histoEtaInputMCWOWeights_LHC17pLowInt_PHO_5TeV2017->SetTitle(Form("Eta_Phojet_LHC17pLowInt_%s_5TeV2017",nameRecCls.Data()));
			histoEtaInputMCWOWeights_LHC17pLowInt_PHO_5TeV2017->Write(Form("Eta_Phojet_LHC17pLowInt_%s_5TeV2017",nameRecCls.Data()),TObject::kOverwrite);
		}

	fMCSpectraInput->Close();

}
