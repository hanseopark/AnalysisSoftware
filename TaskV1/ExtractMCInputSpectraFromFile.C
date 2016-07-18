// provided by Gamma Conversion Group, $ALICE_PHYSICS/PWGGA/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion
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
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

//****************************************************************************
//*********************** Function to scale the MC input yields **************
//*********************** with #Delta Y, 2 #pi and nEvt **********************
//****************************************************************************
void ScaleMCYield( TH1D* histoCorrectedToBeScaled, Double_t deltaRapid, Double_t scaling, Double_t nEvtMC){
    histoCorrectedToBeScaled->Sumw2();
    histoCorrectedToBeScaled->Scale(1./deltaRapid);
    histoCorrectedToBeScaled->Scale(scaling);
    histoCorrectedToBeScaled->Scale(1./nEvtMC);
    for (Int_t i = 1; i < histoCorrectedToBeScaled->GetNbinsX()+1 ; i++){
        Double_t newBinContent = histoCorrectedToBeScaled->GetBinContent(i)/histoCorrectedToBeScaled->GetBinCenter(i)/histoCorrectedToBeScaled->GetBinWidth(i);
        Double_t newBinError = histoCorrectedToBeScaled->GetBinError(i)/histoCorrectedToBeScaled->GetBinCenter(i)/histoCorrectedToBeScaled->GetBinWidth(i);
        histoCorrectedToBeScaled->SetBinContent(i,newBinContent);
        histoCorrectedToBeScaled->SetBinError(i,newBinError);
    }
}
           
           
           
//****************************************************************************
//********* Main function for extraction of MC inputs ************************
//****************************************************************************
void ExtractMCInputSpectraFromFile( TString file                    = "", 
                                    TString cutSelection            = "", 
                                    TString suffix                  = "", 
                                    TString optionEnergy            = "", 
                                    TString optionPeriod            = "",
                                    Int_t mode                      = 10,
                                    TString nameExternalInput       = "ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_2016_07_04.root"
                                ) {
    
    //*****************************************************************************************************************
    //******************************** Default plotting settings ******************************************************
    //*****************************************************************************************************************
    gROOT->Reset();
    StyleSettingsThesis(suffix);
    SetPlotStyle();

    //*****************************************************************************************************************
    //******************************* Set global variables ************************************************************
    //*****************************************************************************************************************
    Int_t fMode                        = mode;
    // mode:   10 // new output merged EMCal
    //         11 // new output merged PHOS
    if (!(fMode == 10 || fMode == 11)){
        cout << "ERROR: You are running the wrong macro, this macro is only designed for merged cluster analysis." << endl;
        return;
    }
    
    TString fEventCutSelection          = "";
    TString fClusterCutSelection        = "";
    TString fClusterMergedCutSelection  = "";
    TString fMesonCutSelection          = "";
    TString fCutSelection               = cutSelection;
    TString dummyString                 = "";
    ReturnSeparatedCutNumberAdvanced( cutSelection, fEventCutSelection, fClusterCutSelection, fClusterMergedCutSelection, dummyString, fMesonCutSelection, mode);

    TString rapidityRange               = "";
    Double_t deltaRapid                 = ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
    
    //Variable defintion
    Double_t scaling                    = 1./(2.*TMath::Pi());

    

    TString outputDir               = Form("%s/ExtractMCInputSpectraFromFile",suffix.Data());
    gSystem->Exec("mkdir -p "+outputDir);
    
    cout<<"Pictures are saved as "<< suffix.Data()<< endl;
    TString fdate                   = ReturnDateString();
    
    //****************************** Specification of collision system ************************************************
    TString fCollisionSystem        = ReturnFullCollisionsSystem(optionEnergy);
    TString fCollisionSystenWrite   = ReturnCollisionEnergyOutputString(optionEnergy);
    TString fGeneratorName          = ReturnGeneratorNameFromMCName(optionPeriod);
    
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
            
    if(cutSelection.Length() == 0){
        cout<<"ERROR: Cut selection is not set, please do!"<<endl;
        return;
    }
    Int_t kCollisionSystem                      = 0; // 0 : pp, 1: PbPb, 2: pPb
    if (optionEnergy.CompareTo("PbPb_2.76TeV") == 0) 
        kCollisionSystem                        = 1;
    if (optionEnergy.CompareTo("pPb_5.023TeV") == 0) 
        kCollisionSystem                        = 2;

    
    //**************************** Determine Centrality *************************************************************
    TString fCentralityString       = GetCentralityString(fEventCutSelection);
    TString fTextCent               = "";
    if (fCentralityString.CompareTo("pp")==0){
        fTextCent                   = "MinBias";
    } else {
        fTextCent                   = Form("%s central", fCentralityString.Data());
    }
    if (fCentralityString.CompareTo("pp")!=0 && !fCentralityString.Contains("0-100%") ){
        fCollisionSystem            = Form("%s %s", fCentralityString.Data(), fCollisionSystem.Data());
    }
                
    //***************************** Load file ***********************************************************************
    TFile f(file.Data());
    
    TString nameMainDir             = "";
    if (mode == 10 || mode == 11){
        nameMainDir                 = "GammaCaloMerged";
    } else {
        cout << "ERROR: Wrong mode aborting here!" << endl;
        return;
    }
    
    TList *TopDir                   = (TList*)f.Get(nameMainDir.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }
    
    TList *HistosGammaConversion    = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelection.Data()));
    if(HistosGammaConversion == NULL){
        cout<<"ERROR: " << Form("Cut Number %s",fCutSelection.Data()) << " not Found in File"<<endl;
        return;
    }

    // ------------------------------- Determine event number -----------------------------------------------
    TList *ESDContainer                 = (TList*) HistosGammaConversion->FindObject(Form("%s ESD histograms",fCutSelection.Data()));
    TH1D* fEventQuality                 = (TH1D*)ESDContainer->FindObject("NEvents");
    Float_t nEvtMC                      = 0;
    if (kCollisionSystem == 1){
        nEvtMC                          = fEventQuality->GetBinContent(1);
    } else {
        nEvtMC                          = GetNEvents(fEventQuality);
        // BinContent 5 - Zvertex-position, BinContent 4 - no Trigger Bit, BinContent 7 - PileUp 
    }
        
    // ------------------------------- Read MC information ---------------------------------------------------    
    TList *MCContainer                  = (TList*)HistosGammaConversion->FindObject(Form("%s MC histograms",fCutSelection.Data()));
    // ----------- read pi0 -----------------------------------
    cout << "reading pi0" << endl;
    TH1D* fHistoMCPi0GGPt               = (TH1D*)MCContainer->FindObject("MC_Pi0_Pt");  
    TH1D* fHistoMCPi0DalitzPt           = (TH1D*)MCContainer->FindObject("MC_Pi0Dalitz_Pt");  
    TH1D* fHistoMCPi0Pt                 = (TH1D*)fHistoMCPi0GGPt->Clone("MC_Pi0_All_Pt");
    fHistoMCPi0Pt->Sumw2();
    
    if (fHistoMCPi0DalitzPt){
        fHistoMCPi0Pt->Add(fHistoMCPi0DalitzPt);
        fHistoMCPi0Pt->Scale(1/(0.98923+0.01174));
    } else {
        fHistoMCPi0Pt->Scale(1/0.98923);
    }    
    // ----------- read eta -----------------------------------
    cout << "reading eta" << endl;
    TH1D* fHistoMCEtaGGPt               = (TH1D*)MCContainer->FindObject("MC_Eta_Pt");  
    TH1D* fHistoMCEtaDalitzPt           = (TH1D*)MCContainer->FindObject("MC_EtaDalitz_Pt");  
    TH1D* fHistoMCEtaPt                 = (TH1D*)fHistoMCEtaGGPt->Clone("MC_Eta_All_Pt");
    fHistoMCEtaPt->Sumw2();
    if (fHistoMCEtaDalitzPt){
        fHistoMCEtaPt->Add(fHistoMCEtaDalitzPt);
        fHistoMCEtaPt->Scale(1/(0.3941+0.000069));
    } else {
        fHistoMCEtaPt->Scale(1/0.3941);
    }    
    // ----------- charged pions ------------------------------
    cout << "reading charged pions" << endl;
    TH1D* fHistoMCPiNegPt               = (TH1D*)MCContainer->FindObject("MC_NegPi_Pt");  
    TH1D* fHistoMCPiPosPt               = (TH1D*)MCContainer->FindObject("MC_PosPi_Pt");
    TH1D* fHistoMCPiPt                  = (TH1D*)fHistoMCPiNegPt->Clone("MC_PiCh_All_Pt");
    fHistoMCPiPt->Sumw2();
    fHistoMCPiPt->Add(fHistoMCPiPosPt);
    fHistoMCPiPt->Scale(0.5);
    // ----------- charged kaons ------------------------------
    cout << "reading charged kaons" << endl;
    TH1D* fHistoMCKNegPt                = (TH1D*)MCContainer->FindObject("MC_NegK_Pt");  
    TH1D* fHistoMCKPosPt                = (TH1D*)MCContainer->FindObject("MC_PosK_Pt");  
    TH1D* fHistoMCKPt                   = (TH1D*)fHistoMCKNegPt->Clone("MC_KCh_All_Pt");
    fHistoMCKPt->Sumw2();
    fHistoMCKPt->Add(fHistoMCKPosPt);
    fHistoMCKPt->Scale(0.5);
    // ----------- neutral kaons ------------------------------
    cout << "reading K0s" << endl;
    TH1D* fHistoMCK0sPt                 = (TH1D*)MCContainer->FindObject("MC_K0s_Pt");  

    ScaleMCYield(fHistoMCPi0Pt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCEtaPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCPiPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCPiNegPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCPiPosPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCKPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCKNegPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCKPosPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCK0sPt,  deltaRapid,  scaling,  nEvtMC );
    
    TH1D* fHistoRatioMCPiDivDataFit     = NULL;
    TH1D* fHistoRatioMCPi0DivDataFit    = NULL;
    TH1D* fHistoRatioMCKDivDataFit      = NULL;
    TH1D* fHistoRatioMCK0sDivDataFit    = NULL;
    TH1D* fHistoRatioDataPiDivDataFit   = NULL;
    TH1D* fHistoRatioDataKDivDataFit    = NULL;
    
    if (nameExternalInput.CompareTo("") != 0){
      TFile* fileDataInput          = new TFile(nameExternalInput.Data());
      TH1D* fHistoChargedPionData   = NULL;
      TH1D* fHistoChargedKaonData   = NULL;
      if (optionEnergy.CompareTo("2.76TeV") == 0){
        fHistoChargedPionData       = (TH1D*)fileDataInput->Get("histoChargedPionSpecPubStat2760GeV");
        TF1* fitChargedPions        = FitObject("l","fitChargedPions","Pi0",fHistoChargedPionData,0.1,20.,NULL,"QNRMEI");
        fHistoRatioDataPiDivDataFit = CalculateHistoRatioToFit (fHistoChargedPionData, fitChargedPions); 
        fHistoRatioMCPiDivDataFit   = CalculateHistoRatioToFit (fHistoMCPiPt, fitChargedPions); 
        fHistoRatioMCPi0DivDataFit  = CalculateHistoRatioToFit (fHistoMCPi0Pt, fitChargedPions); 
        
        fHistoChargedKaonData       = (TH1D*)fileDataInput->Get("histoChargedKaonSpecPubStat2760GeV");
        TF1* fitChargedKaons        =  FitObject("l","ptDistribution","K",fHistoChargedKaonData,0.1,20.,NULL,"QNRMEI");
        
        fHistoRatioDataKDivDataFit  = CalculateHistoRatioToFit (fHistoChargedKaonData, fitChargedKaons); 
        fHistoRatioMCK0sDivDataFit  = CalculateHistoRatioToFit (fHistoMCK0sPt, fitChargedKaons); 
        fHistoRatioMCKDivDataFit    = CalculateHistoRatioToFit (fHistoMCKPt, fitChargedKaons); 
        
        
        
        TCanvas *canvasFitQA = new TCanvas("canvasFitQA","canvasFitQA",1000,800);
        DrawGammaCanvasSettings( canvasFitQA, 0.12, 0.02, 0.02, 0.08);
        canvasFitQA->cd();
        canvasFitQA->SetLogy(1);

        DrawAutoGammaMesonHistos(   fHistoChargedPionData, 
                            "", "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 
                            kTRUE, 10, 1e-10, kFALSE,
                            kFALSE, 0., 0.7, 
                            kFALSE, 0., 10.);
        fHistoChargedPionData->GetYaxis()->SetTitleOffset(1.3);
        DrawGammaSetMarker(fHistoChargedPionData, 20, 1.5, kAzure-6, kAzure-6);
        fHistoChargedPionData->DrawClone("pe");

        DrawGammaSetMarker(fHistoMCPiPt, 20, 1.5, kRed+2, kRed+2);
        fHistoMCPiPt->Draw("same,pe");
        DrawGammaSetMarker(fHistoMCPi0Pt, 24, 1.5, kGreen+2, kGreen+2);
        fHistoMCPi0Pt->Draw("same,pe");

        fitChargedPions->SetLineColor(kAzure-6);
        fitChargedPions->Draw("same");
        fHistoChargedPionData->Draw("same,pe");

        TLegend* legendSpectraPi = GetAndSetLegend2(0.73, 0.70, 0.95, 0.95, 32,1); 
        legendSpectraPi->AddEntry(fHistoChargedPionData,"Data: #frac{#pi^{+}+#pi^{-}}{2}","p");
        legendSpectraPi->AddEntry(fitChargedPions,"Data: fit","l");
        legendSpectraPi->AddEntry(fHistoMCPiPt,"MC: #frac{#pi^{+}+#pi^{-}}{2}","p");
        legendSpectraPi->AddEntry(fHistoMCPi0Pt,"MC: #pi^{0}","p");
        legendSpectraPi->Draw();
        
        TLatex *labelEnergySpectra        = new TLatex(0.15,0.16,Form("%s",fCollisionSystem.Data()));
        SetStyleTLatex( labelEnergySpectra, 0.04,4);
        labelEnergySpectra->Draw();
        TLatex *labelGeneratorSpectra     = new TLatex(0.15,0.12,Form("%s",fGeneratorName.Data()));
        SetStyleTLatex( labelGeneratorSpectra, 0.04,4);
        labelGeneratorSpectra->Draw();

        
        canvasFitQA->SaveAs(Form("%s/Pi_ComparisonMCAndData_%s_%s.%s",outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), suffix.Data()));
        
        DrawAutoGammaMesonHistos(   fHistoChargedKaonData, 
                            "", "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 
                            kTRUE, 10, 1e-10, kFALSE,
                            kFALSE, 0., 0.7, 
                            kFALSE, 0., 10.);
        fHistoChargedKaonData->GetYaxis()->SetTitleOffset(1.3);
        DrawGammaSetMarker(fHistoChargedKaonData, 20, 1.5, kAzure-6, kAzure-6);
        fHistoChargedKaonData->DrawClone("pe");

        DrawGammaSetMarker(fHistoMCKPt, 20, 1.5, kRed+2, kRed+2);
        fHistoMCKPt->Draw("same,pe");
        DrawGammaSetMarker(fHistoMCK0sPt, 24, 1.5, kGreen+2, kGreen+2);
        fHistoMCK0sPt->Draw("same,pe");
        fitChargedKaons->SetLineColor(kAzure-6);
        fitChargedKaons->Draw("same");
        fHistoChargedKaonData->Draw("same,pe");

        TLegend* legendSpectraK = GetAndSetLegend2(0.73, 0.70, 0.95, 0.95, 32,1); 
        legendSpectraK->AddEntry(fHistoChargedKaonData,"Data: #frac{K^{+}+K^{-}}{2}","p");
        legendSpectraK->AddEntry(fitChargedKaons,"Data: fit","l");
        legendSpectraK->AddEntry(fHistoMCKPt,"MC: #frac{K^{+}+K^{-}}{2}","p");
        legendSpectraK->AddEntry(fHistoMCK0sPt,"MC: K^{0}_{s}","p");
        legendSpectraK->Draw();
        labelEnergySpectra->Draw();
        labelGeneratorSpectra->Draw();
        
        canvasFitQA->SaveAs(Form("%s/K_ComparisonMCAndData_%s_%s.%s",outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), suffix.Data()));
    
        TCanvas *canvasRatio = new TCanvas("canvasRatio","canvasRatio",1000,800);
        DrawGammaCanvasSettings( canvasRatio, 0.07, 0.02, 0.02, 0.08);
        canvasRatio->cd();
        canvasRatio->SetLogy(0);
        
        DrawAutoGammaMesonHistos(   fHistoRatioDataPiDivDataFit, 
                            "", "#it{p}_{T} (GeV/#it{c})", "Data/Fit", 
                            kFALSE, 10, 1e-10, kFALSE,
                            kTRUE, 0., 2, 
                            kTRUE, 0., 20.);
        fHistoRatioDataPiDivDataFit->GetYaxis()->SetTitleOffset(0.85);
        DrawGammaSetMarker(fHistoRatioDataPiDivDataFit, 20, 1.5, kAzure-6, kAzure-6);
        fHistoRatioDataPiDivDataFit->DrawClone("pe");
        DrawGammaSetMarker(fHistoRatioMCPiDivDataFit, 20, 1.5, kRed+2, kRed+2);
        fHistoRatioMCPiDivDataFit->Draw("same,pe");
        DrawGammaSetMarker(fHistoRatioMCPi0DivDataFit, 24, 1.5, kGreen+2, kGreen+2);
        fHistoRatioMCPi0DivDataFit->Draw("same,pe");


        TLegend* legendRatioPi = GetAndSetLegend2(0.12, 0.70, 0.32, 0.95, 32,1); 
        legendRatioPi->AddEntry(fHistoRatioDataPiDivDataFit,"Data: #frac{#pi^{+}+#pi^{-}}{2}","p");
        legendRatioPi->AddEntry(fHistoRatioMCPiDivDataFit,"MC: #frac{#pi^{+}+#pi^{-}}{2}","p");
        legendRatioPi->AddEntry(fHistoRatioMCPi0DivDataFit,"MC: #pi^{0}","p");
        legendRatioPi->Draw();
        
        DrawGammaLines(0., 20,1, 1,0.1, kGray+2, 7);
        TLatex *labelEnergyRatio        = new TLatex(0.11,0.16,Form("%s",fCollisionSystem.Data()));
        SetStyleTLatex( labelEnergyRatio, 0.04,4);
        labelEnergyRatio->Draw();
        TLatex *labelGeneratorRatio     = new TLatex(0.11,0.12,Form("%s",fGeneratorName.Data()));
        SetStyleTLatex( labelGeneratorRatio, 0.04,4);
        labelGeneratorRatio->Draw();
        
        canvasRatio->SaveAs(Form("%s/Pi_RatioComparisonMCAndData_%s_%s.%s",outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), suffix.Data()));

        canvasRatio->cd();
        
        DrawAutoGammaMesonHistos(   fHistoRatioDataKDivDataFit, 
                            "", "#it{p}_{T} (GeV/#it{c})", "Data/Fit", 
                            kFALSE, 10, 1e-10, kFALSE,
                            kTRUE, 0., 2, 
                            kTRUE, 0., 20.);
        fHistoRatioDataKDivDataFit->GetYaxis()->SetTitleOffset(0.85);
        DrawGammaSetMarker(fHistoRatioDataKDivDataFit, 20, 1.5, kAzure-6, kAzure-6);
        fHistoRatioDataKDivDataFit->DrawClone("pe");
        
        
        DrawGammaSetMarker(fHistoRatioMCKDivDataFit, 20, 1.5, kRed+2, kRed+2);
        fHistoRatioMCKDivDataFit->Draw("same,pe");
        DrawGammaSetMarker(fHistoRatioMCK0sDivDataFit, 24, 1.5, kGreen+2, kGreen+2);
        fHistoRatioMCK0sDivDataFit->Draw("same,pe");

        DrawGammaLines(0., 20,1, 1,0.1, kGray+2, 7);
        
        TLegend* legendRatioK = GetAndSetLegend2(0.12, 0.70, 0.32, 0.95, 32,1); 
        legendRatioK->AddEntry(fHistoRatioDataKDivDataFit,"Data: #frac{K^{+}+K^{-}}{2}","p");
        legendRatioK->AddEntry(fHistoRatioMCKDivDataFit,"MC: #frac{K^{+}+K^{-}}{2}","p");
        legendRatioK->AddEntry(fHistoRatioMCK0sDivDataFit,"MC: K^{0}_{s}","p");
        legendRatioK->Draw();

        labelEnergyRatio->Draw();
        labelGeneratorRatio->Draw();
        

        
        canvasRatio->SaveAs(Form("%s/K_RatioComparisonMCAndData_%s_%s.%s",outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), suffix.Data()));
        
      }    
    }
    
    //******************************************************************************************
    //************************ Saving histograms for further processing ************************
    //******************************************************************************************                    
    TString nameOutput = Form("%s/MCInputCompilation%s_%s.root", outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data());
    TFile* fOutput2 = new TFile(nameOutput.Data(),"RECREATE");
    cout << "======================================================" << endl;
    cout << nameOutput << endl;
    cout << "======================================================" << endl;

        if (fEventQuality)                          fEventQuality->Write("NEvents");
        
        if (fHistoMCPi0Pt)                          fHistoMCPi0Pt->Write();
        if (fHistoMCEtaPt)                          fHistoMCEtaPt->Write();
        if (fHistoMCPiPt)                           fHistoMCPiPt->Write();
        if (fHistoMCPiNegPt)                        fHistoMCPiNegPt->Write();
        if (fHistoMCPiPosPt)                        fHistoMCPiPosPt->Write();
        if (fHistoMCKPt)                            fHistoMCKPt->Write();
        if (fHistoMCKNegPt)                         fHistoMCKNegPt->Write();
        if (fHistoMCKPosPt)                         fHistoMCKPosPt->Write();
        if (fHistoMCK0sPt)                          fHistoMCK0sPt->Write();
        if (fHistoRatioMCK0sDivDataFit)             fHistoRatioMCK0sDivDataFit->Write("K0sRatioToDataFit");
        if (fHistoRatioMCKDivDataFit)               fHistoRatioMCKDivDataFit->Write("KRatioToDataFit");
        if (fHistoRatioMCPiDivDataFit)              fHistoRatioMCPiDivDataFit->Write("PiRatioToDataFit");
        
    fOutput2->Write();
    fOutput2->Close();

    return;
}

