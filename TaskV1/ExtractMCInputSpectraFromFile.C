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
                                    Int_t mode                      = 10
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
    
    fOutput2->Write();
    fOutput2->Close();

    return;
}

