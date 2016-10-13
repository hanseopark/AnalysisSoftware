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
#include "TSpline.h"
// #include "TSpline3.h"
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
                                    TString nameExternalInput       = "ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_2016_07_04.root",
                                    Double_t minPt                  = 0,
                                    Double_t maxPt                  = 20
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
//     if (!(fMode == 10 || fMode == 11) ){
//         cout << "ERROR: You are running the wrong macro, this macro is only designed for merged cluster analysis." << endl;
//         return;
//     }
    
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
    } else if (mode == 2 ){
        nameMainDir                 = "GammaConvCalo";
    } else if (mode == 4 ){
        nameMainDir                 = "GammaCalo";
    } else if (mode == 0 ){
        nameMainDir                 = "GammaConv";
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
        
    Double_t ptBinning[65]              = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                           1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                           2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
                                           3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,
                                           5.0, 5.4, 5.8, 6.2, 6.6, 7.0, 7.5, 8.0, 8.5, 9.0,
                                           10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 25.0,
                                           30.0, 35.0, 40.0, 45.0, 50.0};
    Int_t nBinsX                        = 64;    
        
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
    TH1D* fHistoMCPi0PtRebinned         = (TH1D*)fHistoMCPi0Pt->Clone("fHistoMCPi0PtRebinned");
    fHistoMCPi0PtRebinned->Rebin(nBinsX,"fHistoMCPi0PtRebinned2",ptBinning);
    fHistoMCPi0PtRebinned               = (TH1D*)gDirectory->Get("fHistoMCPi0PtRebinned2");
    
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
    TH1D* fHistoMCEtaPtRebinned         = (TH1D*)fHistoMCEtaPt->Clone("fHistoMCEtaPtRebinned");
    fHistoMCEtaPtRebinned->Rebin(nBinsX,"fHistoMCEtaPtRebinned2",ptBinning);
    fHistoMCEtaPtRebinned               = (TH1D*)gDirectory->Get("fHistoMCEtaPtRebinned2");

    // ------------ read 2D primary particle histo vs pt --------
    TH2D* fHistoPrimPartYSource         = NULL;
    TH1D* fHistoMCPiNegY                = NULL;
    TH1D* fHistoMCPiPosY                = NULL;
    TH1D* fHistoMCKNegY                 = NULL;
    TH1D* fHistoMCKPosY                 = NULL;
    TH1D* fHistoMCK0sY                  = NULL;
    TH1D* fHistoMCK0lY                  = NULL;
    TH1D* fHistoMCLambdaY               = NULL;
    fHistoPrimPartYSource               = (TH2D*)MCContainer->FindObject("MC_Primary_Y_Source");  
    if (fHistoPrimPartYSource){
        fHistoMCPiPosY                  = (TH1D*)fHistoPrimPartYSource->ProjectionX( "MC_PosPi_Y", fHistoPrimPartYSource->GetYaxis()->FindBin(0.), 
                                                                                    fHistoPrimPartYSource->GetYaxis()->FindBin(0.),"e");     
        fHistoMCPiNegY                  = (TH1D*)fHistoPrimPartYSource->ProjectionX( "MC_NegPi_Y", fHistoPrimPartYSource->GetYaxis()->FindBin(1.), 
                                                                                    fHistoPrimPartYSource->GetYaxis()->FindBin(1.),"e");     
        fHistoMCKPosY                   = (TH1D*)fHistoPrimPartYSource->ProjectionX( "MC_PosK_Y", fHistoPrimPartYSource->GetYaxis()->FindBin(2.), 
                                                                                    fHistoPrimPartYSource->GetYaxis()->FindBin(2.),"e");     
        fHistoMCKNegY                   = (TH1D*)fHistoPrimPartYSource->ProjectionX( "MC_NegK_Y", fHistoPrimPartYSource->GetYaxis()->FindBin(3.), 
                                                                                    fHistoPrimPartYSource->GetYaxis()->FindBin(3.),"e");     
        fHistoMCK0sY                    = (TH1D*)fHistoPrimPartYSource->ProjectionX( "MC_K0s_Y", fHistoPrimPartYSource->GetYaxis()->FindBin(4.), 
                                                                                    fHistoPrimPartYSource->GetYaxis()->FindBin(4.),"e");     
        fHistoMCK0lY                    = (TH1D*)fHistoPrimPartYSource->ProjectionX( "MC_K0l_Y", fHistoPrimPartYSource->GetYaxis()->FindBin(5.), 
                                                                                    fHistoPrimPartYSource->GetYaxis()->FindBin(5.),"e");     
        fHistoMCLambdaY                 = (TH1D*)fHistoPrimPartYSource->ProjectionX( "MC_Lambda_Y", fHistoPrimPartYSource->GetYaxis()->FindBin(6.), 
                                                                                    fHistoPrimPartYSource->GetYaxis()->FindBin(6.),"e");     
    } 
    
    // ------------ read 2D primary particle histo vs pt --------
    TH2D* fHistoPrimPartPtSource        = NULL;
    TH1D* fHistoMCPiNegPt               = NULL;
    TH1D* fHistoMCPiPosPt               = NULL;
    TH1D* fHistoMCKNegPt                = NULL;
    TH1D* fHistoMCKPosPt                = NULL;
    TH1D* fHistoMCK0sPt                 = NULL;
    TH1D* fHistoMCK0lPt                 = NULL;
    TH1D* fHistoMCLambdaPt              = NULL;
    fHistoPrimPartPtSource              = (TH2D*)MCContainer->FindObject("MC_Primary_Pt_Source");  
    if (fHistoPrimPartPtSource){
        fHistoMCPiPosPt               = (TH1D*)fHistoPrimPartPtSource->ProjectionX( "MC_PosPi_Pt", fHistoPrimPartPtSource->GetYaxis()->FindBin(0.), 
                                                                                    fHistoPrimPartPtSource->GetYaxis()->FindBin(0.),"e");     
        fHistoMCPiNegPt               = (TH1D*)fHistoPrimPartPtSource->ProjectionX( "MC_NegPi_Pt", fHistoPrimPartPtSource->GetYaxis()->FindBin(1.), 
                                                                                    fHistoPrimPartPtSource->GetYaxis()->FindBin(1.),"e");     
        fHistoMCKPosPt                = (TH1D*)fHistoPrimPartPtSource->ProjectionX( "MC_PosK_Pt", fHistoPrimPartPtSource->GetYaxis()->FindBin(2.), 
                                                                                    fHistoPrimPartPtSource->GetYaxis()->FindBin(2.),"e");     
        fHistoMCKNegPt                = (TH1D*)fHistoPrimPartPtSource->ProjectionX( "MC_NegK_Pt", fHistoPrimPartPtSource->GetYaxis()->FindBin(3.), 
                                                                                    fHistoPrimPartPtSource->GetYaxis()->FindBin(3.),"e");     
        fHistoMCK0sPt                 = (TH1D*)fHistoPrimPartPtSource->ProjectionX( "MC_K0s_Pt", fHistoPrimPartPtSource->GetYaxis()->FindBin(4.), 
                                                                                    fHistoPrimPartPtSource->GetYaxis()->FindBin(4.),"e");     
        fHistoMCK0lPt                 = (TH1D*)fHistoPrimPartPtSource->ProjectionX( "MC_K0l_Pt", fHistoPrimPartPtSource->GetYaxis()->FindBin(5.), 
                                                                                    fHistoPrimPartPtSource->GetYaxis()->FindBin(5.),"e");     
        fHistoMCLambdaPt              = (TH1D*)fHistoPrimPartPtSource->ProjectionX( "MC_Lambda_Pt", fHistoPrimPartPtSource->GetYaxis()->FindBin(6.), 
                                                                                    fHistoPrimPartPtSource->GetYaxis()->FindBin(6.),"e");     
    } else {
        fHistoMCPiNegPt               = (TH1D*)MCContainer->FindObject("MC_NegPi_Pt");  
        fHistoMCPiPosPt               = (TH1D*)MCContainer->FindObject("MC_PosPi_Pt");
        fHistoMCKNegPt                = (TH1D*)MCContainer->FindObject("MC_NegK_Pt");  
        fHistoMCKPosPt                = (TH1D*)MCContainer->FindObject("MC_PosK_Pt");  
        fHistoMCK0sPt                 = (TH1D*)MCContainer->FindObject("MC_K0s_Pt");  
        fHistoMCK0lPt                 = (TH1D*)MCContainer->FindObject("MC_K0l_Pt");          
    }    
    
    // ----------- charged pions ------------------------------
    cout << "reading charged pions" << endl;
    TH1D* fHistoMCPiPt                  = (TH1D*)fHistoMCPiNegPt->Clone("MC_PiCh_All_Pt");
    fHistoMCPiPt->Sumw2();
    fHistoMCPiPt->Add(fHistoMCPiPosPt);
    fHistoMCPiPt->Scale(0.5);
    TH1D* fHistoMCPiPtRebinned          = (TH1D*)fHistoMCPiPt->Clone("fHistoMCPiPtRebinned");
    fHistoMCPiPtRebinned->Rebin(nBinsX,"fHistoMCPiPtRebinned2",ptBinning);
    fHistoMCPiPtRebinned                = (TH1D*)gDirectory->Get("fHistoMCPiPtRebinned2");
    // ----------- charged kaons ------------------------------
    cout << "reading charged kaons" << endl;
    TH1D* fHistoMCKPt                   = (TH1D*)fHistoMCKNegPt->Clone("MC_KCh_All_Pt");
    fHistoMCKPt->Sumw2();
    fHistoMCKPt->Add(fHistoMCKPosPt);
    fHistoMCKPt->Scale(0.5);
    TH1D* fHistoMCKPtRebinned           = (TH1D*)fHistoMCKPt->Clone("fHistoMCKPtRebinned");
    fHistoMCKPtRebinned->Rebin(nBinsX,"fHistoMCKPtRebinned2",ptBinning);
    fHistoMCKPtRebinned                 = (TH1D*)gDirectory->Get("fHistoMCKPtRebinned2");
    // ----------- neutral kaons ------------------------------
    cout << "reading K0s" << endl;
    TH1D* fHistoMCK0sPtYield            = (TH1D*)fHistoMCK0sPt->Clone("MC_K0s_Pt");
    fHistoMCK0sPtYield->Sumw2();
    TH1D* fHistoMCK0sPtRebinned         = (TH1D*)fHistoMCK0sPt->Clone("fHistoMCK0sPtRebinned");
    fHistoMCK0sPtRebinned->Rebin(nBinsX,"fHistoMCK0sPtRebinned2",ptBinning);
    fHistoMCK0sPtRebinned               = (TH1D*)gDirectory->Get("fHistoMCK0sPtRebinned2");

    // ----------- neutral kaons ------------------------------
    cout << "reading K0l" << endl;
    TH1D* fHistoMCK0lPtYield            = (TH1D*)fHistoMCK0lPt->Clone("MC_K0l_Pt_Yield");
    fHistoMCK0lPtYield->Sumw2();
    TH1D* fHistoMCK0lPtRebinned         = (TH1D*)fHistoMCK0lPt->Clone("fHistoMCK0lPtRebinned");
    fHistoMCK0lPtRebinned->Rebin(nBinsX,"fHistoMCK0lPtRebinned2",ptBinning);
    fHistoMCK0lPtRebinned               = (TH1D*)gDirectory->Get("fHistoMCK0lPtRebinned2");
    
    // ----------- Lambda -------------------------------------
    TH1D* fHistoMCLambdaPtRebinned      = NULL;
    TH1D* fHistoMCLambdaPtYield         = NULL;
    if (fHistoMCLambdaPt){
        fHistoMCLambdaPtYield           = (TH1D*)fHistoMCLambdaPt->Clone("MC_Lambda_Pt_Yield");
        fHistoMCLambdaPtYield->Sumw2();
        fHistoMCLambdaPtRebinned        = (TH1D*)fHistoMCLambdaPt->Clone("fHistoMCLambdaPtRebinned");
        fHistoMCLambdaPtRebinned->Rebin(nBinsX,"fHistoMCLambdaPtRebinned2",ptBinning);
        fHistoMCLambdaPtRebinned        = (TH1D*)gDirectory->Get("fHistoMCLambdaPtRebinned2");
    }
    
    // ----------- secondary pions from K0s, K0l, Lambda ------
    //---> K0s
    TH2D* fHistoMCSecPi0PtSource                = (TH2D*)MCContainer->FindObject("MC_SecPi0_Pt_Source");
    TH1D* fHistoMCSecPi0FromK0sPt               = (TH1D*)fHistoMCSecPi0PtSource->ProjectionX(   "MCSecPi0FromK0s", fHistoMCSecPi0PtSource->GetYaxis()->FindBin(1.), 
                                                                                                fHistoMCSecPi0PtSource->GetYaxis()->FindBin(1.),"e");     
    TH1D* fHistoMCSecPi0FromK0sPtRebinned       = (TH1D*)fHistoMCSecPi0FromK0sPt->Clone("MCSecPi0FromK0sRebinned");
    fHistoMCSecPi0FromK0sPtRebinned->Rebin(nBinsX,"MCSecPi0FromK0sRebinned2",ptBinning);
    fHistoMCSecPi0FromK0sPtRebinned             = (TH1D*)gDirectory->Get("MCSecPi0FromK0sRebinned2");
    TH1D* fHistoRatioPi0FromK0sDivK0s           = (TH1D*)fHistoMCSecPi0FromK0sPtRebinned->Clone("ratioPi0FromK0s");
    fHistoRatioPi0FromK0sDivK0s->Divide(fHistoMCSecPi0FromK0sPtRebinned, fHistoMCK0sPtRebinned);

    //---> K0l
    TH1D* fHistoMCSecPi0FromK0lPt               = (TH1D*)fHistoMCSecPi0PtSource->ProjectionX(   "MCSecPi0FromK0l", fHistoMCSecPi0PtSource->GetYaxis()->FindBin(3.), 
                                                                                                fHistoMCSecPi0PtSource->GetYaxis()->FindBin(3.),"e");     
    TH1D* fHistoMCSecPi0FromK0lPtRebinned       = (TH1D*)fHistoMCSecPi0FromK0lPt->Clone("MCSecPi0FromK0lRebinned");
    fHistoMCSecPi0FromK0lPtRebinned->Rebin(nBinsX,"MCSecPi0FromK0lRebinned2",ptBinning);
    fHistoMCSecPi0FromK0lPtRebinned             = (TH1D*)gDirectory->Get("MCSecPi0FromK0lRebinned2");
    TH1D* fHistoRatioPi0FromK0lDivK0l           = (TH1D*)fHistoMCSecPi0FromK0lPtRebinned->Clone("ratioPi0FromK0l");
    fHistoRatioPi0FromK0lDivK0l->Divide(fHistoMCSecPi0FromK0lPtRebinned, fHistoMCK0lPtRebinned);

    //---> Lambda
    TH1D* fHistoMCSecPi0FromLambdaPt            = (TH1D*)fHistoMCSecPi0PtSource->ProjectionX(   "MCSecPi0FromLambda", fHistoMCSecPi0PtSource->GetYaxis()->FindBin(2.), 
                                                                                                fHistoMCSecPi0PtSource->GetYaxis()->FindBin(2.),"e");     
    TH1D* fHistoMCSecPi0FromLambdaPtRebinned    = (TH1D*)fHistoMCSecPi0FromLambdaPt->Clone("MCSecPi0FromLambdaRebinned");
    fHistoMCSecPi0FromLambdaPtRebinned->Rebin(nBinsX,"MCSecPi0FromLambdaRebinned2",ptBinning);
    TH1D* fHistoRatioPi0FromLambdaDivLambda     = NULL;
    if (fHistoMCLambdaPt){
        fHistoRatioPi0FromLambdaDivLambda       = (TH1D*)fHistoMCSecPi0FromLambdaPtRebinned->Clone("ratioPi0FromLambda");
        fHistoRatioPi0FromLambdaDivLambda->Divide(fHistoMCSecPi0FromLambdaPtRebinned, fHistoMCLambdaPtRebinned);
    }    
    
    // scale yield to per event quantities
    fHistoMCK0sPtYield->Scale(1./nEvtMC);
    fHistoMCK0lPtYield->Scale(1./nEvtMC);
    fHistoMCSecPi0FromK0sPt->Scale(1./nEvtMC);
    fHistoMCSecPi0FromK0lPt->Scale(1./nEvtMC);
    if (fHistoMCLambdaPt){
        fHistoMCLambdaPtYield->Scale(1./nEvtMC);
        fHistoMCSecPi0FromLambdaPt->Scale(1./nEvtMC);
    }
    
    // scale yield to fully invariant numbers
    ScaleMCYield(fHistoMCPi0Pt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCPi0PtRebinned,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCEtaPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCEtaPtRebinned,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCPiPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCPiPtRebinned,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCPiNegPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCPiPosPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCKPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCKPtRebinned,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCKNegPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCKPosPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCK0sPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCK0sPtRebinned,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCK0lPt,  deltaRapid,  scaling,  nEvtMC );
    ScaleMCYield(fHistoMCK0lPtRebinned,  deltaRapid,  scaling,  nEvtMC );
    if (fHistoMCLambdaPt){
        ScaleMCYield(fHistoMCLambdaPt,  deltaRapid,  scaling,  nEvtMC );
        ScaleMCYield(fHistoMCLambdaPtRebinned,  deltaRapid,  scaling,  nEvtMC );
    }
    
    TH1D* fHistoRatioMCPi0DivPi         = NULL;
    TH1D* fHistoRatioMCK0sDivK          = NULL;
    TH1D* fHistoRatioMCK0lDivK          = NULL;
    TH1D* fHistoRatioMCKDivPi           = NULL;
    TH1D* fHistoRatioMCK0sDivPi0        = NULL;
    TH1D* fHistoRatioMCEtaDivPi0        = NULL;
    TCanvas *canvasRatio = new TCanvas("canvasRatio","canvasRatio",1000,800);
    DrawGammaCanvasSettings( canvasRatio, 0.11, 0.02, 0.02, 0.08);
    canvasRatio->cd();
    canvasRatio->SetLogy(0);
    TLatex *labelEnergyRatio        = new TLatex(0.15,0.16,Form("%s",fCollisionSystem.Data()));
    SetStyleTLatex( labelEnergyRatio, 0.04,4);
    TLatex *labelGeneratorRatio     = new TLatex(0.15,0.12,Form("%s",fGeneratorName.Data()));
    SetStyleTLatex( labelGeneratorRatio, 0.04,4);
    
    // build ratio eta/pi^0
    if (fHistoMCPi0Pt && fHistoMCEtaPt){
        fHistoRatioMCEtaDivPi0 = (TH1D*)fHistoMCEtaPtRebinned->Clone("fHistoRatioMCPi0DivPi");
        fHistoRatioMCEtaDivPi0->Divide(fHistoRatioMCEtaDivPi0,fHistoMCPi0PtRebinned);

        Double_t maxPtForFit    = 5;
        Double_t minPtForFit    = 15;
        if (fGeneratorName.Contains("p_{")){
            minPtForFit         = 8;
            maxPtForFit         = 20;
        }    
        TF1* etaToPi0ConstMC    = new TF1("etaToPi0ConstMC","[0]",minPtForFit,maxPtForFit);
        fHistoRatioMCEtaDivPi0->Fit(etaToPi0ConstMC,"QRME0","",minPtForFit,maxPtForFit);
        cout << "***********************************************************************************************************" << endl;
        cout << "***********************************************************************************************************" << endl;
        cout << "***********************************************************************************************************" << endl;
        cout << "high pt eta/pi0 - MC: " << etaToPi0ConstMC->GetParameter(0) << "+-"<< etaToPi0ConstMC->GetParError(0) << endl;
        cout << "***********************************************************************************************************" << endl;
        cout << "***********************************************************************************************************" << endl;
        cout << "***********************************************************************************************************" << endl;

        
        DrawAutoGammaMesonHistos(   fHistoRatioMCEtaDivPi0, 
                            "", "#it{p}_{T} (GeV/#it{c})", "#eta/ #pi^{0}", 
                            kFALSE, 10, 1e-10, kFALSE,
                            kTRUE, 0., 1.1, 
                            kTRUE, minPt, maxPt);
        fHistoRatioMCEtaDivPi0->GetYaxis()->SetTitleOffset(1.2);
        DrawGammaSetMarker(fHistoRatioMCEtaDivPi0, 20, 1.5, kAzure-6, kAzure-6);
        fHistoRatioMCEtaDivPi0->DrawClone("pe");
        
        etaToPi0ConstMC->SetLineColor(kAzure-6);
        etaToPi0ConstMC->SetLineStyle(7);
        etaToPi0ConstMC->Draw("same");
        
//         DrawGammaLines(0., 20,1, 1,0.1, kGray+2, 7);
        labelEnergyRatio->Draw();
        labelGeneratorRatio->Draw();
        
        canvasRatio->SaveAs(Form("%s/EtaToPi0_MC_%s_%s.%s",outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), suffix.Data()));    
    }

    
    // build ratio pi0/pi^\m
    if (fHistoMCPi0Pt && fHistoMCPiPt){
        fHistoRatioMCPi0DivPi = (TH1D*)fHistoMCPi0PtRebinned->Clone("fHistoRatioMCPi0DivPi");
        fHistoRatioMCPi0DivPi->Divide(fHistoRatioMCPi0DivPi,fHistoMCPiPtRebinned);

        DrawAutoGammaMesonHistos(   fHistoRatioMCPi0DivPi, 
                            "", "#it{p}_{T} (GeV/#it{c})", "#pi^{0} / #frac{#pi^{+}+#pi^{-}}{2}", 
                            kFALSE, 10, 1e-10, kFALSE,
                            kTRUE, 0.9, 1.2, 
                            kTRUE, minPt, maxPt);
        fHistoRatioMCPi0DivPi->GetYaxis()->SetTitleOffset(1.2);
        DrawGammaSetMarker(fHistoRatioMCPi0DivPi, 20, 1.5, kAzure-6, kAzure-6);
        fHistoRatioMCPi0DivPi->DrawClone("pe");
        
        DrawGammaLines(0., 20,1, 1,0.1, kGray+2, 7);
        labelEnergyRatio->Draw();
        labelGeneratorRatio->Draw();
        
        canvasRatio->SaveAs(Form("%s/Pi0ToPi_MC_%s_%s.%s",outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), suffix.Data()));    
    }
    // build ratio K0s/K^\pm
    if (fHistoMCK0sPt && fHistoMCKPt){
        fHistoRatioMCK0sDivK = (TH1D*)fHistoMCK0sPtRebinned->Clone("fHistoRatioMCK0sDivK");
        fHistoRatioMCK0sDivK->Divide(fHistoRatioMCK0sDivK,fHistoMCKPtRebinned);
     
        if (fHistoMCK0lPt){
            fHistoRatioMCK0lDivK = (TH1D*)fHistoMCK0lPtRebinned->Clone("fHistoRatioMCK0lDivK");
            fHistoRatioMCK0lDivK->Divide(fHistoRatioMCK0lDivK,fHistoMCKPtRebinned);
        }

        TLegend* legendRatioK0ToK = GetAndSetLegend2(0.15, 0.8, 0.42, 0.95, 32,1); 
        DrawAutoGammaMesonHistos(   fHistoRatioMCK0sDivK, 
                            "", "#it{p}_{T} (GeV/#it{c})", "K^{0} / #frac{K^{+}+K^{-}}{2}", 
                            kFALSE, 10, 1e-10, kFALSE,
                            kTRUE, 0.8, 1.2, 
                            kTRUE, minPt, maxPt);
        fHistoRatioMCK0sDivK->GetYaxis()->SetTitleOffset(1.2);
        DrawGammaSetMarker(fHistoRatioMCK0sDivK, 20, 1.5, kAzure-6, kAzure-6);
        fHistoRatioMCK0sDivK->DrawClone("pe");
        legendRatioK0ToK->AddEntry(fHistoRatioMCK0sDivK,"K^{0}_{S}","p");
        
        if (fHistoRatioMCK0lDivK){
            DrawGammaSetMarker(fHistoRatioMCK0lDivK, 24, 1.6, kGreen+2, kGreen+2);
            fHistoRatioMCK0lDivK->DrawClone("same,pe");
            legendRatioK0ToK->AddEntry(fHistoRatioMCK0lDivK,"K^{0}_{L}","p");
        }
        legendRatioK0ToK->Draw();

        DrawGammaLines(0., 20,1, 1,0.1, kGray+2, 7);
        labelEnergyRatio->Draw();
        labelGeneratorRatio->Draw();
        
        canvasRatio->SaveAs(Form("%s/K0ToK_MC_%s_%s.%s",outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), suffix.Data()));
    }    
    // build ratio K^\pm/pi^\pm
    if ( fHistoMCKPt && fHistoMCPiPt){    
        fHistoRatioMCKDivPi = (TH1D*)fHistoMCKPtRebinned->Clone("fHistoRatioMCKDivPi");
        fHistoRatioMCKDivPi->Divide(fHistoRatioMCKDivPi,fHistoMCPiPtRebinned);

        DrawAutoGammaMesonHistos(   fHistoRatioMCKDivPi, 
                            "", "#it{p}_{T} (GeV/#it{c})", "#frac{K^{+}+K^{-}}{2} / #frac{#pi^{+}+#pi^{-}}{2}", 
                            kFALSE, 10, 1e-10, kFALSE,
                            kFALSE, 0.9, 1.2, 
                            kTRUE, minPt, maxPt);
        fHistoRatioMCKDivPi->GetYaxis()->SetTitleOffset(1.2);
        DrawGammaSetMarker(fHistoRatioMCKDivPi, 20, 1.5, kAzure-6, kAzure-6);
        fHistoRatioMCKDivPi->DrawClone("pe");
        
//         DrawGammaLines(0., 20,1, 1,0.1, kGray+2, 7);
        labelEnergyRatio->Draw();
        labelGeneratorRatio->Draw();
        
    
        canvasRatio->SaveAs(Form("%s/KToPi_MC_%s_%s.%s",outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), suffix.Data()));        
    }   

    // build ratio K0s/pi0
    if ( fHistoMCK0sPt && fHistoMCPi0PtRebinned){    
        fHistoRatioMCK0sDivPi0 = (TH1D*)fHistoMCK0sPtRebinned->Clone("fHistoRatioMCKDivPi");
        fHistoRatioMCK0sDivPi0->Divide(fHistoRatioMCK0sDivPi0,fHistoMCPi0PtRebinned);

        DrawAutoGammaMesonHistos(   fHistoRatioMCK0sDivPi0, 
                            "", "#it{p}_{T} (GeV/#it{c})", "K^{0}_{s} / #pi^{0}", 
                            kFALSE, 10, 1e-10, kFALSE,
                            kFALSE, 0.9, 1.2, 
                            kTRUE, minPt, maxPt);
        fHistoRatioMCK0sDivPi0->GetYaxis()->SetTitleOffset(1.2);
        DrawGammaSetMarker(fHistoRatioMCK0sDivPi0, 20, 1.5, kAzure-6, kAzure-6);
        fHistoRatioMCK0sDivPi0->DrawClone("pe");
        
//         DrawGammaLines(0., 20,1, 1,0.1, kGray+2, 7);
        labelEnergyRatio->Draw();
        labelGeneratorRatio->Draw();
        
        canvasRatio->SaveAs(Form("%s/K0sToPi0_MC_%s_%s.%s",outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), suffix.Data()));        
    }   
    
    // build ratio pi0 from K0s/ K0s
    if ( fHistoRatioPi0FromK0sDivK0s){    

        DrawAutoGammaMesonHistos(   fHistoRatioPi0FromK0sDivK0s, 
                            "", "#it{p}_{T} (GeV/#it{c})", "#pi^{0} from K^{0}_{s} / K^{0}_{s}", 
                            kFALSE, 10, 1e-10, kFALSE,
                            kTRUE, 0, 0.2, 
                            kTRUE, minPt, maxPt);
        fHistoRatioPi0FromK0sDivK0s->GetYaxis()->SetTitleOffset(1.2);
        DrawGammaSetMarker(fHistoRatioPi0FromK0sDivK0s, 20, 1.5, kAzure-6, kAzure-6);
        fHistoRatioPi0FromK0sDivK0s->DrawClone("pe");
                
        
//         DrawGammaLines(0., 20,1, 1,0.1, kGray+2, 7);
        labelEnergyRatio->Draw();
        labelGeneratorRatio->Draw();
        
        canvasRatio->SaveAs(Form("%s/Pi0FromK0sToK0s_MC_%s_%s.%s",outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), suffix.Data()));        
    }   
    
    
    
    // build ratio pi0 from K0l/ K0l
    if ( fHistoRatioPi0FromK0lDivK0l){    
        canvasRatio->cd();
        canvasRatio->SetTopMargin(0.035);
        DrawAutoGammaMesonHistos(   fHistoRatioPi0FromK0lDivK0l, 
                            "", "#it{p}_{T} (GeV/#it{c})", "#pi^{0} from K^{0}_{l} / K^{0}_{l}", 
                            kFALSE, 10, 1e-10, kFALSE,
                            kTRUE, 0, 0.001, 
                            kTRUE, minPt, maxPt);
        fHistoRatioPi0FromK0lDivK0l->GetYaxis()->SetTitleOffset(1.2);
        DrawGammaSetMarker(fHistoRatioPi0FromK0lDivK0l, 20, 1.5, kAzure-6, kAzure-6);
        fHistoRatioPi0FromK0lDivK0l->DrawClone("pe");
                
        
//         DrawGammaLines(0., 20,1, 1,0.1, kGray+2, 7);
        labelEnergyRatio->Draw();
        labelGeneratorRatio->Draw();
        
        canvasRatio->SaveAs(Form("%s/Pi0FromK0lToK0l_MC_%s_%s.%s",outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), suffix.Data()));        
        canvasRatio->SetTopMargin(0.02);
    }   

    // build ratio pi0 from Lambda/ Lambda
    if ( fHistoRatioPi0FromLambdaDivLambda){    
        canvasRatio->cd();
        canvasRatio->SetLogy(1);
        canvasRatio->SetTopMargin(0.02);
        DrawAutoGammaMesonHistos(   fHistoRatioPi0FromLambdaDivLambda, 
                            "", "#it{p}_{T} (GeV/#it{c})", "#pi^{0} from #Lambda / #Lambda", 
                            kFALSE, 10, 1e-10, kFALSE,
                            kTRUE, 0.0001, 100, 
                            kTRUE, minPt, maxPt);
        fHistoRatioPi0FromLambdaDivLambda->GetYaxis()->SetTitleOffset(1.2);
        DrawGammaSetMarker(fHistoRatioPi0FromLambdaDivLambda, 20, 1.5, kAzure-6, kAzure-6);
        fHistoRatioPi0FromLambdaDivLambda->DrawClone("pe");
                
        
//         DrawGammaLines(0., 20,1, 1,0.1, kGray+2, 7);
        labelEnergyRatio->Draw();
        labelGeneratorRatio->Draw();
        
        canvasRatio->SaveAs(Form("%s/Pi0FromLambdaToLambda_MC_%s_%s.%s",outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), suffix.Data()));        
        canvasRatio->SetTopMargin(0.02);
        canvasRatio->SetLogy(0);
    }   
    
    
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
      if (optionEnergy.CompareTo("2.76TeV") == 0 || optionEnergy.CompareTo("7TeV") == 0){
        fHistoChargedPionData       = (TH1D*)fileDataInput->Get("histoChargedPionSpecPubStat2760GeV");
        if( optionEnergy.CompareTo("7TeV") == 0 ) fHistoChargedPionData       = (TH1D*)fileDataInput->Get("histoChargedPionSpecPubStat7TeV");
        TF1* fitChargedPions        = FitObject("l","fitChargedPions","Pi0",fHistoChargedPionData,0.1,20.,NULL,"QNRMEI");
        TSpline5* paramPions        = new TSpline5(fHistoChargedPionData);
        fHistoRatioDataPiDivDataFit = CalculateHistoRatioToSpline (fHistoChargedPionData, paramPions); 
        fHistoRatioMCPiDivDataFit   = CalculateHistoRatioToSpline (fHistoMCPiPtRebinned, paramPions); 
        fHistoRatioMCPi0DivDataFit  = CalculateHistoRatioToSpline (fHistoMCPi0PtRebinned, paramPions); 
        
        fHistoChargedKaonData       = (TH1D*)fileDataInput->Get("histoChargedKaonSpecPubStat2760GeV");
        if( optionEnergy.CompareTo("7TeV") == 0 ) fHistoChargedKaonData       = (TH1D*)fileDataInput->Get("histoChargedKaonSpecPubStat7TeV");
        TF1* fitChargedKaons        =  FitObject("l","ptDistribution","K",fHistoChargedKaonData,0.1,20.,NULL,"QNRMEI");
        TSpline5* paramKaons        = new TSpline5(fHistoChargedKaonData);
        
        fHistoRatioDataKDivDataFit  = CalculateHistoRatioToSpline (fHistoChargedKaonData, paramKaons); 
        fHistoRatioMCK0sDivDataFit  = CalculateHistoRatioToSpline (fHistoMCK0sPtRebinned, paramKaons); 
        fHistoRatioMCKDivDataFit    = CalculateHistoRatioToSpline (fHistoMCKPtRebinned, paramKaons); 
        
        
        
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

        DrawGammaSetMarker(fHistoMCPiPtRebinned, 20, 1.5, kRed+2, kRed+2);
        fHistoMCPiPtRebinned->Draw("same,pe");
        DrawGammaSetMarker(fHistoMCPi0PtRebinned, 24, 1.5, kGreen+2, kGreen+2);
        fHistoMCPi0PtRebinned->Draw("same,pe");

        fitChargedPions->SetLineColor(kAzure-6);
        fitChargedPions->Draw("same");
        paramPions->SetLineColor(kAzure+2);
        paramPions->Draw("same");
        
        fHistoChargedPionData->Draw("same,pe");

        TLegend* legendSpectraPi = GetAndSetLegend2(0.73, 0.70, 0.95, 0.95, 32,1); 
        legendSpectraPi->AddEntry(fHistoChargedPionData,"Data: #frac{#pi^{+}+#pi^{-}}{2}","p");
        legendSpectraPi->AddEntry(fitChargedPions,"Data: fit","l");
        legendSpectraPi->AddEntry(fHistoMCPiPtRebinned,"MC: #frac{#pi^{+}+#pi^{-}}{2}","p");
        legendSpectraPi->AddEntry(fHistoMCPi0PtRebinned,"MC: #pi^{0}","p");
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

        DrawGammaSetMarker(fHistoMCKPtRebinned, 20, 1.5, kRed+2, kRed+2);
        fHistoMCKPtRebinned->Draw("same,pe");
        DrawGammaSetMarker(fHistoMCK0sPtRebinned, 24, 1.5, kGreen+2, kGreen+2);
        fHistoMCK0sPtRebinned->Draw("same,pe");
        fitChargedKaons->SetLineColor(kAzure-6);
        fitChargedKaons->Draw("same");
        fHistoChargedKaonData->Draw("same,pe");
        paramKaons->SetLineColor(kAzure+2);
        paramKaons->Draw("same");

        TLegend* legendSpectraK = GetAndSetLegend2(0.73, 0.70, 0.95, 0.95, 32,1); 
        legendSpectraK->AddEntry(fHistoChargedKaonData,"Data: #frac{K^{+}+K^{-}}{2}","p");
        legendSpectraK->AddEntry(fitChargedKaons,"Data: fit","l");
        legendSpectraK->AddEntry(fHistoMCKPtRebinned,"MC: #frac{K^{+}+K^{-}}{2}","p");
        legendSpectraK->AddEntry(fHistoMCK0sPtRebinned,"MC: K^{0}_{s}","p");
        legendSpectraK->Draw();
        labelEnergySpectra->Draw();
        labelGeneratorSpectra->Draw();
        
        canvasFitQA->SaveAs(Form("%s/K_ComparisonMCAndData_%s_%s.%s",outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), suffix.Data()));
        
        canvasRatio->cd();
        Double_t maxScaleY      = 2;
        Double_t scaleFac       = 0.1;
        if (optionPeriod.CompareTo("LHC15a3a") == 0 || optionPeriod.CompareTo("LHC15g1a") == 0){
            maxScaleY           = 2;            
            fHistoRatioMCPiDivDataFit->Scale(scaleFac);
            fHistoRatioMCPi0DivDataFit->Scale(scaleFac);
            fHistoRatioMCKDivDataFit->Scale(scaleFac);
            fHistoRatioMCK0sDivDataFit->Scale(scaleFac);
            
        }    
        DrawAutoGammaMesonHistos(   fHistoRatioDataPiDivDataFit, 
                            "", "#it{p}_{T} (GeV/#it{c})", "Data/Fit", 
                            kFALSE, 10, 1e-10, kFALSE,
                            kTRUE, 0., maxScaleY, 
                            kTRUE, minPt, maxPt);
        fHistoRatioDataPiDivDataFit->GetYaxis()->SetTitleOffset(0.85);
        DrawGammaSetMarker(fHistoRatioDataPiDivDataFit, 20, 1.5, kAzure-6, kAzure-6);
        fHistoRatioDataPiDivDataFit->DrawClone("pe");
        DrawGammaSetMarker(fHistoRatioMCPiDivDataFit, 20, 1.5, kRed+2, kRed+2);
        fHistoRatioMCPiDivDataFit->Draw("same,pe");
        DrawGammaSetMarker(fHistoRatioMCPi0DivDataFit, 24, 1.5, kGreen+2, kGreen+2);
        fHistoRatioMCPi0DivDataFit->Draw("same,pe");


        TLegend* legendRatioPi = GetAndSetLegend2(0.12, 0.70, 0.32, 0.95, 32,1); 
        legendRatioPi->AddEntry(fHistoRatioDataPiDivDataFit,"Data: #frac{#pi^{+}+#pi^{-}}{2}","p");
        if (optionPeriod.CompareTo("LHC15a3a") == 0 || optionPeriod.CompareTo("LHC15g1a") == 0){
            legendRatioPi->AddEntry(fHistoRatioMCPiDivDataFit,Form("MC: #frac{#pi^{+}+#pi^{-}}{2}, scaled by %1.2f",scaleFac) ,"p");
            legendRatioPi->AddEntry(fHistoRatioMCPi0DivDataFit,Form("MC: #pi^{0}, scaled by %1.2f",scaleFac),"p");
        } else {
            legendRatioPi->AddEntry(fHistoRatioMCPiDivDataFit,"MC: #frac{#pi^{+}+#pi^{-}}{2}","p");
            legendRatioPi->AddEntry(fHistoRatioMCPi0DivDataFit,"MC: #pi^{0}","p");
        }    
        legendRatioPi->Draw();
        
        DrawGammaLines(minPt, maxPt,1, 1,0.1, kGray+2, 7);
        labelEnergyRatio->Draw();
        labelGeneratorRatio->Draw();
        
        canvasRatio->SaveAs(Form("%s/Pi_RatioComparisonMCAndData_%s_%s.%s",outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), suffix.Data()));

        canvasRatio->cd();
        
        DrawAutoGammaMesonHistos(   fHistoRatioDataKDivDataFit, 
                            "", "#it{p}_{T} (GeV/#it{c})", "Data/Fit", 
                            kFALSE, 10, 1e-10, kFALSE,
                            kTRUE, 0., maxScaleY, 
                            kTRUE, minPt, maxPt);
        fHistoRatioDataKDivDataFit->GetYaxis()->SetTitleOffset(0.85);
        DrawGammaSetMarker(fHistoRatioDataKDivDataFit, 20, 1.5, kAzure-6, kAzure-6);
        fHistoRatioDataKDivDataFit->DrawClone("pe");
        
        
        DrawGammaSetMarker(fHistoRatioMCKDivDataFit, 20, 1.5, kRed+2, kRed+2);
        fHistoRatioMCKDivDataFit->Draw("same,pe");
        DrawGammaSetMarker(fHistoRatioMCK0sDivDataFit, 24, 1.5, kGreen+2, kGreen+2);
        fHistoRatioMCK0sDivDataFit->Draw("same,pe");

        DrawGammaLines(minPt, maxPt,1, 1,0.1, kGray+2, 7);
        
        TLegend* legendRatioK = GetAndSetLegend2(0.12, 0.70, 0.32, 0.95, 32,1); 
        legendRatioK->AddEntry(fHistoRatioDataKDivDataFit,"Data: #frac{K^{+}+K^{-}}{2}","p");
        if (optionPeriod.CompareTo("LHC15a3a") == 0 || optionPeriod.CompareTo("LHC15g1a") == 0){
            legendRatioK->AddEntry(fHistoRatioMCKDivDataFit,Form("MC: #frac{K^{+}+K^{-}}{2}, scaled by %1.2f",scaleFac) ,"p");
            legendRatioK->AddEntry(fHistoRatioMCK0sDivDataFit,Form("MC: K^{0}_{s}, scaled by %1.2f",scaleFac) ,"p");
        } else {
            legendRatioK->AddEntry(fHistoRatioMCKDivDataFit,"MC: #frac{K^{+}+K^{-}}{2}","p");
            legendRatioK->AddEntry(fHistoRatioMCK0sDivDataFit,"MC: K^{0}_{s}","p");            
        }    
        legendRatioK->Draw();

        labelEnergyRatio->Draw();
        labelGeneratorRatio->Draw();
        
        canvasRatio->SaveAs(Form("%s/K_RatioComparisonMCAndData_%s_%s.%s",outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), suffix.Data()));
        
      }    
    }
    delete canvasRatio;        
    
    //******************************************************************************************
    //************************ Saving histograms for further processing ************************
    //******************************************************************************************                    
    TString nameOutput = Form("%s/MCInputCompilation%s_%s_%d.root", outputDir.Data(), optionPeriod.Data(), fCollisionSystenWrite.Data(), mode);
    TFile* fOutput2 = new TFile(nameOutput.Data(),"RECREATE");
    cout << "======================================================" << endl;
    cout << nameOutput << endl;
    cout << "======================================================" << endl;

        if (fEventQuality)                          fEventQuality->Write("NEvents");
        
        if (fHistoMCPi0Pt)                          fHistoMCPi0Pt->Write("MC_Pi0_Pt");
        if (fHistoMCPi0PtRebinned)                  fHistoMCPi0PtRebinned->Write("MC_Pi0_Pt_Rebinned");
        if (fHistoMCEtaPt)                          fHistoMCEtaPt->Write("MC_Eta_Pt");
        if (fHistoMCEtaPtRebinned)                  fHistoMCEtaPtRebinned->Write("MC_Eta_Pt_Rebinned");
        if (fHistoMCPiPt)                           fHistoMCPiPt->Write();
        if (fHistoMCPiNegPt)                        fHistoMCPiNegPt->Write();
        if (fHistoMCPiPosPt)                        fHistoMCPiPosPt->Write();
        if (fHistoMCKPt)                            fHistoMCKPt->Write();
        if (fHistoMCKNegPt)                         fHistoMCKNegPt->Write();
        if (fHistoMCKPosPt)                         fHistoMCKPosPt->Write();
        if (fHistoMCK0sPt)                          fHistoMCK0sPt->Write();
        if (fHistoMCK0sPtYield)                     fHistoMCK0sPtYield->Write(Form("MCYield_K0s_Pt_%1.2f",deltaRapid/2.));
        if (fHistoMCK0lPt)                          fHistoMCK0lPt->Write();
        if (fHistoMCK0lPtYield)                     fHistoMCK0lPtYield->Write(Form("MCYield_K0l_Pt_%1.2f",deltaRapid/2.));
        if (fHistoMCK0sPtRebinned)                  fHistoMCK0sPtRebinned->Write("MC_K0s_Pt_Rebinned");
        if (fHistoMCK0lPtRebinned)                  fHistoMCK0lPtRebinned->Write("MC_K0l_Pt_Rebinned");
        if (fHistoMCSecPi0FromK0sPt)                fHistoMCSecPi0FromK0sPt->Write(Form("MCSecPi0FromK0s_%1.2f",deltaRapid/2.));
        if (fHistoMCSecPi0FromK0lPt)                fHistoMCSecPi0FromK0lPt->Write(Form("MCSecPi0FromK0l_%1.2f",deltaRapid/2.));
        if (fHistoRatioMCK0sDivDataFit)             fHistoRatioMCK0sDivDataFit->Write("K0sRatioToDataFit");
        if (fHistoRatioMCKDivDataFit)               fHistoRatioMCKDivDataFit->Write("KRatioToDataFit");
        if (fHistoRatioMCPiDivDataFit)              fHistoRatioMCPiDivDataFit->Write("PiRatioToDataFit");
        if (fHistoRatioPi0FromK0sDivK0s)            fHistoRatioPi0FromK0sDivK0s->Write(Form("MCPi0FromK0sToK0s_%1.2f",deltaRapid/2.));
        if (fHistoRatioPi0FromK0lDivK0l)            fHistoRatioPi0FromK0lDivK0l->Write(Form("MCPi0FromK0lToK0l_%1.2f",deltaRapid/2.));
        if (fHistoRatioMCEtaDivPi0)                 fHistoRatioMCEtaDivPi0->Write("MCEtaToPi0");
        if (fHistoRatioMCK0sDivPi0)                 fHistoRatioMCK0sDivPi0->Write("MCK0sToPi0");
        if (fHistoRatioMCKDivPi)                    fHistoRatioMCKDivPi->Write("MCKToPi");
        
        if (fHistoMCLambdaPt)                       fHistoMCLambdaPt->Write();
        if (fHistoMCLambdaPtYield)                  fHistoMCLambdaPtYield->Write(Form("MCYield_Lambda_Pt_%1.2f",deltaRapid/2.));
        if (fHistoMCLambdaPtRebinned)               fHistoMCK0sPtRebinned->Write("MC_Lambda_Pt_Rebinned");
        if (fHistoMCSecPi0FromLambdaPt)             fHistoMCSecPi0FromLambdaPt->Write(Form("MCSecPi0FromLambda_%1.2f",deltaRapid/2.));
        if (fHistoRatioPi0FromLambdaDivLambda)      fHistoRatioPi0FromLambdaDivLambda->Write(Form("MCPi0FromLambdaToLambda_%1.2f",deltaRapid/2.));

        if (fHistoMCPiNegY)                         fHistoMCPiNegY->Write();
        if (fHistoMCPiPosY)                         fHistoMCPiPosY->Write();
        if (fHistoMCKNegY)                          fHistoMCKNegY->Write();
        if (fHistoMCKPosY)                          fHistoMCKPosY->Write();
        if (fHistoMCK0sY)                           fHistoMCK0sY->Write();        
        if (fHistoMCK0lY)                           fHistoMCK0lY->Write();
        if (fHistoMCLambdaY)                        fHistoMCLambdaY->Write();
        
    fOutput2->Write();
    fOutput2->Close();

    return;
}

