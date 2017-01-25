/****************************************************************************************************************************
******      provided by Gamma Conversion Group, PWGGA,                                                                  *****
******      Friederike Bock, friederike.bock@cern.ch                                                                    *****
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

void RemoveZerosAndPrint(TGraphAsymmErrors* graph){
  while (graph->GetY()[0] == 0. ) graph->RemovePoint(0);
  while (graph->GetY()[graph->GetN()-1] == 0. ) graph->RemovePoint(graph->GetN()-1);
  graph->Print();
}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------- Main function -----------------------------------------------------------------------------
// ----------------------- This macro is used to compile the pPb external input file for data, i.e. charged hadron, pions, kaons ... --------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------------------------------
void PrepareChargedPionDataALICE_pPb(){ 
    
    // *********************************************************************************************************************
    // ************************************** global variable definition ***************************************************
    // *********************************************************************************************************************
    TString dateForOutput                               = ReturnDateStringForOutput();

        // ************************** Low Pt  charged pions pPb************************************************
    TFile *fPionLowPtpPbMinBias = TFile::Open("ExternalInputpPb/20130201.CombinedSpectra_pA_pilot_ycms_0005_minbias.root");
    
    TH1D* histoChargedPionPlusSpecLowPtStatpPb = (TH1D*)fPionLowPtpPbMinBias->Get("stat_mb_pion_plus");
    TH1D* histoChargedPionPlusSpecLowPtSyspPb = (TH1D*)fPionLowPtpPbMinBias->Get("sys_mb_pion_plus");
    TH1D* histoChargedPionMinusSpecLowPtStatpPb = (TH1D*)fPionLowPtpPbMinBias->Get("stat_mb_pion_minus");
    TH1D* histoChargedPionMinusSpecLowPtSyspPb = (TH1D*)fPionLowPtpPbMinBias->Get("sys_mb_pion_minus");
    
    TH1D* histoChargedPionSpecLowPtStatpPb = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb->Clone("histoChargedPionSpecLowPtStatpPb");
    histoChargedPionSpecLowPtStatpPb->Add(histoChargedPionPlusSpecLowPtStatpPb);
    histoChargedPionSpecLowPtStatpPb->Scale(0.5);

    TH1D* histoChargedPionSpecLowPtSyspPb = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb->Clone("histoChargedPionSpecLowPtSyspPb");   
    
    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStatpPb->SetBinContent(i, histoChargedPionSpecLowPtStatpPb->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStatpPb->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb->GetBinError(i)/histoChargedPionSpecLowPtStatpPb->GetBinCenter(i)/(2*TMath::Pi()));
        Double_t error = 0;
        if (histoChargedPionMinusSpecLowPtSyspPb->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb->GetBinContent(i) != 0){
            error= (histoChargedPionMinusSpecLowPtSyspPb->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb->GetBinContent(i);
        } else {
            error = 0;
        }
        histoChargedPionSpecLowPtSyspPb->SetBinContent(i, histoChargedPionSpecLowPtStatpPb->GetBinContent(i));
        histoChargedPionSpecLowPtSyspPb->SetBinError(i,error);      
    }
    cout << "bis hier" << endl;
    // ************************** Low Pt  charged pions pPb depending on centrality************************************************
    TFile *fPionLowPtpPbCent = TFile::Open("ExternalInputpPb/20130201.CombinedSpectra_pA_pilot_ycms_0005.root");
    
    TH1D* histoChargedPionPlusSpecLowPtStatpPb0020 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent0_pion_plus");
    TH1D* histoChargedPionPlusSpecLowPtSyspPb0020 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent0_pion_plus");
    TH1D* histoChargedPionMinusSpecLowPtStatpPb0020 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent0_pion_minus");
    TH1D* histoChargedPionMinusSpecLowPtSyspPb0020 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent0_pion_minus");
    
    TH1D* histoChargedPionSpecLowPtStatpPb0020 = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb0020->Clone("histoChargedPionSpecLowPtStatpPb0020");
    histoChargedPionSpecLowPtStatpPb0020->Add(histoChargedPionPlusSpecLowPtStatpPb0020);
    histoChargedPionSpecLowPtStatpPb0020->Scale(0.5);

    TH1D* histoChargedPionSpecLowPtSyspPb0020 = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb0020->Clone("histoChargedPionSpecLowPtSyspPb0020");   
    
    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb0020->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStatpPb0020->SetBinContent(i, histoChargedPionSpecLowPtStatpPb0020->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb0020->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStatpPb0020->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb0020->GetBinError(i)/histoChargedPionSpecLowPtStatpPb0020->GetBinCenter(i)/(2*TMath::Pi()));
        Double_t error = 0;
        if (histoChargedPionMinusSpecLowPtSyspPb0020->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb0020->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb0020->GetBinContent(i) != 0){
            error= (histoChargedPionMinusSpecLowPtSyspPb0020->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb0020->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb0020->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb0020->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb0020->GetBinContent(i);
        } else {
            error = 0;
        }
        histoChargedPionSpecLowPtSyspPb0020->SetBinContent(i, histoChargedPionSpecLowPtStatpPb0020->GetBinContent(i));
        histoChargedPionSpecLowPtSyspPb0020->SetBinError(i,error);      
    }
    cout << "bis hier" << endl;
    TH1D* histoChargedPionPlusSpecLowPtStatpPb2040 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent1_pion_plus");
    TH1D* histoChargedPionPlusSpecLowPtSyspPb2040 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent1_pion_plus");
    TH1D* histoChargedPionMinusSpecLowPtStatpPb2040 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent1_pion_minus");
    TH1D* histoChargedPionMinusSpecLowPtSyspPb2040 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent1_pion_minus");
    
    TH1D* histoChargedPionSpecLowPtStatpPb2040 = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb2040->Clone("histoChargedPionSpecLowPtStatpPb2040");
    histoChargedPionSpecLowPtStatpPb2040->Add(histoChargedPionPlusSpecLowPtStatpPb2040);
    histoChargedPionSpecLowPtStatpPb2040->Scale(0.5);

    TH1D* histoChargedPionSpecLowPtSyspPb2040 = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb2040->Clone("histoChargedPionSpecLowPtSyspPb2040");   
    
    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb2040->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStatpPb2040->SetBinContent(i, histoChargedPionSpecLowPtStatpPb2040->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb2040->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStatpPb2040->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb2040->GetBinError(i)/histoChargedPionSpecLowPtStatpPb2040->GetBinCenter(i)/(2*TMath::Pi()));
        Double_t error = 0;
        if (histoChargedPionMinusSpecLowPtSyspPb2040->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb2040->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb2040->GetBinContent(i) != 0){
            error= (histoChargedPionMinusSpecLowPtSyspPb2040->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb2040->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb2040->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb2040->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb2040->GetBinContent(i);
        } else {
            error = 0;
        }
        histoChargedPionSpecLowPtSyspPb2040->SetBinContent(i, histoChargedPionSpecLowPtStatpPb2040->GetBinContent(i));
        histoChargedPionSpecLowPtSyspPb2040->SetBinError(i,error);      
    }
    TH1D* histoChargedPionPlusSpecLowPtStatpPb4060 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent2_pion_plus");
    TH1D* histoChargedPionPlusSpecLowPtSyspPb4060 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent2_pion_plus");
    TH1D* histoChargedPionMinusSpecLowPtStatpPb4060 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent2_pion_minus");
    TH1D* histoChargedPionMinusSpecLowPtSyspPb4060 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent2_pion_minus");
    
    TH1D* histoChargedPionSpecLowPtStatpPb4060 = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb4060->Clone("histoChargedPionSpecLowPtStatpPb4060");
    histoChargedPionSpecLowPtStatpPb4060->Add(histoChargedPionPlusSpecLowPtStatpPb4060);
    histoChargedPionSpecLowPtStatpPb4060->Scale(0.5);

    TH1D* histoChargedPionSpecLowPtSyspPb4060 = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb4060->Clone("histoChargedPionSpecLowPtSyspPb4060");   
    
    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb4060->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStatpPb4060->SetBinContent(i, histoChargedPionSpecLowPtStatpPb4060->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb4060->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStatpPb4060->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb4060->GetBinError(i)/histoChargedPionSpecLowPtStatpPb4060->GetBinCenter(i)/(2*TMath::Pi()));
        Double_t error = 0;
        if (histoChargedPionMinusSpecLowPtSyspPb4060->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb4060->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb4060->GetBinContent(i) != 0){
            error= (histoChargedPionMinusSpecLowPtSyspPb4060->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb4060->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb4060->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb4060->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb4060->GetBinContent(i);
        } else {
            error = 0;
        }
        histoChargedPionSpecLowPtSyspPb4060->SetBinContent(i, histoChargedPionSpecLowPtStatpPb4060->GetBinContent(i));
        histoChargedPionSpecLowPtSyspPb4060->SetBinError(i,error);      
    }
    cout << "bis hier" << endl;
    TH1D* histoChargedPionPlusSpecLowPtStatpPb6080 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent3_pion_plus");
    TH1D* histoChargedPionPlusSpecLowPtSyspPb6080 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent3_pion_plus");
    TH1D* histoChargedPionMinusSpecLowPtStatpPb6080 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent3_pion_minus");
    TH1D* histoChargedPionMinusSpecLowPtSyspPb6080 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent3_pion_minus");
    
    TH1D* histoChargedPionSpecLowPtStatpPb6080 = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb6080->Clone("histoChargedPionSpecLowPtStatpPb6080");
    histoChargedPionSpecLowPtStatpPb6080->Add(histoChargedPionPlusSpecLowPtStatpPb6080);
    histoChargedPionSpecLowPtStatpPb6080->Scale(0.5);

    TH1D* histoChargedPionSpecLowPtSyspPb6080 = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb6080->Clone("histoChargedPionSpecLowPtSyspPb6080");   
    
    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb6080->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStatpPb6080->SetBinContent(i, histoChargedPionSpecLowPtStatpPb6080->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb6080->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStatpPb6080->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb6080->GetBinError(i)/histoChargedPionSpecLowPtStatpPb6080->GetBinCenter(i)/(2*TMath::Pi()));
        Double_t error = 0;
        if (histoChargedPionMinusSpecLowPtSyspPb6080->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb6080->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb6080->GetBinContent(i) != 0){
            error= (histoChargedPionMinusSpecLowPtSyspPb6080->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb6080->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb6080->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb6080->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb6080->GetBinContent(i);
        } else {
            error = 0;
        }
        histoChargedPionSpecLowPtSyspPb6080->SetBinContent(i, histoChargedPionSpecLowPtStatpPb6080->GetBinContent(i));
        histoChargedPionSpecLowPtSyspPb6080->SetBinError(i,error);      
    }
    TH1D* histoChargedPionPlusSpecLowPtStatpPb80100 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent4_pion_plus");
    TH1D* histoChargedPionPlusSpecLowPtSyspPb80100 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent4_pion_plus");
    TH1D* histoChargedPionMinusSpecLowPtStatpPb80100 = (TH1D*)fPionLowPtpPbCent->Get("stat_cent4_pion_minus");
    TH1D* histoChargedPionMinusSpecLowPtSyspPb80100 = (TH1D*)fPionLowPtpPbCent->Get("sys_cent4_pion_minus");
    
    TH1D* histoChargedPionSpecLowPtStatpPb80100 = (TH1D*)histoChargedPionMinusSpecLowPtStatpPb80100->Clone("histoChargedPionSpecLowPtStatpPb80100");
    histoChargedPionSpecLowPtStatpPb80100->Add(histoChargedPionPlusSpecLowPtStatpPb80100);
    histoChargedPionSpecLowPtStatpPb80100->Scale(0.5);

    TH1D* histoChargedPionSpecLowPtSyspPb80100 = (TH1D*)histoChargedPionMinusSpecLowPtSyspPb80100->Clone("histoChargedPionSpecLowPtSyspPb80100");   
    
    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyspPb80100->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStatpPb80100->SetBinContent(i, histoChargedPionSpecLowPtStatpPb80100->GetBinContent(i)/histoChargedPionSpecLowPtStatpPb80100->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStatpPb80100->SetBinError(i, i, histoChargedPionSpecLowPtStatpPb80100->GetBinError(i)/histoChargedPionSpecLowPtStatpPb80100->GetBinCenter(i)/(2*TMath::Pi()));
        Double_t error = 0;
        if (histoChargedPionMinusSpecLowPtSyspPb80100->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyspPb80100->GetBinContent(i) != 0 && histoChargedPionSpecLowPtStatpPb80100->GetBinContent(i) != 0){
            error= (histoChargedPionMinusSpecLowPtSyspPb80100->GetBinError(i)/histoChargedPionMinusSpecLowPtSyspPb80100->GetBinContent(i) +  histoChargedPionPlusSpecLowPtSyspPb80100->GetBinError(i) /histoChargedPionPlusSpecLowPtSyspPb80100->GetBinContent(i))/2.*histoChargedPionSpecLowPtStatpPb80100->GetBinContent(i);
        } else {
            error = 0;
        }
        histoChargedPionSpecLowPtSyspPb80100->SetBinContent(i, histoChargedPionSpecLowPtStatpPb80100->GetBinContent(i));
        histoChargedPionSpecLowPtSyspPb80100->SetBinError(i,error);      
    }

    // **************************************************************************************
    // ******************************** Reading Charged Pions *******************************
    // publication: 
    // twiki: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSPECTRALowToHighPtpA
    // **************************************************************************************  
    TFile* fileIdentfiedPub                         = new TFile("ExternalInputpPb/InputRpPb/pPb502.fullpT.INEL.20151204_mb_wo_V0Acorr.root");
    TH1D* histoChargedPionspPbStatErr               = (TH1D*)fileIdentfiedPub->Get("hstat_pPb502_mb_pion_sum");
    histoChargedPionspPbStatErr->Sumw2();
    histoChargedPionspPbStatErr->Scale(0.5);
    TH1D* histoChargedPionspPbSystErr               = (TH1D*)fileIdentfiedPub->Get("hsys_pPb502_mb_pion_sum"); 
    histoChargedPionspPbSystErr->Sumw2();
    histoChargedPionspPbSystErr->Scale(0.5);
    TGraphAsymmErrors* graphChargedPionspPbSyst     = new TGraphAsymmErrors(histoChargedPionspPbStatErr);
    TGraphAsymmErrors* graphChargedPionspPbStat     = new TGraphAsymmErrors(histoChargedPionspPbSystErr);
    TH1D* histoChargedKaonspPbStatErr               = (TH1D*)fileIdentfiedPub->Get("hstat_pPb502_mb_kaon_sum");
    histoChargedKaonspPbStatErr->Sumw2();
    histoChargedKaonspPbStatErr->Scale(0.5);
    TH1D* histoChargedKaonspPbSystErr               = (TH1D*)fileIdentfiedPub->Get("hsys_pPb502_mb_kaon_sum"); 
    histoChargedKaonspPbSystErr->Sumw2();
    histoChargedKaonspPbSystErr->Scale(0.5);
    TGraphAsymmErrors* graphChargedKaonspPbSyst     = new TGraphAsymmErrors(histoChargedKaonspPbStatErr);
    TGraphAsymmErrors* graphChargedKaonspPbStat     = new TGraphAsymmErrors(histoChargedKaonspPbSystErr);
    TH1D* histoChargedProtonspPbStatErr             = (TH1D*)fileIdentfiedPub->Get("hstat_pPb502_mb_proton_sum");
    histoChargedProtonspPbStatErr->Sumw2();
    histoChargedProtonspPbStatErr->Scale(0.5);
    TH1D* histoChargedProtonspPbSystErr             = (TH1D*)fileIdentfiedPub->Get("hsys_pPb502_mb_pion_sum"); 
    histoChargedProtonspPbSystErr->Sumw2();
    histoChargedProtonspPbSystErr->Scale(0.5);
    TGraphAsymmErrors* graphChargedProtonspPbSyst   = new TGraphAsymmErrors(histoChargedProtonspPbSystErr);
    TGraphAsymmErrors* graphChargedProtonspPbStat   = new TGraphAsymmErrors(histoChargedProtonspPbStatErr);

    TH1D* histoChargedPionsStatCent[7]              = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoChargedPionsSystCent[7]              = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoChargedKaonsStatCent[7]              = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoChargedKaonsSystCent[7]              = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoChargedProtonsStatCent[7]            = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoChargedProtonsSystCent[7]            = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TString centStringCharged[7]                    = {"0005","0510","1020","2040","4060","6080","80100"};
    TFile* fileIdentfiedPubCents                    = new TFile("ExternalInputpPb/IdentifiedCharged/pPb502.fullpT.INEL.20151204.root");
    for (Int_t i = 0; i<7; i++){
        histoChargedPionsStatCent[i]                = (TH1D*)fileIdentfiedPubCents->Get(Form("hstat_pPb502_%s_pion_sum",centStringCharged[i].Data()));
        histoChargedPionsStatCent[i]->Sumw2();
        histoChargedPionsStatCent[i]->Scale(0.5);
        histoChargedPionsSystCent[i]                = (TH1D*)fileIdentfiedPubCents->Get(Form("hsys_pPb502_%s_pion_sum",centStringCharged[i].Data()));
        histoChargedPionsSystCent[i]->Sumw2();
        histoChargedPionsSystCent[i]->Scale(0.5);

        histoChargedKaonsStatCent[i]                = (TH1D*)fileIdentfiedPubCents->Get(Form("hstat_pPb502_%s_kaon_sum",centStringCharged[i].Data()));
        histoChargedKaonsStatCent[i]->Sumw2();
        histoChargedKaonsStatCent[i]->Scale(0.5);
        histoChargedKaonsSystCent[i]                = (TH1D*)fileIdentfiedPubCents->Get(Form("hsys_pPb502_%s_kaon_sum",centStringCharged[i].Data()));
        histoChargedKaonsSystCent[i]->Sumw2();
        histoChargedKaonsSystCent[i]->Scale(0.5);

        histoChargedProtonsStatCent[i]              = (TH1D*)fileIdentfiedPubCents->Get(Form("hstat_pPb502_%s_proton_sum",centStringCharged[i].Data()));
        histoChargedProtonsStatCent[i]->Sumw2();
        histoChargedProtonsStatCent[i]->Scale(0.5);
        histoChargedProtonsSystCent[i]              = (TH1D*)fileIdentfiedPubCents->Get(Form("hsys_pPb502_%s_proton_sum",centStringCharged[i].Data()));
        histoChargedProtonsSystCent[i]->Sumw2();
        histoChargedProtonsSystCent[i]->Scale(0.5);
    }
    
    // **************************************************************************************
    // ******************************** Reading Charged Pions *******************************
    // publication: 
    // twiki: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSPECTRALowToHighPtpA
    // **************************************************************************************  
    TFile* fileIdentfiedPubRpPb                     = new TFile("ExternalInputpPb/IdentifiedCharged/RpPb_502_PiKp__20151204.root");
    TH1D* histoChargedPionsRpPbStatErr              = (TH1D*)fileIdentfiedPubRpPb->Get("hstat_RpPb_pion");
    TH1D* histoChargedPionsRpPbSystErr              = (TH1D*)fileIdentfiedPubRpPb->Get("hsys_RpPb_pion"); 
    TH1D* histoChargedKaonsRpPbStatErr              = (TH1D*)fileIdentfiedPubRpPb->Get("hstat_RpPb_kaon");
    TH1D* histoChargedKaonsRpPbSystErr              = (TH1D*)fileIdentfiedPubRpPb->Get("hsys_RpPb_kaon"); 
    TH1D* histoProtonRpPbStatErr                    = (TH1D*)fileIdentfiedPubRpPb->Get("hstat_RpPb_proton");
    TH1D* histoProtonRpPbSystErr                    = (TH1D*)fileIdentfiedPubRpPb->Get("hsys_RpPb_proton"); 
    TH1D* histoChargedHadronRpPbStatErr             = (TH1D*)fileIdentfiedPubRpPb->Get("hstat_RpPb_charged"); // /HepData/8550/d4x1y1
    TH1D* histoChargedHadronRpPbSystErr             = (TH1D*)fileIdentfiedPubRpPb->Get("hsys_RpPb_charged"); // /HepData/8550/d4x1y1
    TGraphAsymmErrors* graphChargedPionsRpPbStatErr = new TGraphAsymmErrors(histoChargedPionsRpPbStatErr);
    TGraphAsymmErrors* graphChargedPionsRpPbSystErr = new TGraphAsymmErrors(histoChargedPionsRpPbSystErr);
    TGraphAsymmErrors* graphChargedKaonsRpPbStatErr = new TGraphAsymmErrors(histoChargedKaonsRpPbStatErr);
    TGraphAsymmErrors* graphChargedKaonsRpPbSystErr = new TGraphAsymmErrors(histoChargedKaonsRpPbSystErr);
    TGraphAsymmErrors* graphProtonRpPbStatErr       = new TGraphAsymmErrors(histoProtonRpPbStatErr);
    TGraphAsymmErrors* graphProtonRpPbSystErr       = new TGraphAsymmErrors(histoProtonRpPbSystErr);
    TGraphAsymmErrors* graphChargedHadronRpPbStatErr= new TGraphAsymmErrors(histoChargedHadronRpPbStatErr);
    TGraphAsymmErrors* graphChargedHadronRpPbSystErr= new TGraphAsymmErrors(histoChargedHadronRpPbSystErr);
    
    RemoveZerosAndPrint(graphChargedPionsRpPbStatErr);
    RemoveZerosAndPrint(graphChargedPionsRpPbSystErr);
    RemoveZerosAndPrint(graphChargedKaonsRpPbStatErr);
    RemoveZerosAndPrint(graphChargedKaonsRpPbSystErr);
    RemoveZerosAndPrint(graphProtonRpPbStatErr);
    RemoveZerosAndPrint(graphProtonRpPbSystErr);
    RemoveZerosAndPrint(graphChargedHadronRpPbStatErr);
    RemoveZerosAndPrint(graphChargedHadronRpPbSystErr);
    
    // *********************************************************************************************************************
    // ********************************** Write Output files ***************************************************************
    // *********************************************************************************************************************    
    TFile fileChargedPionspPb(Form("ExternalInputpPb/IdentifiedCharged/ChargedIdentifiedSpectrapPb_%s.root",dateForOutput.Data()) ,"RECREATE");
        histoChargedPionSpecLowPtSyspPb->Write("histoChargedPionSpecLowPtSyspPb", TObject::kOverwrite);
        histoChargedPionSpecLowPtStatpPb->Write("histoChargedPionSpecLowPtStatpPb", TObject::kOverwrite);
        histoChargedPionSpecLowPtSyspPb0020->Write("histoChargedPionSpecLowPtSyspPb0020", TObject::kOverwrite);
        histoChargedPionSpecLowPtStatpPb0020->Write("histoChargedPionSpecLowPtStatpPb0020", TObject::kOverwrite);
        histoChargedPionSpecLowPtSyspPb2040->Write("histoChargedPionSpecLowPtSyspPb2040", TObject::kOverwrite);
        histoChargedPionSpecLowPtStatpPb2040->Write("histoChargedPionSpecLowPtStatpPb2040", TObject::kOverwrite);
        histoChargedPionSpecLowPtSyspPb4060->Write("histoChargedPionSpecLowPtSyspPb4060", TObject::kOverwrite);
        histoChargedPionSpecLowPtStatpPb4060->Write("histoChargedPionSpecLowPtStatpPb4060", TObject::kOverwrite);
        histoChargedPionSpecLowPtSyspPb6080->Write("histoChargedPionSpecLowPtSyspPb6080", TObject::kOverwrite);
        histoChargedPionSpecLowPtStatpPb6080->Write("histoChargedPionSpecLowPtStatpPb6080", TObject::kOverwrite);
        histoChargedPionSpecLowPtSyspPb80100->Write("histoChargedPionSpecLowPtSyspPb80100", TObject::kOverwrite);
        histoChargedPionSpecLowPtStatpPb80100->Write("histoChargedPionSpecLowPtStatpPb80100", TObject::kOverwrite);
        
        if (histoChargedPionspPbStatErr) histoChargedPionspPbStatErr->Write("histoChargedPionPubStatpPb", TObject::kOverwrite);
        if (histoChargedPionspPbSystErr) histoChargedPionspPbSystErr->Write("histoChargedPionPubSyspPb", TObject::kOverwrite);
        if (histoChargedKaonspPbStatErr) histoChargedKaonspPbStatErr->Write("histoChargedKaonPubStatpPb", TObject::kOverwrite);
        if (histoChargedKaonspPbSystErr) histoChargedKaonspPbSystErr->Write("histoChargedKaonPubSyspPb", TObject::kOverwrite);
        if (histoChargedProtonspPbStatErr) histoChargedProtonspPbStatErr->Write("histoProtonPubStatpPb", TObject::kOverwrite);
        if (histoChargedProtonspPbSystErr) histoChargedProtonspPbSystErr->Write("histoProtonPubSyspPb", TObject::kOverwrite);
        if (graphChargedPionspPbSyst) graphChargedPionspPbSyst->Write("graphChargedPionPubSyspPb", TObject::kOverwrite);
        if (graphChargedPionspPbStat) graphChargedPionspPbStat->Write("graphChargedPionPubStatpPb", TObject::kOverwrite);
        if (graphChargedKaonspPbStat) graphChargedKaonspPbStat->Write("graphChargedKaonPubStatpPb", TObject::kOverwrite);
        if (graphChargedKaonspPbSyst) graphChargedKaonspPbSyst->Write("graphChargedKaonPubSyspPb", TObject::kOverwrite);
        if (graphChargedProtonspPbStat) graphChargedProtonspPbStat->Write("graphProtonPubStatpPb", TObject::kOverwrite);
        if (graphChargedProtonspPbSyst) graphChargedProtonspPbSyst->Write("graphProtonPubSyspPb", TObject::kOverwrite);
        
        for (Int_t i = 0; i < 7; i++){
            if (histoChargedPionsStatCent[i]) histoChargedPionsStatCent[i]->Write(Form("histoChargedPionPubStatpPb%s",centStringCharged[i].Data()), TObject::kOverwrite);
            if (histoChargedPionsSystCent[i]) histoChargedPionsStatCent[i]->Write(Form("histoChargedPionPubSyspPb%s",centStringCharged[i].Data())), TObject::kOverwrite;
            if (histoChargedKaonsStatCent[i]) histoChargedKaonsStatCent[i]->Write(Form("histoChargedKaonPubStatpPb%s",centStringCharged[i].Data()), TObject::kOverwrite);
            if (histoChargedKaonsSystCent[i]) histoChargedKaonsStatCent[i]->Write(Form("histoChargedKaonPubSyspPb%s",centStringCharged[i].Data()), TObject::kOverwrite);
            if (histoChargedProtonsStatCent[i]) histoChargedProtonsStatCent[i]->Write(Form("histoProtonPubStatpPb%s",centStringCharged[i].Data()), TObject::kOverwrite);
            if (histoChargedProtonsSystCent[i]) histoChargedProtonsStatCent[i]->Write(Form("histoProtonPubSyspPb%s",centStringCharged[i].Data()), TObject::kOverwrite);
        }
        
        if (histoChargedPionsRpPbStatErr) histoChargedPionsRpPbStatErr->Write("histoChargedPionPubStatpPb_RpPb");
        if (histoChargedPionsRpPbSystErr) histoChargedPionsRpPbSystErr->Write("histoChargedPionPubSyspPb_RpPb");
        if (histoChargedKaonsRpPbStatErr) histoChargedKaonsRpPbStatErr->Write("histoChargedKaonPubStatpPb_RpPb");
        if (histoChargedKaonsRpPbSystErr) histoChargedKaonsRpPbSystErr->Write("histoChargedKaonPubSyspPb_RpPb");
        if (histoProtonRpPbStatErr) histoProtonRpPbStatErr->Write("histoProtonPubStatpPb_RpPb");
        if (histoProtonRpPbSystErr) histoProtonRpPbSystErr->Write("histoProtonPubSyspPb_RpPb");
        if (histoChargedHadronRpPbStatErr) histoChargedHadronRpPbStatErr->Write("histoChargedHadronPubStatpPb_RpPb");
        if (histoChargedHadronRpPbSystErr) histoChargedHadronRpPbSystErr->Write("histoChargedHadronPubSyspPb_RpPb");
        if (graphChargedPionsRpPbStatErr) graphChargedPionsRpPbStatErr->Write("graphChargedPionPubStatpPb_RpPb");
        if (graphChargedPionsRpPbSystErr) graphChargedPionsRpPbSystErr->Write("graphChargedPionPubSyspPb_RpPb");
        if (graphChargedKaonsRpPbStatErr) graphChargedKaonsRpPbStatErr->Write("graphChargedKaonPubStatpPb_RpPb");
        if (graphChargedKaonsRpPbSystErr) graphChargedKaonsRpPbSystErr->Write("graphChargedKaonPubSyspPb_RpPb");
        if (graphProtonRpPbStatErr) graphProtonRpPbStatErr->Write("graphProtonPubStatpPb_RpPb");
        if (graphProtonRpPbSystErr) graphProtonRpPbSystErr->Write("graphProtonPubSyspPb_RpPb");
        if (graphChargedHadronRpPbStatErr) graphChargedHadronRpPbStatErr->Write("graphChargedHadronPubStatpPb_RpPb");
        if (graphChargedHadronRpPbSystErr) graphChargedHadronRpPbSystErr->Write("graphChargedHadronPubSyspPb_RpPb");
        
    fileChargedPionspPb.Close();
    
}
