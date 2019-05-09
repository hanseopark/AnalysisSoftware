/****************************************************************************************************************************
******      provided by Gamma Conversion Group, PWGGA,                                                   *****
******      Friederike Bock, friederike.bock@cern.ch                                                    *****
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
extern TBenchmark*  gBenchmark;
extern TSystem* gSystem;
extern TMinuit*     gMinuit;
//________________________________________________________________________________________________________________________
TH1D* ConvertYieldHisto(TH1D* input, Bool_t DivideBy2pi, Bool_t DivideByPt, Bool_t MultiplyBy2pi, Bool_t MultiplyByPt){
    if (!input) {
        cout << "Error: Histogram is NULL" << endl;
        return NULL;
    }

    Int_t nBins                 = input->GetNbinsX();
    Double_t newValue           = 0;
    Double_t newErrorValue      = 0;
    Double_t correctionValue    = 1;

    //correct by 2pi if specified
    if (DivideBy2pi) input->Scale(1/(2*TMath::Pi()));
    if (MultiplyBy2pi) input->Scale(2*TMath::Pi());

    for(Int_t i=0;i<nBins;i++){

        //correct by 1/Pt if specified
        if(DivideByPt)    correctionValue  = 1/(input->GetBinCenter(i+1));
        if(MultiplyByPt)  correctionValue  = input->GetBinCenter(i+1);

        //set the value and error of the bin
        input->SetBinContent(i+1,input->GetBinContent(i+1)*correctionValue);
        input->SetBinError(i+1,input->GetBinError(i+1)*correctionValue);
    }

    return input;
}

void PrepareChargedPionDataALICE_PbPb(TString energy = "PbPb_5.02TeV"){
    TDatime now;
    int iDate = now.GetDate();
    int iYear=iDate/10000;
    int iMonth=(iDate%10000)/100;
    int iDay=iDate%100;
    char* cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
            "Jul","Aug","Sep","Oct","Nov","Dec"};
    char cStamp1[25],cStamp2[25];
    sprintf(cStamp1,"%i_%s_%i",iDay, cMonth[iMonth-1], iYear);

    TString dateForOutput                       = ReturnDateStringForOutput();

  if(energy.CompareTo("PbPb_5.02TeV") == 0){
    TFile* filePbPb5TeV    = new TFile("ExternalInputPbPb/Spectra_PbPbLHC15o_Combined_07-03-19_Histograms.root");
    TList* listSummedPionStat = (TList*)filePbPb5TeV->Get("Summed_Pion");
    TList* listSummedPionSyst = (TList*)filePbPb5TeV->Get("Summed_Pion_Sys");

    TH1D*   histoChargedPionSpecStat0005  = (TH1D*)listSummedPionStat->FindObject("hSpectraSummedPion_PbPb_Combined_0.00to5.00");
    TH1D*   histoChargedPionSpecSyst0005  = (TH1D*)listSummedPionSyst->FindObject("hSpectraSummedPion_PbPb_Combined_0.00to5.00");
    histoChargedPionSpecStat0005->Sumw2();
    histoChargedPionSpecSyst0005->Sumw2();
    histoChargedPionSpecStat0005->Scale(0.5);
    histoChargedPionSpecSyst0005->Scale(0.5);
    ConvertYieldHisto(histoChargedPionSpecStat0005, kTRUE, kTRUE, kFALSE, kFALSE);
    ConvertYieldHisto(histoChargedPionSpecSyst0005, kTRUE, kTRUE, kFALSE, kFALSE);

    TH1D*   histoChargedPionSpecStat0510  = (TH1D*)listSummedPionStat->FindObject("hSpectraSummedPion_PbPb_Combined_5.00to10.00");
    TH1D*   histoChargedPionSpecSyst0510  = (TH1D*)listSummedPionSyst->FindObject("hSpectraSummedPion_PbPb_Combined_5.00to10.00");
    histoChargedPionSpecStat0510->Sumw2();
    histoChargedPionSpecSyst0510->Sumw2();
    histoChargedPionSpecStat0510->Scale(0.5);
    histoChargedPionSpecSyst0510->Scale(0.5);
    ConvertYieldHisto(histoChargedPionSpecStat0510, kTRUE, kTRUE, kFALSE, kFALSE);
    ConvertYieldHisto(histoChargedPionSpecSyst0510, kTRUE, kTRUE, kFALSE, kFALSE);

    TH1D* histoChargedPionSpecStat0010 = (TH1D*)histoChargedPionSpecStat0005->Clone("histoChargedPionSpecStat0010");
    TH1D* histoChargedPionSpecSyst0010 = (TH1D*)histoChargedPionSpecSyst0005->Clone("histoChargedPionSpecSyst0010");
    histoChargedPionSpecStat0010->Add(histoChargedPionSpecStat0510);
    histoChargedPionSpecSyst0010->Add(histoChargedPionSpecSyst0510);
    histoChargedPionSpecStat0010->Scale(0.5);
    histoChargedPionSpecSyst0010->Scale(0.5);

    TH1D*   histoChargedPionSpecStat1020  = (TH1D*)listSummedPionStat->FindObject("hSpectraSummedPion_PbPb_Combined_10.00to20.00");
    TH1D*   histoChargedPionSpecSyst1020  = (TH1D*)listSummedPionSyst->FindObject("hSpectraSummedPion_PbPb_Combined_10.00to20.00");
    histoChargedPionSpecStat1020->Sumw2();
    histoChargedPionSpecSyst1020->Sumw2();
    histoChargedPionSpecStat1020->Scale(0.5);
    histoChargedPionSpecSyst1020->Scale(0.5);
    ConvertYieldHisto(histoChargedPionSpecStat1020, kTRUE, kTRUE, kFALSE, kFALSE);
    ConvertYieldHisto(histoChargedPionSpecSyst1020, kTRUE, kTRUE, kFALSE, kFALSE);

    TH1D*   histoChargedPionSpecStat2030  = (TH1D*)listSummedPionStat->FindObject("hSpectraSummedPion_PbPb_Combined_20.00to30.00");
    TH1D*   histoChargedPionSpecSyst2030  = (TH1D*)listSummedPionSyst->FindObject("hSpectraSummedPion_PbPb_Combined_20.00to30.00");
    histoChargedPionSpecStat2030->Sumw2();
    histoChargedPionSpecSyst2030->Sumw2();
    histoChargedPionSpecStat2030->Scale(0.5);
    histoChargedPionSpecSyst2030->Scale(0.5);
    ConvertYieldHisto(histoChargedPionSpecStat2030, kTRUE, kTRUE, kFALSE, kFALSE);
    ConvertYieldHisto(histoChargedPionSpecSyst2030, kTRUE, kTRUE, kFALSE, kFALSE);

    TH1D*   histoChargedPionSpecStat3040  = (TH1D*)listSummedPionStat->FindObject("hSpectraSummedPion_PbPb_Combined_30.00to40.00");
    TH1D*   histoChargedPionSpecSyst3040  = (TH1D*)listSummedPionSyst->FindObject("hSpectraSummedPion_PbPb_Combined_30.00to40.00");
    histoChargedPionSpecStat3040->Sumw2();
    histoChargedPionSpecSyst3040->Sumw2();
    histoChargedPionSpecStat3040->Scale(0.5);
    histoChargedPionSpecSyst3040->Scale(0.5);
    ConvertYieldHisto(histoChargedPionSpecStat3040, kTRUE, kTRUE, kFALSE, kFALSE);
    ConvertYieldHisto(histoChargedPionSpecSyst3040, kTRUE, kTRUE, kFALSE, kFALSE);

    TH1D* histoChargedPionSpecStat2040 = (TH1D*)histoChargedPionSpecStat2030->Clone("histoChargedPionSpecStat2040");
    TH1D* histoChargedPionSpecSyst2040 = (TH1D*)histoChargedPionSpecSyst2030->Clone("histoChargedPionSpecSyst2040");
    histoChargedPionSpecStat2040->Add(histoChargedPionSpecStat3040);
    histoChargedPionSpecSyst2040->Add(histoChargedPionSpecSyst3040);
    histoChargedPionSpecStat2040->Scale(0.5);
    histoChargedPionSpecSyst2040->Scale(0.5);

    TH1D*   histoChargedPionSpecStat4050  = (TH1D*)listSummedPionStat->FindObject("hSpectraSummedPion_PbPb_Combined_40.00to50.00");
    TH1D*   histoChargedPionSpecSyst4050  = (TH1D*)listSummedPionSyst->FindObject("hSpectraSummedPion_PbPb_Combined_40.00to50.00");
    histoChargedPionSpecStat4050->Sumw2();
    histoChargedPionSpecSyst4050->Sumw2();
    histoChargedPionSpecStat4050->Scale(0.5);
    histoChargedPionSpecSyst4050->Scale(0.5);
    ConvertYieldHisto(histoChargedPionSpecStat4050, kTRUE, kTRUE, kFALSE, kFALSE);
    ConvertYieldHisto(histoChargedPionSpecSyst4050, kTRUE, kTRUE, kFALSE, kFALSE);

    TH1D*   histoChargedPionSpecStat5060  = (TH1D*)listSummedPionStat->FindObject("hSpectraSummedPion_PbPb_Combined_50.00to60.00");
    TH1D*   histoChargedPionSpecSyst5060  = (TH1D*)listSummedPionSyst->FindObject("hSpectraSummedPion_PbPb_Combined_50.00to60.00");
    histoChargedPionSpecStat5060->Sumw2();
    histoChargedPionSpecSyst5060->Sumw2();
    histoChargedPionSpecStat5060->Scale(0.5);
    histoChargedPionSpecSyst5060->Scale(0.5);
    ConvertYieldHisto(histoChargedPionSpecStat5060, kTRUE, kTRUE, kFALSE, kFALSE);
    ConvertYieldHisto(histoChargedPionSpecSyst5060, kTRUE, kTRUE, kFALSE, kFALSE);

    TH1D* histoChargedPionSpecStat4060 = (TH1D*)histoChargedPionSpecStat4050->Clone("histoChargedPionSpecStat4060");
    TH1D* histoChargedPionSpecSyst4060 = (TH1D*)histoChargedPionSpecSyst4050->Clone("histoChargedPionSpecSyst4060");
    histoChargedPionSpecStat4060->Add(histoChargedPionSpecStat5060);
    histoChargedPionSpecSyst4060->Add(histoChargedPionSpecSyst5060);
    histoChargedPionSpecStat4060->Scale(0.5);
    histoChargedPionSpecSyst4060->Scale(0.5);

    TH1D*   histoChargedPionSpecStat6070  = (TH1D*)listSummedPionStat->FindObject("hSpectraSummedPion_PbPb_Combined_60.00to70.00");
    TH1D*   histoChargedPionSpecSyst6070  = (TH1D*)listSummedPionSyst->FindObject("hSpectraSummedPion_PbPb_Combined_60.00to70.00");
    histoChargedPionSpecStat6070->Sumw2();
    histoChargedPionSpecSyst6070->Sumw2();
    histoChargedPionSpecStat6070->Scale(0.5);
    histoChargedPionSpecSyst6070->Scale(0.5);
    ConvertYieldHisto(histoChargedPionSpecStat6070, kTRUE, kTRUE, kFALSE, kFALSE);
    ConvertYieldHisto(histoChargedPionSpecSyst6070, kTRUE, kTRUE, kFALSE, kFALSE);

    TH1D*   histoChargedPionSpecStat7080  = (TH1D*)listSummedPionStat->FindObject("hSpectraSummedPion_PbPb_Combined_70.00to80.00");
    TH1D*   histoChargedPionSpecSyst7080  = (TH1D*)listSummedPionSyst->FindObject("hSpectraSummedPion_PbPb_Combined_70.00to80.00");
    histoChargedPionSpecStat7080->Sumw2();
    histoChargedPionSpecSyst7080->Sumw2();
    histoChargedPionSpecStat7080->Scale(0.5);
    histoChargedPionSpecSyst7080->Scale(0.5);
    ConvertYieldHisto(histoChargedPionSpecStat7080, kTRUE, kTRUE, kFALSE, kFALSE);
    ConvertYieldHisto(histoChargedPionSpecSyst7080, kTRUE, kTRUE, kFALSE, kFALSE);

    TH1D* histoChargedPionSpecStat6080 = (TH1D*)histoChargedPionSpecStat6070->Clone("histoChargedPionSpecStat6080");
    TH1D* histoChargedPionSpecSyst6080 = (TH1D*)histoChargedPionSpecSyst6070->Clone("histoChargedPionSpecSyst6080");
    histoChargedPionSpecStat6080->Add(histoChargedPionSpecStat7080);
    histoChargedPionSpecSyst6080->Add(histoChargedPionSpecSyst7080);
    histoChargedPionSpecStat6080->Scale(0.5);
    histoChargedPionSpecSyst6080->Scale(0.5);

    TH1D*   histoChargedPionSpecStat8090  = (TH1D*)listSummedPionStat->FindObject("hSpectraSummedPion_PbPb_Combined_80.00to90.00");
    TH1D*   histoChargedPionSpecSyst8090  = (TH1D*)listSummedPionSyst->FindObject("hSpectraSummedPion_PbPb_Combined_80.00to90.00");
    histoChargedPionSpecStat8090->Sumw2();
    histoChargedPionSpecSyst8090->Sumw2();
    histoChargedPionSpecStat8090->Scale(0.5);
    histoChargedPionSpecSyst8090->Scale(0.5);
    ConvertYieldHisto(histoChargedPionSpecStat8090, kTRUE, kTRUE, kFALSE, kFALSE);
    ConvertYieldHisto(histoChargedPionSpecSyst8090, kTRUE, kTRUE, kFALSE, kFALSE);

    TString outputFileName = Form("ChargedPion_%s_PbPb5TeV.root",cStamp1);
    TFile* fileOutput = new TFile(outputFileName,"RECREATE");

    TString nameCentFolder_0_10 = "0-10%PbPb_5.02TeV";
    fileOutput->mkdir(nameCentFolder_0_10.Data());
    TDirectoryFile* directoryPi0Output_0_10 = (TDirectoryFile*)fileOutput->Get(nameCentFolder_0_10.Data());
    fileOutput->cd(nameCentFolder_0_10.Data());
    if(histoChargedPionSpecStat0010)  histoChargedPionSpecStat0010->Write("histoChargedPionSpecStat");
    if(histoChargedPionSpecSyst0010)  histoChargedPionSpecSyst0010->Write("histoChargedPionSpecSyst");

    TString nameCentFolder_10_20 = "10-20%PbPb_5.02TeV";
    fileOutput->mkdir(nameCentFolder_10_20.Data());
    TDirectoryFile* directoryPi0Output_10_20 = (TDirectoryFile*)fileOutput->Get(nameCentFolder_10_20.Data());
    fileOutput->cd(nameCentFolder_10_20.Data());
    if(histoChargedPionSpecStat1020)  histoChargedPionSpecStat1020->Write("histoChargedPionSpecStat");
    if(histoChargedPionSpecSyst1020)  histoChargedPionSpecSyst1020->Write("histoChargedPionSpecSyst");

    TString nameCentFolder_20_40 = "20-40%PbPb_5.02TeV";
    fileOutput->mkdir(nameCentFolder_20_40.Data());
    TDirectoryFile* directoryPi0Output_20_40 = (TDirectoryFile*)fileOutput->Get(nameCentFolder_20_40.Data());
    fileOutput->cd(nameCentFolder_20_40.Data());
    if(histoChargedPionSpecStat2040)  histoChargedPionSpecStat2040->Write("histoChargedPionSpecStat");
    if(histoChargedPionSpecSyst2040)  histoChargedPionSpecSyst2040->Write("histoChargedPionSpecSyst");

    TString nameCentFolder_40_60 = "40-60%PbPb_5.02TeV";
    fileOutput->mkdir(nameCentFolder_40_60.Data());
    TDirectoryFile* directoryPi0Output_40_60 = (TDirectoryFile*)fileOutput->Get(nameCentFolder_40_60.Data());
    fileOutput->cd(nameCentFolder_40_60.Data());
    if(histoChargedPionSpecStat4060)  histoChargedPionSpecStat4060->Write("histoChargedPionSpecStat");
    if(histoChargedPionSpecSyst4060)  histoChargedPionSpecSyst4060->Write("histoChargedPionSpecSyst");

    TString nameCentFolder_60_80 = "60-80%PbPb_5.02TeV";
    fileOutput->mkdir(nameCentFolder_60_80.Data());
    TDirectoryFile* directoryPi0Output_60_80 = (TDirectoryFile*)fileOutput->Get(nameCentFolder_60_80.Data());
    fileOutput->cd(nameCentFolder_60_80.Data());
    if(histoChargedPionSpecStat6080)  histoChargedPionSpecStat6080->Write("histoChargedPionSpecStat");
    if(histoChargedPionSpecSyst6080)  histoChargedPionSpecSyst6080->Write("histoChargedPionSpecSyst");

    fileOutput->Write();
    fileOutput->Close();

  }else{
    
    //***********************************high pt spectra ********************************************************************
    TFile* fileChargedPionSpectraHighPtFinal    = new TFile("ExternalInputPbPb/IdentifiedCharged/PbPb276.fullpT.INEL.20140329.root");
    
    TH1D*   histoChargedPionSpecHighPtStat0005  = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrum_0_5");
    TH1D*   histoChargedPionSpecHighPtSyst0005  = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrumSyst_0_5");
    histoChargedPionSpecHighPtStat0005->Scale(0.5);
    histoChargedPionSpecHighPtSyst0005->Scale(0.5);

    TH1D*   histoChargedPionSpecHighPtStat0510 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrum_5_10");
    TH1D*   histoChargedPionSpecHighPtSyst0510 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrumSyst_5_10");
    histoChargedPionSpecHighPtStat0510->Scale(0.5);
    histoChargedPionSpecHighPtSyst0510->Scale(0.5);

    TH1D*   histoChargedPionSpecHighPtStat0010 = (TH1D*)histoChargedPionSpecHighPtStat0510->Clone("histoChargedPionSpecHighPtStat0010");
    TH1D*   histoChargedPionSpecHighPtSyst0010 = (TH1D*)histoChargedPionSpecHighPtSyst0510->Clone("histoChargedPionSpecHighPtSyst0010");
    histoChargedPionSpecHighPtStat0010->Add(histoChargedPionSpecHighPtStat0005);
    histoChargedPionSpecHighPtSyst0010->Add(histoChargedPionSpecHighPtSyst0005);
    histoChargedPionSpecHighPtStat0010->Scale(0.5);
    histoChargedPionSpecHighPtSyst0010->Scale(0.5);

    for (Int_t i = 1; i < histoChargedPionSpecHighPtSyst0010->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0;
        if (histoChargedPionSpecHighPtSyst0005->GetBinContent(i) != 0){
            relErrLowerCent= histoChargedPionSpecHighPtSyst0005->GetBinError(i)/histoChargedPionSpecHighPtSyst0005->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent = 0;
        if (histoChargedPionSpecHighPtSyst0510->GetBinContent(i) != 0){
            relErrHigherCent = histoChargedPionSpecHighPtSyst0510->GetBinError(i)/histoChargedPionSpecHighPtSyst0510->GetBinContent(i)*100 ;
        }
        
        if (relErrHigherCent > relErrLowerCent){
            histoChargedPionSpecHighPtSyst0010->SetBinError(i, histoChargedPionSpecHighPtSyst0010->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoChargedPionSpecHighPtSyst0010->SetBinError(i, histoChargedPionSpecHighPtSyst0010->GetBinContent(i)*relErrLowerCent/100);
        }         
    }   


    TH1D*   histoChargedPionSpecHighPtStat1020 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrum_10_20");
    TH1D*   histoChargedPionSpecHighPtSyst1020 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrumSyst_10_20");
    histoChargedPionSpecHighPtStat1020->Scale(0.5);
    histoChargedPionSpecHighPtSyst1020->Scale(0.5);
    
    TH1D*   histoChargedPionSpecHighPtStat2040 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrum_20_40");
    TH1D*   histoChargedPionSpecHighPtSyst2040 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrumSyst_20_40");
    histoChargedPionSpecHighPtStat2040->Scale(0.5);
    histoChargedPionSpecHighPtSyst2040->Scale(0.5);
    
    TH1D*   histoChargedPionSpecHighPtStat4060 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrum_40_60");
    TH1D*   histoChargedPionSpecHighPtSyst4060 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrumSyst_40_60");
    histoChargedPionSpecHighPtStat4060->Scale(0.5);
    histoChargedPionSpecHighPtSyst4060->Scale(0.5);
    
    TH1D*   histoChargedPionSpecHighPtStat6080 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrum_60_80");
    TH1D*   histoChargedPionSpecHighPtSyst6080 = (TH1D*)fileChargedPionSpectraHighPtFinal->Get("hPionSpectrumSyst_60_80");
    histoChargedPionSpecHighPtStat6080->Scale(0.5);
    histoChargedPionSpecHighPtSyst6080->Scale(0.5);
    
    //***********************************high pt spectra ********************************************************************
    TFile* fileChargedKaonSpectraHighPtFinal = new TFile("ExternalInputPbPb/IdentifiedCharged/SpectraHighPtKaonFinal_20131108.root");
        
    TH1D* histoChargedKaonSpecHighPtStat0005 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrum_0_5");
    TH1D* histoChargedKaonSpecHighPtSyst0005 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrumSyst_0_5");
    histoChargedKaonSpecHighPtStat0005->Scale(0.5);
    histoChargedKaonSpecHighPtSyst0005->Scale(0.5);

    TH1D* histoChargedKaonSpecHighPtStat0510 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrum_5_10");
    TH1D* histoChargedKaonSpecHighPtSyst0510 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrumSyst_5_10");
    histoChargedKaonSpecHighPtStat0510->Scale(0.5);
    histoChargedKaonSpecHighPtSyst0510->Scale(0.5);

    TH1D* histoChargedKaonSpecHighPtStat0010 = (TH1D*)histoChargedKaonSpecHighPtStat0510->Clone("histoChargedKaonSpecHighPtStat0010");
    TH1D* histoChargedKaonSpecHighPtSyst0010 = (TH1D*)histoChargedKaonSpecHighPtSyst0510->Clone("histoChargedKaonSpecHighPtSyst0010");
    histoChargedKaonSpecHighPtStat0010->Add(histoChargedKaonSpecHighPtStat0005);
    histoChargedKaonSpecHighPtSyst0010->Add(histoChargedKaonSpecHighPtSyst0005);
    histoChargedKaonSpecHighPtStat0010->Scale(0.5);
    histoChargedKaonSpecHighPtSyst0010->Scale(0.5);

    for (Int_t i = 1; i < histoChargedKaonSpecHighPtSyst0010->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0;
        if (histoChargedKaonSpecHighPtSyst0005->GetBinContent(i) != 0){
            relErrLowerCent= histoChargedKaonSpecHighPtSyst0005->GetBinError(i)/histoChargedKaonSpecHighPtSyst0005->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent = 0;
        if (histoChargedKaonSpecHighPtSyst0510->GetBinContent(i) != 0){
            relErrHigherCent = histoChargedKaonSpecHighPtSyst0510->GetBinError(i)/histoChargedKaonSpecHighPtSyst0510->GetBinContent(i)*100 ;
        }
        
        if (relErrHigherCent > relErrLowerCent){
            histoChargedKaonSpecHighPtSyst0010->SetBinError(i, histoChargedKaonSpecHighPtSyst0010->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoChargedKaonSpecHighPtSyst0010->SetBinError(i, histoChargedKaonSpecHighPtSyst0010->GetBinContent(i)*relErrLowerCent/100);
        }         
    }   

    TH1D* histoChargedKaonSpecHighPtStat1020 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrum_10_20");
    TH1D* histoChargedKaonSpecHighPtSyst1020 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrumSyst_10_20");
    histoChargedKaonSpecHighPtStat1020->Scale(0.5);
    histoChargedKaonSpecHighPtSyst1020->Scale(0.5);
    
    TH1D* histoChargedKaonSpecHighPtStat2040 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrum_20_40");
    TH1D* histoChargedKaonSpecHighPtSyst2040 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrumSyst_20_40");
    histoChargedKaonSpecHighPtStat2040->Scale(0.5);
    histoChargedKaonSpecHighPtSyst2040->Scale(0.5);
    
    TH1D* histoChargedKaonSpecHighPtStat4060 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrum_40_60");
    TH1D* histoChargedKaonSpecHighPtSyst4060 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrumSyst_40_60");
    histoChargedKaonSpecHighPtStat4060->Scale(0.5);
    histoChargedKaonSpecHighPtSyst4060->Scale(0.5);
    
    TH1D* histoChargedKaonSpecHighPtStat6080 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrum_60_80");
    TH1D* histoChargedKaonSpecHighPtSyst6080 = (TH1D*)fileChargedKaonSpectraHighPtFinal->Get("hKaonSpectrumSyst_60_80");
    histoChargedKaonSpecHighPtStat6080->Scale(0.5);
    histoChargedKaonSpecHighPtSyst6080->Scale(0.5);
    
    TFile* fileChargedIndentifiedSpectraLowPtPrelim2012 = new TFile("ExternalInputPbPb/IdentifiedCharged/SPECTRA_COMB_20120809_default.root");
    TH1D*   histoChargedPionMinusSpecLowPtStat0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent0_pion_minus");
    histoChargedPionMinusSpecLowPtStat0005->Sumw2();
    TH1D*   histoChargedPionMinusSpecLowPtSyst0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent0_pion_minus");
    histoChargedPionMinusSpecLowPtSyst0005->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtStat0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent0_pion_plus");
    histoChargedPionPlusSpecLowPtStat0005->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtSyst0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent0_pion_plus");
    histoChargedPionPlusSpecLowPtSyst0005->Sumw2();
    TH1D*   histoChargedPionSpecLowPtStat0005 = (TH1D*)histoChargedPionMinusSpecLowPtStat0005->Clone("histoChargedPionSpecLowPtStat0005");
    TH1D*   histoChargedPionSpecLowPtSyst0005 = (TH1D*)histoChargedPionMinusSpecLowPtSyst0005->Clone("histoChargedPionSpecLowPtSyst0005");
    histoChargedPionSpecLowPtStat0005->Add(histoChargedPionPlusSpecLowPtStat0005);
    histoChargedPionSpecLowPtSyst0005->Add(histoChargedPionPlusSpecLowPtSyst0005);
    histoChargedPionSpecLowPtStat0005->Scale(0.5);
    histoChargedPionSpecLowPtSyst0005->Scale(0.5);
    
    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst0005->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStat0005->SetBinContent(i, histoChargedPionSpecLowPtStat0005->GetBinContent(i)/histoChargedPionSpecLowPtStat0005->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStat0005->SetBinError(i, histoChargedPionSpecLowPtStat0005->GetBinError(i)/histoChargedPionSpecLowPtStat0005->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtSyst0005->SetBinContent(i, histoChargedPionSpecLowPtSyst0005->GetBinContent(i)/(histoChargedPionSpecLowPtSyst0005->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedPionPlusSpecLowPtSyst0005->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst0005->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst0005->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst0005->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst0005->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst0005->GetBinContent(i)*100)/2;
        }
        histoChargedPionSpecLowPtSyst0005->SetBinError(i, histoChargedPionSpecLowPtSyst0005->GetBinContent(i)*fractionalSystematicError/100.);
    }   

    TH1D*   histoChargedPionMinusSpecLowPtStat0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent1_pion_minus");
    histoChargedPionMinusSpecLowPtStat0510->Sumw2();
    TH1D*   histoChargedPionMinusSpecLowPtSyst0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent1_pion_minus");
    histoChargedPionMinusSpecLowPtSyst0510->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtStat0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent1_pion_plus");
    histoChargedPionPlusSpecLowPtStat0510->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtSyst0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent1_pion_plus");
    histoChargedPionPlusSpecLowPtSyst0510->Sumw2();
    TH1D*   histoChargedPionSpecLowPtStat0510 = (TH1D*)histoChargedPionMinusSpecLowPtStat0510->Clone("histoChargedPionSpecLowPtStat0510");
    TH1D*   histoChargedPionSpecLowPtSyst0510 = (TH1D*)histoChargedPionMinusSpecLowPtSyst0510->Clone("histoChargedPionSpecLowPtSyst0510");
    histoChargedPionSpecLowPtStat0510->Add(histoChargedPionPlusSpecLowPtStat0510);
    histoChargedPionSpecLowPtSyst0510->Add(histoChargedPionPlusSpecLowPtSyst0510);
    histoChargedPionSpecLowPtStat0510->Scale(0.5);
    histoChargedPionSpecLowPtSyst0510->Scale(0.5);

    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst0510->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStat0510->SetBinContent(i, histoChargedPionSpecLowPtStat0510->GetBinContent(i)/histoChargedPionSpecLowPtStat0510->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStat0510->SetBinError(i, histoChargedPionSpecLowPtStat0510->GetBinError(i)/histoChargedPionSpecLowPtStat0510->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtSyst0510->SetBinContent(i, histoChargedPionSpecLowPtSyst0510->GetBinContent(i)/(histoChargedPionSpecLowPtSyst0510->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedPionPlusSpecLowPtSyst0510->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst0510->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst0510->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst0510->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst0510->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst0510->GetBinContent(i)*100)/2;
        }
        histoChargedPionSpecLowPtSyst0510->SetBinError(i, histoChargedPionSpecLowPtSyst0510->GetBinContent(i)*fractionalSystematicError/100.);
    }
    
    TH1D*   histoChargedPionSpecLowPtStat0010 = (TH1D*)histoChargedPionSpecLowPtStat0510->Clone("histoChargedPionSpecLowPtStat0010");
    TH1D*   histoChargedPionSpecLowPtSyst0010 = (TH1D*)histoChargedPionSpecLowPtSyst0510->Clone("histoChargedPionSpecLowPtSyst0010");
    histoChargedPionSpecLowPtStat0010->Add(histoChargedPionSpecLowPtStat0005);
    histoChargedPionSpecLowPtSyst0010->Add(histoChargedPionSpecLowPtSyst0005);
    histoChargedPionSpecLowPtStat0010->Scale(0.5);
    histoChargedPionSpecLowPtSyst0010->Scale(0.5);

    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst0010->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0. ;
        if (histoChargedPionSpecLowPtSyst0005->GetBinContent(i) != 0){
            relErrLowerCent = histoChargedPionSpecLowPtSyst0005->GetBinError(i)/histoChargedPionSpecLowPtSyst0005->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent = 0. ;
        if (histoChargedPionSpecLowPtSyst0510->GetBinContent(i) != 0){
            relErrHigherCent = histoChargedPionSpecLowPtSyst0510->GetBinError(i)/histoChargedPionSpecLowPtSyst0510->GetBinContent(i)*100 ;
        }
        if (relErrHigherCent > relErrLowerCent){
            histoChargedPionSpecLowPtSyst0010->SetBinError(i, histoChargedPionSpecLowPtSyst0010->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoChargedPionSpecLowPtSyst0010->SetBinError(i, histoChargedPionSpecLowPtSyst0010->GetBinContent(i)*relErrLowerCent/100);
        }         
    }   

    TH1D*   histoChargedPionMinusSpecLowPtStat1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent2_pion_minus");
    histoChargedPionMinusSpecLowPtStat1020->Sumw2();
    TH1D*   histoChargedPionMinusSpecLowPtSyst1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent2_pion_minus");
    histoChargedPionMinusSpecLowPtSyst1020->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtStat1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent2_pion_plus");
    histoChargedPionPlusSpecLowPtStat1020->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtSyst1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent2_pion_plus");
    histoChargedPionPlusSpecLowPtSyst1020->Sumw2();
    TH1D*   histoChargedPionSpecLowPtStat1020 = (TH1D*)histoChargedPionMinusSpecLowPtStat1020->Clone("histoChargedPionSpecLowPtStat1020");
    TH1D*   histoChargedPionSpecLowPtSyst1020 = (TH1D*)histoChargedPionMinusSpecLowPtSyst1020->Clone("histoChargedPionSpecLowPtSyst1020");
    histoChargedPionSpecLowPtStat1020->Add(histoChargedPionPlusSpecLowPtStat1020);
    histoChargedPionSpecLowPtSyst1020->Add(histoChargedPionPlusSpecLowPtSyst1020);
    histoChargedPionSpecLowPtStat1020->Scale(0.5);
    histoChargedPionSpecLowPtSyst1020->Scale(0.5);
    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst1020->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStat1020->SetBinContent(i, histoChargedPionSpecLowPtStat1020->GetBinContent(i)/histoChargedPionSpecLowPtStat1020->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStat1020->SetBinError(i, histoChargedPionSpecLowPtStat1020->GetBinError(i)/histoChargedPionSpecLowPtStat1020->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtSyst1020->SetBinContent(i, histoChargedPionSpecLowPtSyst1020->GetBinContent(i)/(histoChargedPionSpecLowPtSyst1020->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedPionPlusSpecLowPtSyst1020->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst1020->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst1020->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst1020->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst1020->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst1020->GetBinContent(i)*100)/2;
        }
        histoChargedPionSpecLowPtSyst1020->SetBinError(i, histoChargedPionSpecLowPtSyst1020->GetBinContent(i)*fractionalSystematicError/100.);
    }
    
    TH1D*   histoChargedPionMinusSpecLowPtStat2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent3_pion_minus");
    histoChargedPionMinusSpecLowPtStat2030->Sumw2();
    TH1D*   histoChargedPionMinusSpecLowPtSyst2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent3_pion_minus");
    histoChargedPionMinusSpecLowPtSyst2030->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtStat2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent3_pion_plus");
    histoChargedPionPlusSpecLowPtStat2030->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtSyst2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent3_pion_plus");
    histoChargedPionPlusSpecLowPtSyst2030->Sumw2();
    TH1D*   histoChargedPionSpecLowPtStat2030 = (TH1D*)histoChargedPionMinusSpecLowPtStat2030->Clone("histoChargedPionSpecLowPtStat2030");
    TH1D*   histoChargedPionSpecLowPtSyst2030 = (TH1D*)histoChargedPionMinusSpecLowPtSyst2030->Clone("histoChargedPionSpecLowPtSyst2030");
    histoChargedPionSpecLowPtStat2030->Add(histoChargedPionPlusSpecLowPtStat2030);
    histoChargedPionSpecLowPtSyst2030->Add(histoChargedPionPlusSpecLowPtSyst2030);
    histoChargedPionSpecLowPtStat2030->Scale(0.5);
    histoChargedPionSpecLowPtSyst2030->Scale(0.5);
    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst2030->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStat2030->SetBinContent(i, histoChargedPionSpecLowPtStat2030->GetBinContent(i)/histoChargedPionSpecLowPtStat2030->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStat2030->SetBinError(i, histoChargedPionSpecLowPtStat2030->GetBinError(i)/histoChargedPionSpecLowPtStat2030->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtSyst2030->SetBinContent(i, histoChargedPionSpecLowPtSyst2030->GetBinContent(i)/(histoChargedPionSpecLowPtSyst2030->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedPionPlusSpecLowPtSyst2030->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst2030->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst2030->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst2030->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst2030->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst2030->GetBinContent(i)*100)/2;
        }
        histoChargedPionSpecLowPtSyst2030->SetBinError(i, histoChargedPionSpecLowPtSyst2030->GetBinContent(i)*fractionalSystematicError/100.);
    }
    TH1D*   histoChargedPionMinusSpecLowPtStat3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent4_pion_minus");
    histoChargedPionMinusSpecLowPtStat3040->Sumw2();
    TH1D*   histoChargedPionMinusSpecLowPtSyst3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent4_pion_minus");
    histoChargedPionMinusSpecLowPtSyst3040->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtStat3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent4_pion_plus");
    histoChargedPionPlusSpecLowPtStat3040->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtSyst3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent4_pion_plus");
    histoChargedPionPlusSpecLowPtSyst3040->Sumw2();
    TH1D*   histoChargedPionSpecLowPtStat3040 = (TH1D*)histoChargedPionMinusSpecLowPtStat3040->Clone("histoChargedPionSpecLowPtStat3040");
    TH1D*   histoChargedPionSpecLowPtSyst3040 = (TH1D*)histoChargedPionMinusSpecLowPtSyst3040->Clone("histoChargedPionSpecLowPtSyst3040");
    histoChargedPionSpecLowPtStat3040->Add(histoChargedPionPlusSpecLowPtStat3040);
    histoChargedPionSpecLowPtSyst3040->Add(histoChargedPionPlusSpecLowPtSyst3040);
    histoChargedPionSpecLowPtStat3040->Scale(0.5);
    histoChargedPionSpecLowPtSyst3040->Scale(0.5);
    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst3040->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStat3040->SetBinContent(i, histoChargedPionSpecLowPtStat3040->GetBinContent(i)/histoChargedPionSpecLowPtStat3040->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStat3040->SetBinError(i, histoChargedPionSpecLowPtStat3040->GetBinError(i)/histoChargedPionSpecLowPtStat3040->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtSyst3040->SetBinContent(i, histoChargedPionSpecLowPtSyst3040->GetBinContent(i)/(histoChargedPionSpecLowPtSyst3040->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedPionPlusSpecLowPtSyst3040->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst3040->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst3040->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst3040->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst3040->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst3040->GetBinContent(i)*100)/2;
        }
        histoChargedPionSpecLowPtSyst3040->SetBinError(i, histoChargedPionSpecLowPtSyst3040->GetBinContent(i)*fractionalSystematicError/100.);
    }
    TH1D*   histoChargedPionSpecLowPtStat2040 = (TH1D*)histoChargedPionSpecLowPtStat2030->Clone("histoChargedPionSpecLowPtStat2040");
    TH1D*   histoChargedPionSpecLowPtSyst2040 = (TH1D*)histoChargedPionSpecLowPtSyst2030->Clone("histoChargedPionSpecLowPtSyst2040");
    histoChargedPionSpecLowPtStat2040->Add(histoChargedPionSpecLowPtStat3040);
    histoChargedPionSpecLowPtSyst2040->Add(histoChargedPionSpecLowPtSyst3040);
    histoChargedPionSpecLowPtStat2040->Scale(0.5);
    histoChargedPionSpecLowPtSyst2040->Scale(0.5);

for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst2040->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0;
        if (histoChargedPionSpecLowPtSyst2030->GetBinContent(i) != 0){
            relErrLowerCent = histoChargedPionSpecLowPtSyst2030->GetBinError(i)/histoChargedPionSpecLowPtSyst2030->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent = 0;
        if (histoChargedPionSpecLowPtSyst3040->GetBinContent(i) != 0){
            relErrHigherCent = histoChargedPionSpecLowPtSyst3040->GetBinError(i)/histoChargedPionSpecLowPtSyst3040->GetBinContent(i)*100 ;
        }
        
        if (relErrHigherCent > relErrLowerCent){
            histoChargedPionSpecLowPtSyst2040->SetBinError(i, histoChargedPionSpecLowPtSyst2040->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoChargedPionSpecLowPtSyst2040->SetBinError(i, histoChargedPionSpecLowPtSyst2040->GetBinContent(i)*relErrLowerCent/100);
        }         
    }   
    
    TH1D*   histoChargedPionMinusSpecLowPtStat4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent5_pion_minus");
    histoChargedPionMinusSpecLowPtStat4050->Sumw2();
    TH1D*   histoChargedPionMinusSpecLowPtSyst4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent5_pion_minus");
    histoChargedPionMinusSpecLowPtSyst4050->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtStat4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent5_pion_plus");
    histoChargedPionPlusSpecLowPtStat4050->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtSyst4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent5_pion_plus");
    histoChargedPionPlusSpecLowPtSyst4050->Sumw2();
    TH1D*   histoChargedPionSpecLowPtStat4050 = (TH1D*)histoChargedPionMinusSpecLowPtStat4050->Clone("histoChargedPionSpecLowPtStat4050");
    TH1D*   histoChargedPionSpecLowPtSyst4050 = (TH1D*)histoChargedPionMinusSpecLowPtSyst4050->Clone("histoChargedPionSpecLowPtSyst4050");
    histoChargedPionSpecLowPtStat4050->Add(histoChargedPionPlusSpecLowPtStat4050);
    histoChargedPionSpecLowPtSyst4050->Add(histoChargedPionPlusSpecLowPtSyst4050);
    histoChargedPionSpecLowPtStat4050->Scale(0.5);
    histoChargedPionSpecLowPtSyst4050->Scale(0.5);
    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst4050->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStat4050->SetBinContent(i, histoChargedPionSpecLowPtStat4050->GetBinContent(i)/histoChargedPionSpecLowPtStat4050->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStat4050->SetBinError(i, histoChargedPionSpecLowPtStat4050->GetBinError(i)/histoChargedPionSpecLowPtStat4050->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtSyst4050->SetBinContent(i, histoChargedPionSpecLowPtSyst4050->GetBinContent(i)/(histoChargedPionSpecLowPtSyst4050->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedPionPlusSpecLowPtSyst4050->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst4050->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst4050->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst4050->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst4050->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst4050->GetBinContent(i)*100)/2;
        }
        histoChargedPionSpecLowPtSyst4050->SetBinError(i, histoChargedPionSpecLowPtSyst4050->GetBinContent(i)*fractionalSystematicError/100.);
    }

    TH1D*   histoChargedPionMinusSpecLowPtStat5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent6_pion_minus");
    histoChargedPionMinusSpecLowPtStat5060->Sumw2();
    TH1D*   histoChargedPionMinusSpecLowPtSyst5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent6_pion_minus");
    histoChargedPionMinusSpecLowPtSyst5060->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtStat5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent6_pion_plus");
    histoChargedPionPlusSpecLowPtStat5060->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtSyst5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent6_pion_plus");
    histoChargedPionPlusSpecLowPtSyst5060->Sumw2();
    TH1D*   histoChargedPionSpecLowPtStat5060 = (TH1D*)histoChargedPionMinusSpecLowPtStat5060->Clone("histoChargedPionSpecLowPtStat5060");
    TH1D*   histoChargedPionSpecLowPtSyst5060 = (TH1D*)histoChargedPionMinusSpecLowPtSyst5060->Clone("histoChargedPionSpecLowPtSyst5060");
    histoChargedPionSpecLowPtStat5060->Add(histoChargedPionPlusSpecLowPtStat5060);
    histoChargedPionSpecLowPtSyst5060->Add(histoChargedPionPlusSpecLowPtSyst5060);
    histoChargedPionSpecLowPtStat5060->Scale(0.5);
    histoChargedPionSpecLowPtSyst5060->Scale(0.5);
    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst5060->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStat5060->SetBinContent(i, histoChargedPionSpecLowPtStat5060->GetBinContent(i)/histoChargedPionSpecLowPtStat5060->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStat5060->SetBinError(i, histoChargedPionSpecLowPtStat5060->GetBinError(i)/histoChargedPionSpecLowPtStat5060->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtSyst5060->SetBinContent(i, histoChargedPionSpecLowPtSyst5060->GetBinContent(i)/(histoChargedPionSpecLowPtSyst5060->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedPionPlusSpecLowPtSyst5060->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst5060->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst5060->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst5060->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst5060->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst5060->GetBinContent(i)*100)/2;
        }
        histoChargedPionSpecLowPtSyst5060->SetBinError(i, histoChargedPionSpecLowPtSyst5060->GetBinContent(i)*fractionalSystematicError/100.);
    }
    TH1D*   histoChargedPionSpecLowPtStat4060 = (TH1D*)histoChargedPionSpecLowPtStat4050->Clone("histoChargedPionSpecLowPtStat4060");
    TH1D*   histoChargedPionSpecLowPtSyst4060 = (TH1D*)histoChargedPionSpecLowPtSyst4050->Clone("histoChargedPionSpecLowPtSyst4060");
    histoChargedPionSpecLowPtStat4060->Add(histoChargedPionSpecLowPtStat5060);
    histoChargedPionSpecLowPtSyst4060->Add(histoChargedPionSpecLowPtSyst5060);
    histoChargedPionSpecLowPtStat4060->Scale(0.5);
    histoChargedPionSpecLowPtSyst4060->Scale(0.5);

    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst4060->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0.;
        if (histoChargedPionSpecLowPtSyst4050->GetBinContent(i) != 0){
            relErrLowerCent = histoChargedPionSpecLowPtSyst4050->GetBinError(i)/histoChargedPionSpecLowPtSyst4050->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent = 0.;
        if (histoChargedPionSpecLowPtSyst5060->GetBinContent(i) != 0){
            relErrHigherCent = histoChargedPionSpecLowPtSyst5060->GetBinError(i)/histoChargedPionSpecLowPtSyst5060->GetBinContent(i)*100 ;
        }   
        if (relErrHigherCent > relErrLowerCent){
            histoChargedPionSpecLowPtSyst4060->SetBinError(i, histoChargedPionSpecLowPtSyst4060->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoChargedPionSpecLowPtSyst4060->SetBinError(i, histoChargedPionSpecLowPtSyst4060->GetBinContent(i)*relErrLowerCent/100);
        }         
    }   


    TH1D*   histoChargedPionMinusSpecLowPtStat6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent7_pion_minus");
    histoChargedPionMinusSpecLowPtStat6070->Sumw2();
    TH1D*   histoChargedPionMinusSpecLowPtSyst6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent7_pion_minus");
    histoChargedPionMinusSpecLowPtSyst6070->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtStat6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent7_pion_plus");
    histoChargedPionPlusSpecLowPtStat6070->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtSyst6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent7_pion_plus");
    histoChargedPionPlusSpecLowPtSyst6070->Sumw2();
    TH1D*   histoChargedPionSpecLowPtStat6070 = (TH1D*)histoChargedPionMinusSpecLowPtStat6070->Clone("histoChargedPionSpecLowPtStat6070");
    TH1D*   histoChargedPionSpecLowPtSyst6070 = (TH1D*)histoChargedPionMinusSpecLowPtSyst6070->Clone("histoChargedPionSpecLowPtSyst6070");
    histoChargedPionSpecLowPtStat6070->Add(histoChargedPionPlusSpecLowPtStat6070);
    histoChargedPionSpecLowPtSyst6070->Add(histoChargedPionPlusSpecLowPtSyst6070);
    histoChargedPionSpecLowPtStat6070->Scale(0.5);
    histoChargedPionSpecLowPtSyst6070->Scale(0.5);
    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst6070->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStat6070->SetBinContent(i, histoChargedPionSpecLowPtStat6070->GetBinContent(i)/histoChargedPionSpecLowPtStat6070->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStat6070->SetBinError(i, histoChargedPionSpecLowPtStat6070->GetBinError(i)/histoChargedPionSpecLowPtStat6070->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtSyst6070->SetBinContent(i, histoChargedPionSpecLowPtSyst6070->GetBinContent(i)/(histoChargedPionSpecLowPtSyst6070->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedPionPlusSpecLowPtSyst6070->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst6070->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst6070->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst6070->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst6070->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst6070->GetBinContent(i)*100)/2;
        }
        histoChargedPionSpecLowPtSyst6070->SetBinError(i, histoChargedPionSpecLowPtSyst6070->GetBinContent(i)*fractionalSystematicError/100.);
    }

    TH1D*   histoChargedPionMinusSpecLowPtStat7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent8_pion_minus");
    histoChargedPionMinusSpecLowPtStat7080->Sumw2();
    TH1D*   histoChargedPionMinusSpecLowPtSyst7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent8_pion_minus");
    histoChargedPionMinusSpecLowPtSyst7080->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtStat7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent8_pion_plus");
    histoChargedPionPlusSpecLowPtStat7080->Sumw2();
    TH1D*   histoChargedPionPlusSpecLowPtSyst7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent8_pion_plus");
    histoChargedPionPlusSpecLowPtSyst7080->Sumw2();
    TH1D*   histoChargedPionSpecLowPtStat7080 = (TH1D*)histoChargedPionMinusSpecLowPtStat7080->Clone("histoChargedPionSpecLowPtStat7080");
    TH1D*   histoChargedPionSpecLowPtSyst7080 = (TH1D*)histoChargedPionMinusSpecLowPtSyst7080->Clone("histoChargedPionSpecLowPtSyst7080");
    histoChargedPionSpecLowPtStat7080->Add(histoChargedPionPlusSpecLowPtStat7080);
    histoChargedPionSpecLowPtSyst7080->Add(histoChargedPionPlusSpecLowPtSyst7080);
    histoChargedPionSpecLowPtStat7080->Scale(0.5);
    histoChargedPionSpecLowPtSyst7080->Scale(0.5);
    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst7080->GetNbinsX()+1; i++){
        histoChargedPionSpecLowPtStat7080->SetBinContent(i, histoChargedPionSpecLowPtStat7080->GetBinContent(i)/histoChargedPionSpecLowPtStat7080->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtStat7080->SetBinError(i, histoChargedPionSpecLowPtStat7080->GetBinError(i)/histoChargedPionSpecLowPtStat7080->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedPionSpecLowPtSyst7080->SetBinContent(i, histoChargedPionSpecLowPtSyst7080->GetBinContent(i)/(histoChargedPionSpecLowPtSyst7080->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedPionPlusSpecLowPtSyst7080->GetBinContent(i) != 0 && histoChargedPionMinusSpecLowPtSyst7080->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedPionPlusSpecLowPtSyst7080->GetBinError(i)/histoChargedPionPlusSpecLowPtSyst7080->GetBinContent(i)*100 + histoChargedPionMinusSpecLowPtSyst7080->GetBinError(i)/histoChargedPionMinusSpecLowPtSyst7080->GetBinContent(i)*100)/2;
        }
        histoChargedPionSpecLowPtSyst7080->SetBinError(i, histoChargedPionSpecLowPtSyst7080->GetBinContent(i)*fractionalSystematicError/100.);
    }

    TH1D*   histoChargedPionSpecLowPtStat6080 = (TH1D*)histoChargedPionSpecLowPtStat6070->Clone("histoChargedPionSpecLowPtStat6080");
    TH1D*   histoChargedPionSpecLowPtSyst6080 = (TH1D*)histoChargedPionSpecLowPtSyst6070->Clone("histoChargedPionSpecLowPtSyst6080");
    histoChargedPionSpecLowPtStat6080->Add(histoChargedPionSpecLowPtStat7080);
    histoChargedPionSpecLowPtSyst6080->Add(histoChargedPionSpecLowPtSyst7080);
    histoChargedPionSpecLowPtStat6080->Scale(0.5);
    histoChargedPionSpecLowPtSyst6080->Scale(0.5);

    for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst6080->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0;
        if (histoChargedPionSpecLowPtSyst6070->GetBinContent(i) != 0){
            relErrLowerCent = histoChargedPionSpecLowPtSyst6070->GetBinError(i)/histoChargedPionSpecLowPtSyst6070->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent = 0;
        if (histoChargedPionSpecLowPtSyst7080->GetBinContent(i) != 0){
            relErrHigherCent = histoChargedPionSpecLowPtSyst7080->GetBinError(i)/histoChargedPionSpecLowPtSyst7080->GetBinContent(i)*100 ;
        }
        if (relErrHigherCent > relErrLowerCent){
            histoChargedPionSpecLowPtSyst6080->SetBinError(i, histoChargedPionSpecLowPtSyst6080->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoChargedPionSpecLowPtSyst6080->SetBinError(i, histoChargedPionSpecLowPtSyst6080->GetBinContent(i)*relErrLowerCent/100);
        }         
    }   


    // ****************************************************************************************************************
    // ******************************************charged Kaons ********************************************************
    // ****************************************************************************************************************
    
    TH1D* histoChargedKaonMinusSpecLowPtStat0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent0_kaon_minus");
    histoChargedKaonMinusSpecLowPtStat0005->Sumw2();
    TH1D* histoChargedKaonMinusSpecLowPtSyst0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent0_kaon_minus");
    histoChargedKaonMinusSpecLowPtSyst0005->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtStat0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent0_kaon_plus");
    histoChargedKaonPlusSpecLowPtStat0005->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtSyst0005 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent0_kaon_plus");
    histoChargedKaonPlusSpecLowPtSyst0005->Sumw2();
    TH1D* histoChargedKaonSpecLowPtStat0005 = (TH1D*)histoChargedKaonMinusSpecLowPtStat0005->Clone("histoChargedKaonSpecLowPtStat0005");
    TH1D* histoChargedKaonSpecLowPtSyst0005 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst0005->Clone("histoChargedKaonSpecLowPtSyst0005");
    histoChargedKaonSpecLowPtStat0005->Add(histoChargedKaonPlusSpecLowPtStat0005);
    histoChargedKaonSpecLowPtSyst0005->Add(histoChargedKaonPlusSpecLowPtSyst0005);
    histoChargedKaonSpecLowPtStat0005->Scale(0.5);
    histoChargedKaonSpecLowPtSyst0005->Scale(0.5);
        
    for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst0005->GetNbinsX()+1; i++){
        histoChargedKaonSpecLowPtStat0005->SetBinContent(i, histoChargedKaonSpecLowPtStat0005->GetBinContent(i)/histoChargedKaonSpecLowPtStat0005->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtStat0005->SetBinError(i, histoChargedKaonSpecLowPtStat0005->GetBinError(i)/histoChargedKaonSpecLowPtStat0005->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtSyst0005->SetBinContent(i, histoChargedKaonSpecLowPtSyst0005->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst0005->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedKaonMinusSpecLowPtSyst0005->GetBinContent(i) != 0 && histoChargedKaonPlusSpecLowPtSyst0005->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst0005->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst0005->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst0005->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst0005->GetBinContent(i)*100)/2;
        }
        histoChargedKaonSpecLowPtSyst0005->SetBinError(i, histoChargedKaonSpecLowPtSyst0005->GetBinContent(i)*fractionalSystematicError/100.);
    }  

    TH1D* histoChargedKaonMinusSpecLowPtStat0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent1_kaon_minus");
    histoChargedKaonMinusSpecLowPtStat0510->Sumw2();
    TH1D* histoChargedKaonMinusSpecLowPtSyst0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent1_kaon_minus");
    histoChargedKaonMinusSpecLowPtSyst0510->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtStat0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent1_kaon_plus");
    histoChargedKaonPlusSpecLowPtStat0510->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtSyst0510 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent1_kaon_plus");
    histoChargedKaonPlusSpecLowPtSyst0510->Sumw2();
    TH1D* histoChargedKaonSpecLowPtStat0510 = (TH1D*)histoChargedKaonMinusSpecLowPtStat0510->Clone("histoChargedKaonSpecLowPtStat0510");
    TH1D* histoChargedKaonSpecLowPtSyst0510 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst0510->Clone("histoChargedKaonSpecLowPtSyst0510");
    histoChargedKaonSpecLowPtStat0510->Add(histoChargedKaonPlusSpecLowPtStat0510);
    histoChargedKaonSpecLowPtSyst0510->Add(histoChargedKaonPlusSpecLowPtSyst0510);
    histoChargedKaonSpecLowPtStat0510->Scale(0.5);
    histoChargedKaonSpecLowPtSyst0510->Scale(0.5);
    
    for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst0510->GetNbinsX()+1; i++){
        histoChargedKaonSpecLowPtStat0510->SetBinContent(i, histoChargedKaonSpecLowPtStat0510->GetBinContent(i)/histoChargedKaonSpecLowPtStat0510->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtStat0510->SetBinError(i, histoChargedKaonSpecLowPtStat0510->GetBinError(i)/histoChargedKaonSpecLowPtStat0510->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtSyst0510->SetBinContent(i, histoChargedKaonSpecLowPtSyst0510->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst0510->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedKaonPlusSpecLowPtSyst0510->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst0510->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst0510->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst0510->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst0510->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst0510->GetBinContent(i)*100)/2;
        }
        histoChargedKaonSpecLowPtSyst0510->SetBinError(i, histoChargedKaonSpecLowPtSyst0510->GetBinContent(i)*fractionalSystematicError/100.);
    }
    
    TH1D* histoChargedKaonSpecLowPtStat0010 = (TH1D*)histoChargedKaonSpecLowPtStat0510->Clone("histoChargedKaonSpecLowPtStat0010");
    TH1D* histoChargedKaonSpecLowPtSyst0010 = (TH1D*)histoChargedKaonSpecLowPtSyst0510->Clone("histoChargedKaonSpecLowPtSyst0010");
    histoChargedKaonSpecLowPtStat0010->Add(histoChargedKaonSpecLowPtStat0005);
    histoChargedKaonSpecLowPtSyst0010->Add(histoChargedKaonSpecLowPtSyst0005);
    histoChargedKaonSpecLowPtStat0010->Scale(0.5);
    histoChargedKaonSpecLowPtSyst0010->Scale(0.5);

    for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst0010->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0;
        if (histoChargedKaonSpecLowPtSyst0005->GetBinContent(i) != 0){
            relErrLowerCent = histoChargedKaonSpecLowPtSyst0005->GetBinError(i)/histoChargedKaonSpecLowPtSyst0005->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent = 0;
        if (histoChargedKaonSpecLowPtSyst0510->GetBinContent(i) != 0){
            relErrHigherCent = histoChargedKaonSpecLowPtSyst0510->GetBinError(i)/histoChargedKaonSpecLowPtSyst0510->GetBinContent(i)*100 ;
        }
        if (relErrHigherCent > relErrLowerCent){
            histoChargedKaonSpecLowPtSyst0010->SetBinError(i, histoChargedKaonSpecLowPtSyst0010->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoChargedKaonSpecLowPtSyst0010->SetBinError(i, histoChargedKaonSpecLowPtSyst0010->GetBinContent(i)*relErrLowerCent/100);
        }         
    }   

    TH1D* histoChargedKaonMinusSpecLowPtStat1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent2_kaon_minus");
    histoChargedKaonMinusSpecLowPtStat1020->Sumw2();
    TH1D* histoChargedKaonMinusSpecLowPtSyst1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent2_kaon_minus");
    histoChargedKaonMinusSpecLowPtSyst1020->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtStat1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent2_kaon_plus");
    histoChargedKaonPlusSpecLowPtStat1020->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtSyst1020 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent2_kaon_plus");
    histoChargedKaonPlusSpecLowPtSyst1020->Sumw2();
    TH1D* histoChargedKaonSpecLowPtStat1020 = (TH1D*)histoChargedKaonMinusSpecLowPtStat1020->Clone("histoChargedKaonSpecLowPtStat1020");
    TH1D* histoChargedKaonSpecLowPtSyst1020 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst1020->Clone("histoChargedKaonSpecLowPtSyst1020");
    histoChargedKaonSpecLowPtStat1020->Add(histoChargedKaonPlusSpecLowPtStat1020);
    histoChargedKaonSpecLowPtSyst1020->Add(histoChargedKaonPlusSpecLowPtSyst1020);
    histoChargedKaonSpecLowPtStat1020->Scale(0.5);
    histoChargedKaonSpecLowPtSyst1020->Scale(0.5);
    for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst1020->GetNbinsX()+1; i++){
        histoChargedKaonSpecLowPtStat1020->SetBinContent(i, histoChargedKaonSpecLowPtStat1020->GetBinContent(i)/histoChargedKaonSpecLowPtStat1020->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtStat1020->SetBinError(i, histoChargedKaonSpecLowPtStat1020->GetBinError(i)/histoChargedKaonSpecLowPtStat1020->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtSyst1020->SetBinContent(i, histoChargedKaonSpecLowPtSyst1020->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst1020->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedKaonPlusSpecLowPtSyst1020->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst1020->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst1020->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst1020->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst1020->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst1020->GetBinContent(i)*100)/2;
        }
        histoChargedKaonSpecLowPtSyst1020->SetBinError(i, histoChargedKaonSpecLowPtSyst1020->GetBinContent(i)*fractionalSystematicError/100.);
    }
    
    TH1D* histoChargedKaonMinusSpecLowPtStat2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent3_kaon_minus");
    histoChargedKaonMinusSpecLowPtStat2030->Sumw2();
    TH1D* histoChargedKaonMinusSpecLowPtSyst2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent3_kaon_minus");
    histoChargedKaonMinusSpecLowPtSyst2030->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtStat2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent3_kaon_plus");
    histoChargedKaonPlusSpecLowPtStat2030->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtSyst2030 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent3_kaon_plus");
    histoChargedKaonPlusSpecLowPtSyst2030->Sumw2();
    TH1D* histoChargedKaonSpecLowPtStat2030 = (TH1D*)histoChargedKaonMinusSpecLowPtStat2030->Clone("histoChargedKaonSpecLowPtStat2030");
    TH1D* histoChargedKaonSpecLowPtSyst2030 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst2030->Clone("histoChargedKaonSpecLowPtSyst2030");
    histoChargedKaonSpecLowPtStat2030->Add(histoChargedKaonPlusSpecLowPtStat2030);
    histoChargedKaonSpecLowPtSyst2030->Add(histoChargedKaonPlusSpecLowPtSyst2030);
    histoChargedKaonSpecLowPtStat2030->Scale(0.5);
    histoChargedKaonSpecLowPtSyst2030->Scale(0.5);
    for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst2030->GetNbinsX()+1; i++){
        histoChargedKaonSpecLowPtStat2030->SetBinContent(i, histoChargedKaonSpecLowPtStat2030->GetBinContent(i)/histoChargedKaonSpecLowPtStat2030->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtStat2030->SetBinError(i, histoChargedKaonSpecLowPtStat2030->GetBinError(i)/histoChargedKaonSpecLowPtStat2030->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtSyst2030->SetBinContent(i, histoChargedKaonSpecLowPtSyst2030->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst2030->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedKaonPlusSpecLowPtSyst2030->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst2030->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst2030->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst2030->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst2030->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst2030->GetBinContent(i)*100)/2;
        }
        histoChargedKaonSpecLowPtSyst2030->SetBinError(i, histoChargedKaonSpecLowPtSyst2030->GetBinContent(i)*fractionalSystematicError/100.);
    }
    TH1D* histoChargedKaonMinusSpecLowPtStat3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent4_kaon_minus");
    histoChargedKaonMinusSpecLowPtStat3040->Sumw2();
    TH1D* histoChargedKaonMinusSpecLowPtSyst3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent4_kaon_minus");
    histoChargedKaonMinusSpecLowPtSyst3040->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtStat3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent4_kaon_plus");
    histoChargedKaonPlusSpecLowPtStat3040->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtSyst3040 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent4_kaon_plus");
    histoChargedKaonPlusSpecLowPtSyst3040->Sumw2();
    TH1D* histoChargedKaonSpecLowPtStat3040 = (TH1D*)histoChargedKaonMinusSpecLowPtStat3040->Clone("histoChargedKaonSpecLowPtStat3040");
    TH1D* histoChargedKaonSpecLowPtSyst3040 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst3040->Clone("histoChargedKaonSpecLowPtSyst3040");
    histoChargedKaonSpecLowPtStat3040->Add(histoChargedKaonPlusSpecLowPtStat3040);
    histoChargedKaonSpecLowPtSyst3040->Add(histoChargedKaonPlusSpecLowPtSyst3040);
    histoChargedKaonSpecLowPtStat3040->Scale(0.5);
    histoChargedKaonSpecLowPtSyst3040->Scale(0.5);
    for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst3040->GetNbinsX()+1; i++){
        histoChargedKaonSpecLowPtStat3040->SetBinContent(i, histoChargedKaonSpecLowPtStat3040->GetBinContent(i)/histoChargedKaonSpecLowPtStat3040->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtStat3040->SetBinError(i, histoChargedKaonSpecLowPtStat3040->GetBinError(i)/histoChargedKaonSpecLowPtStat3040->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtSyst3040->SetBinContent(i, histoChargedKaonSpecLowPtSyst3040->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst3040->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedKaonPlusSpecLowPtSyst3040->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst3040->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst3040->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst3040->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst3040->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst3040->GetBinContent(i)*100)/2;
        }
        histoChargedKaonSpecLowPtSyst3040->SetBinError(i, histoChargedKaonSpecLowPtSyst3040->GetBinContent(i)*fractionalSystematicError/100.);
    }
    TH1D* histoChargedKaonSpecLowPtStat2040 = (TH1D*)histoChargedKaonSpecLowPtStat2030->Clone("histoChargedKaonSpecLowPtStat2040");
    TH1D* histoChargedKaonSpecLowPtSyst2040 = (TH1D*)histoChargedKaonSpecLowPtSyst2030->Clone("histoChargedKaonSpecLowPtSyst2040");
    histoChargedKaonSpecLowPtStat2040->Add(histoChargedKaonSpecLowPtStat3040);
    histoChargedKaonSpecLowPtSyst2040->Add(histoChargedKaonSpecLowPtSyst3040);
    histoChargedKaonSpecLowPtStat2040->Scale(0.5);
    histoChargedKaonSpecLowPtSyst2040->Scale(0.5);

    for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst2040->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0;
        if (histoChargedKaonSpecLowPtSyst2030->GetBinContent(i) !=0){
            relErrLowerCent = histoChargedKaonSpecLowPtSyst2030->GetBinError(i)/histoChargedKaonSpecLowPtSyst2030->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent = 0;
        if (histoChargedKaonSpecLowPtSyst3040->GetBinContent(i) != 0){
            relErrHigherCent = histoChargedKaonSpecLowPtSyst3040->GetBinError(i)/histoChargedKaonSpecLowPtSyst3040->GetBinContent(i)*100 ;
        }
        if (relErrHigherCent > relErrLowerCent){
            histoChargedKaonSpecLowPtSyst2040->SetBinError(i, histoChargedKaonSpecLowPtSyst2040->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoChargedKaonSpecLowPtSyst2040->SetBinError(i, histoChargedKaonSpecLowPtSyst2040->GetBinContent(i)*relErrLowerCent/100);
        }         
    }   
    
    TH1D* histoChargedKaonMinusSpecLowPtStat4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent5_kaon_minus");
    histoChargedKaonMinusSpecLowPtStat4050->Sumw2();
    TH1D* histoChargedKaonMinusSpecLowPtSyst4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent5_kaon_minus");
    histoChargedKaonMinusSpecLowPtSyst4050->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtStat4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent5_kaon_plus");
    histoChargedKaonPlusSpecLowPtStat4050->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtSyst4050 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent5_kaon_plus");
    histoChargedKaonPlusSpecLowPtSyst4050->Sumw2();
    TH1D* histoChargedKaonSpecLowPtStat4050 = (TH1D*)histoChargedKaonMinusSpecLowPtStat4050->Clone("histoChargedKaonSpecLowPtStat4050");
    TH1D* histoChargedKaonSpecLowPtSyst4050 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst4050->Clone("histoChargedKaonSpecLowPtSyst4050");
    histoChargedKaonSpecLowPtStat4050->Add(histoChargedKaonPlusSpecLowPtStat4050);
    histoChargedKaonSpecLowPtSyst4050->Add(histoChargedKaonPlusSpecLowPtSyst4050);
    histoChargedKaonSpecLowPtStat4050->Scale(0.5);
    histoChargedKaonSpecLowPtSyst4050->Scale(0.5);
    for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst4050->GetNbinsX()+1; i++){
        histoChargedKaonSpecLowPtStat4050->SetBinContent(i, histoChargedKaonSpecLowPtStat4050->GetBinContent(i)/histoChargedKaonSpecLowPtStat4050->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtStat4050->SetBinError(i, histoChargedKaonSpecLowPtStat4050->GetBinError(i)/histoChargedKaonSpecLowPtStat4050->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtSyst4050->SetBinContent(i, histoChargedKaonSpecLowPtSyst4050->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst4050->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedKaonPlusSpecLowPtSyst4050->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst4050->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst4050->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst4050->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst4050->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst4050->GetBinContent(i)*100)/2;
        }
        histoChargedKaonSpecLowPtSyst4050->SetBinError(i, histoChargedKaonSpecLowPtSyst4050->GetBinContent(i)*fractionalSystematicError/100.);
    }

    TH1D* histoChargedKaonMinusSpecLowPtStat5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent6_kaon_minus");
    histoChargedKaonMinusSpecLowPtStat5060->Sumw2();
    TH1D* histoChargedKaonMinusSpecLowPtSyst5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent6_kaon_minus");
    histoChargedKaonMinusSpecLowPtSyst5060->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtStat5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent6_kaon_plus");
    histoChargedKaonPlusSpecLowPtStat5060->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtSyst5060 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent6_kaon_plus");
    histoChargedKaonPlusSpecLowPtSyst5060->Sumw2();
    TH1D* histoChargedKaonSpecLowPtStat5060 = (TH1D*)histoChargedKaonMinusSpecLowPtStat5060->Clone("histoChargedKaonSpecLowPtStat5060");
    TH1D* histoChargedKaonSpecLowPtSyst5060 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst5060->Clone("histoChargedKaonSpecLowPtSyst5060");
    histoChargedKaonSpecLowPtStat5060->Add(histoChargedKaonPlusSpecLowPtStat5060);
    histoChargedKaonSpecLowPtSyst5060->Add(histoChargedKaonPlusSpecLowPtSyst5060);
    histoChargedKaonSpecLowPtStat5060->Scale(0.5);
    histoChargedKaonSpecLowPtSyst5060->Scale(0.5);
    for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst5060->GetNbinsX()+1; i++){
        histoChargedKaonSpecLowPtStat5060->SetBinContent(i, histoChargedKaonSpecLowPtStat5060->GetBinContent(i)/histoChargedKaonSpecLowPtStat5060->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtStat5060->SetBinError(i, histoChargedKaonSpecLowPtStat5060->GetBinError(i)/histoChargedKaonSpecLowPtStat5060->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtSyst5060->SetBinContent(i, histoChargedKaonSpecLowPtSyst5060->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst5060->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedKaonPlusSpecLowPtSyst5060->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst5060->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst5060->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst5060->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst5060->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst5060->GetBinContent(i)*100)/2;
        }
        histoChargedKaonSpecLowPtSyst5060->SetBinError(i, histoChargedKaonSpecLowPtSyst5060->GetBinContent(i)*fractionalSystematicError/100.);
    }
    TH1D* histoChargedKaonSpecLowPtStat4060 = (TH1D*)histoChargedKaonSpecLowPtStat4050->Clone("histoChargedKaonSpecLowPtStat4060");
    TH1D* histoChargedKaonSpecLowPtSyst4060 = (TH1D*)histoChargedKaonSpecLowPtSyst4050->Clone("histoChargedKaonSpecLowPtSyst4060");
    histoChargedKaonSpecLowPtStat4060->Add(histoChargedKaonSpecLowPtStat5060);
    histoChargedKaonSpecLowPtSyst4060->Add(histoChargedKaonSpecLowPtSyst5060);
    histoChargedKaonSpecLowPtStat4060->Scale(0.5);
    histoChargedKaonSpecLowPtSyst4060->Scale(0.5);

    for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst4060->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0;
        if (histoChargedKaonSpecLowPtSyst4050->GetBinContent(i) !=0 ){
            relErrLowerCent =histoChargedKaonSpecLowPtSyst4050->GetBinError(i)/histoChargedKaonSpecLowPtSyst4050->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent =0;
        if (histoChargedKaonSpecLowPtSyst5060->GetBinContent(i) !=0 ){
            relErrHigherCent = histoChargedKaonSpecLowPtSyst5060->GetBinError(i)/histoChargedKaonSpecLowPtSyst5060->GetBinContent(i)*100 ;
        }
        if (relErrHigherCent > relErrLowerCent){
            histoChargedKaonSpecLowPtSyst4060->SetBinError(i, histoChargedKaonSpecLowPtSyst4060->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoChargedKaonSpecLowPtSyst4060->SetBinError(i, histoChargedKaonSpecLowPtSyst4060->GetBinContent(i)*relErrLowerCent/100);
        }         
    }   

    
    TH1D* histoChargedKaonMinusSpecLowPtStat6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent7_kaon_minus");
    histoChargedKaonMinusSpecLowPtStat6070->Sumw2();
    TH1D* histoChargedKaonMinusSpecLowPtSyst6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent7_kaon_minus");
    histoChargedKaonMinusSpecLowPtSyst6070->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtStat6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent7_kaon_plus");
    histoChargedKaonPlusSpecLowPtStat6070->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtSyst6070 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent7_kaon_plus");
    histoChargedKaonPlusSpecLowPtSyst6070->Sumw2();
    TH1D* histoChargedKaonSpecLowPtStat6070 = (TH1D*)histoChargedKaonMinusSpecLowPtStat6070->Clone("histoChargedKaonSpecLowPtStat6070");
    TH1D* histoChargedKaonSpecLowPtSyst6070 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst6070->Clone("histoChargedKaonSpecLowPtSyst6070");
    histoChargedKaonSpecLowPtStat6070->Add(histoChargedKaonPlusSpecLowPtStat6070);
    histoChargedKaonSpecLowPtSyst6070->Add(histoChargedKaonPlusSpecLowPtSyst6070);
    histoChargedKaonSpecLowPtStat6070->Scale(0.5);
    histoChargedKaonSpecLowPtSyst6070->Scale(0.5);
    for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst6070->GetNbinsX()+1; i++){
        histoChargedKaonSpecLowPtStat6070->SetBinContent(i, histoChargedKaonSpecLowPtStat6070->GetBinContent(i)/histoChargedKaonSpecLowPtStat6070->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtStat6070->SetBinError(i, histoChargedKaonSpecLowPtStat6070->GetBinError(i)/histoChargedKaonSpecLowPtStat6070->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtSyst6070->SetBinContent(i, histoChargedKaonSpecLowPtSyst6070->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst6070->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedKaonPlusSpecLowPtSyst6070->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst6070->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst6070->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst6070->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst6070->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst6070->GetBinContent(i)*100)/2;
        }
        histoChargedKaonSpecLowPtSyst6070->SetBinError(i, histoChargedKaonSpecLowPtSyst6070->GetBinContent(i)*fractionalSystematicError/100.);
    }

    TH1D* histoChargedKaonMinusSpecLowPtStat7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent8_kaon_minus");
    histoChargedKaonMinusSpecLowPtStat7080->Sumw2();
    TH1D* histoChargedKaonMinusSpecLowPtSyst7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent8_kaon_minus");
    histoChargedKaonMinusSpecLowPtSyst7080->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtStat7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("stat_cent8_kaon_plus");
    histoChargedKaonPlusSpecLowPtStat7080->Sumw2();
    TH1D* histoChargedKaonPlusSpecLowPtSyst7080 = (TH1D*)fileChargedIndentifiedSpectraLowPtPrelim2012->Get("sys_cent8_kaon_plus");
    histoChargedKaonPlusSpecLowPtSyst7080->Sumw2();
    TH1D* histoChargedKaonSpecLowPtStat7080 = (TH1D*)histoChargedKaonMinusSpecLowPtStat7080->Clone("histoChargedKaonSpecLowPtStat7080");
    TH1D* histoChargedKaonSpecLowPtSyst7080 = (TH1D*)histoChargedKaonMinusSpecLowPtSyst7080->Clone("histoChargedKaonSpecLowPtSyst7080");
    histoChargedKaonSpecLowPtStat7080->Add(histoChargedKaonPlusSpecLowPtStat7080);
    histoChargedKaonSpecLowPtSyst7080->Add(histoChargedKaonPlusSpecLowPtSyst7080);
    histoChargedKaonSpecLowPtStat7080->Scale(0.5);
    histoChargedKaonSpecLowPtSyst7080->Scale(0.5);
    for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst7080->GetNbinsX()+1; i++){
        histoChargedKaonSpecLowPtStat7080->SetBinContent(i, histoChargedKaonSpecLowPtStat7080->GetBinContent(i)/histoChargedKaonSpecLowPtStat7080->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtStat7080->SetBinError(i, histoChargedKaonSpecLowPtStat7080->GetBinError(i)/histoChargedKaonSpecLowPtStat7080->GetBinCenter(i)/(2*TMath::Pi()));
        histoChargedKaonSpecLowPtSyst7080->SetBinContent(i, histoChargedKaonSpecLowPtSyst7080->GetBinContent(i)/(histoChargedKaonSpecLowPtSyst7080->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoChargedKaonPlusSpecLowPtSyst7080->GetBinContent(i) != 0 && histoChargedKaonMinusSpecLowPtSyst7080->GetBinContent(i) != 0){
            fractionalSystematicError = (histoChargedKaonPlusSpecLowPtSyst7080->GetBinError(i)/histoChargedKaonPlusSpecLowPtSyst7080->GetBinContent(i)*100 + histoChargedKaonMinusSpecLowPtSyst7080->GetBinError(i)/histoChargedKaonMinusSpecLowPtSyst7080->GetBinContent(i)*100)/2;
        }
        histoChargedKaonSpecLowPtSyst7080->SetBinError(i, histoChargedKaonSpecLowPtSyst7080->GetBinContent(i)*fractionalSystematicError/100.);
    }

    TH1D* histoChargedKaonSpecLowPtStat6080 = (TH1D*)histoChargedKaonSpecLowPtStat6070->Clone("histoChargedKaonSpecLowPtStat6080");
    TH1D* histoChargedKaonSpecLowPtSyst6080 = (TH1D*)histoChargedKaonSpecLowPtSyst6070->Clone("histoChargedKaonSpecLowPtSyst6080");
    histoChargedKaonSpecLowPtStat6080->Add(histoChargedKaonSpecLowPtStat7080);
    histoChargedKaonSpecLowPtSyst6080->Add(histoChargedKaonSpecLowPtSyst7080);
    histoChargedKaonSpecLowPtStat6080->Scale(0.5);
    histoChargedKaonSpecLowPtSyst6080->Scale(0.5);

    for (Int_t i = 1; i < histoChargedKaonSpecLowPtSyst6080->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0;
        if (histoChargedKaonSpecLowPtSyst6070->GetBinContent(i) != 0){
            relErrLowerCent = histoChargedKaonSpecLowPtSyst6070->GetBinError(i)/histoChargedKaonSpecLowPtSyst6070->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent = 0;
        if (histoChargedKaonSpecLowPtSyst7080->GetBinContent(i) != 0){
            relErrHigherCent = histoChargedKaonSpecLowPtSyst7080->GetBinError(i)/histoChargedKaonSpecLowPtSyst7080->GetBinContent(i)*100 ;
        }
        if (relErrHigherCent > relErrLowerCent){
            histoChargedKaonSpecLowPtSyst6080->SetBinError(i, histoChargedKaonSpecLowPtSyst6080->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoChargedKaonSpecLowPtSyst6080->SetBinError(i, histoChargedKaonSpecLowPtSyst6080->GetBinContent(i)*relErrLowerCent/100);
        }         
    }   

    TFile* fileK0sFinalPbPb = new TFile("ExternalInputPbPb/NeutralKaon/k0s_lambda_final_spectra_12112013.root");
    TH1F* histoNeutralKaonSpecStat0005 = (TH1F*)fileK0sFinalPbPb->Get("statonly_cent0005_K0s");
    histoNeutralKaonSpecStat0005->Sumw2();
    TH1F* histoNeutralKaonSpecSyst0005 = (TH1F*)fileK0sFinalPbPb->Get("systonly_cent0005_K0s");
    histoNeutralKaonSpecSyst0005->Sumw2();
    for (Int_t i = 1; i < histoNeutralKaonSpecSyst0005->GetNbinsX()+1; i++){
        histoNeutralKaonSpecStat0005->SetBinContent(i, histoNeutralKaonSpecStat0005->GetBinContent(i)/histoNeutralKaonSpecStat0005->GetBinCenter(i)/(2*TMath::Pi()));
        histoNeutralKaonSpecStat0005->SetBinError(i, histoNeutralKaonSpecStat0005->GetBinError(i)/histoNeutralKaonSpecStat0005->GetBinCenter(i)/(2*TMath::Pi()));
        histoNeutralKaonSpecSyst0005->SetBinContent(i, histoNeutralKaonSpecSyst0005->GetBinContent(i)/(histoNeutralKaonSpecSyst0005->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoNeutralKaonSpecSyst0005->GetBinContent(i) != 0){
            fractionalSystematicError = histoNeutralKaonSpecSyst0005->GetBinError(i)/histoNeutralKaonSpecSyst0005->GetBinContent(i)*100;
        }
        histoNeutralKaonSpecSyst0005->SetBinError(i, histoNeutralKaonSpecSyst0005->GetBinContent(i)*fractionalSystematicError/100.);
    }
    
    TH1F* histoNeutralKaonSpecStat0510 = (TH1F*)fileK0sFinalPbPb->Get("statonly_cent0510_K0s");
    histoNeutralKaonSpecStat0510->Sumw2();
    TH1F* histoNeutralKaonSpecSyst0510 = (TH1F*)fileK0sFinalPbPb->Get("systonly_cent0510_K0s");
    histoNeutralKaonSpecSyst0510->Sumw2();
    for (Int_t i = 1; i < histoNeutralKaonSpecSyst0510->GetNbinsX()+1; i++){
        histoNeutralKaonSpecStat0510->SetBinContent(i, histoNeutralKaonSpecStat0510->GetBinContent(i)/histoNeutralKaonSpecStat0510->GetBinCenter(i)/(2*TMath::Pi()));
        histoNeutralKaonSpecStat0510->SetBinError(i, histoNeutralKaonSpecStat0510->GetBinError(i)/histoNeutralKaonSpecStat0510->GetBinCenter(i)/(2*TMath::Pi()));
        histoNeutralKaonSpecSyst0510->SetBinContent(i, histoNeutralKaonSpecSyst0510->GetBinContent(i)/(histoNeutralKaonSpecSyst0510->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoNeutralKaonSpecSyst0510->GetBinContent(i) != 0){
            fractionalSystematicError = histoNeutralKaonSpecSyst0510->GetBinError(i)/histoNeutralKaonSpecSyst0510->GetBinContent(i)*100;
        }
        histoNeutralKaonSpecSyst0510->SetBinError(i, histoNeutralKaonSpecSyst0510->GetBinContent(i)*fractionalSystematicError/100.);
    }

    TH1D* histoNeutralKaonSpecStat0010 = (TH1D*)histoNeutralKaonSpecStat0005->Clone("histoNeutralKaonSpecStat0010");
    TH1D* histoNeutralKaonSpecSyst0010 = (TH1D*)histoNeutralKaonSpecSyst0005->Clone("histoChargedKaonSpecLowPtSyst0010");
    histoNeutralKaonSpecStat0010->Add(histoNeutralKaonSpecStat0510);
    histoNeutralKaonSpecSyst0010->Add(histoNeutralKaonSpecSyst0510);
    histoNeutralKaonSpecStat0010->Scale(0.5);
    histoNeutralKaonSpecSyst0010->Scale(0.5);

    for (Int_t i = 1; i < histoNeutralKaonSpecSyst0010->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0;
        if (histoNeutralKaonSpecSyst0005->GetBinContent(i) != 0){
            relErrLowerCent = histoNeutralKaonSpecSyst0005->GetBinError(i)/histoNeutralKaonSpecSyst0005->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent = 0;
        if (histoNeutralKaonSpecSyst0510->GetBinContent(i) != 0){
            relErrHigherCent = histoNeutralKaonSpecSyst0510->GetBinError(i)/histoNeutralKaonSpecSyst0510->GetBinContent(i)*100 ;
        }
        if (relErrHigherCent > relErrLowerCent){
            histoNeutralKaonSpecSyst0010->SetBinError(i, histoNeutralKaonSpecSyst0010->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoNeutralKaonSpecSyst0010->SetBinError(i, histoNeutralKaonSpecSyst0010->GetBinContent(i)*relErrLowerCent/100);
        }         
    }
    
    TH1F* histoNeutralKaonSpecStat1020 = (TH1F*)fileK0sFinalPbPb->Get("statonly_cent1020_K0s");
    histoNeutralKaonSpecStat1020->Sumw2();
    TH1F* histoNeutralKaonSpecSyst1020 = (TH1F*)fileK0sFinalPbPb->Get("systonly_cent1020_K0s");
    histoNeutralKaonSpecSyst1020->Sumw2();
    for (Int_t i = 1; i < histoNeutralKaonSpecSyst1020->GetNbinsX()+1; i++){
        histoNeutralKaonSpecStat1020->SetBinContent(i, histoNeutralKaonSpecStat1020->GetBinContent(i)/histoNeutralKaonSpecStat1020->GetBinCenter(i)/(2*TMath::Pi()));
        histoNeutralKaonSpecStat1020->SetBinError(i, histoNeutralKaonSpecStat1020->GetBinError(i)/histoNeutralKaonSpecStat1020->GetBinCenter(i)/(2*TMath::Pi()));
        histoNeutralKaonSpecSyst1020->SetBinContent(i, histoNeutralKaonSpecSyst1020->GetBinContent(i)/(histoNeutralKaonSpecSyst1020->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoNeutralKaonSpecSyst1020->GetBinContent(i) != 0){
            fractionalSystematicError = histoNeutralKaonSpecSyst1020->GetBinError(i)/histoNeutralKaonSpecSyst1020->GetBinContent(i)*100;
        }
        histoNeutralKaonSpecSyst1020->SetBinError(i, histoNeutralKaonSpecSyst1020->GetBinContent(i)*fractionalSystematicError/100.);
    }

    TH1F* histoNeutralKaonSpecStat2040 = (TH1F*)fileK0sFinalPbPb->Get("statonly_cent2040_K0s");
    histoNeutralKaonSpecStat2040->Sumw2();
    TH1F* histoNeutralKaonSpecSyst2040 = (TH1F*)fileK0sFinalPbPb->Get("systonly_cent2040_K0s");
    histoNeutralKaonSpecSyst2040->Sumw2();
    for (Int_t i = 1; i < histoNeutralKaonSpecSyst2040->GetNbinsX()+1; i++){
        histoNeutralKaonSpecStat2040->SetBinContent(i, histoNeutralKaonSpecStat2040->GetBinContent(i)/histoNeutralKaonSpecStat2040->GetBinCenter(i)/(2*TMath::Pi()));
        histoNeutralKaonSpecStat2040->SetBinError(i, histoNeutralKaonSpecStat2040->GetBinError(i)/histoNeutralKaonSpecStat2040->GetBinCenter(i)/(2*TMath::Pi()));
        histoNeutralKaonSpecSyst2040->SetBinContent(i, histoNeutralKaonSpecSyst2040->GetBinContent(i)/(histoNeutralKaonSpecSyst2040->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoNeutralKaonSpecSyst2040->GetBinContent(i) != 0){
            fractionalSystematicError = histoNeutralKaonSpecSyst2040->GetBinError(i)/histoNeutralKaonSpecSyst2040->GetBinContent(i)*100;
        }
        histoNeutralKaonSpecSyst2040->SetBinError(i, histoNeutralKaonSpecSyst2040->GetBinContent(i)*fractionalSystematicError/100.);
    }

    TH1F* histoNeutralKaonSpecStat4060 = (TH1F*)fileK0sFinalPbPb->Get("statonly_cent4060_K0s");
    histoNeutralKaonSpecStat4060->Sumw2();
    TH1F* histoNeutralKaonSpecSyst4060 = (TH1F*)fileK0sFinalPbPb->Get("systonly_cent4060_K0s");
    histoNeutralKaonSpecSyst4060->Sumw2();
    for (Int_t i = 1; i < histoNeutralKaonSpecSyst4060->GetNbinsX()+1; i++){
        histoNeutralKaonSpecStat4060->SetBinContent(i, histoNeutralKaonSpecStat4060->GetBinContent(i)/histoNeutralKaonSpecStat4060->GetBinCenter(i)/(2*TMath::Pi()));
        histoNeutralKaonSpecStat4060->SetBinError(i, histoNeutralKaonSpecStat4060->GetBinError(i)/histoNeutralKaonSpecStat4060->GetBinCenter(i)/(2*TMath::Pi()));
        histoNeutralKaonSpecSyst4060->SetBinContent(i, histoNeutralKaonSpecSyst4060->GetBinContent(i)/(histoNeutralKaonSpecSyst4060->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoNeutralKaonSpecSyst4060->GetBinContent(i) != 0){
            fractionalSystematicError = histoNeutralKaonSpecSyst4060->GetBinError(i)/histoNeutralKaonSpecSyst4060->GetBinContent(i)*100;
        }
        histoNeutralKaonSpecSyst4060->SetBinError(i, histoNeutralKaonSpecSyst4060->GetBinContent(i)*fractionalSystematicError/100.);
    }

    TH1F* histoNeutralKaonSpecStat6080 = (TH1F*)fileK0sFinalPbPb->Get("statonly_cent6080_K0s");
    histoNeutralKaonSpecStat6080->Sumw2();
    TH1F* histoNeutralKaonSpecSyst6080 = (TH1F*)fileK0sFinalPbPb->Get("systonly_cent6080_K0s");
    histoNeutralKaonSpecSyst6080->Sumw2();  
    for (Int_t i = 1; i < histoNeutralKaonSpecSyst6080->GetNbinsX()+1; i++){
        histoNeutralKaonSpecStat6080->SetBinContent(i, histoNeutralKaonSpecStat6080->GetBinContent(i)/histoNeutralKaonSpecStat6080->GetBinCenter(i)/(2*TMath::Pi()));
        histoNeutralKaonSpecStat6080->SetBinError(i, histoNeutralKaonSpecStat6080->GetBinError(i)/histoNeutralKaonSpecStat6080->GetBinCenter(i)/(2*TMath::Pi()));
        histoNeutralKaonSpecSyst6080->SetBinContent(i, histoNeutralKaonSpecSyst6080->GetBinContent(i)/(histoNeutralKaonSpecSyst6080->GetBinCenter(i)*2*TMath::Pi()));
        Double_t fractionalSystematicError = 0;
        if (histoNeutralKaonSpecSyst6080->GetBinContent(i) != 0){
            fractionalSystematicError = histoNeutralKaonSpecSyst6080->GetBinError(i)/histoNeutralKaonSpecSyst6080->GetBinContent(i)*100;
        }
        histoNeutralKaonSpecSyst6080->SetBinError(i, histoNeutralKaonSpecSyst6080->GetBinContent(i)*fractionalSystematicError/100.);
    }


    
    TFile* fileChargedPionRAAFinal2013 = new TFile("ExternalInputPbPb/IdentifiedCharged/RAA_Pion_20131108.root");
    TH1D*   histoChargedPionRAAStat2040 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Stat_20_40");
    TH1D*   histoChargedPionRAASyst2040 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Syst_20_40");
    TH1D*   histoChargedPionRAAStat4060 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Stat_40_60");
    TH1D*   histoChargedPionRAASyst4060 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Syst_40_60");
    TH1D*   histoChargedPionRAAStat6080 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Stat_60_80");
    TH1D*   histoChargedPionRAASyst6080 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Syst_60_80");
    TH1D*   histoChargedPionRAAStat1020 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Stat_10_20");
    TH1D*   histoChargedPionRAASyst1020 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Syst_10_20");
    TH1D*   histoChargedPionRAAStat0005 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Stat_0_5");
    TH1D*   histoChargedPionRAASyst0005 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Syst_0_5");
    TH1D*   histoChargedPionRAAStat0510 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Stat_5_10");
    TH1D*   histoChargedPionRAASyst0510 = (TH1D*)fileChargedPionRAAFinal2013->Get("RAAPion_Syst_5_10");


    TFile* fileChargedKaonRAAFinal2013 = new TFile("ExternalInputPbPb/IdentifiedCharged/RAA_Kaon_20131108.root");
    TH1D* histoChargedKaonRAAStat2040 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Stat_20_40");
    TH1D* histoChargedKaonRAASyst2040 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Syst_20_40");
    TH1D* histoChargedKaonRAAStat4060 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Stat_40_60");
    TH1D* histoChargedKaonRAASyst4060 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Syst_40_60");
    TH1D* histoChargedKaonRAAStat6080 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Stat_60_80");
    TH1D* histoChargedKaonRAASyst6080 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Syst_60_80");
    TH1D* histoChargedKaonRAAStat1020 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Stat_10_20");
    TH1D* histoChargedKaonRAASyst1020 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Syst_10_20");
    TH1D* histoChargedKaonRAAStat0005 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Stat_0_5");
    TH1D* histoChargedKaonRAASyst0005 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Syst_0_5");
    TH1D* histoChargedKaonRAAStat0510 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Stat_5_10");
    TH1D* histoChargedKaonRAASyst0510 = (TH1D*)fileChargedKaonRAAFinal2013->Get("RAAKaon_Syst_5_10");

    
    
    cout << "bis hier" << endl;
    //******************************************************************************************************************
    //************************************ Id Charged Xiango ***********************************************************
    //******************************************************************************************************************
    TFile *fIndentifiedChargedXiango = TFile::Open("ExternalInputPbPb/IdentifiedCharged/finalspectra_Xiango_IdentifiedCharged_20130527.root");

    TList* listNegPP2760GeV = (TList*)fIndentifiedChargedXiango->Get("pp2760_NEG");
    TGraphAsymmErrors* graphXiangoChargedPionNegStatPP2760GeV = (TGraphAsymmErrors*)listNegPP2760GeV->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionNegSysPP2760GeV = (TGraphAsymmErrors*)listNegPP2760GeV->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonNegStatPP2760GeV = (TGraphAsymmErrors*)listNegPP2760GeV->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonNegSysPP2760GeV = (TGraphAsymmErrors*)listNegPP2760GeV->FindObject("grsyskaon");
    TList* listPosPP2760GeV = (TList*)fIndentifiedChargedXiango->Get("pp2760_POS");
    TGraphAsymmErrors* graphXiangoChargedPionPosStatPP2760GeV = (TGraphAsymmErrors*)listPosPP2760GeV->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionPosSysPP2760GeV = (TGraphAsymmErrors*)listPosPP2760GeV->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonPosStatPP2760GeV = (TGraphAsymmErrors*)listPosPP2760GeV->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonPosSysPP2760GeV = (TGraphAsymmErrors*)listPosPP2760GeV->FindObject("grsyskaon");   

    TList* listNegPP7TeV = (TList*)fIndentifiedChargedXiango->Get("pp7000_NEG");
    TGraphAsymmErrors* graphXiangoChargedPionNegStatPP7TeV = (TGraphAsymmErrors*)listNegPP7TeV->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionNegSysPP7TeV = (TGraphAsymmErrors*)listNegPP7TeV->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonNegStatPP7TeV = (TGraphAsymmErrors*)listNegPP7TeV->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonNegSysPP7TeV = (TGraphAsymmErrors*)listNegPP7TeV->FindObject("grsyskaon");
    TList* listPosPP7TeV = (TList*)fIndentifiedChargedXiango->Get("pp7000_POS");
    TGraphAsymmErrors* graphXiangoChargedPionPosStatPP7TeV = (TGraphAsymmErrors*)listPosPP7TeV->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionPosSysPP7TeV = (TGraphAsymmErrors*)listPosPP7TeV->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonPosStatPP7TeV = (TGraphAsymmErrors*)listPosPP7TeV->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonPosSysPP7TeV = (TGraphAsymmErrors*)listPosPP7TeV->FindObject("grsyskaon");   

    TList* listNeg0005 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_0_5_NEG");
    TGraphAsymmErrors* graphXiangoChargedPionNegStat0005 = (TGraphAsymmErrors*)listNeg0005->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionNegSys0005 = (TGraphAsymmErrors*)listNeg0005->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonNegStat0005 = (TGraphAsymmErrors*)listNeg0005->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonNegSys0005 = (TGraphAsymmErrors*)listNeg0005->FindObject("grsyskaon");
    TList* listPos0005 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_0_5_POS");
    TGraphAsymmErrors* graphXiangoChargedPionPosStat0005 = (TGraphAsymmErrors*)listPos0005->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionPosSys0005 = (TGraphAsymmErrors*)listPos0005->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonPosStat0005 = (TGraphAsymmErrors*)listPos0005->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonPosSys0005 = (TGraphAsymmErrors*)listPos0005->FindObject("grsyskaon");   
    
    TList* listNeg0510 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_5_10_NEG");
    TGraphAsymmErrors* graphXiangoChargedPionNegStat0510 = (TGraphAsymmErrors*)listNeg0510->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionNegSys0510 = (TGraphAsymmErrors*)listNeg0510->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonNegStat0510 = (TGraphAsymmErrors*)listNeg0510->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonNegSys0510 = (TGraphAsymmErrors*)listNeg0510->FindObject("grsyskaon");
    TList* listPos0510 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_5_10_POS");
    TGraphAsymmErrors* graphXiangoChargedPionPosStat0510 = (TGraphAsymmErrors*)listPos0510->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionPosSys0510 = (TGraphAsymmErrors*)listPos0510->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonPosStat0510 = (TGraphAsymmErrors*)listPos0510->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonPosSys0510 = (TGraphAsymmErrors*)listPos0510->FindObject("grsyskaon");   

    TList* listNeg1020 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_10_20_NEG");
    TGraphAsymmErrors* graphXiangoChargedPionNegStat1020 = (TGraphAsymmErrors*)listNeg1020->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionNegSys1020 = (TGraphAsymmErrors*)listNeg1020->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonNegStat1020 = (TGraphAsymmErrors*)listNeg1020->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonNegSys1020 = (TGraphAsymmErrors*)listNeg1020->FindObject("grsyskaon");
    TList* listPos1020 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_10_20_POS");
    TGraphAsymmErrors* graphXiangoChargedPionPosStat1020 = (TGraphAsymmErrors*)listPos1020->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionPosSys1020 = (TGraphAsymmErrors*)listPos1020->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonPosStat1020 = (TGraphAsymmErrors*)listPos1020->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonPosSys1020 = (TGraphAsymmErrors*)listPos1020->FindObject("grsyskaon");   

    TList* listNeg2040 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_20_40_NEG");
    TGraphAsymmErrors* graphXiangoChargedPionNegStat2040 = (TGraphAsymmErrors*)listNeg2040->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionNegSys2040 = (TGraphAsymmErrors*)listNeg2040->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonNegStat2040 = (TGraphAsymmErrors*)listNeg2040->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonNegSys2040 = (TGraphAsymmErrors*)listNeg2040->FindObject("grsyskaon");
    TList* listPos2040 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_20_40_POS");
    TGraphAsymmErrors* graphXiangoChargedPionPosStat2040 = (TGraphAsymmErrors*)listPos2040->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionPosSys2040 = (TGraphAsymmErrors*)listPos2040->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonPosStat2040 = (TGraphAsymmErrors*)listPos2040->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonPosSys2040 = (TGraphAsymmErrors*)listPos2040->FindObject("grsyskaon");   

    TList* listNeg4060 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_40_60_NEG");
    TGraphAsymmErrors* graphXiangoChargedPionNegStat4060 = (TGraphAsymmErrors*)listNeg4060->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionNegSys4060 = (TGraphAsymmErrors*)listNeg4060->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonNegStat4060 = (TGraphAsymmErrors*)listNeg4060->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonNegSys4060 = (TGraphAsymmErrors*)listNeg4060->FindObject("grsyskaon");
    TList* listPos4060 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_40_60_POS");
    TGraphAsymmErrors* graphXiangoChargedPionPosStat4060 = (TGraphAsymmErrors*)listPos4060->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionPosSys4060 = (TGraphAsymmErrors*)listPos4060->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonPosStat4060 = (TGraphAsymmErrors*)listPos4060->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonPosSys4060 = (TGraphAsymmErrors*)listPos4060->FindObject("grsyskaon");   
    
    TList* listNeg6080 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_60_80_NEG");
    TGraphAsymmErrors* graphXiangoChargedPionNegStat6080 = (TGraphAsymmErrors*)listNeg6080->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionNegSys6080 = (TGraphAsymmErrors*)listNeg6080->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonNegStat6080 = (TGraphAsymmErrors*)listNeg6080->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonNegSys6080 = (TGraphAsymmErrors*)listNeg6080->FindObject("grsyskaon");
    TList* listPos6080 = (TList*)fIndentifiedChargedXiango->Get("PbPb2760_60_80_POS");
    TGraphAsymmErrors* graphXiangoChargedPionPosStat6080 = (TGraphAsymmErrors*)listPos6080->FindObject("grstatpion");
    TGraphAsymmErrors* graphXiangoChargedPionPosSys6080 = (TGraphAsymmErrors*)listPos6080->FindObject("grsyspion");
    TGraphAsymmErrors* graphXiangoChargedKaonPosStat6080 = (TGraphAsymmErrors*)listPos6080->FindObject("grstatkaon");
    TGraphAsymmErrors* graphXiangoChargedKaonPosSys6080 = (TGraphAsymmErrors*)listPos6080->FindObject("grsyskaon");   

    
    
    //--------------------- Write Files--------------------------------------------------------------   
    cout << "bis hier" << endl;
    TFile fileChargedPions(Form("ExternalInputPbPb/IdentifiedCharged/ChargedPionSpectraPbPb_%s.root",cStamp1) ,"RECREATE");
        histoChargedPionSpecLowPtStat0005->Write("histoChargedPionSpecLowPtStat0005");
        histoChargedPionSpecLowPtSyst0005->Write("histoChargedPionSpecLowPtSyst0005");
        histoChargedPionSpecLowPtStat0510->Write("histoChargedPionSpecLowPtStat0510");
        histoChargedPionSpecLowPtSyst0510->Write("histoChargedPionSpecLowPtSyst0510");
        histoChargedPionSpecLowPtStat0010->Write("histoChargedPionSpecLowPtStat0010");
        histoChargedPionSpecLowPtSyst0010->Write("histoChargedPionSpecLowPtSyst0010");
        histoChargedPionSpecLowPtStat1020->Write("histoChargedPionSpecLowPtStat1020");
        histoChargedPionSpecLowPtSyst1020->Write("histoChargedPionSpecLowPtSyst1020");
        histoChargedPionSpecLowPtStat2030->Write("histoChargedPionSpecLowPtStat2030");
        histoChargedPionSpecLowPtSyst2030->Write("histoChargedPionSpecLowPtSyst2030");
        histoChargedPionSpecLowPtStat3040->Write("histoChargedPionSpecLowPtStat3040");
        histoChargedPionSpecLowPtSyst3040->Write("histoChargedPionSpecLowPtSyst3040");
        histoChargedPionSpecLowPtStat2040->Write("histoChargedPionSpecLowPtStat2040");
        histoChargedPionSpecLowPtSyst2040->Write("histoChargedPionSpecLowPtSyst2040");
        histoChargedPionSpecLowPtStat4050->Write("histoChargedPionSpecLowPtStat4050");
        histoChargedPionSpecLowPtSyst4050->Write("histoChargedPionSpecLowPtSyst4050");
        histoChargedPionSpecLowPtStat5060->Write("histoChargedPionSpecLowPtStat5060");
        histoChargedPionSpecLowPtSyst5060->Write("histoChargedPionSpecLowPtSyst5060");
        histoChargedPionSpecLowPtStat4060->Write("histoChargedPionSpecLowPtStat4060");
        histoChargedPionSpecLowPtSyst4060->Write("histoChargedPionSpecLowPtSyst4060");
        histoChargedPionSpecLowPtStat6070->Write("histoChargedPionSpecLowPtStat6070");
        histoChargedPionSpecLowPtSyst6070->Write("histoChargedPionSpecLowPtSyst6070");
        histoChargedPionSpecLowPtStat7080->Write("histoChargedPionSpecLowPtStat7080");
        histoChargedPionSpecLowPtSyst7080->Write("histoChargedPionSpecLowPtSyst7080");
        histoChargedPionSpecLowPtStat6080->Write("histoChargedPionSpecLowPtStat6080");
        histoChargedPionSpecLowPtSyst6080->Write("histoChargedPionSpecLowPtSyst6080");
        
        histoChargedPionSpecHighPtStat0005->Write("histoChargedPionSpecHighPtStat0005");
        histoChargedPionSpecHighPtSyst0005->Write("histoChargedPionSpecHighPtSyst0005");
        histoChargedPionSpecHighPtStat0510->Write("histoChargedPionSpecHighPtStat0510");
        histoChargedPionSpecHighPtSyst0510->Write("histoChargedPionSpecHighPtSyst0510");
        histoChargedPionSpecHighPtStat0010->Write("histoChargedPionSpecHighPtStat0010");
        histoChargedPionSpecHighPtSyst0010->Write("histoChargedPionSpecHighPtSyst0010");
        histoChargedPionSpecHighPtStat1020->Write("histoChargedPionSpecHighPtStat1020");
        histoChargedPionSpecHighPtSyst1020->Write("histoChargedPionSpecHighPtSyst1020");
        histoChargedPionSpecHighPtStat2040->Write("histoChargedPionSpecHighPtStat2040");
        histoChargedPionSpecHighPtSyst2040->Write("histoChargedPionSpecHighPtSyst2040");
        histoChargedPionSpecHighPtStat4060->Write("histoChargedPionSpecHighPtStat4060");
        histoChargedPionSpecHighPtSyst4060->Write("histoChargedPionSpecHighPtSyst4060");
        histoChargedPionSpecHighPtStat6080->Write("histoChargedPionSpecHighPtStat6080");
        histoChargedPionSpecHighPtSyst6080->Write("histoChargedPionSpecHighPtSyst6080");
            
        graphXiangoChargedPionNegStat0005->Write("graphChargedPionNegSpecXiangoStat0005");
        graphXiangoChargedPionNegSys0005->Write("graphChargedPionNegSpecXiangoSyst0005");
        graphXiangoChargedPionPosStat0005->Write("graphChargedPionPosSpecXiangoStat0005");
        graphXiangoChargedPionPosSys0005->Write("graphChargedPionPosSpecXiangoSyst0005");
        graphXiangoChargedPionNegStat0510->Write("graphChargedPionNegSpecXiangoStat0510");
        graphXiangoChargedPionNegSys0510->Write("graphChargedPionNegSpecXiangoSyst0510");
        graphXiangoChargedPionPosStat0510->Write("graphChargedPionPosSpecXiangoStat0510");
        graphXiangoChargedPionPosSys0510->Write("graphChargedPionPosSpecXiangoSyst0510");
        graphXiangoChargedPionNegStat1020->Write("graphChargedPionNegSpecXiangoStat1020");
        graphXiangoChargedPionNegSys1020->Write("graphChargedPionNegSpecXiangoSyst1020");
        graphXiangoChargedPionPosStat1020->Write("graphChargedPionPosSpecXiangoStat1020");
        graphXiangoChargedPionPosSys1020->Write("graphChargedPionPosSpecXiangoSyst1020");
        graphXiangoChargedPionNegStat2040->Write("graphChargedPionNegSpecXiangoStat2040");
        graphXiangoChargedPionNegSys2040->Write("graphChargedPionNegSpecXiangoSyst2040");
        graphXiangoChargedPionPosStat2040->Write("graphChargedPionPosSpecXiangoStat2040");
        graphXiangoChargedPionPosSys2040->Write("graphChargedPionPosSpecXiangoSyst2040");
        graphXiangoChargedPionNegStat4060->Write("graphChargedPionNegSpecXiangoStat4060");
        graphXiangoChargedPionNegSys4060->Write("graphChargedPionNegSpecXiangoSyst4060");
        graphXiangoChargedPionPosStat4060->Write("graphChargedPionPosSpecXiangoStat4060");
        graphXiangoChargedPionPosSys4060->Write("graphChargedPionPosSpecXiangoSyst4060");
        graphXiangoChargedPionNegStat6080->Write("graphChargedPionNegSpecXiangoStat6080");
        graphXiangoChargedPionNegSys6080->Write("graphChargedPionNegSpecXiangoSyst6080");
        graphXiangoChargedPionPosStat6080->Write("graphChargedPionPosSpecXiangoStat6080");
        graphXiangoChargedPionPosSys6080->Write("graphChargedPionPosSpecXiangoSyst6080");

        histoChargedPionRAAStat0005->Write("histoChargedPionRAAStat0005");
        histoChargedPionRAASyst0005->Write("histoChargedPionRAASyst0005");
        histoChargedPionRAAStat0510->Write("histoChargedPionRAAStat0510");
        histoChargedPionRAASyst0510->Write("histoChargedPionRAASyst0510");
        histoChargedPionRAAStat1020->Write("histoChargedPionRAAStat1020");
        histoChargedPionRAASyst1020->Write("histoChargedPionRAASyst1020");
        histoChargedPionRAAStat2040->Write("histoChargedPionRAAStat2040");
        histoChargedPionRAASyst2040->Write("histoChargedPionRAASyst2040");
        histoChargedPionRAAStat4060->Write("histoChargedPionRAAStat4060");
        histoChargedPionRAASyst4060->Write("histoChargedPionRAASyst4060");
        histoChargedPionRAAStat6080->Write("histoChargedPionRAAStat6080");
        histoChargedPionRAASyst6080->Write("histoChargedPionRAASyst6080");
        
    fileChargedPions.Close();

    TFile fileChargedKaons(Form("ExternalInputPbPb/IdentifiedCharged/KaonSpectraPbPb_%s.root",dateForOutput.Data()) ,"RECREATE");
        histoChargedKaonSpecLowPtStat0005->Write("histoChargedKaonSpecLowPtStat0005");
        histoChargedKaonSpecLowPtSyst0005->Write("histoChargedKaonSpecLowPtSyst0005");
        histoChargedKaonSpecLowPtStat0510->Write("histoChargedKaonSpecLowPtStat0510");
        histoChargedKaonSpecLowPtSyst0510->Write("histoChargedKaonSpecLowPtSyst0510");
        histoChargedKaonSpecLowPtStat0010->Write("histoChargedKaonSpecLowPtStat0010");
        histoChargedKaonSpecLowPtSyst0010->Write("histoChargedKaonSpecLowPtSyst0010");
        histoChargedKaonSpecLowPtStat1020->Write("histoChargedKaonSpecLowPtStat1020");
        histoChargedKaonSpecLowPtSyst1020->Write("histoChargedKaonSpecLowPtSyst1020");
        histoChargedKaonSpecLowPtStat2030->Write("histoChargedKaonSpecLowPtStat2030");
        histoChargedKaonSpecLowPtSyst2030->Write("histoChargedKaonSpecLowPtSyst2030");
        histoChargedKaonSpecLowPtStat3040->Write("histoChargedKaonSpecLowPtStat3040");
        histoChargedKaonSpecLowPtSyst3040->Write("histoChargedKaonSpecLowPtSyst3040");
        histoChargedKaonSpecLowPtStat2040->Write("histoChargedKaonSpecLowPtStat2040");
        histoChargedKaonSpecLowPtSyst2040->Write("histoChargedKaonSpecLowPtSyst2040");
        histoChargedKaonSpecLowPtStat4050->Write("histoChargedKaonSpecLowPtStat4050");
        histoChargedKaonSpecLowPtSyst4050->Write("histoChargedKaonSpecLowPtSyst4050");
        histoChargedKaonSpecLowPtStat5060->Write("histoChargedKaonSpecLowPtStat5060");
        histoChargedKaonSpecLowPtSyst5060->Write("histoChargedKaonSpecLowPtSyst5060");
        histoChargedKaonSpecLowPtStat4060->Write("histoChargedKaonSpecLowPtStat4060");
        histoChargedKaonSpecLowPtSyst4060->Write("histoChargedKaonSpecLowPtSyst4060");
        histoChargedKaonSpecLowPtStat6070->Write("histoChargedKaonSpecLowPtStat6070");
        histoChargedKaonSpecLowPtSyst6070->Write("histoChargedKaonSpecLowPtSyst6070");
        histoChargedKaonSpecLowPtStat7080->Write("histoChargedKaonSpecLowPtStat7080");
        histoChargedKaonSpecLowPtSyst7080->Write("histoChargedKaonSpecLowPtSyst7080");
        histoChargedKaonSpecLowPtStat6080->Write("histoChargedKaonSpecLowPtStat6080");
        histoChargedKaonSpecLowPtSyst6080->Write("histoChargedKaonSpecLowPtSyst6080");      
        
        graphXiangoChargedKaonNegStat0005->Write("graphChargedKaonNegSpecXiangoStat0005");
        graphXiangoChargedKaonNegSys0005->Write("graphChargedKaonNegSpecXiangoSyst0005");
        graphXiangoChargedKaonPosStat0005->Write("graphChargedKaonPosSpecXiangoStat0005");
        graphXiangoChargedKaonPosSys0005->Write("graphChargedKaonPosSpecXiangoSyst0005");
        graphXiangoChargedKaonNegStat0510->Write("graphChargedKaonNegSpecXiangoStat0510");
        graphXiangoChargedKaonNegSys0510->Write("graphChargedKaonNegSpecXiangoSyst0510");
        graphXiangoChargedKaonPosStat0510->Write("graphChargedKaonPosSpecXiangoStat0510");
        graphXiangoChargedKaonPosSys0510->Write("graphChargedKaonPosSpecXiangoSyst0510");
        graphXiangoChargedKaonNegStat1020->Write("graphChargedKaonNegSpecXiangoStat1020");
        graphXiangoChargedKaonNegSys1020->Write("graphChargedKaonNegSpecXiangoSyst1020");
        graphXiangoChargedKaonPosStat1020->Write("graphChargedKaonPosSpecXiangoStat1020");
        graphXiangoChargedKaonPosSys1020->Write("graphChargedKaonPosSpecXiangoSyst1020");
        graphXiangoChargedKaonNegStat2040->Write("graphChargedKaonNegSpecXiangoStat2040");
        graphXiangoChargedKaonNegSys2040->Write("graphChargedKaonNegSpecXiangoSyst2040");
        graphXiangoChargedKaonPosStat2040->Write("graphChargedKaonPosSpecXiangoStat2040");
        graphXiangoChargedKaonPosSys2040->Write("graphChargedKaonPosSpecXiangoSyst2040");
        graphXiangoChargedKaonNegStat4060->Write("graphChargedKaonNegSpecXiangoStat4060");
        graphXiangoChargedKaonNegSys4060->Write("graphChargedKaonNegSpecXiangoSyst4060");
        graphXiangoChargedKaonPosStat4060->Write("graphChargedKaonPosSpecXiangoStat4060");
        graphXiangoChargedKaonPosSys4060->Write("graphChargedKaonPosSpecXiangoSyst4060");
        graphXiangoChargedKaonNegStat6080->Write("graphChargedKaonNegSpecXiangoStat6080");
        graphXiangoChargedKaonNegSys6080->Write("graphChargedKaonNegSpecXiangoSyst6080");
        graphXiangoChargedKaonPosStat6080->Write("graphChargedKaonPosSpecXiangoStat6080");
        graphXiangoChargedKaonPosSys6080->Write("graphChargedKaonPosSpecXiangoSyst6080");
    
        histoChargedKaonSpecHighPtStat0005->Write("histoChargedKaonSpecHighPtStat0005");
        histoChargedKaonSpecHighPtSyst0005->Write("histoChargedKaonSpecHighPtSyst0005");
        histoChargedKaonSpecHighPtStat0510->Write("histoChargedKaonSpecHighPtStat0510");
        histoChargedKaonSpecHighPtSyst0510->Write("histoChargedKaonSpecHighPtSyst0510");
        histoChargedKaonSpecHighPtStat0010->Write("histoChargedKaonSpecHighPtStat0010");
        histoChargedKaonSpecHighPtSyst0010->Write("histoChargedKaonSpecHighPtSyst0010");
        histoChargedKaonSpecHighPtStat1020->Write("histoChargedKaonSpecHighPtStat1020");
        histoChargedKaonSpecHighPtSyst1020->Write("histoChargedKaonSpecHighPtSyst1020");
        histoChargedKaonSpecHighPtStat2040->Write("histoChargedKaonSpecHighPtStat2040");
        histoChargedKaonSpecHighPtSyst2040->Write("histoChargedKaonSpecHighPtSyst2040");
        histoChargedKaonSpecHighPtStat4060->Write("histoChargedKaonSpecHighPtStat4060");
        histoChargedKaonSpecHighPtSyst4060->Write("histoChargedKaonSpecHighPtSyst4060");
        histoChargedKaonSpecHighPtStat6080->Write("histoChargedKaonSpecHighPtStat6080");
        histoChargedKaonSpecHighPtSyst6080->Write("histoChargedKaonSpecHighPtSyst6080");

        histoNeutralKaonSpecStat0005->Write("histoNeutralKaonSpecStat0005");
        histoNeutralKaonSpecSyst0005->Write("histoNeutralKaonSpecSyst0005");
        histoNeutralKaonSpecStat0010->Write("histoNeutralKaonSpecStat0010");
        histoNeutralKaonSpecSyst0010->Write("histoNeutralKaonSpecSyst0010");
        histoNeutralKaonSpecStat0510->Write("histoNeutralKaonSpecStat0510");
        histoNeutralKaonSpecSyst0510->Write("histoNeutralKaonSpecSyst0510");
        histoNeutralKaonSpecStat1020->Write("histoNeutralKaonSpecStat1020");
        histoNeutralKaonSpecSyst1020->Write("histoNeutralKaonSpecSyst1020");
        histoNeutralKaonSpecStat2040->Write("histoNeutralKaonSpecStat2040");
        histoNeutralKaonSpecSyst2040->Write("histoNeutralKaonSpecSyst2040");
        histoNeutralKaonSpecStat4060->Write("histoNeutralKaonSpecStat4060");
        histoNeutralKaonSpecSyst4060->Write("histoNeutralKaonSpecSyst4060");
        histoNeutralKaonSpecStat6080->Write("histoNeutralKaonSpecStat6080");
        histoNeutralKaonSpecSyst6080->Write("histoNeutralKaonSpecSyst6080");
        
        histoChargedKaonRAAStat0005->Write("histoChargedKaonRAAStat0005");
        histoChargedKaonRAASyst0005->Write("histoChargedKaonRAASyst0005");
        histoChargedKaonRAAStat0510->Write("histoChargedKaonRAAStat0510");
        histoChargedKaonRAASyst0510->Write("histoChargedKaonRAASyst0510");
        histoChargedKaonRAAStat1020->Write("histoChargedKaonRAAStat1020");
        histoChargedKaonRAASyst1020->Write("histoChargedKaonRAASyst1020");
        histoChargedKaonRAAStat2040->Write("histoChargedKaonRAAStat2040");
        histoChargedKaonRAASyst2040->Write("histoChargedKaonRAASyst2040");
        histoChargedKaonRAAStat4060->Write("histoChargedKaonRAAStat4060");
        histoChargedKaonRAASyst4060->Write("histoChargedKaonRAASyst4060");
        histoChargedKaonRAAStat6080->Write("histoChargedKaonRAAStat6080");
        histoChargedKaonRAASyst6080->Write("histoChargedKaonRAASyst6080");
        
      fileChargedKaons.Close();
    }
}

