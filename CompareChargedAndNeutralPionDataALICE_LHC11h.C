/****************************************************************************************************************************
****** 		provided by Gamma Conversion Group, PWG4, 													*****
******		Ana Marin, marin@physi.uni-heidelberg.de													*****
******	   	Kathrin Koch, kkoch@physi.uni-heidelberg.de 													*****
******		Friederike Bock, friederike.bock@cern.ch													*****
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
#include "CommonHeaders/CombinationFunctions.h"
#include "CompareChargedAndNeutralPionDataALICE_LHC11h.h"

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;

//************************** Convert yield histo ********************************************************************
TH1D* ConvertYieldHisto(TH1D* input){

    // converts yield -> inv. yield

    if (!input) {
        cout << "Error: Histogram is NULL" << endl;
        return NULL;
    }

    Int_t nBins                     = input->GetNbinsX();
    Double_t newValue               = 0;
    Double_t newErrorValue          = 0;
    Double_t correctionValue        = 1;

    // divide py 2*pi
    input->Scale(1/(2*TMath::Pi()));

    // divide by pT
    for(Int_t i=0;i<nBins;i++){
        correctionValue             = 1/(input->GetBinCenter(i+1));
        input->SetBinContent(i+1,   input->GetBinContent(i+1)*correctionValue);
        input->SetBinError(i+1,     input->GetBinError(i+1)*correctionValue);
    }


    return input;
}

void CompareChargedAndNeutralPionDataALICE_LHC11h(
											TString nameFilePbPbLHC11h   = "",
                                            TString nameFilePP           = "ExternalInput/CombNeutralMesons/CombinedResultsPaper_ShiftedX_PP7and900GeV_Pub2012.root",
											TString nameFilePbPb         = "LHC11hExternalInputs/CombinedResultsPbPb_18_Feb_2014.root",
											TString suffix               = "pdf"
  										){

	gROOT->Reset();
	gROOT->SetStyle("Plain");

	StyleSettingsThesis(suffix);
	SetPlotStyle();

	TString dateForOutput 				= ReturnDateStringForOutput();
	TString outputDir 					= Form("%s/%s/ComparisonNeutralAndChargedPions",suffix.Data(),dateForOutput.Data());

    TString fileNameCharged5TeVPbPb 	= "ExternalInputPbPb/IdentifiedCharged/Spectra_PbPbLHC15o_Combined_Histograms.root";
	TString fileNameChargedPionPbPb 	= "ExternalInputPbPb/IdentifiedCharged/ChargedPionSpectraPbPb_4_Apr_2014.root";
	TString fileNameChargedPionPbPbFullPt = "ExternalInputPbPb/IdentifiedCharged/JIRA_PWGLF-258/PbPb276.fullpT.INEL.20140329.root";
	TString fileNameChargedPionPP 		= "ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_4_Apr_2014.root";

    TString fileNameEMCalPion7TeVPP 	= "ExternalInput/EMCAL/7TeV/pi0Spectrum2011EMCALAddedSignalsEffic_7TeV_150323_evi_11cd.root";
	TString fileNameEMCalPion2760GeVPP 	= "ExternalInput/EMCAL/2.76TeV/FinalCombinedXsec_pi0276TeV_25Apr2015_11a13g.root";
	TString fileNameComb2760GeVPP		= "ExternalInputPbPb/NeutralMesonspp276GeVReference/CombinedResultsPaperPP2760GeV_2017_01_26_FrediV2Clusterizer.root";

	gSystem->Exec("mkdir -p "+outputDir);

	cout << "*************************************************************************"<< endl;
	cout << "**********************  Charged Pion 5TeVPbPb full pT *******************"<< endl;
	cout << "*************************************************************************"<< endl;
	TFile* fileChargedInput5TeVPbPbFullPt   = new TFile(fileNameCharged5TeVPbPb.Data());

    TList* folderSummedPionStat             = (TList*)fileChargedInput5TeVPbPbFullPt->Get("Summed_Pion");
	TH1D* histoChargedPionSpec5TeVStat0005 	= (TH1D*)folderSummedPionStat->FindObject("hSpectraSummedPion_PbPb_Combined_0.00to5.00");
    histoChargedPionSpec5TeVStat0005 = ConvertYieldHisto(histoChargedPionSpec5TeVStat0005);
	TH1D* histoChargedPionSpec5TeVStat0510 	= (TH1D*)folderSummedPionStat->FindObject("hSpectraSummedPion_PbPb_Combined_5.00to10.00");
    histoChargedPionSpec5TeVStat0510 = ConvertYieldHisto(histoChargedPionSpec5TeVStat0510);
	TH1D* histoChargedPionSpec5TeVStat2030 	= (TH1D*)folderSummedPionStat->FindObject("hSpectraSummedPion_PbPb_Combined_20.00to30.00");
    histoChargedPionSpec5TeVStat2030 = ConvertYieldHisto(histoChargedPionSpec5TeVStat2030);
	TH1D* histoChargedPionSpec5TeVStat3040 	= (TH1D*)folderSummedPionStat->FindObject("hSpectraSummedPion_PbPb_Combined_30.00to40.00");
    histoChargedPionSpec5TeVStat3040 = ConvertYieldHisto(histoChargedPionSpec5TeVStat3040);
	TH1D* histoChargedPionSpec5TeVStat4050 	= (TH1D*)folderSummedPionStat->FindObject("hSpectraSummedPion_PbPb_Combined_40.00to50.00");
    histoChargedPionSpec5TeVStat4050 = ConvertYieldHisto(histoChargedPionSpec5TeVStat4050);

    TList* folderSummedPionSyst             = (TList*)fileChargedInput5TeVPbPbFullPt->Get("Summed_Pion_Sys");
	TH1D* histoChargedPionSpec5TeVSyst0005 	= (TH1D*)folderSummedPionSyst->FindObject("hSpectraSummedPion_PbPb_Combined_0.00to5.00");
    histoChargedPionSpec5TeVSyst0005 = ConvertYieldHisto(histoChargedPionSpec5TeVSyst0005);
	TH1D* histoChargedPionSpec5TeVSyst0510 	= (TH1D*)folderSummedPionSyst->FindObject("hSpectraSummedPion_PbPb_Combined_5.00to10.00");
    histoChargedPionSpec5TeVSyst0510 = ConvertYieldHisto(histoChargedPionSpec5TeVSyst0510);
	TH1D* histoChargedPionSpec5TeVSyst2030 	= (TH1D*)folderSummedPionSyst->FindObject("hSpectraSummedPion_PbPb_Combined_20.00to30.00");
    histoChargedPionSpec5TeVSyst2030 = ConvertYieldHisto(histoChargedPionSpec5TeVSyst2030);
	TH1D* histoChargedPionSpec5TeVSyst3040 	= (TH1D*)folderSummedPionSyst->FindObject("hSpectraSummedPion_PbPb_Combined_30.00to40.00");
    histoChargedPionSpec5TeVSyst3040 = ConvertYieldHisto(histoChargedPionSpec5TeVSyst3040);
	TH1D* histoChargedPionSpec5TeVSyst4050 	= (TH1D*)folderSummedPionSyst->FindObject("hSpectraSummedPion_PbPb_Combined_40.00to50.00");
    histoChargedPionSpec5TeVSyst4050 = ConvertYieldHisto(histoChargedPionSpec5TeVSyst4050);

    TList* folderSummedKaonStat             = (TList*)fileChargedInput5TeVPbPbFullPt->Get("Summed_Kaon");
	TH1D* histoChargedKaonSpec5TeVStat0005 	= (TH1D*)folderSummedKaonStat->FindObject("hSpectraSummedKaon_PbPb_Combined_0.00to5.00");
    histoChargedKaonSpec5TeVStat0005 = ConvertYieldHisto(histoChargedKaonSpec5TeVStat0005);
	TH1D* histoChargedKaonSpec5TeVStat0510 	= (TH1D*)folderSummedKaonStat->FindObject("hSpectraSummedKaon_PbPb_Combined_5.00to10.00");
    histoChargedKaonSpec5TeVStat0510 = ConvertYieldHisto(histoChargedKaonSpec5TeVStat0510);
	TH1D* histoChargedKaonSpec5TeVStat2030 	= (TH1D*)folderSummedKaonStat->FindObject("hSpectraSummedKaon_PbPb_Combined_20.00to30.00");
    histoChargedKaonSpec5TeVStat2030 = ConvertYieldHisto(histoChargedKaonSpec5TeVStat2030);
	TH1D* histoChargedKaonSpec5TeVStat3040 	= (TH1D*)folderSummedKaonStat->FindObject("hSpectraSummedKaon_PbPb_Combined_30.00to40.00");
    histoChargedKaonSpec5TeVStat3040 = ConvertYieldHisto(histoChargedKaonSpec5TeVStat3040);
	TH1D* histoChargedKaonSpec5TeVStat4050 	= (TH1D*)folderSummedKaonStat->FindObject("hSpectraSummedKaon_PbPb_Combined_40.00to50.00");
    histoChargedKaonSpec5TeVStat4050 = ConvertYieldHisto(histoChargedKaonSpec5TeVStat4050);

    TList* folderSummedKaonSyst             = (TList*)fileChargedInput5TeVPbPbFullPt->Get("Summed_Kaon_Sys");
	TH1D* histoChargedKaonSpec5TeVSyst0005 	= (TH1D*)folderSummedKaonSyst->FindObject("hSpectraSummedKaon_PbPb_Combined_0.00to5.00");
    histoChargedKaonSpec5TeVSyst0005 = ConvertYieldHisto(histoChargedKaonSpec5TeVSyst0005);
	TH1D* histoChargedKaonSpec5TeVSyst0510 	= (TH1D*)folderSummedKaonSyst->FindObject("hSpectraSummedKaon_PbPb_Combined_5.00to10.00");
    histoChargedKaonSpec5TeVSyst0510 = ConvertYieldHisto(histoChargedKaonSpec5TeVSyst0510);
	TH1D* histoChargedKaonSpec5TeVSyst2030 	= (TH1D*)folderSummedKaonSyst->FindObject("hSpectraSummedKaon_PbPb_Combined_20.00to30.00");
    histoChargedKaonSpec5TeVSyst2030 = ConvertYieldHisto(histoChargedKaonSpec5TeVSyst2030);
	TH1D* histoChargedKaonSpec5TeVSyst3040 	= (TH1D*)folderSummedKaonSyst->FindObject("hSpectraSummedKaon_PbPb_Combined_30.00to40.00");
    histoChargedKaonSpec5TeVSyst3040 = ConvertYieldHisto(histoChargedKaonSpec5TeVSyst3040);
	TH1D* histoChargedKaonSpec5TeVSyst4050 	= (TH1D*)folderSummedKaonSyst->FindObject("hSpectraSummedKaon_PbPb_Combined_40.00to50.00");
    histoChargedKaonSpec5TeVSyst4050 = ConvertYieldHisto(histoChargedKaonSpec5TeVSyst4050);


    TFile *fileInvYieldsChargedSpectra5TeV = new TFile("ExternalInputPbPb/IdentifiedCharged/InvYield_PbPbLHC15o_Combined_Histograms.root","RECREATE");
        histoChargedPionSpec5TeVStat0005->Write("InvYield_SummedPion_PbPb_Stat_0005");
        histoChargedPionSpec5TeVStat0510->Write("InvYield_SummedPion_PbPb_Stat_0510");
        histoChargedPionSpec5TeVStat2030->Write("InvYield_SummedPion_PbPb_Stat_2030");
        histoChargedPionSpec5TeVStat3040->Write("InvYield_SummedPion_PbPb_Stat_3040");
        histoChargedPionSpec5TeVStat4050->Write("InvYield_SummedPion_PbPb_Stat_4050");
        histoChargedPionSpec5TeVSyst0005->Write("InvYield_SummedPion_PbPb_Syst_0005");
        histoChargedPionSpec5TeVSyst0510->Write("InvYield_SummedPion_PbPb_Syst_0510");
        histoChargedPionSpec5TeVSyst2030->Write("InvYield_SummedPion_PbPb_Syst_2030");
        histoChargedPionSpec5TeVSyst3040->Write("InvYield_SummedPion_PbPb_Syst_3040");
        histoChargedPionSpec5TeVSyst4050->Write("InvYield_SummedPion_PbPb_Syst_4050");

        histoChargedKaonSpec5TeVStat0005->Write("InvYield_SummedKaon_PbPb_Stat_0005");
        histoChargedKaonSpec5TeVStat0510->Write("InvYield_SummedKaon_PbPb_Stat_0510");
        histoChargedKaonSpec5TeVStat2030->Write("InvYield_SummedKaon_PbPb_Stat_2030");
        histoChargedKaonSpec5TeVStat3040->Write("InvYield_SummedKaon_PbPb_Stat_3040");
        histoChargedKaonSpec5TeVStat4050->Write("InvYield_SummedKaon_PbPb_Stat_4050");
        histoChargedKaonSpec5TeVSyst0005->Write("InvYield_SummedKaon_PbPb_Syst_0005");
        histoChargedKaonSpec5TeVSyst0510->Write("InvYield_SummedKaon_PbPb_Syst_0510");
        histoChargedKaonSpec5TeVSyst2030->Write("InvYield_SummedKaon_PbPb_Syst_2030");
        histoChargedKaonSpec5TeVSyst3040->Write("InvYield_SummedKaon_PbPb_Syst_3040");
        histoChargedKaonSpec5TeVSyst4050->Write("InvYield_SummedKaon_PbPb_Syst_4050");
    fileInvYieldsChargedSpectra5TeV->Close();

    histoChargedPionSpec5TeVStat0005->Scale(0.5);
    histoChargedPionSpec5TeVStat0510->Scale(0.5);
    histoChargedPionSpec5TeVStat2030->Scale(0.5);
    histoChargedPionSpec5TeVStat3040->Scale(0.5);
    histoChargedPionSpec5TeVStat4050->Scale(0.5);

    histoChargedPionSpec5TeVSyst0005->Scale(0.5);
    histoChargedPionSpec5TeVSyst0510->Scale(0.5);
    histoChargedPionSpec5TeVSyst2030->Scale(0.5);
    histoChargedPionSpec5TeVSyst3040->Scale(0.5);
    histoChargedPionSpec5TeVSyst4050->Scale(0.5);

    TH1D* histoChargedPionSpec5TeVStat0010 	= (TH1D*)histoChargedPionSpec5TeVStat0005->Clone("histoChargedPionSpec5TeVStat0010");
    histoChargedPionSpec5TeVStat0010->Add(histoChargedPionSpec5TeVStat0510);
    histoChargedPionSpec5TeVStat0010->Scale(0.5);
    TGraphErrors* graphChargedPionSpec5TeVStat0010 = new TGraphErrors(histoChargedPionSpec5TeVStat0010);

    TH1D* histoChargedPionSpec5TeVSyst0010 	= (TH1D*)histoChargedPionSpec5TeVSyst0005->Clone("histoChargedPionSpec5TeVSyst0010");
    histoChargedPionSpec5TeVSyst0010->Add(histoChargedPionSpec5TeVSyst0510);
	for (Int_t i = 1; i < histoChargedPionSpec5TeVSyst0010->GetNbinsX()+1; i++){
		Double_t relErrLowerCenrHIGH = 0;
		if (histoChargedPionSpec5TeVSyst0005->GetBinContent(i) != 0){
			relErrLowerCenrHIGH= histoChargedPionSpec5TeVSyst0005->GetBinError(i)/histoChargedPionSpec5TeVSyst0005->GetBinContent(i)*100 ;
		}
		Double_t relErrLowerCentHIGH = 0;
		if (histoChargedPionSpec5TeVSyst0510->GetBinContent(i) != 0){
			relErrLowerCentHIGH = histoChargedPionSpec5TeVSyst0510->GetBinError(i)/histoChargedPionSpec5TeVSyst0510->GetBinContent(i)*100 ;
		}

		if (relErrLowerCentHIGH > relErrLowerCenrHIGH){
			histoChargedPionSpec5TeVSyst0010->SetBinError(i, histoChargedPionSpec5TeVSyst0010->GetBinContent(i)*relErrLowerCentHIGH/100);
		} else {
			histoChargedPionSpec5TeVSyst0010->SetBinError(i, histoChargedPionSpec5TeVSyst0010->GetBinContent(i)*relErrLowerCenrHIGH/100);
		}
	}
	histoChargedPionSpec5TeVSyst0010->Scale(0.5);
    TGraphErrors* graphChargedPionSpec5TeVSyst0010 = new TGraphErrors(histoChargedPionSpec5TeVSyst0010);


    TH1D* histoChargedPionSpec5TeVStat2040 	= (TH1D*)histoChargedPionSpec5TeVStat2030->Clone("histoChargedPionSpec5TeVStat2040");
    histoChargedPionSpec5TeVStat2040->Add(histoChargedPionSpec5TeVStat3040);
	histoChargedPionSpec5TeVStat2040->Scale(0.5);

    TH1D* histoChargedPionSpec5TeVStat2050 	= (TH1D*)histoChargedPionSpec5TeVStat2030->Clone("histoChargedPionSpec5TeVStat2050");
    histoChargedPionSpec5TeVStat2050->Add(histoChargedPionSpec5TeVStat3040);
    histoChargedPionSpec5TeVStat2050->Add(histoChargedPionSpec5TeVStat4050);
	histoChargedPionSpec5TeVStat2050->Scale(1./3);
    TGraphErrors* graphChargedPionSpec5TeVStat2050 = new TGraphErrors(histoChargedPionSpec5TeVStat2050);

    TH1D* histoChargedPionSpec5TeVSyst2040 	= (TH1D*)histoChargedPionSpec5TeVSyst2030->Clone("histoChargedPionSpec5TeVSyst2040");
    histoChargedPionSpec5TeVSyst2040->Add(histoChargedPionSpec5TeVSyst3040);
	for (Int_t i = 1; i < histoChargedPionSpec5TeVSyst2040->GetNbinsX()+1; i++){
		Double_t relErrLowerCenrHIGH = 0;
		if (histoChargedPionSpec5TeVSyst2030->GetBinContent(i) != 0){
			relErrLowerCenrHIGH= histoChargedPionSpec5TeVSyst2030->GetBinError(i)/histoChargedPionSpec5TeVSyst2030->GetBinContent(i)*100 ;
		}
		Double_t relErrLowerCentHIGH = 0;
		if (histoChargedPionSpec5TeVSyst3040->GetBinContent(i) != 0){
			relErrLowerCentHIGH = histoChargedPionSpec5TeVSyst3040->GetBinError(i)/histoChargedPionSpec5TeVSyst3040->GetBinContent(i)*100 ;
		}

		if (relErrLowerCentHIGH > relErrLowerCenrHIGH){
			histoChargedPionSpec5TeVSyst2040->SetBinError(i, histoChargedPionSpec5TeVSyst2040->GetBinContent(i)*relErrLowerCentHIGH/100);
		} else {
			histoChargedPionSpec5TeVSyst2040->SetBinError(i, histoChargedPionSpec5TeVSyst2040->GetBinContent(i)*relErrLowerCenrHIGH/100);
		}
	}
	histoChargedPionSpec5TeVSyst2040->Scale(0.5);

    TH1D* histoChargedPionSpec5TeVSyst2050 	= (TH1D*)histoChargedPionSpec5TeVSyst2040->Clone("histoChargedPionSpec5TeVSyst2050");
    histoChargedPionSpec5TeVSyst2050->Add(histoChargedPionSpec5TeVSyst4050);
	histoChargedPionSpec5TeVSyst2050->Scale(1./3);
	for (Int_t i = 1; i < histoChargedPionSpec5TeVSyst2050->GetNbinsX()+1; i++){
		Double_t relErrLowerCenrHIGH = 0;
		if (histoChargedPionSpec5TeVSyst2040->GetBinContent(i) != 0){
			relErrLowerCenrHIGH= histoChargedPionSpec5TeVSyst2040->GetBinError(i)/histoChargedPionSpec5TeVSyst2040->GetBinContent(i)*100 ;
		}
		Double_t relErrLowerCentHIGH = 0;
		if (histoChargedPionSpec5TeVSyst4050->GetBinContent(i) != 0){
			relErrLowerCentHIGH = histoChargedPionSpec5TeVSyst4050->GetBinError(i)/histoChargedPionSpec5TeVSyst4050->GetBinContent(i)*100 ;
		}

		if (relErrLowerCentHIGH > relErrLowerCenrHIGH){
			histoChargedPionSpec5TeVSyst2050->SetBinError(i, histoChargedPionSpec5TeVSyst2050->GetBinContent(i)*relErrLowerCentHIGH/100);
		} else {
			histoChargedPionSpec5TeVSyst2050->SetBinError(i, histoChargedPionSpec5TeVSyst2050->GetBinContent(i)*relErrLowerCenrHIGH/100);
		}
	}
    TGraphErrors* graphChargedPionSpec5TeVSyst2050 = new TGraphErrors(histoChargedPionSpec5TeVSyst2050);


    histoChargedKaonSpec5TeVStat0005->Scale(0.5);
    histoChargedKaonSpec5TeVStat0510->Scale(0.5);
    histoChargedKaonSpec5TeVStat2030->Scale(0.5);
    histoChargedKaonSpec5TeVStat3040->Scale(0.5);
    histoChargedKaonSpec5TeVStat4050->Scale(0.5);

    histoChargedKaonSpec5TeVSyst0005->Scale(0.5);
    histoChargedKaonSpec5TeVSyst0510->Scale(0.5);
    histoChargedKaonSpec5TeVSyst2030->Scale(0.5);
    histoChargedKaonSpec5TeVSyst3040->Scale(0.5);
    histoChargedKaonSpec5TeVSyst4050->Scale(0.5);

    TH1D* histoChargedKaonSpec5TeVStat0010 	= (TH1D*)histoChargedKaonSpec5TeVStat0005->Clone("histoChargedKaonSpec5TeVStat0010");
    histoChargedKaonSpec5TeVStat0010->Add(histoChargedKaonSpec5TeVStat0510);
    histoChargedKaonSpec5TeVStat0010->Scale(0.5);
    TGraphErrors* graphChargedKaonSpec5TeVStat0010 = new TGraphErrors(histoChargedKaonSpec5TeVStat0010);

    TH1D* histoChargedKaonSpec5TeVSyst0010 	= (TH1D*)histoChargedKaonSpec5TeVSyst0005->Clone("histoChargedKaonSpec5TeVSyst0010");
    histoChargedKaonSpec5TeVSyst0010->Add(histoChargedKaonSpec5TeVSyst0510);
	for (Int_t i = 1; i < histoChargedKaonSpec5TeVSyst0010->GetNbinsX()+1; i++){
		Double_t relErrLowerCenrHIGH = 0;
		if (histoChargedKaonSpec5TeVSyst0005->GetBinContent(i) != 0){
			relErrLowerCenrHIGH= histoChargedKaonSpec5TeVSyst0005->GetBinError(i)/histoChargedKaonSpec5TeVSyst0005->GetBinContent(i)*100 ;
		}
		Double_t relErrLowerCentHIGH = 0;
		if (histoChargedKaonSpec5TeVSyst0510->GetBinContent(i) != 0){
			relErrLowerCentHIGH = histoChargedKaonSpec5TeVSyst0510->GetBinError(i)/histoChargedKaonSpec5TeVSyst0510->GetBinContent(i)*100 ;
		}

		if (relErrLowerCentHIGH > relErrLowerCenrHIGH){
			histoChargedKaonSpec5TeVSyst0010->SetBinError(i, histoChargedKaonSpec5TeVSyst0010->GetBinContent(i)*relErrLowerCentHIGH/100);
		} else {
			histoChargedKaonSpec5TeVSyst0010->SetBinError(i, histoChargedKaonSpec5TeVSyst0010->GetBinContent(i)*relErrLowerCenrHIGH/100);
		}
	}
	histoChargedKaonSpec5TeVSyst0010->Scale(0.5);
    TGraphErrors* graphChargedKaonSpec5TeVSyst0010 = new TGraphErrors(histoChargedKaonSpec5TeVSyst0010);


    TH1D* histoChargedKaonSpec5TeVStat2040 	= (TH1D*)histoChargedKaonSpec5TeVStat2030->Clone("histoChargedKaonSpec5TeVStat2040");
    histoChargedKaonSpec5TeVStat2040->Add(histoChargedKaonSpec5TeVStat3040);
	histoChargedKaonSpec5TeVStat2040->Scale(0.5);

    TH1D* histoChargedKaonSpec5TeVStat2050 	= (TH1D*)histoChargedKaonSpec5TeVStat2030->Clone("histoChargedKaonSpec5TeVStat2050");
    histoChargedKaonSpec5TeVStat2050->Add(histoChargedKaonSpec5TeVStat3040);
    histoChargedKaonSpec5TeVStat2050->Add(histoChargedKaonSpec5TeVStat4050);
	histoChargedKaonSpec5TeVStat2050->Scale(1./3);
    TGraphErrors* graphChargedKaonSpec5TeVStat2050 = new TGraphErrors(histoChargedKaonSpec5TeVStat2050);

    TH1D* histoChargedKaonSpec5TeVSyst2040 	= (TH1D*)histoChargedKaonSpec5TeVSyst2030->Clone("histoChargedKaonSpec5TeVSyst2040");
    histoChargedKaonSpec5TeVSyst2040->Add(histoChargedKaonSpec5TeVSyst3040);
	for (Int_t i = 1; i < histoChargedKaonSpec5TeVSyst2040->GetNbinsX()+1; i++){
		Double_t relErrLowerCenrHIGH = 0;
		if (histoChargedKaonSpec5TeVSyst2030->GetBinContent(i) != 0){
			relErrLowerCenrHIGH= histoChargedKaonSpec5TeVSyst2030->GetBinError(i)/histoChargedKaonSpec5TeVSyst2030->GetBinContent(i)*100 ;
		}
		Double_t relErrLowerCentHIGH = 0;
		if (histoChargedKaonSpec5TeVSyst3040->GetBinContent(i) != 0){
			relErrLowerCentHIGH = histoChargedKaonSpec5TeVSyst3040->GetBinError(i)/histoChargedKaonSpec5TeVSyst3040->GetBinContent(i)*100 ;
		}

		if (relErrLowerCentHIGH > relErrLowerCenrHIGH){
			histoChargedKaonSpec5TeVSyst2040->SetBinError(i, histoChargedKaonSpec5TeVSyst2040->GetBinContent(i)*relErrLowerCentHIGH/100);
		} else {
			histoChargedKaonSpec5TeVSyst2040->SetBinError(i, histoChargedKaonSpec5TeVSyst2040->GetBinContent(i)*relErrLowerCenrHIGH/100);
		}
	}
	histoChargedKaonSpec5TeVSyst2040->Scale(0.5);

    TH1D* histoChargedKaonSpec5TeVSyst2050 	= (TH1D*)histoChargedKaonSpec5TeVSyst2040->Clone("histoChargedKaonSpec5TeVSyst2050");
    histoChargedKaonSpec5TeVSyst2050->Add(histoChargedKaonSpec5TeVSyst4050);
	histoChargedKaonSpec5TeVSyst2050->Scale(1./3);
	for (Int_t i = 1; i < histoChargedKaonSpec5TeVSyst2050->GetNbinsX()+1; i++){
		Double_t relErrLowerCenrHIGH = 0;
		if (histoChargedKaonSpec5TeVSyst2040->GetBinContent(i) != 0){
			relErrLowerCenrHIGH= histoChargedKaonSpec5TeVSyst2040->GetBinError(i)/histoChargedKaonSpec5TeVSyst2040->GetBinContent(i)*100 ;
		}
		Double_t relErrLowerCentHIGH = 0;
		if (histoChargedKaonSpec5TeVSyst4050->GetBinContent(i) != 0){
			relErrLowerCentHIGH = histoChargedKaonSpec5TeVSyst4050->GetBinError(i)/histoChargedKaonSpec5TeVSyst4050->GetBinContent(i)*100 ;
		}

		if (relErrLowerCentHIGH > relErrLowerCenrHIGH){
			histoChargedKaonSpec5TeVSyst2050->SetBinError(i, histoChargedKaonSpec5TeVSyst2050->GetBinContent(i)*relErrLowerCentHIGH/100);
		} else {
			histoChargedKaonSpec5TeVSyst2050->SetBinError(i, histoChargedKaonSpec5TeVSyst2050->GetBinContent(i)*relErrLowerCenrHIGH/100);
		}
	}
    TGraphErrors* graphChargedKaonSpec5TeVSyst2050 = new TGraphErrors(histoChargedKaonSpec5TeVSyst2050);



	cout << "*************************************************************************"<< endl;
	cout << "**************************  Charged Pion PbPb full pT *******************"<< endl;
	cout << "*************************************************************************"<< endl;

	TFile* fileChargedPionInputPbPbFullPt 			= new TFile(fileNameChargedPionPbPbFullPt.Data());
	TH1D* histoChargedPionSpecFullPtStat2040 	= (TH1D*)fileChargedPionInputPbPbFullPt->Get("hstat_PbPb276_2040_pion_sum");
	TH1D* histoChargedPionSpecFullPtSyst2040 	= (TH1D*)fileChargedPionInputPbPbFullPt->Get("hsys_PbPb276_2040_pion_sum");
	histoChargedPionSpecFullPtStat2040->Scale(0.5);
	histoChargedPionSpecFullPtSyst2040->Scale(0.5);

	TH1D* histoChargedPionSpecFullPtStat4060 	= (TH1D*)fileChargedPionInputPbPbFullPt->Get("hstat_PbPb276_4060_pion_sum");
	TH1D* histoChargedPionSpecFullPtSyst4060 	= (TH1D*)fileChargedPionInputPbPbFullPt->Get("hsys_PbPb276_4060_pion_sum");
	histoChargedPionSpecFullPtStat4060->Scale(0.5);
	histoChargedPionSpecFullPtSyst4060->Scale(0.5);
	histoChargedPionSpecFullPtStat4060->Scale(0.5);
	histoChargedPionSpecFullPtSyst4060->Scale(0.5);

	TH1D*	histoChargedPionSpecFullPtStat2050 = (TH1D*)histoChargedPionSpecFullPtStat2040->Clone("histoChargedPionSpecFullPtStat2050");
	TH1D*	histoChargedPionSpecFullPtSyst2050 = (TH1D*)histoChargedPionSpecFullPtSyst2040->Clone("histoChargedPionSpecFullPtSyst2050");
	histoChargedPionSpecFullPtStat2050->Add(histoChargedPionSpecFullPtStat4060);
	histoChargedPionSpecFullPtSyst2050->Add(histoChargedPionSpecFullPtSyst4060);
	for (Int_t i = 1; i < histoChargedPionSpecFullPtSyst2050->GetNbinsX()+1; i++){
		Double_t relErrLowerCenrHIGH = 0;
		if (histoChargedPionSpecFullPtSyst4060->GetBinContent(i) != 0){
			relErrLowerCenrHIGH= histoChargedPionSpecFullPtSyst4060->GetBinError(i)/histoChargedPionSpecFullPtSyst4060->GetBinContent(i)*100 ;
		}
		Double_t relErrLowerCentHIGH = 0;
		if (histoChargedPionSpecFullPtSyst2040->GetBinContent(i) != 0){
			relErrLowerCentHIGH = histoChargedPionSpecFullPtSyst2040->GetBinError(i)/histoChargedPionSpecFullPtSyst2040->GetBinContent(i)*100 ;
		}

		if (relErrLowerCentHIGH > relErrLowerCenrHIGH){
			histoChargedPionSpecFullPtSyst2050->SetBinError(i, histoChargedPionSpecFullPtSyst2050->GetBinContent(i)*relErrLowerCentHIGH/100);
		} else {
			histoChargedPionSpecFullPtSyst2050->SetBinError(i, histoChargedPionSpecFullPtSyst2050->GetBinContent(i)*relErrLowerCenrHIGH/100);
		}
	}
	histoChargedPionSpecFullPtStat2050->Scale(3./4);
	histoChargedPionSpecFullPtSyst2050->Scale(3./4);


	TH1D* histoChargedPionSpecOtherPtStat0005 	= (TH1D*)fileChargedPionInputPbPbFullPt->Get("hstat_PbPb276_0005_pion_sum");
	TH1D* histoChargedPionSpecOtherPtSyst0005 	= (TH1D*)fileChargedPionInputPbPbFullPt->Get("hsys_PbPb276_0005_pion_sum");
	histoChargedPionSpecOtherPtStat0005->Scale(0.5);
	histoChargedPionSpecOtherPtSyst0005->Scale(0.5);

	TH1D* histoChargedPionSpecOtherPtStat0510 	= (TH1D*)fileChargedPionInputPbPbFullPt->Get("hstat_PbPb276_0510_pion_sum");
	TH1D* histoChargedPionSpecOtherPtSyst0510 	= (TH1D*)fileChargedPionInputPbPbFullPt->Get("hsys_PbPb276_0510_pion_sum");
	histoChargedPionSpecOtherPtStat0510->Scale(0.5);
	histoChargedPionSpecOtherPtSyst0510->Scale(0.5);

	TH1D*	histoChargedPionSpecFullPtStat0010 = (TH1D*)histoChargedPionSpecOtherPtStat0510->Clone("histoChargedPionSpecFullPtStat0010");
	TH1D*	histoChargedPionSpecFullPtSyst0010 = (TH1D*)histoChargedPionSpecOtherPtSyst0510->Clone("histoChargedPionSpecFullPtSyst0010");
	histoChargedPionSpecFullPtStat0010->Add(histoChargedPionSpecOtherPtStat0005);
	histoChargedPionSpecFullPtSyst0010->Add(histoChargedPionSpecOtherPtSyst0005);
	for (Int_t i = 1; i < histoChargedPionSpecFullPtSyst0010->GetNbinsX()+1; i++){
		Double_t relErrLowerCenrHIGH = 0;
		if (histoChargedPionSpecOtherPtSyst0005->GetBinContent(i) != 0){
			relErrLowerCenrHIGH= histoChargedPionSpecOtherPtSyst0005->GetBinError(i)/histoChargedPionSpecOtherPtSyst0005->GetBinContent(i)*100 ;
		}
		Double_t relErrLowerCentHIGH = 0;
		if (histoChargedPionSpecOtherPtSyst0510->GetBinContent(i) != 0){
			relErrLowerCentHIGH = histoChargedPionSpecOtherPtSyst0510->GetBinError(i)/histoChargedPionSpecOtherPtSyst0510->GetBinContent(i)*100 ;
		}

		if (relErrLowerCentHIGH > relErrLowerCenrHIGH){
			histoChargedPionSpecFullPtSyst0010->SetBinError(i, histoChargedPionSpecFullPtSyst0010->GetBinContent(i)*relErrLowerCentHIGH/100);
		} else {
			histoChargedPionSpecFullPtSyst0010->SetBinError(i, histoChargedPionSpecFullPtSyst0010->GetBinContent(i)*relErrLowerCenrHIGH/100);
		}
	}
	histoChargedPionSpecFullPtStat0010->Scale(0.5);
	histoChargedPionSpecFullPtSyst0010->Scale(0.5);


    TH1D* histoChargedKaonSpecFullPtStat2040  = (TH1D*)fileChargedPionInputPbPbFullPt->Get("hstat_PbPb276_2040_kaon_sum");
    TH1D* histoChargedKaonSpecFullPtSyst2040  = (TH1D*)fileChargedPionInputPbPbFullPt->Get("hsys_PbPb276_2040_kaon_sum");
    histoChargedKaonSpecFullPtStat2040->Scale(0.5);
    histoChargedKaonSpecFullPtSyst2040->Scale(0.5);

	TH1D* histoChargedKaonSpecFullPtStat4060 	= (TH1D*)fileChargedPionInputPbPbFullPt->Get("hstat_PbPb276_4060_kaon_sum");
	TH1D* histoChargedKaonSpecFullPtSyst4060 	= (TH1D*)fileChargedPionInputPbPbFullPt->Get("hsys_PbPb276_4060_kaon_sum");
	histoChargedKaonSpecFullPtStat4060->Scale(0.5);
	histoChargedKaonSpecFullPtSyst4060->Scale(0.5);

	histoChargedKaonSpecFullPtStat4060->Scale(0.5);
	histoChargedKaonSpecFullPtSyst4060->Scale(0.5);
	TH1D*	histoChargedKaonSpecFullPtStat2050 = (TH1D*)histoChargedKaonSpecFullPtStat2040->Clone("histoChargedKaonSpecFullPtStat2050");
	TH1D*	histoChargedKaonSpecFullPtSyst2050 = (TH1D*)histoChargedKaonSpecFullPtSyst2040->Clone("histoChargedKaonSpecFullPtSyst2050");
	histoChargedKaonSpecFullPtStat2050->Add(histoChargedKaonSpecFullPtStat4060);
	histoChargedKaonSpecFullPtSyst2050->Add(histoChargedKaonSpecFullPtSyst4060);
	for (Int_t i = 1; i < histoChargedKaonSpecFullPtSyst2050->GetNbinsX()+1; i++){
		Double_t relErrLowerCenrHIGH = 0;
		if (histoChargedKaonSpecFullPtSyst4060->GetBinContent(i) != 0){
			relErrLowerCenrHIGH= histoChargedKaonSpecFullPtSyst4060->GetBinError(i)/histoChargedKaonSpecFullPtSyst4060->GetBinContent(i)*100 ;
		}
		Double_t relErrLowerCentHIGH = 0;
		if (histoChargedKaonSpecFullPtSyst2040->GetBinContent(i) != 0){
			relErrLowerCentHIGH = histoChargedKaonSpecFullPtSyst2040->GetBinError(i)/histoChargedKaonSpecFullPtSyst2040->GetBinContent(i)*100 ;
		}

		if (relErrLowerCentHIGH > relErrLowerCenrHIGH){
			histoChargedKaonSpecFullPtSyst2050->SetBinError(i, histoChargedKaonSpecFullPtSyst2050->GetBinContent(i)*relErrLowerCentHIGH/100);
		} else {
			histoChargedKaonSpecFullPtSyst2050->SetBinError(i, histoChargedKaonSpecFullPtSyst2050->GetBinContent(i)*relErrLowerCenrHIGH/100);
		}
	}
	histoChargedKaonSpecFullPtStat2050->Scale(3./4);
	histoChargedKaonSpecFullPtSyst2050->Scale(3./4);

    TH1D* histoChargedKaonSpecOtherPtStat0005   = (TH1D*)fileChargedPionInputPbPbFullPt->Get("hstat_PbPb276_0005_kaon_sum");
    TH1D* histoChargedKaonSpecOtherPtSyst0005   = (TH1D*)fileChargedPionInputPbPbFullPt->Get("hsys_PbPb276_0005_kaon_sum");
    histoChargedKaonSpecOtherPtStat0005->Scale(0.5);
    histoChargedKaonSpecOtherPtSyst0005->Scale(0.5);

    TH1D* histoChargedKaonSpecOtherPtStat0510   = (TH1D*)fileChargedPionInputPbPbFullPt->Get("hstat_PbPb276_0510_kaon_sum");
    TH1D* histoChargedKaonSpecOtherPtSyst0510   = (TH1D*)fileChargedPionInputPbPbFullPt->Get("hsys_PbPb276_0510_kaon_sum");
    histoChargedKaonSpecOtherPtStat0510->Scale(0.5);
    histoChargedKaonSpecOtherPtSyst0510->Scale(0.5);

    TH1D*   histoChargedKaonSpecFullPtStat0010 = (TH1D*)histoChargedKaonSpecOtherPtStat0510->Clone("histoChargedKaonSpecFullPtStat0010");
    TH1D*   histoChargedKaonSpecFullPtSyst0010 = (TH1D*)histoChargedKaonSpecOtherPtSyst0510->Clone("histoChargedKaonSpecFullPtSyst0010");
    histoChargedKaonSpecFullPtStat0010->Add(histoChargedKaonSpecOtherPtStat0005);
    histoChargedKaonSpecFullPtSyst0010->Add(histoChargedKaonSpecOtherPtSyst0005);

    for (Int_t i = 1; i < histoChargedKaonSpecFullPtSyst0010->GetNbinsX()+1; i++){
        Double_t relErrLowerCenrHIGH = 0;
        if (histoChargedKaonSpecOtherPtSyst0005->GetBinContent(i) != 0){
            relErrLowerCenrHIGH= histoChargedKaonSpecOtherPtSyst0005->GetBinError(i)/histoChargedKaonSpecOtherPtSyst0005->GetBinContent(i)*100 ;
        }
        Double_t relErrLowerCentHIGH = 0;
        if (histoChargedKaonSpecOtherPtSyst0510->GetBinContent(i) != 0){
            relErrLowerCentHIGH = histoChargedKaonSpecOtherPtSyst0510->GetBinError(i)/histoChargedKaonSpecOtherPtSyst0510->GetBinContent(i)*100 ;
        }

        if (relErrLowerCentHIGH > relErrLowerCenrHIGH){
            histoChargedKaonSpecFullPtSyst0010->SetBinError(i, histoChargedKaonSpecFullPtSyst0010->GetBinContent(i)*relErrLowerCentHIGH/100);
        } else {
            histoChargedKaonSpecFullPtSyst0010->SetBinError(i, histoChargedKaonSpecFullPtSyst0010->GetBinContent(i)*relErrLowerCenrHIGH/100);
        }
    }
    histoChargedKaonSpecFullPtStat0010->Scale(0.5);
    histoChargedKaonSpecFullPtSyst0010->Scale(0.5);

	cout << "*************************************************************************"<< endl;
	cout << "*****************************  charged pion PbPb ************************"<< endl;
	cout << "*************************************************************************"<< endl;
	TFile* fileChargedPionInputPbPb 			= new TFile(fileNameChargedPionPbPb.Data());
	TH1D* histoChargedPionSpecHighPtStat0005 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtStat0005");
	TH1D* histoChargedPionSpecHighPtSyst0005 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtSyst0005");
	TH1D* histoChargedPionSpecHighPtStat0510 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtStat0510");
	TH1D* histoChargedPionSpecHighPtSyst0510 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtSyst0510");
	TH1D* histoChargedPionSpecHighPtStat1020 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtStat1020");
	TH1D* histoChargedPionSpecHighPtSyst1020 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtSyst1020");
	TH1D* histoChargedPionSpecHighPtSyst2040 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtSyst2040");
	TH1D* histoChargedPionSpecHighPtStat2040 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtStat2040");
	TH1D* histoChargedPionSpecHighPtSyst4060 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtSyst4060");
	TH1D* histoChargedPionSpecHighPtStat4060 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtStat4060");
	TH1D* histoChargedPionSpecHighPtSyst6080 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtSyst6080");
	TH1D* histoChargedPionSpecHighPtStat6080 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecHighPtStat6080");
	TH1D* histoChargedPionSpecLowPtStat0005 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtStat0005");
	TH1D* histoChargedPionSpecLowPtSyst0005 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtSyst0005");
	TH1D* histoChargedPionSpecLowPtStat0510 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtStat0510");
	TH1D* histoChargedPionSpecLowPtSyst0510 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtSyst0510");
	TH1D* histoChargedPionSpecLowPtStat1020 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtStat1020");
	TH1D* histoChargedPionSpecLowPtSyst1020 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtSyst1020");
	TH1D* histoChargedPionSpecLowPtStat2040 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtStat2040");
	TH1D* histoChargedPionSpecLowPtSyst2040 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtSyst2040");
	TH1D* histoChargedPionSpecLowPtStat4060 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtStat4060");
	TH1D* histoChargedPionSpecLowPtSyst4060 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtSyst4060");
	TH1D* histoChargedPionSpecLowPtStat6080 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtStat6080");
	TH1D* histoChargedPionSpecLowPtSyst6080 	= (TH1D*)fileChargedPionInputPbPb->Get("histoChargedPionSpecLowPtSyst6080");


	TH1D*	histoChargedPionSpecLowPtStat0010 = (TH1D*)histoChargedPionSpecLowPtStat0510->Clone("histoChargedPionSpecLowPtStat0010");
	TH1D*	histoChargedPionSpecLowPtSyst0010 = (TH1D*)histoChargedPionSpecLowPtSyst0510->Clone("histoChargedPionSpecLowPtSyst0010");
	histoChargedPionSpecLowPtStat0010->Add(histoChargedPionSpecLowPtStat0005);
	histoChargedPionSpecLowPtSyst0010->Add(histoChargedPionSpecLowPtSyst0005);
	histoChargedPionSpecLowPtStat0010->Scale(0.5);
	histoChargedPionSpecLowPtSyst0010->Scale(0.5);
	for (Int_t i = 1; i < histoChargedPionSpecLowPtSyst0010->GetNbinsX()+1; i++){
		Double_t relErrLowerCenrLOW = 0;
		if (histoChargedPionSpecLowPtSyst0005->GetBinContent(i) != 0){
			relErrLowerCenrLOW= histoChargedPionSpecLowPtSyst0005->GetBinError(i)/histoChargedPionSpecLowPtSyst0005->GetBinContent(i)*100 ;
		}
		Double_t relErrLowerCentLOW = 0;
		if (histoChargedPionSpecLowPtSyst0510->GetBinContent(i) != 0){
			relErrLowerCentLOW = histoChargedPionSpecLowPtSyst0510->GetBinError(i)/histoChargedPionSpecLowPtSyst0510->GetBinContent(i)*100 ;
		}

		if (relErrLowerCentLOW > relErrLowerCenrLOW){
			histoChargedPionSpecLowPtSyst0010->SetBinError(i, histoChargedPionSpecLowPtSyst0010->GetBinContent(i)*relErrLowerCentLOW/100);
		} else {
			histoChargedPionSpecLowPtSyst0010->SetBinError(i, histoChargedPionSpecLowPtSyst0010->GetBinContent(i)*relErrLowerCenrLOW/100);
		}
	}

	TH1D*	histoChargedPionSpecHighPtStat0010 = (TH1D*)histoChargedPionSpecHighPtStat0510->Clone("histoChargedPionSpecHighPtStat0010");
	TH1D*	histoChargedPionSpecHighPtSyst0010 = (TH1D*)histoChargedPionSpecHighPtSyst0510->Clone("histoChargedPionSpecHighPtSyst0010");
	histoChargedPionSpecHighPtStat0010->Add(histoChargedPionSpecHighPtStat0005);
	histoChargedPionSpecHighPtSyst0010->Add(histoChargedPionSpecHighPtSyst0005);
	histoChargedPionSpecHighPtStat0010->Scale(0.5);
	histoChargedPionSpecHighPtSyst0010->Scale(0.5);
	for (Int_t i = 1; i < histoChargedPionSpecHighPtSyst0010->GetNbinsX()+1; i++){
		Double_t relErrLowerCenrHIGH = 0;
		if (histoChargedPionSpecHighPtSyst0005->GetBinContent(i) != 0){
			relErrLowerCenrHIGH= histoChargedPionSpecHighPtSyst0005->GetBinError(i)/histoChargedPionSpecHighPtSyst0005->GetBinContent(i)*100 ;
		}
		Double_t relErrLowerCentHIGH = 0;
		if (histoChargedPionSpecHighPtSyst0510->GetBinContent(i) != 0){
			relErrLowerCentHIGH = histoChargedPionSpecHighPtSyst0510->GetBinError(i)/histoChargedPionSpecHighPtSyst0510->GetBinContent(i)*100 ;
		}

		if (relErrLowerCentHIGH > relErrLowerCenrHIGH){
			histoChargedPionSpecHighPtSyst0010->SetBinError(i, histoChargedPionSpecHighPtSyst0010->GetBinContent(i)*relErrLowerCentHIGH/100);
		} else {
			histoChargedPionSpecHighPtSyst0010->SetBinError(i, histoChargedPionSpecHighPtSyst0010->GetBinContent(i)*relErrLowerCenrHIGH/100);
		}
	}

	cout << "*************************************************************************"<< endl;
	cout << "***************************  Pi0 PbPb 2010 ******************************"<< endl;
	cout << "*************************************************************************"<< endl;

	TFile* fCombResults= new TFile(nameFilePbPb.Data());
	TGraphAsymmErrors*	graphYieldCombStatPi02760GeV	 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPComb2760GeV_StatErr");
	TGraphAsymmErrors*	graphYieldCombSysPi02760GeV 		= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPComb2760GeV_SysErr");
	TGraphAsymmErrors*	graphYieldPCMStatPi02760GeV 		= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPPCM2760GeV_StatErr");
	TGraphAsymmErrors*	graphYieldPCMSysPi02760GeV 			= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPPCM2760GeV_SysErr");
	TGraphAsymmErrors*	graphYieldPHOSStatPi02760GeV 		= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPPHOS2760GeV_StatErr");
	TGraphAsymmErrors*	graphYieldPHOSSysPi02760GeV 		= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPPPHOS2760GeV_SysErr");
	TGraphAsymmErrors*	graphYieldPi0CombPbPb0005StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_0005");
	TGraphAsymmErrors*	graphYieldPi0CombPbPb0005SysErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_0005");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb0005StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMStatErr_0005");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb0005SysErr 		= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMSysErr_0005");
	TGraphAsymmErrors*	graphYieldPi0PHOSPbPb0005StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSStatErr_0005");
	TGraphAsymmErrors*	graphYieldPi0PHOSPbPb0005SysErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSSysErr_0005");
	TGraphAsymmErrors*	graphYieldPi0CombPbPb0510StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_0510");
	TGraphAsymmErrors*	graphYieldPi0CombPbPb0510SysErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_0510");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb0510StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMStatErr_0510");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb0510SysErr 		= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMSysErr_0510");
	TGraphAsymmErrors*	graphYieldPi0PHOSPbPb0510StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSStatErr_0510");
	TGraphAsymmErrors*	graphYieldPi0PHOSPbPb0510SysErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSSysErr_0510");

	TGraphAsymmErrors*	graphYieldPi0CombPbPb0010StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_0010");
	TGraphAsymmErrors*	graphYieldPi0CombPbPb0010SysErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_0010");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb0010StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMStatErr_0010");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb0010SysErr 		= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMSysErr_0010");
	TGraphAsymmErrors*	graphYieldPi0PHOSPbPb0010StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSStatErr_0010");
	TGraphAsymmErrors*	graphYieldPi0PHOSPbPb0010SysErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSSysErr_0010");

	TGraphAsymmErrors*	graphYieldPi0CombPbPb1020StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_1020");
	TGraphAsymmErrors*	graphYieldPi0CombPbPb1020SysErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_1020");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb1020StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMStatErr_1020");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb1020SysErr 		= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMSysErr_1020");
	TGraphAsymmErrors*	graphYieldPi0PHOSPbPb1020StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSStatErr_1020");
	TGraphAsymmErrors*	graphYieldPi0PHOSPbPb1020SysErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSSysErr_1020");
	TGraphAsymmErrors*	graphYieldPi0CombPbPb2040StatErr	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_2040");
	TGraphAsymmErrors*	graphYieldPi0CombPbPb2040SysErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_2040");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb2040StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMStatErr_2040");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb2040SysErr 		= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMSysErr_2040");
	TGraphAsymmErrors*	graphYieldPi0PHOSPbPb2040StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSStatErr_2040");
	TGraphAsymmErrors*	graphYieldPi0PHOSPbPb2040SysErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSSysErr_2040");
	TGraphAsymmErrors*	graphYieldPi0CombPbPb4060StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_4060");
	TGraphAsymmErrors*	graphYieldPi0CombPbPb4060SysErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_4060");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb4060StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMStatErr_4060");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb4060SysErr 		= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMSysErr_4060");
	TGraphAsymmErrors*	graphYieldPi0PHOSPbPb4060StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSStatErr_4060");
	TGraphAsymmErrors*	graphYieldPi0PHOSPbPb4060SysErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSSysErr_4060");
	TGraphAsymmErrors*	graphYieldPi0CombPbPb6080StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbStatErr_6080");
	TGraphAsymmErrors*	graphYieldPi0CombPbPb6080SysErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbSysErr_6080");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb6080StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMStatErr_6080");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb6080SysErr 		= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPCMSysErr_6080");
	TGraphAsymmErrors*	graphYieldPi0PHOSPbPb6080StatErr 	= (TGraphAsymmErrors*)fCombResults->Get("InvYieldPbPbPHOSStatErr_6080");

	cout << "*************************************************************************"<< endl;
	cout << "******************************  PbPb 0-10% *******************************"<< endl;
	cout << "*************************************************************************"<< endl;

	TGraphAsymmErrors* graphYieldPi0CombPbPb0010StatErrCopy 	= (TGraphAsymmErrors*) graphYieldPi0CombPbPb0010StatErr->Clone("graphYieldPi0CombPbPb0010StatErrCopy");
	TGraphAsymmErrors* graphYieldPi0CombPbPb0010SysErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0CombPbPb0010SysErr->Clone("graphYieldPi0CombPbPb0010SysErrCopy");
	TGraphAsymmErrors* graphYieldPi0PCMPbPb0010StatErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0PCMPbPb0010StatErr->Clone("graphYieldPi0PCMPbPb0010StatErrCopy");
	TGraphAsymmErrors* graphYieldPi0PCMPbPb0010SysErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0PCMPbPb0010SysErr->Clone("graphYieldPi0PCMPbPb0010SysErrCopy");
	TGraphAsymmErrors* graphYieldPi0PHOSPbPb0010StatErrCopy 	= (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb0010StatErr->Clone("graphYieldPi0PHOSPbPb0010StatErrCopy");
	TGraphAsymmErrors* graphYieldPi0PHOSPbPb0010SysErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb0010SysErr->Clone("graphYieldPi0PHOSPbPb0010SysErrCopy");

	TGraphErrors* graphRatioHighPtChargedPionsComb0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb0010StatErrCopy, graphYieldPi0CombPbPb0010SysErrCopy,
																										  histoChargedPionSpecHighPtStat0010, histoChargedPionSpecHighPtSyst0010,
																									      kTRUE,  kTRUE,
																									      &graphYieldCombStatPi00010RebinnedHighPtComb, &graphYieldCombSysPi00010RebinnedHighPtComb,
																									      &graphChargedPionSpecHighPtStat0010HighPtComb, &graphChargedPionSpecHighPtSyst0010HighPtComb);
// 	graphRatioHighPtChargedPionsComb0010->Print();

	TGraphErrors* graphRatioLowPtChargedPionsComb0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb0010StatErrCopy, graphYieldPi0CombPbPb0010SysErrCopy,
																										 histoChargedPionSpecLowPtStat0010, histoChargedPionSpecLowPtSyst0010,
																									     kTRUE,  kTRUE,
																									     &graphYieldCombStatPi00010RebinnedLowPtComb, &graphYieldCombSysPi00010RebinnedLowPtComb,
																									     &graphChargedPionSpecLowPtStat0010LowPtComb, &graphChargedPionSpecLowPtSyst0010LowPtComb);
// 	graphRatioLowPtChargedPionsComb0010->Print();

	TGraphErrors* graphRatioHighPtChargedPionsPCM0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb0010StatErrCopy, graphYieldPi0PCMPbPb0010SysErrCopy,
																										 histoChargedPionSpecHighPtStat0010, histoChargedPionSpecHighPtSyst0010,
																										 kTRUE,  kTRUE,
																									     &graphYieldPCMStatPi00010RebinnedHighPtPCM, &graphYieldPCMSysPi00010RebinnedHighPtPCM,
																									     &graphChargedPionSpecHighPtStat0010HighPtPCM, &graphChargedPionSpecHighPtSyst0010HighPtPCM);
// 	graphRatioHighPtChargedPionsPCM0010->Print();

	TGraphErrors* graphRatioLowPtChargedPionsPCM0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb0010StatErrCopy, graphYieldPi0PCMPbPb0010SysErrCopy,
																										histoChargedPionSpecLowPtStat0010, histoChargedPionSpecLowPtSyst0010,
																										kTRUE,  kTRUE,
																										&graphYieldPCMStatPi00010RebinnedLowPtPCM, &graphYieldPCMSysPi00010RebinnedLowPtPCM,
																										&graphChargedPionSpecLowPtStat0010LowPtPCM, &graphChargedPionSpecLowPtSyst0010LowPtPCM);
// 	graphRatioLowPtChargedPionsPCM0010->Print();

	graphYieldPi0PHOSPbPb0010StatErrCopy->RemovePoint(0);
	graphYieldPi0PHOSPbPb0010SysErrCopy->RemovePoint(0);
	TGraphErrors* graphRatioHighPtChargedPionsPHOS0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb0010StatErrCopy, graphYieldPi0PHOSPbPb0010SysErrCopy,
																										  histoChargedPionSpecHighPtStat0010, histoChargedPionSpecHighPtSyst0010,
																										  kTRUE,  kTRUE,
																										  &graphYieldPHOSStatPi00010RebinnedHighPtPHOS, &graphYieldPHOSSysPi00010RebinnedHighPtPHOS,
																									      &graphChargedPionSpecHighPtStat0010HighPtPHOS, &graphChargedPionSpecHighPtSyst0010HighPtPHOS);
// 	graphRatioHighPtChargedPionsPHOS0010->Print();

	TGraphErrors* graphRatioLowPtChargedPionsPHOS0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb0010StatErrCopy, graphYieldPi0PHOSPbPb0010SysErrCopy,
																										 histoChargedPionSpecLowPtStat0010, histoChargedPionSpecLowPtSyst0010,
																										 kTRUE,  kTRUE,
																										 &graphYieldPHOSStatPi00010RebinnedLowPtPHOS, &graphYieldPHOSSysPi00010RebinnedLowPtPHOS,
																										 &graphChargedPionSpecLowPtStat0010LowPtPHOS, &graphChargedPionSpecLowPtSyst0010LowPtPHOS);
// 	graphRatioLowPtChargedPionsPHOS0010->Print();

	TGraphErrors* graphRatioFullPtChargedPionsComb0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb0010StatErrCopy, graphYieldPi0CombPbPb0010SysErrCopy,
																										  histoChargedPionSpecFullPtStat0010, histoChargedPionSpecFullPtSyst0010,
																									      kTRUE,  kTRUE,
																									      &graphYieldCombStatPi00010FullPtPtComb, &graphYieldCombSysPi00010FullPtPtComb,
																									      &graphChargedPionSpecFullPtStat0010FullPtComb, &graphChargedPionSpecFullPtSyst0010FullPtComb);
// 	graphRatioFullPtChargedPionsComb0010->Print();

	TGraphErrors* graphRatioFullPtChargedPionsPCM0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb0010StatErrCopy, graphYieldPi0PCMPbPb0010SysErrCopy,
																										 histoChargedPionSpecFullPtStat0010, histoChargedPionSpecFullPtSyst0010,
																										 kTRUE,  kTRUE,
																									     &graphYieldPCMStatPi00010FullPtPtPCM, &graphYieldPCMSysPi00010FullPtPtPCM,
																									     &graphChargedPionSpecFullPtStat0010FullPtPCM, &graphChargedPionSpecFullPtSyst0010FullPtPCM);
// 	graphRatioFullPtChargedPionsPCM0010->Print();

	graphYieldPi0PHOSPbPb0010StatErrCopy->RemovePoint(0);
	graphYieldPi0PHOSPbPb0010SysErrCopy->RemovePoint(0);
	TGraphErrors* graphRatioFullPtChargedPionsPHOS0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb0010StatErrCopy, graphYieldPi0PHOSPbPb0010SysErrCopy,
																										  histoChargedPionSpecFullPtStat0010, histoChargedPionSpecFullPtSyst0010,
																										  kTRUE,  kTRUE,
																										  &graphYieldPHOSStatPi00010FullPtPtPHOS, &graphYieldPHOSSysPi00010FullPtPtPHOS,
																									      &graphChargedPionSpecFullPtStat0010FullPtPHOS, &graphChargedPionSpecFullPtSyst0010FullPtPHOS);
// 	graphRatioFullPtChargedPionsPHOS0010->Print();


	cout << "*************************************************************************"<< endl;
	cout << "******************************  PbPb 20-40% *****************************"<< endl;
	cout << "*************************************************************************"<< endl;

	TGraphAsymmErrors* graphYieldPi0CombPbPb2040StatErrCopy 	= (TGraphAsymmErrors*) graphYieldPi0CombPbPb2040StatErr->Clone("graphYieldPi0CombPbPb2040StatErrCopy");
	TGraphAsymmErrors* graphYieldPi0CombPbPb2040SysErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0CombPbPb2040SysErr->Clone("graphYieldPi0CombPbPb2040SysErrCopy");
	TGraphAsymmErrors* graphYieldPi0PCMPbPb2040StatErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0PCMPbPb2040StatErr->Clone("graphYieldPi0PCMPbPb2040StatErrCopy");
	TGraphAsymmErrors* graphYieldPi0PCMPbPb2040SysErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0PCMPbPb2040SysErr->Clone("graphYieldPi0PCMPbPb2040SysErrCopy");
	TGraphAsymmErrors* graphYieldPi0PHOSPbPb2040StatErrCopy 	= (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb2040StatErr->Clone("graphYieldPi0PHOSPbPb2040StatErrCopy");
	TGraphAsymmErrors* graphYieldPi0PHOSPbPb2040SysErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0PHOSPbPb2040SysErr->Clone("graphYieldPi0PHOSPbPb2040SysErrCopy");

	TGraphErrors* graphYieldCombStatPi02040RebinnedHighPtComb 	= NULL;
	TGraphErrors* graphYieldCombSysPi02040RebinnedHighPtComb 	= NULL;
	TGraphErrors* graphChargedPionSpecHighPtStat2040HighPtComb 	= NULL;
	TGraphErrors* graphChargedPionSpecHighPtSyst2040HighPtComb 	= NULL;
	TGraphErrors* graphRatioHighPtChargedPionsComb2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb2040StatErrCopy, graphYieldPi0CombPbPb2040SysErrCopy,
																										  histoChargedPionSpecHighPtStat2040, histoChargedPionSpecHighPtSyst2040,
																									      kTRUE,  kTRUE,
																									      &graphYieldCombStatPi02040RebinnedHighPtComb, &graphYieldCombSysPi02040RebinnedHighPtComb,
																									      &graphChargedPionSpecHighPtStat2040HighPtComb, &graphChargedPionSpecHighPtSyst2040HighPtComb);
	graphRatioHighPtChargedPionsComb2040->Print();

	TGraphErrors* graphYieldCombStatPi02040RebinnedLowPtComb 	= NULL;
	TGraphErrors* graphYieldCombSysPi02040RebinnedLowPtComb 	= NULL;
	TGraphErrors* graphChargedPionSpecLowPtStat2040LowPtComb 	= NULL;
	TGraphErrors* graphChargedPionSpecLowPtSyst2040LowPtComb 	= NULL;
	TGraphErrors* graphRatioLowPtChargedPionsComb2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb2040StatErrCopy, graphYieldPi0CombPbPb2040SysErrCopy,
																										 histoChargedPionSpecLowPtStat2040, histoChargedPionSpecLowPtSyst2040,
																									     kTRUE,  kTRUE,
																									     &graphYieldCombStatPi02040RebinnedLowPtComb, &graphYieldCombSysPi02040RebinnedLowPtComb,
																										 &graphChargedPionSpecLowPtStat2040LowPtComb, &graphChargedPionSpecLowPtSyst2040LowPtComb);
	graphRatioLowPtChargedPionsComb2040->Print();

	TGraphErrors* graphYieldPCMStatPi02040RebinnedHighPtPCM 	= NULL;
	TGraphErrors* graphYieldPCMSysPi02040RebinnedHighPtPCM 		= NULL;
	TGraphErrors* graphChargedPionSpecHighPtStat2040HighPtPCM 	= NULL;
	TGraphErrors* graphChargedPionSpecHighPtSyst2040HighPtPCM 	= NULL;
	TGraphErrors* graphRatioHighPtChargedPionsPCM2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb2040StatErrCopy, graphYieldPi0PCMPbPb2040SysErrCopy,
																										 histoChargedPionSpecHighPtStat2040, histoChargedPionSpecHighPtSyst2040,
																										 kTRUE,  kTRUE,
																									     &graphYieldPCMStatPi02040RebinnedHighPtPCM, &graphYieldPCMSysPi02040RebinnedHighPtPCM,
																										 &graphChargedPionSpecHighPtStat2040HighPtPCM, &graphChargedPionSpecHighPtSyst2040HighPtPCM);
	graphRatioHighPtChargedPionsPCM2040->Print();

	TGraphErrors* graphYieldPCMStatPi02040RebinnedLowPtPCM 		= NULL;
	TGraphErrors* graphYieldPCMSysPi02040RebinnedLowPtPCM 		= NULL;
	TGraphErrors* graphChargedPionSpecLowPtStat2040LowPtPCM 	= NULL;
	TGraphErrors* graphChargedPionSpecLowPtSyst2040LowPtPCM 	= NULL;
	TGraphErrors* graphRatioLowPtChargedPionsPCM2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb2040StatErrCopy, graphYieldPi0PCMPbPb2040SysErrCopy,
																										histoChargedPionSpecLowPtStat2040, histoChargedPionSpecLowPtSyst2040,
																										kTRUE,  kTRUE,
																										&graphYieldPCMStatPi02040RebinnedLowPtPCM, &graphYieldPCMSysPi02040RebinnedLowPtPCM,
																										&graphChargedPionSpecLowPtStat2040LowPtPCM, &graphChargedPionSpecLowPtSyst2040LowPtPCM);
	graphRatioLowPtChargedPionsPCM2040->Print();

	graphYieldPi0PHOSPbPb2040StatErrCopy->RemovePoint(0);
	graphYieldPi0PHOSPbPb2040SysErrCopy->RemovePoint(0);
	TGraphErrors* graphYieldPHOSStatPi02040RebinnedHighPtPHOS 	= NULL;
	TGraphErrors* graphYieldPHOSSysPi02040RebinnedHighPtPHOS 	= NULL;
	TGraphErrors* graphChargedPionSpecHighPtStat2040HighPtPHOS 	= NULL;
	TGraphErrors* graphChargedPionSpecHighPtSyst2040HighPtPHOS 	= NULL;
	TGraphErrors* graphRatioHighPtChargedPionsPHOS2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb2040StatErrCopy, graphYieldPi0PHOSPbPb2040SysErrCopy,
																										  histoChargedPionSpecHighPtStat2040, histoChargedPionSpecHighPtSyst2040,
																										  kTRUE,  kTRUE,
																										  &graphYieldPHOSStatPi02040RebinnedHighPtPHOS, &graphYieldPHOSSysPi02040RebinnedHighPtPHOS,
																										  &graphChargedPionSpecHighPtStat2040HighPtPHOS, &graphChargedPionSpecHighPtSyst2040HighPtPHOS);
	graphRatioHighPtChargedPionsPHOS2040->Print();

	TGraphErrors* graphYieldPHOSStatPi02040RebinnedLowPtPHOS 	= NULL;
	TGraphErrors* graphYieldPHOSSysPi02040RebinnedLowPtPHOS 	= NULL;
	TGraphErrors* graphChargedPionSpecLowPtStat2040LowPtPHOS 	= NULL;
	TGraphErrors* graphChargedPionSpecLowPtSyst2040LowPtPHOS 	= NULL;
	TGraphErrors* graphRatioLowPtChargedPionsPHOS2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb2040StatErrCopy, graphYieldPi0PHOSPbPb2040SysErrCopy,
																										 histoChargedPionSpecLowPtStat2040, histoChargedPionSpecLowPtSyst2040,
																										 kTRUE,  kTRUE,
																										 &graphYieldPHOSStatPi02040RebinnedLowPtPHOS, &graphYieldPHOSSysPi02040RebinnedLowPtPHOS,
																										 &graphChargedPionSpecLowPtStat2040LowPtPHOS, &graphChargedPionSpecLowPtSyst2040LowPtPHOS);
	graphRatioLowPtChargedPionsPHOS2040->Print();


	TGraphErrors* graphYieldCombStatPi02040FullPtPtComb 	= NULL;
	TGraphErrors* graphYieldCombSysPi02040FullPtPtComb 	= NULL;
	TGraphErrors* graphChargedPionSpecFullPtStat2040FullPtComb 	= NULL;
	TGraphErrors* graphChargedPionSpecFullPtSyst2040FullPtComb 	= NULL;
	TGraphErrors* graphRatioFullPtChargedPionsComb2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0CombPbPb2040StatErrCopy, graphYieldPi0CombPbPb2040SysErrCopy,
																										  histoChargedPionSpecFullPtStat2040, histoChargedPionSpecFullPtSyst2040,
																									      kTRUE,  kTRUE,
																									      &graphYieldCombStatPi02040FullPtPtComb, &graphYieldCombSysPi02040FullPtPtComb,
																									      &graphChargedPionSpecFullPtStat2040FullPtComb, &graphChargedPionSpecFullPtSyst2040FullPtComb);
	graphRatioFullPtChargedPionsComb2040->Print();

	TGraphErrors* graphYieldPCMStatPi02040FullPtPtPCM 	= NULL;
	TGraphErrors* graphYieldPCMSysPi02040FullPtPtPCM 		= NULL;
	TGraphErrors* graphChargedPionSpecFullPtStat2040FullPtPCM 	= NULL;
	TGraphErrors* graphChargedPionSpecFullPtSyst2040FullPtPCM 	= NULL;
	TGraphErrors* graphRatioFullPtChargedPionsPCM2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PCMPbPb2040StatErrCopy, graphYieldPi0PCMPbPb2040SysErrCopy,
																										 histoChargedPionSpecFullPtStat2040, histoChargedPionSpecFullPtSyst2040,
																										 kTRUE,  kTRUE,
																									     &graphYieldPCMStatPi02040FullPtPtPCM, &graphYieldPCMSysPi02040FullPtPtPCM,
																									     &graphChargedPionSpecFullPtStat2040FullPtPCM, &graphChargedPionSpecFullPtSyst2040FullPtPCM);
	graphRatioFullPtChargedPionsPCM2040->Print();

	graphYieldPi0PHOSPbPb2040StatErrCopy->RemovePoint(0);
	graphYieldPi0PHOSPbPb2040SysErrCopy->RemovePoint(0);
	TGraphErrors* graphYieldPHOSStatPi02040FullPtPtPHOS 	= NULL;
	TGraphErrors* graphYieldPHOSSysPi02040FullPtPtPHOS 	= NULL;
	TGraphErrors* graphChargedPionSpecFullPtStat2040FullPtPHOS 	= NULL;
	TGraphErrors* graphChargedPionSpecFullPtSyst2040FullPtPHOS 	= NULL;
	TGraphErrors* graphRatioFullPtChargedPionsPHOS2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphYieldPi0PHOSPbPb2040StatErrCopy, graphYieldPi0PHOSPbPb2040SysErrCopy,
																										  histoChargedPionSpecFullPtStat2040, histoChargedPionSpecFullPtSyst2040,
																										  kTRUE,  kTRUE,
																										  &graphYieldPHOSStatPi02040FullPtPtPHOS, &graphYieldPHOSSysPi02040FullPtPtPHOS,
																									      &graphChargedPionSpecFullPtStat2040FullPtPHOS, &graphChargedPionSpecFullPtSyst2040FullPtPHOS);
	graphRatioFullPtChargedPionsPHOS2040->Print();


	cout << "*************************************************************************"<< endl;
	cout << "***************************  Pi0 PbPb 2011 ******************************"<< endl;
	cout << "*************************************************************************"<< endl;

	TFile* fCombResultsLHC11h = new TFile(nameFilePbPbLHC11h.Data());
    TDirectory* directoryUnshiftedSpectra =             (TDirectory*)fCombResultsLHC11h->Get("UnshiftedSpectra");
    //comb
	TGraphAsymmErrors*	graphYieldPi0CombPbPb0010LHC11hStatErr 		= (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldPi0CombPbPb2760GeVStatErrNoShift_0010");
	TGraphAsymmErrors*	graphYieldPi0CombPbPb0010LHC11hSysErr 		= (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldPi0CombPbPb2760GeVSysErrNoShift_0010");
    TGraphAsymmErrors*  graphYieldEtaCombPbPb0010LHC11hStatErr      = (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldEtaCombPbPb2760GeVStatErrNoShift_0010");
    TGraphAsymmErrors*  graphYieldEtaCombPbPb0010LHC11hSysErr       = (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldEtaCombPbPb2760GeVSysErrNoShift_0010");

    TGraphAsymmErrors*	graphYieldPi0CombPbPb2050LHC11hStatErr 		= (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldPi0CombPbPb2760GeVStatErrNoShift_2050");
	TGraphAsymmErrors*	graphYieldPi0CombPbPb2050LHC11hSysErr 		= (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldPi0CombPbPb2760GeVSysErrNoShift_2050");
    TGraphAsymmErrors*  graphYieldEtaCombPbPb2050LHC11hStatErr      = (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldEtaCombPbPb2760GeVStatErrNoShift_2050");
    TGraphAsymmErrors*  graphYieldEtaCombPbPb2050LHC11hSysErr       = (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldEtaCombPbPb2760GeVSysErrNoShift_2050");

    //PCM
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb0010LHC11hStatErr	 	= (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldPi0PCMPbPb2760GeVStatErrNoShift_0010");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb0010LHC11hSysErr 		= (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldPi0PCMPbPb2760GeVSysErrNoShift_0010");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb2040LHC11hStatErr 		= (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldPi0PCMPbPb2760GeVSysErrNoShift_2040");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb2040LHC11hSysErr 		= (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldPi0PCMPbPb2760GeVSysErrNoShift_2040");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb2050LHC11hStatErr 		= (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldPi0PCMPbPb2760GeVSysErrNoShift_2050");
	TGraphAsymmErrors*	graphYieldPi0PCMPbPb2050LHC11hSysErr 		= (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldPi0PCMPbPb2760GeVSysErrNoShift_2050");

    //PCM
    TGraphAsymmErrors*  graphYieldEtaPCMPbPb0010LHC11hStatErr       = (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldEtaPCMPbPb2760GeVStatErrNoShift_0010");
    TGraphAsymmErrors*  graphYieldEtaPCMPbPb0010LHC11hSysErr        = (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldEtaPCMPbPb2760GeVSysErrNoShift_0010");
    TGraphAsymmErrors*  graphYieldEtaPCMPbPb2040LHC11hStatErr       = (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldEtaPCMPbPb2760GeVSysErrNoShift_2040");
    TGraphAsymmErrors*  graphYieldEtaPCMPbPb2040LHC11hSysErr        = (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldEtaPCMPbPb2760GeVSysErrNoShift_2040");
    TGraphAsymmErrors*  graphYieldEtaPCMPbPb2050LHC11hStatErr       = (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldEtaPCMPbPb2760GeVSysErrNoShift_2050");
    TGraphAsymmErrors*  graphYieldEtaPCMPbPb2050LHC11hSysErr        = (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldEtaPCMPbPb2760GeVSysErrNoShift_2050");

    //EMCal
	TGraphAsymmErrors*	graphYieldPi0EMCalPbPb0010LHC11hStatErr 	= (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldPi0EMCalPbPb2760GeVStatErrNoShift_0010");
	TGraphAsymmErrors*	graphYieldPi0EMCalPbPb0010LHC11hSysErr 		= (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldPi0EMCalPbPb2760GeVSysErrNoShift_0010");
    TGraphAsymmErrors*  graphYieldEtaEMCalPbPb0010LHC11hStatErr     = (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldEtaEMCalPbPb2760GeVStatErrNoShift_0010");
    TGraphAsymmErrors*  graphYieldEtaEMCalPbPb0010LHC11hSysErr      = (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldEtaEMCalPbPb2760GeVSysErrNoShift_0010");
	TGraphAsymmErrors*	graphYieldPi0EMCalPbPb2050LHC11hStatErr 	= (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldPi0EMCalPbPb2760GeVStatErrNoShift_2050");
	TGraphAsymmErrors*	graphYieldPi0EMCalPbPb2050LHC11hSysErr 		= (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldPi0EMCalPbPb2760GeVSysErrNoShift_2050");
    TGraphAsymmErrors*  graphYieldEtaEMCalPbPb2050LHC11hStatErr     = (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldEtaEMCalPbPb2760GeVStatErrNoShift_2050");
    TGraphAsymmErrors*  graphYieldEtaEMCalPbPb2050LHC11hSysErr      = (TGraphAsymmErrors*)directoryUnshiftedSpectra->Get("graphInvYieldEtaEMCalPbPb2760GeVSysErrNoShift_2050");

    TGraphErrors * a = NULL;
    TGraphErrors * b = NULL;
    TGraphErrors * c = NULL;
    TGraphErrors * d = NULL;

	cout << "*************************************************************************"<< endl;
	cout << "**************************  PbPb 0-10% LHC11h ***************************"<< endl;
	cout << "*************************************************************************"<< endl;

	TGraphAsymmErrors* graphYieldPi0CombPbPb0010LHC11hStatErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0CombPbPb0010LHC11hStatErr->Clone("graphYieldPi0CombPbPb0010LHC11hStatErrCopy");
	TGraphAsymmErrors* graphYieldPi0CombPbPb0010LHC11hSysErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0CombPbPb0010LHC11hSysErr->Clone("graphYieldPi0CombPbPb0010LHC11hSysErrCopy");

    TGraphAsymmErrors* graphYieldEtaCombPbPb0010LHC11hStatErrCopy       = (TGraphAsymmErrors*) graphYieldEtaCombPbPb0010LHC11hStatErr->Clone("graphYieldEtaCombPbPb0010LHC11hStatErrCopy");
    TGraphAsymmErrors* graphYieldEtaCombPbPb0010LHC11hSysErrCopy        = (TGraphAsymmErrors*) graphYieldEtaCombPbPb0010LHC11hSysErr->Clone("graphYieldEtaCombPbPb0010LHC11hSysErrCopy");

	TGraphAsymmErrors* graphYieldPi0CombPbPb2050LHC11hStatErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0CombPbPb2050LHC11hStatErr->Clone("graphYieldPi0CombPbPb2050LHC11hStatErrCopy");
	TGraphAsymmErrors* graphYieldPi0CombPbPb2050LHC11hSysErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0CombPbPb2050LHC11hSysErr->Clone("graphYieldPi0CombPbPb2050LHC11hSysErrCopy");

    TGraphAsymmErrors* graphYieldEtaCombPbPb2050LHC11hStatErrCopy       = (TGraphAsymmErrors*) graphYieldEtaCombPbPb2050LHC11hStatErr->Clone("graphYieldEtaCombPbPb2050LHC11hStatErrCopy");
    TGraphAsymmErrors* graphYieldEtaCombPbPb2050LHC11hSysErrCopy        = (TGraphAsymmErrors*) graphYieldEtaCombPbPb2050LHC11hSysErr->Clone("graphYieldEtaCombPbPb2050LHC11hSysErrCopy");

    TGraphAsymmErrors* graphYieldPi0PCMPbPb0010LHC11hStatErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0PCMPbPb0010LHC11hStatErr->Clone("graphYieldPi0PCMPbPb0010LHC11hStatErrCopy");
	TGraphAsymmErrors* graphYieldPi0PCMPbPb0010LHC11hSysErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0PCMPbPb0010LHC11hSysErr->Clone("graphYieldPi0PCMPbPb0010LHC11hSysErrCopy");

    TGraphAsymmErrors* graphYieldEtaPCMPbPb0010LHC11hStatErrCopy        = (TGraphAsymmErrors*) graphYieldEtaPCMPbPb0010LHC11hStatErr->Clone("graphYieldEtaPCMPbPb0010LHC11hStatErrCopy");
    TGraphAsymmErrors* graphYieldEtaPCMPbPb0010LHC11hSysErrCopy         = (TGraphAsymmErrors*) graphYieldEtaPCMPbPb0010LHC11hSysErr->Clone("graphYieldEtaPCMPbPb0010LHC11hSysErrCopy");

	TGraphAsymmErrors* graphYieldPi0PCMPbPb2040LHC11hStatErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0PCMPbPb2040LHC11hStatErr->Clone("graphYieldPi0PCMPbPb2040LHC11hStatErrCopy");
	TGraphAsymmErrors* graphYieldPi0PCMPbPb2040LHC11hSysErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0PCMPbPb2040LHC11hSysErr->Clone("graphYieldPi0PCMPbPb2040LHC11hSysErrCopy");

    TGraphAsymmErrors* graphYieldEtaPCMPbPb2040LHC11hStatErrCopy        = (TGraphAsymmErrors*) graphYieldEtaPCMPbPb2040LHC11hStatErr->Clone("graphYieldEtaPCMPbPb2040LHC11hStatErrCopy");
    TGraphAsymmErrors* graphYieldEtaPCMPbPb2040LHC11hSysErrCopy         = (TGraphAsymmErrors*) graphYieldEtaPCMPbPb2040LHC11hSysErr->Clone("graphYieldEtaPCMPbPb2040LHC11hSysErrCopy");

	TGraphAsymmErrors* graphYieldPi0PCMPbPb2050LHC11hStatErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0PCMPbPb2050LHC11hStatErr->Clone("graphYieldPi0PCMPbPb2050LHC11hStatErrCopy");
	TGraphAsymmErrors* graphYieldPi0PCMPbPb2050LHC11hSysErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0PCMPbPb2050LHC11hSysErr->Clone("graphYieldPi0PCMPbPb2050LHC11hSysErrCopy");

    TGraphAsymmErrors* graphYieldEtaPCMPbPb2050LHC11hStatErrCopy        = (TGraphAsymmErrors*) graphYieldEtaPCMPbPb2050LHC11hStatErr->Clone("graphYieldEtaPCMPbPb2050LHC11hStatErrCopy");
    TGraphAsymmErrors* graphYieldEtaPCMPbPb2050LHC11hSysErrCopy         = (TGraphAsymmErrors*) graphYieldEtaPCMPbPb2050LHC11hSysErr->Clone("graphYieldEtaPCMPbPb2050LHC11hSysErrCopy");

	TGraphAsymmErrors* graphYieldPi0EMCalPbPb0010LHC11hStatErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0EMCalPbPb0010LHC11hStatErr->Clone("graphYieldPi0EMCalPbPb0010LHC11hStatErrCopy");
	TGraphAsymmErrors* graphYieldPi0EMCalPbPb0010LHC11hSysErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0EMCalPbPb0010LHC11hSysErr->Clone("graphYieldPi0EMCalPbPb0010LHC11hSysErrCopy");

    TGraphAsymmErrors* graphYieldEtaEMCalPbPb0010LHC11hStatErrCopy      = (TGraphAsymmErrors*) graphYieldEtaEMCalPbPb0010LHC11hStatErr->Clone("graphYieldEtaEMCalPbPb0010LHC11hStatErrCopy");
    TGraphAsymmErrors* graphYieldEtaEMCalPbPb0010LHC11hSysErrCopy       = (TGraphAsymmErrors*) graphYieldEtaEMCalPbPb0010LHC11hSysErr->Clone("graphYieldEtaEMCalPbPb0010LHC11hSysErrCopy");

	TGraphAsymmErrors* graphYieldPi0EMCalPbPb2050LHC11hStatErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0EMCalPbPb2050LHC11hStatErr->Clone("graphYieldPi0EMCalPbPb2050LHC11hStatErrCopy");
	TGraphAsymmErrors* graphYieldPi0EMCalPbPb2050LHC11hSysErrCopy 		= (TGraphAsymmErrors*) graphYieldPi0EMCalPbPb2050LHC11hSysErr->Clone("graphYieldPi0EMCalPbPb2050LHC11hSysErrCopy");

    TGraphAsymmErrors* graphYieldEtaEMCalPbPb2050LHC11hStatErrCopy      = (TGraphAsymmErrors*) graphYieldEtaEMCalPbPb2050LHC11hStatErr->Clone("graphYieldEtaEMCalPbPb2050LHC11hStatErrCopy");
    TGraphAsymmErrors* graphYieldEtaEMCalPbPb2050LHC11hSysErrCopy       = (TGraphAsymmErrors*) graphYieldEtaEMCalPbPb2050LHC11hSysErr->Clone("graphYieldEtaEMCalPbPb2050LHC11hSysErrCopy");


	TGraphErrors* graphRatioHighPtChargedPionsComb0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0CombPbPb0010LHC11hStatErrCopy, graphYieldPi0CombPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedPionSpecHighPtStat0010, histoChargedPionSpecHighPtSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldCombStatPi00010LHC11hRebinnedHighPtComb, &graphYieldCombSysPi00010LHC11hRebinnedHighPtComb,
                                                                      &graphChargedPionSpecHighPtStat0010LHC11hHighPtComb, &graphChargedPionSpecHighPtSyst0010LHC11hHighPtComb);
// 	graphRatioHighPtChargedPionsComb0010LHC11h->Print();

	TGraphErrors* graphRatioLowPtChargedPionsComb0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0CombPbPb0010LHC11hStatErrCopy, graphYieldPi0CombPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedPionSpecLowPtStat0010, histoChargedPionSpecLowPtSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldCombStatPi00010LHC11hRebinnedLowPtComb, &graphYieldCombSysPi00010LHC11hRebinnedLowPtComb,
                                                                      &graphChargedPionSpecLowPtStat0010LHC11hLowPtComb, &graphChargedPionSpecLowPtSyst0010LHC11hLowPtComb);
// 	graphRatioLowPtChargedPionsComb0010LHC11h->Print();

	TGraphErrors* graphRatioHighPtChargedPionsPCM0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0PCMPbPb0010LHC11hStatErrCopy, graphYieldPi0PCMPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedPionSpecHighPtStat0010, histoChargedPionSpecHighPtSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldPCMStatPi00010LHC11hRebinnedHighPtPCM, &graphYieldPCMSysPi00010LHC11hRebinnedHighPtPCM,
                                                                      &graphChargedPionSpecHighPtStat0010LHC11hHighPtPCM, &graphChargedPionSpecHighPtSyst0010LHC11hHighPtPCM);
// 	graphRatioHighPtChargedPionsPCM0010LHC11h->Print();

	TGraphErrors* graphRatioLowPtChargedPionsPCM0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0PCMPbPb0010LHC11hStatErrCopy, graphYieldPi0PCMPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedPionSpecLowPtStat0010, histoChargedPionSpecLowPtSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldPCMStatPi00010LHC11hRebinnedLowPtPCM, &graphYieldPCMSysPi00010LHC11hRebinnedLowPtPCM,
                                                                      &graphChargedPionSpecLowPtStat0010LHC11hLowPtPCM, &graphChargedPionSpecLowPtSyst0010LHC11hLowPtPCM);
// 	graphRatioLowPtChargedPionsPCM0010LHC11h->Print();

	TGraphErrors* graphRatioHighPtChargedPionsPCM2040LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0PCMPbPb2040LHC11hStatErrCopy, graphYieldPi0PCMPbPb2040LHC11hSysErrCopy,
                                                                      histoChargedPionSpecHighPtStat2040, histoChargedPionSpecHighPtSyst2040,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldPCMStatPi02040LHC11hRebinnedHighPtPCM, &graphYieldPCMSysPi02040LHC11hRebinnedHighPtPCM,
                                                                      &graphChargedPionSpecHighPtStat2040LHC11hHighPtPCM, &graphChargedPionSpecHighPtSyst2040LHC11hHighPtPCM);
// 	graphRatioHighPtChargedPionsPCM2040LHC11h->Print();

	TGraphErrors* graphRatioLowPtChargedPionsPCM2040LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0PCMPbPb2040LHC11hStatErrCopy, graphYieldPi0PCMPbPb2040LHC11hSysErrCopy,
                                                                      histoChargedPionSpecLowPtStat2040, histoChargedPionSpecLowPtSyst2040,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldPCMStatPi02040LHC11hRebinnedLowPtPCM, &graphYieldPCMSysPi02040LHC11hRebinnedLowPtPCM,
                                                                      &graphChargedPionSpecLowPtStat2040LHC11hLowPtPCM, &graphChargedPionSpecLowPtSyst2040LHC11hLowPtPCM);
// 	graphRatioLowPtChargedPionsPCM2040LHC11h->Print();

    TGraphErrors* graphRatioHighPtChargedPionsEMCal0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0EMCalPbPb0010LHC11hStatErrCopy, graphYieldPi0EMCalPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedPionSpecHighPtStat0010, histoChargedPionSpecHighPtSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldEMCalStatPi00010LHC11hRebinnedHighPtEMCal, &graphYieldEMCalSysPi00010LHC11hRebinnedHighPtEMCal,
                                                                      &graphChargedPionSpecHighPtStat0010LHC11hHighPtEMCal, &graphChargedPionSpecHighPtSyst0010LHC11hHighPtEMCal);
//     graphRatioHighPtChargedPionsEMCal0010LHC11h->Print();



	TGraphErrors* graphRatioFullPtChargedPionsComb0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0CombPbPb0010LHC11hStatErrCopy, graphYieldPi0CombPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedPionSpecFullPtStat0010, histoChargedPionSpecFullPtSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldCombStatPi00010LHC11hFullPtPtComb, &graphYieldCombSysPi00010LHC11hFullPtPtComb,
                                                                      &graphChargedPionSpecFullPtStat0010LHC11hFullPtComb, &graphChargedPionSpecFullPtSyst0010LHC11hFullPtComb);
// 	graphRatioFullPtChargedPionsComb0010LHC11h->Print();

	TGraphErrors* graphRatioFullPtChargedPionsPCM0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0PCMPbPb0010LHC11hStatErrCopy, graphYieldPi0PCMPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedPionSpecFullPtStat0010, histoChargedPionSpecFullPtSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldPCMStatPi00010LHC11hFullPtPtPCM, &graphYieldPCMSysPi00010LHC11hFullPtPtPCM,
                                                                      &graphChargedPionSpecFullPtStat0010LHC11hFullPtPCM, &graphChargedPionSpecFullPtSyst0010LHC11hFullPtPCM);
// 	graphRatioFullPtChargedPionsPCM0010LHC11h->Print();

	TGraphErrors* graphRatioFullPtChargedPionsComb2050LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0CombPbPb2050LHC11hStatErrCopy, graphYieldPi0CombPbPb2050LHC11hSysErrCopy,
                                                                      histoChargedPionSpecFullPtStat2050, histoChargedPionSpecFullPtSyst2050,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldCombStatPi02050LHC11hFullPtPtComb, &graphYieldCombSysPi02050LHC11hFullPtPtComb,
                                                                      &graphChargedPionSpecFullPtStat2050LHC11hFullPtComb, &graphChargedPionSpecFullPtSyst2050LHC11hFullPtComb);
// 	graphRatioFullPtChargedPionsComb2050LHC11h->Print();

	TGraphErrors* graphRatioFullPtChargedPionsPCM2050LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0PCMPbPb2050LHC11hStatErrCopy, graphYieldPi0PCMPbPb2050LHC11hSysErrCopy,
                                                                      histoChargedPionSpecFullPtStat2050, histoChargedPionSpecFullPtSyst2050,
                                                                      kTRUE,  kTRUE,
                                                                      &a, &b,
                                                                      &c, &d);
// 	graphRatioFullPtChargedPionsPCM2050LHC11h->Print();


    TGraphErrors* graphRatioFullPtChargedPionsPCM2040LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0PCMPbPb2040LHC11hStatErrCopy, graphYieldPi0PCMPbPb2040LHC11hSysErrCopy,
                                                                      histoChargedPionSpecFullPtStat2040, histoChargedPionSpecFullPtSyst2040,
                                                                      kTRUE,  kTRUE,
                                                                      &a, &b,
                                                                      &c, &d);
// 	graphRatioFullPtChargedPionsPCM2040LHC11h->Print();

	TGraphErrors* graphRatioFullPtChargedPionsEMCal0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0EMCalPbPb0010LHC11hStatErrCopy, graphYieldPi0EMCalPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedPionSpecFullPtStat0010, histoChargedPionSpecFullPtSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &a, &b,
                                                                      &c, &d);
// 	graphRatioFullPtChargedPionsEMCal0010LHC11h->Print();

	TGraphErrors* graphRatioFullPtChargedPionsEMCal2050LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0EMCalPbPb2050LHC11hStatErrCopy, graphYieldPi0EMCalPbPb2050LHC11hSysErrCopy,
                                                                      histoChargedPionSpecFullPtStat2050, histoChargedPionSpecFullPtSyst2050,
                                                                      kTRUE,  kTRUE,
                                                                      &a, &b,
                                                                      &c, &d);
// 	graphRatioFullPtChargedPionsEMCal2050LHC11h->Print();





    TGraphErrors* graphRatioFullPtChargedKaonsComb0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldEtaCombPbPb0010LHC11hStatErrCopy, graphYieldEtaCombPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedKaonSpecFullPtStat0010, histoChargedKaonSpecFullPtSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldCombStatEta0010LHC11hFullPtPtComb, &graphYieldCombSysEta0010LHC11hFullPtPtComb,
                                                                      &graphChargedKaonSpecFullPtStat0010LHC11hFullPtComb, &graphChargedKaonSpecFullPtSyst0010LHC11hFullPtComb);
//     graphRatioFullPtChargedKaonsComb0010LHC11h->Print();

    TGraphErrors* graphRatioFullPtChargedKaonsPCM0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldEtaPCMPbPb0010LHC11hStatErrCopy, graphYieldEtaPCMPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedKaonSpecFullPtStat0010, histoChargedKaonSpecFullPtSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldPCMStatEta0010LHC11hFullPtPtPCM, &graphYieldPCMSysEta0010LHC11hFullPtPtPCM,
                                                                      &graphChargedKaonSpecFullPtStat0010LHC11hFullPtPCM, &graphChargedKaonSpecFullPtSyst0010LHC11hFullPtPCM);
//     graphRatioFullPtChargedKaonsPCM0010LHC11h->Print();

    TGraphErrors* graphRatioFullPtChargedKaonsComb2050LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldEtaCombPbPb2050LHC11hStatErrCopy, graphYieldEtaCombPbPb2050LHC11hSysErrCopy,
                                                                      histoChargedKaonSpecFullPtStat2050, histoChargedKaonSpecFullPtSyst2050,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldCombStatEta2050LHC11hFullPtPtComb, &graphYieldCombSysEta2050LHC11hFullPtPtComb,
                                                                      &graphChargedKaonSpecFullPtStat2050LHC11hFullPtComb, &graphChargedKaonSpecFullPtSyst2050LHC11hFullPtComb);
//     graphRatioFullPtChargedKaonsComb2050LHC11h->Print();

    TGraphErrors* graphRatioFullPtChargedKaonsPCM2050LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldEtaPCMPbPb2050LHC11hStatErrCopy, graphYieldEtaPCMPbPb2050LHC11hSysErrCopy,
                                                                      histoChargedKaonSpecFullPtStat2050, histoChargedKaonSpecFullPtSyst2050,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldPCMStatEta2050LHC11hFullPtPtPCM, &graphYieldPCMSysEta2050LHC11hFullPtPtPCM,
                                                                      &graphChargedKaonSpecFullPtStat2050LHC11hFullPtPCM, &graphChargedKaonSpecFullPtSyst2050LHC11hFullPtPCM);
//     graphRatioFullPtChargedKaonsPCM2050LHC11h->Print();


    TGraphErrors* graphRatioFullPtChargedKaonsPCM2040LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldEtaPCMPbPb2040LHC11hStatErrCopy, graphYieldEtaPCMPbPb2040LHC11hSysErrCopy,
                                                                      histoChargedKaonSpecFullPtStat2040, histoChargedKaonSpecFullPtSyst2040,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldPCMStatEta2040LHC11hFullPtPtPCM, &graphYieldPCMSysEta2040LHC11hFullPtPtPCM,
                                                                      &graphChargedKaonSpecFullPtStat2040LHC11hFullPtPCM, &graphChargedKaonSpecFullPtSyst2040LHC11hFullPtPCM);
//     graphRatioFullPtChargedKaonsPCM2040LHC11h->Print();

    TGraphErrors* graphRatioFullPtChargedKaonsEMCal0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldEtaEMCalPbPb0010LHC11hStatErrCopy, graphYieldEtaEMCalPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedKaonSpecFullPtStat0010, histoChargedKaonSpecFullPtSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldEMCalStatEta0010LHC11hFullPtPtEMCal, &graphYieldEMCalSysEta0010LHC11hFullPtPtEMCal,
                                                                      &graphChargedKaonSpecFullPtStat0010LHC11hFullPtEMCal, &graphChargedKaonSpecFullPtSyst0010LHC11hFullPtEMCal);
//     graphRatioFullPtChargedKaonsEMCal0010LHC11h->Print();

    TGraphErrors* graphRatioFullPtChargedKaonsEMCal2050LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldEtaEMCalPbPb2050LHC11hStatErrCopy, graphYieldEtaEMCalPbPb2050LHC11hSysErrCopy,
                                                                      histoChargedKaonSpecFullPtStat2050, histoChargedKaonSpecFullPtSyst2050,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldEMCalStatEta2050LHC11hFullPtPtEMCal, &graphYieldEMCalSysEta2050LHC11hFullPtPtEMCal,
                                                                      &graphChargedKaonSpecFullPtStat2050LHC11hFullPtEMCal, &graphChargedKaonSpecFullPtSyst2050LHC11hFullPtEMCal);
//     graphRatioFullPtChargedKaonsEMCal2050LHC11h->Print();




	cout << "*************************************************************************"<< endl;
	cout << "*********************  PbPb 5TeV charged comparison *********************"<< endl;
	cout << "*************************************************************************"<< endl;

//     TGraphErrors* graphChargedPionSpec5TeVStat0010Copy = (TGraphErrors*)graphChargedPionSpec5TeVStat0010->Clone();
//     TGraphErrors* graphChargedPionSpec5TeVSyst0010Copy = (TGraphErrors*)graphChargedPionSpec5TeVSyst0010->Clone();
//     TGraphErrors* graphChargedPionSpec5TeVStat2050Copy = (TGraphErrors*)graphChargedPionSpec5TeVStat2050->Clone();
//     TGraphErrors* graphChargedPionSpec5TeVSyst2050Copy = (TGraphErrors*)graphChargedPionSpec5TeVSyst2050->Clone();
//
//     while(graphChargedPionSpec5TeVStat0010Copy->GetY()[0] <=0)graphChargedPionSpec5TeVStat0010Copy->RemovePoint(0);
//     while(graphChargedPionSpec5TeVSyst0010Copy->GetY()[0] <=0)graphChargedPionSpec5TeVSyst0010Copy->RemovePoint(0);
//     while(graphChargedPionSpec5TeVStat0010Copy->GetY()[graphChargedPionSpec5TeVStat0010Copy->GetN()-1] <=0)graphChargedPionSpec5TeVStat0010Copy->RemovePoint(graphChargedPionSpec5TeVStat0010Copy->GetN()-1);
//     while(graphChargedPionSpec5TeVSyst0010Copy->GetY()[graphChargedPionSpec5TeVSyst0010Copy->GetN()-1] <=0)graphChargedPionSpec5TeVSyst0010Copy->RemovePoint(graphChargedPionSpec5TeVSyst0010Copy->GetN()-1);
//
//     while(graphChargedPionSpec5TeVStat2050Copy->GetY()[0] <=0)graphChargedPionSpec5TeVStat2050Copy->RemovePoint(0);
//     while(graphChargedPionSpec5TeVSyst2050Copy->GetY()[0] <=0)graphChargedPionSpec5TeVSyst2050Copy->RemovePoint(0);
//     while(graphChargedPionSpec5TeVStat2050Copy->GetY()[graphChargedPionSpec5TeVStat2050Copy->GetN()-1] <=0)graphChargedPionSpec5TeVStat2050Copy->RemovePoint(graphChargedPionSpec5TeVStat2050Copy->GetN()-1);
//     while(graphChargedPionSpec5TeVSyst2050Copy->GetY()[graphChargedPionSpec5TeVSyst2050Copy->GetN()-1] <=0)graphChargedPionSpec5TeVSyst2050Copy->RemovePoint(graphChargedPionSpec5TeVSyst2050Copy->GetN()-1);


	TGraphErrors* graphRatio5TeVto2760GeVChargedPions0010 = CalculateRatioBetweenSpectraWithDifferentBinning(
																										  histoChargedPionSpecFullPtStat0010, histoChargedPionSpecFullPtSyst0010,
                                                                                                          histoChargedPionSpec5TeVStat0010, histoChargedPionSpec5TeVSyst0010,
																									      kTRUE,  kTRUE,
																									      &a, &b,
																									      &c, &d);
	graphRatio5TeVto2760GeVChargedPions0010->Print();

	TGraphErrors* graphRatio5TeVto2760GeVChargedPions2040 = CalculateRatioBetweenSpectraWithDifferentBinning(
																										  histoChargedPionSpecFullPtStat2040, histoChargedPionSpecFullPtSyst2040,
                                                                                                            histoChargedPionSpec5TeVStat2040, histoChargedPionSpec5TeVSyst2040,
																									      kTRUE,  kTRUE,
																									      &a, &b,
																									      &c, &d);
	graphRatio5TeVto2760GeVChargedPions2040->Print();

	TGraphErrors* graphRatio5TeVto2760GeVChargedPions2050 = CalculateRatioBetweenSpectraWithDifferentBinning(
																										  histoChargedPionSpecFullPtStat2050, histoChargedPionSpecFullPtSyst2050,
                                                                                                          histoChargedPionSpec5TeVStat2050, histoChargedPionSpec5TeVSyst2050,
																									      kTRUE,  kTRUE,
																									      &a, &b,
																									      &c, &d);
	graphRatio5TeVto2760GeVChargedPions2050->Print();

	TGraphErrors* graphRatio5TeVChargedPionsComb0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0CombPbPb0010LHC11hStatErrCopy, graphYieldPi0CombPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedPionSpec5TeVStat0010, histoChargedPionSpec5TeVSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldCombStatPi00010LHC11hFullPtPtComb, &graphYieldCombSysPi00010LHC11hFullPtPtComb,
                                                                      &graphChargedPionSpec5TeVStat0010LHC11hFullPtComb, &graphChargedPionSpec5TeVSyst0010LHC11hFullPtComb);
// 	graphRatio5TeVChargedPionsComb0010LHC11h->Print();

	TGraphErrors* graphRatio5TeVChargedPionsComb2050LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0CombPbPb2050LHC11hStatErrCopy, graphYieldPi0CombPbPb2050LHC11hSysErrCopy,
                                                                      histoChargedPionSpec5TeVStat2050, histoChargedPionSpec5TeVSyst2050,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldCombStatPi02050LHC11hFullPtPtComb, &graphYieldCombSysPi02050LHC11hFullPtPtComb,
                                                                      &graphChargedPionSpec5TeVStat2050LHC11hFullPtComb, &graphChargedPionSpec5TeVSyst2050LHC11hFullPtComb);
// 	graphRatio5TeVChargedPionsComb2050LHC11h->Print();

	TGraphErrors* graphRatio5TeVChargedPionsPCM0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0PCMPbPb0010LHC11hStatErrCopy, graphYieldPi0PCMPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedPionSpec5TeVStat0010, histoChargedPionSpec5TeVSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldPCMStatPi00010LHC11hFullPtPtPCM, &graphYieldPCMSysPi00010LHC11hFullPtPtPCM,
                                                                      &graphChargedPionSpec5TeVStat0010LHC11hFullPtPCM, &graphChargedPionSpec5TeVSyst0010LHC11hFullPtPCM);
// 	graphRatio5TeVChargedPionsPCM0010LHC11h->Print();


    TGraphErrors* graphRatio5TeVChargedPionsPCM2040LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0PCMPbPb2040LHC11hStatErrCopy, graphYieldPi0PCMPbPb2040LHC11hSysErrCopy,
                                                                      histoChargedPionSpec5TeVStat2040, histoChargedPionSpec5TeVSyst2040,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldPCMStatPi02040LHC11hFullPtPtPCM, &graphYieldPCMSysPi02040LHC11hFullPtPtPCM,
                                                                      &graphChargedPionSpec5TeVStat2040LHC11hFullPtPCM, &graphChargedPionSpec5TeVSyst2040LHC11hFullPtPCM);
// 	graphRatio5TeVChargedPionsPCM2040LHC11h->Print();

    TGraphErrors* graphRatio5TeVChargedPionsPCM2050LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0PCMPbPb2050LHC11hStatErrCopy, graphYieldPi0PCMPbPb2050LHC11hSysErrCopy,
                                                                      histoChargedPionSpec5TeVStat2050, histoChargedPionSpec5TeVSyst2050,
                                                                      kTRUE,  kTRUE,
                                                                      &a, &b,
                                                                      &c, &d);
// 	graphRatio5TeVChargedPionsPCM2050LHC11h->Print();


	TGraphErrors* graphRatio5TeVChargedPionsEMCal0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0EMCalPbPb0010LHC11hStatErrCopy, graphYieldPi0EMCalPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedPionSpec5TeVStat0010, histoChargedPionSpec5TeVSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &a, &b,
                                                                      &c, &d);
// 	graphRatio5TeVChargedPionsEMCal0010LHC11h->Print();

	TGraphErrors* graphRatio5TeVChargedPionsEMCal2050LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldPi0EMCalPbPb2050LHC11hStatErrCopy, graphYieldPi0EMCalPbPb2050LHC11hSysErrCopy,
                                                                      histoChargedPionSpec5TeVStat2050, histoChargedPionSpec5TeVSyst2050,
                                                                      kTRUE,  kTRUE,
                                                                      &a, &b,
                                                                      &c, &d);
// 	graphRatio5TeVChargedPionsEMCal2050LHC11h->Print();



	TGraphErrors* graphRatio5TeVto2760GeVChargedKaons0010 = CalculateRatioBetweenSpectraWithDifferentBinning(
																										  histoChargedKaonSpecFullPtStat0010, histoChargedKaonSpecFullPtSyst0010,
                                                                                                          histoChargedKaonSpec5TeVStat0010, histoChargedKaonSpec5TeVSyst0010,
																									      kTRUE,  kTRUE,
																									      &a, &b,
																									      &c, &d);
	graphRatio5TeVto2760GeVChargedKaons0010->Print();

	TGraphErrors* graphRatio5TeVto2760GeVChargedKaons2040 = CalculateRatioBetweenSpectraWithDifferentBinning(
																										  histoChargedKaonSpecFullPtStat2040, histoChargedKaonSpecFullPtSyst2040,
                                                                                                            histoChargedKaonSpec5TeVStat2040, histoChargedKaonSpec5TeVSyst2040,
																									      kTRUE,  kTRUE,
																									      &a, &b,
																									      &c, &d);
	graphRatio5TeVto2760GeVChargedKaons2040->Print();

	TGraphErrors* graphRatio5TeVto2760GeVChargedKaons2050 = CalculateRatioBetweenSpectraWithDifferentBinning(
																										  histoChargedKaonSpecFullPtStat2050, histoChargedKaonSpecFullPtSyst2050,
                                                                                                            histoChargedKaonSpec5TeVStat2050, histoChargedKaonSpec5TeVSyst2050,
																									      kTRUE,  kTRUE,
																									      &a, &b,
																									      &c, &d);
	graphRatio5TeVto2760GeVChargedKaons2050->Print();

	TGraphErrors* graphRatio5TeVChargedKaonsComb0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldEtaCombPbPb0010LHC11hStatErrCopy, graphYieldEtaCombPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedKaonSpec5TeVStat0010, histoChargedKaonSpec5TeVSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldCombStatEta0010LHC11hFullPtPtComb, &graphYieldCombSysEta0010LHC11hFullPtPtComb,
                                                                      &graphChargedKaonSpec5TeVStat0010LHC11hFullPtComb, &graphChargedKaonSpec5TeVSyst0010LHC11hFullPtComb);
// 	graphRatio5TeVChargedKaonsComb0010LHC11h->Print();

	TGraphErrors* graphRatio5TeVChargedKaonsComb2050LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldEtaCombPbPb2050LHC11hStatErrCopy, graphYieldEtaCombPbPb2050LHC11hSysErrCopy,
                                                                      histoChargedKaonSpec5TeVStat2050, histoChargedKaonSpec5TeVSyst2050,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldCombStatEta2050LHC11hFullPtPtComb, &graphYieldCombSysEta2050LHC11hFullPtPtComb,
                                                                      &graphChargedKaonSpec5TeVStat2050LHC11hFullPtComb, &graphChargedKaonSpec5TeVSyst2050LHC11hFullPtComb);
// 	graphRatio5TeVChargedKaonsComb2050LHC11h->Print();

	TGraphErrors* graphRatio5TeVChargedKaonsPCM0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldEtaPCMPbPb0010LHC11hStatErrCopy, graphYieldEtaPCMPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedKaonSpec5TeVStat0010, histoChargedKaonSpec5TeVSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldPCMStatEta0010LHC11hFullPtPtPCM, &graphYieldPCMSysEta0010LHC11hFullPtPtPCM,
                                                                      &graphChargedKaonSpec5TeVStat0010LHC11hFullPtPCM, &graphChargedKaonSpec5TeVSyst0010LHC11hFullPtPCM);
// 	graphRatio5TeVChargedKaonsPCM0010LHC11h->Print();


    TGraphErrors* graphRatio5TeVChargedKaonsPCM2040LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldEtaPCMPbPb2040LHC11hStatErrCopy, graphYieldEtaPCMPbPb2040LHC11hSysErrCopy,
                                                                      histoChargedKaonSpec5TeVStat2040, histoChargedKaonSpec5TeVSyst2040,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldPCMStatEta2040LHC11hFullPtPtPCM, &graphYieldPCMSysEta2040LHC11hFullPtPtPCM,
                                                                      &graphChargedKaonSpec5TeVStat2040LHC11hFullPtPCM, &graphChargedKaonSpec5TeVSyst2040LHC11hFullPtPCM);
// 	graphRatio5TeVChargedKaonsPCM2040LHC11h->Print();

    TGraphErrors* graphRatio5TeVChargedKaonsPCM2050LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldEtaPCMPbPb2050LHC11hStatErrCopy, graphYieldEtaPCMPbPb2050LHC11hSysErrCopy,
                                                                      histoChargedKaonSpec5TeVStat2050, histoChargedKaonSpec5TeVSyst2050,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldPCMStatEta2050LHC11hFullPtPtPCM, &graphYieldPCMSysEta2050LHC11hFullPtPtPCM,
                                                                      &graphChargedKaonSpec5TeVStat2050LHC11hFullPtPCM, &graphChargedKaonSpec5TeVSyst2050LHC11hFullPtPCM);
// 	graphRatio5TeVChargedKaonsPCM2050LHC11h->Print();


	TGraphErrors* graphRatio5TeVChargedKaonsEMCal0010LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldEtaEMCalPbPb0010LHC11hStatErrCopy, graphYieldEtaEMCalPbPb0010LHC11hSysErrCopy,
                                                                      histoChargedKaonSpec5TeVStat0010, histoChargedKaonSpec5TeVSyst0010,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldEMCalStatEta0010LHC11hFullPtPtEMCal, &graphYieldEMCalSysEta0010LHC11hFullPtPtEMCal,
                                                                      &graphChargedKaonSpec5TeVStat0010LHC11hFullPtEMCal, &graphChargedKaonSpec5TeVSyst0010LHC11hFullPtEMCal);
// 	graphRatio5TeVChargedKaonsEMCal0010LHC11h->Print();

	TGraphErrors* graphRatio5TeVChargedKaonsEMCal2050LHC11h = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                      graphYieldEtaEMCalPbPb2050LHC11hStatErrCopy, graphYieldEtaEMCalPbPb2050LHC11hSysErrCopy,
                                                                      histoChargedKaonSpec5TeVStat2050, histoChargedKaonSpec5TeVSyst2050,
                                                                      kTRUE,  kTRUE,
                                                                      &graphYieldEMCalStatEta2050LHC11hFullPtPtEMCal, &graphYieldEMCalSysEta2050LHC11hFullPtPtEMCal,
                                                                      &graphChargedKaonSpec5TeVStat2050LHC11hFullPtEMCal, &graphChargedKaonSpec5TeVSyst2050LHC11hFullPtEMCal);
// 	graphRatio5TeVChargedKaonsEMCal2050LHC11h->Print();




	// ***************************************************************************************************************
	// ********************************* Comparison pi0/pi+-, pi0 combined  PbPb with 2011 ***************************
	// ***************************************************************************************************************
	Double_t textsizeLabels1 = 0;
	Double_t textsizeFac1 = 0;
	Double_t textsizeLabels2 = 0;
	Double_t textsizeFac2 = 0;
	Double_t arrayBoundariesX1_4[3];
	Double_t arrayBoundariesY1_4[2];
	Double_t relativeMarginsX[3];
	Double_t relativeMarginsY[3];
	ReturnCorrectValuesForCanvasScaling(1000,500, 2, 1,0.06, 0.005, 0.005,0.09,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

	TCanvas* canvas6PartCompChargedPionsLHC11h = new TCanvas("canvas6PartCompChargedPionsLHC11h","",0,0,1000,500);  // gives the page size
	DrawGammaCanvasSettings( canvas6PartCompChargedPionsLHC11h,  0.13, 0.02, 0.03, 0.06);
	TPad* pad6PartCompChargedPionsLHC11h1 = new TPad("pad6PartCompChargedPionsLHC11h1", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartCompChargedPionsLHC11h1, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[2]);
	pad6PartCompChargedPionsLHC11h1->Draw();
	TPad* pad6PartCompChargedPionsLHC11h3 = new TPad("pad6PartCompChargedPionsLHC11h3", "", arrayBoundariesX1_4[1], arrayBoundariesY1_4[1], arrayBoundariesX1_4[2], arrayBoundariesY1_4[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartCompChargedPionsLHC11h3, relativeMarginsX[1], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[2]);
	pad6PartCompChargedPionsLHC11h3->Draw();

	Double_t margin = relativeMarginsX[0]*0.8*1000;
	ReturnCorrectValuesTextSize(pad6PartCompChargedPionsLHC11h1,textsizeLabels1, textsizeFac1, 22, margin);
	ReturnCorrectValuesTextSize(pad6PartCompChargedPionsLHC11h3,textsizeLabels2, textsizeFac2, 22, margin);

	TH2F * histo2DCompCombinedRatioLHC11h2 = new TH2F("histo2DCompCombinedRatioLHC11h2","histo2DCompCombinedRatioLHC11h2",1000,0.3,40.,1000,0.2,4.	);
	SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatioLHC11h2, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}",0.85*textsizeLabels1, textsizeLabels1,
								  0.85*textsizeLabels1, textsizeLabels1, 0.8,0.25/(textsizeFac1*margin), 512, 505);
	histo2DCompCombinedRatioLHC11h2->GetYaxis()->SetRangeUser(0.,2.1);
	histo2DCompCombinedRatioLHC11h2->GetXaxis()->SetRangeUser(0.,20.);
	histo2DCompCombinedRatioLHC11h2->GetXaxis()->SetLabelOffset(-0.0105);

	TH2F* histo2DCompCombinedRatioLHC11h = new TH2F("histo2DCompCombinedRatioLHC11h","histo2DCompCombinedRatioLHC11h",1000,0.3,40.,1000,0.2,4.	);
	SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatioLHC11h, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.85*textsizeLabels2, textsizeLabels2,
								  0.85*textsizeLabels2, textsizeLabels2, 0.8,0.25/(textsizeFac2*margin), 512, 505);
	histo2DCompCombinedRatioLHC11h->GetXaxis()->SetLabelOffset(-0.0105);
	histo2DCompCombinedRatioLHC11h->GetYaxis()->SetRangeUser(0.6,2.1);
	histo2DCompCombinedRatioLHC11h->GetXaxis()->SetRangeUser(.0,20.);

	pad6PartCompChargedPionsLHC11h1->cd();
	pad6PartCompChargedPionsLHC11h1->SetLogx();
	histo2DCompCombinedRatioLHC11h2->GetXaxis()->SetRangeUser(0.,30.);
	histo2DCompCombinedRatioLHC11h2->GetYaxis()->SetRangeUser(0.,2.1);
	histo2DCompCombinedRatioLHC11h2->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsComb0010, markerStyleCombHighPt, markerSizeComparison, colorCombHighPt , colorCombHighPt);
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsComb0010, markerStyleCombLowPt, markerSizeComparison, colorCombLowPt , colorCombLowPt);
		graphRatioHighPtChargedPionsComb0010->Draw("E1psame");
		graphRatioLowPtChargedPionsComb0010->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsComb0010LHC11h, markerStyleCombHighPt, markerSizeComparison, kBlue+1 , kBlue+1);
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsComb0010LHC11h, markerStyleCombLowPt, markerSizeComparison, kBlue , kBlue);
		graphRatioHighPtChargedPionsComb0010LHC11h->Draw("E1psame");
		graphRatioLowPtChargedPionsComb0010LHC11h->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsComb0010, 33, markerSizeComparison+0.5, kPink+9 , kPink+9);
		graphRatioFullPtChargedPionsComb0010->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsComb0010LHC11h, 33, markerSizeComparison+0.5, kPink+4,kPink+4);
		graphRatioFullPtChargedPionsComb0010LHC11h->Draw("E1psame");


		TLatex *labelPi0CompChargedPionsPbPb0010 = new TLatex(0.16,0.9,collisionSystemCent10.Data());
		SetStyleTLatex( labelPi0CompChargedPionsPbPb0010, 0.85*textsizeLabels1,4);
		labelPi0CompChargedPionsPbPb0010->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1, kGray, 2);

		TLegend* legendPi0CompChargedPionsPbPb0010 = new TLegend(0.12,0.88-(0.035*3),0.95,0.88);
		legendPi0CompChargedPionsPbPb0010->SetFillColor(0);
		legendPi0CompChargedPionsPbPb0010->SetLineColor(0);
		legendPi0CompChargedPionsPbPb0010->SetNColumns(2);
		legendPi0CompChargedPionsPbPb0010->SetTextSize(0.85*textsizeLabels1);
		legendPi0CompChargedPionsPbPb0010->AddEntry(graphRatioLowPtChargedPionsComb0010,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (2010)","p");
		legendPi0CompChargedPionsPbPb0010->AddEntry(graphRatioHighPtChargedPionsComb0010,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (2010)","p");
		legendPi0CompChargedPionsPbPb0010->AddEntry(graphRatioLowPtChargedPionsComb0010LHC11h,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (2011)","p");
		legendPi0CompChargedPionsPbPb0010->AddEntry(graphRatioHighPtChargedPionsComb0010LHC11h,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (2011)","p");
// 		legendPi0CompChargedPionsPbPb0010->AddEntry((TObject*)0,"charged ref: PWGLF-258","");
// 		legendPi0CompChargedPionsPbPb0010->AddEntry((TObject*)0,"","");
		legendPi0CompChargedPionsPbPb0010->AddEntry(graphRatioFullPtChargedPionsComb0010,"#pi^{0}/#pi^{#pm} (2010)","p");
		legendPi0CompChargedPionsPbPb0010->AddEntry(graphRatioFullPtChargedPionsComb0010LHC11h,"#pi^{0}/#pi^{#pm} (2011)","p");

		legendPi0CompChargedPionsPbPb0010->Draw();
		DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);

	histo2DCompCombinedRatioLHC11h2->Draw("axis,same");
	pad6PartCompChargedPionsLHC11h1->Update();
	pad6PartCompChargedPionsLHC11h3->cd();
	pad6PartCompChargedPionsLHC11h3->SetLogx();
	histo2DCompCombinedRatioLHC11h->GetXaxis()->SetRangeUser(0.,30.);
	histo2DCompCombinedRatioLHC11h->GetYaxis()->SetRangeUser(0.,2.1);
	histo2DCompCombinedRatioLHC11h->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsComb2040, markerStyleCombHighPt, markerSizeComparison, colorCombHighPt , colorCombHighPt);
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsComb2040, markerStyleCombLowPt, markerSizeComparison, colorCombLowPt , colorCombLowPt);
		graphRatioHighPtChargedPionsComb2040->Draw("E1psame");
		graphRatioLowPtChargedPionsComb2040->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsComb2040, 33, markerSizeComparison+0.5, kPink+9 , kPink+9);
		graphRatioFullPtChargedPionsComb2040->Draw("E1psame");

		TLatex *labelPi0CompChargedPionsLHC11hPbPb2040 = new TLatex(0.16-relativeMarginsX[0],0.9,collisionSystemSemiCent.Data());
		SetStyleTLatex( labelPi0CompChargedPionsLHC11hPbPb2040, 0.85*textsizeLabels2,4);
		labelPi0CompChargedPionsLHC11hPbPb2040->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1, kGray, 2);


	histo2DCompCombinedRatioLHC11h->Draw("axis,same");
	pad6PartCompChargedPionsLHC11h3->Update();

	canvas6PartCompChargedPionsLHC11h->Update();
	canvas6PartCompChargedPionsLHC11h->SaveAs(Form("%s/ComparisonChargedToNeutralCombinedMeas_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

	delete pad6PartCompChargedPionsLHC11h1;
	delete pad6PartCompChargedPionsLHC11h3;
	delete canvas6PartCompChargedPionsLHC11h;



	// ***************************************************************************************************************
	// ************************************ Comparison pi0/pi+-, pi0 PCM, PHOS PbPb 2011 *****************************
	// ***************************************************************************************************************
	ReturnCorrectValuesForCanvasScaling(1000,500, 2, 1,0.06, 0.01, 0.01,0.09,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

	TCanvas* canvas6PartCompChargedIndPionsLHC11h = new TCanvas("canvas6PartCompChargedIndPionsLHC11h","",0,0,1000,500);  // gives the page size
	DrawGammaCanvasSettings( canvas6PartCompChargedIndPionsLHC11h,  0.13, 0.04, 0.05, 0.06);

	TPad* pad6PartCompChargedIndPionsLHC11h1 = new TPad("pad6PartCompChargedIndPionsLHC11h1", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartCompChargedIndPionsLHC11h1, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[2]);
	pad6PartCompChargedIndPionsLHC11h1->Draw();

	TPad* pad6PartCompChargedIndPionsLHC11h3 = new TPad("pad6PartCompChargedIndPionsLHC11h3", "", arrayBoundariesX1_4[1], arrayBoundariesY1_4[1], arrayBoundariesX1_4[2], arrayBoundariesY1_4[0],-1, -1, -2);
	DrawGammaPadSettings( pad6PartCompChargedIndPionsLHC11h3, relativeMarginsX[1], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[2]);
	pad6PartCompChargedIndPionsLHC11h3->Draw();

	ReturnCorrectValuesTextSize(pad6PartCompChargedIndPionsLHC11h1,textsizeLabels1, textsizeFac1, 22, margin);
	ReturnCorrectValuesTextSize(pad6PartCompChargedIndPionsLHC11h3,textsizeLabels2, textsizeFac2, 22, margin);

	pad6PartCompChargedIndPionsLHC11h1->cd();
	pad6PartCompChargedIndPionsLHC11h1->SetLogx();
	histo2DCompCombinedRatioLHC11h2->GetXaxis()->SetRangeUser(0.,30.);
	histo2DCompCombinedRatioLHC11h2->GetYaxis()->SetRangeUser(0.,2.1); //2.1
	histo2DCompCombinedRatioLHC11h2->DrawCopy();


		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM0010LHC11h, markerStylePCMHighPt, markerSizeComparison, kBlue+1 , kBlue+1);
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM0010LHC11h, markerStylePCMLowPt, markerSizeComparison, kBlue, kBlue);
		graphRatioHighPtChargedPionsPCM0010LHC11h->Draw("E1psame");
		graphRatioLowPtChargedPionsPCM0010LHC11h->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsEMCal0010LHC11h, markerStylePHOSHighPt, markerSizeComparison, kGreen+2 , kGreen+2);
		graphRatioHighPtChargedPionsEMCal0010LHC11h->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM0010, markerStylePCMHighPt, markerSizeComparison, colorPCMLowPt , colorPCMLowPt);
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM0010, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
		graphRatioHighPtChargedPionsPCM0010->Draw("E1psame");
		graphRatioLowPtChargedPionsPCM0010->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS0010, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS0010, markerStylePHOSLowPt, markerSizeComparison, colorPHOSLowPt, colorPHOSLowPt);
		graphRatioHighPtChargedPionsPHOS0010->Draw("E1psame");
		graphRatioLowPtChargedPionsPHOS0010->Draw("E1psame");

		TLatex *labelPi0CompChargedPionsPbPbLHC11h0010 = new TLatex(0.16,0.92,collisionSystemCent10.Data());
		SetStyleTLatex( labelPi0CompChargedPionsPbPbLHC11h0010, 0.85*textsizeLabels2,4);
		labelPi0CompChargedPionsPbPbLHC11h0010->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray, 2);

		TLegend* legendPi0CompChargedIndPionsPbPb0010 = new TLegend(0.13,0.9-(0.035*6),0.95,0.90);
		legendPi0CompChargedIndPionsPbPb0010->SetFillColor(0);
		legendPi0CompChargedIndPionsPbPb0010->SetLineColor(0);
		legendPi0CompChargedIndPionsPbPb0010->SetNColumns(2);
		legendPi0CompChargedIndPionsPbPb0010->SetTextSize(0.85*textsizeLabels1);
		legendPi0CompChargedIndPionsPbPb0010->AddEntry((TObject*)0,"2010 data:","");
		legendPi0CompChargedIndPionsPbPb0010->AddEntry((TObject*)0,"","");
		legendPi0CompChargedIndPionsPbPb0010->AddEntry(graphRatioLowPtChargedPionsPCM0010,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PCM)","p");
		legendPi0CompChargedIndPionsPbPb0010->AddEntry(graphRatioHighPtChargedPionsPCM0010,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PCM)","p");
		legendPi0CompChargedIndPionsPbPb0010->AddEntry(graphRatioLowPtChargedPionsPHOS0010,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PHOS)","p");
		legendPi0CompChargedIndPionsPbPb0010->AddEntry(graphRatioHighPtChargedPionsPHOS0010,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PHOS)","p");
		legendPi0CompChargedIndPionsPbPb0010->AddEntry((TObject*)0,"2011 data:","");
		legendPi0CompChargedIndPionsPbPb0010->AddEntry((TObject*)0,"","");
		legendPi0CompChargedIndPionsPbPb0010->AddEntry(graphRatioLowPtChargedPionsPCM0010LHC11h,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PCM)","p");
		legendPi0CompChargedIndPionsPbPb0010->AddEntry(graphRatioHighPtChargedPionsPCM0010LHC11h,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PCM)","p");
		legendPi0CompChargedIndPionsPbPb0010->AddEntry((TObject*)0,"","");
		legendPi0CompChargedIndPionsPbPb0010->AddEntry(graphRatioHighPtChargedPionsEMCal0010LHC11h,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (EMCal)","p");
		legendPi0CompChargedIndPionsPbPb0010->Draw();

	histo2DCompCombinedRatioLHC11h2->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h3->cd();
	pad6PartCompChargedIndPionsLHC11h3->SetLogx();
	histo2DCompCombinedRatioLHC11h->GetXaxis()->SetRangeUser(-0.25,15.);
	histo2DCompCombinedRatioLHC11h->GetYaxis()->SetRangeUser(0.,2.1);
	histo2DCompCombinedRatioLHC11h->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM2040LHC11h, markerStylePCMHighPt, markerSizeComparison, kBlue+1 , kBlue+1);
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM2040LHC11h, markerStylePCMLowPt, markerSizeComparison, kBlue, kBlue);
		graphRatioHighPtChargedPionsPCM2040LHC11h->Draw("E1psame");
		graphRatioLowPtChargedPionsPCM2040LHC11h->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM2040, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM2040, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
		graphRatioHighPtChargedPionsPCM2040->Draw("E1psame");
		graphRatioLowPtChargedPionsPCM2040->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS2040, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS2040, markerStylePHOSLowPt, markerSizeComparison, colorPHOSLowPt, colorPHOSLowPt);
		graphRatioHighPtChargedPionsPHOS2040->Draw("E1psame");
		graphRatioLowPtChargedPionsPHOS2040->Draw("E1psame");

		TLatex *labelPi0CompChargedPionsPbPbLHC11h2040 = new TLatex(0.16-relativeMarginsX[0],0.92,collisionSystemSemiCent.Data());
		SetStyleTLatex( labelPi0CompChargedPionsPbPbLHC11h2040, 0.85*textsizeLabels2,4);
		labelPi0CompChargedPionsPbPbLHC11h2040->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray,2 );

	histo2DCompCombinedRatioLHC11h->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h3->Update();

	canvas6PartCompChargedIndPionsLHC11h->Update();
	canvas6PartCompChargedIndPionsLHC11h->SaveAs(Form("%s/ComparisonChargedToNeutralAll_splitPt_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));


    pad6PartCompChargedIndPionsLHC11h1->cd();
	pad6PartCompChargedIndPionsLHC11h1->SetLogx();
	histo2DCompCombinedRatioLHC11h2->GetXaxis()->SetRangeUser(0.,30.);
	histo2DCompCombinedRatioLHC11h2->GetYaxis()->SetRangeUser(0.,2.1); //2.1
	histo2DCompCombinedRatioLHC11h2->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsEMCal0010LHC11h, markerStylePHOSHighPt, markerSizeComparison, kGreen+2 , kGreen+2);
		graphRatioHighPtChargedPionsEMCal0010LHC11h->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS0010, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS0010, markerStylePHOSLowPt, markerSizeComparison, colorPHOSLowPt, colorPHOSLowPt);
		graphRatioHighPtChargedPionsPHOS0010->Draw("E1psame");
		graphRatioLowPtChargedPionsPHOS0010->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsPCM0010, 33, markerSizeComparison+0.5, kPink+9 , kPink+9);
		graphRatioFullPtChargedPionsPCM0010->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsPCM0010LHC11h, 33, markerSizeComparison+0.5, kPink+4,kPink+4);
		graphRatioFullPtChargedPionsPCM0010LHC11h->Draw("E1psame");

        labelPi0CompChargedPionsPbPbLHC11h0010->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray, 2);

		TLegend* legendPi0CompChargedIndPionsPbPb0010_fullPt = new TLegend(0.13,0.67,0.95,0.90);
		legendPi0CompChargedIndPionsPbPb0010_fullPt->SetFillColor(0);
		legendPi0CompChargedIndPionsPbPb0010_fullPt->SetLineColor(0);
		legendPi0CompChargedIndPionsPbPb0010_fullPt->SetNColumns(2);
		legendPi0CompChargedIndPionsPbPb0010_fullPt->SetTextSize(0.85*textsizeLabels1);
// 		legendPi0CompChargedIndPionsPbPb0010_fullPt->AddEntry((TObject*)0,"2010 data:","");
// 		legendPi0CompChargedIndPionsPbPb0010_fullPt->AddEntry((TObject*)0,"","");
		legendPi0CompChargedIndPionsPbPb0010_fullPt->AddEntry(graphRatioFullPtChargedPionsPCM0010,"#pi^{0}/#pi^{#pm} (PCM 2010)","p");
		legendPi0CompChargedIndPionsPbPb0010_fullPt->AddEntry((TObject*)0,"","");
// 		legendPi0CompChargedIndPionsPbPb0010_fullPt->AddEntry((TObject*)0,"2011 data:","");
// 		legendPi0CompChargedIndPionsPbPb0010_fullPt->AddEntry((TObject*)0,"","");
		legendPi0CompChargedIndPionsPbPb0010_fullPt->AddEntry(graphRatioFullPtChargedPionsPCM0010LHC11h,"#pi^{0}/#pi^{#pm} (PCM 2011)","p");
// 		legendPi0CompChargedIndPionsPbPb0010_fullPt->AddEntry((TObject*)0,"","");
		legendPi0CompChargedIndPionsPbPb0010_fullPt->AddEntry(graphRatioLowPtChargedPionsPCM0010LHC11h,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (PCM)","p");
		legendPi0CompChargedIndPionsPbPb0010_fullPt->AddEntry(graphRatioHighPtChargedPionsPCM0010LHC11h,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (PCM)","p");
		legendPi0CompChargedIndPionsPbPb0010_fullPt->AddEntry((TObject*)0,"","");
		legendPi0CompChargedIndPionsPbPb0010_fullPt->AddEntry(graphRatioHighPtChargedPionsEMCal0010LHC11h,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (EMCal)","p");
		legendPi0CompChargedIndPionsPbPb0010_fullPt->Draw();


	histo2DCompCombinedRatioLHC11h2->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h3->cd();
	pad6PartCompChargedIndPionsLHC11h3->SetLogx();
	histo2DCompCombinedRatioLHC11h->GetXaxis()->SetRangeUser(-0.25,15.);
	histo2DCompCombinedRatioLHC11h->GetYaxis()->SetRangeUser(0.,2.1);
	histo2DCompCombinedRatioLHC11h->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPHOS2040, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPHOS2040, markerStylePHOSLowPt, markerSizeComparison, colorPHOSLowPt, colorPHOSLowPt);
		graphRatioHighPtChargedPionsPHOS2040->Draw("E1psame");
		graphRatioLowPtChargedPionsPHOS2040->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsPCM2040, 33, markerSizeComparison+0.5, kPink+9 , kPink+9);
		graphRatioFullPtChargedPionsPCM2040->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsPCM2040LHC11h, 33, markerSizeComparison+0.5, kPink+4,kPink+4);
		graphRatioFullPtChargedPionsPCM2040LHC11h->Draw("E1psame");

		labelPi0CompChargedPionsPbPbLHC11h2040->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray,2 );

	histo2DCompCombinedRatioLHC11h->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h3->Update();
	canvas6PartCompChargedIndPionsLHC11h->Update();
	canvas6PartCompChargedIndPionsLHC11h->SaveAs(Form("%s/ComparisonChargedToNeutralAll_fullPt_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));



    pad6PartCompChargedIndPionsLHC11h1->cd();
	pad6PartCompChargedIndPionsLHC11h1->SetLogx();
	TH2F * histo2DCompCharged5to2760 = new TH2F("histo2DCompCharged5to2760","histo2DCompCharged5to2760",1000,0.,40.,1000,0.1,4.	);
	SetStyleHistoTH2ForGraphs(histo2DCompCharged5to2760, "#it{p}_{T} (GeV/#it{c})","2760 GeV/5 TeV",0.85*textsizeLabels1, textsizeLabels1,
								  0.85*textsizeLabels1, textsizeLabels1, 0.8,0.22/(textsizeFac1*margin), 512, 505);
	histo2DCompCharged5to2760->GetYaxis()->SetRangeUser(0.3,1.3);
	histo2DCompCharged5to2760->GetXaxis()->SetRangeUser(0.3,20.);
	histo2DCompCharged5to2760->GetXaxis()->SetLabelOffset(-0.0105);
    histo2DCompCharged5to2760->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatio5TeVto2760GeVChargedPions0010, markerStylePHOSHighPt, markerSizeComparison, kGray+1 , kGray+1);
		graphRatio5TeVto2760GeVChargedPions0010->Draw("E1psame");
        //Comb to 5TeV
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedPionsComb0010LHC11h, 33, markerSizeComparison+0.5,  kCyan+2 , kCyan+2);
		graphRatio5TeVChargedPionsComb0010LHC11h->Draw("E1psame");
//         //PCM to 5TeV
// 		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedPionsPCM0010LHC11h, 33, markerSizeComparison+0.5,  kCyan+2 , kCyan+2);
// 		graphRatio5TeVChargedPionsPCM0010LHC11h->Draw("E1psame");
        //EMCAL to 5TeV
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedPionsEMCal0010LHC11h, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
		graphRatio5TeVChargedPionsEMCal0010LHC11h->Draw("E1psame");

        labelPi0CompChargedPionsPbPbLHC11h0010->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray, 2);

		TLegend* legendCompCharged5and2760 = new TLegend(0.13,0.77,0.95,0.90);
		legendCompCharged5and2760->SetFillColor(0);
		legendCompCharged5and2760->SetLineColor(0);
		legendCompCharged5and2760->SetTextSize(0.85*textsizeLabels1);
		legendCompCharged5and2760->AddEntry(graphRatio5TeVto2760GeVChargedPions0010,"#pi^{#pm}","p");
		legendCompCharged5and2760->AddEntry(graphRatio5TeVChargedPionsComb0010LHC11h,"#pi^{0}/#pi^{#pm} combined","p");
// 		legendCompCharged5and2760->AddEntry(graphRatio5TeVChargedPionsPCM0010LHC11h,"#pi^{0}/#pi^{#pm} PCM","p");
		legendCompCharged5and2760->AddEntry(graphRatio5TeVChargedPionsEMCal0010LHC11h,"#pi^{0}/#pi^{#pm} EMCal","p");
		legendCompCharged5and2760->Draw();


	histo2DCompCombinedRatioLHC11h2->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h3->cd();
	pad6PartCompChargedIndPionsLHC11h3->SetLogx();
	TH2F* histo2DCompCharged5to27602 = new TH2F("histo2DCompCharged5to27602","histo2DCompCharged5to27602",1000,0.,40.,1000,0.1,4.	);
	SetStyleHistoTH2ForGraphs(histo2DCompCharged5to27602, "#it{p}_{T} (GeV/#it{c})","2760 GeV/5 TeV", 0.85*textsizeLabels2, textsizeLabels2,
								  0.85*textsizeLabels2, textsizeLabels2, 0.8,0.22/(textsizeFac2*margin), 512, 505);
	histo2DCompCharged5to27602->GetXaxis()->SetLabelOffset(-0.0105);
	histo2DCompCharged5to27602->GetYaxis()->SetRangeUser(0.3,1.3);
	histo2DCompCharged5to27602->GetXaxis()->SetRangeUser(0.3,20.);
    histo2DCompCharged5to27602->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatio5TeVto2760GeVChargedPions2040, markerStylePHOSHighPt, markerSizeComparison,kGray , kGray);
		graphRatio5TeVto2760GeVChargedPions2040->Draw("E1psame");
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVto2760GeVChargedPions2050, markerStylePHOSHighPt, markerSizeComparison, kGray+1 , kGray+1);
		graphRatio5TeVto2760GeVChargedPions2050->Draw("E1psame");
        //comb to 5TeV
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedPionsComb2050LHC11h, 33, markerSizeComparison+0.5,  kCyan+2 , kCyan+2);
		graphRatio5TeVChargedPionsComb2050LHC11h->Draw("E1psame");
        //PCM to 5TeV
// 		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedPionsPCM2050LHC11h, 33, markerSizeComparison+0.5,  kCyan+2 , kCyan+2);
// 		graphRatio5TeVChargedPionsPCM2050LHC11h->Draw("E1psame");
        //EMCAL to 5TeV
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedPionsEMCal2050LHC11h, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
		graphRatio5TeVChargedPionsEMCal2050LHC11h->Draw("E1psame");

		TLatex *labelPi0CompChargedPionsPbPbLHC11h2050 = new TLatex(0.16-relativeMarginsX[0],0.92,collisionSystem2050.Data());
		SetStyleTLatex( labelPi0CompChargedPionsPbPbLHC11h2050, 0.85*textsizeLabels2,4);
		labelPi0CompChargedPionsPbPbLHC11h2050->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray,2 );


		TLegend* legendCompCharged5and27602 = new TLegend(0.11,0.8,0.5,0.90);
		legendCompCharged5and27602->SetFillColor(0);
		legendCompCharged5and27602->SetLineColor(0);
		legendCompCharged5and27602->SetTextSize(0.85*textsizeLabels1);
		legendCompCharged5and27602->AddEntry(graphRatio5TeVto2760GeVChargedPions2050,"#pi^{#pm} (20-50%)","p");
		legendCompCharged5and27602->AddEntry(graphRatio5TeVto2760GeVChargedPions2040,"#pi^{#pm} (20-40%)","p");
		legendCompCharged5and27602->Draw();


	histo2DCompCombinedRatioLHC11h->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h3->Update();
	canvas6PartCompChargedIndPionsLHC11h->Update();
	canvas6PartCompChargedIndPionsLHC11h->SaveAs(Form("%s/ComparisonChargedPion5TeVTo2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));


    pad6PartCompChargedIndPionsLHC11h1->cd();
	pad6PartCompChargedIndPionsLHC11h1->SetLogx();
    histo2DCompCharged5to2760->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatio5TeVto2760GeVChargedKaons0010, markerStylePHOSHighPt, markerSizeComparison, kGray+1 , kGray+1);
		graphRatio5TeVto2760GeVChargedKaons0010->Draw("E1psame");
        //comb to 5TeV
        TGraphErrors *graphRatio5TeVChargedKaonsComb0010LHC11h_lowPt = (TGraphErrors*)graphRatio5TeVChargedKaonsComb0010LHC11h->Clone();
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedKaonsComb0010LHC11h_lowPt, 27, markerSizeComparison+0.5,  kCyan+2 , kCyan+2);
		graphRatio5TeVChargedKaonsComb0010LHC11h_lowPt->Draw("E1psame");
        graphRatio5TeVChargedKaonsComb0010LHC11h->RemovePoint(0);
        graphRatio5TeVChargedKaonsComb0010LHC11h->RemovePoint(0);
        graphRatio5TeVChargedKaonsComb0010LHC11h->RemovePoint(0);
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedKaonsComb0010LHC11h, 33, markerSizeComparison+0.5,  kCyan+2 , kCyan+2);
		graphRatio5TeVChargedKaonsComb0010LHC11h->Draw("E1psame");
        //PCM to 5TeV
// 		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedKaonsPCM0010LHC11h, 33, markerSizeComparison+0.5,  kCyan+2 , kCyan+2);
// 		graphRatio5TeVChargedKaonsPCM0010LHC11h->Draw("E1psame");
        //EMCAL to 5TeV
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedKaonsEMCal0010LHC11h, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
		graphRatio5TeVChargedKaonsEMCal0010LHC11h->Draw("E1psame");

        labelPi0CompChargedPionsPbPbLHC11h0010->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray, 2);

		TLegend* legendCompChargedKaon5and2760 = new TLegend(0.11,0.73,0.93,0.90);
		legendCompChargedKaon5and2760->SetFillColor(0);
		legendCompChargedKaon5and2760->SetLineColor(0);
		legendCompChargedKaon5and2760->SetTextSize(0.85*textsizeLabels1);
		legendCompChargedKaon5and2760->AddEntry(graphRatio5TeVto2760GeVChargedKaons0010,"K^{#pm}","p");
		legendCompChargedKaon5and2760->AddEntry(graphRatio5TeVChargedKaonsComb0010LHC11h_lowPt,"#eta/K^{#pm} combined (not to be considered, ","p");
        legendCompChargedKaon5and2760->AddEntry((TObject*)0,"        because of mass difference)","");
		legendCompChargedKaon5and2760->AddEntry(graphRatio5TeVChargedKaonsComb0010LHC11h,"#eta/K^{#pm} combined","p");
// 		legendCompChargedKaon5and2760->AddEntry(graphRatio5TeVChargedKaonsPCM0010LHC11h,"#eta/K^{#pm} PCM","p");
		legendCompChargedKaon5and2760->AddEntry(graphRatio5TeVChargedKaonsEMCal0010LHC11h,"#eta/K^{#pm} EMCal","p");
		legendCompChargedKaon5and2760->Draw();


	histo2DCompCombinedRatioLHC11h2->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h3->cd();
	pad6PartCompChargedIndPionsLHC11h3->SetLogx();
    histo2DCompCharged5to27602->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatio5TeVto2760GeVChargedKaons2040, markerStylePHOSHighPt, markerSizeComparison,kGray , kGray);
		graphRatio5TeVto2760GeVChargedKaons2040->Draw("E1psame");
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVto2760GeVChargedKaons2050, markerStylePHOSHighPt, markerSizeComparison, kGray+1 , kGray+1);
		graphRatio5TeVto2760GeVChargedKaons2050->Draw("E1psame");

        TGraphErrors *graphRatio5TeVChargedKaonsComb2050LHC11h_lowPt = (TGraphErrors*)graphRatio5TeVChargedKaonsComb2050LHC11h->Clone();
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedKaonsComb2050LHC11h_lowPt, 27, markerSizeComparison+0.5,  kCyan+2 , kCyan+2);
		graphRatio5TeVChargedKaonsComb2050LHC11h_lowPt->Draw("E1psame");
        graphRatio5TeVChargedKaonsComb2050LHC11h->RemovePoint(0);
        graphRatio5TeVChargedKaonsComb2050LHC11h->RemovePoint(0);
        graphRatio5TeVChargedKaonsComb2050LHC11h->RemovePoint(0);
        //PCM to 5TeV
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedKaonsComb2050LHC11h, 33, markerSizeComparison+0.5,  kCyan+2 , kCyan+2);
		graphRatio5TeVChargedKaonsComb2050LHC11h->Draw("E1psame");
         //PCM to 5TeV
// 		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedKaonsPCM2050LHC11h, 33, markerSizeComparison+0.5,  kCyan+2 , kCyan+2);
// 		graphRatio5TeVChargedKaonsPCM2050LHC11h->Draw("E1psame");
        //EMCAL to 5TeV
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedKaonsEMCal2050LHC11h, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
		graphRatio5TeVChargedKaonsEMCal2050LHC11h->Draw("E1psame");

		TLatex *labelPi0CompChargedKaonsPbPbLHC11h2050 = new TLatex(0.16-relativeMarginsX[0],0.92,collisionSystem2050.Data());
		SetStyleTLatex( labelPi0CompChargedKaonsPbPbLHC11h2050, 0.85*textsizeLabels2,4);
		labelPi0CompChargedKaonsPbPbLHC11h2050->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray,2 );


		TLegend* legendCompChargedKaon5and27602 = new TLegend(0.13,0.77,0.5,0.90);
		legendCompChargedKaon5and27602->SetFillColor(0);
		legendCompChargedKaon5and27602->SetLineColor(0);
		legendCompChargedKaon5and27602->SetTextSize(0.85*textsizeLabels1);
		legendCompChargedKaon5and27602->AddEntry(graphRatio5TeVto2760GeVChargedKaons2040,"K^{#pm} (20-40%)","p");
		legendCompChargedKaon5and27602->AddEntry(graphRatio5TeVto2760GeVChargedKaons2050,"K^{#pm} (20-50%)","p");
		legendCompChargedKaon5and27602->Draw();


	histo2DCompCombinedRatioLHC11h->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h3->Update();
	canvas6PartCompChargedIndPionsLHC11h->Update();
	canvas6PartCompChargedIndPionsLHC11h->SaveAs(Form("%s/ComparisonChargedKaon5TeVTo2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));


    pad6PartCompChargedIndPionsLHC11h1->cd();
	pad6PartCompChargedIndPionsLHC11h1->SetLogx();
	histo2DCompCombinedRatioLHC11h2->GetXaxis()->SetRangeUser(0.,30.);
	histo2DCompCombinedRatioLHC11h2->GetYaxis()->SetRangeUser(0.,2.1); //2.1
	histo2DCompCombinedRatioLHC11h2->DrawCopy();

        //EMCAL to 2.76TeV
		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsEMCal0010LHC11h, markerStylePHOSHighPt, markerSizeComparison, kGreen+2 , kGreen+2);
		graphRatioFullPtChargedPionsEMCal0010LHC11h->Draw("E1psame");
        //EMCAL to 5TeV
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedPionsEMCal0010LHC11h, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
		graphRatio5TeVChargedPionsEMCal0010LHC11h->Draw("E1psame");
        //PCM to 2.76TeV
		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsPCM0010LHC11h, 33, markerSizeComparison+0.5, kPink+4,kPink+4);
		graphRatioFullPtChargedPionsPCM0010LHC11h->Draw("E1psame");
        //PCM to 5TeV
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedPionsPCM0010LHC11h, 33, markerSizeComparison+0.5,  kCyan+1 , kCyan+1);
		graphRatio5TeVChargedPionsPCM0010LHC11h->Draw("E1psame");

//         //Comb to 5TeV
// 		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedPionsComb0010LHC11h, 33, markerSizeComparison+0.5,  kBlue+1 , kBlue+1);
// 		graphRatio5TeVChargedPionsComb0010LHC11h->Draw("E1psame");

        labelPi0CompChargedPionsPbPbLHC11h0010->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray, 2);

		TLegend* legendPi0CompChargedPions5TeVPbPb0010 = new TLegend(0.13,0.77,0.95,0.90);
		legendPi0CompChargedPions5TeVPbPb0010->SetFillColor(0);
		legendPi0CompChargedPions5TeVPbPb0010->SetLineColor(0);
		legendPi0CompChargedPions5TeVPbPb0010->SetNColumns(2);
		legendPi0CompChargedPions5TeVPbPb0010->SetTextSize(0.85*textsizeLabels1);
		legendPi0CompChargedPions5TeVPbPb0010->AddEntry(graphRatioFullPtChargedPionsPCM0010LHC11h,"#pi^{0}/#pi^{#pm} PCM to 2.76TeV","p");
		legendPi0CompChargedPions5TeVPbPb0010->AddEntry(graphRatio5TeVChargedPionsPCM0010LHC11h,"#pi^{0}/#pi^{#pm} PCM to 5TeV","p");
		legendPi0CompChargedPions5TeVPbPb0010->AddEntry(graphRatioFullPtChargedPionsEMCal0010LHC11h,"#pi^{0}/#pi^{#pm} EMCal to 2.76TeV","p");
		legendPi0CompChargedPions5TeVPbPb0010->AddEntry(graphRatio5TeVChargedPionsEMCal0010LHC11h,"#pi^{0}/#pi^{#pm} EMCal to 5TeV","p");
		legendPi0CompChargedPions5TeVPbPb0010->Draw();


	histo2DCompCombinedRatioLHC11h2->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h3->cd();
	pad6PartCompChargedIndPionsLHC11h3->SetLogx();
	histo2DCompCombinedRatioLHC11h->GetXaxis()->SetRangeUser(-0.25,15.);
	histo2DCompCombinedRatioLHC11h->GetYaxis()->SetRangeUser(0.,2.1);
	histo2DCompCombinedRatioLHC11h->DrawCopy();

        //PCM to 2.76TeV
		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsPCM2040LHC11h, 33, markerSizeComparison+0.5, kPink+4,kPink+4);
		graphRatioFullPtChargedPionsPCM2040LHC11h->Draw("E1psame");
        //PCM to 5TeV
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedPionsPCM2050LHC11h, 33, markerSizeComparison+0.5,  kCyan+1 , kCyan+1);
		graphRatio5TeVChargedPionsPCM2050LHC11h->Draw("E1psame");
        //EMCAL to 2.76TeV
		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsEMCal2050LHC11h, markerStylePHOSHighPt, markerSizeComparison, kGreen+2 , kGreen+2);
		graphRatioFullPtChargedPionsEMCal2050LHC11h->Draw("E1psame");
        //EMCAL to 5TeV
		DrawGammaSetMarkerTGraphErr(graphRatio5TeVChargedPionsEMCal2050LHC11h, markerStylePHOSHighPt, markerSizeComparison, colorPHOSHighPt , colorPHOSHighPt);
		graphRatio5TeVChargedPionsEMCal2050LHC11h->Draw("E1psame");

		TLegend* legendPi0CompChargedPions5TeVPbPb2050 = new TLegend(0.16-relativeMarginsX[0],0.77,0.95,0.90);
		legendPi0CompChargedPions5TeVPbPb2050->SetFillColor(0);
		legendPi0CompChargedPions5TeVPbPb2050->SetLineColor(0);
		legendPi0CompChargedPions5TeVPbPb2050->SetNColumns(2);
		legendPi0CompChargedPions5TeVPbPb2050->SetTextSize(0.85*textsizeLabels1);
		legendPi0CompChargedPions5TeVPbPb2050->AddEntry(graphRatioFullPtChargedPionsPCM2040LHC11h,"#pi^{0}/#pi^{#pm} PCM to 2.76TeV","p");
		legendPi0CompChargedPions5TeVPbPb2050->AddEntry(graphRatio5TeVChargedPionsPCM2050LHC11h,"#pi^{0}/#pi^{#pm} PCM to 5TeV","p");
		legendPi0CompChargedPions5TeVPbPb2050->AddEntry(graphRatioFullPtChargedPionsEMCal2050LHC11h,"#pi^{0}/#pi^{#pm} EMCal to 2.76TeV","p");
		legendPi0CompChargedPions5TeVPbPb2050->AddEntry(graphRatio5TeVChargedPionsEMCal2050LHC11h,"#pi^{0}/#pi^{#pm} EMCal to 5TeV","p");
		legendPi0CompChargedPions5TeVPbPb2050->Draw();

		labelPi0CompChargedPionsPbPbLHC11h2050->Draw();
		DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray,2 );

	histo2DCompCombinedRatioLHC11h->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h3->Update();
	canvas6PartCompChargedIndPionsLHC11h->Update();
	canvas6PartCompChargedIndPionsLHC11h->SaveAs(Form("%s/ComparisonCharged5TeVToNeutralPCMeEMCAL_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));



// 	// ***************************************************************************************************************
// 	// ************************************ Comparison pi0/pi+-, pi0 PCM only (2011) *********************************
// 	// ***************************************************************************************************************
	pad6PartCompChargedIndPionsLHC11h1->cd();
	pad6PartCompChargedIndPionsLHC11h1->SetLogx();
	histo2DCompCombinedRatioLHC11h2->GetXaxis()->SetRangeUser(0.,15.);
	histo2DCompCombinedRatioLHC11h2->GetYaxis()->SetRangeUser(0.,2.1);
	histo2DCompCombinedRatioLHC11h2->DrawCopy();

        //2010
        DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM0010, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM0010, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
		graphRatioHighPtChargedPionsPCM0010->Draw("E1psame");
		graphRatioLowPtChargedPionsPCM0010->Draw("E1psame");

        //2011
		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM0010LHC11h, markerStylePCMHighPt, markerSizeComparison, kRed+1, kRed+1);
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM0010LHC11h, markerStylePCMLowPt, markerSizeComparison, kRed, kRed);
		graphRatioHighPtChargedPionsPCM0010LHC11h->Draw("E1psame");
		graphRatioLowPtChargedPionsPCM0010LHC11h->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsPCM0010, 33, markerSizeComparison+0.5, kPink+9 , kPink+9);
		graphRatioFullPtChargedPionsPCM0010->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsPCM0010LHC11h, 33, markerSizeComparison+0.5, kPink+4,kPink+4);
		graphRatioFullPtChargedPionsPCM0010LHC11h->Draw("E1psame");

		labelPi0CompChargedPionsPbPbLHC11h0010->Draw();
		DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);

		TLegend* legendPi0CompChargedOnlyPCMPionsPbPb0010 = new TLegend(0.15,0.68,0.75,0.9);
		legendPi0CompChargedOnlyPCMPionsPbPb0010->SetFillColor(0);
		legendPi0CompChargedOnlyPCMPionsPbPb0010->SetLineColor(0);
		legendPi0CompChargedOnlyPCMPionsPbPb0010->SetNColumns(2);
		legendPi0CompChargedOnlyPCMPionsPbPb0010->SetTextSize(0.85*textsizeLabels1);
//         legendPi0CompChargedOnlyPCMPionsPbPb0010->SetHeader("PCM only");
     legendPi0CompChargedOnlyPCMPionsPbPb0010->AddEntry((TObject*)0,"#pi^{0}/#pi^{#pm} 2010","");
     legendPi0CompChargedOnlyPCMPionsPbPb0010->AddEntry((TObject*)0,"#pi^{0}/#pi^{#pm} 2011","");
        legendPi0CompChargedOnlyPCMPionsPbPb0010->AddEntry(graphRatioLowPtChargedPionsPCM0010,"low #it{p}_{T}","p");
        legendPi0CompChargedOnlyPCMPionsPbPb0010->AddEntry(graphRatioLowPtChargedPionsPCM0010LHC11h,"low #it{p}_{T}","p");
        legendPi0CompChargedOnlyPCMPionsPbPb0010->AddEntry(graphRatioHighPtChargedPionsPCM0010,"high #it{p}_{T}","p");
        legendPi0CompChargedOnlyPCMPionsPbPb0010->AddEntry(graphRatioHighPtChargedPionsPCM0010LHC11h,"high #it{p}_{T}","p");
// 		legendPi0CompChargedOnlyPCMPionsPbPb0010->AddEntry(graphRatioLowPtChargedPionsPCM0005,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (2010)","p");
// 		legendPi0CompChargedOnlyPCMPionsPbPb0010->AddEntry(graphRatioHighPtChargedPionsPCM0005,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (2010)","p");
// 		legendPi0CompChargedOnlyPCMPionsPbPb0010->AddEntry(graphRatioLowPtChargedPionsPCM0010LHC11h,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (2011)","p");
// 		legendPi0CompChargedOnlyPCMPionsPbPb0010->AddEntry(graphRatioHighPtChargedPionsPCM0010LHC11h,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (2011)","p");
// 		legendPi0CompChargedOnlyPCMPionsPbPb0010->AddEntry((TObject*)0,"charged ref: PWGLF-258","");
// 		legendPi0CompChargedOnlyPCMPionsPbPb0010->AddEntry(graphRatioFullPtChargedPionsPCM0010,"#pi^{0}/#pi^{#pm} (PCM 2010)","p");
		legendPi0CompChargedOnlyPCMPionsPbPb0010->AddEntry(graphRatioFullPtChargedPionsPCM0010LHC11h,"#pi^{0}/#pi^{#pm} (PCM 2011)","p");

		legendPi0CompChargedOnlyPCMPionsPbPb0010->Draw();
		DrawGammaLines(0., 15 , 1, 1 ,1,kGray, 2);

        TLatex *thesisLabel = new TLatex(0.8,0.15,"This thesis");
        SetStyleTLatex( thesisLabel,0.9*textsizeLabels2,4);
//          thesisLabel->Draw();


	histo2DCompCombinedRatioLHC11h2->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h1->Update();
	pad6PartCompChargedIndPionsLHC11h3->cd();
	pad6PartCompChargedIndPionsLHC11h3->SetLogx();
	histo2DCompCombinedRatioLHC11h->GetXaxis()->SetRangeUser(0.5,20.);
	histo2DCompCombinedRatioLHC11h->GetYaxis()->SetRangeUser(0.,2.1);
	histo2DCompCombinedRatioLHC11h->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM2040LHC11h, markerStylePCMHighPt, markerSizeComparison, kRed+1, kRed+1);
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM2040LHC11h, markerStylePCMLowPt, markerSizeComparison, kRed, kRed);
		graphRatioHighPtChargedPionsPCM2040LHC11h->Draw("E1psame");
		graphRatioLowPtChargedPionsPCM2040LHC11h->Draw("E1psame");


		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM2040, markerStylePCMHighPt, markerSizeComparison, colorPCMHighPt , colorPCMHighPt);
		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM2040, markerStylePCMLowPt, markerSizeComparison, colorPCMLowPt, colorPCMLowPt);
		graphRatioHighPtChargedPionsPCM2040->Draw("E1psame");
		graphRatioLowPtChargedPionsPCM2040->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsPCM2040, 33, markerSizeComparison+0.5, kPink+9 , kPink+9);
// 		graphRatioFullPtChargedPionsPCM2040->Draw("E1psame");

		labelPi0CompChargedPionsPbPbLHC11h2040->Draw();
		DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);
//         TLatex *thesisLabel2 = new TLatex(0.78,0.15,"This thesis");
//         SetStyleTLatex( thesisLabel2,0.9*textsizeLabels2,4);
//         thesisLabel2->Draw();

	histo2DCompCombinedRatioLHC11h->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h3->Update();
	canvas6PartCompChargedIndPionsLHC11h->Update();
	canvas6PartCompChargedIndPionsLHC11h->SaveAs(Form("%s/ComparisonChargedToNeutralOnlyPCM_splipT_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));


	pad6PartCompChargedIndPionsLHC11h1->cd();
	pad6PartCompChargedIndPionsLHC11h1->SetLogx();
	histo2DCompCombinedRatioLHC11h2->GetXaxis()->SetRangeUser(0.,15.);
	histo2DCompCombinedRatioLHC11h2->GetYaxis()->SetRangeUser(0.,2.1);
	histo2DCompCombinedRatioLHC11h2->DrawCopy();

        //2010
		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsPCM0010, 33, markerSizeComparison+0.5, kPink+9 , kPink+9);
		graphRatioFullPtChargedPionsPCM0010->Draw("E1psame");

        //2011
		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsPCM0010LHC11h, 33, markerSizeComparison+0.5, kPink+4,kPink+4);
		graphRatioFullPtChargedPionsPCM0010LHC11h->Draw("E1psame");

		labelPi0CompChargedPionsPbPbLHC11h0010->Draw();
		DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);

		TLegend* legendPi0CompChargedOnlyPCMPionsPbPb0010_fullPt = new TLegend(0.15,0.68,0.75,0.9);
		legendPi0CompChargedOnlyPCMPionsPbPb0010_fullPt->SetFillColor(0);
		legendPi0CompChargedOnlyPCMPionsPbPb0010_fullPt->SetLineColor(0);
		legendPi0CompChargedOnlyPCMPionsPbPb0010_fullPt->SetNColumns(2);
		legendPi0CompChargedOnlyPCMPionsPbPb0010_fullPt->SetTextSize(0.85*textsizeLabels1);
        legendPi0CompChargedOnlyPCMPionsPbPb0010_fullPt->AddEntry((TObject*)0,"#pi^{0}/#pi^{#pm} 2010","");
        legendPi0CompChargedOnlyPCMPionsPbPb0010_fullPt->AddEntry((TObject*)0,"#pi^{0}/#pi^{#pm} 2011","");
		legendPi0CompChargedOnlyPCMPionsPbPb0010_fullPt->AddEntry(graphRatioFullPtChargedPionsPCM0010,"#pi^{0}/#pi^{#pm} (PCM 2010)","p");
		legendPi0CompChargedOnlyPCMPionsPbPb0010_fullPt->AddEntry(graphRatioFullPtChargedPionsPCM0010LHC11h,"#pi^{0}/#pi^{#pm} (PCM 2011)","p");
		legendPi0CompChargedOnlyPCMPionsPbPb0010_fullPt->Draw();
		DrawGammaLines(0., 15 , 1, 1 ,1,kGray, 2);


	histo2DCompCombinedRatioLHC11h2->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h1->Update();
	pad6PartCompChargedIndPionsLHC11h3->cd();
	pad6PartCompChargedIndPionsLHC11h3->SetLogx();
	histo2DCompCombinedRatioLHC11h->GetXaxis()->SetRangeUser(0.5,20.);
	histo2DCompCombinedRatioLHC11h->GetYaxis()->SetRangeUser(0.,2.1);
	histo2DCompCombinedRatioLHC11h->DrawCopy();

		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsPCM2040, 33, markerSizeComparison+0.5, kPink+9 , kPink+9);
		graphRatioFullPtChargedPionsPCM2040->Draw("E1psame");

		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsPCM2040LHC11h, 33, markerSizeComparison+0.5, kPink+4,kPink+4);
		graphRatioFullPtChargedPionsPCM2040LHC11h->Draw("E1psame");

		labelPi0CompChargedPionsPbPbLHC11h2040->Draw();
		DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);

	histo2DCompCombinedRatioLHC11h->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h3->Update();
	canvas6PartCompChargedIndPionsLHC11h->Update();
	canvas6PartCompChargedIndPionsLHC11h->SaveAs(Form("%s/ComparisonChargedToNeutralOnlyPCM_fullpT_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));


    pad6PartCompChargedIndPionsLHC11h1->cd();
	pad6PartCompChargedIndPionsLHC11h1->SetLogx();
	histo2DCompCombinedRatioLHC11h2->GetXaxis()->SetRangeUser(0.5,20.);
	histo2DCompCombinedRatioLHC11h2->GetYaxis()->SetRangeUser(0.,2.1);
	histo2DCompCombinedRatioLHC11h2->DrawCopy();

// 		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM0010LHC11h, markerStylePCMHighPt, markerSizeComparison, kRed+1, kRed+1);
// 		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM0010LHC11h, markerStylePCMLowPt, markerSizeComparison, kRed, kRed);
// 		graphRatioHighPtChargedPionsPCM0010LHC11h->Draw("E1psame");
// 		graphRatioLowPtChargedPionsPCM0010LHC11h->Draw("E1psame");
		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsPCM0010LHC11h, markerStylePCMHighPt, markerSizeComparison, kRed+1, kRed+1);
		graphRatioFullPtChargedPionsPCM0010LHC11h->Draw("E1psame");

		labelPi0CompChargedPionsPbPbLHC11h0010->Draw();
		DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);

		TLegend* legendPi0CompChargedOnlyPCMPionsPbPb0010forthesis = new TLegend(0.15,0.88-(0.04*1/*2*/),0.5,0.88);
		legendPi0CompChargedOnlyPCMPionsPbPb0010forthesis->SetFillColor(0);
		legendPi0CompChargedOnlyPCMPionsPbPb0010forthesis->SetLineColor(0);
		legendPi0CompChargedOnlyPCMPionsPbPb0010forthesis->SetNColumns(2);
        legendPi0CompChargedOnlyPCMPionsPbPb0010forthesis->SetTextFont(42);
		legendPi0CompChargedOnlyPCMPionsPbPb0010forthesis->SetTextSize(0.85*textsizeLabels1);
        legendPi0CompChargedOnlyPCMPionsPbPb0010forthesis->AddEntry(graphRatioFullPtChargedPionsPCM0010LHC11h,"#pi^{0}/#pi^{#pm}","p");
//         legendPi0CompChargedOnlyPCMPionsPbPb0010forthesis->AddEntry((TObject*)0,"#pi^{0}/#pi^{#pm} (this thesis)","");
//         legendPi0CompChargedOnlyPCMPionsPbPb0010forthesis->AddEntry((TObject*)0,"","");
//         legendPi0CompChargedOnlyPCMPionsPbPb0010forthesis->AddEntry(graphRatioLowPtChargedPionsPCM0010LHC11h,"low #it{p}_{T}","p");
//         legendPi0CompChargedOnlyPCMPionsPbPb0010forthesis->AddEntry(graphRatioHighPtChargedPionsPCM0010LHC11h,"high #it{p}_{T}","p");

// 		legendPi0CompChargedOnlyPCMPionsPbPb0010forthesis->Draw();
		DrawGammaLines(0., 15 , 1, 1 ,1,kGray, 2);

        TLatex *thesisLabel2 = new TLatex(0.75,0.15,"This thesis");
        SetStyleTLatex( thesisLabel2,0.9*textsizeLabels2,4);

	histo2DCompCombinedRatioLHC11h2->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h1->Update();
	pad6PartCompChargedIndPionsLHC11h3->cd();
	pad6PartCompChargedIndPionsLHC11h3->SetLogx();
	histo2DCompCombinedRatioLHC11h->GetXaxis()->SetRangeUser(0.5,20.);
	histo2DCompCombinedRatioLHC11h->GetYaxis()->SetRangeUser(0.,2.1);
	histo2DCompCombinedRatioLHC11h->DrawCopy();

// 		DrawGammaSetMarkerTGraphErr(graphRatioHighPtChargedPionsPCM2040LHC11h, markerStylePCMHighPt, markerSizeComparison, kRed+1, kRed+1);
// 		DrawGammaSetMarkerTGraphErr(graphRatioLowPtChargedPionsPCM2040LHC11h, markerStylePCMLowPt, markerSizeComparison, kRed, kRed);
// 		graphRatioHighPtChargedPionsPCM2040LHC11h->Draw("E1psame");
// 		graphRatioLowPtChargedPionsPCM2040LHC11h->Draw("E1psame");
		DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedPionsPCM2040LHC11h, markerStylePCMHighPt, markerSizeComparison, kRed+1, kRed+1);
		graphRatioFullPtChargedPionsPCM2040LHC11h->Draw("E1psame");

        labelPi0CompChargedPionsPbPbLHC11h2040->Draw();
		DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);
        thesisLabel2->Draw();

	histo2DCompCombinedRatioLHC11h->Draw("axis,same");
	pad6PartCompChargedIndPionsLHC11h3->Update();
	canvas6PartCompChargedIndPionsLHC11h->Update();
	canvas6PartCompChargedIndPionsLHC11h->SaveAs(Form("%s/ComparisonChargedToNeutralOnlyPCMforThesis_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));


    pad6PartCompChargedIndPionsLHC11h1->cd();
    pad6PartCompChargedIndPionsLHC11h1->SetLogx();

	TH2F * histo2DCompCombinedRatioLHC11hEta2 = new TH2F("histo2DCompCombinedRatioLHC11hEta2","histo2DCompCombinedRatioLHC11hEta2",1000,0.3,40.,1000,0.,4.	);
	SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatioLHC11hEta2, "#it{p}_{T} (GeV/#it{c}) ","#eta/K^{#pm}",0.85*textsizeLabels1, textsizeLabels1,
								  0.85*textsizeLabels1, textsizeLabels1, 0.9,0.25/(textsizeFac1*margin), 512, 505);
	histo2DCompCombinedRatioLHC11hEta2->GetYaxis()->SetRangeUser(0.,2.1);
	histo2DCompCombinedRatioLHC11hEta2->GetXaxis()->SetLabelOffset(-0.0105);

	TH2F* histo2DCompCombinedRatioLHC11hEta = new TH2F("histo2DCompCombinedRatioLHC11hEta","histo2DCompCombinedRatioLHC11hEta",1000,0.3,40.,1000,0.,4.	);
	SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatioLHC11hEta, "#it{p}_{T} (GeV/#it{c}) ","#eta/K^{#pm}", 0.85*textsizeLabels2, textsizeLabels2,
								  0.85*textsizeLabels2, textsizeLabels2, 0.9,0.25/(textsizeFac2*margin), 512, 505);
	histo2DCompCombinedRatioLHC11hEta->GetXaxis()->SetLabelOffset(-0.0105);
	histo2DCompCombinedRatioLHC11hEta->GetYaxis()->SetRangeUser(0.6,2.1);
	histo2DCompCombinedRatioLHC11hEta->GetXaxis()->SetRangeUser(.0,20.);
    histo2DCompCombinedRatioLHC11hEta2->GetXaxis()->SetRangeUser(0.5,15.);
    histo2DCompCombinedRatioLHC11hEta2->GetYaxis()->SetRangeUser(0.,2.2);
    histo2DCompCombinedRatioLHC11hEta2->DrawCopy();


        DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedKaonsPCM0010LHC11h, markerStylePCMHighPt, markerSizeComparison, kRed+1, kRed+1);
        graphRatioFullPtChargedKaonsPCM0010LHC11h->Draw("E1psame");
        labelPi0CompChargedPionsPbPbLHC11h0010->Draw();

        TLegend* legendEtaCompChargedOnlyPCMPionsPbPb0010 = new TLegend(0.15,0.88-(0.04*1),0.5,0.88);
		legendEtaCompChargedOnlyPCMPionsPbPb0010->SetFillColor(0);
		legendEtaCompChargedOnlyPCMPionsPbPb0010->SetLineColor(0);
// 		legendEtaCompChargedOnlyPCMPionsPbPb0010->SetNColumns(2);
		legendEtaCompChargedOnlyPCMPionsPbPb0010->SetTextSize(0.85*textsizeLabels1);
        legendEtaCompChargedOnlyPCMPionsPbPb0010->AddEntry(graphRatioFullPtChargedKaonsPCM0010LHC11h,"#eta/K^{#pm}","p");
//         legendEtaCompChargedOnlyPCMPionsPbPb0010->Draw();
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray, 2);


    histo2DCompCombinedRatioLHC11hEta2->Draw("axis,same");
    pad6PartCompChargedIndPionsLHC11h1->Update();
    pad6PartCompChargedIndPionsLHC11h3->cd();
    pad6PartCompChargedIndPionsLHC11h3->SetLogx();
    histo2DCompCombinedRatioLHC11hEta->GetXaxis()->SetRangeUser(0.5,15.);
    histo2DCompCombinedRatioLHC11hEta->GetYaxis()->SetRangeUser(0.,2.2);
    histo2DCompCombinedRatioLHC11hEta->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioFullPtChargedKaonsPCM2040LHC11h, markerStylePCMHighPt, markerSizeComparison, kRed+1, kRed+1);
        graphRatioFullPtChargedKaonsPCM2040LHC11h->Draw("E1psame");

        labelPi0CompChargedPionsPbPbLHC11h2040->Draw();
        DrawGammaLines(0., 15 , 1, 1 ,1,kGray,2);
        thesisLabel2->Draw();

    histo2DCompCombinedRatioLHC11hEta->Draw("axis,same");
    pad6PartCompChargedIndPionsLHC11h3->Update();
    canvas6PartCompChargedIndPionsLHC11h->Update();
    canvas6PartCompChargedIndPionsLHC11h->SaveAs(Form("%s/ComparisonChargedToEtaOnlyPCM_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    delete pad6PartCompChargedIndPionsLHC11h1;
    delete pad6PartCompChargedIndPionsLHC11h3;
    delete canvas6PartCompChargedIndPionsLHC11h;




}