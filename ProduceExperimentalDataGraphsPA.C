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
#include "CommonHeaders/ConversionFunctions.h"

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;

void ProduceExperimentalDataGraphsPA(){

    // *********************************************************************************************
    // PHENIX published dAu
    // *********************************************************************************************
    TFile* filePHENIXPhoton_dAu = new TFile("ExternalInputPbPb/OtherExperiments/PHENIXHeavyIonData.root");

    TGraphErrors* graphPHENIXdAu200GeVGamma_Stat        = (TGraphErrors*)filePHENIXPhoton_dAu->Get("tEdAu200IntConv00100");
    TGraphErrors* graphPHENIXdAu200GeVGamma_Sys         = (TGraphErrors*)filePHENIXPhoton_dAu->Get("tSdAu200IntConv00100");
    TGraphErrors* graphPHENIXdAu200GeVGamma_Tot         = (TGraphErrors*)filePHENIXPhoton_dAu->Get("tSUMdAu200IntConv00100");
    graphPHENIXdAu200GeVGamma_Stat                      =  ScaleGraph(graphPHENIXdAu200GeVGamma_Stat,1./(5624.336));
    graphPHENIXdAu200GeVGamma_Sys                       = ScaleGraph(graphPHENIXdAu200GeVGamma_Sys,1./(5624.336));
    graphPHENIXdAu200GeVGamma_Tot                       = ScaleGraph(graphPHENIXdAu200GeVGamma_Tot,1./(5624.336));


    // *********************************************************************************************
    // PHENIX preliminary pAu
    // *********************************************************************************************
    TFile* filePHENIXPhoton_pAu = new TFile("ExternalInputpPb/OtherExperiments/PHENIXSmallSysData.root");
    TGraphErrors* graphPHENIXpAu200GeVGamma_StatLoad        = (TGraphErrors*)filePHENIXPhoton_pAu->Get("tEYieldPAu");
    TGraphErrors* graphPHENIXpAu200GeVGamma_SysLoad         = (TGraphErrors*)filePHENIXPhoton_pAu->Get("tSYieldPAu");
    TGraphErrors* graphPHENIXpAu200GeVGamma_Stat            =  ScaleGraph(graphPHENIXpAu200GeVGamma_StatLoad,1./(4520.0144));
    TGraphErrors* graphPHENIXpAu200GeVGamma_Sys             = ScaleGraph(graphPHENIXpAu200GeVGamma_SysLoad,1./(4520.0144));
    TGraphErrors* graphPHENIXpAu200GeVGamma_Tot             = NULL;
    if (graphPHENIXpAu200GeVGamma_Sys && graphPHENIXpAu200GeVGamma_Stat)
        graphPHENIXpAu200GeVGamma_Tot = AddErrorsQuadraticallyTGraph(graphPHENIXpAu200GeVGamma_Stat,graphPHENIXpAu200GeVGamma_Sys);

    // TGraphErrors* graphPHENIXpAu200GeVGamma_Stat         = (TGraphErrors*)graphPHENIXpAu200GeVGamma_StatLoad->Clone("tEYieldPAu2");
    // TGraphErrors* graphPHENIXpAu200GeVGamma_Sys          = (TGraphErrors*)graphPHENIXpAu200GeVGamma_SysLoad->Clone("tSYieldPAu2");
    TGraphErrors* graphPHENIXpAu200GeVGamma0020_StatLoad    = (TGraphErrors*)filePHENIXPhoton_pAu->Get("tEYieldPAuCent");
    TGraphErrors* graphPHENIXpAu200GeVGamma0020_SysLoad     = (TGraphErrors*)filePHENIXPhoton_pAu->Get("tSYieldPAuCent");
    TGraphErrors* graphPHENIXpAu200GeVGamma0020_Stat        = ScaleGraph(graphPHENIXpAu200GeVGamma0020_StatLoad,1./(4520.0144));
    TGraphErrors* graphPHENIXpAu200GeVGamma0020_Sys         = ScaleGraph(graphPHENIXpAu200GeVGamma0020_SysLoad,1./(4520.0144));
    TGraphErrors* graphPHENIXpAu200GeVGamma0020_Tot         = NULL;
    if (graphPHENIXpAu200GeVGamma0020_Sys && graphPHENIXpAu200GeVGamma0020_Stat)
        graphPHENIXpAu200GeVGamma0020_Tot = AddErrorsQuadraticallyTGraph(graphPHENIXpAu200GeVGamma0020_Stat,graphPHENIXpAu200GeVGamma0020_Sys);

    // *********************************************************************************************
    // ALICE work in progress results iso gamma
    // *********************************************************************************************
    TFile* fileALICEIso_pPb5TeV         = new TFile ("ExternalInputpPb/EMCAL/isolated_photon_invariant_yield_20181022.root");
    TGraphErrors* graphALICEpPb5TeVIsoGamma_Stat = (TGraphErrors*)fileALICEIso_pPb5TeV->Get("invariant_yield_EMCNtrChg_EtaBand");

    // *********************************************************************************************
    // ATLAS preliminary results iso gamma
    // *********************************************************************************************
    Int_t npoints_ATLAS_dirgamma_pPb8TeV                    = 17;
    Double_t pT_ATLAS_dirgamma_pPb8TeV[17]                  = { 30, 40, 50, 60, 70, 80, 95, 115, 137.5, 162.5, 187.5, 225, 275, 325, 375, 435, 510};
    // Double_t pTErr_ATLAS_dirgamma_pPb8TeV[17]            = { 5, 5, 5, 5, 5, 5, 10, 10, 12.5, 12.5, 12.5, 25, 25, 25, 25, 35, 40};
    Double_t pTErr_ATLAS_dirgamma_pPb8TeV[17]               = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    Double_t xsection_ATLAS_dirgamma_pPb8TeV[17]            = { 157610, 44823.2, 17190.3, 7723.29, 3746.62, 2025.61, 957.837, 392.197, 158.167, 73.0265, 39.4647, 13.7071, 4.51454, 1.92182, 0.695431, 0.102295, 0.0921991};
    Double_t xsectionErr_ATLAS_dirgamma_pPb8TeV[17]         = { 1836.26, 457.206, 94.5415, 64.3389, 45.6204, 33.7393, 16.9768, 11.108, 6.29703, 4.32938, 3.22617, 1.35363, 0.797988, 0.564711, 0.528743, 0.358239, 0.353141};

    TGraphErrors* graphYieldDirGammapPb8TeV8000GeVATLASTot      = new TGraphErrors( npoints_ATLAS_dirgamma_pPb8TeV, pT_ATLAS_dirgamma_pPb8TeV, xsection_ATLAS_dirgamma_pPb8TeV,
                                                                                    pTErr_ATLAS_dirgamma_pPb8TeV, xsectionErr_ATLAS_dirgamma_pPb8TeV);
    if (graphYieldDirGammapPb8TeV8000GeVATLASTot)
        graphYieldDirGammapPb8TeV8000GeVATLASTot                = ScaleGraph(graphYieldDirGammapPb8TeV8000GeVATLASTot, 1./(70*1e12*2*TMath::Pi())); //70mb=179.7733 GeV/c2, Ncoll pPb8Tev = 7.09

    TFile fileExperimetalSummary("ExternalInputpPb/OtherExperiments/DataCompilationFromOtherEnergiesPA.root","RECREATE");
    fileExperimetalSummary.mkdir("Pi0");
    fileExperimetalSummary.cd("Pi0");

    // NOTE: Meson spectra to be added

    fileExperimetalSummary.mkdir("Eta");
    fileExperimetalSummary.cd("Eta");

    // NOTE: Meson spectra to be added

    fileExperimetalSummary.mkdir("Gamma");
    fileExperimetalSummary.cd("Gamma");
        if (graphPHENIXdAu200GeVGamma_Stat) graphPHENIXdAu200GeVGamma_Stat->Write("graph_InvYieldDirGamma_PHENIX_dAu_200GeV_Stat_0100");
        if (graphPHENIXdAu200GeVGamma_Sys)  graphPHENIXdAu200GeVGamma_Sys->Write("graph_InvYieldDirGamma_PHENIX_dAu_200GeV_Sys_0100");
        if (graphPHENIXdAu200GeVGamma_Tot)  graphPHENIXdAu200GeVGamma_Tot->Write("graph_InvYieldDirGamma_PHENIX_dAu_200GeV_Tot_0100");

        if (graphPHENIXpAu200GeVGamma_Stat) graphPHENIXpAu200GeVGamma_Stat->Write("graph_InvYieldDirGamma_PHENIX_pAu_200GeV_Stat_0100");
        if (graphPHENIXpAu200GeVGamma_Sys)  graphPHENIXpAu200GeVGamma_Sys->Write("graph_InvYieldDirGamma_PHENIX_pAu_200GeV_Sys_0100");
        if (graphPHENIXpAu200GeVGamma_Tot) graphPHENIXpAu200GeVGamma0020_Tot->Write("graph_InvYieldDirGamma_PHENIX_pAu_200GeV_Tot_0100");
        if (graphPHENIXpAu200GeVGamma0020_Stat) graphPHENIXpAu200GeVGamma0020_Stat->Write("graph_InvYieldDirGamma_PHENIX_pAu_200GeV_Stat_0020");
        if (graphPHENIXpAu200GeVGamma0020_Sys)  graphPHENIXpAu200GeVGamma0020_Sys->Write("graph_InvYieldDirGamma_PHENIX_pAu_200GeV_Sys_0020");
        if (graphPHENIXpAu200GeVGamma0020_Tot) graphPHENIXpAu200GeVGamma0020_Tot->Write("graph_InvYieldDirGamma_PHENIX_pAu_200GeV_Tot_0020");
        if (graphALICEpPb5TeVIsoGamma_Stat) graphALICEpPb5TeVIsoGamma_Stat->Write("graph_InvYieldIsoGamma_ALICE_pPb_5TeV_Stat_00100");
        if (graphYieldDirGammapPb8TeV8000GeVATLASTot) graphYieldDirGammapPb8TeV8000GeVATLASTot->Write("graph_InvYieldIsoGamma_ATLAS_pPb_8TeV_Tot_00100");

    fileExperimetalSummary.Close();
}
