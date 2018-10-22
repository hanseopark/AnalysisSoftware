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


    TFile* filePHENIXPhoton_dAu = new TFile("ExternalInputPbPb/OtherExperiments/PHENIXHeavyIonData.root");

    TGraphErrors* graphPHENIXdAu200GeVGamma_Stat        = (TGraphErrors*)filePHENIXPhoton_dAu->Get("tEdAu200IntConv00100");
    TGraphErrors* graphPHENIXdAu200GeVGamma_Sys         = (TGraphErrors*)filePHENIXPhoton_dAu->Get("tSdAu200IntConv00100");
    TGraphErrors* graphPHENIXdAu200GeVGamma_Tot         = (TGraphErrors*)filePHENIXPhoton_dAu->Get("tSUMdAu200IntConv00100");

    TFile* filePHENIXPhoton_pAu = new TFile("ExternalInputpPb/OtherExperiments/PHENIXSmallSysData.root");

    TGraphErrors* graphPHENIXpAu200GeVGamma_StatLoad        = (TGraphErrors*)filePHENIXPhoton_pAu->Get("tEYieldPAu");
    TGraphErrors* graphPHENIXpAu200GeVGamma_SysLoad         = (TGraphErrors*)filePHENIXPhoton_pAu->Get("tSYieldPAu");
    TGraphErrors* graphPHENIXpAu200GeVGamma_Stat        =  ScaleGraph(graphPHENIXpAu200GeVGamma_StatLoad,1);
    TGraphErrors* graphPHENIXpAu200GeVGamma_Sys         = ScaleGraph(graphPHENIXpAu200GeVGamma_SysLoad,1);
    // TGraphErrors* graphPHENIXpAu200GeVGamma_Stat        = (TGraphErrors*)graphPHENIXpAu200GeVGamma_StatLoad->Clone("tEYieldPAu2");
    // TGraphErrors* graphPHENIXpAu200GeVGamma_Sys         = (TGraphErrors*)graphPHENIXpAu200GeVGamma_SysLoad->Clone("tSYieldPAu2");
    TGraphErrors* graphPHENIXpAu200GeVGamma0020_StatLoad    = (TGraphErrors*)filePHENIXPhoton_pAu->Get("tEYieldPAuCent");
    TGraphErrors* graphPHENIXpAu200GeVGamma0020_SysLoad     = (TGraphErrors*)filePHENIXPhoton_pAu->Get("tSYieldPAuCent");
    TGraphErrors* graphPHENIXpAu200GeVGamma0020_Stat    = ScaleGraph(graphPHENIXpAu200GeVGamma0020_StatLoad,1);
    TGraphErrors* graphPHENIXpAu200GeVGamma0020_Sys     = ScaleGraph(graphPHENIXpAu200GeVGamma0020_SysLoad,1);
    // TGraphErrors* graphPHENIXpAu200GeVGamma0020_Stat    = (TGraphErrors*)graphPHENIXpAu200GeVGamma0020_StatLoad->Clone("tEYieldPAuCent2");
    // TGraphErrors* graphPHENIXpAu200GeVGamma0020_Sys     = (TGraphErrors*)graphPHENIXpAu200GeVGamma0020_SysLoad->Clone("tSYieldPAuCent2");





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
        if (graphPHENIXpAu200GeVGamma0020_Stat) graphPHENIXpAu200GeVGamma0020_Stat->Write("graph_InvYieldDirGamma_PHENIX_pAu_200GeV_Stat_0020");
        if (graphPHENIXpAu200GeVGamma0020_Sys)  graphPHENIXpAu200GeVGamma0020_Sys->Write("graph_InvYieldDirGamma_PHENIX_pAu_200GeV_Sys_0020");
    fileExperimetalSummary.Close();
}
