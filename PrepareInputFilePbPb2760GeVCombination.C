/****************************************************************************************************************************
********    provided by Gamma Conversion Group, PWGGA,                                                                  *****
********    Lucia Leardini, leardini@cern.ch                                                                            *****
*****************************************************************************************************************************/
#include <Riostream.h>
#include "TMath.h"
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <TROOT.h>
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
#include "TObject.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
// #include "CommonHeaders/FittingGammaConversion.h"
// #include "CommonHeaders/ConversionFunctions.h"
// #include "CommonHeaders/PlottingGammaConversionHistos.h"
// #include "CommonHeaders/PlottingGammaConversionAdditional.h"
// #include "CommonHeaders/CombinationFunctions.h"
#include "PrepareInputFilePbPb2760GeVCombination.h"


void PrepareInputFilePbPb2760GeVCombination(TString fileNamePCM = "data_PCMResults_PbPb_2.76TeV_16Oct2016.root",
                                            TString suffix = "pdf")
{

    TString dateForOutput           = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;
    TString outputDir                       = Form("%s/%s/CombineMesonMeasurementsPbPb2760GeVX",suffix.Data(),dateForOutput.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    //___________________________________ Declaration of files _____________________________________________
    TString fileNamePP2760GeVpublished      = "ExternalInputPbPb/NeutralMesonspp276GeVReference/CombinedResultsPaperPP2760GeV_2017_01_26_FrediV2Clusterizer.root";
    TString fileNamePPPHOS                  = "ExternalInput/PHOS/2.76TeV/LHC11a_PHOS_pi0_pp2760_noBWCorr_FDcorr_20140218.root";
//     TString fileFinalResultsPi0PPforRAA     = "LHC11hInputFiles/CombinedResultsPaperX_18_Feb_2014.root";
//     TString fileFinalResultsEtaPPforRAA     = "LHC11hInputFiles/CombinedResultsEta2760X_2015_09_01.root";
    TString fileNamePCMPublished            = "ExternalInputPbPb/data_PCMResults_PbPb_2.76TeV_LHC10h.root";
//     "FinalResults/CombinedResultsPbPb_ShiftedX_PaperRAA_13_Aug_2014.root";
    TString fileNamePHOS                    = "ExternalInputPbPb/PHOS/LHC10h_PHOS_pi0_PbPb_06022014.root";
    TString fileNameEMCal                   = "ExternalInputPbPb/EMCAL/LHC11h_EMCal_pi0eta_PbPb_10092015.root";
    TString fileNameChargedSpectra          = "ExternalInputPbPb/IdentifiedCharged/JIRA_PWGLF-258/PbPb276.fullpT.INEL.20140329.root";
    TString fileNameChargedRatios           = "ExternalInputPbPb/IdentifiedCharged/JIRA_PWGLF-258/PbPb276.fullpT.RATIOS.20140329.root";
    TString fileNameChargedPionRAA          = "ExternalInputPbPb/IdentifiedCharged/JIRA_PWGLF-258/RAA_Pion_08052014.root";
    TString fileNameChargedKaonRAA          = "ExternalInputPbPb/IdentifiedCharged/JIRA_PWGLF-258/RAA_Kaon_08052014.root";



    TString rootFiles                       = Form("%s/rootFiles",outputDir.Data());
    gSystem->Exec("mkdir -p "+rootFiles);
    gSystem->Exec(Form("cp %s %s/InputPP.root",         fileNamePP2760GeVpublished.Data(), rootFiles.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCM.root",        fileNamePCM.Data(), rootFiles.Data()));
    gSystem->Exec(Form("cp %s %s/InputPHOS.root",       fileNamePHOS.Data(), rootFiles.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCalFull.root",  fileNameEMCal.Data(), rootFiles.Data()));
    gSystem->Exec(Form("cp %s %s/ChargedSpectra.root",  fileNameChargedSpectra.Data(), rootFiles.Data()));
    gSystem->Exec(Form("cp %s %s/ChargedRatios.root",   fileNameChargedRatios.Data(), rootFiles.Data()));
    gSystem->Exec(Form("cp %s %s/ChargedPionRAA.root",  fileNameChargedPionRAA.Data(), rootFiles.Data()));
    gSystem->Exec(Form("cp %s %s/ChargedKaonRAA.root",  fileNameChargedKaonRAA.Data(), rootFiles.Data()));


    //*********************************************************************************************************//
    //**********************************     Charged pions and kaons     **************************************//
    //*********************************************************************************************************//
    TFile* fileDataALICEChargedSpectra = new TFile(fileNameChargedSpectra);
    TH1D* histoChargedPionSpectraStat0005   = (TH1D*)fileDataALICEChargedSpectra->Get("hstat_PbPb276_0005_pion_sum");
    TH1D* histoChargedPionSpectraSyst0005   = (TH1D*)fileDataALICEChargedSpectra->Get("hsys_PbPb276_0005_pion_sum");
    histoChargedPionSpectraStat0005->Scale(0.5);
    histoChargedPionSpectraSyst0005->Scale(0.5);
    TGraphAsymmErrors* graphChargedPionSpectraStat0005 = new TGraphAsymmErrors(histoChargedPionSpectraStat0005);
    TGraphAsymmErrors* graphChargedPionSpectraSyst0005 = new TGraphAsymmErrors(histoChargedPionSpectraSyst0005);
    TH1D* histoChargedPionSpectraStat0510   = (TH1D*)fileDataALICEChargedSpectra->Get("hstat_PbPb276_0510_pion_sum");
    TH1D* histoChargedPionSpectraSyst0510   = (TH1D*)fileDataALICEChargedSpectra->Get("hsys_PbPb276_0510_pion_sum");
    histoChargedPionSpectraStat0510->Scale(0.5);
    histoChargedPionSpectraSyst0510->Scale(0.5);
    TGraphAsymmErrors* graphChargedPionSpectraStat0510 = new TGraphAsymmErrors(histoChargedPionSpectraStat0510);
    TGraphAsymmErrors* graphChargedPionSpectraSyst0510 = new TGraphAsymmErrors(histoChargedPionSpectraSyst0510);
    TH1D* histoChargedPionSpectraStat0010 = (TH1D*)histoChargedPionSpectraStat0510->Clone("histoChargedPionSpectraStat0010");
    TH1D* histoChargedPionSpectraSyst0010 = (TH1D*)histoChargedPionSpectraStat0510->Clone("histoChargedPionSpectraSyst0010");
    histoChargedPionSpectraStat0010->Add(histoChargedPionSpectraStat0005);
    histoChargedPionSpectraSyst0010->Add(histoChargedPionSpectraSyst0005);
    histoChargedPionSpectraStat0010->Scale(0.5);
    histoChargedPionSpectraSyst0010->Scale(0.5);
    for (Int_t i = 1; i < histoChargedPionSpectraSyst0010->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0;
        if (histoChargedPionSpectraSyst0005->GetBinContent(i) != 0){
            relErrLowerCent= histoChargedPionSpectraSyst0005->GetBinError(i)/histoChargedPionSpectraSyst0005->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent = 0;
        if (histoChargedPionSpectraSyst0510->GetBinContent(i) != 0){
            relErrHigherCent = histoChargedPionSpectraSyst0510->GetBinError(i)/histoChargedPionSpectraSyst0510->GetBinContent(i)*100 ;
        }

        if (relErrHigherCent > relErrLowerCent){
            histoChargedPionSpectraSyst0010->SetBinError(i, histoChargedPionSpectraSyst0010->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoChargedPionSpectraSyst0010->SetBinError(i, histoChargedPionSpectraSyst0010->GetBinContent(i)*relErrLowerCent/100);
        }
    }
    TGraphAsymmErrors* graphChargedPionSpectraStat0010 = new TGraphAsymmErrors(histoChargedPionSpectraStat0010);
    TGraphAsymmErrors* graphChargedPionSpectraSyst0010 = new TGraphAsymmErrors(histoChargedPionSpectraSyst0010);
    TH1D* histoChargedPionSpectraStat2040   = (TH1D*)fileDataALICEChargedSpectra->Get("hstat_PbPb276_2040_pion_sum");
    TH1D* histoChargedPionSpectraSyst2040   = (TH1D*)fileDataALICEChargedSpectra->Get("hsys_PbPb276_2040_pion_sum");
    histoChargedPionSpectraStat2040->Scale(0.5);
    histoChargedPionSpectraSyst2040->Scale(0.5);
    TGraphAsymmErrors* graphChargedPionSpectraStat2040 = new TGraphAsymmErrors(histoChargedPionSpectraStat2040);
    TGraphAsymmErrors* graphChargedPionSpectraSyst2040 = new TGraphAsymmErrors(histoChargedPionSpectraSyst2040);

    TH1D* histoChargedKaonSpectraStat0005   = (TH1D*)fileDataALICEChargedSpectra->Get("hstat_PbPb276_0005_kaon_sum");
    TH1D* histoChargedKaonSpectraSyst0005   = (TH1D*)fileDataALICEChargedSpectra->Get("hsys_PbPb276_0005_kaon_sum");
    histoChargedKaonSpectraStat0005->Scale(0.5);
    histoChargedKaonSpectraSyst0005->Scale(0.5);
    TGraphAsymmErrors* graphChargedKaonSpectraStat0005 = new TGraphAsymmErrors(histoChargedKaonSpectraStat0005);
    TGraphAsymmErrors* graphChargedKaonSpectraSyst0005 = new TGraphAsymmErrors(histoChargedKaonSpectraSyst0005);
    TH1D* histoChargedKaonSpectraStat0510   = (TH1D*)fileDataALICEChargedSpectra->Get("hstat_PbPb276_0510_kaon_sum");
    TH1D* histoChargedKaonSpectraSyst0510   = (TH1D*)fileDataALICEChargedSpectra->Get("hsys_PbPb276_0510_kaon_sum");
    histoChargedKaonSpectraStat0510->Scale(0.5);
    histoChargedKaonSpectraSyst0510->Scale(0.5);
    TGraphAsymmErrors* graphChargedKaonSpectraStat0510 = new TGraphAsymmErrors(histoChargedKaonSpectraStat0510);
    TGraphAsymmErrors* graphChargedKaonSpectraSyst0510 = new TGraphAsymmErrors(histoChargedKaonSpectraSyst0510);
    TH1D* histoChargedKaonSpectraStat0010 = (TH1D*)histoChargedKaonSpectraStat0510->Clone("histoChargedKaonSpectraStat0010");
    TH1D* histoChargedKaonSpectraSyst0010 = (TH1D*)histoChargedKaonSpectraStat0510->Clone("histoChargedKaonSpectraSyst0010");
    histoChargedKaonSpectraStat0010->Add(histoChargedKaonSpectraStat0005);
    histoChargedKaonSpectraSyst0010->Add(histoChargedKaonSpectraSyst0005);
    histoChargedKaonSpectraStat0010->Scale(0.5);
    histoChargedKaonSpectraSyst0010->Scale(0.5);
    for (Int_t i = 1; i < histoChargedKaonSpectraSyst0010->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0;
        if (histoChargedKaonSpectraSyst0005->GetBinContent(i) != 0){
            relErrLowerCent= histoChargedKaonSpectraSyst0005->GetBinError(i)/histoChargedKaonSpectraSyst0005->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent = 0;
        if (histoChargedKaonSpectraSyst0510->GetBinContent(i) != 0){
            relErrHigherCent = histoChargedKaonSpectraSyst0510->GetBinError(i)/histoChargedKaonSpectraSyst0510->GetBinContent(i)*100 ;
        }

        if (relErrHigherCent > relErrLowerCent){
            histoChargedKaonSpectraSyst0010->SetBinError(i, histoChargedKaonSpectraSyst0010->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoChargedKaonSpectraSyst0010->SetBinError(i, histoChargedKaonSpectraSyst0010->GetBinContent(i)*relErrLowerCent/100);
        }
    }
    TGraphAsymmErrors* graphChargedKaonSpectraStat0010 = new TGraphAsymmErrors(histoChargedKaonSpectraStat0010);
    TGraphAsymmErrors* graphChargedKaonSpectraSyst0010 = new TGraphAsymmErrors(histoChargedKaonSpectraSyst0010);
    TH1D* histoChargedKaonSpectraStat2040   = (TH1D*)fileDataALICEChargedSpectra->Get("hstat_PbPb276_2040_kaon_sum");
    TH1D* histoChargedKaonSpectraSyst2040   = (TH1D*)fileDataALICEChargedSpectra->Get("hsys_PbPb276_2040_kaon_sum");
    histoChargedKaonSpectraStat2040->Scale(0.5);
    histoChargedKaonSpectraSyst2040->Scale(0.5);
    TGraphAsymmErrors* graphChargedKaonSpectraStat2040 = new TGraphAsymmErrors(histoChargedKaonSpectraStat2040);
    TGraphAsymmErrors* graphChargedKaonSpectraSyst2040 = new TGraphAsymmErrors(histoChargedKaonSpectraSyst2040);

    TFile* fileDataALICEChargedRatioKaonToPion  = new TFile(fileNameChargedRatios);
    TH1D *histoStatChargedRatioKaonToPion0005 = (TH1D*)fileDataALICEChargedRatioKaonToPion->Get("hstat_PbPb276_0005_kaon_to_pion_sum");
    TH1D *histoSysChargedRatioKaonToPion0005 = (TH1D*)fileDataALICEChargedRatioKaonToPion->Get("hsys_PbPb276_0005_kaon_to_pion_sum");
    TGraphAsymmErrors* graphChargedRatioKaonToPion0005 = new TGraphAsymmErrors(histoStatChargedRatioKaonToPion0005);
    TGraphAsymmErrors* graphChargedRatioKaonToPionSys0005 = new TGraphAsymmErrors(histoSysChargedRatioKaonToPion0005);
    TH1D *histoStatChargedRatioKaonToPion0510 = (TH1D*)fileDataALICEChargedRatioKaonToPion->Get("hstat_PbPb276_0510_kaon_to_pion_sum");
    TH1D *histoSysChargedRatioKaonToPion0510 = (TH1D*)fileDataALICEChargedRatioKaonToPion->Get("hsys_PbPb276_0510_kaon_to_pion_sum");
    TGraphAsymmErrors* graphChargedRatioKaonToPion0510 = new TGraphAsymmErrors(histoStatChargedRatioKaonToPion0510);
    TGraphAsymmErrors* graphChargedRatioKaonToPionSys0510 = new TGraphAsymmErrors(histoSysChargedRatioKaonToPion0510);
    TH1D* histoStatChargedRatioKaonToPion0010 = (TH1D*)histoStatChargedRatioKaonToPion0510->Clone("histoStatChargedRatioKaonToPion0010");
    TH1D* histoSysChargedRatioKaonToPion0010 = (TH1D*)histoSysChargedRatioKaonToPion0510->Clone("histoSystChargedRatioKaonToPion0010");
    histoStatChargedRatioKaonToPion0010->Add(histoStatChargedRatioKaonToPion0005);
    histoSysChargedRatioKaonToPion0010->Add(histoSysChargedRatioKaonToPion0005);
    histoStatChargedRatioKaonToPion0010->Scale(0.5);
    histoSysChargedRatioKaonToPion0010->Scale(0.5);
    for (Int_t i = 1; i < histoSysChargedRatioKaonToPion0010->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0;
        if (histoSysChargedRatioKaonToPion0005->GetBinContent(i) != 0){
            relErrLowerCent= histoSysChargedRatioKaonToPion0005->GetBinError(i)/histoSysChargedRatioKaonToPion0005->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent = 0;
        if (histoSysChargedRatioKaonToPion0510->GetBinContent(i) != 0){
            relErrHigherCent = histoSysChargedRatioKaonToPion0510->GetBinError(i)/histoSysChargedRatioKaonToPion0510->GetBinContent(i)*100 ;
        }

        if (relErrHigherCent > relErrLowerCent){
            histoSysChargedRatioKaonToPion0010->SetBinError(i, histoSysChargedRatioKaonToPion0010->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoSysChargedRatioKaonToPion0010->SetBinError(i, histoSysChargedRatioKaonToPion0010->GetBinContent(i)*relErrLowerCent/100);
        }
    }
    TGraphAsymmErrors* graphChargedRatioKaonToPion0010 = new TGraphAsymmErrors(histoStatChargedRatioKaonToPion0010);
    TGraphAsymmErrors* graphChargedRatioKaonToPionSys0010 = new TGraphAsymmErrors(histoSysChargedRatioKaonToPion0010);
    TH1D *histoStatChargedRatioKaonToPion2040 = (TH1D*)fileDataALICEChargedRatioKaonToPion->Get("hstat_PbPb276_2040_kaon_to_pion_sum");
    TH1D *histoSysChargedRatioKaonToPion2040 = (TH1D*)fileDataALICEChargedRatioKaonToPion->Get("hsys_PbPb276_2040_kaon_to_pion_sum");
    TGraphAsymmErrors* graphChargedRatioKaonToPion2040 = new TGraphAsymmErrors(histoStatChargedRatioKaonToPion2040);
    TGraphAsymmErrors* graphChargedRatioKaonToPionSys2040 = new TGraphAsymmErrors(histoSysChargedRatioKaonToPion2040);

    TFile* fileDataALICEChargedPionRAA = new TFile(fileNameChargedPionRAA);
    TH1D *histoStatChargedPion0005 = (TH1D*)fileDataALICEChargedPionRAA->Get("RAAPion_Stat_0_5");
    TH1D *histoSysChargedPion0005 = (TH1D*)fileDataALICEChargedPionRAA->Get("RAAPion_Syst_0_5");
    TGraphAsymmErrors* graphChargedPionRAA0005 = new TGraphAsymmErrors(histoStatChargedPion0005);
    TGraphAsymmErrors* graphChargedPionRAASys0005 = new TGraphAsymmErrors(histoSysChargedPion0005);
    TH1D *histoStatChargedPion0510 = (TH1D*)fileDataALICEChargedPionRAA->Get("RAAPion_Stat_5_10");
    TH1D *histoSysChargedPion0510 = (TH1D*)fileDataALICEChargedPionRAA->Get("RAAPion_Syst_5_10");
    TGraphAsymmErrors* graphChargedPionRAA0510 = new TGraphAsymmErrors(histoStatChargedPion0510);
    TGraphAsymmErrors* graphChargedPionRAASys0510 = new TGraphAsymmErrors(histoSysChargedPion0510);
    TH1D* histoStatChargedPion0010 = (TH1D*)histoStatChargedPion0510->Clone("histoStatChargedPion0010");
    TH1D* histoSysChargedPion0010 = (TH1D*)histoSysChargedPion0510->Clone("histoSystChargedPion0010");
    histoStatChargedPion0010->Add(histoStatChargedPion0005);
    histoSysChargedPion0010->Add(histoSysChargedPion0005);
    histoStatChargedPion0010->Scale(0.5);
    histoSysChargedPion0010->Scale(0.5);
    for (Int_t i = 1; i < histoSysChargedPion0010->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0;
        if (histoSysChargedPion0005->GetBinContent(i) != 0){
            relErrLowerCent= histoSysChargedPion0005->GetBinError(i)/histoSysChargedPion0005->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent = 0;
        if (histoSysChargedPion0510->GetBinContent(i) != 0){
            relErrHigherCent = histoSysChargedPion0510->GetBinError(i)/histoSysChargedPion0510->GetBinContent(i)*100 ;
        }

        if (relErrHigherCent > relErrLowerCent){
            histoSysChargedPion0010->SetBinError(i, histoSysChargedPion0010->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoSysChargedPion0010->SetBinError(i, histoSysChargedPion0010->GetBinContent(i)*relErrLowerCent/100);
        }
    }
    TGraphAsymmErrors* graphChargedPionRAA0010 = new TGraphAsymmErrors(histoStatChargedPion0010);
    TGraphAsymmErrors* graphChargedPionRAASys0010 = new TGraphAsymmErrors(histoSysChargedPion0010);
    TH1D *histoStatChargedPion2040 = (TH1D*)fileDataALICEChargedPionRAA->Get("RAAPion_Stat_20_40");
    TH1D *histoSysChargedPion2040 = (TH1D*)fileDataALICEChargedPionRAA->Get("RAAPion_Syst_20_40");
    TGraphAsymmErrors* graphChargedPionRAA2040 = new TGraphAsymmErrors(histoStatChargedPion2040);
    TGraphAsymmErrors* graphChargedPionRAASys2040 = new TGraphAsymmErrors(histoSysChargedPion2040);

    TFile* fileDataALICEChargedKaonRAA = new TFile(fileNameChargedKaonRAA);
    TH1D *histoStatChargedKaon0005 = (TH1D*)fileDataALICEChargedKaonRAA->Get("RAAKaon_Stat_0_5");
    TH1D *histoSysChargedKaon0005 = (TH1D*)fileDataALICEChargedKaonRAA->Get("RAAKaon_Syst_0_5");
    TGraphAsymmErrors* graphChargedKaonRAA0005 = new TGraphAsymmErrors(histoStatChargedKaon0005);
    TGraphAsymmErrors* graphChargedKaonRAASys0005 = new TGraphAsymmErrors(histoSysChargedKaon0005);
    TH1D *histoStatChargedKaon0510 = (TH1D*)fileDataALICEChargedKaonRAA->Get("RAAKaon_Stat_5_10");
    TH1D *histoSysChargedKaon0510 = (TH1D*)fileDataALICEChargedKaonRAA->Get("RAAKaon_Syst_5_10");
    TGraphAsymmErrors* graphChargedKaonRAA0510 = new TGraphAsymmErrors(histoStatChargedKaon0510);
    TGraphAsymmErrors* graphChargedKaonRAASys0510 = new TGraphAsymmErrors(histoSysChargedKaon0510);
    TH1D* histoStatChargedKaon0010 = (TH1D*)histoStatChargedKaon0510->Clone("histoStatChargedKaon0010");
    TH1D* histoSysChargedKaon0010 = (TH1D*)histoSysChargedKaon0510->Clone("histoSystChargedKaon0010");
    histoStatChargedKaon0010->Add(histoStatChargedKaon0005);
    histoSysChargedKaon0010->Add(histoSysChargedKaon0005);
    histoStatChargedKaon0010->Scale(0.5);
    histoSysChargedKaon0010->Scale(0.5);
    for (Int_t i = 1; i < histoSysChargedKaon0010->GetNbinsX()+1; i++){
        Double_t relErrLowerCent = 0;
        if (histoSysChargedKaon0005->GetBinContent(i) != 0){
            relErrLowerCent= histoSysChargedKaon0005->GetBinError(i)/histoSysChargedKaon0005->GetBinContent(i)*100 ;
        }
        Double_t relErrHigherCent = 0;
        if (histoSysChargedKaon0510->GetBinContent(i) != 0){
            relErrHigherCent = histoSysChargedKaon0510->GetBinError(i)/histoSysChargedKaon0510->GetBinContent(i)*100 ;
        }

        if (relErrHigherCent > relErrLowerCent){
            histoSysChargedKaon0010->SetBinError(i, histoSysChargedKaon0010->GetBinContent(i)*relErrHigherCent/100);
        } else {
            histoSysChargedKaon0010->SetBinError(i, histoSysChargedKaon0010->GetBinContent(i)*relErrLowerCent/100);
        }
    }
    TGraphAsymmErrors* graphChargedKaonRAA0010 = new TGraphAsymmErrors(histoStatChargedKaon0010);
    TGraphAsymmErrors* graphChargedKaonRAASys0010 = new TGraphAsymmErrors(histoSysChargedKaon0010);
    TH1D *histoStatChargedKaon2040 = (TH1D*)fileDataALICEChargedKaonRAA->Get("RAAKaon_Stat_20_40");
    TH1D *histoSysChargedKaon2040 = (TH1D*)fileDataALICEChargedKaonRAA->Get("RAAKaon_Syst_20_40");
    TGraphAsymmErrors* graphChargedKaonRAA2040 = new TGraphAsymmErrors(histoStatChargedKaon2040);
    TGraphAsymmErrors* graphChargedKaonRAASys2040 = new TGraphAsymmErrors(histoSysChargedKaon2040);


    //*********************************************************************************************************//
    //***************************************   PP 2.76TeV input     ******************************************//
    //*********************************************************************************************************//
    fileFinalResultsPP =        new TFile(fileNamePP2760GeVpublished.Data());
    //***************************************   Pi0 Comb     ******************************************//
    TDirectoryFile* directoryPi0PP2760GeV         = (TDirectoryFile*)fileFinalResultsPP->Get("Pi02.76TeV");
    graphInvSectionCombStatPi02760GeV = (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0Comb2760GeVAStatErr");
    graphInvSectionCombStatPi02760GeV = ScaleGraph(graphInvSectionCombStatPi02760GeV,1./xSection2760GeVppINEL);
    graphInvSectionCombSysPi02760GeV = (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0Comb2760GeVASysErr");
    graphInvSectionCombSysPi02760GeV = ScaleGraph(graphInvSectionCombSysPi02760GeV,1./xSection2760GeVppINEL);
    TF1 *fitInvCrossSectionTsallisPi0Comb2760GeV = (TF1*)directoryPi0PP2760GeV->Get("TsallisFitPi0");
    fitInvCrossSectionTsallisPi0Comb2760GeV->SetParameter(0, fitInvCrossSectionTsallisPi0Comb2760GeV->GetParameter(0)*1./xSection2760GeVppINEL);
    TF1 *fitInvCrossSectionTCMPi0Comb2760GeV = (TF1*)directoryPi0PP2760GeV->Get("TwoComponentModelFitPi0");
    fitInvCrossSectionTCMPi0Comb2760GeV->SetParameter(0, fitInvCrossSectionTCMPi0Comb2760GeV->GetParameter(0)*1./xSection2760GeVppINEL);
    fitInvCrossSectionTCMPi0Comb2760GeV->SetParameter(2, fitInvCrossSectionTCMPi0Comb2760GeV->GetParameter(2)*1./xSection2760GeVppINEL);
    graphInvSectionCombStatPi02760GeVforRAA = (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0Comb2760GeVAStatErr_yShifted");
    graphInvSectionCombStatPi02760GeVforRAA = ScaleGraph(graphInvSectionCombStatPi02760GeVforRAA,1./xSection2760GeVppINEL);
    graphInvSectionCombSysPi02760GeVforRAA = (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0Comb2760GeVASysErr_yShifted");
    graphInvSectionCombSysPi02760GeVforRAA = ScaleGraph(graphInvSectionCombSysPi02760GeVforRAA,1./xSection2760GeVppINEL);

    //***************************************   Pi0 PCM     ******************************************//
    graphInvSectionPCMStatPi02760GeV = (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0PCM2760GeVStatErr");
    graphInvSectionPCMStatPi02760GeV = ScaleGraph(graphInvSectionPCMStatPi02760GeV,1./xSection2760GeVppINEL);
    graphInvSectionPCMSysPi02760GeV = (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0PCM2760GeVSysErr");
    graphInvSectionPCMSysPi02760GeV = ScaleGraph(graphInvSectionPCMSysPi02760GeV,1./xSection2760GeVppINEL);
    graphInvSectionPCMStatPi02760GeVforRAA =       (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0PCM2760GeVStatErr_yShifted");
    graphInvSectionPCMStatPi02760GeVforRAA = ScaleGraph(graphInvSectionPCMStatPi02760GeVforRAA,1./xSection2760GeVppINEL);
    //WITH material budget error
    graphInvSectionPCMSysPi02760GeV_yShifted =       (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0PCM2760GeVSysErr_yShifted");
    graphInvSectionPCMSysPi02760GeV_yShifted = ScaleGraph(graphInvSectionPCMSysPi02760GeV_yShifted,1./xSection2760GeVppINEL);
    // subtraction of material budget error:
    Int_t pi0PCMbins = graphInvSectionPCMSysPi02760GeV_yShifted->GetN();
    Double_t *xWM = graphInvSectionPCMSysPi02760GeV_yShifted->GetX();
    Double_t *yWM = graphInvSectionPCMSysPi02760GeV_yShifted->GetY();
    Double_t xErrWM[pi0PCMbins];
    Double_t yErrWM[pi0PCMbins];
    Double_t yrelErrWM[pi0PCMbins];
    Double_t yErrorWOM[pi0PCMbins];
    Double_t yrelErrorWOM[pi0PCMbins];
    for(Int_t i = 0;i<graphInvSectionPCMSysPi02760GeV_yShifted->GetN();i++){

      xErrWM[i] = graphInvSectionPCMSysPi02760GeV_yShifted->GetErrorXlow(i);
      yrelErrWM[i] = graphInvSectionPCMSysPi02760GeV_yShifted->GetErrorYlow(i);
      yErrWM[i] = (yrelErrWM[i] * 100.)/yWM[i];

      yErrorWOM[i] = TMath::Sqrt( (yErrWM[i]*yErrWM[i]) - (9.*9.) );
      yrelErrorWOM[i] = (yWM[i]*yErrorWOM[i])/100.;
//           cout << "yErrWM[i] " << yErrWM[i] << endl;
//           cout << "yErrorWOM[i] " << yErrorWOM[i] << endl;

    }
    //WITHOUT material budget error
    graphInvSectionPCMSysPi02760GeVforRAA = new TGraphAsymmErrors(pi0PCMbins,xWM,yWM,xErrWM,xErrWM,yrelErrorWOM,yrelErrorWOM);

    //***************************************   Pi0 PHOS     ******************************************//
    graphInvSectionPHOSStatPi02760GeV = (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0PHOS2760GeVStatErr");
    graphInvSectionPHOSStatPi02760GeV = ScaleGraph(graphInvSectionPHOSStatPi02760GeV,1./xSection2760GeVppINEL);
    graphInvSectionPHOSSysPi02760GeV = (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0PHOS2760GeVSysErr");
    graphInvSectionPHOSSysPi02760GeV = ScaleGraph(graphInvSectionPHOSSysPi02760GeV,1./xSection2760GeVppINEL);
    graphInvSectionPHOSStatPi02760GeVforRAA =       (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0PHOS2760GeVStatErr_yShifted");
    graphInvSectionPHOSStatPi02760GeVforRAA = ScaleGraph(graphInvSectionPHOSStatPi02760GeVforRAA,1./xSection2760GeVppINEL);
    //WITH common errors
    graphInvSectionPHOSSysPi02760GeV_yShifted =       (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0PHOS2760GeVSysErr_yShifted");
    graphInvSectionPHOSSysPi02760GeV_yShifted = ScaleGraph(graphInvSectionPHOSSysPi02760GeV_yShifted,1./xSection2760GeVppINEL);
    TFile *filePPPHOS =        new TFile(fileNamePPPHOS.Data());
    TDirectoryFile* directoryPHOSPi02760GeV = (TDirectoryFile*)filePPPHOS->Get("pp2760");
    TH1D *histoPi0PhosSysRAA2760GeV =     (TH1D*)directoryPHOSPi02760GeV->Get("hPi02760GeVSysTypeB"); //already divided by xsection
    TGraphAsymmErrors *graphPi0PhosSysRAA2760GeV = new TGraphAsymmErrors(histoPi0PhosSysRAA2760GeV);
      graphPi0PhosSysRAA2760GeV->RemovePoint(0);
    // taking type B error histo
    Int_t pi0PHOSbins = graphInvSectionPHOSSysPi02760GeV_yShifted->GetN();
    Double_t *xsys = graphInvSectionPHOSSysPi02760GeV_yShifted->GetX();
    Double_t *ysys = graphInvSectionPHOSSysPi02760GeV_yShifted->GetY();
    Double_t xsysErr[pi0PHOSbins];
    Double_t ysysErrtypeB[pi0PHOSbins];
    Double_t percentErrorUp;
    Double_t percentErrorLow;
    //WITHOUT common errors (type B)
    graphInvSectionPHOSSysPi02760GeVforRAA = new TGraphAsymmErrors(pi0PHOSbins);
    for(Int_t i = 0;i<graphInvSectionPHOSSysPi02760GeV_yShifted->GetN();i++){
      graphInvSectionPHOSSysPi02760GeVforRAA->SetPoint(i,xsys[i],ysys[i]);
      percentErrorLow = (graphPi0PhosSysRAA2760GeV->GetErrorYlow(i)*100)/graphPi0PhosSysRAA2760GeV->GetY()[i];
      ysysErrtypeB[i] = percentErrorLow*ysys[i]/100;
      graphInvSectionPHOSSysPi02760GeVforRAA->SetPointError(i,graphInvSectionPHOSSysPi02760GeV_yShifted->GetErrorXlow(i),graphInvSectionPHOSSysPi02760GeV_yShifted->GetErrorXhigh(i),ysysErrtypeB[i],ysysErrtypeB[i]);
    }

    //***************************************   Pi0 EMCal     ******************************************//
    graphInvSectionEMCalStatPi02760GeV = (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0EMCAL2760GeVStatErr");
    graphInvSectionEMCalStatPi02760GeV = ScaleGraph(graphInvSectionEMCalStatPi02760GeV,1./xSection2760GeVppINEL);
    graphInvSectionEMCalSysPi02760GeV = (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0EMCAL2760GeVSysErr");
    graphInvSectionEMCalSysPi02760GeV = ScaleGraph(graphInvSectionEMCalSysPi02760GeV,1./xSection2760GeVppINEL);
    graphInvSectionEMCalStatPi02760GeVforRAA =       (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0PCMEMCAL2760GeVStatErr_yShifted");
    graphInvSectionEMCalStatPi02760GeVforRAA = ScaleGraph(graphInvSectionEMCalStatPi02760GeVforRAA,1./xSection2760GeVppINEL);
    //WITH common errors
    graphInvSectionEMCalSysPi02760GeV_yShifted =       (TGraphAsymmErrors*)directoryPi0PP2760GeV->Get("graphInvCrossSectionPi0PCMEMCAL2760GeVSysErr_yShifted");
    graphInvSectionEMCalSysPi02760GeV_yShifted = ScaleGraph(graphInvSectionEMCalSysPi02760GeV_yShifted,1./xSection2760GeVppINEL);
    // subtraction of common errors (material budget and resolution(?)):
    Int_t pi0EMCalbins = graphInvSectionEMCalSysPi02760GeV_yShifted->GetN();
    Double_t *xWCommonErr = graphInvSectionEMCalSysPi02760GeV_yShifted->GetX();
    Double_t *yWCommonErr = graphInvSectionEMCalSysPi02760GeV_yShifted->GetY();
    Double_t xErrWCommonErr[pi0EMCalbins];
    Double_t yErrWCommonErr[pi0EMCalbins];
    Double_t yrelErrWCommonErr[pi0EMCalbins];
    Double_t yErrorWOCommonErr[pi0EMCalbins];
    Double_t yrelErrorWOCommonErr[pi0EMCalbins];
    for(Int_t i = 0;i<graphInvSectionEMCalSysPi02760GeV_yShifted->GetN();i++){
      xErrWCommonErr[i] = graphInvSectionEMCalSysPi02760GeV_yShifted->GetErrorXlow(i);
      yrelErrWCommonErr[i] = graphInvSectionEMCalSysPi02760GeV_yShifted->GetErrorYlow(i);
      yErrWCommonErr[i] = (yrelErrWCommonErr[i] * 100.)/yWCommonErr[i];
      yErrorWOCommonErr[i] = TMath::Sqrt( (yErrWCommonErr[i]*yErrWCommonErr[i]) - (5.2*5.2) ); // 3% (rougly en. scale) and 4.24 (3% for TRD mat, 3% for TOF mat added in quadrature)
      yrelErrorWOCommonErr[i] = (yWCommonErr[i]*yErrorWOCommonErr[i])/100.;
    }
    //WITHOUT common errors
    graphInvSectionEMCalSysPi02760GeVforRAA = new TGraphAsymmErrors(pi0EMCalbins,xWCommonErr,yWCommonErr,xErrWCommonErr,xErrWCommonErr,yrelErrorWOCommonErr,yrelErrorWOCommonErr);


    //***************************************   Eta Comb     ******************************************//
    TDirectoryFile* directoryEtaPP2760GeV         = (TDirectoryFile*)fileFinalResultsPP->Get("Eta2.76TeV");
    graphInvSectionCombStatEta2760GeV = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaComb2760GeVAStatErr");
    graphInvSectionCombStatEta2760GeV = ScaleGraph(graphInvSectionCombStatEta2760GeV,1./xSection2760GeVppINEL);
    graphInvSectionCombSysEta2760GeV = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaComb2760GeVASysErr");
    graphInvSectionCombSysEta2760GeV = ScaleGraph(graphInvSectionCombSysEta2760GeV,1./xSection2760GeVppINEL);
    TF1 *fitInvCrossSectionTsallisEtaComb2760GeV = (TF1*)directoryEtaPP2760GeV->Get("TsallisFitEta");
    fitInvCrossSectionTsallisEtaComb2760GeV->SetParameter(0, fitInvCrossSectionTsallisEtaComb2760GeV->GetParameter(0)*1./xSection2760GeVppINEL);
    TF1 *fitInvCrossSectionTCMEtaComb2760GeV = (TF1*)directoryEtaPP2760GeV->Get("TwoComponentModelFitEta");
    fitInvCrossSectionTCMEtaComb2760GeV->SetParameter(0, fitInvCrossSectionTCMEtaComb2760GeV->GetParameter(0)*1./xSection2760GeVppINEL);
    fitInvCrossSectionTCMEtaComb2760GeV->SetParameter(2, fitInvCrossSectionTCMEtaComb2760GeV->GetParameter(2)*1./xSection2760GeVppINEL);
    TGraphAsymmErrors* graphInvSectionCombStatEta2760GeV_yShifted = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaComb2760GeVAStatErr_yShifted");
    graphInvSectionCombStatEta2760GeV_yShifted = ScaleGraph(graphInvSectionCombStatEta2760GeV_yShifted,1./xSection2760GeVppINEL);
    TGraphAsymmErrors* graphInvSectionCombSysEta2760GeV_yShifted = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaComb2760GeVASysErr_yShifted");
    graphInvSectionCombSysEta2760GeV_yShifted = ScaleGraph(graphInvSectionCombSysEta2760GeV_yShifted,1./xSection2760GeVppINEL);

    graphRatioEtaToPi0Comb2760GeVStatErr = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphRatioEtaToPi0Comb2760GeVStatErr");
    graphRatioEtaToPi0Comb2760GeVSysErr = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphRatioEtaToPi0Comb2760GeVSysErr");

    //***************************************   Eta PCM     ******************************************//
    graphInvSectionPCMStatEta2760GeV = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaPCM2760GeVStatErr");
    graphInvSectionPCMStatEta2760GeV = ScaleGraph(graphInvSectionPCMStatEta2760GeV,1./xSection2760GeVppINEL);
    graphInvSectionPCMSysEta2760GeV = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaPCM2760GeVSysErr");
    graphInvSectionPCMSysEta2760GeV = ScaleGraph(graphInvSectionPCMSysEta2760GeV,1./xSection2760GeVppINEL);
    graphInvSectionPCMStatEta2760GeVforRAA =       (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaPCM2760GeVStatErr_yShifted");
    graphInvSectionPCMStatEta2760GeVforRAA = ScaleGraph(graphInvSectionPCMStatEta2760GeVforRAA,1./xSection2760GeVppINEL);
    //WITH material budget error
    graphInvSectionPCMSysEta2760GeV_yShifted =       (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaPCM2760GeVSysErr_yShifted");
    graphInvSectionPCMSysEta2760GeV_yShifted = ScaleGraph(graphInvSectionPCMSysEta2760GeV_yShifted,1./xSection2760GeVppINEL);
    // subtraction of material budget error:
    Int_t etaPCMbins = graphInvSectionPCMSysEta2760GeV_yShifted->GetN();
    Double_t *xetaWM = graphInvSectionPCMSysEta2760GeV_yShifted->GetX();
    Double_t *yetaWM = graphInvSectionPCMSysEta2760GeV_yShifted->GetY();
    Double_t xetaErrWM[etaPCMbins];
    Double_t yetaErrWM[etaPCMbins];
    Double_t yetarelErrWM[etaPCMbins];
    Double_t yetaErrorWOM[etaPCMbins];
    Double_t yetarelErrorWOM[etaPCMbins];
    for(Int_t i = 0;i<graphInvSectionPCMSysEta2760GeV_yShifted->GetN();i++){
      xetaErrWM[i] = graphInvSectionPCMSysEta2760GeV_yShifted->GetErrorXlow(i);
      yetarelErrWM[i] = graphInvSectionPCMSysEta2760GeV_yShifted->GetErrorYlow(i);
      yetaErrWM[i] = (yetarelErrWM[i] * 100.)/yetaWM[i];
      yetaErrorWOM[i] = TMath::Sqrt( (yetaErrWM[i]*yetaErrWM[i]) - (9.*9.) );
      yetarelErrorWOM[i] = (yetaWM[i]*yetaErrorWOM[i])/100.;
    }
    //WITHOUT material budget error
    graphInvSectionPCMSysEta2760GeVforRAA = new TGraphAsymmErrors(etaPCMbins,xetaWM,yetaWM,xetaErrWM,xetaErrWM,yetarelErrorWOM,yetarelErrorWOM);

    //***************************************   Eta EMCal     ******************************************//
    graphInvSectionEMCalStatEta2760GeV = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaEMCAL2760GeVStatErr");
    graphInvSectionEMCalStatEta2760GeV = ScaleGraph(graphInvSectionEMCalStatEta2760GeV,1./xSection2760GeVppINEL);
    graphInvSectionEMCalSysEta2760GeV = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaEMCAL2760GeVSysErr");
    graphInvSectionEMCalSysEta2760GeV = ScaleGraph(graphInvSectionEMCalSysEta2760GeV,1./xSection2760GeVppINEL);
    graphInvSectionEMCalStatEta2760GeVforRAA =       (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaEMCAL2760GeVStatErr_yShifted");
      graphInvSectionEMCalStatEta2760GeVforRAA = ScaleGraph(graphInvSectionEMCalStatEta2760GeVforRAA,1./xSection2760GeVppINEL);
    //WITH common errors
    graphInvSectionEMCalSysEta2760GeV_yShifted =       (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphInvCrossSectionEtaEMCAL2760GeVSysErr_yShifted");
    graphInvSectionEMCalSysEta2760GeV_yShifted = ScaleGraph(graphInvSectionEMCalSysEta2760GeV_yShifted,1./xSection2760GeVppINEL);
    // subtraction of common errors (material budget and resolution (?)):
    Int_t etaEMCalbins = graphInvSectionEMCalSysEta2760GeV_yShifted->GetN();
    Double_t *xetaWCommonErr = graphInvSectionEMCalSysEta2760GeV_yShifted->GetX();
    Double_t *yetaWCommonErr = graphInvSectionEMCalSysEta2760GeV_yShifted->GetY();
    Double_t xetaErrWCommonErr[etaEMCalbins];
    Double_t yetaErrWCommonErr[etaEMCalbins];
    Double_t yetarelErrWCommonErr[etaEMCalbins];
    Double_t yetaErrorWOCommonErr[etaEMCalbins];
    Double_t yetarelErrorWOCommonErr[etaEMCalbins];
    for(Int_t i = 0;i<graphInvSectionEMCalSysEta2760GeV_yShifted->GetN();i++){
      xetaErrWCommonErr[i] = graphInvSectionEMCalSysEta2760GeV_yShifted->GetErrorXlow(i);
      yetarelErrWCommonErr[i] = graphInvSectionEMCalSysEta2760GeV_yShifted->GetErrorYlow(i);
      yetaErrWCommonErr[i] = (yetarelErrWCommonErr[i] * 100.)/yetaWCommonErr[i];
      yetaErrorWOCommonErr[i] = TMath::Sqrt( (yetaErrWCommonErr[i]*yetaErrWCommonErr[i]) - (6.5*6.5) ); // 5% (rougly en. scale) and 4.24 (3% for TRD mat, 3% for TOF mat added in quadrature)
      yetarelErrorWOCommonErr[i] = (yetaWCommonErr[i]*yetaErrorWOCommonErr[i])/100.;
    }
    //WITHOUT common errors
    graphInvSectionEMCalSysEta2760GeVforRAA = new TGraphAsymmErrors(etaEMCalbins,xetaWCommonErr,yetaWCommonErr,xetaErrWCommonErr,xetaErrWCommonErr,yetarelErrorWOCommonErr,yetarelErrorWOCommonErr);

	TFile* fileCombinedpp 				= new TFile("FinalResults/CombinedResultsPP_ShiftedX_PaperRAA_16_May_2014.root ");
	TGraphAsymmErrors* graphCombEtaToPi0Ratiopp7TeV =         (TGraphAsymmErrors*)fileCombinedpp->Get("graphEtaToPi0Comb7TeVStat");
	TGraphAsymmErrors* graphCombEtaToPi0RatioSysErrpp7TeV=    (TGraphAsymmErrors*)fileCombinedpp->Get("graphEtaToPi0Comb7TeVSys"); 
	TGraphAsymmErrors* graphCombEtaToPi0Ratiopp7TeVNoXErrors = (TGraphAsymmErrors*)graphCombEtaToPi0Ratiopp7TeV->Clone();


    //*********************************************************************************************************//
    //**************************************   PbPb 2.76TeV input     *****************************************//
    //*********************************************************************************************************//

    //***************************************   Pi0 PCM 2010  ******************************************//
    TFile *filePublished    = new TFile(fileNamePCMPublished.Data());
    if(!filePublished) return;
    TDirectoryFile *Pi0Published0010 = (TDirectoryFile*)filePublished->Get("Pi0_PbPb_2.76TeV_0-10%");
    TH1D* histoPCMPublishedInvYieldStatPbPb2760GeV_0010 = (TH1D*)Pi0Published0010->Get("CorrectedYieldPi0");
    TGraphAsymmErrors* graphPCMPublishedInvYieldStatPbPb2760GeV_0010 = new TGraphAsymmErrors(histoPCMPublishedInvYieldStatPbPb2760GeV_0010);
    graphPCMPublishedInvYieldStatPbPb2760GeV_0010->RemovePoint(0);
    graphPCMPublishedInvYieldStatPbPb2760GeV_0010->RemovePoint(0);
    TGraphAsymmErrors* graphPCMPublishedInvYieldSysPbPb2760GeV_0010  = (TGraphAsymmErrors*)Pi0Published0010->Get("Pi0SystError");

    TDirectoryFile *Pi0Published2040 = (TDirectoryFile*)filePublished->Get("Pi0_PbPb_2.76TeV_20-40%");
    TH1D* histoPCMPublishedInvYieldStatPbPb2760GeV_2040 = (TH1D*)Pi0Published2040->Get("CorrectedYieldPi0");
    TGraphAsymmErrors* graphPCMPublishedInvYieldStatPbPb2760GeV_2040     = new TGraphAsymmErrors(histoPCMPublishedInvYieldStatPbPb2760GeV_2040);
    graphPCMPublishedInvYieldStatPbPb2760GeV_2040->RemovePoint(0);
    graphPCMPublishedInvYieldStatPbPb2760GeV_2040->RemovePoint(0);
    TGraphAsymmErrors* graphPCMPublishedInvYieldSysPbPb2760GeV_2040  = (TGraphAsymmErrors*)Pi0Published2040->Get("Pi0SystError");


    //***************************************   Pi0 PHOS     ******************************************//
    TFile *filePHOSPbPb =       new TFile(fileNamePHOS);
    TDirectory *directoryPHOSPi0PbPb0010 = (TDirectory*)filePHOSPbPb->Get("pi0_PbPb_2760_Centrality_0-10%");
    histoPi0PHOSPbPb0010 =      (TH1D*)directoryPHOSPi0PbPb0010->Get("hPi0_PbPb_cen7_NoBW_Stat");
    histoPi0PHOSSysPbPb0010 =   (TH1D*)directoryPHOSPi0PbPb0010->Get("hPi0_PbPb_cen7_NoBW_Syst");
    TGraphAsymmErrors *graphPHOSPi0InvYieldStatPbPb2760GeV_0010 = new TGraphAsymmErrors(histoPi0PHOSPbPb0010);
    TGraphAsymmErrors *graphPHOSPi0InvYieldSysPbPb2760GeV_0010 = new TGraphAsymmErrors(histoPi0PHOSSysPbPb0010);
    histoPi0PHOSSysRAAPbPb0010 =    (TH1D*)directoryPHOSPi0PbPb0010->Get("hPi0_PbPb_cen7_SystRaa"); //bin shifted
    graphSysErrRAAYieldPi0PHOSPbPb0010 = new TGraphAsymmErrors(histoPi0PHOSSysRAAPbPb0010);
    graphPHOSPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(0);
    graphPHOSPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(0);
    graphPHOSPi0InvYieldSysPbPb2760GeV_0010->RemovePoint(0);
    graphPHOSPi0InvYieldSysPbPb2760GeV_0010->RemovePoint(0);
    graphSysErrRAAYieldPi0PHOSPbPb0010->RemovePoint(0);
    graphSysErrRAAYieldPi0PHOSPbPb0010->RemovePoint(0);
    Double_t error1;
    Double_t error2;
    for(Int_t i=0;i<graphPHOSPi0InvYieldStatPbPb2760GeV_0010->GetN();i++){
      graphSysErrRAAYieldPi0PHOSPbPb0010->SetPoint(i,graphPHOSPi0InvYieldStatPbPb2760GeV_0010->GetX()[i],graphPHOSPi0InvYieldStatPbPb2760GeV_0010->GetY()[i]);
      percentErrorLow = (graphSysErrRAAYieldPi0PHOSPbPb0010->GetErrorYlow(i)*100)/graphSysErrRAAYieldPi0PHOSPbPb0010->GetY()[i];
      graphSysErrRAAYieldPi0PHOSPbPb0010->SetPointError(i,graphPHOSPi0InvYieldStatPbPb2760GeV_0010->GetErrorXlow(i),graphPHOSPi0InvYieldStatPbPb2760GeV_0010->GetErrorXhigh(i),(percentErrorLow*graphPHOSPi0InvYieldStatPbPb2760GeV_0010->GetY()[i]/100.),(percentErrorLow*graphPHOSPi0InvYieldStatPbPb2760GeV_0010->GetY()[i]/100.));
    }
    graphPHOSPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(graphPHOSPi0InvYieldStatPbPb2760GeV_0010->GetN()-1);
    graphPHOSPi0InvYieldSysPbPb2760GeV_0010->RemovePoint(graphPHOSPi0InvYieldSysPbPb2760GeV_0010->GetN()-1);
    graphSysErrRAAYieldPi0PHOSPbPb0010->RemovePoint(graphSysErrRAAYieldPi0PHOSPbPb0010->GetN()-1);
    graphPHOSPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(graphPHOSPi0InvYieldStatPbPb2760GeV_0010->GetN()-1);
    graphPHOSPi0InvYieldSysPbPb2760GeV_0010->RemovePoint(graphPHOSPi0InvYieldSysPbPb2760GeV_0010->GetN()-1);
    graphSysErrRAAYieldPi0PHOSPbPb0010->RemovePoint(graphSysErrRAAYieldPi0PHOSPbPb0010->GetN()-1);

    //***************************************   Pi0 PCM   ******************************************//
    TFile* filePCM = new TFile(fileNamePCM.Data());
    TDirectoryFile* directoryPCMPi0PbPb2760GeV_0010         = (TDirectoryFile*)filePCM->Get("Pi0_PbPb_2.76TeV_0-10%");
    TH1D* histoPCMPi0InvYieldPbPb2760GeV_0010                   = (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("CorrectedYieldPi0");
    TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeV_0010  = new TGraphAsymmErrors(histoPCMPi0InvYieldPbPb2760GeV_0010);
    while(graphPCMPi0InvYieldStatPbPb2760GeV_0010->GetY()[0] < 1e-14) graphPCMPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(0);
    TH1D* histoPCMPi0InvYieldPbPb2760GeVYshifted_0010                   = (TH1D*)directoryPCMPi0PbPb2760GeV_0010->Get("CorrectedYieldBinShiftedPi0");
    TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_0010  = new TGraphAsymmErrors(histoPCMPi0InvYieldPbPb2760GeVYshifted_0010);
    while(graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_0010->GetY()[0] < 1e-14) graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_0010->RemovePoint(0);
    //WITH material budget error
    TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_0010       = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0SystError");
    //WITHOUT material budget error
    TGraphAsymmErrors* graphPCMPi0InvYieldSysWOMat2760GeV_0010          = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0SystErrorA");
    TGraphAsymmErrors* graphPCMPi0RAAStatPbPb2760GeV_0010   = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0RAA");
    TH1D* histoPCMPi0RAAStatPbPb2760GeV_0010 = (TH1D*)GraphAsymErrorsToHist(graphPCMPi0RAAStatPbPb2760GeV_0010,14,"histoPCMPi0RAAStatPbPb2760GeV_0010");
    TGraphAsymmErrors* graphPCMPi0RAASysPbPb2760GeV_0010    = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_0010->Get("Pi0RAASys");

    TDirectoryFile* directoryPCMPi0PbPb2760GeV_2040                 = (TDirectoryFile*)filePCM->Get("Pi0_PbPb_2.76TeV_20-40%");
    TH1D* histoPCMPi0InvYieldPbPb2760GeV_2040                   = (TH1D*)directoryPCMPi0PbPb2760GeV_2040->Get("CorrectedYieldPi0");
    TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeV_2040  = new TGraphAsymmErrors(histoPCMPi0InvYieldPbPb2760GeV_2040);
    while(graphPCMPi0InvYieldStatPbPb2760GeV_2040->GetY()[0] < 1e-14) graphPCMPi0InvYieldStatPbPb2760GeV_2040->RemovePoint(0);
    //WITH material budget error
    TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_2040       = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2040->Get("Pi0SystError");
    //WITHOUT material budget error
    TGraphAsymmErrors* graphPCMPi0InvYieldSysWOMat2760GeV_2040          = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2040->Get("Pi0SystErrorA");
    graphPCMPi0InvYieldSysWOMat2760GeV_2040->RemovePoint(graphPCMPi0InvYieldSysWOMat2760GeV_2040->GetN()-1);
    TGraphAsymmErrors* graphPCMPi0RAAStatPbPb2760GeV_2040   = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2040->Get("Pi0RAA");
    TH1D* histoPCMPi0RAAStatPbPb2760GeV_2040 = (TH1D*)GraphAsymErrorsToHist(graphPCMPi0RAAStatPbPb2760GeV_2040,14,"histoPCMPi0RAAStatPbPb2760GeV_2040");
    TGraphAsymmErrors* graphPCMPi0RAASysPbPb2760GeV_2040    = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2040->Get("Pi0RAASys");

    TDirectory* directoryPCMPi0PbPb2760GeV_2050                     = (TDirectory*)filePCM->Get("Pi0_PbPb_2.76TeV_20-50%");
    TH1D* histoPCMPi0InvYieldPbPb2760GeV_2050                   = (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("CorrectedYieldPi0");
    TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeV_2050  = new TGraphAsymmErrors(histoPCMPi0InvYieldPbPb2760GeV_2050);
    while(graphPCMPi0InvYieldStatPbPb2760GeV_2050->GetY()[0] < 1e-14) graphPCMPi0InvYieldStatPbPb2760GeV_2050->RemovePoint(0);
    TH1D* histoPCMPi0InvYieldPbPb2760GeVYshifted_2050                  = (TH1D*)directoryPCMPi0PbPb2760GeV_2050->Get("CorrectedYieldBinShiftedPi0");
    TGraphAsymmErrors* graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_2050  = new TGraphAsymmErrors(histoPCMPi0InvYieldPbPb2760GeVYshifted_2050);
    while(graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_2050->GetY()[0] < 1e-14) graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_2050->RemovePoint(0);
    //WITH material budget error
    TGraphAsymmErrors* graphPCMPi0InvYieldSysPbPb2760GeV_2050       = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0SystError");
    //WITHOUT material budget error
    TGraphAsymmErrors* graphPCMPi0InvYieldSysWOMat2760GeV_2050          = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0SystErrorA");
    TGraphAsymmErrors* graphPCMPi0RAAStatPbPb2760GeV_2050   = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0RAA");
    TH1D* histoPCMPi0RAAStatPbPb2760GeV_2050 = (TH1D*)GraphAsymErrorsToHist(graphPCMPi0RAAStatPbPb2760GeV_2050,14,"histoPCMPi0RAAStatPbPb2760GeV_2050");
    TGraphAsymmErrors* graphPCMPi0RAASysPbPb2760GeV_2050    = (TGraphAsymmErrors*)directoryPCMPi0PbPb2760GeV_2050->Get("Pi0RAASys");


    //***************************************   Eta PCM   ******************************************//
    TDirectory* directoryPCMEtaPbPb2760GeV_0010                     = (TDirectory*)filePCM->Get("Eta_PbPb_2.76TeV_0-10%");
    TH1D* histoPCMEtaInvYieldPbPb2760GeV_0010                   = (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("CorrectedYieldEta");
    TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeV_0010  = new TGraphAsymmErrors(histoPCMEtaInvYieldPbPb2760GeV_0010);
    graphPCMEtaInvYieldStatPbPb2760GeV_0010->RemovePoint(0);
    graphPCMEtaInvYieldStatPbPb2760GeV_0010->RemovePoint(0);
    TH1D* histoPCMEtaInvYieldPbPb2760GeVYshifted_0010                  = (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("CorrectedYieldBinShiftedEta");
    TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_0010  = new TGraphAsymmErrors(histoPCMEtaInvYieldPbPb2760GeVYshifted_0010);
    graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_0010->RemovePoint(0);
    graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_0010->RemovePoint(0);
    //WITH material budget error
    TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_0010       = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtaSystError");
    //WITHOUT material budget error
    TGraphAsymmErrors* graphPCMEtaInvYieldSysWOMat2760GeV_0010          = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtaSystErrorA");
    TGraphAsymmErrors* graphPCMEtaRAAStatPbPb2760GeV_0010   = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtaRAA");
    TH1D* histoPCMEtaRAAStatPbPb2760GeV_0010 = (TH1D*)GraphAsymErrorsToHist(graphPCMEtaRAAStatPbPb2760GeV_0010,10,"histoPCMEtaRAAStatPbPb2760GeV_0010");
    TGraphAsymmErrors* graphPCMEtaRAASysPbPb2760GeV_0010    = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtaRAASys");
    TH1D* histoPCMEtatoPi0Stat2760GeV_0010                  = (TH1D*)directoryPCMEtaPbPb2760GeV_0010->Get("EtatoPi0Ratio");
    TGraphAsymmErrors* graphPCMEtatoPi0Stat2760GeV_0010 = new TGraphAsymmErrors(histoPCMEtatoPi0Stat2760GeV_0010);
    graphPCMEtatoPi0Stat2760GeV_0010->RemovePoint(0);
    graphPCMEtatoPi0Stat2760GeV_0010->RemovePoint(0);
    TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_0010  = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_0010->Get("EtatoPi0RatioSys");

    TDirectoryFile* directoryPCMEtaPbPb2760GeV_2040                 = (TDirectoryFile*)filePCM->Get("Eta_PbPb_2.76TeV_20-40%");
    TH1D* histoPCMEtaInvYieldPbPb2760GeV_2040                   = (TH1D*)directoryPCMEtaPbPb2760GeV_2040->Get("CorrectedYieldEta");
    TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeV_2040  = new TGraphAsymmErrors(histoPCMEtaInvYieldPbPb2760GeV_2040);
    while(graphPCMEtaInvYieldStatPbPb2760GeV_2040->GetY()[0] < 1e-14) graphPCMEtaInvYieldStatPbPb2760GeV_2040->RemovePoint(0);
    //WITH material budget error
    TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_2040       = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2040->Get("EtaSystError");

    TDirectory* directoryPCMEtaPbPb2760GeV_2050                     = (TDirectory*)filePCM->Get("Eta_PbPb_2.76TeV_20-50%");
    TH1D* histoPCMEtaInvYieldPbPb2760GeV_2050                   = (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("CorrectedYieldEta");
    TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeV_2050  = new TGraphAsymmErrors(histoPCMEtaInvYieldPbPb2760GeV_2050);
    graphPCMEtaInvYieldStatPbPb2760GeV_2050->RemovePoint(0);
    graphPCMEtaInvYieldStatPbPb2760GeV_2050->RemovePoint(0);
    TH1D* histoPCMEtaInvYieldPbPb2760GeVYshifted_2050                  = (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("CorrectedYieldBinShiftedEta");
    TGraphAsymmErrors* graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_2050  = new TGraphAsymmErrors(histoPCMEtaInvYieldPbPb2760GeVYshifted_2050);
    graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_2050->RemovePoint(0);
    graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_2050->RemovePoint(0);
    //WITH material budget error
    TGraphAsymmErrors* graphPCMEtaInvYieldSysPbPb2760GeV_2050       = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtaSystError");
    //WITHOUT material budget error
    TGraphAsymmErrors* graphPCMEtaInvYieldSysWOMat2760GeV_2050  = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtaSystErrorA");
    TGraphAsymmErrors* graphPCMEtaRAAStatPbPb2760GeV_2050   = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtaRAA");
    TH1D* histoPCMEtaRAAStatPbPb2760GeV_2050 = (TH1D*)GraphAsymErrorsToHist(graphPCMEtaRAAStatPbPb2760GeV_2050,10,"histoPCMEtaRAAStatPbPb2760GeV_2050");
    TGraphAsymmErrors* graphPCMEtaRAASysPbPb2760GeV_2050    = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtaRAASys");
    TH1D* histoPCMEtatoPi0Stat2760GeV_2050                  = (TH1D*)directoryPCMEtaPbPb2760GeV_2050->Get("EtatoPi0Ratio");
    TGraphAsymmErrors* graphPCMEtatoPi0Stat2760GeV_2050 = new TGraphAsymmErrors(histoPCMEtatoPi0Stat2760GeV_2050);
    graphPCMEtatoPi0Stat2760GeV_2050->RemovePoint(0);
    graphPCMEtatoPi0Stat2760GeV_2050->RemovePoint(0);
    TGraphAsymmErrors* graphPCMEtatoPi0Sys2760GeV_2050  = (TGraphAsymmErrors*)directoryPCMEtaPbPb2760GeV_2050->Get("EtatoPi0RatioSys");

    //***************************************   Pi0 EMCal   ******************************************//
    TFile* fileEMCal                                = new TFile(fileNameEMCal.Data());
    TDirectory* directoryEMCalPi0PbPb2760GeV                = (TDirectory*)fileEMCal->Get("Pi02.76TeV_PbPb");
    TH1D* histoEMCalPi0InvYieldStatPbPb2760GeV_0010         = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("LuciaBinningInvYieldPbPbStatErrPi_0010");
    TGraphAsymmErrors* graphEMCalPi0InvYieldStatPbPb2760GeV_0010        = new TGraphAsymmErrors(histoEMCalPi0InvYieldStatPbPb2760GeV_0010);
    graphEMCalPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(graphEMCalPi0InvYieldStatPbPb2760GeV_0010->GetN()-1);
    graphEMCalPi0InvYieldStatPbPb2760GeV_0010->RemovePoint(graphEMCalPi0InvYieldStatPbPb2760GeV_0010->GetN()-1);
    TH1D* histoEMCalPi0InvYieldSysPbPb2760GeV_0010  = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("LuciaBinningInvYieldPbPbSysErrPi_0010");
    TGraphAsymmErrors* graphEMCalPi0InvYieldSysPbPb2760GeV_0010     = new TGraphAsymmErrors(histoEMCalPi0InvYieldSysPbPb2760GeV_0010);
    graphEMCalPi0InvYieldSysPbPb2760GeV_0010->RemovePoint(graphEMCalPi0InvYieldSysPbPb2760GeV_0010->GetN()-1);
    graphEMCalPi0InvYieldSysPbPb2760GeV_0010->RemovePoint(graphEMCalPi0InvYieldSysPbPb2760GeV_0010->GetN()-1);
    TH1D* histoEMCalPi0InvYieldSysPbPb2760GeVforRAA_0010  = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("EMCalSysPion010forRAA");
    TGraphAsymmErrors* graphEMCalPi0InvYieldSysPbPb2760GeVforRAA_0010     = new TGraphAsymmErrors(histoEMCalPi0InvYieldSysPbPb2760GeVforRAA_0010);

    TH1D* histoEMCalPi0InvYieldStatPbPb2760GeV_2050         = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("LuciaBinningInvYieldPbPbStatErrPi_2050");
    TGraphAsymmErrors* graphEMCalPi0InvYieldStatPbPb2760GeV_2050        = new TGraphAsymmErrors(histoEMCalPi0InvYieldStatPbPb2760GeV_2050);
    graphEMCalPi0InvYieldStatPbPb2760GeV_2050->RemovePoint(graphEMCalPi0InvYieldStatPbPb2760GeV_2050->GetN()-1);
    graphEMCalPi0InvYieldStatPbPb2760GeV_2050->RemovePoint(graphEMCalPi0InvYieldStatPbPb2760GeV_2050->GetN()-1);
    TH1D* histoEMCalPi0InvYieldSysPbPb2760GeV_2050  = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("LuciaBinningInvYieldPbPbSysErrPi_2050");
    TGraphAsymmErrors* graphEMCalPi0InvYieldSysPbPb2760GeV_2050     = new TGraphAsymmErrors(histoEMCalPi0InvYieldSysPbPb2760GeV_2050);
    graphEMCalPi0InvYieldSysPbPb2760GeV_2050->RemovePoint(graphEMCalPi0InvYieldSysPbPb2760GeV_2050->GetN()-1);
    graphEMCalPi0InvYieldSysPbPb2760GeV_2050->RemovePoint(graphEMCalPi0InvYieldSysPbPb2760GeV_2050->GetN()-1);
    TH1D* histoEMCalPi0InvYieldSysPbPb2760GeVforRAA_2050  = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("EMCalSysPion2050forRAA");
    TGraphAsymmErrors* graphEMCalPi0InvYieldSysPbPb2760GeVforRAA_2050     = new TGraphAsymmErrors(histoEMCalPi0InvYieldSysPbPb2760GeVforRAA_2050);

    //***************************************   Eta EMCal   ******************************************//
    TDirectory* directoryEMCalEtaPbPb2760GeV                = (TDirectory*)fileEMCal->Get("Eta2.76TeV_PbPb");
    TH1D* histoEMCalEtaInvYieldStatPbPb2760GeV_0010         = (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbStatErrEta_0010");
    TGraphAsymmErrors* graphEMCalEtaInvYieldStatPbPb2760GeV_0010        = new TGraphAsymmErrors(histoEMCalEtaInvYieldStatPbPb2760GeV_0010);
    graphEMCalEtaInvYieldStatPbPb2760GeV_0010->RemovePoint(graphEMCalEtaInvYieldStatPbPb2760GeV_0010->GetN()-1);
    graphEMCalEtaInvYieldStatPbPb2760GeV_0010->RemovePoint(graphEMCalEtaInvYieldStatPbPb2760GeV_0010->GetN()-1);
    TH1D*   histoEMCalEtaInvYieldSysPbPb2760GeV_0010    = (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbSysErrEta_0010");
    TGraphAsymmErrors* graphEMCalEtaInvYieldSysPbPb2760GeV_0010     = new TGraphAsymmErrors(histoEMCalEtaInvYieldSysPbPb2760GeV_0010);
    graphEMCalEtaInvYieldSysPbPb2760GeV_0010->RemovePoint(graphEMCalEtaInvYieldSysPbPb2760GeV_0010->GetN()-1);
    graphEMCalEtaInvYieldSysPbPb2760GeV_0010->RemovePoint(graphEMCalEtaInvYieldSysPbPb2760GeV_0010->GetN()-1);
    TH1D* histoEMCalEtaInvYieldSysPbPb2760GeVforRAA_0010  = (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("EMCalSysEta010forRAA");
    TGraphAsymmErrors* graphEMCalEtaInvYieldSysPbPb2760GeVforRAA_0010     = new TGraphAsymmErrors(histoEMCalEtaInvYieldSysPbPb2760GeVforRAA_0010);
    TH1D*   histoEMCalEtatoPi0StatPbPb2760GeV_0010  = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("StatErrEtatoPi0Ratio_0010");
    TGraphAsymmErrors* graphEMCalEtatoPi0Stat2760GeV_0010   = new TGraphAsymmErrors(histoEMCalEtatoPi0StatPbPb2760GeV_0010);
    graphEMCalEtatoPi0Stat2760GeV_0010->RemovePoint(graphEMCalEtatoPi0Stat2760GeV_0010->GetN()-1);
    graphEMCalEtatoPi0Stat2760GeV_0010->RemovePoint(graphEMCalEtatoPi0Stat2760GeV_0010->GetN()-1);
    TH1D*   histoEMCalEtatoPi0SysPbPb2760GeV_0010   = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("SysErrEtatoPi0Ratio_0010");
    TGraphAsymmErrors* graphEMCalEtatoPi0Sys2760GeV_0010    = new TGraphAsymmErrors(histoEMCalEtatoPi0SysPbPb2760GeV_0010);
    graphEMCalEtatoPi0Sys2760GeV_0010->RemovePoint(graphEMCalEtatoPi0Sys2760GeV_0010->GetN()-1);
    graphEMCalEtatoPi0Sys2760GeV_0010->RemovePoint(graphEMCalEtatoPi0Sys2760GeV_0010->GetN()-1);

    TH1D* histoEMCalEtaInvYieldStatPbPb2760GeV_2050         = (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbStatErrEta_2050");
    TGraphAsymmErrors* graphEMCalEtaInvYieldStatPbPb2760GeV_2050        = new TGraphAsymmErrors(histoEMCalEtaInvYieldStatPbPb2760GeV_2050);
    graphEMCalEtaInvYieldStatPbPb2760GeV_2050->RemovePoint(graphEMCalEtaInvYieldStatPbPb2760GeV_2050->GetN()-1);
    graphEMCalEtaInvYieldStatPbPb2760GeV_2050->RemovePoint(graphEMCalEtaInvYieldStatPbPb2760GeV_2050->GetN()-1);
    TH1D*   histoEMCalEtaInvYieldSysPbPb2760GeV_2050    = (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("InvYieldPbPbSysErrEta_2050");
    TGraphAsymmErrors* graphEMCalEtaInvYieldSysPbPb2760GeV_2050     = new TGraphAsymmErrors(histoEMCalEtaInvYieldSysPbPb2760GeV_2050);
    graphEMCalEtaInvYieldSysPbPb2760GeV_2050->RemovePoint(graphEMCalEtaInvYieldSysPbPb2760GeV_2050->GetN()-1);
    graphEMCalEtaInvYieldSysPbPb2760GeV_2050->RemovePoint(graphEMCalEtaInvYieldSysPbPb2760GeV_2050->GetN()-1);
    TH1D* histoEMCalEtaInvYieldSysPbPb2760GeVforRAA_2050  = (TH1D*)directoryEMCalEtaPbPb2760GeV->Get("EMCalSysEta2050forRAA");
    TGraphAsymmErrors* graphEMCalEtaInvYieldSysPbPb2760GeVforRAA_2050     = new TGraphAsymmErrors(histoEMCalEtaInvYieldSysPbPb2760GeVforRAA_2050);
    TH1D*   histoEMCalEtatoPi0StatPbPb2760GeV_2050  = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("StatErrEtatoPi0Ratio_2050");
    TGraphAsymmErrors* graphEMCalEtatoPi0Stat2760GeV_2050   = new TGraphAsymmErrors(histoEMCalEtatoPi0StatPbPb2760GeV_2050);
    graphEMCalEtatoPi0Stat2760GeV_2050->RemovePoint(graphEMCalEtatoPi0Stat2760GeV_2050->GetN()-1);
    graphEMCalEtatoPi0Stat2760GeV_2050->RemovePoint(graphEMCalEtatoPi0Stat2760GeV_2050->GetN()-1);
    TH1D*   histoEMCalEtatoPi0SysPbPb2760GeV_2050   = (TH1D*)directoryEMCalPi0PbPb2760GeV->Get("SysErrEtatoPi0Ratio_2050");
    TGraphAsymmErrors* graphEMCalEtatoPi0Sys2760GeV_2050    = new TGraphAsymmErrors(histoEMCalEtatoPi0SysPbPb2760GeV_2050);
    graphEMCalEtatoPi0Sys2760GeV_2050->RemovePoint(graphEMCalEtatoPi0Sys2760GeV_2050->GetN()-1);
    graphEMCalEtatoPi0Sys2760GeV_2050->RemovePoint(graphEMCalEtatoPi0Sys2760GeV_2050->GetN()-1);

    TDirectory* directoryEMCalRAAPbPb2760GeV                = (TDirectory*)fileEMCal->Get("RAA2.76TeV_PbPb");
    TH1D*   histoEMCalPi0RAAStatPbPb2760GeV_0010    = (TH1D*)directoryEMCalRAAPbPb2760GeV->Get("RAAStatErrPion_0010");
    TGraphAsymmErrors* graphEMCalPi0RAAStatPbPb2760GeV_0010     = new TGraphAsymmErrors(histoEMCalPi0RAAStatPbPb2760GeV_0010);
    TH1D*   histoEMCalPi0RAASysPbPb2760GeV_0010    = (TH1D*)directoryEMCalRAAPbPb2760GeV->Get("RAASystErrPion_0010");
    TGraphAsymmErrors* graphEMCalPi0RAASysPbPb2760GeV_0010  = new TGraphAsymmErrors(histoEMCalPi0RAASysPbPb2760GeV_0010);
    TH1D*   histoEMCalPi0RAAStatPbPb2760GeV_2050    = (TH1D*)directoryEMCalRAAPbPb2760GeV->Get("RAAStatErrPion_2050");
    TGraphAsymmErrors* graphEMCalPi0RAAStatPbPb2760GeV_2050     = new TGraphAsymmErrors(histoEMCalPi0RAAStatPbPb2760GeV_2050);
    TH1D*   histoEMCalPi0RAASysPbPb2760GeV_2050    = (TH1D*)directoryEMCalRAAPbPb2760GeV->Get("RAASystErrPion_2050");
    TGraphAsymmErrors* graphEMCalPi0RAASysPbPb2760GeV_2050  = new TGraphAsymmErrors(histoEMCalPi0RAASysPbPb2760GeV_2050);

    TH1D*   histoEMCalEtaRAAStatPbPb2760GeV_0010    = (TH1D*)directoryEMCalRAAPbPb2760GeV->Get("RAAStatErrEta_0010");
    TGraphAsymmErrors* graphEMCalEtaRAAStatPbPb2760GeV_0010     = new TGraphAsymmErrors(histoEMCalEtaRAAStatPbPb2760GeV_0010);
    TH1D*   histoEMCalEtaRAASysPbPb2760GeV_0010    = (TH1D*)directoryEMCalRAAPbPb2760GeV->Get("RAASystErrEta_0010");
    TGraphAsymmErrors* graphEMCalEtaRAASysPbPb2760GeV_0010  = new TGraphAsymmErrors(histoEMCalEtaRAASysPbPb2760GeV_0010);
    TH1D*   histoEMCalEtaRAAStatPbPb2760GeV_2050    = (TH1D*)directoryEMCalRAAPbPb2760GeV->Get("RAAStatErrEta_2050");
    TGraphAsymmErrors* graphEMCalEtaRAAStatPbPb2760GeV_2050     = new TGraphAsymmErrors(histoEMCalEtaRAAStatPbPb2760GeV_2050);
    TH1D*   histoEMCalEtaRAASysPbPb2760GeV_2050    = (TH1D*)directoryEMCalRAAPbPb2760GeV->Get("RAASystErrEta_2050");
    TGraphAsymmErrors* graphEMCalEtaRAASysPbPb2760GeV_2050  = new TGraphAsymmErrors(histoEMCalEtaRAASysPbPb2760GeV_2050);



  TFile *fInputFilePbPb2760GeV = new TFile(Form("%s/InputALICEResultsPbPb2760GeV_%s.root", outputDir.Data(),dateForOutput.Data()), "RECREATE");
      fInputFilePbPb2760GeV->mkdir("NeutralMesons_PbPb_2.76TeV");
      fInputFilePbPb2760GeV->mkdir("NeutralMesons_PP_2.76TeV");
      fInputFilePbPb2760GeV->mkdir("NeutralMesons_PP_7TeV");
      fInputFilePbPb2760GeV->mkdir("ChargedParticles_PbPb_2.76TeV");

      TDirectoryFile* directoryNeutralMesonPbPb = (TDirectoryFile*)fInputFilePbPb2760GeV->Get("NeutralMesons_PbPb_2.76TeV");
      fInputFilePbPb2760GeV->cd("NeutralMesons_PbPb_2.76TeV");

        graphPCMPublishedInvYieldStatPbPb2760GeV_0010->Write("graphInvYieldPi0PCMPubPbPb2760GeVStatErr_0010");
        graphPCMPublishedInvYieldSysPbPb2760GeV_0010->Write("graphInvYieldPi0PCMPubPbPb2760GeVSysErr_0010");
        graphPCMPublishedInvYieldStatPbPb2760GeV_2040->Write("graphInvYieldPi0PCMPubPbPb2760GeVStatErr_2040");
        graphPCMPublishedInvYieldSysPbPb2760GeV_2040->Write("graphInvYieldPi0PCMPubPbPb2760GeVSysErr_2040");

        histoPi0PHOSPbPb0010->Write("histoInvYieldPi0PHOSPbPb2760GeVStatErr_0010");
        graphPHOSPi0InvYieldStatPbPb2760GeV_0010->Write("graphInvYieldPi0PHOSPbPb2760GeVStatErr_0010");
        graphPHOSPi0InvYieldSysPbPb2760GeV_0010->Write("graphInvYieldPi0PHOSPbPb2760GeVSysErr_0010");
        graphSysErrRAAYieldPi0PHOSPbPb0010->Write("graphInvYieldPi0PHOSPbPb2760GeVSysErr_forRAA_0010");

        histoPCMPi0InvYieldPbPb2760GeV_0010->Write("histoInvYieldPi0PCMPbPb2760GeVStatErr_0010");
        graphPCMPi0InvYieldStatPbPb2760GeV_0010->Write("graphInvYieldPi0PCMPbPb2760GeVStatErr_0010");
        graphPCMPi0InvYieldSysPbPb2760GeV_0010->Write("graphInvYieldPi0PCMPbPb2760GeVSysErr_0010");
        graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_0010->Write("graphInvYieldPi0PCMPbPb2760GeVStatErr_yShifted_0010");
        graphPCMPi0InvYieldSysWOMat2760GeV_0010->Write("graphInvYieldPi0PCMPbPb2760GeVSysErr_forRAA_0010");
        histoPCMPi0RAAStatPbPb2760GeV_0010->Write("histoRAAPi0PCMPbPb2760GeVStatErr_0010");
        graphPCMPi0RAAStatPbPb2760GeV_0010->Write("graphRAAPi0PCMPbPb2760GeVStatErr_0010");
        graphPCMPi0RAASysPbPb2760GeV_0010->Write("graphRAAPi0PCMPbPb2760GeVSysErr_0010");

        histoPCMPi0InvYieldPbPb2760GeV_2040->Write("histoInvYieldPi0PCMPbPb2760GeVStatErr_2040");
        graphPCMPi0InvYieldStatPbPb2760GeV_2040->Write("graphInvYieldPi0PCMPbPb2760GeVStatErr_2040");
        graphPCMPi0InvYieldSysPbPb2760GeV_2040->Write("graphInvYieldPi0PCMPbPb2760GeVSysErr_2040");

        histoPCMPi0InvYieldPbPb2760GeV_2050->Write("histoInvYieldPi0PCMPbPb2760GeVStatErr_2050");
        graphPCMPi0InvYieldStatPbPb2760GeV_2050->Write("graphInvYieldPi0PCMPbPb2760GeVStatErr_2050");
        graphPCMPi0InvYieldSysPbPb2760GeV_2050->Write("graphInvYieldPi0PCMPbPb2760GeVSysErr_2050");
        graphPCMPi0InvYieldStatPbPb2760GeVwithYshift_2050->Write("graphInvYieldPi0PCMPbPb2760GeVStatErr_yShifted_2050");
        graphPCMPi0InvYieldSysWOMat2760GeV_2050->Write("graphInvYieldPi0PCMPbPb2760GeVSysErr_forRAA_2050");
        histoPCMPi0RAAStatPbPb2760GeV_2050->Write("histoRAAPi0PCMPbPb2760GeVStatErr_2050");
        graphPCMPi0RAAStatPbPb2760GeV_2050->Write("graphRAAPi0PCMPbPb2760GeVStatErr_2050");
        graphPCMPi0RAASysPbPb2760GeV_2050->Write("graphRAAPi0PCMPbPb2760GeVSysErr_2050");

        histoPCMEtaInvYieldPbPb2760GeV_0010->Write("histoInvYieldEtaPCMPbPb2760GeVStatErr_0010");
        graphPCMEtaInvYieldStatPbPb2760GeV_0010->Write("graphInvYieldEtaPCMPbPb2760GeVStatErr_0010");
        graphPCMEtaInvYieldSysPbPb2760GeV_0010->Write("graphInvYieldEtaPCMPbPb2760GeVSysErr_0010");
        graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_0010->Write("graphInvYieldEtaPCMPbPb2760GeVStatErr_yShifted_0010");
        graphPCMEtaInvYieldSysWOMat2760GeV_0010->Write("graphInvYieldEtaPCMPbPb2760GeVSysErr_forRAA_0010");
        histoPCMEtaRAAStatPbPb2760GeV_0010->Write("histoRAAEtaPCMPbPb2760GeVStatErr_0010");
        graphPCMEtaRAAStatPbPb2760GeV_0010->Write("graphRAAEtaPCMPbPb2760GeVStatErr_0010");
        graphPCMEtaRAASysPbPb2760GeV_0010->Write("graphRAAEtaPCMPbPb2760GeVSysErr_0010");
        histoPCMEtatoPi0Stat2760GeV_0010->Write("histoEtaToPi0RatioPCMPbPb2760GeVStatErr_0010");
        graphPCMEtatoPi0Stat2760GeV_0010->Write("graphEtaToPi0RatioPCMPbPb2760GeVStatErr_0010");
        graphPCMEtatoPi0Sys2760GeV_0010->Write("graphEtaToPi0RatioPCMPbPb2760GeVSysErr_0010");

        histoPCMEtaInvYieldPbPb2760GeV_2040->Write("histoInvYieldEtaPCMPbPb2760GeVStatErr_2040");
        graphPCMEtaInvYieldStatPbPb2760GeV_2040->Write("graphInvYieldEtaPCMPbPb2760GeVStatErr_2040");
        graphPCMEtaInvYieldSysPbPb2760GeV_2040->Write("graphInvYieldEtaPCMPbPb2760GeVSysErr_2040");

        histoPCMEtaInvYieldPbPb2760GeV_2050->Write("histoInvYieldEtaPCMPbPb2760GeVStatErr_2050");
        graphPCMEtaInvYieldStatPbPb2760GeV_2050->Write("graphInvYieldEtaPCMPbPb2760GeVStatErr_2050");
        graphPCMEtaInvYieldSysPbPb2760GeV_2050->Write("graphInvYieldEtaPCMPbPb2760GeVSysErr_2050");
        graphPCMEtaInvYieldStatPbPb2760GeVwithYshift_2050->Write("graphInvYieldEtaPCMPbPb2760GeVStatErr_yShifted_2050");
        graphPCMEtaInvYieldSysWOMat2760GeV_2050->Write("graphInvYieldEtaPCMPbPb2760GeVSysErr_forRAA_2050");
        histoPCMEtaRAAStatPbPb2760GeV_2050->Write("histoRAAEtaPCMPbPb2760GeVStatErr_2050");
        graphPCMEtaRAAStatPbPb2760GeV_2050->Write("graphRAAEtaPCMPbPb2760GeVStatErr_2050");
        graphPCMEtaRAASysPbPb2760GeV_2050->Write("graphRAAEtaPCMPbPb2760GeVSysErr_2050");
        histoPCMEtatoPi0Stat2760GeV_2050->Write("histoEtaToPi0RatioPCMPbPb2760GeVStatErr_2050");
        graphPCMEtatoPi0Stat2760GeV_2050->Write("graphEtaToPi0RatioPCMPbPb2760GeVStatErr_2050");
        graphPCMEtatoPi0Sys2760GeV_2050->Write("graphEtaToPi0RatioPCMPbPb2760GeVSysErr_2050");

        histoEMCalPi0InvYieldStatPbPb2760GeV_0010->Write("histoInvYieldPi0EMCalPbPb2760GeVStatErr_0010");
        graphEMCalPi0InvYieldStatPbPb2760GeV_0010->Write("graphInvYieldPi0EMCalPbPb2760GeVStatErr_0010");
        graphEMCalPi0InvYieldSysPbPb2760GeV_0010->Write("graphInvYieldPi0EMCalPbPb2760GeVSysErr_0010");
        graphEMCalPi0InvYieldSysPbPb2760GeVforRAA_0010->Write("graphInvYieldPi0EMCalPbPb2760GeVSysErr_forRAA_0010");
        histoEMCalPi0RAAStatPbPb2760GeV_0010->Write("histoRAAPi0EMCalPbPb2760GeVStatErr_0010");
        graphEMCalPi0RAAStatPbPb2760GeV_0010->Write("graphRAAPi0EMCalPbPb2760GeVStatErr_0010");
        graphEMCalPi0RAASysPbPb2760GeV_0010->Write("graphRAAPi0EMCalPbPb2760GeVSysErr_0010");

        histoEMCalPi0InvYieldStatPbPb2760GeV_2050->Write("histoInvYieldPi0EMCalPbPb2760GeVStatErr_2050");
        graphEMCalPi0InvYieldStatPbPb2760GeV_2050->Write("graphInvYieldPi0EMCalPbPb2760GeVStatErr_2050");
        graphEMCalPi0InvYieldSysPbPb2760GeV_2050->Write("graphInvYieldPi0EMCalPbPb2760GeVSysErr_2050");
        graphEMCalPi0InvYieldSysPbPb2760GeVforRAA_2050->Write("graphInvYieldPi0EMCalPbPb2760GeVSysErr_forRAA_2050");
        histoEMCalPi0RAAStatPbPb2760GeV_2050->Write("histoRAAPi0EMCalPbPb2760GeVStatErr_2050");
        graphEMCalPi0RAAStatPbPb2760GeV_2050->Write("graphRAAPi0EMCalPbPb2760GeVStatErr_2050");
        graphEMCalPi0RAASysPbPb2760GeV_2050->Write("graphRAAPi0EMCalPbPb2760GeVSysErr_2050");

        histoEMCalEtaInvYieldStatPbPb2760GeV_0010->Write("histoInvYieldEtaEMCalPbPb2760GeVStatErr_0010");
        graphEMCalEtaInvYieldStatPbPb2760GeV_0010->Write("graphInvYieldEtaEMCalPbPb2760GeVStatErr_0010");
        graphEMCalEtaInvYieldSysPbPb2760GeV_0010->Write("graphInvYieldEtaEMCalPbPb2760GeVSysErr_0010");
        graphEMCalEtaInvYieldSysPbPb2760GeVforRAA_0010->Write("graphInvYieldEtaEMCalPbPb2760GeVSysErr_forRAA_0010");
        histoEMCalEtaRAAStatPbPb2760GeV_0010->Write("histoRAAEtaEMCalPbPb2760GeVStatErr_0010");
        graphEMCalEtaRAAStatPbPb2760GeV_0010->Write("graphRAAEtaEMCalPbPb2760GeVStatErr_0010");
        graphEMCalEtaRAASysPbPb2760GeV_0010->Write("graphRAAEtaEMCalPbPb2760GeVSysErr_0010");
        histoEMCalEtatoPi0StatPbPb2760GeV_0010->Write("histoEtaToPi0RatioEMCalPbPb2760GeVStatErr_0010");
        graphEMCalEtatoPi0Stat2760GeV_0010->Write("graphEtaToPi0RatioEMCalPbPb2760GeVStatErr_0010");
        graphEMCalEtatoPi0Sys2760GeV_0010->Write("graphEtaToPi0RatioEMCalPbPb2760GeVSysErr_0010");

        histoEMCalEtaInvYieldStatPbPb2760GeV_2050->Write("histoInvYieldEtaEMCalPbPb2760GeVStatErr_2050");
        graphEMCalEtaInvYieldStatPbPb2760GeV_2050->Write("graphInvYieldEtaEMCalPbPb2760GeVStatErr_2050");
        graphEMCalEtaInvYieldSysPbPb2760GeV_2050->Write("graphInvYieldEtaEMCalPbPb2760GeVSysErr_2050");
        graphEMCalEtaInvYieldSysPbPb2760GeVforRAA_2050->Write("graphInvYieldEtaEMCalPbPb2760GeVSysErr_forRAA_2050");
        histoEMCalEtaRAAStatPbPb2760GeV_2050->Write("histoRAAEtaEMCalPbPb2760GeVStatErr_2050");
        graphEMCalEtaRAAStatPbPb2760GeV_2050->Write("graphRAAEtaEMCalPbPb2760GeVStatErr_2050");
        graphEMCalEtaRAASysPbPb2760GeV_2050->Write("graphRAAEtaEMCalPbPb2760GeVSysErr_2050");
        histoEMCalEtatoPi0StatPbPb2760GeV_2050->Write("histoEtaToPi0RatioEMCalPbPb2760GeVStatErr_2050");
        graphEMCalEtatoPi0Stat2760GeV_2050->Write("graphEtaToPi0RatioEMCalPbPb2760GeVStatErr_2050");
        graphEMCalEtatoPi0Sys2760GeV_2050->Write("graphEtaToPi0RatioEMCalPbPb2760GeVSysErr_2050");

      TDirectoryFile* directoryNeutralMesonPP7TeV = (TDirectoryFile*)fInputFilePbPb2760GeV->Get("NeutralMesons_PP_7TeV");
      fInputFilePbPb2760GeV->cd("NeutralMesons_PP_7TeV");
      
        graphCombEtaToPi0Ratiopp7TeVNoXErrors->Write("graphCombEtaToPi0Ratiopp7TeVNoXErrors");
        graphCombEtaToPi0RatioSysErrpp7TeV->Write("graphCombEtaToPi0RatioSysErrpp7TeV");
        
      TDirectoryFile* directoryNeutralMesonPP = (TDirectoryFile*)fInputFilePbPb2760GeV->Get("NeutralMesons_PP_2.76TeV");
      fInputFilePbPb2760GeV->cd("NeutralMesons_PP_2.76TeV");

        graphInvSectionCombStatPi02760GeV->Write("graphInvCrossSectionPi0Comb2760GeVAStatErr");
        graphInvSectionCombSysPi02760GeV->Write("graphInvCrossSectionPi0Comb2760GeVASysErr");
        fitInvCrossSectionTsallisPi0Comb2760GeV->Write("TsallisFitPi0");
        fitInvCrossSectionTCMPi0Comb2760GeV->Write("TwoComponentModelFitPi0");
        graphInvSectionCombStatPi02760GeVforRAA->Write("graphInvCrossSectionPi0Comb2760GeVAStatErr_yShifted");
        graphInvSectionCombSysPi02760GeVforRAA->Write("graphInvCrossSectionPi0Comb2760GeVASysErr_yShifted");

        graphInvSectionPCMStatPi02760GeV->Write("graphInvCrossSectionPi0PCM2760GeVStatErr");
        graphInvSectionPCMSysPi02760GeV->Write("graphInvCrossSectionPi0PCM2760GeVSysErr");
        graphInvSectionPCMStatPi02760GeVforRAA->Write("graphInvCrossSectionPi0PCM2760GeVStatErr_yShifted");
        graphInvSectionPCMSysPi02760GeV_yShifted->Write("graphInvCrossSectionPi0PCM2760GeVSysErr_yShifted");
        graphInvSectionPCMSysPi02760GeVforRAA->Write("graphInvCrossSectionPi0PCM2760GeVSysErr_forRAA");

        graphInvSectionPHOSStatPi02760GeV->Write("graphInvCrossSectionPi0PHOS2760GeVStatErr");
        graphInvSectionPHOSSysPi02760GeV->Write("graphInvCrossSectionPi0PHOS2760GeVSysErr");
        graphInvSectionPHOSStatPi02760GeVforRAA->Write("graphInvCrossSectionPi0PHOS2760GeVStatErr_yShifted");
        graphInvSectionPHOSSysPi02760GeV_yShifted->Write("graphInvCrossSectionPi0PHOS2760GeVSysErr_yShifted");
        graphInvSectionPHOSSysPi02760GeVforRAA->Write("graphInvCrossSectionPi0PHOS2760GeVSysErr_forRAA");

        graphInvSectionEMCalStatPi02760GeV->Write("graphInvCrossSectionPi0EMCAL2760GeVStatErr");
        graphInvSectionEMCalSysPi02760GeV->Write("graphInvCrossSectionPi0EMCAL2760GeVSysErr");
        graphInvSectionEMCalStatPi02760GeVforRAA->Write("graphInvCrossSectionPi0EMCAL2760GeVStatErr_yShifted");
        graphInvSectionEMCalSysPi02760GeV_yShifted->Write("graphInvCrossSectionPi0EMCAL2760GeVSysErr_yShifted");
        graphInvSectionEMCalSysPi02760GeVforRAA->Write("graphInvCrossSectionPi0EMCAL2760GeVSysErr_forRAA");

        graphInvSectionCombStatEta2760GeV->Write("graphInvCrossSectionEtaComb2760GeVAStatErr");
        graphInvSectionCombSysEta2760GeV->Write("graphInvCrossSectionEtaComb2760GeVASysErr");
        fitInvCrossSectionTsallisEtaComb2760GeV->Write("TsallisFitEta");
        fitInvCrossSectionTCMEtaComb2760GeV->Write("TwoComponentModelFitEta");
        graphInvSectionCombStatEta2760GeV_yShifted->Write("graphInvCrossSectionEtaComb2760GeVAStatErr_yShifted");
        graphInvSectionCombSysEta2760GeV_yShifted->Write("graphInvCrossSectionEtaComb2760GeVASysErr_yShifted");
        graphRatioEtaToPi0Comb2760GeVStatErr->Write("graphRatioEtaToPi0Comb2760GeVStatErr");
        graphRatioEtaToPi0Comb2760GeVSysErr->Write("graphRatioEtaToPi0Comb2760GeVSysErr");

        graphInvSectionPCMStatEta2760GeV->Write("graphInvCrossSectionEtaPCM2760GeVStatErr");
        graphInvSectionPCMSysEta2760GeV->Write("graphInvCrossSectionEtaPCM2760GeVSysErr");
        graphInvSectionPCMStatEta2760GeVforRAA->Write("graphInvCrossSectionEtaPCM2760GeVStatErr_yShifted");
        graphInvSectionPCMSysEta2760GeV_yShifted->Write("graphInvCrossSectionEtaPCM2760GeVSysErr_yShifted");
        graphInvSectionPCMSysEta2760GeVforRAA->Write("graphInvCrossSectionEtaPCM2760GeVSysErr_forRAA");

        graphInvSectionEMCalStatEta2760GeV->Write("graphInvCrossSectionEtaEMCAL2760GeVStatErr");
        graphInvSectionEMCalSysEta2760GeV->Write("graphInvCrossSectionEtaEMCAL2760GeVSysErr");
        graphInvSectionEMCalStatEta2760GeVforRAA->Write("graphInvCrossSectionEtaEMCAL2760GeVStatErr_yShifted");
        graphInvSectionEMCalSysEta2760GeV_yShifted->Write("graphInvCrossSectionEtaEMCAL2760GeVSysErr_yShifted");
        graphInvSectionEMCalSysEta2760GeVforRAA->Write("graphInvCrossSectionEtaEMCAL2760GeVSysErr_forRAA");

      TDirectoryFile* directoryChargedPbPb = (TDirectoryFile*)fInputFilePbPb2760GeV->Get("ChargedParticles_PbPb_2.76TeV");
      fInputFilePbPb2760GeV->cd("ChargedParticles_PbPb_2.76TeV");

        histoChargedPionSpectraStat0005->Write("histoChargedPionSpectraStat0005");
        histoChargedPionSpectraSyst0005->Write("histoChargedPionSpectraSyst0005");
        histoChargedPionSpectraStat0510->Write("histoChargedPionSpectraStat0510");
        histoChargedPionSpectraSyst0510->Write("histoChargedPionSpectraSyst0510");
        histoChargedPionSpectraStat0010->Write("histoChargedPionSpectraStat0010");
        histoChargedPionSpectraSyst0010->Write("histoChargedPionSpectraSyst0010");
        histoChargedPionSpectraStat2040->Write("histoChargedPionSpectraStat2040");
        histoChargedPionSpectraSyst2040->Write("histoChargedPionSpectraSyst2040");

        histoChargedKaonSpectraStat0005->Write("histoChargedKaonSpectraStat0005");
        histoChargedKaonSpectraSyst0005->Write("histoChargedKaonSpectraSyst0005");
        histoChargedKaonSpectraStat0510->Write("histoChargedKaonSpectraStat0510");
        histoChargedKaonSpectraSyst0510->Write("histoChargedKaonSpectraSyst0510");
        histoChargedKaonSpectraStat0010->Write("histoChargedKaonSpectraStat0010");
        histoChargedKaonSpectraSyst0010->Write("histoChargedKaonSpectraSyst0010");
        histoChargedKaonSpectraStat2040->Write("histoChargedKaonSpectraStat2040");
        histoChargedKaonSpectraSyst2040->Write("histoChargedKaonSpectraSyst2040");

        graphChargedRatioKaonToPion0005->Write("graphChargedRatioKaonToPion0005");
        graphChargedRatioKaonToPionSys0005->Write("graphChargedRatioKaonToPionSys0005");
        graphChargedRatioKaonToPion0510->Write("graphChargedRatioKaonToPion0510");
        graphChargedRatioKaonToPionSys0510->Write("graphChargedRatioKaonToPionSys0510");
        graphChargedRatioKaonToPion0010->Write("graphChargedRatioKaonToPion0010");
        graphChargedRatioKaonToPionSys0010->Write("graphChargedRatioKaonToPionSys0010");
        graphChargedRatioKaonToPion2040->Write("graphChargedRatioKaonToPion2040");
        graphChargedRatioKaonToPionSys2040->Write("graphChargedRatioKaonToPionSys2040");

        graphChargedPionRAA0005->Write("graphChargedPionRAA0005");
        graphChargedPionRAASys0005->Write("graphChargedPionRAASys0005");
        graphChargedPionRAA0510->Write("graphChargedPionRAA0510");
        graphChargedPionRAASys0510->Write("graphChargedPionRAASys0510");
        graphChargedPionRAA0010->Write("graphChargedPionRAA0010");
        graphChargedPionRAASys0010->Write("graphChargedPionRAASys0010");
        graphChargedPionRAA2040->Write("graphChargedPionRAA2040");
        graphChargedPionRAASys2040->Write("graphChargedPionRAASys2040");

        graphChargedKaonRAA0005->Write("graphChargedKaonRAA0005");
        graphChargedKaonRAASys0005->Write("graphChargedKaonRAASys0005");
        graphChargedKaonRAA0510->Write("graphChargedKaonRAA0510");
        graphChargedKaonRAASys0510->Write("graphChargedKaonRAASys0510");
        graphChargedKaonRAA0010->Write("graphChargedKaonRAA0010");
        graphChargedKaonRAASys0010->Write("graphChargedKaonRAASys0010");
        graphChargedKaonRAA2040->Write("graphChargedKaonRAA2040");
        graphChargedKaonRAASys2040->Write("graphChargedKaonRAASys2040");

  fInputFilePbPb2760GeV->Close();

}
