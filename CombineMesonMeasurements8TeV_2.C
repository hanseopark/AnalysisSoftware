/****************************************************************************************************************************
******         provided by Gamma Conversion Group, PWG4,                                                     *****
******        Daniel Mühlheim, d.muehlheim@cern.ch                                                        *****
******        Friederike Bock, friederike.bock@cern.ch                                                    *****
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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ExtractSignalBinning.h"
// #include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

extern TRandom*    gRandom;
extern TBenchmark* gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*    gMinuit;

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    // TString name;
};

//____________________________________________________________________________________________________________________________________________
void CombineMesonMeasurements8TeV_2(    TString fileNamePCM         = "",
                                        TString fileNamePCMEMCAL    = "",
                                        TString fileNameEMCALLow    = "",
                                        TString fileNameEMCALmerged = "",
                                        TString fileNamePHOS        = "",
                                        TString suffix              = "eps",
                                        TString isMC                = "",
                                        TString thesisPlots         = "",
                                        TString bWCorrection        = "X",
                                        TString fileInputCorrFactors= "",
                                        Bool_t plotInvMassBins = kTRUE
                                    ){

    TString date = ReturnDateString();

    gROOT->Reset();
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();
    SetPlotStyle();

    gStyle->SetEndErrorSize(0);

    TString dateForOutput                       = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;
    //___________________________________ Declaration of files _____________________________________________
    TString collisionSystem8TeV              = "pp, #sqrt{#it{s}} = 8 TeV";

    TString fileNameTheory                      = "ExternalInput/Theory/TheoryCompilationPP.root";
    TString fileNameTheoryPythia8               = "ExternalInput/Theory/Pythia/PYTHIA8_Monash2013Tune_8TeV_Pi0_Eta.root";
    TString fileNameTheoryPythia8_4C            = "ExternalInput/Theory/Pythia/PYTHIA8_4CTune_8TeV_Pi0_Eta.root";
    TString fileNameTheoryPi0DSS14              = "ExternalInput/Theory/pp8TeV_DSS14.root";
    TString fileNameEtaToPi07TeV                = "ExternalInput/EtaPi0Ratio7000GeV.root";
    TString fileNameEtaToPi02760GeV             = "ExternalInput/CombinedResultsPaperPP2760GeV_2016_07_05_FrediV2Clusterizer.root";

    TString fileNamePHOSMB                      = "ExternalInput/PHOS/8TeV/data_PHOS-MBResultsFullCorrection_PP_NoBinShifting.root";
    TString fileNamePHOSPHOS                    = "ExternalInput/PHOS/8TeV/data_PHOS-PHOSResultsFullCorrection_PP_NoBinShifting.root";
    TString fileNamePHOSprelim                  = "ExternalInput/PHOS/8TeV/preliminary/pi0_8TeV_PHOS.root";
//    TString fileNamePCMEta                      = "FinalResults/data_PCMResultsFullCorrection_PP_NoBinShifting_usedTosvnEtapaper.root";
    TString fileNameChargedPionPP               = "ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_20_May_2015.root";
    TString fileNameChargedHadronPP             = "ExternalInput/UnidentifiedCharged/ChargedHadrinSpectraPP_20_May_2015.root";
    TString outputDir                           = Form("%s/%s/CombineMesonMeasurements8TeV_merged%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());

    TString nameFinalResDat                     = Form("%s/CombinedResults%s_FitResults.dat",outputDir.Data(),bWCorrection.Data());
    cout << outputDir.Data() << endl;
    cout << fileNamePCM.Data() << endl;

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
//    gSystem->Exec(Form("cp %s %s/InputPCMEta.root", fileNamePCMEta.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMEMCAL.root", fileNamePCMEMCAL.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCALLow.root", fileNameEMCALLow.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCALmerged.root", fileNameEMCALmerged.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/Theory.root", fileNameTheory.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/ChargedPionsPP.root", fileNameChargedPionPP.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/ChargedHadronsPP.root", fileNameChargedHadronPP.Data(), outputDir.Data()));

    Bool_t thesis                               = kFALSE;
    if(thesisPlots.CompareTo("thesis") == 0){
        thesis                                  = kTRUE;
    }

    TString prefix2                             = "";
    if (isMC.CompareTo("kTRUE")==0){
        prefix2                                 = "MC";
    } else {
        prefix2                                 = "Data";
    }

    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();

    Width_t  widthLinesBoxes                    = 1.4;
    Width_t  widthCommonFit                     = 2;

    // Definition of colors, styles and markers sizes
    Color_t  colorComb                          = kMagenta+2;
    Style_t  markerStyleComb                    = 20;
    Size_t   markerSizeComb                     = 2;

    Color_t  colorCombLowPt                     = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kFALSE);
    Color_t  colorCombHighPt                    = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kTRUE);
    Style_t  markerStyleCombLowPt               = 20;
    Style_t  markerStyleCombHighPt              = 20;
    Size_t   markerSizeComparison               = 2;

    Color_t colorTrigg      [10]                = {kBlack, kGray+1, kRed+2, kBlue+2, kGreen+3, kCyan+2, kViolet, kMagenta+2,  kRed-2, kBlue-2};
    Color_t colorTriggShade [10]                = {kGray+1, kGray, kRed-6, kBlue-6, kGreen-8, kCyan-6, kViolet-8, kMagenta-8,  kRed-8, kBlue-8};
    Marker_t markerTrigg    [10]                = {20, 20, 21, 34, 29, 33, 21, 27, 28, 30 };
    Marker_t markerTriggMC  [10]                = {24, 24, 25, 28, 30, 27, 25, 27, 28, 30 };

    Size_t sizeTrigg        [10]                = {1.5, 1.5, 1.5, 2, 2.2, 2., 1.5, 2., 2.5, 1.5 };


    TString  nameMeasGlobal[11]                 = { "PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal",
                                                    "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged",
                                                    "PCMOtherDataset"};
    TString  nameTrigger[3]                     = {"INT7", "EMC7", "EGA"};

    Color_t  colorDet[11];
    Color_t  colorDetMC[11];
    Style_t  markerStyleDet[11];
    Style_t  markerStyleDetMC[11];
    Size_t   markerSizeDet[11];
    Size_t   markerSizeDetMC[11];

    Color_t  colorNLO                           = kAzure-4;
    Style_t  styleMarkerNLOMuHalf               = 24;
    Style_t  styleMarkerNLOMuOne                = 27;
    Style_t  styleMarkerNLOMuTwo                = 30;
    Style_t  styleLineNLOMuHalf                 = 8;
    Style_t  styleLineNLOMuOne                  = 7;
    Style_t  styleLineNLOMuTwo                  = 4;
    Style_t  styleLineNLOMuTwoBKK               = 3;
    Style_t  styleLineNLOMuTwoDSS               = 6;
    Size_t   sizeMarkerNLO                      = 1;
    Width_t  widthLineNLO                       = 2.;

    Int_t    binExInvMass[11]                   = { 7,  0,  12, 0,  7,
                                                    0,  0,  0,  0,  0,
                                                    0 };
    Double_t rangePtExMin[11]                   = { 1.6,    1.4,    2.6,    0.0,    1.6,
                                                    0.0,    0.0,    0.0,    9.0,    0.0,
                                                    20.0 };

    Double_t rangePtExMax[11]                   = { 1.8,    1.6,    3.0,    0.0,    1.8,
                                                    0.0,    0.0,    0.0,    10.,    0.0,
                                                    22.0 };

    for (Int_t i = 0; i < 11; i++){
        colorDet[i]                             = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kFALSE, kFALSE, kTRUE);
        colorDetMC[i]                           = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kTRUE, kFALSE, kTRUE);
        markerStyleDet[i]                       = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
        markerStyleDetMC[i]                     = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kTRUE);
        markerSizeDet[i]                        = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kFALSE)*2;
        markerSizeDetMC[i]                      = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kTRUE)*2;
    }


    //************************** Read data for PCM **************************************************
    TFile* filePCM                                          = new TFile(fileNamePCM.Data());
    TH1D* histoPCMNumberOfEvents                            = (TH1D*)filePCM->Get("histoNumberOfEvents8TeV");
    TDirectory* directoryPCMPi0                             = (TDirectory*)filePCM->Get("Pi08TeV");
        TH1D* histoPCMPi0Mass                               = (TH1D*)directoryPCMPi0->Get("MassPi0");
        TH1D* histoPCMPi0FWHMMeV                            = (TH1D*)directoryPCMPi0->Get("FWHMPi0MeV");
        TH1D* histoPCMPi0TrueMass                           = (TH1D*)directoryPCMPi0->Get("TrueMassPi0");
        TH1D* histoPCMPi0TrueFWHMMeV                        = (TH1D*)directoryPCMPi0->Get("TrueFWHMPi0MeV");
        TH1D* histoPCMPi0Acc                                = (TH1D*)directoryPCMPi0->Get("AcceptancePi0");
        TH1D* histoPCMPi0TrueEffPt                          = (TH1D*)directoryPCMPi0->Get("EfficiencyPi0");
        TH1D* histoPCMPi0InvXSectionStat                    = (TH1D*)directoryPCMPi0->Get("InvCrossSectionPi0");
        histoPCMPi0InvXSectionStat->SetBinContent(1,0);
        histoPCMPi0InvXSectionStat->SetBinError(1,1);
        TGraphAsymmErrors* graphPCMPi0InvXSectionStat       = new TGraphAsymmErrors(histoPCMPi0InvXSectionStat);
        //graphPCMPi0InvXSectionStat->RemovePoint(graphPCMPi0InvXSectionStat->GetN()-1);
        graphPCMPi0InvXSectionStat->RemovePoint(0);
        cout << "Pi0 stat PCM" << endl;
        graphPCMPi0InvXSectionStat->Print();
        TGraphAsymmErrors* graphPCMPi0InvXSectionSysA       = (TGraphAsymmErrors*)directoryPCMPi0->Get("InvCrossSectionPi0SysA");
        TGraphAsymmErrors* graphPCMPi0InvXSectionSys        = (TGraphAsymmErrors*)directoryPCMPi0->Get("InvCrossSectionPi0Sys");
        //graphPCMPi0InvXSectionSys->RemovePoint(graphPCMPi0InvXSectionSys->GetN()-1);
        graphPCMPi0InvXSectionSys->RemovePoint(graphPCMPi0InvXSectionSys->GetN()-1);
        cout << "Pi0 sys PCM" << endl;
        graphPCMPi0InvXSectionSys->Print();
//      TGraphAsymmErrors* graphPCMPi0CorrYieldSysErr    = (TGraphAsymmErrors*)directoryPCMPi0->Get("Pi0SystError");

        TH1D* histoPCMPi0AccTimesEff                        = (TH1D*)histoPCMPi0TrueEffPt->Clone("histoPCMPi0AccTimesEff");
        histoPCMPi0AccTimesEff->Multiply(histoPCMPi0Acc);
        // normalize to full acceptance (delta y and phi)
        histoPCMPi0AccTimesEff->Scale(2*TMath::Pi()*1.6);

        TH1D* histoPi0InvMassSigPlusBGPCM[3];
        TH1D* histoPi0InvMassSigPCM[3];
        TH1D* histoPi0InvMassSigRemBGSubPCM[3];
        TH1D* histoPi0InvMassBGPCM[3];
        TH1D* histoPi0InvMassRemBGPCM[3];
        TH1D* histoPi0InvMassBGTotPCM[3];
        TF1* fitPi0InvMassSigPCM[3];
        TF1* fitPi0InvMassBGPCM[3];
        Bool_t haveAllPi0InvMassPCM[3]                 = {kFALSE, kFALSE, kFALSE};
        if (plotInvMassBins){
            for (Int_t i = 0; i < 1; i++){
                histoPi0InvMassSigPCM[i]               = (TH1D*)directoryPCMPi0->Get("InvMassSig_PtBin07");
                histoPi0InvMassSigPlusBGPCM[i]         = (TH1D*)directoryPCMPi0->Get("InvMassSigPlusBG_PtBin07");
                histoPi0InvMassBGPCM[i]                = (TH1D*)directoryPCMPi0->Get("InvMassBG_PtBin07");
                fitPi0InvMassSigPCM[i]                 = (TF1*)directoryPCMPi0->Get("FitInvMassSig_PtBin07");
                if (histoPi0InvMassSigPCM[i] && histoPi0InvMassSigPlusBGPCM[i] && histoPi0InvMassBGPCM[i] && fitPi0InvMassSigPCM[i]){
                    haveAllPi0InvMassPCM[i]            = kTRUE;
                }


                if (haveAllPi0InvMassPCM[i]){
                    histoPi0InvMassSigPCM[i]->Fit(fitPi0InvMassSigPCM[i],"QRME0");
                    for (Int_t l=0; l < 6; l++){
                        cout << fitPi0InvMassSigPCM[i]->GetParameter(l) << "\t +- " << fitPi0InvMassSigPCM[i]->GetParError(l) << endl;
                    }
                    fitPi0InvMassBGPCM[i]                                  = new TF1("Linearpp","[0]+[1]*x",0.02,0.25);
                    fitPi0InvMassBGPCM[i]->SetParameter(0, fitPi0InvMassSigPCM[i]->GetParameter(4));
                    fitPi0InvMassBGPCM[i]->SetParameter(1, fitPi0InvMassSigPCM[i]->GetParameter(5));
                    TVirtualFitter * fitterPCM                             = TVirtualFitter::GetFitter();
                    Int_t nFreeParPCM                                      = fitPi0InvMassSigPCM[i]->GetNumberFreeParameters();
                    double * covMatrixPCM                                  = fitterPCM->GetCovarianceMatrix();

                    histoPi0InvMassRemBGPCM[i]                             = (TH1D*)histoPi0InvMassBGPCM[i]->Clone(Form("Pi0_InvMassRemBG_Example_%s",nameTrigger[i].Data()));
                    for (Int_t j = 1; j < histoPi0InvMassRemBGPCM[i]->GetNbinsX()+1; j++){
                        histoPi0InvMassRemBGPCM[i]->SetBinContent(j,0);
                        histoPi0InvMassRemBGPCM[i]->SetBinError(j,0);
                    }
                    for (Int_t j = histoPi0InvMassSigPCM[i]->GetXaxis()->FindBin(0.01); j < histoPi0InvMassSigPCM[i]->GetXaxis()->FindBin(0.30)+1; j++){
                        Double_t startBinEdge                                   = histoPi0InvMassSigPCM[i]->GetXaxis()->GetBinLowEdge(j);
                        Double_t endBinEdge                                     = histoPi0InvMassSigPCM[i]->GetXaxis()->GetBinUpEdge(j);
                        Double_t intLinearBack                                  = fitPi0InvMassBGPCM[i]->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
                        Double_t errorLinearBck                                 = pow(( pow( (endBinEdge-startBinEdge)*fitPi0InvMassSigPCM[i]->GetParError(4),2) +
                                                                                        pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPi0InvMassSigPCM[i]->GetParError(5),2)
                                                                                        +2*covMatrixPCM[nFreeParPCM*nFreeParPCM-2]*(endBinEdge-startBinEdge)*0.5*
                                                                                        (endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
//                         cout << j << "\t" << intLinearBack << "\t" << errorLinearBck << endl;
                        histoPi0InvMassRemBGPCM[i]->SetBinContent(j,intLinearBack);
                        histoPi0InvMassRemBGPCM[i]->SetBinError(j,errorLinearBck);
                    }
                    histoPi0InvMassBGTotPCM[i]         = (TH1D*)histoPi0InvMassBGPCM[i]->Clone(Form("Pi0_InvMassTotBG_Example_%s",nameTrigger[i].Data()));
                    histoPi0InvMassBGTotPCM[i]->Sumw2();
                    histoPi0InvMassBGTotPCM[i]->Add(histoPi0InvMassRemBGPCM[i]);
                    histoPi0InvMassSigRemBGSubPCM[i]   = (TH1D*)histoPi0InvMassSigPCM[i]->Clone(Form("Pi0_InvMassSigRemBGSub_Example_%s",nameTrigger[i].Data()));
                    histoPi0InvMassSigRemBGSubPCM[i]->Sumw2();
                    histoPi0InvMassSigRemBGSubPCM[i]->Add(histoPi0InvMassRemBGPCM[i],-1);
                    fitPi0InvMassSigPCM[i]->SetParameter(4, 0);
                    fitPi0InvMassSigPCM[i]->SetParameter(5, 0);
                }
                cout << nameTrigger[i].Data() << "\t" << histoPi0InvMassSigPCM[i] << "\t" << histoPi0InvMassSigPlusBGPCM[i] << "\t" << histoPi0InvMassBGPCM[i] << "\t" << fitPi0InvMassSigPCM[i]
                      << "\t" << haveAllPi0InvMassPCM[i] << endl;
            }
        }

    cout << "here" << endl;
    //************************** Read data for PCMEMCAL **************************************************
    TFile* filePCMEMCAL                                     = new TFile(fileNamePCMEMCAL.Data());
    TH1D* histoPCMEMCALNumberOfEvents                       = (TH1D*)filePCMEMCAL->Get("histoNumberOfEvents8TeV");
    TDirectory* directoryPCMEMCALPi0                        = (TDirectory*)filePCMEMCAL->Get("Pi08TeV");
        TGraphAsymmErrors* graphPCMEMCALPi0Mass             = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("Pi0_Mass_data");
        graphPCMEMCALPi0Mass                                = ScaleGraph(graphPCMEMCALPi0Mass, 1000.);
        TGraphAsymmErrors* graphPCMEMCALPi0FWHM             = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("Pi0_Width_data");
        graphPCMEMCALPi0FWHM                                = ScaleGraph(graphPCMEMCALPi0FWHM, 1000.);
        TGraphAsymmErrors* graphPCMEMCALPi0MassMC           = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("Pi0_Mass_MC");
        graphPCMEMCALPi0MassMC                              = ScaleGraph(graphPCMEMCALPi0MassMC, 1000.);
        TGraphAsymmErrors* graphPCMEMCALPi0FWHMMC           = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("Pi0_Width_MC");
        graphPCMEMCALPi0FWHMMC                              = ScaleGraph(graphPCMEMCALPi0FWHMMC, 1000.);
        TGraphAsymmErrors* graphPCMEMCALPi0Acc              = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("AcceptancePi0");
        TGraphAsymmErrors* graphPCMEMCALPi0EffPt            = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("EfficiencyPi0");
        TH1D* histoPCMEMCALPi0TriggerEff[4];
        for (Int_t i = 0; i < 3; i++){
            histoPCMEMCALPi0TriggerEff[i]                 = (TH1D*)directoryPCMEMCALPi0->Get(Form("TriggerEfficiencyPi0_%s",nameTrigger[i].Data()));
        }
        TGraphAsymmErrors* graphPCMEMCALPi0AccTimesEff      = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("EffTimesAccPi0");
        TGraphAsymmErrors* graphPCMEMCALPi0InvXSectionStat  = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("graphInvCrossSectionPi0");
        TH1D* histoPCMEMCALPi0InvXSectionStat               = (TH1D*)directoryPCMEMCALPi0->Get("InvCrossSectionPi0");
        cout << "Pi0 stat PCM-EMC" << endl;
        graphPCMEMCALPi0InvXSectionStat->Print();
        TGraphAsymmErrors* graphPCMEMCALPi0InvXSectionSys   = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("InvCrossSectionPi0Sys");
        cout << "Pi0 sys PCM-EMC" << endl;
        graphPCMEMCALPi0InvXSectionSys->Print();


        TH1D* histoPi0InvMassSigPlusBGPCMEMCAL[3];
        TH1D* histoPi0InvMassSigPCMEMCAL[3];
        TH1D* histoPi0InvMassSigRemBGSubPCMEMCAL[3];
        TH1D* histoPi0InvMassBGPCMEMCAL[3];
        TH1D* histoPi0InvMassRemBGPCMEMCAL[3];
        TH1D* histoPi0InvMassBGTotPCMEMCAL[3];
        TF1* fitPi0InvMassSigPCMEMCAL[3];
        TF1* fitPi0InvMassBGPCMEMCAL[3];
        Bool_t haveAllPi0InvMassPCMEMCAL[3]                 = {kFALSE, kFALSE, kFALSE};
        if (plotInvMassBins){
            for (Int_t i = 0; i < 3; i++){
                histoPi0InvMassSigPCMEMCAL[i]               = (TH1D*)directoryPCMEMCALPi0->Get(Form("Pi0_InvMassSig_Example_%s",nameTrigger[i].Data()));
                histoPi0InvMassSigPlusBGPCMEMCAL[i]         = (TH1D*)directoryPCMEMCALPi0->Get(Form("Pi0_InvMassSigPlusBG_Example_%s",nameTrigger[i].Data()));
                histoPi0InvMassBGPCMEMCAL[i]                = (TH1D*)directoryPCMEMCALPi0->Get(Form("Pi0_InvMassBG_Example_%s",nameTrigger[i].Data()));
                fitPi0InvMassSigPCMEMCAL[i]                 = (TF1*)directoryPCMEMCALPi0->Get(Form("Pi0_InvMassSigFit_Example_%s",nameTrigger[i].Data()));
                if (histoPi0InvMassSigPCMEMCAL[i] && histoPi0InvMassSigPlusBGPCMEMCAL[i] && histoPi0InvMassBGPCMEMCAL[i] && fitPi0InvMassSigPCMEMCAL[i]){
                    haveAllPi0InvMassPCMEMCAL[i]            = kTRUE;
                }


                if (haveAllPi0InvMassPCMEMCAL[i]){
                    histoPi0InvMassSigPCMEMCAL[i]->Fit(fitPi0InvMassSigPCMEMCAL[i],"QRME0");
                    for (Int_t l=0; l < 6; l++){
                        cout << fitPi0InvMassSigPCMEMCAL[i]->GetParameter(l) << "\t +- " << fitPi0InvMassSigPCMEMCAL[i]->GetParError(l) << endl;
                    }
                    fitPi0InvMassBGPCMEMCAL[i]                                  = new TF1("Linearpp","[0]+[1]*x",0.02,0.25);
                    fitPi0InvMassBGPCMEMCAL[i]->SetParameter(0, fitPi0InvMassSigPCMEMCAL[i]->GetParameter(4));
                    fitPi0InvMassBGPCMEMCAL[i]->SetParameter(1, fitPi0InvMassSigPCMEMCAL[i]->GetParameter(5));
                    TVirtualFitter * fitterPCMEMCAL                             = TVirtualFitter::GetFitter();
                    Int_t nFreeParPCMEMCAL                                      = fitPi0InvMassSigPCMEMCAL[i]->GetNumberFreeParameters();
                    double * covMatrixPCMEMCAL                                  = fitterPCMEMCAL->GetCovarianceMatrix();

                    histoPi0InvMassRemBGPCMEMCAL[i]                             = (TH1D*)histoPi0InvMassBGPCMEMCAL[i]->Clone(Form("Pi0_InvMassRemBG_Example_%s",nameTrigger[i].Data()));
                    for (Int_t j = 1; j < histoPi0InvMassRemBGPCMEMCAL[i]->GetNbinsX()+1; j++){
                        histoPi0InvMassRemBGPCMEMCAL[i]->SetBinContent(j,0);
                        histoPi0InvMassRemBGPCMEMCAL[i]->SetBinError(j,0);
                    }
                    for (Int_t j = histoPi0InvMassSigPCMEMCAL[i]->GetXaxis()->FindBin(0.01); j < histoPi0InvMassSigPCMEMCAL[i]->GetXaxis()->FindBin(0.30)+1; j++){
                        Double_t startBinEdge                                   = histoPi0InvMassSigPCMEMCAL[i]->GetXaxis()->GetBinLowEdge(j);
                        Double_t endBinEdge                                     = histoPi0InvMassSigPCMEMCAL[i]->GetXaxis()->GetBinUpEdge(j);
                        Double_t intLinearBack                                  = fitPi0InvMassBGPCMEMCAL[i]->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
                        Double_t errorLinearBck                                 = pow(( pow( (endBinEdge-startBinEdge)*fitPi0InvMassSigPCMEMCAL[i]->GetParError(4),2) +
                                                                                        pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPi0InvMassSigPCMEMCAL[i]->GetParError(5),2)
                                                                                        +2*covMatrixPCMEMCAL[nFreeParPCMEMCAL*nFreeParPCMEMCAL-2]*(endBinEdge-startBinEdge)*0.5*
                                                                                        (endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
//                         cout << j << "\t" << intLinearBack << "\t" << errorLinearBck << endl;
                        histoPi0InvMassRemBGPCMEMCAL[i]->SetBinContent(j,intLinearBack);
                        histoPi0InvMassRemBGPCMEMCAL[i]->SetBinError(j,errorLinearBck);
                    }
                    histoPi0InvMassBGTotPCMEMCAL[i]         = (TH1D*)histoPi0InvMassBGPCMEMCAL[i]->Clone(Form("Pi0_InvMassTotBG_Example_%s",nameTrigger[i].Data()));
                    histoPi0InvMassBGTotPCMEMCAL[i]->Sumw2();
                    histoPi0InvMassBGTotPCMEMCAL[i]->Add(histoPi0InvMassRemBGPCMEMCAL[i]);
                    histoPi0InvMassSigRemBGSubPCMEMCAL[i]   = (TH1D*)histoPi0InvMassSigPCMEMCAL[i]->Clone(Form("Pi0_InvMassSigRemBGSub_Example_%s",nameTrigger[i].Data()));
                    histoPi0InvMassSigRemBGSubPCMEMCAL[i]->Sumw2();
                    histoPi0InvMassSigRemBGSubPCMEMCAL[i]->Add(histoPi0InvMassRemBGPCMEMCAL[i],-1);
                    fitPi0InvMassSigPCMEMCAL[i]->SetParameter(4, 0);
                    fitPi0InvMassSigPCMEMCAL[i]->SetParameter(5, 0);
                }
                cout << nameTrigger[i].Data() << "\t" << histoPi0InvMassSigPCMEMCAL[i] << "\t" << histoPi0InvMassSigPlusBGPCMEMCAL[i] << "\t" << histoPi0InvMassBGPCMEMCAL[i] << "\t" << fitPi0InvMassSigPCMEMCAL[i]
                      << "\t" << haveAllPi0InvMassPCMEMCAL[i] << endl;
            }
        }


    //************************** Read data for EMCAL ****************************************************
    TFile* fileEMCALLow                                     = new TFile(fileNameEMCALLow.Data());
    TH1D* histoEMCALNumberOfEvents                          = (TH1D*)fileEMCALLow->Get("histoNumberOfEvents8TeV");
    TDirectory* directoryEMCALPi0                           = (TDirectory*)fileEMCALLow->Get("Pi08TeV");
        TGraphAsymmErrors* graphEMCALPi0Mass                = (TGraphAsymmErrors*)directoryEMCALPi0->Get("Pi0_Mass_data");
        graphEMCALPi0Mass                                   = ScaleGraph(graphEMCALPi0Mass, 1000.);
        TGraphAsymmErrors* graphEMCALPi0FWHM                = (TGraphAsymmErrors*)directoryEMCALPi0->Get("Pi0_Width_data");
        graphEMCALPi0FWHM                                   = ScaleGraph(graphEMCALPi0FWHM, 1000.);
        TGraphAsymmErrors* graphEMCALPi0MassMC              = (TGraphAsymmErrors*)directoryEMCALPi0->Get("Pi0_Mass_MC");
        graphEMCALPi0MassMC                                 = ScaleGraph(graphEMCALPi0MassMC, 1000.);
        TGraphAsymmErrors* graphEMCALPi0FWHMMC              = (TGraphAsymmErrors*)directoryEMCALPi0->Get("Pi0_Width_MC");
        graphEMCALPi0FWHMMC                                 = ScaleGraph(graphEMCALPi0FWHMMC, 1000.);
        TGraphAsymmErrors* graphEMCALPi0Acc                 = (TGraphAsymmErrors*)directoryEMCALPi0->Get("AcceptancePi0");
        TGraphAsymmErrors* graphEMCALPi0EffPt               = (TGraphAsymmErrors*)directoryEMCALPi0->Get("EfficiencyPi0");
        TGraphAsymmErrors* graphEMCALPi0AccTimesEff         = (TGraphAsymmErrors*)directoryEMCALPi0->Get("EffTimesAccPi0");
        TH1D* histoEMCALPi0TriggerEff[4];
        for (Int_t i = 2; i < 6; i++){
            histoEMCALPi0TriggerEff[i-2]                    = (TH1D*)directoryEMCALPi0->Get(Form("TriggerEfficiencyPi0_%s",nameTrigger[i].Data()));
        }
        TH1D* histoEMCALPi0InvXSectionStat                  = (TH1D*)directoryEMCALPi0->Get("InvCrossSectionPi0");
        TGraphAsymmErrors* graphEMCALPi0InvXSectionStat     = (TGraphAsymmErrors*)directoryEMCALPi0->Get("graphInvCrossSectionPi0");
        cout << "Pi0 stat EMC-EMC" << endl;
        graphEMCALPi0InvXSectionStat->Print();

        TGraphAsymmErrors* graphEMCALPi0InvXSectionSys      = (TGraphAsymmErrors*)directoryEMCALPi0->Get("InvCrossSectionPi0Sys");
        cout << "Pi0 sys EMC-EMC" << endl;
        graphEMCALPi0InvXSectionSys->Print();


        TH1D* histoPi0InvMassSigPlusBGEMCAL[3];
        TH1D* histoPi0InvMassSigEMCAL[3];
        TH1D* histoPi0InvMassSigRemBGSubEMCAL[3];
        TH1D* histoPi0InvMassBGEMCAL[3];
        TH1D* histoPi0InvMassRemBGEMCAL[3];
        TH1D* histoPi0InvMassBGTotEMCAL[3];
        TF1* fitPi0InvMassSigEMCAL[3];
        TF1* fitPi0InvMassBGEMCAL[3];
        Bool_t haveAllPi0InvMassEMCAL[3]                 = {kFALSE, kFALSE, kFALSE};
        if (plotInvMassBins){
            for (Int_t i = 0; i < 3; i++){
                histoPi0InvMassSigEMCAL[i]                  = (TH1D*)directoryEMCALPi0->Get(Form("Pi0_InvMassSig_Example_%s",nameTrigger[i].Data()));
                histoPi0InvMassSigPlusBGEMCAL[i]            = (TH1D*)directoryEMCALPi0->Get(Form("Pi0_InvMassSigPlusBG_Example_%s",nameTrigger[i].Data()));
                histoPi0InvMassBGEMCAL[i]                   = (TH1D*)directoryEMCALPi0->Get(Form("Pi0_InvMassBG_Example_%s",nameTrigger[i].Data()));
                fitPi0InvMassSigEMCAL[i]                    = (TF1*)directoryEMCALPi0->Get(Form("Pi0_InvMassSigFit_Example_%s",nameTrigger[i].Data()));
                if (histoPi0InvMassSigEMCAL[i] && histoPi0InvMassSigPlusBGEMCAL[i] && histoPi0InvMassBGEMCAL[i] && fitPi0InvMassSigEMCAL[i]){
                    haveAllPi0InvMassEMCAL[i]               = kTRUE;
                }
                if (haveAllPi0InvMassEMCAL[i]){
                    histoPi0InvMassSigEMCAL[i]->Fit(fitPi0InvMassSigEMCAL[i],"QRME0");
                    for (Int_t l=0; l < 6; l++){
                        cout << fitPi0InvMassSigEMCAL[i]->GetParameter(l) << "\t +- " << fitPi0InvMassSigEMCAL[i]->GetParError(l) << endl;
                    }
                    fitPi0InvMassBGEMCAL[i]                                  = new TF1("Linearpp","[0]+[1]*x",0.08,0.30);
                    fitPi0InvMassBGEMCAL[i]->SetParameter(0, fitPi0InvMassSigEMCAL[i]->GetParameter(4));
                    fitPi0InvMassBGEMCAL[i]->SetParameter(1, fitPi0InvMassSigEMCAL[i]->GetParameter(5));
                    TVirtualFitter * fitterEMCAL                             = TVirtualFitter::GetFitter();
                    Int_t nFreeParEMCAL                                      = fitPi0InvMassSigEMCAL[i]->GetNumberFreeParameters();
                    double * covMatrixEMCAL                                  = fitterEMCAL->GetCovarianceMatrix();

                    histoPi0InvMassRemBGEMCAL[i]                             = (TH1D*)histoPi0InvMassBGEMCAL[i]->Clone(Form("Pi0_InvMassRemBG_Example_%s",nameTrigger[i].Data()));
                    for (Int_t j = 1; j < histoPi0InvMassRemBGEMCAL[i]->GetNbinsX()+1; j++){
                        histoPi0InvMassRemBGEMCAL[i]->SetBinContent(j,0);
                        histoPi0InvMassRemBGEMCAL[i]->SetBinError(j,0);
                    }
                    for (Int_t j = histoPi0InvMassSigEMCAL[i]->GetXaxis()->FindBin(0.01); j < histoPi0InvMassSigEMCAL[i]->GetXaxis()->FindBin(0.30)+1; j++){
                        Double_t startBinEdge                                   = histoPi0InvMassSigEMCAL[i]->GetXaxis()->GetBinLowEdge(j);
                        Double_t endBinEdge                                     = histoPi0InvMassSigEMCAL[i]->GetXaxis()->GetBinUpEdge(j);
                        Double_t intLinearBack                                  = fitPi0InvMassBGEMCAL[i]->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
                        Double_t errorLinearBck                                 = pow(( pow( (endBinEdge-startBinEdge)*fitPi0InvMassSigEMCAL[i]->GetParError(4),2) +
                                                                                        pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPi0InvMassSigEMCAL[i]->GetParError(5),2)
                                                                                        +2*covMatrixEMCAL[nFreeParEMCAL*nFreeParEMCAL-2]*(endBinEdge-startBinEdge)*0.5*
                                                                                        (endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
//                         cout << j << "\t" << intLinearBack << "\t" << errorLinearBck << endl;
                        histoPi0InvMassRemBGEMCAL[i]->SetBinContent(j,intLinearBack);
                        histoPi0InvMassRemBGEMCAL[i]->SetBinError(j,errorLinearBck);
                    }
                    histoPi0InvMassBGTotEMCAL[i]         = (TH1D*)histoPi0InvMassBGEMCAL[i]->Clone(Form("Pi0_InvMassTotBG_Example_%s",nameTrigger[i].Data()));
                    histoPi0InvMassBGTotEMCAL[i]->Sumw2();
                    histoPi0InvMassBGTotEMCAL[i]->Add(histoPi0InvMassRemBGEMCAL[i]);
                    histoPi0InvMassSigRemBGSubEMCAL[i]   = (TH1D*)histoPi0InvMassSigEMCAL[i]->Clone(Form("Pi0_InvMassSigRemBGSub_Example_%s",nameTrigger[i].Data()));
                    histoPi0InvMassSigRemBGSubEMCAL[i]->Sumw2();
                    histoPi0InvMassSigRemBGSubEMCAL[i]->Add(histoPi0InvMassRemBGEMCAL[i],-1);
                    fitPi0InvMassSigEMCAL[i]->SetParameter(4, 0);
                    fitPi0InvMassSigEMCAL[i]->SetParameter(5, 0);
                }
                cout << nameTrigger[i].Data() << "\t" << histoPi0InvMassSigEMCAL[i] << "\t" << histoPi0InvMassSigPlusBGEMCAL[i] << "\t" << histoPi0InvMassBGEMCAL[i] << "\t" << fitPi0InvMassSigEMCAL[i]
                      << "\t" << haveAllPi0InvMassEMCAL[i] << endl;
            }
        }



    //************************** Read data for EMCAL merged **************************************************
    TFile* fileEMCALmerged                                  = new TFile(fileNameEMCALmerged.Data());
    TDirectory* directoryEMCALmergedPi0                     = (TDirectory*)fileEMCALmerged->Get("Pi08TeV");
         TH1D* histoEMCALMergedPi0Eff                        = NULL;
         //TH1D* histoEMCALMergedPi0Pur                        = NULL;
        TH1D* histoEMCALMergedPi0AccTimesEff                = NULL;
        TH1D* histoEMCALMergedPi0InvXSectionStat            = NULL;
        TH1D* histoEMCALMergedPi0InvXSectionSys             = NULL;
        TGraphAsymmErrors* graphEMCALMergedPi0AccTimesEff       = NULL;
        TGraphAsymmErrors* graphEMCALMergedPi0Purity            = NULL;
        //TGraphAsymmErrors* graphEMCALMergedPi0AccTimesEffDivPur = NULL;
        TGraphAsymmErrors* graphEMCALMergedPi0InvXSectionStat   = NULL;
        TGraphAsymmErrors* graphEMCALMergedPi0InvXSectionSys    = NULL;
             histoEMCALMergedPi0Eff                          = (TH1D*)directoryEMCALmergedPi0->Get("EfficiencyPi0");
             //histoEMCALMergedPi0Pur                          = (TH1D*)directoryEMCALmergedPi0->Get("PurityPi0");
            graphEMCALMergedPi0AccTimesEff                  = (TGraphAsymmErrors*)directoryEMCALmergedPi0->Get("EffTimesAccPi0_EGA");
            //graphEMCALMergedPi0Purity                       = (TGraphAsymmErrors*)directoryEMCALmergedPi0->Get("PurityPi0");
            //graphEMCALMergedPi0AccTimesEffDivPur            = CalculateGraphErrRatioToOtherTGraphErr(graphEMCALMergedPi0AccTimesEff, graphEMCALMergedPi0Purity,kTRUE);
            graphEMCALMergedPi0InvXSectionStat              = (TGraphAsymmErrors*)directoryEMCALmergedPi0->Get("graphInvCrossSectionPi0");
            graphEMCALMergedPi0InvXSectionSys               = (TGraphAsymmErrors*)directoryEMCALmergedPi0->Get("InvCrossSectionPi0Sys");
                cout << "EMCAL merged stat" << endl;
                for (Int_t i = 0; graphEMCALMergedPi0InvXSectionStat->GetX()[0]< 14; i++){
                    graphEMCALMergedPi0InvXSectionStat->RemovePoint(0);
                }
                graphEMCALMergedPi0InvXSectionStat->Print();
                cout << "EMCAL merged sys" << endl;
                for (Int_t i = 0; graphEMCALMergedPi0InvXSectionSys->GetX()[0]< 14; i++){
                    graphEMCALMergedPi0InvXSectionSys->RemovePoint(0);
                }
                Double_t * yValue               = graphEMCALMergedPi0InvXSectionSys->GetY();
                Double_t* yErrorLow             = graphEMCALMergedPi0InvXSectionSys->GetEYlow();
                Double_t* yErrorHigh            = graphEMCALMergedPi0InvXSectionSys->GetEYhigh();
                Int_t nPoints                   = graphEMCALMergedPi0InvXSectionSys->GetN();
                for (Int_t i = 0; i < nPoints; i++){
                    yErrorLow[i]           = yValue[i]*0.2;
                    yErrorHigh[i]          = yValue[i]*0.2;
                }
                graphEMCALMergedPi0InvXSectionSys->Print();
//         return;

    //************************** Read data for PHOS *****************************************************

        TFile* filePHOS                                     = new TFile(fileNamePHOS.Data());
        TDirectory* directoryPHOSPi0                           = (TDirectory*)filePHOS->Get("Pi08TeV");
            TH1D* histoPHOSPi0Mass                = (TH1D*)directoryPHOSPi0->Get("MassPi0");
            histoPHOSPi0Mass->Scale(1000.);
            //for(Int_t k=0; k<6; k++) histoPHOSPi0Mass->SetBinContent(k,0.);
            //for(Int_t k=histoPHOSPi0Mass->GetNbinsX(), j=k-6; k>=j; k--) histoPHOSPi0Mass->SetBinContent(k,0.);

            TH1D* histoPHOSPi0FWHMMeV                = (TH1D*)directoryPHOSPi0->Get("FWHMPi0MeV");
            //for(Int_t k=0; k<6; k++) histoPHOSPi0FWHMMeV->SetBinContent(k,-1000.);
           // for(Int_t k=histoPHOSPi0FWHMMeV->GetNbinsX(), j=k-6; k>=j; k--) histoPHOSPi0FWHMMeV->SetBinContent(k,-1000.);

            TH1D* histoPHOSPi0TrueMass              = (TH1D*)directoryPHOSPi0->Get("TrueMassPi0");
            histoPHOSPi0TrueMass->Scale(1000.);
           // for(Int_t k=0; k<1; k++) histoPHOSPi0TrueMass->SetBinContent(k,0.);
           // for(Int_t k=histoPHOSPi0TrueMass->GetNbinsX(), j=k-5; k>=j; k--) histoPHOSPi0TrueMass->SetBinContent(k,0.);

            TH1D* histoPHOSPi0TrueFWHMMeV              = (TH1D*)directoryPHOSPi0->Get("TrueFWHMPi0MeV");
           // for(Int_t k=0; k<1; k++) histoPHOSPi0TrueFWHMMeV->SetBinContent(k,-1000.);
           // for(Int_t k=histoPHOSPi0TrueFWHMMeV->GetNbinsX(), j=k-5; k>=j; k--) histoPHOSPi0TrueFWHMMeV->SetBinContent(k,-1000.);

            TH1D* histoPHOSPi0Acc                 = (TH1D*)directoryPHOSPi0->Get("AcceptancePi0");
            //for(Int_t k=histoPHOSPi0Acc->GetNbinsX(), j=k-5; k>=j; k--) histoPHOSPi0Acc->SetBinContent(k,0.);

            TH1D* histoPHOSPi0TrueEffPt               = (TH1D*)directoryPHOSPi0->Get("EfficiencyPi0");
           // for(Int_t k=histoPHOSPi0TrueEffPt->GetNbinsX(), j=k-5; k>=j; k--) histoPHOSPi0TrueEffPt->SetBinContent(k,0.);

            TH1D* histoPHOSPi0AccTimesEff      = (TH1D*)histoPHOSPi0TrueEffPt->Clone("histoPHOSPi0AccTimesEff");
            histoPHOSPi0AccTimesEff->Multiply(histoPHOSPi0Acc);
            histoPHOSPi0AccTimesEff->Scale(2*TMath::Pi());
            //TH1D* graphPHOSPi0AccTimesEff         = (TH1D*)directoryPHOSPi0->Get("EffTimesAccPi0_MB");

            TH1D* histoPHOSPi0InvXSectionStat                   = (TH1D*)directoryPHOSPi0->Get("InvCrossSectionPi0");
            TGraphAsymmErrors* graphPHOSPi0InvXSectionStat      = new TGraphAsymmErrors(histoPHOSPi0InvXSectionStat);
            cout << "Pi0 stat PHOS" << endl;
            graphPHOSPi0InvXSectionStat->Print();

            TGraphErrors* dummySysGraphPHOS                     = (TGraphErrors*)directoryPHOSPi0->Get("InvCrossSectionPi0Sys");
            Double_t* xValuePHOS                                = dummySysGraphPHOS->GetX();
            Double_t* yValuePHOS                                = dummySysGraphPHOS->GetY();
            Double_t* xErrPHOS                                  = dummySysGraphPHOS->GetEX();
            Double_t* yErrPHOS                                  = dummySysGraphPHOS->GetEY();

            TGraphAsymmErrors* graphPHOSPi0InvXSectionSys       = new TGraphAsymmErrors(dummySysGraphPHOS->GetN(), xValuePHOS, yValuePHOS, xErrPHOS, xErrPHOS, yErrPHOS, yErrPHOS);
            cout << "Pi0 sys PHOS" << endl;
            graphPHOSPi0InvXSectionSys->Print();

        TFile* filePHOSprelim                                     = new TFile(fileNamePHOSprelim.Data());

    // *******************************************************************************************************
    // ************************** Loading theory calculations ************************************************
    // *******************************************************************************************************
    TFile* fileTheoryCompilation                            = new TFile(fileNameTheory.Data());
    TFile* fileTheoryPythia8                                = new TFile(fileNameTheoryPythia8.Data());
    TFile* fileTheoryPythia8_4C                             = new TFile(fileNameTheoryPythia8_4C.Data());
    TFile* fileTheorypp8TeVPi0DSS14                         = new TFile(fileNameTheoryPi0DSS14.Data());
    TFile* fileEtaToPi07000GeV                              = new TFile(fileNameEtaToPi07TeV.Data());
    TFile* fileEtaToPi02760GeV                              = new TFile(fileNameEtaToPi02760GeV.Data());
    TDirectory* directoryfileEtaToPi02760GeV                = (TDirectory*)fileEtaToPi02760GeV->Get("Eta2.76TeV");
        TH1F* histoPythia8InvXSection                       = (TH1F*) fileTheoryPythia8->Get("fHistInvXsec_Pi0");
        TH1F* histoPythia8InvXSectionEta                    = (TH1F*) fileTheoryPythia8->Get("fHistInvXsec_Eta");
        TH1F* histoPythia8_4CInvXSection                    = (TH1F*) fileTheoryPythia8_4C->Get("fHistInvXsec_Pi0");
        TH1F* histoPythia8_4CInvXSectionEta                 = (TH1F*) fileTheoryPythia8_4C->Get("fHistInvXsec_Eta");
        TGraphAsymmErrors* graphPi0DSS14                    = (TGraphAsymmErrors*)fileTheorypp8TeVPi0DSS14->Get("fGraphInvXsec_Pi0");
        for (int i=0;i<graphPi0DSS14->GetN();i++) graphPi0DSS14->GetY()[i] /= 2.;
        TGraph* graphEtaToPi07000GeV                        = (TGraph*) fileEtaToPi07000GeV->Get("EtaPi0Ratio7000GeV");
        TGraph* graphEtaToPi02760GeV                        = (TGraph*) directoryfileEtaToPi02760GeV->Get("graphRatioEtaToPi0Comb2760GeVTotErr");
//        TH1F* histoPythia8InvXSection_VarBinning            = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Spec8000GeVVarBinning");
        TGraph* graphNLOCalcPi0MuHalf                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuHalf8000GeV");
        TGraph* graphNLOCalcPi0MuOne                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuOne8000GeV");
        TGraph* graphNLOCalcPi0MuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuTwo8000GeV");
        TGraph* graphNLOCalcEtaMuHalf                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuHalf8000GeV");
        TGraph* graphNLOCalcEtaMuOne                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuOne8000GeV");
        TGraph* graphNLOCalcEtaMuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuTwo8000GeV");
        TGraph* graphNLOEtaToPi0MuHalf                      = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuHalf8000GeV");
        TGraph* graphNLOEtaToPi0MuOne                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuOne8000GeV");
        TGraph* graphNLOEtaToPi0MuTwo                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuTwo8000GeV");
        //TGraph* graphNLODSSCalcMuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcDSSInvSecPi0MuTwo8000GeV");
        //TGraph* graphNLOCGCCalcMuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcCGCInvCrossSec8000GeV");
        //TGraphAsymmErrors* graphNLODSS14Calc                = (TGraphAsymmErrors*)fileTheoryCompilation->Get("graphNLOCalcDSS14InvCrossSec8000GeV");

        while (graphNLOCalcEtaMuHalf->GetX()[graphNLOCalcEtaMuHalf->GetN()-1] > 37. )
            graphNLOCalcEtaMuHalf->RemovePoint(graphNLOCalcEtaMuHalf->GetN()-1);
        while (graphNLOCalcEtaMuOne->GetX()[graphNLOCalcEtaMuOne->GetN()-1] > 37. )
            graphNLOCalcEtaMuOne->RemovePoint(graphNLOCalcEtaMuOne->GetN()-1);
        while (graphNLOCalcEtaMuTwo->GetX()[graphNLOCalcEtaMuTwo->GetN()-1] > 37. )
            graphNLOCalcEtaMuTwo->RemovePoint(graphNLOCalcEtaMuTwo->GetN()-1);

        cout << "mu half" << endl;
        graphNLOEtaToPi0MuHalf->Print();
        cout << "mu one" << endl;
        graphNLOEtaToPi0MuOne->Print();
        cout << "mu two" << endl;
        graphNLOEtaToPi0MuTwo->Print();

    // *******************************************************************************************************
    // ************************** Loading charged pion results ***********************************************
    // *******************************************************************************************************
//    TFile* fileChargedPionInputpp                           = new TFile(fileNameChargedPionPP.Data());
//        TH1D* histoChPiInvYieldPubStatPP                    = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecPubStat7TeV");
//        TH1D* histoChPiInvYieldPubSystPP                    = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecPubSyst7TeV");

//        TH1D* histoChPiInvYieldHighPtStatPP                 = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecHighPtStat7TeV");
//        TH1D* histoChPiInvYieldHighPtSystPP                 = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecHighPtSyst7TeV");
//        TH1D* histoChPiInvYieldLowPtStatCMS                 = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtStat7TeVCMS");
//        TH1D* histoChPiInvYieldLowPtSysCMS                  = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtSys7TeVCMS");
//        TH1D* histoChPiInvYieldLowPtStatPP                  = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtStatPP7TeV");
//        TH1D* histoChPiInvYieldLowPtSysPP                   = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtSysPP7TeV");

    // *******************************************************************************************************
    // ************************** Loading charged hadron results ***********************************************
    // *******************************************************************************************************
//    TFile* fileChargedHadronsInputpp                        = new TFile(fileNameChargedHadronPP.Data());
//        TGraphAsymmErrors* graphChargedHadronsStatPP        = (TGraphAsymmErrors*)fileChargedHadronsInputpp->Get("graphChargedHadronsStatPP8TeV");
//        graphChargedHadronsStatPP                           = ScaleGraph(graphChargedHadronsStatPP,0.5*1e9/xSection8TeVINEL);
//        TGraphAsymmErrors* graphChargedHadronsSysPP         = (TGraphAsymmErrors*)fileChargedHadronsInputpp->Get("graphChargedHadronsSysPP8TeV");
//        graphChargedHadronsSysPP                            = ScaleGraph(graphChargedHadronsSysPP,0.5*1e9/xSection8TeVINEL);

    // *******************************************************************************************************
    // ************************** Combination of different pi0 measurements **********************************
    // *******************************************************************************************************
    // REMARKS:
    //      - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //      - currently only PCM-EMCAL vs others fully implemeted energy independent
    //      - extendable to other energies
    //     - offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented


    // definition of array of histograms (NULL - means we have no measurement at this energy for this rec-method)
    // for statistical error and final value from respective method
    TH1D* statErrorCollectionPi0[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorCollectionPi0[i]   = NULL;
    }
    statErrorCollectionPi0[0]       = (TH1D*)histoPCMPi0InvXSectionStat->Clone("statErrPCMPi0");
    statErrorCollectionPi0[1]       = (TH1D*)histoPHOSPi0InvXSectionStat->Clone("statErrPHOSPi0");
    statErrorCollectionPi0[2]       = (TH1D*)histoEMCALPi0InvXSectionStat->Clone("statErrEMCALPi0");
    statErrorCollectionPi0[4]       = (TH1D*)histoPCMEMCALPi0InvXSectionStat->Clone("statErrPCMEMCALPi0");

    // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
    // for systematic error from respective method
    TGraphAsymmErrors* sysErrorCollectionPi0[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorCollectionPi0[i]    = NULL;
    }
    sysErrorCollectionPi0[0]        = (TGraphAsymmErrors*)graphPCMPi0InvXSectionSys->Clone("sysErrPCMPi0");
    sysErrorCollectionPi0[1]        = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionSys->Clone("sysErrPHOSPi0");
    sysErrorCollectionPi0[2]        = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionSys->Clone("sysErrEMCALPi0");
    sysErrorCollectionPi0[4]        = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionSys->Clone("sysErrPCMEMCALPi0");


    // **********************************************************************************************************************
    // ******************************************* Adding EMCal & EMC merged  ***********************************************
    // **********************************************************************************************************************

    // Definition of final pt binning (has to be set manually)
    Double_t xPtLimitsPi0[100];
    Int_t maxNBinsPi0               = 0;

    maxNBinsPi0                  = GetBinning( xPtLimitsPi0, "Pi0", "8TeV", 11 );
//    maxNBinsPi0--;

    Double_t xPtLimitsPi0WOMerged[70];
    Int_t maxNBinsPi0W0Merged       = GetBinning( xPtLimitsPi0WOMerged, "Pi0", "8TeV", 11 );

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsPi0[11]            = { 0,  6,  0,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t offSetsPi0Sys[11]         = { 1,  6,  7,  0,  5,
                                        0,  0,  0,  0,  36,
                                        0};
    Int_t offSetPi0Shifting[11]     = { 0,  0,  0,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0 };
    Int_t nComBinsPi0Shifting[11]   = { 0,  0,  0,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0 };

    TH1D* histoEMCALMergedPi0InvXSectionStatCorrBin     = new TH1D("histoEMCALMergedPi0InvXSectionStatCorrBin", "", maxNBinsPi0, xPtLimitsPi0);
    Int_t firstBinMergedPi0 = 34;
    cout << graphEMCALMergedPi0InvXSectionStat->GetX()[0] << endl;

    while (histoEMCALMergedPi0InvXSectionStatCorrBin->GetBinCenter(firstBinMergedPi0) < graphEMCALMergedPi0InvXSectionStat->GetX()[0]){
        cout << histoEMCALMergedPi0InvXSectionStatCorrBin->GetBinCenter(firstBinMergedPi0) << endl;
        histoEMCALMergedPi0InvXSectionStatCorrBin->SetBinContent(firstBinMergedPi0, 0);
        histoEMCALMergedPi0InvXSectionStatCorrBin->SetBinError(firstBinMergedPi0, 0);
        firstBinMergedPi0++;
    }
    for (Int_t i = 0; i < graphEMCALMergedPi0InvXSectionStat->GetN(); i++){
        histoEMCALMergedPi0InvXSectionStatCorrBin->SetBinContent(i+firstBinMergedPi0, graphEMCALMergedPi0InvXSectionStat->GetY()[i]);
        histoEMCALMergedPi0InvXSectionStatCorrBin->SetBinError(i+firstBinMergedPi0, graphEMCALMergedPi0InvXSectionStat->GetEYlow()[i]);
    }

    graphEMCALMergedPi0InvXSectionStat->Print();
    for (Int_t i = 1; i < histoEMCALMergedPi0InvXSectionStatCorrBin->GetNbinsX()+1; i++){
        cout << "Bin " << i << "\t" <<  histoEMCALMergedPi0InvXSectionStatCorrBin->GetBinCenter(i) << "\t" << histoEMCALMergedPi0InvXSectionStatCorrBin->GetBinContent(i)
             << "\t" << histoEMCALMergedPi0InvXSectionStatCorrBin->GetBinError(i) << endl;
    }
     //return;

    statErrorCollectionPi0[9]          = (TH1D*)histoEMCALMergedPi0InvXSectionStatCorrBin->Clone("statErrEMCALMergedPi0");
    sysErrorCollectionPi0[9]           = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionSys->Clone("sysErrEMCALMergedPi0");


    // **********************************************************************************************************************
    // ************************ Adding PCM-EMCal measurements to matrix & calculating relativ errors ************************
    // **********************************************************************************************************************


    TGraphAsymmErrors* statErrorRelCollectionPi0[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorRelCollectionPi0[i]        = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (statErrorCollectionPi0[i]){
            statErrorRelCollectionPi0[i]    = new TGraphAsymmErrors(statErrorCollectionPi0[i]);
            statErrorRelCollectionPi0[i]    = CalculateRelErrUpAsymmGraph( statErrorRelCollectionPi0[i], Form("relativeStatErrorPi0_%s", nameMeasGlobal[i].Data()));
        }
    }

    TGraphAsymmErrors* sysErrorRelCollectionPi0[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorRelCollectionPi0[i]         = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (sysErrorCollectionPi0[i]){
          sysErrorRelCollectionPi0[i]       = CalculateRelErrUpAsymmGraph( sysErrorCollectionPi0[i], Form("relativeSysErrorPi0_%s", nameMeasGlobal[i].Data()));
        }
    }

//     return;

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************

    TGraph* graphWeightsPi0A[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsPi0A[i]         = NULL;
    }

    // Declaration & calculation of combined spectrum
    TString fileNamePi0OutputWeightingA                  = Form("%s/Pi0_WeightingMethodA.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombPi0InvXSectionStatA      = NULL;
    TGraphAsymmErrors* graphCombPi0InvXSectionSysA       = NULL;
    TGraphAsymmErrors* graphCombPi0InvXSectionTotA       = CombinePtPointsSpectraFullCorrMat(   statErrorCollectionPi0,    sysErrorCollectionPi0,
                                                                                                xPtLimitsPi0, maxNBinsPi0,
                                                                                                offSetsPi0, offSetsPi0Sys,
                                                                                                graphCombPi0InvXSectionStatA, graphCombPi0InvXSectionSysA,
                                                                                                fileNamePi0OutputWeightingA, "8TeV", "Pi0", kTRUE,
                                                                                                0x0, fileInputCorrFactors
                                                                                            );


    if (graphCombPi0InvXSectionTotA == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }
    // remove bin from 0-0.3 (should have been done automatically in principle)
//    graphCombPi0InvXSectionStatA->RemovePoint(0);
//    graphCombPi0InvXSectionSysA->RemovePoint(0);
//    graphCombPi0InvXSectionTotA->RemovePoint(0);
    cout << __LINE__ << endl;
    graphCombPi0InvXSectionTotA->Print();
//     return;

    // Reading weights from output file for plotting
    ifstream fileWeightsPi0ReadA;
    fileWeightsPi0ReadA.open(fileNamePi0OutputWeightingA,ios_base::in);
    cout << "reading" << fileNamePi0OutputWeightingA << endl;
    Double_t xValuesPi0ReadA[70];
    Double_t weightsPi0ReadA[11][70];
    Int_t availablePi0MeasA[11]        = { -1, -1, -1, -1, -1,
                                        -1, -1, -1, -1, -1,
                                        -1};
    Int_t nMeasSetPi0A                 = 5;
    Int_t nPtBinsPi0ReadA              = 0;
    while(!fileWeightsPi0ReadA.eof() && nPtBinsPi0ReadA < 70){
        TString garbage             = "";
        if (nPtBinsPi0ReadA == 0){
            fileWeightsPi0ReadA >> garbage ;//>> availablePi0Meas[0] >> availablePi0Meas[1] >> availablePi0Meas[2] >> availablePi0Meas[3];
            for (Int_t i = 0; i < nMeasSetPi0A; i++){
                fileWeightsPi0ReadA >> availablePi0MeasA[i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availablePi0MeasA[i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsPi0ReadA >> xValuesPi0ReadA[nPtBinsPi0ReadA-1];
            for (Int_t i = 0; i < nMeasSetPi0A; i++){
                fileWeightsPi0ReadA >> weightsPi0ReadA[availablePi0MeasA[i]][nPtBinsPi0ReadA-1] ;
            }
            cout << "read: "<<  nPtBinsPi0ReadA << "\t"<< xValuesPi0ReadA[nPtBinsPi0ReadA-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetPi0A; i++){
                cout << weightsPi0ReadA[availablePi0MeasA[i]][nPtBinsPi0ReadA-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsPi0ReadA++;
    }
    nPtBinsPi0ReadA                    = nPtBinsPi0ReadA-2 ;
    fileWeightsPi0ReadA.close();

    for (Int_t i = 0; i < nMeasSetPi0A; i++){
        graphWeightsPi0A[availablePi0MeasA[i]]                        = new TGraph(nPtBinsPi0ReadA,xValuesPi0ReadA,weightsPi0ReadA[availablePi0MeasA[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsPi0ReadA; n++){
            if (graphWeightsPi0A[availablePi0MeasA[i]]->GetY()[bin] == 0) graphWeightsPi0A[availablePi0MeasA[i]]->RemovePoint(bin);
            else bin++;
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Plotting weights Method A ************************************************
    // **********************************************************************************************************************
    Int_t textSizeLabelsPixel           = 900*0.04;

    TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
    canvasWeights->SetLogx();

    TH2F * histo2DPi0Weights;
    histo2DPi0Weights = new TH2F("histo2DPi0Weights","histo2DPi0Weights",11000,0.23,110.,1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DPi0Weights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DPi0Weights->GetXaxis()->SetMoreLogLabels();
    histo2DPi0Weights->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DPi0Weights->GetYaxis()->SetRangeUser(-10,10);
    histo2DPi0Weights->Draw("copy");

        TLatex *labelWeightsEnergy      = new TLatex(0.7,0.20,collisionSystem8TeV.Data());
        SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
        labelWeightsEnergy->SetTextFont(43);
        labelWeightsEnergy->Draw();
        TLatex *labelWeightsPi0         = new TLatex(0.7,0.16,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelWeightsPi0, 0.85*textSizeLabelsPixel,4);
        labelWeightsPi0->SetTextFont(43);
        labelWeightsPi0->Draw();

        TLegend* legendWeights   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetPi0A), 32);
        for (Int_t i = 0; i < nMeasSetPi0A; i++){
            DrawGammaSetMarkerTGraph(graphWeightsPi0A[availablePi0MeasA[i]], markerStyleDet[availablePi0MeasA[i]], markerSizeDet[availablePi0MeasA[i]]*0.5, colorDet[availablePi0MeasA[i]] , colorDet[availablePi0MeasA[i]]);
            graphWeightsPi0A[availablePi0MeasA[i]]->Draw("p,same,z");
            legendWeights->AddEntry(graphWeightsPi0A[availablePi0MeasA[i]],nameMeasGlobal[availablePi0MeasA[i]],"p");
        }
        legendWeights->Draw();

        labelWeightsEnergy->Draw();
        labelWeightsPi0->Draw();

//      DrawGammaLines(0.23, 110. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.23, 110. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.23, 110. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.23, 110. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.23, 110. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/Pi0_WeightsA.%s",outputDir.Data(),suffix.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErr->SetLogx();

    TH2F * histo2DRelSysErr;
    histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.23,110.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErr->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelSysErr->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelSysErr->Draw("copy");

        TLegend* legendRelSysErr        = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetPi0A), 0.95, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetPi0A; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionPi0[availablePi0MeasA[i]], markerStyleDet[availablePi0MeasA[i]], markerSizeDet[availablePi0MeasA[i]]*0.5, colorDet[availablePi0MeasA[i]],
                                     colorDet[availablePi0MeasA[i]]);
            sysErrorRelCollectionPi0[availablePi0MeasA[i]]->Draw("p,same,z");
            legendRelSysErr->AddEntry(sysErrorRelCollectionPi0[availablePi0MeasA[i]],nameMeasGlobal[availablePi0MeasA[i]],"p");
        }
        legendRelSysErr->Draw();

        TLatex *labelRelSysErrEnergy    = new TLatex(0.15,0.89,collisionSystem8TeV.Data());
        SetStyleTLatex( labelRelSysErrEnergy, 0.85*textSizeLabelsPixel,4);
        labelRelSysErrEnergy->SetTextFont(43);
        labelRelSysErrEnergy->Draw();
        TLatex *labelRelSysErrPi0       = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelSysErrPi0, 0.85*textSizeLabelsPixel,4);
        labelRelSysErrPi0->SetTextFont(43);
        labelRelSysErrPi0->Draw();

    canvasRelSysErr->SaveAs(Form("%s/Pi0_RelSysErr.%s",outputDir.Data(),suffix.Data()));

    histo2DRelSysErr->GetYaxis()->SetRangeUser(0,30.5);
    histo2DRelSysErr->Draw("copy");

        for (Int_t i = 0; i < nMeasSetPi0A; i++){
            sysErrorRelCollectionPi0[availablePi0MeasA[i]]->Draw("p,same,z");
        }
        legendRelSysErr->Draw();

        labelRelSysErrEnergy->Draw();
        labelRelSysErrPi0->Draw();

    canvasRelSysErr->SaveAs(Form("%s/Pi0_RelSysErrZoomed.%s",outputDir.Data(),suffix.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelStatErr->SetLogx();

    TH2F * histo2DRelStatErr;
    histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.23,110.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErr->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelStatErr->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelStatErr->Draw("copy");
        TLegend* legendRelStatErr       = GetAndSetLegend2(0.14, 0.92-(0.035*nMeasSetPi0A), 0.45, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetPi0A; i++){
            DrawGammaSetMarkerTGraph(statErrorRelCollectionPi0[availablePi0MeasA[i]], markerStyleDet[availablePi0MeasA[i]], markerSizeDet[availablePi0MeasA[i]]*0.5, colorDet[availablePi0MeasA[i]],
                                     colorDet[availablePi0MeasA[i]]);
            statErrorRelCollectionPi0[availablePi0MeasA[i]]->Draw("p,same,z");
            legendRelStatErr->AddEntry(statErrorRelCollectionPi0[availablePi0MeasA[i]],nameMeasGlobal[availablePi0MeasA[i]],"p");
        }
        legendRelStatErr->Draw();

        TLatex *labelRelStatErrEnergy   = new TLatex(0.75,0.89,collisionSystem8TeV.Data());
        SetStyleTLatex( labelRelStatErrEnergy, 0.85*textSizeLabelsPixel,4);
        labelRelStatErrEnergy->SetTextFont(43);
        labelRelStatErrEnergy->Draw();
        TLatex *labelRelStatErrPi0      = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelStatErrPi0, 0.85*textSizeLabelsPixel,4);
        labelRelStatErrPi0->SetTextFont(43);
        labelRelStatErrPi0->Draw();

    canvasRelStatErr->SaveAs(Form("%s/Pi0_RelStatErr.%s",outputDir.Data(),suffix.Data()));

    histo2DRelStatErr->GetYaxis()->SetRangeUser(0,30.5);
    histo2DRelStatErr->Draw("copy");
        for (Int_t i = 0; i < nMeasSetPi0A; i++){
            statErrorRelCollectionPi0[availablePi0MeasA[i]]->Draw("p,same,z");
        }
        legendRelStatErr->Draw();

        labelRelStatErrEnergy->Draw();
        labelRelStatErrPi0->Draw();

    canvasRelStatErr->SaveAs(Form("%s/Pi0_RelStatErrZoomed.%s",outputDir.Data(),suffix.Data()));

    //  *********************************************************************************************************************
    //  ************************ Visualize relative total errors of different combination methods Pi0 ***********************
    //  *********************************************************************************************************************
    TGraphAsymmErrors* graphCombPi0InvXSectionRelStatA       = CalculateRelErrUpAsymmGraph( graphCombPi0InvXSectionStatA, "relativeStatErrorPi0_MethodA");
    TGraphAsymmErrors* graphCombPi0InvXSectionRelSysA        = CalculateRelErrUpAsymmGraph( graphCombPi0InvXSectionSysA, "relativeSysErrorPi0_MethodA");
    TGraphAsymmErrors* graphCombPi0InvXSectionRelTotA        = CalculateRelErrUpAsymmGraph( graphCombPi0InvXSectionTotA, "relativeTotalErrorPi0_MethodA");


    TCanvas* canvasRelTotErr            = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelTotErr->SetLogx();

    TH2F * histo2DRelTotErrPi0;
    histo2DRelTotErrPi0                 = new TH2F("histo2DRelTotErrPi0","histo2DRelTotErrPi0",11000,0.23,110.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrPi0, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrPi0->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelTotErrPi0->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionRelTotA, markerStyleComb+4, markerSizeComb, kBlue+2 , kBlue+2);
        graphCombPi0InvXSectionRelTotA->Draw("p,same,z");

        TLegend* legendRelTotErr = GetAndSetLegend2(0.14, 0.92-(0.035*2), 0.45, 0.92, 32);
        legendRelTotErr->AddEntry(graphCombPi0InvXSectionRelTotA,"All","p");
        legendRelTotErr->Draw();

        TLatex *labelRelTotErrEnergy    = new TLatex(0.75,0.89,collisionSystem8TeV.Data());
        SetStyleTLatex( labelRelTotErrEnergy, 0.85*textSizeLabelsPixel,4);
        labelRelTotErrEnergy->SetTextFont(43);
        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrPi0       = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelTotErrPi0, 0.85*textSizeLabelsPixel,4);
        labelRelTotErrPi0->SetTextFont(43);
        labelRelTotErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/Pi0_RelTotErr.%s",outputDir.Data(),suffix.Data()));


    histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,40.5);
    histo2DRelTotErrPi0->Draw("copy");
        graphCombPi0InvXSectionRelTotA->Draw("p,same,z");

        legendRelTotErr->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/Pi0_RelTotErrZoomed.%s",outputDir.Data(),suffix.Data()));

    histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,30.5);
    histo2DRelTotErrPi0->GetYaxis()->SetTitle("Err (%)");
    histo2DRelTotErrPi0->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionRelTotA, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombPi0InvXSectionRelTotA->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionRelStatA, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombPi0InvXSectionRelStatA->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionRelSysA, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombPi0InvXSectionRelSysA->SetLineStyle(7);
        graphCombPi0InvXSectionRelSysA->Draw("l,x0,same,e1");

        TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
        legendRelTotErr3->AddEntry(graphCombPi0InvXSectionRelTotA,"tot","p");
        legendRelTotErr3->AddEntry(graphCombPi0InvXSectionRelStatA,"stat","l");
        legendRelTotErr3->AddEntry(graphCombPi0InvXSectionRelSysA,"sys","l");
        legendRelTotErr3->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/Pi0_Reldecomp.%s",outputDir.Data(),suffix.Data()));


    // **********************************************************************************************************************
    // ************************************* Calculating bin shifted spectra & fitting **************************************
    // **********************************************************************************************************************

    // Cloning spectra
    TGraphAsymmErrors* graphCombPi0InvXSectionTotAUnshi         = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotA->Clone("Pi0Unshifted");
    TGraphAsymmErrors* graphCombPi0InvXSectionStatAUnshi        = (TGraphAsymmErrors*)graphCombPi0InvXSectionStatA->Clone("Pi0UnshiftedStat");
    TGraphAsymmErrors* graphCombPi0InvXSectionSysAUnshi         = (TGraphAsymmErrors*)graphCombPi0InvXSectionSysA->Clone("Pi0UnshiftedSys");

    TGraphAsymmErrors* graphPCMPi0InvXSectionStatUnshi          = (TGraphAsymmErrors*)graphPCMPi0InvXSectionStat->Clone("Pi0UnshiftedStatPCM");
    TGraphAsymmErrors* graphPCMPi0InvXSectionSysUnshi           = (TGraphAsymmErrors*)graphPCMPi0InvXSectionSys->Clone("Pi0UnshiftedSysPCM");

    TGraphAsymmErrors* graphPHOSPi0InvXSectionStatUnshi         = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionStat->Clone("Pi0UnshiftedStatPHOS");
    TGraphAsymmErrors* graphPHOSPi0InvXSectionSysUnshi          = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionSys->Clone("Pi0UnshiftedSysPHOS");

    TGraphAsymmErrors* graphEMCALPi0InvXSectionStatUnshi        = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionStat->Clone("Pi0UnshiftedStatEMCAL");
    TGraphAsymmErrors* graphEMCALPi0InvXSectionSysUnshi         = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionSys->Clone("Pi0UnshiftedSysEMCAL");

    TGraphAsymmErrors* graphPCMEMCALPi0InvXSectionStatUnshi     = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionStat->Clone("Pi0UnshiftedStatPCMEMCAL");
    TGraphAsymmErrors* graphPCMEMCALPi0InvXSectionSysUnshi      = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionSys->Clone("Pi0UnshiftedSysPCMEMCAL");

    TGraphAsymmErrors* graphEMCALMergedPi0InvXSectionStatUnshi  = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionStat->Clone("Pi0UnshiftedStatEMCALMerged");
    TGraphAsymmErrors* graphEMCALMergedPi0InvXSectionSysUnshi   = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionSys->Clone("Pi0UnshiftedSysEMCALMerged");


    // fitting spectrum with intial parameters
    // Two component model fit from Bylinkin
//    TF1* fitTCMDecomposedL                 = new TF1("twoCompModel_DecLow",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0),0.4,2);
//    fitTCMDecomposedL->SetParameters(graphCombPi0InvXSectionTotA->GetY()[0],0.1);
//    graphCombPi0InvXSectionStatA->Fit(fitTCMDecomposedL,"QNRMEX0+","",0.4,2.);
//    TF1 *fitTCMDecomposedH                 = new TF1("twoCompModel_DecH","[0]/(TMath::Power(1+x*x/([1]*[1]*[2]),[2]))",4,35);
//    fitTCMDecomposedH->SetParameters(graphCombPi0InvXSectionTotA->GetY()[0],0.6, 3);
//    graphCombPi0InvXSectionTotA->Fit(fitTCMDecomposedH,"QNRMEX0+","",4,35);
//    Double_t paramTCMPi0New[5]  = { fitTCMDecomposedL->GetParameter(0),fitTCMDecomposedL->GetParameter(1),
//                                    fitTCMDecomposedH->GetParameter(0),fitTCMDecomposedH->GetParameter(1),fitTCMDecomposedH->GetParameter(2)};
    Double_t paramTCMPi0New[5]  = { graphCombPi0InvXSectionTotA->GetY()[0],0.1,
                                    graphCombPi0InvXSectionTotA->GetY()[3],0.6,3.0};
    TF1* fitTCMInvXSectionPi0   = FitObject("tcm","fitTCMInvCrossSectionPi08TeV","Pi0",graphCombPi0InvXSectionStatA,0.3,100. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);

    TF1* fitTCMInvXSectionPi0highPt   = FitObject("tcm","fitTCMInvXSectionPi0highPt","Pi0",graphCombPi0InvXSectionStatA,8.,100. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);

    // Tsallis fit
    Double_t paramGraphPi0[3]                              = {5e11, 6., 0.13};
    TF1* fitInvXSectionPi0                       = FitObject("l","fitInvCrossSectionPi08TeV","Pi0",histoEMCALPi0InvXSectionStat,0.3,100.,paramGraphPi0,"QNRMEX0+");
    TF1* fitInvXSectionPi0Graph                  = (TF1*)fitInvXSectionPi0->Clone("fitInvCrossSectionPi08TeVGraph");

    // *************************************************************************************************************
    // Shift graphs in X direction if desired
    // *************************************************************************************************************
    if(bWCorrection.Contains("X")){
//         TF1* fitTsallisPi0PtMult                 = FitObject("tcmpt","TsallisMultWithPtPi08TeV","Pi0");

//         fitTsallisPi0PtMult->SetParameters( fitTCMInvXSectionPi0->GetParameter(0),fitTCMInvXSectionPi0->GetParameter(1), fitTCMInvXSectionPi0->GetParameter(2), fitTCMInvXSectionPi0->GetParameter(3),
//                                             fitTCMInvXSectionPi0->GetParameter(4)); // standard parameter optimize if necessary
        TF1* fitTsallisPi0PtMult                 = FitObject("tmpt","TsallisMultWithPtPi08TeV","Pi0");
        fitTsallisPi0PtMult->SetParameters(fitInvXSectionPi0->GetParameter(0),fitInvXSectionPi0->GetParameter(1), fitInvXSectionPi0->GetParameter(2));

        TGraphAsymmErrors* graphCombPi0InvXSectionTotANoShift = (TGraphAsymmErrors*) graphCombPi0InvXSectionTotA->Clone("Pi0_NoShift");

        graphCombPi0InvXSectionTotA              = ApplyXshift(graphCombPi0InvXSectionTotA, fitTsallisPi0PtMult);
        cout << "comb" << endl;
        graphCombPi0InvXSectionStatA->Print();
        graphCombPi0InvXSectionStatA             = ApplyXshiftIndividualSpectra (   graphCombPi0InvXSectionTotA,
                                                                                    graphCombPi0InvXSectionStatA,
                                                                                    fitTsallisPi0PtMult,
                                                                                    0, graphCombPi0InvXSectionStatA->GetN());
        graphCombPi0InvXSectionSysA              = ApplyXshiftIndividualSpectra (   graphCombPi0InvXSectionTotA,
                                                                                    graphCombPi0InvXSectionSysA,
                                                                                    fitTsallisPi0PtMult,
                                                                                    0, graphCombPi0InvXSectionSysA->GetN());
        cout << "PCM" << endl;
        graphPCMPi0InvXSectionStat               = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA,
                                                                                    graphPCMPi0InvXSectionStat,
                                                                                    fitTsallisPi0PtMult,
                                                                                    offSetPi0Shifting[0], nComBinsPi0Shifting[0]);
        graphPCMPi0InvXSectionSys                = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA,
                                                                                    graphPCMPi0InvXSectionSys,
                                                                                    fitTsallisPi0PtMult,
                                                                                    offSetPi0Shifting[0], nComBinsPi0Shifting[0]);
        cout << "PHOS" << endl;
        graphPHOSPi0InvXSectionStat              = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA,
                                                                                    graphPHOSPi0InvXSectionStat,
                                                                                    fitTsallisPi0PtMult,
                                                                                    offSetPi0Shifting[1], nComBinsPi0Shifting[1]);
        graphPHOSPi0InvXSectionSys               = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA,
                                                                                    graphPHOSPi0InvXSectionSys,
                                                                                    fitTsallisPi0PtMult,
                                                                                    offSetPi0Shifting[1], nComBinsPi0Shifting[1]);
        cout << "EMC-EMC" << endl;
        graphEMCALPi0InvXSectionStat             = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA,
                                                                                    graphEMCALPi0InvXSectionStat,
                                                                                    fitTsallisPi0PtMult,
                                                                                    offSetPi0Shifting[2], nComBinsPi0Shifting[2]);
        graphEMCALPi0InvXSectionSys              = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA,
                                                                                    graphEMCALPi0InvXSectionSys,
                                                                                    fitTsallisPi0PtMult,
                                                                                    offSetPi0Shifting[2], nComBinsPi0Shifting[2]);
        cout << "PCM-EMC" << endl;
        graphPCMEMCALPi0InvXSectionStat          = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA,
                                                                                    graphPCMEMCALPi0InvXSectionStat,
                                                                                    fitTsallisPi0PtMult,
                                                                                    offSetPi0Shifting[4], nComBinsPi0Shifting[4]);
        graphPCMEMCALPi0InvXSectionSys           = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA,
                                                                                    graphPCMEMCALPi0InvXSectionSys,
                                                                                    fitTsallisPi0PtMult,
                                                                                    offSetPi0Shifting[4], nComBinsPi0Shifting[4]);

        cout << "EMC merged" << endl;
        graphEMCALMergedPi0InvXSectionStat       = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA,
                                                                                    graphEMCALMergedPi0InvXSectionStat,
                                                                                    fitTsallisPi0PtMult,
                                                                                    offSetPi0Shifting[9], nComBinsPi0Shifting[9]);
        graphEMCALMergedPi0InvXSectionSys        = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA,
                                                                                    graphEMCALMergedPi0InvXSectionSys,
                                                                                    fitTsallisPi0PtMult,
                                                                                    offSetPi0Shifting[9], nComBinsPi0Shifting[9]);

        //***************************************************************************************************************
        //************************************Plotting binshift corrections *********************************************
        //***************************************************************************************************************

        TCanvas* canvasShift = new TCanvas("canvasShift","",0,0,1000,900);// gives the page size
        DrawGammaCanvasSettings( canvasShift, 0.10, 0.017, 0.015, 0.08);
        canvasShift->SetLogx(1);

        Size_t textSizeSpectra          = 0.04;
        TH1F * histoBinShift = new TH1F("histoBinShift","histoBinShift",1000,0.23, 110.);
        SetStyleHistoTH1ForGraphs(histoBinShift, "#it{p}_{T} (GeV/#it{c})","bin shifted (X) / no shift",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 1.1,1.2);
        histoBinShift->GetXaxis()->SetMoreLogLabels(1);
        histoBinShift->GetYaxis()->SetRangeUser(0.95,1.05);
        histoBinShift->DrawCopy();

        Int_t numberPoints   = graphCombPi0InvXSectionTotANoShift->GetN();
        Double_t *xPoint     = graphCombPi0InvXSectionTotANoShift->GetX();
        Double_t* xvalueErrUp  = graphCombPi0InvXSectionTotANoShift->GetEXhigh();
        Double_t* xvalueErrLow = graphCombPi0InvXSectionTotANoShift->GetEXlow();
        Double_t *xPointShift= graphCombPi0InvXSectionTotA->GetX();
        for (Int_t i=0; i<numberPoints; i++) {
          graphCombPi0InvXSectionTotANoShift->SetPoint(i,xPoint[i],xPointShift[i]/xPoint[i]);
          graphCombPi0InvXSectionTotANoShift->SetPointError(i,xvalueErrLow[i],xvalueErrUp[i],0,0);
        }
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionTotANoShift, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombPi0InvXSectionTotANoShift->Draw("p same");

        TLatex *labelRatioToFitBinShift   = new TLatex(0.72, 0.91, collisionSystem8TeV.Data());
        SetStyleTLatex( labelRatioToFitBinShift, textSizeLabelsPixel,4);
        labelRatioToFitBinShift->SetTextFont(43);
        labelRatioToFitBinShift->Draw();
        TLatex *labelRatioToFitALICEBinShift    = new TLatex(0.852, 0.86, "ALICE");
        SetStyleTLatex( labelRatioToFitALICEBinShift, textSizeLabelsPixel,4);
        labelRatioToFitALICEBinShift->SetTextFont(43);
        labelRatioToFitALICEBinShift->Draw();
        TLatex *labelRatioToFitPi0BinShift      = new TLatex(0.826, 0.807, "#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitPi0BinShift, textSizeLabelsPixel,4);
        labelRatioToFitPi0BinShift->SetTextFont(43);
        labelRatioToFitPi0BinShift->Draw();

        canvasShift->Update();
        canvasShift->SaveAs(Form("%s/BinShiftCorrection_Pi0.%s",outputDir.Data(),suffix.Data()));
        canvasShift->SetLogx(0);

        // *************************************************************************************************************
        // Plot control graphs
        // *************************************************************************************************************

        TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
        DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.08);
        canvasDummy2->SetLogy();
        canvasDummy2->SetLogx();
        TH2F * histo2DDummy2;
        histo2DDummy2               = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.23,110.,1000,1e-2,10e12);
        SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 0.8,1.55);
        histo2DDummy2->DrawCopy();


        DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvXSectionSys, markerStyleDet[0] ,markerSizeDet[0]/2, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        graphPCMPi0InvXSectionSys->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphPHOSPi0InvXSectionSys, markerStyleDet[1] ,markerSizeDet[1]/2, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
        graphPHOSPi0InvXSectionSys->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0InvXSectionSys, markerStyleDet[2] ,markerSizeDet[2]/2, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        graphEMCALPi0InvXSectionSys->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0InvXSectionSys, markerStyleDet[4] ,markerSizeDet[4]/2, colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
        graphPCMEMCALPi0InvXSectionSys->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphEMCALMergedPi0InvXSectionSys, markerStyleDet[9] ,markerSizeDet[9]/2, colorDet[9], colorDet[9], widthLinesBoxes, kTRUE);
        graphEMCALMergedPi0InvXSectionSys->Draw("pEsame");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatAUnshi, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
        graphCombPi0InvXSectionStatAUnshi->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatA, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombPi0InvXSectionStatA->Draw("pEsame");

        fitTCMInvXSectionPi0->SetLineColor(kBlue+2);
        fitTCMInvXSectionPi0->Draw("same");

        canvasDummy2->Update();
        canvasDummy2->Print(Form("%s/ComparisonShiftedPi0_8TeV.%s",outputDir.Data(),suffix.Data()));
        delete canvasDummy2;
    }

    // *************************************************************************************************************
    // redo fitting after binshifts
    // *************************************************************************************************************
    // Tsallis function
    graphCombPi0InvXSectionTotA->Fit(fitInvXSectionPi0,"QNRMEX0+","",0.3,100.);
    fitInvXSectionPi0           = FitObject("l","fitInvCrossSectionPi08TeV","Pi0",graphCombPi0InvXSectionTotA,0.3,100.,paramGraphPi0,"QNRMEX0+");
    fitInvXSectionPi0           = FitObject("l","fitInvCrossSectionPi08TeV","Pi0",graphCombPi0InvXSectionTotA,0.3,100.,paramGraphPi0,"QNRMEX0+");
    cout << WriteParameterToFile(fitInvXSectionPi0)<< endl;

    //Two component model from Bylinkin
    fitTCMInvXSectionPi0        = FitObject("tcm","fitTCMInvCrossSectionPi08TeV","Pi0",graphCombPi0InvXSectionTotA,0.3,100. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitTCMInvXSectionPi0)<< endl;

    TF1* fitTCMInvXSectionPi0Plot = new TF1("twoCompModel_plotting",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
    fitTCMInvXSectionPi0Plot->SetRange(0.3,100.);
    fitTCMInvXSectionPi0Plot->SetParameters(fitTCMInvXSectionPi0->GetParameters());
    fitTCMInvXSectionPi0Plot->SetParErrors(fitTCMInvXSectionPi0->GetParErrors());
//     TF1* fitTCMInvXSectionPi02  = FitObject("tcm","fitTCMInvCrossSectionPi08TeV2","Pi0",graphCombPi0InvXSectionTotA,0.4,40.,paramTCMPi0New,"QNRMEX0+","", kFALSE);
//     cout << WriteParameterToFile(fitTCMInvXSectionPi02)<< endl;

    TF1* fitPowInvXSectionPi0   = FitObject("m","fitPowInvXSectionPi08TeV","Pi0",graphCombPi0InvXSectionTotA,5,100. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvXSectionPi0)<< endl;

    TF1* fitPCMTCMInvXSectionPi0    = FitObject("tcm","fitPCMTCMInvCrossSectionPi08TeV","Pi0",graphPCMPi0InvXSectionStat,0.3,8. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
    //cout << WriteParameterToFile(fitPCMTCMInvXSectionPi0)<< endl;

    TF1* fitPHOSTCMInvXSectionPi0    = FitObject("tcm","fitPHOSTCMInvCrossSectionPi08TeV","Pi0",graphPHOSPi0InvXSectionStat,1.0,35. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);

    // *************************************************************************************************************
    // Shift graphs in Y direction as well if desired
    // *************************************************************************************************************
    TGraphAsymmErrors* graphCombPi0InvXSectionTotA_yShifted         = NULL;
    TGraphAsymmErrors* graphCombPi0InvXSectionStatA_yShifted        = NULL;
    TGraphAsymmErrors* graphCombPi0InvXSectionSysA_yShifted         = NULL;
    TGraphAsymmErrors* graphPCMPi0InvXSectionStat_yShifted          = NULL;
    TGraphAsymmErrors* graphPCMPi0InvXSectionSys_yShifted           = NULL;
    TGraphAsymmErrors* graphPHOSPi0InvXSectionStat_yShifted         = NULL;
    TGraphAsymmErrors* graphPHOSPi0InvXSectionSys_yShifted          = NULL;
    TGraphAsymmErrors* graphEMCALPi0InvXSectionStat_yShifted        = NULL;
    TGraphAsymmErrors* graphEMCALPi0InvXSectionSys_yShifted         = NULL;
    TGraphAsymmErrors* graphPCMEMCALPi0InvXSectionStat_yShifted     = NULL;
    TGraphAsymmErrors* graphPCMEMCALPi0InvXSectionSys_yShifted      = NULL;
    TGraphAsymmErrors* graphEMCALMergedPi0InvXSectionStat_yShifted  = NULL;
    TGraphAsymmErrors* graphEMCALMergedPi0InvXSectionSys_yShifted   = NULL;

    if(bWCorrection.Contains("Y") ){
        graphCombPi0InvXSectionTotA_yShifted        = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotAUnshi->Clone("Pi0YShiftedCombTot");
        graphCombPi0InvXSectionTotA_yShifted        =  ApplyYshiftIndividualSpectra( graphCombPi0InvXSectionTotA_yShifted, fitInvXSectionPi0);
        graphCombPi0InvXSectionStatA_yShifted       = (TGraphAsymmErrors*)graphCombPi0InvXSectionStatAUnshi->Clone("Pi0YShiftedCombStat");
        graphCombPi0InvXSectionStatA_yShifted       =  ApplyYshiftIndividualSpectra( graphCombPi0InvXSectionStatA_yShifted, fitInvXSectionPi0);
        graphCombPi0InvXSectionSysA_yShifted        = (TGraphAsymmErrors*)graphCombPi0InvXSectionSysAUnshi->Clone("Pi0YShiftedCombSys");
        graphCombPi0InvXSectionSysA_yShifted        =  ApplyYshiftIndividualSpectra( graphCombPi0InvXSectionSysA_yShifted, fitInvXSectionPi0);

        graphPCMPi0InvXSectionStat_yShifted         = (TGraphAsymmErrors*)graphPCMPi0InvXSectionStatUnshi->Clone("Pi0YShiftedPCMStat");
        graphPCMPi0InvXSectionStat_yShifted         =  ApplyYshiftIndividualSpectra( graphPCMPi0InvXSectionStat_yShifted, fitInvXSectionPi0);
        graphPCMPi0InvXSectionSys_yShifted          = (TGraphAsymmErrors*)graphPCMPi0InvXSectionSysUnshi->Clone("Pi0YShiftedPCMSys");
        graphPCMPi0InvXSectionSys_yShifted          =  ApplyYshiftIndividualSpectra( graphPCMPi0InvXSectionSys_yShifted, fitInvXSectionPi0);

        graphPHOSPi0InvXSectionStat_yShifted         = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionStatUnshi->Clone("Pi0YShiftedPHOSStat");
        graphPHOSPi0InvXSectionStat_yShifted         =  ApplyYshiftIndividualSpectra( graphPHOSPi0InvXSectionStat_yShifted, fitInvXSectionPi0);
        graphPHOSPi0InvXSectionSys_yShifted          = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionSysUnshi->Clone("Pi0YShiftedPHOSSys");
        graphPHOSPi0InvXSectionSys_yShifted          =  ApplyYshiftIndividualSpectra( graphPHOSPi0InvXSectionSys_yShifted, fitInvXSectionPi0);

        graphEMCALPi0InvXSectionStat_yShifted       = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionStatUnshi->Clone("Pi0YShiftedEMCStat");
        graphEMCALPi0InvXSectionStat_yShifted       =  ApplyYshiftIndividualSpectra( graphEMCALPi0InvXSectionStat_yShifted, fitInvXSectionPi0);
        graphEMCALPi0InvXSectionSys_yShifted        = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionSysUnshi->Clone("Pi0YShiftedEMCSys");
        graphEMCALPi0InvXSectionSys_yShifted        =  ApplyYshiftIndividualSpectra( graphEMCALPi0InvXSectionSys_yShifted, fitInvXSectionPi0);

        graphPCMEMCALPi0InvXSectionStat_yShifted    = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionStatUnshi->Clone("Pi0YShiftedPCMEMCStat");
        graphPCMEMCALPi0InvXSectionStat_yShifted    =  ApplyYshiftIndividualSpectra( graphPCMEMCALPi0InvXSectionStat_yShifted, fitInvXSectionPi0);
        graphPCMEMCALPi0InvXSectionSys_yShifted     = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionSysUnshi->Clone("Pi0YShiftedPCMEMCStat");
        graphPCMEMCALPi0InvXSectionSys_yShifted     =  ApplyYshiftIndividualSpectra( graphPCMEMCALPi0InvXSectionSys_yShifted, fitInvXSectionPi0);

        graphEMCALMergedPi0InvXSectionStat_yShifted = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionStatUnshi->Clone("Pi0YShiftedEMCMergedStat");
        graphEMCALMergedPi0InvXSectionStat_yShifted =  ApplyYshiftIndividualSpectra( graphEMCALMergedPi0InvXSectionStat_yShifted, fitInvXSectionPi0);
        graphEMCALMergedPi0InvXSectionSys_yShifted  = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionSysUnshi->Clone("Pi0YShiftedEMCMergedSys");
        graphEMCALMergedPi0InvXSectionSys_yShifted  =  ApplyYshiftIndividualSpectra( graphEMCALMergedPi0InvXSectionSys_yShifted, fitInvXSectionPi0);
    }

    // *************************************************************************************************************
    // Calculate ratios to combined fit
    // *************************************************************************************************************
    TH1D* histoRatioPythia8ToFit                     = (TH1D*) histoPythia8InvXSection->Clone();
    histoRatioPythia8ToFit                           = CalculateHistoRatioToFit (histoRatioPythia8ToFit, fitTCMInvXSectionPi0Plot);
    TH1D* histoRatioPythia8_4CToFit                  = (TH1D*) histoPythia8_4CInvXSection->Clone();
    histoRatioPythia8_4CToFit                        = CalculateHistoRatioToFit (histoRatioPythia8_4CToFit, fitTCMInvXSectionPi0Plot);
//    TH1D* histoRatioPythia8VarBinningToFit           = (TH1D*) histoPythia8InvXSection_VarBinning->Clone();
//    histoRatioPythia8VarBinningToFit                 = CalculateHistoRatioToFit (histoRatioPythia8VarBinningToFit, fitTCMInvXSectionPi0Plot);

    TGraph* graphRatioPi0CombNLOMuHalf               = (TGraph*)graphNLOCalcPi0MuHalf->Clone();cout << __LINE__ << endl;
    TGraph* graphRatioPi0CombNLOMuOne                = (TGraph*)graphNLOCalcPi0MuOne->Clone();cout << __LINE__ << endl;
    TGraph* graphRatioPi0CombNLOMuTwo                = (TGraph*)graphNLOCalcPi0MuTwo->Clone();cout << __LINE__ << endl;
    TGraphAsymmErrors* graphRatioPi0DSS14            = (TGraphAsymmErrors*)graphPi0DSS14->Clone();cout << __LINE__ << endl;
    //TGraphAsymmErrors* graphRatioPi0CombNLODSS14     = (TGraphAsymmErrors*)graphNLODSS14Calc->Clone();cout << __LINE__ << endl;
    graphRatioPi0CombNLOMuHalf                       = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuHalf, fitTCMInvXSectionPi0Plot); cout << __LINE__ << endl;
    graphRatioPi0CombNLOMuOne                        = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuOne, fitTCMInvXSectionPi0Plot); cout << __LINE__ << endl;
    graphRatioPi0CombNLOMuTwo                        = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuTwo, fitTCMInvXSectionPi0Plot); cout << __LINE__ << endl;
    graphRatioPi0DSS14                               = CalculateGraphErrRatioToFit (graphRatioPi0DSS14, fitTCMInvXSectionPi0Plot); cout << __LINE__ << endl;
    //graphRatioPi0CombNLODSS14                        = CalculateGraphErrRatioToFit(graphRatioPi0CombNLODSS14, fitTCMInvXSectionPi0Plot); cout << __LINE__ << endl;

    TGraphAsymmErrors* graphRatioPi0CombCombFitTotA     = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotA->Clone();
    graphRatioPi0CombCombFitTotA                        = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitTotA, fitTCMInvXSectionPi0Plot);
    TGraphAsymmErrors* graphRatioPi0CombCombFitStatA    = (TGraphAsymmErrors*)graphCombPi0InvXSectionStatA->Clone();
    graphRatioPi0CombCombFitStatA                       = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitStatA, fitTCMInvXSectionPi0Plot);
    TGraphAsymmErrors* graphRatioPi0CombCombFitSysA     = (TGraphAsymmErrors*)graphCombPi0InvXSectionSysA->Clone();
    graphRatioPi0CombCombFitSysA                        = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitSysA, fitTCMInvXSectionPi0Plot);
    TGraphAsymmErrors* graphRatioPi0PCMCombFitStat      = (TGraphAsymmErrors*)graphPCMPi0InvXSectionStat->Clone();
    graphRatioPi0PCMCombFitStat                         = CalculateGraphErrRatioToFit(graphRatioPi0PCMCombFitStat, fitTCMInvXSectionPi0Plot);
    TGraphAsymmErrors* graphRatioPi0PCMCombFitSys       = (TGraphAsymmErrors*)graphPCMPi0InvXSectionSys->Clone();
    graphRatioPi0PCMCombFitSys                          = CalculateGraphErrRatioToFit(graphRatioPi0PCMCombFitSys, fitTCMInvXSectionPi0Plot);
    TGraphAsymmErrors* graphRatioPi0PHOSCombFitStat     = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionStat->Clone();
    graphRatioPi0PHOSCombFitStat                        = CalculateGraphErrRatioToFit(graphRatioPi0PHOSCombFitStat, fitTCMInvXSectionPi0Plot);
    TGraphAsymmErrors* graphRatioPi0PHOSCombFitSys      = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionSys->Clone();
    graphRatioPi0PHOSCombFitSys                         = CalculateGraphErrRatioToFit(graphRatioPi0PHOSCombFitSys, fitTCMInvXSectionPi0Plot);
    TGraphAsymmErrors* graphRatioPi0EMCALCombFitStat    = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionStat->Clone();
    graphRatioPi0EMCALCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioPi0EMCALCombFitStat, fitTCMInvXSectionPi0Plot);
    TGraphAsymmErrors* graphRatioPi0EMCALCombFitSys     = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionSys->Clone();
    graphRatioPi0EMCALCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioPi0EMCALCombFitSys, fitTCMInvXSectionPi0Plot);
    TGraphAsymmErrors* graphRatioPi0PCMEMCALCombFitStat = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionStat->Clone();
    graphRatioPi0PCMEMCALCombFitStat                    = CalculateGraphErrRatioToFit(graphRatioPi0PCMEMCALCombFitStat, fitTCMInvXSectionPi0Plot);
    TGraphAsymmErrors* graphRatioPi0PCMEMCALCombFitSys  = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionSys->Clone();
    graphRatioPi0PCMEMCALCombFitSys                     = CalculateGraphErrRatioToFit(graphRatioPi0PCMEMCALCombFitSys, fitTCMInvXSectionPi0Plot);
    TGraphAsymmErrors* graphRatioPi0EMCALMergedCombFitStat  = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionStat->Clone();
    graphRatioPi0EMCALMergedCombFitStat                     = CalculateGraphErrRatioToFit(graphRatioPi0EMCALMergedCombFitStat, fitTCMInvXSectionPi0Plot);
    TGraphAsymmErrors* graphRatioPi0EMCALMergedCombFitSys   = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionSys->Clone();
    graphRatioPi0EMCALMergedCombFitSys                      = CalculateGraphErrRatioToFit(graphRatioPi0EMCALMergedCombFitSys, fitTCMInvXSectionPi0Plot);


    TGraphAsymmErrors* graphRatioPi0PCMCombFitStatH      = (TGraphAsymmErrors*)graphPCMPi0InvXSectionStat->Clone();
    graphRatioPi0PCMCombFitStatH                         = CalculateGraphErrRatioToFit(graphRatioPi0PCMCombFitStatH, fitTCMInvXSectionPi0highPt);
    TGraphAsymmErrors* graphRatioPi0PCMCombFitSysH       = (TGraphAsymmErrors*)graphPCMPi0InvXSectionSys->Clone();
    graphRatioPi0PCMCombFitSysH                          = CalculateGraphErrRatioToFit(graphRatioPi0PCMCombFitSysH, fitTCMInvXSectionPi0highPt);
    TGraphAsymmErrors* graphRatioPi0PHOSCombFitStatH     = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionStat->Clone();
    graphRatioPi0PHOSCombFitStatH                        = CalculateGraphErrRatioToFit(graphRatioPi0PHOSCombFitStatH, fitTCMInvXSectionPi0highPt);
    TGraphAsymmErrors* graphRatioPi0PHOSCombFitSysH      = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionSys->Clone();
    graphRatioPi0PHOSCombFitSysH                         = CalculateGraphErrRatioToFit(graphRatioPi0PHOSCombFitSysH, fitTCMInvXSectionPi0highPt);
    TGraphAsymmErrors* graphRatioPi0EMCALCombFitStatH    = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionStat->Clone();
    graphRatioPi0EMCALCombFitStatH                       = CalculateGraphErrRatioToFit(graphRatioPi0EMCALCombFitStatH, fitTCMInvXSectionPi0highPt);
    TGraphAsymmErrors* graphRatioPi0EMCALCombFitSysH     = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionSys->Clone();
    graphRatioPi0EMCALCombFitSysH                        = CalculateGraphErrRatioToFit(graphRatioPi0EMCALCombFitSysH, fitTCMInvXSectionPi0highPt);
    TGraphAsymmErrors* graphRatioPi0PCMEMCALCombFitStatH = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionStat->Clone();
    graphRatioPi0PCMEMCALCombFitStatH                    = CalculateGraphErrRatioToFit(graphRatioPi0PCMEMCALCombFitStatH, fitTCMInvXSectionPi0highPt);
    TGraphAsymmErrors* graphRatioPi0PCMEMCALCombFitSysH  = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionSys->Clone();
    graphRatioPi0PCMEMCALCombFitSysH                     = CalculateGraphErrRatioToFit(graphRatioPi0PCMEMCALCombFitSysH, fitTCMInvXSectionPi0highPt);
    TGraphAsymmErrors* graphRatioPi0EMCALMergedCombFitStatH  = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionStat->Clone();
    graphRatioPi0EMCALMergedCombFitStatH                     = CalculateGraphErrRatioToFit(graphRatioPi0EMCALMergedCombFitStatH, fitTCMInvXSectionPi0highPt);
    TGraphAsymmErrors* graphRatioPi0EMCALMergedCombFitSysH   = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionSys->Clone();
    graphRatioPi0EMCALMergedCombFitSysH                      = CalculateGraphErrRatioToFit(graphRatioPi0EMCALMergedCombFitSysH, fitTCMInvXSectionPi0highPt);
    // **********************************************************************************************************************
    // ******************************************* Plot Ratio of Comb to Fit ****************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 54;
    TCanvas* canvasRatioToCombFit       = new TCanvas("canvasRatioToCombFit","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToCombFit, 0.08, 0.01, 0.01, 0.125);
    canvasRatioToCombFit->SetLogx();

        Double_t textsizeLabelsPP       = 0;
        Double_t textsizeFacPP          = 0;
        if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
            textsizeLabelsPP            = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
            textsizeFacPP               = (Double_t)1./canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
        } else {
            textsizeLabelsPP            = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
            textsizeFacPP               = (Double_t)1./canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
        }
        cout << textsizeLabelsPP << endl;

    TH2F * histo2DPi0RatioToCombFit;
    histo2DPi0RatioToCombFit               = new TH2F("histo2DPi0RatioToCombFit","histo2DPi0RatioToCombFit",1000,0.23,110.,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DPi0RatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
    histo2DPi0RatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DPi0RatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
    histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
    histo2DPi0RatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioPi0CombCombFitSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioPi0CombCombFitStatA->Draw("p,same,z");

        DrawGammaLines(0.23, 110. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.23, 110. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.23, 110. , 0.9, 0.9,0.1, kGray, 7);

        TLatex *labelRatioToFitEnergy   = new TLatex(0.72, 0.91, collisionSystem8TeV.Data());
        SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel,4);
        labelRatioToFitEnergy->SetTextFont(43);
        labelRatioToFitEnergy->Draw();
        TLatex *labelRatioToFitALICE    = new TLatex(0.852, 0.86, "ALICE");
        SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel,4);
        labelRatioToFitALICE->SetTextFont(43);
        labelRatioToFitALICE->Draw();
        TLatex *labelRatioToFitPi0      = new TLatex(0.826, 0.807, "#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitPi0, textSizeLabelsPixel,4);
        labelRatioToFitPi0->SetTextFont(43);
        labelRatioToFitPi0->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToCombFit_PP8TeV.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // *******************************************Plot Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DPi0RatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitSys, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitStat, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PHOSCombFitSys, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PHOSCombFitStat, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALCombFitSys, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALCombFitStat, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALMergedCombFitSys, markerStyleDet[9] ,markerSizeDet[9]*0.5, colorDet[9], colorDet[9], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALMergedCombFitStat, markerStyleDet[9] ,markerSizeDet[9]*0.5, colorDet[9], colorDet[9]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMEMCALCombFitSys, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMEMCALCombFitStat, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4]);

        graphRatioPi0PCMCombFitSys->Draw("E2same");
        graphRatioPi0PHOSCombFitSys->Draw("E2same");
        graphRatioPi0EMCALCombFitSys->Draw("E2same");
        graphRatioPi0EMCALMergedCombFitSys->Draw("E2same");
        graphRatioPi0PCMEMCALCombFitSys->Draw("E2same");

        graphRatioPi0PCMCombFitStat->Draw("p,same,z");
        graphRatioPi0PHOSCombFitStat->Draw("p,same,z");
        graphRatioPi0EMCALCombFitStat->Draw("p,same,z");
        graphRatioPi0EMCALMergedCombFitStat->Draw("p,same,z");
        graphRatioPi0PCMEMCALCombFitStat->Draw("p,same,z");

        DrawGammaLines(0.23, 110. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.23, 110. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.23, 110. , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitPi0->Draw();

        //****************************** Definition of the Legend ******************************************
        //**************** Row def ************************
        Double_t rowsLegendOnlyPi0Ratio[6]          = {0.91, 0.855, 0.805, 0.755, 0.705, 0.655};
        Double_t rowsLegendOnlyPi0RatioAbs[6]       = {0.91, 2.165, 2.035, 1.905, 1.765, 1.635};
        Double_t columnsLegendOnlyPi0Ratio[3]       = {0.175, 0.445, 0.53};
        Double_t columnsLegendOnlyPi0RatioAbs[3]    = {0.115, 2.9, 4.5};
        Double_t lengthBox                          = 0.2;
        Double_t heightBox                          = 0.08/2;
        //****************** first Column **************************************************
        TLatex *textPCMOnlyRatioPi0                 = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[1],"PCM");
        SetStyleTLatex( textPCMOnlyRatioPi0, textSizeLabelsPixel,4);
        textPCMOnlyRatioPi0->SetTextFont(43);
        textPCMOnlyRatioPi0->Draw();
        TLatex *textPHOSOnlyRatioPi0                = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[2],"PHOS");
        SetStyleTLatex( textPHOSOnlyRatioPi0, textSizeLabelsPixel,4);
        textPHOSOnlyRatioPi0->SetTextFont(43);
        textPHOSOnlyRatioPi0->Draw();
        TLatex *textEMCALOnlyRatioPi0               = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[3],"EMCal");
        SetStyleTLatex( textEMCALOnlyRatioPi0, textSizeLabelsPixel,4);
        textEMCALOnlyRatioPi0->SetTextFont(43);
        textEMCALOnlyRatioPi0->Draw();
        TLatex *textPCMEMCALOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[4],"PCM-EMCal");
        SetStyleTLatex( textPCMEMCALOnlyRatioPi0, textSizeLabelsPixel,4);
        textPCMEMCALOnlyRatioPi0->SetTextFont(43);
        textPCMEMCALOnlyRatioPi0->Draw();
        TLatex *textEMCALMergedOnlyRatioPi0         = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[5],"EMCal merged");
        SetStyleTLatex( textEMCALMergedOnlyRatioPi0, textSizeLabelsPixel,4);
        textEMCALMergedOnlyRatioPi0->SetTextFont(43);
        textEMCALMergedOnlyRatioPi0->Draw();

        //****************** second Column *************************************************
        TLatex *textStatOnlyRatioPi0                = new TLatex(columnsLegendOnlyPi0Ratio[1],rowsLegendOnlyPi0Ratio[0] ,"stat");
        SetStyleTLatex( textStatOnlyRatioPi0, textSizeLabelsPixel,4);
        textStatOnlyRatioPi0->SetTextFont(43);
        textStatOnlyRatioPi0->Draw();
        TLatex *textSysOnlyRatioPi0                 = new TLatex(columnsLegendOnlyPi0Ratio[2] ,rowsLegendOnlyPi0Ratio[0],"syst");
        SetStyleTLatex( textSysOnlyRatioPi0, textSizeLabelsPixel,4);
        textSysOnlyRatioPi0->SetTextFont(43);
        textSysOnlyRatioPi0->Draw();
        TMarker* markerPCMPi0OnlyRatioPi0           = CreateMarkerFromGraph(graphRatioPi0PCMCombFitSys,columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[1],1);
        markerPCMPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[1]);
        TMarker* markerPHOSPi0OnlyRatioPi0          = CreateMarkerFromGraph(graphRatioPi0PHOSCombFitSys, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[2],1);
        markerPHOSPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[2]);
        TMarker* markerEMCALPi0OnlyRatioPi0         = CreateMarkerFromGraph(graphRatioPi0EMCALCombFitSys, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[3],1);
        markerEMCALPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[3]);
        TMarker* markerPCMEMCALPi0OnlyRatioPi0      = CreateMarkerFromGraph(graphRatioPi0PCMEMCALCombFitSys, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[4],1);
        markerPCMEMCALPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[4]);
        TMarker* markerEMCALMergedPi0OnlyRatioPi0   = CreateMarkerFromGraph(graphRatioPi0EMCALMergedCombFitSys, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[5],1);
        markerEMCALMergedPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[5]);

        TBox* boxPCMPi0OnlyRatioPi0                 = CreateBoxFromGraph(graphRatioPi0PCMCombFitSys, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[1]- heightBox,
                                                                         columnsLegendOnlyPi0RatioAbs[2]+ 5*lengthBox, rowsLegendOnlyPi0RatioAbs[1]+ heightBox);
        boxPCMPi0OnlyRatioPi0->Draw("l");
        TBox* boxPHOSPi0OnlyRatioPi0                = CreateBoxFromGraph(graphRatioPi0PHOSCombFitSys, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[2]- heightBox,
                                                                         columnsLegendOnlyPi0RatioAbs[2]+ 5*lengthBox, rowsLegendOnlyPi0RatioAbs[2]+ heightBox);
        boxPHOSPi0OnlyRatioPi0->Draw("l");
        TBox* boxEMCALPi0OnlyRatioPi0               = CreateBoxFromGraph(graphRatioPi0EMCALCombFitSys, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[3]- heightBox,
                                                                         columnsLegendOnlyPi0RatioAbs[2]+ 5*lengthBox, rowsLegendOnlyPi0RatioAbs[3]+ heightBox);
        boxEMCALPi0OnlyRatioPi0->Draw("l");
        TBox* boxPCMEMCALPi0OnlyRatioPi0            = CreateBoxFromGraph(graphRatioPi0PCMEMCALCombFitSys, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[4]- heightBox,
                                                                         columnsLegendOnlyPi0RatioAbs[2]+ 5*lengthBox, rowsLegendOnlyPi0RatioAbs[4]+ heightBox);
        boxPCMEMCALPi0OnlyRatioPi0->Draw("l");
        TBox* boxEMCALMergedPi0OnlyRatioPi0         = CreateBoxFromGraph(graphRatioPi0EMCALMergedCombFitSys, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[5]- heightBox,
                                                                         columnsLegendOnlyPi0RatioAbs[2]+ 5*lengthBox, rowsLegendOnlyPi0RatioAbs[5]+ heightBox);
        boxEMCALMergedPi0OnlyRatioPi0->Draw("l");

    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFit_PP.%s",outputDir.Data(),suffix.Data()));


    // **********************************************************************************************************************
    // *******************************************Plot Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DPi0RatioToCombFit->GetXaxis()->SetRangeUser(7.,79.);
    histo2DPi0RatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitSysH, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitStatH, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PHOSCombFitSysH, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PHOSCombFitStatH, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALCombFitSysH, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALCombFitStatH, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALMergedCombFitSysH, markerStyleDet[9] ,markerSizeDet[9]*0.5, colorDet[9], colorDet[9], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALMergedCombFitStatH, markerStyleDet[9] ,markerSizeDet[9]*0.5, colorDet[9], colorDet[9]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMEMCALCombFitSysH, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMEMCALCombFitStatH, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4]);

        graphRatioPi0PCMCombFitSysH->Draw("E2same");
        graphRatioPi0PHOSCombFitSysH->Draw("E2same");
        graphRatioPi0EMCALCombFitSysH->Draw("E2same");
        graphRatioPi0EMCALMergedCombFitSysH->Draw("E2same");
        graphRatioPi0PCMEMCALCombFitSysH->Draw("E2same");

        graphRatioPi0PCMCombFitStatH->Draw("p,same,z");
        graphRatioPi0PHOSCombFitStatH->Draw("p,same,z");
        graphRatioPi0EMCALCombFitStatH->Draw("p,same,z");
        graphRatioPi0EMCALMergedCombFitStatH->Draw("p,same,z");
        graphRatioPi0PCMEMCALCombFitStatH->Draw("p,same,z");

        DrawGammaLines(7, 110. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(7, 110. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(7, 110. , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitPi0->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFit_PP_highPt.%s",outputDir.Data(),suffix.Data()));

    // *******************************************************************************************************
    // ********************** Ratio to standalone PCM fit ****************************************************
    // *******************************************************************************************************
    TGraphAsymmErrors* graphRatioPi0PCMCombFitPCMStat      = (TGraphAsymmErrors*)graphPCMPi0InvXSectionStat->Clone();
    graphRatioPi0PCMCombFitPCMStat                         = CalculateGraphErrRatioToFit(graphRatioPi0PCMCombFitPCMStat, fitPCMTCMInvXSectionPi0);
    TGraphAsymmErrors* graphRatioPi0PCMCombFitPCMSys       = (TGraphAsymmErrors*)graphPCMPi0InvXSectionSys->Clone();
    graphRatioPi0PCMCombFitPCMSys                          = CalculateGraphErrRatioToFit(graphRatioPi0PCMCombFitPCMSys, fitPCMTCMInvXSectionPi0);

    canvasRatioToCombFit->cd();
    histo2DPi0RatioToCombFit->GetXaxis()->SetRangeUser(0.23,16);
    histo2DPi0RatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitSys, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitStat, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitPCMSys, markerStyleDet[0]+4 ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitPCMStat, markerStyleDet[0]+4 ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);

        graphRatioPi0PCMCombFitSys->Draw("E2same");
        graphRatioPi0PCMCombFitPCMSys->Draw("E2same");

        graphRatioPi0PCMCombFitStat->Draw("p,same,z");
        graphRatioPi0PCMCombFitPCMStat->Draw("p,same,z");

        DrawGammaLines(0.23, 110. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.23, 110. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.23, 110. , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitPi0->Draw();

        TLegend* legendPi0FitPCMStandalone        = GetAndSetLegend2(0.2, 0.92-(0.05*2), 0.45, 0.92, 38);
        legendPi0FitPCMStandalone->AddEntry(graphRatioPi0PCMCombFitSys,"PCM/comb fit");
        legendPi0FitPCMStandalone->AddEntry(graphRatioPi0PCMCombFitPCMSys,"PCM/PCM standalone fit");
        legendPi0FitPCMStandalone->Draw("");

    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfPCMToCombAndPCMStandaloneFit_PP.%s",outputDir.Data(),suffix.Data()));


    // *******************************************************************************************************
    // ********************** Ratio of PHOS prelim/final **************************************************
    // *******************************************************************************************************
    TGraphErrors* graphRatioPHOSprelimFinalMBstat     = (TGraphErrors*)filePHOSprelim->Get("gInvCrossSectionNoBinShift_MB_PHOS");
    graphRatioPHOSprelimFinalMBstat                   = CalculateGraphErrRatioToFit(graphRatioPHOSprelimFinalMBstat, fitPHOSTCMInvXSectionPi0);
    TGraphErrors* graphRatioPHOSprelimFinalMBsys      = (TGraphErrors*)filePHOSprelim->Get("gSysErrInvCrossSectionNoBinShift_MB_PHOS");
    graphRatioPHOSprelimFinalMBsys                    = CalculateGraphErrRatioToFit(graphRatioPHOSprelimFinalMBsys, fitPHOSTCMInvXSectionPi0);
    TGraphErrors* graphRatioPHOSprelimFinalPHOSstat   = (TGraphErrors*)filePHOSprelim->Get("gInvCrossSectionNoBinShift_PHOS_PHOS");
    graphRatioPHOSprelimFinalPHOSstat                 = CalculateGraphErrRatioToFit(graphRatioPHOSprelimFinalPHOSstat, fitPHOSTCMInvXSectionPi0);
    TGraphErrors* graphRatioPHOSprelimFinalPHOSsys    = (TGraphErrors*)filePHOSprelim->Get("gSysErrInvCrossSectionNoBinShift_PHOS_PHOS");
    graphRatioPHOSprelimFinalPHOSsys                  = CalculateGraphErrRatioToFit(graphRatioPHOSprelimFinalPHOSsys, fitPHOSTCMInvXSectionPi0);

    canvasRatioToCombFit->cd();
    histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.5,1.95);
    histo2DPi0RatioToCombFit->GetXaxis()->SetRangeUser(0.8,45);
    histo2DPi0RatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphErr(graphRatioPHOSprelimFinalMBsys, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphErr(graphRatioPHOSprelimFinalMBstat, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);

        DrawGammaSetMarkerTGraphErr(graphRatioPHOSprelimFinalPHOSsys, markerStyleDet[1]+4 ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphErr(graphRatioPHOSprelimFinalPHOSstat, markerStyleDet[1]+4 ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);

        graphRatioPHOSprelimFinalMBsys->Draw("E2same");
        graphRatioPHOSprelimFinalPHOSsys->Draw("E2same");

        graphRatioPHOSprelimFinalMBstat->Draw("p,same,z");
        graphRatioPHOSprelimFinalPHOSstat->Draw("p,same,z");

        DrawGammaLines(0.8, 45. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.8, 45. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.8, 45. , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitPi0->Draw();

        TLegend* legendPi0FitPHOSStandaloneA        = GetAndSetLegend2(0.2, 0.92-(0.05*2), 0.45, 0.92, 38);
        legendPi0FitPHOSStandaloneA->AddEntry(graphRatioPHOSprelimFinalMBsys,"PHOS MB prelim / final");
        legendPi0FitPHOSStandaloneA->AddEntry(graphRatioPHOSprelimFinalPHOSsys,"PHOS PHOS prelim / final");
        legendPi0FitPHOSStandaloneA->Draw("");

    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfPHOSprelimToPHOSfinal_PP.%s",outputDir.Data(),suffix.Data()));

    // *******************************************************************************************************
    // ********************** Ratio of PHOS prelim/comb final **************************************************
    // *******************************************************************************************************
    TGraphErrors* graphRatioPHOSprelimMBstat      = (TGraphErrors*)filePHOSprelim->Get("gInvCrossSectionNoBinShift_MB_PHOS");
    graphRatioPHOSprelimMBstat                         = CalculateGraphErrRatioToFit(graphRatioPHOSprelimMBstat, fitTCMInvXSectionPi0Plot);
    TGraphErrors* graphRatioPHOSprelimMBsys      = (TGraphErrors*)filePHOSprelim->Get("gSysErrInvCrossSectionNoBinShift_MB_PHOS");
    graphRatioPHOSprelimMBsys                         = CalculateGraphErrRatioToFit(graphRatioPHOSprelimMBsys, fitTCMInvXSectionPi0Plot);
    TGraphErrors* graphRatioPHOSprelimPHOSstat    = (TGraphErrors*)filePHOSprelim->Get("gInvCrossSectionNoBinShift_PHOS_PHOS");
    graphRatioPHOSprelimPHOSstat                       = CalculateGraphErrRatioToFit(graphRatioPHOSprelimPHOSstat, fitTCMInvXSectionPi0Plot);
    TGraphErrors* graphRatioPHOSprelimPHOSsys      = (TGraphErrors*)filePHOSprelim->Get("gSysErrInvCrossSectionNoBinShift_PHOS_PHOS");
    graphRatioPHOSprelimPHOSsys                         = CalculateGraphErrRatioToFit(graphRatioPHOSprelimPHOSsys, fitTCMInvXSectionPi0Plot);

    canvasRatioToCombFit->cd();
    histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.5,1.95);
    histo2DPi0RatioToCombFit->GetXaxis()->SetRangeUser(0.8,45);
    histo2DPi0RatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphErr(graphRatioPHOSprelimMBsys, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphErr(graphRatioPHOSprelimMBstat, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);

        DrawGammaSetMarkerTGraphErr(graphRatioPHOSprelimPHOSsys, markerStyleDet[1]+4 ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphErr(graphRatioPHOSprelimPHOSstat, markerStyleDet[1]+4 ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);

        graphRatioPHOSprelimMBsys->Draw("E2same");
        graphRatioPHOSprelimPHOSsys->Draw("E2same");

        graphRatioPHOSprelimMBstat->Draw("p,same,z");
        graphRatioPHOSprelimPHOSstat->Draw("p,same,z");

        DrawGammaLines(0.8, 45. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.8, 45. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.8, 45. , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitPi0->Draw();

        TLegend* legendPi0FitPHOSStandalone        = GetAndSetLegend2(0.2, 0.92-(0.05*2), 0.45, 0.92, 38);
        legendPi0FitPHOSStandalone->AddEntry(graphRatioPHOSprelimMBsys,"PHOS MB prelim / combined final");
        legendPi0FitPHOSStandalone->AddEntry(graphRatioPHOSprelimPHOSsys,"PHOS PHOS prelim / combined final");
        legendPi0FitPHOSStandalone->Draw("");

    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfPHOSprelimToCombFit_PP.%s",outputDir.Data(),suffix.Data()));


    //********************************************************************************************************
    // Plotting simple comparison of data vs fit to pi0 meson spec
    //********************************************************************************************************
    TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.08);
    canvasDummy2->SetLogy();
    canvasDummy2->SetLogx();
    TH2F* histo2DDummy3;
    histo2DDummy3               = new TH2F("histo2DDummy3","histo2DDummy3",1000,0.33,110.,1000,1e-2,1e11);
    SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 0.8,1.55);

    TF1* fitTCMDecomposedPi0L                 = new TF1("twoCompModel_DecLow",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
    fitTCMDecomposedPi0L->SetParameters(fitTCMInvXSectionPi0Plot->GetParameter(0),fitTCMInvXSectionPi0Plot->GetParameter(1));
    fitTCMDecomposedPi0L->SetRange(0.3,110);
    TF1 *fitTCMDecomposedPi0H                 = new TF1("twoCompModel_DecH","[0]/(TMath::Power(1+x*x/([1]*[1]*[2]),[2]))");
   //      graphCombEtaInvXSectionTotA->Fit(fitTCMDecomposedH,"QNRMEX0+","",5,20);
    fitTCMDecomposedPi0H->SetParameters(fitTCMInvXSectionPi0Plot->GetParameter(2),fitTCMInvXSectionPi0Plot->GetParameter(3), fitTCMInvXSectionPi0Plot->GetParameter(4));
    fitTCMDecomposedPi0H->SetRange(0.3,110);

    histo2DDummy3               = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.23,110.,1000,1e-2,9e12);
    SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 0.8,1.55);
    histo2DDummy3->DrawCopy();

    DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatAUnshi, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
    graphCombPi0InvXSectionStatAUnshi->Draw("pEsame");
    DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatA, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
    graphCombPi0InvXSectionStatA->Draw("pEsame");

//    fitInvXSectionPi0->SetLineColor(kBlue+2);
//    fitInvXSectionPi0->Draw("same");
    fitTCMInvXSectionPi0Plot->SetLineColor(kRed+2);
    fitTCMInvXSectionPi0Plot->SetRange(0.3,110.);
    fitTCMInvXSectionPi0Plot->Draw("same");

    fitTCMDecomposedPi0L->SetLineColor(kAzure);
    fitTCMDecomposedPi0L->SetLineStyle(2);
    fitTCMDecomposedPi0L->Draw("same");
    fitTCMDecomposedPi0H->SetLineColor(kGreen+2);
    fitTCMDecomposedPi0H->SetLineStyle(8);
    fitTCMDecomposedPi0H->Draw("same");

    TLatex *labelTCMPi01= new TLatex(0.48, 0.94, Form("TCM low:"));
    TLatex *labelTCMPi02= new TLatex(0.48, 0.90, Form("A_{1}: (%.1e #pm %.1e) - T_{e}: (%.3f #pm %.3f)",fitTCMInvXSectionPi0Plot->GetParameter(0),fitTCMInvXSectionPi0Plot->GetParError(0),fitTCMInvXSectionPi0Plot->GetParameter(1),fitTCMInvXSectionPi0Plot->GetParError(1)));
    TLatex *labelTCMPi03= new TLatex(0.48, 0.86, Form("TCM high:"));
    TLatex *labelTCMPi04= new TLatex(0.48, 0.82, Form("A_{2}: (%.1e #pm %.1e) - T: (%.3f #pm %.3f) - n: (%.3f #pm %.3f)",fitTCMInvXSectionPi0Plot->GetParameter(2),fitTCMInvXSectionPi0Plot->GetParError(2),abs(fitTCMInvXSectionPi0Plot->GetParameter(3)),fitTCMInvXSectionPi0Plot->GetParError(3),fitTCMInvXSectionPi0Plot->GetParameter(4),fitTCMInvXSectionPi0Plot->GetParError(4)));

    TLatex *labelTCMPi05= new TLatex(0.55, 0.75, Form("Bylinkin-Rostovtsev:"));
    TLatex *labelTCMPi06= new TLatex(0.55, 0.71, Form("#it{A}_{1} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}_{2}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{n}"));

    SetStyleTLatex( labelTCMPi01, 0.03,4);
    labelTCMPi01->Draw();
    SetStyleTLatex( labelTCMPi02, 0.02,4);
    labelTCMPi02->Draw();
    SetStyleTLatex( labelTCMPi03, 0.03,4);
    labelTCMPi03->Draw();
    SetStyleTLatex( labelTCMPi04, 0.02,4);
    labelTCMPi04->Draw();
    SetStyleTLatex( labelTCMPi05, 0.03,4);
    labelTCMPi05->Draw();
    SetStyleTLatex( labelTCMPi06, 0.03,4);
    labelTCMPi06->Draw();

    TLatex *labelRelSysErrEnergyC    = new TLatex(0.18,0.94,collisionSystem8TeV.Data());
    SetStyleTLatex( labelRelSysErrEnergyC, 0.85*textSizeLabelsPixel,4);
    labelRelSysErrEnergyC->SetTextFont(43);
    labelRelSysErrEnergyC->Draw();
    TLatex *labelRelSysErrPi0C       = new TLatex(0.18,0.9,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrPi0C, 0.85*textSizeLabelsPixel,4);
    labelRelSysErrPi0C->SetTextFont(43);
    labelRelSysErrPi0C->Draw();

    TLegend* legendWithFit   = GetAndSetLegend2(0.17, 0.14, 0.5, 0.14+(0.035*3), 32);
    legendWithFit->AddEntry(fitTCMDecomposedPi0L,"TCM low","l");
    legendWithFit->AddEntry(fitTCMDecomposedPi0H,"TCM high","l");
    legendWithFit->AddEntry(fitTCMInvXSectionPi0Plot,"Bylinkin-Rostovtsev (TCM)","l");
    legendWithFit->Draw();

    canvasDummy2->Update();
    canvasDummy2->Print(Form("%s/ComparisonWithFitPi0_8TeV.%s",outputDir.Data(),suffix.Data()));
    //********************************************************************************************************
    canvasDummy2->Clear();
    histo2DDummy3->DrawCopy();

    graphCombPi0InvXSectionStatAUnshi->Draw("pEsame");
    graphCombPi0InvXSectionStatA->Draw("pEsame");

    fitInvXSectionPi0->SetLineColor(kRed+2);
    fitInvXSectionPi0->Draw("same");

    TLatex *labelTCMEta20 = new TLatex(0.35, 0.90, Form("dN/dy: (%.1e #pm %.1e) - n: (%.3f #pm %.3f) - T_{Levy} (GeV/c): (%.3f #pm %.3f)",fitInvXSectionPi0->GetParameter(0),fitInvXSectionPi0->GetParError(0),fitInvXSectionPi0->GetParameter(1),fitInvXSectionPi0->GetParError(1),fitInvXSectionPi0->GetParameter(2),fitInvXSectionPi0->GetParError(2)));
    SetStyleTLatex( labelTCMEta20, 0.02,4);
    labelTCMEta20->Draw();

    labelRelSysErrEnergyC->Draw();
    labelRelSysErrPi0C->Draw();

    TLegend* legendWithFitPi02   = GetAndSetLegend2(0.17, 0.14, 0.5, 0.14+(0.035*3), 32);
    legendWithFitPi02->AddEntry(fitInvXSectionPi0,"Levy-Tsallis","l");
    legendWithFitPi02->Draw();

    canvasDummy2->Update();
    canvasDummy2->Print(Form("%s/ComparisonWithFit_Tsallis_Pi0_8TeV.%s",outputDir.Data(),suffix.Data()));

    delete canvasDummy2;
    delete histo2DDummy3;

    // **********************************************************************************************************************
    // ******************************************* Mass and width for pi0 at 8TeV ****************************************
    // **********************************************************************************************************************

    Double_t arrayBoundariesX1_4[2];
    Double_t arrayBoundariesY1_4[3];
    Double_t relativeMarginsX[3];
    Double_t relativeMarginsY[3];
    textSizeLabelsPixel             = 50;
    ReturnCorrectValuesForCanvasScaling(1350,1250, 1, 2,0.09, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

    TCanvas* canvasMassWidthPi0     = new TCanvas("canvasMassWidthPi0","",0,0,1350,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasMassWidthPi0,  0.13, 0.02, 0.03, 0.06);

    TPad* padWidthPi0               = new TPad("padWidthPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthPi0->Draw();

    TPad* padMassPi0                = new TPad("padMassPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
    padMassPi0->Draw();

    TPad* padMassLegend1            = new TPad("padMassLegend1", "", 0.13, 0.36, 0.52, 0.52,-1, -1, -2);
    DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
    padMassLegend1->SetFillStyle(0);
    padMassLegend1->Draw();

    padWidthPi0->cd();
    padWidthPi0->SetLogx();

        Double_t margin                 = relativeMarginsX[0]*2.7*1350;
        Double_t textsizeLabelsWidth    = 0;
        Double_t textsizeFacWidth       = 0;
        if (padWidthPi0->XtoPixel(padWidthPi0->GetX2()) < padWidthPi0->YtoPixel(padWidthPi0->GetY1())){
            textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi0->XtoPixel(padWidthPi0->GetX2()) ;
            textsizeFacWidth            = (Double_t)1./padWidthPi0->XtoPixel(padWidthPi0->GetX2()) ;
        } else {
            textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi0->YtoPixel(padWidthPi0->GetY1());
            textsizeFacWidth            = (Double_t)1./padWidthPi0->YtoPixel(padWidthPi0->GetY1());
        }
        cout << textsizeLabelsWidth << endl;

        TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, 0.23, 50. ,1000., -30, 40);
        SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                                  0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
        histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,32.5);
        histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
        histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
        histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
        histo2DAllPi0FWHM->DrawCopy();

        DrawGammaSetMarker(histoPCMPi0FWHMMeV, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        histoPCMPi0FWHMMeV->Draw("p,same,e");
        DrawGammaSetMarker(histoPCMPi0TrueFWHMMeV, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
        histoPCMPi0TrueFWHMMeV->Draw("p,same,e");

        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0FWHM, markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
        graphEMCALPi0FWHM->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0FWHMMC, markerStyleDetMC[2], markerSizeDetMC[2]*0.55, colorDetMC[2] , colorDetMC[2]);
        graphEMCALPi0FWHMMC->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0FWHM, markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
        graphPCMEMCALPi0FWHM->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0FWHMMC, markerStyleDetMC[4], markerSizeDetMC[4]*0.55, colorDetMC[4] , colorDetMC[4]);
        graphPCMEMCALPi0FWHMMC->Draw("p,same,z");

        TLatex *labelLegendAMass    = new TLatex(0.13,0.06,"a)");
        SetStyleTLatex( labelLegendAMass, textSizeLabelsPixel,4);
        labelLegendAMass->SetTextFont(43);
        labelLegendAMass->Draw();

        TLatex *labelMassPerf       = new TLatex(0.13,0.87,"ALICE performance");
        SetStyleTLatex( labelMassPerf, textSizeLabelsPixel,4);
        labelMassPerf->SetTextFont(43);
        labelMassPerf->Draw();
        TLatex *labelMassEnergy     = new TLatex(0.13,0.78,collisionSystem8TeV.Data());
        SetStyleTLatex( labelMassEnergy, textSizeLabelsPixel,4);
        labelMassEnergy->SetTextFont(43);
        labelMassEnergy->Draw();
        TLatex *labelMassPi0        = new TLatex(0.13,0.69,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelMassPi0, textSizeLabelsPixel,4);
        labelMassPi0->SetTextFont(43);
        labelMassPi0->Draw();

    padMassPi0->cd();
    padMassPi0->SetLogx();

        Double_t textsizeLabelsMass         = 0;
        Double_t textsizeFacMass            = 0;
        if (padMassPi0->XtoPixel(padMassPi0->GetX2()) <padMassPi0->YtoPixel(padMassPi0->GetY1()) ){
            textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
            textsizeFacMass                 = (Double_t)1./padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
        } else {
            textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi0->YtoPixel(padMassPi0->GetY1());
            textsizeFacMass                 = (Double_t)1./padMassPi0->YtoPixel(padMassPi0->GetY1());
        }

        TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, 0.23, 50., 1000., 120., 175);
        SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                                  textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
        histo2DAllPi0Mass->GetYaxis()->SetRangeUser(122.5,171.8);
        histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
        histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
        histo2DAllPi0Mass->GetXaxis()->SetLabelOffset(-0.015);
        histo2DAllPi0Mass->DrawCopy();

        DrawGammaSetMarker(histoPCMPi0Mass, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        histoPCMPi0Mass->Draw("p,same,e");
        DrawGammaSetMarker(histoPCMPi0TrueMass, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
        histoPCMPi0TrueMass->Draw("p,same,e");

        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0Mass, markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
        graphEMCALPi0Mass->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0MassMC, markerStyleDetMC[2], markerSizeDetMC[2]*0.55, colorDetMC[2] , colorDetMC[2]);
        graphEMCALPi0MassMC->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0Mass, markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
        graphPCMEMCALPi0Mass->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0MassMC, markerStyleDetMC[4], markerSizeDetMC[4]*0.55, colorDetMC[4] , colorDetMC[4]);
        graphPCMEMCALPi0MassMC->Draw("p,same,z");


        DrawGammaLines(0.23, 50. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);

        TLatex *labelLegendBMass            = new TLatex(0.13,0.22,"b)");
        SetStyleTLatex( labelLegendBMass, textSizeLabelsPixel,4);
        labelLegendBMass->SetTextFont(43);
        labelLegendBMass->Draw();

        //********************************** Defintion of the Legend **************************************************
        Double_t columnsLegendMass2[3]      = {0.,0.57,0.84};
        Double_t rowsLegendMass2[4]         = {0.75,0.5,0.25,0.01};
        //******************* Offsets ***********************
        Double_t offsetMarkerXMass2         = 0.1;
        Double_t offsetMarkerYMass2         = 0.1;
        //****************** Scale factors ******************
        Double_t scaleMarkerMass2           = 1.2;

        padMassLegend1->cd();
        //****************** first Column **************************************************
        TLatex *textMassPCM                 = new TLatex(columnsLegendMass2[0],rowsLegendMass2[1],"PCM");
        SetStyleTLatex( textMassPCM, textSizeLabelsPixel,4);
        textMassPCM->SetTextFont(43);
        textMassPCM->Draw();
        TLatex *textMassPCMEMCAL            = new TLatex(columnsLegendMass2[0],rowsLegendMass2[2],"PCM-EMCal");
        SetStyleTLatex( textMassPCMEMCAL, textSizeLabelsPixel,4);
        textMassPCMEMCAL->SetTextFont(43);
        textMassPCMEMCAL->Draw();
        TLatex *textMassEMCAL               = new TLatex(columnsLegendMass2[0],rowsLegendMass2[3],"EMCal");
        SetStyleTLatex( textMassEMCAL, textSizeLabelsPixel,4);
        textMassEMCAL->SetTextFont(43);
        textMassEMCAL->Draw();

        //****************** second Column *************************************************
        TLatex *textMassData                = new TLatex(columnsLegendMass2[1],rowsLegendMass2[0] ,"Data");
        SetStyleTLatex( textMassData, textSizeLabelsPixel,4);
        textMassData->SetTextFont(43);
        textMassData->Draw();
        TLatex *textMassMC                  = new TLatex(columnsLegendMass2[2] ,rowsLegendMass2[0],"MC");
        SetStyleTLatex( textMassMC, textSizeLabelsPixel,4);
        textMassMC->SetTextFont(43);
        textMassMC->Draw();

        TMarker* markerPCMPi0Mass        = CreateMarkerFromHisto(histoPCMPi0Mass,columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
        markerPCMPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2);
        TMarker* markerPCMEMCALPi0Mass   = CreateMarkerFromGraph(graphPCMEMCALPi0Mass,columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2 ,scaleMarkerMass2);
        markerPCMEMCALPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2);
        TMarker* markerEMCALPi0Mass      = CreateMarkerFromGraph(graphEMCALPi0Mass,columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2 ,scaleMarkerMass2);
        markerEMCALPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2);

        TMarker* markerPCMPi0MassMC      = CreateMarkerFromHisto(histoPCMPi0TrueMass,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
        markerPCMPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[1]+ offsetMarkerYMass2);
        TMarker* markerPCMEMCALPi0MassMC = CreateMarkerFromGraph(graphPCMEMCALPi0MassMC,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2 ,scaleMarkerMass2);
        markerPCMEMCALPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[2]+ offsetMarkerYMass2);
        TMarker* markerEMCALPi0MassMC    = CreateMarkerFromGraph(graphEMCALPi0MassMC,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2 ,scaleMarkerMass2);
        markerEMCALPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[3]+ offsetMarkerYMass2);

    canvasMassWidthPi0->Update();
    canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // *********************************** Mass and width for pi0 at 8TeV including PHOS *********************************
    // **********************************************************************************************************************

    canvasMassWidthPi0->cd();
    padWidthPi0->Draw();
    padMassPi0->Draw();

    TPad* padMassLegend2                = new TPad("padMassLegend2", "", 0.13, 0.32, 0.48, 0.51,-1, -1, -2);
    DrawGammaPadSettings( padMassLegend2, 0., 0., 0., 0.);
    padMassLegend2->SetFillStyle(0);
    padMassLegend2->Draw();

    padWidthPi0->cd();
    padWidthPi0->SetLogx();

        histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,32.5);
        histo2DAllPi0FWHM->DrawCopy();

        DrawGammaSetMarker(histoPCMPi0FWHMMeV, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        histoPCMPi0FWHMMeV->Draw("p,same,e");
        DrawGammaSetMarker(histoPCMPi0TrueFWHMMeV, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
        histoPCMPi0TrueFWHMMeV->Draw("p,same,e");

        DrawGammaSetMarker(histoPHOSPi0FWHMMeV, markerStyleDet[1], markerSizeDet[1]*0.55, colorDet[1] , colorDet[1]);
        histoPHOSPi0FWHMMeV->Draw("p,same,e");
        DrawGammaSetMarker(histoPHOSPi0TrueFWHMMeV, markerStyleDetMC[1], markerSizeDetMC[1]*0.55, colorDetMC[1] , colorDetMC[1]);
        histoPHOSPi0TrueFWHMMeV->Draw("p,same,e");

        graphEMCALPi0FWHM->Draw("p,same,z");
        graphEMCALPi0FWHMMC->Draw("p,same,z");

        graphPCMEMCALPi0FWHM->Draw("p,same,z");
        graphPCMEMCALPi0FWHMMC->Draw("p,same,z");

        labelMassPerf->Draw();
        labelMassEnergy->Draw();
        labelMassPi0->Draw();
        labelLegendAMass->Draw();

//     //********************************** Defintion of the Legend **************************************************
//         Double_t columnsLegendFWHM[4]       = {0.12, 0.44, 0.525, 0.39};
//         Double_t columnsLegendFWHMAbs[4]    = {4, 1.6, 2.5, 9.8};
//         Double_t rowsLegendFWHM3[5]         = {0.88, 0.8, 0.74, 0.68, 0.62};
//         Double_t rowsLegendFWHMAbs3[5]      = {0.2, 26.8, 24.8, 22.8, 20.8 };
//         //****************** Scale factors ******************
//         Double_t scaleMarkerFWHM            = 1.2;
//
//         //****************** first Column **************************************************
//         TLatex *textFWHMPCM2                = new TLatex(columnsLegendFWHM[0],rowsLegendFWHM3[1],"PCM (FWHM/2.36)");
//         SetStyleTLatex( textFWHMPCM2, textSizeLabelsPixel,4);
//         textFWHMPCM2->SetTextFont(43);
//         textFWHMPCM2->Draw();
//         TLatex *textFWHMPCMEMCAL2           = new TLatex(columnsLegendFWHM[0],rowsLegendFWHM3[2],"PCM-EMCal (FWHM/2.36)");
//         SetStyleTLatex( textFWHMPCMEMCAL2, textSizeLabelsPixel,4);
//         textFWHMPCMEMCAL2->SetTextFont(43);
//         textFWHMPCMEMCAL2->Draw();
//         TLatex *textFWHMEMCAL2              = new TLatex(columnsLegendFWHM[0],rowsLegendFWHM3[3],"EMCal (FWHM/2.36)");
//         SetStyleTLatex( textFWHMEMCAL2, textSizeLabelsPixel,4);
//         textFWHMEMCAL2->SetTextFont(43);
//         textFWHMEMCAL2->Draw();
//         TLatex *textFWHMPHOS                = new TLatex(columnsLegendFWHM[0],rowsLegendFWHM3[4],"PHOS (#sigma)");
//         SetStyleTLatex( textFWHMPHOS, textSizeLabelsPixel,4);
//         textFWHMPHOS->SetTextFont(43);
//         textFWHMPHOS->Draw();
//
//         //****************** second Column *************************************************
//         TLatex *textFWHMData2               = new TLatex(columnsLegendFWHM[1],rowsLegendFWHM3[0] ,"Data");
//         SetStyleTLatex( textFWHMData2, textSizeLabelsPixel ,4);
//         textFWHMData2->SetTextFont(43);
//         textFWHMData2->Draw();
//         TLatex *textFWHMMC2                 = new TLatex(columnsLegendFWHM[2] ,rowsLegendFWHM3[0],"MC");
//         SetStyleTLatex( textFWHMMC2, textSizeLabelsPixel,4);
//         textFWHMMC2->SetTextFont(43);
//         textFWHMMC2->Draw();
//
//         TMarker* markerPCMPi0FWHM        = CreateMarkerFromHisto(histoPCMPi0FWHMMeV,columnsLegendFWHMAbs[1] ,rowsLegendFWHMAbs3[1] ,scaleMarkerFWHM);
//         TMarker* markerPCMEMCALPi0FWHM   = CreateMarkerFromGraph(graphPCMEMCALPi0FWHM,columnsLegendFWHMAbs[1],rowsLegendFWHMAbs3[2],scaleMarkerFWHM);
//         TMarker* markerEMCALPi0FWHM      = CreateMarkerFromGraph(graphEMCALPi0FWHM,columnsLegendFWHMAbs[1],rowsLegendFWHMAbs3[3],scaleMarkerFWHM);
//         TMarker* markerPHOSPi07TeVFWHM          = CreateMarkerFromHisto(histoPHOSPi0FWHMMeV,columnsLegendFWHMAbs[1],rowsLegendFWHMAbs3[3],scaleMarkerFWHM);
//
//         TMarker* markerPCMPi0FWHMMC      = CreateMarkerFromHisto(histoPCMPi0TrueFWHMMeV,columnsLegendFWHMAbs[2],rowsLegendFWHMAbs3[1],scaleMarkerFWHM);
//         TMarker* markerPCMEMCALPi0FWHMMC = CreateMarkerFromGraph(graphPCMEMCALPi0FWHMMC,columnsLegendFWHM[2],rowsLegendFWHMAbs3[2],scaleMarkerFWHM);
//         TMarker* markerEMCALPi0FWHMMC    = CreateMarkerFromGraph(graphEMCALPi0FWHMMC,columnsLegendFWHM[2],rowsLegendFWHMAbs3[3] ,scaleMarkerFWHM);
//         TMarker* markerPHOSPi07TeVFWHMMC        = CreateMarkerFromHisto(histoPHOSPi0TrueFWHMMeV,columnsLegendFWHMAbs[2],rowsLegendFWHMAbs3[3],scaleMarkerFWHM);
//
//         markerPCMPi0FWHM->DrawMarker(columnsLegendFWHMAbs[1] ,rowsLegendFWHMAbs3[1]);
//         markerPCMEMCALPi0FWHM->DrawMarker(columnsLegendFWHMAbs[1] ,rowsLegendFWHMAbs3[2]);
//         markerEMCALPi0FWHM->DrawMarker(columnsLegendFWHMAbs[1] ,rowsLegendFWHMAbs3[3]);
//         markerPHOSPi07TeVFWHM->DrawMarker(columnsLegendFWHMAbs[1] ,rowsLegendFWHMAbs3[4]);
//         markerPCMPi0FWHMMC->DrawMarker(columnsLegendFWHMAbs[2] ,rowsLegendFWHMAbs3[1]);
//         markerPCMEMCALPi0FWHMMC->DrawMarker(columnsLegendFWHMAbs[2],rowsLegendFWHMAbs3[2]);
//         markerEMCALPi0FWHMMC->DrawMarker(columnsLegendFWHMAbs[2],rowsLegendFWHMAbs3[3]);
//         markerPHOSPi07TeVFWHMMC->DrawMarker(columnsLegendFWHMAbs[2] ,rowsLegendFWHMAbs3[4]);


    padMassPi0->cd();
    padMassPi0->SetLogx();

        histo2DAllPi0Mass->DrawCopy();

        histoPCMPi0Mass->Draw("p,same,e");
        histoPCMPi0TrueMass->Draw("p,same,e");

        DrawGammaSetMarker(histoPHOSPi0Mass, markerStyleDet[1], markerSizeDet[1]*0.55, colorDet[1] , colorDet[1]);
        histoPHOSPi0Mass->Draw("p,same,e");
        DrawGammaSetMarker(histoPHOSPi0TrueMass, markerStyleDetMC[1], markerSizeDetMC[1]*0.55, colorDetMC[1] , colorDetMC[1]);
        histoPHOSPi0TrueMass->Draw("p,same,e");

        graphEMCALPi0Mass->Draw("p,same,z");
        graphEMCALPi0MassMC->Draw("p,same,z");

        graphPCMEMCALPi0Mass->Draw("p,same,z");
        graphPCMEMCALPi0MassMC->Draw("p,same,z");


        DrawGammaLines(0.23, 50. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);

        labelLegendBMass->Draw();

        //********************************** Defintion of the Legend **************************************************
        Double_t rowsLegendMass3[5]         = {0.8,0.6,0.4,0.2,0.01};
        Double_t offsetMarkerYMass3         = 0.06;

        padMassLegend2->cd();
        //****************** first Column **************************************************
        TLatex *textMassPCM2                = new TLatex(columnsLegendMass2[0],rowsLegendMass3[1],"PCM");
        SetStyleTLatex( textMassPCM2, textSizeLabelsPixel,4);
        textMassPCM2->SetTextFont(43);
        textMassPCM2->Draw();
        TLatex *textMassPCMEMCAL2           = new TLatex(columnsLegendMass2[0],rowsLegendMass3[2],"PCM-EMCal");
        SetStyleTLatex( textMassPCMEMCAL2, textSizeLabelsPixel,4);
        textMassPCMEMCAL2->SetTextFont(43);
        textMassPCMEMCAL2->Draw();
        TLatex *textMassEMCAL2              = new TLatex(columnsLegendMass2[0],rowsLegendMass3[3],"EMCal");
        SetStyleTLatex( textMassEMCAL2, textSizeLabelsPixel,4);
        textMassEMCAL2->SetTextFont(43);
        textMassEMCAL2->Draw();
        TLatex *textMassPHOS2               = new TLatex(columnsLegendMass2[0],rowsLegendMass3[4],"PHOS");
        SetStyleTLatex( textMassPHOS2, textSizeLabelsPixel,4);
        textMassPHOS2->SetTextFont(43);
        textMassPHOS2->Draw();

        //****************** second Column *************************************************
        TLatex *textMassData2               = new TLatex(columnsLegendMass2[1],rowsLegendMass3[0] ,"Data");
        SetStyleTLatex( textMassData2, textSizeLabelsPixel,4);
        textMassData2->SetTextFont(43);
        textMassData2->Draw();
        TLatex *textMassMC2                 = new TLatex(columnsLegendMass2[2] ,rowsLegendMass3[0],"MC");
        SetStyleTLatex( textMassMC2, textSizeLabelsPixel,4);
        textMassMC2->SetTextFont(43);
        textMassMC2->Draw();

        markerPCMPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass3[1]+ offsetMarkerYMass3);
        markerPCMEMCALPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass3[2]+ offsetMarkerYMass3);
        markerEMCALPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass3[3]+ offsetMarkerYMass3);
        TMarker* markerPHOSPi0Mass   = CreateMarkerFromHisto(histoPHOSPi0Mass,columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass3[4]+ offsetMarkerYMass3 ,scaleMarkerMass2);
        markerPHOSPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass3[4]+ offsetMarkerYMass3);

        markerPCMPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass3[1]+ offsetMarkerYMass3);
        markerPCMEMCALPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass3[2]+ offsetMarkerYMass3);
        markerEMCALPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass3[3]+ offsetMarkerYMass3);
        TMarker* markerPHOSPi0MassMC = CreateMarkerFromHisto(histoPHOSPi0TrueMass,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass3[4]+ offsetMarkerYMass3 ,scaleMarkerMass2);
        markerPHOSPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass3[4]+ offsetMarkerYMass3);

    canvasMassWidthPi0->Update();
    canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth_incPHOS.%s",outputDir.Data(),suffix.Data()));


    // **********************************************************************************************************************
    // ******************************** Cross section for pi0 single measurement 8TeV ************************************
    // **********************************************************************************************************************

    TCanvas* canvasXSectionPi0  = new TCanvas("canvasXSectionPi0","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasXSectionPi0, 0.14, 0.02, 0.02, 0.09);
    canvasXSectionPi0->SetLogx();
    canvasXSectionPi0->SetLogy();

    TH2F * histo2DXSectionPi0;
    histo2DXSectionPi0          = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,0.23,110.,1000,1e-2,9e11);
    SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.45);
    histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
    histo2DXSectionPi0->GetXaxis()->SetLabelOffset(-0.01);
    histo2DXSectionPi0->Draw("copy");


        DrawGammaSetMarkerTGraphAsym(graphPHOSPi0InvXSectionSys, markerStyleDet[1], markerSizeDet[1]*0.75, colorDet[1] , colorDet[1], widthLinesBoxes, kTRUE);
        graphPHOSPi0InvXSectionSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvXSectionSys, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0], widthLinesBoxes, kTRUE);
        graphPCMPi0InvXSectionSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0InvXSectionSys, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2], widthLinesBoxes, kTRUE);
        graphEMCALPi0InvXSectionSys->Draw("E2same");
//      graphPCMEMCALPi0InvXSectionSys->RemovePoint(graphPCMEMCALPi0InvXSectionSys->GetN()-1);
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0InvXSectionSys, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4], widthLinesBoxes, kTRUE);
        graphPCMEMCALPi0InvXSectionSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphEMCALMergedPi0InvXSectionSys, markerStyleDet[9], markerSizeDet[9]*0.75, colorDet[9] , colorDet[9], widthLinesBoxes, kTRUE);
        graphEMCALMergedPi0InvXSectionSys->Draw("E2same");

        DrawGammaSetMarker(histoPHOSPi0InvXSectionStat, markerStyleDet[1], markerSizeDet[1]*0.75, colorDet[1] , colorDet[1]);
        histoPHOSPi0InvXSectionStat->Draw("p,same,e");
//      DrawGammaSetMarker(histoPCMPi0InvXSectionStat, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
//      histoPCMPi0InvXSectionStat->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvXSectionStat,markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
        graphPCMPi0InvXSectionStat->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0InvXSectionStat, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
        graphEMCALPi0InvXSectionStat->Draw("p,same,z");
//      DrawGammaSetMarker(histoEMCALPi0InvXSectionStat, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
//      histoEMCALPi0InvXSectionStat->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0InvXSectionStat, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4]);
        graphPCMEMCALPi0InvXSectionStat->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphEMCALMergedPi0InvXSectionStat, markerStyleDet[9], markerSizeDet[9]*0.75, colorDet[9] , colorDet[9]);
        graphEMCALMergedPi0InvXSectionStat->Draw("p,same,z");


        TLatex *labelEnergyXSectionPi0      = new TLatex(0.64,0.92,collisionSystem8TeV.Data());
        SetStyleTLatex( labelEnergyXSectionPi0, 0.035,4);
        labelEnergyXSectionPi0->Draw();
        TLatex *labelDetSysXSectionPi0      = new TLatex(0.64,0.88,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPi0, 0.035,4);
        labelDetSysXSectionPi0->Draw();

        TLegend* legendXSectionPi0          = new TLegend(0.62,0.62,0.9,0.86);
        legendXSectionPi0->SetFillColor(0);
        legendXSectionPi0->SetLineColor(0);
        legendXSectionPi0->SetTextFont(42);
        legendXSectionPi0->SetTextSize(0.035);
        legendXSectionPi0->AddEntry(graphPCMPi0InvXSectionSys,"PCM","fp");
        legendXSectionPi0->AddEntry(graphPHOSPi0InvXSectionSys,"PHOS","fp");
        legendXSectionPi0->AddEntry(graphEMCALPi0InvXSectionSys,"EMCal","fp");
        legendXSectionPi0->AddEntry(graphPCMEMCALPi0InvXSectionSys,"PCM-EMCal","fp");
        legendXSectionPi0->AddEntry(graphEMCALMergedPi0InvXSectionSys,"EMCal merged","fp");
        legendXSectionPi0->Draw();

    canvasXSectionPi0->SaveAs(Form("%s/Pi0_InvXSectionCompAllSystems.%s",outputDir.Data(),suffix.Data()));

    canvasXSectionPi0->cd();
    histo2DXSectionPi0->Draw("copy");

        graphPHOSPi0InvXSectionSys->Draw("E2same");
        graphPCMPi0InvXSectionSys->Draw("E2same");
        graphEMCALPi0InvXSectionSys->Draw("E2same");
        graphPCMEMCALPi0InvXSectionSys->Draw("E2same");
        graphEMCALMergedPi0InvXSectionSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionSysA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphCombPi0InvXSectionSysA->Draw("E2same");

        histoPHOSPi0InvXSectionStat->Draw("p,same,e");
        graphPCMPi0InvXSectionStat->Draw("p,same,z");
        graphEMCALPi0InvXSectionStat->Draw("p,same,z");
        graphPCMEMCALPi0InvXSectionStat->Draw("p,same,z");
        graphEMCALMergedPi0InvXSectionStat->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatA, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombPi0InvXSectionStatA->Draw("p,same,z");

        labelEnergyXSectionPi0->Draw();
        labelDetSysXSectionPi0->Draw();

        legendXSectionPi0->AddEntry(graphCombPi0InvXSectionSysA,"comb","fp");
        legendXSectionPi0->Draw();

//      DrawGammaSetMarkerTGraphAsym(graphChargedHadronsStatPP, markerStyleDet[1], markerSizeDet[1]*0.75, kBlack , kBlack, widthLinesBoxes);
//      graphChargedHadronsStatPP->Draw("E2same");

   canvasXSectionPi0->SaveAs(Form("%s/Pi0_InvXSectionCompAllSystems_Comb.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************** Acceptance * Efficiency for pi0 single measurement 8TeV **************************
    // **********************************************************************************************************************
    textSizeLabelsPixel             = 55;
    Double_t textSizeLabelsRel      = 55./1200;
    cout << textSizeLabelsRel << endl;

    TCanvas* canvasAcceptanceTimesEff       = new TCanvas("canvasAcceptanceTimesEff", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasAcceptanceTimesEff,  0.1, 0.01, 0.015, 0.095);
    canvasAcceptanceTimesEff->SetLogy(1);
    canvasAcceptanceTimesEff->SetLogx(1);

        TH2F * histo2DAccEff;
        histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, 0.23,  110, 1000, 8e-5, 3 );
        SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec} / #it{P}"),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);//(#times #epsilon_{pur})
        histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
        histo2DAccEff->GetXaxis()->SetLabelOffset(-0.01);
        histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAccEff->DrawCopy();

        DrawGammaSetMarker(histoPCMPi0AccTimesEff, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        histoPCMPi0AccTimesEff->Draw("p,same,e");

        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0AccTimesEff, markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
        graphPCMEMCALPi0AccTimesEff->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0AccTimesEff, markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
        graphEMCALPi0AccTimesEff->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphEMCALMergedPi0AccTimesEff, markerStyleDet[9], markerSizeDet[9]*0.55, colorDet[9] , colorDet[9]);
        graphEMCALMergedPi0AccTimesEff->Draw("p,same,e");

        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.13, 0.83, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffiAccPi0->AddEntry(histoPCMPi0AccTimesEff,"PCM","p");
        legendEffiAccPi0->AddEntry(graphPCMEMCALPi0AccTimesEff,"PCM-EMCal","p");
        legendEffiAccPi0->AddEntry(graphEMCALPi0AccTimesEff,"EMCal","p");
        legendEffiAccPi0->AddEntry(graphEMCALMergedPi0AccTimesEff,"EMCal, merged","p");
        legendEffiAccPi0->Draw();

        TLatex *labelPerfEffi               = new TLatex(0.15,0.92,"ALICE performance");
        SetStyleTLatex( labelPerfEffi, textSizeLabelsRel,4);
        labelPerfEffi->Draw();
        TLatex *labelEnergyEffi             = new TLatex(0.15,0.87,collisionSystem8TeV.Data());
        SetStyleTLatex( labelEnergyEffi, textSizeLabelsRel,4);
        labelEnergyEffi->Draw();
        TLatex *labelDetSysEffiPi0          = new TLatex(0.15,0.82,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysEffiPi0, textSizeLabelsRel,4);
        labelDetSysEffiPi0->Draw();


    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ************************* Acceptance * Efficiency for pi0 single measurement 8TeV inc PHOS ************************
    // **********************************************************************************************************************

    canvasAcceptanceTimesEff->cd();
        histo2DAccEff->DrawCopy();

        DrawGammaSetMarker(histoPHOSPi0AccTimesEff, markerStyleDet[1], markerSizeDet[1]*0.55, colorDet[1] , colorDet[1]);
        histoPHOSPi0AccTimesEff->Draw("p,same,z");

        histoPCMPi0AccTimesEff->Draw("p,same,e");
        graphPCMEMCALPi0AccTimesEff->Draw("p,same,z");
        graphEMCALPi0AccTimesEff->Draw("p,same,z");
        graphEMCALMergedPi0AccTimesEff->Draw("p,same,e");

        TLegend* legendEffiAccPi02          = GetAndSetLegend2(0.55, 0.13, 0.83, 0.13+(5*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffiAccPi02->AddEntry(histoPCMPi0AccTimesEff,"PCM","p");
        legendEffiAccPi02->AddEntry(graphPCMEMCALPi0AccTimesEff,"PCM-EMCal","p");
        legendEffiAccPi02->AddEntry(graphEMCALPi0AccTimesEff,"EMCal","p");
        legendEffiAccPi02->AddEntry(graphEMCALMergedPi0AccTimesEff,"mEMC","p");
        legendEffiAccPi02->AddEntry(histoPHOSPi0AccTimesEff,"PHOS","p");
        legendEffiAccPi02->Draw();

        labelPerfEffi->Draw();
        labelEnergyEffi->Draw();
        labelDetSysEffiPi0->Draw();

    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff_incPHOS.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************************* Comparison to theory calculations Pi0 ************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                     = 48;

    TCanvas* canvasRatioPP                  = new TCanvas("canvasRatioPP","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioPP,  0.12, 0.01, 0.01, 0.11);
    canvasRatioPP->cd();
    canvasRatioPP->SetLogx();

        textsizeLabelsPP                    = 0;
        textsizeFacPP                       = 0;
        if (canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) <canvasRatioPP->YtoPixel(canvasRatioPP->GetY1()) ){
            textsizeLabelsPP                = (Double_t)textSizeLabelsPixel/canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) ;
            textsizeFacPP                   = (Double_t)1./canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) ;
        } else {
        textsizeLabelsPP                    = (Double_t)textSizeLabelsPixel/canvasRatioPP->YtoPixel(canvasRatioPP->GetY1());
        textsizeFacPP                       = (Double_t)1./canvasRatioPP->YtoPixel(canvasRatioPP->GetY1());
        }
        cout << textsizeLabelsPP << endl;

        TH2F * ratio2DTheoryPP       = new TH2F("ratio2DTheoryPP","ratio2DTheoryPP",1000,0.23,110.,1000,0.4,2.85);
        SetStyleHistoTH2ForGraphs(ratio2DTheoryPP, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                                  0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
        ratio2DTheoryPP->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio2DTheoryPP->GetYaxis()->SetNdivisions(505);
        ratio2DTheoryPP->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DTheoryPP->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DTheoryPP->GetXaxis()->SetNoExponent(kTRUE);
//      ratio2DTheoryPP->GetYaxis()->SetLabelOffset(0.01);
        ratio2DTheoryPP->GetXaxis()->SetLabelFont(42);
        ratio2DTheoryPP->GetYaxis()->SetLabelFont(42);
//      ratio2DTheoryPP->GetXaxis()->SetLabelOffset(-0.01);
    //    ratio2DTheoryPP->GetXaxis()->SetTickLength(0.07);
        ratio2DTheoryPP->DrawCopy();

//        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS14, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
//        graphRatioPi0CombNLODSS14->Draw("2,same");

//        DrawGammaNLOTGraph( graphRatioPi0CombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, kGray+1);
//        graphRatioPi0CombNLOMuHalf->Draw("same,c");
//        DrawGammaNLOTGraph( graphRatioPi0CombNLOMuOne, widthCommonFit, styleLineNLOMuOne, kGray+1);
//        graphRatioPi0CombNLOMuOne->Draw("same,c");
//        DrawGammaNLOTGraph( graphRatioPi0CombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, kGray+2);
//        graphRatioPi0CombNLOMuTwo->Draw("same,c");
        graphRatioPi0DSS14->SetLineWidth(widthCommonFit);
        graphRatioPi0DSS14->SetLineColor(colorNLO);
        graphRatioPi0DSS14->SetLineStyle(1);
        graphRatioPi0DSS14->SetFillStyle(1001);
        graphRatioPi0DSS14->SetFillColor(colorNLO);
        graphRatioPi0DSS14->Draw("same,e4");

        DrawGammaSetMarker(histoRatioPythia8ToFit, 24, 1.5, kRed+2 , kRed+2);
        histoRatioPythia8ToFit->SetLineWidth(widthCommonFit);
        histoRatioPythia8ToFit->GetXaxis()->SetRangeUser(0.3,35);
        histoRatioPythia8ToFit->Draw("same,hist,c");
        DrawGammaSetMarker(histoRatioPythia8_4CToFit, 24, 1.5, kGreen+3 , kGreen+3);
        histoRatioPythia8_4CToFit->SetLineWidth(widthCommonFit);
        histoRatioPythia8_4CToFit->SetLineStyle(2);
        histoRatioPythia8_4CToFit->GetXaxis()->SetRangeUser(0.3,35);
        histoRatioPythia8_4CToFit->Draw("same,hist,c");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatA->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA->Draw("2,same");
        graphRatioPi0CombCombFitStatA->Draw("p,same");

        TBox* boxErrorSigmaRatio = CreateBoxConv(kGray+3, 0.3, 1.-(0.0196 ), 0.35, 1.+(0.0196));
        boxErrorSigmaRatio->Draw();
        DrawGammaLines(0.23, 110.,1., 1.,0.1,kGray);

        TLegend* legendRatioTheorypp_3Parted= GetAndSetLegend2(0.35,0.71,0.6,0.97, 0.85* textSizeLabelsPixel);
    //    legendRatioTheorypp_3Parted->SetNColumns(2);
        legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0CombCombFitSysA,"ALICE #pi^{0}","pf");
//        legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0CombNLOMuHalf, "NLO, DSS07 #mu = 0.5 #it{p}_{T}", "l");
//        legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0CombNLOMuOne,  "NLO, DSS07 #mu = #it{p}_{T}", "l");
//        legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0CombNLOMuTwo,  "NLO, DSS07 #mu = 2 #it{p}_{T}", "l");
        legendRatioTheorypp_3Parted->AddEntry(histoRatioPythia8ToFit,  "PYTHIA 8.2, Monash 2013", "l");
        legendRatioTheorypp_3Parted->AddEntry(histoRatioPythia8_4CToFit,  "PYTHIA 8.2, Tune 4C", "l");
        legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0DSS14,  "NLO, PDF:MSTW08 - FF:DSS14", "f");
        legendRatioTheorypp_3Parted->Draw();

        TLatex *labelRatioTheoryPPP   = new TLatex(0.418,0.68,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}");
        SetStyleTLatex( labelRatioTheoryPPP, 0.85*textsizeLabelsPP,4);
        labelRatioTheoryPPP->Draw();

        TLatex *labelRatioTheoryPP   = new TLatex(0.15,0.925,collisionSystem8TeV.Data());
        SetStyleTLatex( labelRatioTheoryPP, 0.85*textsizeLabelsPP,4);
        labelRatioTheoryPP->Draw();


    canvasRatioPP->Update();
    canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PP.%s",outputDir.Data(),suffix.Data()));


    TH2F * ratio2DTheoryPP2       = new TH2F("ratio2DTheoryPP2","ratio2DTheoryPP2",1000,0.23,110.,1000,0.2,5.2);
    SetStyleHistoTH2ForGraphs(ratio2DTheoryPP2, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
    ratio2DTheoryPP2->GetYaxis()->SetMoreLogLabels(kTRUE);
    ratio2DTheoryPP2->GetYaxis()->SetNdivisions(505);
    ratio2DTheoryPP2->GetYaxis()->SetNoExponent(kTRUE);
    ratio2DTheoryPP2->GetXaxis()->SetMoreLogLabels(kTRUE);
    ratio2DTheoryPP2->GetXaxis()->SetNoExponent(kTRUE);
//      ratio2DTheoryPP2->GetYaxis()->SetLabelOffset(0.01);
    ratio2DTheoryPP2->GetXaxis()->SetLabelFont(42);
    ratio2DTheoryPP2->GetYaxis()->SetLabelFont(42);
//      ratio2DTheoryPP2->GetXaxis()->SetLabelOffset(-0.01);
//    ratio2DTheoryPP2->GetXaxis()->SetTickLength(0.07);
    ratio2DTheoryPP2->DrawCopy();

    graphRatioPi0DSS14->SetLineWidth(widthCommonFit);
    graphRatioPi0DSS14->SetLineColor(colorNLO);
    graphRatioPi0DSS14->SetLineStyle(1);
    graphRatioPi0DSS14->SetFillStyle(1001);
    graphRatioPi0DSS14->SetFillColor(colorNLO);
    graphRatioPi0DSS14->Draw("same,e4");

    DrawGammaNLOTGraph( graphRatioPi0CombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, kGray+1);
    graphRatioPi0CombNLOMuHalf->Draw("same,c");
    DrawGammaNLOTGraph( graphRatioPi0CombNLOMuOne, widthCommonFit, styleLineNLOMuOne, kGray+1);
    graphRatioPi0CombNLOMuOne->Draw("same,c");
    DrawGammaNLOTGraph( graphRatioPi0CombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, kGray+2);
    graphRatioPi0CombNLOMuTwo->Draw("same,c");

    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
    graphRatioPi0CombCombFitStatA->SetLineWidth(widthLinesBoxes);
    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
    graphRatioPi0CombCombFitSysA->SetLineWidth(0);
    graphRatioPi0CombCombFitSysA->Draw("2,same");
    graphRatioPi0CombCombFitStatA->Draw("p,same");

    boxErrorSigmaRatio->Draw();
    DrawGammaLines(0.23, 110.,1., 1.,0.1,kGray);

    TLegend* legendRatioTheorypp_3Parted2= GetAndSetLegend2(0.15,0.7,0.4,0.96, 0.85* textSizeLabelsPixel);
    legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0CombCombFitSysA,"ALICE #pi^{0}","pf");
    legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0CombNLOMuHalf, "NLO, PDF:CTEQ6M5 - FF:DSS07 #mu = 0.5 #it{p}_{T}", "l");
    legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0CombNLOMuOne,  "NLO, PDF:CTEQ6M5 - FF:DSS07 #mu = #it{p}_{T}", "l");
    legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0CombNLOMuTwo,  "NLO, PDF:CTEQ6M5 - FF:DSS07 #mu = 2 #it{p}_{T}", "l");
    legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0DSS14,  "NLO, PDF:MSTW08 - FF:DSS14", "f");
    legendRatioTheorypp_3Parted2->Draw();

    TLatex *labelRatioTheoryPPP2   = new TLatex(0.218,0.67,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}");
    SetStyleTLatex( labelRatioTheoryPPP2, 0.85*textsizeLabelsPP,4);
    labelRatioTheoryPPP2->Draw();

    TLatex *labelRatioTheoryPP22   = new TLatex(0.78,0.925,collisionSystem8TeV.Data());
    SetStyleTLatex( labelRatioTheoryPP22, 0.85*textsizeLabelsPP,4);
    labelRatioTheoryPP22->Draw();

    canvasRatioPP->Update();
    canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PP2.%s",outputDir.Data(),suffix.Data()));

    //*************************************************************************************************************
    //***************************** Paper plot X-section and ratios ***********************************************
    //*************************************************************************************************************

    Double_t arrayBoundariesX1_XSec[2];
    Double_t arrayBoundariesY1_XSec[6];
    Double_t relativeMarginsXXSec[3];
    Double_t relativeMarginsYXSec[3];
    textSizeLabelsPixel = 48;
    ReturnCorrectValuesForCanvasScaling(1250,2000, 1, 5,0.135, 0.005, 0.003,0.05,arrayBoundariesX1_XSec,arrayBoundariesY1_XSec,relativeMarginsXXSec,relativeMarginsYXSec);

    TCanvas* canvasInvSectionPaper      = new TCanvas("canvasInvSectionPaper","",0,0,1250,2000);  // gives the page size
    DrawGammaCanvasSettings( canvasInvSectionPaper,  0.13, 0.02, 0.03, 0.06);

    TPad* padInvSectionSpec             = new TPad("padInvSectionSpec", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[3], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[0],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionSpec, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[0], relativeMarginsYXSec[1]);
    padInvSectionSpec->Draw();
    Double_t marginXSec                 = relativeMarginsXXSec[0]*1250;
    Double_t textsizeLabelsXSecUp       = 0;
    Double_t textsizeFacXSecUp          = 0;
    if (padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) < padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1())){
        textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
        textsizeFacXSecUp               = (Double_t)1./padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
    } else {
        textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
        textsizeFacXSecUp               = (Double_t)1./padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
    }

    TPad* padInvSectionNLORatio         = new TPad("padInvSectionNLORatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[4], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[3],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionNLORatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[1]);
    padInvSectionNLORatio->Draw();
    Double_t textsizeLabelsXSecMiddle   = 0;
    Double_t textsizeFacXSecMiddle      = 0;
    if (padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) < padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1())){
        textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) ;
        textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) ;
    } else {
        textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1());
        textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1());
    }

    TPad* padInvSectionPythiaRatio      = new TPad("padInvSectionPythiaRatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[5], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[4],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionPythiaRatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[2]);
    padInvSectionPythiaRatio->Draw();
    Double_t textsizeLabelsXSecDown     = 0;
    Double_t textsizeFacXSecDown        = 0;
    if (padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) < padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1())){
        textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) ;
        textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) ;
    } else {
        textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1());
        textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1());
    }


    padInvSectionSpec->cd();
    padInvSectionSpec->SetLogy(1);
    padInvSectionSpec->SetLogx(1);
        SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",
                                0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.2/(textsizeFacXSecUp*marginXSec));
        histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
        histo2DXSectionPi0->GetXaxis()->SetLabelOffset(+0.01);
        histo2DXSectionPi0->Draw();


        DrawGammaSetMarker(histoPythia8InvXSection, 24, 1.5, kRed+2 , kRed+2);
        histoPythia8InvXSection->SetLineWidth(widthCommonFit);
        histoPythia8InvXSection->GetXaxis()->SetRangeUser(0.3,35);
        histoPythia8InvXSection->Draw("same,hist,c");

        DrawGammaSetMarker(histoPythia8_4CInvXSection, 24, 1.5, kGreen+3 , kGreen+3);
        histoPythia8_4CInvXSection->SetLineWidth(widthCommonFit);
        histoPythia8_4CInvXSection->GetXaxis()->SetRangeUser(0.3,35);
        histoPythia8_4CInvXSection->SetLineStyle(2);
        histoPythia8_4CInvXSection->Draw("same,hist,c");

//        graphNLODSS14Calc->RemovePoint(0);
//        DrawGammaSetMarkerTGraphAsym(graphNLODSS14Calc, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
//        graphNLODSS14Calc->Draw("3,same");

//        DrawGammaNLOTGraph( graphNLOCalcPi0MuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
//        graphNLOCalcPi0MuHalf->Draw("same,c");
//        DrawGammaNLOTGraph( graphNLOCalcPi0MuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
//        graphNLOCalcPi0MuOne->Draw("same,c");
//        DrawGammaNLOTGraph( graphNLOCalcPi0MuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
//        graphNLOCalcPi0MuTwo->Draw("same,c");
        graphPi0DSS14->SetLineWidth(widthCommonFit);
        graphPi0DSS14->SetLineColor(colorNLO);
        graphPi0DSS14->SetLineStyle(1);
        graphPi0DSS14->SetFillStyle(1001);
        graphPi0DSS14->SetFillColor(colorNLO);
        graphPi0DSS14->Draw("same,e4");


        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombPi0InvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatA, markerStyleComb, markerSizeComb, kBlack, kBlack);
        graphCombPi0InvXSectionStatA->Draw("p,same,z");

        DrawGammaSetMarkerTF1( fitTCMInvXSectionPi0Plot, 7, 2, kGray+2);
        fitTCMInvXSectionPi0Plot->Draw("same");

//         DrawGammaSetMarkerTF1( fitPowInvXSectionPi0, 7, 2, kBlue-7);
//         fitPowInvXSectionPi0->Draw("same");


        TLatex *labelEnergyXSectionPaper= new TLatex(0.72, 0.91, collisionSystem8TeV.Data());
        SetStyleTLatex( labelEnergyXSectionPaper, textsizeLabelsXSecUp,4);
        labelEnergyXSectionPaper->Draw();
        TLatex *labelALICEXSectionPaper= new TLatex(0.848,0.87,"ALICE");
        SetStyleTLatex( labelALICEXSectionPaper, textsizeLabelsXSecUp,4);
        labelALICEXSectionPaper->Draw();
        TLatex *labelDetSysXSectionPaper= new TLatex(0.824,0.83,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPaper, textsizeLabelsXSecUp,4);
        labelDetSysXSectionPaper->Draw();

        TLegend* legendXsectionPaper    = GetAndSetLegend2(0.17, 0.18, 0.5, 0.28+0.05*4, textSizeLabelsPixel);
        legendXsectionPaper->SetNColumns(1);
        legendXsectionPaper->SetMargin(0.2);
        legendXsectionPaper->AddEntry(graphCombPi0InvXSectionSysA,"Data","pf");
        legendXsectionPaper->AddEntry((TObject*)0,"norm. unc. 2.9%","");
//        legendXsectionPaper->AddEntry(graphNLODSS14Calc,"NLO, DSS14 ","f");
//        legendXsectionPaper->AddEntry(graphNLOCalcPi0MuHalf, "NLO, DSS07 #mu = 0.5 #it{p}_{T}", "l");
//        legendXsectionPaper->AddEntry(graphNLOCalcPi0MuOne,  "NLO, DSS07 #mu = #it{p}_{T}", "l");
//        legendXsectionPaper->AddEntry(graphNLOCalcPi0MuTwo,  "NLO, DSS07 #mu = 2 #it{p}_{T}", "l");
        legendXsectionPaper->AddEntry(histoPythia8_4CInvXSection,"PYTHIA 8.2, Tune 4C","l");
        legendXsectionPaper->AddEntry(histoPythia8InvXSection,"PYTHIA 8.2, Monash 2013","l");
        legendXsectionPaper->AddEntry(graphPi0DSS14,  "NLO, PDF:MSTW08 - FF:DSS14", "f");
        legendXsectionPaper->Draw();

        TLatex *labelRatioTheoryPP_Paper   = new TLatex(0.24,0.155,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}");
        SetStyleTLatex( labelRatioTheoryPP_Paper, 0.8*textsizeLabelsPP,4);
        labelRatioTheoryPP_Paper->Draw();

        TLegend* legendXsectionPaper2     = GetAndSetLegend2(0.17, 0.06, 0.5, 0.12, textSizeLabelsPixel);
        legendXsectionPaper2->SetMargin(0.2);
        legendXsectionPaper2->AddEntry(fitTCMInvXSectionPi0Plot,"#it{A}_{e} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{n}","l");
        legendXsectionPaper2->Draw();


    padInvSectionNLORatio->cd();
    padInvSectionNLORatio->SetLogx(1);
        TH2F * ratio2DNLO               = new TH2F("ratio2DNLO","ratio2DNLO",1000,0.23,110.,1000,0.55,1.95);
        SetStyleHistoTH2ForGraphs(ratio2DNLO, "#it{p}_{T} (GeV/#it{c})","#frac{NLO, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                  0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
        ratio2DNLO->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio2DNLO->GetYaxis()->SetNdivisions(505);
        ratio2DNLO->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DNLO->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DNLO->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DNLO->GetXaxis()->SetLabelFont(42);
        ratio2DNLO->GetYaxis()->SetLabelFont(42);
        ratio2DNLO->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DNLO->GetXaxis()->SetTickLength(0.07);
        ratio2DNLO->DrawCopy();

//        graphRatioPi0CombNLODSS14->RemovePoint(0);
//        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS14, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
//        graphRatioPi0CombNLODSS14->Draw("3,same");
//        DrawGammaNLOTGraph( graphRatioPi0CombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
//        graphRatioPi0CombNLOMuHalf->Draw("same,c");
//        DrawGammaNLOTGraph( graphRatioPi0CombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
//        graphRatioPi0CombNLOMuOne->Draw("same,c");
//        DrawGammaNLOTGraph( graphRatioPi0CombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
//        graphRatioPi0CombNLOMuTwo->Draw("same,c");
        graphRatioPi0DSS14->SetLineWidth(widthCommonFit);
        graphRatioPi0DSS14->SetLineColor(colorNLO);
        graphRatioPi0DSS14->SetLineStyle(1);
        graphRatioPi0DSS14->SetFillStyle(1001);
        graphRatioPi0DSS14->SetFillColor(colorNLO);
        graphRatioPi0DSS14->Draw("same,e4");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatA->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA->Draw("2,same");
        graphRatioPi0CombCombFitStatA->Draw("p,same");

        boxErrorSigmaRatio->Draw();
        DrawGammaLines(0.23, 110.,1., 1.,0.1,kGray);

    padInvSectionPythiaRatio->cd();
    padInvSectionPythiaRatio->SetLogx(1);
        TH2F * ratio2DPythia            = new TH2F("ratio2DPythia","ratio2DPythia",1000,0.23,110.,1000,0.55,1.95);
        SetStyleHistoTH2ForGraphs(ratio2DPythia, "#it{p}_{T} (GeV/#it{c})","#frac{Pythia, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                                  0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.2/(textsizeFacXSecDown*marginXSec), 510, 505);
        ratio2DPythia->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio2DPythia->GetYaxis()->SetNdivisions(505);
        ratio2DPythia->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DPythia->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DPythia->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DPythia->GetXaxis()->SetLabelFont(42);
        ratio2DPythia->GetYaxis()->SetLabelFont(42);
        ratio2DPythia->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DPythia->GetXaxis()->SetTickLength(0.06);
        ratio2DPythia->GetYaxis()->SetTickLength(0.04);
        ratio2DPythia->DrawCopy();

        DrawGammaSetMarker(histoRatioPythia8ToFit, 24, 1.5, kRed+2 , kRed+2);
        histoRatioPythia8ToFit->SetLineWidth(widthCommonFit);
        histoRatioPythia8ToFit->GetXaxis()->SetRangeUser(0.3,35);
        histoRatioPythia8ToFit->Draw("same,hist,c");

        DrawGammaSetMarker(histoRatioPythia8_4CToFit, 24, 1.5, kGreen+3 , kGreen+3);
        histoRatioPythia8_4CToFit->SetLineWidth(widthCommonFit);
        histoRatioPythia8_4CToFit->SetLineStyle(2);
        histoRatioPythia8_4CToFit->GetXaxis()->SetRangeUser(0.3,35);
        histoRatioPythia8_4CToFit->Draw("same,hist,c");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatA->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA->Draw("2,same");
        graphRatioPi0CombCombFitStatA->Draw("p,same");

        boxErrorSigmaRatio->Draw();
        DrawGammaLines(0.23, 110.,1., 1.,0.1,kGray);

    canvasInvSectionPaper->Print(Form("%s/Pi0_InvXSectionWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));

    // ***************************************************************************************************************
    // ************************** Comparison pi0/pi+-, pi0 pp 8TeV ***************************************************
    // ***************************************************************************************************************

    Double_t p9093_d5x1y1_xval[] = { 0.548, 0.648, 0.7481, 0.8481, 0.9481, 1.048, 1.148, 1.248, 1.348,
                                     1.448, 1.548, 1.648, 1.749, 1.849, 1.949, 2.049, 2.149, 2.249, 2.349,
                                     2.449, 2.617, 2.867, 3.118, 3.368, 3.725, 4.226, 4.729, 5.426, 6.437,
                                     7.444, 8.455, 9.462, 11.7, 16.89, 23.47, 36.99 };
    Double_t p9093_d5x1y1_xerrminus[] = { 0.04800000000000004, 0.04800000000000004, 0.04810000000000003, 0.04809999999999992, 0.04810000000000003, 0.04800000000000004, 0.04799999999999982, 0.04800000000000004, 0.04800000000000004,
                                          0.04800000000000004, 0.04800000000000004, 0.04799999999999982, 0.049000000000000155, 0.04899999999999993, 0.049000000000000155, 0.04899999999999993, 0.04899999999999993, 0.04899999999999993, 0.04900000000000038,
                                          0.04899999999999993, 0.11699999999999999, 0.11699999999999999, 0.11799999999999988, 0.11799999999999988, 0.2250000000000001, 0.22599999999999998, 0.2290000000000001, 0.42600000000000016, 0.4370000000000003,
                                          0.44399999999999995, 0.45500000000000007, 0.46199999999999974, 1.6999999999999993, 1.8900000000000006, 3.469999999999999, 6.990000000000002 };
    Double_t p9093_d5x1y1_xerrplus[] = { 0.051999999999999935, 0.051999999999999935, 0.05190000000000006, 0.05190000000000006, 0.051899999999999946, 0.052000000000000046, 0.052000000000000046, 0.052000000000000046, 0.051999999999999824,
                                         0.052000000000000046, 0.052000000000000046, 0.052000000000000046, 0.050999999999999934, 0.050999999999999934, 0.050999999999999934, 0.051000000000000156, 0.051000000000000156, 0.05099999999999971, 0.05099999999999971,
                                         0.051000000000000156, 0.133, 0.133, 0.13200000000000012, 0.13200000000000012, 0.2749999999999999, 0.274, 0.2709999999999999, 0.5739999999999998, 0.5629999999999997,
                                         0.556, 0.5449999999999999, 0.5380000000000003, 3.3000000000000007, 3.1099999999999994, 6.530000000000001, 13.009999999999998 };
    Double_t p9093_d5x1y1_yval[] = { 1.502, 0.976, 0.6514, 0.4465, 0.3123, 0.2231, 0.162, 0.1194, 0.08913,
                                     0.06728, 0.05137, 0.03968, 0.0308, 0.02412, 0.01897, 0.01509, 0.01204, 0.009758, 0.007881,
                                     0.006396, 0.004563, 0.002838, 0.001825, 0.001188, 6.746E-4, 3.28E-4, 1.701E-4, 7.501E-5, 2.682E-5,
                                     1.132E-5, 5.342E-6, 2.821E-6, 6.721E-7, 7.883E-8, 1.096E-8, 7.32E-10 };
    Double_t p9093_d5x1y1_yerrminus[] = { 0.025470436993502876, 0.01653043049923383, 0.011060362768752208, 0.007561310558421734, 0.005306266194086007, 0.0037212339879265856, 0.002704202827895866, 0.0019991760404976847, 0.0014921546134365566,
                                          0.0011251364141294155, 8.688200064455238E-4, 6.714072460139226E-4, 5.214962707440965E-4, 4.092867835442528E-4, 3.2247847834545484E-4, 2.5637170274622745E-4, 2.0486532144069672E-4, 1.6666000955238184E-4, 1.352552050606556E-4,
                                          1.102505733363777E-4, 7.741420412947484E-5, 4.852996007622508E-5, 3.148570443868137E-5, 2.0772167637490315E-5, 1.169934826218965E-5, 6.232708465025458E-6, 3.3812194309154206E-6, 1.6499859418795059E-6, 7.395709313514154E-7,
                                          3.749271162506121E-7, 2.1025143899626464E-7, 1.330288408579132E-7, 5.421076425397451E-8, 1.0801381786141993E-8, 2.121034360872072E-9, 2.2208906344077367E-10 };
    Double_t p9093_d5x1y1_yerrplus[] = { 0.025470436993502876, 0.01653043049923383, 0.011060362768752208, 0.007561310558421734, 0.005306266194086007, 0.0037212339879265856, 0.002704202827895866, 0.0019991760404976847, 0.0014921546134365566,
                                         0.0011251364141294155, 8.688200064455238E-4, 6.714072460139226E-4, 5.214962707440965E-4, 4.092867835442528E-4, 3.2247847834545484E-4, 2.5637170274622745E-4, 2.0486532144069672E-4, 1.6666000955238184E-4, 1.352552050606556E-4,
                                         1.102505733363777E-4, 7.741420412947484E-5, 4.852996007622508E-5, 3.148570443868137E-5, 2.0772167637490315E-5, 1.169934826218965E-5, 6.232708465025458E-6, 3.3812194309154206E-6, 1.6350133358477538E-6, 7.218146039669743E-7,
                                         3.5748229676446916E-7, 1.9441290492145835E-7, 1.261577682903435E-7, 5.0316249472710104E-8, 8.770346885386005E-9, 1.816219909592448E-9, 1.9310629741155518E-10 };
    Int_t p9093_d5x1y1_numpoints = 36;

    Double_t xSection = ReturnCorrectXSection( "8TeV", 1);
    for(Int_t i=0; i<p9093_d5x1y1_numpoints; i++){
      p9093_d5x1y1_yval[i]     *= xSection*recalcBarn;
      p9093_d5x1y1_yerrminus[i]*= xSection*recalcBarn;
      p9093_d5x1y1_yerrplus[i] *= xSection*recalcBarn;
    }

    TGraphAsymmErrors* p9093_d5x1y1 = new TGraphAsymmErrors(p9093_d5x1y1_numpoints, p9093_d5x1y1_xval, p9093_d5x1y1_yval, p9093_d5x1y1_xerrminus, p9093_d5x1y1_xerrplus, p9093_d5x1y1_yerrminus, p9093_d5x1y1_yerrplus);
    p9093_d5x1y1->SetName("/HepData/9093/d5x1y1");
    p9093_d5x1y1->SetTitle("/HepData/9093/d5x1y1");
    cout << "chargedHadrons:" << endl;
    p9093_d5x1y1->Print();

    TGraphAsymmErrors* graphRatioChargedHadronsInverse = (TGraphAsymmErrors*) p9093_d5x1y1->Clone();
    Double_t * yValueCharged     = graphRatioChargedHadronsInverse->GetY();
    Double_t*  yErrorLowCharged  = graphRatioChargedHadronsInverse->GetEYlow();
    Double_t*  yErrorHighCharged = graphRatioChargedHadronsInverse->GetEYhigh();
    Double_t * yValueChargedInv     = p9093_d5x1y1->GetY();
    Double_t*  yErrorLowChargedInv  = p9093_d5x1y1->GetEYlow();
    Double_t*  yErrorHighChargedInv = p9093_d5x1y1->GetEYhigh();
    for(Int_t i=0; i<p9093_d5x1y1->GetN(); i++){
      yValueCharged[i] = 1/yValueChargedInv[i];
      yErrorLowCharged[i] = yValueCharged[i]*(yErrorLowChargedInv[i]/yValueChargedInv[i]);
      yErrorHighCharged[i] = yValueCharged[i]*(yErrorHighChargedInv[i]/yValueChargedInv[i]);
    }
    TGraphAsymmErrors* graphRatioChargedHadrons = CalculateGraphErrMultiplicationOfFit(graphRatioChargedHadronsInverse, fitTCMInvXSectionPi0Plot);

    // **********************************************************************************************************************
    // **************************Fit ATLAS charged particles and plot ratio to fit ******************************************
    // **********************************************************************************************************************

    Double_t paramTCMATLASChargedNew[5]  = { p9093_d5x1y1->GetY()[0],0.1,
                                             p9093_d5x1y1->GetY()[18],0.6,3.0};
    TF1* fitATLASTCMInvXSectionCharged    = FitObject("tcm","fitATLASTCMInvXSectionCharged","Pi0",p9093_d5x1y1,0.5,40. ,paramTCMATLASChargedNew,"QNRMEX0+","", kFALSE);
    //TF1* fitATLASTCMInvXSectionCharged   = FitObject("m","fitATLASTCMInvXSectionCharged","Pi0",p9093_d5x1y1,5,30. ,NULL,"QNRMEX0+","", kFALSE);
    //Double_t paramGraph[3]                              = {5e9, 6., 0.13};
    //TF1* fitATLASTCMInvXSectionCharged                   = FitObject("l","fitATLASTCMInvXSectionCharged","Pi0",p9093_d5x1y1,2,35.,paramGraph,"QNRMEX0+");
    cout << WriteParameterToFile(fitATLASTCMInvXSectionCharged) << endl;

    TCanvas* canvasDummyATLAS       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasDummyATLAS,  0.15, 0.01, 0.015, 0.08);
    canvasDummyATLAS->SetLogy();
    canvasDummyATLAS->SetLogx();
    TH2F* histo2DDummyATLAS;
    histo2DDummyATLAS               = new TH2F("histo2DDummyATLAS","histo2DDummyATLAS",1000,0.3,110.,1000,1,5e11);
    SetStyleHistoTH2ForGraphs(histo2DDummyATLAS, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 0.8,1.55);
    histo2DDummyATLAS->DrawCopy();

    DrawGammaSetMarkerTGraphAsym(p9093_d5x1y1, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
    p9093_d5x1y1->Draw("pEsame");

//    fitInvXSectionEta->SetLineColor(kBlue+2);
//    fitInvXSectionEta->Draw("same");
    fitATLASTCMInvXSectionCharged->SetLineColor(kRed+2);
    fitATLASTCMInvXSectionCharged->Draw("same");

    TLatex *labelTCCMATLAS1= new TLatex(0.48, 0.94, Form("TCM low:"));
    TLatex *labelTCCMATLAS2= new TLatex(0.48, 0.90, Form("A_{1}: (%.1e #pm %.1e) - T_{e}: (%.3f #pm %.3f)",fitATLASTCMInvXSectionCharged->GetParameter(0),fitATLASTCMInvXSectionCharged->GetParError(0),fitATLASTCMInvXSectionCharged->GetParameter(1),fitATLASTCMInvXSectionCharged->GetParError(1)));
    TLatex *labelTCCMATLAS3= new TLatex(0.48, 0.86, Form("TCM high:"));
    TLatex *labelTCCMATLAS4= new TLatex(0.48, 0.82, Form("A_{2}: (%.1e #pm %.1e) - T: (%.3f #pm %.3f) - n: (%.3f #pm %.3f)",fitATLASTCMInvXSectionCharged->GetParameter(2),fitATLASTCMInvXSectionCharged->GetParError(2),fitATLASTCMInvXSectionCharged->GetParameter(3),fitATLASTCMInvXSectionCharged->GetParError(3),fitATLASTCMInvXSectionCharged->GetParameter(4),fitATLASTCMInvXSectionCharged->GetParError(4)));

    TLatex *labelTCCMATLAS5= new TLatex(0.55, 0.75, Form("Bylinkin-Rostovtsev:"));
    TLatex *labelTCCMATLAS6= new TLatex(0.55, 0.71, Form("#it{A}_{1} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}_{2}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{n}"));

    SetStyleTLatex( labelTCCMATLAS1, 0.03,4);
    labelTCCMATLAS1->Draw();
    SetStyleTLatex( labelTCCMATLAS2, 0.02,4);
    labelTCCMATLAS2->Draw();
    SetStyleTLatex( labelTCCMATLAS3, 0.03,4);
    labelTCCMATLAS3->Draw();
    SetStyleTLatex( labelTCCMATLAS4, 0.02,4);
    labelTCCMATLAS4->Draw();
    SetStyleTLatex( labelTCCMATLAS5, 0.03,4);
    labelTCCMATLAS5->Draw();
    SetStyleTLatex( labelTCCMATLAS6, 0.03,4);
    labelTCCMATLAS6->Draw();

    canvasDummyATLAS->Update();
    canvasDummyATLAS->Print(Form("%s/ATLAS_ComparisonWithFit_8TeV.%s",outputDir.Data(),suffix.Data()));

//-------------------------------------------------------------------------------------------------

    TGraphAsymmErrors* graphRatioATLASChargedFit     = (TGraphAsymmErrors*)p9093_d5x1y1->Clone();
    graphRatioATLASChargedFit                        = CalculateGraphErrRatioToFit(graphRatioATLASChargedFit, fitATLASTCMInvXSectionCharged);

    canvasRatioToCombFit->cd();
    histo2DPi0RatioToCombFit->GetXaxis()->SetRangeUser(0.5,80);
    histo2DPi0RatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioATLASChargedFit, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioATLASChargedFit, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);

        graphRatioATLASChargedFit->Draw("E2same");
        graphRatioATLASChargedFit->Draw("p,same,z");

        DrawGammaLines(0.5, 110. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.5, 110. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.5, 110. , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy->Draw();
        TLatex *labelATLAScharged      = new TLatex(0.82, 0.86, "ATLAS h^{#pm}");
        SetStyleTLatex( labelATLAScharged, textSizeLabelsPixel,4);
        labelATLAScharged->SetTextFont(43);
        labelATLAScharged->Draw();

        TLegend* legendATLAScharged        = GetAndSetLegend2(0.2, 0.92-(0.05*2), 0.45, 0.92, 38);
        legendATLAScharged->AddEntry(graphRatioATLASChargedFit,"ATLAS h^{#pm} / TCM fit");
        legendATLAScharged->Draw("");

    canvasRatioToCombFit->SaveAs(Form("%s/ATLAS_charged_RatioToFit.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ************************** plot ALICE pi0 to ATLAS charged particle comparison ***************************************
    // **********************************************************************************************************************

    textSizeLabelsPixel             = 48;
    TCanvas* canvasCompYieldPPInd   = new TCanvas("canvasCompYieldPPInd","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasCompYieldPPInd,   0.12, 0.01, 0.01, 0.11);
    canvasCompYieldPPInd->SetLogx();

    TH2F * histo2DCompCombinedRatio2;
    histo2DCompCombinedRatio2       = new TH2F("histo2DCompCombinedRatio2","histo2DCompCombinedRatio2",1000,0.23,110.,1000,0.1,4.    );
    SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio2, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
    histo2DCompCombinedRatio2->GetYaxis()->SetNoExponent(kTRUE);
//  histo2DCompCombinedRatio2->GetXaxis()->SetLabelOffset(-0.01);
    histo2DCompCombinedRatio2->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DCompCombinedRatio2->GetXaxis()->SetNoExponent(kTRUE);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.1,1.45);
    histo2DCompCombinedRatio2->GetYaxis()->SetTitle("#pi^{0}/h^{#pm} (a.u.)");
    histo2DCompCombinedRatio2->DrawCopy();

    DrawGammaSetMarkerTGraphAsym(graphRatioChargedHadrons, markerStyleCombHighPt, markerSizeComparison, kBlack , kBlack);
    graphRatioChargedHadrons->Draw("E1psame");

    cout << "graphRatioChargedHadrons:" << endl;
    graphRatioChargedHadrons->Print();
    TLegend* legendPi0CompChargedPionsPP3   = GetAndSetLegend2(0.15, 0.84, 0.9, 0.92, 0.85* textSizeLabelsPixel);
    legendPi0CompChargedPionsPP3->SetNColumns(2);
    legendPi0CompChargedPionsPP3->SetMargin(0.12);
    legendPi0CompChargedPionsPP3->AddEntry(graphRatioChargedHadrons,"#pi^{0} (ALICE) / h^{#pm} (ATLAS)","p");
    legendPi0CompChargedPionsPP3->Draw();

    TLatex *labelRatioTheoryPPA   = new TLatex(0.15,0.8,"#pi^{0}: fit to ALICE combined result");
    SetStyleTLatex( labelRatioTheoryPPA, 0.85*textsizeLabelsPP,4);
    labelRatioTheoryPPA->Draw();
    TLatex *labelRatioTheoryPPAT   = new TLatex(0.15,0.75,"h^{#pm}: measured by ATLAS with N_{ch} >= 1, #it{p}_{T} > 500 MeV, |#eta| < 2.5");
    SetStyleTLatex( labelRatioTheoryPPAT, 0.85*textsizeLabelsPP,4);
    labelRatioTheoryPPAT->Draw();

    labelRatioTheoryPP->Draw();

    DrawGammaLines(0.23, 110 , 1, 1 ,1, kGray, 1);

    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedHadronToNeutralPions_PP8TeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
//---------------------------------------------------------------------------------------------------------------------------------------------------
    TGraphAsymmErrors* graphRatioATLASchargedToALICE     = (TGraphAsymmErrors*)graphCombPi0InvXSectionSysA->Clone();
    graphRatioATLASchargedToALICE                        = CalculateGraphErrRatioToFit(graphRatioATLASchargedToALICE, fitATLASTCMInvXSectionCharged);

    histo2DCompCombinedRatio2->DrawCopy();

    DrawGammaSetMarkerTGraphAsym(graphRatioATLASchargedToALICE, markerStyleCombHighPt, markerSizeComparison, kBlack , kBlack);
    graphRatioATLASchargedToALICE->Draw("E1psame");

    cout << "graphRatioChargedHadrons:" << endl;
    graphRatioATLASchargedToALICE->Print();
    TLegend* legendPi0CompChargedPionsPP4   = GetAndSetLegend2(0.15, 0.84, 0.9, 0.92, 0.85* textSizeLabelsPixel);
    legendPi0CompChargedPionsPP4->SetNColumns(2);
    legendPi0CompChargedPionsPP4->SetMargin(0.12);
    legendPi0CompChargedPionsPP4->AddEntry(graphRatioATLASchargedToALICE,"#pi^{0} (ALICE) / h^{#pm} (ATLAS)","p");
    legendPi0CompChargedPionsPP4->Draw();

    TLatex *labelRatioTheoryPPA2   = new TLatex(0.15,0.8,"#pi^{0}: ALICE combined result");
    SetStyleTLatex( labelRatioTheoryPPA2, 0.85*textsizeLabelsPP,4);
    labelRatioTheoryPPA2->Draw();
    TLatex *labelRatioTheoryPPAT2   = new TLatex(0.15,0.75,"h^{#pm}: fit to ATLAS result with N_{ch} >= 1, #it{p}_{T} > 500 MeV, |#eta| < 2.5");
    SetStyleTLatex( labelRatioTheoryPPAT2, 0.85*textsizeLabelsPP,4);
    labelRatioTheoryPPAT2->Draw();

    labelRatioTheoryPP->Draw();

    DrawGammaLines(0.23, 110 , 1, 1 ,1, kGray, 1);

    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedHadronToNeutralPions2_PP8TeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
//---------------------------------------------------------------------------------------------------------------------------------------------------
    TGraphAsymmErrors* graphRatioATLASchargedToALICEmerged     = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionStat->Clone();
    graphRatioATLASchargedToALICEmerged                        = CalculateGraphErrRatioToFit(graphRatioATLASchargedToALICEmerged, fitATLASTCMInvXSectionCharged);

    histo2DCompCombinedRatio2->DrawCopy();

    DrawGammaSetMarkerTGraphAsym(graphRatioATLASchargedToALICE, markerStyleCombHighPt, markerSizeComparison, kBlack , kBlack);
    graphRatioATLASchargedToALICE->Draw("E1psame");
    DrawGammaSetMarkerTGraphAsym(graphRatioATLASchargedToALICEmerged, markerStyleCombHighPt, markerSizeComparison, colorDet[9] , colorDet[9]);
    graphRatioATLASchargedToALICEmerged->Draw("E1psame");

    cout << "graphRatioChargedHadrons:" << endl;
    graphRatioATLASchargedToALICE->Print();
    TLegend* legendPi0CompChargedPionsPP5   = GetAndSetLegend2(0.15, 0.84, 0.9, 0.92, 0.85* textSizeLabelsPixel);
    legendPi0CompChargedPionsPP5->SetNColumns(2);
    legendPi0CompChargedPionsPP5->SetMargin(0.12);
    legendPi0CompChargedPionsPP5->AddEntry(graphRatioATLASchargedToALICE,"#pi^{0} (ALICE) / h^{#pm} (ATLAS)","p");
    legendPi0CompChargedPionsPP5->AddEntry(graphRatioATLASchargedToALICEmerged,"#pi^{0} (ALICE merged) / h^{#pm} (ATLAS)","p");
    legendPi0CompChargedPionsPP5->Draw();

    labelRatioTheoryPPA2->Draw();
    labelRatioTheoryPPAT2->Draw();

    labelRatioTheoryPP->Draw();

    DrawGammaLines(0.23, 110 , 1, 1 ,1, kGray, 1);

    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedHadronToNeutralPionsWithMerged_PP8TeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // **************************Plot example invariant mass bins ***********************************************************
    // **********************************************************************************************************************

    textSizeLabelsPixel                 = 100*3/5;
    TCanvas* canvasInvMassSamplePlot    = new TCanvas("canvasInvMassSamplePlot","",0,0,1500,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasInvMassSamplePlot,  0.09, 0.01, 0.035, 0.08);

    Style_t markerStyleInvMassSGBG      = 0;
    Size_t markerSizeInvMassSGBG        = 0;
    Color_t markerColorInvMassSGBG      = kBlack;
    Style_t markerStyleInvMassMBG       = 24;
    Size_t markerSizeInvMassMBG         = 1.5;
    Color_t markerColorInvMassMBG       = kGray+2;
    Style_t markerStyleInvMassBG        = 20;
    Size_t markerSizeInvMassBG          = 2;
    Color_t markerColorInvMassBG        = kBlack;
    Style_t markerStyleInvMassSG        = 20;
    Size_t markerSizeInvMassSG          = 3;
    Color_t markerColorInvMassSG        = kRed+2;
    Color_t fitColorInvMassSG           = kAzure+2;

    Double_t marginInvMass          = 0.1*1500;
    Double_t textsizeLabelsInvMass  = 0;
    Double_t textsizeFacInvMass     = 0;
    if (canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) < canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1())){
        textsizeLabelsInvMass       = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;
        textsizeFacInvMass          = (Double_t)1./canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;
    } else {
        textsizeLabelsInvMass       = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1());
        textsizeFacInvMass          = (Double_t)1./canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1());
    }
    cout << textsizeLabelsInvMass << endl;

    TH2F * histo2DPi0InvMassDummy;
    histo2DPi0InvMassDummy             = new TH2F("histo2DPi0InvMassDummy","histo2DPi0InvMassDummy",11000,0.02,0.255,21000,-1000,20000);
    SetStyleHistoTH2ForGraphs(histo2DPi0InvMassDummy, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","Counts",0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,
                            0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,0.88, 0.115/(textsizeFacInvMass*marginInvMass));

    TH2F * histo2DEtaInvMassDummy;
    histo2DEtaInvMassDummy             = new TH2F("histo2DEtaInvMassDummy","histo2DEtaInvMassDummy",11000,0.35,0.695,21000,-1000,20000);
    SetStyleHistoTH2ForGraphs(histo2DEtaInvMassDummy, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","Counts",0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,
                            0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,0.88, 0.115/(textsizeFacInvMass*marginInvMass));

    if (plotInvMassBins){
        for (Int_t i =0 ; i < 3; i++){
            if (haveAllPi0InvMassPCM[i]){
                canvasInvMassSamplePlot->cd();
                histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(0.02,0.255);
                histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSubPCM[i]->GetMinimum(),1.15*histoPi0InvMassSigPlusBGPCM[i]->GetMaximum());
                histo2DPi0InvMassDummy->DrawCopy();

                TLatex *labelInvMassPtRangePCM = new TLatex(0.945,0.9,"#pi^{0}: 0.5 GeV/#it{c} < #it{p}_{T} < 0.6 GeV/#it{c}");

                DrawGammaSetMarker(histoPi0InvMassSigPlusBGPCM[i], markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
                histoPi0InvMassSigPlusBGPCM[i]->SetLineWidth(1);
                histoPi0InvMassSigPlusBGPCM[i]->Draw("hist,e,same");
                DrawGammaSetMarker(histoPi0InvMassBGTotPCM[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
                histoPi0InvMassBGTotPCM[i]->Draw("same");

                DrawGammaSetMarker(histoPi0InvMassSigRemBGSubPCM[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
                histoPi0InvMassSigRemBGSubPCM[i]->Draw("same");
                fitPi0InvMassSigPCM[i]->SetNpx(1000);
                fitPi0InvMassSigPCM[i]->SetRange(0,0.255);
                fitPi0InvMassSigPCM[i]->SetLineColor(fitColorInvMassSG);
                fitPi0InvMassSigPCM[i]->Draw("same");

                //
                TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9,collisionSystem8TeV.Data());
                SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
                labelInvMassEnergy->SetTextFont(43);
                labelInvMassEnergy->Draw();

                TLatex *labelInvMassTrigger      = new TLatex(0.135,0.9-0.8*textsizeLabelsPP,Form("%s triggered",nameTrigger[i].Data()));
                SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4);
                labelInvMassTrigger->SetTextFont(43);
                labelInvMassTrigger->Draw();

                TLatex *labelInvMassRecoPCM  = new TLatex(0.135,0.9-2*0.8*textsizeLabelsPP,"PCM");
                SetStyleTLatex( labelInvMassRecoPCM, 0.85*textSizeLabelsPixel,4);
                labelInvMassRecoPCM->SetTextFont(43);
                labelInvMassRecoPCM->Draw();

                SetStyleTLatex( labelInvMassPtRangePCM, 0.85*textSizeLabelsPixel,4);
                labelInvMassPtRangePCM->SetTextAlign(31);
                labelInvMassPtRangePCM->SetTextFont(43);
                labelInvMassPtRangePCM->Draw();

                TLegend* legendInvMassPCM  = GetAndSetLegend2(0.67, 0.88-5*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
                legendInvMassPCM->SetMargin(0.25);
                legendInvMassPCM->AddEntry(histoPi0InvMassSigPlusBGPCM[i],"Raw real events","l");
                legendInvMassPCM->AddEntry(histoPi0InvMassBGTotPCM[i],"Mixed event +","p");
                legendInvMassPCM->AddEntry((TObject*)0,"corr. BG","");
                legendInvMassPCM->AddEntry(histoPi0InvMassSigRemBGSubPCM[i],"BG subtracted","p");
                legendInvMassPCM->AddEntry(fitPi0InvMassSigPCM[i], "Fit","l");
                legendInvMassPCM->Draw();
                canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBinPCM_%s.%s",outputDir.Data(), nameTrigger[i].Data(), suffix.Data()));
            } else {
                cout << "missing partial input for invariant mass bin for PCM for trigger: " << nameTrigger[i].Data() << endl;
            }

            if (haveAllPi0InvMassPCMEMCAL[i]){
                canvasInvMassSamplePlot->cd();
                histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(0.02,0.255);
                histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSubPCMEMCAL[i]->GetMinimum(),1.15*histoPi0InvMassSigPlusBGPCMEMCAL[i]->GetMaximum());
                histo2DPi0InvMassDummy->DrawCopy();

                TLatex *labelInvMassPtRangePCMEMCAL = new TLatex(0.945,0.9,Form("#pi^{0}: %s", histoPi0InvMassSigPCMEMCAL[i]->GetTitle()));

                DrawGammaSetMarker(histoPi0InvMassSigPlusBGPCMEMCAL[i], markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
                histoPi0InvMassSigPlusBGPCMEMCAL[i]->SetLineWidth(1);
                histoPi0InvMassSigPlusBGPCMEMCAL[i]->Draw("hist,e,same");
                DrawGammaSetMarker(histoPi0InvMassBGTotPCMEMCAL[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
                histoPi0InvMassBGTotPCMEMCAL[i]->Draw("same");

                DrawGammaSetMarker(histoPi0InvMassSigRemBGSubPCMEMCAL[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
                histoPi0InvMassSigRemBGSubPCMEMCAL[i]->Draw("same");
                fitPi0InvMassSigPCMEMCAL[i]->SetRange(0,0.255);
                fitPi0InvMassSigPCMEMCAL[i]->SetLineColor(fitColorInvMassSG);
                fitPi0InvMassSigPCMEMCAL[i]->Draw("same");

                //
                TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9,collisionSystem8TeV.Data());
                SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
                labelInvMassEnergy->SetTextFont(43);
                labelInvMassEnergy->Draw();

                TLatex *labelInvMassTrigger      = new TLatex(0.135,0.9-0.8*textsizeLabelsPP,Form("%s triggered",nameTrigger[i].Data()));
                SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4);
                labelInvMassTrigger->SetTextFont(43);
                labelInvMassTrigger->Draw();

                TLatex *labelInvMassRecoPCMEMC  = new TLatex(0.135,0.9-2*0.8*textsizeLabelsPP,"PCM-EMC");
                SetStyleTLatex( labelInvMassRecoPCMEMC, 0.85*textSizeLabelsPixel,4);
                labelInvMassRecoPCMEMC->SetTextFont(43);
                labelInvMassRecoPCMEMC->Draw();

                SetStyleTLatex( labelInvMassPtRangePCMEMCAL, 0.85*textSizeLabelsPixel,4);
                labelInvMassPtRangePCMEMCAL->SetTextAlign(31);
                labelInvMassPtRangePCMEMCAL->SetTextFont(43);
                labelInvMassPtRangePCMEMCAL->Draw();

                TLegend* legendInvMassPCMEMCAL  = GetAndSetLegend2(0.67, 0.88-5*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
                legendInvMassPCMEMCAL->SetMargin(0.25);
                legendInvMassPCMEMCAL->AddEntry(histoPi0InvMassSigPlusBGPCMEMCAL[i],"Raw real events","l");
                legendInvMassPCMEMCAL->AddEntry(histoPi0InvMassBGTotPCMEMCAL[i],"Mixed event +","p");
                legendInvMassPCMEMCAL->AddEntry((TObject*)0,"corr. BG","");
                legendInvMassPCMEMCAL->AddEntry(histoPi0InvMassSigRemBGSubPCMEMCAL[i],"BG subtracted","p");
                legendInvMassPCMEMCAL->AddEntry(fitPi0InvMassSigPCMEMCAL[i], "Fit","l");
                legendInvMassPCMEMCAL->Draw();
                canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBinPCMEMC_%s.%s",outputDir.Data(), nameTrigger[i].Data(), suffix.Data()));
            } else {
                cout << "missing partial input for invariant mass bin for PCM-EMC for trigger: " << nameTrigger[i].Data() << endl;
            }


            if (haveAllPi0InvMassEMCAL[i]){
                canvasInvMassSamplePlot->cd();
                histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(0.02,0.3);
                if(i==1) histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(0.07,0.3);
                else if(i==2) histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(0.08,0.3);
                histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSubEMCAL[i]->GetMinimum(),1.25*histoPi0InvMassSigPlusBGEMCAL[i]->GetMaximum());
                histo2DPi0InvMassDummy->DrawCopy();

                TLatex *labelInvMassPtRangeEMCAL = new TLatex(0.945,0.9,Form("#pi^{0}: %s", histoPi0InvMassSigEMCAL[i]->GetTitle()));

                DrawGammaSetMarker(histoPi0InvMassSigPlusBGEMCAL[i], markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
                histoPi0InvMassSigPlusBGEMCAL[i]->SetLineWidth(1);
                histoPi0InvMassSigPlusBGEMCAL[i]->Draw("hist,e,same");
                DrawGammaSetMarker(histoPi0InvMassBGTotEMCAL[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
                histoPi0InvMassBGTotEMCAL[i]->Draw("same");

                DrawGammaSetMarker(histoPi0InvMassSigRemBGSubEMCAL[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
                histoPi0InvMassSigRemBGSubEMCAL[i]->Draw("same");
                fitPi0InvMassSigEMCAL[i]->SetRange(0.,0.255);
                fitPi0InvMassSigEMCAL[i]->SetLineColor(fitColorInvMassSG);
                fitPi0InvMassSigEMCAL[i]->Draw("same");

                //
                TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9,collisionSystem8TeV.Data());
                SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
                labelInvMassEnergy->SetTextFont(43);
                labelInvMassEnergy->Draw();

                TLatex *labelInvMassTrigger      = new TLatex(0.135,0.9-0.8*textsizeLabelsPP,Form("%s triggered",nameTrigger[i].Data()));
                SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4);
                labelInvMassTrigger->SetTextFont(43);
                labelInvMassTrigger->Draw();

                TLatex *labelInvMassRecoEMC  = new TLatex(0.135,0.9-2*0.8*textsizeLabelsPP,"EMC");
                SetStyleTLatex( labelInvMassRecoEMC, 0.85*textSizeLabelsPixel,4);
                labelInvMassRecoEMC->SetTextFont(43);
                labelInvMassRecoEMC->Draw();

                SetStyleTLatex( labelInvMassPtRangeEMCAL, 0.85*textSizeLabelsPixel,4);
                labelInvMassPtRangeEMCAL->SetTextAlign(31);
                labelInvMassPtRangeEMCAL->SetTextFont(43);
                labelInvMassPtRangeEMCAL->Draw();

                TLegend* legendInvMassEMCAL  = GetAndSetLegend2(0.67, 0.88-5*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
                legendInvMassEMCAL->SetMargin(0.25);
                legendInvMassEMCAL->AddEntry(histoPi0InvMassSigPlusBGEMCAL[i],"Raw real events","l");
                legendInvMassEMCAL->AddEntry(histoPi0InvMassBGTotEMCAL[i],"Mixed event +","p");
                legendInvMassEMCAL->AddEntry((TObject*)0,"corr. BG","");
                legendInvMassEMCAL->AddEntry(histoPi0InvMassSigRemBGSubEMCAL[i],"BG subtracted","p");
                legendInvMassEMCAL->AddEntry(fitPi0InvMassSigEMCAL[i], "Fit","l");
                legendInvMassEMCAL->Draw();
                canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBinEMC_%s.%s",outputDir.Data(), nameTrigger[i].Data(), suffix.Data()));
            } else {
                cout << "missing partial input for invariant mass bin for PCM-EMC for trigger: " << nameTrigger[i].Data() << endl;
            }
        }
    }

      // **********************************************************************************************************************
      // **************************Plot example invariant mass bin PHOS low pt ***********************************************
      // **********************************************************************************************************************

      TFile* filePHOSMB                                     = new TFile(fileNamePHOSMB.Data());
      TDirectory* directoryPHOSPi0L                         = (TDirectory*)filePHOSMB->Get("Pi08TeV");
      TH1D* histoPHOSSignalPlusBGPi0                        = (TH1D*)directoryPHOSPi0L->Get("InvMassSigPlusBG_PtBin07");
      TH1D* histoPHOSTotalBGPi0                             = (TH1D*)directoryPHOSPi0L->Get("InvMassBG_PtBin07");
      TH1D* histoPHOSSignalPi0                              = (TH1D*)directoryPHOSPi0L->Get("InvMassSig_PtBin07");
      TF1* fitPHOSlow                                       = (TF1*)directoryPHOSPi0L->Get("InvMassFinalFit");

      canvasInvMassSamplePlot->cd();
      histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(0.9*histoPHOSSignalPi0->GetMinimum(),0.9*histoPHOSSignalPi0->GetMaximum());
      histo2DPi0InvMassDummy->DrawCopy();

      TLatex *labelInvMassPtRangePHOSl = new TLatex(0.945,0.9,"#pi^{0}: 1.40 GeV/c < p_{T} < 1.60 GeV/c");

      DrawGammaSetMarker(histoPHOSSignalPlusBGPi0, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
      histoPHOSSignalPlusBGPi0->SetLineWidth(1);
      histoPHOSSignalPlusBGPi0->Draw("hist,e,same");
      DrawGammaSetMarker(histoPHOSTotalBGPi0, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
      histoPHOSTotalBGPi0->Draw("same");

      DrawGammaSetMarker(histoPHOSSignalPi0, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
      histoPHOSSignalPi0->Draw("same");
      fitPHOSlow->SetRange(0,0.255);
      fitPHOSlow->SetLineColor(fitColorInvMassSG);
      fitPHOSlow->Draw("same");

      //
      TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9,collisionSystem8TeV.Data());
      SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
      labelInvMassEnergy->SetTextFont(43);
      labelInvMassEnergy->Draw();

      TLatex *labelInvMassTriggerPHOSl      = new TLatex(0.135,0.9-0.8*textsizeLabelsPP,"INT7 triggered");
      SetStyleTLatex( labelInvMassTriggerPHOSl, 0.85*textSizeLabelsPixel,4);
      labelInvMassTriggerPHOSl->SetTextFont(43);
      labelInvMassTriggerPHOSl->Draw();

      TLatex *labelInvMassRecoPHOS  = new TLatex(0.135,0.9-2*0.8*textsizeLabelsPP,"PHOS");
      SetStyleTLatex( labelInvMassRecoPHOS, 0.85*textSizeLabelsPixel,4);
      labelInvMassRecoPHOS->SetTextFont(43);
      labelInvMassRecoPHOS->Draw();

      SetStyleTLatex( labelInvMassPtRangePHOSl, 0.85*textSizeLabelsPixel,4);
      labelInvMassPtRangePHOSl->SetTextAlign(31);
      labelInvMassPtRangePHOSl->SetTextFont(43);
      labelInvMassPtRangePHOSl->Draw();

      TLegend* legendInvMassPHOSl  = GetAndSetLegend2(0.67, 0.88-5*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
      legendInvMassPHOSl->SetMargin(0.25);
      legendInvMassPHOSl->AddEntry(histoPHOSSignalPlusBGPi0,"Raw real events","l");
      legendInvMassPHOSl->AddEntry(histoPHOSTotalBGPi0,"Mixed event +","p");
      legendInvMassPHOSl->AddEntry((TObject*)0,"corr. BG","");
      legendInvMassPHOSl->AddEntry(histoPHOSSignalPi0,"BG subtracted","p");
      legendInvMassPHOSl->AddEntry(fitPHOSlow, "Fit","l");
      legendInvMassPHOSl->Draw();
      canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBinPHOSlow.%s",outputDir.Data(), suffix.Data()));

      // **********************************************************************************************************************
      // **************************Plot example invariant mass bin PHOS high pt ***********************************************
      // **********************************************************************************************************************

      TFile* filePHOSPHOS                                     = new TFile(fileNamePHOSPHOS.Data());
      TDirectory* directoryPHOSPi0H                           = (TDirectory*)filePHOSPHOS->Get("Pi08TeV");
      TH1D* histoPHOSHighSignalPlusBGPi0                      = (TH1D*)directoryPHOSPi0H->Get("InvMassSigPlusBG_PtBin07");
      TH1D* histoPHOSHighTotalBGPi0                           = (TH1D*)directoryPHOSPi0H->Get("InvMassBG_PtBin07");
      TH1D* histoPHOSHighSignalPi0                            = (TH1D*)directoryPHOSPi0H->Get("InvMassSig_PtBin07");
      TF1* fitPHOShigh                                        = (TF1*)directoryPHOSPi0H->Get("InvMassFinalFit");

      canvasInvMassSamplePlot->cd();
      histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(0.9*histoPHOSHighSignalPi0->GetMinimum(),1.6*histoPHOSHighSignalPi0->GetMaximum());
      histo2DPi0InvMassDummy->DrawCopy();

      TLatex *labelInvMassPtRangePHOSh = new TLatex(0.945,0.9,"#pi^{0}: 20.00 GeV/c < p_{T} < 22.00 GeV/c");

      DrawGammaSetMarker(histoPHOSHighSignalPlusBGPi0, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
      histoPHOSHighSignalPlusBGPi0->SetLineWidth(1);
      histoPHOSHighSignalPlusBGPi0->Draw("hist,e,same");
      DrawGammaSetMarker(histoPHOSHighTotalBGPi0, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
      histoPHOSHighTotalBGPi0->Draw("same");

      DrawGammaSetMarker(histoPHOSHighSignalPi0, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
      histoPHOSHighSignalPi0->Draw("same");
      fitPHOShigh->SetRange(0,0.255);
      fitPHOShigh->SetLineColor(fitColorInvMassSG);
      fitPHOShigh->Draw("same");

      //
      labelInvMassEnergy->Draw();

      TLatex *labelInvMassTriggerPHOSh      = new TLatex(0.135,0.9-0.8*textsizeLabelsPP,"PHOS triggered");
      SetStyleTLatex( labelInvMassTriggerPHOSh, 0.85*textSizeLabelsPixel,4);
      labelInvMassTriggerPHOSh->SetTextFont(43);
      labelInvMassTriggerPHOSh->Draw();

      labelInvMassRecoPHOS->Draw();

      SetStyleTLatex( labelInvMassPtRangePHOSh, 0.85*textSizeLabelsPixel,4);
      labelInvMassPtRangePHOSh->SetTextAlign(31);
      labelInvMassPtRangePHOSh->SetTextFont(43);
      labelInvMassPtRangePHOSh->Draw();

      TLegend* legendInvMassPHOSh  = GetAndSetLegend2(0.67, 0.88-5*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
      legendInvMassPHOSh->SetMargin(0.25);
      legendInvMassPHOSh->AddEntry(histoPHOSHighSignalPlusBGPi0,"Raw real events","l");
      legendInvMassPHOSh->AddEntry(histoPHOSHighTotalBGPi0,"Mixed event +","p");
      legendInvMassPHOSh->AddEntry((TObject*)0,"corr. BG","");
      legendInvMassPHOSh->AddEntry(histoPHOSHighSignalPi0,"BG subtracted","p");
      legendInvMassPHOSh->AddEntry(fitPHOShigh, "Fit","l");
      legendInvMassPHOSh->Draw();
      canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBinPHOShigh.%s",outputDir.Data(), suffix.Data()));

 // **********************************************************************************************************************
 // ************************* Saving of final results ********************************************************************
 // **********************************************************************************************************************

    TString nameOutputCommonFile    = Form("CombinedResultsPaperPP8TeVmerged_%s.root", dateForOutput.Data());
    TFile fCombResults(nameOutputCommonFile.Data(), "RECREATE");

    fCombResults.mkdir("Pi08TeV");
    TDirectoryFile* directoryPi0 = (TDirectoryFile*)fCombResults.Get("Pi08TeV");
    fCombResults.cd("Pi08TeV");
        // PCM component
        graphPCMPi0InvXSectionStat->Write("graphInvCrossSectionPi0PCM8TeVStatErr");
        graphPCMPi0InvXSectionSys->Write("graphInvCrossSectionPi0PCM8TeVSysErr");
        // PHOS component
        graphPHOSPi0InvXSectionStat->Write("graphInvCrossSectionPi0PHOS8TeVStatErr");
        graphPHOSPi0InvXSectionSys->Write("graphInvCrossSectionPi0PHOS8TeVSysErr");
        // EMCAL component
        graphEMCALPi0InvXSectionStat->Write("graphInvCrossSectionPi0EMCAL8TeVStatErr");
        graphEMCALPi0InvXSectionSys->Write("graphInvCrossSectionPi0EMCAL8TeVSysErr");
        // PCM-EMCal component
        graphPCMEMCALPi0InvXSectionStat->Write("graphInvCrossSectionPi0PCMEMCAL8TeVStatErr");
        graphPCMEMCALPi0InvXSectionSys->Write("graphInvCrossSectionPi0PCMEMCAL8TeVSysErr");
        // EMCal merged component
        graphEMCALMergedPi0InvXSectionStat->Write("graphInvCrossSectionPi0EMCALMerged8TeVStatErr");
        graphEMCALMergedPi0InvXSectionSys->Write("graphInvCrossSectionPi0EMCALMerged8TeVSysErr");
        // Final spectrum correlations Method A
        graphCombPi0InvXSectionTotA->Write("graphInvCrossSectionPi0Comb8TeVA");
        graphCombPi0InvXSectionStatA->Write("graphInvCrossSectionPi0Comb8TeVAStatErr");
        graphCombPi0InvXSectionSysA->Write("graphInvCrossSectionPi0Comb8TeVASysErr");

         // fits for eta
        fitInvXSectionPi0->Write("TsallisFitPi0");
        fitTCMInvXSectionPi0Plot->Write("TwoComponentModelFitPi0");

        if (bWCorrection.Contains("Y")){
            if(graphPCMPi0InvXSectionStat_yShifted)graphPCMPi0InvXSectionStat_yShifted->Write("graphInvCrossSectionPi0PCM8TeVStatErr_yShifted");
            if(graphPCMPi0InvXSectionSys_yShifted)graphPCMPi0InvXSectionSys_yShifted->Write("graphInvCrossSectionPi0PCM8TeVSysErr_yShifted");
        // PHOS component
            if(graphPHOSPi0InvXSectionStat_yShifted) graphPHOSPi0InvXSectionStat_yShifted->Write("graphInvCrossSectionPi0PHOS8TeVStatErr_yShifted");
            if(graphPHOSPi0InvXSectionSys_yShifted) graphPHOSPi0InvXSectionSys_yShifted->Write("graphInvCrossSectionPi0PHOS8TeVSysErr_yShifted");
            // EMCAL component
            if(graphEMCALPi0InvXSectionStat_yShifted)graphEMCALPi0InvXSectionStat_yShifted->Write("graphInvCrossSectionPi0EMCAL8TeVStatErr_yShifted");
            if(graphEMCALPi0InvXSectionSys_yShifted)graphEMCALPi0InvXSectionSys_yShifted->Write("graphInvCrossSectionPi0EMCAL8TeVSysErr_yShifted");
            // PCM-EMCal component
            if(graphPCMEMCALPi0InvXSectionStat_yShifted)graphPCMEMCALPi0InvXSectionStat_yShifted->Write("graphInvCrossSectionPi0PCMEMCAL8TeVStatErr_yShifted");
            if(graphPCMEMCALPi0InvXSectionSys_yShifted)graphPCMEMCALPi0InvXSectionSys_yShifted->Write("graphInvCrossSectionPi0PCMEMCAL8TeVSysErr_yShifted");
            // EMCal merged component
            if (graphEMCALMergedPi0InvXSectionStat_yShifted) graphEMCALMergedPi0InvXSectionStat_yShifted->Write("graphInvCrossSectionPi0EMCALMerged8TeVStatErr_yShifted");
            if (graphEMCALMergedPi0InvXSectionSys_yShifted) graphEMCALMergedPi0InvXSectionSys_yShifted->Write("graphInvCrossSectionPi0EMCALMerged8TeVSysErr_yShifted");
            // Final spectrum correlations Method A
            if(graphCombPi0InvXSectionTotA_yShifted)graphCombPi0InvXSectionTotA_yShifted->Write("graphInvCrossSectionPi0Comb8TeVA_yShifted");
            if(graphCombPi0InvXSectionStatA_yShifted)graphCombPi0InvXSectionStatA_yShifted->Write("graphInvCrossSectionPi0Comb8TeVAStatErr_yShifted");
            if(graphCombPi0InvXSectionSysA_yShifted)graphCombPi0InvXSectionSysA_yShifted->Write("graphInvCrossSectionPi0Comb8TeVASysErr_yShifted");

        }
        directoryPi0->mkdir("Supporting");
        directoryPi0->cd("Supporting");
            // Writing full correction factors
            histoPCMPi0AccTimesEff->Write("Pi0CorrectionFactorPCM");
            histoPHOSPi0AccTimesEff->Write("Pi0CorrectionFactorPHOS");
            graphEMCALPi0AccTimesEff->Write("Pi0CorrectionFactorEMCAL");
            graphPCMEMCALPi0AccTimesEff->Write("Pi0CorrectionFactorPCMEMCAL");
            graphEMCALMergedPi0AccTimesEff->Write("Pi0CorrectionFactorEMCALMerged");

            histoPCMPi0Mass->Write("Pi0MassDataPCM");
            histoPCMPi0TrueMass->Write("Pi0MassMCPCM");
            histoPHOSPi0Mass->Write("Pi0MassDataPHOS");
            histoPHOSPi0TrueMass->Write("Pi0MassMCPHOS");
            graphEMCALPi0Mass->Write("Pi0MassDataEMCAL");
            graphEMCALPi0MassMC->Write("Pi0MassMCEMCAL");
            graphPCMEMCALPi0Mass->Write("Pi0MassDataPCMEMCAL");
            graphPCMEMCALPi0MassMC->Write("Pi0MassMCPCMEMCAL");

            histoPCMPi0FWHMMeV->Write("Pi0WidthDataPCM");
            histoPCMPi0TrueFWHMMeV->Write("Pi0WidthMCPCM");
            histoPHOSPi0FWHMMeV->Write("Pi0WidthDataPHOS");
            histoPHOSPi0TrueFWHMMeV->Write("Pi0WidthMCPHOS");
            graphEMCALPi0FWHM->Write("Pi0WidthDataEMCAL");
            graphEMCALPi0FWHMMC->Write("Pi0WidthMCEMCAL");
            graphPCMEMCALPi0FWHM->Write("Pi0WidthDataPCMEMCAL");
            graphPCMEMCALPi0FWHMMC->Write("Pi0WidthMCPCMEMCAL");

    fCombResults.Close();


 // **********************************************************************************************************************
 // ************************* Saving only fits to final results **********************************************************
 // **********************************************************************************************************************

    TString nameOutputCommonFileFitsOnly    = Form("FitsPaperPP8TeVmerged_%s.root", dateForOutput.Data());
    TFile fFitsResults(nameOutputCommonFileFitsOnly.Data(), "RECREATE");

         // fits for pi0
        fitInvXSectionPi0->Write("TsallisFitPi0");
        fitTCMInvXSectionPi0Plot->Write("TwoComponentModelFitPi0");


    fFitsResults.Close();

}
