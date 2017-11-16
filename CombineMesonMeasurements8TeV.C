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
void CombineMesonMeasurements8TeV(      TString fileNamePCM         = "",
                                        TString fileNamePCMEMCAL    = "",
                                        TString fileNameEMCALLow    = "",
                                        TString fileNameEMCALmerged = "",
                                        TString fileNamePHOS        = "",
                                        TString suffix              = "eps",
                                        TString isMC                = "",
                                        TString thesisPlots         = "",
                                        TString bWCorrection        = "X",
                                        TString fileInputCorrFactors= "",
                                        Bool_t plotInvMassBins = kTRUE,
                                        Bool_t plotDate = kFALSE
                                    ){

    TString date = ReturnDateString(kTRUE);

    TString ALICEperfor = "ALICE performance";

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
    TString fileNameEtaToPi0                    = "ExternalInput/WorldDataPi0Eta.root";
    TString fileNamePOWHEG8TeV                  = "ExternalInput/Theory/POWHEG_Pythia/powhegDijetShoweredWithPythia_pi0_eta_invariantXsec_NNPDF23LOas0130.root";

    TString fileName2760GeV                     = "ExternalInput/CombinedResultsPaperPP2760GeV_2016_12_15_FrediV2Clusterizer.root";
    TString fileName7TeV                        = "ExternalInput/CombinedResultsPaperPP7TeV_2017_11_10.root";
    TString fileName7TeVpub                     = "ExternalInput/CombNeutralMesons/CombinedResultsPP_ShiftedX_PaperRAA_16_May_2014_including7TeVand900GeVpublished.root";

    TString fileNamePCMMB                       = "ExternalInput/PCM/8TeV/8TeV_data_PCMResults_InvMassBins.root";

    TString fileNamePHOSMB                      = "ExternalInput/PHOS/8TeV/data_PHOS-MBResultsFullCorrection_PP_NoBinShifting.root";
    TString fileNamePHOSPHOS                    = "ExternalInput/PHOS/8TeV/data_PHOS-PHOSResultsFullCorrection_PP_NoBinShifting.root";
    TString fileNamePHOSprelim                  = "ExternalInput/PHOS/8TeV/preliminary/pi0_8TeV_PHOS.root";
    TString fileNameChargedPionPP               = "ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_2016_08_14.root";
    //TString fileNameChargedHadronPP             = "ExternalInput/UnidentifiedCharged/ChargedHadrinSpectraPP_20_May_2015.root";
    TString outputDir                           = Form("%s/%s/CombineMesonMeasurements8TeV%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
    if(plotDate) outputDir.Append(dateForOutput);

    cout << outputDir.Data() << endl;
    cout << fileNamePCM.Data() << endl;

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMEMCAL.root", fileNamePCMEMCAL.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCALLow.root", fileNameEMCALLow.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCALmerged.root", fileNameEMCALmerged.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/Theory.root", fileNameTheory.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/ChargedPionsPP.root", fileNameChargedPionPP.Data(), outputDir.Data()));
    //gSystem->Exec(Form("cp %s %s/ChargedHadronsPP.root", fileNameChargedHadronPP.Data(), outputDir.Data()));

    fstream fLog;
    fLog.open(Form("%s/CombineMeson8TeV.log",outputDir.Data()), ios::out);
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << dateForOutput.Data() << endl;

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


    TString  nameMeasGlobal[11]                 = { "PCM", "PHOS", "EMC", "PCM-PHOS", "PCM-EMC",
                                                    "PCM-Dalitz", "PHOS-Dalitz", "EMC-Dalitz", "EMC high pT", "mEMC",
                                                    "PCMOtherDataset"};
    TString  nameTrigger[3]                     = {"INT7", "EMC7", "EGA"};
    TString  nameTriggerAlternative[3]          = {"MB trigger", "EMC-L0 trigger", "EMC-L1 trigger"};
    TString  nameSecPi0SourceRead[4]            = {"K0S", "K0L", "Lambda", "Rest"};
    TString  nameSecPi0SourceLabel[4]           = {"K^{0}_{s}", "K^{0}_{l}", "#Lambda", "had. int."};
    Double_t maxSecCorr[4]                      = { 0.05, 0.009, 0.00028, 0.04};

    Color_t  colorDet[11];
    Color_t  colorDetMC[11];
    Style_t  markerStyleDet[11];
    Style_t  markerStyleDetMC[11];
    Size_t   markerSizeDet[11];
    Size_t   markerSizeDetMC[11];

    Color_t  colorCGC                           = kCyan-8;
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

    Bool_t haveEffSecCorr[4][11]                    = { {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE},
                                                        {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE},
                                                        {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE},
                                                        {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE} };
    TGraphAsymmErrors* graphPi0EffSecCorrFromX[4][11]   = { { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                            { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                            { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                            { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL} };

    //************************** Read data for PCM **************************************************
    TFile* filePCM                                          = new TFile(fileNamePCM.Data());
    TH1D* histoPCMNumberOfEvents                            = (TH1D*)filePCM->Get("NEvents_MB");
    TDirectory* directoryPCMPi0                             = (TDirectory*)filePCM->Get("Pi08TeV");
        TH1D* histoPCMPi0Mass                               = (TH1D*)directoryPCMPi0->Get("Pi0_Mass_data_MB");
        histoPCMPi0Mass                                     ->Scale(1000.);
        TH1D* histoPCMPi0FWHMMeV                            = (TH1D*)directoryPCMPi0->Get("Pi0_Width_data_MB");
        histoPCMPi0FWHMMeV                                  ->Scale(1000.);
        TH1D* histoPCMPi0TrueMass                           = (TH1D*)directoryPCMPi0->Get("Pi0_Mass_MC_MB");
        histoPCMPi0TrueMass                                 ->Scale(1000.);
        TH1D* histoPCMPi0TrueFWHMMeV                        = (TH1D*)directoryPCMPi0->Get("Pi0_Width_MC_MB");
        histoPCMPi0TrueFWHMMeV                              ->Scale(1000.);
        TH1D* histoPCMPi0Acc                                = (TH1D*)directoryPCMPi0->Get("AcceptancePi0_MB");
        TH1D* histoPCMPi0TrueEffPt                          = (TH1D*)directoryPCMPi0->Get("EfficiencyPi0_MB");

        TGraphAsymmErrors* graphPCMPi0AccTimesEff      = (TGraphAsymmErrors*)directoryPCMPi0->Get("EffTimesAccPi0");
        for (Int_t k = 0; k < 4; k++){
            graphPi0EffSecCorrFromX[k][0]           = (TGraphAsymmErrors*)directoryPCMPi0->Get(Form("EffectiveSecondaryPi0CorrFrom%s",nameSecPi0SourceRead[k].Data()));
            if (graphPi0EffSecCorrFromX[k][0]){
                cout << nameSecPi0SourceRead[k].Data() << endl;
                graphPi0EffSecCorrFromX[k][0]->Print();
                Int_t nAboveZero                    = 0;
                for (Int_t i = 0; i< graphPi0EffSecCorrFromX[k][0]->GetN(); i++){
                    if(graphPi0EffSecCorrFromX[k][0]->GetY()[i] > 0) nAboveZero++;
                }
                if (nAboveZero>0){
                    haveEffSecCorr[k][0]            = kTRUE;
                } else {
                    graphPi0EffSecCorrFromX[k][0]   = NULL;
                }
            }
        }
        TGraphAsymmErrors* graphPCMPi0InvXSectionStat  = (TGraphAsymmErrors*)directoryPCMPi0->Get("graphInvCrossSectionPi0");
        TH1D* histoPCMPi0InvXSectionStat               = (TH1D*)directoryPCMPi0->Get("InvCrossSectionPi0");
        cout << "Pi0 stat PCM" << endl;
        graphPCMPi0InvXSectionStat->Print();
        TGraphAsymmErrors* graphPCMPi0InvXSectionSys   = (TGraphAsymmErrors*)directoryPCMPi0->Get("InvCrossSectionPi0Sys");
        cout << "Pi0 sys PCM" << endl;
        graphPCMPi0InvXSectionSys->Print();

//        TH1D* histoPCMPi0InvXSectionStat                    = (TH1D*)directoryPCMPi0->Get("InvCrossSectionPi0");
//        histoPCMPi0InvXSectionStat->SetBinContent(1,0);
//        histoPCMPi0InvXSectionStat->SetBinError(1,1);
//        TGraphAsymmErrors* graphPCMPi0InvXSectionStat       = new TGraphAsymmErrors(histoPCMPi0InvXSectionStat);
//        //graphPCMPi0InvXSectionStat->RemovePoint(graphPCMPi0InvXSectionStat->GetN()-1);
//        graphPCMPi0InvXSectionStat->RemovePoint(0);
//        cout << "Pi0 stat PCM" << endl;
//        graphPCMPi0InvXSectionStat->Print();
//        TGraphAsymmErrors* graphPCMPi0InvXSectionSysA       = (TGraphAsymmErrors*)directoryPCMPi0->Get("InvCrossSectionPi0SysA");
//        TGraphAsymmErrors* graphPCMPi0InvXSectionSys        = (TGraphAsymmErrors*)directoryPCMPi0->Get("InvCrossSectionPi0Sys");
//        //graphPCMPi0InvXSectionSys->RemovePoint(graphPCMPi0InvXSectionSys->GetN()-1);
//        graphPCMPi0InvXSectionSys->RemovePoint(graphPCMPi0InvXSectionSys->GetN()-1);
//        cout << "Pi0 sys PCM" << endl;
//        graphPCMPi0InvXSectionSys->Print();
////      TGraphAsymmErrors* graphPCMPi0CorrYieldSysErr    = (TGraphAsymmErrors*)directoryPCMPi0->Get("Pi0SystError");

        TH1D* histoPCMPi0AccTimesEff                        = (TH1D*)histoPCMPi0TrueEffPt->Clone("histoPCMPi0AccTimesEff");
        histoPCMPi0AccTimesEff->Multiply(histoPCMPi0Acc);
        // normalize to full acceptance (delta y and phi)
        histoPCMPi0AccTimesEff->Scale(2*TMath::Pi()*1.6);

        TH1D* histoPi0InvMassSigPlusBGPCM[3];
        TH1D* histoPi0InvMassSigPCM[3];
        TH1D* histoPi0InvMassSigRemBGSubPCM[3];
        TH1D* histoPi0InvMassBGPCM[3];
        TH1D* histoPi0InvMassRemBGPCM[3];
        TH1D* histoPi0InvMassBGMixedPCM[3];
        TF1* fitPi0InvMassSigPCM[3];
        TF1* fitPi0InvMassBGPCM[3];
        Bool_t haveAllPi0InvMassPCM[3]                 = {kFALSE, kFALSE, kFALSE};
        if (plotInvMassBins){
            for (Int_t i = 0; i < 1; i++){
                histoPi0InvMassSigPCM[i]               = (TH1D*)directoryPCMPi0->Get("Pi0_InvMassSig_Example_MB");
                histoPi0InvMassSigPlusBGPCM[i]         = (TH1D*)directoryPCMPi0->Get("Pi0_InvMassSigPlusBG_Example_MB");
                histoPi0InvMassBGPCM[i]                = (TH1D*)directoryPCMPi0->Get("Pi0_InvMassBG_Example_MB");
                fitPi0InvMassSigPCM[i]                 = (TF1*)directoryPCMPi0->Get("Pi0_InvMassSigFit_Example_MB");
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
                    histoPi0InvMassBGMixedPCM[i]         = (TH1D*)histoPi0InvMassBGPCM[i]->Clone(Form("Pi0_InvMassTotBG_Example_%s",nameTrigger[i].Data()));
                    //histoPi0InvMassBGMixedPCM[i]->Sumw2();
                    //histoPi0InvMassBGMixedPCM[i]->Add(histoPi0InvMassRemBGPCM[i]);
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

    TDirectory* directoryPCMEta                             = (TDirectory*)filePCM->Get("Eta8TeV");
        TH1D* histoPCMEtaMass                               = (TH1D*)directoryPCMEta->Get("Eta_Mass_data_MB");
        histoPCMEtaMass                                     ->Scale(1000.);
        TH1D* histoPCMEtaFWHMMeV                            = (TH1D*)directoryPCMEta->Get("Eta_Width_data_MB");
        histoPCMEtaFWHMMeV                                  ->Scale(1000.);
        TH1D* histoPCMEtaTrueMass                           = (TH1D*)directoryPCMEta->Get("Eta_Mass_MC_MB");
        histoPCMEtaTrueMass                                 ->Scale(1000.);
        TH1D* histoPCMEtaTrueFWHMMeV                        = (TH1D*)directoryPCMEta->Get("Eta_Width_MC_MB");
        histoPCMEtaTrueFWHMMeV                              ->Scale(1000.);
        TH1D* histoPCMEtaAcc                                = (TH1D*)directoryPCMEta->Get("AcceptanceEta_MB");
        TH1D* histoPCMEtaTrueEffPt                          = (TH1D*)directoryPCMEta->Get("EfficiencyEta_MB");

        TGraphAsymmErrors* graphPCMEtaAccTimesEff      = (TGraphAsymmErrors*)directoryPCMEta->Get("EffTimesAccEta");
        TH1D* histoPCMEtaInvXSectionStat               = (TH1D*)directoryPCMEta->Get("InvCrossSectionEta");
        TGraphAsymmErrors* graphPCMEtaInvXSectionStat  = (TGraphAsymmErrors*)directoryPCMEta->Get("graphInvCrossSectionEta");
        cout << "Eta stat PCM" << endl;
        graphPCMEtaInvXSectionStat->Print();
        TGraphAsymmErrors* graphPCMEtaInvXSectionSys   = (TGraphAsymmErrors*)directoryPCMEta->Get("InvCrossSectionEtaSys");
        cout << "Eta sys PCM" << endl;
        graphPCMEtaInvXSectionSys->Print();
        TH1D* histoPCMEtaToPi0Stat                     = (TH1D*)directoryPCMEta->Get("EtaToPi0YShiftedStatError");
        TGraphAsymmErrors* graphPCMEtaToPi0Stat        = new TGraphAsymmErrors(histoPCMEtaToPi0Stat);
        graphPCMEtaToPi0Stat->RemovePoint(graphPCMEtaToPi0Stat->GetN()-1);
        TGraphAsymmErrors* graphPCMEtaToPi0Sys         = (TGraphAsymmErrors*)directoryPCMEta->Get("EtaToPi0YShiftedSystError");
//        TH1D* histoPCMEtaInvXSectionStat                    = (TH1D*)directoryPCMEta->Get("InvCrossSectionEta");
//        TGraphAsymmErrors* graphPCMEtaInvXSectionStat       = new TGraphAsymmErrors(histoPCMEtaInvXSectionStat);
//        graphPCMEtaInvXSectionStat->RemovePoint(0);
//        cout << "Eta stat PCM" << endl;
//        graphPCMEtaInvXSectionStat->Print();

        TH1D* histoPCMEtaAccTimesEff                        = (TH1D*)histoPCMEtaTrueEffPt->Clone("histoPCMEtaAccTimesEff");
        histoPCMEtaAccTimesEff->Multiply(histoPCMEtaAcc);
        // normalize to full acceptance (delta y and phi)
        histoPCMEtaAccTimesEff->Scale(2*TMath::Pi()*1.6);

//        //      TGraphAsymmErrors* graphPCMEtaInvXSectionSysA     = (TGraphAsymmErrors*)directoryPCMEta->Get("InvCrossSectionEtaSysA");
//        TGraphAsymmErrors* graphPCMEtaInvXSectionSys        = (TGraphAsymmErrors*)directoryPCMEta->Get("InvCrossSectionEtaSys");
//        cout << "Eta sys PCM" << endl;
//        graphPCMEtaInvXSectionSys->Print();
////      TGraphAsymmErrors* graphPCMEtaCorrYieldSysErr    = (TGraphAsymmErrors*)directoryPCMEta->Get("EtaSystError");
//        TH1D* histoPCMEtaToPi0Stat                          = (TH1D*)directoryPCMEta->Get("EtatoPi0RatioConversionBinShifted");
//        TGraphAsymmErrors* graphPCMEtaToPi0Sys              = (TGraphAsymmErrors*)directoryPCMEta->Get("EtatoPi0RatioConversionBinShiftedSys");

        TH1D* histoEtaInvMassSigPlusBGPCM[3];
        TH1D* histoEtaInvMassSigPCM[3];
        TH1D* histoEtaInvMassSigRemBGSubPCM[3];
        TH1D* histoEtaInvMassBGPCM[3];
        TH1D* histoEtaInvMassRemBGPCM[3];
        TH1D* histoEtaInvMassBGMixedPCM[3];
        TF1* fitEtaInvMassSigPCM[3];
        TF1* fitEtaInvMassBGPCM[3];
        Bool_t haveAllEtaInvMassPCM[3]                 = {kFALSE, kFALSE, kFALSE};
        if (plotInvMassBins){
            for (Int_t i = 0; i < 1; i++){
                histoEtaInvMassSigPCM[i]               = (TH1D*)directoryPCMEta->Get("Eta_InvMassSig_Example_MB");
                histoEtaInvMassSigPlusBGPCM[i]         = (TH1D*)directoryPCMEta->Get("Eta_InvMassSigPlusBG_Example_MB");
                histoEtaInvMassBGPCM[i]                = (TH1D*)directoryPCMEta->Get("Eta_InvMassBG_Example_MB");
                fitEtaInvMassSigPCM[i]                 = (TF1*)directoryPCMEta->Get("Eta_InvMassSigFit_Example_MB");
                if (histoEtaInvMassSigPCM[i] && histoEtaInvMassSigPlusBGPCM[i] && histoEtaInvMassBGPCM[i] && fitEtaInvMassSigPCM[i]){
                    haveAllEtaInvMassPCM[i]            = kTRUE;
                }


                if (haveAllEtaInvMassPCM[i]){
                    histoEtaInvMassSigPCM[i]->Fit(fitEtaInvMassSigPCM[i],"QRME0");
                    for (Int_t l=0; l < 6; l++){
                        cout << fitEtaInvMassSigPCM[i]->GetParameter(l) << "\t +- " << fitEtaInvMassSigPCM[i]->GetParError(l) << endl;
                    }
                    fitEtaInvMassBGPCM[i]                                  = new TF1("Linearpp","[0]+[1]*x",0.0,1.0);
                    fitEtaInvMassBGPCM[i]->SetParameter(0, fitEtaInvMassSigPCM[i]->GetParameter(4));
                    fitEtaInvMassBGPCM[i]->SetParameter(1, fitEtaInvMassSigPCM[i]->GetParameter(5));
                    TVirtualFitter * fitterPCM                             = TVirtualFitter::GetFitter();
                    Int_t nFreeParPCM                                      = fitEtaInvMassSigPCM[i]->GetNumberFreeParameters();
                    double * covMatrixPCM                                  = fitterPCM->GetCovarianceMatrix();

                    histoEtaInvMassRemBGPCM[i]                             = (TH1D*)histoEtaInvMassBGPCM[i]->Clone(Form("Eta_InvMassRemBG_Example_%s",nameTrigger[i].Data()));
                    for (Int_t j = 1; j < histoEtaInvMassRemBGPCM[i]->GetNbinsX()+1; j++){
                        histoEtaInvMassRemBGPCM[i]->SetBinContent(j,0);
                        histoEtaInvMassRemBGPCM[i]->SetBinError(j,0);
                    }
                    for (Int_t j = histoEtaInvMassSigPCM[i]->GetXaxis()->FindBin(0.30); j < histoEtaInvMassSigPCM[i]->GetXaxis()->FindBin(0.70)+1; j++){
                        Double_t startBinEdge                                   = histoEtaInvMassSigPCM[i]->GetXaxis()->GetBinLowEdge(j);
                        Double_t endBinEdge                                     = histoEtaInvMassSigPCM[i]->GetXaxis()->GetBinUpEdge(j);
                        Double_t intLinearBack                                  = fitEtaInvMassBGPCM[i]->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
                        Double_t errorLinearBck                                 = pow(( pow( (endBinEdge-startBinEdge)*fitEtaInvMassSigPCM[i]->GetParError(4),2) +
                                                                                        pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitEtaInvMassSigPCM[i]->GetParError(5),2)
                                                                                        +2*covMatrixPCM[nFreeParPCM*nFreeParPCM-2]*(endBinEdge-startBinEdge)*0.5*
                                                                                        (endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
//                         cout << j << "\t" << intLinearBack << "\t" << errorLinearBck << endl;
                        histoEtaInvMassRemBGPCM[i]->SetBinContent(j,intLinearBack);
                        histoEtaInvMassRemBGPCM[i]->SetBinError(j,errorLinearBck);
                    }
                    histoEtaInvMassBGMixedPCM[i]         = (TH1D*)histoEtaInvMassBGPCM[i]->Clone(Form("Eta_InvMassTotBG_Example_%s",nameTrigger[i].Data()));
                    //histoEtaInvMassBGMixedPCM[i]->Sumw2();
                    //histoEtaInvMassBGMixedPCM[i]->Add(histoEtaInvMassRemBGPCM[i]);
                    histoEtaInvMassSigRemBGSubPCM[i]   = (TH1D*)histoEtaInvMassSigPCM[i]->Clone(Form("Eta_InvMassSigRemBGSub_Example_%s",nameTrigger[i].Data()));
                    histoEtaInvMassSigRemBGSubPCM[i]->Sumw2();
                    histoEtaInvMassSigRemBGSubPCM[i]->Add(histoEtaInvMassRemBGPCM[i],-1);
                    fitEtaInvMassSigPCM[i]->SetParameter(4, 0);
                    fitEtaInvMassSigPCM[i]->SetParameter(5, 0);
                }
                cout << nameTrigger[i].Data() << "\t" << histoEtaInvMassSigPCM[i] << "\t" << histoEtaInvMassSigPlusBGPCM[i] << "\t" << histoEtaInvMassBGPCM[i] << "\t" << fitEtaInvMassSigPCM[i]
                      << "\t" << haveAllEtaInvMassPCM[i] << endl;
            }
        }
    cout << "here" << endl;
    Int_t nEvtPCM                                           = histoPCMNumberOfEvents->GetBinContent(1);
//    cout << "here" << endl;
//    histoPCMPi0Mass->Scale(1000.);
//    histoPCMPi0TrueMass->Scale(1000.);
//    histoPCMEtaMass->Scale(1000.);
//    histoPCMEtaTrueMass->Scale(1000.);
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
        for (Int_t k = 0; k < 4; k++){
            graphPi0EffSecCorrFromX[k][4]                   = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get(Form("EffectiveSecondaryPi0CorrFrom%s",nameSecPi0SourceRead[k].Data()));
            if (graphPi0EffSecCorrFromX[k][4]){
                haveEffSecCorr[k][4]                        = kTRUE;
            }
        }
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
        TH1D* histoPi0InvMassBGMixedPCMEMCAL[3];
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
                    histoPi0InvMassBGMixedPCMEMCAL[i]         = (TH1D*)histoPi0InvMassBGPCMEMCAL[i]->Clone(Form("Pi0_InvMassTotBG_Example_%s",nameTrigger[i].Data()));
                    //histoPi0InvMassBGMixedPCMEMCAL[i]->Sumw2();
                    //histoPi0InvMassBGMixedPCMEMCAL[i]->Add(histoPi0InvMassRemBGPCMEMCAL[i]);
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

    TDirectory* directoryPCMEMCALEta                        = (TDirectory*)filePCMEMCAL->Get("Eta8TeV");
        TGraphAsymmErrors* graphPCMEMCALEtaMass             = (TGraphAsymmErrors*)directoryPCMEMCALEta->Get("Eta_Mass_data");
        graphPCMEMCALEtaMass                                = ScaleGraph(graphPCMEMCALEtaMass, 1000.);
        TGraphAsymmErrors* graphPCMEMCALEtaFWHM             = (TGraphAsymmErrors*)directoryPCMEMCALEta->Get("Eta_Width_data");
        graphPCMEMCALEtaFWHM                                = ScaleGraph(graphPCMEMCALEtaFWHM, 1000.);
        TGraphAsymmErrors* graphPCMEMCALEtaMassMC           = (TGraphAsymmErrors*)directoryPCMEMCALEta->Get("Eta_Mass_MC");
        graphPCMEMCALEtaMassMC                              = ScaleGraph(graphPCMEMCALEtaMassMC, 1000.);
        TGraphAsymmErrors* graphPCMEMCALEtaFWHMMC           = (TGraphAsymmErrors*)directoryPCMEMCALEta->Get("Eta_Width_MC");
        graphPCMEMCALEtaFWHMMC                              = ScaleGraph(graphPCMEMCALEtaFWHMMC, 1000.);
        TGraphAsymmErrors* graphPCMEMCALEtaAcc              = (TGraphAsymmErrors*)directoryPCMEMCALEta->Get("AcceptanceEta");
        TGraphAsymmErrors* graphPCMEMCALEtaEffPt            = (TGraphAsymmErrors*)directoryPCMEMCALEta->Get("EfficiencyEta");
        TGraphAsymmErrors* graphPCMEMCALEtaAccTimesEff      = (TGraphAsymmErrors*)directoryPCMEMCALEta->Get("EffTimesAccEta");
        TH1D* histoPCMEMCALEtaTriggerEff[4];
        for (Int_t i = 2; i < 6; i++){
            histoPCMEMCALEtaTriggerEff[i-2]                 = (TH1D*)directoryPCMEMCALEta->Get(Form("TriggerEfficiencyEta_%s",nameTrigger[i].Data()));
        }
        TH1D* histoPCMEMCALEtaInvXSectionStat               = (TH1D*)directoryPCMEMCALEta->Get("InvCrossSectionEta");
        TGraphAsymmErrors* graphPCMEMCALEtaInvXSectionStat  = (TGraphAsymmErrors*)directoryPCMEMCALEta->Get("graphInvCrossSectionEta");
        cout << "Eta stat PCM-EMC" << endl;
        graphPCMEMCALEtaInvXSectionStat->Print();
        TGraphAsymmErrors* graphPCMEMCALEtaInvXSectionSys   = (TGraphAsymmErrors*)directoryPCMEMCALEta->Get("InvCrossSectionEtaSys");
        cout << "Eta sys PCM-EMC" << endl;
        graphPCMEMCALEtaInvXSectionSys->Print();
        TH1D* histoPCMEMCALEtaToPi0Stat                     = (TH1D*)directoryPCMEMCALEta->Get("EtaToPi0YShiftedStatError");
        TGraphAsymmErrors* graphPCMEMCALEtaToPi0Stat        = new TGraphAsymmErrors(histoPCMEMCALEtaToPi0Stat);
        for(Int_t i=0; i<2; i++) graphPCMEMCALEtaToPi0Stat->RemovePoint(0);

        TGraphAsymmErrors* graphPCMEMCALEtaToPi0Sys         = (TGraphAsymmErrors*)directoryPCMEMCALEta->Get("EtaToPi0YShiftedSystError");


        TH1D* histoEtaInvMassSigPlusBGPCMEMCAL[3];
        TH1D* histoEtaInvMassSigPCMEMCAL[3];
        TH1D* histoEtaInvMassSigRemBGSubPCMEMCAL[3];
        TH1D* histoEtaInvMassBGPCMEMCAL[3];
        TH1D* histoEtaInvMassRemBGPCMEMCAL[3];
        TH1D* histoEtaInvMassBGMixedPCMEMCAL[3];
        TF1* fitEtaInvMassSigPCMEMCAL[3];
        TF1* fitEtaInvMassBGPCMEMCAL[3];
        Bool_t haveAllEtaInvMassPCMEMCAL[3]                 = {kFALSE, kFALSE, kFALSE};
        if (plotInvMassBins){
            for (Int_t i = 0; i < 3; i++){
                histoEtaInvMassSigPCMEMCAL[i]               = (TH1D*)directoryPCMEMCALEta->Get(Form("Eta_InvMassSig_Example_%s",nameTrigger[i].Data()));
                histoEtaInvMassSigPlusBGPCMEMCAL[i]         = (TH1D*)directoryPCMEMCALEta->Get(Form("Eta_InvMassSigPlusBG_Example_%s",nameTrigger[i].Data()));
                histoEtaInvMassBGPCMEMCAL[i]                = (TH1D*)directoryPCMEMCALEta->Get(Form("Eta_InvMassBG_Example_%s",nameTrigger[i].Data()));
                fitEtaInvMassSigPCMEMCAL[i]                 = (TF1*)directoryPCMEMCALEta->Get(Form("Eta_InvMassSigFit_Example_%s",nameTrigger[i].Data()));
                if (histoEtaInvMassSigPCMEMCAL[i] && histoEtaInvMassSigPlusBGPCMEMCAL[i] && histoEtaInvMassBGPCMEMCAL[i] && fitEtaInvMassSigPCMEMCAL[i]){
                    haveAllEtaInvMassPCMEMCAL[i]            = kTRUE;
                }


                if (haveAllEtaInvMassPCMEMCAL[i]){
                    histoEtaInvMassSigPCMEMCAL[i]->Fit(fitEtaInvMassSigPCMEMCAL[i],"QRME0");
                    for (Int_t l=0; l < 6; l++){
                        cout << fitEtaInvMassSigPCMEMCAL[i]->GetParameter(l) << "\t +- " << fitEtaInvMassSigPCMEMCAL[i]->GetParError(l) << endl;
                    }
                    fitEtaInvMassBGPCMEMCAL[i]                                  = new TF1("Linearpp","[0]+[1]*x",0.0,1.0);
                    fitEtaInvMassBGPCMEMCAL[i]->SetParameter(0, fitEtaInvMassSigPCMEMCAL[i]->GetParameter(4));
                    fitEtaInvMassBGPCMEMCAL[i]->SetParameter(1, fitEtaInvMassSigPCMEMCAL[i]->GetParameter(5));
                    TVirtualFitter * fitterPCMEMCAL                             = TVirtualFitter::GetFitter();
                    Int_t nFreeParPCMEMCAL                                      = fitEtaInvMassSigPCMEMCAL[i]->GetNumberFreeParameters();
                    double * covMatrixPCMEMCAL                                  = fitterPCMEMCAL->GetCovarianceMatrix();

                    histoEtaInvMassRemBGPCMEMCAL[i]                             = (TH1D*)histoEtaInvMassBGPCMEMCAL[i]->Clone(Form("Eta_InvMassRemBG_Example_%s",nameTrigger[i].Data()));
                    for (Int_t j = 1; j < histoEtaInvMassRemBGPCMEMCAL[i]->GetNbinsX()+1; j++){
                        histoEtaInvMassRemBGPCMEMCAL[i]->SetBinContent(j,0);
                        histoEtaInvMassRemBGPCMEMCAL[i]->SetBinError(j,0);
                    }
                    for (Int_t j = histoEtaInvMassSigPCMEMCAL[i]->GetXaxis()->FindBin(0.30); j < histoEtaInvMassSigPCMEMCAL[i]->GetXaxis()->FindBin(0.70)+1; j++){
                        Double_t startBinEdge                                   = histoEtaInvMassSigPCMEMCAL[i]->GetXaxis()->GetBinLowEdge(j);
                        Double_t endBinEdge                                     = histoEtaInvMassSigPCMEMCAL[i]->GetXaxis()->GetBinUpEdge(j);
                        Double_t intLinearBack                                  = fitEtaInvMassBGPCMEMCAL[i]->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
                        Double_t errorLinearBck                                 = pow(( pow( (endBinEdge-startBinEdge)*fitEtaInvMassSigPCMEMCAL[i]->GetParError(4),2) +
                                                                                        pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitEtaInvMassSigPCMEMCAL[i]->GetParError(5),2)
                                                                                        +2*covMatrixPCMEMCAL[nFreeParPCMEMCAL*nFreeParPCMEMCAL-2]*(endBinEdge-startBinEdge)*0.5*
                                                                                        (endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
//                         cout << j << "\t" << intLinearBack << "\t" << errorLinearBck << endl;
                        histoEtaInvMassRemBGPCMEMCAL[i]->SetBinContent(j,intLinearBack);
                        histoEtaInvMassRemBGPCMEMCAL[i]->SetBinError(j,errorLinearBck);
                    }
                    histoEtaInvMassBGMixedPCMEMCAL[i]         = (TH1D*)histoEtaInvMassBGPCMEMCAL[i]->Clone(Form("Eta_InvMassTotBG_Example_%s",nameTrigger[i].Data()));
                    //histoEtaInvMassBGMixedPCMEMCAL[i]->Sumw2();
                    //histoEtaInvMassBGMixedPCMEMCAL[i]->Add(histoEtaInvMassRemBGPCMEMCAL[i]);
                    histoEtaInvMassSigRemBGSubPCMEMCAL[i]   = (TH1D*)histoEtaInvMassSigPCMEMCAL[i]->Clone(Form("Eta_InvMassSigRemBGSub_Example_%s",nameTrigger[i].Data()));
                    histoEtaInvMassSigRemBGSubPCMEMCAL[i]->Sumw2();
                    histoEtaInvMassSigRemBGSubPCMEMCAL[i]->Add(histoEtaInvMassRemBGPCMEMCAL[i],-1);
                    fitEtaInvMassSigPCMEMCAL[i]->SetParameter(4, 0);
                    fitEtaInvMassSigPCMEMCAL[i]->SetParameter(5, 0);
                }
                cout << nameTrigger[i].Data() << "\t" << histoEtaInvMassSigPCMEMCAL[i] << "\t" << histoEtaInvMassSigPlusBGPCMEMCAL[i] << "\t" << histoEtaInvMassBGPCMEMCAL[i] << "\t" << fitEtaInvMassSigPCMEMCAL[i]
                      << "\t" << haveAllEtaInvMassPCMEMCAL[i] << endl;
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
        for (Int_t k = 0; k < 4; k++){
            graphPi0EffSecCorrFromX[k][2]                   = (TGraphAsymmErrors*)directoryEMCALPi0->Get(Form("EffectiveSecondaryPi0CorrFrom%s",nameSecPi0SourceRead[k].Data()));
            if (graphPi0EffSecCorrFromX[k][2]){
                haveEffSecCorr[k][2]                        = kTRUE;
            }
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
        TH1D* histoPi0InvMassBGMixedEMCAL[3];
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
                    histoPi0InvMassBGMixedEMCAL[i]         = (TH1D*)histoPi0InvMassBGEMCAL[i]->Clone(Form("Pi0_InvMassTotBG_Example_%s",nameTrigger[i].Data()));
                    //histoPi0InvMassBGMixedEMCAL[i]->Sumw2();
                    //histoPi0InvMassBGMixedEMCAL[i]->Add(histoPi0InvMassRemBGEMCAL[i]);
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


    TDirectory* directoryEMCALEta                           = (TDirectory*)fileEMCALLow->Get("Eta8TeV");
        TGraphAsymmErrors* graphEMCALEtaMass                = (TGraphAsymmErrors*)directoryEMCALEta->Get("Eta_Mass_data");
        graphEMCALEtaMass                                   = ScaleGraph(graphEMCALEtaMass, 1000.);
        TGraphAsymmErrors* graphEMCALEtaFWHM                = (TGraphAsymmErrors*)directoryEMCALEta->Get("Eta_Width_data");
        graphEMCALEtaFWHM                                   = ScaleGraph(graphEMCALEtaFWHM, 1000.);
        TGraphAsymmErrors* graphEMCALEtaMassMC              = (TGraphAsymmErrors*)directoryEMCALEta->Get("Eta_Mass_MC");
        graphEMCALEtaMassMC                                 = ScaleGraph(graphEMCALEtaMassMC, 1000.);
        TGraphAsymmErrors* graphEMCALEtaFWHMMC              = (TGraphAsymmErrors*)directoryEMCALEta->Get("Eta_Width_MC");
        graphEMCALEtaFWHMMC                                 = ScaleGraph(graphEMCALEtaFWHMMC, 1000.);
        TGraphAsymmErrors* graphEMCALEtaAcc                 = (TGraphAsymmErrors*)directoryEMCALEta->Get("AcceptanceEta");
        TGraphAsymmErrors* graphEMCALEtaEffPt               = (TGraphAsymmErrors*)directoryEMCALEta->Get("EfficiencyEta");
        TGraphAsymmErrors* graphEMCALEtaAccTimesEff         = (TGraphAsymmErrors*)directoryEMCALEta->Get("EffTimesAccEta");
        TH1D* histoEMCALEtaTriggerEff[4];
        for (Int_t i = 2; i < 6; i++){
            histoEMCALEtaTriggerEff[i-2]                    = (TH1D*)directoryEMCALEta->Get(Form("TriggerEfficiencyEta_%s",nameTrigger[i].Data()));
        }
        TH1D* histoEMCALEtaInvXSectionStat                  = (TH1D*)directoryEMCALEta->Get("InvCrossSectionEta");
        TGraphAsymmErrors* graphEMCALEtaInvXSectionStat     = (TGraphAsymmErrors*)directoryEMCALEta->Get("graphInvCrossSectionEta");
        cout << "Eta stat EMC-EMC" << endl;
        graphEMCALEtaInvXSectionStat->Print();
        TGraphAsymmErrors* graphEMCALEtaInvXSectionSys      = (TGraphAsymmErrors*)directoryEMCALEta->Get("InvCrossSectionEtaSys");
        cout << "Eta sys EMC-EMC" << endl;
        graphEMCALEtaInvXSectionSys->Print();
        TH1D* histoEMCALEtaToPi0Stat                        = (TH1D*)directoryEMCALEta->Get("EtaToPi0YShiftedStatError");
        TGraphAsymmErrors* graphEMCALEtaToPi0Stat           = new TGraphAsymmErrors(histoEMCALEtaToPi0Stat);
        for(Int_t i=0; i<5; i++) graphEMCALEtaToPi0Stat->RemovePoint(0);
        for(Int_t i=0; i<3; i++) graphEMCALEtaToPi0Stat->RemovePoint(graphEMCALEtaToPi0Stat->GetN()-1);

        TGraphAsymmErrors* graphEMCALEtaToPi0Sys            = (TGraphAsymmErrors*)directoryEMCALEta->Get("EtaToPi0YShiftedSystError");

        TH1D* histoEtaInvMassSigPlusBGEMCAL[3];
        TH1D* histoEtaInvMassSigEMCAL[3];
        TH1D* histoEtaInvMassSigRemBGSubEMCAL[3];
        TH1D* histoEtaInvMassBGEMCAL[3];
        TH1D* histoEtaInvMassRemBGEMCAL[3];
        TH1D* histoEtaInvMassBGMixedEMCAL[3];
        TF1* fitEtaInvMassSigEMCAL[3];
        TF1* fitEtaInvMassBGEMCAL[3];
        Bool_t haveAllEtaInvMassEMCAL[3]                 = {kFALSE, kFALSE, kFALSE};
        if (plotInvMassBins){
            for (Int_t i = 0; i < 3; i++){
                histoEtaInvMassSigEMCAL[i]                  = (TH1D*)directoryEMCALEta->Get(Form("Eta_InvMassSig_Example_%s",nameTrigger[i].Data()));
                histoEtaInvMassSigPlusBGEMCAL[i]            = (TH1D*)directoryEMCALEta->Get(Form("Eta_InvMassSigPlusBG_Example_%s",nameTrigger[i].Data()));
                histoEtaInvMassBGEMCAL[i]                   = (TH1D*)directoryEMCALEta->Get(Form("Eta_InvMassBG_Example_%s",nameTrigger[i].Data()));
                fitEtaInvMassSigEMCAL[i]                    = (TF1*)directoryEMCALEta->Get(Form("Eta_InvMassSigFit_Example_%s",nameTrigger[i].Data()));
                if (histoEtaInvMassSigEMCAL[i] && histoEtaInvMassSigPlusBGEMCAL[i] && histoEtaInvMassBGEMCAL[i] && fitEtaInvMassSigEMCAL[i]){
                    haveAllEtaInvMassEMCAL[i]               = kTRUE;
                }
                if (haveAllEtaInvMassEMCAL[i]){
                    histoEtaInvMassSigEMCAL[i]->Fit(fitEtaInvMassSigEMCAL[i],"QRME0");
                    for (Int_t l=0; l < 6; l++){
                        cout << fitEtaInvMassSigEMCAL[i]->GetParameter(l) << "\t +- " << fitEtaInvMassSigEMCAL[i]->GetParError(l) << endl;
                    }
                    fitEtaInvMassBGEMCAL[i]                                  = new TF1("Linearpp","[0]+[1]*x",0.0,1.0);
                    fitEtaInvMassBGEMCAL[i]->SetParameter(0, fitEtaInvMassSigEMCAL[i]->GetParameter(4));
                    fitEtaInvMassBGEMCAL[i]->SetParameter(1, fitEtaInvMassSigEMCAL[i]->GetParameter(5));
                    TVirtualFitter * fitterEMCAL                             = TVirtualFitter::GetFitter();
                    Int_t nFreeParEMCAL                                      = fitEtaInvMassSigEMCAL[i]->GetNumberFreeParameters();
                    double * covMatrixEMCAL                                  = fitterEMCAL->GetCovarianceMatrix();

                    histoEtaInvMassRemBGEMCAL[i]                             = (TH1D*)histoEtaInvMassBGEMCAL[i]->Clone(Form("Eta_InvMassRemBG_Example_%s",nameTrigger[i].Data()));
                    for (Int_t j = 1; j < histoEtaInvMassRemBGEMCAL[i]->GetNbinsX()+1; j++){
                        histoEtaInvMassRemBGEMCAL[i]->SetBinContent(j,0);
                        histoEtaInvMassRemBGEMCAL[i]->SetBinError(j,0);
                    }
                    for (Int_t j = histoEtaInvMassSigEMCAL[i]->GetXaxis()->FindBin(0.30); j < histoEtaInvMassSigEMCAL[i]->GetXaxis()->FindBin(0.70)+1; j++){
                        Double_t startBinEdge                                   = histoEtaInvMassSigEMCAL[i]->GetXaxis()->GetBinLowEdge(j);
                        Double_t endBinEdge                                     = histoEtaInvMassSigEMCAL[i]->GetXaxis()->GetBinUpEdge(j);
                        Double_t intLinearBack                                  = fitEtaInvMassBGEMCAL[i]->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
                        Double_t errorLinearBck                                 = pow(( pow( (endBinEdge-startBinEdge)*fitEtaInvMassSigEMCAL[i]->GetParError(4),2) +
                                                                                        pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitEtaInvMassSigEMCAL[i]->GetParError(5),2)
                                                                                        +2*covMatrixEMCAL[nFreeParEMCAL*nFreeParEMCAL-2]*(endBinEdge-startBinEdge)*0.5*
                                                                                        (endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
//                         cout << j << "\t" << intLinearBack << "\t" << errorLinearBck << endl;
                        histoEtaInvMassRemBGEMCAL[i]->SetBinContent(j,intLinearBack);
                        histoEtaInvMassRemBGEMCAL[i]->SetBinError(j,errorLinearBck);
                    }
                    histoEtaInvMassBGMixedEMCAL[i]         = (TH1D*)histoEtaInvMassBGEMCAL[i]->Clone(Form("Eta_InvMassTotBG_Example_%s",nameTrigger[i].Data()));
                    //histoEtaInvMassBGMixedEMCAL[i]->Sumw2();
                    //histoEtaInvMassBGMixedEMCAL[i]->Add(histoEtaInvMassRemBGEMCAL[i]);
                    histoEtaInvMassSigRemBGSubEMCAL[i]   = (TH1D*)histoEtaInvMassSigEMCAL[i]->Clone(Form("Eta_InvMassSigRemBGSub_Example_%s",nameTrigger[i].Data()));
                    histoEtaInvMassSigRemBGSubEMCAL[i]->Sumw2();
                    histoEtaInvMassSigRemBGSubEMCAL[i]->Add(histoEtaInvMassRemBGEMCAL[i],-1);
                    fitEtaInvMassSigEMCAL[i]->SetParameter(4, 0);
                    fitEtaInvMassSigEMCAL[i]->SetParameter(5, 0);
                }
                cout << nameTrigger[i].Data() << "\t" << histoEtaInvMassSigEMCAL[i] << "\t" << histoEtaInvMassSigPlusBGEMCAL[i] << "\t" << histoEtaInvMassBGEMCAL[i] << "\t" << fitEtaInvMassSigEMCAL[i]
                      << "\t" << haveAllEtaInvMassEMCAL[i] << endl;
            }
        }

    //************************** Read data for EMCAL merged **************************************************
    TFile* fileEMCALmerged                                  = new TFile(fileNameEMCALmerged.Data());
    TDirectory* directoryEMCALmergedPi0                     = (TDirectory*)fileEMCALmerged->Get("Pi08TeV");
//         TH1D* histoEMCALMergedPi0Eff                        = NULL;
//         TH1D* histoEMCALMergedPi0Pur                        = NULL;
        TH1D* histoEMCALMergedPi0AccTimesEff                = NULL;
        TH1D* histoEMCALMergedPi0InvXSectionStat            = NULL;
        TH1D* histoEMCALMergedPi0InvXSectionSys             = NULL;
        TGraphAsymmErrors* graphEMCALMergedPi0AccTimesEff       = NULL;
        TGraphAsymmErrors* graphEMCALMergedPi0Purity            = NULL;
        TGraphAsymmErrors* graphEMCALMergedPi0AccTimesEffDivPur = NULL;
        TGraphAsymmErrors* graphEMCALMergedPi0InvXSectionStat   = NULL;
        TGraphAsymmErrors* graphEMCALMergedPi0InvXSectionSys    = NULL;
//             histoEMCALMergedPi0Eff                          = (TH1D*)directoryEMCALmergedPi0->Get("EfficiencyPi0");
//             histoEMCALMergedPi0Pur                          = (TH1D*)directoryEMCALmergedPi0->Get("PurityPi0");
            graphEMCALMergedPi0AccTimesEff                  = (TGraphAsymmErrors*)directoryEMCALmergedPi0->Get("EffTimesAccPi0_EGA");
            //graphEMCALMergedPi0Purity                       = (TGraphAsymmErrors*)directoryEMCALmergedPi0->Get("PurityPi0");
            //graphEMCALMergedPi0AccTimesEffDivPur            = CalculateGraphErrRatioToOtherTGraphErr(graphEMCALMergedPi0AccTimesEff, graphEMCALMergedPi0Purity,kTRUE);
            graphEMCALMergedPi0InvXSectionStat              = (TGraphAsymmErrors*)directoryEMCALmergedPi0->Get("graphInvCrossSectionPi0");
            graphEMCALMergedPi0InvXSectionSys               = (TGraphAsymmErrors*)directoryEMCALmergedPi0->Get("InvCrossSectionPi0Sys");
                cout << "EMCAL merged stat" << endl;
                for (Int_t i = 0; graphEMCALMergedPi0InvXSectionStat->GetX()[0]< 10; i++){
                    graphEMCALMergedPi0InvXSectionStat->RemovePoint(0);
                }
                graphEMCALMergedPi0InvXSectionStat->Print();
                cout << "EMCAL merged sys" << endl;
                for (Int_t i = 0; graphEMCALMergedPi0InvXSectionSys->GetX()[0]< 10; i++){
                    graphEMCALMergedPi0InvXSectionSys->RemovePoint(0);
                }
                graphEMCALMergedPi0InvXSectionSys->Print();
//         return;

    //************************** Read data for PHOS *****************************************************

        TFile* filePHOS                                     = new TFile(fileNamePHOS.Data());
        TDirectory* directoryPHOSPi0                           = (TDirectory*)filePHOS->Get("Pi08TeV");
            TH1D* histoPHOSPi0Mass                = (TH1D*)directoryPHOSPi0->Get("MassPi0");
            histoPHOSPi0Mass->Scale(1000.);
              histoPHOSPi0Mass->SetBinContent(histoPHOSPi0Mass->FindBin(32.5),-200.);
            //for(Int_t k=0; k<6; k++) histoPHOSPi0Mass->SetBinContent(k,0.);
            //for(Int_t k=histoPHOSPi0Mass->GetNbinsX(), j=k-6; k>=j; k--) histoPHOSPi0Mass->SetBinContent(k,0.);

            TH1D* histoPHOSPi0FWHMMeV                = (TH1D*)directoryPHOSPi0->Get("FWHMPi0MeV");
              histoPHOSPi0FWHMMeV->SetBinContent(histoPHOSPi0FWHMMeV->FindBin(32.5),-200.);
            //for(Int_t k=0; k<6; k++) histoPHOSPi0FWHMMeV->SetBinContent(k,-1000.);
           // for(Int_t k=histoPHOSPi0FWHMMeV->GetNbinsX(), j=k-6; k>=j; k--) histoPHOSPi0FWHMMeV->SetBinContent(k,-1000.);

            TH1D* histoPHOSPi0TrueMass              = (TH1D*)directoryPHOSPi0->Get("TrueMassPi0");
            histoPHOSPi0TrueMass->Scale(1000.);
              histoPHOSPi0TrueMass->SetBinContent(histoPHOSPi0TrueMass->FindBin(32.5),-200.);
           // for(Int_t k=0; k<1; k++) histoPHOSPi0TrueMass->SetBinContent(k,0.);
           // for(Int_t k=histoPHOSPi0TrueMass->GetNbinsX(), j=k-5; k>=j; k--) histoPHOSPi0TrueMass->SetBinContent(k,0.);

            TH1D* histoPHOSPi0TrueFWHMMeV              = (TH1D*)directoryPHOSPi0->Get("TrueFWHMPi0MeV");
             histoPHOSPi0TrueFWHMMeV->SetBinContent(histoPHOSPi0TrueFWHMMeV->FindBin(32.5),-200.);
           // for(Int_t k=0; k<1; k++) histoPHOSPi0TrueFWHMMeV->SetBinContent(k,-1000.);
           // for(Int_t k=histoPHOSPi0TrueFWHMMeV->GetNbinsX(), j=k-5; k>=j; k--) histoPHOSPi0TrueFWHMMeV->SetBinContent(k,-1000.);

            TH1D* histoPHOSPi0Acc                 = (TH1D*)directoryPHOSPi0->Get("AcceptancePi0");
              histoPHOSPi0Acc->SetBinContent(histoPHOSPi0Acc->FindBin(32.5),-1.);
            //for(Int_t k=histoPHOSPi0Acc->GetNbinsX(), j=k-5; k>=j; k--) histoPHOSPi0Acc->SetBinContent(k,0.);

            TH1D* histoPHOSPi0TrueEffPt               = (TH1D*)directoryPHOSPi0->Get("EfficiencyPi0");
              histoPHOSPi0TrueEffPt->SetBinContent(histoPHOSPi0TrueEffPt->FindBin(32.5),-1.);
           // for(Int_t k=histoPHOSPi0TrueEffPt->GetNbinsX(), j=k-5; k>=j; k--) histoPHOSPi0TrueEffPt->SetBinContent(k,0.);

            TH1D* histoPHOSPi0AccTimesEff      = (TH1D*)histoPHOSPi0TrueEffPt->Clone("histoPHOSPi0AccTimesEff");
              histoPHOSPi0AccTimesEff->SetBinContent(histoPHOSPi0AccTimesEff->FindBin(32.5),-1.);
            histoPHOSPi0AccTimesEff->Multiply(histoPHOSPi0Acc);
            histoPHOSPi0AccTimesEff->Scale(2*TMath::Pi());
            //TH1D* graphPHOSPi0AccTimesEff         = (TH1D*)directoryPHOSPi0->Get("EffTimesAccPi0_MB");

            TH1D* histoPHOSPi0InvXSectionStat                   = (TH1D*)directoryPHOSPi0->Get("InvCrossSectionPi0");
              histoPHOSPi0InvXSectionStat->SetBinContent(histoPHOSPi0InvXSectionStat->FindBin(32.5),0.);
            TGraphAsymmErrors* graphPHOSPi0InvXSectionStat      = new TGraphAsymmErrors(histoPHOSPi0InvXSectionStat);
              graphPHOSPi0InvXSectionStat->RemovePoint(graphPHOSPi0InvXSectionStat->GetN()-1);
            cout << "Pi0 stat PHOS" << endl;
            graphPHOSPi0InvXSectionStat->Print();

            TGraphErrors* dummySysGraphPHOS                     = (TGraphErrors*)directoryPHOSPi0->Get("InvCrossSectionPi0Sys");
            Double_t* xValuePHOS                                = dummySysGraphPHOS->GetX();
            Double_t* yValuePHOS                                = dummySysGraphPHOS->GetY();
            Double_t* xErrPHOS                                  = dummySysGraphPHOS->GetEX();
            Double_t* yErrPHOS                                  = dummySysGraphPHOS->GetEY();

            TGraphAsymmErrors* graphPHOSPi0InvXSectionSys       = new TGraphAsymmErrors(dummySysGraphPHOS->GetN(), xValuePHOS, yValuePHOS, xErrPHOS, xErrPHOS, yErrPHOS, yErrPHOS);
            cout << "Pi0 sys PHOS" << endl;
              graphPHOSPi0InvXSectionSys->RemovePoint(graphPHOSPi0InvXSectionSys->GetN()-1);
            graphPHOSPi0InvXSectionSys->Print();


        TFile* filePHOSprelim                                     = new TFile(fileNamePHOSprelim.Data());

        //***************************************************************************************************************
        //*******************************Plotting trigger rejection factors = fits log scale all in one *****************
        //***************************************************************************************************************

        TFile* filePHOSPHOS                                     = new TFile(fileNamePHOSPHOS.Data());
        TDirectory* directoryPHOSPHOSRejection                  = (TDirectory*)filePHOSPHOS->Get("RejectionFactor");

        TString triggerNameLabel[3] = {"EMC-L0/MB","EMC-L1/EMC-L0","PHOS-L0/MB"};
        TH1F* histoTriggerReject[3];
//        TH1D* triggRejec[3];
        histoTriggerReject[0] = (TH1F*)fileEMCALLow->Get("TriggRejectvsE_EMC7_INT7");
        histoTriggerReject[1] = (TH1F*)fileEMCALLow->Get("TriggRejectvsE_EGA_EMC7");
        histoTriggerReject[1]->GetXaxis()->SetRangeUser(2.3,50.);
        histoTriggerReject[2] = (TH1F*)directoryPHOSPHOSRejection->Get("fHistRejection");

        Double_t minPt[3] = {4.1,12.5,6.};
        Double_t maxPt[3] = {30.,50.,24.};
        Double_t triggRejecFac[3] = {67.0,222.51,12.4E3};
        Double_t triggRejecFacErr[3] = {1.1,4.02,1.5E3};

        Size_t textSizeSpectra2         = 0.0415;
        Int_t textPixelPP               = textSizeSpectra2*1100;
        TCanvas* canvasTriggerReject    = new TCanvas("canvasTriggerReject","",0,0,1500,1100);// gives the page size
        DrawGammaCanvasSettings( canvasTriggerReject, 0.076, 0.015, 0.015, 0.085);
        canvasTriggerReject->SetLogy(1);

        Double_t minTriggReject = 0.1;
        Double_t maxTriggReject = 60000;

        TH2F * histo2DTriggReject = new TH2F("histo2DTriggReject","histo2DTriggReject",1000,0., 50.,10000,minTriggReject, maxTriggReject);
        SetStyleHistoTH2ForGraphs(histo2DTriggReject, "#it{E} (GeV)","#it{RF}", //"#frac{N_{clus,trig A}/N_{Evt, trig A}}{N_{clus,trig B}/N_{Evt,trig B}}",
                                0.85*textSizeSpectra2,textSizeSpectra2, 0.85*textSizeSpectra2,textSizeSpectra2, 0.85,0.85);
        histo2DTriggReject->DrawCopy();

        TLegend* legendTriggReject = GetAndSetLegend2(0.2, 0.12, 0.9, 0.12+(0.9*3*textSizeSpectra2),textPixelPP);
        legendTriggReject->SetMargin(0.02);
        legendTriggReject->SetNColumns(3);
        for (Int_t i = 0; i< 3; i++){
          Color_t color, colorShade = kBlack;
          if(i==0){
            color = kBlue+2;
            colorShade = kBlue-6;
            DrawGammaSetMarker(histoTriggerReject[i], markerTrigg[i], sizeTrigg[4], colorDet[4], colorDet[4]);
            histoTriggerReject[i]->DrawCopy("e,same");
            legendTriggReject->AddEntry(histoTriggerReject[i],Form("   %s",triggerNameLabel[i].Data()),"p");
            legendTriggReject->AddEntry((TObject*)0,Form("       %3.1f < E < %3.1f",minPt[i],maxPt[i]),"");
            legendTriggReject->AddEntry((TObject*)0,Form("               %3.1f #pm %3.1f",triggRejecFac[i], triggRejecFacErr[i]),"");
          }else if(i==1){
            color = kGreen+3;
            colorShade = kGreen-8;
            DrawGammaSetMarker(histoTriggerReject[i], markerTrigg[i], sizeTrigg[4], colorDet[2], colorDet[2]);
            histoTriggerReject[i]->DrawCopy("e,same");
            legendTriggReject->AddEntry(histoTriggerReject[i],Form("   %s",triggerNameLabel[i].Data()),"p");
            legendTriggReject->AddEntry((TObject*)0,Form("     %3.1f < E < %3.1f",minPt[i],maxPt[i]),"");
            legendTriggReject->AddEntry((TObject*)0,Form("             %3.1f #pm %3.1f",triggRejecFac[i], triggRejecFacErr[i]),"");
          }else if(i==2){
            color = kRed+2;
            colorShade = kRed-6;
            DrawGammaSetMarker(histoTriggerReject[i], markerTrigg[i], sizeTrigg[4], colorDet[1], colorDet[1]);
            histoTriggerReject[i]->DrawCopy("e,same");
            legendTriggReject->AddEntry(histoTriggerReject[i],Form("   %s",triggerNameLabel[i].Data()),"p");
            legendTriggReject->AddEntry((TObject*)0,Form("       %3.1f < E < %3.1f",minPt[i],maxPt[i]),"");
            legendTriggReject->AddEntry((TObject*)0,Form("              (%3.1f #pm %3.1f)#times10^{3}",triggRejecFac[i]/1E3, triggRejecFacErr[i]/1E3),"");
          }
          TF1* pol0 = new TF1("pol0","[0]",minPt[i],maxPt[i]); //
          histoTriggerReject[i]->Fit(pol0,"NRME+","",minPt[i],maxPt[i]);

//          triggRejec[i] = (TH1D*)histoTriggerReject[i]->Clone(Form("CL_%i",i));
//          for (Int_t j = 1; j < triggRejec[i]->GetNbinsX()+1; j++){
//              triggRejec[i]->SetBinContent(j,triggRejecFac[i]);
//              triggRejec[i]->SetBinError(j,triggRejecFacErr[i]);
//          }
//          triggRejec[i]->SetStats(kFALSE);
//          triggRejec[i]->SetFillColor(colorShade);
//          triggRejec[i]->SetMarkerSize(0);
//          triggRejec[i]->Draw("e3,same");
          TBox* box = new TBox(0 ,triggRejecFac[i]-triggRejecFacErr[i] , histoTriggerReject[i]->GetXaxis()->GetBinUpEdge(histoTriggerReject[i]->GetNbinsX()), triggRejecFac[i]+triggRejecFacErr[i]);
          box->SetLineColor(colorShade);
          box->SetFillColorAlpha(colorShade,0.1);
          box->Draw();

          pol0->SetParameter(0,triggRejecFac[i]);
          pol0->SetParError(0,triggRejecFacErr[i]);
          pol0->SetLineColor(color);
          pol0->SetLineStyle(7);
          pol0->SetRange(minPt[i],maxPt[i]);
          pol0->Draw("same");
          histoTriggerReject[i]->DrawCopy("e,same");
        }

        legendTriggReject->Draw();
        histo2DTriggReject->Draw("same,axis");
        TLatex *labelPerfTriggRejec = new TLatex(0.7, 0.92,ALICEperfor.Data());
        SetStyleTLatex( labelPerfTriggRejec, textSizeSpectra2,4);
        labelPerfTriggRejec->Draw();

        TLatex *labelWeightsTR      = new TLatex(0.7,0.87,collisionSystem8TeV.Data());
        SetStyleTLatex( labelWeightsTR, textSizeSpectra2,4);
        labelWeightsTR->Draw();

        if(plotDate){
          TLatex *labelDate = new TLatex(0.7, 0.82,date.Data());
          SetStyleTLatex( labelDate, textSizeSpectra2,4);
          labelDate->Draw();
        }

        TLatex *labelPerfTriggFitRange = new TLatex(0.463, 0.125+(0.9*3*textSizeSpectra2)+0.01, "Fit range (GeV)");
        SetStyleTLatex( labelPerfTriggFitRange, textSizeSpectra2,4);
        labelPerfTriggFitRange->Draw();

        TLatex *labelPerfTriggRejecFac = new TLatex(0.723, 0.125+(0.9*3*textSizeSpectra2)+0.01, "Trigger rejection");
        SetStyleTLatex( labelPerfTriggRejecFac, textSizeSpectra2,4);
        labelPerfTriggRejecFac->Draw();

        canvasTriggerReject->Update();
        canvasTriggerReject->SaveAs(Form("%s/TriggerRejectionFactors.%s",outputDir.Data(),suffix.Data()));
        delete canvasTriggerReject;

    // *******************************************************************************************************
    // ************************** Loading theory calculations ************************************************
    // *******************************************************************************************************

    TFile* fileChargedIdentified                   = new TFile(fileNameChargedPionPP.Data());
    TH1F* histoChPion2760GeVStat                   = (TH1F*) fileChargedIdentified->Get("histoChargedPionSpecPubStat2760GeV");
    TH1F* histoChPion2760GeVSys                    = (TH1F*) fileChargedIdentified->Get("histoChargedPionSpecPubSyst2760GeV");
    histoChPion2760GeVStat->Scale(xSection2760GeVINEL);
    histoChPion2760GeVSys->Scale(xSection2760GeVINEL);
    TH1F* histoChPion2760GeV              = new TH1F(*histoChPion2760GeVStat);
    for (Int_t i = 0; i<histoChPion2760GeV->GetNbinsX(); i++){
        histoChPion2760GeV->SetBinError(i, TMath::Sqrt(TMath::Power(histoChPion2760GeVStat->GetBinError(i),2)+TMath::Power(histoChPion2760GeVSys->GetBinError(i),2)));
    }

    TH1F* histoChPion7TeVStat                      = (TH1F*) fileChargedIdentified->Get("histoChargedPionSpecPubStat7TeV");
    TH1F* histoChPion7TeVSys                       = (TH1F*) fileChargedIdentified->Get("histoChargedPionSpecPubSyst7TeV");
    histoChPion7TeVStat->Scale(xSection7TeVINEL*1e12);
    histoChPion7TeVSys->Scale(xSection7TeVINEL*1e12);
    TH1F* histoChPion7TeV              = new TH1F(*histoChPion7TeVStat);
    for (Int_t i = 0; i<histoChPion2760GeV->GetNbinsX(); i++){
        histoChPion7TeV->SetBinError(i,TMath::Sqrt(TMath::Power(histoChPion7TeVStat->GetBinError(i),2)+TMath::Power(histoChPion7TeVSys->GetBinError(i),2)));
    }

    TFile* fileEtaToPi                              = new TFile(fileNameEtaToPi0.Data());
    TGraphErrors *eta2pi0_NA27_275GeV               = (TGraphErrors*)fileEtaToPi->Get("Aguilar400GeV");
    TGraphErrors *eta2pi0_RHIC200GeV                = (TGraphErrors*)fileEtaToPi->Get("Phenix200GeV");

    TFile* fileTheoryCompilation                            = new TFile(fileNameTheory.Data());
    TFile* filePOWHEG8TeV                                   = new TFile(fileNamePOWHEG8TeV.Data());

        // Pythia8 Monash2013:
        TH1F* histoPythia8InvXSection                       = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoPi08TeV");
        histoPythia8InvXSection->GetXaxis()->SetRangeUser(0.3,35);
        TH1F* histoPythia8InvXSectionEta                    = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoEta8TeV");
        histoPythia8InvXSectionEta->GetXaxis()->SetRangeUser(0.4,35);
        TH1F* histoPythia8InvXSectionChPion             = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoChPion8TeV");
        histoPythia8InvXSectionChPion->GetXaxis()->SetRangeUser(0.3,35);
        TGraphErrors* graphPythia8InvXSection               = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoPi08TeV"));
        while(graphPythia8InvXSection->GetX()[0] < 0.3) graphPythia8InvXSection->RemovePoint(0);
        while(graphPythia8InvXSection->GetX()[graphPythia8InvXSection->GetN()-1] > 35.) graphPythia8InvXSection->RemovePoint(graphPythia8InvXSection->GetN()-1);
        TGraphErrors* graphPythia8InvXSectionEta            = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoEta8TeV"));
        while(graphPythia8InvXSectionEta->GetX()[0] < 0.4) graphPythia8InvXSectionEta->RemovePoint(0);
        while(graphPythia8InvXSectionEta->GetX()[graphPythia8InvXSectionEta->GetN()-1] > 35.) graphPythia8InvXSectionEta->RemovePoint(graphPythia8InvXSectionEta->GetN()-1);
        TH1F* histoPythia8EtaToPi0                          = (TH1F*) histoPythia8InvXSectionEta->Clone("Pythia8EtaToPi0");
        histoPythia8EtaToPi0->Divide(histoPythia8InvXSection);
        histoPythia8EtaToPi0->GetXaxis()->SetRangeUser(0.4,25);
        TGraphErrors* graphPythia8EtaToPi0                  = new TGraphErrors(histoPythia8EtaToPi0);
        while(graphPythia8EtaToPi0->GetX()[0] < 0.4) graphPythia8EtaToPi0->RemovePoint(0);
        while(graphPythia8EtaToPi0->GetX()[graphPythia8EtaToPi0->GetN()-1] > 25.) graphPythia8EtaToPi0->RemovePoint(graphPythia8EtaToPi0->GetN()-1);
        // *******************************************************************************************************
        // Pythia8 Tune4C:
        TH1F* histoPythia8_4CInvXSection                    = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Tune4CLegoPi08TeV");
        histoPythia8_4CInvXSection->GetXaxis()->SetRangeUser(0.3,35);
        TH1F* histoPythia8_4CInvXSectionEta                 = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Tune4CLegoEta8TeV");
        histoPythia8_4CInvXSectionEta->GetXaxis()->SetRangeUser(0.4,35);
        TGraphErrors* graphPythia8_4CInvXSection            = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Tune4CLegoPi08TeV"));
        while(graphPythia8_4CInvXSection->GetX()[0] < 0.3) graphPythia8_4CInvXSection->RemovePoint(0);
        while(graphPythia8_4CInvXSection->GetX()[graphPythia8_4CInvXSection->GetN()-1] > 35.) graphPythia8_4CInvXSection->RemovePoint(graphPythia8_4CInvXSection->GetN()-1);
        TGraphErrors* graphPythia8_4CInvXSectionEta         = new TGraphErrors((TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Tune4CLegoEta8TeV"));
        while(graphPythia8_4CInvXSectionEta->GetX()[0] < 0.4) graphPythia8_4CInvXSectionEta->RemovePoint(0);
        while(graphPythia8_4CInvXSectionEta->GetX()[graphPythia8_4CInvXSectionEta->GetN()-1] > 35.) graphPythia8_4CInvXSectionEta->RemovePoint(graphPythia8_4CInvXSectionEta->GetN()-1);
        TH1F* histoPythia8T4CEtaToPi0                       = (TH1F*) histoPythia8_4CInvXSectionEta->Clone("Pythia8T4CEtaToPi0");
        histoPythia8T4CEtaToPi0->Divide(histoPythia8_4CInvXSection);
        histoPythia8T4CEtaToPi0->GetXaxis()->SetRangeUser(0.4,25);
        TGraphErrors* graphPythia8T4CEtaToPi0               = new TGraphErrors(histoPythia8T4CEtaToPi0);
        while(graphPythia8T4CEtaToPi0->GetX()[0] < 0.4) graphPythia8T4CEtaToPi0->RemovePoint(0);
        while(graphPythia8T4CEtaToPi0->GetX()[graphPythia8T4CEtaToPi0->GetN()-1] > 25.) graphPythia8T4CEtaToPi0->RemovePoint(graphPythia8T4CEtaToPi0->GetN()-1);

        // *******************************************************************************************************
        // NLO calc
        TGraphAsymmErrors* graphPi0DSS14                    = (TGraphAsymmErrors*) fileTheoryCompilation->Get("graphNLOCalcDSS14InvSecPi08000GeV");
        while (graphPi0DSS14->GetX()[graphPi0DSS14->GetN()-1] > 42. ) graphPi0DSS14->RemovePoint(graphPi0DSS14->GetN()-1);

        TGraphAsymmErrors* graphEtaAESSS                    = (TGraphAsymmErrors*) fileTheoryCompilation->Get("graphNLOCalcAESSSInvSecEta8000GeV");

        TGraphAsymmErrors* graphEtaToPi07000GeV             = (TGraphAsymmErrors*) fileEtaToPi->Get("Alice7TeV");
        TGraphAsymmErrors* graphEtaToPi02760GeV             = (TGraphAsymmErrors*) fileEtaToPi->Get("Alice2760GeV");
        ProduceGraphAsymmWithoutXErrors(graphEtaToPi02760GeV);

        // *******************************************************************************************************
        TGraphAsymmErrors* graphPOWHEGPi08TeV               = (TGraphAsymmErrors*) filePOWHEG8TeV->Get("pi0");
        while (graphPOWHEGPi08TeV->GetX()[graphPOWHEGPi08TeV->GetN()-1] > 35. ) graphPOWHEGPi08TeV->RemovePoint(graphPOWHEGPi08TeV->GetN()-1);
        while (graphPOWHEGPi08TeV->GetX()[0] < 5. ) graphPOWHEGPi08TeV->RemovePoint(0);
        // *******************************************************************************************************

        TGraph* graphNLOCalcPi0MuHalf                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuHalf8000GeV");
        TGraph* graphNLOCalcPi0MuOne                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuOne8000GeV");
        TGraph* graphNLOCalcPi0MuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuTwo8000GeV");
        TGraph* graphNLOCalcEtaMuHalf                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuHalf8000GeV");
        TGraph* graphNLOCalcEtaMuOne                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuOne8000GeV");
        TGraph* graphNLOCalcEtaMuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuTwo8000GeV");

        while (graphNLOCalcEtaMuHalf->GetX()[graphNLOCalcEtaMuHalf->GetN()-1] > 37. )
            graphNLOCalcEtaMuHalf->RemovePoint(graphNLOCalcEtaMuHalf->GetN()-1);
        while (graphNLOCalcEtaMuOne->GetX()[graphNLOCalcEtaMuOne->GetN()-1] > 37. )
            graphNLOCalcEtaMuOne->RemovePoint(graphNLOCalcEtaMuOne->GetN()-1);
        while (graphNLOCalcEtaMuTwo->GetX()[graphNLOCalcEtaMuTwo->GetN()-1] > 37. )
            graphNLOCalcEtaMuTwo->RemovePoint(graphNLOCalcEtaMuTwo->GetN()-1);

        TGraphAsymmErrors* graphNLOEtaToPi0                = (TGraphAsymmErrors*) fileTheoryCompilation->Get("graphNLOCalcEtaOverPi08000GeV_AESSS_DSS07");
        while (graphNLOEtaToPi0->GetX()[graphNLOEtaToPi0->GetN()-1] > 27. ) graphNLOEtaToPi0->RemovePoint(graphNLOEtaToPi0->GetN()-1);


     TFile* file2760GeV                              = new TFile(fileName2760GeV.Data());
     TGraphAsymmErrors* graph2760GeVPi0              = (TGraphAsymmErrors*) file2760GeV->Get("Pi02.76TeV/graphInvCrossSectionPi0Comb2760GeVATotErr");
     TGraphAsymmErrors* graph2760GeVPi0Stat           = (TGraphAsymmErrors*) file2760GeV->Get("Pi02.76TeV/graphInvCrossSectionPi0Comb2760GeVAStatErr");
     TGraphAsymmErrors* graph2760GeVPi0Sys           = (TGraphAsymmErrors*) file2760GeV->Get("Pi02.76TeV/graphInvCrossSectionPi0Comb2760GeVASysErr");
     TF1* fit2760GeVPi0TCM                           = (TF1*) file2760GeV->Get("Pi02.76TeV/TwoComponentModelFitPi0");
     fit2760GeVPi0TCM->SetName("fit2760GeVPi0TCM");
     TF1* fit2760GeVPi0Tsallis                           = (TF1*) file2760GeV->Get("Pi02.76TeV/TsallisFitPi0");
     fit2760GeVPi0Tsallis->SetName("fit2760GeVPi0Tsallis");

     TGraphAsymmErrors* graph2760GeVEta              = (TGraphAsymmErrors*) file2760GeV->Get("Eta2.76TeV/graphInvCrossSectionEtaComb2760GeVATotErr");
     TGraphAsymmErrors* graph2760GeVEtaStat           = (TGraphAsymmErrors*) file2760GeV->Get("Eta2.76TeV/graphInvCrossSectionEtaComb2760GeVAStatErr");
     TGraphAsymmErrors* graph2760GeVEtaSys           = (TGraphAsymmErrors*) file2760GeV->Get("Eta2.76TeV/graphInvCrossSectionEtaComb2760GeVASysErr");
     TF1* fit2760GeVEtaTCM                           = (TF1*) file2760GeV->Get("Eta2.76TeV/TwoComponentModelFitEta");
     fit2760GeVEtaTCM->SetName("fit2760GeVEtaTCM");
     TF1* fit2760GeVEtaTsallis                       = (TF1*) file2760GeV->Get("Eta2.76TeV/TsallisFitEta");
     fit2760GeVEtaTsallis->SetName("fit2760GeVEtaTsallis");

     TH1F* histoPythia8InvXSection2760GeV                       = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoPi02760GeV");
     histoPythia8InvXSection2760GeV->GetXaxis()->SetRangeUser(0.3,35);
     TH1F* histoPythia8InvXSection2760GeVEta                    = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoEta2760GeV");
     histoPythia8InvXSection2760GeVEta->GetXaxis()->SetRangeUser(0.4,35);
     TH1F* histoPythia8InvXSection2760GeVChPion                 = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoChPion2760GeV");
     histoPythia8InvXSection2760GeVChPion->GetXaxis()->SetRangeUser(0.3,35);

     TFile* file7TeV                              = new TFile(fileName7TeV.Data());
     TGraphAsymmErrors* graph7TeVPi0Stat          = (TGraphAsymmErrors*) file7TeV->Get("Pi07TeV/graphInvCrossSectionPi0Comb7TeVStatErr");
     TGraphAsymmErrors* graph7TeVPi0Sys           = (TGraphAsymmErrors*) file7TeV->Get("Pi07TeV/graphInvCrossSectionPi0Comb7TeVSysErr");
     TGraphAsymmErrors* graph7TeVPi0              = (TGraphAsymmErrors*) file7TeV->Get("Pi07TeV/graphInvCrossSectionPi0Comb7TeV");
     if(!graph7TeVPi0){
       for (Int_t i = 0; i<graph7TeVPi0->GetN(); i++){
           graph7TeVPi0->GetEYlow()[i] = TMath::Sqrt(TMath::Power(graph7TeVPi0Stat->GetEYlow()[i],2)+TMath::Power(graph7TeVPi0Sys->GetEYlow()[i],2));
           graph7TeVPi0->GetEYhigh()[i] = TMath::Sqrt(TMath::Power(graph7TeVPi0Stat->GetEYhigh()[i],2)+TMath::Power(graph7TeVPi0Sys->GetEYhigh()[i],2));
       }
     }

//     TFile* file7TeVpub                              = new TFile(fileName7TeVpub.Data());
//     TGraphAsymmErrors* graph7TeVPi0StatPub          = (TGraphAsymmErrors*) file7TeVpub->Get("graphInvCrossSectionPi0Comb7TeVStatErr");
//     TGraphAsymmErrors* graph7TeVPi0SysPub           = (TGraphAsymmErrors*) file7TeVpub->Get("graphInvCrossSectionPi0Comb7TeVSysErr");

     TF1* fit7TeVPi0TCM                           = (TF1*) file7TeV->Get("Pi07TeV/TwoComponentModelFitPi0");
     fit7TeVPi0TCM->SetName("fit7TeVPi0TCM");
     TF1* fit7TeVPi0Tsallis                           = (TF1*) file7TeV->Get("Pi07TeV/TsallisFitPi0");
     fit7TeVPi0Tsallis->SetName("fit7TeVPi0Tsallis");

     TGraphAsymmErrors* graph7TeVEtaStat          = (TGraphAsymmErrors*) file7TeV->Get("Eta7TeV/graphInvCrossSectionEtaComb7TeVStatErr");
     TGraphAsymmErrors* graph7TeVEtaSys           = (TGraphAsymmErrors*) file7TeV->Get("Eta7TeV/graphInvCrossSectionEtaComb7TeVSysErr");
     TGraphAsymmErrors* graph7TeVEta              = (TGraphAsymmErrors*) file7TeV->Get("Eta7TeV/graphInvCrossSectionEtaComb7TeV");
     if(!graph7TeVEta){
       for (Int_t i = 0; i<graph7TeVEta->GetN(); i++){
           graph7TeVEta->GetEYlow()[i] = TMath::Sqrt(TMath::Power(graph7TeVEtaStat->GetEYlow()[i],2)+TMath::Power(graph7TeVEtaSys->GetEYlow()[i],2));
           graph7TeVEta->GetEYhigh()[i] = TMath::Sqrt(TMath::Power(graph7TeVEtaStat->GetEYhigh()[i],2)+TMath::Power(graph7TeVEtaSys->GetEYhigh()[i],2));
       }
     }
     TF1* fit7TeVEtaTCM                           = (TF1*) file7TeV->Get("Eta7TeV/TwoComponentModelFitEta");
     fit7TeVEtaTCM->SetName("fit7TeVEtaTCM");
     TF1* fit7TeVEtaTsallis                       = (TF1*) file7TeV->Get("Eta7TeV/TsallisFitEta");
     fit7TeVEtaTsallis->SetName("fit7TeVEtaTsallis");

     TH1F* histoPythia8InvXSection7TeV                       = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoPi07TeV");
     histoPythia8InvXSection7TeV->GetXaxis()->SetRangeUser(0.3,35);
     TH1F* histoPythia8InvXSection7TeVEta                    = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoEta7TeV");
     histoPythia8InvXSection7TeVEta->GetXaxis()->SetRangeUser(0.4,35);
     TH1F* histoPythia8InvXSection7TeVChPion                 = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoChPion7TeV");
     histoPythia8InvXSection7TeVChPion->GetXaxis()->SetRangeUser(0.3,35);

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
//    TH1D* histoEMCALMergedPi0InvXSectionStatCorrBin     = new TH1D("histoEMCALMergedPi0InvXSectionStatCorrBin", "", 33, xPtLimitsPi0);
//    Int_t firstBinMergedPi0 = 1;
//    cout << graphEMCALMergedPi0InvXSectionStat->GetX()[0] << endl;

//    while (histoEMCALMergedPi0InvXSectionStatCorrBin->GetBinCenter(firstBinMergedPi0) < graphEMCALMergedPi0InvXSectionStat->GetX()[0]){
//        cout << histoEMCALMergedPi0InvXSectionStatCorrBin->GetBinCenter(firstBinMergedPi0) << endl;
//        histoEMCALMergedPi0InvXSectionStatCorrBin->SetBinContent(firstBinMergedPi0, 0);
//        histoEMCALMergedPi0InvXSectionStatCorrBin->SetBinError(firstBinMergedPi0, 0);
//        firstBinMergedPi0++;
//    }
//    for (Int_t i = 0; i < graphEMCALMergedPi0InvXSectionStat->GetN(); i++){
//        histoEMCALMergedPi0InvXSectionStatCorrBin->SetBinContent(i+firstBinMergedPi0, graphEMCALMergedPi0InvXSectionStat->GetY()[i]);
//        histoEMCALMergedPi0InvXSectionStatCorrBin->SetBinError(i+firstBinMergedPi0, graphEMCALMergedPi0InvXSectionStat->GetEYlow()[i]);
//    }

//    graphEMCALMergedPi0InvXSectionStat->Print();
//    for (Int_t i = 1; i < histoEMCALMergedPi0InvXSectionStatCorrBin->GetNbinsX(); i++){
//        cout << "Bin " << i << "\t" <<  histoEMCALMergedPi0InvXSectionStatCorrBin->GetBinCenter(i) << "\t" << histoEMCALMergedPi0InvXSectionStatCorrBin->GetBinContent(i)
//             << "\t" << histoEMCALMergedPi0InvXSectionStatCorrBin->GetBinError(i) << endl;
//    }

//     return;

//    statErrorCollectionPi0[9]          = (TH1D*)histoEMCALMergedPi0InvXSectionStatCorrBin->Clone("statErrEMCALMergedPi0");
//    sysErrorCollectionPi0[9]           = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionSys->Clone("sysErrEMCALMergedPi0");


    // Definition of final pt binning (has to be set manually)
    Double_t xPtLimitsPi0[100];
    Int_t maxNBinsPi0               = 0;

    maxNBinsPi0                  = GetBinning( xPtLimitsPi0, "Pi0", "8TeV", 2 );
    maxNBinsPi0--;

    Double_t xPtLimitsPi0WOMerged[70];
    Int_t maxNBinsPi0W0Merged       = GetBinning( xPtLimitsPi0WOMerged, "Pi0", "8TeV", 2 );

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsPi0[11]            = { 0,  6,  0,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t offSetsPi0Sys[11]         = { 1,  6,  7,  0,  5,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t offSetPi0Shifting[11]     = { 0,  0,  0,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0 };
    Int_t nComBinsPi0Shifting[11]   = { 0,  0,  0,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0 };

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
    Int_t nMeasSetPi0A                 = 4;
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
    histo2DPi0Weights = new TH2F("histo2DPi0Weights","histo2DPi0Weights",11000,0.23,50.,1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DPi0Weights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DPi0Weights->GetXaxis()->SetMoreLogLabels();
    histo2DPi0Weights->GetXaxis()->SetNoExponent(kTRUE);
    //histo2DPi0Weights->GetXaxis()->SetLabelOffset(-0.01);
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

//      DrawGammaLines(0.23, 50. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.23, 50. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.23, 50. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.23, 50. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.23, 50. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/Pi0_WeightsA.%s",outputDir.Data(),suffix.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErr->SetLogx();

    TH2F * histo2DRelSysErr;
    histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.23,50.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErr->GetXaxis()->SetNoExponent(kTRUE);
    //histo2DRelSysErr->GetXaxis()->SetLabelOffset(-0.01);
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
    histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.23,50.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErr->GetXaxis()->SetNoExponent(kTRUE);
    //histo2DRelStatErr->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelStatErr->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelStatErr->Draw("copy");
        TLegend* legendRelStatErr       = GetAndSetLegend2(0.24, 0.92-(0.035*nMeasSetPi0A), 0.55, 0.92, 32);
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
    histo2DRelTotErrPi0                 = new TH2F("histo2DRelTotErrPi0","histo2DRelTotErrPi0",11000,0.23,50.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrPi0, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrPi0->GetXaxis()->SetNoExponent(kTRUE);
    //histo2DRelTotErrPi0->GetXaxis()->SetLabelOffset(-0.01);
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

    histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,35.5);
    histo2DRelTotErrPi0->GetYaxis()->SetTitle("Err (%)");
    histo2DRelTotErrPi0->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionRelTotA, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombPi0InvXSectionRelTotA->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionRelStatA, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombPi0InvXSectionRelStatA->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionRelSysA, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombPi0InvXSectionRelSysA->SetLineStyle(7);
        graphCombPi0InvXSectionRelSysA->Draw("l,x0,same,e1");

        TLegend* legendRelTotErr3       = GetAndSetLegend2(0.27, 0.92-(0.035*3), 0.58, 0.92, 32);
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
    Double_t paramTCMPi0New[5]  = { graphCombPi0InvXSectionTotA->GetY()[1],0.1,
                                    graphCombPi0InvXSectionTotA->GetY()[4],0.6,3.0};
    TF1* fitTCMInvXSectionPi0   = FitObject("tcm","fitTCMInvCrossSectionPi08TeV","Pi0",graphCombPi0InvXSectionTotA,0.3,35. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);

    // Tsallis fit
    Double_t paramGraphPi0[3]                              = {5e11, 6., 0.13};
    TF1* fitInvXSectionPi0                       = FitObject("l","fitInvCrossSectionPi08TeV","Pi0",graphCombPi0InvXSectionTotA,0.3,35.,paramGraphPi0,"QNRMEX0+");

    // *************************************************************************************************************
    // Shift graphs in X direction if desired
    // *************************************************************************************************************
    if(bWCorrection.Contains("X")){
//         TF1* fitTsallisPi0PtMult                 = FitObject("tcmpt","TsallisMultWithPtPi08TeV","Pi0");

//         fitTsallisPi0PtMult->SetParameters( fitTCMInvXSectionPi0->GetParameter(0),fitTCMInvXSectionPi0->GetParameter(1), fitTCMInvXSectionPi0->GetParameter(2), fitTCMInvXSectionPi0->GetParameter(3),
//                                             fitTCMInvXSectionPi0->GetParameter(4)); // standard parameter optimize if necessary
//        TF1* fitTsallisPi0PtMult                 = FitObject("tmpt","TsallisMultWithPtPi08TeV","Pi0");
//        fitTsallisPi0PtMult->SetRange(0.3,35.);
//        fitTsallisPi0PtMult->SetParameters(fitInvXSectionPi0->GetParameter(0),fitInvXSectionPi0->GetParameter(1), fitInvXSectionPi0->GetParameter(2));

        TF1* fitTsallisPi0PtMult                 = FitObject("tmpt","TsallisMultWithPtPi08TeV","Pi0");
        fitTsallisPi0PtMult->SetRange(0.3,35.);
        fitTsallisPi0PtMult->SetParameters(fitInvXSectionPi0->GetParameter(0),fitInvXSectionPi0->GetParameter(1), fitInvXSectionPi0->GetParameter(2));

        TGraphAsymmErrors* graphCombPi0InvXSectionTotANoShift = (TGraphAsymmErrors*) graphCombPi0InvXSectionTotA->Clone("Pi0_NoShift");

        graphCombPi0InvXSectionTotA              = ApplyXshift(graphCombPi0InvXSectionTotA, fitTsallisPi0PtMult, "Pi0", kTRUE);

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

        TF1* fitTsallisPi0PtMultFromShift                 = FitObject("tmpt","TsallisMultWithPtPi08TeVFromShift","Pi0");
        fitTsallisPi0PtMultFromShift->SetRange(0.3,35.);
        fitTsallisPi0PtMultFromShift->SetParameters(fitTsallisPi0PtMult->GetParameter(0),fitTsallisPi0PtMult->GetParameter(1), fitTsallisPi0PtMult->GetParameter(2));

        TF1* fitTsallisPi0PtMultFromShiftScaled = new TF1("TsallisMultWithPtPi08TeVFromShiftScaled","(1/x)*TsallisMultWithPtPi08TeVFromShift",0.3,35.);

        //***************************************************************************************************************
        //************************************Plotting binshift corrections *********************************************
        //***************************************************************************************************************

        TCanvas* canvasShift = new TCanvas("canvasShift","",0,0,1000,900);// gives the page size
        DrawGammaCanvasSettings( canvasShift, 0.10, 0.017, 0.015, 0.1);
        canvasShift->SetLogx(1);

        Size_t textSizeSpectra          = 0.04;
        TH1F * histoBinShift = new TH1F("histoBinShift","histoBinShift",1000,0.23, 40.);
        SetStyleHistoTH1ForGraphs(histoBinShift, "#it{p}_{T} (GeV/#it{c})","bin shifted (X) / no shift",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 1.1, 1.2);
        histoBinShift->GetXaxis()->SetMoreLogLabels(1);
        histoBinShift->GetXaxis()->SetNoExponent(kTRUE);
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
        histo2DDummy2               = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.23,50.,1000,1,10e12);
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
//        DrawGammaSetMarkerTGraphAsym(graphEMCALMergedPi0InvXSectionSys, markerStyleDet[9] ,markerSizeDet[9]/2, colorDet[9], colorDet[9], widthLinesBoxes, kTRUE);
//        graphEMCALMergedPi0InvXSectionSys->Draw("pEsame");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatAUnshi, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
        graphCombPi0InvXSectionStatAUnshi->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatA, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombPi0InvXSectionStatA->Draw("pEsame");

        fitInvXSectionPi0->SetLineColor(kBlue+2);
        fitInvXSectionPi0->Draw("same");

        fitTsallisPi0PtMultFromShiftScaled->SetLineColor(kRed+2);
        fitTsallisPi0PtMultFromShiftScaled->Draw("same");

        canvasDummy2->Update();
        canvasDummy2->SaveAs(Form("%s/ComparisonShiftedPi0_8TeV.%s",outputDir.Data(),suffix.Data()));
        delete canvasDummy2;
    }

    TGraphAsymmErrors* graphCombPi0InvXSectionStatA_WOXErr = (TGraphAsymmErrors*) graphCombPi0InvXSectionStatA->Clone("graphCombPi0InvXSectionStatA_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphCombPi0InvXSectionStatA_WOXErr);
    TGraphAsymmErrors* graphPCMPi0InvXSectionStat_WOXErr = (TGraphAsymmErrors*) graphPCMPi0InvXSectionStat->Clone("graphPCMPi0InvXSectionStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphPCMPi0InvXSectionStat_WOXErr);
    TGraphAsymmErrors* graphPHOSPi0InvXSectionStat_WOXErr = (TGraphAsymmErrors*) graphPHOSPi0InvXSectionStat->Clone("graphPHOSPi0InvXSectionStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphPHOSPi0InvXSectionStat_WOXErr);
    TGraphAsymmErrors* graphEMCALPi0InvXSectionStat_WOXErr = (TGraphAsymmErrors*) graphEMCALPi0InvXSectionStat->Clone("graphEMCALPi0InvXSectionStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphEMCALPi0InvXSectionStat_WOXErr);
    TGraphAsymmErrors* graphPCMEMCALPi0InvXSectionStat_WOXErr = (TGraphAsymmErrors*) graphPCMEMCALPi0InvXSectionStat->Clone("graphPCMEMCALPi0InvXSectionStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphPCMEMCALPi0InvXSectionStat_WOXErr);
    TGraphAsymmErrors* graphEMCALMergedPi0InvXSectionStat_WOXErr = (TGraphAsymmErrors*) graphEMCALMergedPi0InvXSectionStat->Clone("graphEMCALMergedPi0InvXSectionStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphEMCALMergedPi0InvXSectionStat_WOXErr);

    // *************************************************************************************************************
    // redo fitting after binshifts
    // *************************************************************************************************************
    // Tsallis function
    graphCombPi0InvXSectionTotA->Fit(fitInvXSectionPi0,"QNRMEX0+","",0.3,35.);
    fitInvXSectionPi0           = FitObject("l","fitInvCrossSectionPi08TeV","Pi0",graphCombPi0InvXSectionTotA,0.3,35.,paramGraphPi0,"QNRMEX0+");
    //fitInvXSectionPi0           = FitObject("l","fitInvCrossSectionPi08TeV","Pi0",graphCombPi0InvXSectionTotA,0.3,35.,paramGraphPi0,"QNRMEX0+");
    cout << WriteParameterToFile(fitInvXSectionPi0)<< endl;

    //Two component model from Bylinkin
    fitTCMInvXSectionPi0        = FitObject("tcm","fitTCMInvCrossSectionPi08TeV","Pi0",graphCombPi0InvXSectionTotA,0.3,35. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitTCMInvXSectionPi0)<< endl;

    TF1* fitTCMInvXSectionPi0Plot = new TF1("twoCompModel_plotting",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1]) + [2]/(TMath::Power(1+x*x/([3]*[3]*[4]),[4]) )",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
    fitTCMInvXSectionPi0Plot->SetRange(0.3,40.);
    fitTCMInvXSectionPi0Plot->SetParameters(fitTCMInvXSectionPi0->GetParameters());
    fitTCMInvXSectionPi0Plot->SetParErrors(fitTCMInvXSectionPi0->GetParErrors());
//     TF1* fitTCMInvXSectionPi02  = FitObject("tcm","fitTCMInvCrossSectionPi08TeV2","Pi0",graphCombPi0InvXSectionTotA,0.4,40.,paramTCMPi0New,"QNRMEX0+","", kFALSE);
//     cout << WriteParameterToFile(fitTCMInvXSectionPi02)<< endl;

    Double_t paramPi0Power[3] = {1E11,0.5,6.5};
    TF1* fitPowInvXSectionPi0   = FitObject("powPure","fitPowInvXSectionPi08TeV","Pi0",graphCombPi0InvXSectionTotA,3.5,35. ,paramPi0Power,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvXSectionPi0)<< endl;

    TF1* fitPowInvXSectionPi0Stat   = FitObject("powPure","fitPowInvXSectionPi08TeVStat","Pi0",graphCombPi0InvXSectionStatA,3.5,35. ,paramPi0Power,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvXSectionPi0Stat)<< endl;

    Double_t paramPi0HageDorn[5] = {1E11,0.3,-0.1,0.5,5.95};
    TF1* fitOHagInvYieldPi0Tot   = FitObject("oHag","fitOHagInvYieldPi08TeV","Pi0",graphCombPi0InvXSectionTotA,0.3,35. ,paramPi0HageDorn,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitOHagInvYieldPi0Tot)<< endl;

    TF1* fitPCMTCMInvXSectionPi0    = FitObject("tcm","fitPCMTCMInvCrossSectionPi08TeV","Pi0",graphPCMPi0InvXSectionStat,0.3,8. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
    //cout << WriteParameterToFile(fitPCMTCMInvXSectionPi0)<< endl;

    TF1* fitPHOSTCMInvXSectionPi0    = FitObject("tcm","fitPHOSTCMInvCrossSectionPi08TeV","Pi0",graphPHOSPi0InvXSectionStat,1.0,35. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);

    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "Pi0 - Tsallis" << endl;
    fLog << WriteParameterToFile(fitInvXSectionPi0)<< endl;
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "Pi0 - TCM" << endl;
    fLog << WriteParameterToFile(fitTCMInvXSectionPi0) << endl;
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "Pi0 - Hagedorn" << endl;
    fLog << WriteParameterToFile(fitOHagInvYieldPi0Tot) << endl;
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "Pi0 - PowerLaw" << endl;
    fLog << WriteParameterToFile(fitPowInvXSectionPi0) << endl;
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "Pi0 - PowerLaw - Stat" << endl;
    fLog << WriteParameterToFile(fitPowInvXSectionPi0Stat) << endl;

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

    TGraphErrors* graphRatioPythia8ToFit             = (TGraphErrors*) graphPythia8InvXSection->Clone();
    graphRatioPythia8ToFit                           = CalculateGraphErrRatioToFit (graphRatioPythia8ToFit, fitTCMInvXSectionPi0Plot);
    TGraphErrors* graphRatioPythia8_4CToFit          = (TGraphErrors*) graphPythia8_4CInvXSection->Clone();
    graphRatioPythia8_4CToFit                        = CalculateGraphErrRatioToFit (graphRatioPythia8_4CToFit, fitTCMInvXSectionPi0Plot);
//    TH1D* histoRatioPythia8VarBinningToFit           = (TH1D*) histoPythia8InvXSection_VarBinning->Clone();
//    histoRatioPythia8VarBinningToFit                 = CalculateHistoRatioToFit (histoRatioPythia8VarBinningToFit, fitTCMInvXSectionPi0Plot);

    TGraphAsymmErrors* graphRatioPOWHEG8TeV          = (TGraphAsymmErrors*)graphPOWHEGPi08TeV->Clone();
    graphRatioPOWHEG8TeV                               = CalculateGraphErrRatioToFit (graphRatioPOWHEG8TeV, fitTCMInvXSectionPi0Plot);

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

    TGraphAsymmErrors* graphRatioPi0CombTsallisFitStatA  = (TGraphAsymmErrors*)graphCombPi0InvXSectionStatA->Clone();
    graphRatioPi0CombTsallisFitStatA                     = CalculateGraphErrRatioToFit(graphRatioPi0CombTsallisFitStatA, fitInvXSectionPi0);
    TGraphAsymmErrors* graphRatioPi0CombHagedornFitStatA = (TGraphAsymmErrors*)graphCombPi0InvXSectionStatA->Clone();
    graphRatioPi0CombHagedornFitStatA                    = CalculateGraphErrRatioToFit(graphRatioPi0CombHagedornFitStatA, fitOHagInvYieldPi0Tot);
    TGraphAsymmErrors* graphRatioPi0CombPowerFitStatA    = (TGraphAsymmErrors*)graphCombPi0InvXSectionStatA->Clone();
    graphRatioPi0CombPowerFitStatA                       = CalculateGraphErrRatioToFit(graphRatioPi0CombPowerFitStatA, fitPowInvXSectionPi0);
    TGraphAsymmErrors* graphRatioPi0CombTsallisFitSysA   = (TGraphAsymmErrors*)graphCombPi0InvXSectionSysA->Clone();
    graphRatioPi0CombTsallisFitSysA                      = CalculateGraphErrRatioToFit(graphRatioPi0CombTsallisFitSysA, fitInvXSectionPi0);
    TGraphAsymmErrors* graphRatioPi0CombHagedornFitSysA  = (TGraphAsymmErrors*)graphCombPi0InvXSectionSysA->Clone();
    graphRatioPi0CombHagedornFitSysA                     = CalculateGraphErrRatioToFit(graphRatioPi0CombHagedornFitSysA, fitOHagInvYieldPi0Tot);
    TGraphAsymmErrors* graphRatioPi0CombPowerFitSysA     = (TGraphAsymmErrors*)graphCombPi0InvXSectionSysA->Clone();
    graphRatioPi0CombPowerFitSysA                        = CalculateGraphErrRatioToFit(graphRatioPi0CombPowerFitSysA, fitPowInvXSectionPi0);

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

    TGraphAsymmErrors* graphRatioPi0CombCombFitStatA_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombCombFitStatA->Clone("graphRatioPi0CombCombFitStatA_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatA_WOXErr);

    TGraphAsymmErrors* graphRatioPi0CombTsallisFitStatA_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombTsallisFitStatA->Clone("graphRatioPi0CombTsallisFitStatA_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombTsallisFitStatA_WOXErr);
    TGraphAsymmErrors* graphRatioPi0CombHagedornFitStatA_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombHagedornFitStatA->Clone("graphRatioPi0CombHagedornFitStatA_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombHagedornFitStatA_WOXErr);
    TGraphAsymmErrors* graphRatioPi0CombPowerFitStatA_WOXErr = (TGraphAsymmErrors*) graphRatioPi0CombPowerFitStatA->Clone("graphRatioPi0CombPowerFitStatA_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombPowerFitStatA_WOXErr);

    TGraphAsymmErrors* graphRatioPi0PCMCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioPi0PCMCombFitStat->Clone("graphRatioPi0PCMCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioPi0PCMCombFitStat_WOXErr);
    TGraphAsymmErrors* graphRatioPi0PHOSCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioPi0PHOSCombFitStat->Clone("graphRatioPi0PHOSCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioPi0PHOSCombFitStat_WOXErr);
    TGraphAsymmErrors* graphRatioPi0EMCALCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioPi0EMCALCombFitStat->Clone("graphRatioPi0EMCALCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioPi0EMCALCombFitStat_WOXErr);
    TGraphAsymmErrors* graphRatioPi0PCMEMCALCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioPi0PCMEMCALCombFitStat->Clone("graphRatioPi0PCMEMCALCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioPi0PCMEMCALCombFitStat_WOXErr);
    TGraphAsymmErrors* graphRatioPi0EMCALMergedCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioPi0EMCALMergedCombFitStat->Clone("graphRatioPi0EMCALMergedCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioPi0EMCALMergedCombFitStat_WOXErr);

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
    histo2DPi0RatioToCombFit               = new TH2F("histo2DPi0RatioToCombFit","histo2DPi0RatioToCombFit",1000,0.23,50.,1000,0.2,3.0);
    SetStyleHistoTH2ForGraphs(histo2DPi0RatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/TCM fit", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
    histo2DPi0RatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DPi0RatioToCombFit->GetXaxis()->SetNoExponent(kTRUE);
//  histo2DPi0RatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
    histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DPi0RatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioPi0CombCombFitSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioPi0CombCombFitStatA_WOXErr->Draw("p,same,z");

        DrawGammaLines(0.23, 50. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.23, 50. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.23, 50. , 0.9, 0.9,0.1, kGray, 7);

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
    // *******************************************Plot different ratios to fits *********************************************
    // **********************************************************************************************************************

    histo2DPi0RatioToCombFit->SetYTitle("Data/fit");
    histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
    histo2DPi0RatioToCombFit->Draw("copy");


    while (graphRatioPi0CombPowerFitSysA->GetX()[0] < 1.6) graphRatioPi0CombPowerFitSysA->RemovePoint(0);
    while (graphRatioPi0CombPowerFitStatA_WOXErr->GetX()[0] < 1.6) graphRatioPi0CombPowerFitStatA_WOXErr->RemovePoint(0);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombPowerFitSysA, 33, markerSizeComb*2, kAzure+2 , kAzure+2, widthLinesBoxes, kTRUE);
        graphRatioPi0CombPowerFitSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombPowerFitStatA_WOXErr, 33, markerSizeComb*2, kAzure+2 , kAzure+2);
        graphRatioPi0CombPowerFitStatA_WOXErr->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombTsallisFitSysA, 25, markerSizeComb, colorTrigg[1], colorTrigg[1], widthLinesBoxes, kTRUE);
        graphRatioPi0CombTsallisFitSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombTsallisFitStatA_WOXErr, 25, markerSizeComb, colorTrigg[1] , colorTrigg[1]);
        graphRatioPi0CombTsallisFitStatA_WOXErr->Draw("p,same,z");

        graphRatioPi0CombCombFitSysA->SetMarkerColor(kGray+2);
        graphRatioPi0CombCombFitSysA->Draw("E2same");
        graphRatioPi0CombCombFitStatA_WOXErr->SetMarkerColor(kGray+2);
        graphRatioPi0CombCombFitStatA_WOXErr->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombHagedornFitSysA, 24, markerSizeComb+0.2, colorTrigg[2], colorTrigg[2], widthLinesBoxes, kTRUE);
        graphRatioPi0CombHagedornFitSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombHagedornFitStatA_WOXErr, 24, markerSizeComb+0.2, colorTrigg[2] , colorTrigg[2]);
        graphRatioPi0CombHagedornFitStatA_WOXErr->Draw("p,same,z");

        DrawGammaLines(0.23, 50. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.23, 50. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.23, 50. , 0.9, 0.9,0.1, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitPi0->Draw();

        TLegend* legendRatioPi0Fits= GetAndSetLegend2(0.12,0.95-4*1.05*textsizeLabelsPP,0.37,0.95, textSizeLabelsPixel, 1, "", 43, 0.2);
        legendRatioPi0Fits->AddEntry(graphRatioPi0CombCombFitSysA,"TCM","p");
        legendRatioPi0Fits->AddEntry(graphRatioPi0CombTsallisFitSysA,"Tsallis","p");
        legendRatioPi0Fits->AddEntry(graphRatioPi0CombHagedornFitSysA,"mod. Hagedorn","p");
        legendRatioPi0Fits->AddEntry(graphRatioPi0CombPowerFitSysA,"pure powerlaw, 3.5-35 GeV/#it{c}","p");
        legendRatioPi0Fits->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToDifferentFits_PP8TeV.%s",outputDir.Data(),suffix.Data()));
    histo2DPi0RatioToCombFit->SetYTitle("Data/TCM fit");

    // **********************************************************************************************************************
    // *******************************************Plot Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
    histo2DPi0RatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitSys, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitStat_WOXErr, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PHOSCombFitSys, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PHOSCombFitStat_WOXErr, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALCombFitSys, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALCombFitStat_WOXErr, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
//        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALMergedCombFitSys, markerStyleDet[9] ,markerSizeDet[9]*0.5, colorDet[9], colorDet[9], widthLinesBoxes, kTRUE);
//        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALMergedCombFitStat, markerStyleDet[9] ,markerSizeDet[9]*0.5, colorDet[9], colorDet[9]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMEMCALCombFitSys, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMEMCALCombFitStat_WOXErr, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4]);

        graphRatioPi0PCMCombFitSys->Draw("E2same");
        graphRatioPi0PHOSCombFitSys->Draw("E2same");
        graphRatioPi0EMCALCombFitSys->Draw("E2same");
//        graphRatioPi0EMCALMergedCombFitSys->Draw("E2same");
        graphRatioPi0PCMEMCALCombFitSys->Draw("E2same");

        graphRatioPi0PCMCombFitStat_WOXErr->Draw("p,same,z");
        graphRatioPi0PHOSCombFitStat_WOXErr->Draw("p,same,z");
        graphRatioPi0EMCALCombFitStat_WOXErr->Draw("p,same,z");
//        graphRatioPi0EMCALMergedCombFitStat->Draw("p,same,z");
        graphRatioPi0PCMEMCALCombFitStat_WOXErr->Draw("p,same,z");

        DrawGammaLines(0.23, 50. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.23, 50. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.23, 50. , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitPi0->Draw();

        //****************************** Definition of the Legend ******************************************
        //**************** Row def ************************
        Double_t rowsLegendOnlyPi0Ratio[6]          = {0.91, 0.855, 0.805, 0.755, 0.705, 0.655};
        Double_t rowsLegendOnlyPi0RatioAbs[6]       = {0.91, 2.165, 2.035, 1.905, 1.765, 1.635};
        Double_t columnsLegendOnlyPi0Ratio[3]       = {0.185, 0.385, 0.47};
        Double_t columnsLegendOnlyPi0RatioAbs[3]    = {0.115, 1.8, 2.7};
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
        TLatex *textEMCALOnlyRatioPi0               = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[3],"EMC");
        SetStyleTLatex( textEMCALOnlyRatioPi0, textSizeLabelsPixel,4);
        textEMCALOnlyRatioPi0->SetTextFont(43);
        textEMCALOnlyRatioPi0->Draw();
        TLatex *textPCMEMCALOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[4],"PCM-EMC");
        SetStyleTLatex( textPCMEMCALOnlyRatioPi0, textSizeLabelsPixel,4);
        textPCMEMCALOnlyRatioPi0->SetTextFont(43);
        textPCMEMCALOnlyRatioPi0->Draw();
//        TLatex *textEMCALMergedOnlyRatioPi0         = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[5],"EMC merged");
//        SetStyleTLatex( textEMCALMergedOnlyRatioPi0, textSizeLabelsPixel,4);
//        textEMCALMergedOnlyRatioPi0->SetTextFont(43);
//        textEMCALMergedOnlyRatioPi0->Draw();

        //****************** second Column *************************************************
        TLatex *textStatOnlyRatioPi0                = new TLatex(columnsLegendOnlyPi0Ratio[1]+0.01,rowsLegendOnlyPi0Ratio[0] ,"stat");
        SetStyleTLatex( textStatOnlyRatioPi0, textSizeLabelsPixel,4);
        textStatOnlyRatioPi0->SetTextFont(43);
        textStatOnlyRatioPi0->Draw();
        TLatex *textSysOnlyRatioPi0                 = new TLatex(columnsLegendOnlyPi0Ratio[2]+0.02 ,rowsLegendOnlyPi0Ratio[0],"syst");
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
//        TMarker* markerEMCALMergedPi0OnlyRatioPi0   = CreateMarkerFromGraph(graphRatioPi0EMCALMergedCombFitSys, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[5],1);
//        markerEMCALMergedPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[5]);

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
//        TBox* boxEMCALMergedPi0OnlyRatioPi0         = CreateBoxFromGraph(graphRatioPi0EMCALMergedCombFitSys, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[5]- heightBox,
//                                                                         columnsLegendOnlyPi0RatioAbs[2]+ 5*lengthBox, rowsLegendOnlyPi0RatioAbs[5]+ heightBox);
//        boxEMCALMergedPi0OnlyRatioPi0->Draw("l");

    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFit_PP.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // *******************************************Plot Ratio of Individual meas to Fit with MERGED **************************
    // **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DPi0RatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitSys, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitStat_WOXErr, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PHOSCombFitSys, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PHOSCombFitStat_WOXErr, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALCombFitSys, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALCombFitStat_WOXErr, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALMergedCombFitSys, markerStyleDet[9] ,markerSizeDet[9]*0.5, colorDet[9], colorDet[9], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALMergedCombFitStat_WOXErr, markerStyleDet[9] ,markerSizeDet[9]*0.5, colorDet[9], colorDet[9]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMEMCALCombFitSys, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMEMCALCombFitStat_WOXErr, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4]);

        graphRatioPi0PCMCombFitSys->Draw("E2same");
        graphRatioPi0PHOSCombFitSys->Draw("E2same");
        graphRatioPi0EMCALCombFitSys->Draw("E2same");
        graphRatioPi0EMCALMergedCombFitSys->Draw("E2same");
        graphRatioPi0PCMEMCALCombFitSys->Draw("E2same");

        graphRatioPi0PCMCombFitStat_WOXErr->Draw("p,same,z");
        graphRatioPi0PHOSCombFitStat_WOXErr->Draw("p,same,z");
        graphRatioPi0EMCALCombFitStat_WOXErr->Draw("p,same,z");
        graphRatioPi0EMCALMergedCombFitStat_WOXErr->Draw("p,same,z");
        graphRatioPi0PCMEMCALCombFitStat_WOXErr->Draw("p,same,z");

        DrawGammaLines(0.23, 50. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.23, 50. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.23, 50. , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitPi0->Draw();

        //****************** first Column **************************************************
        textPCMOnlyRatioPi0->Draw();
        textPHOSOnlyRatioPi0->Draw();
        textEMCALOnlyRatioPi0->Draw();
        textPCMEMCALOnlyRatioPi0->Draw();
        TLatex *textEMCALMergedOnlyRatioPi0         = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[5],"mEMC");
        SetStyleTLatex( textEMCALMergedOnlyRatioPi0, textSizeLabelsPixel,4);
        textEMCALMergedOnlyRatioPi0->SetTextFont(43);
        textEMCALMergedOnlyRatioPi0->Draw();

        //****************** second Column *************************************************
        textStatOnlyRatioPi0->Draw();
        textSysOnlyRatioPi0->Draw();
        markerPCMPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[1]);
        markerPHOSPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[2]);
        markerEMCALPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[3]);
        markerPCMEMCALPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[4]);
        TMarker* markerEMCALMergedPi0OnlyRatioPi0   = CreateMarkerFromGraph(graphRatioPi0EMCALMergedCombFitSys, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[5],1);
        markerEMCALMergedPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[5]);

        boxPCMPi0OnlyRatioPi0->Draw("l");
        boxPHOSPi0OnlyRatioPi0->Draw("l");
        boxEMCALPi0OnlyRatioPi0->Draw("l");
        boxPCMEMCALPi0OnlyRatioPi0->Draw("l");
        TBox* boxEMCALMergedPi0OnlyRatioPi0         = CreateBoxFromGraph(graphRatioPi0EMCALMergedCombFitSys, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[5]- heightBox,
                                                                         columnsLegendOnlyPi0RatioAbs[2]+ 5*lengthBox, rowsLegendOnlyPi0RatioAbs[5]+ heightBox);
        boxEMCALMergedPi0OnlyRatioPi0->Draw("l");

    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFitWithMerged_PP.%s",outputDir.Data(),suffix.Data()));

    // *******************************************************************************************************
    // ********************** Ratio to standalone PCM fit ****************************************************
    // *******************************************************************************************************
    TGraphAsymmErrors* graphRatioPi0PCMCombFitPCMStat      = (TGraphAsymmErrors*)graphPCMPi0InvXSectionStat->Clone();
    graphRatioPi0PCMCombFitPCMStat                         = CalculateGraphErrRatioToFit(graphRatioPi0PCMCombFitPCMStat, fitPCMTCMInvXSectionPi0);
    TGraphAsymmErrors* graphRatioPi0PCMCombFitPCMSys       = (TGraphAsymmErrors*)graphPCMPi0InvXSectionSys->Clone();
    graphRatioPi0PCMCombFitPCMSys                          = CalculateGraphErrRatioToFit(graphRatioPi0PCMCombFitPCMSys, fitPCMTCMInvXSectionPi0);

    TGraphAsymmErrors* graphRatioPi0PCMCombFitPCMStat_WOXErr = (TGraphAsymmErrors*) graphRatioPi0PCMCombFitPCMStat->Clone("graphRatioPi0PCMCombFitPCMStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioPi0PCMCombFitPCMStat_WOXErr);

    canvasRatioToCombFit->cd();
    histo2DPi0RatioToCombFit->GetXaxis()->SetRangeUser(0.23,16);
    histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.5,2.4);
    histo2DPi0RatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitSys, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitStat_WOXErr, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitPCMSys, markerStyleDet[0]+4 ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitPCMStat_WOXErr, markerStyleDet[0]+4 ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);

        graphRatioPi0PCMCombFitSys->Draw("E2same");
        graphRatioPi0PCMCombFitPCMSys->Draw("E2same");

        graphRatioPi0PCMCombFitStat_WOXErr->Draw("p,same,z");
        graphRatioPi0PCMCombFitPCMStat_WOXErr->Draw("p,same,z");

        DrawGammaLines(0.23, 16. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.23, 16. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.23, 16. , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitPi0->Draw();

        TLegend* legendPi0FitPCMStandalone        = GetAndSetLegend2(0.2, 0.92-(0.05*2), 0.45, 0.92, 38);
        legendPi0FitPCMStandalone->AddEntry(graphRatioPi0PCMCombFitSys,"PCM/comb fit");
        legendPi0FitPCMStandalone->AddEntry(graphRatioPi0PCMCombFitPCMSys,"PCM/PCM standalone fit");
        legendPi0FitPCMStandalone->Draw("");

    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfPCMToCombAndPCMStandaloneFit_PP.%s",outputDir.Data(),suffix.Data()));


    // *******************************************************************************************************
    // ********************** Ratios to TSALLIS FITS **************************************************
    // *******************************************************************************************************

    TGraphAsymmErrors* graphRatioTsallisPi0CombCombFitTotA     = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotA->Clone();
    graphRatioTsallisPi0CombCombFitTotA                        = CalculateGraphErrRatioToFit(graphRatioTsallisPi0CombCombFitTotA, fitInvXSectionPi0);
    TGraphAsymmErrors* graphRatioTsallisPi0CombCombFitStatA    = (TGraphAsymmErrors*)graphCombPi0InvXSectionStatA->Clone();
    graphRatioTsallisPi0CombCombFitStatA                       = CalculateGraphErrRatioToFit(graphRatioTsallisPi0CombCombFitStatA, fitInvXSectionPi0);
    TGraphAsymmErrors* graphRatioTsallisPi0CombCombFitSysA     = (TGraphAsymmErrors*)graphCombPi0InvXSectionSysA->Clone();
    graphRatioTsallisPi0CombCombFitSysA                        = CalculateGraphErrRatioToFit(graphRatioTsallisPi0CombCombFitSysA, fitInvXSectionPi0);
    TGraphAsymmErrors* graphRatioTsallisPi0PCMCombFitStat      = (TGraphAsymmErrors*)graphPCMPi0InvXSectionStat->Clone();
    graphRatioTsallisPi0PCMCombFitStat                         = CalculateGraphErrRatioToFit(graphRatioTsallisPi0PCMCombFitStat, fitInvXSectionPi0);
    TGraphAsymmErrors* graphRatioTsallisPi0PCMCombFitSys       = (TGraphAsymmErrors*)graphPCMPi0InvXSectionSys->Clone();
    graphRatioTsallisPi0PCMCombFitSys                          = CalculateGraphErrRatioToFit(graphRatioTsallisPi0PCMCombFitSys, fitInvXSectionPi0);
    TGraphAsymmErrors* graphRatioTsallisPi0PHOSCombFitStat     = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionStat->Clone();
    graphRatioTsallisPi0PHOSCombFitStat                        = CalculateGraphErrRatioToFit(graphRatioTsallisPi0PHOSCombFitStat, fitInvXSectionPi0);
    TGraphAsymmErrors* graphRatioTsallisPi0PHOSCombFitSys      = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionSys->Clone();
    graphRatioTsallisPi0PHOSCombFitSys                         = CalculateGraphErrRatioToFit(graphRatioTsallisPi0PHOSCombFitSys, fitInvXSectionPi0);
    TGraphAsymmErrors* graphRatioTsallisPi0EMCALCombFitStat    = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionStat->Clone();
    graphRatioTsallisPi0EMCALCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioTsallisPi0EMCALCombFitStat, fitInvXSectionPi0);
    TGraphAsymmErrors* graphRatioTsallisPi0EMCALCombFitSys     = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionSys->Clone();
    graphRatioTsallisPi0EMCALCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioTsallisPi0EMCALCombFitSys, fitInvXSectionPi0);
    TGraphAsymmErrors* graphRatioTsallisPi0PCMEMCALCombFitStat = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionStat->Clone();
    graphRatioTsallisPi0PCMEMCALCombFitStat                    = CalculateGraphErrRatioToFit(graphRatioTsallisPi0PCMEMCALCombFitStat, fitInvXSectionPi0);
    TGraphAsymmErrors* graphRatioTsallisPi0PCMEMCALCombFitSys  = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionSys->Clone();
    graphRatioTsallisPi0PCMEMCALCombFitSys                     = CalculateGraphErrRatioToFit(graphRatioTsallisPi0PCMEMCALCombFitSys, fitInvXSectionPi0);
    TGraphAsymmErrors* graphRatioTsallisPi0EMCALMergedCombFitStat  = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionStat->Clone();
    graphRatioTsallisPi0EMCALMergedCombFitStat                     = CalculateGraphErrRatioToFit(graphRatioTsallisPi0EMCALMergedCombFitStat, fitInvXSectionPi0);
    TGraphAsymmErrors* graphRatioTsallisPi0EMCALMergedCombFitSys   = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionSys->Clone();
    graphRatioTsallisPi0EMCALMergedCombFitSys                      = CalculateGraphErrRatioToFit(graphRatioTsallisPi0EMCALMergedCombFitSys, fitInvXSectionPi0);

    TGraphAsymmErrors* graphRatioTsallisPi0CombCombFitStatA_WOXErr = (TGraphAsymmErrors*) graphRatioTsallisPi0CombCombFitStatA->Clone("graphRatioTsallisPi0CombCombFitStatA_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioTsallisPi0CombCombFitStatA_WOXErr);
    TGraphAsymmErrors* graphRatioTsallisPi0PCMCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioTsallisPi0PCMCombFitStat->Clone("graphRatioTsallisPi0PCMCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioTsallisPi0PCMCombFitStat_WOXErr);
    TGraphAsymmErrors* graphRatioTsallisPi0PHOSCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioTsallisPi0PHOSCombFitStat->Clone("graphRatioTsallisPi0PHOSCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioTsallisPi0PHOSCombFitStat_WOXErr);
    TGraphAsymmErrors* graphRatioTsallisPi0EMCALCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioTsallisPi0EMCALCombFitStat->Clone("graphRatioTsallisPi0EMCALCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioTsallisPi0EMCALCombFitStat_WOXErr);
    TGraphAsymmErrors* graphRatioTsallisPi0PCMEMCALCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioTsallisPi0PCMEMCALCombFitStat->Clone("graphRatioTsallisPi0PCMEMCALCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioTsallisPi0PCMEMCALCombFitStat_WOXErr);
    TGraphAsymmErrors* graphRatioTsallisPi0EMCALMergedCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioTsallisPi0EMCALMergedCombFitStat->Clone("graphRatioTsallisPi0EMCALMergedCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioTsallisPi0EMCALMergedCombFitStat_WOXErr);

    // **********************************************************************************************************************
    // ******************************************* Plot Ratio of Comb to Fit ****************************************
    // **********************************************************************************************************************

    canvasRatioToCombFit->SetLogx();
    histo2DPi0RatioToCombFit->GetXaxis()->SetRangeUser(0.23,50);
    histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.6,2.1);
    histo2DPi0RatioToCombFit->SetYTitle("Data/Tsallis fit");
    histo2DPi0RatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioTsallisPi0CombCombFitSysA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioTsallisPi0CombCombFitSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioTsallisPi0CombCombFitStatA_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioTsallisPi0CombCombFitStatA_WOXErr->Draw("p,same,z");

        DrawGammaLines(0.23, 50. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.23, 50. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.23, 50. , 0.9, 0.9,0.1, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitPi0->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToCombTsallisFit_PP8TeV.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // *******************************************Plot Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
    histo2DPi0RatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioTsallisPi0PCMCombFitSys, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioTsallisPi0PCMCombFitStat_WOXErr, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
        DrawGammaSetMarkerTGraphAsym(graphRatioTsallisPi0PHOSCombFitSys, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioTsallisPi0PHOSCombFitStat_WOXErr, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);
        DrawGammaSetMarkerTGraphAsym(graphRatioTsallisPi0EMCALCombFitSys, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioTsallisPi0EMCALCombFitStat_WOXErr, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
        DrawGammaSetMarkerTGraphAsym(graphRatioTsallisPi0PCMEMCALCombFitSys, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioTsallisPi0PCMEMCALCombFitStat_WOXErr, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4]);

        graphRatioTsallisPi0PCMCombFitSys->Draw("E2same");
        graphRatioTsallisPi0PHOSCombFitSys->Draw("E2same");
        graphRatioTsallisPi0EMCALCombFitSys->Draw("E2same");
        graphRatioTsallisPi0PCMEMCALCombFitSys->Draw("E2same");

        graphRatioTsallisPi0PCMCombFitStat_WOXErr->Draw("p,same,z");
        graphRatioTsallisPi0PHOSCombFitStat_WOXErr->Draw("p,same,z");
        graphRatioTsallisPi0EMCALCombFitStat_WOXErr->Draw("p,same,z");
        graphRatioTsallisPi0PCMEMCALCombFitStat_WOXErr->Draw("p,same,z");

        DrawGammaLines(0.23, 50. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.23, 50. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.23, 50. , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitPi0->Draw();

        //****************** first Column **************************************************
        textPCMOnlyRatioPi0->Draw();
        textPHOSOnlyRatioPi0->Draw();
        textEMCALOnlyRatioPi0->Draw();
        textPCMEMCALOnlyRatioPi0->Draw();

        //****************** second Column *************************************************
        textStatOnlyRatioPi0->Draw();
        textSysOnlyRatioPi0->Draw();
        markerPCMPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[1]);
        markerPHOSPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[2]);
        markerEMCALPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[3]);
        markerPCMEMCALPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[4]);

        boxPCMPi0OnlyRatioPi0->Draw("l");
        boxPHOSPi0OnlyRatioPi0->Draw("l");
        boxEMCALPi0OnlyRatioPi0->Draw("l");
        boxPCMEMCALPi0OnlyRatioPi0->Draw("l");

    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombTsallisFit_PP.%s",outputDir.Data(),suffix.Data()));
    histo2DPi0RatioToCombFit->SetYTitle("Data/TCM fit");

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
        legendPi0FitPHOSStandaloneA->AddEntry(graphRatioPHOSprelimFinalMBsys,"PHOS MB prelim / PHOS final");
        legendPi0FitPHOSStandaloneA->AddEntry(graphRatioPHOSprelimFinalPHOSsys,"PHOS PHOS prelim / PHOS final");
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

    // *******************************************************************************************************
    // ************************** Combination of different eta measurements **********************************
    // *******************************************************************************************************
    // REMARKS:
    //      - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //      - currently only PCM-EMCAL vs others fully implemeted energy independent
    //      - extendable to other energies
    //     - offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented

    // definition of array of histograms (NULL - means we have no measurement at this energy for this rec-method)
    // for statistical error and final value from respective method
    TH1D* statErrorCollectionEta[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorCollectionEta[i]   = NULL;
    }
    statErrorCollectionEta[0]       = (TH1D*)histoPCMEtaInvXSectionStat->Clone("statErrPCMEta");
    statErrorCollectionEta[2]       = (TH1D*)histoEMCALEtaInvXSectionStat->Clone("statErrEMCALEta");
    statErrorCollectionEta[4]       = (TH1D*)histoPCMEMCALEtaInvXSectionStat->Clone("statErrPCMEMCALEta");

    // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
    // for systematic error from respective method
    TGraphAsymmErrors* sysErrorCollectionEta[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorCollectionEta[i]    = NULL;
    }
    sysErrorCollectionEta[0]        = (TGraphAsymmErrors*)graphPCMEtaInvXSectionSys->Clone("sysErrPCMEta");
    sysErrorCollectionEta[2]        = (TGraphAsymmErrors*)graphEMCALEtaInvXSectionSys->Clone("sysErrEMCALEta");
    sysErrorCollectionEta[4]        = (TGraphAsymmErrors*)graphPCMEMCALEtaInvXSectionSys->Clone("sysErrPCMEMCALEta");

    TH1D* statErrorRelCollectionEta[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorRelCollectionEta[i]        = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (statErrorCollectionEta[i])
            statErrorRelCollectionEta[i]    = CalculateRelErrUpTH1D( statErrorCollectionEta[i], Form("relativeStatErrorEta_%s", nameMeasGlobal[i].Data()));
    }

    TGraphAsymmErrors* sysErrorRelCollectionEta[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorRelCollectionEta[i]         = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (sysErrorCollectionEta[i])
            sysErrorRelCollectionEta[i]     = CalculateRelErrUpAsymmGraph( sysErrorCollectionEta[i], Form("relativeSysErrorEta_%s", nameMeasGlobal[i].Data()));
    }


    // Definition of binning for eta meson, take care that it is the correct one
    Double_t xPtLimitsEtaWOMerged[70];
    Int_t maxNBinsEtaW0Merged       = GetBinning( xPtLimitsEtaWOMerged, "Eta", "8TeV", 2 );
    maxNBinsEtaW0Merged--;
    for (Int_t i = 0; i< maxNBinsEtaW0Merged; i++){
        cout << i << ": "<< xPtLimitsEtaWOMerged[i] <<" - " << xPtLimitsEtaWOMerged[i+1]<< ", " <<endl;
    }

//     return;
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Eta
    Int_t offSetsEta[11]            = { -1,  0,  0,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t offSetsEtaSys[11]         = { 0,  0,  5,  0,  2,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t offSetEtaShifting[11]     = { 0,  0,  0,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t nComBinsEtaShifting[11]   = { 0,  0,  0,  0,  0,
                                        0,  0,  0,  0,  0,
                                        0};

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************

    TGraph* graphWeightsEtaA[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsEtaA[i]         = NULL;
    }

    // Declaration & calculation of combined spectrum
    TString fileNameEtaOutputWeightingA                  = Form("%s/Eta_WeightingMethodA.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombEtaInvXSectionStatA      = NULL;
    TGraphAsymmErrors* graphCombEtaInvXSectionSysA       = NULL;
    TGraphAsymmErrors* graphCombEtaInvXSectionTotA       = CombinePtPointsSpectraFullCorrMat(   statErrorCollectionEta,    sysErrorCollectionEta,
                                                                                                xPtLimitsEtaWOMerged, maxNBinsEtaW0Merged,
//                                                                                                      xPtLimitsPi0, 33,
                                                                                                offSetsEta, offSetsEtaSys,
                                                                                                graphCombEtaInvXSectionStatA, graphCombEtaInvXSectionSysA,
                                                                                                fileNameEtaOutputWeightingA,"8TeV", "Eta", kFALSE,
                                                                                                0x0, fileInputCorrFactors
                                                                                            );
    graphCombEtaInvXSectionStatA->Print();
    if (graphCombEtaInvXSectionStatA == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }

     //return;

    // Reading weights from output file for plotting
    ifstream fileWeightsEtaReadA;
    fileWeightsEtaReadA.open(fileNameEtaOutputWeightingA,ios_base::in);
    cout << "reading" << fileNameEtaOutputWeightingA << endl;
    Double_t xValuesEtaReadA[70];
    Double_t weightsEtaReadA[11][70];
    Int_t availableEtaMeasA[11]        = { -1, -1, -1, -1, -1,
                                        -1, -1, -1, -1, -1,
                                        -1};
    Int_t nMeasSetEtaA                 = 3;
    Int_t nPtBinsEtaReadA              = 0;
    while(!fileWeightsEtaReadA.eof() && nPtBinsEtaReadA < 70){
        TString garbage             = "";
        if (nPtBinsEtaReadA == 0){
            fileWeightsEtaReadA >> garbage ;//>> availableEtaMeas[0] >> availableEtaMeas[1] >> availableEtaMeas[2] >> availableEtaMeas[3];
            for (Int_t i = 0; i < nMeasSetEtaA; i++){
                fileWeightsEtaReadA >> availableEtaMeasA[i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availableEtaMeasA[i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsEtaReadA >> xValuesEtaReadA[nPtBinsEtaReadA-1];
            for (Int_t i = 0; i < nMeasSetEtaA; i++){
                fileWeightsEtaReadA >> weightsEtaReadA[availableEtaMeasA[i]][nPtBinsEtaReadA-1] ;
            }
            cout << "read: "<<  nPtBinsEtaReadA << "\t"<< xValuesEtaReadA[nPtBinsEtaReadA-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetEtaA; i++){
                cout << weightsEtaReadA[availableEtaMeasA[i]][nPtBinsEtaReadA-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsEtaReadA++;
    }
    nPtBinsEtaReadA                    = nPtBinsEtaReadA-2 ;
    fileWeightsEtaReadA.close();

    for (Int_t i = 0; i < nMeasSetEtaA; i++){
        graphWeightsEtaA[availableEtaMeasA[i]]                        = new TGraph(nPtBinsEtaReadA,xValuesEtaReadA,weightsEtaReadA[availableEtaMeasA[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsEtaReadA; n++){
            if (graphWeightsEtaA[availableEtaMeasA[i]]->GetY()[bin] == 0) graphWeightsEtaA[availableEtaMeasA[i]]->RemovePoint(bin);
            else bin++;
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Plotting weights Method A ************************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel           = 900*0.04;

    canvasWeights->cd();

    TH2F * histo2DEtaWeights;
    histo2DEtaWeights = new TH2F("histo2DEtaWeights","histo2DEtaWeights",11000,0.33,50.,1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DEtaWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DEtaWeights->GetXaxis()->SetMoreLogLabels();
    histo2DEtaWeights->GetXaxis()->SetNoExponent(kTRUE);
    //histo2DEtaWeights->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DEtaWeights->GetYaxis()->SetRangeUser(-10,10);
    histo2DEtaWeights->Draw("copy");

        TLegend* legendWeightsEta   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEtaA), 32);
        for (Int_t i = 0; i < nMeasSetEtaA; i++){
            DrawGammaSetMarkerTGraph(graphWeightsEtaA[availableEtaMeasA[i]], markerStyleDet[availableEtaMeasA[i]], markerSizeDet[availableEtaMeasA[i]]*0.5, colorDet[availableEtaMeasA[i]] , colorDet[availableEtaMeasA[i]]);
            graphWeightsEtaA[availableEtaMeasA[i]]->Draw("p,same,z");
            legendWeightsEta->AddEntry(graphWeightsEtaA[availableEtaMeasA[i]],nameMeasGlobal[availableEtaMeasA[i]],"p");
        }
        legendWeightsEta->Draw();

        labelWeightsEnergy->Draw();
        TLatex *labelWeightsEta         = new TLatex(0.7,0.16,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelWeightsEta, 0.85*textSizeLabelsPixel,4);
        labelWeightsEta->SetTextFont(43);
        labelWeightsEta->Draw();

//      DrawGammaLines(0.33, 25. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.33, 50. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.33, 50. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.33, 50. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.33, 50. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/Eta_WeightsA.%s",outputDir.Data(),suffix.Data()));


    // *********************************************************************************************************************
    // ************************************ Visualize relative errors Eta ******************************************************
    // *********************************************************************************************************************

    canvasRelSysErr->cd();

    TH2F * histo2DRelSysErrEta;
    histo2DRelSysErrEta                 = new TH2F("histo2DRelSysErrEta","histo2DRelSysErrEta",11000,0.33,50.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErrEta, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErrEta->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErrEta->GetXaxis()->SetNoExponent(kTRUE);
//    histo2DRelSysErrEta->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelSysErrEta->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelSysErrEta->Draw("copy");

        TLegend* legendRelSysErrEta        = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetEtaA), 0.95, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetEtaA; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionEta[availableEtaMeasA[i]], markerStyleDet[availableEtaMeasA[i]], markerSizeDet[availableEtaMeasA[i]]*0.5, colorDet[availableEtaMeasA[i]],
                                     colorDet[availableEtaMeasA[i]]);
            sysErrorRelCollectionEta[availableEtaMeasA[i]]->Draw("p,same,z");
            legendRelSysErrEta->AddEntry(sysErrorRelCollectionEta[availableEtaMeasA[i]],nameMeasGlobal[availableEtaMeasA[i]],"p");
        }
        legendRelSysErrEta->Draw();

        labelRelSysErrEnergy->Draw();
        TLatex *labelRelSysErrEta       = new TLatex(0.15,0.85,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelSysErrEta, 0.85*textSizeLabelsPixel,4);
        labelRelSysErrEta->SetTextFont(43);
        labelRelSysErrEta->Draw();

    canvasRelSysErr->SaveAs(Form("%s/Eta_RelSysErr.%s",outputDir.Data(),suffix.Data()));

    histo2DRelSysErrEta->GetYaxis()->SetRangeUser(5.,30.5);
    histo2DRelSysErrEta->Draw("copy");

        for (Int_t i = 0; i < nMeasSetEtaA; i++){
            sysErrorRelCollectionEta[availableEtaMeasA[i]]->Draw("p,same,z");
        }
        legendRelSysErrEta->Draw();

        labelRelSysErrEnergy->Draw();
        labelRelSysErrEta->Draw();

    canvasRelSysErr->SaveAs(Form("%s/Eta_RelSysErrZoomed.%s",outputDir.Data(),suffix.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors Eta **************************************************
    //  *********************************************************************************************************************

    canvasRelStatErr->cd();

    TH2F * histo2DRelStatErrEta;
    histo2DRelStatErrEta                = new TH2F("histo2DRelStatErrEta","histo2DRelStatErrEta",11000,0.33,50.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErrEta, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErrEta->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErrEta->GetXaxis()->SetNoExponent(kTRUE);
//    histo2DRelStatErrEta->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRelStatErrEta->GetYaxis()->SetRangeUser(0,60.5);
    histo2DRelStatErrEta->Draw("copy");

        TLegend* legendRelStatErrEta       = GetAndSetLegend2(0.24, 0.92-(0.035*nMeasSetEtaA), 0.55, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetEtaA; i++){
            DrawGammaSetMarker(statErrorRelCollectionEta[availableEtaMeasA[i]], markerStyleDet[availableEtaMeasA[i]], markerSizeDet[availableEtaMeasA[i]]*0.5, colorDet[availableEtaMeasA[i]] ,
                               colorDet[availableEtaMeasA[i]]);
            statErrorRelCollectionEta[availableEtaMeasA[i]]->Draw("p,same,z");
            legendRelStatErrEta->AddEntry(statErrorRelCollectionEta[availableEtaMeasA[i]],nameMeasGlobal[availableEtaMeasA[i]],"p");
        }
        legendRelStatErrEta->Draw();

        labelRelStatErrEnergy->Draw();
        TLatex *labelRelStatErrEta      = new TLatex(0.75,0.85,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelStatErrEta, 0.85*textSizeLabelsPixel,4);
        labelRelStatErrEta->SetTextFont(43);
        labelRelStatErrEta->Draw();

    canvasRelStatErr->SaveAs(Form("%s/Eta_RelStatErr.%s",outputDir.Data(),suffix.Data()));

    histo2DRelStatErrEta->GetYaxis()->SetRangeUser(0,60.5);
    histo2DRelStatErrEta->Draw("copy");
        for (Int_t i = 0; i < nMeasSetEtaA; i++){
            statErrorRelCollectionEta[availableEtaMeasA[i]]->Draw("p,same,z");
        }
        legendRelStatErrEta->Draw();

        labelRelStatErrEnergy->Draw();
        labelRelStatErrEta->Draw();

    canvasRelStatErr->SaveAs(Form("%s/Eta_RelStatErrZoomed.%s",outputDir.Data(),suffix.Data()));


    //  *********************************************************************************************************************
    //  ************************ Visualize relative total errors of different combination methods Eta ***********************
    //  *********************************************************************************************************************
    TGraphAsymmErrors* graphCombEtaInvXSectionRelStatA       = CalculateRelErrUpAsymmGraph( graphCombEtaInvXSectionStatA, "relativeStatErrorEta_MethodA");
    TGraphAsymmErrors* graphCombEtaInvXSectionRelSysA        = CalculateRelErrUpAsymmGraph( graphCombEtaInvXSectionSysA, "relativeSysErrorEta_MethodA");
    TGraphAsymmErrors* graphCombEtaInvXSectionRelTotA        = CalculateRelErrUpAsymmGraph( graphCombEtaInvXSectionTotA, "relativeTotalErrorEta_MethodA");


    canvasRelTotErr->cd();
    TH2F * histo2DRelTotErrEta;
    histo2DRelTotErrEta                 = new TH2F("histo2DRelTotErrEta","histo2DRelTotErrEta",11000,0.33,50.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrEta, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrEta->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrEta->GetXaxis()->SetNoExponent(kTRUE);
//    histo2DRelTotErrEta->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRelTotErrEta->Draw("copy");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionRelTotA, markerStyleComb+4, markerSizeComb, kBlue+2 , kBlue+2);
        graphCombEtaInvXSectionRelTotA->Draw("p,same,z");

        TLegend* legendRelTotErrEta     = GetAndSetLegend2(0.24, 0.92-(0.035*1), 0.55, 0.92, 32);
        legendRelTotErrEta->AddEntry(graphCombEtaInvXSectionRelTotA,"All","p");
        legendRelTotErrEta->Draw();

        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrEta       = new TLatex(0.75,0.85,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelTotErrEta, 0.85*textSizeLabelsPixel,4);
        labelRelTotErrEta->SetTextFont(43);
        labelRelTotErrEta->Draw();

    canvasRelTotErr->SaveAs(Form("%s/Eta_RelTotErr.%s",outputDir.Data(),suffix.Data()));

    histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,65.5);
    histo2DRelTotErrEta->Draw("copy");
        graphCombEtaInvXSectionRelTotA->Draw("p,same,z");

        legendRelTotErrEta->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrEta->Draw();

    canvasRelTotErr->SaveAs(Form("%s/Eta_RelTotErrZoomed.%s",outputDir.Data(),suffix.Data()));

    histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,65.5);
    histo2DRelTotErrEta->GetYaxis()->SetTitle("Err (%)");
    histo2DRelTotErrEta->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionRelTotA, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombEtaInvXSectionRelTotA->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionRelStatA, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombEtaInvXSectionRelStatA->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionRelSysA, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombEtaInvXSectionRelSysA->SetLineStyle(7);
        graphCombEtaInvXSectionRelSysA->Draw("l,x0,same,e1");

        legendRelTotErr3->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrEta->Draw();

    canvasRelTotErr->SaveAs(Form("%s/Eta_RelMethodAdecomp.%s",outputDir.Data(),suffix.Data()));


    // **********************************************************************************************************************
    // ************************************* Calculating bin shifted spectra & fitting **************************************
    // **********************************************************************************************************************

    // Cloning spectra
    TGraphAsymmErrors* graphCombEtaInvXSectionTotAUnshi      = (TGraphAsymmErrors*)graphCombEtaInvXSectionTotA->Clone("EtaUnshifted");
    TGraphAsymmErrors* graphCombEtaInvXSectionStatAUnshi     = (TGraphAsymmErrors*)graphCombEtaInvXSectionStatA->Clone("EtaUnshiftedStat");
    TGraphAsymmErrors* graphCombEtaInvXSectionSysAUnshi      = (TGraphAsymmErrors*)graphCombEtaInvXSectionSysA->Clone("EtaUnshiftedSys");

    TGraphAsymmErrors* graphPCMEtaInvXSectionStatUnshi       = (TGraphAsymmErrors*)graphPCMEtaInvXSectionStat->Clone("EtaUnshiftedStatPCM");
    TGraphAsymmErrors* graphPCMEtaInvXSectionSysUnshi        = (TGraphAsymmErrors*)graphPCMEtaInvXSectionSys->Clone("EtaUnshiftedSysPCM");

    TGraphAsymmErrors* graphEMCALEtaInvXSectionStatUnshi     = (TGraphAsymmErrors*)graphEMCALEtaInvXSectionStat->Clone("EtaUnshiftedStatEMCAL");
    TGraphAsymmErrors* graphEMCALEtaInvXSectionSysUnshi      = (TGraphAsymmErrors*)graphEMCALEtaInvXSectionSys->Clone("EtaUnshiftedSysEMCAL");

    TGraphAsymmErrors* graphPCMEMCALEtaInvXSectionStatUnshi  = (TGraphAsymmErrors*)graphPCMEMCALEtaInvXSectionStat->Clone("EtaUnshiftedStatPCMEMCAL");
    TGraphAsymmErrors* graphPCMEMCALEtaInvXSectionSysUnshi   = (TGraphAsymmErrors*)graphPCMEMCALEtaInvXSectionSys->Clone("EtaUnshiftedSysPCMEMCAL");

    // fitting spectrum with intial parameters
    Double_t paramGraphEta[3]                    = {1.5e9, 6.6, 0.22};
    TF1* fitInvXSectionEta                       = FitObject("l","fitInvCrossSectionEta8TeV","Eta",graphCombEtaInvXSectionTotAUnshi,0.4,35.,paramGraphEta,"QNRMEX0+");

    // *************************************************************************************************************
    // Shift graphs in X direction if desired
    // *************************************************************************************************************
    if(bWCorrection.Contains("X") ){
//        TF1* fitTsallisEtaPtMult                 = FitObject("tmpt","TsallisMultWithPtEta8TeV","Eta");
//        fitTsallisEtaPtMult->SetParameters(paramGraphEta[0],paramGraphEta[1], paramGraphEta[2]) ; // standard parameter optimize if necessary
//         fitTsallisEtaPtMult->SetParameters(paramTCMEta[0],paramTCMEta[1], paramTCMEta[2], paramTCMEta[3], paramTCMEta[4]) ; // standard parameter optimize if necessary

        TF1* fitTsallisEtaPtMult                 = FitObject("tmpt","TsallisMultWithPtEta8TeV","Eta");
        fitTsallisEtaPtMult->SetRange(0.3,35.);
        fitTsallisEtaPtMult->SetParameters(fitInvXSectionEta->GetParameter(0),fitInvXSectionEta->GetParameter(1), fitInvXSectionEta->GetParameter(2));

        TGraphAsymmErrors* graphCombEtaInvXSectionTotANoShift = (TGraphAsymmErrors*) graphCombEtaInvXSectionTotA->Clone("Eta_NoShift");

        graphCombEtaInvXSectionTotA              = ApplyXshift(graphCombEtaInvXSectionTotA, fitTsallisEtaPtMult, "Eta", kTRUE);

        graphCombEtaInvXSectionStatA             = ApplyXshiftIndividualSpectra (   graphCombEtaInvXSectionTotA,
                                                                                    graphCombEtaInvXSectionStatA,
                                                                                    fitTsallisEtaPtMult,
                                                                                    0, graphCombEtaInvXSectionStatA->GetN(),"Eta");
        graphCombEtaInvXSectionSysA              = ApplyXshiftIndividualSpectra (   graphCombEtaInvXSectionTotA,
                                                                                    graphCombEtaInvXSectionSysA,
                                                                                    fitTsallisEtaPtMult,
                                                                                    0, graphCombEtaInvXSectionSysA->GetN(), "Eta");
        graphPCMEtaInvXSectionStat               = ApplyXshiftIndividualSpectra(    graphCombEtaInvXSectionTotA,
                                                                                    graphPCMEtaInvXSectionStat,
                                                                                    fitTsallisEtaPtMult,
                                                                                    offSetEtaShifting[0], nComBinsEtaShifting[0], "Eta");
        graphPCMEtaInvXSectionSys                = ApplyXshiftIndividualSpectra(    graphCombEtaInvXSectionTotA,
                                                                                    graphPCMEtaInvXSectionSys,
                                                                                    fitTsallisEtaPtMult,
                                                                                    offSetEtaShifting[0], nComBinsEtaShifting[0], "Eta");
        graphEMCALEtaInvXSectionStat             = ApplyXshiftIndividualSpectra(    graphCombEtaInvXSectionTotA,
                                                                                    graphEMCALEtaInvXSectionStat,
                                                                                    fitTsallisEtaPtMult,
                                                                                    offSetEtaShifting[2], nComBinsEtaShifting[2], "Eta");
        graphEMCALEtaInvXSectionSys              = ApplyXshiftIndividualSpectra(    graphCombEtaInvXSectionTotA,
                                                                                    graphEMCALEtaInvXSectionSys,
                                                                                    fitTsallisEtaPtMult,
                                                                                    offSetEtaShifting[2], nComBinsEtaShifting[2], "Eta");
        graphPCMEMCALEtaInvXSectionStat          = ApplyXshiftIndividualSpectra(    graphCombEtaInvXSectionTotA,
                                                                                    graphPCMEMCALEtaInvXSectionStat,
                                                                                    fitTsallisEtaPtMult,
                                                                                    offSetEtaShifting[4], nComBinsEtaShifting[4], "Eta");
        graphPCMEMCALEtaInvXSectionSys           = ApplyXshiftIndividualSpectra(    graphCombEtaInvXSectionTotA,
                                                                                    graphPCMEMCALEtaInvXSectionSys,
                                                                                    fitTsallisEtaPtMult,
                                                                                    offSetEtaShifting[4], nComBinsEtaShifting[4], "Eta");

        TF1* fitTsallisEtaPtMultFromShift                 = FitObject("tmpt","TsallisMultWithPtEta8TeVFromShift","Eta");
        fitTsallisEtaPtMultFromShift->SetRange(0.4,35.);
        fitTsallisEtaPtMultFromShift->SetParameters(fitTsallisEtaPtMult->GetParameter(0),fitTsallisEtaPtMult->GetParameter(1), fitTsallisEtaPtMult->GetParameter(2));

        TF1* fitTsallisEtaPtMultFromShiftScaled = new TF1("TsallisMultWithPtEta8TeVFromShiftScaled","(1/x)*TsallisMultWithPtEta8TeVFromShift",0.4,35.);

        //***************************************************************************************************************
        //************************************Plotting binshift corrections *********************************************
        //***************************************************************************************************************

        TCanvas* canvasShift = new TCanvas("canvasShift","",0,0,1000,900);// gives the page size
        DrawGammaCanvasSettings( canvasShift, 0.10, 0.017, 0.015, 0.1);
        canvasShift->SetLogx(1);

        Size_t textSizeSpectra          = 0.04;
        TH1F * histoBinShift = new TH1F("histoBinShift","histoBinShift",1000,0.33, 40.);
        SetStyleHistoTH1ForGraphs(histoBinShift, "#it{p}_{T} (GeV/#it{c})","bin shifted (X) / no shift",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 1.1, 1.2);
        histoBinShift->GetXaxis()->SetMoreLogLabels(1);
        histoBinShift->GetXaxis()->SetNoExponent(kTRUE);
        histoBinShift->GetYaxis()->SetRangeUser(0.95,1.05);
        histoBinShift->DrawCopy();

        Int_t numberPoints   = graphCombEtaInvXSectionTotANoShift->GetN();
        Double_t *xPoint     = graphCombEtaInvXSectionTotANoShift->GetX();
        Double_t* xvalueErrUp  = graphCombEtaInvXSectionTotANoShift->GetEXhigh();
        Double_t* xvalueErrLow = graphCombEtaInvXSectionTotANoShift->GetEXlow();
        Double_t *xPointShift= graphCombEtaInvXSectionTotA->GetX();
        for (Int_t i=0; i<numberPoints; i++) {
          graphCombEtaInvXSectionTotANoShift->SetPoint(i,xPoint[i],xPointShift[i]/xPoint[i]);
          graphCombEtaInvXSectionTotANoShift->SetPointError(i,xvalueErrLow[i],xvalueErrUp[i],0,0);
        }
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionTotANoShift, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombEtaInvXSectionTotANoShift->Draw("p same");

        TLatex *labelRatioToFitBinShift   = new TLatex(0.72, 0.91, collisionSystem8TeV.Data());
        SetStyleTLatex( labelRatioToFitBinShift, textSizeLabelsPixel,4);
        labelRatioToFitBinShift->SetTextFont(43);
        labelRatioToFitBinShift->Draw();
        TLatex *labelRatioToFitALICEBinShift    = new TLatex(0.852, 0.86, "ALICE");
        SetStyleTLatex( labelRatioToFitALICEBinShift, textSizeLabelsPixel,4);
        labelRatioToFitALICEBinShift->SetTextFont(43);
        labelRatioToFitALICEBinShift->Draw();
        TLatex *labelRatioToFitPi0BinShift      = new TLatex(0.826, 0.807, "#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitPi0BinShift, textSizeLabelsPixel,4);
        labelRatioToFitPi0BinShift->SetTextFont(43);
        labelRatioToFitPi0BinShift->Draw();

        canvasShift->Update();
        canvasShift->SaveAs(Form("%s/BinShiftCorrection_Eta.%s",outputDir.Data(),suffix.Data()));
        canvasShift->SetLogx(0);

        // *************************************************************************************************************
        // Plot control graphs
        // *************************************************************************************************************

        TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
        DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.08);
        canvasDummy2->SetLogy();
        canvasDummy2->SetLogx();
        TH2F* histo2DDummy3;
        histo2DDummy3               = new TH2F("histo2DDummy3","histo2DDummy3",1000,0.33,50.,1000,1,1e11);
        SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 0.8,1.55);
        histo2DDummy3->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStatAUnshi, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
        graphCombEtaInvXSectionStatAUnshi->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionTotA, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombEtaInvXSectionTotA->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvXSectionSys, markerStyleDet[0] ,markerSizeDet[0]/2, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        graphPCMEtaInvXSectionSys->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaInvXSectionSys, markerStyleDet[2] ,markerSizeDet[2]/2, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        graphEMCALEtaInvXSectionSys->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaInvXSectionSys, markerStyleDet[4] ,markerSizeDet[4]/2, colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
        graphPCMEMCALEtaInvXSectionSys->Draw("pEsame");

        fitInvXSectionEta->SetLineColor(kBlue+2);
        fitInvXSectionEta->Draw("same");

        fitTsallisEtaPtMultFromShiftScaled->SetLineColor(kRed+2);
        fitTsallisEtaPtMultFromShiftScaled->Draw("same");

        canvasDummy2->Update();
        canvasDummy2->Print(Form("%s/ComparisonShiftedEta_8TeV.%s",outputDir.Data(),suffix.Data()));
        delete canvasDummy2;
        delete histo2DDummy3;
    }

    TGraphAsymmErrors* graphCombEtaInvXSectionStatA_WOXErr = (TGraphAsymmErrors*) graphCombEtaInvXSectionStatA->Clone("graphCombEtaInvXSectionStatA_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphCombEtaInvXSectionStatA_WOXErr);
    TGraphAsymmErrors* graphPCMEtaInvXSectionStat_WOXErr = (TGraphAsymmErrors*) graphPCMEtaInvXSectionStat->Clone("graphPCMEtaInvXSectionStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphPCMEtaInvXSectionStat_WOXErr);
    TGraphAsymmErrors* graphEMCALEtaInvXSectionStat_WOXErr = (TGraphAsymmErrors*) graphEMCALEtaInvXSectionStat->Clone("graphEMCALEtaInvXSectionStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphEMCALEtaInvXSectionStat_WOXErr);
    TGraphAsymmErrors* graphPCMEMCALEtaInvXSectionStat_WOXErr = (TGraphAsymmErrors*) graphPCMEMCALEtaInvXSectionStat->Clone("graphPCMEMCALEtaInvXSectionStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphPCMEMCALEtaInvXSectionStat_WOXErr);

    // *************************************************************************************************************
    // redo fitting after binshifts
    // *************************************************************************************************************
    // Tsallis function
    graphCombEtaInvXSectionTotA->Fit(fitInvXSectionEta,"QNRMEX0+","",0.4,35.);
    fitInvXSectionEta        = FitObject("l","fitInvCrossSectionEta8TeV","Eta",graphCombEtaInvXSectionTotA,0.4,35.,paramGraphEta,"QNRMEX0+");
    //fitInvXSectionEta        = FitObject("l","fitInvCrossSectionEta8TeV","Eta",graphCombEtaInvXSectionTotA,0.4,35.,paramGraphEta,"QNRMEX0+");
    cout << WriteParameterToFile(fitInvXSectionEta)<< endl;

    Double_t paramTCMEta[5]  = {graphCombEtaInvXSectionTotA->GetY()[1],0.2,graphCombEtaInvXSectionTotA->GetY()[3],0.75,3.};
     //Double_t paramTCMEta[5]  = {5E7,0.2,4E9,0.5,3.03};
    // Two component model by Bylinkin
    TF1* fitTCMInvXSectionEta= FitObject("tcm","fitTCMInvCrossSectionEta8TeV","Eta",graphCombEtaInvXSectionTotA,0.4,35.,paramTCMEta,"QNRMEX0+","", kFALSE);
    fitTCMInvXSectionEta     = FitObject("tcm","fitTCMInvCrossSectionEta8TeV","Eta",graphCombEtaInvXSectionTotA,0.4,35.,paramTCMEta,"QNRMEX0+","", kFALSE);

    TF1* fitTCMDecomposedEtaL                 = new TF1("twoCompModel_DecLow",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mesonMassExpectEta,mesonMassExpectEta,mesonMassExpectEta));
    fitTCMDecomposedEtaL->SetParameters(fitTCMInvXSectionEta->GetParameter(0),fitTCMInvXSectionEta->GetParameter(1));
    fitTCMDecomposedEtaL->SetRange(0.3,40);
    TF1 *fitTCMDecomposedEtaH                 = new TF1("twoCompModel_DecH","[0]/(TMath::Power(1+x*x/([1]*[1]*[2]),[2]))");
   //      graphCombEtaInvXSectionTotA->Fit(fitTCMDecomposedH,"QNRMEX0+","",5,20);
    fitTCMDecomposedEtaH->SetParameters(fitTCMInvXSectionEta->GetParameter(2),fitTCMInvXSectionEta->GetParameter(3), fitTCMInvXSectionEta->GetParameter(4));
    fitTCMDecomposedEtaH->SetRange(0.3,40);
    cout << WriteParameterToFile(fitTCMInvXSectionEta)<< endl;

    Double_t paramEtaPower[3] = {1E11,0.5,6.5};
    TF1* fitPowInvXSectionEta   = FitObject("powPure","fitPowInvXSectionEta8TeV","Eta",graphCombEtaInvXSectionTotA,3.5,35. ,paramEtaPower,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvXSectionEta)<< endl;

    Double_t paramEtaHageDorn[5] = {1E11,0.3,-0.1,0.5,5.95};
    TF1* fitOHagInvYieldEtaTot   = FitObject("oHag","fitOHagInvYieldEta8TeV","Eta",graphCombEtaInvXSectionTotA,0.3,35. ,paramEtaHageDorn,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitOHagInvYieldEtaTot)<< endl;

    TF1* fitPowInvXSectionEtaStat   = FitObject("powPure","fitPowInvXSectionEta8TeVStat","Eta",graphCombEtaInvXSectionStatA,3.5,35. ,paramEtaPower,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvXSectionEtaStat)<< endl;

    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "Eta - Tsallis" << endl;
    fLog << WriteParameterToFile(fitInvXSectionEta)<< endl;
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "Eta - TCM" << endl;
    fLog << WriteParameterToFile(fitTCMInvXSectionEta) << endl;
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "Eta - Hagedorn" << endl;
    fLog << WriteParameterToFile(fitOHagInvYieldEtaTot) << endl;
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "Eta - PowerLaw" << endl;
    fLog << WriteParameterToFile(fitPowInvXSectionEta) << endl;
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "Eta - PowerLaw - Stat" << endl;
    fLog << WriteParameterToFile(fitPowInvXSectionEtaStat) << endl;

    // *************************************************************************************************************
    // Shift spectra in Y  direction as well if desired
    // *************************************************************************************************************

    TGraphAsymmErrors* graphCombEtaInvXSectionTotA_yShifted         = NULL;
    TGraphAsymmErrors* graphCombEtaInvXSectionStatA_yShifted        = NULL;
    TGraphAsymmErrors* graphCombEtaInvXSectionSysA_yShifted         = NULL;
    TGraphAsymmErrors* graphPCMEtaInvXSectionStat_yShifted          = NULL;
    TGraphAsymmErrors* graphPCMEtaInvXSectionSys_yShifted           = NULL;
    TGraphAsymmErrors* graphEMCALEtaInvXSectionStat_yShifted        = NULL;
    TGraphAsymmErrors* graphEMCALEtaInvXSectionSys_yShifted         = NULL;
    TGraphAsymmErrors* graphPCMEMCALEtaInvXSectionStat_yShifted     = NULL;
    TGraphAsymmErrors* graphPCMEMCALEtaInvXSectionSys_yShifted      = NULL;

    if(bWCorrection.Contains("Y") ){
        graphCombEtaInvXSectionTotA_yShifted        = (TGraphAsymmErrors*)graphCombEtaInvXSectionTotAUnshi->Clone("EtaYShiftedCombTot");
        graphCombEtaInvXSectionTotA_yShifted        =  ApplyYshiftIndividualSpectra( graphCombEtaInvXSectionTotA_yShifted, fitInvXSectionEta);
        graphCombEtaInvXSectionStatA_yShifted       = (TGraphAsymmErrors*)graphCombEtaInvXSectionStatAUnshi->Clone("EtaYShiftedCombStat");
        graphCombEtaInvXSectionStatA_yShifted       =  ApplyYshiftIndividualSpectra( graphCombEtaInvXSectionStatA_yShifted, fitInvXSectionEta);
        graphCombEtaInvXSectionSysA_yShifted        = (TGraphAsymmErrors*)graphCombEtaInvXSectionSysAUnshi->Clone("EtaYShiftedCombSys");
        graphCombEtaInvXSectionSysA_yShifted        =  ApplyYshiftIndividualSpectra( graphCombEtaInvXSectionSysA_yShifted, fitInvXSectionEta);

        graphPCMEtaInvXSectionStat_yShifted         = (TGraphAsymmErrors*)graphPCMEtaInvXSectionStatUnshi->Clone("EtaYShiftedPCMStat");
        graphPCMEtaInvXSectionStat_yShifted         =  ApplyYshiftIndividualSpectra( graphPCMEtaInvXSectionStat_yShifted, fitInvXSectionEta);
        graphPCMEtaInvXSectionSys_yShifted          = (TGraphAsymmErrors*)graphPCMEtaInvXSectionSysUnshi->Clone("EtaYShiftedPCMSys");
        graphPCMEtaInvXSectionSys_yShifted          =  ApplyYshiftIndividualSpectra( graphPCMEtaInvXSectionSys_yShifted, fitInvXSectionEta);

        graphEMCALEtaInvXSectionStat_yShifted       = (TGraphAsymmErrors*)graphEMCALEtaInvXSectionStatUnshi->Clone("EtaYShiftedEMCStat");
        graphEMCALEtaInvXSectionStat_yShifted       =  ApplyYshiftIndividualSpectra( graphEMCALEtaInvXSectionStat_yShifted, fitInvXSectionEta);
        graphEMCALEtaInvXSectionSys_yShifted        = (TGraphAsymmErrors*)graphEMCALEtaInvXSectionSysUnshi->Clone("EtaYShiftedEMCSys");
        graphEMCALEtaInvXSectionSys_yShifted        =  ApplyYshiftIndividualSpectra( graphEMCALEtaInvXSectionSys_yShifted, fitInvXSectionEta);

        graphPCMEMCALEtaInvXSectionStat_yShifted    = (TGraphAsymmErrors*)graphPCMEMCALEtaInvXSectionStatUnshi->Clone("EtaYShiftedPCMEMCStat");
        graphPCMEMCALEtaInvXSectionStat_yShifted    =  ApplyYshiftIndividualSpectra( graphPCMEMCALEtaInvXSectionStat_yShifted, fitInvXSectionEta);
        graphPCMEMCALEtaInvXSectionSys_yShifted     = (TGraphAsymmErrors*)graphPCMEMCALEtaInvXSectionSysUnshi->Clone("EtaYShiftedPCMEMCStat");
        graphPCMEMCALEtaInvXSectionSys_yShifted     =  ApplyYshiftIndividualSpectra( graphPCMEMCALEtaInvXSectionSys_yShifted, fitInvXSectionEta);
    }

    //********************************************************************************************************
    // Plotting simple comparison of data vs fit to eta meson spec
    //********************************************************************************************************
    TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.08);
    canvasDummy2->SetLogy();
    canvasDummy2->SetLogx();
    TH2F* histo2DDummy3;
    histo2DDummy3               = new TH2F("histo2DDummy3","histo2DDummy3",1000,0.33,50.,1000,1,9e10);
    SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 0.8,1.55);
    histo2DDummy3->DrawCopy();

    DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionTotA, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
    graphCombEtaInvXSectionTotA->Draw("pEsame");

//    fitInvXSectionEta->SetLineColor(kBlue+2);
//    fitInvXSectionEta->Draw("same");
    fitTCMInvXSectionEta->SetLineColor(kRed+2);
    fitTCMInvXSectionEta->Draw("same");

    fitTCMDecomposedEtaL->SetLineColor(kAzure);
    fitTCMDecomposedEtaL->SetLineStyle(2);
    fitTCMDecomposedEtaL->Draw("same");
    fitTCMDecomposedEtaH->SetLineColor(kGreen+2);
    fitTCMDecomposedEtaH->SetLineStyle(8);
    fitTCMDecomposedEtaH->Draw("same");

    TLatex *labelTCMEta1= new TLatex(0.48, 0.94, Form("TCM low:"));
    TLatex *labelTCMEta2= new TLatex(0.48, 0.90, Form("A_{1}: (%.1e #pm %.1e) - T_{e}: (%.3f #pm %.3f)",fitTCMInvXSectionEta->GetParameter(0),fitTCMInvXSectionEta->GetParError(0),fitTCMInvXSectionEta->GetParameter(1),fitTCMInvXSectionEta->GetParError(1)));
    TLatex *labelTCMEta3= new TLatex(0.48, 0.86, Form("TCM high:"));
    TLatex *labelTCMEta4= new TLatex(0.48, 0.82, Form("A_{2}: (%.1e #pm %.1e) - T: (%.3f #pm %.3f) - n: (%.3f #pm %.3f)",fitTCMInvXSectionEta->GetParameter(2),fitTCMInvXSectionEta->GetParError(2),fitTCMInvXSectionEta->GetParameter(3),fitTCMInvXSectionEta->GetParError(3),fitTCMInvXSectionEta->GetParameter(4),fitTCMInvXSectionEta->GetParError(4)));

    TLatex *labelTCMEta5= new TLatex(0.55, 0.75, Form("Bylinkin-Rostovtsev:"));
    TLatex *labelTCMEta6= new TLatex(0.55, 0.71, Form("#it{A}_{1} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}_{2}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{n}"));

    SetStyleTLatex( labelTCMEta1, 0.03,4);
    labelTCMEta1->Draw();
    SetStyleTLatex( labelTCMEta2, 0.02,4);
    labelTCMEta2->Draw();
    SetStyleTLatex( labelTCMEta3, 0.03,4);
    labelTCMEta3->Draw();
    SetStyleTLatex( labelTCMEta4, 0.02,4);
    labelTCMEta4->Draw();
    SetStyleTLatex( labelTCMEta5, 0.03,4);
    labelTCMEta5->Draw();
    SetStyleTLatex( labelTCMEta6, 0.03,4);
    labelTCMEta6->Draw();

    TLatex *labelRelSysErrEnergyC    = new TLatex(0.18,0.94,collisionSystem8TeV.Data());
    SetStyleTLatex( labelRelSysErrEnergyC, 0.85*textSizeLabelsPixel,4);
    labelRelSysErrEnergyC->SetTextFont(43);
    labelRelSysErrEnergyC->Draw();
    TLatex *labelRelSysErrEtaC       = new TLatex(0.18,0.9,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrEtaC, 0.85*textSizeLabelsPixel,4);
    labelRelSysErrEtaC->SetTextFont(43);
    labelRelSysErrEtaC->Draw();

    TLegend* legendWithFitEta   = GetAndSetLegend2(0.17, 0.14, 0.5, 0.14+(0.035*3), 32);
    legendWithFitEta->AddEntry(fitTCMDecomposedEtaL,"TCM low","l");
    legendWithFitEta->AddEntry(fitTCMDecomposedEtaH,"TCM high","l");
    legendWithFitEta->AddEntry(fitTCMInvXSectionEta,"Bylinkin-Rostovtsev (TCM)","l");
    legendWithFitEta->Draw();

    canvasDummy2->Update();
    canvasDummy2->Print(Form("%s/ComparisonWithFitEta_8TeV.%s",outputDir.Data(),suffix.Data()));
//********************************************************************************************************
    canvasDummy2->Clear();
    histo2DDummy3->DrawCopy();

    graphCombEtaInvXSectionTotA->Draw("pEsame");

    fitInvXSectionEta->SetLineColor(kRed+2);
    fitInvXSectionEta->Draw("same");

    TLatex *labelTCMEta10 = new TLatex(0.35, 0.90, Form("dN/dy: (%.1e #pm %.1e) - n: (%.3f #pm %.3f) - T_{Levy} (GeV/c): (%.3f #pm %.3f)",fitInvXSectionEta->GetParameter(0),fitInvXSectionEta->GetParError(0),fitInvXSectionEta->GetParameter(1),fitInvXSectionEta->GetParError(1),fitInvXSectionEta->GetParameter(2),fitInvXSectionEta->GetParError(2)));
    SetStyleTLatex( labelTCMEta10, 0.02,4);
    labelTCMEta10->Draw();

    labelRelSysErrEnergyC->Draw();
    labelRelSysErrEtaC->Draw();

    TLegend* legendWithFitEta2   = GetAndSetLegend2(0.17, 0.14, 0.5, 0.14+(0.035*3), 32);
    legendWithFitEta2->AddEntry(fitInvXSectionEta,"Tsallis","l");
    legendWithFitEta2->Draw();

    canvasDummy2->Update();
    canvasDummy2->Print(Form("%s/ComparisonWithFit_Tsallis_Eta_8TeV.%s",outputDir.Data(),suffix.Data()));

    delete histo2DDummy3;
    canvasDummy2->Clear();

    //********************************************************************************************************
    // Plotting simple comparison of data vs fit to pi0 meson spec
    //********************************************************************************************************

    TF1* fitTCMDecomposedPi0L                 = new TF1("twoCompModel_DecLow",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0));
    fitTCMDecomposedPi0L->SetParameters(fitTCMInvXSectionPi0Plot->GetParameter(0),fitTCMInvXSectionPi0Plot->GetParameter(1));
    fitTCMDecomposedPi0L->SetRange(0.3,40);
    TF1 *fitTCMDecomposedPi0H                 = new TF1("twoCompModel_DecH","[0]/(TMath::Power(1+x*x/([1]*[1]*[2]),[2]))");
   //      graphCombEtaInvXSectionTotA->Fit(fitTCMDecomposedH,"QNRMEX0+","",5,20);
    fitTCMDecomposedPi0H->SetParameters(fitTCMInvXSectionPi0Plot->GetParameter(2),fitTCMInvXSectionPi0Plot->GetParameter(3), fitTCMInvXSectionPi0Plot->GetParameter(4));
    fitTCMDecomposedPi0H->SetRange(0.3,40);

    histo2DDummy3               = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.23,50.,1000,1,9e12);
    SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 0.8,1.55);
    histo2DDummy3->DrawCopy();

    DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatAUnshi, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
    graphCombPi0InvXSectionStatAUnshi->Draw("pEsame");
    DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatA, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
    graphCombPi0InvXSectionStatA->Draw("pEsame");

//    fitInvXSectionPi0->SetLineColor(kBlue+2);
//    fitInvXSectionPi0->Draw("same");
    fitTCMInvXSectionPi0Plot->SetLineColor(kRed+2);
    fitTCMInvXSectionPi0Plot->SetRange(0.3,40.);
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
    legendWithFitPi02->AddEntry(fitInvXSectionPi0,"Tsallis","l");
    legendWithFitPi02->Draw();

    canvasDummy2->Update();
    canvasDummy2->Print(Form("%s/ComparisonWithFit_Tsallis_Pi0_8TeV.%s",outputDir.Data(),suffix.Data()));

    delete canvasDummy2;
    delete histo2DDummy3;


    // *************************************************************************************************************
    // Calculate ratios to combined fit
    // *************************************************************************************************************
    TGraph* graphRatioEtaCombNLOMuHalf                  = (TGraph*)graphNLOCalcEtaMuHalf->Clone();
    TGraph* graphRatioEtaCombNLOMuOne                   = (TGraph*)graphNLOCalcEtaMuOne->Clone();
    TGraph* graphRatioEtaCombNLOMuTwo                   = (TGraph*)graphNLOCalcEtaMuTwo->Clone();
    graphRatioEtaCombNLOMuHalf                          = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuHalf, fitTCMInvXSectionEta);
    graphRatioEtaCombNLOMuOne                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuOne, fitTCMInvXSectionEta);
    graphRatioEtaCombNLOMuTwo                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuTwo, fitTCMInvXSectionEta);

    TGraphAsymmErrors* graphRatioEtaAESSS               = (TGraphAsymmErrors*) graphEtaAESSS->Clone();
    graphRatioEtaAESSS                                  = CalculateGraphErrRatioToFit (graphRatioEtaAESSS, fitTCMInvXSectionEta);

    TGraphAsymmErrors* graphRatioEtaCombCombFitTotA     = (TGraphAsymmErrors*)graphCombEtaInvXSectionTotA->Clone();
    graphRatioEtaCombCombFitTotA                        = CalculateGraphErrRatioToFit(graphRatioEtaCombCombFitTotA, fitTCMInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaCombCombFitStatA    = (TGraphAsymmErrors*)graphCombEtaInvXSectionStatA->Clone();
    graphRatioEtaCombCombFitStatA                       = CalculateGraphErrRatioToFit(graphRatioEtaCombCombFitStatA, fitTCMInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaCombCombFitSysA     = (TGraphAsymmErrors*)graphCombEtaInvXSectionSysA->Clone();
    graphRatioEtaCombCombFitSysA                        = CalculateGraphErrRatioToFit(graphRatioEtaCombCombFitSysA, fitTCMInvXSectionEta);

    TGraphAsymmErrors* graphRatioEtaCombTsallisFitStatA  = (TGraphAsymmErrors*)graphCombEtaInvXSectionStatA->Clone();
    graphRatioEtaCombTsallisFitStatA                     = CalculateGraphErrRatioToFit(graphRatioEtaCombTsallisFitStatA, fitInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaCombHagedornFitStatA = (TGraphAsymmErrors*)graphCombEtaInvXSectionStatA->Clone();
    graphRatioEtaCombHagedornFitStatA                    = CalculateGraphErrRatioToFit(graphRatioEtaCombHagedornFitStatA, fitOHagInvYieldEtaTot);
    TGraphAsymmErrors* graphRatioEtaCombPowerFitStatA    = (TGraphAsymmErrors*)graphCombEtaInvXSectionStatA->Clone();
    graphRatioEtaCombPowerFitStatA                       = CalculateGraphErrRatioToFit(graphRatioEtaCombPowerFitStatA, fitPowInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaCombTsallisFitSysA   = (TGraphAsymmErrors*)graphCombEtaInvXSectionSysA->Clone();
    graphRatioEtaCombTsallisFitSysA                      = CalculateGraphErrRatioToFit(graphRatioEtaCombTsallisFitSysA, fitInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaCombHagedornFitSysA  = (TGraphAsymmErrors*)graphCombEtaInvXSectionSysA->Clone();
    graphRatioEtaCombHagedornFitSysA                     = CalculateGraphErrRatioToFit(graphRatioEtaCombHagedornFitSysA, fitOHagInvYieldEtaTot);
    TGraphAsymmErrors* graphRatioEtaCombPowerFitSysA     = (TGraphAsymmErrors*)graphCombEtaInvXSectionSysA->Clone();
    graphRatioEtaCombPowerFitSysA                        = CalculateGraphErrRatioToFit(graphRatioEtaCombPowerFitSysA, fitPowInvXSectionEta);

    TGraphAsymmErrors* graphRatioEtaPCMCombFitStat      = (TGraphAsymmErrors*)graphPCMEtaInvXSectionStat->Clone();
    graphRatioEtaPCMCombFitStat                         = CalculateGraphErrRatioToFit(graphRatioEtaPCMCombFitStat, fitTCMInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaPCMCombFitSys       = (TGraphAsymmErrors*)graphPCMEtaInvXSectionSys->Clone();
    graphRatioEtaPCMCombFitSys                          = CalculateGraphErrRatioToFit(graphRatioEtaPCMCombFitSys, fitTCMInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaEMCALCombFitStat    = (TGraphAsymmErrors*)graphEMCALEtaInvXSectionStat->Clone();
    graphRatioEtaEMCALCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioEtaEMCALCombFitStat, fitTCMInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaEMCALCombFitSys     = (TGraphAsymmErrors*)graphEMCALEtaInvXSectionSys->Clone();
    graphRatioEtaEMCALCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioEtaEMCALCombFitSys, fitTCMInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaPCMEMCALCombFitStat = (TGraphAsymmErrors*)graphPCMEMCALEtaInvXSectionStat->Clone();
    graphRatioEtaPCMEMCALCombFitStat                    = CalculateGraphErrRatioToFit(graphRatioEtaPCMEMCALCombFitStat, fitTCMInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaPCMEMCALCombFitSys  = (TGraphAsymmErrors*)graphPCMEMCALEtaInvXSectionSys->Clone();
    graphRatioEtaPCMEMCALCombFitSys                     = CalculateGraphErrRatioToFit(graphRatioEtaPCMEMCALCombFitSys, fitTCMInvXSectionEta);

    TH1D* histoRatioPythia8ToFitEta                  = (TH1D*) histoPythia8InvXSectionEta->Clone();
    histoRatioPythia8ToFitEta                        = CalculateHistoRatioToFit (histoRatioPythia8ToFitEta, fitTCMInvXSectionEta);
    histoRatioPythia8ToFitEta->GetXaxis()->SetRangeUser(0.8,35);

    TH1D* histoRatioPythia8_4CToFitEta               = (TH1D*) histoPythia8_4CInvXSectionEta->Clone();
    histoRatioPythia8_4CToFitEta                     = CalculateHistoRatioToFit (histoRatioPythia8_4CToFitEta, fitTCMInvXSectionEta);
    histoRatioPythia8_4CToFitEta->GetXaxis()->SetRangeUser(0.8,35);

    TGraphErrors* graphRatioPythia8ToFitEta             = (TGraphErrors*) graphPythia8InvXSectionEta->Clone();
    graphRatioPythia8ToFitEta                           = CalculateGraphErrRatioToFit (graphRatioPythia8ToFitEta, fitTCMInvXSectionEta);
    while(graphRatioPythia8ToFitEta->GetX()[0] < 0.8) graphRatioPythia8ToFitEta->RemovePoint(0);
    TGraphErrors* graphRatioPythia8_4CToFitEta          = (TGraphErrors*) graphPythia8_4CInvXSectionEta->Clone();
    graphRatioPythia8_4CToFitEta                        = CalculateGraphErrRatioToFit (graphRatioPythia8_4CToFitEta, fitTCMInvXSectionEta);
    while(graphRatioPythia8_4CToFitEta->GetX()[0] < 0.8) graphRatioPythia8_4CToFitEta->RemovePoint(0);

    TGraphAsymmErrors* graphRatioEtaCombCombFitStatA_WOXErr = (TGraphAsymmErrors*) graphRatioEtaCombCombFitStatA->Clone("graphRatioEtaCombCombFitStatA_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombCombFitStatA_WOXErr);

    TGraphAsymmErrors* graphRatioEtaCombTsallisFitStatA_WOXErr = (TGraphAsymmErrors*) graphRatioEtaCombTsallisFitStatA->Clone("graphRatioEtaCombTsallisFitStatA_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombTsallisFitStatA_WOXErr);
    TGraphAsymmErrors* graphRatioEtaCombHagedornFitStatA_WOXErr = (TGraphAsymmErrors*) graphRatioEtaCombHagedornFitStatA->Clone("graphRatioEtaCombHagedornFitStatA_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombHagedornFitStatA_WOXErr);
    TGraphAsymmErrors* graphRatioEtaCombPowerFitStatA_WOXErr = (TGraphAsymmErrors*) graphRatioEtaCombPowerFitStatA->Clone("graphRatioEtaCombPowerFitStatA_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombPowerFitStatA_WOXErr);

    TGraphAsymmErrors* graphRatioEtaPCMCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioEtaPCMCombFitStat->Clone("graphRatioEtaPCMCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioEtaPCMCombFitStat_WOXErr);
    TGraphAsymmErrors* graphRatioEtaEMCALCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioEtaEMCALCombFitStat->Clone("graphRatioEtaEMCALCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioEtaEMCALCombFitStat_WOXErr);
    TGraphAsymmErrors* graphRatioEtaPCMEMCALCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioEtaPCMEMCALCombFitStat->Clone("graphRatioEtaPCMEMCALCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioEtaPCMEMCALCombFitStat_WOXErr);

    // **********************************************************************************************************************
    // ******************************************* Ratio of Comb to Fit ****************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 48;
    canvasRatioToCombFit->cd();
    TH2F * histo2DEtaRatioToCombFit;
    histo2DEtaRatioToCombFit               = new TH2F("histo2DEtaRatioToCombFit","histo2DEtaRatioToCombFit",1000,0.33,50.,1000,0.2,7.    );
    SetStyleHistoTH2ForGraphs(histo2DEtaRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/TCM fit", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
    histo2DEtaRatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DEtaRatioToCombFit->GetXaxis()->SetNoExponent(kTRUE);
//    histo2DEtaRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
    histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.3,1.8);
    histo2DEtaRatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSysA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioEtaCombCombFitSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatA_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioEtaCombCombFitStatA_WOXErr->Draw("p,same,z");

        DrawGammaLines(0.33, 50. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.33, 50. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.33, 50. , 0.9, 0.9,0.5, kGray, 7);

        TLatex *labelRatioToFitEnergy2   = new TLatex(0.73, 0.91, collisionSystem8TeV.Data());
        SetStyleTLatex( labelRatioToFitEnergy2, textSizeLabelsPixel,4);
        labelRatioToFitEnergy2->SetTextFont(43);
        labelRatioToFitEnergy2->Draw();
        TLatex *labelRatioToFitALICE2    = new TLatex(0.85, 0.86, "ALICE");
        SetStyleTLatex( labelRatioToFitALICE2, textSizeLabelsPixel,4);
        labelRatioToFitALICE2->SetTextFont(43);
        labelRatioToFitALICE2->Draw();

        TLatex *labelRatioToFitEta      = new TLatex(0.84,0.82,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitEta, textSizeLabelsPixel,4);
        labelRatioToFitEta->SetTextFont(43);
        labelRatioToFitEta->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfCombToCombFit_PP8TeV.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // *******************************************Plot different ratios to fits *********************************************
    // **********************************************************************************************************************

    histo2DEtaRatioToCombFit->SetYTitle("Data/fit");
    histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.1,2.45);
    histo2DEtaRatioToCombFit->Draw("copy");

    while (graphRatioEtaCombPowerFitSysA->GetX()[0] < 1.6) graphRatioEtaCombPowerFitSysA->RemovePoint(0);
    while (graphRatioEtaCombPowerFitStatA_WOXErr->GetX()[0] < 1.6) graphRatioEtaCombPowerFitStatA_WOXErr->RemovePoint(0);

        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombPowerFitSysA, 33, markerSizeComb*2, kAzure+2 , kAzure+2, widthLinesBoxes, kTRUE);
        graphRatioEtaCombPowerFitSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombPowerFitStatA_WOXErr, 33, markerSizeComb*2, kAzure+2 , kAzure+2);
        graphRatioEtaCombPowerFitStatA_WOXErr->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombTsallisFitSysA, 25, markerSizeComb, colorTrigg[1], colorTrigg[1], widthLinesBoxes, kTRUE);
        graphRatioEtaCombTsallisFitSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombTsallisFitStatA_WOXErr, 25, markerSizeComb, colorTrigg[1] , colorTrigg[1]);
        graphRatioEtaCombTsallisFitStatA_WOXErr->Draw("p,same,z");

        graphRatioEtaCombCombFitSysA->SetMarkerColor(kGray+2);
        graphRatioEtaCombCombFitSysA->Draw("E2same");
        graphRatioEtaCombCombFitStatA_WOXErr->SetMarkerColor(kGray+2);
        graphRatioEtaCombCombFitStatA_WOXErr->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombHagedornFitSysA, 24, markerSizeComb+0.2, colorTrigg[2], colorTrigg[2], widthLinesBoxes, kTRUE);
        graphRatioEtaCombHagedornFitSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombHagedornFitStatA_WOXErr, 24, markerSizeComb+0.2, colorTrigg[2] , colorTrigg[2]);
        graphRatioEtaCombHagedornFitStatA_WOXErr->Draw("p,same,z");

        DrawGammaLines(0.33, 50. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.33, 50. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.33, 50. , 0.9, 0.9,0.1, kGray, 7);

        labelRatioToFitEnergy2->Draw();
        labelRatioToFitALICE2->Draw();
        labelRatioToFitEta->Draw();

        TLegend* legendRatioEtaFits= GetAndSetLegend2(0.12,0.95-4*1.05*textsizeLabelsPP,0.37,0.95, textSizeLabelsPixel, 1, "", 43, 0.2);
        legendRatioEtaFits->AddEntry(graphRatioEtaCombCombFitSysA,"TCM","p");
        legendRatioEtaFits->AddEntry(graphRatioEtaCombTsallisFitSysA,"Tsallis","p");
        legendRatioEtaFits->AddEntry(graphRatioEtaCombHagedornFitSysA,"mod. Hagedorn","p");
        legendRatioEtaFits->AddEntry(graphRatioEtaCombPowerFitSysA,"pure powerlaw, 3.5-35 GeV/#it{c}","p");
        legendRatioEtaFits->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfCombToDifferentFits_PP8TeV.%s",outputDir.Data(),suffix.Data()));
    histo2DEtaRatioToCombFit->SetYTitle("Data/TCM fit");

    // **********************************************************************************************************************
    // ******************************************* Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.55);
    histo2DEtaRatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioEtaPCMCombFitSys, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaPCMCombFitStat_WOXErr, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaEMCALCombFitSys, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaEMCALCombFitStat_WOXErr, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaPCMEMCALCombFitSys, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaPCMEMCALCombFitStat_WOXErr, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4]);

        graphRatioEtaPCMCombFitSys->Draw("E2same");
        graphRatioEtaEMCALCombFitSys->Draw("E2same");
        graphRatioEtaPCMEMCALCombFitSys->Draw("E2same");

        graphRatioEtaPCMCombFitStat_WOXErr->Draw("p,same,z");
        graphRatioEtaEMCALCombFitStat_WOXErr->Draw("p,same,z");
        graphRatioEtaPCMEMCALCombFitStat_WOXErr->Draw("p,same,z");

        DrawGammaLines(0.33, 50. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.33, 50. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.33, 50. , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy2->Draw();
        labelRatioToFitALICE2->Draw();
        labelRatioToFitEta->Draw();

        //****************************** Definition of the Legend ******************************************
        Double_t rowsLegendOnlyEtaRatio[4]          = {0.91, 0.86, 0.81, 0.76};
        Double_t rowsLegendOnlyEtaRatioAbs[4]       = {0.91, 2.245, 2.11, 1.975};
        Double_t columnsLegendOnlyEtaRatio[3]       = {0.115, 0.355, 0.43};
        Double_t columnsLegendOnlyEtaRatioAbs[3]    = {0.115, 1.38, 1.92};

        //****************** first Column **************************************************
        TLatex *textPCMOnlyRatioEta                 = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[1],"PCM");
        SetStyleTLatex( textPCMOnlyRatioEta, textSizeLabelsPixel,4);
        textPCMOnlyRatioEta->SetTextFont(43);
        textPCMOnlyRatioEta->Draw();
        TLatex *textEMCALOnlyRatioEta               = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[2],"EMC");
        SetStyleTLatex( textEMCALOnlyRatioEta, textSizeLabelsPixel,4);
        textEMCALOnlyRatioEta->SetTextFont(43);
        textEMCALOnlyRatioEta->Draw();
        TLatex *textPCMEMCALOnlyRatioEta            = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[3],"PCM-EMC");
        SetStyleTLatex( textPCMEMCALOnlyRatioEta, textSizeLabelsPixel,4);
        textPCMEMCALOnlyRatioEta->SetTextFont(43);
        textPCMEMCALOnlyRatioEta->Draw();

        //****************** second Column *************************************************
        TLatex *textStatOnlyRatioEta                = new TLatex(columnsLegendOnlyEtaRatio[1]-0.05,rowsLegendOnlyEtaRatio[0] ,"stat");
        SetStyleTLatex( textStatOnlyRatioEta, textSizeLabelsPixel,4);
        textStatOnlyRatioEta->SetTextFont(43);
        textStatOnlyRatioEta->Draw();
        TLatex *textSysOnlyRatioEta                 = new TLatex(columnsLegendOnlyEtaRatio[2]-0.05 ,rowsLegendOnlyEtaRatio[0],"syst");
        SetStyleTLatex( textSysOnlyRatioEta, textSizeLabelsPixel,4);
        textSysOnlyRatioEta->SetTextFont(43);
        textSysOnlyRatioEta->Draw();
        TMarker* markerPCMEtaOnlyRatioEta           = CreateMarkerFromGraph(graphRatioEtaPCMCombFitSys,columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[1],1);
        markerPCMEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[1]);
        TMarker* markerEMCALEtaOnlyRatioEta         = CreateMarkerFromGraph(graphRatioEtaEMCALCombFitSys, columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[2],1);
        markerEMCALEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[2]);
        TMarker* markerPCMEMCALEtaOnlyRatioEta      = CreateMarkerFromGraph(graphRatioEtaPCMEMCALCombFitSys, columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[3],1);
        markerPCMEMCALEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[3]);

        TBox* boxPCMEtaOnlyRatioEta                 = CreateBoxFromGraph(graphRatioEtaPCMCombFitSys, columnsLegendOnlyEtaRatioAbs[2]-0.5*lengthBox , rowsLegendOnlyEtaRatioAbs[1]- heightBox,
                                                                         columnsLegendOnlyEtaRatioAbs[2]+ 2*lengthBox, rowsLegendOnlyEtaRatioAbs[1]+ heightBox);
        boxPCMEtaOnlyRatioEta->Draw("l");
        TBox* boxEMCALEtaOnlyRatioEta               = CreateBoxFromGraph(graphRatioEtaEMCALCombFitSys, columnsLegendOnlyEtaRatioAbs[2]-0.5*lengthBox , rowsLegendOnlyEtaRatioAbs[2]- heightBox,
                                                                         columnsLegendOnlyEtaRatioAbs[2]+ 2*lengthBox, rowsLegendOnlyEtaRatioAbs[2]+ heightBox);
        boxEMCALEtaOnlyRatioEta->Draw("l");
        TBox* boxPCMEMCALEtaOnlyRatioEta            = CreateBoxFromGraph(graphRatioEtaPCMEMCALCombFitSys, columnsLegendOnlyEtaRatioAbs[2]-0.5*lengthBox , rowsLegendOnlyEtaRatioAbs[3]- heightBox,
                                                                         columnsLegendOnlyEtaRatioAbs[2]+ 2*lengthBox, rowsLegendOnlyEtaRatioAbs[3]+ heightBox);
        boxPCMEMCALEtaOnlyRatioEta->Draw("l");

    canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfIndividualMeasToCombFit_PP.%s",outputDir.Data(),suffix.Data()));

    // *************************************************************************************************************
    // Calculate ratios to combined TSALLIS fit
    // *************************************************************************************************************

    TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitTotA     = (TGraphAsymmErrors*)graphCombEtaInvXSectionTotA->Clone();
    graphRatioEtaTsallisCombCombFitTotA                        = CalculateGraphErrRatioToFit(graphRatioEtaTsallisCombCombFitTotA, fitInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitStatA    = (TGraphAsymmErrors*)graphCombEtaInvXSectionStatA->Clone();
    graphRatioEtaTsallisCombCombFitStatA                       = CalculateGraphErrRatioToFit(graphRatioEtaTsallisCombCombFitStatA, fitInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitSysA     = (TGraphAsymmErrors*)graphCombEtaInvXSectionSysA->Clone();

    graphRatioEtaTsallisCombCombFitSysA                        = CalculateGraphErrRatioToFit(graphRatioEtaTsallisCombCombFitSysA, fitInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaTsallisPCMCombFitStat      = (TGraphAsymmErrors*)graphPCMEtaInvXSectionStat->Clone();
    graphRatioEtaTsallisPCMCombFitStat                         = CalculateGraphErrRatioToFit(graphRatioEtaTsallisPCMCombFitStat, fitInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaTsallisPCMCombFitSys       = (TGraphAsymmErrors*)graphPCMEtaInvXSectionSys->Clone();
    graphRatioEtaTsallisPCMCombFitSys                          = CalculateGraphErrRatioToFit(graphRatioEtaTsallisPCMCombFitSys, fitInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaTsallisEMCALCombFitStat    = (TGraphAsymmErrors*)graphEMCALEtaInvXSectionStat->Clone();
    graphRatioEtaTsallisEMCALCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioEtaTsallisEMCALCombFitStat, fitInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaTsallisEMCALCombFitSys     = (TGraphAsymmErrors*)graphEMCALEtaInvXSectionSys->Clone();
    graphRatioEtaTsallisEMCALCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioEtaTsallisEMCALCombFitSys, fitInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaTsallisPCMEMCALCombFitStat = (TGraphAsymmErrors*)graphPCMEMCALEtaInvXSectionStat->Clone();
    graphRatioEtaTsallisPCMEMCALCombFitStat                    = CalculateGraphErrRatioToFit(graphRatioEtaTsallisPCMEMCALCombFitStat, fitInvXSectionEta);
    TGraphAsymmErrors* graphRatioEtaTsallisPCMEMCALCombFitSys  = (TGraphAsymmErrors*)graphPCMEMCALEtaInvXSectionSys->Clone();
    graphRatioEtaTsallisPCMEMCALCombFitSys                     = CalculateGraphErrRatioToFit(graphRatioEtaTsallisPCMEMCALCombFitSys, fitInvXSectionEta);

    TGraphAsymmErrors* graphRatioEtaTsallisCombCombFitStatA_WOXErr = (TGraphAsymmErrors*) graphRatioEtaTsallisCombCombFitStatA->Clone("graphRatioEtaTsallisCombCombFitStatA_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioEtaTsallisCombCombFitStatA_WOXErr);
    TGraphAsymmErrors* graphRatioEtaTsallisPCMCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioEtaTsallisPCMCombFitStat->Clone("graphRatioEtaTsallisPCMCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioEtaTsallisPCMCombFitStat_WOXErr);
    TGraphAsymmErrors* graphRatioEtaTsallisEMCALCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioEtaTsallisEMCALCombFitStat->Clone("graphRatioEtaTsallisEMCALCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioEtaTsallisEMCALCombFitStat_WOXErr);
    TGraphAsymmErrors* graphRatioEtaTsallisPCMEMCALCombFitStat_WOXErr = (TGraphAsymmErrors*) graphRatioEtaTsallisPCMEMCALCombFitStat->Clone("graphRatioEtaTsallisPCMEMCALCombFitStat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioEtaTsallisPCMEMCALCombFitStat_WOXErr);

    // **********************************************************************************************************************
    // ******************************************* Ratio of Comb to Fit ****************************************
    // **********************************************************************************************************************
    canvasRatioToCombFit->cd();
    histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.3,1.8);
    histo2DEtaRatioToCombFit->SetYTitle("Data/Tsallis fit");
    histo2DEtaRatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioEtaTsallisCombCombFitSysA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioEtaTsallisCombCombFitSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaTsallisCombCombFitStatA_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioEtaTsallisCombCombFitStatA_WOXErr->Draw("p,same,z");

        DrawGammaLines(0.33, 50. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.33, 50. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.33, 50. , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy2->Draw();
        labelRatioToFitALICE2->Draw();
        labelRatioToFitEta->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfCombToCombTsallisFit_PP8TeV.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************************* Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.55);
    histo2DEtaRatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioEtaTsallisPCMCombFitSys, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaTsallisPCMCombFitStat_WOXErr, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaTsallisEMCALCombFitSys, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaTsallisEMCALCombFitStat_WOXErr, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaTsallisPCMEMCALCombFitSys, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaTsallisPCMEMCALCombFitStat_WOXErr, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4]);

        graphRatioEtaTsallisPCMCombFitSys->Draw("E2same");
        graphRatioEtaTsallisEMCALCombFitSys->Draw("E2same");
        graphRatioEtaTsallisPCMEMCALCombFitSys->Draw("E2same");

        graphRatioEtaTsallisPCMCombFitStat_WOXErr->Draw("p,same,z");
        graphRatioEtaTsallisEMCALCombFitStat_WOXErr->Draw("p,same,z");
        graphRatioEtaTsallisPCMEMCALCombFitStat_WOXErr->Draw("p,same,z");

        DrawGammaLines(0.33, 50. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.33, 50. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.33, 50. , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy2->Draw();
        labelRatioToFitALICE2->Draw();
        labelRatioToFitEta->Draw();

        //****************** first Column **************************************************
        textPCMOnlyRatioEta->Draw();
        textEMCALOnlyRatioEta->Draw();
        textPCMEMCALOnlyRatioEta->Draw();

        //****************** second Column *************************************************
        textStatOnlyRatioEta->Draw();
        textSysOnlyRatioEta->Draw();

        markerPCMEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[1]);
        markerEMCALEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[2]);
        markerPCMEMCALEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[3]);

        boxPCMEtaOnlyRatioEta->Draw("l");
        boxEMCALEtaOnlyRatioEta->Draw("l");
        boxPCMEMCALEtaOnlyRatioEta->Draw("l");

    canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfIndividualMeasToCombTsallisFit_PP.%s",outputDir.Data(),suffix.Data()));

    // *******************************************************************************************************
    // ************************** Combination of different eta/pi0 measurements **********************************
    // *******************************************************************************************************
    // REMARKS:
    //     - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - currently only PCM-EMCAL vs others fully implemeted energy independent
    //     - extendable to other energies
    //     - offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented

    // definition of array of histograms (NULL - means we have no measurement at this energy for this rec-method)
    // for statistical error and final value from respective method
    TH1D* statErrorCollectionEtaToPi0[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorCollectionEtaToPi0[i]   = NULL;
    }
    statErrorCollectionEtaToPi0[0]       = (TH1D*)histoPCMEtaToPi0Stat->Clone("statErrPCMEtaToPi0");
    statErrorCollectionEtaToPi0[2]       = (TH1D*)histoEMCALEtaToPi0Stat->Clone("statErrEMCALEtaToPi0");
    statErrorCollectionEtaToPi0[4]       = (TH1D*)histoPCMEMCALEtaToPi0Stat->Clone("statErrPCMEMCALEtaToPi0");

    // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
    // for systematic error from respective method
    TGraphAsymmErrors* sysErrorCollectionEtaToPi0[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorCollectionEtaToPi0[i]    = NULL;
    }
    sysErrorCollectionEtaToPi0[0]        = (TGraphAsymmErrors*)graphPCMEtaToPi0Sys->Clone("sysErrPCMEtaToPi0");
    sysErrorCollectionEtaToPi0[2]        = (TGraphAsymmErrors*)graphEMCALEtaToPi0Sys->Clone("sysErrEMCALEta");
    sysErrorCollectionEtaToPi0[4]        = (TGraphAsymmErrors*)graphPCMEMCALEtaToPi0Sys->Clone("sysErrPCMEMCALEtaToPi0");

    TH1D* statErrorRelCollectionEtaToPi0[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorRelCollectionEtaToPi0[i]        = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (statErrorCollectionEtaToPi0[i])
            statErrorRelCollectionEtaToPi0[i]    = CalculateRelErrUpTH1D( statErrorCollectionEtaToPi0[i], Form("relativeStatErrorEtaToPi0_%s", nameMeasGlobal[i].Data()));
    }

    TGraphAsymmErrors* sysErrorRelCollectionEtaToPi0[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorRelCollectionEtaToPi0[i]         = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (sysErrorCollectionEtaToPi0[i])
            sysErrorRelCollectionEtaToPi0[i]     = CalculateRelErrUpAsymmGraph( sysErrorCollectionEtaToPi0[i], Form("relativeSysErrorEtaToPi0_%s", nameMeasGlobal[i].Data()));
    }

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for EtaToPi0
    Int_t offSetsEtaToPi0[11]           = { -1,  0,  0,  0,  0,
                                            0,  0,  0,  0,  0,
                                            0};
    Int_t offSetsEtaToPi0Sys[11]        = { 0,  0,  5,  0,  2,
                                            0,  0,  0,  0,  0,
                                            0};
    Int_t offSetEtaToPi0Shifting[11]    = { 0,  0,  0,  0,  0,
                                            0,  0,  0,  0,  0,
                                            0};
    Int_t nComBinsEtaToPi0Shifting[11]  = { 0,  0,  0,  0,  0,
                                            0,  0,  0,  0,  0,
                                            0};

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************

    TGraph* graphWeightsEtaToPi0A[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsEtaToPi0A[i]                    = NULL;
    }

    maxNBinsEtaW0Merged-=2;

    // Declaration & calculation of combined spectrum
    TString fileNameEtaToPi0OutputWeightingA        = Form("%s/EtaToPi0_WeightingMethodA.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombEtaToPi0StatA       = NULL;
    TGraphAsymmErrors* graphCombEtaToPi0SysA        = NULL;
    TGraphAsymmErrors* graphCombEtaToPi0TotA        = CombinePtPointsSpectraFullCorrMat(   statErrorCollectionEtaToPi0,    sysErrorCollectionEtaToPi0,
                                                                                           xPtLimitsEtaWOMerged, maxNBinsEtaW0Merged,
                                                                                           offSetsEtaToPi0, offSetsEtaToPi0Sys,
                                                                                           graphCombEtaToPi0StatA, graphCombEtaToPi0SysA,
                                                                                           fileNameEtaToPi0OutputWeightingA,"8TeV", "EtaToPi0", kFALSE,
                                                                                           0x0, fileInputCorrFactors
                                                                                        );
    graphCombEtaToPi0StatA->Print();
    if (graphCombEtaToPi0StatA == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }

//    graphCombEtaToPi0StatA->RemovePoint(0);
//    graphCombEtaToPi0SysA->RemovePoint(0);
//    graphCombEtaToPi0TotA->RemovePoint(0);

//     return;

    // Reading weights from output file for plotting
    ifstream fileWeightsEtaToPi0ReadA;
    fileWeightsEtaToPi0ReadA.open(fileNameEtaToPi0OutputWeightingA,ios_base::in);
    cout << "reading" << fileNameEtaToPi0OutputWeightingA << endl;
    Double_t xValuesEtaToPi0ReadA[70];
    Double_t weightsEtaToPi0ReadA[11][70];
    Int_t availableEtaToPi0MeasA[11]        = { -1, -1, -1, -1, -1,
                                        -1, -1, -1, -1, -1,
                                        -1};
    Int_t nMeasSetEtaToPi0A                 = 3;
    Int_t nPtBinsEtaToPi0ReadA              = 0;
    while(!fileWeightsEtaToPi0ReadA.eof() && nPtBinsEtaToPi0ReadA < 70){
        TString garbage             = "";
        if (nPtBinsEtaToPi0ReadA == 0){
            fileWeightsEtaToPi0ReadA >> garbage ;//>> availableEtaToPi0Meas[0] >> availableEtaToPi0Meas[1] >> availableEtaToPi0Meas[2] >> availableEtaToPi0Meas[3];
            for (Int_t i = 0; i < nMeasSetEtaToPi0A; i++){
                fileWeightsEtaToPi0ReadA >> availableEtaToPi0MeasA[i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availableEtaToPi0MeasA[i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsEtaToPi0ReadA >> xValuesEtaToPi0ReadA[nPtBinsEtaToPi0ReadA-1];
            for (Int_t i = 0; i < nMeasSetEtaToPi0A; i++){
                fileWeightsEtaToPi0ReadA >> weightsEtaToPi0ReadA[availableEtaToPi0MeasA[i]][nPtBinsEtaToPi0ReadA-1] ;
            }
            cout << "read: "<<  nPtBinsEtaToPi0ReadA << "\t"<< xValuesEtaToPi0ReadA[nPtBinsEtaToPi0ReadA-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetEtaToPi0A; i++){
                cout << weightsEtaToPi0ReadA[availableEtaToPi0MeasA[i]][nPtBinsEtaToPi0ReadA-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsEtaToPi0ReadA++;
    }
    nPtBinsEtaToPi0ReadA                    = nPtBinsEtaToPi0ReadA-2 ;
    fileWeightsEtaToPi0ReadA.close();

    for (Int_t i = 0; i < nMeasSetEtaToPi0A; i++){
        graphWeightsEtaToPi0A[availableEtaToPi0MeasA[i]]                        = new TGraph(nPtBinsEtaToPi0ReadA,xValuesEtaToPi0ReadA,weightsEtaToPi0ReadA[availableEtaToPi0MeasA[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsEtaToPi0ReadA; n++){
            if (graphWeightsEtaToPi0A[availableEtaToPi0MeasA[i]]->GetY()[bin] == 0) graphWeightsEtaToPi0A[availableEtaToPi0MeasA[i]]->RemovePoint(bin);
            else bin++;
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Plotting weights Method A ************************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel           = 900*0.04;

    canvasWeights->cd();

    TH2F * histo2DEtaToPi0Weights;
    histo2DEtaToPi0Weights = new TH2F("histo2DEtaToPi0Weights","histo2DEtaToPi0Weights",11000,0.33,30.,1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DEtaToPi0Weights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DEtaToPi0Weights->GetXaxis()->SetMoreLogLabels();
    histo2DEtaToPi0Weights->GetXaxis()->SetNoExponent(kTRUE);
//    histo2DEtaToPi0Weights->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DEtaToPi0Weights->GetYaxis()->SetRangeUser(-10,10);
    histo2DEtaToPi0Weights->Draw("copy");

        TLegend* legendWeightsEtaToPi0   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEtaToPi0A), 32);
        for (Int_t i = 0; i < nMeasSetEtaToPi0A; i++){
            DrawGammaSetMarkerTGraph(graphWeightsEtaToPi0A[availableEtaToPi0MeasA[i]], markerStyleDet[availableEtaToPi0MeasA[i]], markerSizeDet[availableEtaToPi0MeasA[i]]*0.5, colorDet[availableEtaToPi0MeasA[i]] , colorDet[availableEtaToPi0MeasA[i]]);
            graphWeightsEtaToPi0A[availableEtaToPi0MeasA[i]]->Draw("p,same,z");
            legendWeightsEtaToPi0->AddEntry(graphWeightsEtaToPi0A[availableEtaToPi0MeasA[i]],nameMeasGlobal[availableEtaToPi0MeasA[i]],"p");
        }
        legendWeightsEtaToPi0->Draw();

        labelWeightsEnergy->Draw();
        TLatex *labelWeightsEtaToPi0         = new TLatex(0.7,0.16,"#eta/#pi^{0}");
        SetStyleTLatex( labelWeightsEtaToPi0, 0.85*textSizeLabelsPixel,4);
        labelWeightsEtaToPi0->SetTextFont(43);
        labelWeightsEtaToPi0->Draw();

//      DrawGammaLines(0.33, 25. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.33, 30. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.33, 30. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.33, 30. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.33, 30. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/EtaToPi0_WeightsA.%s",outputDir.Data(),suffix.Data()));

    // ***************************************************************************************************************
    // ******************************** fitting eta/pi0 **************************************************************
    // ***************************************************************************************************************

    TF1 *fitEtaToPi0 = new TF1("fitEtaToPi0","[0]",3.5,25.);
    fitEtaToPi0->SetParameter(0,0.48);

    TGraphAsymmErrors* comEtaPi0 = (TGraphAsymmErrors*) graphCombEtaToPi0StatA->Clone();
    comEtaPi0->Fit(fitEtaToPi0,"QNRMEX0+","",3.5,25.);
    cout << "\n\n\n\n\n++++++++++++++++++++++++++++++++" << endl;
    cout << fitEtaToPi0->GetParameter(0) << ", +- " << fitEtaToPi0->GetParError(0) << endl;
    cout << "++++++++++++++++++++++++++++++++\n\n\n\n\n" << endl;

    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "Eta/Pi0 - Fit pol0: 3.5 < pT < 25.0" << endl;
    fLog << fitEtaToPi0->GetParameter(0) << ", +- " << fitEtaToPi0->GetParError(0) << endl;

    TGraphAsymmErrors* comEtaPi0Tot = (TGraphAsymmErrors*) graphCombEtaToPi0TotA->Clone();
    comEtaPi0Tot->Fit(fitEtaToPi0,"QNRMEX0+","",3.5,25.);
    cout << "\n\n\n\n\n++++++++++++++++++++++++++++++++" << endl;
    cout << fitEtaToPi0->GetParameter(0) << ", +- " << fitEtaToPi0->GetParError(0) << endl;
    cout << "++++++++++++++++++++++++++++++++\n\n\n\n\n" << endl;

    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "Eta/Pi0 - Fit pol0: 3.5 < pT < 25.0" << endl;
    fLog << fitEtaToPi0->GetParameter(0) << ", +- " << fitEtaToPi0->GetParError(0) << endl;


    // *********************************************************************************************************************
    // ************************************ Visualize relative errors EtaToPi0 ******************************************************
    // *********************************************************************************************************************

    canvasRelSysErr->cd();

    TH2F * histo2DRelSysErrEtaToPi0;
    histo2DRelSysErrEtaToPi0                 = new TH2F("histo2DRelSysErrEtaToPi0","histo2DRelSysErrEtaToPi0",11000,0.33,30.,1000,0,40.5);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErrEtaToPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErrEtaToPi0->GetXaxis()->SetNoExponent(kTRUE);
//    histo2DRelSysErrEtaToPi0->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelSysErrEtaToPi0->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelSysErrEtaToPi0->Draw("copy");

        TLegend* legendRelSysErrEtaToPi0        = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetEtaToPi0A), 0.95, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetEtaToPi0A; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionEtaToPi0[availableEtaToPi0MeasA[i]], markerStyleDet[availableEtaToPi0MeasA[i]], markerSizeDet[availableEtaToPi0MeasA[i]]*0.5,
                                     colorDet[availableEtaToPi0MeasA[i]], colorDet[availableEtaToPi0MeasA[i]]);
            sysErrorRelCollectionEtaToPi0[availableEtaToPi0MeasA[i]]->Draw("p,same,z");
            legendRelSysErrEtaToPi0->AddEntry(sysErrorRelCollectionEtaToPi0[availableEtaToPi0MeasA[i]],nameMeasGlobal[availableEtaToPi0MeasA[i]],"p");
        }
        legendRelSysErrEtaToPi0->Draw();

        labelRelSysErrEnergy->Draw();
        TLatex *labelRelSysErrEtaToPi0       = new TLatex(0.15,0.85,"#eta/#pi^{0}");
        SetStyleTLatex( labelRelSysErrEtaToPi0, 0.85*textSizeLabelsPixel,4);
        labelRelSysErrEtaToPi0->SetTextFont(43);
        labelRelSysErrEtaToPi0->Draw();

    canvasRelSysErr->SaveAs(Form("%s/EtaToPi0_RelSysErr.%s",outputDir.Data(),suffix.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors EtaToPi0 **************************************************
    //  *********************************************************************************************************************

    canvasRelStatErr->cd();

    TH2F * histo2DRelStatErrEtaToPi0;
    histo2DRelStatErrEtaToPi0                = new TH2F("histo2DRelStatErrEtaToPi0","histo2DRelStatErrEtaToPi0",11000,0.33,30.,1000,0,65.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErrEtaToPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErrEtaToPi0->GetXaxis()->SetNoExponent(kTRUE);
//    histo2DRelStatErrEtaToPi0->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelStatErrEtaToPi0->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelStatErrEtaToPi0->Draw("copy");

        TLegend* legendRelStatErrEtaToPi0       = GetAndSetLegend2(0.24, 0.92-(0.035*nMeasSetEtaToPi0A), 0.55, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetEtaToPi0A; i++){
            DrawGammaSetMarker(statErrorRelCollectionEtaToPi0[availableEtaToPi0MeasA[i]], markerStyleDet[availableEtaToPi0MeasA[i]], markerSizeDet[availableEtaToPi0MeasA[i]]*0.5,
                               colorDet[availableEtaToPi0MeasA[i]] , colorDet[availableEtaToPi0MeasA[i]]);
            statErrorRelCollectionEtaToPi0[availableEtaToPi0MeasA[i]]->Draw("p,same,z");
            legendRelStatErrEtaToPi0->AddEntry(statErrorRelCollectionEtaToPi0[availableEtaToPi0MeasA[i]],nameMeasGlobal[availableEtaToPi0MeasA[i]],"p");
        }
        legendRelStatErrEtaToPi0->Draw();

        labelRelStatErrEnergy->Draw();
        TLatex *labelRelStatErrEtaToPi0      = new TLatex(0.75,0.85,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelStatErrEtaToPi0, 0.85*textSizeLabelsPixel,4);
        labelRelStatErrEtaToPi0->SetTextFont(43);
        labelRelStatErrEtaToPi0->Draw();

    canvasRelStatErr->SaveAs(Form("%s/EtaToPi0_RelStatErr.%s",outputDir.Data(),suffix.Data()));


    //  *********************************************************************************************************************
    //  ************************ Visualize relative total errors of different combination methods Eta ***********************
    //  *********************************************************************************************************************
    TGraphAsymmErrors* graphCombEtaToPi0RelStatA       = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0StatA, "relativeStatErrorEtaToPi0_MethodA");
    TGraphAsymmErrors* graphCombEtaToPi0RelSysA        = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0SysA, "relativeSysErrorEtaToPi0_MethodA");
    TGraphAsymmErrors* graphCombEtaToPi0RelTotA        = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0TotA, "relativeTotalErrorEtaToPi0_MethodA");


    canvasRelTotErr->cd();
    TH2F * histo2DRelTotErrEtaToPi0;
    histo2DRelTotErrEtaToPi0                 = new TH2F("histo2DRelTotErrEtaToPi0","histo2DRelTotErrEtaToPi0",11000,0.33,30.,1000,0,65.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrEtaToPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrEtaToPi0->GetXaxis()->SetNoExponent(kTRUE);
//    histo2DRelTotErrEtaToPi0->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRelTotErrEtaToPi0->Draw("copy");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelTotA, markerStyleComb+4, markerSizeComb, kBlue+2 , kBlue+2);
        graphCombEtaToPi0RelTotA->Draw("p,same,z");

        TLegend* legendRelTotErrEtaToPi0     = GetAndSetLegend2(0.24, 0.92-(0.035*1), 0.55, 0.92, 32);
        legendRelTotErrEtaToPi0->AddEntry(graphCombEtaToPi0RelTotA,"All","p");
        legendRelTotErrEtaToPi0->Draw();

        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrEtaToPi0       = new TLatex(0.75,0.85,"#eta/#pi^{0}");
        SetStyleTLatex( labelRelTotErrEtaToPi0, 0.85*textSizeLabelsPixel,4);
        labelRelTotErrEtaToPi0->SetTextFont(43);
        labelRelTotErrEtaToPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/EtaToPi0_RelTotErr.%s",outputDir.Data(),suffix.Data()));

    histo2DRelTotErrEtaToPi0->GetYaxis()->SetRangeUser(0,65.5);
    histo2DRelTotErrEtaToPi0->GetYaxis()->SetTitle("Err (%)");
    histo2DRelTotErrEtaToPi0->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelTotA, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombEtaToPi0RelTotA->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelStatA, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombEtaToPi0RelStatA->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelSysA, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombEtaToPi0RelSysA->SetLineStyle(7);
        graphCombEtaToPi0RelSysA->Draw("l,x0,same,e1");

        legendRelTotErr3->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrEtaToPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/EtaToPi0_RelMethodAdecomp.%s",outputDir.Data(),suffix.Data()));


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

        TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, 0.23, 50. ,1000., -30, 60);
        SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                                  0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 510, 505);
        histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.5,31.5);
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

        TLatex *labelLegendAMass    = new TLatex(0.95,0.06,"a)");
        SetStyleTLatex( labelLegendAMass, textSizeLabelsPixel,4);
        labelLegendAMass->SetTextFont(43);
        labelLegendAMass->Draw();

        TLatex *labelMassPerf       = new TLatex(0.13,0.87,ALICEperfor.Data());
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
        if(plotDate){
          TLatex *labelDate = new TLatex(0.13, 0.6, date.Data());
          SetStyleTLatex( labelDate, textSizeLabelsPixel,4);
          labelDate->SetTextFont(43);
          labelDate->Draw();
        }

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

        TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, 0.23, 50., 1000., 110., 190);
        SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                                  textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 510, 505);
        histo2DAllPi0Mass->GetYaxis()->SetRangeUser(119.5,169.5);
        histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
        histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
        histo2DAllPi0Mass->GetXaxis()->SetNoExponent(kTRUE);
        //histo2DAllPi0Mass->GetXaxis()->SetLabelOffset(-0.015);
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

        TLatex *labelLegendBMass            = new TLatex(0.95,0.2,"b)");
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
        TLatex *textMassPCMEMCAL            = new TLatex(columnsLegendMass2[0],rowsLegendMass2[2],"PCM-EMC");
        SetStyleTLatex( textMassPCMEMCAL, textSizeLabelsPixel,4);
        textMassPCMEMCAL->SetTextFont(43);
        textMassPCMEMCAL->Draw();
        TLatex *textMassEMCAL               = new TLatex(columnsLegendMass2[0],rowsLegendMass2[3],"EMC");
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

        histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.5,31.5);
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
        if(plotDate){
          TLatex *labelDate = new TLatex(0.13, 0.6, date.Data());
          SetStyleTLatex( labelDate, textSizeLabelsPixel,4);
          labelDate->SetTextFont(43);
          labelDate->Draw();
        }

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
//         TLatex *textFWHMPCMEMCAL2           = new TLatex(columnsLegendFWHM[0],rowsLegendFWHM3[2],"PCM-EMC (FWHM/2.36)");
//         SetStyleTLatex( textFWHMPCMEMCAL2, textSizeLabelsPixel,4);
//         textFWHMPCMEMCAL2->SetTextFont(43);
//         textFWHMPCMEMCAL2->Draw();
//         TLatex *textFWHMEMCAL2              = new TLatex(columnsLegendFWHM[0],rowsLegendFWHM3[3],"EMC (FWHM/2.36)");
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
        TLatex *textMassPCMEMCAL2           = new TLatex(columnsLegendMass2[0],rowsLegendMass3[2],"PCM-EMC");
        SetStyleTLatex( textMassPCMEMCAL2, textSizeLabelsPixel,4);
        textMassPCMEMCAL2->SetTextFont(43);
        textMassPCMEMCAL2->Draw();
        TLatex *textMassEMCAL2              = new TLatex(columnsLegendMass2[0],rowsLegendMass3[3],"EMC");
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
    // ******************************************* Mass and width for eta at 8TeV ****************************************
    // **********************************************************************************************************************

    TCanvas* canvasMassWidthEta     = new TCanvas("canvasMassWidthEta","",0,0,1350,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasMassWidthEta,  0.13, 0.02, 0.03, 0.06);

    TPad* padWidthEta               = new TPad("padWidthEta", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthEta, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthEta->Draw();

    TPad* padMassEta                = new TPad("padMassEta", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassEta, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
    padMassEta->Draw();

    TPad* padMassLegend3            = new TPad("padMassLegend3", "", 0.62, 0.83, 0.95, 0.98,-1, -1, -2);
    DrawGammaPadSettings( padMassLegend3, 0., 0., 0., 0.);
    padMassLegend3->SetFillStyle(0);
    padMassLegend3->Draw();

    padWidthEta->cd();
    padWidthEta->SetLogx();


        TH2F * histo2DAllEtaFWHM    = new TH2F("histo2DAllEtaFWHM","histo2DAllEtaFWHM", 20, 0.33, 50. ,1000., -3, 79);
        SetStyleHistoTH2ForGraphs(histo2DAllEtaFWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                                  0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.85,0.28/(textsizeFacWidth*margin), 505, 505);
//      histo2DAllEtaFWHM->GetYaxis()->SetRangeUser(-1.,25.5);
        histo2DAllEtaFWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllEtaFWHM->GetXaxis()->SetNoExponent(kTRUE);
        histo2DAllEtaFWHM->GetYaxis()->SetNdivisions(505);
        histo2DAllEtaFWHM->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllEtaFWHM->GetXaxis()->SetTickLength(0.05);
        histo2DAllEtaFWHM->GetYaxis()->SetTickLength(0.026);
        histo2DAllEtaFWHM->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaFWHM, markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
        graphEMCALEtaFWHM->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaFWHMMC, markerStyleDetMC[2], markerSizeDetMC[2]*0.55, colorDetMC[2] , colorDetMC[2]);
        graphEMCALEtaFWHMMC->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaFWHM, markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
        graphPCMEMCALEtaFWHM->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaFWHMMC, markerStyleDetMC[4], markerSizeDetMC[4]*0.55, colorDetMC[4] , colorDetMC[4]);
        graphPCMEMCALEtaFWHMMC->Draw("p,same,z");

        DrawGammaSetMarker(histoPCMEtaFWHMMeV, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        histoPCMEtaFWHMMeV->Draw("p,same,e");
        DrawGammaSetMarker(histoPCMEtaTrueFWHMMeV, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
        histoPCMEtaTrueFWHMMeV->Draw("p,same,e");

        labelMassPerf->Draw();
        labelLegendAMass->Draw();
        labelMassEnergy->Draw();
        TLatex *labelMassEta                = new TLatex(0.13,0.69,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelMassEta, textSizeLabelsPixel,4);
        labelMassEta->SetTextFont(43);
        labelMassEta->Draw();

        if(plotDate){
          TLatex *labelDate = new TLatex(0.13, 0.6, date.Data());
          SetStyleTLatex( labelDate, textSizeLabelsPixel,4);
          labelDate->SetTextFont(43);
          labelDate->Draw();
        }

        padMassLegend3->cd();
        //****************** first Column **************************************************
        textMassPCM->Draw();
        textMassPCMEMCAL->Draw();
        textMassEMCAL->Draw();
        //****************** second Column *************************************************
        textMassData->Draw();
        textMassMC->Draw();
        markerPCMPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2);
        markerPCMEMCALPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2);
        markerEMCALPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2);
        markerPCMPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[1]+ offsetMarkerYMass2);
        markerPCMEMCALPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[2]+ offsetMarkerYMass2);
        markerEMCALPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[3]+ offsetMarkerYMass2);

    padMassEta->cd();
    padMassEta->SetLogx();

        TH2F * histo2DAllEtaMass    = new TH2F("histo2DAllEtaMass","histo2DAllEtaMass",1000, 0.33, 50., 1000., 505., 574);
        SetStyleHistoTH2ForGraphs(histo2DAllEtaMass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                                  textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 505, 505);
        histo2DAllEtaMass->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllEtaMass->GetXaxis()->SetNoExponent(kTRUE);
        histo2DAllEtaMass->GetYaxis()->SetNdivisions(505);
        histo2DAllEtaMass->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllEtaMass->GetXaxis()->SetTickLength(0.05);
//        histo2DAllEtaMass->GetXaxis()->SetLabelOffset(-0.015);
        histo2DAllEtaMass->DrawCopy();

        DrawGammaSetMarker(histoPCMEtaMass, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        histoPCMEtaMass->Draw("p,same,e");
        DrawGammaSetMarker(histoPCMEtaTrueMass, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
        histoPCMEtaTrueMass->Draw("p,same,e");

        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaMass, markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
        graphEMCALEtaMass->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaMassMC, markerStyleDetMC[2], markerSizeDetMC[2]*0.55, colorDetMC[2] , colorDetMC[2]);
        graphEMCALEtaMassMC->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaMass, markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
        graphPCMEMCALEtaMass->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaMassMC, markerStyleDetMC[4], markerSizeDetMC[4]*0.55, colorDetMC[4] , colorDetMC[4]);
        graphPCMEMCALEtaMassMC->Draw("p,same,z");

        DrawGammaLines(0.33, 50. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,0.3, kGray);

        labelLegendBMass->Draw();

    canvasMassWidthEta->Update();
    canvasMassWidthEta->Print(Form("%s/Eta_MassAndWidth.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************** Cross section for pi0 single measurement 8TeV ************************************
    // **********************************************************************************************************************

    TCanvas* canvasXSectionPi0  = new TCanvas("canvasXSectionPi0","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasXSectionPi0, 0.14, 0.02, 0.02, 0.09);
    canvasXSectionPi0->SetLogx();
    canvasXSectionPi0->SetLogy();

    TH2F * histo2DXSectionPi0;
    histo2DXSectionPi0          = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,0.23,50.,1000,6,9e11);
    SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.45);
    histo2DXSectionPi0->GetXaxis()->SetMoreLogLabels();
    histo2DXSectionPi0->GetXaxis()->SetNoExponent(kTRUE);
//    histo2DXSectionPi0->GetXaxis()->SetLabelOffset(-0.01);
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
//        DrawGammaSetMarkerTGraphAsym(graphEMCALMergedPi0InvXSectionSys, markerStyleDet[9], markerSizeDet[9]*0.75, colorDet[9] , colorDet[9], widthLinesBoxes, kTRUE);
//        graphEMCALMergedPi0InvXSectionSys->Draw("E2same");

        DrawGammaSetMarkerTGraphAsym(graphPHOSPi0InvXSectionStat_WOXErr, markerStyleDet[1], markerSizeDet[1]*0.75, colorDet[1] , colorDet[1]);
        graphPHOSPi0InvXSectionStat_WOXErr->Draw("p,same,e");
//      DrawGammaSetMarker(histoPCMPi0InvXSectionStat, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
//      histoPCMPi0InvXSectionStat->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvXSectionStat_WOXErr,markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
        graphPCMPi0InvXSectionStat_WOXErr->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0InvXSectionStat_WOXErr, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
        graphEMCALPi0InvXSectionStat_WOXErr->Draw("p,same,z");
//      DrawGammaSetMarker(histoEMCALPi0InvXSectionStat, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
//      histoEMCALPi0InvXSectionStat->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0InvXSectionStat_WOXErr, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4]);
        graphPCMEMCALPi0InvXSectionStat_WOXErr->Draw("p,same,z");
//        DrawGammaSetMarkerTGraphAsym(graphEMCALMergedPi0InvXSectionStat, markerStyleDet[9], markerSizeDet[9]*0.75, colorDet[9] , colorDet[9]);
//        graphEMCALMergedPi0InvXSectionStat->Draw("p,same,z");


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
        legendXSectionPi0->AddEntry(graphEMCALPi0InvXSectionSys,"EMC","fp");
        legendXSectionPi0->AddEntry(graphPCMEMCALPi0InvXSectionSys,"PCM-EMC","fp");
//        legendXSectionPi0->AddEntry(graphEMCALMergedPi0InvXSectionSys,"EMC merged","fp");
        legendXSectionPi0->Draw();

    canvasXSectionPi0->SaveAs(Form("%s/Pi0_InvXSectionCompAllSystems.%s",outputDir.Data(),suffix.Data()));

    canvasXSectionPi0->cd();
    histo2DXSectionPi0->Draw("copy");

        graphPHOSPi0InvXSectionSys->Draw("E2same");
        graphPCMPi0InvXSectionSys->Draw("E2same");
        graphEMCALPi0InvXSectionSys->Draw("E2same");
        graphPCMEMCALPi0InvXSectionSys->Draw("E2same");
//        graphEMCALMergedPi0InvXSectionSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionSysA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphCombPi0InvXSectionSysA->Draw("E2same");

        graphPHOSPi0InvXSectionStat_WOXErr->Draw("p,same,e");
        graphPCMPi0InvXSectionStat_WOXErr->Draw("p,same,z");
        graphEMCALPi0InvXSectionStat_WOXErr->Draw("p,same,z");
        graphPCMEMCALPi0InvXSectionStat_WOXErr->Draw("p,same,z");
//        graphEMCALMergedPi0InvXSectionStat->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatA_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombPi0InvXSectionStatA_WOXErr->Draw("p,same,z");

        labelEnergyXSectionPi0->Draw();
        labelDetSysXSectionPi0->Draw();

        legendXSectionPi0->AddEntry(graphCombPi0InvXSectionSysA,"comb","fp");
        legendXSectionPi0->Draw();

//      DrawGammaSetMarkerTGraphAsym(graphChargedHadronsStatPP, markerStyleDet[1], markerSizeDet[1]*0.75, kBlack , kBlack, widthLinesBoxes);
//      graphChargedHadronsStatPP->Draw("E2same");

   canvasXSectionPi0->SaveAs(Form("%s/Pi0_InvXSectionCompAllSystems_Comb.%s",outputDir.Data(),suffix.Data()));

    TCanvas* canvasXSectionEta      = new TCanvas("canvasXSectionEta","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasXSectionEta, 0.14, 0.02, 0.02, 0.09);
    canvasXSectionEta->SetLogx();
    canvasXSectionEta->SetLogy();

    TH2F * histo2DXSectionEta;
    histo2DXSectionEta              = new TH2F("histo2DXSectionEta","histo2DXSectionEta",11000,0.23, 50.,1000,6,2e10);
    SetStyleHistoTH2ForGraphs(histo2DXSectionEta, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 1.,1.45);
    histo2DXSectionEta->GetXaxis()->SetMoreLogLabels();
    histo2DXSectionEta->GetXaxis()->SetNoExponent(kTRUE);
//    histo2DXSectionEta->GetXaxis()->SetLabelOffset(-0.01);
    histo2DXSectionEta->Draw("copy");


        DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvXSectionSys, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0], widthLinesBoxes, kTRUE);
        graphPCMEtaInvXSectionSys->Draw("E2same");

        //      graphPCMEMCALEtaInvXSectionSys->RemovePoint(0);
//      graphPCMEMCALEtaInvXSectionSys->RemovePoint(graphPCMEMCALEtaInvXSectionSys->GetN()-1);
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaInvXSectionSys, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4], widthLinesBoxes, kTRUE);
        graphPCMEMCALEtaInvXSectionSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaInvXSectionSys, markerStyleDet[2], markerSizeDet[4]*0.75, colorDet[2] , colorDet[2], widthLinesBoxes, kTRUE);
        graphEMCALEtaInvXSectionSys->Draw("E2same");

        DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvXSectionStat_WOXErr, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
        graphPCMEtaInvXSectionStat_WOXErr->Draw("p,same,e");
//      DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvXSectionStat,markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
//      graphPCMPi0InvXSectionStat->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaInvXSectionStat_WOXErr, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4]);
        graphPCMEMCALEtaInvXSectionStat_WOXErr->Draw("p,same,e");

        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaInvXSectionStat_WOXErr, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
        graphEMCALEtaInvXSectionStat_WOXErr->Draw("p,same,e");


        TLatex *labelEnergyXSectionEta      = new TLatex(0.64,0.92,collisionSystem8TeV.Data());
        SetStyleTLatex( labelEnergyXSectionEta, 0.035,4);
        labelEnergyXSectionEta->Draw();
        TLatex *labelDetSysXSectionEta      = new TLatex(0.64,0.88,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionEta, 0.035,4);
        labelDetSysXSectionEta->Draw();

        TLegend* legendXSectionEta          = new TLegend(0.62,0.72,0.9,0.86);
        legendXSectionEta->SetFillColor(0);
        legendXSectionEta->SetLineColor(0);
        legendXSectionEta->SetTextFont(42);
        legendXSectionEta->SetTextSize(0.035);
        legendXSectionEta->AddEntry(graphPCMEtaInvXSectionSys,"PCM","fp");
//      legendXSectionEta->AddEntry(graphPHOSPi0InvXSectionSys,"PHOS","fp");
        legendXSectionEta->AddEntry(graphEMCALEtaInvXSectionSys,"EMC","fp");
        legendXSectionEta->AddEntry(graphPCMEMCALEtaInvXSectionSys,"PCM-EMC","fp");
        legendXSectionEta->Draw();

    canvasXSectionEta->SaveAs(Form("%s/Eta_InvXSectionCompAllSystems.%s",outputDir.Data(),suffix.Data()));
    histo2DXSectionEta->Draw("copy");

        graphPCMEtaInvXSectionSys->Draw("E2same");
        graphEMCALEtaInvXSectionSys->Draw("E2same");
        graphPCMEMCALEtaInvXSectionSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSysA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphCombEtaInvXSectionSysA->Draw("E2same");

        graphPCMEtaInvXSectionStat_WOXErr->Draw("p,same,e");
        graphEMCALEtaInvXSectionStat_WOXErr->Draw("p,same,e");
        graphPCMEMCALEtaInvXSectionStat_WOXErr->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStatA_WOXErr, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombEtaInvXSectionStatA_WOXErr->Draw("p,same,z");


        labelEnergyXSectionEta->Draw();
        labelDetSysXSectionEta->Draw();

        legendXSectionEta->AddEntry(graphCombEtaInvXSectionSysA,"comb","fp");
        legendXSectionEta->Draw();

//      DrawGammaSetMarkerTGraphAsym(graphChargedHadronsStatPP, markerStyleDet[1], markerSizeDet[1]*0.75, kBlack , kBlack, widthLinesBoxes);
//      graphChargedHadronsStatPP->Draw("E2same");

   canvasXSectionEta->SaveAs(Form("%s/Eta_InvXSectionCompAllSystems_Comb.%s",outputDir.Data(),suffix.Data()));


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
        histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, 0.23, 50, 1000, 2e-6, 2e-0 );
        SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec}"),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);
        histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
        histo2DAccEff->GetXaxis()->SetNoExponent(kTRUE);
        //histo2DAccEff->GetXaxis()->SetLabelOffset(-0.01);
        histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAccEff->DrawCopy();

        DrawGammaSetMarker(histoPCMPi0AccTimesEff, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        histoPCMPi0AccTimesEff->Draw("p,same,e");

        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0AccTimesEff, markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
        graphPCMEMCALPi0AccTimesEff->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0AccTimesEff, markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
        graphEMCALPi0AccTimesEff->Draw("p,same,z");

//        DrawGammaSetMarkerTGraphAsym(graphEMCALMergedPi0AccTimesEffDivPur, markerStyleDet[9], markerSizeDet[9]*0.55, colorDet[9] , colorDet[9]);
//        graphEMCALMergedPi0AccTimesEffDivPur->Draw("p,same,e");

        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.62, 0.13, 0.9, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffiAccPi0->AddEntry(histoPCMPi0AccTimesEff,"PCM","p");
        legendEffiAccPi0->AddEntry(graphPCMEMCALPi0AccTimesEff,"PCM-EMC","p");
        legendEffiAccPi0->AddEntry(graphEMCALPi0AccTimesEff,"EMC","p");
//        legendEffiAccPi0->AddEntry(graphEMCALMergedPi0AccTimesEffDivPur,"EMC, merged","p");
        legendEffiAccPi0->Draw();

        TLatex *labelPerfEffi               = new TLatex(0.15,0.92,"ALICE simulation");
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
//        graphEMCALMergedPi0AccTimesEffDivPur->Draw("p,same,e");

        TLegend* legendEffiAccPi02          = GetAndSetLegend2(0.62, 0.13, 0.9, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffiAccPi02->AddEntry(histoPCMPi0AccTimesEff,"PCM","p");
        legendEffiAccPi02->AddEntry(graphPCMEMCALPi0AccTimesEff,"PCM-EMC","p");
        legendEffiAccPi02->AddEntry(graphEMCALPi0AccTimesEff,"EMC","p");
//        legendEffiAccPi02->AddEntry(graphEMCALMergedPi0AccTimesEffDivPur,"EMC, merged","p");
        legendEffiAccPi02->AddEntry(histoPHOSPi0AccTimesEff,"PHOS","p");
        legendEffiAccPi02->Draw();

        labelPerfEffi->Draw();
        labelEnergyEffi->Draw();
        labelDetSysEffiPi0->Draw();

    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff_incPHOS.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************** Acceptance * Efficiency for eta single measurement 8TeV **************************
    // **********************************************************************************************************************
    canvasAcceptanceTimesEff->cd();

        TH2F * histo2DAccEffEta;
        histo2DAccEffEta                = new TH2F("histo2DAccEffEta", "histo2DAccEffEta",1000, 0.33,  50, 1000, 5e-4, 1.5e-0 );
        SetStyleHistoTH2ForGraphs( histo2DAccEffEta, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec}"),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.);
        histo2DAccEffEta->GetYaxis()->SetLabelOffset(0.001);
//        histo2DAccEffEta->GetXaxis()->SetLabelOffset(-0.01);
        histo2DAccEffEta->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAccEffEta->GetXaxis()->SetNoExponent(kTRUE);
        histo2DAccEffEta->DrawCopy();

        DrawGammaSetMarker(histoPCMEtaAccTimesEff, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        histoPCMEtaAccTimesEff->Draw("p,same,e");

        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaAccTimesEff, markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
        graphPCMEMCALEtaAccTimesEff->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaAccTimesEff, markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
        graphEMCALEtaAccTimesEff->Draw("p,same,z");

        TLegend* legendEffiAccEta           = GetAndSetLegend2(0.62, 0.13, 0.9, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffiAccEta->AddEntry(histoPCMEtaAccTimesEff,"PCM","p");
        legendEffiAccEta->AddEntry(graphPCMEMCALEtaAccTimesEff,"PCM-EMC","p");
        legendEffiAccEta->AddEntry(graphEMCALEtaAccTimesEff,"EMC","p");
        legendEffiAccEta->Draw();

        labelPerfEffi->Draw();
        labelEnergyEffi->Draw();
        TLatex *labelDetSysEffiEta          = new TLatex(0.15,0.82,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysEffiEta, textSizeLabelsRel,4);
        labelDetSysEffiEta->Draw();


    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Eta_AcceptanceTimesEff.%s",outputDir.Data(),suffix.Data()));

//    TCanvas* canvasTriggerEff       = new TCanvas("canvasTriggerEff", "", 200, 10, 1200, 1100);  // gives the page size
//    DrawGammaCanvasSettings( canvasTriggerEff,  0.09, 0.01, 0.015, 0.095);
////  canvasTriggerEff->SetLogx(1);

//        TH2F * histo2DTriggerEffPi0;
//        histo2DTriggerEffPi0                = new TH2F("histo2DTriggerEffPi0", "histo2DTriggerEffPi0",1000, 0.,  21, 1000, 0, 1.15 );
//        SetStyleHistoTH2ForGraphs( histo2DTriggerEffPi0, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{Trigger}",
//                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 0.95);//(#times #epsilon_{pur})
//            histo2DTriggerEffPi0->GetXaxis()->SetMoreLogLabels(kTRUE);
//        histo2DTriggerEffPi0->DrawCopy();

//        DrawGammaSetMarker(histoPCMEMCALPi0TriggerEff[0], markerTriggMC[2], sizeTrigg[2]*1.5, colorTriggShade[2] , colorTriggShade[2]);
//        histoPCMEMCALPi0TriggerEff[0]->Draw("p,same,e");
//        DrawGammaSetMarker(histoPCMEMCALPi0TriggerEff[1], markerTriggMC[3], sizeTrigg[3]*1.5, colorTriggShade[3] , colorTriggShade[3]);
//        histoPCMEMCALPi0TriggerEff[1]->Draw("p,same,e");
//        DrawGammaSetMarker(histoPCMEMCALPi0TriggerEff[3], markerTriggMC[5], sizeTrigg[5]*1.5, colorTriggShade[5] , colorTriggShade[5]);
//        histoPCMEMCALPi0TriggerEff[3]->Draw("p,same,e");
//        DrawGammaSetMarker(histoEMCALPi0TriggerEff[0], markerTrigg[2], sizeTrigg[2]*1.5, colorTrigg[2] , colorTrigg[2]);
//        histoEMCALPi0TriggerEff[0]->Draw("p,same,e");
//        DrawGammaSetMarker(histoEMCALPi0TriggerEff[1], markerTrigg[3], sizeTrigg[3]*1.5, colorTrigg[3] , colorTrigg[3]);
//        histoEMCALPi0TriggerEff[1]->Draw("p,same,e");
//        DrawGammaSetMarker(histoEMCALPi0TriggerEff[3], markerTrigg[5], sizeTrigg[5]*1.5, colorTrigg[5] , colorTrigg[5]);
//        histoEMCALPi0TriggerEff[3]->Draw("p,same,e");


//        TLatex *labelPerfTriggEff               = new TLatex(0.6,0.425,"ALICE simulation");
//        SetStyleTLatex( labelPerfTriggEff, textSizeLabelsRel,4);
//        labelPerfTriggEff->Draw();
//        TLatex *labelEnergyTriggEff             = new TLatex(0.6,0.378,collisionSystem8TeV.Data());
//        SetStyleTLatex( labelEnergyTriggEff, textSizeLabelsRel,4);
//        labelEnergyTriggEff->Draw();
//        TLatex *labelDetSysTriggEffPi0          = new TLatex(0.6,0.33,"#pi^{0} #rightarrow #gamma#gamma");
//        SetStyleTLatex( labelDetSysTriggEffPi0, textSizeLabelsRel,4);
//        labelDetSysTriggEffPi0->Draw();

//        TLegend* legendTriggEffPi0              = GetAndSetLegend2(0.415, 0.13, 0.9, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel,3);
//        legendTriggEffPi0->AddEntry((TObject*)0, Form("%s    ",nameTrigger[2].Data()),"");
//        legendTriggEffPi0->AddEntry(histoPCMEMCALPi0TriggerEff[0],"       ","p");
//        legendTriggEffPi0->AddEntry(histoEMCALPi0TriggerEff[0],"","p");
//        legendTriggEffPi0->AddEntry((TObject*)0, nameTrigger[3].Data(),"");
//        legendTriggEffPi0->AddEntry(histoPCMEMCALPi0TriggerEff[1],"    ","p");
//        legendTriggEffPi0->AddEntry(histoEMCALPi0TriggerEff[1],"","p");
//        legendTriggEffPi0->AddEntry((TObject*)0, nameTrigger[5].Data(),"");
//        legendTriggEffPi0->AddEntry(histoPCMEMCALPi0TriggerEff[3],"    ","p");
//        legendTriggEffPi0->AddEntry(histoEMCALPi0TriggerEff[3],"","p");
//        legendTriggEffPi0->Draw();

//        TLatex *labelPCMEMCALTriggEff           = new TLatex(0.6,0.275,"PCM-EMC");
//        SetStyleTLatex( labelPCMEMCALTriggEff, textSizeLabelsRel,4);
//        labelPCMEMCALTriggEff->Draw();
//        TLatex *labelEMCALTriggEff              = new TLatex(0.84,0.275,"EMC");
//        SetStyleTLatex( labelEMCALTriggEff, textSizeLabelsRel,4);
//        labelEMCALTriggEff->Draw();

//    canvasTriggerEff->Update();
//    canvasTriggerEff->Print(Form("%s/Pi0_TriggerEff.%s",outputDir.Data(),suffix.Data()));

//    canvasTriggerEff->cd();
//    histo2DTriggerEffPi0->DrawCopy();
//        DrawGammaSetMarker(histoPCMEMCALEtaTriggerEff[0], markerTriggMC[2], sizeTrigg[2]*1.5, colorTriggShade[2] , colorTriggShade[2]);
//        histoPCMEMCALEtaTriggerEff[0]->Draw("p,same,e");
//        DrawGammaSetMarker(histoPCMEMCALEtaTriggerEff[1], markerTriggMC[3], sizeTrigg[3]*1.5, colorTriggShade[3] , colorTriggShade[3]);
//        histoPCMEMCALEtaTriggerEff[1]->Draw("p,same,e");
//        DrawGammaSetMarker(histoPCMEMCALEtaTriggerEff[3], markerTriggMC[5], sizeTrigg[5]*1.5, colorTriggShade[5] , colorTriggShade[5]);
//        histoPCMEMCALEtaTriggerEff[3]->Draw("p,same,e");
//        DrawGammaSetMarker(histoEMCALEtaTriggerEff[0], markerTrigg[2], sizeTrigg[2]*1.5, colorTrigg[2] , colorTrigg[2]);
//        histoEMCALEtaTriggerEff[0]->Draw("p,same,e");
//        DrawGammaSetMarker(histoEMCALEtaTriggerEff[1], markerTrigg[3], sizeTrigg[3]*1.5, colorTrigg[3] , colorTrigg[3]);
//        histoEMCALEtaTriggerEff[1]->Draw("p,same,e");
//        DrawGammaSetMarker(histoEMCALEtaTriggerEff[3], markerTrigg[5], sizeTrigg[5]*1.5, colorTrigg[5] , colorTrigg[5]);
//        histoEMCALEtaTriggerEff[3]->Draw("p,same,e");

//        labelPerfTriggEff->Draw();
//        labelEnergyTriggEff->Draw();
//        TLatex *labelDetSysTriggEffEta          = new TLatex(0.6,0.33,"#eta #rightarrow #gamma#gamma");
//        SetStyleTLatex( labelDetSysTriggEffEta, 0.85*textSizeLabelsRel,4);
//        labelDetSysTriggEffEta->Draw();
//        labelPCMEMCALTriggEff->Draw();
//        labelEMCALTriggEff->Draw();
//        legendTriggEffPi0->Draw();


//    canvasTriggerEff->Update();
//    canvasTriggerEff->Print(Form("%s/Eta_TriggerEff.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************** effective secondary correction drawing for different methods ************************
    // **********************************************************************************************************************
    TCanvas* canvasEffectiveSecCorr       = new TCanvas("canvasEffectiveSecCorr", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasEffectiveSecCorr,  0.1, 0.01, 0.04, 0.095);
    canvasEffectiveSecCorr->SetLogx(1);

        TH1F * histo1DEffSecCorr;
        histo1DEffSecCorr                = new TH1F("histo1DEffSecCorr", "histo1DEffSecCorr",1000, 0.23,  31);
        SetStyleHistoTH1ForGraphs( histo1DEffSecCorr, "#it{p}_{T} (GeV/#it{c})","R_{K^{0}_{s}}",
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
        histo1DEffSecCorr->GetYaxis()->SetRangeUser(0, 10 );
        histo1DEffSecCorr->GetXaxis()->SetLabelOffset(-0.01);
        histo1DEffSecCorr->GetXaxis()->SetMoreLogLabels(kTRUE);

        for (Int_t k = 0; k < 4; k++){
            Bool_t plotCorr     = kFALSE;
            Int_t nCorrAvail    = 0;
            for (Int_t i = 0; i < 11; i++){
                if (haveEffSecCorr[k][i]){
                    nCorrAvail++;
                    plotCorr    = kTRUE;
                }
            }
            TLegend* legendEffSecCorrPi0           = GetAndSetLegend2(0.65, 0.925-(nCorrAvail*textSizeLabelsRel), 0.93, 0.925,textSizeLabelsPixel);
            if (plotCorr){
                histo1DEffSecCorr->GetYaxis()->SetTitle(Form("R_{%s}",nameSecPi0SourceLabel[k].Data()));
                histo1DEffSecCorr->GetYaxis()->SetRangeUser(0, maxSecCorr[k]);
                histo1DEffSecCorr->DrawCopy();
                for (Int_t i = 0; i < 11; i++){
                    if (haveEffSecCorr[k][i]){
                        DrawGammaSetMarkerTGraphAsym(graphPi0EffSecCorrFromX[k][i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
                        graphPi0EffSecCorrFromX[k][i]->Draw("p,same,z");
                        legendEffSecCorrPi0->AddEntry(graphPi0EffSecCorrFromX[k][i],nameMeasGlobal[i],"p");
                    }
                }
                legendEffSecCorrPi0->Draw();
                Double_t xAlign = 0.15;
                if(k==2) xAlign += 0.05;
                TLatex *labelPerfSecCorr               = new TLatex(xAlign,0.90,"ALICE performance");
                SetStyleTLatex( labelPerfSecCorr, textSizeLabelsRel,4);
                labelPerfSecCorr->Draw();
                TLatex *labelEnergySecCorr             = new TLatex(xAlign,0.85,collisionSystem8TeV.Data());
                SetStyleTLatex( labelEnergySecCorr, textSizeLabelsRel,4);
                labelEnergySecCorr->Draw();
                TLatex *labelDetSysSecCorrPi0          = new TLatex(xAlign,0.80,"#pi^{0} #rightarrow #gamma#gamma");
                SetStyleTLatex( labelDetSysSecCorrPi0, textSizeLabelsRel,4);
                labelDetSysSecCorrPi0->Draw();

                canvasEffectiveSecCorr->Update();
                canvasEffectiveSecCorr->Print(Form("%s/Pi0_EffectiveSecCorr_%s.%s",outputDir.Data(), nameSecPi0SourceRead[k].Data() , suffix.Data()));
            }
        }

    delete canvasEffectiveSecCorr;

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

        TH2F * ratio2DTheoryPP       = new TH2F("ratio2DTheoryPP","ratio2DTheoryPP",1000,0.23,50.,1000,0.6,2.2);
        SetStyleHistoTH2ForGraphs(ratio2DTheoryPP, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{TCM fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
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

        DrawGammaSetMarkerTGraphErr(graphRatioPythia8_4CToFit, 0, 0, 807 , 807, widthLinesBoxes, kTRUE, 807);
        graphRatioPythia8_4CToFit->Draw("3,same");
        DrawGammaSetMarker(histoRatioPythia8_4CToFit, 24, 1.5, 807 , 807);
        histoRatioPythia8_4CToFit->SetLineWidth(widthCommonFit);
        histoRatioPythia8_4CToFit->Draw("same,hist,l");

        DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFit, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
        graphRatioPythia8ToFit->Draw("3,same");
        DrawGammaSetMarker(histoRatioPythia8ToFit, 24, 1.5, kRed+2 , kRed+2);
        histoRatioPythia8ToFit->SetLineWidth(widthCommonFit);
        histoRatioPythia8ToFit->Draw("same,hist,l");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatA_WOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA->Draw("2,same");
        graphRatioPi0CombCombFitStatA_WOXErr->Draw("p,same");

        TBox* boxErrorSigmaRatio = CreateBoxConv(kGray+2, 0.3, 1.-(0.026 ), 0.35, 1.+(0.026));
        boxErrorSigmaRatio->SetLineWidth(8);
        boxErrorSigmaRatio->Draw();
        DrawGammaLines(0.23, 50.,1., 1.,0.5,kGray+2);

        TLegend* legendRatioTheorypp_3Parted = GetAndSetLegend2(0.2,0.76,0.45,0.96, 0.85* textSizeLabelsPixel);
    //    legendRatioTheorypp_3Parted->SetNColumns(2);
        legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0CombCombFitSysA,"Data","pf");
//        legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0CombNLOMuHalf, "NLO, DSS07 #mu = 0.5 #it{p}_{T}", "l");
//        legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0CombNLOMuOne,  "NLO, DSS07 #mu = #it{p}_{T}", "l");
//        legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0CombNLOMuTwo,  "NLO, DSS07 #mu = 2 #it{p}_{T}", "l");
        legendRatioTheorypp_3Parted->AddEntry(histoRatioPythia8ToFit,  "PYTHIA 8.2, Monash 2013", "l");
        legendRatioTheorypp_3Parted->AddEntry(histoRatioPythia8_4CToFit,  "PYTHIA 8.2, Tune 4C", "l");
        legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0DSS14,  "NLO, PDF:MSTW08 - FF:DSS14", "f");
        legendRatioTheorypp_3Parted->Draw();

        TLegend* legendRatioTheoryNormUnc = GetAndSetLegend2(0.34,0.91,0.59,0.96, 0.85* textSizeLabelsPixel);
        legendRatioTheoryNormUnc->AddEntry(boxErrorSigmaRatio,"norm. unc. 2.6%","l");
        legendRatioTheoryNormUnc->Draw();

        TLatex *labelRatioTheoryPPP   = new TLatex(0.268,0.73,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}");
        SetStyleTLatex( labelRatioTheoryPPP, 0.85*textsizeLabelsPP,4);
        labelRatioTheoryPPP->Draw();

        TLatex *labelRatioTheoryPP   = new TLatex(0.76,0.925,collisionSystem8TeV.Data());
        SetStyleTLatex( labelRatioTheoryPP, 0.85*textsizeLabelsPP,4);
        labelRatioTheoryPP->Draw();
        TLatex *labelRatioTheoryPP1P = new TLatex(0.863,0.875,"ALICE");
        SetStyleTLatex( labelRatioTheoryPP1P, 0.85*textsizeLabelsPP,4);
        labelRatioTheoryPP1P->Draw();
        TLatex *labelRatioTheoryPP2P= new TLatex(0.843,0.83,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioTheoryPP2P, 0.85*textsizeLabelsPP,4);
        labelRatioTheoryPP2P->Draw();


    canvasRatioPP->Update();
    canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PP.%s",outputDir.Data(),suffix.Data()));


    TH2F * ratio2DTheoryPP2       = new TH2F("ratio2DTheoryPP2","ratio2DTheoryPP2",1000,0.23,50.,1000,0.5,3.6);
    SetStyleHistoTH2ForGraphs(ratio2DTheoryPP2, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{TCM fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
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

    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
    graphRatioPi0CombCombFitStatA_WOXErr->SetLineWidth(widthLinesBoxes);
    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
    graphRatioPi0CombCombFitSysA->SetLineWidth(0);
    graphRatioPi0CombCombFitSysA->Draw("2,same");
    graphRatioPi0CombCombFitStatA_WOXErr->Draw("p,same");

    boxErrorSigmaRatio->Draw();
    DrawGammaLines(0.23, 50.,1., 1.,0.5,kGray+2);

    TLegend* legendRatioTheorypp_3Parted2= GetAndSetLegend2(0.15,0.65,0.4,0.96, 0.85* textSizeLabelsPixel);
    legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0CombCombFitSysA,"Data","pf");
    legendRatioTheorypp_3Parted2->AddEntry((TObject*)0,"NLO, PDF:CTEQ6M5 - FF:DSS07", "");
    legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0CombNLOMuHalf, "#mu = 0.5 #it{p}_{T}", "l");
    legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0CombNLOMuOne,  "#mu = #it{p}_{T}", "l");
    legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0CombNLOMuTwo,  "#mu = 2 #it{p}_{T}", "l");
    //legendRatioTheorypp_3Parted2->AddEntry(graphRatioPi0DSS14,  "NLO, PDF:MSTW08 - FF:DSS14", "f");
    legendRatioTheorypp_3Parted2->Draw();

    TLegend* legendRatioTheoryNormUnc2 = GetAndSetLegend2(0.34,0.902,0.59,0.952, 0.85* textSizeLabelsPixel);
    legendRatioTheoryNormUnc2->AddEntry(boxErrorSigmaRatio,"norm. unc. 2.6%","l");
    legendRatioTheoryNormUnc2->Draw();

    TLegend* legendRatioTheorypp_3Parted22= GetAndSetLegend2(0.15,0.14,0.4,0.18, 0.85* textSizeLabelsPixel);
    legendRatioTheorypp_3Parted22->AddEntry(graphRatioPi0DSS14,  "NLO, PDF:MSTW08 - FF:DSS14, 0.5#it{p}_{T} < #mu < 2#it{p}_{T}", "f");
    legendRatioTheorypp_3Parted22->Draw();

//    TLatex *labelRatioTheoryPPP2   = new TLatex(0.218,0.62,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}");
//    SetStyleTLatex( labelRatioTheoryPPP2, 0.85*textsizeLabelsPP,4);
//    labelRatioTheoryPPP2->Draw();

    labelRatioTheoryPP->Draw();
    labelRatioTheoryPP1P->Draw();
    labelRatioTheoryPP2P->Draw();

    canvasRatioPP->Update();
    canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PP2.%s",outputDir.Data(),suffix.Data()));


    TH2F * ratio2DTheoryPP3       = new TH2F("ratio2DTheoryPP3","ratio2DTheoryPP3",1000,0.23,50.,1000,0.4,3.0);
    SetStyleHistoTH2ForGraphs(ratio2DTheoryPP3, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{TCM fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
    ratio2DTheoryPP3->GetYaxis()->SetMoreLogLabels(kTRUE);
    ratio2DTheoryPP3->GetYaxis()->SetNdivisions(505);
    ratio2DTheoryPP3->GetYaxis()->SetNoExponent(kTRUE);
    ratio2DTheoryPP3->GetXaxis()->SetMoreLogLabels(kTRUE);
    ratio2DTheoryPP3->GetXaxis()->SetNoExponent(kTRUE);
//      ratio2DTheoryPP3->GetYaxis()->SetLabelOffset(0.01);
    ratio2DTheoryPP3->GetXaxis()->SetLabelFont(42);
    ratio2DTheoryPP3->GetYaxis()->SetLabelFont(42);
//      ratio2DTheoryPP3->GetXaxis()->SetLabelOffset(-0.01);
//    ratio2DTheoryPP3->GetXaxis()->SetTickLength(0.07);
    ratio2DTheoryPP3->DrawCopy();

    graphRatioPOWHEG8TeV->SetLineWidth(widthCommonFit);
    graphRatioPOWHEG8TeV->SetLineColor(colorNLO);
    graphRatioPOWHEG8TeV->SetLineStyle(1);
    graphRatioPOWHEG8TeV->SetFillStyle(1001);
    graphRatioPOWHEG8TeV->SetFillColor(colorNLO);
    graphRatioPOWHEG8TeV->Draw("same,e4");

    DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFit, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
    graphRatioPythia8ToFit->Draw("3,same");
    DrawGammaSetMarker(histoRatioPythia8ToFit, 24, 1.5, kRed+2 , kRed+2);
    histoRatioPythia8ToFit->SetLineWidth(widthCommonFit);
    histoRatioPythia8ToFit->Draw("same,hist,l");

    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
    graphRatioPi0CombCombFitStatA_WOXErr->SetLineWidth(widthLinesBoxes);
    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
    graphRatioPi0CombCombFitSysA->SetLineWidth(0);
    graphRatioPi0CombCombFitSysA->Draw("2,same");
    graphRatioPi0CombCombFitStatA_WOXErr->Draw("p,same");

    boxErrorSigmaRatio->Draw();
    DrawGammaLines(0.23, 50.,1., 1.,0.5,kGray+2);

    TLegend* legendRatioTheorypp_3Parted3= GetAndSetLegend2(0.15,0.7,0.4,0.96, 0.85* textSizeLabelsPixel);
    legendRatioTheorypp_3Parted3->AddEntry(graphRatioPi0CombCombFitSysA,"Data","pf");
    legendRatioTheorypp_3Parted3->AddEntry(histoRatioPythia8ToFit,  "PYTHIA 8.2, Monash 2013", "l");
    legendRatioTheorypp_3Parted3->AddEntry(graphRatioPOWHEG8TeV,  "POWHEG + PYTHIA 8", "f");
    legendRatioTheorypp_3Parted3->Draw();

    TLegend* legendRatioTheoryNormUnc3 = GetAndSetLegend2(0.34,0.885,0.59,0.935, 0.85* textSizeLabelsPixel);
    legendRatioTheoryNormUnc3->AddEntry(boxErrorSigmaRatio,"norm. unc. 2.6%","l");
    legendRatioTheoryNormUnc3->Draw();

    TLatex *labelRatioTheoryPPP22   = new TLatex(0.218,0.67,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}");
    SetStyleTLatex( labelRatioTheoryPPP22, 0.85*textsizeLabelsPP,4);
    labelRatioTheoryPPP22->Draw();

    labelRatioTheoryPP->Draw();
    labelRatioTheoryPP1P->Draw();
    labelRatioTheoryPP2P->Draw();

    canvasRatioPP->Update();
    canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PP3.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************************* Comparison to theory calculations Eta ************************************
    // **********************************************************************************************************************
    canvasRatioPP->cd();
    canvasRatioPP->SetLogx();

        TH2F * ratio2DTheoryPPEta       = new TH2F("ratio2DTheoryPPEta","ratio2DTheoryPPEta",1000,0.23,50.,1000,0.2,3.9);
        SetStyleHistoTH2ForGraphs(ratio2DTheoryPPEta, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{TCM fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                                  0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
        ratio2DTheoryPPEta->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio2DTheoryPPEta->GetYaxis()->SetNdivisions(505);
        ratio2DTheoryPPEta->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DTheoryPPEta->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DTheoryPPEta->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DTheoryPPEta->GetXaxis()->SetLabelFont(42);
        ratio2DTheoryPPEta->GetYaxis()->SetLabelFont(42);
        ratio2DTheoryPPEta->DrawCopy();

        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
        graphRatioEtaCombNLOMuHalf->Draw("same,c");
        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
        graphRatioEtaCombNLOMuOne->Draw("same,c");
        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
        graphRatioEtaCombNLOMuTwo->Draw("same,c");

        DrawGammaSetMarkerTGraphErr(graphRatioPythia8_4CToFitEta, 0, 0, 807 , 807, widthLinesBoxes, kTRUE, 807);
        graphRatioPythia8_4CToFitEta->Draw("3,same");
        DrawGammaSetMarker(histoRatioPythia8_4CToFitEta, 24, 1.5, 807 , 807);
        histoRatioPythia8_4CToFitEta->SetLineWidth(widthCommonFit);
        histoRatioPythia8_4CToFitEta->Draw("same,hist,l");

        DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFitEta, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
        graphRatioPythia8ToFitEta->Draw("3,same");
        DrawGammaSetMarker(histoRatioPythia8ToFitEta, 24, 1.5, kRed+2 , kRed+2);
        histoRatioPythia8ToFitEta->SetLineWidth(widthCommonFit);
        histoRatioPythia8ToFitEta->Draw("same,hist,l");

        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatA_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioEtaCombCombFitStatA_WOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioEtaCombCombFitSysA->SetLineWidth(0);
        graphRatioEtaCombCombFitSysA->Draw("2,same");
        graphRatioEtaCombCombFitStatA_WOXErr->Draw("p,same");

        boxErrorSigmaRatio->Draw();
        DrawGammaLines(0.23, 50. , 1., 1.,0.5, kGray+2);

        TLegend* legendRatioTheoryppEta_3Parted= GetAndSetLegend2(0.15,0.84-(0.85*textsizeLabelsPP*5),0.40,0.96, 0.85* textSizeLabelsPixel);
        legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaCombCombFitSysA,"Data","pf");
        legendRatioTheoryppEta_3Parted->AddEntry(histoRatioPythia8ToFitEta,  "PYTHIA 8.2, Monash 2013", "l");
        legendRatioTheoryppEta_3Parted->AddEntry(histoRatioPythia8_4CToFitEta,  "PYTHIA 8.2, Tune 4C", "l");
        legendRatioTheoryppEta_3Parted->AddEntry((TObject*)0, "NLO, PDF:CTEQ6M5 - FF:AESSS", "");
        legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaCombNLOMuHalf, "#mu = 0.5 #it{p}_{T}", "l");
        legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaCombNLOMuOne,  "#mu = #it{p}_{T}", "l");
        legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaCombNLOMuTwo,  "#mu = 2 #it{p}_{T}", "l");
        //legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaAESSS, "NLO, PDF:CTEQ6M5 - FF:AESSS", "");
        legendRatioTheoryppEta_3Parted->Draw();

        TLegend* legendRatioTheoryNormUncEta = GetAndSetLegend2(0.31,0.91,0.56,0.96, 0.85* textSizeLabelsPixel);
        legendRatioTheoryNormUncEta->AddEntry(boxErrorSigmaRatio,"norm. unc. 2.6%","l");
        legendRatioTheoryNormUncEta->Draw();

//        TLatex *labelRatioTheoryPP_Paper2   = new TLatex(0.2,0.68,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}");
//        SetStyleTLatex( labelRatioTheoryPP_Paper2, 0.8*textsizeLabelsPP,4);
//        labelRatioTheoryPP_Paper2->Draw();

        TLatex *labelRatioTheoryPP2   = new TLatex(0.78,0.925,collisionSystem8TeV.Data());
        SetStyleTLatex( labelRatioTheoryPP2, 0.85*textsizeLabelsPP,4);
        labelRatioTheoryPP2->Draw();
        TLatex *labelRatioTheoryPP221 = new TLatex(0.883,0.875,"ALICE");
        SetStyleTLatex( labelRatioTheoryPP221, 0.85*textsizeLabelsPP,4);
        labelRatioTheoryPP221->Draw();
        TLatex *labelRatioTheoryPP222= new TLatex(0.873,0.83,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioTheoryPP222, 0.85*textsizeLabelsPP,4);
        labelRatioTheoryPP222->Draw();

    canvasRatioPP->Update();
    canvasRatioPP->Print(Form("%s/Eta_RatioTheoryToData_PP.%s",outputDir.Data(),suffix.Data()));

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

        graphPi0DSS14->SetLineWidth(widthCommonFit);
        graphPi0DSS14->SetLineColor(colorNLO);
        graphPi0DSS14->SetLineStyle(1);
        graphPi0DSS14->SetFillStyle(1001);
        graphPi0DSS14->SetFillColor(colorNLO);
        graphPi0DSS14->Draw("same,e3");

        DrawGammaSetMarkerTGraphErr(graphPythia8_4CInvXSection, 0, 0, 807 , 807, widthLinesBoxes, kTRUE, 807);
        graphPythia8_4CInvXSection->Draw("3,same");
        DrawGammaSetMarker(histoPythia8_4CInvXSection, 24, 1.5, 807 , 807);
        histoPythia8_4CInvXSection->SetLineWidth(widthCommonFit);
        histoPythia8_4CInvXSection->Draw("same,hist,l");

        DrawGammaSetMarkerTGraphErr(graphPythia8InvXSection, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
        graphPythia8InvXSection->Draw("3,same");
        DrawGammaSetMarker(histoPythia8InvXSection, 24, 1.5, kRed+2 , kRed+2);
        histoPythia8InvXSection->SetLineWidth(widthCommonFit);
        histoPythia8InvXSection->Draw("same,hist,l");

        DrawGammaSetMarkerTF1( fitTCMInvXSectionPi0Plot, 7, 2, kGray+2);
        fitTCMInvXSectionPi0Plot->Draw("same");

        DrawGammaSetMarkerTF1( fitInvXSectionPi0, 3, 2, kGray+1);
        fitInvXSectionPi0->Draw("same");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombPi0InvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatA_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
        graphCombPi0InvXSectionStatA_WOXErr->Draw("p,same,z");

        TLatex *labelEnergyXSectionPaper= new TLatex(0.72, 0.91, collisionSystem8TeV.Data());
        SetStyleTLatex( labelEnergyXSectionPaper, textsizeLabelsXSecUp,4);
        labelEnergyXSectionPaper->Draw();
        TLatex *labelALICEXSectionPaper= new TLatex(0.848,0.86,"ALICE");
        SetStyleTLatex( labelALICEXSectionPaper, textsizeLabelsXSecUp,4);
        labelALICEXSectionPaper->Draw();
        TLatex *labelDetSysXSectionPaper= new TLatex(0.824,0.815,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPaper, textsizeLabelsXSecUp,4);
        labelDetSysXSectionPaper->Draw();

        TLegend* legendXsectionPaper    = GetAndSetLegend2(0.17, 0.08, 0.5, 0.18+0.05*6, textSizeLabelsPixel);
        legendXsectionPaper->SetNColumns(1);
        legendXsectionPaper->SetMargin(0.2);
        legendXsectionPaper->AddEntry(graphCombPi0InvXSectionSysA,"Data","pf");
        legendXsectionPaper->AddEntry(boxErrorSigmaRatio,"norm. unc. 2.6%","l");
        legendXsectionPaper->AddEntry(fitTCMInvXSectionPi0Plot,"TCM fit","l");
        legendXsectionPaper->AddEntry(fitInvXSectionPi0,"Tsallis fit","l");
//        legendXsectionPaper->AddEntry(graphNLODSS14Calc,"NLO, DSS14 ","f");
//        legendXsectionPaper->AddEntry(graphNLOCalcPi0MuHalf, "NLO, DSS07 #mu = 0.5 #it{p}_{T}", "l");
//        legendXsectionPaper->AddEntry(graphNLOCalcPi0MuOne,  "NLO, DSS07 #mu = #it{p}_{T}", "l");
//        legendXsectionPaper->AddEntry(graphNLOCalcPi0MuTwo,  "NLO, DSS07 #mu = 2 #it{p}_{T}", "l");
        legendXsectionPaper->AddEntry(histoPythia8InvXSection,"PYTHIA 8.2, Monash 2013","l");
        legendXsectionPaper->AddEntry(histoPythia8_4CInvXSection,"PYTHIA 8.2, Tune 4C","l");
        legendXsectionPaper->AddEntry(graphPi0DSS14,  "NLO, PDF:MSTW08 - FF:DSS14", "f");
        legendXsectionPaper->Draw();

        TLatex *labelRatioTheoryPP_Paper   = new TLatex(0.24,0.055,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}");
        SetStyleTLatex( labelRatioTheoryPP_Paper, 0.8*textsizeLabelsPP,4);
        labelRatioTheoryPP_Paper->Draw();

//        TLegend* legendXsectionPaper2     = GetAndSetLegend2(0.17, 0.06, 0.5, 0.12, textSizeLabelsPixel);
//        legendXsectionPaper2->SetMargin(0.2);
//        legendXsectionPaper2->AddEntry(fitTCMInvXSectionPi0Plot,"TCM fit","l");
//        //legendXsectionPaper2->AddEntry(fitTCMInvXSectionPi0Plot,"#it{A}_{e} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{n}","l");
//        legendXsectionPaper2->Draw();

    padInvSectionNLORatio->cd();
    padInvSectionNLORatio->SetLogx(1);
        TH2F * ratio2DNLO               = new TH2F("ratio2DNLO","ratio2DNLO",1000,0.23,50.,1000,0.6,1.95);
        SetStyleHistoTH2ForGraphs(ratio2DNLO, "#it{p}_{T} (GeV/#it{c})","#frac{NLO, Data}{TCM fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
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

//        TLegend* legendXsectionPaperPi02     = GetAndSetLegend2(0.17, 0.8, 0.4, 0.83+0.05*1, textSizeLabelsPixel*0.8);
//        legendXsectionPaperPi02->AddEntry(boxErrorSigmaRatio,"norm. unc. 2.6%","f");
//        legendXsectionPaperPi02->Draw();

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

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatA_WOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA->Draw("2,same");
        graphRatioPi0CombCombFitStatA_WOXErr->Draw("p,same");

        boxErrorSigmaRatio->Draw();
        DrawGammaLines(0.23, 50. , 1., 1.,0.5, kGray+2);

    padInvSectionPythiaRatio->cd();
    padInvSectionPythiaRatio->SetLogx(1);
        TH2F * ratio2DPythia            = new TH2F("ratio2DPythia","ratio2DPythia",1000,0.23,50.,1000,0.6,1.95);
        SetStyleHistoTH2ForGraphs(ratio2DPythia, "#it{p}_{T} (GeV/#it{c})","#frac{Pythia, Data}{TCM fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
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

        DrawGammaSetMarkerTGraphErr(graphRatioPythia8_4CToFit, 0, 0, 807 , 807, widthLinesBoxes, kTRUE, 807);
        graphRatioPythia8_4CToFit->Draw("3,same");
        DrawGammaSetMarker(histoRatioPythia8_4CToFit, 24, 1.5, 807 , 807);
        histoRatioPythia8_4CToFit->SetLineWidth(widthCommonFit);
        histoRatioPythia8_4CToFit->Draw("same,hist,l");

        DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFit, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
        graphRatioPythia8ToFit->Draw("3,same");
        DrawGammaSetMarker(histoRatioPythia8ToFit, 24, 1.5, kRed+2 , kRed+2);
        histoRatioPythia8ToFit->SetLineWidth(widthCommonFit);
        histoRatioPythia8ToFit->Draw("same,hist,l");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatA_WOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA->Draw("2,same");
        graphRatioPi0CombCombFitStatA_WOXErr->Draw("p,same");

        boxErrorSigmaRatio->Draw();
        DrawGammaLines(0.23, 50. , 1., 1.,0.5, kGray+2);

    canvasInvSectionPaper->Print(Form("%s/Pi0_InvXSectionWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));

    padInvSectionSpec->cd();
    padInvSectionSpec->SetLogy(1);
    padInvSectionSpec->SetLogx(1);
        SetStyleHistoTH2ForGraphs(histo2DXSectionEta, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",
                                0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.2/(textsizeFacXSecUp*marginXSec));
        histo2DXSectionEta->GetXaxis()->SetMoreLogLabels();
        histo2DXSectionEta->GetXaxis()->SetLabelOffset(+0.01);
        histo2DXSectionEta->Draw();

        DrawGammaNLOTGraph( graphNLOCalcEtaMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
        graphNLOCalcEtaMuHalf->Draw("same,c");
        DrawGammaNLOTGraph( graphNLOCalcEtaMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
        graphNLOCalcEtaMuOne->Draw("same,c");
        DrawGammaNLOTGraph( graphNLOCalcEtaMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
        graphNLOCalcEtaMuTwo->Draw("same,c");

        DrawGammaSetMarkerTGraphErr(graphPythia8_4CInvXSectionEta, 0, 0, 807 , 807, widthLinesBoxes, kTRUE, 807);
        graphPythia8_4CInvXSectionEta->Draw("3,same");
        DrawGammaSetMarker(histoPythia8_4CInvXSectionEta, 24, 1.5, 807 , 807);
        histoPythia8_4CInvXSectionEta->SetLineWidth(widthCommonFit);
        histoPythia8_4CInvXSectionEta->Draw("same,hist,l");

        DrawGammaSetMarkerTGraphErr(graphPythia8InvXSectionEta, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
        graphPythia8InvXSectionEta->Draw("3,same");
        DrawGammaSetMarker(histoPythia8InvXSectionEta, 24, 1.5, kRed+2 , kRed+2);
        histoPythia8InvXSectionEta->SetLineWidth(widthCommonFit);
        histoPythia8InvXSectionEta->Draw("same,hist,l");

        DrawGammaSetMarkerTF1( fitTCMInvXSectionEta, 7, 2, kGray+2);
        fitTCMInvXSectionEta->Draw("same");

        DrawGammaSetMarkerTF1( fitInvXSectionEta, 3, 2, kGray+1);
        fitInvXSectionEta->Draw("same");

        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombEtaInvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStatA_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
        graphCombEtaInvXSectionStatA_WOXErr->Draw("p,same,z");

        labelEnergyXSectionPaper->Draw();
        labelALICEXSectionPaper->Draw();
        TLatex *labelDetSysXSectionPaperEta = new TLatex(0.84,0.815,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPaperEta, textsizeLabelsXSecUp,4);
        labelDetSysXSectionPaperEta->Draw();

        TLegend* legendXsectionPaperEta     = GetAndSetLegend2(0.17, 0.03, 0.5, 0.13+0.05*8, textSizeLabelsPixel);
        legendXsectionPaperEta->SetNColumns(1);
        legendXsectionPaperEta->SetMargin(0.2);
        legendXsectionPaperEta->AddEntry(graphCombPi0InvXSectionSysA,"Data","pf");
        legendXsectionPaperEta->AddEntry(boxErrorSigmaRatio, "norm. unc. 2.6%", "l");
        legendXsectionPaperEta->AddEntry(fitTCMInvXSectionEta,"TCM fit","l");
        legendXsectionPaperEta->AddEntry(fitInvXSectionEta,"Tsallis fit","l");
        legendXsectionPaperEta->AddEntry(histoPythia8InvXSectionEta,"PYTHIA 8.2, Monash 2013","l");
        legendXsectionPaperEta->AddEntry(histoPythia8_4CInvXSectionEta,"PYTHIA 8.2, Tune 4C","l");
        legendXsectionPaperEta->AddEntry((TObject*)0, "", "");
        legendXsectionPaperEta->AddEntry(graphNLOCalcEtaMuHalf,"#mu = 0.5 #it{p}_{T}","l");
        legendXsectionPaperEta->AddEntry(graphNLOCalcEtaMuOne,"#mu = #it{p}_{T}","l");
        legendXsectionPaperEta->AddEntry(graphNLOCalcEtaMuTwo,"#mu = 2 #it{p}_{T}","l");
        //legendXsectionPaperEta->AddEntry(graphEtaAESSS,"NLO, PDF:CTEQ6M5 - FF:AESSS", "f");
        legendXsectionPaperEta->Draw();

        TLatex *labelEta = new TLatex(0.175, 0.192,"NLO, PDF:CTEQ6M5 - FF:AESSS");
        SetStyleTLatex( labelEta, 0.75*textsizeLabelsPP,4);
        labelEta->Draw();

        //labelRatioTheoryPP_Paper->Draw();


    padInvSectionNLORatio->cd();
    padInvSectionNLORatio->SetLogx(1);
        TH2F * ratio2DNLOEta                = new TH2F("ratio2DNLOEta","ratio2DNLOEta",1000,0.23,50.,1000,0.35,3.15);
        SetStyleHistoTH2ForGraphs(ratio2DNLOEta, "#it{p}_{T} (GeV/#it{c})","#frac{NLO, Data}{TCM fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                  0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
        ratio2DNLOEta->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio2DNLOEta->GetYaxis()->SetNdivisions(505);
        ratio2DNLOEta->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DNLOEta->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DNLOEta->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DNLOEta->GetXaxis()->SetLabelFont(42);
        ratio2DNLOEta->GetYaxis()->SetLabelFont(42);
        ratio2DNLOEta->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DNLOEta->GetXaxis()->SetTickLength(0.07);
        ratio2DNLOEta->DrawCopy();

//        TLegend* legendXsectionPaperEta2     = GetAndSetLegend2(0.17, 0.8, 0.4, 0.83+0.05*1, textSizeLabelsPixel*0.8);
//        legendXsectionPaperEta2->AddEntry(boxErrorSigmaRatio,"norm. unc. 2.6%","f");
//        legendXsectionPaperEta2->Draw();

        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
        graphRatioEtaCombNLOMuHalf->Draw("same,c");
        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
        graphRatioEtaCombNLOMuOne->Draw("same,c");
        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
        graphRatioEtaCombNLOMuTwo->Draw("same,c");

        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatA_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioEtaCombCombFitStatA_WOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioEtaCombCombFitSysA->SetLineWidth(0);
        graphRatioEtaCombCombFitSysA->Draw("2,same");
        graphRatioEtaCombCombFitStatA_WOXErr->Draw("p,same");

        boxErrorSigmaRatio->Draw();
        DrawGammaLines(0.23, 50. , 1., 1.,0.5, kGray+2);

    padInvSectionPythiaRatio->cd();
    padInvSectionPythiaRatio->SetLogx(1);
        TH2F * ratio2DPythiaEta             = new TH2F("ratio2DPythiaEta","ratio2DPythiaEta",1000,0.23,50.,1000,0.6,1.95);
        SetStyleHistoTH2ForGraphs(ratio2DPythiaEta, "#it{p}_{T} (GeV/#it{c})","#frac{Pythia, Data}{TCM fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                                  0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.2/(textsizeFacXSecDown*marginXSec), 510, 505);
        ratio2DPythiaEta->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio2DPythiaEta->GetYaxis()->SetNdivisions(505);
        ratio2DPythiaEta->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DPythiaEta->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DPythiaEta->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DPythiaEta->GetXaxis()->SetLabelFont(42);
        ratio2DPythiaEta->GetYaxis()->SetLabelFont(42);
        ratio2DPythiaEta->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DPythiaEta->GetXaxis()->SetTickLength(0.06);
        ratio2DPythiaEta->GetYaxis()->SetTickLength(0.04);
        ratio2DPythiaEta->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioPythia8_4CToFitEta, 0, 0, 807 , 807, widthLinesBoxes, kTRUE, 807);
        graphRatioPythia8_4CToFitEta->Draw("3,same");
        DrawGammaSetMarker(histoRatioPythia8_4CToFitEta, 24, 1.5, 807 , 807);
        histoRatioPythia8_4CToFitEta->SetLineWidth(widthCommonFit);
        histoRatioPythia8_4CToFitEta->Draw("same,hist,l");

        DrawGammaSetMarkerTGraphErr(graphRatioPythia8ToFitEta, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
        graphRatioPythia8ToFitEta->Draw("3,same");
        DrawGammaSetMarker(histoRatioPythia8ToFitEta, 24, 1.5, kRed+2 , kRed+2);
        histoRatioPythia8ToFitEta->SetLineWidth(widthCommonFit);
        histoRatioPythia8ToFitEta->Draw("same,hist,l");

        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatA_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioEtaCombCombFitStatA_WOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioEtaCombCombFitSysA->SetLineWidth(0);
        graphRatioEtaCombCombFitSysA->Draw("2,same");
        graphRatioEtaCombCombFitStatA_WOXErr->Draw("p,same");

        boxErrorSigmaRatio->Draw();
        DrawGammaLines(0.23, 50. , 1., 1.,0.5, kGray+2);

    canvasInvSectionPaper->Print(Form("%s/Eta_InvXSectionWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));

    // ***************************************************************************************************************
    // ******************************* pi0+eta combined plot  ********************************************************
    // ***************************************************************************************************************

    canvasXSectionPi0->cd();

    TH2F * histo2DXSectionWithEtaAndPi0;
    histo2DXSectionWithEtaAndPi0          = new TH2F("histo2DXSectionWithEtaAndPi0","histo2DXSectionWithEtaAndPi0",11000,0.23,50.,1000,6e-2,9e11);
    SetStyleHistoTH2ForGraphs(histo2DXSectionWithEtaAndPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.45);
    histo2DXSectionWithEtaAndPi0->GetXaxis()->SetMoreLogLabels();
    histo2DXSectionWithEtaAndPi0->GetXaxis()->SetNoExponent();
    histo2DXSectionWithEtaAndPi0->Draw("copy");

        // scale eta graphs
        Double_t scaleFacEtaForCombPlot                              = 1e-2;
        TGraphAsymmErrors* graphCombEtaInvXSectionStatA_WOXErrCopy   = (TGraphAsymmErrors*) graphCombEtaInvXSectionStatA_WOXErr->Clone("graphCombEtaInvXSectionStatAWOXErrCopy");
        TGraphAsymmErrors* graphCombEtaInvXSectionSysA_Copy          = (TGraphAsymmErrors*) graphCombEtaInvXSectionSysA->Clone("graphCombEtaInvXSectionSysA_Copy");
        graphCombEtaInvXSectionStatA_WOXErrCopy                      = ScaleGraph(graphCombEtaInvXSectionStatA_WOXErrCopy,scaleFacEtaForCombPlot);
        graphCombEtaInvXSectionSysA_Copy                             = ScaleGraph(graphCombEtaInvXSectionSysA_Copy,scaleFacEtaForCombPlot);

        TH1D* histFitTCMInvXSectionEta                               = (TH1D*)fitTCMInvXSectionEta->GetHistogram();
        histFitTCMInvXSectionEta->Scale(scaleFacEtaForCombPlot);

        TF1* tf1FitInvXSectionEta = new TF1("tf1FitInvXSectionEta","(1/100) * fitInvCrossSectionEta8TeV", 0.4, 35.);

        histoPythia8InvXSectionEta->Scale(scaleFacEtaForCombPlot);
        histoPythia8_4CInvXSectionEta->Scale(scaleFacEtaForCombPlot);
        TGraphErrors* graphPythia8InvXSectionEtaScaled              = ScaleGraph(graphPythia8InvXSectionEta,scaleFacEtaForCombPlot);
        TGraphErrors* graphPythia8T4CInvXSectionEtaScaled           = ScaleGraph(graphPythia8_4CInvXSectionEta,scaleFacEtaForCombPlot);

        TGraphAsymmErrors* graphEtaAESSSCopy                        = (TGraphAsymmErrors*)graphEtaAESSS->Clone("graphEtaAESSSCopy");
        graphEtaAESSSCopy                                           = ScaleGraph(graphEtaAESSSCopy,scaleFacEtaForCombPlot);

        // plotting NLO calcs pi0
        graphPi0DSS14->SetLineWidth(widthCommonFit);
        graphPi0DSS14->SetLineColor(colorNLO);
        graphPi0DSS14->SetLineStyle(1);
        graphPi0DSS14->SetFillStyle(1001);
        graphPi0DSS14->SetFillColor(colorNLO);
        graphPi0DSS14->Draw("same,e3");

        // plotting NLO calcs eta
        DrawGammaSetMarkerTGraphAsym(graphEtaAESSSCopy, 0, 0, colorCGC, colorCGC, widthLinesBoxes, kTRUE, colorCGC);
        graphEtaAESSSCopy->Draw("3,same");

        // plotting Pythia 8.2 Monash
        DrawGammaSetMarkerTGraphErr(graphPythia8_4CInvXSection, 0, 0, 807 , 807, widthLinesBoxes, kTRUE, 807);
        graphPythia8_4CInvXSection->Draw("3,same");
        DrawGammaSetMarker(histoPythia8_4CInvXSection, 24, 1.5, 807 , 807);
        histoPythia8_4CInvXSection->SetLineWidth(widthCommonFit);
        histoPythia8_4CInvXSection->Draw("same,hist,l");

        DrawGammaSetMarkerTGraphErr(graphPythia8InvXSection, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
        graphPythia8InvXSection->Draw("3,same");
        DrawGammaSetMarker(histoPythia8InvXSection, 24, 1.5, kRed+2 , kRed+2);
        histoPythia8InvXSection->SetLineWidth(widthCommonFit);
        histoPythia8InvXSection->Draw("same,hist,l");

        // plotting Pythia 8.2 Tune4C
        DrawGammaSetMarkerTGraphErr(graphPythia8T4CInvXSectionEtaScaled, 0, 0, 807 , 807, widthLinesBoxes, kTRUE, 807);
        graphPythia8T4CInvXSectionEtaScaled->Draw("3,same");
        DrawGammaSetMarker(histoPythia8_4CInvXSectionEta, 24, 1.5, 807 , 807);
        histoPythia8_4CInvXSectionEta->SetLineWidth(widthCommonFit);
        histoPythia8_4CInvXSectionEta->Draw("same,hist,l");

        DrawGammaSetMarkerTGraphErr(graphPythia8InvXSectionEtaScaled, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
        graphPythia8InvXSectionEtaScaled->Draw("3,same");
        DrawGammaSetMarker(histoPythia8InvXSectionEta, 24, 1.5, kRed+2 , kRed+2);
        histoPythia8InvXSectionEta->SetLineWidth(widthCommonFit);
        histoPythia8InvXSectionEta->Draw("same,hist,l");

        // plots fits
        fitTCMInvXSectionPi0Plot->Draw("same");
        fitInvXSectionPi0->Draw("same");

        SetStyleHisto(histFitTCMInvXSectionEta, 2, 7, kGray+2);
        histFitTCMInvXSectionEta->Draw("same,c");
        DrawGammaSetMarkerTF1( tf1FitInvXSectionEta, 3, 2, kGray+1);
        tf1FitInvXSectionEta->Draw("same");

        // plot data
        graphCombPi0InvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSysA_Copy, markerStyleComb+4, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
        graphCombEtaInvXSectionSysA_Copy->Draw("E2same");

        graphCombPi0InvXSectionStatA_WOXErr->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStatA_WOXErrCopy, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
        graphCombEtaInvXSectionStatA_WOXErrCopy->Draw("p,same,z");

        // labels lower left corner
        TLegend* legendXsectionPaperAll    = GetAndSetLegend2(0.17, 0.12, 0.5, 0.11+0.04*4, textSizeLabelsPixel, 1, "", 43, 0.2);
        legendXsectionPaperAll->AddEntry(graphCombPi0InvXSectionSysA,"#pi^{0}","pf");
        legendXsectionPaperAll->AddEntry(graphCombEtaInvXSectionSysA_Copy,"#eta (x 10^{-2})","pf");
        legendXsectionPaperAll->AddEntry(fitTCMInvXSectionPi0Plot,"TCM fit","l");
        legendXsectionPaperAll->AddEntry(fitInvXSectionPi0,"Tsallis fit","l");
        legendXsectionPaperAll->Draw();

        TLatex *labelEnergyXSectionPaperAll = new TLatex(0.18, 0.12+0.04*6, collisionSystem8TeV.Data());
        SetStyleTLatex( labelEnergyXSectionPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelEnergyXSectionPaperAll->Draw();
        TLatex *labelALICEXSectionPaperAll  = new TLatex(0.18,0.12+0.04*5,"ALICE");
        SetStyleTLatex( labelALICEXSectionPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelALICEXSectionPaperAll->Draw();
        TLatex *labelALICENormUnPaperAll    = new TLatex(0.18,0.12+0.04*4+0.003,"norm. unc. 2.6%");
        SetStyleTLatex( labelALICENormUnPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelALICENormUnPaperAll->Draw();

        // labels upper right corner
        TLegend* legendXsectionPaperPyBoth  = GetAndSetLegend2(0.5, 0.95-0.04*5, 0.54+0.33, 0.95, textSizeLabelsPixel, 1, "", 43, 0.18);
        legendXsectionPaperPyBoth->AddEntry(histoPythia8InvXSectionEta,"PYTHIA 8.2, Monash 2013","l");
        legendXsectionPaperPyBoth->AddEntry(histoPythia8_4CInvXSectionEta,"PYTHIA 8.2, Tune 4C","l");
        legendXsectionPaperPyBoth->AddEntry(graphPi0DSS14,"#pi^{0} pQCD NLO ","f");
        legendXsectionPaperPyBoth->AddEntry((TObject*)0,"#scale[0.75]{PDF: MSTW, FF: DSS14}","");
        legendXsectionPaperPyBoth->AddEntry(graphEtaAESSSCopy,"#eta pQCD NLO ","f");
        legendXsectionPaperPyBoth->AddEntry((TObject*)0,"#scale[0.75]{PDF: CTEQ6M5, FF: AESSS}","");
        legendXsectionPaperPyBoth->Draw();

    canvasXSectionPi0->SaveAs(Form("%s/InvXSection_Pi0_Eta_Theory.%s",outputDir.Data(),suffix.Data()));

    histo2DXSectionWithEtaAndPi0->Draw("copy");

        // plots fits
        fitTCMInvXSectionPi0Plot->Draw("same");
        fitInvXSectionPi0->Draw("same");

        SetStyleHisto(histFitTCMInvXSectionEta, 2, 7, kGray+2);
        histFitTCMInvXSectionEta->Draw("same,c");
        DrawGammaSetMarkerTF1( tf1FitInvXSectionEta, 3, 2, kGray+1);
        tf1FitInvXSectionEta->Draw("same");

        // plot data
        graphCombPi0InvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSysA_Copy, markerStyleComb+4, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
        graphCombEtaInvXSectionSysA_Copy->Draw("E2same");

        graphCombPi0InvXSectionStatA_WOXErr->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStatA_WOXErrCopy, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
        graphCombEtaInvXSectionStatA_WOXErrCopy->Draw("p,same,z");

        // labels lower left corner
        legendXsectionPaperAll->Draw();

        labelEnergyXSectionPaperAll->Draw();
        labelALICEXSectionPaperAll->Draw();
        labelALICENormUnPaperAll->Draw();

    canvasXSectionPi0->SaveAs(Form("%s/InvXSection_Pi0_Eta_Fits.%s",outputDir.Data(),suffix.Data()));

    histo2DXSectionWithEtaAndPi0->Draw("copy");

        // plot data
        graphCombPi0InvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSysA_Copy, markerStyleComb+4, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
        graphCombEtaInvXSectionSysA_Copy->Draw("E2same");

        graphCombPi0InvXSectionStatA_WOXErr->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStatA_WOXErrCopy, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
        graphCombEtaInvXSectionStatA_WOXErrCopy->Draw("p,same,z");

        // labels lower left corner
        TLegend* legendXsectionPaperAll2    = GetAndSetLegend2(0.17, 0.20, 0.5, 0.19+0.04*2, textSizeLabelsPixel, 1, "", 43, 0.2);
        legendXsectionPaperAll2->AddEntry(graphCombPi0InvXSectionSysA,"#pi^{0}","pf");
        legendXsectionPaperAll2->AddEntry(graphCombEtaInvXSectionSysA_Copy,"#eta (x 10^{-2})","pf");
        legendXsectionPaperAll2->Draw();

        labelEnergyXSectionPaperAll->Draw();
        labelALICEXSectionPaperAll->Draw();
        labelALICENormUnPaperAll->Draw();

    canvasXSectionPi0->SaveAs(Form("%s/InvXSection_Pi0_Eta.%s",outputDir.Data(),suffix.Data()));

    histoPythia8InvXSectionEta->Scale(1/scaleFacEtaForCombPlot);
    histoPythia8_4CInvXSectionEta->Scale(1/scaleFacEtaForCombPlot);

    // ***************************************************************************************************************
    // ******************************* eta/pi0 graphs without x-error  ***********************************************
    // ***************************************************************************************************************

    TGraphAsymmErrors* graphCombEtaToPi0StatA_WOXErr = (TGraphAsymmErrors*) graphCombEtaToPi0StatA->Clone("graphCombEtaToPi0StatA_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphCombEtaToPi0StatA_WOXErr);
    TGraphAsymmErrors* graphPCMEtaToPi0Stat_WOXErr = (TGraphAsymmErrors*) graphPCMEtaToPi0Stat->Clone("graphPCMEtaToPi0Stat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphPCMEtaToPi0Stat_WOXErr);
    TGraphAsymmErrors* graphEMCALEtaToPi0Stat_WOXErr = (TGraphAsymmErrors*) graphEMCALEtaToPi0Stat->Clone("graphEMCALEtaToPi0Stat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphEMCALEtaToPi0Stat_WOXErr);
    TGraphAsymmErrors* graphPCMEMCALEtaToPi0Stat_WOXErr = (TGraphAsymmErrors*) graphPCMEMCALEtaToPi0Stat->Clone("graphPCMEMCALEtaToPi0Stat_WOXErr");
    ProduceGraphAsymmWithoutXErrors(graphPCMEMCALEtaToPi0Stat_WOXErr);

    // ***************************************************************************************************************
    // ******************************* Plotting eta/pi0 ratio for single measurements ********************************
    // ***************************************************************************************************************
    textSizeLabelsPixel                 = 54;
    TCanvas* canvasEtatoPi0combo       = new TCanvas("canvasEtatoPi0combo","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasEtatoPi0combo, 0.09, 0.01, 0.01, 0.125);
    canvasEtatoPi0combo->SetLogx();

        Double_t textsizeLabelsEtaToPi0 = 0;
        Double_t textsizeFacEtaToPi0    = 0;
        if (canvasEtatoPi0combo->XtoPixel(canvasEtatoPi0combo->GetX2()) <canvasEtatoPi0combo->YtoPixel(canvasEtatoPi0combo->GetY1()) ){
            textsizeLabelsEtaToPi0      = (Double_t)textSizeLabelsPixel/canvasEtatoPi0combo->XtoPixel(canvasEtatoPi0combo->GetX2()) ;
            textsizeFacEtaToPi0         = (Double_t)1./canvasEtatoPi0combo->XtoPixel(canvasEtatoPi0combo->GetX2()) ;
        } else {
            textsizeLabelsEtaToPi0      = (Double_t)textSizeLabelsPixel/canvasEtatoPi0combo->YtoPixel(canvasEtatoPi0combo->GetY1());
            textsizeFacEtaToPi0         = (Double_t)1./canvasEtatoPi0combo->YtoPixel(canvasEtatoPi0combo->GetY1());
        }
        cout << textsizeLabelsEtaToPi0 << endl;

        textsizeLabelsEtaToPi0*=0.9;

    TH2F * histo2DEtatoPi0combo;
    histo2DEtatoPi0combo               = new TH2F("histo2DEtatoPi0combo","histo2DEtatoPi0combo",1000,0.33,32.,1000,-0.05,1.05    );
    SetStyleHistoTH2ForGraphs(histo2DEtatoPi0combo, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.85*textsizeLabelsEtaToPi0, textsizeLabelsEtaToPi0,
                              0.85*textsizeLabelsEtaToPi0,textsizeLabelsEtaToPi0, 0.9, 0.65, 510, 510);
    histo2DEtatoPi0combo->GetXaxis()->SetMoreLogLabels();
    histo2DEtatoPi0combo->GetXaxis()->SetNoExponent(kTRUE);
//    histo2DEtatoPi0combo->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(-10,10);
    histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(-0.05,1.05);
    histo2DEtatoPi0combo->Draw("copy");

        // plotting systematics graphs
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaToPi0Sys, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0], widthLinesBoxes, kTRUE);
        graphPCMEtaToPi0Sys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaToPi0Sys, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4], widthLinesBoxes, kTRUE);
        graphPCMEMCALEtaToPi0Sys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaToPi0Sys, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2], widthLinesBoxes, kTRUE);
        graphEMCALEtaToPi0Sys->Draw("E2same");

        // plotting statistics graphs
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaToPi0Stat_WOXErr, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
        graphPCMEtaToPi0Stat_WOXErr->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaToPi0Stat_WOXErr, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4]);
        graphPCMEMCALEtaToPi0Stat_WOXErr->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaToPi0Stat_WOXErr, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
        graphEMCALEtaToPi0Stat_WOXErr->Draw("p,same,e");

        TLegend* legendEtaToPi0 = GetAndSetLegend2(0.47, 0.15, 0.9, 0.15+(textsizeLabelsEtaToPi0*3*0.9), textSizeLabelsPixel);
        legendEtaToPi0->AddEntry(graphPCMEtaToPi0Sys,"PCM","pf");
        legendEtaToPi0->AddEntry(graphPCMEMCALEtaToPi0Sys,"PCM-EMC","pf");
        legendEtaToPi0->AddEntry(graphEMCALEtaToPi0Sys,"EMC","pf");
        legendEtaToPi0->Draw();

        TLatex *labelEnergyEtaToPi0 = new TLatex(0.13, 0.92,collisionSystem8TeV.Data());
        SetStyleTLatex( labelEnergyEtaToPi0, 0.85*textsizeLabelsEtaToPi0,4);
        labelEnergyEtaToPi0->Draw();

        TLatex *labelALICEEtaToPi0 = new TLatex(0.13, 0.92-(1*textsizeLabelsEtaToPi0*0.85),"ALICE");
        SetStyleTLatex( labelALICEEtaToPi0, 0.85*textsizeLabelsEtaToPi0,4);
        labelALICEEtaToPi0->Draw();


    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_differentSystems.%s",outputDir.Data(), suffix.Data()));

    //*************************************************************************************************************
    //*************************************************************************************************************

    histo2DEtatoPi0combo->Draw("copy");

    TLegend* legendXsectionPaperEtaToPi0     = GetAndSetLegend2(0.12, 0.8, 0.45, 0.96, 0.85*textSizeLabelsPixel);
    legendXsectionPaperEtaToPi0->SetNColumns(1);
    legendXsectionPaperEtaToPi0->SetMargin(0.2);
    legendXsectionPaperEtaToPi0->AddEntry(graphCombPi0InvXSectionSysA,"Data","pf");
    legendXsectionPaperEtaToPi0->AddEntry(graphEtaToPi07000GeV,"ALICE pp, #sqrt{#it{s}} = 7 TeV","p");
    legendXsectionPaperEtaToPi0->AddEntry(graphEtaToPi02760GeV,"ALICE pp, #sqrt{#it{s}} = 2.76 TeV","p");
    legendXsectionPaperEtaToPi0->Draw();

    DrawGammaSetMarkerTGraphAsym(graphEtaToPi07000GeV, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[1] , colorDet[1], widthLinesBoxes, kTRUE);
    graphEtaToPi07000GeV->Draw("same,p");
    DrawGammaSetMarkerTGraphAsym(graphEtaToPi02760GeV, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4], widthLinesBoxes, kTRUE);
    graphEtaToPi02760GeV->Draw("same,p");

    // plotting labels
    TLatex *labelEnergyEtaToPi02 = new TLatex(0.75, 0.92,collisionSystem8TeV.Data());
    SetStyleTLatex( labelEnergyEtaToPi02, 0.85*textsizeLabelsEtaToPi0,4);
    labelEnergyEtaToPi02->Draw();

    TLatex *labelALICEEtaToPi02 = new TLatex(0.852, 0.92-(1*textsizeLabelsEtaToPi0*0.85),"ALICE");
    SetStyleTLatex( labelALICEEtaToPi02, 0.85*textsizeLabelsEtaToPi0,4);
    labelALICEEtaToPi02->Draw();

    // plotting data
    graphCombEtaToPi0StatA_WOXErr->Print();
    DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0StatA_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
    graphCombEtaToPi0StatA_WOXErr->SetLineWidth(widthLinesBoxes);
    DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0SysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
    graphCombEtaToPi0SysA->SetLineWidth(0);
    graphCombEtaToPi0SysA->Draw("2,same");
    graphCombEtaToPi0StatA_WOXErr->Draw("p,same");

    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Paper.%s",outputDir.Data(), suffix.Data()));

    //*************************************************************************************************************
    //*************************************************************************************************************

    histo2DEtatoPi0combo->Draw("copy");

    // eta/pi0 mt-scaled
    TH1F *eta2pi0MtScaled = new TH1F("eta2pi0MtScaled","#eta/#pi^{0} from m_{T} scaling",5000,0.4,25.);
    eta2pi0MtScaled->SetLineColor(kBlue+2);
    eta2pi0MtScaled->SetLineWidth(2.);

    TH1F *eta2pi0MtScaledTCM = new TH1F("eta2pi0MtScaledTCM","#eta/#pi^{0} from m_{T} scaling",5000,0.4,25.);
    eta2pi0MtScaledTCM->SetLineColor(kBlue+2);
    eta2pi0MtScaledTCM->SetLineWidth(2.);

    Double_t eta2Pi0Const = 0.454963;
    Double_t mPi0 = 0.134977;
    Double_t mEta = 0.547853;
    for (Int_t i=1; i<=eta2pi0MtScaled->GetNbinsX(); i++) {
      Double_t ptPi0 = eta2pi0MtScaled->GetBinCenter(i);
      if (ptPi0 < 0.3) continue;
      Double_t mtEta = TMath::Sqrt(mEta*mEta + ptPi0*ptPi0);
      Double_t ptEta = TMath::Sqrt(mtEta*mtEta - mPi0*mPi0);
      Double_t Reta2pi0 = fitInvXSectionPi0->Eval(ptEta) / fitInvXSectionPi0->Eval(ptPi0) * eta2Pi0Const;
      eta2pi0MtScaled->SetBinContent(i,Reta2pi0);

      Double_t Reta2pi0TCM = fitTCMInvXSectionPi0Plot->Eval(ptEta) / fitTCMInvXSectionPi0Plot->Eval(ptPi0) * eta2Pi0Const;
      eta2pi0MtScaledTCM->SetBinContent(i,Reta2pi0TCM);
    }

    TGraphAsymmErrors* graphRatioForMt_stat     = (TGraphAsymmErrors*)graphCombEtaToPi0StatA_WOXErr->Clone();
    TGraphAsymmErrors* graphRatioForMt_sys      = (TGraphAsymmErrors*)graphCombEtaToPi0SysA->Clone();

    Int_t n_stat                = graphRatioForMt_stat->GetN();
    Double_t* yValue_stat       = graphRatioForMt_stat->GetY();
    Double_t* yErrorLow_stat    = graphRatioForMt_stat->GetEYlow();
    Double_t* yErrorHigh_stat   = graphRatioForMt_stat->GetEYhigh();
    Double_t* yValue_sys        = graphRatioForMt_sys->GetY();
    Double_t* yErrorLow_sys     = graphRatioForMt_sys->GetEYlow();
    Double_t* yErrorHigh_sys    = graphRatioForMt_sys->GetEYhigh();
    for (Int_t i = 0; i < n_stat; i++){
        Double_t ptPi0 = graphCombEtaToPi0StatA_WOXErr->GetX()[i];
        Double_t mtEta = TMath::Sqrt(mEta*mEta + ptPi0*ptPi0);
        Double_t ptEta = TMath::Sqrt(mtEta*mtEta - mPi0*mPi0);
        Double_t Reta2pi0 = fitInvXSectionPi0->Eval(ptEta) / fitInvXSectionPi0->Eval(ptPi0) * eta2Pi0Const;

        yValue_stat[i]       = yValue_stat[i] / Reta2pi0;
        yErrorLow_stat[i]    = yErrorLow_stat[i] / Reta2pi0;
        yErrorHigh_stat[i]   = yErrorHigh_stat[i] / Reta2pi0;

        yValue_sys[i]        = yValue_sys[i] / Reta2pi0;
        yErrorLow_sys[i]     = yErrorLow_sys[i] / Reta2pi0;
        yErrorHigh_sys[i]    = yErrorHigh_sys[i]/ Reta2pi0;
    }

    textSizeLabelsPixel = 48;
    TLegend* legendXsectionPaperEtaToPi03     = GetAndSetLegend2(0.11, 0.86, 0.96, 0.96, 0.85*textSizeLabelsPixel);
    legendXsectionPaperEtaToPi03->SetNColumns(2);
    legendXsectionPaperEtaToPi03->SetMargin(0.15);
    legendXsectionPaperEtaToPi03->AddEntry(graphCombPi0InvXSectionSysA,"ALICE pp, #sqrt{#it{s}} = 8 TeV","pf");
    legendXsectionPaperEtaToPi03->AddEntry(graphEtaToPi07000GeV,"ALICE pp, #sqrt{#it{s}} = 7 TeV","p");
    legendXsectionPaperEtaToPi03->AddEntry(graphEtaToPi02760GeV,"ALICE pp, #sqrt{#it{s}} = 2.76 TeV","p");
    legendXsectionPaperEtaToPi03->AddEntry(eta2pi0MtScaled,"ALICE pp m_{T}-scaled, #sqrt{#it{s}} = 8 TeV","l");
    legendXsectionPaperEtaToPi03->Draw();

    TLegend* legendXsectionPaperEtaToPi04     = GetAndSetLegend2(0.55, 0.15, 0.93, 0.25, 0.85*textSizeLabelsPixel);
    legendXsectionPaperEtaToPi04->SetNColumns(1);
    legendXsectionPaperEtaToPi04->SetMargin(0.15);
    legendXsectionPaperEtaToPi04->AddEntry(eta2pi0_RHIC200GeV,"PHENIX pp, #sqrt{#it{s}} = 200 GeV","p");
    legendXsectionPaperEtaToPi04->AddEntry(eta2pi0_NA27_275GeV,"NA27 pp, #sqrt{#it{s}} = 27.5 GeV","p");
    legendXsectionPaperEtaToPi04->Draw();

    DrawGammaSetMarkerTGraphErr(eta2pi0_RHIC200GeV, 27, 3., kGray+1 , kGray+1, 2., kFALSE);
    eta2pi0_RHIC200GeV->Draw("same,p");
    DrawGammaSetMarkerTGraphErr(eta2pi0_NA27_275GeV, 25, 2., kGray+2, kGray+2, 2., kFALSE);
    eta2pi0_NA27_275GeV->Draw("p,same");
    DrawGammaSetMarkerTGraphAsym(graphEtaToPi07000GeV, 28, 3., kBlue-6 , kBlue-6, 2., kFALSE);
    graphEtaToPi07000GeV->Draw("same,p");
    DrawGammaSetMarkerTGraphAsym(graphEtaToPi02760GeV, 30, 3., kBlue-6, kBlue-6, 2., kFALSE);
    graphEtaToPi02760GeV->Draw("same,p");
//    DrawGammaSetMarkerTGraphErr(eta2pi0_NA27_275GeV, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[3] , colorDet[3], widthLinesBoxes, kTRUE);
//    eta2pi0_NA27_275GeV->Draw("same,p");
//    DrawGammaSetMarkerTGraph(graphEtaToPi02760GeV, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4], widthLinesBoxes, kTRUE);
//    graphEtaToPi02760GeV->Draw("same,p");

    eta2pi0MtScaled->Draw("][ c same");

    // plotting data
    DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0StatA_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
    graphCombEtaToPi0StatA_WOXErr->SetLineWidth(widthLinesBoxes);
    DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0SysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
    graphCombEtaToPi0SysA->SetLineWidth(0);
    graphCombEtaToPi0SysA->Draw("2,same");
    graphCombEtaToPi0StatA_WOXErr->Draw("p,same");

    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_mT.%s",outputDir.Data(), suffix.Data()));

    //*************************************************************************************************************
    //*************************************************************************************************************

    histo2DEtatoPi0combo->Draw("copy");

    TLegend* legendXsectionPaperEtaToPi03TCM     = GetAndSetLegend2(0.11, 0.86, 0.96, 0.96, 0.85*textSizeLabelsPixel);
    legendXsectionPaperEtaToPi03TCM->SetNColumns(2);
    legendXsectionPaperEtaToPi03TCM->SetMargin(0.15);
    legendXsectionPaperEtaToPi03TCM->AddEntry(graphCombPi0InvXSectionSysA,"ALICE pp, #sqrt{#it{s}} = 8 TeV","pf");
    legendXsectionPaperEtaToPi03TCM->AddEntry(graphEtaToPi07000GeV,"ALICE pp, #sqrt{#it{s}} = 7 TeV","p");
    legendXsectionPaperEtaToPi03TCM->AddEntry(graphEtaToPi02760GeV,"ALICE pp, #sqrt{#it{s}} = 2.76 TeV","p");
    legendXsectionPaperEtaToPi03TCM->AddEntry(eta2pi0MtScaledTCM,"ALICE pp m_{T}-scaled, #sqrt{#it{s}} = 8 TeV","l");
    legendXsectionPaperEtaToPi03TCM->Draw();

    legendXsectionPaperEtaToPi04->Draw();

    eta2pi0_RHIC200GeV->Draw("same,p");
    eta2pi0_NA27_275GeV->Draw("p,same");
    graphEtaToPi07000GeV->Draw("same,p");
    graphEtaToPi02760GeV->Draw("same,p");

    eta2pi0MtScaledTCM->Draw("][ c same");

    // plotting data
    graphCombEtaToPi0SysA->Draw("2,same");
    graphCombEtaToPi0StatA_WOXErr->Draw("p,same");

    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_mT_TCM.%s",outputDir.Data(), suffix.Data()));
    //*************************************************************************************************************
    //*************************************************************************************************************

    TH2F * histo2DEtatoPi0ratio;
    histo2DEtatoPi0ratio               = new TH2F("histo2DEtatoPi0ratio","histo2DEtatoPi0ratio",1000,0.33,32.,1000,-0.05,1.55    );
    SetStyleHistoTH2ForGraphs(histo2DEtatoPi0ratio, "#it{p}_{T} (GeV/#it{c})","ratio", 0.85*textsizeLabelsEtaToPi0, textsizeLabelsEtaToPi0,
                              0.85*textsizeLabelsEtaToPi0,textsizeLabelsEtaToPi0, 0.9, 0.65, 510, 510);
    histo2DEtatoPi0ratio->GetXaxis()->SetMoreLogLabels();
    histo2DEtatoPi0ratio->GetXaxis()->SetNoExponent(kTRUE);
//    histo2DEtatoPi0ratio->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DEtatoPi0ratio->GetYaxis()->SetRangeUser(-10,10);
//    histo2DEtatoPi0ratio->GetYaxis()->SetRangeUser(-0.05,1.05);
    histo2DEtatoPi0ratio->Draw("copy");

    // plotting data
    DrawGammaSetMarkerTGraphAsym(graphRatioForMt_stat, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
    graphRatioForMt_stat->SetLineWidth(widthLinesBoxes);
    DrawGammaSetMarkerTGraphAsym(graphRatioForMt_sys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
    graphRatioForMt_sys->SetLineWidth(0);
    graphRatioForMt_sys->Draw("2,same");
    graphRatioForMt_stat->Draw("p,same");

    TLegend* legendXsectionPaperEtaToPiRatio     = GetAndSetLegend2(0.11, 0.86, 0.96, 0.96, 0.85*textSizeLabelsPixel);
    legendXsectionPaperEtaToPiRatio->SetNColumns(2);
    legendXsectionPaperEtaToPiRatio->SetMargin(0.15);
    legendXsectionPaperEtaToPiRatio->AddEntry(graphRatioForMt_sys,"ALICE pp, #sqrt{#it{s}} = 8 TeV - (#eta/#pi^{0})_{data}/(#eta/#pi^{0})_{m_{T}}","pf");
    legendXsectionPaperEtaToPiRatio->Draw();

    DrawGammaLines(0.33, 32. , 1., 1.,0.5, kGray+2);

    histo2DEtatoPi0ratio->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_mT_ratio.%s",outputDir.Data(), suffix.Data()));

    //*************************************************************************************************************
    //*************************************************************************************************************

    histo2DEtatoPi0combo->Draw("copy");

    textSizeLabelsPixel = 48;
    TLegend* legendXsectionPaperEtaToPi03n     = GetAndSetLegend2(0.11, 0.81, 0.5, 0.96, 0.85*textSizeLabelsPixel);
    legendXsectionPaperEtaToPi03n->SetNColumns(1);
    legendXsectionPaperEtaToPi03n->SetMargin(0.15);
    legendXsectionPaperEtaToPi03n->AddEntry(graphCombPi0InvXSectionSysA,"ALICE pp, #sqrt{#it{s}} = 8 TeV","pf");
    legendXsectionPaperEtaToPi03n->AddEntry(graphEtaToPi07000GeV,"ALICE pp, #sqrt{#it{s}} = 7 TeV","p");
    legendXsectionPaperEtaToPi03n->AddEntry(graphEtaToPi02760GeV,"ALICE pp, #sqrt{#it{s}} = 2.76 TeV","p");
    legendXsectionPaperEtaToPi03n->Draw();

    TLegend* legendXsectionPaperEtaToPi04n     = GetAndSetLegend2(0.53, 0.15, 0.93, 0.25, 0.85*textSizeLabelsPixel);
    legendXsectionPaperEtaToPi04n->SetNColumns(1);
    legendXsectionPaperEtaToPi04n->SetMargin(0.15);
    legendXsectionPaperEtaToPi04n->AddEntry(eta2pi0_RHIC200GeV,"PHENIX pp, #sqrt{#it{s}} = 200 GeV","p");
    legendXsectionPaperEtaToPi04n->AddEntry(eta2pi0_NA27_275GeV,"NA27 pp, #sqrt{#it{s}} = 27.5 GeV","p");
    legendXsectionPaperEtaToPi04n->Draw();

    eta2pi0_RHIC200GeV->Draw("same,p");
    eta2pi0_NA27_275GeV->Draw("p,same");
    graphEtaToPi07000GeV->Draw("same,p");
    graphEtaToPi02760GeV->Draw("same,p");

    // plotting data
    graphCombEtaToPi0SysA->Draw("2,same");
    graphCombEtaToPi0StatA_WOXErr->Draw("p,same");

    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Comparison_no_mT.%s",outputDir.Data(), suffix.Data()));

    //*************************************************************************************************************
    //*************************************************************************************************************

    histo2DEtatoPi0combo->Draw("copy");

    DrawGammaSetMarkerTGraphErr(graphPythia8T4CEtaToPi0, 0, 0, 807 , 807, widthLinesBoxes, kTRUE, 807);
    graphPythia8T4CEtaToPi0->Draw("3,same");
    DrawGammaSetMarker(histoPythia8T4CEtaToPi0, 24, 1.5, 807 , 807);
    histoPythia8T4CEtaToPi0->SetLineWidth(widthCommonFit);
    histoPythia8T4CEtaToPi0->Draw("same,hist,l");

    DrawGammaSetMarkerTGraphErr(graphPythia8EtaToPi0, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
    graphPythia8EtaToPi0->Draw("3,same");
    DrawGammaSetMarker(histoPythia8EtaToPi0, 24, 1.5, kRed+2 , kRed+2);
    histoPythia8EtaToPi0->SetLineWidth(widthCommonFit);
    histoPythia8EtaToPi0->Draw("same,hist,l");

    textSizeLabelsPixel = 48;
    TLegend* legendXsectionPaperEtaToPi02     = GetAndSetLegend2(0.12, 0.69, 0.45, 0.69+0.045*6, 0.85*textSizeLabelsPixel);
    legendXsectionPaperEtaToPi02->SetNColumns(1);
    legendXsectionPaperEtaToPi02->SetMargin(0.2);
    legendXsectionPaperEtaToPi02->AddEntry(graphCombPi0InvXSectionSysA,"Data","pf");
    legendXsectionPaperEtaToPi02->AddEntry(histoPythia8EtaToPi0,"PYTHIA 8.2, Monash 2013","l");
    legendXsectionPaperEtaToPi02->AddEntry(histoPythia8T4CEtaToPi0,"PYTHIA 8.2, Tune 4C","l");
    legendXsectionPaperEtaToPi02->AddEntry(graphNLOEtaToPi0,"NLO, PDF:CTEQ6M5","f");
    legendXsectionPaperEtaToPi02->AddEntry((TObject*)0,"#pi^{0} FF:DSS07, #eta FF:AESSS","");
    legendXsectionPaperEtaToPi02->AddEntry((TObject*)0,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}","");
    //legendXsectionPaperEtaToPi02->AddEntry(graphNLOEtaToPi0MuHalf,"#mu = 0.5 #it{p}_{T}","l");
    //legendXsectionPaperEtaToPi02->AddEntry(graphNLOEtaToPi0MuOne,"#mu = #it{p}_{T}","l");
    //legendXsectionPaperEtaToPi02->AddEntry(graphNLOEtaToPi0MuTwo,"#mu = 2 #it{p}_{T}","l");
    legendXsectionPaperEtaToPi02->Draw();

    // plotting NLO
    graphNLOEtaToPi0->SetLineWidth(widthCommonFit);
    graphNLOEtaToPi0->SetLineColor(colorNLO);
    graphNLOEtaToPi0->SetLineStyle(1);
    graphNLOEtaToPi0->SetFillStyle(1001);
    graphNLOEtaToPi0->SetFillColor(colorNLO);
    graphNLOEtaToPi0->Draw("same,e4");

    // plotting data
    DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0StatA_WOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
    graphCombEtaToPi0StatA_WOXErr->SetLineWidth(widthLinesBoxes);
    DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0SysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
    graphCombEtaToPi0SysA->SetLineWidth(0);
    graphCombEtaToPi0SysA->Draw("2,same");
    graphCombEtaToPi0StatA_WOXErr->Draw("p,same");

    // plotting labels
    labelEnergyEtaToPi02->Draw();
    labelALICEEtaToPi02->Draw();
    //labelPi0EtaToPi02->Draw();

    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Theory_Paper.%s",outputDir.Data(), suffix.Data()));

    //*************************************************************************************************************
    //*************************************************************************************************************

    histo2DEtatoPi0combo->Draw("copy");

    DrawGammaSetMarkerTGraphErr(graphPythia8T4CEtaToPi0, 0, 0, 807 , 807, widthLinesBoxes, kTRUE, 807);
    graphPythia8T4CEtaToPi0->Draw("3,same");
    DrawGammaSetMarker(histoPythia8T4CEtaToPi0, 24, 1.5, 807 , 807);
    histoPythia8T4CEtaToPi0->SetLineWidth(widthCommonFit);
    histoPythia8T4CEtaToPi0->Draw("same,hist,l");

    DrawGammaSetMarkerTGraphErr(graphPythia8EtaToPi0, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
    graphPythia8EtaToPi0->Draw("3,same");
    DrawGammaSetMarker(histoPythia8EtaToPi0, 24, 1.5, kRed+2 , kRed+2);
    histoPythia8EtaToPi0->SetLineWidth(widthCommonFit);
    histoPythia8EtaToPi0->Draw("same,hist,l");

    TLegend* legendXsectionPaperEtaToPi05     = GetAndSetLegend2(0.12, 0.645, 0.45, 0.645+0.045*7, 0.85*textSizeLabelsPixel);
    legendXsectionPaperEtaToPi05->SetNColumns(1);
    legendXsectionPaperEtaToPi05->SetMargin(0.2);
    legendXsectionPaperEtaToPi05->AddEntry(graphCombPi0InvXSectionSysA,"Data","pf");
    legendXsectionPaperEtaToPi05->AddEntry(eta2pi0MtScaled,"ALICE pp, #sqrt{#it{s}} = 8 TeV from m_{T} scaling","l");
    legendXsectionPaperEtaToPi05->AddEntry(histoPythia8EtaToPi0,"PYTHIA 8.2, Monash 2013","l");
    legendXsectionPaperEtaToPi05->AddEntry(histoPythia8T4CEtaToPi0,"PYTHIA 8.2, Tune 4C","l");
    legendXsectionPaperEtaToPi05->AddEntry(graphNLOEtaToPi0,"NLO, PDF:CTEQ6M5","f");
    legendXsectionPaperEtaToPi05->AddEntry((TObject*)0,"#pi^{0} FF:DSS07, #eta FF:AESSS","");
    legendXsectionPaperEtaToPi05->AddEntry((TObject*)0,"0.5#it{p}_{T} < #mu < 2#it{p}_{T}","");
    legendXsectionPaperEtaToPi05->Draw();

    // plotting NLO
    graphNLOEtaToPi0->Draw("same,e4");

    // plotting data
    graphCombEtaToPi0SysA->Draw("2,same");
    graphCombEtaToPi0StatA_WOXErr->Draw("p,same");

    // plotting labels
    labelEnergyEtaToPi02->Draw();
    labelALICEEtaToPi02->Draw();

    eta2pi0MtScaled->Draw("][ c same");

    // plotting data
    graphCombEtaToPi0SysA->Draw("2,same");
    graphCombEtaToPi0StatA_WOXErr->Draw("p,same");

    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Theory_Paper_mT.%s",outputDir.Data(), suffix.Data()));

    //*************************************************************************************************************
    //*************************************************************************************************************

    histo2DEtatoPi0combo->Draw("copy");

    // plotting data
    graphCombEtaToPi0SysA->Draw("2,same");
    graphCombEtaToPi0StatA_WOXErr->Draw("p,same");

    // plotting labels
    labelEnergyEtaToPi02->Draw();
    labelALICEEtaToPi02->Draw();
    //labelPi0EtaToPi02->Draw();

    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Combined.%s",outputDir.Data(), suffix.Data()));

    //*************************************************************************************************************
    //*************************************************************************************************************

    canvasEtatoPi0combo->SetRightMargin(0.02);
    histo2DEtatoPi0combo->GetXaxis()->SetLabelOffset(0.);
    histo2DEtatoPi0combo->GetXaxis()->SetRangeUser(0.,25.0);
    histo2DEtatoPi0combo->Draw("copy");
    legendXsectionPaperEtaToPi02->Draw();

    //plotting MC
    DrawGammaSetMarkerTGraphErr(graphPythia8T4CEtaToPi0, 0, 0, 807 , 807, widthLinesBoxes, kTRUE, 807);
    graphPythia8T4CEtaToPi0->Draw("3,same");
    DrawGammaSetMarker(histoPythia8T4CEtaToPi0, 24, 1.5, 807 , 807);
    histoPythia8T4CEtaToPi0->SetLineWidth(widthCommonFit);
    histoPythia8T4CEtaToPi0->Draw("same,hist,l");

    DrawGammaSetMarkerTGraphErr(graphPythia8EtaToPi0, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
    graphPythia8EtaToPi0->Draw("3,same");
    DrawGammaSetMarker(histoPythia8EtaToPi0, 24, 1.5, kRed+2 , kRed+2);
    histoPythia8EtaToPi0->SetLineWidth(widthCommonFit);
    histoPythia8EtaToPi0->Draw("same,hist,l");

    // plotting NLO
    graphNLOEtaToPi0->Draw("same,e4");

    // plotting data
    graphCombEtaToPi0SysA->Draw("2,same");
    graphCombEtaToPi0StatA_WOXErr->Draw("p,same");

    // plotting labels
    labelEnergyEtaToPi02->Draw();
    labelALICEEtaToPi02->Draw();
    //labelPi0EtaToPi02->Draw();

    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->SetLogx(kFALSE);
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Theory_Paper_LIN.%s",outputDir.Data(), suffix.Data()));

    //*************************************************************************************************************
    //*************************************************************************************************************

    histo2DEtatoPi0combo->Draw("copy");

    // plotting data
    graphCombEtaToPi0SysA->Draw("2,same");
    graphCombEtaToPi0StatA_WOXErr->Draw("p,same");

    // plotting labels
    labelEnergyEtaToPi02->Draw();
    labelALICEEtaToPi02->Draw();
    //labelPi0EtaToPi02->Draw();

    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->SetLogx(kFALSE);
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Combined_LIN.%s",outputDir.Data(), suffix.Data()));

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
    histo2DDummyATLAS               = new TH2F("histo2DDummyATLAS","histo2DDummyATLAS",1000,0.3,70.,1000,1,5e11);
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
    histo2DPi0RatioToCombFit->GetXaxis()->SetRangeUser(0.4,70);
    histo2DPi0RatioToCombFit->GetXaxis()->SetNoExponent(kTRUE);
    histo2DPi0RatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioATLASChargedFit, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioATLASChargedFit, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);

        graphRatioATLASChargedFit->Draw("E2same");
        graphRatioATLASChargedFit->Draw("p,same,z");

        DrawGammaLines(0.5, 70. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.5, 70. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.5, 70. , 0.9, 0.9,0.5, kGray, 7);

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
    histo2DCompCombinedRatio2       = new TH2F("histo2DCompCombinedRatio2","histo2DCompCombinedRatio2",1000,0.23,70.,1000,0.1,4.    );
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

    DrawGammaLines(0.23, 70 , 1, 1 ,1, kGray, 1);

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

    DrawGammaLines(0.23, 70 , 1, 1 ,1, kGray, 1);

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

    DrawGammaLines(0.23, 70 , 1, 1 ,1, kGray, 1);

    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedHadronToNeutralPionsWithMerged_PP8TeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    //*************************************************************************************************************
    //*************************************************************************************************************

    canvasEtatoPi0combo->cd();
    canvasEtatoPi0combo->SetLogx();

    TH2F * histoRatioEnergies;
    histoRatioEnergies               = new TH2F("histoRatioEnergies","histoRatioEnergies",1000,0.3,35.,1000,0.01,6.2    );
    SetStyleHistoTH2ForGraphs(histoRatioEnergies, "#it{p}_{T} (GeV/#it{c})","ratio", 0.85*textsizeLabelsEtaToPi0, textsizeLabelsEtaToPi0,
                              0.85*textsizeLabelsEtaToPi0,textsizeLabelsEtaToPi0, 0.9, 0.65, 510, 510);
    histoRatioEnergies->GetXaxis()->SetMoreLogLabels();
    histoRatioEnergies->GetXaxis()->SetNoExponent(kTRUE);
    histoRatioEnergies->DrawCopy();

    // plotting data 8TeV/2.76TeV
    TH1F *ratioOfEnergies = (TH1F*)histoPythia8InvXSection->Clone("ratioOfEnergies");
    ratioOfEnergies->SetTitle("");
    ratioOfEnergies->SetLineColor(kRed+2);
    ratioOfEnergies->SetMarkerColor(kRed+2);
    ratioOfEnergies->SetMarkerStyle(20);
    ratioOfEnergies->SetLineWidth(2.);

    TH1F *ratioOfEnergiesToPythia = (TH1F*)histoPythia8InvXSection->Clone("ratioOfEnergiesToPythia");
    ratioOfEnergiesToPythia->SetTitle("");
    ratioOfEnergiesToPythia->SetLineColor(kRed+2);
    ratioOfEnergiesToPythia->SetMarkerColor(kRed+2);
    ratioOfEnergiesToPythia->SetMarkerStyle(24);
    ratioOfEnergiesToPythia->SetLineWidth(2.);

    for (Int_t i=1; i<=ratioOfEnergies->GetNbinsX(); i++) {
      Double_t temp = fitTCMInvXSectionPi0->Eval(ratioOfEnergies->GetBinCenter(i))/fit2760GeVPi0TCM->Eval(ratioOfEnergies->GetBinCenter(i));
      ratioOfEnergies->SetBinContent(i,temp);
      ratioOfEnergies->SetBinError(i,0.);
    }

    for (Int_t i=1; i<=ratioOfEnergiesToPythia->GetNbinsX(); i++) {
      Double_t temp = histoPythia8InvXSection->GetBinContent(i)/histoPythia8InvXSection2760GeV->GetBinContent(i);
      ratioOfEnergiesToPythia->SetBinContent(i,temp);
      ratioOfEnergiesToPythia->SetBinError(i,0.);
    }

    ratioOfEnergies->DrawCopy("p,same");
    ratioOfEnergiesToPythia->DrawCopy("p,same");

    // plotting data 8TeV/7TeV
    TH1F *ratioOfEnergies2 = (TH1F*)histoPythia8InvXSection->Clone("ratioOfEnergies2");
    ratioOfEnergies2->SetTitle("");
    ratioOfEnergies2->SetLineColor(kBlue+2);
    ratioOfEnergies2->SetMarkerColor(kBlue+2);
    ratioOfEnergies2->SetMarkerStyle(20);
    ratioOfEnergies2->SetLineWidth(2.);

    TH1F *ratioOfEnergiesToPythia2 = (TH1F*)histoPythia8InvXSection->Clone("ratioOfEnergiesToPythia2");
    ratioOfEnergiesToPythia2->SetTitle("");
    ratioOfEnergiesToPythia2->SetLineColor(kBlue+2);
    ratioOfEnergiesToPythia2->SetMarkerColor(kBlue+2);
    ratioOfEnergiesToPythia2->SetMarkerStyle(24);
    ratioOfEnergiesToPythia2->SetLineWidth(2.);

    for (Int_t i=1; i<=ratioOfEnergies2->GetNbinsX(); i++) {
      Double_t temp = fitTCMInvXSectionPi0->Eval(ratioOfEnergies->GetBinCenter(i))/fit7TeVPi0TCM->Eval(ratioOfEnergies->GetBinCenter(i));
      ratioOfEnergies2->SetBinContent(i,temp);
      ratioOfEnergies2->SetBinError(i,0.);
    }

    for (Int_t i=1; i<=ratioOfEnergiesToPythia2->GetNbinsX(); i++) {
      Double_t temp = histoPythia8InvXSection->GetBinContent(i)/histoPythia8InvXSection7TeV->GetBinContent(i);
      ratioOfEnergiesToPythia2->SetBinContent(i,temp);
      ratioOfEnergiesToPythia2->SetBinError(i,0.);
    }

    ratioOfEnergies2->DrawCopy("p,same");
    ratioOfEnergiesToPythia2->DrawCopy("p,same");


    // plotting data 7TeV/2.76TeV
    TH1F *ratioOfEnergies3 = (TH1F*)histoPythia8InvXSection->Clone("ratioOfEnergies3");
    ratioOfEnergies3->SetTitle("");
    ratioOfEnergies3->SetLineColor(kGreen+2);
    ratioOfEnergies3->SetMarkerColor(kGreen+2);
    ratioOfEnergies3->SetMarkerStyle(20);
    ratioOfEnergies3->SetLineWidth(2.);

    TH1F *ratioOfEnergiesToPythia3 = (TH1F*)histoPythia8InvXSection->Clone("ratioOfEnergiesToPythia3");
    ratioOfEnergiesToPythia3->SetTitle("");
    ratioOfEnergiesToPythia3->SetLineColor(kGreen+2);
    ratioOfEnergiesToPythia3->SetMarkerColor(kGreen+2);
    ratioOfEnergiesToPythia3->SetMarkerStyle(24);
    ratioOfEnergiesToPythia3->SetLineWidth(2.);

    for (Int_t i=1; i<=ratioOfEnergies3->GetNbinsX(); i++) {
      Double_t temp = fit7TeVPi0TCM->Eval(ratioOfEnergies->GetBinCenter(i))/fit2760GeVPi0TCM->Eval(ratioOfEnergies->GetBinCenter(i));
      ratioOfEnergies3->SetBinContent(i,temp);
      ratioOfEnergies3->SetBinError(i,0.);
    }

    for (Int_t i=1; i<=ratioOfEnergiesToPythia3->GetNbinsX(); i++) {
      Double_t temp = histoPythia8InvXSection7TeV->GetBinContent(i)/histoPythia8InvXSection2760GeV->GetBinContent(i);
      ratioOfEnergiesToPythia3->SetBinContent(i,temp);
      ratioOfEnergiesToPythia3->SetBinError(i,0.);
    }

    ratioOfEnergies3->DrawCopy("p,same");
    ratioOfEnergiesToPythia3->DrawCopy("p,same");

    TLegend* legendRatios  = GetAndSetLegend2(0.15, 0.75, 0.4, 0.95, 0.85*textSizeLabelsPixel);
    legendRatios->SetMargin(0.25);
    legendRatios->AddEntry((TObject*)0,"#pi^{0} ALICE - TCM fits","");
    legendRatios->AddEntry(ratioOfEnergies, "8TeV/2.76TeV","p");
    legendRatios->AddEntry(ratioOfEnergies2,"8TeV/7TeV","p");
    legendRatios->AddEntry(ratioOfEnergies3,"7TeV/2.76TeV","p");
    legendRatios->Draw();

    TLegend* legendRatiosMC  = GetAndSetLegend2(0.5, 0.75, 0.75, 0.95, 0.85*textSizeLabelsPixel);
    legendRatiosMC->SetMargin(0.25);
    legendRatiosMC->AddEntry((TObject*)0,"#pi^{0} PYTHIA Monash2013","");
    legendRatiosMC->AddEntry(ratioOfEnergiesToPythia, "8TeV/2.76TeV","p");
    legendRatiosMC->AddEntry(ratioOfEnergiesToPythia2, "8TeV/7TeV","p");
    legendRatiosMC->AddEntry(ratioOfEnergiesToPythia3, "7TeV/2.76TeV","p");
    legendRatiosMC->Draw();

    DrawGammaLines(0.33, 32. , 1., 1.,0.5, kGray+2);

    histoRatioEnergies->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/Pi0_diffEnergy_TCM_ratio.%s",outputDir.Data(), suffix.Data()));

// -----------------------------------------------------------------------------------------------------------------

    histoRatioEnergies->DrawCopy();

    // plotting data 8TeV/2.76TeV
    for (Int_t i=1; i<=ratioOfEnergies->GetNbinsX(); i++) {
      Double_t temp = fitInvXSectionPi0->Eval(ratioOfEnergies->GetBinCenter(i))/fit2760GeVPi0Tsallis->Eval(ratioOfEnergies->GetBinCenter(i));
      ratioOfEnergies->SetBinContent(i,temp);
      ratioOfEnergies->SetBinError(i,0.);
    }

    ratioOfEnergies->DrawCopy("p,same");
    ratioOfEnergiesToPythia->DrawCopy("p,same");

    // plotting data 8TeV/7TeV
    for (Int_t i=1; i<=ratioOfEnergies2->GetNbinsX(); i++) {
      Double_t temp = fitInvXSectionPi0->Eval(ratioOfEnergies->GetBinCenter(i))/fit7TeVPi0Tsallis->Eval(ratioOfEnergies->GetBinCenter(i));
      ratioOfEnergies2->SetBinContent(i,temp);
      ratioOfEnergies2->SetBinError(i,0.);
    }

    ratioOfEnergies2->DrawCopy("p,same");
    ratioOfEnergiesToPythia2->DrawCopy("p,same");

    // plotting data 7TeV/2.76TeV
    for (Int_t i=1; i<=ratioOfEnergies3->GetNbinsX(); i++) {
      Double_t temp = fit7TeVPi0Tsallis->Eval(ratioOfEnergies->GetBinCenter(i))/fit2760GeVPi0Tsallis->Eval(ratioOfEnergies->GetBinCenter(i));
      ratioOfEnergies3->SetBinContent(i,temp);
      ratioOfEnergies3->SetBinError(i,0.);
    }

    ratioOfEnergies3->DrawCopy("p,same");
    ratioOfEnergiesToPythia3->DrawCopy("p,same");

    TLegend* legendRatiosTS  = GetAndSetLegend2(0.15, 0.75, 0.4, 0.95, 0.85*textSizeLabelsPixel);
    legendRatiosTS->SetMargin(0.25);
    legendRatiosTS->AddEntry((TObject*)0,"#pi^{0} ALICE - Tsallis fits","");
    legendRatiosTS->AddEntry(ratioOfEnergies, "8TeV/2.76TeV","p");
    legendRatiosTS->AddEntry(ratioOfEnergies2,"8TeV/7TeV","p");
    legendRatiosTS->AddEntry(ratioOfEnergies3,"7TeV/2.76TeV","p");
    legendRatiosTS->Draw();

    legendRatiosMC->Draw();

    DrawGammaLines(0.33, 32. , 1., 1.,0.5, kGray+2);

    histoRatioEnergies->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/Pi0_diffEnergy_ratio.%s",outputDir.Data(), suffix.Data()));

    // -----------------------------------------------------------------------------------------------------------------

    TH2F * histoRatioEnergiesRa;
    histoRatioEnergiesRa               = new TH2F("histoRatioEnergiesRa","histoRatioEnergiesRa",1000,0.3,30.,1000,0.01,6.9    );
    SetStyleHistoTH2ForGraphs(histoRatioEnergiesRa, "#it{p}_{T} (GeV/#it{c})","ratio", 0.85*textsizeLabelsEtaToPi0, textsizeLabelsEtaToPi0,
                              0.85*textsizeLabelsEtaToPi0,textsizeLabelsEtaToPi0, 0.9, 0.65, 510, 510);
    histoRatioEnergiesRa->GetXaxis()->SetMoreLogLabels();
    histoRatioEnergiesRa->GetXaxis()->SetNoExponent(kTRUE);
    histoRatioEnergiesRa->DrawCopy();

        Bool_t doMaterialError = kFALSE;

        TGraphErrors* graphRatioBinByBin7000_2760AStat = NULL;
        TGraphErrors* graphRatioBinByBin7000_2760ASys  = NULL;
        TGraphErrors* graphRatioBinByBin7000_2760BStat = NULL;
        TGraphErrors* graphRatioBinByBin7000_2760BSys  = NULL;
        TGraphAsymmErrors* AAinputASys = (TGraphAsymmErrors*) graph7TeVPi0Sys->Clone();
        TGraphAsymmErrors* AAinputBSys = (TGraphAsymmErrors*) graph2760GeVPi0Sys->Clone();
        if(doMaterialError){
          for(Int_t i=0; AAinputASys->GetX()[i]<=0.8;i++){
            AAinputASys->GetEYlow()[i] = TMath::Sqrt(TMath::Power(AAinputASys->GetEYlow()[i],2)-TMath::Power((9./100)*AAinputASys->GetY()[i],2));
            AAinputASys->GetEYhigh()[i] = TMath::Sqrt(TMath::Power(AAinputASys->GetEYhigh()[i],2)-TMath::Power((9./100)*AAinputASys->GetY()[i],2));
          }

          for(Int_t i=0; AAinputBSys->GetX()[i]<=0.8;i++){
            AAinputBSys->GetEYlow()[i] = TMath::Sqrt(TMath::Power(AAinputBSys->GetEYlow()[i],2)-TMath::Power((9./100)*AAinputBSys->GetY()[i],2));
            AAinputBSys->GetEYhigh()[i] = TMath::Sqrt(TMath::Power(AAinputBSys->GetEYhigh()[i],2)-TMath::Power((9./100)*AAinputBSys->GetY()[i],2));
          }
        }
        TGraphErrors* graphRatioBinByBin7000_2760 = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                                                            graph7TeVPi0Stat->Clone(), AAinputASys,
                                                                                                            graph2760GeVPi0Stat->Clone(), AAinputBSys,
                                                                                                            kTRUE,  kTRUE,
                                                                                                            &graphRatioBinByBin7000_2760AStat, &graphRatioBinByBin7000_2760ASys,
                                                                                                            &graphRatioBinByBin7000_2760BStat, &graphRatioBinByBin7000_2760BSys )    ;
        TGraphAsymmErrors* AAinputASys_ForStat = (TGraphAsymmErrors*) graph7TeVPi0Sys->Clone();
        TGraphAsymmErrors* AAinputBSys_ForStat = (TGraphAsymmErrors*) graph2760GeVPi0Sys->Clone();
          for(Int_t i=0; i<AAinputASys_ForStat->GetN();i++){
            AAinputASys_ForStat->GetEYlow()[i] = (1./1000.)*AAinputASys_ForStat->GetY()[i];
            AAinputASys_ForStat->GetEYhigh()[i] = (1./1000.)*AAinputASys_ForStat->GetY()[i];
          }

          for(Int_t i=0; i<AAinputBSys_ForStat->GetN();i++){
            AAinputBSys_ForStat->GetEYlow()[i] = (1./1000.)*AAinputBSys_ForStat->GetY()[i];
            AAinputBSys_ForStat->GetEYhigh()[i] = (1./1000.)*AAinputBSys_ForStat->GetY()[i];
          }

        TGraphErrors* graphRatioBinByBin7000_2760stat = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                                                            graph7TeVPi0Stat->Clone(), AAinputASys_ForStat,
                                                                                                            graph2760GeVPi0Stat->Clone(), AAinputBSys_ForStat,
                                                                                                            kTRUE,  kTRUE,
                                                                                                            &graphRatioBinByBin7000_2760AStat, &graphRatioBinByBin7000_2760ASys,
                                                                                                            &graphRatioBinByBin7000_2760BStat, &graphRatioBinByBin7000_2760BSys )    ;


        TGraphErrors* graphRatioBinByBin8000_2760AStat = NULL;
        TGraphErrors* graphRatioBinByBin8000_2760ASys  = NULL;
        TGraphErrors* graphRatioBinByBin8000_2760BStat = NULL;
        TGraphErrors* graphRatioBinByBin8000_2760BSys  = NULL;
        TGraphAsymmErrors* AAAinputASys = (TGraphAsymmErrors*) graphCombPi0InvXSectionSysA->Clone();
        TGraphAsymmErrors* AAAinputBSys = (TGraphAsymmErrors*) graph2760GeVPi0Sys->Clone();
        if(doMaterialError){
          for(Int_t i=0; AAAinputASys->GetX()[i]<=0.8;i++){
            AAAinputASys->GetEYlow()[i] = TMath::Sqrt(TMath::Power(AAAinputASys->GetEYlow()[i],2)-TMath::Power((9./100)*AAAinputASys->GetY()[i],2));
            AAAinputASys->GetEYhigh()[i] = TMath::Sqrt(TMath::Power(AAAinputASys->GetEYhigh()[i],2)-TMath::Power((9./100)*AAAinputASys->GetY()[i],2));
          }

          for(Int_t i=0; AAAinputBSys->GetX()[i]<=0.8;i++){
            AAAinputBSys->GetEYlow()[i] = TMath::Sqrt(TMath::Power(AAAinputBSys->GetEYlow()[i],2)-TMath::Power((9./100)*AAAinputBSys->GetY()[i],2));
            AAAinputBSys->GetEYhigh()[i] = TMath::Sqrt(TMath::Power(AAAinputBSys->GetEYhigh()[i],2)-TMath::Power((9./100)*AAAinputBSys->GetY()[i],2));
          }
        }

        TGraphErrors* graphRatioBinByBin8000_2760 = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                                                            graphCombPi0InvXSectionStatA->Clone(), AAAinputASys,
                                                                                                            graph2760GeVPi0Stat->Clone(), AAAinputBSys,
                                                                                                            kTRUE,  kTRUE,
                                                                                                            &graphRatioBinByBin8000_2760AStat, &graphRatioBinByBin8000_2760ASys,
                                                                                                            &graphRatioBinByBin8000_2760BStat, &graphRatioBinByBin8000_2760BSys )    ;
        TGraphAsymmErrors* AAAinputASys_ForStat = (TGraphAsymmErrors*) graphCombPi0InvXSectionSysA->Clone();
        TGraphAsymmErrors* AAAinputBSys_ForStat = (TGraphAsymmErrors*) graph2760GeVPi0Sys->Clone();
          for(Int_t i=0; i<AAAinputASys_ForStat->GetN();i++){
            AAAinputASys_ForStat->GetEYlow()[i] = (1./1000.)*AAAinputASys_ForStat->GetY()[i];
            AAAinputASys_ForStat->GetEYhigh()[i] = (1./1000.)*AAAinputASys_ForStat->GetY()[i];
          }

          for(Int_t i=0; i<AAAinputBSys_ForStat->GetN();i++){
            AAAinputBSys_ForStat->GetEYlow()[i] = (1./1000.)*AAAinputBSys_ForStat->GetY()[i];
            AAAinputBSys_ForStat->GetEYhigh()[i] = (1./1000.)*AAAinputBSys_ForStat->GetY()[i];
          }
        TGraphErrors* graphRatioBinByBin8000_2760stat = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                                                            graphCombPi0InvXSectionStatA->Clone(), AAAinputASys_ForStat,
                                                                                                            graph2760GeVPi0Stat->Clone(), AAAinputBSys_ForStat,
                                                                                                            kTRUE,  kTRUE,
                                                                                                            &graphRatioBinByBin8000_2760AStat, &graphRatioBinByBin8000_2760ASys,
                                                                                                            &graphRatioBinByBin8000_2760BStat, &graphRatioBinByBin8000_2760BSys )    ;

        TGraphErrors* graphRatioBinByBin8000_7000AStat = NULL;
        TGraphErrors* graphRatioBinByBin8000_7000ASys  = NULL;
        TGraphErrors* graphRatioBinByBin8000_7000BStat = NULL;
        TGraphErrors* graphRatioBinByBin8000_7000BSys  = NULL;
        TGraphAsymmErrors* AAAAinputASys = (TGraphAsymmErrors*) graphCombPi0InvXSectionSysA->Clone();
        TGraphAsymmErrors* AAAAinputBSys = (TGraphAsymmErrors*) graph7TeVPi0Sys->Clone();
        if(doMaterialError){
          for(Int_t i=0; AAAAinputASys->GetX()[i]<=0.8;i++){
            AAAAinputASys->GetEYlow()[i] = TMath::Sqrt(TMath::Power(AAAAinputASys->GetEYlow()[i],2)-TMath::Power((9./100)*AAAAinputASys->GetY()[i],2));
            AAAAinputASys->GetEYhigh()[i] = TMath::Sqrt(TMath::Power(AAAAinputASys->GetEYhigh()[i],2)-TMath::Power((9./100)*AAAAinputASys->GetY()[i],2));
          }

          for(Int_t i=0; AAAAinputBSys->GetX()[i]<=0.8;i++){
            AAAAinputBSys->GetEYlow()[i] = TMath::Sqrt(TMath::Power(AAAAinputBSys->GetEYlow()[i],2)-TMath::Power((9./100)*AAAAinputBSys->GetY()[i],2));
            AAAAinputBSys->GetEYhigh()[i] = TMath::Sqrt(TMath::Power(AAAAinputBSys->GetEYhigh()[i],2)-TMath::Power((9./100)*AAAAinputBSys->GetY()[i],2));
          }
        }
        TGraphErrors* graphRatioBinByBin8000_7000 = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                                                            graphCombPi0InvXSectionStatA->Clone(), AAAAinputASys,
                                                                                                            graph7TeVPi0Stat->Clone(), AAAAinputBSys,
                                                                                                            kTRUE,  kTRUE,
                                                                                                            &graphRatioBinByBin8000_7000AStat, &graphRatioBinByBin8000_7000ASys,
                                                                                                            &graphRatioBinByBin8000_7000BStat, &graphRatioBinByBin8000_7000BSys )    ;
        TGraphAsymmErrors* AAAAinputASys_ForStat = (TGraphAsymmErrors*) graphCombPi0InvXSectionSysA->Clone();
        TGraphAsymmErrors* AAAAinputBSys_ForStat = (TGraphAsymmErrors*) graph7TeVPi0Sys->Clone();
          for(Int_t i=0; i<AAAAinputASys_ForStat->GetN();i++){
            AAAAinputASys_ForStat->GetEYlow()[i] = (1./1000.)*AAAAinputASys_ForStat->GetY()[i];
            AAAAinputASys_ForStat->GetEYhigh()[i] = (1./1000.)*AAAAinputASys_ForStat->GetY()[i];
          }

          for(Int_t i=0; i<AAAAinputBSys_ForStat->GetN();i++){
            AAAAinputBSys_ForStat->GetEYlow()[i] = (1./1000.)*AAAAinputBSys_ForStat->GetY()[i];
            AAAAinputBSys_ForStat->GetEYhigh()[i] = (1./1000.)*AAAAinputBSys_ForStat->GetY()[i];
          }
        TGraphErrors* graphRatioBinByBin8000_7000stat = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                                                            graphCombPi0InvXSectionStatA->Clone(), AAAAinputASys_ForStat,
                                                                                                            graph7TeVPi0Stat->Clone(), AAAAinputBSys_ForStat,
                                                                                                            kTRUE,  kTRUE,
                                                                                                            &graphRatioBinByBin8000_7000AStat, &graphRatioBinByBin8000_7000ASys,
                                                                                                            &graphRatioBinByBin8000_7000BStat, &graphRatioBinByBin8000_7000BSys )    ;

        TBox* box = new TBox(0.3 ,1. , 0.8, 1.2);
        box->SetLineColor(kBlue-6);
        box->SetFillColorAlpha(kBlue-6,0.1);
        box->Draw();

        TBox* box2 = new TBox(0.3 ,1.4 , 0.8, 1.6);
        box2->SetLineColor(kGreen-6);
        box2->SetFillColorAlpha(kGreen-6,0.1);
        box2->Draw();

        TBox* box3 = new TBox(0.3 ,1.5 , 0.8, 1.7);
        box3->SetLineColor(kRed-6);
        box3->SetFillColorAlpha(kRed-6,0.1);
        box3->Draw();

        DrawGammaLines(0.3, 0.8 , 1.5, 1.5,2, kGreen+2);
        DrawGammaLines(0.3, 0.8 , 1.1, 1.1,2, kBlue+2);
        DrawGammaLines(0.3, 0.8 , 1.6, 1.6,2, kRed+2);

        ratioOfEnergiesToPythia->DrawCopy("p,same");
        ratioOfEnergiesToPythia2->DrawCopy("p,same");
        ratioOfEnergiesToPythia3->DrawCopy("p,same");

        graphRatioBinByBin7000_2760->SetMarkerStyle(20);
        graphRatioBinByBin7000_2760->SetMarkerColor(kGreen+2);
        graphRatioBinByBin7000_2760->SetLineColor(kGreen+2);
        graphRatioBinByBin7000_2760->SetMarkerSize(2);
        graphRatioBinByBin7000_2760->Draw("p,same");
        for(Int_t i=0; i<graphRatioBinByBin7000_2760stat->GetN(); i++) graphRatioBinByBin7000_2760stat->GetEX()[i] = 0.;
        graphRatioBinByBin7000_2760stat->SetMarkerStyle(1);
        graphRatioBinByBin7000_2760stat->SetMarkerColor(kGreen+2);
        graphRatioBinByBin7000_2760stat->SetLineColor(kGreen+2);
        graphRatioBinByBin7000_2760stat->SetLineWidth(3);
        graphRatioBinByBin7000_2760stat->Draw("p,same");

        graphRatioBinByBin8000_2760->SetMarkerStyle(20);
        graphRatioBinByBin8000_2760->SetMarkerColor(kRed+2);
        graphRatioBinByBin8000_2760->SetLineColor(kRed+2);
        graphRatioBinByBin8000_2760->SetMarkerSize(2);
        graphRatioBinByBin8000_2760->Draw("p,same");
        for(Int_t i=0; i<graphRatioBinByBin8000_2760stat->GetN(); i++) graphRatioBinByBin8000_2760stat->GetEX()[i] = 0.;
        graphRatioBinByBin8000_2760stat->SetMarkerStyle(1);
        graphRatioBinByBin8000_2760stat->SetMarkerColor(kRed+2);
        graphRatioBinByBin8000_2760stat->SetLineColor(kRed+2);
        graphRatioBinByBin8000_2760stat->SetLineWidth(3);
        graphRatioBinByBin8000_2760stat->Draw("p,same");

        graphRatioBinByBin8000_7000->SetMarkerStyle(20);
        graphRatioBinByBin8000_7000->SetMarkerColor(kBlue+2);
        graphRatioBinByBin8000_7000->SetLineColor(kBlue+2);
        graphRatioBinByBin8000_7000->SetMarkerSize(2);
        graphRatioBinByBin8000_7000->Draw("p,same");
        for(Int_t i=0; i<graphRatioBinByBin8000_7000stat->GetN(); i++) graphRatioBinByBin8000_7000stat->GetEX()[i] = 0.;
        graphRatioBinByBin8000_7000stat->SetMarkerStyle(1);
        graphRatioBinByBin8000_7000stat->SetMarkerColor(kBlue+2);
        graphRatioBinByBin8000_7000stat->SetLineColor(kBlue+2);
        graphRatioBinByBin8000_7000stat->SetLineWidth(3);
        graphRatioBinByBin8000_7000stat->Draw("p,same");

        TLegend* legendRatios_2  = GetAndSetLegend2(0.15, 0.75, 0.4, 0.95, 0.85*textSizeLabelsPixel);
        legendRatios_2->SetMargin(0.25);
        legendRatios_2->AddEntry((TObject*)0,"#pi^{0} ALICE - bin-by-bin ratio","");
        legendRatios_2->AddEntry(graphRatioBinByBin8000_2760, "8TeV/2.76TeV","p");
        legendRatios_2->AddEntry(graphRatioBinByBin8000_7000,"8TeV/7TeV","p");
        legendRatios_2->AddEntry(graphRatioBinByBin7000_2760, "7TeV/2.76TeV","p");
        legendRatios_2->Draw();

        legendRatiosMC->Draw();

        DrawGammaLines(0.33, 32. , 1., 1.,0.5, kGray+2);

        histoRatioEnergiesRa->Draw("axis,same");

        canvasEtatoPi0combo->Update();
        if(doMaterialError) canvasEtatoPi0combo->SaveAs(Form("%s/Pi0_diffEnergy_ratio2_withoutMatErr.%s",outputDir.Data(), suffix.Data()));
        else canvasEtatoPi0combo->SaveAs(Form("%s/Pi0_diffEnergy_ratio2.%s",outputDir.Data(), suffix.Data()));

        histoRatioEnergiesRa->GetXaxis()->SetRangeUser(0.3,3.);
        histoRatioEnergiesRa->GetYaxis()->SetRangeUser(0.5,4.4);
        histoRatioEnergiesRa->DrawCopy();

        box->Draw();
        box2->Draw();
        box3->Draw();
        DrawGammaLines(0.3, 0.8 , 1.1, 1.1,2, kBlue+2);
        DrawGammaLines(0.3, 0.8 , 1.5, 1.5,2, kGreen+2);
        DrawGammaLines(0.3, 0.8 , 1.6, 1.6,2, kRed+2);

        ratioOfEnergiesToPythia->DrawCopy("p,same");
        ratioOfEnergiesToPythia2->DrawCopy("p,same");
        ratioOfEnergiesToPythia3->DrawCopy("p,same");

        graphRatioBinByBin7000_2760->Draw("p,same");
        graphRatioBinByBin7000_2760stat->Draw("p,same");
        graphRatioBinByBin8000_2760->Draw("p,same");
        graphRatioBinByBin8000_2760stat->Draw("p,same");
        graphRatioBinByBin8000_7000->Draw("p,same");
        graphRatioBinByBin8000_7000stat->Draw("p,same");

        legendRatios_2->Draw();

        legendRatiosMC->Draw();
        DrawGammaLines(0.3, 3. , 1., 1.,0.5, kGray+2);
        histoRatioEnergiesRa->Draw("axis,same");

        canvasEtatoPi0combo->Update();
        if(doMaterialError) canvasEtatoPi0combo->SaveAs(Form("%s/Pi0_diffEnergy_ratio2_zoom_withoutMatErr.%s",outputDir.Data(), suffix.Data()));
        else canvasEtatoPi0combo->SaveAs(Form("%s/Pi0_diffEnergy_ratio2_zoom.%s",outputDir.Data(), suffix.Data()));

        // -----------------------------------------------------------------------------------------------------------------

    TH2F * histoRatioEnergies2;
    histoRatioEnergies2               = new TH2F("histoRatioEnergies","histoRatioEnergies",1000,0.3,35.,1000,0.01,7.8    );
    SetStyleHistoTH2ForGraphs(histoRatioEnergies2, "#it{p}_{T} (GeV/#it{c})","ratio", 0.85*textsizeLabelsEtaToPi0, textsizeLabelsEtaToPi0,
                              0.85*textsizeLabelsEtaToPi0,textsizeLabelsEtaToPi0, 0.9, 0.65, 510, 510);
    histoRatioEnergies2->GetXaxis()->SetMoreLogLabels();
    histoRatioEnergies2->GetXaxis()->SetNoExponent(kTRUE);
    histoRatioEnergies2->DrawCopy();

    TH1F *ratioOfEnergiesToPythiaEta = (TH1F*)histoPythia8InvXSectionEta->Clone("ratioOfEnergiesToPythiaEta");
    ratioOfEnergiesToPythiaEta->SetTitle("");
    ratioOfEnergiesToPythiaEta->SetLineColor(kRed+2);
    ratioOfEnergiesToPythiaEta->SetMarkerColor(kRed+2);
    ratioOfEnergiesToPythiaEta->SetMarkerStyle(24);
    ratioOfEnergiesToPythiaEta->SetLineWidth(2.);
    TH1F *ratioOfEnergiesToPythia2Eta = (TH1F*)histoPythia8InvXSectionEta->Clone("ratioOfEnergiesToPythia2Eta");
    ratioOfEnergiesToPythia2Eta->SetTitle("");
    ratioOfEnergiesToPythia2Eta->SetLineColor(kBlue+2);
    ratioOfEnergiesToPythia2Eta->SetMarkerColor(kBlue+2);
    ratioOfEnergiesToPythia2Eta->SetMarkerStyle(24);
    ratioOfEnergiesToPythia2Eta->SetLineWidth(2.);
    TH1F *ratioOfEnergiesToPythia3Eta = (TH1F*)histoPythia8InvXSectionEta->Clone("ratioOfEnergiesToPythia3Eta");
    ratioOfEnergiesToPythia3Eta->SetTitle("");
    ratioOfEnergiesToPythia3Eta->SetLineColor(kGreen+2);
    ratioOfEnergiesToPythia3Eta->SetMarkerColor(kGreen+2);
    ratioOfEnergiesToPythia3Eta->SetMarkerStyle(24);
    ratioOfEnergiesToPythia3Eta->SetLineWidth(2.);

    // eta plotting
    // plotting data 8TeV/2.76TeV
    for (Int_t i=1; i<=ratioOfEnergies->GetNbinsX(); i++) {
      Double_t temp = fitInvXSectionEta->Eval(ratioOfEnergies->GetBinCenter(i))/fit2760GeVEtaTsallis->Eval(ratioOfEnergies->GetBinCenter(i));
      ratioOfEnergies->SetBinContent(i,temp);
      ratioOfEnergies->SetBinError(i,0.);
    }

    for (Int_t i=1; i<=ratioOfEnergiesToPythiaEta->GetNbinsX(); i++) {
      Double_t temp = histoPythia8InvXSectionEta->GetBinContent(i)/histoPythia8InvXSection2760GeVEta->GetBinContent(i);
      ratioOfEnergiesToPythiaEta->SetBinContent(i,temp);
      ratioOfEnergiesToPythiaEta->SetBinError(i,0.);
    }

    ratioOfEnergies->DrawCopy("p,same");
    ratioOfEnergiesToPythiaEta->DrawCopy("p,same");

    // plotting data 8TeV/7TeV
    for (Int_t i=1; i<=ratioOfEnergies2->GetNbinsX(); i++) {
      Double_t temp = fitInvXSectionEta->Eval(ratioOfEnergies->GetBinCenter(i))/fit7TeVEtaTsallis->Eval(ratioOfEnergies->GetBinCenter(i));
      ratioOfEnergies2->SetBinContent(i,temp);
      ratioOfEnergies2->SetBinError(i,0.);
    }

    for (Int_t i=1; i<=ratioOfEnergiesToPythia2Eta->GetNbinsX(); i++) {
      Double_t temp = histoPythia8InvXSectionEta->GetBinContent(i)/histoPythia8InvXSection7TeVEta->GetBinContent(i);
      ratioOfEnergiesToPythia2Eta->SetBinContent(i,temp);
      ratioOfEnergiesToPythia2Eta->SetBinError(i,0.);
    }

    ratioOfEnergies2->DrawCopy("p,same");
    ratioOfEnergiesToPythia2Eta->DrawCopy("p,same");

    // plotting data 7TeV/2.76TeV
    for (Int_t i=1; i<=ratioOfEnergies3->GetNbinsX(); i++) {
      Double_t temp = fit7TeVEtaTsallis->Eval(ratioOfEnergies->GetBinCenter(i))/fit2760GeVEtaTsallis->Eval(ratioOfEnergies->GetBinCenter(i));
      ratioOfEnergies3->SetBinContent(i,temp);
      ratioOfEnergies3->SetBinError(i,0.);
    }

    for (Int_t i=1; i<=ratioOfEnergiesToPythia3Eta->GetNbinsX(); i++) {
      Double_t temp = histoPythia8InvXSection7TeVEta->GetBinContent(i)/histoPythia8InvXSection2760GeVEta->GetBinContent(i);
      ratioOfEnergiesToPythia3Eta->SetBinContent(i,temp);
      ratioOfEnergiesToPythia3Eta->SetBinError(i,0.);
    }

    ratioOfEnergies3->DrawCopy("p,same");
    ratioOfEnergiesToPythia3Eta->DrawCopy("p,same");

    TLegend* legendRatiosEta  = GetAndSetLegend2(0.15, 0.75, 0.4, 0.95, 0.85*textSizeLabelsPixel);
    legendRatiosEta->SetMargin(0.25);
    legendRatiosEta->AddEntry((TObject*)0,"#eta ALICE - Tsallis fits","");
    legendRatiosEta->AddEntry(ratioOfEnergies, "8TeV/2.76TeV","p");
    legendRatiosEta->AddEntry(ratioOfEnergies2,"8TeV/7TeV","p");
    legendRatiosEta->AddEntry(ratioOfEnergies3,"7TeV/2.76TeV","p");
    legendRatiosEta->Draw();

    TLegend* legendRatiosEtaMC  = GetAndSetLegend2(0.5, 0.75, 0.75, 0.95, 0.85*textSizeLabelsPixel);
    legendRatiosEtaMC->SetMargin(0.25);
    legendRatiosEtaMC->AddEntry((TObject*)0,"#eta PYTHIA Monash2013","");
    legendRatiosEtaMC->AddEntry(ratioOfEnergiesToPythiaEta, "8TeV/2.76TeV","p");
    legendRatiosEtaMC->AddEntry(ratioOfEnergiesToPythia2Eta, "8TeV/7TeV","p");
    legendRatiosEtaMC->AddEntry(ratioOfEnergiesToPythia3Eta, "7TeV/2.76TeV","p");
    legendRatiosEtaMC->Draw();

    DrawGammaLines(0.33, 32. , 1., 1.,0.5, kGray+2);

    histoRatioEnergies2->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/Eta_diffEnergy_ratio.%s",outputDir.Data(), suffix.Data()));

    // -----------------------------------------------------------------------------------------------------------------

    histoRatioEnergies2->DrawCopy();

    ratioOfEnergiesToPythiaEta->DrawCopy("p,same");
    ratioOfEnergiesToPythia2Eta->DrawCopy("p,same");
    ratioOfEnergiesToPythia3Eta->DrawCopy("p,same");
graph2760GeVEtaStat->RemovePoint(0);
graph2760GeVEtaSys->RemovePoint(0);
        TGraphErrors* graphRatioBinByBin7000_2760AStatEta = NULL;
        TGraphErrors* graphRatioBinByBin7000_2760ASysEta  = NULL;
        TGraphErrors* graphRatioBinByBin7000_2760BStatEta = NULL;
        TGraphErrors* graphRatioBinByBin7000_2760BSysEta  = NULL;
        TGraphAsymmErrors* AAinputASysEta = (TGraphAsymmErrors*) graph7TeVEtaSys->Clone();
        TGraphAsymmErrors* AAinputBSysEta = (TGraphAsymmErrors*) graph2760GeVEtaSys->Clone();
        if(doMaterialError){
          for(Int_t i=0; AAinputASysEta->GetX()[i]<=0.8;i++){
            AAinputASysEta->GetEYlow()[i] = TMath::Sqrt(TMath::Power(AAinputASysEta->GetEYlow()[i],2)-TMath::Power((9./100)*AAinputASysEta->GetY()[i],2));
            AAinputASysEta->GetEYhigh()[i] = TMath::Sqrt(TMath::Power(AAinputASysEta->GetEYhigh()[i],2)-TMath::Power((9./100)*AAinputASysEta->GetY()[i],2));
          }

          for(Int_t i=0; AAinputBSysEta->GetX()[i]<=0.8;i++){
            AAinputBSysEta->GetEYlow()[i] = TMath::Sqrt(TMath::Power(AAinputBSysEta->GetEYlow()[i],2)-TMath::Power((9./100)*AAinputBSysEta->GetY()[i],2));
            AAinputBSysEta->GetEYhigh()[i] = TMath::Sqrt(TMath::Power(AAinputBSysEta->GetEYhigh()[i],2)-TMath::Power((9./100)*AAinputBSysEta->GetY()[i],2));
          }
        }
        TGraphErrors* graphRatioBinByBin7000_2760Eta = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                                                            graph7TeVEtaStat->Clone(), AAinputASysEta,
                                                                                                            graph2760GeVEtaStat->Clone(), AAinputBSysEta,
                                                                                                            kTRUE,  kTRUE,
                                                                                                            &graphRatioBinByBin7000_2760AStatEta, &graphRatioBinByBin7000_2760ASysEta,
                                                                                                            &graphRatioBinByBin7000_2760BStatEta, &graphRatioBinByBin7000_2760BSysEta )    ;

        TGraphErrors* graphRatioBinByBin8000_2760AStatEta = NULL;
        TGraphErrors* graphRatioBinByBin8000_2760ASysEta  = NULL;
        TGraphErrors* graphRatioBinByBin8000_2760BStatEta = NULL;
        TGraphErrors* graphRatioBinByBin8000_2760BSysEta  = NULL;
        TGraphAsymmErrors* AAAinputASysEta = (TGraphAsymmErrors*) graphCombEtaInvXSectionSysA->Clone();
        TGraphAsymmErrors* AAAinputBSysEta = (TGraphAsymmErrors*) graph2760GeVEtaSys->Clone();
        if(doMaterialError){
          for(Int_t i=0; AAAinputASysEta->GetX()[i]<=0.8;i++){
            AAAinputASysEta->GetEYlow()[i] = TMath::Sqrt(TMath::Power(AAAinputASysEta->GetEYlow()[i],2)-TMath::Power((9./100)*AAAinputASysEta->GetY()[i],2));
            AAAinputASysEta->GetEYhigh()[i] = TMath::Sqrt(TMath::Power(AAAinputASysEta->GetEYhigh()[i],2)-TMath::Power((9./100)*AAAinputASysEta->GetY()[i],2));
          }

          for(Int_t i=0; AAAinputBSysEta->GetX()[i]<=0.8;i++){
            AAAinputBSysEta->GetEYlow()[i] = TMath::Sqrt(TMath::Power(AAAinputBSysEta->GetEYlow()[i],2)-TMath::Power((9./100)*AAAinputBSysEta->GetY()[i],2));
            AAAinputBSysEta->GetEYhigh()[i] = TMath::Sqrt(TMath::Power(AAAinputBSysEta->GetEYhigh()[i],2)-TMath::Power((9./100)*AAAinputBSysEta->GetY()[i],2));
          }
        }

        TGraphErrors* graphRatioBinByBin8000_2760Eta = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                                                            graphCombEtaInvXSectionStatA->Clone(), AAAinputASysEta,
                                                                                                            graph2760GeVEtaStat->Clone(), AAAinputBSysEta,
                                                                                                            kTRUE,  kTRUE,
                                                                                                            &graphRatioBinByBin8000_2760AStatEta, &graphRatioBinByBin8000_2760ASysEta,
                                                                                                            &graphRatioBinByBin8000_2760BStatEta, &graphRatioBinByBin8000_2760BSysEta )    ;
graphRatioBinByBin8000_2760Eta->RemovePoint(0);
        TGraphErrors* graphRatioBinByBin8000_7000AStatEta = NULL;
        TGraphErrors* graphRatioBinByBin8000_7000ASysEta  = NULL;
        TGraphErrors* graphRatioBinByBin8000_7000BStatEta = NULL;
        TGraphErrors* graphRatioBinByBin8000_7000BSysEta  = NULL;
        TGraphAsymmErrors* AAAAinputASysEta = (TGraphAsymmErrors*) graphCombEtaInvXSectionSysA->Clone();
        TGraphAsymmErrors* AAAAinputBSysEta = (TGraphAsymmErrors*) graph7TeVEtaSys->Clone();
        if(doMaterialError){
          for(Int_t i=0; AAAAinputASysEta->GetX()[i]<=0.8;i++){
            AAAAinputASysEta->GetEYlow()[i] = TMath::Sqrt(TMath::Power(AAAAinputASysEta->GetEYlow()[i],2)-TMath::Power((9./100)*AAAAinputASysEta->GetY()[i],2));
            AAAAinputASysEta->GetEYhigh()[i] = TMath::Sqrt(TMath::Power(AAAAinputASysEta->GetEYhigh()[i],2)-TMath::Power((9./100)*AAAAinputASysEta->GetY()[i],2));
          }

          for(Int_t i=0; AAAAinputBSysEta->GetX()[i]<=0.8;i++){
            AAAAinputBSysEta->GetEYlow()[i] = TMath::Sqrt(TMath::Power(AAAAinputBSysEta->GetEYlow()[i],2)-TMath::Power((9./100)*AAAAinputBSysEta->GetY()[i],2));
            AAAAinputBSysEta->GetEYhigh()[i] = TMath::Sqrt(TMath::Power(AAAAinputBSysEta->GetEYhigh()[i],2)-TMath::Power((9./100)*AAAAinputBSysEta->GetY()[i],2));
          }
        }
        TGraphErrors* graphRatioBinByBin8000_7000Eta = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                                                            graphCombEtaInvXSectionStatA->Clone(), AAAAinputASysEta,
                                                                                                            graph7TeVEtaStat->Clone(), AAAAinputBSysEta,
                                                                                                            kTRUE,  kTRUE,
                                                                                                            &graphRatioBinByBin8000_7000AStatEta, &graphRatioBinByBin8000_7000ASysEta,
                                                                                                            &graphRatioBinByBin8000_7000BStatEta, &graphRatioBinByBin8000_7000BSysEta )    ;

        graphRatioBinByBin7000_2760Eta->SetMarkerStyle(20);
        graphRatioBinByBin7000_2760Eta->SetMarkerColor(kGreen+2);
        graphRatioBinByBin7000_2760Eta->SetLineColor(kGreen+2);
        graphRatioBinByBin7000_2760Eta->SetMarkerSize(2);
        graphRatioBinByBin7000_2760Eta->Draw("p,same");

        graphRatioBinByBin8000_2760Eta->SetMarkerStyle(20);
        graphRatioBinByBin8000_2760Eta->SetMarkerColor(kRed+2);
        graphRatioBinByBin8000_2760Eta->SetLineColor(kRed+2);
        graphRatioBinByBin8000_2760Eta->SetMarkerSize(2);
        graphRatioBinByBin8000_2760Eta->Draw("p,same");

        graphRatioBinByBin8000_7000Eta->SetMarkerStyle(20);
        graphRatioBinByBin8000_7000Eta->SetMarkerColor(kBlue+2);
        graphRatioBinByBin8000_7000Eta->SetLineColor(kBlue+2);
        graphRatioBinByBin8000_7000Eta->SetMarkerSize(2);
        graphRatioBinByBin8000_7000Eta->Draw("p,same");

        TLegend* legendRatios_2Eta  = GetAndSetLegend2(0.15, 0.75, 0.4, 0.95, 0.85*textSizeLabelsPixel);
        legendRatios_2Eta->SetMargin(0.25);
        legendRatios_2Eta->AddEntry((TObject*)0,"#eta ALICE - bin-by-bin ratio","");
        legendRatios_2Eta->AddEntry(graphRatioBinByBin8000_2760Eta, "8TeV/2.76TeV","p");
        legendRatios_2Eta->AddEntry(graphRatioBinByBin8000_7000Eta,"8TeV/7TeV","p");
        legendRatios_2Eta->AddEntry(graphRatioBinByBin7000_2760Eta, "7TeV/2.76TeV","p");
        legendRatios_2Eta->Draw();

        legendRatiosEtaMC->Draw();

        DrawGammaLines(0.33, 32. , 1., 1.,0.5, kGray+2);

        histoRatioEnergies2->Draw("axis,same");

        canvasEtatoPi0combo->Update();
        if(doMaterialError) canvasEtatoPi0combo->SaveAs(Form("%s/Eta_diffEnergy_ratio2_withoutMatErr.%s",outputDir.Data(), suffix.Data()));
        else canvasEtatoPi0combo->SaveAs(Form("%s/Eta_diffEnergy_ratio2.%s",outputDir.Data(), suffix.Data()));

        histoRatioEnergiesRa->GetXaxis()->SetRangeUser(0.3,3.);
        histoRatioEnergiesRa->GetYaxis()->SetRangeUser(0.5,4.4);
        histoRatioEnergiesRa->DrawCopy();

        ratioOfEnergiesToPythiaEta->DrawCopy("p,same");
        ratioOfEnergiesToPythia2Eta->DrawCopy("p,same");
        ratioOfEnergiesToPythia3Eta->DrawCopy("p,same");

        graphRatioBinByBin7000_2760Eta->Draw("p,same");
        graphRatioBinByBin8000_2760Eta->Draw("p,same");
        graphRatioBinByBin8000_7000Eta->Draw("p,same");

        legendRatios_2Eta->Draw();

        legendRatiosEtaMC->Draw();
        DrawGammaLines(0.3, 3. , 1., 1.,0.5, kGray+2);
        histoRatioEnergiesRa->Draw("axis,same");

        canvasEtatoPi0combo->Update();
        if(doMaterialError) canvasEtatoPi0combo->SaveAs(Form("%s/Eta_diffEnergy_ratio2_zoom_withoutMatErr.%s",outputDir.Data(), suffix.Data()));
        else canvasEtatoPi0combo->SaveAs(Form("%s/Eta_diffEnergy_ratio2_zoom.%s",outputDir.Data(), suffix.Data()));
// -----------------------------------------------------------------------------------------------------------------


    histoRatioEnergies->DrawCopy();

    // plotting data 8TeV/2.76TeV
    TH1F *ratioOfEnergiesChPionToPythia = (TH1F*)histoPythia8InvXSectionChPion->Clone("ratioOfEnergiesChPionToPythia");
    ratioOfEnergiesChPionToPythia->SetTitle("");
    ratioOfEnergiesChPionToPythia->SetLineColor(kRed+2);
    ratioOfEnergiesChPionToPythia->SetMarkerColor(kRed+2);
    ratioOfEnergiesChPionToPythia->SetMarkerStyle(24);
    ratioOfEnergiesChPionToPythia->SetLineWidth(2.);

    for (Int_t i=1; i<=ratioOfEnergiesChPionToPythia->GetNbinsX(); i++) {
      Double_t temp = histoPythia8InvXSectionChPion->GetBinContent(i)/histoPythia8InvXSection2760GeVChPion->GetBinContent(i);
      ratioOfEnergiesChPionToPythia->SetBinContent(i,temp);
      ratioOfEnergiesChPionToPythia->SetBinError(i,0.);
    }

    ratioOfEnergiesChPionToPythia->DrawCopy("p,same");

    // plotting data 8TeV/7TeV

    TH1F *ratioOfEnergiesChPionToPythia2 = (TH1F*)histoPythia8InvXSectionChPion->Clone("ratioOfEnergiesChPionToPythia2");
    ratioOfEnergiesChPionToPythia2->SetTitle("");
    ratioOfEnergiesChPionToPythia2->SetLineColor(kBlue+2);
    ratioOfEnergiesChPionToPythia2->SetMarkerColor(kBlue+2);
    ratioOfEnergiesChPionToPythia2->SetMarkerStyle(24);
    ratioOfEnergiesChPionToPythia2->SetLineWidth(2.);

    for (Int_t i=1; i<=ratioOfEnergiesChPionToPythia2->GetNbinsX(); i++) {
      Double_t temp = histoPythia8InvXSectionChPion->GetBinContent(i)/histoPythia8InvXSection7TeVChPion->GetBinContent(i);
      ratioOfEnergiesChPionToPythia2->SetBinContent(i,temp);
      ratioOfEnergiesChPionToPythia2->SetBinError(i,0.);
    }

    ratioOfEnergiesChPionToPythia2->DrawCopy("p,same");

    // plotting data 7TeV/2.76TeV
    TH1F *ratioOfEnergiesChPionToPythia3 = (TH1F*)histoPythia8InvXSectionChPion->Clone("ratioOfEnergiesChPionToPythia3");
    ratioOfEnergiesChPionToPythia3->SetTitle("");
    ratioOfEnergiesChPionToPythia3->SetLineColor(kGreen+2);
    ratioOfEnergiesChPionToPythia3->SetMarkerColor(kGreen+2);
    ratioOfEnergiesChPionToPythia3->SetMarkerStyle(24);
    ratioOfEnergiesChPionToPythia3->SetLineWidth(2.);

    for (Int_t i=1; i<=ratioOfEnergiesChPionToPythia3->GetNbinsX(); i++) {
      Double_t temp = histoPythia8InvXSection7TeVChPion->GetBinContent(i)/histoPythia8InvXSection2760GeVChPion->GetBinContent(i);
      ratioOfEnergiesChPionToPythia3->SetBinContent(i,temp);
      ratioOfEnergiesChPionToPythia3->SetBinError(i,0.);
    }

    ratioOfEnergiesChPionToPythia3->DrawCopy("p,same");

    TLegend* legendRatiosMCChPion  = GetAndSetLegend2(0.5, 0.75, 0.75, 0.95, 0.85*textSizeLabelsPixel);
    legendRatiosMCChPion->SetMargin(0.25);
    legendRatiosMCChPion->AddEntry((TObject*)0,"#pi^{+/-} PYTHIA Monash2013","");
    legendRatiosMCChPion->AddEntry(ratioOfEnergiesChPionToPythia, "8TeV/2.76TeV","p");
    legendRatiosMCChPion->AddEntry(ratioOfEnergiesChPionToPythia2, "8TeV/7TeV","p");
    legendRatiosMCChPion->AddEntry(ratioOfEnergiesChPionToPythia3, "7TeV/2.76TeV","p");
    legendRatiosMCChPion->Draw();

    TGraphErrors* graphRatioBinByBin7000_2760ChPion = CalculateRatioBetweenSpectraWithDifferentBinning(
                                                                                                        histoChPion7TeVStat->Clone(), histoChPion7TeVSys->Clone(),
                                                                                                        histoChPion2760GeVStat->Clone(), histoChPion2760GeVSys->Clone(),
                                                                                                        kTRUE,  kTRUE,
                                                                                                        &graphRatioBinByBin8000_2760AStat, &graphRatioBinByBin8000_2760ASys,
                                                                                                        &graphRatioBinByBin8000_2760BStat, &graphRatioBinByBin8000_2760BSys )    ;

    graphRatioBinByBin7000_2760ChPion->SetMarkerStyle(20);
    graphRatioBinByBin7000_2760ChPion->SetMarkerColor(kGreen+2);
    graphRatioBinByBin7000_2760ChPion->SetLineColor(kGreen+2);
    graphRatioBinByBin7000_2760ChPion->SetMarkerSize(2);
    graphRatioBinByBin7000_2760ChPion->Draw("p,same");

    TLegend* legendRatiosChPion  = GetAndSetLegend2(0.15, 0.75, 0.4, 0.95, 0.85*textSizeLabelsPixel);
    legendRatiosChPion->SetMargin(0.25);
    legendRatiosChPion->AddEntry((TObject*)0,"#pi^{+/-} ALICE - bin-by-bin ratio","");
    legendRatiosChPion->AddEntry((TObject*)0, "N/A","p");
    legendRatiosChPion->AddEntry((TObject*)0, "N/A","p");
//    legendRatiosChPion->AddEntry(graphRatioBinByBin8000_2760ChPion, "8TeV/2.76TeV","p");
//    legendRatiosChPion->AddEntry(graphRatioBinByBin8000_7000ChPion,"8TeV/7TeV","p");
    legendRatiosChPion->AddEntry(graphRatioBinByBin7000_2760ChPion, "7TeV/2.76TeV","p");
    legendRatiosChPion->Draw();

    DrawGammaLines(0.33, 32. , 1., 1.,0.5, kGray+2);

    histoRatioEnergies->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/ChPion_diffEnergy_ratio.%s",outputDir.Data(), suffix.Data()));

    // -----------------------------------------------------------------------------------------------------------------

    histoRatioEnergies->DrawCopy();

    ratioOfEnergiesToPythia->SetMarkerStyle(28);
    ratioOfEnergiesToPythia->SetMarkerSize(2);
    ratioOfEnergiesToPythia->DrawCopy("p,same");
    ratioOfEnergiesToPythia2->SetMarkerStyle(28);
    ratioOfEnergiesToPythia2->SetMarkerSize(2);
    ratioOfEnergiesToPythia2->DrawCopy("p,same");
    ratioOfEnergiesToPythia3->SetMarkerStyle(28);
    ratioOfEnergiesToPythia3->SetMarkerSize(2);
    ratioOfEnergiesToPythia3->DrawCopy("p,same");
    legendRatiosMC->SetX1(legendRatios->GetX1());
    legendRatiosMC->SetX2(legendRatios->GetX2());
    legendRatiosMC->SetY1(legendRatios->GetY1());
    legendRatiosMC->SetY2(legendRatios->GetY2());
    legendRatiosMC->Draw();

    ratioOfEnergiesChPionToPythia->SetMarkerStyle(20);
    ratioOfEnergiesChPionToPythia->DrawCopy("p,same");
    ratioOfEnergiesChPionToPythia2->SetMarkerStyle(20);
    ratioOfEnergiesChPionToPythia2->DrawCopy("p,same");
    ratioOfEnergiesChPionToPythia3->SetMarkerStyle(20);
    ratioOfEnergiesChPionToPythia3->DrawCopy("p,same");
    legendRatiosMCChPion->Draw();


    DrawGammaLines(0.3, 32. , 1., 1.,0.5, kGray+2);

    histoRatioEnergies->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/MC_diffEnergy_ratio.%s",outputDir.Data(), suffix.Data()));


    histoRatioEnergiesRa->DrawCopy();

    ratioOfEnergiesToPythia->DrawCopy("p,same");
    ratioOfEnergiesToPythia2->DrawCopy("p,same");
    ratioOfEnergiesToPythia3->DrawCopy("p,same");
    legendRatiosMC->Draw();

    ratioOfEnergiesChPionToPythia->DrawCopy("p,same");
    ratioOfEnergiesChPionToPythia2->DrawCopy("p,same");
    ratioOfEnergiesChPionToPythia3->DrawCopy("p,same");
    legendRatiosMCChPion->Draw();

    DrawGammaLines(0.3, 3. , 1., 1.,0.5, kGray+2);

    histoRatioEnergiesRa->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/MC_diffEnergy_ratio_zoom.%s",outputDir.Data(), suffix.Data()));

    // **********************************************************************************************************************
    // **************************Fit ATLAS charged particles and plot ratio to fit ******************************************
    // **********************************************************************************************************************

    Double_t paramGraphChPi[3]              = {5e11, 6., 0.13};
    TF1* fitInvXSectionChPion7TeV           = FitObject("l","fitInvCrossSectionChPion7TeV","ChargedPi",histoChPion7TeV,0.1,20.,paramGraphChPi,"QNRMEX0+");

    Double_t paramGraphChPi2[3]              = {5e11, 6., 0.13};
    TF1* fitInvXSectionChPion2760GeV           = FitObject("l","fitInvCrossSectionChPion2760GeV","ChargedPi",histoChPion2760GeV,0.1,20.,paramGraphChPi2,"QNRMEX0+");

//    Double_t paramChPion[5]  = { 5e11,0.1,
//                                             5e10,0.6,3.0};
//    TF1* fitInvXSectionChPion7TeV    = FitObject("tcm","fitTCMInvCrossSectionChPion7TeV","ChargedPi",histoChPion7TeV,0.1,20. ,paramChPion,"QNRMEX0+","", kFALSE);

//    Double_t paramChPion2[5]  = { 5e11,0.1,
//                                             5e10,0.6,3.0};
//    TF1* fitInvXSectionChPion2760GeV    = FitObject("tcm","fitTCMInvCrossSectionChPion2760GeV","ChargedPi",histoChPion2760GeV,0.1,20. ,paramChPion2,"QNRMEX0+","", kFALSE);

    cout << "fit charged pions 7 TeV:" << endl;
    cout << WriteParameterToFile(fitInvXSectionChPion7TeV) << endl;
    cout << "fit charged pions 2760 GeV:" << endl;
    cout << WriteParameterToFile(fitInvXSectionChPion2760GeV) << endl;

    canvasDummyATLAS->cd();

    TH2F* histo2DDummyALICE;
    histo2DDummyALICE               = new TH2F("histo2DDummyALICE","histo2DDummyALICE",1000,0.01,20.,1000,1,8e11);
    SetStyleHistoTH2ForGraphs(histo2DDummyALICE, "", 0.032,0.04, 0.04,0.04, 0.8,1.55);
    histo2DDummyALICE->DrawCopy();

    histoChPion7TeVStat->Draw("psame");
    histoChPion2760GeVStat->Draw("psame");

    fitInvXSectionChPion7TeV->SetLineColor(kRed+2);
    fitInvXSectionChPion7TeV->Draw("same");

    fitInvXSectionChPion2760GeV->SetLineColor(kBlue+2);
    fitInvXSectionChPion2760GeV->Draw("same");

    cout << "Integration 7 TeV: 0-20GeV/c (" << fitInvXSectionChPion7TeV->Integral(0.,20.) << "), 0.4-20GeV/c (" << fitInvXSectionChPion7TeV->Integral(0.4,20.) <<") - ratio: " << fitInvXSectionChPion7TeV->Integral(0.,20.)/fitInvXSectionChPion7TeV->Integral(0.4,20.) << endl;
    cout << "Integration 2.76 TeV: 0-20GeV/c (" << fitInvXSectionChPion2760GeV->Integral(0.,20.) << "), 0.4-20GeV/c (" << fitInvXSectionChPion2760GeV->Integral(0.4,20.) <<") - ratio: " << fitInvXSectionChPion2760GeV->Integral(0.,20.)/fitInvXSectionChPion2760GeV->Integral(0.4,20.) << endl;

    cout << "Integration 7 TeV: 0.1-20GeV/c (" << fitInvXSectionChPion7TeV->Integral(0.1,20.) << "), 0.4-20GeV/c (" << fitInvXSectionChPion7TeV->Integral(0.4,20.) <<") - ratio: " << fitInvXSectionChPion7TeV->Integral(0.1,20.)/fitInvXSectionChPion7TeV->Integral(0.4,20.) << endl;
    cout << "Integration 2.76 TeV: 0.1-20GeV/c (" << fitInvXSectionChPion2760GeV->Integral(0.1,20.) << "), 0.4-20GeV/c (" << fitInvXSectionChPion2760GeV->Integral(0.4,20.) <<") - ratio: " << fitInvXSectionChPion2760GeV->Integral(0.1,20.)/fitInvXSectionChPion2760GeV->Integral(0.4,20.) << endl;

    cout << "Integration Histo 7 TeV: 0.1-20GeV/c (" << histoChPion7TeV->Integral(histoChPion7TeV->FindBin(0.11),histoChPion7TeV->FindBin(19.9)) << "), 0.4-20GeV/c (" << histoChPion7TeV->Integral(histoChPion7TeV->FindBin(0.41),histoChPion7TeV->FindBin(19.9)) <<") - ratio: " << histoChPion7TeV->Integral(histoChPion7TeV->FindBin(0.11),histoChPion7TeV->FindBin(19.9))/histoChPion7TeV->Integral(histoChPion7TeV->FindBin(0.41),histoChPion7TeV->FindBin(19.9)) << endl;
    cout << "Integration Histo 2.76 TeV: 0.1-20GeV/c (" << histoChPion2760GeV->Integral(histoChPion2760GeV->FindBin(0.11),histoChPion2760GeV->FindBin(19.9)) << "), 0.4-20GeV/c (" << histoChPion2760GeV->Integral(histoChPion2760GeV->FindBin(0.41),histoChPion2760GeV->FindBin(19.9)) <<") - ratio: " << histoChPion2760GeV->Integral(histoChPion2760GeV->FindBin(0.11),histoChPion2760GeV->FindBin(19.9))/histoChPion2760GeV->Integral(histoChPion2760GeV->FindBin(0.41),histoChPion2760GeV->FindBin(19.9)) << endl;

    canvasDummyATLAS->Update();
    canvasDummyATLAS->Print(Form("%s/ChPion_ComparisonWithFit.%s",outputDir.Data(),suffix.Data()));

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
    Color_t markerColorInvMassMBG1      = kGray+3;
    Color_t markerColorInvMassMBG2      = kGray+1;
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
                DrawGammaSetMarker(histoPi0InvMassBGMixedPCM[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG1, markerColorInvMassMBG1);
                histoPi0InvMassBGMixedPCM[i]->SetLineWidth(markerSizeInvMassMBG);
                histoPi0InvMassBGMixedPCM[i]->Draw("same");
                DrawGammaSetMarker(histoPi0InvMassRemBGPCM[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG2, markerColorInvMassMBG2);
                histoPi0InvMassRemBGPCM[i]->SetLineWidth(markerSizeInvMassMBG);
                histoPi0InvMassRemBGPCM[i]->Draw("same");

                DrawGammaSetMarker(histoPi0InvMassSigRemBGSubPCM[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
                histoPi0InvMassSigRemBGSubPCM[i]->Draw("same");
                fitPi0InvMassSigPCM[i]->SetLineWidth(1);
                fitPi0InvMassSigPCM[i]->SetNpx(1000);
                fitPi0InvMassSigPCM[i]->SetRange(0,0.255);
                fitPi0InvMassSigPCM[i]->SetLineColor(fitColorInvMassSG);
                fitPi0InvMassSigPCM[i]->Draw("same");

                TLatex *labelALICE      = new TLatex(0.135,0.9,ALICEperfor.Data());
                SetStyleTLatex( labelALICE, 0.85*textSizeLabelsPixel,4);
                labelALICE->SetTextFont(43);
                labelALICE->Draw();

                TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9-0.9*0.8*textsizeLabelsPP,collisionSystem8TeV.Data());
                SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
                labelInvMassEnergy->SetTextFont(43);
                labelInvMassEnergy->Draw();

                TLatex *labelInvMassTrigger      = new TLatex(0.135,0.9-0.9*2*0.8*textsizeLabelsPP,Form("%s",nameTriggerAlternative[i].Data()));
                SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4);
                labelInvMassTrigger->SetTextFont(43);
                labelInvMassTrigger->Draw();

                TLatex *labelInvMassRecoPCM  = new TLatex(0.135,0.9-0.9*3*0.8*textsizeLabelsPP,"PCM");
                SetStyleTLatex( labelInvMassRecoPCM, 0.85*textSizeLabelsPixel,4);
                labelInvMassRecoPCM->SetTextFont(43);
                labelInvMassRecoPCM->Draw();

                if(plotDate){
                  TLatex *labelDate = new TLatex(0.135, 0.9-0.9*4*0.8*textsizeLabelsPP,date.Data());
                  SetStyleTLatex( labelDate, 0.85*textSizeLabelsPixel,4);
                  labelDate->SetTextFont(43);
                  labelDate->Draw();
                }

                SetStyleTLatex( labelInvMassPtRangePCM, 0.85*textSizeLabelsPixel,4);
                labelInvMassPtRangePCM->SetTextAlign(31);
                labelInvMassPtRangePCM->SetTextFont(43);
                labelInvMassPtRangePCM->Draw();

                TLegend* legendInvMassPCM  = GetAndSetLegend2(0.67, 0.88-5*0.8*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
                legendInvMassPCM->SetMargin(0.25);
                legendInvMassPCM->AddEntry(histoPi0InvMassSigPlusBGPCM[i],"Raw real events","l");
                legendInvMassPCM->AddEntry(histoPi0InvMassBGMixedPCM[i],"Mixed event BG","p");
                legendInvMassPCM->AddEntry(histoPi0InvMassRemBGPCM[i],"Remain. BG","p");
                legendInvMassPCM->AddEntry(histoPi0InvMassSigRemBGSubPCM[i],"BG subtracted","p");
                legendInvMassPCM->AddEntry(fitPi0InvMassSigPCM[i], "Fit","l");
                legendInvMassPCM->Draw();
                canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBinPCM_%s.%s",outputDir.Data(), nameTrigger[i].Data(), suffix.Data()));
            } else {
                cout << "missing partial input for invariant mass bin for PCM for trigger: " << nameTrigger[i].Data() << endl;
            }

            if (haveAllEtaInvMassPCM[i]){
                canvasInvMassSamplePlot->cd();
                histo2DEtaInvMassDummy->GetYaxis()->SetRangeUser(histoEtaInvMassSigRemBGSubPCM[i]->GetMinimum(),1.1*histoEtaInvMassSigPlusBGPCM[i]->GetMaximum());
                histo2DEtaInvMassDummy->DrawCopy();

                TLatex *labelInvMassPtRangePCM = new TLatex(0.945,0.9,"#eta: 1.1 GeV/#it{c} < #it{p}_{T} < 1.4 GeV/#it{c}");

                DrawGammaSetMarker(histoEtaInvMassSigPlusBGPCM[i], markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
                histoEtaInvMassSigPlusBGPCM[i]->SetLineWidth(1);
                histoEtaInvMassSigPlusBGPCM[i]->Draw("hist,e,same");
                DrawGammaSetMarker(histoEtaInvMassBGMixedPCM[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG1, markerColorInvMassMBG1);
                histoEtaInvMassBGMixedPCM[i]->SetLineWidth(markerSizeInvMassMBG);
                histoEtaInvMassBGMixedPCM[i]->Draw("same");
                DrawGammaSetMarker(histoEtaInvMassRemBGPCM[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG2, markerColorInvMassMBG2);
                histoEtaInvMassRemBGPCM[i]->SetLineWidth(markerSizeInvMassMBG);
                histoEtaInvMassRemBGPCM[i]->Draw("same");

                DrawGammaSetMarker(histoEtaInvMassSigRemBGSubPCM[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
                histoEtaInvMassSigRemBGSubPCM[i]->Draw("same");
                fitEtaInvMassSigPCM[i]->SetLineWidth(1);
                fitEtaInvMassSigPCM[i]->SetNpx(1000);
                fitEtaInvMassSigPCM[i]->SetRange(0.35,0.695);
                fitEtaInvMassSigPCM[i]->SetLineColor(fitColorInvMassSG);
                fitEtaInvMassSigPCM[i]->Draw("same");

                TLatex *labelALICE      = new TLatex(0.135,0.9,ALICEperfor.Data());
                SetStyleTLatex( labelALICE, 0.85*textSizeLabelsPixel,4);
                labelALICE->SetTextFont(43);
                labelALICE->Draw();

                TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9-0.9*0.8*textsizeLabelsPP,collisionSystem8TeV.Data());
                SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
                labelInvMassEnergy->SetTextFont(43);
                labelInvMassEnergy->Draw();

                TLatex *labelInvMassTrigger      = new TLatex(0.135,0.9-0.9*2*0.8*textsizeLabelsPP,Form("%s",nameTriggerAlternative[i].Data()));
                SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4);
                labelInvMassTrigger->SetTextFont(43);
                labelInvMassTrigger->Draw();

                TLatex *labelInvMassRecoPCM  = new TLatex(0.135,0.9-0.9*3*0.8*textsizeLabelsPP,"PCM");
                SetStyleTLatex( labelInvMassRecoPCM, 0.85*textSizeLabelsPixel,4);
                labelInvMassRecoPCM->SetTextFont(43);
                labelInvMassRecoPCM->Draw();

                if(plotDate){
                  TLatex *labelDate = new TLatex(0.135, 0.9-0.9*4*0.8*textsizeLabelsPP,date.Data());
                  SetStyleTLatex( labelDate, 0.85*textSizeLabelsPixel,4);
                  labelDate->SetTextFont(43);
                  labelDate->Draw();
                }

                SetStyleTLatex( labelInvMassPtRangePCM, 0.85*textSizeLabelsPixel,4);
                labelInvMassPtRangePCM->SetTextAlign(31);
                labelInvMassPtRangePCM->SetTextFont(43);
                labelInvMassPtRangePCM->Draw();

                TLegend* legendInvMassPCM  = GetAndSetLegend2(0.67, 0.88-5*0.8*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
                legendInvMassPCM->SetMargin(0.25);
                legendInvMassPCM->AddEntry(histoEtaInvMassSigPlusBGPCM[i],"Raw real events","l");
                legendInvMassPCM->AddEntry(histoEtaInvMassBGMixedPCM[i],"Mixed event BG","p");
                legendInvMassPCM->AddEntry(histoEtaInvMassRemBGPCM[i],"Remain. BG","p");
                legendInvMassPCM->AddEntry(histoEtaInvMassSigRemBGSubPCM[i],"BG subtracted","p");
                legendInvMassPCM->AddEntry(fitEtaInvMassSigPCM[i], "Fit","l");
                legendInvMassPCM->Draw();
                canvasInvMassSamplePlot->SaveAs(Form("%s/Eta_InvMassBinPCM_%s.%s",outputDir.Data(), nameTrigger[i].Data(), suffix.Data()));
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
                DrawGammaSetMarker(histoPi0InvMassBGMixedPCMEMCAL[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG1, markerColorInvMassMBG1);
                histoPi0InvMassBGMixedPCMEMCAL[i]->SetLineWidth(markerSizeInvMassMBG);
                histoPi0InvMassBGMixedPCMEMCAL[i]->Draw("same");
                DrawGammaSetMarker(histoPi0InvMassRemBGPCMEMCAL[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG2, markerColorInvMassMBG2);
                histoPi0InvMassRemBGPCMEMCAL[i]->SetLineWidth(markerSizeInvMassMBG);
                histoPi0InvMassRemBGPCMEMCAL[i]->Draw("same");

                DrawGammaSetMarker(histoPi0InvMassSigRemBGSubPCMEMCAL[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
                histoPi0InvMassSigRemBGSubPCMEMCAL[i]->Draw("same");
                fitPi0InvMassSigPCMEMCAL[i]->SetLineWidth(1);
                fitPi0InvMassSigPCMEMCAL[i]->SetRange(0,0.255);
                fitPi0InvMassSigPCMEMCAL[i]->SetLineColor(fitColorInvMassSG);
                fitPi0InvMassSigPCMEMCAL[i]->Draw("same");

                TLatex *labelALICE      = new TLatex(0.135,0.9,ALICEperfor.Data());
                SetStyleTLatex( labelALICE, 0.85*textSizeLabelsPixel,4);
                labelALICE->SetTextFont(43);
                labelALICE->Draw();

                TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9-0.9*0.8*textsizeLabelsPP,collisionSystem8TeV.Data());
                SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
                labelInvMassEnergy->SetTextFont(43);
                labelInvMassEnergy->Draw();

                TLatex *labelInvMassTrigger      = new TLatex(0.135,0.9-0.9*2*0.8*textsizeLabelsPP,Form("%s",nameTriggerAlternative[i].Data()));
                SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4);
                labelInvMassTrigger->SetTextFont(43);
                labelInvMassTrigger->Draw();

                TLatex *labelInvMassRecoPCMEMC  = new TLatex(0.135,0.9-0.9*3*0.8*textsizeLabelsPP,"PCM-EMC");
                SetStyleTLatex( labelInvMassRecoPCMEMC, 0.85*textSizeLabelsPixel,4);
                labelInvMassRecoPCMEMC->SetTextFont(43);
                labelInvMassRecoPCMEMC->Draw();

                if(plotDate){
                  TLatex *labelDate = new TLatex(0.135, 0.9-0.9*4*0.8*textsizeLabelsPP,date.Data());
                  SetStyleTLatex( labelDate, 0.85*textSizeLabelsPixel,4);
                  labelDate->SetTextFont(43);
                  labelDate->Draw();
                }

                SetStyleTLatex( labelInvMassPtRangePCMEMCAL, 0.85*textSizeLabelsPixel,4);
                labelInvMassPtRangePCMEMCAL->SetTextAlign(31);
                labelInvMassPtRangePCMEMCAL->SetTextFont(43);
                labelInvMassPtRangePCMEMCAL->Draw();

                TLegend* legendInvMassPCMEMCAL  = GetAndSetLegend2(0.67, 0.88-5*0.8*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
                legendInvMassPCMEMCAL->SetMargin(0.25);
                legendInvMassPCMEMCAL->AddEntry(histoPi0InvMassSigPlusBGPCMEMCAL[i],"Raw real events","l");
                legendInvMassPCMEMCAL->AddEntry(histoPi0InvMassBGMixedPCMEMCAL[i],"Mixed event BG","p");
                legendInvMassPCMEMCAL->AddEntry(histoPi0InvMassRemBGPCMEMCAL[i],"Remain. BG","p");
                legendInvMassPCMEMCAL->AddEntry(histoPi0InvMassSigRemBGSubPCMEMCAL[i],"BG subtracted","p");
                legendInvMassPCMEMCAL->AddEntry(fitPi0InvMassSigPCMEMCAL[i], "Fit","l");
                legendInvMassPCMEMCAL->Draw();
                canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBinPCMEMC_%s.%s",outputDir.Data(), nameTrigger[i].Data(), suffix.Data()));
            } else {
                cout << "missing partial input for invariant mass bin for PCM-EMC for trigger: " << nameTrigger[i].Data() << endl;
            }

            if (haveAllEtaInvMassPCMEMCAL[i]){
                canvasInvMassSamplePlot->cd();
                histo2DEtaInvMassDummy->GetYaxis()->SetRangeUser(histoEtaInvMassSigRemBGSubPCMEMCAL[i]->GetMinimum(),1.1*histoEtaInvMassSigPlusBGPCMEMCAL[i]->GetMaximum());
                histo2DEtaInvMassDummy->DrawCopy();

                TLatex *labelInvMassPtRangePCMEMCAL = new TLatex(0.945,0.9,Form("#eta: %s", histoEtaInvMassSigPCMEMCAL[i]->GetTitle()));

                DrawGammaSetMarker(histoEtaInvMassSigPlusBGPCMEMCAL[i], markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
                histoEtaInvMassSigPlusBGPCMEMCAL[i]->SetLineWidth(1);
                histoEtaInvMassSigPlusBGPCMEMCAL[i]->Draw("hist,e,same");
                DrawGammaSetMarker(histoEtaInvMassBGMixedPCMEMCAL[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG1, markerColorInvMassMBG1);
                histoEtaInvMassBGMixedPCMEMCAL[i]->SetLineWidth(markerSizeInvMassMBG);
                histoEtaInvMassBGMixedPCMEMCAL[i]->Draw("same");
                DrawGammaSetMarker(histoEtaInvMassRemBGPCMEMCAL[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG2, markerColorInvMassMBG2);
                histoEtaInvMassRemBGPCMEMCAL[i]->SetLineWidth(markerSizeInvMassMBG);
                histoEtaInvMassRemBGPCMEMCAL[i]->Draw("same");

                DrawGammaSetMarker(histoEtaInvMassSigRemBGSubPCMEMCAL[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
                histoEtaInvMassSigRemBGSubPCMEMCAL[i]->Draw("same");
                fitEtaInvMassSigPCMEMCAL[i]->SetLineWidth(1);
                fitEtaInvMassSigPCMEMCAL[i]->SetRange(0.35,0.695);
                fitEtaInvMassSigPCMEMCAL[i]->SetLineColor(fitColorInvMassSG);
                fitEtaInvMassSigPCMEMCAL[i]->Draw("same");

                TLatex *labelALICE      = new TLatex(0.135,0.9,ALICEperfor.Data());
                SetStyleTLatex( labelALICE, 0.85*textSizeLabelsPixel,4);
                labelALICE->SetTextFont(43);
                labelALICE->Draw();

                TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9-0.9*0.8*textsizeLabelsPP,collisionSystem8TeV.Data());
                SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
                labelInvMassEnergy->SetTextFont(43);
                labelInvMassEnergy->Draw();

                TLatex *labelInvMassTrigger      = new TLatex(0.135,0.9-0.9*2*0.8*textsizeLabelsPP,Form("%s",nameTriggerAlternative[i].Data()));
                SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4);
                labelInvMassTrigger->SetTextFont(43);
                labelInvMassTrigger->Draw();

                TLatex *labelInvMassRecoPCMEMC  = new TLatex(0.135,0.9-0.9*3*0.8*textsizeLabelsPP,"PCM-EMC");
                SetStyleTLatex( labelInvMassRecoPCMEMC, 0.85*textSizeLabelsPixel,4);
                labelInvMassRecoPCMEMC->SetTextFont(43);
                labelInvMassRecoPCMEMC->Draw();

                if(plotDate){
                  TLatex *labelDate = new TLatex(0.135, 0.9-0.9*4*0.8*textsizeLabelsPP,date.Data());
                  SetStyleTLatex( labelDate, 0.85*textSizeLabelsPixel,4);
                  labelDate->SetTextFont(43);
                  labelDate->Draw();
                }

                SetStyleTLatex( labelInvMassPtRangePCMEMCAL, 0.85*textSizeLabelsPixel,4);
                labelInvMassPtRangePCMEMCAL->SetTextAlign(31);
                labelInvMassPtRangePCMEMCAL->SetTextFont(43);
                labelInvMassPtRangePCMEMCAL->Draw();

                TLegend* legendInvMassPCMEMCAL  = GetAndSetLegend2(0.67, 0.88-5*0.8*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
                legendInvMassPCMEMCAL->SetMargin(0.25);
                legendInvMassPCMEMCAL->AddEntry(histoEtaInvMassSigPlusBGPCMEMCAL[i],"Raw real events","l");
                legendInvMassPCMEMCAL->AddEntry(histoEtaInvMassBGMixedPCMEMCAL[i],"Mixed event BG","p");
                legendInvMassPCMEMCAL->AddEntry(histoEtaInvMassRemBGPCMEMCAL[i],"Remain. BG","p");
                legendInvMassPCMEMCAL->AddEntry(histoEtaInvMassSigRemBGSubPCMEMCAL[i],"BG subtracted","p");
                legendInvMassPCMEMCAL->AddEntry(fitEtaInvMassSigPCMEMCAL[i], "Fit","l");
                legendInvMassPCMEMCAL->Draw();
                canvasInvMassSamplePlot->SaveAs(Form("%s/Eta_InvMassBinPCMEMC_%s.%s",outputDir.Data(), nameTrigger[i].Data(), suffix.Data()));
            } else {
                cout << "missing partial input for invariant mass bin for PCM-EMC for trigger: " << nameTrigger[i].Data() << endl;
            }


            if (haveAllPi0InvMassEMCAL[i]){
                canvasInvMassSamplePlot->cd();
                histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(0.02,0.3);
                if(i==1) histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(0.02,0.3);
                else if(i==2) histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(0.08,0.3);
                histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSubEMCAL[i]->GetMinimum(),1.25*histoPi0InvMassSigPlusBGEMCAL[i]->GetMaximum());
                histo2DPi0InvMassDummy->DrawCopy();

                TLatex *labelInvMassPtRangeEMCAL = new TLatex(0.945,0.9,Form("#pi^{0}: %s", histoPi0InvMassSigEMCAL[i]->GetTitle()));

                DrawGammaSetMarker(histoPi0InvMassSigPlusBGEMCAL[i], markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
                histoPi0InvMassSigPlusBGEMCAL[i]->SetLineWidth(1);
                histoPi0InvMassSigPlusBGEMCAL[i]->Draw("hist,e,same");
                DrawGammaSetMarker(histoPi0InvMassBGMixedEMCAL[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG1, markerColorInvMassMBG1);
                histoPi0InvMassBGMixedEMCAL[i]->SetLineWidth(markerSizeInvMassMBG);
                histoPi0InvMassBGMixedEMCAL[i]->Draw("same");
                DrawGammaSetMarker(histoPi0InvMassRemBGEMCAL[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG2, markerColorInvMassMBG2);
                histoPi0InvMassRemBGEMCAL[i]->SetLineWidth(markerSizeInvMassMBG);
                histoPi0InvMassRemBGEMCAL[i]->Draw("same");

                DrawGammaSetMarker(histoPi0InvMassSigRemBGSubEMCAL[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
                histoPi0InvMassSigRemBGSubEMCAL[i]->Draw("same");
                fitPi0InvMassSigEMCAL[i]->SetLineWidth(1);
                fitPi0InvMassSigEMCAL[i]->SetRange(0.,0.255);
                fitPi0InvMassSigEMCAL[i]->SetLineColor(fitColorInvMassSG);
                fitPi0InvMassSigEMCAL[i]->Draw("same");

                TLatex *labelALICE      = new TLatex(0.135,0.9,ALICEperfor.Data());
                SetStyleTLatex( labelALICE, 0.85*textSizeLabelsPixel,4);
                labelALICE->SetTextFont(43);
                labelALICE->Draw();

                TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9-0.9*0.8*textsizeLabelsPP,collisionSystem8TeV.Data());
                SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
                labelInvMassEnergy->SetTextFont(43);
                labelInvMassEnergy->Draw();

                TLatex *labelInvMassTrigger      = new TLatex(0.135,0.9-0.9*0.8*2*textsizeLabelsPP,Form("%s",nameTriggerAlternative[i].Data()));
                SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4);
                labelInvMassTrigger->SetTextFont(43);
                labelInvMassTrigger->Draw();

                TLatex *labelInvMassRecoEMC  = new TLatex(0.135,0.9-0.9*0.8*3*textsizeLabelsPP,"EMC");
                SetStyleTLatex( labelInvMassRecoEMC, 0.85*textSizeLabelsPixel,4);
                labelInvMassRecoEMC->SetTextFont(43);
                labelInvMassRecoEMC->Draw();

                if(plotDate){
                  TLatex *labelDate = new TLatex(0.135, 0.9-0.9*4*0.8*textsizeLabelsPP,date.Data());
                  SetStyleTLatex( labelDate, 0.85*textSizeLabelsPixel,4);
                  labelDate->SetTextFont(43);
                  labelDate->Draw();
                }

                SetStyleTLatex( labelInvMassPtRangeEMCAL, 0.85*textSizeLabelsPixel,4);
                labelInvMassPtRangeEMCAL->SetTextAlign(31);
                labelInvMassPtRangeEMCAL->SetTextFont(43);
                labelInvMassPtRangeEMCAL->Draw();

                TLegend* legendInvMassEMCAL  = GetAndSetLegend2(0.67, 0.88-5*0.8*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
                legendInvMassEMCAL->SetMargin(0.25);
                legendInvMassEMCAL->AddEntry(histoPi0InvMassSigPlusBGEMCAL[i],"Raw real events","l");
                legendInvMassEMCAL->AddEntry(histoPi0InvMassBGMixedEMCAL[i],"Mixed event BG","p");
                legendInvMassEMCAL->AddEntry(histoPi0InvMassRemBGEMCAL[i],"Remain. BG","p");
                legendInvMassEMCAL->AddEntry(histoPi0InvMassSigRemBGSubEMCAL[i],"BG subtracted","p");
                legendInvMassEMCAL->AddEntry(fitPi0InvMassSigEMCAL[i], "Fit","l");
                legendInvMassEMCAL->Draw();
                canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBinEMC_%s.%s",outputDir.Data(), nameTrigger[i].Data(), suffix.Data()));
            } else {
                cout << "missing partial input for invariant mass bin for EMC for trigger: " << nameTrigger[i].Data() << endl;
            }

            if (haveAllEtaInvMassEMCAL[i]){
                canvasInvMassSamplePlot->cd();
                histo2DEtaInvMassDummy->GetYaxis()->SetRangeUser(histoEtaInvMassSigRemBGSubEMCAL[i]->GetMinimum(),1.15*histoEtaInvMassSigPlusBGEMCAL[i]->GetMaximum());
                histo2DEtaInvMassDummy->DrawCopy();

                TLatex *labelInvMassPtRangeEMCAL = new TLatex(0.945,0.9,Form("#eta: %s", histoEtaInvMassSigEMCAL[i]->GetTitle()));

                DrawGammaSetMarker(histoEtaInvMassSigPlusBGEMCAL[i], markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
                histoEtaInvMassSigPlusBGEMCAL[i]->SetLineWidth(1);
                histoEtaInvMassSigPlusBGEMCAL[i]->Draw("hist,e,same");
                DrawGammaSetMarker(histoEtaInvMassBGMixedEMCAL[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG1, markerColorInvMassMBG1);
                histoEtaInvMassBGMixedEMCAL[i]->SetLineWidth(markerSizeInvMassMBG);
                histoEtaInvMassBGMixedEMCAL[i]->Draw("same");
                DrawGammaSetMarker(histoEtaInvMassRemBGEMCAL[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG2, markerColorInvMassMBG2);
                histoEtaInvMassRemBGEMCAL[i]->SetLineWidth(markerSizeInvMassMBG);
                histoEtaInvMassRemBGEMCAL[i]->Draw("same");

                DrawGammaSetMarker(histoEtaInvMassSigRemBGSubEMCAL[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
                histoEtaInvMassSigRemBGSubEMCAL[i]->Draw("same");
                fitEtaInvMassSigEMCAL[i]->SetLineWidth(1);
                fitEtaInvMassSigEMCAL[i]->SetRange(0.35,0.695);
                fitEtaInvMassSigEMCAL[i]->SetLineColor(fitColorInvMassSG);
                fitEtaInvMassSigEMCAL[i]->Draw("same");

                TLatex *labelALICE      = new TLatex(0.135,0.9,ALICEperfor.Data());
                SetStyleTLatex( labelALICE, 0.85*textSizeLabelsPixel,4);
                labelALICE->SetTextFont(43);
                labelALICE->Draw();

                TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9-0.9*0.8*textsizeLabelsPP,collisionSystem8TeV.Data());
                SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
                labelInvMassEnergy->SetTextFont(43);
                labelInvMassEnergy->Draw();

                TLatex *labelInvMassTrigger      = new TLatex(0.135,0.9-0.9*0.8*2*textsizeLabelsPP,Form("%s",nameTriggerAlternative[i].Data()));
                SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4);
                labelInvMassTrigger->SetTextFont(43);
                labelInvMassTrigger->Draw();

                TLatex *labelInvMassRecoEMC  = new TLatex(0.135,0.9-0.9*3*0.8*textsizeLabelsPP,"EMC");
                SetStyleTLatex( labelInvMassRecoEMC, 0.85*textSizeLabelsPixel,4);
                labelInvMassRecoEMC->SetTextFont(43);
                labelInvMassRecoEMC->Draw();

                if(plotDate){
                  TLatex *labelDate = new TLatex(0.135, 0.9-0.9*4*0.8*textsizeLabelsPP,date.Data());
                  SetStyleTLatex( labelDate, 0.85*textSizeLabelsPixel,4);
                  labelDate->SetTextFont(43);
                  labelDate->Draw();
                }

                SetStyleTLatex( labelInvMassPtRangeEMCAL, 0.85*textSizeLabelsPixel,4);
                labelInvMassPtRangeEMCAL->SetTextAlign(31);
                labelInvMassPtRangeEMCAL->SetTextFont(43);
                labelInvMassPtRangeEMCAL->Draw();

                TLegend* legendInvMassEMCAL  = GetAndSetLegend2(0.67, 0.88-5*0.8*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
                legendInvMassEMCAL->SetMargin(0.25);
                legendInvMassEMCAL->AddEntry(histoEtaInvMassSigPlusBGEMCAL[i],"Raw real events","l");
                legendInvMassEMCAL->AddEntry(histoEtaInvMassBGMixedEMCAL[i],"Mixed event BG","p");
                legendInvMassEMCAL->AddEntry(histoEtaInvMassRemBGEMCAL[i],"Remain. BG","p");
                legendInvMassEMCAL->AddEntry(histoEtaInvMassSigRemBGSubEMCAL[i],"BG subtracted","p");
                legendInvMassEMCAL->AddEntry(fitEtaInvMassSigEMCAL[i], "Fit","l");
                legendInvMassEMCAL->Draw();
                canvasInvMassSamplePlot->SaveAs(Form("%s/Eta_InvMassBinEMC_%s.%s",outputDir.Data(), nameTrigger[i].Data(), suffix.Data()));
            } else {
                cout << "missing partial input for invariant mass bin for EMC for trigger: " << nameTrigger[i].Data() << endl;
            }
        }
    }

    // **********************************************************************************************************************
    // **************************Plot example invariant mass bin PCM Pi0 ****************************************************
    // **********************************************************************************************************************

    TFile* filePCMMB                                     = new TFile(fileNamePCMMB.Data());
    TDirectory* directoryPCMPi0L                         = (TDirectory*)filePCMMB->Get("Pi08TeV");

    const Int_t nPCMBins = 1;
    Int_t binN[nPCMBins] = {21};
    Double_t PCMbinLow[nPCMBins]={4.0};
    Double_t PCMbinHigh[nPCMBins]={4.5};

    for(Int_t iBin=0; iBin<nPCMBins; iBin++){
      TH1D* histoPCMSignalPlusBGPi0                        = (TH1D*)directoryPCMPi0L->Get(Form("Pi0_InvMassSigPlusBG_%02i_MB",binN[iBin]));
      TH1D* histoPCMMixedBGPi0                             = (TH1D*)directoryPCMPi0L->Get(Form("Pi0_InvMassBG_Example_%02i_MB",binN[iBin]));
      TH1D* histoPCMSignalPi0                              = (TH1D*)directoryPCMPi0L->Get(Form("Pi0_InvMassSig_%02i_MB",binN[iBin]));
      TF1* fitPCMlow                                       = (TF1*)directoryPCMPi0L->Get(Form("Pi0_InvMassSigFit_%02i_MB",binN[iBin]));

      canvasInvMassSamplePlot->cd();
      histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(0.02,0.3);
      histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(1.2*histoPCMSignalPi0->GetMinimum(),1.2*histoPCMSignalPlusBGPi0->GetMaximum());
      histo2DPi0InvMassDummy->DrawCopy();

      histoPCMSignalPi0->Fit(fitPCMlow,"QRME0");
      TF1* fitPCMlowPi0F                                  = new TF1("Linearpp","[0]+[1]*x",0.0,0.3);
      fitPCMlowPi0F->SetParameter(0, fitPCMlow->GetParameter(4));
      fitPCMlowPi0F->SetParameter(1, fitPCMlow->GetParameter(5));
      TVirtualFitter * fitterPCM                             = TVirtualFitter::GetFitter();
      Int_t nFreeParPCM                                      = fitPCMlow->GetNumberFreeParameters();
      double * covMatrixPCM                                  = fitterPCM->GetCovarianceMatrix();

      TH1D* histoPi0InvMassRemBGPCM                             = (TH1D*)histoPCMMixedBGPi0->Clone(Form("Pi0_InvMassRemBG_Example_%i",iBin));
      for (Int_t j = 1; j < histoPi0InvMassRemBGPCM->GetNbinsX()+1; j++){
          histoPi0InvMassRemBGPCM->SetBinContent(j,0);
          histoPi0InvMassRemBGPCM->SetBinError(j,0);
      }
      for (Int_t j = histoPCMSignalPi0->GetXaxis()->FindBin(0.01); j < histoPCMSignalPi0->GetXaxis()->FindBin(0.30)+1; j++){
          Double_t startBinEdge                                   = histoPCMSignalPi0->GetXaxis()->GetBinLowEdge(j);
          Double_t endBinEdge                                     = histoPCMSignalPi0->GetXaxis()->GetBinUpEdge(j);
          Double_t intLinearBack                                  = fitPCMlowPi0F->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
          Double_t errorLinearBck                                 = pow(( pow( (endBinEdge-startBinEdge)*fitPCMlow->GetParError(4),2) +
                                                                          pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPCMlow->GetParError(5),2)
                                                                          +2*covMatrixPCM[nFreeParPCM*nFreeParPCM-2]*(endBinEdge-startBinEdge)*0.5*
                                                                          (endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
//                         cout << j << "\t" << intLinearBack << "\t" << errorLinearBck << endl;
          histoPi0InvMassRemBGPCM->SetBinContent(j,intLinearBack);
          histoPi0InvMassRemBGPCM->SetBinError(j,errorLinearBck);
      }
      //histoPCMMixedBGPi0->Add(histoPi0InvMassRemBGPCM);
      histoPCMSignalPi0->Add(histoPi0InvMassRemBGPCM,-1);
      fitPCMlow->SetParameter(4, 0);
      fitPCMlow->SetParameter(5, 0);

      TLatex *labelInvMassPtRangePCMl = new TLatex(0.945,0.9,Form("#pi^{0}: %.1f GeV/c < p_{T} < %.1f GeV/c",PCMbinLow[iBin],PCMbinHigh[iBin]));

      DrawGammaSetMarker(histoPCMSignalPlusBGPi0, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
      histoPCMSignalPlusBGPi0->SetLineWidth(1);
      histoPCMSignalPlusBGPi0->Draw("hist,e,same");
      DrawGammaSetMarker(histoPCMMixedBGPi0, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG1, markerColorInvMassMBG1);
      histoPCMMixedBGPi0->SetLineWidth(markerSizeInvMassMBG);
      histoPCMMixedBGPi0->Draw("same");
      DrawGammaSetMarker(histoPi0InvMassRemBGPCM, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG2, markerColorInvMassMBG2);
      histoPi0InvMassRemBGPCM->SetLineWidth(markerSizeInvMassMBG);
      histoPi0InvMassRemBGPCM->Draw("same");

      DrawGammaSetMarker(histoPCMSignalPi0, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
      histoPCMSignalPi0->Draw("same");
      fitPCMlow->SetLineWidth(1);
      fitPCMlow->SetRange(0,0.255);
      fitPCMlow->SetLineColor(fitColorInvMassSG);
      fitPCMlow->Draw("same");

      TLatex *labelALICE      = new TLatex(0.135,0.9,ALICEperfor.Data());
      SetStyleTLatex( labelALICE, 0.85*textSizeLabelsPixel,4);
      labelALICE->SetTextFont(43);
      labelALICE->Draw();

      TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9-0.9*0.8*textsizeLabelsPP,collisionSystem8TeV.Data());
      SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
      labelInvMassEnergy->SetTextFont(43);
      labelInvMassEnergy->Draw();

      TLatex *labelInvMassTriggerPCMl      = new TLatex(0.135,0.9-0.9*2*0.8*textsizeLabelsPP,"MB trigger");
      SetStyleTLatex( labelInvMassTriggerPCMl, 0.85*textSizeLabelsPixel,4);
      labelInvMassTriggerPCMl->SetTextFont(43);
      labelInvMassTriggerPCMl->Draw();

      TLatex *labelInvMassRecoPCM  = new TLatex(0.135,0.9-0.9*3*0.8*textsizeLabelsPP,"PCM");
      SetStyleTLatex( labelInvMassRecoPCM, 0.85*textSizeLabelsPixel,4);
      labelInvMassRecoPCM->SetTextFont(43);
      labelInvMassRecoPCM->Draw();

      if(plotDate){
        TLatex *labelDate = new TLatex(0.135, 0.9-0.9*4*0.8*textsizeLabelsPP,date.Data());
        SetStyleTLatex( labelDate, 0.85*textSizeLabelsPixel,4);
        labelDate->SetTextFont(43);
        labelDate->Draw();
      }

      SetStyleTLatex( labelInvMassPtRangePCMl, 0.85*textSizeLabelsPixel,4);
      labelInvMassPtRangePCMl->SetTextAlign(31);
      labelInvMassPtRangePCMl->SetTextFont(43);
      labelInvMassPtRangePCMl->Draw();

      TLegend* legendInvMassPCMl  = GetAndSetLegend2(0.67, 0.88-5*0.8*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
      legendInvMassPCMl->SetMargin(0.25);
      legendInvMassPCMl->AddEntry(histoPCMSignalPlusBGPi0,"Raw real events","l");
      legendInvMassPCMl->AddEntry(histoPCMMixedBGPi0,"Mixed event BG","p");
      legendInvMassPCMl->AddEntry(histoPi0InvMassRemBGPCM,"Remain. BG","p");
      legendInvMassPCMl->AddEntry(histoPCMSignalPi0,"BG subtracted","p");
      legendInvMassPCMl->AddEntry(fitPCMlow, "Fit","l");
      legendInvMassPCMl->Draw();
      canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBinPCM_Pi0_%i.%s",outputDir.Data(), iBin, suffix.Data()));
    }


    // **********************************************************************************************************************
    // **************************Plot example invariant mass bin PCM Eta ****************************************************
    // **********************************************************************************************************************

    TDirectory* directoryPCMEtaL                         = (TDirectory*)filePCMMB->Get("Eta8TeV");

    const Int_t nPCMEtaBins = 2;
    Int_t binNEta[nPCMEtaBins] = {3,11};
    Double_t PCMEtabinLow[nPCMEtaBins]={1.1,5.0};
    Double_t PCMEtabinHigh[nPCMEtaBins]={1.4,6.0};

    for(Int_t iBin=0; iBin<nPCMEtaBins; iBin++){
      TH1D* histoPCMSignalPlusBGEta                        = (TH1D*)directoryPCMEtaL->Get(Form("Eta_InvMassSigPlusBG_%02i_MB",binNEta[iBin]));
      TH1D* histoPCMMixedBGEta                             = (TH1D*)directoryPCMEtaL->Get(Form("Eta_InvMassBG_%02i_MB",binNEta[iBin]));
      TH1D* histoPCMSignalEta                              = (TH1D*)directoryPCMEtaL->Get(Form("Eta_InvMassSig_%02i_MB",binNEta[iBin]));
      TF1* fitPCMlowEta                                    = (TF1*)directoryPCMEtaL->Get(Form("Eta_InvMassSigFit_%02i_MB",binNEta[iBin]));

      canvasInvMassSamplePlot->cd();
      histo2DEtaInvMassDummy->GetXaxis()->SetRangeUser(0.02,0.3);
      histo2DEtaInvMassDummy->GetYaxis()->SetRangeUser(1.2*histoPCMSignalEta->GetMinimum(),1.*histoPCMSignalPlusBGEta->GetMaximum());
      histo2DEtaInvMassDummy->DrawCopy();


        histoPCMSignalEta->Fit(fitPCMlowEta,"QRME0");
        TF1* fitPCMlowEtaF                                  = new TF1("Linearpp","[0]+[1]*x",0.0,1.0);
        fitPCMlowEtaF->SetParameter(0, fitPCMlowEta->GetParameter(4));
        fitPCMlowEtaF->SetParameter(1, fitPCMlowEta->GetParameter(5));
        TVirtualFitter * fitterPCM                             = TVirtualFitter::GetFitter();
        Int_t nFreeParPCM                                      = fitPCMlowEta->GetNumberFreeParameters();
        double * covMatrixPCM                                  = fitterPCM->GetCovarianceMatrix();

        TH1D* histoEtaInvMassRemBGPCM                             = (TH1D*)histoPCMMixedBGEta->Clone(Form("Eta_InvMassRemBG_Example_%i",iBin));
        for (Int_t j = 1; j < histoEtaInvMassRemBGPCM->GetNbinsX()+1; j++){
            histoEtaInvMassRemBGPCM->SetBinContent(j,0);
            histoEtaInvMassRemBGPCM->SetBinError(j,0);
        }
        for (Int_t j = histoPCMSignalEta->GetXaxis()->FindBin(0.30); j < histoPCMSignalEta->GetXaxis()->FindBin(0.70)+1; j++){
            Double_t startBinEdge                                   = histoPCMSignalEta->GetXaxis()->GetBinLowEdge(j);
            Double_t endBinEdge                                     = histoPCMSignalEta->GetXaxis()->GetBinUpEdge(j);
            Double_t intLinearBack                                  = fitPCMlowEtaF->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
            Double_t errorLinearBck                                 = pow(( pow( (endBinEdge-startBinEdge)*fitPCMlowEta->GetParError(4),2) +
                                                                            pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPCMlowEta->GetParError(5),2)
                                                                            +2*covMatrixPCM[nFreeParPCM*nFreeParPCM-2]*(endBinEdge-startBinEdge)*0.5*
                                                                            (endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
  //                         cout << j << "\t" << intLinearBack << "\t" << errorLinearBck << endl;
            histoEtaInvMassRemBGPCM->SetBinContent(j,intLinearBack);
            histoEtaInvMassRemBGPCM->SetBinError(j,errorLinearBck);
        }
        //histoPCMMixedBGEta->Add(histoEtaInvMassRemBGPCM);
        histoPCMSignalEta->Add(histoEtaInvMassRemBGPCM,-1);
        fitPCMlowEta->SetParameter(4, 0);
        fitPCMlowEta->SetParameter(5, 0);


      TLatex *labelInvMassPtRangePCMl = new TLatex(0.945,0.9,Form("#eta: %.1f GeV/c < p_{T} < %.1f GeV/c",PCMEtabinLow[iBin],PCMEtabinHigh[iBin]));

      DrawGammaSetMarker(histoPCMSignalPlusBGEta, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
      histoPCMSignalPlusBGEta->SetLineWidth(1);
      histoPCMSignalPlusBGEta->Draw("hist,e,same");
      DrawGammaSetMarker(histoPCMMixedBGEta, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG1, markerColorInvMassMBG1);
      histoPCMMixedBGEta->SetLineWidth(markerSizeInvMassMBG);
      histoPCMMixedBGEta->Draw("same");
      DrawGammaSetMarker(histoEtaInvMassRemBGPCM, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG2, markerColorInvMassMBG2);
      histoEtaInvMassRemBGPCM->SetLineWidth(markerSizeInvMassMBG);
      histoEtaInvMassRemBGPCM->Draw("same");

      DrawGammaSetMarker(histoPCMSignalEta, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
      histoPCMSignalEta->Draw("same");
      fitPCMlowEta->SetLineWidth(1);
      fitPCMlowEta->SetRange(0.35,0.7);
      fitPCMlowEta->SetLineColor(fitColorInvMassSG);
      fitPCMlowEta->Draw("same");

      TLatex *labelALICE      = new TLatex(0.135,0.9,ALICEperfor.Data());
      SetStyleTLatex( labelALICE, 0.85*textSizeLabelsPixel,4);
      labelALICE->SetTextFont(43);
      labelALICE->Draw();

      TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9-0.9*0.8*textsizeLabelsPP,collisionSystem8TeV.Data());
      SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
      labelInvMassEnergy->SetTextFont(43);
      labelInvMassEnergy->Draw();

      TLatex *labelInvMassTriggerPCMl      = new TLatex(0.135,0.9-0.9*2*0.8*textsizeLabelsPP,"MB trigger");
      SetStyleTLatex( labelInvMassTriggerPCMl, 0.85*textSizeLabelsPixel,4);
      labelInvMassTriggerPCMl->SetTextFont(43);
      labelInvMassTriggerPCMl->Draw();

      TLatex *labelInvMassRecoPCM  = new TLatex(0.135,0.9-0.9*3*0.8*textsizeLabelsPP,"PCM");
      SetStyleTLatex( labelInvMassRecoPCM, 0.85*textSizeLabelsPixel,4);
      labelInvMassRecoPCM->SetTextFont(43);
      labelInvMassRecoPCM->Draw();

      if(plotDate){
        TLatex *labelDate = new TLatex(0.135, 0.9-0.9*4*0.8*textsizeLabelsPP,date.Data());
        SetStyleTLatex( labelDate, 0.85*textSizeLabelsPixel,4);
        labelDate->SetTextFont(43);
        labelDate->Draw();
      }

      SetStyleTLatex( labelInvMassPtRangePCMl, 0.85*textSizeLabelsPixel,4);
      labelInvMassPtRangePCMl->SetTextAlign(31);
      labelInvMassPtRangePCMl->SetTextFont(43);
      labelInvMassPtRangePCMl->Draw();

      TLegend* legendInvMassPCMl  = GetAndSetLegend2(0.67, 0.88-5*0.8*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
      legendInvMassPCMl->SetMargin(0.25);
      legendInvMassPCMl->AddEntry(histoPCMSignalPlusBGEta,"Raw real events","l");
      legendInvMassPCMl->AddEntry(histoPCMMixedBGEta,"Mixed event BG","p");
      legendInvMassPCMl->AddEntry(histoEtaInvMassRemBGPCM,"Remain. BG","p");
      legendInvMassPCMl->AddEntry(histoPCMSignalEta,"BG subtracted","p");
      legendInvMassPCMl->AddEntry(fitPCMlowEta, "Fit","l");
      legendInvMassPCMl->Draw();
      canvasInvMassSamplePlot->SaveAs(Form("%s/Eta_InvMassBinPCM_Eta_%i.%s",outputDir.Data(), iBin, suffix.Data()));
    }

      // **********************************************************************************************************************
      // **************************Plot example invariant mass bin PHOS low pt ***********************************************
      // **********************************************************************************************************************

      TFile* filePHOSMB                                     = new TFile(fileNamePHOSMB.Data());
      TDirectory* directoryPHOSPi0L                         = (TDirectory*)filePHOSMB->Get("Pi08TeV");

//      const Int_t nBins = 28;
//      Double_t binLow[nBins]={1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,10.0,11.0};
//      Double_t binHigh[nBins]={1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,10.0,11.0,12.0};

      const Int_t nBins = 2;
      Double_t binLow[nBins]={1.6,5.5};
      Double_t binHigh[nBins]={1.8,6.0};

      for(Int_t iBin=0; iBin<nBins; iBin++){
        TH1D* histoPHOSSignalPlusBGPi0                        = (TH1D*)directoryPHOSPi0L->Get(Form("InvMassSigPlusBG_PtBin_%.1fto%.1f",binLow[iBin],binHigh[iBin]));
        TH1D* histoPHOSCombBGPi0                              = (TH1D*)directoryPHOSPi0L->Get(Form("InvMassCombBG_PtBin_%.1fto%.1f",binLow[iBin],binHigh[iBin]));
        TH1D* histoPHOSResiBGPi0                              = (TH1D*)directoryPHOSPi0L->Get(Form("InvMassResiBG_PtBin_%.1fto%.1f",binLow[iBin],binHigh[iBin]));
        TH1D* histoPHOSSignalPi0                              = (TH1D*)directoryPHOSPi0L->Get(Form("InvMassSig_PtBin_%.1fto%.1f",binLow[iBin],binHigh[iBin]));
        TF1* fitPHOSlow                                       = (TF1*)directoryPHOSPi0L->Get(Form("InvMassFinalFit_%.1fto%.1f",binLow[iBin],binHigh[iBin]));

        canvasInvMassSamplePlot->cd();
        histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(0.02,0.3);
        histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(1.2*histoPHOSSignalPi0->GetMinimum(),1.2*histoPHOSSignalPlusBGPi0->GetMaximum());
        histo2DPi0InvMassDummy->DrawCopy();

        TLatex *labelInvMassPtRangePHOSl = new TLatex(0.945,0.9,Form("#pi^{0}: %.1f GeV/c < p_{T} < %.1f GeV/c",binLow[iBin],binHigh[iBin]));

        DrawGammaSetMarker(histoPHOSSignalPlusBGPi0, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
        histoPHOSSignalPlusBGPi0->SetLineWidth(1);
        histoPHOSSignalPlusBGPi0->Draw("hist,e,same");
        DrawGammaSetMarker(histoPHOSCombBGPi0, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG1, markerColorInvMassMBG1);
        histoPHOSCombBGPi0->SetLineWidth(markerSizeInvMassMBG);
        histoPHOSCombBGPi0->Draw("same");
        DrawGammaSetMarker(histoPHOSResiBGPi0, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG2, markerColorInvMassMBG2);
        histoPHOSResiBGPi0->SetLineWidth(markerSizeInvMassMBG);
        histoPHOSResiBGPi0->Draw("same");

        DrawGammaSetMarker(histoPHOSSignalPi0, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
        histoPHOSSignalPi0->SetStats(kFALSE);
        histoPHOSSignalPi0->Draw("same");
        fitPHOSlow->SetLineWidth(1);
        fitPHOSlow->SetRange(0,0.255);
        fitPHOSlow->SetLineColor(fitColorInvMassSG);
        fitPHOSlow->Draw("same");

        TLatex *labelALICE      = new TLatex(0.135,0.9,ALICEperfor.Data());
        SetStyleTLatex( labelALICE, 0.85*textSizeLabelsPixel,4);
        labelALICE->SetTextFont(43);
        labelALICE->Draw();

        TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9-0.9*0.8*textsizeLabelsPP,collisionSystem8TeV.Data());
        SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
        labelInvMassEnergy->SetTextFont(43);
        labelInvMassEnergy->Draw();

        TLatex *labelInvMassTriggerPHOSl      = new TLatex(0.135,0.9-0.9*2*0.8*textsizeLabelsPP,"MB trigger");
        SetStyleTLatex( labelInvMassTriggerPHOSl, 0.85*textSizeLabelsPixel,4);
        labelInvMassTriggerPHOSl->SetTextFont(43);
        labelInvMassTriggerPHOSl->Draw();

        TLatex *labelInvMassRecoPHOS  = new TLatex(0.135,0.9-0.9*3*0.8*textsizeLabelsPP,"PHOS");
        SetStyleTLatex( labelInvMassRecoPHOS, 0.85*textSizeLabelsPixel,4);
        labelInvMassRecoPHOS->SetTextFont(43);
        labelInvMassRecoPHOS->Draw();

        if(plotDate){
          TLatex *labelDate = new TLatex(0.135, 0.9-0.9*4*0.8*textsizeLabelsPP,date.Data());
          SetStyleTLatex( labelDate, 0.85*textSizeLabelsPixel,4);
          labelDate->SetTextFont(43);
          labelDate->Draw();
        }

        SetStyleTLatex( labelInvMassPtRangePHOSl, 0.85*textSizeLabelsPixel,4);
        labelInvMassPtRangePHOSl->SetTextAlign(31);
        labelInvMassPtRangePHOSl->SetTextFont(43);
        labelInvMassPtRangePHOSl->Draw();

        TLegend* legendInvMassPHOSl  = GetAndSetLegend2(0.67, 0.88-5*0.8*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
        legendInvMassPHOSl->SetMargin(0.25);
        legendInvMassPHOSl->AddEntry(histoPHOSSignalPlusBGPi0,"Raw real events","l");
        legendInvMassPHOSl->AddEntry(histoPHOSCombBGPi0,"Mixed event BG","p");
        legendInvMassPHOSl->AddEntry(histoPHOSResiBGPi0,"Remain. BG","p");
        legendInvMassPHOSl->AddEntry(histoPHOSSignalPi0,"BG subtracted","p");
        legendInvMassPHOSl->AddEntry(fitPHOSlow, "Fit","l");
        legendInvMassPHOSl->Draw();
        canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBinPHOS_low_%i.%s",outputDir.Data(), iBin, suffix.Data()));
      }

      // **********************************************************************************************************************
      // **************************Plot example invariant mass bin PHOS high pt ***********************************************
      // **********************************************************************************************************************

      //TFile* filePHOSPHOS                                     = new TFile(fileNamePHOSPHOS.Data());
      TDirectory* directoryPHOSPi0H                           = (TDirectory*)filePHOSPHOS->Get("Pi08TeV");

//      const Int_t nBinsTr = 21;
//      Double_t binTrLow[nBinsTr]={5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,20.0,22.0,26.0,30.0};
//      Double_t binTrHigh[nBinsTr]={6.0,6.5,7.0,7.5,8.0,8.5,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,20.0,22.0,26.0,30.0,35.0};

      const Int_t nBinsTr = 2;
      Double_t binTrLow[nBinsTr]={12.0,22.0};
      Double_t binTrHigh[nBinsTr]={13.0,26.0};

      for(Int_t iBin=0; iBin<nBinsTr; iBin++){
        TH1D* histoPHOSHighSignalPlusBGPi0                      = (TH1D*)directoryPHOSPi0H->Get(Form("InvMassSigPlusBG_PtBin_%.1fto%.1f",binTrLow[iBin],binTrHigh[iBin]));
        TH1D* histoPHOSHighCombBGPi0                            = (TH1D*)directoryPHOSPi0H->Get(Form("InvMassCombBG_PtBin_%.1fto%.1f",binTrLow[iBin],binTrHigh[iBin]));
        TH1D* histoPHOSHighResiBGPi0                            = (TH1D*)directoryPHOSPi0H->Get(Form("InvMassResiBG_PtBin_%.1fto%.1f",binTrLow[iBin],binTrHigh[iBin]));
        TH1D* histoPHOSHighSignalPi0                            = (TH1D*)directoryPHOSPi0H->Get(Form("InvMassSig_PtBin_%.1fto%.1f",binTrLow[iBin],binTrHigh[iBin]));
        TF1* fitPHOShigh                                        = (TF1*)directoryPHOSPi0H->Get(Form("InvMassFinalFit_%.1fto%.1f",binTrLow[iBin],binTrHigh[iBin]));

        canvasInvMassSamplePlot->cd();
        histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(1.2*histoPHOSHighSignalPi0->GetMinimum(),1.3*histoPHOSHighSignalPlusBGPi0->GetMaximum());
        histo2DPi0InvMassDummy->DrawCopy();

        TLatex *labelInvMassPtRangePHOSh = new TLatex(0.945,0.9,Form("#pi^{0}: %.1f GeV/c < p_{T} < %.1f GeV/c",binTrLow[iBin],binTrHigh[iBin]));

        DrawGammaSetMarker(histoPHOSHighSignalPlusBGPi0, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
        histoPHOSHighSignalPlusBGPi0->SetLineWidth(1);
        histoPHOSHighSignalPlusBGPi0->Draw("hist,e,same");
        DrawGammaSetMarker(histoPHOSHighCombBGPi0, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG1, markerColorInvMassMBG1);
        histoPHOSHighCombBGPi0->SetLineWidth(markerSizeInvMassMBG);
        histoPHOSHighCombBGPi0->Draw("same");
        DrawGammaSetMarker(histoPHOSHighResiBGPi0, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG2, markerColorInvMassMBG2);
        histoPHOSHighResiBGPi0->SetLineWidth(markerSizeInvMassMBG);
        histoPHOSHighResiBGPi0->Draw("same");

        DrawGammaSetMarker(histoPHOSHighSignalPi0, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
        histoPHOSHighSignalPi0->SetStats(kFALSE);
        histoPHOSHighSignalPi0->Draw("same");
        fitPHOShigh->SetLineWidth(1);
        fitPHOShigh->SetRange(0,0.255);
        fitPHOShigh->SetLineColor(fitColorInvMassSG);
        fitPHOShigh->Draw("same");

        TLatex *labelALICE      = new TLatex(0.135,0.9,ALICEperfor.Data());
        SetStyleTLatex( labelALICE, 0.85*textSizeLabelsPixel,4);
        labelALICE->SetTextFont(43);
        labelALICE->Draw();

        TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9-0.9*0.8*textsizeLabelsPP,collisionSystem8TeV.Data());
        SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
        labelInvMassEnergy->SetTextFont(43);
        labelInvMassEnergy->Draw();

        TLatex *labelInvMassTriggerPHOSh      = new TLatex(0.135,0.9-0.9*2*0.8*textsizeLabelsPP,"PHOS-L0 trigger");
        SetStyleTLatex( labelInvMassTriggerPHOSh, 0.85*textSizeLabelsPixel,4);
        labelInvMassTriggerPHOSh->SetTextFont(43);
        labelInvMassTriggerPHOSh->Draw();

        TLatex *labelInvMassRecoPHOS  = new TLatex(0.135,0.9-0.9*3*0.8*textsizeLabelsPP,"PHOS");
        SetStyleTLatex( labelInvMassRecoPHOS, 0.85*textSizeLabelsPixel,4);
        labelInvMassRecoPHOS->SetTextFont(43);
        labelInvMassRecoPHOS->Draw();

        if(plotDate){
          TLatex *labelDate = new TLatex(0.135, 0.9-0.9*4*0.8*textsizeLabelsPP,date.Data());
          SetStyleTLatex( labelDate, 0.85*textSizeLabelsPixel,4);
          labelDate->SetTextFont(43);
          labelDate->Draw();
        }

        SetStyleTLatex( labelInvMassPtRangePHOSh, 0.85*textSizeLabelsPixel,4);
        labelInvMassPtRangePHOSh->SetTextAlign(31);
        labelInvMassPtRangePHOSh->SetTextFont(43);
        labelInvMassPtRangePHOSh->Draw();

        TLegend* legendInvMassPHOSh  = GetAndSetLegend2(0.67, 0.88-5*0.8*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
        legendInvMassPHOSh->SetMargin(0.25);
        legendInvMassPHOSh->AddEntry(histoPHOSHighSignalPlusBGPi0,"Raw real events","l");
        legendInvMassPHOSh->AddEntry(histoPHOSHighCombBGPi0,"Mixed event BG","p");
        legendInvMassPHOSh->AddEntry(histoPHOSHighResiBGPi0,"Remain. BG","p");
        legendInvMassPHOSh->AddEntry(histoPHOSHighSignalPi0,"BG subtracted","p");
        legendInvMassPHOSh->AddEntry(fitPHOShigh, "Fit","l");
        legendInvMassPHOSh->Draw();
        canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBinPHOS_high_%i.%s",outputDir.Data(), iBin, suffix.Data()));
      }
 // **********************************************************************************************************************
 // ************************* Saving of final results ********************************************************************
 // **********************************************************************************************************************

    TString nameOutputCommonFile    = Form("CombinedResultsPaperPP8TeV_%s.root", dateForOutput.Data());
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
//        graphEMCALMergedPi0InvXSectionStat->Write("graphInvCrossSectionPi0EMCALMerged8TeVStatErr");
//        graphEMCALMergedPi0InvXSectionSys->Write("graphInvCrossSectionPi0EMCALMerged8TeVSysErr");
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
//            if (graphEMCALMergedPi0InvXSectionStat_yShifted) graphEMCALMergedPi0InvXSectionStat_yShifted->Write("graphInvCrossSectionPi0EMCALMerged8TeVStatErr_yShifted");
//            if (graphEMCALMergedPi0InvXSectionSys_yShifted) graphEMCALMergedPi0InvXSectionSys_yShifted->Write("graphInvCrossSectionPi0EMCALMerged8TeVSysErr_yShifted");
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
//            graphEMCALMergedPi0AccTimesEffDivPur->Write("Pi0CorrectionFactorEMCALMerged");

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

    fCombResults.mkdir("Eta8TeV");
    TDirectoryFile* directoryEta = (TDirectoryFile*)fCombResults.Get("Eta8TeV");
    fCombResults.cd("Eta8TeV");
        // PCM component
        graphPCMEtaInvXSectionStat->Write("graphInvCrossSectionEtaPCM8TeVStatErr");
        graphPCMEtaInvXSectionSys->Write("graphInvCrossSectionEtaPCM8TeVSysErr");
        // EMCAL component
        graphEMCALEtaInvXSectionStat->Write("graphInvCrossSectionEtaEMCAL8TeVStatErr");
        graphEMCALEtaInvXSectionSys->Write("graphInvCrossSectionEtaEMCAL8TeVSysErr");
        // PCM-EMCal component
        graphPCMEMCALEtaInvXSectionStat->Write("graphInvCrossSectionEtaPCMEMCAL8TeVStatErr");
        graphPCMEMCALEtaInvXSectionSys->Write("graphInvCrossSectionEtaPCMEMCAL8TeVSysErr");
        // Final spectrum correlations Method A
        graphCombEtaInvXSectionTotA->Write("graphInvCrossSectionEtaComb8TeVA");
        graphCombEtaInvXSectionStatA->Write("graphInvCrossSectionEtaComb8TeVAStatErr");
        graphCombEtaInvXSectionSysA->Write("graphInvCrossSectionEtaComb8TeVASysErr");

        // fits for eta
        fitInvXSectionEta->Write("TsallisFitEta");
        fitTCMInvXSectionEta->Write("TwoComponentModelFitEta");

        // writing Y shifted graphs in addition
        if (bWCorrection.Contains("Y")){
            if(graphPCMEtaInvXSectionStat_yShifted)graphPCMEtaInvXSectionStat_yShifted->Write("graphInvCrossSectionEtaPCM8TeVStatErr_yShifted");
            if(graphPCMEtaInvXSectionSys_yShifted)graphPCMEtaInvXSectionSys_yShifted->Write("graphInvCrossSectionEtaPCM8TeVSysErr_yShifted");
            // EMCAL component
            if(graphEMCALEtaInvXSectionStat_yShifted)graphEMCALEtaInvXSectionStat_yShifted->Write("graphInvCrossSectionEtaEMCAL8TeVStatErr_yShifted");
            if(graphEMCALEtaInvXSectionSys_yShifted)graphEMCALEtaInvXSectionSys_yShifted->Write("graphInvCrossSectionEtaEMCAL8TeVSysErr_yShifted");
            // PCM-EMCal component
            if(graphPCMEMCALEtaInvXSectionStat_yShifted)graphPCMEMCALEtaInvXSectionStat_yShifted->Write("graphInvCrossSectionEtaPCMEMCAL8TeVStatErr_yShifted");
            if(graphPCMEMCALEtaInvXSectionSys_yShifted)graphPCMEMCALEtaInvXSectionSys_yShifted->Write("graphInvCrossSectionEtaPCMEMCAL8TeVSysErr_yShifted");
            // Final spectrum correlations Method A
            if(graphCombEtaInvXSectionTotA_yShifted)graphCombEtaInvXSectionTotA_yShifted->Write("graphInvCrossSectionEtaComb8TeVA_yShifted");
            if(graphCombEtaInvXSectionStatA_yShifted)graphCombEtaInvXSectionStatA_yShifted->Write("graphInvCrossSectionEtaComb8TeVAStatErr_yShifted");
            if(graphCombEtaInvXSectionSysA_yShifted)graphCombEtaInvXSectionSysA_yShifted->Write("graphInvCrossSectionEtaComb8TeVASysErr_yShifted");
        }

        histoPCMEtaToPi0Stat->Write("histoRatioEtaToPi0PCM8TeVStatErr");
        graphPCMEtaToPi0Sys->Write("graphRatioEtaToPi0PCM8TeVSysErr");
        histoEMCALEtaToPi0Stat->Write("histoRatioEtaToPi0EMCAL8TeVStatErr");
        graphEMCALEtaToPi0Sys->Write("graphRatioEtaToPi0EMCAL8TeVSysErr");
        histoPCMEMCALEtaToPi0Stat->Write("histoRatioEtaToPi0PCMEMCAL8TeVStatErr");
        graphPCMEMCALEtaToPi0Sys->Write("graphRatioEtaToPi0PCMEMCAL8TeVSysErr");
        graphCombEtaToPi0TotA->Write("graphRatioEtaToPi0Comb8TeVTotErr");
        graphCombEtaToPi0StatA->Write("graphRatioEtaToPi0Comb8TeVStatErr");
        graphCombEtaToPi0SysA->Write("graphRatioEtaToPi0Comb8TeVSysErr");

        graphRatioForMt_stat->Write("graphRatioMtScalingToEtaToPi0StatErr");
        graphRatioForMt_sys->Write("graphRatioMtScalingToEtaToPi0SysErr");

        directoryEta->mkdir("Supporting");
        directoryEta->cd("Supporting");
            // Writing full correction factors
            histoPCMEtaAccTimesEff->Write("EtaCorrectionFactorPCM");
            graphEMCALEtaAccTimesEff->Write("EtaCorrectionFactorEMCAL");
            graphPCMEMCALEtaAccTimesEff->Write("EtaCorrectionFactorPCMEMCAL");

            histoPCMEtaMass->Write("EtaMassDataPCM");
            histoPCMEtaTrueMass->Write("EtaMassMCPCM");
            graphEMCALEtaMass->Write("EtaMassDataEMCAL");
            graphEMCALEtaMassMC->Write("EtaMassMCEMCAL");
            graphPCMEMCALEtaMass->Write("EtaMassDataPCMEMCAL");
            graphPCMEMCALEtaMassMC->Write("EtaMassMCPCMEMCAL");

            histoPCMEtaFWHMMeV->Write("EtaWidthDataPCM");
            histoPCMEtaTrueFWHMMeV->Write("EtaWidthMCPCM");
            graphEMCALEtaFWHM->Write("EtaWidthDataEMCAL");
            graphEMCALEtaFWHMMC->Write("EtaWidthMCEMCAL");
            graphPCMEMCALEtaFWHM->Write("EtaWidthDataPCMEMCAL");
            graphPCMEMCALEtaFWHMMC->Write("EtaWidthMCPCMEMCAL");

    fCombResults.Close();


 // **********************************************************************************************************************
 // ************************* Saving only fits to final results **********************************************************
 // **********************************************************************************************************************

    TString nameOutputCommonFileFitsOnly    = Form("FitsPaperPP8TeV_%s.root", dateForOutput.Data());
    TFile fFitsResults(nameOutputCommonFileFitsOnly.Data(), "RECREATE");

         // fits for pi0
        fitInvXSectionPi0->Write("TsallisFitPi0");
        fitTCMInvXSectionPi0Plot->Write("TwoComponentModelFitPi0");

        // fits for eta
        fitInvXSectionEta->Write("TsallisFitEta");
        fitTCMInvXSectionEta->Write("TwoComponentModelFitEta");

    fFitsResults.Close();

}
