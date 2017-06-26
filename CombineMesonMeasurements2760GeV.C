/****************************************************************************************************************************
******         provided by Gamma Conversion Group, PWG4,                                                     *****
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
extern TBenchmark*    gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*      gMinuit;

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    // TString name;
};

void CombineMesonMeasurements2760GeV(   TString fileNamePCM         = "", 
                                        TString fileNamePCMEMCAL    = "", 
                                        TString fileNameEMCALLow    = "",  
                                        TString fileNameEMCALmerged = "",  
                                        TString suffix              = "eps", 
                                        TString isMC                = "", 
                                        TString thesisPlots         = "", 
                                        TString bWCorrection        = "X",
                                        Int_t  flagMerged           = 2,
                                        TString fileNameCorrFactors = "",
                                        Bool_t plotInvMassBins      = kFALSE,
                                        Bool_t flagPCMfile          = 0
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
    TString collisionSystem2760GeV              = "pp, #sqrt{#it{s}} = 2.76 TeV";
    
    TString fileNameTheory                      = "ExternalInput/Theory/TheoryCompilationPP.root";
    TString fileNameEtaToPi0WorldData           = "ExternalInput/WorldDataPi0Eta.root";
    TString fileNamePHOS                        = "ExternalInput/PHOS/2.76TeV/LHC11a_PHOS_pi0_pp2760_noBWCorr_FDcorr_20140218.root";
    TString fileNamePCMEta                      = "FinalResults/data_PCMResultsFullCorrection_PP_NoBinShifting_usedTosvnEtapaper.root";
    TString fileNameChargedPionPP               = "ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_20_May_2015.root";
    TString fileNameChargedHadronPP             = "ExternalInput/UnidentifiedCharged/ChargedHadrinSpectraPP_2016_07_04.root";
    TString fileNameOldPi0Publication           = "FinalResults/CombinedResultsPP_ShiftedX_PaperRAA_16_May_2014.root";
    TString outputDir                           = Form("%s/%s/CombineMesonMeasurements2760GeV%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
    if (flagMerged == 0){
      outputDir = outputDir+"_HaitaoMerged";  
    } else if (flagMerged == 1){
      outputDir = outputDir+"_FrediMergedV1";  
    } else if (flagMerged == 2){
      outputDir = outputDir+"_FrediMergedV2";  
    } else if (flagMerged == 3){
      outputDir = outputDir+"_FrediMergedV1NLM1";  
    } else if (flagMerged == 4){
      outputDir = outputDir+"_FrediMergedV1NLM2";  
    }  

    cout << outputDir.Data() << endl;
    cout << fileNamePCM.Data() << endl;    
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMEta.root", fileNamePCMEta.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMEMCAL.root", fileNamePCMEMCAL.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCALLow.root", fileNameEMCALLow.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCALmerged.root", fileNameEMCALmerged.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/Theory.root", fileNameTheory.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/ChargedPionsPP.root", fileNameChargedPionPP.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/ChargedHadronsPP.root", fileNameChargedHadronPP.Data(), outputDir.Data()));
    

    TString nameFinalResDat                     = Form("%s/CombinedResults%s_FitResults.dat",outputDir.Data(),bWCorrection.Data());
    fstream  fileFitsOutput;
    fileFitsOutput.open(nameFinalResDat.Data(), ios::out);  


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
//     TString  nameMeasGlobalLabel[11]            = { "PCM^{2}", "PHOS^{2}", "EMC^{2}", "PCM-PHOS", "PCM-EMC",
    TString  nameMeasGlobalLabel[11]            = { "PCM", "PHOS", "EMC", "PCM-PHOS", "PCM-EMC",
                                                    "PCM-Dalitz", "PHOS-Dalitz", "EMC-Dalitz", "EMC high pT", "mEMC",
                                                    "PCMOtherDataset"};
    TString  nameTrigger[6]                     = {"INT1", "INT7", "EMC1", "EMC7", "EG2", "EG1"};
    TString  nameSecPi0SourceRead[4]            = {"K0S", "K0L", "Lambda", "Rest"};
    TString  nameSecPi0SourceLabel[4]           = {"K^{0}_{s}", "K^{0}_{l}", "#Lambda", "had. int."};
    Double_t maxSecCorr[4]                      = { 0.15, 0.0020, 0.00035, 0.045};

    
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
    Style_t  styleLineCGC                       = 2;
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
    Double_t rangePtExMin[11]                   = { 1.6,    0.0,    2.6,    0.0,    1.6,
                                                    0.0,    0.0,    0.0,    9.0,    0.0,
                                                    0.0 };
                                            
    Double_t rangePtExMax[11]                   = { 1.8,    0.0,    3.0,    0.0,    1.8,
                                                    0.0,    0.0,    0.0,    10.,    0.0, 
                                                    0.0 };
    
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
    TFile* filePCM                                  = NULL;
    TDirectory* directoryPCMPi0                     = NULL;
    TGraphAsymmErrors* graphPCMPi0Mass              = NULL;
    TGraphAsymmErrors* graphPCMPi0FWHM              = NULL;
    TGraphAsymmErrors* graphPCMPi0MassMC            = NULL;
    TGraphAsymmErrors* graphPCMPi0FWHMMC            = NULL;
    TGraphAsymmErrors* graphPCMPi0Acc               = NULL;
    TGraphAsymmErrors* graphPCMPi0EffPt             = NULL;
    TGraphAsymmErrors* graphPCMPi0InvXSectionStat   = NULL;
    TH1D* histoPCMPi0InvXSectionStat                = NULL;
    TGraphAsymmErrors* graphPCMPi0InvXSectionSys    = NULL;
    TGraphAsymmErrors* graphPCMPi0AccTimesEff       = NULL;
    TDirectory* directoryPCMEta                     = NULL;
    TGraphAsymmErrors* graphPCMEtaMass              = NULL;
    TGraphAsymmErrors* graphPCMEtaFWHM              = NULL;
    TGraphAsymmErrors* graphPCMEtaMassMC            = NULL;
    TGraphAsymmErrors* graphPCMEtaFWHMMC            = NULL;
    TGraphAsymmErrors* graphPCMEtaAcc               = NULL;
    TGraphAsymmErrors* graphPCMEtaEffPt             = NULL;
    TGraphAsymmErrors* graphPCMEtaInvXSectionStat   = NULL;
    TH1D* histoPCMEtaInvXSectionStat                = NULL;
    TGraphAsymmErrors* graphPCMEtaInvXSectionSys    = NULL;
    TGraphAsymmErrors* graphPCMEtaAccTimesEff       = NULL;
    TH1D* histoPCMEtaToPi0Stat                      = NULL;
    TGraphAsymmErrors* graphPCMEtaToPi0Sys          = NULL;

    if (!flagPCMfile){  
        filePCM                                     = new TFile(fileNamePCM.Data());
        directoryPCMPi0                             = (TDirectory*)filePCM->Get("Pi02.76TeV"); 
            TH1D* histoPCMPi0Mass                       = (TH1D*)directoryPCMPi0->Get("MassPi0");
            TH1D* histoPCMPi0FWHMMeV                    = (TH1D*)directoryPCMPi0->Get("FWHMPi0MeV");
            TH1D* histoPCMPi0TrueMass                   = (TH1D*)directoryPCMPi0->Get("TrueMassPi0");
            TH1D* histoPCMPi0TrueFWHMMeV                = (TH1D*)directoryPCMPi0->Get("TrueFWHMPi0MeV");
            TH1D* histoPCMPi0Acc                        = (TH1D*)directoryPCMPi0->Get("AcceptancePi0");
            TH1D* histoPCMPi0TrueEffPt                  = (TH1D*)directoryPCMPi0->Get("EfficiencyPi0");
            histoPCMPi0InvXSectionStat                  = (TH1D*)directoryPCMPi0->Get("InvCrossSectionPi0");
            graphPCMPi0InvXSectionStat                  = new TGraphAsymmErrors(histoPCMPi0InvXSectionStat);
            graphPCMPi0InvXSectionStat->RemovePoint(graphPCMPi0InvXSectionStat->GetN()-1);
            graphPCMPi0InvXSectionStat->RemovePoint(0);
            graphPCMPi0InvXSectionSys                   = (TGraphAsymmErrors*)directoryPCMPi0->Get("InvCrossSectionPi0Sys");
            graphPCMPi0InvXSectionSys->RemovePoint(graphPCMPi0InvXSectionSys->GetN()-1);
            TH1D* histoPCMPi0AccTimesEff                = (TH1D*)histoPCMPi0TrueEffPt->Clone("histoPCMPi0AccTimesEff");
            histoPCMPi0AccTimesEff->Multiply(histoPCMPi0Acc);
            // normalize to full acceptance (delta y and phi)
            histoPCMPi0AccTimesEff->Scale(2*TMath::Pi()*1.6);
            histoPCMPi0Mass->Scale(1000.);
            histoPCMPi0TrueMass->Scale(1000.);

            graphPCMPi0Mass                             = new TGraphAsymmErrors(histoPCMPi0Mass);
            graphPCMPi0FWHM                             = new TGraphAsymmErrors(histoPCMPi0FWHMMeV);
            graphPCMPi0MassMC                           = new TGraphAsymmErrors(histoPCMPi0TrueMass);
            graphPCMPi0FWHMMC                           = new TGraphAsymmErrors(histoPCMPi0TrueFWHMMeV);
            graphPCMPi0Acc                              = new TGraphAsymmErrors(histoPCMPi0Acc);
            graphPCMPi0EffPt                            = new TGraphAsymmErrors(histoPCMPi0TrueEffPt);
            graphPCMPi0AccTimesEff                      = new TGraphAsymmErrors(histoPCMPi0AccTimesEff);
            
            while(graphPCMPi0AccTimesEff->GetY()[0] == 0){
                graphPCMPi0Mass->RemovePoint(0);
                graphPCMPi0FWHM->RemovePoint(0);
                graphPCMPi0MassMC->RemovePoint(0);
                graphPCMPi0FWHMMC->RemovePoint(0);
                graphPCMPi0Acc->RemovePoint(0);
                graphPCMPi0EffPt->RemovePoint(0);
                graphPCMPi0AccTimesEff->RemovePoint(0);
            }
                
        TFile* filePCMEta                           = new TFile(fileNamePCMEta.Data());
        directoryPCMEta                             = (TDirectory*)filePCMEta->Get("Eta2.76TeV"); 
            TH1D* histoPCMEtaMass                       = (TH1D*)directoryPCMEta->Get("MassEta");
            TH1D* histoPCMEtaFWHMMeV                    = (TH1D*)directoryPCMEta->Get("FWHMEtaMeV");
            TH1D* histoPCMEtaTrueMass                   = (TH1D*)directoryPCMEta->Get("TrueMassEta");
            TH1D* histoPCMEtaTrueFWHMMeV                = (TH1D*)directoryPCMEta->Get("TrueFWHMEtaMeV");
            TH1D* histoPCMEtaAcc                        = (TH1D*)directoryPCMEta->Get("AcceptanceEta");
            TH1D* histoPCMEtaTrueEffPt                  = (TH1D*)directoryPCMEta->Get("EfficiencyEta");
            histoPCMEtaInvXSectionStat                  = (TH1D*)directoryPCMEta->Get("InvCrossSectionEta");
            graphPCMEtaInvXSectionStat                  = new TGraphAsymmErrors(histoPCMEtaInvXSectionStat);
            graphPCMEtaInvXSectionStat->RemovePoint(0);

            TH1D* histoPCMEtaAccTimesEff                = (TH1D*)histoPCMEtaTrueEffPt->Clone("histoPCMEtaAccTimesEff");
            histoPCMEtaAccTimesEff->Multiply(histoPCMEtaAcc);
            // normalize to full acceptance (delta y and phi)
            histoPCMEtaAccTimesEff->Scale(2*TMath::Pi()*1.6);
            histoPCMEtaMass->Scale(1000.);
            histoPCMEtaTrueMass->Scale(1000.);
        
            graphPCMEtaMass                             = new TGraphAsymmErrors(histoPCMEtaMass);
            graphPCMEtaFWHM                             = new TGraphAsymmErrors(histoPCMEtaFWHMMeV);
            graphPCMEtaMassMC                           = new TGraphAsymmErrors(histoPCMEtaTrueMass);
            graphPCMEtaFWHMMC                           = new TGraphAsymmErrors(histoPCMEtaTrueFWHMMeV);
            graphPCMEtaAcc                              = new TGraphAsymmErrors(histoPCMEtaAcc);
            graphPCMEtaEffPt                            = new TGraphAsymmErrors(histoPCMEtaTrueEffPt);
            graphPCMEtaAccTimesEff                      = new TGraphAsymmErrors(histoPCMEtaAccTimesEff);

            while(graphPCMEtaAccTimesEff->GetY()[0] == 0){
                graphPCMEtaMass->RemovePoint(0);
                graphPCMEtaFWHM->RemovePoint(0);
                graphPCMEtaMassMC->RemovePoint(0);
                graphPCMEtaFWHMMC->RemovePoint(0);
                graphPCMEtaAcc->RemovePoint(0);
                graphPCMEtaEffPt->RemovePoint(0);
                graphPCMEtaAccTimesEff->RemovePoint(0);
            }
            
            
            graphPCMEtaInvXSectionSys                   = (TGraphAsymmErrors*)directoryPCMEta->Get("InvCrossSectionEtaSys");
            histoPCMEtaToPi0Stat                        = (TH1D*)directoryPCMEta->Get("EtatoPi0RatioConversionBinShifted");
            graphPCMEtaToPi0Sys                         = (TGraphAsymmErrors*)directoryPCMEta->Get("EtatoPi0RatioConversionBinShiftedSys");
    } else {
        filePCM                                     = new TFile(fileNamePCM.Data());
        directoryPCMPi0                             = (TDirectory*)filePCM->Get("Pi02.76TeV"); 
            graphPCMPi0Mass                             = (TGraphAsymmErrors*)directoryPCMPi0->Get("Pi0_Mass_data");
            graphPCMPi0Mass                             = ScaleGraph(graphPCMPi0Mass, 1000.);
            graphPCMPi0FWHM                             = (TGraphAsymmErrors*)directoryPCMPi0->Get("Pi0_Width_data");
            graphPCMPi0FWHM                             = ScaleGraph(graphPCMPi0FWHM, 1000.);
            graphPCMPi0MassMC                           = (TGraphAsymmErrors*)directoryPCMPi0->Get("Pi0_Mass_MC");
            graphPCMPi0MassMC                           = ScaleGraph(graphPCMPi0MassMC, 1000.);
            graphPCMPi0FWHMMC                           = (TGraphAsymmErrors*)directoryPCMPi0->Get("Pi0_Width_MC");
            graphPCMPi0FWHMMC                           = ScaleGraph(graphPCMPi0FWHMMC, 1000.);
            graphPCMPi0Acc                              = (TGraphAsymmErrors*)directoryPCMPi0->Get("AcceptancePi0");
            graphPCMPi0EffPt                            = (TGraphAsymmErrors*)directoryPCMPi0->Get("EfficiencyPi0");
            graphPCMPi0InvXSectionStat                  = (TGraphAsymmErrors*)directoryPCMPi0->Get("graphInvCrossSectionPi0");
            histoPCMPi0InvXSectionStat                  = (TH1D*)directoryPCMPi0->Get("InvCrossSectionPi0");
            cout << "Pi0 stat PCM" << endl;
            graphPCMPi0InvXSectionStat->Print();
            graphPCMPi0InvXSectionSys                   = (TGraphAsymmErrors*)directoryPCMPi0->Get("InvCrossSectionPi0Sys");
            cout << "Pi0 sys PCM" << endl;
            graphPCMPi0InvXSectionSys->Print();
            graphPCMPi0AccTimesEff                      = (TGraphAsymmErrors*)directoryPCMPi0->Get("EffTimesAccPi0");
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

            
        directoryPCMEta                                 = (TDirectory*)filePCM->Get("Eta2.76TeV"); 
            graphPCMEtaMass                             = (TGraphAsymmErrors*)directoryPCMEta->Get("Eta_Mass_data");
            graphPCMEtaMass                             = ScaleGraph(graphPCMEtaMass, 1000.);
            graphPCMEtaFWHM                             = (TGraphAsymmErrors*)directoryPCMEta->Get("Eta_Width_data");
            graphPCMEtaFWHM                             = ScaleGraph(graphPCMEtaFWHM, 1000.);
            graphPCMEtaMassMC                           = (TGraphAsymmErrors*)directoryPCMEta->Get("Eta_Mass_MC");
            graphPCMEtaMassMC                           = ScaleGraph(graphPCMEtaMassMC, 1000.);
            graphPCMEtaFWHMMC                           = (TGraphAsymmErrors*)directoryPCMEta->Get("Eta_Width_MC");
            graphPCMEtaFWHMMC                           = ScaleGraph(graphPCMEtaFWHMMC, 1000.);
            graphPCMEtaAcc                              = (TGraphAsymmErrors*)directoryPCMEta->Get("AcceptanceEta");
            graphPCMEtaEffPt                            = (TGraphAsymmErrors*)directoryPCMEta->Get("EfficiencyEta");
            graphPCMEtaAccTimesEff                      = (TGraphAsymmErrors*)directoryPCMEta->Get("EffTimesAccEta");
            histoPCMEtaInvXSectionStat                  = (TH1D*)directoryPCMEta->Get("InvCrossSectionEta");
            graphPCMEtaInvXSectionStat                  = (TGraphAsymmErrors*)directoryPCMEta->Get("graphInvCrossSectionEta");
            cout << "Eta stat PCM-EMC" << endl;
            graphPCMEtaInvXSectionStat->Print();
            graphPCMEtaInvXSectionSys                   = (TGraphAsymmErrors*)directoryPCMEta->Get("InvCrossSectionEtaSys");
            cout << "Eta sys PCM-EMC" << endl;
            graphPCMEtaInvXSectionSys->Print();        
            histoPCMEtaToPi0Stat                        = (TH1D*)directoryPCMEta->Get("EtaToPi0YShiftedStatError");
            graphPCMEtaToPi0Sys                         = (TGraphAsymmErrors*)directoryPCMEta->Get("EtaToPi0YShiftedSystError");
            
//             histoPCMEtaToPi0Stat                        = (TH1D*)directoryPCMEta->Get("EtaToPi0StatError");
//             graphPCMEtaToPi0Sys                         = (TGraphAsymmErrors*)directoryPCMEta->Get("EtaToPi0SystError");
    }
    
    
    
    //************************** Read data for PCMEMCAL **************************************************
    TFile* filePCMEMCAL                                     = new TFile(fileNamePCMEMCAL.Data());
    TH1D* histoPCMEMCALNumberOfEvents                       = (TH1D*)filePCMEMCAL->Get("histoNumberOfEvents2.76TeV");
    TDirectory* directoryPCMEMCALPi0                        = (TDirectory*)filePCMEMCAL->Get("Pi02.76TeV"); 
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
        for (Int_t i = 2; i < 6; i++){
            histoPCMEMCALPi0TriggerEff[i-2]                 = (TH1D*)directoryPCMEMCALPi0->Get(Form("TriggerEfficiencyPi0_%s",nameTrigger[i].Data()));
        }    
        TGraphAsymmErrors* graphPCMEMCALPi0InvXSectionStat  = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("graphInvCrossSectionPi0");
        TH1D* histoPCMEMCALPi0InvXSectionStat               = (TH1D*)directoryPCMEMCALPi0->Get("InvCrossSectionPi0");
        cout << "Pi0 stat PCM-EMC" << endl;
        graphPCMEMCALPi0InvXSectionStat->Print();
        TGraphAsymmErrors* graphPCMEMCALPi0InvXSectionSys   = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("InvCrossSectionPi0Sys");
        cout << "Pi0 sys PCM-EMC" << endl;
        graphPCMEMCALPi0InvXSectionSys->Print();
        TGraphAsymmErrors* graphPCMEMCALPi0AccTimesEff      = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("EffTimesAccPi0");
        for (Int_t k = 0; k < 4; k++){
            graphPi0EffSecCorrFromX[k][4]                   = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get(Form("EffectiveSecondaryPi0CorrFrom%s",nameSecPi0SourceRead[k].Data()));
            if (graphPi0EffSecCorrFromX[k][4]){
                haveEffSecCorr[k][4]                        = kTRUE;
            }
        }    

        
        TH1D* histoPi0InvMassSigPlusBGPCMEMCAL[6];
        TH1D* histoPi0InvMassSigPCMEMCAL[6];
        TH1D* histoPi0InvMassSigRemBGSubPCMEMCAL[6];
        TH1D* histoPi0InvMassBGPCMEMCAL[6];
        TH1D* histoPi0InvMassRemBGPCMEMCAL[6];
        TH1D* histoPi0InvMassBGTotPCMEMCAL[6];
        TF1* fitPi0InvMassSigPCMEMCAL[6];
        TF1* fitPi0InvMassBGPCMEMCAL[6];
        Bool_t haveAllPi0InvMassPCMEMCAL[6]                 = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
        if (plotInvMassBins){
            for (Int_t i = 0; i < 6; i++){
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
                    fitPi0InvMassBGPCMEMCAL[i]                                  = new TF1("Linearpp","[0]+[1]*x",0.0,0.3);
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

        
    TDirectory* directoryPCMEMCALEta                        = (TDirectory*)filePCMEMCAL->Get("Eta2.76TeV"); 
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
        TGraphAsymmErrors* graphPCMEMCALEtaToPi0Sys         = (TGraphAsymmErrors*)directoryPCMEMCALEta->Get("EtaToPi0YShiftedSystError");

        TH1D* histoEtaInvMassSigPlusBGPCMEMCAL[6];
        TH1D* histoEtaInvMassSigPCMEMCAL[6];
        TH1D* histoEtaInvMassSigRemBGSubPCMEMCAL[6];
        TH1D* histoEtaInvMassBGPCMEMCAL[6];
        TH1D* histoEtaInvMassRemBGPCMEMCAL[6];
        TH1D* histoEtaInvMassBGTotPCMEMCAL[6];
        TF1* fitEtaInvMassSigPCMEMCAL[6];
        TF1* fitEtaInvMassBGPCMEMCAL[6];
        Bool_t haveAllEtaInvMassPCMEMCAL[6]                 = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
        if (plotInvMassBins){
            for (Int_t i = 0; i < 6; i++){
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
                    fitEtaInvMassBGPCMEMCAL[i]                                  = new TF1("Linearpp","[0]+[1]*x",0.0,0.3);
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
                    histoEtaInvMassBGTotPCMEMCAL[i]         = (TH1D*)histoEtaInvMassBGPCMEMCAL[i]->Clone(Form("Eta_InvMassTotBG_Example_%s",nameTrigger[i].Data()));
                    histoEtaInvMassBGTotPCMEMCAL[i]->Sumw2();
                    histoEtaInvMassBGTotPCMEMCAL[i]->Add(histoEtaInvMassRemBGPCMEMCAL[i]);
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
    TDirectory* directoryEMCALPi0                           = (TDirectory*)fileEMCALLow->Get("Pi02.76TeV"); 
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
        for (Int_t k = 0; k < 4; k++){
            graphPi0EffSecCorrFromX[k][2]                   = (TGraphAsymmErrors*)directoryEMCALPi0->Get(Form("EffectiveSecondaryPi0CorrFrom%s",nameSecPi0SourceRead[k].Data()));
            if (graphPi0EffSecCorrFromX[k][2]){
                haveEffSecCorr[k][2]                        = kTRUE;
            }
        }    

        
        TH1D* histoPi0InvMassSigPlusBGEMCAL[6];
        TH1D* histoPi0InvMassSigEMCAL[6];
        TH1D* histoPi0InvMassSigRemBGSubEMCAL[6];
        TH1D* histoPi0InvMassBGEMCAL[6];
        TH1D* histoPi0InvMassRemBGEMCAL[6];
        TH1D* histoPi0InvMassBGTotEMCAL[6];
        TF1* fitPi0InvMassSigEMCAL[6];
        TF1* fitPi0InvMassBGEMCAL[6];
        Bool_t haveAllPi0InvMassEMCAL[6]                 = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
        if (plotInvMassBins){
            for (Int_t i = 0; i < 6; i++){
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
                    fitPi0InvMassBGEMCAL[i]                                  = new TF1("Linearpp","[0]+[1]*x",0.0,0.3);
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
        
    TDirectory* directoryEMCALEta                           = (TDirectory*)fileEMCALLow->Get("Eta2.76TeV"); 
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
        TGraphAsymmErrors* graphEMCALEtaToPi0Sys            = (TGraphAsymmErrors*)directoryEMCALEta->Get("EtaToPi0YShiftedSystError");

        TH1D* histoEtaInvMassSigPlusBGEMCAL[6];
        TH1D* histoEtaInvMassSigEMCAL[6];
        TH1D* histoEtaInvMassSigRemBGSubEMCAL[6];
        TH1D* histoEtaInvMassBGEMCAL[6];
        TH1D* histoEtaInvMassRemBGEMCAL[6];
        TH1D* histoEtaInvMassBGTotEMCAL[6];
        TF1* fitEtaInvMassSigEMCAL[6];
        TF1* fitEtaInvMassBGEMCAL[6];
        Bool_t haveAllEtaInvMassEMCAL[6]                 = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
        if (plotInvMassBins){
            for (Int_t i = 0; i < 6; i++){
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
                    fitEtaInvMassBGEMCAL[i]                                  = new TF1("Linearpp","[0]+[1]*x",0.0,0.3);
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
                    histoEtaInvMassBGTotEMCAL[i]         = (TH1D*)histoEtaInvMassBGEMCAL[i]->Clone(Form("Eta_InvMassTotBG_Example_%s",nameTrigger[i].Data()));
                    histoEtaInvMassBGTotEMCAL[i]->Sumw2();
                    histoEtaInvMassBGTotEMCAL[i]->Add(histoEtaInvMassRemBGEMCAL[i]);
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
    TDirectory* directoryEMCALmergedPi0                     = (TDirectory*)fileEMCALmerged->Get("Pi02.76TeV"); 
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
        if (flagMerged == 0){
            TH1D* histoEMCALMergedPi0Eff                    = (TH1D*)directoryEMCALmergedPi0->Get("Efficiencypi0");
//             histoEMCALMergedPi0Pur                          = (TH1D*)directoryEMCALmergedPi0->Get("PurityPi0");
            Double_t acceptanceEMCALmerged                  = 0.18;
            histoEMCALMergedPi0AccTimesEff                  = (TH1D*)histoEMCALMergedPi0Eff->Clone("histoEMCALMergedPi0AccTimesEff");
//          histoEMCALMergedPi0AccTimesEff->Multiply(histoEMCALMergedPi0Pur);
            histoEMCALMergedPi0AccTimesEff->Scale(acceptanceEMCALmerged);
            graphEMCALMergedPi0AccTimesEffDivPur                  = new TGraphAsymmErrors(histoEMCALMergedPi0AccTimesEff);
            histoEMCALMergedPi0InvXSectionStat              = (TH1D*)directoryEMCALmergedPi0->Get("InvCrossSection_With_stat");
            histoEMCALMergedPi0InvXSectionSys               = (TH1D*)directoryEMCALmergedPi0->Get("InvCrossSection_With_syst");
            graphEMCALMergedPi0InvXSectionStat              = new TGraphAsymmErrors(histoEMCALMergedPi0InvXSectionStat);
            cout << "EMCAL merged stat" << endl;
            graphEMCALMergedPi0InvXSectionStat->Print();
            graphEMCALMergedPi0InvXSectionSys               = new TGraphAsymmErrors(histoEMCALMergedPi0InvXSectionSys);
            cout << "EMCAL merged sys" << endl;
            graphEMCALMergedPi0InvXSectionSys->Print();
            for (Int_t i = 0; i< 8; i++){
                graphEMCALMergedPi0InvXSectionStat->RemovePoint(0);
                graphEMCALMergedPi0InvXSectionSys->RemovePoint(0);
            }
            Int_t l = 1;
//         while (histoEMCALMergedPi0InvXSectionStat->GetBinCenter(l) < 10){
//             histoEMCALMergedPi0InvXSectionStat->SetBinContent(l,0.);
//             histoEMCALMergedPi0InvXSectionStat->SetBinError(l,0.);
//         }    
        } else {
//             histoEMCALMergedPi0Eff                          = (TH1D*)directoryEMCALmergedPi0->Get("EfficiencyPi0");
//             histoEMCALMergedPi0Pur                          = (TH1D*)directoryEMCALmergedPi0->Get("PurityPi0");
            graphEMCALMergedPi0AccTimesEff                  = (TGraphAsymmErrors*)directoryEMCALmergedPi0->Get("EffTimesAccPi0");
            graphEMCALMergedPi0Purity                       = (TGraphAsymmErrors*)directoryEMCALmergedPi0->Get("PurityPi0");
            while (graphEMCALMergedPi0AccTimesEff->GetX()[0] < 16 ) 
                graphEMCALMergedPi0AccTimesEff->RemovePoint(0);
            while (graphEMCALMergedPi0Purity->GetX()[0] < 16 ) 
                graphEMCALMergedPi0Purity->RemovePoint(0);
            graphEMCALMergedPi0AccTimesEffDivPur            = CalculateGraphErrRatioToOtherTGraphErr(graphEMCALMergedPi0AccTimesEff, graphEMCALMergedPi0Purity,kTRUE);
            graphEMCALMergedPi0InvXSectionStat              = (TGraphAsymmErrors*)directoryEMCALmergedPi0->Get("graphInvCrossSectionPi0");
            graphEMCALMergedPi0InvXSectionSys               = (TGraphAsymmErrors*)directoryEMCALmergedPi0->Get("InvCrossSectionPi0Sys");
            if (graphEMCALMergedPi0AccTimesEffDivPur->GetX()[0] < 16 ) graphEMCALMergedPi0AccTimesEffDivPur->RemovePoint(0);
            for (Int_t k = 0; k < 4; k++){
                graphPi0EffSecCorrFromX[k][9]               = (TGraphAsymmErrors*)directoryEMCALmergedPi0->Get(Form("EffectiveSecondaryPi0CorrFrom%s",nameSecPi0SourceRead[k].Data()));
                if (graphPi0EffSecCorrFromX[k][9]){
                    haveEffSecCorr[k][9]                    = kTRUE;
                }
            }    
            
            if (flagMerged == 1 || flagMerged == 4){
                cout << "EMCAL merged stat" << endl;
                for (Int_t i = 0; graphEMCALMergedPi0InvXSectionStat->GetX()[0]< 9; i++){
                    graphEMCALMergedPi0InvXSectionStat->RemovePoint(0);
                }
                graphEMCALMergedPi0InvXSectionStat->Print();
                cout << "EMCAL merged sys" << endl;
                for (Int_t i = 0; graphEMCALMergedPi0InvXSectionSys->GetX()[0]< 9; i++){
                    graphEMCALMergedPi0InvXSectionSys->RemovePoint(0);
                }
                graphEMCALMergedPi0InvXSectionSys->Print();
            } else if (flagMerged == 2){
                cout << "EMCAL merged stat" << endl;
                for (Int_t i = 0; graphEMCALMergedPi0InvXSectionStat->GetX()[0]< 16; i++){
                    graphEMCALMergedPi0InvXSectionStat->RemovePoint(0);
                }
                graphEMCALMergedPi0InvXSectionStat->Print();
                cout << "EMCAL merged sys" << endl;
                for (Int_t i = 0; graphEMCALMergedPi0InvXSectionSys->GetX()[0]< 16; i++){
                    graphEMCALMergedPi0InvXSectionSys->RemovePoint(0);
                }
                graphEMCALMergedPi0InvXSectionSys->Print();    
            } else if (flagMerged == 3){
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
            }    
        }
        
        TH1D* histoDataM02[6];
        TH1D* histoMCrecM02[6];
        TH1D* histoTruePi0M02[6];
        TH1D* histoTrueEtaM02[6];
        TH1D* histoTrueGammaM02[6];
        TH1D* histoTrueElectronM02[6];
        TH1D* histoTrueBGM02[6];
        TH1D* histoTruePi0PureMergedM02[6];
        TH1D* histoTruePi0PConvMergedM02[6];
        TH1D* histoTruePi0OneGammaM02[6];
        TH1D* histoTruePi0OneElectronM02[6];
        Bool_t haveAllPi0M02EMCm[6]                 = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
        if (plotInvMassBins){
            for (Int_t i = 0; i < 6; i++){
                histoDataM02[i]                     = (TH1D*)directoryEMCALmergedPi0->Get(Form("Data_M02_Example_%s",nameTrigger[i].Data()));
                histoMCrecM02[i]                    = (TH1D*)directoryEMCALmergedPi0->Get(Form("MCRec_M02_Example_%s",nameTrigger[i].Data()));
                histoTruePi0M02[i]                  = (TH1D*)directoryEMCALmergedPi0->Get(Form("TruePi0_M02_Example_%s",nameTrigger[i].Data()));
                histoTrueEtaM02[i]                  = (TH1D*)directoryEMCALmergedPi0->Get(Form("TrueEta_M02_Example_%s",nameTrigger[i].Data()));
                histoTrueGammaM02[i]                = (TH1D*)directoryEMCALmergedPi0->Get(Form("TrueGamma_M02_Example_%s",nameTrigger[i].Data()));
                histoTrueElectronM02[i]             = (TH1D*)directoryEMCALmergedPi0->Get(Form("TrueElectron_M02_Example_%s",nameTrigger[i].Data()));
                histoTrueBGM02[i]                   = (TH1D*)directoryEMCALmergedPi0->Get(Form("TrueBG_M02_Example_%s",nameTrigger[i].Data()));
                histoTruePi0PureMergedM02[i]        = (TH1D*)directoryEMCALmergedPi0->Get(Form("TruePi0PureMerged_M02_Example_%s",nameTrigger[i].Data()));
                histoTruePi0PConvMergedM02[i]       = (TH1D*)directoryEMCALmergedPi0->Get(Form("TruePi0PartConvMerged_M02_Example_%s",nameTrigger[i].Data()));
                histoTruePi0OneGammaM02[i]          = (TH1D*)directoryEMCALmergedPi0->Get(Form("TruePi0OneGamma_M02_Example_%s",nameTrigger[i].Data()));
                histoTruePi0OneElectronM02[i]       = (TH1D*)directoryEMCALmergedPi0->Get(Form("TruePi0OneElectron_M02_Example_%s",nameTrigger[i].Data()));

                if (histoDataM02[i] && histoMCrecM02[i] && histoTruePi0M02[i] && histoTrueEtaM02[i]){
                    haveAllPi0M02EMCm[i]            = kTRUE;
                    cout << nameTrigger[i].Data() << "\t found M02 hists" << endl;
                }
            }    
        }    
        
    //************************** Read data for PHOS *****************************************************
    
    TFile* filePHOS                                         = new TFile(fileNamePHOS);
    TDirectory* directoryPHOSPi0                            = (TDirectory*)filePHOS->Get("pp2760"); 
        TH1D* histoPHOSPi0InvXSectionStat                   = (TH1D*)directoryPHOSPi0->Get("hPi02760GeVStat");
        TH1D* histoPHOSPi0InvXSectionSys                    = (TH1D*)directoryPHOSPi0->Get("hPi02760GeVSys");
        histoPHOSPi0InvXSectionStat->Scale(xSection2760GeV*recalcBarn);
        histoPHOSPi0InvXSectionStat->SetBinContent(1,0);
        histoPHOSPi0InvXSectionStat->SetBinError(1,0);
        histoPHOSPi0InvXSectionSys->Scale(xSection2760GeV*recalcBarn);
        TGraphAsymmErrors* graphPHOSPi0InvXSectionStat      = new TGraphAsymmErrors(histoPHOSPi0InvXSectionStat);
        graphPHOSPi0InvXSectionStat->RemovePoint(0);
        graphPHOSPi0InvXSectionStat->Print();
        TGraphAsymmErrors* graphPHOSPi0InvXSectionSys       = new TGraphAsymmErrors(histoPHOSPi0InvXSectionSys);
        graphPHOSPi0InvXSectionSys->RemovePoint(0);
        TH1D* histoPHOSPi0Mass                              = (TH1D*)directoryPHOSPi0->Get("mass2_GS_Data");
        TH1D* histoPHOSPi0FWHMMeV                           = (TH1D*)directoryPHOSPi0->Get("width2_GS_Data");
        TH1D* histoPHOSPi0TrueMass                          = (TH1D*)directoryPHOSPi0->Get("mass2_GS_MC");
        TH1D* histoPHOSPi0TrueFWHMMeV                       = (TH1D*)directoryPHOSPi0->Get("width2_GS_MC");
        histoPHOSPi0Mass->Scale(1000.);
        histoPHOSPi0TrueMass->Scale(1000.);
        histoPHOSPi0FWHMMeV->Scale(1000.);
        histoPHOSPi0TrueFWHMMeV->Scale(1000.);
        histoPHOSPi0Mass->SetBinContent(histoPHOSPi0Mass->GetNbinsX(),0);
        histoPHOSPi0Mass->SetBinContent(1,0);
        histoPHOSPi0TrueMass->SetBinContent(histoPHOSPi0TrueMass->GetNbinsX(),0);
        histoPHOSPi0TrueMass->SetBinContent(1,0);
        histoPHOSPi0FWHMMeV->SetBinContent(histoPHOSPi0FWHMMeV->GetNbinsX(),10000.);
        histoPHOSPi0FWHMMeV->SetBinContent(1,10000.);
        histoPHOSPi0TrueFWHMMeV->SetBinContent(histoPHOSPi0TrueFWHMMeV->GetNbinsX(),10000.);
        histoPHOSPi0TrueFWHMMeV->SetBinContent(1,10000.);
        TH1D* histoPHOSPi0AccTimesEff                       = (TH1D*)directoryPHOSPi0->Get("yield1_int_CB");
        histoPHOSPi0AccTimesEff->Scale(2*TMath::Pi());
        TGraphAsymmErrors* graphPHOSPi0AccTimesEff          = new TGraphAsymmErrors(histoPHOSPi0AccTimesEff);


    // *******************************************************************************************************
    // ************************** Loading theory calculations ************************************************
    // *******************************************************************************************************        
    TFile* fileTheoryCompilation                            = new TFile(fileNameTheory.Data());
        TH1F* histoPythia8InvXSectionPi0                    = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoPi02760GeV");
        TGraphErrors* graphPythia8InvXSectionPi0            = new TGraphErrors(histoPythia8InvXSectionPi0);
        while(graphPythia8InvXSectionPi0->GetX()[0] < 0.3) graphPythia8InvXSectionPi0->RemovePoint(0);
        TH1F* histoPythia8InvXSectionEta                    = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Monash2013LegoEta2760GeV");
        TGraphErrors* graphPythia8InvXSectionEta            = new TGraphErrors(histoPythia8InvXSectionEta);
        while(graphPythia8InvXSectionEta->GetX()[0] < 0.4) graphPythia8InvXSectionEta->RemovePoint(0);
        while(graphPythia8InvXSectionEta->GetX()[graphPythia8InvXSectionEta->GetN()-1] > 25) graphPythia8InvXSectionEta->RemovePoint(graphPythia8InvXSectionEta->GetN()-1);
        TH1F* histoPythia8EtaToPi0                          = (TH1F*) fileTheoryCompilation->Get("histoEtaToPi0RatioPythia8Monash2013Lego2760GeV");
        TGraphErrors* graphPythia8EtaToPi0                  = new TGraphErrors(histoPythia8EtaToPi0);
        while(graphPythia8EtaToPi0->GetX()[0] < 0.4) graphPythia8EtaToPi0->RemovePoint(0);
        while(graphPythia8EtaToPi0->GetX()[graphPythia8EtaToPi0->GetN()-1] > 27) graphPythia8EtaToPi0->RemovePoint(graphPythia8EtaToPi0->GetN()-1);
        TGraph* graphNLOCalcPi0MuHalf                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuHalf2760GeV");
        TGraph* graphNLOCalcPi0MuOne                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuOne2760GeV");
        TGraph* graphNLOCalcPi0MuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuTwo2760GeV");
        TGraph* graphNLOCalcEtaMuHalf                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuHalf2760GeV");
        TGraph* graphNLOCalcEtaMuOne                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuOne2760GeV");
        TGraph* graphNLOCalcEtaMuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuTwo2760GeV");
        TGraphAsymmErrors* graphNLOEtaAESSCalc              = (TGraphAsymmErrors*)fileTheoryCompilation->Get("graphNLOCalcAESSSInvSecEta2760GeV");
        TGraphAsymmErrors* graphNLOEtaToPi0                 = (TGraphAsymmErrors*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi02760GeV_AESSS_DSS07");
        TGraph* graphNLOEtaToPi0MuHalf                      = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuHalf2760GeV");
        TGraph* graphNLOEtaToPi0MuOne                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuOne2760GeV");
        TGraph* graphNLOEtaToPi0MuTwo                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuTwo2760GeV");
        TGraph* graphNLODSSCalcMuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcDSSInvSecPi0MuTwo2760GeV");
        TGraph* graphNLOCGCCalcMuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcCGCInvCrossSec2760GeV");
        TGraphAsymmErrors* graphNLODSS14Calc                = (TGraphAsymmErrors*)fileTheoryCompilation->Get("graphNLOCalcDSS14InvCrossSec2760GeV");
        TGraph* graphCGCPi0                                 = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcCGCInvCrossSec2760GeV");
        
        TDirectory* directoryDirectPhotonTheory             = (TDirectory*)fileTheoryCompilation->Get("DirectPhoton");
        TGraphAsymmErrors* graphNLOCalcDirGam               = (TGraphAsymmErrors*)directoryDirectPhotonTheory->Get("graphDirectPhotonNLOVogelsang_2760GeV");
        
        while (graphNLOCalcEtaMuHalf->GetX()[graphNLOCalcEtaMuHalf->GetN()-1] > 22. )
            graphNLOCalcEtaMuHalf->RemovePoint(graphNLOCalcEtaMuHalf->GetN()-1);
        while (graphNLOCalcEtaMuOne->GetX()[graphNLOCalcEtaMuOne->GetN()-1] > 22. )
            graphNLOCalcEtaMuOne->RemovePoint(graphNLOCalcEtaMuOne->GetN()-1);
        while (graphNLOCalcEtaMuTwo->GetX()[graphNLOCalcEtaMuTwo->GetN()-1] > 22. )
            graphNLOCalcEtaMuTwo->RemovePoint(graphNLOCalcEtaMuTwo->GetN()-1);
            
    // *******************************************************************************************************
    // ************************** Publication pi0 PbPb *******************************************************
    // *******************************************************************************************************
    TFile* filePublishedPi0Old                              = new TFile(fileNameOldPi0Publication);
        TF1* fitPi0OldPub                                   = (TF1*)filePublishedPi0Old->Get("fitInvCrossSectionPi0Comb2760GeV_YShift");
        TGraphAsymmErrors* graphPi0PubOldTot                = (TGraphAsymmErrors*)filePublishedPi0Old->Get("graphInvCrossSectionPi0Comb2760GeV");
        TGraphAsymmErrors* graphPi0PubOldSys                = (TGraphAsymmErrors*)filePublishedPi0Old->Get("graphInvCrossSectionPi0Comb2760GeVSysErr");
        TGraphAsymmErrors* graphPi0PubOldStat               = (TGraphAsymmErrors*)filePublishedPi0Old->Get("graphInvCrossSectionPi0Comb2760GeVStatErr");
        
    // *******************************************************************************************************
    // ************************** Loading eta/pi0 compilation ************************************************
    // *******************************************************************************************************        
    TFile* fileEtaToPi0Compilation                          = new TFile(fileNameEtaToPi0WorldData.Data());
        TGraphErrors* graphPHENIXEtaToPi0200GeV             = (TGraphErrors*)fileEtaToPi0Compilation->Get("Phenix200GeV");
        TGraphAsymmErrors* graphALICEEtaToPi07TeV           = (TGraphAsymmErrors*)fileEtaToPi0Compilation->Get("Alice7TeV");
        
//         return;
    // *******************************************************************************************************
    // ************************** Loading charged pion results ***********************************************
    // *******************************************************************************************************        
    TFile* fileChargedPionInputpp                           = new TFile(fileNameChargedPionPP.Data());
        TH1D* histoChPiInvYieldPubStatPP                    = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecPubStat2760GeV");
        TH1D* histoChPiInvYieldPubSystPP                    = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecPubSyst2760GeV");

        TH1D* histoChPiInvYieldHighPtStatPP                 = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecHighPtStat2760GeV");
        TH1D* histoChPiInvYieldHighPtSystPP                 = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecHighPtSyst2760GeV");
        TH1D* histoChPiInvYieldLowPtStatCMS                 = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtStat2760GeVCMS");
        TH1D* histoChPiInvYieldLowPtSysCMS                  = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtSys2760GeVCMS");
        TH1D* histoChPiInvYieldLowPtStatPP                  = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtStatPP2760GeV");
        TH1D* histoChPiInvYieldLowPtSysPP                   = (TH1D*)fileChargedPionInputpp->Get("histoChargedPionSpecLowPtSysPP2760GeV");

    // *******************************************************************************************************
    // ************************** Loading charged hadron results ***********************************************
    // *******************************************************************************************************        
    TFile* fileChargedHadronsInputpp                        = new TFile(fileNameChargedHadronPP.Data());
        TGraphAsymmErrors* graphChargedHadronsStatPP        = (TGraphAsymmErrors*)fileChargedHadronsInputpp->Get("graphChargedHadronsALICEStatPP2760GeV");
        graphChargedHadronsStatPP                           = ScaleGraph(graphChargedHadronsStatPP,0.5*1e9/xSection2760GeVINEL);
        TGraphAsymmErrors* graphChargedHadronsSysPP         = (TGraphAsymmErrors*)fileChargedHadronsInputpp->Get("graphChargedHadronsALICESysPP2760GeV");
        graphChargedHadronsSysPP                            = ScaleGraph(graphChargedHadronsSysPP,0.5*1e9/xSection2760GeVINEL);
        TGraphAsymmErrors* graphChargedHadronsStatPPCMS     = (TGraphAsymmErrors*)fileChargedHadronsInputpp->Get("graphChargedHadronsCMSStatPP2760GeV");
        graphChargedHadronsStatPPCMS                        = ScaleGraph(graphChargedHadronsStatPPCMS,0.5*1e9/xSection2760GeVINEL);
        TGraphAsymmErrors* graphChargedHadronsSysPPCMS      = (TGraphAsymmErrors*)fileChargedHadronsInputpp->Get("graphChargedHadronsCMSSysPP2760GeV");
        graphChargedHadronsSysPPCMS                         = ScaleGraph(graphChargedHadronsSysPPCMS,0.5*1e9/xSection2760GeVINEL);
        TGraphAsymmErrors* graphChargedHadronsStatPPATLAS   = (TGraphAsymmErrors*)fileChargedHadronsInputpp->Get("graphChargedHadronsATLASStatPP2760GeV");
        graphChargedHadronsStatPPATLAS                      = ScaleGraph(graphChargedHadronsStatPPATLAS,0.5*1e9/xSection2760GeVINEL);
        TGraphAsymmErrors* graphChargedHadronsSysPPATLAS    = (TGraphAsymmErrors*)fileChargedHadronsInputpp->Get("graphChargedHadronsATLASSysPP2760GeV");
        graphChargedHadronsSysPPATLAS                       = ScaleGraph(graphChargedHadronsSysPPATLAS,0.5*1e9/xSection2760GeVINEL);
        
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

    
    // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)    
    // for systematic error from respective method
    TGraphAsymmErrors* sysErrorCollectionPi0[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorCollectionPi0[i]    = NULL;
    }    
    sysErrorCollectionPi0[0]        = (TGraphAsymmErrors*)graphPCMPi0InvXSectionSys->Clone("sysErrPCMPi0");
    sysErrorCollectionPi0[1]        = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionSys->Clone("sysErrPHOSPi0");

    // Definition of final pt binning (has to be set manually)
    Double_t xPtLimitsPi0[100];
    Int_t maxNBinsPi0               = 33;
    
    if (flagMerged == 0) {
        maxNBinsPi0                  = GetBinning( xPtLimitsPi0, "Pi0", "2.76TeV", 20 );
    } else {
        maxNBinsPi0                  = GetBinning( xPtLimitsPi0, "Pi0", "2.76TeV", 10 );
//         if (flagMerged == 1) 
            maxNBinsPi0--;
    }    
    Double_t xPtLimitsPi0WOMerged[50];
    Int_t maxNBinsPi0W0Merged       = GetBinning( xPtLimitsPi0WOMerged, "Pi0", "2.76TeV", 4 );
    
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsPi0[11]            = { 0,  2,  0,  0,  0,
                                        0,  0,  0,  0,  0, 
                                        0};
    Int_t offSetsPi0Sys[11]         = { 1,  3,  6,  0,  3,
                                        0,  0,  0,  0,  21, 
                                        0};
    Int_t offSetPi0Shifting[11]     = { 0,  2,  5,  0,  2,
                                        0,  0,  0,  0,  21,
                                        0 };
    Int_t nComBinsPi0Shifting[11]   = { 16, 16, 24, 0,  20,
                                        0,  0,  0,  0,  31,
                                        0 };
                                        
    if (flagMerged == 1 || flagMerged == 4) {
        offSetsPi0[9]               = 0;
        offSetsPi0Sys[9]            = 20;
        offSetPi0Shifting[9]        = 19;
        nComBinsPi0Shifting[9]      = 30;
    } else if (flagMerged == 2 ) {
        offSetsPi0[9]               = 0;
        offSetsPi0Sys[9]            = 25;   
        offSetPi0Shifting[9]        = 24;
        nComBinsPi0Shifting[9]      = 30;
    } else if (flagMerged == 3) {
        offSetsPi0[9]               = 0;
        offSetsPi0Sys[9]            = 21;
        offSetPi0Shifting[9]        = 20;
        nComBinsPi0Shifting[9]      = 30;
    }
                                        
    // **********************************************************************************************************************
    // ******************************************* Recalculating published spectrum *****************************************
    // **********************************************************************************************************************
    TGraph* graphWeightsPi0Old[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsPi0Old[i]          = NULL;
    }    
                        
    // Declaration & calculation of combined spectrum                            
    TString fileNamePi0OutputWeightingOld                          = Form("%s/Pi0_WeightingOld.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombPi0InvXSectionStatOld    = NULL;
    TGraphAsymmErrors* graphCombPi0InvXSectionSysOld     = NULL;
    TGraphAsymmErrors* graphCombPi0InvXSectionTotOld     = CombinePtPointsSpectraFullCorrMat(   statErrorCollectionPi0,    sysErrorCollectionPi0,     
                                                                                                xPtLimitsPi0WOMerged, 17,
                                                                                                offSetsPi0, offSetsPi0Sys,
                                                                                                graphCombPi0InvXSectionStatOld, graphCombPi0InvXSectionSysOld,
                                                                                                fileNamePi0OutputWeightingOld, "2.76TeV", "Pi0", kTRUE, 
                                                                                                NULL, fileNameCorrFactors
                                                                                             );
    if (graphCombPi0InvXSectionTotOld == NULL) {
        cout << "Aborting: something went wrong during the combination of the old spectra" << endl;
    }    
    graphCombPi0InvXSectionStatOld->Print();
    
    // Reading weights from output file for plotting
    ifstream fileWeightsPi0ReadOld;
    fileWeightsPi0ReadOld.open(fileNamePi0OutputWeightingOld,ios_base::in);
    cout << "reading" << fileNamePi0OutputWeightingOld << endl;
    Double_t xValuesPi0ReadOld[50];
    Double_t weightsPi0ReadOld[11][50];
    Int_t availablePi0MeasOld[11]      = { -1, -1, -1, -1, -1,
                                        -1, -1, -1, -1, -1, 
                                        -1};
    Int_t nMeasSetPi0Old               = 2;
    Int_t nPtBinsPi0ReadOld            = 0;
    while(!fileWeightsPi0ReadOld.eof() && nPtBinsPi0ReadOld < 50){
        TString garbage = "";
        if (nPtBinsPi0ReadOld == 0){
            fileWeightsPi0ReadOld >> garbage ;//>> availablePi0Meas[0] >> availablePi0Meas[1] >> availablePi0Meas[2] >> availablePi0Meas[3];
            for (Int_t i = 0; i < nMeasSetPi0Old; i++){
                fileWeightsPi0ReadOld >> availablePi0MeasOld[i] ;
            }    
            cout << "read following measurements: "; 
            for (Int_t i = 0; i < 11; i++){
                cout << availablePi0MeasOld[i] << "\t" ;
            }    
            cout << endl;
        } else {
            fileWeightsPi0ReadOld >> xValuesPi0ReadOld[nPtBinsPi0ReadOld-1];
            for (Int_t i = 0; i < nMeasSetPi0Old; i++){
                fileWeightsPi0ReadOld >> weightsPi0ReadOld[availablePi0MeasOld[i]][nPtBinsPi0ReadOld-1] ;
            }    
            cout << "read: "<<  nPtBinsPi0ReadOld << "\t"<< xValuesPi0ReadOld[nPtBinsPi0ReadOld-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetPi0Old; i++){
                cout << weightsPi0ReadOld[availablePi0MeasOld[i]][nPtBinsPi0ReadOld-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsPi0ReadOld++;
    }
    nPtBinsPi0ReadOld                  = nPtBinsPi0ReadOld-2 ;
    fileWeightsPi0ReadOld.close();

    for (Int_t i = 0; i < nMeasSetPi0Old; i++){
        graphWeightsPi0Old[availablePi0MeasOld[i]]                        = new TGraph(nPtBinsPi0ReadOld,xValuesPi0ReadOld,weightsPi0ReadOld[availablePi0MeasOld[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsPi0ReadOld; n++){
            if (graphWeightsPi0Old[availablePi0MeasOld[i]]->GetY()[bin] == 0) graphWeightsPi0Old[availablePi0MeasOld[i]]->RemovePoint(bin);
            else bin++;
        }    
    }    

    // **********************************************************************************************************************
    // ******************************************* Plotting weights method only EMC *****************************************
    // **********************************************************************************************************************
    Int_t textSizeLabelsPixel           = 900*0.04;

    TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
    canvasWeights->SetLogx();
   
    TH2F * histo2DPi0Weights;
    histo2DPi0Weights = new TH2F("histo2DPi0Weights","histo2DPi0Weights",11000,0.23,70.,1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DPi0Weights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DPi0Weights->GetXaxis()->SetMoreLogLabels();
    histo2DPi0Weights->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DPi0Weights->GetYaxis()->SetRangeUser(-10,10);
    histo2DPi0Weights->Draw("copy");
    
        TLegend* legendWeightsOld    = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetPi0Old), 32);
        for (Int_t i = 0; i < nMeasSetPi0Old; i++){
            DrawGammaSetMarkerTGraph(graphWeightsPi0Old[availablePi0MeasOld[i]],
                                     markerStyleDet[availablePi0MeasOld[i]], 
                                     markerSizeDet[availablePi0MeasOld[i]]*0.5, 
                                     colorDet[availablePi0MeasOld[i]] , 
                                     colorDet[availablePi0MeasOld[i]]);
            graphWeightsPi0Old[availablePi0MeasOld[i]]->Draw("p,same,z");
            legendWeightsOld->AddEntry(graphWeightsPi0Old[availablePi0MeasOld[i]],nameMeasGlobalLabel[availablePi0MeasOld[i]],"p");
        }    
        legendWeightsOld->Draw();

        TLatex *labelWeightsEnergy      = new TLatex(0.7,0.20,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel, 4, 1, 43);
        labelWeightsEnergy->Draw();
        TLatex *labelWeightsPi0         = new TLatex(0.7,0.16,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelWeightsPi0, 0.85*textSizeLabelsPixel, 4, 1, 43 );
        labelWeightsPi0->Draw();

//      DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
        
    canvasWeights->SaveAs(Form("%s/Pi0_WeightsOldPublished.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************************* Adding EMCal & EMC merged  ***********************************************
    // **********************************************************************************************************************
    TH1D* histoEMCALMergedPi0InvXSectionStatCorrBin     = new TH1D("histoEMCALMergedPi0InvXSectionStatCorrBin", "", 33, xPtLimitsPi0);
    Int_t firstBinMergedPi0 = 1;
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
    for (Int_t i = 1; i < histoEMCALMergedPi0InvXSectionStatCorrBin->GetNbinsX(); i++){
        cout << "Bin " << i << "\t" <<  histoEMCALMergedPi0InvXSectionStatCorrBin->GetBinCenter(i) << "\t" << histoEMCALMergedPi0InvXSectionStatCorrBin->GetBinContent(i) 
             << "\t" << histoEMCALMergedPi0InvXSectionStatCorrBin->GetBinError(i) << endl;
    }    
    
//     return;

    statErrorCollectionPi0[2]          = (TH1D*)histoEMCALPi0InvXSectionStat->Clone("statErrEMCALPi0");
    sysErrorCollectionPi0[2]           = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionSys->Clone("sysErrEMCALPi0");
    statErrorCollectionPi0[9]          = (TH1D*)histoEMCALMergedPi0InvXSectionStatCorrBin->Clone("statErrEMCALMergedPi0");
    sysErrorCollectionPi0[9]           = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionSys->Clone("sysErrEMCALMergedPi0");
    
    // **********************************************************************************************************************
    // ************************ Adding PCM-EMC measurements to matrix & calculating relativ errors ************************
    // **********************************************************************************************************************
    
    statErrorCollectionPi0[4]               = (TH1D*)histoPCMEMCALPi0InvXSectionStat->Clone("statErrPCMEMCALPi0");
    sysErrorCollectionPi0[4]                = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionSys->Clone("sysErrPCMEMCALPi0");

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
        if (sysErrorCollectionPi0[i]) 
            sysErrorRelCollectionPi0[i]     = CalculateRelErrUpAsymmGraph( sysErrorCollectionPi0[i], Form("relativeSysErrorPi0_%s", nameMeasGlobal[i].Data()));
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
                                                                                                fileNamePi0OutputWeightingA, "2.76TeV", "Pi0", kTRUE, 
                                                                                                NULL, fileNameCorrFactors
                                                                                            );
    
    
    if (graphCombPi0InvXSectionTotA == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }
    
    while (graphCombPi0InvXSectionStatA->GetX()[0] < 0.4){
        graphCombPi0InvXSectionStatA->RemovePoint(0);
    }
    while (graphCombPi0InvXSectionTotA->GetX()[0] < 0.4){    
        graphCombPi0InvXSectionTotA->RemovePoint(0);
    }
    while (graphCombPi0InvXSectionSysA->GetX()[0] < 0.4){    
        graphCombPi0InvXSectionSysA->RemovePoint(0);
    }    
    graphCombPi0InvXSectionTotA->Print();
    
    // Reading weights from output file for plotting
    ifstream fileWeightsPi0ReadA;
    fileWeightsPi0ReadA.open(fileNamePi0OutputWeightingA,ios_base::in);
    cout << "reading" << fileNamePi0OutputWeightingA << endl;
    Double_t xValuesPi0ReadA[50];
    Double_t weightsPi0ReadA[11][50];
    Int_t availablePi0MeasA[11]        = { -1, -1, -1, -1, -1,
                                        -1, -1, -1, -1, -1,
                                        -1};
    Int_t nMeasSetPi0A                 = 5;
    Int_t nPtBinsPi0ReadA              = 0;
    while(!fileWeightsPi0ReadA.eof() && nPtBinsPi0ReadA < 50){
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
    canvasWeights->cd();
   
    histo2DPi0Weights->Draw("copy");

        TLegend* legendWeights   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetPi0A), 32);
        for (Int_t i = 0; i < nMeasSetPi0A; i++){
            DrawGammaSetMarkerTGraph(graphWeightsPi0A[availablePi0MeasA[i]], markerStyleDet[availablePi0MeasA[i]], markerSizeDet[availablePi0MeasA[i]]*0.5, colorDet[availablePi0MeasA[i]] , colorDet[availablePi0MeasA[i]]);
            graphWeightsPi0A[availablePi0MeasA[i]]->Draw("p,same,z");
            legendWeights->AddEntry(graphWeightsPi0A[availablePi0MeasA[i]],nameMeasGlobalLabel[availablePi0MeasA[i]],"p");
        }    
        legendWeights->Draw();

        labelWeightsEnergy->Draw();
        labelWeightsPi0->Draw();

//      DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
        
    canvasWeights->SaveAs(Form("%s/Pi0_WeightsA.%s",outputDir.Data(),suffix.Data()));

// 
    // **********************************************************************************************************************
    // ******************************************* Compare new spectrum to old average **************************************
    // **********************************************************************************************************************
    
    TGraphAsymmErrors* graphCombPi0InvXSectionStatAForCompToOld      = (TGraphAsymmErrors*)graphCombPi0InvXSectionStatA->Clone("graphCombPi0InvCrossSectionStat2760GeVAForCompToOld");
    TGraphAsymmErrors* graphCombPi0InvXSectionSysAForCompToOld       = (TGraphAsymmErrors*)graphCombPi0InvXSectionSysA->Clone("graphCombPi0InvCrossSectionSys2760GeVAForCompToOld");
    
    for (Int_t n = graphCombPi0InvXSectionStatAForCompToOld->GetN()-1; graphCombPi0InvXSectionStatAForCompToOld->GetX()[n] > 6; n-- ){
        graphCombPi0InvXSectionStatAForCompToOld->RemovePoint(n);
    }    
    for (Int_t n = graphCombPi0InvXSectionSysAForCompToOld->GetN()-1; graphCombPi0InvXSectionSysAForCompToOld->GetX()[n] > 6; n-- ){
        graphCombPi0InvXSectionSysAForCompToOld->RemovePoint(n);
    }    
    
    graphCombPi0InvXSectionStatAForCompToOld->Print();
    graphCombPi0InvXSectionSysAForCompToOld->Print();
    
    TGraphAsymmErrors*  graphRatioPi0CombNewADivCombOldStat        = CalculateGraphErrRatioToOtherTGraphErr(   graphCombPi0InvXSectionStatAForCompToOld, 
                                                                                                            graphCombPi0InvXSectionStatOld, kTRUE);
    TGraphAsymmErrors*  graphRatioPi0CombNewADivCombOldSys         = CalculateGraphErrRatioToOtherTGraphErr(   graphCombPi0InvXSectionSysAForCompToOld, 
                                                                                                            graphCombPi0InvXSectionSysOld, kTRUE);
    

    // **********************************************************************************************************************
    // ******************************************* Ratio of Comb to Fit ****************************************
    // **********************************************************************************************************************
    
    TCanvas* canvasRatioToOldCombined   = new TCanvas("canvasRatioToOldCombined","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToOldCombined, 0.1, 0.02, 0.035, 0.09);
    canvasRatioToOldCombined->SetLogx();
   
    TH2F * histo2DRatioToOldCombined;
    histo2DRatioToOldCombined           = new TH2F("histo2DRatioToOldCombined","histo2DRatioToOldCombined",11000,0.23,10.,1000,0.75,1.25);
    SetStyleHistoTH2ForGraphs(histo2DRatioToOldCombined, "#it{p}_{T} (GeV/#it{c})","#frac{Comb}{Comb Old}",0.035,0.04, 0.035,0.04, 1.,1.,510,505);
    histo2DRatioToOldCombined->GetXaxis()->SetMoreLogLabels();
    histo2DRatioToOldCombined->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRatioToOldCombined->GetYaxis()->SetRangeUser(-10,10);
    histo2DRatioToOldCombined->Draw("copy");

 
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNewADivCombOldSys, markerStyleComb+4, markerSizeComb, kBlue+2 , kBlue+2, widthLinesBoxes, kTRUE);
        graphRatioPi0CombNewADivCombOldSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNewADivCombOldStat, markerStyleComb+4, markerSizeComb, kBlue+2 , kBlue+2);
        graphRatioPi0CombNewADivCombOldStat->Draw("p,same,z");

        DrawGammaLines(0.23, 10. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.23, 10. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.23, 10. , 0.9, 0.9,0.1, kGray, 7);
        DrawGammaLines(0.23, 10. , 1.05, 1.05,0.1, kGray, 3);
        DrawGammaLines(0.23, 10. , 0.95, 0.95,0.1, kGray, 3);
        
        TLegend* legendRatioToOld       = GetAndSetLegend2(0.67, 0.93-(0.035*1), 0.93, 0.93, 32);
        legendRatioToOld->AddEntry(graphRatioPi0CombNewADivCombOldSys,"All");
        legendRatioToOld->Draw();

        TLatex *labelRatioToOldEnergy   = new TLatex(0.15,0.89,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelRatioToOldEnergy, 0.85*textSizeLabelsPixel, 4, 1, 43);
        labelRatioToOldEnergy->Draw();
        TLatex *labelRatioToOldPi0      = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToOldPi0, 0.85*textSizeLabelsPixel, 4, 1, 43);
        labelRatioToOldPi0->Draw();
        
    canvasRatioToOldCombined->SaveAs(Form("%s/Pi0_RatioOfCombToCombOld_PP2760GeV.%s",outputDir.Data(),suffix.Data()));
    
    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************
    
    TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErr->SetLogx();
   
    TH2F * histo2DRelSysErr;
    histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.23,70.,1000,0,80.5);
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
            legendRelSysErr->AddEntry(sysErrorRelCollectionPi0[availablePi0MeasA[i]],nameMeasGlobalLabel[availablePi0MeasA[i]],"p");
        }    
        legendRelSysErr->Draw();

        TLatex *labelRelSysErrEnergy    = new TLatex(0.15,0.89,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelRelSysErrEnergy, 0.85*textSizeLabelsPixel, 4, 1, 43);
        labelRelSysErrEnergy->Draw();
        TLatex *labelRelSysErrPi0       = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelSysErrPi0, 0.85*textSizeLabelsPixel, 4, 1, 43);
        labelRelSysErrPi0->Draw();
        
    canvasRelSysErr->SaveAs(Form("%s/Pi0_RelSysErr.%s",outputDir.Data(),suffix.Data()));
    
    histo2DRelSysErr->GetYaxis()->SetRangeUser(0,38.5);
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
    histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.23,70.,1000,0,80.5);
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
            legendRelStatErr->AddEntry(statErrorRelCollectionPi0[availablePi0MeasA[i]],nameMeasGlobalLabel[availablePi0MeasA[i]],"p");
        }    
        legendRelStatErr->Draw();

        TLatex *labelRelStatErrEnergy   = new TLatex(0.75,0.89,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelRelStatErrEnergy, 0.85*textSizeLabelsPixel, 4, 1, 43);
        labelRelStatErrEnergy->Draw();
        TLatex *labelRelStatErrPi0      = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelStatErrPi0, 0.85*textSizeLabelsPixel, 4, 1, 43);
        labelRelStatErrPi0->Draw();
        
    canvasRelStatErr->SaveAs(Form("%s/Pi0_RelStatErr.%s",outputDir.Data(),suffix.Data()));

    histo2DRelStatErr->GetYaxis()->SetRangeUser(0,38.5);
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
    TGraphAsymmErrors* graphCombPi0InvXSectionRelTotOld      = CalculateRelErrUpAsymmGraph( graphCombPi0InvXSectionTotOld, "relativeTotalErrorPi0_Old");

    while (graphCombPi0InvXSectionRelStatA->GetX()[graphCombPi0InvXSectionRelStatA->GetN()-1] > 40 )graphCombPi0InvXSectionRelStatA->RemovePoint(graphCombPi0InvXSectionRelStatA->GetN()-1);
    while (graphCombPi0InvXSectionRelSysA->GetX()[graphCombPi0InvXSectionRelSysA->GetN()-1] > 40 )graphCombPi0InvXSectionRelSysA->RemovePoint(graphCombPi0InvXSectionRelSysA->GetN()-1);
    while (graphCombPi0InvXSectionRelTotA->GetX()[graphCombPi0InvXSectionRelTotA->GetN()-1] > 40 )graphCombPi0InvXSectionRelTotA->RemovePoint(graphCombPi0InvXSectionRelTotA->GetN()-1);
    
    TCanvas* canvasRelTotErr            = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelTotErr->SetLogx();
   
    TH2F * histo2DRelTotErrPi0;
    histo2DRelTotErrPi0                 = new TH2F("histo2DRelTotErrPi0","histo2DRelTotErrPi0",11000,0.23,70.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrPi0, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrPi0->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelTotErrPi0->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionRelTotOld, markerStyleComb+1, markerSizeComb, kBlack , kBlack);
        graphCombPi0InvXSectionRelTotOld->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionRelTotA, markerStyleComb+4, markerSizeComb, kBlue+2 , kBlue+2);
        graphCombPi0InvXSectionRelTotA->Draw("p,same,z");

        TLegend* legendRelTotErr = GetAndSetLegend2(0.14, 0.92-(0.035*2), 0.45, 0.92, 32);
        legendRelTotErr->AddEntry(graphCombPi0InvXSectionRelTotOld,"PCM, PHOS","p");
        legendRelTotErr->AddEntry(graphCombPi0InvXSectionRelTotA,"All","p");
        legendRelTotErr->Draw();

        TLatex *labelRelTotErrEnergy    = new TLatex(0.75,0.89,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelRelTotErrEnergy, 0.85*textSizeLabelsPixel, 4, 1, 43);        
        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrPi0       = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelTotErrPi0, 0.85*textSizeLabelsPixel, 4, 1, 43);
        labelRelTotErrPi0->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/Pi0_RelTotErr.%s",outputDir.Data(),suffix.Data()));
    
    
    histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,38.5);
    histo2DRelTotErrPi0->Draw("copy");
        graphCombPi0InvXSectionRelTotOld->Draw("p,same,z");
        graphCombPi0InvXSectionRelTotA->Draw("p,same,z");

        legendRelTotErr->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrPi0->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/Pi0_RelTotErrZoomed.%s",outputDir.Data(),suffix.Data()));

    histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,38.5);
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
// 
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
    TF1* fitTCMDecomposedLPi0   = FitObject("tcmlow","twoCompModelPi0_DecL", "Pi0", NULL, 0.4, 2.);
    TF1* fitTCMDecomposedHPi0   = FitObject("tcmhigh","twoCompModelPi0_DecH", "Pi0", NULL, 4, 50.);
    fitTCMDecomposedLPi0->SetParameters(graphCombPi0InvXSectionTotA->GetY()[2],0.3);
    graphCombPi0InvXSectionStatA->Fit(fitTCMDecomposedLPi0,"QNRMEX0+","",0.4,2.);
    graphCombPi0InvXSectionTotA->Fit(fitTCMDecomposedHPi0,"QNRMEX0+","",4,50);
    fitTCMDecomposedHPi0->SetParameters(graphCombPi0InvXSectionTotA->GetY()[2],0.8, 2);    
    
    Double_t paramTCMPi0New[5]  = { 703678218.2483119965,0.5727060723,
                                    68561118949.5456314087,0.4481417931,3.0869940568};
    TF1* fitTCMInvXSectionPi0   = FitObject("tcm","fitTCMInvCrossSectionPi02760GeV","Pi0",graphCombPi0InvXSectionStatA,0.4,50. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
    
    fitTCMDecomposedLPi0->SetParameter(0, fitTCMInvXSectionPi0->GetParameter(0));
    fitTCMDecomposedLPi0->SetParameter(1, fitTCMInvXSectionPi0->GetParameter(1));
    fitTCMDecomposedHPi0->SetParameter(0, fitTCMInvXSectionPi0->GetParameter(2));
    fitTCMDecomposedHPi0->SetParameter(1, fitTCMInvXSectionPi0->GetParameter(3));
    fitTCMDecomposedHPi0->SetParameter(2, fitTCMInvXSectionPi0->GetParameter(4));
    
    // Tsallis fit 
    Double_t paramGraphPi0[3]                              = {1.0e12, 8., 0.13};
    TF1* fitInvXSectionPi0                       = FitObject("l","fitInvCrossSectionPi02760GeV","Pi0",histoEMCALPi0InvXSectionStat,0.4,50.,paramGraphPi0,"QNRMEX0+");
    TF1* fitInvXSectionPi0Graph                  = (TF1*)fitInvXSectionPi0->Clone("fitInvCrossSectionPi02760GeVGraph"); 

//     return;
    
    // *************************************************************************************************************
    // Shift graphs in X direction if desired
    // *************************************************************************************************************    
    if(bWCorrection.Contains("X")){
//         TF1* fitTsallisPi0PtMult                 = FitObject("tcmpt","TsallisMultWithPtPi02760GeV","Pi0");
        
//         fitTsallisPi0PtMult->SetParameters( fitTCMInvXSectionPi0->GetParameter(0),fitTCMInvXSectionPi0->GetParameter(1), fitTCMInvXSectionPi0->GetParameter(2), fitTCMInvXSectionPi0->GetParameter(3),
//                                             fitTCMInvXSectionPi0->GetParameter(4)); // standard parameter optimize if necessary
        TF1* fitTsallisPi0PtMult                 = FitObject("tmpt","TsallisMultWithPtPi02760GeV","Pi0");
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

        Size_t textSizeSpectra          = 0.04;
        TH1F * histoBinShift = new TH1F("histoBinShift","histoBinShift",1000,0.,40);
        SetStyleHistoTH1ForGraphs(histoBinShift, "#it{p}_{T} (GeV/#it{c})","bin shifted (X) / no shift",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.2);
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

        TLatex *labelRatioToFitBinShift   = new TLatex(0.94, 0.91, collisionSystem2760GeV.Data());
        SetStyleTLatex( labelRatioToFitBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitBinShift->Draw();
        TLatex *labelRatioToFitALICEBinShift    = new TLatex(0.94, 0.86, "ALICE");
        SetStyleTLatex( labelRatioToFitALICEBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitALICEBinShift->Draw();
        TLatex *labelRatioToFitPi0BinShift      = new TLatex(0.94, 0.807, "#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitPi0BinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitPi0BinShift->Draw();

        canvasShift->Update();
        canvasShift->SaveAs(Form("%s/BinShiftCorrection_Pi0.%s",outputDir.Data(),suffix.Data()));

        // *************************************************************************************************************
        // Plot control graphs
        // *************************************************************************************************************

        TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
        DrawGammaCanvasSettings( canvasDummy2,  0.1, 0.01, 0.015, 0.08);
        canvasDummy2->SetLogy();
        canvasDummy2->SetLogx();
        TH2F * histo2DDummy2;
        histo2DDummy2               = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.23,70.,1000,1e-1,10e11);
        SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
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
        canvasDummy2->Print(Form("%s/ComparisonShiftedPi0_2760GeV.%s",outputDir.Data(),suffix.Data()));
        delete canvasDummy2;
    }   

    // *************************************************************************************************************
    // redo fitting after binshifts
    // *************************************************************************************************************
    // Tsallis function    
    graphCombPi0InvXSectionTotA->Fit(fitInvXSectionPi0,"QNRMEX0+","",0.4,50.);
    fitInvXSectionPi0           = FitObject("l","fitInvCrossSectionPi02760GeV","Pi0",graphCombPi0InvXSectionTotA,0.4,40.,paramGraphPi0,"QNRMEX0+");
    fitInvXSectionPi0           = FitObject("l","fitInvCrossSectionPi02760GeV","Pi0",graphCombPi0InvXSectionTotA,0.4,40. ,paramGraphPi0,"QNRMEX0+");
    cout << WriteParameterToFile(fitInvXSectionPi0)<< endl;    
    fileFitsOutput << "Tsallis fit full pt range to tot err"<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitInvXSectionPi0)<< endl;    
    TF1* fitPi0Tsallis1           = FitObject("l","fitPi0Tsallis1","Pi0",graphCombPi0InvXSectionTotA,0.8,40. ,paramGraphPi0,"QNRMEX0+");
    TF1* fitPi0Tsallis2           = FitObject("l","fitPi0Tsallis2","Pi0",graphCombPi0InvXSectionTotA,1.5,40. ,paramGraphPi0,"QNRMEX0+");
    TF1* fitPi0Tsallis3           = FitObject("l","fitPi0Tsallis3","Pi0",graphCombPi0InvXSectionTotA,2.5,40. ,paramGraphPi0,"QNRMEX0+");
    TF1* fitPi0Tsallis4           = FitObject("l","fitPi0Tsallis4","Pi0",graphCombPi0InvXSectionTotA,0.4,10. ,paramGraphPi0,"QNRMEX0+");
    fileFitsOutput << "Tsallis fit 0.8-40 range to tot err" << endl;
    fileFitsOutput <<  WriteParameterToFile(fitPi0Tsallis1)<< endl;    
    fileFitsOutput << "Tsallis fit 1.5-40 range to tot err"<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPi0Tsallis2)<< endl;    
    fileFitsOutput << "Tsallis fit 2.5-40 range to tot err"<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPi0Tsallis3)<< endl;    
    fileFitsOutput << "Tsallis fit 0.4-10 range to tot err"<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPi0Tsallis4)<< endl;    
    
    
    
    //Two component model from Bylinkin
    fitTCMInvXSectionPi0        = FitObject("tcm","fitTCMInvCrossSectionPi02760GeV","Pi0",graphCombPi0InvXSectionTotA,0.4,40. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitTCMInvXSectionPi0)<< endl;
    fileFitsOutput << "TCM fit to tot err"<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitTCMInvXSectionPi0)<< endl;    
//     TF1* fitTCMInvXSectionPi02  = FitObject("tcm","fitTCMInvCrossSectionPi02760GeV2","Pi0",graphCombPi0InvXSectionTotA,0.4,40.,paramTCMPi0New,"QNRMEX0+","", kFALSE);
//     cout << WriteParameterToFile(fitTCMInvXSectionPi02)<< endl;

//     TF1* fitPowInvXSectionPi0   = FitObject("m","fitPowInvXSectionPi02760GeV","Pi0",graphCombPi0InvXSectionTotA,5,40. ,NULL,"QNRMEX0+","", kFALSE);
//     cout << WriteParameterToFile(fitPowInvXSectionPi0)<< endl;
    
    TF1* fitPCMTCMInvXSectionPi0    = FitObject("tcm","fitPCMTCMInvCrossSectionPi02760GeV","Pi0",graphPCMPi0InvXSectionStat,0.4,8. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPCMTCMInvXSectionPi0)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPCMTCMInvXSectionPi0)<< endl;    
//     return;
            
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
    TH1D* histoRatioPi0Pythia8ToFit                     = (TH1D*) histoPythia8InvXSectionPi0->Clone(); 
    histoRatioPi0Pythia8ToFit                           = CalculateHistoRatioToFit (histoRatioPi0Pythia8ToFit, fitTCMInvXSectionPi0); 
    histoRatioPi0Pythia8ToFit->GetXaxis()->SetRangeUser(0.4,40);
    TGraphErrors* graphRatioPi0Pythia8ToFit             = new TGraphErrors(histoRatioPi0Pythia8ToFit);
    while(graphRatioPi0Pythia8ToFit->GetX()[0] < 0.4) graphRatioPi0Pythia8ToFit->RemovePoint(0);
    while(graphRatioPi0Pythia8ToFit->GetX()[graphRatioPi0Pythia8ToFit->GetN()-1] > 40) graphRatioPi0Pythia8ToFit->RemovePoint(graphRatioPi0Pythia8ToFit->GetN()-1);
    
    TGraph* graphRatioPi0CombNLOMuHalf                  = (TGraph*)graphNLOCalcPi0MuHalf->Clone();
    TGraph* graphRatioPi0CombNLOMuOne                   = (TGraph*)graphNLOCalcPi0MuOne->Clone();
    TGraph* graphRatioPi0CombNLOMuTwo                   = (TGraph*)graphNLOCalcPi0MuTwo->Clone();
    TGraphAsymmErrors* graphRatioPi0CombNLODSS14        = (TGraphAsymmErrors*)graphNLODSS14Calc->Clone();
    graphRatioPi0CombNLOMuHalf                          = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuHalf, fitTCMInvXSectionPi0); 
    graphRatioPi0CombNLOMuOne                           = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuOne, fitTCMInvXSectionPi0); 
    graphRatioPi0CombNLOMuTwo                           = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuTwo, fitTCMInvXSectionPi0); 
    graphRatioPi0CombNLODSS14                           = CalculateGraphErrRatioToFit(graphRatioPi0CombNLODSS14, fitTCMInvXSectionPi0); 

    TGraph* graphRatioPi0CombCGC                        = (TGraph*)graphCGCPi0->Clone();
    graphRatioPi0CombCGC                                = CalculateGraphRatioToFit (graphRatioPi0CombCGC, fitTCMInvXSectionPi0); 
    
    TGraphAsymmErrors* graphRatioPi0CombCombFitTotA     = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotA->Clone();
    graphRatioPi0CombCombFitTotA                        = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitTotA, fitTCMInvXSectionPi0); 
    TGraphAsymmErrors* graphRatioPi0CombCombFitStatA    = (TGraphAsymmErrors*)graphCombPi0InvXSectionStatA->Clone();
    graphRatioPi0CombCombFitStatA                       = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitStatA, fitTCMInvXSectionPi0); 
    TGraphAsymmErrors* graphRatioPi0CombCombFitSysA     = (TGraphAsymmErrors*)graphCombPi0InvXSectionSysA->Clone();
    graphRatioPi0CombCombFitSysA                        = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitSysA, fitTCMInvXSectionPi0); 
    TGraphAsymmErrors* graphRatioPi0PCMCombFitStat      = (TGraphAsymmErrors*)graphPCMPi0InvXSectionStat->Clone();
    graphRatioPi0PCMCombFitStat                         = CalculateGraphErrRatioToFit(graphRatioPi0PCMCombFitStat, fitTCMInvXSectionPi0); 
    TGraphAsymmErrors* graphRatioPi0PCMCombFitSys       = (TGraphAsymmErrors*)graphPCMPi0InvXSectionSys->Clone();
    graphRatioPi0PCMCombFitSys                          = CalculateGraphErrRatioToFit(graphRatioPi0PCMCombFitSys, fitTCMInvXSectionPi0); 
    TGraphAsymmErrors* graphRatioPi0PHOSCombFitStat     = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionStat->Clone();
    graphRatioPi0PHOSCombFitStat                        = CalculateGraphErrRatioToFit(graphRatioPi0PHOSCombFitStat, fitTCMInvXSectionPi0); 
    TGraphAsymmErrors* graphRatioPi0PHOSCombFitSys      = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionSys->Clone();
    graphRatioPi0PHOSCombFitSys                         = CalculateGraphErrRatioToFit(graphRatioPi0PHOSCombFitSys, fitTCMInvXSectionPi0); 
    TGraphAsymmErrors* graphRatioPi0EMCALCombFitStat    = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionStat->Clone();
    graphRatioPi0EMCALCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioPi0EMCALCombFitStat, fitTCMInvXSectionPi0); 
    TGraphAsymmErrors* graphRatioPi0EMCALCombFitSys     = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionSys->Clone();
    graphRatioPi0EMCALCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioPi0EMCALCombFitSys, fitTCMInvXSectionPi0); 
    TGraphAsymmErrors* graphRatioPi0PCMEMCALCombFitStat = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionStat->Clone();
    graphRatioPi0PCMEMCALCombFitStat                    = CalculateGraphErrRatioToFit(graphRatioPi0PCMEMCALCombFitStat, fitTCMInvXSectionPi0); 
    TGraphAsymmErrors* graphRatioPi0PCMEMCALCombFitSys  = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionSys->Clone();
    graphRatioPi0PCMEMCALCombFitSys                     = CalculateGraphErrRatioToFit(graphRatioPi0PCMEMCALCombFitSys, fitTCMInvXSectionPi0); 
    TGraphAsymmErrors* graphRatioPi0EMCALMergedCombFitStat  = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionStat->Clone();
    graphRatioPi0EMCALMergedCombFitStat                     = CalculateGraphErrRatioToFit(graphRatioPi0EMCALMergedCombFitStat, fitTCMInvXSectionPi0); 
    TGraphAsymmErrors* graphRatioPi0EMCALMergedCombFitSys   = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionSys->Clone();
    graphRatioPi0EMCALMergedCombFitSys                      = CalculateGraphErrRatioToFit(graphRatioPi0EMCALMergedCombFitSys, fitTCMInvXSectionPi0); 
    

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
    histo2DPi0RatioToCombFit               = new TH2F("histo2DPi0RatioToCombFit","histo2DPi0RatioToCombFit",1000,0.23,70.,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DPi0RatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
    histo2DPi0RatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DPi0RatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
    histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.45,2.05);
    histo2DPi0RatioToCombFit->Draw("copy");

        ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatA);
 
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioPi0CombCombFitSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioPi0CombCombFitStatA->Draw("p,same,z");

        DrawGammaLines(0.23, 70. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.23, 70. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.9, 0.9,0.1, kGray, 7);

        TLatex *labelRatioToFitEnergy   = new TLatex(0.95, 0.92, collisionSystem2760GeV.Data());
        SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitEnergy->Draw();
        TLatex *labelRatioToFitALICE    = new TLatex(0.95, 0.87, "ALICE");
        SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitALICE->Draw();
        TLatex *labelRatioToFitPi0      = new TLatex(0.95, 0.817, "#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitPi0->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToCombFit_PP2760GeV.%s",outputDir.Data(),suffix.Data()));
    
    // **********************************************************************************************************************
    // *******************************************Plot Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************
    
    canvasRatioToCombFit->cd();
    histo2DPi0RatioToCombFit->Draw("copy");
    
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0PCMCombFitStat);
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0PHOSCombFitStat);
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0EMCALCombFitStat);
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0EMCALMergedCombFitStat);
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0PCMEMCALCombFitStat);
        
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

        DrawGammaLines(0.23, 70. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.23, 70. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.9, 0.9,0.5, kGray, 7);
        
        TLatex *labelRatioToFitEnergyInd   = new TLatex(0.115, 0.16, collisionSystem2760GeV.Data());
        SetStyleTLatex( labelRatioToFitEnergyInd, textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
        labelRatioToFitEnergyInd->Draw();
        TLatex *labelRatioToFitALICEInd    = new TLatex(0.95, 0.92, "ALICE");
        SetStyleTLatex( labelRatioToFitALICEInd, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitALICEInd->Draw();
        TLatex *labelRatioToFitPi0Ind      = new TLatex(0.95, 0.87, "#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitPi0Ind, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitPi0Ind->Draw();
    
        //****************************** Definition of the Legend ******************************************
        //**************** Row def ************************
        Double_t rowsLegendOnlyPi0Ratio[6]          = {0.92, 0.867, 0.807, 0.747, 0.687, 0.627};
//         Double_t rowsLegendOnlyPi0RatioAbs[6]       = {0.92, 2.18, 2.03, 1.87, 1.73, 1.56 };
        Double_t rowsLegendOnlyPi0RatioAbs[6]       = {0.92, 1.855, 1.75, 1.635, 1.53, 1.415 };
        
        Double_t columnsLegendOnlyPi0Ratio[3]       = {0.115, 0.325, 0.41};
        Double_t columnsLegendOnlyPi0RatioAbs[3]    = {0.115, 1.3, 1.95};
        Double_t lengthBox                          = 0.15;
        Double_t heightBox                          = 0.07/2;
        //****************** first Column **************************************************
        TLatex *textPCMOnlyRatioPi0                 = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[1],nameMeasGlobalLabel[0]);
        SetStyleTLatex( textPCMOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
        textPCMOnlyRatioPi0->Draw();
        TLatex *textPHOSOnlyRatioPi0                = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[2],nameMeasGlobalLabel[1]);
        SetStyleTLatex( textPHOSOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
        textPHOSOnlyRatioPi0->Draw();
        TLatex *textEMCALOnlyRatioPi0               = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[3],nameMeasGlobalLabel[2]);
        SetStyleTLatex( textEMCALOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
        textEMCALOnlyRatioPi0->Draw();
        TLatex *textPCMEMCALOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[4],nameMeasGlobalLabel[4]);
        SetStyleTLatex( textPCMEMCALOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
        textPCMEMCALOnlyRatioPi0->Draw();
        TLatex *textEMCALMergedOnlyRatioPi0         = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[5],nameMeasGlobalLabel[9]);
        SetStyleTLatex( textEMCALMergedOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
        textEMCALMergedOnlyRatioPi0->Draw();
        
        //****************** second Column *************************************************
        TLatex *textStatOnlyRatioPi0                = new TLatex(columnsLegendOnlyPi0Ratio[1],rowsLegendOnlyPi0Ratio[0] ,"stat");
        SetStyleTLatex( textStatOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
        textStatOnlyRatioPi0->Draw();
        TLatex *textSysOnlyRatioPi0                 = new TLatex(columnsLegendOnlyPi0Ratio[2] ,rowsLegendOnlyPi0Ratio[0],"syst");
        SetStyleTLatex( textSysOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
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

    // *******************************************************************************************************
    // ********************** Ratio to standalone PCM fit ****************************************************
    // *******************************************************************************************************
    TGraphAsymmErrors* graphRatioPi0PCMCombFitPCMStat      = (TGraphAsymmErrors*)graphPCMPi0InvXSectionStat->Clone();
    graphRatioPi0PCMCombFitPCMStat                         = CalculateGraphErrRatioToFit(graphRatioPi0PCMCombFitPCMStat, fitPCMTCMInvXSectionPi0); 
    TGraphAsymmErrors* graphRatioPi0PCMCombFitPCMSys       = (TGraphAsymmErrors*)graphPCMPi0InvXSectionSys->Clone();
    graphRatioPi0PCMCombFitPCMSys                          = CalculateGraphErrRatioToFit(graphRatioPi0PCMCombFitPCMSys, fitPCMTCMInvXSectionPi0); 
 
    canvasRatioToCombFit->cd();
    histo2DPi0RatioToCombFit->Draw("copy");
    
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitSys, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitStat, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitPCMSys, markerStyleDet[0]+4 ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitPCMStat, markerStyleDet[0]+4 ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
        
        graphRatioPi0PCMCombFitSys->Draw("E2same");
        graphRatioPi0PCMCombFitPCMSys->Draw("E2same");
        
        graphRatioPi0PCMCombFitStat->Draw("p,same,z");
        graphRatioPi0PCMCombFitPCMStat->Draw("p,same,z");

        DrawGammaLines(0.23, 70. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.23, 70. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.9, 0.9,0.5, kGray, 7);
        
        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitPi0->Draw();
    
        TLegend* legendPi0FitPCMStandalone        = GetAndSetLegend2(0.12, 0.92-(0.05*2), 0.45, 0.92, 38);
        legendPi0FitPCMStandalone->AddEntry(graphRatioPi0PCMCombFitSys,"PCM/comb fit");
        legendPi0FitPCMStandalone->AddEntry(graphRatioPi0PCMCombFitPCMSys,"PCM/PCM standalone fit");
        legendPi0FitPCMStandalone->Draw("");
        
    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfPCMToCombAndPCMStandaloneFit_PP.%s",outputDir.Data(),suffix.Data()));    
    
    
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
    Double_t xPtLimitsEtaWOMerged[50];
    Int_t maxNBinsEtaW0Merged       = GetBinning( xPtLimitsEtaWOMerged, "Eta", "2.76TeV", 4 );
    for (Int_t i = 0; i< maxNBinsEtaW0Merged; i++){
        cout << i << ": "<< xPtLimitsEtaWOMerged[i] <<" - " << xPtLimitsEtaWOMerged[i+1]<< ", " <<endl;
    }
    cout << endl;
    cout << "Nbins: eta "<< maxNBinsEtaW0Merged << endl;

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Eta                    
    Int_t offSetsEta[11]            = { 0,  0,  0,  0,  0,
                                        0,  0,  0,  0,  0, 
                                        0};
    Int_t offSetsEtaSys[11]         = { 1,  0,  4,  0,  2,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t offSetEtaShifting[11]     = { 0,  0,  3,  0,  1,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t nComBinsEtaShifting[11]   = { 6,  0,  11,  0,  10,
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
                                                                                                fileNameEtaOutputWeightingA,"2.76TeV", "Eta", kFALSE, 
                                                                                                NULL, fileNameCorrFactors
                                                                                            );
    graphCombEtaInvXSectionStatA->Print();
    if (graphCombEtaInvXSectionStatA == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }    
    if (graphCombEtaInvXSectionTotA->GetX()[0] == 0.25){
        graphCombEtaInvXSectionTotA->RemovePoint(0);
        cout << "removed first point from tot-graph" << endl; 
    }    
    if (graphCombEtaInvXSectionStatA->GetX()[0] == 0.25){
        graphCombEtaInvXSectionStatA->RemovePoint(0);
        cout << "removed first point from stat-graph" << endl; 
    }    
    if (graphCombEtaInvXSectionSysA->GetX()[0] == 0.25){
        graphCombEtaInvXSectionSysA->RemovePoint(0);
        cout << "removed first point from sys-graph" << endl; 
    }
//     return;
    
    // Reading weights from output file for plotting
    ifstream fileWeightsEtaReadA;
    fileWeightsEtaReadA.open(fileNameEtaOutputWeightingA,ios_base::in);
    cout << "reading" << fileNameEtaOutputWeightingA << endl;
    Double_t xValuesEtaReadA[50];
    Double_t weightsEtaReadA[11][50];
    Int_t availableEtaMeasA[11]        = { -1, -1, -1, -1, -1,
                                        -1, -1, -1, -1, -1,
                                        -1};
    Int_t nMeasSetEtaA                 = 3;
    Int_t nPtBinsEtaReadA              = 0;
    while(!fileWeightsEtaReadA.eof() && nPtBinsEtaReadA < 50){
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
    histo2DEtaWeights = new TH2F("histo2DEtaWeights","histo2DEtaWeights",11000,0.33, 27.,1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DEtaWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DEtaWeights->GetXaxis()->SetMoreLogLabels();
    histo2DEtaWeights->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DEtaWeights->GetYaxis()->SetRangeUser(-10,10);    
    histo2DEtaWeights->Draw("copy");

        TLegend* legendWeightsEta   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEtaA), 32);
        for (Int_t i = 0; i < nMeasSetEtaA; i++){
            DrawGammaSetMarkerTGraph(graphWeightsEtaA[availableEtaMeasA[i]], markerStyleDet[availableEtaMeasA[i]], markerSizeDet[availableEtaMeasA[i]]*0.5, colorDet[availableEtaMeasA[i]] , colorDet[availableEtaMeasA[i]]);
            graphWeightsEtaA[availableEtaMeasA[i]]->Draw("p,same,z");
            legendWeightsEta->AddEntry(graphWeightsEtaA[availableEtaMeasA[i]],nameMeasGlobalLabel[availableEtaMeasA[i]],"p");
        }    
        legendWeightsEta->Draw();

        labelWeightsEnergy->Draw();
        TLatex *labelWeightsEta         = new TLatex(0.7,0.16,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelWeightsEta, 0.85*textSizeLabelsPixel,4, 1, 43);
        labelWeightsEta->Draw();

//      DrawGammaLines(0.33, 27. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.33, 27. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.33, 27. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.33, 27. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.33, 27. , 0.2, 0.2,0.1, kGray, 3);
        
    canvasWeights->SaveAs(Form("%s/Eta_WeightsA.%s",outputDir.Data(),suffix.Data()));

//     return;
    // *********************************************************************************************************************
    // ************************************ Visualize relative errors Eta ******************************************************
    // *********************************************************************************************************************
    
    canvasRelSysErr->cd();
   
    TH2F * histo2DRelSysErrEta;
    histo2DRelSysErrEta                 = new TH2F("histo2DRelSysErrEta","histo2DRelSysErrEta",11000,0.33, 27.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErrEta, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErrEta->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErrEta->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelSysErrEta->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelSysErrEta->Draw("copy");

        TLegend* legendRelSysErrEta        = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetEtaA), 0.95, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetEtaA; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionEta[availableEtaMeasA[i]], markerStyleDet[availableEtaMeasA[i]], markerSizeDet[availableEtaMeasA[i]]*0.5, colorDet[availableEtaMeasA[i]],
                                     colorDet[availableEtaMeasA[i]]);
            sysErrorRelCollectionEta[availableEtaMeasA[i]]->Draw("p,same,z");
            legendRelSysErrEta->AddEntry(sysErrorRelCollectionEta[availableEtaMeasA[i]],nameMeasGlobalLabel[availableEtaMeasA[i]],"p");
        }    
        legendRelSysErrEta->Draw();

        labelRelSysErrEnergy->Draw();
        TLatex *labelRelSysErrEta       = new TLatex(0.15,0.85,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelSysErrEta, 0.85*textSizeLabelsPixel,4, 1, 43);
        labelRelSysErrEta->Draw();
        
    canvasRelSysErr->SaveAs(Form("%s/Eta_RelSysErr.%s",outputDir.Data(),suffix.Data()));
    
    histo2DRelSysErrEta->GetYaxis()->SetRangeUser(0,50.5);
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
    histo2DRelStatErrEta                = new TH2F("histo2DRelStatErrEta","histo2DRelStatErrEta",11000,0.33, 27.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErrEta, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErrEta->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErrEta->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelStatErrEta->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelStatErrEta->Draw("copy");
        
        TLegend* legendRelStatErrEta       = GetAndSetLegend2(0.14, 0.92-(0.035*nMeasSetEtaA), 0.45, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetEtaA; i++){
            DrawGammaSetMarker(statErrorRelCollectionEta[availableEtaMeasA[i]], markerStyleDet[availableEtaMeasA[i]], markerSizeDet[availableEtaMeasA[i]]*0.5, colorDet[availableEtaMeasA[i]] , 
                               colorDet[availableEtaMeasA[i]]);
            statErrorRelCollectionEta[availableEtaMeasA[i]]->Draw("p,same,z");
            legendRelStatErrEta->AddEntry(statErrorRelCollectionEta[availableEtaMeasA[i]],nameMeasGlobalLabel[availableEtaMeasA[i]],"p");
        }    
        legendRelStatErrEta->Draw();

        labelRelStatErrEnergy->Draw();
        TLatex *labelRelStatErrEta      = new TLatex(0.75,0.85,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelStatErrEta, 0.85*textSizeLabelsPixel,4, 1, 43);
        labelRelStatErrEta->Draw();
        
    canvasRelStatErr->SaveAs(Form("%s/Eta_RelStatErr.%s",outputDir.Data(),suffix.Data()));

    histo2DRelStatErrEta->GetYaxis()->SetRangeUser(0,50.5);
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
    histo2DRelTotErrEta                 = new TH2F("histo2DRelTotErrEta","histo2DRelTotErrEta",11000,0.33, 27.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrEta, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrEta->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrEta->GetXaxis()->SetLabelOffset(-0.01);    
    histo2DRelTotErrEta->Draw("copy");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionRelTotA, markerStyleComb+4, markerSizeComb, kBlue+2 , kBlue+2);
        graphCombEtaInvXSectionRelTotA->Draw("p,same,z");

        TLegend* legendRelTotErrEta     = GetAndSetLegend2(0.14, 0.92-(0.035*1), 0.45, 0.92, 32);
        legendRelTotErrEta->AddEntry(graphCombEtaInvXSectionRelTotA,"All","p");
        legendRelTotErrEta->Draw();

        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrEta       = new TLatex(0.75,0.85,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelTotErrEta, 0.85*textSizeLabelsPixel,4, 1, 43);
        labelRelTotErrEta->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/Eta_RelTotErr.%s",outputDir.Data(),suffix.Data()));
        
    histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,50.5);
    histo2DRelTotErrEta->Draw("copy");
        graphCombEtaInvXSectionRelTotA->Draw("p,same,z");
    
        legendRelTotErr->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrEta->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/Eta_RelTotErrZoomed.%s",outputDir.Data(),suffix.Data()));
    
    histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,50.5);
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
    Double_t paramGraphEta[3]                    = {1.0e10, 7., 0.13};
    TF1* fitInvXSectionEta                       = FitObject("l","fitInvCrossSectionEta2760GeV","Eta",graphCombEtaInvXSectionTotAUnshi,0.5,20.,paramGraphEta,"QNRMEX0+");
    TF1* fitInvXSectionEtaGraph                  = (TF1*)fitInvXSectionEta->Clone("fitInvCrossSectionEta2760GeVGraph"); 
    
    Double_t paramTCMEta[5]  = {graphCombEtaInvXSectionTotA->GetY()[0],0.5,graphCombEtaInvXSectionTotA->GetY()[0]/1000,0.8,3};


    // *************************************************************************************************************
    // Shift graphs in X direction if desired
    // *************************************************************************************************************        
    if(bWCorrection.Contains("X") ){
        TF1* fitTsallisEtaPtMult                 = FitObject("tmpt","TsallisMultWithPtEta2760GeV","Eta");
        fitTsallisEtaPtMult->SetParameters(paramGraphEta[0],paramGraphEta[1], paramGraphEta[2]) ; // standard parameter optimize if necessary
//         fitTsallisEtaPtMult->SetParameters(paramTCMEta[0],paramTCMEta[1], paramTCMEta[2], paramTCMEta[3], paramTCMEta[4]) ; // standard parameter optimize if necessary

        TGraphAsymmErrors* graphCombEtaInvXSectionTotANoShift = (TGraphAsymmErrors*) graphCombEtaInvXSectionTotA->Clone("Eta_NoShift");

        graphCombEtaInvXSectionTotA              = ApplyXshift(graphCombEtaInvXSectionTotA, fitTsallisEtaPtMult);
        
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

        //***************************************************************************************************************
        //************************************Plotting binshift corrections *********************************************
        //***************************************************************************************************************

        TCanvas* canvasShift = new TCanvas("canvasShift","",0,0,1000,900);// gives the page size
        DrawGammaCanvasSettings( canvasShift, 0.10, 0.017, 0.015, 0.08);

        Size_t textSizeSpectra          = 0.04;
        TH1F * histoBinShift = new TH1F("histoBinShift","histoBinShift",1000,0., 20.);
        SetStyleHistoTH1ForGraphs(histoBinShift, "#it{p}_{T} (GeV/#it{c})","bin shifted (X) / no shift",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.2);
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

        TLatex *labelRatioToFitBinShift   = new TLatex(0.94, 0.91, collisionSystem2760GeV.Data());
        SetStyleTLatex( labelRatioToFitBinShift, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
        labelRatioToFitBinShift->Draw();
        TLatex *labelRatioToFitALICEBinShift    = new TLatex(0.94, 0.86, "ALICE");
        SetStyleTLatex( labelRatioToFitALICEBinShift, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
        labelRatioToFitALICEBinShift->Draw();
        TLatex *labelRatioToFitPi0BinShift      = new TLatex(0.94, 0.807, "#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitPi0BinShift, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);;
        labelRatioToFitPi0BinShift->Draw();

        canvasShift->Update();
        canvasShift->SaveAs(Form("%s/BinShiftCorrection_Eta.%s",outputDir.Data(),suffix.Data()));

        // *************************************************************************************************************
        // Plot control graphs
        // *************************************************************************************************************

        TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
        DrawGammaCanvasSettings( canvasDummy2,  0.1, 0.01, 0.015, 0.08);
        canvasDummy2->SetLogy();
        canvasDummy2->SetLogx();
        TH2F* histo2DDummy3;
        histo2DDummy3               = new TH2F("histo2DDummy3","histo2DDummy3",1000,0.23,25.,1000,1e1,10e11);
        SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
        histo2DDummy3->DrawCopy(); 
    
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStatAUnshi, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
        graphCombEtaInvXSectionStatAUnshi->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionTotA, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombEtaInvXSectionTotA->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvXSectionSys, markerStyleDet[0] ,markerSizeDet[0]/2, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        graphPCMEtaInvXSectionSys->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvXSectionStat, markerStyleDet[0] ,markerSizeDet[0]/2, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        graphPCMEtaInvXSectionStat->Draw("psame");
        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaInvXSectionSys, markerStyleDet[2] ,markerSizeDet[2]/2, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        graphEMCALEtaInvXSectionSys->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaInvXSectionSys, markerStyleDet[4] ,markerSizeDet[4]/2, colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
        graphPCMEMCALEtaInvXSectionSys->Draw("pEsame");

        fitInvXSectionEta->SetLineColor(kBlue+2);
        fitInvXSectionEta->Draw("same");
        
//         canvasDummy2->Update();
        canvasDummy2->Print(Form("%s/ComparisonShiftedEta_2760GeV.%s",outputDir.Data(),suffix.Data()));
        delete canvasDummy2;
        delete histo2DDummy3;
    }
    
//     return;
    
    // *************************************************************************************************************
    // redo fitting after binshifts
    // *************************************************************************************************************
    // Tsallis function
    graphCombEtaInvXSectionTotA->Fit(fitInvXSectionEta,"QNRMEX0+","",0.5,20.);
    fitInvXSectionEta        = FitObject("l","fitInvCrossSectionEta2760GeV","Eta",graphCombEtaInvXSectionTotA,0.5,20.,paramGraphEta,"QNRMEX0+");
    fitInvXSectionEta        = FitObject("l","fitInvCrossSectionEta2760GeV","Eta",graphCombEtaInvXSectionTotA,0.5,20. ,paramGraphEta,"QNRMEX0+");
    cout << WriteParameterToFile(fitInvXSectionEta)<< endl;
    fileFitsOutput << "Tsallis eta meson full tot err" << endl;
    fileFitsOutput <<  WriteParameterToFile(fitInvXSectionEta)<< endl;    
    
    // Two component model by Bylinkin
    TF1* fitTCMDecomposedLEta   = FitObject("tcmlow","twoCompModelEta_DecL", "Eta", NULL, 0.5, 20.);
    TF1* fitTCMDecomposedHEta   = FitObject("tcmhigh","twoCompModelEta_DecH", "Eta", NULL, 0.5, 20.);
//     fitTCMDecomposedHEta->SetParameters(graphCombEtaInvXSectionTotA->GetY()[2],0.8, 2);    
//     Double_t paramTCMEtaNew[5]  = { fitTCMDecomposedLEta->GetParameter(0),fitTCMDecomposedLEta->GetParameter(1),
//                                     fitTCMDecomposedHEta->GetParameter(0),fitTCMDecomposedHEta->GetParameter(1),fitTCMDecomposedHEta->GetParameter(2)};
    
    TF1* fitTCMInvXSectionEta= FitObject("tcm","fitTCMInvCrossSectionEta2760GeV","Eta",graphCombEtaInvXSectionTotA,0.5,20.,paramTCMEta,"QNRMEX0+","", kFALSE);
    fitTCMInvXSectionEta     = FitObject("tcm","fitTCMInvCrossSectionEta2760GeV","Eta",graphCombEtaInvXSectionTotA,0.5,20. ,paramTCMEta,"QNRMEX0+","", kFALSE);
    fitTCMDecomposedLEta->SetParameter(0, fitTCMInvXSectionEta->GetParameter(0));
    fitTCMDecomposedLEta->SetParameter(1, fitTCMInvXSectionEta->GetParameter(1));
    fitTCMDecomposedHEta->SetParameter(0, fitTCMInvXSectionEta->GetParameter(2));
    fitTCMDecomposedHEta->FixParameter(1, fitTCMInvXSectionEta->GetParameter(3));
    fitTCMDecomposedHEta->FixParameter(2, fitTCMInvXSectionEta->GetParameter(4));
    graphCombEtaInvXSectionTotA->Fit(fitTCMDecomposedHEta,"QNRMEX0+","",3,20);
    
    cout << WriteParameterToFile(fitTCMInvXSectionEta)<< endl;
    fileFitsOutput << "TCM eta meson full tot err" << endl;
    fileFitsOutput <<  WriteParameterToFile(fitTCMInvXSectionEta)<< endl;    
    
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
    DrawGammaCanvasSettings( canvasDummy2,  0.1, 0.01, 0.015, 0.08);
    canvasDummy2->SetLogy();
    canvasDummy2->SetLogx();
    TH2F* histo2DDummy3;
    histo2DDummy3               = new TH2F("histo2DDummy3","histo2DDummy3",1000,0.23,25.,1000,1e1,10e11);
    SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
    histo2DDummy3->DrawCopy(); 

    DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionTotA, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
    graphCombEtaInvXSectionTotA->Draw("pEsame");

    fitInvXSectionEta->SetLineColor(kBlue+2);
    fitInvXSectionEta->Draw("same");
    fitTCMInvXSectionEta->SetLineColor(kRed+2);
    fitTCMInvXSectionEta->Draw("same");
    
    canvasDummy2->Update();
    canvasDummy2->Print(Form("%s/ComparisonWithFitEta_2760GeV.%s",outputDir.Data(),suffix.Data()));
    delete canvasDummy2;
    delete histo2DDummy3;

//     return;

    // *************************************************************************************************************
    // Calculate ratios to combined fit
    // *************************************************************************************************************    
    TH1D* histoRatioEtaPythia8ToFit                     = (TH1D*) histoPythia8InvXSectionEta->Clone(); 
    histoRatioEtaPythia8ToFit                           = CalculateHistoRatioToFit (histoRatioEtaPythia8ToFit, fitTCMInvXSectionEta); 
    histoRatioEtaPythia8ToFit->GetXaxis()->SetRangeUser(0.6,25);
    TGraphErrors* graphRatioEtaPythia8ToFit             = new TGraphErrors(histoRatioEtaPythia8ToFit); 
    while(graphRatioEtaPythia8ToFit->GetX()[0] < 0.6) graphRatioEtaPythia8ToFit->RemovePoint(0);
    while(graphRatioEtaPythia8ToFit->GetX()[graphRatioEtaPythia8ToFit->GetN()-1] > 25) graphRatioEtaPythia8ToFit->RemovePoint(graphRatioEtaPythia8ToFit->GetN()-1);

    
    TGraph* graphRatioEtaCombNLOMuHalf                  = (TGraph*)graphNLOCalcEtaMuHalf->Clone();
    TGraph* graphRatioEtaCombNLOMuOne                   = (TGraph*)graphNLOCalcEtaMuOne->Clone();
    TGraph* graphRatioEtaCombNLOMuTwo                   = (TGraph*)graphNLOCalcEtaMuTwo->Clone();
    graphRatioEtaCombNLOMuHalf                          = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuHalf, fitTCMInvXSectionEta); 
    graphRatioEtaCombNLOMuOne                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuOne, fitTCMInvXSectionEta); 
    graphRatioEtaCombNLOMuTwo                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuTwo, fitTCMInvXSectionEta); 
    
    TGraphAsymmErrors* graphRatioEtaCombCombFitTotA     = (TGraphAsymmErrors*)graphCombEtaInvXSectionTotA->Clone();
    graphRatioEtaCombCombFitTotA                        = CalculateGraphErrRatioToFit(graphRatioEtaCombCombFitTotA, fitTCMInvXSectionEta); 
    TGraphAsymmErrors* graphRatioEtaCombCombFitStatA    = (TGraphAsymmErrors*)graphCombEtaInvXSectionStatA->Clone();
    graphRatioEtaCombCombFitStatA                       = CalculateGraphErrRatioToFit(graphRatioEtaCombCombFitStatA, fitTCMInvXSectionEta); 
    TGraphAsymmErrors* graphRatioEtaCombCombFitSysA     = (TGraphAsymmErrors*)graphCombEtaInvXSectionSysA->Clone();

    graphRatioEtaCombCombFitSysA                        = CalculateGraphErrRatioToFit(graphRatioEtaCombCombFitSysA, fitTCMInvXSectionEta); 
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
    

    // **********************************************************************************************************************
    // ******************************************* Ratio of Comb to Fit ****************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 54;
    canvasRatioToCombFit->cd();
    TH2F * histo2DEtaRatioToCombFit;
    histo2DEtaRatioToCombFit               = new TH2F("histo2DEtaRatioToCombFit","histo2DEtaRatioToCombFit",1000,0.33, 27.,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DEtaRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
    histo2DEtaRatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DEtaRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
    histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.45,2.05);
    histo2DEtaRatioToCombFit->Draw("copy");

        ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombCombFitStatA);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSysA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioEtaCombCombFitSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatA, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioEtaCombCombFitStatA->Draw("p,same,z");

        DrawGammaLines(0.33, 27. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.33, 27. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.33, 27. , 0.9, 0.9,0.1, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        TLatex *labelRatioToFitEta      = new TLatex(0.852,0.82,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitEta, textSizeLabelsPixel,4, 1, 43);
        labelRatioToFitEta->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfCombToCombFit_PP2760GeV.%s",outputDir.Data(),suffix.Data()));
    
    // **********************************************************************************************************************
    // ******************************************* Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************
    
    canvasRatioToCombFit->cd();
    histo2DEtaRatioToCombFit->Draw("copy");
        ProduceGraphAsymmWithoutXErrors(graphRatioEtaPCMCombFitStat);
        ProduceGraphAsymmWithoutXErrors(graphRatioEtaEMCALCombFitStat);
        ProduceGraphAsymmWithoutXErrors(graphRatioEtaPCMEMCALCombFitStat);
        
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaPCMCombFitSys, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaPCMCombFitStat, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaEMCALCombFitSys, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaEMCALCombFitStat, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaPCMEMCALCombFitSys, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaPCMEMCALCombFitStat, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4]);
        
        graphRatioEtaPCMCombFitSys->Draw("E2same");
        graphRatioEtaEMCALCombFitSys->Draw("E2same");
        graphRatioEtaPCMEMCALCombFitSys->Draw("E2same");
        
        graphRatioEtaPCMCombFitStat->Draw("p,same,z");
        graphRatioEtaEMCALCombFitStat->Draw("p,same,z");
        graphRatioEtaPCMEMCALCombFitStat->Draw("p,same,z");

        DrawGammaLines(0.33, 27. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.33, 27. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.33, 27. , 0.9, 0.9,0.5, kGray, 7);
        
        labelRatioToFitEnergyInd->Draw();
        TLatex *labelRatioToFitALICEInd2   = new TLatex(0.95, 0.21, "ALICE");
        SetStyleTLatex( labelRatioToFitALICEInd2, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitALICEInd2->Draw();
        TLatex *labelRatioToFitEtaInd      = new TLatex(0.95,0.16,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitEtaInd, textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
        labelRatioToFitEtaInd->Draw();

        
        //****************************** Definition of the Legend ******************************************
        Double_t rowsLegendOnlyEtaRatio[4]          = {0.91, 0.85, 0.79, 0.76};
//         Double_t rowsLegendOnlyEtaRatioAbs[4]       = {0.91, 2.245, 2.11, 1.975};
        Double_t rowsLegendOnlyEtaRatioAbs[4]       = {0.91, 1.83, 1.71, 1.975};
        Double_t columnsLegendOnlyEtaRatio[6]       = {0.115, 0.21, 0.285, 0.39, 0.58, 0.655};
//         Double_t columnsLegendOnlyEtaRatioAbs[3]    = {0.115, 1.2, 1.62};
        Double_t columnsLegendOnlyEtaRatioAbs[6]    = {0.115, 0.72, 0.98, 0.39, 4.3, 5.8};
        
        //****************** first Column **************************************************
        TLatex *textPCMOnlyRatioEta                 = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[1],nameMeasGlobalLabel[0]);
        SetStyleTLatex( textPCMOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
        textPCMOnlyRatioEta->Draw();
        TLatex *textEMCALOnlyRatioEta               = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[2],nameMeasGlobalLabel[2]);
        SetStyleTLatex( textEMCALOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
        textEMCALOnlyRatioEta->Draw();
        TLatex *textPCMEMCALOnlyRatioEta            = new TLatex(columnsLegendOnlyEtaRatio[3],rowsLegendOnlyEtaRatio[1],nameMeasGlobalLabel[4]);
        SetStyleTLatex( textPCMEMCALOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
        textPCMEMCALOnlyRatioEta->Draw();
        
        //****************** second Column *************************************************
        TLatex *textStatOnlyRatioEta                = new TLatex(columnsLegendOnlyEtaRatio[1],rowsLegendOnlyEtaRatio[0] ,"stat");
        SetStyleTLatex( textStatOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
        textStatOnlyRatioEta->Draw();
        TLatex *textSysOnlyRatioEta                 = new TLatex(columnsLegendOnlyEtaRatio[2] ,rowsLegendOnlyEtaRatio[0],"syst");
        SetStyleTLatex( textSysOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
        textSysOnlyRatioEta->Draw();
        TLatex *textStatOnlyRatioEta2                = new TLatex(columnsLegendOnlyEtaRatio[4],rowsLegendOnlyEtaRatio[0] ,"stat");
        SetStyleTLatex( textStatOnlyRatioEta2, textSizeLabelsPixel,4, 1, 43);
        textStatOnlyRatioEta2->Draw();
        TLatex *textSysOnlyRatioEta2                 = new TLatex(columnsLegendOnlyEtaRatio[5] ,rowsLegendOnlyEtaRatio[0],"syst");
        SetStyleTLatex( textSysOnlyRatioEta2, textSizeLabelsPixel,4, 1, 43);
        textSysOnlyRatioEta2->Draw();
        TMarker* markerPCMEtaOnlyRatioEta           = CreateMarkerFromGraph(graphRatioEtaPCMCombFitSys,columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[1],1);
        markerPCMEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[1]);
        TMarker* markerEMCALEtaOnlyRatioEta         = CreateMarkerFromGraph(graphRatioEtaEMCALCombFitSys, columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[2],1);
        markerEMCALEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[2]);
        TMarker* markerPCMEMCALEtaOnlyRatioEta      = CreateMarkerFromGraph(graphRatioEtaPCMEMCALCombFitSys, columnsLegendOnlyEtaRatio[4] ,rowsLegendOnlyEtaRatio[1],1);
        markerPCMEMCALEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[4] ,rowsLegendOnlyEtaRatioAbs[1]);

        TBox* boxPCMEtaOnlyRatioEta                 = CreateBoxFromGraph(graphRatioEtaPCMCombFitSys, columnsLegendOnlyEtaRatioAbs[2]-0.5*lengthBox , rowsLegendOnlyEtaRatioAbs[1]- heightBox,
                                                                         columnsLegendOnlyEtaRatioAbs[2]+ 1.5*lengthBox, rowsLegendOnlyEtaRatioAbs[1]+ heightBox);
        boxPCMEtaOnlyRatioEta->Draw("l");
        TBox* boxEMCALEtaOnlyRatioEta               = CreateBoxFromGraph(graphRatioEtaEMCALCombFitSys, columnsLegendOnlyEtaRatioAbs[2]-0.5*lengthBox , rowsLegendOnlyEtaRatioAbs[2]- heightBox,
                                                                         columnsLegendOnlyEtaRatioAbs[2]+ 1.5*lengthBox, rowsLegendOnlyEtaRatioAbs[2]+ heightBox);
        boxEMCALEtaOnlyRatioEta->Draw("l");
        TBox* boxPCMEMCALEtaOnlyRatioEta            = CreateBoxFromGraph(graphRatioEtaPCMEMCALCombFitSys, columnsLegendOnlyEtaRatioAbs[5]-2.7*lengthBox , rowsLegendOnlyEtaRatioAbs[1]- heightBox,
                                                                         columnsLegendOnlyEtaRatioAbs[5]+ 9.5*lengthBox, rowsLegendOnlyEtaRatioAbs[1]+ heightBox);
        boxPCMEMCALEtaOnlyRatioEta->Draw("l");
        
    canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfIndividualMeasToCombFit_PP.%s",outputDir.Data(),suffix.Data()));
    

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
    Int_t offSetsEtaToPi0[11]           = { 0,  0,  0,  0,  0,
                                            0,  0,  0,  0,  0, 
                                            0};
    Int_t offSetsEtaToPi0Sys[11]        = { 1,  0,  4,  0,  2,
                                            0,  0,  0,  0,  0,
                                            0};
    Int_t offSetEtaToPi0Shifting[11]    = { 0,  0,  3,  0,  1,
                                            0,  0,  0,  0,  0,
                                            0};
    Int_t nComBinsEtaToPi0Shifting[11]  = { 6,  0,  11,  0,  10,
                                            0,  0,  0,  0,  0,
                                            0};
                                        
    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************
                                
    TGraph* graphWeightsEtaToPi0A[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsEtaToPi0A[i]                    = NULL;
    }    
                                
    // Declaration & calculation of combined spectrum                            
    TString fileNameEtaToPi0OutputWeightingA        = Form("%s/EtaToPi0_WeightingMethodA.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombEtaToPi0StatA       = NULL;
    TGraphAsymmErrors* graphCombEtaToPi0SysA        = NULL;
    TGraphAsymmErrors* graphCombEtaToPi0TotA        = CombinePtPointsSpectraFullCorrMat(   statErrorCollectionEtaToPi0,    sysErrorCollectionEtaToPi0,     
                                                                                           xPtLimitsEtaWOMerged, maxNBinsEtaW0Merged,
                                                                                           offSetsEtaToPi0, offSetsEtaToPi0Sys,
                                                                                           graphCombEtaToPi0StatA, graphCombEtaToPi0SysA,
                                                                                           fileNameEtaToPi0OutputWeightingA,"2.76TeV", "EtaToPi0", kFALSE,
                                                                                           NULL, fileNameCorrFactors
                                                                                        );

    graphCombEtaToPi0StatA->Print();
    if (graphCombEtaToPi0StatA == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }    

//     return;
    
    // Reading weights from output file for plotting
    ifstream fileWeightsEtaToPi0ReadA;
    fileWeightsEtaToPi0ReadA.open(fileNameEtaToPi0OutputWeightingA,ios_base::in);
    cout << "reading" << fileNameEtaToPi0OutputWeightingA << endl;
    Double_t xValuesEtaToPi0ReadA[50];
    Double_t weightsEtaToPi0ReadA[11][50];
    Int_t availableEtaToPi0MeasA[11]        = { -1, -1, -1, -1, -1,
                                        -1, -1, -1, -1, -1,
                                        -1};
    Int_t nMeasSetEtaToPi0A                 = 3;
    Int_t nPtBinsEtaToPi0ReadA              = 0;
    while(!fileWeightsEtaToPi0ReadA.eof() && nPtBinsEtaToPi0ReadA < 50){
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
    histo2DEtaToPi0Weights = new TH2F("histo2DEtaToPi0Weights","histo2DEtaToPi0Weights",11000,0.33, 27.,1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DEtaToPi0Weights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DEtaToPi0Weights->GetXaxis()->SetMoreLogLabels();
    histo2DEtaToPi0Weights->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DEtaToPi0Weights->GetYaxis()->SetRangeUser(-10,10);    
    histo2DEtaToPi0Weights->Draw("copy");

        TLegend* legendWeightsEtaToPi0   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEtaToPi0A), 32);
        for (Int_t i = 0; i < nMeasSetEtaToPi0A; i++){
            DrawGammaSetMarkerTGraph(graphWeightsEtaToPi0A[availableEtaToPi0MeasA[i]], markerStyleDet[availableEtaToPi0MeasA[i]], markerSizeDet[availableEtaToPi0MeasA[i]]*0.5, colorDet[availableEtaToPi0MeasA[i]] , colorDet[availableEtaToPi0MeasA[i]]);
            graphWeightsEtaToPi0A[availableEtaToPi0MeasA[i]]->Draw("p,same,z");
            legendWeightsEtaToPi0->AddEntry(graphWeightsEtaToPi0A[availableEtaToPi0MeasA[i]],nameMeasGlobalLabel[availableEtaToPi0MeasA[i]],"p");
        }    
        legendWeightsEtaToPi0->Draw();

        labelWeightsEnergy->Draw();
        TLatex *labelWeightsEtaToPi0         = new TLatex(0.7,0.16,"#eta/#pi^{0}");
        SetStyleTLatex( labelWeightsEtaToPi0, 0.85*textSizeLabelsPixel,4, 1, 43);
        labelWeightsEtaToPi0->Draw();

//      DrawGammaLines(0.33, 27. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.33, 27. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.33, 27. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.33, 27. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.33, 27. , 0.2, 0.2,0.1, kGray, 3);
        
    canvasWeights->SaveAs(Form("%s/EtaToPi0_WeightsA.%s",outputDir.Data(),suffix.Data()));


    // *********************************************************************************************************************
    // ************************************ Visualize relative errors EtaToPi0 ******************************************************
    // *********************************************************************************************************************
    
    canvasRelSysErr->cd();
   
    TH2F * histo2DRelSysErrEtaToPi0;
    histo2DRelSysErrEtaToPi0                 = new TH2F("histo2DRelSysErrEtaToPi0","histo2DRelSysErrEtaToPi0",11000,0.33, 27.,1000,0,60.5);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErrEtaToPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErrEtaToPi0->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelSysErrEtaToPi0->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelSysErrEtaToPi0->Draw("copy");

        TLegend* legendRelSysErrEtaToPi0        = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetEtaToPi0A), 0.95, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetEtaToPi0A; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionEtaToPi0[availableEtaToPi0MeasA[i]], markerStyleDet[availableEtaToPi0MeasA[i]], markerSizeDet[availableEtaToPi0MeasA[i]]*0.5, 
                                     colorDet[availableEtaToPi0MeasA[i]], colorDet[availableEtaToPi0MeasA[i]]);
            sysErrorRelCollectionEtaToPi0[availableEtaToPi0MeasA[i]]->Draw("p,same,z");
            legendRelSysErrEtaToPi0->AddEntry(sysErrorRelCollectionEtaToPi0[availableEtaToPi0MeasA[i]],nameMeasGlobalLabel[availableEtaToPi0MeasA[i]],"p");
        }    
        legendRelSysErrEtaToPi0->Draw();

        labelRelSysErrEnergy->Draw();
        TLatex *labelRelSysErrEtaToPi0       = new TLatex(0.15,0.85,"#eta/#pi^{0}");
        SetStyleTLatex( labelRelSysErrEtaToPi0, 0.85*textSizeLabelsPixel,4, 1, 43);
        labelRelSysErrEtaToPi0->Draw();
        
    canvasRelSysErr->SaveAs(Form("%s/EtaToPi0_RelSysErr.%s",outputDir.Data(),suffix.Data()));
        
    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors EtaToPi0 **************************************************
    //  *********************************************************************************************************************
    
    canvasRelStatErr->cd();
   
    TH2F * histo2DRelStatErrEtaToPi0;
    histo2DRelStatErrEtaToPi0                = new TH2F("histo2DRelStatErrEtaToPi0","histo2DRelStatErrEtaToPi0",11000,0.33, 27.,1000,0,60.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErrEtaToPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErrEtaToPi0->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelStatErrEtaToPi0->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelStatErrEtaToPi0->Draw("copy");
        
        TLegend* legendRelStatErrEtaToPi0       = GetAndSetLegend2(0.14, 0.92-(0.035*nMeasSetEtaToPi0A), 0.45, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetEtaToPi0A; i++){
            DrawGammaSetMarker(statErrorRelCollectionEtaToPi0[availableEtaToPi0MeasA[i]], markerStyleDet[availableEtaToPi0MeasA[i]], markerSizeDet[availableEtaToPi0MeasA[i]]*0.5,
                               colorDet[availableEtaToPi0MeasA[i]] , colorDet[availableEtaToPi0MeasA[i]]);
            statErrorRelCollectionEtaToPi0[availableEtaToPi0MeasA[i]]->Draw("p,same,z");
            legendRelStatErrEtaToPi0->AddEntry(statErrorRelCollectionEtaToPi0[availableEtaToPi0MeasA[i]],nameMeasGlobalLabel[availableEtaToPi0MeasA[i]],"p");
        }    
        legendRelStatErrEtaToPi0->Draw();

        labelRelStatErrEnergy->Draw();
        TLatex *labelRelStatErrEtaToPi0      = new TLatex(0.75,0.85,"#eta/#pi^{0}");
        SetStyleTLatex( labelRelStatErrEtaToPi0, 0.85*textSizeLabelsPixel,4, 1, 43);
        labelRelStatErrEtaToPi0->Draw();
        
    canvasRelStatErr->SaveAs(Form("%s/EtaToPi0_RelStatErr.%s",outputDir.Data(),suffix.Data()));


    //  *********************************************************************************************************************
    //  ************************ Visualize relative total errors of different combination methods Eta ***********************
    //  *********************************************************************************************************************
    TGraphAsymmErrors* graphCombEtaToPi0RelStatA       = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0StatA, "relativeStatErrorEtaToPi0_MethodA");
    TGraphAsymmErrors* graphCombEtaToPi0RelSysA        = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0SysA, "relativeSysErrorEtaToPi0_MethodA");
    TGraphAsymmErrors* graphCombEtaToPi0RelTotA        = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0TotA, "relativeTotalErrorEtaToPi0_MethodA");
    
    while(graphCombEtaToPi0RelStatA->GetX()[0] < 0.6) graphCombEtaToPi0RelStatA->RemovePoint(0);
    while(graphCombEtaToPi0RelSysA->GetX()[0] < 0.6) graphCombEtaToPi0RelSysA->RemovePoint(0);
    while(graphCombEtaToPi0RelTotA->GetX()[0] < 0.6) graphCombEtaToPi0RelTotA->RemovePoint(0);

    canvasRelTotErr->cd();
    TH2F * histo2DRelTotErrEtaToPi0;
    histo2DRelTotErrEtaToPi0                 = new TH2F("histo2DRelTotErrEtaToPi0","histo2DRelTotErrEtaToPi0",11000,0.33, 27.,1000,0,60.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrEtaToPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrEtaToPi0->GetXaxis()->SetLabelOffset(-0.01);    
    histo2DRelTotErrEtaToPi0->Draw("copy");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelTotA, markerStyleComb+4, markerSizeComb, kBlue+2 , kBlue+2);
        graphCombEtaToPi0RelTotA->Draw("p,same,z");

        TLegend* legendRelTotErrEtaToPi0     = GetAndSetLegend2(0.14, 0.92-(0.035*1), 0.45, 0.92, 32);
        legendRelTotErrEtaToPi0->AddEntry(graphCombEtaToPi0RelTotA,"All","p");
        legendRelTotErrEtaToPi0->Draw();

        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrEtaToPi0       = new TLatex(0.75,0.85,"#eta/#pi^{0}");
        SetStyleTLatex( labelRelTotErrEtaToPi0, 0.85*textSizeLabelsPixel,4, 1, 43);
        labelRelTotErrEtaToPi0->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/EtaToPi0_RelTotErr.%s",outputDir.Data(),suffix.Data()));
            
    histo2DRelTotErrEtaToPi0->GetYaxis()->SetRangeUser(0,50.5);
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
    // ******************************************* Mass and width for pi0 at 2.76TeV ****************************************
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
        
        TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, 0.23, 25. ,1000., -30, 40);
        SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                                  0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
        histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,28.5);
        histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
        histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
        histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
        histo2DAllPi0FWHM->DrawCopy(); 

        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0FWHM, markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
        graphEMCALPi0FWHM->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0FWHMMC, markerStyleDetMC[2], markerSizeDetMC[2]*0.55, colorDetMC[2] , colorDetMC[2]);
        graphEMCALPi0FWHMMC->Draw("p,same,z");
        
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0FWHM, markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
        graphPCMEMCALPi0FWHM->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0FWHMMC, markerStyleDetMC[4], markerSizeDetMC[4]*0.55, colorDetMC[4] , colorDetMC[4]);
        graphPCMEMCALPi0FWHMMC->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphPCMPi0FWHM, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        graphPCMPi0FWHM->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPCMPi0FWHMMC, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
        graphPCMPi0FWHMMC->Draw("p,same,z");

        TLatex *labelLegendAMass    = new TLatex(0.13,0.06,"a)");
        SetStyleTLatex( labelLegendAMass, textSizeLabelsPixel,4, 1, 43);
        labelLegendAMass->Draw();

        TLatex *labelMassPerf       = new TLatex(0.13,0.87,"ALICE performance");
        SetStyleTLatex( labelMassPerf, textSizeLabelsPixel,4, 1, 43);
        labelMassPerf->Draw();        
        TLatex *labelMassEnergy     = new TLatex(0.13,0.78,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelMassEnergy, textSizeLabelsPixel,4, 1, 43);
        labelMassEnergy->Draw();
        TLatex *labelMassPi0        = new TLatex(0.13,0.69,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelMassPi0, textSizeLabelsPixel,4, 1, 43);
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
        
        TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, 0.23, 25., 1000., 120., 170);
        SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass, 
                                  textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
        histo2DAllPi0Mass->GetYaxis()->SetRangeUser(122.5,161.5);
        histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
        histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
        histo2DAllPi0Mass->GetXaxis()->SetLabelOffset(-0.015);
        histo2DAllPi0Mass->DrawCopy(); 

        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0Mass, markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
        graphEMCALPi0Mass->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0MassMC, markerStyleDetMC[2], markerSizeDetMC[2]*0.55, colorDetMC[2] , colorDetMC[2]);
        graphEMCALPi0MassMC->Draw("p,same,z");
        
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0Mass, markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
        graphPCMEMCALPi0Mass->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0MassMC, markerStyleDetMC[4], markerSizeDetMC[4]*0.55, colorDetMC[4] , colorDetMC[4]);
        graphPCMEMCALPi0MassMC->Draw("p,same,z");
        
        DrawGammaSetMarkerTGraphAsym(graphPCMPi0Mass, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        graphPCMPi0Mass->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPCMPi0MassMC, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
        graphPCMPi0MassMC->Draw("p,same,z");

        DrawGammaLines(0.23, 25. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);

        TLatex *labelLegendBMass            = new TLatex(0.13,0.22,"b)");
        SetStyleTLatex( labelLegendBMass, textSizeLabelsPixel,4, 1, 43);
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
        TLatex *textMassPCM                 = new TLatex(columnsLegendMass2[0],rowsLegendMass2[1],nameMeasGlobalLabel[0]);
        SetStyleTLatex( textMassPCM, textSizeLabelsPixel,4, 1, 43);
        textMassPCM->Draw();
        TLatex *textMassPCMEMCAL            = new TLatex(columnsLegendMass2[0],rowsLegendMass2[2],nameMeasGlobalLabel[4]);
        SetStyleTLatex( textMassPCMEMCAL, textSizeLabelsPixel,4, 1, 43);
        textMassPCMEMCAL->Draw();
        TLatex *textMassEMCAL               = new TLatex(columnsLegendMass2[0],rowsLegendMass2[3],nameMeasGlobalLabel[2]);
        SetStyleTLatex( textMassEMCAL, textSizeLabelsPixel,4, 1, 43);
        textMassEMCAL->Draw();
    
        //****************** second Column *************************************************
        TLatex *textMassData                = new TLatex(columnsLegendMass2[1],rowsLegendMass2[0] ,"Data");
        SetStyleTLatex( textMassData, textSizeLabelsPixel,4, 1, 43);
        textMassData->Draw();
        TLatex *textMassMC                  = new TLatex(columnsLegendMass2[2] ,rowsLegendMass2[0],"MC");
        SetStyleTLatex( textMassMC, textSizeLabelsPixel,4, 1, 43);
        textMassMC->Draw();
        
        TMarker* markerPCMPi0Mass        = CreateMarkerFromGraph(graphPCMPi0Mass,columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
        markerPCMPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2);
        TMarker* markerPCMEMCALPi0Mass   = CreateMarkerFromGraph(graphPCMEMCALPi0Mass,columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2 ,scaleMarkerMass2);
        markerPCMEMCALPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2);
        TMarker* markerEMCALPi0Mass      = CreateMarkerFromGraph(graphEMCALPi0Mass,columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2 ,scaleMarkerMass2);
        markerEMCALPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2);
    
        TMarker* markerPCMPi0MassMC      = CreateMarkerFromGraph(graphPCMPi0MassMC,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
        markerPCMPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[1]+ offsetMarkerYMass2);
        TMarker* markerPCMEMCALPi0MassMC = CreateMarkerFromGraph(graphPCMEMCALPi0MassMC,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2 ,scaleMarkerMass2);
        markerPCMEMCALPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[2]+ offsetMarkerYMass2);
        TMarker* markerEMCALPi0MassMC    = CreateMarkerFromGraph(graphEMCALPi0MassMC,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2 ,scaleMarkerMass2);
        markerEMCALPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[3]+ offsetMarkerYMass2);
        
    canvasMassWidthPi0->Update();
    canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // *********************************** Mass and width for pi0 at 2.76TeV including PHOS *********************************
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

        DrawGammaSetMarker(histoPHOSPi0FWHMMeV, markerStyleDet[1], markerSizeDet[1]*0.55, colorDet[1] , colorDet[1]);
        histoPHOSPi0FWHMMeV->Draw("p,same,e");
        DrawGammaSetMarker(histoPHOSPi0TrueFWHMMeV, markerStyleDetMC[1], markerSizeDetMC[1]*0.55, colorDetMC[1] , colorDetMC[1]);
        histoPHOSPi0TrueFWHMMeV->Draw("p,same,e");
        
        graphEMCALPi0FWHM->Draw("p,same,z");
        graphEMCALPi0FWHMMC->Draw("p,same,z");

        graphPCMEMCALPi0FWHM->Draw("p,same,z");
        graphPCMEMCALPi0FWHMMC->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphPCMPi0FWHM, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        graphPCMPi0FWHM->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPCMPi0FWHMMC, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
        graphPCMPi0FWHMMC->Draw("p,same,z");

        labelMassPerf->Draw();
        labelMassEnergy->Draw();
        labelMassPi0->Draw();        
        labelLegendAMass->Draw();

    padMassPi0->cd();
    padMassPi0->SetLogx();

        histo2DAllPi0Mass->DrawCopy(); 

        DrawGammaSetMarker(histoPHOSPi0Mass, markerStyleDet[1], markerSizeDet[1]*0.55, colorDet[1] , colorDet[1]);
        histoPHOSPi0Mass->Draw("p,same,e");
        DrawGammaSetMarker(histoPHOSPi0TrueMass, markerStyleDetMC[1], markerSizeDetMC[1]*0.55, colorDetMC[1] , colorDetMC[1]);
        histoPHOSPi0TrueMass->Draw("p,same,e");
        
        graphEMCALPi0Mass->Draw("p,same,z");
        graphEMCALPi0MassMC->Draw("p,same,z");
        
        graphPCMEMCALPi0Mass->Draw("p,same,z");
        graphPCMEMCALPi0MassMC->Draw("p,same,z");
        
        graphPCMPi0Mass->Draw("p,same,z");
        graphPCMPi0MassMC->Draw("p,same,z");

        DrawGammaLines(0.23, 25. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);

        labelLegendBMass->Draw();
        
        //********************************** Defintion of the Legend **************************************************    
        Double_t rowsLegendMass3[5]         = {0.8,0.6,0.4,0.2,0.01};
        Double_t offsetMarkerYMass3         = 0.06;

        padMassLegend2->cd();
        //****************** first Column **************************************************
        TLatex *textMassPCM2                = new TLatex(columnsLegendMass2[0],rowsLegendMass3[1],nameMeasGlobalLabel[0]);
        SetStyleTLatex( textMassPCM2, textSizeLabelsPixel,4, 1, 43);
        textMassPCM2->Draw();
        TLatex *textMassPCMEMCAL2           = new TLatex(columnsLegendMass2[0],rowsLegendMass3[2],nameMeasGlobalLabel[4]);
        SetStyleTLatex( textMassPCMEMCAL2, textSizeLabelsPixel,4, 1, 43);
        textMassPCMEMCAL2->Draw();
        TLatex *textMassEMCAL2              = new TLatex(columnsLegendMass2[0],rowsLegendMass3[3],nameMeasGlobalLabel[2]);
        SetStyleTLatex( textMassEMCAL2, textSizeLabelsPixel,4, 1, 43);
        textMassEMCAL2->Draw();
        TLatex *textMassPHOS2               = new TLatex(columnsLegendMass2[0],rowsLegendMass3[4],nameMeasGlobalLabel[1]);
        SetStyleTLatex( textMassPHOS2, textSizeLabelsPixel,4, 1, 43);
        textMassPHOS2->Draw();
    
        //****************** second Column *************************************************
        TLatex *textMassData2               = new TLatex(columnsLegendMass2[1],rowsLegendMass3[0] ,"Data");
        SetStyleTLatex( textMassData2, textSizeLabelsPixel,4, 1, 43);
        textMassData2->Draw();
        TLatex *textMassMC2                 = new TLatex(columnsLegendMass2[2] ,rowsLegendMass3[0],"MC");
        SetStyleTLatex( textMassMC2, textSizeLabelsPixel,4, 1, 43);
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
    // ******************************************* Mass and width for eta at 2.76TeV ****************************************
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

        
        TH2F * histo2DAllEtaFWHM    = new TH2F("histo2DAllEtaFWHM","histo2DAllEtaFWHM", 20, 0.33, 27. ,1000., -5, 92);
        SetStyleHistoTH2ForGraphs(histo2DAllEtaFWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                                  0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.85,0.28/(textsizeFacWidth*margin), 512, 505);
//      histo2DAllEtaFWHM->GetYaxis()->SetRangeUser(-1.,25.5);
        histo2DAllEtaFWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
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

        DrawGammaSetMarkerTGraphAsym(graphPCMEtaFWHM, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        graphPCMEtaFWHM->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaFWHMMC, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
        graphPCMEtaFWHMMC->Draw("p,same,z");
        
        labelMassPerf->Draw();
        labelLegendAMass->Draw();
        labelMassEnergy->Draw();
        TLatex *labelMassEta                = new TLatex(0.13,0.69,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelMassEta, textSizeLabelsPixel,4, 1, 43);
        labelMassEta->Draw();

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
        markerPCMPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.02 ,rowsLegendMass2[1]+ offsetMarkerYMass2);
        markerPCMEMCALPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.02 ,rowsLegendMass2[2]+ offsetMarkerYMass2);
        markerEMCALPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.02 ,rowsLegendMass2[3]+ offsetMarkerYMass2);
        
    padMassEta->cd();
    padMassEta->SetLogx();

        TH2F * histo2DAllEtaMass    = new TH2F("histo2DAllEtaMass","histo2DAllEtaMass",20, 0.33, 27., 1000., 500., 565);
        SetStyleHistoTH2ForGraphs(histo2DAllEtaMass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass, 
                                  textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
        histo2DAllEtaMass->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAllEtaMass->GetYaxis()->SetNdivisions(505);
        histo2DAllEtaMass->GetYaxis()->SetNoExponent(kTRUE);
        histo2DAllEtaMass->GetXaxis()->SetTickLength(0.05);
        histo2DAllEtaMass->GetXaxis()->SetLabelOffset(-0.015);
        histo2DAllEtaMass->DrawCopy(); 

        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaMass, markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
        graphEMCALEtaMass->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaMassMC, markerStyleDetMC[2], markerSizeDetMC[2]*0.55, colorDetMC[2] , colorDetMC[2]);
        graphEMCALEtaMassMC->Draw("p,same,z");
        
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaMass, markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
        graphPCMEMCALEtaMass->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaMassMC, markerStyleDetMC[4], markerSizeDetMC[4]*0.55, colorDetMC[4] , colorDetMC[4]);
        graphPCMEMCALEtaMassMC->Draw("p,same,z");
        
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaMass, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        graphPCMEtaMass->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaMassMC, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
        graphPCMEtaMassMC->Draw("p,same,z");

        DrawGammaLines(0.33, 27. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,0.3, kGray);
        
        labelLegendBMass->Draw();
        
    canvasMassWidthEta->Update();
    canvasMassWidthEta->Print(Form("%s/Eta_MassAndWidth.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************** Cross section for pi0 single measurement 2.76TeV ************************************
    // **********************************************************************************************************************
    
    TCanvas* canvasXSectionPi0  = new TCanvas("canvasXSectionPi0","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasXSectionPi0, 0.14, 0.02, 0.02, 0.08);
    canvasXSectionPi0->SetLogx();
    canvasXSectionPi0->SetLogy();
    
    TH2F * histo2DXSectionPi0;
    histo2DXSectionPi0          = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,0.23,70.,1000,2e-2,10e11);
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

        
        TLatex *labelEnergyXSectionPi0      = new TLatex(0.94,0.92,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelEnergyXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);
        labelEnergyXSectionPi0->Draw();
        TLatex *labelDetSysXSectionPi0      = new TLatex(0.94,0.88,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);
        labelDetSysXSectionPi0->Draw();

                                                
        TLegend* legendXSectionPi0          = GetAndSetLegend2(0.18, 0.13,0.40,0.13+(0.035*6), 0.035, 1, "", 42, 0);
        legendXSectionPi0->AddEntry(graphPCMPi0InvXSectionSys,nameMeasGlobalLabel[0],"fp");
        legendXSectionPi0->AddEntry(graphPHOSPi0InvXSectionSys,nameMeasGlobalLabel[1],"fp");
        legendXSectionPi0->AddEntry(graphEMCALPi0InvXSectionSys,nameMeasGlobalLabel[2],"fp");
        legendXSectionPi0->AddEntry(graphPCMEMCALPi0InvXSectionSys,nameMeasGlobalLabel[4],"fp");
        legendXSectionPi0->AddEntry(graphEMCALMergedPi0InvXSectionSys,nameMeasGlobalLabel[9],"fp");
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
  
    TCanvas* canvasXSectionEta      = new TCanvas("canvasXSectionEta","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasXSectionEta, 0.14, 0.02, 0.02, 0.08);
    canvasXSectionEta->SetLogx();
    canvasXSectionEta->SetLogy();
    
    TH2F * histo2DXSectionEta;
    histo2DXSectionEta              = new TH2F("histo2DXSectionEta","histo2DXSectionEta",11000,0.33, 27.,1000,6e0,10e10);
    SetStyleHistoTH2ForGraphs(histo2DXSectionEta, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.45);
    histo2DXSectionEta->GetXaxis()->SetMoreLogLabels();
    histo2DXSectionEta->GetXaxis()->SetLabelOffset(-0.01);
    histo2DXSectionEta->Draw("copy");

    
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvXSectionSys, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0], widthLinesBoxes, kTRUE);
        graphPCMEtaInvXSectionSys->Draw("E2same");

        //      graphPCMEMCALEtaInvXSectionSys->RemovePoint(0);
//      graphPCMEMCALEtaInvXSectionSys->RemovePoint(graphPCMEMCALEtaInvXSectionSys->GetN()-1);
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaInvXSectionSys, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4], widthLinesBoxes, kTRUE);
        graphPCMEMCALEtaInvXSectionSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaInvXSectionSys, markerStyleDet[2], markerSizeDet[4]*0.75, colorDet[2] , colorDet[2], widthLinesBoxes, kTRUE);
        graphEMCALEtaInvXSectionSys->Draw("E2same");
    
        DrawGammaSetMarker(histoPCMEtaInvXSectionStat, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
        histoPCMEtaInvXSectionStat->Draw("p,same,e");
//      DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvXSectionStat,markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
//      graphPCMPi0InvXSectionStat->Draw("p,same,z");
        DrawGammaSetMarker(histoPCMEMCALEtaInvXSectionStat, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4]);
        histoPCMEMCALEtaInvXSectionStat->Draw("p,same,e");

        DrawGammaSetMarker(histoEMCALEtaInvXSectionStat, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
        histoEMCALEtaInvXSectionStat->Draw("p,same,e");


        TLatex *labelEnergyXSectionEta      = new TLatex(0.94,0.92,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelEnergyXSectionEta, 0.035, 4, 1, 42, kTRUE, 31);
        labelEnergyXSectionEta->Draw();
        TLatex *labelDetSysXSectionEta      = new TLatex(0.94,0.88,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionEta, 0.035, 4, 1, 42, kTRUE, 31);
        labelDetSysXSectionEta->Draw();

        TLegend* legendXSectionEta          = GetAndSetLegend2(0.18, 0.13,0.40,0.13+(0.035*4), 0.035, 1, "", 42, 0);
        legendXSectionEta->AddEntry(graphPCMEtaInvXSectionSys,nameMeasGlobalLabel[0],"fp");
//      legendXSectionEta->AddEntry(graphPHOSPi0InvXSectionSys,nameMeasGlobalLabel[1],"fp");
        legendXSectionEta->AddEntry(graphEMCALEtaInvXSectionSys,nameMeasGlobalLabel[2],"fp");
        legendXSectionEta->AddEntry(graphPCMEMCALEtaInvXSectionSys,nameMeasGlobalLabel[4],"fp");
        legendXSectionEta->Draw();
   
    canvasXSectionEta->SaveAs(Form("%s/Eta_InvXSectionCompAllSystems.%s",outputDir.Data(),suffix.Data()));
    histo2DXSectionEta->Draw("copy");

        graphPCMEtaInvXSectionSys->Draw("E2same");
        graphEMCALEtaInvXSectionSys->Draw("E2same");
        graphPCMEMCALEtaInvXSectionSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSysA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphCombEtaInvXSectionSysA->Draw("E2same");
    
        histoPCMEtaInvXSectionStat->Draw("p,same,e");
        histoEMCALEtaInvXSectionStat->Draw("p,same,e");
        histoPCMEMCALEtaInvXSectionStat->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStatA, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombEtaInvXSectionStatA->Draw("p,same,z");


        labelEnergyXSectionEta->Draw();
        labelDetSysXSectionEta->Draw();

        legendXSectionEta->AddEntry(graphCombEtaInvXSectionSysA,"comb","fp");
        legendXSectionEta->Draw();

//      DrawGammaSetMarkerTGraphAsym(graphChargedHadronsStatPP, markerStyleDet[1], markerSizeDet[1]*0.75, kBlack , kBlack, widthLinesBoxes);
//      graphChargedHadronsStatPP->Draw("E2same");
        
   canvasXSectionEta->SaveAs(Form("%s/Eta_InvXSectionCompAllSystems_Comb.%s",outputDir.Data(),suffix.Data()));


    // **********************************************************************************************************************
    // ******************************** Acceptance * Efficiency for pi0 single measurement 2.76TeV **************************
    // **********************************************************************************************************************
    textSizeLabelsPixel             = 55;
    Double_t textSizeLabelsRel      = 55./1200;
    cout << textSizeLabelsRel << endl;
    
    TCanvas* canvasAcceptanceTimesEff       = new TCanvas("canvasAcceptanceTimesEff", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasAcceptanceTimesEff,  0.1, 0.01, 0.015, 0.095);
    canvasAcceptanceTimesEff->SetLogy(1);
    canvasAcceptanceTimesEff->SetLogx(1);
    
        TH2F * histo2DAccEff;
        histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, 0.23,  70, 1000, 8e-5, 2e-0 );
        SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec} / #it{P}"),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
        histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
        histo2DAccEff->GetXaxis()->SetLabelOffset(-0.01);
        histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAccEff->DrawCopy(); 

        DrawGammaSetMarkerTGraphAsym(graphPCMPi0AccTimesEff, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        graphPCMPi0AccTimesEff->Draw("p,same,z");
        
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0AccTimesEff, markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
        graphPCMEMCALPi0AccTimesEff->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0AccTimesEff, markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
        graphEMCALPi0AccTimesEff->Draw("p,same,z");
        
        DrawGammaSetMarkerTGraphAsym(graphEMCALMergedPi0AccTimesEffDivPur, markerStyleDet[9], markerSizeDet[9]*0.55, colorDet[9] , colorDet[9]);
        graphEMCALMergedPi0AccTimesEffDivPur->Draw("p,same,e");

        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.55, 0.13, 0.83, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffiAccPi0->AddEntry(graphPCMPi0AccTimesEff,nameMeasGlobalLabel[0],"p");
        legendEffiAccPi0->AddEntry(graphPCMEMCALPi0AccTimesEff,nameMeasGlobalLabel[4],"p");
        legendEffiAccPi0->AddEntry(graphEMCALPi0AccTimesEff,nameMeasGlobalLabel[2],"p");
        legendEffiAccPi0->AddEntry(graphEMCALMergedPi0AccTimesEffDivPur,nameMeasGlobalLabel[9],"p");
        legendEffiAccPi0->Draw();

        TLatex *labelPerfEffi               = new TLatex(0.15,0.92,"ALICE performance");
        SetStyleTLatex( labelPerfEffi, textSizeLabelsRel,4);
        labelPerfEffi->Draw();
        TLatex *labelEnergyEffi             = new TLatex(0.15,0.87,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelEnergyEffi, textSizeLabelsRel,4);
        labelEnergyEffi->Draw();
        TLatex *labelDetSysEffiPi0          = new TLatex(0.15,0.82,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysEffiPi0, textSizeLabelsRel,4);
        labelDetSysEffiPi0->Draw();

        
    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ************************* Acceptance * Efficiency for pi0 single measurement 2.76TeV inc PHOS ************************
    // **********************************************************************************************************************
    
    canvasAcceptanceTimesEff->cd();
        histo2DAccEff->DrawCopy();
        
        DrawGammaSetMarkerTGraphAsym(graphPHOSPi0AccTimesEff, markerStyleDet[1], markerSizeDet[1]*0.55, colorDet[1] , colorDet[1]);
        graphPHOSPi0AccTimesEff->Draw("p,same,z");

        graphPCMPi0AccTimesEff->Draw("p,same,z");
        graphPCMEMCALPi0AccTimesEff->Draw("p,same,z");
        graphEMCALPi0AccTimesEff->Draw("p,same,z");
        graphEMCALMergedPi0AccTimesEffDivPur->Draw("p,same,e");

        TLegend* legendEffiAccPi02          = GetAndSetLegend2(0.55, 0.13, 0.83, 0.13+(5*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffiAccPi02->AddEntry(graphPCMPi0AccTimesEff,nameMeasGlobalLabel[0],"p");
        legendEffiAccPi02->AddEntry(graphPCMEMCALPi0AccTimesEff,nameMeasGlobalLabel[4],"p");
        legendEffiAccPi02->AddEntry(graphEMCALPi0AccTimesEff,nameMeasGlobalLabel[2],"p");
        legendEffiAccPi02->AddEntry(graphEMCALMergedPi0AccTimesEffDivPur,nameMeasGlobalLabel[9],"p");
        legendEffiAccPi02->AddEntry(graphPHOSPi0AccTimesEff,nameMeasGlobalLabel[1],"p");
        legendEffiAccPi02->Draw();

        labelPerfEffi->Draw();
        labelEnergyEffi->Draw();
        labelDetSysEffiPi0->Draw();
        
    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff_incPHOS.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************** Acceptance * Efficiency for eta single measurement 2.76TeV **************************
    // **********************************************************************************************************************
    canvasAcceptanceTimesEff->cd();
    
        TH2F * histo2DAccEffEta;
        histo2DAccEffEta                = new TH2F("histo2DAccEffEta", "histo2DAccEffEta",1000, 0.33,  25, 1000, 8e-4, 1.5e-0 );
        SetStyleHistoTH2ForGraphs( histo2DAccEffEta, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec} / #it{P}"),  
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
        histo2DAccEffEta->GetYaxis()->SetLabelOffset(0.001);
        histo2DAccEffEta->GetXaxis()->SetLabelOffset(-0.01);
        histo2DAccEffEta->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAccEffEta->DrawCopy(); 
    
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaAccTimesEff, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        graphPCMEtaAccTimesEff->Draw("p,same,z");
        
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaAccTimesEff, markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
        graphPCMEMCALEtaAccTimesEff->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaAccTimesEff, markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
        graphEMCALEtaAccTimesEff->Draw("p,same,z");
        
        TLegend* legendEffiAccEta           = GetAndSetLegend2(0.62, 0.13, 0.9, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffiAccEta->AddEntry(graphPCMEtaAccTimesEff,nameMeasGlobalLabel[0],"p");
        legendEffiAccEta->AddEntry(graphPCMEMCALEtaAccTimesEff,nameMeasGlobalLabel[4],"p");
        legendEffiAccEta->AddEntry(graphEMCALEtaAccTimesEff,nameMeasGlobalLabel[2],"p");
        legendEffiAccEta->Draw();

        labelPerfEffi->Draw();
        labelEnergyEffi->Draw();
        TLatex *labelDetSysEffiEta          = new TLatex(0.15,0.82,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysEffiEta, textSizeLabelsRel,4);
        labelDetSysEffiEta->Draw();

        
    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Eta_AcceptanceTimesEff.%s",outputDir.Data(),suffix.Data()));

    TCanvas* canvasTriggerEff       = new TCanvas("canvasTriggerEff", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasTriggerEff,  0.09, 0.01, 0.015, 0.095);
//  canvasTriggerEff->SetLogx(1);
    
        TH2F * histo2DTriggerEffPi0;
        histo2DTriggerEffPi0                = new TH2F("histo2DTriggerEffPi0", "histo2DTriggerEffPi0",1000, 0.,  21, 1000, 0, 1.15 );
        SetStyleHistoTH2ForGraphs( histo2DTriggerEffPi0, "#it{p}_{T} (GeV/#it{c})", "#it{#kappa}_{Trig}",  
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 0.95);//(#times #epsilon_{pur})
            histo2DTriggerEffPi0->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DTriggerEffPi0->DrawCopy(); 

        DrawGammaSetMarker(histoPCMEMCALPi0TriggerEff[0], markerTriggMC[2], sizeTrigg[2]*1.5, colorTriggShade[2] , colorTriggShade[2]);
        histoPCMEMCALPi0TriggerEff[0]->Draw("p,same,e");
        DrawGammaSetMarker(histoPCMEMCALPi0TriggerEff[1], markerTriggMC[3], sizeTrigg[3]*1.5, colorTriggShade[3] , colorTriggShade[3]);
        histoPCMEMCALPi0TriggerEff[1]->Draw("p,same,e");
        DrawGammaSetMarker(histoPCMEMCALPi0TriggerEff[3], markerTriggMC[5], sizeTrigg[5]*1.5, colorTriggShade[5] , colorTriggShade[5]);
        histoPCMEMCALPi0TriggerEff[3]->Draw("p,same,e");
        DrawGammaSetMarker(histoEMCALPi0TriggerEff[0], markerTrigg[2], sizeTrigg[2]*1.5, colorTrigg[2] , colorTrigg[2]);
        histoEMCALPi0TriggerEff[0]->Draw("p,same,e");
        DrawGammaSetMarker(histoEMCALPi0TriggerEff[1], markerTrigg[3], sizeTrigg[3]*1.5, colorTrigg[3] , colorTrigg[3]);
        histoEMCALPi0TriggerEff[1]->Draw("p,same,e");
        DrawGammaSetMarker(histoEMCALPi0TriggerEff[3], markerTrigg[5], sizeTrigg[5]*1.5, colorTrigg[5] , colorTrigg[5]);
        histoEMCALPi0TriggerEff[3]->Draw("p,same,e");


        TLatex *labelPerfTriggEff               = new TLatex(0.64,0.425,"ALICE simulation");
        SetStyleTLatex( labelPerfTriggEff, textSizeLabelsRel,4);
        labelPerfTriggEff->Draw();
        TLatex *labelEnergyTriggEff             = new TLatex(0.64,0.378,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelEnergyTriggEff, textSizeLabelsRel,4);
        labelEnergyTriggEff->Draw();
        TLatex *labelDetSysTriggEffPi0          = new TLatex(0.64,0.33,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysTriggEffPi0, textSizeLabelsRel,4);
        labelDetSysTriggEffPi0->Draw();

        TLegend* legendTriggEffPi0              = GetAndSetLegend2(0.475, 0.13, 0.9, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel,3, "", 43, 0.15);
        legendTriggEffPi0->AddEntry((TObject*)0, Form("%s    ",nameTrigger[2].Data()),"");
        legendTriggEffPi0->AddEntry(histoPCMEMCALPi0TriggerEff[0],"       ","p");
        legendTriggEffPi0->AddEntry(histoEMCALPi0TriggerEff[0],"","p");
        legendTriggEffPi0->AddEntry((TObject*)0, nameTrigger[3].Data(),"");
        legendTriggEffPi0->AddEntry(histoPCMEMCALPi0TriggerEff[1],"    ","p");
        legendTriggEffPi0->AddEntry(histoEMCALPi0TriggerEff[1],"","p");
        legendTriggEffPi0->AddEntry((TObject*)0, nameTrigger[5].Data(),"");
        legendTriggEffPi0->AddEntry(histoPCMEMCALPi0TriggerEff[3],"    ","p");
        legendTriggEffPi0->AddEntry(histoEMCALPi0TriggerEff[3],"","p");
        legendTriggEffPi0->Draw();

        TLatex *labelPCMEMCALTriggEff           = new TLatex(0.64,0.275,nameMeasGlobalLabel[4]);
        SetStyleTLatex( labelPCMEMCALTriggEff, textSizeLabelsRel,4);
        labelPCMEMCALTriggEff->Draw();
        TLatex *labelEMCALTriggEff              = new TLatex(0.85,0.275,nameMeasGlobalLabel[2]);
        SetStyleTLatex( labelEMCALTriggEff, textSizeLabelsRel,4);
        labelEMCALTriggEff->Draw();
        
    canvasTriggerEff->Update();
    canvasTriggerEff->Print(Form("%s/Pi0_TriggerEff.%s",outputDir.Data(),suffix.Data()));
    
    canvasTriggerEff->cd();
    histo2DTriggerEffPi0->DrawCopy(); 
        DrawGammaSetMarker(histoPCMEMCALEtaTriggerEff[0], markerTriggMC[2], sizeTrigg[2]*1.5, colorTriggShade[2] , colorTriggShade[2]);
        histoPCMEMCALEtaTriggerEff[0]->Draw("p,same,e");
        DrawGammaSetMarker(histoPCMEMCALEtaTriggerEff[1], markerTriggMC[3], sizeTrigg[3]*1.5, colorTriggShade[3] , colorTriggShade[3]);
        histoPCMEMCALEtaTriggerEff[1]->Draw("p,same,e");
        DrawGammaSetMarker(histoPCMEMCALEtaTriggerEff[3], markerTriggMC[5], sizeTrigg[5]*1.5, colorTriggShade[5] , colorTriggShade[5]);
        histoPCMEMCALEtaTriggerEff[3]->Draw("p,same,e");
        DrawGammaSetMarker(histoEMCALEtaTriggerEff[0], markerTrigg[2], sizeTrigg[2]*1.5, colorTrigg[2] , colorTrigg[2]);
        histoEMCALEtaTriggerEff[0]->Draw("p,same,e");
        DrawGammaSetMarker(histoEMCALEtaTriggerEff[1], markerTrigg[3], sizeTrigg[3]*1.5, colorTrigg[3] , colorTrigg[3]);
        histoEMCALEtaTriggerEff[1]->Draw("p,same,e");
        DrawGammaSetMarker(histoEMCALEtaTriggerEff[3], markerTrigg[5], sizeTrigg[5]*1.5, colorTrigg[5] , colorTrigg[5]);
        histoEMCALEtaTriggerEff[3]->Draw("p,same,e");

        labelPerfTriggEff->Draw();
        labelEnergyTriggEff->Draw();
        TLatex *labelDetSysTriggEffEta          = new TLatex(0.64,0.33,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysTriggEffEta, 0.85*textSizeLabelsRel,4);
        labelDetSysTriggEffEta->Draw();
        labelPCMEMCALTriggEff->Draw();
        labelEMCALTriggEff->Draw();
        legendTriggEffPi0->Draw();

        
    canvasTriggerEff->Update();
    canvasTriggerEff->Print(Form("%s/Eta_TriggerEff.%s",outputDir.Data(),suffix.Data()));
    
    
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
                        legendEffSecCorrPi0->AddEntry(graphPi0EffSecCorrFromX[k][i],nameMeasGlobalLabel[i],"p");
                    }
                }    
                legendEffSecCorrPi0->Draw();
                TLatex *labelPerfSecCorr               = new TLatex(0.15,0.90,"ALICE performance");
                SetStyleTLatex( labelPerfSecCorr, textSizeLabelsRel,4);
                labelPerfSecCorr->Draw();
                TLatex *labelEnergySecCorr             = new TLatex(0.15,0.85,collisionSystem2760GeV.Data());
                SetStyleTLatex( labelEnergySecCorr, textSizeLabelsRel,4);
                labelEnergySecCorr->Draw();
                TLatex *labelDetSysSecCorrPi0          = new TLatex(0.15,0.80,"#pi^{0} #rightarrow #gamma#gamma");
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

        TH2F * ratio2DTheoryPP       = new TH2F("ratio2DTheoryPP","ratio2DTheoryPP",1000,0.23,70.,1000,0.2,3.6);
        SetStyleHistoTH2ForGraphs(ratio2DTheoryPP, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
                                  0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
        ratio2DTheoryPP->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio2DTheoryPP->GetYaxis()->SetNdivisions(505);
        ratio2DTheoryPP->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DTheoryPP->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DTheoryPP->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DTheoryPP->GetXaxis()->SetLabelFont(42);
        ratio2DTheoryPP->GetYaxis()->SetLabelFont(42);
        ratio2DTheoryPP->DrawCopy(); 

        DrawGammaSetMarkerTGraphErr(graphRatioPi0Pythia8ToFit, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2 );
        graphRatioPi0Pythia8ToFit->Draw("3,same");
        
        DrawGammaSetMarker(histoRatioPi0Pythia8ToFit, 24, 1.5, kRed+2 , kRed+2);  
        histoRatioPi0Pythia8ToFit->SetLineWidth(widthCommonFit);
//      histoRatioPi0Pythia8ToFit->GetXaxis()->SetRangeUser(0.5,14);
        histoRatioPi0Pythia8ToFit->Draw("same,hist,l");
        
        graphRatioPi0CombNLODSS14->RemovePoint(0);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS14, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
        graphRatioPi0CombNLODSS14->Draw("3,same");

//         DrawGammaNLOTGraph( graphRatioPi0CombCGC, widthCommonFit, styleLineCGC, colorCGC);
//         graphRatioPi0CombCGC->Draw("same,c");

        DrawGammaNLOTGraph( graphRatioPi0CombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, kGray+1);
        graphRatioPi0CombNLOMuHalf->Draw("same,c");
        DrawGammaNLOTGraph( graphRatioPi0CombNLOMuOne, widthCommonFit, styleLineNLOMuOne, kGray+1);
        graphRatioPi0CombNLOMuOne->Draw("same,c");
        DrawGammaNLOTGraph( graphRatioPi0CombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, kGray+2);
        graphRatioPi0CombNLOMuTwo->Draw("same,c");


        TGraphAsymmErrors* graphRatioPi0CombCombFitStatAWOXErr = (TGraphAsymmErrors*)graphRatioPi0CombCombFitStatA->Clone("graphRatioPi0CombCombFitStatAWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatAWOXErr);
        
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatAWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatAWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA->Draw("2,same");
        graphRatioPi0CombCombFitStatAWOXErr->Draw("p,same");

        TBox* boxErrorSigma2760GeVRatioPi0 = CreateBoxConv(kGray+1, 0.28, 1.-(0.0251), 0.32, 1.+(0.0251));
        boxErrorSigma2760GeVRatioPi0->Draw();
        DrawGammaLines(0.23, 70.,1., 1.,0.1,kGray);

        // upper left corner
        TLatex *labelRatioTheoryNLO   = new TLatex(0.155,0.92,"NLO, PDF: CTEQ6M5 FF: DSS07");
        SetStyleTLatex( labelRatioTheoryNLO, 0.85*textsizeLabelsPP,4);
        labelRatioTheoryNLO->Draw();
        TLegend* legendRatioPi0OldTheo= GetAndSetLegend2(0.15,0.905-3*0.9*textsizeLabelsPP,0.4,0.905, 0.85* textSizeLabelsPixel, 1, "", 43, 0.2);
        legendRatioPi0OldTheo->AddEntry(graphRatioPi0CombNLOMuHalf, "#mu = 0.5 #it{p}_{T}", "l");  
        legendRatioPi0OldTheo->AddEntry(graphRatioPi0CombNLOMuOne,  "#mu = #it{p}_{T}", "l");  
        legendRatioPi0OldTheo->AddEntry(graphRatioPi0CombNLOMuTwo,  "#mu = 2 #it{p}_{T}", "l");  
        legendRatioPi0OldTheo->Draw();
        
        // lower left corner
        TLegend* legendRatioPi0OldTheo3= GetAndSetLegend2(0.15,0.14,0.4,0.15+2*0.85*textsizeLabelsPP, 0.85* textSizeLabelsPixel, 1, "", 43, 0.2);
        legendRatioPi0OldTheo3->AddEntry(graphRatioPi0CombNLODSS14,  "NLO, PDF: MSTW FF: DSS14", "f");  
        legendRatioPi0OldTheo3->AddEntry(histoRatioPi0Pythia8ToFit,  "PYTHIA 8.2, Monash 2013", "l");
        legendRatioPi0OldTheo3->Draw();

        // upper right corner
        TLegend* legendRatioPi0OldTheo2= GetAndSetLegend2(0.72,0.905-1*0.85*textsizeLabelsPP,0.95,0.905, 0.85* textSizeLabelsPixel, 1, "", 43, 0.27);
        legendRatioPi0OldTheo2->AddEntry(graphRatioPi0CombCombFitSysA,"#pi^{0} ALICE","pf");
        legendRatioPi0OldTheo2->Draw();
        
        TLatex *labelRatioTheoryPP2   = new TLatex(0.73,0.92,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelRatioTheoryPP2, 0.85*textsizeLabelsPP,4);
        labelRatioTheoryPP2->Draw();
    
        
    canvasRatioPP->Update();
    canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PP.%s",outputDir.Data(),suffix.Data()));

    canvasRatioPP->cd();
        ratio2DTheoryPP->GetYaxis()->SetNdivisions(512);
        ratio2DTheoryPP->GetYaxis()->SetRangeUser(0.5,1.9);
        ratio2DTheoryPP->DrawCopy(); 

        graphRatioPi0Pythia8ToFit->Draw("3,same");
        histoRatioPi0Pythia8ToFit->Draw("same,hist,l");
        
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS14, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
        graphRatioPi0CombNLODSS14->Draw("3,same");
        

//         DrawGammaNLOTGraph( graphRatioPi0CombCGC, widthCommonFit, styleLineCGC, colorCGC);
//         graphRatioPi0CombCGC->Draw("same,c");
        
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatAWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatAWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA->Draw("2,same");
        graphRatioPi0CombCombFitStatAWOXErr->Draw("p,same");

        boxErrorSigma2760GeVRatioPi0->Draw();
        DrawGammaLines(0.23, 70.,1., 1.,0.1,kGray);

        TLegend* legendRatioTheoryPPNew= GetAndSetLegend2(0.15,0.15,0.4,0.15+0.85*textsizeLabelsPP*3, 0.85*textSizeLabelsPixel);
        legendRatioTheoryPPNew->AddEntry(graphRatioPi0CombCombFitSysA,"#pi^{0} ALICE","pf");
        legendRatioTheoryPPNew->AddEntry(graphRatioPi0CombNLODSS14,"NLO, PDF: MSTW FF: DSS14","f");
//         legendRatioTheoryPPNew->AddEntry((TObject*)0,"Phys.Rev. D91 no. 1, (2015) 014035","");
        legendRatioTheoryPPNew->AddEntry(histoRatioPi0Pythia8ToFit,  "PYTHIA 8.2, Monash 2013", "l");  
        legendRatioTheoryPPNew->Draw();

        TLatex *labelRatioTheoryPP   = new TLatex(0.95,0.92,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelRatioTheoryPP, 0.85*textsizeLabelsPP,4, 1, 42, kTRUE, 31);        
        labelRatioTheoryPP->Draw();
        
    canvasRatioPP->Update();
    canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PP_NewTheory.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************************* Comparison to theory calculations Eta ************************************
    // **********************************************************************************************************************    
    canvasRatioPP->cd();
    canvasRatioPP->SetLogx();

        TH2F * ratio2DTheoryPPEta       = new TH2F("ratio2DTheoryPPEta","ratio2DTheoryPPEta",1000,0.33, 27.,1000,0.2,3.6);
        SetStyleHistoTH2ForGraphs(ratio2DTheoryPPEta, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
                                  0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
        ratio2DTheoryPPEta->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio2DTheoryPPEta->GetYaxis()->SetNdivisions(505);
        ratio2DTheoryPPEta->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DTheoryPPEta->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DTheoryPPEta->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DTheoryPPEta->GetXaxis()->SetLabelFont(42);
        ratio2DTheoryPPEta->GetYaxis()->SetLabelFont(42);
        ratio2DTheoryPPEta->DrawCopy(); 

        DrawGammaSetMarkerTGraphErr(graphRatioEtaPythia8ToFit, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2 );
        graphRatioEtaPythia8ToFit->Draw("3,same");

        DrawGammaSetMarker(histoRatioEtaPythia8ToFit, 24, 1.5, kRed+2 , kRed+2);  
        histoRatioEtaPythia8ToFit->SetLineWidth(widthCommonFit);
//      histoRatioEtaPythia8ToFit->GetXaxis()->SetRangeUser(0.5,14);
        histoRatioEtaPythia8ToFit->Draw("same,hist,l");
        
        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
        graphRatioEtaCombNLOMuHalf->Draw("same,c");
        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
        graphRatioEtaCombNLOMuOne->Draw("same,c");
        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
        graphRatioEtaCombNLOMuTwo->Draw("same,c");
        

        TGraphAsymmErrors* graphRatioEtaCombCombFitStatAWOXErr = (TGraphAsymmErrors*)graphRatioEtaCombCombFitStatA->Clone("graphRatioEtaCombCombFitStatAWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombCombFitStatAWOXErr);
        
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatAWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioEtaCombCombFitStatAWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioEtaCombCombFitSysA->SetLineWidth(0);
        graphRatioEtaCombCombFitSysA->Draw("2,same");
        graphRatioEtaCombCombFitStatAWOXErr->Draw("p,same");

        TBox* boxErrorSigma2760GeVRatioEta = CreateBoxConv(kGray+1, 0.4, 1.-(0.0251), 0.44, 1.+(0.0251));        
        boxErrorSigma2760GeVRatioEta->Draw();
        DrawGammaLines(0.33, 27.,1., 1.,0.1,kGray);

        // upper left corner
        TLatex *labelRatioTheoryNLOEta   = new TLatex(0.155,0.92,"NLO, PDF: CTEQ6M5 FF: AESSS");
        SetStyleTLatex( labelRatioTheoryNLOEta, 0.85*textsizeLabelsPP,4);
        labelRatioTheoryNLOEta->Draw();
        legendRatioPi0OldTheo->Draw();
        
        // lower left corner
        TLegend* legendRatioEtaOldTheo3= GetAndSetLegend2(0.15,0.14,0.4,0.15+1*0.85*textsizeLabelsPP, 0.85* textSizeLabelsPixel, 1, "", 43, 0.2);
        legendRatioEtaOldTheo3->AddEntry(histoRatioPi0Pythia8ToFit,  "PYTHIA 8.2, Monash 2013", "l");  
        legendRatioEtaOldTheo3->Draw();

        // upper right corner
        TLegend* legendRatioEtaOldTheo2= GetAndSetLegend2(0.72,0.905-1*0.85*textsizeLabelsPP,0.95,0.905, 0.85* textSizeLabelsPixel, 1, "", 43, 0.27);
        legendRatioEtaOldTheo2->AddEntry(graphRatioPi0CombCombFitSysA,"#eta ALICE","pf");
        legendRatioEtaOldTheo2->Draw();
        labelRatioTheoryPP2->Draw();
            
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

        DrawGammaSetMarkerTGraphErr(graphPythia8InvXSectionPi0, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2 );
        graphPythia8InvXSectionPi0->Draw("3,same");
        
        DrawGammaSetMarker(histoPythia8InvXSectionPi0, 24, 1.5, kRed+2 , kRed+2);  
        histoPythia8InvXSectionPi0->SetLineWidth(widthCommonFit);
        histoPythia8InvXSectionPi0->Draw("same,hist,l");

        graphNLODSS14Calc->RemovePoint(0);
        DrawGammaSetMarkerTGraphAsym(graphNLODSS14Calc, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
        graphNLODSS14Calc->Draw("3,same");

//         DrawGammaNLOTGraph( graphCGCPi0, widthCommonFit, styleLineCGC, colorCGC);
//         graphCGCPi0->Draw("same,c");


        TGraphAsymmErrors* graphCombPi0InvXSectionStatAWOXErr = (TGraphAsymmErrors*)graphCombPi0InvXSectionStatA->Clone("graphCombPi0InvXSectionStatAWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphCombPi0InvXSectionStatAWOXErr);
        
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombPi0InvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatAWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
        graphCombPi0InvXSectionStatAWOXErr->Draw("p,same,z");

        DrawGammaSetMarkerTF1( fitTCMInvXSectionPi0, 7, 2, kGray+2); 
        fitTCMInvXSectionPi0->Draw("same");        

//         DrawGammaSetMarkerTF1( fitPowInvXSectionPi0, 7, 2, kBlue-7); 
//         fitPowInvXSectionPi0->Draw("same");        
        
        
        TLatex *labelEnergyXSectionPaper= new TLatex(0.95, 0.91, collisionSystem2760GeV.Data());
        SetStyleTLatex( labelEnergyXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
        labelEnergyXSectionPaper->Draw();
        TLatex *labelALICEXSectionPaper= new TLatex(0.95,0.87,"ALICE");
        SetStyleTLatex( labelALICEXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
        labelALICEXSectionPaper->Draw();
        TLatex *labelDetSysXSectionPaper= new TLatex(0.95,0.83,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
        labelDetSysXSectionPaper->Draw();
        
        TGraphErrors* dummyNorm     = new TGraphErrors(1);
        dummyNorm->SetPoint(0,1,0.0251);
        DrawGammaSetMarkerTGraphErr(dummyNorm, 0, 0, kGray+1, kGray+1, widthLinesBoxes, kTRUE, kGray+1);
        
        TLegend* legendXsectionPaper    = GetAndSetLegend2(0.17, 0.03, 0.5, 0.03+0.05*5, textSizeLabelsPixel, 1, "", 43, 0.2);
        legendXsectionPaper->AddEntry(graphCombPi0InvXSectionSysA,"Data","pf");
        legendXsectionPaper->AddEntry(dummyNorm,"Norm. unc. 2.5%","f");
        legendXsectionPaper->AddEntry(fitTCMInvXSectionPi0,"TCM fit","l");
        legendXsectionPaper->AddEntry(graphNLODSS14Calc,"NLO, PDF: MSTW, FF: DSS14","f");
//         legendXsectionPaper->AddEntry((TObject*)0,"Phys.Rev. D91 no. 1, (2015) 014035","");
        legendXsectionPaper->AddEntry(histoPythia8InvXSectionPi0,"PYTHIA 8.2, Monash 2013","l");
        legendXsectionPaper->Draw();
        
    padInvSectionNLORatio->cd();
    padInvSectionNLORatio->SetLogx(1);
        TH2F * ratio2DNLO               = new TH2F("ratio2DNLO","ratio2DNLO",1000,0.23,70.,1000,0.5,1.9);
        SetStyleHistoTH2ForGraphs(ratio2DNLO, "#it{p}_{T} (GeV/#it{c})","#frac{NLO, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle, 
                                  0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 512);
//         ratio2DNLO->GetYaxis()->SetNdivisions(512);
        ratio2DNLO->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DNLO->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DNLO->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DNLO->GetXaxis()->SetLabelFont(42);
        ratio2DNLO->GetYaxis()->SetLabelFont(42);
        ratio2DNLO->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DNLO->GetXaxis()->SetTickLength(0.07);
        ratio2DNLO->DrawCopy();

//         DrawGammaNLOTGraph( graphRatioPi0CombCGC, widthCommonFit, styleLineCGC, colorCGC);
//         graphRatioPi0CombCGC->Draw("same,c");

        boxErrorSigma2760GeVRatioPi0->Draw();
        
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS14, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
        graphRatioPi0CombNLODSS14->Draw("3,same");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatAWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatAWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA->Draw("2,same");
        graphRatioPi0CombCombFitStatAWOXErr->Draw("p,same");
        
        DrawGammaLines(0.23, 70.,1., 1.,0.1,kGray);
        
    padInvSectionPythiaRatio->cd();
    padInvSectionPythiaRatio->SetLogx(1);
        TH2F * ratio2DPythia            = new TH2F("ratio2DPythia","ratio2DPythia",1000,0.23,70.,1000,0.5,1.9);
        SetStyleHistoTH2ForGraphs(ratio2DPythia, "#it{p}_{T} (GeV/#it{c})","#frac{PYTHIA, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown, 
                                  0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.2/(textsizeFacXSecDown*marginXSec), 510, 512);
//         ratio2DPythia->GetYaxis()->SetNdivisions(512);
        ratio2DPythia->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DPythia->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DPythia->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DPythia->GetXaxis()->SetLabelFont(42);
        ratio2DPythia->GetYaxis()->SetLabelFont(42);
        ratio2DPythia->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DPythia->GetXaxis()->SetTickLength(0.06);
        ratio2DPythia->GetYaxis()->SetTickLength(0.04);
        ratio2DPythia->DrawCopy();

        graphRatioPi0Pythia8ToFit->Draw("3,same");
        DrawGammaSetMarker(histoRatioPi0Pythia8ToFit, 24, 1.5, kRed+2 , kRed+2);  
        histoRatioPi0Pythia8ToFit->SetLineWidth(widthCommonFit);
//      histoRatioPi0Pythia8ToFit->GetXaxis()->SetRangeUser(0.5,14);
        histoRatioPi0Pythia8ToFit->Draw("same,hist,l");
        
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatAWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatAWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA->Draw("2,same");
        graphRatioPi0CombCombFitStatAWOXErr->Draw("p,same");
        
        boxErrorSigma2760GeVRatioPi0->Draw();
        
        DrawGammaLines(0.23, 70.,1., 1.,0.1,kGray);
        
    canvasInvSectionPaper->Print(Form("%s/Pi0_InvXSectionWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));
    
    padInvSectionSpec->cd();
    padInvSectionSpec->SetLogy(1);
    padInvSectionSpec->SetLogx(1);
        SetStyleHistoTH2ForGraphs(histo2DXSectionEta, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",
                                0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.2/(textsizeFacXSecUp*marginXSec));
        histo2DXSectionEta->GetXaxis()->SetMoreLogLabels();
        histo2DXSectionEta->GetXaxis()->SetLabelOffset(+0.01);
        histo2DXSectionEta->Draw();

        DrawGammaSetMarkerTGraphErr(graphPythia8InvXSectionEta, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2 );
        graphPythia8InvXSectionEta->Draw("3,same");

        histoPythia8InvXSectionEta->GetXaxis()->SetRangeUser(0.4,25);
        DrawGammaSetMarker(histoPythia8InvXSectionEta, 24, 1.5, kRed+2 , kRed+2);  
        histoPythia8InvXSectionEta->SetLineWidth(widthCommonFit);
        histoPythia8InvXSectionEta->Draw("same,hist,l");

        DrawGammaNLOTGraph( graphNLOCalcEtaMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
        graphNLOCalcEtaMuHalf->Draw("same,c");
        DrawGammaNLOTGraph( graphNLOCalcEtaMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
        graphNLOCalcEtaMuOne->Draw("same,c");
        DrawGammaNLOTGraph( graphNLOCalcEtaMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
        graphNLOCalcEtaMuTwo->Draw("same,c");

        TGraphAsymmErrors* graphCombEtaInvXSectionStatAWOXErr = (TGraphAsymmErrors*)graphCombEtaInvXSectionStatA->Clone("graphCombEtaInvXSectionStatAWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphCombEtaInvXSectionStatAWOXErr);

        
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombEtaInvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStatAWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
        graphCombEtaInvXSectionStatAWOXErr->Draw("p,same,z");

        DrawGammaSetMarkerTF1( fitTCMInvXSectionEta, 7, 2, kGray+2); 
        fitTCMInvXSectionEta->Draw("same");

        // labels lower left corner
        TLegend* legendXsectionPaperEta    = GetAndSetLegend2(0.17, 0.03, 0.5, 0.03+0.05*3, textSizeLabelsPixel, 1, "", 43, 0.2);
        legendXsectionPaperEta->AddEntry(graphCombPi0InvXSectionSysA,"Data","pf");
        legendXsectionPaperEta->AddEntry(dummyNorm,"Norm. unc. 2.5%","f");
        legendXsectionPaperEta->AddEntry(fitTCMInvXSectionPi0,"TCM fit","l");
        legendXsectionPaperEta->Draw();
        
        TLatex *labelEnergyXSectionPaperEta= new TLatex(0.18, 0.03+0.05*5, collisionSystem2760GeV.Data());
        SetStyleTLatex( labelEnergyXSectionPaperEta, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 11);
        labelEnergyXSectionPaperEta->Draw();
        TLatex *labelALICEXSectionPaperEta= new TLatex(0.18,0.03+0.05*4,"ALICE");
        SetStyleTLatex( labelALICEXSectionPaperEta, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 11);
        labelALICEXSectionPaperEta->Draw();
        TLatex *labelDetSysXSectionPaperEta = new TLatex(0.18,0.03+0.055*3,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPaperEta, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 11);
        labelDetSysXSectionPaperEta->Draw();
        

        // labels upper right corner
        TLegend* legendXsectionPaperEtaTheo     = GetAndSetLegend2(0.59, 0.95-0.04*2, 0.59+0.33, 0.95, textSizeLabelsPixel, 1, "", 43, 0.2);
        legendXsectionPaperEtaTheo->AddEntry(histoPythia8InvXSectionEta,"PYTHIA 8.2,","l");
        legendXsectionPaperEtaTheo->AddEntry((TObject*)0,"Monash 2013","");
        legendXsectionPaperEtaTheo->Draw();

        TLatex *labelNLOHeaderEta= new TLatex(0.60, 0.805, "NLO, PDF: CTEQ6M5");
        SetStyleTLatex( labelNLOHeaderEta, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 11);
        labelNLOHeaderEta->Draw();
        TLatex *labelNLOHeaderEta2= new TLatex(0.695, 0.765, "FF: AESSS");
        SetStyleTLatex( labelNLOHeaderEta2, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 11);
        labelNLOHeaderEta2->Draw();
                
        TLegend* legendXsectionPaperEtaTheo2     = GetAndSetLegend2(0.73, 0.755-0.048*3, 0.95, 0.755, textSizeLabelsPixel, 1, "", 43, 0.27);
        legendXsectionPaperEtaTheo2->AddEntry(graphNLOCalcEtaMuHalf,"#mu = 0.5 #it{p}_{T}","l");
        legendXsectionPaperEtaTheo2->AddEntry(graphNLOCalcEtaMuOne,"#mu = #it{p}_{T}","l");
//         legendXsectionPaperEtaTheo2->AddEntry((TObject*)0,"","");
        legendXsectionPaperEtaTheo2->AddEntry(graphNLOCalcEtaMuTwo,"#mu = 2 #it{p}_{T}","l");
        legendXsectionPaperEtaTheo2->Draw();

        
    padInvSectionNLORatio->cd();
    padInvSectionNLORatio->SetLogx(1);
        TH2F * ratio2DNLOEta                = new TH2F("ratio2DNLOEta","ratio2DNLOEta",1000,0.33, 27.,1000,0.2,3.45);
        SetStyleHistoTH2ForGraphs(ratio2DNLOEta, "#it{p}_{T} (GeV/#it{c})","#frac{NLO, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle, 
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

        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
        graphRatioEtaCombNLOMuHalf->Draw("same,c");
        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
        graphRatioEtaCombNLOMuOne->Draw("same,c");
        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
        graphRatioEtaCombNLOMuTwo->Draw("same,c");
        
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatAWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioEtaCombCombFitStatAWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioEtaCombCombFitSysA->SetLineWidth(0);
        graphRatioEtaCombCombFitSysA->Draw("2,same");
        graphRatioEtaCombCombFitStatAWOXErr->Draw("p,same");
        
        boxErrorSigma2760GeVRatioEta->Draw();
        DrawGammaLines(0.33, 27.,1., 1.,0.1,kGray);
        
    padInvSectionPythiaRatio->cd();
    padInvSectionPythiaRatio->SetLogx(1);
        TH2F * ratio2DPythiaEta             = new TH2F("ratio2DPythiaEta","ratio2DPythiaEta",1000,0.33, 27.,1000,0.5,1.9);
        SetStyleHistoTH2ForGraphs(ratio2DPythiaEta, "#it{p}_{T} (GeV/#it{c})","#frac{PYTHIA, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown, 
                                  0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.2/(textsizeFacXSecDown*marginXSec), 510, 512);
        ratio2DPythiaEta->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DPythiaEta->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DPythiaEta->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DPythiaEta->GetXaxis()->SetLabelFont(42);
        ratio2DPythiaEta->GetYaxis()->SetLabelFont(42);
        ratio2DPythiaEta->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DPythiaEta->GetXaxis()->SetTickLength(0.06);
        ratio2DPythiaEta->GetYaxis()->SetTickLength(0.04);
        ratio2DPythiaEta->DrawCopy();

        graphRatioEtaPythia8ToFit->Draw("3,same");
        DrawGammaSetMarker(histoRatioEtaPythia8ToFit, 24, 1.5, kRed+2 , kRed+2);  
        histoRatioEtaPythia8ToFit->SetLineWidth(widthCommonFit);
        histoRatioEtaPythia8ToFit->Draw("same,hist,l");
        
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatAWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioEtaCombCombFitStatAWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioEtaCombCombFitSysA->SetLineWidth(0);
        graphRatioEtaCombCombFitSysA->Draw("2,same");
        graphRatioEtaCombCombFitStatAWOXErr->Draw("p,same");

        boxErrorSigma2760GeVRatioEta->Draw();
        DrawGammaLines(0.33, 27.,1., 1.,0.1,kGray);
        
    canvasInvSectionPaper->Print(Form("%s/Eta_InvXSectionWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));


    // ***************************************************************************************************************
    // ******************************* Plotting eta/pi0 ratio for single measurements ********************************
    // ***************************************************************************************************************
    textSizeLabelsPixel                 = 54;
    TCanvas* canvasEtatoPi0combo       = new TCanvas("canvasEtatoPi0combo","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasEtatoPi0combo, 0.095, 0.002, 0.01, 0.125);
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
    
    TH2F * histo2DEtatoPi0combo;
    histo2DEtatoPi0combo               = new TH2F("histo2DEtatoPi0combo","histo2DEtatoPi0combo",1000,0.33, 27.,1000,0.,1.05 );
    SetStyleHistoTH2ForGraphs(histo2DEtatoPi0combo, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.75*textsizeLabelsEtaToPi0, textsizeLabelsEtaToPi0, 
                              0.75*textsizeLabelsEtaToPi0,1.1*textsizeLabelsEtaToPi0, 0.9, 0.65, 510, 510);
    histo2DEtatoPi0combo->GetXaxis()->SetMoreLogLabels();
    histo2DEtatoPi0combo->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(-10,10);
    histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.0,1.05);
    histo2DEtatoPi0combo->Draw("copy");

        // plotting systematics graphs
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaToPi0Sys, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0], widthLinesBoxes, kTRUE);
        graphPCMEtaToPi0Sys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaToPi0Sys, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4], widthLinesBoxes, kTRUE);
        graphPCMEMCALEtaToPi0Sys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaToPi0Sys, markerStyleDet[2], markerSizeDet[4]*0.75, colorDet[2] , colorDet[2], widthLinesBoxes, kTRUE);
        graphEMCALEtaToPi0Sys->Draw("E2same");

        // plotting statistics graphs
        DrawGammaSetMarker(histoPCMEtaToPi0Stat, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
        histoPCMEtaToPi0Stat->Draw("p,same,e");
        DrawGammaSetMarker(histoPCMEMCALEtaToPi0Stat, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4]);
        histoPCMEMCALEtaToPi0Stat->Draw("p,same,e");
        DrawGammaSetMarker(histoEMCALEtaToPi0Stat, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
        histoEMCALEtaToPi0Stat->Draw("p,same,e");
    
        TLegend* legendEtaToPi0 = GetAndSetLegend2(0.67, 0.15, 0.9, 0.15+(textsizeLabelsEtaToPi0*3*0.9), textSizeLabelsPixel);
        legendEtaToPi0->AddEntry(graphPCMEtaToPi0Sys,nameMeasGlobalLabel[0],"pf");
        legendEtaToPi0->AddEntry(graphPCMEMCALEtaToPi0Sys,nameMeasGlobalLabel[4],"pf");
        legendEtaToPi0->AddEntry(graphEMCALEtaToPi0Sys,nameMeasGlobalLabel[2],"pf");
        legendEtaToPi0->Draw();

        TLatex *labelEnergyEtaToPi0 = new TLatex(0.13, 0.92,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelEnergyEtaToPi0, 0.85*textsizeLabelsEtaToPi0,4, 1, 42, kTRUE, 11);
        labelEnergyEtaToPi0->Draw();
        
        TLatex *labelALICEEtaToPi0 = new TLatex(0.13, 0.92-(1*textsizeLabelsEtaToPi0*0.85),"ALICE");
        SetStyleTLatex( labelALICEEtaToPi0, 0.85*textsizeLabelsEtaToPi0,4, 1, 42, kTRUE, 11);
        labelALICEEtaToPi0->Draw();
        
//         TLatex *labelPi0EtaToPi0 = new TLatex(0.13, 0.92-(2*textsizeLabelsEtaToPi0*0.9),"#eta/#pi^{0}");
//         SetStyleTLatex( labelPi0EtaToPi0, textsizeLabelsEtaToPi0,4, 1, 42, kTRUE, 11);
//         labelPi0EtaToPi0->Draw();

    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_differentSystems.%s",outputDir.Data(), suffix.Data()));

    // ***************************************************************************************************************
    // ******************************* Plotting eta/pi0 ratio for combined measurement *******************************
    // ***************************************************************************************************************    
    histo2DEtatoPi0combo->Draw("copy");

        TGraphAsymmErrors* graphCombEtaToPi0StatAWOXErr = (TGraphAsymmErrors*)graphCombEtaToPi0StatA->Clone("graphCombEtaToPi0StatAWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphCombEtaToPi0StatAWOXErr);
    
        // plotting data
        graphCombEtaToPi0StatA->Print();
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0StatAWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphCombEtaToPi0StatAWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0SysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphCombEtaToPi0SysA->SetLineWidth(0);
        graphCombEtaToPi0SysA->Draw("2,same");
        graphCombEtaToPi0StatAWOXErr->Draw("p,same");
        
        // plotting labels
        labelEnergyEtaToPi0->Draw();
        labelALICEEtaToPi0->Draw();
//         labelPi0EtaToPi0->Draw();
        
 
    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Paper.%s",outputDir.Data(), suffix.Data()));

    
    // ***************************************************************************************************************
    // ********************** Plotting eta/pi0 ratio for combined measurement + Theory *******************************
    // ***************************************************************************************************************    
    histo2DEtatoPi0combo->Draw("copy");
        
        DrawGammaSetMarkerTGraphErr(graphPythia8EtaToPi0, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2);
        graphPythia8EtaToPi0->Draw("3,same");
        DrawGammaSetMarker(histoPythia8EtaToPi0, 24, 1.5, kRed+2 , kRed+2);  
        histoPythia8EtaToPi0->SetLineWidth(widthCommonFit);
        histoPythia8EtaToPi0->Draw("same,hist,l");
        
        DrawGammaSetMarkerTGraphAsym(graphNLOEtaToPi0, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
        graphNLOEtaToPi0->Draw("3,same");

        
        graphCombEtaToPi0SysA->Draw("2,same");
        graphCombEtaToPi0StatAWOXErr->Draw("p,same");
        
        // plotting labels
        labelEnergyEtaToPi0->Draw();
        labelALICEEtaToPi0->Draw();
//         labelPi0EtaToPi0->Draw();

        TLegend* legendEtaToPi0Theory = GetAndSetLegend2(0.48, 0.95, 0.9, 0.95-(textsizeLabelsEtaToPi0*4*0.9), textSizeLabelsPixel*0.85, 1, "", 43, 0.16);
        legendEtaToPi0Theory->AddEntry(graphCombEtaToPi0SysA,"Data","pf");
        legendEtaToPi0Theory->AddEntry(graphNLOEtaToPi0,"NLO, PDF:CTEQ6M5 ","pf");
        legendEtaToPi0Theory->AddEntry((TObject*)0,"#pi^{0} FF: DSS07, #eta FF: AESSS","");
        legendEtaToPi0Theory->AddEntry(histoPythia8EtaToPi0,"PYTHIA 8.2, Monash 2013","l");
        legendEtaToPi0Theory->Draw();
        
    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Theory_Paper.%s",outputDir.Data(), suffix.Data()));

    // ***************************************************************************************************************
    // **************** Plotting eta/pi0 ratio for combined measurement + Theory + other energies ********************
    // ***************************************************************************************************************    
    canvasEtatoPi0combo->cd();
    histo2DEtatoPi0combo->Draw("copy");
        graphPythia8EtaToPi0->Draw("3,same");
        histoPythia8EtaToPi0->Draw("same,hist,l");
        graphNLOEtaToPi0->Draw("3,same");
        
        graphCombEtaToPi0SysA->Draw("2,same");
               
        DrawGammaSetMarkerTGraphErr(graphPHENIXEtaToPi0200GeV, 25, 2., kGray+2, kGray+2, widthLinesBoxes, kFALSE);
        graphPHENIXEtaToPi0200GeV->Draw("p,same");
        DrawGammaSetMarkerTGraphAsym(graphALICEEtaToPi07TeV, 24, 2.2, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphALICEEtaToPi07TeV->Draw("p,same");
        
        graphCombEtaToPi0StatAWOXErr->Draw("p,same");
        
        TLegend* legendEtaToPi0Theory2 = GetAndSetLegend2(0.13, 0.95, 0.46, 0.95-(textsizeLabelsEtaToPi0*4*0.9), textSizeLabelsPixel*0.85, 1, "", 43, 0.16);
        legendEtaToPi0Theory2->AddEntry(graphCombEtaToPi0SysA,Form("ALICE, %s",collisionSystem2760GeV.Data()),"pf");
        legendEtaToPi0Theory2->AddEntry(graphNLOEtaToPi0,"NLO, PDF:CTEQ6M5 ","pf");
        legendEtaToPi0Theory2->AddEntry((TObject*)0,"#pi^{0} FF: DSS07, #eta FF: AESSS","");
        legendEtaToPi0Theory2->AddEntry(histoPythia8EtaToPi0,"PYTHIA 8.2, Monash 2013","l");
        legendEtaToPi0Theory2->Draw();

        TLegend* legendEtaToPi0WorldData = GetAndSetLegend2(0.51, 0.155+(textsizeLabelsEtaToPi0*2*0.9), 0.9, 0.155, textSizeLabelsPixel*0.85, 1, "", 43, 0.16);
        legendEtaToPi0WorldData->AddEntry(graphALICEEtaToPi07TeV,"ALICE, pp, #sqrt{#it{s}} = 7 TeV","p");
        legendEtaToPi0WorldData->AddEntry(graphPHENIXEtaToPi0200GeV,"PHENIX, pp, #sqrt{#it{s}} = 0.2 TeV", "p");
        legendEtaToPi0WorldData->Draw();
        
    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Theory_WorldData_Paper.%s",outputDir.Data(), suffix.Data()));


    // ***************************************************************************************************************
    // **************** Plotting eta/pi0 ratio for combined measurement + other energies ********************
    // ***************************************************************************************************************        
    canvasEtatoPi0combo->cd();
    histo2DEtatoPi0combo->Draw("copy");
        
        graphCombEtaToPi0SysA->Draw("2,same");
               
        DrawGammaSetMarkerTGraphErr(graphPHENIXEtaToPi0200GeV, 25, 2., kGray+2, kGray+2, widthLinesBoxes, kFALSE);
        graphPHENIXEtaToPi0200GeV->Draw("p,same");
        DrawGammaSetMarkerTGraphAsym(graphALICEEtaToPi07TeV, 24, 2.2, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphALICEEtaToPi07TeV->Draw("p,same");
        
        graphCombEtaToPi0StatAWOXErr->Draw("p,same");
        
        TLegend* legendEtaToPi0WorldData2 = GetAndSetLegend2(0.13, 0.95, 0.46, 0.95-(textsizeLabelsEtaToPi0*3*0.9), textSizeLabelsPixel*0.85, 1, "", 43, 0.16);
        legendEtaToPi0WorldData2->AddEntry(graphCombEtaToPi0SysA,Form("ALICE, %s",collisionSystem2760GeV.Data()),"pf");
        legendEtaToPi0WorldData2->AddEntry(graphALICEEtaToPi07TeV,"ALICE, pp, #sqrt{#it{s}} = 7 TeV","p");
        legendEtaToPi0WorldData2->AddEntry(graphPHENIXEtaToPi0200GeV,"PHENIX, pp, #sqrt{#it{s}} = 0.2 TeV", "p");
        legendEtaToPi0WorldData2->Draw();
        
    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_WorldData_Paper.%s",outputDir.Data(), suffix.Data()));
    
    
    // ***************************************************************************************************************
    // ******************************** fitting eta/pi0 **************************************************************
    // ***************************************************************************************************************    
    TF1* etaToPi0ConstData  = new TF1("etaToPi0ConstData","[0]",4,20);
    TF1* etaToPi0ConstMC    = new TF1("etaToPi0ConstMC","[0]",4,20);
    graphCombEtaToPi0StatAWOXErr->Fit(etaToPi0ConstData,"QRME0","",4,20);
    histoPythia8EtaToPi0->Fit(etaToPi0ConstMC,"QRME0","",4,20);
    
    fileFitsOutput << "***********************************************************************************************************" << endl;
    fileFitsOutput << "***********************************************************************************************************" << endl;
    fileFitsOutput << "***********************************************************************************************************" << endl;
    fileFitsOutput << "high pt eta/pi0 - data, stat: " << etaToPi0ConstData->GetParameter(0) << "+-"<< etaToPi0ConstData->GetParError(0) << endl;
    fileFitsOutput << WriteParameterToFile(etaToPi0ConstData) << endl;
    graphCombEtaToPi0TotA->Fit(etaToPi0ConstData,"QRME0","",8,20);
    fileFitsOutput << "high pt eta/pi0 - data, tot: " << etaToPi0ConstData->GetParameter(0) << "+-"<< etaToPi0ConstData->GetParError(0) << endl;
    fileFitsOutput << WriteParameterToFile(etaToPi0ConstData) << endl;
    fileFitsOutput << "high pt eta/pi0 - pythia 8: " << etaToPi0ConstMC->GetParameter(0) << "+-"<< etaToPi0ConstMC->GetParError(0) << endl;
    fileFitsOutput << WriteParameterToFile(etaToPi0ConstMC) << endl;
    fileFitsOutput << "***********************************************************************************************************" << endl;
    fileFitsOutput << "***********************************************************************************************************" << endl;
    fileFitsOutput << "***********************************************************************************************************" << endl;
//     return;
    //*************************************************************************************************************
    //***************************** Comparison to Charged pions ***************************************************
    //*************************************************************************************************************
    TGraphAsymmErrors* graphCombPi0InvYieldTotA      = (TGraphAsymmErrors*) graphCombPi0InvXSectionTotA->Clone("graphCombPi0InvYieldTotA");
    graphCombPi0InvYieldTotA                         = ScaleGraph(graphCombPi0InvYieldTotA,1./xSection2760GeVINEL);
    TGraphAsymmErrors* graphCombPi0InvYieldStatA     = (TGraphAsymmErrors*) graphCombPi0InvXSectionStatA->Clone("graphCombPi0InvYieldStatA");
    graphCombPi0InvYieldStatA                        = ScaleGraph(graphCombPi0InvYieldStatA,1./xSection2760GeVINEL);
    TGraphAsymmErrors* graphCombPi0InvYieldSysA      = (TGraphAsymmErrors*) graphCombPi0InvXSectionSysA->Clone("graphCombPi0InvYieldSysA");
    graphCombPi0InvYieldSysA                         = ScaleGraph(graphCombPi0InvYieldSysA,1./xSection2760GeVINEL);
    TGraphAsymmErrors* graphCombEtaInvYieldTotA      = (TGraphAsymmErrors*) graphCombEtaInvXSectionTotA->Clone("graphCombEtaInvYieldTotA");
    graphCombEtaInvYieldTotA                         = ScaleGraph(graphCombEtaInvYieldTotA,1./xSection2760GeVINEL);
    TGraphAsymmErrors* graphCombEtaInvYieldStatA     = (TGraphAsymmErrors*) graphCombEtaInvXSectionStatA->Clone("graphCombEtaInvYieldStatA");
    graphCombEtaInvYieldStatA                        = ScaleGraph(graphCombEtaInvYieldStatA,1./xSection2760GeVINEL);
    TGraphAsymmErrors* graphCombEtaInvYieldSysA      = (TGraphAsymmErrors*) graphCombEtaInvXSectionSysA->Clone("graphCombEtaInvYieldSysA");
    graphCombEtaInvYieldSysA                         = ScaleGraph(graphCombEtaInvYieldSysA,1./xSection2760GeVINEL);

       cout << "combined Spectrum - high Pt" << endl;
    TGraphErrors* graphChPiInvYieldHighPtStatPPHighPtCombUp         = NULL;
    TGraphErrors* graphChPiInvYieldHighPtSystPPHighPtCombUp         = NULL;
    TGraphErrors* graphCombPi0InvYieldStatRebinnedHighPtComb = NULL;
    TGraphErrors* graphCombPi0InvYieldSysRebinnedHighPtComb  = NULL;
    TGraphErrors* graphRatioPi0HighPtChPisCombA = CalculateRatioBetweenSpectraWithDifferentBinning( graphCombPi0InvYieldStatA, graphCombPi0InvYieldSysA, 
                                                                                                        histoChPiInvYieldHighPtStatPP, histoChPiInvYieldHighPtSystPP,  
                                                                                                        kTRUE,  kTRUE, 
                                                                                                        &graphCombPi0InvYieldStatRebinnedHighPtComb, &graphCombPi0InvYieldSysRebinnedHighPtComb, 
                                                                                                        &graphChPiInvYieldHighPtStatPPHighPtCombUp, &graphChPiInvYieldHighPtSystPPHighPtCombUp )    ;
   
    cout << "combined Spectrum - low Pt" << endl;
    TGraphErrors* graphChPiInvYieldLowPtStatPPLowPtCombUp           = NULL;
    TGraphErrors* graphChPiInvYieldLowPtSystPPLowPtCombUp           = NULL;
    TGraphErrors* graphCombPi0InvYieldStatRebinnedLowPtComb  = NULL;
    TGraphErrors* graphCombPi0InvYieldSysRebinnedLowPtComb   = NULL;
    TGraphErrors* graphRatioPi0LowPtChPisCombA  = CalculateRatioBetweenSpectraWithDifferentBinning( graphCombPi0InvYieldStatA, graphCombPi0InvYieldSysA, 
                                                                                                        histoChPiInvYieldLowPtStatPP, histoChPiInvYieldLowPtSysPP,  
                                                                                                        kTRUE,  kTRUE, 
                                                                                                        &graphCombPi0InvYieldStatRebinnedLowPtComb, &graphCombPi0InvYieldSysRebinnedLowPtComb,
                                                                                                        &graphChPiInvYieldLowPtStatPPLowPtCombUp, &graphChPiInvYieldLowPtSystPPLowPtCombUp )    ;
    
    cout << "combined Spectrum - low Pt CMS" << endl;
   
    TGraphErrors* graphChPiInvYieldCMSStatPPCMSCombUp       = NULL;
    TGraphErrors* graphChPiInvYieldCMSSystPPCMSCombUp       = NULL;
    TGraphErrors* graphCombPi0InvYieldStatRebinnedCMSComb   = NULL;
    TGraphErrors* graphCombPi0InvYieldSysRebinnedCMSComb    = NULL;
    TGraphErrors* graphRatioPi0CMSChPisCombA    = CalculateRatioBetweenSpectraWithDifferentBinning( graphCombPi0InvYieldStatA, graphCombPi0InvYieldSysA, 
                                                                                                        histoChPiInvYieldLowPtStatCMS, histoChPiInvYieldLowPtSysCMS,  
                                                                                                        kTRUE,  kTRUE, 
                                                                                                        &graphCombPi0InvYieldStatRebinnedCMSComb, &graphCombPi0InvYieldSysRebinnedCMSComb,
                                                                                                        &graphChPiInvYieldCMSStatPPCMSCombUp, &graphChPiInvYieldCMSSystPPCMSCombUp ) ;

       cout << "combined Spectrum - published" << endl;
    TGraphErrors* graphChPiInvYieldPubStatPPPubCombUp       = NULL;
    TGraphErrors* graphChPiInvYieldPubSystPPPubCombUp       = NULL;
    TGraphErrors* graphCombPi0InvYieldStatRebinnedPubComb   = NULL;
    TGraphErrors* graphCombPi0InvYieldSysRebinnedPubComb    = NULL;
    TGraphErrors* graphRatioPi0PubChPisCombA    = CalculateRatioBetweenSpectraWithDifferentBinning( graphCombPi0InvYieldStatA, graphCombPi0InvYieldSysA, 
                                                                                                        histoChPiInvYieldPubStatPP, histoChPiInvYieldPubSystPP,  
                                                                                                        kTRUE,  kTRUE, 
                                                                                                        &graphCombPi0InvYieldStatRebinnedPubComb, &graphCombPi0InvYieldSysRebinnedPubComb, 
                                                                                                        &graphChPiInvYieldPubStatPPPubCombUp, &graphChPiInvYieldPubSystPPPubCombUp )    ;

       cout << "combined Spectrum - published" << endl;
    TGraphErrors* graphChHadInvYieldPubStatPPHadCombUp      = NULL;
    TGraphErrors* graphChHadInvYieldPubSystPPHadCombUp      = NULL;
    TGraphErrors* graphCombPi0InvYieldStatRebinnedHadComb   = NULL;
    TGraphErrors* graphCombPi0InvYieldSysRebinnedHadComb    = NULL;
    TGraphErrors* graphRatioPi0PubChHadsCombA   = CalculateRatioBetweenSpectraWithDifferentBinning( graphCombPi0InvYieldStatA, graphCombPi0InvYieldSysA, 
                                                                                                    graphChargedHadronsStatPP, graphChargedHadronsSysPP,  
                                                                                                    kTRUE,  kTRUE, 
                                                                                                    &graphCombPi0InvYieldStatRebinnedHadComb, &graphCombPi0InvYieldSysRebinnedHadComb, 
                                                                                                    &graphChHadInvYieldPubStatPPHadCombUp, &graphChHadInvYieldPubSystPPHadCombUp )    ;

    TF1* fitPowInvYieldChHadALICE       = FitObject("powPure","fitPowInvYieldChHadALICE","pi0",graphChargedHadronsStatPP,7,50. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvYieldChHadALICE)<< endl;                       
    fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldChHadALICE)<< endl;    
    TF1* fitPowInvYieldChHadATLAS       = FitObject("powPure","fitPowInvYieldChHadATLAS","pi0",graphChargedHadronsStatPPATLAS,7,70. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvYieldChHadATLAS)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldChHadATLAS)<< endl;    
    TF1* fitPowInvYieldChHadCMS         = FitObject("powPure","fitPowInvYieldChHadCMS","pi0",graphChargedHadronsStatPPCMS,7,70. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvYieldChHadCMS)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldChHadCMS)<< endl;    

    TGraphAsymmErrors* graphRatioChHadStatToChFitALICE      = CalculateGraphErrRatioToFit(graphChargedHadronsStatPP, fitPowInvYieldChHadALICE);
    TGraphAsymmErrors* graphRatioChHadSysToChFitALICE       = CalculateGraphErrRatioToFit(graphChargedHadronsSysPP, fitPowInvYieldChHadALICE);
    TGraphAsymmErrors* graphRatioPi0StatToChFitALICE        = CalculateGraphErrRatioToFit(graphCombPi0InvYieldStatA, fitPowInvYieldChHadALICE);
    TGraphAsymmErrors* graphRatioPi0SysToChFitALICE         = CalculateGraphErrRatioToFit(graphCombPi0InvYieldSysA, fitPowInvYieldChHadALICE);
    while (graphRatioChHadStatToChFitALICE->GetX()[0] < 7)
        graphRatioChHadStatToChFitALICE->RemovePoint(0);   
    while (graphRatioChHadSysToChFitALICE->GetX()[0] < 7)    
        graphRatioChHadSysToChFitALICE->RemovePoint(0);   
    while (graphRatioPi0StatToChFitALICE->GetX()[0] < 7)
        graphRatioPi0StatToChFitALICE->RemovePoint(0);   
    while (graphRatioPi0SysToChFitALICE->GetX()[0] < 7)
        graphRatioPi0SysToChFitALICE->RemovePoint(0);   
    
    
    TGraphAsymmErrors* graphRatioChHadStatToChFitATLAS      = CalculateGraphErrRatioToFit(graphChargedHadronsStatPPATLAS, fitPowInvYieldChHadATLAS);
    TGraphAsymmErrors* graphRatioChHadSysToChFitATLAS       = CalculateGraphErrRatioToFit(graphChargedHadronsSysPPATLAS, fitPowInvYieldChHadATLAS);
    TGraphAsymmErrors* graphRatioPi0StatToChFitATLAS        = CalculateGraphErrRatioToFit(graphCombPi0InvYieldStatA, fitPowInvYieldChHadATLAS);
    TGraphAsymmErrors* graphRatioPi0SysToChFitATLAS         = CalculateGraphErrRatioToFit(graphCombPi0InvYieldSysA, fitPowInvYieldChHadATLAS);
    while (graphRatioChHadStatToChFitATLAS->GetX()[0] < 7)
        graphRatioChHadStatToChFitATLAS->RemovePoint(0);  
    while (graphRatioChHadSysToChFitATLAS->GetX()[0] < 7)
        graphRatioChHadSysToChFitATLAS->RemovePoint(0);   
    while (graphRatioPi0StatToChFitATLAS->GetX()[0] < 7)
        graphRatioPi0StatToChFitATLAS->RemovePoint(0);   
    while (graphRatioPi0SysToChFitATLAS->GetX()[0] < 7)
        graphRatioPi0SysToChFitATLAS->RemovePoint(0);   
    
    
    TGraphAsymmErrors* graphRatioChHadStatToChFitCMS        = CalculateGraphErrRatioToFit(graphChargedHadronsStatPPCMS, fitPowInvYieldChHadCMS);
    TGraphAsymmErrors* graphRatioChHadSysToChFitCMS         = CalculateGraphErrRatioToFit(graphChargedHadronsSysPPCMS, fitPowInvYieldChHadCMS);
    TGraphAsymmErrors* graphRatioPi0StatToChFitCMS          = CalculateGraphErrRatioToFit(graphCombPi0InvYieldStatA, fitPowInvYieldChHadCMS);
    TGraphAsymmErrors* graphRatioPi0SysToChFitCMS           = CalculateGraphErrRatioToFit(graphCombPi0InvYieldSysA, fitPowInvYieldChHadCMS);
    while (graphRatioChHadStatToChFitCMS->GetX()[0] < 7)
        graphRatioChHadStatToChFitCMS->RemovePoint(0);   
    while (graphRatioChHadSysToChFitCMS->GetX()[0] < 7)
        graphRatioChHadSysToChFitCMS->RemovePoint(0);   
    while (graphRatioPi0StatToChFitCMS->GetX()[0] < 7)
        graphRatioPi0StatToChFitCMS->RemovePoint(0);   
    while (graphRatioPi0SysToChFitCMS->GetX()[0] < 7)
        graphRatioPi0SysToChFitCMS->RemovePoint(0);   
     
    

    // ***************************************************************************************************************
    // ************************** Comparison pi0/pi+-, pi0 updated comb pp 2.76TeV ***********************
    // ***************************************************************************************************************    
    textSizeLabelsPixel             = 48;
    TCanvas* canvasCompYieldPPInd   = new TCanvas("canvasCompYieldPPInd","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasCompYieldPPInd,   0.12, 0.01, 0.01, 0.11);
    canvasCompYieldPPInd->SetLogx();

    TH2F * histo2DCompCombinedRatio2;
    histo2DCompCombinedRatio2       = new TH2F("histo2DCompCombinedRatio2","histo2DCompCombinedRatio2",1000,0.23,70.,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio2, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
    histo2DCompCombinedRatio2->GetYaxis()->SetNoExponent(kTRUE);
//  histo2DCompCombinedRatio2->GetXaxis()->SetLabelOffset(-0.01);
    histo2DCompCombinedRatio2->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DCompCombinedRatio2->GetXaxis()->SetNoExponent(kTRUE);
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.05,2.45);
    histo2DCompCombinedRatio2->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioPi0HighPtChPisCombA, markerStyleCombHighPt, markerSizeComparison, colorCombHighPt , colorCombHighPt);
        graphRatioPi0HighPtChPisCombA->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioPi0LowPtChPisCombA, markerStyleCombLowPt, markerSizeComparison, colorCombLowPt , colorCombLowPt);
        graphRatioPi0LowPtChPisCombA->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioPi0CMSChPisCombA, 22, markerSizeComparison, kRed+1 , kRed+1);
        graphRatioPi0CMSChPisCombA->Draw("E1psame");

        TLegend* legendPi0CompChargedPionsPP    = GetAndSetLegend2(0.15, 0.8, 0.9, 0.90, 0.85* textSizeLabelsPixel, 2, "", 43, 0.12);
        legendPi0CompChargedPionsPP->AddEntry(graphRatioPi0LowPtChPisCombA,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (ALICE)","p");
        legendPi0CompChargedPionsPP->AddEntry(graphRatioPi0HighPtChPisCombA,"#pi^{0}/#pi^{#pm} high #it{p}_{T} (ALICE)","p");
        legendPi0CompChargedPionsPP->AddEntry(graphRatioPi0CMSChPisCombA,"#pi^{0}/#pi^{#pm} low #it{p}_{T} (#pi^{#pm} from CMS)","p");
        legendPi0CompChargedPionsPP->Draw();
    
        labelRatioTheoryPP->Draw();
        
        DrawGammaLines(0., 70 , 1, 1 ,1, kGray, 1);   
   
    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedToNeutralOldCharged_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    canvasCompYieldPPInd->cd();

    histo2DCompCombinedRatio2->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioPi0PubChPisCombA, markerStyleCombHighPt, markerSizeComparison, kBlue+2 , kBlue+2);
        graphRatioPi0PubChPisCombA->Draw("E1psame");
        DrawGammaSetMarkerTGraphErr(graphRatioPi0CMSChPisCombA, 22, markerSizeComparison, kRed+1 , kRed+1);
        graphRatioPi0CMSChPisCombA->Draw("E1psame");

        graphRatioPi0PubChHadsCombA->Print();
        TLegend* legendPi0CompChargedPionsPP2   = GetAndSetLegend2(0.15, 0.86, 0.9, 0.94, 0.85* textSizeLabelsPixel,2, "", 43, 0.12);
        legendPi0CompChargedPionsPP2->AddEntry(graphRatioPi0PubChPisCombA,"#pi^{0}/#pi^{#pm}  (ALICE)","p");
        legendPi0CompChargedPionsPP2->AddEntry(graphRatioPi0CMSChPisCombA,"#pi^{0}/#pi^{#pm}  (#pi^{#pm} from CMS)","p");
        legendPi0CompChargedPionsPP2->Draw();
    
        labelRatioTheoryPP->Draw();
        
        DrawGammaLines(0., 70 , 1, 1 ,1, kGray, 1);   
   
    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedToNeutralPublishedCharged_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    canvasCompYieldPPInd->cd();

    histo2DCompCombinedRatio2->GetYaxis()->SetTitle("#pi^{0}/h^{#pm}");
    histo2DCompCombinedRatio2->GetYaxis()->SetRangeUser(0.05,1.35);    
    histo2DCompCombinedRatio2->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioPi0PubChHadsCombA, markerStyleCombHighPt, markerSizeComparison, kBlack , kBlack);
        graphRatioPi0PubChHadsCombA->Draw("E1psame");

        graphRatioPi0PubChHadsCombA->Print();
        TLegend* legendPi0CompChargedPionsPP3   = GetAndSetLegend2(0.15, 0.86, 0.9, 0.94, 0.85* textSizeLabelsPixel,2, "", 43, 0.12);
        legendPi0CompChargedPionsPP3->AddEntry(graphRatioPi0PubChHadsCombA,"#pi^{0}/h^{#pm}  (ALICE)","p");
//      legendPi0CompChargedPionsPP3->AddEntry(graphRatioPi0CMSChPisCombA,"#pi^{0}/#pi^{#pm}  (#pi^{#pm} from CMS)","p");
        legendPi0CompChargedPionsPP3->Draw();
    
        labelRatioTheoryPP->Draw();
        
        DrawGammaLines(0., 70 , 1, 1 ,1, kGray, 1);   
   
    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedHadronToNeutralPublishedCharged_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    histo2DCompCombinedRatio2->GetYaxis()->SetTitle("#pi^{0}/fit to h^{#pm}");
    histo2DCompCombinedRatio2->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0SysToChFitALICE, markerStyleCombHighPt, markerSizeComparison, kBlack , kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0SysToChFitALICE->Draw("E2psame");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0SysToChFitATLAS, markerStyleCombHighPt, markerSizeComparison, kBlue+1 , kBlue+1, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0SysToChFitATLAS->Draw("E2psame");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0SysToChFitCMS, markerStyleCombHighPt, markerSizeComparison, kRed+1 , kRed+1, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0SysToChFitCMS->Draw("E2psame");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0StatToChFitALICE, markerStyleCombHighPt, markerSizeComparison, kBlack , kBlack);
        graphRatioPi0StatToChFitALICE->Draw("E1psame");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0StatToChFitATLAS, markerStyleCombHighPt, markerSizeComparison, kBlue+1 , kBlue+1);
        graphRatioPi0StatToChFitATLAS->Draw("E1psame");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0StatToChFitCMS, markerStyleCombHighPt, markerSizeComparison, kRed+1 , kRed+1);
        graphRatioPi0StatToChFitCMS->Draw("E1psame");
        
        TLegend* legendPi0CompChargedPionsToFits   = GetAndSetLegend2(0.15, 0.85, 0.6, 0.92, 0.85* textSizeLabelsPixel, 3, "", 43, 0.25);
        legendPi0CompChargedPionsToFits->AddEntry(graphRatioPi0SysToChFitALICE,"ALICE","fp");
        legendPi0CompChargedPionsToFits->AddEntry(graphRatioPi0SysToChFitATLAS,"ATLAS","fp");
        legendPi0CompChargedPionsToFits->AddEntry(graphRatioPi0SysToChFitCMS,"CMS","fp");
        legendPi0CompChargedPionsToFits->Draw();
    
        labelRatioTheoryPP->Draw();
        
        DrawGammaLines(0., 70 , 1, 1 ,1, kGray, 1);   
   
    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedHadronToNeutralPublishedChargedOtherExp_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    
    TH2F * histo2DCompCombinedRatio3;
    histo2DCompCombinedRatio3       = new TH2F("histo2DCompCombinedRatio3","histo2DCompCombinedRatio3",1000,0.23,200.,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio3, "#it{p}_{T} (GeV/#it{c})","h^{#pm}/powerlaw fit to spec", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
    histo2DCompCombinedRatio3->GetYaxis()->SetNoExponent(kTRUE);
//  histo2DCompCombinedRatio3->GetXaxis()->SetLabelOffset(-0.01);
    histo2DCompCombinedRatio3->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DCompCombinedRatio3->GetXaxis()->SetNoExponent(kTRUE);    
    histo2DCompCombinedRatio3->GetYaxis()->SetTitle("");   
    histo2DCompCombinedRatio3->GetYaxis()->SetRangeUser(0.65,1.35);
    histo2DCompCombinedRatio3->GetXaxis()->SetRangeUser(5,170.);
    histo2DCompCombinedRatio3->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphRatioChHadSysToChFitALICE, markerStyleCombHighPt, markerSizeComparison, kBlack , kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioChHadSysToChFitALICE->Draw("E2psame");
        DrawGammaSetMarkerTGraphAsym(graphRatioChHadSysToChFitATLAS, markerStyleCombHighPt, markerSizeComparison, kBlue+1 , kBlue+1, widthLinesBoxes, kTRUE, 0);
        graphRatioChHadSysToChFitATLAS->Draw("E2psame");
        DrawGammaSetMarkerTGraphAsym(graphRatioChHadSysToChFitCMS, markerStyleCombHighPt, markerSizeComparison, kRed+1 , kRed+1, widthLinesBoxes, kTRUE, 0);
        graphRatioChHadSysToChFitCMS->Draw("E2psame");
        DrawGammaSetMarkerTGraphAsym(graphRatioChHadStatToChFitALICE, markerStyleCombHighPt, markerSizeComparison, kBlack , kBlack);
        graphRatioChHadStatToChFitALICE->Draw("E1psame");
        DrawGammaSetMarkerTGraphAsym(graphRatioChHadStatToChFitATLAS, markerStyleCombHighPt, markerSizeComparison, kBlue+1 , kBlue+1);
        graphRatioChHadStatToChFitATLAS->Draw("E1psame");
        DrawGammaSetMarkerTGraphAsym(graphRatioChHadStatToChFitCMS, markerStyleCombHighPt, markerSizeComparison, kRed+1 , kRed+1);
        graphRatioChHadStatToChFitCMS->Draw("E1psame");
        
        TLegend* legendPi0CompChargedPionsFits   = GetAndSetLegend2(0.15, 0.85, 0.6, 0.92, 0.85* textSizeLabelsPixel,3, "", 43, 0.25);
        legendPi0CompChargedPionsFits->AddEntry(graphRatioChHadSysToChFitALICE,"ALICE","fp");
        legendPi0CompChargedPionsFits->AddEntry(graphRatioChHadSysToChFitATLAS,"ATLAS","fp");
        legendPi0CompChargedPionsFits->AddEntry(graphRatioChHadSysToChFitCMS,"CMS","fp");
        legendPi0CompChargedPionsFits->Draw();
    
        labelRatioTheoryPP->Draw();
        
        DrawGammaLines(5., 170 , 1, 1 ,1, kGray, 1);   
   
    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedHadronToFit_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));


    canvasXSectionPi0->cd();
    TH2F * histo2DXSectionWithEtaAndPi0;
    histo2DXSectionWithEtaAndPi0          = new TH2F("histo2DXSectionWithEtaAndPi0","histo2DXSectionWithEtaAndPi0",11000,0.23,70.,1000,2e-3,10e11);
    SetStyleHistoTH2ForGraphs(histo2DXSectionWithEtaAndPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.45);
    histo2DXSectionWithEtaAndPi0->GetXaxis()->SetMoreLogLabels();
    histo2DXSectionWithEtaAndPi0->GetXaxis()->SetLabelOffset(-0.01);
    histo2DXSectionWithEtaAndPi0->Draw("copy");

        // scale eta graphs
        Double_t scaleFacEtaForCombPlot                             = 1e-2;
        TGraphAsymmErrors* graphCombEtaInvXSectionStatAWOXErrCopy   = (TGraphAsymmErrors*) graphCombEtaInvXSectionStatAWOXErr->Clone("graphCombEtaInvXSectionStatAWOXErrCopy");
        TGraphAsymmErrors* graphCombEtaInvXSectionSysACopy          = (TGraphAsymmErrors*) graphCombEtaInvXSectionSysA->Clone("graphCombEtaInvXSectionSysACopy");
        graphCombEtaInvXSectionStatAWOXErrCopy                      = ScaleGraph(graphCombEtaInvXSectionStatAWOXErrCopy,scaleFacEtaForCombPlot);
        graphCombEtaInvXSectionSysACopy                             = ScaleGraph(graphCombEtaInvXSectionSysACopy,scaleFacEtaForCombPlot);
        TH1D* histFitTCMInvXSectionEta                              = (TH1D*)fitTCMInvXSectionEta->GetHistogram();
        histFitTCMInvXSectionEta->Scale(scaleFacEtaForCombPlot);
        histoPythia8InvXSectionEta->Scale(scaleFacEtaForCombPlot);
        TGraphErrors* graphPythia8InvXSectionEtaScaled              = ScaleGraph(graphPythia8InvXSectionEta,scaleFacEtaForCombPlot);
        
        TGraphAsymmErrors* graphNLOEtaAESSCalcCopy                  = (TGraphAsymmErrors*)graphNLOEtaAESSCalc->Clone("graphNLOEtaAESSCalcCopy");
        TGraph* graphNLOCalcEtaMuHalfCopy                           = (TGraph*)graphNLOCalcEtaMuHalf->Clone("graphNLOCalcEtaMuHalfCopy");
        TGraph* graphNLOCalcEtaMuOneCopy                            = (TGraph*)graphNLOCalcEtaMuOne->Clone("graphNLOCalcEtaMuOneCopy");
        TGraph* graphNLOCalcEtaMuTwoCopy                            = (TGraph*)graphNLOCalcEtaMuTwo->Clone("graphNLOCalcEtaMuTwoCopy");
        graphNLOEtaAESSCalcCopy                                     = ScaleGraph(graphNLOEtaAESSCalcCopy,scaleFacEtaForCombPlot);
        graphNLOCalcEtaMuHalfCopy                                   = ScaleGraph(graphNLOCalcEtaMuHalfCopy,scaleFacEtaForCombPlot);
        graphNLOCalcEtaMuOneCopy                                    = ScaleGraph(graphNLOCalcEtaMuOneCopy,scaleFacEtaForCombPlot);
        graphNLOCalcEtaMuTwoCopy                                    = ScaleGraph(graphNLOCalcEtaMuTwoCopy,scaleFacEtaForCombPlot);
        
        // plotting Pythia 8.2 Monash
        graphPythia8InvXSectionPi0->Draw("3,same");
        DrawGammaSetMarker(histoPythia8InvXSectionPi0, 24, 1.5, kRed+2 , kRed+2);  
        histoPythia8InvXSectionPi0->SetLineWidth(widthCommonFit);
        histoPythia8InvXSectionPi0->Draw("same,hist,l");
        DrawGammaSetMarkerTGraphErr(graphPythia8InvXSectionEtaScaled, 0, 0, kRed+2 , kRed+2, widthLinesBoxes, kTRUE, kRed+2 );
        graphPythia8InvXSectionEtaScaled->Draw("3,same");
        DrawGammaSetMarker(histoPythia8InvXSectionEta, 24, 1.5, kRed+2 , kRed+2);  
        histoPythia8InvXSectionEta->SetLineWidth(widthCommonFit);
        histoPythia8InvXSectionEta->Draw("same,hist,l");

        // plotting NLO calcs pi0
        DrawGammaSetMarkerTGraphAsym(graphNLODSS14Calc, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
        graphNLODSS14Calc->Draw("3,same");
        
        // plotting NLO calcs eta
        DrawGammaSetMarkerTGraphAsym(graphNLOEtaAESSCalcCopy, 0, 0, colorCGC, colorCGC, widthLinesBoxes, kTRUE, colorCGC);
        graphNLOEtaAESSCalcCopy->Draw("3,same");       

        // plot data
        graphCombPi0InvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSysACopy, markerStyleComb+4, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
        graphCombEtaInvXSectionSysACopy->Draw("E2same");
        
        graphCombPi0InvXSectionStatAWOXErr->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStatAWOXErrCopy, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
        graphCombEtaInvXSectionStatAWOXErrCopy->Draw("p,same,z");
        
        // plots fits
        fitTCMInvXSectionPi0->Draw("same");                
        SetStyleHisto(histFitTCMInvXSectionEta, 2, 7, kGray+2);
        histFitTCMInvXSectionEta->Draw("same,c");
                        
        // labels lower left corner
        TLegend* legendXsectionPaperAll    = GetAndSetLegend2(0.17, 0.11, 0.5, 0.11+0.04*3, textSizeLabelsPixel, 1, "", 43, 0.2);
        legendXsectionPaperAll->AddEntry(graphCombPi0InvXSectionSysA,"#pi^{0}","pf");
        legendXsectionPaperAll->AddEntry(graphCombEtaInvXSectionSysACopy,"#eta (x 10^{-2})","pf");
//         legendXsectionPaperAll->AddEntry(dummyNorm,"Norm. unc. 2.5%","f");
        legendXsectionPaperAll->AddEntry(fitTCMInvXSectionPi0,"TCM fit","l");
        legendXsectionPaperAll->Draw();
//         TLegend* legendXsectionPaper3     = GetAndSetLegend2(0.17, 0.05+0.08, 0.5, 0.11+0.08, textSizeLabelsPixel, 1, "", 43, 0.2);
//         legendXsectionPaper3->AddEntry(fitTCMInvXSectionPi0,"#it{A}_{e} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{-n}","l");
//         legendXsectionPaper3->Draw();
        
        TLatex *labelEnergyXSectionPaperAll = new TLatex(0.18, 0.11+0.04*5, collisionSystem2760GeV.Data());
        SetStyleTLatex( labelEnergyXSectionPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelEnergyXSectionPaperAll->Draw();
        TLatex *labelALICEXSectionPaperAll  = new TLatex(0.18,0.11+0.04*4,"ALICE");
        SetStyleTLatex( labelALICEXSectionPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelALICEXSectionPaperAll->Draw();
        TLatex *labelALICENormUnPaperAll    = new TLatex(0.18,0.11+0.04*3+0.007,"Norm. unc. 2.5%");
        SetStyleTLatex( labelALICENormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelALICENormUnPaperAll->Draw();
        

        // labels upper right corner
        TLegend* legendXsectionPaperPyBoth  = GetAndSetLegend2(0.54, 0.95-0.04*4-0.04*0.75*2, 0.59+0.33, 0.95, textSizeLabelsPixel, 1, "", 43, 0.18);
        legendXsectionPaperPyBoth->AddEntry(histoPythia8InvXSectionEta,"PYTHIA 8.2","l");
        legendXsectionPaperPyBoth->AddEntry((TObject*)0,"Monash 2013","");
        legendXsectionPaperPyBoth->AddEntry(graphNLODSS14Calc,"#pi^{0} pQCD NLO ","f");
        legendXsectionPaperPyBoth->AddEntry((TObject*)0,"#scale[0.75]{PDF: MSTW, FF: DSS14}","");
        legendXsectionPaperPyBoth->AddEntry(graphNLOEtaAESSCalcCopy,"#eta pQCD NLO ","f");        
        legendXsectionPaperPyBoth->AddEntry((TObject*)0,"#scale[0.75]{PDF: CTEQ6M5, FF: AESSS}","");
        legendXsectionPaperPyBoth->Draw();

    canvasXSectionPi0->SaveAs(Form("%s/InvXSection_Pi0_Eta_Theory.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // *********************************** HFE results for 2.76TeV **********************************************************
    // **********************************************************************************************************************
    
    Double_t xValueHFE[]    = { 0.55,   0.65,   0.75,   0.85,   0.95,   1.05,   1.15,   1.25,   1.35,   1.45,
                                1.625,  1.875,  2.125,  2.375,  2.625,  2.875,  3.25,   3.75,   4.25,   4.75, 
                                5.5,    6.5,    7.5,    9.0,    11.0 };
    Double_t xErrNegHFE[]   = { 0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050, 
                                0.125,  0.125,  0.125,  0.125,  0.125,  0.125,  0.25,   0.25,   0.25,   0.25, 
                                0.5,    0.5,    0.5,    1.0,    1.0 };
    Double_t xErrPosHFE[]   = { 0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050, 
                                0.125,  0.125,  0.125,  0.125,  0.125,  0.125,  0.25,   0.25,   0.25,   0.25, 
                                0.5, 0.5, 0.5, 1.0, 1.0 };
    Double_t yValueHFE[]    = { 19.58,  8.902,  7.438,  5.914,  3.64,   2.333,  1.694,  1.08,   1.237,  0.7166, 
                                0.5536, 0.2903, 0.1756, 0.1035, 0.06608, 0.04384, 0.02174, 0.01383, 0.006796, 0.003306, 
                                0.001664, 7.53E-4, 3.84E-4, 1.6E-4, 4.8E-5 };
    Double_t yErrNegHFE[]   = { 14.00601299442493,      7.11655197409532,       4.035167902330707,      2.3827179858304675,     1.4889949630539387,
                                0.9269843580125827,     0.6120661728930948,     0.41292251089036064,    0.3157483174935379,     0.2211634011313807, 
                                0.12102098991497301,    0.06001666435249463,    0.030992418427738096,   0.019112822920751397,   0.012453513560437472, 
                                0.008021396387163522,   0.0042004285495649135,  0.002537892826736385,   0.0014131514426981985,  8.054321821233617E-4,
                                3.703093301552095E-4,   1.70578427709954E-4,    1.038315944209661E-4,   3.935733730830885E-5,   1.3892443989449805E-5 };
    Double_t yErrPosHFE[]   = { 14.194594041394774,     7.2290926124929396,     4.102352983349921,      2.4218804677357633,     1.5141710603495235, 
                                0.9423226623614653,     0.6225881463696525,     0.41949016675006817,    0.32020618357552055,    0.22401950807909565, 
                                0.12270631605585754,    0.06073483349775481,    0.0313517144666763,     0.019277447963877377,   0.012551099553425588, 
                                0.008072329279706076,   0.00422190715198712,    0.0025442091109026395,  0.0014165313268685586,  8.073295485735673E-4, 
                                3.7126001669988653E-4,  1.70578427709954E-4,    1.038315944209661E-4,   3.935733730830885E-5,   1.3892443989449805E-5 };
    Double_t yErrStatNegHFE[]   = { 1.72,       1.144,      0.628,      0.489,      0.375, 
                                    0.266,      0.18,       0.144,      0.144,      0.1115, 
                                    0.0428,     0.0266,     0.0137,     0.0109,     0.00826, 
                                    0.00552,    0.00294,    0.00197,    0.001169,   2.55E-4, 
                                    1.15E-4,    6.9E-5,     5.0E-5,     1.8E-5,     7.0E-6 };
    Double_t yErrStatPosHFE[]   = { 1.72,       1.144,      0.628,      0.489,      0.375, 
                                    0.266,      0.18,       0.144,      0.144,      0.1115, 
                                    0.0428,     0.0266,     0.0137,     0.0109,     0.00826, 
                                    0.00552,    0.00294,    0.00197,    0.001169,   2.55E-4, 
                                    1.15E-4,    6.9E-5,     5.0E-5,     1.8E-5,     7.0E-6 };
    Int_t nPointsHFE            = 25;
    TGraphAsymmErrors* graphHFE = new TGraphAsymmErrors(nPointsHFE, xValueHFE, yValueHFE, xErrNegHFE, xErrPosHFE, yErrNegHFE, yErrPosHFE);    
    graphHFE                    = ScaleGraph(graphHFE, 1e6);
    TF1* fitPowInvXSectionHFE   = FitObject("powPure","fitPowInvXSectionHFE2760GeV","HFE",graphHFE,3,10. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvXSectionHFE)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPowInvXSectionHFE)<< endl;    

    
    fileFitsOutput<< "Fitting powerlaw to pi0 6-40" << endl;
    TF1* fitPowInvXSectionPi0Tot   = FitObject("powPure","fitPowInvXSectionPi02760GeVTot","Pi0",graphCombPi0InvXSectionTotA,6,40. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvXSectionPi0Tot)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPowInvXSectionPi0Tot)<< endl;    
    TF1* fitPowInvXSectionPi0Stat   = FitObject("powPure","fitPowInvXSectionPi02760GeV","Pi0",graphCombPi0InvXSectionStatA,6,40. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvXSectionPi0Stat)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPowInvXSectionPi0Stat)<< endl;    

    fileFitsOutput<< "Fitting powerlaw to pi0 4-40" << endl;
    TF1* fitPowInvXSectionPi0TotLower   = FitObject("powPure","fitPowInvXSectionPi02760GeVTotLower","Pi0",graphCombPi0InvXSectionTotA,4,40. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvXSectionPi0TotLower)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPowInvXSectionPi0TotLower)<< endl;    
    TF1* fitPowInvXSectionPi0StatLower   = FitObject("powPure","fitPowInvXSectionPi02760GeVStatLower","Pi0",graphCombPi0InvXSectionStatA,4,40. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvXSectionPi0StatLower)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPowInvXSectionPi0StatLower)<< endl;    

    fileFitsOutput<< "Fitting powerlaw to eta 6-20" << endl;
    TF1* fitPowInvXSectionEtaTot   = FitObject("powPure","fitPowInvXSectionEta2760GeVTot","Eta",graphCombEtaInvXSectionTotA,6,20. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvXSectionEtaTot)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPowInvXSectionEtaTot)<< endl;    
    TF1* fitPowInvXSectionEtaStat   = FitObject("powPure","fitPowInvXSectionEta2760GeVStat","Eta",graphCombEtaInvXSectionStatA,6,20. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvXSectionEtaStat)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPowInvXSectionEtaStat)<< endl;    

    fileFitsOutput<< "Fitting powerlaw to eta 4-20" << endl;
    TF1* fitPowInvXSectionEtaTotLower   = FitObject("powPure","fitPowInvXSectionEta2760GeVTotLower","Eta",graphCombEtaInvXSectionTotA,4,20. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvXSectionEtaTotLower)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPowInvXSectionEtaTotLower)<< endl;    
    TF1* fitPowInvXSectionEtaStatLower   = FitObject("powPure","fitPowInvXSectionEta2760GeVStatLower","Eta",graphCombEtaInvXSectionStatA,4,20. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvXSectionEtaStatLower)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPowInvXSectionEtaStatLower)<< endl;    

    
    canvasXSectionPi0->cd();
    TH2F * histo2DXSectionWithHFE;
    histo2DXSectionWithHFE          = new TH2F("histo2DXSectionWithHFE","histo2DXSectionWithHFE",11000,0.23,70.,1000,2e-2,10e11);
    SetStyleHistoTH2ForGraphs(histo2DXSectionWithHFE, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 0.9,1.45);
    histo2DXSectionWithHFE->GetXaxis()->SetMoreLogLabels();
    histo2DXSectionWithHFE->GetXaxis()->SetLabelOffset(-0.01);
    histo2DXSectionWithHFE->Draw("copy");


        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionSysA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphCombPi0InvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSysA, markerStyleComb, markerSizeComb, kGray+2 , kGray+2, widthLinesBoxes, kTRUE);
        graphCombEtaInvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphHFE, markerStyleComb, markerSizeComb, kCyan+2 , kCyan+2, widthLinesBoxes, kTRUE);
        graphHFE->Draw("pE2same");
    
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatA, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombPi0InvXSectionStatA->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStatA, markerStyleComb, markerSizeComb, kGray+2 , kGray+2);
        graphCombEtaInvXSectionStatA->Draw("p,same,z");
        DrawGammaSetMarkerTF1( fitTCMInvXSectionPi0, 7, 2, colorComb); 
        fitPowInvXSectionPi0Tot->SetRange(6,40);
        DrawGammaSetMarkerTF1( fitPowInvXSectionPi0Tot, 1, 2, colorComb-2); 
        fitPowInvXSectionPi0Tot->Draw("same");
        fitTCMInvXSectionPi0->Draw("same");
        
        fitTCMInvXSectionEta->SetRange(0.6,50.);
        DrawGammaSetMarkerTF1( fitTCMInvXSectionEta, 7, 2, kGray+2); 
        
        fitPowInvXSectionEtaTot->SetRange(6,20);
        DrawGammaSetMarkerTF1( fitPowInvXSectionEtaTot, 1, 2, kGray+1); 
        fitPowInvXSectionEtaTot->Draw("same");
        fitTCMInvXSectionEta->Draw("same");

        fitPowInvXSectionHFE->SetRange(5,50.);
        DrawGammaSetMarkerTF1( fitPowInvXSectionHFE, 7, 2, kCyan+2); 
        fitPowInvXSectionHFE->Draw("same");

        DrawGammaSetMarkerTGraphAsym(graphNLOCalcDirGam, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
        graphNLOCalcDirGam->Draw("3,same");
                        
        labelEnergyXSectionPi0->Draw();

        TLegend* legendXSectionWithHFE          = GetAndSetLegend2(0.76, 0.92-0.035*4, 0.98 , 0.92, 32); 
        legendXSectionWithHFE->AddEntry(graphCombPi0InvXSectionSysA,"#pi^{0}","fp");
        legendXSectionWithHFE->AddEntry(graphCombEtaInvXSectionSysA,"#eta","fp");
        legendXSectionWithHFE->AddEntry(graphHFE,"#frac{e^{+}+e^{-}}{2}","fp");
        legendXSectionWithHFE->AddEntry(graphNLOCalcDirGam,"#gamma_{dir} NLO","l");
        legendXSectionWithHFE->Draw();
   
    canvasXSectionPi0->SaveAs(Form("%s/InvXSection_Pi0_Eta_HFE.%s",outputDir.Data(),suffix.Data()));


    histo2DXSectionWithHFE->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionSysA, markerStyleComb, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
        graphCombPi0InvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatA, markerStyleComb, markerSizeComb, kBlack , kBlack);
        graphCombPi0InvXSectionStatA->Draw("p,same,z");

        
        DrawGammaSetMarkerTF1( fitInvXSectionPi0, 9, 2, kGray+2); 
        DrawGammaSetMarkerTF1( fitTCMInvXSectionPi0, 7, 2, kGray+1); 
        DrawGammaSetMarkerTF1( fitPowInvXSectionPi0Tot, 1, 2, kAzure+2); 
        DrawGammaSetMarkerTF1( fitPowInvXSectionPi0TotLower, 3, 2, kAzure+3); 
        DrawGammaSetMarkerTF1( fitTCMDecomposedLPi0, 4, 2, kRed-6); 
        DrawGammaSetMarkerTF1( fitTCMDecomposedHPi0, 5, 2, kRed-8); 
        fitPowInvXSectionPi0Tot->SetRange(6,50);
        fitPowInvXSectionPi0TotLower->SetRange(6,50);
        fitTCMInvXSectionPi0->SetRange(0.2,50);
        fitInvXSectionPi0->SetRange(0.2,50);
        fitTCMDecomposedLPi0->SetRange(0.2,10);
        fitTCMDecomposedHPi0->SetRange(0.2,50);
        fitPowInvXSectionPi0Tot->Draw("same");
        fitPowInvXSectionPi0TotLower->Draw("same");
        fitInvXSectionPi0->Draw("same");
        fitTCMInvXSectionPi0->Draw("same");
        fitTCMDecomposedLPi0->Draw("same");
        fitTCMDecomposedHPi0->Draw("same");
        
        DrawGammaSetMarkerTF1( fitPi0Tsallis1, 7, 2, 807); 
        DrawGammaSetMarkerTF1( fitPi0Tsallis2, 8, 2, kViolet+2); 
        DrawGammaSetMarkerTF1( fitPi0Tsallis3, 3, 2, kPink+2); 
        DrawGammaSetMarkerTF1( fitPi0Tsallis4, 6, 2, kGreen+2); 
        fitPi0Tsallis1->SetRange(0.2,50);
        fitPi0Tsallis2->SetRange(0.2,50);
        fitPi0Tsallis3->SetRange(0.2,50);
        fitPi0Tsallis4->SetRange(0.2,50);
        fitPi0Tsallis1->Draw("same");
        fitPi0Tsallis2->Draw("same");
        fitPi0Tsallis3->Draw("same");
        fitPi0Tsallis4->Draw("same");
        
        labelEnergyXSectionPi0->Draw();
        labelDetSysXSectionPi0->Draw();
        
        TLegend* legendXSectionPi0DiffFits          = GetAndSetLegend2(0.18, 0.13, 0.40 , 0.13+0.03*10, 32); 
        legendXSectionPi0DiffFits->AddEntry(fitInvXSectionPi0,"Levy-Tsallis","l");
        legendXSectionPi0DiffFits->AddEntry(fitPi0Tsallis1,"Levy-Tsallis, 0.8-40","l");
        legendXSectionPi0DiffFits->AddEntry(fitPi0Tsallis2,"Levy-Tsallis, 1.5-40","l");
        legendXSectionPi0DiffFits->AddEntry(fitPi0Tsallis3,"Levy-Tsallis, 2.5-40","l");
        legendXSectionPi0DiffFits->AddEntry(fitPi0Tsallis4,"Levy-Tsallis, 0.4-10","l");
        legendXSectionPi0DiffFits->AddEntry(fitPowInvXSectionPi0Tot,"pure powerlaw 6-40","l");
        legendXSectionPi0DiffFits->AddEntry(fitPowInvXSectionPi0TotLower,"pure powerlaw 4-40","l");
        legendXSectionPi0DiffFits->AddEntry(fitTCMInvXSectionPi0,"full TCM","l");
        legendXSectionPi0DiffFits->AddEntry(fitTCMDecomposedLPi0,"low p_{T} comp. TCM","l");
        legendXSectionPi0DiffFits->AddEntry(fitTCMDecomposedHPi0,"high p_{T} comp. TCM","l");
        legendXSectionPi0DiffFits->Draw();
        
    canvasXSectionPi0->SaveAs(Form("%s/InvXSection_Pi0_WithDiffFits.%s",outputDir.Data(),suffix.Data()));

    TF1* fitTsallisPi0PrevPubRefitted = (TF1*) fitPi0OldPub->Clone("fitTsallisPi0PrevPubRefitted");
    for (Int_t i = 0; i < fitTsallisPi0PrevPubRefitted->GetNpar(); i++){
      fitTsallisPi0PrevPubRefitted->SetParLimits(i, 0.1*fitTsallisPi0PrevPubRefitted->GetParameter(i), 1.9*fitTsallisPi0PrevPubRefitted->GetParameter(i));
    }
    fitTsallisPi0PrevPubRefitted->SetParameter(1,6.3);
    fitTsallisPi0PrevPubRefitted->SetParLimits(1,6.0,6.7);
    fitTsallisPi0PrevPubRefitted->SetRange(0.4,40);
    graphCombPi0InvXSectionTotA->Fit(fitTsallisPi0PrevPubRefitted,"QNRMEX0+","",0.4,40);
    fileFitsOutput <<  WriteParameterToFile(fitTsallisPi0PrevPubRefitted)<< endl;    
    
    histo2DXSectionWithHFE->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionSysA, markerStyleComb, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
        graphCombPi0InvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatA, markerStyleComb, markerSizeComb, kBlack , kBlack);
        graphCombPi0InvXSectionStatA->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphPi0PubOldSys, markerStyleComb+4, markerSizeComb, kGray+1 , kGray+1, widthLinesBoxes, kTRUE);
        graphPi0PubOldSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPi0PubOldStat, markerStyleComb+4, markerSizeComb, kGray+1 , kGray+1);
        graphPi0PubOldStat->Draw("p,same,z");
        
        fitPowInvXSectionPi0Tot->Draw("same");
        
        DrawGammaSetMarkerTF1( fitTCMInvXSectionPi0, 7, 2, kGray+1); 
        fitTCMInvXSectionPi0->Draw("same");
        
        DrawGammaSetMarkerTF1( fitPi0OldPub, 6, 2, kGreen+2); 
        fitPi0OldPub->SetRange(0.2,50);
        fitPi0OldPub->Draw("same");
        
        DrawGammaSetMarkerTF1( fitTsallisPi0PrevPubRefitted, 5, 2, kRed+2); 
        fitTsallisPi0PrevPubRefitted->SetRange(0.2,50);
        fitTsallisPi0PrevPubRefitted->Draw("same");
        
        labelEnergyXSectionPi0->Draw();
        labelDetSysXSectionPi0->Draw();
        
        TLegend* legendXSectionPi0WithOld          = GetAndSetLegend2(0.18, 0.13, 0.40 , 0.13+0.03*6, 32); 
        legendXSectionPi0WithOld->AddEntry(graphCombPi0InvXSectionSysA, "updated","pf");
        legendXSectionPi0WithOld->AddEntry(fitTCMInvXSectionPi0,"full TCM","l");
        legendXSectionPi0WithOld->AddEntry(fitTsallisPi0PrevPubRefitted,"Levy-Tsallis, restr. param.","l");
        legendXSectionPi0WithOld->AddEntry(fitPowInvXSectionPi0Tot,"pure powerlaw 6-40","l");
        legendXSectionPi0WithOld->AddEntry(graphPi0PubOldSys, "old pub","pf");
        legendXSectionPi0WithOld->AddEntry(fitPi0OldPub,"Levy-Tsallis, old pub","l");
        legendXSectionPi0WithOld->Draw();
        
    canvasXSectionPi0->SaveAs(Form("%s/InvXSection_Pi0_WithOld.%s",outputDir.Data(),suffix.Data()));
    
    
    
    TCanvas* canvasRatioPP2                 = new TCanvas("canvasRatioPP2","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioPP2,  0.12, 0.01, 0.01, 0.11);
    canvasRatioPP2->cd();
    canvasRatioPP2->SetLogx();
        textsizeLabelsPP                    = 0;
        textsizeFacPP                       = 0;
        if (canvasRatioPP2->XtoPixel(canvasRatioPP2->GetX2()) <canvasRatioPP2->YtoPixel(canvasRatioPP2->GetY1()) ){
            textsizeLabelsPP                = (Double_t)textSizeLabelsPixel/canvasRatioPP2->XtoPixel(canvasRatioPP2->GetX2()) ;
            textsizeFacPP                   = (Double_t)1./canvasRatioPP2->XtoPixel(canvasRatioPP2->GetX2()) ;
        } else {
        textsizeLabelsPP                    = (Double_t)textSizeLabelsPixel/canvasRatioPP2->YtoPixel(canvasRatioPP2->GetY1());
        textsizeFacPP                       = (Double_t)1./canvasRatioPP2->YtoPixel(canvasRatioPP2->GetY1());
        }
        cout << textsizeLabelsPP << endl;

        TH2F * ratio2DFits       = new TH2F("ratio2DFits","ratio2DFits",1000,0.23,70.,1000,0.2,2.1);
        SetStyleHistoTH2ForGraphs(ratio2DFits, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
                                  0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
        ratio2DFits->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio2DFits->GetYaxis()->SetNdivisions(505);
        ratio2DFits->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DFits->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DFits->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DFits->GetXaxis()->SetLabelFont(42);
        ratio2DFits->GetYaxis()->SetLabelFont(42);

        ratio2DFits->DrawCopy(); 

             
        TGraphAsymmErrors* graphRatioPi0TCMTot      = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotA->Clone();
        graphRatioPi0TCMTot                         = CalculateGraphErrRatioToFit(graphRatioPi0TCMTot, fitTCMInvXSectionPi0); 
        TGraphAsymmErrors* graphRatioPi0OldTCMTot   = (TGraphAsymmErrors*)graphPi0PubOldTot->Clone();
        graphRatioPi0OldTCMTot                      = CalculateGraphErrRatioToFit(graphRatioPi0OldTCMTot, fitTCMInvXSectionPi0); 

        TGraphAsymmErrors* graphRatioPi0TsallisTot  = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotA->Clone();
        graphRatioPi0TsallisTot                     = CalculateGraphErrRatioToFit(graphRatioPi0TsallisTot, fitInvXSectionPi0); 
        TGraphAsymmErrors* graphRatioPi0TsallisOld  = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotA->Clone();
        graphRatioPi0TsallisOld                     = CalculateGraphErrRatioToFit(graphRatioPi0TsallisOld, fitTsallisPi0PrevPubRefitted); 
        TGraphAsymmErrors* graphRatioPi0Tsallis1Tot = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotA->Clone();
        graphRatioPi0Tsallis1Tot                    = CalculateGraphErrRatioToFit(graphRatioPi0Tsallis1Tot, fitPi0Tsallis1); 
        TGraphAsymmErrors* graphRatioPi0Tsallis2Tot = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotA->Clone();
        graphRatioPi0Tsallis2Tot                    = CalculateGraphErrRatioToFit(graphRatioPi0Tsallis2Tot, fitPi0Tsallis2); 
        TGraphAsymmErrors* graphRatioPi0Tsallis3Tot = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotA->Clone();
        graphRatioPi0Tsallis3Tot                    = CalculateGraphErrRatioToFit(graphRatioPi0Tsallis3Tot, fitPi0Tsallis3); 
        TGraphAsymmErrors* graphRatioPi0Tsallis4Tot = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotA->Clone();
        graphRatioPi0Tsallis4Tot                    = CalculateGraphErrRatioToFit(graphRatioPi0Tsallis4Tot, fitPi0Tsallis4); 
        TGraphAsymmErrors* graphRatioPi0PowTot      = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotA->Clone();
        graphRatioPi0PowTot                         = CalculateGraphErrRatioToFit(graphRatioPi0PowTot, fitPowInvXSectionPi0Tot); 
        TGraphAsymmErrors* graphRatioPi0PowTotLow   = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotA->Clone();
        graphRatioPi0PowTotLow                      = CalculateGraphErrRatioToFit(graphRatioPi0PowTotLow, fitPowInvXSectionPi0TotLower); 
        
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0TsallisTot, 25, markerSizeComb, kGray+1, kGray+1, widthLinesBoxes, kFALSE);
        graphRatioPi0TsallisTot->Draw("p,same");
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0Tsallis1Tot, 25, markerSizeComb, 807, 807, widthLinesBoxes, kFALSE);
//         graphRatioPi0Tsallis1Tot->Draw("p,same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0Tsallis2Tot, 25, markerSizeComb, kViolet+2, kViolet+2, widthLinesBoxes, kFALSE);
        graphRatioPi0Tsallis2Tot->Draw("p,same");
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0Tsallis3Tot, 25, markerSizeComb, kPink+2, kPink+2, widthLinesBoxes, kFALSE);
//         graphRatioPi0Tsallis3Tot->Draw("p,same");
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0Tsallis4Tot, 25, markerSizeComb, kGreen+2, kGreen+2, widthLinesBoxes, kFALSE);
//         graphRatioPi0Tsallis4Tot->Draw("p,same");
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0PowTot, 29, markerSizeComb*2, kAzure+2, kAzure+2, widthLinesBoxes, kFALSE);
//         graphRatioPi0PowTot->Draw("p,same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PowTotLow, 33, markerSizeComb*2, kAzure+3, kAzure+3, widthLinesBoxes, kFALSE);
        graphRatioPi0PowTotLow->Draw("p,same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0TsallisOld, 24, markerSizeComb, kRed+2, kRed+2, widthLinesBoxes, kFALSE);
        graphRatioPi0TsallisOld->Draw("p,same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0TCMTot, markerStyleComb, markerSizeComb, kGray+2, kGray+2, widthLinesBoxes, kFALSE);
        graphRatioPi0TCMTot->Draw("p,same");

        DrawGammaLines(0.23, 70. , 1.1, 1.1,1, kGray, 3);
        DrawGammaLines(0.23, 70. , 1, 1,1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.9, 0.9,1, kGray, 3);
        
        TLegend* legendRatioPi0Fits= GetAndSetLegend2(0.15,0.95-5*0.9*textsizeLabelsPP,0.4,0.95, 0.85* textSizeLabelsPixel, 1, "", 43, 0.2);
        legendRatioPi0Fits->AddEntry(graphRatioPi0TCMTot,"full TCM","p");
        legendRatioPi0Fits->AddEntry(graphRatioPi0TsallisTot,"Levy-Tsallis, free","p");
//         legendRatioPi0Fits->AddEntry(graphRatioPi0Tsallis1Tot,"Levy-Tsallis, 0.8-40","p");
        legendRatioPi0Fits->AddEntry(graphRatioPi0Tsallis2Tot,"Levy-Tsallis, 1.5-40","pl");
//         legendRatioPi0Fits->AddEntry(graphRatioPi0Tsallis3Tot,"Levy-Tsallis, 2.5-40","p");
//         legendRatioPi0Fits->AddEntry(graphRatioPi0Tsallis4Tot,"Levy-Tsallis, 0.4-10","p");
        legendRatioPi0Fits->AddEntry(graphRatioPi0TsallisOld,"Levy-Tsallis, restr. param.","p");
//         legendRatioPi0Fits->AddEntry(graphRatioPi0PowTot,"pure powerlaw 6-40","p");
        legendRatioPi0Fits->AddEntry(graphRatioPi0PowTotLow,"pure powerlaw 4-40","p");
        legendRatioPi0Fits->Draw();
        
        // upper right corner
//         legendRatioPi0OldTheo2->Draw();
        labelRatioTheoryPP2->Draw();
        TLatex *labelRatioFitsPi0   = new TLatex(0.73,0.92-0.9*textsizeLabelsPP,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioFitsPi0, 0.85*textsizeLabelsPP,4);
        labelRatioFitsPi0->Draw();

        ratio2DFits->Draw("axis,same");
        
    canvasRatioPP2->Update();
    canvasRatioPP2->Print(Form("%s/Pi0_RatioDiffFits.%s",outputDir.Data(),suffix.Data()));

    
    
    histo2DXSectionWithHFE->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSysA, markerStyleComb, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
        graphCombEtaInvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStatA, markerStyleComb, markerSizeComb, kBlack , kBlack);
        graphCombEtaInvXSectionStatA->Draw("p,same,z");

        DrawGammaSetMarkerTF1( fitInvXSectionEta, 9, 2, kGray+2); 
        DrawGammaSetMarkerTF1( fitTCMInvXSectionEta, 7, 2, kGray+1); 
        DrawGammaSetMarkerTF1( fitPowInvXSectionEtaTot, 1, 2, kAzure+2); 
        DrawGammaSetMarkerTF1( fitTCMDecomposedLEta, 4, 2, kRed-6); 
        DrawGammaSetMarkerTF1( fitTCMDecomposedHEta, 5, 2, kRed-8); 
        fitPowInvXSectionEtaTot->SetRange(6,50);
        fitTCMInvXSectionEta->SetRange(0.2,50);
        fitInvXSectionEta->SetRange(0.2,50);
        fitTCMDecomposedLEta->SetRange(0.2,10);
        fitTCMDecomposedHEta->SetRange(0.2,50);
        fitPowInvXSectionEtaTot->Draw("same");
        fitInvXSectionEta->Draw("same");
        fitTCMInvXSectionEta->Draw("same");
        fitTCMDecomposedLEta->Draw("same");
        fitTCMDecomposedHEta->Draw("same");
        
        labelEnergyXSectionPi0->Draw();
        labelDetSysXSectionEta->Draw();
        
        TLegend* legendXSectionEtaDiffFits          = GetAndSetLegend2(0.18, 0.13, 0.40 , 0.13+0.03*5, 32); 
        legendXSectionEtaDiffFits->AddEntry(fitInvXSectionEta,"Levy-Tsallis","l");
        legendXSectionEtaDiffFits->AddEntry(fitPowInvXSectionEtaTot,"pure powerlaw","l");
        legendXSectionEtaDiffFits->AddEntry(fitTCMInvXSectionEta,"full TCM","l");
        legendXSectionEtaDiffFits->AddEntry(fitTCMDecomposedLEta,"low p_{T} comp. TCM","l");
        legendXSectionEtaDiffFits->AddEntry(fitTCMDecomposedHEta,"high p_{T} comp. TCM","l");
        legendXSectionEtaDiffFits->Draw();
        
    canvasXSectionPi0->SaveAs(Form("%s/InvXSection_Eta_WithDiffFits.%s",outputDir.Data(),suffix.Data()));
    
    if (plotInvMassBins){
        textSizeLabelsPixel                 = 100*3/5;
        TCanvas* canvasInvMassSamplePlot    = new TCanvas("canvasInvMassSamplePlot","",0,0,1500,1500);  // gives the page size
        DrawGammaCanvasSettings( canvasInvMassSamplePlot,  0.09, 0.01, 0.012, 0.08);

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
        histo2DPi0InvMassDummy             = new TH2F("histo2DPi0InvMassDummy","histo2DPi0InvMassDummy",11000,0.05,0.255,21000,-1000,20000);
        SetStyleHistoTH2ForGraphs(histo2DPi0InvMassDummy, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","Counts",0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,
                                0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,0.88, 0.115/(textsizeFacInvMass*marginInvMass));

        TH2F * histo2DEtaInvMassDummy;
        histo2DEtaInvMassDummy             = new TH2F("histo2DEtaInvMassDummy","histo2DEtaInvMassDummy",11000,0.35,0.695,21000,-1000,20000);
        SetStyleHistoTH2ForGraphs(histo2DEtaInvMassDummy, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","Counts",0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,
                                0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,0.88, 0.115/(textsizeFacInvMass*marginInvMass));
        
        for (Int_t i =0 ; i < 6; i++){
            if (haveAllPi0InvMassPCMEMCAL[i]){
                canvasInvMassSamplePlot->cd();
                histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSubPCMEMCAL[i]->GetMinimum(),1.15*histoPi0InvMassSigPlusBGPCMEMCAL[i]->GetMaximum());
                histo2DPi0InvMassDummy->DrawCopy();

                TLatex *labelInvMassPtRangePCMEMCAL = new TLatex(0.955,0.93,Form("#pi^{0}: %s", histoPi0InvMassSigPCMEMCAL[i]->GetTitle()));
            
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
                
                TLatex *labelInvMassPerf      = new TLatex(0.135,0.93,"ALICE performance");
                SetStyleTLatex( labelInvMassPerf, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
                labelInvMassPerf->Draw();

                TLatex *labelInvMassEnergy      = new TLatex(0.135,0.93-0.8*textsizeLabelsPP,collisionSystem2760GeV.Data());
                SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
                labelInvMassEnergy->Draw();
                
                TLatex *labelInvMassTrigger      = new TLatex(0.135,0.93-2*0.8*textsizeLabelsPP,Form("%s triggered",nameTrigger[i].Data()));
                SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
                labelInvMassTrigger->Draw();

                TLatex *labelInvMassRecoPCMEMC  = new TLatex(0.135,0.93-3*0.8*textsizeLabelsPP,nameMeasGlobalLabel[4]);
                SetStyleTLatex( labelInvMassRecoPCMEMC, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
                labelInvMassRecoPCMEMC->Draw();
                
                SetStyleTLatex( labelInvMassPtRangePCMEMCAL, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
                labelInvMassPtRangePCMEMCAL->Draw();

                TLegend* legendInvMassPCMEMCAL  = GetAndSetLegend2(0.67, 0.91-5*0.75*textsizeLabelsPP, 0.9, 0.91, 0.85*textSizeLabelsPixel, 1, "", 43, 0.25);
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

            if (haveAllEtaInvMassPCMEMCAL[i]){
                canvasInvMassSamplePlot->cd();
                histo2DEtaInvMassDummy->GetYaxis()->SetRangeUser(histoEtaInvMassSigRemBGSubPCMEMCAL[i]->GetMinimum()*0.95,1.15*histoEtaInvMassSigPlusBGPCMEMCAL[i]->GetMaximum());
                histo2DEtaInvMassDummy->DrawCopy();

                TLatex *labelInvMassPtRangePCMEMCAL = new TLatex(0.955,0.93,Form("#eta: %s", histoEtaInvMassSigPCMEMCAL[i]->GetTitle()));
            
                DrawGammaSetMarker(histoEtaInvMassSigPlusBGPCMEMCAL[i], markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
                histoEtaInvMassSigPlusBGPCMEMCAL[i]->SetLineWidth(1);
                histoEtaInvMassSigPlusBGPCMEMCAL[i]->Draw("hist,e,same");
                DrawGammaSetMarker(histoEtaInvMassBGTotPCMEMCAL[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
                histoEtaInvMassBGTotPCMEMCAL[i]->Draw("same");


                DrawGammaSetMarker(histoEtaInvMassSigRemBGSubPCMEMCAL[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
                histoEtaInvMassSigRemBGSubPCMEMCAL[i]->Draw("same");
                fitEtaInvMassSigPCMEMCAL[i]->SetRange(0.35,0.695);
                fitEtaInvMassSigPCMEMCAL[i]->SetLineColor(fitColorInvMassSG);
                fitEtaInvMassSigPCMEMCAL[i]->Draw("same");

                TLatex *labelInvMassPerf      = new TLatex(0.135,0.93,"ALICE performance");
                SetStyleTLatex( labelInvMassPerf, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
                labelInvMassPerf->Draw();
                
                TLatex *labelInvMassEnergy      = new TLatex(0.135,0.93-0.8*textsizeLabelsPP,collisionSystem2760GeV.Data());
                SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
                labelInvMassEnergy->Draw();
                
                TLatex *labelInvMassTrigger      = new TLatex(0.135,0.93-2*0.8*textsizeLabelsPP,Form("%s triggered",nameTrigger[i].Data()));
                SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
                labelInvMassTrigger->Draw();

                TLatex *labelInvMassRecoPCMEMC  = new TLatex(0.135,0.93-3*0.8*textsizeLabelsPP,nameMeasGlobalLabel[4]);
                SetStyleTLatex( labelInvMassRecoPCMEMC, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
                labelInvMassRecoPCMEMC->Draw();
                
                SetStyleTLatex( labelInvMassPtRangePCMEMCAL, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
                labelInvMassPtRangePCMEMCAL->Draw();

                TLegend* legendInvMassPCMEMCAL  = GetAndSetLegend2(0.67, 0.91-5*0.75*textsizeLabelsPP, 0.9, 0.91, 0.85*textSizeLabelsPixel, 1, "", 43, 0.25);
                legendInvMassPCMEMCAL->AddEntry(histoEtaInvMassSigPlusBGPCMEMCAL[i],"Raw real events","l");
                legendInvMassPCMEMCAL->AddEntry(histoEtaInvMassBGTotPCMEMCAL[i],"Mixed event +","p");
                legendInvMassPCMEMCAL->AddEntry((TObject*)0,"corr. BG","");
                legendInvMassPCMEMCAL->AddEntry(histoEtaInvMassSigRemBGSubPCMEMCAL[i],"BG subtracted","p");
                legendInvMassPCMEMCAL->AddEntry(fitEtaInvMassSigPCMEMCAL[i], "Fit","l");
                legendInvMassPCMEMCAL->Draw();
                canvasInvMassSamplePlot->SaveAs(Form("%s/Eta_InvMassBinPCMEMC_%s.%s",outputDir.Data(), nameTrigger[i].Data(), suffix.Data()));
            } else {
                cout << "missing partial input for invariant mass bin for PCM-EMC for trigger: " << nameTrigger[i].Data() << endl;
            }
            
            
            if (haveAllPi0InvMassEMCAL[i]){
                canvasInvMassSamplePlot->cd();
                histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSubEMCAL[i]->GetMinimum(),1.15*histoPi0InvMassSigPlusBGEMCAL[i]->GetMaximum());
                histo2DPi0InvMassDummy->DrawCopy();

                TLatex *labelInvMassPtRangeEMCAL = new TLatex(0.955,0.93,Form("#pi^{0}: %s", histoPi0InvMassSigEMCAL[i]->GetTitle()));
            
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
                TLatex *labelInvMassPerf      = new TLatex(0.135,0.93,"ALICE performance");
                SetStyleTLatex( labelInvMassPerf, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
                labelInvMassPerf->Draw();

                TLatex *labelInvMassEnergy      = new TLatex(0.135,0.93-0.8*textsizeLabelsPP,collisionSystem2760GeV.Data());
                SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
                labelInvMassEnergy->Draw();
                
                TLatex *labelInvMassTrigger      = new TLatex(0.135,0.93-0.8*2*textsizeLabelsPP,Form("%s triggered",nameTrigger[i].Data()));
                SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
                labelInvMassTrigger->Draw();

                TLatex *labelInvMassRecoEMC  = new TLatex(0.135,0.93-3*0.8*textsizeLabelsPP,nameMeasGlobalLabel[2]);
                SetStyleTLatex( labelInvMassRecoEMC, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
                labelInvMassRecoEMC->Draw();
                
                SetStyleTLatex( labelInvMassPtRangeEMCAL, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
                labelInvMassPtRangeEMCAL->Draw();

                TLegend* legendInvMassEMCAL  = GetAndSetLegend2(0.67, 0.91-5*0.75*textsizeLabelsPP, 0.9, 0.91, 0.85*textSizeLabelsPixel,1, "", 43, 0.25);
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
            
            if (haveAllEtaInvMassEMCAL[i]){
                canvasInvMassSamplePlot->cd();
                histo2DEtaInvMassDummy->GetYaxis()->SetRangeUser(histoEtaInvMassSigRemBGSubEMCAL[i]->GetMinimum(),1.1*histoEtaInvMassSigPlusBGEMCAL[i]->GetMaximum());
                histo2DEtaInvMassDummy->DrawCopy();

                TLatex *labelInvMassPtRangeEMCAL = new TLatex(0.955,0.93,Form("#eta: %s", histoEtaInvMassSigEMCAL[i]->GetTitle()));
            
                DrawGammaSetMarker(histoEtaInvMassSigPlusBGEMCAL[i], markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
                histoEtaInvMassSigPlusBGEMCAL[i]->SetLineWidth(1);
                histoEtaInvMassSigPlusBGEMCAL[i]->Draw("hist,e,same");
                DrawGammaSetMarker(histoEtaInvMassBGTotEMCAL[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
                histoEtaInvMassBGTotEMCAL[i]->Draw("same");

                DrawGammaSetMarker(histoEtaInvMassSigRemBGSubEMCAL[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
                histoEtaInvMassSigRemBGSubEMCAL[i]->Draw("same");
                fitEtaInvMassSigEMCAL[i]->SetRange(0.35,0.695);
                fitEtaInvMassSigEMCAL[i]->SetLineColor(fitColorInvMassSG);
                fitEtaInvMassSigEMCAL[i]->Draw("same");
                
                TLatex *labelInvMassPerf      = new TLatex(0.135,0.93,"ALICE performance");
                SetStyleTLatex( labelInvMassPerf, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
                labelInvMassPerf->Draw();
                
                TLatex *labelInvMassEnergy      = new TLatex(0.135,0.93-0.8*textsizeLabelsPP,collisionSystem2760GeV.Data());
                SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
                labelInvMassEnergy->Draw();
                
                TLatex *labelInvMassTrigger      = new TLatex(0.135,0.93-2*0.8*textsizeLabelsPP,Form("%s triggered",nameTrigger[i].Data()));
                SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
                labelInvMassTrigger->Draw();

                TLatex *labelInvMassRecoEMC  = new TLatex(0.135,0.93-3*0.8*textsizeLabelsPP,nameMeasGlobalLabel[2]);
                SetStyleTLatex( labelInvMassRecoEMC, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
                labelInvMassRecoEMC->Draw();
                
                SetStyleTLatex( labelInvMassPtRangeEMCAL, 0.85*textSizeLabelsPixel,4, 1, 43, kTRUE, 31);
                labelInvMassPtRangeEMCAL->Draw();

                TLegend* legendInvMassEMCAL  = GetAndSetLegend2(0.67, 0.91-5*0.75*textsizeLabelsPP, 0.9, 0.91, 0.85*textSizeLabelsPixel, 1, "", 43, 0.25);
                legendInvMassEMCAL->AddEntry(histoEtaInvMassSigPlusBGEMCAL[i],"Raw real events","l");
                legendInvMassEMCAL->AddEntry(histoEtaInvMassBGTotEMCAL[i],"Mixed event +","p");
                legendInvMassEMCAL->AddEntry((TObject*)0,"corr. BG","");
                legendInvMassEMCAL->AddEntry(histoEtaInvMassSigRemBGSubEMCAL[i],"BG subtracted","p");
                legendInvMassEMCAL->AddEntry(fitEtaInvMassSigEMCAL[i], "Fit","l");
                legendInvMassEMCAL->Draw();
                canvasInvMassSamplePlot->SaveAs(Form("%s/Eta_InvMassBinEMC_%s.%s",outputDir.Data(), nameTrigger[i].Data(), suffix.Data()));
            } else {
                cout << "missing partial input for invariant mass bin for PCM-EMC for trigger: " << nameTrigger[i].Data() << endl;
            }            
        }        
    }    
    
    // close fit log file
    fileFitsOutput.close();
    
 // **********************************************************************************************************************
 // ************************* Saving of final results ********************************************************************
 // **********************************************************************************************************************
 
    TString nameOutputCommonFile    = "";
    if (flagMerged == 0) 
        nameOutputCommonFile        = Form("CombinedResultsPaperPP2760GeV_%s_Haitao.root", dateForOutput.Data());
    else if (flagMerged == 1) 
        nameOutputCommonFile        = Form("CombinedResultsPaperPP2760GeV_%s_FrediV1Clusterizer.root", dateForOutput.Data());
    else if (flagMerged == 2)
        nameOutputCommonFile        = Form("CombinedResultsPaperPP2760GeV_%s_FrediV2Clusterizer.root", dateForOutput.Data());
    else if (flagMerged == 3) 
        nameOutputCommonFile        = Form("CombinedResultsPaperPP2760GeV_%s_FrediV1ClusterizerNLM1.root", dateForOutput.Data());
    else if (flagMerged == 4) 
        nameOutputCommonFile        = Form("CombinedResultsPaperPP2760GeV_%s_FrediV1ClusterizerNLM2.root", dateForOutput.Data());
    
    TFile fCombResults(nameOutputCommonFile.Data(), "RECREATE");

    fCombResults.mkdir("Pi02.76TeV");
    TDirectoryFile* directoryPi0 = (TDirectoryFile*)fCombResults.Get("Pi02.76TeV"); 
    fCombResults.cd("Pi02.76TeV");
        // PCM component
        graphPCMPi0InvXSectionStat->Write("graphInvCrossSectionPi0PCM2760GeVStatErr");
        graphPCMPi0InvXSectionSys->Write("graphInvCrossSectionPi0PCM2760GeVSysErr");
        // PHOS component
        graphPHOSPi0InvXSectionStat->Write("graphInvCrossSectionPi0PHOS2760GeVStatErr");
        graphPHOSPi0InvXSectionSys->Write("graphInvCrossSectionPi0PHOS2760GeVSysErr");
        // EMCAL component
        graphEMCALPi0InvXSectionStat->Write("graphInvCrossSectionPi0EMCAL2760GeVStatErr");
        graphEMCALPi0InvXSectionSys->Write("graphInvCrossSectionPi0EMCAL2760GeVSysErr");
        // PCM-EMC component
        graphPCMEMCALPi0InvXSectionStat->Write("graphInvCrossSectionPi0PCMEMCAL2760GeVStatErr");
        graphPCMEMCALPi0InvXSectionSys->Write("graphInvCrossSectionPi0PCMEMCAL2760GeVSysErr");
        // EMCal merged component
        graphEMCALMergedPi0InvXSectionStat->Write("graphInvCrossSectionPi0EMCALMerged2760GeVStatErr");
        graphEMCALMergedPi0InvXSectionSys->Write("graphInvCrossSectionPi0EMCALMerged2760GeVSysErr");
        // Final spectrum correlations Method A
        graphCombPi0InvXSectionTotA->Write("graphInvCrossSectionPi0Comb2760GeVATotErr");
        graphCombPi0InvXSectionStatA->Write("graphInvCrossSectionPi0Comb2760GeVAStatErr");
        graphCombPi0InvXSectionSysA->Write("graphInvCrossSectionPi0Comb2760GeVASysErr");  
        // Final inv yield INEL
        graphCombPi0InvYieldTotA->Write("graphInvYieldINELPi0Comb2760GeVATotErr");
        graphCombPi0InvYieldStatA->Write("graphInvYieldINELPi0Comb2760GeVAStatErr");
        graphCombPi0InvYieldSysA->Write("graphInvYieldINELPi0Comb2760GeVASysErr");  

         // fits for eta
        fitInvXSectionPi0->Write("TsallisFitPi0");
        fitTCMInvXSectionPi0->Write("TwoComponentModelFitPi0");
       
        TF1* fitTCMInvYieldSectionPi0 = (TF1*)fitTCMInvXSectionPi0->Clone("fitTCMInvYieldSectionPi0");
        fitTCMInvYieldSectionPi0->SetParameter(0,fitTCMInvXSectionPi0->GetParameter(0)/xSection2760GeVINEL);
        fitTCMInvYieldSectionPi0->SetParError(0,fitTCMInvXSectionPi0->GetParError(0)/xSection2760GeVINEL);
        fitTCMInvYieldSectionPi0->SetParameter(2,fitTCMInvXSectionPi0->GetParameter(2)/xSection2760GeVINEL);
        fitTCMInvYieldSectionPi0->SetParError(2,fitTCMInvXSectionPi0->GetParError(2)/xSection2760GeVINEL);
        fitTCMInvYieldSectionPi0->Write("TwoComponentModelFitYieldPi0");
               
        if (bWCorrection.Contains("Y")){
            if(graphPCMPi0InvXSectionStat_yShifted)graphPCMPi0InvXSectionStat_yShifted->Write("graphInvCrossSectionPi0PCM2760GeVStatErr_yShifted");
            if(graphPCMPi0InvXSectionSys_yShifted)graphPCMPi0InvXSectionSys_yShifted->Write("graphInvCrossSectionPi0PCM2760GeVSysErr_yShifted");
        // PHOS component
            if(graphPHOSPi0InvXSectionStat_yShifted) graphPHOSPi0InvXSectionStat_yShifted->Write("graphInvCrossSectionPi0PHOS2760GeVStatErr_yShifted");
            if(graphPHOSPi0InvXSectionSys_yShifted) graphPHOSPi0InvXSectionSys_yShifted->Write("graphInvCrossSectionPi0PHOS2760GeVSysErr_yShifted");
            // EMCAL component
            if(graphEMCALPi0InvXSectionStat_yShifted)graphEMCALPi0InvXSectionStat_yShifted->Write("graphInvCrossSectionPi0EMCAL2760GeVStatErr_yShifted");
            if(graphEMCALPi0InvXSectionSys_yShifted)graphEMCALPi0InvXSectionSys_yShifted->Write("graphInvCrossSectionPi0EMCAL2760GeVSysErr_yShifted");
            // PCM-EMC component
            if(graphPCMEMCALPi0InvXSectionStat_yShifted)graphPCMEMCALPi0InvXSectionStat_yShifted->Write("graphInvCrossSectionPi0PCMEMCAL2760GeVStatErr_yShifted");
            if(graphPCMEMCALPi0InvXSectionSys_yShifted)graphPCMEMCALPi0InvXSectionSys_yShifted->Write("graphInvCrossSectionPi0PCMEMCAL2760GeVSysErr_yShifted");
            // EMCal merged component
            if (graphEMCALMergedPi0InvXSectionStat_yShifted) graphEMCALMergedPi0InvXSectionStat_yShifted->Write("graphInvCrossSectionPi0EMCALMerged2760GeVStatErr_yShifted");
            if (graphEMCALMergedPi0InvXSectionSys_yShifted) graphEMCALMergedPi0InvXSectionSys_yShifted->Write("graphInvCrossSectionPi0EMCALMerged2760GeVSysErr_yShifted");
            // Final spectrum correlations Method A
            if(graphCombPi0InvXSectionTotA_yShifted)graphCombPi0InvXSectionTotA_yShifted->Write("graphInvCrossSectionPi0Comb2760GeVATotErr_yShifted");
            if(graphCombPi0InvXSectionStatA_yShifted)graphCombPi0InvXSectionStatA_yShifted->Write("graphInvCrossSectionPi0Comb2760GeVAStatErr_yShifted");
            if(graphCombPi0InvXSectionSysA_yShifted)graphCombPi0InvXSectionSysA_yShifted->Write("graphInvCrossSectionPi0Comb2760GeVASysErr_yShifted");  

        }    
        directoryPi0->mkdir("Supporting");
        directoryPi0->cd("Supporting");
            // Writing full correction factors
            graphPCMPi0AccTimesEff->Write("Pi0CorrectionFactorPCM");
            graphPHOSPi0AccTimesEff->Write("Pi0CorrectionFactorPHOS");
            graphEMCALPi0AccTimesEff->Write("Pi0CorrectionFactorEMCAL");
            graphPCMEMCALPi0AccTimesEff->Write("Pi0CorrectionFactorPCMEMCAL");
            graphEMCALMergedPi0AccTimesEffDivPur->Write("Pi0CorrectionFactorEMCALMerged");
        
            graphPCMPi0Mass->Write("Pi0MassDataPCM");
            graphPCMPi0MassMC->Write("Pi0MassMCPCM");
            histoPHOSPi0Mass->Write("Pi0MassDataPHOS");
            histoPHOSPi0TrueMass->Write("Pi0MassMCPHOS");
            graphEMCALPi0Mass->Write("Pi0MassDataEMCAL");
            graphEMCALPi0MassMC->Write("Pi0MassMCEMCAL");
            graphPCMEMCALPi0Mass->Write("Pi0MassDataPCMEMCAL");
            graphPCMEMCALPi0MassMC->Write("Pi0MassMCPCMEMCAL");

            graphPCMPi0FWHM->Write("Pi0WidthDataPCM");
            graphPCMPi0FWHMMC->Write("Pi0WidthMCPCM");
            histoPHOSPi0FWHMMeV->Write("Pi0WidthDataPHOS");
            histoPHOSPi0TrueFWHMMeV->Write("Pi0WidthMCPHOS");
            graphEMCALPi0FWHM->Write("Pi0WidthDataEMCAL");
            graphEMCALPi0FWHMMC->Write("Pi0WidthMCEMCAL");
            graphPCMEMCALPi0FWHM->Write("Pi0WidthDataPCMEMCAL");
            graphPCMEMCALPi0FWHMMC->Write("Pi0WidthMCPCMEMCAL");
            
    fCombResults.mkdir("Eta2.76TeV");
    TDirectoryFile* directoryEta = (TDirectoryFile*)fCombResults.Get("Eta2.76TeV"); 
    fCombResults.cd("Eta2.76TeV");
        // PCM component
        graphPCMEtaInvXSectionStat->Write("graphInvCrossSectionEtaPCM2760GeVStatErr");
        graphPCMEtaInvXSectionSys->Write("graphInvCrossSectionEtaPCM2760GeVSysErr");
        // EMCAL component
        graphEMCALEtaInvXSectionStat->Write("graphInvCrossSectionEtaEMCAL2760GeVStatErr");
        graphEMCALEtaInvXSectionSys->Write("graphInvCrossSectionEtaEMCAL2760GeVSysErr");
        // PCM-EMC component
        graphPCMEMCALEtaInvXSectionStat->Write("graphInvCrossSectionEtaPCMEMCAL2760GeVStatErr");
        graphPCMEMCALEtaInvXSectionSys->Write("graphInvCrossSectionEtaPCMEMCAL2760GeVSysErr");
        // Final spectrum correlations Method A
        graphCombEtaInvXSectionTotA->Write("graphInvCrossSectionEtaComb2760GeVATotErr");
        graphCombEtaInvXSectionStatA->Write("graphInvCrossSectionEtaComb2760GeVAStatErr");
        graphCombEtaInvXSectionSysA->Write("graphInvCrossSectionEtaComb2760GeVASysErr");  

        // Final yield for INEL
        graphCombEtaInvYieldTotA->Write("graphInvYieldINELEtaComb2760GeVATotErr");
        graphCombEtaInvYieldStatA->Write("graphInvYieldINELEtaComb2760GeVAStatErr");
        graphCombEtaInvYieldSysA->Write("graphInvYieldINELEtaComb2760GeVASysErr");  
        
        // fits for eta
        fitInvXSectionEta->Write("TsallisFitEta");
        fitTCMInvXSectionEta->Write("TwoComponentModelFitEta");
        TF1* fitTCMInvYieldSectionEta = (TF1*)fitTCMInvXSectionEta->Clone("fitTCMInvYieldSectionEta");
        fitTCMInvYieldSectionEta->SetParameter(0,fitTCMInvXSectionEta->GetParameter(0)/xSection2760GeVINEL);
        fitTCMInvYieldSectionEta->SetParError(0,fitTCMInvXSectionEta->GetParError(0)/xSection2760GeVINEL);
        fitTCMInvYieldSectionEta->SetParameter(2,fitTCMInvXSectionEta->GetParameter(2)/xSection2760GeVINEL);
        fitTCMInvYieldSectionEta->SetParError(2,fitTCMInvXSectionEta->GetParError(2)/xSection2760GeVINEL);
        fitTCMInvYieldSectionEta->Write("TwoComponentModelFitYieldEta");

        
        // writing Y shifted graphs in addition
        if (bWCorrection.Contains("Y")){
            if(graphPCMEtaInvXSectionStat_yShifted)graphPCMEtaInvXSectionStat_yShifted->Write("graphInvCrossSectionEtaPCM2760GeVStatErr_yShifted");
            if(graphPCMEtaInvXSectionSys_yShifted)graphPCMEtaInvXSectionSys_yShifted->Write("graphInvCrossSectionEtaPCM2760GeVSysErr_yShifted");
            // EMCAL component
            if(graphEMCALEtaInvXSectionStat_yShifted)graphEMCALEtaInvXSectionStat_yShifted->Write("graphInvCrossSectionEtaEMCAL2760GeVStatErr_yShifted");
            if(graphEMCALEtaInvXSectionSys_yShifted)graphEMCALEtaInvXSectionSys_yShifted->Write("graphInvCrossSectionEtaEMCAL2760GeVSysErr_yShifted");
            // PCM-EMC component
            if(graphPCMEMCALEtaInvXSectionStat_yShifted)graphPCMEMCALEtaInvXSectionStat_yShifted->Write("graphInvCrossSectionEtaPCMEMCAL2760GeVStatErr_yShifted");
            if(graphPCMEMCALEtaInvXSectionSys_yShifted)graphPCMEMCALEtaInvXSectionSys_yShifted->Write("graphInvCrossSectionEtaPCMEMCAL2760GeVSysErr_yShifted");
            // Final spectrum correlations Method A
            if(graphCombEtaInvXSectionTotA_yShifted)graphCombEtaInvXSectionTotA_yShifted->Write("graphInvCrossSectionEtaComb2760GeVATotErr_yShifted");
            if(graphCombEtaInvXSectionStatA_yShifted)graphCombEtaInvXSectionStatA_yShifted->Write("graphInvCrossSectionEtaComb2760GeVAStatErr_yShifted");
            if(graphCombEtaInvXSectionSysA_yShifted)graphCombEtaInvXSectionSysA_yShifted->Write("graphInvCrossSectionEtaComb2760GeVASysErr_yShifted");  
        }    

        histoPCMEtaToPi0Stat->Write("histoRatioEtaToPi0PCM2760GeVStatErr");
        graphPCMEtaToPi0Sys->Write("graphRatioEtaToPi0PCM2760GeVSysErr");
        histoEMCALEtaToPi0Stat->Write("histoRatioEtaToPi0EMCAL2760GeVStatErr");
        graphEMCALEtaToPi0Sys->Write("graphRatioEtaToPi0EMCAL2760GeVSysErr");
        histoPCMEMCALEtaToPi0Stat->Write("histoRatioEtaToPi0PCMEMCAL2760GeVStatErr");
        graphPCMEMCALEtaToPi0Sys->Write("graphRatioEtaToPi0PCMEMCAL2760GeVSysErr");
        graphCombEtaToPi0TotA->Write("graphRatioEtaToPi0Comb2760GeVTotErr");
        graphCombEtaToPi0StatA->Write("graphRatioEtaToPi0Comb2760GeVStatErr");
        graphCombEtaToPi0SysA->Write("graphRatioEtaToPi0Comb2760GeVSysErr");

        
        directoryEta->mkdir("Supporting");
        directoryEta->cd("Supporting");
            // Writing full correction factors
            graphPCMEtaAccTimesEff->Write("EtaCorrectionFactorPCM");
            graphEMCALEtaAccTimesEff->Write("EtaCorrectionFactorEMCAL");
            graphPCMEMCALEtaAccTimesEff->Write("EtaCorrectionFactorPCMEMCAL");
        
            graphPCMEtaMass->Write("EtaMassDataPCM");
            graphPCMEtaMassMC->Write("EtaMassMCPCM");
            graphEMCALEtaMass->Write("EtaMassDataEMCAL");
            graphEMCALEtaMassMC->Write("EtaMassMCEMCAL");
            graphPCMEMCALEtaMass->Write("EtaMassDataPCMEMCAL");
            graphPCMEMCALEtaMassMC->Write("EtaMassMCPCMEMCAL");

            graphPCMEtaFWHM->Write("EtaWidthDataPCM");
            graphPCMEtaFWHMMC->Write("EtaWidthMCPCM");
            graphEMCALEtaFWHM->Write("EtaWidthDataEMCAL");
            graphEMCALEtaFWHMMC->Write("EtaWidthMCEMCAL");
            graphPCMEMCALEtaFWHM->Write("EtaWidthDataPCMEMCAL");
            graphPCMEMCALEtaFWHMMC->Write("EtaWidthMCPCMEMCAL");
        
    fCombResults.Close();


 // **********************************************************************************************************************
 // ************************* Saving only fits to final results **********************************************************
 // **********************************************************************************************************************
    
    TString nameOutputCommonFileFitsOnly    = "";
    if (flagMerged == 0) 
        nameOutputCommonFileFitsOnly        = Form("FitsPaperPP2760GeV_%s_Haitao.root", dateForOutput.Data());
    else if (flagMerged == 1) 
        nameOutputCommonFileFitsOnly        = Form("FitsPaperPP2760GeV_%s_FrediV1Clusterizer.root", dateForOutput.Data());
    else if (flagMerged == 2)
        nameOutputCommonFileFitsOnly        = Form("FitsPaperPP2760GeV_%s_FrediV2Clusterizer.root", dateForOutput.Data());
    else if (flagMerged == 3) 
        nameOutputCommonFileFitsOnly        = Form("FitsPaperPP2760GeV_%s_FrediV1ClusterizerNLM1.root", dateForOutput.Data());
    else if (flagMerged == 4) 
        nameOutputCommonFileFitsOnly        = Form("FitsPaperPP2760GeV_%s_FrediV1ClusterizerNLM2.root", dateForOutput.Data());
    
    TFile fFitsResults(nameOutputCommonFileFitsOnly.Data(), "RECREATE");

         // fits for pi0
        fitInvXSectionPi0->Write("TsallisFitPi0");
        fitTCMInvXSectionPi0->Write("TwoComponentModelFitPi0");
           
        // fits for eta
        fitInvXSectionEta->Write("TsallisFitEta");
        fitTCMInvXSectionEta->Write("TwoComponentModelFitEta");
                
    fFitsResults.Close();

 // **********************************************************************************************************************
 // ************************* Saving comparison to comb for diff measurements ********************************************
 // **********************************************************************************************************************
    
    TString nameOutputCommonFileCompOnly    = Form("ComparisonsPaperPP2760GeV_%s.root", dateForOutput.Data());
        
    TFile fCompResults(nameOutputCommonFileCompOnly.Data(), "RECREATE");
        
        graphRatioPi0CombCombFitStatA->Write("Pi0_RatioCombToCombFit_Stat");
        graphRatioPi0CombCombFitSysA->Write("Pi0_RatioCombToCombFit_Syst");
        graphRatioPi0PCMCombFitStat->Write("Pi0_RatioPCMToCombFit_Stat");
        graphRatioPi0PCMCombFitSys->Write("Pi0_RatioPCMToCombFit_Syst");
        graphRatioPi0PHOSCombFitStat->Write("Pi0_RatioPHOSToCombFit_Stat");
        graphRatioPi0PHOSCombFitSys->Write("Pi0_RatioPHOSToCombFit_Syst");
        graphRatioPi0EMCALCombFitStat->Write("Pi0_RatioEMCToCombFit_Stat");
        graphRatioPi0EMCALCombFitSys->Write("Pi0_RatioEMCToCombFit_Syst");
        graphRatioPi0PCMEMCALCombFitStat->Write("Pi0_RatioPCMEMCToCombFit_Stat");
        graphRatioPi0PCMEMCALCombFitSys->Write("Pi0_RatioPCMEMCToCombFit_Syst");
        graphRatioPi0EMCALMergedCombFitStat->Write("Pi0_RatiomEMCToCombFit_Stat");
        graphRatioPi0EMCALMergedCombFitSys->Write("Pi0_RatiomEMCToCombFit_Syst");

        graphRatioEtaCombCombFitStatA->Write("Eta_RatioCombToCombFit_Stat");
        graphRatioEtaCombCombFitSysA->Write("Eta_RatioCombToCombFit_Syst");        
        graphRatioEtaPCMCombFitStat->Write("Eta_RatioPCMToCombFit_Stat");
        graphRatioEtaPCMCombFitSys->Write("Eta_RatioPCMToCombFit_Syst");
        graphRatioEtaEMCALCombFitStat->Write("Eta_RatioEMCToCombFit_Stat");
        graphRatioEtaEMCALCombFitSys->Write("Eta_RatioEMCToCombFit_Syst");
        graphRatioEtaPCMEMCALCombFitStat->Write("Eta_RatioPCMEMCToCombFit_Stat");
        graphRatioEtaPCMEMCALCombFitSys->Write("Eta_RatioPCMEMCToCombFit_Syst");
       
    fCompResults.Close();
    
}
    
