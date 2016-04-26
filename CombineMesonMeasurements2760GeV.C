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

void CombineMesonMeasurements2760GeV(   TString fileNamePCM = "", 
                                        TString fileNamePCMEMCAL = "", 
                                        TString fileNameEMCALLow = "",  
                                        TString fileNameEMCALmerged = "",  
                                        TString suffix = "eps", 
                                        TString isMC= "", 
                                        TString thesisPlots = "", 
                                        TString bWCorrection="X",
                                        Int_t  flagMerged = kFALSE
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
    
    TString fileNameTheory                      = "ExternalInput/TheoryCompilationPP.root";
    TString fileNamePHOS                        = "ExternalInput/PHOS/2.76TeV/LHC11a_PHOS_pi0_pp2760_noBWCorr_FDcorr_20140218.root";
    TString fileNamePCMEta                      = "FinalResults/data_PCMResultsFullCorrection_PP_NoBinShifting_usedTosvnEtapaper.root";
    TString fileNameChargedPionPP               = "ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_20_May_2015.root";
    TString fileNameChargedHadronPP             = "ExternalInput/UnidentifiedCharged/ChargedHadrinSpectraPP_20_May_2015.root";
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
    
    TString nameFinalResDat                     = Form("%s/CombinedResults%s_FitResults.dat",outputDir.Data(),bWCorrection.Data());
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
    TString  nameTrigger[6]                     = {"MB", "INT7", "EMC1", "EMC7", "EG2", "EG1"};
    
    Color_t  colorDet[11];
    Color_t  colorDetMC[11];
    Style_t  markerStyleDet[11];
    Style_t  markerStyleDetMC[11];
    Size_t   markerSizeDet[11];
    Size_t   markerSizeDetMC[11];

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
    

    //************************** Read data for PCM **************************************************
    TFile* filePCM                                          = new TFile(fileNamePCM.Data());
    TH1D* histoPCMNumberOfEvents                            = (TH1D*)filePCM->Get("histoNumberOfEvents2.76TeV");
    TDirectory* directoryPCMPi0                             = (TDirectory*)filePCM->Get("Pi02.76TeV"); 
        TH1D* histoPCMPi0Mass                               = (TH1D*)directoryPCMPi0->Get("MassPi0");
        TH1D* histoPCMPi0FWHMMeV                            = (TH1D*)directoryPCMPi0->Get("FWHMPi0MeV");
        TH1D* histoPCMPi0TrueMass                           = (TH1D*)directoryPCMPi0->Get("TrueMassPi0");
        TH1D* histoPCMPi0TrueFWHMMeV                        = (TH1D*)directoryPCMPi0->Get("TrueFWHMPi0MeV");
        TH1D* histoPCMPi0Acc                                = (TH1D*)directoryPCMPi0->Get("AcceptancePi0");
        TH1D* histoPCMPi0TrueEffPt                          = (TH1D*)directoryPCMPi0->Get("EfficiencyPi0");
        TH1D* histoPCMPi0InvXSectionStat                    = (TH1D*)directoryPCMPi0->Get("InvCrossSectionPi0");
        TGraphAsymmErrors* graphPCMPi0InvXSectionStat       = new TGraphAsymmErrors(histoPCMPi0InvXSectionStat);
        graphPCMPi0InvXSectionStat->RemovePoint(graphPCMPi0InvXSectionStat->GetN()-1);
        graphPCMPi0InvXSectionStat->RemovePoint(0);
        TGraphAsymmErrors* graphPCMPi0InvXSectionSysA       = (TGraphAsymmErrors*)directoryPCMPi0->Get("InvCrossSectionPi0SysA");
        TGraphAsymmErrors* graphPCMPi0InvXSectionSys        = (TGraphAsymmErrors*)directoryPCMPi0->Get("InvCrossSectionPi0Sys");
        graphPCMPi0InvXSectionSys->RemovePoint(graphPCMPi0InvXSectionSys->GetN()-1);
//      TGraphAsymmErrors* graphPCMPi0CorrYieldSysErr    = (TGraphAsymmErrors*)directoryPCMPi0->Get("Pi0SystError");
        
        TH1D* histoPCMPi0AccTimesEff                        = (TH1D*)histoPCMPi0TrueEffPt->Clone("histoPCMPi0AccTimesEff");
        histoPCMPi0AccTimesEff->Multiply(histoPCMPi0Acc);
    
    TFile* filePCMEta                                       = new TFile(fileNamePCMEta.Data());
    TDirectory* directoryPCMEta                             = (TDirectory*)filePCMEta->Get("Eta2.76TeV"); 
        TH1D* histoPCMEtaMass                               = (TH1D*)directoryPCMEta->Get("MassEta");
        TH1D* histoPCMEtaFWHMMeV                            = (TH1D*)directoryPCMEta->Get("FWHMEtaMeV");
        TH1D* histoPCMEtaTrueMass                           = (TH1D*)directoryPCMEta->Get("TrueMassEta");
        TH1D* histoPCMEtaTrueFWHMMeV                        = (TH1D*)directoryPCMEta->Get("TrueFWHMEtaMeV");
        TH1D* histoPCMEtaAcc                                = (TH1D*)directoryPCMEta->Get("AcceptanceEta");
        TH1D* histoPCMEtaTrueEffPt                          = (TH1D*)directoryPCMEta->Get("EfficiencyEta");
        TH1D* histoPCMEtaInvXSectionStat                    = (TH1D*)directoryPCMEta->Get("InvCrossSectionEta");
        TGraphAsymmErrors* graphPCMEtaInvXSectionStat       = new TGraphAsymmErrors(histoPCMEtaInvXSectionStat);
        graphPCMEtaInvXSectionStat->RemovePoint(0);

        TH1D* histoPCMEtaAccTimesEff                        = (TH1D*)histoPCMEtaTrueEffPt->Clone("histoPCMEtaAccTimesEff");
        histoPCMEtaAccTimesEff->Multiply(histoPCMEtaAcc);

        //      TGraphAsymmErrors* graphPCMEtaInvXSectionSysA     = (TGraphAsymmErrors*)directoryPCMEta->Get("InvCrossSectionEtaSysA");
        TGraphAsymmErrors* graphPCMEtaInvXSectionSys        = (TGraphAsymmErrors*)directoryPCMEta->Get("InvCrossSectionEtaSys");
//      TGraphAsymmErrors* graphPCMEtaCorrYieldSysErr    = (TGraphAsymmErrors*)directoryPCMEta->Get("EtaSystError");
        
    cout << "here" << endl;
    Int_t nEvtPCM                                           = histoPCMNumberOfEvents->GetBinContent(1);
    cout << "here" << endl;
    histoPCMPi0Mass->Scale(1000.);
    histoPCMPi0TrueMass->Scale(1000.);
    histoPCMEtaMass->Scale(1000.);
    histoPCMEtaTrueMass->Scale(1000.);
    cout << "here" << endl;
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
               
//      TH1D* histoPCMEMCALSignalPlusBGPi0                          = (TH1D*)directoryPCMEMCALPi0->Get(Form("InvMassSigPlusBG_PtBin%02d",binExInvMass[4]));
//      TH1D* histoPCMEMCALBGPi0                                    = (TH1D*)directoryPCMEMCALPi0->Get(Form("InvMassBG_PtBin%02d",binExInvMass[4]));
//      TH1D* histoPCMEMCALSignalPi0                                = (TH1D*)directoryPCMEMCALPi0->Get(Form("InvMassSig_PtBin%02d",binExInvMass[4]));
//      TH1D* histoPCMEMCALRemainingBGPi0                           = (TH1D*)histoPCMEMCALSignalPi0->Clone("histoPCMEMCALRemainingBGPi0");
//      TF1* fitPCMEMCALSignalPi0                                   = (TF1*)directoryPCMEMCALPi0->Get(Form("FitInvMassSig_PtBin%02d",binExInvMass[4]));
//      histoPCMEMCALSignalPi0->Fit(fitPCMEMCALSignalPi0,"QRME0");
//      for (Int_t i=0; i < 6; i++){
//          cout << fitPCMEMCALSignalPi0->GetParameter(i) << "\t +- " << fitPCMEMCALSignalPi0->GetParError(i) << endl;
//      }     
//      TF1*  fitPCMEMCALLinearBckPi0                               = new TF1("Linearpp","[0]+[1]*x",0.0,0.3);
//      fitPCMEMCALLinearBckPi0->SetParameter(0, fitPCMEMCALSignalPi0->GetParameter(4));
//      fitPCMEMCALLinearBckPi0->SetParameter(1, fitPCMEMCALSignalPi0->GetParameter(5));
//      TVirtualFitter * fitterPCMEMCAL                             = TVirtualFitter::GetFitter();
//      Int_t nFreeParPCMEMCAL                                      = fitPCMEMCALSignalPi0->GetNumberFreeParameters();
//      double * covMatrixPCMEMCAL                                  = fitterPCMEMCAL->GetCovarianceMatrix();
//      for (Int_t i = 1; i < histoPCMEMCALSignalPi0->GetXaxis()->FindBin(0.3); i++){
//          Double_t startBinEdge                                   = histoPCMEMCALSignalPi0->GetXaxis()->GetBinLowEdge(i);
//          Double_t endBinEdge                                     = histoPCMEMCALSignalPi0->GetXaxis()->GetBinUpEdge(i);
//          Double_t intLinearBack                                  = fitPCMEMCALLinearBckPi0->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
//          Double_t errorLinearBck                                 = pow(( pow( (endBinEdge-startBinEdge)*fitPCMEMCALSignalPi0->GetParError(4),2) + 
//                                                                          pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPCMEMCALSignalPi0->GetParError(5),2)
//                                                                          +2*covMatrixPCMEMCAL[nFreeParPCMEMCAL*nFreeParPCMEMCAL-2]*(endBinEdge-startBinEdge)*0.5*
//                                                                          (endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
//          histoPCMEMCALRemainingBGPi0->SetBinContent(i,intLinearBack);
//          histoPCMEMCALRemainingBGPi0->SetBinError(i,errorLinearBck);
//          cout << fitPCMEMCALLinearBckPi0->Eval(startBinEdge) << "\t" <<fitPCMEMCALLinearBckPi0->Eval(endBinEdge) << "\t" 
//               << histoPCMEMCALRemainingBGPi0->GetBinContent(i) << "\t" <<histoPCMEMCALSignalPi0->GetBinContent(i) << endl;
//      }
//      histoPCMEMCALSignalPi0->Add(histoPCMEMCALRemainingBGPi0,-1.);
//      fitPCMEMCALSignalPi0->SetParameter(4,0.);
//      fitPCMEMCALSignalPi0->SetParameter(5,0.);
        
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

        
        
//      TDirectory* directoryEMCALPi0Add                     = (TDirectory*)directoryEMCALPi0->Get("Slices_in_pT_Data"); 
//      TH1D* histoEMCALSignalPlusBGPi0                             = (TH1D*)directoryEMCALPi0Add->Get(Form("h1_FG_%d",binExInvMass[2]));
//      TH1D* histoEMCALTotalBGPi0                                  = (TH1D*)directoryEMCALPi0Add->Get(Form("h1_TotalBG_%d",binExInvMass[2]));
//      TH1D* histoEMCALBGPi0                                       = (TH1D*)directoryEMCALPi0Add->Get(Form("h1_BG_%d",binExInvMass[2]));
//      TH1D* histoEMCALSignalPi0                                   = (TH1D*)directoryEMCALPi0Add->Get(Form("h1_Su_%d",binExInvMass[2]));
//      TH1D* histoEMCALRemainingBGPi0                              = (TH1D*)directoryEMCALPi0Add->Get(Form("h1_CB_%d",binExInvMass[2]));
//      TF1* fitEMCALSignalPi0                                      = (TF1*)directoryEMCALPi0Add->Get(Form("f_Cryst_%d",binExInvMass[2]));
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
        
    //************************** Read data for EMCAL merged **************************************************    
    TFile* fileEMCALmerged                                  = new TFile(fileNameEMCALmerged.Data());
    TDirectory* directoryEMCALmergedPi0                     = (TDirectory*)fileEMCALmerged->Get("Pi02.76TeV"); 
//         TH1D* histoEMCALMergedPi0Eff                        = NULL;
//         TH1D* histoEMCALMergedPi0Pur                        = NULL;
        TH1D* histoEMCALMergedPi0AccTimesEff                = NULL;
        TH1D* histoEMCALMergedPi0InvXSectionStat            = NULL;
        TH1D* histoEMCALMergedPi0InvXSectionSys             = NULL;
        TGraphAsymmErrors* graphEMCALMergedPi0AccTimesEff       = NULL;
        TGraphAsymmErrors* graphEMCALMergedPi0InvXSectionStat   = NULL;
        TGraphAsymmErrors* graphEMCALMergedPi0InvXSectionSys    = NULL;
        if (flagMerged == 0){
            TH1D* histoEMCALMergedPi0Eff                    = (TH1D*)directoryEMCALmergedPi0->Get("Efficiencypi0");
//             histoEMCALMergedPi0Pur                          = (TH1D*)directoryEMCALmergedPi0->Get("PurityPi0");
            Double_t acceptanceEMCALmerged                  = 0.18;
            histoEMCALMergedPi0AccTimesEff                  = (TH1D*)histoEMCALMergedPi0Eff->Clone("histoEMCALMergedPi0AccTimesEff");
//          histoEMCALMergedPi0AccTimesEff->Multiply(histoEMCALMergedPi0Pur);
            histoEMCALMergedPi0AccTimesEff->Scale(acceptanceEMCALmerged);
            graphEMCALMergedPi0AccTimesEff                  = new TGraphAsymmErrors(histoEMCALMergedPi0AccTimesEff);
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
            graphEMCALMergedPi0InvXSectionStat              = (TGraphAsymmErrors*)directoryEMCALmergedPi0->Get("graphInvCrossSectionPi0");
            graphEMCALMergedPi0InvXSectionSys               = (TGraphAsymmErrors*)directoryEMCALmergedPi0->Get("InvCrossSectionPi0Sys");
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
            } else if (flagMerged == 2 || flagMerged == 3){
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
//         return;

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
        TGraphAsymmErrors* graphPHOSPi0AccTimesEff          = new TGraphAsymmErrors(histoPHOSPi0AccTimesEff);


    // *******************************************************************************************************
    // ************************** Loading theory calculations ************************************************
    // *******************************************************************************************************        
    TFile* fileTheoryCompilation                            = new TFile(fileNameTheory.Data());
        TH1F* histoPythia8InvXSection                       = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Spec2760GeV");
        TH1F* histoPythia8InvXSection_VarBinning            = (TH1F*) fileTheoryCompilation->Get("histoInvSecPythia8Spec2760GeVVarBinning");
        TGraph* graphNLOCalcPi0MuHalf                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuHalf2760GeV");
        TGraph* graphNLOCalcPi0MuOne                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuOne2760GeV");
        TGraph* graphNLOCalcPi0MuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecPi0MuTwo2760GeV");
        TGraph* graphNLOCalcEtaMuHalf                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuHalf2760GeV");
        TGraph* graphNLOCalcEtaMuOne                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuOne2760GeV");
        TGraph* graphNLOCalcEtaMuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcInvSecEtaMuTwo2760GeV");
        TGraph* graphEtaToPi0NLOMuHalf                      = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuHalf2760GeV");
        TGraph* graphEtaToPi0NLOMuOne                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuOne2760GeV");
        TGraph* graphEtaToPi0NLOMuTwo                       = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcEtaOverPi0MuTwo2760GeV");
        TGraph* graphNLODSSCalcMuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcDSSInvSecPi0MuTwo2760GeV");
        TGraph* graphNLOCGCCalcMuTwo                        = (TGraph*)fileTheoryCompilation->Get("graphNLOCalcCGCInvCrossSec2760GeV");
        TGraphAsymmErrors* graphNLODSS14Calc                = (TGraphAsymmErrors*)fileTheoryCompilation->Get("graphNLOCalcDSS14InvCrossSec2760GeV");
    
        while (graphNLOCalcEtaMuHalf->GetX()[graphNLOCalcEtaMuHalf->GetN()-1] > 22. )
            graphNLOCalcEtaMuHalf->RemovePoint(graphNLOCalcEtaMuHalf->GetN()-1);
        while (graphNLOCalcEtaMuOne->GetX()[graphNLOCalcEtaMuOne->GetN()-1] > 22. )
            graphNLOCalcEtaMuOne->RemovePoint(graphNLOCalcEtaMuOne->GetN()-1);
        while (graphNLOCalcEtaMuTwo->GetX()[graphNLOCalcEtaMuTwo->GetN()-1] > 22. )
            graphNLOCalcEtaMuTwo->RemovePoint(graphNLOCalcEtaMuTwo->GetN()-1);
            
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
        TGraphAsymmErrors* graphChargedHadronsStatPP        = (TGraphAsymmErrors*)fileChargedHadronsInputpp->Get("graphChargedHadronsStatPP2760GeV");
        graphChargedHadronsStatPP                           = ScaleGraph(graphChargedHadronsStatPP,0.5*1e9/xSection2760GeVINEL);
        TGraphAsymmErrors* graphChargedHadronsSysPP         = (TGraphAsymmErrors*)fileChargedHadronsInputpp->Get("graphChargedHadronsSysPP2760GeV");
        graphChargedHadronsSysPP                            = ScaleGraph(graphChargedHadronsSysPP,0.5*1e9/xSection2760GeVINEL);
        
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
    Int_t offSetsPi0[11]            = { 0,  2,  0,  0,  0,
                                        0,  0,  0,  0,  0, 
                                        0};
    Int_t offSetsPi0Sys[11]         = { 1,  3,  6,  0,  3,
                                        0,  0,  0,  0,  21, 
                                        0};
    if (flagMerged == 1 || flagMerged == 4) {
        offSetsPi0[9]               = 0;
        offSetsPi0Sys[9]            = 20;
    } else if (flagMerged == 2 || flagMerged == 3) {
        offSetsPi0[9]               = 0;
        offSetsPi0Sys[9]            = 21;
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
                                                                                                fileNamePi0OutputWeightingOld,1
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
    
        TLegend* legendAccWeightsOld    = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetPi0Old), 32);
        for (Int_t i = 0; i < nMeasSetPi0Old; i++){
            DrawGammaSetMarkerTGraph(graphWeightsPi0Old[availablePi0MeasOld[i]],
                                     markerStyleDet[availablePi0MeasOld[i]], 
                                     markerSizeDet[availablePi0MeasOld[i]]*0.5, 
                                     colorDet[availablePi0MeasOld[i]] , 
                                     colorDet[availablePi0MeasOld[i]]);
            graphWeightsPi0Old[availablePi0MeasOld[i]]->Draw("p,same,z");
            legendAccWeightsOld->AddEntry(graphWeightsPi0Old[availablePi0MeasOld[i]],nameMeasGlobal[availablePi0MeasOld[i]],"p");
        }    
        legendAccWeightsOld->Draw();

        TLatex *labelWeightsEnergy      = new TLatex(0.7,0.20,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
        labelWeightsEnergy->SetTextFont(43);
        labelWeightsEnergy->Draw();
        TLatex *labelWeightsPi0         = new TLatex(0.7,0.16,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelWeightsPi0, 0.85*textSizeLabelsPixel,4);
        labelWeightsPi0->SetTextFont(43);
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
    // ************************ Adding PCM-EMCal measurements to matrix & calculating relativ errors ************************
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
                                                                                                fileNamePi0OutputWeightingA, "2.76TeV", "Pi0", kTRUE
                                                                                            );
    
    if (graphCombPi0InvXSectionTotA == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
    }    
    graphCombPi0InvXSectionStatA->RemovePoint(0);
    graphCombPi0InvXSectionSysA->RemovePoint(0);
    graphCombPi0InvXSectionTotA->RemovePoint(0);
    graphCombPi0InvXSectionTotA->Print();
    
//     return;
    
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

        TLegend* legendAccWeights   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetPi0A), 32);
        for (Int_t i = 0; i < nMeasSetPi0A; i++){
            DrawGammaSetMarkerTGraph(graphWeightsPi0A[availablePi0MeasA[i]], markerStyleDet[availablePi0MeasA[i]], markerSizeDet[availablePi0MeasA[i]]*0.5, colorDet[availablePi0MeasA[i]] , colorDet[availablePi0MeasA[i]]);
            graphWeightsPi0A[availablePi0MeasA[i]]->Draw("p,same,z");
            legendAccWeights->AddEntry(graphWeightsPi0A[availablePi0MeasA[i]],nameMeasGlobal[availablePi0MeasA[i]],"p");
        }    
        legendAccWeights->Draw();

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
        SetStyleTLatex( labelRatioToOldEnergy, 0.85*textSizeLabelsPixel,4);
        labelRatioToOldEnergy->SetTextFont(43);
        labelRatioToOldEnergy->Draw();
        TLatex *labelRatioToOldPi0      = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToOldPi0, 0.85*textSizeLabelsPixel,4);
        labelRatioToOldPi0->SetTextFont(43);
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
            legendRelSysErr->AddEntry(sysErrorRelCollectionPi0[availablePi0MeasA[i]],nameMeasGlobal[availablePi0MeasA[i]],"p");
        }    
        legendRelSysErr->Draw();

        TLatex *labelRelSysErrEnergy    = new TLatex(0.15,0.89,collisionSystem2760GeV.Data());
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
            sysErrorRelCollectionPi0[availablePi0MeasA[i]]->Draw("p,same,z");
            legendRelStatErr->AddEntry(statErrorRelCollectionPi0[availablePi0MeasA[i]],nameMeasGlobal[availablePi0MeasA[i]],"p");
        }    
        legendRelStatErr->Draw();

        TLatex *labelRelStatErrEnergy   = new TLatex(0.75,0.89,collisionSystem2760GeV.Data());
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
    TGraphAsymmErrors* graphCombPi0InvXSectionRelTotOld      = CalculateRelErrUpAsymmGraph( graphCombPi0InvXSectionTotOld, "relativeTotalErrorPi0_Old");

    
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
        SetStyleTLatex( labelRelTotErrEnergy, 0.85*textSizeLabelsPixel,4);
        labelRelTotErrEnergy->SetTextFont(43);
        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrPi0       = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelTotErrPi0, 0.85*textSizeLabelsPixel,4);
        labelRelTotErrPi0->SetTextFont(43);
        labelRelTotErrPi0->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/Pi0_RelTotErr.%s",outputDir.Data(),suffix.Data()));
    
    
    histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,30.5);
    histo2DRelTotErrPi0->Draw("copy");
        graphCombPi0InvXSectionRelTotOld->Draw("p,same,z");
        graphCombPi0InvXSectionRelTotA->Draw("p,same,z");

        legendRelTotErr->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrPi0->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/Pi0_RelTotErrZoomed.%s",outputDir.Data(),suffix.Data()));

    histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,80.5);
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
    TGraphAsymmErrors* graphCombPi0InvXSectionTotAUnShifted     = (TGraphAsymmErrors*)graphCombPi0InvXSectionTotA->Clone("Pi0Unshifted"); 
    TGraphAsymmErrors* graphCombPi0InvXSectionStatAUnShifted    = (TGraphAsymmErrors*)graphCombPi0InvXSectionStatA->Clone("Pi0UnshiftedStat"); 
    TGraphAsymmErrors* graphCombPi0InvXSectionSysAUnShifted     = (TGraphAsymmErrors*)graphCombPi0InvXSectionSysA->Clone("Pi0UnshiftedSys"); 

    TGraphAsymmErrors* graphPCMPi0InvXSectionStatUnShifted      = (TGraphAsymmErrors*)graphPCMPi0InvXSectionStat->Clone("Pi0UnshiftedStatPCM"); 
    TGraphAsymmErrors* graphPCMPi0InvXSectionSysUnShifted       = (TGraphAsymmErrors*)graphPCMPi0InvXSectionSys->Clone("Pi0UnshiftedSysPCM"); 

    TGraphAsymmErrors* graphPHOSPi0InvXSectionStatUnshi         = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionStat->Clone("Pi0UnshiftedStatPHOS"); 
    TGraphAsymmErrors* graphPHOSPi0InvXSectionSysUnshi          = (TGraphAsymmErrors*)graphPHOSPi0InvXSectionSys->Clone("Pi0UnshiftedSysPHOS"); 
    
    TGraphAsymmErrors* graphEMCALPi0InvXSectionStatUnshi        = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionStat->Clone("Pi0UnshiftedStatEMCAL"); 
    TGraphAsymmErrors* graphEMCALPi0InvXSectionSysUnshi         = (TGraphAsymmErrors*)graphEMCALPi0InvXSectionSys->Clone("Pi0UnshiftedSysEMCAL"); 

    TGraphAsymmErrors* graphPCMEMCALPi0InvXSectionStatUnshi     = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionStat->Clone("Pi0UnshiftedStatPCMEMCAL"); 
    TGraphAsymmErrors* graphPCMEMCALPi0InvXSectionSysUnshi      = (TGraphAsymmErrors*)graphPCMEMCALPi0InvXSectionSys->Clone("Pi0UnshiftedSysPCMEMCAL"); 

    TGraphAsymmErrors* graphEMCALMergedPi0InvXSectionStatUnshi  = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionStat->Clone("Pi0UnshiftedStatEMCALMerged"); 
    TGraphAsymmErrors* graphEMCALMergedPi0InvXSectionSysUnshi   = (TGraphAsymmErrors*)graphEMCALMergedPi0InvXSectionSys->Clone("Pi0UnshiftedSysEMCALMerged"); 

    
    // Calculating binshifts
    Double_t paramGraphPi0[3]                              = {1.0e12, 8., 0.13};
    TF1* fitInvXSectionPi0                       = FitObject("l","fitInvCrossSectionPi02760GeV","Pi0",histoEMCALPi0InvXSectionStat,0.4,50.,paramGraphPi0,"QNRMEX0+");
    TF1* fitInvXSectionPi0Graph                  = (TF1*)fitInvXSectionPi0->Clone("fitInvCrossSectionPi02760GeVGraph"); 
    
    if(bWCorrection.CompareTo("X")==0 ){
        TF1* fitTsallisPi0PtMult                 = FitObject("tmpt","TsallisMultWithPtPi02760GeV","Pi0");
        fitTsallisPi0PtMult->SetParameters(paramGraphPi0[0],paramGraphPi0[1], paramGraphPi0[2]) ; // standard parameter optimize if necessary
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
        graphPCMPi0InvXSectionStat               = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA,
                                                                                    graphPCMPi0InvXSectionStat,
                                                                                    fitTsallisPi0PtMult, 
                                                                                    0, 16);
        graphPCMPi0InvXSectionSys                = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA, 
                                                                                    graphPCMPi0InvXSectionSys, 
                                                                                    fitTsallisPi0PtMult, 
                                                                                    0, 16);
        graphPHOSPi0InvXSectionStat              = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA, 
                                                                                    graphPHOSPi0InvXSectionStat, 
                                                                                    fitTsallisPi0PtMult,
                                                                                    2, 16);
        graphPHOSPi0InvXSectionSys               = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA, 
                                                                                    graphPHOSPi0InvXSectionSys, 
                                                                                    fitTsallisPi0PtMult,
                                                                                    2, 16);
        graphEMCALPi0InvXSectionStat             = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA, 
                                                                                    graphEMCALPi0InvXSectionStat, 
                                                                                    fitTsallisPi0PtMult,
                                                                                    5, 16);
        
        graphEMCALPi0InvXSectionSys              = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA, 
                                                                                    graphEMCALPi0InvXSectionSys, 
                                                                                    fitTsallisPi0PtMult,
                                                                                    5, 16);
        graphPCMEMCALPi0InvXSectionStat          = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA,
                                                                                    graphPCMEMCALPi0InvXSectionStat,
                                                                                    fitTsallisPi0PtMult, 
                                                                                    2, 17);
        graphPCMEMCALPi0InvXSectionSys           = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA, 
                                                                                    graphPCMEMCALPi0InvXSectionSys, 
                                                                                    fitTsallisPi0PtMult, 
                                                                                    2, 17);
        if (flagMerged == 0){
            graphEMCALMergedPi0InvXSectionStat       = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA, 
                                                                                        graphEMCALMergedPi0InvXSectionStat, 
                                                                                        fitTsallisPi0PtMult,
                                                                                        21, 12);
            graphEMCALMergedPi0InvXSectionSys        = ApplyXshiftIndividualSpectra(    graphCombPi0InvXSectionTotA, 
                                                                                        graphEMCALMergedPi0InvXSectionSys, 
                                                                                        fitTsallisPi0PtMult,
                                                                                        21, 12);
        }
                
        TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
        DrawGammaCanvasSettings( canvasDummy2,  0.1, 0.01, 0.015, 0.08);
        canvasDummy2->SetLogy();
        canvasDummy2->SetLogx();
        TH2F * histo2DDummy2;
        histo2DDummy2               = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.23,70.,1000,1e-1,10e11);
        SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
        histo2DDummy2->DrawCopy(); 
    
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatAUnShifted, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
        graphCombPi0InvXSectionStatAUnShifted->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionTotA, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombPi0InvXSectionTotA->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvXSectionSys, markerStyleDet[0] ,markerSizeDet[0], colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        graphPCMPi0InvXSectionSys->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphPHOSPi0InvXSectionSys, markerStyleDet[1] ,markerSizeDet[1], colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
        graphPHOSPi0InvXSectionSys->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0InvXSectionSys, markerStyleDet[2] ,markerSizeDet[2], colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        graphEMCALPi0InvXSectionSys->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0InvXSectionSys, markerStyleDet[4] ,markerSizeDet[4], colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
        graphPCMEMCALPi0InvXSectionSys->Draw("pEsame");

        fitInvXSectionPi0->SetLineColor(kBlue+2);
        fitInvXSectionPi0->Draw("same");
        
        canvasDummy2->Update();
        canvasDummy2->Print(Form("%s/ComparisonShiftedPi0_2760GeV.%s",outputDir.Data(),suffix.Data()));
        delete canvasDummy2;
    }
        
//     return;    
    graphCombPi0InvXSectionTotA->Fit(fitInvXSectionPi0,"QNRMEX0+","",0.4,50.);
    
    fitInvXSectionPi0        = FitObject("l","fitInvCrossSectionPi02760GeV","Pi0",graphCombPi0InvXSectionTotA,0.4,50.,paramGraphPi0,"QNRMEX0+");
    fitInvXSectionPi0        = FitObject("l","fitInvCrossSectionPi02760GeV","Pi0",graphCombPi0InvXSectionTotA,0.4,50. ,paramGraphPi0,"QNRMEX0+");

    cout << WriteParameterToFile(fitInvXSectionPi0)<< endl;
    
    Double_t paramTCMPi0[5]  = {graphCombPi0InvXSectionTotA->GetY()[2],0.3,graphCombPi0InvXSectionTotA->GetY()[2],0.8,2};
    cout << paramTCMPi0[0] << "\t" << paramTCMPi0[1] <<  "\t" << paramTCMPi0[2] << "\t" << paramTCMPi0[3] << "\t" << paramTCMPi0[4] << endl;
    TF1* fitTCMInvXSectionPi0       = FitObject("tcm","fitTCMInvCrossSectionPi02760GeV","Pi0",graphCombPi0InvXSectionStatA,0.4,50.,paramTCMPi0,"QNRMEX0+","", kFALSE);
//     fitTCMInvXSectionPi0->SetParLimits(1,0.2,0.4);
//     fitTCMInvXSectionPi0->SetParLimits(3,0.6,1.);
//     fitTCMInvXSectionPi0->SetParLimits(4,1.,4.);
//     graphCombPi0InvXSectionStatA->Fit(fitTCMInvXSectionPi0,"QNRMEX0+","",0.4,50.);
//     fitTCMInvXSectionPi0     = FitObject("tcm","fitTCMInvCrossSectionPi02760GeV","Pi0",graphCombPi0InvXSectionStatA,0.4,50. ,paramTCMPi0,"QNRMEX0+");
        
    TF1* fitTCMDecomposedL                 = new TF1("twoCompModel_DecLow",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mesonMassExpectPi0,mesonMassExpectPi0,mesonMassExpectPi0),0.4,2);
    fitTCMDecomposedL->SetParameters(graphCombPi0InvXSectionTotA->GetY()[2],0.3);
    graphCombPi0InvXSectionStatA->Fit(fitTCMDecomposedL,"QNRMEX0+","",0.4,2.);
    TF1 *fitTCMDecomposedH                 = new TF1("twoCompModel_DecH","[0]/(TMath::Power(1+x*x/([1]*[1]*[2]),[2]))",4,50);
    graphCombPi0InvXSectionTotA->Fit(fitTCMDecomposedH,"QNRMEX0+","",4,50);
    fitTCMDecomposedH->SetParameters(graphCombPi0InvXSectionTotA->GetY()[2],0.8, 2);
        
//     cout << WriteParameterToFile(fitTCMInvXSectionPi0)<< endl;
//     cout << WriteParameterToFile(fitTCMDecomposedL)<< endl;
//     cout << WriteParameterToFile(fitTCMDecomposedH)<< endl;
    
    Double_t paramTCMPi0New[5]  = { fitTCMDecomposedL->GetParameter(0),fitTCMDecomposedL->GetParameter(1),
                                    fitTCMDecomposedH->GetParameter(0),fitTCMDecomposedH->GetParameter(1),fitTCMDecomposedH->GetParameter(2)};
    fitTCMInvXSectionPi0        = FitObject("tcm","fitTCMInvCrossSectionPi02760GeV","Pi0",graphCombPi0InvXSectionStatA,0.4,50. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);

    TF1* fitTCMInvXSectionPi02  = FitObject("tcm","fitTCMInvCrossSectionPi02760GeV2","Pi0",graphCombPi0InvXSectionStatA,0.4,50.,paramTCMPi0New,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitTCMInvXSectionPi02)<< endl;

    cout << WriteParameterToFile(fitTCMInvXSectionPi0)<< endl;
    
//     return;
//  TString forOutput= WriteParameterToFile(fitInvXSectionPi0);
//  cout<< forOutput.Data()<< endl;
//  

    TH1D* histoRatioPythia8ToFit                     = (TH1D*) histoPythia8InvXSection->Clone(); 
    histoRatioPythia8ToFit                           = CalculateHistoRatioToFit (histoRatioPythia8ToFit, fitTCMInvXSectionPi0); 
    TH1D* histoRatioPythia8VarBinningToFit           = (TH1D*) histoPythia8InvXSection_VarBinning->Clone(); 
    histoRatioPythia8VarBinningToFit                 = CalculateHistoRatioToFit (histoRatioPythia8VarBinningToFit, fitTCMInvXSectionPi0); 

    TGraph* graphRatioPi0CombNLOMuHalf               = (TGraph*)graphNLOCalcPi0MuHalf->Clone();
    TGraph* graphRatioPi0CombNLOMuOne                = (TGraph*)graphNLOCalcPi0MuOne->Clone();
    TGraph* graphRatioPi0CombNLOMuTwo                = (TGraph*)graphNLOCalcPi0MuTwo->Clone();
    TGraphAsymmErrors* graphRatioPi0CombNLODSS14     = (TGraphAsymmErrors*)graphNLODSS14Calc->Clone();
    graphRatioPi0CombNLOMuHalf                       = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuHalf, fitTCMInvXSectionPi0); 
    graphRatioPi0CombNLOMuOne                        = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuOne, fitTCMInvXSectionPi0); 
    graphRatioPi0CombNLOMuTwo                        = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuTwo, fitTCMInvXSectionPi0); 
    graphRatioPi0CombNLODSS14                        = CalculateGraphErrRatioToFit(graphRatioPi0CombNLODSS14, fitTCMInvXSectionPi0); 
    
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
    // ******************************************* Ratio of Comb to Fit ****************************************
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
    histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
    histo2DPi0RatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioPi0CombCombFitSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioPi0CombCombFitStatA->Draw("p,same,z");

        DrawGammaLines(0.23, 70. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.23, 70. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.9, 0.9,0.1, kGray, 7);

        TLatex *labelRatioToFitEnergy   = new TLatex(0.665, 0.91, collisionSystem2760GeV.Data());
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

    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToCombFit_PP2760GeV.%s",outputDir.Data(),suffix.Data()));
    
    // **********************************************************************************************************************
    // ******************************************* Ratio of Individual meas to Fit ******************************************
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

        DrawGammaLines(0.23, 70. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.23, 70. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.9, 0.9,0.5, kGray, 7);
        
        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitPi0->Draw();
    
        //****************************** Definition of the Legend ******************************************
        //**************** Row def ************************
        Double_t rowsLegendOnlyPi0Ratio[6]          = {0.91, 0.855, 0.805, 0.755, 0.705, 0.655};
        Double_t rowsLegendOnlyPi0RatioAbs[6]       = {0.91, 2.165, 2.035, 1.905, 1.765, 1.635};
        Double_t columnsLegendOnlyPi0Ratio[3]       = {0.115, 0.375, 0.46};
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
//  for (Int_t i = 0; i< maxNBinsEtaW0Merged; i++){
//      cout << xPtLimitsEtaWOMerged[i] << ", ";
//  }
//  cout << endl;
    
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Eta                    
    Int_t offSetsEta[11]            = { 0,  0,  0,  0,  0,
                                        0,  0,  0,  0,  0, 
                                        0};
    Int_t offSetsEtaSys[11]         = { 1,  0,  4,  0,  2,
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
                                                                                                fileNameEtaOutputWeightingA,"2.76TeV", "Eta", kFALSE
                                                                                            );
    
    
//    
    
    graphCombEtaInvXSectionStatA->RemovePoint(0);
    graphCombEtaInvXSectionSysA->RemovePoint(0);
    graphCombEtaInvXSectionTotA->RemovePoint(0);
    graphCombEtaInvXSectionStatA->RemovePoint(graphCombEtaInvXSectionStatA->GetN()-1);
    graphCombEtaInvXSectionSysA->RemovePoint(graphCombEtaInvXSectionSysA->GetN()-1);
    graphCombEtaInvXSectionTotA->RemovePoint(graphCombEtaInvXSectionTotA->GetN()-1);

    graphCombEtaInvXSectionStatA->Print();
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
    histo2DEtaWeights = new TH2F("histo2DEtaWeights","histo2DEtaWeights",11000,0.33,25.,1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DEtaWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DEtaWeights->GetXaxis()->SetMoreLogLabels();
    histo2DEtaWeights->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DEtaWeights->GetYaxis()->SetRangeUser(-10,10);    
    histo2DEtaWeights->Draw("copy");

        TLegend* legendAccWeightsEta   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEtaA), 32);
        for (Int_t i = 0; i < nMeasSetEtaA; i++){
            DrawGammaSetMarkerTGraph(graphWeightsEtaA[availableEtaMeasA[i]], markerStyleDet[availableEtaMeasA[i]], markerSizeDet[availableEtaMeasA[i]]*0.5, colorDet[availableEtaMeasA[i]] , colorDet[availableEtaMeasA[i]]);
            graphWeightsEtaA[availableEtaMeasA[i]]->Draw("p,same,z");
            legendAccWeightsEta->AddEntry(graphWeightsEtaA[availableEtaMeasA[i]],nameMeasGlobal[availableEtaMeasA[i]],"p");
        }    
        legendAccWeightsEta->Draw();

        labelWeightsEnergy->Draw();
        TLatex *labelWeightsEta         = new TLatex(0.7,0.16,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelWeightsEta, 0.85*textSizeLabelsPixel,4);
        labelWeightsEta->SetTextFont(43);
        labelWeightsEta->Draw();

//      DrawGammaLines(0.33, 25. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.33, 25. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.33, 25. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.33, 25. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.33, 25. , 0.2, 0.2,0.1, kGray, 3);
        
    canvasWeights->SaveAs(Form("%s/Eta_WeightsA.%s",outputDir.Data(),suffix.Data()));


    // *********************************************************************************************************************
    // ************************************ Visualize relative errors Eta ******************************************************
    // *********************************************************************************************************************
    
    canvasRelSysErr->cd();
   
    TH2F * histo2DRelSysErrEta;
    histo2DRelSysErrEta                 = new TH2F("histo2DRelSysErrEta","histo2DRelSysErrEta",11000,0.33,25.,1000,0,80.5);
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
            legendRelSysErrEta->AddEntry(sysErrorRelCollectionEta[availableEtaMeasA[i]],nameMeasGlobal[availableEtaMeasA[i]],"p");
        }    
        legendRelSysErrEta->Draw();

        labelRelSysErrEnergy->Draw();
        TLatex *labelRelSysErrEta       = new TLatex(0.15,0.85,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelSysErrEta, 0.85*textSizeLabelsPixel,4);
        labelRelSysErrEta->SetTextFont(43);
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
    histo2DRelStatErrEta                = new TH2F("histo2DRelStatErrEta","histo2DRelStatErrEta",11000,0.33,25.,1000,0,80.5);
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
            legendRelStatErrEta->AddEntry(statErrorRelCollectionEta[availableEtaMeasA[i]],nameMeasGlobal[availableEtaMeasA[i]],"p");
        }    
        legendRelStatErrEta->Draw();

        labelRelStatErrEnergy->Draw();
        TLatex *labelRelStatErrEta      = new TLatex(0.75,0.85,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelStatErrEta, 0.85*textSizeLabelsPixel,4);
        labelRelStatErrEta->SetTextFont(43);
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
    histo2DRelTotErrEta                 = new TH2F("histo2DRelTotErrEta","histo2DRelTotErrEta",11000,0.33,25.,1000,0,80.5);
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
        SetStyleTLatex( labelRelTotErrEta, 0.85*textSizeLabelsPixel,4);
        labelRelTotErrEta->SetTextFont(43);
        labelRelTotErrEta->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/Eta_RelTotErr.%s",outputDir.Data(),suffix.Data()));
        
    histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,50.5);
    histo2DRelTotErrEta->Draw("copy");
        graphCombEtaInvXSectionRelTotA->Draw("p,same,z");
    
        legendRelTotErr->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrEta->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/Eta_RelTotErrZoomed.%s",outputDir.Data(),suffix.Data()));
    
    histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,80.5);
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
    
    // Calculating binshifts
    Double_t paramGraphEta[3]                    = {1.0e10, 7., 0.13};
    TF1* fitInvXSectionEta                       = FitObject("l","fitInvCrossSectionEta2760GeV","Eta",graphCombEtaInvXSectionTotAUnshi,0.5,20.,paramGraphEta,"QNRMEX0+");
    TF1* fitInvXSectionEtaGraph                  = (TF1*)fitInvXSectionEta->Clone("fitInvCrossSectionEta2760GeVGraph"); 
    
    Double_t paramTCMEta[5]  = {graphCombEtaInvXSectionTotA->GetY()[0],0.5,graphCombEtaInvXSectionTotA->GetY()[0]/1000,0.8,3};
    
    if(bWCorrection.CompareTo("X")==0 ){
        TF1* fitTsallisEtaPtMult                 = FitObject("tmpt","TsallisMultWithPtEta2760GeV","Eta");
        fitTsallisEtaPtMult->SetParameters(paramGraphEta[0],paramGraphEta[1], paramGraphEta[2]) ; // standard parameter optimize if necessary
//         fitTsallisEtaPtMult->SetParameters(paramTCMEta[0],paramTCMEta[1], paramTCMEta[2], paramTCMEta[3], paramTCMEta[4]) ; // standard parameter optimize if necessary
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
                                                                                    0, 6, "Eta");
        graphPCMEtaInvXSectionSys                = ApplyXshiftIndividualSpectra(    graphCombEtaInvXSectionTotA, 
                                                                                    graphPCMEtaInvXSectionSys, 
                                                                                    fitTsallisEtaPtMult, 
                                                                                    0, 6, "Eta");
        graphEMCALEtaInvXSectionStat             = ApplyXshiftIndividualSpectra(    graphCombEtaInvXSectionTotA, 
                                                                                    graphEMCALEtaInvXSectionStat, 
                                                                                    fitTsallisEtaPtMult,
                                                                                    3, 8, "Eta");
        graphEMCALEtaInvXSectionSys              = ApplyXshiftIndividualSpectra(    graphCombEtaInvXSectionTotA, 
                                                                                    graphEMCALEtaInvXSectionSys, 
                                                                                    fitTsallisEtaPtMult,
                                                                                    3, 8, "Eta");
        graphPCMEMCALEtaInvXSectionStat          = ApplyXshiftIndividualSpectra(    graphCombEtaInvXSectionTotA,
                                                                                    graphPCMEMCALEtaInvXSectionStat,
                                                                                    fitTsallisEtaPtMult, 
                                                                                    1, 9, "Eta");
        graphPCMEMCALEtaInvXSectionSys           = ApplyXshiftIndividualSpectra(    graphCombEtaInvXSectionTotA, 
                                                                                    graphPCMEMCALEtaInvXSectionSys, 
                                                                                    fitTsallisEtaPtMult, 
                                                                                    1, 9, "Eta");

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
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvXSectionSys, markerStyleDet[0] ,markerSizeDet[0], colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        graphPCMEtaInvXSectionSys->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaInvXSectionSys, markerStyleDet[2] ,markerSizeDet[2], colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        graphEMCALEtaInvXSectionSys->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaInvXSectionSys, markerStyleDet[4] ,markerSizeDet[4], colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
        graphPCMEMCALEtaInvXSectionSys->Draw("pEsame");

        fitInvXSectionEta->SetLineColor(kBlue+2);
        fitInvXSectionEta->Draw("same");
        
        canvasDummy2->Update();
        canvasDummy2->Print(Form("%s/ComparisonShiftedEta_2760GeV.%s",outputDir.Data(),suffix.Data()));
        delete canvasDummy2;
        delete histo2DDummy3;
    }
//     
    graphCombEtaInvXSectionTotA->Fit(fitInvXSectionEta,"QNRMEX0+","",0.5,20.);
    
    fitInvXSectionEta        = FitObject("l","fitInvCrossSectionEta2760GeV","Eta",graphCombEtaInvXSectionTotA,0.5,20.,paramGraphEta,"QNRMEX0+");
    fitInvXSectionEta        = FitObject("l","fitInvCrossSectionEta2760GeV","Eta",graphCombEtaInvXSectionTotA,0.5,20. ,paramGraphEta,"QNRMEX0+");

    cout << WriteParameterToFile(fitInvXSectionEta)<< endl;
    
//     Double_t paramTCMEta[5]  = {graphCombEtaInvXSectionTotA->GetY()[0],0.5,graphCombEtaInvXSectionTotA->GetY()[0]/1000,0.8,3};
    TF1* fitTCMInvXSectionEta= FitObject("tcm","fitTCMInvCrossSectionEta2760GeV","Eta",graphCombEtaInvXSectionTotA,0.5,20.,paramTCMEta,"QNRMEX0+","", kFALSE);
    fitTCMInvXSectionEta     = FitObject("tcm","fitTCMInvCrossSectionEta2760GeV","Eta",graphCombEtaInvXSectionTotA,0.5,20. ,paramTCMEta,"QNRMEX0+","", kFALSE);
//  TF1* fitTCMDecomposedL                 = new TF1("twoCompModel_DecLow",Form("[0]*exp(-(TMath::Sqrt(x*x+%.10f*%.10f)-%.10f)/[1])",mesonMassExpectEta,mesonMassExpectEta,mesonMassExpectEta));
//  fitTCMDecomposedL->SetParameters(fitTCMInvXSectionEta->GetParameter(0),fitTCMInvXSectionEta->GetParameter(1));
//  TF1 *fitTCMDecomposedH                 = new TF1("twoCompModel_DecH","[0]/(1 + x*x/TMath::Power(1+x*x/([1]*[1]*[2]),-[2]) )");
// //      graphCombEtaInvXSectionTotA->Fit(fitTCMDecomposedH,"QNRMEX0+","",5,20);
//  fitTCMDecomposedH->SetParameters(fitTCMInvXSectionEta->GetParameter(2),fitTCMInvXSectionEta->GetParameter(3), fitTCMInvXSectionEta->GetParameter(4));
    cout << WriteParameterToFile(fitTCMInvXSectionEta)<< endl;
//     return;
    
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
    textSizeLabelsPixel                 = 48;
    canvasRatioToCombFit->cd();
    TH2F * histo2DEtaRatioToCombFit;
    histo2DEtaRatioToCombFit               = new TH2F("histo2DEtaRatioToCombFit","histo2DEtaRatioToCombFit",1000,0.33,25.,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DEtaRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
    histo2DEtaRatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DEtaRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
    histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.55);
    histo2DEtaRatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSysA, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioEtaCombCombFitSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatA, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioEtaCombCombFitStatA->Draw("p,same,z");

        DrawGammaLines(0.33, 25. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.33, 25. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.33, 25. , 0.9, 0.9,0.1, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        TLatex *labelRatioToFitEta      = new TLatex(0.852,0.82,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitEta, textSizeLabelsPixel,4);
        labelRatioToFitEta->SetTextFont(43);
        labelRatioToFitEta->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfCombToCombFit_PP2760GeV.%s",outputDir.Data(),suffix.Data()));
    
    // **********************************************************************************************************************
    // ******************************************* Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************
    
    canvasRatioToCombFit->cd();
    histo2DEtaRatioToCombFit->Draw("copy");
    
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

        DrawGammaLines(0.33, 25. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.33, 25. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.33, 25. , 0.9, 0.9,0.5, kGray, 7);
        
        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitEta->Draw();
    
        //****************************** Definition of the Legend ******************************************
        Double_t rowsLegendOnlyEtaRatio[4]          = {0.91, 0.86, 0.81, 0.76};
        Double_t rowsLegendOnlyEtaRatioAbs[4]       = {0.91, 2.245, 2.11, 1.975};
        Double_t columnsLegendOnlyEtaRatio[3]       = {0.115, 0.325, 0.40};
        Double_t columnsLegendOnlyEtaRatioAbs[3]    = {0.115, 1.2, 1.62};

        //****************** first Column **************************************************
        TLatex *textPCMOnlyRatioEta                 = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[1],"PCM");
        SetStyleTLatex( textPCMOnlyRatioEta, textSizeLabelsPixel,4);
        textPCMOnlyRatioEta->SetTextFont(43);
        textPCMOnlyRatioEta->Draw();
        TLatex *textEMCALOnlyRatioEta               = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[2],"EMCal");
        SetStyleTLatex( textEMCALOnlyRatioEta, textSizeLabelsPixel,4);
        textEMCALOnlyRatioEta->SetTextFont(43);
        textEMCALOnlyRatioEta->Draw();
        TLatex *textPCMEMCALOnlyRatioEta            = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[3],"PCM-EMCal");
        SetStyleTLatex( textPCMEMCALOnlyRatioEta, textSizeLabelsPixel,4);
        textPCMEMCALOnlyRatioEta->SetTextFont(43);
        textPCMEMCALOnlyRatioEta->Draw();
        
        //****************** second Column *************************************************
        TLatex *textStatOnlyRatioEta                = new TLatex(columnsLegendOnlyEtaRatio[1],rowsLegendOnlyEtaRatio[0] ,"stat");
        SetStyleTLatex( textStatOnlyRatioEta, textSizeLabelsPixel,4);
        textStatOnlyRatioEta->SetTextFont(43);
        textStatOnlyRatioEta->Draw();
        TLatex *textSysOnlyRatioEta                 = new TLatex(columnsLegendOnlyEtaRatio[2] ,rowsLegendOnlyEtaRatio[0],"syst");
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

        DrawGammaSetMarker(histoPCMPi0FWHMMeV, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        histoPCMPi0FWHMMeV->Draw("p,same,e");
        DrawGammaSetMarker(histoPCMPi0TrueFWHMMeV, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
        histoPCMPi0TrueFWHMMeV->Draw("p,same,e");

        TLatex *labelLegendAMass    = new TLatex(0.13,0.06,"a)");
        SetStyleTLatex( labelLegendAMass, textSizeLabelsPixel,4);
        labelLegendAMass->SetTextFont(43);
        labelLegendAMass->Draw();

        TLatex *labelMassPerf       = new TLatex(0.13,0.87,"ALICE performance");
        SetStyleTLatex( labelMassPerf, textSizeLabelsPixel,4);
        labelMassPerf->SetTextFont(43);
        labelMassPerf->Draw();        
        TLatex *labelMassEnergy     = new TLatex(0.13,0.78,collisionSystem2760GeV.Data());
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
        
        DrawGammaSetMarker(histoPCMPi0Mass, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        histoPCMPi0Mass->Draw("p,same,e");
        DrawGammaSetMarker(histoPCMPi0TrueMass, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
        histoPCMPi0TrueMass->Draw("p,same,e");

        DrawGammaLines(0.23, 25. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);

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

        DrawGammaSetMarker(histoPCMPi0FWHMMeV, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        histoPCMPi0FWHMMeV->Draw("p,same,e");
        DrawGammaSetMarker(histoPCMPi0TrueFWHMMeV, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
        histoPCMPi0TrueFWHMMeV->Draw("p,same,e");

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

        DrawGammaSetMarker(histoPHOSPi0Mass, markerStyleDet[1], markerSizeDet[1]*0.55, colorDet[1] , colorDet[1]);
        histoPHOSPi0Mass->Draw("p,same,e");
        DrawGammaSetMarker(histoPHOSPi0TrueMass, markerStyleDetMC[1], markerSizeDetMC[1]*0.55, colorDetMC[1] , colorDetMC[1]);
        histoPHOSPi0TrueMass->Draw("p,same,e");
        
        graphEMCALPi0Mass->Draw("p,same,z");
        graphEMCALPi0MassMC->Draw("p,same,z");
        
        graphPCMEMCALPi0Mass->Draw("p,same,z");
        graphPCMEMCALPi0MassMC->Draw("p,same,z");
        
        histoPCMPi0Mass->Draw("p,same,e");
        histoPCMPi0TrueMass->Draw("p,same,e");

        DrawGammaLines(0.23, 25. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);

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

        
        TH2F * histo2DAllEtaFWHM    = new TH2F("histo2DAllEtaFWHM","histo2DAllEtaFWHM", 20, 0.33, 25. ,1000., -5, 92);
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

        TH2F * histo2DAllEtaMass    = new TH2F("histo2DAllEtaMass","histo2DAllEtaMass",20, 0.33, 25., 1000., 500., 565);
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
        
        DrawGammaSetMarker(histoPCMEtaMass, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        histoPCMEtaMass->Draw("p,same,e");
        DrawGammaSetMarker(histoPCMEtaTrueMass, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
        histoPCMEtaTrueMass->Draw("p,same,e");

        DrawGammaLines(0.33, 25. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,0.3, kGray);
        
        labelLegendBMass->Draw();
        
    canvasMassWidthEta->Update();
    canvasMassWidthEta->Print(Form("%s/Eta_MassAndWidth.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************** Cross section for pi0 single measurement 2.76TeV ************************************
    // **********************************************************************************************************************
    
    TCanvas* canvasXSectionPi0  = new TCanvas("canvasXSectionPi0","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasXSectionPi0, 0.14, 0.02, 0.02, 0.09);
    canvasXSectionPi0->SetLogx();
    canvasXSectionPi0->SetLogy();
    
    TH2F * histo2DXSectionPi0;
    histo2DXSectionPi0          = new TH2F("histo2DXSectionPi0","histo2DXSectionPi0",11000,0.23,70.,1000,2e-2,10e11);
    SetStyleHistoTH2ForGraphs(histo2DXSectionPi0, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 1.,1.45);
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

        
        TLatex *labelEnergyXSectionPi0      = new TLatex(0.64,0.92,collisionSystem2760GeV.Data());
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
  
    TCanvas* canvasXSectionEta      = new TCanvas("canvasXSectionEta","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasXSectionEta, 0.14, 0.02, 0.02, 0.09);
    canvasXSectionEta->SetLogx();
    canvasXSectionEta->SetLogy();
    
    TH2F * histo2DXSectionEta;
    histo2DXSectionEta              = new TH2F("histo2DXSectionEta","histo2DXSectionEta",11000,0.33,25.,1000,2e1,10e10);
    SetStyleHistoTH2ForGraphs(histo2DXSectionEta, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )",0.035,0.04, 0.035,0.04, 1.,1.45);
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


        TLatex *labelEnergyXSectionEta      = new TLatex(0.64,0.92,collisionSystem2760GeV.Data());
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
        legendXSectionEta->AddEntry(graphEMCALEtaInvXSectionSys,"EMCal","fp");
        legendXSectionEta->AddEntry(graphPCMEMCALEtaInvXSectionSys,"PCM-EMCal","fp");
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
        histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, 0.23,  70, 1000, 8e-6, 2e-1 );
        SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{T} (GeV/#it{c})", "A #times #epsilon_{eff} ",  
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.1);//(#times #epsilon_{pur})
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

        histoPCMPi0AccTimesEff->Draw("p,same,e");
        graphPCMEMCALPi0AccTimesEff->Draw("p,same,z");
        graphEMCALPi0AccTimesEff->Draw("p,same,z");
        graphEMCALMergedPi0AccTimesEff->Draw("p,same,e");

        TLegend* legendEffiAccPi02          = GetAndSetLegend2(0.55, 0.13, 0.83, 0.13+(5*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffiAccPi02->AddEntry(histoPCMPi0AccTimesEff,"PCM","p");
        legendEffiAccPi02->AddEntry(graphPCMEMCALPi0AccTimesEff,"PCM-EMCal","p");
        legendEffiAccPi02->AddEntry(graphEMCALPi0AccTimesEff,"EMCal","p");
        legendEffiAccPi02->AddEntry(graphEMCALMergedPi0AccTimesEff,"EMCal, merged","p");
        legendEffiAccPi02->AddEntry(graphPHOSPi0AccTimesEff,"PHOS","p");
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
        histo2DAccEffEta                = new TH2F("histo2DAccEffEta", "histo2DAccEffEta",1000, 0.33,  25, 1000, 8e-5, 1.5e-1 );
        SetStyleHistoTH2ForGraphs( histo2DAccEffEta, "#it{p}_{T} (GeV/#it{c})", "A #times #epsilon_{eff} ",  
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.1);//(#times #epsilon_{pur})
        histo2DAccEffEta->GetYaxis()->SetLabelOffset(0.001);
        histo2DAccEffEta->GetXaxis()->SetLabelOffset(-0.01);
        histo2DAccEffEta->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DAccEffEta->DrawCopy(); 
    
        DrawGammaSetMarker(histoPCMEtaAccTimesEff, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        histoPCMEtaAccTimesEff->Draw("p,same,e");
        
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaAccTimesEff, markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
        graphPCMEMCALEtaAccTimesEff->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaAccTimesEff, markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
        graphEMCALEtaAccTimesEff->Draw("p,same,z");
        
        TLegend* legendEffiAccEta           = GetAndSetLegend2(0.62, 0.13, 0.9, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffiAccEta->AddEntry(histoPCMEtaAccTimesEff,"PCM","p");
        legendEffiAccEta->AddEntry(graphPCMEMCALEtaAccTimesEff,"PCM-EMCal","p");
        legendEffiAccEta->AddEntry(graphEMCALEtaAccTimesEff,"EMCal","p");
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
        SetStyleHistoTH2ForGraphs( histo2DTriggerEffPi0, "#it{p}_{T} (GeV/#it{c})", "#epsilon_{Trigger}",  
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


        TLatex *labelPerfTriggEff               = new TLatex(0.6,0.425,"ALICE simulation");
        SetStyleTLatex( labelPerfTriggEff, textSizeLabelsRel,4);
        labelPerfTriggEff->Draw();
        TLatex *labelEnergyTriggEff             = new TLatex(0.6,0.378,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelEnergyTriggEff, textSizeLabelsRel,4);
        labelEnergyTriggEff->Draw();
        TLatex *labelDetSysTriggEffPi0          = new TLatex(0.6,0.33,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysTriggEffPi0, textSizeLabelsRel,4);
        labelDetSysTriggEffPi0->Draw();

        TLegend* legendTriggEffPi0              = GetAndSetLegend2(0.415, 0.13, 0.9, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel,3);
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

        TLatex *labelPCMEMCALTriggEff           = new TLatex(0.6,0.275,"PCM-EMCal");
        SetStyleTLatex( labelPCMEMCALTriggEff, textSizeLabelsRel,4);
        labelPCMEMCALTriggEff->Draw();
        TLatex *labelEMCALTriggEff              = new TLatex(0.84,0.275,"EMCal");
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
        TLatex *labelDetSysTriggEffEta          = new TLatex(0.6,0.33,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysTriggEffEta, 0.85*textSizeLabelsRel,4);
        labelDetSysTriggEffEta->Draw();
        labelPCMEMCALTriggEff->Draw();
        labelEMCALTriggEff->Draw();
        legendTriggEffPi0->Draw();

        
    canvasTriggerEff->Update();
    canvasTriggerEff->Print(Form("%s/Eta_TriggerEff.%s",outputDir.Data(),suffix.Data()));
    
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
//      ratio2DTheoryPP->GetYaxis()->SetLabelOffset(0.01);
        ratio2DTheoryPP->GetXaxis()->SetLabelFont(42);
        ratio2DTheoryPP->GetYaxis()->SetLabelFont(42);
//      ratio2DTheoryPP->GetXaxis()->SetLabelOffset(-0.01);
    //    ratio2DTheoryPP->GetXaxis()->SetTickLength(0.07);
        ratio2DTheoryPP->DrawCopy(); 

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS14, 0, 0, kGreen-6, kGreen-6, widthLinesBoxes, kTRUE, kGreen-6);
        graphRatioPi0CombNLODSS14->Draw("2,same");
        
        DrawGammaSetMarker(histoRatioPythia8ToFit, 24, 1.5, kRed+2 , kRed+2);  
        histoRatioPythia8ToFit->SetLineWidth(widthCommonFit);
//      histoRatioPythia8ToFit->GetXaxis()->SetRangeUser(0.5,14);
        histoRatioPythia8ToFit->Draw("same,hist,c");
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

        
        DrawGammaLines(0.23, 70.,1., 1.,0.1,kGray);

        TLegend* legendRatioTheorypp_3Parted= GetAndSetLegend2(0.15,0.66,0.4,0.92, 0.85* textSizeLabelsPixel);
    //    legendRatioTheorypp_3Parted->SetNColumns(2);
        legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0CombCombFitSysA,"#pi^{0} ALICE","pf");
        legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0CombNLOMuHalf, "NLO, DSS07 #mu = 0.5 #it{p}_{T}", "l");  
        legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0CombNLOMuOne,  "NLO, DSS07 #mu = #it{p}_{T}", "l");  
        legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0CombNLOMuTwo,  "NLO, DSS07 #mu = 2 #it{p}_{T}", "l");  
        legendRatioTheorypp_3Parted->AddEntry(graphRatioPi0CombNLODSS14,  "NLO, DSS14", "f");  
        legendRatioTheorypp_3Parted->AddEntry(histoRatioPythia8ToFit,  "Pythia 8, Tune 4C", "l");  
        legendRatioTheorypp_3Parted->Draw();
        
        TLatex *labelRatioTheoryPP   = new TLatex(0.15,0.93,collisionSystem2760GeV.Data());
        SetStyleTLatex( labelRatioTheoryPP, 0.85*textsizeLabelsPP,4);
        labelRatioTheoryPP->Draw();
    
        
    canvasRatioPP->Update();
    canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PP.%s",outputDir.Data(),suffix.Data()));

    canvasRatioPP->cd();
        ratio2DTheoryPP->GetYaxis()->SetRangeUser(0.05,2.45);
        ratio2DTheoryPP->DrawCopy(); 

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS14, 0, 0, kGreen-6, kGreen-6, widthLinesBoxes, kTRUE, kGreen-6);
        graphRatioPi0CombNLODSS14->Draw("2,same");
        
        DrawGammaSetMarker(histoRatioPythia8ToFit, 24, 1.5, kRed+2 , kRed+2);  
        histoRatioPythia8ToFit->SetLineWidth(widthCommonFit);
//      histoRatioPythia8ToFit->GetXaxis()->SetRangeUser(0.5,14);
        histoRatioPythia8ToFit->Draw("same,hist,c");
        
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatA->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA->Draw("2,same");
        graphRatioPi0CombCombFitStatA->Draw("p,same");

        DrawGammaLines(0.23, 70.,1., 1.,0.1,kGray);

        TLegend* legendRatioTheoryPPNew= GetAndSetLegend2(0.15,0.76,0.4,0.92, 0.85*textSizeLabelsPixel);
        legendRatioTheoryPPNew->AddEntry(graphRatioPi0CombCombFitSysA,"#pi^{0} ALICE","pf");
        legendRatioTheoryPPNew->AddEntry(graphRatioPi0CombNLODSS14,  "NLO, DSS14", "f");  
        legendRatioTheoryPPNew->AddEntry(histoRatioPythia8ToFit,  "Pythia 8, Tune 4C", "l");  
        legendRatioTheoryPPNew->Draw();
        
        labelRatioTheoryPP->Draw();

        
    canvasRatioPP->Update();
    canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PP_NewTheory.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************************* Comparison to theory calculations Eta ************************************
    // **********************************************************************************************************************    
    canvasRatioPP->cd();
    canvasRatioPP->SetLogx();

        TH2F * ratio2DTheoryPPEta       = new TH2F("ratio2DTheoryPPEta","ratio2DTheoryPPEta",1000,0.33,25.,1000,0.2,3.6);
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

        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, kGreen-6);
        graphRatioEtaCombNLOMuHalf->Draw("same,c");
        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuOne, widthCommonFit, styleLineNLOMuOne, kGreen-6);
        graphRatioEtaCombNLOMuOne->Draw("same,c");
        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, kGreen-6);
        graphRatioEtaCombNLOMuTwo->Draw("same,c");
        
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioEtaCombCombFitStatA->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioEtaCombCombFitSysA->SetLineWidth(0);
        graphRatioEtaCombCombFitSysA->Draw("2,same");
        graphRatioEtaCombCombFitStatA->Draw("p,same");

        DrawGammaLines(0.23, 70.,1., 1.,0.1,kGray);

        TLegend* legendRatioTheoryppEta_3Parted= GetAndSetLegend2(0.15,0.92-(0.85*textsizeLabelsPP*4),0.4,0.92, 0.85* textSizeLabelsPixel);
        legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaCombCombFitSysA,"#eta ALICE","pf");
        legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaCombNLOMuHalf, "NLO, AESS #mu = 0.5 #it{p}_{T}", "l");  
        legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaCombNLOMuOne,  "NLO, AESS #mu = #it{p}_{T}", "l");  
        legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaCombNLOMuTwo,  "NLO, AESS #mu = 2 #it{p}_{T}", "l");  
//         legendRatioTheoryppEta_3Parted->AddEntry(graphRatioEtaCombNLODSS14,  "NLO, DSS14", "f");  
//         legendRatioTheoryppEta_3Parted->AddEntry(histoRatioPythia8ToFit,  "Pythia 8, Tune 4C", "l");  
        legendRatioTheoryppEta_3Parted->Draw();
        labelRatioTheoryPP->Draw();
            
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

        
        DrawGammaSetMarker(histoPythia8InvXSection, 24, 1.5, kRed+2 , kRed+2);  
        histoPythia8InvXSection->SetLineWidth(widthCommonFit);
        histoPythia8InvXSection->Draw("same,hist,c");

        graphNLODSS14Calc->RemovePoint(0);
        DrawGammaSetMarkerTGraphAsym(graphNLODSS14Calc, 0, 0, kGreen-6, kGreen-6, widthLinesBoxes, kTRUE, kGreen-6);
        graphNLODSS14Calc->Draw("3,same");
        
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombPi0InvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvXSectionStatA, markerStyleComb, markerSizeComb, kBlack, kBlack);
        graphCombPi0InvXSectionStatA->Draw("p,same,z");

        fitTCMInvXSectionPi0->SetLineColor(kGray+2);
        fitTCMInvXSectionPi0->SetLineStyle(7);
        fitTCMInvXSectionPi0->SetLineWidth(2);
        fitTCMInvXSectionPi0->Draw("same");        
                
        TLatex *labelEnergyXSectionPaper= new TLatex(0.665, 0.91, collisionSystem2760GeV.Data());
        SetStyleTLatex( labelEnergyXSectionPaper, textsizeLabelsXSecUp,4);
        labelEnergyXSectionPaper->Draw();
        TLatex *labelALICEXSectionPaper= new TLatex(0.848,0.87,"ALICE");
        SetStyleTLatex( labelALICEXSectionPaper, textsizeLabelsXSecUp,4);
        labelALICEXSectionPaper->Draw();
        TLatex *labelDetSysXSectionPaper= new TLatex(0.824,0.83,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPaper, textsizeLabelsXSecUp,4);
        labelDetSysXSectionPaper->Draw();
        
        TLegend* legendXsectionPaper    = GetAndSetLegend2(0.17, 0.14, 0.5, 0.14+0.05*3, textSizeLabelsPixel);
        legendXsectionPaper->SetNColumns(1);
        legendXsectionPaper->SetMargin(0.2);
        legendXsectionPaper->AddEntry(graphCombPi0InvXSectionSysA,"Data","pf");
        legendXsectionPaper->AddEntry(graphNLODSS14Calc,"NLO, DSS14 ","f");
        legendXsectionPaper->AddEntry(histoPythia8InvXSection,"Pythia 8.1, Tune 4C","l");
        legendXsectionPaper->Draw();

        TLegend* legendXsectionPaper2     = GetAndSetLegend2(0.17, 0.06, 0.5, 0.12, textSizeLabelsPixel);        
        legendXsectionPaper2->SetMargin(0.2);
        legendXsectionPaper2->AddEntry(fitTCMInvXSectionPi0,"#it{A}_{e} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{n}","l");
        legendXsectionPaper2->Draw();
        
        
    padInvSectionNLORatio->cd();
    padInvSectionNLORatio->SetLogx(1);
        TH2F * ratio2DNLO               = new TH2F("ratio2DNLO","ratio2DNLO",1000,0.23,70.,1000,0.2,2.45);
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

        graphRatioPi0CombNLODSS14->RemovePoint(0);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS14, 0, 0, kGreen-6, kGreen-6, widthLinesBoxes, kTRUE, kGreen-6);
        graphRatioPi0CombNLODSS14->Draw("3,same");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatA->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA->Draw("2,same");
        graphRatioPi0CombCombFitStatA->Draw("p,same");
        
        DrawGammaLines(0.23, 70.,1., 1.,0.1,kGray);
        
    padInvSectionPythiaRatio->cd();
    padInvSectionPythiaRatio->SetLogx(1);
        TH2F * ratio2DPythia            = new TH2F("ratio2DPythia","ratio2DPythia",1000,0.23,70.,1000,0.2,2.45);
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
//      histoRatioPythia8ToFit->GetXaxis()->SetRangeUser(0.5,14);
        histoRatioPythia8ToFit->Draw("same,hist,c");
        
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatA->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSysA->SetLineWidth(0);
        graphRatioPi0CombCombFitSysA->Draw("2,same");
        graphRatioPi0CombCombFitStatA->Draw("p,same");

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

//         DrawGammaSetMarker(histoPythia8InvXSection, 24, 1.5, kRed+2 , kRed+2);  
//         histoPythia8InvXSection->SetLineWidth(widthCommonFit);
//         histoPythia8InvXSection->Draw("same,hist,c");

        DrawGammaNLOTGraph( graphNLOCalcEtaMuHalf, widthCommonFit, styleLineNLOMuHalf, kGreen-6);
        graphNLOCalcEtaMuHalf->Draw("same,c");
        DrawGammaNLOTGraph( graphNLOCalcEtaMuOne, widthCommonFit, styleLineNLOMuOne, kGreen-6);
        graphNLOCalcEtaMuOne->Draw("same,c");
        DrawGammaNLOTGraph( graphNLOCalcEtaMuTwo, widthCommonFit, styleLineNLOMuTwo, kGreen-6);
        graphNLOCalcEtaMuTwo->Draw("same,c");
        
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombEtaInvXSectionSysA->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvXSectionStatA, markerStyleComb, markerSizeComb, kBlack, kBlack);
        graphCombEtaInvXSectionStatA->Draw("p,same,z");

        fitTCMInvXSectionEta->SetLineColor(kGray+2);
        fitTCMInvXSectionEta->SetLineStyle(7);
        fitTCMInvXSectionEta->SetLineWidth(2);
        fitTCMInvXSectionEta->Draw("same");
        
        labelEnergyXSectionPaper->Draw();
        labelALICEXSectionPaper->Draw();
        TLatex *labelDetSysXSectionPaperEta = new TLatex(0.84,0.83,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPaperEta, textsizeLabelsXSecUp,4);
        labelDetSysXSectionPaperEta->Draw();
        
        TLegend* legendXsectionPaperEta     = GetAndSetLegend2(0.17, 0.14, 0.5, 0.14+0.05*4, textSizeLabelsPixel);
        legendXsectionPaperEta->SetNColumns(1);
        legendXsectionPaperEta->SetMargin(0.2);
        legendXsectionPaperEta->AddEntry(graphCombPi0InvXSectionSysA,"Data","pf");
        legendXsectionPaperEta->AddEntry(graphNLOCalcEtaMuHalf,"NLO, AESS, #mu = 0.5 #it{p}_{T}","l");
        legendXsectionPaperEta->AddEntry(graphNLOCalcEtaMuOne,"NLO, AESS, #mu = #it{p}_{T}","l");
        legendXsectionPaperEta->AddEntry(graphNLOCalcEtaMuTwo,"NLO, AESS, #mu = 2 #it{p}_{T}","l");
        legendXsectionPaperEta->Draw();

        legendXsectionPaper2->Draw();
        
    padInvSectionNLORatio->cd();
    padInvSectionNLORatio->SetLogx(1);
        TH2F * ratio2DNLOEta                = new TH2F("ratio2DNLOEta","ratio2DNLOEta",1000,0.33,25.,1000,0.2,3.45);
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

        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, kGreen-6);
        graphRatioEtaCombNLOMuHalf->Draw("same,c");
        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuOne, widthCommonFit, styleLineNLOMuOne, kGreen-6);
        graphRatioEtaCombNLOMuOne->Draw("same,c");
        DrawGammaNLOTGraph( graphRatioEtaCombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, kGreen-6);
        graphRatioEtaCombNLOMuTwo->Draw("same,c");
        
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioEtaCombCombFitStatA->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioEtaCombCombFitSysA->SetLineWidth(0);
        graphRatioEtaCombCombFitSysA->Draw("2,same");
        graphRatioEtaCombCombFitStatA->Draw("p,same");
        
        DrawGammaLines(0.33, 25.,1., 1.,0.1,kGray);
        
    padInvSectionPythiaRatio->cd();
    padInvSectionPythiaRatio->SetLogx(1);
        TH2F * ratio2DPythiaEta             = new TH2F("ratio2DPythiaEta","ratio2DPythiaEta",1000,0.33,25.,1000,0.2,3.45);
        SetStyleHistoTH2ForGraphs(ratio2DPythiaEta, "#it{p}_{T} (GeV/#it{c})","#frac{Pythia, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown, 
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

//         DrawGammaSetMarker(histoRatioPythia8ToFit, 24, 1.5, kRed+2 , kRed+2);  
//         histoRatioPythia8ToFit->SetLineWidth(widthCommonFit);
//         histoRatioPythia8ToFit->Draw("same,hist,c");
        
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioEtaCombCombFitStatA->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSysA, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioEtaCombCombFitSysA->SetLineWidth(0);
        graphRatioEtaCombCombFitSysA->Draw("2,same");
        graphRatioEtaCombCombFitStatA->Draw("p,same");

        DrawGammaLines(0.33, 25.,1., 1.,0.1,kGray);
        
    canvasInvSectionPaper->Print(Form("%s/Eta_InvXSectionWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));
    
    //*************************************************************************************************************
    //***************************** Comparison to Charged pions ***************************************************
    //*************************************************************************************************************

    TGraphAsymmErrors* graphCombPi0InvYieldStatA     = (TGraphAsymmErrors*) graphCombPi0InvXSectionStatA->Clone("graphCombPi0InvYieldStatA");
    graphCombPi0InvYieldStatA                        = ScaleGraph(graphCombPi0InvYieldStatA,1./xSection2760GeVINEL);
    TGraphAsymmErrors* graphCombPi0InvYieldSysA      = (TGraphAsymmErrors*) graphCombPi0InvXSectionSysA->Clone("graphCombPi0InvYieldSysA");
    graphCombPi0InvYieldSysA                         = ScaleGraph(graphCombPi0InvYieldSysA,1./xSection2760GeVINEL);

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
   
    TGraphErrors* graphChPiInvYieldCMSStatPPCMSCombUp               = NULL;
    TGraphErrors* graphChPiInvYieldCMSSystPPCMSCombUp               = NULL;
    TGraphErrors* graphCombPi0InvYieldStatRebinnedCMSComb    = NULL;
    TGraphErrors* graphCombPi0InvYieldSysRebinnedCMSComb     = NULL;
    TGraphErrors* graphRatioPi0CMSChPisCombA    = CalculateRatioBetweenSpectraWithDifferentBinning( graphCombPi0InvYieldStatA, graphCombPi0InvYieldSysA, 
                                                                                                        histoChPiInvYieldLowPtStatCMS, histoChPiInvYieldLowPtSysCMS,  
                                                                                                        kTRUE,  kTRUE, 
                                                                                                        &graphCombPi0InvYieldStatRebinnedCMSComb, &graphCombPi0InvYieldSysRebinnedCMSComb,
                                                                                                        &graphChPiInvYieldCMSStatPPCMSCombUp, &graphChPiInvYieldCMSSystPPCMSCombUp ) ;

       cout << "combined Spectrum - published" << endl;
    TGraphErrors* graphChPiInvYieldPubStatPPPubCombUp               = NULL;
    TGraphErrors* graphChPiInvYieldPubSystPPPubCombUp               = NULL;
    TGraphErrors* graphCombPi0InvYieldStatRebinnedPubComb    = NULL;
    TGraphErrors* graphCombPi0InvYieldSysRebinnedPubComb     = NULL;
    TGraphErrors* graphRatioPi0PubChPisCombA    = CalculateRatioBetweenSpectraWithDifferentBinning( graphCombPi0InvYieldStatA, graphCombPi0InvYieldSysA, 
                                                                                                        histoChPiInvYieldPubStatPP, histoChPiInvYieldPubSystPP,  
                                                                                                        kTRUE,  kTRUE, 
                                                                                                        &graphCombPi0InvYieldStatRebinnedPubComb, &graphCombPi0InvYieldSysRebinnedPubComb, 
                                                                                                        &graphChPiInvYieldPubStatPPPubCombUp, &graphChPiInvYieldPubSystPPPubCombUp )    ;

       cout << "combined Spectrum - published" << endl;
    TGraphErrors* graphChHadInvYieldPubStatPPHadCombUp              = NULL;
    TGraphErrors* graphChHadInvYieldPubSystPPHadCombUp              = NULL;
    TGraphErrors* graphCombPi0InvYieldStatRebinnedHadComb    = NULL;
    TGraphErrors* graphCombPi0InvYieldSysRebinnedHadComb     = NULL;
    TGraphErrors* graphRatioPi0PubChHadsCombA   = CalculateRatioBetweenSpectraWithDifferentBinning( graphCombPi0InvYieldStatA, graphCombPi0InvYieldSysA, 
                                                                                                    graphChargedHadronsStatPP, graphChargedHadronsSysPP,  
                                                                                                    kTRUE,  kTRUE, 
                                                                                                    &graphCombPi0InvYieldStatRebinnedHadComb, &graphCombPi0InvYieldSysRebinnedHadComb, 
                                                                                                    &graphChHadInvYieldPubStatPPHadCombUp, &graphChHadInvYieldPubSystPPHadCombUp )    ;
                                                                                                    
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

        TLegend* legendPi0CompChargedPionsPP    = GetAndSetLegend2(0.15, 0.8, 0.9, 0.90, 0.85* textSizeLabelsPixel);
        legendPi0CompChargedPionsPP->SetNColumns(2);
        legendPi0CompChargedPionsPP->SetMargin(0.12);
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
        TLegend* legendPi0CompChargedPionsPP2   = GetAndSetLegend2(0.15, 0.86, 0.9, 0.94, 0.85* textSizeLabelsPixel);
        legendPi0CompChargedPionsPP2->SetNColumns(2);
        legendPi0CompChargedPionsPP2->SetMargin(0.12);
        legendPi0CompChargedPionsPP2->AddEntry(graphRatioPi0PubChPisCombA,"#pi^{0}/#pi^{#pm}  (ALICE)","p");
        legendPi0CompChargedPionsPP2->AddEntry(graphRatioPi0CMSChPisCombA,"#pi^{0}/#pi^{#pm}  (#pi^{#pm} from CMS)","p");
        legendPi0CompChargedPionsPP2->Draw();
    
        labelRatioTheoryPP->Draw();
        
        DrawGammaLines(0., 70 , 1, 1 ,1, kGray, 1);   
   
    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedToNeutralPublishedCharged_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    canvasCompYieldPPInd->cd();

    histo2DCompCombinedRatio2->GetYaxis()->SetTitle("#pi^{0}/h^{#pm}");
    histo2DCompCombinedRatio2->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioPi0PubChHadsCombA, markerStyleCombHighPt, markerSizeComparison, kBlack , kBlack);
        graphRatioPi0PubChHadsCombA->Draw("E1psame");

        graphRatioPi0PubChHadsCombA->Print();
        TLegend* legendPi0CompChargedPionsPP3   = GetAndSetLegend2(0.15, 0.86, 0.9, 0.94, 0.85* textSizeLabelsPixel);
        legendPi0CompChargedPionsPP3->SetNColumns(2);
        legendPi0CompChargedPionsPP3->SetMargin(0.12);
        legendPi0CompChargedPionsPP3->AddEntry(graphRatioPi0PubChHadsCombA,"#pi^{0}/h^{#pm}  (ALICE)","p");
//      legendPi0CompChargedPionsPP3->AddEntry(graphRatioPi0CMSChPisCombA,"#pi^{0}/#pi^{#pm}  (#pi^{#pm} from CMS)","p");
        legendPi0CompChargedPionsPP3->Draw();
    
        labelRatioTheoryPP->Draw();
        
        DrawGammaLines(0., 70 , 1, 1 ,1, kGray, 1);   
   
    canvasCompYieldPPInd->Update();
    canvasCompYieldPPInd->Print(Form("%s/ComparisonChargedHadronToNeutralPublishedCharged_PP2760GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
    
//  // **********************************************************************************************************************
//  // **************************Plot example invariant mass bin PCM -EMCAL *************************************************
//  // **********************************************************************************************************************
//  
//  Double_t arrayBoundariesX1_InvMass[2];
//  Double_t arrayBoundariesY1_InvMass[3];
//  Double_t relativeMarginsXInvMass[3];
//  Double_t relativeMarginsYInvMass[3];
//  textSizeLabelsPixel                 = 90;
//  ReturnCorrectValuesForCanvasScaling(2500,3200, 1, 2,0.1, 0.005, 0.005,0.06,arrayBoundariesX1_InvMass,arrayBoundariesY1_InvMass,relativeMarginsXInvMass,relativeMarginsYInvMass);
// 
//  TCanvas* canvasInvMassSamplePlot    = new TCanvas("canvasInvMassSamplePlot","",0,0,2500,3200);  // gives the page size
//  DrawGammaCanvasSettings( canvasInvMassSamplePlot,  0.13, 0.02, 0.03, 0.06);
// 
//  TPad* padInvMassSampleWithBG        = new TPad("padInvMassSampleWithBG", "", arrayBoundariesX1_InvMass[0], arrayBoundariesY1_InvMass[1], arrayBoundariesX1_InvMass[1], arrayBoundariesY1_InvMass[0],-1, -1, -2);
//  DrawGammaPadSettings( padInvMassSampleWithBG, relativeMarginsXInvMass[0], relativeMarginsXInvMass[2], relativeMarginsYInvMass[0], relativeMarginsYInvMass[1]);
//  padInvMassSampleWithBG->Draw();
// 
//  TPad* padInvMassSampleSignal        = new TPad("padInvMassSampleSignal", "", arrayBoundariesX1_InvMass[0], arrayBoundariesY1_InvMass[2], arrayBoundariesX1_InvMass[1], arrayBoundariesY1_InvMass[1],-1, -1, -2);
//  DrawGammaPadSettings( padInvMassSampleSignal, relativeMarginsXInvMass[0], relativeMarginsXInvMass[2], relativeMarginsYInvMass[1], relativeMarginsYInvMass[2]);
//  padInvMassSampleSignal->Draw();
// 
//  Style_t markerStyleInvMassSGBG      = 20;
//  Style_t markerStyleInvMassSGBG2     = 24;
//  Size_t markerSizeInvMassSGBG        = 3;
//  Color_t markerColorInvMassSGBG      = kBlue+2;
//  Color_t markerColorInvMassSGBG2     = kBlue-6;
//  Style_t markerStyleInvMassMBG       = 24;
//  Size_t markerSizeInvMassMBG         = 1.5;
//  Color_t markerColorInvMassMBG       = kRed+2;
//  Style_t markerStyleInvMassCBG       = 24;
//  Size_t markerSizeInvMassCBG         = 1.5;
//  Color_t markerColorInvMassCBG       = kGreen+2;
//  Style_t markerStyleInvMassBG        = 20;
//  Size_t markerSizeInvMassBG          = 2;
//  Color_t markerColorInvMassBG        = kBlack;
//  Style_t markerStyleInvMassSG        = 20;
//  Size_t markerSizeInvMassSG          = 3;
//  Color_t markerColorInvMassSG        = kBlack;
//  Color_t fitColorInvMassSG           = kAzure+1;
//      
//  padInvMassSampleWithBG->cd();
//      Double_t marginInvMass          = relativeMarginsXInvMass[0]*2500;
//      Double_t textsizeLabelsInvMassUp = 0;
//      Double_t textsizeFacInvMassUp   = 0;
//      if (padInvMassSampleWithBG->XtoPixel(padInvMassSampleWithBG->GetX2()) < padInvMassSampleWithBG->YtoPixel(padInvMassSampleWithBG->GetY1())){
//          textsizeLabelsInvMassUp     = (Double_t)textSizeLabelsPixel/padInvMassSampleWithBG->XtoPixel(padInvMassSampleWithBG->GetX2()) ;
//          textsizeFacInvMassUp        = (Double_t)1./padInvMassSampleWithBG->XtoPixel(padInvMassSampleWithBG->GetX2()) ;
//      } else {
//          textsizeLabelsInvMassUp     = (Double_t)textSizeLabelsPixel/padInvMassSampleWithBG->YtoPixel(padInvMassSampleWithBG->GetY1());
//          textsizeFacInvMassUp        = (Double_t)1./padInvMassSampleWithBG->YtoPixel(padInvMassSampleWithBG->GetY1());
//      }
//      cout << textsizeLabelsInvMassUp << endl;
//      
//      TH2F * histo2DInvMassDummyUp;
//      histo2DInvMassDummyUp           = new TH2F("histo2DInvMassDummyUp","histo2DInvMassDummyUp",11000,0.05,0.225,20100,-100,20000);
//      SetStyleHistoTH2ForGraphs(histo2DInvMassDummyUp, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","counts",0.85*textsizeLabelsInvMassUp, textsizeLabelsInvMassUp,
//                                0.85*textsizeLabelsInvMassUp, textsizeLabelsInvMassUp,1., 0.14/(textsizeFacInvMassUp*marginInvMass));
// //      histo2DInvMassDummyUp->GetXaxis()->SetLabelOffset(-0.01);
//      histo2DInvMassDummyUp->GetYaxis()->SetRangeUser(-20.,345);
//      histo2DInvMassDummyUp->GetXaxis()->SetTickLength(0.04);
//      histo2DInvMassDummyUp->GetYaxis()->SetTickLength(0.02);
//      histo2DInvMassDummyUp->DrawCopy();
//      
//      DrawGammaSetMarker(histoPCMEMCALSignalPlusBGPi0, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
//      histoPCMEMCALSignalPlusBGPi0->Draw("same");
//      TH1D* histoPCMEMCALTotalBGPi0 = (TH1D*)histoPCMEMCALBGPi0->Clone("histoPCMEMCALTotalBGPi0");
//      histoPCMEMCALTotalBGPi0->Add(histoPCMEMCALRemainingBGPi0);
//      DrawGammaSetMarker(histoPCMEMCALTotalBGPi0, markerStyleInvMassBG, markerSizeInvMassBG, markerColorInvMassBG, markerColorInvMassBG);
//      histoPCMEMCALTotalBGPi0->Draw("same");
//      DrawGammaSetMarker(histoPCMEMCALBGPi0, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
//      histoPCMEMCALBGPi0->Draw("same");
//      DrawGammaSetMarker(histoPCMEMCALRemainingBGPi0, markerStyleInvMassCBG, markerSizeInvMassCBG, markerColorInvMassCBG, markerColorInvMassCBG);
//      histoPCMEMCALRemainingBGPi0->Draw("same");
//  
//      TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9,collisionSystem2760GeV.Data());
//      SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
//      labelInvMassEnergy->SetTextFont(43);
//      labelInvMassEnergy->Draw();
// 
//      TLatex *labelInvMassRecoPCMEMC  = new TLatex(0.135,0.84,"#gamma's rec. PCM, EMCal");
//      SetStyleTLatex( labelInvMassRecoPCMEMC, 0.85*textSizeLabelsPixel,4);
//      labelInvMassRecoPCMEMC->SetTextFont(43);
//      labelInvMassRecoPCMEMC->Draw();
//      
//      TLatex *labelInvMassPtRangePCMEMCAL = new TLatex(0.135,0.78,Form("%3.1f GeV/#it{c} < #it{p}_{T} < %3.1f GeV/#it{c}",rangePtExMin[4],rangePtExMax[4]));
//      SetStyleTLatex( labelInvMassPtRangePCMEMCAL, 0.85*textSizeLabelsPixel,4);
//      labelInvMassPtRangePCMEMCAL->SetTextFont(43);
//      labelInvMassPtRangePCMEMCAL->Draw();
// 
//      TLegend* legendInvMassPCMEMCAL  = GetAndSetLegend2(0.64, 0.74, 0.9, 0.98, 0.85*textSizeLabelsPixel);
//      legendInvMassPCMEMCAL->SetMargin(0.25);
//      legendInvMassPCMEMCAL->AddEntry(histoPCMEMCALSignalPlusBGPi0,"raw real events","p");
//      legendInvMassPCMEMCAL->AddEntry(histoPCMEMCALTotalBGPi0,"estimated total bkg","p");
//      legendInvMassPCMEMCAL->AddEntry(histoPCMEMCALBGPi0,"mixed event bkg","p");
//      legendInvMassPCMEMCAL->AddEntry(histoPCMEMCALRemainingBGPi0,"correlated bkg","p");
//      legendInvMassPCMEMCAL->Draw();
//      
//      
//  padInvMassSampleSignal->cd();
//      Double_t textsizeLabelsInvMassDown  = 0;
//      Double_t textsizeFacInvMassDown     = 0;
//      if (padInvMassSampleSignal->XtoPixel(padInvMassSampleSignal->GetX2()) < padInvMassSampleSignal->YtoPixel(padInvMassSampleSignal->GetY1())){
//          textsizeLabelsInvMassDown       = (Double_t)textSizeLabelsPixel/padInvMassSampleSignal->XtoPixel(padInvMassSampleSignal->GetX2()) ;
//          textsizeFacInvMassDown          = (Double_t)1./padInvMassSampleSignal->XtoPixel(padInvMassSampleSignal->GetX2()) ;
//      } else {
//          textsizeLabelsInvMassDown       = (Double_t)textSizeLabelsPixel/padInvMassSampleSignal->YtoPixel(padInvMassSampleSignal->GetY1());
//          textsizeFacInvMassDown          = (Double_t)1./padInvMassSampleSignal->YtoPixel(padInvMassSampleSignal->GetY1());
//      }
//      cout << textsizeLabelsInvMassDown << endl;
//      
//      TH2F * histo2DInvMassDummyDown;
//      histo2DInvMassDummyDown             = new TH2F("histo2DInvMassDummyDown","histo2DInvMassDummyDown",11000,0.05,0.225,20100,-100,20000);
//      SetStyleHistoTH2ForGraphs(histo2DInvMassDummyDown, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","counts",0.85*textsizeLabelsInvMassDown, textsizeLabelsInvMassDown,
//                                0.85*textsizeLabelsInvMassDown, textsizeLabelsInvMassDown,0.84,0.14/(textsizeFacInvMassDown*marginInvMass));
// //      histo2DInvMassDummyDown->GetXaxis()->SetLabelOffset(-0.01);
//      histo2DInvMassDummyDown->GetYaxis()->SetRangeUser(-20.,345);
//      histo2DInvMassDummyDown->GetXaxis()->SetTickLength(0.04);
//      histo2DInvMassDummyDown->GetYaxis()->SetTickLength(0.02);
//      histo2DInvMassDummyDown->DrawCopy();
// 
//      cout << 0.1/(textsizeFacInvMassDown*marginInvMass) << endl;
//      
//      TH1D* histoPCMEMCALSignalPlusBGPi0Clone     = (TH1D*)histoPCMEMCALSignalPlusBGPi0->Clone("histoPCMEMCALSignalPlusBGPi0Plot");
//      DrawGammaSetMarker(histoPCMEMCALSignalPlusBGPi0Clone, markerStyleInvMassSGBG2, markerSizeInvMassSGBG, markerColorInvMassSGBG2, markerColorInvMassSGBG2);
//      histoPCMEMCALSignalPlusBGPi0Clone->Draw("same");
//      DrawGammaSetMarker(histoPCMEMCALSignalPi0, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
//      histoPCMEMCALSignalPi0->Draw("same");
//      fitPCMEMCALSignalPi0->SetRange(0,0.25);
//      fitPCMEMCALSignalPi0->SetLineColor(fitColorInvMassSG);
//      fitPCMEMCALSignalPi0->Draw("same");
// 
//      TLegend* legendInvMassSigPCMEMCAL   = GetAndSetLegend2(0.64, 0.8, 0.9, 0.98, 0.85*textSizeLabelsPixel);
//      legendInvMassSigPCMEMCAL->SetMargin(0.25);
//      legendInvMassSigPCMEMCAL->AddEntry(histoPCMEMCALSignalPlusBGPi0Clone,"raw real events","p");
//      legendInvMassSigPCMEMCAL->AddEntry(histoPCMEMCALSignalPi0,"bkg subtracted","p");
//      legendInvMassSigPCMEMCAL->AddEntry(fitPCMEMCALSignalPi0,"fit","l");
//      legendInvMassSigPCMEMCAL->Draw();
//      
//  canvasInvMassSamplePlot->Print(Form("%s/InvMassSampleBinPCMEMCAL.%s",outputDir.Data(),suffix.Data()));
// 
//  // **********************************************************************************************************************
//  // **************************Plot example invariant mass bin EMCAL low pt ***********************************************
//  // **********************************************************************************************************************
//  
//  canvasInvMassSamplePlot->cd();
//  padInvMassSampleWithBG->cd();
// 
//      histo2DInvMassDummyUp->GetXaxis()->SetRangeUser(0.05,0.225);
//      histo2DInvMassDummyUp->GetYaxis()->SetRangeUser(-20.,410);
//      histo2DInvMassDummyUp->DrawCopy();
//      
//      DrawGammaSetMarker(histoEMCALSignalPlusBGPi0, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
//      histoEMCALSignalPlusBGPi0->Draw("same");
//      
//      DrawGammaSetMarker(histoEMCALTotalBGPi0, markerStyleInvMassBG, markerSizeInvMassBG, markerColorInvMassBG, markerColorInvMassBG);
//      histoEMCALTotalBGPi0->Draw("same");
//      DrawGammaSetMarker(histoEMCALBGPi0, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
//      histoEMCALBGPi0->Draw("same");
//      DrawGammaSetMarker(histoEMCALRemainingBGPi0, markerStyleInvMassCBG, markerSizeInvMassCBG, markerColorInvMassCBG, markerColorInvMassCBG);
//      histoEMCALRemainingBGPi0->Draw("same");
//  
//      labelInvMassEnergy->Draw();
// 
//      TLatex *labelInvMassRecoEMC         = new TLatex(0.135,0.84,"#gamma's rec. EMCal");
//      SetStyleTLatex( labelInvMassRecoEMC, 0.85*textSizeLabelsPixel,4);
//      labelInvMassRecoEMC->SetTextFont(43);
//      labelInvMassRecoEMC->Draw();
//      
//      TLatex *labelInvMassPtRangeEMCAL    = new TLatex(0.135,0.78,Form("%3.1f GeV/#it{c} < #it{p}_{T} < %3.1f GeV/#it{c}",rangePtExMin[2],rangePtExMax[2]));
//      SetStyleTLatex( labelInvMassPtRangeEMCAL, 0.85*textSizeLabelsPixel,4);
//      labelInvMassPtRangeEMCAL->SetTextFont(43);
//      labelInvMassPtRangeEMCAL->Draw();
// 
//      TLegend* legendInvMassEMCAL         = GetAndSetLegend2(0.64, 0.74, 0.9, 0.98, 0.85*textSizeLabelsPixel);
//      legendInvMassEMCAL->SetMargin(0.25);
//      legendInvMassEMCAL->AddEntry(histoEMCALSignalPlusBGPi0,"raw real events","p");
//      legendInvMassEMCAL->AddEntry(histoEMCALTotalBGPi0,"estimated total bkg","p");
//      legendInvMassEMCAL->AddEntry(histoEMCALBGPi0,"mixed event bkg","p");
//      legendInvMassEMCAL->AddEntry(histoEMCALRemainingBGPi0,"correlated bkg","p");
//      legendInvMassEMCAL->Draw();
//      
//  padInvMassSampleSignal->cd();
//      
//      histo2DInvMassDummyDown->GetXaxis()->SetRangeUser(0.05,0.245);
//      histo2DInvMassDummyDown->GetYaxis()->SetRangeUser(-20.,410);
//      histo2DInvMassDummyDown->DrawCopy();
//      
//      TH1D* histoEMCALSignalPlusBGPi0Clone    = (TH1D*)histoEMCALSignalPlusBGPi0->Clone("histoEMCALSignalPlusBGPi0Clone");
//      DrawGammaSetMarker(histoEMCALSignalPlusBGPi0Clone, markerStyleInvMassSGBG2, markerSizeInvMassSGBG, markerColorInvMassSGBG2, markerColorInvMassSGBG2);
//      histoEMCALSignalPlusBGPi0Clone->Draw("same");
//      DrawGammaSetMarker(histoEMCALSignalPi0, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
//      histoEMCALSignalPi0->Draw("same");
//      fitEMCALSignalPi0->SetRange(0,0.25);
//      fitEMCALSignalPi0->SetLineColor(fitColorInvMassSG);
//      fitEMCALSignalPi0->Draw("same");
// 
//      TLegend* legendInvMassSigEMCAL      = GetAndSetLegend2(0.64, 0.8, 0.9, 0.98, 0.85*textSizeLabelsPixel);
//      legendInvMassSigEMCAL->SetMargin(0.25);
//      legendInvMassSigEMCAL->AddEntry(histoEMCALSignalPlusBGPi0Clone,"raw real events","p");
//      legendInvMassSigEMCAL->AddEntry(histoEMCALSignalPi0,"bkg subtracted","p");
//      legendInvMassSigEMCAL->AddEntry(fitEMCALSignalPi0,"fit","l");
//      legendInvMassSigEMCAL->Draw();
// 
//      
//  canvasInvMassSamplePlot->Print(Form("%s/InvMassSampleBinEMCALLow.%s",outputDir.Data(),suffix.Data()));
// 
//  // **********************************************************************************************************************
//  // **************************Plot example invariant mass bin EMCAL high pt ***********************************************
//  // **********************************************************************************************************************
//  
//  canvasInvMassSamplePlot->cd();
//  padInvMassSampleWithBG->cd();
// 
//      histo2DInvMassDummyUp->GetXaxis()->SetRangeUser(0.05,0.225);
//      histo2DInvMassDummyUp->GetYaxis()->SetRangeUser(-20.,410);
//      histo2DInvMassDummyUp->DrawCopy();
//      
//      DrawGammaSetMarker(histoEMCALHighSignalPlusBGPi0, markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
//      histoEMCALHighSignalPlusBGPi0->Draw("same");
//      
//      DrawGammaSetMarker(histoEMCALHighTotalBGPi0, markerStyleInvMassBG, markerSizeInvMassBG, markerColorInvMassBG, markerColorInvMassBG);
//      histoEMCALHighTotalBGPi0->Draw("same");
//      DrawGammaSetMarker(histoEMCALHighBGPi0, markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
//      histoEMCALHighBGPi0->Draw("same");
//      DrawGammaSetMarker(histoEMCALHighRemainingBGPi0, markerStyleInvMassCBG, markerSizeInvMassCBG, markerColorInvMassCBG, markerColorInvMassCBG);
//      histoEMCALHighRemainingBGPi0->Draw("same");
//  
//      labelInvMassEnergy->Draw();
//      labelInvMassRecoEMC->Draw();
//      
//      TLatex *labelInvMassPtRangeEMCALH   = new TLatex(0.135,0.78,Form("%3.1f GeV/#it{c} < #it{p}_{T} < %3.1f GeV/#it{c}",rangePtExMin[8],rangePtExMax[8]));
//      SetStyleTLatex( labelInvMassPtRangeEMCALH, 0.85*textSizeLabelsPixel,4);
//      labelInvMassPtRangeEMCALH->SetTextFont(43);
//      labelInvMassPtRangeEMCALH->Draw();
// 
//      TLegend* legendInvMassEMCALH        = GetAndSetLegend2(0.64, 0.74, 0.9, 0.98, 0.85*textSizeLabelsPixel);
//      legendInvMassEMCALH->SetMargin(0.25);
//      legendInvMassEMCALH->AddEntry(histoEMCALHighSignalPlusBGPi0,"raw real events","p");
//      legendInvMassEMCALH->AddEntry(histoEMCALHighTotalBGPi0,"estimated total bkg","p");
//      legendInvMassEMCALH->AddEntry(histoEMCALHighBGPi0,"mixed event bkg","p");
//      legendInvMassEMCALH->AddEntry(histoEMCALHighRemainingBGPi0,"correlated bkg","p");
//      legendInvMassEMCALH->Draw();
//      
//  padInvMassSampleSignal->cd();
//      
//      histo2DInvMassDummyDown->GetXaxis()->SetRangeUser(0.05,0.245);
//      histo2DInvMassDummyDown->GetYaxis()->SetRangeUser(-20.,410);
//      histo2DInvMassDummyDown->DrawCopy();
//      
//      TH1D* histoEMCALHighSignalPlusBGPi0Clone    = (TH1D*)histoEMCALHighSignalPlusBGPi0->Clone("histoEMCALHighSignalPlusBGPi0Clone");
//      DrawGammaSetMarker(histoEMCALHighSignalPlusBGPi0Clone, markerStyleInvMassSGBG2, markerSizeInvMassSGBG, markerColorInvMassSGBG2, markerColorInvMassSGBG2);
//      histoEMCALHighSignalPlusBGPi0Clone->Draw("same");
//      DrawGammaSetMarker(histoEMCALHighSignalPi0, markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
//      histoEMCALHighSignalPi0->Draw("same");
// //      histoEMCALHighSignalPi0Templ->SetRange(0,0.25);
//      histoEMCALHighSignalPi0Templ->SetLineColor(fitColorInvMassSG);
//      histoEMCALHighSignalPi0Templ->Draw("same,hist");
// 
//      TLegend* legendInvMassSigEMCALH             = GetAndSetLegend2(0.64, 0.8, 0.9, 0.98, 0.85*textSizeLabelsPixel);
//      legendInvMassSigEMCALH->SetMargin(0.25);
//      legendInvMassSigEMCALH->AddEntry(histoEMCALHighSignalPlusBGPi0Clone,"raw real events","p");
//      legendInvMassSigEMCALH->AddEntry(histoEMCALHighSignalPi0,"bkg subtracted","p");
//      legendInvMassSigEMCALH->AddEntry(histoEMCALHighSignalPi0Templ,"template fit","l");
//      legendInvMassSigEMCALH->Draw();
// 
//      
//  canvasInvMassSamplePlot->Print(Form("%s/InvMassSampleBinEMCALHigh.%s",outputDir.Data(),suffix.Data()));
    
    
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
        // PCM-EMCal component
        graphPCMEMCALPi0InvXSectionStat->Write("graphInvCrossSectionPi0PCMEMCAL2760GeVStatErr");
        graphPCMEMCALPi0InvXSectionSys->Write("graphInvCrossSectionPi0PCMEMCAL2760GeVSysErr");
        // EMCal merged component
        graphEMCALMergedPi0InvXSectionStat->Write("graphInvCrossSectionPi0EMCALMerged2760GeVStatErr");
        graphEMCALMergedPi0InvXSectionSys->Write("graphInvCrossSectionPi0EMCALMerged2760GeVSysErr");
        // Final spectrum correlations Method A
        graphCombPi0InvXSectionTotA->Write("graphInvCrossSectionPi0Comb2760GeVA");
        graphCombPi0InvXSectionStatA->Write("graphInvCrossSectionPi0Comb2760GeVAStatErr");
        graphCombPi0InvXSectionSysA->Write("graphInvCrossSectionPi0Comb2760GeVASysErr");  
    
    fCombResults.mkdir("Eta2.76TeV");
    TDirectoryFile* directoryEta = (TDirectoryFile*)fCombResults.Get("Eta2.76TeV"); 
    fCombResults.cd("Eta2.76TeV");
        // PCM component
        graphPCMEtaInvXSectionStat->Write("graphInvCrossSectionEtaPCM2760GeVStatErr");
        graphPCMEtaInvXSectionSys->Write("graphInvCrossSectionEtaPCM2760GeVSysErr");
        // EMCAL component
        graphEMCALEtaInvXSectionStat->Write("graphInvCrossSectionEtaEMCAL2760GeVStatErr");
        graphEMCALEtaInvXSectionSys->Write("graphInvCrossSectionEtaEMCAL2760GeVSysErr");
        // PCM-EMCal component
        graphPCMEMCALEtaInvXSectionStat->Write("graphInvCrossSectionEtaPCMEMCAL2760GeVStatErr");
        graphPCMEMCALEtaInvXSectionSys->Write("graphInvCrossSectionEtaPCMEMCAL2760GeVSysErr");
        // Final spectrum correlations Method A
        graphCombEtaInvXSectionTotA->Write("graphInvCrossSectionEtaComb2760GeVA");
        graphCombEtaInvXSectionStatA->Write("graphInvCrossSectionEtaComb2760GeVAStatErr");
        graphCombEtaInvXSectionSysA->Write("graphInvCrossSectionEtaComb2760GeVASysErr");  
    
    fCombResults.Close();

 

}
    
