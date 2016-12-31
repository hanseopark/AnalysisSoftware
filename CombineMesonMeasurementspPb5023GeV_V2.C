/****************************************************************************************************************************
******        provided by Gamma Conversion Group, PWGGA,                                                     *****
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

void CombineMesonMeasurementspPb5023GeV_V2(     TString fileNamePCM         = "", 
                                                Bool_t flagPCMfile          = 0,
                                                TString fileNameEMCAL       = "",  
                                                Bool_t flagEMCfile          = 0,
                                                TString fileNamePHOS        = "",
                                                Bool_t flagPHOSfile         = 0,
                                                TString fileNameDalitz      = "",
                                                Bool_t flagDalitzfile       = 0,
                                                TString fileNamePCMEMCAL    = "",
                                                TString suffix              = "eps", 
                                                TString isMC                = "", 
                                                TString bWCorrection        = "X",
                                                TString fileNameCorrFactors = "",
                                                Bool_t isNSD                = kTRUE
                                                
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
    TString collisionSystempPb                  = "p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";    
    TString outputDir                           = Form("%s/%s/CombineMesonMeasurementspPb5023GeV%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
    
    TString nameFinalResDat                     = Form("%s/CombinedResults%s_FitResults.dat",outputDir.Data(),bWCorrection.Data());
    fstream  fileFitsOutput;
    fileFitsOutput.open(nameFinalResDat.Data(), ios::out);  

    TString fileNameChargedPion                 = "ExternalInputpPb/IdentifiedCharged/ChargedIdentifiedSpectrapPb_2016_12_22.root";
    TString fileNameEtaToPi0WorldData           = "ExternalInput/WorldDataPi0Eta.root";
    
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMEMCAL.root", fileNamePCMEMCAL.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCAL.root", fileNameEMCAL.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputDalitz.root", fileNameDalitz.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputIdentifiedCharged.root", fileNameChargedPion.Data(), outputDir.Data()));    
    
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
                                                    "PCM-Dal", "PHOS-Dal", "EMC-Dal", "EMC high pT", "mEMC",
                                                    "PCMOtherDataset"};
    TString  nameTrigger[6]                     = {"INT1", "INT7", "EMC1", "EMC7", "EG2", "EG1"};
    
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
    Double_t scalingToNSD                       = 0.964;
    if (!isNSD)                                                    
        scalingToNSD                            = 1;
    
    for (Int_t i = 0; i < 11; i++){
        colorDet[i]                             = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kFALSE, kFALSE, kTRUE);
        colorDetMC[i]                           = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kTRUE, kFALSE, kTRUE);
        markerStyleDet[i]                       = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
        markerStyleDetMC[i]                     = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kTRUE);
        markerSizeDet[i]                        = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kFALSE)*2;
        markerSizeDetMC[i]                      = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kTRUE)*2;
    }    
    

    //************************** Read data for PCM **************************************************
    TFile* filePCM                                  = NULL;
    TDirectory* directoryPCMPi0                     = NULL;
    TGraphAsymmErrors* graphPCMPi0Mass              = NULL;
    TGraphAsymmErrors* graphPCMPi0FWHM              = NULL;
    TGraphAsymmErrors* graphPCMPi0MassMC            = NULL;
    TGraphAsymmErrors* graphPCMPi0FWHMMC            = NULL;
    TGraphAsymmErrors* graphPCMPi0cc                = NULL;
    TGraphAsymmErrors* graphPCMPi0EffPt             = NULL;
    TGraphAsymmErrors* graphPCMPi0InvYieldStat      = NULL;
    TH1D* histoPCMPi0InvYieldStat                   = NULL;
    TGraphAsymmErrors* graphPCMPi0InvYieldSys       = NULL;
    TGraphAsymmErrors* graphPCMPi0AccTimesEff        = NULL;
    TDirectory* directoryPCMEta                     = NULL;
    TGraphAsymmErrors* graphPCMEtaMass              = NULL;
    TGraphAsymmErrors* graphPCMEtaFWHM              = NULL;
    TGraphAsymmErrors* graphPCMEtaMassMC            = NULL;
    TGraphAsymmErrors* graphPCMEtaFWHMMC            = NULL;
    TGraphAsymmErrors* graphPCMEtaAcc               = NULL;
    TGraphAsymmErrors* graphPCMEtaEffPt             = NULL;
    TGraphAsymmErrors* graphPCMEtaInvYieldStat      = NULL;
    TH1D* histoPCMEtaInvYieldStat                   = NULL;
    TGraphAsymmErrors* graphPCMEtaInvYieldSys       = NULL;
    TGraphAsymmErrors* graphPCMEtaAccTimesEff       = NULL;
    TH1D* histoPCMEtaToPi0Stat                      = NULL;
    TGraphAsymmErrors* graphPCMEtaToPi0Sys          = NULL;

    if (!flagPCMfile){  
        filePCM                                     = new TFile(fileNamePCM.Data());
        directoryPCMPi0                             = (TDirectory*)filePCM->Get("Pi0_pPb_5.023TeV_0-100%"); 
            TH1D* histoPCMPi0Mass                       = (TH1D*)directoryPCMPi0->Get("MassPi0");
            TH1D* histoPCMPi0FWHMMeV                    = (TH1D*)directoryPCMPi0->Get("FWHMPi0MeV");
            TH1D* histoPCMPi0TrueMass                   = (TH1D*)directoryPCMPi0->Get("TrueMassPi0");
            TH1D* histoPCMPi0TrueFWHMMeV                = (TH1D*)directoryPCMPi0->Get("TrueFWHMPi0MeV");
            TH1D* histoPCMPi0cc                         = (TH1D*)directoryPCMPi0->Get("AcceptancePi0");
            TH1D* histoPCMPi0TrueEffPt                  = (TH1D*)directoryPCMPi0->Get("EfficiencyPi0");
            histoPCMPi0InvYieldStat                     = (TH1D*)directoryPCMPi0->Get("CorrectedYieldPi0");
            graphPCMPi0InvYieldStat                     = new TGraphAsymmErrors(histoPCMPi0InvYieldStat);
            while (graphPCMPi0InvYieldStat->GetY()[0] == 0) graphPCMPi0InvYieldStat->RemovePoint(0);
            graphPCMPi0InvYieldSys                      = (TGraphAsymmErrors*)directoryPCMPi0->Get("Pi0SystError");
            while (graphPCMPi0InvYieldSys->GetY()[0] == 0) graphPCMPi0InvYieldSys->RemovePoint(0);
            if (isNSD) {
                histoPCMPi0InvYieldStat->Scale(scalingToNSD);
                graphPCMPi0InvYieldStat                 = ScaleGraph(graphPCMPi0InvYieldStat,scalingToNSD);
                graphPCMPi0InvYieldSys                  = ScaleGraph(graphPCMPi0InvYieldSys,scalingToNSD);
            }
            TH1D* histoPCMPi0AccTimesEff                 = (TH1D*)histoPCMPi0TrueEffPt->Clone("histoPCMPi0AccTimesEff");
            histoPCMPi0AccTimesEff->Multiply(histoPCMPi0cc);
            // normalize to full acceptance (delta y and phi)
            histoPCMPi0AccTimesEff->Scale(2*TMath::Pi()*1.6);
            histoPCMPi0Mass->Scale(1000.);
            histoPCMPi0TrueMass->Scale(1000.);

            graphPCMPi0Mass                             = new TGraphAsymmErrors(histoPCMPi0Mass);
            graphPCMPi0FWHM                             = new TGraphAsymmErrors(histoPCMPi0FWHMMeV);
            graphPCMPi0MassMC                           = new TGraphAsymmErrors(histoPCMPi0TrueMass);
            graphPCMPi0FWHMMC                           = new TGraphAsymmErrors(histoPCMPi0TrueFWHMMeV);
            graphPCMPi0cc                               = new TGraphAsymmErrors(histoPCMPi0cc);
            graphPCMPi0EffPt                            = new TGraphAsymmErrors(histoPCMPi0TrueEffPt);
            graphPCMPi0AccTimesEff                       = new TGraphAsymmErrors(histoPCMPi0AccTimesEff);
            
            while(graphPCMPi0AccTimesEff->GetY()[0] == 0){
                graphPCMPi0Mass->RemovePoint(0);
                graphPCMPi0FWHM->RemovePoint(0);
                graphPCMPi0MassMC->RemovePoint(0);
                graphPCMPi0FWHMMC->RemovePoint(0);
                graphPCMPi0cc->RemovePoint(0);
                graphPCMPi0EffPt->RemovePoint(0);
                graphPCMPi0AccTimesEff->RemovePoint(0);
            }
                
        directoryPCMEta                             = (TDirectory*)filePCM->Get("Eta_pPb_5.023TeV_0-100%"); 
            TH1D* histoPCMEtaMass                       = (TH1D*)directoryPCMEta->Get("MassEta");
            TH1D* histoPCMEtaFWHMMeV                    = (TH1D*)directoryPCMEta->Get("FWHMEtaMeV");
            TH1D* histoPCMEtaTrueMass                   = (TH1D*)directoryPCMEta->Get("TrueMassEta");
            TH1D* histoPCMEtaTrueFWHMMeV                = (TH1D*)directoryPCMEta->Get("TrueFWHMEtaMeV");
            TH1D* histoPCMEtaAcc                        = (TH1D*)directoryPCMEta->Get("AcceptanceEta");
            TH1D* histoPCMEtaTrueEffPt                  = (TH1D*)directoryPCMEta->Get("EfficiencyEta");
            histoPCMEtaInvYieldStat                     = (TH1D*)directoryPCMEta->Get("CorrectedYieldEta");
            graphPCMEtaInvYieldStat                     = new TGraphAsymmErrors(histoPCMEtaInvYieldStat);
            while (graphPCMEtaInvYieldStat->GetY()[0] == 0)graphPCMEtaInvYieldStat->RemovePoint(0);
            graphPCMEtaInvYieldSys                      = (TGraphAsymmErrors*)directoryPCMEta->Get("EtaSystError");
            while (graphPCMEtaInvYieldSys->GetY()[0] == 0)graphPCMEtaInvYieldSys->RemovePoint(0);
            if (isNSD) {
                histoPCMEtaInvYieldStat->Scale(scalingToNSD);
                graphPCMEtaInvYieldStat                 = ScaleGraph(graphPCMEtaInvYieldStat,scalingToNSD);
                graphPCMEtaInvYieldSys                  = ScaleGraph(graphPCMEtaInvYieldSys,scalingToNSD);
            }

            histoPCMEtaToPi0Stat                        = (TH1D*)directoryPCMEta->Get("EtatoPi0Ratio");
            graphPCMEtaToPi0Sys                         = (TGraphAsymmErrors*)directoryPCMEta->Get("EtatoPi0RatioSys");
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
    } else {
        filePCM                                     = new TFile(fileNamePCM.Data());
        directoryPCMPi0                             = (TDirectory*)filePCM->Get("Pi0pPb_5.023TeV"); 
            graphPCMPi0Mass                             = (TGraphAsymmErrors*)directoryPCMPi0->Get("Pi0_Mass_data");
            graphPCMPi0Mass                             = ScaleGraph(graphPCMPi0Mass, 1000.);
            graphPCMPi0FWHM                             = (TGraphAsymmErrors*)directoryPCMPi0->Get("Pi0_Width_data");
            graphPCMPi0FWHM                             = ScaleGraph(graphPCMPi0FWHM, 1000.);
            graphPCMPi0MassMC                           = (TGraphAsymmErrors*)directoryPCMPi0->Get("Pi0_Mass_MC");
            graphPCMPi0MassMC                           = ScaleGraph(graphPCMPi0MassMC, 1000.);
            graphPCMPi0FWHMMC                           = (TGraphAsymmErrors*)directoryPCMPi0->Get("Pi0_Width_MC");
            graphPCMPi0FWHMMC                           = ScaleGraph(graphPCMPi0FWHMMC, 1000.);
            graphPCMPi0cc                               = (TGraphAsymmErrors*)directoryPCMPi0->Get("AcceptancePi0");
            graphPCMPi0EffPt                            = (TGraphAsymmErrors*)directoryPCMPi0->Get("EfficiencyPi0");
            graphPCMPi0InvYieldStat                     = (TGraphAsymmErrors*)directoryPCMPi0->Get("graphCorrectedYieldPi0");
            histoPCMPi0InvYieldStat                     = (TH1D*)directoryPCMPi0->Get("CorrectedYieldPi0");
            cout << "Pi0 stat PCM" << endl;
            graphPCMPi0InvYieldStat->Print();
            graphPCMPi0InvYieldSys                      = (TGraphAsymmErrors*)directoryPCMPi0->Get("Pi0SystError");
            cout << "Pi0 sys PCM" << endl;
            graphPCMPi0InvYieldSys->Print();
            if (isNSD) {
                histoPCMPi0InvYieldStat->Scale(scalingToNSD);
                graphPCMPi0InvYieldStat                 = ScaleGraph(graphPCMPi0InvYieldStat,scalingToNSD);
                graphPCMPi0InvYieldSys                  = ScaleGraph(graphPCMPi0InvYieldSys,scalingToNSD);
            }
            graphPCMPi0AccTimesEff                      = (TGraphAsymmErrors*)directoryPCMPi0->Get("EffTimesAccPi0");

        directoryPCMEta                                 = (TDirectory*)filePCM->Get("EtapPb_5.023TeV"); 
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
            histoPCMEtaInvYieldStat                     = (TH1D*)directoryPCMEta->Get("CorrectedYieldEta");
            
            graphPCMEtaInvYieldStat                     = (TGraphAsymmErrors*)directoryPCMEta->Get("graphCorrectedYieldEta");
            cout << "Eta stat PCM-EMC" << endl;
            graphPCMEtaInvYieldStat->Print();
            graphPCMEtaInvYieldSys                      = (TGraphAsymmErrors*)directoryPCMEta->Get("EtaSystError");
            cout << "Eta sys PCM-EMC" << endl;
            graphPCMEtaInvYieldSys->Print();        
            if (isNSD) {
                histoPCMEtaInvYieldStat->Scale(scalingToNSD);
                graphPCMEtaInvYieldStat                 = ScaleGraph(graphPCMEtaInvYieldStat,scalingToNSD);
                graphPCMEtaInvYieldSys                  = ScaleGraph(graphPCMEtaInvYieldSys,scalingToNSD);
            }
            histoPCMEtaToPi0Stat                        = (TH1D*)directoryPCMEta->Get("EtaToPi0YShiftedStatError");
            graphPCMEtaToPi0Sys                         = (TGraphAsymmErrors*)directoryPCMEta->Get("EtaToPi0YShiftedSystError");            
            if (!histoPCMEtaToPi0Stat){
                histoPCMEtaToPi0Stat                        = (TH1D*)directoryPCMEta->Get("EtaToPi0StatError");
            } 
            if (!graphPCMEtaToPi0Sys){
                graphPCMEtaToPi0Sys                         = (TGraphAsymmErrors*)directoryPCMEta->Get("EtaToPi0SystError");
            }    
    }
    graphPCMPi0AccTimesEff                              = ScaleGraph(graphPCMPi0AccTimesEff,0.98823);
    
    
    //************************** Read data for PCMEMCAL **************************************************
    TFile* filePCMEMCAL                                     = new TFile(fileNamePCMEMCAL.Data());
    TDirectory* directoryPCMEMCALPi0                        = (TDirectory*)filePCMEMCAL->Get("Pi0pPb_5.023TeV"); 
        TGraphAsymmErrors* graphPCMEMCALPi0Mass             = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("Pi0_Mass_data");
        graphPCMEMCALPi0Mass                                = ScaleGraph(graphPCMEMCALPi0Mass, 1000.);
        TGraphAsymmErrors* graphPCMEMCALPi0FWHM             = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("Pi0_Width_data");
        graphPCMEMCALPi0FWHM                                = ScaleGraph(graphPCMEMCALPi0FWHM, 1000.);
        TGraphAsymmErrors* graphPCMEMCALPi0MassMC           = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("Pi0_Mass_MC");
        graphPCMEMCALPi0MassMC                              = ScaleGraph(graphPCMEMCALPi0MassMC, 1000.);
        TGraphAsymmErrors* graphPCMEMCALPi0FWHMMC           = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("Pi0_Width_MC");
        graphPCMEMCALPi0FWHMMC                              = ScaleGraph(graphPCMEMCALPi0FWHMMC, 1000.);
        TGraphAsymmErrors* graphPCMEMCALPi0cc               = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("AcceptancePi0");
        TGraphAsymmErrors* graphPCMEMCALPi0EffPt            = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("EfficiencyPi0");
        TGraphAsymmErrors* graphPCMEMCALPi0InvYieldStat     = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("graphCorrectedYieldPi0");
        TH1D* histoPCMEMCALPi0InvYieldStat                  = (TH1D*)directoryPCMEMCALPi0->Get("CorrectedYieldPi0");
        cout << "Pi0 stat PCM-EMC" << endl;
        graphPCMEMCALPi0InvYieldStat->Print();
        TGraphAsymmErrors* graphPCMEMCALPi0InvYieldSys      = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("Pi0SystError");
        cout << "Pi0 sys PCM-EMC" << endl;
        graphPCMEMCALPi0InvYieldSys->Print();
        if (isNSD) {
            histoPCMEMCALPi0InvYieldStat->Scale(scalingToNSD);
            graphPCMEMCALPi0InvYieldStat                    = ScaleGraph(graphPCMEMCALPi0InvYieldStat,scalingToNSD);
            graphPCMEMCALPi0InvYieldSys                     = ScaleGraph(graphPCMEMCALPi0InvYieldSys,scalingToNSD);
        }
        TGraphAsymmErrors* graphPCMEMCALPi0AccTimesEff      = (TGraphAsymmErrors*)directoryPCMEMCALPi0->Get("EffTimesAccPi0");
        graphPCMEMCALPi0AccTimesEff                         = ScaleGraph(graphPCMEMCALPi0AccTimesEff,0.98823);
        
    TDirectory* directoryPCMEMCALEta                        = (TDirectory*)filePCMEMCAL->Get("EtapPb_5.023TeV"); 
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
        TH1D* histoPCMEMCALEtaInvYieldStat                  = (TH1D*)directoryPCMEMCALEta->Get("CorrectedYieldEta");
        TGraphAsymmErrors* graphPCMEMCALEtaInvYieldStat     = (TGraphAsymmErrors*)directoryPCMEMCALEta->Get("graphCorrectedYieldEta");
        cout << "Eta stat PCM-EMC" << endl;
        graphPCMEMCALEtaInvYieldStat->Print();
        TGraphAsymmErrors* graphPCMEMCALEtaInvYieldSys   = (TGraphAsymmErrors*)directoryPCMEMCALEta->Get("EtaSystError");
        cout << "Eta sys PCM-EMC" << endl;
        graphPCMEMCALEtaInvYieldSys->Print();        
        if (isNSD) {
            histoPCMEMCALEtaInvYieldStat->Scale(scalingToNSD);
            graphPCMEMCALEtaInvYieldStat                    = ScaleGraph(graphPCMEMCALEtaInvYieldStat,scalingToNSD);
            graphPCMEMCALEtaInvYieldSys                     = ScaleGraph(graphPCMEMCALEtaInvYieldSys,scalingToNSD);
        }
        TH1D* histoPCMEMCALEtaToPi0Stat                     = (TH1D*)directoryPCMEMCALEta->Get("EtaToPi0YShiftedStatError");
        TGraphAsymmErrors* graphPCMEMCALEtaToPi0Sys         = (TGraphAsymmErrors*)directoryPCMEMCALEta->Get("EtaToPi0YShiftedSystError");
        if (!histoPCMEMCALEtaToPi0Stat){
            histoPCMEMCALEtaToPi0Stat                       = (TH1D*)directoryPCMEMCALEta->Get("EtaToPi0StatError");
        } 
        if (!graphPCMEMCALEtaToPi0Sys){
            graphPCMEMCALEtaToPi0Sys                        = (TGraphAsymmErrors*)directoryPCMEMCALEta->Get("EtaToPi0SystError");
        } 
        
        
    //************************** Read data for EMCAL ****************************************************
    TFile* fileEMCAL                                    = NULL;
    TDirectory* directoryEMCALPi0                       = NULL;
    TGraphAsymmErrors* graphEMCALPi0Mass                = NULL;
    TGraphAsymmErrors* graphEMCALPi0FWHM                = NULL;
    TGraphAsymmErrors* graphEMCALPi0MassMC              = NULL;
    TGraphAsymmErrors* graphEMCALPi0FWHMMC              = NULL;
    TGraphAsymmErrors* graphEMCALPi0cc                  = NULL;
    TGraphAsymmErrors* graphEMCALPi0EffPt               = NULL;
    TGraphAsymmErrors* graphEMCALPi0AccTimesEff          = NULL;
    TH1D* histoEMCALPi0InvYieldStat                     = NULL;
    TGraphAsymmErrors* graphEMCALPi0InvYieldStat        = NULL;
    TGraphAsymmErrors* graphEMCALPi0InvYieldSys         = NULL;
    TDirectory* directoryEMCALEta                       = NULL;
    TGraphAsymmErrors* graphEMCALEtaMass                = NULL;
    TGraphAsymmErrors* graphEMCALEtaFWHM                = NULL;
    TGraphAsymmErrors* graphEMCALEtaMassMC              = NULL;
    TGraphAsymmErrors* graphEMCALEtaFWHMMC              = NULL;
    TGraphAsymmErrors* graphEMCALEtaAcc                 = NULL;
    TGraphAsymmErrors* graphEMCALEtaEffPt               = NULL;
    TGraphAsymmErrors* graphEMCALEtaAccTimesEff         = NULL;
    TH1D* histoEMCALEtaInvYieldStat                     = NULL;
    TGraphAsymmErrors* graphEMCALEtaInvYieldStat        = NULL;
    TGraphAsymmErrors* graphEMCALEtaInvYieldSys         = NULL;
    TH1D* histoEMCALEtaToPi0Stat                        = NULL;
    TGraphAsymmErrors* graphEMCALEtaToPi0Sys            = NULL;

    if (!flagEMCfile){
        fileEMCAL                                   = new TFile(fileNameEMCAL.Data());
        directoryEMCALPi0                           = (TDirectory*)fileEMCAL->Get("Pi0_pPb_5.023TeV_0-100%"); 
            TH1D* histoEMCALPi0Mass                     = (TH1D*)directoryEMCALPi0->Get("MassPi0");
            TH1D* histoEMCALPi0FWHMMeV                  = (TH1D*)directoryEMCALPi0->Get("FWHMPi0MeV");
            TH1D* histoEMCALPi0TrueMass                 = (TH1D*)directoryEMCALPi0->Get("TrueMassPi0");
            TH1D* histoEMCALPi0TrueFWHMMeV              = (TH1D*)directoryEMCALPi0->Get("TrueFWHMPi0MeV");
            TH1D* histoEMCALPi0cc                       = (TH1D*)directoryEMCALPi0->Get("AcceptancePi0");
            TH1D* histoEMCALPi0TrueEffPt                = (TH1D*)directoryEMCALPi0->Get("EfficiencyPi0");
            histoEMCALPi0InvYieldStat                   = (TH1D*)directoryEMCALPi0->Get("CorrectedYieldPi0");
            graphEMCALPi0InvYieldStat                   = new TGraphAsymmErrors(histoEMCALPi0InvYieldStat);
            while (graphEMCALPi0InvYieldStat->GetY()[0] == 0) graphEMCALPi0InvYieldStat->RemovePoint(0);
            graphEMCALPi0InvYieldSys                    = (TGraphAsymmErrors*)directoryEMCALPi0->Get("Pi0SystError");
            while (graphEMCALPi0InvYieldSys->GetY()[0] == 0) graphEMCALPi0InvYieldSys->RemovePoint(0);
            if (isNSD) {
                histoEMCALPi0InvYieldStat->Scale(scalingToNSD);
                graphEMCALPi0InvYieldStat               = ScaleGraph(graphEMCALPi0InvYieldStat,scalingToNSD);
                graphEMCALPi0InvYieldSys                = ScaleGraph(graphEMCALPi0InvYieldSys,scalingToNSD);
            }
            TH1D* histoEMCALPi0AccTimesEff              = (TH1D*)histoEMCALPi0TrueEffPt->Clone("histoEMCALPi0AccTimesEff");
            histoEMCALPi0AccTimesEff->Multiply(histoEMCALPi0cc);
            // normalize to full acceptance (delta y and phi)
            histoEMCALPi0AccTimesEff->Scale(2*TMath::Pi()*1.6);
            histoEMCALPi0Mass->Scale(1000.);
            histoEMCALPi0TrueMass->Scale(1000.);

            graphEMCALPi0Mass                           = new TGraphAsymmErrors(histoEMCALPi0Mass);
            graphEMCALPi0FWHM                           = new TGraphAsymmErrors(histoEMCALPi0FWHMMeV);
            graphEMCALPi0MassMC                         = new TGraphAsymmErrors(histoEMCALPi0TrueMass);
            graphEMCALPi0FWHMMC                         = new TGraphAsymmErrors(histoEMCALPi0TrueFWHMMeV);
            graphEMCALPi0cc                             = new TGraphAsymmErrors(histoEMCALPi0cc);
            graphEMCALPi0EffPt                          = new TGraphAsymmErrors(histoEMCALPi0TrueEffPt);
            graphEMCALPi0AccTimesEff                     = new TGraphAsymmErrors(histoEMCALPi0AccTimesEff);
            
            while(graphEMCALPi0AccTimesEff->GetY()[0] == 0){
                graphEMCALPi0Mass->RemovePoint(0);
                graphEMCALPi0FWHM->RemovePoint(0);
                graphEMCALPi0MassMC->RemovePoint(0);
                graphEMCALPi0FWHMMC->RemovePoint(0);
                graphEMCALPi0cc->RemovePoint(0);
                graphEMCALPi0EffPt->RemovePoint(0);
                graphEMCALPi0AccTimesEff->RemovePoint(0);
            }
                
        directoryEMCALEta                           = (TDirectory*)fileEMCAL->Get("Eta_pPb_5.023TeV_0-100%"); 
            TH1D* histoEMCALEtaMass                     = (TH1D*)directoryEMCALEta->Get("MassEta");
            TH1D* histoEMCALEtaFWHMMeV                  = (TH1D*)directoryEMCALEta->Get("FWHMEtaMeV");
            TH1D* histoEMCALEtaTrueMass                 = (TH1D*)directoryEMCALEta->Get("TrueMassEta");
            TH1D* histoEMCALEtaTrueFWHMMeV              = (TH1D*)directoryEMCALEta->Get("TrueFWHMEtaMeV");
            TH1D* histoEMCALEtaAcc                      = (TH1D*)directoryEMCALEta->Get("AcceptanceEta");
            TH1D* histoEMCALEtaTrueEffPt                = (TH1D*)directoryEMCALEta->Get("EfficiencyEta");
            histoEMCALEtaInvYieldStat                   = (TH1D*)directoryEMCALEta->Get("CorrectedYieldEta");
            graphEMCALEtaInvYieldStat                   = new TGraphAsymmErrors(histoEMCALEtaInvYieldStat);
            while (graphEMCALEtaInvYieldStat->GetY()[0] == 0)graphEMCALEtaInvYieldStat->RemovePoint(0);
            graphEMCALEtaInvYieldSys                    = (TGraphAsymmErrors*)directoryEMCALEta->Get("EtaSystError");
            while (graphEMCALEtaInvYieldSys->GetY()[0] == 0)graphEMCALEtaInvYieldSys->RemovePoint(0);
            histoEMCALEtaToPi0Stat                      = (TH1D*)directoryEMCALEta->Get("EtatoPi0Ratio");
            graphEMCALEtaToPi0Sys                       = (TGraphAsymmErrors*)directoryEMCALEta->Get("EtatoPi0RatioSys");
            if (isNSD) {
                histoEMCALEtaInvYieldStat->Scale(scalingToNSD);
                graphEMCALEtaInvYieldStat               = ScaleGraph(graphEMCALEtaInvYieldStat,scalingToNSD);
                graphEMCALEtaInvYieldSys                = ScaleGraph(graphEMCALEtaInvYieldSys,scalingToNSD);
            }
            TH1D* histoEMCALEtaAccTimesEff              = (TH1D*)histoEMCALEtaTrueEffPt->Clone("histoEMCALEtaAccTimesEff");
            histoEMCALEtaAccTimesEff->Multiply(histoEMCALEtaAcc);
            // normalize to full acceptance (delta y and phi)
            histoEMCALEtaAccTimesEff->Scale(2*TMath::Pi()*1.6);
            histoEMCALEtaMass->Scale(1000.);
            histoEMCALEtaTrueMass->Scale(1000.);
        
            graphEMCALEtaMass                           = new TGraphAsymmErrors(histoEMCALEtaMass);
            graphEMCALEtaFWHM                           = new TGraphAsymmErrors(histoEMCALEtaFWHMMeV);
            graphEMCALEtaMassMC                         = new TGraphAsymmErrors(histoEMCALEtaTrueMass);
            graphEMCALEtaFWHMMC                         = new TGraphAsymmErrors(histoEMCALEtaTrueFWHMMeV);
            graphEMCALEtaAcc                            = new TGraphAsymmErrors(histoEMCALEtaAcc);
            graphEMCALEtaEffPt                          = new TGraphAsymmErrors(histoEMCALEtaTrueEffPt);
            graphEMCALEtaAccTimesEff                    = new TGraphAsymmErrors(histoEMCALEtaAccTimesEff);

            while(graphEMCALEtaAccTimesEff->GetY()[0] == 0){
                graphEMCALEtaMass->RemovePoint(0);
                graphEMCALEtaFWHM->RemovePoint(0);
                graphEMCALEtaMassMC->RemovePoint(0);
                graphEMCALEtaFWHMMC->RemovePoint(0);
                graphEMCALEtaAcc->RemovePoint(0);
                graphEMCALEtaEffPt->RemovePoint(0);
                graphEMCALEtaAccTimesEff->RemovePoint(0);
            }
    } else {    
        fileEMCAL                                   = new TFile(fileNameEMCAL.Data());
        directoryEMCALPi0                           = (TDirectory*)fileEMCAL->Get("Pi0pPb_5.023TeV"); 
            graphEMCALPi0Mass                           = (TGraphAsymmErrors*)directoryEMCALPi0->Get("Pi0_Mass_data");
            graphEMCALPi0Mass                           = ScaleGraph(graphEMCALPi0Mass, 1000.);
            graphEMCALPi0FWHM                           = (TGraphAsymmErrors*)directoryEMCALPi0->Get("Pi0_Width_data");
            graphEMCALPi0FWHM                           = ScaleGraph(graphEMCALPi0FWHM, 1000.);
            graphEMCALPi0MassMC                         = (TGraphAsymmErrors*)directoryEMCALPi0->Get("Pi0_Mass_MC");
            graphEMCALPi0MassMC                         = ScaleGraph(graphEMCALPi0MassMC, 1000.);
            graphEMCALPi0FWHMMC                         = (TGraphAsymmErrors*)directoryEMCALPi0->Get("Pi0_Width_MC");
            graphEMCALPi0FWHMMC                         = ScaleGraph(graphEMCALPi0FWHMMC, 1000.);
            graphEMCALPi0cc                             = (TGraphAsymmErrors*)directoryEMCALPi0->Get("AcceptancePi0");
            graphEMCALPi0EffPt                          = (TGraphAsymmErrors*)directoryEMCALPi0->Get("EfficiencyPi0");
            graphEMCALPi0AccTimesEff                     = (TGraphAsymmErrors*)directoryEMCALPi0->Get("EffTimesAccPi0");
            histoEMCALPi0InvYieldStat                   = (TH1D*)directoryEMCALPi0->Get("CorrectedYieldPi0");
            graphEMCALPi0InvYieldStat                   = (TGraphAsymmErrors*)directoryEMCALPi0->Get("graphCorrectedYieldPi0");    
            graphEMCALPi0InvYieldSys                    = (TGraphAsymmErrors*)directoryEMCALPi0->Get("Pi0SystError");
            if (isNSD) {
                histoEMCALPi0InvYieldStat->Scale(scalingToNSD);
                graphEMCALPi0InvYieldStat               = ScaleGraph(graphEMCALPi0InvYieldStat,scalingToNSD);
                graphEMCALPi0InvYieldSys                = ScaleGraph(graphEMCALPi0InvYieldSys,scalingToNSD);
            }
            
        directoryEMCALEta                           = (TDirectory*)fileEMCAL->Get("EtapPb_5.023TeV"); 
            graphEMCALEtaMass                           = (TGraphAsymmErrors*)directoryEMCALEta->Get("Eta_Mass_data");
            graphEMCALEtaMass                           = ScaleGraph(graphEMCALEtaMass, 1000.);
            graphEMCALEtaFWHM                           = (TGraphAsymmErrors*)directoryEMCALEta->Get("Eta_Width_data");
            graphEMCALEtaFWHM                           = ScaleGraph(graphEMCALEtaFWHM, 1000.);
            graphEMCALEtaMassMC                         = (TGraphAsymmErrors*)directoryEMCALEta->Get("Eta_Mass_MC");
            graphEMCALEtaMassMC                         = ScaleGraph(graphEMCALEtaMassMC, 1000.);
            graphEMCALEtaFWHMMC                         = (TGraphAsymmErrors*)directoryEMCALEta->Get("Eta_Width_MC");
            graphEMCALEtaFWHMMC                         = ScaleGraph(graphEMCALEtaFWHMMC, 1000.);
            graphEMCALEtaAcc                            = (TGraphAsymmErrors*)directoryEMCALEta->Get("AcceptanceEta");
            graphEMCALEtaEffPt                          = (TGraphAsymmErrors*)directoryEMCALEta->Get("EfficiencyEta");
            graphEMCALEtaAccTimesEff                    = (TGraphAsymmErrors*)directoryEMCALEta->Get("EffTimesAccEta");
            histoEMCALEtaInvYieldStat                   = (TH1D*)directoryEMCALEta->Get("CorrectedYieldEta");
            graphEMCALEtaInvYieldStat                   = (TGraphAsymmErrors*)directoryEMCALEta->Get("graphCorrectedYieldEta");
            graphEMCALEtaInvYieldSys                    = (TGraphAsymmErrors*)directoryEMCALEta->Get("EtaSystError");
            if (isNSD) {
                histoEMCALEtaInvYieldStat->Scale(scalingToNSD);
                graphEMCALEtaInvYieldStat               = ScaleGraph(graphEMCALEtaInvYieldStat,scalingToNSD);
                graphEMCALEtaInvYieldSys                = ScaleGraph(graphEMCALEtaInvYieldSys,scalingToNSD);
            }
            histoEMCALEtaToPi0Stat                      = (TH1D*)directoryEMCALEta->Get("EtaToPi0YShiftedStatError");
            graphEMCALEtaToPi0Sys                       = (TGraphAsymmErrors*)directoryEMCALEta->Get("EtaToPi0YShiftedSystError");
            if (!histoEMCALEtaToPi0Stat){
                histoEMCALEtaToPi0Stat                  = (TH1D*)directoryEMCALEta->Get("EtaToPi0StatError");
            } 
            if (!graphEMCALEtaToPi0Sys){
                graphEMCALEtaToPi0Sys                   = (TGraphAsymmErrors*)directoryEMCALEta->Get("EtaToPi0SystError");
            } 
    }    
    graphEMCALPi0AccTimesEff                            = ScaleGraph(graphEMCALPi0AccTimesEff,0.98823);
    //************************** Read data for PHOS *****************************************************
    TFile* filePHOS                                 = new TFile(fileNamePHOS);
        TH1D* histoPHOSPi0InvYieldStat                  = (TH1D*)filePHOS->Get("hCor_stat");
        TH1D* histoPHOSPi0InvYieldSys                   = (TH1D*)filePHOS->Get("hCor_syst");
        TGraphAsymmErrors* graphPHOSPi0InvYieldStat     = new TGraphAsymmErrors(histoPHOSPi0InvYieldStat);
        graphPHOSPi0InvYieldStat->RemovePoint(0);
        graphPHOSPi0InvYieldStat->Print();
        TGraphAsymmErrors* graphPHOSPi0InvYieldSys      = new TGraphAsymmErrors(histoPHOSPi0InvYieldSys);
        graphPHOSPi0InvYieldSys->RemovePoint(0);
        if (isNSD) {
            histoPHOSPi0InvYieldStat->Scale(scalingToNSD);
            graphPHOSPi0InvYieldStat                    = ScaleGraph(graphPHOSPi0InvYieldStat,scalingToNSD);
            graphPHOSPi0InvYieldSys                     = ScaleGraph(graphPHOSPi0InvYieldSys,scalingToNSD);
        }        
        TH1D* histoPHOSPi0Mass                          = (TH1D*)filePHOS->Get("Gmean_Real");
        TH1D* histoPHOSPi0FWHMMeV                       = (TH1D*)filePHOS->Get("Gsigma_Real");
        TH1D* histoPHOSPi0TrueMass                      = (TH1D*)filePHOS->Get("Gmean_MC");
        TH1D* histoPHOSPi0TrueFWHMMeV                   = (TH1D*)filePHOS->Get("Gsigma_MC");
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
//         TH1D* histoPHOSPi0AccTimesEff                       = (TH1D*)directoryPHOSPi0->Get("yield1_int_CB");
//         histoPHOSPi0AccTimesEff->Scale(2*TMath::Pi());
//         TGraphAsymmErrors* graphPHOSPi0AccTimesEff          = new TGraphAsymmErrors(histoPHOSPi0AccTimesEff);


    //************************** Read data for Dalitz *************************************************
    TFile* fileDalitz                               = NULL;
    TDirectory* directoryDalitzPi0                  = NULL;
    TGraphAsymmErrors* graphDalitzPi0Mass           = NULL;
    TGraphAsymmErrors* graphDalitzPi0FWHM           = NULL;
    TGraphAsymmErrors* graphDalitzPi0MassMC         = NULL;
    TGraphAsymmErrors* graphDalitzPi0FWHMMC         = NULL;
    TGraphAsymmErrors* graphDalitzPi0cc             = NULL;
    TGraphAsymmErrors* graphDalitzPi0EffPt          = NULL;
    TGraphAsymmErrors* graphDalitzPi0InvYieldStat   = NULL;
    TH1D* histoDalitzPi0InvYieldStat                = NULL;
    TGraphAsymmErrors* graphDalitzPi0InvYieldSys    = NULL;
    TGraphAsymmErrors* graphDalitzPi0AccTimesEff     = NULL;

    if (!flagDalitzfile){  
        fileDalitz                                  = new TFile(fileNameDalitz.Data());
        directoryDalitzPi0                          = (TDirectory*)fileDalitz->Get("Pi0_pPb_5.023TeV_0-100%"); 
            TH1D* histoDalitzPi0Mass                    = (TH1D*)directoryDalitzPi0->Get("MassPi0");
            TH1D* histoDalitzPi0FWHMMeV                 = (TH1D*)directoryDalitzPi0->Get("FWHMPi0MeV");
            TH1D* histoDalitzPi0TrueMass                = (TH1D*)directoryDalitzPi0->Get("TrueMassPi0");
            TH1D* histoDalitzPi0TrueFWHMMeV             = (TH1D*)directoryDalitzPi0->Get("TrueFWHMPi0MeV");
            TH1D* histoDalitzPi0cc                      = (TH1D*)directoryDalitzPi0->Get("AcceptancePi0");
            TH1D* histoDalitzPi0TrueEffPt               = (TH1D*)directoryDalitzPi0->Get("EfficiencyPi0");
            histoDalitzPi0InvYieldStat                  = (TH1D*)directoryDalitzPi0->Get("CorrectedYieldPi0");
            graphDalitzPi0InvYieldStat                  = new TGraphAsymmErrors(histoDalitzPi0InvYieldStat);
            while (graphDalitzPi0InvYieldStat->GetY()[0] == 0) graphDalitzPi0InvYieldStat->RemovePoint(0);
            graphDalitzPi0InvYieldSys                   = (TGraphAsymmErrors*)directoryDalitzPi0->Get("Pi0SystError");
            while (graphDalitzPi0InvYieldSys->GetY()[0] == 0) graphDalitzPi0InvYieldSys->RemovePoint(0);
            if (isNSD) {
                histoDalitzPi0InvYieldStat->Scale(scalingToNSD);
                graphDalitzPi0InvYieldStat              = ScaleGraph(graphDalitzPi0InvYieldStat,scalingToNSD);
                graphDalitzPi0InvYieldSys               = ScaleGraph(graphDalitzPi0InvYieldSys,scalingToNSD);
            }

            TH1D* histoDalitzPi0AccTimesEff              = (TH1D*)histoDalitzPi0TrueEffPt->Clone("histoDalitzPi0AccTimesEff");
            histoDalitzPi0AccTimesEff->Multiply(histoDalitzPi0cc);
            // normalize to full acceptance (delta y and phi)
            histoDalitzPi0AccTimesEff->Scale(2*TMath::Pi()*1.6);
            histoDalitzPi0Mass->Scale(1000.);
            histoDalitzPi0TrueMass->Scale(1000.);

            graphDalitzPi0Mass                          = new TGraphAsymmErrors(histoDalitzPi0Mass);
            graphDalitzPi0FWHM                          = new TGraphAsymmErrors(histoDalitzPi0FWHMMeV);
            graphDalitzPi0MassMC                        = new TGraphAsymmErrors(histoDalitzPi0TrueMass);
            graphDalitzPi0FWHMMC                        = new TGraphAsymmErrors(histoDalitzPi0TrueFWHMMeV);
            graphDalitzPi0cc                            = new TGraphAsymmErrors(histoDalitzPi0cc);
            graphDalitzPi0EffPt                         = new TGraphAsymmErrors(histoDalitzPi0TrueEffPt);
            graphDalitzPi0AccTimesEff                    = new TGraphAsymmErrors(histoDalitzPi0AccTimesEff);
            
            while(graphDalitzPi0AccTimesEff->GetY()[0] == 0){
                graphDalitzPi0Mass->RemovePoint(0);
                graphDalitzPi0FWHM->RemovePoint(0);
                graphDalitzPi0MassMC->RemovePoint(0);
                graphDalitzPi0FWHMMC->RemovePoint(0);
                graphDalitzPi0cc->RemovePoint(0);
                graphDalitzPi0EffPt->RemovePoint(0);
                graphDalitzPi0AccTimesEff->RemovePoint(0);
            }
                
    } else {
        fileDalitz                                  = new TFile(fileNameDalitz.Data());
        directoryDalitzPi0                          = (TDirectory*)fileDalitz->Get("Pi0pPb_5.023TeV"); 
            graphDalitzPi0Mass                          = (TGraphAsymmErrors*)directoryDalitzPi0->Get("Pi0_Mass_data");
            graphDalitzPi0Mass                          = ScaleGraph(graphDalitzPi0Mass, 1000.);
            graphDalitzPi0FWHM                          = (TGraphAsymmErrors*)directoryDalitzPi0->Get("Pi0_Width_data");
            graphDalitzPi0FWHM                          = ScaleGraph(graphDalitzPi0FWHM, 1000.);
            graphDalitzPi0MassMC                        = (TGraphAsymmErrors*)directoryDalitzPi0->Get("Pi0_Mass_MC");
            graphDalitzPi0MassMC                        = ScaleGraph(graphDalitzPi0MassMC, 1000.);
            graphDalitzPi0FWHMMC                        = (TGraphAsymmErrors*)directoryDalitzPi0->Get("Pi0_Width_MC");
            graphDalitzPi0FWHMMC                        = ScaleGraph(graphDalitzPi0FWHMMC, 1000.);
            graphDalitzPi0cc                            = (TGraphAsymmErrors*)directoryDalitzPi0->Get("AcceptancePi0");
            graphDalitzPi0EffPt                         = (TGraphAsymmErrors*)directoryDalitzPi0->Get("EfficiencyPi0");
            graphDalitzPi0InvYieldStat                  = (TGraphAsymmErrors*)directoryDalitzPi0->Get("graphCorrectedYieldPi0");
            histoDalitzPi0InvYieldStat                  = (TH1D*)directoryDalitzPi0->Get("CorrectedYieldPi0");
            cout << "Pi0 stat Dalitz" << endl;
            graphDalitzPi0InvYieldStat->Print();
            graphDalitzPi0InvYieldSys                   = (TGraphAsymmErrors*)directoryDalitzPi0->Get("Pi0SystError");
            cout << "Pi0 sys Dalitz" << endl;
            graphDalitzPi0InvYieldSys->Print();
            if (isNSD) {
                histoDalitzPi0InvYieldStat->Scale(scalingToNSD);
                graphDalitzPi0InvYieldStat              = ScaleGraph(graphDalitzPi0InvYieldStat,scalingToNSD);
                graphDalitzPi0InvYieldSys               = ScaleGraph(graphDalitzPi0InvYieldSys,scalingToNSD);
            }
            graphDalitzPi0AccTimesEff                   = (TGraphAsymmErrors*)directoryDalitzPi0->Get("EffTimesAccPi0");

    }
    graphDalitzPi0AccTimesEff                           = ScaleGraph(graphDalitzPi0AccTimesEff,0.01174);    
 
    
    // *******************************************************************************************************
    // ************************** Loading charged pion results ***********************************************
    // *******************************************************************************************************        
    TFile* fileChargedPionInput                         = new TFile(fileNameChargedPion.Data());
        TH1D* histoChPiInvYieldPubStat                  = (TH1D*)fileChargedPionInput->Get("histoChargedPionPubStatpPb");
        TH1D* histoChPiInvYieldPubSyst                  = (TH1D*)fileChargedPionInput->Get("graphChargedPionPubStatpPb");

    // *******************************************************************************************************
    // ************************** Loading eta/pi0 compilation ************************************************
    // *******************************************************************************************************        
    TFile* fileEtaToPi0Compilation                          = new TFile(fileNameEtaToPi0WorldData.Data());
        TGraphErrors* graphPHENIXEtaToPi0200GeV             = (TGraphErrors*)fileEtaToPi0Compilation->Get("Phenix200GeV");
        TGraphAsymmErrors* graphALICEEtaToPi07TeV           = (TGraphAsymmErrors*)fileEtaToPi0Compilation->Get("Alice7TeV");
        TGraphAsymmErrors* graphALICEEtaToPi02760GeV        = (TGraphAsymmErrors*)fileEtaToPi0Compilation->Get("Alice2760GeV");
        
    
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
    statErrorCollectionPi0[0]       = (TH1D*)histoPCMPi0InvYieldStat->Clone("statErrPCMPi0");
    statErrorCollectionPi0[1]       = (TH1D*)histoPHOSPi0InvYieldStat->Clone("statErrPHOSPi0");
    statErrorCollectionPi0[2]       = (TH1D*)histoEMCALPi0InvYieldStat->Clone("statErrEMCALPi0");
    statErrorCollectionPi0[5]       = (TH1D*)histoDalitzPi0InvYieldStat->Clone("statErrDalitzPi0");
    
    // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)    
    // for systematic error from respective method
    TGraphAsymmErrors* sysErrorCollectionPi0[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorCollectionPi0[i]    = NULL;
    }    
    sysErrorCollectionPi0[0]        = (TGraphAsymmErrors*)graphPCMPi0InvYieldSys->Clone("sysErrPCMPi0");
    sysErrorCollectionPi0[1]        = (TGraphAsymmErrors*)graphPHOSPi0InvYieldSys->Clone("sysErrPHOSPi0");
    sysErrorCollectionPi0[2]        = (TGraphAsymmErrors*)graphEMCALPi0InvYieldSys->Clone("sysErrEMCALPi0");
    sysErrorCollectionPi0[5]        = (TGraphAsymmErrors*)graphDalitzPi0InvYieldSys->Clone("sysErrDalitzPi0");

    // Definition of final pt binning (has to be set manually)
    Double_t xPtLimitsPi0[100];
    Int_t maxNBinsPi0               = GetBinning( xPtLimitsPi0, "Pi0", "pPb_5.023TeV", 20 );
    for (Int_t i = 0; i< maxNBinsPi0; i++){
        cout << i << ": "<< xPtLimitsPi0[i] <<" - " << xPtLimitsPi0[i+1]<< ", " <<endl;
    }
        
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsPi0[11]            = { 0,  6,  0,  0,  0,
                                        3,  0,  0,  0,  0, 
                                        0};
    Int_t offSetsPi0Sys[11]         = { 1,  7,  9,  0,  6,
                                        4,  0,  0,  0,  21, 
                                        0};
    Int_t offSetPi0Shifting[11]     = { 0,  6,  8,  0,  5,
                                        3,  0,  0,  0,  21,
                                        0 };
    Int_t nComBinsPi0Shifting[11]   = { 29, 31, 30, 0,  31,
                                        17,  0,  0,  0,  0,
                                        0 };
                                        

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

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************
                                
    TGraph* graphWeightsWOPEPi0[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsWOPEPi0[i]         = NULL;
    }    
                                
    // Declaration & calculation of combined spectrum                            
    TString fileNamePi0OutputWeightingWOPE           = Form("%s/Pi0_WeightingMethodWOPE.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombPi0InvYieldWOPEStat  = NULL;
    TGraphAsymmErrors* graphCombPi0InvYieldWOPESys   = NULL;
    TGraphAsymmErrors* graphCombPi0InvYieldWOPETot   = CombinePtPointsSpectraFullCorrMat(   statErrorCollectionPi0,    sysErrorCollectionPi0,     
                                                                                            xPtLimitsPi0, maxNBinsPi0,
                                                                                            offSetsPi0, offSetsPi0Sys,
                                                                                            graphCombPi0InvYieldWOPEStat, graphCombPi0InvYieldWOPESys,
                                                                                            fileNamePi0OutputWeightingWOPE, "pPb_5.023TeV", "Pi0", kTRUE, 
                                                                                            NULL, fileNameCorrFactors
                                                                                        );
    
    
    if (graphCombPi0InvYieldWOPETot == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }
    
//     return;
    
    while (graphCombPi0InvYieldWOPEStat->GetX()[0] < 0.3){
        graphCombPi0InvYieldWOPEStat->RemovePoint(0);
    }
    while (graphCombPi0InvYieldWOPETot->GetX()[0] < 0.3){    
        graphCombPi0InvYieldWOPETot->RemovePoint(0);
    }
    while (graphCombPi0InvYieldWOPESys->GetX()[0] < 0.3){    
        graphCombPi0InvYieldWOPESys->RemovePoint(0);
    }    
    graphCombPi0InvYieldWOPETot->Print();

    // Reading weights from output file for plotting
    ifstream fileWeightsPi0ReadWOPE;
    fileWeightsPi0ReadWOPE.open(fileNamePi0OutputWeightingWOPE,ios_base::in);
    cout << "reading" << fileNamePi0OutputWeightingWOPE << endl;
    Double_t xValuesPi0ReadWOPE[50];
    Double_t weightsPi0ReadWOPE[11][50];
    Int_t availablePi0MeasWOPE[11]    = {   -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1,
                                            -1};
    Int_t nMeasSetPi0WOPE             = 4;
    Int_t nPtBinsPi0ReadWOPE          = 0;
    while(!fileWeightsPi0ReadWOPE.eof() && nPtBinsPi0ReadWOPE < 50){
        TString garbage             = "";
        if (nPtBinsPi0ReadWOPE == 0){
            fileWeightsPi0ReadWOPE >> garbage ;//>> availablePi0MeasWOPE[0] >> availablePi0MeasWOPE[1] >> availablePi0MeasWOPE[2] >> availablePi0MeasWOPE[3];
            for (Int_t i = 0; i < nMeasSetPi0WOPE; i++){
                fileWeightsPi0ReadWOPE >> availablePi0MeasWOPE[i] ;
            }    
            cout << "read following measurements: "; 
            for (Int_t i = 0; i < 11; i++){
                cout << availablePi0MeasWOPE[i] << "\t" ;
            }    
            cout << endl;
        } else {
            fileWeightsPi0ReadWOPE >> xValuesPi0ReadWOPE[nPtBinsPi0ReadWOPE-1];
            for (Int_t i = 0; i < nMeasSetPi0WOPE; i++){
                fileWeightsPi0ReadWOPE >> weightsPi0ReadWOPE[availablePi0MeasWOPE[i]][nPtBinsPi0ReadWOPE-1] ;
            }    
            cout << "read: "<<  nPtBinsPi0ReadWOPE << "\t"<< xValuesPi0ReadWOPE[nPtBinsPi0ReadWOPE-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetPi0WOPE; i++){
                cout << weightsPi0ReadWOPE[availablePi0MeasWOPE[i]][nPtBinsPi0ReadWOPE-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsPi0ReadWOPE++;
    }
    nPtBinsPi0ReadWOPE                  = nPtBinsPi0ReadWOPE-2 ;
    fileWeightsPi0ReadWOPE.close();

    for (Int_t i = 0; i < nMeasSetPi0WOPE; i++){
        graphWeightsWOPEPi0[availablePi0MeasWOPE[i]]                        = new TGraph(nPtBinsPi0ReadWOPE,xValuesPi0ReadWOPE,weightsPi0ReadWOPE[availablePi0MeasWOPE[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsPi0ReadWOPE; n++){
            if (graphWeightsWOPEPi0[availablePi0MeasWOPE[i]]->GetY()[bin] == 0) graphWeightsWOPEPi0[availablePi0MeasWOPE[i]]->RemovePoint(bin);
            else bin++;
        }    
    }    

    // **********************************************************************************************************************
    // ************************ Adding PCM-EMC measurements to matrix & calculating relativ errors ************************
    // **********************************************************************************************************************
    
    statErrorCollectionPi0[4]               = (TH1D*)histoPCMEMCALPi0InvYieldStat->Clone("statErrPCMEMCALPi0");
    sysErrorCollectionPi0[4]                = (TGraphAsymmErrors*)graphPCMEMCALPi0InvYieldSys->Clone("sysErrPCMEMCALPi0");

    if (statErrorCollectionPi0[4]){
        statErrorRelCollectionPi0[4]        = new TGraphAsymmErrors(statErrorCollectionPi0[4]);
        statErrorRelCollectionPi0[4]        = CalculateRelErrUpAsymmGraph( statErrorRelCollectionPi0[4], Form("relativeStatErrorPi0_%s", nameMeasGlobal[4].Data()));
    }        
    if (sysErrorCollectionPi0[4]) 
        sysErrorRelCollectionPi0[4]         = CalculateRelErrUpAsymmGraph( sysErrorCollectionPi0[4], Form("relativeSysErrorPi0_%s", nameMeasGlobal[4].Data()));

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************
                                
    TGraph* graphWeightsPi0[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsPi0[i]         = NULL;
    }    
                                
    // Declaration & calculation of combined spectrum                            
    TString fileNamePi0OutputWeighting           = Form("%s/Pi0_WeightingMethod.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombPi0InvYieldStat  = NULL;
    TGraphAsymmErrors* graphCombPi0InvYieldSys   = NULL;
    TGraphAsymmErrors* graphCombPi0InvYieldTot   = CombinePtPointsSpectraFullCorrMat(   statErrorCollectionPi0,    sysErrorCollectionPi0,     
                                                                                        xPtLimitsPi0, maxNBinsPi0,
                                                                                        offSetsPi0, offSetsPi0Sys,
                                                                                        graphCombPi0InvYieldStat, graphCombPi0InvYieldSys,
                                                                                        fileNamePi0OutputWeighting, "pPb_5.023TeV", "Pi0", kTRUE, 
                                                                                        NULL, fileNameCorrFactors
                                                                                    );
    
    
    if (graphCombPi0InvYieldTot == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }
    
//     return;
    
    while (graphCombPi0InvYieldStat->GetX()[0] < 0.3){
        graphCombPi0InvYieldStat->RemovePoint(0);
    }
    while (graphCombPi0InvYieldTot->GetX()[0] < 0.3){    
        graphCombPi0InvYieldTot->RemovePoint(0);
    }
    while (graphCombPi0InvYieldSys->GetX()[0] < 0.3){    
        graphCombPi0InvYieldSys->RemovePoint(0);
    }    
    graphCombPi0InvYieldTot->Print();

    // Reading weights from output file for plotting
    ifstream fileWeightsPi0Read;
    fileWeightsPi0Read.open(fileNamePi0OutputWeighting,ios_base::in);
    cout << "reading" << fileNamePi0OutputWeighting << endl;
    Double_t xValuesPi0Read[50];
    Double_t weightsPi0Read[11][50];
    Int_t availablePi0Meas[11]      = { -1, -1, -1, -1, -1,
                                        -1, -1, -1, -1, -1,
                                        -1};
    Int_t nMeasSetPi0               = 5;
    Int_t nPtBinsPi0Read            = 0;
    while(!fileWeightsPi0Read.eof() && nPtBinsPi0Read < 50){
        TString garbage             = "";
        if (nPtBinsPi0Read == 0){
            fileWeightsPi0Read >> garbage ;//>> availablePi0Meas[0] >> availablePi0Meas[1] >> availablePi0Meas[2] >> availablePi0Meas[3];
            for (Int_t i = 0; i < nMeasSetPi0; i++){
                fileWeightsPi0Read >> availablePi0Meas[i] ;
            }    
            cout << "read following measurements: "; 
            for (Int_t i = 0; i < 11; i++){
                cout << availablePi0Meas[i] << "\t" ;
            }    
            cout << endl;
        } else {
            fileWeightsPi0Read >> xValuesPi0Read[nPtBinsPi0Read-1];
            for (Int_t i = 0; i < nMeasSetPi0; i++){
                fileWeightsPi0Read >> weightsPi0Read[availablePi0Meas[i]][nPtBinsPi0Read-1] ;
            }    
            cout << "read: "<<  nPtBinsPi0Read << "\t"<< xValuesPi0Read[nPtBinsPi0Read-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetPi0; i++){
                cout << weightsPi0Read[availablePi0Meas[i]][nPtBinsPi0Read-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsPi0Read++;
    }
    nPtBinsPi0Read                  = nPtBinsPi0Read-2 ;
    fileWeightsPi0Read.close();

    for (Int_t i = 0; i < nMeasSetPi0; i++){
        graphWeightsPi0[availablePi0Meas[i]]                        = new TGraph(nPtBinsPi0Read,xValuesPi0Read,weightsPi0Read[availablePi0Meas[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsPi0Read; n++){
            if (graphWeightsPi0[availablePi0Meas[i]]->GetY()[bin] == 0) graphWeightsPi0[availablePi0Meas[i]]->RemovePoint(bin);
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
    histo2DPi0Weights = new TH2F("histo2DPi0Weights","histo2DPi0Weights",11000,0.23,31.,1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DPi0Weights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DPi0Weights->GetXaxis()->SetMoreLogLabels();
    histo2DPi0Weights->GetXaxis()->SetLabelOffset(-0.01);
    histo2DPi0Weights->Draw("copy");
    
        TLegend* legendWeights   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetPi0WOPE), 32);
        for (Int_t i = 0; i < nMeasSetPi0WOPE; i++){
            DrawGammaSetMarkerTGraph(graphWeightsWOPEPi0[availablePi0MeasWOPE[i]], markerStyleDet[availablePi0MeasWOPE[i]], markerSizeDet[availablePi0MeasWOPE[i]]*0.5, colorDet[availablePi0MeasWOPE[i]] , colorDet[availablePi0MeasWOPE[i]]);
            graphWeightsWOPEPi0[availablePi0MeasWOPE[i]]->Draw("p,same,z");
            legendWeights->AddEntry(graphWeightsWOPEPi0[availablePi0MeasWOPE[i]],nameMeasGlobalLabel[availablePi0MeasWOPE[i]],"p");
        }    
        legendWeights->Draw();

        TLatex *labelWeightsEnergy      = new TLatex(0.7,0.20,collisionSystempPb.Data());
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
        
    canvasWeights->SaveAs(Form("%s/Pi0_Weights_WOPCMEMC.%s",outputDir.Data(),suffix.Data()));
    
    histo2DPi0Weights->Draw("copy");
        TLegend* legendWeights2   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetPi0), 32);
        for (Int_t i = 0; i < nMeasSetPi0; i++){
            DrawGammaSetMarkerTGraph(graphWeightsPi0[availablePi0Meas[i]], markerStyleDet[availablePi0Meas[i]], markerSizeDet[availablePi0Meas[i]]*0.5, colorDet[availablePi0Meas[i]] , colorDet[availablePi0Meas[i]]);
            graphWeightsPi0[availablePi0Meas[i]]->Draw("p,same,z");
            legendWeights2->AddEntry(graphWeightsPi0[availablePi0Meas[i]],nameMeasGlobalLabel[availablePi0Meas[i]],"p");
        }    
        legendWeights2->Draw();

        labelWeightsEnergy->Draw();
        labelWeightsPi0->Draw();

//      DrawGammaLines(0.23, 70. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);
        
    canvasWeights->SaveAs(Form("%s/Pi0_Weights.%s",outputDir.Data(),suffix.Data()));
    
    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************
    
    TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErr->SetLogx();
   
    TH2F * histo2DRelSysErr;
    histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.23,31.,1000,0,45.5);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErr->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRelSysErr->Draw("copy");

        TLegend* legendRelSysErr        = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetPi0), 0.95, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetPi0WOPE; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionPi0[availablePi0MeasWOPE[i]], markerStyleDet[availablePi0MeasWOPE[i]], markerSizeDet[availablePi0MeasWOPE[i]]*0.5, colorDet[availablePi0MeasWOPE[i]],
                                     colorDet[availablePi0MeasWOPE[i]]);
            sysErrorRelCollectionPi0[availablePi0MeasWOPE[i]]->Draw("p,same,z");
            legendRelSysErr->AddEntry(sysErrorRelCollectionPi0[availablePi0MeasWOPE[i]],nameMeasGlobalLabel[availablePi0MeasWOPE[i]],"p");
        }    
        legendRelSysErr->Draw();

        TLatex *labelRelSysErrEnergy    = new TLatex(0.15,0.89,collisionSystempPb.Data());
        SetStyleTLatex( labelRelSysErrEnergy, 0.85*textSizeLabelsPixel, 4, 1, 43);
        labelRelSysErrEnergy->Draw();
        TLatex *labelRelSysErrPi0       = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelSysErrPi0, 0.85*textSizeLabelsPixel, 4, 1, 43);
        labelRelSysErrPi0->Draw();
        
    canvasRelSysErr->SaveAs(Form("%s/Pi0_RelSysErr_WOPCMEMC.%s",outputDir.Data(),suffix.Data()));
    
    histo2DRelSysErr->Draw("copy");
        TLegend* legendRelSysErr2       = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetPi0), 0.95, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetPi0; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionPi0[availablePi0Meas[i]], markerStyleDet[availablePi0Meas[i]], markerSizeDet[availablePi0Meas[i]]*0.5, colorDet[availablePi0Meas[i]],
                                     colorDet[availablePi0Meas[i]]);
            sysErrorRelCollectionPi0[availablePi0Meas[i]]->Draw("p,same,z");
            legendRelSysErr2->AddEntry(sysErrorRelCollectionPi0[availablePi0Meas[i]],nameMeasGlobalLabel[availablePi0Meas[i]],"p");
        }    
        legendRelSysErr2->Draw();

        for (Int_t i = 0; i < nMeasSetPi0; i++){
            sysErrorRelCollectionPi0[availablePi0Meas[i]]->Draw("p,same,z");
        }    
        legendRelSysErr2->Draw();

        labelRelSysErrEnergy->Draw();
        labelRelSysErrPi0->Draw();
        
    canvasRelSysErr->SaveAs(Form("%s/Pi0_RelSysErr.%s",outputDir.Data(),suffix.Data()));
    
    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************
    
    TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelStatErr->SetLogx();
   
    TH2F * histo2DRelStatErr;
    histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.23,31.,1000,0,45.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErr->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRelStatErr->Draw("copy");
        TLegend* legendRelStatErr       = GetAndSetLegend2(0.34, 0.92-(0.035*nMeasSetPi0WOPE), 0.45, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetPi0WOPE; i++){
            DrawGammaSetMarkerTGraph(statErrorRelCollectionPi0[availablePi0MeasWOPE[i]], markerStyleDet[availablePi0MeasWOPE[i]], markerSizeDet[availablePi0MeasWOPE[i]]*0.5, colorDet[availablePi0MeasWOPE[i]],
                                     colorDet[availablePi0MeasWOPE[i]]);
            statErrorRelCollectionPi0[availablePi0MeasWOPE[i]]->Draw("p,same,z");
            legendRelStatErr->AddEntry(statErrorRelCollectionPi0[availablePi0MeasWOPE[i]],nameMeasGlobalLabel[availablePi0MeasWOPE[i]],"p");
        }    
        legendRelStatErr->Draw();

        TLatex *labelRelStatErrEnergy   = new TLatex(0.75,0.89,collisionSystempPb.Data());
        SetStyleTLatex( labelRelStatErrEnergy, 0.85*textSizeLabelsPixel, 4, 1, 43);
        labelRelStatErrEnergy->Draw();
        TLatex *labelRelStatErrPi0      = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelStatErrPi0, 0.85*textSizeLabelsPixel, 4, 1, 43);
        labelRelStatErrPi0->Draw();
        
    canvasRelStatErr->SaveAs(Form("%s/Pi0_RelStatErr_WOPCMEMC.%s",outputDir.Data(),suffix.Data()));

    histo2DRelStatErr->Draw("copy");
        TLegend* legendRelStatErr2       = GetAndSetLegend2(0.34, 0.92-(0.035*nMeasSetPi0), 0.45, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetPi0; i++){
            DrawGammaSetMarkerTGraph(statErrorRelCollectionPi0[availablePi0Meas[i]], markerStyleDet[availablePi0Meas[i]], markerSizeDet[availablePi0Meas[i]]*0.5, colorDet[availablePi0Meas[i]],
                                     colorDet[availablePi0Meas[i]]);
            statErrorRelCollectionPi0[availablePi0Meas[i]]->Draw("p,same,z");
            legendRelStatErr2->AddEntry(statErrorRelCollectionPi0[availablePi0Meas[i]],nameMeasGlobalLabel[availablePi0Meas[i]],"p");
        }    
        legendRelStatErr2->Draw();

        labelRelStatErrEnergy->Draw();
        labelRelStatErrPi0->Draw();
        
    canvasRelStatErr->SaveAs(Form("%s/Pi0_RelStatErr.%s",outputDir.Data(),suffix.Data()));
    
    //  *********************************************************************************************************************
    //  ************************ Visualize relative total errors of different combination methods Pi0 ***********************
    //  *********************************************************************************************************************
    TGraphAsymmErrors* graphCombPi0InvYieldWOPERelStat  = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldWOPEStat, "relativeStatErrorPi0_MethodWOPE");
    TGraphAsymmErrors* graphCombPi0InvYieldWOPERelSys   = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldWOPESys, "relativeSysErrorPi0_MethodWOPE");
    TGraphAsymmErrors* graphCombPi0InvYieldWOPERelTot   = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldWOPETot, "relativeTotalErrorPi0_MethodWOPE");
    TGraphAsymmErrors* graphCombPi0InvYieldRelStat      = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldStat, "relativeStatErrorPi0_Method");
    TGraphAsymmErrors* graphCombPi0InvYieldRelSys       = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldSys, "relativeSysErrorPi0_Method");
    TGraphAsymmErrors* graphCombPi0InvYieldRelTot       = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldTot, "relativeTotalErrorPi0_Method");
    
    TCanvas* canvasRelTotErr            = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelTotErr->SetLogx();
   
    TH2F * histo2DRelTotErrPi0;
    histo2DRelTotErrPi0                 = new TH2F("histo2DRelTotErrPi0","histo2DRelTotErrPi0",11000,0.23,31.,1000,0,45.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrPi0, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrPi0->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRelTotErrPi0->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldWOPERelTot, markerStyleComb, markerSizeComb, kBlue+2 , kBlue+2);
        graphCombPi0InvYieldWOPERelTot->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot, markerStyleComb+4, markerSizeComb, kRed+2 , kRed+2);
        graphCombPi0InvYieldRelTot->Draw("p,same,z");

        TLegend* legendRelTotErr1       = GetAndSetLegend2(0.20, 0.92-(0.035*2), 0.45, 0.92, 32);
        legendRelTotErr1->AddEntry(graphCombPi0InvYieldWOPERelTot,"PCM, Dalitz, EMC, PHOS","p");
        legendRelTotErr1->AddEntry(graphCombPi0InvYieldRelTot,"All","p");
        legendRelTotErr1->Draw();

        TLatex *labelRelTotErrEnergy    = new TLatex(0.75,0.89,collisionSystempPb.Data());
        SetStyleTLatex( labelRelTotErrEnergy, 0.85*textSizeLabelsPixel, 4, 1, 43);        
        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrPi0       = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelTotErrPi0, 0.85*textSizeLabelsPixel, 4, 1, 43);
        labelRelTotErrPi0->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/Pi0_TotErr_Comp.%s",outputDir.Data(),suffix.Data()));
    histo2DRelTotErrPi0->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldWOPERelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombPi0InvYieldWOPERelTot->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldWOPERelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombPi0InvYieldWOPERelStat->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldWOPERelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombPi0InvYieldWOPERelSys->SetLineStyle(7);
        graphCombPi0InvYieldWOPERelSys->Draw("l,x0,same,e1");

        TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
        legendRelTotErr3->AddEntry(graphCombPi0InvYieldWOPERelTot,"tot","p");
        legendRelTotErr3->AddEntry(graphCombPi0InvYieldWOPERelStat,"stat","l");
        legendRelTotErr3->AddEntry(graphCombPi0InvYieldWOPERelSys,"sys","l");
        legendRelTotErr3->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrPi0->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/Pi0_Reldecomp_WOPCMEMC.%s",outputDir.Data(),suffix.Data()));
    
    histo2DRelTotErrPi0->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombPi0InvYieldRelTot->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombPi0InvYieldRelStat->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombPi0InvYieldRelSys->SetLineStyle(7);
        graphCombPi0InvYieldRelSys->Draw("l,x0,same,e1");
        legendRelTotErr3->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrPi0->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/Pi0_Reldecomp.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ************************************* Calculating bin shifted spectra & fitting **************************************
    // **********************************************************************************************************************
    
    // Cloning spectra
    TGraphAsymmErrors* graphCombPi0InvYieldTotUnshi         = (TGraphAsymmErrors*)graphCombPi0InvYieldTot->Clone("Pi0Unshifted"); 
    TGraphAsymmErrors* graphCombPi0InvYieldStatUnshi        = (TGraphAsymmErrors*)graphCombPi0InvYieldStat->Clone("Pi0UnshiftedStat"); 
    TGraphAsymmErrors* graphCombPi0InvYieldSysUnshi         = (TGraphAsymmErrors*)graphCombPi0InvYieldSys->Clone("Pi0UnshiftedSys"); 

    TGraphAsymmErrors* graphPCMPi0InvYieldStatUnshi         = (TGraphAsymmErrors*)graphPCMPi0InvYieldStat->Clone("Pi0UnshiftedStatPCM"); 
    TGraphAsymmErrors* graphPCMPi0InvYieldSysUnshi          = (TGraphAsymmErrors*)graphPCMPi0InvYieldSys->Clone("Pi0UnshiftedSysPCM"); 

    TGraphAsymmErrors* graphPHOSPi0InvYieldStatUnshi        = (TGraphAsymmErrors*)graphPHOSPi0InvYieldStat->Clone("Pi0UnshiftedStatPHOS"); 
    TGraphAsymmErrors* graphPHOSPi0InvYieldSysUnshi         = (TGraphAsymmErrors*)graphPHOSPi0InvYieldSys->Clone("Pi0UnshiftedSysPHOS"); 
    
    TGraphAsymmErrors* graphEMCALPi0InvYieldStatUnshi       = (TGraphAsymmErrors*)graphEMCALPi0InvYieldStat->Clone("Pi0UnshiftedStatEMCAL"); 
    TGraphAsymmErrors* graphEMCALPi0InvYieldSysUnshi        = (TGraphAsymmErrors*)graphEMCALPi0InvYieldSys->Clone("Pi0UnshiftedSysEMCAL"); 

    TGraphAsymmErrors* graphPCMEMCALPi0InvYieldStatUnshi    = (TGraphAsymmErrors*)graphPCMEMCALPi0InvYieldStat->Clone("Pi0UnshiftedStatPCMEMCAL"); 
    TGraphAsymmErrors* graphPCMEMCALPi0InvYieldSysUnshi     = (TGraphAsymmErrors*)graphPCMEMCALPi0InvYieldSys->Clone("Pi0UnshiftedSysPCMEMCAL"); 

    TGraphAsymmErrors* graphDalitzPi0InvYieldStatUnshi      = (TGraphAsymmErrors*)graphDalitzPi0InvYieldStat->Clone("Pi0UnshiftedStatDalitz"); 
    TGraphAsymmErrors* graphDalitzPi0InvYieldSysUnshi       = (TGraphAsymmErrors*)graphDalitzPi0InvYieldSys->Clone("Pi0UnshiftedSysDalitz"); 
   
    // fitting spectrum with intial parameters
    // Two component model fit from Bylinkin
    TF1* fitTCMDecomposedLPi0           = FitObject("tcmlow","twoCompModelPi0_DecL", "Pi0", NULL, 0.3, 2.);
    TF1* fitTCMDecomposedHPi0           = FitObject("tcmhigh","twoCompModelPi0_DecH", "Pi0", NULL, 4, 50.);
    fitTCMDecomposedLPi0->SetParameters(graphCombPi0InvYieldTot->GetY()[2],0.3);
    graphCombPi0InvYieldStat->Fit(fitTCMDecomposedLPi0,"QNRMEX0+","",0.4,2.);
    graphCombPi0InvYieldTot->Fit(fitTCMDecomposedHPi0,"QNRMEX0+","",4,50);
    fitTCMDecomposedHPi0->SetParameters(graphCombPi0InvYieldTot->GetY()[2],0.8, 2);    
    
    cout << WriteParameterToFile(fitTCMDecomposedLPi0)<< endl;    
    fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedLPi0)<< endl;    
    cout << WriteParameterToFile(fitTCMDecomposedHPi0)<< endl;    
    fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedHPi0)<< endl;    
    
    Double_t paramTCMPi0New[5]          = { 5.265,0.33,
                                            1.9,0.46,3.1};
    TF1* fitTCMInvYieldPi0              = FitObject("tcm","fitTCMInvYieldPi02760GeV","Pi0",graphCombPi0InvYieldStat,0.3,20. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
    
    fitTCMDecomposedLPi0->SetParameter(0, fitTCMInvYieldPi0->GetParameter(0));
    fitTCMDecomposedLPi0->SetParameter(1, fitTCMInvYieldPi0->GetParameter(1));
    fitTCMDecomposedHPi0->SetParameter(0, fitTCMInvYieldPi0->GetParameter(2));
    fitTCMDecomposedHPi0->SetParameter(1, fitTCMInvYieldPi0->GetParameter(3));
    fitTCMDecomposedHPi0->SetParameter(2, fitTCMInvYieldPi0->GetParameter(4));
    
    // Tsallis fit 
    Double_t paramGraphPi0[3]           = {1.0e12, 8., 0.13};
    TF1* fitInvYieldPi0                 = FitObject("l","fitInvYieldPi02760GeV","Pi0",histoEMCALPi0InvYieldStat,0.3,20.,paramGraphPi0,"QNRMEX0+");
    TF1* fitInvYieldPi0Graph            = (TF1*)fitInvYieldPi0->Clone("fitInvYieldPi02760GeVGraph"); 

    // *************************************************************************************************************
    // Shift graphs in X direction if desired
    // *************************************************************************************************************    
    if(bWCorrection.Contains("X")){
//         TF1* fitShiftingPi0            = FitObject("tmpt","ShiftingPi0","Pi0");
//         fitShiftingPi0->SetParameters(fitInvYieldPi0->GetParameter(0),fitInvYieldPi0->GetParameter(1), fitInvYieldPi0->GetParameter(2));
        TF1* fitShiftingPi0            = FitObject("tcmpt","ShiftingPi0","Pi0");
        fitShiftingPi0->SetParameters(fitTCMInvYieldPi0->GetParameter(0),fitTCMInvYieldPi0->GetParameter(1), fitTCMInvYieldPi0->GetParameter(2), fitTCMInvYieldPi0->GetParameter(3),fitTCMInvYieldPi0->GetParameter(4));
        
        TGraphAsymmErrors* graphCombPi0InvYieldTotNoShift = (TGraphAsymmErrors*) graphCombPi0InvYieldTot->Clone("Pi0_NoShift");

        graphCombPi0InvYieldTot             = ApplyXshift(graphCombPi0InvYieldTot, fitShiftingPi0);
        cout << "comb" << endl;
        graphCombPi0InvYieldStat->Print();
        graphCombPi0InvYieldStat            = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot, 
                                                                            graphCombPi0InvYieldStat, 
                                                                            fitShiftingPi0,
                                                                            0, graphCombPi0InvYieldStat->GetN());
        graphCombPi0InvYieldSys             = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot, 
                                                                            graphCombPi0InvYieldSys, 
                                                                            fitShiftingPi0, 
                                                                            0, graphCombPi0InvYieldSys->GetN());
        cout << "PCM" << endl;        
        graphPCMPi0InvYieldStat             = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot,
                                                                            graphPCMPi0InvYieldStat,
                                                                            fitShiftingPi0, 
                                                                            offSetPi0Shifting[0], nComBinsPi0Shifting[0]);
        graphPCMPi0InvYieldSys              = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot, 
                                                                            graphPCMPi0InvYieldSys, 
                                                                            fitShiftingPi0, 
                                                                            offSetPi0Shifting[0], nComBinsPi0Shifting[0]);
        cout << "PHOS" << endl;        
        graphPHOSPi0InvYieldStat            = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot, 
                                                                            graphPHOSPi0InvYieldStat, 
                                                                            fitShiftingPi0,
                                                                            offSetPi0Shifting[1], nComBinsPi0Shifting[1]);
        graphPHOSPi0InvYieldSys             = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot, 
                                                                            graphPHOSPi0InvYieldSys, 
                                                                            fitShiftingPi0,
                                                                            offSetPi0Shifting[1], nComBinsPi0Shifting[1]);

        cout << "EMC-EMC" << endl;        
        graphEMCALPi0InvYieldStat           = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot, 
                                                                            graphEMCALPi0InvYieldStat, 
                                                                            fitShiftingPi0,
                                                                            offSetPi0Shifting[2], nComBinsPi0Shifting[2]);
        graphEMCALPi0InvYieldSys            = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot, 
                                                                            graphEMCALPi0InvYieldSys, 
                                                                            fitShiftingPi0,
                                                                            offSetPi0Shifting[2], nComBinsPi0Shifting[2]);
        cout << "PCM-EMC" << endl;        
        graphPCMEMCALPi0InvYieldStat        = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot,
                                                                            graphPCMEMCALPi0InvYieldStat,
                                                                            fitShiftingPi0, 
                                                                            offSetPi0Shifting[4], nComBinsPi0Shifting[4]);
        graphPCMEMCALPi0InvYieldSys         = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot, 
                                                                            graphPCMEMCALPi0InvYieldSys, 
                                                                            fitShiftingPi0, 
                                                                            offSetPi0Shifting[4], nComBinsPi0Shifting[4]);
    
        cout << "Dalitz" << endl;        
        graphDalitzPi0InvYieldStat          = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot,
                                                                            graphDalitzPi0InvYieldStat,
                                                                            fitShiftingPi0, 
                                                                            offSetPi0Shifting[5], nComBinsPi0Shifting[5]);
        graphDalitzPi0InvYieldSys           = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot, 
                                                                            graphDalitzPi0InvYieldSys, 
                                                                            fitShiftingPi0, 
                                                                            offSetPi0Shifting[5], nComBinsPi0Shifting[5]);

        //***************************************************************************************************************
        //************************************Plotting binshift corrections *********************************************
        //***************************************************************************************************************

        TCanvas* canvasShift = new TCanvas("canvasShift","",0,0,1000,900);// gives the page size
        DrawGammaCanvasSettings( canvasShift, 0.12, 0.017, 0.015, 0.08);

        Size_t textSizeSpectra          = 0.04;
        TH1F * histoBinShift = new TH1F("histoBinShift","histoBinShift",1000,0.,20);
        SetStyleHistoTH1ForGraphs(histoBinShift, "#it{p}_{T} (GeV/#it{c})","bin shifted (X) / no shift",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.2);
        histoBinShift->GetYaxis()->SetRangeUser(0.95,1.05);
        histoBinShift->DrawCopy();

            Int_t numberPoints      = graphCombPi0InvYieldTotNoShift->GetN();
            Double_t *xPoint        = graphCombPi0InvYieldTotNoShift->GetX();
            Double_t* xvalueErrUp   = graphCombPi0InvYieldTotNoShift->GetEXhigh();
            Double_t* xvalueErrLow  = graphCombPi0InvYieldTotNoShift->GetEXlow();
            Double_t *xPointShift= graphCombPi0InvYieldTot->GetX();
            for (Int_t i=0; i<numberPoints; i++) {
            graphCombPi0InvYieldTotNoShift->SetPoint(i,xPoint[i],xPointShift[i]/xPoint[i]);
            graphCombPi0InvYieldTotNoShift->SetPointError(i,xvalueErrLow[i],xvalueErrUp[i],0,0);
            }
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldTotNoShift, markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphCombPi0InvYieldTotNoShift->Draw("p same");

            TLatex *labelRatioToFitBinShift   = new TLatex(0.94, 0.91, collisionSystempPb.Data());
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
        DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.09);
        canvasDummy2->SetLogy();
        canvasDummy2->SetLogx();
        TH2F * histo2DDummy2;
        histo2DDummy2               = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.53,31.,1000,1e-9,10e1);
        SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
        histo2DDummy2->DrawCopy(); 

       
            DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvYieldSys, markerStyleDet[0] ,markerSizeDet[0]/2, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
            graphPCMPi0InvYieldSys->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphPHOSPi0InvYieldSys, markerStyleDet[1] ,markerSizeDet[1]/2, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
            graphPHOSPi0InvYieldSys->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphEMCALPi0InvYieldSys, markerStyleDet[2] ,markerSizeDet[2]/2, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
            graphEMCALPi0InvYieldSys->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0InvYieldSys, markerStyleDet[4] ,markerSizeDet[4]/2, colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
            graphPCMEMCALPi0InvYieldSys->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphDalitzPi0InvYieldSys, markerStyleDet[5] ,markerSizeDet[5]/2, colorDet[5], colorDet[5], widthLinesBoxes, kTRUE);
            graphDalitzPi0InvYieldSys->Draw("pEsame");

            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStatUnshi, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombPi0InvYieldStatUnshi->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombPi0InvYieldStat->Draw("pEsame");

            fitTCMInvYieldPi0->SetLineColor(kBlue+2);
            fitTCMInvYieldPi0->Draw("same");
        
        canvasDummy2->Update();
        canvasDummy2->Print(Form("%s/ComparisonShiftedPi0_2760GeV.%s",outputDir.Data(),suffix.Data()));
        delete canvasDummy2;
    }   

    // *************************************************************************************************************
    // redo fitting after binshifts
    // *************************************************************************************************************
    // Tsallis function    
//     graphCombPi0InvYieldTot->Fit(fitInvYieldPi0,"QNRMEX0+","",0.3,20.);
//     fitInvYieldPi0           = FitObject("l","fitInvYieldPi0pPb5023GeV","Pi0",graphCombPi0InvYieldTot,0.3,20.,paramGraphPi0,"QNRMEX0+");
//     fitInvYieldPi0           = FitObject("l","fitInvYieldPi0pPb5023GeV","Pi0",graphCombPi0InvYieldTot,0.3,20. ,paramGraphPi0,"QNRMEX0+");
    cout << WriteParameterToFile(fitInvYieldPi0)<< endl;    
    fileFitsOutput <<  WriteParameterToFile(fitInvYieldPi0)<< endl;    
    //Two component model from Bylinkin
    fitTCMInvYieldPi0        = FitObject("tcm","fitTCMInvYieldPi0pPb5023GeV","Pi0",graphCombPi0InvYieldTot,0.3,20. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitTCMInvYieldPi0)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitTCMInvYieldPi0)<< endl;    
//     TF1* fitTCMInvYieldPi02  = FitObject("tcm","fitTCMInvYieldPi0pPb5023GeV2","Pi0",graphCombPi0InvYieldTot,0.4,40.,paramTCMPi0New,"QNRMEX0+","", kFALSE);
//     cout << WriteParameterToFile(fitTCMInvYieldPi02)<< endl;

    TF1* fitPowInvYieldPi0   = FitObject("m","fitPowInvYieldPi0pPb5023GeV","Pi0",graphCombPi0InvYieldTot,5,40. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvYieldPi0)<< endl;
                
    // *************************************************************************************************************
    // Shift graphs in Y direction as well if desired
    // *************************************************************************************************************
    TGraphAsymmErrors* graphCombPi0InvYieldTot_yShifted         = NULL;
    TGraphAsymmErrors* graphCombPi0InvYieldStat_yShifted        = NULL;
    TGraphAsymmErrors* graphCombPi0InvYieldSys_yShifted         = NULL;
    TGraphAsymmErrors* graphPCMPi0InvYieldStat_yShifted         = NULL;
    TGraphAsymmErrors* graphPCMPi0InvYieldSys_yShifted          = NULL;
    TGraphAsymmErrors* graphPHOSPi0InvYieldStat_yShifted        = NULL;
    TGraphAsymmErrors* graphPHOSPi0InvYieldSys_yShifted         = NULL;
    TGraphAsymmErrors* graphEMCALPi0InvYieldStat_yShifted       = NULL;
    TGraphAsymmErrors* graphEMCALPi0InvYieldSys_yShifted        = NULL;
    TGraphAsymmErrors* graphPCMEMCALPi0InvYieldStat_yShifted    = NULL;
    TGraphAsymmErrors* graphPCMEMCALPi0InvYieldSys_yShifted     = NULL;
    TGraphAsymmErrors* graphDalitzPi0InvYieldStat_yShifted      = NULL;
    TGraphAsymmErrors* graphDalitzPi0InvYieldSys_yShifted       = NULL;
        
    if(bWCorrection.Contains("Y") ){
        graphCombPi0InvYieldTot_yShifted        = (TGraphAsymmErrors*)graphCombPi0InvYieldTotUnshi->Clone("Pi0YShiftedCombTot");
        graphCombPi0InvYieldTot_yShifted        =  ApplyYshiftIndividualSpectra( graphCombPi0InvYieldTot_yShifted, fitInvYieldPi0);
        graphCombPi0InvYieldStat_yShifted       = (TGraphAsymmErrors*)graphCombPi0InvYieldStatUnshi->Clone("Pi0YShiftedCombStat");
        graphCombPi0InvYieldStat_yShifted       =  ApplyYshiftIndividualSpectra( graphCombPi0InvYieldStat_yShifted, fitInvYieldPi0);
        graphCombPi0InvYieldSys_yShifted        = (TGraphAsymmErrors*)graphCombPi0InvYieldSysUnshi->Clone("Pi0YShiftedCombSys");
        graphCombPi0InvYieldSys_yShifted        =  ApplyYshiftIndividualSpectra( graphCombPi0InvYieldSys_yShifted, fitInvYieldPi0);
        
        graphPCMPi0InvYieldStat_yShifted        = (TGraphAsymmErrors*)graphPCMPi0InvYieldStatUnshi->Clone("Pi0YShiftedPCMStat");
        graphPCMPi0InvYieldStat_yShifted        =  ApplyYshiftIndividualSpectra( graphPCMPi0InvYieldStat_yShifted, fitInvYieldPi0);
        graphPCMPi0InvYieldSys_yShifted         = (TGraphAsymmErrors*)graphPCMPi0InvYieldSysUnshi->Clone("Pi0YShiftedPCMSys");
        graphPCMPi0InvYieldSys_yShifted         =  ApplyYshiftIndividualSpectra( graphPCMPi0InvYieldSys_yShifted, fitInvYieldPi0);

        graphPHOSPi0InvYieldStat_yShifted       = (TGraphAsymmErrors*)graphPHOSPi0InvYieldStatUnshi->Clone("Pi0YShiftedPHOSStat");
        graphPHOSPi0InvYieldStat_yShifted       =  ApplyYshiftIndividualSpectra( graphPHOSPi0InvYieldStat_yShifted, fitInvYieldPi0);
        graphPHOSPi0InvYieldSys_yShifted        = (TGraphAsymmErrors*)graphPHOSPi0InvYieldSysUnshi->Clone("Pi0YShiftedPHOSSys");
        graphPHOSPi0InvYieldSys_yShifted        =  ApplyYshiftIndividualSpectra( graphPHOSPi0InvYieldSys_yShifted, fitInvYieldPi0);
        
        graphEMCALPi0InvYieldStat_yShifted      = (TGraphAsymmErrors*)graphEMCALPi0InvYieldStatUnshi->Clone("Pi0YShiftedEMCStat");
        graphEMCALPi0InvYieldStat_yShifted      =  ApplyYshiftIndividualSpectra( graphEMCALPi0InvYieldStat_yShifted, fitInvYieldPi0);
        graphEMCALPi0InvYieldSys_yShifted       = (TGraphAsymmErrors*)graphEMCALPi0InvYieldSysUnshi->Clone("Pi0YShiftedEMCSys");
        graphEMCALPi0InvYieldSys_yShifted       =  ApplyYshiftIndividualSpectra( graphEMCALPi0InvYieldSys_yShifted, fitInvYieldPi0);
        
        graphPCMEMCALPi0InvYieldStat_yShifted   = (TGraphAsymmErrors*)graphPCMEMCALPi0InvYieldStatUnshi->Clone("Pi0YShiftedPCMEMCStat");
        graphPCMEMCALPi0InvYieldStat_yShifted   =  ApplyYshiftIndividualSpectra( graphPCMEMCALPi0InvYieldStat_yShifted, fitInvYieldPi0);
        graphPCMEMCALPi0InvYieldSys_yShifted    = (TGraphAsymmErrors*)graphPCMEMCALPi0InvYieldSysUnshi->Clone("Pi0YShiftedPCMEMCStat");
        graphPCMEMCALPi0InvYieldSys_yShifted    =  ApplyYshiftIndividualSpectra( graphPCMEMCALPi0InvYieldSys_yShifted, fitInvYieldPi0);

        graphDalitzPi0InvYieldStat_yShifted     = (TGraphAsymmErrors*)graphDalitzPi0InvYieldStatUnshi->Clone("Pi0YShiftedPCMEMCStat");
        graphDalitzPi0InvYieldStat_yShifted     =  ApplyYshiftIndividualSpectra( graphDalitzPi0InvYieldStat_yShifted, fitInvYieldPi0);
        graphDalitzPi0InvYieldSys_yShifted      = (TGraphAsymmErrors*)graphDalitzPi0InvYieldSysUnshi->Clone("Pi0YShiftedPCMEMCStat");
        graphDalitzPi0InvYieldSys_yShifted      =  ApplyYshiftIndividualSpectra( graphDalitzPi0InvYieldSys_yShifted, fitInvYieldPi0);
    }

    // *************************************************************************************************************
    // Calculate ratios to combined fit
    // *************************************************************************************************************
//     TH1D* histoRatioPi0Pythia8ToFit                     = (TH1D*) histoPythia8InvYieldPi0->Clone(); 
//     histoRatioPi0Pythia8ToFit                           = CalculateHistoRatioToFit (histoRatioPi0Pythia8ToFit, fitTCMInvYieldPi0); 
//     histoRatioPi0Pythia8ToFit->GetXaxis()->SetRangeUser(0.4,40);
//     
//     TGraph* graphRatioPi0CombNLOMuHalf                  = (TGraph*)graphNLOCalcPi0MuHalf->Clone();
//     TGraph* graphRatioPi0CombNLOMuOne                   = (TGraph*)graphNLOCalcPi0MuOne->Clone();
//     TGraph* graphRatioPi0CombNLOMuTwo                   = (TGraph*)graphNLOCalcPi0MuTwo->Clone();
//     TGraphAsymmErrors* graphRatioPi0CombNLODSS14        = (TGraphAsymmErrors*)graphNLODSS14Calc->Clone();
//     graphRatioPi0CombNLOMuHalf                          = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuHalf, fitTCMInvYieldPi0); 
//     graphRatioPi0CombNLOMuOne                           = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuOne, fitTCMInvYieldPi0); 
//     graphRatioPi0CombNLOMuTwo                           = CalculateGraphRatioToFit (graphRatioPi0CombNLOMuTwo, fitTCMInvYieldPi0); 
//     graphRatioPi0CombNLODSS14                           = CalculateGraphErrRatioToFit(graphRatioPi0CombNLODSS14, fitTCMInvYieldPi0); 
// 
//     TGraph* graphRatioPi0CombCGC                        = (TGraph*)graphCGCPi0->Clone();
//     graphRatioPi0CombCGC                                = CalculateGraphRatioToFit (graphRatioPi0CombCGC, fitTCMInvYieldPi0); 
    
    TGraphAsymmErrors* graphRatioPi0CombCombFitTot      = (TGraphAsymmErrors*)graphCombPi0InvYieldTot->Clone();
    graphRatioPi0CombCombFitTot                         = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitTot, fitTCMInvYieldPi0); 
    TGraphAsymmErrors* graphRatioPi0CombCombFitStat     = (TGraphAsymmErrors*)graphCombPi0InvYieldStat->Clone();
    graphRatioPi0CombCombFitStat                        = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitStat, fitTCMInvYieldPi0); 
    TGraphAsymmErrors* graphRatioPi0CombCombFitSys      = (TGraphAsymmErrors*)graphCombPi0InvYieldSys->Clone();
    graphRatioPi0CombCombFitSys                         = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitSys, fitTCMInvYieldPi0); 
    TGraphAsymmErrors* graphRatioPi0PCMCombFitStat      = (TGraphAsymmErrors*)graphPCMPi0InvYieldStat->Clone();
    graphRatioPi0PCMCombFitStat                         = CalculateGraphErrRatioToFit(graphRatioPi0PCMCombFitStat, fitTCMInvYieldPi0); 
    TGraphAsymmErrors* graphRatioPi0PCMCombFitSys       = (TGraphAsymmErrors*)graphPCMPi0InvYieldSys->Clone();
    graphRatioPi0PCMCombFitSys                          = CalculateGraphErrRatioToFit(graphRatioPi0PCMCombFitSys, fitTCMInvYieldPi0); 
    TGraphAsymmErrors* graphRatioPi0PHOSCombFitStat     = (TGraphAsymmErrors*)graphPHOSPi0InvYieldStat->Clone();
    graphRatioPi0PHOSCombFitStat                        = CalculateGraphErrRatioToFit(graphRatioPi0PHOSCombFitStat, fitTCMInvYieldPi0); 
    TGraphAsymmErrors* graphRatioPi0PHOSCombFitSys      = (TGraphAsymmErrors*)graphPHOSPi0InvYieldSys->Clone();
    graphRatioPi0PHOSCombFitSys                         = CalculateGraphErrRatioToFit(graphRatioPi0PHOSCombFitSys, fitTCMInvYieldPi0); 
    TGraphAsymmErrors* graphRatioPi0EMCALCombFitStat    = (TGraphAsymmErrors*)graphEMCALPi0InvYieldStat->Clone();
    graphRatioPi0EMCALCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioPi0EMCALCombFitStat, fitTCMInvYieldPi0); 
    TGraphAsymmErrors* graphRatioPi0EMCALCombFitSys     = (TGraphAsymmErrors*)graphEMCALPi0InvYieldSys->Clone();
    graphRatioPi0EMCALCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioPi0EMCALCombFitSys, fitTCMInvYieldPi0); 
    TGraphAsymmErrors* graphRatioPi0PCMEMCALCombFitStat = (TGraphAsymmErrors*)graphPCMEMCALPi0InvYieldStat->Clone();
    graphRatioPi0PCMEMCALCombFitStat                    = CalculateGraphErrRatioToFit(graphRatioPi0PCMEMCALCombFitStat, fitTCMInvYieldPi0); 
    TGraphAsymmErrors* graphRatioPi0PCMEMCALCombFitSys  = (TGraphAsymmErrors*)graphPCMEMCALPi0InvYieldSys->Clone();
    graphRatioPi0PCMEMCALCombFitSys                     = CalculateGraphErrRatioToFit(graphRatioPi0PCMEMCALCombFitSys, fitTCMInvYieldPi0); 
    TGraphAsymmErrors* graphRatioPi0DalitzCombFitStat   = (TGraphAsymmErrors*)graphDalitzPi0InvYieldStat->Clone();
    graphRatioPi0DalitzCombFitStat                      = CalculateGraphErrRatioToFit(graphRatioPi0DalitzCombFitStat, fitTCMInvYieldPi0); 
    TGraphAsymmErrors* graphRatioPi0DalitzCombFitSys    = (TGraphAsymmErrors*)graphDalitzPi0InvYieldSys->Clone();
    graphRatioPi0DalitzCombFitSys                       = CalculateGraphErrRatioToFit(graphRatioPi0DalitzCombFitSys, fitTCMInvYieldPi0); 

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
    histo2DPi0RatioToCombFit               = new TH2F("histo2DPi0RatioToCombFit","histo2DPi0RatioToCombFit",1000,0.23,31.,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DPi0RatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
    histo2DPi0RatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DPi0RatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
    histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.45,2.75);
    histo2DPi0RatioToCombFit->Draw("copy");

        ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStat);
 
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioPi0CombCombFitSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStat, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioPi0CombCombFitStat->Draw("p,same,z");

        DrawGammaLines(0.23, 70. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.23, 70. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.9, 0.9,0.1, kGray, 7);

        TLatex *labelRatioToFitEnergy   = new TLatex(0.95, 0.92, collisionSystempPb.Data());
        SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitEnergy->Draw();
        TLatex *labelRatioToFitALICE    = new TLatex(0.95, 0.87, "ALICE");
        SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitALICE->Draw();
        TLatex *labelRatioToFitPi0      = new TLatex(0.95, 0.817, "#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitPi0->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToCombFit_PPb5023GeV.%s",outputDir.Data(),suffix.Data()));
    
    // **********************************************************************************************************************
    // *******************************************Plot Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************
    
    canvasRatioToCombFit->cd();
    histo2DPi0RatioToCombFit->Draw("copy");
    
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0PCMCombFitStat);
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0PHOSCombFitStat);
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0EMCALCombFitStat);
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0PCMEMCALCombFitStat);
        
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitSys, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMCombFitStat, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0DalitzCombFitSys, markerStyleDet[5] ,markerSizeDet[5]*0.5, colorDet[5], colorDet[5], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0DalitzCombFitStat, markerStyleDet[5] ,markerSizeDet[5]*0.5, colorDet[5], colorDet[5]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PHOSCombFitSys, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PHOSCombFitStat, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALCombFitSys, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0EMCALCombFitStat, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMEMCALCombFitSys, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0PCMEMCALCombFitStat, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4]);
        
        graphRatioPi0PCMCombFitSys->Draw("E2same");
        graphRatioPi0DalitzCombFitSys->Draw("E2same");
        graphRatioPi0PHOSCombFitSys->Draw("E2same");
        graphRatioPi0EMCALCombFitSys->Draw("E2same");
        graphRatioPi0PCMEMCALCombFitSys->Draw("E2same");
        
        graphRatioPi0PCMCombFitStat->Draw("p,same,z");
        graphRatioPi0DalitzCombFitStat->Draw("p,same,z");
        graphRatioPi0PHOSCombFitStat->Draw("p,same,z");
        graphRatioPi0EMCALCombFitStat->Draw("p,same,z");
        graphRatioPi0PCMEMCALCombFitStat->Draw("p,same,z");

        DrawGammaLines(0.23, 70. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.23, 70. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.9, 0.9,0.5, kGray, 7);
        
        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitPi0->Draw();
    
        //****************************** Definition of the Legend ******************************************
        //**************** Row def ************************
        Double_t rowsLegendOnlyPi0Ratio[6]          = {0.92, 0.867, 0.807, 0.747, 0.687, 0.627};
        Double_t rowsLegendOnlyPi0RatioAbs[6]       = {0.92, 2.48, 2.32, 2.16, 2.00, 1.84 };
//         Double_t rowsLegendOnlyPi0RatioAbs[6]       = {0.92, 2.18, 2.03, 1.87, 1.73, 1.56 };
        Double_t columnsLegendOnlyPi0Ratio[3]       = {0.115, 0.325, 0.41};
        Double_t columnsLegendOnlyPi0RatioAbs[3]    = {0.115, 1., 1.43};
        Double_t lengthBox                          = 0.1;
        Double_t heightBox                          = 0.08/2;
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
        TLatex *textDalitzOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[5],nameMeasGlobalLabel[5]);
        SetStyleTLatex( textDalitzOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
        textDalitzOnlyRatioPi0->Draw();
        
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
        TMarker* markerDalitzPi0OnlyRatioPi0        = CreateMarkerFromGraph(graphRatioPi0DalitzCombFitSys, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[5],1);
        markerDalitzPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[5]);

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
        TBox* boxDalitzPi0OnlyRatioPi0              = CreateBoxFromGraph(graphRatioPi0DalitzCombFitSys, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[5]- heightBox,
                                                                         columnsLegendOnlyPi0RatioAbs[2]+ 5*lengthBox, rowsLegendOnlyPi0RatioAbs[5]+ heightBox);
        boxDalitzPi0OnlyRatioPi0->Draw("l");
        
    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFit_PP.%s",outputDir.Data(),suffix.Data()));
    
    
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
    statErrorCollectionEta[0]       = (TH1D*)histoPCMEtaInvYieldStat->Clone("statErrPCMEta");
    statErrorCollectionEta[2]       = (TH1D*)histoEMCALEtaInvYieldStat->Clone("statErrEMCALEta");
//     statErrorCollectionEta[4]       = (TH1D*)histoPCMEMCALEtaInvYieldStat->Clone("statErrPCMEMCALEta");
   
    // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)    
    // for systematic error from respective method
    TGraphAsymmErrors* sysErrorCollectionEta[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorCollectionEta[i]    = NULL;
    }    
    sysErrorCollectionEta[0]        = (TGraphAsymmErrors*)graphPCMEtaInvYieldSys->Clone("sysErrPCMEta");
    sysErrorCollectionEta[2]        = (TGraphAsymmErrors*)graphEMCALEtaInvYieldSys->Clone("sysErrEMCALEta");
//     sysErrorCollectionEta[4]        = (TGraphAsymmErrors*)graphPCMEMCALEtaInvYieldSys->Clone("sysErrPCMEMCALEta");

    TGraphAsymmErrors* statErrorRelCollectionEta[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorRelCollectionEta[i]        = NULL;
    }    
    for (Int_t i = 0; i < 11; i++){
        if (statErrorCollectionEta[i]){
            statErrorRelCollectionEta[i]    = new TGraphAsymmErrors(statErrorCollectionEta[i]);
            while (statErrorRelCollectionEta[i]->GetY()[0] == 0) 
                statErrorRelCollectionEta[i]->RemovePoint(0);
            statErrorRelCollectionEta[i]    = CalculateRelErrUpAsymmGraph( statErrorRelCollectionEta[i], Form("relativeStatErrorEta_%s", nameMeasGlobal[i].Data()));
        }   
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
    Double_t xPtLimitsEta[50];
    Int_t maxNBinsEta       = GetBinning( xPtLimitsEta, "Eta", "pPb_5.023TeV", 20 );
    for (Int_t i = 0; i< maxNBinsEta; i++){
        cout << i << ": "<< xPtLimitsEta[i] <<" - " << xPtLimitsEta[i+1]<< ", " <<endl;
    }
    cout << endl;
    cout << "Nbins: eta "<< maxNBinsEta << endl;

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Eta                    
    Int_t offSetsEta[11]            = { 0,  0,  0,  0,  0,
                                        0,  0,  0,  0,  0, 
                                        0};
    Int_t offSetsEtaSys[11]         = { 3,  0,  8,  0,  5,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t offSetEtaShifting[11]     = { 0,  0,  5,  0,  2,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t nComBinsEtaShifting[11]   = { 13,  0,  16,  0,  15,
                                        0,  0,  0,  0,  0,
                                        0};                                        
    
    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************
                                
                                        
    // Combination w/o PCM-EMC
    TGraph* graphWeightsEtaWOPE[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsEtaWOPE[i]         = NULL;
    }    
                                
    // Declaration & calculation of combined spectrum                            
    TString fileNameEtaOutputWeightingWOPE              = Form("%s/Eta_WeightingMethodWOPE.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombEtaInvYieldWOPEStat     = NULL;
    TGraphAsymmErrors* graphCombEtaInvYieldWOPESys      = NULL;
    TGraphAsymmErrors* graphCombEtaInvYieldWOPETot      = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionEta,    sysErrorCollectionEta,     
                                                                                                xPtLimitsEta, maxNBinsEta,
                                                                                                offSetsEta, offSetsEtaSys,
                                                                                                graphCombEtaInvYieldWOPEStat, graphCombEtaInvYieldWOPESys,
                                                                                                fileNameEtaOutputWeightingWOPE,"pPb_5.023TeV", "Eta", kFALSE, 
                                                                                                NULL, fileNameCorrFactors
                                                                                            );
    if (graphCombEtaInvYieldWOPEStat == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }    
    graphCombEtaInvYieldWOPEStat->Print();
    while (graphCombEtaInvYieldWOPETot->GetX()[0] < 0.7){
        graphCombEtaInvYieldWOPETot->RemovePoint(0);
        cout << "removed first point from tot-graph" << endl; 
    }    
    while (graphCombEtaInvYieldWOPEStat->GetX()[0] < 0.7){
        graphCombEtaInvYieldWOPEStat->RemovePoint(0);
        cout << "removed first point from stat-graph" << endl; 
    }    
    while (graphCombEtaInvYieldWOPESys->GetX()[0] < 0.7){
        graphCombEtaInvYieldWOPESys->RemovePoint(0);
        cout << "removed first point from sys-graph" << endl; 
    }
    
    // Reading weights from output file for plotting
    ifstream fileWeightsEtaReadWOPE;
    fileWeightsEtaReadWOPE.open(fileNameEtaOutputWeightingWOPE,ios_base::in);
    cout << "reading" << fileNameEtaOutputWeightingWOPE << endl;
    Double_t xValuesEtaReadWOPE[50];
    Double_t weightsEtaReadWOPE[11][50];
    Int_t availableEtaMeasWOPE[11]      = { -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1,
                                            -1};
    Int_t nMeasSetEtaWOPE               = 2;
    Int_t nPtBinsEtaReadWOPE            = 0;
    while(!fileWeightsEtaReadWOPE.eof() && nPtBinsEtaReadWOPE < 50){
        TString garbage             = "";
        if (nPtBinsEtaReadWOPE == 0){
            fileWeightsEtaReadWOPE >> garbage ;//>> availableEtaMeasWOPE[0] >> availableEtaMeasWOPE[1] >> availableEtaMeasWOPE[2] >> availableEtaMeasWOPE[3];
            for (Int_t i = 0; i < nMeasSetEtaWOPE; i++){
                fileWeightsEtaReadWOPE >> availableEtaMeasWOPE[i] ;
            }    
            cout << "read following measurements: "; 
            for (Int_t i = 0; i < 11; i++){
                cout << availableEtaMeasWOPE[i] << "\t" ;
            }    
            cout << endl;
        } else {
            fileWeightsEtaReadWOPE >> xValuesEtaReadWOPE[nPtBinsEtaReadWOPE-1];
            for (Int_t i = 0; i < nMeasSetEtaWOPE; i++){
                fileWeightsEtaReadWOPE >> weightsEtaReadWOPE[availableEtaMeasWOPE[i]][nPtBinsEtaReadWOPE-1] ;
            }    
            cout << "read: "<<  nPtBinsEtaReadWOPE << "\t"<< xValuesEtaReadWOPE[nPtBinsEtaReadWOPE-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetEtaWOPE; i++){
                cout << weightsEtaReadWOPE[availableEtaMeasWOPE[i]][nPtBinsEtaReadWOPE-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsEtaReadWOPE++;
    }
    nPtBinsEtaReadWOPE                    = nPtBinsEtaReadWOPE-2 ;
    fileWeightsEtaReadWOPE.close();

    for (Int_t i = 0; i < nMeasSetEtaWOPE; i++){
        graphWeightsEtaWOPE[availableEtaMeasWOPE[i]]                        = new TGraph(nPtBinsEtaReadWOPE,xValuesEtaReadWOPE,weightsEtaReadWOPE[availableEtaMeasWOPE[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsEtaReadWOPE; n++){
            if (graphWeightsEtaWOPE[availableEtaMeasWOPE[i]]->GetY()[bin] == 0) graphWeightsEtaWOPE[availableEtaMeasWOPE[i]]->RemovePoint(bin);
            else bin++;
        }    
    }    

    // **********************************************************************************************************************
    // ************************ Adding PCM-EMC measurements to matrix & calculating relativ errors ************************
    // **********************************************************************************************************************
    
    statErrorCollectionEta[4]               = (TH1D*)histoPCMEMCALEtaInvYieldStat->Clone("statErrPCMEMCALEta");
    sysErrorCollectionEta[4]                = (TGraphAsymmErrors*)graphPCMEMCALEtaInvYieldSys->Clone("sysErrPCMEMCALEta");

    if (statErrorCollectionEta[4]){
        statErrorRelCollectionEta[4]        = new TGraphAsymmErrors(statErrorCollectionEta[4]);
        while (statErrorRelCollectionEta[4]->GetY()[0] == 0) 
            statErrorRelCollectionEta[4]->RemovePoint(0);
        statErrorRelCollectionEta[4]        = CalculateRelErrUpAsymmGraph( statErrorRelCollectionEta[4], Form("relativeStatErrorEta_%s", nameMeasGlobal[4].Data()));
    }        
    if (sysErrorCollectionEta[4]) 
        sysErrorRelCollectionEta[4]         = CalculateRelErrUpAsymmGraph( sysErrorCollectionEta[4], Form("relativeSysErrorEta_%s", nameMeasGlobal[4].Data()));

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation incl PCM-EMC ********************************
    // **********************************************************************************************************************
                                
    TGraph* graphWeightsEta[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsEta[i]         = NULL;
    }    
                                
    // Declaration & calculation of combined spectrum                            
    TString fileNameEtaOutputWeighting              = Form("%s/Eta_WeightingMethod.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombEtaInvYieldStat     = NULL;
    TGraphAsymmErrors* graphCombEtaInvYieldSys      = NULL;
    TGraphAsymmErrors* graphCombEtaInvYieldTot      = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionEta,    sysErrorCollectionEta,     
                                                                                                xPtLimitsEta, maxNBinsEta,
                                                                                                offSetsEta, offSetsEtaSys,
                                                                                                graphCombEtaInvYieldStat, graphCombEtaInvYieldSys,
                                                                                                fileNameEtaOutputWeighting,"pPb_5.023TeV", "Eta", kFALSE, 
                                                                                                NULL, fileNameCorrFactors
                                                                                            );
    graphCombEtaInvYieldStat->Print();
    if (graphCombEtaInvYieldStat == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }    
    while (graphCombEtaInvYieldTot->GetX()[0] < 0.7){
        graphCombEtaInvYieldTot->RemovePoint(0);
        cout << "removed first point from tot-graph" << endl; 
    }    
    while (graphCombEtaInvYieldStat->GetX()[0] < 0.7){
        graphCombEtaInvYieldStat->RemovePoint(0);
        cout << "removed first point from stat-graph" << endl; 
    }    
    while (graphCombEtaInvYieldSys->GetX()[0] < 0.7){
        graphCombEtaInvYieldSys->RemovePoint(0);
        cout << "removed first point from sys-graph" << endl; 
    }
    
    // Reading weights from output file for plotting
    ifstream fileWeightsEtaRead;
    fileWeightsEtaRead.open(fileNameEtaOutputWeighting,ios_base::in);
    cout << "reading" << fileNameEtaOutputWeighting << endl;
    Double_t xValuesEtaRead[50];
    Double_t weightsEtaRead[11][50];
    Int_t availableEtaMeas[11]      = { -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1,
                                            -1};
    Int_t nMeasSetEta               = 3;
    Int_t nPtBinsEtaRead            = 0;
    while(!fileWeightsEtaRead.eof() && nPtBinsEtaRead < 50){
        TString garbage             = "";
        if (nPtBinsEtaRead == 0){
            fileWeightsEtaRead >> garbage ;//>> availableEtaMeas[0] >> availableEtaMeas[1] >> availableEtaMeas[2] >> availableEtaMeas[3];
            for (Int_t i = 0; i < nMeasSetEta; i++){
                fileWeightsEtaRead >> availableEtaMeas[i] ;
            }    
            cout << "read following measurements: "; 
            for (Int_t i = 0; i < 11; i++){
                cout << availableEtaMeas[i] << "\t" ;
            }    
            cout << endl;
        } else {
            fileWeightsEtaRead >> xValuesEtaRead[nPtBinsEtaRead-1];
            for (Int_t i = 0; i < nMeasSetEta; i++){
                fileWeightsEtaRead >> weightsEtaRead[availableEtaMeas[i]][nPtBinsEtaRead-1] ;
            }    
            cout << "read: "<<  nPtBinsEtaRead << "\t"<< xValuesEtaRead[nPtBinsEtaRead-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetEta; i++){
                cout << weightsEtaRead[availableEtaMeas[i]][nPtBinsEtaRead-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsEtaRead++;
    }
    nPtBinsEtaRead                    = nPtBinsEtaRead-2 ;
    fileWeightsEtaRead.close();

    for (Int_t i = 0; i < nMeasSetEta; i++){
        graphWeightsEta[availableEtaMeas[i]]                        = new TGraph(nPtBinsEtaRead,xValuesEtaRead,weightsEtaRead[availableEtaMeas[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsEtaRead; n++){
            if (graphWeightsEta[availableEtaMeas[i]]->GetY()[bin] == 0) graphWeightsEta[availableEtaMeas[i]]->RemovePoint(bin);
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

        TLegend* legendWeightsEtaWOPE   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEtaWOPE), 32);
        for (Int_t i = 0; i < nMeasSetEtaWOPE; i++){
            DrawGammaSetMarkerTGraph(graphWeightsEtaWOPE[availableEtaMeasWOPE[i]], markerStyleDet[availableEtaMeasWOPE[i]], markerSizeDet[availableEtaMeasWOPE[i]]*0.5, colorDet[availableEtaMeasWOPE[i]] , colorDet[availableEtaMeasWOPE[i]]);
            graphWeightsEtaWOPE[availableEtaMeasWOPE[i]]->Draw("p,same,z");
            legendWeightsEtaWOPE->AddEntry(graphWeightsEtaWOPE[availableEtaMeasWOPE[i]],nameMeasGlobalLabel[availableEtaMeasWOPE[i]],"p");
        }    
        legendWeightsEtaWOPE->Draw();

        labelWeightsEnergy->Draw();
        TLatex *labelWeightsEta         = new TLatex(0.7,0.16,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelWeightsEta, 0.85*textSizeLabelsPixel,4, 1, 43);
        labelWeightsEta->Draw();

//      DrawGammaLines(0.33, 25. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.33, 25. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.33, 25. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.33, 25. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.33, 25. , 0.2, 0.2,0.1, kGray, 3);
        
    canvasWeights->SaveAs(Form("%s/Eta_Weights_WOPCMEMC.%s",outputDir.Data(),suffix.Data()));
    
    histo2DEtaWeights->Draw("copy");
        TLegend* legendWeightsEta   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEta), 32);
        for (Int_t i = 0; i < nMeasSetEta; i++){
            DrawGammaSetMarkerTGraph(graphWeightsEta[availableEtaMeas[i]], markerStyleDet[availableEtaMeas[i]], markerSizeDet[availableEtaMeas[i]]*0.5, colorDet[availableEtaMeas[i]] , colorDet[availableEtaMeas[i]]);
            graphWeightsEta[availableEtaMeas[i]]->Draw("p,same,z");
            legendWeightsEta->AddEntry(graphWeightsEta[availableEtaMeas[i]],nameMeasGlobalLabel[availableEtaMeas[i]],"p");
        }    
        legendWeightsEta->Draw();

        labelWeightsEnergy->Draw();
        labelWeightsEta->Draw();

//      DrawGammaLines(0.33, 25. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.33, 25. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.33, 25. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.33, 25. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.33, 25. , 0.2, 0.2,0.1, kGray, 3);
        
    canvasWeights->SaveAs(Form("%s/Eta_Weights.%s",outputDir.Data(),suffix.Data()));

    
//     return;
    // *********************************************************************************************************************
    // ************************************ Visualize relative errors Eta ******************************************************
    // *********************************************************************************************************************
    
    canvasRelSysErr->cd();
   
    TH2F * histo2DRelSysErrEta;
    histo2DRelSysErrEta                 = new TH2F("histo2DRelSysErrEta","histo2DRelSysErrEta",11000,0.33,25.,1000,0,60.5);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErrEta, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErrEta->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErrEta->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelSysErrEta->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelSysErrEta->Draw("copy");

        TLegend* legendRelSysErrEtaWOPE        = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetEtaWOPE), 0.95, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetEtaWOPE; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionEta[availableEtaMeasWOPE[i]], markerStyleDet[availableEtaMeasWOPE[i]], markerSizeDet[availableEtaMeasWOPE[i]]*0.5, colorDet[availableEtaMeasWOPE[i]],
                                     colorDet[availableEtaMeasWOPE[i]]);
            sysErrorRelCollectionEta[availableEtaMeasWOPE[i]]->Draw("p,same,z");
            legendRelSysErrEtaWOPE->AddEntry(sysErrorRelCollectionEta[availableEtaMeasWOPE[i]],nameMeasGlobalLabel[availableEtaMeasWOPE[i]],"p");
        }    
        legendRelSysErrEtaWOPE->Draw();

        labelRelSysErrEnergy->Draw();
        TLatex *labelRelSysErrEta       = new TLatex(0.15,0.85,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelSysErrEta, 0.85*textSizeLabelsPixel,4, 1, 43);
        labelRelSysErrEta->Draw();
        
    canvasRelSysErr->SaveAs(Form("%s/Eta_RelSysErr_WOPE.%s",outputDir.Data(),suffix.Data()));
    histo2DRelSysErrEta->Draw("copy");

        TLegend* legendRelSysErrEta        = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetEta), 0.95, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetEta; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionEta[availableEtaMeas[i]], markerStyleDet[availableEtaMeas[i]], markerSizeDet[availableEtaMeas[i]]*0.5, colorDet[availableEtaMeas[i]],
                                     colorDet[availableEtaMeas[i]]);
            sysErrorRelCollectionEta[availableEtaMeas[i]]->Draw("p,same,z");
            legendRelSysErrEta->AddEntry(sysErrorRelCollectionEta[availableEtaMeas[i]],nameMeasGlobalLabel[availableEtaMeas[i]],"p");
        }    
        legendRelSysErrEta->Draw();

        labelRelSysErrEnergy->Draw();
        labelRelSysErrEta->Draw();
        
    canvasRelSysErr->SaveAs(Form("%s/Eta_RelSysErr.%s",outputDir.Data(),suffix.Data()));
        
    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors Eta **************************************************
    //  *********************************************************************************************************************
    
    canvasRelStatErr->cd();
   
    TH2F * histo2DRelStatErrEta;
    histo2DRelStatErrEta                = new TH2F("histo2DRelStatErrEta","histo2DRelStatErrEta",11000,0.33,25.,1000,0,60.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErrEta, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErrEta->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErrEta->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelStatErrEta->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelStatErrEta->Draw("copy");
        
        TLegend* legendRelStatErrEtaWOPE       = GetAndSetLegend2(0.14, 0.92-(0.035*nMeasSetEtaWOPE), 0.45, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetEtaWOPE; i++){
            DrawGammaSetMarkerTGraph(statErrorRelCollectionEta[availableEtaMeasWOPE[i]], markerStyleDet[availableEtaMeasWOPE[i]], markerSizeDet[availableEtaMeasWOPE[i]]*0.5, colorDet[availableEtaMeasWOPE[i]] , 
                               colorDet[availableEtaMeasWOPE[i]]);
            statErrorRelCollectionEta[availableEtaMeasWOPE[i]]->Draw("p,same,z");
            legendRelStatErrEtaWOPE->AddEntry(statErrorRelCollectionEta[availableEtaMeasWOPE[i]],nameMeasGlobalLabel[availableEtaMeasWOPE[i]],"p");
        }    
        legendRelStatErrEtaWOPE->Draw();

        labelRelStatErrEnergy->Draw();
        TLatex *labelRelStatErrEta      = new TLatex(0.75,0.85,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelStatErrEta, 0.85*textSizeLabelsPixel,4, 1, 43);
        labelRelStatErrEta->Draw();
        
    canvasRelStatErr->SaveAs(Form("%s/Eta_RelStatErr_WOPE.%s",outputDir.Data(),suffix.Data()));
    histo2DRelStatErrEta->Draw("copy");
        
        TLegend* legendRelStatErrEta       = GetAndSetLegend2(0.14, 0.92-(0.035*nMeasSetEta), 0.45, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetEta; i++){
            DrawGammaSetMarkerTGraph(statErrorRelCollectionEta[availableEtaMeas[i]], markerStyleDet[availableEtaMeas[i]], markerSizeDet[availableEtaMeas[i]]*0.5, colorDet[availableEtaMeas[i]] , 
                               colorDet[availableEtaMeas[i]]);
            statErrorRelCollectionEta[availableEtaMeas[i]]->Draw("p,same,z");
            legendRelStatErrEta->AddEntry(statErrorRelCollectionEta[availableEtaMeas[i]],nameMeasGlobalLabel[availableEtaMeas[i]],"p");
        }    
        legendRelStatErrEta->Draw();

        labelRelStatErrEnergy->Draw();
        labelRelStatErrEta->Draw();
        
    canvasRelStatErr->SaveAs(Form("%s/Eta_RelStatErr.%s",outputDir.Data(),suffix.Data()));


    
    //  *********************************************************************************************************************
    //  ************************ Visualize relative total errors of different combination methods Eta ***********************
    //  *********************************************************************************************************************
    TGraphAsymmErrors* graphCombEtaInvYieldWOPERelStat   = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldWOPEStat, "relativeStatErrorEta_Method");
    TGraphAsymmErrors* graphCombEtaInvYieldRelWOPESys    = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldWOPESys, "relativeSysErrorEta_Method");
    TGraphAsymmErrors* graphCombEtaInvYieldWOPERelTot    = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldWOPETot, "relativeTotalErrorEta_Method");
    TGraphAsymmErrors* graphCombEtaInvYieldRelStat       = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldStat, "relativeStatErrorEta_Method");
    TGraphAsymmErrors* graphCombEtaInvYieldRelSys        = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldSys, "relativeSysErrorEta_Method");
    TGraphAsymmErrors* graphCombEtaInvYieldRelTot        = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldTot, "relativeTotalErrorEta_Method");
    
    
    canvasRelTotErr->cd();
    TH2F * histo2DRelTotErrEta;
    histo2DRelTotErrEta                 = new TH2F("histo2DRelTotErrEta","histo2DRelTotErrEta",11000,0.33,25.,1000,0,60.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrEta, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrEta->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrEta->GetXaxis()->SetLabelOffset(-0.01);    
    histo2DRelTotErrEta->Draw("copy");
    
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldWOPERelTot, markerStyleComb, markerSizeComb, kBlue+2 , kBlue+2);
        graphCombEtaInvYieldWOPERelTot->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelTot, markerStyleComb+4, markerSizeComb, kRed+2 , kRed+2);
        graphCombEtaInvYieldRelTot->Draw("p,same,z");

        TLegend* legendRelTotErrEta     = GetAndSetLegend2(0.14, 0.92-(0.035*2), 0.45, 0.92, 32);
        legendRelTotErrEta->AddEntry(graphCombEtaInvYieldWOPERelTot,"PCM, EMC","p");
        legendRelTotErrEta->AddEntry(graphCombEtaInvYieldRelTot,"All","p");
        legendRelTotErrEta->Draw();

        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrEta       = new TLatex(0.75,0.85,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelTotErrEta, 0.85*textSizeLabelsPixel,4, 1, 43);
        labelRelTotErrEta->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/Eta_RelTotErr.%s",outputDir.Data(),suffix.Data()));
    
    histo2DRelTotErrEta->GetYaxis()->SetTitle("Err (%)");
    histo2DRelTotErrEta->Draw("copy");
        
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldWOPERelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombEtaInvYieldWOPERelTot->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldWOPERelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombEtaInvYieldWOPERelStat->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelWOPESys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombEtaInvYieldRelWOPESys->SetLineStyle(7);
        graphCombEtaInvYieldRelWOPESys->Draw("l,x0,same,e1");

        legendRelTotErr3->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrEta->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/Eta_RelMethoddecomp_WOPE.%s",outputDir.Data(),suffix.Data()));

    histo2DRelTotErrEta->Draw("copy");
        
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombEtaInvYieldRelTot->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombEtaInvYieldRelStat->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombEtaInvYieldRelSys->SetLineStyle(7);
        graphCombEtaInvYieldRelSys->Draw("l,x0,same,e1");

        legendRelTotErr3->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrEta->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/Eta_RelMethoddecomp.%s",outputDir.Data(),suffix.Data()));
    
        

    // **********************************************************************************************************************
    // ************************************* Calculating bin shifted spectra & fitting **************************************
    // **********************************************************************************************************************
    
    // Cloning spectra
    TGraphAsymmErrors* graphCombEtaInvYieldTotUnshi      = (TGraphAsymmErrors*)graphCombEtaInvYieldTot->Clone("EtaUnshifted"); 
    TGraphAsymmErrors* graphCombEtaInvYieldStatUnshi     = (TGraphAsymmErrors*)graphCombEtaInvYieldStat->Clone("EtaUnshiftedStat"); 
    TGraphAsymmErrors* graphCombEtaInvYieldSysUnshi      = (TGraphAsymmErrors*)graphCombEtaInvYieldSys->Clone("EtaUnshiftedSys"); 

    TGraphAsymmErrors* graphPCMEtaInvYieldStatUnshi       = (TGraphAsymmErrors*)graphPCMEtaInvYieldStat->Clone("EtaUnshiftedStatPCM"); 
    TGraphAsymmErrors* graphPCMEtaInvYieldSysUnshi        = (TGraphAsymmErrors*)graphPCMEtaInvYieldSys->Clone("EtaUnshiftedSysPCM"); 
    
    TGraphAsymmErrors* graphEMCALEtaInvYieldStatUnshi     = (TGraphAsymmErrors*)graphEMCALEtaInvYieldStat->Clone("EtaUnshiftedStatEMCAL"); 
    TGraphAsymmErrors* graphEMCALEtaInvYieldSysUnshi      = (TGraphAsymmErrors*)graphEMCALEtaInvYieldSys->Clone("EtaUnshiftedSysEMCAL"); 

    TGraphAsymmErrors* graphPCMEMCALEtaInvYieldStatUnshi  = (TGraphAsymmErrors*)graphPCMEMCALEtaInvYieldStat->Clone("EtaUnshiftedStatPCMEMCAL"); 
    TGraphAsymmErrors* graphPCMEMCALEtaInvYieldSysUnshi   = (TGraphAsymmErrors*)graphPCMEMCALEtaInvYieldSys->Clone("EtaUnshiftedSysPCMEMCAL"); 
    
    // fitting spectrum with intial parameters
    Double_t paramGraphEta[3]                    = {1.0e10, 7., 0.13};
    TF1* fitInvYieldEta                         = FitObject("l","fitInvYieldEta2760GeV","Eta",graphCombEtaInvYieldTotUnshi,0.7,20.,paramGraphEta,"QNRMEX0+");
    TF1* fitInvYieldEtaGraph                    = (TF1*)fitInvYieldEta->Clone("fitInvYieldEta2760GeVGraph"); 
    
    Double_t paramTCMEta[5]                     = {graphCombEtaInvYieldTot->GetY()[0],0.5,graphCombEtaInvYieldTot->GetY()[0]/1000,0.8,3};


    // *************************************************************************************************************
    // Shift graphs in X direction if desired
    // *************************************************************************************************************        
    if(bWCorrection.Contains("X") ){
        TF1* fitShiftingEta                     = FitObject("tmpt","ShiftingEta","Eta");
        fitShiftingEta->SetParameters(paramGraphEta[0],paramGraphEta[1], paramGraphEta[2]) ; // standard parameter optimize if necessary
        
//         TF1* fitShiftingEta                     = FitObject("tcmpt","ShiftingEta","Eta");
//         fitShiftingEta->SetParameters(paramTCMEta[0],paramTCMEta[1], paramTCMEta[2], paramTCMEta[3], paramTCMEta[4]) ; // standard parameter optimize if necessary

        TGraphAsymmErrors* graphCombEtaInvYieldTotNoShift   = (TGraphAsymmErrors*) graphCombEtaInvYieldTot->Clone("Eta_NoShift");
        
        TGraphAsymmErrors* graphCombEtaInvYieldTotMultPt    = (TGraphAsymmErrors*) graphCombEtaInvYieldTot->Clone("Eta_MultPt");
        for (Int_t i = 0; i< graphCombEtaInvYieldTotMultPt->GetN(); i++){
            graphCombEtaInvYieldTotMultPt->SetPoint(i,graphCombEtaInvYieldTotMultPt->GetX()[i],graphCombEtaInvYieldTotMultPt->GetY()[i]*graphCombEtaInvYieldTotMultPt->GetX()[i]);
            graphCombEtaInvYieldTotMultPt->SetPointEYhigh(i,graphCombEtaInvYieldTotMultPt->GetErrorYhigh(i)*graphCombEtaInvYieldTotMultPt->GetX()[i]);
            graphCombEtaInvYieldTotMultPt->SetPointEYlow(i,graphCombEtaInvYieldTotMultPt->GetErrorYlow(i)*graphCombEtaInvYieldTotMultPt->GetX()[i]);
        }    
//         fitShiftingEta                          = FitObject("tmpt","ShiftingEta","Eta",graphCombEtaInvYieldTotMultPt,0.7,20.,paramGraphEta,"QNRMEX0+");
        fitShiftingEta                          = FitObject("tcmpt","ShiftingEta","Eta",graphCombEtaInvYieldTotMultPt,0.7,20.,paramTCMEta,"QNRMEX0+");
        
        graphCombEtaInvYieldTot->Print();
        
        
        graphCombEtaInvYieldTot                 = ApplyXshift(graphCombEtaInvYieldTot, fitShiftingEta);
        
        
        graphCombEtaInvYieldStat                = ApplyXshiftIndividualSpectra (   graphCombEtaInvYieldTot, 
                                                                                    graphCombEtaInvYieldStat, 
                                                                                    fitShiftingEta,
                                                                                    0, graphCombEtaInvYieldStat->GetN(),"Eta");
        graphCombEtaInvYieldSys                 = ApplyXshiftIndividualSpectra (   graphCombEtaInvYieldTot, 
                                                                                    graphCombEtaInvYieldSys, 
                                                                                    fitShiftingEta, 
                                                                                    0, graphCombEtaInvYieldSys->GetN(), "Eta");
        graphPCMEtaInvYieldStat                 = ApplyXshiftIndividualSpectra(    graphCombEtaInvYieldTot,
                                                                                    graphPCMEtaInvYieldStat,
                                                                                    fitShiftingEta, 
                                                                                    offSetEtaShifting[0], nComBinsEtaShifting[0], "Eta");
        graphPCMEtaInvYieldSys                  = ApplyXshiftIndividualSpectra(    graphCombEtaInvYieldTot, 
                                                                                    graphPCMEtaInvYieldSys, 
                                                                                    fitShiftingEta, 
                                                                                    offSetEtaShifting[0], nComBinsEtaShifting[0], "Eta");
        graphEMCALEtaInvYieldStat               = ApplyXshiftIndividualSpectra(    graphCombEtaInvYieldTot, 
                                                                                    graphEMCALEtaInvYieldStat, 
                                                                                    fitShiftingEta,
                                                                                    offSetEtaShifting[2], nComBinsEtaShifting[2], "Eta");
        graphEMCALEtaInvYieldSys                = ApplyXshiftIndividualSpectra(    graphCombEtaInvYieldTot, 
                                                                                    graphEMCALEtaInvYieldSys, 
                                                                                    fitShiftingEta,
                                                                                    offSetEtaShifting[2], nComBinsEtaShifting[2], "Eta");
        graphPCMEMCALEtaInvYieldStat            = ApplyXshiftIndividualSpectra(    graphCombEtaInvYieldTot,
                                                                                    graphPCMEMCALEtaInvYieldStat,
                                                                                    fitShiftingEta, 
                                                                                    offSetEtaShifting[4], nComBinsEtaShifting[4], "Eta");
        graphPCMEMCALEtaInvYieldSys             = ApplyXshiftIndividualSpectra(    graphCombEtaInvYieldTot, 
                                                                                    graphPCMEMCALEtaInvYieldSys, 
                                                                                    fitShiftingEta, 
                                                                                    offSetEtaShifting[4], nComBinsEtaShifting[4], "Eta");

        //***************************************************************************************************************
        //************************************Plotting binshift corrections *********************************************
        //***************************************************************************************************************

        TCanvas* canvasShift = new TCanvas("canvasShift","",0,0,1000,900);// gives the page size
        DrawGammaCanvasSettings( canvasShift, 0.12, 0.017, 0.015, 0.08);

        Size_t textSizeSpectra          = 0.04;
        TH1F * histoBinShift = new TH1F("histoBinShift","histoBinShift",1000,0., 20.);
        SetStyleHistoTH1ForGraphs(histoBinShift, "#it{p}_{T} (GeV/#it{c})","bin shifted (X) / no shift",
                                0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.2);
        histoBinShift->GetYaxis()->SetRangeUser(0.95,1.05);
        histoBinShift->DrawCopy();

        Int_t numberPoints   = graphCombEtaInvYieldTotNoShift->GetN();
        Double_t *xPoint     = graphCombEtaInvYieldTotNoShift->GetX();
        Double_t* xvalueErrUp  = graphCombEtaInvYieldTotNoShift->GetEXhigh();
        Double_t* xvalueErrLow = graphCombEtaInvYieldTotNoShift->GetEXlow();
        Double_t *xPointShift= graphCombEtaInvYieldTot->GetX();
        for (Int_t i=0; i<numberPoints; i++) {
          graphCombEtaInvYieldTotNoShift->SetPoint(i,xPoint[i],xPointShift[i]/xPoint[i]);
          graphCombEtaInvYieldTotNoShift->SetPointError(i,xvalueErrLow[i],xvalueErrUp[i],0,0);
        }
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldTotNoShift, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombEtaInvYieldTotNoShift->Draw("p same");

        TLatex *labelRatioToFitBinShift   = new TLatex(0.94, 0.91, collisionSystempPb.Data());
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
        DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.1);
        canvasDummy2->SetLogy();
        canvasDummy2->SetLogx();
        TH2F* histo2DDummy3;
        histo2DDummy3               = new TH2F("histo2DDummy3","histo2DDummy3",1000,0.23,25.,1000,1e-10,1e1);
        SetStyleHistoTH2ForGraphs(histo2DDummy3, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
        histo2DDummy3->DrawCopy(); 
    
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStatUnshi, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
        graphCombEtaInvYieldStatUnshi->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldTot, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombEtaInvYieldTot->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvYieldSys, markerStyleDet[0] ,markerSizeDet[0]/2, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        graphPCMEtaInvYieldSys->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvYieldStat, markerStyleDet[0] ,markerSizeDet[0]/2, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        graphPCMEtaInvYieldStat->Draw("psame");
        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaInvYieldSys, markerStyleDet[2] ,markerSizeDet[2]/2, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        graphEMCALEtaInvYieldSys->Draw("pEsame");


        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaInvYieldSys, markerStyleDet[4] ,markerSizeDet[4]/2, colorDet[4], colorDet[4], widthLinesBoxes, kTRUE);
        graphPCMEMCALEtaInvYieldSys->Draw("pEsame");

//         fitShiftingEta->SetLineColor(kRed+2);
//         fitShiftingEta->SetLineStyle(7);
        fitInvYieldEta->SetLineColor(kBlue+2);
        fitInvYieldEta->Draw("same");
//         fitShiftingEta->Draw("same");
        
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
    graphCombEtaInvYieldTot->Fit(fitInvYieldEta,"QNRMEX0+","",0.7, 20.);
    fitInvYieldEta        = FitObject("l","fitInvYieldEta2760GeV","Eta",graphCombEtaInvYieldTot,0.7, 20.,paramGraphEta,"QNRMEX0+");
    fitInvYieldEta        = FitObject("l","fitInvYieldEta2760GeV","Eta",graphCombEtaInvYieldTot,0.7, 20. ,paramGraphEta,"QNRMEX0+");
    cout << WriteParameterToFile(fitInvYieldEta)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitInvYieldEta)<< endl;    
    
    // Two component model by Bylinkin
    TF1* fitTCMDecomposedLEta   = FitObject("tcmlow","twoCompModelEta_DecL", "Eta", NULL, 0.7, 5.);
    TF1* fitTCMDecomposedHEta   = FitObject("tcmhigh","twoCompModelEta_DecH", "Eta", NULL, 0.7, 20.);
//     fitTCMDecomposedHEta->SetParameters(graphCombEtaInvYieldTot->GetY()[2],0.8, 2);    
//     Double_t paramTCMEtaNew[5]  = { fitTCMDecomposedLEta->GetParameter(0),fitTCMDecomposedLEta->GetParameter(1),
//                                     fitTCMDecomposedHEta->GetParameter(0),fitTCMDecomposedHEta->GetParameter(1),fitTCMDecomposedHEta->GetParameter(2)};
    
    TF1* fitTCMInvYieldEta= FitObject("tcm","fitTCMInvYieldEta2760GeV","Eta",graphCombEtaInvYieldTot,0.7, 20.,paramTCMEta,"QNRMEX0+","", kFALSE);
    fitTCMInvYieldEta     = FitObject("tcm","fitTCMInvYieldEta2760GeV","Eta",graphCombEtaInvYieldTot,0.7, 20. ,paramTCMEta,"QNRMEX0+","", kFALSE);
    fitTCMDecomposedLEta->SetParameter(0, fitTCMInvYieldEta->GetParameter(0));
    fitTCMDecomposedLEta->SetParameter(1, fitTCMInvYieldEta->GetParameter(1));
    fitTCMDecomposedHEta->SetParameter(0, fitTCMInvYieldEta->GetParameter(2));
    fitTCMDecomposedHEta->FixParameter(1, fitTCMInvYieldEta->GetParameter(3));
    fitTCMDecomposedHEta->FixParameter(2, fitTCMInvYieldEta->GetParameter(4));
    graphCombEtaInvYieldTot->Fit(fitTCMDecomposedHEta,"QNRMEX0+","",3,20);
    
    cout << WriteParameterToFile(fitTCMInvYieldEta)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitTCMInvYieldEta)<< endl;    
    
    // *************************************************************************************************************
    // Shift spectra in Y  direction as well if desired
    // *************************************************************************************************************

    TGraphAsymmErrors* graphCombEtaInvYieldTot_yShifted         = NULL;
    TGraphAsymmErrors* graphCombEtaInvYieldStat_yShifted        = NULL;
    TGraphAsymmErrors* graphCombEtaInvYieldSys_yShifted         = NULL;
    TGraphAsymmErrors* graphPCMEtaInvYieldStat_yShifted          = NULL;
    TGraphAsymmErrors* graphPCMEtaInvYieldSys_yShifted           = NULL;
    TGraphAsymmErrors* graphEMCALEtaInvYieldStat_yShifted        = NULL;
    TGraphAsymmErrors* graphEMCALEtaInvYieldSys_yShifted         = NULL;
    TGraphAsymmErrors* graphPCMEMCALEtaInvYieldStat_yShifted     = NULL;
    TGraphAsymmErrors* graphPCMEMCALEtaInvYieldSys_yShifted      = NULL;
        
    if(bWCorrection.Contains("Y") ){
        graphCombEtaInvYieldTot_yShifted        = (TGraphAsymmErrors*)graphCombEtaInvYieldTotUnshi->Clone("EtaYShiftedCombTot");
        graphCombEtaInvYieldTot_yShifted        =  ApplyYshiftIndividualSpectra( graphCombEtaInvYieldTot_yShifted, fitInvYieldEta);
        graphCombEtaInvYieldStat_yShifted       = (TGraphAsymmErrors*)graphCombEtaInvYieldStatUnshi->Clone("EtaYShiftedCombStat");
        graphCombEtaInvYieldStat_yShifted       =  ApplyYshiftIndividualSpectra( graphCombEtaInvYieldStat_yShifted, fitInvYieldEta);
        graphCombEtaInvYieldSys_yShifted        = (TGraphAsymmErrors*)graphCombEtaInvYieldSysUnshi->Clone("EtaYShiftedCombSys");
        graphCombEtaInvYieldSys_yShifted        =  ApplyYshiftIndividualSpectra( graphCombEtaInvYieldSys_yShifted, fitInvYieldEta);
        
        graphPCMEtaInvYieldStat_yShifted         = (TGraphAsymmErrors*)graphPCMEtaInvYieldStatUnshi->Clone("EtaYShiftedPCMStat");
        graphPCMEtaInvYieldStat_yShifted         =  ApplyYshiftIndividualSpectra( graphPCMEtaInvYieldStat_yShifted, fitInvYieldEta);
        graphPCMEtaInvYieldSys_yShifted          = (TGraphAsymmErrors*)graphPCMEtaInvYieldSysUnshi->Clone("EtaYShiftedPCMSys");
        graphPCMEtaInvYieldSys_yShifted          =  ApplyYshiftIndividualSpectra( graphPCMEtaInvYieldSys_yShifted, fitInvYieldEta);
        
        graphEMCALEtaInvYieldStat_yShifted       = (TGraphAsymmErrors*)graphEMCALEtaInvYieldStatUnshi->Clone("EtaYShiftedEMCStat");
        graphEMCALEtaInvYieldStat_yShifted       =  ApplyYshiftIndividualSpectra( graphEMCALEtaInvYieldStat_yShifted, fitInvYieldEta);
        graphEMCALEtaInvYieldSys_yShifted        = (TGraphAsymmErrors*)graphEMCALEtaInvYieldSysUnshi->Clone("EtaYShiftedEMCSys");
        graphEMCALEtaInvYieldSys_yShifted        =  ApplyYshiftIndividualSpectra( graphEMCALEtaInvYieldSys_yShifted, fitInvYieldEta);
        
        graphPCMEMCALEtaInvYieldStat_yShifted    = (TGraphAsymmErrors*)graphPCMEMCALEtaInvYieldStatUnshi->Clone("EtaYShiftedPCMEMCStat");
        graphPCMEMCALEtaInvYieldStat_yShifted    =  ApplyYshiftIndividualSpectra( graphPCMEMCALEtaInvYieldStat_yShifted, fitInvYieldEta);
        graphPCMEMCALEtaInvYieldSys_yShifted     = (TGraphAsymmErrors*)graphPCMEMCALEtaInvYieldSysUnshi->Clone("EtaYShiftedPCMEMCStat");
        graphPCMEMCALEtaInvYieldSys_yShifted     =  ApplyYshiftIndividualSpectra( graphPCMEMCALEtaInvYieldSys_yShifted, fitInvYieldEta);
    }

    // *************************************************************************************************************
    // Calculate ratios to combined fit
    // *************************************************************************************************************    
//     TH1D* histoRatioEtaPythia8ToFit                     = (TH1D*) histoPythia8InvYieldEta->Clone(); 
//     histoRatioEtaPythia8ToFit                           = CalculateHistoRatioToFit (histoRatioEtaPythia8ToFit, fitTCMInvYieldEta); 
//     histoRatioEtaPythia8ToFit->GetXaxis()->SetRangeUser(0.6,25);
// 
//     TGraph* graphRatioEtaCombNLOMuHalf                  = (TGraph*)graphNLOCalcEtaMuHalf->Clone();
//     TGraph* graphRatioEtaCombNLOMuOne                   = (TGraph*)graphNLOCalcEtaMuOne->Clone();
//     TGraph* graphRatioEtaCombNLOMuTwo                   = (TGraph*)graphNLOCalcEtaMuTwo->Clone();
//     graphRatioEtaCombNLOMuHalf                          = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuHalf, fitTCMInvYieldEta); 
//     graphRatioEtaCombNLOMuOne                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuOne, fitTCMInvYieldEta); 
//     graphRatioEtaCombNLOMuTwo                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuTwo, fitTCMInvYieldEta); 
    
    TGraphAsymmErrors* graphRatioEtaCombCombFitTot     = (TGraphAsymmErrors*)graphCombEtaInvYieldTot->Clone();
    graphRatioEtaCombCombFitTot                        = CalculateGraphErrRatioToFit(graphRatioEtaCombCombFitTot, fitTCMInvYieldEta); 
    TGraphAsymmErrors* graphRatioEtaCombCombFitStat    = (TGraphAsymmErrors*)graphCombEtaInvYieldStat->Clone();
    graphRatioEtaCombCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioEtaCombCombFitStat, fitTCMInvYieldEta); 
    TGraphAsymmErrors* graphRatioEtaCombCombFitSys     = (TGraphAsymmErrors*)graphCombEtaInvYieldSys->Clone();

    graphRatioEtaCombCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioEtaCombCombFitSys, fitTCMInvYieldEta); 
    TGraphAsymmErrors* graphRatioEtaPCMCombFitStat      = (TGraphAsymmErrors*)graphPCMEtaInvYieldStat->Clone();
    graphRatioEtaPCMCombFitStat                         = CalculateGraphErrRatioToFit(graphRatioEtaPCMCombFitStat, fitTCMInvYieldEta); 
    TGraphAsymmErrors* graphRatioEtaPCMCombFitSys       = (TGraphAsymmErrors*)graphPCMEtaInvYieldSys->Clone();
    graphRatioEtaPCMCombFitSys                          = CalculateGraphErrRatioToFit(graphRatioEtaPCMCombFitSys, fitTCMInvYieldEta); 
    TGraphAsymmErrors* graphRatioEtaEMCALCombFitStat    = (TGraphAsymmErrors*)graphEMCALEtaInvYieldStat->Clone();
    graphRatioEtaEMCALCombFitStat                       = CalculateGraphErrRatioToFit(graphRatioEtaEMCALCombFitStat, fitTCMInvYieldEta); 
    TGraphAsymmErrors* graphRatioEtaEMCALCombFitSys     = (TGraphAsymmErrors*)graphEMCALEtaInvYieldSys->Clone();
    graphRatioEtaEMCALCombFitSys                        = CalculateGraphErrRatioToFit(graphRatioEtaEMCALCombFitSys, fitTCMInvYieldEta); 
    TGraphAsymmErrors* graphRatioEtaPCMEMCALCombFitStat = (TGraphAsymmErrors*)graphPCMEMCALEtaInvYieldStat->Clone();
    graphRatioEtaPCMEMCALCombFitStat                    = CalculateGraphErrRatioToFit(graphRatioEtaPCMEMCALCombFitStat, fitTCMInvYieldEta); 
    TGraphAsymmErrors* graphRatioEtaPCMEMCALCombFitSys  = (TGraphAsymmErrors*)graphPCMEMCALEtaInvYieldSys->Clone();
    graphRatioEtaPCMEMCALCombFitSys                     = CalculateGraphErrRatioToFit(graphRatioEtaPCMEMCALCombFitSys, fitTCMInvYieldEta); 
    

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

        ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombCombFitStat);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioEtaCombCombFitSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStat, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioEtaCombCombFitStat->Draw("p,same,z");

        DrawGammaLines(0.33, 25. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.33, 25. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.33, 25. , 0.9, 0.9,0.1, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        TLatex *labelRatioToFitEta      = new TLatex(0.852,0.82,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitEta, textSizeLabelsPixel,4, 1, 43);
        labelRatioToFitEta->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfCombToCombFit_PPb5023GeV.%s",outputDir.Data(),suffix.Data()));
    
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
        TLatex *textPCMOnlyRatioEta                 = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[1],nameMeasGlobalLabel[0]);
        SetStyleTLatex( textPCMOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
        textPCMOnlyRatioEta->Draw();
        TLatex *textEMCALOnlyRatioEta               = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[2],nameMeasGlobalLabel[2]);
        SetStyleTLatex( textEMCALOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
        textEMCALOnlyRatioEta->Draw();
        TLatex *textPCMEMCALOnlyRatioEta            = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[3],nameMeasGlobalLabel[4]);
        SetStyleTLatex( textPCMEMCALOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
        textPCMEMCALOnlyRatioEta->Draw();
        
        //****************** second Column *************************************************
        TLatex *textStatOnlyRatioEta                = new TLatex(columnsLegendOnlyEtaRatio[1],rowsLegendOnlyEtaRatio[0] ,"stat");
        SetStyleTLatex( textStatOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
        textStatOnlyRatioEta->Draw();
        TLatex *textSysOnlyRatioEta                 = new TLatex(columnsLegendOnlyEtaRatio[2] ,rowsLegendOnlyEtaRatio[0],"syst");
        SetStyleTLatex( textSysOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
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
   
    // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)    
    // for systematic error from respective method
    TGraphAsymmErrors* sysErrorCollectionEtaToPi0[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorCollectionEtaToPi0[i]    = NULL;
    }    
    sysErrorCollectionEtaToPi0[0]        = (TGraphAsymmErrors*)graphPCMEtaToPi0Sys->Clone("sysErrPCMEtaToPi0");
    sysErrorCollectionEtaToPi0[2]        = (TGraphAsymmErrors*)graphEMCALEtaToPi0Sys->Clone("sysErrEMCALEtaToPi0");

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
    Int_t offSetsEtaToPi0Sys[11]        = { 3,  0,  8,  0,  5,
                                            0,  0,  0,  0,  0,
                                            0};

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************
                                
    TGraph* graphWeightsEtaToPi0WOPE[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsEtaToPi0WOPE[i]                    = NULL;
    }    
                                
    // Declaration & calculation of combined spectrum                            
    TString fileNameEtaToPi0OutputWeightingWOPE    = Form("%s/EtaToPi0_WeightingWOPE.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombEtaToPi0StatWOPE   = NULL;
    TGraphAsymmErrors* graphCombEtaToPi0SysWOPE    = NULL;
    TGraphAsymmErrors* graphCombEtaToPi0TotWOPE    = CombinePtPointsSpectraFullCorrMat(   statErrorCollectionEtaToPi0,    sysErrorCollectionEtaToPi0,     
                                                                                           xPtLimitsEta, maxNBinsEta-1,
                                                                                           offSetsEtaToPi0, offSetsEtaToPi0Sys,
                                                                                           graphCombEtaToPi0StatWOPE, graphCombEtaToPi0SysWOPE,
                                                                                           fileNameEtaToPi0OutputWeightingWOPE,"2.76TeV", "EtaToPi0", kFALSE,
                                                                                           NULL, fileNameCorrFactors
                                                                                        );

    
    if (graphCombEtaToPi0StatWOPE == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    } else {    
        while(graphCombEtaToPi0StatWOPE->GetX()[0] < 0.8){
            graphCombEtaToPi0StatWOPE->RemovePoint(0);
            graphCombEtaToPi0SysWOPE->RemovePoint(0);
            graphCombEtaToPi0TotWOPE->RemovePoint(0);
        }
        graphCombEtaToPi0StatWOPE->Print();
    }    
    
    // Reading weights from output file for plotting
    ifstream fileWeightsEtaToPi0ReadWOPE;
    fileWeightsEtaToPi0ReadWOPE.open(fileNameEtaToPi0OutputWeightingWOPE,ios_base::in);
    cout << "reading" << fileNameEtaToPi0OutputWeightingWOPE << endl;
    Double_t xValuesEtaToPi0ReadWOPE[50];
    Double_t weightsEtaToPi0ReadWOPE[11][50];
    Int_t availableEtaToPi0MeasWOPE[11]         = { -1, -1, -1, -1, -1,
                                                    -1, -1, -1, -1, -1,
                                                    -1};
    Int_t nMeasSetEtaToPi0WOPE                  = 2;
    Int_t nPtBinsEtaToPi0ReadWOPE               = 0;
    while(!fileWeightsEtaToPi0ReadWOPE.eof() && nPtBinsEtaToPi0ReadWOPE < 50){
        TString garbage             = "";
        if (nPtBinsEtaToPi0ReadWOPE == 0){
            fileWeightsEtaToPi0ReadWOPE >> garbage ;//>> availableEtaToPi0MeasWOPE[0] >> availableEtaToPi0MeasWOPE[1] >> availableEtaToPi0MeasWOPE[2] >> availableEtaToPi0MeasWOPE[3];
            for (Int_t i = 0; i < nMeasSetEtaToPi0WOPE; i++){
                fileWeightsEtaToPi0ReadWOPE >> availableEtaToPi0MeasWOPE[i] ;
            }    
            cout << "read following measurements: "; 
            for (Int_t i = 0; i < 11; i++){
                cout << availableEtaToPi0MeasWOPE[i] << "\t" ;
            }    
            cout << endl;
        } else {
            fileWeightsEtaToPi0ReadWOPE >> xValuesEtaToPi0ReadWOPE[nPtBinsEtaToPi0ReadWOPE-1];
            for (Int_t i = 0; i < nMeasSetEtaToPi0WOPE; i++){
                fileWeightsEtaToPi0ReadWOPE >> weightsEtaToPi0ReadWOPE[availableEtaToPi0MeasWOPE[i]][nPtBinsEtaToPi0ReadWOPE-1] ;
            }    
            cout << "read: "<<  nPtBinsEtaToPi0ReadWOPE << "\t"<< xValuesEtaToPi0ReadWOPE[nPtBinsEtaToPi0ReadWOPE-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetEtaToPi0WOPE; i++){
                cout << weightsEtaToPi0ReadWOPE[availableEtaToPi0MeasWOPE[i]][nPtBinsEtaToPi0ReadWOPE-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsEtaToPi0ReadWOPE++;
    }
    nPtBinsEtaToPi0ReadWOPE                    = nPtBinsEtaToPi0ReadWOPE-2 ;
    fileWeightsEtaToPi0ReadWOPE.close();

    for (Int_t i = 0; i < nMeasSetEtaToPi0WOPE; i++){
        graphWeightsEtaToPi0WOPE[availableEtaToPi0MeasWOPE[i]]          = new TGraph(nPtBinsEtaToPi0ReadWOPE,xValuesEtaToPi0ReadWOPE,weightsEtaToPi0ReadWOPE[availableEtaToPi0MeasWOPE[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsEtaToPi0ReadWOPE; n++){
            if (graphWeightsEtaToPi0WOPE[availableEtaToPi0MeasWOPE[i]]->GetY()[bin] == 0) graphWeightsEtaToPi0WOPE[availableEtaToPi0MeasWOPE[i]]->RemovePoint(bin);
            else bin++;
        }    
    }    
    // **********************************************************************************************************************
    // ************************ Adding PCM-EMC measurements to matrix & calculating relativ errors ************************
    // **********************************************************************************************************************
    
    statErrorCollectionEtaToPi0[4]          = (TH1D*)histoPCMEMCALEtaToPi0Stat->Clone("statErrPCMEMCALEtaToPi0");
    sysErrorCollectionEtaToPi0[4]           = (TGraphAsymmErrors*)graphPCMEMCALEtaToPi0Sys->Clone("sysErrPCMEMCALEtaToPi0");

    if (statErrorCollectionEtaToPi0[4]){
        statErrorRelCollectionEtaToPi0[4]    = CalculateRelErrUpTH1D( statErrorCollectionEtaToPi0[4], Form("relativeStatErrorEtaToPi0_%s", nameMeasGlobal[4].Data()));
//         statErrorRelCollectionEta[4]        = new TGraphAsymmErrors(statErrorCollectionEta[4]);
//         while (statErrorRelCollectionEta[4]->GetY()[0] == 0) 
//             statErrorRelCollectionEta[4]->RemovePoint(0);
//         statErrorRelCollectionEta[4]        = CalculateRelErrUpAsymmGraph( statErrorRelCollectionEta[4], Form("relativeStatErrorEta_%s", nameMeasGlobal[4].Data()));
    }        
    if (sysErrorCollectionEta[4]) 
        sysErrorRelCollectionEtaToPi0[4]     = CalculateRelErrUpAsymmGraph( sysErrorCollectionEtaToPi0[4], Form("relativeSysErrorEtaToPi0_%s", nameMeasGlobal[4].Data()));
                                            
    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************
                                
    TGraph* graphWeightsEtaToPi0[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeightsEtaToPi0[i]                    = NULL;
    }    
                                
    // Declaration & calculation of combined spectrum                            
    TString fileNameEtaToPi0OutputWeighting        = Form("%s/EtaToPi0_WeightingMethod.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombEtaToPi0Stat       = NULL;
    TGraphAsymmErrors* graphCombEtaToPi0Sys        = NULL;
    TGraphAsymmErrors* graphCombEtaToPi0Tot        = CombinePtPointsSpectraFullCorrMat(   statErrorCollectionEtaToPi0,    sysErrorCollectionEtaToPi0,     
                                                                                           xPtLimitsEta, maxNBinsEta-1,
                                                                                           offSetsEtaToPi0, offSetsEtaToPi0Sys,
                                                                                           graphCombEtaToPi0Stat, graphCombEtaToPi0Sys,
                                                                                           fileNameEtaToPi0OutputWeighting,"2.76TeV", "EtaToPi0", kFALSE,
                                                                                           NULL, fileNameCorrFactors
                                                                                        );

    
    if (graphCombEtaToPi0Stat == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    } else {    
        while(graphCombEtaToPi0Stat->GetX()[0] < 0.8){
            graphCombEtaToPi0Stat->RemovePoint(0);
            graphCombEtaToPi0Sys->RemovePoint(0);
            graphCombEtaToPi0Tot->RemovePoint(0);
        }
        graphCombEtaToPi0Stat->Print();
    }    
    
    // Reading weights from output file for plotting
    ifstream fileWeightsEtaToPi0Read;
    fileWeightsEtaToPi0Read.open(fileNameEtaToPi0OutputWeighting,ios_base::in);
    cout << "reading" << fileNameEtaToPi0OutputWeighting << endl;
    Double_t xValuesEtaToPi0Read[50];
    Double_t weightsEtaToPi0Read[11][50];
    Int_t availableEtaToPi0Meas[11]         = { -1, -1, -1, -1, -1,
                                                -1, -1, -1, -1, -1,
                                                -1};
    Int_t nMeasSetEtaToPi0                  = 3;
    Int_t nPtBinsEtaToPi0Read               = 0;
    while(!fileWeightsEtaToPi0Read.eof() && nPtBinsEtaToPi0Read < 50){
        TString garbage             = "";
        if (nPtBinsEtaToPi0Read == 0){
            fileWeightsEtaToPi0Read >> garbage ;//>> availableEtaToPi0Meas[0] >> availableEtaToPi0Meas[1] >> availableEtaToPi0Meas[2] >> availableEtaToPi0Meas[3];
            for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                fileWeightsEtaToPi0Read >> availableEtaToPi0Meas[i] ;
            }    
            cout << "read following measurements: "; 
            for (Int_t i = 0; i < 11; i++){
                cout << availableEtaToPi0Meas[i] << "\t" ;
            }    
            cout << endl;
        } else {
            fileWeightsEtaToPi0Read >> xValuesEtaToPi0Read[nPtBinsEtaToPi0Read-1];
            for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                fileWeightsEtaToPi0Read >> weightsEtaToPi0Read[availableEtaToPi0Meas[i]][nPtBinsEtaToPi0Read-1] ;
            }    
            cout << "read: "<<  nPtBinsEtaToPi0Read << "\t"<< xValuesEtaToPi0Read[nPtBinsEtaToPi0Read-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                cout << weightsEtaToPi0Read[availableEtaToPi0Meas[i]][nPtBinsEtaToPi0Read-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsEtaToPi0Read++;
    }
    nPtBinsEtaToPi0Read                    = nPtBinsEtaToPi0Read-2 ;
    fileWeightsEtaToPi0Read.close();

    for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
        graphWeightsEtaToPi0[availableEtaToPi0Meas[i]]                        = new TGraph(nPtBinsEtaToPi0Read,xValuesEtaToPi0Read,weightsEtaToPi0Read[availableEtaToPi0Meas[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsEtaToPi0Read; n++){
            if (graphWeightsEtaToPi0[availableEtaToPi0Meas[i]]->GetY()[bin] == 0) graphWeightsEtaToPi0[availableEtaToPi0Meas[i]]->RemovePoint(bin);
            else bin++;
        }    
    }    

    // **********************************************************************************************************************
    // ******************************************* Plotting weights Method **************************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel           = 900*0.04;

    canvasWeights->cd();
   
    TH2F * histo2DEtaToPi0Weights;
    histo2DEtaToPi0Weights = new TH2F("histo2DEtaToPi0Weights","histo2DEtaToPi0Weights",11000,0.33,25.,1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DEtaToPi0Weights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DEtaToPi0Weights->GetXaxis()->SetMoreLogLabels();
    histo2DEtaToPi0Weights->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DEtaToPi0Weights->GetYaxis()->SetRangeUser(-10,10);    
    histo2DEtaToPi0Weights->Draw("copy");
        TLegend* legendWeightsEtaToPi0WOPE= GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEtaToPi0WOPE), 32);
        for (Int_t i = 0; i < nMeasSetEtaToPi0WOPE; i++){
            DrawGammaSetMarkerTGraph(graphWeightsEtaToPi0WOPE[availableEtaToPi0MeasWOPE[i]], markerStyleDet[availableEtaToPi0MeasWOPE[i]], markerSizeDet[availableEtaToPi0MeasWOPE[i]]*0.5, colorDet[availableEtaToPi0MeasWOPE[i]] , colorDet[availableEtaToPi0Meas[i]]);
            graphWeightsEtaToPi0WOPE[availableEtaToPi0MeasWOPE[i]]->Draw("p,same,z");
            legendWeightsEtaToPi0WOPE->AddEntry(graphWeightsEtaToPi0[availableEtaToPi0MeasWOPE[i]],nameMeasGlobalLabel[availableEtaToPi0MeasWOPE[i]],"p");
        }    
        legendWeightsEtaToPi0WOPE->Draw();

        labelWeightsEnergy->Draw();
        TLatex *labelWeightsEtaToPi0         = new TLatex(0.7,0.16,"#eta/#pi^{0}");
        SetStyleTLatex( labelWeightsEtaToPi0, 0.85*textSizeLabelsPixel,4, 1, 43);
        labelWeightsEtaToPi0->Draw();

//      DrawGammaLines(0.33, 25. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.33, 25. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.33, 25. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.33, 25. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.33, 25. , 0.2, 0.2,0.1, kGray, 3);
        
    canvasWeights->SaveAs(Form("%s/EtaToPi0_Weights_WOPCMEMC.%s",outputDir.Data(),suffix.Data()));

        TLegend* legendWeightsEtaToPi0   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEtaToPi0), 32);
        for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
            DrawGammaSetMarkerTGraph(graphWeightsEtaToPi0[availableEtaToPi0Meas[i]], markerStyleDet[availableEtaToPi0Meas[i]], markerSizeDet[availableEtaToPi0Meas[i]]*0.5, colorDet[availableEtaToPi0Meas[i]] , colorDet[availableEtaToPi0Meas[i]]);
            graphWeightsEtaToPi0[availableEtaToPi0Meas[i]]->Draw("p,same,z");
            legendWeightsEtaToPi0->AddEntry(graphWeightsEtaToPi0[availableEtaToPi0Meas[i]],nameMeasGlobalLabel[availableEtaToPi0Meas[i]],"p");
        }    
        legendWeightsEtaToPi0->Draw();

        labelWeightsEnergy->Draw();
        labelWeightsEtaToPi0->Draw();

        DrawGammaLines(0.33, 25. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.33, 25. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.33, 25. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.33, 25. , 0.2, 0.2,0.1, kGray, 3);
        
    canvasWeights->SaveAs(Form("%s/EtaToPi0_Weights.%s",outputDir.Data(),suffix.Data()));


    // *********************************************************************************************************************
    // ************************************ Visualize relative errors EtaToPi0 ******************************************************
    // *********************************************************************************************************************
    
    canvasRelSysErr->cd();
   
    TH2F * histo2DRelSysErrEtaToPi0;
    histo2DRelSysErrEtaToPi0                 = new TH2F("histo2DRelSysErrEtaToPi0","histo2DRelSysErrEtaToPi0",11000,0.33,25.,1000,0,60.5);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErrEtaToPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErrEtaToPi0->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelSysErrEtaToPi0->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelSysErrEtaToPi0->Draw("copy");

        TLegend* legendRelSysErrEtaToPi0        = GetAndSetLegend2(0.62, 0.92-(0.035*nMeasSetEtaToPi0), 0.95, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionEtaToPi0[availableEtaToPi0Meas[i]], markerStyleDet[availableEtaToPi0Meas[i]], markerSizeDet[availableEtaToPi0Meas[i]]*0.5, 
                                     colorDet[availableEtaToPi0Meas[i]], colorDet[availableEtaToPi0Meas[i]]);
            sysErrorRelCollectionEtaToPi0[availableEtaToPi0Meas[i]]->Draw("p,same,z");
            legendRelSysErrEtaToPi0->AddEntry(sysErrorRelCollectionEtaToPi0[availableEtaToPi0Meas[i]],nameMeasGlobalLabel[availableEtaToPi0Meas[i]],"p");
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
    histo2DRelStatErrEtaToPi0                = new TH2F("histo2DRelStatErrEtaToPi0","histo2DRelStatErrEtaToPi0",11000,0.33,25.,1000,0,60.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErrEtaToPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErrEtaToPi0->GetXaxis()->SetLabelOffset(-0.01);
//  histo2DRelStatErrEtaToPi0->GetYaxis()->SetRangeUser(-10,10);
    histo2DRelStatErrEtaToPi0->Draw("copy");
        
        TLegend* legendRelStatErrEtaToPi0       = GetAndSetLegend2(0.14, 0.92-(0.035*nMeasSetEtaToPi0), 0.45, 0.92, 32);
        for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
            DrawGammaSetMarker(statErrorRelCollectionEtaToPi0[availableEtaToPi0Meas[i]], markerStyleDet[availableEtaToPi0Meas[i]], markerSizeDet[availableEtaToPi0Meas[i]]*0.5,
                               colorDet[availableEtaToPi0Meas[i]] , colorDet[availableEtaToPi0Meas[i]]);
            statErrorRelCollectionEtaToPi0[availableEtaToPi0Meas[i]]->Draw("p,same,z");
            legendRelStatErrEtaToPi0->AddEntry(statErrorRelCollectionEtaToPi0[availableEtaToPi0Meas[i]],nameMeasGlobalLabel[availableEtaToPi0Meas[i]],"p");
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
    TGraphAsymmErrors* graphCombEtaToPi0RelStat       = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0Stat, "relativeStatErrorEtaToPi0_Method");
    TGraphAsymmErrors* graphCombEtaToPi0RelSys        = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0Sys, "relativeSysErrorEtaToPi0_Method");
    TGraphAsymmErrors* graphCombEtaToPi0RelTot        = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0Tot, "relativeTotalErrorEtaToPi0_Method");
    TGraphAsymmErrors* graphCombEtaToPi0RelStatWOPE   = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0StatWOPE, "relativeStatErrorEtaToPi0_MethodWOPE");
    TGraphAsymmErrors* graphCombEtaToPi0RelSysWOPE    = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0SysWOPE, "relativeSysErrorEtaToPi0_MethodWOPE");
    TGraphAsymmErrors* graphCombEtaToPi0RelTotWOPE    = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0TotWOPE, "relativeTotalErrorEtaToPi0_MethodWOPE");
    
    while(graphCombEtaToPi0RelStat->GetX()[0] < 0.8) graphCombEtaToPi0RelStat->RemovePoint(0);
    while(graphCombEtaToPi0RelSys->GetX()[0] < 0.8) graphCombEtaToPi0RelSys->RemovePoint(0);
    while(graphCombEtaToPi0RelTot->GetX()[0] < 0.8) graphCombEtaToPi0RelTot->RemovePoint(0);
    while(graphCombEtaToPi0RelStatWOPE->GetX()[0] < 0.8) graphCombEtaToPi0RelStatWOPE->RemovePoint(0);
    while(graphCombEtaToPi0RelSysWOPE->GetX()[0] < 0.8) graphCombEtaToPi0RelSysWOPE->RemovePoint(0);
    while(graphCombEtaToPi0RelTotWOPE->GetX()[0] < 0.8) graphCombEtaToPi0RelTotWOPE->RemovePoint(0);

    canvasRelTotErr->cd();
    TH2F * histo2DRelTotErrEtaToPi0;
    histo2DRelTotErrEtaToPi0                 = new TH2F("histo2DRelTotErrEtaToPi0","histo2DRelTotErrEtaToPi0",11000,0.33,25.,1000,0,60.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrEtaToPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrEtaToPi0->GetXaxis()->SetLabelOffset(-0.01);    
    histo2DRelTotErrEtaToPi0->Draw("copy");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelTotWOPE, markerStyleComb, markerSizeComb, kBlue+2 , kBlue+2);
        graphCombEtaToPi0RelTotWOPE->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelTot, markerStyleComb+4, markerSizeComb, kRed+2 , kRed+2);
        graphCombEtaToPi0RelTot->Draw("p,same,z");

        TLegend* legendRelTotErrEtaToPi0     = GetAndSetLegend2(0.14, 0.92-(0.035*2), 0.45, 0.92, 32);
        legendRelTotErrEtaToPi0->AddEntry(graphCombEtaToPi0RelTotWOPE,"PCM, EMC","p");
        legendRelTotErrEtaToPi0->AddEntry(graphCombEtaToPi0RelTot,"All","p");
        legendRelTotErrEtaToPi0->Draw();

        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrEtaToPi0       = new TLatex(0.75,0.85,"#eta/#pi^{0}");
        SetStyleTLatex( labelRelTotErrEtaToPi0, 0.85*textSizeLabelsPixel,4, 1, 43);
        labelRelTotErrEtaToPi0->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/EtaToPi0_RelTotErr.%s",outputDir.Data(),suffix.Data()));
            
    histo2DRelTotErrEtaToPi0->GetYaxis()->SetRangeUser(0,50.5);
    histo2DRelTotErrEtaToPi0->GetYaxis()->SetTitle("Err (%)");
    histo2DRelTotErrEtaToPi0->Draw("copy");
        
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombEtaToPi0RelTot->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombEtaToPi0RelStat->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombEtaToPi0RelSys->SetLineStyle(7);
        graphCombEtaToPi0RelSys->Draw("l,x0,same,e1");

        legendRelTotErr3->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrEtaToPi0->Draw();
        
    canvasRelTotErr->SaveAs(Form("%s/EtaToPi0_RelMethoddecomp.%s",outputDir.Data(),suffix.Data()));

    
    // **********************************************************************************************************************
    // ******************************************* Mass and width for pi0 ***************************************************
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
    
    TPad* padMassLegend1            = new TPad("padMassLegend1", "", 0.13, 0.36, 0.82, 0.52,-1, -1, -2);
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
        histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,20.5);
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
        TLatex *labelMassEnergy     = new TLatex(0.13,0.78,collisionSystempPb.Data());
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
        histo2DAllPi0Mass->GetYaxis()->SetRangeUser(126.5,150.5);
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
        Double_t columnsLegendMass2[6]      = {0.,0.25,0.38, 0.5, 0.75, 0.9};
        Double_t rowsLegendMass2[4]         = {0.75,0.5,0.25,0.01};
        //******************* Offsets ***********************
        Double_t offsetMarkerXMass2         = 0.05;
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
        markerPCMPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.01 ,rowsLegendMass2[1]+ offsetMarkerYMass2);
        TMarker* markerPCMEMCALPi0MassMC = CreateMarkerFromGraph(graphPCMEMCALPi0MassMC,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[2]+ offsetMarkerYMass2 ,scaleMarkerMass2);
        markerPCMEMCALPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.01 ,rowsLegendMass2[2]+ offsetMarkerYMass2);
        TMarker* markerEMCALPi0MassMC    = CreateMarkerFromGraph(graphEMCALPi0MassMC,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[3]+ offsetMarkerYMass2 ,scaleMarkerMass2);
        markerEMCALPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.01 ,rowsLegendMass2[3]+ offsetMarkerYMass2);
        
    canvasMassWidthPi0->Update();
    canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // *********************************** Mass and width for pi0 at 2.76TeV including PHOS *********************************
    // **********************************************************************************************************************
    
    canvasMassWidthPi0->cd();
    padWidthPi0->Draw();
    padMassPi0->Draw();
    
    TPad* padMassLegend2                = new TPad("padMassLegend2", "", 0.13, 0.36, 0.82, 0.52,-1, -1, -2);
    DrawGammaPadSettings( padMassLegend2, 0., 0., 0., 0.);
    padMassLegend2->SetFillStyle(0);
    padMassLegend2->Draw();
    
    padWidthPi0->cd();
    padWidthPi0->SetLogx(); 

        histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,20.5);
        histo2DAllPi0FWHM->DrawCopy(); 

        DrawGammaSetMarker(histoPHOSPi0FWHMMeV, markerStyleDet[1], markerSizeDet[1]*0.55, colorDet[1] , colorDet[1]);
        histoPHOSPi0FWHMMeV->Draw("p,same,e");
        DrawGammaSetMarker(histoPHOSPi0TrueFWHMMeV, markerStyleDetMC[1], markerSizeDetMC[1]*0.55, colorDetMC[1] , colorDetMC[1]);
        histoPHOSPi0TrueFWHMMeV->Draw("p,same,e");

        DrawGammaSetMarkerTGraphAsym(graphDalitzPi0FWHM, markerStyleDet[5], markerSizeDet[5]*0.55, colorDet[5] , colorDet[5]);
        graphDalitzPi0FWHM->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphDalitzPi0FWHMMC, markerStyleDetMC[5], markerSizeDetMC[5]*0.55, colorDetMC[5] , colorDetMC[5]);
        graphDalitzPi0FWHMMC->Draw("p,same,z");
        
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

        DrawGammaSetMarkerTGraphAsym(graphDalitzPi0Mass, markerStyleDet[5], markerSizeDet[5]*0.55, colorDet[5] , colorDet[5]);
        graphDalitzPi0Mass->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphDalitzPi0MassMC, markerStyleDetMC[5], markerSizeDetMC[5]*0.55, colorDetMC[5] , colorDetMC[5]);
        graphDalitzPi0MassMC->Draw("p,same,z");
        
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
        TLatex *textMassPHOS2               = new TLatex(columnsLegendMass2[3],rowsLegendMass3[1],nameMeasGlobalLabel[1]);
        SetStyleTLatex( textMassPHOS2, textSizeLabelsPixel,4, 1, 43);
        textMassPHOS2->Draw();
        TLatex *textMassDalitz2               = new TLatex(columnsLegendMass2[3],rowsLegendMass3[2],nameMeasGlobalLabel[5]);
        SetStyleTLatex( textMassDalitz2, textSizeLabelsPixel,4, 1, 43);
        textMassDalitz2->Draw();
    
        //****************** second Column *************************************************
        TLatex *textMassData3                = new TLatex(columnsLegendMass2[1],rowsLegendMass3[0] ,"Data");
        SetStyleTLatex( textMassData3, textSizeLabelsPixel,4, 1, 43);
        textMassData3->Draw();
        TLatex *textMassMC3                  = new TLatex(columnsLegendMass2[2] ,rowsLegendMass3[0],"MC");
        SetStyleTLatex( textMassMC3, textSizeLabelsPixel,4, 1, 43);
        textMassMC3->Draw();
        
        TLatex *textMassData2               = new TLatex(columnsLegendMass2[4],rowsLegendMass3[0] ,"Data");
        SetStyleTLatex( textMassData2, textSizeLabelsPixel,4, 1, 43);
        textMassData2->Draw();
        TLatex *textMassMC2                 = new TLatex(columnsLegendMass2[5] ,rowsLegendMass3[0],"MC");
        SetStyleTLatex( textMassMC2, textSizeLabelsPixel,4, 1, 43);
        textMassMC2->Draw();
        
        markerPCMPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass3[1]+ offsetMarkerYMass3);
        markerPCMEMCALPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass3[2]+ offsetMarkerYMass3);
        markerEMCALPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass3[3]+ offsetMarkerYMass3);
        TMarker* markerPHOSPi0Mass   = CreateMarkerFromHisto(histoPHOSPi0Mass,columnsLegendMass2[4]+ offsetMarkerXMass2 ,rowsLegendMass3[1]+ offsetMarkerYMass3 ,scaleMarkerMass2);
        markerPHOSPi0Mass->DrawMarker(columnsLegendMass2[4]+ offsetMarkerXMass2 ,rowsLegendMass3[1]+ offsetMarkerYMass3);
        TMarker* markerDalitzPi0Mass   = CreateMarkerFromGraph(graphDalitzPi0Mass,columnsLegendMass2[4]+ offsetMarkerXMass2 ,rowsLegendMass3[2]+ offsetMarkerYMass3 ,scaleMarkerMass2);
        markerDalitzPi0Mass->DrawMarker(columnsLegendMass2[4]+ offsetMarkerXMass2 ,rowsLegendMass3[2]+ offsetMarkerYMass3);
    
        markerPCMPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.01 ,rowsLegendMass3[1]+ offsetMarkerYMass3);
        markerPCMEMCALPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.01 ,rowsLegendMass3[2]+ offsetMarkerYMass3);
        markerEMCALPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.01 ,rowsLegendMass3[3]+ offsetMarkerYMass3);
        TMarker* markerPHOSPi0MassMC = CreateMarkerFromHisto(histoPHOSPi0TrueMass,columnsLegendMass2[5]+ offsetMarkerXMass2 ,rowsLegendMass3[1]+ offsetMarkerYMass3 ,scaleMarkerMass2);
        markerPHOSPi0MassMC->DrawMarker(columnsLegendMass2[5]+ offsetMarkerXMass2-0.01 ,rowsLegendMass3[1]+ offsetMarkerYMass3);
        TMarker* markerDalitzPi0MassMC = CreateMarkerFromGraph(graphDalitzPi0MassMC,columnsLegendMass2[5]+ offsetMarkerXMass2 ,rowsLegendMass3[2]+ offsetMarkerYMass3 ,scaleMarkerMass2);
        markerDalitzPi0MassMC->DrawMarker(columnsLegendMass2[5]+ offsetMarkerXMass2-0.01 ,rowsLegendMass3[2]+ offsetMarkerYMass3);

    canvasMassWidthPi0->Update();
    canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth_incPHOSandDalitz.%s",outputDir.Data(),suffix.Data()));
    
    
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
        Double_t columnsLegendMass3[3]      = {0.,0.57,0.84};
        //****************** first Column **************************************************
        textMassPCM->Draw();
        textMassPCMEMCAL->Draw();
        textMassEMCAL->Draw();
        //****************** second Column *************************************************
        TLatex *textMassData4                = new TLatex(columnsLegendMass3[1],rowsLegendMass2[0] ,"Data");
        SetStyleTLatex( textMassData4, textSizeLabelsPixel,4, 1, 43);
        textMassData4->Draw();
        TLatex *textMassMC4                  = new TLatex(columnsLegendMass3[2] ,rowsLegendMass2[0],"MC");
        SetStyleTLatex( textMassMC4, textSizeLabelsPixel,4, 1, 43);
        textMassMC4->Draw();

        markerPCMPi0Mass->DrawMarker(columnsLegendMass3[1]+ offsetMarkerXMass2*2 ,rowsLegendMass2[1]+ offsetMarkerYMass2);
        markerPCMEMCALPi0Mass->DrawMarker(columnsLegendMass3[1]+ offsetMarkerXMass2*2 ,rowsLegendMass2[2]+ offsetMarkerYMass2);
        markerEMCALPi0Mass->DrawMarker(columnsLegendMass3[1]+ offsetMarkerXMass2*2 ,rowsLegendMass2[3]+ offsetMarkerYMass2);
        markerPCMPi0MassMC->DrawMarker(columnsLegendMass3[2]+ offsetMarkerXMass2*2-0.03 ,rowsLegendMass2[1]+ offsetMarkerYMass2);
        markerPCMEMCALPi0MassMC->DrawMarker(columnsLegendMass3[2]+ offsetMarkerXMass2*2-0.03 ,rowsLegendMass2[2]+ offsetMarkerYMass2);
        markerEMCALPi0MassMC->DrawMarker(columnsLegendMass3[2]+ offsetMarkerXMass2*2-0.03 ,rowsLegendMass2[3]+ offsetMarkerYMass2);
        
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
        
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaMass, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        graphPCMEtaMass->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaMassMC, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
        graphPCMEtaMassMC->Draw("p,same,z");

        DrawGammaLines(0.33, 25. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,0.3, kGray);
        
        labelLegendBMass->Draw();
        
    canvasMassWidthEta->Update();
    canvasMassWidthEta->Print(Form("%s/Eta_MassAndWidth.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************** Cross section for pi0 single measurement 2.76TeV ************************************
    // **********************************************************************************************************************
    
    TCanvas* canvasXSectionPi0  = new TCanvas("canvasXSectionPi0","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasXSectionPi0, 0.16, 0.02, 0.02, 0.08);
    canvasXSectionPi0->SetLogx();
    canvasXSectionPi0->SetLogy();
    
    TH2F * histo2DYieldPi0;
    histo2DYieldPi0          = new TH2F("histo2DYieldPi0","histo2DYieldPi0",11000,0.23,31.,1000,1e-8,10e1);
    SetStyleHistoTH2ForGraphs(histo2DYieldPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.65);
    histo2DYieldPi0->GetXaxis()->SetMoreLogLabels();
    histo2DYieldPi0->GetXaxis()->SetLabelOffset(-0.01);
    histo2DYieldPi0->Draw("copy");


        DrawGammaSetMarkerTGraphAsym(graphPHOSPi0InvYieldSys, markerStyleDet[1], markerSizeDet[1]*0.75, colorDet[1] , colorDet[1], widthLinesBoxes, kTRUE);
        graphPHOSPi0InvYieldSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphDalitzPi0InvYieldSys, markerStyleDet[5], markerSizeDet[5]*0.75, colorDet[5] , colorDet[5], widthLinesBoxes, kTRUE);
        graphDalitzPi0InvYieldSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvYieldSys, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0], widthLinesBoxes, kTRUE);
        graphPCMPi0InvYieldSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0InvYieldSys, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2], widthLinesBoxes, kTRUE);
        graphEMCALPi0InvYieldSys->Draw("E2same");
//      graphPCMEMCALPi0InvYieldSys->RemovePoint(graphPCMEMCALPi0InvYieldSys->GetN()-1);
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0InvYieldSys, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4], widthLinesBoxes, kTRUE);
        graphPCMEMCALPi0InvYieldSys->Draw("E2same");

        DrawGammaSetMarker(histoPHOSPi0InvYieldStat, markerStyleDet[1], markerSizeDet[1]*0.75, colorDet[1] , colorDet[1]);
        histoPHOSPi0InvYieldStat->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphDalitzPi0InvYieldStat,markerStyleDet[5], markerSizeDet[5]*0.75, colorDet[5] , colorDet[5]);
        graphDalitzPi0InvYieldStat->Draw("p,same,z");
//      DrawGammaSetMarker(histoPCMPi0InvYieldStat, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
//      histoPCMPi0InvYieldStat->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphPCMPi0InvYieldStat,markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
        graphPCMPi0InvYieldStat->Draw("p,same,z");
            
        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0InvYieldStat, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
        graphEMCALPi0InvYieldStat->Draw("p,same,z");
//      DrawGammaSetMarker(histoEMCALPi0InvYieldStat, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
//      histoEMCALPi0InvYieldStat->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0InvYieldStat, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4]);
        graphPCMEMCALPi0InvYieldStat->Draw("p,same,z");

        
        TLatex *labelEnergyXSectionPi0      = new TLatex(0.94,0.92,collisionSystempPb.Data());
        SetStyleTLatex( labelEnergyXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);
        labelEnergyXSectionPi0->Draw();
        TLatex *labelDetSysXSectionPi0      = new TLatex(0.94,0.88,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);
        labelDetSysXSectionPi0->Draw();

                                                
        TLegend* legendXSectionPi0          = GetAndSetLegend2(0.20, 0.13,0.40,0.13+(0.035*7), 0.035, 1, "", 42, 0);
        legendXSectionPi0->AddEntry(graphPCMPi0InvYieldSys,nameMeasGlobalLabel[0],"fp");
        legendXSectionPi0->AddEntry(graphPHOSPi0InvYieldSys,nameMeasGlobalLabel[1],"fp");
        legendXSectionPi0->AddEntry(graphEMCALPi0InvYieldSys,nameMeasGlobalLabel[2],"fp");
        legendXSectionPi0->AddEntry(graphPCMEMCALPi0InvYieldSys,nameMeasGlobalLabel[4],"fp");
        legendXSectionPi0->AddEntry(graphDalitzPi0InvYieldSys,nameMeasGlobalLabel[5],"fp");
        legendXSectionPi0->Draw();
   
    canvasXSectionPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems.%s",outputDir.Data(),suffix.Data()));
    
    canvasXSectionPi0->cd();
    histo2DYieldPi0->Draw("copy");

        graphPHOSPi0InvYieldSys->Draw("E2same");
        graphDalitzPi0InvYieldSys->Draw("E2same");
        graphPCMPi0InvYieldSys->Draw("E2same");
        graphEMCALPi0InvYieldSys->Draw("E2same");
        graphPCMEMCALPi0InvYieldSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphCombPi0InvYieldSys->Draw("E2same");
    
        histoPHOSPi0InvYieldStat->Draw("p,same,e");
        graphDalitzPi0InvYieldStat->Draw("p,same,z");
        graphPCMPi0InvYieldStat->Draw("p,same,z");
        graphEMCALPi0InvYieldStat->Draw("p,same,z");
        graphPCMEMCALPi0InvYieldStat->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombPi0InvYieldStat->Draw("p,same,z");

        labelEnergyXSectionPi0->Draw();
        labelDetSysXSectionPi0->Draw();

        legendXSectionPi0->AddEntry(graphCombPi0InvYieldSys,"comb","fp");
        legendXSectionPi0->Draw();

//      DrawGammaSetMarkerTGraphAsym(graphChargedHadronsStatPP, markerStyleDet[1], markerSizeDet[1]*0.75, kBlack , kBlack, widthLinesBoxes);
//      graphChargedHadronsStatPP->Draw("E2same");
        
    canvasXSectionPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems_Comb.%s",outputDir.Data(),suffix.Data()));
  
    TCanvas* canvasXSectionEta      = new TCanvas("canvasXSectionEta","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasXSectionEta, 0.16, 0.02, 0.02, 0.08);
    canvasXSectionEta->SetLogx();
    canvasXSectionEta->SetLogy();
    
    TH2F * histo2DInvYieldEta;
    histo2DInvYieldEta              = new TH2F("histo2DInvYieldEta","histo2DInvYieldEta",11000,0.33,31.,1000,1e-9,1e-0);
    SetStyleHistoTH2ForGraphs(histo2DInvYieldEta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.65);
    histo2DInvYieldEta->GetXaxis()->SetMoreLogLabels();
    histo2DInvYieldEta->GetXaxis()->SetLabelOffset(-0.01);
    histo2DInvYieldEta->Draw("copy");

    
        DrawGammaSetMarkerTGraphAsym(graphPCMEtaInvYieldSys, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0], widthLinesBoxes, kTRUE);
        graphPCMEtaInvYieldSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALEtaInvYieldSys, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4], widthLinesBoxes, kTRUE);
        graphPCMEMCALEtaInvYieldSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphEMCALEtaInvYieldSys, markerStyleDet[2], markerSizeDet[4]*0.75, colorDet[2] , colorDet[2], widthLinesBoxes, kTRUE);
        graphEMCALEtaInvYieldSys->Draw("E2same");
    
        DrawGammaSetMarker(histoPCMEtaInvYieldStat, markerStyleDet[0], markerSizeDet[0]*0.75, colorDet[0] , colorDet[0]);
        histoPCMEtaInvYieldStat->Draw("p,same,e");
        DrawGammaSetMarker(histoPCMEMCALEtaInvYieldStat, markerStyleDet[4], markerSizeDet[4]*0.75, colorDet[4] , colorDet[4]);
        histoPCMEMCALEtaInvYieldStat->Draw("p,same,e");

        DrawGammaSetMarker(histoEMCALEtaInvYieldStat, markerStyleDet[2], markerSizeDet[2]*0.75, colorDet[2] , colorDet[2]);
        histoEMCALEtaInvYieldStat->Draw("p,same,e");


        TLatex *labelEnergyXSectionEta      = new TLatex(0.94,0.92,collisionSystempPb.Data());
        SetStyleTLatex( labelEnergyXSectionEta, 0.035, 4, 1, 42, kTRUE, 31);
        labelEnergyXSectionEta->Draw();
        TLatex *labelDetSysXSectionEta      = new TLatex(0.94,0.88,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionEta, 0.035, 4, 1, 42, kTRUE, 31);
        labelDetSysXSectionEta->Draw();

        TLegend* legendXSectionEta          = GetAndSetLegend2(0.20, 0.13,0.40,0.13+(0.035*4), 0.035, 1, "", 42, 0);
        legendXSectionEta->AddEntry(graphPCMEtaInvYieldSys,nameMeasGlobalLabel[0],"fp");
//      legendXSectionEta->AddEntry(graphPHOSPi0InvYieldSys,nameMeasGlobalLabel[1],"fp");
        legendXSectionEta->AddEntry(graphEMCALEtaInvYieldSys,nameMeasGlobalLabel[2],"fp");
        legendXSectionEta->AddEntry(graphPCMEMCALEtaInvYieldSys,nameMeasGlobalLabel[4],"fp");
        legendXSectionEta->Draw();
   
    canvasXSectionEta->SaveAs(Form("%s/Eta_InvYieldCompAllSystems.%s",outputDir.Data(),suffix.Data()));
    histo2DInvYieldEta->Draw("copy");

        graphPCMEtaInvYieldSys->Draw("E2same");
        graphEMCALEtaInvYieldSys->Draw("E2same");
        graphPCMEMCALEtaInvYieldSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphCombEtaInvYieldSys->Draw("E2same");
    
        histoPCMEtaInvYieldStat->Draw("p,same,e");
        histoEMCALEtaInvYieldStat->Draw("p,same,e");
        histoPCMEMCALEtaInvYieldStat->Draw("p,same,e");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombEtaInvYieldStat->Draw("p,same,z");


        labelEnergyXSectionEta->Draw();
        labelDetSysXSectionEta->Draw();

        legendXSectionEta->AddEntry(graphCombEtaInvYieldSys,"comb","fp");
        legendXSectionEta->Draw();
        
    canvasXSectionEta->SaveAs(Form("%s/Eta_InvYieldCompAllSystems_Comb.%s",outputDir.Data(),suffix.Data()));


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
    
        TH1F * histo1DAccEff;
        histo1DAccEff                = new TH1F("histo1DAccEff", "histo1DAccEff",1000, 0.23,  31);
        SetStyleHistoTH1ForGraphs( histo1DAccEff, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec}#upoint#it{BR}"),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
        histo1DAccEff->GetYaxis()->SetRangeUser(1e-5, 2e-0 );
        histo1DAccEff->GetYaxis()->SetLabelOffset(0.001);
        histo1DAccEff->GetXaxis()->SetLabelOffset(-0.01);
        histo1DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo1DAccEff->DrawCopy(); 

        DrawGammaSetMarkerTGraphAsym(graphPCMPi0AccTimesEff, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
        graphPCMPi0AccTimesEff->Draw("p,same,z");
        
        DrawGammaSetMarkerTGraphAsym(graphPCMEMCALPi0AccTimesEff, markerStyleDet[4], markerSizeDet[4]*0.55, colorDet[4] , colorDet[4]);
        graphPCMEMCALPi0AccTimesEff->Draw("p,same,z");

        DrawGammaSetMarkerTGraphAsym(graphEMCALPi0AccTimesEff, markerStyleDet[2], markerSizeDet[2]*0.55, colorDet[2] , colorDet[2]);
        graphEMCALPi0AccTimesEff->Draw("p,same,z");
        
        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(4*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffiAccPi0->AddEntry(graphPCMPi0AccTimesEff,nameMeasGlobalLabel[0],"p");
        legendEffiAccPi0->AddEntry(graphPCMEMCALPi0AccTimesEff,nameMeasGlobalLabel[4],"p");
        legendEffiAccPi0->AddEntry(graphEMCALPi0AccTimesEff,nameMeasGlobalLabel[2],"p");
        legendEffiAccPi0->Draw();

        TLatex *labelPerfEffi               = new TLatex(0.15,0.92,"ALICE performance");
        SetStyleTLatex( labelPerfEffi, textSizeLabelsRel,4);
        labelPerfEffi->Draw();
        TLatex *labelEnergyEffi             = new TLatex(0.15,0.87,collisionSystempPb.Data());
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
        
        histo1DAccEff->DrawCopy();
        
//         DrawGammaSetMarkerTGraphAsym(graphPHOSPi0AccTimesEff, markerStyleDet[1], markerSizeDet[1]*0.55, colorDet[1] , colorDet[1]);
//         graphPHOSPi0AccTimesEff->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphDalitzPi0AccTimesEff, markerStyleDet[5], markerSizeDet[5]*0.55, colorDet[5] , colorDet[5]);
        graphDalitzPi0AccTimesEff->Draw("p,same,z");

        graphPCMPi0AccTimesEff->Draw("p,same,z");
        graphPCMEMCALPi0AccTimesEff->Draw("p,same,z");
        graphEMCALPi0AccTimesEff->Draw("p,same,z");

        TLegend* legendEffiAccPi02          = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(5*textSizeLabelsRel),textSizeLabelsPixel);
        legendEffiAccPi02->AddEntry(graphPCMPi0AccTimesEff,nameMeasGlobalLabel[0],"p");
        legendEffiAccPi02->AddEntry(graphPCMEMCALPi0AccTimesEff,nameMeasGlobalLabel[4],"p");
        legendEffiAccPi02->AddEntry(graphEMCALPi0AccTimesEff,nameMeasGlobalLabel[2],"p");
//         legendEffiAccPi02->AddEntry(graphPHOSPi0AccTimesEff,nameMeasGlobalLabel[1],"p");
        legendEffiAccPi02->AddEntry(graphDalitzPi0AccTimesEff,nameMeasGlobalLabel[5],"p");
        legendEffiAccPi02->Draw();

        labelPerfEffi->Draw();
        labelEnergyEffi->Draw();
        labelDetSysEffiPi0->Draw();
        
    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff_incPHOSAndDalitz.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************** Acceptance * Efficiency for eta single measurement **********************************
    // **********************************************************************************************************************
    canvasAcceptanceTimesEff->cd();
        
        TH1F * histo1DAccEffEta;
        histo1DAccEffEta                = new TH1F("histo1DAccEffEta", "histo1DAccEffEta",1000, 0.33, 31);
        SetStyleHistoTH1ForGraphs( histo1DAccEffEta, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec}"),  
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
        histo1DAccEffEta->GetYaxis()->SetLabelOffset(0.001);
        histo1DAccEffEta->GetXaxis()->SetLabelOffset(-0.01);
        histo1DAccEffEta->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo1DAccEffEta->GetYaxis()->SetRangeUser(8e-4, 2.0e-0 );
        histo1DAccEffEta->DrawCopy(); 
        
        
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

// //     
// //     // **********************************************************************************************************************
// //     // ******************************************* Comparison to theory calculations Pi0 ************************************
// //     // **********************************************************************************************************************    
// //     textSizeLabelsPixel                     = 48;
// //     
// //     TCanvas* canvasRatioPP                  = new TCanvas("canvasRatioPP","",200,10,1350,900);  // gives the page size
// //     DrawGammaCanvasSettings( canvasRatioPP,  0.12, 0.01, 0.01, 0.11);
// //     canvasRatioPP->cd();
// //     canvasRatioPP->SetLogx();
// // 
// //         textsizeLabelsPP                    = 0;
// //         textsizeFacPP                       = 0;
// //         if (canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) <canvasRatioPP->YtoPixel(canvasRatioPP->GetY1()) ){
// //             textsizeLabelsPP                = (Double_t)textSizeLabelsPixel/canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) ;
// //             textsizeFacPP                   = (Double_t)1./canvasRatioPP->XtoPixel(canvasRatioPP->GetX2()) ;
// //         } else {
// //         textsizeLabelsPP                    = (Double_t)textSizeLabelsPixel/canvasRatioPP->YtoPixel(canvasRatioPP->GetY1());
// //         textsizeFacPP                       = (Double_t)1./canvasRatioPP->YtoPixel(canvasRatioPP->GetY1());
// //         }
// //         cout << textsizeLabelsPP << endl;
// // 
// //         TH2F * ratio2DTheoryPP       = new TH2F("ratio2DTheoryPP","ratio2DTheoryPP",1000,0.23,31.,1000,0.2,3.6);
// //         SetStyleHistoTH2ForGraphs(ratio2DTheoryPP, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
// //                                   0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
// //         ratio2DTheoryPP->GetYaxis()->SetMoreLogLabels(kTRUE);
// //         ratio2DTheoryPP->GetYaxis()->SetNdivisions(505);
// //         ratio2DTheoryPP->GetYaxis()->SetNoExponent(kTRUE);
// //         ratio2DTheoryPP->GetXaxis()->SetMoreLogLabels(kTRUE);
// //         ratio2DTheoryPP->GetXaxis()->SetNoExponent(kTRUE);
// //         ratio2DTheoryPP->GetXaxis()->SetLabelFont(42);
// //         ratio2DTheoryPP->GetYaxis()->SetLabelFont(42);
// //         ratio2DTheoryPP->DrawCopy(); 
// // 
// //         graphRatioPi0CombNLODSS14->RemovePoint(0);
// //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS14, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
// //         graphRatioPi0CombNLODSS14->Draw("3,same");
// // 
// // //         DrawGammaNLOTGraph( graphRatioPi0CombCGC, widthCommonFit, styleLineCGC, colorCGC);
// // //         graphRatioPi0CombCGC->Draw("same,c");
// //         
// //         DrawGammaSetMarker(histoRatioPi0Pythia8ToFit, 24, 1.5, kRed+2 , kRed+2);  
// //         histoRatioPi0Pythia8ToFit->SetLineWidth(widthCommonFit);
// // //      histoRatioPi0Pythia8ToFit->GetXaxis()->SetRangeUser(0.5,14);
// //         histoRatioPi0Pythia8ToFit->Draw("same,hist,l");
// //         DrawGammaNLOTGraph( graphRatioPi0CombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, kGray+1);
// //         graphRatioPi0CombNLOMuHalf->Draw("same,c");
// //         DrawGammaNLOTGraph( graphRatioPi0CombNLOMuOne, widthCommonFit, styleLineNLOMuOne, kGray+1);
// //         graphRatioPi0CombNLOMuOne->Draw("same,c");
// //         DrawGammaNLOTGraph( graphRatioPi0CombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, kGray+2);
// //         graphRatioPi0CombNLOMuTwo->Draw("same,c");
// // 
// // 
// //         TGraphAsymmErrors* graphRatioPi0CombCombFitStatWOXErr = (TGraphAsymmErrors*)graphRatioPi0CombCombFitStat->Clone("graphRatioPi0CombCombFitStatWOXErr");
// //         ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatWOXErr);
// //         
// //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
// //         graphRatioPi0CombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
// //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
// //         graphRatioPi0CombCombFitSys->SetLineWidth(0);
// //         graphRatioPi0CombCombFitSys->Draw("2,same");
// //         graphRatioPi0CombCombFitStatWOXErr->Draw("p,same");
// // 
// //         TBox* boxErrorSigma2760GeVRatioPi0 = CreateBoxConv(kGray+1, 0.28, 1.-(0.0544 ), 0.32, 1.+(0.0544));
// //         boxErrorSigma2760GeVRatioPi0->Draw();
// //         DrawGammaLines(0.23, 70.,1., 1.,0.1,kGray);
// // 
// //         // upper left corner
// //         TLatex *labelRatioTheoryNLO   = new TLatex(0.155,0.93,"NLO, PDF: CTEQ6M5 FF: DSS07");
// //         SetStyleTLatex( labelRatioTheoryNLO, 0.85*textsizeLabelsPP,4);
// //         labelRatioTheoryNLO->Draw();
// //         TLegend* legendRatioPi0OldTheo= GetAndSetLegend2(0.15,0.91-3*0.9*textsizeLabelsPP,0.4,0.91, 0.85* textSizeLabelsPixel, 1, "", 43, 0.2);
// //         legendRatioPi0OldTheo->AddEntry(graphRatioPi0CombNLOMuHalf, "#mu = 0.5 #it{p}_{T}", "l");  
// //         legendRatioPi0OldTheo->AddEntry(graphRatioPi0CombNLOMuOne,  "#mu = #it{p}_{T}", "l");  
// //         legendRatioPi0OldTheo->AddEntry(graphRatioPi0CombNLOMuTwo,  "#mu = 2 #it{p}_{T}", "l");  
// //         legendRatioPi0OldTheo->Draw();
// //         
// //         // lower left corner
// //         TLegend* legendRatioPi0OldTheo3= GetAndSetLegend2(0.15,0.14,0.4,0.15+2*0.85*textsizeLabelsPP, 0.85* textSizeLabelsPixel, 1, "", 43, 0.2);
// //         legendRatioPi0OldTheo3->AddEntry(graphRatioPi0CombNLODSS14,  "NLO, PDF: MSTW FF: DSS14", "f");  
// //         legendRatioPi0OldTheo3->AddEntry(histoRatioPi0Pythia8ToFit,  "PYTHIA 8.2, Monash 2013", "l");
// //         legendRatioPi0OldTheo3->Draw();
// // 
// //         // upper right corner
// //         TLegend* legendRatioPi0OldTheo2= GetAndSetLegend2(0.72,0.91-1*0.85*textsizeLabelsPP,0.95,0.91, 0.85* textSizeLabelsPixel, 1, "", 43, 0.27);
// //         legendRatioPi0OldTheo2->AddEntry(graphRatioPi0CombCombFitSys,"#pi^{0} ALICE","pf");
// //         legendRatioPi0OldTheo2->Draw();
// //         
// //         TLatex *labelRatioTheory2   = new TLatex(0.73,0.93,collisionSystempPb.Data());
// //         SetStyleTLatex( labelRatioTheory2, 0.85*textsizeLabelsPP,4);
// //         labelRatioTheory2->Draw();
// //     
// //         
// //     canvasRatioPP->Update();
// //     canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PP.%s",outputDir.Data(),suffix.Data()));
// // 
// //     canvasRatioPP->cd();
// //         ratio2DTheoryPP->GetYaxis()->SetNdivisions(512);
// //         ratio2DTheoryPP->GetYaxis()->SetRangeUser(0.5,1.9);
// //         ratio2DTheoryPP->DrawCopy(); 
// // 
// //         
// //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS14, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
// //         graphRatioPi0CombNLODSS14->Draw("3,same");
// //         
// //         DrawGammaSetMarker(histoRatioPi0Pythia8ToFit, 24, 1.5, kRed+2 , kRed+2);  
// //         histoRatioPi0Pythia8ToFit->SetLineWidth(widthCommonFit);
// // //      histoRatioPi0Pythia8ToFit->GetXaxis()->SetRangeUser(0.5,14);
// //         histoRatioPi0Pythia8ToFit->Draw("same,hist,l");
// // 
// // //         DrawGammaNLOTGraph( graphRatioPi0CombCGC, widthCommonFit, styleLineCGC, colorCGC);
// // //         graphRatioPi0CombCGC->Draw("same,c");
// //         
// //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
// //         graphRatioPi0CombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
// //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
// //         graphRatioPi0CombCombFitSys->SetLineWidth(0);
// //         graphRatioPi0CombCombFitSys->Draw("2,same");
// //         graphRatioPi0CombCombFitStatWOXErr->Draw("p,same");
// // 
// //         boxErrorSigma2760GeVRatioPi0->Draw();
// //         DrawGammaLines(0.23, 70.,1., 1.,0.1,kGray);
// // 
// //         TLegend* legendRatioTheoryPPNew= GetAndSetLegend2(0.15,0.92-0.85*textsizeLabelsPP*3,0.4,0.92, 0.85*textSizeLabelsPixel);
// //         legendRatioTheoryPPNew->AddEntry(graphRatioPi0CombCombFitSys,"#pi^{0} ALICE","pf");
// //         legendRatioTheoryPPNew->AddEntry(graphRatioPi0CombNLODSS14,"NLO, PDF: MSTW FF: DSS14","f");
// // //         legendRatioTheoryPPNew->AddEntry((TObject*)0,"Phys.Rev. D91 no. 1, (2015) 014035","");
// //         legendRatioTheoryPPNew->AddEntry(histoRatioPi0Pythia8ToFit,  "PYTHIA 8.2, Monash 2013", "l");  
// //         legendRatioTheoryPPNew->Draw();
// // 
        TLatex *labelRatioTheory   = new TLatex(0.15,0.925,collisionSystempPb.Data());
        SetStyleTLatex( labelRatioTheory, 0.85*textsizeLabelsPP,4);        
        labelRatioTheory->Draw();
// //         
// //     canvasRatioPP->Update();
// //     canvasRatioPP->Print(Form("%s/Pi0_RatioTheoryToData_PP_NewTheory.%s",outputDir.Data(),suffix.Data()));
// // 
// //     // **********************************************************************************************************************
// //     // ******************************************* Comparison to theory calculations Eta ************************************
// //     // **********************************************************************************************************************    
// //     canvasRatioPP->cd();
// //     canvasRatioPP->SetLogx();
// // 
// //         TH2F * ratio2DTheoryPPEta       = new TH2F("ratio2DTheoryPPEta","ratio2DTheoryPPEta",1000,0.33,25.,1000,0.2,3.6);
// //         SetStyleHistoTH2ForGraphs(ratio2DTheoryPPEta, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
// //                                   0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
// //         ratio2DTheoryPPEta->GetYaxis()->SetMoreLogLabels(kTRUE);
// //         ratio2DTheoryPPEta->GetYaxis()->SetNdivisions(505);
// //         ratio2DTheoryPPEta->GetYaxis()->SetNoExponent(kTRUE);
// //         ratio2DTheoryPPEta->GetXaxis()->SetMoreLogLabels(kTRUE);
// //         ratio2DTheoryPPEta->GetXaxis()->SetNoExponent(kTRUE);
// //         ratio2DTheoryPPEta->GetXaxis()->SetLabelFont(42);
// //         ratio2DTheoryPPEta->GetYaxis()->SetLabelFont(42);
// //         ratio2DTheoryPPEta->DrawCopy(); 
// // 
// //         DrawGammaNLOTGraph( graphRatioEtaCombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
// //         graphRatioEtaCombNLOMuHalf->Draw("same,c");
// //         DrawGammaNLOTGraph( graphRatioEtaCombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
// //         graphRatioEtaCombNLOMuOne->Draw("same,c");
// //         DrawGammaNLOTGraph( graphRatioEtaCombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
// //         graphRatioEtaCombNLOMuTwo->Draw("same,c");
// //         
// //         DrawGammaSetMarker(histoRatioEtaPythia8ToFit, 24, 1.5, kRed+2 , kRed+2);  
// //         histoRatioEtaPythia8ToFit->SetLineWidth(widthCommonFit);
// // //      histoRatioEtaPythia8ToFit->GetXaxis()->SetRangeUser(0.5,14);
// //         histoRatioEtaPythia8ToFit->Draw("same,hist,l");
// // 
// //         TGraphAsymmErrors* graphRatioEtaCombCombFitStatWOXErr = (TGraphAsymmErrors*)graphRatioEtaCombCombFitStat->Clone("graphRatioEtaCombCombFitStatWOXErr");
// //         ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombCombFitStatWOXErr);
// //         
// //         DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
// //         graphRatioEtaCombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
// //         DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
// //         graphRatioEtaCombCombFitSys->SetLineWidth(0);
// //         graphRatioEtaCombCombFitSys->Draw("2,same");
// //         graphRatioEtaCombCombFitStatWOXErr->Draw("p,same");
// // 
// //         TBox* boxErrorSigma2760GeVRatioEta = CreateBoxConv(kGray+1, 0.40, 1.-(0.0544), 0.44, 1.+(0.0544));        
// //         boxErrorSigma2760GeVRatioEta->Draw();
// //         DrawGammaLines(0.33, 25.,1., 1.,0.1,kGray);
// // 
// //         // upper left corner
// //         TLatex *labelRatioTheoryNLOEta   = new TLatex(0.155,0.93,"NLO, PDF: CTEQ6M5 FF: AESSS");
// //         SetStyleTLatex( labelRatioTheoryNLOEta, 0.85*textsizeLabelsPP,4);
// //         labelRatioTheoryNLOEta->Draw();
// //         legendRatioPi0OldTheo->Draw();
// //         
// //         // lower left corner
// //         TLegend* legendRatioEtaOldTheo3= GetAndSetLegend2(0.15,0.14,0.4,0.15+1*0.85*textsizeLabelsPP, 0.85* textSizeLabelsPixel, 1, "", 43, 0.2);
// //         legendRatioEtaOldTheo3->AddEntry(histoRatioPi0Pythia8ToFit,  "PYTHIA 8.2, Monash 2013", "l");  
// //         legendRatioEtaOldTheo3->Draw();
// // 
// //         // upper right corner
// //         TLegend* legendRatioEtaOldTheo2= GetAndSetLegend2(0.72,0.91-1*0.85*textsizeLabelsPP,0.95,0.91, 0.85* textSizeLabelsPixel, 1, "", 43, 0.27);
// //         legendRatioEtaOldTheo2->AddEntry(graphRatioPi0CombCombFitSys,"#eta ALICE","pf");
// //         legendRatioEtaOldTheo2->Draw();
// //         labelRatioTheory2->Draw();
// //             
// //     canvasRatioPP->Update();
// //     canvasRatioPP->Print(Form("%s/Eta_RatioTheoryToData_PP.%s",outputDir.Data(),suffix.Data()));
// //     
// //     //*************************************************************************************************************
// //     //***************************** Paper plot X-section and ratios ***********************************************
// //     //*************************************************************************************************************
// //     
// //     Double_t arrayBoundariesX1_XSec[2];
// //     Double_t arrayBoundariesY1_XSec[6];
// //     Double_t relativeMarginsXXSec[3];
// //     Double_t relativeMarginsYXSec[3];
// //     textSizeLabelsPixel = 48;
// //     ReturnCorrectValuesForCanvasScaling(1250,2000, 1, 5,0.135, 0.005, 0.003,0.05,arrayBoundariesX1_XSec,arrayBoundariesY1_XSec,relativeMarginsXXSec,relativeMarginsYXSec);
// //     
// //     TCanvas* canvasInvSectionPaper      = new TCanvas("canvasInvSectionPaper","",0,0,1250,2000);  // gives the page size
// //     DrawGammaCanvasSettings( canvasInvSectionPaper,  0.13, 0.02, 0.03, 0.06);
// // 
// //     TPad* padInvSectionSpec             = new TPad("padInvSectionSpec", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[3], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[0],-1, -1, -2);
// //     DrawGammaPadSettings( padInvSectionSpec, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[0], relativeMarginsYXSec[1]);
// //     padInvSectionSpec->Draw();
// //     Double_t marginXSec                 = relativeMarginsXXSec[0]*1250;
// //     Double_t textsizeLabelsXSecUp       = 0;
// //     Double_t textsizeFacXSecUp          = 0;
// //     if (padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) < padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1())){
// //         textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
// //         textsizeFacXSecUp               = (Double_t)1./padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
// //     } else {
// //         textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
// //         textsizeFacXSecUp                   = (Double_t)1./padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
// //     }
// //     
// //     TPad* padInvSectionNLORatio         = new TPad("padInvSectionNLORatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[4], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[3],-1, -1, -2);
// //     DrawGammaPadSettings( padInvSectionNLORatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[1]);
// //     padInvSectionNLORatio->Draw();
// //     Double_t textsizeLabelsXSecMiddle   = 0;
// //     Double_t textsizeFacXSecMiddle      = 0;
// //     if (padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) < padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1())){
// //         textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) ;
// //         textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) ;
// //     } else {
// //         textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1());
// //         textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1());
// //     }
// //     
// //     TPad* padInvSectionPythiaRatio      = new TPad("padInvSectionPythiaRatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[5], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[4],-1, -1, -2);
// //     DrawGammaPadSettings( padInvSectionPythiaRatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[2]);
// //     padInvSectionPythiaRatio->Draw();
// //     Double_t textsizeLabelsXSecDown     = 0;
// //     Double_t textsizeFacXSecDown        = 0;
// //     if (padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) < padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1())){
// //         textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) ;
// //         textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) ;
// //     } else {
// //         textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1());
// //         textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1());
// //     }
// // 
// //     
// //     padInvSectionSpec->cd();
// //     padInvSectionSpec->SetLogy(1);
// //     padInvSectionSpec->SetLogx(1);
// //         SetStyleHistoTH2ForGraphs(histo2DYieldPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
// //                                 0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.2/(textsizeFacXSecUp*marginXSec));
// //         histo2DYieldPi0->GetXaxis()->SetMoreLogLabels();
// //         histo2DYieldPi0->GetXaxis()->SetLabelOffset(+0.01);
// //         histo2DYieldPi0->Draw();
// // 
// //         
// //         DrawGammaSetMarker(histoPythia8InvYieldPi0, 24, 1.5, kRed+2 , kRed+2);  
// //         histoPythia8InvYieldPi0->SetLineWidth(widthCommonFit);
// //         histoPythia8InvYieldPi0->Draw("same,hist,l");
// // 
// //         graphNLODSS14Calc->RemovePoint(0);
// //         DrawGammaSetMarkerTGraphAsym(graphNLODSS14Calc, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
// //         graphNLODSS14Calc->Draw("3,same");
// // 
// // //         DrawGammaNLOTGraph( graphCGCPi0, widthCommonFit, styleLineCGC, colorCGC);
// // //         graphCGCPi0->Draw("same,c");
// // 
// // 
        TGraphAsymmErrors* graphCombPi0InvYieldStatWOXErr = (TGraphAsymmErrors*)graphCombPi0InvYieldStat->Clone("graphCombPi0InvYieldStatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphCombPi0InvYieldStatWOXErr);
        
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombPi0InvYieldSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
        graphCombPi0InvYieldStatWOXErr->Draw("p,same,z");

        DrawGammaSetMarkerTF1( fitTCMInvYieldPi0, 7, 2, kGray+2); 
        fitTCMInvYieldPi0->Draw("same");        

        DrawGammaSetMarkerTF1( fitPowInvYieldPi0, 7, 2, kBlue-7); 
        fitPowInvYieldPi0->Draw("same");        
// //         
// //         
// //         TLatex *labelEnergyXSectionPaper= new TLatex(0.95, 0.91, collisionSystempPb.Data());
// //         SetStyleTLatex( labelEnergyXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
// //         labelEnergyXSectionPaper->Draw();
// //         TLatex *labelALICEXSectionPaper= new TLatex(0.95,0.87,"ALICE");
// //         SetStyleTLatex( labelALICEXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
// //         labelALICEXSectionPaper->Draw();
// //         TLatex *labelDetSysXSectionPaper= new TLatex(0.95,0.83,"#pi^{0} #rightarrow #gamma#gamma");
// //         SetStyleTLatex( labelDetSysXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
// //         labelDetSysXSectionPaper->Draw();
// //         
// //         TLegend* legendXsectionPaper    = GetAndSetLegend2(0.17, 0.14, 0.5, 0.14+0.05*4, textSizeLabelsPixel, 1, "", 43, 0.2);
// //         legendXsectionPaper->AddEntry(graphCombPi0InvYieldSys,"Data,","pf");
// //         legendXsectionPaper->AddEntry((TObject*)0,"Norm. unc. 5.4%","");
// //         legendXsectionPaper->AddEntry(graphNLODSS14Calc,"NLO, PDF: MSTW, FF: DSS14","f");
// // //         legendXsectionPaper->AddEntry((TObject*)0,"Phys.Rev. D91 no. 1, (2015) 014035","");
// //         legendXsectionPaper->AddEntry(histoPythia8InvYieldPi0,"PYTHIA 8.2, Monash 2013","l");
// //         legendXsectionPaper->Draw();
// // 
// //         TLegend* legendXsectionPaper2     = GetAndSetLegend2(0.17, 0.05, 0.5, 0.11, textSizeLabelsPixel, 1, "", 43, 0.2);
// //         legendXsectionPaper2->AddEntry(fitTCMInvYieldPi0,"#it{A}_{e} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{n}","l");
// //         legendXsectionPaper2->Draw();
// //         
// //         
// //     padInvSectionNLORatio->cd();
// //     padInvSectionNLORatio->SetLogx(1);
// //         TH2F * ratio2DNLO               = new TH2F("ratio2DNLO","ratio2DNLO",1000,0.23,31.,1000,0.5,1.9);
// //         SetStyleHistoTH2ForGraphs(ratio2DNLO, "#it{p}_{T} (GeV/#it{c})","#frac{NLO, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle, 
// //                                   0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 512);
// // //         ratio2DNLO->GetYaxis()->SetNdivisions(512);
// //         ratio2DNLO->GetYaxis()->SetNoExponent(kTRUE);
// //         ratio2DNLO->GetXaxis()->SetMoreLogLabels(kTRUE);
// //         ratio2DNLO->GetXaxis()->SetNoExponent(kTRUE);
// //         ratio2DNLO->GetXaxis()->SetLabelFont(42);
// //         ratio2DNLO->GetYaxis()->SetLabelFont(42);
// //         ratio2DNLO->GetYaxis()->SetLabelOffset(+0.01);
// //         ratio2DNLO->GetXaxis()->SetTickLength(0.07);
// //         ratio2DNLO->DrawCopy();
// // 
// // //         DrawGammaNLOTGraph( graphRatioPi0CombCGC, widthCommonFit, styleLineCGC, colorCGC);
// // //         graphRatioPi0CombCGC->Draw("same,c");
// // 
// //         boxErrorSigma2760GeVRatioPi0->Draw();
// //         
// //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS14, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
// //         graphRatioPi0CombNLODSS14->Draw("3,same");
// // 
// //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
// //         graphRatioPi0CombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
// //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
// //         graphRatioPi0CombCombFitSys->SetLineWidth(0);
// //         graphRatioPi0CombCombFitSys->Draw("2,same");
// //         graphRatioPi0CombCombFitStatWOXErr->Draw("p,same");
// //         
// //         DrawGammaLines(0.23, 70.,1., 1.,0.1,kGray);
// //         
// //     padInvSectionPythiaRatio->cd();
// //     padInvSectionPythiaRatio->SetLogx(1);
// //         TH2F * ratio2DPythia            = new TH2F("ratio2DPythia","ratio2DPythia",1000,0.23,31.,1000,0.5,1.9);
// //         SetStyleHistoTH2ForGraphs(ratio2DPythia, "#it{p}_{T} (GeV/#it{c})","#frac{PYTHIA, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown, 
// //                                   0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.2/(textsizeFacXSecDown*marginXSec), 510, 512);
// // //         ratio2DPythia->GetYaxis()->SetNdivisions(512);
// //         ratio2DPythia->GetYaxis()->SetNoExponent(kTRUE);
// //         ratio2DPythia->GetXaxis()->SetMoreLogLabels(kTRUE);
// //         ratio2DPythia->GetXaxis()->SetNoExponent(kTRUE);
// //         ratio2DPythia->GetXaxis()->SetLabelFont(42);
// //         ratio2DPythia->GetYaxis()->SetLabelFont(42);
// //         ratio2DPythia->GetYaxis()->SetLabelOffset(+0.01);
// //         ratio2DPythia->GetXaxis()->SetTickLength(0.06);
// //         ratio2DPythia->GetYaxis()->SetTickLength(0.04);
// //         ratio2DPythia->DrawCopy();
// // 
// //         DrawGammaSetMarker(histoRatioPi0Pythia8ToFit, 24, 1.5, kRed+2 , kRed+2);  
// //         histoRatioPi0Pythia8ToFit->SetLineWidth(widthCommonFit);
// // //      histoRatioPi0Pythia8ToFit->GetXaxis()->SetRangeUser(0.5,14);
// //         histoRatioPi0Pythia8ToFit->Draw("same,hist,l");
// //         
// //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
// //         graphRatioPi0CombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
// //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
// //         graphRatioPi0CombCombFitSys->SetLineWidth(0);
// //         graphRatioPi0CombCombFitSys->Draw("2,same");
// //         graphRatioPi0CombCombFitStatWOXErr->Draw("p,same");
// //         
// //         boxErrorSigma2760GeVRatioPi0->Draw();
// //         
// //         DrawGammaLines(0.23, 70.,1., 1.,0.1,kGray);
// //         
// //     canvasInvSectionPaper->Print(Form("%s/Pi0_InvYieldWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));
// //     
// //     padInvSectionSpec->cd();
// //     padInvSectionSpec->SetLogy(1);
// //     padInvSectionSpec->SetLogx(1);
// //         SetStyleHistoTH2ForGraphs(histo2DInvYieldEta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
// //                                 0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.2/(textsizeFacXSecUp*marginXSec));
// //         histo2DInvYieldEta->GetXaxis()->SetMoreLogLabels();
// //         histo2DInvYieldEta->GetXaxis()->SetLabelOffset(+0.01);
// //         histo2DInvYieldEta->Draw();
// // 
// //         DrawGammaSetMarker(histoPythia8InvYieldEta, 24, 1.5, kRed+2 , kRed+2);  
// //         histoPythia8InvYieldEta->SetLineWidth(widthCommonFit);
// //         histoPythia8InvYieldEta->Draw("same,hist,l");
// // 
// //         DrawGammaNLOTGraph( graphNLOCalcEtaMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
// //         graphNLOCalcEtaMuHalf->Draw("same,c");
// //         DrawGammaNLOTGraph( graphNLOCalcEtaMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
// //         graphNLOCalcEtaMuOne->Draw("same,c");
// //         DrawGammaNLOTGraph( graphNLOCalcEtaMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
// //         graphNLOCalcEtaMuTwo->Draw("same,c");
// // 
        TGraphAsymmErrors* graphCombEtaInvYieldStatWOXErr = (TGraphAsymmErrors*)graphCombEtaInvYieldStat->Clone("graphCombEtaInvYieldStatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphCombEtaInvYieldStatWOXErr);
// // 
// //         
// //         DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
// //         graphCombEtaInvYieldSys->Draw("E2same");
// //         DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
// //         graphCombEtaInvYieldStatWOXErr->Draw("p,same,z");
// // 
// //         DrawGammaSetMarkerTF1( fitTCMInvYieldEta, 7, 2, kGray+2); 
// //         fitTCMInvYieldEta->Draw("same");
// // 
// //         // labels lower left corner
// //         TLegend* legendXsectionPaperEta    = GetAndSetLegend2(0.17, 0.11, 0.5, 0.11+0.05*2, textSizeLabelsPixel, 1, "", 43, 0.2);
// //         legendXsectionPaperEta->AddEntry(graphCombPi0InvYieldSys,"Data,","pf");
// //         legendXsectionPaperEta->AddEntry((TObject*)0,"Norm. unc. 5.4%","");
// //         legendXsectionPaperEta->Draw();
// //         legendXsectionPaper2->Draw();
// //         
// //         TLatex *labelEnergyXSectionPaperEta= new TLatex(0.18, 0.13+0.05*4, collisionSystempPb.Data());
// //         SetStyleTLatex( labelEnergyXSectionPaperEta, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 11);
// //         labelEnergyXSectionPaperEta->Draw();
// //         TLatex *labelALICEXSectionPaperEta= new TLatex(0.18,0.13+0.05*3,"ALICE");
// //         SetStyleTLatex( labelALICEXSectionPaperEta, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 11);
// //         labelALICEXSectionPaperEta->Draw();
// //         TLatex *labelDetSysXSectionPaperEta = new TLatex(0.18,0.13+0.05*2,"#eta #rightarrow #gamma#gamma");
// //         SetStyleTLatex( labelDetSysXSectionPaperEta, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 11);
// //         labelDetSysXSectionPaperEta->Draw();
// //         
// // 
// //         // labels upper right corner
// //         TLegend* legendXsectionPaperEtaTheo     = GetAndSetLegend2(0.59, 0.95-0.04*2, 0.59+0.33, 0.95, textSizeLabelsPixel, 1, "", 43, 0.2);
// //         legendXsectionPaperEtaTheo->AddEntry(histoPythia8InvYieldEta,"PYTHIA 8.2,","l");
// //         legendXsectionPaperEtaTheo->AddEntry((TObject*)0,"Monash 2013","");
// //         legendXsectionPaperEtaTheo->Draw();
// // 
// //         TLatex *labelNLOHeaderEta= new TLatex(0.60, 0.805, "NLO, PDF: CTEQ6M5");
// //         SetStyleTLatex( labelNLOHeaderEta, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 11);
// //         labelNLOHeaderEta->Draw();
// //         TLatex *labelNLOHeaderEta2= new TLatex(0.695, 0.765, "FF: AESSS");
// //         SetStyleTLatex( labelNLOHeaderEta2, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 11);
// //         labelNLOHeaderEta2->Draw();
// //                 
// //         TLegend* legendXsectionPaperEtaTheo2     = GetAndSetLegend2(0.73, 0.755-0.048*3, 0.95, 0.755, textSizeLabelsPixel, 1, "", 43, 0.27);
// //         legendXsectionPaperEtaTheo2->AddEntry(graphNLOCalcEtaMuHalf,"#mu = 0.5 #it{p}_{T}","l");
// //         legendXsectionPaperEtaTheo2->AddEntry(graphNLOCalcEtaMuOne,"#mu = #it{p}_{T}","l");
// // //         legendXsectionPaperEtaTheo2->AddEntry((TObject*)0,"","");
// //         legendXsectionPaperEtaTheo2->AddEntry(graphNLOCalcEtaMuTwo,"#mu = 2 #it{p}_{T}","l");
// //         legendXsectionPaperEtaTheo2->Draw();
// // 
// //         
// //     padInvSectionNLORatio->cd();
// //     padInvSectionNLORatio->SetLogx(1);
// //         TH2F * ratio2DNLOEta                = new TH2F("ratio2DNLOEta","ratio2DNLOEta",1000,0.33,25.,1000,0.2,3.45);
// //         SetStyleHistoTH2ForGraphs(ratio2DNLOEta, "#it{p}_{T} (GeV/#it{c})","#frac{NLO, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle, 
// //                                   0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.2/(textsizeFacXSecMiddle*marginXSec), 510, 505);
// //         ratio2DNLOEta->GetYaxis()->SetMoreLogLabels(kTRUE);
// //         ratio2DNLOEta->GetYaxis()->SetNdivisions(505);
// //         ratio2DNLOEta->GetYaxis()->SetNoExponent(kTRUE);
// //         ratio2DNLOEta->GetXaxis()->SetMoreLogLabels(kTRUE);
// //         ratio2DNLOEta->GetXaxis()->SetNoExponent(kTRUE);
// //         ratio2DNLOEta->GetXaxis()->SetLabelFont(42);
// //         ratio2DNLOEta->GetYaxis()->SetLabelFont(42);
// //         ratio2DNLOEta->GetYaxis()->SetLabelOffset(+0.01);
// //         ratio2DNLOEta->GetXaxis()->SetTickLength(0.07);
// //         ratio2DNLOEta->DrawCopy();
// // 
// //         DrawGammaNLOTGraph( graphRatioEtaCombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorNLO);
// //         graphRatioEtaCombNLOMuHalf->Draw("same,c");
// //         DrawGammaNLOTGraph( graphRatioEtaCombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorNLO);
// //         graphRatioEtaCombNLOMuOne->Draw("same,c");
// //         DrawGammaNLOTGraph( graphRatioEtaCombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorNLO);
// //         graphRatioEtaCombNLOMuTwo->Draw("same,c");
// //         
// //         DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
// //         graphRatioEtaCombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
// //         DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
// //         graphRatioEtaCombCombFitSys->SetLineWidth(0);
// //         graphRatioEtaCombCombFitSys->Draw("2,same");
// //         graphRatioEtaCombCombFitStatWOXErr->Draw("p,same");
// //         
// //         boxErrorSigma2760GeVRatioEta->Draw();
// //         DrawGammaLines(0.33, 25.,1., 1.,0.1,kGray);
// //         
// //     padInvSectionPythiaRatio->cd();
// //     padInvSectionPythiaRatio->SetLogx(1);
// //         TH2F * ratio2DPythiaEta             = new TH2F("ratio2DPythiaEta","ratio2DPythiaEta",1000,0.33,25.,1000,0.5,1.9);
// //         SetStyleHistoTH2ForGraphs(ratio2DPythiaEta, "#it{p}_{T} (GeV/#it{c})","#frac{PYTHIA, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown, 
// //                                   0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.2/(textsizeFacXSecDown*marginXSec), 510, 512);
// //         ratio2DPythiaEta->GetYaxis()->SetNoExponent(kTRUE);
// //         ratio2DPythiaEta->GetXaxis()->SetMoreLogLabels(kTRUE);
// //         ratio2DPythiaEta->GetXaxis()->SetNoExponent(kTRUE);
// //         ratio2DPythiaEta->GetXaxis()->SetLabelFont(42);
// //         ratio2DPythiaEta->GetYaxis()->SetLabelFont(42);
// //         ratio2DPythiaEta->GetYaxis()->SetLabelOffset(+0.01);
// //         ratio2DPythiaEta->GetXaxis()->SetTickLength(0.06);
// //         ratio2DPythiaEta->GetYaxis()->SetTickLength(0.04);
// //         ratio2DPythiaEta->DrawCopy();
// // 
// //         DrawGammaSetMarker(histoRatioEtaPythia8ToFit, 24, 1.5, kRed+2 , kRed+2);  
// //         histoRatioEtaPythia8ToFit->SetLineWidth(widthCommonFit);
// //         histoRatioEtaPythia8ToFit->Draw("same,hist,l");
// //         
// //         DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
// //         graphRatioEtaCombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
// //         DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
// //         graphRatioEtaCombCombFitSys->SetLineWidth(0);
// //         graphRatioEtaCombCombFitSys->Draw("2,same");
// //         graphRatioEtaCombCombFitStatWOXErr->Draw("p,same");
// // 
// //         boxErrorSigma2760GeVRatioEta->Draw();
// //         DrawGammaLines(0.33, 25.,1., 1.,0.1,kGray);
// //         
// //     canvasInvSectionPaper->Print(Form("%s/Eta_InvYieldWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));
// // 
// // 
    // ***************************************************************************************************************
    // ******************************* Plotting eta/pi0 ratio for single measurements ********************************
    // ***************************************************************************************************************
    textSizeLabelsPixel                 = 54;
    TCanvas* canvasEtatoPi0combo       = new TCanvas("canvasEtatoPi0combo","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasEtatoPi0combo, 0.1, 0.01, 0.01, 0.125);
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
    histo2DEtatoPi0combo               = new TH2F("histo2DEtatoPi0combo","histo2DEtatoPi0combo",1000,0.33,25.,1000,0.,1.05    );
    SetStyleHistoTH2ForGraphs(histo2DEtatoPi0combo, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.85*textsizeLabelsEtaToPi0, textsizeLabelsEtaToPi0, 
                              0.85*textsizeLabelsEtaToPi0,1.1*textsizeLabelsEtaToPi0, 0.9, 0.65, 510, 510);
    histo2DEtatoPi0combo->GetXaxis()->SetMoreLogLabels();
    histo2DEtatoPi0combo->GetXaxis()->SetLabelOffset(-0.01);
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

        TLatex *labelEnergyEtaToPi0 = new TLatex(0.13, 0.92,collisionSystempPb.Data());
        SetStyleTLatex( labelEnergyEtaToPi0, 0.85*textsizeLabelsEtaToPi0,4, 1, 42, kTRUE, 11);
        labelEnergyEtaToPi0->Draw();
        
        TLatex *labelALICEEtaToPi0 = new TLatex(0.13, 0.92-(1*textsizeLabelsEtaToPi0*0.85),"ALICE");
        SetStyleTLatex( labelALICEEtaToPi0, 0.85*textsizeLabelsEtaToPi0,4, 1, 42, kTRUE, 11);
        labelALICEEtaToPi0->Draw();
        
        TLatex *labelPi0EtaToPi0 = new TLatex(0.13, 0.92-(2*textsizeLabelsEtaToPi0*0.9),"#eta/#pi^{0}");
        SetStyleTLatex( labelPi0EtaToPi0, textsizeLabelsEtaToPi0,4, 1, 42, kTRUE, 11);
        labelPi0EtaToPi0->Draw();

    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_differentSystems.%s",outputDir.Data(), suffix.Data()));

    // ***************************************************************************************************************
    // ******************************* Plotting eta/pi0 ratio for combined measurement *******************************
    // ***************************************************************************************************************    
    histo2DEtatoPi0combo->Draw("copy");

        TGraphAsymmErrors* graphCombEtaToPi0StatWOXErr = (TGraphAsymmErrors*)graphCombEtaToPi0Stat->Clone("graphCombEtaToPi0StatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphCombEtaToPi0StatWOXErr);
    
        // plotting data
        graphCombEtaToPi0Stat->Print();
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0StatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphCombEtaToPi0StatWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphCombEtaToPi0Sys->SetLineWidth(0);
        graphCombEtaToPi0Sys->Draw("2,same");
        graphCombEtaToPi0StatWOXErr->Draw("p,same");
        
        // plotting labels
        labelEnergyEtaToPi0->Draw();
        labelALICEEtaToPi0->Draw();
        labelPi0EtaToPi0->Draw();
        
 
    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Paper.%s",outputDir.Data(), suffix.Data()));

// //     
// // //     // ***************************************************************************************************************
// // //     // ********************** Plotting eta/pi0 ratio for combined measurement + Theory *******************************
// // //     // ***************************************************************************************************************    
// // //     histo2DEtatoPi0combo->Draw("copy");
// // //         DrawGammaSetMarkerTGraphAsym(graphNLOEtaToPi0, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
// // //         graphNLOEtaToPi0->Draw("3,same");
// // // 
// // //         DrawGammaSetMarker(histoPythia8EtaToPi0, 24, 1.5, kRed+2 , kRed+2);  
// // //         histoPythia8EtaToPi0->SetLineWidth(widthCommonFit);
// // //         histoPythia8EtaToPi0->Draw("same,hist,l");
// // // 
// // //         graphCombEtaToPi0Sys->Draw("2,same");
// // //         graphCombEtaToPi0StatWOXErr->Draw("p,same");
// // //         
// // //         // plotting labels
// // //         labelEnergyEtaToPi0->Draw();
// // //         labelALICEEtaToPi0->Draw();
// // //         labelPi0EtaToPi0->Draw();
// // // 
// // //         TLegend* legendEtaToPi0Theory = GetAndSetLegend2(0.48, 0.95, 0.9, 0.95-(textsizeLabelsEtaToPi0*4*0.9), textSizeLabelsPixel*0.85, 1, "", 43, 0.16);
// // //         legendEtaToPi0Theory->AddEntry(graphCombEtaToPi0Sys,"Data","pf");
// // //         legendEtaToPi0Theory->AddEntry(graphNLOEtaToPi0,"NLO, PDF:CTEQ6M5 ","pf");
// // //         legendEtaToPi0Theory->AddEntry((TObject*)0,"#pi^{0} FF: DSS07, #eta FF: AESSS","");
// // //         legendEtaToPi0Theory->AddEntry(histoPythia8EtaToPi0,"PYTHIA 8.2, Monash 2013","l");
// // //         legendEtaToPi0Theory->Draw();
// // //         
// // //     histo2DEtatoPi0combo->Draw("axis,same");
// // // 
// // //     canvasEtatoPi0combo->Update();
// // //     canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Theory_Paper.%s",outputDir.Data(), suffix.Data()));
// // // 
    // ***************************************************************************************************************
    // **************** Plotting eta/pi0 ratio for combined measurement + Theory + other energies ********************
    // ***************************************************************************************************************    
    canvasEtatoPi0combo->cd();
    histo2DEtatoPi0combo->Draw("copy");
//         graphNLOEtaToPi0->Draw("3,same");
//         histoPythia8EtaToPi0->Draw("same,hist,l");

        graphCombEtaToPi0Sys->Draw("2,same");
               
        DrawGammaSetMarkerTGraphErr(graphPHENIXEtaToPi0200GeV, 25, 2., kGray+2, kGray+2, widthLinesBoxes, kFALSE);
        graphPHENIXEtaToPi0200GeV->Draw("p,same");
        DrawGammaSetMarkerTGraphAsym(graphALICEEtaToPi07TeV, 24, 2.2, kBlue-6, kBlue-6, widthLinesBoxes, kFALSE);
        graphALICEEtaToPi07TeV->Draw("p,same");
        
        TGraphAsymmErrors* graphALICEEtaToPi02760GeVWOXErr = (TGraphAsymmErrors*)graphALICEEtaToPi02760GeV->Clone("graphALICEEtaToPi02760GeVWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphALICEEtaToPi02760GeVWOXErr);
    
        DrawGammaSetMarkerTGraphAsym(graphALICEEtaToPi02760GeVWOXErr, 25, 2.2, kBlue-6, kBlue-6, widthLinesBoxes, kFALSE);
        graphALICEEtaToPi02760GeVWOXErr->Draw("p,same");
        
        graphCombEtaToPi0StatWOXErr->Draw("p,same");
        
//         TLegend* legendEtaToPi0Theory2 = GetAndSetLegend2(0.13, 0.95, 0.46, 0.95-(textsizeLabelsEtaToPi0*4*0.9), textSizeLabelsPixel*0.85, 1, "", 43, 0.16);
        TLegend* legendEtaToPi0Theory2 = GetAndSetLegend2(0.13, 0.95, 0.46, 0.95-(textsizeLabelsEtaToPi0*1*0.9), textSizeLabelsPixel*0.85, 1, "", 43, 0.16);
        legendEtaToPi0Theory2->AddEntry(graphCombEtaToPi0Sys,Form("ALICE, %s",collisionSystempPb.Data()),"pf");
//         legendEtaToPi0Theory2->AddEntry(graphNLOEtaToPi0,"NLO, PDF:CTEQ6M5 ","pf");
//         legendEtaToPi0Theory2->AddEntry((TObject*)0,"#pi^{0} FF: DSS07, #eta FF: AESSS","");
//         legendEtaToPi0Theory2->AddEntry(histoPythia8EtaToPi0,"PYTHIA 8.2, Monash 2013","l");
        legendEtaToPi0Theory2->Draw();

        TLegend* legendEtaToPi0WorldData = GetAndSetLegend2(0.53, 0.145+(textsizeLabelsEtaToPi0*3*0.9), 0.9, 0.145, textSizeLabelsPixel*0.85, 1, "", 43, 0.16);
        legendEtaToPi0WorldData->AddEntry(graphALICEEtaToPi07TeV,"ALICE, pp, #sqrt{#it{s}} = 7 TeV","p");
        legendEtaToPi0WorldData->AddEntry(graphALICEEtaToPi02760GeVWOXErr,"ALICE, pp, #sqrt{#it{s}} = 2.76 TeV","p");
        legendEtaToPi0WorldData->AddEntry(graphPHENIXEtaToPi0200GeV,"PHENIX, pp, #sqrt{#it{s}} = 0.2 TeV", "p");
        legendEtaToPi0WorldData->Draw();
        
    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Theory_WorldData_Paper.%s",outputDir.Data(), suffix.Data()));
// // // 
// // // 
// // //     // ***************************************************************************************************************
// // //     // **************** Plotting eta/pi0 ratio for combined measurement + other energies ********************
// // //     // ***************************************************************************************************************        
// // //     canvasEtatoPi0combo->cd();
// // //     histo2DEtatoPi0combo->Draw("copy");
// // //         
// // //         graphCombEtaToPi0Sys->Draw("2,same");
// // //                
// // //         DrawGammaSetMarkerTGraphErr(graphPHENIXEtaToPi0200GeV, 25, 2., kGray+2, kGray+2, widthLinesBoxes, kFALSE);
// // //         graphPHENIXEtaToPi0200GeV->Draw("p,same");
// // //         DrawGammaSetMarkerTGraphAsym(graphALICEEtaToPi07TeV, 24, 2.2, kBlack, kBlack, widthLinesBoxes, kFALSE);
// // //         graphALICEEtaToPi07TeV->Draw("p,same");
// // //         
// // //         graphCombEtaToPi0StatWOXErr->Draw("p,same");
// // //         
// // //         TLegend* legendEtaToPi0WorldData2 = GetAndSetLegend2(0.13, 0.95, 0.46, 0.95-(textsizeLabelsEtaToPi0*3*0.9), textSizeLabelsPixel*0.85, 1, "", 43, 0.16);
// // //         legendEtaToPi0WorldData2->AddEntry(graphCombEtaToPi0Sys,Form("ALICE, %s",collisionSystempPb.Data()),"pf");
// // //         legendEtaToPi0WorldData2->AddEntry(graphALICEEtaToPi07TeV,"ALICE, pp, #sqrt{#it{s}} = 7 TeV","p");
// // //         legendEtaToPi0WorldData2->AddEntry(graphPHENIXEtaToPi0200GeV,"PHENIX, pp, #sqrt{#it{s}} = 0.2 TeV", "p");
// // //         legendEtaToPi0WorldData2->Draw();
// // //         
// // //     histo2DEtatoPi0combo->Draw("axis,same");
// // // 
// // //     canvasEtatoPi0combo->Update();
// // //     canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_WorldData_Paper.%s",outputDir.Data(), suffix.Data()));
// //     
// //     
    // ***************************************************************************************************************
    // ******************************** fitting eta/pi0 **************************************************************
    // ***************************************************************************************************************    
    TF1* etaToPi0ConstData  = new TF1("etaToPi0ConstData","[0]",6,20);
//     TF1* etaToPi0ConstMC    = new TF1("etaToPi0ConstMC","[0]",8,20);
    graphCombEtaToPi0StatWOXErr->Fit(etaToPi0ConstData,"QRME0","",6,20);
//     histoPythia8EtaToPi0->Fit(etaToPi0ConstMC,"QRME0","",8,20);
    
    cout << "***********************************************************************************************************" << endl;
    cout << "***********************************************************************************************************" << endl;
    cout << "***********************************************************************************************************" << endl;
    cout << "high pt eta/pi0 - data, stat: " << etaToPi0ConstData->GetParameter(0) << "+-"<< etaToPi0ConstData->GetParError(0) << endl;
    graphCombEtaToPi0Tot->Fit(etaToPi0ConstData,"QRME0","",6,20);
    cout << "high pt eta/pi0 - data, tot: " << etaToPi0ConstData->GetParameter(0) << "+-"<< etaToPi0ConstData->GetParError(0) << endl;
//     cout << "high pt eta/pi0 - pythia 8: " << etaToPi0ConstMC->GetParameter(0) << "+-"<< etaToPi0ConstMC->GetParError(0) << endl;
    cout << "***********************************************************************************************************" << endl;
    cout << "***********************************************************************************************************" << endl;
    cout << "***********************************************************************************************************" << endl;


    canvasXSectionPi0->cd();
    TH2F * histo2DXSectionWithEtaAndPi0;
    histo2DXSectionWithEtaAndPi0          = new TH2F("histo2DXSectionWithEtaAndPi0","histo2DXSectionWithEtaAndPi0",11000,0.23,31.,1000,2e-12,10e1);
    SetStyleHistoTH2ForGraphs(histo2DXSectionWithEtaAndPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.65);
    histo2DXSectionWithEtaAndPi0->GetXaxis()->SetMoreLogLabels();
    histo2DXSectionWithEtaAndPi0->GetXaxis()->SetLabelOffset(-0.01);
    histo2DXSectionWithEtaAndPi0->Draw("copy");

        // scale eta graphs
        Double_t scaleFacEtaForCombPlot                             = 1e-2;
        TGraphAsymmErrors* graphCombEtaInvYieldStatWOXErrCopy   = (TGraphAsymmErrors*) graphCombEtaInvYieldStatWOXErr->Clone("graphCombEtaInvYieldStatWOXErrCopy");
        TGraphAsymmErrors* graphCombEtaInvYieldSysCopy          = (TGraphAsymmErrors*) graphCombEtaInvYieldSys->Clone("graphCombEtaInvYieldSysCopy");
        graphCombEtaInvYieldStatWOXErrCopy                      = ScaleGraph(graphCombEtaInvYieldStatWOXErrCopy,scaleFacEtaForCombPlot);
        graphCombEtaInvYieldSysCopy                             = ScaleGraph(graphCombEtaInvYieldSysCopy,scaleFacEtaForCombPlot);
        TH1D* histFitTCMInvYieldEta                              = (TH1D*)fitTCMInvYieldEta->GetHistogram();
        histFitTCMInvYieldEta->Scale(scaleFacEtaForCombPlot);
//         histoPythia8InvYieldEta->Scale(scaleFacEtaForCombPlot);
//         TGraphAsymmErrors* graphNLOEtaAESSCalcCopy                  = (TGraphAsymmErrors*)graphNLOEtaAESSCalc->Clone("graphNLOEtaAESSCalcCopy");
//         TGraph* graphNLOCalcEtaMuHalfCopy                           = (TGraph*)graphNLOCalcEtaMuHalf->Clone("graphNLOCalcEtaMuHalfCopy");
//         TGraph* graphNLOCalcEtaMuOneCopy                            = (TGraph*)graphNLOCalcEtaMuOne->Clone("graphNLOCalcEtaMuOneCopy");
//         TGraph* graphNLOCalcEtaMuTwoCopy                            = (TGraph*)graphNLOCalcEtaMuTwo->Clone("graphNLOCalcEtaMuTwoCopy");
//         graphNLOEtaAESSCalcCopy                                     = ScaleGraph(graphNLOEtaAESSCalcCopy,scaleFacEtaForCombPlot);
//         graphNLOCalcEtaMuHalfCopy                                   = ScaleGraph(graphNLOCalcEtaMuHalfCopy,scaleFacEtaForCombPlot);
//         graphNLOCalcEtaMuOneCopy                                    = ScaleGraph(graphNLOCalcEtaMuOneCopy,scaleFacEtaForCombPlot);
//         graphNLOCalcEtaMuTwoCopy                                    = ScaleGraph(graphNLOCalcEtaMuTwoCopy,scaleFacEtaForCombPlot);
        
        // plotting Pythia 8.2 Monash
//         DrawGammaSetMarker(histoPythia8InvYieldPi0, 24, 1.5, kRed+2 , kRed+2);  
//         histoPythia8InvYieldPi0->SetLineWidth(widthCommonFit);
//         histoPythia8InvYieldPi0->Draw("same,hist,l");
//         DrawGammaSetMarker(histoPythia8InvYieldEta, 24, 1.5, kRed+2 , kRed+2);  
//         histoPythia8InvYieldEta->SetLineWidth(widthCommonFit);
//         histoPythia8InvYieldEta->Draw("same,hist,l");

//         // plotting NLO calcs pi0
//         DrawGammaSetMarkerTGraphAsym(graphNLODSS14Calc, 0, 0, colorNLO, colorNLO, widthLinesBoxes, kTRUE, colorNLO);
//         graphNLODSS14Calc->Draw("3,same");
        
//         // plotting NLO calcs eta
//         DrawGammaSetMarkerTGraphAsym(graphNLOEtaAESSCalcCopy, 0, 0, colorCGC, colorCGC, widthLinesBoxes, kTRUE, colorCGC);
//         graphNLOEtaAESSCalcCopy->Draw("3,same");       

        // plot data
        graphCombPi0InvYieldSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSysCopy, markerStyleComb+4, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
        graphCombEtaInvYieldSysCopy->Draw("E2same");
        
        graphCombPi0InvYieldStatWOXErr->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStatWOXErrCopy, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
        graphCombEtaInvYieldStatWOXErrCopy->Draw("p,same,z");
        
        // plots fits
        fitTCMInvYieldPi0->Draw("same");                
        SetStyleHisto(histFitTCMInvYieldEta, 2, 7, kGray+2);
        histFitTCMInvYieldEta->Draw("same,c");
                        
        // labels lower left corner
        TLegend* legendXsectionPaperAll    = GetAndSetLegend2(0.19, 0.19, 0.5, 0.19+0.04*1, textSizeLabelsPixel, 2, "", 43, 0.4);
        legendXsectionPaperAll->AddEntry(graphCombPi0InvYieldSys,"#pi^{0}","pf");
        legendXsectionPaperAll->AddEntry(graphCombEtaInvYieldSysCopy,"#eta (x 10^{-2})","pf");
        legendXsectionPaperAll->Draw();
        TLegend* legendXsectionPaper3     = GetAndSetLegend2(0.19, 0.05+0.08, 0.5, 0.11+0.08, textSizeLabelsPixel, 1, "", 43, 0.2);
        legendXsectionPaper3->AddEntry(fitTCMInvYieldPi0,"#it{A}_{e} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{n}","l");
        legendXsectionPaper3->Draw();
        
        TLatex *labelEnergyXSectionPaperAll = new TLatex(0.20, 0.20+0.04*3, collisionSystempPb.Data());
        SetStyleTLatex( labelEnergyXSectionPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelEnergyXSectionPaperAll->Draw();
        TLatex *labelALICEXSectionPaperAll  = new TLatex(0.20,0.20+0.04*2,"ALICE");
        SetStyleTLatex( labelALICEXSectionPaperAll, textSizeLabelsPixel,4, 1, 43, kTRUE, 11);
        labelALICEXSectionPaperAll->Draw();
        TLatex *labelALICENormUnPaperAll    = new TLatex(0.20,0.20+0.05*1,"Norm. unc. 5.4%");
        SetStyleTLatex( labelALICENormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelALICENormUnPaperAll->Draw();
        

//         // labels upper right corner
//         TLegend* legendXsectionPaperPyBoth  = GetAndSetLegend2(0.54, 0.95-0.04*4-0.04*0.75*2, 0.59+0.33, 0.95, textSizeLabelsPixel, 1, "", 43, 0.18);
//         legendXsectionPaperPyBoth->AddEntry(histoPythia8InvYieldEta,"PYTHIA 8.2","l");
//         legendXsectionPaperPyBoth->AddEntry((TObject*)0,"Monash 2013","");
//         legendXsectionPaperPyBoth->AddEntry(graphNLODSS14Calc,"#pi^{0} pQCD NLO ","f");
//         legendXsectionPaperPyBoth->AddEntry((TObject*)0,"#scale[0.75]{PDF: MSTW, FF: DSS14}","");
//         legendXsectionPaperPyBoth->AddEntry(graphNLOEtaAESSCalcCopy,"#eta pQCD NLO ","f");        
//         legendXsectionPaperPyBoth->AddEntry((TObject*)0,"#scale[0.75]{PDF: CTEQ6M5, FF: AESSS}","");
//         legendXsectionPaperPyBoth->Draw();

    canvasXSectionPi0->SaveAs(Form("%s/InvYield_Pi0_Eta_Theory.%s",outputDir.Data(),suffix.Data()));

    
    TF1* fitPowInvYieldPi0Tot   = FitObject("powPure","fitPowInvYieldPi02760GeVTot","Pi0",graphCombPi0InvYieldTot,6,40. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvYieldPi0Tot)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldPi0Tot)<< endl;    
    TF1* fitPowInvYieldPi0Stat   = FitObject("powPure","fitPowInvYieldPi02760GeV","Pi0",graphCombPi0InvYieldStat,6,40. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvYieldPi0Stat)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldPi0Stat)<< endl;    
    TF1* fitPowInvYieldEtaTot   = FitObject("powPure","fitPowInvYieldEta2760GeVTot","Eta",graphCombEtaInvYieldTot,6,20. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvYieldEtaTot)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldEtaTot)<< endl;    
    TF1* fitPowInvYieldEtaStat   = FitObject("powPure","fitPowInvYieldEta2760GeVStat","Eta",graphCombEtaInvYieldStat,6,20. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvYieldEtaStat)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldEtaStat)<< endl;    
    
    canvasXSectionPi0->cd();
    histo2DYieldPi0->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys, markerStyleComb, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
        graphCombPi0InvYieldSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat, markerStyleComb, markerSizeComb, kBlack , kBlack);
        graphCombPi0InvYieldStat->Draw("p,same,z");
        
        DrawGammaSetMarkerTF1( fitInvYieldPi0, 9, 2, kGray+2); 
        DrawGammaSetMarkerTF1( fitTCMInvYieldPi0, 7, 2, kGray+1); 
        DrawGammaSetMarkerTF1( fitPowInvYieldPi0Tot, 1, 2, kAzure+2); 
        DrawGammaSetMarkerTF1( fitTCMDecomposedLPi0, 4, 2, kRed-6); 
        DrawGammaSetMarkerTF1( fitTCMDecomposedHPi0, 5, 2, kRed-8); 
        fitPowInvYieldPi0Tot->SetRange(6,50);
        fitTCMInvYieldPi0->SetRange(0.2,50);
        fitInvYieldPi0->SetRange(0.2,50);
        fitTCMDecomposedLPi0->SetRange(0.2,10);
        fitTCMDecomposedHPi0->SetRange(0.2,50);
        fitPowInvYieldPi0Tot->Draw("same");
        fitInvYieldPi0->Draw("same");
        fitTCMInvYieldPi0->Draw("same");
        fitTCMDecomposedLPi0->Draw("same");
        fitTCMDecomposedHPi0->Draw("same");
        
        labelEnergyXSectionPi0->Draw();
        labelDetSysXSectionPi0->Draw();
        
        TLegend* legendXSectionPi0DiffFits          = GetAndSetLegend2(0.2, 0.13, 0.40 , 0.13+0.03*5, 32); 
        legendXSectionPi0DiffFits->AddEntry(fitInvYieldPi0,"Levy-Tsallis","l");
        legendXSectionPi0DiffFits->AddEntry(fitPowInvYieldPi0Tot,"pure powerlaw","l");
        legendXSectionPi0DiffFits->AddEntry(fitTCMInvYieldPi0,"full TCM","l");
        legendXSectionPi0DiffFits->AddEntry(fitTCMDecomposedLPi0,"low p_{T} comp. TCM","l");
        legendXSectionPi0DiffFits->AddEntry(fitTCMDecomposedHPi0,"high p_{T} comp. TCM","l");
        legendXSectionPi0DiffFits->Draw();
        
    canvasXSectionPi0->SaveAs(Form("%s/InvYield_Pi0_WithDiffFits.%s",outputDir.Data(),suffix.Data()));


    histo2DInvYieldEta->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys, markerStyleComb, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
        graphCombEtaInvYieldSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat, markerStyleComb, markerSizeComb, kBlack , kBlack);
        graphCombEtaInvYieldStat->Draw("p,same,z");

        DrawGammaSetMarkerTF1( fitInvYieldEta, 9, 2, kGray+2); 
        DrawGammaSetMarkerTF1( fitTCMInvYieldEta, 7, 2, kGray+1); 
        DrawGammaSetMarkerTF1( fitPowInvYieldEtaTot, 1, 2, kAzure+2); 
        DrawGammaSetMarkerTF1( fitTCMDecomposedLEta, 4, 2, kRed-6); 
        DrawGammaSetMarkerTF1( fitTCMDecomposedHEta, 5, 2, kRed-8); 
        fitPowInvYieldEtaTot->SetRange(6,50);
        fitTCMInvYieldEta->SetRange(0.2,50);
        fitInvYieldEta->SetRange(0.2,50);
        fitTCMDecomposedLEta->SetRange(0.2,10);
        fitTCMDecomposedHEta->SetRange(0.2,50);
        fitPowInvYieldEtaTot->Draw("same");
        fitInvYieldEta->Draw("same");
        fitTCMInvYieldEta->Draw("same");
        fitTCMDecomposedLEta->Draw("same");
        fitTCMDecomposedHEta->Draw("same");
        
        labelEnergyXSectionPi0->Draw();
        labelDetSysXSectionEta->Draw();
        
        TLegend* legendXSectionEtaDiffFits          = GetAndSetLegend2(0.2, 0.13, 0.40 , 0.13+0.03*5, 32); 
        legendXSectionEtaDiffFits->AddEntry(fitInvYieldEta,"Levy-Tsallis","l");
        legendXSectionEtaDiffFits->AddEntry(fitPowInvYieldEtaTot,"pure powerlaw","l");
        legendXSectionEtaDiffFits->AddEntry(fitTCMInvYieldEta,"full TCM","l");
        legendXSectionEtaDiffFits->AddEntry(fitTCMDecomposedLEta,"low p_{T} comp. TCM","l");
        legendXSectionEtaDiffFits->AddEntry(fitTCMDecomposedHEta,"high p_{T} comp. TCM","l");
        legendXSectionEtaDiffFits->Draw();
        
    canvasXSectionPi0->SaveAs(Form("%s/InvYield_Eta_WithDiffFits.%s",outputDir.Data(),suffix.Data()));
    
    
    //*************************************************************************************************************
    //***************************** Comparison to Charged pions ***************************************************
    //*************************************************************************************************************    
    cout << "combined Spectrum - published" << endl;
    TGraphErrors* graphChPiInvYieldPubStatPubCombUp         = NULL;
    TGraphErrors* graphChPiInvYieldPubSystPubCombUp         = NULL;
    TGraphErrors* graphCombPi0InvYieldStatRebinnedPubComb   = NULL;
    TGraphErrors* graphCombPi0InvYieldSysRebinnedPubComb    = NULL;
    TGraphErrors* graphRatioPi0PubChPiComb                  = CalculateRatioBetweenSpectraWithDifferentBinning( graphCombPi0InvYieldStat, graphCombPi0InvYieldSys, 
                                                                                                                histoChPiInvYieldPubStat, histoChPiInvYieldPubSyst,  
                                                                                                                kTRUE,  kTRUE, 
                                                                                                                &graphCombPi0InvYieldStatRebinnedPubComb, &graphCombPi0InvYieldSysRebinnedPubComb, 
                                                                                                                &graphChPiInvYieldPubStatPubCombUp, &graphChPiInvYieldPubSystPubCombUp )    ;
    cout << "PHOS Spectrum - published" << endl;
    TGraphErrors* graphChPiInvYieldPubStatPubPHOS           = NULL;
    TGraphErrors* graphChPiInvYieldPubSystPubPHOS           = NULL;
    TGraphErrors* graphPHOSPi0InvYieldStatRebinnedPubPHOS   = NULL;
    TGraphErrors* graphPHOSPi0InvYieldSysRebinnedPubPHOS    = NULL;
    TGraphErrors* graphRatioPi0PubChPiPHOS                  = CalculateRatioBetweenSpectraWithDifferentBinning( graphPHOSPi0InvYieldStat, graphPHOSPi0InvYieldSys, 
                                                                                                                histoChPiInvYieldPubStat, histoChPiInvYieldPubSyst,  
                                                                                                                kTRUE,  kTRUE, 
                                                                                                                &graphPHOSPi0InvYieldStatRebinnedPubPHOS, &graphPHOSPi0InvYieldSysRebinnedPubPHOS, 
                                                                                                                &graphChPiInvYieldPubStatPubPHOS, &graphChPiInvYieldPubSystPubPHOS )    ;
    cout << "PCM Spectrum - published" << endl;
    TGraphErrors* graphChPiInvYieldPubStatPubPCM            = NULL;
    TGraphErrors* graphChPiInvYieldPubSystPubPCM            = NULL;
    TGraphErrors* graphPCMPi0InvYieldStatRebinnedPubPCM     = NULL;
    TGraphErrors* graphPCMPi0InvYieldSysRebinnedPubPCM      = NULL;
    TGraphErrors* graphRatioPi0PubChPiPCM                   = CalculateRatioBetweenSpectraWithDifferentBinning(  graphPCMPi0InvYieldStat, graphPCMPi0InvYieldSys, 
                                                                                                                histoChPiInvYieldPubStat, histoChPiInvYieldPubSyst,  
                                                                                                                kTRUE,  kTRUE, 
                                                                                                                &graphPCMPi0InvYieldStatRebinnedPubPCM, &graphPCMPi0InvYieldSysRebinnedPubPCM, 
                                                                                                                &graphChPiInvYieldPubStatPubPCM, &graphChPiInvYieldPubSystPubPCM )    ;
    cout << "PCM-EMC Spectrum - published" << endl;
    TGraphErrors* graphChPiInvYieldPubStatPubPCMEMC             = NULL;
    TGraphErrors* graphChPiInvYieldPubSystPubPCMEMC             = NULL;
    TGraphErrors* graphPCMEMCPi0InvYieldStatRebinnedPubPCMEMC   = NULL;
    TGraphErrors* graphPCMEMCPi0InvYieldSysRebinnedPubPCMEMC    = NULL;
    TGraphErrors* graphRatioPi0PubChPiPCMEMC                    = CalculateRatioBetweenSpectraWithDifferentBinning( graphPCMEMCALPi0InvYieldStat, graphPCMEMCALPi0InvYieldSys, 
                                                                                                                histoChPiInvYieldPubStat, histoChPiInvYieldPubSyst,  
                                                                                                                kTRUE,  kTRUE, 
                                                                                                                &graphPCMEMCPi0InvYieldStatRebinnedPubPCMEMC, &graphPCMEMCPi0InvYieldSysRebinnedPubPCMEMC, 
                                                                                                                &graphChPiInvYieldPubStatPubPCMEMC, &graphChPiInvYieldPubSystPubPCMEMC )    ;
    cout << "EMC Spectrum - published" << endl;
    TGraphErrors* graphChPiInvYieldPubStatPubEMC            = NULL;
    TGraphErrors* graphChPiInvYieldPubSystPubEMC            = NULL;
    TGraphErrors* graphEMCPi0InvYieldStatRebinnedPubEMC     = NULL;
    TGraphErrors* graphEMCPi0InvYieldSysRebinnedPubEMC      = NULL;
    TGraphErrors* graphRatioPi0PubChPiEMC                   = CalculateRatioBetweenSpectraWithDifferentBinning( graphEMCALPi0InvYieldStat, graphEMCALPi0InvYieldSys, 
                                                                                                                histoChPiInvYieldPubStat, histoChPiInvYieldPubSyst,  
                                                                                                                kTRUE,  kTRUE, 
                                                                                                                &graphEMCPi0InvYieldStatRebinnedPubEMC, &graphEMCPi0InvYieldSysRebinnedPubEMC, 
                                                                                                                &graphChPiInvYieldPubStatPubEMC, &graphChPiInvYieldPubSystPubEMC )    ;

    cout << "Dalitz Spectrum - published" << endl;
    TGraphErrors* graphChPiInvYieldPubStatPubDal            = NULL;
    TGraphErrors* graphChPiInvYieldPubSystPubDal            = NULL;
    TGraphErrors* graphDalPi0InvYieldStatRebinnedPubDal     = NULL;
    TGraphErrors* graphDalPi0InvYieldSysRebinnedPubDal      = NULL;
    TGraphErrors* graphRatioPi0PubChPiDal                   = CalculateRatioBetweenSpectraWithDifferentBinning( graphDalitzPi0InvYieldStat, graphDalitzPi0InvYieldSys, 
                                                                                                                histoChPiInvYieldPubStat, histoChPiInvYieldPubSyst,  
                                                                                                                kTRUE,  kTRUE, 
                                                                                                                &graphDalPi0InvYieldStatRebinnedPubDal, &graphDalPi0InvYieldSysRebinnedPubDal, 
                                                                                                                &graphChPiInvYieldPubStatPubDal, &graphChPiInvYieldPubSystPubDal )    ;

    // ***************************************************************************************************************
    // ************************** Comparison pi0/pi+-, pi0 updated comb **********************************************
    // ***************************************************************************************************************    
    textSizeLabelsPixel             = 48;
    TCanvas* canvasCompYield   = new TCanvas("canvasCompYield","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasCompYield,   0.12, 0.01, 0.01, 0.11);
    canvasCompYield->SetLogx();

    TH2F * histo2DCompCombinedRatio;
    histo2DCompCombinedRatio       = new TH2F("histo2DCompCombinedRatio","histo2DCompCombinedRatio",1000,0.23,30.,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.85*textsizeLabelsPP, textsizeLabelsPP, 
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
    histo2DCompCombinedRatio->GetYaxis()->SetNoExponent(kTRUE);
//  histo2DCompCombinedRatio->GetXaxis()->SetLabelOffset(-0.01);
    histo2DCompCombinedRatio->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DCompCombinedRatio->GetXaxis()->SetNoExponent(kTRUE);
    histo2DCompCombinedRatio->GetYaxis()->SetRangeUser(0.05,2.45);
    histo2DCompCombinedRatio->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioPi0PubChPiComb, markerStyleCombLowPt, markerSizeComparison, colorCombLowPt , colorCombLowPt);
        graphRatioPi0PubChPiComb->Draw("E1psame");

        TLegend* legendPi0CompChargedPionsComb    = GetAndSetLegend2(0.15, 0.8, 0.9, 0.90, 0.85* textSizeLabelsPixel, 2, "", 43, 0.12);
        legendPi0CompChargedPionsComb->AddEntry(graphRatioPi0PubChPiComb,"#pi^{0}/#pi^{#pm}","p");
        legendPi0CompChargedPionsComb->Draw();
    
        labelRatioTheory->Draw();
        
        DrawGammaLines(0., 30 , 1, 1 ,1, kGray, 1);   
   
    canvasCompYield->Update();
    canvasCompYield->Print(Form("%s/ComparisonChargedToNeutralComb_PPb5023GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));

    // ***************************************************************************************************************
    // ************************** Comparison pi0/pi+-, pi0 updated individual ****************************************
    // ***************************************************************************************************************    
    histo2DCompCombinedRatio->DrawCopy();

        DrawGammaSetMarkerTGraphErr(graphRatioPi0PubChPiPCM, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
        DrawGammaSetMarkerTGraphErr(graphRatioPi0PubChPiDal, markerStyleDet[5] ,markerSizeDet[5]*0.5, colorDet[5], colorDet[5]);
        DrawGammaSetMarkerTGraphErr(graphRatioPi0PubChPiPHOS, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);
        DrawGammaSetMarkerTGraphErr(graphRatioPi0PubChPiEMC, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);
        DrawGammaSetMarkerTGraphErr(graphRatioPi0PubChPiPCMEMC, markerStyleDet[4] ,markerSizeDet[4]*0.5, colorDet[4], colorDet[4]);        
        
        graphRatioPi0PubChPiPCM->Draw("p,same,z");
        graphRatioPi0PubChPiDal->Draw("p,same,z");
        graphRatioPi0PubChPiPHOS->Draw("p,same,z");
        graphRatioPi0PubChPiEMC->Draw("p,same,z");
        graphRatioPi0PubChPiPCMEMC->Draw("p,same,z");
   
        TLegend* legendPi0CompChargedPionsInd    = GetAndSetLegend2(0.15, 0.15, 0.9, 0.15+3*0.035, 0.85* textSizeLabelsPixel, 2, "", 43, 0.12);
        legendPi0CompChargedPionsInd->AddEntry(graphRatioPi0PubChPiPCM,"#pi^{0}/#pi^{#pm} PCM","p");
        legendPi0CompChargedPionsInd->AddEntry(graphRatioPi0PubChPiPHOS,"#pi^{0}/#pi^{#pm} PHOS","p");
        legendPi0CompChargedPionsInd->AddEntry(graphRatioPi0PubChPiEMC,"#pi^{0}/#pi^{#pm} EMC","p");
        legendPi0CompChargedPionsInd->AddEntry(graphRatioPi0PubChPiDal,"#pi^{0}/#pi^{#pm} PCM-Dal","p");
        legendPi0CompChargedPionsInd->AddEntry(graphRatioPi0PubChPiPCMEMC,"#pi^{0}/#pi^{#pm} PCM-EMC","p");
        legendPi0CompChargedPionsInd->Draw();
   
        labelRatioTheory->Draw();
        
        DrawGammaLines(0., 30 , 1, 1 ,1, kGray, 1);   
   
    canvasCompYield->Update();
    canvasCompYield->Print(Form("%s/ComparisonChargedToNeutralIndividual_PPb5023GeV_%s.%s",outputDir.Data(),dateForOutput.Data(),suffix.Data()));
    
    // close fit log file
    fileFitsOutput.close();
    
 // **********************************************************************************************************************
 // ************************* Saving of final results ********************************************************************
 // **********************************************************************************************************************
 
    TString nameOutputCommonFile    = Form("CombinedResultsPaperPPb5023GeV_%s", dateForOutput.Data());
    
    TFile fCombResults(nameOutputCommonFile.Data(), "RECREATE");

    fCombResults.mkdir("Pi0pPb_5.023TeV");
    TDirectoryFile* directoryPi0 = (TDirectoryFile*)fCombResults.Get("Pi0pPb_5.023TeV"); 
    fCombResults.cd("Pi0pPb_5.023TeV");
        // PCM component
        graphPCMPi0InvYieldStat->Write("graphInvYieldPi0PCMpPb5023GeVStatErr");
        graphPCMPi0InvYieldSys->Write("graphInvYieldPi0PCMpPb5023GeVSysErr");
        // PHOS component
        graphPHOSPi0InvYieldStat->Write("graphInvYieldPi0PHOSpPb5023GeVStatErr");
        graphPHOSPi0InvYieldSys->Write("graphInvYieldPi0PHOSpPb5023GeVSysErr");
        // EMCAL component
        graphEMCALPi0InvYieldStat->Write("graphInvYieldPi0EMCALpPb5023GeVStatErr");
        graphEMCALPi0InvYieldSys->Write("graphInvYieldPi0EMCALpPb5023GeVSysErr");
        // PCM-EMC component
        graphPCMEMCALPi0InvYieldStat->Write("graphInvYieldPi0PCMEMCALpPb5023GeVStatErr");
        graphPCMEMCALPi0InvYieldSys->Write("graphInvYieldPi0PCMEMCALpPb5023GeVSysErr");
        // Dalitz component
        graphDalitzPi0InvYieldStat->Write("graphInvYieldPi0DalitzpPb5023GeVStatErr");
        graphDalitzPi0InvYieldSys->Write("graphInvYieldPi0DalitzpPb5023GeVSysErr");
        // Final spectrum correlations Method A
        graphCombPi0InvYieldTot->Write("graphInvYieldPi0CombpPb5023GeVATotErr");
        graphCombPi0InvYieldStat->Write("graphInvYieldPi0CombpPb5023GeVAStatErr");
        graphCombPi0InvYieldSys->Write("graphInvYieldPi0CombpPb5023GeVASysErr");  
        // Final inv yield INEL
        graphCombPi0InvYieldTot->Write("graphInvYieldINELPi0CombpPb5023GeVATotErr");
        graphCombPi0InvYieldStat->Write("graphInvYieldINELPi0CombpPb5023GeVAStatErr");
        graphCombPi0InvYieldSys->Write("graphInvYieldINELPi0CombpPb5023GeVASysErr");  

         // fits for eta
        fitInvYieldPi0->Write("TsallisFitPi0");
        fitTCMInvYieldPi0->Write("TwoComponentModelFitPi0");
                      
        if (bWCorrection.Contains("Y")){
            if(graphPCMPi0InvYieldStat_yShifted)graphPCMPi0InvYieldStat_yShifted->Write("graphInvYieldPi0PCMpPb5023GeVStatErr_yShifted");
            if(graphPCMPi0InvYieldSys_yShifted)graphPCMPi0InvYieldSys_yShifted->Write("graphInvYieldPi0PCMpPb5023GeVSysErr_yShifted");
            // PHOS component
            if(graphPHOSPi0InvYieldStat_yShifted) graphPHOSPi0InvYieldStat_yShifted->Write("graphInvYieldPi0PHOSpPb5023GeVStatErr_yShifted");
            if(graphPHOSPi0InvYieldSys_yShifted) graphPHOSPi0InvYieldSys_yShifted->Write("graphInvYieldPi0PHOSpPb5023GeVSysErr_yShifted");
            // EMCAL component
            if(graphEMCALPi0InvYieldStat_yShifted)graphEMCALPi0InvYieldStat_yShifted->Write("graphInvYieldPi0EMCALpPb5023GeVStatErr_yShifted");
            if(graphEMCALPi0InvYieldSys_yShifted)graphEMCALPi0InvYieldSys_yShifted->Write("graphInvYieldPi0EMCALpPb5023GeVSysErr_yShifted");
            // PCM-EMC component
            if(graphPCMEMCALPi0InvYieldStat_yShifted)graphPCMEMCALPi0InvYieldStat_yShifted->Write("graphInvYieldPi0PCMEMCALpPb5023GeVStatErr_yShifted");
            if(graphPCMEMCALPi0InvYieldSys_yShifted)graphPCMEMCALPi0InvYieldSys_yShifted->Write("graphInvYieldPi0PCMEMCALpPb5023GeVSysErr_yShifted");
            // Dalitz component
            if(graphDalitzPi0InvYieldStat_yShifted)graphDalitzPi0InvYieldStat_yShifted->Write("graphInvYieldPi0DalitzpPb5023GeVStatErr_yShifted");
            if(graphDalitzPi0InvYieldSys_yShifted)graphDalitzPi0InvYieldSys_yShifted->Write("graphInvYieldPi0DalitzpPb5023GeVSysErr_yShifted");
            // Final spectrum correlations Method A
            if(graphCombPi0InvYieldTot_yShifted)graphCombPi0InvYieldTot_yShifted->Write("graphInvYieldPi0CombpPb5023GeVATotErr_yShifted");
            if(graphCombPi0InvYieldStat_yShifted)graphCombPi0InvYieldStat_yShifted->Write("graphInvYieldPi0CombpPb5023GeVAStatErr_yShifted");
            if(graphCombPi0InvYieldSys_yShifted)graphCombPi0InvYieldSys_yShifted->Write("graphInvYieldPi0CombpPb5023GeVASysErr_yShifted");  

        }    
        directoryPi0->mkdir("Supporting");
        directoryPi0->cd("Supporting");
            // Writing full correction factors
            graphPCMPi0AccTimesEff->Write("Pi0CorrectionFactorPCM");
//             graphPHOSPi0AccTimesEff->Write("Pi0CorrectionFactorPHOS");
            graphEMCALPi0AccTimesEff->Write("Pi0CorrectionFactorEMCAL");
            graphPCMEMCALPi0AccTimesEff->Write("Pi0CorrectionFactorPCMEMCAL");
            graphPCMEMCALPi0AccTimesEff->Write("Pi0CorrectionFactorDalitz");
            
            graphPCMPi0Mass->Write("Pi0MassDataPCM");
            graphPCMPi0MassMC->Write("Pi0MassMCPCM");
            histoPHOSPi0Mass->Write("Pi0MassDataPHOS");
            histoPHOSPi0TrueMass->Write("Pi0MassMCPHOS");
            graphEMCALPi0Mass->Write("Pi0MassDataEMCAL");
            graphEMCALPi0MassMC->Write("Pi0MassMCEMCAL");
            graphPCMEMCALPi0Mass->Write("Pi0MassDataPCMEMCAL");
            graphPCMEMCALPi0MassMC->Write("Pi0MassMCPCMEMCAL");
            graphDalitzPi0Mass->Write("Pi0MassDataDalitz");
            graphDalitzPi0MassMC->Write("Pi0MassMCDalitz");

            graphPCMPi0FWHM->Write("Pi0WidthDataPCM");
            graphPCMPi0FWHMMC->Write("Pi0WidthMCPCM");
            histoPHOSPi0FWHMMeV->Write("Pi0WidthDataPHOS");
            histoPHOSPi0TrueFWHMMeV->Write("Pi0WidthMCPHOS");
            graphEMCALPi0FWHM->Write("Pi0WidthDataEMCAL");
            graphEMCALPi0FWHMMC->Write("Pi0WidthMCEMCAL");
            graphPCMEMCALPi0FWHM->Write("Pi0WidthDataPCMEMCAL");
            graphPCMEMCALPi0FWHMMC->Write("Pi0WidthMCPCMEMCAL");
            graphDalitzPi0FWHM->Write("Pi0WidthDataDalitz");
            graphDalitzPi0FWHMMC->Write("Pi0WidthMCDalitz");
            
    fCombResults.mkdir("EtapPb_5.023TeV");
    TDirectoryFile* directoryEta = (TDirectoryFile*)fCombResults.Get("EtapPb_5.023TeV"); 
    fCombResults.cd("EtapPb_5.023TeV");
        // PCM component
        graphPCMEtaInvYieldStat->Write("graphInvYieldEtaPCMpPb5023GeVStatErr");
        graphPCMEtaInvYieldSys->Write("graphInvYieldEtaPCMpPb5023GeVSysErr");
        // EMCAL component
        graphEMCALEtaInvYieldStat->Write("graphInvYieldEtaEMCALpPb5023GeVStatErr");
        graphEMCALEtaInvYieldSys->Write("graphInvYieldEtaEMCALpPb5023GeVSysErr");
        // PCM-EMC component
        graphPCMEMCALEtaInvYieldStat->Write("graphInvYieldEtaPCMEMCALpPb5023GeVStatErr");
        graphPCMEMCALEtaInvYieldSys->Write("graphInvYieldEtaPCMEMCALpPb5023GeVSysErr");
        // Final spectrum correlations Method A
        graphCombEtaInvYieldTot->Write("graphInvYieldEtaCombpPb5023GeVATotErr");
        graphCombEtaInvYieldStat->Write("graphInvYieldEtaCombpPb5023GeVAStatErr");
        graphCombEtaInvYieldSys->Write("graphInvYieldEtaCombpPb5023GeVASysErr");  

        // Final yield for INEL
        graphCombEtaInvYieldTot->Write("graphInvYieldINELEtaCombpPb5023GeVATotErr");
        graphCombEtaInvYieldStat->Write("graphInvYieldINELEtaCombpPb5023GeVAStatErr");
        graphCombEtaInvYieldSys->Write("graphInvYieldINELEtaCombpPb5023GeVASysErr");  
        
        // fits for eta
        fitInvYieldEta->Write("TsallisFitEta");
        fitTCMInvYieldEta->Write("TwoComponentModelFitEta");
        
        // writing Y shifted graphs in addition
        if (bWCorrection.Contains("Y")){
            if(graphPCMEtaInvYieldStat_yShifted)graphPCMEtaInvYieldStat_yShifted->Write("graphInvYieldEtaPCMpPb5023GeVStatErr_yShifted");
            if(graphPCMEtaInvYieldSys_yShifted)graphPCMEtaInvYieldSys_yShifted->Write("graphInvYieldEtaPCMpPb5023GeVSysErr_yShifted");
            // EMCAL component
            if(graphEMCALEtaInvYieldStat_yShifted)graphEMCALEtaInvYieldStat_yShifted->Write("graphInvYieldEtaEMCALpPb5023GeVStatErr_yShifted");
            if(graphEMCALEtaInvYieldSys_yShifted)graphEMCALEtaInvYieldSys_yShifted->Write("graphInvYieldEtaEMCALpPb5023GeVSysErr_yShifted");
            // PCM-EMC component
            if(graphPCMEMCALEtaInvYieldStat_yShifted)graphPCMEMCALEtaInvYieldStat_yShifted->Write("graphInvYieldEtaPCMEMCALpPb5023GeVStatErr_yShifted");
            if(graphPCMEMCALEtaInvYieldSys_yShifted)graphPCMEMCALEtaInvYieldSys_yShifted->Write("graphInvYieldEtaPCMEMCALpPb5023GeVSysErr_yShifted");
            // Final spectrum correlations Method A
            if(graphCombEtaInvYieldTot_yShifted)graphCombEtaInvYieldTot_yShifted->Write("graphInvYieldEtaCombpPb5023GeVATotErr_yShifted");
            if(graphCombEtaInvYieldStat_yShifted)graphCombEtaInvYieldStat_yShifted->Write("graphInvYieldEtaCombpPb5023GeVAStatErr_yShifted");
            if(graphCombEtaInvYieldSys_yShifted)graphCombEtaInvYieldSys_yShifted->Write("graphInvYieldEtaCombpPb5023GeVASysErr_yShifted");  
        }    

        histoPCMEtaToPi0Stat->Write("histoRatioEtaToPi0PCMpPb5023GeVStatErr");
        graphPCMEtaToPi0Sys->Write("graphRatioEtaToPi0PCMpPb5023GeVSysErr");
        histoEMCALEtaToPi0Stat->Write("histoRatioEtaToPi0EMCALpPb5023GeVStatErr");
        graphEMCALEtaToPi0Sys->Write("graphRatioEtaToPi0EMCALpPb5023GeVSysErr");
        histoPCMEMCALEtaToPi0Stat->Write("histoRatioEtaToPi0PCMEMCALpPb5023GeVStatErr");
        graphPCMEMCALEtaToPi0Sys->Write("graphRatioEtaToPi0PCMEMCALpPb5023GeVSysErr");
        graphCombEtaToPi0Tot->Write("graphRatioEtaToPi0CombpPb5023GeVTotErr");
        graphCombEtaToPi0Stat->Write("graphRatioEtaToPi0CombpPb5023GeVStatErr");
        graphCombEtaToPi0Sys->Write("graphRatioEtaToPi0CombpPb5023GeVSysErr");

        
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


    //  **********************************************************************************************************************
    //  ************************* Saving only fits to final results **********************************************************
    //  **********************************************************************************************************************
    
    TString nameOutputCommonFileFitsOnly    = Form("FitsPaperPPpPb5023GeV_%s.root", dateForOutput.Data());
    TFile fFitsResults(nameOutputCommonFileFitsOnly.Data(), "RECREATE");

        // fits for pi0
        fitInvYieldPi0->Write("TsallisFitPi0");
        fitTCMInvYieldPi0->Write("TwoComponentModelFitPi0");
           
        // fits for eta
        fitInvYieldEta->Write("TsallisFitEta");
        fitTCMInvYieldEta->Write("TwoComponentModelFitEta");
                
    fFitsResults.Close();

    // **********************************************************************************************************************
    // ************************* Saving comparison to comb for diff measurements ********************************************
    // **********************************************************************************************************************
        
    TString nameOutputCommonFileCompOnly    = Form("ComparisonsPaperPPb5023GeV_%s.root", dateForOutput.Data());
        
    TFile fCompResults(nameOutputCommonFileCompOnly.Data(), "RECREATE");
        
        graphRatioPi0CombCombFitStat->Write("Pi0_RatioCombToCombFit_Stat");
        graphRatioPi0CombCombFitSys->Write("Pi0_RatioCombToCombFit_Syst");
        graphRatioPi0PCMCombFitStat->Write("Pi0_RatioPCMToCombFit_Stat");
        graphRatioPi0PCMCombFitSys->Write("Pi0_RatioPCMToCombFit_Syst");
        graphRatioPi0PHOSCombFitStat->Write("Pi0_RatioPHOSToCombFit_Stat");
        graphRatioPi0PHOSCombFitSys->Write("Pi0_RatioPHOSToCombFit_Syst");
        graphRatioPi0EMCALCombFitStat->Write("Pi0_RatioEMCToCombFit_Stat");
        graphRatioPi0EMCALCombFitSys->Write("Pi0_RatioEMCToCombFit_Syst");
        graphRatioPi0PCMEMCALCombFitStat->Write("Pi0_RatioPCMEMCToCombFit_Stat");
        graphRatioPi0PCMEMCALCombFitSys->Write("Pi0_RatioPCMEMCToCombFit_Syst");
        graphRatioPi0DalitzCombFitStat->Write("Pi0_RatioDalitzToCombFit_Stat");
        graphRatioPi0DalitzCombFitSys->Write("Pi0_RatioDalitzToCombFit_Syst");

        graphRatioEtaCombCombFitStat->Write("Eta_RatioCombToCombFit_Stat");
        graphRatioEtaCombCombFitSys->Write("Eta_RatioCombToCombFit_Syst");        
        graphRatioEtaPCMCombFitStat->Write("Eta_RatioPCMToCombFit_Stat");
        graphRatioEtaPCMCombFitSys->Write("Eta_RatioPCMToCombFit_Syst");
        graphRatioEtaEMCALCombFitStat->Write("Eta_RatioEMCToCombFit_Stat");
        graphRatioEtaEMCALCombFitSys->Write("Eta_RatioEMCToCombFit_Syst");
        graphRatioEtaPCMEMCALCombFitStat->Write("Eta_RatioPCMEMCToCombFit_Stat");
        graphRatioEtaPCMEMCALCombFitSys->Write("Eta_RatioPCMEMCToCombFit_Syst");
       
    fCompResults.Close();
    
}
    
