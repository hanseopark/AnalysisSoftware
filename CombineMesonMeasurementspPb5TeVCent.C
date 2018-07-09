/****************************************************************************************************************************
******        provided by Gamma Conversion Group, PWGGA,                                                                *****
******        Friederike Bock, friederike.bock@cern.ch                                                                  *****
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

void CombineMesonMeasurementspPb5TeVCent(   TString fileNamePCM             = "",
                                            TString fileNameEMCAL           = "",
                                            TString fileNamePHOS            = "",
                                            TString fileNamePCMEMCAL        = "",
                                            TString fileNamePCMPHOS         = "",
                                            TString fileNameDalitz         = "",
                                            TString fileNamePCMR2           = "",
                                            TString fileNameEMCALR2         = "",
                                            TString fileNamePHOSR2          = "",
                                            TString fileNamePCMEMCALR2      = "",
                                            TString fileNamePCMPHOSR2       = "",
                                            TString suffix                  = "eps",
                                            Bool_t isMC                     = kFALSE,
                                            TString bWCorrection            = "X",
                                            TString fileNameCorrFactors     = "",
                                            TString fileNameInterpolation   = "",
                                            TString fileConfigRpPbErr      = ""
                                        ){

    TString date = ReturnDateString();

    gROOT->Reset();
    gROOT->SetStyle("Plain");
//
    StyleSettingsThesis();
    SetPlotStyle();

    gStyle->SetEndErrorSize(0);

    TString dateForOutput                       = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;

    TString textALICE                           = "ALICE";
    //___________________________________ Declaration of files _____________________________________________
    TString collisionSystempPb                  = "p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
    TString collisionSystempp                   = "pp, #sqrt{#it{s}} = 5.02 TeV";
    TString outputDir                           = Form("%s/%s/CombineMesonMeasurements%spPb5TeVCent", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());
    TString outputDirFile                       = Form("%s/%s/CombineMesonMeasurements%spPb5TeVCent/Inputs", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());
    TString outputDirSupportComb                = Form("%s/%s/CombineMesonMeasurements%spPb5TeVCent/SupportingCombination", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());
    TString outputDirSupport                    = Form("%s/%s/CombineMesonMeasurements%spPb5TeVCent/Supporting", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());

    TString nameFinalResDat                     = Form("%s/CombinedResults%s_FitResults.dat", outputDir.Data(), bWCorrection.Data());
    fstream  fileFitsOutput;
    fileFitsOutput.open(nameFinalResDat.Data(), ios::out);

    TString fileNameTheory                      = "ExternalInputPPb/Theory/TheoryCompilationpPb.root";

    gSystem->Exec("mkdir -p "+outputDirFile);
    gSystem->Exec("mkdir -p "+outputDirSupportComb);
    gSystem->Exec("mkdir -p "+outputDirSupport);
    if (fileNamePCM.CompareTo("")!=0 )          gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDirFile.Data()));
    if (fileNamePCMEMCAL.CompareTo("")!=0 )     gSystem->Exec(Form("cp %s %s/InputPCMEMCAL.root", fileNamePCMEMCAL.Data(), outputDirFile.Data()));
    if (fileNamePCMPHOS.CompareTo("")!=0 )      gSystem->Exec(Form("cp %s %s/InputPCMPHOS.root", fileNamePCMPHOS.Data(), outputDirFile.Data()));
    if (fileNamePHOS.CompareTo("")!=0 )         gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDirFile.Data()));
    if (fileNameEMCAL.CompareTo("")!=0 )        gSystem->Exec(Form("cp %s %s/InputEMCAL.root", fileNameEMCAL.Data(), outputDirFile.Data()));
    if (fileNameDalitz.CompareTo("")!=0 )       gSystem->Exec(Form("cp %s %s/InputDalitz.root", fileNameDalitz.Data(), outputDirFile.Data()));
    if (fileNameInterpolation.CompareTo("")!=0 )gSystem->Exec(Form("cp %s %s/InputInterpolation.root", fileNameInterpolation.Data(), outputDirFile.Data()));
    if (fileNameCorrFactors.CompareTo("")!=0 )  gSystem->Exec(Form("cp %s %s/InputCorrFactors.root", fileNameCorrFactors.Data(), outputDirFile.Data()));
    if (fileConfigRpPbErr.CompareTo("")!=0 )    gSystem->Exec(Form("cp %s %s/configFileRpA.root", fileConfigRpPbErr.Data(), outputDirFile.Data()));

    TString prefix2                             = "";
    if (isMC){
        prefix2                                 = "MC";
    } else {
        prefix2                                 = "Data";
    }

    Int_t maxCentRun1                           = 5;
    Int_t maxCentRun2                           = 4;
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
    TString  centArray[16]                      = { "0-20%", "20-40%", "40-60%", "60-100%", "0-100%",   "0-5%", "0-20%", "0-100%", "60-100%", "10-20%",
                                                    "80-100%", "20-40%", "40-60%", "5-10%", "0-10%",    "60-80%"};
    TString  centArrayOutput[16]                = { "0020", "2040", "4060", "60100", "00100",           "0005", "0020", "00100",  "60100", "1020",
                                                    "80100", "2040", "4060", "0510", "0010",            "6080"};
    TString  centArrayCorr[16]                  = { "0020/", "2040/", "4060/", "60100/", "00100/",      "", "", "",  "", "",
                                                    "", "", "", "", "",                                 ""};
    TString  runArray[16]                       = { "", "", "", "", "",                              "Run2", "Run2", "Run2", "Run2", "Run2",
                                                    "Run2", "Run2", "Run2", "Run2", "Run2",          "Run2" };
    Bool_t  enableCent[16]                      = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,  kTRUE, kTRUE, kTRUE, kTRUE, kFALSE,
                                                    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
    Bool_t  enableCentComb[16]                  = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,  kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                                    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};

    TString  nameMeasGlobal[11]                 = { "PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal",
                                                    "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged",
                                                    "PCMOtherDataset"};
    TString  nameMeasGlobalLabel[11]            = { "PCM", "PHOS", "EMC", "PCM-PHOS", "PCM-EMC",
                                                    "PCM-Dal", "PHOS-Dal", "EMC-Dal", "EMChigh", "mEMC",
                                                    "PCMOtherDataset"};
    TString  nameSecPi0SourceRead[4]            = {"K0S", "K0L", "Lambda", "Rest"};
    TString  nameSecPi0SourceLabel[4]           = {"K^{0}_{s}", "K^{0}_{l}", "#Lambda", "had. int."};
    Double_t maxSecCorr[4]                      = { 0.05, 0.007, 0.0003, 0.03};
    Double_t scalingToNSD                       = 0.964;
    Bool_t scaleNSD[16][11]                     = {
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-20
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 20-40
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 40-60
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 60-100
                                                    { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 },  // MB R1
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-5
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-20
                                                    { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1  },  // MB R2
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 60-100 R2
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 10-20 R2
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 80-100 R2
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 20-40 R2
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 40-60 R2
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 5-10 R2
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-10 R2
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }  // 60-80 R2
                                                  };
    Double_t branchingRatioPi0[11]              =  { 0.98823,  0.98823,  0.98823,  0.98823,  0.98823,
                                                     0.01174,  0.01174,  0.01174,  0.98823,  0.98823,
                                                     0.98823 };
    Double_t branchingRatioEta[11]              =  { 0.3931,  0.3931,  0.3931,  0.3931,  0.3931,
                                                     0.000068,  0.000068,  0.000068,  0.3931,  0.3931,
                                                     0.3931 };

    Double_t minPtPi0Plotting                   = 0.23;
    Double_t maxPtPi0Plotting                   = 51.;
    Double_t minPtEtaPlotting                   = 0.43;
    Double_t maxPtEtaPlotting                   = 31.;

    Color_t  colorDet[11];
    Color_t  colorDetMC[11];
    Style_t  markerStyleDet[11];
    Style_t  markerStyleDetMC[11];
    Size_t   markerSizeDet[11];
    Size_t   markerSizeDetMC[11];
    for (Int_t meth = 0; meth< 11; meth++){
        colorDet[meth]                          = GetDefaultColorDiffDetectors(nameMeasGlobal[meth].Data(), kFALSE, kFALSE, kTRUE);
        colorDetMC[meth]                        = GetDefaultColorDiffDetectors(nameMeasGlobal[meth].Data(), kTRUE, kFALSE, kTRUE);
        markerStyleDet[meth]                    = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[meth].Data(), kFALSE);
        markerStyleDetMC[meth]                  = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[meth].Data(), kTRUE);
        markerSizeDet[meth]                     = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[meth].Data(), kFALSE)*2;
        markerSizeDetMC[meth]                   = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[meth].Data(), kTRUE)*2;
    }

    Color_t  colorCent[16];
    Color_t  colorCentMC[16];
    Style_t  markerStyleCent[16];
    Style_t  markerStyleCentMC[16];
    Size_t   markerSizeCent[16];
    Size_t   markerSizeCentMC[16];
    Double_t nCollpPb[16];
    Double_t nCollErrpPb[16];
    Double_t tpPb[16];
    Double_t tpPbErr[16];
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        colorCent[cent]                         = GetColorDefaultColor("pPb_5.023TeV", "", centArray[cent]);
        colorCentMC[cent]                       = GetColorDefaultColor("pPb_5.023TeV", "HIJING", centArray[cent]);
        markerStyleCent[cent]                   = GetDefaultMarkerStyle("pPb_5.023TeV", "", centArray[cent]);
        markerStyleCentMC[cent]                 = GetDefaultMarkerStyle("pPb_5.023TeV", "HIJING", centArray[cent]);
        markerSizeCent[cent]                    = GetDefaultMarkerSize("pPb_5.023TeV", "", centArray[cent])*2;
        markerSizeCentMC[cent]                  = GetDefaultMarkerSize("pPb_5.023TeV", "HIJING", centArray[cent])*2;
        nCollpPb[cent]                          = GetNCollFromName(centArray[cent], "pPb_5.023TeV");
        nCollErrpPb[cent]                       = GetNCollErrFromName(centArray[cent], "pPb_5.023TeV");
        tpPb[cent]                              = GetTAAFromName(centArray[cent], "pPb_5.023TeV");
        tpPbErr[cent]                           = GetTAAErrFromName(centArray[cent], "pPb_5.023TeV");
    }
    Double_t xSection5TeV                       = ReturnCorrectXSection("5TeV", 1);; // option is wrong fix when possible
    Double_t xSection5TeVErr                    = xSection5023GeVINELErr*1e-3;

    Color_t  colorCGC                           = kGreen+2;
    Color_t  colorDSS                           = kAzure+2;
    Color_t  colorDSSBand                       = kAzure-9;
    Color_t  colorDSSnPDF                       = kOrange+4;
    Color_t  colorDSSnPDFBand                   = kOrange+5;
    Color_t  colorDSSnPDFEPPS                   = kViolet+3;
    Color_t  colorDSSnPDFEPPSBand               = kViolet-8;
    Color_t  colorEPOS3                         = kGreen-2;
    Color_t  colorMcGill                        = kRed-6;
    Color_t  colorDPMJet                        = kBlue+2;
    Color_t  colorHIJING                        = kBlue-7;
    Style_t  styleLineDPMJet                    = 7;
    Style_t  styleLineHIJING                    = 8;
    Style_t  styleLineEPOS3                     = 1;
    Style_t  styleMarkerDSS                     = 24;
    Style_t  styleMarkerNLOMuOne                = 27;
    Style_t  styleMarkerNLOMuTwo                = 30;
    Style_t  styleLineCGC                       = 2;
    Style_t  styleLineDSS                       = 5;
    Style_t  styleLineDSSnPDF                   = 8;
    Style_t  styleLineDSSnPDFEPPS               = 6;
    Style_t  styleLineNLOMuOne                  = 7;
    Style_t  styleLineNLOMuTwo                  = 4;
    Size_t   sizeMarkerNLO                      = 1;
    Width_t  widthLineNLO                       = 2.;

    //***********************************************************************************************
    //************************** Definition of final pt binning (has to be set manually) ************
    //***********************************************************************************************
    cout << "Setting Pi0 binning" << endl;
    Double_t xPtLimitsPi0[16][100];
    Int_t maxNBinsPi0Abs                            = 0;
    Int_t maxNBinsPi0[16]                           = {0};

    cout << "Setting Eta binning" << endl;
    Double_t xPtLimitsEta[16][100];
    Int_t maxNBinsEtaAbs                            = 0;
    Int_t maxNBinsEta[16]                           = {0};

    Int_t nTotMeasPi0[16]                           = {0};
    Int_t nTotMeasEta[16]                           = {0};
    TString fileNamesMethod[11]                     = {"", "", "", "", "", "", "", "", "", "", ""};
    TString fileNamespPbPi0DetailedSys[11]          = {"", "", "", "", "", "", "", "", "", "", ""};
    TString fileNamesppPi0DetailedSys[11]           = {"", "", "", "", "", "", "", "", "", "", ""};
    TString fileNamesRpPbPi0DetailedSys[11]         = {"", "", "", "", "", "", "", "", "", "", ""};
    TString fileNamespPbEtaDetailedSys[11]          = {"", "", "", "", "", "", "", "", "", "", ""};
    TString fileNamesppEtaDetailedSys[11]           = {"", "", "", "", "", "", "", "", "", "", ""};
    TString fileNamesRpPbEtaDetailedSys[11]         = {"", "", "", "", "", "", "", "", "", "", ""};
    Bool_t havePi0SysDetailedpPb[11]                = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
    Bool_t havePi0SysDetailedpp[11]                 = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
    Bool_t haveEtaSysDetailedpPb[11]                = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
    Bool_t haveEtaSysDetailedpp[11]                 = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
    Bool_t haveEffSecCorr[16][4][11];
    vector<TString>* ptSysRemNames                  = new vector<TString>[11];
    TFile* fileMethod[11];
    TDirectory* directoryPi0[16][11];
    TGraphAsymmErrors* graphPi0EffSecCorrFromX[16][4][11];
    TGraphAsymmErrors* graphPi0Eff[16][11];
    TGraphAsymmErrors* graphPi0Acc[16][11];
    TGraphAsymmErrors* graphPi0EffTimesAcc[16][11];
    TGraphAsymmErrors* graphPi0Mass[16][11];
    TGraphAsymmErrors* graphPi0MassMC[16][11];
    TGraphAsymmErrors* graphPi0Width[16][11];
    TGraphAsymmErrors* graphPi0WidthMC[16][11];
    TH1D* histoPi0InvYieldStat[16][11];
    TGraphAsymmErrors* graphPi0InvYieldStat[16][11];
    TGraphAsymmErrors* graphPi0InvYieldSys[16][11];
    TDirectory* directoryEta[16][11];
    TGraphAsymmErrors* graphEtaEff[16][11];
    TGraphAsymmErrors* graphEtaAcc[16][11];
    TGraphAsymmErrors* graphEtaEffTimesAcc[16][11];
    TGraphAsymmErrors* graphEtaMass[16][11];
    TGraphAsymmErrors* graphEtaMassMC[16][11];
    TGraphAsymmErrors* graphEtaWidth[16][11];
    TGraphAsymmErrors* graphEtaWidthMC[16][11];
    TH1D* histoEtaInvYieldStat[16][11];
    TGraphAsymmErrors* graphEtaInvYieldStat[16][11];
    TGraphAsymmErrors* graphEtaInvYieldSys[16][11];
    TH1D* histoEtaToPi0Stat[16][11];
    TGraphAsymmErrors* graphEtaToPi0Stat[16][11];
    TGraphAsymmErrors* graphEtaToPi0Sys[16][11];

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        maxNBinsPi0[cent]                           = GetBinning( xPtLimitsPi0[cent], maxNBinsPi0Abs, "Pi0", Form("pPb_5.023TeV%s",runArray[cent].Data()), 20, -1, kFALSE, centArray[cent]);
        cout << "Pi0 "<< centArray[cent] << ": ";
        for (Int_t i = 0; i< maxNBinsPi0[cent]+1; i++){
            cout << xPtLimitsPi0[cent][i] << ", " ;
        }
        cout << endl;
        for (Int_t meth = 0; meth < 11; meth++){
            if (cent == 0) fileMethod[meth]             = NULL;
            directoryPi0[cent][meth]                    = NULL;
            graphPi0Eff[cent][meth]                     = NULL;
            graphPi0Acc[cent][meth]                     = NULL;
            graphPi0EffTimesAcc[cent][meth]             = NULL;
            graphPi0Mass[cent][meth]                    = NULL;
            graphPi0MassMC[cent][meth]                  = NULL;
            graphPi0Width[cent][meth]                   = NULL;
            graphPi0WidthMC[cent][meth]                 = NULL;
            histoPi0InvYieldStat[cent][meth]            = NULL;
            graphPi0InvYieldStat[cent][meth]            = NULL;
            graphPi0InvYieldSys[cent][meth]             = NULL;
            for (Int_t k = 0; k < 4; k++){
                graphPi0EffSecCorrFromX[cent][k][meth]  = NULL;
                haveEffSecCorr[cent][k][meth]           = kFALSE;
            }
        }
        maxNBinsEta[cent]                           = GetBinning( xPtLimitsEta[cent], maxNBinsEtaAbs, "Eta", Form("pPb_5.023TeV%s",runArray[cent].Data()), 20, -1, kFALSE, centArray[cent]);
        cout << "Eta "<< centArray[cent] << ": ";
        for (Int_t i = 0; i< maxNBinsEta[cent]+1; i++){
            cout << xPtLimitsEta[cent][i] << ", " ;
        }
        cout << endl;
        for (Int_t meth = 0; meth < 11; meth++){
            if (cent == 0) fileMethod[meth]             = NULL;
            directoryEta[cent][meth]                    = NULL;
            graphEtaEff[cent][meth]                     = NULL;
            graphEtaAcc[cent][meth]                     = NULL;
            graphEtaEffTimesAcc[cent][meth]             = NULL;
            graphEtaMass[cent][meth]                    = NULL;
            graphEtaMassMC[cent][meth]                  = NULL;
            graphEtaWidth[cent][meth]                   = NULL;
            graphEtaWidthMC[cent][meth]                 = NULL;
            histoEtaInvYieldStat[cent][meth]            = NULL;
            graphEtaInvYieldStat[cent][meth]            = NULL;
            graphEtaInvYieldSys[cent][meth]             = NULL;
            histoEtaToPi0Stat[cent][meth]               = NULL;
            graphEtaToPi0Stat[cent][meth]               = NULL;
            graphEtaToPi0Sys[cent][meth]                = NULL;
        }
    }

    // *******************************************************************************************************
    // ********************* Set file names for inputs *******************************************************
    // *******************************************************************************************************
    if (fileNamePCM.CompareTo("")!=0)       fileNamesMethod[0]    = fileNamePCM;
    if (fileNamePHOS.CompareTo("")!=0)      fileNamesMethod[1]    = fileNamePHOS;
    if (fileNameEMCAL.CompareTo("")!=0)     fileNamesMethod[2]    = fileNameEMCAL;
    if (fileNamePCMPHOS.CompareTo("")!=0)   fileNamesMethod[3]    = fileNamePCMPHOS;
    if (fileNamePCMEMCAL.CompareTo("")!=0)  fileNamesMethod[4]    = fileNamePCMEMCAL;
    if (fileNameDalitz.CompareTo("")!=0)    fileNamesMethod[5]    = fileNameDalitz;

    if (fileConfigRpPbErr.CompareTo("") != 0){
        ifstream inPbConfig(fileConfigRpPbErr.Data());
        Int_t nReadSys  = 0;
        while(!inPbConfig.eof() && nReadSys<11 ){
            cout << "read line:" <<  nReadSys << "\t"<< nameMeasGlobalLabel[nReadSys].Data() << endl;
            TString nameCurrentSys  = "";
            inPbConfig >> nameCurrentSys ;
            if (nameCurrentSys.CompareTo(nameMeasGlobalLabel[nReadSys].Data()) != 0){
                cout << "wrong order in configuration file, was expecting " <<  nameMeasGlobalLabel[nReadSys].Data() << endl;
                return;
            } else {
                inPbConfig >> fileNamespPbPi0DetailedSys[nReadSys] >> fileNamesppPi0DetailedSys[nReadSys] >> fileNamespPbEtaDetailedSys[nReadSys] >> fileNamesppEtaDetailedSys[nReadSys];
                if (fileNamespPbPi0DetailedSys[nReadSys].CompareTo("bla") != 0){
                    cout << "pi0 pPb sys detailed: " << fileNamespPbPi0DetailedSys[nReadSys] << endl;
                    havePi0SysDetailedpPb[nReadSys]         = kTRUE;
                } else {
                    fileNamespPbPi0DetailedSys[nReadSys]    = "";
                }
                if (fileNamesppPi0DetailedSys[nReadSys].CompareTo("bla") != 0){
                    cout << "pi0 pp sys detailed: " << fileNamesppPi0DetailedSys[nReadSys] << endl;
                    havePi0SysDetailedpp[nReadSys]          = kTRUE;
                } else {
                    fileNamesppPi0DetailedSys[nReadSys]     = "";
                }
                if (fileNamespPbEtaDetailedSys[nReadSys].CompareTo("bla") != 0){
                    cout << "eta pPb sys detailed: " << fileNamespPbEtaDetailedSys[nReadSys] << endl;
                    haveEtaSysDetailedpPb[nReadSys]         = kTRUE;
                } else {
                    fileNamespPbEtaDetailedSys[nReadSys]    = "";
                }
                if (fileNamesppEtaDetailedSys[nReadSys].CompareTo("bla") != 0){
                    cout << "eta pp sys detailed: " << fileNamesppEtaDetailedSys[nReadSys] << endl;
                    haveEtaSysDetailedpp[nReadSys]          = kTRUE;
                } else {
                    fileNamesppEtaDetailedSys[nReadSys]     = "";
                }
                if (!(havePi0SysDetailedpPb[nReadSys] && havePi0SysDetailedpp[nReadSys])){
                    havePi0SysDetailedpPb[nReadSys]         = kTRUE;
                    havePi0SysDetailedpp[nReadSys]          = kTRUE;
                    fileNamespPbPi0DetailedSys[nReadSys]    = "";
                    fileNamesppPi0DetailedSys[nReadSys]     = "";
                } else {
                    fileNamesRpPbPi0DetailedSys[nReadSys]  = Form("%s/Pi0RpPb_%s_detailedSys.dat", outputDir.Data(),nameMeasGlobalLabel[nReadSys].Data());
                }
                if (!(haveEtaSysDetailedpPb[nReadSys] && haveEtaSysDetailedpp[nReadSys])){
                    haveEtaSysDetailedpPb[nReadSys]         = kTRUE;
                    haveEtaSysDetailedpp[nReadSys]          = kTRUE;
                    fileNamespPbEtaDetailedSys[nReadSys]    = "";
                    fileNamesppEtaDetailedSys[nReadSys]     = "";
                } else {
                    fileNamesRpPbEtaDetailedSys[nReadSys]  = Form("%s/EtaRpPb_%s_detailedSys.dat", outputDir.Data(),nameMeasGlobalLabel[nReadSys].Data());
                }
                TString currentString = "";
                inPbConfig >> currentString;
                Int_t counter = 0;
                while (currentString.CompareTo("STOP") != 0 && counter < 20 ){
                    ptSysRemNames[nReadSys].push_back(currentString);
                    counter++;
                    inPbConfig >> currentString;
                }
                cout << "need to take out "<< ptSysRemNames[nReadSys].size() << " systematic errors" << endl;
            }
            nReadSys++;
        }
    }


    // *******************************************************************************************************
    // ***************************** Read in data for different methods **************************************
    // *******************************************************************************************************
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent ++){
        if (!enableCent[cent]) continue;
        cout << "**********************************************************************" << endl;
        cout << "Reading in cent " << centArray[cent].Data() << endl;
        cout << "**********************************************************************" << endl;
        if (cent == maxCentRun1){
            cout << "**********************************************************************" << endl;
            cout << "Resetting files for run2" << endl;
            cout << "**********************************************************************" << endl;
            for (Int_t meth = 0; meth < 11; meth++)   fileNamesMethod[meth]   = "";
            if (fileNamePCMR2.CompareTo("")!=0)       fileNamesMethod[0]    = fileNamePCMR2;
            if (fileNamePHOSR2.CompareTo("")!=0)      fileNamesMethod[1]    = fileNamePHOSR2;
            if (fileNameEMCALR2.CompareTo("")!=0)     fileNamesMethod[2]    = fileNameEMCALR2;
            if (fileNamePCMPHOSR2.CompareTo("")!=0)   fileNamesMethod[3]    = fileNamePCMPHOSR2;
            if (fileNamePCMEMCALR2.CompareTo("")!=0)  fileNamesMethod[4]    = fileNamePCMEMCALR2;
            if (fileNamePCMR2.CompareTo("")!=0)  gSystem->Exec(Form("cp %s %s/InputPCMR2.root", fileNamePCMR2.Data(), outputDirFile.Data()));
            if (fileNamePCMEMCALR2.CompareTo("")!=0)  gSystem->Exec(Form("cp %s %s/InputPCMEMCALR2.root", fileNamePCMEMCALR2.Data(), outputDirFile.Data()));
            if (fileNamePCMPHOSR2.CompareTo("")!=0)  gSystem->Exec(Form("cp %s %s/InputPCMPHOSR2.root", fileNamePCMPHOSR2.Data(), outputDirFile.Data()));
            if (fileNamePHOSR2.CompareTo("")!=0)  gSystem->Exec(Form("cp %s %s/InputPHOSR2.root", fileNamePHOSR2.Data(), outputDirFile.Data()));
            if (fileNameEMCALR2.CompareTo("")!=0)  gSystem->Exec(Form("cp %s %s/InputEMCALR2.root", fileNameEMCALR2.Data(), outputDirFile.Data()));
        }

        for (Int_t meth = 0; meth < 11; meth++){
            if (fileNamesMethod[meth].CompareTo("") != 0){
                cout << "Reading in " << nameMeasGlobalLabel[meth].Data() << " from " << fileNamesMethod[meth].Data() << endl;
                if (!fileMethod[meth] || cent == maxCentRun1)
                    fileMethod[meth]                               = new TFile(fileNamesMethod[meth].Data());
                if (!fileMethod[meth]) {
                    cout << "file " << fileNamesMethod[meth].Data() << " not found! Skipping " << nameMeasGlobalLabel[meth].Data()  << endl;
                    continue;
                }
                directoryPi0[cent][meth]                                = (TDirectory*)fileMethod[meth]->Get(Form("Pi0%spPb_5.023TeV%s",centArray[cent].Data(), runArray[cent].Data()));
//                 if (!directoryPi0[cent][meth]) directoryPi0[cent][meth] = (TDirectory*)fileMethod[meth]->Get(Form("Pi0%spPb_5.023TeV%s e",centArray[cent].Data(), runArray[cent].Data()));
                if (!directoryPi0[cent][meth]) {
                    cout << "File doesn't contain directory " << Form("Pi0%spPb_5.023TeV%s",centArray[cent].Data(), runArray[cent].Data()) << "! Skipping " << nameMeasGlobalLabel[meth].Data()  << endl;
                    continue;
                } else {
                    nTotMeasPi0[cent]++;
                    // reading supporting figures
                    graphPi0Mass[cent][meth]                             = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("Pi0_Mass_data");
                    if (graphPi0Mass[cent][meth]) graphPi0Mass[cent][meth]                             = ScaleGraph(graphPi0Mass[cent][meth], 1000.);
                    graphPi0Width[cent][meth]                            = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("Pi0_Width_data");
                    if (graphPi0Width[cent][meth]) graphPi0Width[cent][meth]                            = ScaleGraph(graphPi0Width[cent][meth], 1000.);
                    graphPi0MassMC[cent][meth]                           = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("Pi0_Mass_MC");
                    if (graphPi0MassMC[cent][meth]) graphPi0MassMC[cent][meth]                           = ScaleGraph(graphPi0MassMC[cent][meth], 1000.);
                    graphPi0WidthMC[cent][meth]                          = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("Pi0_Width_MC");
                    if (graphPi0WidthMC[cent][meth]) graphPi0WidthMC[cent][meth]                          = ScaleGraph(graphPi0WidthMC[cent][meth], 1000.);
                    graphPi0Acc[cent][meth]                              = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("AcceptancePi0");
                    graphPi0Eff[cent][meth]                              = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("EfficiencyPi0");
                    graphPi0EffTimesAcc[cent][meth]                      = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("EffTimesAccPi0");
                    if (graphPi0EffTimesAcc[cent][meth]) graphPi0EffTimesAcc[cent][meth] = ScaleGraph( graphPi0EffTimesAcc[cent][meth], branchingRatioPi0[meth]);
                    for (Int_t k = 0; k < 4; k++){
                        graphPi0EffSecCorrFromX[cent][k][meth]           = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get(Form("EffectiveSecondaryPi0CorrFrom%s",nameSecPi0SourceRead[k].Data()));
                        if (graphPi0EffSecCorrFromX[cent][k][meth]){
//                             cout << nameSecPi0SourceRead[k].Data() << endl;
//                             graphPi0EffSecCorrFromX[cent][k][meth]->Print();
                            Int_t nAboveZero                    = 0;
//                             cout << "k " << k << ", method " << meth << endl;
                            for (Int_t j = 0; j< graphPi0EffSecCorrFromX[cent][k][meth]->GetN(); j++){
                                if(graphPi0EffSecCorrFromX[cent][k][meth]->GetY()[j] > 0) nAboveZero++;
                            }
//                             cout << nAboveZero << endl;
                            if (nAboveZero>0){
                                haveEffSecCorr[cent][k][meth]            = kTRUE;
                            } else {
                                graphPi0EffSecCorrFromX[cent][k][meth]   = NULL;
                            }
                        }
                    }
                    // reading yields
                    graphPi0InvYieldStat[cent][meth]                     = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("graphCorrectedYieldPi0");
                    histoPi0InvYieldStat[cent][meth]                     = (TH1D*)directoryPi0[cent][meth]->Get("CorrectedYieldPi0");
                    cout << "Pi0 stat error" << endl;
                    graphPi0InvYieldStat[cent][meth]->Print();
                    graphPi0InvYieldSys[cent][meth]                      = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("Pi0SystError");
                    cout << "Pi0 sys error" << endl;
                    if (graphPi0InvYieldSys[cent][meth]) graphPi0InvYieldSys[cent][meth]->Print();

                    if (scaleNSD[cent][meth]){
                        if (histoPi0InvYieldStat[cent][meth]) histoPi0InvYieldStat[cent][meth]->Scale(scalingToNSD);
                        if (graphPi0InvYieldStat[cent][meth]) graphPi0InvYieldStat[cent][meth]                 = ScaleGraph(graphPi0InvYieldStat[cent][meth],scalingToNSD);
                        if (graphPi0InvYieldSys[cent][meth]) graphPi0InvYieldSys[cent][meth]                  = ScaleGraph(graphPi0InvYieldSys[cent][meth],scalingToNSD);
                    }
                }
                directoryEta[cent][meth]                             = (TDirectory*)fileMethod[meth]->Get(Form("Eta%spPb_5.023TeV%s",centArray[cent].Data(), runArray[cent].Data()));
                if (!directoryEta[cent][meth]) {
                    cout << "File doesn't contain directory " << Form("Eta%spPb_5.023TeV%s",centArray[cent].Data(), runArray[cent].Data()) << "! Skipping " << nameMeasGlobalLabel[meth].Data()  << endl;
                    continue;
                } else {
                    nTotMeasEta[cent]++;
                    // reading supporting figures
                    graphEtaMass[cent][meth]                             = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("Eta_Mass_data");
                    if (graphEtaMass[cent][meth]) graphEtaMass[cent][meth]                             = ScaleGraph(graphEtaMass[cent][meth], 1000.);
                    graphEtaWidth[cent][meth]                            = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("Eta_Width_data");
                    if (graphEtaWidth[cent][meth]) graphEtaWidth[cent][meth]                            = ScaleGraph(graphEtaWidth[cent][meth], 1000.);
                    graphEtaMassMC[cent][meth]                           = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("Eta_Mass_MC");
                    if (graphEtaMassMC[cent][meth]) graphEtaMassMC[cent][meth]                           = ScaleGraph(graphEtaMassMC[cent][meth], 1000.);
                    graphEtaWidthMC[cent][meth]                          = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("Eta_Width_MC");
                    if (graphEtaWidthMC[cent][meth]) graphEtaWidthMC[cent][meth]                          = ScaleGraph(graphEtaWidthMC[cent][meth], 1000.);
                    graphEtaAcc[cent][meth]                              = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("AcceptanceEta");
                    graphEtaEff[cent][meth]                              = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("EfficiencyEta");
                    graphEtaEffTimesAcc[cent][meth]                      = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("EffTimesAccEta");
                    if (graphEtaEffTimesAcc[cent][meth]) graphEtaEffTimesAcc[cent][meth] = ScaleGraph( graphEtaEffTimesAcc[cent][meth], branchingRatioEta[meth]);
                    // reading yields
                    graphEtaInvYieldStat[cent][meth]                     = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("graphCorrectedYieldEta");
                    histoEtaInvYieldStat[cent][meth]                     = (TH1D*)directoryEta[cent][meth]->Get("CorrectedYieldEta");
                    cout << "Eta stat error" << endl;
                    graphEtaInvYieldStat[cent][meth]->Print();
                    graphEtaInvYieldSys[cent][meth]                      = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("EtaSystError");
                    if (scaleNSD[cent][meth]){
                        if (histoEtaInvYieldStat[cent][meth]) histoEtaInvYieldStat[cent][meth]->Scale(scalingToNSD);
                        if (graphEtaInvYieldStat[cent][meth]) graphEtaInvYieldStat[cent][meth]                 = ScaleGraph(graphEtaInvYieldStat[cent][meth],scalingToNSD);
                        if (graphEtaInvYieldSys[cent][meth]) graphEtaInvYieldSys[cent][meth]                  = ScaleGraph(graphEtaInvYieldSys[cent][meth],scalingToNSD);
                    }
                    if (meth == 1){
                        Int_t ptBin = 1;
                        while (histoEtaInvYieldStat[cent][meth]->GetBinCenter(ptBin) < 4){
                            histoEtaInvYieldStat[cent][meth]->SetBinContent(ptBin, 0);
                            ptBin++;
                        }
                    }
                    cout << "Eta sys error" << endl;
                    if (graphEtaInvYieldSys[cent][meth]) graphEtaInvYieldSys[cent][meth]->Print();
                    histoEtaToPi0Stat[cent][meth]                        = (TH1D*)directoryEta[cent][meth]->Get("EtaToPi0StatError");
                    graphEtaToPi0Stat[cent][meth]                        = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("graphEtaToPi0StatError");
                    graphEtaToPi0Sys[cent][meth]                         = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("EtaToPi0SystError");
                    if (meth == 1){
                        Int_t ptBin = 1;
                        if (histoEtaInvYieldStat[cent][meth]){
                            while (histoEtaInvYieldStat[cent][meth]->GetBinCenter(ptBin) < 4 && ptBin < histoEtaInvYieldStat[cent][meth]->GetNbinsX()){
                                histoEtaInvYieldStat[cent][meth]->SetBinContent(ptBin, 0);
                                ptBin++;
                            }
                        }
                        ptBin = 1;
                        if (histoEtaToPi0Stat[cent][meth]){
                            while (histoEtaToPi0Stat[cent][meth]->GetBinCenter(ptBin) < 4 && ptBin < histoEtaToPi0Stat[cent][meth]->GetNbinsX()){
                                histoEtaToPi0Stat[cent][meth]->SetBinContent(ptBin, 0);
                                ptBin++;
                            }
                        }
                    }
                }
            }
        }
    }


    // *******************************************************************************************************
    // ************************** Loading interpolated spectra ***********************************************
    // *******************************************************************************************************
    TGraphAsymmErrors* statErrorCollectionPi0PP[11];
    TGraphAsymmErrors* systErrorCollectionPi0PP[11];
    TGraphAsymmErrors* systErrorUnCorrCollectionPi0PP[11];
    TGraphAsymmErrors* systErrorInterCollectionPi0PP[11];
    TGraphAsymmErrors* statErrorCollectionEtaPP[11];
    TGraphAsymmErrors* systErrorCollectionEtaPP[11];
    TGraphAsymmErrors* systErrorUnCorrCollectionEtaPP[11];
    TGraphAsymmErrors* systErrorInterCollectionEtaPP[11];
    Bool_t haveRefPPPi0[11]                                 = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                                                kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                                                kFALSE };
    Bool_t haveRefPPEta[11]                                 = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                                                kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                                                kFALSE };
    for (Int_t meth = 0; meth< 11; meth++){
        statErrorCollectionPi0PP[meth]         = NULL;
        systErrorCollectionPi0PP[meth]         = NULL;
        systErrorUnCorrCollectionPi0PP[meth]   = NULL;
        systErrorInterCollectionPi0PP[meth]    = NULL;
        statErrorCollectionEtaPP[meth]         = NULL;
        systErrorCollectionEtaPP[meth]         = NULL;
        systErrorUnCorrCollectionEtaPP[meth]   = NULL;
        systErrorInterCollectionEtaPP[meth]    = NULL;
    }

    TFile* fileInterpolation                        = new TFile(fileNameInterpolation.Data());
        for (Int_t meth = 0; meth< 11; meth++){
            statErrorCollectionPi0PP[meth]                 = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionStatErr%s_Pi0_5.023TeV",nameMeasGlobalLabel[meth].Data()));
            systErrorCollectionPi0PP[meth]                 = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionSystErr%s_Pi0_5.023TeV",nameMeasGlobalLabel[meth].Data()));
            systErrorUnCorrCollectionPi0PP[meth]           = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionUnCorrSystErr%s_Pi0_5.023TeV",nameMeasGlobalLabel[meth].Data()));
            systErrorInterCollectionPi0PP[meth]            = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionInterpolSystErr%s_Pi0_5.023TeV",nameMeasGlobalLabel[meth].Data()));
            if (statErrorCollectionPi0PP[meth] && systErrorCollectionPi0PP[meth] && systErrorUnCorrCollectionPi0PP[meth] && systErrorInterCollectionPi0PP[meth])
                haveRefPPPi0[meth]                         = kTRUE;
            if (haveRefPPPi0[meth])
                cout << "found pi0 pp reference for " << nameMeasGlobalLabel[meth].Data() << endl;
            statErrorCollectionEtaPP[meth]                 = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionStatErr%s_Eta_5.023TeV",nameMeasGlobalLabel[meth].Data()));
            systErrorCollectionEtaPP[meth]                 = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionSystErr%s_Eta_5.023TeV",nameMeasGlobalLabel[meth].Data()));
            systErrorUnCorrCollectionEtaPP[meth]           = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionUnCorrSystErr%s_Eta_5.023TeV",nameMeasGlobalLabel[meth].Data()));
            systErrorInterCollectionEtaPP[meth]            = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionInterpolSystErr%s_Eta_5.023TeV",nameMeasGlobalLabel[meth].Data()));
            if (statErrorCollectionEtaPP[meth] && systErrorCollectionEtaPP[meth] && systErrorUnCorrCollectionEtaPP[meth] && systErrorInterCollectionEtaPP[meth])
                haveRefPPEta[meth]                         = kTRUE;
            if (haveRefPPEta[meth])
                cout << "found eta pp reference for " << nameMeasGlobalLabel[meth].Data() << endl;
        }

        TGraphAsymmErrors* graphPPCombPi0Stat           = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionStatErrComb_Pi0_5.023TeV");
        TGraphAsymmErrors* graphPPCombPi0UncorrSys      = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionUnCorrSystErrComb_Pi0_5.023TeV");
        TGraphAsymmErrors* graphPPCombPi0FullSys        = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionSystErrComb_Pi0_5.023TeV");
        TGraphAsymmErrors* graphPPCombPi0InterSys       = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionInterpolSystErrComb_Pi0_5.023TeV");

        TGraphAsymmErrors* graphPPInvYieldCombPi0Stat   = (TGraphAsymmErrors*)graphPPCombPi0Stat->Clone("graphInvYieldStatErrComb_Pi0_5.023TeV");
        TGraphAsymmErrors* graphPPInvYieldCombPi0Sys    = (TGraphAsymmErrors*)graphPPCombPi0FullSys->Clone("graphInvYieldSystErrComb_Pi0_5.023TeV");
        graphPPInvYieldCombPi0Stat                      = ScaleGraph(graphPPCombPi0Stat,1/(xSection5TeV*recalcBarn));
        graphPPInvYieldCombPi0Sys                       = ScaleGraph(graphPPInvYieldCombPi0Sys,1/(xSection5TeV*recalcBarn));
        TGraphAsymmErrors* graphPPCombEtaStat           = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionStatErrComb_Eta_5.023TeV");
        TGraphAsymmErrors* graphPPCombEtaUncorrSys      = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionUnCorrSystErrComb_Eta_5.023TeV");
        TGraphAsymmErrors* graphPPCombEtaFullSys        = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionSystErrComb_Eta_5.023TeV");
        TGraphAsymmErrors* graphPPCombEtaInterSys       = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionInterpolSystErrComb_Eta_5.023TeV");

        TGraphAsymmErrors* graphPPInvYieldCombEtaStat   = (TGraphAsymmErrors*)graphPPCombEtaStat->Clone("graphInvYieldStatErrComb_Eta_5.023TeV");
        TGraphAsymmErrors* graphPPInvYieldCombEtaSys    = (TGraphAsymmErrors*)graphPPCombEtaFullSys->Clone("graphInvYieldSystErrComb_Eta_5.023TeV");
        graphPPInvYieldCombEtaStat                      = ScaleGraph(graphPPCombEtaStat,1/(xSection5TeV*recalcBarn));
        graphPPInvYieldCombEtaSys                       = ScaleGraph(graphPPInvYieldCombEtaSys,1/(xSection5TeV*recalcBarn));



    // *******************************************************************************************************
    // ************************** Comparison of different pp interpolated pi0 spectra ************************
    // *******************************************************************************************************
    // fitting spectrum with intial parameters
    // Two component model fit from Bylinkin
    TF1* fitTCMDecomposedLPi0PP                             = FitObject("tcmlow","twoCompModelPi0_DecL_PP", "Pi0", NULL, 0.3, 2.);
    TF1* fitTCMDecomposedHPi0PP                             = FitObject("tcmhigh","twoCompModelPi0_DecH_PP", "Pi0", NULL, 4, 50.);
    fitTCMDecomposedLPi0PP->SetParameters(graphPPCombPi0Stat->GetY()[2],0.3);
    graphPPCombPi0Stat->Fit(fitTCMDecomposedLPi0PP,"QNRMEX0+","",0.3,0.8);
    graphPPCombPi0Stat->Fit(fitTCMDecomposedHPi0PP,"QNRMEX0+","",3,20);
    fitTCMDecomposedHPi0PP->SetParameters(graphPPCombPi0Stat->GetY()[2],0.8, 2);

    cout << WriteParameterToFile(fitTCMDecomposedLPi0PP)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedLPi0PP)<< endl;
    cout << WriteParameterToFile(fitTCMDecomposedHPi0PP)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedHPi0PP)<< endl;

    Double_t paramTCMPi0NewPP[5]                            = { fitTCMDecomposedLPi0PP->GetParameter(0),fitTCMDecomposedLPi0PP->GetParameter(1),
                                                                fitTCMDecomposedHPi0PP->GetParameter(0),fitTCMDecomposedHPi0PP->GetParameter(1),fitTCMDecomposedHPi0PP->GetParameter(2)};

    //Two component model from Bylinkin
    TF1* fitTCMInvXSetionPi0PP                              = FitObject("tcm","fitTCMInvXSetionPi0PP5TeV","Pi0",graphPPCombPi0Stat,0.3,20. ,paramTCMPi0NewPP,"QNRMEX0+","", kFALSE);
    cout << "fitting pp pi0 spectra" << endl;
    cout << WriteParameterToFile(fitTCMInvXSetionPi0PP)<< endl;
    TF1* fitTsallisInvXSectionPi0PP                         = FitObject("l","fitTsallisInvXSectionPi0PP5TeV","Pi0",graphPPCombPi0Stat,0.3,20. ,NULL,"QNRMEX0+","", kFALSE);
    cout << "fitting pp pi0 spectra" << endl;
    cout << WriteParameterToFile(fitTsallisInvXSectionPi0PP)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitTsallisInvXSectionPi0PP)<< endl;

    fileFitsOutput <<  WriteParameterToFile(fitTCMInvXSetionPi0PP)<< endl;
    TF1* fitPPTCMInvYieldPi0                                = FitObject("tcm","fitTCMInvYieldPi0PP5TeV","Pi0",NULL,0.3,20. ,NULL,"QNRMEX0+","", kFALSE);
    cout << "fitting pp pi0 spectra" << endl;
    fitPPTCMInvYieldPi0->SetParameter(0,graphPPInvYieldCombPi0Stat->GetY()[2]);
    fitPPTCMInvYieldPi0->SetParameter(1,fitTCMInvXSetionPi0PP->GetParameter(1));
    fitPPTCMInvYieldPi0->SetParameter(2,graphPPInvYieldCombPi0Stat->GetY()[2]);
    fitPPTCMInvYieldPi0->SetParameter(3,fitTCMInvXSetionPi0PP->GetParameter(3));
    fitPPTCMInvYieldPi0->SetParameter(4,fitTCMInvXSetionPi0PP->GetParameter(4));
    graphPPInvYieldCombPi0Stat->Fit(fitPPTCMInvYieldPi0,"QNRMEX0+","",0.3,20);
    cout << WriteParameterToFile(fitPPTCMInvYieldPi0)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitPPTCMInvYieldPi0)<< endl;
    fitPPTCMInvYieldPi0->SetRange(0.3,20);

//     TGraphAsymmErrors* graphRatioPPPi0NLODSS14          = (TGraphAsymmErrors*)graphNLODSS14Pi0PP->Clone("graphRatioPPPi0NLODSS14ToFit");
//     graphRatioPPPi0NLODSS14                             = CalculateGraphErrRatioToFit(graphRatioPPPi0NLODSS14, fitPPTCMInvYieldPi0);
//     TGraph* graphRatioPPPi0NLODSS14Center               = (TGraph*)graphNLODSS14Pi0PPCenter->Clone("graphRatioPPPi0NLODSS14CenterToFit");
//     graphRatioPPPi0NLODSS14Center                       = CalculateGraphRatioToFit(graphRatioPPPi0NLODSS14Center, fitPPTCMInvYieldPi0);


    // *************************************************************************************************************
    // Shift graphs in Y direction as well if desired
    // *************************************************************************************************************

    TGraphAsymmErrors* graphRatioPi0IndCombFitStatPP[11]    = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphRatioPi0IndCombFitSysPP[11]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

    TGraphAsymmErrors* graphRatioPi0CombToFitStatPP         = (TGraphAsymmErrors*)graphPPCombPi0Stat->Clone();
    graphRatioPi0CombToFitStatPP                            = CalculateGraphErrRatioToFit(graphRatioPi0CombToFitStatPP, fitTCMInvXSetionPi0PP);
    TGraphAsymmErrors* graphRatioPi0CombToFitSysPP          = (TGraphAsymmErrors*)graphPPCombPi0FullSys->Clone();
    graphRatioPi0CombToFitSysPP                             = CalculateGraphErrRatioToFit(graphRatioPi0CombToFitSysPP, fitTCMInvXSetionPi0PP);

    for (Int_t meth = 0; meth< 11; meth++){
        if (statErrorCollectionPi0PP[meth]){
            graphRatioPi0IndCombFitStatPP[meth]                = (TGraphAsymmErrors*)statErrorCollectionPi0PP[meth]->Clone(Form("RatioPi0%sToCombFitStatPP", nameMeasGlobalLabel[meth].Data()));
            graphRatioPi0IndCombFitStatPP[meth]                = CalculateGraphErrRatioToFit(graphRatioPi0IndCombFitStatPP[meth], fitTCMInvXSetionPi0PP);
        }
        if (systErrorCollectionPi0PP[meth]){
            graphRatioPi0IndCombFitSysPP[meth]                 = (TGraphAsymmErrors*)systErrorCollectionPi0PP[meth]->Clone(Form("RatioPi0%sToCombFitSystPP", nameMeasGlobalLabel[meth].Data()));
            graphRatioPi0IndCombFitSysPP[meth]                 = CalculateGraphErrRatioToFit(graphRatioPi0IndCombFitSysPP[meth], fitTCMInvXSetionPi0PP);
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Plot Ratio of Comb to Fit for pp *****************************************
    // **********************************************************************************************************************
    Int_t textSizeLabelsPixel               = 52;
    TCanvas* canvasRatioToCombFitPP         = new TCanvas("canvasRatioToCombFitPP","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToCombFitPP, 0.08, 0.01, 0.01, 0.125);
    canvasRatioToCombFitPP->SetLogx();

    Double_t textsizeLabelsPP      = 0;
    if (canvasRatioToCombFitPP->XtoPixel(canvasRatioToCombFitPP->GetX2()) <canvasRatioToCombFitPP->YtoPixel(canvasRatioToCombFitPP->GetY1()) ){
        textsizeLabelsPP           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFitPP->XtoPixel(canvasRatioToCombFitPP->GetX2()) ;
    } else {
        textsizeLabelsPP           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFitPP->YtoPixel(canvasRatioToCombFitPP->GetY1());
    }
    cout << textsizeLabelsPP << endl;

    TH2F * histo2DPi0RatioToCombFitPP;
    histo2DPi0RatioToCombFitPP               = new TH2F("histo2DPi0RatioToCombFitPP","histo2DPi0RatioToCombFitPP",1000,minPtPi0Plotting, maxPtPi0Plotting,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DPi0RatioToCombFitPP, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
    histo2DPi0RatioToCombFitPP->GetXaxis()->SetMoreLogLabels();
    histo2DPi0RatioToCombFitPP->GetXaxis()->SetLabelOffset(-0.01);
    //  histo2DPi0RatioToCombFitPP->GetYaxis()->SetRangeUser(-10,10);
    histo2DPi0RatioToCombFitPP->GetYaxis()->SetRangeUser(0.2,1.82);
    histo2DPi0RatioToCombFitPP->Draw("copy");

    ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombToFitStatPP);

    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombToFitSysPP, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
    graphRatioPi0CombToFitSysPP->Draw("E2same");
    DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombToFitStatPP, markerStyleComb, markerSizeComb, colorComb , colorComb);
    graphRatioPi0CombToFitStatPP->Draw("p,same,z");

    DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1., 1.,0.1, kGray+2);
    DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.1, 1.1,0.1, kGray, 7);
    DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.9, 0.9,0.1, kGray, 7);

    TLatex *labelRatioToFitEnergyPP   = new TLatex(0.95, 0.92, collisionSystempp.Data());
    SetStyleTLatex( labelRatioToFitEnergyPP, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitEnergyPP->Draw();
    TLatex *labelRatioToFitALICEPP    = new TLatex(0.95, 0.86, textALICE.Data());
    SetStyleTLatex( labelRatioToFitALICEPP, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitALICEPP->Draw();
    TLatex *labelRatioToFitPi0PP      = new TLatex(0.12, 0.92, "#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRatioToFitPi0PP, textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
    labelRatioToFitPi0PP->Draw();

    canvasRatioToCombFitPP->SaveAs(Form("%s/Pi0_RatioOfCombToCombFit_PP5TeV.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // *******************************************Plot Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************

    canvasRatioToCombFitPP->cd();
    histo2DPi0RatioToCombFitPP->Draw("copy");

    for (Int_t meth = 10; meth > -1 ; meth--){
        if (graphRatioPi0IndCombFitSysPP[meth]){
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0IndCombFitSysPP[meth], markerStyleDetMC[meth] ,markerSizeDetMC[meth]*0.5, colorDet[meth], colorDet[meth], widthLinesBoxes, kTRUE);
            graphRatioPi0IndCombFitSysPP[meth]->Draw("E2same");
        }
        if (graphRatioPi0IndCombFitStatPP[meth]){
            ProduceGraphAsymmWithoutXErrors(graphRatioPi0IndCombFitStatPP[meth]);
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0IndCombFitStatPP[meth], markerStyleDetMC[meth] ,markerSizeDetMC[meth]*0.5, colorDet[meth], colorDet[meth]);
            graphRatioPi0IndCombFitStatPP[meth]->Draw("p,same,z");
        }
    }
    if (graphRatioPi0IndCombFitStatPP[4])graphRatioPi0IndCombFitStatPP[4]->Draw("p,same,z");

    DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1., 1.,0.5, kGray+2);
    DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.1, 1.1,0.5, kGray, 7);
    DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.9, 0.9,0.5, kGray, 7);
    DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.2, 1.2,0.5, kGray, 9);
    DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.8, 0.8,0.5, kGray, 9);

    labelRatioToFitEnergyPP->Draw();
    labelRatioToFitALICEPP->Draw();
    labelRatioToFitPi0PP->Draw();
    histo2DPi0RatioToCombFitPP->Draw("same,axis");

    //****************************** Definition of the Legend ******************************************
    //**************** Row def ************************
    Double_t rowsLegendOnlyPi0RatioPP[4]          = {0.30, 0.25, 0.20, 0.15};
    Double_t rowsLegendOnlyPi0RatioPPAbs[4]       = {0.92, 0.47, 0.38, 0.29 };
    Double_t columnsLegendOnlyPi0RatioPP[6]       = {0.14, 0.26, 0.35, 0.48, 0.7, 0.8};
    Double_t columnsLegendOnlyPi0RatioPPAbs[6]    = {0.215, 0.69, 0.95, 2, 6.6, 9.6};
    Double_t lengthBoxPP                          = 0.2;
    Double_t heightBoxPP                          = 0.06/2;
    //****************** first Column **************************************************
    TLatex *textPCMOnlyRatioPi0PP                 = new TLatex(columnsLegendOnlyPi0RatioPP[0],rowsLegendOnlyPi0RatioPP[1],nameMeasGlobalLabel[0]);
    SetStyleTLatex( textPCMOnlyRatioPi0PP, textSizeLabelsPixel,4, 1, 43);
    textPCMOnlyRatioPi0PP->Draw();
    TLatex *textPHOSOnlyRatioPi0PP                = new TLatex(columnsLegendOnlyPi0RatioPP[0],rowsLegendOnlyPi0RatioPP[2],nameMeasGlobalLabel[1]);
    SetStyleTLatex( textPHOSOnlyRatioPi0PP, textSizeLabelsPixel,4, 1, 43);
    textPHOSOnlyRatioPi0PP->Draw();
    TLatex *textEMCALOnlyRatioPi0PP               = new TLatex(columnsLegendOnlyPi0RatioPP[0],rowsLegendOnlyPi0RatioPP[3],nameMeasGlobalLabel[2]);
    SetStyleTLatex( textEMCALOnlyRatioPi0PP, textSizeLabelsPixel,4, 1, 43);
    textEMCALOnlyRatioPi0PP->Draw();
    TLatex *textPCMEMCALOnlyRatioPi0PP            = new TLatex(columnsLegendOnlyPi0RatioPP[3],rowsLegendOnlyPi0RatioPP[1],nameMeasGlobalLabel[4]);
    SetStyleTLatex( textPCMEMCALOnlyRatioPi0PP, textSizeLabelsPixel,4, 1, 43);
    textPCMEMCALOnlyRatioPi0PP->Draw();
    TLatex *textPCMPHOSOnlyRatioPi0PP            = new TLatex(columnsLegendOnlyPi0RatioPP[3],rowsLegendOnlyPi0RatioPP[2],nameMeasGlobalLabel[3]);
    SetStyleTLatex( textPCMPHOSOnlyRatioPi0PP, textSizeLabelsPixel,4, 1, 43);
    textPCMPHOSOnlyRatioPi0PP->Draw();
    TLatex *textDalitzOnlyRatioPi0PP            = new TLatex(columnsLegendOnlyPi0RatioPP[3],rowsLegendOnlyPi0RatioPP[3],nameMeasGlobalLabel[5]);
    SetStyleTLatex( textDalitzOnlyRatioPi0PP, textSizeLabelsPixel,4, 1, 43);
    textDalitzOnlyRatioPi0PP->Draw();

    //****************** second Column *************************************************
    TLatex *textStatOnlyRatioPi0PP                = new TLatex(columnsLegendOnlyPi0RatioPP[1],rowsLegendOnlyPi0RatioPP[0] ,"stat");
    SetStyleTLatex( textStatOnlyRatioPi0PP, textSizeLabelsPixel,4, 1, 43);
    textStatOnlyRatioPi0PP->Draw();
    TLatex *textSysOnlyRatioPi0PP                 = new TLatex(columnsLegendOnlyPi0RatioPP[2] ,rowsLegendOnlyPi0RatioPP[0],"syst");
    SetStyleTLatex( textSysOnlyRatioPi0PP, textSizeLabelsPixel,4, 1, 43);
    textSysOnlyRatioPi0PP->Draw();
    TLatex *textStatOnlyRatioPi0PP2               = new TLatex(columnsLegendOnlyPi0RatioPP[4],rowsLegendOnlyPi0RatioPP[0] ,"stat");
    SetStyleTLatex( textStatOnlyRatioPi0PP2, textSizeLabelsPixel,4, 1, 43);
    textStatOnlyRatioPi0PP2->Draw();
    TLatex *textSysOnlyRatioPi0PP2                = new TLatex(columnsLegendOnlyPi0RatioPP[5] ,rowsLegendOnlyPi0RatioPP[0],"syst");
    SetStyleTLatex( textSysOnlyRatioPi0PP2, textSizeLabelsPixel,4, 1, 43);
    textSysOnlyRatioPi0PP2->Draw();

    TMarker* markerPCMPi0OnlyRatioPi0PP             = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[0],columnsLegendOnlyPi0RatioPP[1] ,rowsLegendOnlyPi0RatioPP[1],1);
    markerPCMPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[1] ,rowsLegendOnlyPi0RatioPPAbs[1]);
    TMarker* markerPHOSPi0OnlyRatioPi0PP            = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[1], columnsLegendOnlyPi0RatioPP[1] ,rowsLegendOnlyPi0RatioPP[2],1);
    markerPHOSPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[1] ,rowsLegendOnlyPi0RatioPPAbs[2]);
    TMarker* markerEMCALPi0OnlyRatioPi0PP           = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[2], columnsLegendOnlyPi0RatioPP[1] ,rowsLegendOnlyPi0RatioPP[3],1);
    markerEMCALPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[1] ,rowsLegendOnlyPi0RatioPPAbs[3]);
    TMarker* markerPCMEMCALPi0OnlyRatioPi0PP        = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[4], columnsLegendOnlyPi0RatioPP[3] ,rowsLegendOnlyPi0RatioPP[1],1);
    markerPCMEMCALPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[4] ,rowsLegendOnlyPi0RatioPPAbs[1]);
    TMarker* markerPCMPHOSPi0OnlyRatioPi0PP         = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[3], columnsLegendOnlyPi0RatioPP[3] ,rowsLegendOnlyPi0RatioPP[2],1);
    markerPCMPHOSPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[4] ,rowsLegendOnlyPi0RatioPPAbs[2]);
    TMarker* markerDalitzPi0OnlyRatioPi0PP          = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[5], columnsLegendOnlyPi0RatioPP[3] ,rowsLegendOnlyPi0RatioPP[3],1);
    markerDalitzPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[4] ,rowsLegendOnlyPi0RatioPPAbs[3]);

    TBox* boxPCMPi0OnlyRatioPi0PP                 = CreateBoxFromGraph(graphRatioPi0IndCombFitSysPP[0], columnsLegendOnlyPi0RatioPPAbs[2]-0.1*lengthBoxPP , rowsLegendOnlyPi0RatioPPAbs[1]- heightBoxPP,
                                                                     columnsLegendOnlyPi0RatioPPAbs[2]+ 1.7*lengthBoxPP, rowsLegendOnlyPi0RatioPPAbs[1]+ heightBoxPP);
    boxPCMPi0OnlyRatioPi0PP->Draw("l");
    TBox* boxPHOSPi0OnlyRatioPi0PP                = CreateBoxFromGraph(graphRatioPi0IndCombFitSysPP[1], columnsLegendOnlyPi0RatioPPAbs[2]-0.1*lengthBoxPP , rowsLegendOnlyPi0RatioPPAbs[2]- heightBoxPP,
                                                                     columnsLegendOnlyPi0RatioPPAbs[2]+ 1.7*lengthBoxPP, rowsLegendOnlyPi0RatioPPAbs[2]+ heightBoxPP);
    boxPHOSPi0OnlyRatioPi0PP->Draw("l");
    TBox* boxEMCALPi0OnlyRatioPi0PP               = CreateBoxFromGraph(graphRatioPi0IndCombFitSysPP[2], columnsLegendOnlyPi0RatioPPAbs[2]-0.1*lengthBoxPP , rowsLegendOnlyPi0RatioPPAbs[3]- heightBoxPP,
                                                                     columnsLegendOnlyPi0RatioPPAbs[2]+ 1.7*lengthBoxPP, rowsLegendOnlyPi0RatioPPAbs[3]+ heightBoxPP);
    boxEMCALPi0OnlyRatioPi0PP->Draw("l");
    TBox* boxPCMEMCALPi0OnlyRatioPi0PP            = CreateBoxFromGraph(graphRatioPi0IndCombFitSysPP[4], columnsLegendOnlyPi0RatioPPAbs[5]-0.5*lengthBoxPP , rowsLegendOnlyPi0RatioPPAbs[1]- heightBoxPP,
                                                                     columnsLegendOnlyPi0RatioPPAbs[5]+ 17.2*lengthBoxPP, rowsLegendOnlyPi0RatioPPAbs[1]+ heightBoxPP);
    boxPCMEMCALPi0OnlyRatioPi0PP->Draw("l");
    TBox* boxPCMPHOSPi0OnlyRatioPi0PP             = CreateBoxFromGraph(graphRatioPi0IndCombFitSysPP[3], columnsLegendOnlyPi0RatioPPAbs[5]-0.5*lengthBoxPP , rowsLegendOnlyPi0RatioPPAbs[2]- heightBoxPP,
                                                                       columnsLegendOnlyPi0RatioPPAbs[5]+ 17.2*lengthBoxPP, rowsLegendOnlyPi0RatioPPAbs[2]+ heightBoxPP);
    boxPCMPHOSPi0OnlyRatioPi0PP->Draw("l");
    TBox* boxDalitzPi0OnlyRatioPi0PP             = CreateBoxFromGraph(graphRatioPi0IndCombFitSysPP[5], columnsLegendOnlyPi0RatioPPAbs[5]-0.5*lengthBoxPP , rowsLegendOnlyPi0RatioPPAbs[3]- heightBoxPP,
                                                                      columnsLegendOnlyPi0RatioPPAbs[5]+ 17.2*lengthBoxPP, rowsLegendOnlyPi0RatioPPAbs[3]+ heightBoxPP);
    boxDalitzPi0OnlyRatioPi0PP->Draw("l");

    canvasRatioToCombFitPP->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFit_PP5TeV.%s",outputDir.Data(),suffix.Data()));

    // *******************************************************************************************************
    // ************************** Combination of different pi0 measurements **********************************
    // *******************************************************************************************************
    // REMARKS:
    //     - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - extendable to other energies
    //     - offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented

    TH1D* statErrorCollectionPi0[16][11];
    TGraphAsymmErrors* statErrorGraphCollectionPi0[16][11];
    TGraphAsymmErrors* sysErrorCollectionPi0[16][11];
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCent[cent]) continue;
        for (Int_t meth = 0; meth< 11; meth++){
            // definition of array of histograms (NULL - means we have no measurement at this energy for this rec-method)
            // for statistical error and final value from respective method
            statErrorCollectionPi0[cent][meth]          = NULL;
            cout <<cent << "\t"<< meth << "\t" <<  histoPi0InvYieldStat[cent][meth] << endl;
            if (histoPi0InvYieldStat[cent][meth]) statErrorCollectionPi0[cent][meth]      = (TH1D*)histoPi0InvYieldStat[cent][meth]->Clone(Form("statErr%i_%sPi0",cent,nameMeasGlobalLabel[meth].Data()));
            statErrorGraphCollectionPi0[cent][meth]     = NULL;
            if (graphPi0InvYieldStat[cent][meth]) statErrorGraphCollectionPi0[cent][meth] = (TGraphAsymmErrors*)graphPi0InvYieldStat[cent][meth]->Clone(Form("statErrGraph%i_%sPi0", cent,
                                                                                                                                                             nameMeasGlobalLabel[meth].Data()));
            // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
            // for systematic error from respective method
            sysErrorCollectionPi0[cent][meth]           = NULL;
            if (graphPi0InvYieldSys[cent][meth]) sysErrorCollectionPi0[cent][meth]        = (TGraphAsymmErrors*)graphPi0InvYieldSys[cent][meth]->Clone(Form("sysErr%i_%sPi0", cent,
                                                                                                                                                            nameMeasGlobalLabel[meth].Data()));
        }
    }

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsPi0[16][11]             = {
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 0-20
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 20-40
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 40-60
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 60-100
                                            { 0,  6,  0,  0,  0,  3,  0,  0,  0,  0,  0 },  // MB R1
                                            { 0,  6,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 0-5
                                            { 0,  6,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 0-20
                                            { 0,  6, 16,  0,  8,  3,  0,  0,  0,  0,  0 },  // MB R2
                                            { 0,  6,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 60-100 R2
                                            { 0,  6,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 10-20 R2
                                            { 0,  6,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 80-100 R2
                                            { 0,  6,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 20-40 R2
                                            { 0,  6,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 40-60 R2
                                            { 0,  6,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 5-10 R2
                                            { 0,  6,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 0-10 R2
                                            { 0,  6,  1,  1,  1,  3,  0,  0,  0,  0,  0 }  // 60-80 R2
                                           };
    Int_t offSetsPi0Sys[16][11]          = {
                                            { 1,  5,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 0-20
                                            { 1,  5,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 20-40
                                            { 1,  5,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 40-60
                                            { 1,  5,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 60-100
                                            { 1,  6,  9,  3,  6,  4,  0,  0,  0,  21, 0 },  // MB R1
                                            { 1,  7,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 0-5
                                            { 1,  7,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 0-20
                                            { 1,  7, 23,  5, 10,  0,  0,  0,  0,  21, 0 },  // MB R2
                                            { 1,  7,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 60-100
                                            { 1,  7,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 10-20
                                            { 1,  7,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 80-100
                                            { 1,  7,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 20-40
                                            { 1,  7,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 40-60
                                            { 1,  7,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 5-10
                                            { 1,  7,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 0-10
                                            { 1,  7,  8,  2,  5,  0,  0,  0,  0,  21, 0 }   // 60-80
                                           };
    Int_t offSetPi0Shifting[16][11]      = {
                                            { 0,  4,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 0-20
                                            { 0,  4,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 20-40
                                            { 0,  4,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 40-60
                                            { 0,  4,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 60-100
                                            { 0,  5,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // MB R1
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 0-5 R2
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 0-20 R2
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // MB R2
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 60-100 R2
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 10-20
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 80-100
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 20-40
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 40-60
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 5-10
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 0-10
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 }   // 60-80
                                           };
    Int_t nComBinsPi0Shifting[16][11]    = {
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 0-20
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 20-40
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 40-60
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 60-100
                                            { 30, 10, 30-6-7, 20-3,  31, 17,  0,  0,  0,  0, 0 },  // MB R1
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 0-5
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 0-20
                                            { 0, 0, 0, 0,  0, 0,  0,  0,  0,  0, 0 },  // MB R2
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 60-100
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 10-20
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 80-100
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 20-40
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 40-60
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 5-10
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 0-10
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 }  // 60-80
    };
    Double_t minPtPi0[16]                = { 0.4, 0.4, 0.4, 0.4, 0.3,   0.4, 0.4, 0.3, 0.4, 0.4,
                                             0.4, 0.4, 0.4, 0.4, 0.4,   0.4};

    TGraphAsymmErrors* statErrorRelCollectionPi0[16][11];
    TGraphAsymmErrors* sysErrorRelCollectionPi0[16][11];
    TGraph* graphWeightsPi0[16][11];
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        for (Int_t meth = 0; meth< 11; meth++){
            graphWeightsPi0[cent][meth]                 = NULL;
            statErrorRelCollectionPi0[cent][meth]       = NULL;
            sysErrorRelCollectionPi0[cent][meth]        = NULL;
        }
    }
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCent[cent]) continue;
        for (Int_t meth = 0; meth< 11; meth++){
            if (statErrorCollectionPi0[cent][meth]){
                statErrorRelCollectionPi0[cent][meth]   = new TGraphAsymmErrors(statErrorCollectionPi0[cent][meth]);
                statErrorRelCollectionPi0[cent][meth]   = CalculateRelErrUpAsymmGraph( statErrorRelCollectionPi0[cent][meth], Form("relativeStatErrorPi0_%i_%s", cent,nameMeasGlobal[meth].Data()));
            }
            if (sysErrorCollectionPi0[cent][meth])
                sysErrorRelCollectionPi0[cent][meth]    = CalculateRelErrUpAsymmGraph( sysErrorCollectionPi0[cent][meth], Form("relativeSysErrorPi0%i_%s", cent, nameMeasGlobal[meth].Data()));
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************

    TGraphAsymmErrors* graphCombPi0InvYieldStat[16]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldStatWOXErr[16]       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldSys[16]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldTot[16]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldTotUnshi[16]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldStatUnshi[16]        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldSysUnshi[16]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldRelStat[16]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldRelSys[16]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldRelTot[16]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldStat_yShifted[16]    = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldSys_yShifted[16]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldTot_yShifted[16]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioPi0CombCombFitTot[16]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioPi0CombCombFitStat[16]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioPi0CombCombFitSys[16]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioPi0CombCombFitStatWOXErr[16]   = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };

    TGraphAsymmErrors* graphIndPi0InvYieldStatUnshi[16][11];
    TGraphAsymmErrors* graphIndPi0InvYieldSysUnshi[16][11];
    TGraphAsymmErrors* graphIndPi0InvYieldStat[16][11];
    TGraphAsymmErrors* graphIndPi0InvYieldSys[16][11];
    TGraphAsymmErrors* graphIndPi0InvYieldStat_yShifted[16][11];
    TGraphAsymmErrors* graphIndPi0InvYieldSys_yShifted[16][11];
    TGraphAsymmErrors* graphRatioPi0IndCombFitStat[16][11];
    TGraphAsymmErrors* graphRatioPi0IndCombFitSys[16][11];
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        for (Int_t meth = 0; meth< 11; meth++){
            graphIndPi0InvYieldStatUnshi[cent][meth]                = NULL;
            graphIndPi0InvYieldSysUnshi[cent][meth]                 = NULL;
            graphIndPi0InvYieldStat[cent][meth]                     = NULL;
            graphIndPi0InvYieldSys[cent][meth]                      = NULL;
            graphIndPi0InvYieldStat_yShifted[cent][meth]            = NULL;
            graphIndPi0InvYieldSys_yShifted[cent][meth]             = NULL;
            graphRatioPi0IndCombFitStat[cent][meth]                 = NULL;
            graphRatioPi0IndCombFitSys[cent][meth]                  = NULL;
        }
    }
    Double_t paramTCMPi0New[5]                          = { 20.72,0.167, 0.454, 0.629, 3.163};
    Double_t paramGraphPi0[3]                           = {1.0e12, 8., 0.13};
    TF1* fitTCMDecomposedLPi0[16]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitTCMDecomposedHPi0[16]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitTCMInvYieldPi0[16]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitInvYieldPi0[16]                             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitInvYieldPi0Graph[16]                        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitPowInvYieldPi0[16]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };

    // *************************************************************************************************************
    // ************************************** Define plotting environment ******************************************
    // *************************************************************************************************************
    textSizeLabelsPixel                 = 900*0.04;

    TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
    canvasWeights->SetLogx();
    TH2F * histo2DPi0Weights;
    histo2DPi0Weights = new TH2F("histo2DPi0Weights","histo2DPi0Weights",11000,minPtPi0Plotting, maxPtPi0Plotting,1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DPi0Weights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DPi0Weights->GetXaxis()->SetMoreLogLabels();
    histo2DPi0Weights->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelWeightsPi0         = new TLatex(0.95,0.15,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelWeightsPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelWeightsEnergy      = new TLatex(0.95,0.20,collisionSystempPb.Data());
    SetStyleTLatex( labelWeightsEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);


    TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErr->SetLogx();
    TH2F * histo2DRelSysErr;
    histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,minPtPi0Plotting, maxPtPi0Plotting,1000,0,50.0);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErr->GetYaxis()->SetRangeUser(0,39.5);
    histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErr->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelRelSysErrEnergy    = new TLatex(0.15,0.89,collisionSystempPb.Data());
    SetStyleTLatex( labelRelSysErrEnergy, textSizeLabelsPixel, 4, 1, 43);
    TLatex *labelRelSysErrPi0       = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrPi0, textSizeLabelsPixel, 4, 1, 43);

    TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelStatErr->SetLogx();
    TH2F * histo2DRelStatErr;
    histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,minPtPi0Plotting, maxPtPi0Plotting,1000,0,50.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErr->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelRelStatErrEnergy   = new TLatex(0.95,0.89,collisionSystempPb.Data());
    SetStyleTLatex( labelRelStatErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelStatErrPi0      = new TLatex(0.95,0.85,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelStatErrPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    TCanvas* canvasRelTotErr            = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelTotErr->SetLogx();
    TH2F * histo2DRelTotErrPi0;
    histo2DRelTotErrPi0                 = new TH2F("histo2DRelTotErrPi0","histo2DRelTotErrPi0",11000,minPtPi0Plotting, maxPtPi0Plotting,1000,0,50.0);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrPi0, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,39.5);
    histo2DRelTotErrPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrPi0->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelRelTotErrEnergy    = new TLatex(0.95,0.89,collisionSystempPb.Data());
    SetStyleTLatex( labelRelTotErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelTotErrPi0       = new TLatex(0.95,0.85,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCentComb[cent]) continue;
        cout << "start combining pi0 for " << centArray[cent].Data() << " " << runArray[cent].Data() << endl;
        // Declaration & calculation of combined spectrum
        TString fileNamePi0OutputWeighting      = Form("%s/Pi0_WeightingMethod_%s%s.dat",outputDir.Data(), centArrayOutput[cent].Data(), runArray[cent].Data());
        graphCombPi0InvYieldTot[cent]           = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionPi0[cent],    sysErrorCollectionPi0[cent],
                                                                                        xPtLimitsPi0[cent], maxNBinsPi0[cent],
                                                                                        offSetsPi0[cent], offSetsPi0Sys[cent],
                                                                                        graphCombPi0InvYieldStat[cent], graphCombPi0InvYieldSys[cent],
                                                                                        fileNamePi0OutputWeighting, "pPb_5.023TeV", "Pi0", kTRUE,
                                                                                        NULL, fileNameCorrFactors, centArrayCorr[cent]
                                                                                    );


        if (graphCombPi0InvYieldTot[cent] == NULL) {
            cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
            return;
        }
        while (graphCombPi0InvYieldStat[cent]->GetX()[0] < minPtPi0[cent]){
            graphCombPi0InvYieldStat[cent]->RemovePoint(0);
        }
        while (graphCombPi0InvYieldTot[cent]->GetX()[0] < minPtPi0[cent]){
            graphCombPi0InvYieldTot[cent]->RemovePoint(0);
        }
        while (graphCombPi0InvYieldSys[cent]->GetX()[0] < minPtPi0[cent]){
            graphCombPi0InvYieldSys[cent]->RemovePoint(0);
        }
        graphCombPi0InvYieldTot[cent]->Print();

        // Reading weights from output file for plotting
        ifstream fileWeightsPi0Read;
        fileWeightsPi0Read.open(fileNamePi0OutputWeighting,ios_base::in);
        cout << "reading" << fileNamePi0OutputWeighting << endl;
        Double_t xValuesPi0Read[100];
        Double_t weightsPi0Read[11][100];
        Int_t availablePi0Meas[11]    = {   -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1,
                                            -1};
        Int_t nMeasSetPi0             = nTotMeasPi0[cent];
        Int_t nPtBinsPi0Read          = 0;
        while(!fileWeightsPi0Read.eof() && nPtBinsPi0Read < 100){
            TString garbage             = "";
            if (nPtBinsPi0Read == 0){
                fileWeightsPi0Read >> garbage ;//>> availablePi0Meas[0] >> availablePi0Meas[1] >> availablePi0Meas[2] >> availablePi0Meas[3];
                for (Int_t i = 0; i < nMeasSetPi0; i++){
                    fileWeightsPi0Read >> availablePi0Meas[i] ;
                }
                cout << "read following measurements: ";
                for (Int_t i = 0; i < nMeasSetPi0; i++){
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
            graphWeightsPi0[cent][availablePi0Meas[i]]                        = new TGraph(nPtBinsPi0Read,xValuesPi0Read,weightsPi0Read[availablePi0Meas[i]]);
            Int_t bin = 0;
            for (Int_t n = 0; n< nPtBinsPi0Read; n++){
                if (graphWeightsPi0[cent][availablePi0Meas[i]]->GetY()[bin] == 0) graphWeightsPi0[cent][availablePi0Meas[i]]->RemovePoint(bin);
                else bin++;
            }
        }

        // **********************************************************************************************************************
        // ******************************************* Plotting weights method only EMC *****************************************
        // **********************************************************************************************************************
        textSizeLabelsPixel                 = 900*0.04;
        canvasWeights->cd();
        histo2DPi0Weights->Draw("copy");

            TLegend* legendWeights   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetPi0), 32);
            for (Int_t i = 0; i < nMeasSetPi0; i++){
                DrawGammaSetMarkerTGraph(graphWeightsPi0[cent][availablePi0Meas[i]], markerStyleDet[availablePi0Meas[i]], markerSizeDet[availablePi0Meas[i]]*0.5, colorDet[availablePi0Meas[i]] , colorDet[availablePi0Meas[i]]);
                graphWeightsPi0[cent][availablePi0Meas[i]]->Draw("p,same,z");
                legendWeights->AddEntry(graphWeightsPi0[cent][availablePi0Meas[i]],nameMeasGlobalLabel[availablePi0Meas[i]],"p");
            }
            legendWeights->Draw();

            labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
            labelWeightsEnergy->Draw();
            labelWeightsPi0->Draw();

            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.2, 0.2,0.1, kGray, 3);

        canvasWeights->SaveAs(Form("%s/Pi0_Weights_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));

        //  *********************************************************************************************************************
        //  ************************************ Visualize relative errors ******************************************************
        //  *********************************************************************************************************************

        canvasRelSysErr->cd();
        histo2DRelSysErr->GetYaxis()->SetRangeUser(0,39.5);
        histo2DRelSysErr->Draw("copy");

            cout << "\n*\n*\n*\n*\n*\n*\n*\n*\n*\n*\n* \n text label size:" << textSizeLabelsPixel << endl;
            TLegend* legendRelSysErr        = GetAndSetLegend2(0.62, 0.92-(0.04*nMeasSetPi0), 0.95, 0.92, textSizeLabelsPixel);
            for (Int_t i = 0; i < nMeasSetPi0; i++){
                DrawGammaSetMarkerTGraph(sysErrorRelCollectionPi0[cent][availablePi0Meas[i]], markerStyleDet[availablePi0Meas[i]], markerSizeDet[availablePi0Meas[i]]*0.5, colorDet[availablePi0Meas[i]],
                                         colorDet[availablePi0Meas[i]]);
                sysErrorRelCollectionPi0[cent][availablePi0Meas[i]]->Draw("p,same,z");
                legendRelSysErr->AddEntry(sysErrorRelCollectionPi0[cent][availablePi0Meas[i]],nameMeasGlobalLabel[availablePi0Meas[i]],"p");
            }
            legendRelSysErr->Draw();


            labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
            labelRelSysErrEnergy->Draw();
            labelRelSysErrPi0->Draw();

        canvasRelSysErr->SaveAs(Form("%s/Pi0_RelSysErr_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));
        delete legendRelSysErr;
        //  *********************************************************************************************************************
        //  ************************************ Visualize relative errors ******************************************************
        //  *********************************************************************************************************************

        canvasRelStatErr->cd();
        histo2DRelStatErr->GetYaxis()->SetRangeUser(0,39.5);
        histo2DRelStatErr->Draw("copy");
            TLegend* legendRelStatErr       = GetAndSetLegend2(0.14, 0.92-(0.04*nMeasSetPi0), 0.45, 0.92, textSizeLabelsPixel);
            for (Int_t i = 0; i < nMeasSetPi0; i++){
                DrawGammaSetMarkerTGraph(statErrorRelCollectionPi0[cent][availablePi0Meas[i]], markerStyleDet[availablePi0Meas[i]], markerSizeDet[availablePi0Meas[i]]*0.5, colorDet[availablePi0Meas[i]],
                                         colorDet[availablePi0Meas[i]]);
                statErrorRelCollectionPi0[cent][availablePi0Meas[i]]->Draw("p,same,z");
                legendRelStatErr->AddEntry(statErrorRelCollectionPi0[cent][availablePi0Meas[i]],nameMeasGlobalLabel[availablePi0Meas[i]],"p");
                cout << Form("%s stat: %s", centArray[cent].Data(), nameMeasGlobalLabel[availablePi0Meas[i]].Data() ) << endl;
                statErrorRelCollectionPi0[cent][availablePi0Meas[i]]->Print();
            }
            legendRelStatErr->Draw();

            labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
            labelRelStatErrEnergy->Draw();
            labelRelStatErrPi0->Draw();

        canvasRelStatErr->SaveAs(Form("%s/Pi0_RelStatErr_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));

    //     //  *********************************************************************************************************************
    //     //  ************************ Visualize relative total errors of different combination methods Pi0 ***********************
    //     //  *********************************************************************************************************************
    //
        graphCombPi0InvYieldRelStat[cent]   = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldStat[cent], "relativeStatErrorPi0_Method");
        graphCombPi0InvYieldRelSys[cent]    = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldSys[cent], "relativeSysErrorPi0_Method");
        graphCombPi0InvYieldRelTot[cent]    = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldTot[cent], "relativeTotalErrorPi0_Method");

        canvasRelTotErr->cd();
        histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,39.5);
        histo2DRelTotErrPi0->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot[cent], markerStyleComb, markerSizeComb, kBlue+2 , kBlue+2);
            graphCombPi0InvYieldRelTot[cent]->Draw("p,same,z");

            TLegend* legendRelTotErr1       = GetAndSetLegend2(0.20, 0.92-(0.035), 0.45, 0.92, 32);
            legendRelTotErr1->AddEntry(graphCombPi0InvYieldRelTot[cent],"All","p");
            legendRelTotErr1->Draw();


            labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
            labelRelTotErrEnergy->Draw();
            labelRelTotErrPi0->Draw();

        canvasRelTotErr->SaveAs(Form("%s/Pi0_TotErr_Comp_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));
        histo2DRelTotErrPi0->GetYaxis()->SetTitle("Err (%)");
        histo2DRelTotErrPi0->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphCombPi0InvYieldRelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelStat[cent], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
            graphCombPi0InvYieldRelStat[cent]->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelSys[cent], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
            graphCombPi0InvYieldRelSys[cent]->SetLineStyle(7);
            graphCombPi0InvYieldRelSys[cent]->Draw("l,x0,same,e1");

            TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
            legendRelTotErr3->AddEntry(graphCombPi0InvYieldRelTot[cent],"tot","p");
            legendRelTotErr3->AddEntry(graphCombPi0InvYieldRelStat[cent],"stat","l");
            legendRelTotErr3->AddEntry(graphCombPi0InvYieldRelSys[cent],"sys","l");
            legendRelTotErr3->Draw();

            labelRelTotErrEnergy->Draw();
            labelRelTotErrPi0->Draw();

        canvasRelTotErr->SaveAs(Form("%s/Pi0_Reldecomp_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));


        // **********************************************************************************************************************
        // ************************************* Calculating bin shifted spectra & fitting **************************************
        // **********************************************************************************************************************

        // Cloning spectra
        graphCombPi0InvYieldTotUnshi[cent]          = (TGraphAsymmErrors*)graphCombPi0InvYieldTot[cent]->Clone("Pi0Unshifted");
        graphCombPi0InvYieldStatUnshi[cent]         = (TGraphAsymmErrors*)graphCombPi0InvYieldStat[cent]->Clone("Pi0UnshiftedStat");
        graphCombPi0InvYieldSysUnshi[cent]          = (TGraphAsymmErrors*)graphCombPi0InvYieldSys[cent]->Clone("Pi0UnshiftedSys");

        for (Int_t meth = 0; meth< 11; meth++){
            if (statErrorGraphCollectionPi0[cent][meth]){
                graphIndPi0InvYieldStatUnshi[cent][meth]                 = (TGraphAsymmErrors*)statErrorGraphCollectionPi0[cent][meth]->Clone(Form("Pi0UnshiftedStat%s",nameMeasGlobalLabel[meth].Data()));
                graphIndPi0InvYieldStat[cent][meth]                      = (TGraphAsymmErrors*)statErrorGraphCollectionPi0[cent][meth]->Clone(Form("Pi0Stat%s",nameMeasGlobalLabel[meth].Data()));
                graphIndPi0InvYieldStat_yShifted[cent][meth]             = (TGraphAsymmErrors*)statErrorGraphCollectionPi0[cent][meth]->Clone(Form("Pi0YShiftedStat%s",nameMeasGlobalLabel[meth].Data()));
            }
            if (sysErrorCollectionPi0[cent][meth]){
                graphIndPi0InvYieldSysUnshi[cent][meth]                  = (TGraphAsymmErrors*)sysErrorCollectionPi0[cent][meth]->Clone(Form("Pi0UnshiftedSys%s",nameMeasGlobalLabel[meth].Data()));
                graphIndPi0InvYieldSys[cent][meth]                       = (TGraphAsymmErrors*)sysErrorCollectionPi0[cent][meth]->Clone(Form("Pi0Sys%s",nameMeasGlobalLabel[meth].Data()));
                graphIndPi0InvYieldSys_yShifted[cent][meth]              = (TGraphAsymmErrors*)sysErrorCollectionPi0[cent][meth]->Clone(Form("Pi0YShiftedSys%s",nameMeasGlobalLabel[meth].Data()));
            }
        }

        // fitting spectrum with intial parameters
        // Two component model fit from Bylinkin
        fitTCMDecomposedLPi0[cent]           = FitObject("tcmlow",Form("twoCompModelPi0_DecL_%s", centArrayOutput[cent].Data()), "Pi0", NULL, minPtPi0[cent], 2.);
        fitTCMDecomposedHPi0[cent]           = FitObject("tcmhigh",Form("twoCompModelPi0_DecH_%s", centArrayOutput[cent].Data()), "Pi0", NULL, 4, 50.);
        fitTCMDecomposedLPi0[cent]->SetParameters(graphCombPi0InvYieldTot[cent]->GetY()[2],minPtPi0[cent]);
        graphCombPi0InvYieldStat[cent]->Fit(fitTCMDecomposedLPi0[cent],"QNRMEX0+","",minPtPi0[cent],0.8);
        graphCombPi0InvYieldStat[cent]->Fit(fitTCMDecomposedHPi0[cent],"QNRMEX0+","",3,20);

        cout << WriteParameterToFile(fitTCMDecomposedLPi0[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedLPi0[cent])<< endl;
        cout << WriteParameterToFile(fitTCMDecomposedHPi0[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedHPi0[cent])<< endl;

        fitTCMInvYieldPi0[cent]              = FitObject("tcm",Form("fitTCMInvYieldPi0_%s", centArrayOutput[cent].Data()),"Pi0",graphCombPi0InvYieldStat[cent], minPtPi0[cent], 20. , paramTCMPi0New,
                                                         "QNRMEX0+", "", kFALSE);
        fitTCMDecomposedLPi0[cent]->SetParameter(0, fitTCMInvYieldPi0[cent]->GetParameter(0));
        fitTCMDecomposedLPi0[cent]->SetParameter(1, fitTCMInvYieldPi0[cent]->GetParameter(1));
        fitTCMDecomposedHPi0[cent]->SetParameter(0, fitTCMInvYieldPi0[cent]->GetParameter(2));
        fitTCMDecomposedHPi0[cent]->SetParameter(1, fitTCMInvYieldPi0[cent]->GetParameter(3));
        fitTCMDecomposedHPi0[cent]->SetParameter(2, fitTCMInvYieldPi0[cent]->GetParameter(4));

        // Tsallis fit
        fitInvYieldPi0[cent]                 = FitObject("l","fitInvYieldPi0","Pi0",histoPi0InvYieldStat[cent][2],0.3,20.,paramGraphPi0,"QNRMEX0+");
        fitInvYieldPi0Graph[cent]            = (TF1*)fitInvYieldPi0[cent]->Clone("fitInvYieldPi0Graph");

        // *************************************************************************************************************
        // Shift graphs in X direction if desired
        // *************************************************************************************************************
        if(bWCorrection.Contains("X")){
            TF1* fitShiftingPi0            = FitObject("tmpt","ShiftingPi0","Pi0");
            fitShiftingPi0->SetParameters(fitInvYieldPi0[cent]->GetParameter(0),fitInvYieldPi0[cent]->GetParameter(1), fitInvYieldPi0[cent]->GetParameter(2));
    //         TF1* fitShiftingPi0                 = FitObject("tcmpt","ShiftingPi0","Pi0");
    //         fitShiftingPi0->SetParameters(fitTCMInvYieldPi0[cent]->GetParameter(0),fitTCMInvYieldPi0[cent]->GetParameter(1), fitTCMInvYieldPi0[cent]->GetParameter(2), fitTCMInvYieldPi0[cent]->GetParameter(3),fitTCMInvYieldPi0[cent]->GetParameter(4));

            TGraphAsymmErrors* graphCombPi0InvYieldTotNoShift = (TGraphAsymmErrors*) graphCombPi0InvYieldTot[cent]->Clone("Pi0_NoShift");

            graphCombPi0InvYieldTot[cent]          = ApplyXshift(graphCombPi0InvYieldTot[cent], fitShiftingPi0);
            cout << "comb" << endl;
            graphCombPi0InvYieldStat[cent]->Print();
            graphCombPi0InvYieldStat[cent]         = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot[cent],
                                                                                graphCombPi0InvYieldStat[cent],
                                                                                fitShiftingPi0,
                                                                                0, graphCombPi0InvYieldStat[cent]->GetN());
            graphCombPi0InvYieldSys[cent]          = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot[cent],
                                                                                graphCombPi0InvYieldSys[cent],
                                                                                fitShiftingPi0,
                                                                                0, graphCombPi0InvYieldSys[cent]->GetN());
            for (Int_t meth = 0; meth< 11; meth++){
                if (graphIndPi0InvYieldStat[cent][meth]){
                    cout << "shiting stat err of " << nameMeasGlobalLabel[meth].Data();
                    graphIndPi0InvYieldStat[cent][meth]  = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot[cent],
                                                                                graphIndPi0InvYieldStat[cent][meth],
                                                                                fitShiftingPi0,
                                                                                offSetPi0Shifting[cent][meth], nComBinsPi0Shifting[cent][meth]);

                }
                if (graphIndPi0InvYieldSys[cent][meth]){
                    cout << "shiting sys err of " << nameMeasGlobalLabel[meth].Data();
                    graphIndPi0InvYieldSys[cent][meth]   = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot[cent],
                                                                                graphIndPi0InvYieldSys[cent][meth],
                                                                                fitShiftingPi0,
                                                                                offSetPi0Shifting[cent][meth], nComBinsPi0Shifting[cent][meth]);
                }
            }

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
                Double_t *xPointShift= graphCombPi0InvYieldTot[cent]->GetX();
                for (Int_t i=0; i<numberPoints; i++) {
                    graphCombPi0InvYieldTotNoShift->SetPoint(i,xPoint[i],xPointShift[i]/xPoint[i]);
                    graphCombPi0InvYieldTotNoShift->SetPointError(i,xvalueErrLow[i],xvalueErrUp[i],0,0);
                }
                DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldTotNoShift, markerStyleComb, markerSizeComb, colorComb , colorComb);
                graphCombPi0InvYieldTotNoShift->Draw("p same");

                TLatex *labelRatioToFitBinShift   = new TLatex(0.94, 0.91, Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
                SetStyleTLatex( labelRatioToFitBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitBinShift->Draw();
                TLatex *labelRatioToFitALICEBinShift    = new TLatex(0.94, 0.86, textALICE.Data());
                SetStyleTLatex( labelRatioToFitALICEBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitALICEBinShift->Draw();
                TLatex *labelRatioToFitPi0BinShift      = new TLatex(0.94, 0.807, "#pi^{0} #rightarrow #gamma#gamma");
                SetStyleTLatex( labelRatioToFitPi0BinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitPi0BinShift->Draw();

            canvasShift->Update();
            canvasShift->SaveAs(Form("%s/BinShiftCorrection_Pi0_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));

            // *************************************************************************************************************
            // Plot control graphs
            // *************************************************************************************************************

            TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
            DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.09);
            canvasDummy2->SetLogy();
            canvasDummy2->SetLogx();
            TH2F * histo2DDummy2;
            histo2DDummy2               = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.33,31.,1000,1e-9,10e1);
            SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
            histo2DDummy2->DrawCopy();

                for (Int_t meth = 0; meth< 11; meth++){
                    if (graphIndPi0InvYieldSys[cent][meth]){
                        DrawGammaSetMarkerTGraphAsym(graphIndPi0InvYieldSys[cent][meth], markerStyleDet[meth] ,markerSizeDet[meth]/2, colorDet[meth], colorDet[meth], widthLinesBoxes, kTRUE);
                        graphIndPi0InvYieldSys[cent][meth]->Draw("pEsame");
                    }
                }

                DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStatUnshi[cent], 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
                graphCombPi0InvYieldStatUnshi[cent]->Draw("pEsame");
                DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat[cent], 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
                graphCombPi0InvYieldStat[cent]->Draw("pEsame");

                fitTCMInvYieldPi0[cent]->SetLineColor(kBlue+2);
                fitTCMInvYieldPi0[cent]->Draw("same");

            canvasDummy2->Update();
            canvasDummy2->Print(Form("%s/ComparisonShiftedPi0_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete canvasDummy2;
        }

        // *************************************************************************************************************
        // redo fitting after binshifts
        // *************************************************************************************************************
        // Tsallis function
        cout << WriteParameterToFile(fitInvYieldPi0[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitInvYieldPi0[cent])<< endl;
        //Two component model from Bylinkin
        fitTCMInvYieldPi0[cent]        = FitObject("tcm","fitTCMInvYieldPi0pPb5TeVCent","Pi0",graphCombPi0InvYieldTot[cent],0.3,20. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
        cout << WriteParameterToFile(fitTCMInvYieldPi0[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMInvYieldPi0[cent])<< endl;

        fitPowInvYieldPi0[cent]        = FitObject("m","fitPowInvYieldPi0pPb5TeVCent","Pi0",graphCombPi0InvYieldTot[cent],5,40. ,NULL,"QNRMEX0+","", kFALSE);
        cout << WriteParameterToFile(fitPowInvYieldPi0[cent])<< endl;

        // *************************************************************************************************************
        // Shift graphs in Y direction as well if desired
        // *************************************************************************************************************

        if(bWCorrection.Contains("Y") ){
            graphCombPi0InvYieldTot_yShifted[cent]        = (TGraphAsymmErrors*)graphCombPi0InvYieldTotUnshi[cent]->Clone("Pi0YShiftedCombTot");
            graphCombPi0InvYieldTot_yShifted[cent]        =  ApplyYshiftIndividualSpectra( graphCombPi0InvYieldTot_yShifted[cent], fitInvYieldPi0[cent]);
            graphCombPi0InvYieldStat_yShifted[cent]       = (TGraphAsymmErrors*)graphCombPi0InvYieldStatUnshi[cent]->Clone("Pi0YShiftedCombStat");
            graphCombPi0InvYieldStat_yShifted[cent]       =  ApplyYshiftIndividualSpectra( graphCombPi0InvYieldStat_yShifted[cent], fitInvYieldPi0[cent]);
            graphCombPi0InvYieldSys_yShifted[cent]        = (TGraphAsymmErrors*)graphCombPi0InvYieldSysUnshi[cent]->Clone("Pi0YShiftedCombSys");
            graphCombPi0InvYieldSys_yShifted[cent]        =  ApplyYshiftIndividualSpectra( graphCombPi0InvYieldSys_yShifted[cent], fitInvYieldPi0[cent]);

            for (Int_t meth = 0; meth< 11; meth++){
                if (graphIndPi0InvYieldStat_yShifted[cent][meth]){
                    graphIndPi0InvYieldStat_yShifted[cent][meth] = ApplyYshiftIndividualSpectra( graphIndPi0InvYieldStat_yShifted[cent][meth], fitInvYieldPi0[cent]);
                }
                if (graphIndPi0InvYieldSys_yShifted[cent][meth]){
                    graphIndPi0InvYieldSys_yShifted[cent][meth]  = ApplyYshiftIndividualSpectra( graphIndPi0InvYieldSys_yShifted[cent][meth], fitInvYieldPi0[cent]);
                }
            }
        }

        // *************************************************************************************************************
        // Calculate ratios to combined fit
        // *************************************************************************************************************
    //     TH1D* histoRatioPi0DPMJetToFit                      = (TH1D*) histoDPMJetPi0->Clone("histoRatioPi0DPMJetToFit");
    //     histoRatioPi0DPMJetToFit                            = CalculateHistoRatioToFit (histoRatioPi0DPMJetToFit, fitTCMInvYieldPi0[cent]);
    //     histoRatioPi0DPMJetToFit->GetXaxis()->SetRangeUser(0.3, 20);
    //     TH1D* histoRatioPi0HIJINGToFit                      = (TH1D*) histoHIJINGPi0->Clone("histoRatioPi0HIJINGToFit");
    //     histoRatioPi0HIJINGToFit                            = CalculateHistoRatioToFit (histoRatioPi0HIJINGToFit, fitTCMInvYieldPi0[cent]);
    //     histoRatioPi0HIJINGToFit->GetXaxis()->SetRangeUser(0.3, 20);
    //     TH1D* histoRatioPi0EPOS3ToFit                       = (TH1D*) histoEPOS3Pi0->Clone("histoRatioPi0EPOS3ToFit");
    //     histoRatioPi0EPOS3ToFit                             = CalculateHistoRatioToFit (histoRatioPi0EPOS3ToFit, fitTCMInvYieldPi0[cent]);
    //     histoRatioPi0EPOS3ToFit->GetXaxis()->SetRangeUser(0.3, 20);
    //     TGraphErrors* graphRatioPi0EPOS3ToFit               = new TGraphErrors(histoRatioPi0EPOS3ToFit);
    //     while(graphRatioPi0EPOS3ToFit->GetX()[0] < 0.3)
    //         graphRatioPi0EPOS3ToFit->RemovePoint(0);
    //     while(graphRatioPi0EPOS3ToFit->GetX()[graphRatioPi0EPOS3ToFit->GetN()-1] > 20)
    //         graphRatioPi0EPOS3ToFit->RemovePoint(graphRatioPi0EPOS3ToFit->GetN()-1);
    //     TGraphErrors* graphRatioPi0McGillToFit              = (TGraphErrors*)graphMcGillPi0->Clone("graphRatioPi0McGillToFit");
    //     graphRatioPi0McGillToFit                            = CalculateGraphErrRatioToFit(graphRatioPi0McGillToFit, fitTCMInvYieldPi0[cent]);
    //     while(graphRatioPi0McGillToFit->GetX()[0] < 0.3)
    //         graphRatioPi0McGillToFit->RemovePoint(0);


    //     TGraph* graphRatioPi0CGC                            = (TGraph*)graphPi0CGC5TeV->Clone("graphRatioPi0CGC");
    //     graphRatioPi0CGC                                    = CalculateGraphRatioToFit(graphRatioPi0CGC, fitTCMInvYieldPi0[cent]);
    //     TGraphAsymmErrors* graphRatioPi0NLODSS14            = (TGraphAsymmErrors*)graphNLODSS14Pi0->Clone("graphRatioPi0NLODSS14ToFit");
    //     graphRatioPi0NLODSS14                               = CalculateGraphErrRatioToFit(graphRatioPi0NLODSS14, fitTCMInvYieldPi0[cent]);
    //     TGraph* graphRatioPi0NLODSS14Center                 = (TGraph*)graphNLODSS14Pi0Center->Clone("graphRatioPi0NLODSS14CenterToFit");
    //     graphRatioPi0NLODSS14Center                         = CalculateGraphRatioToFit(graphRatioPi0NLODSS14Center, fitTCMInvYieldPi0[cent]);
    //     TGraphAsymmErrors* graphRatioPi0NLODSS14nPDF        = (TGraphAsymmErrors*)graphNLODSS14nPDFPi0->Clone("graphRatioPi0NLODSS14nPDFToFit");
    //     graphRatioPi0NLODSS14nPDF                           = CalculateGraphErrRatioToFit(graphRatioPi0NLODSS14nPDF, fitTCMInvYieldPi0[cent]);
    //     TGraph* graphRatioPi0NLODSS14nPDFCenter             = (TGraph*)graphNLODSS14nPDFPi0Center->Clone("graphRatioPi0NLODSS14nPDFCenterToFit");
    //     graphRatioPi0NLODSS14nPDFCenter                     = CalculateGraphRatioToFit(graphRatioPi0NLODSS14nPDFCenter, fitTCMInvYieldPi0[cent]);
    //     TGraphAsymmErrors* graphRatioPi0NLODSS14nPDFEPPS16  = (TGraphAsymmErrors*)graphNLODSS14nPDFEPPS16Pi0->Clone("graphRatioPi0NLODSS14nPDFEPPS16ToFit");
    //     graphRatioPi0NLODSS14nPDFEPPS16                     = CalculateGraphErrRatioToFit(graphRatioPi0NLODSS14nPDFEPPS16, fitTCMInvYieldPi0[cent]);
    //     TGraph* graphRatioPi0NLODSS14nPDFEPPS16Center       = (TGraph*)graphNLODSS14nPDFEPPS16Pi0Center->Clone("graphRatioPi0NLODSS14nPDFEPPS16CenterToFit");
    //     graphRatioPi0NLODSS14nPDFEPPS16Center               = CalculateGraphRatioToFit(graphRatioPi0NLODSS14nPDFEPPS16Center, fitTCMInvYieldPi0[cent]);
    //
    //
        graphRatioPi0CombCombFitTot[cent]                       = (TGraphAsymmErrors*)graphCombPi0InvYieldTot[cent]->Clone();
        graphRatioPi0CombCombFitTot[cent]                       = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitTot[cent], fitTCMInvYieldPi0[cent]);
        graphRatioPi0CombCombFitStat[cent]                      = (TGraphAsymmErrors*)graphCombPi0InvYieldStat[cent]->Clone();
        graphRatioPi0CombCombFitStat[cent]                      = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitStat[cent], fitTCMInvYieldPi0[cent]);
        graphRatioPi0CombCombFitSys[cent]                       = (TGraphAsymmErrors*)graphCombPi0InvYieldSys[cent]->Clone();
        graphRatioPi0CombCombFitSys[cent]                       = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitSys[cent], fitTCMInvYieldPi0[cent]);

        for (Int_t meth = 0; meth< 11; meth++){
            if (graphIndPi0InvYieldStat[cent][meth]){
                graphRatioPi0IndCombFitStat[cent][meth]              = (TGraphAsymmErrors*)graphIndPi0InvYieldStat[cent][meth]->Clone(Form("RatioPi0%sToCombFitStat", nameMeasGlobalLabel[meth].Data()));
                graphRatioPi0IndCombFitStat[cent][meth]              = CalculateGraphErrRatioToFit(graphRatioPi0IndCombFitStat[cent][meth], fitTCMInvYieldPi0[cent]);
            }
            if (graphIndPi0InvYieldSys[cent][meth]){
                graphRatioPi0IndCombFitSys[cent][meth]              = (TGraphAsymmErrors*)graphIndPi0InvYieldSys[cent][meth]->Clone(Form("RatioPi0%sToCombFitSyst", nameMeasGlobalLabel[meth].Data()));
                graphRatioPi0IndCombFitSys[cent][meth]              = CalculateGraphErrRatioToFit(graphRatioPi0IndCombFitSys[cent][meth], fitTCMInvYieldPi0[cent]);
            }
        }

        // **********************************************************************************************************************
        // ******************************************* Plot Ratio of Comb to Fit ************************************************
        // **********************************************************************************************************************
        textSizeLabelsPixel                 = 54;
        TCanvas* canvasRatioToCombFit       = new TCanvas("canvasRatioToCombFit","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRatioToCombFit, 0.08, 0.01, 0.01, 0.125);
        canvasRatioToCombFit->SetLogx();

            Double_t textsizeLabelspPb      = 0;
            if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
                textsizeLabelspPb           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
            } else {
                textsizeLabelspPb           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
            }
            cout << textsizeLabelspPb << endl;

        TH2F * histo2DPi0RatioToCombFit;
        histo2DPi0RatioToCombFit               = new TH2F("histo2DPi0RatioToCombFit","histo2DPi0RatioToCombFit",1000,minPtPi0Plotting, maxPtPi0Plotting,1000,0.2,4.    );
        SetStyleHistoTH2ForGraphs(histo2DPi0RatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelspPb, textsizeLabelspPb,
                                  0.85*textsizeLabelspPb,textsizeLabelspPb, 0.9, 0.65, 510, 505);
        histo2DPi0RatioToCombFit->GetXaxis()->SetMoreLogLabels();
        histo2DPi0RatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
        histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.2,1.82);
        histo2DPi0RatioToCombFit->Draw("copy");

            ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStat[cent]);

            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys[cent], markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
            graphRatioPi0CombCombFitSys[cent]->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStat[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphRatioPi0CombCombFitStat[cent]->Draw("p,same,z");

            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1., 1.,0.1, kGray+2);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.1, 1.1,0.1, kGray, 7);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.9, 0.9,0.1, kGray, 7);

            TLatex *labelRatioToFitEnergy   = new TLatex(0.95, 0.92, Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
            SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
            labelRatioToFitEnergy->Draw();
            TLatex *labelRatioToFitALICE    = new TLatex(0.95, 0.86, textALICE.Data());
            SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
            labelRatioToFitALICE->Draw();
            TLatex *labelRatioToFitPi0      = new TLatex(0.12, 0.92, "#pi^{0} #rightarrow #gamma#gamma");
            SetStyleTLatex( labelRatioToFitPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
            labelRatioToFitPi0->Draw();

        canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToCombFit_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));

        // **********************************************************************************************************************
        // *******************************************Plot Ratio of Individual meas to Fit ******************************************
        // **********************************************************************************************************************

        canvasRatioToCombFit->cd();
        histo2DPi0RatioToCombFit->Draw("copy");

            for (Int_t meth = 10; meth > -1 ; meth--){
                if (graphRatioPi0IndCombFitSys[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRatioPi0IndCombFitSys[cent][meth], markerStyleDet[meth] ,markerSizeDet[meth]*0.5, colorDet[meth], colorDet[meth], widthLinesBoxes, kTRUE);
                    graphRatioPi0IndCombFitSys[cent][meth]->Draw("E2same");
                }
                if (graphRatioPi0IndCombFitStat[cent][meth]){
                    ProduceGraphAsymmWithoutXErrors(graphRatioPi0IndCombFitStat[cent][meth]);
                    DrawGammaSetMarkerTGraphAsym(graphRatioPi0IndCombFitStat[cent][meth], markerStyleDet[meth] ,markerSizeDet[meth]*0.5, colorDet[meth], colorDet[meth]);
                    graphRatioPi0IndCombFitStat[cent][meth]->Draw("p,same,z");
                }
            }
            if (graphRatioPi0IndCombFitStat[cent][4]) graphRatioPi0IndCombFitStat[cent][4]->Draw("p,same,z");

            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1., 1.,0.5, kGray+2);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.1, 1.1,0.5, kGray, 7);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.9, 0.9,0.5, kGray, 7);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.2, 1.2,0.5, kGray, 9);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.8, 0.8,0.5, kGray, 9);

            labelRatioToFitEnergy->Draw();
            labelRatioToFitALICE->Draw();
            labelRatioToFitPi0->Draw();
            histo2DPi0RatioToCombFit->Draw("same,axis");

            //****************************** Definition of the Legend ******************************************
            //**************** Row def ************************
            Double_t rowsLegendOnlyPi0Ratio[4]          = {0.30, 0.25, 0.20, 0.15};
            Double_t rowsLegendOnlyPi0RatioAbs[4]       = {0.92, 0.47, 0.38, 0.29 };
            Double_t columnsLegendOnlyPi0Ratio[6]       = {0.14, 0.26, 0.35, 0.48, 0.7, 0.8};
            Double_t columnsLegendOnlyPi0RatioAbs[6]    = {0.215, 0.77, 1.15, 2, 9.5, 14.4};
            Double_t lengthBox                          = 0.2;
            Double_t heightBox                          = 0.06/2;
            //****************** first Column **************************************************

            TLatex *textPCMOnlyRatioPi0                 = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[1],nameMeasGlobalLabel[0]);
            SetStyleTLatex( textPCMOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);

            if (graphRatioPi0IndCombFitSys[cent][0]) textPCMOnlyRatioPi0->Draw();
            TLatex *textPHOSOnlyRatioPi0                = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[2],nameMeasGlobalLabel[1]);
            SetStyleTLatex( textPHOSOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
            if (graphRatioPi0IndCombFitSys[cent][1]) textPHOSOnlyRatioPi0->Draw();
            TLatex *textEMCALOnlyRatioPi0               = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[3],nameMeasGlobalLabel[2]);
            SetStyleTLatex( textEMCALOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
            if (graphRatioPi0IndCombFitSys[cent][2]) textEMCALOnlyRatioPi0->Draw();
            TLatex *textPCMEMCALOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[3],rowsLegendOnlyPi0Ratio[1],nameMeasGlobalLabel[4]);
            SetStyleTLatex( textPCMEMCALOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
            if (graphRatioPi0IndCombFitSys[cent][4]) textPCMEMCALOnlyRatioPi0->Draw();
            TLatex *textPCMPHOSOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[3],rowsLegendOnlyPi0Ratio[2],nameMeasGlobalLabel[3]);
            SetStyleTLatex( textPCMPHOSOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
            if (graphRatioPi0IndCombFitSys[cent][3]) textPCMPHOSOnlyRatioPi0->Draw();
            TLatex *textDalitzOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[3],rowsLegendOnlyPi0Ratio[3],nameMeasGlobalLabel[5]);
            SetStyleTLatex( textDalitzOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
            if (graphRatioPi0IndCombFitSys[cent][5] )textDalitzOnlyRatioPi0->Draw();

            //****************** second Column *************************************************
            TLatex *textStatOnlyRatioPi0                = new TLatex(columnsLegendOnlyPi0Ratio[1],rowsLegendOnlyPi0Ratio[0] ,"stat");
            SetStyleTLatex( textStatOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
            textStatOnlyRatioPi0->Draw();
            TLatex *textSysOnlyRatioPi0                 = new TLatex(columnsLegendOnlyPi0Ratio[2] ,rowsLegendOnlyPi0Ratio[0],"syst");
            SetStyleTLatex( textSysOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
            textSysOnlyRatioPi0->Draw();
            TLatex *textStatOnlyRatioPi02               = new TLatex(columnsLegendOnlyPi0Ratio[4],rowsLegendOnlyPi0Ratio[0] ,"stat");
            SetStyleTLatex( textStatOnlyRatioPi02, textSizeLabelsPixel,4, 1, 43);
            textStatOnlyRatioPi02->Draw();
            TLatex *textSysOnlyRatioPi02                = new TLatex(columnsLegendOnlyPi0Ratio[5] ,rowsLegendOnlyPi0Ratio[0],"syst");
            SetStyleTLatex( textSysOnlyRatioPi02, textSizeLabelsPixel,4, 1, 43);
            textSysOnlyRatioPi02->Draw();

            if (graphRatioPi0IndCombFitSys[cent][0]){
                TMarker* markerPCMPi0OnlyRatioPi0           = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[cent][0],columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[1],1);
                markerPCMPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[1]);
            }
            if (graphRatioPi0IndCombFitSys[cent][1]){
                TMarker* markerPHOSPi0OnlyRatioPi0          = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[cent][1], columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[2],1);
                markerPHOSPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[2]);
            }
            if (graphRatioPi0IndCombFitSys[cent][2]){
                TMarker* markerEMCALPi0OnlyRatioPi0         = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[cent][2], columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[3],1);
                markerEMCALPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[3]);
            }
            if (graphRatioPi0IndCombFitSys[cent][4]){
                TMarker* markerPCMEMCALPi0OnlyRatioPi0      = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[cent][4], columnsLegendOnlyPi0Ratio[3] ,rowsLegendOnlyPi0Ratio[1],1);
                markerPCMEMCALPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[4] ,rowsLegendOnlyPi0RatioAbs[1]);
            }
            if (graphRatioPi0IndCombFitSys[cent][3]){
                TMarker* markerPCMPHOSPi0OnlyRatioPi0      = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[cent][3], columnsLegendOnlyPi0Ratio[3] ,rowsLegendOnlyPi0Ratio[2],1);
                markerPCMPHOSPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[4] ,rowsLegendOnlyPi0RatioAbs[2]);
            }
            if (graphRatioPi0IndCombFitSys[cent][5]){
                TMarker* markerDalitzPi0OnlyRatioPi0      = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[cent][5], columnsLegendOnlyPi0Ratio[3] ,rowsLegendOnlyPi0Ratio[3],1);
                markerDalitzPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[4] ,rowsLegendOnlyPi0RatioAbs[3]);
            }

            if (graphRatioPi0IndCombFitSys[cent][0]){
                TBox* boxPCMPi0OnlyRatioPi0                 = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[cent][0], columnsLegendOnlyPi0RatioAbs[2]-0.2*lengthBox , rowsLegendOnlyPi0RatioAbs[1]- heightBox,
                                                                                columnsLegendOnlyPi0RatioAbs[2]+ 1.8*lengthBox, rowsLegendOnlyPi0RatioAbs[1]+ heightBox);
                boxPCMPi0OnlyRatioPi0->Draw("l");
            }
            if (graphRatioPi0IndCombFitSys[cent][1]){
                TBox* boxPHOSPi0OnlyRatioPi0                = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[cent][1], columnsLegendOnlyPi0RatioAbs[2]-0.2*lengthBox , rowsLegendOnlyPi0RatioAbs[2]- heightBox,
                                                                             columnsLegendOnlyPi0RatioAbs[2]+ 1.8*lengthBox, rowsLegendOnlyPi0RatioAbs[2]+ heightBox);
                boxPHOSPi0OnlyRatioPi0->Draw("l");
            }
            if (graphRatioPi0IndCombFitSys[cent][2]){
                TBox* boxEMCALPi0OnlyRatioPi0               = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[cent][2], columnsLegendOnlyPi0RatioAbs[2]-0.2*lengthBox , rowsLegendOnlyPi0RatioAbs[3]- heightBox,
                                                                                columnsLegendOnlyPi0RatioAbs[2]+ 1.8*lengthBox, rowsLegendOnlyPi0RatioAbs[3]+ heightBox);
                boxEMCALPi0OnlyRatioPi0->Draw("l");
            }
            if (graphRatioPi0IndCombFitSys[cent][4]){
                TBox* boxPCMEMCALPi0OnlyRatioPi0            = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[cent][4], columnsLegendOnlyPi0RatioAbs[5]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[1]- heightBox,
                                                                                columnsLegendOnlyPi0RatioAbs[5]+ 27*lengthBox, rowsLegendOnlyPi0RatioAbs[1]+ heightBox);
                boxPCMEMCALPi0OnlyRatioPi0->Draw("l");
            }
            if (graphRatioPi0IndCombFitSys[cent][3]){
                TBox* boxPCMPHOSPi0OnlyRatioPi0             = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[cent][3], columnsLegendOnlyPi0RatioAbs[5]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[2]- heightBox,
                                                                                columnsLegendOnlyPi0RatioAbs[5]+ 27*lengthBox, rowsLegendOnlyPi0RatioAbs[2]+ heightBox);
                boxPCMPHOSPi0OnlyRatioPi0->Draw("l");
            }
            if (graphRatioPi0IndCombFitSys[cent][5]){
                TBox* boxDalitzPi0OnlyRatioPi0             = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[cent][5], columnsLegendOnlyPi0RatioAbs[5]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[3]- heightBox,
                                                                                 columnsLegendOnlyPi0RatioAbs[5]+ 27*lengthBox, rowsLegendOnlyPi0RatioAbs[3]+ heightBox);
                boxDalitzPi0OnlyRatioPi0->Draw("l");
            }
            canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFit_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));
    }

    // *******************************************************************************************************
    // ************************** Combination of different eta measurements **********************************
    // *******************************************************************************************************
    // REMARKS:
    //     - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - extendable to other energies
    //     - offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented

    TH1D* statErrorCollectionEta[16][11];
    TGraphAsymmErrors* statErrorGraphCollectionEta[16][11];
    TGraphAsymmErrors* sysErrorCollectionEta[16][11];
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCent[cent]) continue;
        for (Int_t meth = 0; meth< 11; meth++){
            // definition of array of histograms (NULL - means we have no measurement at this energy for this rec-method)
            // for statistical error and final value from respective method
            statErrorCollectionEta[cent][meth]          = NULL;
            cout <<cent << "\t"<< meth << "\t" <<  histoEtaInvYieldStat[cent][meth] << endl;
            if (histoEtaInvYieldStat[cent][meth]) statErrorCollectionEta[cent][meth]      = (TH1D*)histoEtaInvYieldStat[cent][meth]->Clone(Form("statErr%i_%sEta",cent,nameMeasGlobalLabel[meth].Data()));
            statErrorGraphCollectionEta[cent][meth]     = NULL;
            if (graphEtaInvYieldStat[cent][meth]) statErrorGraphCollectionEta[cent][meth] = (TGraphAsymmErrors*)graphEtaInvYieldStat[cent][meth]->Clone(Form("statErrGraph%i_%sEta", cent,
                                                                                                                                                             nameMeasGlobalLabel[meth].Data()));
            // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
            // for systematic error from respective method
            sysErrorCollectionEta[cent][meth]           = NULL;
            if (graphEtaInvYieldSys[cent][meth]) sysErrorCollectionEta[cent][meth]        = (TGraphAsymmErrors*)graphEtaInvYieldSys[cent][meth]->Clone(Form("sysErr%i_%sEta", cent,
                                                                                                                                                            nameMeasGlobalLabel[meth].Data()));
        }
    }

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Eta
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsEta[16][11]             = {
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-20
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 20-40
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 40-60
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 60-100
                                            { 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0 },  // MB R1
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-5 R2
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-20 R2
                                            { 0,  7,  10,  5,  7,  0,  0,  0,  0,  0,  0 },  // MB R2
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 60-100 R2
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 10-20
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 80-100
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 20-40
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 40-60
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 5-10
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-10
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 }  // 60-80
    };
    Int_t offSetsEtaSys[16][11]          = {
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 },  // 0-20
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 },  // 20-40
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 },  // 40-60
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 },  // 60-100
                                            { 3,  0,  7,  5,  5,  0,  0,  0,  0,  0, 0 },  // MB R1
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 },  // 0-5 R2
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 },  // 0-20 R2
                                            { 2,  11,  13,  8,  10,  0,  0,  0,  0,  0, 0 },  // MB R2
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 },  // 60-100 R2
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 },  // 10-20
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 },  // 80-100
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 },  // 20-40
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 },  // 40-60
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 },  // 5-10
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 },  // 0-10
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 }  // 60-80
                                            };
    Int_t offSetEtaShifting[16][11]      = {
                                            { 0,  0,  3,  1,  2,  0,  0,  0,  0,  21, 0 },  // 0-20
                                            { 0,  0,  3,  1,  2,  0,  0,  0,  0,  21, 0 },  // 20-40
                                            { 0,  0,  3,  1,  2,  0,  0,  0,  0,  21, 0 },  // 40-60
                                            { 0,  0,  3,  1,  2,  0,  0,  0,  0,  21, 0 },  // 60-100
                                            { 0,  0,  4,  2,  2,  0,  0,  0,  0,  0,  0 },  // MB R1
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 0-5
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 0-20
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },   // MB R2
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 60-100
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 10-20
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 80-100
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 20-40
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 40-60
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 5-10
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 0-10
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 }  // 60-80
                                            };
    Int_t nComBinsEtaShifting[16][11]    = {
                                            { 9, 0, 10, 8,  7, 0,  0,  0,  0,  0, 0 },  // 0-20
                                            { 9, 0, 10, 8,  7, 0,  0,  0,  0,  0, 0 },  // 20-40
                                            { 9, 0, 10, 8,  7, 0,  0,  0,  0,  0, 0 },  // 40-60
                                            { 9, 0, 10, 8,  7, 0,  0,  0,  0,  0, 0 },  // 60-100
                                            { 13, 0,  17, 11, 15,  0,  0,  0,  0,  0, 0 },  // MB R1
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 0-5
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 0-20
                                            { 0, 0, 0, 0,  0, 0,  0,  0,  0,  0, 0 },  // MB R2
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 60-100
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 10-20
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 80-100
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 20-40
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 40-60
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 5-10
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 0-10
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 60-80
                                           };
    Double_t minPtEta[16]                = { 0.7, 0.7, 0.7, 0.7, 0.7,   0.4, 0.4, 0.4, 0.4, 0.4,
                                             0.4, 0.4, 0.4, 0.4, 0.3,   0.4};

    TGraphAsymmErrors* statErrorRelCollectionEta[16][11];
    TGraphAsymmErrors* sysErrorRelCollectionEta[16][11];
    TGraph* graphWeightsEta[16][11];
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        for (Int_t meth = 0; meth< 11; meth++){
            graphWeightsEta[cent][meth]                 = NULL;
            statErrorRelCollectionEta[cent][meth]       = NULL;
            sysErrorRelCollectionEta[cent][meth]        = NULL;
        }
    }
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCent[cent]) continue;
        for (Int_t meth = 0; meth< 11; meth++){
            if (statErrorCollectionEta[cent][meth]){
                statErrorRelCollectionEta[cent][meth]   = new TGraphAsymmErrors(statErrorCollectionEta[cent][meth]);
                statErrorRelCollectionEta[cent][meth]   = CalculateRelErrUpAsymmGraph( statErrorRelCollectionEta[cent][meth], Form("relativeStatErrorEta_%i_%s", cent,nameMeasGlobal[meth].Data()));
            }
            if (sysErrorCollectionEta[cent][meth])
                sysErrorRelCollectionEta[cent][meth]    = CalculateRelErrUpAsymmGraph( sysErrorCollectionEta[cent][meth], Form("relativeSysErrorEta%i_%s", cent, nameMeasGlobal[meth].Data()));
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************

    TGraphAsymmErrors* graphCombEtaInvYieldStat[16]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldStatWOXErr[16]       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldSys[16]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldTot[16]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldTotUnshi[16]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldStatUnshi[16]        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldSysUnshi[16]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldRelStat[16]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldRelSys[16]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldRelTot[16]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldStat_yShifted[16]    = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldSys_yShifted[16]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldTot_yShifted[16]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioEtaCombCombFitTot[16]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioEtaCombCombFitStat[16]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioEtaCombCombFitStatWOXErr[16]   = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioEtaCombCombFitSys[16]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };

    TGraphAsymmErrors* graphIndEtaInvYieldStatUnshi[16][11];
    TGraphAsymmErrors* graphIndEtaInvYieldSysUnshi[16][11];
    TGraphAsymmErrors* graphIndEtaInvYieldStat[16][11];
    TGraphAsymmErrors* graphIndEtaInvYieldSys[16][11];
    TGraphAsymmErrors* graphIndEtaInvYieldStat_yShifted[16][11];
    TGraphAsymmErrors* graphIndEtaInvYieldSys_yShifted[16][11];
    TGraphAsymmErrors* graphRatioEtaIndCombFitStat[16][11];
    TGraphAsymmErrors* graphRatioEtaIndCombFitSys[16][11];
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        for (Int_t meth = 0; meth< 11; meth++){
            graphIndEtaInvYieldStatUnshi[cent][meth]                = NULL;
            graphIndEtaInvYieldSysUnshi[cent][meth]                 = NULL;
            graphIndEtaInvYieldStat[cent][meth]                     = NULL;
            graphIndEtaInvYieldSys[cent][meth]                      = NULL;
            graphIndEtaInvYieldStat_yShifted[cent][meth]            = NULL;
            graphIndEtaInvYieldSys_yShifted[cent][meth]             = NULL;
            graphRatioEtaIndCombFitStat[cent][meth]                 = NULL;
            graphRatioEtaIndCombFitSys[cent][meth]                  = NULL;
        }
    }
    Double_t paramTCMEtaNew[5]                          = { 20.72,0.167, 0.454, 0.629, 3.163};
    Double_t paramGraphEta[3]                           = {1.0e12, 8., 0.13};
    TF1* fitTCMDecomposedLEta[16]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitTCMDecomposedHEta[16]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitTCMInvYieldEta[16]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitInvYieldEta[16]                             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitInvYieldEtaGraph[16]                        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitPowInvYieldEta[16]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };

    // *************************************************************************************************************
    // ************************************** Define plotting environment ******************************************
    // *************************************************************************************************************
    textSizeLabelsPixel                 = 900*0.04;
    TH2F * histo2DEtaWeights;
    histo2DEtaWeights = new TH2F("histo2DEtaWeights","histo2DEtaWeights",11000,minPtEtaPlotting, maxPtEtaPlotting,1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DEtaWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DEtaWeights->GetXaxis()->SetMoreLogLabels();
    histo2DEtaWeights->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelWeightsEta         = new TLatex(0.95,0.15,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelWeightsEta, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    TLatex *labelRelSysErrEta       = new TLatex(0.15,0.85,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrEta, textSizeLabelsPixel, 4, 1, 43);

    TLatex *labelRelStatErrEta      = new TLatex(0.95,0.85,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelStatErrEta, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    TH2F * histo2DRelTotErrEta;
    histo2DRelTotErrEta                 = new TH2F("histo2DRelTotErrEta","histo2DRelTotErrEta",11000,minPtEtaPlotting, maxPtEtaPlotting,1000,0,50.0);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrEta, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,39.5);
    histo2DRelTotErrEta->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrEta->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelRelTotErrEta       = new TLatex(0.95,0.85,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrEta, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCentComb[cent]) continue;
        // Declaration & calculation of combined spectrum
        cout << "start combining eta for " << centArray[cent].Data() << " " << runArray[cent].Data() << endl;
        TString fileNameEtaOutputWeighting      = Form("%s/Eta_WeightingMethod_%s%s.dat",outputDir.Data(), centArrayOutput[cent].Data(), runArray[cent].Data());
        graphCombEtaInvYieldTot[cent]           = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionEta[cent],    sysErrorCollectionEta[cent],
                                                                                        xPtLimitsEta[cent], maxNBinsEta[cent],
                                                                                        offSetsEta[cent], offSetsEtaSys[cent],
                                                                                        graphCombEtaInvYieldStat[cent], graphCombEtaInvYieldSys[cent],
                                                                                        fileNameEtaOutputWeighting, "pPb_5.023TeV", "Eta", kTRUE,
                                                                                        NULL, fileNameCorrFactors, centArrayCorr[cent]
                                                                                    );


        if (graphCombEtaInvYieldTot[cent] == NULL) {
            cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
            return;
        }
        while (graphCombEtaInvYieldStat[cent]->GetX()[0] < minPtEta[cent]){
            graphCombEtaInvYieldStat[cent]->RemovePoint(0);
        }
        while (graphCombEtaInvYieldTot[cent]->GetX()[0] < minPtEta[cent]){
            graphCombEtaInvYieldTot[cent]->RemovePoint(0);
        }
        while (graphCombEtaInvYieldSys[cent]->GetX()[0] < minPtEta[cent]){
            graphCombEtaInvYieldSys[cent]->RemovePoint(0);
        }
        graphCombEtaInvYieldTot[cent]->Print();

        // Reading weights from output file for plotting
        ifstream fileWeightsEtaRead;
        fileWeightsEtaRead.open(fileNameEtaOutputWeighting,ios_base::in);
        cout << "reading" << fileNameEtaOutputWeighting << endl;
        Double_t xValuesEtaRead[50];
        Double_t weightsEtaRead[11][50];
        Int_t availableEtaMeas[11]    = {   -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1,
                                            -1};
        Int_t nMeasSetEta             = nTotMeasEta[cent];
        Int_t nPtBinsEtaRead          = 0;
        while(!fileWeightsEtaRead.eof() && nPtBinsEtaRead < 50){
            TString garbage             = "";
            if (nPtBinsEtaRead == 0){
                fileWeightsEtaRead >> garbage ;//>> availableEtaMeas[0] >> availableEtaMeas[1] >> availableEtaMeas[2] >> availableEtaMeas[3];
                for (Int_t i = 0; i < nMeasSetEta; i++){
                    fileWeightsEtaRead >> availableEtaMeas[i] ;
                }
                cout << "read following measurements: ";
                for (Int_t i = 0; i < nMeasSetEta; i++){
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
        nPtBinsEtaRead                  = nPtBinsEtaRead-2 ;
        fileWeightsEtaRead.close();

        for (Int_t i = 0; i < nMeasSetEta; i++){
            graphWeightsEta[cent][availableEtaMeas[i]]                        = new TGraph(nPtBinsEtaRead,xValuesEtaRead,weightsEtaRead[availableEtaMeas[i]]);
            Int_t bin = 0;
            for (Int_t n = 0; n< nPtBinsEtaRead; n++){
                if (graphWeightsEta[cent][availableEtaMeas[i]]->GetY()[bin] == 0) graphWeightsEta[cent][availableEtaMeas[i]]->RemovePoint(bin);
                else bin++;
            }
        }

        // **********************************************************************************************************************
        // ******************************************* Plotting weights method only EMC *****************************************
        // **********************************************************************************************************************
        textSizeLabelsPixel                 = 900*0.04;
        canvasWeights->cd();
        histo2DEtaWeights->Draw("copy");

            TLegend* legendWeights   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEta), 32);
            for (Int_t i = 0; i < nMeasSetEta; i++){
                DrawGammaSetMarkerTGraph(graphWeightsEta[cent][availableEtaMeas[i]], markerStyleDet[availableEtaMeas[i]], markerSizeDet[availableEtaMeas[i]]*0.5, colorDet[availableEtaMeas[i]] , colorDet[availableEtaMeas[i]]);
                graphWeightsEta[cent][availableEtaMeas[i]]->Draw("p,same,z");
                legendWeights->AddEntry(graphWeightsEta[cent][availableEtaMeas[i]],nameMeasGlobalLabel[availableEtaMeas[i]],"p");
            }
            legendWeights->Draw();

            labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
            labelWeightsEnergy->Draw();
            labelWeightsEta->Draw();

            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.2, 0.2,0.1, kGray, 3);

        canvasWeights->SaveAs(Form("%s/Eta_Weights_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));

        //  *********************************************************************************************************************
        //  ************************************ Visualize relative errors ******************************************************
        //  *********************************************************************************************************************

        canvasRelSysErr->cd();
        histo2DRelSysErr->GetYaxis()->SetRangeUser(0,39.5);
        histo2DRelSysErr->Draw("copy");

            TLegend* legendRelSysErrEta        = GetAndSetLegend2(0.62, 0.92-(0.04*nMeasSetEta), 0.95, 0.92, textSizeLabelsPixel);
            for (Int_t i = 0; i < nMeasSetEta; i++){
                DrawGammaSetMarkerTGraph(sysErrorRelCollectionEta[cent][availableEtaMeas[i]], markerStyleDet[availableEtaMeas[i]], markerSizeDet[availableEtaMeas[i]]*0.5, colorDet[availableEtaMeas[i]],
                                         colorDet[availableEtaMeas[i]]);
                sysErrorRelCollectionEta[cent][availableEtaMeas[i]]->Draw("p,same,z");
                legendRelSysErrEta->AddEntry(sysErrorRelCollectionEta[cent][availableEtaMeas[i]],nameMeasGlobalLabel[availableEtaMeas[i]],"p");
            }
            legendRelSysErrEta->Draw();

            labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
            labelRelSysErrEnergy->Draw();
            labelRelSysErrEta->Draw();

        canvasRelSysErr->SaveAs(Form("%s/Eta_RelSysErr_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));

        //  *********************************************************************************************************************
        //  ************************************ Visualize relative errors ******************************************************
        //  *********************************************************************************************************************

        canvasRelStatErr->cd();
        histo2DRelStatErr->GetYaxis()->SetRangeUser(0,39.5);
        histo2DRelStatErr->Draw("copy");
            TLegend* legendRelStatErr       = GetAndSetLegend2(0.14, 0.92-(0.04*nMeasSetEta), 0.45, 0.92, textSizeLabelsPixel);
            for (Int_t i = 0; i < nMeasSetEta; i++){
                DrawGammaSetMarkerTGraph(statErrorRelCollectionEta[cent][availableEtaMeas[i]], markerStyleDet[availableEtaMeas[i]], markerSizeDet[availableEtaMeas[i]]*0.5, colorDet[availableEtaMeas[i]],
                                         colorDet[availableEtaMeas[i]]);
                statErrorRelCollectionEta[cent][availableEtaMeas[i]]->Draw("p,same,z");
                legendRelStatErr->AddEntry(statErrorRelCollectionEta[cent][availableEtaMeas[i]],nameMeasGlobalLabel[availableEtaMeas[i]],"p");
            }
            legendRelStatErr->Draw();

            labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
            labelRelStatErrEnergy->Draw();
            labelRelStatErrEta->Draw();

        canvasRelStatErr->SaveAs(Form("%s/Eta_RelStatErr_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));

    //     //  *********************************************************************************************************************
    //     //  ************************ Visualize relative total errors of different combination methods Eta ***********************
    //     //  *********************************************************************************************************************
    //
        graphCombEtaInvYieldRelStat[cent]   = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldStat[cent], "relativeStatErrorEta_Method");
        graphCombEtaInvYieldRelSys[cent]    = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldSys[cent], "relativeSysErrorEta_Method");
        graphCombEtaInvYieldRelTot[cent]    = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldTot[cent], "relativeTotalErrorEta_Method");

        canvasRelTotErr->cd();
        histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,39.5);
        histo2DRelTotErrEta->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelTot[cent], markerStyleComb, markerSizeComb, kBlue+2 , kBlue+2);
            graphCombEtaInvYieldRelTot[cent]->Draw("p,same,z");

            TLegend* legendRelTotErr1       = GetAndSetLegend2(0.20, 0.92-(0.035), 0.45, 0.92, 32);
            legendRelTotErr1->AddEntry(graphCombEtaInvYieldRelTot[cent],"All","p");
            legendRelTotErr1->Draw();


            labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
            labelRelTotErrEnergy->Draw();
            labelRelTotErrEta->Draw();

        canvasRelTotErr->SaveAs(Form("%s/Eta_TotErr_Comp_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));
        histo2DRelTotErrEta->GetYaxis()->SetTitle("Err (%)");
        histo2DRelTotErrEta->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphCombEtaInvYieldRelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelStat[cent], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
            graphCombEtaInvYieldRelStat[cent]->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelSys[cent], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
            graphCombEtaInvYieldRelSys[cent]->SetLineStyle(7);
            graphCombEtaInvYieldRelSys[cent]->Draw("l,x0,same,e1");

            TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
            legendRelTotErr3->AddEntry(graphCombEtaInvYieldRelTot[cent],"tot","p");
            legendRelTotErr3->AddEntry(graphCombEtaInvYieldRelStat[cent],"stat","l");
            legendRelTotErr3->AddEntry(graphCombEtaInvYieldRelSys[cent],"sys","l");
            legendRelTotErr3->Draw();

            labelRelTotErrEnergy->Draw();
            labelRelTotErrEta->Draw();

        canvasRelTotErr->SaveAs(Form("%s/Eta_Reldecomp_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));


        // **********************************************************************************************************************
        // ************************************* Calculating bin shifted spectra & fitting **************************************
        // **********************************************************************************************************************

        // Cloning spectra
        graphCombEtaInvYieldTotUnshi[cent]          = (TGraphAsymmErrors*)graphCombEtaInvYieldTot[cent]->Clone("EtaUnshifted");
        graphCombEtaInvYieldStatUnshi[cent]         = (TGraphAsymmErrors*)graphCombEtaInvYieldStat[cent]->Clone("EtaUnshiftedStat");
        graphCombEtaInvYieldSysUnshi[cent]          = (TGraphAsymmErrors*)graphCombEtaInvYieldSys[cent]->Clone("EtaUnshiftedSys");

        for (Int_t meth = 0; meth< 11; meth++){
            if (statErrorGraphCollectionEta[cent][meth]){
                graphIndEtaInvYieldStatUnshi[cent][meth]                 = (TGraphAsymmErrors*)statErrorGraphCollectionEta[cent][meth]->Clone(Form("EtaUnshiftedStat%s",nameMeasGlobalLabel[meth].Data()));
                graphIndEtaInvYieldStat[cent][meth]                      = (TGraphAsymmErrors*)statErrorGraphCollectionEta[cent][meth]->Clone(Form("EtaStat%s",nameMeasGlobalLabel[meth].Data()));
                graphIndEtaInvYieldStat_yShifted[cent][meth]             = (TGraphAsymmErrors*)statErrorGraphCollectionEta[cent][meth]->Clone(Form("EtaYShiftedStat%s",nameMeasGlobalLabel[meth].Data()));
            }
            if (sysErrorCollectionEta[cent][meth]){
                graphIndEtaInvYieldSysUnshi[cent][meth]                  = (TGraphAsymmErrors*)sysErrorCollectionEta[cent][meth]->Clone(Form("EtaUnshiftedSys%s",nameMeasGlobalLabel[meth].Data()));
                graphIndEtaInvYieldSys[cent][meth]                       = (TGraphAsymmErrors*)sysErrorCollectionEta[cent][meth]->Clone(Form("EtaSys%s",nameMeasGlobalLabel[meth].Data()));
                graphIndEtaInvYieldSys_yShifted[cent][meth]              = (TGraphAsymmErrors*)sysErrorCollectionEta[cent][meth]->Clone(Form("EtaYShiftedSys%s",nameMeasGlobalLabel[meth].Data()));
            }
        }

        // fitting spectrum with intial parameters
        // Two component model fit from Bylinkin
        fitTCMDecomposedLEta[cent]           = FitObject("tcmlow",Form("twoCompModelEta_DecL_%s", centArrayOutput[cent].Data()), "Eta", NULL, minPtEta[cent], 2.);
        fitTCMDecomposedHEta[cent]           = FitObject("tcmhigh",Form("twoCompModelEta_DecH_%s", centArrayOutput[cent].Data()), "Eta", NULL, 4, 50.);
        fitTCMDecomposedLEta[cent]->SetParameters(graphCombEtaInvYieldTot[cent]->GetY()[2],minPtEta[cent]);
        graphCombEtaInvYieldStat[cent]->Fit(fitTCMDecomposedLEta[cent],"QNRMEX0+","",minPtEta[cent],0.8);
        graphCombEtaInvYieldStat[cent]->Fit(fitTCMDecomposedHEta[cent],"QNRMEX0+","",3,20);

        cout << WriteParameterToFile(fitTCMDecomposedLEta[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedLEta[cent])<< endl;
        cout << WriteParameterToFile(fitTCMDecomposedHEta[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedHEta[cent])<< endl;

        fitTCMInvYieldEta[cent]              = FitObject("tcm",Form("fitTCMInvYieldEta_%s", centArrayOutput[cent].Data()),"Eta",graphCombEtaInvYieldStat[cent], minPtEta[cent], 20. , paramTCMEtaNew,
                                                         "QNRMEX0+", "", kFALSE);
        fitTCMDecomposedLEta[cent]->SetParameter(0, fitTCMInvYieldEta[cent]->GetParameter(0));
        fitTCMDecomposedLEta[cent]->SetParameter(1, fitTCMInvYieldEta[cent]->GetParameter(1));
        fitTCMDecomposedHEta[cent]->SetParameter(0, fitTCMInvYieldEta[cent]->GetParameter(2));
        fitTCMDecomposedHEta[cent]->SetParameter(1, fitTCMInvYieldEta[cent]->GetParameter(3));
        fitTCMDecomposedHEta[cent]->SetParameter(2, fitTCMInvYieldEta[cent]->GetParameter(4));

        // Tsallis fit
        fitInvYieldEta[cent]                 = FitObject("l","fitInvYieldEta","Eta",histoEtaInvYieldStat[cent][2],0.3,20.,paramGraphEta,"QNRMEX0+");
        fitInvYieldEtaGraph[cent]            = (TF1*)fitInvYieldEta[cent]->Clone("fitInvYieldEtaGraph");

        // *************************************************************************************************************
        // Shift graphs in X direction if desired
        // *************************************************************************************************************
        if(bWCorrection.Contains("X")){
            TF1* fitShiftingEta            = FitObject("tmpt","ShiftingEta","Eta");
            fitShiftingEta->SetParameters(fitInvYieldEta[cent]->GetParameter(0),fitInvYieldEta[cent]->GetParameter(1), fitInvYieldEta[cent]->GetParameter(2));
    //         TF1* fitShiftingEta                 = FitObject("tcmpt","ShiftingEta","Eta");
    //         fitShiftingEta->SetParameters(fitTCMInvYieldEta[cent]->GetParameter(0),fitTCMInvYieldEta[cent]->GetParameter(1), fitTCMInvYieldEta[cent]->GetParameter(2), fitTCMInvYieldEta[cent]->GetParameter(3),fitTCMInvYieldEta[cent]->GetParameter(4));

            TGraphAsymmErrors* graphCombEtaInvYieldTotNoShift = (TGraphAsymmErrors*) graphCombEtaInvYieldTot[cent]->Clone("Eta_NoShift");

            graphCombEtaInvYieldTot[cent]          = ApplyXshift(graphCombEtaInvYieldTot[cent], fitShiftingEta);
            cout << "comb" << endl;
            graphCombEtaInvYieldStat[cent]->Print();
            graphCombEtaInvYieldStat[cent]         = ApplyXshiftIndividualSpectra( graphCombEtaInvYieldTot[cent],
                                                                                graphCombEtaInvYieldStat[cent],
                                                                                fitShiftingEta,
                                                                                0, graphCombEtaInvYieldStat[cent]->GetN());
            graphCombEtaInvYieldSys[cent]          = ApplyXshiftIndividualSpectra( graphCombEtaInvYieldTot[cent],
                                                                                graphCombEtaInvYieldSys[cent],
                                                                                fitShiftingEta,
                                                                                0, graphCombEtaInvYieldSys[cent]->GetN());
            for (Int_t meth = 0; meth< 11; meth++){
                if (graphIndEtaInvYieldStat[cent][meth]){
                    cout << "shiting stat err of " << nameMeasGlobalLabel[meth].Data();
                    graphIndEtaInvYieldStat[cent][meth]  = ApplyXshiftIndividualSpectra( graphCombEtaInvYieldTot[cent],
                                                                                graphIndEtaInvYieldStat[cent][meth],
                                                                                fitShiftingEta,
                                                                                offSetEtaShifting[cent][meth], nComBinsEtaShifting[cent][meth]);

                }
                if (graphIndEtaInvYieldSys[cent][meth]){
                    cout << "shiting sys err of " << nameMeasGlobalLabel[meth].Data();
                    graphIndEtaInvYieldSys[cent][meth]   = ApplyXshiftIndividualSpectra( graphCombEtaInvYieldTot[cent],
                                                                                graphIndEtaInvYieldSys[cent][meth],
                                                                                fitShiftingEta,
                                                                                offSetEtaShifting[cent][meth], nComBinsEtaShifting[cent][meth]);
                }
            }

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

                Int_t numberPoints      = graphCombEtaInvYieldTotNoShift->GetN();
                Double_t *xPoint        = graphCombEtaInvYieldTotNoShift->GetX();
                Double_t* xvalueErrUp   = graphCombEtaInvYieldTotNoShift->GetEXhigh();
                Double_t* xvalueErrLow  = graphCombEtaInvYieldTotNoShift->GetEXlow();
                Double_t *xPointShift= graphCombEtaInvYieldTot[cent]->GetX();
                for (Int_t i=0; i<numberPoints; i++) {
                    graphCombEtaInvYieldTotNoShift->SetPoint(i,xPoint[i],xPointShift[i]/xPoint[i]);
                    graphCombEtaInvYieldTotNoShift->SetPointError(i,xvalueErrLow[i],xvalueErrUp[i],0,0);
                }
                DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldTotNoShift, markerStyleComb, markerSizeComb, colorComb , colorComb);
                graphCombEtaInvYieldTotNoShift->Draw("p same");

                TLatex *labelRatioToFitBinShift   = new TLatex(0.94, 0.91, Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
                SetStyleTLatex( labelRatioToFitBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitBinShift->Draw();
                TLatex *labelRatioToFitALICEBinShift    = new TLatex(0.94, 0.86, textALICE.Data());
                SetStyleTLatex( labelRatioToFitALICEBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitALICEBinShift->Draw();
                TLatex *labelRatioToFitEtaBinShift      = new TLatex(0.94, 0.807, "#eta #rightarrow #gamma#gamma");
                SetStyleTLatex( labelRatioToFitEtaBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitEtaBinShift->Draw();

            canvasShift->Update();
            canvasShift->SaveAs(Form("%s/BinShiftCorrection_Eta_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));

            // *************************************************************************************************************
            // Plot control graphs
            // *************************************************************************************************************

            TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
            DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.09);
            canvasDummy2->SetLogy();
            canvasDummy2->SetLogx();
            TH2F * histo2DDummy2;
            histo2DDummy2               = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.33,31.,1000,1e-9,10e1);
            SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
            histo2DDummy2->DrawCopy();

                for (Int_t meth = 0; meth< 11; meth++){
                    if (graphIndEtaInvYieldSys[cent][meth]){
                        DrawGammaSetMarkerTGraphAsym(graphIndEtaInvYieldSys[cent][meth], markerStyleDet[meth] ,markerSizeDet[meth]/2, colorDet[meth], colorDet[meth], widthLinesBoxes, kTRUE);
                        graphIndEtaInvYieldSys[cent][meth]->Draw("pEsame");
                    }
                }

                DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStatUnshi[cent], 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
                graphCombEtaInvYieldStatUnshi[cent]->Draw("pEsame");
                DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat[cent], 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
                graphCombEtaInvYieldStat[cent]->Draw("pEsame");

                fitTCMInvYieldEta[cent]->SetLineColor(kBlue+2);
                fitTCMInvYieldEta[cent]->Draw("same");

            canvasDummy2->Update();
            canvasDummy2->Print(Form("%s/ComparisonShiftedEta_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete canvasDummy2;
        }

        // *************************************************************************************************************
        // redo fitting after binshifts
        // *************************************************************************************************************
        // Tsallis function
        cout << WriteParameterToFile(fitInvYieldEta[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitInvYieldEta[cent])<< endl;
        //Two component model from Bylinkin
        fitTCMInvYieldEta[cent]        = FitObject("tcm","fitTCMInvYieldEtapPb5TeVCent","Eta",graphCombEtaInvYieldTot[cent],0.3,20. ,paramTCMEtaNew,"QNRMEX0+","", kFALSE);
        cout << WriteParameterToFile(fitTCMInvYieldEta[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMInvYieldEta[cent])<< endl;

        fitPowInvYieldEta[cent]        = FitObject("m","fitPowInvYieldEtapPb5TeVCent","Eta",graphCombEtaInvYieldTot[cent],5,40. ,NULL,"QNRMEX0+","", kFALSE);
        cout << WriteParameterToFile(fitPowInvYieldEta[cent])<< endl;

        // *************************************************************************************************************
        // Shift graphs in Y direction as well if desired
        // *************************************************************************************************************

        if(bWCorrection.Contains("Y") ){
            graphCombEtaInvYieldTot_yShifted[cent]        = (TGraphAsymmErrors*)graphCombEtaInvYieldTotUnshi[cent]->Clone("EtaYShiftedCombTot");
            graphCombEtaInvYieldTot_yShifted[cent]        =  ApplyYshiftIndividualSpectra( graphCombEtaInvYieldTot_yShifted[cent], fitInvYieldEta[cent]);
            graphCombEtaInvYieldStat_yShifted[cent]       = (TGraphAsymmErrors*)graphCombEtaInvYieldStatUnshi[cent]->Clone("EtaYShiftedCombStat");
            graphCombEtaInvYieldStat_yShifted[cent]       =  ApplyYshiftIndividualSpectra( graphCombEtaInvYieldStat_yShifted[cent], fitInvYieldEta[cent]);
            graphCombEtaInvYieldSys_yShifted[cent]        = (TGraphAsymmErrors*)graphCombEtaInvYieldSysUnshi[cent]->Clone("EtaYShiftedCombSys");
            graphCombEtaInvYieldSys_yShifted[cent]        =  ApplyYshiftIndividualSpectra( graphCombEtaInvYieldSys_yShifted[cent], fitInvYieldEta[cent]);

            for (Int_t meth = 0; meth< 11; meth++){
                if (graphIndEtaInvYieldStat_yShifted[cent][meth]){
                    graphIndEtaInvYieldStat_yShifted[cent][meth] = ApplyYshiftIndividualSpectra( graphIndEtaInvYieldStat_yShifted[cent][meth], fitInvYieldEta[cent]);
                }
                if (graphIndEtaInvYieldSys_yShifted[cent][meth]){
                    graphIndEtaInvYieldSys_yShifted[cent][meth]  = ApplyYshiftIndividualSpectra( graphIndEtaInvYieldSys_yShifted[cent][meth], fitInvYieldEta[cent]);
                }
            }
        }

        // *************************************************************************************************************
        // Calculate ratios to combined fit
        // *************************************************************************************************************
    //     TH1D* histoRatioEtaDPMJetToFit                      = (TH1D*) histoDPMJetEta->Clone("histoRatioEtaDPMJetToFit");
    //     histoRatioEtaDPMJetToFit                            = CalculateHistoRatioToFit (histoRatioEtaDPMJetToFit, fitTCMInvYieldEta[cent]);
    //     histoRatioEtaDPMJetToFit->GetXaxis()->SetRangeUser(0.3, 20);
    //     TH1D* histoRatioEtaHIJINGToFit                      = (TH1D*) histoHIJINGEta->Clone("histoRatioEtaHIJINGToFit");
    //     histoRatioEtaHIJINGToFit                            = CalculateHistoRatioToFit (histoRatioEtaHIJINGToFit, fitTCMInvYieldEta[cent]);
    //     histoRatioEtaHIJINGToFit->GetXaxis()->SetRangeUser(0.3, 20);
    //     TH1D* histoRatioEtaEPOS3ToFit                       = (TH1D*) histoEPOS3Eta->Clone("histoRatioEtaEPOS3ToFit");
    //     histoRatioEtaEPOS3ToFit                             = CalculateHistoRatioToFit (histoRatioEtaEPOS3ToFit, fitTCMInvYieldEta[cent]);
    //     histoRatioEtaEPOS3ToFit->GetXaxis()->SetRangeUser(0.3, 20);
    //     TGraphErrors* graphRatioEtaEPOS3ToFit               = new TGraphErrors(histoRatioEtaEPOS3ToFit);
    //     while(graphRatioEtaEPOS3ToFit->GetX()[0] < 0.3)
    //         graphRatioEtaEPOS3ToFit->RemovePoint(0);
    //     while(graphRatioEtaEPOS3ToFit->GetX()[graphRatioEtaEPOS3ToFit->GetN()-1] > 20)
    //         graphRatioEtaEPOS3ToFit->RemovePoint(graphRatioEtaEPOS3ToFit->GetN()-1);
    //     TGraphErrors* graphRatioEtaMcGillToFit              = (TGraphErrors*)graphMcGillEta->Clone("graphRatioEtaMcGillToFit");
    //     graphRatioEtaMcGillToFit                            = CalculateGraphErrRatioToFit(graphRatioEtaMcGillToFit, fitTCMInvYieldEta[cent]);
    //     while(graphRatioEtaMcGillToFit->GetX()[0] < 0.3)
    //         graphRatioEtaMcGillToFit->RemovePoint(0);


    //     TGraph* graphRatioEtaCGC                            = (TGraph*)graphEtaCGC5TeV->Clone("graphRatioEtaCGC");
    //     graphRatioEtaCGC                                    = CalculateGraphRatioToFit(graphRatioEtaCGC, fitTCMInvYieldEta[cent]);
    //     TGraphAsymmErrors* graphRatioEtaNLODSS14            = (TGraphAsymmErrors*)graphNLODSS14Eta->Clone("graphRatioEtaNLODSS14ToFit");
    //     graphRatioEtaNLODSS14                               = CalculateGraphErrRatioToFit(graphRatioEtaNLODSS14, fitTCMInvYieldEta[cent]);
    //     TGraph* graphRatioEtaNLODSS14Center                 = (TGraph*)graphNLODSS14EtaCenter->Clone("graphRatioEtaNLODSS14CenterToFit");
    //     graphRatioEtaNLODSS14Center                         = CalculateGraphRatioToFit(graphRatioEtaNLODSS14Center, fitTCMInvYieldEta[cent]);
    //     TGraphAsymmErrors* graphRatioEtaNLODSS14nPDF        = (TGraphAsymmErrors*)graphNLODSS14nPDFEta->Clone("graphRatioEtaNLODSS14nPDFToFit");
    //     graphRatioEtaNLODSS14nPDF                           = CalculateGraphErrRatioToFit(graphRatioEtaNLODSS14nPDF, fitTCMInvYieldEta[cent]);
    //     TGraph* graphRatioEtaNLODSS14nPDFCenter             = (TGraph*)graphNLODSS14nPDFEtaCenter->Clone("graphRatioEtaNLODSS14nPDFCenterToFit");
    //     graphRatioEtaNLODSS14nPDFCenter                     = CalculateGraphRatioToFit(graphRatioEtaNLODSS14nPDFCenter, fitTCMInvYieldEta[cent]);
    //     TGraphAsymmErrors* graphRatioEtaNLODSS14nPDFEPPS16  = (TGraphAsymmErrors*)graphNLODSS14nPDFEPPS16Eta->Clone("graphRatioEtaNLODSS14nPDFEPPS16ToFit");
    //     graphRatioEtaNLODSS14nPDFEPPS16                     = CalculateGraphErrRatioToFit(graphRatioEtaNLODSS14nPDFEPPS16, fitTCMInvYieldEta[cent]);
    //     TGraph* graphRatioEtaNLODSS14nPDFEPPS16Center       = (TGraph*)graphNLODSS14nPDFEPPS16EtaCenter->Clone("graphRatioEtaNLODSS14nPDFEPPS16CenterToFit");
    //     graphRatioEtaNLODSS14nPDFEPPS16Center               = CalculateGraphRatioToFit(graphRatioEtaNLODSS14nPDFEPPS16Center, fitTCMInvYieldEta[cent]);
    //
    //
        graphRatioEtaCombCombFitTot[cent]                       = (TGraphAsymmErrors*)graphCombEtaInvYieldTot[cent]->Clone();
        graphRatioEtaCombCombFitTot[cent]                       = CalculateGraphErrRatioToFit(graphRatioEtaCombCombFitTot[cent], fitTCMInvYieldEta[cent]);
        graphRatioEtaCombCombFitStat[cent]                      = (TGraphAsymmErrors*)graphCombEtaInvYieldStat[cent]->Clone();
        graphRatioEtaCombCombFitStat[cent]                      = CalculateGraphErrRatioToFit(graphRatioEtaCombCombFitStat[cent], fitTCMInvYieldEta[cent]);
        graphRatioEtaCombCombFitSys[cent]                       = (TGraphAsymmErrors*)graphCombEtaInvYieldSys[cent]->Clone();
        graphRatioEtaCombCombFitSys[cent]                       = CalculateGraphErrRatioToFit(graphRatioEtaCombCombFitSys[cent], fitTCMInvYieldEta[cent]);

        for (Int_t meth = 0; meth< 11; meth++){
            if (graphIndEtaInvYieldStat[cent][meth]){
                graphRatioEtaIndCombFitStat[cent][meth]              = (TGraphAsymmErrors*)graphIndEtaInvYieldStat[cent][meth]->Clone(Form("RatioEta%sToCombFitStat", nameMeasGlobalLabel[meth].Data()));
                graphRatioEtaIndCombFitStat[cent][meth]              = CalculateGraphErrRatioToFit(graphRatioEtaIndCombFitStat[cent][meth], fitTCMInvYieldEta[cent]);
            }
            if (graphIndEtaInvYieldSys[cent][meth]){
                graphRatioEtaIndCombFitSys[cent][meth]              = (TGraphAsymmErrors*)graphIndEtaInvYieldSys[cent][meth]->Clone(Form("RatioEta%sToCombFitSyst", nameMeasGlobalLabel[meth].Data()));
                graphRatioEtaIndCombFitSys[cent][meth]              = CalculateGraphErrRatioToFit(graphRatioEtaIndCombFitSys[cent][meth], fitTCMInvYieldEta[cent]);
            }
        }

        // **********************************************************************************************************************
        // ******************************************* Plot Ratio of Comb to Fit ************************************************
        // **********************************************************************************************************************
        textSizeLabelsPixel                 = 54;
        TCanvas* canvasRatioToCombFit       = new TCanvas("canvasRatioToCombFit","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRatioToCombFit, 0.08, 0.01, 0.01, 0.125);
        canvasRatioToCombFit->SetLogx();

            Double_t textsizeLabelspPb      = 0;
            if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
                textsizeLabelspPb           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
            } else {
                textsizeLabelspPb           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
            }
            cout << textsizeLabelspPb << endl;

        TH2F * histo2DEtaRatioToCombFit;
        histo2DEtaRatioToCombFit               = new TH2F("histo2DEtaRatioToCombFit","histo2DEtaRatioToCombFit",1000,minPtEtaPlotting, maxPtEtaPlotting,1000,0.2,4.    );
        SetStyleHistoTH2ForGraphs(histo2DEtaRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelspPb, textsizeLabelspPb,
                                  0.85*textsizeLabelspPb,textsizeLabelspPb, 0.9, 0.65, 510, 505);
        histo2DEtaRatioToCombFit->GetXaxis()->SetMoreLogLabels();
        histo2DEtaRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
        histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.2,1.82);
        histo2DEtaRatioToCombFit->Draw("copy");

            ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombCombFitStat[cent]);

            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys[cent], markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
            graphRatioEtaCombCombFitSys[cent]->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStat[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphRatioEtaCombCombFitStat[cent]->Draw("p,same,z");

            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1., 1.,0.1, kGray+2);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1.1, 1.1,0.1, kGray, 7);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.9, 0.9,0.1, kGray, 7);

            TLatex *labelRatioToFitEnergy   = new TLatex(0.95, 0.92, Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
            SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
            labelRatioToFitEnergy->Draw();
            TLatex *labelRatioToFitALICE    = new TLatex(0.95, 0.86, textALICE.Data());
            SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
            labelRatioToFitALICE->Draw();
            TLatex *labelRatioToFitEta      = new TLatex(0.12, 0.92, "#eta #rightarrow #gamma#gamma");
            SetStyleTLatex( labelRatioToFitEta, textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
            labelRatioToFitEta->Draw();

        canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfCombToCombFit_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));

        // **********************************************************************************************************************
        // *******************************************Plot Ratio of Individual meas to Fit ******************************************
        // **********************************************************************************************************************

        canvasRatioToCombFit->cd();
        histo2DEtaRatioToCombFit->Draw("copy");

            for (Int_t meth = 10; meth > -1 ; meth--){
                if (graphRatioEtaIndCombFitSys[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRatioEtaIndCombFitSys[cent][meth], markerStyleDet[meth] ,markerSizeDet[meth]*0.5, colorDet[meth], colorDet[meth], widthLinesBoxes, kTRUE);
                    graphRatioEtaIndCombFitSys[cent][meth]->Draw("E2same");
                }
                if (graphRatioEtaIndCombFitStat[cent][meth]){
                    ProduceGraphAsymmWithoutXErrors(graphRatioEtaIndCombFitStat[cent][meth]);
                    DrawGammaSetMarkerTGraphAsym(graphRatioEtaIndCombFitStat[cent][meth], markerStyleDet[meth] ,markerSizeDet[meth]*0.5, colorDet[meth], colorDet[meth]);
                    graphRatioEtaIndCombFitStat[cent][meth]->Draw("p,same,z");
                }
            }
            if (graphRatioEtaIndCombFitStat[cent][4]){
                graphRatioEtaIndCombFitStat[cent][4]->Draw("p,same,z");
            }

            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1., 1.,0.5, kGray+2);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1.1, 1.1,0.5, kGray, 7);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.9, 0.9,0.5, kGray, 7);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1.2, 1.2,0.5, kGray, 9);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.8, 0.8,0.5, kGray, 9);

            labelRatioToFitEnergy->Draw();
            labelRatioToFitALICE->Draw();
            labelRatioToFitEta->Draw();
            histo2DEtaRatioToCombFit->Draw("same,axis");

            //****************************** Definition of the Legend ******************************************
            //**************** Row def ************************
            Double_t rowsLegendOnlyEtaRatio[4]          = {0.30, 0.25, 0.20, 0.15};
            Double_t rowsLegendOnlyEtaRatioAbs[4]       = {0.92, 0.47, 0.38, 0.29 };
            Double_t columnsLegendOnlyEtaRatio[6]       = {0.14, 0.26, 0.35, 0.48, 0.7, 0.8};
            Double_t columnsLegendOnlyEtaRatioAbs[6]    = {0.215, 1.19, 1.60, 2, 9.4, 13.1};
            Double_t lengthBox                          = 0.2;
            Double_t heightBox                          = 0.06/2;
            //****************** first Column **************************************************
            TLatex *textPCMOnlyRatioEta                 = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[1],nameMeasGlobalLabel[0]);
            SetStyleTLatex( textPCMOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
            textPCMOnlyRatioEta->Draw();
            TLatex *textPHOSOnlyRatioEta                = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[2],nameMeasGlobalLabel[1]);
            SetStyleTLatex( textPHOSOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
            textPHOSOnlyRatioEta->Draw();
            TLatex *textEMCALOnlyRatioEta               = new TLatex(columnsLegendOnlyEtaRatio[0],rowsLegendOnlyEtaRatio[3],nameMeasGlobalLabel[2]);
            SetStyleTLatex( textEMCALOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
            textEMCALOnlyRatioEta->Draw();
            TLatex *textPCMEMCALOnlyRatioEta            = new TLatex(columnsLegendOnlyEtaRatio[3],rowsLegendOnlyEtaRatio[1],nameMeasGlobalLabel[4]);
            SetStyleTLatex( textPCMEMCALOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
            textPCMEMCALOnlyRatioEta->Draw();
            TLatex *textPCMPHOSOnlyRatioEta            = new TLatex(columnsLegendOnlyEtaRatio[3],rowsLegendOnlyEtaRatio[2],nameMeasGlobalLabel[3]);
            SetStyleTLatex( textPCMPHOSOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
            textPCMPHOSOnlyRatioEta->Draw();

            //****************** second Column *************************************************
            TLatex *textStatOnlyRatioEta                = new TLatex(columnsLegendOnlyEtaRatio[1],rowsLegendOnlyEtaRatio[0] ,"stat");
            SetStyleTLatex( textStatOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
            textStatOnlyRatioEta->Draw();
            TLatex *textSysOnlyRatioEta                 = new TLatex(columnsLegendOnlyEtaRatio[2] ,rowsLegendOnlyEtaRatio[0],"syst");
            SetStyleTLatex( textSysOnlyRatioEta, textSizeLabelsPixel,4, 1, 43);
            textSysOnlyRatioEta->Draw();
            TLatex *textStatOnlyRatioEta2               = new TLatex(columnsLegendOnlyEtaRatio[4],rowsLegendOnlyEtaRatio[0] ,"stat");
            SetStyleTLatex( textStatOnlyRatioEta2, textSizeLabelsPixel,4, 1, 43);
            textStatOnlyRatioEta2->Draw();
            TLatex *textSysOnlyRatioEta2                = new TLatex(columnsLegendOnlyEtaRatio[5] ,rowsLegendOnlyEtaRatio[0],"syst");
            SetStyleTLatex( textSysOnlyRatioEta2, textSizeLabelsPixel,4, 1, 43);
            textSysOnlyRatioEta2->Draw();

            if (graphRatioEtaIndCombFitSys[cent][0]){
                TMarker* markerPCMEtaOnlyRatioEta           = CreateMarkerFromGraph(graphRatioEtaIndCombFitSys[cent][0],columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[1],1);
                markerPCMEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[1]);
            }
            if (graphRatioEtaIndCombFitSys[cent][1]){
                TMarker* markerPHOSEtaOnlyRatioEta          = CreateMarkerFromGraph(graphRatioEtaIndCombFitSys[cent][1], columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[2],1);
                markerPHOSEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[2]);
            }
            if (graphRatioEtaIndCombFitSys[cent][2]){
                TMarker* markerEMCALEtaOnlyRatioEta         = CreateMarkerFromGraph(graphRatioEtaIndCombFitSys[cent][2], columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[3],1);
                markerEMCALEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[3]);
            }
            if (graphRatioEtaIndCombFitSys[cent][4]){
                TMarker* markerPCMEMCALEtaOnlyRatioEta      = CreateMarkerFromGraph(graphRatioEtaIndCombFitSys[cent][4], columnsLegendOnlyEtaRatio[3] ,rowsLegendOnlyEtaRatio[1],1);
                markerPCMEMCALEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[4] ,rowsLegendOnlyEtaRatioAbs[1]);
            }
            if (graphRatioEtaIndCombFitSys[cent][3]){
                TMarker* markerPCMPHOSEtaOnlyRatioEta      = CreateMarkerFromGraph(graphRatioEtaIndCombFitSys[cent][3], columnsLegendOnlyEtaRatio[3] ,rowsLegendOnlyEtaRatio[2],1);
                markerPCMPHOSEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[4] ,rowsLegendOnlyEtaRatioAbs[2]);
            }
            if (graphRatioEtaIndCombFitSys[cent][0]){
                TBox* boxPCMEtaOnlyRatioEta                 = CreateBoxFromGraph(graphRatioEtaIndCombFitSys[cent][0], columnsLegendOnlyEtaRatioAbs[2]-0.1*lengthBox , rowsLegendOnlyEtaRatioAbs[1]- heightBox,
                                                                                columnsLegendOnlyEtaRatioAbs[2]+ 2.3*lengthBox, rowsLegendOnlyEtaRatioAbs[1]+ heightBox);
                boxPCMEtaOnlyRatioEta->Draw("l");
            }
            if (graphRatioEtaIndCombFitSys[cent][1]){
                TBox* boxPHOSEtaOnlyRatioEta                = CreateBoxFromGraph(graphRatioEtaIndCombFitSys[cent][1], columnsLegendOnlyEtaRatioAbs[2]-0.1*lengthBox , rowsLegendOnlyEtaRatioAbs[2]- heightBox,
                                                                             columnsLegendOnlyEtaRatioAbs[2]+ 2.3*lengthBox, rowsLegendOnlyEtaRatioAbs[2]+ heightBox);
                boxPHOSEtaOnlyRatioEta->Draw("l");
            }
            if (graphRatioEtaIndCombFitSys[cent][2]){
                TBox* boxEMCALEtaOnlyRatioEta               = CreateBoxFromGraph(graphRatioEtaIndCombFitSys[cent][2], columnsLegendOnlyEtaRatioAbs[2]-0.1*lengthBox , rowsLegendOnlyEtaRatioAbs[3]- heightBox,
                                                                                columnsLegendOnlyEtaRatioAbs[2]+ 2.3*lengthBox, rowsLegendOnlyEtaRatioAbs[3]+ heightBox);
                boxEMCALEtaOnlyRatioEta->Draw("l");
            }
            if (graphRatioEtaIndCombFitSys[cent][4]){
                TBox* boxPCMEMCALEtaOnlyRatioEta            = CreateBoxFromGraph(graphRatioEtaIndCombFitSys[cent][4], columnsLegendOnlyEtaRatioAbs[5]-0.5*lengthBox , rowsLegendOnlyEtaRatioAbs[1]- heightBox,
                                                                                columnsLegendOnlyEtaRatioAbs[5]+ 19*lengthBox, rowsLegendOnlyEtaRatioAbs[1]+ heightBox);
                boxPCMEMCALEtaOnlyRatioEta->Draw("l");
            }
            if (graphRatioEtaIndCombFitSys[cent][3]){
                TBox* boxPCMPHOSEtaOnlyRatioEta             = CreateBoxFromGraph(graphRatioEtaIndCombFitSys[cent][3], columnsLegendOnlyEtaRatioAbs[5]-0.5*lengthBox , rowsLegendOnlyEtaRatioAbs[2]- heightBox,
                                                                                columnsLegendOnlyEtaRatioAbs[5]+ 19*lengthBox, rowsLegendOnlyEtaRatioAbs[2]+ heightBox);
                boxPCMPHOSEtaOnlyRatioEta->Draw("l");
            }
        canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfIndividualMeasToCombFit_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), runArray[cent].Data() ,suffix.Data()));
    }

    // **********************************************************************************************************************
    // ******************************************* Mass and width for pi0 ***************************************************
    // **********************************************************************************************************************

    Double_t arrayBoundariesX1_4[2];
    Double_t arrayBoundariesY1_4[3];
    Double_t relativeMarginsX[3];
    Double_t relativeMarginsY[3];
    textSizeLabelsPixel                 = 50;
    Double_t textSizeLabelsRel          = ((Double_t)textSizeLabelsPixel)/1250;
    ReturnCorrectValuesForCanvasScaling(1350,1250, 1, 2,0.09, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

    TCanvas* canvasMassWidthPi0         = new TCanvas("canvasMassWidthPi0","",0,0,1350,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasMassWidthPi0,  0.13, 0.02, 0.03, 0.06);

    TPad* padWidthPi0                   = new TPad("padWidthPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthPi0->Draw();

    TPad* padMassPi0                    = new TPad("padMassPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
    padMassPi0->Draw();

    TPad* padMassLegendPi0                 = new TPad("padMassLegendPi0", "", 0.13, 0.36, 0.82, 0.52,-1, -1, -2);
    DrawGammaPadSettings( padMassLegendPi0, 0., 0., 0., 0.);
    padMassLegendPi0->SetFillStyle(0);
    padMassLegendPi0->Draw();

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

    TH2F * histo2DAllPi0FWHM        = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, minPtPi0Plotting, maxPtPi0Plotting ,1000., -30, 40);
    SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                              0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
    histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,30.5);
    histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
    histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
    histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);

    Double_t textsizeLabelsMass     = 0;
    Double_t textsizeFacMass        = 0;
    if (padMassPi0->XtoPixel(padMassPi0->GetX2()) <padMassPi0->YtoPixel(padMassPi0->GetY1()) ){
        textsizeLabelsMass          = (Double_t)textSizeLabelsPixel/padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
        textsizeFacMass             = (Double_t)1./padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
    } else {
        textsizeLabelsMass          = (Double_t)textSizeLabelsPixel/padMassPi0->YtoPixel(padMassPi0->GetY1());
        textsizeFacMass             = (Double_t)1./padMassPi0->YtoPixel(padMassPi0->GetY1());
    }

    TH2F * histo2DAllPi0Mass        = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, minPtPi0Plotting, maxPtPi0Plotting, 1000., 120., 170);
    SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                              textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
    histo2DAllPi0Mass->GetYaxis()->SetRangeUser(125.5,162.8);
    histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
    histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
    histo2DAllPi0Mass->GetXaxis()->SetLabelOffset(-0.015);

    TLatex *labelLegendAMass        = new TLatex(0.13,0.06,"a)");
    SetStyleTLatex( labelLegendAMass, textSizeLabelsPixel,4, 1, 43);
    TLatex *labelMassPerf           = new TLatex(0.13,0.875,"ALICE performance");
    SetStyleTLatex( labelMassPerf, textSizeLabelsPixel,4, 1, 43);
    TLatex *labelMassEnergy         = new TLatex(0.13,0.775,collisionSystempPb.Data());
    SetStyleTLatex( labelMassEnergy, textSizeLabelsPixel,4, 1, 43);
    TLatex *labelMassPi0            = new TLatex(0.13,0.69,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelMassPi0, textSizeLabelsPixel,4, 1, 43);
    TLatex *labelLegendBMass        = new TLatex(0.13,0.22,"b)");
    SetStyleTLatex( labelLegendBMass, textSizeLabelsPixel,4, 1, 43);
    //********************************** Defintion of the Legend **************************************************
    Double_t columnsLegendMass[6]  = {0.,0.25,0.38, 0.5, 0.78, 0.9};
    Double_t rowsLegendMass[5]     = {0.8,0.6,0.4,0.2,0.01};
    Double_t offsetMarkerYMass     = 0.06;
    Double_t offsetMarkerXMass     = 0.05;
    Double_t scaleMarkerMass       = 1.2;

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCent[cent]) continue;
        padWidthPi0->cd();
        padWidthPi0->SetLogx();

            histo2DAllPi0FWHM->DrawCopy();

            for (Int_t meth = 10; meth > -1; meth--){
                if (graphPi0Width[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphPi0Width[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.55, colorDet[meth] , colorDet[meth]);
                    graphPi0Width[cent][meth]->Draw("p,same,z");
                }
                if (graphPi0WidthMC[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphPi0WidthMC[cent][meth], markerStyleDetMC[meth], markerSizeDetMC[meth]*0.55, colorDetMC[meth] , colorDetMC[meth]);
                    graphPi0WidthMC[cent][meth]->Draw("p,same,z");
                }
            }

            labelLegendAMass->Draw();
            labelMassPerf->Draw();
            labelMassEnergy->SetText(0.13,0.775,Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
            labelMassEnergy->Draw();
            labelMassPi0->Draw();

        padMassPi0->cd();
        padMassPi0->SetLogx();

            histo2DAllPi0Mass->DrawCopy();

            for (Int_t meth = 10; meth > -1; meth--){
                if (graphPi0Mass[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphPi0Mass[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.55, colorDet[meth] , colorDet[meth]);
                    graphPi0Mass[cent][meth]->Draw("p,same,z");
                }
                if (graphPi0MassMC[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphPi0MassMC[cent][meth], markerStyleDetMC[meth], markerSizeDetMC[meth]*0.55, colorDetMC[meth] , colorDetMC[meth]);
                    graphPi0MassMC[cent][meth]->Draw("p,same,z");
                }
            }

            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);
            labelLegendBMass->Draw();

            padMassLegendPi0->cd();
            padMassLegendPi0->Clear();
            //****************** first Column **************************************************
            TLatex *textMassPCM2            = new TLatex(columnsLegendMass[0],rowsLegendMass[1],nameMeasGlobalLabel[0]);
            SetStyleTLatex( textMassPCM2, textSizeLabelsPixel,4, 1, 43);
            if (graphPi0Mass[cent][0]) textMassPCM2->Draw();
            TLatex *textMassPCMEMCAL       = new TLatex(columnsLegendMass[0],rowsLegendMass[2],nameMeasGlobalLabel[4]);
            SetStyleTLatex( textMassPCMEMCAL, textSizeLabelsPixel,4, 1, 43);
            if (graphPi0Mass[cent][4]) textMassPCMEMCAL->Draw();
            TLatex *textMassEMCAL          = new TLatex(columnsLegendMass[0],rowsLegendMass[3],nameMeasGlobalLabel[2]);
            SetStyleTLatex( textMassEMCAL, textSizeLabelsPixel,4, 1, 43);
            if (graphPi0Mass[cent][2]) textMassEMCAL->Draw();
            TLatex *textMassPHOS           = new TLatex(columnsLegendMass[3],rowsLegendMass[1],nameMeasGlobalLabel[1]);
            SetStyleTLatex( textMassPHOS, textSizeLabelsPixel,4, 1, 43);
            if (graphPi0Mass[cent][1]) textMassPHOS->Draw();
            TLatex *textMassPCMPHOS             = new TLatex(columnsLegendMass[3],rowsLegendMass[2],nameMeasGlobalLabel[3]);
            SetStyleTLatex( textMassPCMPHOS, textSizeLabelsPixel,4, 1, 43);
            if (graphPi0Mass[cent][3]) textMassPCMPHOS->Draw();
            TLatex *textMassDalitz             = new TLatex(columnsLegendMass[3],rowsLegendMass[2],nameMeasGlobalLabel[5]);
            SetStyleTLatex( textMassDalitz, textSizeLabelsPixel,4, 1, 43);
            if (graphPi0Mass[cent][5]) textMassDalitz->Draw();

            //****************** second Column *************************************************
            TLatex *textMassData1           = new TLatex(columnsLegendMass[1],rowsLegendMass[0] ,"Data");
            SetStyleTLatex( textMassData1, textSizeLabelsPixel,4, 1, 43);
            if (graphPi0Mass[cent][2] || graphPi0Mass[cent][4] || graphPi0Mass[cent][0]) textMassData1->Draw();
            TLatex *textMassMC1             = new TLatex(columnsLegendMass[2] ,rowsLegendMass[0],"MC");
            SetStyleTLatex( textMassMC1, textSizeLabelsPixel,4, 1, 43);
            if (graphPi0Mass[cent][2] || graphPi0Mass[cent][4] || graphPi0Mass[cent][0]) textMassMC1->Draw();

            TLatex *textMassData2           = new TLatex(columnsLegendMass[4],rowsLegendMass[0] ,"Data");
            SetStyleTLatex( textMassData2, textSizeLabelsPixel,4, 1, 43);
            if (graphPi0Mass[cent][1] || graphPi0Mass[cent][3] ) textMassData2->Draw();
            TLatex *textMassMC2             = new TLatex(columnsLegendMass[5] ,rowsLegendMass[0],"MC");
            SetStyleTLatex( textMassMC2, textSizeLabelsPixel,4, 1, 43);
            if (graphPi0Mass[cent][1] || graphPi0Mass[cent][3] ) textMassMC2->Draw();

            TMarker* markerPi0Mass[11]      = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
            TMarker* markerPi0MassMC[11]    = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

            for (Int_t meth = 0; meth< 11; meth++){
                if (graphPi0Mass[cent][meth]) markerPi0Mass[meth]            = CreateMarkerFromGraph(graphPi0Mass[cent][meth], 0.1, 0.1, scaleMarkerMass);
                if (graphPi0MassMC[cent][meth]) markerPi0MassMC[meth]        = CreateMarkerFromGraph(graphPi0MassMC[cent][meth], 0.1, 0.1, scaleMarkerMass);
            }
            if (markerPi0Mass[0] ) markerPi0Mass[0]->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[1]+ offsetMarkerYMass);
            if (markerPi0Mass[4] ) markerPi0Mass[4]->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
            if (markerPi0Mass[2] ) markerPi0Mass[2]->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
            if (markerPi0Mass[1] ) markerPi0Mass[1]->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[1]+ offsetMarkerYMass);
            if (markerPi0Mass[3] ) markerPi0Mass[3]->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
            if (markerPi0Mass[5] ) markerPi0Mass[5]->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);

            if (markerPi0MassMC[0]) markerPi0MassMC[0]->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass-0.01 ,rowsLegendMass[1]+ offsetMarkerYMass);
            if (markerPi0MassMC[4]) markerPi0MassMC[4]->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass-0.01 ,rowsLegendMass[2]+ offsetMarkerYMass);
            if (markerPi0MassMC[2]) markerPi0MassMC[2]->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass-0.01 ,rowsLegendMass[3]+ offsetMarkerYMass);
            if (markerPi0MassMC[1]) markerPi0MassMC[1]->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass-0.01 ,rowsLegendMass[1]+ offsetMarkerYMass);
            if (markerPi0MassMC[3]) markerPi0MassMC[3]->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass-0.01 ,rowsLegendMass[2]+ offsetMarkerYMass);
            if (markerPi0MassMC[5]) markerPi0MassMC[5]->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass-0.01 ,rowsLegendMass[3]+ offsetMarkerYMass);

        canvasMassWidthPi0->Update();
        canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth_%s%s.%s",outputDir.Data(),centArrayOutput[cent].Data(), runArray[cent].Data(), suffix.Data()));

        for (Int_t meth = 0; meth< 11; meth++){
            delete markerPi0Mass[meth];
            delete markerPi0MassMC[meth];
        }

        delete textMassPCM2;
        delete textMassPCMEMCAL;
        delete textMassEMCAL;
        delete textMassPHOS;
        delete textMassPCMPHOS;
        delete textMassData1;
        delete textMassMC1;
        delete textMassData2;
        delete textMassMC2;
    }

    Float_t xRangeMass[2][11] = {   {0.23, 0.73, 0.99, 0.23, 0.63, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23},
                                    {20., 20., 30., 20., 30., 20., 30., 30., 30., 30. ,30.} };
    Float_t yRangeMass[2][11] = {   {129.5, 129.5, 125.5, 129.5, 125.5, 125.5, 125.5, 125.5, 125.5, 125.5, 125.5},
                                    {144.5, 151.5, 162.8, 148.5, 162.8, 162.8, 162.8, 162.8, 162.8, 162.8 ,162.8} };

    DrawGammaCanvasSettings( canvasMassWidthPi0,  0.09, 0.02, 0.01, 0.082);
    TH2F * histo2DCentPi0Mass        = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",2000, minPtPi0Plotting, maxPtPi0Plotting, 1000., 120., 170);
    SetStyleHistoTH2ForGraphs(histo2DCentPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel,
                              textSizeLabelsRel, 0.9, 1.1, 512, 505);
    histo2DCentPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DCentPi0Mass->GetYaxis()->SetNdivisions(505);
    histo2DCentPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
    histo2DCentPi0Mass->GetXaxis()->SetTickLength(0.025);
    histo2DCentPi0Mass->GetXaxis()->SetLabelOffset(-0.008);
    TLatex *labelPerfMass2           = new TLatex(0.15,0.92,"ALICE performance");
    SetStyleTLatex( labelPerfMass2, textSizeLabelsRel,4);
    TLatex *labelEnergyMass2         = new TLatex(0.15,0.87,collisionSystempPb.Data());
    SetStyleTLatex( labelEnergyMass2, textSizeLabelsRel,4);
    TLatex *labelDetSysPi0Mass      = new TLatex(0.15,0.82,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysPi0Mass, textSizeLabelsRel,4);


    for (Int_t meth = 0; meth < 11; meth++){
        if (!fileNamesMethod[meth].CompareTo("")) continue;

        canvasMassWidthPi0->cd();
        canvasMassWidthPi0->SetLogx(1);
        histo2DCentPi0Mass->GetXaxis()->SetRangeUser(xRangeMass[0][meth],xRangeMass[1][meth]);
        histo2DCentPi0Mass->GetYaxis()->SetRangeUser(yRangeMass[0][meth],yRangeMass[1][meth]);
        histo2DCentPi0Mass->DrawCopy();

        TLegend* legendMass           = GetAndSetLegend2(0.65, 0.95-(5*textSizeLabelsRel), 0.93, 0.95,textSizeLabelsPixel);
        for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
            if (!enableCent[cent]) continue;
            if (graphPi0Mass[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphPi0Mass[cent][meth], markerStyleCent[cent], markerSizeCent[cent]*0.55, colorCent[cent] , colorCent[cent]);
                graphPi0Mass[cent][meth]->Draw("p,same,z");
                legendMass->AddEntry(graphPi0Mass[cent][meth],centArray[cent],"p");
            }
        }
        legendMass->Draw();
        labelPerfMass2->Draw();
        labelEnergyMass2->SetText(0.15,0.87,Form("%s", collisionSystempPb.Data()));
        labelEnergyMass2->Draw();
        labelDetSysPi0Mass->SetText(0.15,0.82,Form("#pi^{0} #rightarrow #gamma#gamma, %s", nameMeasGlobalLabel[meth].Data()));
        labelDetSysPi0Mass->Draw();

        canvasMassWidthPi0->Update();
        canvasMassWidthPi0->Print(Form("%s/Pi0_MassPositionData_%s.%s",outputDirSupport.Data(), nameMeasGlobalLabel[meth].Data(), suffix.Data()));

        for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
            if (!enableCent[cent]) continue;
            if (graphPi0MassMC[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphPi0MassMC[cent][meth], markerStyleCentMC[cent], markerSizeCent[cent]*0.55, colorCentMC[cent] , colorCentMC[cent]);
                graphPi0MassMC[cent][meth]->Draw("p,same,z");
            }
        }

        canvasMassWidthPi0->Update();
        canvasMassWidthPi0->Print(Form("%s/Pi0_MassPositionDataAndMC_%s.%s",outputDirSupport.Data(), nameMeasGlobalLabel[meth].Data(), suffix.Data()));
    }



    // **********************************************************************************************************************
    // ******************************************* Mass and width for pi0 ***************************************************
    // **********************************************************************************************************************

    TCanvas* canvasMassWidthEta         = new TCanvas("canvasMassWidthEta","",0,0,1350,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasMassWidthEta,  0.13, 0.02, 0.03, 0.06);

    TPad* padWidthEta                   = new TPad("padWidthEta", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthEta, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthEta->Draw();

    TPad* padMassEta                    = new TPad("padMassEta", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassEta, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
    padMassEta->Draw();

    TPad* padMassLegendEta                 = new TPad("padMassLegendEta", "", 0.13, 0.36, 0.82, 0.52,-1, -1, -2);
    DrawGammaPadSettings( padMassLegendEta, 0., 0., 0., 0.);
    padMassLegendEta->SetFillStyle(0);
    padMassLegendEta->Draw();

    TH2F * histo2DAllEtaFWHM        = new TH2F("histo2DAllEtaFWHM","histo2DAllEtaFWHM", 100, minPtEtaPlotting, maxPtEtaPlotting ,1000., -5, 92);
    SetStyleHistoTH2ForGraphs(histo2DAllEtaFWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                              0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
    histo2DAllEtaFWHM->GetYaxis()->SetRangeUser(-5.,92);
    histo2DAllEtaFWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllEtaFWHM->GetYaxis()->SetNdivisions(505);
    histo2DAllEtaFWHM->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllEtaFWHM->GetXaxis()->SetTickLength(0.05);
    histo2DAllEtaFWHM->GetYaxis()->SetTickLength(0.026);

    TH2F * histo2DAllEtaMass        = new TH2F("histo2DAllEtaMass","histo2DAllEtaMass",100, minPtEtaPlotting, maxPtEtaPlotting, 1000., 505., 590);
    SetStyleHistoTH2ForGraphs(histo2DAllEtaMass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                              textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
    histo2DAllEtaMass->GetYaxis()->SetRangeUser(125.5,162.8);
    histo2DAllEtaMass->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllEtaMass->GetYaxis()->SetNdivisions(505);
    histo2DAllEtaMass->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllEtaMass->GetXaxis()->SetTickLength(0.05);
    histo2DAllEtaMass->GetXaxis()->SetLabelOffset(-0.015);

    TLatex *labelMassEta            = new TLatex(0.13,0.69,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelMassEta, textSizeLabelsPixel,4, 1, 43);

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCent[cent]) continue;
        padWidthEta->cd();
        padWidthEta->SetLogx();

            histo2DAllEtaFWHM->DrawCopy();

            for (Int_t meth = 10; meth > -1; meth--){
                if (graphEtaWidth[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphEtaWidth[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.55, colorDet[meth] , colorDet[meth]);
                    graphEtaWidth[cent][meth]->Draw("p,same,z");
                }
                if (graphEtaWidthMC[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphEtaWidthMC[cent][meth], markerStyleDetMC[meth], markerSizeDetMC[meth]*0.55, colorDetMC[meth] , colorDetMC[meth]);
                    graphEtaWidthMC[cent][meth]->Draw("p,same,z");
                }
            }

            labelLegendAMass->Draw();
            labelMassPerf->Draw();
            labelMassEnergy->SetText(0.13,0.775,Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
            labelMassEnergy->Draw();
            labelMassEta->Draw();

        padMassEta->cd();
        padMassEta->SetLogx();

            histo2DAllEtaMass->DrawCopy();

            for (Int_t meth = 10; meth > -1; meth--){
                if (graphEtaMass[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphEtaMass[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.55, colorDet[meth] , colorDet[meth]);
                    graphEtaMass[cent][meth]->Draw("p,same,z");
                }
                if (graphEtaMassMC[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphEtaMassMC[cent][meth], markerStyleDetMC[meth], markerSizeDetMC[meth]*0.55, colorDetMC[meth] , colorDetMC[meth]);
                    graphEtaMassMC[cent][meth]->Draw("p,same,z");
                }
            }

            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,0.1, kGray);
            labelLegendBMass->Draw();

            padMassLegendEta->cd();
            padMassLegendEta->Clear();
            //****************** first Column **************************************************
            TLatex *textMassPCM2            = new TLatex(columnsLegendMass[0],rowsLegendMass[1],nameMeasGlobalLabel[0]);
            SetStyleTLatex( textMassPCM2, textSizeLabelsPixel,4, 1, 43);
            if (graphEtaMass[cent][0]) textMassPCM2->Draw();
            TLatex *textMassPCMEMCAL       = new TLatex(columnsLegendMass[0],rowsLegendMass[2],nameMeasGlobalLabel[4]);
            SetStyleTLatex( textMassPCMEMCAL, textSizeLabelsPixel,4, 1, 43);
            if (graphEtaMass[cent][4]) textMassPCMEMCAL->Draw();
            TLatex *textMassEMCAL          = new TLatex(columnsLegendMass[0],rowsLegendMass[3],nameMeasGlobalLabel[2]);
            SetStyleTLatex( textMassEMCAL, textSizeLabelsPixel,4, 1, 43);
            if (graphEtaMass[cent][2]) textMassEMCAL->Draw();
            TLatex *textMassPHOS           = new TLatex(columnsLegendMass[3],rowsLegendMass[1],nameMeasGlobalLabel[1]);
            SetStyleTLatex( textMassPHOS, textSizeLabelsPixel,4, 1, 43);
            if (graphEtaMass[cent][1]) textMassPHOS->Draw();
            TLatex *textMassPCMPHOS             = new TLatex(columnsLegendMass[3],rowsLegendMass[2],nameMeasGlobalLabel[3]);
            SetStyleTLatex( textMassPCMPHOS, textSizeLabelsPixel,4, 1, 43);
            if (graphEtaMass[cent][3]) textMassPCMPHOS->Draw();

            //****************** second Column *************************************************
            TLatex *textMassData1           = new TLatex(columnsLegendMass[1],rowsLegendMass[0] ,"Data");
            SetStyleTLatex( textMassData1, textSizeLabelsPixel,4, 1, 43);
            if (graphEtaMass[cent][2] || graphEtaMass[cent][4] || graphEtaMass[cent][0]) textMassData1->Draw();
            TLatex *textMassMC1             = new TLatex(columnsLegendMass[2] ,rowsLegendMass[0],"MC");
            SetStyleTLatex( textMassMC1, textSizeLabelsPixel,4, 1, 43);
            if (graphEtaMass[cent][2] || graphEtaMass[cent][4] || graphEtaMass[cent][0]) textMassMC1->Draw();

            TLatex *textMassData2           = new TLatex(columnsLegendMass[4],rowsLegendMass[0] ,"Data");
            SetStyleTLatex( textMassData2, textSizeLabelsPixel,4, 1, 43);
            if (graphEtaMass[cent][1] || graphEtaMass[cent][3] ) textMassData2->Draw();
            TLatex *textMassMC2             = new TLatex(columnsLegendMass[5] ,rowsLegendMass[0],"MC");
            SetStyleTLatex( textMassMC2, textSizeLabelsPixel,4, 1, 43);
            if (graphEtaMass[cent][1] || graphEtaMass[cent][3] ) textMassMC2->Draw();

            TMarker* markerEtaMass[11]      = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
            TMarker* markerEtaMassMC[11]    = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

            for (Int_t meth = 0; meth< 11; meth++){
                if (graphEtaMass[cent][meth]) markerEtaMass[meth]            = CreateMarkerFromGraph(graphEtaMass[cent][meth], 0.1, 0.1, scaleMarkerMass);
                if (graphEtaMassMC[cent][meth]) markerEtaMassMC[meth]        = CreateMarkerFromGraph(graphEtaMassMC[cent][meth], 0.1, 0.1, scaleMarkerMass);
            }
            if (markerEtaMass[0]) markerEtaMass[0]->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[1]+ offsetMarkerYMass);
            if (markerEtaMass[4]) markerEtaMass[4]->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
            if (markerEtaMass[2]) markerEtaMass[2]->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
            if (markerEtaMass[1]) markerEtaMass[1]->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[1]+ offsetMarkerYMass);
            if (markerEtaMass[3]) markerEtaMass[3]->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);

            if (markerEtaMassMC[0]) markerEtaMassMC[0]->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass-0.01 ,rowsLegendMass[1]+ offsetMarkerYMass);
            if (markerEtaMassMC[4]) markerEtaMassMC[4]->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass-0.01 ,rowsLegendMass[2]+ offsetMarkerYMass);
            if (markerEtaMassMC[2]) markerEtaMassMC[2]->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass-0.01 ,rowsLegendMass[3]+ offsetMarkerYMass);
            if (markerEtaMassMC[1]) markerEtaMassMC[1]->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass-0.01 ,rowsLegendMass[1]+ offsetMarkerYMass);
            if (markerEtaMassMC[3]) markerEtaMassMC[3]->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass-0.01 ,rowsLegendMass[2]+ offsetMarkerYMass);

        canvasMassWidthEta->Update();
        canvasMassWidthEta->Print(Form("%s/Eta_MassAndWidth_%s%s.%s",outputDir.Data(),centArrayOutput[cent].Data(), runArray[cent].Data(), suffix.Data()));

        for (Int_t meth = 0; meth< 11; meth++){
            delete markerEtaMass[meth];
            delete markerEtaMassMC[meth];
        }

        delete textMassPCM2;
        delete textMassPCMEMCAL;
        delete textMassEMCAL;
        delete textMassPHOS;
        delete textMassPCMPHOS;
        delete textMassData1;
        delete textMassMC1;
        delete textMassData2;
        delete textMassMC2;
    }

    Float_t xRangeMassEta[2][11] = {   {0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.43, 0.23, 0.43, 0.43, 0.43},
                                    {25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 25.} };
    Float_t yRangeMassEta[2][11] = {   {505, 505, 505, 505, 505, 505, 505, 505, 505, 505, 505},
                                    {590, 590, 590, 590, 590, 590, 590, 590, 590, 590 ,590} };

    DrawGammaCanvasSettings( canvasMassWidthEta,  0.09, 0.02, 0.01, 0.082);
    TH2F * histo2DCentEtaMass        = new TH2F("histo2DAllEtaMass","histo2DAllEtaMass",2000, minPtEtaPlotting, maxPtEtaPlotting, 1000., 505, 590);
    SetStyleHistoTH2ForGraphs(histo2DCentEtaMass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel,
                              textSizeLabelsRel, 0.9, 1.1, 512, 505);
    histo2DCentEtaMass->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DCentEtaMass->GetYaxis()->SetNdivisions(505);
    histo2DCentEtaMass->GetYaxis()->SetNoExponent(kTRUE);
    histo2DCentEtaMass->GetXaxis()->SetTickLength(0.025);
    histo2DCentEtaMass->GetXaxis()->SetLabelOffset(-0.008);
    TLatex *labelDetSysEtaMass      = new TLatex(0.15,0.82,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysEtaMass, textSizeLabelsRel,4);


    for (Int_t meth = 0; meth < 11; meth++){
        if (!fileNamesMethod[meth].CompareTo("")) continue;

        canvasMassWidthEta->cd();
        canvasMassWidthEta->SetLogx(1);
        histo2DCentEtaMass->GetXaxis()->SetRangeUser(xRangeMassEta[0][meth],xRangeMassEta[1][meth]);
        histo2DCentEtaMass->GetYaxis()->SetRangeUser(yRangeMassEta[0][meth],yRangeMassEta[1][meth]);
        histo2DCentEtaMass->DrawCopy();

        TLegend* legendMass           = GetAndSetLegend2(0.65, 0.95-(5*textSizeLabelsRel), 0.93, 0.95,textSizeLabelsPixel);
        for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
            if (!enableCent[cent]) continue;
            if (graphEtaMass[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphEtaMass[cent][meth], markerStyleCent[cent], markerSizeCent[cent]*0.55, colorCent[cent] , colorCent[cent]);
                graphEtaMass[cent][meth]->Draw("p,same,z");
                legendMass->AddEntry(graphEtaMass[cent][meth],centArray[cent],"p");
            }
        }
        legendMass->Draw();
        labelPerfMass2->Draw();
        labelEnergyMass2->SetText(0.15,0.87,Form("%s", collisionSystempPb.Data()));
        labelEnergyMass2->Draw();
        labelDetSysEtaMass->SetText(0.15,0.82,Form("#eta #rightarrow #gamma#gamma, %s", nameMeasGlobalLabel[meth].Data()));
        labelDetSysEtaMass->Draw();

        canvasMassWidthEta->Update();
        canvasMassWidthEta->Print(Form("%s/Eta_MassPositionData_%s.%s",outputDirSupport.Data(), nameMeasGlobalLabel[meth].Data(), suffix.Data()));

        for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
            if (!enableCent[cent]) continue;
            if (graphEtaMassMC[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphEtaMassMC[cent][meth], markerStyleCentMC[cent], markerSizeCent[cent]*0.55, colorCentMC[cent] , colorCentMC[cent]);
                graphEtaMassMC[cent][meth]->Draw("p,same,z");
            }
        }

        canvasMassWidthEta->Update();
        canvasMassWidthEta->Print(Form("%s/Eta_MassPositionDataAndMC_%s.%s",outputDirSupport.Data(), nameMeasGlobalLabel[meth].Data(), suffix.Data()));
    }


    // **********************************************************************************************************************
    // ******************************** Cross section for pi0 single measurement 2.76TeV ************************************
    // **********************************************************************************************************************
    TCanvas* canvasInvYieldPi0          = new TCanvas("canvasInvYieldPi0","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasInvYieldPi0, 0.16, 0.02, 0.02, 0.08);
    canvasInvYieldPi0->SetLogx();
    canvasInvYieldPi0->SetLogy();

    TH2F * histo2DYieldPi0              = new TH2F("histo2DYieldPi0","histo2DYieldPi0",11000,minPtPi0Plotting, maxPtPi0Plotting,1000,6e-11,3e1);
    SetStyleHistoTH2ForGraphs(histo2DYieldPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.65);
    histo2DYieldPi0->GetXaxis()->SetMoreLogLabels();
    histo2DYieldPi0->GetXaxis()->SetLabelOffset(-0.01);

    TLatex *labelEnergyXSectionPi0  = new TLatex(0.94,0.92,collisionSystempPb.Data());
    SetStyleTLatex( labelEnergyXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);
    TLatex *labelDetSysXSectionPi0  = new TLatex(0.94,0.88,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCentComb[cent]) continue;
        canvasInvYieldPi0->cd();
        histo2DYieldPi0->Draw("copy");

            for (Int_t meth = 10; meth> -1; meth--){
                if (graphIndPi0InvYieldSys[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphIndPi0InvYieldSys[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.75, colorDet[meth] , colorDet[meth], widthLinesBoxes, kTRUE);
                    graphIndPi0InvYieldSys[cent][meth]->Draw("E2same");
                }
            }
            for (Int_t meth = 10; meth> -1; meth--){
                if (graphIndPi0InvYieldStat[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphIndPi0InvYieldStat[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.75, colorDet[meth] , colorDet[meth]);
                    graphIndPi0InvYieldStat[cent][meth]->Draw("p,same,z");
                }
            }
            labelEnergyXSectionPi0->SetText(0.94,0.92,Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
            labelEnergyXSectionPi0->Draw();
            labelDetSysXSectionPi0->Draw();


            TLegend* legendXSectionPi0      = GetAndSetLegend2(0.20, 0.13,0.40,0.13+(0.035*nTotMeasPi0[cent]), 0.035, 1, "", 42, 0);
            for (Int_t meth = 0; meth < 11; meth++){
                if (graphIndPi0InvYieldSys[cent][meth]) legendXSectionPi0->AddEntry(graphIndPi0InvYieldSys[cent][meth],nameMeasGlobalLabel[meth],"fp");
            }
            legendXSectionPi0->Draw();

        canvasInvYieldPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), runArray[cent].Data(), suffix.Data()));

        canvasInvYieldPi0->cd();
        histo2DYieldPi0->Draw("copy");
            for (Int_t meth = 10; meth> -1; meth--){
                        if (graphIndPi0InvYieldSys[cent][meth]){
                            graphIndPi0InvYieldSys[cent][meth]->Draw("E2same");
                }
            }
            if (graphCombPi0InvYieldSys[cent]){
                DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys[cent], markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
                graphCombPi0InvYieldSys[cent]->Draw("E2same");
                legendXSectionPi0->AddEntry(graphCombPi0InvYieldSys[cent],"comb","fp");
            }
            for (Int_t meth = 10; meth> -1; meth--){
                if (graphIndPi0InvYieldSys[cent][meth]){
                    graphIndPi0InvYieldStat[cent][meth]->Draw("p,same,z");
                }
            }
            if (graphCombPi0InvYieldStat[cent]){
                DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
                graphCombPi0InvYieldStat[cent]->Draw("p,same,z");
            }
            labelEnergyXSectionPi0->Draw();
            labelDetSysXSectionPi0->Draw();

            legendXSectionPi0->Draw();

        canvasInvYieldPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems_Comb_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), runArray[cent].Data(), suffix.Data()));
    }

    // **********************************************************************************************************************
    // ******************************** Cross section for eta single measurement 2.76TeV ************************************
    // **********************************************************************************************************************
    TCanvas* canvasInvYieldEta          = new TCanvas("canvasInvYieldEta","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasInvYieldEta, 0.16, 0.02, 0.02, 0.08);
    canvasInvYieldEta->SetLogx();
    canvasInvYieldEta->SetLogy();

    TH2F * histo2DYieldEta              = new TH2F("histo2DYieldEta","histo2DYieldEta",11000,minPtEtaPlotting, maxPtEtaPlotting,1000,6e-10,3e-0);
    SetStyleHistoTH2ForGraphs(histo2DYieldEta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.65);
    histo2DYieldEta->GetXaxis()->SetMoreLogLabels();
    histo2DYieldEta->GetXaxis()->SetLabelOffset(-0.01);

    TLatex *labelEnergyXSectionEta  = new TLatex(0.94,0.92,collisionSystempPb.Data());
    SetStyleTLatex( labelEnergyXSectionEta, 0.035,4, 1, 42, kTRUE, 31);
    TLatex *labelDetSysXSectionEta  = new TLatex(0.94,0.88,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysXSectionEta, 0.035,4, 1, 42, kTRUE, 31);

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCentComb[cent]) continue;
        canvasInvYieldEta->cd();
        histo2DYieldEta->Draw("copy");

        for (Int_t meth = 10; meth> -1; meth--){
            if (graphIndEtaInvYieldSys[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphIndEtaInvYieldSys[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.75, colorDet[meth] , colorDet[meth], widthLinesBoxes, kTRUE);
                graphIndEtaInvYieldSys[cent][meth]->Draw("E2same");
            }
        }
        for (Int_t meth = 10; meth> -1; meth--){
            if (graphIndEtaInvYieldStat[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphIndEtaInvYieldStat[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.75, colorDet[meth] , colorDet[meth]);
                graphIndEtaInvYieldStat[cent][meth]->Draw("p,same,z");
            }
        }
        labelEnergyXSectionEta->SetText(0.94,0.92,Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
        labelEnergyXSectionEta->Draw();
        labelDetSysXSectionEta->Draw();


        TLegend* legendXSectionEta      = GetAndSetLegend2(0.20, 0.13,0.40,0.13+(0.035*nTotMeasEta[cent]), 0.035, 1, "", 42, 0);
        for (Int_t meth = 0; meth < 11; meth++){
            if (graphIndEtaInvYieldSys[cent][meth]) legendXSectionEta->AddEntry(graphIndEtaInvYieldSys[cent][meth],nameMeasGlobalLabel[meth],"fp");
        }
        legendXSectionEta->Draw();

        canvasInvYieldEta->SaveAs(Form("%s/Eta_InvYieldCompAllSystems_%s%s.%s",outputDir.Data(),centArrayOutput[cent].Data(), runArray[cent].Data(), suffix.Data()));

        canvasInvYieldEta->cd();
        histo2DYieldEta->Draw("copy");
            for (Int_t meth = 10; meth> -1; meth--){
                        if (graphIndEtaInvYieldSys[cent][meth]){
                            graphIndEtaInvYieldSys[cent][meth]->Draw("E2same");
                }
            }
            if (graphCombEtaInvYieldSys[cent]){
                DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys[cent], markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
                graphCombEtaInvYieldSys[cent]->Draw("E2same");
                legendXSectionEta->AddEntry(graphCombEtaInvYieldSys[cent],"comb","fp");
            }
            for (Int_t meth = 10; meth> -1; meth--){
                if (graphIndEtaInvYieldSys[cent][meth]){
                    graphIndEtaInvYieldStat[cent][meth]->Draw("p,same,z");
                }
            }
            if (graphCombEtaInvYieldStat[cent]){
                DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
                graphCombEtaInvYieldStat[cent]->Draw("p,same,z");
            }
            labelEnergyXSectionEta->Draw();
            labelDetSysXSectionEta->Draw();

            legendXSectionEta->Draw();

        canvasInvYieldEta->SaveAs(Form("%s/Eta_InvYieldCompAllSystems_Comb_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), runArray[cent].Data(), suffix.Data()));

    }

    // **********************************************************************************************************************
    // ******************************** Acceptance * Efficiency for pi0 single measurement **********************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 55;
    textSizeLabelsRel                   = ((Double_t)textSizeLabelsPixel)/1200;
    cout << textSizeLabelsRel << endl;

    TCanvas* canvasAcceptanceTimesEff   = new TCanvas("canvasAcceptanceTimesEff", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasAcceptanceTimesEff,  0.1, 0.01, 0.015, 0.095);
    canvasAcceptanceTimesEff->SetLogy(1);
    canvasAcceptanceTimesEff->SetLogx(1);

    TH1F * histo1DAccEff            = new TH1F("histo1DAccEff", "histo1DAccEff",1000, minPtPi0Plotting,  maxPtPi0Plotting);
    SetStyleHistoTH1ForGraphs(  histo1DAccEff, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec}#upoint#it{BR}"),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
    histo1DAccEff->GetYaxis()->SetRangeUser(7e-6, 2e-0 );
    histo1DAccEff->GetYaxis()->SetLabelOffset(0.001);
    histo1DAccEff->GetXaxis()->SetLabelOffset(-0.01);
    histo1DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
    TLatex *labelPerfEffi           = new TLatex(0.137,0.92,"ALICE performance");
    SetStyleTLatex( labelPerfEffi, textSizeLabelsRel,4);
    TLatex *labelEnergyEffi         = new TLatex(0.137,0.87,collisionSystempPb.Data());
    SetStyleTLatex( labelEnergyEffi, textSizeLabelsRel,4);
    TLatex *labelDetSysEffiPi0      = new TLatex(0.137,0.82,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysEffiPi0, textSizeLabelsRel,4);

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCent[cent]) continue;
        canvasAcceptanceTimesEff->cd();
        histo1DAccEff->DrawCopy();

            for (Int_t meth = 0; meth < 11; meth++){
                if (graphPi0EffTimesAcc[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphPi0EffTimesAcc[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.55, colorDet[meth] , colorDet[meth]);
                    graphPi0EffTimesAcc[cent][meth]->Draw("p,same,z");
                }
            }

            TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(nTotMeasPi0[cent]*textSizeLabelsRel),textSizeLabelsPixel);
            if (graphPi0EffTimesAcc[cent][0]) legendEffiAccPi0->AddEntry(graphPi0EffTimesAcc[cent][0],nameMeasGlobalLabel[0],"p");
            if (graphPi0EffTimesAcc[cent][4]) legendEffiAccPi0->AddEntry(graphPi0EffTimesAcc[cent][4],nameMeasGlobalLabel[4],"p");
            if (graphPi0EffTimesAcc[cent][2]) legendEffiAccPi0->AddEntry(graphPi0EffTimesAcc[cent][2],nameMeasGlobalLabel[2],"p");
            if (graphPi0EffTimesAcc[cent][1]) legendEffiAccPi0->AddEntry(graphPi0EffTimesAcc[cent][1],nameMeasGlobalLabel[1],"p");
            if (graphPi0EffTimesAcc[cent][3]) legendEffiAccPi0->AddEntry(graphPi0EffTimesAcc[cent][3],nameMeasGlobalLabel[3],"p");
            if (graphPi0EffTimesAcc[cent][5]) legendEffiAccPi0->AddEntry(graphPi0EffTimesAcc[cent][5],nameMeasGlobalLabel[5],"p");
            legendEffiAccPi0->Draw();

            labelPerfEffi->Draw();
            labelEnergyEffi->SetText(0.137,0.87,Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
            labelEnergyEffi->Draw();
            labelDetSysEffiPi0->Draw();

        canvasAcceptanceTimesEff->Update();
        canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff_%s%s.%s",outputDir.Data(),centArrayOutput[cent].Data(), runArray[cent].Data(), suffix.Data()));
    }

    Double_t rangMinEffi[11]                = { 7e-5, 1e-3, 8e-4, 9e-5, 9e-5, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4 };
    Double_t rangMaxEffi[11]                = { 5e-2, 3e-1, 2e-0, 5e-2, 3e-1, 2e-0, 2e-0, 2e-0, 2e-0, 2e-0, 2e-0 };

    for (Int_t meth = 0; meth < 11; meth++){
        if (!fileNamesMethod[meth].CompareTo("")) continue;
        canvasAcceptanceTimesEff->cd();
        histo1DAccEff->GetYaxis()->SetRangeUser(rangMinEffi[meth], rangMaxEffi[meth] );
        histo1DAccEff->DrawCopy();

        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(5*textSizeLabelsRel),textSizeLabelsPixel);
        for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
            if (!enableCent[cent]) continue;
            if (graphPi0EffTimesAcc[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphPi0EffTimesAcc[cent][meth], markerStyleCent[cent], markerSizeCent[cent]*0.55, colorCent[cent] , colorCent[cent]);
                graphPi0EffTimesAcc[cent][meth]->Draw("p,same,z");
                legendEffiAccPi0->AddEntry(graphPi0EffTimesAcc[cent][meth],centArray[cent],"p");
            }
        }
        legendEffiAccPi0->Draw();

        labelPerfEffi->Draw();
        labelEnergyEffi->SetText(0.137,0.87,Form("%s", collisionSystempPb.Data()));
        labelEnergyEffi->Draw();
        labelDetSysEffiPi0->SetText(0.137,0.82,Form("#pi^{0} #rightarrow #gamma#gamma, %s", nameMeasGlobalLabel[meth].Data()));
        labelDetSysEffiPi0->Draw();

        canvasAcceptanceTimesEff->Update();
        canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff_%s.%s",outputDirSupport.Data(), nameMeasGlobalLabel[meth].Data(), suffix.Data()));

    }

    TH1F * histo1DAccEffEta            = new TH1F("histo1DAccEffEta", "histo1DAccEffEta",1000, minPtEtaPlotting,  maxPtEtaPlotting);
    SetStyleHistoTH1ForGraphs(  histo1DAccEffEta, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec}#upoint#it{BR}"),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
    histo1DAccEffEta->GetYaxis()->SetRangeUser(9e-5, 2e-0 );
    histo1DAccEffEta->GetYaxis()->SetLabelOffset(0.001);
    histo1DAccEffEta->GetXaxis()->SetLabelOffset(-0.01);
    histo1DAccEffEta->GetXaxis()->SetMoreLogLabels(kTRUE);

    TLatex *labelDetSysEffiEta      = new TLatex(0.137,0.82,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysEffiEta, textSizeLabelsRel,4);

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCent[cent]) continue;
        canvasAcceptanceTimesEff->cd();
        histo1DAccEffEta->GetYaxis()->SetRangeUser(9e-6, 2e-0 );
        histo1DAccEffEta->DrawCopy();

        for (Int_t meth = 0; meth < 11; meth++){
            if (graphEtaEffTimesAcc[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphEtaEffTimesAcc[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.55, colorDet[meth] , colorDet[meth]);
                graphEtaEffTimesAcc[cent][meth]->Draw("p,same,z");
            }
        }

        TLegend* legendEffiAccEta           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(nTotMeasEta[cent]*textSizeLabelsRel),textSizeLabelsPixel);
        if (graphEtaEffTimesAcc[cent][0]) legendEffiAccEta->AddEntry(graphEtaEffTimesAcc[cent][0],nameMeasGlobalLabel[0],"p");
        if (graphEtaEffTimesAcc[cent][4]) legendEffiAccEta->AddEntry(graphEtaEffTimesAcc[cent][4],nameMeasGlobalLabel[4],"p");
        if (graphEtaEffTimesAcc[cent][2]) legendEffiAccEta->AddEntry(graphEtaEffTimesAcc[cent][2],nameMeasGlobalLabel[2],"p");
        if (graphEtaEffTimesAcc[cent][1]) legendEffiAccEta->AddEntry(graphEtaEffTimesAcc[cent][1],nameMeasGlobalLabel[1],"p");
        if (graphEtaEffTimesAcc[cent][3]) legendEffiAccEta->AddEntry(graphEtaEffTimesAcc[cent][3],nameMeasGlobalLabel[3],"p");
        legendEffiAccEta->Draw();

        labelPerfEffi->Draw();
        labelEnergyEffi->SetText(0.137,0.87,Form("%s %s", centArray[cent].Data(), collisionSystempPb.Data()));
        labelEnergyEffi->Draw();
        labelDetSysEffiEta->Draw();

        canvasAcceptanceTimesEff->Update();
        canvasAcceptanceTimesEff->Print(Form("%s/Eta_AcceptanceTimesEff_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), runArray[cent].Data(), suffix.Data()));
    }

    Double_t rangMinEffiEta[11]                = { 9e-4, 1e-3, 8e-4, 9e-4, 9e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4 };
    Double_t rangMaxEffiEta[11]                = { 5e-2, 3e-1, 2e-0, 5e-2, 3e-1, 2e-0, 2e-0, 2e-0, 2e-0, 2e-0, 2e-0 };

    for (Int_t meth = 0; meth < 11; meth++){
        if (!fileNamesMethod[meth].CompareTo("")) continue;
        canvasAcceptanceTimesEff->cd();
        histo1DAccEffEta->GetYaxis()->SetRangeUser(rangMinEffiEta[meth], rangMaxEffiEta[meth] );
        histo1DAccEffEta->DrawCopy();

        TLegend* legendEffiAccEta           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(5*textSizeLabelsRel),textSizeLabelsPixel);
        for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
            if (!enableCent[cent]) continue;
            if (graphEtaEffTimesAcc[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphEtaEffTimesAcc[cent][meth], markerStyleCent[cent], markerSizeCent[cent]*0.55, colorCent[cent] , colorCent[cent]);
                graphEtaEffTimesAcc[cent][meth]->Draw("p,same,z");
                legendEffiAccEta->AddEntry(graphEtaEffTimesAcc[cent][meth],centArray[cent],"p");
            }
        }
        legendEffiAccEta->Draw();

        labelPerfEffi->Draw();
        labelEnergyEffi->SetText(0.137,0.87,Form("%s", collisionSystempPb.Data()));
        labelEnergyEffi->Draw();
        labelDetSysEffiEta->SetText(0.137,0.82,Form("#pi^{0} #rightarrow #gamma#gamma, %s", nameMeasGlobalLabel[meth].Data()));
        labelDetSysEffiEta->Draw();

        canvasAcceptanceTimesEff->Update();
        canvasAcceptanceTimesEff->Print(Form("%s/Eta_AcceptanceTimesEff_%s.%s",outputDirSupport.Data(), nameMeasGlobalLabel[meth].Data(), suffix.Data()));

    }

    delete labelPerfEffi;
    delete labelEnergyEffi;
    delete labelDetSysEffiPi0;
    delete labelDetSysEffiEta;
    delete canvasAcceptanceTimesEff;

    // **********************************************************************************************************************
    // ******************************** effective secondary correction drawing for different methods ************************
    // **********************************************************************************************************************
    TCanvas* canvasEffectiveSecCorr     = new TCanvas("canvasEffectiveSecCorr", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasEffectiveSecCorr,  0.1, 0.01, 0.04, 0.095);
    canvasEffectiveSecCorr->SetLogx(1);

        TH1F * histo1DEffSecCorr;
        histo1DEffSecCorr               = new TH1F("histo1DEffSecCorr", "histo1DEffSecCorr",1000, minPtPi0Plotting,  maxPtPi0Plotting);
        SetStyleHistoTH1ForGraphs(  histo1DEffSecCorr, "#it{p}_{T} (GeV/#it{c})","R_{K^{0}_{s}}",
                                    0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
        histo1DEffSecCorr->GetYaxis()->SetRangeUser(0, 10 );
        histo1DEffSecCorr->GetXaxis()->SetLabelOffset(-0.01);
        histo1DEffSecCorr->GetXaxis()->SetMoreLogLabels(kTRUE);
        TLatex *labelPerfSecCorr        = new TLatex(0.15,0.90,"ALICE performance");
        SetStyleTLatex( labelPerfSecCorr, textSizeLabelsRel,4);
        TLatex *labelEnergySecCorr      = new TLatex(0.15,0.85,collisionSystempPb.Data());
        SetStyleTLatex( labelEnergySecCorr, textSizeLabelsRel,4);
        TLatex *labelDetSysSecCorrPi0   = new TLatex(0.15,0.80,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysSecCorrPi0, textSizeLabelsRel,4);

        for (Int_t cent= 0; cent < maxCentRun1+maxCentRun2; cent++ ){
            if (!enableCent[cent]) continue;
            for (Int_t k = 0; k < 4; k++){
                Bool_t plotCorr             = kFALSE;
                Int_t nCorrAvail            = 0;
                for (Int_t meth = 0; meth < 11; meth++){
                    if (haveEffSecCorr[cent][k][meth]){
                        nCorrAvail++;
                        plotCorr            = kTRUE;
                    }
                }
                TLegend* legendEffSecCorrPi0        = GetAndSetLegend2(0.65, 0.925-(nCorrAvail*textSizeLabelsRel), 0.93, 0.925,textSizeLabelsPixel);
                if (plotCorr){
                    histo1DEffSecCorr->GetYaxis()->SetTitle(Form("R_{%s}",nameSecPi0SourceLabel[k].Data()));
                    histo1DEffSecCorr->GetYaxis()->SetRangeUser(0, maxSecCorr[k]);
                    histo1DEffSecCorr->DrawCopy();
                    for (Int_t meth = 0; meth < 11; meth++){
                        if (haveEffSecCorr[cent][k][meth]){
                            DrawGammaSetMarkerTGraphAsym(graphPi0EffSecCorrFromX[cent][k][meth], markerStyleDet[meth], markerSizeDet[meth]*0.55, colorDet[meth] , colorDet[meth]);
                            graphPi0EffSecCorrFromX[cent][k][meth]->Draw("p,same,z");
                            legendEffSecCorrPi0->AddEntry(graphPi0EffSecCorrFromX[cent][k][meth],nameMeasGlobalLabel[meth],"p");
                        }
                    }
                    legendEffSecCorrPi0->Draw();
                    labelPerfSecCorr->Draw();
                    labelEnergySecCorr->Draw();
                    labelDetSysSecCorrPi0->Draw();

                    canvasEffectiveSecCorr->Update();
                    canvasEffectiveSecCorr->Print(Form("%s/Pi0_EffectiveSecCorr_%s%s_%s.%s",outputDir.Data(), nameSecPi0SourceRead[k].Data(), centArrayOutput[cent].Data(), runArray[cent].Data(),
                                                       suffix.Data()));
                }
            }
        }
    delete canvasEffectiveSecCorr;

    fileFitsOutput << "*******************************************************************************************" << endl;
    fileFitsOutput << "****************************** Power law fit pi0 ******************************************" << endl;
    fileFitsOutput << "*******************************************************************************************" << endl;
    TF1* fitPowInvYieldPi0Tot[16]   = { NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL,
                                        NULL, NULL, NULL, NULL, NULL,   NULL };
    TF1* fitPowInvYieldPi0Stat[16]   = { NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,
                                        NULL, NULL, NULL, NULL, NULL,   NULL };
    TF1* fitOHagInvYieldPi0Tot[16]   = { NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,
                                        NULL, NULL, NULL, NULL, NULL,   NULL };
    TF1* fitPowInvYieldEtaTot[16]   = { NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL,
                                        NULL, NULL, NULL, NULL, NULL,   NULL };
    TF1* fitPowInvYieldEtaStat[16]   = { NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,
                                        NULL, NULL, NULL, NULL, NULL,   NULL };
    TF1* fitOHagInvYieldEtaTot[16]   = { NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,
                                        NULL, NULL, NULL, NULL, NULL,   NULL };

    for (Int_t cent= 0; cent < maxCentRun1+maxCentRun2; cent++ ){
        if (!enableCentComb[cent]) continue;
        cout << centArray[cent].Data() << ":"  << endl;
        fileFitsOutput << centArray[cent].Data() << ":"  << endl;
        fitPowInvYieldPi0Tot[cent]   = FitObject("powPure","fitPowInvYieldPi0Tot","Pi0",graphCombPi0InvYieldTot[cent],4,20. ,NULL,"QNRMEX0+","", kFALSE);
        cout << WriteParameterToFile(fitPowInvYieldPi0Tot[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldPi0Tot[cent])<< endl;
        fitPowInvYieldPi0Stat[cent]   = FitObject("powPure","fitPowInvYieldPi0","Pi0",graphCombPi0InvYieldStat[cent],4,20. ,NULL,"QNRMEX0+","", kFALSE);
        cout << WriteParameterToFile(fitPowInvYieldPi0Stat[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldPi0Stat[cent])<< endl;
        fitOHagInvYieldPi0Tot[cent]   = FitObject("oHag","fitOHagInvYieldPi0","Pi0",graphCombPi0InvYieldTot[cent],0.3,20. ,NULL,"QNRMEX0+","", kFALSE);
        cout << WriteParameterToFile(fitOHagInvYieldPi0Tot[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitOHagInvYieldPi0Tot[cent])<< endl;
        graphCombPi0InvYieldStatWOXErr[cent]        = (TGraphAsymmErrors*)graphCombPi0InvYieldStat[cent]->Clone("graphCombPi0InvYieldStatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphCombPi0InvYieldStatWOXErr[cent]);
        graphRatioPi0CombCombFitStatWOXErr[cent]    = (TGraphAsymmErrors*)graphRatioPi0CombCombFitStat[cent]->Clone("graphRatioPi0CombCombFitStatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatWOXErr[cent]);

        fitPowInvYieldEtaTot[cent]   = FitObject("powPure","fitPowInvYieldEtaTot","Eta",graphCombEtaInvYieldTot[cent],4,20. ,NULL,"QNRMEX0+","", kFALSE);
        cout << WriteParameterToFile(fitPowInvYieldEtaTot[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldEtaTot[cent])<< endl;
        fitPowInvYieldEtaStat[cent]   = FitObject("powPure","fitPowInvYieldEta","Eta",graphCombEtaInvYieldStat[cent],4,20. ,NULL,"QNRMEX0+","", kFALSE);
        cout << WriteParameterToFile(fitPowInvYieldEtaStat[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldEtaStat[cent])<< endl;
        fitOHagInvYieldEtaTot[cent]   = FitObject("oHag","fitOHagInvYieldEta","Eta",graphCombEtaInvYieldTot[cent],0.3,20. ,NULL,"QNRMEX0+","", kFALSE);
        cout << WriteParameterToFile(fitOHagInvYieldEtaTot[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitOHagInvYieldEtaTot[cent])<< endl;
        graphCombEtaInvYieldStatWOXErr[cent]        = (TGraphAsymmErrors*)graphCombEtaInvYieldStat[cent]->Clone("graphCombEtaInvYieldStatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphCombEtaInvYieldStatWOXErr[cent]);
        graphRatioEtaCombCombFitStatWOXErr[cent]    = (TGraphAsymmErrors*)graphRatioEtaCombCombFitStat[cent]->Clone("graphRatioEtaCombCombFitStatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombCombFitStatWOXErr[cent]);
    }

//     fileFitsOutput << "NLO pp" << endl;
//     TF1* fitPowInvYieldPi0NLOPP         = FitObject("powPure","fitPowInvYieldPi0NLOPP","Pi0",graphNLODSS14Pi0PP,4,20. ,NULL,"QNRMEX0+","", kFALSE);
//     fileFitsOutput << WriteParameterToFile(fitPowInvYieldPi0NLOPP)<< endl;
//
//     fileFitsOutput << "NLO pPb nCTEQ" << endl;
//     TF1* fitPowInvYieldPi0NLOpPbnCTEQ   = FitObject("powPure","fitPowInvYieldPi0NLOpPbnCTEQ","Pi0",graphNLODSS14nPDFPi0,4,20. ,NULL,"QNRMEX0+","", kFALSE);
//     fileFitsOutput << WriteParameterToFile(fitPowInvYieldPi0NLOpPbnCTEQ)<< endl;
//
//     fileFitsOutput << "NLO pp EPPS16" << endl;
//     TF1* fitPowInvYieldPi0NLOpPbEPPS    = FitObject("powPure","fitPowInvYieldPi0NLOpPbEPPS","Pi0",graphNLODSS14nPDFEPPS16Pi0,4,20. ,NULL,"QNRMEX0+","", kFALSE);
//     fileFitsOutput << WriteParameterToFile(fitPowInvYieldPi0NLOpPbEPPS)<< endl;
//

//     // **********************************************************************************************************************
//     // ******************************************* Comparison to theory calculations Pi0 ************************************
//     // **********************************************************************************************************************
//     textSizeLabelsPixel                     = 48;
//
//     TCanvas* canvasRatiopPb                 = new TCanvas("canvasRatiopPb","",200,10,1350,900);  // gives the page size
//     DrawGammaCanvasSettings( canvasRatiopPb,  0.12, 0.01, 0.01, 0.11);
//     canvasRatiopPb->cd();
//     canvasRatiopPb->SetLogx();
//
//         textsizeLabelspPb                   = 0;
//         if (canvasRatiopPb->XtoPixel(canvasRatiopPb->GetX2()) <canvasRatiopPb->YtoPixel(canvasRatiopPb->GetY1()) ){
//             textsizeLabelspPb               = (Double_t)textSizeLabelsPixel/canvasRatiopPb->XtoPixel(canvasRatiopPb->GetX2()) ;
//         } else {
//             textsizeLabelspPb               = (Double_t)textSizeLabelsPixel/canvasRatiopPb->YtoPixel(canvasRatiopPb->GetY1());
//         }
//         cout << textsizeLabelspPb << endl;
//
//         TH2F * ratio2DTheorypPb             = new TH2F("ratio2DTheorypPb","ratio2DTheorypPb",1000,minPtPi0Plotting, 25.,1000,0.2,2.95);
//         SetStyleHistoTH2ForGraphs(  ratio2DTheorypPb, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelspPb, textsizeLabelspPb,
//                                     0.85*textsizeLabelspPb,textsizeLabelspPb, 0.9, 0.95, 510, 505);
//         ratio2DTheorypPb->GetYaxis()->SetMoreLogLabels(kTRUE);
//         ratio2DTheorypPb->GetYaxis()->SetNdivisions(510);
//         ratio2DTheorypPb->GetYaxis()->SetNoExponent(kTRUE);
//         ratio2DTheorypPb->GetXaxis()->SetMoreLogLabels(kTRUE);
//         ratio2DTheorypPb->GetXaxis()->SetNoExponent(kTRUE);
//         ratio2DTheorypPb->GetXaxis()->SetLabelFont(42);
//         ratio2DTheorypPb->GetYaxis()->SetLabelFont(42);
//         ratio2DTheorypPb->DrawCopy();
//
//
// //         SetStyleHisto(histoRatioPi0DPMJetToFit, widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );
// //         SetStyleHisto(histoRatioPi0HIJINGToFit, widthCommonFit*1.5, styleLineHIJING, colorHIJING );
// //         SetStyleHisto(histoRatioPi0EPOS3ToFit, widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );
// //         DrawGammaSetMarkerTGraphErr(graphRatioPi0EPOS3ToFit, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
// //         DrawGammaSetMarkerTGraphErr(graphRatioPi0McGillToFit, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
//
// //         histoRatioPi0EPOS3ToFit->Draw("same,hist,l");
// //         graphRatioPi0EPOS3ToFit->Draw("same,3");
// //         graphRatioPi0McGillToFit->Draw("same,3");
// //         histoRatioPi0DPMJetToFit->Draw("same,hist,l");
// //         histoRatioPi0HIJINGToFit->Draw("same,hist,l");
//
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr, markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
//         graphRatioPi0CombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
//         graphRatioPi0CombCombFitSys[cent]->SetLineWidth(0);
//         graphRatioPi0CombCombFitSys[cent]->Draw("2,same");
//         graphRatioPi0CombCombFitStatWOXErr->Draw("p,same");
//
        TBox* boxErrorNormRatioPi0          = CreateBoxConv(kGray+1, 0.28, 1.-(0.031 ), 0.32, 1.+(0.031));
        TBox* boxErrorNormRatioEta          = CreateBoxConv(kGray+1, 0.48, 1.-(0.031 ), 0.52, 1.+(0.031));
//         boxErrorNormRatioPi0->Draw();
//         DrawGammaLines(minPtPi0Plotting, 25.,1., 1.,0.1,kGray);
//
//         TLegend* legendRatioTheorypPbNew= GetAndSetLegend2(0.69,0.91-0.85*textsizeLabelspPb*6,0.94,0.91, 0.85*textSizeLabelsPixel);
//         legendRatioTheorypPbNew->AddEntry(graphRatioPi0CombCombFitSys[cent],"#pi^{0} ALICE","pf");
// //         legendRatioTheorypPbNew->AddEntry(histoRatioPi0DPMJetToFit,  "DPMJet", "l");
// //         legendRatioTheorypPbNew->AddEntry(histoRatioPi0HIJINGToFit,  "HIJING", "l");
// //         legendRatioTheorypPbNew->AddEntry(histoRatioPi0EPOS3ToFit,  "EPOS", "l");
// //         legendRatioTheorypPbNew->AddEntry(graphRatioPi0McGillToFit,  "Shen #it{et al.}", "f");
//         legendRatioTheorypPbNew->Draw();
//
//         TLatex *labelRatioTheory            = new TLatex(0.97,0.925,collisionSystempPb.Data());
//         SetStyleTLatex( labelRatioTheory, 0.85*textsizeLabelspPb,4, 1, 42, kTRUE, 31);
//         labelRatioTheory->Draw();
//
//         ratio2DTheorypPb->Draw("same,axis");
//     canvasRatiopPb->Update();
//     canvasRatiopPb->Print(Form("%s/Pi0_RatioTheoryToData_pPb.%s",outputDir.Data(),suffix.Data()));
//
//
    //*************************************************************************************************************
    //***************************** Paper plot invyield and ratios ************************************************
    //*************************************************************************************************************

    Double_t arrayBoundariesX1_XSec[2];
    Double_t arrayBoundariesY1_XSec[6];
    Double_t relativeMarginsXXSec[3];
    Double_t relativeMarginsYXSec[3];
    textSizeLabelsPixel = 48;
    ReturnCorrectValuesForCanvasScaling(1250,2000, 1, 5,0.15, 0.005, 0.003,0.05,arrayBoundariesX1_XSec,arrayBoundariesY1_XSec,relativeMarginsXXSec,relativeMarginsYXSec);

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
        textsizeFacXSecUp                   = (Double_t)1./padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
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

    TLegend* legendXsectionPaper    = GetAndSetLegend2(0.19, 0.03, 0.7, 0.03+0.05*(maxCentRun1+maxCentRun2), textSizeLabelsPixel, 2, "", 43, 0.25);
    TGraphErrors* dummyNorm     = new TGraphErrors(1);
    dummyNorm->SetPoint(0,1,0.0251);
    DrawGammaSetMarkerTGraphErr(dummyNorm, 0, 0, kGray+1, kGray+1, widthLinesBoxes, kTRUE, kGray+1);

    //*************************************************************************************************************
    //***************************** Paper plot invyield pi0 with soft-physics *************************************
    //*************************************************************************************************************
    padInvSectionSpec->cd();
    padInvSectionSpec->SetLogy(1);
    padInvSectionSpec->SetLogx(1);
        SetStyleHistoTH2ForGraphs(histo2DYieldPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                                0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.257/(textsizeFacXSecUp*marginXSec));
        histo2DYieldPi0->GetXaxis()->SetMoreLogLabels();
        histo2DYieldPi0->GetXaxis()->SetLabelOffset(+0.01);
        histo2DYieldPi0->Draw();

//         SetStyleHisto(histoEPOS3Pi0, widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );
//         DrawGammaSetMarkerTGraphErr(graphEPOS3Pi0, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
//         histoEPOS3Pi0->GetXaxis()->SetRangeUser(0.3,20);
//         histoEPOS3Pi0->Draw("same,hist,l");
//         while(graphEPOS3Pi0->GetX()[0] < 0.3)
//             graphEPOS3Pi0->RemovePoint(0);
//         while(graphEPOS3Pi0->GetX()[graphEPOS3Pi0->GetN()-1] > 20)
//             graphEPOS3Pi0->RemovePoint(graphEPOS3Pi0->GetN()-1);
//         graphEPOS3Pi0->Draw("same,3");

//         DrawGammaSetMarkerTGraphErr(graphMcGillPi0, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
//         while(graphMcGillPi0->GetX()[0] < 0.3)
//             graphMcGillPi0->RemovePoint(0);
//         graphMcGillPi0->Draw("same,3");

//         SetStyleHisto(histoDPMJetPi0, widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );
//         histoDPMJetPi0->GetXaxis()->SetRangeUser(0.3,20);
//         histoDPMJetPi0->Draw("same,hist,l");
//         SetStyleHisto(histoHIJINGPi0, widthCommonFit*1.5, styleLineHIJING, colorHIJING );
//         histoHIJINGPi0->GetXaxis()->SetRangeUser(0.3,20);
//         histoHIJINGPi0->Draw("same,hist,l");

        for (Int_t cent= 0; cent < maxCentRun1+maxCentRun2; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
            graphCombPi0InvYieldSys[cent]->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent]);
            graphCombPi0InvYieldStatWOXErr[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTF1( fitTCMInvYieldPi0[cent], 7, 2, colorCentMC[cent]);
            fitTCMInvYieldPi0[cent]->Draw("same");
            legendXsectionPaper->AddEntry(graphCombPi0InvYieldSys[cent],centArray[cent].Data(),"pf");
            legendXsectionPaper->AddEntry(fitTCMInvYieldPi0[cent],"TCM fit","l");
        }

        TLatex *labelEnergyXSectionPaper= new TLatex(0.95, 0.91, collisionSystempPb.Data());
        SetStyleTLatex( labelEnergyXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
        labelEnergyXSectionPaper->Draw();
        TLatex *labelALICEXSectionPaper= new TLatex(0.95,0.867,textALICE.Data());
        SetStyleTLatex( labelALICEXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
        labelALICEXSectionPaper->Draw();
        TLatex *labelDetSysXSectionPaper= new TLatex(0.95,0.83,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
        labelDetSysXSectionPaper->Draw();

//         legendXsectionPaper->AddEntry(dummyNorm,"Norm. unc. 3.1%","f");

//         legendXsectionPaper->AddEntry(histoDPMJetPi0,"DPMJet","l");
//         legendXsectionPaper->AddEntry(histoHIJINGPi0,"HIJING","l");
//         legendXsectionPaper->AddEntry(histoEPOS3Pi0,"EPOS3","l");
//         legendXsectionPaper->AddEntry(graphMcGillPi0,"Shen #it{et al.}","f");
//         legendXsectionPaper->AddEntry(graphMcGillPi0,"iEBE-VISHNU","f");
        legendXsectionPaper->Draw();


    padInvSectionNLORatio->cd();
    padInvSectionNLORatio->SetLogx(1);
        TH2F * ratio2DEPOS               = new TH2F("ratio2DEPOS","ratio2DEPOS",1000,minPtPi0Plotting, maxPtPi0Plotting,1000,0.3,2.2);
        SetStyleHistoTH2ForGraphs(ratio2DEPOS, "#it{p}_{T} (GeV/#it{c})","#frac{EPOS, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                  0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.257/(textsizeFacXSecMiddle*marginXSec), 510, 508);
        ratio2DEPOS->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DEPOS->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DEPOS->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DEPOS->GetXaxis()->SetLabelFont(42);
        ratio2DEPOS->GetYaxis()->SetLabelFont(42);
        ratio2DEPOS->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DEPOS->GetXaxis()->SetTickLength(0.07);
        ratio2DEPOS->DrawCopy();

        boxErrorNormRatioPi0->Draw();

//         SetStyleHisto(histoRatioPi0EPOS3ToFit, widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );
//         DrawGammaSetMarkerTGraphErr(graphRatioPi0EPOS3ToFit, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
//         histoRatioPi0EPOS3ToFit->Draw("same,hist,l");
//         graphRatioPi0EPOS3ToFit->Draw("same,3");
//         DrawGammaSetMarkerTGraphErr(graphRatioPi0McGillToFit, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
//         graphRatioPi0McGillToFit->Draw("same,3");

        for (Int_t cent= 0; cent < maxCentRun1+maxCentRun2; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
            graphRatioPi0CombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioPi0CombCombFitSys[cent]->SetLineWidth(0);
            graphRatioPi0CombCombFitSys[cent]->Draw("2,same");
            graphRatioPi0CombCombFitStatWOXErr[cent]->Draw("p,same");
        }
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting,1., 1.,0.1,kGray);

    padInvSectionPythiaRatio->cd();
    padInvSectionPythiaRatio->SetLogx(1);
        TH2F * ratio2DGenerator            = new TH2F("ratio2DGenerator","ratio2DGenerator",1000,minPtPi0Plotting, maxPtPi0Plotting,1000,0.3,2.2);
        SetStyleHistoTH2ForGraphs(ratio2DGenerator, "#it{p}_{T} (GeV/#it{c})","#frac{Generator, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                                  0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.257/(textsizeFacXSecDown*marginXSec), 510, 508);
        ratio2DGenerator->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DGenerator->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DGenerator->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DGenerator->GetXaxis()->SetLabelFont(42);
        ratio2DGenerator->GetYaxis()->SetLabelFont(42);
        ratio2DGenerator->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DGenerator->GetXaxis()->SetTickLength(0.06);
        ratio2DGenerator->GetYaxis()->SetTickLength(0.04);
        ratio2DGenerator->DrawCopy();

//         SetStyleHisto(histoRatioPi0DPMJetToFit, widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );
//         histoRatioPi0DPMJetToFit->GetXaxis()->SetRangeUser(0.5,14);
//         histoRatioPi0DPMJetToFit->Draw("same,hist,l");
//         SetStyleHisto(histoRatioPi0HIJINGToFit, widthCommonFit*1.5, styleLineHIJING, colorHIJING );
//         histoRatioPi0HIJINGToFit->GetXaxis()->SetRangeUser(0.5,14);
//         histoRatioPi0HIJINGToFit->Draw("same,hist,l");
        for (Int_t cent= 0; cent < maxCentRun1+maxCentRun2; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
            graphRatioPi0CombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioPi0CombCombFitSys[cent]->SetLineWidth(0);
            graphRatioPi0CombCombFitSys[cent]->Draw("2,same");
            graphRatioPi0CombCombFitStatWOXErr[cent]->Draw("p,same");
        }
        boxErrorNormRatioPi0->Draw();

        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting,1., 1.,0.1,kGray);

    canvasInvSectionPaper->Print(Form("%s/Pi0_InvYieldWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));
//
//     //*************************************************************************************************************
//     //***************************** Paper plot invyield pi0 with hard-physics *************************************
//     //*************************************************************************************************************
//     cout << __LINE__ << endl;
//     padInvSectionSpec->cd();
//     padInvSectionSpec->SetLogy(1);
//     padInvSectionSpec->SetLogx(1);
//         histo2DYieldPi0->DrawCopy();
//
// //         DrawGammaSetMarkerTGraphAsym(graphNLODSS14Pi0, 0, 0, colorDSSBand, colorDSSBand, widthLinesBoxes, kTRUE, colorDSSBand);
// //         graphNLODSS14Pi0->Draw("3,same");
// //         DrawGammaSetMarkerTGraphAsym(graphNLODSS14Pi0PP, 0, 0, colorDSSBand, colorDSSBand, widthLinesBoxes, kTRUE, colorDSSBand);
// //         graphNLODSS14Pi0PP->Draw("3,same");
// //         DrawGammaSetMarkerTGraphAsym(graphNLODSS14nPDFPi0, 0, 0, colorDSSnPDFBand, colorDSSnPDFBand, widthLinesBoxes, kTRUE, colorDSSnPDFBand, kTRUE);
// //         graphNLODSS14nPDFPi0->Draw("3,same");
// //         DrawGammaSetMarkerTGraphAsym(graphNLODSS14nPDFEPPS16Pi0, 0, 0, colorDSSnPDFEPPSBand, colorDSSnPDFEPPSBand, widthLinesBoxes, kTRUE, colorDSSnPDFEPPSBand, kTRUE);
// //         graphNLODSS14nPDFEPPS16Pi0->Draw("3,same");
// //
// //         DrawGammaNLOTGraph( graphNLODSS14Pi0Center, widthCommonFit, styleLineDSS, colorDSS);
// //         graphNLODSS14Pi0Center->Draw("same,c");
// //         DrawGammaNLOTGraph( graphPi0CGC5TeV, widthCommonFit, styleLineCGC, colorCGC);
// //         graphPi0CGC5TeV->Draw("same,c");
// //         DrawGammaNLOTGraph( graphNLODSS14Pi0PPCenter, widthCommonFit, styleLineDSS, colorDSS);
// //         graphNLODSS14Pi0PPCenter->Draw("same,c");
// //         DrawGammaNLOTGraph( graphNLODSS14nPDFEPPS16Pi0Center, widthCommonFit, styleLineDSSnPDFEPPS, colorDSSnPDFEPPS);
// //         graphNLODSS14nPDFEPPS16Pi0Center->Draw("same,c");
// //         DrawGammaNLOTGraph( graphNLODSS14nPDFPi0Center, widthCommonFit, styleLineDSSnPDF, colorDSSnPDF);
// //         graphNLODSS14nPDFPi0Center->Draw("same,c");
//
//         DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
//         graphCombPi0InvYieldSys[cent]->Draw("E2same");
//         DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombPi0Sys, markerStyleComb+4, markerSizeComb+0.2, kGray+1, kGray+1, widthLinesBoxes, kTRUE);
//         graphPPInvYieldCombPi0Sys->Draw("E2same");
//         DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent]);
//         graphCombPi0InvYieldStatWOXErr[cent]->Draw("p,same,z");
//
//         TGraphAsymmErrors* graphPPInvYieldCombPi0StatWOXErr = (TGraphAsymmErrors*)graphPPInvYieldCombPi0Stat->Clone("graphPPInvYieldCombPi0StatWOXErr");
//         ProduceGraphAsymmWithoutXErrors(graphPPInvYieldCombPi0StatWOXErr);
//
//         DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombPi0StatWOXErr, markerStyleComb+4, markerSizeComb+0.2, kGray+1, kGray+1);
//         graphPPInvYieldCombPi0StatWOXErr->Draw("p,same,z");
//
// //         DrawGammaSetMarkerTF1( fitPowInvYieldPi0NLOpPbnCTEQ, 7, 4, kOrange);
// //         fitPowInvYieldPi0NLOpPbnCTEQ->Draw("same");
//
//         DrawGammaSetMarkerTF1( fitTCMInvYieldPi0[cent], 7, 2, kGray+2);
//         fitTCMInvYieldPi0[cent]->Draw("same");
//         DrawGammaSetMarkerTF1( fitPPTCMInvYieldPi0, 4, 2, kGray+1);
//         fitPPTCMInvYieldPi0->Draw("same");
//
//         TLatex *labelALICENLOPaper= new TLatex(0.19,0.92,textALICE.Data());
//         SetStyleTLatex( labelALICENLOPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 11);
//         labelALICENLOPaper->Draw();
//
//         TLegend* legendPaperPi0WithNLO    = GetAndSetLegend2(0.19, 0.03, 0.5, 0.03+0.05*8, textSizeLabelsPixel, 1, collisionSystempPb.Data(), 43, 0.2);
// //         legendPaperPi0WithNLO->AddEntry(graphCombPi0InvYieldSys[cent],Form("%s %s","#pi^{0}", collisionSystempPb.Data()),"pf");
//         legendPaperPi0WithNLO->AddEntry(graphCombPi0InvYieldSys[cent],"#pi^{0}","pf");
//         legendPaperPi0WithNLO->AddEntry(dummyNorm,"Norm. unc. 3.1%","f");
//         legendPaperPi0WithNLO->AddEntry(fitTCMInvYieldPi0[cent],"TCM fit","l");
// //         legendPaperPi0WithNLO->AddEntry(graphPi0CGC5TeV,"CGC","l");
// //         legendPaperPi0WithNLO->AddEntry(graphNLODSS14Pi0,"NLO, PDF: CT10, FF: DSS14","f");
// //         legendPaperPi0WithNLO->AddEntry(graphNLODSS14nPDFPi0,"NLO, nPDF: nCTEQ, FF: DSS14","f");
// //         legendPaperPi0WithNLO->AddEntry(graphNLODSS14nPDFEPPS16Pi0,"NLO, nPDF: EPPS16, FF: DSS14","f");
//         legendPaperPi0WithNLO->Draw();
//
//         DrawGammaLines(0.30, 0.385, 9.85e-8, 9.85e-8, widthCommonFit, colorDSS, styleLineDSS);
//         DrawGammaLines(0.30, 0.385, 2.85e-8, 2.85e-8, widthCommonFit, colorDSSnPDF, styleLineDSSnPDF);
//         DrawGammaLines(0.30, 0.385, 8.0e-9, 8.0e-9, widthCommonFit, colorDSSnPDFEPPS, styleLineDSSnPDFEPPS);
//
//         TGraphErrors* dummyNorm2     = new TGraphErrors(1);
//         dummyNorm2->SetPoint(0,1,0.0);
//         DrawGammaSetMarkerTGraphErr(dummyNorm2, 0, 0, kGray+2, kGray+2, widthLinesBoxes, kTRUE, kGray+2);
//
//         TLegend* legendPaperPi0WithNLOPP    = GetAndSetLegend2(0.60, 0.955-0.05*4, 0.95, 0.955, textSizeLabelsPixel, 1, collisionSystempp.Data(), 43, 0.177);
//         legendPaperPi0WithNLOPP->SetTextAlign(12);
//         legendPaperPi0WithNLOPP->AddEntry(graphPPInvYieldCombPi0Sys,"#pi^{0} interpolated","pf");
// //         legendPaperPi0WithNLOPP->AddEntry(graphPPInvYieldCombPi0Sys,Form("%s %s","#pi^{0}", collisionSystempp.Data()),"pf");
//         legendPaperPi0WithNLOPP->AddEntry(dummyNorm2,"Norm. unc. 0.89%","f");
//         legendPaperPi0WithNLOPP->AddEntry(fitPPTCMInvYieldPi0,"TCM fit","l");
//         legendPaperPi0WithNLOPP->Draw();
//
//     padInvSectionNLORatio->cd();
//     padInvSectionNLORatio->SetLogx(1);
//     TH2F * ratio2DNLOpPb               = new TH2F("ratio2DNLOpPb","ratio2DNLOpPb",1000,minPtPi0Plotting, 25.,1000,0.3,2.2);
//     SetStyleHistoTH2ForGraphs(ratio2DNLOpPb, "#it{p}_{T} (GeV/#it{c})","#frac{NLO, p-Pb Data}{p-Pb fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
//                               0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.257/(textsizeFacXSecMiddle*marginXSec), 510, 508);
//     //         ratio2DNLOpPb->GetYaxis()->SetNdivisions(512);
//     ratio2DNLOpPb->GetYaxis()->SetNoExponent(kTRUE);
//     ratio2DNLOpPb->GetXaxis()->SetMoreLogLabels(kTRUE);
//     ratio2DNLOpPb->GetXaxis()->SetNoExponent(kTRUE);
//     ratio2DNLOpPb->GetXaxis()->SetLabelFont(42);
//     ratio2DNLOpPb->GetYaxis()->SetLabelFont(42);
//     ratio2DNLOpPb->GetYaxis()->SetLabelOffset(+0.01);
//     ratio2DNLOpPb->GetXaxis()->SetTickLength(0.07);
//     ratio2DNLOpPb->DrawCopy();
//
//         boxErrorNormRatioPi0->Draw();
//
// //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0NLODSS14, 0, 0, colorDSSBand, colorDSSBand, widthLinesBoxes, kTRUE, colorDSSBand);
// //         graphRatioPi0NLODSS14->Draw("3,same");
// //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0NLODSS14nPDF, 0, 0, colorDSSnPDFBand, colorDSSnPDFBand, widthLinesBoxes, kTRUE, colorDSSnPDFBand, kTRUE);
// //         graphRatioPi0NLODSS14nPDF->Draw("3,same");
// //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0NLODSS14nPDFEPPS16, 0, 0, colorDSSnPDFEPPSBand, colorDSSnPDFEPPSBand, widthLinesBoxes, kTRUE, colorDSSnPDFEPPSBand, kTRUE);
// //         graphRatioPi0NLODSS14nPDFEPPS16->Draw("3,same");
//
// //         DrawGammaNLOTGraph( graphRatioPi0NLODSS14Center, widthCommonFit, styleLineDSS, colorDSS);
// //         graphRatioPi0NLODSS14Center->Draw("same,c");
// //         DrawGammaNLOTGraph( graphRatioPi0CGC, widthCommonFit, styleLineCGC, colorCGC);
// //         DrawGammaNLOTGraph( graphRatioPi0NLODSS14nPDFCenter, widthCommonFit, styleLineDSSnPDF, colorDSSnPDF);
// //         DrawGammaNLOTGraph( graphRatioPi0NLODSS14nPDFEPPS16Center, widthCommonFit, styleLineDSSnPDFEPPS, colorDSSnPDFEPPS);
//
//
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
//         graphRatioPi0CombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
//         graphRatioPi0CombCombFitSys[cent]->SetLineWidth(0);
//         graphRatioPi0CombCombFitSys[cent]->Draw("2,same");
//         graphRatioPi0CombCombFitStatWOXErr[cent]->Draw("p,same");
//
// //         graphRatioPi0CGC->Draw("same,c");
// //         graphRatioPi0NLODSS14nPDFEPPS16Center->Draw("same,c");
// //         graphRatioPi0NLODSS14nPDFCenter->Draw("same,c");
//
//         TLatex *labelEnergyRatioPaperNLOpPb= new TLatex(0.19, 1-textsizeLabelsXSecMiddle*1.5, collisionSystempPb.Data());
//         SetStyleTLatex( labelEnergyRatioPaperNLOpPb, textsizeLabelsXSecMiddle,4, 1, 42, kTRUE, 11);
//         labelEnergyRatioPaperNLOpPb->Draw();
//
//         ratio2DNLOpPb->Draw("same,axis");
//         DrawGammaLines(minPtPi0Plotting, 25.,1., 1.,0.1,kGray);
//
//     cout << __LINE__ << endl;
//
//     padInvSectionPythiaRatio->cd();
//     padInvSectionPythiaRatio->SetLogx(1);
//     TH2F * ratio2DNLOpp            = new TH2F("ratio2DNLOpp","ratio2DNLOpp",1000,minPtPi0Plotting, 25.,1000,0.3,2.2);
//     SetStyleHistoTH2ForGraphs(ratio2DNLOpp, "#it{p}_{T} (GeV/#it{c})","#frac{NLO, pp Data}{pp fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
//                               0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.257/(textsizeFacXSecDown*marginXSec), 510, 508);
//     //         ratio2DNLOpp->GetYaxis()->SetNdivisions(512);
//     ratio2DNLOpp->GetYaxis()->SetNoExponent(kTRUE);
//     ratio2DNLOpp->GetXaxis()->SetMoreLogLabels(kTRUE);
//     ratio2DNLOpp->GetXaxis()->SetNoExponent(kTRUE);
//     ratio2DNLOpp->GetXaxis()->SetLabelFont(42);
//     ratio2DNLOpp->GetYaxis()->SetLabelFont(42);
//     ratio2DNLOpp->GetYaxis()->SetLabelOffset(+0.01);
//     ratio2DNLOpp->GetXaxis()->SetTickLength(0.06);
//     ratio2DNLOpp->GetYaxis()->SetTickLength(0.04);
//     ratio2DNLOpp->DrawCopy();
//
//         TBox* boxErrorNormRatioPi0PP          = CreateBoxConv(kGray+2, 0.28, 1.-(xSection5TeVErr/xSection5TeV ), 0.32, 1.+(xSection5TeVErr/xSection5TeV));
//         boxErrorNormRatioPi0PP->Draw();
//
// //         DrawGammaSetMarkerTGraphAsym(graphRatioPPPi0NLODSS14, 0, 0, colorDSSBand, colorDSSBand, widthLinesBoxes, kTRUE, colorDSSBand);
// //         graphRatioPPPi0NLODSS14->Draw("3,same");
// //         DrawGammaNLOTGraph( graphRatioPPPi0NLODSS14Center, widthCommonFit, styleLineDSS, colorDSS);
// //         graphRatioPPPi0NLODSS14Center->Draw("same,c");
//
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombToFitStatPP, markerStyleComb+4, markerSizeComb+0.2, kGray+2, kGray+2, widthLinesBoxes, kFALSE);
//         graphRatioPi0CombToFitStatPP->SetLineWidth(widthLinesBoxes);
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombToFitSysPP, markerStyleComb+4, markerSizeComb+0.2, kGray+2, kGray+2, widthLinesBoxes, kTRUE, 0);
//         graphRatioPi0CombToFitSysPP->SetLineWidth(0);
//         graphRatioPi0CombToFitSysPP->Draw("2,same");
//         graphRatioPi0CombToFitStatPP->Draw("p,same");
//
//         TLatex *labelEnergyRatioPaperNLOPP= new TLatex(0.19, 1-textsizeLabelsXSecDown*1.5, collisionSystempp.Data());
//         SetStyleTLatex( labelEnergyRatioPaperNLOPP, textsizeLabelsXSecDown,4, 1, 42, kTRUE, 11);
//         labelEnergyRatioPaperNLOPP->Draw();
//
//         ratio2DNLOpp->Draw("same,axis");
//         DrawGammaLines(minPtPi0Plotting, 25.,1., 1.,0.1,kGray);
//
//     canvasInvSectionPaper->Print(Form("%s/Pi0_InvYieldWithRatios_Paper_NLOs.%s",outputDir.Data(),suffix.Data()));
//
//     cout << __LINE__ << endl;
//
//
//
//     //*****************************************************************************************************************
//     // Plotting ratio to different fitting functions for pi0
//     //*****************************************************************************************************************
//     textSizeLabelsPixel                         = 52;
//     TCanvas* canvasRatioDiffFits                = new TCanvas("canvasRatioDiffFits","",200,10,1350,900);  // gives the page size
//     DrawGammaCanvasSettings( canvasRatioDiffFits,  0.082, 0.01, 0.01, 0.12);
//     canvasRatioDiffFits->cd();
//     canvasRatioDiffFits->SetLogx();
//     Double_t textsizeLabelsFits                 = 0;
//     if (canvasRatioDiffFits->XtoPixel(canvasRatioDiffFits->GetX2()) <canvasRatioDiffFits->YtoPixel(canvasRatioDiffFits->GetY1()) ){
//         textsizeLabelsFits                      = (Double_t)textSizeLabelsPixel/canvasRatioDiffFits->XtoPixel(canvasRatioDiffFits->GetX2()) ;
//     } else {
//         textsizeLabelsFits                      = (Double_t)textSizeLabelsPixel/canvasRatioDiffFits->YtoPixel(canvasRatioDiffFits->GetY1());
//     }
//     cout << textsizeLabelsFits << endl;
//
//     // calculating ratios to fits
//     TGraphAsymmErrors* graphRatioPi0TCMTot      = (TGraphAsymmErrors*)graphCombPi0InvYieldTot[cent]->Clone();
//     graphRatioPi0TCMTot                         = CalculateGraphErrRatioToFit(graphRatioPi0TCMTot, fitTCMInvYieldPi0[cent]);
//     TGraphAsymmErrors* graphRatioPi0TsallisTot  = (TGraphAsymmErrors*)graphCombPi0InvYieldTot[cent]->Clone();
//     graphRatioPi0TsallisTot                     = CalculateGraphErrRatioToFit(graphRatioPi0TsallisTot, fitInvYieldPi0[cent]);
//     TGraphAsymmErrors* graphRatioPi0PowTot      = (TGraphAsymmErrors*)graphCombPi0InvYieldTot[cent]->Clone();
//     graphRatioPi0PowTot                         = CalculateGraphErrRatioToFit(graphRatioPi0PowTot, fitPowInvYieldPi0Tot);
//     TGraphAsymmErrors* graphRatioPi0OHagTot     = (TGraphAsymmErrors*)graphCombPi0InvYieldTot[cent]->Clone();
//     graphRatioPi0OHagTot                        = CalculateGraphErrRatioToFit(graphRatioPi0OHagTot, fitOHagInvYieldPi0Tot);
//     // restricting pT range of power law fit drawing
//     while (graphRatioPi0PowTot->GetX()[0] < 2.) graphRatioPi0PowTot->RemovePoint(0);
//
//
//     TH2F * ratio2DFits       = new TH2F("ratio2DFits","ratio2DFits",1000,minPtPi0Plotting, 25.,1000,0.2,2.3);
//     SetStyleHistoTH2ForGraphs(ratio2DFits, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsFits, textsizeLabelsFits,
//                               0.85*textsizeLabelsFits,textsizeLabelsFits, 0.88, 0.7, 510, 505);
//     ratio2DFits->GetYaxis()->SetMoreLogLabels(kTRUE);
//     ratio2DFits->GetYaxis()->SetNdivisions(505);
//     ratio2DFits->GetYaxis()->SetNoExponent(kTRUE);
//     ratio2DFits->GetXaxis()->SetMoreLogLabels(kTRUE);
//     ratio2DFits->GetXaxis()->SetNoExponent(kTRUE);
//     ratio2DFits->GetXaxis()->SetLabelFont(42);
//     ratio2DFits->GetYaxis()->SetLabelFont(42);
//     ratio2DFits->DrawCopy();
//
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0TsallisTot, 25, markerSizeComb, kGray+1, kGray+1, widthLinesBoxes, kFALSE);
//         graphRatioPi0TsallisTot->Draw("p,same");
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0PowTot, 33, markerSizeComb*2, kAzure+2, kAzure+2, widthLinesBoxes, kFALSE);
//         graphRatioPi0PowTot->Draw("p,same");
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0TCMTot, markerStyleComb, markerSizeComb, kGray+2, kGray+2, widthLinesBoxes, kFALSE);
//         graphRatioPi0TCMTot->Draw("p,same");
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0OHagTot, 24, markerSizeComb+0.2, kRed+2, kRed+2, widthLinesBoxes, kFALSE);
//         graphRatioPi0OHagTot->Draw("p,same");
//
//         DrawGammaLines(minPtPi0Plotting, 25. , 1.1, 1.1,1, kGray, 3);
//         DrawGammaLines(minPtPi0Plotting, 25. , 1, 1,1, kGray, 7);
//         DrawGammaLines(minPtPi0Plotting, 25. , 0.9, 0.9,1, kGray, 3);
//
//         TLegend* legendRatioPi0Fits= GetAndSetLegend2(0.12,0.95-4*1.05*textsizeLabelsFits,0.37,0.95, textSizeLabelsPixel, 1, "", 43, 0.2);
//         legendRatioPi0Fits->AddEntry(graphRatioPi0TCMTot,"TCM","p");
//         legendRatioPi0Fits->AddEntry(graphRatioPi0TsallisTot,"Levy-Tsallis","p");
//         legendRatioPi0Fits->AddEntry(graphRatioPi0OHagTot,"mod. Hagedorn","p");
//         legendRatioPi0Fits->AddEntry(graphRatioPi0PowTot,"pure powerlaw, 4-20 GeV/#it{c}","p");
//         legendRatioPi0Fits->Draw();
//
//         TLatex *labelRatioDiffFitsEnergy    = new TLatex(0.95,0.91,collisionSystempPb.Data());
//         SetStyleTLatex( labelRatioDiffFitsEnergy, textsizeLabelsFits,4,1, 42, kTRUE, 31);
//         labelRatioDiffFitsEnergy->Draw();
//         TLatex *labelRatioFitsPi0           = new TLatex(0.95,0.91-1.05*textsizeLabelsFits,"#pi^{0} #rightarrow #gamma#gamma");
//         SetStyleTLatex( labelRatioFitsPi0, textsizeLabelsFits,4, 1, 42, kTRUE, 31);
//         labelRatioFitsPi0->Draw();
//         TLatex *labelRatioFitsThesis        = new TLatex(0.13,0.16,textALICE);
//         SetStyleTLatex( labelRatioFitsThesis, textsizeLabelsFits,4, 1, 42, kTRUE, 11);
//         labelRatioFitsThesis->Draw();
//
//
//         ratio2DFits->Draw("axis,same");
//
//     canvasRatioDiffFits->Update();
//     canvasRatioDiffFits->Print(Form("%s/Pi0_RatioDiffFits.%s",outputDir.Data(),suffix.Data()));

    //*************************************************************************************************************
    //***************************** Paper plot invyield eta with soft-physics *************************************
    //*************************************************************************************************************
    TLegend* legendXsectionPaperEta    = GetAndSetLegend2(0.19, 0.03, 0.7, 0.03+0.05*(maxCentRun1+maxCentRun2), textSizeLabelsPixel, 2, "", 43, 0.25);

    canvasInvSectionPaper->cd();
    padInvSectionSpec->cd();
    padInvSectionSpec->SetLogy(1);
    padInvSectionSpec->SetLogx(1);
        SetStyleHistoTH2ForGraphs(histo2DYieldEta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                                0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.257/(textsizeFacXSecUp*marginXSec));
        histo2DYieldEta->GetXaxis()->SetMoreLogLabels();
        histo2DYieldEta->GetXaxis()->SetLabelOffset(+0.01);
        histo2DYieldEta->Draw();

//         SetStyleHisto(histoEPOS3Eta, widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );
//         DrawGammaSetMarkerTGraphErr(graphEPOS3Eta, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
//         histoEPOS3Eta->GetXaxis()->SetRangeUser(0.3,20);
//         histoEPOS3Eta->Draw("same,hist,l");
//         while(graphEPOS3Eta->GetX()[0] < 0.3)
//             graphEPOS3Eta->RemovePoint(0);
//         while(graphEPOS3Eta->GetX()[graphEPOS3Eta->GetN()-1] > 20)
//             graphEPOS3Eta->RemovePoint(graphEPOS3Eta->GetN()-1);
//         graphEPOS3Eta->Draw("same,3");

//         DrawGammaSetMarkerTGraphErr(graphMcGillEta, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
//         while(graphMcGillEta->GetX()[0] < 0.3)
//             graphMcGillEta->RemovePoint(0);
//         graphMcGillEta->Draw("same,3");

//         SetStyleHisto(histoDPMJetEta, widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );
//         histoDPMJetEta->GetXaxis()->SetRangeUser(0.3,20);
//         histoDPMJetEta->Draw("same,hist,l");
//         SetStyleHisto(histoHIJINGEta, widthCommonFit*1.5, styleLineHIJING, colorHIJING );
//         histoHIJINGEta->GetXaxis()->SetRangeUser(0.3,20);
//         histoHIJINGEta->Draw("same,hist,l");

        for (Int_t cent= 0; cent < maxCentRun1+maxCentRun2; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
            graphCombEtaInvYieldSys[cent]->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent]);
            graphCombEtaInvYieldStatWOXErr[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTF1( fitTCMInvYieldEta[cent], 7, 2, colorCentMC[cent]);
            fitTCMInvYieldEta[cent]->Draw("same");
            legendXsectionPaperEta->AddEntry(graphCombEtaInvYieldSys[cent],centArray[cent].Data(),"pf");
            legendXsectionPaperEta->AddEntry(fitTCMInvYieldEta[cent],"TCM fit","l");
        }

        labelEnergyXSectionPaper->Draw();
        labelALICEXSectionPaper->Draw();
        TLatex *labelDetSysXSectionPaperEta= new TLatex(0.95,0.83,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPaperEta, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
        labelDetSysXSectionPaperEta->Draw();

//         legendXsectionPaper->AddEntry(dummyNorm,"Norm. unc. 3.1%","f");

//         legendXsectionPaper->AddEntry(histoDPMJetEta,"DPMJet","l");
//         legendXsectionPaper->AddEntry(histoHIJINGEta,"HIJING","l");
//         legendXsectionPaper->AddEntry(histoEPOS3Eta,"EPOS3","l");
//         legendXsectionPaper->AddEntry(graphMcGillEta,"Shen #it{et al.}","f");
//         legendXsectionPaper->AddEntry(graphMcGillEta,"iEBE-VISHNU","f");
        legendXsectionPaper->Draw();


    padInvSectionNLORatio->cd();
    padInvSectionNLORatio->SetLogx(1);
        TH2F * ratio2DEPOSEta               = new TH2F("ratio2DEPOSEta","ratio2DEPOSEta",1000, minPtEtaPlotting, maxPtEtaPlotting,1000,0.3,2.2);
        SetStyleHistoTH2ForGraphs(ratio2DEPOSEta, "#it{p}_{T} (GeV/#it{c})","#frac{EPOS, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                  0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.257/(textsizeFacXSecMiddle*marginXSec), 510, 508);
        ratio2DEPOSEta->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DEPOSEta->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DEPOSEta->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DEPOSEta->GetXaxis()->SetLabelFont(42);
        ratio2DEPOSEta->GetYaxis()->SetLabelFont(42);
        ratio2DEPOSEta->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DEPOSEta->GetXaxis()->SetTickLength(0.07);
        ratio2DEPOSEta->DrawCopy();

        boxErrorNormRatioEta->Draw();

//         SetStyleHisto(histoRatioEtaEPOS3ToFit, widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );
//         DrawGammaSetMarkerTGraphErr(graphRatioEtaEPOS3ToFit, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
//         histoRatioEtaEPOS3ToFit->Draw("same,hist,l");
//         graphRatioEtaEPOS3ToFit->Draw("same,3");
//         DrawGammaSetMarkerTGraphErr(graphRatioEtaMcGillToFit, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
//         graphRatioEtaMcGillToFit->Draw("same,3");

        for (Int_t cent= 0; cent < maxCentRun1+maxCentRun2; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
            graphRatioEtaCombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioEtaCombCombFitSys[cent]->SetLineWidth(0);
            graphRatioEtaCombCombFitSys[cent]->Draw("2,same");
            graphRatioEtaCombCombFitStatWOXErr[cent]->Draw("p,same");
        }
        DrawGammaLines(minPtEtaPlotting, 31.,1., 1.,0.1,kGray);

    padInvSectionPythiaRatio->cd();
    padInvSectionPythiaRatio->SetLogx(1);
        TH2F * ratio2DGeneratorEta            = new TH2F("ratio2DGeneratorEta","ratio2DGeneratorEta",1000, minPtEtaPlotting, maxPtEtaPlotting,1000,0.3,2.2);
        SetStyleHistoTH2ForGraphs(ratio2DGeneratorEta, "#it{p}_{T} (GeV/#it{c})","#frac{Generator, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                                  0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.257/(textsizeFacXSecDown*marginXSec), 510, 508);
        ratio2DGeneratorEta->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DGeneratorEta->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DGeneratorEta->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DGeneratorEta->GetXaxis()->SetLabelFont(42);
        ratio2DGeneratorEta->GetYaxis()->SetLabelFont(42);
        ratio2DGeneratorEta->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DGeneratorEta->GetXaxis()->SetTickLength(0.06);
        ratio2DGeneratorEta->GetYaxis()->SetTickLength(0.04);
        ratio2DGeneratorEta->DrawCopy();

//         SetStyleHisto(histoRatioEtaDPMJetToFit, widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );
//         histoRatioEtaDPMJetToFit->GetXaxis()->SetRangeUser(0.5,14);
//         histoRatioEtaDPMJetToFit->Draw("same,hist,l");
//         SetStyleHisto(histoRatioEtaHIJINGToFit, widthCommonFit*1.5, styleLineHIJING, colorHIJING );
//         histoRatioEtaHIJINGToFit->GetXaxis()->SetRangeUser(0.5,14);
//         histoRatioEtaHIJINGToFit->Draw("same,hist,l");
        for (Int_t cent= 0; cent < maxCentRun1+maxCentRun2; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
            graphRatioEtaCombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioEtaCombCombFitSys[cent]->SetLineWidth(0);
            graphRatioEtaCombCombFitSys[cent]->Draw("2,same");
            graphRatioEtaCombCombFitStatWOXErr[cent]->Draw("p,same");
        }
        boxErrorNormRatioEta->Draw();

        DrawGammaLines(minPtEtaPlotting, 31.,1., 1.,0.1,kGray);

    canvasInvSectionPaper->Print(Form("%s/Eta_InvYieldWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));



    // **********************************************************************************************************************
    // ************************* Saving of final results ********************************************************************
    // **********************************************************************************************************************

    TString nameOutputCommonFile    = Form("CombinedResultsPaperpPb5TeVCent_%s.root", dateForOutput.Data());

    TFile fCombResults(nameOutputCommonFile.Data(), "UPDATE");

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCent[cent]) continue;
        TDirectoryFile* directoryPi0Output  = NULL;
        directoryPi0Output                  = (TDirectoryFile*)fCombResults.Get(Form("Pi0pPb5TeV_%s%s", centArray[cent].Data(), runArray[cent].Data()));
        if (!directoryPi0Output){
            fCombResults.mkdir(Form("Pi0pPb5TeV_%s%s", centArray[cent].Data(), runArray[cent].Data()));
            directoryPi0Output              = (TDirectoryFile*)fCombResults.Get(Form("Pi0pPb5TeV_%s%s", centArray[cent].Data(), runArray[cent].Data()));
        }
        fCombResults.cd(Form("Pi0pPb5TeV_%s%s", centArray[cent].Data(), runArray[cent].Data()));

            // Final spectrum
            if (graphCombPi0InvYieldTot[cent]) graphCombPi0InvYieldTot[cent]->Write("graphInvYieldPi0CombTotErr",TObject::kOverwrite);
            if (graphCombPi0InvYieldStat[cent]) graphCombPi0InvYieldStat[cent]->Write("graphInvYieldPi0CombStatErr",TObject::kOverwrite);
            if (graphCombPi0InvYieldSys[cent]) graphCombPi0InvYieldSys[cent]->Write("graphInvYieldPi0CombSysErr",TObject::kOverwrite);
             // fits for pi0
            if (fitInvYieldPi0[cent]) fitInvYieldPi0[cent]->Write("TsallisFitPi0");
            if (fitTCMInvYieldPi0[cent]) fitTCMInvYieldPi0[cent]->Write("TwoComponentModelFitPi0",TObject::kOverwrite);

            // Final inv yield INEL
            if (bWCorrection.Contains("Y")){
                // Final spectrum correlations Method A
                if(graphCombPi0InvYieldTot_yShifted[cent])graphCombPi0InvYieldTot_yShifted[cent]->Write("graphInvYieldPi0CombTotErr_yShifted",TObject::kOverwrite);
                if(graphCombPi0InvYieldStat_yShifted[cent])graphCombPi0InvYieldStat_yShifted[cent]->Write("graphInvYieldPi0CombStatErr_yShifted",TObject::kOverwrite);
                if(graphCombPi0InvYieldSys_yShifted[cent])graphCombPi0InvYieldSys_yShifted[cent]->Write("graphInvYieldPi0CombSysErr_yShifted",TObject::kOverwrite);
            }
            // writing individual measurements
            for (Int_t meth = 0; meth< 11; meth++){
                if (graphIndPi0InvYieldStat[cent]){
                    if (graphIndPi0InvYieldStat[cent][meth]) graphIndPi0InvYieldStat[cent][meth]->Write(Form("graphInvYieldPi0%sStatErr", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                    if (graphIndPi0InvYieldSys[cent][meth]) graphIndPi0InvYieldSys[cent][meth]->Write(Form("graphInvYieldPi0%sSysErr", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                    if (bWCorrection.Contains("Y")){
                        if (graphIndPi0InvYieldStat_yShifted[cent][meth]) graphIndPi0InvYieldStat_yShifted[cent][meth]->Write(Form("graphInvYieldPi0%sStatErr_yShifted", nameMeasGlobalLabel[meth].Data()),
                                                                                                                            TObject::kOverwrite);
                        if (graphIndPi0InvYieldSys_yShifted[cent][meth]) graphIndPi0InvYieldSys_yShifted[cent][meth]->Write(Form("graphInvYieldPi0%sSysErr_yShifted", nameMeasGlobalLabel[meth].Data()),
                                                                                                                            TObject::kOverwrite);
                    }
                }
            }

            TDirectoryFile* supporting  = NULL;
            supporting                  = (TDirectoryFile*)directoryPi0Output->Get("Supporting");
            if (!supporting)
                directoryPi0Output->mkdir("Supporting");
            directoryPi0Output->cd("Supporting");
                // Writing full correction factors
                for (Int_t meth = 0; meth< 11; meth++){
                    if (graphPi0EffTimesAcc[cent][meth]) graphPi0EffTimesAcc[cent][meth]->Write(Form("Pi0CorrectionFactor%s", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                    if (graphPi0Mass[cent][meth]) graphPi0Mass[cent][meth]->Write(Form("Pi0MassData%s", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                    if (graphPi0MassMC[cent][meth]) graphPi0MassMC[cent][meth]->Write(Form("Pi0MassMC%s", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                    if (graphPi0Width[cent][meth]) graphPi0Width[cent][meth]->Write(Form("Pi0WidthData%s", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                    if (graphPi0WidthMC[cent][meth]) graphPi0WidthMC[cent][meth]->Write(Form("Pi0WidthMC%s", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                }

        TDirectoryFile* directoryEtaOutput  = NULL;
        directoryEtaOutput                  = (TDirectoryFile*)fCombResults.Get(Form("EtapPb5TeV_%s%s", centArray[cent].Data(), runArray[cent].Data()));
        if (!directoryEtaOutput){
            fCombResults.mkdir(Form("EtapPb5TeV_%s%s", centArray[cent].Data(), runArray[cent].Data()));
            directoryEtaOutput              = (TDirectoryFile*)fCombResults.Get(Form("EtapPb5TeV_%s%s", centArray[cent].Data(), runArray[cent].Data()));
        }
        fCombResults.cd(Form("EtapPb5TeV_%s%s", centArray[cent].Data(), runArray[cent].Data()));

            // Final spectrum
            if (graphCombEtaInvYieldTot[cent])  graphCombEtaInvYieldTot[cent]->Write("graphInvYieldEtaCombTotErr",TObject::kOverwrite);
            if (graphCombEtaInvYieldStat[cent]) graphCombEtaInvYieldStat[cent]->Write("graphInvYieldEtaCombStatErr",TObject::kOverwrite);
            if (graphCombEtaInvYieldSys[cent])  graphCombEtaInvYieldSys[cent]->Write("graphInvYieldEtaCombSysErr",TObject::kOverwrite);
                // fits for pi0
            if (fitInvYieldEta[cent])           fitInvYieldEta[cent]->Write("TsallisFitEta");
            if (fitTCMInvYieldEta[cent])        fitTCMInvYieldEta[cent]->Write("TwoComponentModelFitEta",TObject::kOverwrite);

            // Final inv yield INEL
            if (bWCorrection.Contains("Y")){
                // Final spectrum correlations Method A
                if(graphCombEtaInvYieldTot_yShifted[cent])graphCombEtaInvYieldTot_yShifted[cent]->Write("graphInvYieldEtaCombTotErr_yShifted",TObject::kOverwrite);
                if(graphCombEtaInvYieldStat_yShifted[cent])graphCombEtaInvYieldStat_yShifted[cent]->Write("graphInvYieldEtaCombStatErr_yShifted",TObject::kOverwrite);
                if(graphCombEtaInvYieldSys_yShifted[cent])graphCombEtaInvYieldSys_yShifted[cent]->Write("graphInvYieldEtaCombSysErr_yShifted",TObject::kOverwrite);
            }
            // writing individual measurements
            for (Int_t meth = 0; meth< 11; meth++){
                if (graphIndEtaInvYieldStat[cent][meth]) graphIndEtaInvYieldStat[cent][meth]->Write(Form("graphInvYieldEta%sStatErr", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                if (graphIndEtaInvYieldSys[cent][meth]) graphIndEtaInvYieldSys[cent][meth]->Write(Form("graphInvYieldEta%sSysErr", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                if (bWCorrection.Contains("Y")){
                    if (graphIndEtaInvYieldStat_yShifted[cent][meth]) graphIndEtaInvYieldStat_yShifted[cent][meth]->Write(Form("graphInvYieldEta%sStatErr_yShifted", nameMeasGlobalLabel[meth].Data()),
                        TObject::kOverwrite);
                    if (graphIndEtaInvYieldSys_yShifted[cent][meth]) graphIndEtaInvYieldSys_yShifted[cent][meth]->Write(Form("graphInvYieldEta%sSysErr_yShifted", nameMeasGlobalLabel[meth].Data()),
                        TObject::kOverwrite);
                }
                if (graphEtaToPi0Stat[cent][meth]) graphEtaToPi0Stat[cent][meth]->Write(Form("graphEtaToPi0%sStatErr", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                if (graphEtaToPi0Sys[cent][meth]) graphEtaToPi0Sys[cent][meth]->Write(Form("graphEtaToPi0%sSysErr", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
            }

            TDirectoryFile* supportingEta  = NULL;
            supportingEta                  = (TDirectoryFile*)directoryEtaOutput->Get("Supporting");
            if (!supportingEta)
                directoryEtaOutput->mkdir("Supporting");
            directoryEtaOutput->cd("Supporting");
                // Writing full correction factors
                for (Int_t meth = 0; meth< 11; meth++){
                    if (graphEtaEffTimesAcc[cent][meth]) graphEtaEffTimesAcc[cent][meth]->Write(Form("EtaCorrectionFactor%s", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                    if (graphEtaMass[cent][meth]) graphEtaMass[cent][meth]->Write(Form("EtaMassData%s", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                    if (graphEtaMassMC[cent][meth]) graphEtaMassMC[cent][meth]->Write(Form("EtaMassMC%s", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                    if (graphEtaWidth[cent][meth]) graphEtaWidth[cent][meth]->Write(Form("EtaWidthData%s", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                    if (graphEtaWidthMC[cent][meth]) graphEtaWidthMC[cent][meth]->Write(Form("EtaWidthMC%s", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                }
    }
    fCombResults.Close();


    //  **********************************************************************************************************************
    //  ************************* Saving only fits to final results **********************************************************
    //  **********************************************************************************************************************

    TString nameOutputCommonFileFitsOnly    = Form("FitsPaperpPb5TeVCent_%s.root", dateForOutput.Data());
    TFile fFitsResults(nameOutputCommonFileFitsOnly.Data(), "RECREATE");

        for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
            if (!enableCentComb[cent]) continue;
            // fits for pi0
            fitInvYieldPi0[cent]->Write(Form("TsallisFitPi0_%s%s",centArray[cent].Data(), runArray[cent].Data()));
            fitTCMInvYieldPi0[cent]->Write(Form("TwoComponentModelFitPi0_%s%s",centArray[cent].Data(), runArray[cent].Data()));
            // fits for eta
            fitInvYieldEta[cent]->Write(Form("TsallisFitEta_%s%s",centArray[cent].Data(), runArray[cent].Data()));
            fitTCMInvYieldEta[cent]->Write(Form("TwoComponentModelFitEta_%s%s",centArray[cent].Data(), runArray[cent].Data()));
        }
    fFitsResults.Close();

//     // **********************************************************************************************************************
//     // ************************* Saving comparison to comb for diff measurements ********************************************
//     // **********************************************************************************************************************
//
//     TString nameOutputCommonFileCompOnly    = Form("ComparisonsPaperpPb5TeVCent_%s.root", dateForOutput.Data());
//     TFile fCompResults(nameOutputCommonFileCompOnly.Data(), "RECREATE");
//
//         graphRatioPi0CombCombFitStat[cent]->Write("Pi0_RatioCombToCombFit_Stat");
//         graphRatioPi0CombCombFitSys[cent]->Write("Pi0_RatioCombToCombFit_Syst");
//         for (Int_t meth = 0; meth < 11; meth++){
//             if (graphRatioPi0IndCombFitStat[cent][meth]) graphRatioPi0IndCombFitStat[cent][meth]->Write();
//             if (graphRatioPi0IndCombFitSys[cent][meth]) graphRatioPi0IndCombFitSys[cent][meth]->Write();
//         }
//         graphRatioEtaCombCombFitStat->Write("Eta_RatioCombToCombFit_Stat");
//         graphRatioEtaCombCombFitSys->Write("Eta_RatioCombToCombFit_Syst");
//         for (Int_t meth = 0; meth < 11; meth++){
//             if (graphRatioEtaIndCombFitStat[meth]) graphRatioEtaIndCombFitStat[meth]->Write();
//             if (graphRatioEtaIndCombFitSys[meth]) graphRatioEtaIndCombFitSys[meth]->Write();
//         }
//     fCompResults.Close();
}

