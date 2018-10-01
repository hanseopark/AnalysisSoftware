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
                                            TString pathConfigsRpPbErr      = "",
                                            TString pathConfigsRCPErr       = "",
                                            TString pathConfigsRMBErr       = ""
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

    TString textALICE                           = "ALICE preliminary";
    //___________________________________ Declaration of files _____________________________________________
    TString collisionSystempPb                  = "p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
    TString collisionSystempp                   = "pp, #sqrt{#it{s}} = 5.02 TeV";
    TString outputDir                           = Form("%s/%s/CombineMesonMeasurements%spPb5TeVCent", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());
    TString outputDirFile                       = Form("%s/%s/CombineMesonMeasurements%spPb5TeVCent/Inputs", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());
    TString outputDirSupportComb                = Form("%s/%s/CombineMesonMeasurements%spPb5TeVCent/SupportingCombination", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());
    TString outputDirSupport                    = Form("%s/%s/CombineMesonMeasurements%spPb5TeVCent/Supporting", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());
    TString outputDirSupportPaper               = Form("%s/%s/CombineMesonMeasurements%spPb5TeVCent/SupportingPaper", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());

    TString nameFinalResDat                     = Form("%s/CombinedResults%s_FitResults.dat", outputDir.Data(), bWCorrection.Data());

    fstream  fileFitsOutput;
    fileFitsOutput.open(nameFinalResDat.Data(), ios::out);

    TString fileNameTheory                      = "ExternalInputpPb/Theory/TheoryCompilationPPb.root";
    TString fileNameOtherParticleInput          = "ExternalInputpPb/IdentifiedCharged/ChargedIdentifiedSpectrapPb_2018_09_20.root";
    TString fileNameOlderPPbMB                  = "ExternalInputpPb/CombNeutralMesons/CombinedResultsPaper_pPb_5023GeV_2018_05_25_Published.root";

    Int_t maxCentRun1                           = 14;
    Int_t maxCentRun2                           = 0;
    Int_t intRefMB                              = 1;
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
    TString  centArray[17]                      = { "0-100%", "0-100%", "0-20%", "20-40%", "40-60%",        "60-100%", "0-20%", "20-40%", "40-60%", "60-100%",
                                                    "0-20%", "20-40%", "40-60%", "60-100%", "0-100%",   "0-10%", "60-80%"};
    TString  centArray2[17]                     = { "0-100%", "0-100%", "0-20% V0A", "20-40% V0A", "40-60% V0A", "60-100% V0A",  "0-20% CL1", "20-40% CL1", "40-60% CL1", "60-100% CL1",
                                                    "0-20% ZNA", "20-40% ZNA", "40-60% ZNA", "60-100% ZNA", "0-100%",   "0-10%", "60-80%"};
    TString  centArray3[17]                     = { "0-100%", "0-100%", "0-20% V0A", "20-40% V0A",          "40-60% V0A", "60-100% V0A",  "0-20% CL1", "20-40% CL1", "40-60% CL1", "60-100% CL1",
                                                    "0-20% #it{N}_{coll} Pb-side,", "20-40% #it{N}_{coll} Pb-side,", "40-60% #it{N}_{coll} Pb-side,", "60-100% #it{N}_{coll} Pb-side,", "0-100%",   "0-10%", "60-80%"};
    TString  centArrayOutput[17]                = { "00100", "00100", "0020", "2040", "4060",           "60100", "0020", "2040", "4060",  "60100",
                                                    "0020", "2040", "4060", "60100", "00100",           "0010", "6080"};
    TString  centArrayCorr[17]                  = { "00100/", "00100Cent/", "0020/", "2040/", "4060/",  "60100/",  "CL1_0020/", "CL1_2040/", "CL1_4060/",  "CL1_60100/",
                                                    "ZNA_0020/", "ZNA_2040/", "ZNA_4060/", "ZNA_60100/", "", "", ""};
    TString  runArray[17]                       = { "", "Cent","", "", "",                               "", "", "", "", "",
                                                    "", "", "", "", "Run2",                             "Run2","Run2" };
    TString addCentString[17]                    = {"_V0A", "_V0A", "_V0A","_V0A","_V0A",               "_V0A", "_CL1", "_CL1", "_CL1", "_CL1",
                                                    "_ZNA", "_ZNA", "_ZNA", "_ZNA", "",                  "", ""};
//     TString  centArray[17]                      = { "0-20%", "20-40%", "40-60%", "60-100%", "0-100%",   "0-100%", "0-5%", "0-20%", "0-100%", "60-100%",
//                                                     "10-20%", "80-100%", "20-40%", "40-60%", "5-10%",   "0-10%", "60-80%"};
//     TString  centArrayOutput[17]                = { "0020", "2040", "4060", "60100", "00100",           "00100", "0005", "0020", "00100",  "60100",
//                                                     "1020", "80100", "2040", "4060", "0510",            "0010", "6080"};
//     TString  centArrayCorr[17]                  = { "0020/", "2040/", "4060/", "60100/", "00100/",       "00100Cent/", "", "", "",  "",
//                                                     "", "", "", "", "",                                 "", ""};
//     TString  runArray[17]                       = { "", "", "", "", "",                                 "Cent", "Run2", "Run2", "Run2", "Run2", "Run2",
//                                                     "Run2", "Run2", "Run2", "Run2", "Run2",             "Run2" };
//     TString addCentString[17]                    = {"_V0A", "_V0A", "_V0A","_V0A","_V0A",               "_V0A", "", "", "", "",
//                                                     "", "", "", "", "",                                 "", ""};
    Bool_t  enableCent[17]                      = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,  kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
                                                    kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE};
    Bool_t  enableCentComb[17]                  = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
                                                    kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE};
    Bool_t  enableCentRpPb[17]                  = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
                                                    kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kFALSE};
    Bool_t  enableCentRCP[17]                   = { kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kFALSE, kTRUE, kTRUE, kTRUE, kFALSE,
                                                    kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE};
    Bool_t  enableCentRMB[17]                   = { kFALSE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
                                                    kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kFALSE};
    Int_t ppRefInt[17]                          = { 0, 1, 1, 1, 1,  1, 1, 1, 1, 1,
                                                    1, 1, 1, 1, 2,  2, 2 };
    Int_t rCPRefInt[17]                         = { -1, -1, 5, 5, 5,   -1, 9, 9, 9, -1,
                                                    13, 13, 13, -1, -1,  -1, -1 };
    Int_t rMBRefInt[17]                         = { -1, -1, 1, 1, 1, 1, 1, 1, 1, 1,
                                                    1, 1, 1, 1, -1,  -1, -1 };
    Double_t rCPNColl[17]                       = { 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0,  0, 0, };
    Double_t rCPNCollErr[17]                    = { 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0,  0, 0, };
    Double_t rMBNColl[17]                       = { 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0,  0, 0, };
    Double_t rMBNCollErr[17]                    = { 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0,  0, 0, };
    TString nameCentEst[3]                      = { "V0A", "CL1", "ZNA"};
    TString nameCentEstRatios[3]                = { "V0A, ", "CL1, ", "#it{N}_{coll} Pb-side, "};


    TString  nameMeasGlobal[11]                 = { "PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal",
                                                    "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged",
                                                    "PCMOtherDataset"};
    TString  nameMeasGlobalLabel[11]            = { "PCM", "PHOS", "EMC", "PCM-PHOS", "PCM-EMC",
                                                    "PCM-Dal", "PHOS-Dal", "EMC-Dal", "EMChigh", "mEMC",
                                                    "PCMOtherDataset"};
    TString  nameConfigFileRpPbErr[17]          = { "", "", "", "", "",  "", "", "", "", "",
                                                    "", "", "", "", "",  "", ""};
    TString  nameConfigFileRCPErr[17]           = { "", "", "", "", "",  "", "", "", "", "",
                                                    "", "", "", "", "",  "", ""};
    TString  nameConfigFileRMBErr[17]           = { "", "", "", "", "",  "", "", "", "", "",
                                                    "", "", "", "", "",  "", ""};
    // **********************************************************************************************
    // ********************** Copy inputs into 1 directory for bookeeping ***************************
    // **********************************************************************************************
    gSystem->Exec("mkdir -p "+outputDirFile);
    gSystem->Exec("mkdir -p "+outputDirSupportComb);
    gSystem->Exec("mkdir -p "+outputDirSupport);
    gSystem->Exec("mkdir -p "+outputDirSupportPaper);
    if (fileNamePCM.CompareTo("")!=0 )          gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDirFile.Data()));
    if (fileNamePCMEMCAL.CompareTo("")!=0 )     gSystem->Exec(Form("cp %s %s/InputPCMEMCAL.root", fileNamePCMEMCAL.Data(), outputDirFile.Data()));
    if (fileNamePCMPHOS.CompareTo("")!=0 )      gSystem->Exec(Form("cp %s %s/InputPCMPHOS.root", fileNamePCMPHOS.Data(), outputDirFile.Data()));
    if (fileNamePHOS.CompareTo("")!=0 )         gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDirFile.Data()));
    if (fileNameEMCAL.CompareTo("")!=0 )        gSystem->Exec(Form("cp %s %s/InputEMCAL.root", fileNameEMCAL.Data(), outputDirFile.Data()));
    if (fileNameDalitz.CompareTo("")!=0 )       gSystem->Exec(Form("cp %s %s/InputDalitz.root", fileNameDalitz.Data(), outputDirFile.Data()));
    if (fileNameInterpolation.CompareTo("")!=0 )gSystem->Exec(Form("cp %s %s/InputInterpolation.root", fileNameInterpolation.Data(), outputDirFile.Data()));
    if (fileNameCorrFactors.CompareTo("")!=0 )  gSystem->Exec(Form("cp %s %s/InputCorrFactors.root", fileNameCorrFactors.Data(), outputDirFile.Data()));
    if (fileNameTheory.CompareTo("")!= 0)       gSystem->Exec(Form("cp %s %s/InputTheory.root", fileNameTheory.Data(), outputDirFile.Data()));
    if (pathConfigsRpPbErr.CompareTo("")!=0 ){
        for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
            if (!enableCentRpPb[cent]) continue;
            nameConfigFileRpPbErr[cent]               = Form("%s/configFileRpA_%s%s%s.txt", pathConfigsRpPbErr.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data());
            gSystem->Exec(Form("cp %s %s/configFileRpA_%s%s%s.txt", nameConfigFileRpPbErr[cent].Data(), outputDirFile.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data()));
            cout << "copied: " << nameConfigFileRpPbErr[cent].Data() << " if available" << endl;
        }
    }
    if (pathConfigsRCPErr.CompareTo("")!=0 ){
        for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
            if (rCPRefInt[cent] == -1)
                enableCentRCP[cent] = kFALSE;
            if (!enableCentRCP[cent]) continue;
            nameConfigFileRCPErr[cent]               = Form("%s/configFileRCP_%s%s%s.txt", pathConfigsRCPErr.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data());
            gSystem->Exec(Form("cp %s %s/configFileRCP_%s%s%s.txt", nameConfigFileRCPErr[cent].Data(), outputDirFile.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data()));
            cout << "copied: " << nameConfigFileRCPErr[cent].Data() << " if available" << endl;
        }
    }
    if (pathConfigsRMBErr.CompareTo("")!=0 ){
        for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
            if (rMBRefInt[cent] == -1)
                enableCentRMB[cent] = kFALSE;
            if (!enableCentRMB[cent]) continue;
            nameConfigFileRMBErr[cent]               = Form("%s/configFileRMB_%s%s%s.txt", pathConfigsRMBErr.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data());
            gSystem->Exec(Form("cp %s %s/configFileRMB_%s%s%s.txt", nameConfigFileRMBErr[cent].Data(), outputDirFile.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data()));
            cout << "copied: " << nameConfigFileRMBErr[cent].Data() << " if available" << endl;
        }
    }

    TString prefix2                             = "";
    if (isMC){
        prefix2                                 = "MC";
    } else {
        prefix2                                 = "Data";
    }

    TString  nameSecPi0SourceRead[4]            = {"K0S", "K0L", "Lambda", "Rest"};
    TString  nameSecPi0SourceLabel[4]           = {"K^{0}_{s}", "K^{0}_{l}", "#Lambda", "had. int."};
    Double_t maxSecCorr[4]                      = { 0.05, 0.007, 0.0003, 0.03};
    // scale factors from https://twiki.cern.ch/twiki/pub/ALICE/PAPaperCentrality/normalization.pdf
    Double_t scalingToNSD[17]                   = { 0.964, 0.964, 1, 1, 1,  0.97, 1, 1, 1, 0.97,   1, 1, 0.995, 0.9725, 0.964,     1, 1};
    Int_t modeToBeUsed[17]                      = { 21, 20, 20, 20, 20,     20, 20, 20, 20, 20,
                                                    20, 20, 20, 20, 20,     20, 20 };
    Bool_t scaleNSD[17][11]                     = {
                                                    { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 },  // MB R1
                                                    { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 },  // MB R1 cent binning
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-20
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 20-40
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 40-60
                                                    { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 },  // 60-100
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-20 CL1
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 20-40 CL1
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 40-60 CL1
                                                    { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 },  // 60-100 CL1
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-20 ZNA
                                                    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 20-40 ZNA
                                                    { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 },  // 40-60 ZNA
                                                    { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 },  // 60-100 ZNA
                                                    { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1  },  // MB R2
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
    Double_t maxPtEtaPlotting                   = 40.;
    Double_t minPtEtaToPi0Plotting              = 0.43;
    Double_t maxPtEtaToPi0Plotting              = 40.;

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

    Color_t  colorCent[17];
    Color_t  colorCentMC[17];
    Color_t  colorCentBox[17];
    Style_t  markerStyleCent[17];
    Style_t  markerStyleCentMC[17];
    Size_t   markerSizeCent[17];
    Size_t   markerSizeCentMC[17];
    Double_t nCollpPb[17];
    Double_t nCollErrpPb[17];
    Double_t tpPb[17];
    Double_t tpPbErr[17];
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        colorCent[cent]                         = GetColorDefaultColor("pPb_5.023TeV", "", centArray[cent]);
        colorCentMC[cent]                       = GetColorDefaultColor("pPb_5.023TeV", "HIJING", centArray[cent]);
        colorCentBox[cent]                      = GetColorDefaultColor("pPb_5.023TeV", "", centArray[cent],kTRUE);
        markerStyleCent[cent]                   = GetDefaultMarkerStyle("pPb_5.023TeV", "", centArray[cent]);
        markerStyleCentMC[cent]                 = GetDefaultMarkerStyle("pPb_5.023TeV", "HIJING", centArray[cent]);
        markerSizeCent[cent]                    = GetDefaultMarkerSize("pPb_5.023TeV", "", centArray[cent])*2;
        markerSizeCentMC[cent]                  = GetDefaultMarkerSize("pPb_5.023TeV", "HIJING", centArray[cent])*2;
        nCollpPb[cent]                          = GetNCollFromName(centArrayOutput[cent]+addCentString[cent], "pPb_5.023TeV");
        nCollErrpPb[cent]                       = GetNCollErrFromName(centArrayOutput[cent]+addCentString[cent], "pPb_5.023TeV");
        tpPb[cent]                              = GetTAAFromName(centArrayOutput[cent]+addCentString[cent], "pPb_5.023TeV");
        tpPbErr[cent]                           = GetTAAErrFromName(centArrayOutput[cent]+addCentString[cent], "pPb_5.023TeV");
        clog << cent << "\t" << centArrayOutput[cent]+addCentString[cent] << "\t"<< nCollpPb[cent] << "\t" << nCollErrpPb[cent] << "\t"<< tpPb[cent] << "\t"<< tpPbErr[cent] << endl;;
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
    Double_t xPtLimitsPi0[17][100];
    Int_t maxNBinsPi0Abs                            = 0;
    Int_t maxNBinsPi0[17]                           = {0};

    cout << "Setting Eta binning" << endl;
    Double_t xPtLimitsEta[17][100];
    Double_t xPtLimitsEtaToPi0[17][100];
    Int_t maxNBinsEtaAbs                            = 0;
    Int_t maxNBinsEta[17]                           = {0};
    Int_t maxNBinsEtaToPi0Abs                       = 0;
    Int_t maxNBinsEtaToPi0[17]                      = {0};
    Int_t nTotMeasPi0[17]                           = {0};
    Int_t nTotMeasEta[17]                           = {0};
    Int_t nTotMeasEtaToPi0[17]                      = {0};
    TString fileNamesMethod[11]                     = {"", "", "", "", "", "", "", "", "", "", ""};
    TString fileNamespPbPi0DetailedSys[17][11];
    TString fileNamesppPi0DetailedSys[17][11];
    TString fileNamesperiPi0DetailedSys[17][11];
    TString fileNamesminPi0DetailedSys[17][11];
    TString fileNamesRpPbPi0DetailedSys[17][11];
    TString fileNamesRCPPi0DetailedSys[17][11];
    TString fileNamesRMBPi0DetailedSys[17][11];
    TString fileNamespPbEtaDetailedSys[17][11];
    TString fileNamesppEtaDetailedSys[17][11];
    TString fileNamesperiEtaDetailedSys[17][11];
    TString fileNamesminEtaDetailedSys[17][11];
    TString fileNamesRpPbEtaDetailedSys[17][11];
    TString fileNamesRCPEtaDetailedSys[17][11];
    TString fileNamesRMBEtaDetailedSys[17][11];
    Bool_t havePi0SysDetailedpPb[17][11];
    Bool_t havePi0SysDetailedpp[17][11];
    Bool_t havePi0SysDetailedperi[17][11];
    Bool_t havePi0SysDetailedmin[17][11];
    Bool_t haveEtaSysDetailedpPb[17][11];
    Bool_t haveEtaSysDetailedpp[17][11];
    Bool_t haveEtaSysDetailedperi[17][11];
    Bool_t haveEtaSysDetailedmin[17][11];
    vector<TString> ptSysRemNames[17][11];
    vector<TString> ptSysRemRCPNames[17][11];
    vector<TString> ptSysRemRMBNames[17][11];
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        for (Int_t meth = 0; meth < 11; meth++){
            fileNamespPbPi0DetailedSys[cent][meth]      = "";
            fileNamesppPi0DetailedSys[cent][meth]       = "";
            fileNamesperiPi0DetailedSys[cent][meth]     = "";
            fileNamesminPi0DetailedSys[cent][meth]      = "";
            fileNamesRpPbPi0DetailedSys[cent][meth]     = "";
            fileNamesRCPPi0DetailedSys[cent][meth]      = "";
            fileNamesRMBPi0DetailedSys[cent][meth]      = "";
            fileNamespPbEtaDetailedSys[cent][meth]      = "";
            fileNamesppEtaDetailedSys[cent][meth]       = "";
            fileNamesperiEtaDetailedSys[cent][meth]       = "";
            fileNamesminEtaDetailedSys[cent][meth]       = "";
            fileNamesRpPbEtaDetailedSys[cent][meth]     = "";
            fileNamesRCPEtaDetailedSys[cent][meth]      = "";
            fileNamesRMBEtaDetailedSys[cent][meth]      = "";
            havePi0SysDetailedpPb[cent][meth]           = kFALSE;
            havePi0SysDetailedpp[cent][meth]            = kFALSE;
            havePi0SysDetailedperi[cent][meth]          = kFALSE;
            havePi0SysDetailedmin[cent][meth]           = kFALSE;
            haveEtaSysDetailedpPb[cent][meth]           = kFALSE;
            haveEtaSysDetailedpp[cent][meth]            = kFALSE;
            haveEtaSysDetailedperi[cent][meth]          = kFALSE;
            haveEtaSysDetailedmin[cent][meth]           = kFALSE;
        }
    }

    Bool_t haveEffSecCorr[17][4][11];
    TFile* fileMethod[11];
    TDirectory* directoryPi0[17][11];
    TGraphAsymmErrors* graphPi0EffSecCorrFromX[17][4][11];
    TGraphAsymmErrors* graphPi0Eff[17][11];
    TGraphAsymmErrors* graphPi0Acc[17][11];
    TGraphAsymmErrors* graphPi0EffTimesAcc[17][11];
    TGraphAsymmErrors* graphPi0Mass[17][11];
    TGraphAsymmErrors* graphPi0MassMC[17][11];
    TGraphAsymmErrors* graphPi0Width[17][11];
    TGraphAsymmErrors* graphPi0WidthMC[17][11];
    TH1D* histoPi0InvYieldStat[17][11];
    TGraphAsymmErrors* graphPi0InvYieldStat[17][11];
    TGraphAsymmErrors* graphPi0InvYieldSys[17][11];
    TDirectory* directoryEta[17][11];
    TGraphAsymmErrors* graphEtaEff[17][11];
    TGraphAsymmErrors* graphEtaAcc[17][11];
    TGraphAsymmErrors* graphEtaEffTimesAcc[17][11];
    TGraphAsymmErrors* graphEtaMass[17][11];
    TGraphAsymmErrors* graphEtaMassMC[17][11];
    TGraphAsymmErrors* graphEtaWidth[17][11];
    TGraphAsymmErrors* graphEtaWidthMC[17][11];
    TH1D* histoEtaInvYieldStat[17][11];
    TGraphAsymmErrors* graphEtaInvYieldStat[17][11];
    TGraphAsymmErrors* graphEtaInvYieldSys[17][11];
    TH1D* histoEtaToPi0Stat[17][11];
    TGraphAsymmErrors* graphEtaToPi0Stat[17][11];
    TGraphAsymmErrors* graphEtaToPi0Sys[17][11];

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        maxNBinsPi0[cent]                           = GetBinning( xPtLimitsPi0[cent], maxNBinsPi0Abs, "Pi0", Form("pPb_5.023TeV%s",runArray[cent].Data()), modeToBeUsed[cent], -1, kFALSE, centArray[cent]);
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
        maxNBinsEta[cent]                           = GetBinning( xPtLimitsEta[cent], maxNBinsEtaAbs, "Eta", Form("pPb_5.023TeV%s",runArray[cent].Data()), modeToBeUsed[cent], -1, kFALSE, centArray[cent]);
        maxNBinsEtaToPi0[cent]                      = GetBinning( xPtLimitsEtaToPi0[cent], maxNBinsEtaToPi0Abs, "Eta", Form("pPb_5.023TeV%s",runArray[cent].Data()), modeToBeUsed[cent], -1, kFALSE, centArray[cent]);
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

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){

        if (nameConfigFileRpPbErr[cent].CompareTo("") != 0){
            cout << "===================================================================================================================" << endl;
            cout << "RpPb details for " << centArray[cent].Data() << "\t" << addCentString[cent].Data() << "\t" << runArray[cent].Data() << endl;
            ifstream inPbConfig(nameConfigFileRpPbErr[cent].Data());
            Int_t nReadSys  = 0;
            while(!inPbConfig.eof() && nReadSys<11 ){
                cout << "read line:" <<  nReadSys << "\t"<< nameMeasGlobalLabel[nReadSys].Data() << endl;
                TString nameCurrentSys  = "";
                inPbConfig >> nameCurrentSys ;
                if (nameCurrentSys.CompareTo(nameMeasGlobalLabel[nReadSys].Data()) != 0){
                    cout << "wrong order in configuration file, was expecting " <<  nameMeasGlobalLabel[nReadSys].Data() << endl;
                    return;
                } else {
                    inPbConfig >> fileNamespPbPi0DetailedSys[cent][nReadSys] >> fileNamesppPi0DetailedSys[cent][nReadSys] >> fileNamespPbEtaDetailedSys[cent][nReadSys] >> fileNamesppEtaDetailedSys[cent][nReadSys];
                    if (fileNamespPbPi0DetailedSys[cent][nReadSys].CompareTo("bla") != 0){
                        cout << "pi0 pPb sys detailed: " << fileNamespPbPi0DetailedSys[cent][nReadSys] << endl;
                        havePi0SysDetailedpPb[cent][nReadSys]         = kTRUE;
                    } else {
                        fileNamespPbPi0DetailedSys[cent][nReadSys]    = "";
                    }
                    if (fileNamesppPi0DetailedSys[cent][nReadSys].CompareTo("bla") != 0){
                        cout << "pi0 pp sys detailed: " << fileNamesppPi0DetailedSys[cent][nReadSys] << endl;
                        havePi0SysDetailedpp[cent][nReadSys]          = kTRUE;
                    } else {
                        fileNamesppPi0DetailedSys[cent][nReadSys]     = "";
                    }
                    if (fileNamespPbEtaDetailedSys[cent][nReadSys].CompareTo("bla") != 0){
                        cout << "eta pPb sys detailed: " << fileNamespPbEtaDetailedSys[cent][nReadSys] << endl;
                        haveEtaSysDetailedpPb[cent][nReadSys]         = kTRUE;
                    } else {
                        fileNamespPbEtaDetailedSys[cent][nReadSys]    = "";
                    }
                    if (fileNamesppEtaDetailedSys[cent][nReadSys].CompareTo("bla") != 0){
                        cout << "eta pp sys detailed: " << fileNamesppEtaDetailedSys[cent][nReadSys] << endl;
                        haveEtaSysDetailedpp[cent][nReadSys]          = kTRUE;
                    } else {
                        fileNamesppEtaDetailedSys[cent][nReadSys]     = "";
                    }
                    if (!(havePi0SysDetailedpPb[cent][nReadSys] && havePi0SysDetailedpp[cent][nReadSys])){
                        havePi0SysDetailedpPb[cent][nReadSys]         = kTRUE;
                        havePi0SysDetailedpp[cent][nReadSys]          = kTRUE;
                        fileNamespPbPi0DetailedSys[cent][nReadSys]    = "";
                        fileNamesppPi0DetailedSys[cent][nReadSys]     = "";
                    } else {
                        fileNamesRpPbPi0DetailedSys[cent][nReadSys]  = Form("%s/Pi0RpPb_%s_detailedSys%s%s%s.dat", outputDirSupportComb.Data(),nameMeasGlobalLabel[nReadSys].Data(),
                                                                            centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data());
                    }
                    if (!(haveEtaSysDetailedpPb[cent][nReadSys] && haveEtaSysDetailedpp[cent][nReadSys])){
                        haveEtaSysDetailedpPb[cent][nReadSys]         = kTRUE;
                        haveEtaSysDetailedpp[cent][nReadSys]          = kTRUE;
                        fileNamespPbEtaDetailedSys[cent][nReadSys]    = "";
                        fileNamesppEtaDetailedSys[cent][nReadSys]     = "";
                    } else {
                        fileNamesRpPbEtaDetailedSys[cent][nReadSys]  = Form("%s/EtaRpPb_%s_detailedSys%s%s%s.dat", outputDirSupportComb.Data(),nameMeasGlobalLabel[nReadSys].Data(),
                                                                            centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data());
                    }
                    TString currentString = "";
                    inPbConfig >> currentString;
                    Int_t counter = 0;
                    while (currentString.CompareTo("STOP") != 0 && counter < 20 ){
                        ptSysRemNames[cent][nReadSys].push_back(currentString);
                        counter++;
                        inPbConfig >> currentString;
                    }
                    cout << "need to take out "<< ptSysRemNames[cent][nReadSys].size() << " systematic errors" << endl;
                }
                nReadSys++;
            }
            cout << "===================================================================================================================" << endl;
        }
        if (nameConfigFileRCPErr[cent].CompareTo("") != 0){
            cout << "===================================================================================================================" << endl;
            cout << "RCP details for " << centArray[cent].Data() << "\t" << addCentString[cent].Data() << "\t" << runArray[cent].Data() << endl;
            ifstream inPbConfig(nameConfigFileRCPErr[cent].Data());
            Int_t nReadSys  = 0;
            while(!inPbConfig.eof() && nReadSys<11 ){
                cout << "read line:" <<  nReadSys << "\t"<< nameMeasGlobalLabel[nReadSys].Data() << endl;
                TString nameCurrentSys  = "";
                inPbConfig >> nameCurrentSys ;
                if (nameCurrentSys.CompareTo(nameMeasGlobalLabel[nReadSys].Data()) != 0){
                    cout << "wrong order in configuration file, was expecting " <<  nameMeasGlobalLabel[nReadSys].Data() << endl;
                    return;
                } else {
                    inPbConfig >> fileNamesperiPi0DetailedSys[cent][nReadSys] >>  fileNamesperiEtaDetailedSys[cent][nReadSys];
                    if (fileNamesperiPi0DetailedSys[cent][nReadSys].CompareTo("bla") != 0){
                        cout << "pi0 peri sys detailed: " << fileNamesperiPi0DetailedSys[cent][nReadSys] << endl;
                        havePi0SysDetailedperi[cent][nReadSys]          = kTRUE;
                    } else {
                        fileNamesperiPi0DetailedSys[cent][nReadSys]     = "";
                    }
                    if (fileNamesperiEtaDetailedSys[cent][nReadSys].CompareTo("bla") != 0){
                        cout << "eta peri sys detailed: " << fileNamesperiEtaDetailedSys[cent][nReadSys] << endl;
                        haveEtaSysDetailedperi[cent][nReadSys]          = kTRUE;
                    } else {
                        fileNamesperiEtaDetailedSys[cent][nReadSys]     = "";
                    }
                    if (!(havePi0SysDetailedpPb[cent][nReadSys] && havePi0SysDetailedperi[cent][nReadSys])){
                        havePi0SysDetailedperi[cent][nReadSys]          = kTRUE;
                        fileNamesperiPi0DetailedSys[cent][nReadSys]     = "";
                    } else {
                        fileNamesRCPPi0DetailedSys[cent][nReadSys]  = Form("%s/Pi0RCP_%s_detailedSys%s%s%s.dat", outputDirSupportComb.Data(),nameMeasGlobalLabel[nReadSys].Data(),
                                                                            centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data());
                    }
                    if (!(haveEtaSysDetailedpPb[cent][nReadSys] && haveEtaSysDetailedperi[cent][nReadSys])){
                        haveEtaSysDetailedperi[cent][nReadSys]          = kTRUE;
                        fileNamesperiEtaDetailedSys[cent][nReadSys]     = "";
                    } else {
                        fileNamesRCPEtaDetailedSys[cent][nReadSys]  = Form("%s/EtaRCP_%s_detailedSys%s%s%s.dat", outputDirSupportComb.Data(),nameMeasGlobalLabel[nReadSys].Data(),
                                                                            centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data());
                    }
                    TString currentString = "";
                    inPbConfig >> currentString;
                    Int_t counter = 0;
                    while (currentString.CompareTo("STOP") != 0 && counter < 30 ){
                        ptSysRemRCPNames[cent][nReadSys].push_back(currentString);
                        counter++;
                        inPbConfig >> currentString;
                    }
                    cout << "need to take out "<< ptSysRemRCPNames[cent][nReadSys].size() << " systematic errors" << endl;
                }
                nReadSys++;
            }
            cout << "===================================================================================================================" << endl;
        }
        if (nameConfigFileRMBErr[cent].CompareTo("") != 0){
            cout << "===================================================================================================================" << endl;
            cout << "RMB details for " << centArray[cent].Data() << "\t" << addCentString[cent].Data() << "\t" << runArray[cent].Data() << endl;
            ifstream inPbConfig(nameConfigFileRMBErr[cent].Data());
            Int_t nReadSys  = 0;
            while(!inPbConfig.eof() && nReadSys<11 ){
                cout << "read line:" <<  nReadSys << "\t"<< nameMeasGlobalLabel[nReadSys].Data() << endl;
                TString nameCurrentSys  = "";
                inPbConfig >> nameCurrentSys ;
                if (nameCurrentSys.CompareTo(nameMeasGlobalLabel[nReadSys].Data()) != 0){
                    cout << "wrong order in configuration file, was expecting " <<  nameMeasGlobalLabel[nReadSys].Data() << endl;
                    return;
                } else {
                    inPbConfig >> fileNamesminPi0DetailedSys[cent][nReadSys] >>  fileNamesminEtaDetailedSys[cent][nReadSys];
                    if (fileNamesminPi0DetailedSys[cent][nReadSys].CompareTo("bla") != 0){
                        cout << "pi0 min sys detailed: " << fileNamesminPi0DetailedSys[cent][nReadSys] << endl;
                        havePi0SysDetailedmin[cent][nReadSys]          = kTRUE;
                    } else {
                        fileNamesminPi0DetailedSys[cent][nReadSys]     = "";
                    }
                    if (fileNamesminEtaDetailedSys[cent][nReadSys].CompareTo("bla") != 0){
                        cout << "eta min sys detailed: " << fileNamesminEtaDetailedSys[cent][nReadSys] << endl;
                        haveEtaSysDetailedmin[cent][nReadSys]          = kTRUE;
                    } else {
                        fileNamesminEtaDetailedSys[cent][nReadSys]     = "";
                    }
                    if (!(havePi0SysDetailedpPb[cent][nReadSys] && havePi0SysDetailedmin[cent][nReadSys])){
                        havePi0SysDetailedmin[cent][nReadSys]          = kTRUE;
                        fileNamesminPi0DetailedSys[cent][nReadSys]     = "";
                    } else {
                        fileNamesRMBPi0DetailedSys[cent][nReadSys]  = Form("%s/Pi0RMB_%s_detailedSys%s%s%s.dat", outputDirSupportComb.Data(),nameMeasGlobalLabel[nReadSys].Data(),
                                                                           centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data());
                    }
                    if (!(haveEtaSysDetailedpPb[cent][nReadSys] && haveEtaSysDetailedmin[cent][nReadSys])){
                        haveEtaSysDetailedmin[cent][nReadSys]          = kTRUE;
                        fileNamesminEtaDetailedSys[cent][nReadSys]     = "";
                    } else {
                        fileNamesRMBEtaDetailedSys[cent][nReadSys]  = Form("%s/EtaRMB_%s_detailedSys%s%s%s.dat", outputDirSupportComb.Data(),nameMeasGlobalLabel[nReadSys].Data(),
                                                                           centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data());
                    }
                    TString currentString = "";
                    inPbConfig >> currentString;
                    Int_t counter = 0;
                    while (currentString.CompareTo("STOP") != 0 && counter < 30 ){
                        ptSysRemRMBNames[cent][nReadSys].push_back(currentString);
                        counter++;
                        inPbConfig >> currentString;
                    }
                    cout << "need to take out "<< ptSysRemRMBNames[cent][nReadSys].size() << " systematic errors" << endl;
                }
                nReadSys++;
            }
            cout << "===================================================================================================================" << endl;
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
                directoryPi0[cent][meth]                                = (TDirectory*)fileMethod[meth]->Get(Form("Pi0%spPb_5.023TeV%s%s",centArray[cent].Data(), runArray[cent].Data(),
                                                                                                                  addCentString[cent].Data()));
                if (runArray[cent].CompareTo("Cent")==0 && meth == 1)
                    directoryPi0[cent][meth]                                = (TDirectory*)fileMethod[meth]->Get(Form("Pi0%sCentDeppPb_5.023TeV%s",centArray[cent].Data(),
                                                                                                                      addCentString[cent].Data()));

                if (!directoryPi0[cent][meth]) {
                    cout << "File doesn't contain directory " << Form("Pi0%spPb_5.023TeV%s%s",centArray[cent].Data(), runArray[cent].Data(), addCentString[cent].Data()) << "! Skipping " << nameMeasGlobalLabel[meth].Data()  << endl;
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

                    if (graphPi0Width[cent][meth]->GetY()[0] == 0 && graphPi0Width[cent][meth]->GetN()>0){
                        graphPi0Width[cent][meth]->RemovePoint(0);
                        graphPi0Mass[cent][meth]->RemovePoint(0);
                        graphPi0MassMC[cent][meth]->RemovePoint(0);
                        graphPi0WidthMC[cent][meth]->RemovePoint(0);
                        graphPi0EffTimesAcc[cent][meth]->RemovePoint(0);
                    }
                    while (graphPi0Width[cent][meth]->GetY()[graphPi0Width[cent][meth]->GetN()-1] == 0 && graphPi0Width[cent][meth]->GetN()>0){
                        graphPi0Width[cent][meth]->RemovePoint(graphPi0Width[cent][meth]->GetN()-1);
                        graphPi0Mass[cent][meth]->RemovePoint(graphPi0Mass[cent][meth]->GetN()-1);
                        graphPi0MassMC[cent][meth]->RemovePoint(graphPi0MassMC[cent][meth]->GetN()-1);
                        graphPi0WidthMC[cent][meth]->RemovePoint(graphPi0WidthMC[cent][meth]->GetN()-1);
                        graphPi0EffTimesAcc[cent][meth]->RemovePoint(graphPi0EffTimesAcc[cent][meth]->GetN()-1);
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
                        if (histoPi0InvYieldStat[cent][meth]) histoPi0InvYieldStat[cent][meth]->Scale(scalingToNSD[cent]);
                        if (graphPi0InvYieldStat[cent][meth]) graphPi0InvYieldStat[cent][meth]      = ScaleGraph(graphPi0InvYieldStat[cent][meth],scalingToNSD[cent]);
                        if (graphPi0InvYieldSys[cent][meth]) graphPi0InvYieldSys[cent][meth]        = ScaleGraph(graphPi0InvYieldSys[cent][meth],scalingToNSD[cent]);
                    }
                }
                directoryEta[cent][meth]                             = (TDirectory*)fileMethod[meth]->Get(Form("Eta%spPb_5.023TeV%s%s",centArray[cent].Data(), runArray[cent].Data(),
                                                                                                               addCentString[cent].Data()));
                if (runArray[cent].CompareTo("Cent") == 0 && meth == 1)
                    directoryEta[cent][meth]                         = (TDirectory*)fileMethod[meth]->Get(Form("Eta%sCentDeppPb_5.023TeV%s",centArray[cent].Data(),
                                                                                                                      addCentString[cent].Data()));
                if (!directoryEta[cent][meth]) {
                    cout << "File doesn't contain directory " << Form("Eta%spPb_5.023TeV%s%s",centArray[cent].Data(), runArray[cent].Data(), addCentString[cent].Data()) << "! Skipping " << nameMeasGlobalLabel[meth].Data()  << endl;
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
                    // reading eta yields
                    graphEtaInvYieldStat[cent][meth]                     = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("graphCorrectedYieldEta");
                    histoEtaInvYieldStat[cent][meth]                     = (TH1D*)directoryEta[cent][meth]->Get("CorrectedYieldEta");
                    cout << "Eta stat error" << endl;
                    graphEtaInvYieldStat[cent][meth]->Print();
                    graphEtaInvYieldSys[cent][meth]                      = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("EtaSystError");
                    if (scaleNSD[cent][meth]){
                        if (histoEtaInvYieldStat[cent][meth]) histoEtaInvYieldStat[cent][meth]->Scale(scalingToNSD[cent]);
                        if (graphEtaInvYieldStat[cent][meth]) graphEtaInvYieldStat[cent][meth]                 = ScaleGraph(graphEtaInvYieldStat[cent][meth],scalingToNSD[cent]);
                        if (graphEtaInvYieldSys[cent][meth]) graphEtaInvYieldSys[cent][meth]                  = ScaleGraph(graphEtaInvYieldSys[cent][meth],scalingToNSD[cent]);
                    }

                    // remove 0 entries from effi, acc, width, mass
                    while (graphEtaWidth[cent][meth]->GetY()[0] == 0 && graphEtaWidth[cent][meth]->GetN()>0){
                        graphEtaWidth[cent][meth]->RemovePoint(0);
                        graphEtaMass[cent][meth]->RemovePoint(0);
                        graphEtaMassMC[cent][meth]->RemovePoint(0);
                        graphEtaWidthMC[cent][meth]->RemovePoint(0);
                        graphEtaEffTimesAcc[cent][meth]->RemovePoint(0);
                    }
                    while (graphEtaWidth[cent][meth]->GetY()[graphEtaWidth[cent][meth]->GetN()-1] == 0 && graphEtaWidth[cent][meth]->GetN()>0){
                        graphEtaWidth[cent][meth]->RemovePoint(graphEtaWidth[cent][meth]->GetN()-1);
                        graphEtaMass[cent][meth]->RemovePoint(graphEtaMass[cent][meth]->GetN()-1);
                        graphEtaMassMC[cent][meth]->RemovePoint(graphEtaMassMC[cent][meth]->GetN()-1);
                        graphEtaWidthMC[cent][meth]->RemovePoint(graphEtaWidthMC[cent][meth]->GetN()-1);
                        graphEtaEffTimesAcc[cent][meth]->RemovePoint(graphEtaEffTimesAcc[cent][meth]->GetN()-1);
                    }

                    // loading eta/pi0 hists/graphs
                    cout << "Eta sys error" << endl;
                    if (graphEtaInvYieldSys[cent][meth]) graphEtaInvYieldSys[cent][meth]->Print();
                    histoEtaToPi0Stat[cent][meth]                        = (TH1D*)directoryEta[cent][meth]->Get("EtaToPi0StatError");
                    if (!histoEtaToPi0Stat[cent][meth])
                        histoEtaToPi0Stat[cent][meth]                    = (TH1D*)directoryEta[cent][meth]->Get("EtaToPi0YShiftedStatError");
                    graphEtaToPi0Stat[cent][meth]                        = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("graphEtaToPi0StatError");
                    if (!graphEtaToPi0Stat[cent][meth])
                        graphEtaToPi0Stat[cent][meth]                    = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("graphEtaToPi0YShiftedStatError");
                    graphEtaToPi0Sys[cent][meth]                         = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("EtaToPi0SystError");
                    if (!graphEtaToPi0Sys[cent][meth])
                        graphEtaToPi0Sys[cent][meth]                     = (TGraphAsymmErrors*)directoryEta[cent][meth]->Get("EtaToPi0YShiftedSystError");
                    if (graphEtaToPi0Stat[cent][meth])
                        nTotMeasEtaToPi0[cent]++;
                    // manually set bins below 4 GeV to 0 for PHOS for Run1
                    if (meth == 1 && cent < 5){
                        Int_t ptBin = 1;
                        while (histoEtaInvYieldStat[cent][meth]->GetBinCenter(ptBin) < 4){
                            histoEtaInvYieldStat[cent][meth]->SetBinContent(ptBin, 0);
                            ptBin++;
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
    // ************************** Loading theory spectra *****************************************************
    // *******************************************************************************************************
    TFile* fileTheory                                   = new TFile(fileNameTheory.Data());
    TDirectory* directoryTheory[17]                     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1F* histoDPMJetPi0[17]                            = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1F* histoDPMJetEta[17]                            = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1F* histoDPMJetEtaToPi0[17]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1F* histoDPMJetPi0ToPiCh[17]                      = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1F* histoDPMJetEtaToKCh[17]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1F* histoEPOSLHCPi0[17]                           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1F* histoEPOSLHCEta[17]                           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1F* histoEPOSLHCEtaToPi0[17]                      = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1F* histoEPOSLHCPi0ToPiCh[17]                     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1F* histoEPOSLHCEtaToKCh[17]                      = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1F* histoHIJINGPi0[17]                            = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1F* histoHIJINGEta[17]                            = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1F* histoHIJINGEtaToPi0[17]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1F* histoHIJINGPi0ToPiCh[17]                      = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1F* histoHIJINGEtaToKCh[17]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent ++){
        if (!enableCent[cent] || cent > maxCentRun1) continue;
        cout << "**********************************************************************" << endl;
        cout << "Reading in cent " << centArray[cent].Data() << endl;
        cout << "**********************************************************************" << endl;
        TString currentCent     =  centArray[cent].Data();
        TString addTheo         = "";
        if (centArray[cent].CompareTo("0-100%") == 0) currentCent = "";
        directoryTheory[cent]                               = (TDirectory*)fileTheory->Get(Form("%spPb_5.023TeV",currentCent.Data()));
        if (directoryTheory[cent]){
            if (centArray[cent].CompareTo("0-100%") == 0){
                histoDPMJetPi0[cent]                                = (TH1F*) directoryTheory[cent]->Get("histoPi0SpecDPMJet_MCGen");
                histoDPMJetEta[cent]                                = (TH1F*) directoryTheory[cent]->Get("histoEtaSpecDPMJet_MCGen");
                histoDPMJetEtaToPi0[cent]                           = (TH1F*) directoryTheory[cent]->Get("histoEtaToPi0DPMJet_MCGen");
                histoDPMJetPi0ToPiCh[cent]                          = (TH1F*) directoryTheory[cent]->Get("histoPi0ToPiChDPMJet_MCGen");
                histoEPOSLHCPi0[cent]                               = (TH1F*) directoryTheory[cent]->Get("histoPi0SpecEPOSLHC_MCGen");
                histoEPOSLHCEta[cent]                               = (TH1F*) directoryTheory[cent]->Get("histoEtaSpecEPOSLHC_MCGen");
                histoEPOSLHCEtaToPi0[cent]                          = (TH1F*) directoryTheory[cent]->Get("histoEtaToPi0EPOSLHC_MCGen");
                histoEPOSLHCPi0ToPiCh[cent]                         = (TH1F*) directoryTheory[cent]->Get("histoPi0ToPiChEPOSLHC_MCGen");
                histoHIJINGPi0[cent]                                = (TH1F*) directoryTheory[cent]->Get("histoPi0SpecHIJING_MCGen");
                histoHIJINGEta[cent]                                = (TH1F*) directoryTheory[cent]->Get("histoEtaSpecHIJING_MCGen");
                histoHIJINGEtaToPi0[cent]                           = (TH1F*) directoryTheory[cent]->Get("histoEtaToPi0HIJING_MCGen");
                histoHIJINGPi0ToPiCh[cent]                          = (TH1F*) directoryTheory[cent]->Get("histoPi0ToPiChHIJING_MCGen");
                histoHIJINGEtaToKCh[cent]                           = (TH1F*) directoryTheory[cent]->Get("histoEtaToKChHIJING");
            } else {
                histoDPMJetPi0[cent]                                = (TH1F*) directoryTheory[cent]->Get("histoPi0SpecDPMJet_Reb");
                histoDPMJetEta[cent]                                = (TH1F*) directoryTheory[cent]->Get("histoEtaSpecDPMJet_Reb");
                histoDPMJetEtaToPi0[cent]                           = (TH1F*) directoryTheory[cent]->Get("histoEtaToPi0DPMJet");
                histoDPMJetPi0ToPiCh[cent]                          = (TH1F*) directoryTheory[cent]->Get("histoPi0ToPiChDPMJet");
                histoEPOSLHCPi0[cent]                               = (TH1F*) directoryTheory[cent]->Get("histoPi0SpecEPOSLHC_Reb");
                histoEPOSLHCEta[cent]                               = (TH1F*) directoryTheory[cent]->Get("histoEtaSpecEPOSLHC_Reb");
                histoEPOSLHCEtaToPi0[cent]                          = (TH1F*) directoryTheory[cent]->Get("histoEtaToPi0EPOSLHC");
                histoEPOSLHCPi0ToPiCh[cent]                         = (TH1F*) directoryTheory[cent]->Get("histoPi0ToPiChEPOSLHC");
            }
            histoDPMJetEtaToKCh[cent]                           = (TH1F*) directoryTheory[cent]->Get("histoEtaToKChDPMJet");
            histoEPOSLHCEtaToKCh[cent]                          = (TH1F*) directoryTheory[cent]->Get("histoEtaToKChEPOSLHC");
        }
    }
    // *******************************************************************************************************
    // ************************** Loading interpolated spectra ***********************************************
    // *******************************************************************************************************
    TString addNamePP[2]                                    = {"", "Cent"};
    TGraphAsymmErrors* statErrorCollectionPi0PP[2][11];
    TGraphAsymmErrors* systErrorCollectionPi0PP[2][11];
    TGraphAsymmErrors* systErrorUnCorrCollectionPi0PP[2][11];
    TGraphAsymmErrors* systErrorInterCollectionPi0PP[2][11];
    TGraphAsymmErrors* systErrorUnCorrInterCollectionPi0PP[2][11];
    TGraphAsymmErrors* statErrorCollectionInvYieldPi0PP[2][11];
    TGraphAsymmErrors* systErrorCollectionInvYieldPi0PP[2][11];
    TGraphAsymmErrors* systErrorUnCorrCollectionInvYieldPi0PP[2][11];
    TGraphAsymmErrors* systErrorInterCollectionInvYieldPi0PP[2][11];
    TGraphAsymmErrors* systErrorUnCorrInterCollectionInvYieldPi0PP[2][11];

    TGraphAsymmErrors* statErrorCollectionEtaPP[2][11];
    TGraphAsymmErrors* systErrorCollectionEtaPP[2][11];
    TGraphAsymmErrors* systErrorUnCorrCollectionEtaPP[2][11];
    TGraphAsymmErrors* systErrorInterCollectionEtaPP[2][11];
    TGraphAsymmErrors* systErrorUnCorrInterCollectionEtaPP[2][11];
    TGraphAsymmErrors* statErrorCollectionInvYieldEtaPP[2][11];
    TGraphAsymmErrors* systErrorCollectionInvYieldEtaPP[2][11];
    TGraphAsymmErrors* systErrorUnCorrCollectionInvYieldEtaPP[2][11];
    TGraphAsymmErrors* systErrorInterCollectionInvYieldEtaPP[2][11];
    TGraphAsymmErrors* systErrorUnCorrInterCollectionInvYieldEtaPP[2][11];
    TGraphAsymmErrors* graphRatioPi0IndCombFitStatPP[2][11];
    TGraphAsymmErrors* graphRatioPi0IndCombFitSysPP[2][11];
    TGraphAsymmErrors* graphRatioEtaIndCombFitStatPP[2][11];
    TGraphAsymmErrors* graphRatioEtaIndCombFitSysPP[2][11];
    TGraphAsymmErrors* graphPPCombPi0Stat[2]                = {NULL, NULL};
    TGraphAsymmErrors* graphPPCombPi0UncorrSys[2]           = {NULL, NULL};
    TGraphAsymmErrors* graphPPCombPi0FullSys[2]             = {NULL, NULL};
    TGraphAsymmErrors* graphPPCombPi0InterSys[2]            = {NULL, NULL};
    TGraphAsymmErrors* graphPPInvYieldCombPi0Stat[2]        = {NULL, NULL};
    TGraphAsymmErrors* graphPPInvYieldCombPi0StatW0XErr[2]  = {NULL, NULL};
    TGraphAsymmErrors* graphPPInvYieldCombPi0Sys[2]         = {NULL, NULL};
    TGraphAsymmErrors* graphPPCombEtaStat[2]                = {NULL, NULL};
    TGraphAsymmErrors* graphPPCombEtaUncorrSys[2]           = {NULL, NULL};
    TGraphAsymmErrors* graphPPCombEtaFullSys[2]             = {NULL, NULL};
    TGraphAsymmErrors* graphPPCombEtaInterSys[2]            = {NULL, NULL};
    TGraphAsymmErrors* graphPPInvYieldCombEtaStat[2]        = {NULL, NULL};
    TGraphAsymmErrors* graphPPInvYieldCombEtaStatW0XErr[2]  = {NULL, NULL};
    TGraphAsymmErrors* graphPPInvYieldCombEtaSys[2]         = {NULL, NULL};
    TGraphAsymmErrors* graphRatioPi0CombToFitStatPP[2]      = {NULL, NULL};
    TGraphAsymmErrors* graphRatioPi0CombToFitSysPP[2]       = {NULL, NULL};
    TGraphAsymmErrors* graphRatioEtaCombToFitStatPP[2]      = {NULL, NULL};
    TGraphAsymmErrors* graphRatioEtaCombToFitSysPP[2]       = {NULL, NULL};
    TF1* fitPPTCMInvYieldPi0[2]                             = {NULL, NULL};
    TF1* fitTCMInvXSetionPi0PP[2]                           = {NULL, NULL};
    TF1* fitTsallisInvXSectionPi0PP[2]                      = {NULL, NULL};
    TF1* fitPPTCMInvYieldEta[2]                             = {NULL, NULL};
    TF1* fitTCMInvXSetionEtaPP[2]                           = {NULL, NULL};
    TF1* fitTsallisInvXSectionEtaPP[2]                      = {NULL, NULL};


    Bool_t haveRefPPPi0[2][11]                              = { { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE },
                                                                { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE } };
    Bool_t haveRefPPEta[2][11]                              = { { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE },
                                                                { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE } };
    for (Int_t ppEst = 0; ppEst < 2; ppEst++){
        for (Int_t meth = 0; meth< 11; meth++){
            statErrorCollectionPi0PP[ppEst][meth]                   = NULL;
            systErrorCollectionPi0PP[ppEst][meth]                   = NULL;
            systErrorUnCorrCollectionPi0PP[ppEst][meth]             = NULL;
            systErrorInterCollectionPi0PP[ppEst][meth]              = NULL;
            systErrorUnCorrInterCollectionPi0PP[ppEst][meth]        = NULL;
            statErrorCollectionInvYieldPi0PP[ppEst][meth]           = NULL;
            systErrorCollectionInvYieldPi0PP[ppEst][meth]           = NULL;
            systErrorUnCorrCollectionInvYieldPi0PP[ppEst][meth]     = NULL;
            systErrorInterCollectionInvYieldPi0PP[ppEst][meth]      = NULL;
            systErrorUnCorrInterCollectionInvYieldPi0PP[ppEst][meth]= NULL;
            statErrorCollectionEtaPP[ppEst][meth]                   = NULL;
            systErrorCollectionEtaPP[ppEst][meth]                   = NULL;
            systErrorUnCorrCollectionEtaPP[ppEst][meth]             = NULL;
            systErrorInterCollectionEtaPP[ppEst][meth]              = NULL;
            systErrorUnCorrInterCollectionEtaPP[ppEst][meth]        = NULL;
            statErrorCollectionInvYieldEtaPP[ppEst][meth]           = NULL;
            systErrorCollectionInvYieldEtaPP[ppEst][meth]           = NULL;
            systErrorUnCorrCollectionInvYieldEtaPP[ppEst][meth]     = NULL;
            systErrorInterCollectionInvYieldEtaPP[ppEst][meth]      = NULL;
            systErrorUnCorrInterCollectionInvYieldEtaPP[ppEst][meth]= NULL;
            graphRatioPi0IndCombFitStatPP[ppEst][meth]              = NULL;
            graphRatioPi0IndCombFitSysPP[ppEst][meth]               = NULL;
            graphRatioEtaIndCombFitStatPP[ppEst][meth]              = NULL;
            graphRatioEtaIndCombFitSysPP[ppEst][meth]               = NULL;
        }
    }

    TFile* fileInterpolation                        = new TFile(fileNameInterpolation.Data());
    for (Int_t ppEst = 0; ppEst < 2; ppEst++){
        for (Int_t meth = 0; meth< 11; meth++){
            statErrorCollectionPi0PP[ppEst][meth]                 = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionStatErr%s_Pi0_5.023TeV%s",
                                                                                                                    nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
            systErrorCollectionPi0PP[ppEst][meth]                 = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionSystErr%s_Pi0_5.023TeV%s",
                                                                                                                    nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
            systErrorUnCorrCollectionPi0PP[ppEst][meth]           = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionUnCorrSystErr%s_Pi0_5.023TeV%s",
                                                                                                                    nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
            systErrorInterCollectionPi0PP[ppEst][meth]            = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionInterpolSystErr%s_Pi0_5.023TeV%s",
                                                                                                                    nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
            if (statErrorCollectionPi0PP[ppEst][meth] && systErrorCollectionPi0PP[ppEst][meth] && systErrorUnCorrCollectionPi0PP[ppEst][meth] && systErrorInterCollectionPi0PP[ppEst][meth]){
                haveRefPPPi0[ppEst][meth]                         = kTRUE;
                systErrorUnCorrInterCollectionPi0PP[ppEst][meth]  = (TGraphAsymmErrors*)systErrorInterCollectionPi0PP[ppEst][meth]->Clone(Form("graphPP%sPi0InterUncorrSyst%s",
                                                                                                                    nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
                systErrorUnCorrInterCollectionPi0PP[ppEst][meth]  = AddErrorsQuadraticallyTGraph(systErrorInterCollectionPi0PP[ppEst][meth], systErrorUnCorrCollectionPi0PP[ppEst][meth]);

                statErrorCollectionInvYieldPi0PP[ppEst][meth]           = (TGraphAsymmErrors*)statErrorCollectionPi0PP[ppEst][meth]->Clone(Form("graphInvYieldStatErr%s_Pi0_5.023TeV%s",
                                                                                                                        nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
                systErrorCollectionInvYieldPi0PP[ppEst][meth]           = (TGraphAsymmErrors*)systErrorCollectionPi0PP[ppEst][meth]->Clone(Form("graphInvYieldSystErr%s_Pi0_5.023TeV%s",
                                                                                                                        nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
                systErrorUnCorrCollectionInvYieldPi0PP[ppEst][meth]     = (TGraphAsymmErrors*)systErrorUnCorrCollectionPi0PP[ppEst][meth]->Clone(Form("graphInvYieldUnCorrSystErr%s_Pi0_5.023TeV%s",
                                                                                                                        nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
                systErrorInterCollectionInvYieldPi0PP[ppEst][meth]      = (TGraphAsymmErrors*)systErrorInterCollectionPi0PP[ppEst][meth]->Clone(Form("graphInvYieldInterpolSystErr%s_Pi0_5.023TeV%s",
                                                                                                                        nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
                systErrorUnCorrInterCollectionInvYieldPi0PP[ppEst][meth]= (TGraphAsymmErrors*)systErrorUnCorrInterCollectionPi0PP[ppEst][meth]->Clone(Form("graphInvYieldInterpolSystErr%s_Pi0_5.023TeV%s",
                                                                                                                        nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));

                statErrorCollectionInvYieldPi0PP[ppEst][meth]           = ScaleGraph(statErrorCollectionInvYieldPi0PP[ppEst][meth],1/(xSection5TeV*recalcBarn));
                systErrorCollectionInvYieldPi0PP[ppEst][meth]           = ScaleGraph(systErrorCollectionInvYieldPi0PP[ppEst][meth],1/(xSection5TeV*recalcBarn));
                systErrorUnCorrCollectionInvYieldPi0PP[ppEst][meth]     = ScaleGraph(systErrorUnCorrCollectionInvYieldPi0PP[ppEst][meth],1/(xSection5TeV*recalcBarn));
                systErrorInterCollectionInvYieldPi0PP[ppEst][meth]      = ScaleGraph(systErrorInterCollectionInvYieldPi0PP[ppEst][meth],1/(xSection5TeV*recalcBarn));
                systErrorUnCorrInterCollectionInvYieldPi0PP[ppEst][meth]= ScaleGraph(systErrorUnCorrInterCollectionInvYieldPi0PP[ppEst][meth],1/(xSection5TeV*recalcBarn));

            }
            if (haveRefPPPi0[ppEst][meth])
                cout << "found pi0 pp reference for " << nameMeasGlobalLabel[meth].Data() << endl;



            statErrorCollectionEtaPP[ppEst][meth]                 = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionStatErr%s_Eta_5.023TeV%s",
                                                                                                                    nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
            systErrorCollectionEtaPP[ppEst][meth]                 = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionSystErr%s_Eta_5.023TeV%s",
                                                                                                                    nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
            systErrorUnCorrCollectionEtaPP[ppEst][meth]           = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionUnCorrSystErr%s_Eta_5.023TeV%s",
                                                                                                                    nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
            systErrorInterCollectionEtaPP[ppEst][meth]            = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionInterpolSystErr%s_Eta_5.023TeV%s",
                                                                                                                    nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
            if (statErrorCollectionEtaPP[ppEst][meth] && systErrorCollectionEtaPP[ppEst][meth] && systErrorUnCorrCollectionEtaPP[ppEst][meth] && systErrorInterCollectionEtaPP[ppEst][meth]){
                haveRefPPEta[ppEst][meth]                         = kTRUE;
                systErrorUnCorrInterCollectionEtaPP[ppEst][meth]        = (TGraphAsymmErrors*)systErrorInterCollectionEtaPP[ppEst][meth]->Clone(Form("graphPP%sEtaInterUncorrSyst%s",
                                                                                                                                               nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
                systErrorUnCorrInterCollectionEtaPP[ppEst][meth]        = AddErrorsQuadraticallyTGraph(systErrorInterCollectionEtaPP[ppEst][meth], systErrorUnCorrCollectionEtaPP[ppEst][meth]);

                statErrorCollectionInvYieldEtaPP[ppEst][meth]           = (TGraphAsymmErrors*)statErrorCollectionEtaPP[ppEst][meth]->Clone(Form("graphInvYieldStatErr%s_Eta_5.023TeV%s",
                                                                                                                            nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
                systErrorCollectionInvYieldEtaPP[ppEst][meth]           = (TGraphAsymmErrors*)systErrorCollectionEtaPP[ppEst][meth]->Clone(Form("graphInvYieldSystErr%s_Eta_5.023TeV%s",
                                                                                                                            nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
                systErrorUnCorrCollectionInvYieldEtaPP[ppEst][meth]     = (TGraphAsymmErrors*)systErrorUnCorrCollectionEtaPP[ppEst][meth]->Clone(Form("graphInvYieldUnCorrSystErr%s_Eta_5.023TeV%s",
                                                                                                                            nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
                systErrorInterCollectionInvYieldEtaPP[ppEst][meth]      = (TGraphAsymmErrors*)systErrorInterCollectionEtaPP[ppEst][meth]->Clone(Form("graphInvYieldInterpolSystErr%s_Eta_5.023TeV%s",
                                                                                                                            nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));
                systErrorUnCorrInterCollectionInvYieldEtaPP[ppEst][meth]= (TGraphAsymmErrors*)systErrorUnCorrInterCollectionEtaPP[ppEst][meth]->Clone(Form("graphInvYieldInterpolSystErr%s_Eta_5.023TeV%s",
                                                                                                                            nameMeasGlobalLabel[meth].Data(), addNamePP[ppEst].Data()));

                statErrorCollectionInvYieldEtaPP[ppEst][meth]           = ScaleGraph(statErrorCollectionInvYieldEtaPP[ppEst][meth],1/(xSection5TeV*recalcBarn));
                systErrorCollectionInvYieldEtaPP[ppEst][meth]           = ScaleGraph(systErrorCollectionInvYieldEtaPP[ppEst][meth],1/(xSection5TeV*recalcBarn));
                systErrorUnCorrCollectionInvYieldEtaPP[ppEst][meth]     = ScaleGraph(systErrorUnCorrCollectionInvYieldEtaPP[ppEst][meth],1/(xSection5TeV*recalcBarn));
                systErrorInterCollectionInvYieldEtaPP[ppEst][meth]      = ScaleGraph(systErrorInterCollectionInvYieldEtaPP[ppEst][meth],1/(xSection5TeV*recalcBarn));
                systErrorUnCorrInterCollectionInvYieldEtaPP[ppEst][meth]= ScaleGraph(systErrorUnCorrInterCollectionInvYieldEtaPP[ppEst][meth],1/(xSection5TeV*recalcBarn));
            }
            if (haveRefPPEta[ppEst][meth])
                cout << "found eta pp reference for " << nameMeasGlobalLabel[meth].Data() << endl;
        }

        graphPPCombPi0Stat[ppEst]                       = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionStatErrComb_Pi0_5.023TeV%s", addNamePP[ppEst].Data()));
        graphPPCombPi0UncorrSys[ppEst]                  = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionUnCorrSystErrComb_Pi0_5.023TeV%s", addNamePP[ppEst].Data()));
        graphPPCombPi0FullSys[ppEst]                    = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionSystErrComb_Pi0_5.023TeV%s", addNamePP[ppEst].Data()));
        graphPPCombPi0InterSys[ppEst]                   = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionInterpolSystErrComb_Pi0_5.023TeV%s", addNamePP[ppEst].Data()));

        graphPPInvYieldCombPi0Stat[ppEst]               = (TGraphAsymmErrors*)graphPPCombPi0Stat[ppEst]->Clone(Form("graphInvYieldStatErrComb_Pi0_5.023TeV%s", addNamePP[ppEst].Data()));
        graphPPInvYieldCombPi0Sys[ppEst]                = (TGraphAsymmErrors*)graphPPCombPi0FullSys[ppEst]->Clone(Form("graphInvYieldSystErrComb_Pi0_5.023TeV%s", addNamePP[ppEst].Data()));
        graphPPInvYieldCombPi0Stat[ppEst]               = ScaleGraph(graphPPCombPi0Stat[ppEst],1/(xSection5TeV*recalcBarn));
        graphPPInvYieldCombPi0Sys[ppEst]                = ScaleGraph(graphPPInvYieldCombPi0Sys[ppEst],1/(xSection5TeV*recalcBarn));
        graphPPInvYieldCombPi0StatW0XErr[ppEst]         = (TGraphAsymmErrors*)graphPPInvYieldCombPi0Stat[ppEst]->Clone(Form("graphInvYieldStatErrCombW0XErr_Pi0_5.023TeV%s", addNamePP[ppEst].Data()));
        ProduceGraphAsymmWithoutXErrors(graphPPInvYieldCombPi0StatW0XErr[ppEst]);


        graphPPCombEtaStat[ppEst]                       = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionStatErrComb_Eta_5.023TeV%s", addNamePP[ppEst].Data()));
        graphPPCombEtaUncorrSys[ppEst]                  = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionUnCorrSystErrComb_Eta_5.023TeV%s", addNamePP[ppEst].Data()));
        graphPPCombEtaFullSys[ppEst]                    = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionSystErrComb_Eta_5.023TeV%s", addNamePP[ppEst].Data()));
        graphPPCombEtaInterSys[ppEst]                   = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionInterpolSystErrComb_Eta_5.023TeV%s", addNamePP[ppEst].Data()));

        graphPPInvYieldCombEtaStat[ppEst]               = (TGraphAsymmErrors*)graphPPCombEtaStat[ppEst]->Clone(Form("graphInvYieldStatErrComb_Eta_5.023TeV%s", addNamePP[ppEst].Data()));
        graphPPInvYieldCombEtaSys[ppEst]                = (TGraphAsymmErrors*)graphPPCombEtaFullSys[ppEst]->Clone(Form("graphInvYieldSystErrComb_Eta_5.023TeV%s", addNamePP[ppEst].Data()));
        graphPPInvYieldCombEtaStat[ppEst]               = ScaleGraph(graphPPCombEtaStat[ppEst],1/(xSection5TeV*recalcBarn));
        graphPPInvYieldCombEtaSys[ppEst]                = ScaleGraph(graphPPInvYieldCombEtaSys[ppEst],1/(xSection5TeV*recalcBarn));
        graphPPInvYieldCombEtaStatW0XErr[ppEst]         = (TGraphAsymmErrors*)graphPPInvYieldCombEtaStat[ppEst]->Clone(Form("graphInvYieldStatErrCombW0XErr_Eta_5.023TeV%s", addNamePP[ppEst].Data()));
        ProduceGraphAsymmWithoutXErrors(graphPPInvYieldCombEtaStatW0XErr[ppEst]);
    }

    TFile* fileOtherParticle                            = new TFile(fileNameOtherParticleInput.Data());
    TDirectory* directoryOtherPart[17]                  = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphDRpAStat[17]                = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphDRpASys[17]                 = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphDRpATot[17]                 = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphChHadRpAStat[17]            = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphChHadRpASys[17]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphChHadRpATot[17]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphChKToPiStat[17]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphChKToPiSys[17]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphChKToPiTot[17]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    for (Int_t cent = 0; cent < 17; cent++){

//         if (centArray[cent].BeginsWith("0-100%")) continue;
        directoryOtherPart[cent]                               = (TDirectory*)fileOtherParticle->Get(Form("pPb5TeV_%s%s",centArray[cent].Data(), addCentString[cent].Data()));
        if (directoryOtherPart[cent]){
            graphDRpAStat[cent]                                 = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("DMesonQpPbStatW0XErr");
            graphDRpASys[cent]                                  = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("DMesonQpPbSys");
            graphDRpATot[cent]                                  = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("DMesonQpPbTot");
            graphChHadRpAStat[cent]                             = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("ChargedHadQpPbStatW0XErr");
            graphChHadRpASys[cent]                              = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("ChargedHadQpPbSys");
            graphChHadRpATot[cent]                              = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("ChargedHadQpPbTot");
            graphChKToPiStat[cent]                              = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("graphChargedKaonToPionPubStatW0XErr");
            graphChKToPiSys[cent]                               = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("graphChargedKaonToPionPubSys");
            graphChKToPiTot[cent]                               = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("graphChargedKaonToPionPubTot");
        }
    }

    TFile* fileOlderMB                                  = new TFile(fileNameOlderPPbMB);
    TGraphAsymmErrors* graphCombPi0OlderMBStat          = (TGraphAsymmErrors*)fileOlderMB->Get("Pi0pPb/graphInvYieldPi0CombpPb5023GeVStatErr_NSD");
    TGraphAsymmErrors* graphCombPi0OlderMBSys           = (TGraphAsymmErrors*)fileOlderMB->Get("Pi0pPb/graphInvYieldPi0CombpPb5023GeVSystErr_NSD");
    TGraphAsymmErrors* graphCombPi0OlderMBTot           = (TGraphAsymmErrors*)fileOlderMB->Get("Pi0pPb/graphInvYieldPi0CombpPb5023GeVTotErr_NSD");
    TGraphAsymmErrors* graphCombEtaOlderMBStat          = (TGraphAsymmErrors*)fileOlderMB->Get("EtapPb/graphInvYieldEtaCombpPb5023GeVStatErr_NSD");
    TGraphAsymmErrors* graphCombEtaOlderMBSys           = (TGraphAsymmErrors*)fileOlderMB->Get("EtapPb/graphInvYieldEtaCombpPb5023GeVSystErr_NSD");
    TGraphAsymmErrors* graphCombEtaOlderMBTot           = (TGraphAsymmErrors*)fileOlderMB->Get("EtapPb/graphInvYieldEtaCombpPb5023GeVTotErr_NSD");
    TGraphAsymmErrors* graphCombPi0RpPbOlderMBStat      = (TGraphAsymmErrors*)fileOlderMB->Get("Pi0RpPb/CombinedPi0RpPbStatErr");
    TGraphAsymmErrors* graphCombPi0RpPbOlderMBSys       = (TGraphAsymmErrors*)fileOlderMB->Get("Pi0RpPb/CombinedPi0RpPbSystErr");
    TGraphAsymmErrors* graphCombPi0RpPbOlderMBTot       = AddErrorsOfGraphsQuadratically(graphCombPi0RpPbOlderMBStat,graphCombPi0RpPbOlderMBSys);
    TGraphAsymmErrors* graphCombEtaRpPbOlderMBStat      = (TGraphAsymmErrors*)fileOlderMB->Get("EtaRpPb/CombinedEtaRpPbStatErr");
    TGraphAsymmErrors* graphCombEtaRpPbOlderMBSys       = (TGraphAsymmErrors*)fileOlderMB->Get("EtaRpPb/CombinedEtaRpPbSystErr");
    TGraphAsymmErrors* graphCombEtaRpPbOlderMBTot       = AddErrorsOfGraphsQuadratically(graphCombEtaRpPbOlderMBStat,graphCombEtaRpPbOlderMBSys);
    TGraphAsymmErrors* graphCombEtaToPi0OlderMBStat     = (TGraphAsymmErrors*)fileOlderMB->Get("EtapPb/EtaPi0pPbStat");
    TGraphAsymmErrors* graphCombEtaToPi0OlderMBSys      = (TGraphAsymmErrors*)fileOlderMB->Get("EtapPb/EtaPi0pPbSyst");
    TGraphAsymmErrors* graphCombEtaToPi0OlderMBTot      = AddErrorsOfGraphsQuadratically(graphCombEtaToPi0OlderMBStat,graphCombEtaToPi0OlderMBSys);
    TF1* fitCombPi0OlderMB                              = (TF1*)fileOlderMB->Get("Pi0pPb/FitCombPi0");
    TF1* fitCombEtaOlderMB                              = (TF1*)fileOlderMB->Get("EtapPb/FitCombEta");

    TGraphAsymmErrors* graphCombPi0OlderMBRelErrStat   = (TGraphAsymmErrors*)graphCombPi0OlderMBStat->Clone("relativeStatErrorPi0OldMB");
    RemoveZerosAtBeginningAndEndFromGraph(graphCombPi0OlderMBRelErrStat);
    graphCombPi0OlderMBRelErrStat   = CalculateRelErrUpAsymmGraph( graphCombPi0OlderMBRelErrStat,"relativeStatErrorPi0OldMB");
    TGraphAsymmErrors* graphCombPi0OlderMBRelErrSys   = (TGraphAsymmErrors*)graphCombPi0OlderMBSys->Clone("relativeSysErrorPi0OldMB");
    RemoveZerosAtBeginningAndEndFromGraph(graphCombPi0OlderMBRelErrSys);
    graphCombPi0OlderMBRelErrSys   = CalculateRelErrUpAsymmGraph( graphCombPi0OlderMBRelErrSys,"relativeSysErrorPi0OldMB");
    TGraphAsymmErrors* graphCombPi0OlderMBRelErrTot   = (TGraphAsymmErrors*)graphCombPi0OlderMBTot->Clone("relativeTotErrorPi0OldMB");
    RemoveZerosAtBeginningAndEndFromGraph(graphCombPi0OlderMBRelErrTot);
    graphCombPi0OlderMBRelErrTot   = CalculateRelErrUpAsymmGraph( graphCombPi0OlderMBRelErrTot,"relativeTotErrorPi0OldMB");
    TGraphAsymmErrors* graphCombEtaOlderMBRelErrStat   = (TGraphAsymmErrors*)graphCombEtaOlderMBStat->Clone("relativeStatErrorEtaOldMB");
    RemoveZerosAtBeginningAndEndFromGraph(graphCombEtaOlderMBRelErrStat);
    graphCombEtaOlderMBRelErrStat   = CalculateRelErrUpAsymmGraph( graphCombEtaOlderMBRelErrStat,"relativeStatErrorEtaOldMB");
    TGraphAsymmErrors* graphCombEtaOlderMBRelErrSys   = (TGraphAsymmErrors*)graphCombEtaOlderMBSys->Clone("relativeSysErrorEtaOldMB");
    RemoveZerosAtBeginningAndEndFromGraph(graphCombEtaOlderMBRelErrSys);
    graphCombEtaOlderMBRelErrSys   = CalculateRelErrUpAsymmGraph( graphCombEtaOlderMBRelErrSys,"relativeSysErrorEtaOldMB");
    TGraphAsymmErrors* graphCombEtaOlderMBRelErrTot   = (TGraphAsymmErrors*)graphCombEtaOlderMBTot->Clone("relativeTotErrorEtaOldMB");
    RemoveZerosAtBeginningAndEndFromGraph(graphCombEtaOlderMBRelErrTot);
    graphCombEtaOlderMBRelErrTot   = CalculateRelErrUpAsymmGraph( graphCombEtaOlderMBRelErrTot,"relativeTotErrorEtaOldMB");
    TGraphAsymmErrors* graphCombPi0RpPbOlderMBRelErrStat   = (TGraphAsymmErrors*)graphCombPi0RpPbOlderMBStat->Clone("relativeStatErrorPi0RpPbOldMB");
    RemoveZerosAtBeginningAndEndFromGraph(graphCombPi0RpPbOlderMBRelErrStat);
    graphCombPi0RpPbOlderMBRelErrStat   = CalculateRelErrUpAsymmGraph( graphCombPi0RpPbOlderMBRelErrStat,"relativeStatErrorPi0RpPbOldMB");
    TGraphAsymmErrors* graphCombPi0RpPbOlderMBRelErrSys   = (TGraphAsymmErrors*)graphCombPi0RpPbOlderMBSys->Clone("relativeSysErrorPi0RpPbOldMB");
    RemoveZerosAtBeginningAndEndFromGraph(graphCombPi0RpPbOlderMBRelErrSys);
    graphCombPi0RpPbOlderMBRelErrSys   = CalculateRelErrUpAsymmGraph( graphCombPi0RpPbOlderMBRelErrSys,"relativeSysErrorPi0RpPbOldMB");
    TGraphAsymmErrors* graphCombPi0RpPbOlderMBRelErrTot   = (TGraphAsymmErrors*)graphCombPi0RpPbOlderMBTot->Clone("relativeTotErrorPi0RpPbOldMB");
    RemoveZerosAtBeginningAndEndFromGraph(graphCombPi0RpPbOlderMBRelErrTot);
    graphCombPi0RpPbOlderMBRelErrTot   = CalculateRelErrUpAsymmGraph( graphCombPi0RpPbOlderMBRelErrTot,"relativeTotErrorPi0RpPbOldMB");
    TGraphAsymmErrors* graphCombEtaRpPbOlderMBRelErrStat   = (TGraphAsymmErrors*)graphCombEtaRpPbOlderMBStat->Clone("relativeStatErrorEtaRpPbOldMB");
    RemoveZerosAtBeginningAndEndFromGraph(graphCombEtaRpPbOlderMBRelErrStat);
    graphCombEtaRpPbOlderMBRelErrStat   = CalculateRelErrUpAsymmGraph( graphCombEtaRpPbOlderMBRelErrStat,"relativeStatErrorEtaRpPbOldMB");
    TGraphAsymmErrors* graphCombEtaRpPbOlderMBRelErrSys   = (TGraphAsymmErrors*)graphCombEtaRpPbOlderMBSys->Clone("relativeSysErrorEtaRpPbOldMB");
    RemoveZerosAtBeginningAndEndFromGraph(graphCombEtaRpPbOlderMBRelErrSys);
    graphCombEtaRpPbOlderMBRelErrSys   = CalculateRelErrUpAsymmGraph( graphCombEtaRpPbOlderMBRelErrSys,"relativeSysErrorEtaRpPbOldMB");
    TGraphAsymmErrors* graphCombEtaRpPbOlderMBRelErrTot   = (TGraphAsymmErrors*)graphCombEtaRpPbOlderMBTot->Clone("relativeTotErrorEtaRpPbOldMB");
    RemoveZerosAtBeginningAndEndFromGraph(graphCombEtaRpPbOlderMBRelErrTot);
    graphCombEtaRpPbOlderMBRelErrTot   = CalculateRelErrUpAsymmGraph( graphCombEtaRpPbOlderMBRelErrTot,"relativeTotErrorEtaRpPbOldMB");

    // *******************************************************************************************************
    // ******************************** Setup for pp plotting of reference ***********************************
    // *******************************************************************************************************
    Int_t textSizeLabelsPixel               = 52;
    TLatex *labelRatioToFitEnergyPP   = new TLatex(0.95, 0.92, collisionSystempp.Data());
    SetStyleTLatex( labelRatioToFitEnergyPP, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitEnergyPP->Draw();
    TLatex *labelRatioToFitALICEPP    = new TLatex(0.95, 0.86, textALICE.Data());
    SetStyleTLatex( labelRatioToFitALICEPP, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    labelRatioToFitALICEPP->Draw();
    TLatex *labelRatioToFitPi0PP      = new TLatex(0.12, 0.92, "#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRatioToFitPi0PP, textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
    labelRatioToFitPi0PP->Draw();
    TLatex *labelRatioToFitEtaPP      = new TLatex(0.12, 0.92, "#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRatioToFitEtaPP, textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
    labelRatioToFitEtaPP->Draw();

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
    histo2DPi0RatioToCombFitPP->GetYaxis()->SetRangeUser(0.2,1.82);

    TH2F * histo2DEtaRatioToCombFitPP;
    histo2DEtaRatioToCombFitPP               = new TH2F("histo2DEtaRatioToCombFitPP","histo2DEtaRatioToCombFitPP",1000,minPtEtaPlotting, maxPtEtaPlotting,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DEtaRatioToCombFitPP, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                              0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.65, 510, 505);
    histo2DEtaRatioToCombFitPP->GetXaxis()->SetMoreLogLabels();
    histo2DEtaRatioToCombFitPP->GetXaxis()->SetLabelOffset(-0.01);
    histo2DEtaRatioToCombFitPP->GetYaxis()->SetRangeUser(0.2,1.82);

    //****************************** Definition of the Legend ******************************************
    //**************** Row def ************************
    Double_t rowsLegendOnlyPi0RatioPP[4]          = {0.30, 0.25, 0.20, 0.15};
    Double_t rowsLegendOnlyPi0RatioPPAbs[4]       = {0.92, 0.47, 0.38, 0.29 };
    Double_t columnsLegendOnlyPi0RatioPP[6]       = {0.14, 0.26, 0.35, 0.48, 0.7, 0.8};
    Double_t columnsLegendOnlyPi0RatioPPAbs[6]    = {0.215, 0.8, 1.2, 2, 11., 17.0};
    Double_t lengthBoxPP                          = 0.3;
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
    //****************************** Definition of the Legend ******************************************
    //**************** Row def ************************
    Double_t rowsLegendOnlyEtaRatioPP[4]          = {0.30, 0.25, 0.20, 0.15};
    Double_t rowsLegendOnlyEtaRatioPPAbs[4]       = {0.92, 0.47, 0.38, 0.29 };
    Double_t columnsLegendOnlyEtaRatioPP[6]       = {0.14, 0.26, 0.35, 0.48, 0.7, 0.8};
    Double_t columnsLegendOnlyEtaRatioPPAbs[6]    = {0.215, 1.2, 1.65, 2, 10.2, 14.2};
    Double_t lengthBoxPPEta                       = 0.3;
    //****************** first Column **************************************************
    TLatex *textPCMOnlyRatioEtaPP                 = new TLatex(columnsLegendOnlyEtaRatioPP[0],rowsLegendOnlyEtaRatioPP[1],nameMeasGlobalLabel[0]);
    SetStyleTLatex( textPCMOnlyRatioEtaPP, textSizeLabelsPixel,4, 1, 43);
    textPCMOnlyRatioEtaPP->Draw();
    TLatex *textPHOSOnlyRatioEtaPP                = new TLatex(columnsLegendOnlyEtaRatioPP[0],rowsLegendOnlyEtaRatioPP[2],nameMeasGlobalLabel[1]);
    SetStyleTLatex( textPHOSOnlyRatioEtaPP, textSizeLabelsPixel,4, 1, 43);
    textPHOSOnlyRatioEtaPP->Draw();
    TLatex *textEMCALOnlyRatioEtaPP               = new TLatex(columnsLegendOnlyEtaRatioPP[0],rowsLegendOnlyEtaRatioPP[3],nameMeasGlobalLabel[2]);
    SetStyleTLatex( textEMCALOnlyRatioEtaPP, textSizeLabelsPixel,4, 1, 43);
    textEMCALOnlyRatioEtaPP->Draw();
    TLatex *textPCMEMCALOnlyRatioEtaPP            = new TLatex(columnsLegendOnlyEtaRatioPP[3],rowsLegendOnlyEtaRatioPP[1],nameMeasGlobalLabel[4]);
    SetStyleTLatex( textPCMEMCALOnlyRatioEtaPP, textSizeLabelsPixel,4, 1, 43);
    textPCMEMCALOnlyRatioEtaPP->Draw();
    TLatex *textPCMPHOSOnlyRatioEtaPP            = new TLatex(columnsLegendOnlyEtaRatioPP[3],rowsLegendOnlyEtaRatioPP[2],nameMeasGlobalLabel[3]);
    SetStyleTLatex( textPCMPHOSOnlyRatioEtaPP, textSizeLabelsPixel,4, 1, 43);
    textPCMPHOSOnlyRatioEtaPP->Draw();

    //****************** second Column *************************************************
    TLatex *textStatOnlyRatioEtaPP                = new TLatex(columnsLegendOnlyEtaRatioPP[1],rowsLegendOnlyEtaRatioPP[0] ,"stat");
    SetStyleTLatex( textStatOnlyRatioEtaPP, textSizeLabelsPixel,4, 1, 43);
    textStatOnlyRatioEtaPP->Draw();
    TLatex *textSysOnlyRatioEtaPP                 = new TLatex(columnsLegendOnlyEtaRatioPP[2] ,rowsLegendOnlyEtaRatioPP[0],"syst");
    SetStyleTLatex( textSysOnlyRatioEtaPP, textSizeLabelsPixel,4, 1, 43);
    textSysOnlyRatioEtaPP->Draw();
    TLatex *textStatOnlyRatioEtaPP2               = new TLatex(columnsLegendOnlyEtaRatioPP[4],rowsLegendOnlyEtaRatioPP[0] ,"stat");
    SetStyleTLatex( textStatOnlyRatioEtaPP2, textSizeLabelsPixel,4, 1, 43);
    textStatOnlyRatioEtaPP2->Draw();
    TLatex *textSysOnlyRatioEtaPP2                = new TLatex(columnsLegendOnlyEtaRatioPP[5] ,rowsLegendOnlyEtaRatioPP[0],"syst");
    SetStyleTLatex( textSysOnlyRatioEtaPP2, textSizeLabelsPixel,4, 1, 43);
    textSysOnlyRatioEtaPP2->Draw();


    for (Int_t ppEst = 0; ppEst < 2; ppEst++){
        // *******************************************************************************************************
        // ************************** Comparison of different pp interpolated pi0 spectra ************************
        // *******************************************************************************************************
        // fitting spectrum with intial parameters
        // Two component model fit from Bylinkin
        TF1* fitTCMDecomposedLPi0PP                             = FitObject("tcmlow",Form("twoCompModelPi0_DecL_PP%s", addNamePP[ppEst].Data()), "Pi0", NULL, 0.3, 2.);
        TF1* fitTCMDecomposedHPi0PP                             = FitObject("tcmhigh",Form("twoCompModelPi0_DecH_PP%s", addNamePP[ppEst].Data()), "Pi0", NULL, 4, 50.);
        fitTCMDecomposedLPi0PP->SetParameters(graphPPCombPi0Stat[ppEst]->GetY()[2],0.3);
        graphPPCombPi0Stat[ppEst]->Fit(fitTCMDecomposedLPi0PP,"QNRMEX0+","",0.3,0.8);
        graphPPCombPi0Stat[ppEst]->Fit(fitTCMDecomposedHPi0PP,"QNRMEX0+","",3, 20.);
        fitTCMDecomposedHPi0PP->SetParameters(graphPPCombPi0Stat[ppEst]->GetY()[2],0.8, 2);

        cout << WriteParameterToFile(fitTCMDecomposedLPi0PP)<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedLPi0PP)<< endl;
        cout << WriteParameterToFile(fitTCMDecomposedHPi0PP)<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedHPi0PP)<< endl;

        Double_t paramTCMPi0NewPP[5]                            = { fitTCMDecomposedLPi0PP->GetParameter(0),fitTCMDecomposedLPi0PP->GetParameter(1),
                                                                    fitTCMDecomposedHPi0PP->GetParameter(0),fitTCMDecomposedHPi0PP->GetParameter(1),fitTCMDecomposedHPi0PP->GetParameter(2)};

        //Two component model from Bylinkin
        fitTCMInvXSetionPi0PP[ppEst]                            = FitObject("tcm",Form("fitTCMInvXSetionPi0PP5TeV%s",addNamePP[ppEst].Data()),"Pi0",graphPPCombPi0Stat[ppEst], 0.3, 20. ,paramTCMPi0NewPP,"QNRMEX0+","", kFALSE);
        cout << "fitting pp pi0 spectra" << endl;
        cout << WriteParameterToFile(fitTCMInvXSetionPi0PP[ppEst])<< endl;
        fitTsallisInvXSectionPi0PP[ppEst]                       = FitObject("l",Form("fitTsallisInvXSectionPi0PP5TeV%s",addNamePP[ppEst].Data()),"Pi0",graphPPCombPi0Stat[ppEst], 0.3, 20. ,
                                                                              NULL,"QNRMEX0+","", kFALSE);
        cout << "fitting pp pi0 spectra" << endl;
        cout << WriteParameterToFile(fitTsallisInvXSectionPi0PP[ppEst])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTsallisInvXSectionPi0PP[ppEst])<< endl;

        fileFitsOutput <<  WriteParameterToFile(fitTCMInvXSetionPi0PP[ppEst])<< endl;
        fitPPTCMInvYieldPi0[ppEst]                              = FitObject("tcm",Form("fitTCMInvYieldPi0PP5TeV%s",addNamePP[ppEst].Data()),"Pi0",NULL, 0.3, 20. ,NULL,"QNRMEX0+","", kFALSE);
        cout << "fitting pp pi0 spectra" << endl;
        fitPPTCMInvYieldPi0[ppEst]->SetParameter(0,graphPPInvYieldCombPi0Stat[ppEst]->GetY()[2]);
        fitPPTCMInvYieldPi0[ppEst]->SetParameter(1,fitTCMInvXSetionPi0PP[ppEst]->GetParameter(1));
        fitPPTCMInvYieldPi0[ppEst]->SetParameter(2,graphPPInvYieldCombPi0Stat[ppEst]->GetY()[2]);
        fitPPTCMInvYieldPi0[ppEst]->SetParameter(3,fitTCMInvXSetionPi0PP[ppEst]->GetParameter(3));
        fitPPTCMInvYieldPi0[ppEst]->SetParameter(4,fitTCMInvXSetionPi0PP[ppEst]->GetParameter(4));
        graphPPInvYieldCombPi0Stat[ppEst]->Fit(fitPPTCMInvYieldPi0[ppEst],"QNRMEX0+","", 0.3, 20.);
        cout << WriteParameterToFile(fitPPTCMInvYieldPi0[ppEst])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitPPTCMInvYieldPi0[ppEst])<< endl;
        fitPPTCMInvYieldPi0[ppEst]->SetRange(0.3, 20.);

    //     TGraphAsymmErrors* graphRatioPPPi0NLODSS14          = (TGraphAsymmErrors*)graphNLODSS14Pi0PP->Clone("graphRatioPPPi0NLODSS14ToFit");
    //     graphRatioPPPi0NLODSS14                             = CalculateGraphErrRatioToFit(graphRatioPPPi0NLODSS14, fitPPTCMInvYieldPi0[ppEst]);
    //     TGraph* graphRatioPPPi0NLODSS14Center            mei   = (TGraph*)graphNLODSS14Pi0PPCenter->Clone("graphRatioPPPi0NLODSS14CenterToFit");
    //     graphRatioPPPi0NLODSS14Center                       = CalculateGraphRatioToFit(graphRatioPPPi0NLODSS14Center, fitPPTCMInvYieldPi0[ppEst]);


        // *************************************************************************************************************
        // Shift graphs in Y direction as well if desired
        // *************************************************************************************************************

        graphRatioPi0CombToFitStatPP[ppEst]                         = (TGraphAsymmErrors*)graphPPCombPi0Stat[ppEst]->Clone();
        graphRatioPi0CombToFitStatPP[ppEst]                         = CalculateGraphErrRatioToFit(graphRatioPi0CombToFitStatPP[ppEst], fitTCMInvXSetionPi0PP[ppEst]);
        graphRatioPi0CombToFitSysPP[ppEst]                          = (TGraphAsymmErrors*)graphPPCombPi0FullSys[ppEst]->Clone();
        graphRatioPi0CombToFitSysPP[ppEst]                          = CalculateGraphErrRatioToFit(graphRatioPi0CombToFitSysPP[ppEst], fitTCMInvXSetionPi0PP[ppEst]);

        for (Int_t meth = 0; meth< 11; meth++){
            if (statErrorCollectionPi0PP[ppEst][meth]){
                graphRatioPi0IndCombFitStatPP[ppEst][meth]                = (TGraphAsymmErrors*)statErrorCollectionPi0PP[ppEst][meth]->Clone(Form("RatioPi0%sToCombFitStatPP", nameMeasGlobalLabel[meth].Data()));
                graphRatioPi0IndCombFitStatPP[ppEst][meth]                = CalculateGraphErrRatioToFit(graphRatioPi0IndCombFitStatPP[ppEst][meth], fitTCMInvXSetionPi0PP[ppEst]);
            }
            if (systErrorCollectionPi0PP[ppEst][meth]){
                graphRatioPi0IndCombFitSysPP[ppEst][meth]                 = (TGraphAsymmErrors*)systErrorCollectionPi0PP[ppEst][meth]->Clone(Form("RatioPi0%sToCombFitSystPP", nameMeasGlobalLabel[meth].Data()));
                graphRatioPi0IndCombFitSysPP[ppEst][meth]                 = CalculateGraphErrRatioToFit(graphRatioPi0IndCombFitSysPP[ppEst][meth], fitTCMInvXSetionPi0PP[ppEst]);
            }
        }

        // **********************************************************************************************************************
        // ******************************************* Plot Ratio of Comb to Fit for pp *****************************************
        // **********************************************************************************************************************
        canvasRatioToCombFitPP->cd();
        histo2DPi0RatioToCombFitPP->Draw("copy");

        ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombToFitStatPP[ppEst]);

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombToFitSysPP[ppEst], markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioPi0CombToFitSysPP[ppEst]->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombToFitStatPP[ppEst], markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioPi0CombToFitStatPP[ppEst]->Draw("p,same,z");

        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1., 1.,0.1, kGray+2);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.9, 0.9,0.1, kGray, 7);

        labelRatioToFitEnergyPP->Draw();
        labelRatioToFitALICEPP->Draw();
        labelRatioToFitPi0PP->Draw();

        canvasRatioToCombFitPP->SaveAs(Form("%s/Pi0_RatioOfCombToCombFit_PP5TeV%s.%s",outputDirSupportPaper.Data(), addNamePP[ppEst].Data(), suffix.Data()));

        // **********************************************************************************************************************
        // *******************************************Plot Ratio of Individual meas to Fit ******************************************
        // **********************************************************************************************************************

        canvasRatioToCombFitPP->cd();
        histo2DPi0RatioToCombFitPP->Draw("copy");

        for (Int_t meth = 10; meth > -1 ; meth--){
            if (graphRatioPi0IndCombFitSysPP[ppEst][meth]){
                DrawGammaSetMarkerTGraphAsym(graphRatioPi0IndCombFitSysPP[ppEst][meth], markerStyleDetMC[meth] ,markerSizeDetMC[meth]*0.5, colorDet[meth], colorDet[meth], widthLinesBoxes, kTRUE);
                graphRatioPi0IndCombFitSysPP[ppEst][meth]->Draw("E2same");
            }
            if (graphRatioPi0IndCombFitStatPP[ppEst][meth]){
                ProduceGraphAsymmWithoutXErrors(graphRatioPi0IndCombFitStatPP[ppEst][meth]);
                DrawGammaSetMarkerTGraphAsym(graphRatioPi0IndCombFitStatPP[ppEst][meth], markerStyleDetMC[meth] ,markerSizeDetMC[meth]*0.5, colorDet[meth], colorDet[meth]);
                graphRatioPi0IndCombFitStatPP[ppEst][meth]->Draw("p,same,z");
            }
        }
        if (graphRatioPi0IndCombFitStatPP[ppEst][4])graphRatioPi0IndCombFitStatPP[ppEst][4]->Draw("p,same,z");

        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1., 1.,0.5, kGray+2);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.9, 0.9,0.5, kGray, 7);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.2, 1.2,0.5, kGray, 9);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.8, 0.8,0.5, kGray, 9);

        labelRatioToFitEnergyPP->Draw();
        labelRatioToFitALICEPP->Draw();
        labelRatioToFitPi0PP->Draw();
        histo2DPi0RatioToCombFitPP->Draw("same,axis");

        //****************** first Column *************************************************
        textPCMOnlyRatioPi0PP->Draw();
        textPHOSOnlyRatioPi0PP->Draw();
        textEMCALOnlyRatioPi0PP->Draw();
        textPCMEMCALOnlyRatioPi0PP->Draw();
        textPCMPHOSOnlyRatioPi0PP->Draw();
        textDalitzOnlyRatioPi0PP->Draw();

        //****************** second Column *************************************************
        textStatOnlyRatioPi0PP->Draw();
        textSysOnlyRatioPi0PP->Draw();
        textStatOnlyRatioPi0PP2->Draw();
        textSysOnlyRatioPi0PP2->Draw();

        TMarker* markerPCMPi0OnlyRatioPi0PP             = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[ppEst][0],columnsLegendOnlyPi0RatioPP[1] ,rowsLegendOnlyPi0RatioPP[1],1);
        markerPCMPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[1] ,rowsLegendOnlyPi0RatioPPAbs[1]);
        TMarker* markerPHOSPi0OnlyRatioPi0PP            = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[ppEst][1], columnsLegendOnlyPi0RatioPP[1] ,rowsLegendOnlyPi0RatioPP[2],1);
        markerPHOSPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[1] ,rowsLegendOnlyPi0RatioPPAbs[2]);
        TMarker* markerEMCALPi0OnlyRatioPi0PP           = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[ppEst][2], columnsLegendOnlyPi0RatioPP[1] ,rowsLegendOnlyPi0RatioPP[3],1);
        markerEMCALPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[1] ,rowsLegendOnlyPi0RatioPPAbs[3]);
        TMarker* markerPCMEMCALPi0OnlyRatioPi0PP        = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[ppEst][4], columnsLegendOnlyPi0RatioPP[3] ,rowsLegendOnlyPi0RatioPP[1],1);
        markerPCMEMCALPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[4] ,rowsLegendOnlyPi0RatioPPAbs[1]);
        TMarker* markerPCMPHOSPi0OnlyRatioPi0PP         = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[ppEst][3], columnsLegendOnlyPi0RatioPP[3] ,rowsLegendOnlyPi0RatioPP[2],1);
        markerPCMPHOSPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[4] ,rowsLegendOnlyPi0RatioPPAbs[2]);

        TBox* boxPCMPi0OnlyRatioPi0PP                 = CreateBoxFromGraph(graphRatioPi0IndCombFitSysPP[ppEst][0], columnsLegendOnlyPi0RatioPPAbs[2]-0.1*lengthBoxPP , rowsLegendOnlyPi0RatioPPAbs[1]- heightBoxPP,
                                                                        columnsLegendOnlyPi0RatioPPAbs[2]+ 1.7*lengthBoxPP, rowsLegendOnlyPi0RatioPPAbs[1]+ heightBoxPP);
        boxPCMPi0OnlyRatioPi0PP->Draw("l");
        TBox* boxPHOSPi0OnlyRatioPi0PP                = CreateBoxFromGraph(graphRatioPi0IndCombFitSysPP[ppEst][1], columnsLegendOnlyPi0RatioPPAbs[2]-0.1*lengthBoxPP , rowsLegendOnlyPi0RatioPPAbs[2]- heightBoxPP,
                                                                        columnsLegendOnlyPi0RatioPPAbs[2]+ 1.7*lengthBoxPP, rowsLegendOnlyPi0RatioPPAbs[2]+ heightBoxPP);
        boxPHOSPi0OnlyRatioPi0PP->Draw("l");
        TBox* boxEMCALPi0OnlyRatioPi0PP               = CreateBoxFromGraph(graphRatioPi0IndCombFitSysPP[ppEst][2], columnsLegendOnlyPi0RatioPPAbs[2]-0.1*lengthBoxPP , rowsLegendOnlyPi0RatioPPAbs[3]- heightBoxPP,
                                                                        columnsLegendOnlyPi0RatioPPAbs[2]+ 1.7*lengthBoxPP, rowsLegendOnlyPi0RatioPPAbs[3]+ heightBoxPP);
        boxEMCALPi0OnlyRatioPi0PP->Draw("l");
        TBox* boxPCMEMCALPi0OnlyRatioPi0PP            = CreateBoxFromGraph(graphRatioPi0IndCombFitSysPP[ppEst][4], columnsLegendOnlyPi0RatioPPAbs[5]-0.5*lengthBoxPP , rowsLegendOnlyPi0RatioPPAbs[1]- heightBoxPP,
                                                                        columnsLegendOnlyPi0RatioPPAbs[5]+ 17.2*lengthBoxPP, rowsLegendOnlyPi0RatioPPAbs[1]+ heightBoxPP);
        boxPCMEMCALPi0OnlyRatioPi0PP->Draw("l");
        TBox* boxPCMPHOSPi0OnlyRatioPi0PP             = CreateBoxFromGraph(graphRatioPi0IndCombFitSysPP[ppEst][3], columnsLegendOnlyPi0RatioPPAbs[5]-0.5*lengthBoxPP , rowsLegendOnlyPi0RatioPPAbs[2]- heightBoxPP,
                                                                        columnsLegendOnlyPi0RatioPPAbs[5]+ 17.2*lengthBoxPP, rowsLegendOnlyPi0RatioPPAbs[2]+ heightBoxPP);
        boxPCMPHOSPi0OnlyRatioPi0PP->Draw("l");

        if(graphRatioPi0IndCombFitSysPP[ppEst][5]) {
            TMarker* markerDalitzPi0OnlyRatioPi0PP      = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[ppEst][5], columnsLegendOnlyPi0RatioPP[3] ,rowsLegendOnlyPi0RatioPP[3],1);
            markerDalitzPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[4] ,rowsLegendOnlyPi0RatioPPAbs[3]);
            TBox* boxDalitzPi0OnlyRatioPi0PP             = CreateBoxFromGraph(graphRatioPi0IndCombFitSysPP[ppEst][5], columnsLegendOnlyPi0RatioPPAbs[5]-0.5*lengthBoxPP , rowsLegendOnlyPi0RatioPPAbs[3]- heightBoxPP,
                                                                              columnsLegendOnlyPi0RatioPPAbs[5]+ 17.2*lengthBoxPP, rowsLegendOnlyPi0RatioPPAbs[3]+ heightBoxPP);
            boxDalitzPi0OnlyRatioPi0PP->Draw("l");
        }

        canvasRatioToCombFitPP->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFit_PP5TeV%s.%s",outputDirSupportPaper.Data(), addNamePP[ppEst].Data(),suffix.Data()));

        // *******************************************************************************************************
        // ************************** Comparison of different pp interpolated eta spectra ************************
        // *******************************************************************************************************
        // fitting spectrum with intial parameters
        // Two component model fit from Bylinkin
        TF1* fitTCMDecomposedLEtaPP                             = FitObject("tcmlow",Form("twoCompModelEta_DecL_PP%s", addNamePP[ppEst].Data()), "Eta", NULL, 0.3, 2.);
        TF1* fitTCMDecomposedHEtaPP                             = FitObject("tcmhigh",Form("twoCompModelEta_DecH_PP%s", addNamePP[ppEst].Data()), "Eta", NULL, 4, 50.);
        fitTCMDecomposedLEtaPP->SetParameters(graphPPCombEtaStat[ppEst]->GetY()[2],0.3);
        graphPPCombEtaStat[ppEst]->Fit(fitTCMDecomposedLEtaPP,"QNRMEX0+","",0.3,0.8);
        graphPPCombEtaStat[ppEst]->Fit(fitTCMDecomposedHEtaPP,"QNRMEX0+","",3, 20.);
        fitTCMDecomposedHEtaPP->SetParameters(graphPPCombEtaStat[ppEst]->GetY()[2],0.8, 2);

        cout << WriteParameterToFile(fitTCMDecomposedLEtaPP)<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedLEtaPP)<< endl;
        cout << WriteParameterToFile(fitTCMDecomposedHEtaPP)<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedHEtaPP)<< endl;

        Double_t paramTCMEtaNewPP[5]                            = { fitTCMDecomposedLEtaPP->GetParameter(0),fitTCMDecomposedLEtaPP->GetParameter(1),
                                                                    fitTCMDecomposedHEtaPP->GetParameter(0),fitTCMDecomposedHEtaPP->GetParameter(1),fitTCMDecomposedHEtaPP->GetParameter(2)};

        //Two component model from Bylinkin
        fitTCMInvXSetionEtaPP[ppEst]                            = FitObject("tcm",Form("fitTCMInvXSetionEtaPP5TeV%s",addNamePP[ppEst].Data()),"Eta",graphPPCombEtaStat[ppEst], 0.3, 20. ,paramTCMEtaNewPP,"QNRMEX0+","", kFALSE);
        cout << "fitting pp eta spectra" << endl;
        cout << WriteParameterToFile(fitTCMInvXSetionEtaPP[ppEst])<< endl;
        fitTsallisInvXSectionEtaPP[ppEst]                       = FitObject("l",Form("fitTsallisInvXSectionEtaPP5TeV%s",addNamePP[ppEst].Data()),"Eta",graphPPCombEtaStat[ppEst], 0.3, 20. ,
                                                                            NULL,"QNRMEX0+","", kFALSE);
        cout << "fitting pp eta spectra" << endl;
        cout << WriteParameterToFile(fitTsallisInvXSectionEtaPP[ppEst])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTsallisInvXSectionEtaPP[ppEst])<< endl;

        fileFitsOutput <<  WriteParameterToFile(fitTCMInvXSetionEtaPP[ppEst])<< endl;
        fitPPTCMInvYieldEta[ppEst]                              = FitObject("tcm",Form("fitTCMInvYieldEtaPP5TeV%s",addNamePP[ppEst].Data()),"Eta",NULL, 0.3, 20. ,NULL,"QNRMEX0+","", kFALSE);
        cout << "fitting pp eta spectra" << endl;
        fitPPTCMInvYieldEta[ppEst]->SetParameter(0,graphPPInvYieldCombEtaStat[ppEst]->GetY()[2]);
        fitPPTCMInvYieldEta[ppEst]->SetParameter(1,fitTCMInvXSetionEtaPP[ppEst]->GetParameter(1));
        fitPPTCMInvYieldEta[ppEst]->SetParameter(2,graphPPInvYieldCombEtaStat[ppEst]->GetY()[2]);
        fitPPTCMInvYieldEta[ppEst]->SetParameter(3,fitTCMInvXSetionEtaPP[ppEst]->GetParameter(3));
        fitPPTCMInvYieldEta[ppEst]->SetParameter(4,fitTCMInvXSetionEtaPP[ppEst]->GetParameter(4));
        graphPPInvYieldCombEtaStat[ppEst]->Fit(fitPPTCMInvYieldEta[ppEst],"QNRMEX0+","", 0.3, 20.);
        cout << WriteParameterToFile(fitPPTCMInvYieldEta[ppEst])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitPPTCMInvYieldEta[ppEst])<< endl;
        fitPPTCMInvYieldEta[ppEst]->SetRange(0.3, 20.);

        //     TGraphAsymmErrors* graphRatioPPEtaNLODSS14          = (TGraphAsymmErrors*)graphNLODSS14EtaPP->Clone("graphRatioPPEtaNLODSS14ToFit");
        //     graphRatioPPEtaNLODSS14                             = CalculateGraphErrRatioToFit(graphRatioPPEtaNLODSS14, fitPPTCMInvYieldEta[ppEst]);
        //     TGraph* graphRatioPPEtaNLODSS14Center            mei   = (TGraph*)graphNLODSS14EtaPPCenter->Clone("graphRatioPPEtaNLODSS14CenterToFit");
        //     graphRatioPPEtaNLODSS14Center                       = CalculateGraphRatioToFit(graphRatioPPEtaNLODSS14Center, fitPPTCMInvYieldEta[ppEst]);


        // *************************************************************************************************************
        // Shift graphs in Y direction as well if desired
        // *************************************************************************************************************

        graphRatioEtaCombToFitStatPP[ppEst]                         = (TGraphAsymmErrors*)graphPPCombEtaStat[ppEst]->Clone();
        graphRatioEtaCombToFitStatPP[ppEst]                         = CalculateGraphErrRatioToFit(graphRatioEtaCombToFitStatPP[ppEst], fitTCMInvXSetionEtaPP[ppEst]);
        graphRatioEtaCombToFitSysPP[ppEst]                          = (TGraphAsymmErrors*)graphPPCombEtaFullSys[ppEst]->Clone();
        graphRatioEtaCombToFitSysPP[ppEst]                          = CalculateGraphErrRatioToFit(graphRatioEtaCombToFitSysPP[ppEst], fitTCMInvXSetionEtaPP[ppEst]);

        for (Int_t meth = 0; meth< 11; meth++){
            if (statErrorCollectionEtaPP[ppEst][meth]){
                graphRatioEtaIndCombFitStatPP[ppEst][meth]                = (TGraphAsymmErrors*)statErrorCollectionEtaPP[ppEst][meth]->Clone(Form("RatioEta%sToCombFitStatPP", nameMeasGlobalLabel[meth].Data()));
                graphRatioEtaIndCombFitStatPP[ppEst][meth]                = CalculateGraphErrRatioToFit(graphRatioEtaIndCombFitStatPP[ppEst][meth], fitTCMInvXSetionEtaPP[ppEst]);
            }
            if (systErrorCollectionEtaPP[ppEst][meth]){
                graphRatioEtaIndCombFitSysPP[ppEst][meth]                 = (TGraphAsymmErrors*)systErrorCollectionEtaPP[ppEst][meth]->Clone(Form("RatioEta%sToCombFitSystPP", nameMeasGlobalLabel[meth].Data()));
                graphRatioEtaIndCombFitSysPP[ppEst][meth]                 = CalculateGraphErrRatioToFit(graphRatioEtaIndCombFitSysPP[ppEst][meth], fitTCMInvXSetionEtaPP[ppEst]);
            }
        }

        // **********************************************************************************************************************
        // ******************************************* Plot Ratio of Comb to Fit for pp *****************************************
        // **********************************************************************************************************************
        canvasRatioToCombFitPP->cd();
        histo2DEtaRatioToCombFitPP->Draw("copy");

        ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombToFitStatPP[ppEst]);

        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombToFitSysPP[ppEst], markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioEtaCombToFitSysPP[ppEst]->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombToFitStatPP[ppEst], markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioEtaCombToFitStatPP[ppEst]->Draw("p,same,z");

        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1., 1.,0.1, kGray+2);
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.9, 0.9,0.1, kGray, 7);

        labelRatioToFitEnergyPP->Draw();
        labelRatioToFitALICEPP->Draw();
        labelRatioToFitEtaPP->Draw();

        canvasRatioToCombFitPP->SaveAs(Form("%s/Eta_RatioOfCombToCombFit_PP5TeV%s.%s",outputDirSupportPaper.Data(), addNamePP[ppEst].Data(), suffix.Data()));

        // **********************************************************************************************************************
        // *******************************************Plot Ratio of Individual meas to Fit ******************************************
        // **********************************************************************************************************************

        canvasRatioToCombFitPP->cd();
        histo2DEtaRatioToCombFitPP->Draw("copy");

        for (Int_t meth = 10; meth > -1 ; meth--){
            if (graphRatioEtaIndCombFitSysPP[ppEst][meth]){
                DrawGammaSetMarkerTGraphAsym(graphRatioEtaIndCombFitSysPP[ppEst][meth], markerStyleDetMC[meth] ,markerSizeDetMC[meth]*0.5, colorDet[meth], colorDet[meth], widthLinesBoxes, kTRUE);
                graphRatioEtaIndCombFitSysPP[ppEst][meth]->Draw("E2same");
            }
            if (graphRatioEtaIndCombFitStatPP[ppEst][meth]){
                ProduceGraphAsymmWithoutXErrors(graphRatioEtaIndCombFitStatPP[ppEst][meth]);
                DrawGammaSetMarkerTGraphAsym(graphRatioEtaIndCombFitStatPP[ppEst][meth], markerStyleDetMC[meth] ,markerSizeDetMC[meth]*0.5, colorDet[meth], colorDet[meth]);
                graphRatioEtaIndCombFitStatPP[ppEst][meth]->Draw("p,same,z");
            }
        }
        if (graphRatioEtaIndCombFitStatPP[ppEst][4])graphRatioEtaIndCombFitStatPP[ppEst][4]->Draw("p,same,z");

        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1., 1.,0.5, kGray+2);
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.9, 0.9,0.5, kGray, 7);
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1.2, 1.2,0.5, kGray, 9);
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.8, 0.8,0.5, kGray, 9);

        labelRatioToFitEnergyPP->Draw();
        labelRatioToFitALICEPP->Draw();
        labelRatioToFitEtaPP->Draw();
        histo2DEtaRatioToCombFitPP->Draw("same,axis");

        //****************** first Column *************************************************
        textPCMOnlyRatioEtaPP->Draw();
        textPHOSOnlyRatioEtaPP->Draw();
        textEMCALOnlyRatioEtaPP->Draw();
        textPCMEMCALOnlyRatioEtaPP->Draw();
        textPCMPHOSOnlyRatioEtaPP->Draw();

        //****************** second Column *************************************************
        textStatOnlyRatioEtaPP->Draw();
        textSysOnlyRatioEtaPP->Draw();
        textStatOnlyRatioEtaPP2->Draw();
        textSysOnlyRatioEtaPP2->Draw();

        TMarker* markerPCMEtaOnlyRatioEtaPP             = CreateMarkerFromGraph(graphRatioEtaIndCombFitSysPP[ppEst][0],columnsLegendOnlyEtaRatioPP[1] ,rowsLegendOnlyEtaRatioPP[1],1);
        markerPCMEtaOnlyRatioEtaPP->DrawMarker(columnsLegendOnlyEtaRatioPPAbs[1] ,rowsLegendOnlyEtaRatioPPAbs[1]);
        TMarker* markerPHOSEtaOnlyRatioEtaPP            = CreateMarkerFromGraph(graphRatioEtaIndCombFitSysPP[ppEst][1], columnsLegendOnlyEtaRatioPP[1] ,rowsLegendOnlyEtaRatioPP[2],1);
        markerPHOSEtaOnlyRatioEtaPP->DrawMarker(columnsLegendOnlyEtaRatioPPAbs[1] ,rowsLegendOnlyEtaRatioPPAbs[2]);
        TMarker* markerEMCALEtaOnlyRatioEtaPP           = CreateMarkerFromGraph(graphRatioEtaIndCombFitSysPP[ppEst][2], columnsLegendOnlyEtaRatioPP[1] ,rowsLegendOnlyEtaRatioPP[3],1);
        markerEMCALEtaOnlyRatioEtaPP->DrawMarker(columnsLegendOnlyEtaRatioPPAbs[1] ,rowsLegendOnlyEtaRatioPPAbs[3]);
        TMarker* markerPCMEMCALEtaOnlyRatioEtaPP        = CreateMarkerFromGraph(graphRatioEtaIndCombFitSysPP[ppEst][4], columnsLegendOnlyEtaRatioPP[3] ,rowsLegendOnlyEtaRatioPP[1],1);
        markerPCMEMCALEtaOnlyRatioEtaPP->DrawMarker(columnsLegendOnlyEtaRatioPPAbs[4] ,rowsLegendOnlyEtaRatioPPAbs[1]);
        TMarker* markerPCMPHOSEtaOnlyRatioEtaPP         = CreateMarkerFromGraph(graphRatioEtaIndCombFitSysPP[ppEst][3], columnsLegendOnlyEtaRatioPP[3] ,rowsLegendOnlyEtaRatioPP[2],1);
        markerPCMPHOSEtaOnlyRatioEtaPP->DrawMarker(columnsLegendOnlyEtaRatioPPAbs[4] ,rowsLegendOnlyEtaRatioPPAbs[2]);

        TBox* boxPCMEtaOnlyRatioEtaPP                 = CreateBoxFromGraph(graphRatioEtaIndCombFitSysPP[ppEst][0], columnsLegendOnlyEtaRatioPPAbs[2]-0.1*lengthBoxPPEta , rowsLegendOnlyEtaRatioPPAbs[1]- heightBoxPP,
                                                                            columnsLegendOnlyEtaRatioPPAbs[2]+ 1.7*lengthBoxPPEta, rowsLegendOnlyEtaRatioPPAbs[1]+ heightBoxPP);
        boxPCMEtaOnlyRatioEtaPP->Draw("l");
        TBox* boxPHOSEtaOnlyRatioEtaPP                = CreateBoxFromGraph(graphRatioEtaIndCombFitSysPP[ppEst][1], columnsLegendOnlyEtaRatioPPAbs[2]-0.1*lengthBoxPPEta , rowsLegendOnlyEtaRatioPPAbs[2]- heightBoxPP,
                                                                            columnsLegendOnlyEtaRatioPPAbs[2]+ 1.7*lengthBoxPPEta, rowsLegendOnlyEtaRatioPPAbs[2]+ heightBoxPP);
        boxPHOSEtaOnlyRatioEtaPP->Draw("l");
        TBox* boxEMCALEtaOnlyRatioEtaPP               = CreateBoxFromGraph(graphRatioEtaIndCombFitSysPP[ppEst][2], columnsLegendOnlyEtaRatioPPAbs[2]-0.1*lengthBoxPPEta , rowsLegendOnlyEtaRatioPPAbs[3]- heightBoxPP,
                                                                            columnsLegendOnlyEtaRatioPPAbs[2]+ 1.7*lengthBoxPPEta, rowsLegendOnlyEtaRatioPPAbs[3]+ heightBoxPP);
        boxEMCALEtaOnlyRatioEtaPP->Draw("l");
        TBox* boxPCMEMCALEtaOnlyRatioEtaPP            = CreateBoxFromGraph(graphRatioEtaIndCombFitSysPP[ppEst][4], columnsLegendOnlyEtaRatioPPAbs[5]-0.5*lengthBoxPPEta , rowsLegendOnlyEtaRatioPPAbs[1]- heightBoxPP,
                                                                            columnsLegendOnlyEtaRatioPPAbs[5]+ 17.*lengthBoxPPEta, rowsLegendOnlyEtaRatioPPAbs[1]+ heightBoxPP);
        boxPCMEMCALEtaOnlyRatioEtaPP->Draw("l");
        TBox* boxPCMPHOSEtaOnlyRatioEtaPP             = CreateBoxFromGraph(graphRatioEtaIndCombFitSysPP[ppEst][3], columnsLegendOnlyEtaRatioPPAbs[5]-0.5*lengthBoxPPEta , rowsLegendOnlyEtaRatioPPAbs[2]- heightBoxPP,
                                                                            columnsLegendOnlyEtaRatioPPAbs[5]+ 17.*lengthBoxPPEta, rowsLegendOnlyEtaRatioPPAbs[2]+ heightBoxPP);
        boxPCMPHOSEtaOnlyRatioEtaPP->Draw("l");

        canvasRatioToCombFitPP->SaveAs(Form("%s/Eta_RatioOfIndividualMeasToCombFit_PP5TeV%s.%s",outputDirSupportPaper.Data(), addNamePP[ppEst].Data(), suffix.Data()));
    }

    // *******************************************************************************************************
    // ************************** Combination of different pi0 measurements **********************************
    // *******************************************************************************************************
    // REMARKS:
    //     - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - extendable to other energies
    //     - offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented

    TH1D* statErrorCollectionPi0[17][11];
    TGraphAsymmErrors* statErrorGraphCollectionPi0[17][11];
    TGraphAsymmErrors* sysErrorCollectionPi0[17][11];
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
            if (statErrorGraphCollectionPi0[cent][meth]) RemoveZerosAtBeginningAndEndFromGraph(statErrorGraphCollectionPi0[cent][meth]);
            // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
            // for systematic error from respective method
            sysErrorCollectionPi0[cent][meth]           = NULL;
            if (graphPi0InvYieldSys[cent][meth]) sysErrorCollectionPi0[cent][meth]        = (TGraphAsymmErrors*)graphPi0InvYieldSys[cent][meth]->Clone(Form("sysErr%i_%sPi0", cent,
                                                                                                                                                            nameMeasGlobalLabel[meth].Data()));
            if (sysErrorCollectionPi0[cent][meth]) RemoveZerosAtBeginningAndEndFromGraph(sysErrorCollectionPi0[cent][meth]);
        }
    }

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsPi0[17][11]             = {
                                            { 0,  6,  0,  0,  0,  3,  0,  0,  0,  0,  0 },  // MB R1
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // MB R1 cent binning
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 0-20
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 20-40
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 40-60
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 60-100
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 0-20 CL1
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 20-40 CL1
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 40-60 CL1
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 60-100 CL1
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 0-20 ZNA
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 20-40 ZNA
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 40-60 ZNA
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 60-100 ZNA
                                            { 0,  6, 16,  0,  8,  3,  0,  0,  0,  0,  0 },  // MB R2
                                            { 0,  6,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 0-10 R2
                                            { 0,  6,  1,  1,  1,  3,  0,  0,  0,  0,  0 }  // 60-80 R2
                                           };
    Int_t offSetsPi0Sys[17][11]          = {
                                            { 1,  7,  9,  3,  6,  4,  0,  0,  0,  21, 0 },  // MB R1
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // MB R1 cent binning
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 0-20
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 20-40
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 40-60
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 60-100
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 0-20 CL1
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 20-40 CL1
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 40-60 CL1
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 60-100 CL1
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 0-20 ZNA
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 20-40 ZNA
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 40-60 ZNA
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 60-100 ZNA
                                            { 1,  7, 23,  5, 10,  0,  0,  0,  0,  21, 0 },  // MB R2
                                            { 1,  7,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 0-10
                                            { 1,  7,  8,  2,  5,  0,  0,  0,  0,  21, 0 }   // 60-80
                                           };
    Int_t offSetPi0Shifting[17][11]      = {
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // MB R1
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // MB R1 cent binning
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 0-20
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 20-40
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 40-60
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 60-100
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 0-20 CL1
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 20-40 CL1
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 40-60 CL1
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 60-100 CL1
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 0-20 ZNA
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 20-40 ZNA
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 40-60 ZNA
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 60-100 ZNA
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // MB R2
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 0-10
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 }   // 60-80
                                           };
    Int_t nComBinsPi0Shifting[17][11]    = {
                                            { 29, 10, 17, 17,  29, 17,  0,  0,  0,  0, 0 },  // MB R1
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // MB R1 cent binning
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 0-20
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 20-40
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 40-60
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 60-100
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 0-20 CL1
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 20-40 CL1
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 40-60 CL1
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 60-100 CL1
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 0-20 ZNA
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 20-40 ZNA
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 40-60 ZNA
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 60-100 ZNA
                                            { 0, 0, 0, 0,  0, 0,  0,  0,  0,  0, 0 },  // MB R2
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 0-10
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 }  // 60-80
    };
    Double_t minPtPi0[17]                = { 0.3, 0.4, 0.4, 0.4, 0.4,   0.4, 0.4, 0.4, 0.4, 0.4,
                                             0.4, 0.4, 0.4, 0.4, 0.3,   0.4, 0.4};

    TGraphAsymmErrors* statErrorRelCollectionPi0[17][11];
    TGraphAsymmErrors* sysErrorRelCollectionPi0[17][11];
    TGraph* graphWeightsPi0[17][11];
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
                RemoveZerosAtBeginningAndEndFromGraph(statErrorRelCollectionPi0[cent][meth]);
                statErrorRelCollectionPi0[cent][meth]   = CalculateRelErrUpAsymmGraph( statErrorRelCollectionPi0[cent][meth], Form("relativeStatErrorPi0_%i_%s", cent,nameMeasGlobal[meth].Data()));
            }
            if (sysErrorCollectionPi0[cent][meth]){
                sysErrorRelCollectionPi0[cent][meth]    = (TGraphAsymmErrors*)sysErrorCollectionPi0[cent][meth]->Clone(Form("relativeSysErrorPi0%i_%s", cent, nameMeasGlobal[meth].Data()));
                RemoveZerosAtBeginningAndEndFromGraph(sysErrorRelCollectionPi0[cent][meth]);
                sysErrorRelCollectionPi0[cent][meth]    = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionPi0[cent][meth], Form("relativeSysErrorPi0%i_%s", cent, nameMeasGlobal[meth].Data()));
            }
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************

    TGraphAsymmErrors* graphCombPi0InvYieldStat[17]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldStatWOXErr[17]       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldSys[17]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldTot[17]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldTotUnshi[17]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldStatUnshi[17]        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldSysUnshi[17]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldRelStat[17]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldRelSys[17]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldRelTot[17]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldStat_yShifted[17]    = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldSys_yShifted[17]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombPi0InvYieldTot_yShifted[17]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioPi0CombCombFitTot[17]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioPi0CombCombFitStat[17]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioPi0CombCombFitSys[17]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioPi0CombCombFitStatWOXErr[17]   = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1D* histoRatioPi0DPMJetToFit[17]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1D* histoRatioPi0HIJINGToFit[17]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1D* histoRatioPi0EPOSLHCToFit[17]                         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };

    TGraphAsymmErrors* graphIndPi0InvYieldStatUnshi[17][11];
    TGraphAsymmErrors* graphIndPi0InvYieldSysUnshi[17][11];
    TGraphAsymmErrors* graphIndPi0InvYieldStat[17][11];
    TGraphAsymmErrors* graphIndPi0InvYieldSys[17][11];
    TGraphAsymmErrors* graphIndPi0InvYieldStat_yShifted[17][11];
    TGraphAsymmErrors* graphIndPi0InvYieldSys_yShifted[17][11];
    TGraphAsymmErrors* graphRatioPi0IndCombFitStat[17][11];
    TGraphAsymmErrors* graphRatioPi0IndCombFitSys[17][11];
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
    TF1* fitTCMDecomposedLPi0[17]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitTCMDecomposedHPi0[17]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitTCMInvYieldPi0[17]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitInvYieldPi0[17]                             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitInvYieldPi0Graph[17]                        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitPowInvYieldPi0[17]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };

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
    TLatex *labelWeightsPi0RpPb     = new TLatex(0.95,0.15,"#it{R}_{pPb}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelWeightsPi0RpPb, 0.85*textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelWeightsPi0RCP     = new TLatex(0.95,0.15,"#it{R}_{CP}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelWeightsPi0RCP, 0.85*textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelWeightsPi0RMB     = new TLatex(0.95,0.15,"#it{R}_{MB}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelWeightsPi0RMB, 0.85*textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);


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
    TLatex *labelRelSysErrPi0RpPb    = new TLatex(0.15,0.85,"#it{R}_{pPb}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrPi0RpPb, textSizeLabelsPixel, 4, 1, 43);
    TLatex *labelRelSysErrPi0RCP    = new TLatex(0.15,0.85,"#it{R}_{CP}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrPi0RCP, textSizeLabelsPixel, 4, 1, 43);
    TLatex *labelRelSysErrPi0RMB    = new TLatex(0.15,0.85,"#it{R}_{MB}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrPi0RMB, textSizeLabelsPixel, 4, 1, 43);


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
    TLatex *labelRelStatErrPi0RpPb  = new TLatex(0.95,0.85,"#it{R}_{pPb}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelStatErrPi0RpPb, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelStatErrPi0RCP  = new TLatex(0.95,0.85,"#it{R}_{CP}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelStatErrPi0RCP, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelStatErrPi0RMB  = new TLatex(0.95,0.85,"#it{R}_{MB}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelStatErrPi0RMB, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

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
    TLatex *labelRelTotErrPi0RpPb       = new TLatex(0.95,0.85,"#it{R}_{pPb}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrPi0RpPb, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelTotErrPi0RCP       = new TLatex(0.95,0.85,"#it{R}_{CP}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrPi0RCP, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelTotErrPi0RMB       = new TLatex(0.95,0.85,"#it{R}_{MB}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrPi0RMB, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);


    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCentComb[cent]) continue;
        cout << "start combining pi0 for " << centArray[cent].Data() << " " << runArray[cent].Data() << endl;
        // Declaration & calculation of combined spectrum
        TString fileNamePi0OutputWeighting      = Form("%s/Pi0_WeightingMethod_%s%s%s.dat",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data());
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
//         graphCombPi0InvYieldTot[cent]->Print();

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
//             cout << cent << "\t" << availablePi0Meas[i] << endl;
//             graphWeightsPi0[cent][availablePi0Meas[i]]->Print();
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

            labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelWeightsEnergy->Draw();
            labelWeightsPi0->Draw();

            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.2, 0.2,0.1, kGray, 3);

        canvasWeights->SaveAs(Form("%s/Pi0_Weights_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

        //  *********************************************************************************************************************
        //  ************************************ Visualize relative errors ******************************************************
        //  *********************************************************************************************************************

        canvasRelSysErr->cd();
        histo2DRelSysErr->GetYaxis()->SetRangeUser(0,39.5);
        histo2DRelSysErr->Draw("copy");

            cout << "sys error pi0" << endl;
            TLegend* legendRelSysErr        = GetAndSetLegend2(0.62, 0.92-(0.04*nMeasSetPi0), 0.95, 0.92, textSizeLabelsPixel);
            for (Int_t i = 0; i < nMeasSetPi0; i++){
                DrawGammaSetMarkerTGraph(sysErrorRelCollectionPi0[cent][availablePi0Meas[i]], markerStyleDet[availablePi0Meas[i]], markerSizeDet[availablePi0Meas[i]]*0.5, colorDet[availablePi0Meas[i]],
                                         colorDet[availablePi0Meas[i]]);
                sysErrorRelCollectionPi0[cent][availablePi0Meas[i]]->Draw("p,same,z");
                legendRelSysErr->AddEntry(sysErrorRelCollectionPi0[cent][availablePi0Meas[i]],nameMeasGlobalLabel[availablePi0Meas[i]],"p");
//                 cout << cent << "\t" << availablePi0Meas[i] << endl;
//                 sysErrorRelCollectionPi0[cent][availablePi0Meas[i]]->Print();

            }
            legendRelSysErr->Draw();


            labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelSysErrEnergy->Draw();
            labelRelSysErrPi0->Draw();

        canvasRelSysErr->SaveAs(Form("%s/Pi0_RelSysErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
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
//                 cout << Form("%s stat: %s", centArray[cent].Data(), nameMeasGlobalLabel[availablePi0Meas[i]].Data() ) << endl;
//                 statErrorRelCollectionPi0[cent][availablePi0Meas[i]]->Print();
            }
            legendRelStatErr->Draw();

            labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelStatErrEnergy->Draw();
            labelRelStatErrPi0->Draw();

        canvasRelStatErr->SaveAs(Form("%s/Pi0_RelStatErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

        //  *********************************************************************************************************************
        //  ************************ Visualize relative total errors of different combination methods Pi0 ***********************
        //  *********************************************************************************************************************

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


            labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelTotErrEnergy->Draw();
            labelRelTotErrPi0->Draw();

        canvasRelTotErr->SaveAs(Form("%s/Pi0_TotErr_Comp_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
        histo2DRelTotErrPi0->GetYaxis()->SetTitle("Err (%)");
        histo2DRelTotErrPi0->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphCombPi0InvYieldRelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelStat[cent], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
            graphCombPi0InvYieldRelStat[cent]->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelSys[cent], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
            graphCombPi0InvYieldRelSys[cent]->SetLineStyle(7);
            graphCombPi0InvYieldRelSys[cent]->Draw("l,x0,same,e1");

            if (centArray[cent].BeginsWith("0-100%")){
                DrawGammaSetMarkerTGraphAsym(graphCombPi0OlderMBRelErrTot, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
                graphCombPi0OlderMBRelErrTot->Draw("p,same,z");
                DrawGammaSetMarkerTGraphAsym(graphCombPi0OlderMBRelErrStat, markerStyleComb+4, markerSizeComb, kGray+2 , kGray+2);
                graphCombPi0OlderMBRelErrStat->Draw("l,x0,same,e1");
                DrawGammaSetMarkerTGraphAsym(graphCombPi0OlderMBRelErrSys, markerStyleComb+4, markerSizeComb, kGray+1, kGray+1);
                graphCombPi0OlderMBRelErrSys->SetLineStyle(7);
                graphCombPi0OlderMBRelErrSys->Draw("l,x0,same,e1");

                TLegend* legendRelTotErr4       = GetAndSetLegend2(0.35, 0.92-(0.035*3), 0.65, 0.92, 32);
                legendRelTotErr4->AddEntry(graphCombPi0OlderMBRelErrTot,"pub tot","p");
                legendRelTotErr4->AddEntry(graphCombPi0OlderMBRelErrStat,"pub stat","l");
                legendRelTotErr4->AddEntry(graphCombPi0OlderMBRelErrSys,"pub sys","l");
                legendRelTotErr4->Draw();

            }

            TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.35, 0.92, 32);
            legendRelTotErr3->AddEntry(graphCombPi0InvYieldRelTot[cent],"tot","p");
            legendRelTotErr3->AddEntry(graphCombPi0InvYieldRelStat[cent],"stat","l");
            legendRelTotErr3->AddEntry(graphCombPi0InvYieldRelSys[cent],"sys","l");
            legendRelTotErr3->Draw();

            labelRelTotErrEnergy->Draw();
            labelRelTotErrPi0->Draw();

        canvasRelTotErr->SaveAs(Form("%s/Pi0_Reldecomp_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

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
        graphCombPi0InvYieldStat[cent]->Fit(fitTCMDecomposedHPi0[cent],"QNRMEX0+","",3, xPtLimitsPi0[cent][maxNBinsPi0[cent]]);

        cout << WriteParameterToFile(fitTCMDecomposedLPi0[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedLPi0[cent])<< endl;
        cout << WriteParameterToFile(fitTCMDecomposedHPi0[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedHPi0[cent])<< endl;

        fitTCMInvYieldPi0[cent]              = FitObject("tcm",Form("fitTCMInvYieldPi0_%s", centArrayOutput[cent].Data()),"Pi0",graphCombPi0InvYieldStat[cent],
                                                         minPtPi0[cent], xPtLimitsPi0[cent][maxNBinsPi0[cent]] , paramTCMPi0New,
                                                         "QNRMEX0+", "", kFALSE);
        fitTCMDecomposedLPi0[cent]->SetParameter(0, fitTCMInvYieldPi0[cent]->GetParameter(0));
        fitTCMDecomposedLPi0[cent]->SetParameter(1, fitTCMInvYieldPi0[cent]->GetParameter(1));
        fitTCMDecomposedHPi0[cent]->SetParameter(0, fitTCMInvYieldPi0[cent]->GetParameter(2));
        fitTCMDecomposedHPi0[cent]->SetParameter(1, fitTCMInvYieldPi0[cent]->GetParameter(3));
        fitTCMDecomposedHPi0[cent]->SetParameter(2, fitTCMInvYieldPi0[cent]->GetParameter(4));

        // Tsallis fit
        fitInvYieldPi0[cent]                 = FitObject("l","fitInvYieldPi0","Pi0",histoPi0InvYieldStat[cent][2], minPtPi0[cent], xPtLimitsPi0[cent][maxNBinsPi0[cent]], paramGraphPi0,"QNRMEX0+");
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
//             cout << "comb" << endl;
//             graphCombPi0InvYieldStat[cent]->Print();
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

                TLatex *labelRatioToFitBinShift   = new TLatex(0.94, 0.91, Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
                SetStyleTLatex( labelRatioToFitBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitBinShift->Draw();
                TLatex *labelRatioToFitALICEBinShift    = new TLatex(0.94, 0.86, textALICE.Data());
                SetStyleTLatex( labelRatioToFitALICEBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitALICEBinShift->Draw();
                TLatex *labelRatioToFitPi0BinShift      = new TLatex(0.94, 0.807, "#pi^{0} #rightarrow #gamma#gamma");
                SetStyleTLatex( labelRatioToFitPi0BinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitPi0BinShift->Draw();

            canvasShift->Update();
            canvasShift->SaveAs(Form("%s/BinShiftCorrection_Pi0_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

            // *************************************************************************************************************
            // Plot control graphs
            // *************************************************************************************************************

            TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
            DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.09);
            canvasDummy2->SetLogy();
            canvasDummy2->SetLogx();
            TH2F * histo2DDummy2;
            histo2DDummy2               = new TH2F("histo2DDummy2","histo2DDummy2",1000,minPtPi0Plotting,maxPtPi0Plotting,1000,1e-9,10e1);
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
            canvasDummy2->Print(Form("%s/ComparisonShiftedPi0_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete canvasDummy2;
        }

        // *************************************************************************************************************
        // redo fitting after binshifts
        // *************************************************************************************************************
        // Tsallis function
        cout << WriteParameterToFile(fitInvYieldPi0[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitInvYieldPi0[cent])<< endl;
        //Two component model from Bylinkin
        fitTCMInvYieldPi0[cent]        = FitObject("tcm","fitTCMInvYieldPi0pPb5TeVCent","Pi0",graphCombPi0InvYieldTot[cent], minPtPi0[cent], xPtLimitsPi0[cent][maxNBinsPi0[cent]],
                                                   paramTCMPi0New,"QNRMEX0+","", kFALSE);
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
        if (histoDPMJetPi0[cent]){
            histoRatioPi0DPMJetToFit[cent]                     = (TH1D*) histoDPMJetPi0[cent]->Clone(Form("histoRatioPi0DPMJetToFit%s%s",centArray[cent].Data(), runArray[cent].Data() ));
            histoRatioPi0DPMJetToFit[cent]                     = CalculateHistoRatioToFit (histoRatioPi0DPMJetToFit[cent], fitTCMInvYieldPi0[cent]);
            histoRatioPi0DPMJetToFit[cent]->GetXaxis()->SetRangeUser(minPtPi0[cent], xPtLimitsPi0[cent][maxNBinsPi0[cent]]);
        }
        if (histoHIJINGPi0[cent]){
            histoRatioPi0HIJINGToFit[cent]                     = (TH1D*) histoHIJINGPi0[cent]->Clone(Form("histoRatioPi0HIJINGToFit%s%s",centArray[cent].Data(), runArray[cent].Data() ));
            histoRatioPi0HIJINGToFit[cent]                     = CalculateHistoRatioToFit (histoRatioPi0HIJINGToFit[cent], fitTCMInvYieldPi0[cent]);
            histoRatioPi0HIJINGToFit[cent]->GetXaxis()->SetRangeUser(minPtPi0[cent], xPtLimitsPi0[cent][maxNBinsPi0[cent]]);
        }
        if (histoEPOSLHCPi0[cent]){
            histoRatioPi0EPOSLHCToFit[cent]                    = (TH1D*) histoEPOSLHCPi0[cent]->Clone(Form("histoRatioPi0EPOSLHCToFit%s%s",centArray[cent].Data(), runArray[cent].Data() ));
            histoRatioPi0EPOSLHCToFit[cent]                    = CalculateHistoRatioToFit (histoRatioPi0EPOSLHCToFit[cent], fitTCMInvYieldPi0[cent]);
            histoRatioPi0EPOSLHCToFit[cent]->GetXaxis()->SetRangeUser(minPtPi0[cent], xPtLimitsPi0[cent][maxNBinsPi0[cent]]);
        }

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
        graphCombPi0InvYieldStatWOXErr[cent]                    = (TGraphAsymmErrors*)graphCombPi0InvYieldStat[cent]->Clone("graphCombPi0InvYieldStatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphCombPi0InvYieldStatWOXErr[cent]);
        graphRatioPi0CombCombFitStatWOXErr[cent]                = (TGraphAsymmErrors*)graphRatioPi0CombCombFitStat[cent]->Clone("graphRatioPi0CombCombFitStatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatWOXErr[cent]);

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

            TLatex *labelRatioToFitEnergy   = new TLatex(0.95, 0.92, Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
            labelRatioToFitEnergy->Draw();
            TLatex *labelRatioToFitALICE    = new TLatex(0.95, 0.86, textALICE.Data());
            SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
            labelRatioToFitALICE->Draw();
            TLatex *labelRatioToFitPi0      = new TLatex(0.12, 0.92, "#pi^{0} #rightarrow #gamma#gamma");
            SetStyleTLatex( labelRatioToFitPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
            labelRatioToFitPi0->Draw();

        canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToCombFit_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

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
            Double_t columnsLegendOnlyPi0RatioAbs[6]    = {0.215, 0.8, 1.17, 2, 11.2, 16.8};
            Double_t columnsLegendOnlyPi0RatioAbs2[6]   = {0.215, 0.8, 1.68, 2, 11.2, 24.4};
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
                TBox* boxPCMPi0OnlyRatioPi0                 = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[cent][0], columnsLegendOnlyPi0RatioAbs[2] , rowsLegendOnlyPi0RatioAbs[1]- heightBox,
                                                                                columnsLegendOnlyPi0RatioAbs2[2], rowsLegendOnlyPi0RatioAbs[1]+ heightBox);
                boxPCMPi0OnlyRatioPi0->Draw("l");
            }
            if (graphRatioPi0IndCombFitSys[cent][1]){
                TBox* boxPHOSPi0OnlyRatioPi0                = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[cent][1], columnsLegendOnlyPi0RatioAbs[2] , rowsLegendOnlyPi0RatioAbs[2]- heightBox,
                                                                             columnsLegendOnlyPi0RatioAbs2[2], rowsLegendOnlyPi0RatioAbs[2]+ heightBox);
                boxPHOSPi0OnlyRatioPi0->Draw("l");
            }
            if (graphRatioPi0IndCombFitSys[cent][2]){
                TBox* boxEMCALPi0OnlyRatioPi0               = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[cent][2], columnsLegendOnlyPi0RatioAbs[2] , rowsLegendOnlyPi0RatioAbs[3]- heightBox,
                                                                                columnsLegendOnlyPi0RatioAbs2[2], rowsLegendOnlyPi0RatioAbs[3]+ heightBox);
                boxEMCALPi0OnlyRatioPi0->Draw("l");
            }
            if (graphRatioPi0IndCombFitSys[cent][4]){
                TBox* boxPCMEMCALPi0OnlyRatioPi0            = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[cent][4], columnsLegendOnlyPi0RatioAbs[5] , rowsLegendOnlyPi0RatioAbs[1]- heightBox,
                                                                                columnsLegendOnlyPi0RatioAbs2[5], rowsLegendOnlyPi0RatioAbs[1]+ heightBox);
                boxPCMEMCALPi0OnlyRatioPi0->Draw("l");
            }
            if (graphRatioPi0IndCombFitSys[cent][3]){
                TBox* boxPCMPHOSPi0OnlyRatioPi0             = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[cent][3], columnsLegendOnlyPi0RatioAbs[5], rowsLegendOnlyPi0RatioAbs[2]- heightBox,
                                                                                columnsLegendOnlyPi0RatioAbs2[5], rowsLegendOnlyPi0RatioAbs[2]+ heightBox);
                boxPCMPHOSPi0OnlyRatioPi0->Draw("l");
            }
            if (graphRatioPi0IndCombFitSys[cent][5]){
                TBox* boxDalitzPi0OnlyRatioPi0             = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[cent][5], columnsLegendOnlyPi0RatioAbs[5], rowsLegendOnlyPi0RatioAbs[3]- heightBox,
                                                                                 columnsLegendOnlyPi0RatioAbs2[5], rowsLegendOnlyPi0RatioAbs[3]+ heightBox);
                boxDalitzPi0OnlyRatioPi0->Draw("l");
            }
            canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFit_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(),
                                              runArray[cent].Data() ,suffix.Data()));
    }

    // *******************************************************************************************************
    // ************************** Combination of different eta measurements **********************************
    // *******************************************************************************************************
    // REMARKS:
    //     - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - extendable to other energies
    //     - offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented

    TH1D* statErrorCollectionEta[17][11];
    TGraphAsymmErrors* statErrorGraphCollectionEta[17][11];
    TGraphAsymmErrors* sysErrorCollectionEta[17][11];
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
            if (statErrorGraphCollectionEta[cent][meth]) RemoveZerosAtBeginningAndEndFromGraph(statErrorGraphCollectionEta[cent][meth]);
            // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
            // for systematic error from respective method
            sysErrorCollectionEta[cent][meth]           = NULL;
            if (graphEtaInvYieldSys[cent][meth]) sysErrorCollectionEta[cent][meth]        = (TGraphAsymmErrors*)graphEtaInvYieldSys[cent][meth]->Clone(Form("sysErr%i_%sEta", cent,
                                                                                                                                                            nameMeasGlobalLabel[meth].Data()));
            if (sysErrorCollectionEta[cent][meth]) RemoveZerosAtBeginningAndEndFromGraph(sysErrorCollectionEta[cent][meth]);
        }
    }

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Eta
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsEta[17][11]             = {
                                            { 0,  9,  0,  1,  0,  0,  0,  0,  0,  0,  0 },  // MB R1
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // MB R1 cent binning
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-20
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 20-40
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 40-60
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 60-100
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-20 CL1
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 20-40 CL1
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 40-60 CL1
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 60-100 CL1
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-20 ZNA
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 20-40 ZNA
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 40-60 ZNA
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 60-100 ZNA
                                            { 0,  7,  10,  5,  7,  0,  0,  0,  0,  0,  0 },  // MB R2
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-10
                                            { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 }  // 60-80
    };
    Int_t offSetsEtaSys[17][11]          = {
                                            { 3,  12,  7,  5,  5,  0,  0,  0,  0,  0, 0 },  // MB R1
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // MB R1 cent binning
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 0-20
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 20-40
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 40-60
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 60-100
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 0-20 CL1
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 20-40 CL1
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 40-60 CL1
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 60-100 CL1
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 0-20 ZNA
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 20-40 ZNA
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 40-60 ZNA
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 60-100 ZNA
                                            { 2,  11,  13,  8,  10,  0,  0,  0,  0,  0, 0 },  // MB R2
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 },  // 0-10
                                            { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 }  // 60-80
                                            };
    Int_t offSetEtaShifting[17][11]      = {
                                            { 0,  0,  4,  2,  2,  0,  0,  0,  0,  0,  0 },  // MB R1
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 },  // MB R1 cent binning
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 },  // 0-20
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 },  // 20-40
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 },  // 40-60
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 },  // 60-100
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 },  // 0-20 CL1
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 },  // 20-40 CL1
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 },  // 40-60 CL1
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 },  // 60-100 CL1
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 },  // 0-20 ZNA
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 },  // 20-40 ZNA
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 },  // 40-60 ZNA
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 },  // 60-100 ZNA
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },   // MB R2
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 0-10
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 }  // 60-80
                                            };
    Int_t nComBinsEtaShifting[17][11]    = {
                                            { 13, 0,  17, 11, 15,  0,  0,  0,  0,  0, 0 },  // MB R1
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // MB R1 cent binning
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // 0-20
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // 20-40
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // 40-60
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // 60-100
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // 0-20 CL1
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // 20-40 CL1
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // 40-60 CL1
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // 60-100 CL1
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // 0-20 ZNA
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // 20-40 ZNA
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // 40-60 ZNA
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // 60-100 ZNA
                                            { 0, 0, 0, 0,  0, 0,  0,  0,  0,  0, 0 },  // MB R2
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 0-10
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },  // 60-80
                                           };
    Double_t minPtEta[17]                = { 0.7, 0.7, 0.7, 0.7, 0.7,   0.7, 0.7, 0.7, 0.7, 0.7,
                                             0.7, 0.7, 0.7, 0.7, 0.4,   0.4, 0.4};

    TGraphAsymmErrors* statErrorRelCollectionEta[17][11];
    TGraphAsymmErrors* sysErrorRelCollectionEta[17][11];
    TGraph* graphWeightsEta[17][11];
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
                RemoveZerosAtBeginningAndEndFromGraph(statErrorRelCollectionEta[cent][meth]);
                statErrorRelCollectionEta[cent][meth]   = CalculateRelErrUpAsymmGraph( statErrorRelCollectionEta[cent][meth], Form("relativeStatErrorEta_%i_%s", cent,nameMeasGlobal[meth].Data()));
            }
            if (sysErrorCollectionEta[cent][meth]){
                sysErrorRelCollectionEta[cent][meth]    = (TGraphAsymmErrors*)sysErrorCollectionEta[cent][meth]->Clone(Form("relativeSysErrorEta%i_%s", cent, nameMeasGlobal[meth].Data()));
                RemoveZerosAtBeginningAndEndFromGraph(sysErrorRelCollectionEta[cent][meth]);
                sysErrorRelCollectionEta[cent][meth]    = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionEta[cent][meth], Form("relativeSysErrorEta%i_%s", cent, nameMeasGlobal[meth].Data()));
            }
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************

    TGraphAsymmErrors* graphCombEtaInvYieldStat[17]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldStatWOXErr[17]       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldSys[17]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldTot[17]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldTotUnshi[17]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldStatUnshi[17]        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldSysUnshi[17]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldRelStat[17]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldRelSys[17]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldRelTot[17]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldStat_yShifted[17]    = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldSys_yShifted[17]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaInvYieldTot_yShifted[17]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioEtaCombCombFitTot[17]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioEtaCombCombFitStat[17]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioEtaCombCombFitStatWOXErr[17]   = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphRatioEtaCombCombFitSys[17]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1D* histoRatioEtaDPMJetToFit[17]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1D* histoRatioEtaHIJINGToFit[17]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TH1D* histoRatioEtaEPOSLHCToFit[17]                         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };

    TGraphAsymmErrors* graphIndEtaInvYieldStatUnshi[17][11];
    TGraphAsymmErrors* graphIndEtaInvYieldSysUnshi[17][11];
    TGraphAsymmErrors* graphIndEtaInvYieldStat[17][11];
    TGraphAsymmErrors* graphIndEtaInvYieldSys[17][11];
    TGraphAsymmErrors* graphIndEtaInvYieldStat_yShifted[17][11];
    TGraphAsymmErrors* graphIndEtaInvYieldSys_yShifted[17][11];
    TGraphAsymmErrors* graphRatioEtaIndCombFitStat[17][11];
    TGraphAsymmErrors* graphRatioEtaIndCombFitSys[17][11];
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
    TF1* fitTCMDecomposedLEta[17]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitTCMDecomposedHEta[17]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitTCMInvYieldEta[17]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitInvYieldEta[17]                             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitInvYieldEtaGraph[17]                        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TF1* fitPowInvYieldEta[17]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };

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
    TLatex *labelWeightsEtaRpPb         = new TLatex(0.95,0.15,"#it{R}_{pA}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelWeightsEtaRpPb, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelWeightsEtaRCP         = new TLatex(0.95,0.15,"#it{R}_{CP}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelWeightsEtaRCP, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelWeightsEtaRMB         = new TLatex(0.95,0.15,"#it{R}_{MB}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelWeightsEtaRMB, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    TLatex *labelRelSysErrEta       = new TLatex(0.15,0.85,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrEta, textSizeLabelsPixel, 4, 1, 43);
    TLatex *labelRelSysErrEtaRpPb       = new TLatex(0.15,0.85,"#it{R}_{pA}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrEtaRpPb, textSizeLabelsPixel, 4, 1, 43);
    TLatex *labelRelSysErrEtaRCP       = new TLatex(0.15,0.85,"#it{R}_{CP}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrEtaRCP, textSizeLabelsPixel, 4, 1, 43);
    TLatex *labelRelSysErrEtaRMB       = new TLatex(0.15,0.85,"#it{R}_{MB}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrEtaRMB, textSizeLabelsPixel, 4, 1, 43);

    TLatex *labelRelStatErrEta      = new TLatex(0.95,0.85,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelStatErrEta, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelStatErrEtaRpPb      = new TLatex(0.95,0.85,"#it{R}_{pA}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelStatErrEtaRpPb, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelStatErrEtaRCP      = new TLatex(0.95,0.85,"#it{R}_{CP}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelStatErrEtaRCP, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelStatErrEtaRMB      = new TLatex(0.95,0.85,"#it{R}_{MB}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelStatErrEtaRMB, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    TH2F * histo2DRelTotErrEta;
    histo2DRelTotErrEta                 = new TH2F("histo2DRelTotErrEta","histo2DRelTotErrEta",11000,minPtEtaPlotting, maxPtEtaPlotting,1000,0,50.0);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrEta, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,39.5);
    histo2DRelTotErrEta->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrEta->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelRelTotErrEta       = new TLatex(0.95,0.85,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrEta, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelTotErrEtaRpPb   = new TLatex(0.95,0.85,"#it{R}_{pA}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrEtaRpPb, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelTotErrEtaRCP   = new TLatex(0.95,0.85,"#it{R}_{CP}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrEtaRCP, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelTotErrEtaRMB   = new TLatex(0.95,0.85,"#it{R}_{MB}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrEtaRMB, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCentComb[cent]) continue;
        // Declaration & calculation of combined spectrum
        cout << "start combining eta for " << centArray[cent].Data() << " " << runArray[cent].Data() << endl;
        TString fileNameEtaOutputWeighting      = Form("%s/Eta_WeightingMethod_%s%s%s.dat",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data());
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
//         graphCombEtaInvYieldTot[cent]->Print();

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

            labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelWeightsEnergy->Draw();
            labelWeightsEta->Draw();

            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.2, 0.2,0.1, kGray, 3);

        canvasWeights->SaveAs(Form("%s/Eta_Weights_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

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

            labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelSysErrEnergy->Draw();
            labelRelSysErrEta->Draw();

        canvasRelSysErr->SaveAs(Form("%s/Eta_RelSysErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
        delete legendRelSysErrEta;
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

            labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelStatErrEnergy->Draw();
            labelRelStatErrEta->Draw();

        canvasRelStatErr->SaveAs(Form("%s/Eta_RelStatErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

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


            labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelTotErrEnergy->Draw();
            labelRelTotErrEta->Draw();

        canvasRelTotErr->SaveAs(Form("%s/Eta_TotErr_Comp_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
        histo2DRelTotErrEta->GetYaxis()->SetTitle("Err (%)");
        histo2DRelTotErrEta->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphCombEtaInvYieldRelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelStat[cent], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
            graphCombEtaInvYieldRelStat[cent]->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelSys[cent], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
            graphCombEtaInvYieldRelSys[cent]->SetLineStyle(7);
            graphCombEtaInvYieldRelSys[cent]->Draw("l,x0,same,e1");

            if (centArray[cent].BeginsWith("0-100%")){
                DrawGammaSetMarkerTGraphAsym(graphCombEtaOlderMBRelErrTot, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
                graphCombEtaOlderMBRelErrTot->Draw("p,same,z");
                DrawGammaSetMarkerTGraphAsym(graphCombEtaOlderMBRelErrStat, markerStyleComb+4, markerSizeComb, kGray+2 , kGray+2);
                graphCombEtaOlderMBRelErrStat->Draw("l,x0,same,e1");
                DrawGammaSetMarkerTGraphAsym(graphCombEtaOlderMBRelErrSys, markerStyleComb+4, markerSizeComb, kGray+1, kGray+1);
                graphCombEtaOlderMBRelErrSys->SetLineStyle(7);
                graphCombEtaOlderMBRelErrSys->Draw("l,x0,same,e1");

                TLegend* legendRelTotErr4       = GetAndSetLegend2(0.35, 0.92-(0.035*3), 0.65, 0.92, 32);
                legendRelTotErr4->AddEntry(graphCombEtaOlderMBRelErrTot,"pub tot","p");
                legendRelTotErr4->AddEntry(graphCombEtaOlderMBRelErrStat,"pub stat","l");
                legendRelTotErr4->AddEntry(graphCombEtaOlderMBRelErrSys,"pub sys","l");
                legendRelTotErr4->Draw();

            }


            TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
            legendRelTotErr3->AddEntry(graphCombEtaInvYieldRelTot[cent],"tot","p");
            legendRelTotErr3->AddEntry(graphCombEtaInvYieldRelStat[cent],"stat","l");
            legendRelTotErr3->AddEntry(graphCombEtaInvYieldRelSys[cent],"sys","l");
            legendRelTotErr3->Draw();

            labelRelTotErrEnergy->Draw();
            labelRelTotErrEta->Draw();

        canvasRelTotErr->SaveAs(Form("%s/Eta_Reldecomp_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

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
        graphCombEtaInvYieldStat[cent]->Fit(fitTCMDecomposedLEta[cent],"QNRMEX0+","",minPtEta[cent],1.5);
        graphCombEtaInvYieldStat[cent]->Fit(fitTCMDecomposedHEta[cent],"QNRMEX0+","",3, xPtLimitsEta[cent][maxNBinsEta[cent]]);

        cout << WriteParameterToFile(fitTCMDecomposedLEta[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedLEta[cent])<< endl;
        cout << WriteParameterToFile(fitTCMDecomposedHEta[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedHEta[cent])<< endl;

        if (centArray[cent].Contains("40-60%") || centArray[cent].Contains("60-100%")){
            paramTCMEtaNew[0]                = 0.504415230;
            paramTCMEtaNew[1]                = 0.269;
            paramTCMEtaNew[2]                = 0.06;
            paramTCMEtaNew[3]                = 1.1379990140;
            paramTCMEtaNew[4]                = 3.37;
        }

        fitTCMInvYieldEta[cent]              = FitObject("tcm",Form("fitTCMInvYieldEta_%s", centArrayOutput[cent].Data()),"Eta",graphCombEtaInvYieldStat[cent],
                                                         minPtEta[cent], xPtLimitsEta[cent][maxNBinsEta[cent]] , paramTCMEtaNew,
                                                         "QNRMEX0+", "", kFALSE);
        if (centArray[cent].Contains("0-20%") ){
            fitTCMInvYieldEta[cent]->SetParameter(0, 2);
            fitTCMInvYieldEta[cent]->SetParLimits(0, 0, 10);
            fitTCMInvYieldEta[cent]->SetParameter(1, 0.12);
            fitTCMInvYieldEta[cent]->SetParLimits(1, 0.09, 0.3);
            graphCombEtaInvYieldStat[cent]->Fit(fitTCMInvYieldEta[cent],"QNRMEX0+","",minPtEta[cent]+0.2, xPtLimitsEta[cent][maxNBinsEta[cent]] );
        } else if (centArray[cent].Contains("20-40%") ){
            fitTCMInvYieldEta[cent]->SetParameter(0, 2);
            fitTCMInvYieldEta[cent]->SetParLimits(0, 0, 5);
            fitTCMInvYieldEta[cent]->SetParameter(1, 0.12);
            fitTCMInvYieldEta[cent]->SetParLimits(1, 0.1, 0.3);
            graphCombEtaInvYieldStat[cent]->Fit(fitTCMInvYieldEta[cent],"QNRMEX0+","",minPtEta[cent]+0.2, xPtLimitsEta[cent][maxNBinsEta[cent]] );
        } else if (centArray[cent].Contains("40-60%") ){
            fitTCMInvYieldEta[cent]->SetParameter(0, 1);
            fitTCMInvYieldEta[cent]->SetParLimits(0, 0, 2);
            fitTCMInvYieldEta[cent]->SetParameter(1, 0.12);
            fitTCMInvYieldEta[cent]->SetParLimits(1, 0.1, 0.3);
            graphCombEtaInvYieldStat[cent]->Fit(fitTCMInvYieldEta[cent],"QNRMEX0+","",minPtEta[cent]+0.2, xPtLimitsEta[cent][maxNBinsEta[cent]] );
        } else if (centArray[cent].Contains("60-100%") ){
            fitTCMInvYieldEta[cent]->SetParameter(0, 0.5);
            fitTCMInvYieldEta[cent]->SetParLimits(0, 0, 1);
            fitTCMInvYieldEta[cent]->SetParameter(1, 0.12);
            fitTCMInvYieldEta[cent]->SetParLimits(1, 0.1, 0.3);
            graphCombEtaInvYieldStat[cent]->Fit(fitTCMInvYieldEta[cent],"QNRMEX0+","",minPtEta[cent]+0.2, xPtLimitsEta[cent][maxNBinsEta[cent]] );
        }


        fitTCMDecomposedLEta[cent]->SetParameter(0, fitTCMInvYieldEta[cent]->GetParameter(0));
        fitTCMDecomposedLEta[cent]->SetParameter(1, fitTCMInvYieldEta[cent]->GetParameter(1));
        fitTCMDecomposedHEta[cent]->SetParameter(0, fitTCMInvYieldEta[cent]->GetParameter(2));
        fitTCMDecomposedHEta[cent]->SetParameter(1, fitTCMInvYieldEta[cent]->GetParameter(3));
        fitTCMDecomposedHEta[cent]->SetParameter(2, fitTCMInvYieldEta[cent]->GetParameter(4));

        // Tsallis fit
        fitInvYieldEta[cent]                 = FitObject("l","fitInvYieldEta","Eta",histoEtaInvYieldStat[cent][2], 0.3, xPtLimitsEta[cent][maxNBinsEta[cent]], paramGraphEta,"QNRMEX0+");
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
//             cout << "comb" << endl;
//             graphCombEtaInvYieldStat[cent]->Print();
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

                TLatex *labelRatioToFitBinShift   = new TLatex(0.94, 0.91, Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
                SetStyleTLatex( labelRatioToFitBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitBinShift->Draw();
                TLatex *labelRatioToFitALICEBinShift    = new TLatex(0.94, 0.86, textALICE.Data());
                SetStyleTLatex( labelRatioToFitALICEBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitALICEBinShift->Draw();
                TLatex *labelRatioToFitEtaBinShift      = new TLatex(0.94, 0.807, "#eta #rightarrow #gamma#gamma");
                SetStyleTLatex( labelRatioToFitEtaBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitEtaBinShift->Draw();

            canvasShift->Update();
            canvasShift->SaveAs(Form("%s/BinShiftCorrection_Eta_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

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
            canvasDummy2->Print(Form("%s/ComparisonShiftedEta_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete canvasDummy2;
        }

        // *************************************************************************************************************
        // redo fitting after binshifts
        // *************************************************************************************************************
        // Tsallis function
        cout << WriteParameterToFile(fitInvYieldEta[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitInvYieldEta[cent])<< endl;
        //Two component model from Bylinkin
        fitTCMInvYieldEta[cent]        = FitObject("tcm","fitTCMInvYieldEtapPb5TeVCent","Eta",graphCombEtaInvYieldTot[cent],
                                                   minPtEta[cent], xPtLimitsEta[cent][maxNBinsEta[cent]] ,paramTCMEtaNew,"QNRMEX0+","", kFALSE);
        if (centArray[cent].Contains("0-20%") ){
            fitTCMInvYieldEta[cent]->SetParameter(0, 2);
            fitTCMInvYieldEta[cent]->SetParLimits(0, 0, 10);
            fitTCMInvYieldEta[cent]->SetParameter(1, 0.12);
            fitTCMInvYieldEta[cent]->SetParLimits(1, 0.09, 0.3);
        } else if (centArray[cent].Contains("20-40%") ){
            fitTCMInvYieldEta[cent]->SetParameter(0, 2);
            fitTCMInvYieldEta[cent]->SetParLimits(0, 0, 5);
            fitTCMInvYieldEta[cent]->SetParameter(1, 0.12);
            fitTCMInvYieldEta[cent]->SetParLimits(1, 0.1, 0.3);
        } else if (centArray[cent].Contains("40-60%") ){
            fitTCMInvYieldEta[cent]->SetParameter(0, 1);
            fitTCMInvYieldEta[cent]->SetParLimits(0, 0, 2);
            fitTCMInvYieldEta[cent]->SetParameter(1, 0.12);
            fitTCMInvYieldEta[cent]->SetParLimits(1, 0.1, 0.3);
        } else if (centArray[cent].Contains("60-100%") ){
            fitTCMInvYieldEta[cent]->SetParameter(0, 0.5);
            fitTCMInvYieldEta[cent]->SetParLimits(0, 0, 1);
            fitTCMInvYieldEta[cent]->SetParameter(1, 0.12);
            fitTCMInvYieldEta[cent]->SetParLimits(1, 0.1, 0.3);

        }
        graphCombEtaInvYieldStat[cent]->Fit(fitTCMInvYieldEta[cent],"QNRMEX0+" );

        cout << WriteParameterToFile(fitTCMInvYieldEta[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMInvYieldEta[cent])<< endl;

        fitPowInvYieldEta[cent]        = FitObject("m","fitPowInvYieldEtapPb5TeVCent","Eta",graphCombEtaInvYieldTot[cent], 5, 40. ,NULL,"QNRMEX0+","", kFALSE);
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
        if (histoDPMJetEta[cent]){
            histoRatioEtaDPMJetToFit[cent]                     = (TH1D*) histoDPMJetEta[cent]->Clone(Form("histoRatioEtaDPMJetToFit%s%s",centArray[cent].Data(), runArray[cent].Data() ));
            histoRatioEtaDPMJetToFit[cent]                     = CalculateHistoRatioToFit (histoRatioEtaDPMJetToFit[cent], fitTCMInvYieldEta[cent]);
            histoRatioEtaDPMJetToFit[cent]->GetXaxis()->SetRangeUser(minPtEta[cent], xPtLimitsEta[cent][maxNBinsEta[cent]]);
        }
        if (histoHIJINGEta[cent]){
            histoRatioEtaHIJINGToFit[cent]                     = (TH1D*) histoHIJINGEta[cent]->Clone(Form("histoRatioEtaHIJINGToFit%s%s",centArray[cent].Data(), runArray[cent].Data() ));
            histoRatioEtaHIJINGToFit[cent]                     = CalculateHistoRatioToFit (histoRatioEtaHIJINGToFit[cent], fitTCMInvYieldEta[cent]);
            histoRatioEtaHIJINGToFit[cent]->GetXaxis()->SetRangeUser(minPtEta[cent], xPtLimitsEta[cent][maxNBinsEta[cent]]);
        }
        if (histoEPOSLHCEta[cent]){
            histoRatioEtaEPOSLHCToFit[cent]                    = (TH1D*) histoEPOSLHCEta[cent]->Clone(Form("histoRatioEtaEPOSLHCToFit%s%s",centArray[cent].Data(), runArray[cent].Data() ));
            histoRatioEtaEPOSLHCToFit[cent]                    = CalculateHistoRatioToFit (histoRatioEtaEPOSLHCToFit[cent], fitTCMInvYieldEta[cent]);
            histoRatioEtaEPOSLHCToFit[cent]->GetXaxis()->SetRangeUser(minPtEta[cent], xPtLimitsEta[cent][maxNBinsEta[cent]]);
        }

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
        graphCombEtaInvYieldStatWOXErr[cent]                    = (TGraphAsymmErrors*)graphCombEtaInvYieldStat[cent]->Clone("graphCombEtaInvYieldStatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphCombEtaInvYieldStatWOXErr[cent]);
        graphRatioEtaCombCombFitStatWOXErr[cent]                = (TGraphAsymmErrors*)graphRatioEtaCombCombFitStat[cent]->Clone("graphRatioEtaCombCombFitStatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombCombFitStatWOXErr[cent]);

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

            TLatex *labelRatioToFitEnergy   = new TLatex(0.95, 0.92, Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
            labelRatioToFitEnergy->Draw();
            TLatex *labelRatioToFitALICE    = new TLatex(0.95, 0.86, textALICE.Data());
            SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
            labelRatioToFitALICE->Draw();
            TLatex *labelRatioToFitEta      = new TLatex(0.12, 0.92, "#eta #rightarrow #gamma#gamma");
            SetStyleTLatex( labelRatioToFitEta, textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
            labelRatioToFitEta->Draw();

        canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfCombToCombFit_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

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
            Double_t columnsLegendOnlyEtaRatioAbs[6]    = {0.215, 1.22, 1.68, 2, 11.1, 15.6};
            Double_t columnsLegendOnlyEtaRatioAbs2[6]   = {0.215, 1.22, 2.25, 2, 11.1, 21.3};
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
                TBox* boxPCMEtaOnlyRatioEta                 = CreateBoxFromGraph(graphRatioEtaIndCombFitSys[cent][0], columnsLegendOnlyEtaRatioAbs[2] , rowsLegendOnlyEtaRatioAbs[1]- heightBox,
                                                                                columnsLegendOnlyEtaRatioAbs2[2], rowsLegendOnlyEtaRatioAbs[1]+ heightBox);
                boxPCMEtaOnlyRatioEta->Draw("l");
            }
            if (graphRatioEtaIndCombFitSys[cent][1]){
                TBox* boxPHOSEtaOnlyRatioEta                = CreateBoxFromGraph(graphRatioEtaIndCombFitSys[cent][1], columnsLegendOnlyEtaRatioAbs[2], rowsLegendOnlyEtaRatioAbs[2]- heightBox,
                                                                                columnsLegendOnlyEtaRatioAbs2[2], rowsLegendOnlyEtaRatioAbs[2]+ heightBox);
                boxPHOSEtaOnlyRatioEta->Draw("l");
            }
            if (graphRatioEtaIndCombFitSys[cent][2]){
                TBox* boxEMCALEtaOnlyRatioEta               = CreateBoxFromGraph(graphRatioEtaIndCombFitSys[cent][2], columnsLegendOnlyEtaRatioAbs[2], rowsLegendOnlyEtaRatioAbs[3]- heightBox,
                                                                                columnsLegendOnlyEtaRatioAbs2[2], rowsLegendOnlyEtaRatioAbs[3]+ heightBox);
                boxEMCALEtaOnlyRatioEta->Draw("l");
            }
            if (graphRatioEtaIndCombFitSys[cent][4]){
                TBox* boxPCMEMCALEtaOnlyRatioEta            = CreateBoxFromGraph(graphRatioEtaIndCombFitSys[cent][4], columnsLegendOnlyEtaRatioAbs[5], rowsLegendOnlyEtaRatioAbs[1]- heightBox,
                                                                                columnsLegendOnlyEtaRatioAbs2[5], rowsLegendOnlyEtaRatioAbs[1]+ heightBox);
                boxPCMEMCALEtaOnlyRatioEta->Draw("l");
            }
            if (graphRatioEtaIndCombFitSys[cent][3]){
                TBox* boxPCMPHOSEtaOnlyRatioEta             = CreateBoxFromGraph(graphRatioEtaIndCombFitSys[cent][3], columnsLegendOnlyEtaRatioAbs[5], rowsLegendOnlyEtaRatioAbs[2]- heightBox,
                                                                                columnsLegendOnlyEtaRatioAbs2[5], rowsLegendOnlyEtaRatioAbs[2]+ heightBox);
                boxPCMPHOSEtaOnlyRatioEta->Draw("l");
            }
        canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfIndividualMeasToCombFit_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(),
                                          runArray[cent].Data() ,suffix.Data()));
    }


        // *******************************************************************************************************
    // ************************** Combination of different eta measurements **********************************
    // *******************************************************************************************************
    // REMARKS:
    //     - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - extendable to other energies
    //     - offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented

    TH1D* statErrorCollectionEtaToPi0[17][11];
    TGraphAsymmErrors* statErrorGraphCollectionEtaToPi0[17][11];
    TGraphAsymmErrors* sysErrorCollectionEtaToPi0[17][11];
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCent[cent]) continue;
        for (Int_t meth = 0; meth< 11; meth++){
            // definition of array of histograms (NULL - means we have no measurement at this energy for this rec-method)
            // for statistical error and final value from respective method
            statErrorCollectionEtaToPi0[cent][meth]          = NULL;
            cout <<cent << "\t"<< meth << "\t" <<  histoEtaToPi0Stat[cent][meth] << endl;
            if (histoEtaToPi0Stat[cent][meth]) statErrorCollectionEtaToPi0[cent][meth]      = (TH1D*)histoEtaToPi0Stat[cent][meth]->Clone(  Form("statErr%i_%sEtaToPi0",
                                                                                                                                                    cent,nameMeasGlobalLabel[meth].Data()));
            statErrorGraphCollectionEtaToPi0[cent][meth]     = NULL;
            if (graphEtaToPi0Stat[cent][meth]) statErrorGraphCollectionEtaToPi0[cent][meth] = (TGraphAsymmErrors*)graphEtaToPi0Stat[cent][meth]->Clone( Form("statErrGraph%i_%sEtaToPi0", cent,
                                                                                                                                                                nameMeasGlobalLabel[meth].Data()));
            // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
            // for systematic error from respective method
            sysErrorCollectionEtaToPi0[cent][meth]           = NULL;
            if (graphEtaToPi0Sys[cent][meth]) sysErrorCollectionEtaToPi0[cent][meth]        = (TGraphAsymmErrors*)graphEtaToPi0Sys[cent][meth]->Clone(  Form("sysErr%i_%sEtaToPi0", cent,
                                                                                                                                                                nameMeasGlobalLabel[meth].Data()));
            if (graphEtaToPi0Stat[cent][meth] && graphEtaToPi0Sys[cent][meth]){
                while (graphEtaToPi0Stat[cent][meth]->GetY()[0] == 0 && graphEtaToPi0Stat[cent][meth]->GetN() > 0){
                    graphEtaToPi0Stat[cent][meth]->RemovePoint(0);
                    graphEtaToPi0Sys[cent][meth]->RemovePoint(0);
                }
                while (graphEtaToPi0Stat[cent][meth]->GetY()[graphEtaToPi0Stat[cent][meth]->GetN()-1] == 0 && graphEtaToPi0Stat[cent][meth]->GetN() > 0){
                    graphEtaToPi0Stat[cent][meth]->RemovePoint(graphEtaToPi0Stat[cent][meth]->GetN()-1);
                    graphEtaToPi0Sys[cent][meth]->RemovePoint(graphEtaToPi0Sys[cent][meth]->GetN()-1);
                }
            }
        }
    }

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for EtaToPi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsEtaToPi0[17][11]           = {
                                                { 0,  9,  0,  1,  0,  0,  0,  0,  0,  0,  0 },  // MB R1
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // MB R1 cent binning
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-20
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 20-40
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 40-60
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 60-100
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-20 CL1
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 20-40 CL1
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 40-60 CL1
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 60-100 CL1
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-20 ZNA
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 20-40 ZNA
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 40-60 ZNA
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 60-100 ZNA
                                                { 0,  7,  10,  5,  7,  0,  0,  0,  0,  0,  0 },  // MB R2
                                                { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-10
                                                { 0,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0 }  // 60-80
                                              };
    Int_t offSetsEtaToPi0Sys[17][11]        = {
                                                { 3,  9,  7,  5,  5,  0,  0,  0,  0,  0, 0 },  // MB R1
                                                { 3,  7,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // MB R1 cent binning
                                                { 3,  7,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 0-20
                                                { 3,  7,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 20-40
                                                { 3,  7,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 40-60
                                                { 3,  7,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 60-100
                                                { 3,  7,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 0-20 CL1
                                                { 3,  7,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 20-40 CL1
                                                { 3,  7,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 40-60 CL1
                                                { 3,  7,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 60-100 CL1
                                                { 3,  7,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 0-20 ZNA
                                                { 3,  7,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 20-40 ZNA
                                                { 3,  7,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 40-60 ZNA
                                                { 3,  7,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 60-100 ZNA
                                                { 2,  11,  13,  8,  10,  0,  0,  0,  0,  0, 0 },  // MB R2
                                                { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 },  // 0-10
                                                { 3,  0,  6,  4,  5,  0,  0,  0,  0,  0, 0 }  // 60-80
                                              };
    Double_t minPtEtaToPi0[17]              = { 0.7, 0.7, 0.7, 0.7, 0.7,   0.7, 0.7, 0.7, 0.7, 0.7,
                                                0.7, 0.7, 0.7, 0.7, 0.4,   0.4, 0.4};

    TGraphAsymmErrors* statErrorRelCollectionEtaToPi0[17][11];
    TGraphAsymmErrors* sysErrorRelCollectionEtaToPi0[17][11];
    TGraph* graphWeightsEtaToPi0[17][11];
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        for (Int_t meth = 0; meth< 11; meth++){
            graphWeightsEtaToPi0[cent][meth]                 = NULL;
            statErrorRelCollectionEtaToPi0[cent][meth]       = NULL;
            sysErrorRelCollectionEtaToPi0[cent][meth]        = NULL;
        }
    }

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCent[cent]) continue;
        for (Int_t meth = 0; meth< 11; meth++){
            if (statErrorCollectionEtaToPi0[cent][meth]){
                statErrorRelCollectionEtaToPi0[cent][meth]   = new TGraphAsymmErrors(statErrorCollectionEtaToPi0[cent][meth]);
                RemoveZerosAtBeginningAndEndFromGraph(statErrorRelCollectionEtaToPi0[cent][meth]);                statErrorRelCollectionEtaToPi0[cent][meth]   = CalculateRelErrUpAsymmGraph( statErrorRelCollectionEtaToPi0[cent][meth], Form("relativeStatErrorEtaToPi0_%i_%s", cent,nameMeasGlobal[meth].Data()));
            }

            if (sysErrorCollectionEtaToPi0[cent][meth]){
                sysErrorRelCollectionEtaToPi0[cent][meth]   = (TGraphAsymmErrors*)sysErrorCollectionEtaToPi0[cent][meth]->Clone(Form("relativeSysErrorEtaToPi0%i_%s", cent, nameMeasGlobal[meth].Data()));
                RemoveZerosAtBeginningAndEndFromGraph(sysErrorRelCollectionEtaToPi0[cent][meth]);
                sysErrorRelCollectionEtaToPi0[cent][meth]    = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionEtaToPi0[cent][meth], Form("relativeSysErrorEtaToPi0%i_%s", cent, nameMeasGlobal[meth].Data()));
            }
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************

    TGraphAsymmErrors* graphCombEtaToPi0Stat[17]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaToPi0StatWOXErr[17]       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaToPi0Sys[17]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaToPi0Tot[17]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaToPi0RelStat[17]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaToPi0RelSys[17]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };
    TGraphAsymmErrors* graphCombEtaToPi0RelTot[17]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL  };

    TGraphAsymmErrors* graphIndEtaToPi0Stat[17][11];
    TGraphAsymmErrors* graphIndEtaToPi0Sys[17][11];
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        for (Int_t meth = 0; meth< 11; meth++){
            graphIndEtaToPi0Stat[cent][meth]                     = NULL;
            graphIndEtaToPi0Sys[cent][meth]                      = NULL;
        }
    }

    // *************************************************************************************************************
    // ************************************** Define plotting environment ******************************************
    // *************************************************************************************************************
    textSizeLabelsPixel                 = 900*0.04;
    TH2F * histo2DEtaToPi0Weights;
    histo2DEtaToPi0Weights = new TH2F("histo2DEtaToPi0Weights","histo2DEtaToPi0Weights",11000,minPtEtaToPi0Plotting, maxPtEtaToPi0Plotting,1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DEtaToPi0Weights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DEtaToPi0Weights->GetXaxis()->SetMoreLogLabels();
    histo2DEtaToPi0Weights->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelWeightsEtaToPi0         = new TLatex(0.95,0.15,"#eta/#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelWeightsEtaToPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    TLatex *labelRelSysErrEtaToPi0       = new TLatex(0.15,0.85,"#eta/#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrEtaToPi0, textSizeLabelsPixel, 4, 1, 43);

    TLatex *labelRelStatErrEtaToPi0      = new TLatex(0.95,0.85,"#eta/#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelStatErrEtaToPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    TH2F * histo2DRelTotErrEtaToPi0;
    histo2DRelTotErrEtaToPi0                 = new TH2F("histo2DRelTotErrEtaToPi0","histo2DRelTotErrEtaToPi0",11000,minPtEtaPlotting, maxPtEtaPlotting,1000,0,50.0);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrEtaToPi0->GetYaxis()->SetRangeUser(0,39.5);
    histo2DRelTotErrEtaToPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrEtaToPi0->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelRelTotErrEtaToPi0       = new TLatex(0.95,0.85,"#eta/#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrEtaToPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCentComb[cent]) continue;
        // Declaration & calculation of combined spectrum
        cout << "start combining eta for " << centArray[cent].Data() << " " << runArray[cent].Data() << endl;
        TString fileNameEtaToPi0OutputWeighting         = Form("%s/EtaToPi0_WeightingMethod_%s%s%s.dat",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data());
        graphCombEtaToPi0Tot[cent]                      = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionEtaToPi0[cent],    sysErrorCollectionEtaToPi0[cent],
                                                                                                xPtLimitsEtaToPi0[cent], maxNBinsEtaToPi0[cent],
                                                                                                offSetsEtaToPi0[cent], offSetsEtaToPi0Sys[cent],
                                                                                                graphCombEtaToPi0Stat[cent], graphCombEtaToPi0Sys[cent],
                                                                                                fileNameEtaToPi0OutputWeighting, "pPb_5.023TeV", "EtaToPi0", kTRUE,
                                                                                                NULL, fileNameCorrFactors, centArrayCorr[cent]
                                                                                            );


        if (graphCombEtaToPi0Tot[cent] == NULL) {
            cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
            return;
        }
        while (graphCombEtaToPi0Stat[cent]->GetX()[0] < minPtEtaToPi0[cent]){
            graphCombEtaToPi0Stat[cent]->RemovePoint(0);
        }
        while (graphCombEtaToPi0Tot[cent]->GetX()[0] < minPtEtaToPi0[cent]){
            graphCombEtaToPi0Tot[cent]->RemovePoint(0);
        }
        while (graphCombEtaToPi0Sys[cent]->GetX()[0] < minPtEtaToPi0[cent]){
            graphCombEtaToPi0Sys[cent]->RemovePoint(0);
        }
//         graphCombEtaToPi0Tot[cent]->Print();

        // ***************************************************************************************************************
        // ******************************** removing X-errors from graphs ************************************************
        // ***************************************************************************************************************
        for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
            if (!enableCentComb[cent]) continue;
            if (graphCombEtaToPi0Stat[cent]){
                graphCombEtaToPi0StatWOXErr[cent]           = (TGraphAsymmErrors*)graphCombEtaToPi0Stat[cent]->Clone(Form("EtaToPi0StatW0XErr%s%s", centArray[cent].Data(), runArray[cent].Data()));
                ProduceGraphAsymmWithoutXErrors(graphCombEtaToPi0StatWOXErr[cent]);
            }
        }

        // Reading weights from output file for plotting
        ifstream fileWeightsEtaToPi0Read;
        fileWeightsEtaToPi0Read.open(fileNameEtaToPi0OutputWeighting,ios_base::in);
        cout << "reading: " << fileNameEtaToPi0OutputWeighting << endl;
        Double_t xValuesEtaToPi0Read[50];
        Double_t weightsEtaToPi0Read[11][50];
        Int_t availableEtaToPi0Meas[11]    = {   -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1,
                                            -1};
        Int_t nMeasSetEtaToPi0             = nTotMeasEtaToPi0[cent];
        Int_t nPtBinsEtaToPi0Read          = 0;
        while(!fileWeightsEtaToPi0Read.eof() && nPtBinsEtaToPi0Read < 50){
            TString garbage             = "";
            if (nPtBinsEtaToPi0Read == 0){
                fileWeightsEtaToPi0Read >> garbage ;//>> availableEtaToPi0Meas[0] >> availableEtaToPi0Meas[1] >> availableEtaToPi0Meas[2] >> availableEtaToPi0Meas[3];
                for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                    fileWeightsEtaToPi0Read >> availableEtaToPi0Meas[i] ;
                }
                cout << "read following measurements: ";
                for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
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
        nPtBinsEtaToPi0Read                  = nPtBinsEtaToPi0Read-2 ;
        fileWeightsEtaToPi0Read.close();

        for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
            graphWeightsEtaToPi0[cent][availableEtaToPi0Meas[i]]                        = new TGraph(nPtBinsEtaToPi0Read,xValuesEtaToPi0Read,weightsEtaToPi0Read[availableEtaToPi0Meas[i]]);
            Int_t bin = 0;
            for (Int_t n = 0; n< nPtBinsEtaToPi0Read; n++){
                if (graphWeightsEtaToPi0[cent][availableEtaToPi0Meas[i]]->GetY()[bin] == 0) graphWeightsEtaToPi0[cent][availableEtaToPi0Meas[i]]->RemovePoint(bin);
                else bin++;
            }
//             cout << "weights" << cent << "\t" << availableEtaToPi0Meas[i] << endl;
//             graphWeightsEtaToPi0[cent][availableEtaToPi0Meas[i]]->Print();
        }


        // **********************************************************************************************************************
        // ******************************************* Plotting weights method only EMC *****************************************
        // **********************************************************************************************************************
        textSizeLabelsPixel                 = 900*0.04;
        canvasWeights->cd();
        histo2DEtaToPi0Weights->Draw("copy");

            TLegend* legendWeights   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEtaToPi0), 32);
            for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                DrawGammaSetMarkerTGraph(graphWeightsEtaToPi0[cent][availableEtaToPi0Meas[i]], markerStyleDet[availableEtaToPi0Meas[i]], markerSizeDet[availableEtaToPi0Meas[i]]*0.5, colorDet[availableEtaToPi0Meas[i]] , colorDet[availableEtaToPi0Meas[i]]);
                graphWeightsEtaToPi0[cent][availableEtaToPi0Meas[i]]->Draw("p,same,z");
                legendWeights->AddEntry(graphWeightsEtaToPi0[cent][availableEtaToPi0Meas[i]],nameMeasGlobalLabel[availableEtaToPi0Meas[i]],"p");
            }
            legendWeights->Draw();

            labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelWeightsEnergy->Draw();
            labelWeightsEtaToPi0->Draw();

            DrawGammaLines(minPtEtaToPi0Plotting, maxPtEtaToPi0Plotting , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(minPtEtaToPi0Plotting, maxPtEtaToPi0Plotting , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(minPtEtaToPi0Plotting, maxPtEtaToPi0Plotting , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(minPtEtaToPi0Plotting, maxPtEtaToPi0Plotting , 0.2, 0.2,0.1, kGray, 3);

        canvasWeights->SaveAs(Form("%s/EtaToPi0_Weights_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

        //  *********************************************************************************************************************
        //  ************************************ Visualize relative errors ******************************************************
        //  *********************************************************************************************************************

        canvasRelSysErr->cd();
        histo2DRelSysErr->GetYaxis()->SetRangeUser(0,39.5);
        histo2DRelSysErr->Draw("copy");

            TLegend* legendRelSysErrEtaToPi0        = GetAndSetLegend2(0.62, 0.92-(0.04*nMeasSetEtaToPi0), 0.95, 0.92, textSizeLabelsPixel);
            for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                DrawGammaSetMarkerTGraph(sysErrorRelCollectionEtaToPi0[cent][availableEtaToPi0Meas[i]], markerStyleDet[availableEtaToPi0Meas[i]], markerSizeDet[availableEtaToPi0Meas[i]]*0.5, colorDet[availableEtaToPi0Meas[i]],
                                         colorDet[availableEtaToPi0Meas[i]]);
                sysErrorRelCollectionEtaToPi0[cent][availableEtaToPi0Meas[i]]->Draw("p,same,z");
                legendRelSysErrEtaToPi0->AddEntry(sysErrorRelCollectionEtaToPi0[cent][availableEtaToPi0Meas[i]],nameMeasGlobalLabel[availableEtaToPi0Meas[i]],"p");
            }
            legendRelSysErrEtaToPi0->Draw();

            labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelSysErrEnergy->Draw();
            labelRelSysErrEtaToPi0->Draw();

        canvasRelSysErr->SaveAs(Form("%s/EtaToPi0_RelSysErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
        delete legendRelSysErrEtaToPi0;
        //  *********************************************************************************************************************
        //  ************************************ Visualize relative errors ******************************************************
        //  *********************************************************************************************************************

        canvasRelStatErr->cd();
        histo2DRelStatErr->GetYaxis()->SetRangeUser(0,39.5);
        histo2DRelStatErr->Draw("copy");
            TLegend* legendRelStatErr       = GetAndSetLegend2(0.14, 0.92-(0.04*nMeasSetEtaToPi0), 0.45, 0.92, textSizeLabelsPixel);
            for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                DrawGammaSetMarkerTGraph(statErrorRelCollectionEtaToPi0[cent][availableEtaToPi0Meas[i]], markerStyleDet[availableEtaToPi0Meas[i]], markerSizeDet[availableEtaToPi0Meas[i]]*0.5, colorDet[availableEtaToPi0Meas[i]],
                                         colorDet[availableEtaToPi0Meas[i]]);
                statErrorRelCollectionEtaToPi0[cent][availableEtaToPi0Meas[i]]->Draw("p,same,z");
                legendRelStatErr->AddEntry(statErrorRelCollectionEtaToPi0[cent][availableEtaToPi0Meas[i]],nameMeasGlobalLabel[availableEtaToPi0Meas[i]],"p");
            }
            legendRelStatErr->Draw();

            labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelStatErrEnergy->Draw();
            labelRelStatErrEtaToPi0->Draw();

        canvasRelStatErr->SaveAs(Form("%s/EtaToPi0_RelStatErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

        //  *********************************************************************************************************************
        //  ************************ Visualize relative total errors of different combination methods EtaToPi0 ***********************
        //  *********************************************************************************************************************

        graphCombEtaToPi0RelStat[cent]   = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0Stat[cent], "relativeStatErrorEtaToPi0_Method");
        graphCombEtaToPi0RelSys[cent]    = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0Sys[cent], "relativeSysErrorEtaToPi0_Method");
        graphCombEtaToPi0RelTot[cent]    = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0Tot[cent], "relativeTotalErrorEtaToPi0_Method");

        canvasRelTotErr->cd();
        histo2DRelTotErrEtaToPi0->GetYaxis()->SetRangeUser(0,39.5);
        histo2DRelTotErrEtaToPi0->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelTot[cent], markerStyleComb, markerSizeComb, kBlue+2 , kBlue+2);
            graphCombEtaToPi0RelTot[cent]->Draw("p,same,z");

            TLegend* legendRelTotErr1       = GetAndSetLegend2(0.20, 0.92-(0.035), 0.45, 0.92, 32);
            legendRelTotErr1->AddEntry(graphCombEtaToPi0RelTot[cent],"All","p");
            legendRelTotErr1->Draw();


            labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelTotErrEnergy->Draw();
            labelRelTotErrEtaToPi0->Draw();

        canvasRelTotErr->SaveAs(Form("%s/EtaToPi0_TotErr_Comp_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
        histo2DRelTotErrEtaToPi0->GetYaxis()->SetTitle("Err (%)");
        histo2DRelTotErrEtaToPi0->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphCombEtaToPi0RelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelStat[cent], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
            graphCombEtaToPi0RelStat[cent]->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelSys[cent], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
            graphCombEtaToPi0RelSys[cent]->SetLineStyle(7);
            graphCombEtaToPi0RelSys[cent]->Draw("l,x0,same,e1");

            TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
            legendRelTotErr3->AddEntry(graphCombEtaToPi0RelTot[cent],"tot","p");
            legendRelTotErr3->AddEntry(graphCombEtaToPi0RelStat[cent],"stat","l");
            legendRelTotErr3->AddEntry(graphCombEtaToPi0RelSys[cent],"sys","l");
            legendRelTotErr3->Draw();

            labelRelTotErrEnergy->Draw();
            labelRelTotErrEtaToPi0->Draw();

        canvasRelTotErr->SaveAs(Form("%s/EtaToPi0_Reldecomp_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

    }

    TGraphAsymmErrors* graphRpPbIndStatPi0[17][11];
    TGraphAsymmErrors* graphRpPbIndStatPi0WOXErr[17][11];
    TH1D* histRpPbIndStatPi0[17][11];
    TGraphAsymmErrors* graphRpPbIndSystPi0[17][11];
    TGraphAsymmErrors* graphRpPbIndCombPi0[17][11];
    TGraph* graphWeightsPi0RpPb[17][11];
    TGraphAsymmErrors* statErrorRelCollectionPi0RpA[17][11];
    TGraphAsymmErrors* sysErrorRelCollectionPi0RpA[17][11];
    TGraphAsymmErrors* graphRpPbCombStatPi0[17];
    TGraphAsymmErrors* graphRpPbCombSystPi0[17];
    TGraphAsymmErrors* graphRpPbCombStatPi0WOXErr[17];
    TGraphAsymmErrors* graphRpPbCombCombPi0[17];
    TGraphAsymmErrors* graphCombPi0RpPbRelStat[17];
    TGraphAsymmErrors* graphCombPi0RpPbRelSys[17];
    TGraphAsymmErrors* graphCombPi0RpPbRelTot[17];

    // *******************************************************************************************************
    // *************************** RpPb calculation for pi0 for individual measurements **********************
    // *******************************************************************************************************
    cout << "starting RpPb calculation for pi0" << endl;

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if ( !enableCentRpPb[cent] ) continue;
        // *******************************************************************************************************
        // *************************** RpPb calculation for pi0 for individual measurements **********************
        // *******************************************************************************************************
        for (Int_t meth = 0; meth < 11; meth++){
            cout << "RpPb calc for pi0 " << centArray[cent].Data() << addCentString[cent].Data() << "\t" << runArray[cent].Data() << "\t" << nameMeasGlobal[meth].Data() << endl;
            graphWeightsPi0RpPb[cent][meth]                              = NULL;
            graphRpPbIndStatPi0[cent][meth]                              = NULL;
            graphRpPbIndStatPi0WOXErr[cent][meth]                        = NULL;
            histRpPbIndStatPi0[cent][meth]                               = NULL;
            graphRpPbIndSystPi0[cent][meth]                              = NULL;
            graphRpPbIndCombPi0[cent][meth]                              = NULL;
            statErrorRelCollectionPi0RpA[cent][meth]                     = NULL;
            sysErrorRelCollectionPi0RpA[cent][meth]                      = NULL;
            if (haveRefPPPi0[ppRefInt[cent]][meth]  ){ //&& !( i == 3 )
                if(bWCorrection.Contains("Y") ){
                    TGraphAsymmErrors* tempPPStat               = (TGraphAsymmErrors*)statErrorCollectionPi0PP[ppRefInt[cent]][meth]->Clone("tempStat");
                    TGraphAsymmErrors* tempPPSys                = (TGraphAsymmErrors*)systErrorUnCorrInterCollectionPi0PP[ppRefInt[cent]][meth]->Clone("tempSys");
                    while (tempPPStat->GetX()[0] < graphIndPi0InvYieldStat_yShifted[cent][meth]->GetX()[0])
                        tempPPStat->RemovePoint(0);
                    while (tempPPSys->GetX()[0] < graphIndPi0InvYieldSys_yShifted[cent][meth]->GetX()[0])
                        tempPPSys->RemovePoint(0);
                    while (tempPPStat->GetX()[tempPPStat->GetN()-1] > graphIndPi0InvYieldStat_yShifted[cent][meth]->GetX()[graphIndPi0InvYieldStat_yShifted[cent][meth]->GetN()-1])
                        tempPPStat->RemovePoint(tempPPStat->GetN()-1);
                    while (tempPPSys->GetX()[tempPPSys->GetN()-1] > graphIndPi0InvYieldSys_yShifted[cent][meth]->GetX()[graphIndPi0InvYieldSys_yShifted[cent][meth]->GetN()-1])
                        tempPPSys->RemovePoint(tempPPSys->GetN()-1);

                    cout << "stat pp" << endl;
                    tempPPStat->Print();
                    cout << "sys pp" << endl;
                    tempPPSys->Print();
                    cout << "stat pPb" << endl;
                    graphIndPi0InvYieldStat_yShifted[cent][meth]->Print();
                    cout << "sys pPb" << endl;
                    graphIndPi0InvYieldSys_yShifted[cent][meth]->Print();


                    graphRpPbIndCombPi0[cent][meth]                      = CalcRpPbV2(  tempPPStat,
                                                                                        tempPPSys,
                                                                                        graphIndPi0InvYieldStat_yShifted[cent][meth],
                                                                                        graphIndPi0InvYieldSys_yShifted[cent][meth],
                                                                                        &graphRpPbIndStatPi0[cent][meth],
                                                                                        &graphRpPbIndSystPi0[cent][meth],
                                                                                        tpPb[cent], tpPbErr[cent] ,
                                                                                        ptSysRemNames[cent][meth], fileNamespPbPi0DetailedSys[cent][meth], fileNamesppPi0DetailedSys[cent][meth], fileNamesRpPbPi0DetailedSys[cent][meth]
                    );
                    graphRpPbIndCombPi0[cent][meth]->Print();
                    histRpPbIndStatPi0[cent][meth]                       = (TH1D*)statErrorCollectionPi0[cent][meth]->Clone(Form("histRpPb%sStatPi0",nameMeasGlobalLabel[meth].Data()));
                    for (Int_t j = 0; j< histRpPbIndStatPi0[cent][meth]->GetNbinsX()+1; j++){
                        histRpPbIndStatPi0[cent][meth]->SetBinContent(j,0);
                        histRpPbIndStatPi0[cent][meth]->SetBinError(j,0);
                    }
                    for (Int_t j = 0; j < graphRpPbIndStatPi0[cent][meth]->GetN(); j++){
                        Int_t bin                               = histRpPbIndStatPi0[cent][meth]->GetXaxis()->FindBin(graphRpPbIndStatPi0[cent][meth]->GetX()[j]);
                        histRpPbIndStatPi0[cent][meth]->SetBinContent(bin,graphRpPbIndStatPi0[cent][meth]->GetY()[j]);
                        histRpPbIndStatPi0[cent][meth]->SetBinError(bin,graphRpPbIndStatPi0[cent][meth]->GetEYlow()[j]);
                    }
                }

                if (histRpPbIndStatPi0[cent][meth]){
                    statErrorRelCollectionPi0RpA[cent][meth]   = (TGraphAsymmErrors*)graphRpPbIndStatPi0[cent][meth]->Clone(Form("relativeStatErrorPi0RpPb%i_%s", cent, nameMeasGlobal[meth].Data()));
                    RemoveZerosAtBeginningAndEndFromGraph(statErrorRelCollectionPi0RpA[cent][meth]);
                    statErrorRelCollectionPi0RpA[cent][meth]   = CalculateRelErrUpAsymmGraph( statErrorRelCollectionPi0RpA[cent][meth], Form("relativeStatErrorPi0RpPb_%i_%s", cent,nameMeasGlobal[meth].Data()));
                    statErrorRelCollectionPi0RpA[cent][meth]->Print();
                    graphRpPbIndStatPi0WOXErr[cent][meth]      = (TGraphAsymmErrors*)graphRpPbIndStatPi0[cent][meth]->Clone(Form("graphRpPb%sStatPi0WOXErr%i", nameMeasGlobalLabel[meth].Data(), cent));

                }
                if (graphRpPbIndSystPi0[cent][meth]){
                    sysErrorRelCollectionPi0RpA[cent][meth]    = (TGraphAsymmErrors*)graphRpPbIndSystPi0[cent][meth]->Clone(Form("relativeSysErrorPi0RpPb%i_%s", cent, nameMeasGlobal[meth].Data()));
                    RemoveZerosAtBeginningAndEndFromGraph(sysErrorRelCollectionPi0RpA[cent][meth]);
                    sysErrorRelCollectionPi0RpA[cent][meth]    = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionPi0RpA[cent][meth], Form("relativeSysErrorPi0RpPb%i_%s", cent, nameMeasGlobal[meth].Data()));
                }
            }
        }

        // *******************************************************************************************************
        // *************************** RpPb calculation for pi0 **************************************************
        // *******************************************************************************************************
        graphRpPbCombStatPi0[cent]                      = NULL;
        graphRpPbCombStatPi0[cent]                      = NULL;
        graphRpPbCombStatPi0[cent]                      = NULL;

        if (!enableCentComb[cent]) continue;
        if(bWCorrection.Contains("Y") ){
            TString fileNamePi0RpPbOutputWeighting      = Form("%s/Pi0RpPb_WeightingMethod_%s%s%s.dat", outputDirSupportComb.Data(), centArrayOutput[cent].Data(),
                                                                addCentString[cent].Data(), runArray[cent].Data());
            graphRpPbCombCombPi0[cent]                  = CombinePtPointsSpectraFullCorrMat(    histRpPbIndStatPi0[cent],    graphRpPbIndSystPi0[cent],
                                                                                                xPtLimitsPi0[cent], maxNBinsPi0[cent],
                                                                                                offSetsPi0[cent], offSetsPi0Sys[cent],
                                                                                                graphRpPbCombStatPi0[cent], graphRpPbCombSystPi0[cent],
                                                                                                fileNamePi0RpPbOutputWeighting, "pPb_5.023TeV", "Pi0RpPb", kTRUE,
                                                                                                NULL, fileNameCorrFactors, centArrayCorr[cent]
                                                                                            );
            if (graphRpPbCombCombPi0 == NULL) {
                cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
                return;
            }
            while (graphRpPbCombStatPi0[cent]->GetX()[0] < minPtPi0[cent]){
                graphRpPbCombStatPi0[cent]->RemovePoint(0);
            }
            while (graphRpPbCombSystPi0[cent]->GetX()[0] < minPtPi0[cent]){
                graphRpPbCombSystPi0[cent]->RemovePoint(0);
            }
            while (graphRpPbCombCombPi0[cent]->GetX()[0] < minPtPi0[cent]){
                graphRpPbCombCombPi0[cent]->RemovePoint(0);
            }
            graphRpPbCombCombPi0[cent]->Print();

            graphRpPbCombStatPi0WOXErr[cent]            = (TGraphAsymmErrors*)graphRpPbCombStatPi0[cent]->Clone("graphRpPbCombStatPi0WOXErr");
            ProduceGraphAsymmWithoutXErrors(graphRpPbCombStatPi0WOXErr[cent]);


            // Reading weights from output file for plotting
            ifstream fileWeightsPi0RpPbRead;
            fileWeightsPi0RpPbRead.open(fileNamePi0RpPbOutputWeighting,ios_base::in);
            cout << "reading" << fileNamePi0RpPbOutputWeighting << endl;
            Double_t xValuesPi0RpPbRead[50];
            Double_t weightsPi0RpPbRead[11][50];
            Int_t availablePi0RpPbMeas[11]      = { -1, -1, -1, -1, -1,
                                                -1, -1, -1, -1, -1,
                                                -1};
            Int_t nMeasSetPi0RpPb               = nTotMeasPi0[cent];
            Int_t nPtBinsPi0RpPbRead            = 0;
            while(!fileWeightsPi0RpPbRead.eof() && nPtBinsPi0RpPbRead < 50){
                TString garbage             = "";
                if (nPtBinsPi0RpPbRead == 0){
                    fileWeightsPi0RpPbRead >> garbage ;//>> availablePi0RpPbMeas[0] >> availablePi0RpPbMeas[1] >> availablePi0RpPbMeas[2] >> availablePi0RpPbMeas[3];
                    for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
                        fileWeightsPi0RpPbRead >> availablePi0RpPbMeas[i] ;
                    }
                    cout << "read following measurements: ";
                    for (Int_t meth= 0; meth < 11; meth++){
                        cout << availablePi0RpPbMeas[meth] << "\t" ;
                    }
                    cout << endl;
                } else {
                    fileWeightsPi0RpPbRead >> xValuesPi0RpPbRead[nPtBinsPi0RpPbRead-1];
                    for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
                        fileWeightsPi0RpPbRead >> weightsPi0RpPbRead[availablePi0RpPbMeas[i]][nPtBinsPi0RpPbRead-1] ;
                    }
                    cout << "read: "<<  nPtBinsPi0RpPbRead << "\t"<< xValuesPi0RpPbRead[nPtBinsPi0RpPbRead-1] << "\t" ;
                    for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
                        cout << weightsPi0RpPbRead[availablePi0RpPbMeas[i]][nPtBinsPi0RpPbRead-1] << "\t";
                    }
                    cout << endl;
                }
                nPtBinsPi0RpPbRead++;
            }
            nPtBinsPi0RpPbRead                  = nPtBinsPi0RpPbRead-2 ;
            fileWeightsPi0RpPbRead.close();

            for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
                graphWeightsPi0RpPb[cent][availablePi0RpPbMeas[i]]      = new TGraph(nPtBinsPi0RpPbRead,xValuesPi0RpPbRead,weightsPi0RpPbRead[availablePi0RpPbMeas[i]]);
                Int_t bin = 0;
                for (Int_t n = 0; n< nPtBinsPi0RpPbRead; n++){
                    if (graphWeightsPi0RpPb[cent][availablePi0RpPbMeas[i]]->GetY()[bin] == 0) graphWeightsPi0RpPb[cent][availablePi0RpPbMeas[i]]->RemovePoint(bin);
                    else bin++;
                }
            }

            // **********************************************************************************************************************
            // ******************************************* Plotting weights method only EMC *****************************************
            // **********************************************************************************************************************
            textSizeLabelsPixel           = 900*0.04;
            canvasWeights->cd();

            histo2DPi0Weights->Draw("copy");
                TLegend* legendWeightsPi0RpPb   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetPi0RpPb+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
                for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
                    DrawGammaSetMarkerTGraph(graphWeightsPi0RpPb[cent][availablePi0RpPbMeas[i]], markerStyleDet[availablePi0RpPbMeas[i]], markerSizeDet[availablePi0RpPbMeas[i]]*0.5, colorDet[availablePi0RpPbMeas[i]] , colorDet[availablePi0RpPbMeas[i]]);
                    graphWeightsPi0RpPb[cent][availablePi0RpPbMeas[i]]->Draw("p,same,z");
                    legendWeightsPi0RpPb->AddEntry(graphWeightsPi0RpPb[cent][availablePi0RpPbMeas[i]],nameMeasGlobalLabel[availablePi0RpPbMeas[i]],"p");
                }
                legendWeightsPi0RpPb->Draw();

                labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
                labelWeightsEnergy->Draw();
                labelWeightsPi0RpPb->Draw();

                DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
                DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
                DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
                DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);

            canvasWeights->SaveAs(Form("%s/Pi0RpPb_Weights_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelSysErr->cd();

                histo2DRelSysErr->GetYaxis()->SetRangeUser(0,49.5);
                histo2DRelSysErr->Draw("copy");
                TLegend* legendRelSysErrRpPb       = GetAndSetLegend2(0.60, 0.92-(0.04*(nMeasSetPi0RpPb)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
                for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
                    if (!sysErrorRelCollectionPi0RpA[cent][availablePi0RpPbMeas[i]]) continue;
                    DrawGammaSetMarkerTGraph(sysErrorRelCollectionPi0RpA[cent][availablePi0RpPbMeas[i]], markerStyleDet[availablePi0RpPbMeas[i]], markerSizeDet[availablePi0RpPbMeas[i]]*0.5,
                                            colorDet[availablePi0RpPbMeas[i]], colorDet[availablePi0RpPbMeas[i]]);
                    sysErrorRelCollectionPi0RpA[cent][availablePi0RpPbMeas[i]]->Draw("p,same,z");
                    legendRelSysErrRpPb->AddEntry(sysErrorRelCollectionPi0RpA[cent][availablePi0RpPbMeas[i]],nameMeasGlobalLabel[availablePi0RpPbMeas[i]],"p");
                }
                legendRelSysErrRpPb->Draw();

                labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
                labelRelSysErrEnergy->Draw();
                labelRelSysErrPi0RpPb->Draw();

            canvasRelSysErr->SaveAs(Form("%s/Pi0RpPb_RelSysErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelSysErrRpPb;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelStatErr->cd();

                histo2DRelStatErr->GetYaxis()->SetRangeUser(0,49.5);
                histo2DRelStatErr->Draw("copy");
                TLegend* legendRelStatErrRpPb       = GetAndSetLegend2(0.12, 0.92-(0.04*(nMeasSetPi0RpPb)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
                for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
                    if (!statErrorRelCollectionPi0RpA[cent][availablePi0RpPbMeas[i]]) continue;
                    DrawGammaSetMarkerTGraph(statErrorRelCollectionPi0RpA[cent][availablePi0RpPbMeas[i]], markerStyleDet[availablePi0RpPbMeas[i]], markerSizeDet[availablePi0RpPbMeas[i]]*0.5,
                                            colorDet[availablePi0RpPbMeas[i]], colorDet[availablePi0RpPbMeas[i]]);
                    statErrorRelCollectionPi0RpA[cent][availablePi0RpPbMeas[i]]->Draw("p,same,z");
                    legendRelStatErrRpPb->AddEntry(statErrorRelCollectionPi0RpA[cent][availablePi0RpPbMeas[i]],nameMeasGlobalLabel[availablePi0RpPbMeas[i]],"p");
                }
                legendRelStatErrRpPb->Draw();

                labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
                labelRelStatErrEnergy->Draw();
                labelRelStatErrPi0RpPb->Draw();

            canvasRelStatErr->SaveAs(Form("%s/Pi0RpPb_RelStatErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelStatErrRpPb;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            graphCombPi0RpPbRelStat[cent]      = CalculateRelErrUpAsymmGraph( graphRpPbCombStatPi0[cent], "relativeStatErrorPi0RpPb_Method");
            graphCombPi0RpPbRelSys[cent]       = CalculateRelErrUpAsymmGraph( graphRpPbCombSystPi0[cent], "relativeSysErrorPi0RpPb_Method");
            graphCombPi0RpPbRelTot[cent]       = CalculateRelErrUpAsymmGraph( graphRpPbCombCombPi0[cent], "relativeTotalErrorPi0RpPb_Method");

            canvasRelTotErr->cd();

                histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,49.5);
                histo2DRelTotErrPi0->Draw("copy");
                DrawGammaSetMarkerTGraphAsym(graphCombPi0RpPbRelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
                graphCombPi0RpPbRelTot[cent]->Draw("p,same,z");
                DrawGammaSetMarkerTGraphAsym(graphCombPi0RpPbRelStat[cent], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
                graphCombPi0RpPbRelStat[cent]->Draw("l,x0,same,e1");
                DrawGammaSetMarkerTGraphAsym(graphCombPi0RpPbRelSys[cent], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
                graphCombPi0RpPbRelSys[cent]->SetLineStyle(7);
                graphCombPi0RpPbRelSys[cent]->Draw("l,x0,same,e1");

                if (centArray[cent].BeginsWith("0-100%")){
                    DrawGammaSetMarkerTGraphAsym(graphCombPi0RpPbOlderMBRelErrTot, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
                    graphCombPi0RpPbOlderMBRelErrTot->Draw("p,same,z");
                    DrawGammaSetMarkerTGraphAsym(graphCombPi0RpPbOlderMBRelErrStat, markerStyleComb+4, markerSizeComb, kGray+2 , kGray+2);
                    graphCombPi0RpPbOlderMBRelErrStat->Draw("l,x0,same,e1");
                    DrawGammaSetMarkerTGraphAsym(graphCombPi0RpPbOlderMBRelErrSys, markerStyleComb+4, markerSizeComb, kGray+1, kGray+1);
                    graphCombPi0RpPbOlderMBRelErrSys->SetLineStyle(7);
                    graphCombPi0RpPbOlderMBRelErrSys->Draw("l,x0,same,e1");

                    TLegend* legendRelTotErr4       = GetAndSetLegend2(0.35, 0.92-(0.035*3), 0.65, 0.92, 32);
                    legendRelTotErr4->AddEntry(graphCombPi0RpPbOlderMBRelErrTot,"pub tot","p");
                    legendRelTotErr4->AddEntry(graphCombPi0RpPbOlderMBRelErrStat,"pub stat","l");
                    legendRelTotErr4->AddEntry(graphCombPi0RpPbOlderMBRelErrSys,"pub sys","l");
                    legendRelTotErr4->Draw();

                }

                TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.35, 0.92, 32);
                legendRelTotErr3->AddEntry(graphCombPi0RpPbRelTot[cent],"tot","p");
                legendRelTotErr3->AddEntry(graphCombPi0RpPbRelStat[cent],"stat","l");
                legendRelTotErr3->AddEntry(graphCombPi0RpPbRelSys[cent],"sys","l");
                legendRelTotErr3->Draw();



                labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
                labelRelTotErrEnergy->Draw();
                labelRelTotErrPi0RpPb->Draw();

            canvasRelTotErr->SaveAs(Form("%s/Pi0RpPb_Reldecomp_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelTotErr3;
        }
    }

    TGraphAsymmErrors* graphRpPbIndStatEta[17][11];
    TGraphAsymmErrors* graphRpPbIndStatEtaWOXErr[17][11];
    TH1D* histRpPbIndStatEta[17][11];
    TGraphAsymmErrors* graphRpPbIndSystEta[17][11];
    TGraphAsymmErrors* graphRpPbIndCombEta[17][11];
    TGraph* graphWeightsEtaRpPb[17][11];
    TGraphAsymmErrors* statErrorRelCollectionEtaRpA[17][11];
    TGraphAsymmErrors* sysErrorRelCollectionEtaRpA[17][11];
    TGraphAsymmErrors* graphRpPbCombStatEta[17];
    TGraphAsymmErrors* graphRpPbCombSystEta[17];
    TGraphAsymmErrors* graphRpPbCombStatEtaWOXErr[17];
    TGraphAsymmErrors* graphRpPbCombCombEta[17];
    TGraphAsymmErrors* graphCombEtaRpPbRelStat[17];
    TGraphAsymmErrors* graphCombEtaRpPbRelSys[17];
    TGraphAsymmErrors* graphCombEtaRpPbRelTot[17];

    // *******************************************************************************************************
    // *************************** RpPb calculation for eta for individual measurements **********************
    // *******************************************************************************************************
    cout << "starting RpPb calculation for eta" << endl;

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if ( !enableCentRpPb[cent] ) continue;
        // *******************************************************************************************************
        // *************************** RpPb calculation for eta for individual measurements **********************
        // *******************************************************************************************************
        for (Int_t meth = 0; meth < 11; meth++){
            cout << "RpPb calc for eta " << centArray[cent].Data() << addCentString[cent].Data() << "\t" << runArray[cent].Data() << "\t" << nameMeasGlobal[meth].Data() << endl;
            graphWeightsEtaRpPb[cent][meth]                              = NULL;
            graphRpPbIndStatEta[cent][meth]                              = NULL;
            graphRpPbIndStatEtaWOXErr[cent][meth]                        = NULL;
            histRpPbIndStatEta[cent][meth]                               = NULL;
            graphRpPbIndSystEta[cent][meth]                              = NULL;
            graphRpPbIndCombEta[cent][meth]                              = NULL;
            statErrorRelCollectionEtaRpA[cent][meth]                     = NULL;
            sysErrorRelCollectionEtaRpA[cent][meth]                      = NULL;
            if (haveRefPPEta[ppRefInt[ppRefInt[cent]]][meth]  ){ //&& !( i == 3 )
                if(bWCorrection.Contains("Y") ){
                    TGraphAsymmErrors* tempPPStat               = (TGraphAsymmErrors*)statErrorCollectionEtaPP[ppRefInt[cent]][meth]->Clone("tempStat");
                    TGraphAsymmErrors* tempPPSys                = (TGraphAsymmErrors*)systErrorUnCorrInterCollectionEtaPP[ppRefInt[cent]][meth]->Clone("tempSys");
                    while (tempPPStat->GetX()[0] < graphIndEtaInvYieldStat_yShifted[cent][meth]->GetX()[0])
                        tempPPStat->RemovePoint(0);
                    while (tempPPSys->GetX()[0] < graphIndEtaInvYieldSys_yShifted[cent][meth]->GetX()[0])
                        tempPPSys->RemovePoint(0);
                    while (tempPPStat->GetX()[tempPPStat->GetN()-1] > graphIndEtaInvYieldStat_yShifted[cent][meth]->GetX()[graphIndEtaInvYieldStat_yShifted[cent][meth]->GetN()-1])
                        tempPPStat->RemovePoint(tempPPStat->GetN()-1);
                    while (tempPPSys->GetX()[tempPPSys->GetN()-1] > graphIndEtaInvYieldSys_yShifted[cent][meth]->GetX()[graphIndEtaInvYieldSys_yShifted[cent][meth]->GetN()-1])
                        tempPPSys->RemovePoint(tempPPSys->GetN()-1);

                    cout << "stat pp" << endl;
                    tempPPStat->Print();
                    cout << "sys pp" << endl;
                    tempPPSys->Print();
                    cout << "stat pPb" << endl;
                    graphIndEtaInvYieldStat_yShifted[cent][meth]->Print();
                    cout << "sys pPb" << endl;
                    graphIndEtaInvYieldSys_yShifted[cent][meth]->Print();

                    graphRpPbIndCombEta[cent][meth]                      = CalcRpPbV2(  tempPPStat,
                                                                                        tempPPSys,
                                                                                        graphIndEtaInvYieldStat_yShifted[cent][meth],
                                                                                        graphIndEtaInvYieldSys_yShifted[cent][meth],
                                                                                        &graphRpPbIndStatEta[cent][meth],
                                                                                        &graphRpPbIndSystEta[cent][meth],
                                                                                        tpPb[cent], tpPbErr[cent] ,
                                                                                        ptSysRemNames[cent][meth], fileNamespPbEtaDetailedSys[cent][meth], fileNamesppEtaDetailedSys[cent][meth], fileNamesRpPbEtaDetailedSys[cent][meth]
                    );
                    graphRpPbIndCombEta[cent][meth]->Print();
                    histRpPbIndStatEta[cent][meth]                       = (TH1D*)statErrorCollectionEta[cent][meth]->Clone(Form("histRpPb%sStatEta",nameMeasGlobalLabel[meth].Data()));
                    for (Int_t j = 0; j< histRpPbIndStatEta[cent][meth]->GetNbinsX()+1; j++){
                        histRpPbIndStatEta[cent][meth]->SetBinContent(j,0);
                        histRpPbIndStatEta[cent][meth]->SetBinError(j,0);
                    }
                    for (Int_t j = 0; j < graphRpPbIndStatEta[cent][meth]->GetN(); j++){
                        Int_t bin                               = histRpPbIndStatEta[cent][meth]->GetXaxis()->FindBin(graphRpPbIndStatEta[cent][meth]->GetX()[j]);
                        histRpPbIndStatEta[cent][meth]->SetBinContent(bin,graphRpPbIndStatEta[cent][meth]->GetY()[j]);
                        histRpPbIndStatEta[cent][meth]->SetBinError(bin,graphRpPbIndStatEta[cent][meth]->GetEYlow()[j]);
                    }
                }

                if (histRpPbIndStatEta[cent][meth]){
                    statErrorRelCollectionEtaRpA[cent][meth]   = new TGraphAsymmErrors(histRpPbIndStatEta[cent][meth]);
                    RemoveZerosAtBeginningAndEndFromGraph(statErrorRelCollectionEtaRpA[cent][meth]);
                    statErrorRelCollectionEtaRpA[cent][meth]   = CalculateRelErrUpAsymmGraph( statErrorRelCollectionEtaRpA[cent][meth], Form("relativeStatErrorEtaRpPb_%i_%s", cent,nameMeasGlobal[meth].Data()));

                    graphRpPbIndStatEtaWOXErr[cent][meth]      = (TGraphAsymmErrors*)graphRpPbIndStatEta[cent][meth]->Clone(Form("graphRpPb%sStatEtaWOXErr%i", nameMeasGlobalLabel[meth].Data(), cent));


                }
                if (graphRpPbIndSystEta[cent][meth]){
                    sysErrorRelCollectionEtaRpA[cent][meth]    = (TGraphAsymmErrors*)graphRpPbIndSystEta[cent][meth]->Clone(Form("relativeSysErrorEtaRpPb%i_%s", cent, nameMeasGlobal[meth].Data()));
                    RemoveZerosAtBeginningAndEndFromGraph(sysErrorRelCollectionEtaRpA[cent][meth]);
                    sysErrorRelCollectionEtaRpA[cent][meth]    = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionEtaRpA[cent][meth], Form("relativeSysErrorEtaRpPb%i_%s", cent, nameMeasGlobal[meth].Data()));
                }
            }
        }

        // *******************************************************************************************************
        // *************************** RpPb calculation for pi0 **************************************************
        // *******************************************************************************************************
        graphRpPbCombStatEta[cent]                      = NULL;
        graphRpPbCombStatEta[cent]                      = NULL;
        graphRpPbCombStatEta[cent]                      = NULL;

        if (!enableCentComb[cent]) continue;
        if(bWCorrection.Contains("Y") ){
            TString fileNameEtaRpPbOutputWeighting      = Form("%s/EtaRpPb_WeightingMethod_%s%s%s.dat", outputDirSupportComb.Data(), centArrayOutput[cent].Data(),
                                                                addCentString[cent].Data(), runArray[cent].Data());
            graphRpPbCombCombEta[cent]                  = CombinePtPointsSpectraFullCorrMat(    histRpPbIndStatEta[cent],    graphRpPbIndSystEta[cent],
                                                                                                xPtLimitsEta[cent], maxNBinsEta[cent],
                                                                                                offSetsEta[cent], offSetsEtaSys[cent],
                                                                                                graphRpPbCombStatEta[cent], graphRpPbCombSystEta[cent],
                                                                                                fileNameEtaRpPbOutputWeighting, "pPb_5.023TeV", "EtaRpPb", kTRUE,
                                                                                                NULL, fileNameCorrFactors, centArrayCorr[cent]
                                                                                            );
            if (graphRpPbCombCombEta == NULL) {
                cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
                return;
            }
            while (graphRpPbCombStatEta[cent]->GetX()[0] < minPtEta[cent]){
                graphRpPbCombStatEta[cent]->RemovePoint(0);
            }
            while (graphRpPbCombSystEta[cent]->GetX()[0] < minPtEta[cent]){
                graphRpPbCombSystEta[cent]->RemovePoint(0);
            }
            while (graphRpPbCombCombEta[cent]->GetX()[0] < minPtEta[cent]){
                graphRpPbCombCombEta[cent]->RemovePoint(0);
            }
            graphRpPbCombCombEta[cent]->Print();

            graphRpPbCombStatEtaWOXErr[cent]            = (TGraphAsymmErrors*)graphRpPbCombStatEta[cent]->Clone("graphRpPbCombStatEtaWOXErr");
            ProduceGraphAsymmWithoutXErrors(graphRpPbCombStatEtaWOXErr[cent]);


            // Reading weights from output file for plotting
            ifstream fileWeightsEtaRpPbRead;
            fileWeightsEtaRpPbRead.open(fileNameEtaRpPbOutputWeighting,ios_base::in);
            cout << "reading" << fileNameEtaRpPbOutputWeighting << endl;
            Double_t xValuesEtaRpPbRead[50];
            Double_t weightsEtaRpPbRead[11][50];
            Int_t availableEtaRpPbMeas[11]      = { -1, -1, -1, -1, -1,
                                                -1, -1, -1, -1, -1,
                                                -1};
            Int_t nMeasSetEtaRpPb               = nTotMeasEta[cent];
            Int_t nPtBinsEtaRpPbRead            = 0;
            while(!fileWeightsEtaRpPbRead.eof() && nPtBinsEtaRpPbRead < 50){
                TString garbage             = "";
                if (nPtBinsEtaRpPbRead == 0){
                    fileWeightsEtaRpPbRead >> garbage ;//>> availableEtaRpPbMeas[0] >> availableEtaRpPbMeas[1] >> availableEtaRpPbMeas[2] >> availableEtaRpPbMeas[3];
                    for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
                        fileWeightsEtaRpPbRead >> availableEtaRpPbMeas[i] ;
                    }
                    cout << "read following measurements: ";
                    for (Int_t meth= 0; meth < 11; meth++){
                        cout << availableEtaRpPbMeas[meth] << "\t" ;
                    }
                    cout << endl;
                } else {
                    fileWeightsEtaRpPbRead >> xValuesEtaRpPbRead[nPtBinsEtaRpPbRead-1];
                    for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
                        fileWeightsEtaRpPbRead >> weightsEtaRpPbRead[availableEtaRpPbMeas[i]][nPtBinsEtaRpPbRead-1] ;
                    }
                    cout << "read: "<<  nPtBinsEtaRpPbRead << "\t"<< xValuesEtaRpPbRead[nPtBinsEtaRpPbRead-1] << "\t" ;
                    for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
                        cout << weightsEtaRpPbRead[availableEtaRpPbMeas[i]][nPtBinsEtaRpPbRead-1] << "\t";
                    }
                    cout << endl;
                }
                nPtBinsEtaRpPbRead++;
            }
            nPtBinsEtaRpPbRead                  = nPtBinsEtaRpPbRead-2 ;
            fileWeightsEtaRpPbRead.close();

            for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
                graphWeightsEtaRpPb[cent][availableEtaRpPbMeas[i]]      = new TGraph(nPtBinsEtaRpPbRead,xValuesEtaRpPbRead,weightsEtaRpPbRead[availableEtaRpPbMeas[i]]);
                Int_t bin = 0;
                for (Int_t n = 0; n< nPtBinsEtaRpPbRead; n++){
                    if (graphWeightsEtaRpPb[cent][availableEtaRpPbMeas[i]]->GetY()[bin] == 0) graphWeightsEtaRpPb[cent][availableEtaRpPbMeas[i]]->RemovePoint(bin);
                    else bin++;
                }
            }

            // **********************************************************************************************************************
            // ******************************************* Plotting weights method only EMC *****************************************
            // **********************************************************************************************************************
            textSizeLabelsPixel           = 900*0.04;
            canvasWeights->cd();

            histo2DEtaWeights->Draw("copy");
                TLegend* legendWeightsEtaRpPb   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetEtaRpPb+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
                for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
                    DrawGammaSetMarkerTGraph(graphWeightsEtaRpPb[cent][availableEtaRpPbMeas[i]], markerStyleDet[availableEtaRpPbMeas[i]], markerSizeDet[availableEtaRpPbMeas[i]]*0.5, colorDet[availableEtaRpPbMeas[i]] , colorDet[availableEtaRpPbMeas[i]]);
                    graphWeightsEtaRpPb[cent][availableEtaRpPbMeas[i]]->Draw("p,same,z");
                    legendWeightsEtaRpPb->AddEntry(graphWeightsEtaRpPb[cent][availableEtaRpPbMeas[i]],nameMeasGlobalLabel[availableEtaRpPbMeas[i]],"p");
                }
                legendWeightsEtaRpPb->Draw();

                labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
                labelWeightsEnergy->Draw();
                labelWeightsEtaRpPb->Draw();

                DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
                DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
                DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
                DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);

            canvasWeights->SaveAs(Form("%s/EtaRpPb_Weights_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelSysErr->cd();

                histo2DRelSysErr->GetYaxis()->SetRangeUser(0,49.5);
                histo2DRelSysErr->Draw("copy");
                TLegend* legendRelSysErrRpPb       = GetAndSetLegend2(0.60, 0.92-(0.04*(nMeasSetEtaRpPb)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
                for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
                    if (!sysErrorRelCollectionEtaRpA[cent][availableEtaRpPbMeas[i]]) continue;
                    DrawGammaSetMarkerTGraph(sysErrorRelCollectionEtaRpA[cent][availableEtaRpPbMeas[i]], markerStyleDet[availableEtaRpPbMeas[i]], markerSizeDet[availableEtaRpPbMeas[i]]*0.5,
                                            colorDet[availableEtaRpPbMeas[i]], colorDet[availableEtaRpPbMeas[i]]);
                    sysErrorRelCollectionEtaRpA[cent][availableEtaRpPbMeas[i]]->Draw("p,same,z");
                    legendRelSysErrRpPb->AddEntry(sysErrorRelCollectionEtaRpA[cent][availableEtaRpPbMeas[i]],nameMeasGlobalLabel[availableEtaRpPbMeas[i]],"p");
                }
                legendRelSysErrRpPb->Draw();

                labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
                labelRelSysErrEnergy->Draw();
                labelRelSysErrEtaRpPb->Draw();

            canvasRelSysErr->SaveAs(Form("%s/EtaRpPb_RelSysErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelSysErrRpPb;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelStatErr->cd();

                histo2DRelStatErr->GetYaxis()->SetRangeUser(0,49.5);
                histo2DRelStatErr->Draw("copy");
                TLegend* legendRelStatErrRpPb       = GetAndSetLegend2(0.12, 0.92-(0.04*(nMeasSetEtaRpPb)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
                for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
                    if (!statErrorRelCollectionEtaRpA[cent][availableEtaRpPbMeas[i]]) continue;
                    DrawGammaSetMarkerTGraph(statErrorRelCollectionEtaRpA[cent][availableEtaRpPbMeas[i]], markerStyleDet[availableEtaRpPbMeas[i]], markerSizeDet[availableEtaRpPbMeas[i]]*0.5,
                                            colorDet[availableEtaRpPbMeas[i]], colorDet[availableEtaRpPbMeas[i]]);
                    statErrorRelCollectionEtaRpA[cent][availableEtaRpPbMeas[i]]->Draw("p,same,z");
                    legendRelStatErrRpPb->AddEntry(statErrorRelCollectionEtaRpA[cent][availableEtaRpPbMeas[i]],nameMeasGlobalLabel[availableEtaRpPbMeas[i]],"p");
                }
                legendRelStatErrRpPb->Draw();

                labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
                labelRelStatErrEnergy->Draw();
                labelRelStatErrEtaRpPb->Draw();

            canvasRelStatErr->SaveAs(Form("%s/EtaRpPb_RelStatErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelStatErrRpPb;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            graphCombEtaRpPbRelStat[cent]      = CalculateRelErrUpAsymmGraph( graphRpPbCombStatEta[cent], "relativeStatErrorEtaRpPb_Method");
            graphCombEtaRpPbRelSys[cent]       = CalculateRelErrUpAsymmGraph( graphRpPbCombSystEta[cent], "relativeSysErrorEtaRpPb_Method");
            graphCombEtaRpPbRelTot[cent]       = CalculateRelErrUpAsymmGraph( graphRpPbCombCombEta[cent], "relativeTotalErrorEtaRpPb_Method");

            canvasRelTotErr->cd();

                histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,49.5);
                histo2DRelTotErrEta->Draw("copy");
                DrawGammaSetMarkerTGraphAsym(graphCombEtaRpPbRelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
                graphCombEtaRpPbRelTot[cent]->Draw("p,same,z");
                DrawGammaSetMarkerTGraphAsym(graphCombEtaRpPbRelStat[cent], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
                graphCombEtaRpPbRelStat[cent]->Draw("l,x0,same,e1");
                DrawGammaSetMarkerTGraphAsym(graphCombEtaRpPbRelSys[cent], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
                graphCombEtaRpPbRelSys[cent]->SetLineStyle(7);
                graphCombEtaRpPbRelSys[cent]->Draw("l,x0,same,e1");

                if (centArray[cent].BeginsWith("0-100%")){
                    DrawGammaSetMarkerTGraphAsym(graphCombEtaRpPbOlderMBRelErrTot, markerStyleComb+4, markerSizeComb, kBlack , kBlack);
                    graphCombEtaRpPbOlderMBRelErrTot->Draw("p,same,z");
                    DrawGammaSetMarkerTGraphAsym(graphCombEtaRpPbOlderMBRelErrStat, markerStyleComb+4, markerSizeComb, kGray+2 , kGray+2);
                    graphCombEtaRpPbOlderMBRelErrStat->Draw("l,x0,same,e1");
                    DrawGammaSetMarkerTGraphAsym(graphCombEtaRpPbOlderMBRelErrSys, markerStyleComb+4, markerSizeComb, kGray+1, kGray+1);
                    graphCombEtaRpPbOlderMBRelErrSys->SetLineStyle(7);
                    graphCombEtaRpPbOlderMBRelErrSys->Draw("l,x0,same,e1");

                    TLegend* legendRelTotErr4       = GetAndSetLegend2(0.35, 0.92-(0.035*3), 0.65, 0.92, 32);
                    legendRelTotErr4->AddEntry(graphCombEtaRpPbOlderMBRelErrTot,"pub tot","p");
                    legendRelTotErr4->AddEntry(graphCombEtaRpPbOlderMBRelErrStat,"pub stat","l");
                    legendRelTotErr4->AddEntry(graphCombEtaRpPbOlderMBRelErrSys,"pub sys","l");
                    legendRelTotErr4->Draw();
                }

                TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.35, 0.92, 32);
                legendRelTotErr3->AddEntry(graphCombEtaRpPbRelTot[cent],"tot","p");
                legendRelTotErr3->AddEntry(graphCombEtaRpPbRelStat[cent],"stat","l");
                legendRelTotErr3->AddEntry(graphCombEtaRpPbRelSys[cent],"sys","l");
                legendRelTotErr3->Draw();

                labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
                labelRelTotErrEnergy->Draw();
                labelRelTotErrEtaRpPb->Draw();

            canvasRelTotErr->SaveAs(Form("%s/EtaRpPb_Reldecomp_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelTotErr3;
        }
    }

    TGraphAsymmErrors* graphRCPIndStatPi0[17][11];
    TGraphAsymmErrors* graphRCPIndStatPi0WOXErr[17][11];
    TH1D* histRCPIndStatPi0[17][11];
    TGraphAsymmErrors* graphRCPIndSystPi0[17][11];
    TGraphAsymmErrors* graphRCPIndCombPi0[17][11];
    TGraph* graphWeightsPi0RCP[17][11];
    TGraphAsymmErrors* statErrorRelCollectionPi0RCP[17][11];
    TGraphAsymmErrors* sysErrorRelCollectionPi0RCP[17][11];
    TGraphAsymmErrors* graphRCPCombStatPi0[17];
    TGraphAsymmErrors* graphRCPCombSystPi0[17];
    TGraphAsymmErrors* graphRCPCombStatPi0WOXErr[17];
    TGraphAsymmErrors* graphRCPCombCombPi0[17];
    TGraphAsymmErrors* graphCombPi0RCPRelStat[17];
    TGraphAsymmErrors* graphCombPi0RCPRelSys[17];
    TGraphAsymmErrors* graphCombPi0RCPRelTot[17];


    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if ( !enableCentRCP[cent] ) continue;
        // *******************************************************************************************************
        // *************************** RCP calculation for pi0 for individual measurements **********************
        // *******************************************************************************************************
        rCPNColl[cent]      = nCollpPb[cent]/(nCollpPb[rCPRefInt[cent]]);
        rCPNCollErr[cent]   = rCPNColl[cent]* TMath::Sqrt(TMath::Power(nCollErrpPb[cent]/nCollpPb[cent],2) + TMath::Power(nCollErrpPb[rCPRefInt[cent]]/nCollpPb[rCPRefInt[cent]],2));
        clog << nCollpPb[cent] << "\t" << nCollpPb[rCPRefInt[cent]] << endl;
        clog << rCPNColl[cent] << "\t" << rCPNCollErr[cent] << endl;
        for (Int_t meth = 0; meth < 11; meth++){
            cout << "RCP calc for pi0 " << centArray[cent].Data() << addCentString[cent].Data() << "\t" << runArray[cent].Data() << "\t" << nameMeasGlobal[meth].Data() << endl;
            graphWeightsPi0RCP[cent][meth]                              = NULL;
            graphRCPIndStatPi0[cent][meth]                              = NULL;
            graphRCPIndStatPi0WOXErr[cent][meth]                        = NULL;
            histRCPIndStatPi0[cent][meth]                               = NULL;
            graphRCPIndSystPi0[cent][meth]                              = NULL;
            graphRCPIndCombPi0[cent][meth]                              = NULL;
            statErrorRelCollectionPi0RCP[cent][meth]                     = NULL;
            sysErrorRelCollectionPi0RCP[cent][meth]                      = NULL;
            if (havePi0SysDetailedperi[cent][meth] && graphIndPi0InvYieldStat_yShifted[rCPRefInt[cent]][meth] && graphIndPi0InvYieldStat_yShifted[cent][meth] ){ //&& !( i == 3 )
                if(bWCorrection.Contains("Y")  ){
                    TGraphAsymmErrors* tempGraphPeriStat                = (TGraphAsymmErrors*)graphIndPi0InvYieldStat_yShifted[rCPRefInt[cent]][meth]->Clone("tempPeriStat");
                    TGraphAsymmErrors* tempGraphPeriSys                 = (TGraphAsymmErrors*)graphIndPi0InvYieldSys_yShifted[rCPRefInt[cent]][meth]->Clone("tempPeriSys");
                    tempGraphPeriStat                                   = ScaleGraph(tempGraphPeriStat, 1./tpPb[rCPRefInt[cent]]);
                    tempGraphPeriSys                                    = ScaleGraph(tempGraphPeriSys, 1./tpPb[rCPRefInt[cent]]);
                    RecalculateErrorsBasedOnDetailedInputFile( tempGraphPeriSys, fileNamesperiPi0DetailedSys[cent][meth] );
                    TGraphAsymmErrors* tempGraphCentStat                = (TGraphAsymmErrors*)graphIndPi0InvYieldStat_yShifted[cent][meth]->Clone("tempCentStat");
                    TGraphAsymmErrors* tempGraphCentSys                 = (TGraphAsymmErrors*)graphIndPi0InvYieldSys_yShifted[cent][meth]->Clone("tempCentSys");
                    tempGraphCentStat                                   = ScaleGraph(tempGraphCentStat, 1./tpPb[cent]);
                    tempGraphCentSys                                    = ScaleGraph(tempGraphCentSys, 1./tpPb[cent]);
                    if (tempGraphPeriStat->GetX()[0] != tempGraphCentStat->GetX()[0])
                        continue;

                    graphRCPIndCombPi0[cent][meth]                      = CalcRpPbV2(   tempGraphPeriStat,
                                                                                        tempGraphPeriSys,
                                                                                        tempGraphCentStat,
                                                                                        tempGraphCentSys,
                                                                                        &graphRCPIndStatPi0[cent][meth],
                                                                                        &graphRCPIndSystPi0[cent][meth],
                                                                                        1, 0 ,
                                                                                        ptSysRemRCPNames[cent][meth], fileNamespPbPi0DetailedSys[cent][meth], fileNamesperiPi0DetailedSys[cent][meth], fileNamesRCPPi0DetailedSys[cent][meth]
                                                                                    );
                    graphRCPIndCombPi0[cent][meth]->Print();
                    histRCPIndStatPi0[cent][meth]                       = (TH1D*)statErrorCollectionPi0[cent][meth]->Clone(Form("histRCP%sStatPi0",nameMeasGlobalLabel[meth].Data()));
                    for (Int_t j = 0; j< histRCPIndStatPi0[cent][meth]->GetNbinsX()+1; j++){
                        histRCPIndStatPi0[cent][meth]->SetBinContent(j,0);
                        histRCPIndStatPi0[cent][meth]->SetBinError(j,0);
                    }
                    for (Int_t j = 0; j < graphRCPIndStatPi0[cent][meth]->GetN(); j++){
                        Int_t bin                               = histRCPIndStatPi0[cent][meth]->GetXaxis()->FindBin(graphRCPIndStatPi0[cent][meth]->GetX()[j]);
                        histRCPIndStatPi0[cent][meth]->SetBinContent(bin,graphRCPIndStatPi0[cent][meth]->GetY()[j]);
                        histRCPIndStatPi0[cent][meth]->SetBinError(bin,graphRCPIndStatPi0[cent][meth]->GetEYlow()[j]);
                    }
                }

                if (histRCPIndStatPi0[cent][meth]){
                    statErrorRelCollectionPi0RCP[cent][meth]   = (TGraphAsymmErrors*)graphRCPIndStatPi0[cent][meth]->Clone(Form("relativeStatErrorPi0RCP%i_%s", cent, nameMeasGlobal[meth].Data()));
                    RemoveZerosAtBeginningAndEndFromGraph(statErrorRelCollectionPi0RCP[cent][meth]);
                    statErrorRelCollectionPi0RCP[cent][meth]   = CalculateRelErrUpAsymmGraph( statErrorRelCollectionPi0RCP[cent][meth], Form("relativeStatErrorPi0RCP_%i_%s", cent,nameMeasGlobal[meth].Data()));
                    statErrorRelCollectionPi0RCP[cent][meth]->Print();
                    graphRCPIndStatPi0WOXErr[cent][meth]      = (TGraphAsymmErrors*)graphRCPIndStatPi0[cent][meth]->Clone(Form("graphRCP%sStatPi0WOXErr%i", nameMeasGlobalLabel[meth].Data(), cent));

                }
                if (graphRCPIndSystPi0[cent][meth]){
                    sysErrorRelCollectionPi0RCP[cent][meth]    = (TGraphAsymmErrors*)graphRCPIndSystPi0[cent][meth]->Clone(Form("relativeSysErrorPi0RCP%i_%s", cent, nameMeasGlobal[meth].Data()));
                    RemoveZerosAtBeginningAndEndFromGraph(sysErrorRelCollectionPi0RCP[cent][meth]);
                    sysErrorRelCollectionPi0RCP[cent][meth]    = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionPi0RCP[cent][meth], Form("relativeSysErrorPi0RCP%i_%s", cent, nameMeasGlobal[meth].Data()));
                }
            }
        }
        // *******************************************************************************************************
        // *************************** RCP calculation for pi0 **************************************************
        // *******************************************************************************************************
        graphRCPCombStatPi0[cent]                      = NULL;
        graphRCPCombStatPi0[cent]                      = NULL;
        graphRCPCombStatPi0[cent]                      = NULL;

        if (!enableCentComb[cent]) continue;
        if(bWCorrection.Contains("Y") ){
            TString fileNamePi0RCPOutputWeighting      = Form("%s/Pi0RCP_WeightingMethod_%s%s%s.dat", outputDirSupportComb.Data(), centArrayOutput[cent].Data(),
                                                               addCentString[cent].Data(), runArray[cent].Data());
            graphRCPCombCombPi0[cent]                  = CombinePtPointsSpectraFullCorrMat(    histRCPIndStatPi0[cent],    graphRCPIndSystPi0[cent],
                                                                                                xPtLimitsPi0[cent], maxNBinsPi0[cent],
                                                                                                offSetsPi0[cent], offSetsPi0Sys[cent],
                                                                                                graphRCPCombStatPi0[cent], graphRCPCombSystPi0[cent],
                                                                                                fileNamePi0RCPOutputWeighting, "pPb_5.023TeV", "Pi0RCP", kTRUE,
                                                                                                NULL, fileNameCorrFactors, centArrayCorr[cent]
            );
            if (graphRCPCombCombPi0 == NULL) {
                cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
                return;
            }
            while (graphRCPCombStatPi0[cent]->GetX()[0] < minPtPi0[cent]){
                graphRCPCombStatPi0[cent]->RemovePoint(0);
            }
            while (graphRCPCombSystPi0[cent]->GetX()[0] < minPtPi0[cent]){
                graphRCPCombSystPi0[cent]->RemovePoint(0);
            }
            while (graphRCPCombCombPi0[cent]->GetX()[0] < minPtPi0[cent]){
                graphRCPCombCombPi0[cent]->RemovePoint(0);
            }
            graphRCPCombCombPi0[cent]->Print();

            graphRCPCombStatPi0WOXErr[cent]            = (TGraphAsymmErrors*)graphRCPCombStatPi0[cent]->Clone("graphRCPCombStatPi0WOXErr");
            ProduceGraphAsymmWithoutXErrors(graphRCPCombStatPi0WOXErr[cent]);


            // Reading weights from output file for plotting
            ifstream fileWeightsPi0RCPRead;
            fileWeightsPi0RCPRead.open(fileNamePi0RCPOutputWeighting,ios_base::in);
            cout << "reading" << fileNamePi0RCPOutputWeighting << endl;
            Double_t xValuesPi0RCPRead[50];
            Double_t weightsPi0RCPRead[11][50];
            Int_t availablePi0RCPMeas[11]      = { -1, -1, -1, -1, -1,
                                                    -1, -1, -1, -1, -1,
                                                    -1};
            Int_t nMeasSetPi0RCP               = nTotMeasPi0[cent];
            Int_t nPtBinsPi0RCPRead            = 0;
            while(!fileWeightsPi0RCPRead.eof() && nPtBinsPi0RCPRead < 50){
                TString garbage             = "";
                if (nPtBinsPi0RCPRead == 0){
                    fileWeightsPi0RCPRead >> garbage ;//>> availablePi0RCPMeas[0] >> availablePi0RCPMeas[1] >> availablePi0RCPMeas[2] >> availablePi0RCPMeas[3];
                    for (Int_t i = 0; i < nMeasSetPi0RCP; i++){
                        fileWeightsPi0RCPRead >> availablePi0RCPMeas[i] ;
                    }
                    cout << "read following measurements: ";
                    for (Int_t meth= 0; meth < 11; meth++){
                        cout << availablePi0RCPMeas[meth] << "\t" ;
                    }
                    cout << endl;
                } else {
                    fileWeightsPi0RCPRead >> xValuesPi0RCPRead[nPtBinsPi0RCPRead-1];
                    for (Int_t i = 0; i < nMeasSetPi0RCP; i++){
                        fileWeightsPi0RCPRead >> weightsPi0RCPRead[availablePi0RCPMeas[i]][nPtBinsPi0RCPRead-1] ;
                    }
                    cout << "read: "<<  nPtBinsPi0RCPRead << "\t"<< xValuesPi0RCPRead[nPtBinsPi0RCPRead-1] << "\t" ;
                    for (Int_t i = 0; i < nMeasSetPi0RCP; i++){
                        cout << weightsPi0RCPRead[availablePi0RCPMeas[i]][nPtBinsPi0RCPRead-1] << "\t";
                    }
                    cout << endl;
                }
                nPtBinsPi0RCPRead++;
            }
            nPtBinsPi0RCPRead                  = nPtBinsPi0RCPRead-2 ;
            fileWeightsPi0RCPRead.close();

            for (Int_t i = 0; i < nMeasSetPi0RCP; i++){
                graphWeightsPi0RCP[cent][availablePi0RCPMeas[i]]      = new TGraph(nPtBinsPi0RCPRead,xValuesPi0RCPRead,weightsPi0RCPRead[availablePi0RCPMeas[i]]);
                Int_t bin = 0;
                for (Int_t n = 0; n< nPtBinsPi0RCPRead; n++){
                    if (graphWeightsPi0RCP[cent][availablePi0RCPMeas[i]]->GetY()[bin] == 0) graphWeightsPi0RCP[cent][availablePi0RCPMeas[i]]->RemovePoint(bin);
                    else bin++;
                }
            }

            // **********************************************************************************************************************
            // ******************************************* Plotting weights method only EMC *****************************************
            // **********************************************************************************************************************
            textSizeLabelsPixel           = 900*0.04;
            canvasWeights->cd();

            histo2DPi0Weights->Draw("copy");
            TLegend* legendWeightsPi0RCP   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetPi0RCP+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
            for (Int_t i = 0; i < nMeasSetPi0RCP; i++){
                DrawGammaSetMarkerTGraph(graphWeightsPi0RCP[cent][availablePi0RCPMeas[i]], markerStyleDet[availablePi0RCPMeas[i]], markerSizeDet[availablePi0RCPMeas[i]]*0.5, colorDet[availablePi0RCPMeas[i]] , colorDet[availablePi0RCPMeas[i]]);
                graphWeightsPi0RCP[cent][availablePi0RCPMeas[i]]->Draw("p,same,z");
                legendWeightsPi0RCP->AddEntry(graphWeightsPi0RCP[cent][availablePi0RCPMeas[i]],nameMeasGlobalLabel[availablePi0RCPMeas[i]],"p");
            }
            legendWeightsPi0RCP->Draw();

            labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelWeightsEnergy->Draw();
            labelWeightsPi0RCP->Draw();

            DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);

            canvasWeights->SaveAs(Form("%s/Pi0RCP_Weights_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelSysErr->cd();

            histo2DRelSysErr->GetYaxis()->SetRangeUser(0,49.5);
            histo2DRelSysErr->Draw("copy");
            TLegend* legendRelSysErrRCP       = GetAndSetLegend2(0.60, 0.92-(0.04*(nMeasSetPi0RCP)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
            for (Int_t i = 0; i < nMeasSetPi0RCP; i++){
                if (!sysErrorRelCollectionPi0RCP[cent][availablePi0RCPMeas[i]]) continue;
                DrawGammaSetMarkerTGraph(sysErrorRelCollectionPi0RCP[cent][availablePi0RCPMeas[i]], markerStyleDet[availablePi0RCPMeas[i]], markerSizeDet[availablePi0RCPMeas[i]]*0.5,
                                        colorDet[availablePi0RCPMeas[i]], colorDet[availablePi0RCPMeas[i]]);
                sysErrorRelCollectionPi0RCP[cent][availablePi0RCPMeas[i]]->Draw("p,same,z");
                legendRelSysErrRCP->AddEntry(sysErrorRelCollectionPi0RCP[cent][availablePi0RCPMeas[i]],nameMeasGlobalLabel[availablePi0RCPMeas[i]],"p");
            }
            legendRelSysErrRCP->Draw();

            labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelSysErrEnergy->Draw();
            labelRelSysErrPi0RpPb->Draw();

            canvasRelSysErr->SaveAs(Form("%s/Pi0RCP_RelSysErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelSysErrRCP;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelStatErr->cd();

            histo2DRelStatErr->GetYaxis()->SetRangeUser(0,49.5);
            histo2DRelStatErr->Draw("copy");
            TLegend* legendRelStatErrRCP       = GetAndSetLegend2(0.12, 0.92-(0.04*(nMeasSetPi0RCP)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
            for (Int_t i = 0; i < nMeasSetPi0RCP; i++){
                if (!statErrorRelCollectionPi0RCP[cent][availablePi0RCPMeas[i]]) continue;
                DrawGammaSetMarkerTGraph(statErrorRelCollectionPi0RCP[cent][availablePi0RCPMeas[i]], markerStyleDet[availablePi0RCPMeas[i]], markerSizeDet[availablePi0RCPMeas[i]]*0.5,
                                        colorDet[availablePi0RCPMeas[i]], colorDet[availablePi0RCPMeas[i]]);
                statErrorRelCollectionPi0RCP[cent][availablePi0RCPMeas[i]]->Draw("p,same,z");
                legendRelStatErrRCP->AddEntry(statErrorRelCollectionPi0RCP[cent][availablePi0RCPMeas[i]],nameMeasGlobalLabel[availablePi0RCPMeas[i]],"p");
            }
            legendRelStatErrRCP->Draw();

            labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelStatErrEnergy->Draw();
            labelRelStatErrPi0RCP->Draw();

            canvasRelStatErr->SaveAs(Form("%s/Pi0RCP_RelStatErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelStatErrRCP;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            graphCombPi0RCPRelStat[cent]      = CalculateRelErrUpAsymmGraph( graphRCPCombStatPi0[cent], "relativeStatErrorPi0RCP_Method");
            graphCombPi0RCPRelSys[cent]       = CalculateRelErrUpAsymmGraph( graphRCPCombSystPi0[cent], "relativeSysErrorPi0RCP_Method");
            graphCombPi0RCPRelTot[cent]       = CalculateRelErrUpAsymmGraph( graphRCPCombCombPi0[cent], "relativeTotalErrorPi0RCP_Method");

            canvasRelTotErr->cd();

            histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,49.5);
            histo2DRelTotErrPi0->Draw("copy");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0RCPRelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphCombPi0RCPRelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0RCPRelStat[cent], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
            graphCombPi0RCPRelStat[cent]->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0RCPRelSys[cent], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
            graphCombPi0RCPRelSys[cent]->SetLineStyle(7);
            graphCombPi0RCPRelSys[cent]->Draw("l,x0,same,e1");

            TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
            legendRelTotErr3->AddEntry(graphCombPi0RCPRelTot[cent],"tot","p");
            legendRelTotErr3->AddEntry(graphCombPi0RCPRelStat[cent],"stat","l");
            legendRelTotErr3->AddEntry(graphCombPi0RCPRelSys[cent],"sys","l");
            legendRelTotErr3->Draw();

            labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelTotErrEnergy->Draw();
            labelRelTotErrPi0RCP->Draw();

            canvasRelTotErr->SaveAs(Form("%s/Pi0RCP_Reldecomp_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelTotErr3;
        }
    }

    TGraphAsymmErrors* graphRCPIndStatEta[17][11];
    TGraphAsymmErrors* graphRCPIndStatEtaWOXErr[17][11];
    TH1D* histRCPIndStatEta[17][11];
    TGraphAsymmErrors* graphRCPIndSystEta[17][11];
    TGraphAsymmErrors* graphRCPIndCombEta[17][11];
    TGraph* graphWeightsEtaRCP[17][11];
    TGraphAsymmErrors* statErrorRelCollectionEtaRCP[17][11];
    TGraphAsymmErrors* sysErrorRelCollectionEtaRCP[17][11];
    TGraphAsymmErrors* graphRCPCombStatEta[17];
    TGraphAsymmErrors* graphRCPCombSystEta[17];
    TGraphAsymmErrors* graphRCPCombStatEtaWOXErr[17];
    TGraphAsymmErrors* graphRCPCombCombEta[17];
    TGraphAsymmErrors* graphCombEtaRCPRelStat[17];
    TGraphAsymmErrors* graphCombEtaRCPRelSys[17];
    TGraphAsymmErrors* graphCombEtaRCPRelTot[17];


    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if ( !enableCentRCP[cent] ) continue;
        // *******************************************************************************************************
        // *************************** RCP calculation for pi0 for individual measurements **********************
        // *******************************************************************************************************
        for (Int_t meth = 0; meth < 11; meth++){
            cout << "RCP calc for pi0 " << centArray[cent].Data() << addCentString[cent].Data() << "\t" << runArray[cent].Data() << "\t" << nameMeasGlobal[meth].Data() << endl;
            graphWeightsEtaRCP[cent][meth]                              = NULL;
            graphRCPIndStatEta[cent][meth]                              = NULL;
            graphRCPIndStatEtaWOXErr[cent][meth]                        = NULL;
            histRCPIndStatEta[cent][meth]                               = NULL;
            graphRCPIndSystEta[cent][meth]                              = NULL;
            graphRCPIndCombEta[cent][meth]                              = NULL;
            statErrorRelCollectionEtaRCP[cent][meth]                     = NULL;
            sysErrorRelCollectionEtaRCP[cent][meth]                      = NULL;
            if (haveEtaSysDetailedperi[cent][meth] && graphIndEtaInvYieldStat_yShifted[rCPRefInt[cent]][meth] && graphIndEtaInvYieldStat_yShifted[cent][meth] ){ //&& !( i == 3 )
                if(bWCorrection.Contains("Y")  ){
                    TGraphAsymmErrors* tempGraphPeriStat                = (TGraphAsymmErrors*)graphIndEtaInvYieldStat_yShifted[rCPRefInt[cent]][meth]->Clone("tempPeriStat");
                    TGraphAsymmErrors* tempGraphPeriSys                 = (TGraphAsymmErrors*)graphIndEtaInvYieldSys_yShifted[rCPRefInt[cent]][meth]->Clone("tempPeriSys");
                    tempGraphPeriStat                                   = ScaleGraph(tempGraphPeriStat, 1./tpPb[rCPRefInt[cent]]);
                    tempGraphPeriSys                                    = ScaleGraph(tempGraphPeriSys, 1./tpPb[rCPRefInt[cent]]);
                    RecalculateErrorsBasedOnDetailedInputFile( tempGraphPeriSys, fileNamesperiEtaDetailedSys[cent][meth] );
                    TGraphAsymmErrors* tempGraphCentStat                = (TGraphAsymmErrors*)graphIndEtaInvYieldStat_yShifted[cent][meth]->Clone("tempCentStat");
                    TGraphAsymmErrors* tempGraphCentSys                 = (TGraphAsymmErrors*)graphIndEtaInvYieldSys_yShifted[cent][meth]->Clone("tempCentSys");
                    tempGraphCentStat                                   = ScaleGraph(tempGraphCentStat, 1./tpPb[cent]);
                    tempGraphCentSys                                    = ScaleGraph(tempGraphCentSys, 1./tpPb[cent]);

                    if (tempGraphPeriStat->GetX()[0] != tempGraphCentStat->GetX()[0])
                        continue;

                    graphRCPIndCombEta[cent][meth]                      = CalcRpPbV2(   tempGraphPeriStat,
                                                                                        tempGraphPeriSys,
                                                                                        tempGraphCentStat,
                                                                                        tempGraphCentSys,
                                                                                        &graphRCPIndStatEta[cent][meth],
                                                                                        &graphRCPIndSystEta[cent][meth],
                                                                                        1, 0 ,
                                                                                        ptSysRemRCPNames[cent][meth], fileNamespPbEtaDetailedSys[cent][meth], fileNamesperiEtaDetailedSys[cent][meth], fileNamesRCPEtaDetailedSys[cent][meth]
                    );
                    graphRCPIndCombEta[cent][meth]->Print();
                    histRCPIndStatEta[cent][meth]                       = (TH1D*)statErrorCollectionEta[cent][meth]->Clone(Form("histRCP%sStatEta",nameMeasGlobalLabel[meth].Data()));
                    for (Int_t j = 0; j< histRCPIndStatEta[cent][meth]->GetNbinsX()+1; j++){
                        histRCPIndStatEta[cent][meth]->SetBinContent(j,0);
                        histRCPIndStatEta[cent][meth]->SetBinError(j,0);
                    }
                    for (Int_t j = 0; j < graphRCPIndStatEta[cent][meth]->GetN(); j++){
                        Int_t bin                               = histRCPIndStatEta[cent][meth]->GetXaxis()->FindBin(graphRCPIndStatEta[cent][meth]->GetX()[j]);
                        histRCPIndStatEta[cent][meth]->SetBinContent(bin,graphRCPIndStatEta[cent][meth]->GetY()[j]);
                        histRCPIndStatEta[cent][meth]->SetBinError(bin,graphRCPIndStatEta[cent][meth]->GetEYlow()[j]);
                    }
                }

                if (histRCPIndStatEta[cent][meth]){
                    statErrorRelCollectionEtaRCP[cent][meth]   = (TGraphAsymmErrors*)graphRCPIndStatEta[cent][meth]->Clone(Form("relativeStatErrorEtaRCP%i_%s", cent, nameMeasGlobal[meth].Data()));
                    RemoveZerosAtBeginningAndEndFromGraph(statErrorRelCollectionEtaRCP[cent][meth]);
                    statErrorRelCollectionEtaRCP[cent][meth]   = CalculateRelErrUpAsymmGraph( statErrorRelCollectionEtaRCP[cent][meth], Form("relativeStatErrorEtaRCP_%i_%s", cent,nameMeasGlobal[meth].Data()));
                    statErrorRelCollectionEtaRCP[cent][meth]->Print();
                    graphRCPIndStatEtaWOXErr[cent][meth]      = (TGraphAsymmErrors*)graphRCPIndStatEta[cent][meth]->Clone(Form("graphRCP%sStatEtaWOXErr%i", nameMeasGlobalLabel[meth].Data(), cent));

                }
                if (graphRCPIndSystEta[cent][meth]){
                    sysErrorRelCollectionEtaRCP[cent][meth]    = (TGraphAsymmErrors*)graphRCPIndSystEta[cent][meth]->Clone(Form("relativeSysErrorEtaRCP%i_%s", cent, nameMeasGlobal[meth].Data()));
                    RemoveZerosAtBeginningAndEndFromGraph(sysErrorRelCollectionEtaRCP[cent][meth]);
                    sysErrorRelCollectionEtaRCP[cent][meth]    = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionEtaRCP[cent][meth], Form("relativeSysErrorEtaRCP%i_%s", cent, nameMeasGlobal[meth].Data()));
                }
            }
        }
        // *******************************************************************************************************
        // *************************** RCP calculation for eta **************************************************
        // *******************************************************************************************************
        graphRCPCombStatEta[cent]                      = NULL;
        graphRCPCombStatEta[cent]                      = NULL;
        graphRCPCombStatEta[cent]                      = NULL;

        if (!enableCentComb[cent]) continue;
        if(bWCorrection.Contains("Y") ){
            TString fileNameEtaRCPOutputWeighting      = Form("%s/EtaRCP_WeightingMethod_%s%s%s.dat", outputDirSupportComb.Data(), centArrayOutput[cent].Data(),
                                                               addCentString[cent].Data(), runArray[cent].Data());
            graphRCPCombCombEta[cent]                  = CombinePtPointsSpectraFullCorrMat(    histRCPIndStatEta[cent],    graphRCPIndSystEta[cent],
                                                                                                xPtLimitsEta[cent], maxNBinsEta[cent],
                                                                                                offSetsEta[cent], offSetsEtaSys[cent],
                                                                                                graphRCPCombStatEta[cent], graphRCPCombSystEta[cent],
                                                                                                fileNameEtaRCPOutputWeighting, "pPb_5.023TeV", "EtaRCP", kTRUE,
                                                                                                NULL, fileNameCorrFactors, centArrayCorr[cent]
            );
            if (graphRCPCombCombEta == NULL) {
                cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
                return;
            }
            while (graphRCPCombStatEta[cent]->GetX()[0] < minPtEta[cent]){
                graphRCPCombStatEta[cent]->RemovePoint(0);
            }
            while (graphRCPCombSystEta[cent]->GetX()[0] < minPtEta[cent]){
                graphRCPCombSystEta[cent]->RemovePoint(0);
            }
            while (graphRCPCombCombEta[cent]->GetX()[0] < minPtEta[cent]){
                graphRCPCombCombEta[cent]->RemovePoint(0);
            }
            graphRCPCombCombEta[cent]->Print();

            graphRCPCombStatEtaWOXErr[cent]            = (TGraphAsymmErrors*)graphRCPCombStatEta[cent]->Clone("graphRCPCombStatEtaWOXErr");
            ProduceGraphAsymmWithoutXErrors(graphRCPCombStatEtaWOXErr[cent]);


            // Reading weights from output file for plotting
            ifstream fileWeightsEtaRCPRead;
            fileWeightsEtaRCPRead.open(fileNameEtaRCPOutputWeighting,ios_base::in);
            cout << "reading" << fileNameEtaRCPOutputWeighting << endl;
            Double_t xValuesEtaRCPRead[50];
            Double_t weightsEtaRCPRead[11][50];
            Int_t availableEtaRCPMeas[11]      = { -1, -1, -1, -1, -1,
                                                    -1, -1, -1, -1, -1,
                                                    -1};
            Int_t nMeasSetEtaRCP               = nTotMeasEta[cent];
            Int_t nPtBinsEtaRCPRead            = 0;
            while(!fileWeightsEtaRCPRead.eof() && nPtBinsEtaRCPRead < 50){
                TString garbage             = "";
                if (nPtBinsEtaRCPRead == 0){
                    fileWeightsEtaRCPRead >> garbage ;//>> availableEtaRCPMeas[0] >> availableEtaRCPMeas[1] >> availableEtaRCPMeas[2] >> availableEtaRCPMeas[3];
                    for (Int_t i = 0; i < nMeasSetEtaRCP; i++){
                        fileWeightsEtaRCPRead >> availableEtaRCPMeas[i] ;
                    }
                    cout << "read following measurements: ";
                    for (Int_t meth= 0; meth < 11; meth++){
                        cout << availableEtaRCPMeas[meth] << "\t" ;
                    }
                    cout << endl;
                } else {
                    fileWeightsEtaRCPRead >> xValuesEtaRCPRead[nPtBinsEtaRCPRead-1];
                    for (Int_t i = 0; i < nMeasSetEtaRCP; i++){
                        fileWeightsEtaRCPRead >> weightsEtaRCPRead[availableEtaRCPMeas[i]][nPtBinsEtaRCPRead-1] ;
                    }
                    cout << "read: "<<  nPtBinsEtaRCPRead << "\t"<< xValuesEtaRCPRead[nPtBinsEtaRCPRead-1] << "\t" ;
                    for (Int_t i = 0; i < nMeasSetEtaRCP; i++){
                        cout << weightsEtaRCPRead[availableEtaRCPMeas[i]][nPtBinsEtaRCPRead-1] << "\t";
                    }
                    cout << endl;
                }
                nPtBinsEtaRCPRead++;
            }
            nPtBinsEtaRCPRead                  = nPtBinsEtaRCPRead-2 ;
            fileWeightsEtaRCPRead.close();

            for (Int_t i = 0; i < nMeasSetEtaRCP; i++){
                graphWeightsEtaRCP[cent][availableEtaRCPMeas[i]]      = new TGraph(nPtBinsEtaRCPRead,xValuesEtaRCPRead,weightsEtaRCPRead[availableEtaRCPMeas[i]]);
                Int_t bin = 0;
                for (Int_t n = 0; n< nPtBinsEtaRCPRead; n++){
                    if (graphWeightsEtaRCP[cent][availableEtaRCPMeas[i]]->GetY()[bin] == 0) graphWeightsEtaRCP[cent][availableEtaRCPMeas[i]]->RemovePoint(bin);
                    else bin++;
                }
            }

            // **********************************************************************************************************************
            // ******************************************* Plotting weights method only EMC *****************************************
            // **********************************************************************************************************************
            textSizeLabelsPixel           = 900*0.04;
            canvasWeights->cd();

            histo2DEtaWeights->Draw("copy");
            TLegend* legendWeightsEtaRCP   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetEtaRCP+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
            for (Int_t i = 0; i < nMeasSetEtaRCP; i++){
                DrawGammaSetMarkerTGraph(graphWeightsEtaRCP[cent][availableEtaRCPMeas[i]], markerStyleDet[availableEtaRCPMeas[i]], markerSizeDet[availableEtaRCPMeas[i]]*0.5, colorDet[availableEtaRCPMeas[i]] , colorDet[availableEtaRCPMeas[i]]);
                graphWeightsEtaRCP[cent][availableEtaRCPMeas[i]]->Draw("p,same,z");
                legendWeightsEtaRCP->AddEntry(graphWeightsEtaRCP[cent][availableEtaRCPMeas[i]],nameMeasGlobalLabel[availableEtaRCPMeas[i]],"p");
            }
            legendWeightsEtaRCP->Draw();

            labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelWeightsEnergy->Draw();
            labelWeightsEtaRCP->Draw();

            DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);

            canvasWeights->SaveAs(Form("%s/EtaRCP_Weights_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelSysErr->cd();

            histo2DRelSysErr->GetYaxis()->SetRangeUser(0,49.5);
            histo2DRelSysErr->Draw("copy");
            TLegend* legendRelSysErrRCP       = GetAndSetLegend2(0.60, 0.92-(0.04*(nMeasSetEtaRCP)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
            for (Int_t i = 0; i < nMeasSetEtaRCP; i++){
                if (!sysErrorRelCollectionEtaRCP[cent][availableEtaRCPMeas[i]]) continue;
                DrawGammaSetMarkerTGraph(sysErrorRelCollectionEtaRCP[cent][availableEtaRCPMeas[i]], markerStyleDet[availableEtaRCPMeas[i]], markerSizeDet[availableEtaRCPMeas[i]]*0.5,
                                        colorDet[availableEtaRCPMeas[i]], colorDet[availableEtaRCPMeas[i]]);
                sysErrorRelCollectionEtaRCP[cent][availableEtaRCPMeas[i]]->Draw("p,same,z");
                legendRelSysErrRCP->AddEntry(sysErrorRelCollectionEtaRCP[cent][availableEtaRCPMeas[i]],nameMeasGlobalLabel[availableEtaRCPMeas[i]],"p");
            }
            legendRelSysErrRCP->Draw();

            labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelSysErrEnergy->Draw();
            labelRelSysErrEtaRpPb->Draw();

            canvasRelSysErr->SaveAs(Form("%s/EtaRCP_RelSysErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelSysErrRCP;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelStatErr->cd();

            histo2DRelStatErr->GetYaxis()->SetRangeUser(0,49.5);
            histo2DRelStatErr->Draw("copy");
            TLegend* legendRelStatErrRCP       = GetAndSetLegend2(0.12, 0.92-(0.04*(nMeasSetEtaRCP)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
            for (Int_t i = 0; i < nMeasSetEtaRCP; i++){
                if (!statErrorRelCollectionEtaRCP[cent][availableEtaRCPMeas[i]]) continue;
                DrawGammaSetMarkerTGraph(statErrorRelCollectionEtaRCP[cent][availableEtaRCPMeas[i]], markerStyleDet[availableEtaRCPMeas[i]], markerSizeDet[availableEtaRCPMeas[i]]*0.5,
                                        colorDet[availableEtaRCPMeas[i]], colorDet[availableEtaRCPMeas[i]]);
                statErrorRelCollectionEtaRCP[cent][availableEtaRCPMeas[i]]->Draw("p,same,z");
                legendRelStatErrRCP->AddEntry(statErrorRelCollectionEtaRCP[cent][availableEtaRCPMeas[i]],nameMeasGlobalLabel[availableEtaRCPMeas[i]],"p");
            }
            legendRelStatErrRCP->Draw();

            labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelStatErrEnergy->Draw();
            labelRelStatErrEtaRCP->Draw();

            canvasRelStatErr->SaveAs(Form("%s/EtaRCP_RelStatErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelStatErrRCP;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            graphCombEtaRCPRelStat[cent]      = CalculateRelErrUpAsymmGraph( graphRCPCombStatEta[cent], "relativeStatErrorEtaRCP_Method");
            graphCombEtaRCPRelSys[cent]       = CalculateRelErrUpAsymmGraph( graphRCPCombSystEta[cent], "relativeSysErrorEtaRCP_Method");
            graphCombEtaRCPRelTot[cent]       = CalculateRelErrUpAsymmGraph( graphRCPCombCombEta[cent], "relativeTotalErrorEtaRCP_Method");

            canvasRelTotErr->cd();

            histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,49.5);
            histo2DRelTotErrEta->Draw("copy");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaRCPRelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphCombEtaRCPRelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaRCPRelStat[cent], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
            graphCombEtaRCPRelStat[cent]->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaRCPRelSys[cent], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
            graphCombEtaRCPRelSys[cent]->SetLineStyle(7);
            graphCombEtaRCPRelSys[cent]->Draw("l,x0,same,e1");

            TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
            legendRelTotErr3->AddEntry(graphCombEtaRCPRelTot[cent],"tot","p");
            legendRelTotErr3->AddEntry(graphCombEtaRCPRelStat[cent],"stat","l");
            legendRelTotErr3->AddEntry(graphCombEtaRCPRelSys[cent],"sys","l");
            legendRelTotErr3->Draw();

            labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelTotErrEnergy->Draw();
            labelRelTotErrEtaRCP->Draw();

            canvasRelTotErr->SaveAs(Form("%s/EtaRCP_Reldecomp_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelTotErr3;
        }
    }

    TGraphAsymmErrors* graphRMBIndStatPi0[17][11];
    TGraphAsymmErrors* graphRMBIndStatPi0WOXErr[17][11];
    TH1D* histRMBIndStatPi0[17][11];
    TGraphAsymmErrors* graphRMBIndSystPi0[17][11];
    TGraphAsymmErrors* graphRMBIndCombPi0[17][11];
    TGraph* graphWeightsPi0RMB[17][11];
    TGraphAsymmErrors* statErrorRelCollectionPi0RMB[17][11];
    TGraphAsymmErrors* sysErrorRelCollectionPi0RMB[17][11];
    TGraphAsymmErrors* graphRMBCombStatPi0[17];
    TGraphAsymmErrors* graphRMBCombSystPi0[17];
    TGraphAsymmErrors* graphRMBCombStatPi0WOXErr[17];
    TGraphAsymmErrors* graphRMBCombCombPi0[17];
    TGraphAsymmErrors* graphCombPi0RMBRelStat[17];
    TGraphAsymmErrors* graphCombPi0RMBRelSys[17];
    TGraphAsymmErrors* graphCombPi0RMBRelTot[17];


    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if ( !enableCentRMB[cent] ) continue;
        rMBNColl[cent]      = nCollpPb[cent]/(nCollpPb[rMBRefInt[cent]]);
        rMBNCollErr[cent]   = rMBNColl[cent]* TMath::Sqrt(TMath::Power(nCollErrpPb[cent]/nCollpPb[cent],2) + TMath::Power(nCollErrpPb[rMBRefInt[cent]]/nCollpPb[rMBRefInt[cent]],2));
        clog << nCollpPb[cent] << "\t" << nCollpPb[rMBRefInt[cent]] << endl;
        clog << rMBNColl[cent] << "\t" << rMBNCollErr[cent] << endl;

        // *******************************************************************************************************
        // *************************** RMB calculation for pi0 for individual measurements **********************
        // *******************************************************************************************************
        for (Int_t meth = 0; meth < 11; meth++){
            cout << "RMB calc for pi0 " << centArray[cent].Data() << addCentString[cent].Data() << "\t" << runArray[cent].Data() << "\t" << nameMeasGlobal[meth].Data() << endl;
            graphWeightsPi0RMB[cent][meth]                              = NULL;
            graphRMBIndStatPi0[cent][meth]                              = NULL;
            graphRMBIndStatPi0WOXErr[cent][meth]                        = NULL;
            histRMBIndStatPi0[cent][meth]                               = NULL;
            graphRMBIndSystPi0[cent][meth]                              = NULL;
            graphRMBIndCombPi0[cent][meth]                              = NULL;
            statErrorRelCollectionPi0RMB[cent][meth]                     = NULL;
            sysErrorRelCollectionPi0RMB[cent][meth]                      = NULL;
            if (havePi0SysDetailedmin[cent][meth] && graphIndPi0InvYieldStat_yShifted[rMBRefInt[cent]][meth] && graphIndPi0InvYieldStat_yShifted[cent][meth] ){ //&& !( i == 3 )
                if(bWCorrection.Contains("Y")  ){
                    TGraphAsymmErrors* tempGraphPeriStat                = (TGraphAsymmErrors*)graphIndPi0InvYieldStat_yShifted[rMBRefInt[cent]][meth]->Clone("tempPeriStat");
                    TGraphAsymmErrors* tempGraphPeriSys                 = (TGraphAsymmErrors*)graphIndPi0InvYieldSys_yShifted[rMBRefInt[cent]][meth]->Clone("tempPeriSys");
                    tempGraphPeriStat                                   = ScaleGraph(tempGraphPeriStat, 1./tpPb[rMBRefInt[cent]]);
                    tempGraphPeriSys                                    = ScaleGraph(tempGraphPeriSys, 1./tpPb[rMBRefInt[cent]]);
                    RecalculateErrorsBasedOnDetailedInputFile( tempGraphPeriSys, fileNamesminPi0DetailedSys[cent][meth] );
                    TGraphAsymmErrors* tempGraphCentStat                = (TGraphAsymmErrors*)graphIndPi0InvYieldStat_yShifted[cent][meth]->Clone("tempCentStat");
                    TGraphAsymmErrors* tempGraphCentSys                 = (TGraphAsymmErrors*)graphIndPi0InvYieldSys_yShifted[cent][meth]->Clone("tempCentSys");
                    tempGraphCentStat                                   = ScaleGraph(tempGraphCentStat, 1./tpPb[cent]);
                    tempGraphCentSys                                    = ScaleGraph(tempGraphCentSys, 1./tpPb[cent]);
                    if (tempGraphPeriStat->GetX()[0] != tempGraphCentStat->GetX()[0])
                        continue;

                    graphRMBIndCombPi0[cent][meth]                      = CalcRpPbV2(   tempGraphPeriStat,
                                                                                        tempGraphPeriSys,
                                                                                        tempGraphCentStat,
                                                                                        tempGraphCentSys,
                                                                                        &graphRMBIndStatPi0[cent][meth],
                                                                                        &graphRMBIndSystPi0[cent][meth],
                                                                                        1, 0 ,
                                                                                        ptSysRemRMBNames[cent][meth], fileNamespPbPi0DetailedSys[cent][meth],
                                                                                        fileNamesminPi0DetailedSys[cent][meth], fileNamesRMBPi0DetailedSys[cent][meth]
                                                                                    );
                    graphRMBIndCombPi0[cent][meth]->Print();
                    histRMBIndStatPi0[cent][meth]                       = (TH1D*)statErrorCollectionPi0[cent][meth]->Clone(Form("histRMB%sStatPi0",nameMeasGlobalLabel[meth].Data()));
                    for (Int_t j = 0; j< histRMBIndStatPi0[cent][meth]->GetNbinsX()+1; j++){
                        histRMBIndStatPi0[cent][meth]->SetBinContent(j,0);
                        histRMBIndStatPi0[cent][meth]->SetBinError(j,0);
                    }
                    for (Int_t j = 0; j < graphRMBIndStatPi0[cent][meth]->GetN(); j++){
                        Int_t bin                               = histRMBIndStatPi0[cent][meth]->GetXaxis()->FindBin(graphRMBIndStatPi0[cent][meth]->GetX()[j]);
                        histRMBIndStatPi0[cent][meth]->SetBinContent(bin,graphRMBIndStatPi0[cent][meth]->GetY()[j]);
                        histRMBIndStatPi0[cent][meth]->SetBinError(bin,graphRMBIndStatPi0[cent][meth]->GetEYlow()[j]);
                    }
                }

                if (histRMBIndStatPi0[cent][meth]){
                    statErrorRelCollectionPi0RMB[cent][meth]   = (TGraphAsymmErrors*)graphRMBIndStatPi0[cent][meth]->Clone(Form("relativeStatErrorPi0RMB%i_%s", cent, nameMeasGlobal[meth].Data()));
                    RemoveZerosAtBeginningAndEndFromGraph(statErrorRelCollectionPi0RMB[cent][meth]);
                    statErrorRelCollectionPi0RMB[cent][meth]   = CalculateRelErrUpAsymmGraph( statErrorRelCollectionPi0RMB[cent][meth], Form("relativeStatErrorPi0RMB_%i_%s", cent,nameMeasGlobal[meth].Data()));
                    statErrorRelCollectionPi0RMB[cent][meth]->Print();
                    graphRMBIndStatPi0WOXErr[cent][meth]      = (TGraphAsymmErrors*)graphRMBIndStatPi0[cent][meth]->Clone(Form("graphRMB%sStatPi0WOXErr%i", nameMeasGlobalLabel[meth].Data(), cent));

                }
                if (graphRMBIndSystPi0[cent][meth]){
                    sysErrorRelCollectionPi0RMB[cent][meth]    = (TGraphAsymmErrors*)graphRMBIndSystPi0[cent][meth]->Clone(Form("relativeSysErrorPi0RMB%i_%s", cent, nameMeasGlobal[meth].Data()));
                    RemoveZerosAtBeginningAndEndFromGraph(sysErrorRelCollectionPi0RMB[cent][meth]);
                    sysErrorRelCollectionPi0RMB[cent][meth]    = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionPi0RMB[cent][meth], Form("relativeSysErrorPi0RMB%i_%s", cent, nameMeasGlobal[meth].Data()));
                }
            }
        }
        // *******************************************************************************************************
        // *************************** RMB calculation for pi0 **************************************************
        // *******************************************************************************************************
        graphRMBCombStatPi0[cent]                      = NULL;
        graphRMBCombStatPi0[cent]                      = NULL;
        graphRMBCombStatPi0[cent]                      = NULL;

        if (!enableCentComb[cent]) continue;
        if(bWCorrection.Contains("Y") ){
            TString fileNamePi0RMBOutputWeighting      = Form("%s/Pi0RMB_WeightingMethod_%s%s%s.dat", outputDirSupportComb.Data(), centArrayOutput[cent].Data(),
                                                               addCentString[cent].Data(), runArray[cent].Data());
            graphRMBCombCombPi0[cent]                  = CombinePtPointsSpectraFullCorrMat(    histRMBIndStatPi0[cent],    graphRMBIndSystPi0[cent],
                                                                                                xPtLimitsPi0[cent], maxNBinsPi0[cent],
                                                                                                offSetsPi0[cent], offSetsPi0Sys[cent],
                                                                                                graphRMBCombStatPi0[cent], graphRMBCombSystPi0[cent],
                                                                                                fileNamePi0RMBOutputWeighting, "pPb_5.023TeV", "Pi0RMB", kTRUE,
                                                                                                NULL, fileNameCorrFactors, centArrayCorr[cent]
            );
            if (graphRMBCombCombPi0 == NULL) {
                cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
                return;
            }
            while (graphRMBCombStatPi0[cent]->GetX()[0] < minPtPi0[cent]){
                graphRMBCombStatPi0[cent]->RemovePoint(0);
            }
            while (graphRMBCombSystPi0[cent]->GetX()[0] < minPtPi0[cent]){
                graphRMBCombSystPi0[cent]->RemovePoint(0);
            }
            while (graphRMBCombCombPi0[cent]->GetX()[0] < minPtPi0[cent]){
                graphRMBCombCombPi0[cent]->RemovePoint(0);
            }
            graphRMBCombCombPi0[cent]->Print();

            graphRMBCombStatPi0WOXErr[cent]            = (TGraphAsymmErrors*)graphRMBCombStatPi0[cent]->Clone("graphRMBCombStatPi0WOXErr");
            ProduceGraphAsymmWithoutXErrors(graphRMBCombStatPi0WOXErr[cent]);


            // Reading weights from output file for plotting
            ifstream fileWeightsPi0RMBRead;
            fileWeightsPi0RMBRead.open(fileNamePi0RMBOutputWeighting,ios_base::in);
            cout << "reading" << fileNamePi0RMBOutputWeighting << endl;
            Double_t xValuesPi0RMBRead[50];
            Double_t weightsPi0RMBRead[11][50];
            Int_t availablePi0RMBMeas[11]      = { -1, -1, -1, -1, -1,
                                                    -1, -1, -1, -1, -1,
                                                    -1};
            Int_t nMeasSetPi0RMB               = nTotMeasPi0[cent];
            Int_t nPtBinsPi0RMBRead            = 0;
            while(!fileWeightsPi0RMBRead.eof() && nPtBinsPi0RMBRead < 50){
                TString garbage             = "";
                if (nPtBinsPi0RMBRead == 0){
                    fileWeightsPi0RMBRead >> garbage ;//>> availablePi0RMBMeas[0] >> availablePi0RMBMeas[1] >> availablePi0RMBMeas[2] >> availablePi0RMBMeas[3];
                    for (Int_t i = 0; i < nMeasSetPi0RMB; i++){
                        fileWeightsPi0RMBRead >> availablePi0RMBMeas[i] ;
                    }
                    cout << "read following measurements: ";
                    for (Int_t meth= 0; meth < 11; meth++){
                        cout << availablePi0RMBMeas[meth] << "\t" ;
                    }
                    cout << endl;
                } else {
                    fileWeightsPi0RMBRead >> xValuesPi0RMBRead[nPtBinsPi0RMBRead-1];
                    for (Int_t i = 0; i < nMeasSetPi0RMB; i++){
                        fileWeightsPi0RMBRead >> weightsPi0RMBRead[availablePi0RMBMeas[i]][nPtBinsPi0RMBRead-1] ;
                    }
                    cout << "read: "<<  nPtBinsPi0RMBRead << "\t"<< xValuesPi0RMBRead[nPtBinsPi0RMBRead-1] << "\t" ;
                    for (Int_t i = 0; i < nMeasSetPi0RMB; i++){
                        cout << weightsPi0RMBRead[availablePi0RMBMeas[i]][nPtBinsPi0RMBRead-1] << "\t";
                    }
                    cout << endl;
                }
                nPtBinsPi0RMBRead++;
            }
            nPtBinsPi0RMBRead                  = nPtBinsPi0RMBRead-2 ;
            fileWeightsPi0RMBRead.close();

            for (Int_t i = 0; i < nMeasSetPi0RMB; i++){
                graphWeightsPi0RMB[cent][availablePi0RMBMeas[i]]      = new TGraph(nPtBinsPi0RMBRead,xValuesPi0RMBRead,weightsPi0RMBRead[availablePi0RMBMeas[i]]);
                Int_t bin = 0;
                for (Int_t n = 0; n< nPtBinsPi0RMBRead; n++){
                    if (graphWeightsPi0RMB[cent][availablePi0RMBMeas[i]]->GetY()[bin] == 0) graphWeightsPi0RMB[cent][availablePi0RMBMeas[i]]->RemovePoint(bin);
                    else bin++;
                }
            }

            // **********************************************************************************************************************
            // ******************************************* Plotting weights method only EMC *****************************************
            // **********************************************************************************************************************
            textSizeLabelsPixel           = 900*0.04;
            canvasWeights->cd();

            histo2DPi0Weights->Draw("copy");
            TLegend* legendWeightsPi0RMB   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetPi0RMB+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
            for (Int_t i = 0; i < nMeasSetPi0RMB; i++){
                DrawGammaSetMarkerTGraph(graphWeightsPi0RMB[cent][availablePi0RMBMeas[i]], markerStyleDet[availablePi0RMBMeas[i]], markerSizeDet[availablePi0RMBMeas[i]]*0.5, colorDet[availablePi0RMBMeas[i]] , colorDet[availablePi0RMBMeas[i]]);
                graphWeightsPi0RMB[cent][availablePi0RMBMeas[i]]->Draw("p,same,z");
                legendWeightsPi0RMB->AddEntry(graphWeightsPi0RMB[cent][availablePi0RMBMeas[i]],nameMeasGlobalLabel[availablePi0RMBMeas[i]],"p");
            }
            legendWeightsPi0RMB->Draw();

            labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelWeightsEnergy->Draw();
            labelWeightsPi0RMB->Draw();

            DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);

            canvasWeights->SaveAs(Form("%s/Pi0RMB_Weights_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelSysErr->cd();

            histo2DRelSysErr->GetYaxis()->SetRangeUser(0,49.5);
            histo2DRelSysErr->Draw("copy");
            TLegend* legendRelSysErrRMB       = GetAndSetLegend2(0.60, 0.92-(0.04*(nMeasSetPi0RMB)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
            for (Int_t i = 0; i < nMeasSetPi0RMB; i++){
                if (!sysErrorRelCollectionPi0RMB[cent][availablePi0RMBMeas[i]]) continue;
                DrawGammaSetMarkerTGraph(sysErrorRelCollectionPi0RMB[cent][availablePi0RMBMeas[i]], markerStyleDet[availablePi0RMBMeas[i]], markerSizeDet[availablePi0RMBMeas[i]]*0.5,
                                        colorDet[availablePi0RMBMeas[i]], colorDet[availablePi0RMBMeas[i]]);
                sysErrorRelCollectionPi0RMB[cent][availablePi0RMBMeas[i]]->Draw("p,same,z");
                legendRelSysErrRMB->AddEntry(sysErrorRelCollectionPi0RMB[cent][availablePi0RMBMeas[i]],nameMeasGlobalLabel[availablePi0RMBMeas[i]],"p");
            }
            legendRelSysErrRMB->Draw();

            labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelSysErrEnergy->Draw();
            labelRelSysErrPi0RpPb->Draw();

            canvasRelSysErr->SaveAs(Form("%s/Pi0RMB_RelSysErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelSysErrRMB;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelStatErr->cd();

            histo2DRelStatErr->GetYaxis()->SetRangeUser(0,49.5);
            histo2DRelStatErr->Draw("copy");
            TLegend* legendRelStatErrRMB       = GetAndSetLegend2(0.12, 0.92-(0.04*(nMeasSetPi0RMB)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
            for (Int_t i = 0; i < nMeasSetPi0RMB; i++){
                if (!statErrorRelCollectionPi0RMB[cent][availablePi0RMBMeas[i]]) continue;
                DrawGammaSetMarkerTGraph(statErrorRelCollectionPi0RMB[cent][availablePi0RMBMeas[i]], markerStyleDet[availablePi0RMBMeas[i]], markerSizeDet[availablePi0RMBMeas[i]]*0.5,
                                        colorDet[availablePi0RMBMeas[i]], colorDet[availablePi0RMBMeas[i]]);
                statErrorRelCollectionPi0RMB[cent][availablePi0RMBMeas[i]]->Draw("p,same,z");
                legendRelStatErrRMB->AddEntry(statErrorRelCollectionPi0RMB[cent][availablePi0RMBMeas[i]],nameMeasGlobalLabel[availablePi0RMBMeas[i]],"p");
            }
            legendRelStatErrRMB->Draw();

            labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelStatErrEnergy->Draw();
            labelRelStatErrPi0RMB->Draw();

            canvasRelStatErr->SaveAs(Form("%s/Pi0RMB_RelStatErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelStatErrRMB;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            graphCombPi0RMBRelStat[cent]      = CalculateRelErrUpAsymmGraph( graphRMBCombStatPi0[cent], "relativeStatErrorPi0RMB_Method");
            graphCombPi0RMBRelSys[cent]       = CalculateRelErrUpAsymmGraph( graphRMBCombSystPi0[cent], "relativeSysErrorPi0RMB_Method");
            graphCombPi0RMBRelTot[cent]       = CalculateRelErrUpAsymmGraph( graphRMBCombCombPi0[cent], "relativeTotalErrorPi0RMB_Method");

            canvasRelTotErr->cd();

            histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,49.5);
            histo2DRelTotErrPi0->Draw("copy");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0RMBRelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphCombPi0RMBRelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0RMBRelStat[cent], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
            graphCombPi0RMBRelStat[cent]->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0RMBRelSys[cent], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
            graphCombPi0RMBRelSys[cent]->SetLineStyle(7);
            graphCombPi0RMBRelSys[cent]->Draw("l,x0,same,e1");

            TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
            legendRelTotErr3->AddEntry(graphCombPi0RMBRelTot[cent],"tot","p");
            legendRelTotErr3->AddEntry(graphCombPi0RMBRelStat[cent],"stat","l");
            legendRelTotErr3->AddEntry(graphCombPi0RMBRelSys[cent],"sys","l");
            legendRelTotErr3->Draw();

            labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelTotErrEnergy->Draw();
            labelRelTotErrPi0RMB->Draw();

            canvasRelTotErr->SaveAs(Form("%s/Pi0RMB_Reldecomp_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelTotErr3;
        }
    }

    TGraphAsymmErrors* graphRMBIndStatEta[17][11];
    TGraphAsymmErrors* graphRMBIndStatEtaWOXErr[17][11];
    TH1D* histRMBIndStatEta[17][11];
    TGraphAsymmErrors* graphRMBIndSystEta[17][11];
    TGraphAsymmErrors* graphRMBIndCombEta[17][11];
    TGraph* graphWeightsEtaRMB[17][11];
    TGraphAsymmErrors* statErrorRelCollectionEtaRMB[17][11];
    TGraphAsymmErrors* sysErrorRelCollectionEtaRMB[17][11];
    TGraphAsymmErrors* graphRMBCombStatEta[17];
    TGraphAsymmErrors* graphRMBCombSystEta[17];
    TGraphAsymmErrors* graphRMBCombStatEtaWOXErr[17];
    TGraphAsymmErrors* graphRMBCombCombEta[17];
    TGraphAsymmErrors* graphCombEtaRMBRelStat[17];
    TGraphAsymmErrors* graphCombEtaRMBRelSys[17];
    TGraphAsymmErrors* graphCombEtaRMBRelTot[17];


    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if ( !enableCentRMB[cent] ) continue;
        // *******************************************************************************************************
        // *************************** RMB calculation for pi0 for individual measurements **********************
        // *******************************************************************************************************
        for (Int_t meth = 0; meth < 11; meth++){
            cout << "RMB calc for pi0 " << centArray[cent].Data() << addCentString[cent].Data() << "\t" << runArray[cent].Data() << "\t" << nameMeasGlobal[meth].Data() << endl;
            graphWeightsEtaRMB[cent][meth]                              = NULL;
            graphRMBIndStatEta[cent][meth]                              = NULL;
            graphRMBIndStatEtaWOXErr[cent][meth]                        = NULL;
            histRMBIndStatEta[cent][meth]                               = NULL;
            graphRMBIndSystEta[cent][meth]                              = NULL;
            graphRMBIndCombEta[cent][meth]                              = NULL;
            statErrorRelCollectionEtaRMB[cent][meth]                     = NULL;
            sysErrorRelCollectionEtaRMB[cent][meth]                      = NULL;
            if (haveEtaSysDetailedmin[cent][meth] && graphIndEtaInvYieldStat_yShifted[rMBRefInt[cent]][meth] && graphIndEtaInvYieldStat_yShifted[cent][meth] ){ //&& !( i == 3 )
                if(bWCorrection.Contains("Y")  ){
                    TGraphAsymmErrors* tempGraphPeriStat                = (TGraphAsymmErrors*)graphIndEtaInvYieldStat_yShifted[rMBRefInt[cent]][meth]->Clone("tempPeriStat");
                    TGraphAsymmErrors* tempGraphPeriSys                 = (TGraphAsymmErrors*)graphIndEtaInvYieldSys_yShifted[rMBRefInt[cent]][meth]->Clone("tempPeriSys");
                    tempGraphPeriStat                                   = ScaleGraph(tempGraphPeriStat, 1./tpPb[rMBRefInt[cent]]);
                    tempGraphPeriSys                                    = ScaleGraph(tempGraphPeriSys, 1./tpPb[rMBRefInt[cent]]);
                    RecalculateErrorsBasedOnDetailedInputFile( tempGraphPeriSys, fileNamesminEtaDetailedSys[cent][meth] );
                    TGraphAsymmErrors* tempGraphCentStat                = (TGraphAsymmErrors*)graphIndEtaInvYieldStat_yShifted[cent][meth]->Clone("tempCentStat");
                    TGraphAsymmErrors* tempGraphCentSys                 = (TGraphAsymmErrors*)graphIndEtaInvYieldSys_yShifted[cent][meth]->Clone("tempCentSys");
                    tempGraphCentStat                                   = ScaleGraph(tempGraphCentStat, 1./tpPb[cent]);
                    tempGraphCentSys                                    = ScaleGraph(tempGraphCentSys, 1./tpPb[cent]);
                    if (tempGraphPeriStat->GetX()[0] != tempGraphCentStat->GetX()[0])
                        continue;

                    graphRMBIndCombEta[cent][meth]                      = CalcRpPbV2(   tempGraphPeriStat,
                                                                                        tempGraphPeriSys,
                                                                                        tempGraphCentStat,
                                                                                        tempGraphCentSys,
                                                                                        &graphRMBIndStatEta[cent][meth],
                                                                                        &graphRMBIndSystEta[cent][meth],
                                                                                        1, 0 ,
                                                                                        ptSysRemRMBNames[cent][meth], fileNamespPbEtaDetailedSys[cent][meth],
                                                                                        fileNamesminEtaDetailedSys[cent][meth], fileNamesRMBEtaDetailedSys[cent][meth]
                    );
                    graphRMBIndCombEta[cent][meth]->Print();
                    histRMBIndStatEta[cent][meth]                       = (TH1D*)statErrorCollectionEta[cent][meth]->Clone(Form("histRMB%sStatEta",nameMeasGlobalLabel[meth].Data()));
                    for (Int_t j = 0; j< histRMBIndStatEta[cent][meth]->GetNbinsX()+1; j++){
                        histRMBIndStatEta[cent][meth]->SetBinContent(j,0);
                        histRMBIndStatEta[cent][meth]->SetBinError(j,0);
                    }
                    for (Int_t j = 0; j < graphRMBIndStatEta[cent][meth]->GetN(); j++){
                        Int_t bin                               = histRMBIndStatEta[cent][meth]->GetXaxis()->FindBin(graphRMBIndStatEta[cent][meth]->GetX()[j]);
                        histRMBIndStatEta[cent][meth]->SetBinContent(bin,graphRMBIndStatEta[cent][meth]->GetY()[j]);
                        histRMBIndStatEta[cent][meth]->SetBinError(bin,graphRMBIndStatEta[cent][meth]->GetEYlow()[j]);
                    }
                }

                if (histRMBIndStatEta[cent][meth]){
                    statErrorRelCollectionEtaRMB[cent][meth]   = (TGraphAsymmErrors*)graphRMBIndStatEta[cent][meth]->Clone(Form("relativeStatErrorEtaRMB%i_%s", cent, nameMeasGlobal[meth].Data()));
                    RemoveZerosAtBeginningAndEndFromGraph(statErrorRelCollectionEtaRMB[cent][meth]);
                    statErrorRelCollectionEtaRMB[cent][meth]   = CalculateRelErrUpAsymmGraph( statErrorRelCollectionEtaRMB[cent][meth], Form("relativeStatErrorEtaRMB_%i_%s", cent,nameMeasGlobal[meth].Data()));
                    statErrorRelCollectionEtaRMB[cent][meth]->Print();
                    graphRMBIndStatEtaWOXErr[cent][meth]      = (TGraphAsymmErrors*)graphRMBIndStatEta[cent][meth]->Clone(Form("graphRMB%sStatEtaWOXErr%i", nameMeasGlobalLabel[meth].Data(), cent));

                }
                if (graphRMBIndSystEta[cent][meth]){
                    sysErrorRelCollectionEtaRMB[cent][meth]    = (TGraphAsymmErrors*)graphRMBIndSystEta[cent][meth]->Clone(Form("relativeSysErrorEtaRMB%i_%s", cent, nameMeasGlobal[meth].Data()));
                    RemoveZerosAtBeginningAndEndFromGraph(sysErrorRelCollectionEtaRMB[cent][meth]);
                    sysErrorRelCollectionEtaRMB[cent][meth]    = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionEtaRMB[cent][meth], Form("relativeSysErrorEtaRMB%i_%s", cent, nameMeasGlobal[meth].Data()));
                }
            }
        }
        // *******************************************************************************************************
        // *************************** RMB calculation for eta **************************************************
        // *******************************************************************************************************
        graphRMBCombStatEta[cent]                      = NULL;
        graphRMBCombStatEta[cent]                      = NULL;
        graphRMBCombStatEta[cent]                      = NULL;

        if (!enableCentComb[cent]) continue;
        if(bWCorrection.Contains("Y") ){
            TString fileNameEtaRMBOutputWeighting      = Form("%s/EtaRMB_WeightingMethod_%s%s%s.dat", outputDirSupportComb.Data(), centArrayOutput[cent].Data(),
                                                               addCentString[cent].Data(), runArray[cent].Data());
            graphRMBCombCombEta[cent]                  = CombinePtPointsSpectraFullCorrMat(    histRMBIndStatEta[cent],    graphRMBIndSystEta[cent],
                                                                                                xPtLimitsEta[cent], maxNBinsEta[cent],
                                                                                                offSetsEta[cent], offSetsEtaSys[cent],
                                                                                                graphRMBCombStatEta[cent], graphRMBCombSystEta[cent],
                                                                                                fileNameEtaRMBOutputWeighting, "pPb_5.023TeV", "EtaRMB", kTRUE,
                                                                                                NULL, fileNameCorrFactors, centArrayCorr[cent]
            );
            if (graphRMBCombCombEta == NULL) {
                cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
                return;
            }
            while (graphRMBCombStatEta[cent]->GetX()[0] < minPtEta[cent]){
                graphRMBCombStatEta[cent]->RemovePoint(0);
            }
            while (graphRMBCombSystEta[cent]->GetX()[0] < minPtEta[cent]){
                graphRMBCombSystEta[cent]->RemovePoint(0);
            }
            while (graphRMBCombCombEta[cent]->GetX()[0] < minPtEta[cent]){
                graphRMBCombCombEta[cent]->RemovePoint(0);
            }
            graphRMBCombCombEta[cent]->Print();

            graphRMBCombStatEtaWOXErr[cent]            = (TGraphAsymmErrors*)graphRMBCombStatEta[cent]->Clone("graphRMBCombStatEtaWOXErr");
            ProduceGraphAsymmWithoutXErrors(graphRMBCombStatEtaWOXErr[cent]);


            // Reading weights from output file for plotting
            ifstream fileWeightsEtaRMBRead;
            fileWeightsEtaRMBRead.open(fileNameEtaRMBOutputWeighting,ios_base::in);
            cout << "reading" << fileNameEtaRMBOutputWeighting << endl;
            Double_t xValuesEtaRMBRead[50];
            Double_t weightsEtaRMBRead[11][50];
            Int_t availableEtaRMBMeas[11]      = { -1, -1, -1, -1, -1,
                                                    -1, -1, -1, -1, -1,
                                                    -1};
            Int_t nMeasSetEtaRMB               = nTotMeasEta[cent];
            Int_t nPtBinsEtaRMBRead            = 0;
            while(!fileWeightsEtaRMBRead.eof() && nPtBinsEtaRMBRead < 50){
                TString garbage             = "";
                if (nPtBinsEtaRMBRead == 0){
                    fileWeightsEtaRMBRead >> garbage ;//>> availableEtaRMBMeas[0] >> availableEtaRMBMeas[1] >> availableEtaRMBMeas[2] >> availableEtaRMBMeas[3];
                    for (Int_t i = 0; i < nMeasSetEtaRMB; i++){
                        fileWeightsEtaRMBRead >> availableEtaRMBMeas[i] ;
                    }
                    cout << "read following measurements: ";
                    for (Int_t meth= 0; meth < 11; meth++){
                        cout << availableEtaRMBMeas[meth] << "\t" ;
                    }
                    cout << endl;
                } else {
                    fileWeightsEtaRMBRead >> xValuesEtaRMBRead[nPtBinsEtaRMBRead-1];
                    for (Int_t i = 0; i < nMeasSetEtaRMB; i++){
                        fileWeightsEtaRMBRead >> weightsEtaRMBRead[availableEtaRMBMeas[i]][nPtBinsEtaRMBRead-1] ;
                    }
                    cout << "read: "<<  nPtBinsEtaRMBRead << "\t"<< xValuesEtaRMBRead[nPtBinsEtaRMBRead-1] << "\t" ;
                    for (Int_t i = 0; i < nMeasSetEtaRMB; i++){
                        cout << weightsEtaRMBRead[availableEtaRMBMeas[i]][nPtBinsEtaRMBRead-1] << "\t";
                    }
                    cout << endl;
                }
                nPtBinsEtaRMBRead++;
            }
            nPtBinsEtaRMBRead                  = nPtBinsEtaRMBRead-2 ;
            fileWeightsEtaRMBRead.close();

            for (Int_t i = 0; i < nMeasSetEtaRMB; i++){
                graphWeightsEtaRMB[cent][availableEtaRMBMeas[i]]      = new TGraph(nPtBinsEtaRMBRead,xValuesEtaRMBRead,weightsEtaRMBRead[availableEtaRMBMeas[i]]);
                Int_t bin = 0;
                for (Int_t n = 0; n< nPtBinsEtaRMBRead; n++){
                    if (graphWeightsEtaRMB[cent][availableEtaRMBMeas[i]]->GetY()[bin] == 0) graphWeightsEtaRMB[cent][availableEtaRMBMeas[i]]->RemovePoint(bin);
                    else bin++;
                }
            }

            // **********************************************************************************************************************
            // ******************************************* Plotting weights method only EMC *****************************************
            // **********************************************************************************************************************
            textSizeLabelsPixel           = 900*0.04;
            canvasWeights->cd();

            histo2DEtaWeights->Draw("copy");
            TLegend* legendWeightsEtaRMB   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetEtaRMB+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
            for (Int_t i = 0; i < nMeasSetEtaRMB; i++){
                DrawGammaSetMarkerTGraph(graphWeightsEtaRMB[cent][availableEtaRMBMeas[i]], markerStyleDet[availableEtaRMBMeas[i]], markerSizeDet[availableEtaRMBMeas[i]]*0.5, colorDet[availableEtaRMBMeas[i]] , colorDet[availableEtaRMBMeas[i]]);
                graphWeightsEtaRMB[cent][availableEtaRMBMeas[i]]->Draw("p,same,z");
                legendWeightsEtaRMB->AddEntry(graphWeightsEtaRMB[cent][availableEtaRMBMeas[i]],nameMeasGlobalLabel[availableEtaRMBMeas[i]],"p");
            }
            legendWeightsEtaRMB->Draw();

            labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelWeightsEnergy->Draw();
            labelWeightsEtaRMB->Draw();

            DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);

            canvasWeights->SaveAs(Form("%s/EtaRMB_Weights_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelSysErr->cd();

            histo2DRelSysErr->GetYaxis()->SetRangeUser(0,49.5);
            histo2DRelSysErr->Draw("copy");
            TLegend* legendRelSysErrRMB       = GetAndSetLegend2(0.60, 0.92-(0.04*(nMeasSetEtaRMB)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
            for (Int_t i = 0; i < nMeasSetEtaRMB; i++){
                if (!sysErrorRelCollectionEtaRMB[cent][availableEtaRMBMeas[i]]) continue;
                DrawGammaSetMarkerTGraph(sysErrorRelCollectionEtaRMB[cent][availableEtaRMBMeas[i]], markerStyleDet[availableEtaRMBMeas[i]], markerSizeDet[availableEtaRMBMeas[i]]*0.5,
                                        colorDet[availableEtaRMBMeas[i]], colorDet[availableEtaRMBMeas[i]]);
                sysErrorRelCollectionEtaRMB[cent][availableEtaRMBMeas[i]]->Draw("p,same,z");
                legendRelSysErrRMB->AddEntry(sysErrorRelCollectionEtaRMB[cent][availableEtaRMBMeas[i]],nameMeasGlobalLabel[availableEtaRMBMeas[i]],"p");
            }
            legendRelSysErrRMB->Draw();

            labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelSysErrEnergy->Draw();
            labelRelSysErrEtaRpPb->Draw();

            canvasRelSysErr->SaveAs(Form("%s/EtaRMB_RelSysErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelSysErrRMB;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelStatErr->cd();

            histo2DRelStatErr->GetYaxis()->SetRangeUser(0,49.5);
            histo2DRelStatErr->Draw("copy");
            TLegend* legendRelStatErrRMB       = GetAndSetLegend2(0.12, 0.92-(0.04*(nMeasSetEtaRMB)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
            for (Int_t i = 0; i < nMeasSetEtaRMB; i++){
                if (!statErrorRelCollectionEtaRMB[cent][availableEtaRMBMeas[i]]) continue;
                DrawGammaSetMarkerTGraph(statErrorRelCollectionEtaRMB[cent][availableEtaRMBMeas[i]], markerStyleDet[availableEtaRMBMeas[i]], markerSizeDet[availableEtaRMBMeas[i]]*0.5,
                                        colorDet[availableEtaRMBMeas[i]], colorDet[availableEtaRMBMeas[i]]);
                statErrorRelCollectionEtaRMB[cent][availableEtaRMBMeas[i]]->Draw("p,same,z");
                legendRelStatErrRMB->AddEntry(statErrorRelCollectionEtaRMB[cent][availableEtaRMBMeas[i]],nameMeasGlobalLabel[availableEtaRMBMeas[i]],"p");
            }
            legendRelStatErrRMB->Draw();

            labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelStatErrEnergy->Draw();
            labelRelStatErrEtaRMB->Draw();

            canvasRelStatErr->SaveAs(Form("%s/EtaRMB_RelStatErr_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelStatErrRMB;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            graphCombEtaRMBRelStat[cent]      = CalculateRelErrUpAsymmGraph( graphRMBCombStatEta[cent], "relativeStatErrorEtaRMB_Method");
            graphCombEtaRMBRelSys[cent]       = CalculateRelErrUpAsymmGraph( graphRMBCombSystEta[cent], "relativeSysErrorEtaRMB_Method");
            graphCombEtaRMBRelTot[cent]       = CalculateRelErrUpAsymmGraph( graphRMBCombCombEta[cent], "relativeTotalErrorEtaRMB_Method");

            canvasRelTotErr->cd();

            histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,49.5);
            histo2DRelTotErrEta->Draw("copy");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaRMBRelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphCombEtaRMBRelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaRMBRelStat[cent], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
            graphCombEtaRMBRelStat[cent]->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaRMBRelSys[cent], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
            graphCombEtaRMBRelSys[cent]->SetLineStyle(7);
            graphCombEtaRMBRelSys[cent]->Draw("l,x0,same,e1");

            TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
            legendRelTotErr3->AddEntry(graphCombEtaRMBRelTot[cent],"tot","p");
            legendRelTotErr3->AddEntry(graphCombEtaRMBRelStat[cent],"stat","l");
            legendRelTotErr3->AddEntry(graphCombEtaRMBRelSys[cent],"sys","l");
            legendRelTotErr3->Draw();

            labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelRelTotErrEnergy->Draw();
            labelRelTotErrEtaRMB->Draw();

            canvasRelTotErr->SaveAs(Form("%s/EtaRMB_Reldecomp_%s%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data() ,suffix.Data()));
            delete legendRelTotErr3;
        }
    }

    //===================================================================================================================
    //===================================================================================================================
    //==================== STARTING MAIN PLOTTING FOR PAPER =============================================================
    //===================================================================================================================
    //===================================================================================================================


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
//     cout << textsizeLabelsWidth << endl;

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
            labelMassEnergy->SetText(0.13,0.775,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
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
        canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth_%s%s%s.%s",outputDirSupportPaper.Data(),centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));

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
                legendMass->AddEntry(graphPi0Mass[cent][meth],centArray2[cent],"p");
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

    TLatex *labelMassEta            = new TLatex(0.13,0.69,"#eta #rightarrow #gamma#gamma");
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
            labelMassEnergy->SetText(0.13,0.775,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
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
        canvasMassWidthEta->Print(Form("%s/Eta_MassAndWidth_%s%s%s.%s",outputDirSupportPaper.Data(),centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));

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
                legendMass->AddEntry(graphEtaMass[cent][meth],centArray2[cent],"p");
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
    SetStyleHistoTH2ForGraphs(histo2DYieldPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.68);
    histo2DYieldPi0->GetXaxis()->SetMoreLogLabels();
    histo2DYieldPi0->GetXaxis()->SetLabelOffset(-0.01);

    TLatex *labelEnergyXSectionPi0  = new TLatex(0.935,0.92,collisionSystempPb.Data());
    SetStyleTLatex( labelEnergyXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);
    TLatex *labelALICEXSectionPi0  = new TLatex(0.935,0.882,textALICE.Data());
    SetStyleTLatex( labelALICEXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);
    TLatex *labelDetSysXSectionPi0  = new TLatex(0.935,0.844,"#pi^{0} #rightarrow #gamma#gamma");
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
            labelEnergyXSectionPi0->SetText(0.94,0.92,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelEnergyXSectionPi0->Draw();
            labelALICEXSectionPi0->Draw();
            labelDetSysXSectionPi0->Draw();


            TLegend* legendXSectionPi0      = GetAndSetLegend2(0.20, 0.13,0.40,0.13+(0.035*nTotMeasPi0[cent]), 0.035, 1, "", 42, 0);
            for (Int_t meth = 0; meth < 11; meth++){
                if (graphIndPi0InvYieldSys[cent][meth]) legendXSectionPi0->AddEntry(graphIndPi0InvYieldSys[cent][meth],nameMeasGlobalLabel[meth],"fp");
            }
            legendXSectionPi0->Draw();

        canvasInvYieldPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));

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
            labelALICEXSectionPi0->Draw();
            labelDetSysXSectionPi0->Draw();

            legendXSectionPi0->Draw();

        canvasInvYieldPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems_Comb_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(),  runArray[cent].Data(), suffix.Data()));

        Int_t nLinesLegendTheo              = 3;
        if (centArray[cent].BeginsWith("0-100%"))
            nLinesLegendTheo                = 4;

        TLegend* legendXSectionPi0Theo      = GetAndSetLegend2(0.20, 0.12,0.40,0.12+(0.035*nLinesLegendTheo), 0.035, 1, "", 42, 0);


        histo2DYieldPi0->Draw("copy");
            if (graphCombPi0InvYieldSys[cent]){
                DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphCombPi0InvYieldSys[cent]->Draw("E2same");
                legendXSectionPi0Theo->AddEntry(graphCombPi0InvYieldSys[cent],"Data","fp");
            }
            if (graphCombPi0InvYieldStat[cent]){
                DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphCombPi0InvYieldStat[cent]->Draw("p,same,z");
            }

            if (histoDPMJetPi0[cent]){
                SetStyleHisto(histoDPMJetPi0[cent], widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );
                histoDPMJetPi0[cent]->GetXaxis()->SetRangeUser(minPtPi0[cent],xPtLimitsPi0[cent][maxNBinsPi0[cent]]);
                histoDPMJetPi0[cent]->Draw("same,hist,l");
                legendXSectionPi0Theo->AddEntry(histoDPMJetPi0[cent],"DPMJet","l");
            }
            if (histoHIJINGPi0[cent]){
                SetStyleHisto(histoHIJINGPi0[cent], widthCommonFit*1.5, styleLineHIJING, colorHIJING );
                histoHIJINGPi0[cent]->GetXaxis()->SetRangeUser(minPtPi0[cent],xPtLimitsPi0[cent][maxNBinsPi0[cent]]);
                histoHIJINGPi0[cent]->Draw("same,hist,l");
                legendXSectionPi0Theo->AddEntry(histoHIJINGPi0[cent],"HIJING","l");
            }
            if (histoEPOSLHCPi0[cent]){
                SetStyleHisto(histoEPOSLHCPi0[cent], widthCommonFit*1.5,  styleLineEPOS3, colorEPOS3 );
                histoEPOSLHCPi0[cent]->GetXaxis()->SetRangeUser(minPtPi0[cent],xPtLimitsPi0[cent][maxNBinsPi0[cent]]);
                histoEPOSLHCPi0[cent]->Draw("same,hist,l");
                legendXSectionPi0Theo->AddEntry(histoEPOSLHCPi0[cent],"EPOS-LHC","l");
            }

            labelEnergyXSectionPi0->Draw();
            labelALICEXSectionPi0->Draw();
            labelDetSysXSectionPi0->Draw();

            legendXSectionPi0Theo->Draw();

        canvasInvYieldPi0->SaveAs(Form("%s/Pi0_InvYield_CombWithTheory_%s%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));
    }

    for (Int_t est = 0; est < 3; est++){

        Int_t nCurrEst = 0;
        for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
            if ( (centArray2[cent].Contains(nameCentEst[est].Data()) || centArray2[cent].CompareTo("0-100%") == 0) && cent != 5  && enableCentComb[cent])
                nCurrEst++;
        }

        histo2DYieldPi0->Draw("copy");
        TLegend* legendPi0YieldPaper    = GetAndSetLegend2(0.19, 0.10, 0.7, 0.10+0.035*(nCurrEst), textSizeLabelsPixel, 2, "", 43, 0.25);

            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                //         for (Int_t cent= 0; cent < maxCentRun1+maxCentRun2; cent++ ){
                if (!( centArray2[cent].Contains(nameCentEst[est].Data()) || centArray2[cent].CompareTo("0-100%") == 0) ) continue;
                if (!enableCentComb[cent] || cent == 1) continue;
                DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphCombPi0InvYieldSys[cent]->Draw("E2same");
                DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent]);
                graphCombPi0InvYieldStatWOXErr[cent]->Draw("p,same,z");
                DrawGammaSetMarkerTF1( fitTCMInvYieldPi0[cent], 7, 2, colorCentMC[cent]);
                fitTCMInvYieldPi0[cent]->Draw("same");
                legendPi0YieldPaper->AddEntry(graphCombPi0InvYieldSys[cent],centArray[cent].Data(),"pf");
                legendPi0YieldPaper->AddEntry(fitTCMInvYieldPi0[cent],"TCM fit","l");

            }

            if (enableCentComb[0] ){
                graphCombPi0InvYieldSys[0]->Draw("E2same");
                graphCombPi0InvYieldStatWOXErr[0]->Draw("p,same,z");
                fitTCMInvYieldPi0[0]->Draw("same");
            }

            TLatex *labelEnergyYieldPaper= new TLatex(0.935, 0.92, Form("%s %s", nameCentEst[est].Data(), collisionSystempPb.Data()));
            SetStyleTLatex( labelEnergyYieldPaper, 0.035, 4, 1, 42, kTRUE, 31);
            labelEnergyYieldPaper->Draw();
            TLatex *labelALICEYieldPaper= new TLatex(0.935,0.882,textALICE.Data());
            SetStyleTLatex( labelALICEYieldPaper, 0.035, 4, 1, 42, kTRUE, 31);
            labelALICEYieldPaper->Draw();
            TLatex *labelDetSysYieldPaper= new TLatex(0.935,0.844,"#pi^{0} #rightarrow #gamma#gamma");
            SetStyleTLatex( labelDetSysYieldPaper, 0.035, 4, 1, 42, kTRUE, 31);
            labelDetSysYieldPaper->Draw();
            legendPi0YieldPaper->Draw();
            canvasInvYieldPi0->SaveAs(Form("%s/Pi0_InvYieldAllCent_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));

        histo2DYieldPi0->Draw("copy");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0OlderMBSys, markerStyleCent[0], markerSizeCent[0]*0.75, colorCent[0], colorCent[0], widthLinesBoxes, kTRUE);
            graphCombPi0OlderMBSys->Draw("E2same");

            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                //         for (Int_t cent= 0; cent < maxCentRun1+maxCentRun2; cent++ ){
                if (!( centArray2[cent].Contains(nameCentEst[est].Data())  )) continue;
                if (!enableCentComb[cent] ) continue;
                graphCombPi0InvYieldSys[cent]->Draw("E2same");
                graphCombPi0InvYieldStatWOXErr[cent]->Draw("p,same,z");
                DrawGammaSetMarkerTF1( fitTCMInvYieldPi0[cent], 7, 2, colorCentMC[cent]);
                fitTCMInvYieldPi0[cent]->Draw("same");
            }

            ProduceGraphAsymmWithoutXErrors(graphCombPi0OlderMBStat);
            DrawGammaSetMarkerTGraphAsym(graphCombPi0OlderMBStat, markerStyleCent[0], markerSizeCent[0]*0.75, colorCent[0], colorCent[0]);
            graphCombPi0OlderMBStat->Draw("p,same,z");
            DrawGammaSetMarkerTF1( fitCombPi0OlderMB, 7, 2, colorCentMC[0]);
            fitCombPi0OlderMB->Draw("same");

            labelEnergyYieldPaper->Draw();
            labelALICEYieldPaper->Draw();
            labelDetSysYieldPaper->Draw();
            legendPi0YieldPaper->Draw();
            canvasInvYieldPi0->SaveAs(Form("%s/Pi0_InvYieldAllCent_OldMB_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));

        delete legendPi0YieldPaper;
    }

    // **********************************************************************************************************************
    // ******************************** Cross section for eta single measurement 2.76TeV ************************************
    // **********************************************************************************************************************
    TCanvas* canvasInvYieldEta          = new TCanvas("canvasInvYieldEta","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasInvYieldEta, 0.16, 0.02, 0.02, 0.08);
    canvasInvYieldEta->SetLogx();
    canvasInvYieldEta->SetLogy();

    TH2F * histo2DYieldEta              = new TH2F("histo2DYieldEta","histo2DYieldEta",11000,minPtEtaPlotting, maxPtEtaPlotting,1000,6e-11,3e-0);
    SetStyleHistoTH2ForGraphs(histo2DYieldEta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.68);
    histo2DYieldEta->GetXaxis()->SetMoreLogLabels();
    histo2DYieldEta->GetXaxis()->SetLabelOffset(-0.01);

    TLatex *labelEnergyXSectionEta  = new TLatex(0.935,0.92,collisionSystempPb.Data());
    SetStyleTLatex( labelEnergyXSectionEta, 0.035,4, 1, 42, kTRUE, 31);
    TLatex *labelALICEXSectionEta   = new TLatex(0.935,0.882,textALICE.Data());
    SetStyleTLatex( labelALICEXSectionEta, 0.035,4, 1, 42, kTRUE, 31);
    TLatex *labelDetSysXSectionEta  = new TLatex(0.935,0.844,"#eta #rightarrow #gamma#gamma");
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
        labelEnergyXSectionEta->SetText(0.94,0.92,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
        labelEnergyXSectionEta->Draw();
        labelALICEXSectionEta->Draw();
        labelDetSysXSectionEta->Draw();


        TLegend* legendXSectionEta      = GetAndSetLegend2(0.20, 0.13,0.40,0.13+(0.035*nTotMeasEta[cent]), 0.035, 1, "", 42, 0);
        for (Int_t meth = 0; meth < 11; meth++){
            if (graphIndEtaInvYieldSys[cent][meth]) legendXSectionEta->AddEntry(graphIndEtaInvYieldSys[cent][meth],nameMeasGlobalLabel[meth],"fp");
        }
        legendXSectionEta->Draw();

        canvasInvYieldEta->SaveAs(Form("%s/Eta_InvYieldCompAllSystems_%s%s%s.%s",outputDirSupportPaper.Data(),centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));

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
            labelALICEXSectionEta->Draw();
            labelDetSysXSectionEta->Draw();
            legendXSectionEta->Draw();

        canvasInvYieldEta->SaveAs(Form("%s/Eta_InvYieldCompAllSystems_Comb_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(),
                                       runArray[cent].Data(), suffix.Data()));

        Int_t nLinesLegendTheo              = 3;
        if (centArray[cent].BeginsWith("0-100%"))
            nLinesLegendTheo                = 4;

        TLegend* legendXSectionEtaTheo      = GetAndSetLegend2(0.20, 0.12,0.40,0.12+(0.035*nLinesLegendTheo), 0.035, 1, "", 42, 0);

        histo2DYieldEta->Draw("copy");
            if (graphCombEtaInvYieldSys[cent]){
                DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphCombEtaInvYieldSys[cent]->Draw("E2same");
                legendXSectionEtaTheo->AddEntry(graphCombEtaInvYieldSys[cent],"Data","fp");
            }
            if (graphCombEtaInvYieldStat[cent]){
                DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStat[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphCombEtaInvYieldStat[cent]->Draw("p,same,z");
            }

            if (histoDPMJetEta[cent]){
                SetStyleHisto(histoDPMJetEta[cent], widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );
                histoDPMJetEta[cent]->GetXaxis()->SetRangeUser(minPtEta[cent],xPtLimitsEta[cent][maxNBinsEta[cent]]);
                histoDPMJetEta[cent]->Draw("same,hist,l");
                legendXSectionEtaTheo->AddEntry(histoDPMJetEta[cent],"DPMJet","l");
            }
            if (histoHIJINGEta[cent]){
                SetStyleHisto(histoHIJINGEta[cent], widthCommonFit*1.5, styleLineHIJING, colorHIJING );
                histoHIJINGEta[cent]->GetXaxis()->SetRangeUser(minPtEta[cent],xPtLimitsEta[cent][maxNBinsEta[cent]]);
                histoHIJINGEta[cent]->Draw("same,hist,l");
                legendXSectionEtaTheo->AddEntry(histoHIJINGEta[cent],"HIJING","l");
            }
            if (histoEPOSLHCEta[cent]){
                SetStyleHisto(histoEPOSLHCEta[cent], widthCommonFit*1.5,  styleLineEPOS3, colorEPOS3 );
                histoEPOSLHCEta[cent]->GetXaxis()->SetRangeUser(minPtEta[cent],xPtLimitsEta[cent][maxNBinsEta[cent]]);
                histoEPOSLHCEta[cent]->Draw("same,hist,l");
                legendXSectionEtaTheo->AddEntry(histoEPOSLHCEta[cent],"EPOS-LHC","l");
            }

            labelEnergyXSectionEta->Draw();
            labelDetSysXSectionEta->Draw();
            labelALICEXSectionEta->Draw();

            legendXSectionEtaTheo->Draw();

        canvasInvYieldEta->SaveAs(Form("%s/Eta_InvYield_CombWithTheory_%s%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));
    }
    canvasInvYieldEta->cd();

    for (Int_t est = 0; est < 3; est++){

        Int_t nCurrEst = 0;
        for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
            if ( (centArray2[cent].Contains(nameCentEst[est].Data()) || centArray2[cent].CompareTo("0-100%") == 0) && cent != 5  && enableCentComb[cent])
                nCurrEst++;
        }

        histo2DYieldEta->Draw("copy");
        TLegend* legendEtaYieldPaper    = GetAndSetLegend2(0.19, 0.10, 0.7, 0.10+0.035*(nCurrEst), textSizeLabelsPixel, 2, "", 43, 0.25);

        for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
            //         for (Int_t cent= 0; cent < maxCentRun1+maxCentRun2; cent++ ){
            if (!( centArray2[cent].Contains(nameCentEst[est].Data()) || centArray2[cent].CompareTo("0-100%") == 0) ) continue;
            if (!enableCentComb[cent] || cent == 1) continue;
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
            graphCombEtaInvYieldSys[cent]->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent]);
            graphCombEtaInvYieldStatWOXErr[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTF1( fitTCMInvYieldEta[cent], 7, 2, colorCentMC[cent]);
            fitTCMInvYieldEta[cent]->Draw("same");
            legendEtaYieldPaper->AddEntry(graphCombEtaInvYieldSys[cent],centArray[cent].Data(),"pf");
            legendEtaYieldPaper->AddEntry(fitTCMInvYieldEta[cent],"TCM fit","l");

        }

        if (enableCentComb[0] ){
            graphCombEtaInvYieldSys[0]->Draw("E2same");
            graphCombEtaInvYieldStatWOXErr[0]->Draw("p,same,z");
            fitTCMInvYieldEta[0]->Draw("same");
        }

        TLatex *labelEnergyYieldPaper= new TLatex(0.935, 0.92, Form("%s %s", nameCentEst[est].Data(), collisionSystempPb.Data()));
        SetStyleTLatex( labelEnergyYieldPaper, 0.035, 4, 1, 42, kTRUE, 31);
        labelEnergyYieldPaper->Draw();
        TLatex *labelALICEYieldPaper= new TLatex(0.935,0.882,textALICE.Data());
        SetStyleTLatex( labelALICEYieldPaper, 0.035, 4, 1, 42, kTRUE, 31);
        labelALICEYieldPaper->Draw();
        TLatex *labelDetSysYieldPaper= new TLatex(0.935,0.844,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysYieldPaper, 0.035, 4, 1, 42, kTRUE, 31);
        labelDetSysYieldPaper->Draw();
        legendEtaYieldPaper->Draw();
        canvasInvYieldEta->SaveAs(Form("%s/Eta_InvYieldAllCent_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));

        histo2DYieldEta->Draw("copy");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaOlderMBSys, markerStyleCent[0], markerSizeCent[0]*0.75, colorCent[0], colorCent[0], widthLinesBoxes, kTRUE);
            graphCombEtaOlderMBSys->Draw("E2same");

            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                //         for (Int_t cent= 0; cent < maxCentRun1+maxCentRun2; cent++ ){
                if (!( centArray2[cent].Contains(nameCentEst[est].Data())  )) continue;
                if (!enableCentComb[cent] ) continue;
                graphCombEtaInvYieldSys[cent]->Draw("E2same");
                graphCombEtaInvYieldStatWOXErr[cent]->Draw("p,same,z");
                DrawGammaSetMarkerTF1( fitTCMInvYieldEta[cent], 7, 2, colorCentMC[cent]);
                fitTCMInvYieldEta[cent]->Draw("same");
            }

            ProduceGraphAsymmWithoutXErrors(graphCombEtaOlderMBStat);
            DrawGammaSetMarkerTGraphAsym(graphCombEtaOlderMBStat, markerStyleCent[0], markerSizeCent[0]*0.75, colorCent[0], colorCent[0]);
            graphCombEtaOlderMBStat->Draw("p,same,z");
            DrawGammaSetMarkerTF1( fitCombEtaOlderMB, 7, 2, colorCentMC[0]);
            fitCombEtaOlderMB->Draw("same");

            labelEnergyYieldPaper->Draw();
            labelALICEYieldPaper->Draw();
            labelDetSysYieldPaper->Draw();
            legendEtaYieldPaper->Draw();
            canvasInvYieldEta->SaveAs(Form("%s/Eta_InvYieldAllCent_OldMB_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));

        delete legendEtaYieldPaper;
    }

    // **********************************************************************************************************************
    // ******************************** Cross section for pi0 & eta single measurement **************************************
    // **********************************************************************************************************************
    TH2F * histo2DYieldPi0Eta              = new TH2F("histo2DYieldPi0Eta","histo2DYieldPi0Eta",11000,minPtPi0Plotting, maxPtPi0Plotting,1000,6e-13,3e1);
    SetStyleHistoTH2ForGraphs(histo2DYieldPi0Eta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.68);
    histo2DYieldPi0Eta->GetXaxis()->SetMoreLogLabels();
    histo2DYieldPi0Eta->GetXaxis()->SetLabelOffset(-0.01);

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCentComb[cent]) continue;
        canvasInvYieldPi0->cd();
        histo2DYieldPi0Eta->Draw("copy");

        histo2DYieldPi0Eta->Draw("copy");
        if (graphCombPi0InvYieldSys[cent]){
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
            graphCombPi0InvYieldSys[cent]->Draw("E2same");
        }
        if (graphCombPi0InvYieldStat[cent]){
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
            graphCombPi0InvYieldStat[cent]->Draw("p,same,z");
        }

        if (graphCombEtaInvYieldSys[cent]){
            TGraphAsymmErrors* graphCombEtaInvYieldSysScaled = ScaleGraph( graphCombEtaInvYieldSys[cent], 0.01);
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSysScaled, markerStyleCentMC[cent], markerSizeCentMC[cent]*0.75, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys[cent], markerStyleCentMC[cent], markerSizeCentMC[cent]*0.75, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
            graphCombEtaInvYieldSysScaled->Draw("E2same");
        }
        if (graphCombEtaInvYieldStat[cent]){
            TGraphAsymmErrors* graphCombEtaInvYieldStatScaled = ScaleGraph( graphCombEtaInvYieldStat[cent], 0.01);
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStatScaled, markerStyleCentMC[cent], markerSizeCentMC[cent]*0.75, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
            graphCombEtaInvYieldStatScaled->Draw("p,same,z");
        }

        Int_t nTheory = 0;
        if (histoDPMJetPi0[cent]){
            SetStyleHisto(histoDPMJetPi0[cent], widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );
            histoDPMJetPi0[cent]->GetXaxis()->SetRangeUser(minPtPi0[cent],xPtLimitsPi0[cent][maxNBinsPi0[cent]]);
            histoDPMJetPi0[cent]->Draw("same,hist,l");

            histoDPMJetEta[cent]->Scale(0.01);
            SetStyleHisto(histoDPMJetEta[cent], widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );
            histoDPMJetEta[cent]->GetXaxis()->SetRangeUser(minPtEta[cent],xPtLimitsEta[cent][maxNBinsEta[cent]]);
            histoDPMJetEta[cent]->Draw("same,hist,l");
            nTheory++;
        }
        if (histoHIJINGPi0[cent]){
            SetStyleHisto(histoHIJINGPi0[cent], widthCommonFit*1.5, styleLineHIJING, colorHIJING );
            histoHIJINGPi0[cent]->GetXaxis()->SetRangeUser(minPtPi0[cent],xPtLimitsPi0[cent][maxNBinsPi0[cent]]);
            histoHIJINGPi0[cent]->Draw("same,hist,l");
            histoHIJINGEta[cent]->Scale(0.01);
            SetStyleHisto(histoHIJINGEta[cent], widthCommonFit*1.5, styleLineHIJING, colorHIJING );
            histoHIJINGEta[cent]->GetXaxis()->SetRangeUser(minPtEta[cent],xPtLimitsEta[cent][maxNBinsEta[cent]]);
            histoHIJINGEta[cent]->Draw("same,hist,l");
            nTheory++;
        }
        if (histoEPOSLHCPi0[cent]){
            SetStyleHisto(histoEPOSLHCPi0[cent], widthCommonFit*1.5,  styleLineEPOS3, colorEPOS3 );
            histoEPOSLHCPi0[cent]->GetXaxis()->SetRangeUser(minPtPi0[cent],xPtLimitsPi0[cent][maxNBinsPi0[cent]]);
            histoEPOSLHCPi0[cent]->Draw("same,hist,l");
            histoEPOSLHCEta[cent]->Scale(0.01);
            SetStyleHisto(histoEPOSLHCEta[cent], widthCommonFit*1.5,  styleLineEPOS3, colorEPOS3 );
            histoEPOSLHCEta[cent]->GetXaxis()->SetRangeUser(minPtEta[cent],xPtLimitsEta[cent][maxNBinsEta[cent]]);
            histoEPOSLHCEta[cent]->Draw("same,hist,l");
//
            nTheory++;
        }

        TLegend* legendXSectionPi0Eta      = GetAndSetLegend2(0.20, 0.12,0.40,0.12+(0.035*(2+nTheory)), 0.035, 1, "", 42, 0);
            if (graphCombPi0InvYieldSys[cent]) legendXSectionPi0Eta->AddEntry(graphCombPi0InvYieldSys[cent],"#pi^{0}","fp");
            if (graphCombEtaInvYieldSys[cent]) legendXSectionPi0Eta->AddEntry(graphCombEtaInvYieldSys[cent],"#eta #times 10^{-2}","fp");
            if (histoDPMJetPi0[cent]) legendXSectionPi0Eta->AddEntry(histoDPMJetPi0[cent],"DPMJet","l");
            if (histoHIJINGPi0[cent]) legendXSectionPi0Eta->AddEntry(histoHIJINGPi0[cent],"HIJING","l");
            if (histoEPOSLHCPi0[cent]) legendXSectionPi0Eta->AddEntry(histoEPOSLHCPi0[cent],"EPOS-LHC","l");
        legendXSectionPi0Eta->Draw();
        labelEnergyXSectionPi0->SetText(0.94,0.92,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
        labelEnergyXSectionPi0->Draw();
        labelALICEXSectionPi0->Draw();

        canvasInvYieldPi0->SaveAs(Form("%s/Pi0AndEta_InvYield_CombWithTheory_%s%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));

        if (histoHIJINGEta[cent])histoHIJINGEta[cent]->Scale(100);
        if (histoEPOSLHCEta[cent])histoEPOSLHCEta[cent]->Scale(100);
        if (histoDPMJetEta[cent])histoDPMJetEta[cent]->Scale(100);
    }



    // **********************************************************************************************************************
    // ******************************** Acceptance * Efficiency for pi0 single measurement **********************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 55;
    textSizeLabelsRel                   = ((Double_t)textSizeLabelsPixel)/1200;
//     cout << textSizeLabelsRel << endl;

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
            labelEnergyEffi->SetText(0.137,0.87,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
            labelEnergyEffi->Draw();
            labelDetSysEffiPi0->Draw();

        canvasAcceptanceTimesEff->Update();
        canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff_%s%s%s.%s",outputDirSupportPaper.Data(),centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));
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
                legendEffiAccPi0->AddEntry(graphPi0EffTimesAcc[cent][meth],centArray2[cent],"p");
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
        canvasAcceptanceTimesEff->Print(Form("%s/Eta_AcceptanceTimesEff_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));
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
                legendEffiAccEta->AddEntry(graphEtaEffTimesAcc[cent][meth],centArray2[cent],"p");
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
                    canvasEffectiveSecCorr->Print(Form("%s/Pi0_EffectiveSecCorr_%s%s%s_%s.%s",outputDirSupportPaper.Data(), nameSecPi0SourceRead[k].Data(),
                                                       centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));
                }
            }
        }
    delete canvasEffectiveSecCorr;

    fileFitsOutput << "*******************************************************************************************" << endl;
    fileFitsOutput << "****************************** Power law fit pi0 ******************************************" << endl;
    fileFitsOutput << "*******************************************************************************************" << endl;
    TF1* fitPowInvYieldPi0Tot[17]   = { NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL,
                                        NULL, NULL, NULL, NULL, NULL,   NULL, NULL };
    TF1* fitPowInvYieldPi0Stat[17]   = { NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,
                                        NULL, NULL, NULL, NULL, NULL,   NULL, NULL };
    TF1* fitOHagInvYieldPi0Tot[17]   = { NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,
                                        NULL, NULL, NULL, NULL, NULL,   NULL, NULL };
    TF1* fitPowInvYieldEtaTot[17]   = { NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL, NULL,
                                        NULL, NULL, NULL, NULL, NULL,   NULL, NULL };
    TF1* fitPowInvYieldEtaStat[17]   = { NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,
                                        NULL, NULL, NULL, NULL, NULL,   NULL, NULL };
    TF1* fitOHagInvYieldEtaTot[17]   = { NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,
                                        NULL, NULL, NULL, NULL, NULL,   NULL, NULL };

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

        fitPowInvYieldEtaTot[cent]   = FitObject("powPure","fitPowInvYieldEtaTot","Eta",graphCombEtaInvYieldTot[cent],4,20. ,NULL,"QNRMEX0+","", kFALSE);
        cout << WriteParameterToFile(fitPowInvYieldEtaTot[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldEtaTot[cent])<< endl;
        fitPowInvYieldEtaStat[cent]   = FitObject("powPure","fitPowInvYieldEta","Eta",graphCombEtaInvYieldStat[cent],4,20. ,NULL,"QNRMEX0+","", kFALSE);
        cout << WriteParameterToFile(fitPowInvYieldEtaStat[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldEtaStat[cent])<< endl;
        fitOHagInvYieldEtaTot[cent]   = FitObject("oHag","fitOHagInvYieldEta","Eta",graphCombEtaInvYieldTot[cent],0.3,20. ,NULL,"QNRMEX0+","", kFALSE);
        cout << WriteParameterToFile(fitOHagInvYieldEtaTot[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitOHagInvYieldEtaTot[cent])<< endl;
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
    Double_t arrayBoundariesY1_XSec[7];
    Double_t relativeMarginsXXSec[3];
    Double_t relativeMarginsYXSec[3];
    textSizeLabelsPixel = 48;
    ReturnCorrectValuesForCanvasScaling(1250,2400, 1, 6,0.15, 0.005, 0.003,0.05,arrayBoundariesX1_XSec,arrayBoundariesY1_XSec,relativeMarginsXXSec,relativeMarginsYXSec);

    TCanvas* canvasInvSectionPaper      = new TCanvas("canvasInvSectionPaper","",0,0,1250,2400);  // gives the page size
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

    TPad* padInvSectionMBRatio         = new TPad("padInvSectionMBRatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[4], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[3],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionMBRatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[1]);
    padInvSectionMBRatio->Draw();
    Double_t textsizeLabelsXSecMiddle   = 0;
    Double_t textsizeFacXSecMiddle      = 0;
    if (padInvSectionMBRatio->XtoPixel(padInvSectionMBRatio->GetX2()) < padInvSectionMBRatio->YtoPixel(padInvSectionMBRatio->GetY1())){
        textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionMBRatio->XtoPixel(padInvSectionMBRatio->GetX2()) ;
        textsizeFacXSecMiddle           = (Double_t)1./padInvSectionMBRatio->XtoPixel(padInvSectionMBRatio->GetX2()) ;
    } else {
        textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionMBRatio->YtoPixel(padInvSectionMBRatio->GetY1());
        textsizeFacXSecMiddle           = (Double_t)1./padInvSectionMBRatio->YtoPixel(padInvSectionMBRatio->GetY1());
    }

    TPad* padInvSectionCentRatio1      = new TPad("padInvSectionCentRatio1", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[5], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[4],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionCentRatio1, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[1]);
    padInvSectionCentRatio1->Draw();
    Double_t textsizeLabelsXSecDown     = 0;
    Double_t textsizeFacXSecDown        = 0;
    if (padInvSectionCentRatio1->XtoPixel(padInvSectionCentRatio1->GetX2()) < padInvSectionCentRatio1->YtoPixel(padInvSectionCentRatio1->GetY1())){
        textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionCentRatio1->XtoPixel(padInvSectionCentRatio1->GetX2()) ;
        textsizeFacXSecDown             = (Double_t)1./padInvSectionCentRatio1->XtoPixel(padInvSectionCentRatio1->GetX2()) ;
    } else {
        textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionCentRatio1->YtoPixel(padInvSectionCentRatio1->GetY1());
        textsizeFacXSecDown             = (Double_t)1./padInvSectionCentRatio1->YtoPixel(padInvSectionCentRatio1->GetY1());
    }

    TPad* padInvSectionCentRatio2      = new TPad("padInvSectionCentRatio2", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[6], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[5],-1, -1, -2);
    DrawGammaPadSettings( padInvSectionCentRatio2, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[2]);
    padInvSectionCentRatio2->Draw();
    Double_t textsizeLabelsXSecDown2    = 0;
    Double_t textsizeFacXSecDown2       = 0;
    if (padInvSectionCentRatio2->XtoPixel(padInvSectionCentRatio2->GetX2()) < padInvSectionCentRatio2->YtoPixel(padInvSectionCentRatio2->GetY1())){
        textsizeLabelsXSecDown2         = (Double_t)textSizeLabelsPixel/padInvSectionCentRatio2->XtoPixel(padInvSectionCentRatio2->GetX2()) ;
        textsizeFacXSecDown2            = (Double_t)1./padInvSectionCentRatio2->XtoPixel(padInvSectionCentRatio2->GetX2()) ;
    } else {
        textsizeLabelsXSecDown2         = (Double_t)textSizeLabelsPixel/padInvSectionCentRatio2->YtoPixel(padInvSectionCentRatio2->GetY1());
        textsizeFacXSecDown2            = (Double_t)1./padInvSectionCentRatio2->YtoPixel(padInvSectionCentRatio2->GetY1());
    }

//     TLegend* legendXsectionPaper    = GetAndSetLegend2(0.19, 0.03, 0.7, 0.03+0.05*(maxCentRun1+maxCentRun2), textSizeLabelsPixel, 2, "", 43, 0.25);
    TLegend* legendXsectionPaper    = GetAndSetLegend2(0.19, 0.03, 0.7, 0.03+0.05*(maxCentRun1+1), textSizeLabelsPixel, 1, "", 43, 0.15);
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

        for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
            if (!enableCentComb[cent] || cent == 1) continue;
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
            graphCombPi0InvYieldSys[cent]->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent]);
            graphCombPi0InvYieldStatWOXErr[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTF1( fitTCMInvYieldPi0[cent], 7, 2, colorCentMC[cent]);
            fitTCMInvYieldPi0[cent]->Draw("same");
            legendXsectionPaper->AddEntry(graphCombPi0InvYieldSys[cent],centArray2[cent].Data(),"pf");
//             legendXsectionPaper->AddEntry(fitTCMInvYieldPi0[cent],"TCM fit","l");
        }

        DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombPi0Sys[0], 20, markerSizeCent[0]*0.75, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphPPInvYieldCombPi0Sys[0]->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombPi0StatW0XErr[0], 20, markerSizeCent[0]*0.75, kBlack, kBlack);
        graphPPInvYieldCombPi0StatW0XErr[0]->Draw("p,same,z");
        DrawGammaSetMarkerTF1( fitPPTCMInvYieldPi0[0], 7, 2, kGray+2);
        fitPPTCMInvYieldPi0[0]->Draw("same");
        legendXsectionPaper->AddEntry(graphPPInvYieldCombPi0Sys[0],collisionSystempp,"pf");
        legendXsectionPaper->AddEntry(fitPPTCMInvYieldPi0[0],"TCM fits","l");

        TLatex *labelEnergyXSectionPaper= new TLatex(0.935, 0.91, collisionSystempPb.Data());
        SetStyleTLatex( labelEnergyXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
        labelEnergyXSectionPaper->Draw();
        TLatex *labelALICEXSectionPaper= new TLatex(0.935,0.867,textALICE.Data());
        SetStyleTLatex( labelALICEXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
        labelALICEXSectionPaper->Draw();
        TLatex *labelDetSysXSectionPaper= new TLatex(0.935,0.83,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
        labelDetSysXSectionPaper->Draw();
        legendXsectionPaper->Draw();

    padInvSectionMBRatio->cd();
    padInvSectionMBRatio->SetLogx(1);
        TH2F * ratio2DData1               = new TH2F("ratio2DData1","ratio2DData1",1000,minPtPi0Plotting, maxPtPi0Plotting,1000,0.65,1.55);
        SetStyleHistoTH2ForGraphs(ratio2DData1, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                  0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.257/(textsizeFacXSecMiddle*marginXSec), 510, 508);
        ratio2DData1->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DData1->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DData1->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DData1->GetXaxis()->SetLabelFont(42);
        ratio2DData1->GetYaxis()->SetLabelFont(42);
        ratio2DData1->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DData1->GetXaxis()->SetTickLength(0.07);
        ratio2DData1->DrawCopy();

        boxErrorNormRatioPi0->Draw();

//             if (histoRatioPi0EPOSLHCToFit[cent]){
//                 SetStyleHisto(histoRatioPi0EPOSLHCToFit[cent], widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );
//                 histoRatioPi0EPOSLHCToFit[cent]->Draw("same,hist,l");
//             }
        if (graphRatioPi0CombCombFitStatWOXErr[intRefMB]){
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr[intRefMB], markerStyleCent[intRefMB], markerSizeCent[intRefMB], colorCent[intRefMB], colorCent[intRefMB], widthLinesBoxes, kFALSE);
            graphRatioPi0CombCombFitStatWOXErr[intRefMB]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys[intRefMB], markerStyleCent[intRefMB], markerSizeCent[intRefMB], colorCent[intRefMB], colorCent[intRefMB], widthLinesBoxes, kTRUE, 0);
            graphRatioPi0CombCombFitSys[intRefMB]->SetLineWidth(0);
            graphRatioPi0CombCombFitSys[intRefMB]->Draw("2,same");
            graphRatioPi0CombCombFitStatWOXErr[intRefMB]->Draw("p,same");
        }
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting,1., 1.,0.1,kGray+1);

    padInvSectionCentRatio1->cd();
    padInvSectionCentRatio1->SetLogx(1);
        TH2F * ratioData2            = new TH2F("ratioData2","ratioData2",1000,minPtPi0Plotting, maxPtPi0Plotting,1000,0.65,1.55);
        SetStyleHistoTH2ForGraphs(ratioData2, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                                  0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.257/(textsizeFacXSecDown*marginXSec), 510, 508);
        ratioData2->GetYaxis()->SetNoExponent(kTRUE);
        ratioData2->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratioData2->GetXaxis()->SetNoExponent(kTRUE);
        ratioData2->GetXaxis()->SetLabelFont(42);
        ratioData2->GetYaxis()->SetLabelFont(42);
        ratioData2->GetYaxis()->SetLabelOffset(+0.01);
        ratioData2->GetXaxis()->SetTickLength(0.06);
        ratioData2->GetYaxis()->SetTickLength(0.04);
        ratioData2->DrawCopy();

        for (Int_t cent= 0; cent < 2; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
            graphRatioPi0CombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioPi0CombCombFitSys[cent]->SetLineWidth(0);
            graphRatioPi0CombCombFitSys[cent]->Draw("2,same");
            graphRatioPi0CombCombFitStatWOXErr[cent]->Draw("p,same");
        }
        boxErrorNormRatioPi0->Draw();
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting,1., 1.,0.1,kGray+1);

    padInvSectionCentRatio2->cd();
    padInvSectionCentRatio2->SetLogx(1);
        TH2F * ratioData3            = new TH2F("ratioData3","ratioData3",1000,minPtPi0Plotting, maxPtPi0Plotting,1000,0.65,1.55);
        SetStyleHistoTH2ForGraphs(ratioData3, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}", 0.85*textsizeLabelsXSecDown2, textsizeLabelsXSecDown2,
                                0.85*textsizeLabelsXSecDown2,textsizeLabelsXSecDown2, 0.9,0.257/(textsizeFacXSecDown2*marginXSec), 510, 508);
        ratioData3->GetYaxis()->SetNoExponent(kTRUE);
        ratioData3->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratioData3->GetXaxis()->SetNoExponent(kTRUE);
        ratioData3->GetXaxis()->SetLabelFont(42);
        ratioData3->GetYaxis()->SetLabelFont(42);
        ratioData3->GetYaxis()->SetLabelOffset(+0.01);
        ratioData3->GetXaxis()->SetTickLength(0.06);
        ratioData3->GetYaxis()->SetTickLength(0.04);
        ratioData3->DrawCopy();

        for (Int_t cent= 2; cent < 4; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
            graphRatioPi0CombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioPi0CombCombFitSys[cent]->SetLineWidth(0);
            graphRatioPi0CombCombFitSys[cent]->Draw("2,same");
            graphRatioPi0CombCombFitStatWOXErr[cent]->Draw("p,same");
        }

        boxErrorNormRatioPi0->Draw();
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting,1., 1.,0.1,kGray+1);

    canvasInvSectionPaper->Print(Form("%s/Pi0_InvYieldWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));

    //*************************************************************************************************************
    //***************************** Paper plot invyield eta with soft-physics *************************************
    //*************************************************************************************************************
    TLegend* legendXsectionPaperEta    = GetAndSetLegend2(0.19, 0.03, 0.7, 0.03+0.05*(maxCentRun1+1), textSizeLabelsPixel, 1, "", 43, 0.15);
//     TLegend* legendXsectionPaperEta    = GetAndSetLegend2(0.19, 0.03, 0.7, 0.03+0.05*(maxCentRun1+maxCentRun2), textSizeLabelsPixel, 2, "", 43, 0.25);

    canvasInvSectionPaper->cd();
    padInvSectionSpec->cd();
    padInvSectionSpec->SetLogy(1);
    padInvSectionSpec->SetLogx(1);
        SetStyleHistoTH2ForGraphs(histo2DYieldEta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                                0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.257/(textsizeFacXSecUp*marginXSec));
        histo2DYieldEta->GetXaxis()->SetMoreLogLabels();
        histo2DYieldEta->GetXaxis()->SetLabelOffset(+0.01);
        histo2DYieldEta->Draw();

        for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
            if (!enableCentComb[cent] || cent == 1 ) continue;
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
            graphCombEtaInvYieldSys[cent]->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent]);
            graphCombEtaInvYieldStatWOXErr[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTF1( fitTCMInvYieldEta[cent], 7, 2, colorCentMC[cent]);
            fitTCMInvYieldEta[cent]->Draw("same");
            legendXsectionPaperEta->AddEntry(graphCombEtaInvYieldSys[cent],centArray2[cent].Data(),"pf");
//             legendXsectionPaperEta->AddEntry(fitTCMInvYieldEta[cent],"TCM fit","l");
        }
        DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombEtaSys[0], 20, markerSizeCent[0]*0.75, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphPPInvYieldCombEtaSys[0]->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombEtaStatW0XErr[0], 20, markerSizeCent[0]*0.75, kBlack, kBlack);
        graphPPInvYieldCombEtaStatW0XErr[0]->Draw("p,same,z");
        DrawGammaSetMarkerTF1( fitPPTCMInvYieldEta[0], 7, 2, kGray+2);
        fitPPTCMInvYieldEta[0]->Draw("same");
        legendXsectionPaperEta->AddEntry(graphPPInvYieldCombEtaSys[0],collisionSystempp,"pf");
        legendXsectionPaperEta->AddEntry(fitPPTCMInvYieldEta[0],"TCM fits","l");

        labelEnergyXSectionPaper->Draw();
        labelALICEXSectionPaper->Draw();
        TLatex *labelDetSysXSectionPaperEta= new TLatex(0.935,0.83,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPaperEta, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
        labelDetSysXSectionPaperEta->Draw();
        legendXsectionPaper->Draw();


    padInvSectionMBRatio->cd();
    padInvSectionMBRatio->SetLogx(1);
        TH2F * ratio2DData1Eta               = new TH2F("ratio2DData1Eta","ratio2DData1Eta",1000, minPtEtaPlotting, maxPtEtaPlotting,1000,0.3,2.2);
        SetStyleHistoTH2ForGraphs(ratio2DData1Eta, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
                                  0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.257/(textsizeFacXSecMiddle*marginXSec), 510, 508);
        ratio2DData1Eta->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DData1Eta->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DData1Eta->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DData1Eta->GetXaxis()->SetLabelFont(42);
        ratio2DData1Eta->GetYaxis()->SetLabelFont(42);
        ratio2DData1Eta->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DData1Eta->GetXaxis()->SetTickLength(0.07);
        ratio2DData1Eta->DrawCopy();

        boxErrorNormRatioEta->Draw();
        if (graphRatioEtaCombCombFitStatWOXErr[intRefMB]){
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatWOXErr[intRefMB], markerStyleCent[intRefMB], markerSizeCent[intRefMB], colorCent[intRefMB], colorCent[intRefMB], widthLinesBoxes, kFALSE);
            graphRatioEtaCombCombFitStatWOXErr[intRefMB]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys[intRefMB], markerStyleCent[intRefMB], markerSizeCent[intRefMB], colorCent[intRefMB], colorCent[intRefMB], widthLinesBoxes, kTRUE, 0);
            graphRatioEtaCombCombFitSys[intRefMB]->SetLineWidth(0);
            graphRatioEtaCombCombFitSys[intRefMB]->Draw("2,same");
            graphRatioEtaCombCombFitStatWOXErr[intRefMB]->Draw("p,same");
        }

        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting,1., 1.,0.1,kGray);

    padInvSectionCentRatio1->cd();
    padInvSectionCentRatio1->SetLogx(1);
        TH2F * ratioData2Eta            = new TH2F("ratioData2Eta","ratioData2Eta",1000, minPtEtaPlotting, maxPtEtaPlotting,1000,0.3,2.2);
        SetStyleHistoTH2ForGraphs(ratioData2Eta, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
                                  0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.257/(textsizeFacXSecDown*marginXSec), 510, 508);
        ratioData2Eta->GetYaxis()->SetNoExponent(kTRUE);
        ratioData2Eta->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratioData2Eta->GetXaxis()->SetNoExponent(kTRUE);
        ratioData2Eta->GetXaxis()->SetLabelFont(42);
        ratioData2Eta->GetYaxis()->SetLabelFont(42);
        ratioData2Eta->GetYaxis()->SetLabelOffset(+0.01);
        ratioData2Eta->GetXaxis()->SetTickLength(0.06);
        ratioData2Eta->GetYaxis()->SetTickLength(0.04);
        ratioData2Eta->DrawCopy();

        for (Int_t cent= 0; cent < 2; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
            graphRatioEtaCombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioEtaCombCombFitSys[cent]->SetLineWidth(0);
            graphRatioEtaCombCombFitSys[cent]->Draw("2,same");
            graphRatioEtaCombCombFitStatWOXErr[cent]->Draw("p,same");
        }
        boxErrorNormRatioEta->Draw();

        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting,1., 1.,0.1,kGray);

    padInvSectionCentRatio2->cd();
    padInvSectionCentRatio2->SetLogx(1);
        TH2F * ratioData3Eta            = new TH2F("ratioData3Eta","ratioData3Eta",1000, minPtEtaPlotting, maxPtEtaPlotting,1000,0.3,2.2);
        SetStyleHistoTH2ForGraphs(ratioData3Eta, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{fit}", 0.85*textsizeLabelsXSecDown2, textsizeLabelsXSecDown2,
                                0.85*textsizeLabelsXSecDown2,textsizeLabelsXSecDown2, 0.9,0.257/(textsizeFacXSecDown2*marginXSec), 510, 508);
        ratioData3Eta->GetYaxis()->SetNoExponent(kTRUE);
        ratioData3Eta->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratioData3Eta->GetXaxis()->SetNoExponent(kTRUE);
        ratioData3Eta->GetXaxis()->SetLabelFont(42);
        ratioData3Eta->GetYaxis()->SetLabelFont(42);
        ratioData3Eta->GetYaxis()->SetLabelOffset(+0.01);
        ratioData3Eta->GetXaxis()->SetTickLength(0.06);
        ratioData3Eta->GetYaxis()->SetTickLength(0.04);
        ratioData3Eta->DrawCopy();

        for (Int_t cent= 2; cent < 4; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
            graphRatioEtaCombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioEtaCombCombFitSys[cent]->SetLineWidth(0);
            graphRatioEtaCombCombFitSys[cent]->Draw("2,same");
            graphRatioEtaCombCombFitStatWOXErr[cent]->Draw("p,same");
        }
        boxErrorNormRatioEta->Draw();
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting,1., 1.,0.1,kGray);

    canvasInvSectionPaper->Print(Form("%s/Eta_InvYieldWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));



    //*************************************************************************************************************
    //***************************** Paper plot invyield and ratios ************************************************
    //*************************************************************************************************************
    Double_t arrayBoundariesX1_XSec2[2];
    Double_t arrayBoundariesY1_XSec2[6];
    Double_t relativeMarginsXXSec2[3];
    Double_t relativeMarginsYXSec2[3];
    textSizeLabelsPixel = 48;
    ReturnCorrectValuesForCanvasScaling(1250,2000, 1, 5,0.15, 0.005, 0.003,0.05,arrayBoundariesX1_XSec2,arrayBoundariesY1_XSec2,relativeMarginsXXSec2,relativeMarginsYXSec2);

    TCanvas* canvasRatioGeneratorsPaper      = new TCanvas("canvasRatioGeneratorsPaper","",0,0,1250,2000);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioGeneratorsPaper,  0.13, 0.02, 0.03, 0.06);
    canvasRatioGeneratorsPaper->cd();

    TPad* padRatioGenerators[5]             = {NULL, NULL, NULL, NULL, NULL};
    Double_t textSizeLablesGenerators[5]    = {0};
    Double_t textsizeFacXSecGenerators[5]   = {0};
    marginXSec                              = relativeMarginsXXSec2[0]*1250;
    for (Int_t i = 0; i < 5; i++){
        padRatioGenerators[i]               = new TPad(Form("padRatioGen%i",i), "", arrayBoundariesX1_XSec2[0], arrayBoundariesY1_XSec2[i+1], arrayBoundariesX1_XSec2[1], arrayBoundariesY1_XSec2[i],-1, -1, -2);
        if (i == 0) DrawGammaPadSettings( padRatioGenerators[i], relativeMarginsXXSec2[0], relativeMarginsXXSec2[2], relativeMarginsYXSec2[0], relativeMarginsYXSec2[1]);
        else if (i == 4) DrawGammaPadSettings( padRatioGenerators[i], relativeMarginsXXSec2[0], relativeMarginsXXSec2[2], relativeMarginsYXSec2[1], relativeMarginsYXSec2[2]);
        else DrawGammaPadSettings( padRatioGenerators[i], relativeMarginsXXSec2[0], relativeMarginsXXSec2[2], relativeMarginsYXSec2[1], relativeMarginsYXSec2[1]);
        padRatioGenerators[i]->Draw();

        if (padRatioGenerators[i]->XtoPixel(padRatioGenerators[i]->GetX2()) < padRatioGenerators[i]->YtoPixel(padRatioGenerators[i]->GetY1())){
            textSizeLablesGenerators[i]            = (Double_t)textSizeLabelsPixel/padRatioGenerators[i]->XtoPixel(padRatioGenerators[i]->GetX2()) ;
            textsizeFacXSecGenerators[i]           = (Double_t)1./padRatioGenerators[i]->XtoPixel(padRatioGenerators[i]->GetX2()) ;
        } else {
            textSizeLablesGenerators[i]            = (Double_t)textSizeLabelsPixel/padRatioGenerators[i]->YtoPixel(padRatioGenerators[i]->GetY1());
            textsizeFacXSecGenerators[i]           = (Double_t)1./padRatioGenerators[i]->YtoPixel(padRatioGenerators[i]->GetY1());
        }
    }

    Double_t   yPosLabel[5]         = {0.84, 0.85, 0.85, 0.85, 0.89};

    TLatex *labelEnergyRatioPaper[5]   = {NULL, NULL, NULL, NULL, NULL};
    TLatex *labelDetSysRatioPaperPi0= new TLatex(0.19,yPosLabel[0],"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysRatioPaperPi0, textSizeLablesGenerators[0],4, 1, 42, kTRUE, 11);
    TLatex *labelDetSysRatioPaperEta= new TLatex(0.19,yPosLabel[0],"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysRatioPaperEta, textSizeLablesGenerators[0],4, 1, 42, kTRUE, 11);

    TLegend* legendTheoRatio = GetAndSetLegend2(0.18, yPosLabel[0]-0.03-(1*textSizeLablesGenerators[0]), 0.84, yPosLabel[0]-0.03, textSizeLabelsPixel, 3, "", 43, 0.22);
    legendTheoRatio->AddEntry(histoDPMJetPi0[0], "DPMJet", "l");
    legendTheoRatio->AddEntry(histoEPOSLHCPi0[0], "EPOS-LHC", "l");
    legendTheoRatio->AddEntry(histoHIJINGPi0[0], "HIJING", "l");



    for (Int_t cent = 0; cent < 5; cent++ ){
        padRatioGenerators[cent]->cd();
        padRatioGenerators[cent]->SetLogx(1);

        TH2F* ratio2DDataAndTheo               = new TH2F(Form ("ratio2DDataAndTheo%i",cent),Form ("ratio2DDataAndTheo%i",cent),1000,minPtPi0Plotting, maxPtPi0Plotting,1000,0.3,2.2);
        SetStyleHistoTH2ForGraphs(ratio2DDataAndTheo, "#it{p}_{T} (GeV/#it{c})","#frac{Data, Generator}{fit}", 0.85*textSizeLablesGenerators[cent], textSizeLablesGenerators[cent],
                                  0.85*textSizeLablesGenerators[cent],textSizeLablesGenerators[cent], 1,0.257/(textsizeFacXSecGenerators[cent]*marginXSec), 510, 508);
        ratio2DDataAndTheo->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DDataAndTheo->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DDataAndTheo->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DDataAndTheo->GetXaxis()->SetLabelFont(42);
        ratio2DDataAndTheo->GetYaxis()->SetLabelFont(42);
        ratio2DDataAndTheo->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DDataAndTheo->GetXaxis()->SetTickLength(0.07);
        ratio2DDataAndTheo->DrawCopy();

        labelEnergyRatioPaper[cent]     = new TLatex(0.935, yPosLabel[cent],Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
        SetStyleTLatex( labelEnergyRatioPaper[cent], textSizeLablesGenerators[cent],4, 1, 42, kTRUE, 31);
            boxErrorNormRatioPi0->Draw();
            if (cent == 0) legendTheoRatio->Draw();

            if (histoRatioPi0EPOSLHCToFit[cent]){
                SetStyleHisto(histoRatioPi0EPOSLHCToFit[cent], widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );
                histoRatioPi0EPOSLHCToFit[cent]->Draw("same,hist,l");
            }
            if (histoRatioPi0HIJINGToFit[cent]){
                SetStyleHisto(histoRatioPi0HIJINGToFit[cent], widthCommonFit*1.5, styleLineHIJING, colorHIJING );
                histoRatioPi0HIJINGToFit[cent]->Draw("same,hist,l");
            }
            if (histoRatioPi0DPMJetToFit[cent]){
                SetStyleHisto(histoRatioPi0DPMJetToFit[cent], widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );
                histoRatioPi0DPMJetToFit[cent]->Draw("same,hist,l");
            }
            if (graphRatioPi0CombCombFitStatWOXErr[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
                graphRatioPi0CombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
                DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
                graphRatioPi0CombCombFitSys[cent]->SetLineWidth(0);
                graphRatioPi0CombCombFitSys[cent]->Draw("2,same");
                graphRatioPi0CombCombFitStatWOXErr[cent]->Draw("p,same");
            }
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting,1., 1.,0.1,kGray+1);
            labelEnergyRatioPaper[cent]->Draw();
            if (cent == 0) labelDetSysRatioPaperPi0->Draw();
            ratio2DDataAndTheo->DrawCopy("axis,same");
        delete ratio2DDataAndTheo;
    }
    canvasRatioGeneratorsPaper->Print(Form("%s/Pi0_RatioToGenerator_Paper.%s",outputDir.Data(),suffix.Data()));

    canvasRatioGeneratorsPaper->cd();
    for (Int_t cent = 0; cent < 5; cent++ ){
        padRatioGenerators[cent]->cd();
        padRatioGenerators[cent]->SetLogx(1);

        TH2F* ratio2DDataAndTheo               = new TH2F(Form ("ratio2DDataAndTheo%i",cent),Form ("ratio2DDataAndTheo%i",cent),1000,minPtEtaPlotting, maxPtEtaPlotting,1000,0.3,2.2);
        SetStyleHistoTH2ForGraphs(ratio2DDataAndTheo, "#it{p}_{T} (GeV/#it{c})","#frac{Data, Generator}{fit}", 0.85*textSizeLablesGenerators[cent], textSizeLablesGenerators[cent],
                                  0.85*textSizeLablesGenerators[cent],textSizeLablesGenerators[cent], 1,0.257/(textsizeFacXSecGenerators[cent]*marginXSec), 510, 508);
        ratio2DDataAndTheo->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DDataAndTheo->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DDataAndTheo->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DDataAndTheo->GetXaxis()->SetLabelFont(42);
        ratio2DDataAndTheo->GetYaxis()->SetLabelFont(42);
        ratio2DDataAndTheo->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DDataAndTheo->GetXaxis()->SetTickLength(0.07);
        ratio2DDataAndTheo->DrawCopy();

        boxErrorNormRatioEta->Draw();
        if (cent == 0) legendTheoRatio->Draw();

        if (histoRatioEtaEPOSLHCToFit[cent]){
            SetStyleHisto(histoRatioEtaEPOSLHCToFit[cent], widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );
            histoRatioEtaEPOSLHCToFit[cent]->Draw("same,hist,l");
        }
        if (histoRatioEtaHIJINGToFit[cent]){
            SetStyleHisto(histoRatioEtaHIJINGToFit[cent], widthCommonFit*1.5, styleLineHIJING, colorHIJING );
            histoRatioEtaHIJINGToFit[cent]->Draw("same,hist,l");
        }
        if (histoRatioEtaDPMJetToFit[cent]){
            SetStyleHisto(histoRatioEtaDPMJetToFit[cent], widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );
            histoRatioEtaDPMJetToFit[cent]->Draw("same,hist,l");
        }
        if (graphRatioEtaCombCombFitStatWOXErr[cent]){
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatWOXErr[cent], markerStyleCentMC[cent], markerSizeCentMC[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
            graphRatioEtaCombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys[cent], markerStyleCentMC[cent], markerSizeCentMC[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioEtaCombCombFitSys[cent]->SetLineWidth(0);
            graphRatioEtaCombCombFitSys[cent]->Draw("2,same");
            graphRatioEtaCombCombFitStatWOXErr[cent]->Draw("p,same");
        }
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting,1., 1.,0.1,kGray+1);
        if (labelEnergyRatioPaper[cent]) labelEnergyRatioPaper[cent]->Draw();
        if (cent == 0) labelDetSysRatioPaperEta->Draw();
        ratio2DDataAndTheo->DrawCopy("axis,same");
        delete ratio2DDataAndTheo;
    }
    canvasRatioGeneratorsPaper->Print(Form("%s/Eta_RatioToGenerator_Paper.%s",outputDir.Data(),suffix.Data()));

    // ***************************************************************************************************************
    // ******************************** fitting eta/pi0 **************************************************************
    // ***************************************************************************************************************

    TF1* etaToPi0ConstData[17];
    TF1* etaToPi0ConstDataStat[17];
    TF1* etaToPi0ConstMC[17];
    TF1* etaToPi0ConstMC2[17];
    TF1* etaToPi0ConstMC3[17];
    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCentComb) continue;
        etaToPi0ConstData[cent]         = new TF1(Form("etaToPi0ConstData_%d",cent),"[0]",4,xPtLimitsEtaToPi0[cent][maxNBinsEtaToPi0[cent]]);
        etaToPi0ConstDataStat[cent]     = new TF1(Form("etaToPi0ConstDataStat_%d",cent),"[0]",4,xPtLimitsEtaToPi0[cent][maxNBinsEtaToPi0[cent]]);
        etaToPi0ConstMC[cent]           = new TF1(Form("etaToPi0ConstMCStat_%d",cent),"[0]",4,xPtLimitsEtaToPi0[cent][maxNBinsEtaToPi0[cent]]);
        etaToPi0ConstMC2[cent]          = new TF1(Form("etaToPi0ConstMC2_%d",cent),"[0]",4,xPtLimitsEtaToPi0[cent][maxNBinsEtaToPi0[cent]]);
        etaToPi0ConstMC3[cent]          = new TF1(Form("etaToPi0ConstMC3_%d",cent),"[0]",4,xPtLimitsEtaToPi0[cent][maxNBinsEtaToPi0[cent]]);

        if (graphCombEtaToPi0StatWOXErr[cent])  graphCombEtaToPi0StatWOXErr[cent]->Fit(etaToPi0ConstDataStat[cent],"QRME0","",4,xPtLimitsEtaToPi0[cent][maxNBinsEtaToPi0[cent]]);
        if (graphCombEtaToPi0Tot[cent])  graphCombEtaToPi0Tot[cent]->Fit(etaToPi0ConstData[cent],"QRME0","",4,xPtLimitsEtaToPi0[cent][maxNBinsEtaToPi0[cent]]);
        if (histoDPMJetEtaToPi0[cent])  histoDPMJetEtaToPi0[cent]->Fit(etaToPi0ConstMC[cent],"QRME0","",4,xPtLimitsEtaToPi0[cent][maxNBinsEtaToPi0[cent]]);
        if (histoHIJINGEtaToPi0[cent])  histoHIJINGEtaToPi0[cent]->Fit(etaToPi0ConstMC2[cent],"QRME0","",4,xPtLimitsEtaToPi0[cent][maxNBinsEtaToPi0[cent]]);
        if (histoEPOSLHCEtaToPi0[cent])  histoEPOSLHCEtaToPi0[cent]->Fit(etaToPi0ConstMC3[cent],"QRME0","",4,xPtLimitsEtaToPi0[cent][maxNBinsEtaToPi0[cent]]);

        cout << "***********************************************************************************************************" << endl;
        cout << "***********************************************************************************************************" << endl;
        cout << "***********************************************************************************************************" << endl;
        cout << centArray[cent].Data() << "\t" << addCentString[cent].Data() << "\t" << runArray[cent].Data() << endl;
        cout << "high pt eta/pi0 - data, stat: " << etaToPi0ConstDataStat[cent]->GetParameter(0) << "+-"<< etaToPi0ConstDataStat[cent]->GetParError(0) << endl;
        cout << "high pt eta/pi0 - data, tot: " << etaToPi0ConstData[cent]->GetParameter(0) << "+-"<< etaToPi0ConstData[cent]->GetParError(0) << endl;
        if(histoDPMJetEtaToPi0[cent]) cout << "high pt eta/pi0 - DPMJet: " << etaToPi0ConstMC[cent]->GetParameter(0) << "+-"<< etaToPi0ConstMC[cent]->GetParError(0) << endl;
        if(histoHIJINGEtaToPi0[cent]) cout << "high pt eta/pi0 - HIJING: " << etaToPi0ConstMC2[cent]->GetParameter(0) << "+-"<< etaToPi0ConstMC2[cent]->GetParError(0) << endl;
        if(histoEPOSLHCEtaToPi0[cent]) cout << "high pt eta/pi0 - EPOS-LHC: " << etaToPi0ConstMC3[cent]->GetParameter(0) << "+-"<< etaToPi0ConstMC3[cent]->GetParError(0) << endl;
        cout << "***********************************************************************************************************" << endl;
        cout << "***********************************************************************************************************" << endl;
        cout << "***********************************************************************************************************" << endl;

        fileFitsOutput << "***********************************************************************************************************" << endl;
        fileFitsOutput << "***********************************************************************************************************" << endl;
        fileFitsOutput << "***********************************************************************************************************" << endl;
        fileFitsOutput << centArray[cent].Data() << "\t" << addCentString[cent].Data() << "\t" << runArray[cent].Data() << endl;
        fileFitsOutput << "high pt eta/pi0 - data, stat: " << etaToPi0ConstDataStat[cent]->GetParameter(0) << "+-"<< etaToPi0ConstDataStat[cent]->GetParError(0) << endl;
        fileFitsOutput << "high pt eta/pi0 - data, tot: " << etaToPi0ConstData[cent]->GetParameter(0) << "+-"<< etaToPi0ConstData[cent]->GetParError(0) << endl;
        if(histoDPMJetEtaToPi0[cent]) fileFitsOutput << "high pt eta/pi0 - DPMJet: " << etaToPi0ConstMC[cent]->GetParameter(0) << "+-"<< etaToPi0ConstMC[cent]->GetParError(0) << endl;
        if(histoHIJINGEtaToPi0[cent]) fileFitsOutput << "high pt eta/pi0 - HIJING: " << etaToPi0ConstMC2[cent]->GetParameter(0) << "+-"<< etaToPi0ConstMC2[cent]->GetParError(0) << endl;
        if(histoEPOSLHCEtaToPi0[cent]) fileFitsOutput << "high pt eta/pi0 - EPOS-LHC: " << etaToPi0ConstMC3[cent]->GetParameter(0) << "+-"<< etaToPi0ConstMC3[cent]->GetParError(0) << endl;
        fileFitsOutput << "***********************************************************************************************************" << endl;
        fileFitsOutput << "***********************************************************************************************************" << endl;
        fileFitsOutput << "***********************************************************************************************************" << endl;

    }
    TF1* etaToPi0OlderData              = new TF1("etaToPi0ConstOlderData","[0]",4,20);
    TF1* etaToPi0OlderDataStat          = new TF1("etaToPi0ConstOlderDataStat","[0]",4,20);
    graphCombEtaToPi0OlderMBTot->Fit(etaToPi0OlderData, "QRME0","",4,20);
    graphCombEtaToPi0OlderMBStat->Fit(etaToPi0OlderDataStat, "QRME0","",4,20);

    // ***************************************************************************************************************
    // ******************************* Plotting eta/pi0 ratio for single measurements ********************************
    // ***************************************************************************************************************
    textSizeLabelsPixel                 = 54;
    TCanvas* canvasEtatoPi0combo       = new TCanvas("canvasEtatoPi0combo","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasEtatoPi0combo, 0.1, 0.03, 0.01, 0.125);
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

    TH2F * histo2DEtatoPi0combo     = new TH2F("histo2DEtatoPi0combo","histo2DEtatoPi0combo",1000,minPtEtaToPi0Plotting, maxPtEtaToPi0Plotting,1000,0.,1.05    );
    SetStyleHistoTH2ForGraphs(histo2DEtatoPi0combo, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.85*textsizeLabelsEtaToPi0, textsizeLabelsEtaToPi0,
                              0.85*textsizeLabelsEtaToPi0,1.1*textsizeLabelsEtaToPi0, 0.9, 0.65, 510, 510);
    histo2DEtatoPi0combo->GetXaxis()->SetMoreLogLabels();
    histo2DEtatoPi0combo->GetXaxis()->SetLabelOffset(-0.01);
    histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.0,1.01);
    histo2DEtatoPi0combo->Draw("copy");
    TH2F * histo2DParticleRatios     = new TH2F("histo2DParticleRatios","histo2DParticleRatios",1000,minPtEtaToPi0Plotting, maxPtEtaToPi0Plotting,1000,0.,1.05    );
    SetStyleHistoTH2ForGraphs(histo2DParticleRatios, "#it{p}_{T} (GeV/#it{c})","particle ratio", 0.85*textsizeLabelsEtaToPi0, textsizeLabelsEtaToPi0,
                              0.85*textsizeLabelsEtaToPi0,1.1*textsizeLabelsEtaToPi0, 0.9, 0.65, 510, 510);
    histo2DParticleRatios->GetXaxis()->SetMoreLogLabels();
    histo2DParticleRatios->GetXaxis()->SetLabelOffset(-0.01);
    histo2DParticleRatios->GetYaxis()->SetRangeUser(0.0,1.01);
    histo2DParticleRatios->Draw("copy");

    TLatex *labelEnergyEtaToPi0 = new TLatex(0.13, 0.92,collisionSystempPb.Data());
    SetStyleTLatex( labelEnergyEtaToPi0, textsizeLabelsEtaToPi0,4, 1, 42, kTRUE, 11);
    labelEnergyEtaToPi0->Draw();
    TLatex *labelALICEEtaToPi0 = new TLatex(0.13, 0.92-(1*textsizeLabelsEtaToPi0),textALICE.Data());
    SetStyleTLatex( labelALICEEtaToPi0, textsizeLabelsEtaToPi0,4, 1, 42, kTRUE, 11);
    labelALICEEtaToPi0->Draw();

    Double_t totErrRelEtaToPi0MB          = etaToPi0ConstData[0]->GetParError(0)/etaToPi0ConstData[0]->GetParameter(0);
    Double_t totErrEtaToPi0MB             = totErrRelEtaToPi0MB*etaToPi0ConstDataStat[0]->GetParameter(0);
    Double_t statErrRelEtaToPi0MB         = etaToPi0ConstDataStat[0]->GetParError(0)/etaToPi0ConstDataStat[0]->GetParameter(0);
    Double_t statErrEtaToPi0MB            = etaToPi0ConstDataStat[0]->GetParError(0);
    Double_t sysErrRelEtaToPi0MB          = TMath::Sqrt(pow(totErrRelEtaToPi0MB,2)-pow(statErrRelEtaToPi0MB,2));
    Double_t sysErrEtaToPi0MB             = sysErrRelEtaToPi0MB*etaToPi0ConstDataStat[0]->GetParameter(0);

    TBox* boxHighPtFitMB                  = CreateBoxConv(kGray, 4,etaToPi0ConstDataStat[0]->GetParameter(0)-totErrEtaToPi0MB, maxPtEtaToPi0Plotting, etaToPi0ConstDataStat[0]->GetParameter(0)+totErrEtaToPi0MB);
    TGraphErrors* dummyHighPtMB           = new TGraphErrors(1);
    dummyHighPtMB->SetPoint(0,1,0.0251);
    DrawGammaSetMarkerTGraphErr(dummyHighPtMB, 0, 0, kGray, kGray, widthLinesBoxes, kTRUE, kGray);
    dummyHighPtMB->SetLineStyle(7);
    dummyHighPtMB->SetLineColor(kGray+2);
    dummyHighPtMB->SetLineWidth(2);

    Double_t totErrRelEtaToPi0OldMB          = etaToPi0OlderData->GetParError(0)/etaToPi0OlderData->GetParameter(0);
    Double_t totErrEtaToPi0OldMB             = totErrRelEtaToPi0OldMB*etaToPi0OlderDataStat->GetParameter(0);
    Double_t statErrRelEtaToPi0OldMB         = etaToPi0OlderDataStat->GetParError(0)/etaToPi0OlderDataStat->GetParameter(0);
    Double_t statErrEtaToPi0OldMB            = etaToPi0OlderDataStat->GetParError(0);
    Double_t sysErrRelEtaToPi0OldMB          = TMath::Sqrt(pow(totErrRelEtaToPi0OldMB,2)-pow(statErrRelEtaToPi0OldMB,2));
    Double_t sysErrEtaToPi0OldMB             = sysErrRelEtaToPi0OldMB*etaToPi0OlderDataStat->GetParameter(0);

    TBox* boxHighPtFitOldMB                  = CreateBoxConv(kGray, 4,etaToPi0OlderDataStat->GetParameter(0)-totErrEtaToPi0OldMB, maxPtEtaToPi0Plotting, etaToPi0OlderDataStat->GetParameter(0)+totErrEtaToPi0OldMB);
    TGraphErrors* dummyHighPtOldMB           = new TGraphErrors(1);
    dummyHighPtOldMB->SetPoint(0,1,0.0251);
    DrawGammaSetMarkerTGraphErr(dummyHighPtOldMB, 0, 0, kGray, kGray, widthLinesBoxes, kTRUE, kGray);
    dummyHighPtOldMB->SetLineStyle(7);
    dummyHighPtOldMB->SetLineColor(kGray+2);
    dummyHighPtOldMB->SetLineWidth(2);


    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCentComb[cent] ) continue;
        histo2DEtatoPi0combo->Draw("copy");
        for (Int_t meth = 10; meth > -1 ; meth--){
            if (!fileNamesMethod[meth].CompareTo("")) continue;
            if (graphEtaToPi0Sys[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Sys[cent][meth], markerStyleDet[meth] ,markerSizeDet[meth]*0.5, colorDet[meth], colorDet[meth], widthLinesBoxes, kTRUE);
                graphEtaToPi0Sys[cent][meth]->Draw("E2same");
            }
            if (graphEtaToPi0Stat[cent][meth]){
                ProduceGraphAsymmWithoutXErrors(graphEtaToPi0Stat[cent][meth]);
                DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Stat[cent][meth], markerStyleDet[meth] ,markerSizeDet[meth]*0.5, colorDet[meth], colorDet[meth]);
                graphEtaToPi0Stat[cent][meth]->Draw("p,same,z");
            }
        }
        if (graphEtaToPi0Stat[cent][4]){
            graphEtaToPi0Stat[cent][4]->Draw("p,same,z");
        }
        TLegend* legendEtaToPi0 = GetAndSetLegend2(0.13, 0.92-(1*textsizeLabelsEtaToPi0)-0.01, 0.6, 0.92-(4*textsizeLabelsEtaToPi0)-0.01, textSizeLabelsPixel,2, "", 43, 0.25);
        if (graphEtaToPi0Sys[cent][0]) legendEtaToPi0->AddEntry(graphEtaToPi0Sys[cent][0],nameMeasGlobalLabel[0],"pf");
        if (graphEtaToPi0Sys[cent][4]) legendEtaToPi0->AddEntry(graphEtaToPi0Sys[cent][4],nameMeasGlobalLabel[4],"pf");
        if (graphEtaToPi0Sys[cent][1]) legendEtaToPi0->AddEntry(graphEtaToPi0Sys[cent][1],nameMeasGlobalLabel[1],"pf");
        if (graphEtaToPi0Sys[cent][3]) legendEtaToPi0->AddEntry(graphEtaToPi0Sys[cent][3],nameMeasGlobalLabel[3],"pf");
        if (graphEtaToPi0Sys[cent][2]) legendEtaToPi0->AddEntry(graphEtaToPi0Sys[cent][2],nameMeasGlobalLabel[2],"pf");
        legendEtaToPi0->Draw();

        labelEnergyEtaToPi0->SetText(0.13, 0.92,Form("%s %s", centArray2[cent].Data(), collisionSystempPb.Data()));
        labelEnergyEtaToPi0->Draw();
        labelALICEEtaToPi0->Draw();

        histo2DEtatoPi0combo->Draw("same,axis");

        canvasEtatoPi0combo->Update();
        canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_differentSystems_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));


        histo2DEtatoPi0combo->Draw("copy");
        for (Int_t meth = 10; meth > -1 ; meth--){
            if (!fileNamesMethod[meth].CompareTo("")) continue;
            if (graphEtaToPi0Sys[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Sys[cent][meth], markerStyleDet[meth] ,markerSizeDet[meth]*0.5, colorDet[meth], colorDet[meth], widthLinesBoxes, kTRUE);
                graphEtaToPi0Sys[cent][meth]->Draw("E2same");
            }
            if (graphEtaToPi0Stat[cent][meth]){
                ProduceGraphAsymmWithoutXErrors(graphEtaToPi0Stat[cent][meth]);
                DrawGammaSetMarkerTGraphAsym(graphEtaToPi0Stat[cent][meth], markerStyleDet[meth] ,markerSizeDet[meth]*0.5, colorDet[meth], colorDet[meth]);
                graphEtaToPi0Stat[cent][meth]->Draw("p,same,z");
            }
        }

        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
        graphCombEtaToPi0Sys[cent]->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0StatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent]);
        graphCombEtaToPi0StatWOXErr[cent]->Draw("p,same,z");
        if (graphCombEtaToPi0Sys[cent]) legendEtaToPi0->AddEntry(graphCombEtaToPi0Sys[cent],"comb","pf");

        legendEtaToPi0->Draw();
        labelEnergyEtaToPi0->Draw();
        labelALICEEtaToPi0->Draw();

        canvasEtatoPi0combo->Update();
        canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_differentSystemsWComb_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));

        histo2DEtatoPi0combo->Draw("copy");

            // eta/pi0 mt-scaled
            TH1F *eta2pi0MtScaledTCM = new TH1F("eta2pi0MtScaledTCM","#eta/#pi^{0} from m_{T} scaling",5000,0.4,maxPtEtaToPi0Plotting);
            eta2pi0MtScaledTCM->SetLineColor(kBlue+2);
            eta2pi0MtScaledTCM->SetLineWidth(2.);

            Double_t eta2Pi0Const   = etaToPi0ConstDataStat[cent]->GetParameter(0);
            Double_t mPi0           = 0.134977;
            Double_t mEta           = 0.547853;
            for (Int_t i=1; i<=eta2pi0MtScaledTCM->GetNbinsX(); i++) {
                Double_t ptPi0          = eta2pi0MtScaledTCM->GetBinCenter(i);
                if (ptPi0 < 0.3) continue;
                Double_t mtEta          = TMath::Sqrt(mEta*mEta + ptPi0*ptPi0);
                Double_t ptEta          = TMath::Sqrt(mtEta*mtEta - mPi0*mPi0);
                Double_t Reta2pi0TCM    = fitTCMInvYieldPi0[cent]->Eval(ptEta) / fitTCMInvYieldPi0[cent]->Eval(ptPi0) * eta2Pi0Const;
                eta2pi0MtScaledTCM->SetBinContent(i,Reta2pi0TCM);
            }


            Double_t totErrRelEtaToPi0          = etaToPi0ConstData[cent]->GetParError(0)/etaToPi0ConstData[cent]->GetParameter(0);
            Double_t totErrEtaToPi0             = totErrRelEtaToPi0*etaToPi0ConstDataStat[cent]->GetParameter(0);
            Double_t statErrRelEtaToPi0         = etaToPi0ConstDataStat[cent]->GetParError(0)/etaToPi0ConstDataStat[cent]->GetParameter(0);
            Double_t statErrEtaToPi0            = etaToPi0ConstDataStat[cent]->GetParError(0);
            Double_t sysErrRelEtaToPi0          = TMath::Sqrt(pow(totErrRelEtaToPi0,2)-pow(statErrRelEtaToPi0,2));
            Double_t sysErrEtaToPi0             = sysErrRelEtaToPi0*etaToPi0ConstDataStat[cent]->GetParameter(0);

            TBox* boxHighPtFit                  = CreateBoxConv(kGray, 4,etaToPi0ConstDataStat[cent]->GetParameter(0)-totErrEtaToPi0, maxPtEtaToPi0Plotting, etaToPi0ConstDataStat[cent]->GetParameter(0)+totErrEtaToPi0);
            boxHighPtFit->Draw();
            TGraphErrors* dummyHighPt           = new TGraphErrors(1);
            dummyHighPt->SetPoint(0,1,0.0251);
            DrawGammaSetMarkerTGraphErr(dummyHighPt, 0, 0, kGray, kGray, widthLinesBoxes, kTRUE, kGray);
            dummyHighPt->SetLineStyle(7);
            dummyHighPt->SetLineColor(kGray+2);
            dummyHighPt->SetLineWidth(2);
//             DrawGammaSetMarkerTGraphErr(graphMcGillEtaToPi0, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
//             graphMcGillEtaToPi0->Draw("3,same");

            Int_t nTheory = 0;
            if (histoDPMJetEtaToPi0[cent]){
                SetStyleHisto(histoDPMJetEtaToPi0[cent], widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );
                histoDPMJetEtaToPi0[cent]->GetXaxis()->SetRangeUser(minPtPi0[cent],xPtLimitsPi0[cent][maxNBinsPi0[cent]]);
                histoDPMJetEtaToPi0[cent]->Draw("same,hist,l");
                nTheory++;
            }
            if (histoHIJINGEtaToPi0[cent]){
                SetStyleHisto(histoHIJINGEtaToPi0[cent], widthCommonFit*1.5, styleLineHIJING, colorHIJING );
                histoHIJINGEtaToPi0[cent]->GetXaxis()->SetRangeUser(minPtPi0[cent],xPtLimitsPi0[cent][maxNBinsPi0[cent]]);
                histoHIJINGEtaToPi0[cent]->Draw("same,hist,l");
                nTheory++;
            }
            if (histoEPOSLHCEtaToPi0[cent]){
                SetStyleHisto(histoEPOSLHCEtaToPi0[cent], widthCommonFit*1.5,  styleLineEPOS3, colorEPOS3 );
                histoEPOSLHCEtaToPi0[cent]->GetXaxis()->SetRangeUser(minPtPi0[cent],xPtLimitsPi0[cent][maxNBinsPi0[cent]]);
                histoEPOSLHCEtaToPi0[cent]->Draw("same,hist,l");
                nTheory++;
            }

            SetStyleHisto(eta2pi0MtScaledTCM, widthCommonFit*1.5, 1, kRed+2);
            eta2pi0MtScaledTCM->Draw("same,hist,l");

            DrawGammaSetMarkerTF1( etaToPi0ConstDataStat[cent], 7, 2, kGray+2);
            etaToPi0ConstDataStat[cent]->SetRange(4,maxPtEtaToPi0Plotting);
            etaToPi0ConstDataStat[cent]->Draw("same");


            graphCombEtaToPi0Sys[cent]->Draw("E2same");
            graphCombEtaToPi0StatWOXErr[cent]->Draw("p,same,z");

            // plotting labels
            labelEnergyEtaToPi0->Draw();
            labelALICEEtaToPi0->Draw();

            TLegend* legendEtaToPi0TheoryHighpt = GetAndSetLegend2(0.13, 0.905-(4*textsizeLabelsEtaToPi0*1.05), 0.53, 0.905-(1*textsizeLabelsEtaToPi0), textSizeLabelsPixel*0.85, 1, "", 43, 0.14);
            legendEtaToPi0TheoryHighpt->AddEntry(dummyHighPt, "high #it{p}_{T} average","fl");
            legendEtaToPi0TheoryHighpt->AddEntry((TObject*)0, Form ("%2.3f #pm %1.3f^{stat} #pm %1.3f^{sys}", etaToPi0ConstDataStat[cent]->GetParameter(0), statErrEtaToPi0, sysErrRelEtaToPi0),"");
            legendEtaToPi0TheoryHighpt->AddEntry(eta2pi0MtScaledTCM, "m_{T}-scaled #eta/#pi^{0}","l");
            legendEtaToPi0TheoryHighpt->Draw();


            TLegend* legendEtaToPi0Theory = GetAndSetLegend2(0.68, 0.15+(textsizeLabelsEtaToPi0*nTheory*0.9), 0.88, 0.15, textSizeLabelsPixel*0.85, 1, "", 43, 0.26);
//             legendEtaToPi0Theory->AddEntry((TObject*)0,"","");
            if (histoDPMJetEtaToPi0[cent])legendEtaToPi0Theory->AddEntry(histoDPMJetEtaToPi0[cent],"DPMJet","l");
            if (histoEPOSLHCEtaToPi0[cent])legendEtaToPi0Theory->AddEntry(histoEPOSLHCEtaToPi0[cent],"EPOS-LHC","l");
            if (histoHIJINGEtaToPi0[cent])legendEtaToPi0Theory->AddEntry(histoHIJINGEtaToPi0[cent],"HIJING","l");
//             legendEtaToPi0Theory->AddEntry(graphMcGillEtaToPi0,"iEBE-VISHNU","f");
    //         legendEtaToPi0Theory->AddEntry(graphMcGillEtaToPi0,"Shen #it{et al.}","f");
            legendEtaToPi0Theory->Draw();

        histo2DEtatoPi0combo->Draw("axis,same");

        canvasEtatoPi0combo->Update();
        canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Theory_Paper_%s%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(),suffix.Data()));

        if (graphChKToPiTot[cent]){
            histo2DParticleRatios->Draw("copy");

                graphCombEtaToPi0Sys[cent]->Draw("E2same");
                graphCombEtaToPi0StatWOXErr[cent]->Draw("p,same,z");

//                 DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
//                 graphCombEtaToPi0Sys[cent]->Draw("E2same");
                DrawGammaSetMarkerTGraphAsym(graphChKToPiTot[cent], 25, markerSizeCent[cent]*0.5, kGray+1, kGray+1);
                graphChKToPiTot[cent]->Draw("p,same,z");

                // plotting labels
                labelEnergyEtaToPi0->Draw();
                labelALICEEtaToPi0->Draw();

                TLegend* legendParticleRatios = GetAndSetLegend2(0.13, 0.905-(3*textsizeLabelsEtaToPi0*1.05), 0.53, 0.905-(1*textsizeLabelsEtaToPi0), textSizeLabelsPixel*0.85, 1, "", 43, 0.14);
                legendParticleRatios->AddEntry(graphCombEtaToPi0Sys[cent], "#eta/#pi^{0}","fp");
                legendParticleRatios->AddEntry(graphChKToPiTot[cent], "K^{#pm}/#pi^{#pm}","ep");
                legendParticleRatios->Draw();

            canvasEtatoPi0combo->Update();
            canvasEtatoPi0combo->SaveAs(Form("%s/ParticleRatios_%s%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(),suffix.Data()));
            delete legendParticleRatios;

        }

        delete legendEtaToPi0TheoryHighpt;
        delete dummyHighPt;
        delete boxHighPtFit;
        delete eta2pi0MtScaledTCM;
        delete legendEtaToPi0;
        delete legendEtaToPi0Theory;
    }


    for (Int_t est = 0; est < 3; est++){

        Int_t nCurrEst      = 0;
        Bool_t haveCharged  = kFALSE;
        for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
            if ( (centArray2[cent].Contains(nameCentEst[est].Data()) || centArray2[cent].CompareTo("0-100%") == 0) && cent != 5  && enableCentComb[cent]){
                nCurrEst++;
                if ( graphChKToPiTot[cent] )
                    haveCharged = kTRUE;
            }
        }

        histo2DEtatoPi0combo->Draw("copy");
        TLegend* legendEtaToPi0Comb = GetAndSetLegend2(0.13, 0.91-(1*textsizeLabelsEtaToPi0*1.)-0.01, 0.6, 0.91-(nCurrEst*textsizeLabelsEtaToPi0*1.01)-0.01, textSizeLabelsPixel,1, "", 43, 0.125);

            for (Int_t cent = 0; cent < maxCentRun1; cent++){
                if (! (centArray2[cent].Contains(nameCentEst[est].Data()) || centArray2[cent].CompareTo("0-100%") == 0)) continue;
                if (!enableCentComb[cent] || cent == 1) continue;
                DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys[cent], markerStyleCent[cent], markerSizeCent[cent]*0.55, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphCombEtaToPi0Sys[cent]->Draw("E2same");
                DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0StatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.55, colorCent[cent], colorCent[cent]);
                graphCombEtaToPi0StatWOXErr[cent]->Draw("p,same,z");
                legendEtaToPi0Comb->AddEntry(graphCombEtaToPi0Sys[cent],centArray[cent].Data(),"pf");
            }

            if (enableCentComb[0] ){
                graphCombEtaToPi0Sys[0]->Draw("E2same");
                graphCombEtaToPi0StatWOXErr[0]->Draw("p,same,z");
            }

            // plotting labels
            labelEnergyEtaToPi0->SetText(0.13, 0.92,Form("%s %s", nameCentEst[est].Data(), collisionSystempPb.Data()));
            labelEnergyEtaToPi0->Draw();
            labelALICEEtaToPi0->Draw();
            legendEtaToPi0Comb->Draw();
            histo2DEtatoPi0combo->Draw("axis,same");

        canvasEtatoPi0combo->Update();
        canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Paper_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));


        histo2DEtatoPi0combo->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0OlderMBSys, markerStyleCent[0], markerSizeCent[0]*0.55, colorCent[0], colorCent[0], widthLinesBoxes, kTRUE);
            graphCombEtaToPi0OlderMBSys->Draw("E2same");

            for (Int_t cent = 0; cent < maxCentRun1; cent++){
                if (! (centArray2[cent].Contains(nameCentEst[est].Data()) )) continue;
                if (!enableCentComb[cent] ) continue;
                DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys[cent], markerStyleCent[cent], markerSizeCent[cent]*0.55, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphCombEtaToPi0Sys[cent]->Draw("E2same");
                DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0StatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.55, colorCent[cent], colorCent[cent]);
                graphCombEtaToPi0StatWOXErr[cent]->Draw("p,same,z");
            }

            ProduceGraphAsymmWithoutXErrors(graphCombEtaToPi0OlderMBStat);
            DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0OlderMBStat, markerStyleCent[0], markerSizeCent[0]*0.55, colorCent[0], colorCent[0], widthLinesBoxes);
            graphCombEtaToPi0OlderMBStat->Draw("p,same,z");

            // plotting labels
            labelEnergyEtaToPi0->SetText(0.13, 0.92,Form("%s %s", nameCentEst[est].Data(), collisionSystempPb.Data()));
            labelEnergyEtaToPi0->Draw();
            labelALICEEtaToPi0->Draw();
            legendEtaToPi0Comb->Draw();
            histo2DEtatoPi0combo->Draw("axis,same");

        canvasEtatoPi0combo->Update();
        canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Paper_OldMB_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));

        histo2DEtatoPi0combo->Draw("copy");
        TLegend* legendEtaToPi0Comb3 = GetAndSetLegend2(0.45, 0.14, 0.9, 0.14+(2*textsizeLabelsEtaToPi0*1.), textSizeLabelsPixel,2, "", 43, 0.25);
        boxHighPtFitMB->Draw();


        for (Int_t cent = 0; cent < maxCentRun1; cent++){
            if (! (centArray2[cent].Contains(nameCentEst[est].Data()) )) continue;
            if (!enableCentComb[cent] ) continue;
            graphCombEtaToPi0Sys[cent]->Draw("E2same");
            graphCombEtaToPi0StatWOXErr[cent]->Draw("p,same,z");
            legendEtaToPi0Comb3->AddEntry(graphCombEtaToPi0Sys[cent],centArray[cent].Data(),"pf");
        }

        if (enableCentComb[0] ){
            graphCombEtaToPi0Sys[0]->Draw("E2same");
            graphCombEtaToPi0StatWOXErr[0]->Draw("p,same,z");
            DrawGammaSetMarkerTF1( etaToPi0ConstDataStat[0], 7, 2, kGray+2);
            etaToPi0ConstDataStat[0]->SetRange(4,maxPtEtaToPi0Plotting);
            etaToPi0ConstDataStat[0]->Draw("same");
            TLegend* legendEtaToPi0Highpt2 = GetAndSetLegend2(0.13, 0.905-(4*textsizeLabelsEtaToPi0*1.05), 0.53, 0.905-(1*textsizeLabelsEtaToPi0), textSizeLabelsPixel*0.85, 1, "", 43, 0.14);
            legendEtaToPi0Highpt2->AddEntry(graphCombEtaToPi0Sys[0], centArray[0].Data(),"pf");
            legendEtaToPi0Highpt2->AddEntry(dummyHighPtMB, Form ("high #it{p}_{T} average %s", centArray[0].Data() ),"fl");
            legendEtaToPi0Highpt2->AddEntry((TObject*)0, Form ("%2.3f #pm %1.3f^{stat} #pm %1.3f^{sys}", etaToPi0ConstDataStat[0]->GetParameter(0), statErrEtaToPi0MB, sysErrRelEtaToPi0MB),"");
            legendEtaToPi0Highpt2->Draw();
        }

        // plotting labels
        labelEnergyEtaToPi0->Draw();
        labelALICEEtaToPi0->Draw();
        legendEtaToPi0Comb3->Draw();


        histo2DEtatoPi0combo->Draw("axis,same");

        canvasEtatoPi0combo->Update();
        canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Paper_WithMBAverage_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));


        histo2DEtatoPi0combo->Draw("copy");
        boxHighPtFitOldMB->Draw();


        for (Int_t cent = 0; cent < maxCentRun1; cent++){
            if (! (centArray2[cent].Contains(nameCentEst[est].Data()) )) continue;
            if (!enableCentComb[cent] ) continue;
            graphCombEtaToPi0Sys[cent]->Draw("E2same");
            graphCombEtaToPi0StatWOXErr[cent]->Draw("p,same,z");
        }

        graphCombEtaToPi0OlderMBSys->Draw("E2same");
        graphCombEtaToPi0OlderMBStat->Draw("p,same,z");
        DrawGammaSetMarkerTF1( etaToPi0OlderDataStat, 7, 2, kGray+2);
        etaToPi0OlderDataStat->SetRange(4,maxPtEtaToPi0Plotting);
        etaToPi0OlderDataStat->Draw("same");
        TLegend* legendEtaToPi0Highpt3 = GetAndSetLegend2(0.13, 0.905-(4*textsizeLabelsEtaToPi0*1.05), 0.53, 0.905-(1*textsizeLabelsEtaToPi0), textSizeLabelsPixel*0.85, 1, "", 43, 0.14);
        legendEtaToPi0Highpt3->AddEntry(graphCombEtaToPi0OlderMBSys, centArray[0].Data(),"pf");
        legendEtaToPi0Highpt3->AddEntry(dummyHighPtOldMB, Form ("high #it{p}_{T} average %s", centArray[0].Data() ),"fl");
        legendEtaToPi0Highpt3->AddEntry((TObject*)0, Form ("%2.3f #pm %1.3f^{stat} #pm %1.3f^{sys}", etaToPi0OlderDataStat->GetParameter(0), statErrEtaToPi0OldMB, sysErrRelEtaToPi0OldMB),"");
        legendEtaToPi0Highpt3->Draw();

        // plotting labels
        labelEnergyEtaToPi0->Draw();
        labelALICEEtaToPi0->Draw();
        legendEtaToPi0Comb3->Draw();


        histo2DEtatoPi0combo->Draw("axis,same");

        canvasEtatoPi0combo->Update();
        canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Paper_WithMBAverageOldMB_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));


        if (haveCharged){

            histo2DParticleRatios->Draw("copy");

            TLegend* legendParticleRatios = GetAndSetLegend2(0.13, 0.905-(3*textsizeLabelsEtaToPi0*1.05), 0.63, 0.905-(1*textsizeLabelsEtaToPi0), textSizeLabelsPixel*0.85, 2, "", 43, 0.2);
            for (Int_t cent = 0; cent < maxCentRun1; cent++){
                if (! (centArray2[cent].Contains(nameCentEst[est].Data())))  continue;
                if (!(centArray[cent].CompareTo("0-20%") == 0 ||centArray[cent].CompareTo("60-100%") == 0 )) continue;
                if (!enableCentComb[cent] || cent == 1) continue;
                legendParticleRatios->AddEntry(graphCombEtaToPi0Sys[cent], "#eta/#pi^{0} "+centArray[cent],"fp");

                graphCombEtaToPi0Sys[cent]->Draw("E2same");
                graphCombEtaToPi0StatWOXErr[cent]->Draw("p,same,z");
            }
            for (Int_t cent = 0; cent < maxCentRun1; cent++){
                if (! (centArray2[cent].Contains(nameCentEst[est].Data()))) continue;
                if (! (centArray[cent].CompareTo("0-20%") == 0 ||centArray[cent].CompareTo("60-100%") == 0 )) continue;
                if (!enableCentComb[cent] || cent == 1) continue;
                DrawGammaSetMarkerTGraphAsym(graphChKToPiTot[cent], 25, markerSizeCent[cent]*0.5, colorCentBox[cent], colorCentBox[cent]);
                graphChKToPiTot[cent]->Draw("p,same,z");
                legendParticleRatios->AddEntry(graphChKToPiTot[cent], "K^{#pm}/#pi^{#pm} "+centArray[cent],"ep");

            //                 DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
            //                 graphCombEtaToPi0Sys[cent]->Draw("E2same");
            }
            legendParticleRatios->Draw();
            // plotting labels
            labelEnergyEtaToPi0->Draw();
            labelALICEEtaToPi0->Draw();


            canvasEtatoPi0combo->Update();
            canvasEtatoPi0combo->SaveAs(Form("%s/ParticleRatios_CentAndPeri_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));
            delete legendParticleRatios;
        }
    }


    //*****************************************************************************************************************
    // Plotting RpA s for different mesons
    //*****************************************************************************************************************
    textSizeLabelsPixel                         = 54;
    TCanvas* canvasRpPb = new TCanvas("canvasRpPb","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasRpPb,  0.085, 0.01, 0.015, 0.08);
    canvasRpPb->SetLogx();
    TH2F * histo2DRpPb  = new TH2F("histo2DRpPb","histo2DRpPb",1000,minPtPi0Plotting,maxPtPi0Plotting,1000,0.0,1.95);
    SetStyleHistoTH2ForGraphs(histo2DRpPb, "#it{p}_{T} (GeV/#it{c})","#it{Q}_{pA}", 0.035,0.04, 0.035,0.04, 0.8,1., 512, 505);
    histo2DRpPb->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRpPb->DrawCopy();
    TH2F * histo2DRpPbEta  = new TH2F("histo2DRpPbEta","histo2DRpPbEta",1000,minPtEtaPlotting,maxPtEtaPlotting,1000,0.0,1.75);
    SetStyleHistoTH2ForGraphs(histo2DRpPbEta, "#it{p}_{T} (GeV/#it{c})","#it{Q}_{pA}", 0.035,0.04, 0.035,0.04, 0.8,1., 512, 505);
    histo2DRpPbEta->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRpPbEta->DrawCopy();

    TLatex *labelALICERpPb  = new TLatex(0.12,0.95-0.04*2.1,textALICE.Data());
    SetStyleTLatex( labelALICERpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
    TLatex *labelEnergyRpPb     = new TLatex(0.12, 0.95-0.04*1, collisionSystempPb.Data());
    SetStyleTLatex( labelEnergyRpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
    TLatex *labelPi0RpPb  = new TLatex(0.12,0.95-0.04*3,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelPi0RpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
    TLatex *labelEtaRpPb  = new TLatex(0.12,0.95-0.04*3,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelEtaRpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);


    for (Int_t est = 0; est < 3; est++){

        Int_t nCurrEst = 0;
        for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
            if ( (centArray2[cent].Contains(nameCentEst[est].Data()) || centArray2[cent].CompareTo("0-100%") == 0) && cent != 5 && enableCentComb[cent])
                nCurrEst++;
        }


        canvasRpPb->cd();
        histo2DRpPb->Draw("copy");
        TLegend* legendPi0RpAPaper    = GetAndSetLegend2(0.35, 0.10, 0.95, 0.10+0.04*1.25*(nCurrEst+1)/3, textSizeLabelsPixel, 3, "", 43, 0.25);
            Int_t plot= 0;
            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if (!(centArray2[cent].Contains(nameCentEst[est].Data()) || centArray2[cent].CompareTo("0-100%") == 0)) continue;
                if (!enableCentComb[cent] || !enableCentRpPb[cent] || cent == 1) continue;
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystPi0[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRpPbCombSystPi0[cent]->Draw("E2same");
                legendPi0RpAPaper->AddEntry(graphRpPbCombSystPi0[cent],centArray[cent].Data(),"pf");
                TBox* boxErrorNormRpAPi0          = CreateBoxConv(colorCentBox[cent], 0.24+plot*0.01, 1.-(nCollErrpPb[cent]/nCollpPb[cent]), 0.25+plot*0.01, 1.+(nCollErrpPb[cent]/nCollpPb[cent]));
                boxErrorNormRpAPi0->Draw();
                plot++;
            }
            DrawGammaLines(minPtPi0Plotting,maxPtPi0Plotting , 1, 1 ,1, kGray, 7);
            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if (!(centArray2[cent].Contains(nameCentEst[est].Data()) || centArray2[cent].CompareTo("0-100%") == 0)) continue;
                if (!enableCentComb[cent] || !enableCentRpPb[cent] || cent == 1) continue;
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatPi0WOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent]);
                graphRpPbCombStatPi0WOXErr[cent]->Draw("p,same,z");
            }

            labelALICERpPb->Draw();
            labelEnergyRpPb->SetText(0.12, 0.95-0.04*1,Form("%s %s", nameCentEstRatios[est].Data(), collisionSystempPb.Data()));
            labelEnergyRpPb->Draw();
            labelPi0RpPb->Draw();
            legendPi0RpAPaper->Draw();
        canvasRpPb->SaveAs(Form("%s/Pi0_QpPb_AllCent_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));


        canvasRpPb->cd();
        histo2DRpPb->Draw("copy");
            plot= 0;
            TBox* boxErrorNormRpAPi02          = CreateBoxConv(colorCentBox[0], 0.24+plot*0.01, 1.-(nCollErrpPb[0]/nCollpPb[0]), 0.25+plot*0.01, 1.+(nCollErrpPb[0]/nCollpPb[0]));
            boxErrorNormRpAPi02->Draw();
            plot++;
            DrawGammaSetMarkerTGraphAsym(graphCombPi0RpPbOlderMBSys, markerStyleCent[0], markerSizeCent[0]*0.5, colorCent[0], colorCent[0], widthLinesBoxes, kTRUE);
            graphCombPi0RpPbOlderMBSys->Draw("E2same");


            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if (!(centArray2[cent].Contains(nameCentEst[est].Data()) )) continue;
                if (!enableCentComb[cent] || !enableCentRpPb[cent] ) continue;
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystPi0[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRpPbCombSystPi0[cent]->Draw("E2same");
                TBox* boxErrorNormRpAPi0          = CreateBoxConv(colorCentBox[cent], 0.24+plot*0.01, 1.-(nCollErrpPb[cent]/nCollpPb[cent]), 0.25+plot*0.01, 1.+(nCollErrpPb[cent]/nCollpPb[cent]));
                boxErrorNormRpAPi0->Draw();
                plot++;
            }
            DrawGammaLines(minPtPi0Plotting,maxPtPi0Plotting , 1, 1 ,1, kGray, 7);
            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if (!(centArray2[cent].Contains(nameCentEst[est].Data()) )) continue;
                if (!enableCentComb[cent] || !enableCentRpPb[cent] ) continue;
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatPi0WOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent]);
                graphRpPbCombStatPi0WOXErr[cent]->Draw("p,same,z");
            }

            ProduceGraphAsymmWithoutXErrors(graphCombPi0RpPbOlderMBStat);
            DrawGammaSetMarkerTGraphAsym(graphCombPi0RpPbOlderMBStat, markerStyleCent[0], markerSizeCent[0]*0.5, colorCent[0], colorCent[0]);
            graphCombPi0RpPbOlderMBStat->Draw("p,same,z");


            labelALICERpPb->Draw();
            labelEnergyRpPb->Draw();
            labelPi0RpPb->Draw();
            legendPi0RpAPaper->Draw();
        canvasRpPb->SaveAs(Form("%s/Pi0_QpPb_AllCent_oldMB_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));

        canvasRpPb->cd();
        histo2DRpPbEta->Draw("copy");
        TLegend* legendEtaRpAPaper    = GetAndSetLegend2(0.3, 0.10, 0.9, 0.10+0.04*1.25*(nCurrEst+1)/3, textSizeLabelsPixel, 3, "", 43, 0.25);

            plot= 0;
            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if (!(centArray2[cent].Contains(nameCentEst[est].Data()) || centArray2[cent].CompareTo("0-100%") == 0)) continue;
                if (!enableCentComb[cent] || !enableCentRpPb[cent] ||  cent == 1) continue;
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystEta[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRpPbCombSystEta[cent]->Draw("E2same");
                legendEtaRpAPaper->AddEntry(graphRpPbCombSystEta[cent],centArray[cent].Data(),"pf");
                TBox* boxErrorNormRpAEta          = CreateBoxConv(colorCentBox[cent], 0.45+plot*0.01, 1.-(nCollErrpPb[cent]/nCollpPb[cent]), 0.46+plot*0.01, 1.+(nCollErrpPb[cent]/nCollpPb[cent]));
                boxErrorNormRpAEta->Draw();
                plot++;
            }
            DrawGammaLines(minPtEtaPlotting,maxPtEtaPlotting , 1, 1 ,1, kGray, 7);

            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if (!(centArray2[cent].Contains(nameCentEst[est].Data()) || centArray2[cent].CompareTo("0-100%") == 0)) continue;
                if (!enableCentComb[cent] || !enableCentRpPb[cent] ||  cent == 1) continue;
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatEtaWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent]);
                graphRpPbCombStatEtaWOXErr[cent]->Draw("p,same,z");
            }

            labelALICERpPb->Draw();
            labelEnergyRpPb->Draw();
            labelEtaRpPb->Draw();
            legendEtaRpAPaper->Draw();
        canvasRpPb->SaveAs(Form("%s/Eta_QpPb_AllCent_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));

        histo2DRpPbEta->Draw("copy");

            plot= 0;

            TBox* boxErrorNormRpAEta2          = CreateBoxConv(colorCentBox[0], 0.45+plot*0.01, 1.-(nCollErrpPb[0]/nCollpPb[0]), 0.46+plot*0.01, 1.+(nCollErrpPb[0]/nCollpPb[0]));
            boxErrorNormRpAEta2->Draw();
            plot++;
            DrawGammaSetMarkerTGraphAsym(graphCombEtaRpPbOlderMBSys, markerStyleCent[0], markerSizeCent[0]*0.5, colorCent[0], colorCent[0], widthLinesBoxes, kTRUE);
            graphCombEtaRpPbOlderMBSys->Draw("E2same");


            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if (!(centArray2[cent].Contains(nameCentEst[est].Data()) )) continue;
                if (!enableCentComb[cent] || !enableCentRpPb[cent] ) continue;
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystEta[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRpPbCombSystEta[cent]->Draw("E2same");
                TBox* boxErrorNormRpAEta          = CreateBoxConv(colorCentBox[cent], 0.45+plot*0.01, 1.-(nCollErrpPb[cent]/nCollpPb[cent]), 0.46+plot*0.01, 1.+(nCollErrpPb[cent]/nCollpPb[cent]));
                boxErrorNormRpAEta->Draw();
                plot++;
            }
            DrawGammaLines(minPtEtaPlotting,maxPtEtaPlotting , 1, 1 ,1, kGray, 7);
            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if (!(centArray2[cent].Contains(nameCentEst[est].Data()) )) continue;
                if (!enableCentComb[cent] || !enableCentRpPb[cent] ) continue;
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatEtaWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent]);
                graphRpPbCombStatEtaWOXErr[cent]->Draw("p,same,z");
            }

            ProduceGraphAsymmWithoutXErrors(graphCombEtaRpPbOlderMBStat);
            DrawGammaSetMarkerTGraphAsym(graphCombEtaRpPbOlderMBStat, markerStyleCent[0], markerSizeCent[0]*0.5, colorCent[0], colorCent[0]);
            graphCombEtaRpPbOlderMBStat->Draw("p,same,z");


        labelALICERpPb->Draw();
        labelEnergyRpPb->Draw();
        labelEtaRpPb->Draw();
        legendEtaRpAPaper->Draw();
        canvasRpPb->SaveAs(Form("%s/Eta_QpPb_AllCent_oldMB_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));
    }


    if (bWCorrection.Contains("Y")){
        for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
            if (!enableCentRpPb[cent]) continue;

            TBox* boxErrorNormRpAPi0          = CreateBoxConv(colorCentBox[cent], 0.24, 1.-(nCollErrpPb[cent]/nCollpPb[cent]), 0.25, 1.+(nCollErrpPb[cent]/nCollpPb[cent]));
            TBox* boxErrorNormRpAEta          = CreateBoxConv(colorCentBox[cent], 0.45, 1.-(nCollErrpPb[cent]/nCollpPb[cent]), 0.46, 1.+(nCollErrpPb[cent]/nCollpPb[cent]));

            if (!centArray[cent].CompareTo("0-100%")){
                histo2DRpPb->GetYaxis()->SetTitle("#it{R}_{pA}");
                histo2DRpPbEta->GetYaxis()->SetTitle("#it{R}_{pA}");
            } else {
                histo2DRpPb->GetYaxis()->SetTitle("#it{Q}_{pA}");
                histo2DRpPbEta->GetYaxis()->SetTitle("#it{Q}_{pA}");
            }
            canvasRpPb->cd();
            histo2DRpPb->DrawCopy();
            boxErrorNormRpAPi0->Draw();

            if (graphRpPbCombSystPi0[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystPi0[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRpPbCombSystPi0[cent]->Draw("E2same");
            }
            if (graphRpPbCombSystEta[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystEta[cent], markerStyleCentMC[cent], markerSizeCentMC[cent]*0.5, colorCentMC[cent], colorCentMC[cent], widthLinesBoxes, kTRUE);
                graphRpPbCombSystEta[cent]->Draw("E2same");
            }
            DrawGammaLines(minPtPi0Plotting,maxPtPi0Plotting , 1, 1 ,1, kGray, 7);
            if (graphRpPbCombStatPi0WOXErr[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatPi0WOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRpPbCombStatPi0WOXErr[cent]->Draw("p,same,z");
            }

            if (graphRpPbCombStatEtaWOXErr[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatEtaWOXErr[cent], markerStyleCentMC[cent], markerSizeCentMC[cent]*0.5, colorCentMC[cent], colorCentMC[cent], widthLinesBoxes, kTRUE);
                graphRpPbCombStatEtaWOXErr[cent]->Draw("p,same,z");
            }

            labelEnergyRpPb->SetText(0.12, 0.95-0.04*1,Form("%s %s", centArray3[cent].Data(), collisionSystempPb.Data()));
            labelEnergyRpPb->Draw();
            labelALICERpPb->Draw();

            TLegend* legendRpPbComb     = GetAndSetLegend2(0.12, 0.95-0.04*1.25*4, 0.32 , 0.95-0.04*1.25*2, textSizeLabelsPixel*0.85,1,"",43,0.3);
            legendRpPbComb->AddEntry(graphRpPbCombSystPi0[cent],"#pi^{0}","pf");
            legendRpPbComb->AddEntry(graphRpPbCombSystEta[cent],"#eta","pf");
            legendRpPbComb->Draw();

            canvasRpPb->Update();
            canvasRpPb->Print(Form("%s/Pi0AndEta_RpPb_%s%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));
            delete legendRpPbComb;

            boxErrorNormRpAPi0->Draw();
//             if (graphDRpASys[cent]){
//                 DrawGammaSetMarkerTGraphAsym(graphDRpASys[cent], 24, markerSizeCentMC[cent]*0.5, kGray+2, kGray+2, widthLinesBoxes, kTRUE);
//                 graphDRpASys[cent]->Draw("E2same");
//             }
//             if (graphChHadRpASys[cent]){
//                 DrawGammaSetMarkerTGraphAsym(graphChHadRpASys[cent], 25, markerSizeCentMC[cent]*0.5, kGray+1, kGray+1, widthLinesBoxes, kTRUE);
//                 graphChHadRpASys[cent]->Draw("E2same");
//             }

            if (graphRpPbCombSystPi0[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystPi0[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRpPbCombSystPi0[cent]->Draw("E2same");
            }
            if (graphRpPbCombSystEta[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystEta[cent], markerStyleCentMC[cent], markerSizeCentMC[cent]*0.5, colorCentMC[cent], colorCentMC[cent], widthLinesBoxes, kTRUE);
                graphRpPbCombSystEta[cent]->Draw("E2same");
            }
            DrawGammaLines(minPtPi0Plotting,maxPtPi0Plotting , 1, 1 ,1, kGray, 7);
            if (graphChHadRpATot[cent]){
                gStyle->SetEndErrorSize(2);
                DrawGammaSetMarkerTGraphAsym(graphChHadRpATot[cent], 21, markerSizeCentMC[cent]*0.25, kGray+1, kGray+1, widthLinesBoxes, kTRUE);
                graphChHadRpATot[cent]->Draw("p,same");
                gStyle->SetEndErrorSize(0);
            }
            if (graphRpPbCombStatPi0WOXErr[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatPi0WOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRpPbCombStatPi0WOXErr[cent]->Draw("p,same,z");
            }
            if (graphDRpATot[cent]){
                gStyle->SetEndErrorSize(2);
                DrawGammaSetMarkerTGraphAsym(graphDRpATot[cent], 20, markerSizeCentMC[cent]*0.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
                graphDRpATot[cent]->Draw("p,same");
                gStyle->SetEndErrorSize(0);
            }
            if (graphRpPbCombStatEtaWOXErr[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatEtaWOXErr[cent], markerStyleCentMC[cent], markerSizeCentMC[cent]*0.5, colorCentMC[cent], colorCentMC[cent], widthLinesBoxes, kTRUE);
                graphRpPbCombStatEtaWOXErr[cent]->Draw("p,same,z");
            }



            labelEnergyRpPb->SetText(0.12, 0.95-0.04*1,Form("%s %s", centArray3[cent].Data(), collisionSystempPb.Data()));
            labelEnergyRpPb->Draw();
            labelALICERpPb->Draw();

            TLegend* legendRpPbComb2     = GetAndSetLegend2(0.12, 0.95-0.04*1.25*4, 0.42 , 0.95-0.04*1.25*2, textSizeLabelsPixel*0.85,2,"",43,0.3);
            legendRpPbComb2->AddEntry(graphRpPbCombSystPi0[cent],"#pi^{0}","pf");
            legendRpPbComb2->AddEntry(graphRpPbCombSystEta[cent],"#eta","pf");
            if (graphDRpATot[cent])legendRpPbComb2->AddEntry(graphDRpATot[cent],"D^{0}","pe");
            if (graphChHadRpATot[cent])legendRpPbComb2->AddEntry(graphChHadRpATot[cent],"h^{#pm}","pe");
            legendRpPbComb2->Draw();

            canvasRpPb->Update();
            canvasRpPb->Print(Form("%s/Pi0AndEtaPlusOther_RpPb_%s%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));
            delete legendRpPbComb2;

            //*****************************************************************************************************************
            // Plotting RpA pi0 diff rec techniques
            //*****************************************************************************************************************
            histo2DRpPb->DrawCopy();

            boxErrorNormRpAPi0->Draw();
            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRpPbIndSystPi0[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRpPbIndSystPi0[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth], widthLinesBoxes, kTRUE);
                    graphRpPbIndSystPi0[cent][meth]->Draw("E2same");
                }
            }

            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRpPbIndStatPi0WOXErr[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRpPbIndStatPi0WOXErr[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth]);
                    graphRpPbIndStatPi0WOXErr[cent][meth]->Draw("p,same,z");
                }
            }

            DrawGammaLines(minPtPi0Plotting,maxPtPi0Plotting , 1, 1 ,1, kGray, 7);

            labelEnergyRpPb->Draw();
            labelALICERpPb->Draw();
            labelPi0RpPb->Draw();

            TLegend* legendRpPbInd     = GetAndSetLegend2(0.45, 0.13, 0.95 , 0.13+0.04*3, textSizeLabelsPixel*0.85,2, "", 43, 0.2);
            for (Int_t meth = 0; meth< 11; meth++){
                if (graphRpPbIndSystPi0[cent][meth])legendRpPbInd->AddEntry(graphRpPbIndSystPi0[cent][meth],nameMeasGlobalLabel[meth].Data(),"pf");
            }
            legendRpPbInd->Draw();

            canvasRpPb->Update();
            canvasRpPb->Print(Form("%s/Pi0_RpPb_IndividualMeasurements_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));
            delete legendRpPbInd;

            //*****************************************************************************************************************
            // Plotting RpA eta diff rec techniques
            //*****************************************************************************************************************
            histo2DRpPbEta->DrawCopy();
            boxErrorNormRpAEta->Draw();
            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRpPbIndSystEta[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRpPbIndSystEta[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth], widthLinesBoxes, kTRUE);
                    graphRpPbIndSystEta[cent][meth]->Draw("E2same");
                }
            }
            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRpPbIndStatEtaWOXErr[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRpPbIndStatEtaWOXErr[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth]);
                    graphRpPbIndStatEtaWOXErr[cent][meth]->Draw("p,same,z");
                }
            }

            DrawGammaLines(minPtEtaPlotting,maxPtEtaPlotting , 1, 1 ,1, kGray, 7);

            labelEnergyRpPb->Draw();
            labelALICERpPb->Draw();
            labelEtaRpPb->Draw();

            TLegend* legendRpPbIndEta     = GetAndSetLegend2(0.45, 0.13, 0.95 , 0.13+0.04*3, textSizeLabelsPixel*0.85,2, "", 43, 0.2);
            for (Int_t meth = 0; meth< 11; meth++){
                if (graphRpPbIndSystEta[cent][meth])legendRpPbIndEta->AddEntry(graphRpPbIndSystEta[cent][meth],nameMeasGlobalLabel[meth].Data(),"pf");
            }
            legendRpPbIndEta->Draw();

            canvasRpPb->Update();
            canvasRpPb->Print(Form("%s/Eta_RpPb_IndividualMeasurements_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));
            delete legendRpPbIndEta;
        }
    }

    //*****************************************************************************************************************
    // Plotting RCP s for different mesons
    //*****************************************************************************************************************
    textSizeLabelsPixel                         = 54;
    TCanvas* canvasRCP = new TCanvas("canvasRCP","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasRCP,  0.085, 0.01, 0.015, 0.08);
    canvasRCP->SetLogx();
    TH2F * histo2DRCP  = new TH2F("histo2DRCP","histo2DRCP",1000,minPtPi0Plotting,maxPtPi0Plotting,1000,0.0,4.5);
    SetStyleHistoTH2ForGraphs(histo2DRCP, "#it{p}_{T} (GeV/#it{c})","#it{Q}_{CP}", 0.035,0.04, 0.035,0.04, 0.8,1., 512, 505);
    histo2DRCP->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRCP->DrawCopy();
    TH2F * histo2DRCPEta  = new TH2F("histo2DRCPEta","histo2DRCPEta",1000,minPtEtaPlotting,maxPtEtaPlotting,1000,0.0,4.5);
    SetStyleHistoTH2ForGraphs(histo2DRCPEta, "#it{p}_{T} (GeV/#it{c})","#it{Q}_{CP}", 0.035,0.04, 0.035,0.04, 0.8,1., 512, 505);
    histo2DRCPEta->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRCPEta->DrawCopy();

    TLatex *labelALICERCP  = new TLatex(0.12,0.95-0.04*2.1,textALICE.Data());
    SetStyleTLatex( labelALICERCP, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
    TLatex *labelEnergyRCP     = new TLatex(0.12, 0.95-0.04*1, collisionSystempPb.Data());
    SetStyleTLatex( labelEnergyRCP, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
    TLatex *labelPi0RCP  = new TLatex(0.12,0.95-0.04*3,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelPi0RCP, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
    TLatex *labelEtaRCP  = new TLatex(0.12,0.95-0.04*3,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelEtaRCP, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);

    canvasRCP->cd();


    if (bWCorrection.Contains("Y")){

        for (Int_t est = 0; est < 3; est++){

            Int_t nCurrEst = 0;
            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if ( (centArray2[cent].Contains(nameCentEst[est].Data()) ) && cent != 5 && enableCentComb[cent] && enableCentRCP[cent])
                    nCurrEst++;
            }
            if (nCurrEst ==  0) continue;

            histo2DRCP->Draw("copy");
            TLegend* legendPi0RCPPaper    = GetAndSetLegend2(0.35, 0.10, 0.95, 0.10+0.04*1.25*(nCurrEst)/3, textSizeLabelsPixel, 3, "", 43, 0.25);

            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if ( !centArray2[cent].Contains(nameCentEst[est].Data()) ) continue;
                if (!enableCentComb[cent] || !enableCentRCP[cent] || cent == 1) continue;
                DrawGammaSetMarkerTGraphAsym(graphRCPCombSystPi0[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRCPCombSystPi0[cent]->Draw("E2same");
                legendPi0RCPPaper->AddEntry(graphRCPCombSystPi0[cent],centArray[cent].Data(),"pf");
    //             DrawGammaLines(10, maxPtPi0Plotting , rCPNColl[cent], rCPNColl[cent],2, colorCentMC[cent]);
            }
            DrawGammaLines(minPtPi0Plotting,maxPtPi0Plotting , 1, 1 ,1, kGray, 7);

            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if ( !centArray2[cent].Contains(nameCentEst[est].Data()) ) continue;
                if (!enableCentComb[cent] || !enableCentRCP[cent] || cent == 1) continue;
                DrawGammaSetMarkerTGraphAsym(graphRCPCombStatPi0WOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent]);
                graphRCPCombStatPi0WOXErr[cent]->Draw("p,same,z");
            }

            labelALICERCP->Draw();
            labelEnergyRCP->SetText(0.12, 0.95-0.04*1,Form("%s %s", nameCentEstRatios[est].Data(), collisionSystempPb.Data()));
            labelEnergyRCP->Draw();
            labelPi0RCP->Draw();
            legendPi0RCPPaper->Draw();
            canvasRCP->SaveAs(Form("%s/Pi0_QCP_AllCent_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));

            canvasRCP->cd();
            histo2DRCPEta->Draw("copy");
            TLegend* legendEtaRCPPaper    = GetAndSetLegend2(0.3, 0.10, 0.9, 0.10+0.04*1.25*(nCurrEst)/3, textSizeLabelsPixel, 3, "", 43, 0.25);

            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if ( !centArray2[cent].Contains(nameCentEst[est].Data()) ) continue;
                if (!enableCentComb[cent] || !enableCentRCP[cent] ||  cent == 1) continue;
                DrawGammaSetMarkerTGraphAsym(graphRCPCombSystEta[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRCPCombSystEta[cent]->Draw("E2same");
                legendEtaRCPPaper->AddEntry(graphRCPCombSystEta[cent],centArray[cent].Data(),"pf");
    //             DrawGammaLines(10, maxPtPi0Plotting , rCPNColl[cent], rCPNColl[cent],0.1, colorCentMC[cent]);
            }
            DrawGammaLines(minPtEtaPlotting,maxPtEtaPlotting , 1, 1 ,1, kGray, 7);

            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if ( !centArray2[cent].Contains(nameCentEst[est].Data()) ) continue;
                if (!enableCentComb[cent] || !enableCentRCP[cent] ||  cent == 1) continue;
                DrawGammaSetMarkerTGraphAsym(graphRCPCombStatEtaWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent]);
                graphRCPCombStatEtaWOXErr[cent]->Draw("p,same,z");
            }

            labelALICERCP->Draw();
            labelEnergyRCP->Draw();
            labelEtaRCP->Draw();
            legendEtaRCPPaper->Draw();
            canvasRCP->SaveAs(Form("%s/Eta_QCP_AllCent_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));
        }

        for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
            if (!enableCentRCP[cent]) continue;

            Double_t errNormRCP               = TMath::Sqrt(TMath::Power(nCollErrpPb[cent]/nCollpPb[cent],2)+TMath::Power(nCollErrpPb[rCPRefInt[cent]]/nCollpPb[rCPRefInt[cent]],2));
            TBox* boxErrorNormRCPPi0          = CreateBoxConv(colorCentBox[cent], 0.24, 1.-(errNormRCP), 0.25, 1.+(errNormRCP));
            TBox* boxErrorNormRCPEta          = CreateBoxConv(colorCentBox[cent], 0.45, 1.-(errNormRCP), 0.46, 1.+(errNormRCP));

            canvasRCP->cd();
            histo2DRCP->DrawCopy();

            boxErrorNormRCPPi0->Draw();
            if (graphRCPCombSystPi0[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRCPCombSystPi0[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRCPCombSystPi0[cent]->Draw("E2same");
            }
            if (graphRCPCombSystEta[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRCPCombSystEta[cent], markerStyleCentMC[cent], markerSizeCentMC[cent]*0.5, colorCentMC[cent], colorCentMC[cent], widthLinesBoxes, kTRUE);
                graphRCPCombSystEta[cent]->Draw("E2same");
            }
            DrawGammaLines(minPtPi0Plotting,maxPtPi0Plotting , 1, 1 ,1, kGray, 7);
            if (graphRCPCombStatPi0WOXErr[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRCPCombStatPi0WOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRCPCombStatPi0WOXErr[cent]->Draw("p,same,z");
            }
            if (graphRCPCombStatEtaWOXErr[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRCPCombStatEtaWOXErr[cent], markerStyleCentMC[cent], markerSizeCentMC[cent]*0.5, colorCentMC[cent], colorCentMC[cent], widthLinesBoxes, kTRUE);
                graphRCPCombStatEtaWOXErr[cent]->Draw("p,same,z");
            }

            labelEnergyRCP->SetText(0.12, 0.95-0.04*1,Form("%s %s", centArray3[cent].Data(), collisionSystempPb.Data()));
            labelEnergyRCP->Draw();
            labelALICERCP->Draw();

            TLegend* legendRCPComb     = GetAndSetLegend2(0.12, 0.95-0.04*1.25*4, 0.32 , 0.95-0.04*1.25*2, textSizeLabelsPixel*0.85,1,"",43,0.3);
            legendRCPComb->AddEntry(graphRCPCombSystPi0[cent],"#pi^{0}","pf");
            legendRCPComb->AddEntry(graphRCPCombSystEta[cent],"#eta","pf");
            legendRCPComb->Draw();
            histo2DRCP->Draw("same,axis");

            canvasRCP->Update();
            canvasRCP->Print(Form("%s/Pi0AndEta_RCP_%s%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));
            delete legendRCPComb;

            //*****************************************************************************************************************
            // Plotting RCP pi0 diff rec techniques
            //*****************************************************************************************************************
            histo2DRCP->DrawCopy();

            boxErrorNormRCPPi0->Draw();

            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRCPIndSystPi0[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRCPIndSystPi0[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth], widthLinesBoxes, kTRUE);
                    graphRCPIndSystPi0[cent][meth]->Draw("E2same");
                }
            }

            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRCPIndStatPi0WOXErr[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRCPIndStatPi0WOXErr[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth]);
                    graphRCPIndStatPi0WOXErr[cent][meth]->Draw("p,same,z");
                }
            }


            labelEnergyRCP->Draw();
            labelALICERCP->Draw();
            labelPi0RCP->Draw();

            TLegend* legendRCPInd     = GetAndSetLegend2(0.45, 0.13, 0.95 , 0.13+0.04*3, textSizeLabelsPixel*0.85,2, "", 43, 0.2);
            for (Int_t meth = 0; meth< 11; meth++){
                if (graphRCPIndSystPi0[cent][meth])legendRCPInd->AddEntry(graphRCPIndSystPi0[cent][meth],nameMeasGlobalLabel[meth].Data(),"pf");
            }
            legendRCPInd->Draw();
            histo2DRCP->Draw("same,axis");

            canvasRCP->Update();
            canvasRCP->Print(Form("%s/Pi0_RCP_IndividualMeasurements_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));
            delete legendRCPInd;

            //*****************************************************************************************************************
            // Plotting RCP eta diff rec techniques
            //*****************************************************************************************************************
            histo2DRCPEta->DrawCopy();
            boxErrorNormRCPEta->Draw();
            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRCPIndSystEta[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRCPIndSystEta[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth], widthLinesBoxes, kTRUE);
                    graphRCPIndSystEta[cent][meth]->Draw("E2same");
                }
            }
            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRCPIndStatEtaWOXErr[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRCPIndStatEtaWOXErr[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth]);
                    graphRCPIndStatEtaWOXErr[cent][meth]->Draw("p,same,z");
                }
            }

            labelEnergyRCP->Draw();
            labelALICERCP->Draw();
            labelEtaRCP->Draw();

            TLegend* legendRCPIndEta     = GetAndSetLegend2(0.45, 0.13, 0.95 , 0.13+0.04*3, textSizeLabelsPixel*0.85,2, "", 43, 0.2);
            for (Int_t meth = 0; meth< 11; meth++){
                if (graphRCPIndSystEta[cent][meth])legendRCPIndEta->AddEntry(graphRCPIndSystEta[cent][meth],nameMeasGlobalLabel[meth].Data(),"pf");
            }
            legendRCPIndEta->Draw();

            canvasRCP->Update();
            canvasRCP->Print(Form("%s/Eta_RCP_IndividualMeasurements_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));
            delete legendRCPIndEta;
        }
    }

    //*****************************************************************************************************************
    // Plotting RMB s for different mesons
    //*****************************************************************************************************************
    textSizeLabelsPixel                         = 54;
    TCanvas* canvasRMB = new TCanvas("canvasRMB","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasRMB,  0.085, 0.01, 0.015, 0.08);
    canvasRMB->SetLogx();
    TH2F * histo2DRMB  = new TH2F("histo2DRMB","histo2DRMB",1000,minPtPi0Plotting,maxPtPi0Plotting,1000,0.0,2.01);
    SetStyleHistoTH2ForGraphs(histo2DRMB, "#it{p}_{T} (GeV/#it{c})","#it{Q}_{MB}", 0.035,0.04, 0.035,0.04, 0.8,1., 512, 505);
    histo2DRMB->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRMB->DrawCopy();
    TH2F * histo2DRMBEta  = new TH2F("histo2DRMBEta","histo2DRMBEta",1000,minPtEtaPlotting,maxPtEtaPlotting,1000,0.0,2.01);
    SetStyleHistoTH2ForGraphs(histo2DRMBEta, "#it{p}_{T} (GeV/#it{c})","#it{Q}_{MB}", 0.035,0.04, 0.035,0.04, 0.8,1., 512, 505);
    histo2DRMBEta->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRMBEta->DrawCopy();

    TLatex *labelALICERMB  = new TLatex(0.12,0.95-0.04*2.1,textALICE.Data());
    SetStyleTLatex( labelALICERMB, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
    TLatex *labelEnergyRMB     = new TLatex(0.12, 0.95-0.04*1, collisionSystempPb.Data());
    SetStyleTLatex( labelEnergyRMB, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
    TLatex *labelPi0RMB  = new TLatex(0.12,0.95-0.04*3,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelPi0RMB, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
    TLatex *labelEtaRMB  = new TLatex(0.12,0.95-0.04*3,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelEtaRMB, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);

    if (bWCorrection.Contains("Y")){

        for (Int_t est = 0; est < 3; est++){

            Int_t nCurrEst = 0;
            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if ( (centArray2[cent].Contains(nameCentEst[est].Data()) ) && cent != 5 && enableCentComb[cent] && enableCentRMB[cent])
                    nCurrEst++;
            }
            if (nCurrEst ==  0) continue;
            canvasRMB->cd();
            histo2DRMB->Draw("copy");
            TLegend* legendPi0RMBPaper    = GetAndSetLegend2(0.45, 0.10, 0.95, 0.10+0.04*1.25*(nCurrEst)/2, textSizeLabelsPixel, 2, "", 43, 0.25);

            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if ( !centArray2[cent].Contains(nameCentEst[est].Data()) ) continue;
                if (!enableCentComb[cent] || !enableCentRMB[cent] || cent == 1) continue;
                DrawGammaSetMarkerTGraphAsym(graphRMBCombSystPi0[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRMBCombSystPi0[cent]->Draw("E2same");
                legendPi0RMBPaper->AddEntry(graphRMBCombSystPi0[cent],centArray[cent].Data(),"pf");
//                 DrawGammaLines(15, maxPtPi0Plotting , rMBNColl[cent], rMBNColl[cent],2, colorCent[cent], 7);
            }
            DrawGammaLines(minPtPi0Plotting,maxPtPi0Plotting , 1, 1 ,1, kGray, 7);
            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if ( !centArray2[cent].Contains(nameCentEst[est].Data()) ) continue;
                if (!enableCentComb[cent] || !enableCentRMB[cent] || cent == 1) continue;
                DrawGammaSetMarkerTGraphAsym(graphRMBCombStatPi0WOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent]);
                graphRMBCombStatPi0WOXErr[cent]->Draw("p,same,z");
            }

            labelALICERMB->Draw();
            labelEnergyRMB->SetText(0.12, 0.95-0.04*1,Form("%s %s", nameCentEstRatios[est].Data(), collisionSystempPb.Data()));
            labelEnergyRMB->Draw();
            labelPi0RMB->Draw();
            legendPi0RMBPaper->Draw();
            canvasRMB->SaveAs(Form("%s/Pi0_QMB_AllCent_%s.%s",outputDir.Data(), nameCentEst[est].Data() , suffix.Data()));

            canvasRMB->cd();
            histo2DRMBEta->Draw("copy");
            TLegend* legendEtaRMBPaper    = GetAndSetLegend2(0.4, 0.10, 0.9, 0.10+0.04*1.25*(nCurrEst)/2, textSizeLabelsPixel, 2, "", 43, 0.25);

            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if ( !centArray2[cent].Contains(nameCentEst[est].Data()) ) continue;
                if (!enableCentComb[cent] || !enableCentRMB[cent] ||  cent == 1) continue;
                DrawGammaSetMarkerTGraphAsym(graphRMBCombSystEta[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRMBCombSystEta[cent]->Draw("E2same");
                legendEtaRMBPaper->AddEntry(graphRMBCombSystEta[cent],centArray[cent].Data(),"pf");
//                 DrawGammaLines(15, maxPtEtaPlotting , rMBNColl[cent], rMBNColl[cent],2, colorCent[cent], 7);
            }
            DrawGammaLines(minPtEtaPlotting,maxPtEtaPlotting , 1, 1 ,1, kGray, 7);
            for (Int_t cent= 0; cent < maxCentRun1; cent++ ){
                if ( !centArray2[cent].Contains(nameCentEst[est].Data()) ) continue;
                if (!enableCentComb[cent] || !enableCentRMB[cent] ||  cent == 1) continue;
                DrawGammaSetMarkerTGraphAsym(graphRMBCombStatEtaWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent]);
                graphRMBCombStatEtaWOXErr[cent]->Draw("p,same,z");
            }

            labelALICERMB->Draw();
            labelEnergyRMB->Draw();
            labelEtaRMB->Draw();
            legendEtaRMBPaper->Draw();
            canvasRMB->SaveAs(Form("%s/Eta_QMB_AllCent_%s.%s",outputDir.Data(), nameCentEst[est].Data(), suffix.Data()));
        }

        for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
            if (!enableCentRMB[cent]) continue;

            Double_t errNormRMB               = TMath::Sqrt(TMath::Power(nCollErrpPb[cent]/nCollpPb[cent],2)+TMath::Power(nCollErrpPb[rMBRefInt[cent]]/nCollpPb[rMBRefInt[cent]],2));
            TBox* boxErrorNormRMBPi0          = CreateBoxConv(colorCentBox[cent], 0.24, 1.-(errNormRMB), 0.25, 1.+(errNormRMB));
            TBox* boxErrorNormRMBEta          = CreateBoxConv(colorCentBox[cent], 0.45, 1.-(errNormRMB), 0.46, 1.+(errNormRMB));

            canvasRMB->cd();
            histo2DRMB->DrawCopy();
            boxErrorNormRMBPi0->Draw();
            if (graphRMBCombSystPi0[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRMBCombSystPi0[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRMBCombSystPi0[cent]->Draw("E2same");
            }
            if (graphRMBCombSystEta[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRMBCombSystEta[cent], markerStyleCentMC[cent], markerSizeCentMC[cent]*0.5, colorCentMC[cent], colorCentMC[cent], widthLinesBoxes, kTRUE);
                graphRMBCombSystEta[cent]->Draw("E2same");
            }
            DrawGammaLines(minPtPi0Plotting,maxPtPi0Plotting , 1, 1 ,1, kGray, 7);
            if (graphRMBCombStatPi0WOXErr[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRMBCombStatPi0WOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRMBCombStatPi0WOXErr[cent]->Draw("p,same,z");
            }

            if (graphRMBCombStatEtaWOXErr[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRMBCombStatEtaWOXErr[cent], markerStyleCentMC[cent], markerSizeCentMC[cent]*0.5, colorCentMC[cent], colorCentMC[cent], widthLinesBoxes, kTRUE);
                graphRMBCombStatEtaWOXErr[cent]->Draw("p,same,z");
            }

            labelEnergyRMB->SetText(0.12, 0.95-0.04*1,Form("%s %s", centArray3[cent].Data(), collisionSystempPb.Data()));
            labelEnergyRMB->Draw();
            labelALICERMB->Draw();

            TLegend* legendRMBComb     = GetAndSetLegend2(0.12, 0.95-0.04*1.25*4, 0.32 , 0.95-0.04*1.25*2, textSizeLabelsPixel*0.85,1,"",43,0.3);
            legendRMBComb->AddEntry(graphRMBCombSystPi0[cent],"#pi^{0}","pf");
            legendRMBComb->AddEntry(graphRMBCombSystEta[cent],"#eta","pf");
            legendRMBComb->Draw();

            canvasRMB->Update();
            canvasRMB->Print(Form("%s/Pi0AndEta_RMB_%s%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));
            delete legendRMBComb;

            //*****************************************************************************************************************
            // Plotting RMB pi0 diff rec techniques
            //*****************************************************************************************************************
            histo2DRMB->DrawCopy();

            boxErrorNormRMBPi0->Draw();
            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRMBIndSystPi0[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRMBIndSystPi0[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth], widthLinesBoxes, kTRUE);
                    graphRMBIndSystPi0[cent][meth]->Draw("E2same");
                }
            }

            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRMBIndStatPi0WOXErr[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRMBIndStatPi0WOXErr[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth]);
                    graphRMBIndStatPi0WOXErr[cent][meth]->Draw("p,same,z");
                }
            }

            labelEnergyRMB->Draw();
            labelALICERMB->Draw();
            labelPi0RMB->Draw();

            TLegend* legendRMBInd     = GetAndSetLegend2(0.45, 0.13, 0.95 , 0.13+0.04*3, textSizeLabelsPixel*0.85,2, "", 43, 0.2);
            for (Int_t meth = 0; meth< 11; meth++){
                if (graphRMBIndSystPi0[cent][meth])legendRMBInd->AddEntry(graphRMBIndSystPi0[cent][meth],nameMeasGlobalLabel[meth].Data(),"pf");
            }
            legendRMBInd->Draw();

            canvasRMB->Update();
            canvasRMB->Print(Form("%s/Pi0_RMB_IndividualMeasurements_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));
            delete legendRMBInd;

            //*****************************************************************************************************************
            // Plotting RMB eta diff rec techniques
            //*****************************************************************************************************************
            histo2DRMBEta->DrawCopy();
            boxErrorNormRMBEta->Draw();
            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRMBIndSystEta[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRMBIndSystEta[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth], widthLinesBoxes, kTRUE);
                    graphRMBIndSystEta[cent][meth]->Draw("E2same");
                }
            }
            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRMBIndStatEtaWOXErr[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRMBIndStatEtaWOXErr[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth]);
                    graphRMBIndStatEtaWOXErr[cent][meth]->Draw("p,same,z");
                }
            }
            labelEnergyRMB->Draw();
            labelALICERMB->Draw();
            labelEtaRMB->Draw();

            TLegend* legendRMBIndEta     = GetAndSetLegend2(0.45, 0.13, 0.95 , 0.13+0.04*3, textSizeLabelsPixel*0.85,2, "", 43, 0.2);
            for (Int_t meth = 0; meth< 11; meth++){
                if (graphRMBIndSystEta[cent][meth])legendRMBIndEta->AddEntry(graphRMBIndSystEta[cent][meth],nameMeasGlobalLabel[meth].Data(),"pf");
            }
            legendRMBIndEta->Draw();

            canvasRMB->Update();
            canvasRMB->Print(Form("%s/Eta_RMB_IndividualMeasurements_%s%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), runArray[cent].Data(), suffix.Data()));
            delete legendRMBIndEta;
        }
    }

    // **********************************************************************************************************************
    // ************************* Saving of final results ********************************************************************
    // **********************************************************************************************************************

    TString nameOutputCommonFile    = Form("CombinedResultsPaperpPb5TeVCent_%s.root", dateForOutput.Data());

    TFile fCombResults(nameOutputCommonFile.Data(), "UPDATE");

    for (Int_t cent = 0; cent < maxCentRun1+maxCentRun2; cent++){
        if (!enableCent[cent]) continue;
        TDirectoryFile* directoryPi0Output  = NULL;
        directoryPi0Output                  = (TDirectoryFile*)fCombResults.Get(Form("Pi0pPb5TeV_%s%s%s", centArray[cent].Data(), addCentString[cent].Data(), runArray[cent].Data()));
        if (!directoryPi0Output){
            fCombResults.mkdir(Form("Pi0pPb5TeV_%s%s%s", centArray[cent].Data(), addCentString[cent].Data(), runArray[cent].Data()));
            directoryPi0Output              = (TDirectoryFile*)fCombResults.Get(Form("Pi0pPb5TeV_%s%s%s", centArray[cent].Data(), addCentString[cent].Data(), runArray[cent].Data()));
        }
        fCombResults.cd(Form("Pi0pPb5TeV_%s%s%s", centArray[cent].Data(), addCentString[cent].Data(), runArray[cent].Data()));

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
        directoryEtaOutput                  = (TDirectoryFile*)fCombResults.Get(Form("EtapPb5TeV_%s%s%s", centArray[cent].Data(), addCentString[cent].Data(), runArray[cent].Data()));
        if (!directoryEtaOutput){
            fCombResults.mkdir(Form("EtapPb5TeV_%s%s%s", centArray[cent].Data(), addCentString[cent].Data(), runArray[cent].Data()));
            directoryEtaOutput              = (TDirectoryFile*)fCombResults.Get(Form("EtapPb5TeV_%s%s%s", centArray[cent].Data(), addCentString[cent].Data(), runArray[cent].Data()));
        }
        fCombResults.cd(Form("EtapPb5TeV_%s%s%s", centArray[cent].Data(), addCentString[cent].Data(), runArray[cent].Data()));

            // Final spectrum
            if (graphCombEtaInvYieldTot[cent])  graphCombEtaInvYieldTot[cent]->Write("graphInvYieldEtaCombTotErr",TObject::kOverwrite);
            if (graphCombEtaInvYieldStat[cent]) graphCombEtaInvYieldStat[cent]->Write("graphInvYieldEtaCombStatErr",TObject::kOverwrite);
            if (graphCombEtaInvYieldSys[cent])  graphCombEtaInvYieldSys[cent]->Write("graphInvYieldEtaCombSysErr",TObject::kOverwrite);
                // fits for pi0
            if (fitInvYieldEta[cent])           fitInvYieldEta[cent]->Write("TsallisFitEta",TObject::kOverwrite);
            if (fitTCMInvYieldEta[cent])        fitTCMInvYieldEta[cent]->Write("TwoComponentModelFitEta",TObject::kOverwrite);

            if (graphCombEtaToPi0Tot[cent])     graphCombEtaToPi0Sys[cent]->Write("graphEtaToPi0CombTotErr",TObject::kOverwrite);
            if (graphCombEtaToPi0Stat[cent])    graphCombEtaToPi0Stat[cent]->Write("graphEtaToPi0CombStatErr",TObject::kOverwrite);
            if (graphCombEtaToPi0Sys[cent])     graphCombEtaToPi0Sys[cent]->Write("graphEtaToPi0CombSysErr",TObject::kOverwrite);


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
                if (histoEtaToPi0Stat[cent][meth]) histoEtaToPi0Stat[cent][meth]->Write(Form("histoEtaToPi0%sStatErr", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
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
            fitInvYieldPi0[cent]->Write(Form("TsallisFitPi0_%s%s%s",centArray[cent].Data(), addCentString[cent].Data(), runArray[cent].Data()));
            fitTCMInvYieldPi0[cent]->Write(Form("TwoComponentModelFitPi0_%s%s%s",centArray[cent].Data(), addCentString[cent].Data(), runArray[cent].Data()));
            // fits for eta
            fitInvYieldEta[cent]->Write(Form("TsallisFitEta_%s%s%s",centArray[cent].Data(), addCentString[cent].Data(), runArray[cent].Data()));
            fitTCMInvYieldEta[cent]->Write(Form("TwoComponentModelFitEta_%s%s%s",centArray[cent].Data(), addCentString[cent].Data(), runArray[cent].Data()));
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

