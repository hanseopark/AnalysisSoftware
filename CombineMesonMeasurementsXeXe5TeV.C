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

void CombineMesonMeasurementsXeXe5TeV(  TString fileNamePCM             = "",
                                        TString fileNameEMCAL           = "",
                                        TString fileNamePHOS            = "",
                                        TString fileNamePCMEMCAL        = "",
                                        TString fileNamePCMPHOS         = "",
                                        TString suffix                  = "eps",
                                        Bool_t isMC                     = kFALSE,
                                        TString bWCorrection            = "X",
                                        TString fileNameCorrFactors     = "",
                                        TString fileNameInterpolation   = "",
                                        TString fileConfigRXeXeErr      = ""
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
    TString collisionSystemXeXe                 = "Xe-Xe, #sqrt{#it{s}_{_{NN}}} = 5.44 TeV";
    TString collisionSystempp                   = "pp, #sqrt{#it{s}} = 5.02 TeV";
    TString outputDir                           = Form("%s/%s/CombineMesonMeasurements%sXeXe5TeV", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());

    TString nameFinalResDat                     = Form("%s/CombinedResults%s_FitResults.dat", outputDir.Data(), bWCorrection.Data());
    fstream  fileFitsOutput;
    fileFitsOutput.open(nameFinalResDat.Data(), ios::out);

    TString fileNameTheory                      = "ExternalInputPbPb/Theory/TheoryCompilationPbPb.root";

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMEMCAL.root", fileNamePCMEMCAL.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMPHOS.root", fileNamePCMPHOS.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCAL.root", fileNameEMCAL.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputInterpolation.root", fileNameInterpolation.Data(), outputDir.Data()));

    TString prefix2                             = "";
    if (isMC){
        prefix2                                 = "MC";
    } else {
        prefix2                                 = "Data";
    }

    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();

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
    TString  centArray[5]                       = {"0-20%", "20-40%", "40-80%", "0-40%", "0-80%"};
    TString  centArrayOutput[5]                 = {"0020", "2040", "4080", "0040", "0080"};

    TString  nameMeasGlobal[11]                 = { "PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal",
                                                    "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged",
                                                    "PCMOtherDataset"};
    TString  nameMeasGlobalLabel[11]            = { "PCM", "PHOS", "EMC", "PCM-PHOS", "PCM-EMC",
                                                    "PCM-Dal", "PHOS-Dal", "EMC-Dal", "EMChigh", "mEMC",
                                                    "PCMOtherDataset"};
    TString  nameSecPi0SourceRead[4]            = {"K0S", "K0L", "Lambda", "Rest"};
    TString  nameSecPi0SourceLabel[4]           = {"K^{0}_{s}", "K^{0}_{l}", "#Lambda", "had. int."};
    Double_t maxSecCorr[4]                      = { 0.05, 0.007, 0.0003, 0.03};

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

    Color_t  colorCent[5];
    Color_t  colorCentMC[5];
    Style_t  markerStyleCent[5];
    Style_t  markerStyleCentMC[5];
    Size_t   markerSizeCent[5];
    Size_t   markerSizeCentMC[5];
    Double_t nCollXeXe[5];
    Double_t nCollErrXeXe[5];
    Double_t tXeXe[5];
    Double_t tXeXeErr[5];
    for (Int_t cent = 0; cent< 5; cent++){
        colorCent[cent]                         = GetColorDefaultColor("XeXe_5.44TeV", "", centArray[cent]);
        colorCentMC[cent]                       = GetColorDefaultColor("XeXe_5.44TeV", "HIJING", centArray[cent]);
        markerStyleCent[cent]                   = GetDefaultMarkerStyle("XeXe_5.44TeV", "", centArray[cent]);
        markerStyleCentMC[cent]                 = GetDefaultMarkerStyle("XeXe_5.44TeV", "HIJING", centArray[cent]);
        markerSizeCent[cent]                    = GetDefaultMarkerSize("XeXe_5.44TeV", "", centArray[cent])*2;
        markerSizeCentMC[cent]                  = GetDefaultMarkerSize("XeXe_5.44TeV", "HIJING", centArray[cent])*2;
        nCollXeXe[cent]                         = GetNCollFromName(centArray[cent], "XeXe_5.44TeV");
        nCollErrXeXe[cent]                      = GetNCollErrFromName(centArray[cent], "XeXe_5.44TeV");
        tXeXe[cent]                             = GetTAAFromName(centArray[cent], "XeXe_5.44TeV");
        tXeXeErr[cent]                          = GetTAAErrFromName(centArray[cent], "XeXe_5.44TeV");
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
    Double_t xPtLimitsPi0[5][100];
    Int_t maxNBinsPi0Abs                        = 0;
    Int_t maxNBinsPi0[5]                        = {0};

    Int_t nTotMeas[5]                               = {0};
    TString fileNamesMethod[11]                     = {"", "", "", "", "", "", "", "", "", "", ""};
    TString fileNamesXeXePi0DetailedSys[11]         = {"", "", "", "", "", "", "", "", "", "", ""};
    TString fileNamesppPi0DetailedSys[11]           = {"", "", "", "", "", "", "", "", "", "", ""};
    TString fileNamesRXeXePi0DetailedSys[11]        = {"", "", "", "", "", "", "", "", "", "", ""};
    Bool_t havePi0SysDetailedXeXe[11]               = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
    Bool_t havePi0SysDetailedpp[11]                 = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
    Bool_t haveEffSecCorr[5][4][11];
    vector<TString>* ptSysRemNames                  = new vector<TString>[11];
    TFile* fileMethod[11];
    TDirectory* directoryPi0[5][11];
    TGraphAsymmErrors* graphPi0EffSecCorrFromX[5][4][11];
    TGraphAsymmErrors* graphPi0Eff[5][11];
    TGraphAsymmErrors* graphPi0Acc[5][11];
    TGraphAsymmErrors* graphPi0EffTimesAcc[5][11];
    TGraphAsymmErrors* graphPi0Mass[5][11];
    TGraphAsymmErrors* graphPi0MassMC[5][11];
    TGraphAsymmErrors* graphPi0Width[5][11];
    TGraphAsymmErrors* graphPi0WidthMC[5][11];
    TH1D* histoPi0InvYieldStat[5][11];
    TGraphAsymmErrors* graphPi0InvYieldStat[5][11];
    TGraphAsymmErrors* graphPi0InvYieldSys[5][11];

    for (Int_t cent = 0; cent < 5; cent++){
        maxNBinsPi0[cent]                           = GetBinning( xPtLimitsPi0[cent], maxNBinsPi0Abs, "Pi0", "XeXe_5.44TeV", 20, -1, kFALSE, centArray[cent]);
        for (Int_t i = 0; i< maxNBinsPi0[cent]; i++){
            cout << i << ": "<< xPtLimitsPi0[cent][i] <<" - " << xPtLimitsPi0[cent][i+1]<< ", " <<endl;
        }
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
    }

    // *******************************************************************************************************
    // ********************* Set file names for inputs *******************************************************
    // *******************************************************************************************************
    if (fileNamePCM.CompareTo("")!=0)       fileNamesMethod[0]    = fileNamePCM;
    if (fileNamePHOS.CompareTo("")!=0)      fileNamesMethod[1]    = fileNamePHOS;
    if (fileNameEMCAL.CompareTo("")!=0)     fileNamesMethod[2]    = fileNameEMCAL;
    if (fileNamePCMPHOS.CompareTo("")!=0)   fileNamesMethod[3]    = fileNamePCMPHOS;
    if (fileNamePCMEMCAL.CompareTo("")!=0)  fileNamesMethod[4]    = fileNamePCMEMCAL;

    if (fileConfigRXeXeErr.CompareTo("") != 0){
        ifstream inPbConfig(fileConfigRXeXeErr.Data());
        Int_t nReadSys  = 0;
        while(!inPbConfig.eof() && nReadSys<11 ){
            cout << "read line:" <<  nReadSys << "\t"<< nameMeasGlobalLabel[nReadSys].Data() << endl;
            TString nameCurrentSys  = "";
            inPbConfig >> nameCurrentSys ;
            if (nameCurrentSys.CompareTo(nameMeasGlobalLabel[nReadSys].Data()) != 0){
                cout << "wrong order in configuration file, was expecting " <<  nameMeasGlobalLabel[nReadSys].Data() << endl;
                return;
            } else {
                inPbConfig >> fileNamesXeXePi0DetailedSys[nReadSys] >> fileNamesppPi0DetailedSys[nReadSys] ;
                if (fileNamesXeXePi0DetailedSys[nReadSys].CompareTo("bla") != 0){
                    cout << "pi0 XeXe sys detailed: " << fileNamesXeXePi0DetailedSys[nReadSys] << endl;
                    havePi0SysDetailedXeXe[nReadSys]         = kTRUE;
                } else {
                    fileNamesXeXePi0DetailedSys[nReadSys]    = "";
                }
                if (fileNamesppPi0DetailedSys[nReadSys].CompareTo("bla") != 0){
                    cout << "pi0 pp sys detailed: " << fileNamesppPi0DetailedSys[nReadSys] << endl;
                    havePi0SysDetailedpp[nReadSys]          = kTRUE;
                } else {
                    fileNamesppPi0DetailedSys[nReadSys]     = "";
                }
                if (!(havePi0SysDetailedXeXe[nReadSys] && havePi0SysDetailedpp[nReadSys])){
                    havePi0SysDetailedXeXe[nReadSys]        = kTRUE;
                    havePi0SysDetailedpp[nReadSys]          = kTRUE;
                    fileNamesXeXePi0DetailedSys[nReadSys]   = "";
                    fileNamesppPi0DetailedSys[nReadSys]     = "";
                } else {
                    fileNamesRXeXePi0DetailedSys[nReadSys]  = Form("%s/Pi0RXeXe_%s_detailedSys.dat", outputDir.Data(),nameMeasGlobalLabel[nReadSys].Data());
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
    for (Int_t cent = 0; cent < 5; cent ++){
        cout << "Reading in cent " << centArray[cent].Data() << endl;
        for (Int_t meth = 0; meth < 11; meth++){
            if (fileNamesMethod[meth].CompareTo("") != 0){
                cout << "Reading in " << nameMeasGlobalLabel[meth].Data() << " from " << fileNamesMethod[meth].Data() << endl;
                if (!fileMethod[meth])
                    fileMethod[meth]                               = new TFile(fileNamesMethod[meth].Data());
                if (!fileMethod[meth]) {
                    cout << "file " << fileNamesMethod[meth].Data() << " not found! Skipping " << nameMeasGlobalLabel[meth].Data()  << endl;
                    continue;
                }
                directoryPi0[cent][meth]                             = (TDirectory*)fileMethod[meth]->Get(Form("Pi0%sXeXe_5.44TeV",centArray[cent].Data()));
                if (!directoryPi0[cent][meth]) {
                    cout << "File doesn't contain directory " << Form("Pi0%sXeXe_5.44TeV",centArray[cent].Data()) << "! Skipping " << nameMeasGlobalLabel[meth].Data()  << endl;
                    continue;
                }
                nTotMeas[cent]++;
                // reading supporting figures
                    graphPi0Mass[cent][meth]                             = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("Pi0_Mass_data");
                    graphPi0Mass[cent][meth]                             = ScaleGraph(graphPi0Mass[cent][meth], 1000.);
                    graphPi0Width[cent][meth]                            = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("Pi0_Width_data");
                    graphPi0Width[cent][meth]                            = ScaleGraph(graphPi0Width[cent][meth], 1000.);
                    graphPi0MassMC[cent][meth]                           = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("Pi0_Mass_MC");
                    graphPi0MassMC[cent][meth]                           = ScaleGraph(graphPi0MassMC[cent][meth], 1000.);
                    graphPi0WidthMC[cent][meth]                          = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("Pi0_Width_MC");
                    graphPi0WidthMC[cent][meth]                          = ScaleGraph(graphPi0WidthMC[cent][meth], 1000.);
                    graphPi0Acc[cent][meth]                              = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("AcceptancePi0");
                    graphPi0Eff[cent][meth]                              = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("EfficiencyPi0");
                    graphPi0EffTimesAcc[cent][meth]                      = (TGraphAsymmErrors*)directoryPi0[cent][meth]->Get("EffTimesAccPi0");
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
                    graphPi0InvYieldSys[cent][meth]->Print();

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
    Bool_t haveRefPPPi0[11]                                 = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                                                kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                                                kFALSE };
    for (Int_t meth = 0; meth< 11; meth++){
        statErrorCollectionPi0PP[meth]         = NULL;
        systErrorCollectionPi0PP[meth]         = NULL;
        systErrorUnCorrCollectionPi0PP[meth]   = NULL;
        systErrorInterCollectionPi0PP[meth]    = NULL;
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
        }

        TGraphAsymmErrors* graphPPCombPi0Stat           = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionStatErrComb_Pi0_5.023TeV");
        TGraphAsymmErrors* graphPPCombPi0UncorrSys      = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionUnCorrSystErrComb_Pi0_5.023TeV");
        TGraphAsymmErrors* graphPPCombPi0FullSys        = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionSystErrComb_Pi0_5.023TeV");
        TGraphAsymmErrors* graphPPCombPi0InterSys       = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionInterpolSystErrComb_Pi0_5.023TeV");

        TGraphAsymmErrors* graphPPInvYieldCombPi0Stat   = (TGraphAsymmErrors*)graphPPCombPi0Stat->Clone("graphInvYieldStatErrComb_Pi0_5.023TeV");
        TGraphAsymmErrors* graphPPInvYieldCombPi0Sys    = (TGraphAsymmErrors*)graphPPCombPi0FullSys->Clone("graphInvYieldSystErrComb_Pi0_5.023TeV");
        graphPPInvYieldCombPi0Stat                      = ScaleGraph(graphPPCombPi0Stat,1/(xSection5TeV*recalcBarn));
        graphPPInvYieldCombPi0Sys                       = ScaleGraph(graphPPInvYieldCombPi0Sys,1/(xSection5TeV*recalcBarn));



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
    histo2DPi0RatioToCombFitPP               = new TH2F("histo2DPi0RatioToCombFitPP","histo2DPi0RatioToCombFitPP",1000,0.23, 25.,1000,0.2,4.    );
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

    DrawGammaLines(0.23, 25. , 1., 1.,0.1, kGray+2);
    DrawGammaLines(0.23, 25. , 1.1, 1.1,0.1, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.9, 0.9,0.1, kGray, 7);

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

    DrawGammaLines(0.23, 25. , 1., 1.,0.5, kGray+2);
    DrawGammaLines(0.23, 25. , 1.1, 1.1,0.5, kGray, 7);
    DrawGammaLines(0.23, 25. , 0.9, 0.9,0.5, kGray, 7);
    DrawGammaLines(0.23, 25. , 1.2, 1.2,0.5, kGray, 9);
    DrawGammaLines(0.23, 25. , 0.8, 0.8,0.5, kGray, 9);

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

    TMarker* markerPCMPi0OnlyRatioPi0PP           = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[0],columnsLegendOnlyPi0RatioPP[1] ,rowsLegendOnlyPi0RatioPP[1],1);
    markerPCMPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[1] ,rowsLegendOnlyPi0RatioPPAbs[1]);
    TMarker* markerPHOSPi0OnlyRatioPi0PP          = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[1], columnsLegendOnlyPi0RatioPP[1] ,rowsLegendOnlyPi0RatioPP[2],1);
    markerPHOSPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[1] ,rowsLegendOnlyPi0RatioPPAbs[2]);
    TMarker* markerEMCALPi0OnlyRatioPi0PP         = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[2], columnsLegendOnlyPi0RatioPP[1] ,rowsLegendOnlyPi0RatioPP[3],1);
    markerEMCALPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[1] ,rowsLegendOnlyPi0RatioPPAbs[3]);
    TMarker* markerPCMEMCALPi0OnlyRatioPi0PP      = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[4], columnsLegendOnlyPi0RatioPP[3] ,rowsLegendOnlyPi0RatioPP[1],1);
    markerPCMEMCALPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[4] ,rowsLegendOnlyPi0RatioPPAbs[1]);
    TMarker* markerPCMPHOSPi0OnlyRatioPi0PP      = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[3], columnsLegendOnlyPi0RatioPP[3] ,rowsLegendOnlyPi0RatioPP[2],1);
    markerPCMPHOSPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[4] ,rowsLegendOnlyPi0RatioPPAbs[2]);

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
                                                                     columnsLegendOnlyPi0RatioPPAbs[5]+ 18*lengthBoxPP, rowsLegendOnlyPi0RatioPPAbs[1]+ heightBoxPP);
    boxPCMEMCALPi0OnlyRatioPi0PP->Draw("l");
    TBox* boxPCMPHOSPi0OnlyRatioPi0PP             = CreateBoxFromGraph(graphRatioPi0IndCombFitSysPP[3], columnsLegendOnlyPi0RatioPPAbs[5]-0.5*lengthBoxPP , rowsLegendOnlyPi0RatioPPAbs[2]- heightBoxPP,
                                                                     columnsLegendOnlyPi0RatioPPAbs[5]+ 18*lengthBoxPP, rowsLegendOnlyPi0RatioPPAbs[2]+ heightBoxPP);
    boxPCMPHOSPi0OnlyRatioPi0PP->Draw("l");

    canvasRatioToCombFitPP->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFit_PP5TeV.%s",outputDir.Data(),suffix.Data()));




    // *******************************************************************************************************
    // ************************** Combination of different pi0 measurements **********************************
    // *******************************************************************************************************
    // REMARKS:
    //     - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - extendable to other energies
    //     - offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented


    TH1D* statErrorCollectionPi0[5][11];
    TGraphAsymmErrors* statErrorGraphCollectionPi0[5][11];
    TGraphAsymmErrors* sysErrorCollectionPi0[5][11];
    for (Int_t cent = 0; cent < 5; cent++){
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
    Int_t offSetsPi0[5][11]             = { { 0,  6,  0,  0,  0,  3,  0,  0,  0,  0,  0 },
                                            { 0,  6,  0,  0,  0,  3,  0,  0,  0,  0,  0 },
                                            { 0,  6,  0,  0,  0,  3,  0,  0,  0,  0,  0 },
                                            { 0,  6,  0,  0,  0,  3,  0,  0,  0,  0,  0 },
                                            { 0,  6,  0,  0,  0,  3,  0,  0,  0,  0,  0 } };
    Int_t offSetsPi0Sys[5][11]          = { { 1,  7,  9,  3,  6,  4,  0,  0,  0,  21, 0 },
                                            { 1,  7,  9,  3,  6,  4,  0,  0,  0,  21, 0 },
                                            { 1,  7,  9,  3,  6,  4,  0,  0,  0,  21, 0 },
                                            { 1,  7,  9,  3,  6,  4,  0,  0,  0,  21, 0 },
                                            { 1,  7,  9,  3,  6,  4,  0,  0,  0,  21, 0 } };
    Int_t offSetPi0Shifting[5][11]      = { { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 } };
    Int_t nComBinsPi0Shifting[5][11]    = { { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 },
                                            { 30, 31, 30, 0,  31, 17,  0,  0,  0,  0, 0 } };

    TGraphAsymmErrors* statErrorRelCollectionPi0[5][11];
    TGraphAsymmErrors* sysErrorRelCollectionPi0[5][11];
    TGraph* graphWeightsPi0[5][11];
    for (Int_t cent = 0; cent < 5; cent++){
        for (Int_t meth = 0; meth< 11; meth++){
            graphWeightsPi0[cent][meth]                 = NULL;
            statErrorRelCollectionPi0[cent][meth]       = NULL;
            sysErrorRelCollectionPi0[cent][meth]        = NULL;
        }
    }
    for (Int_t cent = 0; cent < 5; cent++){
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

    TGraphAsymmErrors* graphCombPi0InvYieldStat[5]  = { NULL, NULL, NULL, NULL, NULL };
    TGraphAsymmErrors* graphCombPi0InvYieldSys[5]   = { NULL, NULL, NULL, NULL, NULL };
    TGraphAsymmErrors* graphCombPi0InvYieldTot[5]   = { NULL, NULL, NULL, NULL, NULL };
    TGraphAsymmErrors* graphIndPi0InvYieldStatUnshi[5][11];
    TGraphAsymmErrors* graphIndPi0InvYieldSysUnshi[5][11];
    TGraphAsymmErrors* graphIndPi0InvYieldStat[5][11];
    TGraphAsymmErrors* graphIndPi0InvYieldSys[5][11];
    TGraphAsymmErrors* graphIndPi0InvYieldStat_yShifted[5][11];
    TGraphAsymmErrors* graphIndPi0InvYieldSys_yShifted[5][11];
    for (Int_t cent = 0; cent < 5; cent++){
        for (Int_t meth = 0; meth< 11; meth++){
            graphIndPi0InvYieldStatUnshi[cent][meth]                = NULL;
            graphIndPi0InvYieldSysUnshi[cent][meth]                 = NULL;
            graphIndPi0InvYieldStat[cent][meth]                     = NULL;
            graphIndPi0InvYieldSys[cent][meth]                      = NULL;
            graphIndPi0InvYieldStat_yShifted[cent][meth]            = NULL;
            graphIndPi0InvYieldSys_yShifted[cent][meth]             = NULL;
        }
    }
    for (Int_t cent = 0; cent < 5; cent++){
    //     // Declaration & calculation of combined spectrum
    //     TString fileNamePi0OutputWeighting           = Form("%s/Pi0_WeightingMethod.dat",outputDir.Data());
    //     TGraphAsymmErrors* graphCombPi0InvYieldStat  = NULL;
    //     TGraphAsymmErrors* graphCombPi0InvYieldSys   = NULL;
    //     TGraphAsymmErrors* graphCombPi0InvYieldTot   = CombinePtPointsSpectraFullCorrMat(   statErrorCollectionPi0,    sysErrorCollectionPi0,
    //                                                                                         xPtLimitsPi0, maxNBinsPi0,
    //                                                                                         offSetsPi0, offSetsPi0Sys,
    //                                                                                         graphCombPi0InvYieldStat, graphCombPi0InvYieldSys,
    //                                                                                         fileNamePi0OutputWeighting, "XeXe_5.44TeV", "Pi0", kTRUE,
    //                                                                                         NULL, fileNameCorrFactors
    //                                                                                     );
    //
    //
    //     if (graphCombPi0InvYieldTot == NULL) {
    //         cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
    //         return;
    //     }
    //     while (graphCombPi0InvYieldStat->GetX()[0] < 0.3){
    //         graphCombPi0InvYieldStat->RemovePoint(0);
    //     }
    //     while (graphCombPi0InvYieldTot->GetX()[0] < 0.3){
    //         graphCombPi0InvYieldTot->RemovePoint(0);
    //     }
    //     while (graphCombPi0InvYieldSys->GetX()[0] < 0.3){
    //         graphCombPi0InvYieldSys->RemovePoint(0);
    //     }
    //     graphCombPi0InvYieldTot->Print();
    //
    //     // Reading weights from output file for plotting
    //     ifstream fileWeightsPi0Read;
    //     fileWeightsPi0Read.open(fileNamePi0OutputWeighting,ios_base::in);
    //     cout << "reading" << fileNamePi0OutputWeighting << endl;
    //     Double_t xValuesPi0Read[50];
    //     Double_t weightsPi0Read[11][50];
    //     Int_t availablePi0Meas[11]    = {   -1, -1, -1, -1, -1,
    //                                             -1, -1, -1, -1, -1,
    //                                             -1};
    //     Int_t nMeasSetPi0             = 4;
    //     Int_t nPtBinsPi0Read          = 0;
    //     while(!fileWeightsPi0Read.eof() && nPtBinsPi0Read < 50){
    //         TString garbage             = "";
    //         if (nPtBinsPi0Read == 0){
    //             fileWeightsPi0Read >> garbage ;//>> availablePi0Meas[0] >> availablePi0Meas[1] >> availablePi0Meas[2] >> availablePi0Meas[3];
    //             for (Int_t i = 0; i < nMeasSetPi0; i++){
    //                 fileWeightsPi0Read >> availablePi0Meas[i] ;
    //             }
    //             cout << "read following measurements: ";
    //             for (Int_t meth = 0; meth < 11; meth++){
    //                 cout << availablePi0Meas[i] << "\t" ;
    //             }
    //             cout << endl;
    //         } else {
    //             fileWeightsPi0Read >> xValuesPi0Read[nPtBinsPi0Read-1];
    //             for (Int_t i = 0; i < nMeasSetPi0; i++){
    //                 fileWeightsPi0Read >> weightsPi0Read[availablePi0Meas[i]][nPtBinsPi0Read-1] ;
    //             }
    //             cout << "read: "<<  nPtBinsPi0Read << "\t"<< xValuesPi0Read[nPtBinsPi0Read-1] << "\t" ;
    //             for (Int_t i = 0; i < nMeasSetPi0; i++){
    //                 cout << weightsPi0Read[availablePi0Meas[i]][nPtBinsPi0Read-1] << "\t";
    //             }
    //             cout << endl;
    //         }
    //         nPtBinsPi0Read++;
    //     }
    //     nPtBinsPi0Read                  = nPtBinsPi0Read-2 ;
    //     fileWeightsPi0Read.close();
    //
    //     for (Int_t i = 0; i < nMeasSetPi0; i++){
    //         graphWeightsPi0[availablePi0Meas[i]]                        = new TGraph(nPtBinsPi0Read,xValuesPi0Read,weightsPi0Read[availablePi0Meas[i]]);
    //         Int_t bin = 0;
    //         for (Int_t n = 0; n< nPtBinsPi0Read; n++){
    //             if (graphWeightsPi0[availablePi0Meas[i]]->GetY()[bin] == 0) graphWeightsPi0[availablePi0Meas[i]]->RemovePoint(bin);
    //             else bin++;
    //         }
    //     }
    //
    //
    //     // **********************************************************************************************************************
    //     // ******************************************* Plotting weights method only EMC *****************************************
    //     // **********************************************************************************************************************
    //     textSizeLabelsPixel                 = 900*0.04;
    //
    //     TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
    //     DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
    //     canvasWeights->SetLogx();
    //
    //     TH2F * histo2DPi0Weights;
    //     histo2DPi0Weights = new TH2F("histo2DPi0Weights","histo2DPi0Weights",11000,0.23, 25.,1000,-0.7,1.3);
    //     SetStyleHistoTH2ForGraphs(histo2DPi0Weights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    //     histo2DPi0Weights->GetXaxis()->SetMoreLogLabels();
    //     histo2DPi0Weights->GetXaxis()->SetLabelOffset(-0.01);
    //     histo2DPi0Weights->Draw("copy");
    //
    //         TLegend* legendWeights   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetPi0), 32);
    //         for (Int_t i = 0; i < nMeasSetPi0; i++){
    //             DrawGammaSetMarkerTGraph(graphWeightsPi0[availablePi0Meas[i]], markerStyleDet[availablePi0Meas[i]], markerSizeDet[availablePi0Meas[i]]*0.5, colorDet[availablePi0Meas[i]] , colorDet[availablePi0Meas[i]]);
    //             graphWeightsPi0[availablePi0Meas[i]]->Draw("p,same,z");
    //             legendWeights->AddEntry(graphWeightsPi0[availablePi0Meas[i]],nameMeasGlobalLabel[availablePi0Meas[i]],"p");
    //         }
    //         legendWeights->Draw();
    //
    //         TLatex *labelWeightsEnergy      = new TLatex(0.95,0.20,collisionSystemXeXe.Data());
    //         SetStyleTLatex( labelWeightsEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    //         labelWeightsEnergy->Draw();
    //         TLatex *labelWeightsPi0         = new TLatex(0.95,0.15,"#pi^{0} #rightarrow #gamma#gamma");
    //         SetStyleTLatex( labelWeightsPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    //         labelWeightsPi0->Draw();
    //
    // //      DrawGammaLines(0.23, 25. , 0.8, 0.8,0.1, kGray, 3);
    //         DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
    //         DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
    //         DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
    //         DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);
    //
    //     canvasWeights->SaveAs(Form("%s/Pi0_Weights.%s",outputDir.Data(),suffix.Data()));
    //
    //     //  *********************************************************************************************************************
    //     //  ************************************ Visualize relative errors ******************************************************
    //     //  *********************************************************************************************************************
    //
    //     TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
    //     DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    //     canvasRelSysErr->SetLogx();
    //
    //     TH2F * histo2DRelSysErr;
    //     histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.23, 25.,1000,0,50.0);
    //     SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    //     histo2DRelSysErr->GetYaxis()->SetRangeUser(0,39.5);
    //     histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
    //     histo2DRelSysErr->GetXaxis()->SetLabelOffset(-0.01);
    //     histo2DRelSysErr->Draw("copy");
    //
    //         TLegend* legendRelSysErr        = GetAndSetLegend2(0.62, 0.92-(0.04*nMeasSetPi0), 0.95, 0.92, textSizeLabelsPixel);
    //         for (Int_t i = 0; i < nMeasSetPi0; i++){
    //             DrawGammaSetMarkerTGraph(sysErrorRelCollectionPi0[availablePi0Meas[i]], markerStyleDet[availablePi0Meas[i]], markerSizeDet[availablePi0Meas[i]]*0.5, colorDet[availablePi0Meas[i]],
    //                                      colorDet[availablePi0Meas[i]]);
    //             sysErrorRelCollectionPi0[availablePi0Meas[i]]->Draw("p,same,z");
    //             legendRelSysErr->AddEntry(sysErrorRelCollectionPi0[availablePi0Meas[i]],nameMeasGlobalLabel[availablePi0Meas[i]],"p");
    //         }
    //         legendRelSysErr->Draw();
    //
    //         TLatex *labelRelSysErrEnergy    = new TLatex(0.15,0.89,collisionSystemXeXe.Data());
    //         SetStyleTLatex( labelRelSysErrEnergy, textSizeLabelsPixel, 4, 1, 43);
    //         labelRelSysErrEnergy->Draw();
    //         TLatex *labelRelSysErrPi0       = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
    //         SetStyleTLatex( labelRelSysErrPi0, textSizeLabelsPixel, 4, 1, 43);
    //         labelRelSysErrPi0->Draw();
    //
    //     canvasRelSysErr->SaveAs(Form("%s/Pi0_RelSysErr.%s",outputDir.Data(),suffix.Data()));
    //
    //     //  *********************************************************************************************************************
    //     //  ************************************ Visualize relative errors ******************************************************
    //     //  *********************************************************************************************************************
    //
    //     TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
    //     DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
    //     canvasRelStatErr->SetLogx();
    //
    //     TH2F * histo2DRelStatErr;
    //     histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.23, 25.,1000,0,50.5);
    //     SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    //     histo2DRelStatErr->GetYaxis()->SetRangeUser(0,39.5);
    //     histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
    //     histo2DRelStatErr->GetXaxis()->SetLabelOffset(-0.01);
    //     histo2DRelStatErr->Draw("copy");
    //         TLegend* legendRelStatErr       = GetAndSetLegend2(0.14, 0.92-(0.04*nMeasSetPi0), 0.45, 0.92, textSizeLabelsPixel);
    //         for (Int_t i = 0; i < nMeasSetPi0; i++){
    //             DrawGammaSetMarkerTGraph(statErrorRelCollectionPi0[availablePi0Meas[i]], markerStyleDet[availablePi0Meas[i]], markerSizeDet[availablePi0Meas[i]]*0.5, colorDet[availablePi0Meas[i]],
    //                                      colorDet[availablePi0Meas[i]]);
    //             statErrorRelCollectionPi0[availablePi0Meas[i]]->Draw("p,same,z");
    //             legendRelStatErr->AddEntry(statErrorRelCollectionPi0[availablePi0Meas[i]],nameMeasGlobalLabel[availablePi0Meas[i]],"p");
    //         }
    //         legendRelStatErr->Draw();
    //
    //         TLatex *labelRelStatErrEnergy   = new TLatex(0.95,0.89,collisionSystemXeXe.Data());
    //         SetStyleTLatex( labelRelStatErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    //         labelRelStatErrEnergy->Draw();
    //         TLatex *labelRelStatErrPi0      = new TLatex(0.95,0.85,"#pi^{0} #rightarrow #gamma#gamma");
    //         SetStyleTLatex( labelRelStatErrPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    //         labelRelStatErrPi0->Draw();
    //
    //     canvasRelStatErr->SaveAs(Form("%s/Pi0_RelStatErr.%s",outputDir.Data(),suffix.Data()));
    //
    //     //  *********************************************************************************************************************
    //     //  ************************ Visualize relative total errors of different combination methods Pi0 ***********************
    //     //  *********************************************************************************************************************
    //     TGraphAsymmErrors* graphCombPi0InvYieldRelStat  = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldStat, "relativeStatErrorPi0_Method");
    //     TGraphAsymmErrors* graphCombPi0InvYieldRelSys   = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldSys, "relativeSysErrorPi0_Method");
    //     TGraphAsymmErrors* graphCombPi0InvYieldRelTot   = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldTot, "relativeTotalErrorPi0_Method");
    //
    //     TCanvas* canvasRelTotErr            = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
    //     DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
    //     canvasRelTotErr->SetLogx();
    //
    //     TH2F * histo2DRelTotErrPi0;
    //     histo2DRelTotErrPi0                 = new TH2F("histo2DRelTotErrPi0","histo2DRelTotErrPi0",11000,0.23, 25.,1000,0,50.0);
    //     SetStyleHistoTH2ForGraphs(histo2DRelTotErrPi0, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    //     histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,39.5);
    //     histo2DRelTotErrPi0->GetXaxis()->SetMoreLogLabels();
    //     histo2DRelTotErrPi0->GetXaxis()->SetLabelOffset(-0.01);
    //     histo2DRelTotErrPi0->Draw("copy");
    //
    //         DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot, markerStyleComb, markerSizeComb, kBlue+2 , kBlue+2);
    //         graphCombPi0InvYieldRelTot->Draw("p,same,z");
    //
    //         TLegend* legendRelTotErr1       = GetAndSetLegend2(0.20, 0.92-(0.035), 0.45, 0.92, 32);
    //         legendRelTotErr1->AddEntry(graphCombPi0InvYieldRelTot,"All","p");
    //         legendRelTotErr1->Draw();
    //
    //         TLatex *labelRelTotErrEnergy    = new TLatex(0.95,0.89,collisionSystemXeXe.Data());
    //         SetStyleTLatex( labelRelTotErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    //         labelRelTotErrEnergy->Draw();
    //         TLatex *labelRelTotErrPi0       = new TLatex(0.95,0.85,"#pi^{0} #rightarrow #gamma#gamma");
    //         SetStyleTLatex( labelRelTotErrPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    //         labelRelTotErrPi0->Draw();
    //
    //     canvasRelTotErr->SaveAs(Form("%s/Pi0_TotErr_Comp.%s",outputDir.Data(),suffix.Data()));
    //     histo2DRelTotErrPi0->GetYaxis()->SetTitle("Err (%)");
    //     histo2DRelTotErrPi0->Draw("copy");
    //
    //         DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
    //         graphCombPi0InvYieldRelTot->Draw("p,same,z");
    //         DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
    //         graphCombPi0InvYieldRelStat->Draw("l,x0,same,e1");
    //         DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
    //         graphCombPi0InvYieldRelSys->SetLineStyle(7);
    //         graphCombPi0InvYieldRelSys->Draw("l,x0,same,e1");
    //
    //         TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
    //         legendRelTotErr3->AddEntry(graphCombPi0InvYieldRelTot,"tot","p");
    //         legendRelTotErr3->AddEntry(graphCombPi0InvYieldRelStat,"stat","l");
    //         legendRelTotErr3->AddEntry(graphCombPi0InvYieldRelSys,"sys","l");
    //         legendRelTotErr3->Draw();
    //
    //         labelRelTotErrEnergy->Draw();
    //         labelRelTotErrPi0->Draw();
    //
    //     canvasRelTotErr->SaveAs(Form("%s/Pi0_Reldecomp.%s",outputDir.Data(),suffix.Data()));
    //
    //
    //     // **********************************************************************************************************************
    //     // ************************************* Calculating bin shifted spectra & fitting **************************************
    //     // **********************************************************************************************************************
    //
    //     // Cloning spectra
    //     TGraphAsymmErrors* graphCombPi0InvYieldTotUnshi         = (TGraphAsymmErrors*)graphCombPi0InvYieldTot->Clone("Pi0Unshifted");
    //     TGraphAsymmErrors* graphCombPi0InvYieldStatUnshi        = (TGraphAsymmErrors*)graphCombPi0InvYieldStat->Clone("Pi0UnshiftedStat");
    //     TGraphAsymmErrors* graphCombPi0InvYieldSysUnshi         = (TGraphAsymmErrors*)graphCombPi0InvYieldSys->Clone("Pi0UnshiftedSys");
    //
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
    //
    //     // fitting spectrum with intial parameters
    //     // Two component model fit from Bylinkin
    //     TF1* fitTCMDecomposedLPi0           = FitObject("tcmlow","twoCompModelPi0_DecL", "Pi0", NULL, 0.3, 2.);
    //     TF1* fitTCMDecomposedHPi0           = FitObject("tcmhigh","twoCompModelPi0_DecH", "Pi0", NULL, 4, 50.);
    //     fitTCMDecomposedLPi0->SetParameters(graphCombPi0InvYieldTot->GetY()[2],0.3);
    //     graphCombPi0InvYieldStat->Fit(fitTCMDecomposedLPi0,"QNRMEX0+","",0.3,0.8);
    //     graphCombPi0InvYieldStat->Fit(fitTCMDecomposedHPi0,"QNRMEX0+","",3,20);
    // //     fitTCMDecomposedHPi0->SetParameters(graphCombPi0InvYieldTot->GetY()[2],0.8, 2);
    //
    //     cout << WriteParameterToFile(fitTCMDecomposedLPi0)<< endl;
    //     fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedLPi0)<< endl;
    //     cout << WriteParameterToFile(fitTCMDecomposedHPi0)<< endl;
    //     fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedHPi0)<< endl;
    //
    //     Double_t paramTCMPi0New[5]          = { 20.72,0.167, 0.454, 0.629, 3.163};
    //
    // //     Double_t paramTCMPi0New[5]          = { 5.265,0.33,
    // //                                             1.9,0.46,3.1};
    //     TF1* fitTCMInvYieldPi0              = FitObject("tcm","fitTCMInvYieldPi0","Pi0",graphCombPi0InvYieldStat,0.3,20. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
    //
    //     fitTCMDecomposedLPi0->SetParameter(0, fitTCMInvYieldPi0->GetParameter(0));
    //     fitTCMDecomposedLPi0->SetParameter(1, fitTCMInvYieldPi0->GetParameter(1));
    //     fitTCMDecomposedHPi0->SetParameter(0, fitTCMInvYieldPi0->GetParameter(2));
    //     fitTCMDecomposedHPi0->SetParameter(1, fitTCMInvYieldPi0->GetParameter(3));
    //     fitTCMDecomposedHPi0->SetParameter(2, fitTCMInvYieldPi0->GetParameter(4));
    //
    //     // Tsallis fit
    //     Double_t paramGraphPi0[3]           = {1.0e12, 8., 0.13};
    //     TF1* fitInvYieldPi0                 = FitObject("l","fitInvYieldPi0","Pi0",histoPi0InvYieldStat[2],0.3,20.,paramGraphPi0,"QNRMEX0+");
    //     TF1* fitInvYieldPi0Graph            = (TF1*)fitInvYieldPi0->Clone("fitInvYieldPi0Graph");
    //
    //     // *************************************************************************************************************
    //     // Shift graphs in X direction if desired
    //     // *************************************************************************************************************
    //     if(bWCorrection.Contains("X")){
    //         TF1* fitShiftingPi0            = FitObject("tmpt","ShiftingPi0","Pi0");
    //         fitShiftingPi0->SetParameters(fitInvYieldPi0->GetParameter(0),fitInvYieldPi0->GetParameter(1), fitInvYieldPi0->GetParameter(2));
    // //         TF1* fitShiftingPi0                 = FitObject("tcmpt","ShiftingPi0","Pi0");
    // //         fitShiftingPi0->SetParameters(fitTCMInvYieldPi0->GetParameter(0),fitTCMInvYieldPi0->GetParameter(1), fitTCMInvYieldPi0->GetParameter(2), fitTCMInvYieldPi0->GetParameter(3),fitTCMInvYieldPi0->GetParameter(4));
    //
    //         TGraphAsymmErrors* graphCombPi0InvYieldTotNoShift = (TGraphAsymmErrors*) graphCombPi0InvYieldTot->Clone("Pi0_NoShift");
    //
    //         graphCombPi0InvYieldTot          = ApplyXshift(graphCombPi0InvYieldTot, fitShiftingPi0);
    //         cout << "comb" << endl;
    //         graphCombPi0InvYieldStat->Print();
    //         graphCombPi0InvYieldStat         = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot,
    //                                                                             graphCombPi0InvYieldStat,
    //                                                                             fitShiftingPi0,
    //                                                                             0, graphCombPi0InvYieldStat->GetN());
    //         graphCombPi0InvYieldSys          = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot,
    //                                                                             graphCombPi0InvYieldSys,
    //                                                                             fitShiftingPi0,
    //                                                                             0, graphCombPi0InvYieldSys->GetN());
    //         for (Int_t meth = 0; meth< 11; meth++){
    //             if (graphIndPi0InvYieldStat[meth]){
    //                 cout << "shiting stat err of " << nameMeasGlobalLabel[meth].Data();
    //                 graphIndPi0InvYieldStat[meth]  = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot,
    //                                                                             graphIndPi0InvYieldStat[meth],
    //                                                                             fitShiftingPi0,
    //                                                                             offSetPi0Shifting[meth], nComBinsPi0Shifting[meth]);
    //
    //             }
    //             if (graphIndPi0InvYieldSys[meth]){
    //                 cout << "shiting sys err of " << nameMeasGlobalLabel[meth].Data();
    //                 graphIndPi0InvYieldSys[meth]   = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot,
    //                                                                             graphIndPi0InvYieldSys[meth],
    //                                                                             fitShiftingPi0,
    //                                                                             offSetPi0Shifting[meth], nComBinsPi0Shifting[meth]);
    //             }
    //         }
    //
    //         //***************************************************************************************************************
    //         //************************************Plotting binshift corrections *********************************************
    //         //***************************************************************************************************************
    //
    //         TCanvas* canvasShift = new TCanvas("canvasShift","",0,0,1000,900);// gives the page size
    //         DrawGammaCanvasSettings( canvasShift, 0.12, 0.017, 0.015, 0.08);
    //
    //         Size_t textSizeSpectra          = 0.04;
    //         TH1F * histoBinShift = new TH1F("histoBinShift","histoBinShift",1000,0.,20);
    //         SetStyleHistoTH1ForGraphs(histoBinShift, "#it{p}_{T} (GeV/#it{c})","bin shifted (X) / no shift",
    //                                 0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, 0.85,1.2);
    //         histoBinShift->GetYaxis()->SetRangeUser(0.95,1.05);
    //         histoBinShift->DrawCopy();
    //
    //             Int_t numberPoints      = graphCombPi0InvYieldTotNoShift->GetN();
    //             Double_t *xPoint        = graphCombPi0InvYieldTotNoShift->GetX();
    //             Double_t* xvalueErrUp   = graphCombPi0InvYieldTotNoShift->GetEXhigh();
    //             Double_t* xvalueErrLow  = graphCombPi0InvYieldTotNoShift->GetEXlow();
    //             Double_t *xPointShift= graphCombPi0InvYieldTot->GetX();
    //             for (Int_t i=0; i<numberPoints; i++) {
    //                 graphCombPi0InvYieldTotNoShift->SetPoint(i,xPoint[i],xPointShift[i]/xPoint[i]);
    //                 graphCombPi0InvYieldTotNoShift->SetPointError(i,xvalueErrLow[i],xvalueErrUp[i],0,0);
    //             }
    //             DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldTotNoShift, markerStyleComb, markerSizeComb, colorComb , colorComb);
    //             graphCombPi0InvYieldTotNoShift->Draw("p same");
    //
    //             TLatex *labelRatioToFitBinShift   = new TLatex(0.94, 0.91, collisionSystemXeXe.Data());
    //             SetStyleTLatex( labelRatioToFitBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    //             labelRatioToFitBinShift->Draw();
    //             TLatex *labelRatioToFitALICEBinShift    = new TLatex(0.94, 0.86, textALICE.Data());
    //             SetStyleTLatex( labelRatioToFitALICEBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    //             labelRatioToFitALICEBinShift->Draw();
    //             TLatex *labelRatioToFitPi0BinShift      = new TLatex(0.94, 0.807, "#pi^{0} #rightarrow #gamma#gamma");
    //             SetStyleTLatex( labelRatioToFitPi0BinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    //             labelRatioToFitPi0BinShift->Draw();
    //
    //         canvasShift->Update();
    //         canvasShift->SaveAs(Form("%s/BinShiftCorrection_Pi0.%s",outputDir.Data(),suffix.Data()));
    //
    //         // *************************************************************************************************************
    //         // Plot control graphs
    //         // *************************************************************************************************************
    //
    //         TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
    //         DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.09);
    //         canvasDummy2->SetLogy();
    //         canvasDummy2->SetLogx();
    //         TH2F * histo2DDummy2;
    //         histo2DDummy2               = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.33,31.,1000,1e-9,10e1);
    //         SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
    //         histo2DDummy2->DrawCopy();
    //
    //             for (Int_t meth = 0; meth< 11; meth++){
    //                 if (graphIndPi0InvYieldSys[meth]){
    //                     DrawGammaSetMarkerTGraphAsym(graphIndPi0InvYieldSys[meth], markerStyleDet[meth] ,markerSizeDet[meth]/2, colorDet[meth], colorDet[meth], widthLinesBoxes, kTRUE);
    //                     graphIndPi0InvYieldSys[meth]->Draw("pEsame");
    //                 }
    //             }
    //
    //             DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStatUnshi, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
    //             graphCombPi0InvYieldStatUnshi->Draw("pEsame");
    //             DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
    //             graphCombPi0InvYieldStat->Draw("pEsame");
    //
    //             fitTCMInvYieldPi0->SetLineColor(kBlue+2);
    //             fitTCMInvYieldPi0->Draw("same");
    //
    //         canvasDummy2->Update();
    //         canvasDummy2->Print(Form("%s/ComparisonShiftedPi0_PP5TeV.%s",outputDir.Data(),suffix.Data()));
    //         delete canvasDummy2;
    //     }
    //
    //     // *************************************************************************************************************
    //     // redo fitting after binshifts
    //     // *************************************************************************************************************
    //     // Tsallis function
    //     cout << WriteParameterToFile(fitInvYieldPi0)<< endl;
    //     fileFitsOutput <<  WriteParameterToFile(fitInvYieldPi0)<< endl;
    //     //Two component model from Bylinkin
    //     fitTCMInvYieldPi0        = FitObject("tcm","fitTCMInvYieldPi0XeXe5TeV","Pi0",graphCombPi0InvYieldTot,0.3,20. ,paramTCMPi0New,"QNRMEX0+","", kFALSE);
    //     cout << WriteParameterToFile(fitTCMInvYieldPi0)<< endl;
    //     fileFitsOutput <<  WriteParameterToFile(fitTCMInvYieldPi0)<< endl;
    //
    //     TF1* fitPowInvYieldPi0   = FitObject("m","fitPowInvYieldPi0XeXe5TeV","Pi0",graphCombPi0InvYieldTot,5,40. ,NULL,"QNRMEX0+","", kFALSE);
    //     cout << WriteParameterToFile(fitPowInvYieldPi0)<< endl;
    //
    //     // *************************************************************************************************************
    //     // Shift graphs in Y direction as well if desired
    //     // *************************************************************************************************************
    //     TGraphAsymmErrors* graphCombPi0InvYieldTot_yShifted         = NULL;
    //     TGraphAsymmErrors* graphCombPi0InvYieldStat_yShifted        = NULL;
    //     TGraphAsymmErrors* graphCombPi0InvYieldSys_yShifted         = NULL;
    //
    //     if(bWCorrection.Contains("Y") ){
    //         graphCombPi0InvYieldTot_yShifted        = (TGraphAsymmErrors*)graphCombPi0InvYieldTotUnshi->Clone("Pi0YShiftedCombTot");
    //         graphCombPi0InvYieldTot_yShifted        =  ApplyYshiftIndividualSpectra( graphCombPi0InvYieldTot_yShifted, fitInvYieldPi0);
    //         graphCombPi0InvYieldStat_yShifted       = (TGraphAsymmErrors*)graphCombPi0InvYieldStatUnshi->Clone("Pi0YShiftedCombStat");
    //         graphCombPi0InvYieldStat_yShifted       =  ApplyYshiftIndividualSpectra( graphCombPi0InvYieldStat_yShifted, fitInvYieldPi0);
    //         graphCombPi0InvYieldSys_yShifted        = (TGraphAsymmErrors*)graphCombPi0InvYieldSysUnshi->Clone("Pi0YShiftedCombSys");
    //         graphCombPi0InvYieldSys_yShifted        =  ApplyYshiftIndividualSpectra( graphCombPi0InvYieldSys_yShifted, fitInvYieldPi0);
    //
    //         for (Int_t meth = 0; meth< 11; meth++){
    //             if (graphIndPi0InvYieldStat_yShifted[meth]){
    //                 graphIndPi0InvYieldStat_yShifted[meth] = ApplyYshiftIndividualSpectra( graphIndPi0InvYieldStat_yShifted[meth], fitInvYieldPi0);
    //             }
    //             if (graphIndPi0InvYieldSys_yShifted[meth]){
    //                 graphIndPi0InvYieldSys_yShifted[meth]  = ApplyYshiftIndividualSpectra( graphIndPi0InvYieldSys_yShifted[meth], fitInvYieldPi0);
    //             }
    //         }
    //     }
    //
    //     // *************************************************************************************************************
    //     // Calculate ratios to combined fit
    //     // *************************************************************************************************************
    // //     TH1D* histoRatioPi0DPMJetToFit                      = (TH1D*) histoDPMJetPi0->Clone("histoRatioPi0DPMJetToFit");
    // //     histoRatioPi0DPMJetToFit                            = CalculateHistoRatioToFit (histoRatioPi0DPMJetToFit, fitTCMInvYieldPi0);
    // //     histoRatioPi0DPMJetToFit->GetXaxis()->SetRangeUser(0.3, 20);
    // //     TH1D* histoRatioPi0HIJINGToFit                      = (TH1D*) histoHIJINGPi0->Clone("histoRatioPi0HIJINGToFit");
    // //     histoRatioPi0HIJINGToFit                            = CalculateHistoRatioToFit (histoRatioPi0HIJINGToFit, fitTCMInvYieldPi0);
    // //     histoRatioPi0HIJINGToFit->GetXaxis()->SetRangeUser(0.3, 20);
    // //     TH1D* histoRatioPi0EPOS3ToFit                       = (TH1D*) histoEPOS3Pi0->Clone("histoRatioPi0EPOS3ToFit");
    // //     histoRatioPi0EPOS3ToFit                             = CalculateHistoRatioToFit (histoRatioPi0EPOS3ToFit, fitTCMInvYieldPi0);
    // //     histoRatioPi0EPOS3ToFit->GetXaxis()->SetRangeUser(0.3, 20);
    // //     TGraphErrors* graphRatioPi0EPOS3ToFit               = new TGraphErrors(histoRatioPi0EPOS3ToFit);
    // //     while(graphRatioPi0EPOS3ToFit->GetX()[0] < 0.3)
    // //         graphRatioPi0EPOS3ToFit->RemovePoint(0);
    // //     while(graphRatioPi0EPOS3ToFit->GetX()[graphRatioPi0EPOS3ToFit->GetN()-1] > 20)
    // //         graphRatioPi0EPOS3ToFit->RemovePoint(graphRatioPi0EPOS3ToFit->GetN()-1);
    // //     TGraphErrors* graphRatioPi0McGillToFit              = (TGraphErrors*)graphMcGillPi0->Clone("graphRatioPi0McGillToFit");
    // //     graphRatioPi0McGillToFit                            = CalculateGraphErrRatioToFit(graphRatioPi0McGillToFit, fitTCMInvYieldPi0);
    // //     while(graphRatioPi0McGillToFit->GetX()[0] < 0.3)
    // //         graphRatioPi0McGillToFit->RemovePoint(0);
    //
    //
    // //     TGraph* graphRatioPi0CGC                            = (TGraph*)graphPi0CGC5TeV->Clone("graphRatioPi0CGC");
    // //     graphRatioPi0CGC                                    = CalculateGraphRatioToFit(graphRatioPi0CGC, fitTCMInvYieldPi0);
    // //     TGraphAsymmErrors* graphRatioPi0NLODSS14            = (TGraphAsymmErrors*)graphNLODSS14Pi0->Clone("graphRatioPi0NLODSS14ToFit");
    // //     graphRatioPi0NLODSS14                               = CalculateGraphErrRatioToFit(graphRatioPi0NLODSS14, fitTCMInvYieldPi0);
    // //     TGraph* graphRatioPi0NLODSS14Center                 = (TGraph*)graphNLODSS14Pi0Center->Clone("graphRatioPi0NLODSS14CenterToFit");
    // //     graphRatioPi0NLODSS14Center                         = CalculateGraphRatioToFit(graphRatioPi0NLODSS14Center, fitTCMInvYieldPi0);
    // //     TGraphAsymmErrors* graphRatioPi0NLODSS14nPDF        = (TGraphAsymmErrors*)graphNLODSS14nPDFPi0->Clone("graphRatioPi0NLODSS14nPDFToFit");
    // //     graphRatioPi0NLODSS14nPDF                           = CalculateGraphErrRatioToFit(graphRatioPi0NLODSS14nPDF, fitTCMInvYieldPi0);
    // //     TGraph* graphRatioPi0NLODSS14nPDFCenter             = (TGraph*)graphNLODSS14nPDFPi0Center->Clone("graphRatioPi0NLODSS14nPDFCenterToFit");
    // //     graphRatioPi0NLODSS14nPDFCenter                     = CalculateGraphRatioToFit(graphRatioPi0NLODSS14nPDFCenter, fitTCMInvYieldPi0);
    // //     TGraphAsymmErrors* graphRatioPi0NLODSS14nPDFEPPS16  = (TGraphAsymmErrors*)graphNLODSS14nPDFEPPS16Pi0->Clone("graphRatioPi0NLODSS14nPDFEPPS16ToFit");
    // //     graphRatioPi0NLODSS14nPDFEPPS16                     = CalculateGraphErrRatioToFit(graphRatioPi0NLODSS14nPDFEPPS16, fitTCMInvYieldPi0);
    // //     TGraph* graphRatioPi0NLODSS14nPDFEPPS16Center       = (TGraph*)graphNLODSS14nPDFEPPS16Pi0Center->Clone("graphRatioPi0NLODSS14nPDFEPPS16CenterToFit");
    // //     graphRatioPi0NLODSS14nPDFEPPS16Center               = CalculateGraphRatioToFit(graphRatioPi0NLODSS14nPDFEPPS16Center, fitTCMInvYieldPi0);
    //
    //
    //     TGraphAsymmErrors* graphRatioPi0IndCombFitStat[11]  = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    //     TGraphAsymmErrors* graphRatioPi0IndCombFitSys[11]   = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    //
    //     TGraphAsymmErrors* graphRatioPi0CombCombFitTot      = (TGraphAsymmErrors*)graphCombPi0InvYieldTot->Clone();
    //     graphRatioPi0CombCombFitTot                         = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitTot, fitTCMInvYieldPi0);
    //     TGraphAsymmErrors* graphRatioPi0CombCombFitStat     = (TGraphAsymmErrors*)graphCombPi0InvYieldStat->Clone();
    //     graphRatioPi0CombCombFitStat                        = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitStat, fitTCMInvYieldPi0);
    //     TGraphAsymmErrors* graphRatioPi0CombCombFitSys      = (TGraphAsymmErrors*)graphCombPi0InvYieldSys->Clone();
    //     graphRatioPi0CombCombFitSys                         = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitSys, fitTCMInvYieldPi0);
    //
    //     for (Int_t meth = 0; meth< 11; meth++){
    //         if (graphIndPi0InvYieldStat[meth]){
    //             graphRatioPi0IndCombFitStat[meth]              = (TGraphAsymmErrors*)graphIndPi0InvYieldStat[meth]->Clone(Form("RatioPi0%sToCombFitStat", nameMeasGlobalLabel[meth].Data()));
    //             graphRatioPi0IndCombFitStat[meth]              = CalculateGraphErrRatioToFit(graphRatioPi0IndCombFitStat[meth], fitTCMInvYieldPi0);
    //         }
    //         if (graphIndPi0InvYieldSys[meth]){
    //             graphRatioPi0IndCombFitSys[meth]              = (TGraphAsymmErrors*)graphIndPi0InvYieldSys[meth]->Clone(Form("RatioPi0%sToCombFitSyst", nameMeasGlobalLabel[meth].Data()));
    //             graphRatioPi0IndCombFitSys[meth]              = CalculateGraphErrRatioToFit(graphRatioPi0IndCombFitSys[meth], fitTCMInvYieldPi0);
    //         }
    //     }
    //
    //     // **********************************************************************************************************************
    //     // ******************************************* Plot Ratio of Comb to Fit ************************************************
    //     // **********************************************************************************************************************
    //     textSizeLabelsPixel                 = 54;
    //     TCanvas* canvasRatioToCombFit       = new TCanvas("canvasRatioToCombFit","",200,10,1350,900);  // gives the page size
    //     DrawGammaCanvasSettings( canvasRatioToCombFit, 0.08, 0.01, 0.01, 0.125);
    //     canvasRatioToCombFit->SetLogx();
    //
    //         Double_t textsizeLabelsXeXe      = 0;
    //         if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
    //             textsizeLabelsXeXe           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
    //         } else {
    //             textsizeLabelsXeXe           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
    //         }
    //         cout << textsizeLabelsXeXe << endl;
    //
    //     TH2F * histo2DPi0RatioToCombFit;
    //     histo2DPi0RatioToCombFit               = new TH2F("histo2DPi0RatioToCombFit","histo2DPi0RatioToCombFit",1000,0.23, 25.,1000,0.2,4.    );
    //     SetStyleHistoTH2ForGraphs(histo2DPi0RatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsXeXe, textsizeLabelsXeXe,
    //                               0.85*textsizeLabelsXeXe,textsizeLabelsXeXe, 0.9, 0.65, 510, 505);
    //     histo2DPi0RatioToCombFit->GetXaxis()->SetMoreLogLabels();
    //     histo2DPi0RatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
    // //  histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(-10,10);
    //     histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.2,1.82);
    //     histo2DPi0RatioToCombFit->Draw("copy");
    //
    //         ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStat);
    //
    //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
    //         graphRatioPi0CombCombFitSys->Draw("E2same");
    //         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStat, markerStyleComb, markerSizeComb, colorComb , colorComb);
    //         graphRatioPi0CombCombFitStat->Draw("p,same,z");
    //
    //         DrawGammaLines(0.23, 25. , 1., 1.,0.1, kGray+2);
    //         DrawGammaLines(0.23, 25. , 1.1, 1.1,0.1, kGray, 7);
    //         DrawGammaLines(0.23, 25. , 0.9, 0.9,0.1, kGray, 7);
    //
    //         TLatex *labelRatioToFitEnergy   = new TLatex(0.95, 0.92, collisionSystemXeXe.Data());
    //         SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    //         labelRatioToFitEnergy->Draw();
    //         TLatex *labelRatioToFitALICE    = new TLatex(0.95, 0.86, textALICE.Data());
    //         SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    //         labelRatioToFitALICE->Draw();
    //         TLatex *labelRatioToFitPi0      = new TLatex(0.12, 0.92, "#pi^{0} #rightarrow #gamma#gamma");
    //         SetStyleTLatex( labelRatioToFitPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
    //         labelRatioToFitPi0->Draw();
    //
    //     canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToCombFit_XeXe5TeV.%s",outputDir.Data(),suffix.Data()));
    //
    //     // **********************************************************************************************************************
    //     // *******************************************Plot Ratio of Individual meas to Fit ******************************************
    //     // **********************************************************************************************************************
    //
    //     canvasRatioToCombFit->cd();
    //     histo2DPi0RatioToCombFit->Draw("copy");
    //
    //         for (Int_t meth = 10; meth > -1 ; meth--){
    //             if (graphRatioPi0IndCombFitSys[meth]){
    //                 DrawGammaSetMarkerTGraphAsym(graphRatioPi0IndCombFitSys[meth], markerStyleDet[meth] ,markerSizeDet[meth]*0.5, colorDet[meth], colorDet[meth], widthLinesBoxes, kTRUE);
    //                 graphRatioPi0IndCombFitSys[meth]->Draw("E2same");
    //             }
    //             if (graphRatioPi0IndCombFitStat[meth]){
    //                 ProduceGraphAsymmWithoutXErrors(graphRatioPi0IndCombFitStat[meth]);
    //                 DrawGammaSetMarkerTGraphAsym(graphRatioPi0IndCombFitStat[meth], markerStyleDet[meth] ,markerSizeDet[meth]*0.5, colorDet[meth], colorDet[meth]);
    //                 graphRatioPi0IndCombFitStat[meth]->Draw("p,same,z");
    //             }
    //         }
    //         graphRatioPi0IndCombFitStat[4]->Draw("p,same,z");
    //
    //         DrawGammaLines(0.23, 25. , 1., 1.,0.5, kGray+2);
    //         DrawGammaLines(0.23, 25. , 1.1, 1.1,0.5, kGray, 7);
    //         DrawGammaLines(0.23, 25. , 0.9, 0.9,0.5, kGray, 7);
    //         DrawGammaLines(0.23, 25. , 1.2, 1.2,0.5, kGray, 9);
    //         DrawGammaLines(0.23, 25. , 0.8, 0.8,0.5, kGray, 9);
    //
    //         labelRatioToFitEnergy->Draw();
    //         labelRatioToFitALICE->Draw();
    //         labelRatioToFitPi0->Draw();
    //         histo2DPi0RatioToCombFit->Draw("same,axis");
    //
    //         //****************************** Definition of the Legend ******************************************
    //         //**************** Row def ************************
    //         Double_t rowsLegendOnlyPi0Ratio[4]          = {0.30, 0.25, 0.20, 0.15};
    //         Double_t rowsLegendOnlyPi0RatioAbs[4]       = {0.92, 0.47, 0.38, 0.29 };
    //         Double_t columnsLegendOnlyPi0Ratio[6]       = {0.14, 0.26, 0.35, 0.48, 0.7, 0.8};
    //         Double_t columnsLegendOnlyPi0RatioAbs[6]    = {0.215, 0.69, 0.95, 2, 6.6, 9.6};
    //         Double_t lengthBox                          = 0.2;
    //         Double_t heightBox                          = 0.06/2;
    //         //****************** first Column **************************************************
    //         TLatex *textPCMOnlyRatioPi0                 = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[1],nameMeasGlobalLabel[0]);
    //         SetStyleTLatex( textPCMOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
    //         textPCMOnlyRatioPi0->Draw();
    //         TLatex *textPHOSOnlyRatioPi0                = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[2],nameMeasGlobalLabel[1]);
    //         SetStyleTLatex( textPHOSOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
    //         textPHOSOnlyRatioPi0->Draw();
    //         TLatex *textEMCALOnlyRatioPi0               = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[3],nameMeasGlobalLabel[2]);
    //         SetStyleTLatex( textEMCALOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
    //         textEMCALOnlyRatioPi0->Draw();
    //         TLatex *textPCMEMCALOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[3],rowsLegendOnlyPi0Ratio[1],nameMeasGlobalLabel[4]);
    //         SetStyleTLatex( textPCMEMCALOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
    //         textPCMEMCALOnlyRatioPi0->Draw();
    //         TLatex *textPCMPHOSOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[3],rowsLegendOnlyPi0Ratio[2],nameMeasGlobalLabel[3]);
    //         SetStyleTLatex( textPCMPHOSOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
    //         textPCMPHOSOnlyRatioPi0->Draw();
    //
    //         //****************** second Column *************************************************
    //         TLatex *textStatOnlyRatioPi0                = new TLatex(columnsLegendOnlyPi0Ratio[1],rowsLegendOnlyPi0Ratio[0] ,"stat");
    //         SetStyleTLatex( textStatOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
    //         textStatOnlyRatioPi0->Draw();
    //         TLatex *textSysOnlyRatioPi0                 = new TLatex(columnsLegendOnlyPi0Ratio[2] ,rowsLegendOnlyPi0Ratio[0],"syst");
    //         SetStyleTLatex( textSysOnlyRatioPi0, textSizeLabelsPixel,4, 1, 43);
    //         textSysOnlyRatioPi0->Draw();
    //         TLatex *textStatOnlyRatioPi02               = new TLatex(columnsLegendOnlyPi0Ratio[4],rowsLegendOnlyPi0Ratio[0] ,"stat");
    //         SetStyleTLatex( textStatOnlyRatioPi02, textSizeLabelsPixel,4, 1, 43);
    //         textStatOnlyRatioPi02->Draw();
    //         TLatex *textSysOnlyRatioPi02                = new TLatex(columnsLegendOnlyPi0Ratio[5] ,rowsLegendOnlyPi0Ratio[0],"syst");
    //         SetStyleTLatex( textSysOnlyRatioPi02, textSizeLabelsPixel,4, 1, 43);
    //         textSysOnlyRatioPi02->Draw();
    //
    //         TMarker* markerPCMPi0OnlyRatioPi0           = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[0],columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[1],1);
    //         markerPCMPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[1]);
    //         TMarker* markerPHOSPi0OnlyRatioPi0          = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[1], columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[2],1);
    //         markerPHOSPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[2]);
    //         TMarker* markerEMCALPi0OnlyRatioPi0         = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[2], columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[3],1);
    //         markerEMCALPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[3]);
    //         TMarker* markerPCMEMCALPi0OnlyRatioPi0      = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[4], columnsLegendOnlyPi0Ratio[3] ,rowsLegendOnlyPi0Ratio[1],1);
    //         markerPCMEMCALPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[4] ,rowsLegendOnlyPi0RatioAbs[1]);
    //         TMarker* markerPCMPHOSPi0OnlyRatioPi0      = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[3], columnsLegendOnlyPi0Ratio[3] ,rowsLegendOnlyPi0Ratio[2],1);
    //         markerPCMPHOSPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[4] ,rowsLegendOnlyPi0RatioAbs[2]);
    //
    //         TBox* boxPCMPi0OnlyRatioPi0                 = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[0], columnsLegendOnlyPi0RatioAbs[2]-0.1*lengthBox , rowsLegendOnlyPi0RatioAbs[1]- heightBox,
    //                                                                          columnsLegendOnlyPi0RatioAbs[2]+ 1.7*lengthBox, rowsLegendOnlyPi0RatioAbs[1]+ heightBox);
    //         boxPCMPi0OnlyRatioPi0->Draw("l");
    //         TBox* boxPHOSPi0OnlyRatioPi0                = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[1], columnsLegendOnlyPi0RatioAbs[2]-0.1*lengthBox , rowsLegendOnlyPi0RatioAbs[2]- heightBox,
    //                                                                          columnsLegendOnlyPi0RatioAbs[2]+ 1.7*lengthBox, rowsLegendOnlyPi0RatioAbs[2]+ heightBox);
    //         boxPHOSPi0OnlyRatioPi0->Draw("l");
    //         TBox* boxEMCALPi0OnlyRatioPi0               = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[2], columnsLegendOnlyPi0RatioAbs[2]-0.1*lengthBox , rowsLegendOnlyPi0RatioAbs[3]- heightBox,
    //                                                                          columnsLegendOnlyPi0RatioAbs[2]+ 1.7*lengthBox, rowsLegendOnlyPi0RatioAbs[3]+ heightBox);
    //         boxEMCALPi0OnlyRatioPi0->Draw("l");
    //         TBox* boxPCMEMCALPi0OnlyRatioPi0            = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[4], columnsLegendOnlyPi0RatioAbs[5]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[1]- heightBox,
    //                                                                          columnsLegendOnlyPi0RatioAbs[5]+ 18*lengthBox, rowsLegendOnlyPi0RatioAbs[1]+ heightBox);
    //         boxPCMEMCALPi0OnlyRatioPi0->Draw("l");
    //         TBox* boxPCMPHOSPi0OnlyRatioPi0             = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[3], columnsLegendOnlyPi0RatioAbs[5]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[2]- heightBox,
    //                                                                          columnsLegendOnlyPi0RatioAbs[5]+ 18*lengthBox, rowsLegendOnlyPi0RatioAbs[2]+ heightBox);
    //         boxPCMPHOSPi0OnlyRatioPi0->Draw("l");
    //
    //     canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFit.%s",outputDir.Data(),suffix.Data()));
    }

    // **********************************************************************************************************************
    // ******************************************* Mass and width for pi0 ***************************************************
    // **********************************************************************************************************************

    Double_t arrayBoundariesX1_4[2];
    Double_t arrayBoundariesY1_4[3];
    Double_t relativeMarginsX[3];
    Double_t relativeMarginsY[3];
    textSizeLabelsPixel                 = 50;
    ReturnCorrectValuesForCanvasScaling(1350,1250, 1, 2,0.09, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

    TCanvas* canvasMassWidthPi0         = new TCanvas("canvasMassWidthPi0","",0,0,1350,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasMassWidthPi0,  0.13, 0.02, 0.03, 0.06);

    TPad* padWidthPi0                   = new TPad("padWidthPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthPi0->Draw();

    TPad* padMassPi0                    = new TPad("padMassPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
    padMassPi0->Draw();

    TPad* padMassLegend                 = new TPad("padMassLegend", "", 0.13, 0.36, 0.82, 0.52,-1, -1, -2);
    DrawGammaPadSettings( padMassLegend, 0., 0., 0., 0.);
    padMassLegend->SetFillStyle(0);
    padMassLegend->Draw();

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

    TH2F * histo2DAllPi0FWHM        = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, 0.23, 25. ,1000., -30, 40);
    SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                              0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
    histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,20.5);
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

    TH2F * histo2DAllPi0Mass        = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, 0.23, 25., 1000., 120., 170);
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
    TLatex *labelMassEnergy         = new TLatex(0.13,0.775,collisionSystemXeXe.Data());
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


    for (Int_t cent = 0; cent < 5; cent++){
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
            labelMassEnergy->SetText(0.13,0.775,Form("%s %s", centArray[cent].Data(), collisionSystemXeXe.Data()));
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

            DrawGammaLines(0.23, 25. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);
            labelLegendBMass->Draw();

            padMassLegend->cd();
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
            if (markerPi0Mass[0]) markerPi0Mass[0]->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[1]+ offsetMarkerYMass);
            if (markerPi0Mass[4]) markerPi0Mass[4]->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
            if (markerPi0Mass[2]) markerPi0Mass[2]->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
            if (markerPi0Mass[1]) markerPi0Mass[1]->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[1]+ offsetMarkerYMass);
            if (markerPi0Mass[3]) markerPi0Mass[3]->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);

            if (markerPi0MassMC[0]) markerPi0MassMC[0]->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass-0.01 ,rowsLegendMass[1]+ offsetMarkerYMass);
            if (markerPi0MassMC[4]) markerPi0MassMC[4]->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass-0.01 ,rowsLegendMass[2]+ offsetMarkerYMass);
            if (markerPi0MassMC[2]) markerPi0MassMC[2]->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass-0.01 ,rowsLegendMass[3]+ offsetMarkerYMass);
            if (markerPi0MassMC[1]) markerPi0MassMC[1]->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass-0.01 ,rowsLegendMass[1]+ offsetMarkerYMass);
            if (markerPi0MassMC[3]) markerPi0MassMC[3]->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass-0.01 ,rowsLegendMass[2]+ offsetMarkerYMass);

        canvasMassWidthPi0->Update();
        canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth_%s.%s",outputDir.Data(),centArrayOutput[cent].Data(), suffix.Data()));

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


    // **********************************************************************************************************************
    // ******************************** Cross section for pi0 single measurement 2.76TeV ************************************
    // **********************************************************************************************************************
    TCanvas* canvasInvYieldPi0          = new TCanvas("canvasInvYieldPi0","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasInvYieldPi0, 0.16, 0.02, 0.02, 0.08);
    canvasInvYieldPi0->SetLogx();
    canvasInvYieldPi0->SetLogy();

    TH2F * histo2DYieldPi0              = new TH2F("histo2DYieldPi0","histo2DYieldPi0",11000,0.23, 25.,1000,6e-7,3e2);
    SetStyleHistoTH2ForGraphs(histo2DYieldPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.65);
    histo2DYieldPi0->GetXaxis()->SetMoreLogLabels();
    histo2DYieldPi0->GetXaxis()->SetLabelOffset(-0.01);

    TLatex *labelEnergyXSectionPi0  = new TLatex(0.94,0.92,collisionSystemXeXe.Data());
    SetStyleTLatex( labelEnergyXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);
    TLatex *labelDetSysXSectionPi0  = new TLatex(0.94,0.88,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);

    for (Int_t cent = 0; cent < 5; cent++){
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
            labelEnergyXSectionPi0->SetText(0.94,0.92,Form("%s %s", centArray[cent].Data(), collisionSystemXeXe.Data()));
            labelEnergyXSectionPi0->Draw();
            labelDetSysXSectionPi0->Draw();


            TLegend* legendXSectionPi0      = GetAndSetLegend2(0.20, 0.13,0.40,0.13+(0.035*nTotMeas[cent]), 0.035, 1, "", 42, 0);
            for (Int_t meth = 0; meth < 11; meth++){
                if (graphIndPi0InvYieldSys[cent][meth]) legendXSectionPi0->AddEntry(graphIndPi0InvYieldSys[cent][meth],nameMeasGlobalLabel[meth],"fp");
            }
            legendXSectionPi0->Draw();

        canvasInvYieldPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems_%s.%s", outputDir.Data(), centArrayOutput[cent].Data(), suffix.Data()));

//         canvasInvYieldPi0->cd();
//         histo2DYieldPi0->Draw("copy");
//             for (Int_t meth = 10; meth> -1; meth--){
//                         if (graphIndPi0InvYieldSys[cent][meth]){
//                             graphIndPi0InvYieldSys[cent][meth]->Draw("E2same");
//                 }
//             }
//             DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys[cent], markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
//             graphCombPi0InvYieldSys[cent]->Draw("E2same");
//             for (Int_t meth = 10; meth> -1; meth--){
//                 if (graphIndPi0InvYieldSys[cent][meth]){
//                     graphIndPi0InvYieldStat[cent][meth]->Draw("p,same,z");
//                 }
//             }
//             DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
//             graphCombPi0InvYieldStat[cent]->Draw("p,same,z");
//
//             labelEnergyXSectionPi0->Draw();
//             labelDetSysXSectionPi0->Draw();
//
//             legendXSectionPi0->AddEntry(graphCombPi0InvYieldSys[cent],"comb","fp");
//             legendXSectionPi0->Draw();
//
//         canvasInvYieldPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems_Comb_%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), suffix.Data()));

    }
    // **********************************************************************************************************************
    // ******************************** Acceptance * Efficiency for pi0 single measurement **********************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 55;
    Double_t textSizeLabelsRel          = 55./1200;
    cout << textSizeLabelsRel << endl;

    TCanvas* canvasAcceptanceTimesEff   = new TCanvas("canvasAcceptanceTimesEff", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasAcceptanceTimesEff,  0.1, 0.01, 0.015, 0.095);
    canvasAcceptanceTimesEff->SetLogy(1);
    canvasAcceptanceTimesEff->SetLogx(1);

    TH1F * histo1DAccEff            = new TH1F("histo1DAccEff", "histo1DAccEff",1000, 0.23,  31);
    SetStyleHistoTH1ForGraphs(  histo1DAccEff, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec}"),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
    histo1DAccEff->GetYaxis()->SetRangeUser(3e-4, 2e-0 );
    histo1DAccEff->GetYaxis()->SetLabelOffset(0.001);
    histo1DAccEff->GetXaxis()->SetLabelOffset(-0.01);
    histo1DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
    TLatex *labelPerfEffi           = new TLatex(0.15,0.92,"ALICE performance");
    SetStyleTLatex( labelPerfEffi, textSizeLabelsRel,4);
    TLatex *labelEnergyEffi         = new TLatex(0.15,0.87,collisionSystemXeXe.Data());
    SetStyleTLatex( labelEnergyEffi, textSizeLabelsRel,4);
    TLatex *labelDetSysEffiPi0      = new TLatex(0.15,0.82,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysEffiPi0, textSizeLabelsRel,4);

    for (Int_t cent = 0; cent < 5; cent++){
        canvasAcceptanceTimesEff->cd();
        histo1DAccEff->DrawCopy();

            for (Int_t meth = 0; meth < 11; meth++){
                if (graphPi0EffTimesAcc[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphPi0EffTimesAcc[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.55, colorDet[meth] , colorDet[meth]);
                    graphPi0EffTimesAcc[cent][meth]->Draw("p,same,z");
                }
            }

            TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(nTotMeas[cent]*textSizeLabelsRel),textSizeLabelsPixel);
            if (graphPi0EffTimesAcc[cent][0]) legendEffiAccPi0->AddEntry(graphPi0EffTimesAcc[cent][0],nameMeasGlobalLabel[0],"p");
            if (graphPi0EffTimesAcc[cent][4]) legendEffiAccPi0->AddEntry(graphPi0EffTimesAcc[cent][4],nameMeasGlobalLabel[4],"p");
            if (graphPi0EffTimesAcc[cent][2]) legendEffiAccPi0->AddEntry(graphPi0EffTimesAcc[cent][2],nameMeasGlobalLabel[2],"p");
            if (graphPi0EffTimesAcc[cent][1]) legendEffiAccPi0->AddEntry(graphPi0EffTimesAcc[cent][1],nameMeasGlobalLabel[1],"p");
            if (graphPi0EffTimesAcc[cent][3]) legendEffiAccPi0->AddEntry(graphPi0EffTimesAcc[cent][3],nameMeasGlobalLabel[3],"p");
            legendEffiAccPi0->Draw();

            labelPerfEffi->Draw();
            labelEnergyEffi->SetText(0.15,0.87,Form("%s %s", centArray[cent].Data(), collisionSystemXeXe.Data()));
            labelEnergyEffi->Draw();
            labelDetSysEffiPi0->Draw();

        canvasAcceptanceTimesEff->Update();
        canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff_%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), suffix.Data()));
    }

    Double_t rangMinEffi[11]                = { 3e-4, 1e-3, 8e-4, 3e-4, 3e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4 };
    Double_t rangMaxEffi[11]                = { 5e-2, 3e-1, 2e-0, 5e-2, 3e-1, 2e-0, 2e-0, 2e-0, 2e-0, 2e-0, 2e-0 };

    for (Int_t meth = 0; meth < 11; meth++){
        if (!fileNamesMethod[meth].CompareTo("")) continue;
        canvasAcceptanceTimesEff->cd();
        histo1DAccEff->GetYaxis()->SetRangeUser(rangMinEffi[meth], rangMaxEffi[meth] );
        histo1DAccEff->DrawCopy();

        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(5*textSizeLabelsRel),textSizeLabelsPixel);
        for (Int_t cent = 0; cent < 5; cent++){
            if (graphPi0EffTimesAcc[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphPi0EffTimesAcc[cent][meth], markerStyleCent[cent], markerSizeCent[cent]*0.55, colorCent[cent] , colorCent[cent]);
                graphPi0EffTimesAcc[cent][meth]->Draw("p,same,z");
                legendEffiAccPi0->AddEntry(graphPi0EffTimesAcc[cent][meth],centArray[cent],"p");
            }
        }
        legendEffiAccPi0->Draw();

        labelPerfEffi->Draw();
        labelEnergyEffi->SetText(0.15,0.87,Form("%s", collisionSystemXeXe.Data()));
        labelEnergyEffi->Draw();
        labelDetSysEffiPi0->SetText(0.15,0.82,Form("#pi^{0} #rightarrow #gamma#gamma, %s", nameMeasGlobalLabel[meth].Data()));
        labelDetSysEffiPi0->Draw();

        canvasAcceptanceTimesEff->Update();
        canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff_%s.%s",outputDir.Data(), nameMeasGlobalLabel[meth].Data(), suffix.Data()));

    }
    delete labelPerfEffi;
    delete labelEnergyEffi;
    delete labelDetSysEffiPi0;
    delete canvasAcceptanceTimesEff;

    // **********************************************************************************************************************
    // ******************************** effective secondary correction drawing for different methods ************************
    // **********************************************************************************************************************
    TCanvas* canvasEffectiveSecCorr     = new TCanvas("canvasEffectiveSecCorr", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasEffectiveSecCorr,  0.1, 0.01, 0.04, 0.095);
    canvasEffectiveSecCorr->SetLogx(1);

        TH1F * histo1DEffSecCorr;
        histo1DEffSecCorr               = new TH1F("histo1DEffSecCorr", "histo1DEffSecCorr",1000, 0.23,  25);
        SetStyleHistoTH1ForGraphs(  histo1DEffSecCorr, "#it{p}_{T} (GeV/#it{c})","R_{K^{0}_{s}}",
                                    0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
        histo1DEffSecCorr->GetYaxis()->SetRangeUser(0, 10 );
        histo1DEffSecCorr->GetXaxis()->SetLabelOffset(-0.01);
        histo1DEffSecCorr->GetXaxis()->SetMoreLogLabels(kTRUE);
        TLatex *labelPerfSecCorr        = new TLatex(0.15,0.90,"ALICE performance");
        SetStyleTLatex( labelPerfSecCorr, textSizeLabelsRel,4);
        TLatex *labelEnergySecCorr      = new TLatex(0.15,0.85,collisionSystemXeXe.Data());
        SetStyleTLatex( labelEnergySecCorr, textSizeLabelsRel,4);
        TLatex *labelDetSysSecCorrPi0   = new TLatex(0.15,0.80,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysSecCorrPi0, textSizeLabelsRel,4);

        for (Int_t cent= 0; cent < 5; cent++ ){
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
                    canvasEffectiveSecCorr->Print(Form("%s/Pi0_EffectiveSecCorr_%s_%s.%s",outputDir.Data(), nameSecPi0SourceRead[k].Data(), centArrayOutput[cent].Data(), suffix.Data()));
                }
            }
        }
    delete canvasEffectiveSecCorr;

//     fileFitsOutput << "*******************************************************************************************" << endl;
//     fileFitsOutput << "****************************** Power law fit pi0 ******************************************" << endl;
//     fileFitsOutput << "*******************************************************************************************" << endl;
//     TF1* fitPowInvYieldPi0Tot   = FitObject("powPure","fitPowInvYieldPi0Tot","Pi0",graphCombPi0InvYieldTot,4,20. ,NULL,"QNRMEX0+","", kFALSE);
//     cout << WriteParameterToFile(fitPowInvYieldPi0Tot)<< endl;
//     fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldPi0Tot)<< endl;
//     TF1* fitPowInvYieldPi0Stat   = FitObject("powPure","fitPowInvYieldPi0","Pi0",graphCombPi0InvYieldStat,4,20. ,NULL,"QNRMEX0+","", kFALSE);
//     cout << WriteParameterToFile(fitPowInvYieldPi0Stat)<< endl;
//     fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldPi0Stat)<< endl;
//     TF1* fitOHagInvYieldPi0Tot   = FitObject("oHag","fitOHagInvYieldPi0","Pi0",graphCombPi0InvYieldTot,0.3,20. ,NULL,"QNRMEX0+","", kFALSE);
//     cout << WriteParameterToFile(fitOHagInvYieldPi0Tot)<< endl;
//     fileFitsOutput <<  WriteParameterToFile(fitOHagInvYieldPi0Tot)<< endl;
//
//     fileFitsOutput << "NLO pp" << endl;
//     TF1* fitPowInvYieldPi0NLOPP         = FitObject("powPure","fitPowInvYieldPi0NLOPP","Pi0",graphNLODSS14Pi0PP,4,20. ,NULL,"QNRMEX0+","", kFALSE);
//     fileFitsOutput << WriteParameterToFile(fitPowInvYieldPi0NLOPP)<< endl;

//     fileFitsOutput << "NLO XeXe nCTEQ" << endl;
//     TF1* fitPowInvYieldPi0NLOXeXenCTEQ   = FitObject("powPure","fitPowInvYieldPi0NLOXeXenCTEQ","Pi0",graphNLODSS14nPDFPi0,4,20. ,NULL,"QNRMEX0+","", kFALSE);
//     fileFitsOutput << WriteParameterToFile(fitPowInvYieldPi0NLOXeXenCTEQ)<< endl;
//
//     fileFitsOutput << "NLO pp EPPS16" << endl;
//     TF1* fitPowInvYieldPi0NLOXeXeEPPS    = FitObject("powPure","fitPowInvYieldPi0NLOXeXeEPPS","Pi0",graphNLODSS14nPDFEPPS16Pi0,4,20. ,NULL,"QNRMEX0+","", kFALSE);
//     fileFitsOutput << WriteParameterToFile(fitPowInvYieldPi0NLOXeXeEPPS)<< endl;
//

//     // **********************************************************************************************************************
//     // ******************************************* Comparison to theory calculations Pi0 ************************************
//     // **********************************************************************************************************************
//     textSizeLabelsPixel                     = 48;
//
//     TCanvas* canvasRatioXeXe                 = new TCanvas("canvasRatioXeXe","",200,10,1350,900);  // gives the page size
//     DrawGammaCanvasSettings( canvasRatioXeXe,  0.12, 0.01, 0.01, 0.11);
//     canvasRatioXeXe->cd();
//     canvasRatioXeXe->SetLogx();
//
//         textsizeLabelsXeXe                   = 0;
//         if (canvasRatioXeXe->XtoPixel(canvasRatioXeXe->GetX2()) <canvasRatioXeXe->YtoPixel(canvasRatioXeXe->GetY1()) ){
//             textsizeLabelsXeXe               = (Double_t)textSizeLabelsPixel/canvasRatioXeXe->XtoPixel(canvasRatioXeXe->GetX2()) ;
//         } else {
//             textsizeLabelsXeXe               = (Double_t)textSizeLabelsPixel/canvasRatioXeXe->YtoPixel(canvasRatioXeXe->GetY1());
//         }
//         cout << textsizeLabelsXeXe << endl;
//
//         TH2F * ratio2DTheoryXeXe             = new TH2F("ratio2DTheoryXeXe","ratio2DTheoryXeXe",1000,0.23, 25.,1000,0.2,2.95);
//         SetStyleHistoTH2ForGraphs(  ratio2DTheoryXeXe, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsXeXe, textsizeLabelsXeXe,
//                                     0.85*textsizeLabelsXeXe,textsizeLabelsXeXe, 0.9, 0.95, 510, 505);
//         ratio2DTheoryXeXe->GetYaxis()->SetMoreLogLabels(kTRUE);
//         ratio2DTheoryXeXe->GetYaxis()->SetNdivisions(510);
//         ratio2DTheoryXeXe->GetYaxis()->SetNoExponent(kTRUE);
//         ratio2DTheoryXeXe->GetXaxis()->SetMoreLogLabels(kTRUE);
//         ratio2DTheoryXeXe->GetXaxis()->SetNoExponent(kTRUE);
//         ratio2DTheoryXeXe->GetXaxis()->SetLabelFont(42);
//         ratio2DTheoryXeXe->GetYaxis()->SetLabelFont(42);
//         ratio2DTheoryXeXe->DrawCopy();
//
//         TGraphAsymmErrors* graphRatioPi0CombCombFitStatWOXErr   = (TGraphAsymmErrors*)graphRatioPi0CombCombFitStat->Clone("graphRatioPi0CombCombFitStatWOXErr");
//         ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatWOXErr);
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
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
//         graphRatioPi0CombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
//         graphRatioPi0CombCombFitSys->SetLineWidth(0);
//         graphRatioPi0CombCombFitSys->Draw("2,same");
//         graphRatioPi0CombCombFitStatWOXErr->Draw("p,same");
//
//         TBox* boxErrorNormRatioPi0          = CreateBoxConv(kGray+1, 0.28, 1.-(0.031 ), 0.32, 1.+(0.031));
//         boxErrorNormRatioPi0->Draw();
//         DrawGammaLines(0.23, 25.,1., 1.,0.1,kGray);
//
//         TLegend* legendRatioTheoryXeXeNew= GetAndSetLegend2(0.69,0.91-0.85*textsizeLabelsXeXe*6,0.94,0.91, 0.85*textSizeLabelsPixel);
//         legendRatioTheoryXeXeNew->AddEntry(graphRatioPi0CombCombFitSys,"#pi^{0} ALICE","pf");
// //         legendRatioTheoryXeXeNew->AddEntry(histoRatioPi0DPMJetToFit,  "DPMJet", "l");
// //         legendRatioTheoryXeXeNew->AddEntry(histoRatioPi0HIJINGToFit,  "HIJING", "l");
// //         legendRatioTheoryXeXeNew->AddEntry(histoRatioPi0EPOS3ToFit,  "EPOS", "l");
// //         legendRatioTheoryXeXeNew->AddEntry(graphRatioPi0McGillToFit,  "Shen #it{et al.}", "f");
//         legendRatioTheoryXeXeNew->Draw();
//
//         TLatex *labelRatioTheory            = new TLatex(0.97,0.925,collisionSystemXeXe.Data());
//         SetStyleTLatex( labelRatioTheory, 0.85*textsizeLabelsXeXe,4, 1, 42, kTRUE, 31);
//         labelRatioTheory->Draw();
//
//         ratio2DTheoryXeXe->Draw("same,axis");
//     canvasRatioXeXe->Update();
//     canvasRatioXeXe->Print(Form("%s/Pi0_RatioTheoryToData_XeXe.%s",outputDir.Data(),suffix.Data()));
//
//
//     //*************************************************************************************************************
//     //***************************** Paper plot invyield and ratios ************************************************
//     //*************************************************************************************************************
//
//     Double_t arrayBoundariesX1_XSec[2];
//     Double_t arrayBoundariesY1_XSec[6];
//     Double_t relativeMarginsXXSec[3];
//     Double_t relativeMarginsYXSec[3];
//     textSizeLabelsPixel = 48;
//     ReturnCorrectValuesForCanvasScaling(1250,2000, 1, 5,0.15, 0.005, 0.003,0.05,arrayBoundariesX1_XSec,arrayBoundariesY1_XSec,relativeMarginsXXSec,relativeMarginsYXSec);
//
//     TCanvas* canvasInvSectionPaper      = new TCanvas("canvasInvSectionPaper","",0,0,1250,2000);  // gives the page size
//     DrawGammaCanvasSettings( canvasInvSectionPaper,  0.13, 0.02, 0.03, 0.06);
//
//     TPad* padInvSectionSpec             = new TPad("padInvSectionSpec", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[3], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[0],-1, -1, -2);
//     DrawGammaPadSettings( padInvSectionSpec, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[0], relativeMarginsYXSec[1]);
//     padInvSectionSpec->Draw();
//     Double_t marginXSec                 = relativeMarginsXXSec[0]*1250;
//     Double_t textsizeLabelsXSecUp       = 0;
//     Double_t textsizeFacXSecUp          = 0;
//     if (padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) < padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1())){
//         textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
//         textsizeFacXSecUp               = (Double_t)1./padInvSectionSpec->XtoPixel(padInvSectionSpec->GetX2()) ;
//     } else {
//         textsizeLabelsXSecUp            = (Double_t)textSizeLabelsPixel/padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
//         textsizeFacXSecUp                   = (Double_t)1./padInvSectionSpec->YtoPixel(padInvSectionSpec->GetY1());
//     }
//
//     TPad* padInvSectionNLORatio         = new TPad("padInvSectionNLORatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[4], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[3],-1, -1, -2);
//     DrawGammaPadSettings( padInvSectionNLORatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[1]);
//     padInvSectionNLORatio->Draw();
//     Double_t textsizeLabelsXSecMiddle   = 0;
//     Double_t textsizeFacXSecMiddle      = 0;
//     if (padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) < padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1())){
//         textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) ;
//         textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatio->XtoPixel(padInvSectionNLORatio->GetX2()) ;
//     } else {
//         textsizeLabelsXSecMiddle        = (Double_t)textSizeLabelsPixel/padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1());
//         textsizeFacXSecMiddle           = (Double_t)1./padInvSectionNLORatio->YtoPixel(padInvSectionNLORatio->GetY1());
//     }
//
//     TPad* padInvSectionPythiaRatio      = new TPad("padInvSectionPythiaRatio", "", arrayBoundariesX1_XSec[0], arrayBoundariesY1_XSec[5], arrayBoundariesX1_XSec[1], arrayBoundariesY1_XSec[4],-1, -1, -2);
//     DrawGammaPadSettings( padInvSectionPythiaRatio, relativeMarginsXXSec[0], relativeMarginsXXSec[2], relativeMarginsYXSec[1], relativeMarginsYXSec[2]);
//     padInvSectionPythiaRatio->Draw();
//     Double_t textsizeLabelsXSecDown     = 0;
//     Double_t textsizeFacXSecDown        = 0;
//     if (padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) < padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1())){
//         textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) ;
//         textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatio->XtoPixel(padInvSectionPythiaRatio->GetX2()) ;
//     } else {
//         textsizeLabelsXSecDown          = (Double_t)textSizeLabelsPixel/padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1());
//         textsizeFacXSecDown             = (Double_t)1./padInvSectionPythiaRatio->YtoPixel(padInvSectionPythiaRatio->GetY1());
//     }
//
//
//     //*************************************************************************************************************
//     //***************************** Paper plot invyield pi0 with soft-physics *************************************
//     //*************************************************************************************************************
//     padInvSectionSpec->cd();
//     padInvSectionSpec->SetLogy(1);
//     padInvSectionSpec->SetLogx(1);
//         SetStyleHistoTH2ForGraphs(histo2DYieldPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
//                                 0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.257/(textsizeFacXSecUp*marginXSec));
//         histo2DYieldPi0->GetXaxis()->SetMoreLogLabels();
//         histo2DYieldPi0->GetXaxis()->SetLabelOffset(+0.01);
//         histo2DYieldPi0->Draw();
//
// //         SetStyleHisto(histoEPOS3Pi0, widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );
// //         DrawGammaSetMarkerTGraphErr(graphEPOS3Pi0, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
// //         histoEPOS3Pi0->GetXaxis()->SetRangeUser(0.3,20);
// //         histoEPOS3Pi0->Draw("same,hist,l");
// //         while(graphEPOS3Pi0->GetX()[0] < 0.3)
// //             graphEPOS3Pi0->RemovePoint(0);
// //         while(graphEPOS3Pi0->GetX()[graphEPOS3Pi0->GetN()-1] > 20)
// //             graphEPOS3Pi0->RemovePoint(graphEPOS3Pi0->GetN()-1);
// //         graphEPOS3Pi0->Draw("same,3");
//
// //         DrawGammaSetMarkerTGraphErr(graphMcGillPi0, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
// //         while(graphMcGillPi0->GetX()[0] < 0.3)
// //             graphMcGillPi0->RemovePoint(0);
// //         graphMcGillPi0->Draw("same,3");
//
// //         SetStyleHisto(histoDPMJetPi0, widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );
// //         histoDPMJetPi0->GetXaxis()->SetRangeUser(0.3,20);
// //         histoDPMJetPi0->Draw("same,hist,l");
// //         SetStyleHisto(histoHIJINGPi0, widthCommonFit*1.5, styleLineHIJING, colorHIJING );
// //         histoHIJINGPi0->GetXaxis()->SetRangeUser(0.3,20);
// //         histoHIJINGPi0->Draw("same,hist,l");
//
//         TGraphAsymmErrors* graphCombPi0InvYieldStatWOXErr = (TGraphAsymmErrors*)graphCombPi0InvYieldStat->Clone("graphCombPi0InvYieldStatWOXErr");
//         ProduceGraphAsymmWithoutXErrors(graphCombPi0InvYieldStatWOXErr);
//
//         DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
//         graphCombPi0InvYieldSys->Draw("E2same");
//         DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
//         graphCombPi0InvYieldStatWOXErr->Draw("p,same,z");
//
//         DrawGammaSetMarkerTF1( fitTCMInvYieldPi0, 7, 2, kGray+2);
//         fitTCMInvYieldPi0->Draw("same");
//
// //         DrawGammaSetMarkerTF1( fitPowInvYieldPi0, 7, 2, kBlue-7);
// //         fitPowInvYieldPi0->Draw("same");
//
//         TLatex *labelEnergyXSectionPaper= new TLatex(0.95, 0.91, collisionSystemXeXe.Data());
//         SetStyleTLatex( labelEnergyXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
//         labelEnergyXSectionPaper->Draw();
//         TLatex *labelALICEXSectionPaper= new TLatex(0.95,0.867,textALICE.Data());
//         SetStyleTLatex( labelALICEXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
//         labelALICEXSectionPaper->Draw();
//         TLatex *labelDetSysXSectionPaper= new TLatex(0.95,0.83,"#pi^{0} #rightarrow #gamma#gamma");
//         SetStyleTLatex( labelDetSysXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
//         labelDetSysXSectionPaper->Draw();
//
//         TGraphErrors* dummyNorm     = new TGraphErrors(1);
//         dummyNorm->SetPoint(0,1,0.0251);
//         DrawGammaSetMarkerTGraphErr(dummyNorm, 0, 0, kGray+1, kGray+1, widthLinesBoxes, kTRUE, kGray+1);
//
//         TLegend* legendXsectionPaper    = GetAndSetLegend2(0.19, 0.03, 0.5, 0.03+0.05*7, textSizeLabelsPixel, 1, "", 43, 0.2);
//         legendXsectionPaper->AddEntry(graphCombPi0InvYieldSys,"Data","pf");
//         legendXsectionPaper->AddEntry(dummyNorm,"Norm. unc. 3.1%","f");
//         legendXsectionPaper->AddEntry(fitTCMInvYieldPi0,"TCM fit","l");
// //         legendXsectionPaper->AddEntry(histoDPMJetPi0,"DPMJet","l");
// //         legendXsectionPaper->AddEntry(histoHIJINGPi0,"HIJING","l");
// //         legendXsectionPaper->AddEntry(histoEPOS3Pi0,"EPOS3","l");
// //         legendXsectionPaper->AddEntry(graphMcGillPi0,"Shen #it{et al.}","f");
// //         legendXsectionPaper->AddEntry(graphMcGillPi0,"iEBE-VISHNU","f");
//         legendXsectionPaper->Draw();
//
//
//     padInvSectionNLORatio->cd();
//     padInvSectionNLORatio->SetLogx(1);
//         TH2F * ratio2DEPOS               = new TH2F("ratio2DEPOS","ratio2DEPOS",1000,0.23, 25.,1000,0.3,2.2);
//         SetStyleHistoTH2ForGraphs(ratio2DEPOS, "#it{p}_{T} (GeV/#it{c})","#frac{EPOS, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
//                                   0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.257/(textsizeFacXSecMiddle*marginXSec), 510, 508);
//         ratio2DEPOS->GetYaxis()->SetNoExponent(kTRUE);
//         ratio2DEPOS->GetXaxis()->SetMoreLogLabels(kTRUE);
//         ratio2DEPOS->GetXaxis()->SetNoExponent(kTRUE);
//         ratio2DEPOS->GetXaxis()->SetLabelFont(42);
//         ratio2DEPOS->GetYaxis()->SetLabelFont(42);
//         ratio2DEPOS->GetYaxis()->SetLabelOffset(+0.01);
//         ratio2DEPOS->GetXaxis()->SetTickLength(0.07);
//         ratio2DEPOS->DrawCopy();
//
//         boxErrorNormRatioPi0->Draw();
//
// //         SetStyleHisto(histoRatioPi0EPOS3ToFit, widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );
// //         DrawGammaSetMarkerTGraphErr(graphRatioPi0EPOS3ToFit, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
// //         histoRatioPi0EPOS3ToFit->Draw("same,hist,l");
// //         graphRatioPi0EPOS3ToFit->Draw("same,3");
// //         DrawGammaSetMarkerTGraphErr(graphRatioPi0McGillToFit, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
// //         graphRatioPi0McGillToFit->Draw("same,3");
//
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
//         graphRatioPi0CombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
//         graphRatioPi0CombCombFitSys->SetLineWidth(0);
//         graphRatioPi0CombCombFitSys->Draw("2,same");
//         graphRatioPi0CombCombFitStatWOXErr->Draw("p,same");
//
//         DrawGammaLines(0.23, 25.,1., 1.,0.1,kGray);
//
//     padInvSectionPythiaRatio->cd();
//     padInvSectionPythiaRatio->SetLogx(1);
//         TH2F * ratio2DGenerator            = new TH2F("ratio2DGenerator","ratio2DGenerator",1000,0.23, 25.,1000,0.3,2.2);
//         SetStyleHistoTH2ForGraphs(ratio2DGenerator, "#it{p}_{T} (GeV/#it{c})","#frac{Generator, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown,
//                                   0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.257/(textsizeFacXSecDown*marginXSec), 510, 508);
//         ratio2DGenerator->GetYaxis()->SetNoExponent(kTRUE);
//         ratio2DGenerator->GetXaxis()->SetMoreLogLabels(kTRUE);
//         ratio2DGenerator->GetXaxis()->SetNoExponent(kTRUE);
//         ratio2DGenerator->GetXaxis()->SetLabelFont(42);
//         ratio2DGenerator->GetYaxis()->SetLabelFont(42);
//         ratio2DGenerator->GetYaxis()->SetLabelOffset(+0.01);
//         ratio2DGenerator->GetXaxis()->SetTickLength(0.06);
//         ratio2DGenerator->GetYaxis()->SetTickLength(0.04);
//         ratio2DGenerator->DrawCopy();
//
// //         SetStyleHisto(histoRatioPi0DPMJetToFit, widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );
// //         histoRatioPi0DPMJetToFit->GetXaxis()->SetRangeUser(0.5,14);
// //         histoRatioPi0DPMJetToFit->Draw("same,hist,l");
// //         SetStyleHisto(histoRatioPi0HIJINGToFit, widthCommonFit*1.5, styleLineHIJING, colorHIJING );
// //         histoRatioPi0HIJINGToFit->GetXaxis()->SetRangeUser(0.5,14);
// //         histoRatioPi0HIJINGToFit->Draw("same,hist,l");
//
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
//         graphRatioPi0CombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
//         graphRatioPi0CombCombFitSys->SetLineWidth(0);
//         graphRatioPi0CombCombFitSys->Draw("2,same");
//         graphRatioPi0CombCombFitStatWOXErr->Draw("p,same");
//
//         boxErrorNormRatioPi0->Draw();
//
//         DrawGammaLines(0.23, 25.,1., 1.,0.1,kGray);
//
//     canvasInvSectionPaper->Print(Form("%s/Pi0_InvYieldWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));
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
//         DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
//         graphCombPi0InvYieldSys->Draw("E2same");
//         DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombPi0Sys, markerStyleComb+4, markerSizeComb+0.2, kGray+1, kGray+1, widthLinesBoxes, kTRUE);
//         graphPPInvYieldCombPi0Sys->Draw("E2same");
//         DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
//         graphCombPi0InvYieldStatWOXErr->Draw("p,same,z");
//
//         TGraphAsymmErrors* graphPPInvYieldCombPi0StatWOXErr = (TGraphAsymmErrors*)graphPPInvYieldCombPi0Stat->Clone("graphPPInvYieldCombPi0StatWOXErr");
//         ProduceGraphAsymmWithoutXErrors(graphPPInvYieldCombPi0StatWOXErr);
//
//         DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombPi0StatWOXErr, markerStyleComb+4, markerSizeComb+0.2, kGray+1, kGray+1);
//         graphPPInvYieldCombPi0StatWOXErr->Draw("p,same,z");
//
// //         DrawGammaSetMarkerTF1( fitPowInvYieldPi0NLOXeXenCTEQ, 7, 4, kOrange);
// //         fitPowInvYieldPi0NLOXeXenCTEQ->Draw("same");
//
//         DrawGammaSetMarkerTF1( fitTCMInvYieldPi0, 7, 2, kGray+2);
//         fitTCMInvYieldPi0->Draw("same");
//         DrawGammaSetMarkerTF1( fitPPTCMInvYieldPi0, 4, 2, kGray+1);
//         fitPPTCMInvYieldPi0->Draw("same");
//
//         TLatex *labelALICENLOPaper= new TLatex(0.19,0.92,textALICE.Data());
//         SetStyleTLatex( labelALICENLOPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 11);
//         labelALICENLOPaper->Draw();
//
//         TLegend* legendPaperPi0WithNLO    = GetAndSetLegend2(0.19, 0.03, 0.5, 0.03+0.05*8, textSizeLabelsPixel, 1, collisionSystemXeXe.Data(), 43, 0.2);
// //         legendPaperPi0WithNLO->AddEntry(graphCombPi0InvYieldSys,Form("%s %s","#pi^{0}", collisionSystemXeXe.Data()),"pf");
//         legendPaperPi0WithNLO->AddEntry(graphCombPi0InvYieldSys,"#pi^{0}","pf");
//         legendPaperPi0WithNLO->AddEntry(dummyNorm,"Norm. unc. 3.1%","f");
//         legendPaperPi0WithNLO->AddEntry(fitTCMInvYieldPi0,"TCM fit","l");
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
//     TH2F * ratio2DNLOXeXe               = new TH2F("ratio2DNLOXeXe","ratio2DNLOXeXe",1000,0.23, 25.,1000,0.3,2.2);
//     SetStyleHistoTH2ForGraphs(ratio2DNLOXeXe, "#it{p}_{T} (GeV/#it{c})","#frac{NLO, p-Pb Data}{p-Pb fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle,
//                               0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.257/(textsizeFacXSecMiddle*marginXSec), 510, 508);
//     //         ratio2DNLOXeXe->GetYaxis()->SetNdivisions(512);
//     ratio2DNLOXeXe->GetYaxis()->SetNoExponent(kTRUE);
//     ratio2DNLOXeXe->GetXaxis()->SetMoreLogLabels(kTRUE);
//     ratio2DNLOXeXe->GetXaxis()->SetNoExponent(kTRUE);
//     ratio2DNLOXeXe->GetXaxis()->SetLabelFont(42);
//     ratio2DNLOXeXe->GetYaxis()->SetLabelFont(42);
//     ratio2DNLOXeXe->GetYaxis()->SetLabelOffset(+0.01);
//     ratio2DNLOXeXe->GetXaxis()->SetTickLength(0.07);
//     ratio2DNLOXeXe->DrawCopy();
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
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
//         graphRatioPi0CombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
//         graphRatioPi0CombCombFitSys->SetLineWidth(0);
//         graphRatioPi0CombCombFitSys->Draw("2,same");
//         graphRatioPi0CombCombFitStatWOXErr->Draw("p,same");
//
// //         graphRatioPi0CGC->Draw("same,c");
// //         graphRatioPi0NLODSS14nPDFEPPS16Center->Draw("same,c");
// //         graphRatioPi0NLODSS14nPDFCenter->Draw("same,c");
//
//         TLatex *labelEnergyRatioPaperNLOXeXe= new TLatex(0.19, 1-textsizeLabelsXSecMiddle*1.5, collisionSystemXeXe.Data());
//         SetStyleTLatex( labelEnergyRatioPaperNLOXeXe, textsizeLabelsXSecMiddle,4, 1, 42, kTRUE, 11);
//         labelEnergyRatioPaperNLOXeXe->Draw();
//
//         ratio2DNLOXeXe->Draw("same,axis");
//         DrawGammaLines(0.23, 25.,1., 1.,0.1,kGray);
//
//     cout << __LINE__ << endl;
//
//     padInvSectionPythiaRatio->cd();
//     padInvSectionPythiaRatio->SetLogx(1);
//     TH2F * ratio2DNLOpp            = new TH2F("ratio2DNLOpp","ratio2DNLOpp",1000,0.23, 25.,1000,0.3,2.2);
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
//         DrawGammaLines(0.23, 25.,1., 1.,0.1,kGray);
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
//     TGraphAsymmErrors* graphRatioPi0TCMTot      = (TGraphAsymmErrors*)graphCombPi0InvYieldTot->Clone();
//     graphRatioPi0TCMTot                         = CalculateGraphErrRatioToFit(graphRatioPi0TCMTot, fitTCMInvYieldPi0);
//     TGraphAsymmErrors* graphRatioPi0TsallisTot  = (TGraphAsymmErrors*)graphCombPi0InvYieldTot->Clone();
//     graphRatioPi0TsallisTot                     = CalculateGraphErrRatioToFit(graphRatioPi0TsallisTot, fitInvYieldPi0);
//     TGraphAsymmErrors* graphRatioPi0PowTot      = (TGraphAsymmErrors*)graphCombPi0InvYieldTot->Clone();
//     graphRatioPi0PowTot                         = CalculateGraphErrRatioToFit(graphRatioPi0PowTot, fitPowInvYieldPi0Tot);
//     TGraphAsymmErrors* graphRatioPi0OHagTot     = (TGraphAsymmErrors*)graphCombPi0InvYieldTot->Clone();
//     graphRatioPi0OHagTot                        = CalculateGraphErrRatioToFit(graphRatioPi0OHagTot, fitOHagInvYieldPi0Tot);
//     // restricting pT range of power law fit drawing
//     while (graphRatioPi0PowTot->GetX()[0] < 2.) graphRatioPi0PowTot->RemovePoint(0);
//
//
//     TH2F * ratio2DFits       = new TH2F("ratio2DFits","ratio2DFits",1000,0.23, 25.,1000,0.2,2.3);
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
//         DrawGammaLines(0.23, 25. , 1.1, 1.1,1, kGray, 3);
//         DrawGammaLines(0.23, 25. , 1, 1,1, kGray, 7);
//         DrawGammaLines(0.23, 25. , 0.9, 0.9,1, kGray, 3);
//
//         TLegend* legendRatioPi0Fits= GetAndSetLegend2(0.12,0.95-4*1.05*textsizeLabelsFits,0.37,0.95, textSizeLabelsPixel, 1, "", 43, 0.2);
//         legendRatioPi0Fits->AddEntry(graphRatioPi0TCMTot,"TCM","p");
//         legendRatioPi0Fits->AddEntry(graphRatioPi0TsallisTot,"Levy-Tsallis","p");
//         legendRatioPi0Fits->AddEntry(graphRatioPi0OHagTot,"mod. Hagedorn","p");
//         legendRatioPi0Fits->AddEntry(graphRatioPi0PowTot,"pure powerlaw, 4-20 GeV/#it{c}","p");
//         legendRatioPi0Fits->Draw();
//
//         TLatex *labelRatioDiffFitsEnergy    = new TLatex(0.95,0.91,collisionSystemXeXe.Data());
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


    // **********************************************************************************************************************
    // ************************* Saving of final results ********************************************************************
    // **********************************************************************************************************************

    TString nameOutputCommonFile    = Form("CombinedResultsPaperXeXe5TeV_%s.root", dateForOutput.Data());

    TFile fCombResults(nameOutputCommonFile.Data(), "UPDATE");

    for (Int_t cent = 0; cent< 5; cent++){
        TDirectoryFile* directoryPi0Output  = NULL;
        directoryPi0Output                  = (TDirectoryFile*)fCombResults.Get(Form("Pi0XeXe5TeV_%s", centArray[cent].Data()));
        if (!directoryPi0Output){
            fCombResults.mkdir(Form("Pi0XeXe5TeV_%s", centArray[cent].Data()));
            directoryPi0Output              = (TDirectoryFile*)fCombResults.Get(Form("Pi0XeXe5TeV_%s", centArray[cent].Data()));
        }
        fCombResults.cd(Form("Pi0XeXe5TeV_%s", centArray[cent].Data()));

    //         // Final spectrum
    //         graphCombPi0InvYieldTot[cent]->Write("graphInvYieldPi0CombXeXe5TeVTotErr",TObject::kOverwrite);
    //         graphCombPi0InvYieldStat[cent]->Write("graphInvYieldPi0CombXeXe5TeVStatErr",TObject::kOverwrite);
    //         graphCombPi0InvYieldSys[cent]->Write("graphInvYieldPi0CombXeXe5TeVSysErr",TObject::kOverwrite);
    //          // fits for pi0
    //         fitInvYieldPi0[cent]->Write("TsallisFitPi0");
    //         fitTCMInvYieldPi0[cent]->Write("TwoComponentModelFitPi0",TObject::kOverwrite);
    //
    //         // Final inv yield INEL
    //         if (bWCorrection.Contains("Y")){
    //             // Final spectrum correlations Method A
    //             if(graphCombPi0InvYieldTot_yShifted[cent])graphCombPi0InvYieldTot_yShifted[cent]->Write("graphInvYieldPi0CombXeXe5TeVTotErr_yShifted",TObject::kOverwrite);
    //             if(graphCombPi0InvYieldStat_yShifted[cent])graphCombPi0InvYieldStat_yShifted[cent]->Write("graphInvYieldPi0CombXeXe5TeVStatErr_yShifted",TObject::kOverwrite);
    //             if(graphCombPi0InvYieldSys_yShifted[cent])graphCombPi0InvYieldSys_yShifted[cent]->Write("graphInvYieldPi0CombXeXe5TeVSysErr_yShifted",TObject::kOverwrite);
    //         }
            // writing individual measurements
            for (Int_t meth = 0; meth< 11; meth++){
                if (graphIndPi0InvYieldStat[cent][meth]) graphIndPi0InvYieldStat[cent][meth]->Write(Form("graphInvYieldPi0%sStatErr", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                if (graphIndPi0InvYieldSys[cent][meth]) graphIndPi0InvYieldSys[cent][meth]->Write(Form("graphInvYieldPi0%sSysErr", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
                if (bWCorrection.Contains("Y")){
                    if (graphIndPi0InvYieldStat_yShifted[cent][meth]) graphIndPi0InvYieldStat_yShifted[cent][meth]->Write(Form("graphInvYieldPi0%sStatErr_yShifted", nameMeasGlobalLabel[meth].Data()),
                                                                                                                          TObject::kOverwrite);
                    if (graphIndPi0InvYieldSys_yShifted[cent][meth]) graphIndPi0InvYieldSys_yShifted[cent][meth]->Write(Form("graphInvYieldPi0%sSysErr_yShifted", nameMeasGlobalLabel[meth].Data()),
                                                                                                                        TObject::kOverwrite);
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
        }
    fCombResults.Close();


    //  **********************************************************************************************************************
    //  ************************* Saving only fits to final results **********************************************************
    //  **********************************************************************************************************************

//     TString nameOutputCommonFileFitsOnly    = Form("FitsPaperXeXe5TeV_%s.root", dateForOutput.Data());
//     TFile fFitsResults(nameOutputCommonFileFitsOnly.Data(), "RECREATE");
//
//         // fits for pi0
//         fitInvYieldPi0->Write("TsallisFitPi0");
//         fitTCMInvYieldPi0->Write("TwoComponentModelFitPi0");
//
//     fFitsResults.Close();

//     // **********************************************************************************************************************
//     // ************************* Saving comparison to comb for diff measurements ********************************************
//     // **********************************************************************************************************************
//
//     TString nameOutputCommonFileCompOnly    = Form("ComparisonsPaperXeXe5TeV_%s.root", dateForOutput.Data());
//     TFile fCompResults(nameOutputCommonFileCompOnly.Data(), "RECREATE");
//
//         graphRatioPi0CombCombFitStat->Write("Pi0_RatioCombToCombFit_Stat");
//         graphRatioPi0CombCombFitSys->Write("Pi0_RatioCombToCombFit_Syst");
//         for (Int_t meth = 0; meth < 11; meth++){
//             if (graphRatioPi0IndCombFitStat[meth]) graphRatioPi0IndCombFitStat[meth]->Write();
//             if (graphRatioPi0IndCombFitSys[meth]) graphRatioPi0IndCombFitSys[meth]->Write();
//         }
//     fCompResults.Close();
}

