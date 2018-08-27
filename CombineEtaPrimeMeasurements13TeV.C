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

void CombineEtaPrimeMeasurements13TeV(  TString fileNamePCM             = "",
                                        TString fileNameEMCAL           = "",
                                        TString fileNamePHOS            = "",
                                        TString fileNamePCMEMCAL        = "",
                                        TString fileNamePCMPHOS         = "",
                                        TString suffix                  = "eps",
                                        Bool_t isMC                     = kFALSE,
                                        TString bWCorrection            = "X",
                                        TString fileNameCorrFactors     = ""
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
    TString collisionSystem                     = "pp, #sqrt{#it{s}} = 13 TeV";
    TString outputDir                           = Form("%s/%s/CombineEtaPrimeMeasurements%s13TeV", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());
    TString outputDirFile                       = Form("%s/%s/CombineEtaPrimeMeasurements%s13TeV/Inputs", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());
    TString outputDirSupportComb                = Form("%s/%s/CombineEtaPrimeMeasurements%s13TeV/SupportingCombination", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());
    TString outputDirSupport                    = Form("%s/%s/CombineEtaPrimeMeasurements%s13TeV/Supporting", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());
    TString outputDirSupportPaper               = Form("%s/%s/CombineEtaPrimeMeasurements%s13TeV/SupportingPaper", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());

    TString nameFinalResDat                     = Form("%s/CombinedResults%s_FitResults.dat", outputDir.Data(), bWCorrection.Data());
    fstream  fileFitsOutput;
    fileFitsOutput.open(nameFinalResDat.Data(), ios::out);

    TString fileNameTheory                      = "ExternalInput/Theory/TheoryCompilation.root";

    gSystem->Exec("mkdir -p "+outputDirFile);
    gSystem->Exec("mkdir -p "+outputDirSupportComb);
    gSystem->Exec("mkdir -p "+outputDirSupport);
    gSystem->Exec("mkdir -p "+outputDirSupportPaper);
    if (fileNamePCM.CompareTo("")!=0 )          gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDirFile.Data()));
    if (fileNamePCMEMCAL.CompareTo("")!=0 )     gSystem->Exec(Form("cp %s %s/InputPCMEMCAL.root", fileNamePCMEMCAL.Data(), outputDirFile.Data()));
    if (fileNamePCMPHOS.CompareTo("")!=0 )      gSystem->Exec(Form("cp %s %s/InputPCMPHOS.root", fileNamePCMPHOS.Data(), outputDirFile.Data()));
    if (fileNamePHOS.CompareTo("")!=0 )         gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDirFile.Data()));
    if (fileNameEMCAL.CompareTo("")!=0 )        gSystem->Exec(Form("cp %s %s/InputEMCAL.root", fileNameEMCAL.Data(), outputDirFile.Data()));
    if (fileNameCorrFactors.CompareTo("")!=0 )  gSystem->Exec(Form("cp %s %s/InputCorrFactors.root", fileNameCorrFactors.Data(), outputDirFile.Data()));
    if (fileNameTheory.CompareTo("")!= 0)       gSystem->Exec(Form("cp %s %s/InputTheory.root", fileNameTheory.Data(), outputDirFile.Data()));

    TString prefix2                             = "";
    if (isMC){
        prefix2                                 = "MC";
    } else {
        prefix2                                 = "Data";
    }

    Double_t mesonMassExpectEtaPrime            = TDatabasePDG::Instance()->GetParticle(311)->Mass();

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

    TString  nameMeasGlobal[11]                 = { "PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal",
                                                    "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged",
                                                    "PCMOtherDataset"};
    TString  nameMeasGlobalLabel[11]            = { "PCM", "PHOS", "EMC", "PCM-PHOS", "PCM-EMC",
                                                    "PCM-Dal", "PHOS-Dal", "EMC-Dal", "EMChigh", "mEMC",
                                                    "PCMOtherDataset"};
    TString  nameSecEtaPrimeSourceRead[4]            = {"K0S", "K0L", "Lambda", "Rest"};
    TString  nameSecEtaPrimeSourceLabel[4]           = {"K^{0}_{s}", "K^{0}_{l}", "#Lambda", "had. int."};
    Double_t maxSecCorr[4]                      = { 0.05, 0.007, 0.0003, 0.03};
    Double_t branchingRatioEtaPrime[11]         = { 0.0218,  0.0218,  0.0218,  0.0218,  0.0218,
                                                    0.00001,  0.00001,  0.00001,  0.0218,  0.0218,
                                                     0.0218 };
    Double_t minPtEtaPrimePlotting              = 0.23;
    Double_t maxPtEtaPrimePlotting              = 51.;

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

    Color_t  color                              = GetColorDefaultColor("13TeV", "", "");
    Color_t  colorMC                            = GetColorDefaultColor("13TeV", "Pythia8", "");
    Style_t  markerStyle                        = GetDefaultMarkerStyle("13TeV", "", "");
    Style_t  markerStyleMC                      = GetDefaultMarkerStyle("13TeV", "Pythia8", "");
    Size_t   markerSize                         = GetDefaultMarkerSize("13TeV", "", "")*2;
    Size_t   markerSizeMC                       = GetDefaultMarkerSize("13TeV", "Pythia8", "");
    Double_t xSection13TeV                      = ReturnCorrectXSection("13TeV", 1);; // option is wrong fix when possible
    Double_t xSection13TeVErr                   = xSection13TeVINELErr*1e-3;

    //***********************************************************************************************
    //************************** Definition of final pt binning (has to be set manually) ************
    //***********************************************************************************************
    cout << "Setting EtaPrime binning" << endl;
    Double_t xPtLimitsEtaPrime[100];
    Int_t maxNBinsEtaPrimeAbs                            = 0;
    Int_t maxNBinsEtaPrime                               = 0;

    Int_t nTotMeasEtaPrime                               = 0;
    TString fileNamesMethod[11]                          = {"", "", "", "", "", "", "", "", "", "", ""};
    TString fileNamesEtaPrimeDetailedSys[11]             = {"", "", "", "", "", "", "", "", "", "", ""};
    Bool_t haveEtaPrimeSysDetailed[11]                   = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
    vector<TString>* ptSysRemNames                       = new vector<TString>[11];
    TFile* fileMethod[11];
    TDirectory* directoryEtaPrime[11];
    TGraphAsymmErrors* graphEtaPrimeEff[11];
    TGraphAsymmErrors* graphEtaPrimeAcc[11];
    TGraphAsymmErrors* graphEtaPrimeEffTimesAcc[11];
    TGraphAsymmErrors* graphEtaPrimeMass[11];
    TGraphAsymmErrors* graphEtaPrimeMassMC[11];
    TGraphAsymmErrors* graphEtaPrimeWidth[11];
    TGraphAsymmErrors* graphEtaPrimeWidthMC[11];
    TH1D* histoEtaPrimeInvYieldStat[11];
    TGraphAsymmErrors* graphEtaPrimeInvYieldStat[11];
    TGraphAsymmErrors* graphEtaPrimeInvYieldSys[11];

    maxNBinsEtaPrime                           = GetBinning( xPtLimitsEtaPrime, maxNBinsEtaPrimeAbs, "EtaPrime", "13TeV", 20, -1, kFALSE, "");
    cout << "EtaPrime "<< ": ";
    for (Int_t i = 0; i< maxNBinsEtaPrime+1; i++){
        cout << xPtLimitsEtaPrime[i] << ", " ;
    }
    cout << endl;
    for (Int_t meth = 0; meth < 11; meth++){
        fileMethod[meth]                      = NULL;
        directoryEtaPrime[meth]                    = NULL;
        graphEtaPrimeEff[meth]                     = NULL;
        graphEtaPrimeAcc[meth]                     = NULL;
        graphEtaPrimeEffTimesAcc[meth]             = NULL;
        graphEtaPrimeMass[meth]                    = NULL;
        graphEtaPrimeMassMC[meth]                  = NULL;
        graphEtaPrimeWidth[meth]                   = NULL;
        graphEtaPrimeWidthMC[meth]                 = NULL;
        histoEtaPrimeInvYieldStat[meth]            = NULL;
        graphEtaPrimeInvYieldStat[meth]            = NULL;
        graphEtaPrimeInvYieldSys[meth]             = NULL;
    }

    // *******************************************************************************************************
    // ********************* Set file names for inputs *******************************************************
    // *******************************************************************************************************
    if (fileNamePCM.CompareTo("")!=0)       fileNamesMethod[0]    = fileNamePCM;
    if (fileNamePHOS.CompareTo("")!=0)      fileNamesMethod[1]    = fileNamePHOS;
    if (fileNameEMCAL.CompareTo("")!=0)     fileNamesMethod[2]    = fileNameEMCAL;
    if (fileNamePCMPHOS.CompareTo("")!=0)   fileNamesMethod[3]    = fileNamePCMPHOS;
    if (fileNamePCMEMCAL.CompareTo("")!=0)  fileNamesMethod[4]    = fileNamePCMEMCAL;

    // *******************************************************************************************************
    // ***************************** Read in data for different methods **************************************
    // *******************************************************************************************************
    for (Int_t meth = 0; meth < 11; meth++){
        if (fileNamesMethod[meth].CompareTo("") != 0){
            cout << "Reading in " << nameMeasGlobalLabel[meth].Data() << " from " << fileNamesMethod[meth].Data() << endl;
            if (!fileMethod[meth] )
                fileMethod[meth]                               = new TFile(fileNamesMethod[meth].Data());
            if (!fileMethod[meth]) {
                cout << "file " << fileNamesMethod[meth].Data() << " not found! Skipping " << nameMeasGlobalLabel[meth].Data()  << endl;
                continue;
            }
            directoryEtaPrime[meth]                                = (TDirectory*)fileMethod[meth]->Get("EtaPrime13TeV");
            if (!directoryEtaPrime[meth]) {
                cout << "File doesn't contain directory " << "EtaPrime13TeV" << "! Skipping " << nameMeasGlobalLabel[meth].Data()  << endl;
                continue;
            } else {
                nTotMeasEtaPrime++;
                // reading supporting figures
                graphEtaPrimeMass[meth]                             = (TGraphAsymmErrors*)directoryEtaPrime[meth]->Get("EtaPrime_Mass_data");
                if (graphEtaPrimeMass[meth]) graphEtaPrimeMass[meth]                             = ScaleGraph(graphEtaPrimeMass[meth], 1000.);
                graphEtaPrimeWidth[meth]                            = (TGraphAsymmErrors*)directoryEtaPrime[meth]->Get("EtaPrime_Width_data");
                if (graphEtaPrimeWidth[meth]) graphEtaPrimeWidth[meth]                            = ScaleGraph(graphEtaPrimeWidth[meth], 1000.);
                graphEtaPrimeMassMC[meth]                           = (TGraphAsymmErrors*)directoryEtaPrime[meth]->Get("EtaPrime_Mass_MC");
                if (graphEtaPrimeMassMC[meth]) graphEtaPrimeMassMC[meth]                           = ScaleGraph(graphEtaPrimeMassMC[meth], 1000.);
                graphEtaPrimeWidthMC[meth]                          = (TGraphAsymmErrors*)directoryEtaPrime[meth]->Get("EtaPrime_Width_MC");
                if (graphEtaPrimeWidthMC[meth]) graphEtaPrimeWidthMC[meth]                          = ScaleGraph(graphEtaPrimeWidthMC[meth], 1000.);
                graphEtaPrimeAcc[meth]                              = (TGraphAsymmErrors*)directoryEtaPrime[meth]->Get("AcceptanceEtaPrime");
                graphEtaPrimeEff[meth]                              = (TGraphAsymmErrors*)directoryEtaPrime[meth]->Get("EfficiencyEtaPrime");
                graphEtaPrimeEffTimesAcc[meth]                      = (TGraphAsymmErrors*)directoryEtaPrime[meth]->Get("EffTimesAccEtaPrime");
                if (graphEtaPrimeEffTimesAcc[meth]) graphEtaPrimeEffTimesAcc[meth] = ScaleGraph( graphEtaPrimeEffTimesAcc[meth], branchingRatioEtaPrime[meth]);

                if (graphEtaPrimeWidth[meth]->GetY()[0] == 0 && graphEtaPrimeWidth[meth]->GetN()>0){
                    graphEtaPrimeWidth[meth]->RemovePoint(0);
                    graphEtaPrimeMass[meth]->RemovePoint(0);
                    graphEtaPrimeMassMC[meth]->RemovePoint(0);
                    graphEtaPrimeWidthMC[meth]->RemovePoint(0);
                    graphEtaPrimeEffTimesAcc[meth]->RemovePoint(0);
                }
                while (graphEtaPrimeWidth[meth]->GetY()[graphEtaPrimeWidth[meth]->GetN()-1] == 0 && graphEtaPrimeWidth[meth]->GetN()>0){
                    graphEtaPrimeWidth[meth]->RemovePoint(graphEtaPrimeWidth[meth]->GetN()-1);
                    graphEtaPrimeMass[meth]->RemovePoint(graphEtaPrimeMass[meth]->GetN()-1);
                    graphEtaPrimeMassMC[meth]->RemovePoint(graphEtaPrimeMassMC[meth]->GetN()-1);
                    graphEtaPrimeWidthMC[meth]->RemovePoint(graphEtaPrimeWidthMC[meth]->GetN()-1);
                    graphEtaPrimeEffTimesAcc[meth]->RemovePoint(graphEtaPrimeEffTimesAcc[meth]->GetN()-1);
                }

                // reading yields
                graphEtaPrimeInvYieldStat[meth]                     = (TGraphAsymmErrors*)directoryEtaPrime[meth]->Get("graphCorrectedYieldEtaPrime");
                histoEtaPrimeInvYieldStat[meth]                     = (TH1D*)directoryEtaPrime[meth]->Get("CorrectedYieldEtaPrime");
                cout << "EtaPrime stat error" << endl;
                graphEtaPrimeInvYieldStat[meth]->Print();
                graphEtaPrimeInvYieldSys[meth]                      = (TGraphAsymmErrors*)directoryEtaPrime[meth]->Get("EtaPrimeSystError");
                cout << "EtaPrime sys error" << endl;
                if (graphEtaPrimeInvYieldSys[meth]) graphEtaPrimeInvYieldSys[meth]->Print();

            }
        }
    }

/*
    // *******************************************************************************************************
    // ************************** Combination of different etaprime measurements **********************************
    // *******************************************************************************************************
    // REMARKS:
    //     - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - extendable to other energies
    //     - offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented

    TH1D* statErrorCollectionEtaPrime[11];
    TGraphAsymmErrors* statErrorGraphCollectionEtaPrime[11];
    TGraphAsymmErrors* sysErrorCollectionEtaPrime[11];
    for (Int_t meth = 0; meth< 11; meth++){
        // definition of array of histograms (NULL - means we have no measurement at this energy for this rec-method)
        // for statistical error and final value from respective method
        statErrorCollectionEtaPrime[meth]          = NULL;
        cout << meth << "\t" <<  histoEtaPrimeInvYieldStat[meth] << endl;
        if (histoEtaPrimeInvYieldStat[meth]) statErrorCollectionEtaPrime[meth]      = (TH1D*)histoEtaPrimeInvYieldStat[meth]->Clone(Form("statErr_%sEtaPrime",nameMeasGlobalLabel[meth].Data()));
        statErrorGraphCollectionEtaPrime[meth]     = NULL;
        if (graphEtaPrimeInvYieldStat[meth]) statErrorGraphCollectionEtaPrime[meth] = (TGraphAsymmErrors*)graphEtaPrimeInvYieldStat[meth]->Clone(Form("statErrGraph_%sEtaPrime",
                                                                                                                                                            nameMeasGlobalLabel[meth].Data()));
        // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
        // for systematic error from respective method
        sysErrorCollectionEtaPrime[meth]           = NULL;
        if (graphEtaPrimeInvYieldSys[meth]) sysErrorCollectionEtaPrime[meth]        = (TGraphAsymmErrors*)graphEtaPrimeInvYieldSys[meth]->Clone(Form("sysErr_%sEtaPrime",
                                                                                                                                                        nameMeasGlobalLabel[meth].Data()));
    }

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for EtaPrime
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsEtaPrime[11]             = { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 };
    Int_t offSetsEtaPrimeSys[11]          = { 1,  5,  8,  2,  5,  0,  0,  0,  0,  21, 0 };
    Int_t offSetEtaPrimeShifting[11]      = { 0,  4,  7,  1,  4,  0,  0,  0,  0,  21, 0 };
    Int_t nComBinsEtaPrimeShifting[11]    = { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 };
    Double_t minPtEtaPrime                = 0.4;

    TGraphAsymmErrors* statErrorRelCollectionEtaPrime[11];
    TGraphAsymmErrors* sysErrorRelCollectionEtaPrime[11];
    TGraph* graphWeightsEtaPrime[11];
    for (int_t meth = 0; meth< 11; meth++){
        graphweightsetaprime[meth]                 = null;
        staterrorrelcollectionetaprime[meth]       = null;
        syserrorrelcollectionetaprime[meth]        = null;
    }
    for (Int_t meth = 0; meth< 11; meth++){
        if (statErrorCollectionEtaPrime[meth]){
            statErrorRelCollectionEtaPrime[meth]   = new TGraphAsymmErrors(statErrorCollectionEtaPrime[meth]);
            RemoveZerosAtBeginningAndEndFromGraph(statErrorRelCollectionEtaPrime[meth]);
            statErrorRelCollectionEtaPrime[meth]   = CalculateRelErrUpAsymmGraph( statErrorRelCollectionEtaPrime[meth], Form("relativeStatErrorEtaPrime_%s", nameMeasGlobal[meth].Data()));
        }
        if (sysErrorCollectionEtaPrime[meth]){
            sysErrorRelCollectionEtaPrime[meth]    = (TGraphAsymmErrors*)sysErrorCollectionEtaPrime[meth]->Clone(Form("relativeSysErrorEtaPrime_%s", nameMeasGlobal[meth].Data()));
            RemoveZerosAtBeginningAndEndFromGraph(sysErrorRelCollectionEtaPrime[meth]);
            sysErrorRelCollectionEtaPrime[meth]    = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionEtaPrime[meth], Form("relativeSysErrorEtaPrime_%s", nameMeasGlobal[meth].Data()));
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Assuming maximal correlation *********************************************
    // **********************************************************************************************************************

    TGraphAsymmErrors* graphCombEtaPrimeInvYieldStat             = NULL;
    TGraphAsymmErrors* graphCombEtaPrimeInvYieldStatWOXErr       = NULL;
    TGraphAsymmErrors* graphCombEtaPrimeInvYieldSys              = NULL;
    TGraphAsymmErrors* graphCombEtaPrimeInvYieldTot              = NULL;
    TGraphAsymmErrors* graphCombEtaPrimeInvYieldTotUnshi         = NULL;
    TGraphAsymmErrors* graphCombEtaPrimeInvYieldStatUnshi        = NULL;
    TGraphAsymmErrors* graphCombEtaPrimeInvYieldSysUnshi         = NULL;
    TGraphAsymmErrors* graphCombEtaPrimeInvYieldRelStat          = NULL;
    TGraphAsymmErrors* graphCombEtaPrimeInvYieldRelSys           = NULL;
    TGraphAsymmErrors* graphCombEtaPrimeInvYieldRelTot           = NULL;
    TGraphAsymmErrors* graphCombEtaPrimeInvYieldStat_yShifted    = NULL;
    TGraphAsymmErrors* graphCombEtaPrimeInvYieldSys_yShifted     = NULL;
    TGraphAsymmErrors* graphCombEtaPrimeInvYieldTot_yShifted     = NULL;
    TGraphAsymmErrors* graphRatioEtaPrimeCombCombFitTot          = NULL;
    TGraphAsymmErrors* graphRatioEtaPrimeCombCombFitStat         = NULL;
    TGraphAsymmErrors* graphRatioEtaPrimeCombCombFitSys          = NULL;
    TGraphAsymmErrors* graphRatioEtaPrimeCombCombFitStatWOXErr   = NULL;

    TGraphAsymmErrors* graphIndEtaPrimeInvYieldStatUnshi[11];
    TGraphAsymmErrors* graphIndEtaPrimeInvYieldSysUnshi[11];
    TGraphAsymmErrors* graphIndEtaPrimeInvYieldStat[11];
    TGraphAsymmErrors* graphIndEtaPrimeInvYieldSys[11];
    TGraphAsymmErrors* graphIndEtaPrimeInvYieldStat_yShifted[11];
    TGraphAsymmErrors* graphIndEtaPrimeInvYieldSys_yShifted[11];
    TGraphAsymmErrors* graphRatioEtaPrimeIndCombFitStat[11];
    TGraphAsymmErrors* graphRatioEtaPrimeIndCombFitSys[11];
    for (Int_t meth = 0; meth< 11; meth++){
        graphIndEtaPrimeInvYieldStatUnshi[meth]                = NULL;
        graphIndEtaPrimeInvYieldSysUnshi[meth]                 = NULL;
        graphIndEtaPrimeInvYieldStat[meth]                     = NULL;
        graphIndEtaPrimeInvYieldSys[meth]                      = NULL;
        graphIndEtaPrimeInvYieldStat_yShifted[meth]            = NULL;
        graphIndEtaPrimeInvYieldSys_yShifted[meth]             = NULL;
        graphRatioEtaPrimeIndCombFitStat[meth]                 = NULL;
        graphRatioEtaPrimeIndCombFitSys[meth]                  = NULL;
    }
    Double_t paramTCMEtaPrimeNew[5]                          = { 20.72,0.167, 0.454, 0.629, 3.163};
    Double_t paramGraphEtaPrime[3]                           = {1.0e12, 8., 0.13};
    TF1* fitTCMDecomposedLEtaPrime                       = NULL;
    TF1* fitTCMDecomposedHEtaPrime                       = NULL;
    TF1* fitTCMInvYieldEtaPrime                          = NULL;
    TF1* fitInvYieldEtaPrime                             = NULL;
    TF1* fitInvYieldEtaPrimeGraph                        = NULL;
    TF1* fitPowInvYieldEtaPrime                          = NULL;

    // *************************************************************************************************************
    // ************************************** Define plotting environment ******************************************
    *************************************************************************************************************
    textSizeLabelsPixel                 = 900*0.04;

    TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
    canvasWeights->SetLogx();
    TH2F * histo2DEtaPrimeWeights;
    histo2DEtaPrimeWeights = new TH2F("histo2DEtaPrimeWeights","histo2DEtaPrimeWeights",11000,minPtEtaPrimePlotting, maxPtEtaPrimePlotting,1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DEtaPrimeWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DEtaPrimeWeights->GetXaxis()->SetMoreLogLabels();
    histo2DEtaPrimeWeights->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelWeightsEtaPrime         = new TLatex(0.95,0.15,"#eta' #rightarrow #gamma#gamma");
    SetStyleTLatex( labelWeightsEtaPrime, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelWeightsEnergy      = new TLatex(0.95,0.20,collisionSystem.Data());
    SetStyleTLatex( labelWeightsEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);


    TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErr->SetLogx();
    TH2F * histo2DRelSysErr;
    histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,minPtEtaPrimePlotting, maxPtEtaPrimePlotting,1000,0,50.0);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErr->GetYaxis()->SetRangeUser(0,39.5);
    histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErr->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelRelSysErrEnergy    = new TLatex(0.15,0.89,collisionSystem.Data());
    SetStyleTLatex( labelRelSysErrEnergy, textSizeLabelsPixel, 4, 1, 43);
    TLatex *labelRelSysErrEtaPrime       = new TLatex(0.15,0.85,"#eta' #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrEtaPrime, textSizeLabelsPixel, 4, 1, 43);

    TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelStatErr->SetLogx();
    TH2F * histo2DRelStatErr;
    histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,minPtEtaPrimePlotting, maxPtEtaPrimePlotting,1000,0,50.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErr->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelRelStatErrEnergy   = new TLatex(0.95,0.89,collisionSystem.Data());
    SetStyleTLatex( labelRelStatErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelStatErrEtaPrime      = new TLatex(0.95,0.85,"#eta' #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelStatErrEtaPrime, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    TCanvas* canvasRelTotErr            = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelTotErr->SetLogx();
    TH2F * histo2DRelTotErrEtaPrime;
    histo2DRelTotErrEtaPrime                 = new TH2F("histo2DRelTotErrEtaPrime","histo2DRelTotErrEtaPrime",11000,minPtEtaPrimePlotting, maxPtEtaPrimePlotting,1000,0,50.0);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrEtaPrime, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrEtaPrime->GetYaxis()->SetRangeUser(0,39.5);
    histo2DRelTotErrEtaPrime->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrEtaPrime->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelRelTotErrEnergy    = new TLatex(0.95,0.89,collisionSystem.Data());
    SetStyleTLatex( labelRelTotErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelTotErrEtaPrime       = new TLatex(0.95,0.85,"#eta' #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrEtaPrime, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    // Declaration & calculation of combined spectrum
    TString fileNameEtaPrimeOutputWeighting      = Form("%s/EtaPrime_WeightingMethod.dat",outputDirSupportComb.Data());
    graphCombEtaPrimeInvYieldTot                 = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionEtaPrime,    sysErrorCollectionEtaPrime,
                                                                                    xPtLimitsEtaPrime, maxNBinsEtaPrime,
                                                                                    offSetsEtaPrime, offSetsEtaPrimeSys,
                                                                                    graphCombEtaPrimeInvYieldStat, graphCombEtaPrimeInvYieldSys,
                                                                                    fileNameEtaPrimeOutputWeighting, "13TeV", "EtaPrime", kTRUE,
                                                                                    NULL, fileNameCorrFactors, ""
                                                                                );


    if (graphCombEtaPrimeInvYieldTot == NULL) {
        cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
        return;
    }
    while (graphCombEtaPrimeInvYieldStat->GetX()[0] < minPtEtaPrime){
        graphCombEtaPrimeInvYieldStat->RemovePoint(0);
    }
    while (graphCombEtaPrimeInvYieldTot->GetX()[0] < minPtEtaPrime){
        graphCombEtaPrimeInvYieldTot->RemovePoint(0);
    }
    while (graphCombEtaPrimeInvYieldSys->GetX()[0] < minPtEtaPrime){
        graphCombEtaPrimeInvYieldSys->RemovePoint(0);
    }
//         graphCombEtaPrimeInvYieldTot->Print();

    // Reading weights from output file for plotting
    ifstream fileWeightsEtaPrimeRead;
    fileWeightsEtaPrimeRead.open(fileNameEtaPrimeOutputWeighting,ios_base::in);
    cout << "reading" << fileNameEtaPrimeOutputWeighting << endl;
    Double_t xValuesEtaPrimeRead[100];
    Double_t weightsEtaPrimeRead[11][100];
    Int_t availableEtaPrimeMeas[11]    = {   -1, -1, -1, -1, -1,
                                        -1, -1, -1, -1, -1,
                                        -1};
    Int_t nMeasSetEtaPrime             = nTotMeasEtaPrime;
    Int_t nPtBinsEtaPrimeRead          = 0;
    while(!fileWeightsEtaPrimeRead.eof() && nPtBinsEtaPrimeRead < 100){
        TString garbage             = "";
        if (nPtBinsEtaPrimeRead == 0){
            fileWeightsEtaPrimeRead >> garbage ;//>> availableEtaPrimeMeas[0] >> availableEtaPrimeMeas[1] >> availableEtaPrimeMeas[2] >> availableEtaPrimeMeas[3];
            for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
                fileWeightsEtaPrimeRead >> availableEtaPrimeMeas[i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
                cout << availableEtaPrimeMeas[i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsEtaPrimeRead >> xValuesEtaPrimeRead[nPtBinsEtaPrimeRead-1];
            for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
                fileWeightsEtaPrimeRead >> weightsEtaPrimeRead[availableEtaPrimeMeas[i]][nPtBinsEtaPrimeRead-1] ;
            }
            cout << "read: "<<  nPtBinsEtaPrimeRead << "\t"<< xValuesEtaPrimeRead[nPtBinsEtaPrimeRead-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
                cout << weightsEtaPrimeRead[availableEtaPrimeMeas[i]][nPtBinsEtaPrimeRead-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsEtaPrimeRead++;
    }
    nPtBinsEtaPrimeRead                  = nPtBinsEtaPrimeRead-2 ;
    fileWeightsEtaPrimeRead.close();

    for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
        graphWeightsEtaPrime[availableEtaPrimeMeas[i]]                        = new TGraph(nPtBinsEtaPrimeRead,xValuesEtaPrimeRead,weightsEtaPrimeRead[availableEtaPrimeMeas[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsEtaPrimeRead; n++){
            if (graphWeightsEtaPrime[availableEtaPrimeMeas[i]]->GetY()[bin] == 0) graphWeightsEtaPrime[availableEtaPrimeMeas[i]]->RemovePoint(bin);
            else bin++;
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Plotting weights method only EMC *****************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 900*0.04;
    canvasWeights->cd();
    histo2DEtaPrimeWeights->Draw("copy");

        TLegend* legendWeights   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEtaPrime), 32);
        for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
            DrawGammaSetMarkerTGraph(graphWeightsEtaPrime[availableEtaPrimeMeas[i]], markerStyleDet[availableEtaPrimeMeas[i]], markerSizeDet[availableEtaPrimeMeas[i]]*0.5, colorDet[availableEtaPrimeMeas[i]] , colorDet[availableEtaPrimeMeas[i]]);
            graphWeightsEtaPrime[availableEtaPrimeMeas[i]]->Draw("p,same,z");
            legendWeights->AddEntry(graphWeightsEtaPrime[availableEtaPrimeMeas[i]],nameMeasGlobalLabel[availableEtaPrimeMeas[i]],"p");
        }
        legendWeights->Draw();

        labelWeightsEnergy->SetText(0.95,0.20,collisionSystem.Data());
        labelWeightsEnergy->Draw();
        labelWeightsEtaPrime->Draw();

        DrawGammaLines(minPtEtaPrimePlotting, maxPtEtaPrimePlotting , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(minPtEtaPrimePlotting, maxPtEtaPrimePlotting , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(minPtEtaPrimePlotting, maxPtEtaPrimePlotting , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(minPtEtaPrimePlotting, maxPtEtaPrimePlotting , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/EtaPrime_Weights.%s",outputDirSupportComb.Data(), suffix.Data()));

    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    canvasRelSysErr->cd();
    histo2DRelSysErr->GetYaxis()->SetRangeUser(0,39.5);
    histo2DRelSysErr->Draw("copy");

        cout << "sys error etaprime" << endl;
        TLegend* legendRelSysErr        = GetAndSetLegend2(0.62, 0.92-(0.04*nMeasSetEtaPrime), 0.95, 0.92, textSizeLabelsPixel);
        for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollectionEtaPrime[availableEtaPrimeMeas[i]], markerStyleDet[availableEtaPrimeMeas[i]], markerSizeDet[availableEtaPrimeMeas[i]]*0.5, colorDet[availableEtaPrimeMeas[i]],
                                        colorDet[availableEtaPrimeMeas[i]]);
            sysErrorRelCollectionEtaPrime[availableEtaPrimeMeas[i]]->Draw("p,same,z");
            legendRelSysErr->AddEntry(sysErrorRelCollectionEtaPrime[availableEtaPrimeMeas[i]],nameMeasGlobalLabel[availableEtaPrimeMeas[i]],"p");
        }
        legendRelSysErr->Draw();


        labelRelSysErrEnergy->SetText(0.15,0.89,collisionSystem.Data());
        labelRelSysErrEnergy->Draw();
        labelRelSysErrEtaPrime->Draw();

    canvasRelSysErr->SaveAs(Form("%s/EtaPrime_RelSysErr.%s",outputDirSupportComb.Data(), suffix.Data()));
    delete legendRelSysErr;
    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************

    canvasRelStatErr->cd();
    histo2DRelStatErr->GetYaxis()->SetRangeUser(0,39.5);
    histo2DRelStatErr->Draw("copy");
        TLegend* legendRelStatErr       = GetAndSetLegend2(0.14, 0.92-(0.04*nMeasSetEtaPrime), 0.45, 0.92, textSizeLabelsPixel);
        for (Int_t i = 0; i < nMeasSetEtaPrime; i++){
            DrawGammaSetMarkerTGraph(statErrorRelCollectionEtaPrime[availableEtaPrimeMeas[i]], markerStyleDet[availableEtaPrimeMeas[i]], markerSizeDet[availableEtaPrimeMeas[i]]*0.5, colorDet[availableEtaPrimeMeas[i]],
                                        colorDet[availableEtaPrimeMeas[i]]);
            statErrorRelCollectionEtaPrime[availableEtaPrimeMeas[i]]->Draw("p,same,z");
            legendRelStatErr->AddEntry(statErrorRelCollectionEtaPrime[availableEtaPrimeMeas[i]],nameMeasGlobalLabel[availableEtaPrimeMeas[i]],"p");
        }
        legendRelStatErr->Draw();

        labelRelStatErrEnergy->SetText(0.95,0.89,collisionSystem.Data());
        labelRelStatErrEnergy->Draw();
        labelRelStatErrEtaPrime->Draw();

    canvasRelStatErr->SaveAs(Form("%s/EtaPrime_RelStatErr.%s",outputDirSupportComb.Data(), suffix.Data()));

    //  *********************************************************************************************************************
    //  ************************ Visualize relative total errors of different combination methods EtaPrime ***********************
    //  *********************************************************************************************************************

    graphCombEtaPrimeInvYieldRelStat   = CalculateRelErrUpAsymmGraph( graphCombEtaPrimeInvYieldStat, "relativeStatErrorEtaPrime_Method");
    graphCombEtaPrimeInvYieldRelSys    = CalculateRelErrUpAsymmGraph( graphCombEtaPrimeInvYieldSys, "relativeSysErrorEtaPrime_Method");
    graphCombEtaPrimeInvYieldRelTot    = CalculateRelErrUpAsymmGraph( graphCombEtaPrimeInvYieldTot, "relativeTotalErrorEtaPrime_Method");

    canvasRelTotErr->cd();
    histo2DRelTotErrEtaPrime->GetYaxis()->SetRangeUser(0,39.5);
    histo2DRelTotErrEtaPrime->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombEtaPrimeInvYieldRelTot, markerStyleComb, markerSizeComb, kBlue+2 , kBlue+2);
        graphCombEtaPrimeInvYieldRelTot->Draw("p,same,z");

        TLegend* legendRelTotErr1       = GetAndSetLegend2(0.20, 0.92-(0.035), 0.45, 0.92, 32);
        legendRelTotErr1->AddEntry(graphCombEtaPrimeInvYieldRelTot,"All","p");
        legendRelTotErr1->Draw();


        labelRelTotErrEnergy->SetText(0.95,0.89,collisionSystem.Data());
        labelRelTotErrEnergy->Draw();
        labelRelTotErrEtaPrime->Draw();

    canvasRelTotErr->SaveAs(Form("%s/EtaPrime_TotErr_Comp.%s",outputDirSupportComb.Data(), suffix.Data()));
    histo2DRelTotErrEtaPrime->GetYaxis()->SetTitle("Err (%)");
    histo2DRelTotErrEtaPrime->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombEtaPrimeInvYieldRelTot, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombEtaPrimeInvYieldRelTot->Draw("p,same,z");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaPrimeInvYieldRelStat, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombEtaPrimeInvYieldRelStat->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaPrimeInvYieldRelSys, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombEtaPrimeInvYieldRelSys->SetLineStyle(7);
        graphCombEtaPrimeInvYieldRelSys->Draw("l,x0,same,e1");

        TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.45, 0.92, 32);
        legendRelTotErr3->AddEntry(graphCombEtaPrimeInvYieldRelTot,"tot","p");
        legendRelTotErr3->AddEntry(graphCombEtaPrimeInvYieldRelStat,"stat","l");
        legendRelTotErr3->AddEntry(graphCombEtaPrimeInvYieldRelSys,"sys","l");
        legendRelTotErr3->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrEtaPrime->Draw();

    canvasRelTotErr->SaveAs(Form("%s/EtaPrime_Reldecomp.%s",outputDirSupportComb.Data(), suffix.Data()));

    // **********************************************************************************************************************
    // ************************************* Calculating bin shifted spectra & fitting **************************************
    // **********************************************************************************************************************

    // Cloning spectra
    graphCombEtaPrimeInvYieldTotUnshi          = (TGraphAsymmErrors*)graphCombEtaPrimeInvYieldTot->Clone("EtaPrimeUnshifted");
    graphCombEtaPrimeInvYieldStatUnshi         = (TGraphAsymmErrors*)graphCombEtaPrimeInvYieldStat->Clone("EtaPrimeUnshiftedStat");
    graphCombEtaPrimeInvYieldSysUnshi          = (TGraphAsymmErrors*)graphCombEtaPrimeInvYieldSys->Clone("EtaPrimeUnshiftedSys");

    for (Int_t meth = 0; meth< 11; meth++){
        if (statErrorGraphCollectionEtaPrime[meth]){
            graphIndEtaPrimeInvYieldStatUnshi[meth]                 = (TGraphAsymmErrors*)statErrorGraphCollectionEtaPrime[meth]->Clone(Form("EtaPrimeUnshiftedStat%s",nameMeasGlobalLabel[meth].Data()));
            graphIndEtaPrimeInvYieldStat[meth]                      = (TGraphAsymmErrors*)statErrorGraphCollectionEtaPrime[meth]->Clone(Form("EtaPrimeStat%s",nameMeasGlobalLabel[meth].Data()));
            graphIndEtaPrimeInvYieldStat_yShifted[meth]             = (TGraphAsymmErrors*)statErrorGraphCollectionEtaPrime[meth]->Clone(Form("EtaPrimeYShiftedStat%s",nameMeasGlobalLabel[meth].Data()));
        }
        if (sysErrorCollectionEtaPrime[meth]){
            graphIndEtaPrimeInvYieldSysUnshi[meth]                  = (TGraphAsymmErrors*)sysErrorCollectionEtaPrime[meth]->Clone(Form("EtaPrimeUnshiftedSys%s",nameMeasGlobalLabel[meth].Data()));
            graphIndEtaPrimeInvYieldSys[meth]                       = (TGraphAsymmErrors*)sysErrorCollectionEtaPrime[meth]->Clone(Form("EtaPrimeSys%s",nameMeasGlobalLabel[meth].Data()));
            graphIndEtaPrimeInvYieldSys_yShifted[meth]              = (TGraphAsymmErrors*)sysErrorCollectionEtaPrime[meth]->Clone(Form("EtaPrimeYShiftedSys%s",nameMeasGlobalLabel[meth].Data()));
        }
    }

    // fitting spectrum with intial parameters
    // Two component model fit from Bylinkin
    fitTCMDecomposedLEtaPrime           = FitObject("tcmlow","twoCompModelEtaPrime_DecL", "EtaPrime", NULL, minPtEtaPrime, 2.);
    fitTCMDecomposedHEtaPrime           = FitObject("tcmhigh","twoCompModelEtaPrime_DecH", "EtaPrime", NULL, 4, 50.);
    fitTCMDecomposedLEtaPrime->SetParameters(graphCombEtaPrimeInvYieldTot->GetY()[2],minPtEtaPrime);
    graphCombEtaPrimeInvYieldStat->Fit(fitTCMDecomposedLEtaPrime,"QNRMEX0+","",minPtEtaPrime,0.8);
    graphCombEtaPrimeInvYieldStat->Fit(fitTCMDecomposedHEtaPrime,"QNRMEX0+","",3, xPtLimitsEtaPrime[maxNBinsEtaPrime]);

    cout << WriteParameterToFile(fitTCMDecomposedLEtaPrime)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedLEtaPrime)<< endl;
    cout << WriteParameterToFile(fitTCMDecomposedHEtaPrime)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedHEtaPrime)<< endl;

    fitTCMInvYieldEtaPrime              = FitObject("tcm","fitTCMInvYieldEtaPrime","EtaPrime",graphCombEtaPrimeInvYieldStat,
                                                        minPtEtaPrime, xPtLimitsEtaPrime[maxNBinsEtaPrime] , paramTCMEtaPrimeNew,
                                                        "QNRMEX0+", "", kFALSE);
    fitTCMDecomposedLEtaPrime->SetParameter(0, fitTCMInvYieldEtaPrime->GetParameter(0));
    fitTCMDecomposedLEtaPrime->SetParameter(1, fitTCMInvYieldEtaPrime->GetParameter(1));
    fitTCMDecomposedHEtaPrime->SetParameter(0, fitTCMInvYieldEtaPrime->GetParameter(2));
    fitTCMDecomposedHEtaPrime->SetParameter(1, fitTCMInvYieldEtaPrime->GetParameter(3));
    fitTCMDecomposedHEtaPrime->SetParameter(2, fitTCMInvYieldEtaPrime->GetParameter(4));

    // Tsallis fit
    fitInvYieldEtaPrime                 = FitObject("l","fitInvYieldEtaPrime","EtaPrime",histoEtaPrimeInvYieldStat[2], minPtEtaPrime, xPtLimitsEtaPrime[maxNBinsEtaPrime], paramGraphEtaPrime,"QNRMEX0+");
    fitInvYieldEtaPrimeGraph            = (TF1*)fitInvYieldEtaPrime->Clone("fitInvYieldEtaPrimeGraph");

    // *************************************************************************************************************
    // Shift graphs in X direction if desired
    // *************************************************************************************************************
    if(bWCorrection.Contains("X")){
        TF1* fitShiftingEtaPrime            = FitObject("tmpt","ShiftingEtaPrime","EtaPrime");
        fitShiftingEtaPrime->SetParameters(fitInvYieldEtaPrime->GetParameter(0),fitInvYieldEtaPrime->GetParameter(1), fitInvYieldEtaPrime->GetParameter(2));
//         TF1* fitShiftingEtaPrime                 = FitObject("tcmpt","ShiftingEtaPrime","EtaPrime");
//         fitShiftingEtaPrime->SetParameters(fitTCMInvYieldEtaPrime->GetParameter(0),fitTCMInvYieldEtaPrime->GetParameter(1), fitTCMInvYieldEtaPrime->GetParameter(2), fitTCMInvYieldEtaPrime->GetParameter(3),fitTCMInvYieldEtaPrime->GetParameter(4));

        TGraphAsymmErrors* graphCombEtaPrimeInvYieldTotNoShift = (TGraphAsymmErrors*) graphCombEtaPrimeInvYieldTot->Clone("EtaPrime_NoShift");

        graphCombEtaPrimeInvYieldTot          = ApplyXshift(graphCombEtaPrimeInvYieldTot, fitShiftingEtaPrime);
//             cout << "comb" << endl;
//             graphCombEtaPrimeInvYieldStat->Print();
        graphCombEtaPrimeInvYieldStat         = ApplyXshiftIndividualSpectra( graphCombEtaPrimeInvYieldTot,
                                                                            graphCombEtaPrimeInvYieldStat,
                                                                            fitShiftingEtaPrime,
                                                                            0, graphCombEtaPrimeInvYieldStat->GetN());
        graphCombEtaPrimeInvYieldSys          = ApplyXshiftIndividualSpectra( graphCombEtaPrimeInvYieldTot,
                                                                            graphCombEtaPrimeInvYieldSys,
                                                                            fitShiftingEtaPrime,
                                                                            0, graphCombEtaPrimeInvYieldSys->GetN());
        for (Int_t meth = 0; meth< 11; meth++){
            if (graphIndEtaPrimeInvYieldStat[meth]){
                cout << "shiting stat err of " << nameMeasGlobalLabel[meth].Data();
                graphIndEtaPrimeInvYieldStat[meth]  = ApplyXshiftIndividualSpectra( graphCombEtaPrimeInvYieldTot,
                                                                            graphIndEtaPrimeInvYieldStat[meth],
                                                                            fitShiftingEtaPrime,
                                                                            offSetEtaPrimeShifting[meth], nComBinsEtaPrimeShifting[meth]);

            }
            if (graphIndEtaPrimeInvYieldSys[meth]){
                cout << "shiting sys err of " << nameMeasGlobalLabel[meth].Data();
                graphIndEtaPrimeInvYieldSys[meth]   = ApplyXshiftIndividualSpectra( graphCombEtaPrimeInvYieldTot,
                                                                            graphIndEtaPrimeInvYieldSys[meth],
                                                                            fitShiftingEtaPrime,
                                                                            offSetEtaPrimeShifting[meth], nComBinsEtaPrimeShifting[meth]);
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

            Int_t numberPoints      = graphCombEtaPrimeInvYieldTotNoShift->GetN();
            Double_t *xPoint        = graphCombEtaPrimeInvYieldTotNoShift->GetX();
            Double_t* xvalueErrUp   = graphCombEtaPrimeInvYieldTotNoShift->GetEXhigh();
            Double_t* xvalueErrLow  = graphCombEtaPrimeInvYieldTotNoShift->GetEXlow();
            Double_t *xPointShift= graphCombEtaPrimeInvYieldTot->GetX();
            for (Int_t i=0; i<numberPoints; i++) {
                graphCombEtaPrimeInvYieldTotNoShift->SetPoint(i,xPoint[i],xPointShift[i]/xPoint[i]);
                graphCombEtaPrimeInvYieldTotNoShift->SetPointError(i,xvalueErrLow[i],xvalueErrUp[i],0,0);
            }
            DrawGammaSetMarkerTGraphAsym(graphCombEtaPrimeInvYieldTotNoShift, markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphCombEtaPrimeInvYieldTotNoShift->Draw("p same");

            TLatex *labelRatioToFitBinShift   = new TLatex(0.94, 0.91, collisionSystem.Data());
            SetStyleTLatex( labelRatioToFitBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
            labelRatioToFitBinShift->Draw();
            TLatex *labelRatioToFitALICEBinShift    = new TLatex(0.94, 0.86, textALICE.Data());
            SetStyleTLatex( labelRatioToFitALICEBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
            labelRatioToFitALICEBinShift->Draw();
            TLatex *labelRatioToFitEtaPrimeBinShift      = new TLatex(0.94, 0.807, "#eta' #rightarrow #gamma#gamma");
            SetStyleTLatex( labelRatioToFitEtaPrimeBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
            labelRatioToFitEtaPrimeBinShift->Draw();

        canvasShift->Update();
        canvasShift->SaveAs(Form("%s/BinShiftCorrection_EtaPrime.%s",outputDirSupportComb.Data(), suffix.Data()));

        // *************************************************************************************************************
        // Plot control graphs
        // *************************************************************************************************************

        TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
        DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.09);
        canvasDummy2->SetLogy();
        canvasDummy2->SetLogx();
        TH2F * histo2DDummy2;
        histo2DDummy2               = new TH2F("histo2DDummy2","histo2DDummy2",1000,minPtEtaPrimePlotting,maxPtEtaPrimePlotting,1000,1e-9,10e1);
        SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
        histo2DDummy2->DrawCopy();

            for (Int_t meth = 0; meth< 11; meth++){
                if (graphIndEtaPrimeInvYieldSys[meth]){
                    DrawGammaSetMarkerTGraphAsym(graphIndEtaPrimeInvYieldSys[meth], markerStyleDet[meth] ,markerSizeDet[meth]/2, colorDet[meth], colorDet[meth], widthLinesBoxes, kTRUE);
                    graphIndEtaPrimeInvYieldSys[meth]->Draw("pEsame");
                }
            }

            DrawGammaSetMarkerTGraphAsym(graphCombEtaPrimeInvYieldStatUnshi, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombEtaPrimeInvYieldStatUnshi->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaPrimeInvYieldStat, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombEtaPrimeInvYieldStat->Draw("pEsame");

            fitTCMInvYieldEtaPrime->SetLineColor(kBlue+2);
            fitTCMInvYieldEtaPrime->Draw("same");

        canvasDummy2->Update();
        canvasDummy2->Print(Form("%s/ComparisonShiftedEtaPrime.%s",outputDirSupportComb.Data(), suffix.Data()));
        delete canvasDummy2;
    }

    // *************************************************************************************************************
    // redo fitting after binshifts
    // *************************************************************************************************************
    // Tsallis function
    cout << WriteParameterToFile(fitInvYieldEtaPrime)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitInvYieldEtaPrime)<< endl;
    //Two component model from Bylinkin
    fitTCMInvYieldEtaPrime        = FitObject("tcm","fitTCMInvYieldEtaPrime13TeVCent","EtaPrime",graphCombEtaPrimeInvYieldTot, minPtEtaPrime, xPtLimitsEtaPrime[maxNBinsEtaPrime],
                                                paramTCMEtaPrimeNew,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitTCMInvYieldEtaPrime)<< endl;
    fileFitsOutput <<  WriteParameterToFile(fitTCMInvYieldEtaPrime)<< endl;

    fitPowInvYieldEtaPrime        = FitObject("m","fitPowInvYieldEtaPrime13TeVCent","EtaPrime",graphCombEtaPrimeInvYieldTot,5,40. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvYieldEtaPrime)<< endl;

    // *************************************************************************************************************
    // Shift graphs in Y direction as well if desired
    // *************************************************************************************************************

    if(bWCorrection.Contains("Y") ){
        graphCombEtaPrimeInvYieldTot_yShifted        = (TGraphAsymmErrors*)graphCombEtaPrimeInvYieldTotUnshi->Clone("EtaPrimeYShiftedCombTot");
        graphCombEtaPrimeInvYieldTot_yShifted        =  ApplyYshiftIndividualSpectra( graphCombEtaPrimeInvYieldTot_yShifted, fitInvYieldEtaPrime);
        graphCombEtaPrimeInvYieldStat_yShifted       = (TGraphAsymmErrors*)graphCombEtaPrimeInvYieldStatUnshi->Clone("EtaPrimeYShiftedCombStat");
        graphCombEtaPrimeInvYieldStat_yShifted       =  ApplyYshiftIndividualSpectra( graphCombEtaPrimeInvYieldStat_yShifted, fitInvYieldEtaPrime);
        graphCombEtaPrimeInvYieldSys_yShifted        = (TGraphAsymmErrors*)graphCombEtaPrimeInvYieldSysUnshi->Clone("EtaPrimeYShiftedCombSys");
        graphCombEtaPrimeInvYieldSys_yShifted        =  ApplyYshiftIndividualSpectra( graphCombEtaPrimeInvYieldSys_yShifted, fitInvYieldEtaPrime);

        for (Int_t meth = 0; meth< 11; meth++){
            if (graphIndEtaPrimeInvYieldStat_yShifted[meth]){
                graphIndEtaPrimeInvYieldStat_yShifted[meth] = ApplyYshiftIndividualSpectra( graphIndEtaPrimeInvYieldStat_yShifted[meth], fitInvYieldEtaPrime);
            }
            if (graphIndEtaPrimeInvYieldSys_yShifted[meth]){
                graphIndEtaPrimeInvYieldSys_yShifted[meth]  = ApplyYshiftIndividualSpectra( graphIndEtaPrimeInvYieldSys_yShifted[meth], fitInvYieldEtaPrime);
            }
        }
    }


    graphRatioEtaPrimeCombCombFitTot                       = (TGraphAsymmErrors*)graphCombEtaPrimeInvYieldTot->Clone();
    graphRatioEtaPrimeCombCombFitTot                       = CalculateGraphErrRatioToFit(graphRatioEtaPrimeCombCombFitTot, fitTCMInvYieldEtaPrime);
    graphRatioEtaPrimeCombCombFitStat                      = (TGraphAsymmErrors*)graphCombEtaPrimeInvYieldStat->Clone();
    graphRatioEtaPrimeCombCombFitStat                      = CalculateGraphErrRatioToFit(graphRatioEtaPrimeCombCombFitStat, fitTCMInvYieldEtaPrime);
    graphRatioEtaPrimeCombCombFitSys                       = (TGraphAsymmErrors*)graphCombEtaPrimeInvYieldSys->Clone();
    graphRatioEtaPrimeCombCombFitSys                       = CalculateGraphErrRatioToFit(graphRatioEtaPrimeCombCombFitSys, fitTCMInvYieldEtaPrime);
    graphCombEtaPrimeInvYieldStatWOXErr                    = (TGraphAsymmErrors*)graphCombEtaPrimeInvYieldStat->Clone("graphCombEtaPrimeInvYieldStatWOXErr");
    ProduceGraphAsymmWithoutXErrors(graphCombEtaPrimeInvYieldStatWOXErr);
    graphRatioEtaPrimeCombCombFitStatWOXErr                = (TGraphAsymmErrors*)graphRatioEtaPrimeCombCombFitStat->Clone("graphRatioEtaPrimeCombCombFitStatWOXErr");
    ProduceGraphAsymmWithoutXErrors(graphRatioEtaPrimeCombCombFitStatWOXErr);

    for (Int_t meth = 0; meth< 11; meth++){
        if (graphIndEtaPrimeInvYieldStat[meth]){
            graphRatioEtaPrimeIndCombFitStat[meth]              = (TGraphAsymmErrors*)graphIndEtaPrimeInvYieldStat[meth]->Clone(Form("RatioEtaPrime%sToCombFitStat", nameMeasGlobalLabel[meth].Data()));
            graphRatioEtaPrimeIndCombFitStat[meth]              = CalculateGraphErrRatioToFit(graphRatioEtaPrimeIndCombFitStat[meth], fitTCMInvYieldEtaPrime);
        }
        if (graphIndEtaPrimeInvYieldSys[meth]){
            graphRatioEtaPrimeIndCombFitSys[meth]              = (TGraphAsymmErrors*)graphIndEtaPrimeInvYieldSys[meth]->Clone(Form("RatioEtaPrime%sToCombFitSyst", nameMeasGlobalLabel[meth].Data()));
            graphRatioEtaPrimeIndCombFitSys[meth]              = CalculateGraphErrRatioToFit(graphRatioEtaPrimeIndCombFitSys[meth], fitTCMInvYieldEtaPrime);
        }
    }

    // **********************************************************************************************************************
    // ******************************************* Plot Ratio of Comb to Fit ************************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 54;
    TCanvas* canvasRatioToCombFit       = new TCanvas("canvasRatioToCombFit","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToCombFit, 0.08, 0.01, 0.01, 0.125);
    canvasRatioToCombFit->SetLogx();

        Double_t textsizeLabels      = 0;
        if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
            textsizeLabels           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
        } else {
            textsizeLabels           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
        }
        cout << textsizeLabels << endl;

    TH2F * histo2DEtaPrimeRatioToCombFit;
    histo2DEtaPrimeRatioToCombFit               = new TH2F("histo2DEtaPrimeRatioToCombFit","histo2DEtaPrimeRatioToCombFit",1000,minPtEtaPrimePlotting, maxPtEtaPrimePlotting,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DEtaPrimeRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabels, textsizeLabels,
                                0.85*textsizeLabels,textsizeLabels, 0.9, 0.65, 510, 505);
    histo2DEtaPrimeRatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DEtaPrimeRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
    histo2DEtaPrimeRatioToCombFit->GetYaxis()->SetRangeUser(0.2,1.82);
    histo2DEtaPrimeRatioToCombFit->Draw("copy");

        ProduceGraphAsymmWithoutXErrors(graphRatioEtaPrimeCombCombFitStat);

        DrawGammaSetMarkerTGraphAsym(graphRatioEtaPrimeCombCombFitSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioEtaPrimeCombCombFitSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaPrimeCombCombFitStat, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioEtaPrimeCombCombFitStat->Draw("p,same,z");

        DrawGammaLines(minPtEtaPrimePlotting, maxPtEtaPrimePlotting , 1., 1.,0.1, kGray+2);
        DrawGammaLines(minPtEtaPrimePlotting, maxPtEtaPrimePlotting , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(minPtEtaPrimePlotting, maxPtEtaPrimePlotting , 0.9, 0.9,0.1, kGray, 7);

        TLatex *labelRatioToFitEnergy   = new TLatex(0.95, 0.92, collisionSystem.Data());
        SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitEnergy->Draw();
        TLatex *labelRatioToFitALICE    = new TLatex(0.95, 0.86, textALICE.Data());
        SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitALICE->Draw();
        TLatex *labelRatioToFitEtaPrime      = new TLatex(0.12, 0.92, "#eta' #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitEtaPrime, textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
        labelRatioToFitEtaPrime->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/EtaPrime_RatioOfCombToCombFit.%s",outputDirSupportPaper.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // *******************************************Plot Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DEtaPrimeRatioToCombFit->Draw("copy");

        for (Int_t meth = 10; meth > -1 ; meth--){
            if (graphRatioEtaPrimeIndCombFitSys[meth]){
                DrawGammaSetMarkerTGraphAsym(graphRatioEtaPrimeIndCombFitSys[meth], markerStyleDet[meth] ,markerSizeDet[meth]*0.5, colorDet[meth], colorDet[meth], widthLinesBoxes, kTRUE);
                graphRatioEtaPrimeIndCombFitSys[meth]->Draw("E2same");
            }
            if (graphRatioEtaPrimeIndCombFitStat[meth]){
                ProduceGraphAsymmWithoutXErrors(graphRatioEtaPrimeIndCombFitStat[meth]);
                DrawGammaSetMarkerTGraphAsym(graphRatioEtaPrimeIndCombFitStat[meth], markerStyleDet[meth] ,markerSizeDet[meth]*0.5, colorDet[meth], colorDet[meth]);
                graphRatioEtaPrimeIndCombFitStat[meth]->Draw("p,same,z");
            }
        }
        if (graphRatioEtaPrimeIndCombFitStat[4]) graphRatioEtaPrimeIndCombFitStat[4]->Draw("p,same,z");

        DrawGammaLines(minPtEtaPrimePlotting, maxPtEtaPrimePlotting , 1., 1.,0.5, kGray+2);
        DrawGammaLines(minPtEtaPrimePlotting, maxPtEtaPrimePlotting , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(minPtEtaPrimePlotting, maxPtEtaPrimePlotting , 0.9, 0.9,0.5, kGray, 7);
        DrawGammaLines(minPtEtaPrimePlotting, maxPtEtaPrimePlotting , 1.2, 1.2,0.5, kGray, 9);
        DrawGammaLines(minPtEtaPrimePlotting, maxPtEtaPrimePlotting , 0.8, 0.8,0.5, kGray, 9);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitEtaPrime->Draw();
        histo2DEtaPrimeRatioToCombFit->Draw("same,axis");

        //****************************** Definition of the Legend ******************************************
        //**************** Row def ************************
        Double_t rowsLegendOnlyEtaPrimeRatio[4]          = {0.30, 0.25, 0.20, 0.15};
        Double_t rowsLegendOnlyEtaPrimeRatioAbs[4]       = {0.92, 0.47, 0.38, 0.29 };
        Double_t columnsLegendOnlyEtaPrimeRatio[6]       = {0.14, 0.26, 0.35, 0.48, 0.7, 0.8};
        Double_t columnsLegendOnlyEtaPrimeRatioAbs[6]    = {0.215, 0.8, 1.17, 2, 11.2, 16.8};
        Double_t columnsLegendOnlyEtaPrimeRatioAbs2[6]   = {0.215, 0.8, 1.68, 2, 11.2, 24.4};
        Double_t lengthBox                          = 0.2;
        Double_t heightBox                          = 0.06/2;
        //****************** first Column **************************************************

        TLatex *textPCMOnlyRatioEtaPrime                 = new TLatex(columnsLegendOnlyEtaPrimeRatio[0],rowsLegendOnlyEtaPrimeRatio[1],nameMeasGlobalLabel[0]);
        SetStyleTLatex( textPCMOnlyRatioEtaPrime, textSizeLabelsPixel,4, 1, 43);

        if (graphRatioEtaPrimeIndCombFitSys[0]) textPCMOnlyRatioEtaPrime->Draw();
        TLatex *textPHOSOnlyRatioEtaPrime                = new TLatex(columnsLegendOnlyEtaPrimeRatio[0],rowsLegendOnlyEtaPrimeRatio[2],nameMeasGlobalLabel[1]);
        SetStyleTLatex( textPHOSOnlyRatioEtaPrime, textSizeLabelsPixel,4, 1, 43);
        if (graphRatioEtaPrimeIndCombFitSys[1]) textPHOSOnlyRatioEtaPrime->Draw();
        TLatex *textEMCALOnlyRatioEtaPrime               = new TLatex(columnsLegendOnlyEtaPrimeRatio[0],rowsLegendOnlyEtaPrimeRatio[3],nameMeasGlobalLabel[2]);
        SetStyleTLatex( textEMCALOnlyRatioEtaPrime, textSizeLabelsPixel,4, 1, 43);
        if (graphRatioEtaPrimeIndCombFitSys[2]) textEMCALOnlyRatioEtaPrime->Draw();
        TLatex *textPCMEMCALOnlyRatioEtaPrime            = new TLatex(columnsLegendOnlyEtaPrimeRatio[3],rowsLegendOnlyEtaPrimeRatio[1],nameMeasGlobalLabel[4]);
        SetStyleTLatex( textPCMEMCALOnlyRatioEtaPrime, textSizeLabelsPixel,4, 1, 43);
        if (graphRatioEtaPrimeIndCombFitSys[4]) textPCMEMCALOnlyRatioEtaPrime->Draw();
        TLatex *textPCMPHOSOnlyRatioEtaPrime            = new TLatex(columnsLegendOnlyEtaPrimeRatio[3],rowsLegendOnlyEtaPrimeRatio[2],nameMeasGlobalLabel[3]);
        SetStyleTLatex( textPCMPHOSOnlyRatioEtaPrime, textSizeLabelsPixel,4, 1, 43);
        if (graphRatioEtaPrimeIndCombFitSys[3]) textPCMPHOSOnlyRatioEtaPrime->Draw();
        TLatex *textDalitzOnlyRatioEtaPrime            = new TLatex(columnsLegendOnlyEtaPrimeRatio[3],rowsLegendOnlyEtaPrimeRatio[3],nameMeasGlobalLabel[5]);
        SetStyleTLatex( textDalitzOnlyRatioEtaPrime, textSizeLabelsPixel,4, 1, 43);
        if (graphRatioEtaPrimeIndCombFitSys[5] )textDalitzOnlyRatioEtaPrime->Draw();

        //****************** second Column *************************************************
        TLatex *textStatOnlyRatioEtaPrime                = new TLatex(columnsLegendOnlyEtaPrimeRatio[1],rowsLegendOnlyEtaPrimeRatio[0] ,"stat");
        SetStyleTLatex( textStatOnlyRatioEtaPrime, textSizeLabelsPixel,4, 1, 43);
        textStatOnlyRatioEtaPrime->Draw();
        TLatex *textSysOnlyRatioEtaPrime                 = new TLatex(columnsLegendOnlyEtaPrimeRatio[2] ,rowsLegendOnlyEtaPrimeRatio[0],"syst");
        SetStyleTLatex( textSysOnlyRatioEtaPrime, textSizeLabelsPixel,4, 1, 43);
        textSysOnlyRatioEtaPrime->Draw();
        TLatex *textStatOnlyRatioEtaPrime2               = new TLatex(columnsLegendOnlyEtaPrimeRatio[4],rowsLegendOnlyEtaPrimeRatio[0] ,"stat");
        SetStyleTLatex( textStatOnlyRatioEtaPrime2, textSizeLabelsPixel,4, 1, 43);
        textStatOnlyRatioEtaPrime2->Draw();
        TLatex *textSysOnlyRatioEtaPrime2                = new TLatex(columnsLegendOnlyEtaPrimeRatio[5] ,rowsLegendOnlyEtaPrimeRatio[0],"syst");
        SetStyleTLatex( textSysOnlyRatioEtaPrime2, textSizeLabelsPixel,4, 1, 43);
        textSysOnlyRatioEtaPrime2->Draw();

        if (graphRatioEtaPrimeIndCombFitSys[0]){
            TMarker* markerPCMEtaPrimeOnlyRatioEtaPrime           = CreateMarkerFromGraph(graphRatioEtaPrimeIndCombFitSys[0],columnsLegendOnlyEtaPrimeRatio[1] ,rowsLegendOnlyEtaPrimeRatio[1],1);
            markerPCMEtaPrimeOnlyRatioEtaPrime->DrawMarker(columnsLegendOnlyEtaPrimeRatioAbs[1] ,rowsLegendOnlyEtaPrimeRatioAbs[1]);
        }
        if (graphRatioEtaPrimeIndCombFitSys[1]){
            TMarker* markerPHOSEtaPrimeOnlyRatioEtaPrime          = CreateMarkerFromGraph(graphRatioEtaPrimeIndCombFitSys[1], columnsLegendOnlyEtaPrimeRatio[1] ,rowsLegendOnlyEtaPrimeRatio[2],1);
            markerPHOSEtaPrimeOnlyRatioEtaPrime->DrawMarker(columnsLegendOnlyEtaPrimeRatioAbs[1] ,rowsLegendOnlyEtaPrimeRatioAbs[2]);
        }
        if (graphRatioEtaPrimeIndCombFitSys[2]){
            TMarker* markerEMCALEtaPrimeOnlyRatioEtaPrime         = CreateMarkerFromGraph(graphRatioEtaPrimeIndCombFitSys[2], columnsLegendOnlyEtaPrimeRatio[1] ,rowsLegendOnlyEtaPrimeRatio[3],1);
            markerEMCALEtaPrimeOnlyRatioEtaPrime->DrawMarker(columnsLegendOnlyEtaPrimeRatioAbs[1] ,rowsLegendOnlyEtaPrimeRatioAbs[3]);
        }
        if (graphRatioEtaPrimeIndCombFitSys[4]){
            TMarker* markerPCMEMCALEtaPrimeOnlyRatioEtaPrime      = CreateMarkerFromGraph(graphRatioEtaPrimeIndCombFitSys[4], columnsLegendOnlyEtaPrimeRatio[3] ,rowsLegendOnlyEtaPrimeRatio[1],1);
            markerPCMEMCALEtaPrimeOnlyRatioEtaPrime->DrawMarker(columnsLegendOnlyEtaPrimeRatioAbs[4] ,rowsLegendOnlyEtaPrimeRatioAbs[1]);
        }
        if (graphRatioEtaPrimeIndCombFitSys[3]){
            TMarker* markerPCMPHOSEtaPrimeOnlyRatioEtaPrime      = CreateMarkerFromGraph(graphRatioEtaPrimeIndCombFitSys[3], columnsLegendOnlyEtaPrimeRatio[3] ,rowsLegendOnlyEtaPrimeRatio[2],1);
            markerPCMPHOSEtaPrimeOnlyRatioEtaPrime->DrawMarker(columnsLegendOnlyEtaPrimeRatioAbs[4] ,rowsLegendOnlyEtaPrimeRatioAbs[2]);
        }
        if (graphRatioEtaPrimeIndCombFitSys[5]){
            TMarker* markerDalitzEtaPrimeOnlyRatioEtaPrime      = CreateMarkerFromGraph(graphRatioEtaPrimeIndCombFitSys[5], columnsLegendOnlyEtaPrimeRatio[3] ,rowsLegendOnlyEtaPrimeRatio[3],1);
            markerDalitzEtaPrimeOnlyRatioEtaPrime->DrawMarker(columnsLegendOnlyEtaPrimeRatioAbs[4] ,rowsLegendOnlyEtaPrimeRatioAbs[3]);
        }

        if (graphRatioEtaPrimeIndCombFitSys[0]){
            TBox* boxPCMEtaPrimeOnlyRatioEtaPrime                 = CreateBoxFromGraph(graphRatioEtaPrimeIndCombFitSys[0], columnsLegendOnlyEtaPrimeRatioAbs[2] , rowsLegendOnlyEtaPrimeRatioAbs[1]- heightBox,
                                                                            columnsLegendOnlyEtaPrimeRatioAbs2[2], rowsLegendOnlyEtaPrimeRatioAbs[1]+ heightBox);
            boxPCMEtaPrimeOnlyRatioEtaPrime->Draw("l");
        }
        if (graphRatioEtaPrimeIndCombFitSys[1]){
            TBox* boxPHOSEtaPrimeOnlyRatioEtaPrime                = CreateBoxFromGraph(graphRatioEtaPrimeIndCombFitSys[1], columnsLegendOnlyEtaPrimeRatioAbs[2] , rowsLegendOnlyEtaPrimeRatioAbs[2]- heightBox,
                                                                            columnsLegendOnlyEtaPrimeRatioAbs2[2], rowsLegendOnlyEtaPrimeRatioAbs[2]+ heightBox);
            boxPHOSEtaPrimeOnlyRatioEtaPrime->Draw("l");
        }
        if (graphRatioEtaPrimeIndCombFitSys[2]){
            TBox* boxEMCALEtaPrimeOnlyRatioEtaPrime               = CreateBoxFromGraph(graphRatioEtaPrimeIndCombFitSys[2], columnsLegendOnlyEtaPrimeRatioAbs[2] , rowsLegendOnlyEtaPrimeRatioAbs[3]- heightBox,
                                                                            columnsLegendOnlyEtaPrimeRatioAbs2[2], rowsLegendOnlyEtaPrimeRatioAbs[3]+ heightBox);
            boxEMCALEtaPrimeOnlyRatioEtaPrime->Draw("l");
        }
        if (graphRatioEtaPrimeIndCombFitSys[4]){
            TBox* boxPCMEMCALEtaPrimeOnlyRatioEtaPrime            = CreateBoxFromGraph(graphRatioEtaPrimeIndCombFitSys[4], columnsLegendOnlyEtaPrimeRatioAbs[5] , rowsLegendOnlyEtaPrimeRatioAbs[1]- heightBox,
                                                                            columnsLegendOnlyEtaPrimeRatioAbs2[5], rowsLegendOnlyEtaPrimeRatioAbs[1]+ heightBox);
            boxPCMEMCALEtaPrimeOnlyRatioEtaPrime->Draw("l");
        }
        if (graphRatioEtaPrimeIndCombFitSys[3]){
            TBox* boxPCMPHOSEtaPrimeOnlyRatioEtaPrime             = CreateBoxFromGraph(graphRatioEtaPrimeIndCombFitSys[3], columnsLegendOnlyEtaPrimeRatioAbs[5], rowsLegendOnlyEtaPrimeRatioAbs[2]- heightBox,
                                                                            columnsLegendOnlyEtaPrimeRatioAbs2[5], rowsLegendOnlyEtaPrimeRatioAbs[2]+ heightBox);
            boxPCMPHOSEtaPrimeOnlyRatioEtaPrime->Draw("l");
        }
        if (graphRatioEtaPrimeIndCombFitSys[5]){
            TBox* boxDalitzEtaPrimeOnlyRatioEtaPrime             = CreateBoxFromGraph(graphRatioEtaPrimeIndCombFitSys[5], columnsLegendOnlyEtaPrimeRatioAbs[5], rowsLegendOnlyEtaPrimeRatioAbs[3]- heightBox,
                                                                                columnsLegendOnlyEtaPrimeRatioAbs2[5], rowsLegendOnlyEtaPrimeRatioAbs[3]+ heightBox);
            boxDalitzEtaPrimeOnlyRatioEtaPrime->Draw("l");
        }
        canvasRatioToCombFit->SaveAs(Form("%s/EtaPrime_RatioOfIndividualMeasToCombFit.%s",outputDirSupportPaper.Data(), suffix.Data()));*/

    // **********************************************************************************************************************
    // ******************************************* Mass and width for etaprime ***************************************************
    // **********************************************************************************************************************

    Double_t arrayBoundariesX1_4[2];
    Double_t arrayBoundariesY1_4[3];
    Double_t relativeMarginsX[3];
    Double_t relativeMarginsY[3];
    Int_t textSizeLabelsPixel           = 50;
    Double_t textSizeLabelsRel          = ((Double_t)textSizeLabelsPixel)/1250;
    ReturnCorrectValuesForCanvasScaling(1350,1250, 1, 2,0.1, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

    TCanvas* canvasMassWidthEtaPrime         = new TCanvas("canvasMassWidthEtaPrime","",0,0,1350,1250);  // gives the page size
    DrawGammaCanvasSettings( canvasMassWidthEtaPrime,  0.13, 0.02, 0.03, 0.06);

    TPad* padWidthEtaPrime                   = new TPad("padWidthEtaPrime", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
    DrawGammaPadSettings( padWidthEtaPrime, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
    padWidthEtaPrime->Draw();

    TPad* padMassEtaPrime                    = new TPad("padMassEtaPrime", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
    DrawGammaPadSettings( padMassEtaPrime, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
    padMassEtaPrime->Draw();

    TPad* padMassLegendEtaPrime                 = new TPad("padMassLegendEtaPrime", "", 0.13, 0.36, 0.82, 0.52,-1, -1, -2);
    DrawGammaPadSettings( padMassLegendEtaPrime, 0., 0., 0., 0.);
    padMassLegendEtaPrime->SetFillStyle(0);
    padMassLegendEtaPrime->Draw();

    Double_t margin                 = relativeMarginsX[0]*2*1350;
    Double_t textsizeLabelsWidth    = 0;
    Double_t textsizeFacWidth       = 0;
    if (padWidthEtaPrime->XtoPixel(padWidthEtaPrime->GetX2()) < padWidthEtaPrime->YtoPixel(padWidthEtaPrime->GetY1())){
        textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthEtaPrime->XtoPixel(padWidthEtaPrime->GetX2()) ;
        textsizeFacWidth            = (Double_t)1./padWidthEtaPrime->XtoPixel(padWidthEtaPrime->GetX2()) ;
    } else {
        textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthEtaPrime->YtoPixel(padWidthEtaPrime->GetY1());
        textsizeFacWidth            = (Double_t)1./padWidthEtaPrime->YtoPixel(padWidthEtaPrime->GetY1());
    }
//     cout << textsizeLabelsWidth << endl;

    TH2F * histo2DAllEtaPrimeFWHM        = new TH2F("histo2DAllEtaPrimeFWHM","histo2DAllEtaPrimeFWHM", 20, minPtEtaPrimePlotting, maxPtEtaPrimePlotting ,1000., -30, 40);
    SetStyleHistoTH2ForGraphs(histo2DAllEtaPrimeFWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                              0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
    histo2DAllEtaPrimeFWHM->GetYaxis()->SetRangeUser(-1.,30.5);
    histo2DAllEtaPrimeFWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllEtaPrimeFWHM->GetYaxis()->SetNdivisions(505);
    histo2DAllEtaPrimeFWHM->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllEtaPrimeFWHM->GetXaxis()->SetTickLength(0.05);
    histo2DAllEtaPrimeFWHM->GetYaxis()->SetTickLength(0.026);

    Double_t textsizeLabelsMass     = 0;
    Double_t textsizeFacMass        = 0;
    if (padMassEtaPrime->XtoPixel(padMassEtaPrime->GetX2()) <padMassEtaPrime->YtoPixel(padMassEtaPrime->GetY1()) ){
        textsizeLabelsMass          = (Double_t)textSizeLabelsPixel/padMassEtaPrime->XtoPixel(padMassEtaPrime->GetX2()) ;
        textsizeFacMass             = (Double_t)1./padMassEtaPrime->XtoPixel(padMassEtaPrime->GetX2()) ;
    } else {
        textsizeLabelsMass          = (Double_t)textSizeLabelsPixel/padMassEtaPrime->YtoPixel(padMassEtaPrime->GetY1());
        textsizeFacMass             = (Double_t)1./padMassEtaPrime->YtoPixel(padMassEtaPrime->GetY1());
    }

    TH2F * histo2DAllEtaPrimeMass        = new TH2F("histo2DAllEtaPrimeMass","histo2DAllEtaPrimeMass",20, minPtEtaPrimePlotting, maxPtEtaPrimePlotting, 1000., 750., 1350);
    SetStyleHistoTH2ForGraphs(histo2DAllEtaPrimeMass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                              textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
    histo2DAllEtaPrimeMass->GetYaxis()->SetRangeUser(125.5,162.8);
    histo2DAllEtaPrimeMass->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllEtaPrimeMass->GetYaxis()->SetNdivisions(505);
    histo2DAllEtaPrimeMass->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllEtaPrimeMass->GetXaxis()->SetTickLength(0.05);
    histo2DAllEtaPrimeMass->GetXaxis()->SetLabelOffset(-0.015);

    TLatex *labelLegendAMass        = new TLatex(0.13,0.06,"a)");
    SetStyleTLatex( labelLegendAMass, textSizeLabelsPixel,4, 1, 43);
    TLatex *labelMassPerf           = new TLatex(0.13,0.875,"ALICE performance");
    SetStyleTLatex( labelMassPerf, textSizeLabelsPixel,4, 1, 43);
    TLatex *labelMassEnergy         = new TLatex(0.13,0.775,collisionSystem.Data());
    SetStyleTLatex( labelMassEnergy, textSizeLabelsPixel,4, 1, 43);
    TLatex *labelMassEtaPrime            = new TLatex(0.13,0.69,"#eta' #rightarrow #gamma#gamma");
    SetStyleTLatex( labelMassEtaPrime, textSizeLabelsPixel,4, 1, 43);
    TLatex *labelLegendBMass        = new TLatex(0.13,0.22,"b)");
    SetStyleTLatex( labelLegendBMass, textSizeLabelsPixel,4, 1, 43);
    //********************************** Defintion of the Legend **************************************************
    Double_t columnsLegendMass[6]  = {0.,0.25,0.38, 0.5, 0.78, 0.9};
    Double_t rowsLegendMass[5]     = {0.8,0.6,0.4,0.2,0.01};
    Double_t offsetMarkerYMass     = 0.06;
    Double_t offsetMarkerXMass     = 0.05;
    Double_t scaleMarkerMass       = 1.2;

    padWidthEtaPrime->cd();
    padWidthEtaPrime->SetLogx();

        histo2DAllEtaPrimeFWHM->DrawCopy();

        for (Int_t meth = 10; meth > -1; meth--){
            if (graphEtaPrimeWidth[meth]){
                DrawGammaSetMarkerTGraphAsym(graphEtaPrimeWidth[meth], markerStyleDet[meth], markerSizeDet[meth]*0.55, colorDet[meth] , colorDet[meth]);
                graphEtaPrimeWidth[meth]->Draw("p,same,z");
            }
            if (graphEtaPrimeWidthMC[meth]){
                DrawGammaSetMarkerTGraphAsym(graphEtaPrimeWidthMC[meth], markerStyleDetMC[meth], markerSizeDetMC[meth]*0.55, colorDetMC[meth] , colorDetMC[meth]);
                graphEtaPrimeWidthMC[meth]->Draw("p,same,z");
            }
        }

        labelLegendAMass->Draw();
        labelMassPerf->Draw();
        labelMassEnergy->SetText(0.13,0.775,collisionSystem.Data());
        labelMassEnergy->Draw();
        labelMassEtaPrime->Draw();

    padMassEtaPrime->cd();
    padMassEtaPrime->SetLogx();

        histo2DAllEtaPrimeMass->DrawCopy();

        for (Int_t meth = 10; meth > -1; meth--){
            if (graphEtaPrimeMass[meth]){
                DrawGammaSetMarkerTGraphAsym(graphEtaPrimeMass[meth], markerStyleDet[meth], markerSizeDet[meth]*0.55, colorDet[meth] , colorDet[meth]);
                graphEtaPrimeMass[meth]->Draw("p,same,z");
            }
            if (graphEtaPrimeMassMC[meth]){
                DrawGammaSetMarkerTGraphAsym(graphEtaPrimeMassMC[meth], markerStyleDetMC[meth], markerSizeDetMC[meth]*0.55, colorDetMC[meth] , colorDetMC[meth]);
                graphEtaPrimeMassMC[meth]->Draw("p,same,z");
            }
        }

        DrawGammaLines(minPtEtaPrimePlotting, maxPtEtaPrimePlotting , mesonMassExpectEtaPrime*1000., mesonMassExpectEtaPrime*1000.,0.1, kGray);
        labelLegendBMass->Draw();

        padMassLegendEtaPrime->cd();
        padMassLegendEtaPrime->Clear();
        //****************** first Column **************************************************
        TLatex *textMassPCM2            = new TLatex(columnsLegendMass[0],rowsLegendMass[1],nameMeasGlobalLabel[0]);
        SetStyleTLatex( textMassPCM2, textSizeLabelsPixel,4, 1, 43);
        if (graphEtaPrimeMass[0]) textMassPCM2->Draw();
        TLatex *textMassPCMEMCAL       = new TLatex(columnsLegendMass[0],rowsLegendMass[2],nameMeasGlobalLabel[4]);
        SetStyleTLatex( textMassPCMEMCAL, textSizeLabelsPixel,4, 1, 43);
        if (graphEtaPrimeMass[4]) textMassPCMEMCAL->Draw();
        TLatex *textMassEMCAL          = new TLatex(columnsLegendMass[0],rowsLegendMass[3],nameMeasGlobalLabel[2]);
        SetStyleTLatex( textMassEMCAL, textSizeLabelsPixel,4, 1, 43);
        if (graphEtaPrimeMass[2]) textMassEMCAL->Draw();
        TLatex *textMassPHOS           = new TLatex(columnsLegendMass[3],rowsLegendMass[1],nameMeasGlobalLabel[1]);
        SetStyleTLatex( textMassPHOS, textSizeLabelsPixel,4, 1, 43);
        if (graphEtaPrimeMass[1]) textMassPHOS->Draw();
        TLatex *textMassPCMPHOS             = new TLatex(columnsLegendMass[3],rowsLegendMass[2],nameMeasGlobalLabel[3]);
        SetStyleTLatex( textMassPCMPHOS, textSizeLabelsPixel,4, 1, 43);
        if (graphEtaPrimeMass[3]) textMassPCMPHOS->Draw();
        TLatex *textMassDalitz             = new TLatex(columnsLegendMass[3],rowsLegendMass[2],nameMeasGlobalLabel[5]);
        SetStyleTLatex( textMassDalitz, textSizeLabelsPixel,4, 1, 43);
        if (graphEtaPrimeMass[5]) textMassDalitz->Draw();

        //****************** second Column *************************************************
        TLatex *textMassData1           = new TLatex(columnsLegendMass[1],rowsLegendMass[0] ,"Data");
        SetStyleTLatex( textMassData1, textSizeLabelsPixel,4, 1, 43);
        if (graphEtaPrimeMass[2] || graphEtaPrimeMass[4] || graphEtaPrimeMass[0]) textMassData1->Draw();
        TLatex *textMassMC1             = new TLatex(columnsLegendMass[2] ,rowsLegendMass[0],"MC");
        SetStyleTLatex( textMassMC1, textSizeLabelsPixel,4, 1, 43);
        if (graphEtaPrimeMass[2] || graphEtaPrimeMass[4] || graphEtaPrimeMass[0]) textMassMC1->Draw();

        TLatex *textMassData2           = new TLatex(columnsLegendMass[4],rowsLegendMass[0] ,"Data");
        SetStyleTLatex( textMassData2, textSizeLabelsPixel,4, 1, 43);
        if (graphEtaPrimeMass[1] || graphEtaPrimeMass[3] ) textMassData2->Draw();
        TLatex *textMassMC2             = new TLatex(columnsLegendMass[5] ,rowsLegendMass[0],"MC");
        SetStyleTLatex( textMassMC2, textSizeLabelsPixel,4, 1, 43);
        if (graphEtaPrimeMass[1] || graphEtaPrimeMass[3] ) textMassMC2->Draw();

        TMarker* markerEtaPrimeMass[11]      = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
        TMarker* markerEtaPrimeMassMC[11]    = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

        for (Int_t meth = 0; meth< 11; meth++){
            if (graphEtaPrimeMass[meth]) markerEtaPrimeMass[meth]            = CreateMarkerFromGraph(graphEtaPrimeMass[meth], 0.1, 0.1, scaleMarkerMass);
            if (graphEtaPrimeMassMC[meth]) markerEtaPrimeMassMC[meth]        = CreateMarkerFromGraph(graphEtaPrimeMassMC[meth], 0.1, 0.1, scaleMarkerMass);
        }
        if (markerEtaPrimeMass[0] ) markerEtaPrimeMass[0]->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[1]+ offsetMarkerYMass);
        if (markerEtaPrimeMass[4] ) markerEtaPrimeMass[4]->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
        if (markerEtaPrimeMass[2] ) markerEtaPrimeMass[2]->DrawMarker(columnsLegendMass[1]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);
        if (markerEtaPrimeMass[1] ) markerEtaPrimeMass[1]->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[1]+ offsetMarkerYMass);
        if (markerEtaPrimeMass[3] ) markerEtaPrimeMass[3]->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[2]+ offsetMarkerYMass);
        if (markerEtaPrimeMass[5] ) markerEtaPrimeMass[5]->DrawMarker(columnsLegendMass[4]+ offsetMarkerXMass ,rowsLegendMass[3]+ offsetMarkerYMass);

        if (markerEtaPrimeMassMC[0]) markerEtaPrimeMassMC[0]->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass-0.01 ,rowsLegendMass[1]+ offsetMarkerYMass);
        if (markerEtaPrimeMassMC[4]) markerEtaPrimeMassMC[4]->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass-0.01 ,rowsLegendMass[2]+ offsetMarkerYMass);
        if (markerEtaPrimeMassMC[2]) markerEtaPrimeMassMC[2]->DrawMarker(columnsLegendMass[2]+ offsetMarkerXMass-0.01 ,rowsLegendMass[3]+ offsetMarkerYMass);
        if (markerEtaPrimeMassMC[1]) markerEtaPrimeMassMC[1]->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass-0.01 ,rowsLegendMass[1]+ offsetMarkerYMass);
        if (markerEtaPrimeMassMC[3]) markerEtaPrimeMassMC[3]->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass-0.01 ,rowsLegendMass[2]+ offsetMarkerYMass);
        if (markerEtaPrimeMassMC[5]) markerEtaPrimeMassMC[5]->DrawMarker(columnsLegendMass[5]+ offsetMarkerXMass-0.01 ,rowsLegendMass[3]+ offsetMarkerYMass);

    canvasMassWidthEtaPrime->Update();
    canvasMassWidthEtaPrime->Print(Form("%s/EtaPrime_MassAndWidth.%s",outputDirSupportPaper.Data(), suffix.Data()));

    for (Int_t meth = 0; meth< 11; meth++){
        delete markerEtaPrimeMass[meth];
        delete markerEtaPrimeMassMC[meth];
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


//     // **********************************************************************************************************************
//     // ******************************** Cross section for etaprime single measurement 2.76TeV ************************************
//     // **********************************************************************************************************************
//     TCanvas* canvasInvYieldEtaPrime          = new TCanvas("canvasInvYieldEtaPrime","",200,10,1350,1350*1.15);  // gives the page size
//     DrawGammaCanvasSettings( canvasInvYieldEtaPrime, 0.16, 0.02, 0.02, 0.08);
//     canvasInvYieldEtaPrime->SetLogx();
//     canvasInvYieldEtaPrime->SetLogy();
//
//     TH2F * histo2DYieldEtaPrime              = new TH2F("histo2DYieldEtaPrime","histo2DYieldEtaPrime",11000,minPtEtaPrimePlotting, maxPtEtaPrimePlotting,1000,6e-11,3e1);
//     SetStyleHistoTH2ForGraphs(histo2DYieldEtaPrime, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.65);
//     histo2DYieldEtaPrime->GetXaxis()->SetMoreLogLabels();
//     histo2DYieldEtaPrime->GetXaxis()->SetLabelOffset(-0.01);
//
//     TLatex *labelEnergyXSectionEtaPrime  = new TLatex(0.94,0.92,collisionSystem.Data());
//     SetStyleTLatex( labelEnergyXSectionEtaPrime, 0.035,4, 1, 42, kTRUE, 31);
//     TLatex *labelDetSysXSectionEtaPrime  = new TLatex(0.94,0.88,"#eta' #rightarrow #gamma#gamma");
//     SetStyleTLatex( labelDetSysXSectionEtaPrime, 0.035,4, 1, 42, kTRUE, 31);
//
//     canvasInvYieldEtaPrime->cd();
//     histo2DYieldEtaPrime->Draw("copy");
//
//         for (Int_t meth = 10; meth> -1; meth--){
//             if (graphIndEtaPrimeInvYieldSys[meth]){
//                 DrawGammaSetMarkerTGraphAsym(graphIndEtaPrimeInvYieldSys[meth], markerStyleDet[meth], markerSizeDet[meth]*0.75, colorDet[meth] , colorDet[meth], widthLinesBoxes, kTRUE);
//                 graphIndEtaPrimeInvYieldSys[meth]->Draw("E2same");
//             }
//         }
//         for (Int_t meth = 10; meth> -1; meth--){
//             if (graphIndEtaPrimeInvYieldStat[meth]){
//                 DrawGammaSetMarkerTGraphAsym(graphIndEtaPrimeInvYieldStat[meth], markerStyleDet[meth], markerSizeDet[meth]*0.75, colorDet[meth] , colorDet[meth]);
//                 graphIndEtaPrimeInvYieldStat[meth]->Draw("p,same,z");
//             }
//         }
//         labelEnergyXSectionEtaPrime->SetText(0.94,0.92,collisionSystem.Data());
//         labelEnergyXSectionEtaPrime->Draw();
//         labelDetSysXSectionEtaPrime->Draw();
//
//
//         TLegend* legendXSectionEtaPrime      = GetAndSetLegend2(0.20, 0.13,0.40,0.13+(0.035*nTotMeasEtaPrime), 0.035, 1, "", 42, 0);
//         for (Int_t meth = 0; meth < 11; meth++){
//             if (graphIndEtaPrimeInvYieldSys[meth]) legendXSectionEtaPrime->AddEntry(graphIndEtaPrimeInvYieldSys[meth],nameMeasGlobalLabel[meth],"fp");
//         }
//         legendXSectionEtaPrime->Draw();
//
//     canvasInvYieldEtaPrime->SaveAs(Form("%s/EtaPrime_InvYieldCompAllSystems.%s",outputDirSupportPaper.Data(), suffix.Data()));
//
//     canvasInvYieldEtaPrime->cd();
//     histo2DYieldEtaPrime->Draw("copy");
//         for (Int_t meth = 10; meth> -1; meth--){
//             if (graphIndEtaPrimeInvYieldSys[meth]){
//                 graphIndEtaPrimeInvYieldSys[meth]->Draw("E2same");
//             }
//         }
//         if (graphCombEtaPrimeInvYieldSys){
//             DrawGammaSetMarkerTGraphAsym(graphCombEtaPrimeInvYieldSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
//             graphCombEtaPrimeInvYieldSys->Draw("E2same");
//             legendXSectionEtaPrime->AddEntry(graphCombEtaPrimeInvYieldSys,"comb","fp");
//         }
//         for (Int_t meth = 10; meth> -1; meth--){
//             if (graphIndEtaPrimeInvYieldSys[meth]){
//                 graphIndEtaPrimeInvYieldStat[meth]->Draw("p,same,z");
//             }
//         }
//         if (graphCombEtaPrimeInvYieldStat){
//             DrawGammaSetMarkerTGraphAsym(graphCombEtaPrimeInvYieldStat, markerStyleComb, markerSizeComb, colorComb , colorComb);
//             graphCombEtaPrimeInvYieldStat->Draw("p,same,z");
//         }
//         labelEnergyXSectionEtaPrime->Draw();
//         labelDetSysXSectionEtaPrime->Draw();
//
//         legendXSectionEtaPrime->Draw();
//
//     canvasInvYieldEtaPrime->SaveAs(Form("%s/EtaPrime_InvYieldCompAllSystems_Comb.%s",outputDirSupportPaper.Data(), suffix.Data()));
//
//
    // **********************************************************************************************************************
    // ******************************** Acceptance * Efficiency for etaprime single measurement **********************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 55;
    textSizeLabelsRel                   = ((Double_t)textSizeLabelsPixel)/1200;
//     cout << textSizeLabelsRel << endl;

    TCanvas* canvasAcceptanceTimesEff   = new TCanvas("canvasAcceptanceTimesEff", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasAcceptanceTimesEff,  0.1, 0.01, 0.015, 0.095);
    canvasAcceptanceTimesEff->SetLogy(1);
    canvasAcceptanceTimesEff->SetLogx(1);

    TH1F * histo1DAccEff            = new TH1F("histo1DAccEff", "histo1DAccEff",1000, minPtEtaPrimePlotting,  maxPtEtaPrimePlotting);
    SetStyleHistoTH1ForGraphs(  histo1DAccEff, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec}#upoint#it{BR}"),
                                0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1.04);//(#times #epsilon_{pur})
    histo1DAccEff->GetYaxis()->SetRangeUser(7e-7, 7e-1 );
    histo1DAccEff->GetYaxis()->SetLabelOffset(0.001);
    histo1DAccEff->GetXaxis()->SetLabelOffset(-0.01);
    histo1DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
    TLatex *labelPerfEffi           = new TLatex(0.137,0.92,"ALICE performance");
    SetStyleTLatex( labelPerfEffi, textSizeLabelsRel,4);
    TLatex *labelEnergyEffi         = new TLatex(0.137,0.87,collisionSystem.Data());
    SetStyleTLatex( labelEnergyEffi, textSizeLabelsRel,4);
    TLatex *labelDetSysEffiEtaPrime      = new TLatex(0.137,0.82,"#eta' #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysEffiEtaPrime, textSizeLabelsRel,4);

    canvasAcceptanceTimesEff->cd();
    histo1DAccEff->DrawCopy();

        for (Int_t meth = 0; meth < 11; meth++){
            if (graphEtaPrimeEffTimesAcc[meth]){
                DrawGammaSetMarkerTGraphAsym(graphEtaPrimeEffTimesAcc[meth], markerStyleDet[meth], markerSizeDet[meth]*0.55, colorDet[meth] , colorDet[meth]);
                graphEtaPrimeEffTimesAcc[meth]->Draw("p,same,z");
            }
        }

        TLegend* legendEffiAccEtaPrime           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(nTotMeasEtaPrime*textSizeLabelsRel),textSizeLabelsPixel);
        if (graphEtaPrimeEffTimesAcc[0]) legendEffiAccEtaPrime->AddEntry(graphEtaPrimeEffTimesAcc[0],nameMeasGlobalLabel[0],"p");
        if (graphEtaPrimeEffTimesAcc[4]) legendEffiAccEtaPrime->AddEntry(graphEtaPrimeEffTimesAcc[4],nameMeasGlobalLabel[4],"p");
        if (graphEtaPrimeEffTimesAcc[2]) legendEffiAccEtaPrime->AddEntry(graphEtaPrimeEffTimesAcc[2],nameMeasGlobalLabel[2],"p");
        if (graphEtaPrimeEffTimesAcc[1]) legendEffiAccEtaPrime->AddEntry(graphEtaPrimeEffTimesAcc[1],nameMeasGlobalLabel[1],"p");
        if (graphEtaPrimeEffTimesAcc[3]) legendEffiAccEtaPrime->AddEntry(graphEtaPrimeEffTimesAcc[3],nameMeasGlobalLabel[3],"p");
        if (graphEtaPrimeEffTimesAcc[5]) legendEffiAccEtaPrime->AddEntry(graphEtaPrimeEffTimesAcc[5],nameMeasGlobalLabel[5],"p");
        legendEffiAccEtaPrime->Draw();

        labelPerfEffi->Draw();
        labelEnergyEffi->SetText(0.137,0.87,collisionSystem.Data());
        labelEnergyEffi->Draw();
        labelDetSysEffiEtaPrime->Draw();

    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/EtaPrime_AcceptanceTimesEff.%s",outputDirSupportPaper.Data(), suffix.Data()));



//     fileFitsOutput << "*******************************************************************************************" << endl;
//     fileFitsOutput << "****************************** Power law fit etaprime ******************************************" << endl;
//     fileFitsOutput << "*******************************************************************************************" << endl;
//     TF1* fitPowInvYieldEtaPrimeTot       = NULL;
//     TF1* fitPowInvYieldEtaPrimeStat      = NULL;
//     TF1* fitOHagInvYieldEtaPrimeTot      = NULL;
//
//     fitPowInvYieldEtaPrimeTot   = FitObject("powPure","fitPowInvYieldEtaPrimeTot","EtaPrime",graphCombEtaPrimeInvYieldTot,4,20. ,NULL,"QNRMEX0+","", kFALSE);
//     cout << WriteParameterToFile(fitPowInvYieldEtaPrimeTot)<< endl;
//     fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldEtaPrimeTot)<< endl;
//     fitPowInvYieldEtaPrimeStat   = FitObject("powPure","fitPowInvYieldEtaPrime","EtaPrime",graphCombEtaPrimeInvYieldStat,4,20. ,NULL,"QNRMEX0+","", kFALSE);
//     cout << WriteParameterToFile(fitPowInvYieldEtaPrimeStat)<< endl;
//     fileFitsOutput <<  WriteParameterToFile(fitPowInvYieldEtaPrimeStat)<< endl;
//     fitOHagInvYieldEtaPrimeTot   = FitObject("oHag","fitOHagInvYieldEtaPrime","EtaPrime",graphCombEtaPrimeInvYieldTot,0.3,20. ,NULL,"QNRMEX0+","", kFALSE);
//     cout << WriteParameterToFile(fitOHagInvYieldEtaPrimeTot)<< endl;
//     fileFitsOutput <<  WriteParameterToFile(fitOHagInvYieldEtaPrimeTot)<< endl;
//
//
//
//     // **********************************************************************************************************************
//     // ************************* Saving of final results ********************************************************************
//     // **********************************************************************************************************************
//
//     TString nameOutputCommonFile    = Form("CombinedResultsPaper13TeV%s.root", dateForOutput.Data());
//
//     TFile fCombResults(nameOutputCommonFile.Data(), "UPDATE");
//
//     TDirectoryFile* directoryEtaPrimeOutput  = NULL;
//     directoryEtaPrimeOutput                  = (TDirectoryFile*)fCombResults.Get(Form("EtaPrime13TeV"));
//     if (!directoryEtaPrimeOutput){
//         fCombResults.mkdir(Form("EtaPrime13TeV"));
//         directoryEtaPrimeOutput              = (TDirectoryFile*)fCombResults.Get(Form("EtaPrime13TeV"));
//     }
//     fCombResults.cd(Form("EtaPrime13TeV"));
//
//         // Final spectrum
//         if (graphCombEtaPrimeInvYieldTot) graphCombEtaPrimeInvYieldTot->Write("graphInvYieldEtaPrimeCombTotErr",TObject::kOverwrite);
//         if (graphCombEtaPrimeInvYieldStat) graphCombEtaPrimeInvYieldStat->Write("graphInvYieldEtaPrimeCombStatErr",TObject::kOverwrite);
//         if (graphCombEtaPrimeInvYieldSys) graphCombEtaPrimeInvYieldSys->Write("graphInvYieldEtaPrimeCombSysErr",TObject::kOverwrite);
//             // fits for etaprime
//         if (fitInvYieldEtaPrime) fitInvYieldEtaPrime->Write("TsallisFitEtaPrime");
//         if (fitTCMInvYieldEtaPrime) fitTCMInvYieldEtaPrime->Write("TwoComponentModelFitEtaPrime",TObject::kOverwrite);
//
//         // Final inv yield INEL
//         if (bWCorrection.Contains("Y")){
//             // Final spectrum correlations Method A
//             if(graphCombEtaPrimeInvYieldTot_yShifted)graphCombEtaPrimeInvYieldTot_yShifted->Write("graphInvYieldEtaPrimeCombTotErr_yShifted",TObject::kOverwrite);
//             if(graphCombEtaPrimeInvYieldStat_yShifted)graphCombEtaPrimeInvYieldStat_yShifted->Write("graphInvYieldEtaPrimeCombStatErr_yShifted",TObject::kOverwrite);
//             if(graphCombEtaPrimeInvYieldSys_yShifted)graphCombEtaPrimeInvYieldSys_yShifted->Write("graphInvYieldEtaPrimeCombSysErr_yShifted",TObject::kOverwrite);
//         }
//         // writing individual measurements
//         for (Int_t meth = 0; meth< 11; meth++){
//             if (graphIndEtaPrimeInvYieldStat){
//                 if (graphIndEtaPrimeInvYieldStat[meth]) graphIndEtaPrimeInvYieldStat[meth]->Write(Form("graphInvYieldEtaPrime%sStatErr", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
//                 if (graphIndEtaPrimeInvYieldSys[meth]) graphIndEtaPrimeInvYieldSys[meth]->Write(Form("graphInvYieldEtaPrime%sSysErr", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
//                 if (bWCorrection.Contains("Y")){
//                     if (graphIndEtaPrimeInvYieldStat_yShifted[meth]) graphIndEtaPrimeInvYieldStat_yShifted[meth]->Write(Form("graphInvYieldEtaPrime%sStatErr_yShifted", nameMeasGlobalLabel[meth].Data()),
//                                                                                                                         TObject::kOverwrite);
//                     if (graphIndEtaPrimeInvYieldSys_yShifted[meth]) graphIndEtaPrimeInvYieldSys_yShifted[meth]->Write(Form("graphInvYieldEtaPrime%sSysErr_yShifted", nameMeasGlobalLabel[meth].Data()),
//                                                                                                                         TObject::kOverwrite);
//                 }
//             }
//         }
//
//         TDirectoryFile* supporting  = NULL;
//         supporting                  = (TDirectoryFile*)directoryEtaPrimeOutput->Get("Supporting");
//         if (!supporting)
//             directoryEtaPrimeOutput->mkdir("Supporting");
//         directoryEtaPrimeOutput->cd("Supporting");
//             // Writing full correction factors
//             for (Int_t meth = 0; meth< 11; meth++){
//                 if (graphEtaPrimeEffTimesAcc[meth]) graphEtaPrimeEffTimesAcc[meth]->Write(Form("EtaPrimeCorrectionFactor%s", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
//                 if (graphEtaPrimeMass[meth]) graphEtaPrimeMass[meth]->Write(Form("EtaPrimeMassData%s", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
//                 if (graphEtaPrimeMassMC[meth]) graphEtaPrimeMassMC[meth]->Write(Form("EtaPrimeMassMC%s", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
//                 if (graphEtaPrimeWidth[meth]) graphEtaPrimeWidth[meth]->Write(Form("EtaPrimeWidthData%s", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
//                 if (graphEtaPrimeWidthMC[meth]) graphEtaPrimeWidthMC[meth]->Write(Form("EtaPrimeWidthMC%s", nameMeasGlobalLabel[meth].Data()), TObject::kOverwrite);
//             }
//
//     fCombResults.Close();
//
//
//     //  **********************************************************************************************************************
//     //  ************************* Saving only fits to final results **********************************************************
//     //  **********************************************************************************************************************
//
//     TString nameOutputCommonFileFitsOnly    = Form("FitsPaper13TeVCent_%s.root", dateForOutput.Data());
//     TFile fFitsResults(nameOutputCommonFileFitsOnly.Data(), "RECREATE");
//
//         // fits for etaprime
//         fitInvYieldEtaPrime->Write("TsallisFitEtaPrime");
//         fitTCMInvYieldEtaPrime->Write("TwoComponentModelFitEtaPrime");
//     fFitsResults.Close();
}

