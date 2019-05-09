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

void CombineMesonMeasurementsPbPb5TeVCent(  TString fileNamePCM             = "",
                                            TString fileNameEMCAL           = "",
                                            TString fileNamePHOS            = "",
                                            TString fileNamePCMEMCAL        = "",
                                            TString fileNamePCMPHOS         = "",
                                            TString suffix                  = "eps",
                                            Bool_t isMC                     = kFALSE,
                                            TString bWCorrection            = "X",
                                            TString fileNameCorrFactors     = "",
                                            TString fileNameReference   = "",
                                            TString pathConfigsRAAErr       = ""
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
    TString collisionSystemPbPb                 = "Pb-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
    TString collisionSystempp                   = "pp, #sqrt{#it{s}} = 5.02 TeV";
    TString outputDir                           = Form("%s/%s/CombineMesonMeasurements%sPbPb5TeVCent", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());
    TString outputDirFile                       = Form("%s/%s/CombineMesonMeasurements%sPbPb5TeVCent/Inputs", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());
    TString outputDirSupportComb                = Form("%s/%s/CombineMesonMeasurements%sPbPb5TeVCent/SupportingCombination", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());
    TString outputDirSupport                    = Form("%s/%s/CombineMesonMeasurements%sPbPb5TeVCent/Supporting", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());
    TString outputDirSupportPaper               = Form("%s/%s/CombineMesonMeasurements%sPbPb5TeVCent/SupportingPaper", suffix.Data(), dateForOutput.Data(), bWCorrection.Data());

    TString nameFinalResDat                     = Form("%s/CombinedResults%s_FitResults.dat", outputDir.Data(), bWCorrection.Data());

    fstream  fileFitsOutput;
    fileFitsOutput.open(nameFinalResDat.Data(), ios::out);

    TString fileNameTheory                      = "ExternalInputPbPb/Theory/TheoryCompilationPbPb.root";
    TString fileNameOtherParticleInput          = "ExternalInputPbPb/IdentifiedCharged/ChargedPion_29_Apr_2019_PbPb5TeV.root";

    Int_t maxCent                               = 9;
    Int_t intRefMB                              = 1;
    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();

    Width_t  widthLinesBoxes                    = 1;
    Width_t  widthCommonFit                     = 2;

    // Definition of colors, styles and markers sizes
    Color_t  colorComb                          = kMagenta+2;
    Style_t  markerStyleComb                    = 20;
    Size_t   markerSizeComb                     = 2;

    Color_t  colorCombLowPt                     = 807;
    Color_t  colorCombHighPt                    = kGreen+2;
    Style_t  markerStyleCombLowPt               = 20;
    Style_t  markerStyleCombHighPt              = 20;
    Size_t   markerSizeComparison               = 2;

    TString  centArray[9]                       = { "0-10%", "10-20%", "20-40%", "40-60%", "60-80%", "0-5%", "5-10%", "0-20%", "20-50%"};
    TString  centArray2[9]                      = { "0-10%", "10-20%", "20-40%", "40-60%", "60-80%", "0-5%", "5-10%", "0-20%", "20-50%"};
    TString  centArray3[9]                      = { "0-10%", "10-20%", "20-40%", "40-60%", "60-80%", "0-5%", "5-10%", "0-20%", "20-50%"};
    TString  centArrayOutput[9]                 = { "0010", "1020", "2040", "4060", "6080",  "0005", "0510", "0020", "2050"};
    TString  centArrayCorr[9]                   = { "0010/", "1020/", "2040/", "4060/",  "6080/",  "0005/", "0510/", "0020/", "2050/"};
    TString addCentString[9]                    = {"_V0M", "_V0M", "_V0M","_V0M","_V0M",               "_V0M", "_V0M", "_V0M", "_V0M"};
    Bool_t  enableCent[9]                       = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,   kFALSE, kFALSE, kFALSE, kFALSE};
    Bool_t  enableCentComb[9]                   = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,   kFALSE, kFALSE, kFALSE, kFALSE};
    Bool_t  enableCentRAA[9]                    = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,   kFALSE, kFALSE, kFALSE, kFALSE};
    TString nameCentEst                         =  "V0M";

    TString  nameMeasGlobal[11]                 = { "PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal",
                                                    "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged",
                                                    "PCMOtherDataset"};
    TString  nameMeasGlobalLabel[11]            = { "PCM", "PHOS", "EMC", "PCM-PHOS", "PCM-EMC",
                                                    "PCM-Dal", "PHOS-Dal", "EMC-Dal", "EMChigh", "mEMC",
                                                    "PCMOtherDataset"};
    TString  nameConfigFileRAAErr[9]            = { "", "", "", "", "",  "", "", "", ""};

    TString nameTheory[3]                       = {"Paquet", "SHM-EQ", "SHM-NEQ"};
    TString nameTheoryRead[3]                   = {"Paquet", "SHM_EQ", "SHM_NEQ"};
    TString nameTheoryRAA[3]                    = {"Djordjevic", "Djordjevic T_{cont}", "Vitev"};
    TString nameTheoryRAARead[3]                = {"DjordjevicBjorken", "DjordjevicConstTemp", "Vitev"};
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
    if (fileNameReference.CompareTo("")!=0 )    gSystem->Exec(Form("cp %s %s/InputReference.root", fileNameReference.Data(), outputDirFile.Data()));
    if (fileNameCorrFactors.CompareTo("")!=0 )  gSystem->Exec(Form("cp %s %s/InputCorrFactors.root", fileNameCorrFactors.Data(), outputDirFile.Data()));
    if (fileNameTheory.CompareTo("")!= 0)       gSystem->Exec(Form("cp %s %s/InputTheory.root", fileNameTheory.Data(), outputDirFile.Data()));
    if (pathConfigsRAAErr.CompareTo("")!=0 ){
        for (Int_t cent = 0; cent < maxCent; cent++){
            if (!enableCentRAA[cent]) continue;
            nameConfigFileRAAErr[cent]               = Form("%s/configFileRAA_%s%s.txt", pathConfigsRAAErr.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data());
            gSystem->Exec(Form("cp %s %s/configFileRAA_%s%s.txt", nameConfigFileRAAErr[cent].Data(), outputDirFile.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data()));
            cout << "copied: " << nameConfigFileRAAErr[cent].Data() << " if available" << endl;
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
    Int_t modeToBeUsed[9]                       = { 20, 20, 20, 20, 20,     20, 20, 20, 20};
    Double_t branchingRatioPi0[11]              =  { 0.98823,  0.98823,  0.98823,  0.98823,  0.98823,
                                                     0.01174,  0.01174,  0.01174,  0.98823,  0.98823,
                                                     0.98823 };
    Double_t branchingRatioEta[11]              =  { 0.3931,  0.3931,  0.3931,  0.3931,  0.3931,
                                                     0.000068,  0.000068,  0.000068,  0.3931,  0.3931,
                                                     0.3931 };

    Double_t minPtPi0Plotting                   = 0.23;
    Double_t maxPtPi0Plotting                   = 51.;
    Double_t minPtEtaPlotting                   = 0.33;
    Double_t maxPtEtaPlotting                   = 40.;
    Double_t minPtEtaToPi0Plotting              = 0.33;
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

    Color_t  colorCent[9];
    Color_t  colorCentMC[9];
    Color_t  colorCentBox[9];
    Style_t  markerStyleCent[9];
    Style_t  markerStyleCentMC[9];
    Size_t   markerSizeCent[9];
    Size_t   markerSizeCentMC[9];
    Double_t nCollPbPb[9];
    Double_t nCollErrPbPb[9];
    Double_t tPbPb[9];
    Double_t tPbPbErr[9];
    for (Int_t cent = 0; cent < maxCent; cent++){
        colorCent[cent]                         = GetColorDefaultColor("PbPb_5.02TeV", "", centArray[cent]);
        colorCentMC[cent]                       = GetColorDefaultColor("PbPb_5.02TeV", "HIJING", centArray[cent]);
        colorCentBox[cent]                      = GetColorDefaultColor("PbPb_5.02TeV", "", centArray[cent],kTRUE);
        markerStyleCent[cent]                   = GetDefaultMarkerStyle("PbPb_5.02TeV", "", centArray[cent]);
        markerStyleCentMC[cent]                 = GetDefaultMarkerStyle("PbPb_5.02TeV", "HIJING", centArray[cent]);

        markerSizeCent[cent]                    = GetDefaultMarkerSize("PbPb_5.02TeV", "", centArray[cent])*2;
        markerSizeCentMC[cent]                  = GetDefaultMarkerSize("PbPb_5.02TeV", "HIJING", centArray[cent])*2;
        nCollPbPb[cent]                         = GetNCollFromName(centArrayOutput[cent], "PbPb_5.02TeV");
        nCollErrPbPb[cent]                      = GetNCollErrFromName(centArrayOutput[cent], "PbPb_5.02TeV");
        tPbPb[cent]                             = GetTAAFromName(centArrayOutput[cent], "PbPb_5.02TeV")*1e3*(1/recalcBarn);
        tPbPbErr[cent]                          = GetTAAErrFromName(centArrayOutput[cent], "PbPb_5.02TeV")*1e3*(1/recalcBarn);
        cout << markerStyleCent[cent] << "\t" << markerStyleCentMC[cent] << "\t"<<markerSizeCent[cent] << "\t" << markerSizeCentMC[cent] << endl;
        clog << cent << "\t" << centArrayOutput[cent]+addCentString[cent] << "\t"<< nCollPbPb[cent] << "\t" << nCollErrPbPb[cent] << "\t"<< tPbPb[cent] << "\t"<< tPbPbErr[cent] << endl;;
    }

    Double_t xSection5TeV                       = ReturnCorrectXSection("5TeV", 1);
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
    Double_t xPtLimitsPi0[9][100];
    Int_t maxNBinsPi0Abs                            = 0;
    Int_t maxNBinsPi0[9]                           = {0};

    cout << "Setting Eta binning" << endl;
    Double_t xPtLimitsEta[9][100];
    Double_t xPtLimitsEtaToPi0[9][100];
    Int_t maxNBinsEtaAbs                            = 0;
    Int_t maxNBinsEta[9]                           = {0};
    Int_t maxNBinsEtaToPi0Abs                       = 0;
    Int_t maxNBinsEtaToPi0[9]                      = {0};
    Int_t nTotMeasPi0[9]                           = {0};
    Int_t nTotMeasEta[9]                           = {0};
    Int_t nTotMeasEtaToPi0[9]                      = {0};
    TString fileNamesMethod[11]                     = {"", "", "", "", "", "", "", "", "", "", ""};
    TString fileNamesPbPbPi0DetailedSys[9][11];
    TString fileNamesppPi0DetailedSys[9][11];
    TString fileNamesRAAPi0DetailedSys[9][11];
    TString fileNamesPbPbEtaDetailedSys[9][11];
    TString fileNamesppEtaDetailedSys[9][11];
    TString fileNamesRAAEtaDetailedSys[9][11];
    Bool_t havePi0SysDetailedPbPb[9][11];
    Bool_t havePi0SysDetailedpp[9][11];
    Bool_t haveEtaSysDetailedPbPb[9][11];
    Bool_t haveEtaSysDetailedpp[9][11];
    vector<TString> ptSysRemNames[9][11];
    for (Int_t cent = 0; cent < maxCent; cent++){
        for (Int_t meth = 0; meth < 11; meth++){
            fileNamesPbPbPi0DetailedSys[cent][meth]      = "";
            fileNamesppPi0DetailedSys[cent][meth]       = "";
            fileNamesRAAPi0DetailedSys[cent][meth]     = "";
            fileNamesPbPbEtaDetailedSys[cent][meth]      = "";
            fileNamesppEtaDetailedSys[cent][meth]       = "";
            fileNamesRAAEtaDetailedSys[cent][meth]     = "";
            havePi0SysDetailedPbPb[cent][meth]           = kFALSE;
            havePi0SysDetailedpp[cent][meth]            = kFALSE;
            haveEtaSysDetailedPbPb[cent][meth]           = kFALSE;
            haveEtaSysDetailedpp[cent][meth]            = kFALSE;
        }
    }

    Bool_t haveEffSecCorr[9][4][11];
    TFile* fileMethod[11];
    TDirectory* directoryPi0[9][11];
    TGraphAsymmErrors* graphPi0EffSecCorrFromX[9][4][11];
    TGraphAsymmErrors* graphPi0Eff[9][11];
    TGraphAsymmErrors* graphPi0Acc[9][11];
    TGraphAsymmErrors* graphPi0EffTimesAcc[9][11];
    TGraphAsymmErrors* graphPi0Mass[9][11];
    TGraphAsymmErrors* graphPi0MassMC[9][11];
    TGraphAsymmErrors* graphPi0Width[9][11];
    TGraphAsymmErrors* graphPi0WidthMC[9][11];
    TH1D* histoPi0InvYieldStat[9][11];
    TGraphAsymmErrors* graphPi0InvYieldStat[9][11];
    TGraphAsymmErrors* graphPi0InvYieldSys[9][11];
    TDirectory* directoryEta[9][11];
    TGraphAsymmErrors* graphEtaEff[9][11];
    TGraphAsymmErrors* graphEtaAcc[9][11];
    TGraphAsymmErrors* graphEtaEffTimesAcc[9][11];
    TGraphAsymmErrors* graphEtaMass[9][11];
    TGraphAsymmErrors* graphEtaMassMC[9][11];
    TGraphAsymmErrors* graphEtaWidth[9][11];
    TGraphAsymmErrors* graphEtaWidthMC[9][11];
    TH1D* histoEtaInvYieldStat[9][11];
    TGraphAsymmErrors* graphEtaInvYieldStat[9][11];
    TGraphAsymmErrors* graphEtaInvYieldSys[9][11];
    TH1D* histoEtaToPi0Stat[9][11];
    TGraphAsymmErrors* graphEtaToPi0Stat[9][11];
    TGraphAsymmErrors* graphEtaToPi0Sys[9][11];

    for (Int_t cent = 0; cent < maxCent; cent++){
        maxNBinsPi0[cent]                           = GetBinning( xPtLimitsPi0[cent], maxNBinsPi0Abs, "Pi0", Form("PbPb_5.02TeV"), modeToBeUsed[cent], -1, kFALSE, centArray[cent]);
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
        maxNBinsEta[cent]                           = GetBinning( xPtLimitsEta[cent], maxNBinsEtaAbs, "Eta", Form("PbPb_5.02TeV"), modeToBeUsed[cent], -1, kFALSE, centArray[cent]);
        maxNBinsEtaToPi0[cent]                      = GetBinning( xPtLimitsEtaToPi0[cent], maxNBinsEtaToPi0Abs, "Eta", Form("PbPb_5.02TeV"), modeToBeUsed[cent], -1, kFALSE, centArray[cent]);
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

    for (Int_t cent = 0; cent < maxCent; cent++){

        if (nameConfigFileRAAErr[cent].CompareTo("") != 0){
            cout << "===================================================================================================================" << endl;
            cout << "RAA details for " << centArray[cent].Data() << "\t" << addCentString[cent].Data() << "\t" << endl;
            ifstream inPbConfig(nameConfigFileRAAErr[cent].Data());
            Int_t nReadSys  = 0;
            while(!inPbConfig.eof() && nReadSys<11 ){
                cout << "read line:" <<  nReadSys << "\t"<< nameMeasGlobalLabel[nReadSys].Data() << endl;
                TString nameCurrentSys  = "";
                inPbConfig >> nameCurrentSys ;
                if (nameCurrentSys.CompareTo(nameMeasGlobalLabel[nReadSys].Data()) != 0){
                    cout << "wrong order in configuration file, was expecting " <<  nameMeasGlobalLabel[nReadSys].Data() << endl;
                    return;
                } else {
                    inPbConfig >> fileNamesPbPbPi0DetailedSys[cent][nReadSys] >> fileNamesppPi0DetailedSys[cent][nReadSys] >> fileNamesPbPbEtaDetailedSys[cent][nReadSys] >> fileNamesppEtaDetailedSys[cent][nReadSys];
                    if (fileNamesPbPbPi0DetailedSys[cent][nReadSys].CompareTo("bla") != 0){
                        cout << "pi0 PbPb sys detailed: " << fileNamesPbPbPi0DetailedSys[cent][nReadSys] << endl;
                        havePi0SysDetailedPbPb[cent][nReadSys]         = kTRUE;
                    } else {
                        fileNamesPbPbPi0DetailedSys[cent][nReadSys]    = "";
                    }
                    if (fileNamesppPi0DetailedSys[cent][nReadSys].CompareTo("bla") != 0){
                        cout << "pi0 pp sys detailed: " << fileNamesppPi0DetailedSys[cent][nReadSys] << endl;
                        havePi0SysDetailedpp[cent][nReadSys]          = kTRUE;
                    } else {
                        fileNamesppPi0DetailedSys[cent][nReadSys]     = "";
                    }
                    if (fileNamesPbPbEtaDetailedSys[cent][nReadSys].CompareTo("bla") != 0){
                        cout << "eta PbPb sys detailed: " << fileNamesPbPbEtaDetailedSys[cent][nReadSys] << endl;
                        haveEtaSysDetailedPbPb[cent][nReadSys]         = kTRUE;
                    } else {
                        fileNamesPbPbEtaDetailedSys[cent][nReadSys]    = "";
                    }
                    if (fileNamesppEtaDetailedSys[cent][nReadSys].CompareTo("bla") != 0){
                        cout << "eta pp sys detailed: " << fileNamesppEtaDetailedSys[cent][nReadSys] << endl;
                        haveEtaSysDetailedpp[cent][nReadSys]          = kTRUE;
                    } else {
                        fileNamesppEtaDetailedSys[cent][nReadSys]     = "";
                    }
                    if (!(havePi0SysDetailedPbPb[cent][nReadSys] && havePi0SysDetailedpp[cent][nReadSys])){
                        havePi0SysDetailedPbPb[cent][nReadSys]         = kTRUE;
                        havePi0SysDetailedpp[cent][nReadSys]          = kTRUE;
                        fileNamesPbPbPi0DetailedSys[cent][nReadSys]    = "";
                        fileNamesppPi0DetailedSys[cent][nReadSys]     = "";
                    } else {
                        fileNamesRAAPi0DetailedSys[cent][nReadSys]  = Form("%s/Pi0RAA_%s_detailedSys%s%s.dat", outputDirSupportComb.Data(),nameMeasGlobalLabel[nReadSys].Data(),
                                                                            centArrayOutput[cent].Data(), addCentString[cent].Data());
                    }
                    if (!(haveEtaSysDetailedPbPb[cent][nReadSys] && haveEtaSysDetailedpp[cent][nReadSys])){
                        haveEtaSysDetailedPbPb[cent][nReadSys]         = kTRUE;
                        haveEtaSysDetailedpp[cent][nReadSys]          = kTRUE;
                        fileNamesPbPbEtaDetailedSys[cent][nReadSys]    = "";
                        fileNamesppEtaDetailedSys[cent][nReadSys]     = "";
                    } else {
                        fileNamesRAAEtaDetailedSys[cent][nReadSys]  = Form("%s/EtaRAA_%s_detailedSys%s%s.dat", outputDirSupportComb.Data(),nameMeasGlobalLabel[nReadSys].Data(),
                                                                            centArrayOutput[cent].Data(), addCentString[cent].Data());
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
    }


    // *******************************************************************************************************
    // ***************************** Read in data for different methods **************************************
    // *******************************************************************************************************
    for (Int_t cent = 0; cent < maxCent; cent ++){
        if (!enableCent[cent]) continue;
        cout << "**********************************************************************" << endl;
        cout << "Reading in cent " << centArray[cent].Data() << endl;
        cout << "**********************************************************************" << endl;

        for (Int_t meth = 0; meth < 11; meth++){
            if (fileNamesMethod[meth].CompareTo("") != 0){
                cout << "Reading in " << nameMeasGlobalLabel[meth].Data() << " from " << fileNamesMethod[meth].Data() << endl;
                if (!fileMethod[meth] || cent == maxCent)
                    fileMethod[meth]                               = new TFile(fileNamesMethod[meth].Data());
                if (!fileMethod[meth]) {
                    cout << "file " << fileNamesMethod[meth].Data() << " not found! Skipping " << nameMeasGlobalLabel[meth].Data()  << endl;
                    continue;
                }
                TString nameDirPi0                                      = Form("Pi0%sPbPb_5.02TeV%s",centArray[cent].Data(), addCentString[cent].Data());
//                 if (meth == 1)
//                     nameDirPi0                                          = Form("Pi05TeV_%s%s",centArrayReadPHOS[cent].Data(), addCentString[cent].Data());
                directoryPi0[cent][meth]                                = (TDirectory*)fileMethod[meth]->Get(nameDirPi0.Data());

                if (!directoryPi0[cent][meth]) {
                    cout << "File doesn't contain directory " << nameDirPi0.Data() << "! Skipping " << nameMeasGlobalLabel[meth].Data()  << endl;
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
                    if (graphPi0WidthMC[cent][meth]) graphPi0WidthMC[cent][meth]                         = ScaleGraph(graphPi0WidthMC[cent][meth], 1000.);
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
                    cout << graphPi0InvYieldSys[cent][meth] << endl;
                    cout << "Pi0 sys error" << endl;
                    if (graphPi0InvYieldSys[cent][meth]) graphPi0InvYieldSys[cent][meth]->Print();

                }

                TString nameDirEta                                      = Form("Eta%sPbPb_5.02TeV%s",centArray[cent].Data(), addCentString[cent].Data());
//                 if (meth == 1)
//                     nameDirEta                                          = Form("Eta5TeV_%s%s",centArrayReadPHOS[cent].Data(), addCentString[cent].Data());
                directoryEta[cent][meth]                                = (TDirectory*)fileMethod[meth]->Get(nameDirEta.Data());

                if (!directoryEta[cent][meth]) {
                    cout << "File doesn't contain directory " << Form("Eta%sPbPb_5.02TeV%s",centArray[cent].Data(), addCentString[cent].Data()) << "! Skipping " << nameMeasGlobalLabel[meth].Data()  << endl;
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
    TH1F* histoDPMJetPi0[9]                             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1F* histoDPMJetEta[9]                             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1F* histoDPMJetEtaToPi0[9]                        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1F* histoDPMJetPi0ToPiCh[9]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1F* histoDPMJetEtaToKCh[9]                        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1F* histoEPOSLHCPi0[9]                            = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1F* histoEPOSLHCEta[9]                            = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1F* histoEPOSLHCEtaToPi0[9]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1F* histoEPOSLHCPi0ToPiCh[9]                      = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1F* histoEPOSLHCEtaToKCh[9]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1F* histoHIJINGPi0[9]                             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1F* histoHIJINGEta[9]                             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1F* histoHIJINGEtaToPi0[9]                        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1F* histoHIJINGPi0ToPiCh[9]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1F* histoHIJINGEtaToKCh[9]                        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraph* graphTheoryPi0[3][9]                        = { { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                            { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                            { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL} };
    TGraph* graphTheoryEta[3][9]                        = { { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                            { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                            { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL} };
    TGraph* graphTheoryEtaToPi0[3][9]                   = { { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                            { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                            { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL} };
    TGraphErrors* graphTheoryRAAPi0[3][9]                     = { { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                            { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                            { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL} };
    TGraphErrors* graphTheoryRAAEta[3][9]               = { { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                            { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                                            { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL} };

    TDirectory* directoryTheoryPi0                      = (TDirectory*)fileTheory->Get("Pi0_PbPb_5.02TeV");
    TDirectory* directoryTheoryEta                      = (TDirectory*)fileTheory->Get("Eta_PbPb_5.02TeV");
    if (!directoryTheoryPi0 && !directoryTheoryEta){
        cout << "ATTENTION: no theory curves available!" << endl;
    }

    for (Int_t cent = 0; cent < maxCent; cent ++){
        if (!enableCent[cent] || cent > maxCent) continue;
        cout << "**********************************************************************" << endl;
        cout << "Reading in cent " << centArray[cent].Data() << endl;
        cout << "**********************************************************************" << endl;
        TString currentCent                             =  centArrayOutput[cent].Data();
        for (Int_t theo = 0; theo < 3; theo++){
            graphTheoryPi0[theo][cent]                          = (TGraph*)fileTheory->Get(Form("Pi0_PbPb_5.02TeV/Spectra_%s_%s", nameTheoryRead[theo].Data(), currentCent.Data()));
            graphTheoryEta[theo][cent]                          = (TGraph*)fileTheory->Get(Form("Eta_PbPb_5.02TeV/Spectra_%s_%s", nameTheoryRead[theo].Data(), currentCent.Data()));
            graphTheoryEtaToPi0[theo][cent]                     = (TGraph*)fileTheory->Get(Form("Eta_PbPb_5.02TeV/EtaToPi0_%s_%s", nameTheoryRead[theo].Data(), currentCent.Data()));
            if (graphTheoryPi0[theo][cent] || graphTheoryEta[theo][cent] || graphTheoryEtaToPi0[theo][cent]){
                cout << "Found input for theory " << nameTheoryRead[theo].Data() << " for centrality " << centArray[cent].Data() << endl;
                if (graphTheoryPi0[theo][cent]) cout << "pi0 available \t" ;
                if (graphTheoryEta[theo][cent]) cout << "eta available \t" ;
                if (graphTheoryEtaToPi0[theo][cent]) cout << "eta/pi0 available \t" ;
                cout << endl;
            }
            graphTheoryRAAPi0[theo][cent]                       = (TGraphErrors*)fileTheory->Get(Form("Pi0_PbPb_5.02TeV/RAA_%s_%s", nameTheoryRAARead[theo].Data(), currentCent.Data()));
            graphTheoryRAAEta[theo][cent]                       = (TGraphErrors*)fileTheory->Get(Form("Eta_PbPb_5.02TeV/RAA_%s_%s", nameTheoryRAARead[theo].Data(), currentCent.Data()));
            if (graphTheoryRAAPi0[theo][cent] || graphTheoryRAAEta[theo][cent]){
                cout << "Found input for RAA theory " << nameTheoryRAARead[theo].Data() << " for centrality " << centArray[cent].Data() << endl;
                if (graphTheoryRAAPi0[theo][cent]) cout << "pi0 available \t" ;
                if (graphTheoryRAAEta[theo][cent]) cout << "eta available \t" ;
                cout << endl;
            }
        }

//         if (centArray[cent].CompareTo("0-100%") == 0){
//             histoDPMJetPi0[cent]                                = (TH1F*) directoryTheory[cent]->Get("histoPi0SpecDPMJet_MCGen");
//             histoDPMJetEta[cent]                                = (TH1F*) directoryTheory[cent]->Get("histoEtaSpecDPMJet_MCGen");
//             histoDPMJetEtaToPi0[cent]                           = (TH1F*) directoryTheory[cent]->Get("histoEtaToPi0DPMJet_MCGen");
//             histoDPMJetPi0ToPiCh[cent]                          = (TH1F*) directoryTheory[cent]->Get("histoPi0ToPiChDPMJet_MCGen");
//             histoEPOSLHCPi0[cent]                               = (TH1F*) directoryTheory[cent]->Get("histoPi0SpecEPOSLHC_MCGen");
//             histoEPOSLHCEta[cent]                               = (TH1F*) directoryTheory[cent]->Get("histoEtaSpecEPOSLHC_MCGen");
//             histoEPOSLHCEtaToPi0[cent]                          = (TH1F*) directoryTheory[cent]->Get("histoEtaToPi0EPOSLHC_MCGen");
//             histoEPOSLHCPi0ToPiCh[cent]                         = (TH1F*) directoryTheory[cent]->Get("histoPi0ToPiChEPOSLHC_MCGen");
//             histoHIJINGPi0[cent]                                = (TH1F*) directoryTheory[cent]->Get("histoPi0SpecHIJING_MCGen");
//             histoHIJINGEta[cent]                                = (TH1F*) directoryTheory[cent]->Get("histoEtaSpecHIJING_MCGen");
//             histoHIJINGEtaToPi0[cent]                           = (TH1F*) directoryTheory[cent]->Get("histoEtaToPi0HIJING_MCGen");
//             histoHIJINGPi0ToPiCh[cent]                          = (TH1F*) directoryTheory[cent]->Get("histoPi0ToPiChHIJING_MCGen");
//             histoHIJINGEtaToKCh[cent]                           = (TH1F*) directoryTheory[cent]->Get("histoEtaToKChHIJING");
//         } else {
//             histoDPMJetPi0[cent]                                = (TH1F*) directoryTheory[cent]->Get("histoPi0SpecDPMJet_Reb");
//             histoDPMJetEta[cent]                                = (TH1F*) directoryTheory[cent]->Get("histoEtaSpecDPMJet_Reb");
//             histoDPMJetEtaToPi0[cent]                           = (TH1F*) directoryTheory[cent]->Get("histoEtaToPi0DPMJet");
//             histoDPMJetPi0ToPiCh[cent]                          = (TH1F*) directoryTheory[cent]->Get("histoPi0ToPiChDPMJet");
//             histoEPOSLHCPi0[cent]                               = (TH1F*) directoryTheory[cent]->Get("histoPi0SpecEPOSLHC_Reb");
//             histoEPOSLHCEta[cent]                               = (TH1F*) directoryTheory[cent]->Get("histoEtaSpecEPOSLHC_Reb");
//             histoEPOSLHCEtaToPi0[cent]                          = (TH1F*) directoryTheory[cent]->Get("histoEtaToPi0EPOSLHC");
//             histoEPOSLHCPi0ToPiCh[cent]                         = (TH1F*) directoryTheory[cent]->Get("histoPi0ToPiChEPOSLHC");
//         }
//         histoDPMJetEtaToKCh[cent]                           = (TH1F*) directoryTheory[cent]->Get("histoEtaToKChDPMJet");
//         histoEPOSLHCEtaToKCh[cent]                          = (TH1F*) directoryTheory[cent]->Get("histoEtaToKChEPOSLHC");
    }

    // *******************************************************************************************************
    // ************************** Loading interpolated spectra ***********************************************
    // *******************************************************************************************************
    TGraphAsymmErrors* statErrorCollectionPi0PP[11];
    TGraphAsymmErrors* systErrorCollectionPi0PP[11];
    TGraphAsymmErrors* systErrorUnCorrCollectionPi0PP[11];
    TGraphAsymmErrors* statErrorCollectionInvYieldPi0PP[11];
    TGraphAsymmErrors* systErrorCollectionInvYieldPi0PP[11];
    TGraphAsymmErrors* systErrorUnCorrCollectionInvYieldPi0PP[11];

    TGraphAsymmErrors* statErrorCollectionEtaPP[11];
    TGraphAsymmErrors* systErrorCollectionEtaPP[11];
    TGraphAsymmErrors* systErrorUnCorrCollectionEtaPP[11];
    TGraphAsymmErrors* statErrorCollectionInvYieldEtaPP[11];
    TGraphAsymmErrors* systErrorCollectionInvYieldEtaPP[11];
    TGraphAsymmErrors* systErrorUnCorrCollectionInvYieldEtaPP[11];
    TGraphAsymmErrors* graphRatioPi0IndCombFitStatPP[11];
    TGraphAsymmErrors* graphRatioPi0IndCombFitSysPP[11];
    TGraphAsymmErrors* graphRatioEtaIndCombFitStatPP[11];
    TGraphAsymmErrors* graphRatioEtaIndCombFitSysPP[11];
    TGraphAsymmErrors* graphPPCombPi0Stat                = NULL;
    TGraphAsymmErrors* graphPPCombPi0UncorrSys           = NULL;
    TGraphAsymmErrors* graphPPCombPi0FullSys             = NULL;
    TGraphAsymmErrors* graphPPCombPi0InterSys            = NULL;
    TGraphAsymmErrors* graphPPInvYieldCombPi0Stat        = NULL;
    TGraphAsymmErrors* graphPPInvYieldCombPi0StatW0XErr  = NULL;
    TGraphAsymmErrors* graphPPInvYieldCombPi0Sys         = NULL;
    TGraphAsymmErrors* graphPPCombEtaStat                = NULL;
    TGraphAsymmErrors* graphPPCombEtaUncorrSys           = NULL;
    TGraphAsymmErrors* graphPPCombEtaFullSys             = NULL;
    TGraphAsymmErrors* graphPPCombEtaInterSys            = NULL;
    TGraphAsymmErrors* graphPPInvYieldCombEtaStat        = NULL;
    TGraphAsymmErrors* graphPPInvYieldCombEtaStatW0XErr  = NULL;
    TGraphAsymmErrors* graphPPInvYieldCombEtaSys         = NULL;
    TGraphAsymmErrors* graphRatioPi0CombToFitStatPP      = NULL;
    TGraphAsymmErrors* graphRatioPi0CombToFitSysPP       = NULL;
    TGraphAsymmErrors* graphRatioEtaCombToFitStatPP      = NULL;
    TGraphAsymmErrors* graphRatioEtaCombToFitSysPP       = NULL;
    TF1* fitPPTCMInvYieldPi0                             = NULL;
    TF1* fitTCMInvXSetionPi0PP                           = NULL;
    TF1* fitTCMInvXSetionPi0PPOrg                        = NULL;
    TF1* fitTsallisInvXSectionPi0PP                      = NULL;
    TF1* fitPPTCMInvYieldEta                             = NULL;
    TF1* fitTCMInvXSetionEtaPP                           = NULL;
    TF1* fitTCMInvXSetionEtaPPOrg                        = NULL;
    TF1* fitTsallisInvXSectionEtaPP                      = NULL;

    TGraphAsymmErrors* graphPPCombEtaToPi0Stat           = NULL;
    TGraphAsymmErrors* graphPPCombEtaToPi0Sys            = NULL;
    TGraphAsymmErrors* graphPPCombEtaToPi0Tot            = NULL;
    TGraphAsymmErrors* graphPPCombEtaToPi0StatW0XErr     = NULL;

    Bool_t haveRefPPPi0[11]                              = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE } ;
    Bool_t haveRefPPEta[11]                              = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE } ;
    for (Int_t meth = 0; meth< 11; meth++){
        statErrorCollectionPi0PP[meth]                   = NULL;
        systErrorCollectionPi0PP[meth]                   = NULL;
        systErrorUnCorrCollectionPi0PP[meth]             = NULL;
        statErrorCollectionInvYieldPi0PP[meth]           = NULL;
        systErrorCollectionInvYieldPi0PP[meth]           = NULL;
        systErrorUnCorrCollectionInvYieldPi0PP[meth]     = NULL;
        statErrorCollectionEtaPP[meth]                   = NULL;
        systErrorCollectionEtaPP[meth]                   = NULL;
        systErrorUnCorrCollectionEtaPP[meth]             = NULL;
        statErrorCollectionInvYieldEtaPP[meth]           = NULL;
        systErrorCollectionInvYieldEtaPP[meth]           = NULL;
        systErrorUnCorrCollectionInvYieldEtaPP[meth]     = NULL;
        graphRatioPi0IndCombFitStatPP[meth]              = NULL;
        graphRatioPi0IndCombFitSysPP[meth]               = NULL;
        graphRatioEtaIndCombFitStatPP[meth]              = NULL;
        graphRatioEtaIndCombFitSysPP[meth]               = NULL;
    }

    TFile* fileReference                        = new TFile(fileNameReference.Data());

    for (Int_t meth = 0; meth< 11; meth++){
        TString currMethodRead                          = nameMeasGlobalLabel[meth].Data();
        TString currMethod                              = nameMeasGlobalLabel[meth].Data();
        statErrorCollectionPi0PP[meth]                 = (TGraphAsymmErrors*)fileReference->Get(Form("Pi05TeVPbPbRef/graphInvCrossSectionPi0%s5TeVStatErr_yShifted",
                                                                                                                currMethodRead.Data()));
        systErrorCollectionPi0PP[meth]                 = (TGraphAsymmErrors*)fileReference->Get(Form("Pi05TeVPbPbRef/graphInvCrossSectionPi0%s5TeVSysErr_yShifted",
                                                                                                                currMethodRead.Data()));
        systErrorUnCorrCollectionPi0PP[meth]           = (TGraphAsymmErrors*)fileReference->Get(Form("Pi05TeVPbPbRef/graphInvCrossSectionPi0%s5TeVSysErr_yShifted",
                                                                                                         currMethodRead.Data()));
        if (statErrorCollectionPi0PP[meth] && systErrorCollectionPi0PP[meth] && systErrorUnCorrCollectionPi0PP[meth] ){
            haveRefPPPi0[meth]                              = kTRUE;
            statErrorCollectionInvYieldPi0PP[meth]          = (TGraphAsymmErrors*)statErrorCollectionPi0PP[meth]->Clone(Form("graphInvYieldPi0%s5TeVStatErr",
                                                                                                                    currMethod.Data()));
            systErrorCollectionInvYieldPi0PP[meth]          = (TGraphAsymmErrors*)systErrorCollectionPi0PP[meth]->Clone(Form("graphInvYieldPi0%s5TeVSysErr",
                                                                                                                    currMethod.Data()));
            systErrorUnCorrCollectionInvYieldPi0PP[meth]    = (TGraphAsymmErrors*)systErrorUnCorrCollectionPi0PP[meth]->Clone(Form("graphInvYieldPi0%s5TeVUnCorrSysErr",
                                                                                                                    currMethod.Data()));

            statErrorCollectionInvYieldPi0PP[meth]          = ScaleGraph(statErrorCollectionInvYieldPi0PP[meth],1/(xSection5TeV*recalcBarn));
            systErrorCollectionInvYieldPi0PP[meth]          = ScaleGraph(systErrorCollectionInvYieldPi0PP[meth],1/(xSection5TeV*recalcBarn));
            systErrorUnCorrCollectionInvYieldPi0PP[meth]    = ScaleGraph(systErrorUnCorrCollectionInvYieldPi0PP[meth],1/(xSection5TeV*recalcBarn));

        }
        if (haveRefPPPi0[meth])
            cout << "found pi0 pp reference for " << currMethod.Data() << endl;


        if (meth == 1){
            currMethodRead                              = "Comb";
        }
        statErrorCollectionEtaPP[meth]                  = (TGraphAsymmErrors*)fileReference->Get(Form("Eta5TeVPbPbRef/graphInvCrossSectionEta%s5TeVStatErr_yShifted",
                                                                                                                currMethodRead.Data()));
        systErrorCollectionEtaPP[meth]                  = (TGraphAsymmErrors*)fileReference->Get(Form("Eta5TeVPbPbRef/graphInvCrossSectionEta%s5TeVSysErr_yShifted",
                                                                                                                currMethodRead.Data()));
        systErrorUnCorrCollectionEtaPP[meth]            = (TGraphAsymmErrors*)fileReference->Get(Form("Eta5TeVPbPbRef/graphInvCrossSectionEta%s5TeVSysErr_yShifted",
                                                                                                         currMethodRead.Data()));
        if (statErrorCollectionEtaPP[meth] && systErrorCollectionEtaPP[meth] && systErrorUnCorrCollectionEtaPP[meth] ){
            haveRefPPEta[meth]                          = kTRUE;
            statErrorCollectionInvYieldEtaPP[meth]           = (TGraphAsymmErrors*)statErrorCollectionEtaPP[meth]->Clone(Form("graphInvYieldEta%s5TeVStatErr",
                                                                                                                        currMethod.Data()));
            systErrorCollectionInvYieldEtaPP[meth]           = (TGraphAsymmErrors*)systErrorCollectionEtaPP[meth]->Clone(Form("graphInvYieldEta%s5TeVSysErr",
                                                                                                                        currMethod.Data()));
            systErrorUnCorrCollectionInvYieldEtaPP[meth]     = (TGraphAsymmErrors*)systErrorUnCorrCollectionEtaPP[meth]->Clone(Form("graphInvYieldEta%s5TeVUnCorrSysErr",
                                                                                                                        currMethod.Data()));

            statErrorCollectionInvYieldEtaPP[meth]           = ScaleGraph(statErrorCollectionInvYieldEtaPP[meth],1/(xSection5TeV*recalcBarn));
            systErrorCollectionInvYieldEtaPP[meth]           = ScaleGraph(systErrorCollectionInvYieldEtaPP[meth],1/(xSection5TeV*recalcBarn));
            systErrorUnCorrCollectionInvYieldEtaPP[meth]     = ScaleGraph(systErrorUnCorrCollectionInvYieldEtaPP[meth],1/(xSection5TeV*recalcBarn));
        }
        if (haveRefPPEta[meth])
            cout << "found eta pp reference for " << currMethod.Data() << endl;
    }

    graphPPCombPi0Stat                      = (TGraphAsymmErrors*)fileReference->Get(Form("Pi05TeV/graphInvCrossSectionPi0Comb5TeVStatErr"));
    graphPPCombPi0UncorrSys                 = (TGraphAsymmErrors*)fileReference->Get(Form("Pi05TeV/graphInvCrossSectionPi0Comb5TeVSysErr"));
    graphPPCombPi0FullSys                   = (TGraphAsymmErrors*)fileReference->Get(Form("Pi05TeV/graphInvCrossSectionPi0Comb5TeVSysErr"));
    fitTCMInvXSetionPi0PPOrg                = (TF1*)fileReference->Get("Pi05TeV/TwoComponentModelFitPi0");

    graphPPInvYieldCombPi0Stat              = (TGraphAsymmErrors*)graphPPCombPi0Stat->Clone(Form("graphInvYieldPi0Comb5TeVStatErr"));
    graphPPInvYieldCombPi0Sys               = (TGraphAsymmErrors*)graphPPCombPi0FullSys->Clone(Form("graphInvYieldPi0Comb5TeVSysErr"));
    graphPPInvYieldCombPi0Stat              = ScaleGraph(graphPPCombPi0Stat,1/(xSection5TeV*recalcBarn));
    graphPPInvYieldCombPi0Sys               = ScaleGraph(graphPPInvYieldCombPi0Sys,1/(xSection5TeV*recalcBarn));
    graphPPInvYieldCombPi0StatW0XErr        = (TGraphAsymmErrors*)graphPPInvYieldCombPi0Stat->Clone(Form("graphInvYieldPi0Comb5TeVStatErrW0XErr"));
    ProduceGraphAsymmWithoutXErrors(graphPPInvYieldCombPi0StatW0XErr);

    graphPPCombEtaStat                      = (TGraphAsymmErrors*)fileReference->Get(Form("Eta5TeV/graphInvCrossSectionEtaComb5TeVStatErr"));
    graphPPCombEtaUncorrSys                 = (TGraphAsymmErrors*)fileReference->Get(Form("Eta5TeV/graphInvCrossSectionEtaComb5TeVSysErr"));
    graphPPCombEtaFullSys                   = (TGraphAsymmErrors*)fileReference->Get(Form("Eta5TeV/graphInvCrossSectionEtaComb5TeVSysErr"));
    fitTCMInvXSetionEtaPPOrg                = (TF1*)fileReference->Get("Eta5TeV/TwoComponentModelFitEta");
    graphPPCombEtaToPi0Stat                 = (TGraphAsymmErrors*)fileReference->Get(Form("Eta5TeV/graphRatioEtaToPi0Comb5TeVStatErr"));
    graphPPCombEtaToPi0Sys                  = (TGraphAsymmErrors*)fileReference->Get(Form("Eta5TeV/graphRatioEtaToPi0Comb5TeVSysErr"));
    graphPPCombEtaToPi0Tot                  = (TGraphAsymmErrors*)fileReference->Get(Form("Eta5TeV/graphRatioEtaToPi0Comb5TeVTotErr"));
    graphPPCombEtaToPi0StatW0XErr           = (TGraphAsymmErrors*)graphPPCombEtaToPi0Stat->Clone(Form("graphEtaToPi0PPComb5TeVStatErrW0XErr"));
    ProduceGraphAsymmWithoutXErrors(graphPPCombEtaToPi0StatW0XErr);

    graphPPInvYieldCombEtaStat              = (TGraphAsymmErrors*)graphPPCombEtaStat->Clone(Form("graphInvYieldEtaComb5TeVStatErr"));
    graphPPInvYieldCombEtaSys               = (TGraphAsymmErrors*)graphPPCombEtaFullSys->Clone(Form("graphInvYieldEtaComb5TeVSysErr"));
    graphPPInvYieldCombEtaStat              = ScaleGraph(graphPPCombEtaStat,1/(xSection5TeV*recalcBarn));
    graphPPInvYieldCombEtaSys               = ScaleGraph(graphPPInvYieldCombEtaSys,1/(xSection5TeV*recalcBarn));
    graphPPInvYieldCombEtaStatW0XErr        = (TGraphAsymmErrors*)graphPPInvYieldCombEtaStat->Clone(Form("graphInvYieldEtaComb5TeVStatErrW0XErr"));
    ProduceGraphAsymmWithoutXErrors(graphPPInvYieldCombEtaStatW0XErr);


    TFile* fileOtherParticle                           = new TFile(fileNameOtherParticleInput.Data());
    TDirectory* directoryOtherPart[9]                  = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histChPiSpecStat[9]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histChPiSpecSyst[9]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChPiSpecStatErr[9]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChPiSpecSystErr[9]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChPiSpecRelStatErr[9]      = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChPiSpecRelSystErr[9]      = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphChPiSpecRelTotErr[9]       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
//     TGraphAsymmErrors* graphDRAAStat[9]                = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
//     TGraphAsymmErrors* graphDRAASys[9]                 = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
//     TGraphAsymmErrors* graphDRAATot[9]                 = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
//     TGraphAsymmErrors* graphChHadRAAStat[9]            = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
//     TGraphAsymmErrors* graphChHadRAASys[9]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
//     TGraphAsymmErrors* graphChHadRAATot[9]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
//     TGraphAsymmErrors* graphChKToPiStat[9]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
//     TGraphAsymmErrors* graphChKToPiSys[9]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
//     TGraphAsymmErrors* graphChKToPiTot[9]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    for (Int_t cent = 0; cent < maxCent; cent++){
        TString folderNameOtherPart                     = Form("%sPbPb_5.02TeV",centArray[cent].Data());
        directoryOtherPart[cent]                               = (TDirectory*)fileOtherParticle->Get(folderNameOtherPart.Data());
        cout << "trying to find " << folderNameOtherPart.Data() << endl;
        if (directoryOtherPart[cent]){
            cout << "found it" << endl;
            histChPiSpecStat[cent]                              = (TH1D*)directoryOtherPart[cent]->Get("histoChargedPionSpecStat");
            histChPiSpecSyst[cent]                              = (TH1D*)directoryOtherPart[cent]->Get("histoChargedPionSpecSyst");
            cout << histChPiSpecStat[cent] << "\t" << histChPiSpecSyst[cent] << endl;
            if (histChPiSpecStat[cent] && histChPiSpecSyst[cent]){
                graphChPiSpecSystErr[cent]                          = new TGraphAsymmErrors(histChPiSpecSyst[cent]);
                graphChPiSpecStatErr[cent]                          = new TGraphAsymmErrors(histChPiSpecStat[cent]);
                graphChPiSpecRelTotErr[cent]                        = AddErrorsOfGraphsQuadratically(graphChPiSpecStatErr[cent],graphChPiSpecSystErr[cent]);
                graphChPiSpecRelStatErr[cent]                       = CalculateRelErrUpAsymmGraph( graphChPiSpecStatErr[cent], Form("relativeStatErrorPiCh_%s",centArray[cent].Data()));
                graphChPiSpecRelSystErr[cent]                       = CalculateRelErrUpAsymmGraph( graphChPiSpecSystErr[cent], Form("relativeSysErrorPiCh_%s",centArray[cent].Data()));
                graphChPiSpecRelTotErr[cent]                        = CalculateRelErrUpAsymmGraph( graphChPiSpecRelTotErr[cent], Form("relativeTotErrorPiCh_%s",centArray[cent].Data()));
            }

//             graphDRAAStat[cent]                                 = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("DMesonQPbPbStatW0XErr");
//             graphDRAASys[cent]                                  = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("DMesonQPbPbSys");
//             graphDRAATot[cent]                                  = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("DMesonQPbPbTot");
//             graphChHadRAAStat[cent]                             = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("ChargedHadQPbPbStatW0XErr");
//             graphChHadRAASys[cent]                              = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("ChargedHadQPbPbSys");
//             graphChHadRAATot[cent]                              = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("ChargedHadQPbPbTot");
//             graphChKToPiStat[cent]                              = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("graphChargedKaonToPionPubStatW0XErr");
//             graphChKToPiSys[cent]                               = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("graphChargedKaonToPionPubSys");
//             graphChKToPiTot[cent]                               = (TGraphAsymmErrors*)directoryOtherPart[cent]->Get("graphChargedKaonToPionPubTot");
        }
    }

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
        TF1* fitTCMDecomposedLPi0PP                             = FitObject("tcmlow",Form("twoCompModelPi0_DecL_PP"), "Pi0", NULL, 0.3, 2.);
        TF1* fitTCMDecomposedHPi0PP                             = FitObject("tcmhigh",Form("twoCompModelPi0_DecH_PP"), "Pi0", NULL, 4, 50.);
        fitTCMDecomposedLPi0PP->SetParameters(graphPPCombPi0Stat->GetY()[2],0.3);
        graphPPCombPi0Stat->Fit(fitTCMDecomposedLPi0PP,"QNRMEX0+","",0.3,0.8);
        graphPPCombPi0Stat->Fit(fitTCMDecomposedHPi0PP,"QNRMEX0+","",3, 30.);
        fitTCMDecomposedHPi0PP->SetParameters(graphPPCombPi0Stat->GetY()[2],0.8, 2);

        cout << WriteParameterToFile(fitTCMDecomposedLPi0PP)<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedLPi0PP)<< endl;
        cout << WriteParameterToFile(fitTCMDecomposedHPi0PP)<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedHPi0PP)<< endl;

        Double_t paramTCMPi0NewPP[5]                            = { fitTCMInvXSetionPi0PPOrg->GetParameter(0),fitTCMInvXSetionPi0PPOrg->GetParameter(1),
                                                                    fitTCMInvXSetionPi0PPOrg->GetParameter(2),fitTCMInvXSetionPi0PPOrg->GetParameter(3),fitTCMInvXSetionPi0PPOrg->GetParameter(4)};

        //Two component model from Bylinkin
        fitTCMInvXSetionPi0PP                            = FitObject("tcm",Form("fitTCMInvXSetionPi0PP5TeV"),"Pi0",graphPPCombPi0Stat, 0.3, 30. ,paramTCMPi0NewPP,"QNRMEX0+","", kFALSE);
        cout << "fitting pp pi0 spectra" << endl;
        cout << WriteParameterToFile(fitTCMInvXSetionPi0PP)<< endl;
        fitTsallisInvXSectionPi0PP                       = FitObject("l",Form("fitTsallisInvXSectionPi0PP5TeV"),"Pi0",graphPPCombPi0Stat, 0.3, 30. ,
                                                                              NULL,"QNRMEX0+","", kFALSE);
        cout << "fitting pp pi0 spectra" << endl;
        cout << WriteParameterToFile(fitTsallisInvXSectionPi0PP)<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTsallisInvXSectionPi0PP)<< endl;

        fileFitsOutput <<  WriteParameterToFile(fitTCMInvXSetionPi0PP)<< endl;
        fitPPTCMInvYieldPi0                              = FitObject("tcm",Form("fitTCMInvYieldPi0PP5TeV"),"Pi0",NULL, 0.3, 30. ,NULL,"QNRMEX0+","", kFALSE);
        cout << "fitting pp pi0 spectra" << endl;
        fitPPTCMInvYieldPi0->SetParameter(0,fitTCMInvXSetionPi0PP->GetParameter(0)/(xSection5TeV*recalcBarn));
        fitPPTCMInvYieldPi0->SetParameter(1,fitTCMInvXSetionPi0PP->GetParameter(1));
        fitPPTCMInvYieldPi0->SetParameter(2,fitTCMInvXSetionPi0PP->GetParameter(2)/(xSection5TeV*recalcBarn));
        fitPPTCMInvYieldPi0->SetParameter(3,fitTCMInvXSetionPi0PP->GetParameter(3));
        fitPPTCMInvYieldPi0->SetParameter(4,fitTCMInvXSetionPi0PP->GetParameter(4));
        graphPPInvYieldCombPi0Stat->Fit(fitPPTCMInvYieldPi0,"QNRMEX0+","", 0.3, 30.);
        cout << WriteParameterToFile(fitPPTCMInvYieldPi0)<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitPPTCMInvYieldPi0)<< endl;
        fitPPTCMInvYieldPi0->SetRange(0.3, 30.);

        // *************************************************************************************************************
        // Shift graphs in Y direction as well if desired
        // *************************************************************************************************************

        graphRatioPi0CombToFitStatPP                         = (TGraphAsymmErrors*)graphPPCombPi0Stat->Clone();
        graphRatioPi0CombToFitStatPP                         = CalculateGraphErrRatioToFit(graphRatioPi0CombToFitStatPP, fitTCMInvXSetionPi0PP);
        graphRatioPi0CombToFitSysPP                          = (TGraphAsymmErrors*)graphPPCombPi0FullSys->Clone();
        graphRatioPi0CombToFitSysPP                          = CalculateGraphErrRatioToFit(graphRatioPi0CombToFitSysPP, fitTCMInvXSetionPi0PP);

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
        canvasRatioToCombFitPP->cd();
        histo2DPi0RatioToCombFitPP->Draw("copy");

        ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombToFitStatPP);

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombToFitSysPP, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioPi0CombToFitSysPP->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombToFitStatPP, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioPi0CombToFitStatPP->Draw("p,same,z");

        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1., 1.,1, kGray+2);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.1, 1.1,1, kGray, 7);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.9, 0.9,1, kGray, 7);

        labelRatioToFitEnergyPP->Draw();
        labelRatioToFitALICEPP->Draw();
        labelRatioToFitPi0PP->Draw();

        canvasRatioToCombFitPP->SaveAs(Form("%s/Pi0_RatioOfCombToCombFit_PP5TeV.%s",outputDirSupportPaper.Data(), suffix.Data()));

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

        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1., 1.,1, kGray+2);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.1, 1.1,1, kGray, 7);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.9, 0.9,1, kGray, 7);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.2, 1.2,1, kGray, 9);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.8, 0.8,1, kGray, 9);

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

        if(graphRatioPi0IndCombFitSysPP[5]) {
            TMarker* markerDalitzPi0OnlyRatioPi0PP      = CreateMarkerFromGraph(graphRatioPi0IndCombFitSysPP[5], columnsLegendOnlyPi0RatioPP[3] ,rowsLegendOnlyPi0RatioPP[3],1);
            markerDalitzPi0OnlyRatioPi0PP->DrawMarker(columnsLegendOnlyPi0RatioPPAbs[4] ,rowsLegendOnlyPi0RatioPPAbs[3]);
            TBox* boxDalitzPi0OnlyRatioPi0PP             = CreateBoxFromGraph(graphRatioPi0IndCombFitSysPP[5], columnsLegendOnlyPi0RatioPPAbs[5]-0.5*lengthBoxPP , rowsLegendOnlyPi0RatioPPAbs[3]- heightBoxPP,
                                                                              columnsLegendOnlyPi0RatioPPAbs[5]+ 17.2*lengthBoxPP, rowsLegendOnlyPi0RatioPPAbs[3]+ heightBoxPP);
            boxDalitzPi0OnlyRatioPi0PP->Draw("l");
        }

        canvasRatioToCombFitPP->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFit_PP5TeV.%s",outputDirSupportPaper.Data(),suffix.Data()));

        // *******************************************************************************************************
        // ************************** Comparison of different pp interpolated eta spectra ************************
        // *******************************************************************************************************
        // fitting spectrum with intial parameters
        // Two component model fit from Bylinkin
        TF1* fitTCMDecomposedLEtaPP                             = FitObject("tcmlow",Form("twoCompModelEta_DecL_PP"), "Eta", NULL, 0.8, 2.);
        TF1* fitTCMDecomposedHEtaPP                             = FitObject("tcmhigh",Form("twoCompModelEta_DecH_PP"), "Eta", NULL, 4, 50.);
        fitTCMDecomposedLEtaPP->SetParameters(graphPPCombEtaStat->GetY()[2],0.3);
        graphPPCombEtaStat->Fit(fitTCMDecomposedLEtaPP,"QNRMEX0+","",0.8,2);
        graphPPCombEtaStat->Fit(fitTCMDecomposedHEtaPP,"QNRMEX0+","",3, 20.);
        fitTCMDecomposedHEtaPP->SetParameters(graphPPCombEtaStat->GetY()[2],0.8, 2);

        cout << WriteParameterToFile(fitTCMDecomposedLEtaPP)<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedLEtaPP)<< endl;
        cout << WriteParameterToFile(fitTCMDecomposedHEtaPP)<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedHEtaPP)<< endl;


        Double_t paramTCMEtaNewPP[5]                            = { fitTCMInvXSetionEtaPPOrg->GetParameter(0),fitTCMInvXSetionEtaPPOrg->GetParameter(1),
                                                                    fitTCMInvXSetionEtaPPOrg->GetParameter(2),fitTCMInvXSetionEtaPPOrg->GetParameter(3),fitTCMInvXSetionEtaPPOrg->GetParameter(4)};

        //Two component model from Bylinkin
        fitTCMInvXSetionEtaPP                            = FitObject("tcm",Form("fitTCMInvXSetionEtaPP5TeV"),"Eta",graphPPCombEtaStat, 0.5, 30. ,paramTCMEtaNewPP,"QNRMEX0+","", kFALSE);
        cout << "fitting pp eta spectra" << endl;
        cout << WriteParameterToFile(fitTCMInvXSetionEtaPP)<< endl;
        fitTsallisInvXSectionEtaPP                       = FitObject("l",Form("fitTsallisInvXSectionEtaPP5TeV"),"Eta",graphPPCombEtaStat, 0.5, 30. ,
                                                                            NULL,"QNRMEX0+","", kFALSE);
        cout << "fitting pp eta spectra" << endl;
        cout << WriteParameterToFile(fitTsallisInvXSectionEtaPP)<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTsallisInvXSectionEtaPP)<< endl;

        fileFitsOutput <<  WriteParameterToFile(fitTCMInvXSetionEtaPP)<< endl;
        fitPPTCMInvYieldEta                              = FitObject("tcm",Form("fitTCMInvYieldEtaPP5TeV"),"Eta",NULL, 0.5, 30. ,NULL,"QNRMEX0+","", kFALSE);
        cout << "fitting pp eta spectra" << endl;
        fitPPTCMInvYieldEta->SetParameter(0,fitTCMInvXSetionEtaPP->GetParameter(0)/(xSection5TeV*recalcBarn));
        fitPPTCMInvYieldEta->SetParameter(1,fitTCMInvXSetionEtaPP->GetParameter(1));
        fitPPTCMInvYieldEta->SetParameter(2,fitTCMInvXSetionEtaPP->GetParameter(2)/(xSection5TeV*recalcBarn));
        fitPPTCMInvYieldEta->SetParameter(3,fitTCMInvXSetionEtaPP->GetParameter(3));
        fitPPTCMInvYieldEta->SetParameter(4,fitTCMInvXSetionEtaPP->GetParameter(4));
        graphPPInvYieldCombEtaStat->Fit(fitPPTCMInvYieldEta,"QNRMEX0+","", 0.5, 30.);
        cout << WriteParameterToFile(fitPPTCMInvYieldEta)<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitPPTCMInvYieldEta)<< endl;
        fitPPTCMInvYieldEta->SetRange(0.5, 30.);

        // *************************************************************************************************************
        // Shift graphs in Y direction as well if desired
        // *************************************************************************************************************

        graphRatioEtaCombToFitStatPP                         = (TGraphAsymmErrors*)graphPPCombEtaStat->Clone();
        graphRatioEtaCombToFitStatPP                         = CalculateGraphErrRatioToFit(graphRatioEtaCombToFitStatPP, fitTCMInvXSetionEtaPP);
        graphRatioEtaCombToFitSysPP                          = (TGraphAsymmErrors*)graphPPCombEtaFullSys->Clone();
        graphRatioEtaCombToFitSysPP                          = CalculateGraphErrRatioToFit(graphRatioEtaCombToFitSysPP, fitTCMInvXSetionEtaPP);

        for (Int_t meth = 0; meth< 11; meth++){
            if (statErrorCollectionEtaPP[meth]){
                graphRatioEtaIndCombFitStatPP[meth]                = (TGraphAsymmErrors*)statErrorCollectionEtaPP[meth]->Clone(Form("RatioEta%sToCombFitStatPP", nameMeasGlobalLabel[meth].Data()));
                graphRatioEtaIndCombFitStatPP[meth]                = CalculateGraphErrRatioToFit(graphRatioEtaIndCombFitStatPP[meth], fitTCMInvXSetionEtaPP);
            }
            if (systErrorCollectionEtaPP[meth]){
                graphRatioEtaIndCombFitSysPP[meth]                 = (TGraphAsymmErrors*)systErrorCollectionEtaPP[meth]->Clone(Form("RatioEta%sToCombFitSystPP", nameMeasGlobalLabel[meth].Data()));
                graphRatioEtaIndCombFitSysPP[meth]                 = CalculateGraphErrRatioToFit(graphRatioEtaIndCombFitSysPP[meth], fitTCMInvXSetionEtaPP);
            }
        }

        // **********************************************************************************************************************
        // ******************************************* Plot Ratio of Comb to Fit for pp *****************************************
        // **********************************************************************************************************************
        canvasRatioToCombFitPP->cd();
        histo2DEtaRatioToCombFitPP->Draw("copy");

        ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombToFitStatPP);

        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombToFitSysPP, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioEtaCombToFitSysPP->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombToFitStatPP, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioEtaCombToFitStatPP->Draw("p,same,z");

        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1., 1.,1, kGray+2);
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1.1, 1.1,1, kGray, 7);
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.9, 0.9,1, kGray, 7);

        labelRatioToFitEnergyPP->Draw();
        labelRatioToFitALICEPP->Draw();
        labelRatioToFitEtaPP->Draw();

        canvasRatioToCombFitPP->SaveAs(Form("%s/Eta_RatioOfCombToCombFit_PP5TeV.%s",outputDirSupportPaper.Data(), suffix.Data()));

        // **********************************************************************************************************************
        // *******************************************Plot Ratio of Individual meas to Fit ******************************************
        // **********************************************************************************************************************

        canvasRatioToCombFitPP->cd();
        histo2DEtaRatioToCombFitPP->Draw("copy");

        for (Int_t meth = 10; meth > -1 ; meth--){
            if (graphRatioEtaIndCombFitSysPP[meth]){
                DrawGammaSetMarkerTGraphAsym(graphRatioEtaIndCombFitSysPP[meth], markerStyleDetMC[meth] ,markerSizeDetMC[meth]*0.5, colorDet[meth], colorDet[meth], widthLinesBoxes, kTRUE);
                graphRatioEtaIndCombFitSysPP[meth]->Draw("E2same");
            }
            if (graphRatioEtaIndCombFitStatPP[meth]){
                ProduceGraphAsymmWithoutXErrors(graphRatioEtaIndCombFitStatPP[meth]);
                DrawGammaSetMarkerTGraphAsym(graphRatioEtaIndCombFitStatPP[meth], markerStyleDetMC[meth] ,markerSizeDetMC[meth]*0.5, colorDet[meth], colorDet[meth]);
                graphRatioEtaIndCombFitStatPP[meth]->Draw("p,same,z");
            }
        }
        if (graphRatioEtaIndCombFitStatPP[4])graphRatioEtaIndCombFitStatPP[4]->Draw("p,same,z");

        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1., 1.,1, kGray+2);
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1.1, 1.1,1, kGray, 7);
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.9, 0.9,1, kGray, 7);
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1.2, 1.2,1, kGray, 9);
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.8, 0.8,1, kGray, 9);

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

        TMarker* markerPCMEtaOnlyRatioEtaPP             = CreateMarkerFromGraph(graphRatioEtaIndCombFitSysPP[0],columnsLegendOnlyEtaRatioPP[1] ,rowsLegendOnlyEtaRatioPP[1],1);
        markerPCMEtaOnlyRatioEtaPP->DrawMarker(columnsLegendOnlyEtaRatioPPAbs[1] ,rowsLegendOnlyEtaRatioPPAbs[1]);
        TMarker* markerPHOSEtaOnlyRatioEtaPP            = CreateMarkerFromGraph(graphRatioEtaIndCombFitSysPP[1], columnsLegendOnlyEtaRatioPP[1] ,rowsLegendOnlyEtaRatioPP[2],1);
        markerPHOSEtaOnlyRatioEtaPP->DrawMarker(columnsLegendOnlyEtaRatioPPAbs[1] ,rowsLegendOnlyEtaRatioPPAbs[2]);
        TMarker* markerEMCALEtaOnlyRatioEtaPP           = CreateMarkerFromGraph(graphRatioEtaIndCombFitSysPP[2], columnsLegendOnlyEtaRatioPP[1] ,rowsLegendOnlyEtaRatioPP[3],1);
        markerEMCALEtaOnlyRatioEtaPP->DrawMarker(columnsLegendOnlyEtaRatioPPAbs[1] ,rowsLegendOnlyEtaRatioPPAbs[3]);
        TMarker* markerPCMEMCALEtaOnlyRatioEtaPP        = CreateMarkerFromGraph(graphRatioEtaIndCombFitSysPP[4], columnsLegendOnlyEtaRatioPP[3] ,rowsLegendOnlyEtaRatioPP[1],1);
        markerPCMEMCALEtaOnlyRatioEtaPP->DrawMarker(columnsLegendOnlyEtaRatioPPAbs[4] ,rowsLegendOnlyEtaRatioPPAbs[1]);
        TMarker* markerPCMPHOSEtaOnlyRatioEtaPP         = CreateMarkerFromGraph(graphRatioEtaIndCombFitSysPP[3], columnsLegendOnlyEtaRatioPP[3] ,rowsLegendOnlyEtaRatioPP[2],1);
        markerPCMPHOSEtaOnlyRatioEtaPP->DrawMarker(columnsLegendOnlyEtaRatioPPAbs[4] ,rowsLegendOnlyEtaRatioPPAbs[2]);

        TBox* boxPCMEtaOnlyRatioEtaPP                 = CreateBoxFromGraph(graphRatioEtaIndCombFitSysPP[0], columnsLegendOnlyEtaRatioPPAbs[2]-0.1*lengthBoxPPEta , rowsLegendOnlyEtaRatioPPAbs[1]- heightBoxPP,
                                                                            columnsLegendOnlyEtaRatioPPAbs[2]+ 1.7*lengthBoxPPEta, rowsLegendOnlyEtaRatioPPAbs[1]+ heightBoxPP);
        boxPCMEtaOnlyRatioEtaPP->Draw("l");
        TBox* boxPHOSEtaOnlyRatioEtaPP                = CreateBoxFromGraph(graphRatioEtaIndCombFitSysPP[1], columnsLegendOnlyEtaRatioPPAbs[2]-0.1*lengthBoxPPEta , rowsLegendOnlyEtaRatioPPAbs[2]- heightBoxPP,
                                                                            columnsLegendOnlyEtaRatioPPAbs[2]+ 1.7*lengthBoxPPEta, rowsLegendOnlyEtaRatioPPAbs[2]+ heightBoxPP);
        boxPHOSEtaOnlyRatioEtaPP->Draw("l");
        TBox* boxEMCALEtaOnlyRatioEtaPP               = CreateBoxFromGraph(graphRatioEtaIndCombFitSysPP[2], columnsLegendOnlyEtaRatioPPAbs[2]-0.1*lengthBoxPPEta , rowsLegendOnlyEtaRatioPPAbs[3]- heightBoxPP,
                                                                            columnsLegendOnlyEtaRatioPPAbs[2]+ 1.7*lengthBoxPPEta, rowsLegendOnlyEtaRatioPPAbs[3]+ heightBoxPP);
        boxEMCALEtaOnlyRatioEtaPP->Draw("l");
        TBox* boxPCMEMCALEtaOnlyRatioEtaPP            = CreateBoxFromGraph(graphRatioEtaIndCombFitSysPP[4], columnsLegendOnlyEtaRatioPPAbs[5]-0.5*lengthBoxPPEta , rowsLegendOnlyEtaRatioPPAbs[1]- heightBoxPP,
                                                                            columnsLegendOnlyEtaRatioPPAbs[5]+ 17.*lengthBoxPPEta, rowsLegendOnlyEtaRatioPPAbs[1]+ heightBoxPP);
        boxPCMEMCALEtaOnlyRatioEtaPP->Draw("l");
        TBox* boxPCMPHOSEtaOnlyRatioEtaPP             = CreateBoxFromGraph(graphRatioEtaIndCombFitSysPP[3], columnsLegendOnlyEtaRatioPPAbs[5]-0.5*lengthBoxPPEta , rowsLegendOnlyEtaRatioPPAbs[2]- heightBoxPP,
                                                                            columnsLegendOnlyEtaRatioPPAbs[5]+ 17.*lengthBoxPPEta, rowsLegendOnlyEtaRatioPPAbs[2]+ heightBoxPP);
        boxPCMPHOSEtaOnlyRatioEtaPP->Draw("l");

        canvasRatioToCombFitPP->SaveAs(Form("%s/Eta_RatioOfIndividualMeasToCombFit_PP5TeV.%s",outputDirSupportPaper.Data(), suffix.Data()));
    }

    // *******************************************************************************************************
    // ************************** Combination of different pi0 measurements **********************************
    // *******************************************************************************************************
    // REMARKS:
    //     - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - extendable to other energies
    //     - offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented

    TH1D* statErrorCollectionPi0[9][11];
    TGraphAsymmErrors* statErrorGraphCollectionPi0[9][11];
    TGraphAsymmErrors* sysErrorCollectionPi0[9][11];
    TH1D* statErrorCollectionOnlyIndepMethodsPi0[9][11];
    TGraphAsymmErrors* statErrorGraphCollectionOnlyIndepMethodsPi0[9][11];
    TGraphAsymmErrors* sysErrorCollectionOnlyIndepMethodsPi0[9][11];

    for (Int_t cent = 0; cent < maxCent; cent++){
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
            if (meth < 3){
                statErrorCollectionOnlyIndepMethodsPi0[cent][meth]          = NULL;
                if (histoPi0InvYieldStat[cent][meth]) statErrorCollectionOnlyIndepMethodsPi0[cent][meth]      = (TH1D*)histoPi0InvYieldStat[cent][meth]->Clone(Form("statErrIndepMeth%i_%sPi0",cent,nameMeasGlobalLabel[meth].Data()));
                statErrorGraphCollectionOnlyIndepMethodsPi0[cent][meth]     = NULL;
                if (graphPi0InvYieldStat[cent][meth]) statErrorGraphCollectionOnlyIndepMethodsPi0[cent][meth] = (TGraphAsymmErrors*)graphPi0InvYieldStat[cent][meth]->Clone(Form("statErrGraphIndepMeth%i_%sPi0", cent,
                    nameMeasGlobalLabel[meth].Data()));
                if (statErrorGraphCollectionOnlyIndepMethodsPi0[cent][meth]) RemoveZerosAtBeginningAndEndFromGraph(statErrorGraphCollectionOnlyIndepMethodsPi0[cent][meth]);
                // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
                // for systematic error from respective method
                sysErrorCollectionOnlyIndepMethodsPi0[cent][meth]           = NULL;
                if (graphPi0InvYieldSys[cent][meth]) sysErrorCollectionOnlyIndepMethodsPi0[cent][meth]        = (TGraphAsymmErrors*)graphPi0InvYieldSys[cent][meth]->Clone(Form("sysErrIndepMeth%i_%sPi0", cent,
                    nameMeasGlobalLabel[meth].Data()));
            } else {
                statErrorCollectionOnlyIndepMethodsPi0[cent][meth]          = NULL;
                statErrorGraphCollectionOnlyIndepMethodsPi0[cent][meth]     = NULL;
                sysErrorCollectionOnlyIndepMethodsPi0[cent][meth]           = NULL;
            }
        }
    }

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Pi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsPi0[9][11]             = {
                                            { 0,  1,  0,  0,  0,  3,  0,  0,  0,  0,  0 },  // 0-10
                                            { 0,  1,  0,  0,  0,  3,  0,  0,  0,  0,  0 },  // 10-20
                                            { 0,  1,  0,  0,  0,  3,  0,  0,  0,  0,  0 },  // 20-40
                                            { 0,  1,  0,  0,  0,  3,  0,  0,  0,  0,  0 },  // 40-60
                                            { 0,  1,  0,  0,  0,  3,  0,  0,  0,  0,  0 },  // 60-80
                                            { 0,  6,  0,  0,  0,  3,  0,  0,  0,  0,  0 },  // 0-5%
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 5-10%
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 },  // 0-20
                                            { 0,  5,  1,  1,  1,  3,  0,  0,  0,  0,  0 }   // 20-50
                                           };
    Int_t offSetsPi0Sys[9][11]          = {
                                            { 1,  1,  23,  3,  10,  0,  0,  0,  0,  21, 0 },  // 0-10
                                            { 1,  1,  19,  3,  10,  0,  0,  0,  0,  21, 0 },  // 10-20
                                            { 1,  1,  6,  3,  4,  0,  0,  0,  0,  21, 0 },  // 20-40
                                            { 1,  1,  6,  3,  4,  0,  0,  0,  0,  21, 0 },  // 40-60
                                            { 1,  1,  6,  3,  4,  0,  0,  0,  0,  21, 0 },  // 60-80
                                            { 1,  7,  9,  3,  6,  4,  0,  0,  0,  21, 0 },  // 0-5
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 5-10
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 },  // 0-20
                                            { 1,  6,  8,  2,  5,  0,  0,  0,  0,  21, 0 }   // 20-50
                                           };
    Int_t offSetPi0Shifting[9][11]      = {
                                            { 0,  0,  22,  2,  9,  0,  0,  0,  0,  21, 0 },  // 0-10
                                            { 0,  0,  18,  2,  9,  0,  0,  0,  0,  21, 0 },  // 10-20
                                            { 0,  0,  5,  2,  3,  0,  0,  0,  0,  21, 0 },  // 20-40
                                            { 0,  0,  5,  2,  3,  0,  0,  0,  0,  21, 0 },  // 40-60
                                            { 0,  0,  5,  2,  3,  0,  0,  0,  0,  21, 0 },  // 60-80
                                            { 0,  6,  8,  2,  5,  3,  0,  0,  0,  21, 0 },  // 0-5
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 5-10
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 },  // 0-20
                                            { 0,  5,  7,  1,  4,  0,  0,  0,  0,  21, 0 }   // 20-50
                                           };
    Int_t nComBinsPi0Shifting[9][11]    = {
                                            { 27, 33, 30, 27,  30, 0,  0,  0,  0,  0, 0 },  // 0-10
                                            { 27, 33, 30, 27,  30, 0,  0,  0,  0,  0, 0 },  // 10-20
                                            { 27, 32, 30, 26,  30, 0,  0,  0,  0,  0, 0 },  // 20-40
                                            { 23, 31, 30, 26,  30, 0,  0,  0,  0,  0, 0 },  // 40-60
                                            { 23, 31, 30, 26,  30, 0,  0,  0,  0,  0, 0 },  // 60-80
                                            { 29, 10, 17, 17,  29, 17,  0,  0,  0,  0, 0 },  // 0-5
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 5-10
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 },  // 0-20
                                            { 23, 10, 22, 20,  22, 0,  0,  0,  0,  0, 0 }   // 20-50
                                           };
    Double_t minPtPi0[9]                = { 0.4, 0.4, 0.4, 0.4, 0.4,   0.4, 0.4, 0.4, 0.4};

    TGraphAsymmErrors* statErrorRelCollectionPi0[9][11];
    TGraphAsymmErrors* sysErrorRelCollectionPi0[9][11];
    TGraph* graphWeightsPi0[9][11];
    TGraph* graphWeightsIndepMethodsPi0[9][11];

    for (Int_t cent = 0; cent < maxCent; cent++){
        for (Int_t meth = 0; meth< 11; meth++){
            graphWeightsPi0[cent][meth]                 = NULL;
            graphWeightsIndepMethodsPi0[cent][meth]     = NULL;
            statErrorRelCollectionPi0[cent][meth]       = NULL;
            sysErrorRelCollectionPi0[cent][meth]        = NULL;
        }
    }
    for (Int_t cent = 0; cent < maxCent; cent++){
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

    TGraphAsymmErrors* graphCombPi0InvYieldStat[9]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombPi0InvYieldStatWOXErr[9]       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombPi0InvYieldSys[9]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombPi0InvYieldTot[9]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombPi0InvYieldTotUnshi[9]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombPi0InvYieldStatUnshi[9]        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombPi0InvYieldSysUnshi[9]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombPi0InvYieldRelStat[9]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombPi0InvYieldRelSys[9]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombPi0InvYieldRelTot[9]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombPi0InvYieldStat_yShifted[9]    = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombPi0InvYieldSys_yShifted[9]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombPi0InvYieldTot_yShifted[9]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphRatioPi0CombCombFitTot[9]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphRatioPi0CombCombFitStat[9]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphRatioPi0CombCombFitSys[9]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphRatioPi0CombCombFitStatWOXErr[9]   = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoRatioPi0DPMJetToFit[9]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoRatioPi0HIJINGToFit[9]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoRatioPi0EPOSLHCToFit[9]                         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

    TGraphAsymmErrors* graphCombIndepMethodsPi0InvYieldStat[9]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombIndepMethodsPi0InvYieldSys[9]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombIndepMethodsPi0InvYieldTot[9]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombIndepMethodsPi0InvYieldRelStat[9]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombIndepMethodsPi0InvYieldRelSys[9]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombIndepMethodsPi0InvYieldRelTot[9]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

    TGraphAsymmErrors* graphIndPi0InvYieldStatUnshi[9][11];
    TGraphAsymmErrors* graphIndPi0InvYieldSysUnshi[9][11];
    TGraphAsymmErrors* graphIndPi0InvYieldStat[9][11];
    TGraphAsymmErrors* graphIndPi0InvYieldSys[9][11];
    TGraphAsymmErrors* graphIndPi0InvYieldStat_yShifted[9][11];
    TGraphAsymmErrors* graphIndPi0InvYieldSys_yShifted[9][11];
    TGraphAsymmErrors* graphRatioPi0IndCombFitStat[9][11];
    TGraphAsymmErrors* graphRatioPi0IndCombFitSys[9][11];
    for (Int_t cent = 0; cent < maxCent; cent++){
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
    TF1* fitTCMDecomposedLPi0[9]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TF1* fitTCMDecomposedHPi0[9]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TF1* fitTCMInvYieldPi0[9]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TF1* fitInvYieldPi0[9]                             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TF1* fitInvYieldPi0Graph[9]                        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TF1* fitPowInvYieldPi0[9]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

    // *************************************************************************************************************
    // ************************************** Define plotting environment ******************************************
    // *************************************************************************************************************
    textSizeLabelsPixel                 = 900*0.04;

    TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
    canvasWeights->SetLogx();
    TH2F * histo2DPi0Weights = new TH2F("histo2DPi0Weights","histo2DPi0Weights",11000,minPtPi0Plotting, maxPtPi0Plotting,1000,-0.7,1.3);
    SetStyleHistoTH2ForGraphs(histo2DPi0Weights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DPi0Weights->GetXaxis()->SetMoreLogLabels();
    histo2DPi0Weights->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelWeightsPi0         = new TLatex(0.95,0.15,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelWeightsPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelWeightsEnergy      = new TLatex(0.95,0.20,collisionSystemPbPb.Data());
    SetStyleTLatex( labelWeightsEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelWeightsPi0RAA     = new TLatex(0.95,0.15,"#it{R}_{PbPb}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelWeightsPi0RAA, 0.85*textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErr->SetLogx();
    TH2F * histo2DRelSysErr;
    histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,minPtPi0Plotting, maxPtPi0Plotting,1000,0,70.5);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErr->GetYaxis()->SetRangeUser(0,39.5);
    histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErr->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelRelSysErrEnergy    = new TLatex(0.15,0.89,collisionSystemPbPb.Data());
    SetStyleTLatex( labelRelSysErrEnergy, textSizeLabelsPixel, 4, 1, 43);
    TLatex *labelRelSysErrPi0       = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrPi0, textSizeLabelsPixel, 4, 1, 43);
    TLatex *labelRelSysErrPi0RAA    = new TLatex(0.15,0.85,"#it{R}_{PbPb}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrPi0RAA, textSizeLabelsPixel, 4, 1, 43);


    TCanvas* canvasRelStatErr           = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelStatErr->SetLogx();
    TH2F * histo2DRelStatErr;
    histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,minPtPi0Plotting, maxPtPi0Plotting,1000,0,70.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErr->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelRelStatErrEnergy   = new TLatex(0.95,0.89,collisionSystemPbPb.Data());
    SetStyleTLatex( labelRelStatErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelStatErrPi0      = new TLatex(0.95,0.85,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelStatErrPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelStatErrPi0RAA  = new TLatex(0.95,0.85,"#it{R}_{PbPb}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelStatErrPi0RAA, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    TCanvas* canvasRelTotErr            = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelTotErr->SetLogx();
    TH2F * histo2DRelTotErrPi0;
    histo2DRelTotErrPi0                 = new TH2F("histo2DRelTotErrPi0","histo2DRelTotErrPi0",11000,minPtPi0Plotting, maxPtPi0Plotting,1000,0,70.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrPi0, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,39.5);
    histo2DRelTotErrPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrPi0->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelRelTotErrEnergy    = new TLatex(0.95,0.89,collisionSystemPbPb.Data());
    SetStyleTLatex( labelRelTotErrEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelTotErrPi0       = new TLatex(0.95,0.85,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelTotErrPi0RAA       = new TLatex(0.95,0.85,"#it{R}_{PbPb}: #pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrPi0RAA, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);


    for (Int_t cent = 0; cent < maxCent; cent++){
        if (!enableCentComb[cent]) continue;
        cout << "start combining pi0 for " << centArray[cent].Data() << " " << endl;
        // Declaration & calculation of combined spectrum with all available techniques
        TString fileNamePi0OutputWeighting      = Form("%s/Pi0_WeightingMethod_%s%s.dat",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data());
        graphCombPi0InvYieldTot[cent]           = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionPi0[cent],    sysErrorCollectionPi0[cent],
                                                                                        xPtLimitsPi0[cent], maxNBinsPi0[cent],
                                                                                        offSetsPi0[cent], offSetsPi0Sys[cent],
                                                                                        graphCombPi0InvYieldStat[cent], graphCombPi0InvYieldSys[cent],
                                                                                        fileNamePi0OutputWeighting, "PbPb_5.02TeV", "Pi0", kTRUE,
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
        }


        // Declaration & calculation of combined spectrum with only fully independent techniques
        TString fileNamePi0OutputWeightingIndepMethod      = Form("%s/Pi0_WeightingMethodIndepMethod_%s%s.dat",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data());
        graphCombIndepMethodsPi0InvYieldTot[cent]          = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionOnlyIndepMethodsPi0[cent],    sysErrorCollectionOnlyIndepMethodsPi0[cent],
                                                                                                   xPtLimitsPi0[cent], maxNBinsPi0[cent],
                                                                                                   offSetsPi0[cent], offSetsPi0Sys[cent],
                                                                                                   graphCombIndepMethodsPi0InvYieldStat[cent], graphCombIndepMethodsPi0InvYieldSys[cent],
                                                                                                   fileNamePi0OutputWeightingIndepMethod, "PbPb_5.02TeV", "Pi0", kTRUE,
                                                                                                   NULL, fileNameCorrFactors, centArrayCorr[cent]
        );


        if (graphCombIndepMethodsPi0InvYieldTot[cent] == NULL) {
            cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
            return;
        }
        while (graphCombIndepMethodsPi0InvYieldStat[cent]->GetX()[0] < minPtPi0[cent]){
            graphCombIndepMethodsPi0InvYieldStat[cent]->RemovePoint(0);
        }
        while (graphCombIndepMethodsPi0InvYieldTot[cent]->GetX()[0] < minPtPi0[cent]){
            graphCombIndepMethodsPi0InvYieldTot[cent]->RemovePoint(0);
        }
        while (graphCombIndepMethodsPi0InvYieldSys[cent]->GetX()[0] < minPtPi0[cent]){
            graphCombIndepMethodsPi0InvYieldSys[cent]->RemovePoint(0);
        }

        // Reading weights from output file for plotting
        ifstream fileWeightsPi0ReadIndepMethods;
        fileWeightsPi0ReadIndepMethods.open(fileNamePi0OutputWeightingIndepMethod,ios_base::in);
        cout << "reading" << fileNamePi0OutputWeightingIndepMethod << endl;
        Double_t xValuesPi0ReadIndepMethods[100];
        Double_t weightsPi0ReadIndepMethods[11][100];
        Int_t availablePi0MeasIndepMethods[11]    = {   -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1,
                                            -1};
        Int_t nMeasSetPi0IndepMethods       = 3;
        Int_t nPtBinsPi0ReadIndepMethods    = 0;
        while(!fileWeightsPi0ReadIndepMethods.eof() && nPtBinsPi0ReadIndepMethods < 100){
            TString garbage             = "";
            if (nPtBinsPi0ReadIndepMethods == 0){
                fileWeightsPi0ReadIndepMethods >> garbage ;//>> availablePi0MeasIndepMethods[0] >> availablePi0MeasIndepMethods[1] >> availablePi0MeasIndepMethods[2] >> availablePi0MeasIndepMethods[3];
                for (Int_t i = 0; i < nMeasSetPi0IndepMethods; i++){
                    fileWeightsPi0ReadIndepMethods >> availablePi0MeasIndepMethods[i] ;
                }
                cout << "read following measurements: ";
                for (Int_t i = 0; i < nMeasSetPi0IndepMethods; i++){
                    cout << availablePi0MeasIndepMethods[i] << "\t" ;
                }
                cout << endl;
            } else {
                fileWeightsPi0ReadIndepMethods >> xValuesPi0ReadIndepMethods[nPtBinsPi0ReadIndepMethods-1];
                for (Int_t i = 0; i < nMeasSetPi0IndepMethods; i++){
                    fileWeightsPi0ReadIndepMethods >> weightsPi0ReadIndepMethods[availablePi0MeasIndepMethods[i]][nPtBinsPi0ReadIndepMethods-1] ;
                }
                cout << "read: "<<  nPtBinsPi0ReadIndepMethods << "\t"<< xValuesPi0ReadIndepMethods[nPtBinsPi0ReadIndepMethods-1] << "\t" ;
                for (Int_t i = 0; i < nMeasSetPi0IndepMethods; i++){
                    cout << weightsPi0ReadIndepMethods[availablePi0MeasIndepMethods[i]][nPtBinsPi0ReadIndepMethods-1] << "\t";
                }
                cout << endl;
            }
            nPtBinsPi0ReadIndepMethods++;
        }
        nPtBinsPi0ReadIndepMethods                  = nPtBinsPi0ReadIndepMethods-2 ;
        fileWeightsPi0ReadIndepMethods.close();

        for (Int_t i = 0; i < nMeasSetPi0IndepMethods; i++){
            graphWeightsIndepMethodsPi0[cent][availablePi0MeasIndepMethods[i]]                        = new TGraph(nPtBinsPi0ReadIndepMethods,xValuesPi0ReadIndepMethods,weightsPi0ReadIndepMethods[availablePi0MeasIndepMethods[i]]);
            Int_t bin = 0;
            for (Int_t n = 0; n< nPtBinsPi0ReadIndepMethods; n++){
                if (graphWeightsIndepMethodsPi0[cent][availablePi0MeasIndepMethods[i]]->GetY()[bin] == 0) graphWeightsIndepMethodsPi0[cent][availablePi0MeasIndepMethods[i]]->RemovePoint(bin);
                else bin++;
            }
        }


        // **********************************************************************************************************************
        // ********************************* Plotting weights method all available Methods **************************************
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

            labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            labelWeightsEnergy->Draw();
            labelWeightsPi0->Draw();

            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.2, 0.2,0.1, kGray, 3);

        canvasWeights->SaveAs(Form("%s/Pi0_Weights_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

        // **********************************************************************************************************************
        // ********************************* Plotting weights method all available Methods **************************************
        // **********************************************************************************************************************
        canvasWeights->cd();
        histo2DPi0Weights->Draw("copy");

        TLegend* legendWeightsOnlyIndepMethods   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetPi0IndepMethods), 32);
        for (Int_t i = 0; i < nMeasSetPi0IndepMethods; i++){
            DrawGammaSetMarkerTGraph(graphWeightsIndepMethodsPi0[cent][availablePi0MeasIndepMethods[i]], markerStyleDet[availablePi0MeasIndepMethods[i]], markerSizeDet[availablePi0MeasIndepMethods[i]]*0.5, colorDet[availablePi0MeasIndepMethods[i]] , colorDet[availablePi0MeasIndepMethods[i]]);
            graphWeightsIndepMethodsPi0[cent][availablePi0MeasIndepMethods[i]]->Draw("p,same,z");
            legendWeightsOnlyIndepMethods->AddEntry(graphWeightsIndepMethodsPi0[cent][availablePi0MeasIndepMethods[i]],nameMeasGlobalLabel[availablePi0MeasIndepMethods[i]],"p");
        }
        legendWeightsOnlyIndepMethods->Draw();

        labelWeightsPi0->Draw();

        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.2, 0.2,0.1, kGray, 3);

        canvasWeights->SaveAs(Form("%s/Pi0_WeightsOnlyIndepMethods_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

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

            }
            legendRelSysErr->Draw();


            labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            labelRelSysErrEnergy->Draw();
            labelRelSysErrPi0->Draw();

        canvasRelSysErr->SaveAs(Form("%s/Pi0_RelSysErr_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
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
            }
            legendRelStatErr->Draw();

            labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            labelRelStatErrEnergy->Draw();
            labelRelStatErrPi0->Draw();

        canvasRelStatErr->SaveAs(Form("%s/Pi0_RelStatErr_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

        //  *********************************************************************************************************************
        //  ************************ Visualize relative total errors of different combination methods Pi0 ***********************
        //  *********************************************************************************************************************

        graphCombPi0InvYieldRelStat[cent]               = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldStat[cent], "relativeStatErrorPi0_Method");
        graphCombPi0InvYieldRelSys[cent]                = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldSys[cent], "relativeSysErrorPi0_Method");
        graphCombPi0InvYieldRelTot[cent]                = CalculateRelErrUpAsymmGraph( graphCombPi0InvYieldTot[cent], "relativeTotalErrorPi0_Method");
        graphCombIndepMethodsPi0InvYieldRelStat[cent]   = CalculateRelErrUpAsymmGraph( graphCombIndepMethodsPi0InvYieldStat[cent], "relativeStatErrorPi0_IndepMethods");
        graphCombIndepMethodsPi0InvYieldRelSys[cent]    = CalculateRelErrUpAsymmGraph( graphCombIndepMethodsPi0InvYieldSys[cent], "relativeSysErrorPi0_IndepMethods");
        graphCombIndepMethodsPi0InvYieldRelTot[cent]    = CalculateRelErrUpAsymmGraph( graphCombIndepMethodsPi0InvYieldTot[cent], "relativeTotalErrorPi0_IndepMethods");

        canvasRelTotErr->cd();
        histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,39.5);
        histo2DRelTotErrPi0->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphCombPi0InvYieldRelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombIndepMethodsPi0InvYieldRelTot[cent], markerStyleComb+4, markerSizeComb, colorCombHighPt , colorCombHighPt);
            graphCombIndepMethodsPi0InvYieldRelTot[cent]->Draw("p,same,z");

            TLegend* legendRelTotErr1       = GetAndSetLegend2(0.20, 0.92-(0.035*2), 0.45, 0.92, 32);
            legendRelTotErr1->AddEntry(graphCombPi0InvYieldRelTot[cent],"All methods","p");
            legendRelTotErr1->AddEntry(graphCombIndepMethodsPi0InvYieldRelTot[cent],"Independent methods","p");
            legendRelTotErr1->Draw();


            labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            labelRelTotErrEnergy->Draw();
            labelRelTotErrPi0->Draw();

        canvasRelTotErr->SaveAs(Form("%s/Pi0_TotErr_Comp_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
        histo2DRelTotErrPi0->GetYaxis()->SetTitle("Err (%)");
        histo2DRelTotErrPi0->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphCombPi0InvYieldRelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelStat[cent], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
            graphCombPi0InvYieldRelStat[cent]->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelSys[cent], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
            graphCombPi0InvYieldRelSys[cent]->SetLineStyle(7);
            graphCombPi0InvYieldRelSys[cent]->Draw("l,x0,same,e1");

            TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.35, 0.92, 32);
            legendRelTotErr3->AddEntry(graphCombPi0InvYieldRelTot[cent],"tot","p");
            legendRelTotErr3->AddEntry(graphCombPi0InvYieldRelStat[cent],"stat","l");
            legendRelTotErr3->AddEntry(graphCombPi0InvYieldRelSys[cent],"sys","l");
            legendRelTotErr3->Draw();

            labelRelTotErrEnergy->Draw();
            labelRelTotErrPi0->Draw();

        canvasRelTotErr->SaveAs(Form("%s/Pi0_Reldecomp_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

        if (histChPiSpecStat[cent] && histChPiSpecSyst[cent]){
            histo2DRelTotErrPi0->GetYaxis()->SetTitle("Err (%)");
            histo2DRelTotErrPi0->Draw("copy");
                DrawGammaSetMarkerTGraphAsym(graphChPiSpecRelTotErr[cent], markerStyleComb, markerSizeComb, kBlue+1 , kBlue+1);
                graphChPiSpecRelTotErr[cent]->Draw("p,same,z");
                DrawGammaSetMarkerTGraphAsym(graphChPiSpecRelStatErr[cent], markerStyleComb, markerSizeComb, kBlue-4 , kBlue-4);
                graphChPiSpecRelStatErr[cent]->Draw("l,x0,same,e1");
                DrawGammaSetMarkerTGraphAsym(graphChPiSpecRelSystErr[cent], markerStyleComb, markerSizeComb, kBlue+3, kBlue+3);
                graphChPiSpecRelSystErr[cent]->SetLineStyle(8);
                graphChPiSpecRelSystErr[cent]->Draw("l,x0,same,e1");

                DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
                graphCombPi0InvYieldRelTot[cent]->Draw("p,same,z");
                DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelStat[cent], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
                graphCombPi0InvYieldRelStat[cent]->Draw("l,x0,same,e1");
                DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldRelSys[cent], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
                graphCombPi0InvYieldRelSys[cent]->SetLineStyle(7);
                graphCombPi0InvYieldRelSys[cent]->Draw("l,x0,same,e1");


                TLegend* legendRelTotErrWithCh      = GetAndSetLegend2(0.24, 0.92-(0.035*3), 0.55, 0.92, 32,2, "",43,0.3);
                legendRelTotErrWithCh->AddEntry(graphCombPi0InvYieldRelTot[cent],"#pi^{0} tot","p");
                legendRelTotErrWithCh->AddEntry(graphChPiSpecRelTotErr[cent],"#pi^{#pm} tot","p");
                legendRelTotErrWithCh->AddEntry(graphCombPi0InvYieldRelStat[cent],"#pi^{0} stat","l");
                legendRelTotErrWithCh->AddEntry(graphChPiSpecRelStatErr[cent],"#pi^{#pm} stat","l");
                legendRelTotErrWithCh->AddEntry(graphCombPi0InvYieldRelSys[cent],"#pi^{0} sys","l");
                legendRelTotErrWithCh->AddEntry(graphChPiSpecRelSystErr[cent],"#pi^{#pm} sys","l");
                legendRelTotErrWithCh->Draw();

                labelRelTotErrEnergy->Draw();
                labelRelTotErrPi0->Draw();

            canvasRelTotErr->SaveAs(Form("%s/Pi0_ReldecompWithChargedPions_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
        }

        histo2DRelTotErrPi0->GetYaxis()->SetTitle("Err (%)");
        histo2DRelTotErrPi0->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombIndepMethodsPi0InvYieldRelTot[cent], markerStyleComb, markerSizeComb, colorCombHighPt , colorCombHighPt);
            graphCombIndepMethodsPi0InvYieldRelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombIndepMethodsPi0InvYieldRelStat[cent], markerStyleComb, markerSizeComb, colorCombHighPt-6 , colorCombHighPt-6);
            graphCombIndepMethodsPi0InvYieldRelStat[cent]->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphCombIndepMethodsPi0InvYieldRelSys[cent], markerStyleComb, markerSizeComb, colorCombHighPt+2, colorCombHighPt+2);
            graphCombIndepMethodsPi0InvYieldRelSys[cent]->SetLineStyle(7);
            graphCombIndepMethodsPi0InvYieldRelSys[cent]->Draw("l,x0,same,e1");

            TLegend* legendRelTotErr7       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.35, 0.92, 32);
            legendRelTotErr7->AddEntry(graphCombIndepMethodsPi0InvYieldRelTot[cent],"tot","p");
            legendRelTotErr7->AddEntry(graphCombIndepMethodsPi0InvYieldRelStat[cent],"stat","l");
            legendRelTotErr7->AddEntry(graphCombIndepMethodsPi0InvYieldRelSys[cent],"sys","l");
            legendRelTotErr7->Draw();

            labelRelTotErrEnergy->Draw();
            labelRelTotErrPi0->Draw();

        canvasRelTotErr->SaveAs(Form("%s/Pi0_ReldecompOnlyIndepMethods_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

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
            graphCombPi0InvYieldStat[cent]         = ApplyXshiftIndividualSpectra(  graphCombPi0InvYieldTot[cent],
                                                                                    graphCombPi0InvYieldStat[cent],
                                                                                    fitShiftingPi0,
                                                                                    0, graphCombPi0InvYieldStat[cent]->GetN());
            graphCombPi0InvYieldSys[cent]          = ApplyXshiftIndividualSpectra(  graphCombPi0InvYieldTot[cent],
                                                                                    graphCombPi0InvYieldSys[cent],
                                                                                    fitShiftingPi0,
                                                                                    0, graphCombPi0InvYieldSys[cent]->GetN());
            for (Int_t meth = 0; meth< 11; meth++){
                if (graphIndPi0InvYieldStat[cent][meth]){
                    cout << "shiting stat err of " << nameMeasGlobalLabel[meth].Data();
                    graphIndPi0InvYieldStat[cent][meth]  = ApplyXshiftIndividualSpectra(    graphCombPi0InvYieldTot[cent],
                                                                                            graphIndPi0InvYieldStat[cent][meth],
                                                                                            fitShiftingPi0,
                                                                                            offSetPi0Shifting[cent][meth], nComBinsPi0Shifting[cent][meth]);

                }
                if (graphIndPi0InvYieldSys[cent][meth]){
                    cout << "shiting sys err of " << nameMeasGlobalLabel[meth].Data();
                    graphIndPi0InvYieldSys[cent][meth]   = ApplyXshiftIndividualSpectra(    graphCombPi0InvYieldTot[cent],
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

                TLatex *labelRatioToFitBinShift   = new TLatex(0.94, 0.91, Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
                SetStyleTLatex( labelRatioToFitBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitBinShift->Draw();
                TLatex *labelRatioToFitALICEBinShift    = new TLatex(0.94, 0.86, textALICE.Data());
                SetStyleTLatex( labelRatioToFitALICEBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitALICEBinShift->Draw();
                TLatex *labelRatioToFitPi0BinShift      = new TLatex(0.94, 0.807, "#pi^{0} #rightarrow #gamma#gamma");
                SetStyleTLatex( labelRatioToFitPi0BinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitPi0BinShift->Draw();

            canvasShift->Update();
            canvasShift->SaveAs(Form("%s/BinShiftCorrection_Pi0_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

            // *************************************************************************************************************
            // Plot control graphs
            // *************************************************************************************************************

            TCanvas* canvasDummy2       = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
            DrawGammaCanvasSettings( canvasDummy2,  0.15, 0.01, 0.015, 0.09);
            canvasDummy2->SetLogy();
            canvasDummy2->SetLogx();
            TH2F * histo2DDummy2;
            histo2DDummy2               = new TH2F("histo2DDummy2","histo2DDummy2",1000,minPtPi0Plotting,maxPtPi0Plotting,1000,1e-9,10e4);
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
            canvasDummy2->Print(Form("%s/ComparisonShiftedPi0_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
            delete canvasDummy2;
        }

        // *************************************************************************************************************
        // redo fitting after binshifts
        // *************************************************************************************************************
        // Tsallis function
        cout << WriteParameterToFile(fitInvYieldPi0[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitInvYieldPi0[cent])<< endl;
        //Two component model from Bylinkin
        fitTCMInvYieldPi0[cent]        = FitObject("tcm","fitTCMInvYieldPi0PbPb5TeVCent","Pi0",graphCombPi0InvYieldTot[cent], minPtPi0[cent], xPtLimitsPi0[cent][maxNBinsPi0[cent]],
                                                   paramTCMPi0New,"QNRMEX0+","", kFALSE);
        cout << WriteParameterToFile(fitTCMInvYieldPi0[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitTCMInvYieldPi0[cent])<< endl;

        fitPowInvYieldPi0[cent]        = FitObject("m","fitPowInvYieldPi0PbPb5TeVCent","Pi0",graphCombPi0InvYieldTot[cent],5,40. ,NULL,"QNRMEX0+","", kFALSE);
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
            histoRatioPi0DPMJetToFit[cent]                     = (TH1D*) histoDPMJetPi0[cent]->Clone(Form("histoRatioPi0DPMJetToFit%s",centArray[cent].Data()));
            histoRatioPi0DPMJetToFit[cent]                     = CalculateHistoRatioToFit (histoRatioPi0DPMJetToFit[cent], fitTCMInvYieldPi0[cent]);
            histoRatioPi0DPMJetToFit[cent]->GetXaxis()->SetRangeUser(minPtPi0[cent], xPtLimitsPi0[cent][maxNBinsPi0[cent]]);
        }
        if (histoHIJINGPi0[cent]){
            histoRatioPi0HIJINGToFit[cent]                     = (TH1D*) histoHIJINGPi0[cent]->Clone(Form("histoRatioPi0HIJINGToFit%s",centArray[cent].Data() ));
            histoRatioPi0HIJINGToFit[cent]                     = CalculateHistoRatioToFit (histoRatioPi0HIJINGToFit[cent], fitTCMInvYieldPi0[cent]);
            histoRatioPi0HIJINGToFit[cent]->GetXaxis()->SetRangeUser(minPtPi0[cent], xPtLimitsPi0[cent][maxNBinsPi0[cent]]);
        }
        if (histoEPOSLHCPi0[cent]){
            histoRatioPi0EPOSLHCToFit[cent]                    = (TH1D*) histoEPOSLHCPi0[cent]->Clone(Form("histoRatioPi0EPOSLHCToFit%s",centArray[cent].Data() ));
            histoRatioPi0EPOSLHCToFit[cent]                    = CalculateHistoRatioToFit (histoRatioPi0EPOSLHCToFit[cent], fitTCMInvYieldPi0[cent]);
            histoRatioPi0EPOSLHCToFit[cent]->GetXaxis()->SetRangeUser(minPtPi0[cent], xPtLimitsPi0[cent][maxNBinsPi0[cent]]);
        }

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

            Double_t textsizeLabelsPbPb      = 0;
            if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
                textsizeLabelsPbPb           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
            } else {
                textsizeLabelsPbPb           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
            }
            cout << textsizeLabelsPbPb << endl;

        TH2F * histo2DPi0RatioToCombFit;
        histo2DPi0RatioToCombFit               = new TH2F("histo2DPi0RatioToCombFit","histo2DPi0RatioToCombFit",1000,minPtPi0Plotting, maxPtPi0Plotting,1000,0.2,4.    );
        SetStyleHistoTH2ForGraphs(histo2DPi0RatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsPbPb, textsizeLabelsPbPb,
                                  0.85*textsizeLabelsPbPb,textsizeLabelsPbPb, 0.9, 0.65, 510, 505);
        histo2DPi0RatioToCombFit->GetXaxis()->SetMoreLogLabels();
        histo2DPi0RatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
        histo2DPi0RatioToCombFit->GetYaxis()->SetRangeUser(0.2,1.82);
        histo2DPi0RatioToCombFit->Draw("copy");

            ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStat[cent]);

            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys[cent], markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
            graphRatioPi0CombCombFitSys[cent]->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStat[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphRatioPi0CombCombFitStat[cent]->Draw("p,same,z");

            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1., 1.,1, kGray+2);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.1, 1.1,1, kGray, 7);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.9, 0.9,1, kGray, 7);

            TLatex *labelRatioToFitEnergy   = new TLatex(0.95, 0.92, Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
            labelRatioToFitEnergy->Draw();
            TLatex *labelRatioToFitALICE    = new TLatex(0.95, 0.86, textALICE.Data());
            SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
            labelRatioToFitALICE->Draw();
            TLatex *labelRatioToFitPi0      = new TLatex(0.12, 0.92, "#pi^{0} #rightarrow #gamma#gamma");
            SetStyleTLatex( labelRatioToFitPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
            labelRatioToFitPi0->Draw();

        canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfCombToCombFit_%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

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

            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1., 1.,1, kGray+2);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.1, 1.1,1, kGray, 7);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.9, 0.9,1, kGray, 7);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.2, 1.2,1, kGray, 9);
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.8, 0.8,1, kGray, 9);

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
            canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFit_%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(),
                                              suffix.Data()));
    }

    // *******************************************************************************************************
    // ************************** Combination of different eta measurements **********************************
    // *******************************************************************************************************
    // REMARKS:
    //     - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - extendable to other energies
    //     - offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented

    TH1D* statErrorCollectionEta[9][11];
    TGraphAsymmErrors* statErrorGraphCollectionEta[9][11];
    TGraphAsymmErrors* sysErrorCollectionEta[9][11];
    TH1D* statErrorCollectionOnlyIndepMethodsEta[9][11];
    TGraphAsymmErrors* statErrorGraphCollectionOnlyIndepMethodsEta[9][11];
    TGraphAsymmErrors* sysErrorCollectionOnlyIndepMethodsEta[9][11];
    for (Int_t cent = 0; cent < maxCent; cent++){
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
            if (meth < 3){
                statErrorCollectionOnlyIndepMethodsEta[cent][meth]          = NULL;
                if (histoEtaInvYieldStat[cent][meth]) statErrorCollectionOnlyIndepMethodsEta[cent][meth]      = (TH1D*)histoEtaInvYieldStat[cent][meth]->Clone(Form("statErrIndepMeth%i_%sEta",cent,nameMeasGlobalLabel[meth].Data()));
                statErrorGraphCollectionOnlyIndepMethodsEta[cent][meth]     = NULL;
                if (graphEtaInvYieldStat[cent][meth]) statErrorGraphCollectionOnlyIndepMethodsEta[cent][meth] = (TGraphAsymmErrors*)graphEtaInvYieldStat[cent][meth]->Clone(Form("statErrGraphIndepMeth%i_%sEta", cent,
                    nameMeasGlobalLabel[meth].Data()));
                if (statErrorGraphCollectionOnlyIndepMethodsEta[cent][meth]) RemoveZerosAtBeginningAndEndFromGraph(statErrorGraphCollectionOnlyIndepMethodsEta[cent][meth]);
                // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
                // for systematic error from respective method
                sysErrorCollectionOnlyIndepMethodsEta[cent][meth]           = NULL;
                if (graphEtaInvYieldSys[cent][meth]) sysErrorCollectionOnlyIndepMethodsEta[cent][meth]        = (TGraphAsymmErrors*)graphEtaInvYieldSys[cent][meth]->Clone(Form("sysErrIndepMeth%i_%sEta", cent,
                    nameMeasGlobalLabel[meth].Data()));
            } else {
                statErrorCollectionOnlyIndepMethodsEta[cent][meth]          = NULL;
                statErrorGraphCollectionOnlyIndepMethodsEta[cent][meth]     = NULL;
                sysErrorCollectionOnlyIndepMethodsEta[cent][meth]           = NULL;
            }
        }
    }

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Eta
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsEta[9][11]             = {
                                            { 0,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-10
                                            { 0,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 10-20
                                            { 0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 20-40
                                            { 0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 40-60
                                            { 0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 60-80
                                            { 0,  9,  0,  1,  0,  0,  0,  0,  0,  0,  0 },  // 0-5
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 5-10
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-20
                                            { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 }  // 20-50
                                          };
    Int_t offSetsEtaSys[9][11]          = {
                                            { 1,  4,  6,  2,  5,  0,  0,  0,  0,  0, 0 },  // 0-10
                                            { 1,  4,  5,  2,  3,  0,  0,  0,  0,  0, 0 },  // 10-20
                                            { 1,  3,  3,  2,  2,  0,  0,  0,  0,  0, 0 },  // 20-40
                                            { 1,  3,  3,  2,  2,  0,  0,  0,  0,  0, 0 },  // 40-60
                                            { 1,  3,  3,  2,  2,  0,  0,  0,  0,  0, 0 },  // 60-80
                                            { 3,  12,  7,  5,  5,  0,  0,  0,  0,  0, 0 },  // 0-5
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 5-10
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 0-20
                                            { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 }  // 20-50
                                            };
    Int_t offSetEtaShifting[9][11]      = {
                                            { 0,  3,  5,  1,  4,  0,  0,  0,  0,  21, 0 },  // 0-10
                                            { 0,  3,  4,  1,  2,  0,  0,  0,  0,  21, 0 },  // 10-20
                                            { 0,  2,  2,  1,  1,  0,  0,  0,  0,  21, 0 },  // 20-40
                                            { 0,  2,  2,  1,  1,  0,  0,  0,  0,  21, 0 },  // 40-60
                                            { 0,  2,  2,  1,  1,  0,  0,  0,  0,  21, 0 },  // 60-80
                                            { 0,  0,  4,  2,  2,  0,  0,  0,  0,  0,  0 },  // 0-5
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 },  // 5-10
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 },  // 0-20
                                            { 0,  0,  3,  2,  2,  0,  0,  0,  0,  21, 0 }  // 20-50
                                            };
    Int_t nComBinsEtaShifting[9][11]    = {
                                            { 8, 10, 10, 6,  7, 0,  0,  0,  0,  0, 0 },  // 0-10
                                            { 8, 10, 10, 6,  7, 0,  0,  0,  0,  0, 0 },  // 10-20
                                            { 7, 10, 10, 6,  9, 0,  0,  0,  0,  0, 0 },  // 20-40
                                            { 7, 9, 10, 6,  9, 0,  0,  0,  0,  0, 0 },  // 40-60
                                            { 7, 9, 10, 6,  9, 0,  0,  0,  0,  0, 0 },  // 60-80
                                            { 13, 0, 17, 11, 15,  0,  0,  0,  0,  0, 0 },  // 0-5
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // 5-10
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // 0-20
                                            { 9, 0, 10, 0,  7, 0,  0,  0,  0,  0, 0 },  // 20-50
                                           };
    Double_t minPtEta[9]                = { 0.8, 0.8, 0.8, 0.8, 0.8,   0.7, 0.7, 0.7, 0.7};

    TGraphAsymmErrors* statErrorRelCollectionEta[9][11];
    TGraphAsymmErrors* sysErrorRelCollectionEta[9][11];
    TGraph* graphWeightsEta[9][11];
    TGraph* graphWeightsIndepMethodsEta[9][11];
    for (Int_t cent = 0; cent < maxCent; cent++){
        for (Int_t meth = 0; meth< 11; meth++){
            graphWeightsEta[cent][meth]                 = NULL;
            graphWeightsIndepMethodsEta[cent][meth]     = NULL;
            statErrorRelCollectionEta[cent][meth]       = NULL;
            sysErrorRelCollectionEta[cent][meth]        = NULL;
        }
    }
    for (Int_t cent = 0; cent < maxCent; cent++){
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

    TGraphAsymmErrors* graphCombEtaInvYieldStat[9]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaInvYieldStatWOXErr[9]       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaInvYieldSys[9]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaInvYieldTot[9]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaInvYieldTotUnshi[9]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaInvYieldStatUnshi[9]        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaInvYieldSysUnshi[9]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaInvYieldRelStat[9]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaInvYieldRelSys[9]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaInvYieldRelTot[9]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaInvYieldStat_yShifted[9]    = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaInvYieldSys_yShifted[9]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaInvYieldTot_yShifted[9]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphRatioEtaCombCombFitTot[9]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphRatioEtaCombCombFitStat[9]         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphRatioEtaCombCombFitStatWOXErr[9]   = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphRatioEtaCombCombFitSys[9]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoRatioEtaDPMJetToFit[9]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoRatioEtaHIJINGToFit[9]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D* histoRatioEtaEPOSLHCToFit[9]                         = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

    TGraphAsymmErrors* graphCombIndepMethodsEtaInvYieldStat[9]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombIndepMethodsEtaInvYieldSys[9]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombIndepMethodsEtaInvYieldTot[9]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombIndepMethodsEtaInvYieldRelStat[9]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombIndepMethodsEtaInvYieldRelSys[9]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombIndepMethodsEtaInvYieldRelTot[9]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

    TGraphAsymmErrors* graphIndEtaInvYieldStatUnshi[9][11];
    TGraphAsymmErrors* graphIndEtaInvYieldSysUnshi[9][11];
    TGraphAsymmErrors* graphIndEtaInvYieldStat[9][11];
    TGraphAsymmErrors* graphIndEtaInvYieldSys[9][11];
    TGraphAsymmErrors* graphIndEtaInvYieldStat_yShifted[9][11];
    TGraphAsymmErrors* graphIndEtaInvYieldSys_yShifted[9][11];
    TGraphAsymmErrors* graphRatioEtaIndCombFitStat[9][11];
    TGraphAsymmErrors* graphRatioEtaIndCombFitSys[9][11];
    for (Int_t cent = 0; cent < maxCent; cent++){
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
    TF1* fitTCMDecomposedLEta[9]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TF1* fitTCMDecomposedHEta[9]                       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TF1* fitTCMInvYieldEta[9]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TF1* fitInvYieldEta[9]                             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TF1* fitInvYieldEtaGraph[9]                        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TF1* fitPowInvYieldEta[9]                          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

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
    TLatex *labelWeightsEtaRAA         = new TLatex(0.95,0.15,"#it{R}_{pA}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelWeightsEtaRAA, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    TLatex *labelRelSysErrEta       = new TLatex(0.15,0.85,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrEta, textSizeLabelsPixel, 4, 1, 43);
    TLatex *labelRelSysErrEtaRAA       = new TLatex(0.15,0.85,"#it{R}_{pA}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelSysErrEtaRAA, textSizeLabelsPixel, 4, 1, 43);

    TLatex *labelRelStatErrEta      = new TLatex(0.95,0.85,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelStatErrEta, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelStatErrEtaRAA      = new TLatex(0.95,0.85,"#it{R}_{pA}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelStatErrEtaRAA, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    TH2F * histo2DRelTotErrEta;
    histo2DRelTotErrEta                 = new TH2F("histo2DRelTotErrEta","histo2DRelTotErrEta",11000,minPtEtaPlotting, maxPtEtaPlotting,1000,0,70.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrEta, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,39.5);
    histo2DRelTotErrEta->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrEta->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelRelTotErrEta       = new TLatex(0.95,0.85,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrEta, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
    TLatex *labelRelTotErrEtaRAA   = new TLatex(0.95,0.85,"#it{R}_{pA}: #eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrEtaRAA, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    for (Int_t cent = 0; cent < maxCent; cent++){
        if (!enableCentComb[cent]) continue;
        // Declaration & calculation of combined spectrum
        cout << "start combining eta for " << centArray[cent].Data()  << endl;
        TString fileNameEtaOutputWeighting      = Form("%s/Eta_WeightingMethod_%s%s.dat",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data());
        graphCombEtaInvYieldTot[cent]           = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionEta[cent],    sysErrorCollectionEta[cent],
                                                                                        xPtLimitsEta[cent], maxNBinsEta[cent],
                                                                                        offSetsEta[cent], offSetsEtaSys[cent],
                                                                                        graphCombEtaInvYieldStat[cent], graphCombEtaInvYieldSys[cent],
                                                                                        fileNameEtaOutputWeighting, "PbPb_5.02TeV", "Eta", kTRUE,
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

        // Declaration & calculation of combined spectrum with only fully independent techniques
        TString fileNameEtaOutputWeightingIndepMethod      = Form("%s/Eta_WeightingMethodIndepMethod_%s%s.dat",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data());
        graphCombIndepMethodsEtaInvYieldTot[cent]          = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionOnlyIndepMethodsEta[cent],    sysErrorCollectionOnlyIndepMethodsEta[cent],
                                                                                                   xPtLimitsEta[cent], maxNBinsEta[cent],
                                                                                                   offSetsEta[cent], offSetsEtaSys[cent],
                                                                                                   graphCombIndepMethodsEtaInvYieldStat[cent], graphCombIndepMethodsEtaInvYieldSys[cent],
                                                                                                   fileNameEtaOutputWeightingIndepMethod, "PbPb_5.02TeV", "Eta", kTRUE,
                                                                                                   NULL, fileNameCorrFactors, centArrayCorr[cent]
        );


        if (graphCombIndepMethodsEtaInvYieldTot[cent] == NULL) {
            cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
            return;
        }
        while (graphCombIndepMethodsEtaInvYieldStat[cent]->GetX()[0] < minPtEta[cent]){
            graphCombIndepMethodsEtaInvYieldStat[cent]->RemovePoint(0);
        }
        while (graphCombIndepMethodsEtaInvYieldTot[cent]->GetX()[0] < minPtEta[cent]){
            graphCombIndepMethodsEtaInvYieldTot[cent]->RemovePoint(0);
        }
        while (graphCombIndepMethodsEtaInvYieldSys[cent]->GetX()[0] < minPtEta[cent]){
            graphCombIndepMethodsEtaInvYieldSys[cent]->RemovePoint(0);
        }

        // Reading weights from output file for plotting
        ifstream fileWeightsEtaReadIndepMethods;
        fileWeightsEtaReadIndepMethods.open(fileNameEtaOutputWeightingIndepMethod,ios_base::in);
        cout << "reading" << fileNameEtaOutputWeightingIndepMethod << endl;
        Double_t xValuesEtaReadIndepMethods[100];
        Double_t weightsEtaReadIndepMethods[11][100];
        Int_t availableEtaMeasIndepMethods[11]    = {   -1, -1, -1, -1, -1,
                                                        -1, -1, -1, -1, -1,
                                                        -1};
        Int_t nMeasSetEtaIndepMethods       = 3;
        Int_t nPtBinsEtaReadIndepMethods    = 0;
        while(!fileWeightsEtaReadIndepMethods.eof() && nPtBinsEtaReadIndepMethods < 100){
            TString garbage             = "";
            if (nPtBinsEtaReadIndepMethods == 0){
                fileWeightsEtaReadIndepMethods >> garbage ;//>> availableEtaMeasIndepMethods[0] >> availableEtaMeasIndepMethods[1] >> availableEtaMeasIndepMethods[2] >> availableEtaMeasIndepMethods[3];
                for (Int_t i = 0; i < nMeasSetEtaIndepMethods; i++){
                    fileWeightsEtaReadIndepMethods >> availableEtaMeasIndepMethods[i] ;
                }
                cout << "read following measurements: ";
                for (Int_t i = 0; i < nMeasSetEtaIndepMethods; i++){
                    cout << availableEtaMeasIndepMethods[i] << "\t" ;
                }
                cout << endl;
            } else {
                fileWeightsEtaReadIndepMethods >> xValuesEtaReadIndepMethods[nPtBinsEtaReadIndepMethods-1];
                for (Int_t i = 0; i < nMeasSetEtaIndepMethods; i++){
                    fileWeightsEtaReadIndepMethods >> weightsEtaReadIndepMethods[availableEtaMeasIndepMethods[i]][nPtBinsEtaReadIndepMethods-1] ;
                }
                cout << "read: "<<  nPtBinsEtaReadIndepMethods << "\t"<< xValuesEtaReadIndepMethods[nPtBinsEtaReadIndepMethods-1] << "\t" ;
                for (Int_t i = 0; i < nMeasSetEtaIndepMethods; i++){
                    cout << weightsEtaReadIndepMethods[availableEtaMeasIndepMethods[i]][nPtBinsEtaReadIndepMethods-1] << "\t";
                }
                cout << endl;
            }
            nPtBinsEtaReadIndepMethods++;
        }
        nPtBinsEtaReadIndepMethods                  = nPtBinsEtaReadIndepMethods-2 ;
        fileWeightsEtaReadIndepMethods.close();

        for (Int_t i = 0; i < nMeasSetEtaIndepMethods; i++){
            graphWeightsIndepMethodsEta[cent][availableEtaMeasIndepMethods[i]]                        = new TGraph(nPtBinsEtaReadIndepMethods,xValuesEtaReadIndepMethods,weightsEtaReadIndepMethods[availableEtaMeasIndepMethods[i]]);
            Int_t bin = 0;
            for (Int_t n = 0; n< nPtBinsEtaReadIndepMethods; n++){
                if (graphWeightsIndepMethodsEta[cent][availableEtaMeasIndepMethods[i]]->GetY()[bin] == 0) graphWeightsIndepMethodsEta[cent][availableEtaMeasIndepMethods[i]]->RemovePoint(bin);
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

            labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            labelWeightsEnergy->Draw();
            labelWeightsEta->Draw();

            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.2, 0.2,0.1, kGray, 3);

        canvasWeights->SaveAs(Form("%s/Eta_Weights_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

        canvasWeights->cd();
        histo2DEtaWeights->Draw("copy");

            TLegend* legendWeightsOnlyIndepMethods   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEtaIndepMethods), 32);
            for (Int_t i = 0; i < nMeasSetEtaIndepMethods; i++){
                DrawGammaSetMarkerTGraph(graphWeightsIndepMethodsEta[cent][availableEtaMeasIndepMethods[i]], markerStyleDet[availableEtaMeasIndepMethods[i]], markerSizeDet[availableEtaMeasIndepMethods[i]]*0.5, colorDet[availableEtaMeasIndepMethods[i]] , colorDet[availableEtaMeasIndepMethods[i]]);
                graphWeightsIndepMethodsEta[cent][availableEtaMeasIndepMethods[i]]->Draw("p,same,z");
                legendWeightsOnlyIndepMethods->AddEntry(graphWeightsIndepMethodsEta[cent][availableEtaMeasIndepMethods[i]],nameMeasGlobalLabel[availableEtaMeasIndepMethods[i]],"p");
            }
            legendWeightsOnlyIndepMethods->Draw();

            labelWeightsEta->Draw();

            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.2, 0.2,0.1, kGray, 3);

        canvasWeights->SaveAs(Form("%s/Eta_WeightsOnlyIndepMethods_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

        //  *********************************************************************************************************************
        //  ************************************ Visualize relative errors ******************************************************
        //  *********************************************************************************************************************


        canvasRelSysErr->cd();
        histo2DRelSysErr->GetYaxis()->SetRangeUser(0,65.5);
        histo2DRelSysErr->GetXaxis()->SetRangeUser(minPtEtaPlotting,maxPtEtaPlotting);
        histo2DRelSysErr->Draw("copy");

            TLegend* legendRelSysErrEta        = GetAndSetLegend2(0.62, 0.92-(0.04*nMeasSetEta), 0.95, 0.92, textSizeLabelsPixel);
            for (Int_t i = 0; i < nMeasSetEta; i++){
                DrawGammaSetMarkerTGraph(sysErrorRelCollectionEta[cent][availableEtaMeas[i]], markerStyleDet[availableEtaMeas[i]], markerSizeDet[availableEtaMeas[i]]*0.5, colorDet[availableEtaMeas[i]],
                                         colorDet[availableEtaMeas[i]]);
                sysErrorRelCollectionEta[cent][availableEtaMeas[i]]->Draw("p,same,z");
                legendRelSysErrEta->AddEntry(sysErrorRelCollectionEta[cent][availableEtaMeas[i]],nameMeasGlobalLabel[availableEtaMeas[i]],"p");
            }
            legendRelSysErrEta->Draw();

            labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            labelRelSysErrEnergy->Draw();
            labelRelSysErrEta->Draw();

        canvasRelSysErr->SaveAs(Form("%s/Eta_RelSysErr_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
        delete legendRelSysErrEta;
        //  *********************************************************************************************************************
        //  ************************************ Visualize relative errors ******************************************************
        //  *********************************************************************************************************************

        canvasRelStatErr->cd();
        histo2DRelStatErr->GetYaxis()->SetRangeUser(0,65.5);
        histo2DRelStatErr->GetXaxis()->SetRangeUser(minPtEtaPlotting,maxPtEtaPlotting);
        histo2DRelStatErr->Draw("copy");
            TLegend* legendRelStatErr       = GetAndSetLegend2(0.14, 0.92-(0.04*nMeasSetEta), 0.45, 0.92, textSizeLabelsPixel);
            for (Int_t i = 0; i < nMeasSetEta; i++){
                DrawGammaSetMarkerTGraph(statErrorRelCollectionEta[cent][availableEtaMeas[i]], markerStyleDet[availableEtaMeas[i]], markerSizeDet[availableEtaMeas[i]]*0.5, colorDet[availableEtaMeas[i]],
                                         colorDet[availableEtaMeas[i]]);
                statErrorRelCollectionEta[cent][availableEtaMeas[i]]->Draw("p,same,z");
                legendRelStatErr->AddEntry(statErrorRelCollectionEta[cent][availableEtaMeas[i]],nameMeasGlobalLabel[availableEtaMeas[i]],"p");
            }
            legendRelStatErr->Draw();

            labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            labelRelStatErrEnergy->Draw();
            labelRelStatErrEta->Draw();

        canvasRelStatErr->SaveAs(Form("%s/Eta_RelStatErr_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

        //  *********************************************************************************************************************
        //  ************************ Visualize relative total errors of different combination methods Eta ***********************
        //  *********************************************************************************************************************


        graphCombEtaInvYieldRelStat[cent]               = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldStat[cent], "relativeStatErrorEta_Method");
        graphCombEtaInvYieldRelSys[cent]                = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldSys[cent], "relativeSysErrorEta_Method");
        graphCombEtaInvYieldRelTot[cent]                = CalculateRelErrUpAsymmGraph( graphCombEtaInvYieldTot[cent], "relativeTotalErrorEta_Method");
        graphCombIndepMethodsEtaInvYieldRelStat[cent]   = CalculateRelErrUpAsymmGraph( graphCombIndepMethodsEtaInvYieldStat[cent], "relativeStatErrorEta_IndepMethods");
        graphCombIndepMethodsEtaInvYieldRelSys[cent]    = CalculateRelErrUpAsymmGraph( graphCombIndepMethodsEtaInvYieldSys[cent], "relativeSysErrorEta_IndepMethods");
        graphCombIndepMethodsEtaInvYieldRelTot[cent]    = CalculateRelErrUpAsymmGraph( graphCombIndepMethodsEtaInvYieldTot[cent], "relativeTotalErrorEta_IndepMethods");

        canvasRelTotErr->cd();
        histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,65.5);
        histo2DRelTotErrEta->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldRelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphCombEtaInvYieldRelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombIndepMethodsEtaInvYieldRelTot[cent], markerStyleComb+4, markerSizeComb, colorCombHighPt , colorCombHighPt);
            graphCombIndepMethodsEtaInvYieldRelTot[cent]->Draw("p,same,z");

            TLegend* legendRelTotErr1       = GetAndSetLegend2(0.20, 0.92-(0.035*2), 0.45, 0.92, 32);
            legendRelTotErr1->AddEntry(graphCombEtaInvYieldRelTot[cent],"All methods","p");
            legendRelTotErr1->AddEntry(graphCombIndepMethodsEtaInvYieldRelTot[cent],"Independent methods","p");
            legendRelTotErr1->Draw();


            labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            labelRelTotErrEnergy->Draw();
            labelRelTotErrEta->Draw();

        canvasRelTotErr->SaveAs(Form("%s/Eta_TotErr_Comp_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

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

        canvasRelTotErr->SaveAs(Form("%s/Eta_Reldecomp_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
        histo2DRelTotErrEta->GetYaxis()->SetTitle("Err (%)");
        histo2DRelTotErrEta->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombIndepMethodsEtaInvYieldRelTot[cent], markerStyleComb, markerSizeComb, colorCombHighPt , colorCombHighPt);
            graphCombIndepMethodsEtaInvYieldRelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombIndepMethodsEtaInvYieldRelStat[cent], markerStyleComb, markerSizeComb, colorCombHighPt-6 , colorCombHighPt-6);
            graphCombIndepMethodsEtaInvYieldRelStat[cent]->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphCombIndepMethodsEtaInvYieldRelSys[cent], markerStyleComb, markerSizeComb, colorCombHighPt+2, colorCombHighPt+2);
            graphCombIndepMethodsEtaInvYieldRelSys[cent]->SetLineStyle(7);
            graphCombIndepMethodsEtaInvYieldRelSys[cent]->Draw("l,x0,same,e1");

            TLegend* legendRelTotErr7       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.35, 0.92, 32);
            legendRelTotErr7->AddEntry(graphCombIndepMethodsEtaInvYieldRelTot[cent],"tot","p");
            legendRelTotErr7->AddEntry(graphCombIndepMethodsEtaInvYieldRelStat[cent],"stat","l");
            legendRelTotErr7->AddEntry(graphCombIndepMethodsEtaInvYieldRelSys[cent],"sys","l");
            legendRelTotErr7->Draw();

            labelRelTotErrEnergy->Draw();
            labelRelTotErrEta->Draw();

        canvasRelTotErr->SaveAs(Form("%s/Eta_ReldecompOnlyIndepMethods_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

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

                TLatex *labelRatioToFitBinShift   = new TLatex(0.94, 0.91, Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
                SetStyleTLatex( labelRatioToFitBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitBinShift->Draw();
                TLatex *labelRatioToFitALICEBinShift    = new TLatex(0.94, 0.86, textALICE.Data());
                SetStyleTLatex( labelRatioToFitALICEBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitALICEBinShift->Draw();
                TLatex *labelRatioToFitEtaBinShift      = new TLatex(0.94, 0.807, "#eta #rightarrow #gamma#gamma");
                SetStyleTLatex( labelRatioToFitEtaBinShift, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
                labelRatioToFitEtaBinShift->Draw();

            canvasShift->Update();
            canvasShift->SaveAs(Form("%s/BinShiftCorrection_Eta_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

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
            canvasDummy2->Print(Form("%s/ComparisonShiftedEta_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
            delete canvasDummy2;
        }

        // *************************************************************************************************************
        // redo fitting after binshifts
        // *************************************************************************************************************
        // Tsallis function
        cout << WriteParameterToFile(fitInvYieldEta[cent])<< endl;
        fileFitsOutput <<  WriteParameterToFile(fitInvYieldEta[cent])<< endl;
        //Two component model from Bylinkin
        fitTCMInvYieldEta[cent]        = FitObject("tcm","fitTCMInvYieldEtaPbPb5TeVCent","Eta",graphCombEtaInvYieldTot[cent],
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

        fitPowInvYieldEta[cent]        = FitObject("m","fitPowInvYieldEtaPbPb5TeVCent","Eta",graphCombEtaInvYieldTot[cent], 5, 40. ,NULL,"QNRMEX0+","", kFALSE);
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
            histoRatioEtaDPMJetToFit[cent]                     = (TH1D*) histoDPMJetEta[cent]->Clone(Form("histoRatioEtaDPMJetToFit%s",centArray[cent].Data() ));
            histoRatioEtaDPMJetToFit[cent]                     = CalculateHistoRatioToFit (histoRatioEtaDPMJetToFit[cent], fitTCMInvYieldEta[cent]);
            histoRatioEtaDPMJetToFit[cent]->GetXaxis()->SetRangeUser(minPtEta[cent], xPtLimitsEta[cent][maxNBinsEta[cent]]);
        }
        if (histoHIJINGEta[cent]){
            histoRatioEtaHIJINGToFit[cent]                     = (TH1D*) histoHIJINGEta[cent]->Clone(Form("histoRatioEtaHIJINGToFit%s",centArray[cent].Data()));
            histoRatioEtaHIJINGToFit[cent]                     = CalculateHistoRatioToFit (histoRatioEtaHIJINGToFit[cent], fitTCMInvYieldEta[cent]);
            histoRatioEtaHIJINGToFit[cent]->GetXaxis()->SetRangeUser(minPtEta[cent], xPtLimitsEta[cent][maxNBinsEta[cent]]);
        }
        if (histoEPOSLHCEta[cent]){
            histoRatioEtaEPOSLHCToFit[cent]                    = (TH1D*) histoEPOSLHCEta[cent]->Clone(Form("histoRatioEtaEPOSLHCToFit%s",centArray[cent].Data() ));
            histoRatioEtaEPOSLHCToFit[cent]                    = CalculateHistoRatioToFit (histoRatioEtaEPOSLHCToFit[cent], fitTCMInvYieldEta[cent]);
            histoRatioEtaEPOSLHCToFit[cent]->GetXaxis()->SetRangeUser(minPtEta[cent], xPtLimitsEta[cent][maxNBinsEta[cent]]);
        }

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

        Double_t textsizeLabelsPbPb      = 0;
        if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
            textsizeLabelsPbPb           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
        } else {
            textsizeLabelsPbPb           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
        }

        TH2F * histo2DEtaRatioToCombFit;
        histo2DEtaRatioToCombFit               = new TH2F("histo2DEtaRatioToCombFit","histo2DEtaRatioToCombFit",1000,minPtEtaPlotting, maxPtEtaPlotting,1000,0.2,4.    );
        SetStyleHistoTH2ForGraphs(histo2DEtaRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsPbPb, textsizeLabelsPbPb,
                                  0.85*textsizeLabelsPbPb,textsizeLabelsPbPb, 0.9, 0.65, 510, 505);
        histo2DEtaRatioToCombFit->GetXaxis()->SetMoreLogLabels();
        histo2DEtaRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
        histo2DEtaRatioToCombFit->GetYaxis()->SetRangeUser(0.2,1.82);
        histo2DEtaRatioToCombFit->Draw("copy");

            ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombCombFitStat[cent]);

            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys[cent], markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
            graphRatioEtaCombCombFitSys[cent]->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStat[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphRatioEtaCombCombFitStat[cent]->Draw("p,same,z");

            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1., 1.,1, kGray+2);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1.1, 1.1,1, kGray, 7);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.9, 0.9,1, kGray, 7);

            TLatex *labelRatioToFitEnergy   = new TLatex(0.95, 0.92, Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
            labelRatioToFitEnergy->Draw();
            TLatex *labelRatioToFitALICE    = new TLatex(0.95, 0.86, textALICE.Data());
            SetStyleTLatex( labelRatioToFitALICE, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
            labelRatioToFitALICE->Draw();
            TLatex *labelRatioToFitEta      = new TLatex(0.12, 0.92, "#eta #rightarrow #gamma#gamma");
            SetStyleTLatex( labelRatioToFitEta, textSizeLabelsPixel, 4, 1, 43, kTRUE, 11);
            labelRatioToFitEta->Draw();

        canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfCombToCombFit_%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data() ,suffix.Data()));

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

            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1., 1.,1, kGray+2);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1.1, 1.1,1, kGray, 7);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.9, 0.9,1, kGray, 7);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 1.2, 1.2,1, kGray, 9);
            DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting , 0.8, 0.8,1, kGray, 9);

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
        canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfIndividualMeasToCombFit_%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(),
                                          suffix.Data()));
    }


        // *******************************************************************************************************
    // ************************** Combination of different eta measurements **********************************
    // *******************************************************************************************************
    // REMARKS:
    //     - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //     - extendable to other energies
    //     - offsets have to be determined manually, see cout's in shell from combination function, more can be uncommented

    TH1D* statErrorCollectionEtaToPi0[9][11];
    TGraphAsymmErrors* statErrorGraphCollectionEtaToPi0[9][11];
    TGraphAsymmErrors* sysErrorCollectionEtaToPi0[9][11];
    TH1D* statErrorCollectionOnlyIndepMethodsEtaToPi0[9][11];
    TGraphAsymmErrors* statErrorGraphCollectionOnlyIndepMethodsEtaToPi0[9][11];
    TGraphAsymmErrors* sysErrorCollectionOnlyIndepMethodsEtaToPi0[9][11];
    for (Int_t cent = 0; cent < maxCent; cent++){
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
            if (meth < 3){
                statErrorCollectionOnlyIndepMethodsEtaToPi0[cent][meth]          = NULL;
                if (histoEtaToPi0Stat[cent][meth]) statErrorCollectionOnlyIndepMethodsEtaToPi0[cent][meth]      = (TH1D*)histoEtaToPi0Stat[cent][meth]->Clone(Form("statErrIndepMeth%i_%sEtaToPi0",cent,nameMeasGlobalLabel[meth].Data()));
                statErrorGraphCollectionOnlyIndepMethodsEtaToPi0[cent][meth]     = NULL;
                if (graphEtaToPi0Stat[cent][meth]) statErrorGraphCollectionOnlyIndepMethodsEtaToPi0[cent][meth] = (TGraphAsymmErrors*)graphEtaToPi0Stat[cent][meth]->Clone(Form("statErrGraphIndepMeth%i_%sEtaToPi0", cent,
                    nameMeasGlobalLabel[meth].Data()));
                if (statErrorGraphCollectionOnlyIndepMethodsEtaToPi0[cent][meth]) RemoveZerosAtBeginningAndEndFromGraph(statErrorGraphCollectionOnlyIndepMethodsEtaToPi0[cent][meth]);
                // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
                // for systematic error from respective method
                sysErrorCollectionOnlyIndepMethodsEtaToPi0[cent][meth]           = NULL;
                if (graphEtaToPi0Sys[cent][meth]) sysErrorCollectionOnlyIndepMethodsEtaToPi0[cent][meth]        = (TGraphAsymmErrors*)graphEtaToPi0Sys[cent][meth]->Clone(Form("sysErrIndepMeth%i_%sEtaToPi0", cent,
                    nameMeasGlobalLabel[meth].Data()));
            } else {
                statErrorCollectionOnlyIndepMethodsEtaToPi0[cent][meth]          = NULL;
                statErrorGraphCollectionOnlyIndepMethodsEtaToPi0[cent][meth]     = NULL;
                sysErrorCollectionOnlyIndepMethodsEtaToPi0[cent][meth]           = NULL;
            }
        }
    }

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for EtaToPi0
    // {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMC", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "spare", "EMCAL merged","PCMOtherDataset"};
    Int_t offSetsEtaToPi0[9][11]           = {
                                                { 0,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-10
                                                { 0,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 10-20
                                                { 0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 20-40
                                                { 0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 40-60
                                                { 0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 60-80
                                                { 0,  9,  0,  1,  0,  0,  0,  0,  0,  0,  0 },  // 0-5
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 5-10
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 },  // 0-20
                                                { 0,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0 }   // 20-50
                                              };
    Int_t offSetsEtaToPi0Sys[9][11]        = {
                                                { 1,  4,  6,  2,  5,  0,  0,  0,  0,  0, 0  },  // 0-10
                                                { 1,  4,  5,  2,  3,  0,  0,  0,  0,  0, 0 },  // 10-20
                                                { 1,  3,  3,  2,  2,  0,  0,  0,  0,  0, 0 },  // 20-40
                                                { 1,  3,  3,  2,  2,  0,  0,  0,  0,  0, 0 },  // 40-60
                                                { 1,  3,  3,  2,  2,  0,  0,  0,  0,  0, 0 },  // 60-80
                                                { 3,  12,  7,  5,  5,  0,  0,  0,  0,  0, 0 },  // 0-5
                                                { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 5-10
                                                { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 },  // 0-20
                                                { 3,  10,  6,  5,  5,  0,  0,  0,  0,  0, 0 }  // 20-50
                                              };
    Double_t minPtEtaToPi0[9]              = { 0.8, 0.8, 0.8, 0.8, 0.8,   0.7, 0.7, 0.7, 0.7};

    TGraphAsymmErrors* statErrorRelCollectionEtaToPi0[9][11];
    TGraphAsymmErrors* sysErrorRelCollectionEtaToPi0[9][11];
    TGraph* graphWeightsEtaToPi0[9][11];
    TGraph* graphWeightsIndepMethodsEtaToPi0[9][11];
    for (Int_t cent = 0; cent < maxCent; cent++){
        for (Int_t meth = 0; meth< 11; meth++){
            graphWeightsEtaToPi0[cent][meth]                 = NULL;
            graphWeightsIndepMethodsEtaToPi0[cent][meth]     = NULL;
            statErrorRelCollectionEtaToPi0[cent][meth]       = NULL;
            sysErrorRelCollectionEtaToPi0[cent][meth]        = NULL;
        }
    }

    for (Int_t cent = 0; cent < maxCent; cent++){
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

    TGraphAsymmErrors* graphCombEtaToPi0Stat[9]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaToPi0StatWOXErr[9]       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaToPi0Sys[9]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaToPi0Tot[9]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaToPi0RelStat[9]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaToPi0RelSys[9]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombEtaToPi0RelTot[9]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

    TGraphAsymmErrors* graphCombIndepMethodsEtaToPi0Stat[9]             = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombIndepMethodsEtaToPi0Sys[9]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombIndepMethodsEtaToPi0Tot[9]              = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombIndepMethodsEtaToPi0RelStat[9]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombIndepMethodsEtaToPi0RelSys[9]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphCombIndepMethodsEtaToPi0RelTot[9]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

    TGraphAsymmErrors* graphIndEtaToPi0Stat[9][11];
    TGraphAsymmErrors* graphIndEtaToPi0Sys[9][11];
    for (Int_t cent = 0; cent < maxCent; cent++){
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
    histo2DRelTotErrEtaToPi0                 = new TH2F("histo2DRelTotErrEtaToPi0","histo2DRelTotErrEtaToPi0",11000,minPtEtaPlotting, maxPtEtaPlotting,1000,0,70.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErrEtaToPi0, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErrEtaToPi0->GetYaxis()->SetRangeUser(0,39.5);
    histo2DRelTotErrEtaToPi0->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErrEtaToPi0->GetXaxis()->SetLabelOffset(-0.01);
    TLatex *labelRelTotErrEtaToPi0       = new TLatex(0.95,0.85,"#eta/#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRelTotErrEtaToPi0, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);

    for (Int_t cent = 0; cent < maxCent; cent++){
        if (!enableCentComb[cent]) continue;
        // Declaration & calculation of combined spectrum
        cout << "start combining eta for " << centArray[cent].Data()  << endl;
        TString fileNameEtaToPi0OutputWeighting         = Form("%s/EtaToPi0_WeightingMethod_%s%s.dat",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data());
        graphCombEtaToPi0Tot[cent]                      = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionEtaToPi0[cent],    sysErrorCollectionEtaToPi0[cent],
                                                                                                xPtLimitsEtaToPi0[cent], maxNBinsEtaToPi0[cent],
                                                                                                offSetsEtaToPi0[cent], offSetsEtaToPi0Sys[cent],
                                                                                                graphCombEtaToPi0Stat[cent], graphCombEtaToPi0Sys[cent],
                                                                                                fileNameEtaToPi0OutputWeighting, "PbPb_5.02TeV", "EtaToPi0", kTRUE,
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
        for (Int_t cent = 0; cent < maxCent; cent++){
            if (!enableCentComb[cent]) continue;
            if (graphCombEtaToPi0Stat[cent]){
                graphCombEtaToPi0StatWOXErr[cent]           = (TGraphAsymmErrors*)graphCombEtaToPi0Stat[cent]->Clone(Form("EtaToPi0StatW0XErr%s", centArray[cent].Data()));
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

        // Declaration & calculation of combined spectrum with only fully independent techniques
        TString fileNameEtaToPi0OutputWeightingIndepMethod      = Form("%s/EtaToPi0_WeightingMethodIndepMethod_%s%s.dat",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data());
        graphCombIndepMethodsEtaToPi0Tot[cent]          = CombinePtPointsSpectraFullCorrMat(    statErrorCollectionOnlyIndepMethodsEtaToPi0[cent],    sysErrorCollectionOnlyIndepMethodsEtaToPi0[cent],
                                                                                                   xPtLimitsEtaToPi0[cent], maxNBinsEtaToPi0[cent],
                                                                                                   offSetsEtaToPi0[cent], offSetsEtaToPi0Sys[cent],
                                                                                                   graphCombIndepMethodsEtaToPi0Stat[cent], graphCombIndepMethodsEtaToPi0Sys[cent],
                                                                                                   fileNameEtaToPi0OutputWeightingIndepMethod, "PbPb_5.02TeV", "EtaToPi0", kTRUE,
                                                                                                   NULL, fileNameCorrFactors, centArrayCorr[cent]
        );


        if (graphCombIndepMethodsEtaToPi0Tot[cent] == NULL) {
            cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
            return;
        }
        while (graphCombIndepMethodsEtaToPi0Stat[cent]->GetX()[0] < minPtEtaToPi0[cent]){
            graphCombIndepMethodsEtaToPi0Stat[cent]->RemovePoint(0);
        }
        while (graphCombIndepMethodsEtaToPi0Tot[cent]->GetX()[0] < minPtEtaToPi0[cent]){
            graphCombIndepMethodsEtaToPi0Tot[cent]->RemovePoint(0);
        }
        while (graphCombIndepMethodsEtaToPi0Sys[cent]->GetX()[0] < minPtEtaToPi0[cent]){
            graphCombIndepMethodsEtaToPi0Sys[cent]->RemovePoint(0);
        }

        // Reading weights from output file for plotting
        ifstream fileWeightsEtaToPi0ReadIndepMethods;
        fileWeightsEtaToPi0ReadIndepMethods.open(fileNameEtaToPi0OutputWeightingIndepMethod,ios_base::in);
        cout << "reading" << fileNameEtaToPi0OutputWeightingIndepMethod << endl;
        Double_t xValuesEtaToPi0ReadIndepMethods[100];
        Double_t weightsEtaToPi0ReadIndepMethods[11][100];
        Int_t availableEtaToPi0MeasIndepMethods[11]    = {   -1, -1, -1, -1, -1,
                                                        -1, -1, -1, -1, -1,
                                                        -1};
        Int_t nMeasSetEtaToPi0IndepMethods       = 3;
        Int_t nPtBinsEtaToPi0ReadIndepMethods    = 0;
        while(!fileWeightsEtaToPi0ReadIndepMethods.eof() && nPtBinsEtaToPi0ReadIndepMethods < 100){
            TString garbage             = "";
            if (nPtBinsEtaToPi0ReadIndepMethods == 0){
                fileWeightsEtaToPi0ReadIndepMethods >> garbage ;//>> availableEtaToPi0MeasIndepMethods[0] >> availableEtaToPi0MeasIndepMethods[1] >> availableEtaToPi0MeasIndepMethods[2] >> availableEtaToPi0MeasIndepMethods[3];
                for (Int_t i = 0; i < nMeasSetEtaToPi0IndepMethods; i++){
                    fileWeightsEtaToPi0ReadIndepMethods >> availableEtaToPi0MeasIndepMethods[i] ;
                }
                cout << "read following measurements: ";
                for (Int_t i = 0; i < nMeasSetEtaToPi0IndepMethods; i++){
                    cout << availableEtaToPi0MeasIndepMethods[i] << "\t" ;
                }
                cout << endl;
            } else {
                fileWeightsEtaToPi0ReadIndepMethods >> xValuesEtaToPi0ReadIndepMethods[nPtBinsEtaToPi0ReadIndepMethods-1];
                for (Int_t i = 0; i < nMeasSetEtaToPi0IndepMethods; i++){
                    fileWeightsEtaToPi0ReadIndepMethods >> weightsEtaToPi0ReadIndepMethods[availableEtaToPi0MeasIndepMethods[i]][nPtBinsEtaToPi0ReadIndepMethods-1] ;
                }
                cout << "read: "<<  nPtBinsEtaToPi0ReadIndepMethods << "\t"<< xValuesEtaToPi0ReadIndepMethods[nPtBinsEtaToPi0ReadIndepMethods-1] << "\t" ;
                for (Int_t i = 0; i < nMeasSetEtaToPi0IndepMethods; i++){
                    cout << weightsEtaToPi0ReadIndepMethods[availableEtaToPi0MeasIndepMethods[i]][nPtBinsEtaToPi0ReadIndepMethods-1] << "\t";
                }
                cout << endl;
            }
            nPtBinsEtaToPi0ReadIndepMethods++;
        }
        nPtBinsEtaToPi0ReadIndepMethods                  = nPtBinsEtaToPi0ReadIndepMethods-2 ;
        fileWeightsEtaToPi0ReadIndepMethods.close();

        for (Int_t i = 0; i < nMeasSetEtaToPi0IndepMethods; i++){
            graphWeightsIndepMethodsEtaToPi0[cent][availableEtaToPi0MeasIndepMethods[i]]                        = new TGraph(nPtBinsEtaToPi0ReadIndepMethods,xValuesEtaToPi0ReadIndepMethods,weightsEtaToPi0ReadIndepMethods[availableEtaToPi0MeasIndepMethods[i]]);
            Int_t bin = 0;
            for (Int_t n = 0; n< nPtBinsEtaToPi0ReadIndepMethods; n++){
                if (graphWeightsIndepMethodsEtaToPi0[cent][availableEtaToPi0MeasIndepMethods[i]]->GetY()[bin] == 0) graphWeightsIndepMethodsEtaToPi0[cent][availableEtaToPi0MeasIndepMethods[i]]->RemovePoint(bin);
                else bin++;
            }
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

            labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            labelWeightsEnergy->Draw();
            labelWeightsEtaToPi0->Draw();

            DrawGammaLines(minPtEtaToPi0Plotting, maxPtEtaToPi0Plotting , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(minPtEtaToPi0Plotting, maxPtEtaToPi0Plotting , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(minPtEtaToPi0Plotting, maxPtEtaToPi0Plotting , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(minPtEtaToPi0Plotting, maxPtEtaToPi0Plotting , 0.2, 0.2,0.1, kGray, 3);

        canvasWeights->SaveAs(Form("%s/EtaToPi0_Weights_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(),suffix.Data()));

        canvasWeights->cd();
        histo2DEtaToPi0Weights->Draw("copy");

            TLegend* legendWeightsOnlyIndepMethods   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEtaToPi0IndepMethods), 32);
            for (Int_t i = 0; i < nMeasSetEtaToPi0IndepMethods; i++){
                DrawGammaSetMarkerTGraph(graphWeightsIndepMethodsEtaToPi0[cent][availableEtaToPi0MeasIndepMethods[i]], markerStyleDet[availableEtaToPi0MeasIndepMethods[i]], markerSizeDet[availableEtaToPi0MeasIndepMethods[i]]*0.5, colorDet[availableEtaToPi0MeasIndepMethods[i]] , colorDet[availableEtaToPi0MeasIndepMethods[i]]);
                graphWeightsIndepMethodsEtaToPi0[cent][availableEtaToPi0MeasIndepMethods[i]]->Draw("p,same,z");
                legendWeightsOnlyIndepMethods->AddEntry(graphWeightsIndepMethodsEtaToPi0[cent][availableEtaToPi0MeasIndepMethods[i]],nameMeasGlobalLabel[availableEtaToPi0MeasIndepMethods[i]],"p");
            }
            legendWeightsOnlyIndepMethods->Draw();

            labelWeightsEtaToPi0->Draw();

            DrawGammaLines(minPtEtaToPi0Plotting, maxPtEtaToPi0Plotting , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(minPtEtaToPi0Plotting, maxPtEtaToPi0Plotting , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(minPtEtaToPi0Plotting, maxPtEtaToPi0Plotting , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(minPtEtaToPi0Plotting, maxPtEtaToPi0Plotting , 0.2, 0.2,0.1, kGray, 3);

        canvasWeights->SaveAs(Form("%s/EtaToPi0_WeightsOnlyIndepMethods_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

        //  *********************************************************************************************************************
        //  ************************************ Visualize relative errors ******************************************************
        //  *********************************************************************************************************************

        canvasRelSysErr->cd();
        histo2DRelSysErr->GetYaxis()->SetRangeUser(0,65.5);
        histo2DRelSysErr->Draw("copy");

            TLegend* legendRelSysErrEtaToPi0        = GetAndSetLegend2(0.62, 0.92-(0.04*nMeasSetEtaToPi0), 0.95, 0.92, textSizeLabelsPixel);
            for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                DrawGammaSetMarkerTGraph(sysErrorRelCollectionEtaToPi0[cent][availableEtaToPi0Meas[i]], markerStyleDet[availableEtaToPi0Meas[i]], markerSizeDet[availableEtaToPi0Meas[i]]*0.5, colorDet[availableEtaToPi0Meas[i]],
                                         colorDet[availableEtaToPi0Meas[i]]);
                sysErrorRelCollectionEtaToPi0[cent][availableEtaToPi0Meas[i]]->Draw("p,same,z");
                legendRelSysErrEtaToPi0->AddEntry(sysErrorRelCollectionEtaToPi0[cent][availableEtaToPi0Meas[i]],nameMeasGlobalLabel[availableEtaToPi0Meas[i]],"p");
            }
            legendRelSysErrEtaToPi0->Draw();

            labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            labelRelSysErrEnergy->Draw();
            labelRelSysErrEtaToPi0->Draw();

        canvasRelSysErr->SaveAs(Form("%s/EtaToPi0_RelSysErr_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data() ,suffix.Data()));
        delete legendRelSysErrEtaToPi0;
        //  *********************************************************************************************************************
        //  ************************************ Visualize relative errors ******************************************************
        //  *********************************************************************************************************************

        canvasRelStatErr->cd();
        histo2DRelStatErr->GetYaxis()->SetRangeUser(0,65.5);
        histo2DRelStatErr->Draw("copy");
            TLegend* legendRelStatErr       = GetAndSetLegend2(0.14, 0.92-(0.04*nMeasSetEtaToPi0), 0.45, 0.92, textSizeLabelsPixel);
            for (Int_t i = 0; i < nMeasSetEtaToPi0; i++){
                DrawGammaSetMarkerTGraph(statErrorRelCollectionEtaToPi0[cent][availableEtaToPi0Meas[i]], markerStyleDet[availableEtaToPi0Meas[i]], markerSizeDet[availableEtaToPi0Meas[i]]*0.5, colorDet[availableEtaToPi0Meas[i]],
                                         colorDet[availableEtaToPi0Meas[i]]);
                statErrorRelCollectionEtaToPi0[cent][availableEtaToPi0Meas[i]]->Draw("p,same,z");
                legendRelStatErr->AddEntry(statErrorRelCollectionEtaToPi0[cent][availableEtaToPi0Meas[i]],nameMeasGlobalLabel[availableEtaToPi0Meas[i]],"p");
            }
            legendRelStatErr->Draw();

            labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            labelRelStatErrEnergy->Draw();
            labelRelStatErrEtaToPi0->Draw();

        canvasRelStatErr->SaveAs(Form("%s/EtaToPi0_RelStatErr_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

        //  *********************************************************************************************************************
        //  ************************ Visualize relative total errors of different combination methods EtaToPi0 ***********************
        //  *********************************************************************************************************************
        graphCombEtaToPi0RelStat[cent]   = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0Stat[cent], "relativeStatErrorEtaToPi0_Method");
        graphCombEtaToPi0RelSys[cent]    = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0Sys[cent], "relativeSysErrorEtaToPi0_Method");
        graphCombEtaToPi0RelTot[cent]    = CalculateRelErrUpAsymmGraph( graphCombEtaToPi0Tot[cent], "relativeTotalErrorEtaToPi0_Method");
        graphCombIndepMethodsEtaToPi0RelStat[cent]   = CalculateRelErrUpAsymmGraph( graphCombIndepMethodsEtaToPi0Stat[cent], "relativeStatErrorEtaToPi0_IndepMethod");
        graphCombIndepMethodsEtaToPi0RelSys[cent]    = CalculateRelErrUpAsymmGraph( graphCombIndepMethodsEtaToPi0Sys[cent], "relativeSysErrorEtaToPi0_IndepMethod");
        graphCombIndepMethodsEtaToPi0RelTot[cent]    = CalculateRelErrUpAsymmGraph( graphCombIndepMethodsEtaToPi0Tot[cent], "relativeTotalErrorEtaToPi0_IndepMethod");

        canvasRelTotErr->cd();
        histo2DRelTotErrEtaToPi0->GetYaxis()->SetRangeUser(0,65.5);
        histo2DRelTotErrEtaToPi0->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0RelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
            graphCombEtaToPi0RelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombIndepMethodsEtaToPi0RelTot[cent], markerStyleComb+4, markerSizeComb, colorCombHighPt , colorCombHighPt);
            graphCombIndepMethodsEtaToPi0RelTot[cent]->Draw("p,same,z");

            TLegend* legendRelTotErr1       = GetAndSetLegend2(0.20, 0.92-(0.035*2), 0.45, 0.92, 32);
            legendRelTotErr1->AddEntry(graphCombEtaToPi0RelTot[cent],"All methods","p");
            legendRelTotErr1->AddEntry(graphCombIndepMethodsEtaToPi0RelTot[cent],"Independent methods","p");
            legendRelTotErr1->Draw();


            labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            labelRelTotErrEnergy->Draw();
            labelRelTotErrEtaToPi0->Draw();

        canvasRelTotErr->SaveAs(Form("%s/EtaToPi0_TotErr_Comp_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

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

        canvasRelTotErr->SaveAs(Form("%s/EtaToPi0_Reldecomp_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
        histo2DRelTotErrEtaToPi0->GetYaxis()->SetTitle("Err (%)");
        histo2DRelTotErrEtaToPi0->Draw("copy");

            DrawGammaSetMarkerTGraphAsym(graphCombIndepMethodsEtaToPi0RelTot[cent], markerStyleComb, markerSizeComb, colorCombHighPt , colorCombHighPt);
            graphCombIndepMethodsEtaToPi0RelTot[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTGraphAsym(graphCombIndepMethodsEtaToPi0RelStat[cent], markerStyleComb, markerSizeComb, colorCombHighPt-6 , colorCombHighPt-6);
            graphCombIndepMethodsEtaToPi0RelStat[cent]->Draw("l,x0,same,e1");
            DrawGammaSetMarkerTGraphAsym(graphCombIndepMethodsEtaToPi0RelSys[cent], markerStyleComb, markerSizeComb, colorCombHighPt+2, colorCombHighPt+2);
            graphCombIndepMethodsEtaToPi0RelSys[cent]->SetLineStyle(7);
            graphCombIndepMethodsEtaToPi0RelSys[cent]->Draw("l,x0,same,e1");

            TLegend* legendRelTotErr7       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.35, 0.92, 32);
            legendRelTotErr7->AddEntry(graphCombIndepMethodsEtaToPi0RelTot[cent],"tot","p");
            legendRelTotErr7->AddEntry(graphCombIndepMethodsEtaToPi0RelStat[cent],"stat","l");
            legendRelTotErr7->AddEntry(graphCombIndepMethodsEtaToPi0RelSys[cent],"sys","l");
            legendRelTotErr7->Draw();

            labelRelTotErrEnergy->Draw();
            labelRelTotErrEtaToPi0->Draw();

        canvasRelTotErr->SaveAs(Form("%s/EtaToPi0_ReldecompOnlyIndepMethods_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
    }

    TGraphAsymmErrors* graphRAAIndStatPi0[9][11];
    TGraphAsymmErrors* graphRAAIndStatPi0WOXErr[9][11];
    TH1D* histRAAIndStatPi0[9][11];
    TGraphAsymmErrors* graphRAAIndSystPi0[9][11];
    TGraphAsymmErrors* graphRAAIndCombPi0[9][11];
    TGraph* graphWeightsPi0RAA[9][11];
    TGraphAsymmErrors* statErrorRelCollectionPi0RAA[9][11];
    TGraphAsymmErrors* sysErrorRelCollectionPi0RAA[9][11];
    TGraphAsymmErrors* graphRAACombStatPi0[9];
    TGraphAsymmErrors* graphRAACombSystPi0[9];
    TGraphAsymmErrors* graphRAACombStatPi0WOXErr[9];
    TGraphAsymmErrors* graphRAACombCombPi0[9];
    TGraphAsymmErrors* graphCombPi0RAARelStat[9];
    TGraphAsymmErrors* graphCombPi0RAARelSys[9];
    TGraphAsymmErrors* graphCombPi0RAARelTot[9];

    // *******************************************************************************************************
    // *************************** RAA calculation for pi0 for individual measurements **********************
    // *******************************************************************************************************
    cout << "starting RAA calculation for pi0" << endl;

    for (Int_t cent = 0; cent < maxCent; cent++){
        if ( !enableCentRAA[cent] ) continue;
        // *******************************************************************************************************
        // *************************** RAA calculation for pi0 for individual measurements **********************
        // *******************************************************************************************************
        for (Int_t meth = 0; meth < 11; meth++){
            cout << "RAA calc for pi0 " << centArray[cent].Data() << addCentString[cent].Data() << "\t" << "\t" << nameMeasGlobal[meth].Data() << endl;
            graphWeightsPi0RAA[cent][meth]                              = NULL;
            graphRAAIndStatPi0[cent][meth]                              = NULL;
            graphRAAIndStatPi0WOXErr[cent][meth]                        = NULL;
            histRAAIndStatPi0[cent][meth]                               = NULL;
            graphRAAIndSystPi0[cent][meth]                              = NULL;
            graphRAAIndCombPi0[cent][meth]                              = NULL;
            statErrorRelCollectionPi0RAA[cent][meth]                     = NULL;
            sysErrorRelCollectionPi0RAA[cent][meth]                      = NULL;
            if (haveRefPPPi0[meth]  ){ //&& !( i == 3 )
                if(bWCorrection.Contains("Y") ){
                    TGraphAsymmErrors* tempPPStat               = (TGraphAsymmErrors*)statErrorCollectionPi0PP[meth]->Clone("tempStat");
                    TGraphAsymmErrors* tempPPSys                = (TGraphAsymmErrors*)systErrorUnCorrCollectionPi0PP[meth]->Clone("tempSys");
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
                    cout << "stat PbPb" << endl;
                    graphIndPi0InvYieldStat_yShifted[cent][meth]->Print();
                    cout << "sys PbPb" << endl;
                    graphIndPi0InvYieldSys_yShifted[cent][meth]->Print();


                    graphRAAIndCombPi0[cent][meth]                      = CalcRAAV2(  tempPPStat,
                                                                                        tempPPSys,
                                                                                        graphIndPi0InvYieldStat_yShifted[cent][meth],
                                                                                        graphIndPi0InvYieldSys_yShifted[cent][meth],
                                                                                        &graphRAAIndStatPi0[cent][meth],
                                                                                        &graphRAAIndSystPi0[cent][meth],
                                                                                        tPbPb[cent], tPbPbErr[cent] ,
                                                                                        ptSysRemNames[cent][meth], fileNamesPbPbPi0DetailedSys[cent][meth], fileNamesppPi0DetailedSys[cent][meth], fileNamesRAAPi0DetailedSys[cent][meth]
                    );
                    graphRAAIndCombPi0[cent][meth]->Print();
                    histRAAIndStatPi0[cent][meth]                       = (TH1D*)statErrorCollectionPi0[cent][meth]->Clone(Form("histRAA%sStatPi0",nameMeasGlobalLabel[meth].Data()));
                    for (Int_t j = 0; j< histRAAIndStatPi0[cent][meth]->GetNbinsX()+1; j++){
                        histRAAIndStatPi0[cent][meth]->SetBinContent(j,0);
                        histRAAIndStatPi0[cent][meth]->SetBinError(j,0);
                    }
                    for (Int_t j = 0; j < graphRAAIndStatPi0[cent][meth]->GetN(); j++){
                        Int_t bin                               = histRAAIndStatPi0[cent][meth]->GetXaxis()->FindBin(graphRAAIndStatPi0[cent][meth]->GetX()[j]);
                        histRAAIndStatPi0[cent][meth]->SetBinContent(bin,graphRAAIndStatPi0[cent][meth]->GetY()[j]);
                        histRAAIndStatPi0[cent][meth]->SetBinError(bin,graphRAAIndStatPi0[cent][meth]->GetEYlow()[j]);
                    }
                }

                if (histRAAIndStatPi0[cent][meth]){
                    statErrorRelCollectionPi0RAA[cent][meth]   = (TGraphAsymmErrors*)graphRAAIndStatPi0[cent][meth]->Clone(Form("relativeStatErrorPi0RAA%i_%s", cent, nameMeasGlobal[meth].Data()));
                    RemoveZerosAtBeginningAndEndFromGraph(statErrorRelCollectionPi0RAA[cent][meth]);
                    statErrorRelCollectionPi0RAA[cent][meth]   = CalculateRelErrUpAsymmGraph( statErrorRelCollectionPi0RAA[cent][meth], Form("relativeStatErrorPi0RAA_%i_%s", cent,nameMeasGlobal[meth].Data()));
                    statErrorRelCollectionPi0RAA[cent][meth]->Print();
                    graphRAAIndStatPi0WOXErr[cent][meth]      = (TGraphAsymmErrors*)graphRAAIndStatPi0[cent][meth]->Clone(Form("graphRAA%sStatPi0WOXErr%i", nameMeasGlobalLabel[meth].Data(), cent));

                }
                if (graphRAAIndSystPi0[cent][meth]){
                    sysErrorRelCollectionPi0RAA[cent][meth]    = (TGraphAsymmErrors*)graphRAAIndSystPi0[cent][meth]->Clone(Form("relativeSysErrorPi0RAA%i_%s", cent, nameMeasGlobal[meth].Data()));
                    RemoveZerosAtBeginningAndEndFromGraph(sysErrorRelCollectionPi0RAA[cent][meth]);
                    sysErrorRelCollectionPi0RAA[cent][meth]    = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionPi0RAA[cent][meth], Form("relativeSysErrorPi0RAA%i_%s", cent, nameMeasGlobal[meth].Data()));
                }
            }
        }

        // *******************************************************************************************************
        // *************************** RAA calculation for pi0 **************************************************
        // *******************************************************************************************************
        graphRAACombStatPi0[cent]                      = NULL;
        graphRAACombStatPi0[cent]                      = NULL;
        graphRAACombStatPi0[cent]                      = NULL;

        if (!enableCentComb[cent]) continue;
        if(bWCorrection.Contains("Y") ){
            TString fileNamePi0RAAOutputWeighting      = Form("%s/Pi0RAA_WeightingMethod_%s%s.dat", outputDirSupportComb.Data(), centArrayOutput[cent].Data(),
                                                                addCentString[cent].Data());
            graphRAACombCombPi0[cent]                  = CombinePtPointsSpectraFullCorrMat(    histRAAIndStatPi0[cent],    graphRAAIndSystPi0[cent],
                                                                                                xPtLimitsPi0[cent], maxNBinsPi0[cent],
                                                                                                offSetsPi0[cent], offSetsPi0Sys[cent],
                                                                                                graphRAACombStatPi0[cent], graphRAACombSystPi0[cent],
                                                                                                fileNamePi0RAAOutputWeighting, "PbPb_5.02TeV", "Pi0RAA", kTRUE,
                                                                                                NULL, fileNameCorrFactors, centArrayCorr[cent]
                                                                                            );
            if (graphRAACombCombPi0 == NULL) {
                cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
                return;
            }
            while (graphRAACombStatPi0[cent]->GetX()[0] < minPtPi0[cent]){
                graphRAACombStatPi0[cent]->RemovePoint(0);
            }
            while (graphRAACombSystPi0[cent]->GetX()[0] < minPtPi0[cent]){
                graphRAACombSystPi0[cent]->RemovePoint(0);
            }
            while (graphRAACombCombPi0[cent]->GetX()[0] < minPtPi0[cent]){
                graphRAACombCombPi0[cent]->RemovePoint(0);
            }
            graphRAACombCombPi0[cent]->Print();

            graphRAACombStatPi0WOXErr[cent]            = (TGraphAsymmErrors*)graphRAACombStatPi0[cent]->Clone("graphRAACombStatPi0WOXErr");
            ProduceGraphAsymmWithoutXErrors(graphRAACombStatPi0WOXErr[cent]);


            // Reading weights from output file for plotting
            ifstream fileWeightsPi0RAARead;
            fileWeightsPi0RAARead.open(fileNamePi0RAAOutputWeighting,ios_base::in);
            cout << "reading" << fileNamePi0RAAOutputWeighting << endl;
            Double_t xValuesPi0RAARead[50];
            Double_t weightsPi0RAARead[11][50];
            Int_t availablePi0RAAMeas[11]      = { -1, -1, -1, -1, -1,
                                                -1, -1, -1, -1, -1,
                                                -1};
            Int_t nMeasSetPi0RAA               = nTotMeasPi0[cent];
            Int_t nPtBinsPi0RAARead            = 0;
            while(!fileWeightsPi0RAARead.eof() && nPtBinsPi0RAARead < 50){
                TString garbage             = "";
                if (nPtBinsPi0RAARead == 0){
                    fileWeightsPi0RAARead >> garbage ;//>> availablePi0RAAMeas[0] >> availablePi0RAAMeas[1] >> availablePi0RAAMeas[2] >> availablePi0RAAMeas[3];
                    for (Int_t i = 0; i < nMeasSetPi0RAA; i++){
                        fileWeightsPi0RAARead >> availablePi0RAAMeas[i] ;
                    }
                    cout << "read following measurements: ";
                    for (Int_t meth= 0; meth < 11; meth++){
                        cout << availablePi0RAAMeas[meth] << "\t" ;
                    }
                    cout << endl;
                } else {
                    fileWeightsPi0RAARead >> xValuesPi0RAARead[nPtBinsPi0RAARead-1];
                    for (Int_t i = 0; i < nMeasSetPi0RAA; i++){
                        fileWeightsPi0RAARead >> weightsPi0RAARead[availablePi0RAAMeas[i]][nPtBinsPi0RAARead-1] ;
                    }
                    cout << "read: "<<  nPtBinsPi0RAARead << "\t"<< xValuesPi0RAARead[nPtBinsPi0RAARead-1] << "\t" ;
                    for (Int_t i = 0; i < nMeasSetPi0RAA; i++){
                        cout << weightsPi0RAARead[availablePi0RAAMeas[i]][nPtBinsPi0RAARead-1] << "\t";
                    }
                    cout << endl;
                }
                nPtBinsPi0RAARead++;
            }
            nPtBinsPi0RAARead                  = nPtBinsPi0RAARead-2 ;
            fileWeightsPi0RAARead.close();

            for (Int_t i = 0; i < nMeasSetPi0RAA; i++){
                graphWeightsPi0RAA[cent][availablePi0RAAMeas[i]]      = new TGraph(nPtBinsPi0RAARead,xValuesPi0RAARead,weightsPi0RAARead[availablePi0RAAMeas[i]]);
                Int_t bin = 0;
                for (Int_t n = 0; n< nPtBinsPi0RAARead; n++){
                    if (graphWeightsPi0RAA[cent][availablePi0RAAMeas[i]]->GetY()[bin] == 0) graphWeightsPi0RAA[cent][availablePi0RAAMeas[i]]->RemovePoint(bin);
                    else bin++;
                }
            }

            // **********************************************************************************************************************
            // ******************************************* Plotting weights method only EMC *****************************************
            // **********************************************************************************************************************
            textSizeLabelsPixel           = 900*0.04;
            canvasWeights->cd();

            histo2DPi0Weights->Draw("copy");
                TLegend* legendWeightsPi0RAA   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetPi0RAA+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
                for (Int_t i = 0; i < nMeasSetPi0RAA; i++){
                    DrawGammaSetMarkerTGraph(graphWeightsPi0RAA[cent][availablePi0RAAMeas[i]], markerStyleDet[availablePi0RAAMeas[i]], markerSizeDet[availablePi0RAAMeas[i]]*0.5, colorDet[availablePi0RAAMeas[i]] , colorDet[availablePi0RAAMeas[i]]);
                    graphWeightsPi0RAA[cent][availablePi0RAAMeas[i]]->Draw("p,same,z");
                    legendWeightsPi0RAA->AddEntry(graphWeightsPi0RAA[cent][availablePi0RAAMeas[i]],nameMeasGlobalLabel[availablePi0RAAMeas[i]],"p");
                }
                legendWeightsPi0RAA->Draw();

                labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
                labelWeightsEnergy->Draw();
                labelWeightsPi0RAA->Draw();

                DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
                DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
                DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
                DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);

            canvasWeights->SaveAs(Form("%s/Pi0RAA_Weights_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data() ,suffix.Data()));

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelSysErr->cd();

                histo2DRelSysErr->GetYaxis()->SetRangeUser(0,49.5);
                histo2DRelSysErr->Draw("copy");
                TLegend* legendRelSysErrRAA       = GetAndSetLegend2(0.60, 0.92-(0.04*(nMeasSetPi0RAA)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
                for (Int_t i = 0; i < nMeasSetPi0RAA; i++){
                    if (!sysErrorRelCollectionPi0RAA[cent][availablePi0RAAMeas[i]]) continue;
                    DrawGammaSetMarkerTGraph(sysErrorRelCollectionPi0RAA[cent][availablePi0RAAMeas[i]], markerStyleDet[availablePi0RAAMeas[i]], markerSizeDet[availablePi0RAAMeas[i]]*0.5,
                                            colorDet[availablePi0RAAMeas[i]], colorDet[availablePi0RAAMeas[i]]);
                    sysErrorRelCollectionPi0RAA[cent][availablePi0RAAMeas[i]]->Draw("p,same,z");
                    legendRelSysErrRAA->AddEntry(sysErrorRelCollectionPi0RAA[cent][availablePi0RAAMeas[i]],nameMeasGlobalLabel[availablePi0RAAMeas[i]],"p");
                }
                legendRelSysErrRAA->Draw();

                labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
                labelRelSysErrEnergy->Draw();
                labelRelSysErrPi0RAA->Draw();

            canvasRelSysErr->SaveAs(Form("%s/Pi0RAA_RelSysErr_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
            delete legendRelSysErrRAA;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelStatErr->cd();

                histo2DRelStatErr->GetYaxis()->SetRangeUser(0,49.5);
                histo2DRelStatErr->Draw("copy");
                TLegend* legendRelStatErrRAA       = GetAndSetLegend2(0.12, 0.92-(0.04*(nMeasSetPi0RAA)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
                for (Int_t i = 0; i < nMeasSetPi0RAA; i++){
                    if (!statErrorRelCollectionPi0RAA[cent][availablePi0RAAMeas[i]]) continue;
                    DrawGammaSetMarkerTGraph(statErrorRelCollectionPi0RAA[cent][availablePi0RAAMeas[i]], markerStyleDet[availablePi0RAAMeas[i]], markerSizeDet[availablePi0RAAMeas[i]]*0.5,
                                            colorDet[availablePi0RAAMeas[i]], colorDet[availablePi0RAAMeas[i]]);
                    statErrorRelCollectionPi0RAA[cent][availablePi0RAAMeas[i]]->Draw("p,same,z");
                    legendRelStatErrRAA->AddEntry(statErrorRelCollectionPi0RAA[cent][availablePi0RAAMeas[i]],nameMeasGlobalLabel[availablePi0RAAMeas[i]],"p");
                }
                legendRelStatErrRAA->Draw();

                labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
                labelRelStatErrEnergy->Draw();
                labelRelStatErrPi0RAA->Draw();

            canvasRelStatErr->SaveAs(Form("%s/Pi0RAA_RelStatErr_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
            delete legendRelStatErrRAA;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            graphCombPi0RAARelStat[cent]      = CalculateRelErrUpAsymmGraph( graphRAACombStatPi0[cent], "relativeStatErrorPi0RAA_Method");
            graphCombPi0RAARelSys[cent]       = CalculateRelErrUpAsymmGraph( graphRAACombSystPi0[cent], "relativeSysErrorPi0RAA_Method");
            graphCombPi0RAARelTot[cent]       = CalculateRelErrUpAsymmGraph( graphRAACombCombPi0[cent], "relativeTotalErrorPi0RAA_Method");

            canvasRelTotErr->cd();

                histo2DRelTotErrPi0->GetYaxis()->SetRangeUser(0,49.5);
                histo2DRelTotErrPi0->Draw("copy");
                DrawGammaSetMarkerTGraphAsym(graphCombPi0RAARelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
                graphCombPi0RAARelTot[cent]->Draw("p,same,z");
                DrawGammaSetMarkerTGraphAsym(graphCombPi0RAARelStat[cent], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
                graphCombPi0RAARelStat[cent]->Draw("l,x0,same,e1");
                DrawGammaSetMarkerTGraphAsym(graphCombPi0RAARelSys[cent], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
                graphCombPi0RAARelSys[cent]->SetLineStyle(7);
                graphCombPi0RAARelSys[cent]->Draw("l,x0,same,e1");

                TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.35, 0.92, 32);
                legendRelTotErr3->AddEntry(graphCombPi0RAARelTot[cent],"tot","p");
                legendRelTotErr3->AddEntry(graphCombPi0RAARelStat[cent],"stat","l");
                legendRelTotErr3->AddEntry(graphCombPi0RAARelSys[cent],"sys","l");
                legendRelTotErr3->Draw();

                labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
                labelRelTotErrEnergy->Draw();
                labelRelTotErrPi0RAA->Draw();

            canvasRelTotErr->SaveAs(Form("%s/Pi0RAA_Reldecomp_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
            delete legendRelTotErr3;
        }
    }

    TGraphAsymmErrors* graphRAAIndStatEta[9][11];
    TGraphAsymmErrors* graphRAAIndStatEtaWOXErr[9][11];
    TH1D* histRAAIndStatEta[9][11];
    TGraphAsymmErrors* graphRAAIndSystEta[9][11];
    TGraphAsymmErrors* graphRAAIndCombEta[9][11];
    TGraph* graphWeightsEtaRAA[9][11];
    TGraphAsymmErrors* statErrorRelCollectionEtaRAA[9][11];
    TGraphAsymmErrors* sysErrorRelCollectionEtaRAA[9][11];
    TGraphAsymmErrors* graphRAACombStatEta[9];
    TGraphAsymmErrors* graphRAACombSystEta[9];
    TGraphAsymmErrors* graphRAACombStatEtaWOXErr[9];
    TGraphAsymmErrors* graphRAACombCombEta[9];
    TGraphAsymmErrors* graphCombEtaRAARelStat[9];
    TGraphAsymmErrors* graphCombEtaRAARelSys[9];
    TGraphAsymmErrors* graphCombEtaRAARelTot[9];

    // *******************************************************************************************************
    // *************************** RAA calculation for eta for individual measurements **********************
    // *******************************************************************************************************
    cout << "starting RAA calculation for eta" << endl;

    for (Int_t cent = 0; cent < maxCent; cent++){
        if ( !enableCentRAA[cent] ) continue;
        // *******************************************************************************************************
        // *************************** RAA calculation for eta for individual measurements **********************
        // *******************************************************************************************************
        for (Int_t meth = 0; meth < 11; meth++){
            cout << "RAA calc for eta " << centArray[cent].Data() << addCentString[cent].Data() << "\t" << nameMeasGlobal[meth].Data() << endl;
            graphWeightsEtaRAA[cent][meth]                              = NULL;
            graphRAAIndStatEta[cent][meth]                              = NULL;
            graphRAAIndStatEtaWOXErr[cent][meth]                        = NULL;
            histRAAIndStatEta[cent][meth]                               = NULL;
            graphRAAIndSystEta[cent][meth]                              = NULL;
            graphRAAIndCombEta[cent][meth]                              = NULL;
            statErrorRelCollectionEtaRAA[cent][meth]                     = NULL;
            sysErrorRelCollectionEtaRAA[cent][meth]                      = NULL;
            if (haveRefPPEta[meth]  ){ //&& !( i == 3 )
                if(bWCorrection.Contains("Y") ){
                    TGraphAsymmErrors* tempPPStat               = (TGraphAsymmErrors*)statErrorCollectionEtaPP[meth]->Clone("tempStat");
                    TGraphAsymmErrors* tempPPSys                = (TGraphAsymmErrors*)systErrorUnCorrCollectionEtaPP[meth]->Clone("tempSys");
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
                    cout << "stat PbPb" << endl;
                    graphIndEtaInvYieldStat_yShifted[cent][meth]->Print();
                    cout << "sys PbPb" << endl;
                    graphIndEtaInvYieldSys_yShifted[cent][meth]->Print();

                    graphRAAIndCombEta[cent][meth]                      = CalcRAAV2(  tempPPStat,
                                                                                        tempPPSys,
                                                                                        graphIndEtaInvYieldStat_yShifted[cent][meth],
                                                                                        graphIndEtaInvYieldSys_yShifted[cent][meth],
                                                                                        &graphRAAIndStatEta[cent][meth],
                                                                                        &graphRAAIndSystEta[cent][meth],
                                                                                        tPbPb[cent], tPbPbErr[cent] ,
                                                                                        ptSysRemNames[cent][meth], fileNamesPbPbEtaDetailedSys[cent][meth], fileNamesppEtaDetailedSys[cent][meth], fileNamesRAAEtaDetailedSys[cent][meth]
                    );
                    graphRAAIndCombEta[cent][meth]->Print();
                    histRAAIndStatEta[cent][meth]                       = (TH1D*)statErrorCollectionEta[cent][meth]->Clone(Form("histRAA%sStatEta",nameMeasGlobalLabel[meth].Data()));
                    for (Int_t j = 0; j< histRAAIndStatEta[cent][meth]->GetNbinsX()+1; j++){
                        histRAAIndStatEta[cent][meth]->SetBinContent(j,0);
                        histRAAIndStatEta[cent][meth]->SetBinError(j,0);
                    }
                    for (Int_t j = 0; j < graphRAAIndStatEta[cent][meth]->GetN(); j++){
                        Int_t bin                               = histRAAIndStatEta[cent][meth]->GetXaxis()->FindBin(graphRAAIndStatEta[cent][meth]->GetX()[j]);
                        histRAAIndStatEta[cent][meth]->SetBinContent(bin,graphRAAIndStatEta[cent][meth]->GetY()[j]);
                        histRAAIndStatEta[cent][meth]->SetBinError(bin,graphRAAIndStatEta[cent][meth]->GetEYlow()[j]);
                    }
                }

                if (histRAAIndStatEta[cent][meth]){
                    statErrorRelCollectionEtaRAA[cent][meth]   = new TGraphAsymmErrors(histRAAIndStatEta[cent][meth]);
                    RemoveZerosAtBeginningAndEndFromGraph(statErrorRelCollectionEtaRAA[cent][meth]);
                    statErrorRelCollectionEtaRAA[cent][meth]   = CalculateRelErrUpAsymmGraph( statErrorRelCollectionEtaRAA[cent][meth], Form("relativeStatErrorEtaRAA_%i_%s", cent,nameMeasGlobal[meth].Data()));

                    graphRAAIndStatEtaWOXErr[cent][meth]      = (TGraphAsymmErrors*)graphRAAIndStatEta[cent][meth]->Clone(Form("graphRAA%sStatEtaWOXErr%i", nameMeasGlobalLabel[meth].Data(), cent));


                }
                if (graphRAAIndSystEta[cent][meth]){
                    sysErrorRelCollectionEtaRAA[cent][meth]    = (TGraphAsymmErrors*)graphRAAIndSystEta[cent][meth]->Clone(Form("relativeSysErrorEtaRAA%i_%s", cent, nameMeasGlobal[meth].Data()));
                    RemoveZerosAtBeginningAndEndFromGraph(sysErrorRelCollectionEtaRAA[cent][meth]);
                    sysErrorRelCollectionEtaRAA[cent][meth]    = CalculateRelErrUpAsymmGraph( sysErrorRelCollectionEtaRAA[cent][meth], Form("relativeSysErrorEtaRAA%i_%s", cent, nameMeasGlobal[meth].Data()));
                }
            }
        }

        // *******************************************************************************************************
        // *************************** RAA calculation for pi0 **************************************************
        // *******************************************************************************************************
        graphRAACombStatEta[cent]                      = NULL;
        graphRAACombStatEta[cent]                      = NULL;
        graphRAACombStatEta[cent]                      = NULL;

        if (!enableCentComb[cent]) continue;
        if(bWCorrection.Contains("Y") ){
            TString fileNameEtaRAAOutputWeighting      = Form("%s/EtaRAA_WeightingMethod_%s%s.dat", outputDirSupportComb.Data(), centArrayOutput[cent].Data(),
                                                                addCentString[cent].Data());
            graphRAACombCombEta[cent]                  = CombinePtPointsSpectraFullCorrMat(    histRAAIndStatEta[cent],    graphRAAIndSystEta[cent],
                                                                                                xPtLimitsEta[cent], maxNBinsEta[cent],
                                                                                                offSetsEta[cent], offSetsEtaSys[cent],
                                                                                                graphRAACombStatEta[cent], graphRAACombSystEta[cent],
                                                                                                fileNameEtaRAAOutputWeighting, "PbPb_5.02TeV", "EtaRAA", kTRUE,
                                                                                                NULL, fileNameCorrFactors, centArrayCorr[cent]
                                                                                            );
            if (graphRAACombCombEta == NULL) {
                cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
                return;
            }
            while (graphRAACombStatEta[cent]->GetX()[0] < minPtEta[cent]){
                graphRAACombStatEta[cent]->RemovePoint(0);
            }
            while (graphRAACombSystEta[cent]->GetX()[0] < minPtEta[cent]){
                graphRAACombSystEta[cent]->RemovePoint(0);
            }
            while (graphRAACombCombEta[cent]->GetX()[0] < minPtEta[cent]){
                graphRAACombCombEta[cent]->RemovePoint(0);
            }
            graphRAACombCombEta[cent]->Print();

            graphRAACombStatEtaWOXErr[cent]            = (TGraphAsymmErrors*)graphRAACombStatEta[cent]->Clone("graphRAACombStatEtaWOXErr");
            ProduceGraphAsymmWithoutXErrors(graphRAACombStatEtaWOXErr[cent]);


            // Reading weights from output file for plotting
            ifstream fileWeightsEtaRAARead;
            fileWeightsEtaRAARead.open(fileNameEtaRAAOutputWeighting,ios_base::in);
            cout << "reading" << fileNameEtaRAAOutputWeighting << endl;
            Double_t xValuesEtaRAARead[50];
            Double_t weightsEtaRAARead[11][50];
            Int_t availableEtaRAAMeas[11]      = { -1, -1, -1, -1, -1,
                                                -1, -1, -1, -1, -1,
                                                -1};
            Int_t nMeasSetEtaRAA               = nTotMeasEta[cent];
            Int_t nPtBinsEtaRAARead            = 0;
            while(!fileWeightsEtaRAARead.eof() && nPtBinsEtaRAARead < 50){
                TString garbage             = "";
                if (nPtBinsEtaRAARead == 0){
                    fileWeightsEtaRAARead >> garbage ;//>> availableEtaRAAMeas[0] >> availableEtaRAAMeas[1] >> availableEtaRAAMeas[2] >> availableEtaRAAMeas[3];
                    for (Int_t i = 0; i < nMeasSetEtaRAA; i++){
                        fileWeightsEtaRAARead >> availableEtaRAAMeas[i] ;
                    }
                    cout << "read following measurements: ";
                    for (Int_t meth= 0; meth < 11; meth++){
                        cout << availableEtaRAAMeas[meth] << "\t" ;
                    }
                    cout << endl;
                } else {
                    fileWeightsEtaRAARead >> xValuesEtaRAARead[nPtBinsEtaRAARead-1];
                    for (Int_t i = 0; i < nMeasSetEtaRAA; i++){
                        fileWeightsEtaRAARead >> weightsEtaRAARead[availableEtaRAAMeas[i]][nPtBinsEtaRAARead-1] ;
                    }
                    cout << "read: "<<  nPtBinsEtaRAARead << "\t"<< xValuesEtaRAARead[nPtBinsEtaRAARead-1] << "\t" ;
                    for (Int_t i = 0; i < nMeasSetEtaRAA; i++){
                        cout << weightsEtaRAARead[availableEtaRAAMeas[i]][nPtBinsEtaRAARead-1] << "\t";
                    }
                    cout << endl;
                }
                nPtBinsEtaRAARead++;
            }
            nPtBinsEtaRAARead                  = nPtBinsEtaRAARead-2 ;
            fileWeightsEtaRAARead.close();

            for (Int_t i = 0; i < nMeasSetEtaRAA; i++){
                graphWeightsEtaRAA[cent][availableEtaRAAMeas[i]]      = new TGraph(nPtBinsEtaRAARead,xValuesEtaRAARead,weightsEtaRAARead[availableEtaRAAMeas[i]]);
                Int_t bin = 0;
                for (Int_t n = 0; n< nPtBinsEtaRAARead; n++){
                    if (graphWeightsEtaRAA[cent][availableEtaRAAMeas[i]]->GetY()[bin] == 0) graphWeightsEtaRAA[cent][availableEtaRAAMeas[i]]->RemovePoint(bin);
                    else bin++;
                }
            }

            // **********************************************************************************************************************
            // ******************************************* Plotting weights method only EMC *****************************************
            // **********************************************************************************************************************
            textSizeLabelsPixel           = 900*0.04;
            canvasWeights->cd();

            histo2DEtaWeights->Draw("copy");
                TLegend* legendWeightsEtaRAA   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.04*(nMeasSetEtaRAA+1)/2), textSizeLabelsPixel, 2, "", 43, 0);
                for (Int_t i = 0; i < nMeasSetEtaRAA; i++){
                    DrawGammaSetMarkerTGraph(graphWeightsEtaRAA[cent][availableEtaRAAMeas[i]], markerStyleDet[availableEtaRAAMeas[i]], markerSizeDet[availableEtaRAAMeas[i]]*0.5, colorDet[availableEtaRAAMeas[i]] , colorDet[availableEtaRAAMeas[i]]);
                    graphWeightsEtaRAA[cent][availableEtaRAAMeas[i]]->Draw("p,same,z");
                    legendWeightsEtaRAA->AddEntry(graphWeightsEtaRAA[cent][availableEtaRAAMeas[i]],nameMeasGlobalLabel[availableEtaRAAMeas[i]],"p");
                }
                legendWeightsEtaRAA->Draw();

                labelWeightsEnergy->SetText(0.95,0.20,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
                labelWeightsEnergy->Draw();
                labelWeightsEtaRAA->Draw();

                DrawGammaLines(0.23, 25. , 0.5, 0.5,1, kGray, 7);
                DrawGammaLines(0.23, 25. , 0.4, 0.4,1, kGray, 1);
                DrawGammaLines(0.23, 25. , 0.3, 0.3,1, kGray, 7);
                DrawGammaLines(0.23, 25. , 0.2, 0.2,1, kGray, 3);

            canvasWeights->SaveAs(Form("%s/EtaRAA_Weights_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data() ,suffix.Data()));

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelSysErr->cd();

                histo2DRelSysErr->GetYaxis()->SetRangeUser(0,49.5);
                histo2DRelSysErr->Draw("copy");
                TLegend* legendRelSysErrRAA       = GetAndSetLegend2(0.60, 0.92-(0.04*(nMeasSetEtaRAA)/2), 0.95, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
                for (Int_t i = 0; i < nMeasSetEtaRAA; i++){
                    if (!sysErrorRelCollectionEtaRAA[cent][availableEtaRAAMeas[i]]) continue;
                    DrawGammaSetMarkerTGraph(sysErrorRelCollectionEtaRAA[cent][availableEtaRAAMeas[i]], markerStyleDet[availableEtaRAAMeas[i]], markerSizeDet[availableEtaRAAMeas[i]]*0.5,
                                            colorDet[availableEtaRAAMeas[i]], colorDet[availableEtaRAAMeas[i]]);
                    sysErrorRelCollectionEtaRAA[cent][availableEtaRAAMeas[i]]->Draw("p,same,z");
                    legendRelSysErrRAA->AddEntry(sysErrorRelCollectionEtaRAA[cent][availableEtaRAAMeas[i]],nameMeasGlobalLabel[availableEtaRAAMeas[i]],"p");
                }
                legendRelSysErrRAA->Draw();

                labelRelSysErrEnergy->SetText(0.15,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
                labelRelSysErrEnergy->Draw();
                labelRelSysErrEtaRAA->Draw();

            canvasRelSysErr->SaveAs(Form("%s/EtaRAA_RelSysErr_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data() ,suffix.Data()));
            delete legendRelSysErrRAA;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            canvasRelStatErr->cd();

                histo2DRelStatErr->GetYaxis()->SetRangeUser(0,49.5);
                histo2DRelStatErr->Draw("copy");
                TLegend* legendRelStatErrRAA       = GetAndSetLegend2(0.12, 0.92-(0.04*(nMeasSetEtaRAA)/2), 0.45, 0.92, textSizeLabelsPixel, 2, "", 43, 0);
                for (Int_t i = 0; i < nMeasSetEtaRAA; i++){
                    if (!statErrorRelCollectionEtaRAA[cent][availableEtaRAAMeas[i]]) continue;
                    DrawGammaSetMarkerTGraph(statErrorRelCollectionEtaRAA[cent][availableEtaRAAMeas[i]], markerStyleDet[availableEtaRAAMeas[i]], markerSizeDet[availableEtaRAAMeas[i]]*0.5,
                                            colorDet[availableEtaRAAMeas[i]], colorDet[availableEtaRAAMeas[i]]);
                    statErrorRelCollectionEtaRAA[cent][availableEtaRAAMeas[i]]->Draw("p,same,z");
                    legendRelStatErrRAA->AddEntry(statErrorRelCollectionEtaRAA[cent][availableEtaRAAMeas[i]],nameMeasGlobalLabel[availableEtaRAAMeas[i]],"p");
                }
                legendRelStatErrRAA->Draw();

                labelRelStatErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
                labelRelStatErrEnergy->Draw();
                labelRelStatErrEtaRAA->Draw();

            canvasRelStatErr->SaveAs(Form("%s/EtaRAA_RelStatErr_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
            delete legendRelStatErrRAA;

            //  *********************************************************************************************************************
            //  ************************************ Visualize relative errors ******************************************************
            //  *********************************************************************************************************************

            graphCombEtaRAARelStat[cent]      = CalculateRelErrUpAsymmGraph( graphRAACombStatEta[cent], "relativeStatErrorEtaRAA_Method");
            graphCombEtaRAARelSys[cent]       = CalculateRelErrUpAsymmGraph( graphRAACombSystEta[cent], "relativeSysErrorEtaRAA_Method");
            graphCombEtaRAARelTot[cent]       = CalculateRelErrUpAsymmGraph( graphRAACombCombEta[cent], "relativeTotalErrorEtaRAA_Method");

            canvasRelTotErr->cd();

                histo2DRelTotErrEta->GetYaxis()->SetRangeUser(0,49.5);
                histo2DRelTotErrEta->Draw("copy");
                DrawGammaSetMarkerTGraphAsym(graphCombEtaRAARelTot[cent], markerStyleComb, markerSizeComb, colorComb , colorComb);
                graphCombEtaRAARelTot[cent]->Draw("p,same,z");
                DrawGammaSetMarkerTGraphAsym(graphCombEtaRAARelStat[cent], markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
                graphCombEtaRAARelStat[cent]->Draw("l,x0,same,e1");
                DrawGammaSetMarkerTGraphAsym(graphCombEtaRAARelSys[cent], markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
                graphCombEtaRAARelSys[cent]->SetLineStyle(7);
                graphCombEtaRAARelSys[cent]->Draw("l,x0,same,e1");

                TLegend* legendRelTotErr3       = GetAndSetLegend2(0.14, 0.92-(0.035*3), 0.35, 0.92, 32);
                legendRelTotErr3->AddEntry(graphCombEtaRAARelTot[cent],"tot","p");
                legendRelTotErr3->AddEntry(graphCombEtaRAARelStat[cent],"stat","l");
                legendRelTotErr3->AddEntry(graphCombEtaRAARelSys[cent],"sys","l");
                legendRelTotErr3->Draw();

                labelRelTotErrEnergy->SetText(0.95,0.89,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
                labelRelTotErrEnergy->Draw();
                labelRelTotErrEtaRAA->Draw();

            canvasRelTotErr->SaveAs(Form("%s/EtaRAA_Reldecomp_%s%s.%s",outputDirSupportComb.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
            delete legendRelTotErr3;
        }
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
    TLatex *labelMassEnergy         = new TLatex(0.13,0.775,collisionSystemPbPb.Data());
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

    for (Int_t cent = 0; cent < maxCent; cent++){
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
            labelMassEnergy->SetText(0.13,0.775,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
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
        canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth_%s%s.%s",outputDirSupportPaper.Data(),centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

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
    TLatex *labelEnergyMass2         = new TLatex(0.15,0.87,collisionSystemPbPb.Data());
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
        for (Int_t cent = 0; cent < maxCent; cent++){
            if (!enableCent[cent]) continue;
            if (graphPi0Mass[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphPi0Mass[cent][meth], markerStyleCent[cent], markerSizeCent[cent]*0.55, colorCent[cent] , colorCent[cent]);
                graphPi0Mass[cent][meth]->Draw("p,same,z");
                legendMass->AddEntry(graphPi0Mass[cent][meth],centArray2[cent],"p");
            }
        }
        legendMass->Draw();
        labelPerfMass2->Draw();
        labelEnergyMass2->SetText(0.15,0.87,Form("%s", collisionSystemPbPb.Data()));
        labelEnergyMass2->Draw();
        labelDetSysPi0Mass->SetText(0.15,0.82,Form("#pi^{0} #rightarrow #gamma#gamma, %s", nameMeasGlobalLabel[meth].Data()));
        labelDetSysPi0Mass->Draw();

        canvasMassWidthPi0->Update();
        canvasMassWidthPi0->Print(Form("%s/Pi0_MassPositionData_%s.%s",outputDirSupport.Data(), nameMeasGlobalLabel[meth].Data(), suffix.Data()));

        for (Int_t cent = 0; cent < maxCent; cent++){
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

    for (Int_t cent = 0; cent < maxCent; cent++){
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
            labelMassEnergy->SetText(0.13,0.775,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
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
        canvasMassWidthEta->Print(Form("%s/Eta_MassAndWidth_%s%s.%s",outputDirSupportPaper.Data(),centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

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
        for (Int_t cent = 0; cent < maxCent; cent++){
            if (!enableCent[cent]) continue;
            if (graphEtaMass[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphEtaMass[cent][meth], markerStyleCent[cent], markerSizeCent[cent]*0.55, colorCent[cent] , colorCent[cent]);
                graphEtaMass[cent][meth]->Draw("p,same,z");
                legendMass->AddEntry(graphEtaMass[cent][meth],centArray2[cent],"p");
            }
        }
        legendMass->Draw();
        labelPerfMass2->Draw();
        labelEnergyMass2->SetText(0.15,0.87,Form("%s", collisionSystemPbPb.Data()));
        labelEnergyMass2->Draw();
        labelDetSysEtaMass->SetText(0.15,0.82,Form("#eta #rightarrow #gamma#gamma, %s", nameMeasGlobalLabel[meth].Data()));
        labelDetSysEtaMass->Draw();

        canvasMassWidthEta->Update();
        canvasMassWidthEta->Print(Form("%s/Eta_MassPositionData_%s.%s",outputDirSupport.Data(), nameMeasGlobalLabel[meth].Data(), suffix.Data()));

        for (Int_t cent = 0; cent < maxCent; cent++){
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

    TH2F * histo2DYieldPi0              = new TH2F("histo2DYieldPi0","histo2DYieldPi0",11000,minPtPi0Plotting, maxPtPi0Plotting,1000,6e-11,3e3);
    SetStyleHistoTH2ForGraphs(histo2DYieldPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.68);
    histo2DYieldPi0->GetXaxis()->SetMoreLogLabels();
    histo2DYieldPi0->GetXaxis()->SetLabelOffset(-0.01);

    TLatex *labelEnergyXSectionPi0  = new TLatex(0.935,0.92,collisionSystemPbPb.Data());
    SetStyleTLatex( labelEnergyXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);
    TLatex *labelALICEXSectionPi0  = new TLatex(0.935,0.882,textALICE.Data());
    SetStyleTLatex( labelALICEXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);
    TLatex *labelDetSysXSectionPi0  = new TLatex(0.935,0.844,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);


    for (Int_t cent = 0; cent < maxCent; cent++){
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
            labelEnergyXSectionPi0->SetText(0.94,0.92,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            labelEnergyXSectionPi0->Draw();
            labelALICEXSectionPi0->Draw();
            labelDetSysXSectionPi0->Draw();


            TLegend* legendXSectionPi0      = GetAndSetLegend2(0.20, 0.13,0.40,0.13+(0.035*nTotMeasPi0[cent]), 0.035, 1, "", 42, 0);
            for (Int_t meth = 0; meth < 11; meth++){
                if (graphIndPi0InvYieldSys[cent][meth]) legendXSectionPi0->AddEntry(graphIndPi0InvYieldSys[cent][meth],nameMeasGlobalLabel[meth],"fp");
            }
            legendXSectionPi0->Draw();

        canvasInvYieldPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems_%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

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

        canvasInvYieldPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems_Comb_%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(),  suffix.Data()));

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

        canvasInvYieldPi0->SaveAs(Form("%s/Pi0_InvYield_CombWithTheory_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
    }

    histo2DYieldPi0->Draw("copy");
    Int_t nCurrEst  = 5;
    TLegend* legendPi0YieldPaper    = GetAndSetLegend2(0.19, 0.10, 0.7, 0.10+0.035*(nCurrEst), textSizeLabelsPixel, 2, "", 43, 0.25);

        for (Int_t cent= 0; cent < maxCent; cent++ ){
            //         for (Int_t cent= 0; cent < maxCent; cent++ ){
            if (!enableCentComb[cent] ) continue;
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

        TLatex *labelEnergyYieldPaper= new TLatex(0.935, 0.92, Form("%s",  collisionSystemPbPb.Data()));
        SetStyleTLatex( labelEnergyYieldPaper, 0.035, 4, 1, 42, kTRUE, 31);
        labelEnergyYieldPaper->Draw();
        TLatex *labelALICEYieldPaper= new TLatex(0.935,0.882,textALICE.Data());
        SetStyleTLatex( labelALICEYieldPaper, 0.035, 4, 1, 42, kTRUE, 31);
        labelALICEYieldPaper->Draw();
        TLatex *labelDetSysYieldPaper= new TLatex(0.935,0.844,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysYieldPaper, 0.035, 4, 1, 42, kTRUE, 31);
        labelDetSysYieldPaper->Draw();
        legendPi0YieldPaper->Draw();
    canvasInvYieldPi0->SaveAs(Form("%s/Pi0_InvYieldAllCent_%s.%s",outputDir.Data(), nameCentEst.Data(), suffix.Data()));
    delete legendPi0YieldPaper;

    histo2DYieldPi0->Draw("copy");
    TLegend* legendPi0YieldPaperWPP    = GetAndSetLegend2(0.19, 0.10, 0.7, 0.10+0.035*(nCurrEst+1), textSizeLabelsPixel, 2, "", 43, 0.25);

    for (Int_t cent= 0; cent < maxCent; cent++ ){
        //         for (Int_t cent= 0; cent < maxCent; cent++ ){
        if (!enableCentComb[cent] ) continue;
        graphCombPi0InvYieldSys[cent]->Draw("E2same");
        graphCombPi0InvYieldStatWOXErr[cent]->Draw("p,same,z");
        fitTCMInvYieldPi0[cent]->Draw("same");
        legendPi0YieldPaperWPP->AddEntry(graphCombPi0InvYieldSys[cent],centArray[cent].Data(),"pf");
        legendPi0YieldPaperWPP->AddEntry(fitTCMInvYieldPi0[cent],"TCM fit","l");

    }

    if (enableCentComb[0] ){
        graphCombPi0InvYieldSys[0]->Draw("E2same");
        graphCombPi0InvYieldStatWOXErr[0]->Draw("p,same,z");
        fitTCMInvYieldPi0[0]->Draw("same");
    }

    if (graphPPInvYieldCombPi0Sys && graphPPInvYieldCombPi0StatW0XErr){
        DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombPi0Sys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphPPInvYieldCombPi0Sys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombPi0StatW0XErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
        graphPPInvYieldCombPi0StatW0XErr->Draw("p,same,z");
        DrawGammaSetMarkerTF1( fitPPTCMInvYieldPi0, 7, 2, kGray+2);
        fitPPTCMInvYieldPi0->Draw("same");
        legendPi0YieldPaperWPP->AddEntry(graphPPInvYieldCombPi0Sys,collisionSystempp.Data(),"pf");
        legendPi0YieldPaperWPP->AddEntry(fitPPTCMInvYieldPi0,"TCM fit","l");
    }
    labelEnergyYieldPaper->Draw();
    labelALICEYieldPaper->Draw();
    labelDetSysYieldPaper->Draw();
    legendPi0YieldPaperWPP->Draw();
    canvasInvYieldPi0->SaveAs(Form("%s/Pi0_InvYieldAllCentwithPP_%s.%s",outputDir.Data(), nameCentEst.Data(), suffix.Data()));
    delete legendPi0YieldPaperWPP;


    // **********************************************************************************************************************
    // ******************************** Cross section for eta single measurement 2.76TeV ************************************
    // **********************************************************************************************************************
    TCanvas* canvasInvYieldEta          = new TCanvas("canvasInvYieldEta","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasInvYieldEta, 0.16, 0.02, 0.02, 0.08);
    canvasInvYieldEta->SetLogx();
    canvasInvYieldEta->SetLogy();

    TH2F * histo2DYieldEta              = new TH2F("histo2DYieldEta","histo2DYieldEta",11000,minPtEtaPlotting, maxPtEtaPlotting,1000,6e-11,3e2);
    SetStyleHistoTH2ForGraphs(histo2DYieldEta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.68);
    histo2DYieldEta->GetXaxis()->SetMoreLogLabels();
    histo2DYieldEta->GetXaxis()->SetLabelOffset(-0.01);

    TLatex *labelEnergyXSectionEta  = new TLatex(0.935,0.92,collisionSystemPbPb.Data());
    SetStyleTLatex( labelEnergyXSectionEta, 0.035,4, 1, 42, kTRUE, 31);
    TLatex *labelALICEXSectionEta   = new TLatex(0.935,0.882,textALICE.Data());
    SetStyleTLatex( labelALICEXSectionEta, 0.035,4, 1, 42, kTRUE, 31);
    TLatex *labelDetSysXSectionEta  = new TLatex(0.935,0.844,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysXSectionEta, 0.035,4, 1, 42, kTRUE, 31);

    for (Int_t cent = 0; cent < maxCent; cent++){
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
        labelEnergyXSectionEta->SetText(0.94,0.92,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
        labelEnergyXSectionEta->Draw();
        labelALICEXSectionEta->Draw();
        labelDetSysXSectionEta->Draw();


        TLegend* legendXSectionEta      = GetAndSetLegend2(0.20, 0.13,0.40,0.13+(0.035*nTotMeasEta[cent]), 0.035, 1, "", 42, 0);
        for (Int_t meth = 0; meth < 11; meth++){
            if (graphIndEtaInvYieldSys[cent][meth]) legendXSectionEta->AddEntry(graphIndEtaInvYieldSys[cent][meth],nameMeasGlobalLabel[meth],"fp");
        }
        legendXSectionEta->Draw();

        canvasInvYieldEta->SaveAs(Form("%s/Eta_InvYieldCompAllSystems_%s%s.%s",outputDirSupportPaper.Data(),centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

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

        canvasInvYieldEta->SaveAs(Form("%s/Eta_InvYieldCompAllSystems_Comb_%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(),
                                       suffix.Data()));

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

        canvasInvYieldEta->SaveAs(Form("%s/Eta_InvYield_CombWithTheory_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
    }
    canvasInvYieldEta->cd();

    histo2DYieldEta->Draw("copy");
    TLegend* legendEtaYieldPaper    = GetAndSetLegend2(0.19, 0.10, 0.7, 0.10+0.035*(nCurrEst), textSizeLabelsPixel, 2, "", 43, 0.25);

    for (Int_t cent= 0; cent < maxCent; cent++ ){
        //         for (Int_t cent= 0; cent < maxCent; cent++ ){
        if (!enableCentComb[cent] ) continue;
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

//     TLatex *labelEnergyYieldPaper= new TLatex(0.935, 0.92, Form("%s",  collisionSystemPbPb.Data()));
    SetStyleTLatex( labelEnergyYieldPaper, 0.035, 4, 1, 42, kTRUE, 31);
    labelEnergyYieldPaper->Draw();
//     TLatex *labelALICEYieldPaper= new TLatex(0.935,0.882,textALICE.Data());
    SetStyleTLatex( labelALICEYieldPaper, 0.035, 4, 1, 42, kTRUE, 31);
    labelALICEYieldPaper->Draw();
    TLatex *labelDetSysYieldPaperEta= new TLatex(0.935,0.844,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysYieldPaperEta, 0.035, 4, 1, 42, kTRUE, 31);
    labelDetSysYieldPaperEta->Draw();
    legendEtaYieldPaper->Draw();
    canvasInvYieldEta->SaveAs(Form("%s/Eta_InvYieldAllCent_%s.%s",outputDir.Data(), nameCentEst.Data(), suffix.Data()));
    delete legendEtaYieldPaper;

    histo2DYieldEta->Draw("copy");
    TLegend* legendEtaYieldPaperWPP    = GetAndSetLegend2(0.19, 0.10, 0.7, 0.10+0.035*(nCurrEst+1), textSizeLabelsPixel, 2, "", 43, 0.25);

    for (Int_t cent= 0; cent < maxCent; cent++ ){
        //         for (Int_t cent= 0; cent < maxCent; cent++ ){
        if (!enableCentComb[cent] ) continue;
        graphCombEtaInvYieldSys[cent]->Draw("E2same");
        graphCombEtaInvYieldStatWOXErr[cent]->Draw("p,same,z");
        fitTCMInvYieldEta[cent]->Draw("same");
        legendEtaYieldPaperWPP->AddEntry(graphCombEtaInvYieldSys[cent],centArray[cent].Data(),"pf");
        legendEtaYieldPaperWPP->AddEntry(fitTCMInvYieldEta[cent],"TCM fit","l");

    }

    if (enableCentComb[0] ){
        graphCombEtaInvYieldSys[0]->Draw("E2same");
        graphCombEtaInvYieldStatWOXErr[0]->Draw("p,same,z");
        fitTCMInvYieldEta[0]->Draw("same");
    }

    if (graphPPInvYieldCombEtaSys && graphPPInvYieldCombEtaStatW0XErr){
        DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombEtaSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphPPInvYieldCombEtaSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombEtaStatW0XErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
        graphPPInvYieldCombEtaStatW0XErr->Draw("p,same,z");
        DrawGammaSetMarkerTF1( fitPPTCMInvYieldEta, 7, 2, kGray+2);
        fitPPTCMInvYieldEta->Draw("same");
        legendEtaYieldPaperWPP->AddEntry(graphPPInvYieldCombEtaSys,collisionSystempp.Data(),"pf");
        legendEtaYieldPaperWPP->AddEntry(fitPPTCMInvYieldEta,"TCM fit","l");
    }
    labelEnergyYieldPaper->Draw();
    labelALICEYieldPaper->Draw();
    labelDetSysYieldPaperEta->Draw();
    legendEtaYieldPaperWPP->Draw();
    canvasInvYieldEta->SaveAs(Form("%s/Eta_InvYieldAllCentwithPP_%s.%s",outputDir.Data(), nameCentEst.Data(), suffix.Data()));
    delete legendEtaYieldPaperWPP;


    // **********************************************************************************************************************
    // ******************************** Cross section for pi0 & eta single measurement **************************************
    // **********************************************************************************************************************
    TH2F * histo2DYieldPi0Eta              = new TH2F("histo2DYieldPi0Eta","histo2DYieldPi0Eta",11000,minPtPi0Plotting, maxPtPi0Plotting,1000,6e-13,3e4);
    SetStyleHistoTH2ForGraphs(histo2DYieldPi0Eta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.68);
    histo2DYieldPi0Eta->GetXaxis()->SetMoreLogLabels();
    histo2DYieldPi0Eta->GetXaxis()->SetLabelOffset(-0.01);

    for (Int_t cent = 0; cent < maxCent; cent++){
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
        labelEnergyXSectionPi0->SetText(0.94,0.92,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
        labelEnergyXSectionPi0->Draw();
        labelALICEXSectionPi0->Draw();

        canvasInvYieldPi0->SaveAs(Form("%s/Pi0AndEta_InvYield_CombWithTheory_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

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
    TLatex *labelEnergyEffi         = new TLatex(0.137,0.87,collisionSystemPbPb.Data());
    SetStyleTLatex( labelEnergyEffi, textSizeLabelsRel,4);
    TLatex *labelDetSysEffiPi0      = new TLatex(0.137,0.82,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelDetSysEffiPi0, textSizeLabelsRel,4);

    for (Int_t cent = 0; cent < maxCent; cent++){
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
            labelEnergyEffi->SetText(0.137,0.87,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
            labelEnergyEffi->Draw();
            labelDetSysEffiPi0->Draw();

        canvasAcceptanceTimesEff->Update();
        canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff_%s%s.%s",outputDirSupportPaper.Data(),centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
    }

    Double_t rangMinEffi[11]                = { 7e-5, 1e-3, 8e-4, 9e-5, 9e-5, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4 };
    Double_t rangMaxEffi[11]                = { 5e-2, 3e-1, 2e-0, 5e-2, 3e-1, 2e-0, 2e-0, 2e-0, 2e-0, 2e-0, 2e-0 };

    for (Int_t meth = 0; meth < 11; meth++){
        if (!fileNamesMethod[meth].CompareTo("")) continue;
        canvasAcceptanceTimesEff->cd();
        histo1DAccEff->GetYaxis()->SetRangeUser(rangMinEffi[meth], rangMaxEffi[meth] );
        histo1DAccEff->DrawCopy();

        TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(5*textSizeLabelsRel),textSizeLabelsPixel);
        for (Int_t cent = 0; cent < maxCent; cent++){
            if (!enableCent[cent]) continue;
            if (graphPi0EffTimesAcc[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphPi0EffTimesAcc[cent][meth], markerStyleCent[cent], markerSizeCent[cent]*0.55, colorCent[cent] , colorCent[cent]);
                graphPi0EffTimesAcc[cent][meth]->Draw("p,same,z");
                legendEffiAccPi0->AddEntry(graphPi0EffTimesAcc[cent][meth],centArray2[cent],"p");
            }
        }
        legendEffiAccPi0->Draw();

        labelPerfEffi->Draw();
        labelEnergyEffi->SetText(0.137,0.87,Form("%s", collisionSystemPbPb.Data()));
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

    for (Int_t cent = 0; cent < maxCent; cent++){
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
        labelEnergyEffi->SetText(0.137,0.87,Form("%s %s", centArray[cent].Data(), collisionSystemPbPb.Data()));
        labelEnergyEffi->Draw();
        labelDetSysEffiEta->Draw();

        canvasAcceptanceTimesEff->Update();
        canvasAcceptanceTimesEff->Print(Form("%s/Eta_AcceptanceTimesEff_%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
    }

    Double_t rangMinEffiEta[11]                = { 9e-4, 1e-3, 8e-4, 9e-4, 9e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4 };
    Double_t rangMaxEffiEta[11]                = { 5e-2, 3e-1, 2e-0, 5e-2, 3e-1, 2e-0, 2e-0, 2e-0, 2e-0, 2e-0, 2e-0 };

    for (Int_t meth = 0; meth < 11; meth++){
        if (!fileNamesMethod[meth].CompareTo("")) continue;
        canvasAcceptanceTimesEff->cd();
        histo1DAccEffEta->GetYaxis()->SetRangeUser(rangMinEffiEta[meth], rangMaxEffiEta[meth] );
        histo1DAccEffEta->DrawCopy();

        TLegend* legendEffiAccEta           = GetAndSetLegend2(0.65, 0.13, 0.93, 0.13+(5*textSizeLabelsRel),textSizeLabelsPixel);
        for (Int_t cent = 0; cent < maxCent; cent++){
            if (!enableCent[cent]) continue;
            if (graphEtaEffTimesAcc[cent][meth]){
                DrawGammaSetMarkerTGraphAsym(graphEtaEffTimesAcc[cent][meth], markerStyleCent[cent], markerSizeCent[cent]*0.55, colorCent[cent] , colorCent[cent]);
                graphEtaEffTimesAcc[cent][meth]->Draw("p,same,z");
                legendEffiAccEta->AddEntry(graphEtaEffTimesAcc[cent][meth],centArray2[cent],"p");
            }
        }
        legendEffiAccEta->Draw();

        labelPerfEffi->Draw();
        labelEnergyEffi->SetText(0.137,0.87,Form("%s", collisionSystemPbPb.Data()));
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
        TLatex *labelEnergySecCorr      = new TLatex(0.15,0.85,collisionSystemPbPb.Data());
        SetStyleTLatex( labelEnergySecCorr, textSizeLabelsRel,4);
        TLatex *labelDetSysSecCorrPi0   = new TLatex(0.15,0.80,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysSecCorrPi0, textSizeLabelsRel,4);

        for (Int_t cent= 0; cent < maxCent; cent++ ){
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
                    canvasEffectiveSecCorr->Print(Form("%s/Pi0_EffectiveSecCorr_%s%s_%s.%s",outputDirSupportPaper.Data(), nameSecPi0SourceRead[k].Data(),
                                                       centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
                }
            }
        }
    delete canvasEffectiveSecCorr;

    fileFitsOutput << "*******************************************************************************************" << endl;
    fileFitsOutput << "****************************** Power law fit pi0 ******************************************" << endl;
    fileFitsOutput << "*******************************************************************************************" << endl;
    TF1* fitPowInvYieldPi0Tot[9]   = { NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL};
    TF1* fitPowInvYieldPi0Stat[9]   = { NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL};
    TF1* fitOHagInvYieldPi0Tot[9]   = { NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL};
    TF1* fitPowInvYieldEtaTot[9]   = { NULL, NULL, NULL, NULL, NULL,   NULL, NULL, NULL, NULL};
    TF1* fitPowInvYieldEtaStat[9]   = { NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL};
    TF1* fitOHagInvYieldEtaTot[9]   = { NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL};

    for (Int_t cent= 0; cent < maxCent; cent++ ){
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

    TBox* boxErrorNormRatioPi0          = CreateBoxConv(kGray+1, 0.28, 1.-(0.031 ), 0.32, 1.+(0.031));
    TBox* boxErrorNormRatioEta          = CreateBoxConv(kGray+1, 0.48, 1.-(0.031 ), 0.52, 1.+(0.031));

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

//     TLegend* legendXsectionPaper    = GetAndSetLegend2(0.19, 0.03, 0.7, 0.03+0.05*(maxCent), textSizeLabelsPixel, 2, "", 43, 0.25);
    TLegend* legendXsectionPaper    = GetAndSetLegend2(0.19, 0.03, 0.7, 0.03+0.05*(maxCent+1), textSizeLabelsPixel, 1, "", 43, 0.15);
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

        for (Int_t cent= 0; cent < maxCent; cent++ ){
            if (!enableCentComb[cent] ) continue;
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
            graphCombPi0InvYieldSys[cent]->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.75, colorCent[cent], colorCent[cent]);
            graphCombPi0InvYieldStatWOXErr[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTF1( fitTCMInvYieldPi0[cent], 7, 2, colorCentMC[cent]);
            fitTCMInvYieldPi0[cent]->Draw("same");
            legendXsectionPaper->AddEntry(graphCombPi0InvYieldSys[cent],centArray2[cent].Data(),"pf");
//             legendXsectionPaper->AddEntry(fitTCMInvYieldPi0[cent],"TCM fit","l");
        }

        DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombPi0Sys, 20, markerSizeCent[0]*0.75, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphPPInvYieldCombPi0Sys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombPi0StatW0XErr, 20, markerSizeCent[0]*0.75, kBlack, kBlack);
        graphPPInvYieldCombPi0StatW0XErr->Draw("p,same,z");
        DrawGammaSetMarkerTF1( fitPPTCMInvYieldPi0, 7, 2, kGray+2);
        fitPPTCMInvYieldPi0->Draw("same");
        legendXsectionPaper->AddEntry(graphPPInvYieldCombPi0Sys,collisionSystempp,"pf");
        legendXsectionPaper->AddEntry(fitPPTCMInvYieldPi0,"TCM fits","l");

        TLatex *labelEnergyXSectionPaper= new TLatex(0.935, 0.91, collisionSystemPbPb.Data());
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

        for (Int_t cent= 0; cent < 1; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
            graphRatioPi0CombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioPi0CombCombFitSys[cent]->SetLineWidth(0);
            graphRatioPi0CombCombFitSys[cent]->Draw("2,same");
            graphRatioPi0CombCombFitStatWOXErr[cent]->Draw("p,same");
        }
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting,1., 1.,1,kGray+1);

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

        for (Int_t cent= 1; cent < 3; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
            graphRatioPi0CombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioPi0CombCombFitSys[cent]->SetLineWidth(0);
            graphRatioPi0CombCombFitSys[cent]->Draw("2,same");
            graphRatioPi0CombCombFitStatWOXErr[cent]->Draw("p,same");
        }
        boxErrorNormRatioPi0->Draw();
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting,1., 1.,1,kGray+1);

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

        for (Int_t cent= 3; cent < 5; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
            graphRatioPi0CombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioPi0CombCombFitSys[cent]->SetLineWidth(0);
            graphRatioPi0CombCombFitSys[cent]->Draw("2,same");
            graphRatioPi0CombCombFitStatWOXErr[cent]->Draw("p,same");
        }

        boxErrorNormRatioPi0->Draw();
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting,1., 1.,1,kGray+1);

    canvasInvSectionPaper->Print(Form("%s/Pi0_InvYieldWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));

    //*************************************************************************************************************
    //***************************** Paper plot invyield eta with soft-physics *************************************
    //*************************************************************************************************************
    TLegend* legendXsectionPaperEta    = GetAndSetLegend2(0.19, 0.03, 0.7, 0.03+0.05*(maxCent+1), textSizeLabelsPixel, 1, "", 43, 0.15);
//     TLegend* legendXsectionPaperEta    = GetAndSetLegend2(0.19, 0.03, 0.7, 0.03+0.05*(maxCent), textSizeLabelsPixel, 2, "", 43, 0.25);

    canvasInvSectionPaper->cd();
    padInvSectionSpec->cd();
    padInvSectionSpec->SetLogy(1);
    padInvSectionSpec->SetLogx(1);
        SetStyleHistoTH2ForGraphs(histo2DYieldEta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                                0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.257/(textsizeFacXSecUp*marginXSec));
        histo2DYieldEta->GetXaxis()->SetMoreLogLabels();
        histo2DYieldEta->GetXaxis()->SetLabelOffset(+0.01);
        histo2DYieldEta->Draw();

        for (Int_t cent= 0; cent < maxCent; cent++ ){
            if (!enableCentComb[cent] ) continue;
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
            graphCombEtaInvYieldSys[cent]->Draw("E2same");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent]);
            graphCombEtaInvYieldStatWOXErr[cent]->Draw("p,same,z");
            DrawGammaSetMarkerTF1( fitTCMInvYieldEta[cent], 7, 2, colorCentMC[cent]);
            fitTCMInvYieldEta[cent]->Draw("same");
            legendXsectionPaperEta->AddEntry(graphCombEtaInvYieldSys[cent],centArray2[cent].Data(),"pf");
//             legendXsectionPaperEta->AddEntry(fitTCMInvYieldEta[cent],"TCM fit","l");
        }
        DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombEtaSys, 20, markerSizeCent[0]*0.75, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphPPInvYieldCombEtaSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphPPInvYieldCombEtaStatW0XErr, 20, markerSizeCent[0]*0.75, kBlack, kBlack);
        graphPPInvYieldCombEtaStatW0XErr->Draw("p,same,z");
        DrawGammaSetMarkerTF1( fitPPTCMInvYieldEta, 7, 2, kGray+2);
        fitPPTCMInvYieldEta->Draw("same");
        legendXsectionPaperEta->AddEntry(graphPPInvYieldCombEtaSys,collisionSystempp,"pf");
        legendXsectionPaperEta->AddEntry(fitPPTCMInvYieldEta,"TCM fits","l");

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
        for (Int_t cent= 0; cent < 1; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
            graphRatioEtaCombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioEtaCombCombFitSys[cent]->SetLineWidth(0);
            graphRatioEtaCombCombFitSys[cent]->Draw("2,same");
            graphRatioEtaCombCombFitStatWOXErr[cent]->Draw("p,same");
        }

        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting,1., 1.,1,kGray);

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

        for (Int_t cent= 1; cent < 3; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
            graphRatioEtaCombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioEtaCombCombFitSys[cent]->SetLineWidth(0);
            graphRatioEtaCombCombFitSys[cent]->Draw("2,same");
            graphRatioEtaCombCombFitStatWOXErr[cent]->Draw("p,same");
        }
        boxErrorNormRatioEta->Draw();

        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting,1., 1.,1,kGray);

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

        for (Int_t cent= 3; cent < 5; cent++ ){
            if (!enableCentComb[cent]) continue;
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kFALSE);
            graphRatioEtaCombCombFitStatWOXErr[cent]->SetLineWidth(widthLinesBoxes);
            DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioEtaCombCombFitSys[cent]->SetLineWidth(0);
            graphRatioEtaCombCombFitSys[cent]->Draw("2,same");
            graphRatioEtaCombCombFitStatWOXErr[cent]->Draw("p,same");
        }
        boxErrorNormRatioEta->Draw();
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting,1., 1.,1,kGray);

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

        labelEnergyRatioPaper[cent]     = new TLatex(0.935, yPosLabel[cent],Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
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
            DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting,1., 1.,1,kGray+1);
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
        DrawGammaLines(minPtEtaPlotting, maxPtEtaPlotting,1., 1.,1,kGray+1);
        if (labelEnergyRatioPaper[cent]) labelEnergyRatioPaper[cent]->Draw();
        if (cent == 0) labelDetSysRatioPaperEta->Draw();
        ratio2DDataAndTheo->DrawCopy("axis,same");
        delete ratio2DDataAndTheo;
    }
    canvasRatioGeneratorsPaper->Print(Form("%s/Eta_RatioToGenerator_Paper.%s",outputDir.Data(),suffix.Data()));

    // ***************************************************************************************************************
    // ******************************** fitting eta/pi0 **************************************************************
    // ***************************************************************************************************************

    TF1* etaToPi0ConstData[9];
    TF1* etaToPi0ConstDataStat[9];
    TF1* etaToPi0ConstMC[9];
    TF1* etaToPi0ConstMC2[9];
    TF1* etaToPi0ConstMC3[9];
    for (Int_t cent = 0; cent < maxCent; cent++){
        if (!enableCentComb[cent]) continue;
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
        cout << centArray[cent].Data() << "\t" << addCentString[cent].Data()  << endl;
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
        fileFitsOutput << centArray[cent].Data() << "\t" << addCentString[cent].Data() << endl;
        fileFitsOutput << "high pt eta/pi0 - data, stat: " << etaToPi0ConstDataStat[cent]->GetParameter(0) << "+-"<< etaToPi0ConstDataStat[cent]->GetParError(0) << endl;
        fileFitsOutput << "high pt eta/pi0 - data, tot: " << etaToPi0ConstData[cent]->GetParameter(0) << "+-"<< etaToPi0ConstData[cent]->GetParError(0) << endl;
        if(histoDPMJetEtaToPi0[cent]) fileFitsOutput << "high pt eta/pi0 - DPMJet: " << etaToPi0ConstMC[cent]->GetParameter(0) << "+-"<< etaToPi0ConstMC[cent]->GetParError(0) << endl;
        if(histoHIJINGEtaToPi0[cent]) fileFitsOutput << "high pt eta/pi0 - HIJING: " << etaToPi0ConstMC2[cent]->GetParameter(0) << "+-"<< etaToPi0ConstMC2[cent]->GetParError(0) << endl;
        if(histoEPOSLHCEtaToPi0[cent]) fileFitsOutput << "high pt eta/pi0 - EPOS-LHC: " << etaToPi0ConstMC3[cent]->GetParameter(0) << "+-"<< etaToPi0ConstMC3[cent]->GetParError(0) << endl;
        fileFitsOutput << "***********************************************************************************************************" << endl;
        fileFitsOutput << "***********************************************************************************************************" << endl;
        fileFitsOutput << "***********************************************************************************************************" << endl;

    }

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

    TLatex *labelEnergyEtaToPi0 = new TLatex(0.13, 0.92,collisionSystemPbPb.Data());
    SetStyleTLatex( labelEnergyEtaToPi0, textsizeLabelsEtaToPi0,4, 1, 42, kTRUE, 11);
    labelEnergyEtaToPi0->Draw();
    TLatex *labelALICEEtaToPi0 = new TLatex(0.13, 0.92-(1*textsizeLabelsEtaToPi0),textALICE.Data());
    SetStyleTLatex( labelALICEEtaToPi0, textsizeLabelsEtaToPi0,4, 1, 42, kTRUE, 11);
    labelALICEEtaToPi0->Draw();

    for (Int_t cent = 0; cent < maxCent; cent++){
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

        labelEnergyEtaToPi0->SetText(0.13, 0.92,Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
        labelEnergyEtaToPi0->Draw();
        labelALICEEtaToPi0->Draw();

        histo2DEtatoPi0combo->Draw("same,axis");

        canvasEtatoPi0combo->Update();
        canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_differentSystems_%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));


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
        canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_differentSystemsWComb_%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

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
        canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Theory_Paper_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

        histo2DEtatoPi0combo->Draw("copy");

            boxHighPtFit->Draw();
            eta2pi0MtScaledTCM->Draw("same,hist,l");
            etaToPi0ConstDataStat[cent]->Draw("same");

            graphCombEtaToPi0Sys[cent]->Draw("E2same");
            graphCombEtaToPi0StatWOXErr[cent]->Draw("p,same,z");

            // plotting labels
            labelEnergyEtaToPi0->Draw();
            labelALICEEtaToPi0->Draw();

            legendEtaToPi0TheoryHighpt->Draw();

            if (graphPPCombEtaToPi0StatW0XErr && graphPPCombEtaToPi0Sys){
                DrawGammaSetMarkerTGraphAsym(graphPPCombEtaToPi0Sys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
                graphPPCombEtaToPi0Sys->Draw("E2same");
                DrawGammaSetMarkerTGraphAsym(graphPPCombEtaToPi0StatW0XErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
                graphPPCombEtaToPi0StatW0XErr->Draw("p,same,z");
            }

            TLegend* legendEtaToPi0PP = GetAndSetLegend2(0.58, 0.15+(textsizeLabelsEtaToPi0*1*0.9), 0.88, 0.15, textSizeLabelsPixel*0.85, 1, "", 43, 0.2);
            legendEtaToPi0PP->AddEntry(graphPPCombEtaToPi0Sys,collisionSystempp.Data(),"pf");
            legendEtaToPi0PP->Draw();

        histo2DEtatoPi0combo->Draw("axis,same");

        canvasEtatoPi0combo->Update();
        canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_WithPP_Paper_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));



//         if (graphChKToPiTot[cent]){
//             histo2DParticleRatios->Draw("copy");
//
//                 graphCombEtaToPi0Sys[cent]->Draw("E2same");
//                 graphCombEtaToPi0StatWOXErr[cent]->Draw("p,same,z");
//
// //                 DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
// //                 graphCombEtaToPi0Sys[cent]->Draw("E2same");
//                 DrawGammaSetMarkerTGraphAsym(graphChKToPiTot[cent], 25, markerSizeCent[cent]*0.5, kGray+1, kGray+1);
//                 graphChKToPiTot[cent]->Draw("p,same,z");
//
//                 // plotting labels
//                 labelEnergyEtaToPi0->Draw();
//                 labelALICEEtaToPi0->Draw();
//
//                 TLegend* legendParticleRatios = GetAndSetLegend2(0.13, 0.905-(3*textsizeLabelsEtaToPi0*1.05), 0.53, 0.905-(1*textsizeLabelsEtaToPi0), textSizeLabelsPixel*0.85, 1, "", 43, 0.14);
//                 legendParticleRatios->AddEntry(graphCombEtaToPi0Sys[cent], "#eta/#pi^{0}","fp");
//                 legendParticleRatios->AddEntry(graphChKToPiTot[cent], "K^{#pm}/#pi^{#pm}","ep");
//                 legendParticleRatios->Draw();
//
//             canvasEtatoPi0combo->Update();
//             canvasEtatoPi0combo->SaveAs(Form("%s/ParticleRatios_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
//             delete legendParticleRatios;
//
//         }

        delete legendEtaToPi0TheoryHighpt;
        delete dummyHighPt;
        delete boxHighPtFit;
        delete eta2pi0MtScaledTCM;
        delete legendEtaToPi0;
        delete legendEtaToPi0Theory;
    }
    Bool_t haveCharged  = kFALSE;

    histo2DEtatoPi0combo->Draw("copy");
    TLegend* legendEtaToPi0Comb = GetAndSetLegend2(0.13, 0.91-(1*textsizeLabelsEtaToPi0*1.)-0.01, 0.6, 0.91-(nCurrEst*textsizeLabelsEtaToPi0*1.01)-0.01, textSizeLabelsPixel,1, "", 43, 0.125);

        for (Int_t cent = 0; cent < maxCent; cent++){
            if (!enableCentComb[cent] ) continue;
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
        labelEnergyEtaToPi0->SetText(0.13, 0.92,Form("%s",  collisionSystemPbPb.Data()));
        labelEnergyEtaToPi0->Draw();
        labelALICEEtaToPi0->Draw();
        legendEtaToPi0Comb->Draw();
        histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Paper_%s.%s",outputDir.Data(), nameCentEst.Data(), suffix.Data()));
    histo2DEtatoPi0combo->Draw("copy");

        for (Int_t cent = 0; cent < maxCent; cent++){
            if (!enableCentComb[cent] ) continue;
            graphCombEtaToPi0Sys[cent]->Draw("E2same");
            graphCombEtaToPi0StatWOXErr[cent]->Draw("p,same,z");
        }

        if (enableCentComb[0] ){
            graphCombEtaToPi0Sys[0]->Draw("E2same");
            graphCombEtaToPi0StatWOXErr[0]->Draw("p,same,z");
        }

        // plotting labels
        labelEnergyEtaToPi0->Draw();
        labelALICEEtaToPi0->Draw();
        legendEtaToPi0Comb->Draw();

        if (graphPPCombEtaToPi0StatW0XErr && graphPPCombEtaToPi0Sys){
            graphPPCombEtaToPi0Sys->Draw("E2same");
            graphPPCombEtaToPi0StatW0XErr->Draw("p,same,z");
        }

        TLegend* legendEtaToPi0PP2 = GetAndSetLegend2(0.58, 0.15+(textsizeLabelsEtaToPi0*1*0.9), 0.88, 0.15, textSizeLabelsPixel*0.85, 1, "", 43, 0.2);
        legendEtaToPi0PP2->AddEntry(graphPPCombEtaToPi0Sys,collisionSystempp.Data(),"pf");
        legendEtaToPi0PP2->Draw();

        histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_WithPP_Paper_%s.%s",outputDir.Data(), nameCentEst.Data(), suffix.Data()));

//     if (haveCharged){
//
//         histo2DParticleRatios->Draw("copy");
//
//         TLegend* legendParticleRatios = GetAndSetLegend2(0.13, 0.905-(3*textsizeLabelsEtaToPi0*1.05), 0.63, 0.905-(1*textsizeLabelsEtaToPi0), textSizeLabelsPixel*0.85, 2, "", 43, 0.2);
//         for (Int_t cent = 0; cent < maxCent; cent++){
//             if (!(centArray[cent].CompareTo("0-20%") == 0 ||centArray[cent].CompareTo("60-100%") == 0 )) continue;
//             if (!enableCentComb[cent] ) continue;
//             legendParticleRatios->AddEntry(graphCombEtaToPi0Sys[cent], "#eta/#pi^{0} "+centArray[cent],"fp");
//
//             graphCombEtaToPi0Sys[cent]->Draw("E2same");
//             graphCombEtaToPi0StatWOXErr[cent]->Draw("p,same,z");
//         }
//         for (Int_t cent = 0; cent < maxCent; cent++){
//             if (! (centArray[cent].CompareTo("0-20%") == 0 ||centArray[cent].CompareTo("60-100%") == 0 )) continue;
//             if (!enableCentComb[cent] ) continue;
//             DrawGammaSetMarkerTGraphAsym(graphChKToPiTot[cent], 25, markerSizeCent[cent]*0.5, colorCentBox[cent], colorCentBox[cent]);
//             graphChKToPiTot[cent]->Draw("p,same,z");
//             legendParticleRatios->AddEntry(graphChKToPiTot[cent], "K^{#pm}/#pi^{#pm} "+centArray[cent],"ep");
//
//         //                 DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys[cent], markerStyleCent[cent], markerSizeCent[cent], colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
//         //                 graphCombEtaToPi0Sys[cent]->Draw("E2same");
//         }
//         legendParticleRatios->Draw();
//         // plotting labels
//         labelEnergyEtaToPi0->Draw();
//         labelALICEEtaToPi0->Draw();
//
//
//         canvasEtatoPi0combo->Update();
//         canvasEtatoPi0combo->SaveAs(Form("%s/ParticleRatios_CentAndPeri_%s.%s",outputDir.Data(), nameCentEst.Data(), suffix.Data()));
//         delete legendParticleRatios;
//     }

    //*****************************************************************************************************************
    // Plotting RAA s for different mesons
    //*****************************************************************************************************************
    textSizeLabelsPixel                         = 54;
    TCanvas* canvasRAA = new TCanvas("canvasRAA","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasRAA,  0.085, 0.01, 0.015, 0.08);
    canvasRAA->SetLogx();
    TH2F * histo2DRAA  = new TH2F("histo2DRAA","histo2DRAA",1000,minPtPi0Plotting,maxPtPi0Plotting,1000,0.0,1.41);
    SetStyleHistoTH2ForGraphs(histo2DRAA, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.035,0.04, 0.035,0.04, 0.8,1., 512, 505);
    histo2DRAA->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRAA->DrawCopy();
    TH2F * histo2DRAAEta  = new TH2F("histo2DRAAEta","histo2DRAAEta",1000,minPtEtaPlotting,maxPtEtaPlotting,1000,0.0,1.41);
    SetStyleHistoTH2ForGraphs(histo2DRAAEta, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.035,0.04, 0.035,0.04, 0.8,1., 512, 505);
    histo2DRAAEta->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRAAEta->DrawCopy();

    TLatex *labelALICERAA  = new TLatex(0.12,0.95-0.04*2.1,textALICE.Data());
    SetStyleTLatex( labelALICERAA, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
    TLatex *labelEnergyRAA     = new TLatex(0.12, 0.95-0.04*1, collisionSystemPbPb.Data());
    SetStyleTLatex( labelEnergyRAA, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
    TLatex *labelPi0RAA  = new TLatex(0.12,0.95-0.04*3,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelPi0RAA, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
    TLatex *labelEtaRAA  = new TLatex(0.12,0.95-0.04*3,"#eta #rightarrow #gamma#gamma");
    SetStyleTLatex( labelEtaRAA, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);

    canvasRAA->cd();
    histo2DRAA->Draw("copy");
    TLegend* legendPi0RAAPaper    = GetAndSetLegend2(0.95-0.4, 0.95-0.04*1.25*(nCurrEst+1)/2, 0.95, 0.95, textSizeLabelsPixel, 2, "", 43, 0.25);
        Int_t plot= 0;
        for (Int_t cent= 0; cent < maxCent; cent++ ){
            if (!enableCentComb[cent] || !enableCentRAA[cent] ) continue;
            DrawGammaSetMarkerTGraphAsym(graphRAACombSystPi0[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
            graphRAACombSystPi0[cent]->Draw("E2same");
            legendPi0RAAPaper->AddEntry(graphRAACombSystPi0[cent],centArray[cent].Data(),"pf");
            TBox* boxErrorNormRAAPi0          = CreateBoxConv(colorCentBox[cent], 0.24+plot*0.01, 1.-(nCollErrPbPb[cent]/nCollPbPb[cent]), 0.25+plot*0.01, 1.+(nCollErrPbPb[cent]/nCollPbPb[cent]));
            boxErrorNormRAAPi0->Draw();
            plot++;
        }
        DrawGammaLines(minPtPi0Plotting,maxPtPi0Plotting , 1, 1 ,1, kGray, 7);
        for (Int_t cent= 0; cent < maxCent; cent++ ){
            if (!enableCentComb[cent] || !enableCentRAA[cent] ) continue;
            DrawGammaSetMarkerTGraphAsym(graphRAACombStatPi0WOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent]);
            graphRAACombStatPi0WOXErr[cent]->Draw("p,same,z");
        }

        labelALICERAA->Draw();
        labelEnergyRAA->SetText(0.12, 0.95-0.04*1,Form("%s", collisionSystemPbPb.Data()));
        labelEnergyRAA->Draw();
        labelPi0RAA->Draw();
        legendPi0RAAPaper->Draw();
    canvasRAA->SaveAs(Form("%s/Pi0_RAA_AllCent_%s.%s",outputDir.Data(), nameCentEst.Data(), suffix.Data()));


    canvasRAA->cd();
    histo2DRAAEta->Draw("copy");
    TLegend* legendEtaRAAPaper    = GetAndSetLegend2(0.95-0.4, 0.95-0.04*1.25*(nCurrEst+1)/2, 0.95, 0.95, textSizeLabelsPixel, 2, "", 43, 0.25);

        plot= 0;
        for (Int_t cent= 0; cent < maxCent; cent++ ){
            if (!enableCentComb[cent] || !enableCentRAA[cent] ) continue;
            DrawGammaSetMarkerTGraphAsym(graphRAACombSystEta[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
            graphRAACombSystEta[cent]->Draw("E2same");
            legendEtaRAAPaper->AddEntry(graphRAACombSystEta[cent],centArray[cent].Data(),"pf");
            TBox* boxErrorNormRAAEta          = CreateBoxConv(colorCentBox[cent], 0.45+plot*0.01, 1.-(nCollErrPbPb[cent]/nCollPbPb[cent]), 0.46+plot*0.01, 1.+(nCollErrPbPb[cent]/nCollPbPb[cent]));
            boxErrorNormRAAEta->Draw();
            plot++;
        }
        DrawGammaLines(minPtEtaPlotting,maxPtEtaPlotting , 1, 1 ,1, kGray, 7);

        for (Int_t cent= 0; cent < maxCent; cent++ ){
            if (!enableCentComb[cent] || !enableCentRAA[cent] ) continue;
            DrawGammaSetMarkerTGraphAsym(graphRAACombStatEtaWOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent]);
            graphRAACombStatEtaWOXErr[cent]->Draw("p,same,z");
        }

        labelALICERAA->Draw();
        labelEnergyRAA->Draw();
        labelEtaRAA->Draw();
        legendEtaRAAPaper->Draw();
    canvasRAA->SaveAs(Form("%s/Eta_RAA_AllCent_%s.%s",outputDir.Data(), nameCentEst.Data(), suffix.Data()));


    if (bWCorrection.Contains("Y")){
        for (Int_t cent = 0; cent < maxCent; cent++){
            if (!enableCentRAA[cent]) continue;

            TBox* boxErrorNormRAAPi0          = CreateBoxConv(colorCentBox[cent], 0.24, 1.-(nCollErrPbPb[cent]/nCollPbPb[cent]), 0.25, 1.+(nCollErrPbPb[cent]/nCollPbPb[cent]));
            TBox* boxErrorNormRAAEta          = CreateBoxConv(colorCentBox[cent], 0.45, 1.-(nCollErrPbPb[cent]/nCollPbPb[cent]), 0.46, 1.+(nCollErrPbPb[cent]/nCollPbPb[cent]));

            histo2DRAA->GetYaxis()->SetTitle("#it{R}_{AA}");
            histo2DRAAEta->GetYaxis()->SetTitle("#it{R}_{AA}");
            canvasRAA->cd();
            histo2DRAA->DrawCopy();
            boxErrorNormRAAPi0->Draw();

            if (graphRAACombSystPi0[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRAACombSystPi0[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRAACombSystPi0[cent]->Draw("E2same");
            }
            if (graphRAACombSystEta[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRAACombSystEta[cent], markerStyleCentMC[cent], markerSizeCentMC[cent]*0.5, colorCentMC[cent], colorCentMC[cent], widthLinesBoxes, kTRUE);
                graphRAACombSystEta[cent]->Draw("E2same");
            }
            DrawGammaLines(minPtPi0Plotting,maxPtPi0Plotting , 1, 1 ,1, kGray, 7);
            if (graphRAACombStatPi0WOXErr[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRAACombStatPi0WOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRAACombStatPi0WOXErr[cent]->Draw("p,same,z");
            }

            if (graphRAACombStatEtaWOXErr[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRAACombStatEtaWOXErr[cent], markerStyleCentMC[cent], markerSizeCentMC[cent]*0.5, colorCentMC[cent], colorCentMC[cent], widthLinesBoxes, kTRUE);
                graphRAACombStatEtaWOXErr[cent]->Draw("p,same,z");
            }

            labelEnergyRAA->SetText(0.12, 0.95-0.04*1,Form("%s %s", centArray3[cent].Data(), collisionSystemPbPb.Data()));
            labelEnergyRAA->Draw();
            labelALICERAA->Draw();

            TLegend* legendRAAComb     = GetAndSetLegend2(0.12, 0.95-0.04*1.25*4, 0.32 , 0.95-0.04*1.25*2, textSizeLabelsPixel*0.85,1,"",43,0.3);
            legendRAAComb->AddEntry(graphRAACombSystPi0[cent],"#pi^{0}","pf");
            legendRAAComb->AddEntry(graphRAACombSystEta[cent],"#eta","pf");
            legendRAAComb->Draw();

            canvasRAA->Update();
            canvasRAA->Print(Form("%s/Pi0AndEta_RAA_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
            delete legendRAAComb;

            boxErrorNormRAAPi0->Draw();
//             if (graphDRAASys[cent]){
//                 DrawGammaSetMarkerTGraphAsym(graphDRAASys[cent], 24, markerSizeCentMC[cent]*0.5, kGray+2, kGray+2, widthLinesBoxes, kTRUE);
//                 graphDRAASys[cent]->Draw("E2same");
//             }
//             if (graphChHadRAASys[cent]){
//                 DrawGammaSetMarkerTGraphAsym(graphChHadRAASys[cent], 25, markerSizeCentMC[cent]*0.5, kGray+1, kGray+1, widthLinesBoxes, kTRUE);
//                 graphChHadRAASys[cent]->Draw("E2same");
//             }

            if (graphRAACombSystPi0[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRAACombSystPi0[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRAACombSystPi0[cent]->Draw("E2same");
            }
            if (graphRAACombSystEta[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRAACombSystEta[cent], markerStyleCentMC[cent], markerSizeCentMC[cent]*0.5, colorCentMC[cent], colorCentMC[cent], widthLinesBoxes, kTRUE);
                graphRAACombSystEta[cent]->Draw("E2same");
            }
            DrawGammaLines(minPtPi0Plotting,maxPtPi0Plotting , 1, 1 ,1, kGray, 7);
//             if (graphChHadRAATot[cent]){
//                 gStyle->SetEndErrorSize(2);
//                 DrawGammaSetMarkerTGraphAsym(graphChHadRAATot[cent], 21, markerSizeCentMC[cent]*0.25, kGray+1, kGray+1, widthLinesBoxes, kTRUE);
//                 graphChHadRAATot[cent]->Draw("p,same");
//                 gStyle->SetEndErrorSize(0);
//             }
            if (graphRAACombStatPi0WOXErr[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRAACombStatPi0WOXErr[cent], markerStyleCent[cent], markerSizeCent[cent]*0.5, colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE);
                graphRAACombStatPi0WOXErr[cent]->Draw("p,same,z");
            }
//             if (graphDRAATot[cent]){
//                 gStyle->SetEndErrorSize(2);
//                 DrawGammaSetMarkerTGraphAsym(graphDRAATot[cent], 20, markerSizeCentMC[cent]*0.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
//                 graphDRAATot[cent]->Draw("p,same");
//                 gStyle->SetEndErrorSize(0);
//             }
            if (graphRAACombStatEtaWOXErr[cent]){
                DrawGammaSetMarkerTGraphAsym(graphRAACombStatEtaWOXErr[cent], markerStyleCentMC[cent], markerSizeCentMC[cent]*0.5, colorCentMC[cent], colorCentMC[cent], widthLinesBoxes, kTRUE);
                graphRAACombStatEtaWOXErr[cent]->Draw("p,same,z");
            }



            labelEnergyRAA->SetText(0.12, 0.95-0.04*1,Form("%s %s", centArray3[cent].Data(), collisionSystemPbPb.Data()));
            labelEnergyRAA->Draw();
            labelALICERAA->Draw();

            TLegend* legendRAAComb2     = GetAndSetLegend2(0.12, 0.95-0.04*1.25*4, 0.42 , 0.95-0.04*1.25*2, textSizeLabelsPixel*0.85,2,"",43,0.3);
            legendRAAComb2->AddEntry(graphRAACombSystPi0[cent],"#pi^{0}","pf");
            legendRAAComb2->AddEntry(graphRAACombSystEta[cent],"#eta","pf");
//             if (graphDRAATot[cent])legendRAAComb2->AddEntry(graphDRAATot[cent],"D^{0}","pe");
//             if (graphChHadRAATot[cent])legendRAAComb2->AddEntry(graphChHadRAATot[cent],"h^{#pm}","pe");
            legendRAAComb2->Draw();

            canvasRAA->Update();
            canvasRAA->Print(Form("%s/Pi0AndEtaPlusOther_RAA_%s%s.%s",outputDir.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
            delete legendRAAComb2;

            //*****************************************************************************************************************
            // Plotting RAA pi0 diff rec techniques
            //*****************************************************************************************************************
            histo2DRAA->DrawCopy();

            boxErrorNormRAAPi0->Draw();
            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRAAIndSystPi0[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRAAIndSystPi0[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth], widthLinesBoxes, kTRUE);
                    graphRAAIndSystPi0[cent][meth]->Draw("E2same");
                }
            }

            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRAAIndStatPi0WOXErr[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRAAIndStatPi0WOXErr[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth]);
                    graphRAAIndStatPi0WOXErr[cent][meth]->Draw("p,same,z");
                }
            }

            DrawGammaLines(minPtPi0Plotting,maxPtPi0Plotting , 1, 1 ,1, kGray, 7);

            labelEnergyRAA->Draw();
            labelALICERAA->Draw();
            labelPi0RAA->Draw();

            TLegend* legendRAAInd     = GetAndSetLegend2(0.45, 0.13, 0.95 , 0.13+0.04*3, textSizeLabelsPixel*0.85,2, "", 43, 0.2);
            for (Int_t meth = 0; meth< 11; meth++){
                if (graphRAAIndSystPi0[cent][meth])legendRAAInd->AddEntry(graphRAAIndSystPi0[cent][meth],nameMeasGlobalLabel[meth].Data(),"pf");
            }
            legendRAAInd->Draw();

            canvasRAA->Update();
            canvasRAA->Print(Form("%s/Pi0_RAA_IndividualMeasurements_%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
            delete legendRAAInd;

            //*****************************************************************************************************************
            // Plotting RAA eta diff rec techniques
            //*****************************************************************************************************************
            histo2DRAAEta->DrawCopy();
            boxErrorNormRAAEta->Draw();
            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRAAIndSystEta[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRAAIndSystEta[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth], widthLinesBoxes, kTRUE);
                    graphRAAIndSystEta[cent][meth]->Draw("E2same");
                }
            }
            for (Int_t meth = 10; meth> -1; meth--){
                if (graphRAAIndStatEtaWOXErr[cent][meth]){
                    DrawGammaSetMarkerTGraphAsym(graphRAAIndStatEtaWOXErr[cent][meth], markerStyleDet[meth], markerSizeDet[meth]*0.5, colorDet[meth] , colorDet[meth]);
                    graphRAAIndStatEtaWOXErr[cent][meth]->Draw("p,same,z");
                }
            }

            DrawGammaLines(minPtEtaPlotting,maxPtEtaPlotting , 1, 1 ,1, kGray, 7);

            labelEnergyRAA->Draw();
            labelALICERAA->Draw();
            labelEtaRAA->Draw();

            TLegend* legendRAAIndEta     = GetAndSetLegend2(0.45, 0.13, 0.95 , 0.13+0.04*3, textSizeLabelsPixel*0.85,2, "", 43, 0.2);
            for (Int_t meth = 0; meth< 11; meth++){
                if (graphRAAIndSystEta[cent][meth])legendRAAIndEta->AddEntry(graphRAAIndSystEta[cent][meth],nameMeasGlobalLabel[meth].Data(),"pf");
            }
            legendRAAIndEta->Draw();

            canvasRAA->Update();
            canvasRAA->Print(Form("%s/Eta_RAA_IndividualMeasurements_%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));
            delete legendRAAIndEta;
        }
    }



    cout << "starting charged pion Comparison" << endl;
    TGraphErrors* graphChPiStatRebinnedPi0Comb[9]       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphErrors* graphChPiSystRebinnedPi0Comb[9]       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphErrors* graphCombPi0StatRebinnedChPi[9]       = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphErrors* graphCombPi0SysRebinnedChPi[9]        = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphErrors* graphRatioCombPi0ChPiComb[9]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphErrors* graphRatioIndPi0ChPiComb[9][11];
    TGraphErrors* graphChPiStatRebinnedPi0Ind[9][11];
    TGraphErrors* graphChPiSystRebinnedPi0Ind[9][11];
    TGraphErrors* graphIndPi0StatRebinnedChPi[9][11];
    TGraphErrors* graphIndPi0SysRebinnedChPi[9][11];
    for (Int_t cent = 0; cent < maxCent; cent++){
        for (Int_t meth = 0; meth < 11; meth++){
            graphRatioIndPi0ChPiComb[cent][meth]        = NULL;
            graphChPiStatRebinnedPi0Ind[cent][meth]     = NULL;
            graphChPiSystRebinnedPi0Ind[cent][meth]     = NULL;
            graphIndPi0StatRebinnedChPi[cent][meth]     = NULL;
            graphIndPi0SysRebinnedChPi[cent][meth]      = NULL;
        }
    }

    for (Int_t cent = 0; cent < maxCent; cent++){
        if (!enableCent[cent]) continue;
        for (Int_t meth = 0; meth < 11; meth++){
            if (graphIndPi0InvYieldStat[cent][meth] && graphIndPi0InvYieldSys[cent][meth])
                graphRatioIndPi0ChPiComb[cent][meth]     = CalculateRatioBetweenSpectraWithDifferentBinning( graphIndPi0InvYieldStat[cent][meth], graphIndPi0InvYieldSys[cent][meth],
                                                                                                        histChPiSpecStat[cent], histChPiSpecSyst[cent],
                                                                                                        kTRUE,  kTRUE,
                                                                                                        &graphIndPi0StatRebinnedChPi[cent][meth], &graphIndPi0SysRebinnedChPi[cent][meth],
                                                                                                        &graphChPiStatRebinnedPi0Ind[cent][meth], &graphChPiSystRebinnedPi0Ind[cent][meth] )    ;
        }

        if (!enableCentComb[cent]) continue;
        cout << "charged pion comp for: " << centArray[cent].Data() << endl;
        cout << graphCombPi0InvYieldStat[cent] << "\t" << graphCombPi0InvYieldSys[cent] << "\t" << histChPiSpecStat[cent] << "\t" << histChPiSpecSyst[cent] << endl;
        graphRatioCombPi0ChPiComb[cent]     = CalculateRatioBetweenSpectraWithDifferentBinning( graphCombPi0InvYieldStat[cent], graphCombPi0InvYieldSys[cent],
                                                                                                histChPiSpecStat[cent], histChPiSpecSyst[cent],
                                                                                                kTRUE,  kTRUE,
                                                                                                &graphCombPi0StatRebinnedChPi[cent], &graphChPiSystRebinnedPi0Comb[cent],
                                                                                                &graphChPiStatRebinnedPi0Comb[cent], &graphChPiSystRebinnedPi0Comb[cent] )    ;
    }


    TCanvas* canvasRatioChargedPion      = new TCanvas("canvasRatioChargedPion","",0,0,1250,2000);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioChargedPion,  0.13, 0.02, 0.03, 0.06);
    canvasRatioChargedPion->cd();

    TPad* padRatioChargedPion[5]             = {NULL, NULL, NULL, NULL, NULL};
    Double_t textSizeLablesChPionRatio[5]    = {0};
    Double_t textsizeFacXSecChPionRatio[5]   = {0};
    marginXSec                              = relativeMarginsXXSec2[0]*1250;
    for (Int_t i = 0; i < 5; i++){
        padRatioChargedPion[i]               = new TPad(Form("padRatioGen%i",i), "", arrayBoundariesX1_XSec2[0], arrayBoundariesY1_XSec2[i+1], arrayBoundariesX1_XSec2[1], arrayBoundariesY1_XSec2[i],-1, -1, -2);
        if (i == 0) DrawGammaPadSettings( padRatioChargedPion[i], relativeMarginsXXSec2[0], relativeMarginsXXSec2[2], relativeMarginsYXSec2[0], relativeMarginsYXSec2[1]);
        else if (i == 4) DrawGammaPadSettings( padRatioChargedPion[i], relativeMarginsXXSec2[0], relativeMarginsXXSec2[2], relativeMarginsYXSec2[1], relativeMarginsYXSec2[2]);
        else DrawGammaPadSettings( padRatioChargedPion[i], relativeMarginsXXSec2[0], relativeMarginsXXSec2[2], relativeMarginsYXSec2[1], relativeMarginsYXSec2[1]);
        padRatioChargedPion[i]->Draw();

        if (padRatioChargedPion[i]->XtoPixel(padRatioChargedPion[i]->GetX2()) < padRatioChargedPion[i]->YtoPixel(padRatioChargedPion[i]->GetY1())){
            textSizeLablesChPionRatio[i]            = (Double_t)textSizeLabelsPixel/padRatioChargedPion[i]->XtoPixel(padRatioChargedPion[i]->GetX2()) ;
            textsizeFacXSecChPionRatio[i]           = (Double_t)1./padRatioChargedPion[i]->XtoPixel(padRatioChargedPion[i]->GetX2()) ;
        } else {
            textSizeLablesChPionRatio[i]            = (Double_t)textSizeLabelsPixel/padRatioChargedPion[i]->YtoPixel(padRatioChargedPion[i]->GetY1());
            textsizeFacXSecChPionRatio[i]           = (Double_t)1./padRatioChargedPion[i]->YtoPixel(padRatioChargedPion[i]->GetY1());
        }
    }


    for (Int_t cent = 0; cent < 5; cent++ ){
        padRatioChargedPion[cent]->cd();
        padRatioChargedPion[cent]->SetLogx(1);

        TH2F* ratio2DNPiToChPi               = new TH2F(Form ("ratio2DNPiToChPi%i",cent),Form ("ratio2DNPiToChPi%i",cent),1000,minPtPi0Plotting, maxPtPi0Plotting,1000,0.653,1.95);
        SetStyleHistoTH2ForGraphs(ratio2DNPiToChPi, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.85*textSizeLablesChPionRatio[cent], textSizeLablesChPionRatio[cent],
                                  0.85*textSizeLablesChPionRatio[cent],textSizeLablesChPionRatio[cent], 1,0.257/(textsizeFacXSecChPionRatio[cent]*marginXSec), 510, 508);
        ratio2DNPiToChPi->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DNPiToChPi->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DNPiToChPi->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DNPiToChPi->GetXaxis()->SetLabelFont(42);
        ratio2DNPiToChPi->GetYaxis()->SetLabelFont(42);
        ratio2DNPiToChPi->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DNPiToChPi->GetXaxis()->SetTickLength(0.07);
        ratio2DNPiToChPi->DrawCopy();

        labelEnergyRatioPaper[cent]     = new TLatex(0.935, yPosLabel[cent],Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
        SetStyleTLatex( labelEnergyRatioPaper[cent], textSizeLablesChPionRatio[cent],4, 1, 42, kTRUE, 31);

        if (graphRatioCombPi0ChPiComb[cent]){
            DrawGammaSetMarkerTGraphErr(graphRatioCombPi0ChPiComb[cent], markerStyleCent[cent], markerSizeCent[cent]/2., colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioCombPi0ChPiComb[cent]->Draw("pE,same");
        }
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting,1., 1.,1,kGray+1);
        labelEnergyRatioPaper[cent]->Draw();
        ratio2DNPiToChPi->DrawCopy("axis,same");
        delete ratio2DNPiToChPi;
    }
    canvasRatioChargedPion->Print(Form("%s/NeutralToChargedPion.%s",outputDir.Data(),suffix.Data()));



    for (Int_t cent = 0; cent < maxCent; cent++){
        if (!enableCent[cent]) continue;
        TCanvas* canvasRatioToChargedPi       = new TCanvas("canvasRatioToChargedPi","",200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRatioToChargedPi, 0.08, 0.01, 0.01, 0.125);
        canvasRatioToChargedPi->SetLogx();

        Double_t textsizeLabelsPbPb      = 0;
        if (canvasRatioToChargedPi->XtoPixel(canvasRatioToChargedPi->GetX2()) <canvasRatioToChargedPi->YtoPixel(canvasRatioToChargedPi->GetY1()) ){
            textsizeLabelsPbPb           = (Double_t)textSizeLabelsPixel/canvasRatioToChargedPi->XtoPixel(canvasRatioToChargedPi->GetX2()) ;
        } else {
            textsizeLabelsPbPb           = (Double_t)textSizeLabelsPixel/canvasRatioToChargedPi->YtoPixel(canvasRatioToChargedPi->GetY1());
        }
        cout << textsizeLabelsPbPb << endl;

        TH2F * histo2DPi0RatioToChargedPion;
        histo2DPi0RatioToChargedPion               = new TH2F("histo2DPi0RatioToChargedPion","histo2DPi0RatioToChargedPion",1000,minPtPi0Plotting, maxPtPi0Plotting,1000,0.653,1.95    );
        SetStyleHistoTH2ForGraphs(histo2DPi0RatioToChargedPion, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.85*textsizeLabelsPbPb, textsizeLabelsPbPb,
                                0.85*textsizeLabelsPbPb,textsizeLabelsPbPb, 0.9, 0.65, 510, 505);
        histo2DPi0RatioToChargedPion->GetXaxis()->SetMoreLogLabels();
        histo2DPi0RatioToChargedPion->GetXaxis()->SetLabelOffset(-0.01);
        histo2DPi0RatioToChargedPion->GetYaxis()->SetRangeUser(0.653,1.95);
        histo2DPi0RatioToChargedPion->Draw("copy");

        if (graphRatioCombPi0ChPiComb[cent]){
            ProduceGraphErrWithoutXErrors(graphRatioCombPi0ChPiComb[cent]);
            DrawGammaSetMarkerTGraphErr(graphRatioCombPi0ChPiComb[cent], markerStyleCent[cent], markerSizeCent[cent]/2., colorCent[cent], colorCent[cent], widthLinesBoxes, kTRUE, 0);
            graphRatioCombPi0ChPiComb[cent]->Draw("pE,same");
        }

        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1., 1.,1, kGray+2);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.1, 1.1,1, kGray, 7);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.9, 0.9,1, kGray, 7);

        TLatex *labelRatioToFitEnergy   = new TLatex(0.95, 0.92, Form("%s %s", centArray2[cent].Data(), collisionSystemPbPb.Data()));
        SetStyleTLatex( labelRatioToFitEnergy, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitEnergy->Draw();

        canvasRatioToChargedPi->SaveAs(Form("%s/CombNeutralToChargedPion_%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(), suffix.Data()));

        // **********************************************************************************************************************
        // *******************************************Plot Ratio of Individual meas to Fit ******************************************
        // **********************************************************************************************************************

        canvasRatioToChargedPi->cd();
        histo2DPi0RatioToChargedPion->Draw("copy");

        for (Int_t meth = 10; meth > -1 ; meth--){
            if (graphRatioIndPi0ChPiComb[cent][meth]){
                ProduceGraphErrWithoutXErrors(graphRatioIndPi0ChPiComb[cent][meth]);
                DrawGammaSetMarkerTGraphErr(graphRatioIndPi0ChPiComb[cent][meth], markerStyleDet[meth] ,markerSizeDet[meth]*0.5, colorDet[meth], colorDet[meth]);
                graphRatioIndPi0ChPiComb[cent][meth]->Draw("p,same,z");
            }
        }

        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1., 1.,1, kGray+2);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 1.1, 1.1,1, kGray, 7);
        DrawGammaLines(minPtPi0Plotting, maxPtPi0Plotting , 0.9, 0.9,1, kGray, 7);

        labelRatioToFitEnergy->Draw();
        histo2DPi0RatioToChargedPion->Draw("same,axis");

        canvasRatioToChargedPi->SaveAs(Form("%s/IndNeutralToChargedPion_%s%s.%s",outputDirSupportPaper.Data(), centArrayOutput[cent].Data(), addCentString[cent].Data(),
                                        suffix.Data()));
    }
    // **********************************************************************************************************************
    // ************************* Saving of final results ********************************************************************
    // **********************************************************************************************************************

    TString nameOutputCommonFile    = Form("CombinedResultsPaperPbPb5TeVCent_%s.root", dateForOutput.Data());

    TFile fCombResults(nameOutputCommonFile.Data(), "UPDATE");

    for (Int_t cent = 0; cent < maxCent; cent++){
        if (!enableCent[cent]) continue;
        TDirectoryFile* directoryPi0Output  = NULL;
        directoryPi0Output                  = (TDirectoryFile*)fCombResults.Get(Form("Pi0PbPb5TeV_%s%s", centArray[cent].Data(), addCentString[cent].Data()));
        if (!directoryPi0Output){
            fCombResults.mkdir(Form("Pi0PbPb5TeV_%s%s", centArray[cent].Data(), addCentString[cent].Data()));
            directoryPi0Output              = (TDirectoryFile*)fCombResults.Get(Form("Pi0PbPb5TeV_%s%s", centArray[cent].Data(), addCentString[cent].Data()));
        }
        fCombResults.cd(Form("Pi0PbPb5TeV_%s%s", centArray[cent].Data(), addCentString[cent].Data()));

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
        directoryEtaOutput                  = (TDirectoryFile*)fCombResults.Get(Form("EtaPbPb5TeV_%s%s", centArray[cent].Data(), addCentString[cent].Data()));
        if (!directoryEtaOutput){
            fCombResults.mkdir(Form("EtaPbPb5TeV_%s%s", centArray[cent].Data(), addCentString[cent].Data()));
            directoryEtaOutput              = (TDirectoryFile*)fCombResults.Get(Form("EtaPbPb5TeV_%s%s", centArray[cent].Data(), addCentString[cent].Data()));
        }
        fCombResults.cd(Form("EtaPbPb5TeV_%s%s", centArray[cent].Data(), addCentString[cent].Data()));

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

    TString nameOutputCommonFileFitsOnly    = Form("FitsPaperPbPb5TeVCent_%s.root", dateForOutput.Data());
    TFile fFitsResults(nameOutputCommonFileFitsOnly.Data(), "RECREATE");

        for (Int_t cent = 0; cent < maxCent; cent++){
            if (!enableCentComb[cent]) continue;
            // fits for pi0
            fitInvYieldPi0[cent]->Write(Form("TsallisFitPi0_%s%s",centArray[cent].Data(), addCentString[cent].Data()));
            fitTCMInvYieldPi0[cent]->Write(Form("TwoComponentModelFitPi0_%s%s",centArray[cent].Data(), addCentString[cent].Data()));
            // fits for eta
            fitInvYieldEta[cent]->Write(Form("TsallisFitEta_%s%s",centArray[cent].Data(), addCentString[cent].Data()));
            fitTCMInvYieldEta[cent]->Write(Form("TwoComponentModelFitEta_%s%s",centArray[cent].Data(), addCentString[cent].Data()));
        }
    fFitsResults.Close();

//     // **********************************************************************************************************************
//     // ************************* Saving comparison to comb for diff measurements ********************************************
//     // **********************************************************************************************************************
//
//     TString nameOutputCommonFileCompOnly    = Form("ComparisonsPaperPbPb5TeVCent_%s.root", dateForOutput.Data());
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

