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
extern TBenchmark*    gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*      gMinuit;

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    // TString name;
};

void CombineMesonMeasurementspPb5023GeV_V2(     TString fileNamePCM             = "", 
                                                Bool_t flagPCMfile              = 0,
                                                TString fileNameEMCAL           = "",  
                                                Bool_t flagEMCfile              = 0,
                                                TString fileNamePHOS            = "",
                                                Bool_t flagPHOSfile             = 0,
                                                TString fileNameDalitz          = "",
                                                Bool_t flagDalitzfile           = 0,
                                                TString fileNamePCMEMCAL        = "",
                                                TString suffix                  = "eps", 
                                                TString isMC                    = "", 
                                                TString bWCorrection            = "X",
                                                TString fileNameCorrFactors     = "",
                                                TString fileNameInterpolation   = "",
                                                TString fileConfigRpPbErr       = "",
                                                Bool_t isNSD                    = kTRUE
                                                
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

    TString fileNameTheory                      = "ExternalInputpPb/Theory/TheoryCompilationPPb.root";
    TString fileNameChargedPion                 = "ExternalInputpPb/IdentifiedCharged/ChargedIdentifiedSpectrapPb_2017_01_22.root";
    TString fileNameEtaToPi0WorldDataPP         = "ExternalInput/WorldDataPi0Eta.root";
    TString fileNameEtaToPi0WorldDataPPb        = "ExternalInputpPb/WorldDataPi0EtapPb.root";
    
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMEMCAL.root", fileNamePCMEMCAL.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCAL.root", fileNameEMCAL.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputDalitz.root", fileNameDalitz.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputIdentifiedCharged.root", fileNameChargedPion.Data(), outputDir.Data()));    
    gSystem->Exec(Form("cp %s %s/InputInterpolation.root", fileNameInterpolation.Data(), outputDir.Data()));    
    
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
                                                    "PCM-Dal", "PHOS-Dal", "EMC-Dal", "EMChigh", "mEMC",
                                                    "PCMOtherDataset"};
    TString  nameTrigger[6]                     = {"INT1", "INT7", "EMC1", "EMC7", "EG2", "EG1"};
    
    Color_t  colorDet[11];
    Color_t  colorDetMC[11];
    Style_t  markerStyleDet[11];
    Style_t  markerStyleDetMC[11];
    Size_t   markerSizeDet[11];
    Size_t   markerSizeDetMC[11];

    Double_t nCollpPb                           = GetNCollFromName("", "pPb_5.023TeV");
    Double_t nCollErrpPb                        = GetNCollErrFromName("", "pPb_5.023TeV");
    Double_t tpPb                               = GetTAAFromName("", "pPb_5.023TeV");
    Double_t tpPbErr                            = GetTAAErrFromName("", "pPb_5.023TeV");
                            
    Color_t  colorCGC                           = kRed+2;
    Color_t  colorDSS                           = kAzure+2;
    Color_t  colorDSSBand                       = kAzure-9;
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
    Style_t  styleLineDSS                       = 8;
    Style_t  styleLineNLOMuOne                  = 7;
    Style_t  styleLineNLOMuTwo                  = 4;
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

    //***********************************************************************************************
    //************************** Definition of final pt binning (has to be set manually) ************
    //***********************************************************************************************
    cout << "Setting Pi0 binning" << endl;
    Double_t xPtLimitsPi0[100];
    Int_t maxNBinsPi0               = GetBinning( xPtLimitsPi0, "Pi0", "pPb_5.023TeV", 20 );
    for (Int_t i = 0; i< maxNBinsPi0; i++){
        cout << i << ": "<< xPtLimitsPi0[i] <<" - " << xPtLimitsPi0[i+1]<< ", " <<endl;
    }
    cout << "Setting Eta binning" << endl;
    Double_t xPtLimitsEta[50];
    Int_t maxNBinsEta       = GetBinning( xPtLimitsEta, "Eta", "pPb_5.023TeV", 20 );
//     xPtLimitsEta[maxNBinsEta+1] = 30;   //extend binning by 1 large bin
    for (Int_t i = 0; i< maxNBinsEta; i++){
        cout << i << ": "<< xPtLimitsEta[i] <<" - " << xPtLimitsEta[i+1]<< ", " <<endl;
    }
    cout << endl;
//     return;
    
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
    vector<TString>* ptSysRemNames                  = new vector<TString>[11];
    
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
                if (!(havePi0SysDetailedpPb[nReadSys] && havePi0SysDetailedpp[nReadSys])){
                    havePi0SysDetailedpPb[nReadSys]         = kTRUE;
                    havePi0SysDetailedpp[nReadSys]          = kTRUE;
                    fileNamespPbPi0DetailedSys[nReadSys]    = "";
                    fileNamesppPi0DetailedSys[nReadSys]     = "";
                } else {
                    fileNamesRpPbPi0DetailedSys[nReadSys]   = Form("%s/Pi0RpPb_%s_detailedSys.dat", outputDir.Data(),nameMeasGlobalLabel[nReadSys].Data());
                }    
                if (fileNamespPbEtaDetailedSys[nReadSys].CompareTo("bla") != 0){
                    cout << "pPb sys detailed: " << fileNamespPbEtaDetailedSys[nReadSys] << endl;
                    haveEtaSysDetailedpPb[nReadSys]         = kTRUE;
                } else {
                    fileNamespPbEtaDetailedSys[nReadSys]    = "";
                }    
                if (fileNamesppEtaDetailedSys[nReadSys].CompareTo("bla") != 0){
                    cout << "pPb sys detailed: " << fileNamesppEtaDetailedSys[nReadSys] << endl;
                    haveEtaSysDetailedpp[nReadSys]          = kTRUE;
                } else {
                    fileNamesppEtaDetailedSys[nReadSys]     = "";
                }    
                if (!(haveEtaSysDetailedpPb[nReadSys] && haveEtaSysDetailedpp[nReadSys])){
                    haveEtaSysDetailedpPb[nReadSys]         = kTRUE;
                    haveEtaSysDetailedpp[nReadSys]          = kTRUE;
                    fileNamespPbEtaDetailedSys[nReadSys]    = "";
                    fileNamesppEtaDetailedSys[nReadSys]     = "";
                } else {
                    fileNamesRpPbEtaDetailedSys[nReadSys]   = Form("%s/EtaRpPb_%s_detailedSys.dat", outputDir.Data(),nameMeasGlobalLabel[nReadSys].Data());
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
//     return;
    
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
    TGraphAsymmErrors* graphPCMPi0AccTimesEff       = NULL;
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
    TGraphAsymmErrors* graphDalitzPi0AccTimesEff    = NULL;

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
    for (Int_t i = 0; i< 11; i++){
        statErrorCollectionPi0PP[i]         = NULL;
        systErrorCollectionPi0PP[i]         = NULL;
        systErrorUnCorrCollectionPi0PP[i]   = NULL;
        systErrorInterCollectionPi0PP[i]    = NULL;
        statErrorCollectionEtaPP[i]         = NULL;
        systErrorCollectionEtaPP[i]         = NULL;
        systErrorUnCorrCollectionEtaPP[i]   = NULL;
        systErrorInterCollectionEtaPP[i]    = NULL;
    }    
    
    TFile* fileInterpolation                        = new TFile(fileNameInterpolation.Data());
        for (Int_t i = 0; i< 11; i++){
            statErrorCollectionPi0PP[i]                 = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionStatErr%s_Pi0_5.023TeV",nameMeasGlobalLabel[i].Data()));
            systErrorCollectionPi0PP[i]                 = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionSystErr%s_Pi0_5.023TeV",nameMeasGlobalLabel[i].Data()));
            systErrorUnCorrCollectionPi0PP[i]           = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionUnCorrSystErr%s_Pi0_5.023TeV",nameMeasGlobalLabel[i].Data()));
            systErrorInterCollectionPi0PP[i]            = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionInterpolSystErr%s_Pi0_5.023TeV",nameMeasGlobalLabel[i].Data()));
            if (statErrorCollectionPi0PP[i] && systErrorCollectionPi0PP[i] && systErrorUnCorrCollectionPi0PP[i] && systErrorInterCollectionPi0PP[i])
                haveRefPPPi0[i]                         = kTRUE;
            statErrorCollectionEtaPP[i]                 = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionStatErr%s_Eta_5.023TeV",nameMeasGlobalLabel[i].Data()));
            systErrorCollectionEtaPP[i]                 = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionSystErr%s_Eta_5.023TeV",nameMeasGlobalLabel[i].Data()));
            systErrorUnCorrCollectionEtaPP[i]           = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionUnCorrSystErr%s_Eta_5.023TeV",nameMeasGlobalLabel[i].Data()));
            systErrorInterCollectionEtaPP[i]            = (TGraphAsymmErrors*)fileInterpolation->Get(Form("graphInvXSectionInterpolSystErr%s_Eta_5.023TeV",nameMeasGlobalLabel[i].Data()));
            if (statErrorCollectionEtaPP[i] && systErrorCollectionEtaPP[i] && systErrorUnCorrCollectionEtaPP[i] && systErrorInterCollectionEtaPP[i])
                haveRefPPEta[i]                         = kTRUE;
            if (haveRefPPPi0[i])
                cout << "found pi0 pp reference for " << nameMeasGlobalLabel[i].Data() << endl; 
            if (haveRefPPEta[i])
                cout << "found eta pp reference for " << nameMeasGlobalLabel[i].Data() << endl; 
        }    
        
        TGraphAsymmErrors* graphPPCombPi0Stat           = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionStatErrComb_Pi0_5.023TeV");
        TGraphAsymmErrors* graphPPCombPi0UncorrSys      = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionUnCorrSystErrComb_Pi0_5.023TeV");
        TGraphAsymmErrors* graphPPCombPi0FullSys        = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionSystErrComb_Pi0_5.023TeV");
        TGraphAsymmErrors* graphPPCombPi0InterSys       = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionInterpolSystErrComb_Pi0_5.023TeV");
        TGraphAsymmErrors* graphPPCombEtaStat           = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionStatErrComb_Eta_5.023TeV");
        TGraphAsymmErrors* graphPPCombEtaUncorrSys      = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionUnCorrSystErrComb_Eta_5.023TeV");
        TGraphAsymmErrors* graphPPCombEtaFullSys        = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionSystErrComb_Eta_5.023TeV");
        TGraphAsymmErrors* graphPPCombEtaInterSys       = (TGraphAsymmErrors*)fileInterpolation->Get("graphInvXSectionInterpolSystErrComb_Eta_5.023TeV");
        
    // *******************************************************************************************************
    // ************************** Loading charged pion results ***********************************************
    // *******************************************************************************************************        
    TFile* fileChargedPionInput                         = new TFile(fileNameChargedPion.Data());
        TH1D* histoChPiInvYieldPubStat                  = (TH1D*)fileChargedPionInput->Get("histoChargedPionPubStatpPb");
        TH1D* histoChPiInvYieldPubSyst                  = (TH1D*)fileChargedPionInput->Get("graphChargedPionPubStatpPb");
        TGraphAsymmErrors* graphChPiRpPBPubStat         = (TGraphAsymmErrors*)fileChargedPionInput->Get("graphChargedPionPubStatpPb_RpPb");
        TGraphAsymmErrors* graphChPiRpPBPubSyst         = (TGraphAsymmErrors*)fileChargedPionInput->Get("graphChargedPionPubSyspPb_RpPb");
        TGraphAsymmErrors* graphChKaRpPBPubStat         = (TGraphAsymmErrors*)fileChargedPionInput->Get("graphChargedKaonPubStatpPb_RpPb");
        TGraphAsymmErrors* graphChKaRpPBPubSyst         = (TGraphAsymmErrors*)fileChargedPionInput->Get("graphChargedKaonPubSyspPb_RpPb");
        TGraphAsymmErrors* graphChHadRpPBPubStat        = (TGraphAsymmErrors*)fileChargedPionInput->Get("graphChargedHadronPubStatpPb_RpPb");
        TGraphAsymmErrors* graphChHadRpPBPubSyst        = (TGraphAsymmErrors*)fileChargedPionInput->Get("graphChargedHadronPubSyspPb_RpPb");

    // *******************************************************************************************************
    // ************************** Loading eta/pi0 compilation ************************************************
    // *******************************************************************************************************        
    TFile* fileEtaToPi0CompilationPP                        = new TFile(fileNameEtaToPi0WorldDataPP.Data());
        TGraphErrors* graphPHENIXEtaToPi0200GeV             = (TGraphErrors*)fileEtaToPi0CompilationPP->Get("Phenix200GeV");
        TGraphAsymmErrors* graphALICEEtaToPi07TeV           = (TGraphAsymmErrors*)fileEtaToPi0CompilationPP->Get("Alice7TeV");
        TGraphAsymmErrors* graphALICEEtaToPi02760GeV        = (TGraphAsymmErrors*)fileEtaToPi0CompilationPP->Get("Alice2760GeV");
    TFile* fileEtaToPi0CompilationPPb                       = new TFile(fileNameEtaToPi0WorldDataPPb.Data());
        TGraphErrors* graphCERESEtaToPi0pAu29100MeV         = (TGraphErrors*)fileEtaToPi0CompilationPPb->Get("AgakishievpAu29100MeV");
        TGraphErrors* graphCERESEtaToPi0pBe29100MeV         = (TGraphErrors*)fileEtaToPi0CompilationPPb->Get("AgakishievpBe29100MeV");
        TGraphErrors* graphPHENIXEtaToPi0dAu200GeV          = (TGraphErrors*)fileEtaToPi0CompilationPPb->Get("PhenixdAu200GeV");
      
    // *******************************************************************************************************
    // ************************** Load cocktail inputs *******************************************************
    // *******************************************************************************************************
    TFile* fileCocktail                                     = new TFile("CocktailInput/cocktail_PbPb_0020c_dmtsallis_MTEta.root");
        TDirectory* directoryCocktailMtScaledEta            = (TDirectory*) fileCocktail->Get("cocktail_PbPb_0020c_dmtsallis_MTEta");
        TH1D* cocktailPi0                                   = (TH1D* )directoryCocktailMtScaledEta->Get("ptPi0");
        TH1D* cocktailEta_MtScaled                          = (TH1D* )directoryCocktailMtScaledEta->Get("ptEta");
        cocktailPi0->Sumw2();
        cocktailEta_MtScaled->Sumw2();

        // rebin cocktail input
        TH1D* cocktailPi0_Reb                               = (TH1D* )cocktailPi0->Rebin(maxNBinsEta,"ptPionReb",xPtLimitsEta);
        TH1D* cocktailEta_MtScaledReb                       = (TH1D* )cocktailEta_MtScaled->Rebin(maxNBinsEta,"ptEtaMTScaledRebinned",xPtLimitsEta);

        // calculate eta/pi0 ratio
        TH1D* cocktailEtaToPi0_MtScaled                     = (TH1D* )cocktailEta_MtScaled->Clone("EtaToPi0_MtScaled");
        cocktailEtaToPi0_MtScaled->Sumw2();
        cocktailEtaToPi0_MtScaled->Divide(cocktailEtaToPi0_MtScaled,cocktailPi0);
        TH1D* cocktailEtaToPi0_MtScaledReb                  = (TH1D* )cocktailEta_MtScaledReb->Clone("EtaToPi0_MtScaledRebinned");
        cocktailEtaToPi0_MtScaledReb->Sumw2();
        cocktailEtaToPi0_MtScaledReb->Divide(cocktailEtaToPi0_MtScaledReb,cocktailPi0_Reb);
        
    // *******************************************************************************************************
    // ************************** Loading theory calculations ************************************************
    // *******************************************************************************************************        
    TFile* fileTheoryCompilation                            = new TFile(fileNameTheory.Data());
        TDirectory* directoryMB                             = (TDirectory*) fileTheoryCompilation->Get("pPb_5.023TeV");
        TH1F* histoDPMJetPi0                                = (TH1F*) directoryMB->Get("histoPi0SpecDPMJet5023GeV_Reb");
        TH1F* histoDPMJetEta                                = (TH1F*) directoryMB->Get("histoEtaSpecDPMJet5023GeV_Reb");
        TH1F* histoDPMJetEtaToPi0                           = (TH1F*) directoryMB->Get("histoEtaToPi0DPMJet5023GeV");
        TH1F* histoHIJINGPi0                                = (TH1F*) directoryMB->Get("histoPi0SpecHIJING5023GeV_Reb");
        TH1F* histoHIJINGEta                                = (TH1F*) directoryMB->Get("histoEtaSpecHIJING5023GeV_Reb");
        TH1F* histoHIJINGEtaToPi0                           = (TH1F*) directoryMB->Get("histoEtaToPi0HIJING5023GeV");
        TH1D* histoEPOS3Pi0                                 = (TH1D*) directoryMB->Get("histoPi0SpecEPOS35023GeV_Reb");
        TGraphErrors* graphEPOS3Pi0                         = new TGraphErrors(histoEPOS3Pi0);
        graphEPOS3Pi0->RemovePoint(0);
        TH1D* histoEPOS3Eta                                 = (TH1D*) directoryMB->Get("histoEtaSpecEPOS35023GeV_Reb");
        TGraphErrors* graphEPOS3Eta                         = new TGraphErrors(histoEPOS3Eta);
        graphEPOS3Eta->RemovePoint(0);
        TH1D* histoEPOS3EtaToPi0_Reb                        = (TH1D* )directoryMB->Get("histoEtaToPi0EPOS35023GeV_Reb");
        TGraphErrors* graphEPOS3EtaToPi0_Reb                = new TGraphErrors(histoEPOS3EtaToPi0_Reb);
        graphEPOS3EtaToPi0_Reb->RemovePoint(0);
        TGraphErrors* graphMcGillPi0                        = (TGraphErrors*) directoryMB->Get("graphPi0SpecMcGill5023GeV");
        TGraphErrors* graphMcGillEta                        = (TGraphErrors*) directoryMB->Get("graphEtaSpecMcGill5023GeV");
        TGraphErrors* graphMcGillEtaToPi0                   = (TGraphErrors*) directoryMB->Get("graphEtaToPi0McGill5023GeV");
        
        TGraph* graphPi0RpACGC5023GeV                       = (TGraph*) directoryMB->Get("graphPi0RpACGC5023GeV");  
        TGraph* graphPi0RpAEPS09sKKP5023GeV                 = (TGraph*) directoryMB->Get("graphPi0RpAEPS09sKKP5023GeV");  
        TGraph* graphPi0RpAEPS09sAKK5023GeV                 = (TGraph*) directoryMB->Get("graphPi0RpAEPS09sAKK5023GeV");  
        TGraph* graphPi0RpAEPS09sDSS5023GeV                 = (TGraph*) directoryMB->Get("graphPi0RpAEPS09sDSS5023GeV");
        TGraphAsymmErrors* graphPi0RpAErrEPS09sDSS5023GeV   = (TGraphAsymmErrors*)directoryMB->Get("graphPi0RpAAsymmErrEPS09sDSS5023GeV");
        
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


    TGraphAsymmErrors* statErrorGraphCollectionPi0[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorGraphCollectionPi0[i]   = NULL;
    }    
    statErrorGraphCollectionPi0[0]  = (TGraphAsymmErrors*)graphPCMPi0InvYieldStat->Clone("statErrGraphPCMPi0");
    statErrorGraphCollectionPi0[1]  = (TGraphAsymmErrors*)graphPHOSPi0InvYieldStat->Clone("statErrGraphPHOSPi0");
    statErrorGraphCollectionPi0[2]  = (TGraphAsymmErrors*)graphEMCALPi0InvYieldStat->Clone("statErrGraphEMCALPi0");
    statErrorGraphCollectionPi0[5]  = (TGraphAsymmErrors*)graphDalitzPi0InvYieldStat->Clone("statErrGraphDalitzPi0");
    
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
    Int_t nComBinsPi0Shifting[11]   = { 30, 31, 30, 0,  31,
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
    statErrorGraphCollectionPi0[4]          = (TGraphAsymmErrors*)graphPCMEMCALPi0InvYieldStat->Clone("statErrGraphPCMEMCALPi0");

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
    histo2DPi0Weights = new TH2F("histo2DPi0Weights","histo2DPi0Weights",11000,0.23, 25.,1000,-0.7,1.3);
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

//      DrawGammaLines(0.23, 25. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);
        
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

//      DrawGammaLines(0.23, 25. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);
        
    canvasWeights->SaveAs(Form("%s/Pi0_Weights.%s",outputDir.Data(),suffix.Data()));
    
    //  *********************************************************************************************************************
    //  ************************************ Visualize relative errors ******************************************************
    //  *********************************************************************************************************************
    
    TCanvas* canvasRelSysErr            = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErr->SetLogx();
   
    TH2F * histo2DRelSysErr;
    histo2DRelSysErr                    = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.23, 25.,1000,0,45.5);
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
    histo2DRelStatErr                   = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.23, 25.,1000,0,45.5);
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
    histo2DRelTotErrPi0                 = new TH2F("histo2DRelTotErrPi0","histo2DRelTotErrPi0",11000,0.23, 25.,1000,0,45.5);
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

    TGraphAsymmErrors* graphIndPi0InvYieldStatUnshi[11]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndPi0InvYieldSysUnshi[11]      = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndPi0InvYieldStat[11]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndPi0InvYieldSys[11]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndPi0InvYieldStat_yShifted[11] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndPi0InvYieldSys_yShifted[11]  = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    
    for (Int_t i = 0; i< 11; i++){
        if (statErrorGraphCollectionPi0[i]){
            graphIndPi0InvYieldStatUnshi[i]                 = (TGraphAsymmErrors*)statErrorGraphCollectionPi0[i]->Clone(Form("Pi0UnshiftedStat%s",nameMeasGlobalLabel[i].Data()));
            graphIndPi0InvYieldStat[i]                      = (TGraphAsymmErrors*)statErrorGraphCollectionPi0[i]->Clone(Form("Pi0Stat%s",nameMeasGlobalLabel[i].Data()));
            graphIndPi0InvYieldStat_yShifted[i]             = (TGraphAsymmErrors*)statErrorGraphCollectionPi0[i]->Clone(Form("Pi0YShiftedStat%s",nameMeasGlobalLabel[i].Data()));
        } 
        if (sysErrorCollectionPi0[i]){
            graphIndPi0InvYieldSysUnshi[i]                  = (TGraphAsymmErrors*)sysErrorCollectionPi0[i]->Clone(Form("Pi0UnshiftedSys%s",nameMeasGlobalLabel[i].Data()));
            graphIndPi0InvYieldSys[i]                       = (TGraphAsymmErrors*)sysErrorCollectionPi0[i]->Clone(Form("Pi0Sys%s",nameMeasGlobalLabel[i].Data()));
            graphIndPi0InvYieldSys_yShifted[i]              = (TGraphAsymmErrors*)sysErrorCollectionPi0[i]->Clone(Form("Pi0YShiftedSys%s",nameMeasGlobalLabel[i].Data()));
        }
    }
   
    // fitting spectrum with intial parameters
    // Two component model fit from Bylinkin
    TF1* fitTCMDecomposedLPi0           = FitObject("tcmlow","twoCompModelPi0_DecL", "Pi0", NULL, 0.3, 2.);
    TF1* fitTCMDecomposedHPi0           = FitObject("tcmhigh","twoCompModelPi0_DecH", "Pi0", NULL, 4, 50.);
    fitTCMDecomposedLPi0->SetParameters(graphCombPi0InvYieldTot->GetY()[2],0.3);
    graphCombPi0InvYieldStat->Fit(fitTCMDecomposedLPi0,"QNRMEX0+","",0.3,0.8);
    graphCombPi0InvYieldStat->Fit(fitTCMDecomposedHPi0,"QNRMEX0+","",3,20);
//     fitTCMDecomposedHPi0->SetParameters(graphCombPi0InvYieldTot->GetY()[2],0.8, 2);    
    
    cout << WriteParameterToFile(fitTCMDecomposedLPi0)<< endl;    
    fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedLPi0)<< endl;    
    cout << WriteParameterToFile(fitTCMDecomposedHPi0)<< endl;    
    fileFitsOutput <<  WriteParameterToFile(fitTCMDecomposedHPi0)<< endl;    

    Double_t paramTCMPi0New[5]          = { fitTCMDecomposedLPi0->GetParameter(0),fitTCMDecomposedLPi0->GetParameter(1),
                                            fitTCMDecomposedHPi0->GetParameter(0),fitTCMDecomposedHPi0->GetParameter(1),fitTCMDecomposedHPi0->GetParameter(2)};
    
//     Double_t paramTCMPi0New[5]          = { 5.265,0.33,
//                                             1.9,0.46,3.1};
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
        TF1* fitShiftingPi0            = FitObject("tmpt","ShiftingPi0","Pi0");
        fitShiftingPi0->SetParameters(fitInvYieldPi0->GetParameter(0),fitInvYieldPi0->GetParameter(1), fitInvYieldPi0->GetParameter(2));
//         TF1* fitShiftingPi0                 = FitObject("tcmpt","ShiftingPi0","Pi0");
//         fitShiftingPi0->SetParameters(fitTCMInvYieldPi0->GetParameter(0),fitTCMInvYieldPi0->GetParameter(1), fitTCMInvYieldPi0->GetParameter(2), fitTCMInvYieldPi0->GetParameter(3),fitTCMInvYieldPi0->GetParameter(4));
        
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
        for (Int_t i = 0; i< 11; i++){
            if (graphIndPi0InvYieldStat[i]){
                cout << "shiting stat err of " << nameMeasGlobalLabel[i].Data();
                graphIndPi0InvYieldStat[i]  = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot,
                                                                            graphIndPi0InvYieldStat[i],
                                                                            fitShiftingPi0, 
                                                                            offSetPi0Shifting[i], nComBinsPi0Shifting[i]);
                
            }
            if (graphIndPi0InvYieldSys[i]){
                cout << "shiting sys err of " << nameMeasGlobalLabel[i].Data();
                graphIndPi0InvYieldSys[i]   = ApplyXshiftIndividualSpectra( graphCombPi0InvYieldTot,
                                                                            graphIndPi0InvYieldSys[i],
                                                                            fitShiftingPi0, 
                                                                            offSetPi0Shifting[i], nComBinsPi0Shifting[i]);
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
        histo2DDummy2               = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.33,31.,1000,1e-9,10e1);
        SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
        histo2DDummy2->DrawCopy(); 

            for (Int_t i = 0; i < 11; i++){
                if (graphIndPi0InvYieldSys[i]){
                    DrawGammaSetMarkerTGraphAsym(graphIndPi0InvYieldSys[i], markerStyleDet[i] ,markerSizeDet[i]/2, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
                    graphIndPi0InvYieldSys[i]->Draw("pEsame");
                }    
            }    

            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStatUnshi, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombPi0InvYieldStatUnshi->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombPi0InvYieldStat->Draw("pEsame");

            fitTCMInvYieldPi0->SetLineColor(kBlue+2);
            fitTCMInvYieldPi0->Draw("same");
        
        canvasDummy2->Update();
        canvasDummy2->Print(Form("%s/ComparisonShiftedPi0_PPb5023GeV.%s",outputDir.Data(),suffix.Data()));
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
//     TF1* fitTCMInvYieldPi02  = FitObject("tcm","fitTCMInvYieldPi0pPb5023GeV2","Pi0",graphCombPi0InvYieldTot,0.2, 20.,paramTCMPi0New,"QNRMEX0+","", kFALSE);
//     cout << WriteParameterToFile(fitTCMInvYieldPi02)<< endl;

    TF1* fitPowInvYieldPi0   = FitObject("m","fitPowInvYieldPi0pPb5023GeV","Pi0",graphCombPi0InvYieldTot,5,40. ,NULL,"QNRMEX0+","", kFALSE);
    cout << WriteParameterToFile(fitPowInvYieldPi0)<< endl;
                
    // *************************************************************************************************************
    // Shift graphs in Y direction as well if desired
    // *************************************************************************************************************
    TGraphAsymmErrors* graphCombPi0InvYieldTot_yShifted         = NULL;
    TGraphAsymmErrors* graphCombPi0InvYieldStat_yShifted        = NULL;
    TGraphAsymmErrors* graphCombPi0InvYieldSys_yShifted         = NULL;
        
    if(bWCorrection.Contains("Y") ){
        graphCombPi0InvYieldTot_yShifted        = (TGraphAsymmErrors*)graphCombPi0InvYieldTotUnshi->Clone("Pi0YShiftedCombTot");
        graphCombPi0InvYieldTot_yShifted        =  ApplyYshiftIndividualSpectra( graphCombPi0InvYieldTot_yShifted, fitInvYieldPi0);
        graphCombPi0InvYieldStat_yShifted       = (TGraphAsymmErrors*)graphCombPi0InvYieldStatUnshi->Clone("Pi0YShiftedCombStat");
        graphCombPi0InvYieldStat_yShifted       =  ApplyYshiftIndividualSpectra( graphCombPi0InvYieldStat_yShifted, fitInvYieldPi0);
        graphCombPi0InvYieldSys_yShifted        = (TGraphAsymmErrors*)graphCombPi0InvYieldSysUnshi->Clone("Pi0YShiftedCombSys");
        graphCombPi0InvYieldSys_yShifted        =  ApplyYshiftIndividualSpectra( graphCombPi0InvYieldSys_yShifted, fitInvYieldPi0);
        
        for (Int_t i = 0; i< 11; i++){
            if (graphIndPi0InvYieldStat_yShifted[i]){
                graphIndPi0InvYieldStat_yShifted[i] = ApplyYshiftIndividualSpectra( graphIndPi0InvYieldStat_yShifted[i], fitInvYieldPi0);
            }    
            if (graphIndPi0InvYieldSys_yShifted[i]){
                graphIndPi0InvYieldSys_yShifted[i]  = ApplyYshiftIndividualSpectra( graphIndPi0InvYieldSys_yShifted[i], fitInvYieldPi0);
            }
        }    
    }

    // *************************************************************************************************************
    // Calculate ratios to combined fit
    // *************************************************************************************************************
    TH1D* histoRatioPi0DPMJetToFit                      = (TH1D*) histoDPMJetPi0->Clone("histoRatioPi0DPMJetToFit"); 
    histoRatioPi0DPMJetToFit                            = CalculateHistoRatioToFit (histoRatioPi0DPMJetToFit, fitTCMInvYieldPi0); 
    histoRatioPi0DPMJetToFit->GetXaxis()->SetRangeUser(0.3, 20);
    TH1D* histoRatioPi0HIJINGToFit                      = (TH1D*) histoHIJINGPi0->Clone("histoRatioPi0HIJINGToFit"); 
    histoRatioPi0HIJINGToFit                            = CalculateHistoRatioToFit (histoRatioPi0HIJINGToFit, fitTCMInvYieldPi0); 
    histoRatioPi0HIJINGToFit->GetXaxis()->SetRangeUser(0.3, 20);
    TH1D* histoRatioPi0EPOS3ToFit                       = (TH1D*) histoEPOS3Pi0->Clone("histoRatioPi0EPOS3ToFit"); 
    histoRatioPi0EPOS3ToFit                             = CalculateHistoRatioToFit (histoRatioPi0EPOS3ToFit, fitTCMInvYieldPi0); 
    histoRatioPi0EPOS3ToFit->GetXaxis()->SetRangeUser(0.3, 20);
    TGraphErrors* graphRatioPi0EPOS3ToFit               = new TGraphErrors(histoRatioPi0EPOS3ToFit);
    while(graphRatioPi0EPOS3ToFit->GetX()[0] < 0.3)
        graphRatioPi0EPOS3ToFit->RemovePoint(0);
    while(graphRatioPi0EPOS3ToFit->GetX()[graphRatioPi0EPOS3ToFit->GetN()-1] > 20)
        graphRatioPi0EPOS3ToFit->RemovePoint(graphRatioPi0EPOS3ToFit->GetN()-1);
    TGraphErrors* graphRatioPi0McGillToFit                   = (TGraphErrors*)graphMcGillPi0->Clone("graphRatioPi0McGillToFit");
    graphRatioPi0McGillToFit                                 = CalculateGraphErrRatioToFit(graphRatioPi0McGillToFit, fitTCMInvYieldPi0); 
    while(graphRatioPi0McGillToFit->GetX()[0] < 0.3)
        graphRatioPi0McGillToFit->RemovePoint(0);

// 
    TGraphAsymmErrors* graphRatioPi0IndCombFitStat[11]  = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphRatioPi0IndCombFitSys[11]   = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    
    TGraphAsymmErrors* graphRatioPi0CombCombFitTot      = (TGraphAsymmErrors*)graphCombPi0InvYieldTot->Clone();
    graphRatioPi0CombCombFitTot                         = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitTot, fitTCMInvYieldPi0); 
    TGraphAsymmErrors* graphRatioPi0CombCombFitStat     = (TGraphAsymmErrors*)graphCombPi0InvYieldStat->Clone();
    graphRatioPi0CombCombFitStat                        = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitStat, fitTCMInvYieldPi0); 
    TGraphAsymmErrors* graphRatioPi0CombCombFitSys      = (TGraphAsymmErrors*)graphCombPi0InvYieldSys->Clone();
    graphRatioPi0CombCombFitSys                         = CalculateGraphErrRatioToFit(graphRatioPi0CombCombFitSys, fitTCMInvYieldPi0); 
    
    for (Int_t i= 0; i< 11; i++){
        if (graphIndPi0InvYieldStat[i]){
            graphRatioPi0IndCombFitStat[i]              = (TGraphAsymmErrors*)graphIndPi0InvYieldStat[i]->Clone(Form("RatioPi0%sToCombFitStat", nameMeasGlobalLabel[i].Data()));
            graphRatioPi0IndCombFitStat[i]              = CalculateGraphErrRatioToFit(graphRatioPi0IndCombFitStat[i], fitTCMInvYieldPi0); 
        }    
        if (graphIndPi0InvYieldSys[i]){
            graphRatioPi0IndCombFitSys[i]              = (TGraphAsymmErrors*)graphIndPi0InvYieldSys[i]->Clone(Form("RatioPi0%sToCombFitSyst", nameMeasGlobalLabel[i].Data()));
            graphRatioPi0IndCombFitSys[i]              = CalculateGraphErrRatioToFit(graphRatioPi0IndCombFitSys[i], fitTCMInvYieldPi0); 
        }    
    }

    // **********************************************************************************************************************
    // ******************************************* Plot Ratio of Comb to Fit ****************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 54;
    TCanvas* canvasRatioToCombFit       = new TCanvas("canvasRatioToCombFit","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToCombFit, 0.08, 0.01, 0.01, 0.125);
    canvasRatioToCombFit->SetLogx();

        Double_t textsizeLabelsPPb      = 0;
        Double_t textsizeFacPPb         = 0;
        if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
            textsizeLabelsPPb           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
            textsizeFacPPb              = (Double_t)1./canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
        } else {
            textsizeLabelsPPb           = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
            textsizeFacPPb              = (Double_t)1./canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
        }
        cout << textsizeLabelsPPb << endl;
    
    TH2F * histo2DPi0RatioToCombFit;
    histo2DPi0RatioToCombFit               = new TH2F("histo2DPi0RatioToCombFit","histo2DPi0RatioToCombFit",1000,0.23, 25.,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DPi0RatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsPPb, textsizeLabelsPPb, 
                              0.85*textsizeLabelsPPb,textsizeLabelsPPb, 0.9, 0.65, 510, 505);
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

        DrawGammaLines(0.23, 25. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.23, 25. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.23, 25. , 0.9, 0.9,0.1, kGray, 7);

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
    
        for (Int_t i = 10; i > -1 ; i--){
            if (graphRatioPi0IndCombFitSys[i]){
                DrawGammaSetMarkerTGraphAsym(graphRatioPi0IndCombFitSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
                graphRatioPi0IndCombFitSys[i]->Draw("E2same");
            } 
            if (graphRatioPi0IndCombFitStat[i]){
                ProduceGraphAsymmWithoutXErrors(graphRatioPi0IndCombFitStat[i]);
                DrawGammaSetMarkerTGraphAsym(graphRatioPi0IndCombFitStat[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i]);
                graphRatioPi0IndCombFitStat[i]->Draw("p,same,z");
            }
        }
        graphRatioPi0IndCombFitStat[4]->Draw("p,same,z");

        DrawGammaLines(0.23, 25. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.23, 25. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.23, 25. , 0.9, 0.9,0.5, kGray, 7);
        
        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitPi0->Draw();
    
        //****************************** Definition of the Legend ******************************************
        //**************** Row def ************************
        Double_t rowsLegendOnlyPi0Ratio[6]          = {0.92, 0.867, 0.807, 0.747, 0.687, 0.627};
        Double_t rowsLegendOnlyPi0RatioAbs[6]       = {0.92, 2.48, 2.32, 2.16, 2.00, 1.84 };
//         Double_t rowsLegendOnlyPi0RatioAbs[6]       = {0.92, 2.18, 2.03, 1.87, 1.73, 1.56 };
        Double_t columnsLegendOnlyPi0Ratio[3]       = {0.205, 0.39, 0.47};
        Double_t columnsLegendOnlyPi0RatioAbs[3]    = {0.215, 1.35, 1.83};
//         Double_t columnsLegendOnlyPi0Ratio[3]       = {0.115, 0.325, 0.41};
//         Double_t columnsLegendOnlyPi0RatioAbs[3]    = {0.115, 1., 1.43};
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
        TMarker* markerPCMPi0OnlyRatioPi0           = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[0],columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[1],1);
        markerPCMPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[1]);
        TMarker* markerPHOSPi0OnlyRatioPi0          = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[1], columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[2],1);
        markerPHOSPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[2]);
        TMarker* markerEMCALPi0OnlyRatioPi0         = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[2], columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[3],1);
        markerEMCALPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[3]);
        TMarker* markerPCMEMCALPi0OnlyRatioPi0      = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[4], columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[4],1);
        markerPCMEMCALPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[4]);
        TMarker* markerDalitzPi0OnlyRatioPi0        = CreateMarkerFromGraph(graphRatioPi0IndCombFitSys[5], columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[5],1);
        markerDalitzPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[5]);

        TBox* boxPCMPi0OnlyRatioPi0                 = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[0], columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[1]- heightBox,
                                                                         columnsLegendOnlyPi0RatioAbs[2]+ 5*lengthBox, rowsLegendOnlyPi0RatioAbs[1]+ heightBox);
        boxPCMPi0OnlyRatioPi0->Draw("l");
        TBox* boxPHOSPi0OnlyRatioPi0                = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[1], columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[2]- heightBox,
                                                                         columnsLegendOnlyPi0RatioAbs[2]+ 5*lengthBox, rowsLegendOnlyPi0RatioAbs[2]+ heightBox);
        boxPHOSPi0OnlyRatioPi0->Draw("l");
        TBox* boxEMCALPi0OnlyRatioPi0               = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[2], columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[3]- heightBox,
                                                                         columnsLegendOnlyPi0RatioAbs[2]+ 5*lengthBox, rowsLegendOnlyPi0RatioAbs[3]+ heightBox);
        boxEMCALPi0OnlyRatioPi0->Draw("l");
        TBox* boxPCMEMCALPi0OnlyRatioPi0            = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[4], columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[4]- heightBox,
                                                                         columnsLegendOnlyPi0RatioAbs[2]+ 5*lengthBox, rowsLegendOnlyPi0RatioAbs[4]+ heightBox);
        boxPCMEMCALPi0OnlyRatioPi0->Draw("l");
        TBox* boxDalitzPi0OnlyRatioPi0              = CreateBoxFromGraph(graphRatioPi0IndCombFitSys[5], columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[5]- heightBox,
                                                                         columnsLegendOnlyPi0RatioAbs[2]+ 5*lengthBox, rowsLegendOnlyPi0RatioAbs[5]+ heightBox);
        boxDalitzPi0OnlyRatioPi0->Draw("l");
        
    canvasRatioToCombFit->SaveAs(Form("%s/Pi0_RatioOfIndividualMeasToCombFit_PP.%s",outputDir.Data(),suffix.Data()));
    
    // *******************************************************************************************************
    // *************************** RpPb calculation for pi0 for individual measurements **********************
    // *******************************************************************************************************    
    TGraphAsymmErrors* graphRpPbIndStatPi0[11];
    TH1D* histRpPbIndStatPi0[11];
    TGraphAsymmErrors* graphRpPbIndSystPi0[11];
    TGraphAsymmErrors* systErrorUnCorrInterCollectionPi0PP[11];
    TGraphAsymmErrors* graphRpPbIndCombPi0[11];
    TGraph* graphWeightsPi0RpPb[11];

    for (Int_t i = 0; i < 11; i++){
        graphWeightsPi0RpPb[i]                              = NULL;
        graphRpPbIndStatPi0[i]                              = NULL;
        histRpPbIndStatPi0[i]                               = NULL;
        graphRpPbIndSystPi0[i]                              = NULL;
        graphRpPbIndCombPi0[i]                              = NULL;
        systErrorUnCorrInterCollectionPi0PP[i]              = NULL;
        if (haveRefPPPi0[i] ){
            systErrorUnCorrInterCollectionPi0PP[i]          = (TGraphAsymmErrors*)systErrorInterCollectionPi0PP[i]->Clone(Form("graphPP%sPi0InterUncorrSyst",nameMeasGlobalLabel[i].Data()));
            systErrorUnCorrInterCollectionPi0PP[i]          = AddErrorsQuadraticallyTGraph(systErrorInterCollectionPi0PP[i], systErrorUnCorrCollectionPi0PP[i]);
            if(bWCorrection.Contains("Y") ){
                graphRpPbIndCombPi0[i]                      = CalcRpPbV2(   statErrorCollectionPi0PP[i], 
                                                                            systErrorUnCorrInterCollectionPi0PP[i],
                                                                            graphIndPi0InvYieldStat_yShifted[i],
                                                                            graphIndPi0InvYieldSys_yShifted[i],
                                                                            &graphRpPbIndStatPi0[i],
                                                                            &graphRpPbIndSystPi0[i],
                                                                            tpPb, tpPbErr , 
                                                                            ptSysRemNames[i], fileNamespPbPi0DetailedSys[i], fileNamesppPi0DetailedSys[i], fileNamesRpPbPi0DetailedSys[i]
                                                                        );
                graphRpPbIndCombPi0[i]->Print();
                histRpPbIndStatPi0[i]                       = (TH1D*)statErrorCollectionPi0[i]->Clone(Form("histRpPb%sStatPi0",nameMeasGlobalLabel[i].Data()));
                for (Int_t j = 0; j< histRpPbIndStatPi0[i]->GetNbinsX()+1; j++){
                    histRpPbIndStatPi0[i]->SetBinContent(j,0);
                    histRpPbIndStatPi0[i]->SetBinError(j,0);
                }    
                for (Int_t j = 0; j < graphRpPbIndStatPi0[i]->GetN(); j++){
                    Int_t bin                               = histRpPbIndStatPi0[i]->GetXaxis()->FindBin(graphRpPbIndStatPi0[i]->GetX()[j]);
                    histRpPbIndStatPi0[i]->SetBinContent(bin,graphRpPbIndStatPi0[i]->GetY()[j]);
                    histRpPbIndStatPi0[i]->SetBinError(bin,graphRpPbIndStatPi0[i]->GetEYlow()[j]);
                }
            }    
        }
    }    

    // *******************************************************************************************************
    // *************************** RpPb calculation for pi0 **************************************************
    // *******************************************************************************************************
    
//     TGraphAsymmErrors* graphPPCombPi0InterUncorrSys = (TGraphAsymmErrors*)graphPPCombPi0InterSys->Clone("graphPPCombPi0InterUncorrSys");
//     graphPPCombPi0InterUncorrSys                    = AddErrorsQuadraticallyTGraph(graphPPCombPi0InterSys, graphPPCombPi0UncorrSys);
    TGraphAsymmErrors* graphRpPbCombStatPi0         = NULL;
    TGraphAsymmErrors* graphRpPbCombSystPi0         = NULL;
    TGraphAsymmErrors* graphRpPbCombCombPi0         = NULL;
//     if(bWCorrection.Contains("Y") ){
//         graphRpPbCombCombPi0                        = CalcRpPbV2(   graphPPCombPi0Stat, 
//                                                                     graphPPCombPi0InterUncorrSys, 
//                                                                     graphCombPi0InvYieldStat_yShifted, 
//                                                                     graphCombPi0InvYieldSys_yShifted,
//                                                                     &graphRpPbCombStatPi0,
//                                                                     &graphRpPbCombSystPi0,
//                                                                     tpPb, tpPbErr, ptSysRemNames[0] );
//         graphRpPbCombStatPi0->Print();
//     }    

    if(bWCorrection.Contains("Y") ){
        TString fileNamePi0RpPbOutputWeighting       = Form("%s/Pi0RpPb_WeightingMethod.dat",outputDir.Data());
        graphRpPbCombCombPi0                         = CombinePtPointsSpectraFullCorrMat(   histRpPbIndStatPi0,    graphRpPbIndSystPi0,     
                                                                                            xPtLimitsPi0, maxNBinsPi0,
                                                                                            offSetsPi0, offSetsPi0Sys,
                                                                                            graphRpPbCombStatPi0, graphRpPbCombSystPi0,
                                                                                            fileNamePi0RpPbOutputWeighting, "pPb_5.023TeV", "Pi0RpPb", kTRUE, 
                                                                                            NULL, fileNameCorrFactors
                                                                                            );
        if (graphRpPbCombCombPi0 == NULL) {
            cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
            return;
        }
        
    //     return;
        
        while (graphRpPbCombStatPi0->GetX()[0] < 0.3){
            graphRpPbCombStatPi0->RemovePoint(0);
        }
        while (graphRpPbCombSystPi0->GetX()[0] < 0.3){    
            graphRpPbCombSystPi0->RemovePoint(0);
        }
        while (graphCombPi0InvYieldSys->GetX()[0] < 0.3){    
            graphCombPi0InvYieldSys->RemovePoint(0);
        }    
        graphRpPbCombCombPi0->Print();

        // Reading weights from output file for plotting
        ifstream fileWeightsPi0RpPbRead;
        fileWeightsPi0RpPbRead.open(fileNamePi0RpPbOutputWeighting,ios_base::in);
        cout << "reading" << fileNamePi0RpPbOutputWeighting << endl;
        Double_t xValuesPi0RpPbRead[50];
        Double_t weightsPi0RpPbRead[11][50];
        Int_t availablePi0RpPbMeas[11]      = { -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1,
                                            -1};
        Int_t nMeasSetPi0RpPb               = 5;
        Int_t nPtBinsPi0RpPbRead            = 0;
        while(!fileWeightsPi0RpPbRead.eof() && nPtBinsPi0RpPbRead < 50){
            TString garbage             = "";
            if (nPtBinsPi0RpPbRead == 0){
                fileWeightsPi0RpPbRead >> garbage ;//>> availablePi0RpPbMeas[0] >> availablePi0RpPbMeas[1] >> availablePi0RpPbMeas[2] >> availablePi0RpPbMeas[3];
                for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
                    fileWeightsPi0RpPbRead >> availablePi0RpPbMeas[i] ;
                }    
                cout << "read following measurements: "; 
                for (Int_t i = 0; i < 11; i++){
                    cout << availablePi0RpPbMeas[i] << "\t" ;
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
            graphWeightsPi0RpPb[availablePi0RpPbMeas[i]]                        = new TGraph(nPtBinsPi0RpPbRead,xValuesPi0RpPbRead,weightsPi0RpPbRead[availablePi0RpPbMeas[i]]);
            Int_t bin = 0;
            for (Int_t n = 0; n< nPtBinsPi0RpPbRead; n++){
                if (graphWeightsPi0RpPb[availablePi0RpPbMeas[i]]->GetY()[bin] == 0) graphWeightsPi0RpPb[availablePi0RpPbMeas[i]]->RemovePoint(bin);
                else bin++;
            }    
        }    

        // **********************************************************************************************************************
        // ******************************************* Plotting weights method only EMC *****************************************
        // **********************************************************************************************************************
        textSizeLabelsPixel           = 900*0.04;
        canvasWeights->cd();
        
        histo2DPi0Weights->Draw("copy");
            TLegend* legendWeightsPi0RpPb   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetPi0RpPb), 32);
            for (Int_t i = 0; i < nMeasSetPi0RpPb; i++){
                DrawGammaSetMarkerTGraph(graphWeightsPi0RpPb[availablePi0RpPbMeas[i]], markerStyleDet[availablePi0RpPbMeas[i]], markerSizeDet[availablePi0RpPbMeas[i]]*0.5, colorDet[availablePi0RpPbMeas[i]] , colorDet[availablePi0RpPbMeas[i]]);
                graphWeightsPi0RpPb[availablePi0RpPbMeas[i]]->Draw("p,same,z");
                legendWeightsPi0RpPb->AddEntry(graphWeightsPi0RpPb[availablePi0RpPbMeas[i]],nameMeasGlobalLabel[availablePi0RpPbMeas[i]],"p");
            }    
            legendWeightsPi0RpPb->Draw();

            labelWeightsEnergy->Draw();
            TLatex *labelWeightsPi0RpPb         = new TLatex(0.7,0.16,"R_{pPb}: #pi^{0} #rightarrow #gamma#gamma");
            SetStyleTLatex( labelWeightsPi0RpPb, 0.85*textSizeLabelsPixel, 4, 1, 43 );
            labelWeightsPi0RpPb->Draw();

    //      DrawGammaLines(0.23, 25. , 0.8, 0.8,0.1, kGray, 3);
            DrawGammaLines(0.23, 25. , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(0.23, 25. , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(0.23, 25. , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(0.23, 25. , 0.2, 0.2,0.1, kGray, 3);
            
        canvasWeights->SaveAs(Form("%s/Pi0RpPb_Weights.%s",outputDir.Data(),suffix.Data()));
        
    
    }

    
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

    TGraphAsymmErrors* statErrorGraphCollectionEta[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorGraphCollectionEta[i]   = NULL;
    }    
    statErrorGraphCollectionEta[0]  = (TGraphAsymmErrors*)graphPCMEtaInvYieldStat->Clone("statErrGraphPCMEta");
    statErrorGraphCollectionEta[2]  = (TGraphAsymmErrors*)graphEMCALEtaInvYieldStat->Clone("statErrGraphEMCALEta");
   
    
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

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match for Eta                    
    Int_t offSetsEta[11]            = { 0,  0,  0,  0,  0,
                                        0,  0,  0,  0,  0, 
                                        0};
    Int_t offSetsEtaSys[11]         = { 3,  0,  7,  0,  5,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t offSetEtaShifting[11]     = { 0,  0,  4,  0,  2,
                                        0,  0,  0,  0,  0,
                                        0};
    Int_t nComBinsEtaShifting[11]   = { 13,  0,  17,  0,  15,
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
    statErrorGraphCollectionEta[4]          = (TGraphAsymmErrors*)graphPCMEMCALEtaInvYieldStat->Clone("statErrGraphPCMEMCALEta");
    
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
    histo2DEtaWeights = new TH2F("histo2DEtaWeights","histo2DEtaWeights",11000,0.43, 25.,1000,-0.7,1.3);
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

//      DrawGammaLines(0.43, 25. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.43, 25. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.43, 25. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.43, 25. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.43, 25. , 0.2, 0.2,0.1, kGray, 3);
        
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

//      DrawGammaLines(0.43, 25. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.43, 25. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.43, 25. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.43, 25. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.43, 25. , 0.2, 0.2,0.1, kGray, 3);
        
    canvasWeights->SaveAs(Form("%s/Eta_Weights.%s",outputDir.Data(),suffix.Data()));

    
//     return;
    // *********************************************************************************************************************
    // ************************************ Visualize relative errors Eta ******************************************************
    // *********************************************************************************************************************
    
    canvasRelSysErr->cd();
   
    TH2F * histo2DRelSysErrEta;
    histo2DRelSysErrEta                 = new TH2F("histo2DRelSysErrEta","histo2DRelSysErrEta",11000,0.43, 25.,1000,0,60.5);
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
    histo2DRelStatErrEta                = new TH2F("histo2DRelStatErrEta","histo2DRelStatErrEta",11000,0.43, 25.,1000,0,60.5);
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
    histo2DRelTotErrEta                 = new TH2F("histo2DRelTotErrEta","histo2DRelTotErrEta",11000,0.43, 25.,1000,0,60.5);
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

    TGraphAsymmErrors* graphIndEtaInvYieldStatUnshi[11]     = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndEtaInvYieldSysUnshi[11]      = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndEtaInvYieldStat[11]          = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndEtaInvYieldSys[11]           = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndEtaInvYieldStat_yShifted[11] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphIndEtaInvYieldSys_yShifted[11]  = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

    for (Int_t i = 0; i< 11; i++){
        if (statErrorGraphCollectionEta[i]){
            graphIndEtaInvYieldStatUnshi[i]                 = (TGraphAsymmErrors*)statErrorGraphCollectionEta[i]->Clone(Form("EtaUnshiftedStat%s",nameMeasGlobalLabel[i].Data()));
            graphIndEtaInvYieldStat[i]                      = (TGraphAsymmErrors*)statErrorGraphCollectionEta[i]->Clone(Form("EtaStat%s",nameMeasGlobalLabel[i].Data()));
            graphIndEtaInvYieldStat_yShifted[i]             = (TGraphAsymmErrors*)statErrorGraphCollectionEta[i]->Clone(Form("EtaYShiftedStat%s",nameMeasGlobalLabel[i].Data()));
        } 
        if (sysErrorCollectionEta[i]){
            graphIndEtaInvYieldSysUnshi[i]                  = (TGraphAsymmErrors*)sysErrorCollectionEta[i]->Clone(Form("EtaUnshiftedSys%s",nameMeasGlobalLabel[i].Data()));
            graphIndEtaInvYieldSys[i]                       = (TGraphAsymmErrors*)sysErrorCollectionEta[i]->Clone(Form("EtaSys%s",nameMeasGlobalLabel[i].Data()));
            graphIndEtaInvYieldSys_yShifted[i]              = (TGraphAsymmErrors*)sysErrorCollectionEta[i]->Clone(Form("EtaYShiftedSys%s",nameMeasGlobalLabel[i].Data()));
        }
    }
   
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
        
        for (Int_t i = 0; i< 11; i++){
            if (graphIndEtaInvYieldStat[i]){
                cout << "shiting stat err of " << nameMeasGlobalLabel[i].Data();
                graphIndEtaInvYieldStat[i]  = ApplyXshiftIndividualSpectra( graphCombEtaInvYieldTot,
                                                                            graphIndEtaInvYieldStat[i],
                                                                            fitShiftingEta, 
                                                                            offSetEtaShifting[i], nComBinsEtaShifting[i], "Eta");
                
            }
            if (graphIndEtaInvYieldSys[i]){
                cout << "shiting sys err of " << nameMeasGlobalLabel[i].Data();
                graphIndEtaInvYieldSys[i]   = ApplyXshiftIndividualSpectra( graphCombEtaInvYieldTot,
                                                                            graphIndEtaInvYieldSys[i],
                                                                            fitShiftingEta, 
                                                                            offSetEtaShifting[i], nComBinsEtaShifting[i], "Eta");
            }
        }
        
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
    
            for (Int_t i = 0; i < 11; i++){
                if (graphIndEtaInvYieldSys[i]){
                    DrawGammaSetMarkerTGraphAsym(graphIndEtaInvYieldSys[i], markerStyleDet[i] ,markerSizeDet[i]/2, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
                    graphIndEtaInvYieldSys[i]->Draw("pEsame");
                }    
            }    

            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStatUnshi, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
            graphCombEtaInvYieldStatUnshi->Draw("pEsame");
            DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldTot, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
            graphCombEtaInvYieldTot->Draw("pEsame");
            
//         fitShiftingEta->SetLineColor(kRed+2);
//         fitShiftingEta->SetLineStyle(7);
            fitInvYieldEta->SetLineColor(kBlue+2);
            fitInvYieldEta->Draw("same");
//         fitShiftingEta->Draw("same");
        
//         canvasDummy2->Update();
        canvasDummy2->Print(Form("%s/ComparisonShiftedEta_PPb5023GeV.%s",outputDir.Data(),suffix.Data()));
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
        
    if(bWCorrection.Contains("Y") ){
        graphCombEtaInvYieldTot_yShifted        = (TGraphAsymmErrors*)graphCombEtaInvYieldTotUnshi->Clone("EtaYShiftedCombTot");
        graphCombEtaInvYieldTot_yShifted        =  ApplyYshiftIndividualSpectra( graphCombEtaInvYieldTot_yShifted, fitInvYieldEta);
        graphCombEtaInvYieldStat_yShifted       = (TGraphAsymmErrors*)graphCombEtaInvYieldStatUnshi->Clone("EtaYShiftedCombStat");
        graphCombEtaInvYieldStat_yShifted       =  ApplyYshiftIndividualSpectra( graphCombEtaInvYieldStat_yShifted, fitInvYieldEta);
        graphCombEtaInvYieldSys_yShifted        = (TGraphAsymmErrors*)graphCombEtaInvYieldSysUnshi->Clone("EtaYShiftedCombSys");
        graphCombEtaInvYieldSys_yShifted        =  ApplyYshiftIndividualSpectra( graphCombEtaInvYieldSys_yShifted, fitInvYieldEta);

        for (Int_t i = 0; i< 11; i++){
            if (graphIndEtaInvYieldStat_yShifted[i]){
                graphIndEtaInvYieldStat_yShifted[i] = ApplyYshiftIndividualSpectra( graphIndEtaInvYieldStat_yShifted[i], fitInvYieldEta);
            }    
            if (graphIndEtaInvYieldSys_yShifted[i]){
                graphIndEtaInvYieldSys_yShifted[i]  = ApplyYshiftIndividualSpectra( graphIndEtaInvYieldSys_yShifted[i], fitInvYieldEta);
            }
        }
    }

    // *************************************************************************************************************
    // Calculate ratios to combined fit
    // *************************************************************************************************************    
    TH1D* histoRatioEtaDPMJetToFit                      = (TH1D*) histoDPMJetEta->Clone("histoRatioEtaDPMJetToFit"); 
    histoRatioEtaDPMJetToFit                            = CalculateHistoRatioToFit (histoRatioEtaDPMJetToFit, fitTCMInvYieldEta); 
    histoRatioEtaDPMJetToFit->GetXaxis()->SetRangeUser(0.6,20);
    TH1D* histoRatioEtaHIJINGToFit                      = (TH1D*) histoHIJINGEta->Clone("histoRatioEtaHIJINGToFit"); 
    histoRatioEtaHIJINGToFit                            = CalculateHistoRatioToFit (histoRatioEtaHIJINGToFit, fitTCMInvYieldEta); 
    histoRatioEtaHIJINGToFit->GetXaxis()->SetRangeUser(0.6,20);
    TH1D* histoRatioEtaEPOS3ToFit                       = (TH1D*) histoEPOS3Eta->Clone("histoRatioEtaEPOS3ToFit"); 
    histoRatioEtaEPOS3ToFit                             = CalculateHistoRatioToFit (histoRatioEtaEPOS3ToFit, fitTCMInvYieldEta); 
    histoRatioEtaEPOS3ToFit->GetXaxis()->SetRangeUser(0.6,20);
    TGraphErrors* graphRatioEtaEPOS3ToFit               = new TGraphErrors(histoRatioEtaEPOS3ToFit);
    while(graphRatioEtaEPOS3ToFit->GetX()[0] < 0.6)
        graphRatioEtaEPOS3ToFit->RemovePoint(0);
    while(graphRatioEtaEPOS3ToFit->GetX()[graphRatioEtaEPOS3ToFit->GetN()-1] > 20)
        graphRatioEtaEPOS3ToFit->RemovePoint(graphRatioEtaEPOS3ToFit->GetN()-1);
    TGraphErrors* graphRatioEtaMcGillToFit                   = (TGraphErrors*)graphMcGillEta->Clone("graphRatioEtaMcGillToFit");
    graphRatioEtaMcGillToFit                                 = CalculateGraphErrRatioToFit(graphRatioEtaMcGillToFit, fitTCMInvYieldEta); 
    while(graphRatioEtaMcGillToFit->GetX()[0] < 0.6)
        graphRatioEtaMcGillToFit->RemovePoint(0);
    
    
//     TGraph* graphRatioEtaCombNLOMuHalf                  = (TGraph*)graphNLOCalcEtaMuHalf->Clone();
//     TGraph* graphRatioEtaCombNLOMuOne                   = (TGraph*)graphNLOCalcEtaMuOne->Clone();
//     TGraph* graphRatioEtaCombNLOMuTwo                   = (TGraph*)graphNLOCalcEtaMuTwo->Clone();
//     graphRatioEtaCombNLOMuHalf                          = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuHalf, fitTCMInvYieldEta); 
//     graphRatioEtaCombNLOMuOne                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuOne, fitTCMInvYieldEta); 
//     graphRatioEtaCombNLOMuTwo                           = CalculateGraphRatioToFit (graphRatioEtaCombNLOMuTwo, fitTCMInvYieldEta); 
    
    TGraphAsymmErrors* graphRatioEtaCombCombFitTot      = (TGraphAsymmErrors*)graphCombEtaInvYieldTot->Clone();
    graphRatioEtaCombCombFitTot                         = CalculateGraphErrRatioToFit(graphRatioEtaCombCombFitTot, fitTCMInvYieldEta); 
    TGraphAsymmErrors* graphRatioEtaCombCombFitStat     = (TGraphAsymmErrors*)graphCombEtaInvYieldStat->Clone();
    graphRatioEtaCombCombFitStat                        = CalculateGraphErrRatioToFit(graphRatioEtaCombCombFitStat, fitTCMInvYieldEta); 
    TGraphAsymmErrors* graphRatioEtaCombCombFitSys      = (TGraphAsymmErrors*)graphCombEtaInvYieldSys->Clone();
    graphRatioEtaCombCombFitSys                         = CalculateGraphErrRatioToFit(graphRatioEtaCombCombFitSys, fitTCMInvYieldEta); 
    
    TGraphAsymmErrors* graphRatioEtaIndCombFitStat[11]  = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphRatioEtaIndCombFitSys[11]   = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    
    for (Int_t i= 0; i< 11; i++){
        if (graphIndEtaInvYieldStat[i]){
            graphRatioEtaIndCombFitStat[i]              = (TGraphAsymmErrors*)graphIndEtaInvYieldStat[i]->Clone(Form("RatioEta%sToCombFitStat", nameMeasGlobalLabel[i].Data()));
            graphRatioEtaIndCombFitStat[i]              = CalculateGraphErrRatioToFit(graphRatioEtaIndCombFitStat[i], fitTCMInvYieldEta); 
        }    
        if (graphIndEtaInvYieldSys[i]){
            graphRatioEtaIndCombFitSys[i]              = (TGraphAsymmErrors*)graphIndEtaInvYieldSys[i]->Clone(Form("RatioEta%sToCombFitSyst", nameMeasGlobalLabel[i].Data()));
            graphRatioEtaIndCombFitSys[i]              = CalculateGraphErrRatioToFit(graphRatioEtaIndCombFitSys[i], fitTCMInvYieldEta); 
        }    
    }

    // **********************************************************************************************************************
    // ******************************************* Ratio of Comb to Fit ****************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 48;
    canvasRatioToCombFit->cd();
    TH2F * histo2DEtaRatioToCombFit;
    histo2DEtaRatioToCombFit               = new TH2F("histo2DEtaRatioToCombFit","histo2DEtaRatioToCombFit",1000,0.43, 25.,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DEtaRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.85*textsizeLabelsPPb, textsizeLabelsPPb, 
                              0.85*textsizeLabelsPPb,textsizeLabelsPPb, 0.9, 0.65, 510, 505);
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

        DrawGammaLines(0.43, 25. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.43, 25. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.43, 25. , 0.9, 0.9,0.1, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        TLatex *labelRatioToFitEta      = new TLatex(0.95, 0.82, "#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitEta, textSizeLabelsPixel, 4, 1, 43, kTRUE, 31);
        labelRatioToFitEta->Draw();

    canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfCombToCombFit_PPb5023GeV.%s",outputDir.Data(),suffix.Data()));
    
    // **********************************************************************************************************************
    // ******************************************* Ratio of Individual meas to Fit ******************************************
    // **********************************************************************************************************************
    
    canvasRatioToCombFit->cd();
    histo2DEtaRatioToCombFit->Draw("copy");

        for (Int_t i = 0; i < 11 ; i++){
            if (graphRatioEtaIndCombFitSys[i]){
                DrawGammaSetMarkerTGraphAsym(graphRatioEtaIndCombFitSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
                graphRatioEtaIndCombFitSys[i]->Draw("E2same");
            } 
            if (graphRatioEtaIndCombFitStat[i]){
                ProduceGraphAsymmWithoutXErrors(graphRatioEtaIndCombFitStat[i]);
                DrawGammaSetMarkerTGraphAsym(graphRatioEtaIndCombFitStat[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i]);
                graphRatioEtaIndCombFitStat[i]->Draw("p,same,z");
            }
        }

        DrawGammaLines(0.43, 25. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.43, 25. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.43, 25. , 0.9, 0.9,0.5, kGray, 7);
        
        labelRatioToFitEnergy->Draw();
        labelRatioToFitALICE->Draw();
        labelRatioToFitEta->Draw();
    
        //****************************** Definition of the Legend ******************************************
        Double_t rowsLegendOnlyEtaRatio[4]          = {0.91, 0.86, 0.81, 0.76};
        Double_t rowsLegendOnlyEtaRatioAbs[4]       = {0.91, 2.245, 2.11, 1.975};
        Double_t columnsLegendOnlyEtaRatio[3]       = {0.115, 0.3, 0.375};
        Double_t columnsLegendOnlyEtaRatioAbs[3]    = {0.115, 1.3, 1.68};

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
        TMarker* markerPCMEtaOnlyRatioEta           = CreateMarkerFromGraph(graphRatioEtaIndCombFitSys[0],columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[1],1);
        markerPCMEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[1]);
        TMarker* markerEMCALEtaOnlyRatioEta         = CreateMarkerFromGraph(graphRatioEtaIndCombFitSys[2], columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[2],1);
        markerEMCALEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[2]);
        TMarker* markerPCMEMCALEtaOnlyRatioEta      = CreateMarkerFromGraph(graphRatioEtaIndCombFitSys[4], columnsLegendOnlyEtaRatio[1] ,rowsLegendOnlyEtaRatio[3],1);
        markerPCMEMCALEtaOnlyRatioEta->DrawMarker(columnsLegendOnlyEtaRatioAbs[1] ,rowsLegendOnlyEtaRatioAbs[3]);

        TBox* boxPCMEtaOnlyRatioEta                 = CreateBoxFromGraph(graphRatioEtaIndCombFitSys[0], columnsLegendOnlyEtaRatioAbs[2]-0.5*lengthBox , rowsLegendOnlyEtaRatioAbs[1]- heightBox,
                                                                         columnsLegendOnlyEtaRatioAbs[2]+ 4*lengthBox, rowsLegendOnlyEtaRatioAbs[1]+ heightBox);
        boxPCMEtaOnlyRatioEta->Draw("l");
        TBox* boxEMCALEtaOnlyRatioEta               = CreateBoxFromGraph(graphRatioEtaIndCombFitSys[2], columnsLegendOnlyEtaRatioAbs[2]-0.5*lengthBox , rowsLegendOnlyEtaRatioAbs[2]- heightBox,
                                                                         columnsLegendOnlyEtaRatioAbs[2]+ 4*lengthBox, rowsLegendOnlyEtaRatioAbs[2]+ heightBox);
        boxEMCALEtaOnlyRatioEta->Draw("l");
        TBox* boxPCMEMCALEtaOnlyRatioEta            = CreateBoxFromGraph(graphRatioEtaIndCombFitSys[4], columnsLegendOnlyEtaRatioAbs[2]-0.5*lengthBox , rowsLegendOnlyEtaRatioAbs[3]- heightBox,
                                                                         columnsLegendOnlyEtaRatioAbs[2]+ 4*lengthBox, rowsLegendOnlyEtaRatioAbs[3]+ heightBox);
        boxPCMEMCALEtaOnlyRatioEta->Draw("l");
        
    canvasRatioToCombFit->SaveAs(Form("%s/Eta_RatioOfIndividualMeasToCombFit_PP.%s",outputDir.Data(),suffix.Data()));
    

    
    // *******************************************************************************************************
    // *************************** RpPb calculation for eta for individual measurements **********************
    // *******************************************************************************************************    
    TGraphAsymmErrors* graphRpPbIndStatEta[11];
    TH1D* histRpPbIndStatEta[11];
    TGraphAsymmErrors* graphRpPbIndSystEta[11];
    TGraphAsymmErrors* systErrorUnCorrInterCollectionEtaPP[11];
    TGraphAsymmErrors* graphRpPbIndCombEta[11];
    TGraph* graphWeightsEtaRpPb[11];
    for (Int_t i = 0; i < 11; i++){
        graphWeightsEtaRpPb[i]                              = NULL;
        graphRpPbIndStatEta[i]                              = NULL;
        graphRpPbIndSystEta[i]                              = NULL;
        graphRpPbIndCombEta[i]                              = NULL;
        systErrorUnCorrInterCollectionEtaPP[i]              = NULL;
        histRpPbIndStatEta[i]                               = NULL;
        if (haveRefPPEta[i] && graphIndEtaInvYieldStat_yShifted[i]){
            cout << "running eta for: " << nameMeasGlobalLabel[i].Data() << endl;
            systErrorUnCorrInterCollectionEtaPP[i]          = (TGraphAsymmErrors*)systErrorInterCollectionEtaPP[i]->Clone(Form("graphPP%sEtaInterUncorrSyst",nameMeasGlobalLabel[i].Data()));
            systErrorUnCorrInterCollectionEtaPP[i]          = AddErrorsQuadraticallyTGraph(systErrorInterCollectionEtaPP[i], systErrorUnCorrCollectionEtaPP[i]);
            if(bWCorrection.Contains("Y") ){
                graphRpPbIndCombEta[i]                      = CalcRpPbV2(   statErrorCollectionEtaPP[i], 
                                                                            systErrorUnCorrInterCollectionEtaPP[i],
                                                                            graphIndEtaInvYieldStat_yShifted[i],
                                                                            graphIndEtaInvYieldSys_yShifted[i],
                                                                            &graphRpPbIndStatEta[i],
                                                                            &graphRpPbIndSystEta[i],
                                                                            tpPb, tpPbErr,
                                                                            ptSysRemNames[i], fileNamespPbEtaDetailedSys[i], fileNamesppEtaDetailedSys[i], fileNamesRpPbEtaDetailedSys[i]
                                                                         );
                graphRpPbIndCombEta[i]->Print();
                histRpPbIndStatEta[i]                       = (TH1D*)statErrorCollectionEta[i]->Clone(Form("histRpPb%sStatEta",nameMeasGlobalLabel[i].Data()));
                for (Int_t j = 0; j< histRpPbIndStatEta[i]->GetNbinsX()+1; j++){
                    histRpPbIndStatEta[i]->SetBinContent(j,0);
                    histRpPbIndStatEta[i]->SetBinError(j,0);
                }    
                for (Int_t j = 0; j < graphRpPbIndStatEta[i]->GetN(); j++){
                    Int_t bin                               = histRpPbIndStatEta[i]->GetXaxis()->FindBin(graphRpPbIndStatEta[i]->GetX()[j]);
                    histRpPbIndStatEta[i]->SetBinContent(bin,graphRpPbIndStatEta[i]->GetY()[j]);
                    histRpPbIndStatEta[i]->SetBinError(bin,graphRpPbIndStatEta[i]->GetEYlow()[j]);
                }
            }
            
        }
    }
    // *******************************************************************************************************
    // *************************** RpPb calculation for eta **************************************************
    // *******************************************************************************************************
//     TGraphAsymmErrors* graphPPCombEtaInterUncorrSys = (TGraphAsymmErrors*)graphPPCombEtaInterSys->Clone("graphPPCombEtaInterUncorrSys");
//     graphPPCombEtaInterUncorrSys                    = AddErrorsQuadraticallyTGraph(graphPPCombEtaInterSys, graphPPCombEtaUncorrSys);
    TGraphAsymmErrors* graphRpPbCombStatEta         = NULL;
    TGraphAsymmErrors* graphRpPbCombSystEta         = NULL;
    TGraphAsymmErrors* graphRpPbCombCombEta         = NULL;
//     if(bWCorrection.Contains("Y") ){
//         graphRpPbCombCombEta                            = CalcRpPbV2(   graphPPCombEtaStat, 
//                                                                         graphPPCombEtaInterUncorrSys, 
//                                                                         graphCombEtaInvYieldStat_yShifted, 
//                                                                         graphCombEtaInvYieldSys_yShifted,
//                                                                         &graphRpPbCombStatEta,
//                                                                         &graphRpPbCombSystEta,
//                                                                         tpPb, tpPbErr, ptSysRemNames[0] );
//         graphRpPbCombStatEta->Print();
//     }    

    if(bWCorrection.Contains("Y") ){
        TString fileNameEtaRpPbOutputWeighting       = Form("%s/EtaRpPb_WeightingMethod.dat",outputDir.Data());
        graphRpPbCombCombEta                         = CombinePtPointsSpectraFullCorrMat(   histRpPbIndStatEta,    graphRpPbIndSystEta,     
                                                                                            xPtLimitsEta, maxNBinsEta,
                                                                                            offSetsEta, offSetsEtaSys,
                                                                                            graphRpPbCombStatEta, graphRpPbCombSystEta,
                                                                                            fileNameEtaRpPbOutputWeighting, "pPb_5.023TeV", "EtaRpPb",  kFALSE, 
                                                                                            NULL, fileNameCorrFactors
                                                                                        );
        if (graphRpPbCombCombEta == NULL) {
            cout << "Aborting: something went wrong during the combination of the new spectra" << endl;
            return;
        }
        
    //     return;
        
        while (graphRpPbCombStatEta->GetX()[0] < 0.7){
            graphRpPbCombStatEta->RemovePoint(0);
        }
        while (graphRpPbCombSystEta->GetX()[0] < 0.7){    
            graphRpPbCombSystEta->RemovePoint(0);
        }
        while (graphCombEtaInvYieldSys->GetX()[0] < 0.7){    
            graphCombEtaInvYieldSys->RemovePoint(0);
        }    
        graphRpPbCombCombEta->Print();

        // Reading weights from output file for plotting
        ifstream fileWeightsEtaRpPbRead;
        fileWeightsEtaRpPbRead.open(fileNameEtaRpPbOutputWeighting,ios_base::in);
        cout << "reading" << fileNameEtaRpPbOutputWeighting << endl;
        Double_t xValuesEtaRpPbRead[50];
        Double_t weightsEtaRpPbRead[11][50];
        Int_t availableEtaRpPbMeas[11]      = { -1, -1, -1, -1, -1,
                                            -1, -1, -1, -1, -1,
                                            -1};
        Int_t nMeasSetEtaRpPb               = 3;
        Int_t nPtBinsEtaRpPbRead            = 0;
        while(!fileWeightsEtaRpPbRead.eof() && nPtBinsEtaRpPbRead < 50){
            TString garbage             = "";
            if (nPtBinsEtaRpPbRead == 0){
                fileWeightsEtaRpPbRead >> garbage ;//>> availableEtaRpPbMeas[0] >> availableEtaRpPbMeas[1] >> availableEtaRpPbMeas[2] >> availableEtaRpPbMeas[3];
                for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
                    fileWeightsEtaRpPbRead >> availableEtaRpPbMeas[i] ;
                }    
                cout << "read following measurements: "; 
                for (Int_t i = 0; i < 11; i++){
                    cout << availableEtaRpPbMeas[i] << "\t" ;
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
            graphWeightsEtaRpPb[availableEtaRpPbMeas[i]]                        = new TGraph(nPtBinsEtaRpPbRead,xValuesEtaRpPbRead,weightsEtaRpPbRead[availableEtaRpPbMeas[i]]);
            Int_t bin = 0;
            for (Int_t n = 0; n< nPtBinsEtaRpPbRead; n++){
                if (graphWeightsEtaRpPb[availableEtaRpPbMeas[i]]->GetY()[bin] == 0) graphWeightsEtaRpPb[availableEtaRpPbMeas[i]]->RemovePoint(bin);
                else bin++;
            }    
        }    

        // **********************************************************************************************************************
        // ******************************************* Plotting weights method only EMC *****************************************
        // **********************************************************************************************************************
        textSizeLabelsPixel           = 900*0.04;
        canvasWeights->cd();
        
        histo2DEtaWeights->Draw("copy");
            TLegend* legendWeightsEtaRpPb   = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSetEtaRpPb), 32);
            for (Int_t i = 0; i < nMeasSetEtaRpPb; i++){
                DrawGammaSetMarkerTGraph(graphWeightsEtaRpPb[availableEtaRpPbMeas[i]], markerStyleDet[availableEtaRpPbMeas[i]], markerSizeDet[availableEtaRpPbMeas[i]]*0.5, colorDet[availableEtaRpPbMeas[i]] , colorDet[availableEtaRpPbMeas[i]]);
                graphWeightsEtaRpPb[availableEtaRpPbMeas[i]]->Draw("p,same,z");
                legendWeightsEtaRpPb->AddEntry(graphWeightsEtaRpPb[availableEtaRpPbMeas[i]],nameMeasGlobalLabel[availableEtaRpPbMeas[i]],"p");
            }    
            legendWeightsEtaRpPb->Draw();

            labelWeightsEnergy->Draw();
            TLatex *labelWeightsEtaRpPb         = new TLatex(0.7,0.16,"R_{pPb}: #eta #rightarrow #gamma#gamma");
            SetStyleTLatex( labelWeightsEtaRpPb, 0.85*textSizeLabelsPixel, 4, 1, 43 );
            labelWeightsEtaRpPb->Draw();

    //      DrawGammaLines(0.23, 25. , 0.8, 0.8,0.1, kGray, 3);
            DrawGammaLines(0.43, 25. , 0.5, 0.5,0.1, kGray, 7);
            DrawGammaLines(0.43, 25. , 0.4, 0.4,0.1, kGray, 1);
            DrawGammaLines(0.43, 25. , 0.3, 0.3,0.1, kGray, 7);
            DrawGammaLines(0.43, 25. , 0.2, 0.2,0.1, kGray, 3);
            
        canvasWeights->SaveAs(Form("%s/EtaRpPb_Weights.%s",outputDir.Data(),suffix.Data()));
    }
    
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
    Int_t offSetsEtaToPi0Sys[11]        = { 3,  0,  7,  0,  5,
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
    histo2DEtaToPi0Weights = new TH2F("histo2DEtaToPi0Weights","histo2DEtaToPi0Weights",11000,0.43, 25.,1000,-0.7,1.3);
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

//      DrawGammaLines(0.43, 25. , 0.8, 0.8,0.1, kGray, 3);
        DrawGammaLines(0.43, 25. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.43, 25. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.43, 25. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.43, 25. , 0.2, 0.2,0.1, kGray, 3);
        
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

        DrawGammaLines(0.43, 25. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.43, 25. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.43, 25. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.43, 25. , 0.2, 0.2,0.1, kGray, 3);
        
    canvasWeights->SaveAs(Form("%s/EtaToPi0_Weights.%s",outputDir.Data(),suffix.Data()));


    // *********************************************************************************************************************
    // ************************************ Visualize relative errors EtaToPi0 ******************************************************
    // *********************************************************************************************************************
    
    canvasRelSysErr->cd();
   
    TH2F * histo2DRelSysErrEtaToPi0;
    histo2DRelSysErrEtaToPi0                 = new TH2F("histo2DRelSysErrEtaToPi0","histo2DRelSysErrEtaToPi0",11000,0.43, 25.,1000,0,60.5);
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
    histo2DRelStatErrEtaToPi0                = new TH2F("histo2DRelStatErrEtaToPi0","histo2DRelStatErrEtaToPi0",11000,0.43, 25.,1000,0,60.5);
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
    histo2DRelTotErrEtaToPi0                 = new TH2F("histo2DRelTotErrEtaToPi0","histo2DRelTotErrEtaToPi0",11000,0.43, 25.,1000,0,60.5);
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

        TLatex *labelMassPerf       = new TLatex(0.13,0.875,"ALICE performance");
        SetStyleTLatex( labelMassPerf, textSizeLabelsPixel,4, 1, 43);
        labelMassPerf->Draw();        
        TLatex *labelMassEnergy     = new TLatex(0.13,0.775,collisionSystempPb.Data());
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
        histo2DAllPi0Mass->GetYaxis()->SetRangeUser(126.5,150.8);
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

        
        TH2F * histo2DAllEtaFWHM    = new TH2F("histo2DAllEtaFWHM","histo2DAllEtaFWHM", 20, 0.43, 25. ,1000., -5, 92);
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

        TH2F * histo2DAllEtaMass    = new TH2F("histo2DAllEtaMass","histo2DAllEtaMass",20, 0.43, 25., 1000., 505., 570);
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

        DrawGammaLines(0.43, 25. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,0.3, kGray);
        
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
    histo2DYieldPi0          = new TH2F("histo2DYieldPi0","histo2DYieldPi0",11000,0.23, 25.,1000,2e-9,10e1);
    SetStyleHistoTH2ForGraphs(histo2DYieldPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.65);
    histo2DYieldPi0->GetXaxis()->SetMoreLogLabels();
    histo2DYieldPi0->GetXaxis()->SetLabelOffset(-0.01);
    histo2DYieldPi0->Draw("copy");

        for (Int_t i = 10; i> -1; i--){
            if (graphIndPi0InvYieldSys[i]){
                DrawGammaSetMarkerTGraphAsym(graphIndPi0InvYieldSys[i], markerStyleDet[i], markerSizeDet[i]*0.75, colorDet[i] , colorDet[i], widthLinesBoxes, kTRUE);
                graphIndPi0InvYieldSys[i]->Draw("E2same");
            }    
        }    
        for (Int_t i = 10; i> -1; i--){
            if (graphIndPi0InvYieldStat[i]){
                DrawGammaSetMarkerTGraphAsym(graphIndPi0InvYieldStat[i], markerStyleDet[i], markerSizeDet[i]*0.75, colorDet[i] , colorDet[i]);
                graphIndPi0InvYieldStat[i]->Draw("p,same,z");
            }    
        }    
        TLatex *labelEnergyXSectionPi0      = new TLatex(0.94,0.92,collisionSystempPb.Data());
        SetStyleTLatex( labelEnergyXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);
        labelEnergyXSectionPi0->Draw();
        TLatex *labelDetSysXSectionPi0      = new TLatex(0.94,0.88,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPi0, 0.035,4, 1, 42, kTRUE, 31);
        labelDetSysXSectionPi0->Draw();

                                                
        TLegend* legendXSectionPi0          = GetAndSetLegend2(0.20, 0.13,0.40,0.13+(0.035*7), 0.035, 1, "", 42, 0);
        for (Int_t i = 0; i < 11; i++){
            if (graphIndPi0InvYieldSys[i]) legendXSectionPi0->AddEntry(graphIndPi0InvYieldSys[i],nameMeasGlobalLabel[i],"fp");
        }    
        legendXSectionPi0->Draw();
   
    canvasXSectionPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems.%s",outputDir.Data(),suffix.Data()));
    
    canvasXSectionPi0->cd();
    histo2DYieldPi0->Draw("copy");
        for (Int_t i = 10; i> -1; i--){
            if (graphIndPi0InvYieldSys[i]){
                graphIndPi0InvYieldSys[i]->Draw("E2same");
            }    
        }    
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphCombPi0InvYieldSys->Draw("E2same");
        for (Int_t i = 10; i> -1; i--){
            if (graphIndPi0InvYieldSys[i]){
                graphIndPi0InvYieldStat[i]->Draw("p,same,z");
            }    
        }    
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStat, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphCombPi0InvYieldStat->Draw("p,same,z");

        labelEnergyXSectionPi0->Draw();
        labelDetSysXSectionPi0->Draw();

        legendXSectionPi0->AddEntry(graphCombPi0InvYieldSys,"comb","fp");
        legendXSectionPi0->Draw();
        
    canvasXSectionPi0->SaveAs(Form("%s/Pi0_InvYieldCompAllSystems_Comb.%s",outputDir.Data(),suffix.Data()));
  
    TCanvas* canvasXSectionEta      = new TCanvas("canvasXSectionEta","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasXSectionEta, 0.16, 0.02, 0.02, 0.08);
    canvasXSectionEta->SetLogx();
    canvasXSectionEta->SetLogy();
    
    TH2F * histo2DInvYieldEta;
    histo2DInvYieldEta              = new TH2F("histo2DInvYieldEta","histo2DInvYieldEta",11000,0.43, 25.,1000,5e-10,1e-0);
    SetStyleHistoTH2ForGraphs(histo2DInvYieldEta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",0.035,0.04, 0.035,0.04, 0.9,1.65);
    histo2DInvYieldEta->GetXaxis()->SetMoreLogLabels();
    histo2DInvYieldEta->GetXaxis()->SetLabelOffset(-0.01);
    histo2DInvYieldEta->Draw("copy");

        for (Int_t i = 10; i> -1; i--){
            if (graphIndEtaInvYieldSys[i]){
                DrawGammaSetMarkerTGraphAsym(graphIndEtaInvYieldSys[i], markerStyleDet[i], markerSizeDet[i]*0.75, colorDet[i] , colorDet[i], widthLinesBoxes, kTRUE);
                graphIndEtaInvYieldSys[i]->Draw("E2same");
            }    
        }    
        for (Int_t i = 10; i> -1; i--){
            if (graphIndEtaInvYieldStat[i]){
                DrawGammaSetMarkerTGraphAsym(graphIndEtaInvYieldStat[i], markerStyleDet[i], markerSizeDet[i]*0.75, colorDet[i] , colorDet[i]);
                graphIndEtaInvYieldStat[i]->Draw("p,same,z");
            }    
        }

        TLatex *labelEnergyXSectionEta      = new TLatex(0.94,0.92,collisionSystempPb.Data());
        SetStyleTLatex( labelEnergyXSectionEta, 0.035, 4, 1, 42, kTRUE, 31);
        labelEnergyXSectionEta->Draw();
        TLatex *labelDetSysXSectionEta      = new TLatex(0.94,0.88,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionEta, 0.035, 4, 1, 42, kTRUE, 31);
        labelDetSysXSectionEta->Draw();

        TLegend* legendXSectionEta          = GetAndSetLegend2(0.20, 0.13,0.40,0.13+(0.035*4), 0.035, 1, "", 42, 0);
        for (Int_t i = 0; i < 11; i++){
            if (graphIndEtaInvYieldSys[i]) legendXSectionEta->AddEntry(graphIndEtaInvYieldSys[i],nameMeasGlobalLabel[i],"fp");
        }
        legendXSectionEta->Draw();
   
    canvasXSectionEta->SaveAs(Form("%s/Eta_InvYieldCompAllSystems.%s",outputDir.Data(),suffix.Data()));
    histo2DInvYieldEta->Draw("copy");

        for (Int_t i = 10; i> -1; i--){
            if (graphIndEtaInvYieldSys[i]){
                graphIndEtaInvYieldSys[i]->Draw("E2same");
            }    
        }    
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphCombEtaInvYieldSys->Draw("E2same");
        for (Int_t i = 10; i> -1; i--){
            if (graphIndEtaInvYieldStat[i]){
                graphIndEtaInvYieldStat[i]->Draw("p,same,z");
            }    
        }
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

    
    // **********************************************************************************************************************
    // ******************************************* Comparison to theory calculations Pi0 ************************************
    // **********************************************************************************************************************    
    textSizeLabelsPixel                     = 48;
    
    TCanvas* canvasRatioPPb                 = new TCanvas("canvasRatioPPb","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioPPb,  0.12, 0.01, 0.01, 0.11);
    canvasRatioPPb->cd();
    canvasRatioPPb->SetLogx();

        textsizeLabelsPPb                    = 0;
        textsizeFacPPb                       = 0;
        if (canvasRatioPPb->XtoPixel(canvasRatioPPb->GetX2()) <canvasRatioPPb->YtoPixel(canvasRatioPPb->GetY1()) ){
            textsizeLabelsPPb                = (Double_t)textSizeLabelsPixel/canvasRatioPPb->XtoPixel(canvasRatioPPb->GetX2()) ;
            textsizeFacPPb                   = (Double_t)1./canvasRatioPPb->XtoPixel(canvasRatioPPb->GetX2()) ;
        } else {
            textsizeLabelsPPb                = (Double_t)textSizeLabelsPixel/canvasRatioPPb->YtoPixel(canvasRatioPPb->GetY1());
            textsizeFacPPb                   = (Double_t)1./canvasRatioPPb->YtoPixel(canvasRatioPPb->GetY1());
        }
        cout << textsizeLabelsPPb << endl;

        TH2F * ratio2DTheoryPPb       = new TH2F("ratio2DTheoryPPb","ratio2DTheoryPPb",1000,0.23, 25.,1000,0.2,2.95);
        SetStyleHistoTH2ForGraphs(ratio2DTheoryPPb, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsPPb, textsizeLabelsPPb, 
                                  0.85*textsizeLabelsPPb,textsizeLabelsPPb, 0.9, 0.95, 510, 505);
        ratio2DTheoryPPb->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio2DTheoryPPb->GetYaxis()->SetNdivisions(510);
        ratio2DTheoryPPb->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DTheoryPPb->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DTheoryPPb->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DTheoryPPb->GetXaxis()->SetLabelFont(42);
        ratio2DTheoryPPb->GetYaxis()->SetLabelFont(42);
        ratio2DTheoryPPb->DrawCopy(); 

        TGraphAsymmErrors* graphRatioPi0CombCombFitStatWOXErr = (TGraphAsymmErrors*)graphRatioPi0CombCombFitStat->Clone("graphRatioPi0CombCombFitStatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRatioPi0CombCombFitStatWOXErr);
        
        
//         graphRatioPi0CombNLODSS14->RemovePoint(0);        
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS14, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
//         graphRatioPi0CombNLODSS14->Draw("3,same");
        SetStyleHisto(histoRatioPi0DPMJetToFit, widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );  
        SetStyleHisto(histoRatioPi0HIJINGToFit, widthCommonFit*1.5, styleLineHIJING, colorHIJING );  
        SetStyleHisto(histoRatioPi0EPOS3ToFit, widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );  
        DrawGammaSetMarkerTGraphErr(graphRatioPi0EPOS3ToFit, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
        DrawGammaSetMarkerTGraphErr(graphRatioPi0McGillToFit, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
        
        histoRatioPi0EPOS3ToFit->Draw("same,hist,l");
        graphRatioPi0EPOS3ToFit->Draw("same,3");
        graphRatioPi0McGillToFit->Draw("same,3");
        histoRatioPi0DPMJetToFit->Draw("same,hist,l");
        histoRatioPi0HIJINGToFit->Draw("same,hist,l");
        
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSys->SetLineWidth(0);
        graphRatioPi0CombCombFitSys->Draw("2,same");
        graphRatioPi0CombCombFitStatWOXErr->Draw("p,same");

        TBox* boxErrorNormRatioPi0 = CreateBoxConv(kGray+1, 0.28, 1.-(0.031 ), 0.32, 1.+(0.031));
        boxErrorNormRatioPi0->Draw();
        DrawGammaLines(0.23, 25.,1., 1.,0.1,kGray);

        TLegend* legendRatioTheoryPPbNew= GetAndSetLegend2(0.69,0.91-0.85*textsizeLabelsPPb*5,0.94,0.91, 0.85*textSizeLabelsPixel);
        legendRatioTheoryPPbNew->AddEntry(graphRatioPi0CombCombFitSys,"#pi^{0} ALICE","pf");
//         legendRatioTheoryPPbNew->AddEntry(graphRatioPi0CombNLODSS14,"NLO, PDF: MSTW FF: DSS14","f");
//         legendRatioTheoryPPbNew->AddEntry((TObject*)0,"Phys.Rev. D91 no. 1, (2015) 014035","");
        legendRatioTheoryPPbNew->AddEntry(histoRatioPi0DPMJetToFit,  "DPMJet", "l");  
        legendRatioTheoryPPbNew->AddEntry(histoRatioPi0HIJINGToFit,  "HIJING", "l");  
        legendRatioTheoryPPbNew->AddEntry(histoRatioPi0EPOS3ToFit,  "EPOS", "l");  
        legendRatioTheoryPPbNew->AddEntry(graphRatioPi0McGillToFit,  "Shen #it{et al.}", "f");  
        legendRatioTheoryPPbNew->Draw();

        TLatex *labelRatioTheory   = new TLatex(0.97,0.925,collisionSystempPb.Data());
        SetStyleTLatex( labelRatioTheory, 0.85*textsizeLabelsPPb,4, 1, 42, kTRUE, 31);
        labelRatioTheory->Draw();
        
        ratio2DTheoryPPb->Draw("same,axis");
    canvasRatioPPb->Update();
    canvasRatioPPb->Print(Form("%s/Pi0_RatioTheoryToData_PPb.%s",outputDir.Data(),suffix.Data()));

    // **********************************************************************************************************************
    // ******************************************* Comparison to theory calculations Eta ************************************
    // **********************************************************************************************************************    
    canvasRatioPPb->cd();
    canvasRatioPPb->SetLogx();

        TH2F * ratio2DTheoryPPbEta       = new TH2F("ratio2DTheoryPPbEta","ratio2DTheoryPPbEta",1000,0.43, 25.,1000,0.3,3.2);
        SetStyleHistoTH2ForGraphs(ratio2DTheoryPPbEta, "#it{p}_{T} (GeV/#it{c})","#frac{Theory, Data}{fit}", 0.85*textsizeLabelsPPb, textsizeLabelsPPb, 
                                  0.85*textsizeLabelsPPb,textsizeLabelsPPb, 0.9, 0.95, 510, 505);
        ratio2DTheoryPPbEta->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio2DTheoryPPbEta->GetYaxis()->SetNdivisions(505);
        ratio2DTheoryPPbEta->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DTheoryPPbEta->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DTheoryPPbEta->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DTheoryPPbEta->GetXaxis()->SetLabelFont(42);
        ratio2DTheoryPPbEta->GetYaxis()->SetLabelFont(42);
        ratio2DTheoryPPbEta->DrawCopy(); 

//         DrawGammaNLOTGraph( graphRatioEtaCombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorEPOS3);
//         graphRatioEtaCombNLOMuHalf->Draw("same,c");
//         DrawGammaNLOTGraph( graphRatioEtaCombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorEPOS3);
//         graphRatioEtaCombNLOMuOne->Draw("same,c");
//         DrawGammaNLOTGraph( graphRatioEtaCombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorEPOS3);
//         graphRatioEtaCombNLOMuTwo->Draw("same,c");
        
        SetStyleHisto(histoRatioEtaDPMJetToFit, widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );  
        SetStyleHisto(histoRatioEtaHIJINGToFit, widthCommonFit*1.5, styleLineHIJING, colorHIJING );  
        SetStyleHisto(histoRatioEtaEPOS3ToFit, widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );  
        DrawGammaSetMarkerTGraphErr(graphRatioEtaEPOS3ToFit, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
        DrawGammaSetMarkerTGraphErr(graphRatioEtaMcGillToFit, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
        
        histoRatioEtaEPOS3ToFit->Draw("same,hist,l");
        graphRatioEtaEPOS3ToFit->Draw("same,3");
        graphRatioEtaMcGillToFit->Draw("same,3");
        histoRatioEtaDPMJetToFit->Draw("same,hist,l");
        histoRatioEtaHIJINGToFit->Draw("same,hist,l");
        
        TGraphAsymmErrors* graphRatioEtaCombCombFitStatWOXErr = (TGraphAsymmErrors*)graphRatioEtaCombCombFitStat->Clone("graphRatioEtaCombCombFitStatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphRatioEtaCombCombFitStatWOXErr);
        
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioEtaCombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioEtaCombCombFitSys->SetLineWidth(0);
        graphRatioEtaCombCombFitSys->Draw("2,same");
        graphRatioEtaCombCombFitStatWOXErr->Draw("p,same");

        TBox* boxErrorNormRatioEta = CreateBoxConv(kGray+1, 0.55, 1.-(0.031), 0.6, 1.+(0.031));        
        boxErrorNormRatioEta->Draw();
        DrawGammaLines(0.43, 25.,1., 1.,0.1,kGray);

        // upper left corner
//         TLatex *labelRatioTheoryNLOEta   = new TLatex(0.155,0.93,"NLO, PDF: CTEQ6M5 FF: AESSS");
//         SetStyleTLatex( labelRatioTheoryNLOEta, 0.85*textsizeLabelsPPb,4);
//         labelRatioTheoryNLOEta->Draw();
//         legendRatioPi0OldTheo->Draw();
        
        // lower left corner
        TLegend* legendRatioEtaOldTheo3= GetAndSetLegend2(0.35,0.91-4*0.85*textsizeLabelsPPb,0.6,0.91, 0.85* textSizeLabelsPixel, 1, "", 43, 0.27);
        legendRatioEtaOldTheo3->AddEntry(histoRatioPi0DPMJetToFit,  "DPMJet", "l");  
        legendRatioEtaOldTheo3->AddEntry(histoRatioPi0HIJINGToFit,  "HIJING", "l");  
        legendRatioEtaOldTheo3->AddEntry(histoRatioPi0EPOS3ToFit,  "EPOS3", "l");  
        legendRatioEtaOldTheo3->AddEntry(graphRatioEtaMcGillToFit,  "Shen #it{et al.}", "f");  
        legendRatioEtaOldTheo3->Draw();

        // upper right corner
        TLegend* legendRatioEtaOldTheo2= GetAndSetLegend2(0.15,0.91-1*0.85*textsizeLabelsPPb,0.4,0.91, 0.85* textSizeLabelsPixel, 1, "", 43, 0.27);
        legendRatioEtaOldTheo2->AddEntry(graphRatioPi0CombCombFitSys,"#eta ALICE","pf");
        legendRatioEtaOldTheo2->Draw();
        TLatex *labelRatioTheory2   = new TLatex(0.16,0.925,collisionSystempPb.Data());
        SetStyleTLatex( labelRatioTheory2, 0.85*textsizeLabelsPPb,4);
        labelRatioTheory2->Draw();
            
        ratio2DTheoryPPbEta->Draw("same,axis");
    canvasRatioPPb->Update();
    canvasRatioPPb->Print(Form("%s/Eta_RatioTheoryToData_PPb.%s",outputDir.Data(),suffix.Data()));
    
    //*************************************************************************************************************
    //***************************** Paper plot X-section and ratios ***********************************************
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

    
    padInvSectionSpec->cd();
    padInvSectionSpec->SetLogy(1);
    padInvSectionSpec->SetLogx(1);
        SetStyleHistoTH2ForGraphs(histo2DYieldPi0, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                                0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.26/(textsizeFacXSecUp*marginXSec));
        histo2DYieldPi0->GetXaxis()->SetMoreLogLabels();
        histo2DYieldPi0->GetXaxis()->SetLabelOffset(+0.01);
        histo2DYieldPi0->Draw();
        
        SetStyleHisto(histoEPOS3Pi0, widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );  
        DrawGammaSetMarkerTGraphErr(graphEPOS3Pi0, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
        histoEPOS3Pi0->GetXaxis()->SetRangeUser(0.3,20);
        histoEPOS3Pi0->Draw("same,hist,l");
        while(graphEPOS3Pi0->GetX()[0] < 0.3)
            graphEPOS3Pi0->RemovePoint(0);
        while(graphEPOS3Pi0->GetX()[graphEPOS3Pi0->GetN()-1] > 20)
            graphEPOS3Pi0->RemovePoint(graphEPOS3Pi0->GetN()-1);
        graphEPOS3Pi0->Draw("same,3");
        
        DrawGammaSetMarkerTGraphErr(graphMcGillPi0, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
        while(graphMcGillPi0->GetX()[0] < 0.3)
            graphMcGillPi0->RemovePoint(0);
        graphMcGillPi0->Draw("same,3");
        
        SetStyleHisto(histoDPMJetPi0, widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );  
        histoDPMJetPi0->GetXaxis()->SetRangeUser(0.3,20);
        histoDPMJetPi0->Draw("same,hist,l");
        SetStyleHisto(histoHIJINGPi0, widthCommonFit*1.5, styleLineHIJING, colorHIJING );  
        histoHIJINGPi0->GetXaxis()->SetRangeUser(0.3,20);
        histoHIJINGPi0->Draw("same,hist,l");

//         graphNLODSS14Calc->RemovePoint(0);
//         DrawGammaSetMarkerTGraphAsym(graphNLODSS14Calc, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
//         graphNLODSS14Calc->Draw("3,same");

        TGraphAsymmErrors* graphCombPi0InvYieldStatWOXErr = (TGraphAsymmErrors*)graphCombPi0InvYieldStat->Clone("graphCombPi0InvYieldStatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphCombPi0InvYieldStatWOXErr);
        
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombPi0InvYieldSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvYieldStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
        graphCombPi0InvYieldStatWOXErr->Draw("p,same,z");

        DrawGammaSetMarkerTF1( fitTCMInvYieldPi0, 7, 2, kGray+2); 
        fitTCMInvYieldPi0->Draw("same");        

//         DrawGammaSetMarkerTF1( fitPowInvYieldPi0, 7, 2, kBlue-7); 
//         fitPowInvYieldPi0->Draw("same");        
                
        TLatex *labelEnergyXSectionPaper= new TLatex(0.95, 0.91, collisionSystempPb.Data());
        SetStyleTLatex( labelEnergyXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
        labelEnergyXSectionPaper->Draw();
        TLatex *labelALICEXSectionPaper= new TLatex(0.95,0.87,"ALICE");
        SetStyleTLatex( labelALICEXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
        labelALICEXSectionPaper->Draw();
        TLatex *labelDetSysXSectionPaper= new TLatex(0.95,0.83,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPaper, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 31);
        labelDetSysXSectionPaper->Draw();
        
        TLegend* legendXsectionPaper    = GetAndSetLegend2(0.19, 0.14, 0.5, 0.14+0.05*6, textSizeLabelsPixel, 1, "", 43, 0.2);
        legendXsectionPaper->AddEntry(graphCombPi0InvYieldSys,"Data,","pf");
        legendXsectionPaper->AddEntry((TObject*)0,"Norm. unc. 3.1%","");
//         legendXsectionPaper->AddEntry(graphNLODSS14Calc,"NLO, PDF: MSTW, FF: DSS14","f");
//         legendXsectionPaper->AddEntry((TObject*)0,"Phys.Rev. D91 no. 1, (2015) 014035","");
        legendXsectionPaper->AddEntry(histoDPMJetPi0,"DPMJet","l");
        legendXsectionPaper->AddEntry(histoHIJINGPi0,"HIJING","l");
        legendXsectionPaper->AddEntry(histoEPOS3Pi0,"EPOS3","l");
        legendXsectionPaper->AddEntry(graphMcGillPi0,"Shen #it{et al.}","f");
        legendXsectionPaper->Draw();

        TLegend* legendXsectionPaper2     = GetAndSetLegend2(0.19, 0.05, 0.5, 0.11, textSizeLabelsPixel, 1, "", 43, 0.2);
        legendXsectionPaper2->AddEntry(fitTCMInvYieldPi0,"#it{A}_{e} exp(-#it{E}_{T, kin}/#it{T}_{e}) + #it{A}/#(){1 + #frac{#it{p}_{T}^{2}}{#it{T}^{2}#upoint n}}^{n}","l");
        legendXsectionPaper2->Draw();
        
        
    padInvSectionNLORatio->cd();
    padInvSectionNLORatio->SetLogx(1);
        TH2F * ratio2DEPOS               = new TH2F("ratio2DEPOS","ratio2DEPOS",1000,0.23, 25.,1000,0.3,2.2);
        SetStyleHistoTH2ForGraphs(ratio2DEPOS, "#it{p}_{T} (GeV/#it{c})","#frac{EPOS, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle, 
                                  0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.26/(textsizeFacXSecMiddle*marginXSec), 510, 508);
//         ratio2DEPOS->GetYaxis()->SetNdivisions(512);
        ratio2DEPOS->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DEPOS->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DEPOS->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DEPOS->GetXaxis()->SetLabelFont(42);
        ratio2DEPOS->GetYaxis()->SetLabelFont(42);
        ratio2DEPOS->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DEPOS->GetXaxis()->SetTickLength(0.07);
        ratio2DEPOS->DrawCopy();

        boxErrorNormRatioPi0->Draw();
        
//         DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombNLODSS14, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
//         graphRatioPi0CombNLODSS14->Draw("3,same");
        SetStyleHisto(histoRatioPi0EPOS3ToFit, widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );  
        DrawGammaSetMarkerTGraphErr(graphRatioPi0EPOS3ToFit, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
        histoRatioPi0EPOS3ToFit->Draw("same,hist,l");
        graphRatioPi0EPOS3ToFit->Draw("same,3");
        DrawGammaSetMarkerTGraphErr(graphRatioPi0McGillToFit, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
        graphRatioPi0McGillToFit->Draw("same,3");

        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSys->SetLineWidth(0);
        graphRatioPi0CombCombFitSys->Draw("2,same");
        graphRatioPi0CombCombFitStatWOXErr->Draw("p,same");
        
        DrawGammaLines(0.23, 25.,1., 1.,0.1,kGray);
        
    padInvSectionPythiaRatio->cd();
    padInvSectionPythiaRatio->SetLogx(1);
        TH2F * ratio2DGenerator            = new TH2F("ratio2DGenerator","ratio2DGenerator",1000,0.23, 25.,1000,0.3,2.2);
        SetStyleHistoTH2ForGraphs(ratio2DGenerator, "#it{p}_{T} (GeV/#it{c})","#frac{Generator, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown, 
                                  0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.26/(textsizeFacXSecDown*marginXSec), 510, 508);
//         ratio2DGenerator->GetYaxis()->SetNdivisions(512);
        ratio2DGenerator->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DGenerator->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DGenerator->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DGenerator->GetXaxis()->SetLabelFont(42);
        ratio2DGenerator->GetYaxis()->SetLabelFont(42);
        ratio2DGenerator->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DGenerator->GetXaxis()->SetTickLength(0.06);
        ratio2DGenerator->GetYaxis()->SetTickLength(0.04);
        ratio2DGenerator->DrawCopy();

        SetStyleHisto(histoRatioPi0DPMJetToFit, widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );  
//         histoRatioPi0DPMJetToFit->GetXaxis()->SetRangeUser(0.5,14);
        histoRatioPi0DPMJetToFit->Draw("same,hist,l");
        SetStyleHisto(histoRatioPi0HIJINGToFit, widthCommonFit*1.5, styleLineHIJING, colorHIJING );  
//         histoRatioPi0HIJINGToFit->GetXaxis()->SetRangeUser(0.5,14);
        histoRatioPi0HIJINGToFit->Draw("same,hist,l");
        
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioPi0CombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioPi0CombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioPi0CombCombFitSys->SetLineWidth(0);
        graphRatioPi0CombCombFitSys->Draw("2,same");
        graphRatioPi0CombCombFitStatWOXErr->Draw("p,same");
        
        boxErrorNormRatioPi0->Draw();
        
        DrawGammaLines(0.23, 25.,1., 1.,0.1,kGray);
        
    canvasInvSectionPaper->Print(Form("%s/Pi0_InvYieldWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));
    
    padInvSectionSpec->cd();
    padInvSectionSpec->SetLogy(1);
    padInvSectionSpec->SetLogx(1);
        SetStyleHistoTH2ForGraphs(histo2DInvYieldEta, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                                0.85*textsizeLabelsXSecUp,textsizeLabelsXSecUp, 0.85*textsizeLabelsXSecUp, textsizeLabelsXSecUp, 1,0.26/(textsizeFacXSecUp*marginXSec));
        histo2DInvYieldEta->GetXaxis()->SetMoreLogLabels();
        histo2DInvYieldEta->GetXaxis()->SetLabelOffset(+0.01);
        histo2DInvYieldEta->Draw();

        SetStyleHisto(histoEPOS3Eta, widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );  
        DrawGammaSetMarkerTGraphErr(graphEPOS3Eta, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
        histoEPOS3Eta->GetXaxis()->SetRangeUser(0.6,20);
        histoEPOS3Eta->Draw("same,hist,l");
        while(graphEPOS3Eta->GetX()[0] < 0.6)
            graphEPOS3Eta->RemovePoint(0);
        while(graphEPOS3Eta->GetX()[graphEPOS3Eta->GetN()-1] > 20)
            graphEPOS3Eta->RemovePoint(graphEPOS3Eta->GetN()-1);
        graphEPOS3Eta->Draw("same,3");

        DrawGammaSetMarkerTGraphErr(graphMcGillEta, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
        while(graphMcGillEta->GetX()[0] < 0.6)
            graphMcGillEta->RemovePoint(0);
        graphMcGillEta->Draw("same,3");
        
        SetStyleHisto(histoDPMJetEta, widthCommonFit, styleLineDPMJet, colorDPMJet);          
        histoDPMJetEta->GetXaxis()->SetRangeUser(0.6,20);
        histoDPMJetEta->Draw("same,hist,l");
        SetStyleHisto(histoHIJINGEta, widthCommonFit, styleLineHIJING, colorHIJING);          
        histoHIJINGEta->GetXaxis()->SetRangeUser(0.6,20);
        histoHIJINGEta->Draw("same,hist,l");
        SetStyleHisto(histoEPOS3Eta, widthCommonFit, styleLineEPOS3, colorEPOS3);          
        histoEPOS3Eta->GetXaxis()->SetRangeUser(0.6,20);
        histoEPOS3Eta->Draw("same,hist,l");

        TGraphAsymmErrors* graphCombEtaInvYieldStatWOXErr = (TGraphAsymmErrors*)graphCombEtaInvYieldStat->Clone("graphCombEtaInvYieldStatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphCombEtaInvYieldStatWOXErr);

        
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombEtaInvYieldSys->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphCombEtaInvYieldStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack);
        graphCombEtaInvYieldStatWOXErr->Draw("p,same,z");

        DrawGammaSetMarkerTF1( fitTCMInvYieldEta, 7, 2, kGray+2); 
        fitTCMInvYieldEta->Draw("same");

        // labels lower left corner
        TLegend* legendXsectionPaperEta    = GetAndSetLegend2(0.19, 0.11, 0.5, 0.11+0.05*2, textSizeLabelsPixel, 1, "", 43, 0.2);
        legendXsectionPaperEta->AddEntry(graphCombPi0InvYieldSys,"Data,","pf");
        legendXsectionPaperEta->AddEntry((TObject*)0,"Norm. unc. 3.1%","");
        legendXsectionPaperEta->Draw();
        legendXsectionPaper2->Draw();
        
        TLatex *labelEnergyXSectionPaperEta= new TLatex(0.20, 0.13+0.05*4, collisionSystempPb.Data());
        SetStyleTLatex( labelEnergyXSectionPaperEta, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 11);
        labelEnergyXSectionPaperEta->Draw();
        TLatex *labelALICEXSectionPaperEta= new TLatex(0.20,0.13+0.05*3,"ALICE");
        SetStyleTLatex( labelALICEXSectionPaperEta, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 11);
        labelALICEXSectionPaperEta->Draw();
        TLatex *labelDetSysXSectionPaperEta = new TLatex(0.20,0.13+0.05*2,"#eta #rightarrow #gamma#gamma");
        SetStyleTLatex( labelDetSysXSectionPaperEta, textsizeLabelsXSecUp,4, 1, 42, kTRUE, 11);
        labelDetSysXSectionPaperEta->Draw();
        

        // labels upper right corner
        TLegend* legendXsectionPaperEtaTheo     = GetAndSetLegend2(0.68, 0.95-0.04*4, 0.75+0.33, 0.95, textSizeLabelsPixel, 1, "", 43, 0.2);
        legendXsectionPaperEtaTheo->AddEntry(histoDPMJetEta,"DPMJet","l");
        legendXsectionPaperEtaTheo->AddEntry(histoHIJINGEta,"HIJING","l");
        legendXsectionPaperEtaTheo->AddEntry(histoEPOS3Eta,"EPOS3","l");
        legendXsectionPaperEtaTheo->AddEntry(graphMcGillEta,"Shen #it{et al.}","f");
        legendXsectionPaperEtaTheo->Draw();
        
    padInvSectionNLORatio->cd();
    padInvSectionNLORatio->SetLogx(1);
        TH2F * ratio2DEPOSEta                = new TH2F("ratio2DEPOSEta","ratio2DEPOSEta",1000,0.43, 25.,1000,0.3,3.25);
        SetStyleHistoTH2ForGraphs(ratio2DEPOSEta, "#it{p}_{T} (GeV/#it{c})","#frac{EPOS, Data}{fit}", 0.85*textsizeLabelsXSecMiddle, textsizeLabelsXSecMiddle, 
                                  0.85*textsizeLabelsXSecMiddle,textsizeLabelsXSecMiddle, 1,0.26/(textsizeFacXSecMiddle*marginXSec), 510, 510);
        ratio2DEPOSEta->GetYaxis()->SetMoreLogLabels(kTRUE);
        ratio2DEPOSEta->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DEPOSEta->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DEPOSEta->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DEPOSEta->GetXaxis()->SetLabelFont(42);
        ratio2DEPOSEta->GetYaxis()->SetLabelFont(42);
        ratio2DEPOSEta->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DEPOSEta->GetXaxis()->SetTickLength(0.07);
        ratio2DEPOSEta->DrawCopy();

//         DrawGammaNLOTGraph( graphRatioEtaCombNLOMuHalf, widthCommonFit, styleLineNLOMuHalf, colorEPOS3);
//         graphRatioEtaCombNLOMuHalf->Draw("same,c");
//         DrawGammaNLOTGraph( graphRatioEtaCombNLOMuOne, widthCommonFit, styleLineNLOMuOne, colorEPOS3);
//         graphRatioEtaCombNLOMuOne->Draw("same,c");
//         DrawGammaNLOTGraph( graphRatioEtaCombNLOMuTwo, widthCommonFit, styleLineNLOMuTwo, colorEPOS3);
//         graphRatioEtaCombNLOMuTwo->Draw("same,c");
        SetStyleHisto(histoRatioEtaEPOS3ToFit, widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );  
        DrawGammaSetMarkerTGraphErr(graphRatioEtaEPOS3ToFit, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
        histoRatioEtaEPOS3ToFit->Draw("same,hist,l");
        graphRatioEtaEPOS3ToFit->Draw("same,3");

        DrawGammaSetMarkerTGraphErr(graphRatioEtaMcGillToFit, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
        graphRatioEtaMcGillToFit->Draw("same,3");
        
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioEtaCombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioEtaCombCombFitSys->SetLineWidth(0);
        graphRatioEtaCombCombFitSys->Draw("2,same");
        graphRatioEtaCombCombFitStatWOXErr->Draw("p,same");
        
        boxErrorNormRatioEta->Draw();
        DrawGammaLines(0.43, 25.,1., 1.,0.1,kGray);
        
    padInvSectionPythiaRatio->cd();
    padInvSectionPythiaRatio->SetLogx(1);
        TH2F * ratio2DGeneratorEta             = new TH2F("ratio2DGeneratorEta","ratio2DGeneratorEta",1000,0.43, 25.,1000,0.3,2.95);
        SetStyleHistoTH2ForGraphs(ratio2DGeneratorEta, "#it{p}_{T} (GeV/#it{c})","#frac{Generator, Data}{fit}", 0.85*textsizeLabelsXSecDown, textsizeLabelsXSecDown, 
                                  0.85*textsizeLabelsXSecDown,textsizeLabelsXSecDown, 0.9,0.26/(textsizeFacXSecDown*marginXSec), 510, 510);
        ratio2DGeneratorEta->GetYaxis()->SetNoExponent(kTRUE);
        ratio2DGeneratorEta->GetXaxis()->SetMoreLogLabels(kTRUE);
        ratio2DGeneratorEta->GetXaxis()->SetNoExponent(kTRUE);
        ratio2DGeneratorEta->GetXaxis()->SetLabelFont(42);
        ratio2DGeneratorEta->GetYaxis()->SetLabelFont(42);
        ratio2DGeneratorEta->GetYaxis()->SetLabelOffset(+0.01);
        ratio2DGeneratorEta->GetXaxis()->SetTickLength(0.06);
        ratio2DGeneratorEta->GetYaxis()->SetTickLength(0.04);
        ratio2DGeneratorEta->DrawCopy();

        SetStyleHisto(histoRatioEtaDPMJetToFit, widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );  
        histoRatioEtaDPMJetToFit->Draw("same,hist,l");
        SetStyleHisto(histoRatioEtaHIJINGToFit, widthCommonFit*1.5, styleLineHIJING, colorHIJING );  
        histoRatioEtaHIJINGToFit->Draw("same,hist,l");
        
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitStatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphRatioEtaCombCombFitStatWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaCombCombFitSys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kTRUE, 0);
        graphRatioEtaCombCombFitSys->SetLineWidth(0);
        graphRatioEtaCombCombFitSys->Draw("2,same");
        graphRatioEtaCombCombFitStatWOXErr->Draw("p,same");

        boxErrorNormRatioEta->Draw();
        DrawGammaLines(0.43, 25.,1., 1.,0.1,kGray);
        
    canvasInvSectionPaper->Print(Form("%s/Eta_InvYieldWithRatios_Paper.%s",outputDir.Data(),suffix.Data()));


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
    histo2DEtatoPi0combo                = new TH2F("histo2DEtatoPi0combo","histo2DEtatoPi0combo",1000,0.43, 25.,1000,0.,1.05    );
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

        TGraphAsymmErrors* graphCombEtaToPi0StatWOXErr = (TGraphAsymmErrors*)graphCombEtaToPi0Stat->Clone("graphCombEtaToPi0StatWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphCombEtaToPi0StatWOXErr);
    
        // plotting data
        graphCombEtaToPi0Stat->Print();
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0StatWOXErr, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes, kFALSE);
        graphCombEtaToPi0StatWOXErr->SetLineWidth(widthLinesBoxes);
        DrawGammaSetMarkerTGraphAsym(graphCombEtaToPi0Sys, markerStyleComb, markerSizeComb, kBlack, kBlack, widthLinesBoxes*2, kTRUE, 0);
        graphCombEtaToPi0Sys->SetLineWidth(0);
        graphCombEtaToPi0Sys->Draw("2,same");
        graphCombEtaToPi0StatWOXErr->Draw("p,same");
        
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
//         DrawGammaSetMarkerTGraphAsym(graphNLOEtaToPi0, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
//         graphNLOEtaToPi0->Draw("3,same");

        DrawGammaSetMarkerTGraphErr(graphMcGillEtaToPi0, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
        graphMcGillEtaToPi0->Draw("3,same");
    
        SetStyleHisto(histoEPOS3EtaToPi0_Reb, widthCommonFit*1.5, 1, colorEPOS3);          
        DrawGammaSetMarkerTGraphErr(graphEPOS3EtaToPi0_Reb, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
        histoEPOS3EtaToPi0_Reb->Draw("same,hist,l");
        graphEPOS3EtaToPi0_Reb->Draw("3,same");

        
        SetStyleHisto(cocktailEtaToPi0_MtScaledReb, widthCommonFit*1.5, 1, kRed+2);  
        cocktailEtaToPi0_MtScaledReb->Draw("same,hist,l");

        SetStyleHisto(histoDPMJetEtaToPi0, widthCommonFit*1.5, styleLineDPMJet, colorDPMJet );  
        histoDPMJetEtaToPi0->Draw("same,hist,l");

        SetStyleHisto(histoHIJINGEtaToPi0, widthCommonFit*1.5, 8, colorHIJING);  
        histoHIJINGEtaToPi0->Draw("same,hist,l");
        
        
        graphCombEtaToPi0Sys->Draw("2,same");
        graphCombEtaToPi0StatWOXErr->Draw("p,same");
        
        // plotting labels
        labelEnergyEtaToPi0->Draw();
        labelALICEEtaToPi0->Draw();
//         labelPi0EtaToPi0->Draw();
        
        TLegend* legendEtaToPi0Theory = GetAndSetLegend2(0.48, 0.15+(textsizeLabelsEtaToPi0*3*0.9), 0.9, 0.15, textSizeLabelsPixel*0.85, 2, "", 43, 0.25);
//         legendEtaToPi0Theory->AddEntry(graphNLOEtaToPi0,"NLO, PDF:CTEQ6M5 ","pf");
//         legendEtaToPi0Theory->AddEntry((TObject*)0,"#pi^{0} FF: DSS07, #eta FF: AESSS","");
        legendEtaToPi0Theory->AddEntry(cocktailEtaToPi0_MtScaledReb,"#eta from #it{m}_{_{T}} scaled #pi^{0}","l");
        legendEtaToPi0Theory->AddEntry((TObject*)0,"","");
        legendEtaToPi0Theory->AddEntry(histoDPMJetEtaToPi0,"DPMJet","l");
        legendEtaToPi0Theory->AddEntry(histoEPOS3EtaToPi0_Reb,"EPOS3","l");
        legendEtaToPi0Theory->AddEntry(histoHIJINGEtaToPi0,"HIJING","l");
        legendEtaToPi0Theory->AddEntry(graphMcGillEtaToPi0,"Shen #it{et al.}","f");
        legendEtaToPi0Theory->Draw();
        
    histo2DEtatoPi0combo->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_Theory_Paper.%s",outputDir.Data(), suffix.Data()));

    // ***************************************************************************************************************
    // **************** Plotting eta/pi0 ratio for combined measurement + Theory + other energies ********************
    // ***************************************************************************************************************    
    canvasEtatoPi0combo->cd();
    TH2F * histo2DEtatoPi0WD           = new TH2F("histo2DEtatoPi0WD","histo2DEtatoPi0WD",1000,0.04,25.,1000,-0.18,1.02    );
    SetStyleHistoTH2ForGraphs(histo2DEtatoPi0WD, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.85*textsizeLabelsEtaToPi0, textsizeLabelsEtaToPi0, 
                              0.85*textsizeLabelsEtaToPi0,1.1*textsizeLabelsEtaToPi0, 0.9, 0.65, 510, 510);
    histo2DEtatoPi0WD->GetXaxis()->SetMoreLogLabels();
    histo2DEtatoPi0WD->GetXaxis()->SetLabelOffset(-0.01);
    histo2DEtatoPi0WD->Draw("copy");

        graphCombEtaToPi0Sys->Draw("2,same");
               
        DrawGammaSetMarkerTGraphErr(graphPHENIXEtaToPi0200GeV, 25, 2., kGray+1, kGray+1, widthLinesBoxes, kFALSE);
        graphPHENIXEtaToPi0200GeV->Draw("p,same");
        DrawGammaSetMarkerTGraphAsym(graphALICEEtaToPi07TeV, 24, 2.2, kBlue-6, kBlue-6, widthLinesBoxes, kFALSE);
        graphALICEEtaToPi07TeV->Draw("p,same");
        
        TGraphAsymmErrors* graphALICEEtaToPi02760GeVWOXErr = (TGraphAsymmErrors*)graphALICEEtaToPi02760GeV->Clone("graphALICEEtaToPi02760GeVWOXErr");
        ProduceGraphAsymmWithoutXErrors(graphALICEEtaToPi02760GeVWOXErr);
        DrawGammaSetMarkerTGraphAsym(graphALICEEtaToPi02760GeVWOXErr, 25, 2.2, kBlue-6, kBlue-6, widthLinesBoxes, kFALSE);
        graphALICEEtaToPi02760GeVWOXErr->Draw("p,same");
        
        DrawGammaSetMarkerTGraphErr(graphCERESEtaToPi0pAu29100MeV, 24, 2.5, kGreen-6, kGreen-6, widthLinesBoxes, kFALSE);
        graphCERESEtaToPi0pAu29100MeV->Draw("p,same");
        DrawGammaSetMarkerTGraphErr(graphCERESEtaToPi0pBe29100MeV, 25, 2., kGreen-6, kGreen-6, widthLinesBoxes, kFALSE);
        graphCERESEtaToPi0pBe29100MeV->Draw("p,same");
        DrawGammaSetMarkerTGraphErr(graphPHENIXEtaToPi0dAu200GeV, 24, 2., kGray+1, kGray+1, widthLinesBoxes, kFALSE);
        graphPHENIXEtaToPi0dAu200GeV->Draw("p,same");
        
        graphCombEtaToPi0StatWOXErr->Draw("p,same");
        
        TLegend* legendEtaToPi0Theory2 = GetAndSetLegend2(0.13, 0.95, 0.46, 0.95-(textsizeLabelsEtaToPi0*4*0.9), textSizeLabelsPixel*0.85, 1, "", 43, 0.16);
        legendEtaToPi0Theory2->AddEntry(graphCombEtaToPi0Sys,Form("ALICE, %s",collisionSystempPb.Data()),"pf");
        legendEtaToPi0Theory2->AddEntry(graphPHENIXEtaToPi0dAu200GeV,"PHENIX, d-Au, #sqrt{#it{s}_{_{NN}}} = 0.2 TeV","p");
        legendEtaToPi0Theory2->AddEntry(graphCERESEtaToPi0pAu29100MeV,"TAPS/CERES, p-Au, #sqrt{#it{s}_{_{NN}}} = 29.1 GeV","p");
        legendEtaToPi0Theory2->AddEntry(graphCERESEtaToPi0pBe29100MeV,"TAPS/CERES, p-Be, #sqrt{#it{s}_{_{NN}}} = 29.1 GeV","p");
        legendEtaToPi0Theory2->Draw();

        TLegend* legendEtaToPi0WorldData = GetAndSetLegend2(0.53, 0.145+(textsizeLabelsEtaToPi0*3*0.9), 0.9, 0.145, textSizeLabelsPixel*0.85, 1, "", 43, 0.16);
        legendEtaToPi0WorldData->AddEntry(graphALICEEtaToPi07TeV,"ALICE, pp, #sqrt{#it{s}} = 7 TeV","p");
        legendEtaToPi0WorldData->AddEntry(graphALICEEtaToPi02760GeVWOXErr,"ALICE, pp, #sqrt{#it{s}} = 2.76 TeV","p");
        legendEtaToPi0WorldData->AddEntry(graphPHENIXEtaToPi0200GeV,"PHENIX, pp, #sqrt{#it{s}} = 0.2 TeV", "p");
        legendEtaToPi0WorldData->Draw();
        
    histo2DEtatoPi0WD->Draw("axis,same");

    canvasEtatoPi0combo->Update();
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0_WorldData_Paper.%s",outputDir.Data(), suffix.Data()));

    // ***************************************************************************************************************
    // ******************************** fitting eta/pi0 **************************************************************
    // ***************************************************************************************************************    
    TF1* etaToPi0ConstData  = new TF1("etaToPi0ConstData","[0]",6,20);
    TF1* etaToPi0ConstMC    = new TF1("etaToPi0ConstMC","[0]",8,20);
    TF1* etaToPi0ConstMC2   = new TF1("etaToPi0ConstMC","[0]",8,20);
    graphCombEtaToPi0StatWOXErr->Fit(etaToPi0ConstData,"QRME0","",6,20);
    histoDPMJetEtaToPi0->Fit(etaToPi0ConstMC,"QRME0","",8,20);
    histoHIJINGEtaToPi0->Fit(etaToPi0ConstMC2,"QRME0","",8,20);
    
    cout << "***********************************************************************************************************" << endl;
    cout << "***********************************************************************************************************" << endl;
    cout << "***********************************************************************************************************" << endl;
    cout << "high pt eta/pi0 - data, stat: " << etaToPi0ConstData->GetParameter(0) << "+-"<< etaToPi0ConstData->GetParError(0) << endl;
    graphCombEtaToPi0Tot->Fit(etaToPi0ConstData,"QRME0","",6,20);
    cout << "high pt eta/pi0 - data, tot: " << etaToPi0ConstData->GetParameter(0) << "+-"<< etaToPi0ConstData->GetParError(0) << endl;
    cout << "high pt eta/pi0 - DPMJet: " << etaToPi0ConstMC->GetParameter(0) << "+-"<< etaToPi0ConstMC->GetParError(0) << endl;
    cout << "high pt eta/pi0 - HIJING: " << etaToPi0ConstMC2->GetParameter(0) << "+-"<< etaToPi0ConstMC2->GetParError(0) << endl;
    cout << "***********************************************************************************************************" << endl;
    cout << "***********************************************************************************************************" << endl;
    cout << "***********************************************************************************************************" << endl;

    fileFitsOutput << "***********************************************************************************************************" << endl;
    fileFitsOutput << "***********************************************************************************************************" << endl;
    fileFitsOutput << "***********************************************************************************************************" << endl;
    graphCombEtaToPi0StatWOXErr->Fit(etaToPi0ConstData,"QRME0","",6,20);
    fileFitsOutput << "high pt eta/pi0 - data, stat: " << etaToPi0ConstData->GetParameter(0) << "+-"<< etaToPi0ConstData->GetParError(0) << endl;
    graphCombEtaToPi0Tot->Fit(etaToPi0ConstData,"QRME0","",6,20);
    fileFitsOutput << "high pt eta/pi0 - data, tot: " << etaToPi0ConstData->GetParameter(0) << "+-"<< etaToPi0ConstData->GetParError(0) << endl;
    fileFitsOutput << "high pt eta/pi0 - DPMJet: " << etaToPi0ConstMC->GetParameter(0) << "+-"<< etaToPi0ConstMC->GetParError(0) << endl;
    fileFitsOutput << "high pt eta/pi0 - HIJING: " << etaToPi0ConstMC2->GetParameter(0) << "+-"<< etaToPi0ConstMC2->GetParError(0) << endl;
    fileFitsOutput << "***********************************************************************************************************" << endl;
    fileFitsOutput << "***********************************************************************************************************" << endl;
    fileFitsOutput << "***********************************************************************************************************" << endl;

    canvasXSectionPi0->cd();
    TH2F * histo2DXSectionWithEtaAndPi0;
    histo2DXSectionWithEtaAndPi0          = new TH2F("histo2DXSectionWithEtaAndPi0","histo2DXSectionWithEtaAndPi0",11000,0.23, 25.,1000,2e-12,10e1);
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
        TH1D* histFitTCMInvYieldEta                             = (TH1D*)fitTCMInvYieldEta->GetHistogram();
        histFitTCMInvYieldEta->Scale(scaleFacEtaForCombPlot);
        histoDPMJetEta->Scale(scaleFacEtaForCombPlot);
        histoHIJINGEta->Scale(scaleFacEtaForCombPlot);
        histoEPOS3Eta->Scale(scaleFacEtaForCombPlot);
        TGraphErrors* graphEPOS3EtaCopy                         = (TGraphErrors*) graphEPOS3Eta->Clone("graphEPOS3EtaCopy");
        graphEPOS3EtaCopy                                       = ScaleGraph(graphEPOS3EtaCopy,scaleFacEtaForCombPlot);
        TGraphErrors* graphMcGillEtaCopy                        = (TGraphErrors*) graphMcGillEta->Clone("graphMcGillEtaCopy");
        graphMcGillEtaCopy                                      = ScaleGraph(graphMcGillEtaCopy,scaleFacEtaForCombPlot);
//         TGraphAsymmErrors* graphNLOEtaAESSCalcCopy                  = (TGraphAsymmErrors*)graphNLOEtaAESSCalc->Clone("graphNLOEtaAESSCalcCopy");
//         TGraph* graphNLOCalcEtaMuHalfCopy                           = (TGraph*)graphNLOCalcEtaMuHalf->Clone("graphNLOCalcEtaMuHalfCopy");
//         TGraph* graphNLOCalcEtaMuOneCopy                            = (TGraph*)graphNLOCalcEtaMuOne->Clone("graphNLOCalcEtaMuOneCopy");
//         TGraph* graphNLOCalcEtaMuTwoCopy                            = (TGraph*)graphNLOCalcEtaMuTwo->Clone("graphNLOCalcEtaMuTwoCopy");
//         graphNLOEtaAESSCalcCopy                                     = ScaleGraph(graphNLOEtaAESSCalcCopy,scaleFacEtaForCombPlot);
//         graphNLOCalcEtaMuHalfCopy                                   = ScaleGraph(graphNLOCalcEtaMuHalfCopy,scaleFacEtaForCombPlot);
//         graphNLOCalcEtaMuOneCopy                                    = ScaleGraph(graphNLOCalcEtaMuOneCopy,scaleFacEtaForCombPlot);
//         graphNLOCalcEtaMuTwoCopy                                    = ScaleGraph(graphNLOCalcEtaMuTwoCopy,scaleFacEtaForCombPlot);
        
        // plotting Generators
        histoEPOS3Pi0->Draw("same,hist,l");
        graphEPOS3Pi0->Draw("same,3");
        graphMcGillPi0->Draw("same,3");
        SetStyleHisto(histoEPOS3Eta, widthCommonFit*1.5, styleLineEPOS3, colorEPOS3 );  

        DrawGammaSetMarkerTGraphErr(graphEPOS3EtaCopy, 0, 0, colorEPOS3, colorEPOS3, widthLinesBoxes, kTRUE, colorEPOS3);
        histoEPOS3Eta->Draw("same,hist,l");
        graphEPOS3EtaCopy->Draw("same,3");        
        DrawGammaSetMarkerTGraphErr(graphMcGillEtaCopy, 0, 0, colorMcGill, colorMcGill, widthLinesBoxes, kTRUE, colorMcGill);
        graphMcGillEtaCopy->Draw("same,3");        
        SetStyleHisto(histoDPMJetPi0, widthCommonFit, styleLineDPMJet, colorDPMJet);          
        histoDPMJetPi0->Draw("same,hist,l");
        SetStyleHisto(histoDPMJetEta, widthCommonFit, styleLineDPMJet, colorDPMJet);          
        histoDPMJetEta->Draw("same,hist,l");
        DrawGammaSetMarker(histoHIJINGPi0, 24, 1.5, colorHIJING , colorHIJING);  
        histoHIJINGPi0->SetLineWidth(widthCommonFit);
        histoHIJINGPi0->Draw("same,hist,l");
        DrawGammaSetMarker(histoHIJINGEta, 24, 1.5, colorHIJING , colorHIJING);  
        histoHIJINGEta->SetLineWidth(widthCommonFit);
        histoHIJINGEta->Draw("same,hist,l");

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
        TLatex *labelALICENormUnPaperAll    = new TLatex(0.20,0.20+0.05*1,"Norm. unc. 3.1%");
        SetStyleTLatex( labelALICENormUnPaperAll, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
        labelALICENormUnPaperAll->Draw();
        

        // labels upper right corner
        TLegend* legendXsectionPaperPyBoth  = GetAndSetLegend2(0.71, 0.95-0.04*3, 0.74+0.33, 0.95, textSizeLabelsPixel, 1, "", 43, 0.18);
        legendXsectionPaperPyBoth->AddEntry(histoDPMJetEta,"DPMJet","l");
        legendXsectionPaperPyBoth->AddEntry(histoHIJINGEta,"HIJING","l");
        legendXsectionPaperPyBoth->AddEntry(histoEPOS3Eta,"EPOS3","l");
        legendXsectionPaperPyBoth->AddEntry(graphMcGillPi0,"Shen #it{et al.}","f");
        legendXsectionPaperPyBoth->Draw();

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


    TCanvas* canvasRpPb = new TCanvas("canvasRpPb","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasRpPb,  0.1, 0.01, 0.015, 0.08);
    canvasRpPb->SetLogx();
    
    if (bWCorrection.Contains("Y")){
        TH2F * histo2DRpPb  = new TH2F("histo2DRpPb","histo2DRpPb",1000,0.23,25.,1000,0,1.8);
        SetStyleHistoTH2ForGraphs(histo2DRpPb, "#it{p}_{T} (GeV/#it{c})","#it{R}_{pA}", 0.035,0.04, 0.035,0.04, 0.8,1., 512, 505);
        histo2DRpPb->GetXaxis()->SetLabelOffset(-0.01);
        histo2DRpPb->GetYaxis()->SetRangeUser(0.3,1.62);
        histo2DRpPb->DrawCopy(); 

            DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystPi0, markerStyleComb, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
            graphRpPbCombSystPi0->Draw("E2same");
            graphRpPbCombSystEta->RemovePoint(graphRpPbCombSystEta->GetN()-1);
            DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystEta, markerStyleComb+4, markerSizeComb, kGray+1 , kGray+1, widthLinesBoxes, kTRUE);
            graphRpPbCombSystEta->Draw("E2same");

            DrawGammaLines(0.23, 25 , 1, 1 ,1, kGray, 7);   

            
            TGraphAsymmErrors* graphRpPbCombStatPi0WOXErr   = (TGraphAsymmErrors*)graphRpPbCombStatPi0->Clone("graphRpPbCombStatPi0WOXErr");
            ProduceGraphAsymmWithoutXErrors(graphRpPbCombStatPi0WOXErr);
            DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatPi0WOXErr, markerStyleComb, markerSizeComb, kBlack , kBlack);
            graphRpPbCombStatPi0WOXErr->Draw("p,same,z");
            TGraphAsymmErrors* graphRpPbCombStatEtaWOXErr   = (TGraphAsymmErrors*)graphRpPbCombStatEta->Clone("graphRpPbCombStatEtaWOXErr");
            ProduceGraphAsymmWithoutXErrors(graphRpPbCombStatEtaWOXErr);
            graphRpPbCombStatEtaWOXErr->RemovePoint(graphRpPbCombStatEtaWOXErr->GetN()-1);
            DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatEtaWOXErr, markerStyleComb+4, markerSizeComb, kGray+1 , kGray+1);
            graphRpPbCombStatEtaWOXErr->Draw("p,same,z");


            TLatex *labelEnergyRpPb     = new TLatex(0.15, 0.95-0.04*1, collisionSystempPb.Data());
            SetStyleTLatex( labelEnergyRpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
            labelEnergyRpPb->Draw();
            TLatex *labelALICERpPb  = new TLatex(0.15,0.95-0.04*2,"ALICE");
            SetStyleTLatex( labelALICERpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
            labelALICERpPb->Draw();

            
            TLegend* legendRpPbComb     = GetAndSetLegend2(0.15, 0.95-0.04*1.25*4, 0.35 , 0.95-0.04*1.25*2, textSizeLabelsPixel*0.85); 
            legendRpPbComb->AddEntry(graphRpPbCombSystPi0,"#pi^{0}","pf");
            legendRpPbComb->AddEntry(graphRpPbCombSystEta,"#eta","pf");
            legendRpPbComb->Draw();
            
        canvasRpPb->Update();
        canvasRpPb->Print(Form("%s/Pi0AndEta_RpPb.%s",outputDir.Data(),suffix.Data()));

        histo2DRpPb->DrawCopy(); 

            

            DrawGammaSetMarkerTGraphAsym(graphPi0RpAErrEPS09sDSS5023GeV, widthCommonFit, styleLineDSS, colorDSSBand, colorDSSBand, widthLinesBoxes, kTRUE, colorDSSBand);
            graphPi0RpAErrEPS09sDSS5023GeV->Draw("same,3");

            DrawGammaLines(0.23, 25 , 1, 1 ,1, kGray, 7);   
            DrawGammaNLOTGraph( graphPi0RpACGC5023GeV, widthCommonFit, styleLineCGC, colorCGC);
            DrawGammaNLOTGraph( graphPi0RpAEPS09sDSS5023GeV, widthCommonFit, styleLineDSS, colorDSS);
            
            graphRpPbCombSystPi0->Draw("E2same");
            graphRpPbCombStatPi0WOXErr->Draw("p,same,z");
            graphPi0RpACGC5023GeV->Draw("same,c");
            graphPi0RpAEPS09sDSS5023GeV->Draw("same,c");
            

            labelEnergyRpPb->Draw();
            labelALICERpPb->Draw();

            
            TLegend* legendRpPbPi0          = GetAndSetLegend2(0.15, 0.95-0.04*1.25*3, 0.35 , 0.95-0.04*1.25*2, textSizeLabelsPixel*0.85); 
            legendRpPbPi0->AddEntry(graphRpPbCombSystPi0,"#pi^{0} Data","pf");
//             legendRpPbPi0->AddEntry(graphRpPbCombSystEta,"#eta","pf");
            legendRpPbPi0->Draw();

            TLegend* legendRpPbPi0Theory    = GetAndSetLegend2(0.55, 0.13, 0.95 , 0.13+0.04*5, textSizeLabelsPixel*0.85,1, "", 43, 0.15);
            legendRpPbPi0Theory->AddEntry(graphPi0RpAEPS09sDSS5023GeV,"EPS09s fDSS NLO","l");
            legendRpPbPi0Theory->AddEntry(graphPi0RpAErrEPS09sDSS5023GeV,"fDSS errors","f");
            legendRpPbPi0Theory->AddEntry((TObject*)0,"#scale[0.75]{JHEP 1207 (2012) 073}","");
            legendRpPbPi0Theory->AddEntry(graphPi0RpACGC5023GeV,"CGC","l");
            legendRpPbPi0Theory->AddEntry((TObject*)0,"#scale[0.75]{Phys.Rev. D88 (2013)114020}","");
            legendRpPbPi0Theory->Draw();
            
            histo2DRpPb->Draw("same,axis");
        canvasRpPb->Update();
        canvasRpPb->Print(Form("%s/Pi0WithTheory_RpPb.%s",outputDir.Data(),suffix.Data()));


        histo2DRpPb->DrawCopy(); 

            TGraphAsymmErrors* graphChPiRpPBPubComb     = AddErrorsQuadraticallyTGraph(graphChPiRpPBPubSyst, graphChPiRpPBPubStat);
            TGraphAsymmErrors* graphChKaRpPBPubComb     = AddErrorsQuadraticallyTGraph(graphChKaRpPBPubSyst, graphChKaRpPBPubStat);
            TGraphAsymmErrors* graphChHadRpPBPubComb    = AddErrorsQuadraticallyTGraph(graphChHadRpPBPubSyst, graphChHadRpPBPubStat);
            
            ProduceGraphAsymmWithoutXErrors(graphChPiRpPBPubComb);
            DrawGammaSetMarkerTGraphAsym(graphChPiRpPBPubComb, markerStyleComb+5, markerSizeComb, kBlue-6 , kBlue-6);
            ProduceGraphAsymmWithoutXErrors(graphChHadRpPBPubComb);
            DrawGammaSetMarkerTGraphAsym(graphChHadRpPBPubComb, markerStyleComb+4, markerSizeComb, kGray+1 , kGray+1);
            
            graphRpPbCombSystPi0->Draw("E2same");
            DrawGammaLines(0.23, 25 , 1, 1 ,1, kGray, 7);   

            graphChPiRpPBPubComb->Draw("p,same,z");
            graphChHadRpPBPubComb->Draw("p,same,z");
            graphRpPbCombStatPi0WOXErr->Draw("p,same,z");
            

            labelEnergyRpPb->Draw();
            labelALICERpPb->Draw();

            TLegend* legendRpPbCombCh     = GetAndSetLegend2(0.15, 0.95-0.04*1.25*5, 0.35 , 0.95-0.04*1.25*2, textSizeLabelsPixel*0.85); 
            legendRpPbCombCh->AddEntry(graphRpPbCombSystPi0,"#pi^{0}","pf");
            legendRpPbCombCh->AddEntry(graphChPiRpPBPubComb,"#pi^{#pm}","pe");
            legendRpPbCombCh->AddEntry(graphChHadRpPBPubComb,"h^{#pm}","pe");
            legendRpPbCombCh->Draw();
            
        canvasRpPb->Update();
        canvasRpPb->Print(Form("%s/Pi0WithChargedPionsAndHadron_RpPb.%s",outputDir.Data(),suffix.Data()));

        histo2DRpPb->DrawCopy(); 
            DrawGammaSetMarkerTGraphAsym(graphRpPbCombSystEta, markerStyleComb, markerSizeComb, kBlack , kBlack, widthLinesBoxes, kTRUE);
            DrawGammaSetMarkerTGraphAsym(graphRpPbCombStatEtaWOXErr, markerStyleComb, markerSizeComb, kBlack , kBlack);

            
            ProduceGraphAsymmWithoutXErrors(graphChKaRpPBPubComb);
            DrawGammaSetMarkerTGraphAsym(graphChKaRpPBPubComb, markerStyleComb+5, markerSizeComb, kBlue-6 , kBlue-6);
            DrawGammaSetMarkerTGraphAsym(graphChHadRpPBPubComb,  markerStyleComb+4, markerSizeComb, kGray+1 , kGray+1);
            
            graphRpPbCombSystEta->Draw("E2same");
            DrawGammaLines(0.23, 25 , 1, 1 ,1, kGray, 7);   

            graphChKaRpPBPubComb->Draw("p,same,z");
            graphChHadRpPBPubComb->Draw("p,same,z");
            graphRpPbCombStatEtaWOXErr->Draw("p,same,z");
            

            labelEnergyRpPb->Draw();
            labelALICERpPb->Draw();

            TLegend* legendRpPbEtaCombCh     = GetAndSetLegend2(0.15, 0.95-0.04*1.25*5, 0.35 , 0.95-0.04*1.25*2, textSizeLabelsPixel*0.85); 
            legendRpPbEtaCombCh->AddEntry(graphRpPbCombSystEta,"#eta","pf");
            legendRpPbEtaCombCh->AddEntry(graphChKaRpPBPubComb,"K^{#pm}","pe");
            legendRpPbEtaCombCh->AddEntry(graphChHadRpPBPubComb,"h^{#pm}","pe");
            legendRpPbEtaCombCh->Draw();
            
        canvasRpPb->Update();
        canvasRpPb->Print(Form("%s/EtaWithChargedKaonsAndHadron_RpPb.%s",outputDir.Data(),suffix.Data()));

        
        histo2DRpPb->GetYaxis()->SetRangeUser(0,1.8);
        histo2DRpPb->DrawCopy(); 
        
            for (Int_t i = 10; i> -1; i--){
                if (graphRpPbIndSystPi0[i]){
                    DrawGammaSetMarkerTGraphAsym(graphRpPbIndSystPi0[i], markerStyleDet[i], markerSizeDet[i]*0.5, colorDet[i] , colorDet[i], widthLinesBoxes, kTRUE);
                    graphRpPbIndSystPi0[i]->Draw("E2same");
                }    
            }    
            TGraphAsymmErrors* graphRpPbIndStatPi0WOXErr[11]    = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

            for (Int_t i = 10; i> -1; i--){
                if (graphRpPbIndStatPi0[i]){
                    graphRpPbIndStatPi0WOXErr[i]                = (TGraphAsymmErrors*)graphRpPbIndStatPi0[i]->Clone(Form("graphRpPb%sStatPi0WOXErr", nameMeasGlobalLabel[i].Data()));
                    ProduceGraphAsymmWithoutXErrors(graphRpPbIndStatPi0WOXErr[i]);
                    DrawGammaSetMarkerTGraphAsym(graphRpPbIndStatPi0WOXErr[i], markerStyleDet[i], markerSizeDet[i]*0.5, colorDet[i] , colorDet[i]);
                    graphRpPbIndStatPi0WOXErr[i]->Draw("p,same,z");
                }    
            }    
            graphRpPbIndStatPi0WOXErr[4]->Draw("p,same,z");

            DrawGammaLines(0.23, 25 , 1, 1 ,1, kGray, 7);   

            labelEnergyRpPb->Draw();
            labelALICERpPb->Draw();
            TLatex *labelPi0RpPb  = new TLatex(0.15,0.95-0.04*3,"#pi^{0} #rightarrow #gamma#gamma");
            SetStyleTLatex( labelPi0RpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
            labelPi0RpPb->Draw();

            
            TLegend* legendRpPbInd     = GetAndSetLegend2(0.45, 0.13, 0.95 , 0.13+0.04*3, textSizeLabelsPixel*0.85,2, "", 43, 0.2);
            for (Int_t i = 0; i< 11; i++){
                if (graphRpPbIndSystPi0[i])legendRpPbInd->AddEntry(graphRpPbIndSystPi0[i],nameMeasGlobalLabel[i].Data(),"pf");
            }
            legendRpPbInd->Draw();
            
        canvasRpPb->Update();
        canvasRpPb->Print(Form("%s/Pi0_RpPb_IndividualMeasurements.%s",outputDir.Data(),suffix.Data()));

        for (Int_t i = 0; i < 11; i++){
            if (graphRpPbIndSystPi0[i]){
                histo2DRpPb->DrawCopy();
                    graphRpPbIndSystPi0[i]->Draw("E2same");
                    graphRpPbIndStatPi0WOXErr[i]->Draw("p,same,z");
                    
                    DrawGammaLines(0.23, 25 , 1, 1 ,1, kGray, 7);   

                    labelEnergyRpPb->Draw();
                    labelALICERpPb->Draw();
                    labelPi0RpPb->Draw();
                    
                canvasRpPb->Update();
                canvasRpPb->Print(Form("%s/Pi0_RpPb_%s.%s",outputDir.Data(), nameMeasGlobalLabel[i].Data(), suffix.Data()));
            }
        }    
        
        histo2DRpPb->DrawCopy(); 
        
            for (Int_t i = 10; i> -1; i--){
                if (graphRpPbIndSystEta[i]){
                    if (i==2)graphRpPbIndSystEta[i]->RemovePoint(graphRpPbIndSystEta[i]->GetN()-1);
                    DrawGammaSetMarkerTGraphAsym(graphRpPbIndSystEta[i], markerStyleDet[i], markerSizeDet[i]*0.5, colorDet[i] , colorDet[i], widthLinesBoxes, kTRUE);
                    graphRpPbIndSystEta[i]->Draw("E2same");
                }    
            }    
            TGraphAsymmErrors* graphRpPbIndStatEtaWOXErr[11]    = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

            for (Int_t i = 10; i> -1; i--){
                if (graphRpPbIndStatEta[i]){
                    if (i==2)graphRpPbIndStatEta[i]->RemovePoint(graphRpPbIndStatEta[i]->GetN()-1);
                    graphRpPbIndStatEtaWOXErr[i]                = (TGraphAsymmErrors*)graphRpPbIndStatEta[i]->Clone(Form("graphRpPb%sStatEtaWOXErr", nameMeasGlobalLabel[i].Data()));
                    ProduceGraphAsymmWithoutXErrors(graphRpPbIndStatEtaWOXErr[i]);
                    DrawGammaSetMarkerTGraphAsym(graphRpPbIndStatEtaWOXErr[i], markerStyleDet[i], markerSizeDet[i]*0.5, colorDet[i] , colorDet[i]);
                    graphRpPbIndStatEtaWOXErr[i]->Draw("p,same,z");
                }    
            }    
            graphRpPbIndStatEtaWOXErr[4]->Draw("p,same,z");
            
            DrawGammaLines(0.23, 25 , 1, 1 ,1, kGray, 7);   

            labelEnergyRpPb->Draw();
            labelALICERpPb->Draw();
            TLatex *labelEtaRpPb  = new TLatex(0.15,0.95-0.04*3,"#eta #rightarrow #gamma#gamma");
            SetStyleTLatex( labelEtaRpPb, textSizeLabelsPixel*0.85,4, 1, 43, kTRUE, 11);
            labelEtaRpPb->Draw();

            
            TLegend* legendRpPbIndEta     = GetAndSetLegend2(0.75, 0.13, 0.95 , 0.13+0.04*3, textSizeLabelsPixel*0.85,1, "", 43, 0.2);
            for (Int_t i = 0; i< 11; i++){
                if (graphRpPbIndSystEta[i])legendRpPbIndEta->AddEntry(graphRpPbIndSystEta[i],nameMeasGlobalLabel[i].Data(),"pf");
            }
            legendRpPbIndEta->Draw();
            
        canvasRpPb->Update();
        canvasRpPb->Print(Form("%s/Eta_RpPb_IndividualMeasurements.%s",outputDir.Data(),suffix.Data()));

        for (Int_t i = 0; i < 11; i++){
            if (graphRpPbIndSystEta[i]){
                histo2DRpPb->DrawCopy();
                    graphRpPbIndSystEta[i]->Draw("E2same");
                    graphRpPbIndStatEtaWOXErr[i]->Draw("p,same,z");
                    
                    DrawGammaLines(0.23, 25 , 1, 1 ,1, kGray, 7);   

                    labelEnergyRpPb->Draw();
                    labelALICERpPb->Draw();
                    labelEtaRpPb->Draw();
                    
                canvasRpPb->Update();
                canvasRpPb->Print(Form("%s/Eta_RpPb_%s.%s",outputDir.Data(), nameMeasGlobalLabel[i].Data(), suffix.Data()));
            }
        }    
    }
    
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
    TGraphErrors* graphRatioPi0PubChPiPHOS                  = CalculateRatioBetweenSpectraWithDifferentBinning( graphIndPi0InvYieldStat[1], graphIndPi0InvYieldSys[1], 
                                                                                                                histoChPiInvYieldPubStat, histoChPiInvYieldPubSyst,  
                                                                                                                kTRUE,  kTRUE, 
                                                                                                                &graphPHOSPi0InvYieldStatRebinnedPubPHOS, &graphPHOSPi0InvYieldSysRebinnedPubPHOS, 
                                                                                                                &graphChPiInvYieldPubStatPubPHOS, &graphChPiInvYieldPubSystPubPHOS )    ;
    cout << "PCM Spectrum - published" << endl;
    TGraphErrors* graphChPiInvYieldPubStatPubPCM            = NULL;
    TGraphErrors* graphChPiInvYieldPubSystPubPCM            = NULL;
    TGraphErrors* graphPCMPi0InvYieldStatRebinnedPubPCM     = NULL;
    TGraphErrors* graphPCMPi0InvYieldSysRebinnedPubPCM      = NULL;
    TGraphErrors* graphRatioPi0PubChPiPCM                   = CalculateRatioBetweenSpectraWithDifferentBinning(  graphIndPi0InvYieldStat[0], graphIndPi0InvYieldSys[0], 
                                                                                                                histoChPiInvYieldPubStat, histoChPiInvYieldPubSyst,  
                                                                                                                kTRUE,  kTRUE, 
                                                                                                                &graphPCMPi0InvYieldStatRebinnedPubPCM, &graphPCMPi0InvYieldSysRebinnedPubPCM, 
                                                                                                                &graphChPiInvYieldPubStatPubPCM, &graphChPiInvYieldPubSystPubPCM )    ;
    cout << "PCM-EMC Spectrum - published" << endl;
    TGraphErrors* graphChPiInvYieldPubStatPubPCMEMC             = NULL;
    TGraphErrors* graphChPiInvYieldPubSystPubPCMEMC             = NULL;
    TGraphErrors* graphPCMEMCPi0InvYieldStatRebinnedPubPCMEMC   = NULL;
    TGraphErrors* graphPCMEMCPi0InvYieldSysRebinnedPubPCMEMC    = NULL;
    TGraphErrors* graphRatioPi0PubChPiPCMEMC                    = CalculateRatioBetweenSpectraWithDifferentBinning( graphIndPi0InvYieldStat[4], graphIndPi0InvYieldSys[4], 
                                                                                                                histoChPiInvYieldPubStat, histoChPiInvYieldPubSyst,  
                                                                                                                kTRUE,  kTRUE, 
                                                                                                                &graphPCMEMCPi0InvYieldStatRebinnedPubPCMEMC, &graphPCMEMCPi0InvYieldSysRebinnedPubPCMEMC, 
                                                                                                                &graphChPiInvYieldPubStatPubPCMEMC, &graphChPiInvYieldPubSystPubPCMEMC )    ;
    cout << "EMC Spectrum - published" << endl;
    TGraphErrors* graphChPiInvYieldPubStatPubEMC            = NULL;
    TGraphErrors* graphChPiInvYieldPubSystPubEMC            = NULL;
    TGraphErrors* graphEMCPi0InvYieldStatRebinnedPubEMC     = NULL;
    TGraphErrors* graphEMCPi0InvYieldSysRebinnedPubEMC      = NULL;
    TGraphErrors* graphRatioPi0PubChPiEMC                   = CalculateRatioBetweenSpectraWithDifferentBinning( graphIndPi0InvYieldStat[2], graphIndPi0InvYieldSys[2], 
                                                                                                                histoChPiInvYieldPubStat, histoChPiInvYieldPubSyst,  
                                                                                                                kTRUE,  kTRUE, 
                                                                                                                &graphEMCPi0InvYieldStatRebinnedPubEMC, &graphEMCPi0InvYieldSysRebinnedPubEMC, 
                                                                                                                &graphChPiInvYieldPubStatPubEMC, &graphChPiInvYieldPubSystPubEMC )    ;

    cout << "Dalitz Spectrum - published" << endl;
    TGraphErrors* graphChPiInvYieldPubStatPubDal            = NULL;
    TGraphErrors* graphChPiInvYieldPubSystPubDal            = NULL;
    TGraphErrors* graphDalPi0InvYieldStatRebinnedPubDal     = NULL;
    TGraphErrors* graphDalPi0InvYieldSysRebinnedPubDal      = NULL;
    TGraphErrors* graphRatioPi0PubChPiDal                   = CalculateRatioBetweenSpectraWithDifferentBinning( graphIndPi0InvYieldStat[5], graphIndPi0InvYieldSys[5], 
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
    SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio, "#it{p}_{T} (GeV/#it{c})","#pi^{0}/#pi^{#pm}", 0.85*textsizeLabelsPPb, textsizeLabelsPPb, 
                              0.85*textsizeLabelsPPb,textsizeLabelsPPb, 0.9, 0.95, 510, 505);
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
 
    TString nameOutputCommonFile    = Form("CombinedResultsPaperPPb5023GeV_%s.root", dateForOutput.Data());
    
    TFile fCombResults(nameOutputCommonFile.Data(), "RECREATE");

    fCombResults.mkdir("Pi0pPb_5.023TeV");
    TDirectoryFile* directoryPi0 = (TDirectoryFile*)fCombResults.Get("Pi0pPb_5.023TeV"); 
    fCombResults.cd("Pi0pPb_5.023TeV");
        
        // Final spectrum 
        graphCombPi0InvYieldTot->Write("graphInvYieldPi0CombpPb5023GeVTotErr");
        graphCombPi0InvYieldStat->Write("graphInvYieldPi0CombpPb5023GeVStatErr");
        graphCombPi0InvYieldSys->Write("graphInvYieldPi0CombpPb5023GeVSysErr");  
         // fits for pi0
        fitInvYieldPi0->Write("TsallisFitPi0");
        fitTCMInvYieldPi0->Write("TwoComponentModelFitPi0");

        // Final inv yield INEL
        if (bWCorrection.Contains("Y")){
            // Final spectrum correlations Method A
            if(graphCombPi0InvYieldTot_yShifted)graphCombPi0InvYieldTot_yShifted->Write("graphInvYieldPi0CombpPb5023GeVTotErr_yShifted");
            if(graphCombPi0InvYieldStat_yShifted)graphCombPi0InvYieldStat_yShifted->Write("graphInvYieldPi0CombpPb5023GeVStatErr_yShifted");
            if(graphCombPi0InvYieldSys_yShifted)graphCombPi0InvYieldSys_yShifted->Write("graphInvYieldPi0CombpPb5023GeVSysErr_yShifted");  
        }
        // writing individual measurements
        for (Int_t i = 0; i< 11; i++){
            if (graphIndPi0InvYieldStat[i]) graphIndPi0InvYieldStat[i]->Write(Form("graphInvYieldPi0%spPb5023GeVStatErr",nameMeasGlobalLabel[i].Data()));
            if (graphIndPi0InvYieldSys[i]) graphIndPi0InvYieldSys[i]->Write(Form("graphInvYieldPi0%spPb5023GeVSysErr",nameMeasGlobalLabel[i].Data()));
            if (bWCorrection.Contains("Y")){
                if (graphIndPi0InvYieldStat_yShifted[i]) graphIndPi0InvYieldStat_yShifted[i]->Write(Form("graphInvYieldPi0%spPb5023GeVStatErr_yShifted",nameMeasGlobalLabel[i].Data()));
                if (graphIndPi0InvYieldSys_yShifted[i]) graphIndPi0InvYieldSys_yShifted[i]->Write(Form("graphInvYieldPi0%spPb5023GeVSysErr_yShifted",nameMeasGlobalLabel[i].Data()));                
            }    
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

        // Final spectrum 
        graphCombEtaInvYieldTot->Write("graphInvYieldEtaCombpPb5023GeVTotErr");
        graphCombEtaInvYieldStat->Write("graphInvYieldEtaCombpPb5023GeVStatErr");
        graphCombEtaInvYieldSys->Write("graphInvYieldEtaCombpPb5023GeVSysErr");  
        // writing Y shifted graphs in addition
        if (bWCorrection.Contains("Y")){
            // Final spectrum 
            if(graphCombEtaInvYieldTot_yShifted)graphCombEtaInvYieldTot_yShifted->Write("graphInvYieldEtaCombpPb5023GeVTotErr_yShifted");
            if(graphCombEtaInvYieldStat_yShifted)graphCombEtaInvYieldStat_yShifted->Write("graphInvYieldEtaCombpPb5023GeVStatErr_yShifted");
            if(graphCombEtaInvYieldSys_yShifted)graphCombEtaInvYieldSys_yShifted->Write("graphInvYieldEtaCombpPb5023GeVSysErr_yShifted");  
        }    
        // fits for eta
        fitInvYieldEta->Write("TsallisFitEta");
        fitTCMInvYieldEta->Write("TwoComponentModelFitEta");
    
        // writing individual measurements
        for (Int_t i = 0; i< 11; i++){
            if (graphIndEtaInvYieldStat[i]) graphIndEtaInvYieldStat[i]->Write(Form("graphInvYieldEta%spPb5023GeVStatErr",nameMeasGlobalLabel[i].Data()));
            if (graphIndEtaInvYieldSys[i]) graphIndEtaInvYieldSys[i]->Write(Form("graphInvYieldEta%spPb5023GeVSysErr",nameMeasGlobalLabel[i].Data()));
            if (bWCorrection.Contains("Y")){
                if (graphIndEtaInvYieldStat_yShifted[i]) graphIndEtaInvYieldStat_yShifted[i]->Write(Form("graphInvYieldEta%spPb5023GeVStatErr_yShifted",nameMeasGlobalLabel[i].Data()));
                if (graphIndEtaInvYieldSys_yShifted[i]) graphIndEtaInvYieldSys_yShifted[i]->Write(Form("graphInvYieldEta%spPb5023GeVSysErr_yShifted",nameMeasGlobalLabel[i].Data()));                
            }    
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
        graphRatioEtaCombCombFitStat->Write("Eta_RatioCombToCombFit_Stat");
        graphRatioEtaCombCombFitSys->Write("Eta_RatioCombToCombFit_Syst");                
        for (Int_t i = 0; i < 11; i++){
            if (graphRatioPi0IndCombFitStat[i]) graphRatioPi0IndCombFitStat[i]->Write();
            if (graphRatioPi0IndCombFitSys[i]) graphRatioPi0IndCombFitSys[i]->Write();
            if (graphRatioEtaIndCombFitStat[i]) graphRatioEtaIndCombFitStat[i]->Write();
            if (graphRatioEtaIndCombFitSys[i]) graphRatioEtaIndCombFitSys[i]->Write();
        }    
        
       
    fCompResults.Close();
    
}
    
