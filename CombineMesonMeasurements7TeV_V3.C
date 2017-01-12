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
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

extern TRandom*    gRandom;
extern TBenchmark* gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*    gMinuit;

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    //    TString name;
};

void CombineMesonMeasurements7TeV_V3(   TString fileNamePCM     = "CombinationInput7TeV/data_PCMResultsFullCorrection_PP_NoBinShifting_v2.root",
                                        TString fileNameEMCAL   = "CombinationInput7TeV/data_EMCalResultsFullCorrection_PP_NoBinShifting.root",
                                        TString fileNamePHOS    = "CombinationInput7TeV/pp7TeV_pass4_ppareek_PHOSResultsFullCorrection_13122016_v2.root",
                                        TString suffix          = "pdf",
                                        TString isMC            = "",
                                        TString thesisPlots     = "",
                                        TString bWCorrection    =""
                                    ){

    TString date                                = ReturnDateString();

    gROOT->Reset();
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();
    SetPlotStyle();

    TString dateForOutput                       = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;
    //___________________________________ Declaration of files _____________________________________________
    TString collisionSystem7TeV                 = "pp, #sqrt{#it{s}} = 7 TeV";

    TString fileNameTheory                      = "ExternalInput/Theory/TheoryCompilationPP.root";
    TString fileNameChargedPionPP               = "ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_20_May_2015.root";
    TString fileNameChargedHadronPP             = "ExternalInput/UnidentifiedCharged/ChargedHadrinSpectraPP_20_May_2015.root";
    TString outputDir                           = Form("%s/%s/CombineMesonMeasurements7TeV%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
    TString nameFinalResDat                     = Form("%s/CombinedResults%s_FitResults.dat",outputDir.Data(),bWCorrection.Data());
    cout << outputDir.Data() << endl;
    cout << fileNamePCM.Data() << endl;

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCALLow.root", fileNameEMCAL.Data(), outputDir.Data()));
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

    Double_t xSection7TeV                       = 62.22*1e-3;
    Double_t xSection7TeVV0AND                  = 54.31*1e-3;
    Double_t xSection7TeVErrUp                  = 2.18;
    Double_t xSection7TeVErrDown                = 2.18;
    Double_t xSection7TeVppINEL                 = 73.2*1e9;
    Double_t recalcBarn                         = 1e12; //NLO in pbarn!!!!

    Width_t widthLinesBoxes                     = 1.4;
    Width_t widthCommonFit                      = 2;

    // Definition of colors, styles and markers sizes
    Color_t colorComb                           = kBlue+2;
    Style_t markerStyleComb                     = 20;
    Size_t  markerSizeComb                      = 2;

    Color_t colorCombLowPt                      = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kFALSE);
    Color_t colorCombHighPt                     = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kTRUE);
    Style_t markerStyleCombLowPt                = 20;
    Style_t markerStyleCombHighPt               = 20;
    Size_t  markerSizeComparison                = 2;

    TString nameMeasGlobal[11]                  = {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged", "PCMOtherDataset"};
    Color_t colorDet[11];
    Color_t colorDetMC[11];
    Style_t markerStyleDet[11];
    Style_t markerStyleDetMC[11];
    Size_t  markerSizeDet[11];
    Size_t  markerSizeDetMC[11];

    Style_t styleMarkerNLOMuHalf                = 24;
    Style_t styleMarkerNLOMuOne                 = 27;
    Style_t styleMarkerNLOMuTwo                 = 30;
    Style_t styleLineNLOMuHalf                  = 8;
    Style_t styleLineNLOMuOne                   = 7;
    Style_t styleLineNLOMuTwo                   = 4;
    Style_t styleLineNLOMuTwoBKK                = 3;
    Style_t styleLineNLOMuTwoDSS                = 6;
    Size_t  sizeMarkerNLO                       = 1;
    Width_t widthLineNLO                        = 2.;

    for (Int_t i = 0; i < 11; i++){
        colorDet[i]                             = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kFALSE, kFALSE, kTRUE);
        colorDetMC[i]                           = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kTRUE, kFALSE, kTRUE);
        markerStyleDet[i]                       = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
        markerStyleDetMC[i]                     = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kTRUE);
        markerSizeDet[i]                        = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kFALSE)*2;
        markerSizeDetMC[i]                      = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kTRUE)*2;
    }

    //************************** Read data for PCM **************************************************
    TFile* filePCM                              = new TFile(fileNamePCM.Data());
    TH1D* histoPCMNumberOfEvents7TeV            = (TH1D*)filePCM->Get("histoNumberOfEvents7TeV");
    TDirectory* directoryPCMPi07TeV             = (TDirectory*)filePCM->Get("Pi07TeV");

        TH1D* histoPCMPi0InvCrossSection7TeV                    = (TH1D*)directoryPCMPi07TeV->Get("InvCrossSectionPi0");
        TGraphAsymmErrors* graphPCMPi0InvCrossSectionSys7TeV    = (TGraphAsymmErrors*)directoryPCMPi07TeV->Get("InvCrossSectionPi0Sys");
        cout << "PCM sys" << endl;
        graphPCMPi0InvCrossSectionSys7TeV->Print();

        TGraphAsymmErrors* graphPCMPi0InvCrossSectionStat7TeV   = new TGraphAsymmErrors(histoPCMPi0InvCrossSection7TeV);
        cout << "PCM stat" << endl;
        graphPCMPi0InvCrossSectionStat7TeV->RemovePoint(0);
        graphPCMPi0InvCrossSectionStat7TeV->Print();

    //************************** Read data for EMCAL ****************************************************
    TFile* fileEMCAL                            = new TFile(fileNameEMCAL.Data());
    TDirectory* directoryEMCALPi07TeV           = (TDirectory*)fileEMCAL->Get("Pi07TeV");

        TH1D* histALPi0InvCrossSection7TeV                  = (TH1D*)directoryEMCALPi07TeV->Get("InvCrossSectionPi0");
        TGraphAsymmErrors* graphEMCALPi0InvCrossSectionStat7TeV = new TGraphAsymmErrors(histALPi0InvCrossSection7TeV);
        cout << "EMCAL stat" <<endl;
        graphEMCALPi0InvCrossSectionStat7TeV->Print();

        TGraphAsymmErrors* graphEMCALPi0InvCrossSectionSys7TeV  = (TGraphAsymmErrors*)directoryEMCALPi07TeV->Get("InvCrossSectionPi0Sys");
        cout << "EMCAL sys" <<endl;
        graphEMCALPi0InvCrossSectionSys7TeV->Print();

    //************************** Read data for PHOS *****************************************************
    TFile* filePHOS                             = new TFile(fileNamePHOS);
        TDirectory* directoryPHOSPi07TeV        = (TDirectory*)filePHOS->Get("Pi07TeV");

        TH1D* histoPHOSPi0InvCrossSection7TeV                   = (TH1D*)directoryPHOSPi07TeV->Get("InvCrossSectionPi0");
        TGraphAsymmErrors* graphPHOSPi0InvCrossSectionStat7TeV  = new TGraphAsymmErrors(histoPHOSPi0InvCrossSection7TeV);
        cout << "PHOS stat" <<endl;
        graphPHOSPi0InvCrossSectionStat7TeV->Print();

        TGraphAsymmErrors* graphPHOSPi0InvCrossSectionSys7TeV   = (TGraphAsymmErrors*)directoryPHOSPi07TeV->Get("InvCrossSectionPi0Sys");
        cout << "PHOS sys" <<endl;
        graphPHOSPi0InvCrossSectionSys7TeV->Print();


    // *******************************************************************************************************
    // ************************** Combination of different measurements **************************************
    // *******************************************************************************************************
    // REMARKS:
    //       - order of measurements defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //       - correlations are defined in CombinePtPointsSpectraFullCorrMat from CombinationFunctions.h
    //       - currently only PCM-EMCAL vs others fully implemeted energy independent
    //       - extendable to other energies
    //       - offsets have to be determined manually, see cout's in shell from combination function


    // definition of array of histograms (NULL - means we have no measurement at this energy for this rec-method)
    // for statistical error and final value from respective method
    TH1D* statErrorCollection[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorCollection[i]                  = NULL;
    }
    statErrorCollection[0]                      = (TH1D*)histoPCMPi0InvCrossSection7TeV->Clone("statErrPCMPi0");
    statErrorCollection[1]                      = (TH1D*)histoPHOSPi0InvCrossSection7TeV->Clone("statErrPHOSPi0");
    statErrorCollection[2]                      = (TH1D*)histALPi0InvCrossSection7TeV->Clone("statErrEMCALPi0");

    // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
    // for systematic error from respective method
    TGraphAsymmErrors* sysErrorCollection[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorCollection[i]                   = NULL;
    }
    sysErrorCollection[0]                       = (TGraphAsymmErrors*)graphPCMPi0InvCrossSectionSys7TeV->Clone("sysErrPCMPi0");
    sysErrorCollection[1]                       = (TGraphAsymmErrors*)graphPHOSPi0InvCrossSectionSys7TeV->Clone("sysErrPHOSPi0");
    sysErrorCollection[2]                       = (TGraphAsymmErrors*)graphEMCALPi0InvCrossSectionSys7TeV->Clone("sysErrEMCALPi0");


    // Definition of final pt binning (has to be set manually)
    Double_t xPtLimits[40]                      =  {0.0, 0.3, 0.4, 0.5, 0.6,
                                                    0.8, 1.0, 1.2, 1.4, 1.6,
                                                    1.8, 2.0, 2.2, 2.4, 2.6,
                                                    2.8, 3.0, 3.2, 3.4, 3.6,
                                                    3.8, 4.0, 4.5, 5.0, 5.5,
                                                    6.0, 7.0, 8.0, 9.0, 10.,
                                                    11., 12., 13., 14., 16.,
                                                    18., 20., 25.
                                                    };

    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match
    Int_t offSets[11]                           =  {0, 5, 4, 0, 0,
                                                    0, 0, 0, 0, 0, 0};
    Int_t offSetsSys[11]                        =  {1, 5, 4, 0, 2,
                                                    0, 0, 0, 0, 0, 0};



    //    **********************************************************************************************************************
    //    ******************************************* Calculation of spectrum including EMCal only *****************************
    //    **********************************************************************************************************************

    TH1D* statErrorRelCollection[11];
    for (Int_t i = 0; i< 11; i++){
        statErrorRelCollection[i] = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (statErrorCollection[i]) statErrorRelCollection[i] = CalculateRelErrUpTH1D( statErrorCollection[i], Form("relativeStatError_%s", nameMeasGlobal[i].Data()));
    }

    TGraphAsymmErrors* sysErrorRelCollection[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorRelCollection[i] = NULL;
    }
    for (Int_t i = 0; i < 11; i++){
        if (sysErrorCollection[i]) sysErrorRelCollection[i] = CalculateRelErrUpAsymmGraph( sysErrorCollection[i], Form("relativeSysError_%s", nameMeasGlobal[i].Data()));
    }

    TGraph* graphWeights[11];
    for (Int_t i = 0; i< 11; i++){
        graphWeights[i] = NULL;
    }

    // Declaration & calculation of combined spectrum
    TString fileNameOutputWeighting                       = Form("%s/Weighting.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombPi0InvCrossSectionStat7TeVPCMEMCPHOS= NULL;
    TGraphAsymmErrors* graphCombPi0InvCrossSectionSys7TeVPCMEMCPHOS = NULL;
    TGraphAsymmErrors* graphCombPi0InvCrossSectionTot7TeV = CombinePtPointsSpectraFullCorrMat( statErrorCollection, sysErrorCollection,
                                                                                                   xPtLimits, 37,
                                                                                                   offSets, offSetsSys,
                                                                                                   graphCombPi0InvCrossSectionStat7TeVPCMEMCPHOS, graphCombPi0InvCrossSectionSys7TeVPCMEMCPHOS,
                                                                                                   fileNameOutputWeighting,1
                                                                                                );
//     return;

    graphCombPi0InvCrossSectionStat7TeVPCMEMCPHOS->Print();

    // Reading weights from output file for plotting
    ifstream fileWeightsRead;
    fileWeightsRead.open(fileNameOutputWeighting,ios_base::in);
    cout << "reading" << fileNameOutputWeighting << endl;
    Double_t xValuesRead[50];
    Double_t weightsRead[11][50];
    Int_t availableMeas[11]                 = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    Int_t nMeasSet                          = 3;
    Int_t nPtBinsRead                       = 0;
    while(!fileWeightsRead.eof() && nPtBinsRead < 50){
        TString garbage                         = "";
        if (nPtBinsRead == 0){
            fileWeightsRead >> garbage ;
            for (Int_t i = 0; i < nMeasSet; i++){
                fileWeightsRead >> availableMeas[i] ;
            }
            cout << "read following measurements: ";
            for (Int_t i = 0; i < 11; i++){
                cout << availableMeas[i] << "\t" ;
            }
            cout << endl;
        } else {
            fileWeightsRead >> xValuesRead[nPtBinsRead-1];
            for (Int_t i = 0; i < nMeasSet; i++){
                fileWeightsRead >> weightsRead[availableMeas[i]][nPtBinsRead-1] ;
            }
            cout << "read: "<<  nPtBinsRead << "\t"<< xValuesRead[nPtBinsRead-1] << "\t" ;
            for (Int_t i = 0; i < nMeasSet; i++){
                cout << weightsRead[availableMeas[i]][nPtBinsRead-1] << "\t";
            }
            cout << endl;
        }
        nPtBinsRead++;
    }
    nPtBinsRead = nPtBinsRead-2 ;
    fileWeightsRead.close();

    for (Int_t i = 0; i < nMeasSet; i++){
        graphWeights[availableMeas[i]]  = new TGraph(nPtBinsRead,xValuesRead,weightsRead[availableMeas[i]]);
        Int_t bin = 0;
        for (Int_t n = 0; n< nPtBinsRead; n++){
            if (graphWeights[availableMeas[i]]->GetY()[bin] == 0) graphWeights[availableMeas[i]]->RemovePoint(bin);
            else bin++;
        }
    }

    //    **********************************************************************************************************************
    //    ******************************************* Plotting weights method only EMC *****************************************
    //    **********************************************************************************************************************
    Int_t textSizeLabelsPixel                   = 900*0.04;

    TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
    canvasWeights->SetLogx();

    TH2F * histo2DWeights;
    histo2DWeights                              = new TH2F("histo2DWeights","histo2DWeights",11000,0.23,70.,1000,-0.5,1.1);
    SetStyleHistoTH2ForGraphs(histo2DWeights, "#it{p}_{T} (GeV/#it{c})","#omega_{a} for BLUE",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DWeights->GetXaxis()->SetMoreLogLabels();
    histo2DWeights->GetXaxis()->SetLabelOffset(-0.01);
    canvasWeights->cd();
    histo2DWeights->Draw("copy");

        TLegend* legendAccWeights           = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSet*1.35), 32);
        for (Int_t i = 0; i < nMeasSet; i++){
            DrawGammaSetMarkerTGraph(graphWeights[availableMeas[i]],
                                    markerStyleDet[availableMeas[i]],
                                    markerSizeDet[availableMeas[i]]*0.5,
                                    colorDet[availableMeas[i]] ,
                                    colorDet[availableMeas[i]]);
            graphWeights[availableMeas[i]]->Draw("p,same,e1");
            legendAccWeights->AddEntry(graphWeights[availableMeas[i]],nameMeasGlobal[availableMeas[i]],"p");
        }
        legendAccWeights->Draw();
        TLatex *labelWeightsEnergy              = new TLatex(0.7,0.20,collisionSystem7TeV.Data());
        SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
        labelWeightsEnergy->SetTextFont(43);
        labelWeightsEnergy->Draw();
        TLatex *labelWeightsPi0                 = new TLatex(0.7,0.16,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelWeightsPi0, 0.85*textSizeLabelsPixel,4);
        labelWeightsPi0->SetTextFont(43);
        labelWeightsPi0->Draw();

        DrawGammaLines(0.23, 70. , 0.5, 0.5,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.4, 0.4,0.1, kGray, 1);
        DrawGammaLines(0.23, 70. , 0.3, 0.3,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.2, 0.2,0.1, kGray, 3);

    canvasWeights->SaveAs(Form("%s/WeightsMethods.%s",outputDir.Data(),suffix.Data()));



    //    **********************************************************************************************************************
    //    ******************************************* Ratio of Comb to Fit ****************************************
    //    **********************************************************************************************************************

    TCanvas* canvasRatioToOldCombined           = new TCanvas("canvasRatioToOldCombined","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToOldCombined, 0.1, 0.02, 0.035, 0.09);
    canvasRatioToOldCombined->SetLogx();

    TH2F * histo2DRatioToOldCombined;
    histo2DRatioToOldCombined                   = new TH2F("histo2DRatioToOldCombined","histo2DRatioToOldCombined",11000,0.23,10.,1000,0.75,1.25);
    SetStyleHistoTH2ForGraphs(histo2DRatioToOldCombined, "#it{p}_{T} (GeV/#it{c})","#frac{Comb}{Comb Old}",0.035,0.04, 0.035,0.04, 1.,1.,510,505);
    histo2DRatioToOldCombined->GetXaxis()->SetMoreLogLabels();
    histo2DRatioToOldCombined->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRatioToOldCombined->Draw("copy");

        DrawGammaLines(0.23, 10. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.23, 10. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.23, 10. , 0.9, 0.9,0.1, kGray, 7);
        DrawGammaLines(0.23, 10. , 1.05, 1.05,0.1, kGray, 3);
        DrawGammaLines(0.23, 10. , 0.95, 0.95,0.1, kGray, 3);

        TLatex *labelRatioToOldEnergy           = new TLatex(0.15,0.89,collisionSystem7TeV.Data());
        SetStyleTLatex( labelRatioToOldEnergy, 0.85*textSizeLabelsPixel,4);
        labelRatioToOldEnergy->SetTextFont(43);
        labelRatioToOldEnergy->Draw();
        TLatex *labelRatioToOldPi0              = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToOldPi0, 0.85*textSizeLabelsPixel,4);
        labelRatioToOldPi0->SetTextFont(43);
        labelRatioToOldPi0->Draw();

    canvasRatioToOldCombined->SaveAs(Form("%s/RatioOfCombToCombOld_PP7TeV.%s",outputDir.Data(),suffix.Data()));


    //     *********************************************************************************************************************
    //     ************************************ Visualize relative errors ******************************************************
    //     *********************************************************************************************************************

    TCanvas* canvasRelSysErr                    = new TCanvas("canvasRelSysErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelSysErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelSysErr->SetLogx();

    TH2F * histo2DRelSysErr;
    histo2DRelSysErr                            = new TH2F("histo2DRelSysErr","histo2DRelSysErr",11000,0.23,70.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelSysErr, "#it{p}_{T} (GeV/#it{c})","sys Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelSysErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelSysErr->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRelSysErr->GetYaxis()->SetRangeUser(0,45.5);
    histo2DRelSysErr->Draw("copy");

        TLegend* legendRelSysErr                = GetAndSetLegend2(0.62, 0.94-(0.035*nMeasSet*1.35), 0.95, 0.94, 32);
        for (Int_t i = 0; i < nMeasSet; i++){
            DrawGammaSetMarkerTGraph(sysErrorRelCollection[availableMeas[i]], markerStyleDet[availableMeas[i]], markerSizeDet[availableMeas[i]]*0.5, colorDet[availableMeas[i]],
                                    colorDet[availableMeas[i]]);
            sysErrorRelCollection[availableMeas[i]]->Draw("p,same,e1");
            legendRelSysErr->AddEntry(sysErrorRelCollection[availableMeas[i]],nameMeasGlobal[availableMeas[i]],"p");
        }
        legendRelSysErr->Draw();

        TLatex *labelRelSysErrEnergy            = new TLatex(0.15,0.89,collisionSystem7TeV.Data());
        SetStyleTLatex( labelRelSysErrEnergy, 0.85*textSizeLabelsPixel,4);
        labelRelSysErrEnergy->SetTextFont(43);
        labelRelSysErrEnergy->Draw();
        TLatex *labelRelSysErrPi0               = new TLatex(0.15,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelSysErrPi0, 0.85*textSizeLabelsPixel,4);
        labelRelSysErrPi0->SetTextFont(43);
        labelRelSysErrPi0->Draw();

    canvasRelSysErr->SaveAs(Form("%s/RelSysErr.%s",outputDir.Data(),suffix.Data()));

    //     *********************************************************************************************************************
    //     ************************************ Visualize relative errors ******************************************************
    //     *********************************************************************************************************************

    TCanvas* canvasRelStatErr                   = new TCanvas("canvasRelStatErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelStatErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelStatErr->SetLogx();

    TH2F * histo2DRelStatErr;
    histo2DRelStatErr                           = new TH2F("histo2DRelStatErr","histo2DRelStatErr",11000,0.23,70.,1000,0,80.5);
    SetStyleHistoTH2ForGraphs(histo2DRelStatErr, "#it{p}_{T} (GeV/#it{c})","stat Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelStatErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelStatErr->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRelStatErr->GetYaxis()->SetRangeUser(0,45.5);
    histo2DRelStatErr->Draw("copy");

        TLegend* legendRelStatErr               = GetAndSetLegend2(0.14, 0.94-(0.035*nMeasSet*1.35), 0.45, 0.94, 32);
        for (Int_t i = 0; i < nMeasSet; i++){
            if (availableMeas[i]== 2){
                DrawGammaSetMarker(statErrorRelCollection[availableMeas[i]], markerStyleDet[availableMeas[i]], markerSizeDet[availableMeas[i]]*0.5, colorDet[availableMeas[i]] ,
                            colorDet[availableMeas[i]]);
                TGraphAsymmErrors* graphDummy     = new TGraphAsymmErrors(statErrorRelCollection[availableMeas[i]]);
                DrawGammaSetMarkerTGraphAsym(graphDummy, markerStyleDet[availableMeas[i]], markerSizeDet[availableMeas[i]]*0.5, colorDet[availableMeas[i]],
                                    colorDet[availableMeas[i]]);
                graphDummy->Draw("same,p,x0");
                legendRelStatErr->AddEntry(graphDummy,nameMeasGlobal[availableMeas[i]],"p");

            } else {
                DrawGammaSetMarker(statErrorRelCollection[availableMeas[i]], markerStyleDet[availableMeas[i]], markerSizeDet[availableMeas[i]]*0.5, colorDet[availableMeas[i]] ,
                            colorDet[availableMeas[i]]);
                statErrorRelCollection[availableMeas[i]]->Draw("p,same,e1");
                legendRelStatErr->AddEntry(statErrorRelCollection[availableMeas[i]],nameMeasGlobal[availableMeas[i]],"p");

            }
        }
        legendRelStatErr->Draw();

        TLatex *labelRelStatErrEnergy           = new TLatex(0.75,0.89,collisionSystem7TeV.Data());
        SetStyleTLatex( labelRelStatErrEnergy, 0.85*textSizeLabelsPixel,4);
        labelRelStatErrEnergy->SetTextFont(43);
        labelRelStatErrEnergy->Draw();
        TLatex *labelRelStatErrPi0              = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelStatErrPi0, 0.85*textSizeLabelsPixel,4);
        labelRelStatErrPi0->SetTextFont(43);
        labelRelStatErrPi0->Draw();

    canvasRelStatErr->SaveAs(Form("%s/RelStatErr.%s",outputDir.Data(),suffix.Data()));


    //************************************************************************************************************************
    //************************************** Comparison sys and stat for new and old combined ********************************
    //************************************************************************************************************************
    TGraphAsymmErrors* graphCombPi0InvCrossSectionRelStat7TeV     = CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionStat7TeVPCMEMCPHOS, "relativeStatError_");
    TGraphAsymmErrors* graphCombPi0InvCrossSectionRelSys7TeV      = CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionSys7TeVPCMEMCPHOS, "relativeSysError_");
    TGraphAsymmErrors* graphCombPi0InvCrossSectionRelTot7TeV      = CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionTot7TeV, "relativeTotalError_");

    TCanvas* canvasRelTotErr                    = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelTotErr->SetLogx();

    TH2F * histo2DRelTotErr;
    histo2DRelTotErr                            = new TH2F("histo2DRelTotErr","histo2DRelTotErr",11000,0.23,70.,1000,0,57.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErr, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErr->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRelTotErr->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelTot7TeV, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
        graphCombPi0InvCrossSectionRelTot7TeV->Draw("p,same,e1");

        TLegend* legendRelTotErr                = GetAndSetLegend2(0.14, 0.94-(0.035*4*1.35), 0.45, 0.94, 32);
        legendRelTotErr->AddEntry(graphCombPi0InvCrossSectionRelTot7TeV,"PCM, PHOS, EMCAL","p");
        legendRelTotErr->Draw();

        TLatex *labelRelTotErrEnergy            = new TLatex(0.75,0.89,collisionSystem7TeV.Data());
        SetStyleTLatex( labelRelTotErrEnergy, 0.85*textSizeLabelsPixel,4);
        labelRelTotErrEnergy->SetTextFont(43);
        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrPi0               = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelTotErrPi0, 0.85*textSizeLabelsPixel,4);
        labelRelTotErrPi0->SetTextFont(43);
        labelRelTotErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/RelTotErr.%s",outputDir.Data(),suffix.Data()));

    histo2DRelSysErr->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelSys7TeV, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
        graphCombPi0InvCrossSectionRelSys7TeV->Draw("p,same,e1");

        legendRelTotErr->Draw();
        labelRelTotErrEnergy->Draw();
        labelRelTotErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/RelSysErr_diffComb.%s",outputDir.Data(),suffix.Data()));

    histo2DRelStatErr->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelStat7TeV, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
        graphCombPi0InvCrossSectionRelStat7TeV->Draw("p,same,e1");

        legendRelTotErr->Draw();
        labelRelTotErrEnergy->Draw();
        labelRelTotErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/RelStatErr_diffComb.%s",outputDir.Data(),suffix.Data()));

    histo2DRelTotErr->GetYaxis()->SetRangeUser(0,57.5);
    histo2DRelTotErr->GetYaxis()->SetTitle("Err (%)");
        histo2DRelTotErr->Draw("copy");

        graphCombPi0InvCrossSectionRelTot7TeV->Draw("p,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelStat7TeV, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombPi0InvCrossSectionRelStat7TeV->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelSys7TeV, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombPi0InvCrossSectionRelSys7TeV->SetLineStyle(7);
        graphCombPi0InvCrossSectionRelSys7TeV->Draw("l,x0,same,e1");

        TLegend* legendRelTotErr2 = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
        legendRelTotErr2->AddEntry(graphCombPi0InvCrossSectionRelTot7TeV,"tot","p");
        legendRelTotErr2->AddEntry(graphCombPi0InvCrossSectionRelStat7TeV,"stat","l");
        legendRelTotErr2->AddEntry(graphCombPi0InvCrossSectionRelSys7TeV,"sys","l");
        legendRelTotErr2->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/Reldecomp.%s",outputDir.Data(),suffix.Data()));



    //    **********************************************************************************************************************
    //    ************************************* Calculating bin shifted spectra & fitting **************************************
    //    **********************************************************************************************************************

    // Cloning spectra
    TGraphAsymmErrors* graphCombPi0InvCrossSectionTot7TeVUnShifted              = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionTot7TeV->Clone("Unshifted");
    TGraphAsymmErrors* graphCombPi0InvCrossSectionStat7TeVPCMEMCPHOSUnShifted   = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionStat7TeVPCMEMCPHOS->Clone("UnshiftedStat");
    TGraphAsymmErrors* graphCombPi0InvCrossSectionSys7TeVPCMEMCPHOSUnShifted    = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionSys7TeVPCMEMCPHOS->Clone("UnshiftedSys");

    TGraphAsymmErrors* graphPCMPi0InvCrossSectionStat7TeVUnShifted              = (TGraphAsymmErrors*)graphPCMPi0InvCrossSectionStat7TeV->Clone("UnshiftedStatPCM");
    TGraphAsymmErrors* graphPCMPi0InvCrossSectionSys7TeVUnShifted               = (TGraphAsymmErrors*)graphPCMPi0InvCrossSectionSys7TeV->Clone("UnshiftedSysPCM");

    TGraphAsymmErrors* graphPHOSPi0InvCrossSectionStat7TeVUnshifted             = (TGraphAsymmErrors*)graphPHOSPi0InvCrossSectionStat7TeV->Clone("UnshiftedStatPHOS");
    TGraphAsymmErrors* graphPHOSPi0InvCrossSectionSys7TeVUnshifted              = (TGraphAsymmErrors*)graphPHOSPi0InvCrossSectionSys7TeV->Clone("UnshiftedSysPHOS");

    TGraphAsymmErrors* graphEMCALPi0InvCrossSectionStat7TeVUnshifted            = (TGraphAsymmErrors*)graphEMCALPi0InvCrossSectionStat7TeV->Clone("UnshiftedStatEMCAL");
    TGraphAsymmErrors* graphEMCALPi0InvCrossSectionSys7TeVUnshifted             = (TGraphAsymmErrors*)graphEMCALPi0InvCrossSectionSys7TeV->Clone("UnshiftedSysEMCAL");

    // Calculating binshifts
    Double_t paramGraph[3]                      = {1.0e12, 8., 0.13};
    TF1* fitInvCrossSectionPi07TeV              = FitObject("l","fitInvCrossSectionPi07TeV","Pi0",histALPi0InvCrossSection7TeV,0.3,25.,paramGraph,"QNRMEX0+");
    TF1* fitInvCrossSectionPi07TeVGraph         = (TF1*)fitInvCrossSectionPi07TeV->Clone("fitInvCrossSectionPi07TeVGraph");

    if(bWCorrection.CompareTo("X")==0 ){
        TF1* fitTsallisPi07TeVPtMult                    = FitObject("tmpt","TsallisMultWithPtPi07TeV","Pi0");
        fitTsallisPi07TeVPtMult->SetParameters(paramGraph[0],paramGraph[1], paramGraph[2]) ; // standard parameter optimize if necessary
        graphCombPi0InvCrossSectionTot7TeV              = ApplyXshift(graphCombPi0InvCrossSectionTot7TeV, fitTsallisPi07TeVPtMult);
        graphCombPi0InvCrossSectionStat7TeVPCMEMCPHOS   = ApplyXshiftIndividualSpectra (graphCombPi0InvCrossSectionTot7TeV,
                                                                                        graphCombPi0InvCrossSectionStat7TeVPCMEMCPHOS,
                                                                                        fitTsallisPi07TeVPtMult,
                                                                                        0, graphCombPi0InvCrossSectionStat7TeVPCMEMCPHOS->GetN());
        graphCombPi0InvCrossSectionSys7TeVPCMEMCPHOS    = ApplyXshiftIndividualSpectra (graphCombPi0InvCrossSectionTot7TeV,
                                                                                        graphCombPi0InvCrossSectionSys7TeVPCMEMCPHOS,
                                                                                        fitTsallisPi07TeVPtMult,
                                                                                        0, graphCombPi0InvCrossSectionSys7TeVPCMEMCPHOS->GetN());
        graphPCMPi0InvCrossSectionStat7TeV              = ApplyXshiftIndividualSpectra( graphCombPi0InvCrossSectionTot7TeV,
                                                                                        graphPCMPi0InvCrossSectionStat7TeV,
                                                                                        fitTsallisPi07TeVPtMult,
                                                                                        0, 25);
        graphPCMPi0InvCrossSectionSys7TeV               = ApplyXshiftIndividualSpectra( graphCombPi0InvCrossSectionTot7TeV,
                                                                                        graphPCMPi0InvCrossSectionSys7TeV,
                                                                                        fitTsallisPi07TeVPtMult,
                                                                                        0, 25);
        graphPHOSPi0InvCrossSectionStat7TeV             = ApplyXshiftIndividualSpectra( graphCombPi0InvCrossSectionTot7TeV,
                                                                                        graphPHOSPi0InvCrossSectionStat7TeV,
                                                                                        fitTsallisPi07TeVPtMult,
                                                                                        4, 25);
        graphPHOSPi0InvCrossSectionSys7TeV              = ApplyXshiftIndividualSpectra( graphCombPi0InvCrossSectionTot7TeV,
                                                                                        graphPHOSPi0InvCrossSectionSys7TeV,
                                                                                        fitTsallisPi07TeVPtMult,
                                                                                        4, 25);
        graphEMCALPi0InvCrossSectionStat7TeV            = ApplyXshiftIndividualSpectra( graphCombPi0InvCrossSectionTot7TeV,
                                                                                        graphEMCALPi0InvCrossSectionStat7TeV,
                                                                                        fitTsallisPi07TeVPtMult,
                                                                                        3, graphCombPi0InvCrossSectionTot7TeV->GetN());
        graphEMCALPi0InvCrossSectionSys7TeV             = ApplyXshiftIndividualSpectra( graphCombPi0InvCrossSectionTot7TeV,
                                                                                        graphEMCALPi0InvCrossSectionSys7TeV,
                                                                                        fitTsallisPi07TeVPtMult,
                                                                                        3, graphCombPi0InvCrossSectionTot7TeV->GetN());


        TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
        DrawGammaCanvasSettings( canvasDummy2,  0.13, 0.01, 0.015, 0.08);
        canvasDummy2->SetLogy();
        canvasDummy2->SetLogx();
        TH2F * histo2DDummy2;
        histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.23,70.,1000,1e-1,10e11);
        SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
        histo2DDummy2->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionStat7TeVPCMEMCPHOS, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
        graphCombPi0InvCrossSectionStat7TeVPCMEMCPHOS->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionTot7TeV, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombPi0InvCrossSectionTot7TeV->Draw("pEsame");

        fitInvCrossSectionPi07TeV->SetLineColor(kBlue+2);
        fitInvCrossSectionPi07TeV->Draw("same");

        canvasDummy2->Update();
        canvasDummy2->Print(Form("%s/ComparisonShiftedPi0_7TeV.%s",outputDir.Data(),suffix.Data()));
    }

    graphCombPi0InvCrossSectionTot7TeV->Fit(fitInvCrossSectionPi07TeV,"QNRMEX0+","",0.4,50.);

    fitInvCrossSectionPi07TeV = FitObject("l","fitInvCrossSectionPi07TeV","Pi0",graphCombPi0InvCrossSectionTot7TeV,0.4,50.,paramGraph,"QNRMEX0+");
    fitInvCrossSectionPi07TeV = FitObject("l","fitInvCrossSectionPi07TeV","Pi0",graphCombPi0InvCrossSectionTot7TeV,0.4,50. ,paramGraph,"QNRMEX0+");

    cout << WriteParameterToFile(fitInvCrossSectionPi07TeV)<< endl;
    Double_t paramTCM[5] = {graphCombPi0InvCrossSectionTot7TeV->GetY()[0],0.3,graphCombPi0InvCrossSectionTot7TeV->GetY()[0]/10000,0.8,3};
    TF1* fitTCMInvCrossSectionPi07TeV = FitObject("tcm","fitTCMInvCrossSectionPi07TeV","Pi0",graphCombPi0InvCrossSectionTot7TeV,0.4,50.,paramTCM,"QNRMEX0+");
    fitTCMInvCrossSectionPi07TeV = FitObject("l","fitPtLevy","Pi0",graphCombPi0InvCrossSectionTot7TeV,0.4,25,NULL,"QNRME+");

    cout << WriteParameterToFile(fitTCMInvCrossSectionPi07TeV)<< endl;

    TString forOutput= WriteParameterToFile(fitInvCrossSectionPi07TeV);
    cout<< forOutput.Data()<< endl;

    
    TGraphAsymmErrors* graphRatioCombCombFitTot7TeV     = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionTot7TeV->Clone();
    graphRatioCombCombFitTot7TeV                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitTot7TeV, fitTCMInvCrossSectionPi07TeV);
    TGraphAsymmErrors* graphRatioCombCombFitStat7TeV    = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionStat7TeVPCMEMCPHOS->Clone();
    graphRatioCombCombFitStat7TeV                       = CalculateGraphErrRatioToFit(graphRatioCombCombFitStat7TeV, fitTCMInvCrossSectionPi07TeV);
    TGraphAsymmErrors* graphRatioCombCombFitSys7TeV     = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionSys7TeVPCMEMCPHOS->Clone();
    graphRatioCombCombFitSys7TeV                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitSys7TeV, fitTCMInvCrossSectionPi07TeV);

    TGraphAsymmErrors* graphRatioPCMCombFitStat7TeV     = (TGraphAsymmErrors*)graphPCMPi0InvCrossSectionStat7TeV->Clone();
    graphRatioPCMCombFitStat7TeV                        = CalculateGraphErrRatioToFit(graphRatioPCMCombFitStat7TeV, fitTCMInvCrossSectionPi07TeV);
    TGraphAsymmErrors* graphRatioPCMCombFitSys7TeV      = (TGraphAsymmErrors*)graphPCMPi0InvCrossSectionSys7TeV->Clone();
    graphRatioPCMCombFitSys7TeV                         = CalculateGraphErrRatioToFit(graphRatioPCMCombFitSys7TeV, fitTCMInvCrossSectionPi07TeV);
    TGraphAsymmErrors* graphRatioPHOSCombFitStat7TeV    = (TGraphAsymmErrors*)graphPHOSPi0InvCrossSectionStat7TeV->Clone();
    graphRatioPHOSCombFitStat7TeV                       = CalculateGraphErrRatioToFit(graphRatioPHOSCombFitStat7TeV, fitTCMInvCrossSectionPi07TeV);
    TGraphAsymmErrors* graphRatioPHOSCombFitSys7TeV     = (TGraphAsymmErrors*)graphPHOSPi0InvCrossSectionSys7TeV->Clone();
    graphRatioPHOSCombFitSys7TeV                        = CalculateGraphErrRatioToFit(graphRatioPHOSCombFitSys7TeV, fitTCMInvCrossSectionPi07TeV);
    TGraphAsymmErrors* graphRatiALCombFitStat7TeV       = (TGraphAsymmErrors*)graphEMCALPi0InvCrossSectionStat7TeV->Clone();
    graphRatiALCombFitStat7TeV                          = CalculateGraphErrRatioToFit(graphRatiALCombFitStat7TeV, fitTCMInvCrossSectionPi07TeV);
    TGraphAsymmErrors* graphRatiALCombFitSys7TeV        = (TGraphAsymmErrors*)graphEMCALPi0InvCrossSectionSys7TeV->Clone();
    graphRatiALCombFitSys7TeV                           = CalculateGraphErrRatioToFit(graphRatiALCombFitSys7TeV, fitTCMInvCrossSectionPi07TeV);


    //  **********************************************************************************************************************
    //  ******************************************* Plot with Fit ****************************************
    //  **********************************************************************************************************************
    textSizeLabelsPixel                         = 48;
    TCanvas* canvasPlotwithFit                  = new TCanvas("canvasPlotwithFit","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasPlotwithFit, 0.12, 0.01, 0.01, 0.11);
    canvasPlotwithFit->SetLogx();
    canvasPlotwithFit->SetLogy();

        Double_t textsizeLabelsPP2              = 0;
        Double_t textsizeFacPP2                 = 0;
        if (canvasPlotwithFit->XtoPixel(canvasPlotwithFit->GetX2()) <canvasPlotwithFit->YtoPixel(canvasPlotwithFit->GetY1()) ){
            textsizeLabelsPP2 = (Double_t)textSizeLabelsPixel/canvasPlotwithFit->XtoPixel(canvasPlotwithFit->GetX2()) ;
            textsizeFacPP2 = (Double_t)1./canvasPlotwithFit->XtoPixel(canvasPlotwithFit->GetX2()) ;
        } else {
            textsizeLabelsPP2 = (Double_t)textSizeLabelsPixel/canvasPlotwithFit->YtoPixel(canvasPlotwithFit->GetY1());
            textsizeFacPP2 = (Double_t)1./canvasPlotwithFit->YtoPixel(canvasPlotwithFit->GetY1());
        }
        cout << textsizeLabelsPP2 << endl;

    TH2F * histo2DCombWithFit;
    histo2DCombWithFit                          = new TH2F("histo2DCombWithFit","histo2DCombWithFit",1000,0.23,70.,1000,0.01,9e13  );
    SetStyleHistoTH2ForGraphs(histo2DCombWithFit, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.85*textsizeLabelsPP2, textsizeLabelsPP2,
                            0.85*textsizeLabelsPP2,textsizeLabelsPP2, 0.9, 0.95, 510, 505);
    histo2DCombWithFit->GetXaxis()->SetMoreLogLabels();
    histo2DCombWithFit->GetXaxis()->SetLabelOffset(-0.01);
    histo2DCombWithFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionTot7TeV, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphCombPi0InvCrossSectionTot7TeV->Draw("E2same");
        fitTCMInvCrossSectionPi07TeV->Draw("p,same,e1");

    canvasPlotwithFit->SaveAs(Form("%s/CombWithFit_PP7TeV.%s",outputDir.Data(),suffix.Data()));

    //    **********************************************************************************************************************
    //    ******************************************* Ratio of Comb to Fit ****************************************
    //    **********************************************************************************************************************
    textSizeLabelsPixel                         = 48;
    TCanvas* canvasRatioToCombFit               = new TCanvas("canvasRatioToCombFit","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatioToCombFit, 0.12, 0.01, 0.01, 0.11);
    canvasRatioToCombFit->SetLogx();

        Double_t textsizeLabelsPP = 0;
        Double_t textsizeFacPP= 0;
        if (canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) <canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1()) ){
            textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
            textsizeFacPP = (Double_t)1./canvasRatioToCombFit->XtoPixel(canvasRatioToCombFit->GetX2()) ;
        } else {
            textsizeLabelsPP = (Double_t)textSizeLabelsPixel/canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
            textsizeFacPP = (Double_t)1./canvasRatioToCombFit->YtoPixel(canvasRatioToCombFit->GetY1());
        }
        cout << textsizeLabelsPP << endl;

    TH2F * histo2DRatioToCombFit;
    histo2DRatioToCombFit                       = new TH2F("histo2DRatioToCombFit","histo2DRatioToCombFit",1000,0.23,70.,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Comb Fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                            0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
    histo2DRatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
    histo2DRatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys7TeV, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
        graphRatioCombCombFitSys7TeV->Draw("E2same");
        DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat7TeV, markerStyleComb, markerSizeComb, colorComb , colorComb);
        graphRatioCombCombFitStat7TeV->Draw("p,same,e1");

        DrawGammaLines(0.23, 70. , 1., 1.,0.1, kGray+2);
        DrawGammaLines(0.23, 70. , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.9, 0.9,0.1, kGray, 7);

        TLatex *labelRatioToFitEnergy           = new TLatex(0.73,0.92,collisionSystem7TeV.Data());
        SetStyleTLatex( labelRatioToFitEnergy, 0.85*textSizeLabelsPixel,4);
        labelRatioToFitEnergy->SetTextFont(43);
        labelRatioToFitEnergy->Draw();
        TLatex *labelRatioToFitPi0              = new TLatex(0.73,0.87,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRatioToFitPi0, 0.85*textSizeLabelsPixel,4);
        labelRatioToFitPi0->SetTextFont(43);
        labelRatioToFitPi0->Draw();


    canvasRatioToCombFit->SaveAs(Form("%s/RatioOfCombToCombFit_PP7TeV.%s",outputDir.Data(),suffix.Data()));

    //    **********************************************************************************************************************
    //    ******************************************* Ratio of Individual meas to Fit ******************************************
    //    **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DRatioToCombFit->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitSys7TeV, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPCMCombFitStat7TeV, markerStyleDet[0] ,markerSizeDet[0]*0.5, colorDet[0], colorDet[0]);
        DrawGammaSetMarkerTGraphAsym(graphRatioPHOSCombFitSys7TeV, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioPHOSCombFitStat7TeV, markerStyleDet[1] ,markerSizeDet[1]*0.5, colorDet[1], colorDet[1]);
        DrawGammaSetMarkerTGraphAsym(graphRatiALCombFitSys7TeV, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2], widthLinesBoxes, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatiALCombFitStat7TeV, markerStyleDet[2] ,markerSizeDet[2]*0.5, colorDet[2], colorDet[2]);

        graphRatioPCMCombFitSys7TeV->Draw("E2same");
        graphRatioPHOSCombFitSys7TeV->Draw("E2same");
        graphRatiALCombFitSys7TeV->Draw("E2same");

        graphRatioPCMCombFitStat7TeV->Draw("p,same,e");
        graphRatioPHOSCombFitStat7TeV->Draw("p,same,e");
        graphRatiALCombFitStat7TeV->Draw("p,same,e");

        DrawGammaLines(0.23, 70. , 1., 1.,0.5, kGray+2);
        DrawGammaLines(0.23, 70. , 1.1, 1.1,0.5, kGray, 7);
        DrawGammaLines(0.23, 70. , 0.9, 0.9,0.5, kGray, 7);

        labelRatioToFitEnergy->Draw();
        labelRatioToFitPi0->Draw();

        //****************************** Definition of the Legend ******************************************
        //**************** Row def ************************
        Double_t rowsLegendOnlyPi0Ratio[5]      = {0.92,0.88,0.84,0.80,0.76};
        Double_t rowsLegendOnlyPi0RatioAbs[5]   = {0.91,2.2,2.1,2.0,1.9};
        Double_t columnsLegendOnlyPi0Ratio[3]   = {0.15,0.32, 0.38};
        Double_t columnsLegendOnlyPi0RatioAbs[3]= {0.15,1.04, 1.37};
        Double_t lengthBox                      = 0.2/2;
        Double_t heightBox                      = 0.08/2;
        //****************** first Column **************************************************
        TLatex *textPCMOnlyRatioPi0             = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[1],"PCM");
        SetStyleTLatex( textPCMOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
        textPCMOnlyRatioPi0->SetTextFont(43);
        textPCMOnlyRatioPi0->Draw();
        TLatex *textPHOSOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[2],"PHOS");
        SetStyleTLatex( textPHOSOnlyRatioPi0,  0.85*textSizeLabelsPixel,4);
        textPHOSOnlyRatioPi0->SetTextFont(43);
        textPHOSOnlyRatioPi0->Draw();
        TLatex *textEMCALOnlyRatioPi0           = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[3],"EMCal");
        SetStyleTLatex( textEMCALOnlyRatioPi0,  0.85*textSizeLabelsPixel,4);
        textEMCALOnlyRatioPi0->SetTextFont(43);
        textEMCALOnlyRatioPi0->Draw();

        //****************** second Column *************************************************
        TLatex *textStatOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[1],rowsLegendOnlyPi0Ratio[0] ,"stat");
        SetStyleTLatex( textStatOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
        textStatOnlyRatioPi0->SetTextFont(43);
        textStatOnlyRatioPi0->Draw();
        TLatex *textSysOnlyRatioPi0             = new TLatex(columnsLegendOnlyPi0Ratio[2] ,rowsLegendOnlyPi0Ratio[0],"syst");
        SetStyleTLatex( textSysOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
        textSysOnlyRatioPi0->SetTextFont(43);
        textSysOnlyRatioPi0->Draw();
        TMarker* markerPCMPi0OnlyRatioPi0       = CreateMarkerFromGraph(graphRatioPCMCombFitSys7TeV,columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[1],1);
        markerPCMPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[1]);
        TMarker* markerPHOSPi0OnlyRatioPi0      = CreateMarkerFromGraph(graphRatioPHOSCombFitSys7TeV, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[2],1);
        markerPHOSPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[2]);
        TMarker* markerEMCALPi0OnlyRatioPi0     = CreateMarkerFromGraph(graphRatiALCombFitSys7TeV, columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[3],1);
        markerEMCALPi0OnlyRatioPi0->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[3]);

        TBox* boxPCMPi0OnlyRatioPi0 = CreateBoxFromGraph(graphRatioPCMCombFitSys7TeV, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[1]- heightBox,
                                                        columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox, rowsLegendOnlyPi0RatioAbs[1]+ heightBox);
        boxPCMPi0OnlyRatioPi0->Draw("l");
        TBox* boxPHOSPi0OnlyRatioPi0 = CreateBoxFromGraph(graphRatioPHOSCombFitSys7TeV, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[2]- heightBox,
                                                        columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox, rowsLegendOnlyPi0RatioAbs[2]+ heightBox);
        boxPHOSPi0OnlyRatioPi0->Draw("l");
        TBox* boxEMCALPi0OnlyRatioPi0 = CreateBoxFromGraph(graphRatiALCombFitSys7TeV, columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[3]- heightBox,
                                                        columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox, rowsLegendOnlyPi0RatioAbs[3]+ heightBox);
        boxEMCALPi0OnlyRatioPi0->Draw("l");

    canvasRatioToCombFit->SaveAs(Form("%s/RatioOfIndividualMeasToCombFit_PP7TeV.%s",outputDir.Data(),suffix.Data()));

}

