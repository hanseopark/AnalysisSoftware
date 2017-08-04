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

void drawLatexAdd(TString latextext, Double_t textcolumn, Double_t textrow, Double_t textSizePixel,Bool_t setFont = kFALSE){
    TLatex *latexDummy                  = new TLatex(textcolumn ,textrow,latextext);
    SetStyleTLatex( latexDummy, textSizePixel,4);
    if(setFont)
        latexDummy->SetTextFont(43);
    latexDummy->Draw();
}

// void CombineMesonMeasurements900GeV_V3(   TString fileNamePCM     = "CombinationInput900GeV/data_PCMResultsFullCorrection_PP_NoBinShifting_newComb.root",
void CombineMesonMeasurements900GeV(    TString fileNamePCM     = "/home/nschmidt/AnalysisResults/pp/900GeV/PCM/data_PCMResultsFullCorrection_PP_20170225.root",
                                        TString fileNamePHOS    = "/home/nschmidt/AnalysisResults/pp/900GeV/PHOS/PHOS_pp_pi0_900GeV_noBWcorr_K0Scorr_20111206.root",
                                        TString fileNameEMCal   = "/home/nschmidt/AnalysisResults/pp/900GeV/EMCal/data_EMCAL-EMCALResultsFullCorrection_PP_20170420.root",
                                        TString fileNamePCMEMCal= "/home/nschmidt/AnalysisResults/pp/900GeV/PCM-EMC/data_PCM-EMCALResultsFullCorrection_PP_20170420.root",
                                        TString suffix          = "pdf",
                                        TString isMC            = "",
                                        TString thesisPlots     = "",
                                        TString bWCorrection    ="",
                                        Int_t numbersofmeas     = 5,
                                        Bool_t useDanielmeas    = kTRUE
                                    ){

    TString date                                = ReturnDateString();

    gROOT->Reset();
    gROOT->SetStyle("Plain");

    StyleSettingsThesis();
    SetPlotStyle();

    TString dateForOutput                       = ReturnDateStringForOutput();
    cout << dateForOutput.Data() << endl;
    //___________________________________ Declaration of files _____________________________________________
    TString collisionSystem900GeV               = "pp, #sqrt{#it{s}} = 900 GeV";
    TString outputDir                           = Form("%s/%s/CombineMesonMeasurements900GeV%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
    cout << outputDir.Data() << endl;
    cout << fileNamePCM.Data() << endl;

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputEMCal.root", fileNameEMCal.Data(), outputDir.Data()));
    gSystem->Exec(Form("cp %s %s/InputPCMEMCal.root", fileNamePCMEMCal.Data(), outputDir.Data()));

    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();

    Width_t widthLinesBoxes                     = 1.4;

    // Definition of colors, styles and markers sizes
    Color_t colorComb                           = kBlue+2;
    Style_t markerStyleComb                     = 20;
    Size_t  markerSizeComb                      = 2;

    TString nameMeasGlobal[11]                  = {"PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged", "PCMOtherDataset"};
    Color_t colorDet[11];
    Color_t colorDetMC[11];
    Style_t markerStyleDet[11];
    Style_t markerStyleDetMC[11];
    Size_t  markerSizeDet[11];
    Size_t  markerSizeDetMC[11];

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

    TFile* inputFile[10];
        inputFile[0]                            = new TFile(fileNamePCM.Data());
        inputFile[1]                            = new TFile(fileNamePHOS.Data());
        inputFile[2]                            = new TFile(fileNameEMCal.Data());
        inputFile[4]                            = new TFile(fileNamePCMEMCal.Data());

    TDirectory* directoryPi0[10];
    TDirectory* directoryEta[10];
            directoryPi0[0]                     = (TDirectory*)inputFile[0]->Get("Pi0900GeV");
            directoryEta[0]                     = (TDirectory*)inputFile[0]->Get("Eta900GeV");
            directoryPi0[2]                     = (TDirectory*)inputFile[2]->Get("Pi0900GeV");
            directoryPi0[4]                     = (TDirectory*)inputFile[4]->Get("Pi0900GeV");

    TH1D* histoNumberOfEvents[10];
    TH1D* histoPi0Mass[10];
    TH1D* histoPi0InvCrossSectionSys[10];
    TH1D* histoPi0FWHMMeV[10];
    TH1D* histoPi0TrueMass[10];
    TH1D* histoPi0TrueFWHMMeV[10];
    TH1D* histoEtaMass[10];
    TH1D* histoEtaFWHMMeV[10];
    TH1D* histoEtaTrueMass[10];
    TH1D* histoEtaTrueFWHMMeV[10];
    TH1D* histoPi0Acc[10];
    TH1D* histoPi0TrueEffPt[10];
    TH1D* histoPi0AccTimesEff[10];
    TH1D* histoPi0InvCrossSection[10];
    TH1D* histoEtaInvCrossSection[10];
    TGraphAsymmErrors* graphPi0InvCrossSectionSys[10];
    TGraphAsymmErrors* graphPi0InvCrossSectionStat[10];
    TGraphAsymmErrors* graphEtaInvCrossSectionStat[10];
    TGraphAsymmErrors* graphEtaInvCrossSectionSys[10];
    TGraphAsymmErrors* graphEtaToPi0PCMStat;
    TGraphAsymmErrors* graphEtaToPi0PCMSys;
    TH1D* histoPi0InvMassSigPlusBG[10];
    TH1D* histoPi0InvMassSig[10];
    TH1D* histoPi0InvMassSigRemBGSub[10];
    TH1D* histoPi0InvMassBG[10];
    TH1D* histoPi0InvMassRemBG[10];
    TH1D* histoPi0InvMassBGTot[10];
    TF1* fitPi0InvMassSig[10];
    TF1* fitPi0InvMassBG[10];
    Bool_t haveAllPi0InvMass[10]                = {kFALSE, kFALSE, kFALSE,kFALSE,kFALSE};
    TString strInvMassBin[10]                   = {"04", "22", "3To3_2", "04",""};
    TH1D* histoEtaToPi0Stat[10];
    TGraphAsymmErrors* graphEtaToPi0Stat[10]    = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
    TGraphAsymmErrors* graphEtaToPi0Sys[10]     = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

    Double_t xSection900GeV                     = 47.78*1e9;
    Double_t rapidityMeas[10]                   = {1.6, 1,1, 1.6,1,1,1};

    for (Int_t i = 0; i < 5; i++){
        if(i!=3){
            // LOAD PHOS INPUT
            if(i==1){
                // load cross section systematics and datapoints
                histoPi0InvCrossSectionSys[i]   = (TH1D*)inputFile[i]->Get("hPi0900GeVSys");
                graphPi0InvCrossSectionSys[i]   = new TGraphAsymmErrors(histoPi0InvCrossSectionSys[i]);
                cout << nameMeasGlobal[i].Data() << " sys:" << endl;
                //                 graphPi0InvCrossSectionSys[i]->Print();

                histoPi0InvCrossSection[i]      = (TH1D*)inputFile[i]->Get("hPi0900GeVStat");
                graphPi0InvCrossSectionStat[i]  = new TGraphAsymmErrors(histoPi0InvCrossSection[i]);
                cout << nameMeasGlobal[i].Data() << " stat:" << endl;
                //                 graphPi0InvCrossSectionStat[i]->Print();
                histoPi0InvCrossSection[i]->Scale(xSection900GeV);
                histoPi0InvCrossSectionSys[i]->Scale(xSection900GeV);
                graphPi0InvCrossSectionStat[i]  = ScaleGraph (graphPi0InvCrossSectionStat[i], xSection900GeV);
                graphPi0InvCrossSectionSys[i]   = ScaleGraph (graphPi0InvCrossSectionSys[i], xSection900GeV);

            } else {
                //______________________________ Neutral pion inputs
                graphPi0InvCrossSectionSys[i]   = (TGraphAsymmErrors*)directoryPi0[i]->Get("InvCrossSectionPi0Sys");
                cout << nameMeasGlobal[i].Data() << " sys:" << endl;
                //                 graphPi0InvCrossSectionSys[i]->Print();

                histoPi0InvCrossSection[i]      = (TH1D*)directoryPi0[i]->Get("InvCrossSectionPi0");
                graphPi0InvCrossSectionStat[i]  = new TGraphAsymmErrors(histoPi0InvCrossSection[i]);
                cout << nameMeasGlobal[i].Data() << " stat:" << endl;
                //                 graphPi0InvCrossSectionStat[i]->Print();

                histoPi0Mass[i]                     = (TH1D*)directoryPi0[i]->Get("Pi0_Mass_data_INT1");
                histoPi0FWHMMeV[i]                  = (TH1D*)directoryPi0[i]->Get("Pi0_Width_data_INT1");
                histoPi0TrueMass[i]                 = (TH1D*)directoryPi0[i]->Get("Pi0_Mass_MC_INT1");
                histoPi0TrueFWHMMeV[i]              = (TH1D*)directoryPi0[i]->Get("Pi0_Width_MC_INT1");
                histoPi0Mass[i]                     ->Scale(1000);
                histoPi0FWHMMeV[i]                  ->Scale(1000);
                histoPi0TrueMass[i]                 ->Scale(1000);
                histoPi0TrueFWHMMeV[i]              ->Scale(1000);
                // load acceptance and efficiency and calculate acc*eff*y*2pi
                histoPi0Acc[i]                      = (TH1D*)directoryPi0[i]->Get("AcceptancePi0_INT1");
                histoPi0TrueEffPt[i]                = (TH1D*)directoryPi0[i]->Get("EfficiencyPi0_INT1");
                histoPi0AccTimesEff[i]              = (TH1D*)histoPi0TrueEffPt[i]->Clone(Form("histoPi0AccTimesEff%s",nameMeasGlobal[i].Data()));
                histoPi0AccTimesEff[i]->Multiply(histoPi0Acc[i]);
                histoPi0AccTimesEff[i]->Scale(2*TMath::Pi()*rapidityMeas[i]);

                // load invariant mass example bins
                histoPi0InvMassSig[i]               = (TH1D*)directoryPi0[i]->Get("Pi0_InvMassSig_Example_INT1");
                histoPi0InvMassSigPlusBG[i]         = (TH1D*)directoryPi0[i]->Get("Pi0_InvMassSigPlusBG_Example_INT1");
                histoPi0InvMassBG[i]                = (TH1D*)directoryPi0[i]->Get("Pi0_InvMassBG_Example_INT1");
                fitPi0InvMassSig[i]                 = (TF1*)directoryPi0[i]->Get("Pi0_InvMassSigFit_Example_INT1");
                if (histoPi0InvMassSig[i] && histoPi0InvMassSigPlusBG[i] && histoPi0InvMassBG[i] && fitPi0InvMassSig[i]){
                    haveAllPi0InvMass[i]            = kTRUE;
                    histoPi0InvMassBGTot[i]         = (TH1D*)histoPi0InvMassBG[i]->Clone(Form("Pi0_InvMassTotBG_Example_%s",nameMeasGlobal[i].Data()));
                    histoPi0InvMassSigRemBGSub[i]   = (TH1D*)histoPi0InvMassSig[i]->Clone(Form("Pi0_InvMassSigRemBGSub_Example_%s",nameMeasGlobal[i].Data()));
                }
                //______________________________ Eta meson inputs
                if(i==0){
                    histoEtaInvCrossSection[i]      = (TH1D*)directoryEta[i]->Get("InvCrossSectionEta");
                    graphEtaInvCrossSectionStat[i]  = (TGraphAsymmErrors*)directoryEta[i]->Get("graphInvCrossSectionEta");
                    graphEtaInvCrossSectionSys[i]   = (TGraphAsymmErrors*)directoryEta[i]->Get("InvCrossSectionEtaSys");
                    graphEtaToPi0PCMStat            = (TGraphAsymmErrors*)directoryEta[i]->Get("graphEtaToPi0StatError");
                    graphEtaToPi0PCMSys             = (TGraphAsymmErrors*)directoryEta[i]->Get("EtaToPi0SystError");
                }

            }
        }
    }

    cout << "Finished loading inputs" << endl;
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
    for (Int_t i = 0; i< numbersofmeas; i++){
        if(i!=3)
        statErrorCollection[i]                  = (TH1D*)histoPi0InvCrossSection[i]->Clone(Form("statErr%sPi0",nameMeasGlobal[i].Data()));
    }

    // definition of array of TGraphAsymmErrors (NULL - means we have no measurement at this energy for this rec-method)
    // for systematic error from respective method
    TGraphAsymmErrors* sysErrorCollection[11];
    for (Int_t i = 0; i< 11; i++){
        sysErrorCollection[i]                   = NULL;
    }
    for (Int_t i = 0; i< numbersofmeas; i++){
        if(i!=3)
        sysErrorCollection[i]                   = (TGraphAsymmErrors*)graphPi0InvCrossSectionSys[i]->Clone(Form("sysErr%sPi0",nameMeasGlobal[i].Data()));
    }


    // Definition of final pt binning (has to be set manually)
    Double_t xPtLimits[51]                      =  {0.0, 0.4, 0.6, 0.8, 1.0,
                                                    1.2, 1.4, 1.6, 2.0, 2.5,
                                                    3.0, 3.5, 4.0, 5.0, 7.0,
                                                    10.
                                                    };
    // Definition of offsets for stat & sys see output of function in shell, make sure pt bins match
    Int_t offSets[11]                           =  {0,   2,    5,      2,        3,           0,   0,     0,      0,        0, 0};
    Int_t offSetsSys[11]                        =  {1,   2,    6,      3,        4,           0,   0,     0,      2,        0, 0};
                                                // pcm, phos, emcal, pcmphos, pcmemcal
                                                //  0    1      2       3         4


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
    TString fileNameOutputWeighting             = Form("%s/Weighting.dat",outputDir.Data());
    TGraphAsymmErrors* graphCombinedPi0InvCrossSectionStat= NULL;
    TGraphAsymmErrors* graphCombinedPi0InvCrossSectionSys = NULL;
    TGraphAsymmErrors* graphCombPi0InvCrossSectionTot = CombinePtPointsSpectraFullCorrMat( statErrorCollection, sysErrorCollection,
                                                                                                   xPtLimits, 15,
                                                                                                   offSets, offSetsSys,
                                                                                                   graphCombinedPi0InvCrossSectionStat, graphCombinedPi0InvCrossSectionSys,
                                                                                                   fileNameOutputWeighting,1
                                                                                                );
//     return;
    graphCombinedPi0InvCrossSectionStat->Print();

    // Reading weights from output file for plotting
    ifstream fileWeightsRead;
    fileWeightsRead.open(fileNameOutputWeighting,ios_base::in);
    cout << "reading" << fileNameOutputWeighting << endl;
    Double_t xValuesRead[50];
    Double_t weightsRead[11][50];
    Int_t availableMeas[11]                     = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    Int_t nMeasSet                              = 4;
    Int_t nPtBinsRead                           = 0;
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

        TLegend* legendAccWeights               = GetAndSetLegend2(0.12, 0.14, 0.45, 0.14+(0.035*nMeasSet*1.35), 32);
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
        TLatex *labelWeightsEnergy              = new TLatex(0.7,0.20,collisionSystem900GeV.Data());
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

        TLatex *labelRelSysErrEnergy            = new TLatex(0.15,0.89,collisionSystem900GeV.Data());
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
                TGraphAsymmErrors* graphDummy   = new TGraphAsymmErrors(statErrorRelCollection[availableMeas[i]]);
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

        TLatex *labelRelStatErrEnergy           = new TLatex(0.75,0.89,collisionSystem900GeV.Data());
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
    TGraphAsymmErrors* graphCombPi0InvCrossSectionRelStat900GeV     = CalculateRelErrUpAsymmGraph( graphCombinedPi0InvCrossSectionStat, "relativeStatError_");
    TGraphAsymmErrors* graphCombPi0InvCrossSectionRelSys900GeV      = CalculateRelErrUpAsymmGraph( graphCombinedPi0InvCrossSectionSys, "relativeSysError_");
    TGraphAsymmErrors* graphCombPi0InvCrossSectionRelTot900GeV      = CalculateRelErrUpAsymmGraph( graphCombPi0InvCrossSectionTot, "relativeTotalError_");

    TCanvas* canvasRelTotErr                    = new TCanvas("canvasRelTotErr","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRelTotErr, 0.08, 0.02, 0.035, 0.09);
    canvasRelTotErr->SetLogx();

    TH2F * histo2DRelTotErr;
    histo2DRelTotErr                            = new TH2F("histo2DRelTotErr","histo2DRelTotErr",11000,0.23,70.,1000,0,57.5);
    SetStyleHistoTH2ForGraphs(histo2DRelTotErr, "#it{p}_{T} (GeV/#it{c})","tot Err (%)",0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DRelTotErr->GetXaxis()->SetMoreLogLabels();
    histo2DRelTotErr->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRelTotErr->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelTot900GeV, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
        graphCombPi0InvCrossSectionRelTot900GeV->Draw("p,same,e1");

        TLegend* legendRelTotErr                = GetAndSetLegend2(0.14, 0.94-(0.035*4*1.35), 0.45, 0.94, 32);
        legendRelTotErr->AddEntry(graphCombPi0InvCrossSectionRelTot900GeV,"PCM, PHOS, EMCAL, PCM-EMCal","p");
        legendRelTotErr->Draw();

        TLatex *labelRelTotErrEnergy            = new TLatex(0.75,0.89,collisionSystem900GeV.Data());
        SetStyleTLatex( labelRelTotErrEnergy, 0.85*textSizeLabelsPixel,4);
        labelRelTotErrEnergy->SetTextFont(43);
        labelRelTotErrEnergy->Draw();
        TLatex *labelRelTotErrPi0               = new TLatex(0.75,0.85,"#pi^{0} #rightarrow #gamma#gamma");
        SetStyleTLatex( labelRelTotErrPi0, 0.85*textSizeLabelsPixel,4);
        labelRelTotErrPi0->SetTextFont(43);
        labelRelTotErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/RelTotErr.%s",outputDir.Data(),suffix.Data()));

    histo2DRelSysErr->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelSys900GeV, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
        graphCombPi0InvCrossSectionRelSys900GeV->Draw("p,same,e1");

        legendRelTotErr->Draw();
        labelRelTotErrEnergy->Draw();
        labelRelTotErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/RelSysErr_diffComb.%s",outputDir.Data(),suffix.Data()));

    histo2DRelStatErr->Draw("copy");

        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelStat900GeV, markerStyleComb+5, markerSizeComb, kGreen+2 , kGreen+2);
        graphCombPi0InvCrossSectionRelStat900GeV->Draw("p,same,e1");

        legendRelTotErr->Draw();
        labelRelTotErrEnergy->Draw();
        labelRelTotErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/RelStatErr_diffComb.%s",outputDir.Data(),suffix.Data()));

    histo2DRelTotErr->GetYaxis()->SetRangeUser(0,57.5);
    histo2DRelTotErr->GetYaxis()->SetTitle("Err (%)");
        histo2DRelTotErr->Draw("copy");

        graphCombPi0InvCrossSectionRelTot900GeV->Draw("p,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelStat900GeV, markerStyleComb, markerSizeComb, colorComb-6 , colorComb-6);
        graphCombPi0InvCrossSectionRelStat900GeV->Draw("l,x0,same,e1");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionRelSys900GeV, markerStyleComb, markerSizeComb, colorComb+2, colorComb+2);
        graphCombPi0InvCrossSectionRelSys900GeV->SetLineStyle(7);
        graphCombPi0InvCrossSectionRelSys900GeV->Draw("l,x0,same,e1");

        TLegend* legendRelTotErr2 = GetAndSetLegend2(0.14, 0.94-(0.035*3*1.35), 0.45, 0.94, 32);
        legendRelTotErr2->AddEntry(graphCombPi0InvCrossSectionRelTot900GeV,"tot","p");
        legendRelTotErr2->AddEntry(graphCombPi0InvCrossSectionRelStat900GeV,"stat","l");
        legendRelTotErr2->AddEntry(graphCombPi0InvCrossSectionRelSys900GeV,"sys","l");
        legendRelTotErr2->Draw();

        labelRelTotErrEnergy->Draw();
        labelRelTotErrPi0->Draw();

    canvasRelTotErr->SaveAs(Form("%s/Reldecomp.%s",outputDir.Data(),suffix.Data()));



    //    **********************************************************************************************************************
    //    ************************************* Calculating bin shifted spectra & fitting **************************************
    //    **********************************************************************************************************************

    // Cloning spectra
    TGraphAsymmErrors* graphCombPi0InvCrossSectionTotUnShifted              = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionTot->Clone("Unshifted");
    TGraphAsymmErrors* graphCombinedPi0InvCrossSectionStatUnShifted   = (TGraphAsymmErrors*)graphCombinedPi0InvCrossSectionStat->Clone("UnshiftedStat");
    TGraphAsymmErrors* graphCombinedPi0InvCrossSectionSysUnShifted    = (TGraphAsymmErrors*)graphCombinedPi0InvCrossSectionSys->Clone("UnshiftedSys");

    TGraphAsymmErrors* graphPi0InvCrossSectionStatUnShifted[10];
    TGraphAsymmErrors* graphPi0InvCrossSectionSysUnShifted[10];
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3){
        graphPi0InvCrossSectionStatUnShifted[i] = (TGraphAsymmErrors*)graphPi0InvCrossSectionStat[i]->Clone(Form("UnshiftedStat%s",nameMeasGlobal[i].Data()));
        graphPi0InvCrossSectionSysUnShifted[i]  = (TGraphAsymmErrors*)graphPi0InvCrossSectionSys[i] ->Clone(Form("UnshiftedSys%s",nameMeasGlobal[i].Data()));
        }
    }
cout << __LINE__ << endl;

    // Calculating binshifts
    Double_t paramGraph[3]                      = {1.0e12, 8., 0.13};
    TF1* fitInvCrossSectionPi0900GeV              = FitObject("l","fitInvCrossSectionPi0900GeV","Pi0",histoPi0InvCrossSection[0],0.3,25.,paramGraph,"QNRMEX0+");
    TF1* fitInvCrossSectionPi0900GeVGraph         = (TF1*)histoPi0InvCrossSection[0]->Clone("fitInvCrossSectionPi0900GeVGraph");

    if(bWCorrection.CompareTo("X")==0 ){
        TF1* fitTsallisPi0900GeVPtMult            = FitObject("tmpt","TsallisMultWithPtPi0900GeV","Pi0");
        fitTsallisPi0900GeVPtMult->SetParameters(paramGraph[0],paramGraph[1], paramGraph[2]) ; // standard parameter optimize if necessary
        graphCombPi0InvCrossSectionTot          = ApplyXshift(graphCombPi0InvCrossSectionTot, fitTsallisPi0900GeVPtMult);
        for (Int_t i = 0; i < numbersofmeas; i++){
            graphPi0InvCrossSectionStat[i]      = ApplyXshiftIndividualSpectra (graphCombPi0InvCrossSectionTot,
                                                                                        graphPi0InvCrossSectionStat[i],
                                                                                        fitTsallisPi0900GeVPtMult,
                                                                                        0, graphPi0InvCrossSectionStat[i]->GetN());
            graphPi0InvCrossSectionSys[i]       = ApplyXshiftIndividualSpectra (graphCombPi0InvCrossSectionTot,
                                                                                        graphPi0InvCrossSectionSys[i],
                                                                                        fitTsallisPi0900GeVPtMult,
                                                                                        0, graphPi0InvCrossSectionSys[i]->GetN());
        }

        TCanvas* canvasDummy2 = new TCanvas("canvasDummy2","",200,10,1200,1100);  // gives the page size
        DrawGammaCanvasSettings( canvasDummy2,  0.13, 0.01, 0.015, 0.08);
        canvasDummy2->SetLogy();
        canvasDummy2->SetLogx();
        TH2F * histo2DDummy2;
        histo2DDummy2 = new TH2F("histo2DDummy2","histo2DDummy2",1000,0.23,19.9,1000,1e-1,10e11);
        SetStyleHistoTH2ForGraphs(histo2DDummy2, "#it{p}_{T} (GeV/#it{c})","#it{E} #frac{d^{3}#sigma}{d#it{p}^{3}} (pb GeV^{-2} #it{c}^{3} )", 0.032,0.04, 0.04,0.04, 1,1.55);
        histo2DDummy2->DrawCopy();

        DrawGammaSetMarkerTGraphAsym(graphCombinedPi0InvCrossSectionStat, 20, 1.5, kRed, kRed, widthLinesBoxes, kTRUE);
        graphCombinedPi0InvCrossSectionStat->Draw("pEsame");
        DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionTot, 24, 1.5, kBlack, kBlack, widthLinesBoxes, kTRUE);
        graphCombPi0InvCrossSectionTot->Draw("pEsame");

        fitInvCrossSectionPi0900GeV->SetLineColor(kBlue+2);
        fitInvCrossSectionPi0900GeV->Draw("same");

        canvasDummy2->Update();
        canvasDummy2->Print(Form("%s/ComparisonShiftedPi0_900GeV.%s",outputDir.Data(),suffix.Data()));
    }

    graphCombPi0InvCrossSectionTot->Fit(fitInvCrossSectionPi0900GeV,"QNRMEX0+","",0.4,50.);
    fitInvCrossSectionPi0900GeV = FitObject("l","fitInvCrossSectionPi0900GeV","Pi0",graphCombPi0InvCrossSectionTot,0.4,50. ,paramGraph,"QNRMEX0+");

    cout << WriteParameterToFile(fitInvCrossSectionPi0900GeV)<< endl;

    TString forOutput= WriteParameterToFile(fitInvCrossSectionPi0900GeV);
    cout<< forOutput.Data()<< endl;


    TGraphAsymmErrors* graphRatioCombCombFitTot900GeV     = (TGraphAsymmErrors*)graphCombPi0InvCrossSectionTot->Clone();
    graphRatioCombCombFitTot900GeV                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitTot900GeV, fitInvCrossSectionPi0900GeV);
    TGraphAsymmErrors* graphRatioCombCombFitStat900GeV    = (TGraphAsymmErrors*)graphCombinedPi0InvCrossSectionStat->Clone();
    graphRatioCombCombFitStat900GeV                       = CalculateGraphErrRatioToFit(graphRatioCombCombFitStat900GeV, fitInvCrossSectionPi0900GeV);
    TGraphAsymmErrors* graphRatioCombCombFitSys900GeV     = (TGraphAsymmErrors*)graphCombinedPi0InvCrossSectionSys->Clone();
    graphRatioCombCombFitSys900GeV                        = CalculateGraphErrRatioToFit(graphRatioCombCombFitSys900GeV, fitInvCrossSectionPi0900GeV);

    TGraphAsymmErrors* graphRatioCombFitStat[10];
    TGraphAsymmErrors* graphRatioCombFitSys[10];
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3){
            graphRatioCombFitStat[i]                = (TGraphAsymmErrors*)graphPi0InvCrossSectionStat[i]->Clone();
            graphRatioCombFitStat[i]                = CalculateGraphErrRatioToFit(graphRatioCombFitStat[i], fitInvCrossSectionPi0900GeV);
            graphRatioCombFitSys[i]                 = (TGraphAsymmErrors*)graphPi0InvCrossSectionSys[i]->Clone();
            graphRatioCombFitSys[i]                 = CalculateGraphErrRatioToFit(graphRatioCombFitSys[i], fitInvCrossSectionPi0900GeV);
        }
    }

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
    histo2DCombWithFit                          = new TH2F("histo2DCombWithFit","histo2DCombWithFit",1000,0.23,19.9,1000,0.01,9e13  );
    SetStyleHistoTH2ForGraphs(histo2DCombWithFit, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.85*textsizeLabelsPP2, textsizeLabelsPP2,
                            0.85*textsizeLabelsPP2,textsizeLabelsPP2, 0.9, 0.95, 510, 505);
    histo2DCombWithFit->GetXaxis()->SetMoreLogLabels();
    histo2DCombWithFit->GetXaxis()->SetLabelOffset(-0.01);
    histo2DCombWithFit->Draw("copy");

    DrawGammaSetMarkerTGraphAsym(graphCombPi0InvCrossSectionTot, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
    graphCombPi0InvCrossSectionTot->Draw("E2same");
    fitInvCrossSectionPi0900GeV->Draw("p,same,e1");

    canvasPlotwithFit->SaveAs(Form("%s/CombWithFit_PP900GeV.%s",outputDir.Data(),suffix.Data()));

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
    histo2DRatioToCombFit                       = new TH2F("histo2DRatioToCombFit","histo2DRatioToCombFit",1000,0.23,19.9,1000,0.2,4.    );
    SetStyleHistoTH2ForGraphs(histo2DRatioToCombFit, "#it{p}_{T} (GeV/#it{c})","#frac{Data}{Comb Fit}", 0.85*textsizeLabelsPP, textsizeLabelsPP,
                            0.85*textsizeLabelsPP,textsizeLabelsPP, 0.9, 0.95, 510, 505);
    histo2DRatioToCombFit->GetXaxis()->SetMoreLogLabels();
    histo2DRatioToCombFit->GetXaxis()->SetLabelOffset(-0.01);
    histo2DRatioToCombFit->GetYaxis()->SetRangeUser(0.05,2.45);
    histo2DRatioToCombFit->Draw("copy");

    DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitSys900GeV, markerStyleComb, markerSizeComb, colorComb , colorComb, widthLinesBoxes, kTRUE);
    graphRatioCombCombFitSys900GeV->Draw("E2same");
    DrawGammaSetMarkerTGraphAsym(graphRatioCombCombFitStat900GeV, markerStyleComb, markerSizeComb, colorComb , colorComb);
    graphRatioCombCombFitStat900GeV->Draw("p,same,e1");

    DrawGammaLines(0.23, 70. , 1., 1.,0.1, kGray+2);
    DrawGammaLines(0.23, 70. , 1.1, 1.1,0.1, kGray, 7);
    DrawGammaLines(0.23, 70. , 0.9, 0.9,0.1, kGray, 7);

    TLatex *labelRatioToFitEnergy           = new TLatex(0.73,0.92,collisionSystem900GeV.Data());
    SetStyleTLatex( labelRatioToFitEnergy, 0.85*textSizeLabelsPixel,4);
    labelRatioToFitEnergy->SetTextFont(43);
    labelRatioToFitEnergy->Draw();
    TLatex *labelRatioToFitPi0              = new TLatex(0.73,0.87,"#pi^{0} #rightarrow #gamma#gamma");
    SetStyleTLatex( labelRatioToFitPi0, 0.85*textSizeLabelsPixel,4);
    labelRatioToFitPi0->SetTextFont(43);
    labelRatioToFitPi0->Draw();


    canvasRatioToCombFit->SaveAs(Form("%s/RatioOfCombToCombFit_PP900GeV.%s",outputDir.Data(),suffix.Data()));

    //    **********************************************************************************************************************
    //    ******************************************* Ratio of Individual meas to Fit ******************************************
    //    **********************************************************************************************************************

    canvasRatioToCombFit->cd();
    histo2DRatioToCombFit->Draw("copy");

    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3){
            DrawGammaSetMarkerTGraphAsym(graphRatioCombFitSys[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i], widthLinesBoxes, kTRUE);
            DrawGammaSetMarkerTGraphAsym(graphRatioCombFitStat[i], markerStyleDet[i] ,markerSizeDet[i]*0.5, colorDet[i], colorDet[i]);
        }
    }
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3)
            graphRatioCombFitSys[i]->Draw("E2same");
    }
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3)
            graphRatioCombFitStat[i]->Draw("p,same,e");
    }

    DrawGammaLines(0.23, 70. , 1., 1.,0.5, kGray+2);
    DrawGammaLines(0.23, 70. , 1.1, 1.1,0.5, kGray, 7);
    DrawGammaLines(0.23, 70. , 0.9, 0.9,0.5, kGray, 7);

    labelRatioToFitEnergy->Draw();
    labelRatioToFitPi0->Draw();

    //****************************** Definition of the Legend ******************************************
    //**************** Row def ************************
    Double_t rowsLegendOnlyPi0Ratio[6]      = {0.92,0.88,0.84,0.80,0.79,0.76};
    Double_t rowsLegendOnlyPi0RatioAbs[6]   = {0.91,2.2,2.1,2.0,1.95,1.9};
    Double_t columnsLegendOnlyPi0Ratio[3]   = {0.15,0.32, 0.38};
    Double_t columnsLegendOnlyPi0RatioAbs[3]= {0.15,.72, 0.9};
    Double_t lengthBox                      = 0.11/2;
    Double_t heightBox                      = 0.08/2;
    //****************** first Column **************************************************
    TLatex *textSingleMeasRatioPi0[10];
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3){
            textSingleMeasRatioPi0[i]           = new TLatex(columnsLegendOnlyPi0Ratio[0],rowsLegendOnlyPi0Ratio[i+1],nameMeasGlobal[i].Data());
            SetStyleTLatex( textSingleMeasRatioPi0[i], 0.85*textSizeLabelsPixel,4);
            textSingleMeasRatioPi0[i]->SetTextFont(43);
            textSingleMeasRatioPi0[i]->Draw();
        }
    }

    //****************** second Column *************************************************
    TLatex *textStatOnlyRatioPi0            = new TLatex(columnsLegendOnlyPi0Ratio[1],rowsLegendOnlyPi0Ratio[0] ,"stat");
    SetStyleTLatex( textStatOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
    textStatOnlyRatioPi0->SetTextFont(43);
    textStatOnlyRatioPi0->Draw();
    TLatex *textSysOnlyRatioPi0             = new TLatex(columnsLegendOnlyPi0Ratio[2] ,rowsLegendOnlyPi0Ratio[0],"syst");
    SetStyleTLatex( textSysOnlyRatioPi0, 0.85*textSizeLabelsPixel,4);
    textSysOnlyRatioPi0->SetTextFont(43);
    textSysOnlyRatioPi0->Draw();

    TMarker* markerPi0OnlyRatio[10];
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3){
            markerPi0OnlyRatio[i]               = CreateMarkerFromGraph(graphRatioCombFitSys[i],columnsLegendOnlyPi0Ratio[1] ,rowsLegendOnlyPi0Ratio[i+1],1);

            markerPi0OnlyRatio[i]->DrawMarker(columnsLegendOnlyPi0RatioAbs[1] ,rowsLegendOnlyPi0RatioAbs[i+1]);
        }
    }

    TBox* boxPi0OnlyRatio[10];
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3){
            boxPi0OnlyRatio[i]                  = CreateBoxFromGraph(graphRatioCombFitSys[i], columnsLegendOnlyPi0RatioAbs[2]-0.5*lengthBox , rowsLegendOnlyPi0RatioAbs[i+1]- heightBox, columnsLegendOnlyPi0RatioAbs[2]+ 3*lengthBox, rowsLegendOnlyPi0RatioAbs[i+1]+ heightBox);
            boxPi0OnlyRatio[i]->Draw("l");
        }
    }

    canvasRatioToCombFit->SaveAs(Form("%s/RatioOfIndividualMeasToCombFit_PP900GeV.%s",outputDir.Data(),suffix.Data()));





    // **********************************************************************************************************************
    // ******************************************* Mass and width for pi0 at 900GeV ****************************************
    // **********************************************************************************************************************
    cout << "plotting pi0 mass" << endl;
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

    TPad* padMassLegend1            = new TPad("padMassLegend1", "", 0.13, 0.32, 0.52, 0.52,-1, -1, -2);
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

    TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, 0.23, 19.9 ,1000., -30, 40);
    SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
                            0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
    histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,24.5);
    histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
    histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
    histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
    histo2DAllPi0FWHM->DrawCopy();

    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoPi0FWHMMeV[i] && histoPi0TrueFWHMMeV[i]&&i!=3&&i!=1){
            DrawGammaSetMarker(histoPi0FWHMMeV[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
            histoPi0FWHMMeV[i]->Draw("p,same,e");
            DrawGammaSetMarker(histoPi0TrueFWHMMeV[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
            histoPi0TrueFWHMMeV[i]->Draw("p,same,e");
        }
    }

    TLatex *labelLegendAMass    = new TLatex(0.13,0.06,"a)");
    SetStyleTLatex( labelLegendAMass, textSizeLabelsPixel,4);
    labelLegendAMass->SetTextFont(43);
    labelLegendAMass->Draw();

    TLatex *labelMassPerf       = new TLatex(0.13,0.87,"ALICE performance");
    SetStyleTLatex( labelMassPerf, textSizeLabelsPixel,4);
    labelMassPerf->SetTextFont(43);
    labelMassPerf->Draw();
    TLatex *labelMassEnergy     = new TLatex(0.13,0.78,collisionSystem900GeV.Data());
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

    TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, 0.23, 19.9, 1000., 125.1, 155.9);
    SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,
                            textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
    histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
    histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
    histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
    histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
    histo2DAllPi0Mass->GetXaxis()->SetLabelOffset(-0.015);
    histo2DAllPi0Mass->DrawCopy();

    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoPi0Mass[i] && histoPi0TrueMass[i]&&i!=3&&i!=1){
            DrawGammaSetMarker(histoPi0Mass[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
            histoPi0Mass[i]->Draw("p,same,e");
            DrawGammaSetMarker(histoPi0TrueMass[i], markerStyleDetMC[i], markerSizeDetMC[i]*0.55, colorDetMC[i] , colorDetMC[i]);
            histoPi0TrueMass[i]->Draw("p,same,e");
        }
    }

    DrawGammaLines(0.23, 50. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);

    TLatex *labelLegendBMass            = new TLatex(0.13,0.22,"b)");
    SetStyleTLatex( labelLegendBMass, textSizeLabelsPixel,4);
    labelLegendBMass->SetTextFont(43);
    labelLegendBMass->Draw();

    //********************************** Defintion of the Legend **************************************************
    Double_t columnsLegendMass2[3]      = {0.,0.57,0.84};
    //         Double_t rowsLegendMass2[5] = {0.8,0.6,0.4,0.2,0.01};
    //         Double_t rowsLegendMass2[6] = {0.84,0.66,0.50,0.33,0.16,0.01};
    //           Double_t  rowsLegendMass2[7]= {0.84,0.66,0.50,0.33,0.01,0.16};
    Double_t  rowsLegendMass2[9]= {0.84,0.66,0.51,0.50,0.331,0.33,0.01,0.16}; //setting for use without PHOS and PCM-PHOS peak positions
    //******************* Offsets ***********************
    Double_t offsetMarkerXMass2         = 0.1;
    Double_t offsetMarkerYMass2         = 0.1;
    //****************** Scale factors ******************
    Double_t scaleMarkerMass2           = 1.2;

    padMassLegend1->cd();
    //****************** first Column **************************************************
    TLatex *textMassPCM[10];
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoPi0Mass[i] && histoPi0TrueMass[i] && histoPi0FWHMMeV[i] && histoPi0TrueFWHMMeV[i]&&i!=3&&i!=1){
            textMassPCM[i]                  = new TLatex(columnsLegendMass2[0],rowsLegendMass2[i+1],nameMeasGlobal[i].Data());
            SetStyleTLatex( textMassPCM[i], textSizeLabelsPixel,4);
            textMassPCM[i]->SetTextFont(43);
            textMassPCM[i]->Draw();
        }
    }
    //****************** second Column *************************************************
    TLatex *textMassData                = new TLatex(columnsLegendMass2[1],rowsLegendMass2[0] ,"Data");
    SetStyleTLatex( textMassData, textSizeLabelsPixel,4);
    textMassData->SetTextFont(43);
    textMassData->Draw();
    TLatex *textMassMC                  = new TLatex(columnsLegendMass2[2] ,rowsLegendMass2[0],"MC");
    SetStyleTLatex( textMassMC, textSizeLabelsPixel,4);
    textMassMC->SetTextFont(43);
    textMassMC->Draw();

    TMarker* markerPCMPi0Mass[10];
    TMarker* markerPCMPi0MassMC[10];
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoPi0Mass[i] && histoPi0TrueMass[i]&&i!=3&&i!=1){
            markerPCMPi0Mass[i]             = CreateMarkerFromHisto(histoPi0Mass[i],columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
            markerPCMPi0Mass[i]->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
            markerPCMPi0MassMC[i]           = CreateMarkerFromHisto(histoPi0TrueMass[i],columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
            markerPCMPi0MassMC[i]->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[i+1]+ offsetMarkerYMass2);
        }
    }

    canvasMassWidthPi0->Update();
    canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth.%s",outputDir.Data(),suffix.Data()));


    // **********************************************************************************************************************
    // ******************************** Acceptance * Efficiency for pi0 single measurement  *********************************
    // **********************************************************************************************************************
    textSizeLabelsPixel             = 55;
    Double_t textSizeLabelsRel      = 55./1200;
    cout << textSizeLabelsRel << endl;

    TCanvas* canvasAcceptanceTimesEff       = new TCanvas("canvasAcceptanceTimesEff", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasAcceptanceTimesEff,  0.1, 0.01, 0.015, 0.095);
    canvasAcceptanceTimesEff->SetLogy(1);
    canvasAcceptanceTimesEff->SetLogx(1);

    TH2F * histo2DAccEff;
    histo2DAccEff                = new TH2F("histo2DAccEff", "histo2DAccEff",1000, 0.23,  19.9, 1000, 8e-5, 3 );
    SetStyleHistoTH2ForGraphs( histo2DAccEff, "#it{p}_{T} (GeV/#it{c})", Form("%s%s","#it{#varepsilon} = 2#pi#upoint#Delta","#it{y}#upoint#it{A}#upoint#it{#varepsilon}_{rec} / #it{P}"),
                            0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 1);//(#times #epsilon_{pur})
                            histo2DAccEff->GetYaxis()->SetLabelOffset(0.001);
                            histo2DAccEff->GetXaxis()->SetLabelOffset(-0.01);
                            histo2DAccEff->GetXaxis()->SetMoreLogLabels(kTRUE);
                            histo2DAccEff->DrawCopy();

    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoPi0AccTimesEff[i]&&i!=3&&i!=1){
            DrawGammaSetMarker(histoPi0AccTimesEff[i], markerStyleDet[i], markerSizeDet[i]*0.55, colorDet[i] , colorDet[i]);
            histoPi0AccTimesEff[i]->Draw("p,same,e");
        }
    }

    TLegend* legendEffiAccPi0           = GetAndSetLegend2(0.59, 0.13, 0.93, 0.13+(3*textSizeLabelsRel),textSizeLabelsPixel);
    for (Int_t i = 0; i < numbersofmeas; i++){
        if(histoPi0AccTimesEff[i]&&i!=3&&i!=1){
            legendEffiAccPi0->AddEntry(histoPi0AccTimesEff[i],nameMeasGlobal[i].Data(),"p");
        }
    }
    legendEffiAccPi0->Draw();

    drawLatexAdd("ALICE performance",0.15,0.92,textSizeLabelsRel,kFALSE);
    drawLatexAdd(collisionSystem900GeV.Data(),0.15,0.87,textSizeLabelsRel,kFALSE);
    drawLatexAdd("#pi^{0} #rightarrow #gamma#gamma",0.15,0.82,textSizeLabelsRel,kFALSE);

    canvasAcceptanceTimesEff->Update();
    canvasAcceptanceTimesEff->Print(Form("%s/Pi0_AcceptanceTimesEff.%s",outputDir.Data(),suffix.Data()));


    // **********************************************************************************************************************
    // **************************Plot example invariant mass bins ***********************************************************
    // **********************************************************************************************************************

    textSizeLabelsPixel                         = 100*3/5;
    TCanvas* canvasInvMassSamplePlot            = new TCanvas("canvasInvMassSamplePlot","",0,0,1500,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasInvMassSamplePlot,  0.09, 0.01, 0.035, 0.08);

    Style_t markerStyleInvMassSGBG              = 0;
    Size_t markerSizeInvMassSGBG                = 0;
    Color_t markerColorInvMassSGBG              = kBlack;
    Style_t markerStyleInvMassMBG               = 24;
    Size_t markerSizeInvMassMBG                 = 1.5;
    Color_t markerColorInvMassMBG               = kGray+2;
    Style_t markerStyleInvMassBG                = 20;
    Size_t markerSizeInvMassBG                  = 2;
    Color_t markerColorInvMassBG                = kBlack;
    Style_t markerStyleInvMassSG                = 20;
    Size_t markerSizeInvMassSG                  = 3;
    Color_t markerColorInvMassSG                = kRed+2;
    Color_t fitColorInvMassSG                   = kAzure+2;

    Double_t marginInvMass                      = 0.1*1500;
    Double_t textsizeLabelsInvMass              = 0;
    Double_t textsizeFacInvMass                 = 0;

    TString strLowerEdgeExamplePi0[10]          = {"0.6","6.0","1.7","0.8","1.6"};
    TString strUpperEdgeExamplePi0[10]          = {"0.8","7.0","1.8","0.9","1.8"};
    if (canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) < canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1())){
        textsizeLabelsInvMass                   = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;
        textsizeFacInvMass                      = (Double_t)1./canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;
    } else {
        textsizeLabelsInvMass                   = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1());
        textsizeFacInvMass                      = (Double_t)1./canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1());
    }
    cout << textsizeLabelsInvMass << endl;

    TH2F * histo2DPi0InvMassDummy;
    histo2DPi0InvMassDummy                      = new TH2F("histo2DPi0InvMassDummy","histo2DPi0InvMassDummy",11000,0.05,0.235,21000,-1000,20000);
    SetStyleHistoTH2ForGraphs(histo2DPi0InvMassDummy, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","Counts",0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,
                                0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,0.88, 0.115/(textsizeFacInvMass*marginInvMass));

    TH2F * histo2DEtaInvMassDummy;
    histo2DEtaInvMassDummy                      = new TH2F("histo2DEtaInvMassDummy","histo2DEtaInvMassDummy",11000,0.35,0.695,21000,-1000,20000);
    SetStyleHistoTH2ForGraphs(histo2DEtaInvMassDummy, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","Counts",0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,
                                0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,0.88, 0.115/(textsizeFacInvMass*marginInvMass));

    for (Int_t i =0 ; i < numbersofmeas; i++){
        if (haveAllPi0InvMass[i]){
            canvasInvMassSamplePlot->cd();
            histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(0.02,0.255);
            histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSub[i]->GetMinimum(),1.15*histoPi0InvMassSigPlusBG[i]->GetMaximum());
            if(i==2)
                histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSub[i]->GetMinimum(),1.35*histoPi0InvMassSigPlusBG[i]->GetMaximum());
            histo2DPi0InvMassDummy->DrawCopy();

            TLatex *labelInvMassPtRange = new TLatex(0.945,0.9,Form("#pi^{0}: %s GeV/#it{c} < #it{p}_{T} < %s GeV/#it{c}",strLowerEdgeExamplePi0[i].Data(),strUpperEdgeExamplePi0[i].Data()));

            DrawGammaSetMarker(histoPi0InvMassSigPlusBG[i], markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
            histoPi0InvMassSigPlusBG[i]->SetLineWidth(1);
            histoPi0InvMassSigPlusBG[i]->Draw("hist,e,same");
            DrawGammaSetMarker(histoPi0InvMassBGTot[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
            histoPi0InvMassBGTot[i]->Draw("same");

            DrawGammaSetMarker(histoPi0InvMassSigRemBGSub[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
            histoPi0InvMassSigRemBGSub[i]->Draw("same");
            fitPi0InvMassSig[i]->SetNpx(1000);
            fitPi0InvMassSig[i]->SetRange(0,0.255);
            fitPi0InvMassSig[i]->SetLineColor(fitColorInvMassSG);
            fitPi0InvMassSig[i]->SetLineWidth(1);
            fitPi0InvMassSig[i]->Draw("same");

            TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9,collisionSystem900GeV.Data());
            SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
            labelInvMassEnergy->SetTextFont(43);
            labelInvMassEnergy->Draw();

            TLatex *labelInvMassTrigger      = new TLatex(0.135,0.9-0.8*textsizeLabelsPP,"MinBias");
            SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4);
            labelInvMassTrigger->SetTextFont(43);
            labelInvMassTrigger->Draw();

            TLatex *labelInvMassReco  = new TLatex(0.135,0.9-2*0.8*textsizeLabelsPP,Form("%s",nameMeasGlobal[i].Data()));
            SetStyleTLatex( labelInvMassReco, 0.85*textSizeLabelsPixel,4);
            labelInvMassReco->SetTextFont(43);
            labelInvMassReco->Draw();

            SetStyleTLatex( labelInvMassPtRange, 0.85*textSizeLabelsPixel,4);
            labelInvMassPtRange->SetTextAlign(31);
            labelInvMassPtRange->SetTextFont(43);
            labelInvMassPtRange->Draw();

            TLegend* legendInvMass  = GetAndSetLegend2(0.67, 0.88-5*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
            legendInvMass->SetMargin(0.25);
            legendInvMass->AddEntry(histoPi0InvMassSigPlusBG[i],"Raw real events","l");
            if(i!=2){
                legendInvMass->AddEntry(histoPi0InvMassBGTot[i],"Mixed event +","p");
                legendInvMass->AddEntry((TObject*)0,"corr. BG","");
            } else{
                legendInvMass->AddEntry(histoPi0InvMassBGTot[i],"Mixed event","p");
            }
            legendInvMass->AddEntry(histoPi0InvMassSigRemBGSub[i],"BG subtracted","p");
            legendInvMass->AddEntry(fitPi0InvMassSig[i], "Fit","l");
            legendInvMass->Draw();
            canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBin_%s.%s",outputDir.Data(), nameMeasGlobal[i].Data(), suffix.Data()));
        } else {
            cout << "missing partial input for invariant mass bin for  for trigger: " << nameMeasGlobal[i].Data() << endl;
        }
    }
    
        TString nameMeasGlobal22[11]                  = {"PCM", "PHOS", "EMCAL", "PCMPHOS", "PCMEMCAL", "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged", "PCMOtherDataset"};

    
    TString nameOutputCommonFile                = Form("%s/CombinedResultsPaperPP900GeV_%s.root", outputDir.Data(), dateForOutput.Data());
    TFile fCombResults(nameOutputCommonFile.Data(), "RECREATE");

    fCombResults.mkdir("Pi0900GeV");
    TDirectoryFile* directoryPi02               = (TDirectoryFile*)fCombResults.Get("Pi0900GeV");
    fCombResults.cd("Pi0900GeV");
    // PCM component
//     graphCombinedPi0InvCrossSectionStat         ->Write("graphInvYieldPi0CombStat");
//     graphCombinedPi0InvCrossSectionSys          ->Write("graphInvYieldPi0CombSys");
    graphCombPi0InvCrossSectionTot         ->Write("graphInvCrossSectionPi0Comb900GeVA");
    graphCombinedPi0InvCrossSectionStat         ->Write("graphInvCrossSectionPi0Comb900GeVAStatErr");
    graphCombinedPi0InvCrossSectionSys          ->Write("graphInvCrossSectionPi0Comb900GeVASysErr");

    TGraphAsymmErrors* graphxsecdummystat;
    TGraphAsymmErrors* graphxsecdummysys;

    for (Int_t i = 0; i < numbersofmeas; i++){
        if(i!=3){
//             graphPi0InvCrossSectionStat[i]      ->Write(Form("graphInvYieldPi0%sStat",nameMeasGlobal[i].Data()));
//             graphPi0InvCrossSectionSys[i]       ->Write(Form("graphInvYieldPi0%sSys",nameMeasGlobal[i].Data()));
            graphPi0InvCrossSectionStat[i]      ->Write(Form("graphInvCrossSectionPi0%s900GeVStatErr",nameMeasGlobal22[i].Data()));
            graphPi0InvCrossSectionSys[i]       ->Write(Form("graphInvCrossSectionPi0%s900GeVSysErr",nameMeasGlobal22[i].Data()));
        }
    }
    fCombResults.mkdir("Eta900GeV");
    TDirectoryFile* directoryEta2               = (TDirectoryFile*)fCombResults.Get("Eta900GeV");
    fCombResults.cd("Eta900GeV");

//     Double_t yValue = 0;
//     for(Int_t i = 0; i< graphEtaInvCrossSectionStat->Get
//     yValue                       = histoStat[meas]->GetBinContent(binCounters[meas]+1-startOffsets[meas]);
//     yStatErr[meas]                          = histoStat[meas]->GetBinError(binCounters[meas]+1-startOffsets[meas]);
//     cout << binCounters[meas]-sysOffsets[meas] << "\t"<< graphSyst[meas]->GetX()[binCounters[meas]-sysOffsets[meas]] 
//         << "\t"<< graphSyst[meas]->GetY()[binCounters[meas]-sysOffsets[meas]]
//         << "\t"<< graphSyst[meas]->GetErrorYlow(binCounters[meas]-sysOffsets[meas]) << endl;
//     ySysErr[meas]                           = graphSyst[meas]->GetErrorYhigh(binCounters[meas]-sysOffsets[meas]);
//     yTotErr[meas]                           = TMath::Sqrt(yStatErr[meas]*yStatErr[meas]+ySysErr[meas]*ySysErr[meas]);
//     
//     
    TGraphAsymmErrors * graphEtaInvCrossSectionTot = (TGraphAsymmErrors*) graphEtaInvCrossSectionStat[0]->Clone("graphInvCrossSectionEtaComb900GeVA");
    for (int j=0;j<graphEtaInvCrossSectionTot->GetN();j++){
                graphEtaInvCrossSectionTot->GetEYhigh()[j] = TMath::Sqrt( graphEtaInvCrossSectionStat[0]->GetEYhigh()[j]*graphEtaInvCrossSectionStat[0]->GetEYhigh()[j] + graphEtaInvCrossSectionSys[0]->GetEYhigh()[j] * graphEtaInvCrossSectionSys[0]->GetEYhigh()[j]);
            }
    
    
    if(graphEtaInvCrossSectionStat[0]){
//         graphEtaInvCrossSectionStat[0]          ->Write("graphInvYieldEtaPCMStat");
//         graphxsecdummystat                      = (TGraphAsymmErrors*)graphEtaInvCrossSectionStat[0]->Clone("graphInvCrossSectionEtaPCM900GeVStatErr");
//         graphxsecdummystat                      = ScaleGraph (graphxsecdummystat, xSection900GeV);
        graphEtaInvCrossSectionTot          ->Write("graphInvCrossSectionEtaComb900GeVA");
        graphEtaInvCrossSectionStat[0]          ->Write("graphInvCrossSectionEtaPCM900GeVStatErr");
    }
    if(graphEtaInvCrossSectionSys[0]){
//         graphEtaInvCrossSectionSys[0]           ->Write("graphInvYieldEtaPCMSys");
//         graphxsecdummysys                       = (TGraphAsymmErrors*)graphEtaInvCrossSectionSys[0]->Clone("graphInvCrossSectionEtaPCM900GeVSysErr");
//         graphxsecdummysys                       = ScaleGraph (graphxsecdummysys, xSection900GeV);
        graphEtaInvCrossSectionSys[0]           ->Write("graphInvCrossSectionEtaPCM900GeVSysErr");
    }

    if(graphEtaToPi0PCMStat) graphEtaToPi0PCMStat   ->Write("graphEtaToPi0PCMStat");
    if(graphEtaToPi0PCMSys) graphEtaToPi0PCMSys     ->Write("graphEtaToPi0PCMSys");
    if(graphEtaToPi0PCMStat) graphEtaToPi0PCMStat   ->Write("graphRatioEtaToPi0PCM900GeVStatErr");
    if(graphEtaToPi0PCMSys) graphEtaToPi0PCMSys     ->Write("graphRatioEtaToPi0PCM900GeVSysErr");

    fCombResults.Close();

}