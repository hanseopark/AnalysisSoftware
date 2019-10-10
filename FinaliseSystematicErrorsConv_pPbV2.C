#include <Riostream.h>
#include <fstream>
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
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"

void FinaliseSystematicErrorsConv_pPbV2(    TString nameDataFileErrors          = "",
                                            TString energy                      = "",
                                            TString meson                       = "",
                                            Int_t numberOfPtBins                = 1 ,
                                            Int_t numberCutStudies              = 1,
                                            Double_t startPtSys                 = 0,
                                            TString additionalName              = "",
                                            TString additionalNameOutput        = "",
                                            TString suffix = "eps"
                                     ){

    StyleSettingsThesis();
    SetPlotStyle();

    // ***************************************************************************************************
    // ***************************** labeling and color settings *****************************************
    // ***************************************************************************************************
    TString date                            = ReturnDateString();
    TString dateForOutput                   = ReturnDateStringForOutput();
    TString collisionSystem                 = ReturnFullCollisionsSystem(energy);
    TString energyForOutput                 = energy;
    energyForOutput.ReplaceAll(".","_");

    // ***************************************************************************************************
    // ******************************* general variable definition  **************************************
    // ***************************************************************************************************
    const Int_t nPtBins                     = numberOfPtBins;
    const Int_t nCuts                       = numberCutStudies;
    Double_t* ptBins                        = 0;
    Double_t* ptBinsErr                     = 0;
    TString nameCutVariation[16];
    TString nameCutVariationSC[16];

    TString nameCutVariationSCCurrent[16]   = { "YieldExtraction"   ,"BGEstimate"   ,"dEdxE"    ,"dEdxPi"       ,"TPCCluster"   ,
                                                "SinglePt"          ,"Qt"           ,"Alpha"    ,"Chi2PsiPair"  ,"CosPoint"     ,
                                                "BG"                ,"MCSmearing"   ,"BGEstimateIterations", "Efficiency"      ,"Eta"     ,
                                                ""
                                              };
    Color_t color[20];
    Color_t markerStyle[20];
    for (Int_t k =0; k<16; k++ ){
        cout << "variation: " << nameCutVariationSCCurrent[k].Data() << endl;
        color[k]                            = GetColorSystematics( nameCutVariationSCCurrent[k] );
        markerStyle[k]                      = GetMarkerStyleSystematics( nameCutVariationSCCurrent[k] );
        nameCutVariation[k]                 = GetSystematicsName(nameCutVariationSCCurrent[k]);
        cout << "Color: " << color[k] << " \t marker: " << markerStyle[k] << endl;
        cout << "name for writing: " << nameCutVariation[k].Data() << endl;
    }

    for (Int_t i = 0; i < numberCutStudies; i++){
        nameCutVariationSC[i]               = nameCutVariationSCCurrent[i];
    }

    if (meson.CompareTo("EtaToPi0") == 0){
        nameCutVariation[0]                 = "Yield extraction #eta";
        nameCutVariation[1]                 = "Yield extraction #pi^{0}";
        nameCutVariationSC[1]               = "YieldExtractionPi0";
        color[1]                            = GetColorSystematics( "YieldExtractionPi0" );
        markerStyle[1]                      = GetMarkerStyleSystematics( "YieldExtractionPi0" );
    }
    if (meson.CompareTo("Pi0Ratio") == 0){
        nameCutVariation[0]                 = "Yield extraction pPb";
        nameCutVariation[1]                 = "Yield extraction pp";
        nameCutVariationSC[1]               = "YieldExtraction";
        color[1]                            = GetColorSystematics( "YieldExtractionPi0" );
        markerStyle[1]                      = GetMarkerStyleSystematics( "YieldExtractionPi0" );
    }

    // ***************************************************************************************************
    // ******************************** Booleans for smoothing *******************************************
    // ***************************************************************************************************
    Bool_t bsmooth[16]                      = { 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothMBPi05TeV[16]             = { 0, 0, 1, 1, 1,  1, 1, 1, 1, 1,
                                                1, 1, 0, 1, 0,  0 };
    Bool_t bsmoothMBEta5TeV[16]             = { 0, 0, 1, 1, 1,  1, 1, 1, 1, 1,
                                                1, 1, 0, 1, 0,  0 };
    Bool_t bsmoothMBEtaToPi05TeV[16]        = { 0, 0, 1, 1, 1,  1, 1, 1, 1, 1,
                                                1, 1, 0, 1, 0,  0 };
    Bool_t bsmoothCentPi05TeV[16]           = { 0, 0, 1, 1, 1,  1, 1, 1, 1, 1,
                                                1, 1, 0, 1, 0,  0 };
    Bool_t bsmoothCentEta5TeV[16]           = { 0, 0, 1, 1, 1,  1, 1, 1, 1, 1,
                                                1, 1, 0, 1, 0,  0 };
    Bool_t bsmoothCentEtaToPi05TeV[16]      = { 0, 0, 1, 1, 1,  1, 1, 1, 1, 1,
                                                1, 1, 0, 1, 0,  0 };

    for (Int_t i = 0; i < numberCutStudies; i++){
        if (energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0 ){
            if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Pi0")==0){
                bsmooth[i]                      = bsmoothMBPi05TeV[i];
            } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Eta")==0){
                bsmooth[i]                      = bsmoothMBEta5TeV[i];
            } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("EtaToPi0")==0){
                bsmooth[i]                      = bsmoothMBEtaToPi05TeV[i];
            // for multiplicity bins
            } else if (additionalNameOutput.CompareTo("") && meson.CompareTo("Pi0")==0){
                bsmooth[i]                      = bsmoothCentPi05TeV[i];
            } else if (additionalNameOutput.CompareTo("") && meson.CompareTo("Eta")==0){
                bsmooth[i]                      = bsmoothCentEta5TeV[i];
            } else if (additionalNameOutput.CompareTo("") && meson.CompareTo("EtaToPi0")==0){
                bsmooth[i]                      = bsmoothCentEtaToPi05TeV[i];
            }
        }
    }

    //****************************************************************************************************
    //**************************** Create output directory ***********************************************
    //****************************************************************************************************
    gSystem->Exec("mkdir -p SystematicErrorsCalculatedConv");

    // ***************************************************************************************************
    // ****************************** Initialize error vectors & graphs **********************************
    // ***************************************************************************************************
    Double_t* errorsNeg[nCuts];
    Double_t errorsNegCorr[nCuts][nPtBins];
    Double_t errorsNegSummed[nPtBins];
    Double_t errorsNegCorrSummed[nPtBins];
    Double_t errorsNegCorrMatSummed[nPtBins];

    Double_t* errorsNegErr[nCuts];
    Double_t errorsNegErrCorr[nCuts][nPtBins];
    Double_t errorsNegErrSummed[nPtBins];
    Double_t errorsNegErrCorrSummed[nPtBins];

    Double_t* errorsPos[nCuts];
    Double_t errorsPosCorr[nCuts][nPtBins];
    Double_t errorsPosSummed[nPtBins];
    Double_t errorsPosCorrSummed[nPtBins];
    Double_t errorsPosCorrMatSummed[nPtBins];

    Double_t* errorsPosErr[nCuts];
    Double_t errorsPosErrSummed[nPtBins];
    Double_t errorsPosErrCorr[nCuts][nPtBins];
    Double_t errorsPosErrCorrSummed[nPtBins];

    Double_t errorsMean[nCuts][nPtBins];
    Double_t errorsMeanCorr[nCuts][nPtBins];
    Double_t errorsMeanSummed[nPtBins];
    Double_t errorsMeanCorrSummed[nPtBins];
    Double_t errorsMeanCorrMatSummed[nPtBins];

    Double_t errorsMeanErr[nCuts][nPtBins];
    Double_t errorsMeanErrCorr[nCuts][nPtBins];
    Double_t errorsMeanErrSummed[nPtBins];
    Double_t errorsMeanErrCorrSummed[nPtBins];
    Double_t errorsMeanErrCorrMatSummed[nPtBins];

    TGraphErrors* negativeErrors[nCuts];
    TGraphErrors* negativeErrorsSummed;
    TGraphErrors* positiveErrors[nCuts];
    TGraphErrors* positiveErrorsSummed;
    TGraphErrors* negativeErrorsCorr[nCuts];
    TGraphErrors* negativeErrorsCorrSummed;
    TGraphErrors* positiveErrorsCorr[nCuts];
    TGraphErrors* positiveErrorsCorrSummed;
    TGraphErrors* meanErrors[nCuts];
    TGraphErrors* meanErrorsSummed;
    TGraphErrors* meanErrorsCorr[nCuts];
    TGraphErrors* meanErrorsCorrSummed;
    TGraphErrors* meanErrorsCorrSummedIncMat;

    for (Int_t l = 0; l < nPtBins; l++){
        errorsPosSummed[l]                  = 0.;
        errorsNegSummed[l]                  = 0.;
        errorsMeanSummed[l]                 = 0.;
        errorsPosCorrSummed[l]              = 0.;
        errorsNegCorrSummed[l]              = 0.;
        errorsMeanCorrSummed[l]             = 0.;
    }

    // ***************************************************************************************************
    // ****************************** Read & process data from file **************************************
    // ***************************************************************************************************
    TFile* fileErrorInput                   = new TFile(nameDataFileErrors);
    for (Int_t i = 0; i < nCuts; i++){
        TGraphAsymmErrors* graphPosErrors   = NULL;
        TGraphAsymmErrors* graphNegErrors   = NULL;
        if ( nameCutVariationSC[i].Contains("BGEstimate") && meson.CompareTo("Pi0Ratio")  ){
            TString nameGraphPos            = Form("%s_SystErrorRel_%s_%s",meson.Data(),nameCutVariationSC[i].Data(), additionalName.Data() );
            TString nameGraphNeg            = Form("%s_SystErrorRel_%s_%s",meson.Data(),nameCutVariationSC[i].Data(), additionalName.Data() );
            if (meson.CompareTo("EtaToPi0")==0){
                nameGraphPos                = Form("Eta_SystErrorRel_%s_%s",nameCutVariationSC[i].Data() , additionalName.Data() );
                nameGraphNeg                = Form("Eta_SystErrorRel_%s_%s",nameCutVariationSC[i].Data() , additionalName.Data() );
            }
            cout << "trying to read: " << nameGraphPos.Data() << "\t" << nameGraphNeg.Data() << endl;
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else if ( nameCutVariationSC[i].CompareTo("YieldExtraction")==0 && meson.CompareTo("EtaToPi0") ){
            TString nameGraphPos            = Form("%s_SystErrorRelPos_%s_%s",meson.Data(), nameCutVariationSC[i].Data(), additionalName.Data() );
            TString nameGraphNeg            = Form("%s_SystErrorRelNeg_%s_%s",meson.Data(), nameCutVariationSC[i].Data(), additionalName.Data() );
            cout << "trying to read: " << nameGraphPos.Data() << "\t" << nameGraphNeg.Data() << endl;
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());

        } else if ( nameCutVariationSC[i].CompareTo("YieldExtraction")==0 && !meson.CompareTo("EtaToPi0") &&  i == 0){
            TString nameGraphPos            = Form("Eta_SystErrorRelPos_%s_%s",nameCutVariationSC[i].Data(), additionalName.Data() );
            TString nameGraphNeg            = Form("Eta_SystErrorRelNeg_%s_%s",nameCutVariationSC[i].Data(), additionalName.Data() );
            cout << "trying to read: " << nameGraphPos.Data() << "\t" << nameGraphNeg.Data() << endl;
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());

        } else if ( nameCutVariationSC[i].CompareTo("YieldExtractionPi0")==0 && !meson.CompareTo("EtaToPi0") &&  i == 1){
            TString nameGraphPos            = Form("Pi0EtaBinning_SystErrorRelPos_%s_%s","YieldExtraction", additionalName.Data() );
            TString nameGraphNeg            = Form("Pi0EtaBinning_SystErrorRelNeg_%s_%s","YieldExtraction", additionalName.Data() );
            cout << "trying to read: " << nameGraphPos.Data() << "\t" << nameGraphNeg.Data() << endl;
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else if ( nameCutVariationSC[i].CompareTo("YieldExtraction")==0 && !meson.CompareTo("Pi0Ratio") &&  i == 0){
            TString nameGraphPos            = Form("Eta_SystErrorRelPos_%s_%s",nameCutVariationSC[i].Data(), additionalName.Data() );
            TString nameGraphNeg            = Form("Eta_SystErrorRelNeg_%s_%s",nameCutVariationSC[i].Data(), additionalName.Data() );
            cout << "trying to read: " << nameGraphPos.Data() << "\t" << nameGraphNeg.Data() << endl;
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());

        } else if ( nameCutVariationSC[i].CompareTo("YieldExtractionPi0")==0 && !meson.CompareTo("Pi0Ratio") &&  i == 1){
            TString nameGraphPos            = Form("Pi0EtaBinning_SystErrorRelPos_%s_%s","YieldExtraction", "pp" );
            TString nameGraphNeg            = Form("Pi0EtaBinning_SystErrorRelNeg_%s_%s","YieldExtraction", "pp" );
            cout << "trying to read: " << nameGraphPos.Data() << "\t" << nameGraphNeg.Data() << endl;
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());

        } else {
            TString nameGraphPos            = Form("%s_SystErrorRelPos_%s%s",meson.Data(),nameCutVariationSC[i].Data(), additionalName.Data() );
            TString nameGraphNeg            = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),nameCutVariationSC[i].Data(), additionalName.Data() );
            cout << "trying to read: " << nameGraphPos.Data() << "\t" << nameGraphNeg.Data() << endl;
            if(!meson.CompareTo("Pi0Ratio") ){
                nameGraphPos            = Form("%s_SystErrorRelPos_%s",meson.Data(),nameCutVariationSC[i].Data());
                nameGraphNeg            = Form("%s_SystErrorRelNeg_%s",meson.Data(),nameCutVariationSC[i].Data() );
            }
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
            if ( graphPosErrors == NULL ){
                cout << "systematic wasn't contained, setting it to 0" << endl;
                TString nameGraphPos            = Form("%s_SystErrorRelPos_%s_%s", meson.Data(),nameCutVariationSC[0].Data(), additionalName.Data() );
                TString nameGraphNeg            = Form("%s_SystErrorRelNeg_%s_%s", meson.Data(),nameCutVariationSC[0].Data(), additionalName.Data() );
                if (meson.CompareTo("EtaToPi0") == 0){
                    nameGraphPos                = Form("Eta_SystErrorRelPos_%s_%s", nameCutVariationSC[0].Data(), additionalName.Data() );
                    nameGraphNeg                = Form("Eta_SystErrorRelNeg_%s_%s", nameCutVariationSC[0].Data(), additionalName.Data() );
                }
                cout << "trying to read: " << nameGraphPos.Data() << "\t" << nameGraphNeg.Data() << endl;
                graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
                graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
                for (Int_t k = 0; k< graphPosErrors->GetN(); k++){
                    graphPosErrors->SetPoint(k, graphPosErrors->GetX()[k],0);
                    graphPosErrors->SetPointEYhigh (k, 0);
                    graphPosErrors->SetPointEYlow (k, 0);
                    graphNegErrors->SetPoint(k, graphNegErrors->GetX()[k],0);
                    graphNegErrors->SetPointEYhigh (k, 0);
                    graphNegErrors->SetPointEYlow (k, 0);
                }
            }
        }
        if (!graphPosErrors) cout << "positive graph, aborting" << endl;
        if (!graphNegErrors) cout << "negative graph, aborting" << endl;
        if (!graphNegErrors || !graphPosErrors)
            return;

        while (graphPosErrors->GetX()[0] < startPtSys){
            graphPosErrors->RemovePoint(0);
            graphNegErrors->RemovePoint(0);
        }
        if (i == 0) {
            ptBins      = graphNegErrors->GetX();
            ptBinsErr   = graphNegErrors->GetEXhigh();
        }
        errorsNeg[i]    = graphNegErrors->GetY();
        errorsNegErr[i] = graphNegErrors->GetEYhigh();
        errorsPos[i]    = graphPosErrors->GetY();
        errorsPosErr[i] = graphPosErrors->GetEYhigh();

        cout << nameCutVariationSC[i].Data() << endl;
        CalculateMeanSysErr(errorsMean[i], errorsMeanErr[i], errorsPos[i], errorsNeg[i], nPtBins);
        CorrectSystematicErrorsWithMean(errorsPos[i],errorsPosErr[i], errorsPosCorr[i], errorsPosErrCorr[i], nPtBins);
        CorrectSystematicErrorsWithMean(errorsNeg[i],errorsNegErr[i], errorsNegCorr[i], errorsNegErrCorr[i], nPtBins);
        CorrectSystematicErrorsWithMean(errorsMean[i], errorsMeanErr[i], errorsMeanCorr[i], errorsMeanErrCorr[i], nPtBins);


        // Routine for manual smoothing of systematic errors
        // ATTTENTION! you have to do this manually for each data set/trigger never trust the values mentioned here

        if (bsmooth[i]){
            Double_t minPt      = -10;
            Double_t errorReset = -10000;

            // Yield extraction - cutstudies nr 0/1
            if (nameCutVariationSC[i].CompareTo("YieldExtraction")==0 && (meson.CompareTo("Eta")==0 || meson.CompareTo("EtaToPi0")==0)){
                for (Int_t k = 0; k < nPtBins; k++){
                    if (errorsMean[i][k] > 20.0){
                        errorsMean[i][k]        = 12.5;
                        errorsMeanErr[i][k]     = 0.125;
                        errorsMeanCorr[i][k]    = 12.5;
                        errorsMeanErrCorr[i][k] = 0.125;
                    }
                }
            }

            // BGEstimate - cutstudies nr 1
            if (nameCutVariationSC[i].CompareTo("BGEstimate")==0  ){
                minPt           = 0;
                errorReset      = 0.3;
            }

            // dEdxE - cutstudies nr 2
            if (nameCutVariationSC[i].CompareTo("dEdxE")==0 ){
                if (meson.CompareTo("Pi0")==0){
                    for (Int_t k = 1; k < nPtBins; k++){
                        errorReset              = errorsMean[i][k];
                        errorReset              = 0.5 + 6.51640e+00/pow(6.76502e+00,ptBins[k]);
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = 0.01*errorReset;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = 0.01*errorReset;
                    }
                } else if (meson.CompareTo("Eta")==0 || meson.CompareTo("EtaToPi0")==0){
                    for (Int_t k = 0; k < nPtBins; k++){
                        errorReset              = errorsMean[i][k];
                        errorReset              = 13.5*pow(0.17,ptBins[k])+0.5+pow(ptBins[k],2)*0.025;
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = 0.01*errorReset;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = 0.01*errorReset;
                    }
                }
            }
            // dEdxPi - cutstudies nr 3
            if (nameCutVariationSC[i].CompareTo("dEdxPi")==0 ){
                if (meson.CompareTo("Pi0")==0){
                    for (Int_t k = 0; k < nPtBins; k++){
                        errorReset              = 1.05+pow(ptBins[k],4)*0.0002;
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = 0.01*errorReset;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = 0.01*errorReset;
                    }
                } else if (meson.CompareTo("Eta")==0 || meson.CompareTo("EtaToPi0")==0 ){
                    for (Int_t k = 0; k < nPtBins; k++){
                        errorReset              = 1.45+pow(ptBins[k],2)*0.03;
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = 0.01*errorReset;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = 0.01*errorReset;
                    }
                }
            }

            // TPCCluster - cutstudies nr 4
            if (nameCutVariationSC[i].CompareTo("TPCCluster")==0 ){
                minPt           = 0;
                if (meson.CompareTo("Pi0")==0 ){
                    errorReset      = 0.05;
                } else {
                    errorReset      = 0.02;
                }
            }

            // SinglePt - cutstudies nr 5
            if (nameCutVariationSC[i].CompareTo("SinglePt")==0 ){
                if (meson.CompareTo("Pi0")==0 ){
                    for (Int_t k = 0; k < nPtBins; k++){
                        errorReset              = 20*pow(0.01,ptBins[k])+0.3+pow(ptBins[k],2)*0.003;
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = 0.01*errorReset;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = 0.01*errorReset;
                    }
                } else if (meson.CompareTo("Eta")==0 ){
                    for (Int_t k = 1; k < nPtBins; k++){
                        errorReset              = 40*pow(0.05,ptBins[k])+0.5+pow(ptBins[k],2)*0.001;
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = 0.01*errorReset;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = 0.01*errorReset;
                    }
                } else if (meson.CompareTo("EtaToPi0")==0 ){
                    for (Int_t k = 1; k < nPtBins; k++){
                        errorReset              = 40*pow(0.05,ptBins[k])+0.5+pow(ptBins[k],2)*0.0025;
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = 0.01*errorReset;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = 0.01*errorReset;
                    }
                }
        }

            // Qt - cutstudies nr 6
            if (nameCutVariationSC[i].CompareTo("Qt")==0){
                if (meson.CompareTo("Pi0")==0 ){
                    for (Int_t k = 3; k < nPtBins; k++){
                        errorReset              = 0.15+pow(ptBins[k],2)*0.024+ptBins[k]*0.02;
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = 0.01*errorReset;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = 0.01*errorReset;
                    }
                } else if ( meson.CompareTo("Eta")==0  ||  meson.CompareTo("EtaToPi0")==0  ){
                    for (Int_t k = 2; k < nPtBins; k++){
                        errorReset              = 0.8+pow(ptBins[k],2)*0.024+ptBins[k]*0.02;
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = 0.01*errorReset;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = 0.01*errorReset;
                    }
                }
            }


            // Alpha - cutstudies nr 7
            if (nameCutVariationSC[i].CompareTo("Alpha")==0){
                minPt           = 0;
                errorReset      = 0.3;
                if (meson.CompareTo("EtaToPi0")==0  )
                    errorReset      = 0.45;
            }

            // Chi2/PsiPair - cutstudies nr 8
            if (nameCutVariationSC[i].CompareTo("Chi2PsiPair")==0 ){
                if ( meson.CompareTo("Pi0")==0){
                    for (Int_t k = 2; k < nPtBins; k++){
                        errorReset              = 0.45+pow(ptBins[k],2)*0.01+ptBins[k]*0.01;
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = 0.01*errorReset;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = 0.01*errorReset;
                    }
                } else if (meson.CompareTo("Eta")==0 || meson.CompareTo("EtaToPi0")==0){
                    for (Int_t k = 2; k < nPtBins; k++){
                        errorReset              = 1.65+pow(ptBins[k],2)*0.03+ptBins[k]*0.02;
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = 0.01*errorReset;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = 0.01*errorReset;
                    }
                }
            }

            // CosPoint - cutstudies nr 9
            if (nameCutVariationSC[i].CompareTo("CosPoint")==0  ){
                minPt           = 0;
                errorReset      = 0.;
            }

            // BG - cutstudies nr 10
            if (nameCutVariationSC[i].CompareTo("BG")==0  ){
                minPt           = 0;
                errorReset      = 0.15;
                if (meson.CompareTo("EtaToPi0")==0 )
                    errorReset      = 0.25;
            }

            // MCSmearing - cutstudies nr 11
            if (nameCutVariationSC[i].CompareTo("MCSmearing")==0  ){
                if (meson.CompareTo("Pi0")==0 ){
                    minPt           = 0;
                    errorReset      = 0.5;
                } else if (meson.CompareTo("Eta")==0 ){
                    minPt           = 0;
                    errorReset      = 1.2;
                } else if (meson.CompareTo("EtaToPi0")==0 ){
                    minPt           = 0;
                    errorReset      = 0.8;
                }
            }

            // BGEstimateIterations - cutstudies nr 12


            // Eta - cutstudies nr 13
            if (nameCutVariationSC[i].CompareTo("Efficiency")==0  && energy.CompareTo("pPb_5.023TeVCent")){
                if ( meson.CompareTo("Pi0") == 0){
                    minPt           = 0;
                    errorReset      = 2.0;
                } else if ( meson.CompareTo("Eta")==0  ){
                    minPt           = 0;
                    errorReset      = 2.0;
                }
            }
            // Eta - cutstudies nr 14
            if (nameCutVariationSC[i].CompareTo("Eta")==0  ){
                if ( meson.Contains("Pi0") ){
                    minPt           = 0;
                    errorReset      = 0.25;
                } else if ( meson.CompareTo("Eta")==0  ){
                    minPt           = 0;
                    errorReset      = 3.0;
                }
            }


            // reset to a constant if enabled
            if (minPt != -10 && errorReset != -10000 ){
                for (Int_t k = 0; k < nPtBins; k++){
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }
            }
        }

        for (Int_t l = 0; l < nPtBins; l++){
            errorsPosSummed[l]          = errorsPosSummed[l]+pow(errorsPos[i][l],2);
            errorsNegSummed[l]          = errorsNegSummed[l]+ pow(errorsNeg[i][l],2);
            errorsMeanSummed[l]         = errorsMeanSummed[l]+ pow(errorsMean[i][l],2);
            errorsPosCorrSummed[l]      = errorsPosCorrSummed[l]+pow(errorsPosCorr[i][l],2);
            errorsNegCorrSummed[l]      = errorsNegCorrSummed[l] +pow(errorsNegCorr[i][l],2);
            errorsMeanCorrSummed[l]     = errorsMeanCorrSummed[l]+ pow(errorsMeanCorr[i][l],2);
        }
        negativeErrors[i]           = new TGraphErrors(nPtBins,ptBins ,errorsNeg[i] ,ptBinsErr ,errorsNegErr[i] );
        meanErrors[i]               = new TGraphErrors(nPtBins,ptBins ,errorsMean[i] ,ptBinsErr ,errorsMeanErr[i] );
        positiveErrors[i]           = new TGraphErrors(nPtBins,ptBins ,errorsPos[i] ,ptBinsErr ,errorsPosErr[i] );
        negativeErrorsCorr[i]       = new TGraphErrors(nPtBins,ptBins ,errorsNegCorr[i] ,ptBinsErr ,errorsNegErrCorr[i] );
        meanErrorsCorr[i]           = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorr[i] ,ptBinsErr ,errorsMeanErrCorr[i] );
        positiveErrorsCorr[i]       = new TGraphErrors(nPtBins,ptBins ,errorsPosCorr[i] ,ptBinsErr ,errorsPosErrCorr[i] );

    }

    Double_t errorMaterial      = 4.50;
    if (meson.CompareTo("EtaToPi0") == 0 || meson.CompareTo("Pi0Ratio") == 0)
        errorMaterial           = 0.0;

    for (Int_t l = 0; l < nPtBins; l++){
        errorsPosSummed[l]              = pow(errorsPosSummed[l],0.5);
        errorsMeanSummed[l]             = pow(errorsMeanSummed[l],0.5);
        errorsPosErrSummed[l]           = errorsPosSummed[l]*0.001;
        errorsMeanErrSummed[l]          = errorsMeanSummed[l]*0.001;
        errorsNegSummed[l]              = -pow(errorsNegSummed[l],0.5);
        errorsNegErrSummed[l]           = errorsNegSummed[l]*0.001;
        errorsPosCorrMatSummed[l]       = pow(errorsPosCorrSummed[l]+ pow(2*errorMaterial ,2.),0.5);
        errorsMeanCorrMatSummed[l]      = pow(errorsMeanCorrSummed[l]+ pow(2*errorMaterial ,2.),0.5);
        errorsNegCorrMatSummed[l]       = -pow(errorsNegCorrSummed[l]+ pow(2*errorMaterial ,2.),0.5);

        errorsPosCorrSummed[l]          = pow(errorsPosCorrSummed[l],0.5);
        errorsMeanCorrSummed[l]         = pow(errorsMeanCorrSummed[l],0.5);
        errorsPosErrCorrSummed[l]       = errorsPosCorrSummed[l]*0.001;
        errorsMeanErrCorrSummed[l]      = errorsMeanCorrSummed[l]*0.001;
        errorsMeanErrCorrMatSummed[l]   = errorsMeanCorrMatSummed[l]*0.001;
        errorsNegCorrSummed[l]          = -pow(errorsNegCorrSummed[l],0.5);
        errorsNegErrCorrSummed[l]       = errorsNegCorrSummed[l]*0.001;
// 		cout << errorsMeanSummed[l] << "\t" << errorsMeanCorrMatSummed[l]<< endl;
    }

    Double_t errorsMat[nPtBins];
    for (Int_t l = 0; l < nPtBins; l++){
        errorsMat[l]            = 2*errorMaterial;

    }
    TGraphErrors* graphMaterialError    = new TGraphErrors(nPtBins,ptBins ,errorsMat ,ptBinsErr ,errorsMeanErrSummed );
    negativeErrorsSummed                = new TGraphErrors(nPtBins,ptBins ,errorsNegSummed ,ptBinsErr ,errorsNegErrSummed );
    negativeErrorsCorrSummed            = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrSummed ,ptBinsErr ,errorsNegErrCorrSummed );
    positiveErrorsSummed                = new TGraphErrors(nPtBins,ptBins ,errorsPosSummed ,ptBinsErr ,errorsPosErrSummed );
    positiveErrorsCorrSummed            = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrSummed ,ptBinsErr ,errorsPosErrCorrSummed );
    meanErrorsSummed                    = new TGraphErrors(nPtBins,ptBins ,errorsMeanSummed ,ptBinsErr ,errorsMeanErrSummed );
    meanErrorsCorrSummed                = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSummed ,ptBinsErr ,errorsMeanErrCorrSummed );
    meanErrorsCorrSummedIncMat          = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrMatSummed ,ptBinsErr ,errorsMeanErrCorrMatSummed );


    // **************************************************************************************************************************************************
    // ********************************************** Start Plotting of systematics *********************************************************************
    // **************************************************************************************************************************************************

    // Give legend position for plotting
    Double_t minXLegend     = 0.12;
    Double_t maxYLegend     = 0.95;
    if (meson.CompareTo("Eta") == 0){
        minXLegend          = 0.16;
    }
    Double_t widthLegend    = 0.25;
    if (numberCutStudies> 7)
        widthLegend         = 0.5;
    Double_t heightLegend   = 1.05* 0.035 * (numberCutStudies+3);
    if (numberCutStudies> 7)
        heightLegend        = 1.05* 0.035 * (numberCutStudies/2+1);
    if (numberCutStudies> 7 && meson.CompareTo("Eta") == 0)
        heightLegend        = 1.05* 0.035 * (numberCutStudies/2+1);
    if (numberCutStudies> 7 && meson.CompareTo("EtaToPi0") == 0)
        heightLegend        = 1.05* 0.035 * (numberCutStudies/2);

    // ***************************************************************************************************
    // ****************************** Plot all mean erros separately *************************************
    // ***************************************************************************************************
    TCanvas* canvasSysErrMean = new TCanvas("canvasSysErrMean","",200,10,1350,900);// gives the page size
    DrawGammaCanvasSettings( canvasSysErrMean, 0.08, 0.01, 0.015, 0.09);

    // create dummy histo
    TH2D *histo2DSysErrMean ;
    if (meson.Contains("Pi0") ){
        histo2DSysErrMean = new TH2D("histo2DSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,0.,25.);
    } else {
        histo2DSysErrMean = new TH2D("histo2DSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,0.,25.);
    }
    SetStyleHistoTH2ForGraphs( histo2DSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean systematic Err %", 0.03, 0.04, 0.03, 0.04,
                               1,0.9, 510, 510);
    histo2DSysErrMean->Draw();

    // create legend
    TLegend* legendMean = GetAndSetLegend2(minXLegend,maxYLegend-heightLegend,minXLegend+widthLegend,maxYLegend, 30);
    if (numberCutStudies> 7) legendMean->SetNColumns(2);

    for(Int_t i = 0;i< numberCutStudies ;i++){
        DrawGammaSetMarkerTGraphErr(meanErrors[i], markerStyle[i], 1.,color[i],color[i]);
        meanErrors[i]->Draw("pE0,csame");
        legendMean->AddEntry(meanErrors[i],nameCutVariation[i].Data(),"p");
    }
    legendMean->Draw();

    // plot labeling
    TLatex *labelMeson;
    if (meson.CompareTo("EtaToPi0") == 0){
        labelMeson= new TLatex(0.95,0.89,Form("#eta/#pi^{0} rec. PCM"));
    } else if (meson.Contains("Pi0Ratio")){
        labelMeson= new TLatex(0.95,0.89,Form("#it{R}_{pA} PCM"));
    } else if (meson.Contains("Pi0")){
        labelMeson= new TLatex(0.95,0.89,Form("#pi^{0} PCM"));
    } else {
        labelMeson= new TLatex(0.95,0.89,Form("#eta PCM"));
    }
    SetStyleTLatex( labelMeson, 0.038,4,1,42,kTRUE, 31);
    labelMeson->Draw();

    TLatex *labelCentrality = new TLatex(0.95,0.93,Form("%s %s", additionalName.Data(), collisionSystem.Data()));
    SetStyleTLatex( labelCentrality, 0.038,4,1,42,kTRUE, 31);
    labelCentrality->Draw();

    canvasSysErrMean->Update();
    canvasSysErrMean->SaveAs(Form("SystematicErrorsCalculatedConv/SysMean_%s_%s%s_%s.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
    delete canvasSysErrMean;

    // ***************************************************************************************************
    // ********************* Plot all mean erros separately after smoothing ******************************
    // ***************************************************************************************************
    TCanvas* canvasNewSysErrMean = new TCanvas("canvasNewSysErrMean","",200,10,1350,900);// gives the page size
    DrawGammaCanvasSettings( canvasNewSysErrMean, 0.08, 0.01, 0.015, 0.09);

    // create dummy histo
    TH2D *histo2DNewSysErrMean ;
    if (meson.Contains("Pi0")){
        histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,25.);
    } else {
        histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,25.);
    }
    SetStyleHistoTH2ForGraphs( histo2DNewSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %", 0.03, 0.04, 0.03, 0.04,
                               1,0.9, 510, 510);
    histo2DNewSysErrMean->Draw();

    // create legend
    TLegend* legendMeanNew = GetAndSetLegend2(minXLegend,maxYLegend-heightLegend,minXLegend+widthLegend,maxYLegend, 30);
    legendMeanNew->SetMargin(0.1);
    if (numberCutStudies> 7) legendMeanNew->SetNColumns(2);

    for(Int_t i = 0;i< numberCutStudies ;i++){
        cout << i << "\t"<< additionalNameOutput.Data() << endl;
        DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], markerStyle[i], 1.,color[i],color[i]);
        meanErrorsCorr[i]->Draw("pX0,csame");
        legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");
    }
    if (meson.CompareTo("EtaToPi0") && meson.CompareTo("Pi0Ratio")){
        DrawGammaSetMarkerTGraphErr(graphMaterialError, GetMarkerStyleSystematics( "Material"), 1., GetColorSystematics( "Material" ), GetColorSystematics( "Material"));
        graphMaterialError->Draw("p,csame");
        legendMeanNew->AddEntry(graphMaterialError,"Material","p");
    }
    DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
    meanErrorsCorrSummedIncMat->Draw("p,csame");
    legendMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
    legendMeanNew->Draw();

    // labeling
    labelMeson->Draw();
    labelCentrality->Draw();

    canvasNewSysErrMean->Update();
    canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedConv/SysMeanNewWithMean_%s_%s%s_%s.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));



    // ***************************************************************************************************
    // ********************* Plot unsmoothed errors with fits ********************************************
    // ***************************************************************************************************
    for (Int_t cut =0 ;cut < numberCutStudies;cut++ ){

        canvasNewSysErrMean->cd();
        histo2DNewSysErrMean->GetYaxis()->SetRangeUser(-0.5,15);
        histo2DNewSysErrMean->Draw();


        if (!bsmooth[cut]){ ;
            cout <<endl << endl<<  "variation: " << cut << " \t"<< nameCutVariation[cut].Data() << endl;
            Double_t minPt = startPtSys;
            Double_t maxPt = ptBins[nPtBins-2]+1;

            TF1* pol0 = new TF1("pol0","[0]",minPt,maxPt);//
            TF1* pol1 = new TF1("pol1","[0]+[1]*x",minPt,maxPt);//
            TF1* pol2 = new TF1("pol2","[0]+[1]*x+[2]*x*x",minPt,maxPt);//
            TF1* pol4 = new TF1("pol4","[0]+[1]*x+[2]*x*x+[3]*x*x*x*x",minPt,maxPt);//
            TF1* bla  = new TF1("bla","[0]+[1]/pow([2],x)",minPt,maxPt);
            //             TF1* bla  = new TF1("bla","1.5+50/pow(5,x)",minPt,maxPt);
            bla->SetParameter(0,1.5);
            bla->SetParameter(1,50);
            bla->SetParameter(2,5);
            pol4->SetParLimits(3,0,10);

            meanErrorsCorr[cut]->Fit(pol4,"NRMEX0+","",minPt,maxPt);
            meanErrorsCorr[cut]->Fit(pol2,"NRMEX0+","",minPt,maxPt);
            meanErrorsCorr[cut]->Fit(pol1,"NRMEX0+","",minPt,maxPt);
            meanErrorsCorr[cut]->Fit(pol0,"NRMEX0+","",minPt,maxPt);
            meanErrorsCorr[cut]->Fit(bla,"NRMEX0+","",minPt,maxPt);
            pol4->SetLineColor(kRed+2);
            pol2->SetLineColor(kBlue+2);
            pol1->SetLineColor(kCyan+2);
            pol0->SetLineColor(kBlack);
            bla->SetLineColor(kMagenta+2);

            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[cut], markerStyle[cut], 1.,color[cut],color[cut]);
            meanErrorsCorr[cut]->Draw("p,csame");
            pol4->Draw("same");
            pol2->Draw("same");
            pol1->Draw("same");
            pol0->Draw("same");
            bla->Draw("same");
            cout <<"plotting: " << color[cut] << "\t" << markerStyle[cut] <<endl;

            meanErrorsCorr[cut]->Print();
            canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedConv/SysMeanNewWithMeanSingle_%s_%s%s_%s_Variation%d.%s",meson.Data(), energyForOutput.Data(), additionalNameOutput.Data(), dateForOutput.Data(),cut, suffix.Data()));
        } else {
            cout <<endl << endl<<  "variation: " << cut << " \t"<< nameCutVariation[cut].Data() << endl;

            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[cut], markerStyle[cut], 1.,color[cut],color[cut]);
            meanErrorsCorr[cut]->Draw("p,csame");

            canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedConv/SysMeanNewWithMeanSingleSmoothed_%s_%s%s_%s_Variation%d.%s",meson.Data(), energyForOutput.Data(), additionalNameOutput.Data(), dateForOutput.Data(),cut, suffix.Data()));
        }
    }

    // ***************************************************************************************************
    // ********************* Create output files with errors *********************************************
    // ***************************************************************************************************
    const char *SysErrDatname = Form("SystematicErrorsCalculatedConv/SystematicErrorPCM_%s_%s%s_%s.dat", meson.Data(), energyForOutput.Data(), additionalNameOutput.Data(), dateForOutput.Data());
    fstream SysErrDat;
    cout << SysErrDatname << endl;
    SysErrDat.open(SysErrDatname, ios::out);
    for (Int_t l=0; l< nPtBins-1; l++){
        SysErrDat << ptBins[l] <<"\t" << errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
    }
    SysErrDat.close();

    const char *SysErrDatnameMean = Form("SystematicErrorsCalculatedConv/SystematicErrorAveragedPCM_%s_%s%s_%s.dat", meson.Data(), energyForOutput.Data(), additionalNameOutput.Data(), dateForOutput.Data());
    fstream SysErrDatAver;
    cout << SysErrDatnameMean << endl;
    SysErrDatAver.open(SysErrDatnameMean, ios::out);
    for (Int_t l=0; l< nPtBins; l++){
        SysErrDatAver  << ptBins[l] <<"\t"<< "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
        // 			SysErrDatAver << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
    }
    SysErrDatAver.close();

    const char *SysErrDatnameMeanSingleErr = Form("SystematicErrorsCalculatedConv/SystematicErrorAveragedSinglePCM_%s_%s%s_%s.dat", meson.Data(), energyForOutput.Data(), additionalNameOutput.Data(), dateForOutput.Data());
    fstream SysErrDatAverSingle;
    cout << SysErrDatnameMeanSingleErr << endl;
    SysErrDatAverSingle.open(SysErrDatnameMeanSingleErr, ios::out);
    SysErrDatAverSingle << "Pt bin\t" ;
    for (Int_t i= 0; i< numberCutStudies; i++){
//         if(!meson.CompareTo("EtaToPi0")&&i==1)
//             continue;
        SysErrDatAverSingle << nameCutVariationSC[i] << "\t";
    }
    if(meson.CompareTo("EtaToPi0")  &&  meson.CompareTo("Pi0Ratio"))
        SysErrDatAverSingle << "Material";
    SysErrDatAverSingle << endl;
    for (Int_t l=0;l< nPtBins;l++){
        SysErrDatAverSingle << ptBins[l] << "\t";
        for (Int_t i= 0; i< numberCutStudies; i++){
//             if(!meson.CompareTo("EtaToPi0")&&i==1)
//                 continue;
            SysErrDatAverSingle << errorsMeanCorr[i][l] << "\t";
        }
        if(meson.CompareTo("EtaToPi0") &&  meson.CompareTo("Pi0Ratio"))
            SysErrDatAverSingle << 9 << "\t";
        SysErrDatAverSingle << errorsMeanCorrMatSummed[l] << endl;
    }
    SysErrDatAverSingle.close();


    // ***************************************************************************************************
    // ********************* Group errors according to topic *********************************************
    // ***************************************************************************************************
    Double_t errorsMeanCorrPID[nPtBins];
    Double_t errorsMeanCorrSignalExtraction[nPtBins];
    Double_t errorsMeanCorrTrackReco[nPtBins];
    Double_t errorsMeanCorrPhotonReco[nPtBins];
    Double_t errorsMeanCorrPileup[nPtBins];

    for (Int_t l=0; l< nPtBins; l++){

        //0 - YieldExtraction
        //1 - BGEstimate & YieldExtractionPi0 for eta/pi0
        //2 - dEdxE
        //3 - dEdxPi
        //4 - TPCCluster
        //5 - SinglePt
        //6 - Qt
        //7 - Alpha
        //8 - Chi2PsiPair
        //9 - CosPoint
        //10 - BG
        //11 - MCSmearing
        //12 - BGEstimateIterations
        //13 - Eta

        // Signal extraction: Yield extraction, Alpha ,BG, MC Smearing
        errorsMeanCorrSignalExtraction[l]       =   TMath::Sqrt(pow(errorsMeanCorr[0][l],2)+    // Yield extraction
                                                                pow(errorsMeanCorr[7][l],2)+    // Alpha
                                                                pow(errorsMeanCorr[9][l],2)+    // BG
                                                                pow(errorsMeanCorr[11][l],2));  // MCSmearing
        if (!meson.CompareTo("EtaToPi0") || !meson.CompareTo("Pi0Ratio")){
            errorsMeanCorrSignalExtraction[l]   =   TMath::Sqrt(pow(errorsMeanCorr[0][l],2)+    // Yield extraction eta
                                                                pow(errorsMeanCorr[7][l],2)+    // Alpha
                                                                pow(errorsMeanCorr[9][l],2)+    // BG
                                                                pow(errorsMeanCorr[1][l],2)+    // Yield extraction pi0
                                                                pow(errorsMeanCorr[11][l],2)  );// MCSmearing

        }
        // PID: dEdxE, dEdxPi
        errorsMeanCorrPID[l]                    =   TMath::Sqrt(pow(errorsMeanCorr[2][l],2)+    // dEdxE
                                                                pow(errorsMeanCorr[3][l],2));   // dEdxPi

        // photon reco: Chi2, Qt, PsiPair, CosPoint, MinR
        errorsMeanCorrPhotonReco[l]             =   TMath::Sqrt(pow(errorsMeanCorr[6][l],2)+    // Chi2PsiPair
                                                                pow(errorsMeanCorr[8][l],2));   // Qt

        // track reconstruction: TPCCluster, Single pT, Eta
        errorsMeanCorrTrackReco[l]              =   TMath::Sqrt(pow(errorsMeanCorr[4][l],2)+    // TPCCluster
                                                                pow(errorsMeanCorr[5][l],2));   // Single pT
        // pileup
        if (numberCutStudies > 12){
            errorsMeanCorrPileup[l]                 =   TMath::Sqrt(pow(errorsMeanCorr[1][l],2)+    // out of bunch methods
            pow(errorsMeanCorr[12][l],2));  // out of bunch iterations
        } else {
            errorsMeanCorrPileup[l]                 =   TMath::Sqrt(pow(errorsMeanCorr[1][l],2)); // out of bunch methods
        }

    }
    TGraphErrors* meanErrorsPID                 = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPID ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsPhotonReco          = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPhotonReco ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsSignalExtraction    = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSignalExtraction ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsTrackReco           = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrTrackReco ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsPileup              = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPileup ,ptBinsErr ,errorsMeanErrCorrSummed );

    // ***************************************************************************************************
    // ********************* Plot grouped errors for better understanding ********************************
    // ***************************************************************************************************
    Double_t minXLegend2 = 0.13;
    Double_t maxYLegend2 = 0.95;
    Double_t widthLegend2 = 0.52;
    Double_t heightLegend2 = 0.15;

    TCanvas* canvasSummedErrMean = new TCanvas("canvasSummedErrMean","",200,10,1350,900);// gives the page size
    DrawGammaCanvasSettings( canvasSummedErrMean, 0.08, 0.01, 0.015, 0.09);

    // create dummy histo
    TH2D *histo2DSummedErrMean ;
    if (meson.Contains("Pi0") ){
        histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,25.);
    } else {
        histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,25.);
    }
    SetStyleHistoTH2ForGraphs( histo2DSummedErrMean, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %", 0.03, 0.04, 0.03, 0.04,
                               1,0.9, 510, 510);
    histo2DSummedErrMean->Draw();

    // create legend
    TLegend* legendSummedMeanNew = GetAndSetLegend2(minXLegend2,maxYLegend2-heightLegend2,minXLegend2+widthLegend2,maxYLegend2, 30);
    legendSummedMeanNew->SetNColumns(2);
    legendSummedMeanNew->SetMargin(0.1);

    // Signal extraction error
    DrawGammaSetMarkerTGraphErr(meanErrorsSignalExtraction, 20, 1.,color[0],color[0]);
    meanErrorsSignalExtraction->Draw("p,csame");
    legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"Signal Extraction","p");
    cout << "here" << endl;
    DrawGammaSetMarkerTGraphErr(meanErrorsPID, markerStyle[2], 1.,color[2],color[2]);
    meanErrorsPID->Draw("p,csame");
    legendSummedMeanNew->AddEntry(meanErrorsPID,"Electron PID","p");
    cout << "here" << endl;
    DrawGammaSetMarkerTGraphErr(meanErrorsTrackReco, markerStyle[5], 1.,color[5],color[5]);
    meanErrorsTrackReco->Draw("p,csame");
    legendSummedMeanNew->AddEntry(meanErrorsTrackReco,"Track Reconstruction","p");
    cout << "here" << endl;
    DrawGammaSetMarkerTGraphErr(meanErrorsPhotonReco, markerStyle[8], 1.,color[6],color[6]);
    meanErrorsPhotonReco->Draw("p,csame");
    legendSummedMeanNew->AddEntry(meanErrorsPhotonReco,"Photon Reconstruction","p");
    cout << "here" << endl;
    if (meson.CompareTo("EtaToPi0") &&  meson.CompareTo("Pi0Ratio")){
        DrawGammaSetMarkerTGraphErr(meanErrorsPileup, GetMarkerStyleSystematics( "BGEstimate"), 1.,GetColorSystematics( "BGEstimate"),GetColorSystematics( "BGEstimate"));
        meanErrorsPileup->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsPileup,"Pileup Estimate","p");
        cout << "here" << endl;
        DrawGammaSetMarkerTGraphErr(graphMaterialError, GetMarkerStyleSystematics( "InnerMaterial"), 1.,GetColorSystematics( "InnerMaterial"),GetColorSystematics( "InnerMaterial"));
        graphMaterialError->Draw("p,csame");
        legendSummedMeanNew->AddEntry(graphMaterialError,"Material","p");
    }

    DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
    meanErrorsCorrSummedIncMat->Draw("p,csame");
    legendSummedMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
    legendSummedMeanNew->Draw();

    labelMeson->Draw();
    labelCentrality->Draw();

    canvasSummedErrMean->Update();
    canvasSummedErrMean->SaveAs(Form("SystematicErrorsCalculatedConv/SysErrorSummedVisu_%s_%s%s_%s.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));

    delete canvasSummedErrMean;

    const char *SysErrDatnameMeanPaper = Form("SystematicErrorsCalculatedConv/SystematicErrorAveragedPCMPaper_%s_%s%s_%s.dat",meson.Data(),energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data());
    fstream SysErrDatAverPaper;
    cout << SysErrDatnameMeanPaper << endl;
    SysErrDatAverPaper.open(SysErrDatnameMeanPaper, ios::out);
    SysErrDatAverPaper  << "p_{T}" << "\t Material \t Yield Extraction \t PID \t photon reco \t track recon \t summed" <<  endl;
    for (Int_t l=0; l< nPtBins; l++){
        SysErrDatAverPaper << ptBins[l] <<"\t" << errorMaterial*2 << "\t" << errorsMeanCorrSignalExtraction[l] << "\t" << errorsMeanCorrPID[l]<< "\t" << errorsMeanCorrPhotonReco[l]<< "\t" <<errorsMeanCorrTrackReco[l] <<"\t" << errorsMeanCorrPileup[l] << "\t" << errorsMeanCorrMatSummed[l]<< endl;
    }
    SysErrDatAverPaper.close();
}