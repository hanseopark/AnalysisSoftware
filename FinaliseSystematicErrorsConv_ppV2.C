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

void FinaliseSystematicErrorsConv_ppV2( TString nameDataFileErrors      = "",
                                        TString energy                  = "",
                                        TString meson                   = "",
                                        Int_t numberOfPtBins            = 1 ,
                                        Int_t numberCutStudies          = 1,
                                        Double_t startPtSys             = 0,
                                        TString additionalName          = "",
                                        TString additionalNameOutput    = "",
                                        TString suffix                  = "eps",
                                        Bool_t showBeforeSmoothing      = kFALSE
                                      ){

    // ***************************************************************************************************
    // ****************************** General style settings *********************************************
    // ***************************************************************************************************

    StyleSettingsThesis();
    SetPlotStyle();

    // ***************************************************************************************************
    // ****************************** Create output directory ********************************************
    // ***************************************************************************************************
    gSystem->Exec("mkdir -p SystematicErrorsCalculatedConv");

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

    TString nameCutVariationSCCurrent[16]   = { "YieldExtraction_pp", "BGEstimate_pp", "dEdxE", "dEdxPi", "TPCCluster",
                                                "SinglePt", "Chi2" , "Qt" , "Alpha" , "BG" ,
                                                "Periods" , "SPD" , "CosPoint", "PsiPair" ,"Eta" ,
                                                "MCSmearing" };
    if(!energy.CompareTo("900GeV"))
        nameCutVariationSCCurrent[10] = "SPD";
    if(!energy.CompareTo("7TeV"))
        nameCutVariationSCCurrent[10] = "BGEstimateIterations_pp";
    if(!energy.CompareTo("8TeV"))
        nameCutVariationSCCurrent[12] = "PileupDCA";
    if(!energy.CompareTo("5TeV2017")){
        nameCutVariationSCCurrent[9] = "BGEstimateIterations_pp";
        nameCutVariationSCCurrent[10] = "DoubleCount";
        if (meson.CompareTo("EtaToPi0"))
            nameCutVariationSCCurrent[11] = "CosPoint";

    }
    Color_t color[20];
    Color_t markerStyle[20];
    for (Int_t k =0; k<16; k++ ){
        cout << "variation: " << nameCutVariationSCCurrent[k].Data() << endl;
        color[k]                            = GetColorSystematics( nameCutVariationSCCurrent[k] );
        markerStyle[k]                      = GetMarkerStyleSystematics( nameCutVariationSCCurrent[k] );
        nameCutVariation[k]                 = GetSystematicsName(nameCutVariationSCCurrent[k]);
        cout << "name for writing: " << nameCutVariation[k].Data() << endl;
    }

    for (Int_t i = 0; i < numberCutStudies; i++){
        nameCutVariationSC[i]               = nameCutVariationSCCurrent[i];
    }

    if (meson.CompareTo("EtaToPi0") == 0){
        nameCutVariation[0]                 = "Yield extraction #eta";
        nameCutVariation[1]                 = "Yield extraction #pi^{0}";
        nameCutVariationSC[1]               = "YieldExtraction_pp";
        color[1]                            = GetColorSystematics( "YieldExtractionPi0" );
        markerStyle[1]                      = GetMarkerStyleSystematics( "YieldExtractionPi0" );

        if (energy.CompareTo("2.76TeV") == 0){
            nameCutVariation[10]            = "pile-up";
            nameCutVariationSC[10]          = "BGEstimate_pp";
            color[10]                       = GetColorSystematics( nameCutVariationSC[10] );
            markerStyle[10]                 = GetMarkerStyleSystematics( nameCutVariationSC[10] );
        }
        if (energy.CompareTo("5TeV2017") == 0){
            nameCutVariation[11]            = "pile-up";
            nameCutVariationSC[11]          = "BGEstimate_pp";
            color[11]                       = GetColorSystematics( nameCutVariationSC[10] );
            markerStyle[11]                 = GetMarkerStyleSystematics( nameCutVariationSC[10] );
        }
    }

    // ***************************************************************************************************
    // ******************************** Booleans for smoothing *******************************************
    // ***************************************************************************************************
    Bool_t bsmooth[16]                      = { 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothMBPi07TeV[16]             = { 1, 0, 1, 1, 1,  1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,  0 };
    Bool_t bsmoothMBEta7TeV[16]             = { 1, 0, 1, 1, 1,  1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,  0 };
    Bool_t bsmoothMBEtaToPi07TeV[16]        = { 1, 1, 1, 1, 1,  1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,  0 };
    Bool_t bsmoothMBPi08TeV[16]             = { 1, 0, 1, 1, 1,  1, 1, 1, 1, 1,
                                                1, 1, 0, 0, 0,  0 };
    Bool_t bsmoothMBEta8TeV[16]             = { 1, 1, 1, 1, 1,  1, 1, 1, 1, 1,
                                                1, 1, 0, 0, 0,  0 };
    Bool_t bsmoothMBEtaToPi08TeV[16]        = { 1, 1, 1, 1, 1,  1, 1, 1, 1, 1,
                                                1, 1, 0, 0, 0,  0 };
    Bool_t bsmoothMBPi02760GeV[16]          = { 1, 0, 1, 1, 1,  1, 1, 1, 1, 1,
                                                0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothMBEta2760GeV[16]          = { 1, 0, 1, 1, 1,  1, 1, 1, 1, 1,
                                                0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothMBEtaToPi02760GeV[16]     = { 1, 0, 1, 1, 1,  1, 1, 1, 1, 1,
                                                0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothMBPi0900GeV[16]           = { 1, 1, 1, 0, 0,  1, 1, 1, 1, 0,
                                                0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothMBEta900GeV[16]           = { 1, 1, 0, 0, 0,  1, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothMBEtaToPi0900GeV[16]      = { 0, 0, 0, 0, 0,  1, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothMBPi05TeV2017[16]         = { 1, 0, 1, 1, 1,  1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,  0 };
    Bool_t bsmoothMBEta5TeV2017[16]         = { 1, 0, 1, 1, 1,  1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,  0 };
    Bool_t bsmoothMBEtaToPi05TeV2017[16]    = { 1, 0, 1, 1, 1,  1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,  0 };

    for (Int_t i = 0; i < numberCutStudies; i++){
        if (energy.CompareTo("900GeV") == 0){
            if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Pi0")==0){
                bsmooth[i]                      = bsmoothMBPi0900GeV[i];
            } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Eta")==0){
                bsmooth[i]                      = bsmoothMBEta900GeV[i];
            } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("EtaToPi0")==0){
                bsmooth[i]                      = bsmoothMBEtaToPi0900GeV[i];
            }
        } else if (energy.CompareTo("2.76TeV") == 0){
            if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Pi0")==0){
                bsmooth[i]                      = bsmoothMBPi02760GeV[i];
            } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Eta")==0){
                bsmooth[i]                      = bsmoothMBEta2760GeV[i];
            } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("EtaToPi0")==0){
                bsmooth[i]                      = bsmoothMBEtaToPi02760GeV[i];
            }
        } else if (energy.CompareTo("7TeV") == 0){
            if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Pi0")==0){
                bsmooth[i]                      = bsmoothMBPi07TeV[i];
            } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Eta")==0){
                bsmooth[i]                      = bsmoothMBEta7TeV[i];
            } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("EtaToPi0")==0){
                bsmooth[i]                      = bsmoothMBEtaToPi07TeV[i];
            }
        } else if (energy.CompareTo("8TeV") == 0){
            if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Pi0")==0){
                bsmooth[i]                      = bsmoothMBPi08TeV[i];
            } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Eta")==0){
                bsmooth[i]                      = bsmoothMBEta8TeV[i];
            } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("EtaToPi0")==0){
                bsmooth[i]                      = bsmoothMBEtaToPi08TeV[i];
            }
        } else if (energy.CompareTo("5TeV2017") == 0){
            if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Pi0")==0){
                bsmooth[i]                      = bsmoothMBPi05TeV2017[i];
            } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Eta")==0){
                bsmooth[i]                      = bsmoothMBEta5TeV2017[i];
            } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("EtaToPi0")==0){
                bsmooth[i]                      = bsmoothMBEtaToPi05TeV2017[i];
            }
        }
    }

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
    Double_t errorsMeanCorrBefore[nCuts][nPtBins];
    Double_t errorsMeanSummed[nPtBins];
    Double_t errorsMeanCorrSummed[nPtBins];
    Double_t errorsMeanCorrMatSummed[nPtBins];

    Double_t errorsMeanErr[nCuts][nPtBins];
    Double_t errorsMeanErrCorr[nCuts][nPtBins];
    Double_t errorsMeanErrCorrBefore[nCuts][nPtBins];
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
    TGraphErrors* meanErrorsCorrBefore[nCuts];
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
        if ( nameCutVariationSC[i].CompareTo("BGEstimate_pp")==0  ){
            TString nameGraphPos            = Form("%s_SystErrorRel_%s",meson.Data(),nameCutVariationSC[i].Data()  );
            TString nameGraphNeg            = Form("%s_SystErrorRel_%s",meson.Data(),nameCutVariationSC[i].Data()  );
            if (meson.CompareTo("EtaToPi0")==0){
                nameGraphPos                = Form("Eta_SystErrorRel_%s",nameCutVariationSC[i].Data()  );
                nameGraphNeg                = Form("Eta_SystErrorRel_%s",nameCutVariationSC[i].Data()  );
            }
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else if ( nameCutVariationSC[i].CompareTo("BGEstimateIterations_pp")==0  ){
            TString nameGraphPos            = Form("%s_SystErrorRel_%s",meson.Data(),nameCutVariationSC[i].Data()  );
            TString nameGraphNeg            = Form("%s_SystErrorRel_%s",meson.Data(),nameCutVariationSC[i].Data()  );
            if (meson.CompareTo("EtaToPi0")==0){
                nameGraphPos                = Form("Eta_SystErrorRel_%s",nameCutVariationSC[i].Data()  );
                nameGraphNeg                = Form("Eta_SystErrorRel_%s",nameCutVariationSC[i].Data()  );
            }
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else if ( nameCutVariationSC[i].CompareTo("YieldExtraction_pp")==0  && (!meson.CompareTo("Pi0") || !meson.CompareTo("Eta") )){
          TString nameGraphPos            = Form("%s_SystErrorRelPos_%s",meson.Data(),nameCutVariationSC[i].Data() );
          TString nameGraphNeg            = Form("%s_SystErrorRelNeg_%s",meson.Data(),nameCutVariationSC[i].Data() );
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());

        } else if ( nameCutVariationSC[i].CompareTo("YieldExtraction_pp")==0 && !meson.CompareTo("EtaToPi0") &&  i == 0){

            TString nameGraphPos            = Form("Eta_SystErrorRelPos_%s",nameCutVariationSC[i].Data()  );
            TString nameGraphNeg            = Form("Eta_SystErrorRelNeg_%s",nameCutVariationSC[i].Data()  );
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else if ( nameCutVariationSC[i].CompareTo("YieldExtraction_pp")==0 && !meson.CompareTo("EtaToPi0") &&  i == 1){

            TString nameGraphPos            = Form("Pi0EtaBinning_SystErrorRelPos_%s",nameCutVariationSC[i].Data()  );
            TString nameGraphNeg            = Form("Pi0EtaBinning_SystErrorRelNeg_%s",nameCutVariationSC[i].Data()  );
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
//         } else if ( nameCutVariationSC[i].CompareTo("Periods")==0  ){
//             cout << "including period uncertainty for 7 and 8 TeV pp" << endl;
//             TString nameGraphPos            = Form("%s_SystErrorRelPos_%s",meson.Data(),nameCutVariationSC[i-1].Data()  );
//             TString nameGraphNeg            = Form("%s_SystErrorRelNeg_%s",meson.Data(),nameCutVariationSC[i-1].Data()  );
//             cout << "Cutstudies " << i << "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
//             graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
//             graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else {
            TString nameGraphPos            = Form("%s_SystErrorRelPos_%spp",meson.Data(),nameCutVariationSC[i].Data() );
            TString nameGraphNeg            = Form("%s_SystErrorRelNeg_%spp",meson.Data(),nameCutVariationSC[i].Data() );
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
            if ( graphPosErrors == NULL ){
                cout << "systematic wasn't contained, setting it to 0" << endl;
                TString nameGraphPos            = Form("%s_SystErrorRelPos_%s",meson.Data(),nameCutVariationSC[0].Data() );
                TString nameGraphNeg            = Form("%s_SystErrorRelNeg_%s",meson.Data(),nameCutVariationSC[0].Data() );
                if (meson.CompareTo("EtaToPi0") == 0){
                    nameGraphPos                = Form("Eta_SystErrorRelPos_%s",nameCutVariationSC[0].Data() );
                    nameGraphNeg                = Form("Eta_SystErrorRelNeg_%s",nameCutVariationSC[0].Data() );
                }
                cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
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
        while (graphPosErrors->GetX()[0] < startPtSys){
            graphPosErrors->RemovePoint(0);
            graphNegErrors->RemovePoint(0);
        }


        if (i == 0) {
            ptBins                          = graphNegErrors->GetX();
            ptBinsErr                       = graphNegErrors->GetEXhigh();
        }
        errorsNeg[i]                        = graphNegErrors->GetY();
        errorsNegErr[i]                     = graphNegErrors->GetEYhigh();
        errorsPos[i]                        = graphPosErrors->GetY();
        errorsPosErr[i]                     = graphPosErrors->GetEYhigh();

        cout << nameCutVariationSC[i].Data() << endl;
        CalculateMeanSysErr(            errorsMean[i],  errorsMeanErr[i],   errorsPos[i],       errorsNeg[i],           nPtBins);
        CorrectSystematicErrorsWithMean(errorsPos[i],   errorsPosErr[i],    errorsPosCorr[i],   errorsPosErrCorr[i],    nPtBins);
        CorrectSystematicErrorsWithMean(errorsNeg[i],   errorsNegErr[i],    errorsNegCorr[i],   errorsNegErrCorr[i],    nPtBins);
        CorrectSystematicErrorsWithMean(errorsMean[i],  errorsMeanErr[i],   errorsMeanCorr[i],  errorsMeanErrCorr[i],   nPtBins);
        CorrectSystematicErrorsWithMean(errorsMean[i],  errorsMeanErr[i],   errorsMeanCorrBefore[i],  errorsMeanErrCorrBefore[i],   nPtBins);


        // Routine for manual smoothing of systematic errors
        // ATTTENTION! you have to do this manually for each data set/trigger never trust the values mentioned here

        if (bsmooth[i]){
            Double_t minPt      = -10;
            Double_t errorReset = -10000;
            showBeforeSmoothing = kTRUE;

            // period variation
            if (nameCutVariationSC[i].Contains("Periods") || nameCutVariationSC[i].Contains("periods")){
                minPt       = 0;
                if (!energy.CompareTo("8TeV"))
                    errorReset  = 3.0;
                else if (!energy.CompareTo("7TeV"))
                    errorReset  = 1.5;
            }

            // pileup from out of bunch
            if (nameCutVariationSC[i].CompareTo("BGEstimate_pp")==0 ){
                minPt       = 0;
                if (!energy.CompareTo("8TeV"))
                    errorReset  = 0.7;
                else if (!energy.CompareTo("7TeV"))
                    errorReset  = 0.3;
                else if (!energy.CompareTo("900GeV"))
                    errorReset  = 0.1;
            }

            // SPD pileup
            if (nameCutVariationSC[i].CompareTo("SPD")==0){
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("8TeV")){
                        errorReset              = 0.5+pow(ptBins[k],2)*0.01;
                    } else if (!energy.CompareTo("7TeV")){
                        errorReset              = 1.0;
                    }
                    errorsMean[i][k]            = errorReset;
                    errorsMeanErr[i][k]         = errorReset*0.01;
                    errorsMeanCorr[i][k]        = errorReset;
                    errorsMeanErrCorr[i][k]     = errorReset*0.01;
                }
            }

            // BG variations
            if (nameCutVariationSC[i].CompareTo("BG")==0){
                for (Int_t k = 0; k < nPtBins; k++){
                  if (!energy.CompareTo("8TeV"))
                      errorReset  = 0.3;
                  else if (!energy.CompareTo("7TeV"))
                      errorReset  = 0.2;
                      if(meson.Contains("Eta"))
                        errorReset = 0.0471559+50.2245/pow(96.9644,ptBins[k]);
                  else if (!energy.CompareTo("2.76TeV"))
                      errorReset  = 0.2;
                  errorsMean[i][k]        = errorReset;
                  errorsMeanErr[i][k]     = errorReset*0.01;
                  errorsMeanCorr[i][k]    = errorReset;
                  errorsMeanErrCorr[i][k] = errorReset*0.01;
              }
            }

            // TPCCluster
            if (nameCutVariationSC[i].CompareTo("TPCCluster")==0 ){
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("2.76TeV"))
                        errorReset  = 0.5;
                    if (!energy.CompareTo("7TeV"))
                        errorReset  = 0.2;
                    if (!energy.CompareTo("8TeV")){
                        if (!meson.CompareTo("Pi0"))
                            errorReset  = 0.05;
                        else if (!meson.CompareTo("Eta")||!meson.CompareTo("EtaToPi0")){
                            if (ptBins[k] < 1.3)
                                errorReset      = 0.04+pow(ptBins[k]-1.3,2)*4;
                            else
                                errorReset      = 0.04;
                        }
                    }
                    if (!energy.CompareTo("5TeV2017")){
                        if (!meson.CompareTo("Pi0"))
                            errorReset  = 0.05;
                        else if (!meson.CompareTo("Eta")||!meson.CompareTo("EtaToPi0")){
                            errorReset      = 0.05;
                        }
                    }
                    errorsMean[i][k]        = errorReset;
                    errorsMeanErr[i][k]     = errorReset*0.01;
                    errorsMeanCorr[i][k]    = errorReset;
                    errorsMeanErrCorr[i][k] = errorReset*0.01;
                }
            }

            // dEdx e variations
            if (nameCutVariationSC[i].CompareTo("dEdxE")==0 ){
                if (!meson.CompareTo("Pi0")){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (!energy.CompareTo("8TeV")){
                            errorReset          = 1.0+pow(ptBins[k]-5,2)*0.065;
                        } else if (!energy.CompareTo("7TeV")){
                            if (ptBins[k] < 2.0)
                                errorReset      = 0.1+pow(ptBins[k]-2,2)*1.1;
                            else
                                errorReset      = 0.119514-0.0542626*ptBins[k]+0.0126631*pow(ptBins[k],2);
                        } else if (!energy.CompareTo("2.76TeV")){
                            errorReset          = 0.9+pow(ptBins[k],2)*0.005;
                        } else if (!energy.CompareTo("900GeV")){
                            errorReset          = 2.-1.9*ptBins[k]+0.45*pow(ptBins[k],2);
                        } else if (!energy.CompareTo("5TeV2017")){
                            if (ptBins[k] < 4.0)
                                errorReset      = 0.5;
                            else
                            errorReset          = 0.5+pow(ptBins[k]-5,2)*0.065;
                        }
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = errorReset*0.01;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = errorReset*0.01;
                    }
                } else if (meson.Contains("Eta")){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (!energy.CompareTo("8TeV")){
                            errorReset         = 1.2+pow(ptBins[k]-4,2)*0.265;
                        } else if (!energy.CompareTo("7TeV")){
                            if (ptBins[k] < 2.0)
                                errorReset          = 0.5+pow(ptBins[k]-1.5,2)*2.3;
                            else
                                errorReset          = 0.5;
                        } else if (!energy.CompareTo("2.76TeV")){
                            errorReset          = 2*(0.9+pow(ptBins[k],2)*0.005);
                        } else if (!energy.CompareTo("5TeV2017")){
                            errorReset          = 2*(0.55+pow(ptBins[k],2)*0.005);
                        }
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = errorReset*0.01;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = errorReset*0.01;
                    }
                }
            }

            // dEdx Pi variation
            if (nameCutVariationSC[i].CompareTo("dEdxPi")==0 ){
                if (!meson.CompareTo("Pi0") ){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (!energy.CompareTo("8TeV"))
                            errorReset          = 0.15+pow(ptBins[k],2)*0.01;
                        else if (!energy.CompareTo("7TeV"))
                            errorReset          = 0.0126348+-0.024018*ptBins[k]+0.0157711*pow(ptBins[k],2)+5.02481e-06*pow(ptBins[k],4);
                        else if (!energy.CompareTo("2.76TeV"))
                            errorReset          = 0.95+pow(ptBins[k],2)*0.016;
                        else if (!energy.CompareTo("5TeV2017"))
                            errorReset          = 0.08+0.14*ptBins[k]+0.023*pow(ptBins[k],2);
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = errorReset*0.01;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = errorReset*0.01;
                    }
                } else if (meson.Contains("Eta")){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (!energy.CompareTo("8TeV"))
                            errorReset          = 1.3;
                        else if (!energy.CompareTo("7TeV"))
                            errorReset          = 0.0565787+-0.0274106*ptBins[k]+0.0143598*pow(ptBins[k],2);
                        else if (!energy.CompareTo("2.76TeV"))
                            errorReset          = 0.95+pow(ptBins[k],2)*0.046;
                        else if (!energy.CompareTo("5TeV2017"))
                            errorReset          = 1.2*(1.78039+-0.35*ptBins[k]+0.048*pow(ptBins[k],2));
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = errorReset*0.01;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = errorReset*0.01;
                    }
                }
            }

            // qt variations
            if (nameCutVariationSC[i].CompareTo("Qt")==0 ){
                if ( !meson.CompareTo("Pi0") ){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (!energy.CompareTo("8TeV"))
                            errorReset          = 0.5+pow(ptBins[k],2)*0.028;
                        else if (!energy.CompareTo("7TeV"))
                            errorReset          = 0.2+pow(ptBins[k],2)*0.028;
                        else if (!energy.CompareTo("2.76TeV"))
                            errorReset          = 0.8+pow(ptBins[k],2)*0.03;
                        else if (!energy.CompareTo("900GeV"))
                            errorReset          = 1.72-2.28*ptBins[k]+.859*pow(ptBins[k],2);
                        else if (!energy.CompareTo("5TeV2017"))
                            errorReset          = 0.3+0.02*ptBins[k]+0.078*pow(ptBins[k],1.55);
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = errorReset*0.01;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = errorReset*0.01;
                    }
                } else if (meson.Contains("Eta")){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (!energy.CompareTo("8TeV"))
                            errorReset          = 2+pow(ptBins[k],2)*0.09;
                        else if (!energy.CompareTo("7TeV"))
                            errorReset          = 0.8+pow(ptBins[k],2)*0.035;
                        else if (!energy.CompareTo("2.76TeV"))
                            errorReset          = 2*(0.8+pow(ptBins[k],2)*0.03);
                        else if (!energy.CompareTo("5TeV2017"))
                            errorReset          = 0.8+pow(ptBins[k],2)*0.035;
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = errorReset*0.01;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = errorReset*0.01;
                    }
                }
            }

            // psipair and chi2 variation
            if ((nameCutVariationSC[i].CompareTo("PsiPair")==0 || nameCutVariationSC[i].CompareTo("Chi2")==0)){
                if (!meson.CompareTo("Pi0")){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (!energy.CompareTo("8TeV"))
                            errorReset          = 1.5+pow(ptBins[k]-3,2)*0.03;
                        else if (!energy.CompareTo("7TeV"))
                            errorReset          = 0.75+pow(ptBins[k],2)*0.015;
                        else if (!energy.CompareTo("2.76TeV"))
                            errorReset          = 1.25+pow(ptBins[k],2)*0.075;
                        else if (!energy.CompareTo("900GeV"))
                            errorReset          = 6.-5.4*ptBins[k]+1.85*pow(ptBins[k],2);
                        else if (!energy.CompareTo("5TeV2017"))
                            errorReset          = 1.8+8.52726/pow(14.8707,ptBins[k]);
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = errorReset*0.01;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = errorReset*0.01;
                    }
                } else if (meson.Contains("Eta")){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (!energy.CompareTo("8TeV")){
                            if (ptBins[k] > 2.5)
                                errorReset      = 2.5+pow(ptBins[k]-3,2)*0.275;
                            else
                                errorReset      = 2.5+pow(ptBins[k]-2.5,2)*1.1;
                        } else if (!energy.CompareTo("7TeV")){
                            errorReset          = 1.2+pow(ptBins[k]-3,2)*0.15;
                        } else if (!energy.CompareTo("2.76TeV")){
                            errorReset          = 3.2+pow(ptBins[k]-3,2)*0.09;
                        } else if (!energy.CompareTo("5TeV2017")){
                            errorReset          = 3.2+pow(ptBins[k]-3,2)*0.05;
                        }
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = errorReset*0.01;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = errorReset*0.01;
                    }
                }
            }

            // single pt variation
            if (nameCutVariationSC[i].CompareTo("SinglePt")==0 ){
                if (!meson.CompareTo("Pi0")){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (!energy.CompareTo("8TeV")){
                            errorReset          = 1.0;
                        } else if (!energy.CompareTo("7TeV")){
                          if (ptBins[k] > 1)
                              errorReset      = 0.5;
                          else
                              errorReset      = 0.5+pow(ptBins[k]-1,2)*2.35;
                        } else if (!energy.CompareTo("2.76TeV")){
                            if (ptBins[k] > 3)
                                errorReset      = 0.75;
                            else
                                errorReset      = 0.75+pow(ptBins[k]-3,2)*0.35;
                        } else if (!energy.CompareTo("900GeV")){
                            errorReset          = 1.8-1.4*ptBins[k]+0.46*pow(ptBins[k],2);
                        } else if (!energy.CompareTo("5TeV2017")){
                            errorReset          = 0.35+-0.065*ptBins[k]+0.0141391*pow(ptBins[k],2);
                        }
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = errorReset*0.01;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = errorReset*0.01;
                    }
                } else if (meson.Contains("Eta")){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (!energy.CompareTo("8TeV")){
                            if (ptBins[k] > 1.6)
                                errorReset      = 1.5+pow(ptBins[k]-1.5,2)*0.03;
                            else
                                errorReset      = 1.5+pow(ptBins[k]-1.6,2)*3;
                        } else if (!energy.CompareTo("7TeV")){
                            errorReset          = 1.76187+-0.69099*ptBins[k]+0.071933*pow(ptBins[k],2);
                        } else if (!energy.CompareTo("2.76TeV")){
                            if (ptBins[k] > 3)
                                errorReset      = 2*0.75;
                            else
                                errorReset      = 2*(0.75+pow(ptBins[k]-3,2)*0.35);
                        } else if (!energy.CompareTo("900GeV")){
                            errorReset          = 3.;
                        } else if (!energy.CompareTo("5TeV2017")){
                            if (ptBins[k] > 2)
                                errorReset      = 0.36;
                            else
                                errorReset      = 1.7*(0.085+1.6/pow(3.9,ptBins[k]));
                        }
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = errorReset*0.01;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = errorReset*0.01;
                    }
                }
            }

            // alpha variation
            if (nameCutVariationSC[i].CompareTo("Alpha")==0 ){
                if (!meson.CompareTo("Pi0")){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (!energy.CompareTo("8TeV"))
                            if (ptBins[k] > 5.0)
                                errorReset      = 1.2;
                            else
                                errorReset      = 0.1+pow(ptBins[k],2)*0.035;
                        else if (!energy.CompareTo("7TeV"))
                            errorReset          = 0.1+pow(ptBins[k],2)*0.035;
                        else if (!energy.CompareTo("2.76TeV"))
                            errorReset          = 0.4+pow(ptBins[k],2)*0.015;
                        else if (!energy.CompareTo("900GeV")){
                            if (ptBins[k] > 1.5)
                                errorReset      = 0.1+pow(ptBins[k]-1.5,2)*0.7;
                            else
                                errorReset      = 0.1;
                        } else if (!energy.CompareTo("5TeV2017")){
                            errorReset          = 0.25+-0.0292745*ptBins[k]+0.018*pow(ptBins[k],2);
                        }
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = errorReset*0.01;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = errorReset*0.01;
                    }
                } else if (meson.Contains("Eta")){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (!energy.CompareTo("8TeV"))
                            errorReset          = 0.4+pow(ptBins[k],2)*0.045;
                        else if (!energy.CompareTo("7TeV"))
                            errorReset          = 0.1+pow(ptBins[k],2)*0.045;
                        else if (!energy.CompareTo("2.76TeV"))
                            errorReset          = 0.4+pow(ptBins[k],2)*0.015;
                        else if (!energy.CompareTo("5TeV2017"))
                            errorReset          = 1.8*(0.926229+-0.241193*ptBins[k]+0.0278201*pow(ptBins[k],2));
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = errorReset*0.01;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = errorReset*0.01;
                    }
                }
            }

            // Double counting
            if (nameCutVariationSC[i].CompareTo("DoubleCount")==0 ){
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("5TeV2017")){
                        if (!meson.CompareTo("Pi0"))
                            errorReset  = 0.25;
                        else if (meson.Contains("Eta")){
                            errorReset  = 4*0.25;
                        }
                    }
                    errorsMean[i][k]        = errorReset;
                    errorsMeanErr[i][k]     = errorReset*0.01;
                    errorsMeanCorr[i][k]    = errorReset;
                    errorsMeanErrCorr[i][k] = errorReset*0.01;
                }
            }

            // yield variation
            if (nameCutVariationSC[i].CompareTo("YieldExtraction_pp")==0 ) {
                // smoothing for eta normal and etat for eta/pi0 ratio
                if( meson.CompareTo("Eta")==0 || (!meson.CompareTo("EtaToPi0") && i==0) ){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (!energy.CompareTo("8TeV")){
                            if (ptBins[k] > 2.0)
                                errorReset      = 5.+pow(ptBins[k]-2,2)*0.17;
                            else
                                errorReset      = 5.+pow(ptBins[k]-2,2)*10;
                        } else if (!energy.CompareTo("7TeV")){
                            if (ptBins[k] > 3.0)
                                errorReset      = 2.5+pow(ptBins[k]-3,2)*0.15;
                            else
                                errorReset      = 2.5+pow(ptBins[k]-3,2)*1.7;
                        } else if (!energy.CompareTo("2.76TeV")){
                            errorReset          = 6.0;
                        } else if (!energy.CompareTo("900GeV")){
                            errorReset          = 14.;
                        } else if (!energy.CompareTo("5TeV2017")){
                            errorReset      = 6.88139+-2.50537*ptBins[k]+0.275235*pow(ptBins[k],1.9);
                        }
                        errorsMean[i][k]        = errorReset;
                        errorsMeanErr[i][k]     = errorReset*0.01;
                        errorsMeanCorr[i][k]    = errorReset;
                        errorsMeanErrCorr[i][k] = errorReset*0.01;
                    }
                // smoothing for pi0 and pi0 for eta/pi0 ratio
                } else if (meson.CompareTo("Pi0")==0 || (!meson.CompareTo("EtaToPi0") && i==1) ) {
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (!energy.CompareTo("8TeV")){
                            // original Pi0 binning
                            if (ptBins[k] > 1.0)
                                errorReset          = 3.+pow(ptBins[k]-1,2)*0.04; //
                            else
                                errorReset          = 3.+pow(ptBins[k]-1,2)*18;
                            // or if Pi0 in eta binning prefered (no pt distinction):
                            //  errorReset          = 2.+pow(ptBins[k]-2.5,2)*0.18;
                            errorsMean[i][k]        = errorReset;
                            errorsMeanErr[i][k]     = errorReset*0.01;
                            errorsMeanCorr[i][k]    = errorReset;
                            errorsMeanErrCorr[i][k] = errorReset*0.01;
                            // by adding the pi0etabinning error, set the eta yield extraction error to the total eta/pi0 yield error
                            errorsMean[0][k]            = pow(pow(errorsMean[0][k],2)+pow(errorsMean[i][k],2),0.5);
                            errorsMeanErr[0][k]         = pow(pow(errorsMeanErr[0][k]*0.01,2)+pow(errorsMeanErr[i][k]*0.01,2),0.5);
                            errorsMeanCorr[0][k]        = pow(pow(errorsMeanCorr[0][k],2)+pow(errorsMeanCorr[i][k],2),0.5);
                            errorsMeanErrCorr[0][k]     = pow(pow(errorsMeanErrCorr[0][k]*0.01,2)+pow(errorsMeanErrCorr[i][k]*0.01,2),0.5);
                        } else if (!energy.CompareTo("7TeV")){
                            if (ptBins[k] > 2.0)
                                errorReset      = 1.5+pow(ptBins[k]-2,2)*0.03;
                            else
                                errorReset      = 1.5+pow(ptBins[k]-2,2)*1.5;
                            errorsMean[i][k]        = errorReset;
                            errorsMeanErr[i][k]     = errorReset*0.01;
                            errorsMeanCorr[i][k]    = errorReset;
                            errorsMeanErrCorr[i][k] = errorReset*0.01;
                        } else if (!energy.CompareTo("2.76TeV")){
                            if (errorsMean[i][k] > 7.0){
                                errorReset              = 5;
                                errorsMean[i][k]        = errorReset;
                                errorsMeanErr[i][k]     = errorReset*0.01;
                                errorsMeanCorr[i][k]    = errorReset;
                                errorsMeanErrCorr[i][k] = errorReset*0.01;
                            }
                        } else if (!energy.CompareTo("900GeV")){
                            errorReset              = 18.-18.2*ptBins[k]+4.74*pow(ptBins[k],2);
                            errorsMean[i][k]        = errorReset;
                            errorsMeanErr[i][k]     = errorReset*0.01;
                            errorsMeanCorr[i][k]    = errorReset;
                            errorsMeanErrCorr[i][k] = errorReset*0.01;
                        }
                    }
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
        } // end smoothing

        for (Int_t l = 0; l < nPtBins; l++){
            errorsPosSummed[l]      = errorsPosSummed[l]+pow(errorsPos[i][l],2);
            errorsNegSummed[l]      = errorsNegSummed[l]+ pow(errorsNeg[i][l],2);
            errorsMeanSummed[l]     = errorsMeanSummed[l]+ pow(errorsMean[i][l],2);
            errorsPosCorrSummed[l]  = errorsPosCorrSummed[l]+pow(errorsPosCorr[i][l],2);
            errorsNegCorrSummed[l]  = errorsNegCorrSummed[l] +pow(errorsNegCorr[i][l],2);
            errorsMeanCorrSummed[l] = errorsMeanCorrSummed[l]+ pow(errorsMeanCorr[i][l],2);
        }
        negativeErrors[i]           = new TGraphErrors(nPtBins,ptBins ,errorsNeg[i] ,ptBinsErr ,errorsNegErr[i] );
        meanErrors[i]               = new TGraphErrors(nPtBins,ptBins ,errorsMean[i] ,ptBinsErr ,errorsMeanErr[i] );
        positiveErrors[i]           = new TGraphErrors(nPtBins,ptBins ,errorsPos[i] ,ptBinsErr ,errorsPosErr[i] );
        negativeErrorsCorr[i]       = new TGraphErrors(nPtBins,ptBins ,errorsNegCorr[i] ,ptBinsErr ,errorsNegErrCorr[i] );
        meanErrorsCorr[i]           = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorr[i] ,ptBinsErr ,errorsMeanErrCorr[i] );
        meanErrorsCorrBefore[i]     = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrBefore[i] ,ptBinsErr ,errorsMeanErrCorrBefore[i] );
        positiveErrorsCorr[i]       = new TGraphErrors(nPtBins,ptBins ,errorsPosCorr[i] ,ptBinsErr ,errorsPosErrCorr[i] );

        if ( !meson.CompareTo("EtaToPi0") && !nameCutVariationSC[i].CompareTo("YieldExtraction_pp") && i==1){
            // set the pi0etabinning yield error to zero as both yield errors are already combined in the eta yield error
            for (Int_t k = 0; k < nPtBins; k++){
                if (ptBins[k] > -10){
                    errorsMean[i][k]            = 0;
                    errorsMeanErr[i][k]         = 0;
                    errorsMeanCorr[i][k]        = 0;
                    errorsMeanErrCorr[i][k]     = 0;
                }
            }
        }
    }

    Double_t errorMaterial = 4.50;
    if (meson.CompareTo("EtaToPi0") == 0)
        errorMaterial       = 0.;


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
    }

    Double_t errorsMat[nPtBins];
    Double_t errorsErrMat[nPtBins];
    for (Int_t l = 0; l < nPtBins; l++){
        errorsMat[l]    = 2*errorMaterial;
        errorsErrMat[l] = 2*errorMaterial/100;
    }

    TGraphErrors* graphMaterialError    = new TGraphErrors(nPtBins,ptBins ,errorsMat ,ptBinsErr ,errorsErrMat );
    negativeErrorsSummed                = new TGraphErrors(nPtBins,ptBins ,errorsNegSummed ,ptBinsErr ,errorsNegErrSummed );
    negativeErrorsCorrSummed            = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrSummed ,ptBinsErr ,errorsNegErrCorrSummed );
    positiveErrorsSummed                = new TGraphErrors(nPtBins,ptBins ,errorsPosSummed ,ptBinsErr ,errorsPosErrSummed );
    positiveErrorsCorrSummed            = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrSummed ,ptBinsErr ,errorsPosErrCorrSummed );
    meanErrorsSummed                    = new TGraphErrors(nPtBins,ptBins ,errorsMeanSummed ,ptBinsErr ,errorsMeanErrSummed );
    meanErrorsCorrSummed                = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSummed ,ptBinsErr ,errorsMeanErrCorrSummed );
    meanErrorsCorrSummedIncMat          = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrMatSummed ,ptBinsErr ,errorsMeanErrCorrMatSummed );


    // Give legend position for plotting
    Double_t minXLegend     = 0.12;
    Double_t maxYLegend     = 0.95;
    if (meson.CompareTo("Eta") == 0){
        minXLegend          = 0.23;
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
            labelMeson= new TLatex(0.75,0.89,Form("#eta/#pi^{0} rec. #gamma_{conv}"));
        } else if (meson.Contains("Pi0")){
            labelMeson= new TLatex(0.75,0.89,Form("#pi^{0} #rightarrow #gamma_{conv}#gamma_{conv}"));
        } else {
            labelMeson= new TLatex(0.75,0.89,Form("#eta #rightarrow #gamma_{conv}#gamma_{conv}"));
        }
        SetStyleTLatex( labelMeson, 0.038,4,1,42,kTRUE, 311);
        labelMeson->Draw();

        TLatex *labelCentrality = new TLatex(0.75,0.93,Form("%s",collisionSystem.Data() ));
        SetStyleTLatex( labelCentrality, 0.038,4,1,42,kTRUE, 11);
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
        if (meson.CompareTo("EtaToPi0")){
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
        histo2DNewSysErrMean->GetYaxis()->SetRangeUser(0.,10.);
        histo2DNewSysErrMean->Draw();

            if (!showBeforeSmoothing) continue;
//             if (bsmooth[cut]) continue;
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

            meanErrorsCorrBefore[cut]->Fit(pol4,"NRMEX0+","",minPt,maxPt);
            meanErrorsCorrBefore[cut]->Fit(pol2,"NRMEX0+","",minPt,maxPt);
            meanErrorsCorrBefore[cut]->Fit(pol1,"NRMEX0+","",minPt,maxPt);
            meanErrorsCorrBefore[cut]->Fit(pol0,"NRMEX0+","",minPt,maxPt);
            meanErrorsCorrBefore[cut]->Fit(bla,"NRMEX0+","",minPt,maxPt);

            cout << "The following functional forms can be used for smoothing of sys variation " << cut << endl;
            cout << "blk (pol0):\t " << pol0->GetParameter(0) << endl;
            cout << "cyan(pol1):\t " << pol1->GetParameter(0) << "+" << pol1->GetParameter(1) << "*ptBins[k]" << endl;
            cout << "blue(pol2):\t " << pol2->GetParameter(0) << "+" << pol2->GetParameter(1) << "*ptBins[k]+" << pol2->GetParameter(2) << "*pow(ptBins[k],2)" << endl;
            cout << "red (pol4):\t " << pol4->GetParameter(0) <<"+"<< pol4->GetParameter(1)<< "*ptBins[k]+"<< pol4->GetParameter(2)<< "*pow(ptBins[k],2)+"<< pol4->GetParameter(3)<< "*pow(ptBins[k],4)" << endl;
            cout << "pink (pow):\t " << bla->GetParameter(0) << "+" << bla->GetParameter(1) << "/pow(" << bla->GetParameter(2) << ",ptBins[k])" << endl;

            pol4->SetLineColor(kRed+2);
            pol2->SetLineColor(kBlue+2);
            pol1->SetLineColor(kCyan+2);
            pol0->SetLineColor(kBlack);
            bla->SetLineColor(kMagenta+2);

            DrawGammaSetMarkerTGraphErr(meanErrorsCorrBefore[cut], 20+cut, 1.,kGray+1,kGray+1);
            meanErrorsCorrBefore[cut]->Draw("p,csame");
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[cut], 20+cut, 1.,color[cut],color[cut]);
            meanErrorsCorr[cut]->Draw("p,csame");
            pol4->Draw("same");
            pol2->Draw("same");
            pol1->Draw("same");
            pol0->Draw("same");
            bla->Draw("same");


            TLegend* legendMeanNewSingle = GetAndSetLegend2(minXLegend,maxYLegend-heightLegend,minXLegend+widthLegend,maxYLegend, 30);
            legendMeanNewSingle->SetMargin(0.1);
            legendMeanNewSingle->AddEntry(meanErrorsCorr[cut],nameCutVariation[cut].Data(),"p");
            legendMeanNewSingle->Draw();

        canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedConv/SysMeanNewWithMeanSingle_%s_%s%s_%s_Variation%d.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),cut,suffix.Data()));
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
    for (Int_t l=0; l< nPtBins-1; l++){
//         SysErrDatAver  << ptBins[l] <<"\t"<< "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
        SysErrDatAver  << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
        // 			SysErrDatAver << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
    }
    SysErrDatAver.close();

    const char *SysErrDatnameMeanSingleErr = Form("SystematicErrorsCalculatedConv/SystematicErrorAveragedSinglePCM_%s_%s%s_%s.dat", meson.Data(), energyForOutput.Data(), additionalNameOutput.Data(), dateForOutput.Data());
    fstream SysErrDatAverSingle;
    cout << SysErrDatnameMeanSingleErr << endl;
    SysErrDatAverSingle.open(SysErrDatnameMeanSingleErr, ios::out);
    SysErrDatAverSingle << "Pt bin\t" ;
    for (Int_t i= 0; i< numberCutStudies; i++){
        if(!meson.CompareTo("EtaToPi0")&&i==1)
            continue;
        SysErrDatAverSingle << nameCutVariationSC[i] << "\t";
    }
    if(meson.CompareTo("EtaToPi0"))
        SysErrDatAverSingle << "Material";
    SysErrDatAverSingle << endl;
    for (Int_t l=0;l< nPtBins-1;l++){
        SysErrDatAverSingle << ptBins[l] << "\t";
        for (Int_t i= 0; i< numberCutStudies; i++){
            if(!meson.CompareTo("EtaToPi0")&&i==1)
                continue;
            SysErrDatAverSingle << errorsMeanCorr[i][l] << "\t";
        }
        if(meson.CompareTo("EtaToPi0"))
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
        //"0 YieldExtraction_pp"
        //"1 BGEstimate_pp"
        //"2 dEdxE"
        //"3 dEdxPi"
        //"4 TPCCluster"
        //"5 SinglePt"
        //"6 Chi2"
        //"7 Qt"
        //"8 Alpha"
        //"9 BG"
        //"10 Periods"
        //"11 SPD"

        // Signal extraction: Yield extraction, Alpha ,BG, MC Smearing
        errorsMeanCorrSignalExtraction[l]   =   TMath::Sqrt(pow(errorsMeanCorr[0][l],2)+    // Yield extraction
                                                            pow(errorsMeanCorr[8][l],2)+    // Alpha
                                                            pow(errorsMeanCorr[9][l],2));   // BG

        if (meson.CompareTo("EtaToPi0") == 0){
            errorsMeanCorrSignalExtraction[l]   =   TMath::Sqrt(pow(errorsMeanCorr[0][l],2)+    // Yield extraction eta
                                                                pow(errorsMeanCorr[8][l],2)+    // Alpha
                                                                pow(errorsMeanCorr[9][l],2)+    // BG
                                                                pow(errorsMeanCorr[1][l],2));   // Yield extraction pi0

        }
        // PID: dEdxE, dEdxPi
        errorsMeanCorrPID[l]                =   TMath::Sqrt(pow(errorsMeanCorr[2][l],2)+    // dEdxE
                                                            pow(errorsMeanCorr[3][l],2));   // dEdxPi

        // photon reco: Chi2, Qt, PsiPair, CosPoint, MinR
        errorsMeanCorrPhotonReco[l]         =   TMath::Sqrt(pow(errorsMeanCorr[6][l],2)+    // Chi2
                                                            pow(errorsMeanCorr[7][l],2));   // Qt

        // track reconstruction: TPCCluster, Single pT, Eta
        errorsMeanCorrTrackReco[l]          =   TMath::Sqrt(pow(errorsMeanCorr[4][l],2)+    // TPCCluster
                                                            pow(errorsMeanCorr[5][l],2));   // Single pT
        // pileup
        if(!energy.CompareTo("8TeV")||!energy.CompareTo("7TeV"))
            errorsMeanCorrPileup[l]         =   TMath::Sqrt(pow(errorsMeanCorr[1][l],2)+ // out of bunch methods
                                                            pow(errorsMeanCorr[11][l],2)+// SPD
                                                            pow(errorsMeanCorr[12][l],2));   // DCA out of bunch

    }
    TGraphErrors* meanErrorsPID                 = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPID ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsPhotonReco          = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPhotonReco ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsSignalExtraction    = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSignalExtraction ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsTrackReco           = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrTrackReco ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsPileup              = NULL;
    if(!energy.CompareTo("8TeV")||!energy.CompareTo("7TeV"))
        meanErrorsPileup              = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPileup ,ptBinsErr ,errorsMeanErrCorrSummed );


    // ***************************************************************************************************
    // ********************* Plot grouped errors for better understanding ********************************
    // ***************************************************************************************************
    Double_t minXLegend2 = 0.13;
    Double_t maxYLegend2 = 0.95;
    if (meson.CompareTo("Eta") == 0){
        minXLegend2 = 0.20;
    }
    Double_t widthLegend2 = 0.52;
    Double_t heightLegend2 = 0.15;

    TCanvas* canvasSummedErrMean = new TCanvas("canvasSummedErrMean","",200,10,1350,900);// gives the page size
    DrawGammaCanvasSettings( canvasSummedErrMean, 0.08, 0.01, 0.015, 0.09);

        // create dummy histo
        TH2D *histo2DSummedErrMean ;
        if (meson.Contains("Pi0") ){
            histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,25.);
        } else {
            histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,35.);
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
        DrawGammaSetMarkerTGraphErr(meanErrorsPID, 21, 1.,color[1],color[1]);
        meanErrorsPID->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsPID,"Electron PID","p");
        cout << "here" << endl;
        DrawGammaSetMarkerTGraphErr(meanErrorsTrackReco, 22, 1.,color[2],color[2]);
        meanErrorsTrackReco->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsTrackReco,"Track Reconstruction","p");
        cout << "here" << endl;
        DrawGammaSetMarkerTGraphErr(meanErrorsPhotonReco, 23, 1.,color[3],color[3]);
        meanErrorsPhotonReco->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsPhotonReco,"Photon Reconstruction","p");
        cout << "here" << endl;
        if (meson.CompareTo("EtaToPi0")){
            if(!energy.CompareTo("8TeV")||!energy.CompareTo("7TeV")) {
                DrawGammaSetMarkerTGraphErr(meanErrorsPileup, 25, 1.,color[5],color[5]);
                meanErrorsPileup->Draw("p,csame");
                legendSummedMeanNew->AddEntry(meanErrorsPileup,"Pileup Estimate","p");
            }else{
                DrawGammaSetMarkerTGraphErr(meanErrorsCorr[1], 25, 1.,color[5],color[5]);
                meanErrorsCorr[1]->Draw("p,csame");
                legendSummedMeanNew->AddEntry(meanErrorsCorr[1],"Pileup Estimate","p");
            }
            cout << "here" << endl;
            DrawGammaSetMarkerTGraphErr(graphMaterialError, 24, 1.,color[4],color[4]);
            graphMaterialError->Draw("p,csame");
            legendSummedMeanNew->AddEntry(graphMaterialError,"Material","p");
        } else {
            if (meanErrorsCorr[10]&&!energy.CompareTo("900GeV")){
                DrawGammaSetMarkerTGraphErr(meanErrorsCorr[10], 25, 1.,color[5],color[5]);
                meanErrorsCorr[10]->Draw("p,csame");
                legendSummedMeanNew->AddEntry(meanErrorsCorr[10],"Pileup Estimate","p");
                cout << "here900" << endl;
            } else if (meanErrorsCorr[11]){
                DrawGammaSetMarkerTGraphErr(meanErrorsCorr[11], 25, 1.,color[5],color[5]);
                meanErrorsCorr[11]->Draw("p,csame");
                legendSummedMeanNew->AddEntry(meanErrorsCorr[11],"Pileup Estimate","p");
                cout << "here" << endl;
            }
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
        SysErrDatAverPaper << ptBins[l] <<"\t" << errorMaterial*2 << "\t" << errorsMeanCorrSignalExtraction[l] << "\t" << errorsMeanCorrPID[l]<< "\t" << errorsMeanCorrPhotonReco[l]<< "\t" <<errorsMeanCorrTrackReco[l] <<"\t" << errorsMeanCorrMatSummed[l]<< endl;
    }
    SysErrDatAverPaper.close();

}