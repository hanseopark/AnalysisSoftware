// ------------------------------------------------
// This code finalizes the systematic errors for --
// the direct photon analysis. Gamma spectrum as --
// well as double ratio systematics are produced --
// ------------------------------------------------
// Code is maintained by:                        --
//          - Friederike Bock                    --
// for the Photon Conversion Group               --
// ------------------------------------------------
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

//*************************************************************************************************************
//*********** Main function to calculate gamma errors for calo measurements in pp collisions ******************
//*************************************************************************************************************
void FinaliseSystematicErrorsCalo_Gammas_pp(    TString nameDataFileErrors      = "",
                                                TString energy                  = "",
                                                TString spectrumName            = "",
                                                Int_t numberCutStudies          = 1,
                                                Double_t minPt                  = 0,
                                                Double_t maxPt                  = 10,
                                                TString suffix                  = "pdf",
                                                Int_t mode                      = 0
                                            ){

    //**************************************************************
    // Set plotting style and labels
    StyleSettingsThesis();
    SetPlotStyle();

    Double_t textSizeSpectra                    = 0.04;
    TString date                                = ReturnDateString();
    TString collisionSystem                     = ReturnFullCollisionsSystem(energy);
    TString detectionProcess                    = ReturnFullTextReconstructionProcess(4);


    TString dateForOutput                       = ReturnDateStringForOutput();
    TString energyForOutput                     = ReturnCollisionEnergyOutputString(energy);

    TLatex *labelGamma                          = new TLatex(0.95,0.88,detectionProcess);
    SetStyleTLatex( labelGamma, 0.038,4);
    labelGamma->SetTextAlign(31);

    Double_t  posYLabel                         = 0.84;

    TLatex *labelEnergy                         = new TLatex(0.95,0.92,collisionSystem);
    SetStyleTLatex( labelEnergy, 0.038,4);
    labelEnergy->SetTextAlign(31);
    TLatex *labelSpectrum = 0x0;
    if(!spectrumName.CompareTo("IncRatio"))
        labelSpectrum                           = new TLatex(0.95,posYLabel,"#gamma_{inc}/#pi^{0}");
    if(!spectrumName.CompareTo("DoubleRatio"))
        labelSpectrum                           = new TLatex(0.95,posYLabel,"R_{#gamma}");
    if(!spectrumName.CompareTo("Gamma"))
        labelSpectrum                           = new TLatex(0.95,posYLabel,"#gamma_{inc}");
    SetStyleTLatex( labelSpectrum, 0.038,4);
    labelSpectrum->SetTextAlign(31);


    //**************************************************************
    // setup to read different error sources from root file
    Int_t numberOfEntriesPos                    = 0;
    Int_t numberOfEntriesNeg                    = 0;

    // read root file with different cutvariations
    TFile* fileErrorInput                       = new TFile(nameDataFileErrors.Data());

    const Int_t nPtBins                         = 60;
    Int_t nPtBinsActive                         = 60;
    const Int_t nCuts                           = numberCutStudies;
    Int_t nCutsActive                           = nCuts;
    Double_t* ptBins                            = 0;
    Double_t* ptBinsErr                         = 0;

    Double_t yRangesSysPlotting[2]              = {-0.5,24.9};
    if (energy.CompareTo("2.76TeV") == 0 ){
        if (spectrumName.Contains("Gamma"))
            yRangesSysPlotting[1]               = 20.5;
        else if (spectrumName.Contains("IncRatio"))
            yRangesSysPlotting[1]               = 20.5;
        else
            yRangesSysPlotting[1]               = 20.5;
    }else if (energy.CompareTo("8TeV") == 0 ){
      if (spectrumName.Contains("Gamma"))
          yRangesSysPlotting[1]               = 15.5;
      else if (spectrumName.Contains("IncRatio"))
          yRangesSysPlotting[1]               = 15.5;
      else
          yRangesSysPlotting[1]               = 15.5;
    }

    // Set names of cut variations for file input
    TString nameCutVariationSC[15]              = { "ClusterNonLinearity", "ClusterM02", "ClusterMinEnergy", "ClusterTrackMatchingCalo", "ClusterTiming",
                                                    "ClusterNCells", "ClusterMaterialTRD",  "SPD", "BG","Alpha" ,
                                                    "Cocktail", "IntRange", "Efficiency", "Rapidity", "OpeningAngle"
                                                  };

    // Set colors and markers
    Color_t color[15];
    Color_t markerStyle[15];
    // Set names of cut variations for legends
    TString nameCutVariation[15];

    for (Int_t k =0; k<nCuts; k++ ){
        cout << "variation: " << nameCutVariationSC[k].Data() << endl;
        color[k]                                = GetColorSystematics( nameCutVariationSC[k] );
        markerStyle[k]                          = GetMarkerStyleSystematics( nameCutVariationSC[k] );
        nameCutVariation[k]                     = GetSystematicsName(nameCutVariationSC[k]);
        cout << "name for writing: " << nameCutVariation[k].Data() << "\t"<< color[k]  << "\t" << markerStyle[k] << endl;
    }

    if(!energy.CompareTo("7TeV")){
        nameCutVariationSC[12]                   = "7TeVEfficiency";
    }else if(!energy.CompareTo("8TeV")){
        nameCutVariationSC[12]                   = "8TeVEfficiency";
    }

    // Create output folder
    gSystem->Exec("mkdir -p GammaSystematicErrorsCalculated");

    // ***************************************************************************************************
    // ******************************** Booleans for enabling systematic errors **************************
    // ***************************************************************************************************
    Bool_t benable[15]                          = { 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0 };
    Bool_t benableIncGammapp2760GeV[15]         = { 1, 1, 1, 1, 1,  1, 1, 1, 0, 0,
                                                    0, 0, 1, 0, 0 };
    Bool_t benableIncRatiopp2760GeV[15]         = { 1, 1, 1, 1, 1,  1, 1, 0, 0, 1,
                                                    0, 1, 1, 1, 1 };
    Bool_t benableDRpp2760GeV[15]               = { 1, 1, 1, 1, 1,  1, 1, 0, 0, 1,
                                                    1, 1, 1, 1, 1 };

    Bool_t benableIncGammapp8TeV[15]            = { 1, 1, 1, 1, 1,  1, 1, 1, 0, 0,
                                                    0, 0, 1, 0, 0 };
    Bool_t benableIncRatiopp8TeV[15]            = { 1, 1, 1, 1, 1,  1, 1, 0, 0, 1,
                                                    0, 1, 1, 0, 1 };
    Bool_t benableDRpp8TeV[15]                  = { 1, 1, 1, 1, 1,  1, 1, 0, 0, 1,
                                                    1, 1, 1, 0, 1 };
    // ***************************************************************************************************
    // ******************************** Booleans for smoothing *******************************************
    // ***************************************************************************************************
    Bool_t bsmooth[15]                          = { 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0 };
    Bool_t bsmoothIncGammapp2760GeV[15]         = { 1, 0, 1, 1, 1,  1, 1, 1, 0, 0,
                                                    0, 0, 1, 0, 0 };
    Bool_t bsmoothIncRatiopp2760GeV[15]         = { 1, 0, 1, 1, 1,  1, 1, 1, 0, 1,
                                                    0, 0, 1, 1, 0 };
    Bool_t bsmoothDRpp2760GeV[15]               = { 1, 0, 1, 1, 1,  1, 1, 1, 0, 1,
                                                    0, 0, 1, 1, 0 };

    Bool_t bsmoothIncGammapp8TeV[15]            = { 1, 1, 1, 1, 1,  1, 1, 1, 0, 1,
                                                    1, 0, 1, 0, 0 };
    Bool_t bsmoothIncRatiopp8TeV[15]            = { 1, 1, 1, 1, 1,  1, 1, 1, 0, 1,
                                                    1, 1, 1, 0, 1 };
    Bool_t bsmoothDRpp8TeV[15]                  = { 1, 1, 1, 1, 1,  1, 1, 1, 0, 1,
                                                    1, 1, 1, 0, 1 };

    for (Int_t i = 0; i < numberCutStudies; i++){
        if (energy.CompareTo("2.76TeV") == 0){
            if(!spectrumName.CompareTo("IncRatio")){
                bsmooth[i]                      = bsmoothIncRatiopp2760GeV[i];
                benable[i]                      = benableIncRatiopp2760GeV[i];
            } else if(!spectrumName.CompareTo("DoubleRatio")){
                bsmooth[i]                      = bsmoothDRpp2760GeV[i];
                benable[i]                      = benableDRpp2760GeV[i];
            } else if(!spectrumName.CompareTo("Gamma")){
                bsmooth[i]                      = bsmoothIncGammapp2760GeV[i];
                benable[i]                      = benableIncGammapp2760GeV[i];
            }
        }else if (energy.CompareTo("8TeV") == 0){
          if(!spectrumName.CompareTo("IncRatio")){
              bsmooth[i]                      = bsmoothIncRatiopp8TeV[i];
              benable[i]                      = benableIncRatiopp8TeV[i];
          } else if(!spectrumName.CompareTo("DoubleRatio")){
              bsmooth[i]                      = bsmoothDRpp8TeV[i];
              benable[i]                      = benableDRpp8TeV[i];
          } else if(!spectrumName.CompareTo("Gamma")){
              bsmooth[i]                      = bsmoothIncGammapp8TeV[i];
              benable[i]                      = benableIncGammapp8TeV[i];
          }
        }
        if (!benable[i]) nCutsActive--;
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
        errorsPosSummed[l]                      = 0.;
        errorsNegSummed[l]                      = 0.;
        errorsMeanSummed[l]                     = 0.;
        errorsPosCorrSummed[l]                  = 0.;
        errorsNegCorrSummed[l]                  = 0.;
        errorsMeanCorrSummed[l]                 = 0.;
    }

    for (Int_t i = 0; i < nCuts; i++){
        if (!benable[i]) {
            cout << "*****************************************************************" << endl;
            cout << "skipping: " << nameCutVariationSC[i].Data() << endl;
            cout << "*****************************************************************" << endl;
            for (Int_t l = 0; l < nPtBins; l++){
                errorsMean[i][l]                = 0;
                errorsMeanErr[i][l]             = 0.0;
                errorsMeanCorr[i][l]            = 0;
                errorsMeanErrCorr[i][l]         = 0.0;
            }
            continue;
        }

        TGraphAsymmErrors* graphPosErrors       = NULL;
        TGraphAsymmErrors* graphNegErrors       = NULL;

        // Set currently undetermined uncertainties
        if (    nameCutVariationSC[i].CompareTo("SPD")==0  || nameCutVariationSC[i].CompareTo("Efficiency")==0 || nameCutVariationSC[i].CompareTo("ClusterTiming")==0 ||
                nameCutVariationSC[i].CompareTo("Alpha")==0 || nameCutVariationSC[i].CompareTo("Rapidity")==0){
            TString nameGraphPos;
            TString nameGraphNeg;
            nameGraphPos                        = Form("%s_SystErrorRelPos_%s_pp",spectrumName.Data(),nameCutVariationSC[0].Data());
            nameGraphNeg                        = Form("%s_SystErrorRelNeg_%s_pp",spectrumName.Data(),nameCutVariationSC[0].Data());
            cout << "Cutstudies " << i << "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<< "\t fixing" <<  endl;
            graphPosErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else {
            // Load input graphs from systematics file
            TString nameGraphPos                = Form("%s_SystErrorRelPos_%s_pp",spectrumName.Data(),nameCutVariationSC[i].Data());
            TString nameGraphNeg                = Form("%s_SystErrorRelNeg_%s_pp",spectrumName.Data(),nameCutVariationSC[i].Data());
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        }

        // Remove first points depending on chosen offset
        while (graphPosErrors->GetX()[0] < minPt ){
            graphPosErrors->RemovePoint(0);
            graphNegErrors->RemovePoint(0);
        }
        while (graphPosErrors->GetX()[graphPosErrors->GetN()-1] > maxPt){
            graphPosErrors->RemovePoint(graphPosErrors->GetN()-1);
            graphNegErrors->RemovePoint(graphNegErrors->GetN()-1);
        }
        if (i == 0) {
            nPtBinsActive                       = graphNegErrors->GetN();
            ptBins                              = graphNegErrors->GetX();
            ptBinsErr                           = graphNegErrors->GetEXhigh();
        }
        errorsNeg[i]                            = graphNegErrors->GetY();
        errorsNegErr[i]                         = graphNegErrors->GetEYhigh();
        errorsPos[i]                            = graphPosErrors->GetY();
        errorsPosErr[i]                         = graphPosErrors->GetEYhigh();

        // Calculate systematic error from input spectrum
        cout << nameCutVariationSC[i].Data() << endl;
        CalculateMeanSysErr(            errorsMean[i],  errorsMeanErr[i],   errorsPos[i],       errorsNeg[i],           nPtBinsActive);
        CorrectSystematicErrorsWithMean(errorsPos[i],   errorsPosErr[i],    errorsPosCorr[i],   errorsPosErrCorr[i],    nPtBinsActive);
        CorrectSystematicErrorsWithMean(errorsNeg[i],   errorsNegErr[i],    errorsNegCorr[i],   errorsNegErrCorr[i],    nPtBinsActive);
        CorrectSystematicErrorsWithMean(errorsMean[i],  errorsMeanErr[i],   errorsMeanCorr[i],  errorsMeanErrCorr[i],   nPtBinsActive);

        // ***************************************************************************************************
        // ************************ Adjust errors if requested to fixed values *******************************
        // ***************************************************************************************************
        if (bsmooth[i]){
            cout << "smoothing on" << endl;
            Double_t errorFixed                 = -1;
            Bool_t adjustPtDependent            = kFALSE;

            // 0 "ClusterNonLinearity"
            if (!nameCutVariationSC[i].CompareTo("ClusterNonLinearity")){
                if (spectrumName.Contains("Ratio")){
                    adjustPtDependent           = kTRUE;
                    for (Int_t k = 0; k < nPtBins; k++){
                        if(!energy.CompareTo("2.76TeV")){
                          errorFixed              = 1.2+(0.01)*ptBins[k]+(0.01)*ptBins[k]*ptBins[k]; // non lin only
                          errorFixed              = errorFixed+0.1*7;
                        }else if(!energy.CompareTo("8TeV")){
                          errorFixed   = 1.8+25./pow(10,ptBins[k]);
                        }
                        if (errorFixed != -1){
                            errorsMean[i][k]            = errorFixed;
                            errorsMeanErr[i][k]         = errorFixed*0.01;
                            errorsMeanCorr[i][k]        = errorFixed;
                            errorsMeanErrCorr[i][k]     = errorFixed*0.01;
                        }
                    }
                } else {
                    adjustPtDependent           = kTRUE;
                    for (Int_t k = 0; k < nPtBins; k++){
                        if(!energy.CompareTo("2.76TeV")){
                          errorFixed              = 0.65+(0.02)*ptBins[k]+(0.01)*ptBins[k]*ptBins[k];
                          errorFixed              = errorFixed+0.1*7;
                        }else if(!energy.CompareTo("8TeV")){
                          errorFixed   = 0.9+25./pow(10,ptBins[k]);
                        }
                        if (errorFixed != -1){
                            errorsMean[i][k]            = errorFixed;
                            errorsMeanErr[i][k]         = errorFixed*0.01;
                            errorsMeanCorr[i][k]        = errorFixed;
                            errorsMeanErrCorr[i][k]     = errorFixed*0.01;
                        }
                    }
                }
            }

            // 1 "ClusterM02"
            if (!nameCutVariationSC[i].CompareTo("ClusterM02")){
                if (spectrumName.Contains("Ratio")){
                    adjustPtDependent           = kTRUE;
                    for (Int_t k = 0; k < nPtBins; k++){
                        if(!energy.CompareTo("2.76TeV")){
                          if (ptBins[k] > 2.8)
                          errorFixed              =  1.45+(0.06)*ptBins[k]+(0.04)*ptBins[k]*ptBins[k];
                        }else if(!energy.CompareTo("8TeV")){
                          errorFixed              =  0.9;
                          if(ptBins[k] > 5.)
                            errorFixed           += 0.1*(ptBins[k]-5.)+0.01*(ptBins[k]-5.)*(ptBins[k]-5.);
                        }
                        if (errorFixed != -1){
                            errorsMean[i][k]            = errorFixed;
                            errorsMeanErr[i][k]         = errorFixed*0.01;
                            errorsMeanCorr[i][k]        = errorFixed;
                            errorsMeanErrCorr[i][k]     = errorFixed*0.01;
                        }
                    }
                } else {
                    adjustPtDependent           = kTRUE;
                    for (Int_t k = 0; k < nPtBins; k++){
                        if(!energy.CompareTo("2.76TeV")){
                          if (ptBins[k] > 5.5)
                              errorFixed              =  10.5;
                        }else if(!energy.CompareTo("8TeV")){
                          errorFixed              =  1.5;
                          if(ptBins[k] > 6.)
                            errorFixed           += 0.1*(ptBins[k]-6.);
                        }
                        if (errorFixed != -1){
                            errorsMean[i][k]            = errorFixed;
                            errorsMeanErr[i][k]         = errorFixed*0.01;
                            errorsMeanCorr[i][k]        = errorFixed;
                            errorsMeanErrCorr[i][k]     = errorFixed*0.01;
                        }
                    }

                }
            }

            // 2 "ClusterMinEnergy"
            if (!nameCutVariationSC[i].CompareTo("ClusterMinEnergy")){
                if (spectrumName.Contains("Ratio")){
                    adjustPtDependent           = kTRUE;
                    for (Int_t k = 0; k < nPtBins; k++){
                        if(!energy.CompareTo("2.76TeV")){
                          errorFixed                  = 0.8+55/pow(6.4,ptBins[k]);   ;
                        }else if(!energy.CompareTo("8TeV")){
                          errorFixed                  = 0.5+30/pow(7.4,ptBins[k]);   ;
                        }
                        if (errorFixed != -1){
                            errorsMean[i][k]            = errorFixed;
                            errorsMeanErr[i][k]         = errorFixed*0.01;
                            errorsMeanCorr[i][k]        = errorFixed;
                            errorsMeanErrCorr[i][k]     = errorFixed*0.01;
                        }
                    }
                } else {
                    if(!energy.CompareTo("2.76TeV")){
                      errorFixed      = 0.5;
                    }else if(!energy.CompareTo("8TeV")){
                      errorFixed      = 0.1;
                    }
                }
            }

            // 3 "ClusterTrackMatching"
            if (!nameCutVariationSC[i].CompareTo("ClusterTrackMatchingCalo")){
                if (spectrumName.Contains("Ratio")){
                    adjustPtDependent           = kTRUE;
                    for (Int_t k = 0; k < nPtBins; k++){
                        if(!energy.CompareTo("2.76TeV")){
                          if (ptBins[k] > 1.6)
                            errorFixed                  = 0.65+(0.003)*ptBins[k]+(0.014)*ptBins[k]*ptBins[k]; // parametrisation
                        }else if(!energy.CompareTo("8TeV")){
                            errorFixed              =  1.0;
                            if(ptBins[k] < 3.)
                              errorFixed           += 0.3*(ptBins[k]-3.)*(ptBins[k]-3.);
                            if(ptBins[k] > 7.)
                              errorFixed           += 0.015*(ptBins[k]-7.)*(ptBins[k]-7.)+0.1*(ptBins[k]-7.);
                        }
                        if (errorFixed != -1){
                            errorsMean[i][k]            = errorFixed;
                            errorsMeanErr[i][k]         = errorFixed*0.01;
                            errorsMeanCorr[i][k]        = errorFixed;
                            errorsMeanErrCorr[i][k]     = errorFixed*0.01;
                        }
                    }
                } else {
                    adjustPtDependent           = kTRUE;
                    for (Int_t k = 0; k < nPtBins; k++){
                        if(!energy.CompareTo("2.76TeV")){
                          if (ptBins[k] > 1.6)
                            errorFixed              =  0.65;
                        }else if(!energy.CompareTo("8TeV")){
                            errorFixed              =  1.5;
                            if(ptBins[k] < 3.)
                              errorFixed           += 0.2*(ptBins[k]-3.)*(ptBins[k]-3.);
                            if(ptBins[k] > 7.)
                              errorFixed           += 0.015*(ptBins[k]-7.)*(ptBins[k]-7.)+0.05*(ptBins[k]-7.);
                        }
                        if (errorFixed != -1){
                            errorsMean[i][k]            = errorFixed;
                            errorsMeanErr[i][k]         = errorFixed*0.01;
                            errorsMeanCorr[i][k]        = errorFixed;
                            errorsMeanErrCorr[i][k]     = errorFixed*0.01;
                        }
                    }

                }
            }

            // fix ClusterTiming pileup sys #4
            if (!nameCutVariationSC[i].CompareTo("ClusterTiming")){
                if(!energy.CompareTo("2.76TeV")){
                  errorFixed                  = 0.3;
                }else if(!energy.CompareTo("8TeV")){
                  errorFixed                  = 0.5;
                }
            }


            // 5 "ClusterNCells"
            if (!nameCutVariationSC[i].CompareTo("ClusterNCells")){
                if(!energy.CompareTo("2.76TeV")){
                  errorFixed                  = 0.23;
                }else if(!energy.CompareTo("8TeV")){
                  errorFixed                  = 0.23;
                }
            }

            // fix ClusterMaterialTRD sys #6
            if (!nameCutVariationSC[i].CompareTo("ClusterMaterialTRD")){
                if (spectrumName.Contains("Ratio")){
                    if(!energy.CompareTo("2.76TeV")){
                      errorFixed                      = 4.2;
                    }else if(!energy.CompareTo("8TeV")){
                      errorFixed                      = 3.0;
                    }
                }else{
                    if(!energy.CompareTo("2.76TeV")){
                        errorFixed                   = 2.8;
                    }else if(!energy.CompareTo("8TeV")){
                        errorFixed                   = 2.1;
                    }
                }
            }

            // fix SPD pileup sys #7
            if (!nameCutVariationSC[i].CompareTo("SPD")){
                if(!energy.CompareTo("2.76TeV")){
                    errorFixed                  = 0.25;
                }else if(!energy.CompareTo("8TeV")){
                    errorFixed                  = 0.1;
                }
            }

            // fix BG sys #8
            if (!nameCutVariationSC[i].CompareTo("BG")){
                adjustPtDependent               = kTRUE;
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if(!energy.CompareTo("2.76TeV")){
                      errorFixed                  = 0.3+pow(ptBins[k],2)*0.0035;
                    }
                    if (errorFixed != -1){
                        errorsMean[i][k]        = errorFixed;
                        errorsMeanErr[i][k]     = errorFixed*0.01;
                        errorsMeanCorr[i][k]    = errorFixed;
                        errorsMeanErrCorr[i][k] = errorFixed*0.01;
                    }
                }
            }

            // fix alpha uncertainties #9
            if (!nameCutVariationSC[i].CompareTo("Alpha")){
                if(!energy.CompareTo("2.76TeV")){
                    errorFixed                  = 0.25;
                }else if(!energy.CompareTo("8TeV")){
                    errorFixed                  = 0.25;
                }
            }

            // fix "Cocktail" #10
            if (!nameCutVariationSC[i].CompareTo("Cocktail")){
                if(!energy.CompareTo("8TeV")){
                    errorFixed                  = 3.5;
                }
            }

            // fix IntRange sys #11
            if (!nameCutVariationSC[i].CompareTo("IntRange")){
                if (spectrumName.Contains("Ratio")){
                    adjustPtDependent           = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                        if(!energy.CompareTo("2.76TeV")){
                          if (ptBins[k] > 1.6)
                              errorFixed          = 1.6+pow(ptBins[k],2)*0.018;
                        }else if(!energy.CompareTo("8TeV")){
                          errorFixed              = 1.25;
                          if(ptBins[k]>=6.) errorFixed += 0.02*(ptBins[k]-6.)*(ptBins[k]-6.);
                          if(ptBins[k]<=2.) errorFixed += -0.45+45/pow(10.,ptBins[k]);
                        }
                        if (errorFixed != -1){
                            errorsMean[i][k]        = errorFixed;
                            errorsMeanErr[i][k]     = errorFixed*0.01;
                            errorsMeanCorr[i][k]    = errorFixed;
                            errorsMeanErrCorr[i][k] = errorFixed*0.01;
                        }
                    }
                }
            }

            // fix Efficiency sys #12
            if (!nameCutVariationSC[i].CompareTo("Efficiency")||!nameCutVariationSC[i].CompareTo("8TeVEfficiency")){
                if (spectrumName.Contains("Ratio")){
                    if(!energy.CompareTo("2.76TeV")){
                        errorFixed          = TMath::Sqrt(2*2+1.5*1.5);
                    }else if(!energy.CompareTo("8TeV")){
                        errorFixed          = TMath::Sqrt(1.5*1.5+2*2+1.0*1.0);
                    }
                } else {
                    if(!energy.CompareTo("2.76TeV")){
                        errorFixed          = 1.5;
                    }else if(!energy.CompareTo("8TeV")){
                        errorFixed          = TMath::Sqrt(1.5*1.5+1.0*1.0);
                    }
                }
            }

            // 13 "Rapidity"
            if (!nameCutVariationSC[i].CompareTo("Rapidity")){
                if (spectrumName.Contains("Ratio")){
                    errorFixed          = 2.05;
                }
            }
            // 14 "OpeningAngle"
            if (!nameCutVariationSC[i].CompareTo("OpeningAngle")){
                if (spectrumName.Contains("Ratio")){
                    adjustPtDependent           = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                        if(!energy.CompareTo("2.76TeV")){
                          errorFixed          = 0.1+pow(ptBins[k],2)*0.004;
                        }else if(!energy.CompareTo("8TeV")){
                            errorFixed           = 0.25;
                            if(ptBins[k]>=6.) errorFixed += (0.15)*(ptBins[k]-6.);
                        }
                        if (errorFixed != -1){
                            errorsMean[i][k]        = errorFixed;
                            errorsMeanErr[i][k]     = errorFixed*0.01;
                            errorsMeanCorr[i][k]    = errorFixed;
                            errorsMeanErrCorr[i][k] = errorFixed*0.01;
                        }
                    }

                }
            }


            // put fixed values for pt independent errors, which were adjusted
            if (!adjustPtDependent && errorFixed != -1){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    errorsMean[i][k]        = errorFixed;
                    errorsMeanErr[i][k]     = errorFixed*0.01;
                    errorsMeanCorr[i][k]    = errorFixed;
                    errorsMeanErrCorr[i][k] = errorFixed*0.01;
                }
            }

        }


        // Add systematic error contribution from current cutvariation to total summed error
        for (Int_t l = 0; l < nPtBinsActive; l++){
            errorsPosSummed[l]                  = errorsPosSummed[l]+pow(errorsPos[i][l],2);
            errorsNegSummed[l]                  = errorsNegSummed[l]+ pow(errorsNeg[i][l],2);
            errorsMeanSummed[l]                 = errorsMeanSummed[l]+ pow(errorsMean[i][l],2);
            errorsPosCorrSummed[l]              = errorsPosCorrSummed[l]+pow(errorsPosCorr[i][l],2);
            errorsNegCorrSummed[l]              = errorsNegCorrSummed[l] +pow(errorsNegCorr[i][l],2);
            errorsMeanCorrSummed[l]             = errorsMeanCorrSummed[l]+ pow(errorsMeanCorr[i][l],2);
        }
        negativeErrors[i]                       = new TGraphErrors(nPtBinsActive,ptBins ,errorsNeg[i] ,ptBinsErr ,errorsNegErr[i] );
        meanErrors[i]                           = new TGraphErrors(nPtBinsActive,ptBins ,errorsMean[i] ,ptBinsErr ,errorsMeanErr[i] );
        positiveErrors[i]                       = new TGraphErrors(nPtBinsActive,ptBins ,errorsPos[i] ,ptBinsErr ,errorsPosErr[i] );
        negativeErrorsCorr[i]                   = new TGraphErrors(nPtBinsActive,ptBins ,errorsNegCorr[i] ,ptBinsErr ,errorsNegErrCorr[i] );
        meanErrorsCorr[i]                       = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorr[i] ,ptBinsErr ,errorsMeanErrCorr[i] );
        positiveErrorsCorr[i]                   = new TGraphErrors(nPtBinsActive,ptBins ,errorsPosCorr[i] ,ptBinsErr ,errorsPosErrCorr[i] );

    }

    for (Int_t l = 0; l < nPtBinsActive; l++){
        errorsPosSummed[l]                      = pow(errorsPosSummed[l],0.5);
        errorsMeanSummed[l]                     = pow(errorsMeanSummed[l],0.5);
        errorsPosErrSummed[l]                   = errorsPosSummed[l]*0.001;
        errorsMeanErrSummed[l]                  = errorsMeanSummed[l]*0.001;
        errorsNegSummed[l]                      = -pow(errorsNegSummed[l],0.5);
        errorsNegErrSummed[l]                   = errorsNegSummed[l]*0.001;
        errorsPosCorrSummed[l]                  = pow(errorsPosCorrSummed[l],0.5);
        errorsMeanCorrSummed[l]                 = pow(errorsMeanCorrSummed[l],0.5);
        errorsPosErrCorrSummed[l]               = errorsPosCorrSummed[l]*0.001;
        errorsMeanErrCorrSummed[l]              = errorsMeanCorrSummed[l]*0.001;
        errorsNegCorrSummed[l]                  = -pow(errorsNegCorrSummed[l],0.5);
        errorsNegErrCorrSummed[l]               = errorsNegCorrSummed[l]*0.001;
    }

    negativeErrorsSummed                        = new TGraphErrors(nPtBinsActive,ptBins ,errorsNegSummed ,ptBinsErr ,errorsNegErrSummed );
    negativeErrorsCorrSummed                    = new TGraphErrors(nPtBinsActive,ptBins ,errorsNegCorrSummed ,ptBinsErr ,errorsNegErrCorrSummed );
    positiveErrorsSummed                        = new TGraphErrors(nPtBinsActive,ptBins ,errorsPosSummed ,ptBinsErr ,errorsPosErrSummed );
    positiveErrorsCorrSummed                    = new TGraphErrors(nPtBinsActive,ptBins ,errorsPosCorrSummed ,ptBinsErr ,errorsPosErrCorrSummed );
    meanErrorsSummed                            = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanSummed ,ptBinsErr ,errorsMeanErrSummed );
    meanErrorsCorrSummed                        = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrSummed ,ptBinsErr ,errorsMeanErrCorrSummed );

    //++++++++++++++++++++++++++++++ PLOTTING OF SYSMEAN +++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TCanvas* canvasNewSysErrMean       = new TCanvas("canvasNewSysErrMean", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasNewSysErrMean,  0.075, 0.01, 0.015, 0.095);
    canvasNewSysErrMean->SetLogx(1);

    Double_t minXLegend     = 0.11;
    Double_t maxYLegend     = 0.95;
    Double_t widthLegend    = 0.23;
    Double_t heightLegend   = 1.12* 0.035*0.85 * (nCutsActive+1);
    Int_t nColumnsLegend    = 1;
    if ( nCutsActive+1 > 7){
        widthLegend         = 0.5;
        heightLegend        = 1.05* 0.035*0.85 * (nCutsActive/1.5+1);
        nColumnsLegend      = 2;
    }

    Double_t textSizeLabelsPixelmean        = 55;
    Double_t textSizeLabelsRelmean          = 55./1200;

        TH2F * histo2DNewSysErrMean;
        histo2DNewSysErrMean                = new TH2F("histo2DNewSysErrMean", "histo2DNewSysErrMean",1000, ptBins[0]-(ptBins[1]-ptBins[0])/2-0.07,1.2*(ptBins[nPtBinsActive-1]+ptBinsErr[nPtBinsActive-1]), 1000.,
                                                       yRangesSysPlotting[0], yRangesSysPlotting[1]);
        SetStyleHistoTH2ForGraphs( histo2DNewSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %",
                                0.85*textSizeLabelsRelmean, textSizeLabelsRelmean, 0.85*textSizeLabelsRelmean, textSizeLabelsRelmean, 0.9, 0.75);//(#times #epsilon_{pur})
        histo2DNewSysErrMean->GetYaxis()->SetLabelOffset(0.001);
        //histo2DNewSysErrMean->GetXaxis()->SetLabelOffset(-0.01);
        histo2DNewSysErrMean->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DNewSysErrMean->GetXaxis()->SetNoExponent();
        histo2DNewSysErrMean->DrawCopy();

        // create legend
        TLegend* legendMeanNew              = GetAndSetLegend2(minXLegend, maxYLegend-heightLegend, minXLegend+widthLegend, maxYLegend, 40*0.85, nColumnsLegend, "", 43, 0.1);
        for(Int_t i = 0;i< nCuts ;i++){
            if (benable[i]){
                DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], markerStyle[i], 1.,color[i],color[i]);
                meanErrorsCorr[i]->Draw("pX0,csame");
                legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");
            }
        }

        DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummed, 20, 1.,kBlack,kBlack);
        meanErrorsCorrSummed->Draw("p,csame");
        legendMeanNew->AddEntry(meanErrorsCorrSummed,"qd. sum","p");
        legendMeanNew->Draw();

        // labeling
        labelGamma->Draw();
        labelEnergy->Draw();
        labelSpectrum->Draw();

    canvasNewSysErrMean->Update();
    canvasNewSysErrMean->SaveAs(Form("GammaSystematicErrorsCalculated/SysMeanNewWithMeanEMC_%s_%s_%s.%s",spectrumName.Data(), energyForOutput.Data(), dateForOutput.Data(),suffix.Data()));

    //+++++++++++++++++++++++++ SAVING SYSTEMATICS TO DAT FILE +++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatname                   = Form("GammaSystematicErrorsCalculated/SystematicErrorEMC_%s_%s_%s.dat",spectrumName.Data(), energyForOutput.Data(), dateForOutput.Data());
    fstream SysErrDat;
    cout << SysErrDatname << endl;
    SysErrDat.open(SysErrDatname, ios::out);
    for (Int_t l=0; l< nPtBinsActive; l++){
        SysErrDat << ptBins[l] << "\t" <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
    }
    SysErrDat.close();

    //+++++++++++++++++++++++++ SAVING AVERAGE SYSTEMATICS TO DAT ++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatnameMean               = Form("GammaSystematicErrorsCalculated/SystematicErrorAveragedEMC_%s_%s_%s.dat", spectrumName.Data(), energyForOutput.Data(), dateForOutput.Data());
    fstream SysErrDatAver;
    cout << SysErrDatnameMean << endl;
    SysErrDatAver.open(SysErrDatnameMean, ios::out);
    for (Int_t l=0; l< nPtBinsActive; l++){
        SysErrDatAver  << ptBins[l] << "\t" << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
    }
    SysErrDatAver.close();

    //+++++++++++++++++++++++++ SAVING DETAILED SYSTEMATICS ++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatnameMeanSingleErr      = Form("GammaSystematicErrorsCalculated/SystematicErrorAveragedSingleEMC_%s_%s_%s.dat", spectrumName.Data(),energyForOutput.Data(), dateForOutput.Data());
    fstream SysErrDatAverSingle;
    cout << SysErrDatnameMeanSingleErr << endl;
    SysErrDatAverSingle.open(SysErrDatnameMeanSingleErr, ios::out);
    SysErrDatAverSingle << "Pt bin\t" ;
    for (Int_t i= 0; i< numberCutStudies; i++){
        if (benable[i]) SysErrDatAverSingle << nameCutVariationSC[i] << "\t";
    }
    SysErrDatAverSingle << endl;
    for (Int_t l=0;l< nPtBinsActive;l++){
        SysErrDatAverSingle << ptBins[l] << "\t";
        for (Int_t i= 0; i< numberCutStudies; i++){
            if (benable[i]) SysErrDatAverSingle << errorsMeanCorr[i][l] << "\t";
        }
        SysErrDatAverSingle << errorsMeanCorrSummed[l] << endl;
    }
    SysErrDatAverSingle.close();

    Double_t errorsMeanCorrSignalExtraction[nPtBinsActive];
    Double_t errorsMeanCorrClusterProp[nPtBinsActive];

    for (Int_t l=0; l< nPtBinsActive; l++){
        // 0 "ClusterNonLinearity"
        // 1 "ClusterM02"
        // 2 "ClusterMinEnergy"
        // 3 "ClusterTrackMatching"
        // 4 "ClusterTiming"
        // 5 "ClusterNCells"
        // 6 "ClusterMaterialTRD"
        // 7 "SPD"
        // 8 "BG"
        // 9 "Alpha"
        // 10 "Cocktail"
        // 11 "IntRange"
        // 12 "Efficiency"
        // 13 "Rapidity"
        // 14 "OpeningAngle"

        // Signal extraction: Pi0Yield extraction, Alpha ,BG, Openangle, Rapidity
        errorsMeanCorrSignalExtraction[l]       =   TMath::Sqrt(pow(errorsMeanCorr[9][l],2)+   // Alpha
                                                                pow(errorsMeanCorr[8][l],2)+   // BG
                                                                pow(errorsMeanCorr[13][l],2)+  // Rapidity
                                                                pow(errorsMeanCorr[14][l],2)+  // OpeningAngle
                                                                pow(errorsMeanCorr[11][l],2)); // Pi0 Yield Extraction


        errorsMeanCorrClusterProp[l]            =   TMath::Sqrt(pow(errorsMeanCorr[2][l],2)+    // ClusterMinEnergy
                                                                pow(errorsMeanCorr[5][l],2)+    // ClusterNCells
                                                                pow(errorsMeanCorr[1][l],2)+    // ClusterM02
                                                                pow(errorsMeanCorr[4][l],2));   // ClusterTiming

    }
    TGraphErrors* meanErrorsSignalExtraction    = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrSignalExtraction ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsClusterProp         = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrClusterProp ,ptBinsErr ,errorsMeanErrCorrSummed );

    //++++++++++++++++++++++++ PLOTTING OF SYSERRORSUMMEDVISU ++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    Double_t textSizeLabelsPixel        = 55;
    Double_t textSizeLabelsRel          = 55./1200;

    TCanvas* canvasSummedErrMean        = new TCanvas("canvasSummedErrMean", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasSummedErrMean,  0.075, 0.01, 0.015, 0.095);
    canvasSummedErrMean->SetLogx(1);

        TH2F * histo2DSummedErrMean;
        histo2DSummedErrMean            = new TH2F("histo2DSummedErrMean", "histo2DSummedErrMean",1000, ptBins[0]-(ptBins[1]-ptBins[0])/2-0.07, 1.2*(ptBins[nPtBinsActive-1]+ptBinsErr[nPtBinsActive-1]),
                                                   1000.,yRangesSysPlotting[0],yRangesSysPlotting[1] );
        SetStyleHistoTH2ForGraphs(  histo2DSummedErrMean, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %",
                                    0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 0.75);//(#times #epsilon_{pur})
        histo2DSummedErrMean->GetYaxis()->SetLabelOffset(0.001);
        //histo2DSummedErrMean->GetXaxis()->SetLabelOffset(-0.01);
        histo2DSummedErrMean->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DSummedErrMean->GetXaxis()->SetNoExponent();
        histo2DSummedErrMean->DrawCopy();

        Double_t minXLegend2        = 0.11;
        Double_t maxYLegend2        = 0.95;
        Double_t widthLegend2       = 0.5;
        Double_t heightLegend2      = 4*0.04;

        // create legend
        TLegend* legendSummedMeanNew    = GetAndSetLegend2(minXLegend2, maxYLegend2-heightLegend2, minXLegend2+widthLegend2, maxYLegend2, 40, 2, "", 43, 0.1);
        Size_t markersizeSummed     = 1.3;
        // 0 "ClusterNonLinearity"
        // 1 "ClusterM02"
        // 2 "ClusterMinEnergy"
        // 3 "ClusterTrackMatching"
        // 4 "ClusterTiming"
        // 5 "ClusterNCells"
        // 6 "ClusterMaterialTRD"
        // 7 "SPD"
        // 8 "BG"
        // 9 "Alpha"
        // 10 "Cocktail"
        // 11 "IntRange"
        // 12 "Efficiency"
        // 13 "Rapidity"
        // 14 "OpeningAngle"

        // Signal extraction error
        if (benable[8] || benable[9] || benable[11] || benable[13] || benable[14]){
            DrawGammaSetMarkerTGraphErr(meanErrorsSignalExtraction, GetMarkerStyleSystematics("IntRange"), markersizeSummed, GetColorSystematics("IntRange"),GetColorSystematics("IntRange"));
            meanErrorsSignalExtraction->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"signal ext. #pi^{0}","p");
        }
        if (benable[7] ){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[7], GetMarkerStyleSystematics("Pileup"), markersizeSummed, GetColorSystematics("Pileup"),GetColorSystematics("Pileup"));
            meanErrorsCorr[7]->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsCorr[7],"pileup","p");
        }
        if (benable[1] || benable[2] || benable[4] || benable[5] ){
            DrawGammaSetMarkerTGraphErr(meanErrorsClusterProp, GetMarkerStyleSystematics("ClusterM02"), markersizeSummed, GetColorSystematics("ClusterM02"),GetColorSystematics("ClusterM02"));
            meanErrorsClusterProp->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsClusterProp,"cluster prop.","p");
        }
        if (benable[0]){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[0], markerStyle[0], markersizeSummed,color[0],color[0]);
            meanErrorsCorr[0]->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsCorr[0],"cl. energy scale","p");
        }
        if (benable[3]){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[3], markerStyle[3], markersizeSummed,color[3],color[3]);
            meanErrorsCorr[3]->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsCorr[3],"track match. cl.","p");
        }
        if (benable[12]){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[12], markerStyle[12], markersizeSummed,color[12],color[12]);
            meanErrorsCorr[12]->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsCorr[12],"efficiency","p");
        }
        if (benable[10]){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[10], markerStyle[10], markersizeSummed,color[10],color[10]);
            meanErrorsCorr[10]->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsCorr[10],"cocktail","p");
        }
        if (benable[6]){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[6], markerStyle[6], markersizeSummed,color[6],color[6]);
            meanErrorsCorr[6]->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsCorr[6],"outer material","p");
        }

        DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummed, 20, markersizeSummed,kBlack,kBlack);
        meanErrorsCorrSummed->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsCorrSummed,"qd. sum","p");
        legendSummedMeanNew->Draw();

        labelGamma->Draw();
        labelEnergy->Draw();
        labelSpectrum->Draw();

    canvasSummedErrMean->Update();
    canvasSummedErrMean->SaveAs(Form("GammaSystematicErrorsCalculated/SysErrorSummedVisuEMC_%s_%s_%s.%s",spectrumName.Data(), energyForOutput.Data(), dateForOutput.Data(),suffix.Data()));
    delete canvasSummedErrMean;


    //+++++++++++++++++++++++++ SAVING DETAILED SYSTEMATICS ++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatnameMeanSinglePaperErr      = Form("GammaSystematicErrorsCalculated/SystematicErrorAveragedSinglePaperEMC_%s_%s_%s.dat", spectrumName.Data(),
                                                            energyForOutput.Data(), dateForOutput.Data());
    fstream SysErrDatAverSinglePaper;
    cout << SysErrDatnameMeanSinglePaperErr << endl;
    SysErrDatAverSinglePaper.open(SysErrDatnameMeanSinglePaperErr, ios::out);
    SysErrDatAverSinglePaper << "Pt bin\t" ;

    if (benable[8] || benable[9] || benable[11] || benable[13] || benable[14])  SysErrDatAverSinglePaper << "SignalExtraction" << "\t";
    if (benable[7])                                                             SysErrDatAverSinglePaper << "Pileup" << "\t";
    if (benable[1] || benable[2] || benable[4] || benable[5] )                  SysErrDatAverSinglePaper << "ClusterProp" << "\t";
    if (benable[0])                                                             SysErrDatAverSinglePaper << "ClusterEnergyScale" << "\t";
    if (benable[3])                                                             SysErrDatAverSinglePaper << "ClusterTM" << "\t";
    if (benable[12])                                                            SysErrDatAverSinglePaper << "Efficiency" << "\t";
    if (benable[10])                                                            SysErrDatAverSinglePaper << "Cocktail" << "\t";
    if (benable[6])                                                             SysErrDatAverSinglePaper << "OuterMaterial" << "\t";
    SysErrDatAverSinglePaper << endl;

    for (Int_t l=0;l< nPtBinsActive;l++){
        SysErrDatAverSinglePaper << ptBins[l] << "\t";
        if (benable[8] || benable[9] || benable[11] || benable[13] || benable[14])  SysErrDatAverSinglePaper << errorsMeanCorrSignalExtraction[l] << "\t";
        if (benable[7])                                                             SysErrDatAverSinglePaper << errorsMeanCorr[7][l] << "\t";
        if (benable[1] || benable[2] || benable[4] || benable[5] )                  SysErrDatAverSinglePaper << errorsMeanCorrClusterProp[l] << "\t";
        if (benable[0])                                                             SysErrDatAverSinglePaper << errorsMeanCorr[0][l] << "\t";
        if (benable[3])                                                             SysErrDatAverSinglePaper << errorsMeanCorr[3][l] << "\t";
        if (benable[12])                                                            SysErrDatAverSinglePaper << errorsMeanCorr[12][l] << "\t";
        if (benable[10])                                                            SysErrDatAverSinglePaper << errorsMeanCorr[10][l] << "\t";
        if (benable[6])                                                             SysErrDatAverSinglePaper << errorsMeanCorr[6][l] << "\t";
        SysErrDatAverSinglePaper << errorsMeanCorrSummed[l] << endl;
    }
    SysErrDatAverSinglePaper.close();

}
