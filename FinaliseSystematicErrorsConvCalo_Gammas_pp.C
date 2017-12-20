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
//*********** Main function to calculate gamma errors for conversion measurements in pp collisions ************
//*************************************************************************************************************
void FinaliseSystematicErrorsConvCalo_Gammas_pp (   TString nameDataFileErrors      = "",
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
    TString detectionProcess                    = ReturnFullTextReconstructionProcess(0);
    TString detectionProcessPi0                 = ReturnTextReconstructionProcess(41);


    TString dateForOutput                       = ReturnDateStringForOutput();
    TString energyForOutput                     = ReturnCollisionEnergyOutputString(energy);

    TLatex *labelGamma                          = new TLatex(0.95,0.88,detectionProcess);
    SetStyleTLatex( labelGamma, 0.038,4);
    labelGamma->SetTextAlign(31);
    TLatex *labelPi0                            = new TLatex(0.95,0.84,detectionProcessPi0);
    SetStyleTLatex( labelPi0, 0.038,4);
    labelPi0->SetTextAlign(31);

    Double_t  posYLabel                         = 0.84;
    if(!spectrumName.CompareTo("IncRatio") || !spectrumName.CompareTo("DoubleRatio"))
        posYLabel                               = posYLabel-0.04;

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
    if (energy.CompareTo("pPb_5.023TeV") == 0 ){
        if (spectrumName.Contains("Gamma"))
            yRangesSysPlotting[1]               = 15.5;
        else if (spectrumName.Contains("IncRatio"))
            yRangesSysPlotting[1]               = 20.5;
        else
            yRangesSysPlotting[1]               = 20.5;
    }

    // Set names of cut variations for file input
    TString nameCutVariationSC[24]              = { "dEdxE", "dEdxPi","TPCCluster", "SinglePt", "Chi2",
                                                    "Qt" , "CosPoint", "ToCloseV0s", "ConvPhi", "ClusterNonLinearity",
                                                    "ClusterM02", "ClusterMinEnergy", "ClusterTrackMatching", "ClusterTiming", "ClusterNCells",
                                                    "ClusterMaterialTRD",  "OOBPileupGamma", "SPD", "BG","Alpha" ,
                                                    "Cocktail", "Material", "IntRange", "Efficiency"
                                                  };

    // Set colors and markers
    Color_t color[24];
    Color_t markerStyle[24];
    // Set names of cut variations for legends
    TString nameCutVariation[24];

    for (Int_t k =0; k<nCuts; k++ ){
        cout << "variation: " << nameCutVariationSC[k].Data() << endl;
        color[k]                                = GetColorSystematics( nameCutVariationSC[k] );
        markerStyle[k]                          = GetMarkerStyleSystematics( nameCutVariationSC[k] );
        nameCutVariation[k]                     = GetSystematicsName(nameCutVariationSC[k]);
        cout << "name for writing: " << nameCutVariation[k].Data() << "\t"<< color[k]  << "\t" << markerStyle[k] << endl;
    }

    // Create output folder
    gSystem->Exec("mkdir -p GammaSystematicErrorsCalculated");

    // ***************************************************************************************************
    // ******************************** Booleans for enabling systematic errors **************************
    // ***************************************************************************************************
    Bool_t benable[24]                          = { 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                    0, 0, 0, 0 };
    Bool_t benableIncGamma2760GeV[24]           = { 1, 1, 1, 1, 1,  1, 0, 0, 1, 0,
                                                    0, 0, 0, 0, 0,  0, 1, 1, 0, 0,
                                                    0, 1, 0, 0 };
    Bool_t benableIncRatio2760GeV[24]           = { 1, 1, 1, 1, 1,  1, 0, 0, 1, 1,
                                                    1, 1, 1, 1, 1,  1, 1, 0, 0, 1,
                                                    0, 0, 1, 1 };
    Bool_t benableDR2760GeV[24]                 = { 1, 1, 1, 1, 1,  1, 0, 0, 1, 1,
                                                    1, 1, 1, 1, 1,  1, 1, 0, 0, 1,
                                                    1, 0, 1, 1 };

    Bool_t benableIncGamma8TeV[24]              = { 1, 1, 1, 1, 1,  1, 0, 0, 1, 0,
                                                    0, 0, 0, 0, 0,  0, 1, 1, 0, 0,
                                                    1, 1, 0, 1 };
    Bool_t benableIncRatio8TeV[24]              = { 1, 1, 1, 1, 1,  1, 0, 0, 1, 1,
                                                    1, 1, 1, 1, 1,  1, 1, 1, 0, 1,
                                                    1, 0, 1, 1 };
    Bool_t benableDR8TeV[24]                    = { 1, 1, 1, 1, 1,  1, 0, 0, 1, 1,
                                                    1, 1, 1, 1, 1,  1, 1, 1, 0, 1,
                                                    1, 0, 1, 1 };

    // ***************************************************************************************************
    // ******************************** Booleans for smoothing *******************************************
    // ***************************************************************************************************
    Bool_t bsmooth[24]                          = { 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                    0, 0, 0, 0 };
    Bool_t bsmoothIncGamma2760GeV[24]           = { 1, 1, 0, 0, 1,  1, 0, 0, 1, 0,
                                                    0, 0, 0, 0, 0,  0, 1, 1, 0, 0,
                                                    0, 1, 0, 1 };
    Bool_t bsmoothIncRatio2760GeV[24]           = { 1, 1, 0, 1, 1,  1, 0, 0, 1, 1,
                                                    1, 1, 1, 1, 1,  1, 1, 1, 0, 1,
                                                    0, 0, 0, 1 };
    Bool_t bsmoothDR2760GeV[24]                 = { 1, 1, 0, 1, 1,  1, 0, 0, 1, 1,
                                                    1, 1, 1, 1, 1,  1, 1, 1, 0, 1,
                                                    0, 0, 0, 1 };

    Bool_t bsmoothIncGamma8TeV[24]              = { 1, 1, 1, 1, 1,  1, 0, 0, 1, 0,
                                                    0, 0, 0, 0, 0,  0, 1, 1, 0, 0,
                                                    0, 1, 0, 1 };
    Bool_t bsmoothIncRatio8TeV[24]              = { 1, 1, 1, 1, 1,  1, 0, 0, 1, 1,
                                                    1, 1, 1, 1, 1,  1, 1, 1, 0, 1,
                                                    0, 0, 1, 1 };
    Bool_t bsmoothDR8TeV[24]                    = { 1, 1, 1, 1, 1,  1, 0, 0, 1, 1,
                                                    1, 1, 1, 1, 1,  1, 1, 1, 0, 1,
                                                    1, 0, 1, 1 };

    for (Int_t i = 0; i < numberCutStudies; i++){
        if (energy.CompareTo("2.76TeV") == 0){
            if(!spectrumName.CompareTo("IncRatio")){
                bsmooth[i]                      = bsmoothIncRatio2760GeV[i];
                benable[i]                      = benableIncRatio2760GeV[i];
            } else if(!spectrumName.CompareTo("DoubleRatio")){
                bsmooth[i]                      = bsmoothDR2760GeV[i];
                benable[i]                      = benableDR2760GeV[i];
            } else if(!spectrumName.CompareTo("Gamma")){
                bsmooth[i]                      = bsmoothIncGamma2760GeV[i];
                benable[i]                      = benableIncGamma2760GeV[i];
            }
        } else if (energy.CompareTo("8TeV") == 0){
          if(!spectrumName.CompareTo("IncRatio")){
              bsmooth[i]                      = bsmoothIncRatio8TeV[i];
              benable[i]                      = benableIncRatio8TeV[i];
          } else if(!spectrumName.CompareTo("DoubleRatio")){
              bsmooth[i]                      = bsmoothDR8TeV[i];
              benable[i]                      = benableDR8TeV[i];
          } else if(!spectrumName.CompareTo("Gamma")){
              bsmooth[i]                      = bsmoothIncGamma8TeV[i];
              benable[i]                      = benableIncGamma8TeV[i];
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
        if ( nameCutVariationSC[i].CompareTo("SPD")==0  || nameCutVariationSC[i].CompareTo("Material")==0 || nameCutVariationSC[i].CompareTo("Efficiency")==0){
            TString nameGraphPos;
            TString nameGraphNeg;
            nameGraphPos                        = Form("%s_SystErrorRelPos_%s_pp",spectrumName.Data(),nameCutVariationSC[0].Data());
            nameGraphNeg                        = Form("%s_SystErrorRelNeg_%s_pp",spectrumName.Data(),nameCutVariationSC[0].Data());
            cout << "Cutstudies " << i << "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<< "\t fixing" <<  endl;
            graphPosErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else if ( nameCutVariationSC[i].CompareTo("OOBPileupGamma")==0 ){
            TString nameGraphPos                = Form("Gamma_SystErrorRelPos_%s_pp",nameCutVariationSC[i].Data());
            TString nameGraphNeg                = Form("Gamma_SystErrorRelNeg_%s_pp",nameCutVariationSC[i].Data());
            cout << "Cutstudies " << i << "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
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

            // fix dEdx e error #0
            if (!nameCutVariationSC[i].CompareTo("dEdxE")){
                adjustPtDependent               = kTRUE;
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if (spectrumName.Contains("Ratio")){
                      if(!energy.CompareTo("2.76TeV")){
                        if (ptBins[k]>1.5){
                            errorFixed      = 0.2+pow(ptBins[k],1.9)*0.04;
                        }
                      }else if(!energy.CompareTo("8TeV")){
                        errorFixed = 0.1;
                      }
                    } else {
                      if(!energy.CompareTo("2.76TeV")){
                        if (ptBins[k]>2.0)
                            errorFixed      = 0.05+pow(ptBins[k],1.8)*0.04;
                      }else if(!energy.CompareTo("8TeV")){
                            errorFixed = 0.05;
                      }
                    }

                    if (errorFixed != -1){
                        errorsMean[i][k]        = errorFixed;
                        errorsMeanErr[i][k]     = errorFixed*0.01;
                        errorsMeanCorr[i][k]    = errorFixed;
                        errorsMeanErrCorr[i][k] = errorFixed*0.01;
                    }
                }
            }

            // fix dEdxPi sys #1
            if (!nameCutVariationSC[i].CompareTo("dEdxPi")){
                if (spectrumName.Contains("Ratio")){
                    adjustPtDependent               = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                        if(!energy.CompareTo("2.76TeV")){
                          errorFixed                  = 0.15+pow(ptBins[k],2.6)*0.02+0.01*ptBins[k];
                        }else if(!energy.CompareTo("8TeV")){
                          if(ptBins[k]>=1.8) errorFixed = 0.8;
                        }
                        if (errorFixed != -1){
                            errorsMean[i][k]        = errorFixed;
                            errorsMeanErr[i][k]     = errorFixed*0.01;
                            errorsMeanCorr[i][k]    = errorFixed;
                            errorsMeanErrCorr[i][k] = errorFixed*0.01;
                        }
                    }
                } else {
                    adjustPtDependent               = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                          if(!energy.CompareTo("2.76TeV")){
                            if (ptBins[k] > 1.6){
                                errorFixed                  = 0.3+pow(ptBins[k],2.2)*0.012+0.001*ptBins[k];
                            }
                          }else if(!energy.CompareTo("8TeV")){
                            if(ptBins[k]>=1.8) errorFixed = 0.6;
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

            // fix Single pt sys #3
            if (!nameCutVariationSC[i].CompareTo("TPCCluster")){
                if (spectrumName.Contains("Ratio")){
                    adjustPtDependent               = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                        if(!energy.CompareTo("2.76TeV")){
                          errorFixed          = 0.05+pow(ptBins[k],1.5)*0.02+0.1*ptBins[k];
                        }else if(!energy.CompareTo("8TeV")){
                          errorFixed = 0.05;
                        }

                        if (errorFixed != -1){
                            errorsMean[i][k]        = errorFixed;
                            errorsMeanErr[i][k]     = errorFixed*0.01;
                            errorsMeanCorr[i][k]    = errorFixed;
                            errorsMeanErrCorr[i][k] = errorFixed*0.01;
                        }
                    }
                } else {
                  if(!energy.CompareTo("2.76TeV")){
                    errorFixed                  = 0.05;
                  }else if(!energy.CompareTo("8TeV")){
                    errorFixed = 0.05;
                  }
                }
            }

            // fix Single pt sys #2
            if (!nameCutVariationSC[i].CompareTo("SinglePt")){
                if (spectrumName.Contains("Ratio")){
                    adjustPtDependent               = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                        if(!energy.CompareTo("2.76TeV")){
                          if (ptBins[k] > 1.2)
                              errorFixed              = 1.3;
                        }else if(!energy.CompareTo("8TeV")){
                              errorFixed = 0.2;
                        }
                        if (errorFixed != -1){
                            errorsMean[i][k]        = errorFixed;
                            errorsMeanErr[i][k]     = errorFixed*0.01;
                            errorsMeanCorr[i][k]    = errorFixed;
                            errorsMeanErrCorr[i][k] = errorFixed*0.01;
                        }
                    }
                } else {
                  if(!energy.CompareTo("2.76TeV")){
                    errorFixed                  = 0.15;
                  }else if(!energy.CompareTo("8TeV")){
                    errorFixed                  = 0.1;
                  }
                }
            }

            // fix Chi2/psi pair uncertainties #4
            if (!nameCutVariationSC[i].CompareTo("Chi2")){
                if (spectrumName.Contains("Ratio")){
                    adjustPtDependent           = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                        if(!energy.CompareTo("2.76TeV")){
                          if (ptBins[k] > 1.2)
                              errorFixed          = 0.5+pow(ptBins[k],2.1)*0.012+0.02*ptBins[k];
                        }else if(!energy.CompareTo("8TeV")){
                          errorFixed              = 0.8;
                          if(ptBins[k]>4.) errorFixed += 0.1*(ptBins[k]-4.);
                        }
                        if (errorFixed != -1){
                            errorsMean[i][k]        = errorFixed;
                            errorsMeanErr[i][k]     = errorFixed*0.01;
                            errorsMeanCorr[i][k]    = errorFixed;
                            errorsMeanErrCorr[i][k] = errorFixed*0.01;
                        }
                    }
                } else {
                    adjustPtDependent           = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                        if(!energy.CompareTo("2.76TeV")){
                          errorFixed              = 0.2+pow(ptBins[k],1.8)*0.07+0.03*ptBins[k];
                        }else if(!energy.CompareTo("8TeV")){
                          errorFixed              = 1.2;
                          if(ptBins[k]>4.) errorFixed += 0.05*(ptBins[k]-4.);
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

            // fix Qt sys #5
            if (!nameCutVariationSC[i].CompareTo("Qt")){
                adjustPtDependent               = kTRUE;
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if (spectrumName.Contains("Ratio")){
                      if(!energy.CompareTo("2.76TeV")){
                        errorFixed      = 0.25+pow(ptBins[k],1.5)*0.08+0.02*ptBins[k]*ptBins[k];
                      }else if(!energy.CompareTo("8TeV")){
                        errorFixed              = 0.1;
                        if(ptBins[k]>3.) errorFixed += 0.01*(ptBins[k]-3.)*(ptBins[k]-3.) + 0.02*(ptBins[k]-3.);
                      }
                    } else {
                      if(!energy.CompareTo("2.76TeV")){
                        errorFixed      = 0.15+pow(ptBins[k],2.3)*0.02+0.04*ptBins[k];
                      }else if(!energy.CompareTo("8TeV")){
                        errorFixed              = 0.1;
                        if(ptBins[k]>3.) errorFixed += 0.02*(ptBins[k]-3.)*(ptBins[k]-3.) + 0.04*(ptBins[k]-3.);
                      }
                    }
                    if (errorFixed != -1){
                        errorsMean[i][k]        = errorFixed;
                        errorsMeanErr[i][k]     = errorFixed*0.01;
                        errorsMeanCorr[i][k]    = errorFixed;
                        errorsMeanErrCorr[i][k] = errorFixed*0.01;
                    }
                }
            }

            // 6 "CosPoint"

            // 7 "ToCloseV0s"

            // fix ConvPhi pileup sys #8
            if (!nameCutVariationSC[i].CompareTo("ConvPhi")){
              if(!energy.CompareTo("2.76TeV")){
                errorFixed                  = 0.8;
              }else if(!energy.CompareTo("8TeV")){
                errorFixed                  = 0.5;
              }
            }

            // 9 "ClusterNonLinearity"
            if (!nameCutVariationSC[i].CompareTo("ClusterNonLinearity")){
                adjustPtDependent           = kTRUE;
                for (Int_t k = 0; k < nPtBins; k++){
                    if(!energy.CompareTo("2.76TeV")){
                      errorFixed              = 2.+(0.01)*ptBins[k]+(0.01)*ptBins[k]*ptBins[k];
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
            }

            // 10 "ClusterM02"
            if (!nameCutVariationSC[i].CompareTo("ClusterM02")){
                adjustPtDependent           = kTRUE;
                for (Int_t k = 0; k < nPtBins; k++){
                    if(!energy.CompareTo("2.76TeV")){
                      errorFixed              =  1.8+(0.01)*ptBins[k]+(0.07)*ptBins[k]*ptBins[k];
                    }else if(!energy.CompareTo("8TeV")){
                      errorFixed              =  2.0;
                      if(ptBins[k]>3.) errorFixed += (0.1)*(ptBins[k]-3.);
                    }
                    if (errorFixed != -1){
                        errorsMean[i][k]            = errorFixed;
                        errorsMeanErr[i][k]         = errorFixed*0.01;
                        errorsMeanCorr[i][k]        = errorFixed;
                        errorsMeanErrCorr[i][k]     = errorFixed*0.01;
                    }
                }

            }
            // 11 "ClusterMinEnergy"
            if (!nameCutVariationSC[i].CompareTo("ClusterMinEnergy")){
              adjustPtDependent           = kTRUE;
              for (Int_t k = 0; k < nPtBins; k++){
                if(!energy.CompareTo("2.76TeV")){
                  errorFixed                  = 1.4;
                }else if(!energy.CompareTo("8TeV")){
                  errorFixed                  = 0.75+20/pow(10,ptBins[k]);;
                }
                if (errorFixed != -1){
                    errorsMean[i][k]            = errorFixed;
                    errorsMeanErr[i][k]         = errorFixed*0.01;
                    errorsMeanCorr[i][k]        = errorFixed;
                    errorsMeanErrCorr[i][k]     = errorFixed*0.01;
                }
              }
            }

            // 12 "ClusterTrackMatching"
            if (!nameCutVariationSC[i].CompareTo("ClusterTrackMatching")){
                adjustPtDependent           = kTRUE;
                for (Int_t k = 0; k < nPtBins; k++){
                    if(!energy.CompareTo("2.76TeV")){
                      if (ptBins[k] > 1.4)
                          errorFixed              =  1.2+(0.014)*ptBins[k]+(0.09)*ptBins[k]*ptBins[k]; // parametrisation
                    }else if(!energy.CompareTo("8TeV")){
                          errorFixed = 0.2;
                         if(ptBins[k]>=4) errorFixed += 2E-3*(ptBins[k]-4)*(ptBins[k]-4)+0.25*(ptBins[k]-4);
                    }
                    if (errorFixed != -1){
                        errorsMean[i][k]            = errorFixed;
                        errorsMeanErr[i][k]         = errorFixed*0.01;
                        errorsMeanCorr[i][k]        = errorFixed;
                        errorsMeanErrCorr[i][k]     = errorFixed*0.01;
                    }
                }
            }

            // fix ClusterTiming pileup sys #13
            if (!nameCutVariationSC[i].CompareTo("ClusterTiming")){
              if(!energy.CompareTo("2.76TeV")){
                errorFixed                  = 0.6;
              }else if(!energy.CompareTo("8TeV")){
                errorFixed                  = 0.7;
              }
            }

            // 14 "ClusterNCells"
            if (!nameCutVariationSC[i].CompareTo("ClusterNCells")){
              if(!energy.CompareTo("2.76TeV")){
                errorFixed                  = 0.5;
              }else if(!energy.CompareTo("8TeV")){
                errorFixed                  = 0.5;
              }
            }


            // fix ClusterMaterialTRD sys #15
            if (!nameCutVariationSC[i].CompareTo("ClusterMaterialTRD")){
              if(!energy.CompareTo("2.76TeV")){
                errorFixed                  = 2.1;
              }else if(!energy.CompareTo("8TeV")){
                errorFixed                  = 2.1;
              }
            }

            // fix out-of-bunch gamma sys #16
            if (!nameCutVariationSC[i].CompareTo("OOBPileupGamma")){
                adjustPtDependent           = kTRUE;
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if(!energy.CompareTo("2.76TeV")){
                      if (ptBins[k] > 1.5)
                          errorFixed              = 0.9;
                    }else if(!energy.CompareTo("8TeV")){
                      errorFixed                  = 2.25;
                      if(ptBins[k] < 3.) errorFixed += pow(ptBins[k]-3,2)*0.25;
                      if(ptBins[k] > 3.) errorFixed += 0.1*(ptBins[k]-3.);
                    }
                    if (errorFixed != -1){
                        errorsMean[i][k]        = errorFixed;
                        errorsMeanErr[i][k]     = errorFixed*0.01;
                        errorsMeanCorr[i][k]    = errorFixed;
                        errorsMeanErrCorr[i][k] = errorFixed*0.01;
                    }
                }
            }

            // fix SPD pileup sys #17
            if (!nameCutVariationSC[i].CompareTo("SPD")){
              if(!energy.CompareTo("2.76TeV")){
                errorFixed                  = 0.25;
              }else if(!energy.CompareTo("8TeV")){
                errorFixed                  = 0.1;
              }
            }

            // fix BG sys #18
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


            // fix alpha uncertainties #19
            if (!nameCutVariationSC[i].CompareTo("Alpha")){
              if(!energy.CompareTo("2.76TeV")){
                errorFixed                  = 0.25;
              }else if(!energy.CompareTo("8TeV")){
                errorFixed                  = 0.5;
              }
            }
            // fix "Cocktail" #20
            if (!nameCutVariationSC[i].CompareTo("Cocktail")){
              if (!spectrumName.CompareTo("DoubleRatio")){
                adjustPtDependent           = kTRUE;
                for (Int_t k = 0; k < nPtBinsActive; k++){
                  if(!energy.CompareTo("8TeV")){
                    if(ptBins[k]<=3.) errorFixed = 1. ;
                    else              errorFixed = -1;
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
            // fix Material sys #21
            if (!nameCutVariationSC[i].CompareTo("Material")){
                if (!spectrumName.Contains("Ratio")){
                    errorFixed                  = 4.5;
                } else {
                    errorFixed                  = 0.;
                }
            }

            // fix IntRange sys #22
            if (!nameCutVariationSC[i].CompareTo("IntRange")){
                if (spectrumName.Contains("Ratio")){
                    adjustPtDependent           = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                        if(!energy.CompareTo("2.76TeV")){
                          if (ptBins[k] > 1.6)
                              errorFixed          = 1.6+pow(ptBins[k],2)*0.018;
                        }else if(!energy.CompareTo("8TeV")){
                              errorFixed          = 1.75;
                              if(ptBins[k]>=2.) errorFixed         += 0.15*(ptBins[k]-2) ;
                              if(ptBins[k]==0.9) errorFixed = 3.;
                              else if(ptBins[k]==1.1) errorFixed = 2.;
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
            // fix Efficiency uncertainties #23
            if (!nameCutVariationSC[i].CompareTo("Efficiency")){
              if(!energy.CompareTo("2.76TeV")){
                errorFixed                  = 2.0;
              }else if(!energy.CompareTo("8TeV")){
                if(spectrumName.Contains("Ratio")){
                   errorFixed                  = TMath::Sqrt(2.0*2.0 + 0.5*0.5); // pi0_eff - gamma_eff
                }else{
                   errorFixed                  = TMath::Sqrt(0.5*0.5); // gamma_eff
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
        if (spectrumName.Contains("Ratio")) labelPi0->Draw();
        labelEnergy->Draw();
        labelSpectrum->Draw();

    canvasNewSysErrMean->Update();
    canvasNewSysErrMean->SaveAs(Form("GammaSystematicErrorsCalculated/SysMeanNewWithMeanPCMEMC_%s_%s_%s.%s",spectrumName.Data(), energyForOutput.Data(), dateForOutput.Data(),suffix.Data()));

    //+++++++++++++++++++++++++ SAVING SYSTEMATICS TO DAT FILE +++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatname                   = Form("GammaSystematicErrorsCalculated/SystematicErrorPCMEMC_%s_%s_%s.dat",spectrumName.Data(), energyForOutput.Data(), dateForOutput.Data());
    fstream SysErrDat;
    cout << SysErrDatname << endl;
    SysErrDat.open(SysErrDatname, ios::out);
    for (Int_t l=0; l< nPtBinsActive; l++){
        SysErrDat << ptBins[l] << "\t" <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
    }
    SysErrDat.close();

    //+++++++++++++++++++++++++ SAVING AVERAGE SYSTEMATICS TO DAT ++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatnameMean               = Form("GammaSystematicErrorsCalculated/SystematicErrorAveragedPCMEMC_%s_%s_%s.dat", spectrumName.Data(), energyForOutput.Data(), dateForOutput.Data());
    fstream SysErrDatAver;
    cout << SysErrDatnameMean << endl;
    SysErrDatAver.open(SysErrDatnameMean, ios::out);
    for (Int_t l=0; l< nPtBinsActive; l++){
        SysErrDatAver  << ptBins[l] << "\t" << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
    }
    SysErrDatAver.close();

    //+++++++++++++++++++++++++ SAVING DETAILED SYSTEMATICS ++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatnameMeanSingleErr      = Form("GammaSystematicErrorsCalculated/SystematicErrorAveragedSinglePCMEMC_%s_%s_%s.dat", spectrumName.Data(),energyForOutput.Data(), dateForOutput.Data());
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

    Double_t errorsMeanCorrPID[nPtBinsActive];
    Double_t errorsMeanCorrSignalExtraction[nPtBinsActive];
    Double_t errorsMeanCorrPileup[nPtBinsActive];
    Double_t errorsMeanCorrTrackReco[nPtBinsActive];
    Double_t errorsMeanCorrPhotonReco[nPtBinsActive];
    Double_t errorsMeanCorrClusterProp[nPtBinsActive];

    for (Int_t l=0; l< nPtBinsActive; l++){
        // 0 "dEdxE"
        // 1 "dEdxPi"
        // 2 "TPCCluster"
        // 3 "SinglePt"
        // 4 "Chi2"
        // 5 "Qt"
        // 6 "CosPoint"
        // 7 "ToCloseV0s"
        // 8 "ConvPhi"
        // 9 "ClusterNonLinearity"
        // 10 "ClusterM02"
        // 11 "ClusterMinEnergy"
        // 12 "ClusterTrackMatching"
        // 13 "ClusterTiming"
        // 14 "ClusterNCells"
        // 15 "ClusterMaterialTRD"
        // 16 "OOBPileupGamma"
        // 17 "SPD"
        // 18 "BG"
        // 19 "Alpha"
        // 20 "Cocktail"
        // 21 "Material"
        // 22 "IntRange"
        // 23 "Efficiency"

        // Signal extraction: Pi0Yield extraction, Alpha ,BG
        errorsMeanCorrSignalExtraction[l]       =   TMath::Sqrt(pow(errorsMeanCorr[19][l],2)+   // Alpha
                                                                pow(errorsMeanCorr[18][l],2)+   // BG
                                                                pow(errorsMeanCorr[22][l],2));  // Pi0 Yield Extraction

        // Pileup:  SPD, Out-of-bunch
        errorsMeanCorrPileup[l]                 =   TMath::Sqrt(pow(errorsMeanCorr[16][l],2)+   // OOBPileupGamma
                                                                pow(errorsMeanCorr[17][l],2));  // SPD

        // PID: dEdxE, dEdxPi
        errorsMeanCorrPID[l]                    =   TMath::Sqrt(pow(errorsMeanCorr[0][l],2)+    // dEdxE
                                                                pow(errorsMeanCorr[1][l],2));   // dEdxPi

        // photon reco: Chi2+PsiPair, Qt, CosPoint, MinR
        errorsMeanCorrPhotonReco[l]             =   TMath::Sqrt(pow(errorsMeanCorr[4][l],2)+    // Chi2 PsiPair
                                                                pow(errorsMeanCorr[5][l],2)+    // Qt
                                                                pow(errorsMeanCorr[7][l],2)+    // DoubleCount
                                                                pow(errorsMeanCorr[8][l],2)+    // ConvPhi
                                                                pow(errorsMeanCorr[6][l],2));   // CosPoint

        // track reconstruction: TPCCluster, Single pT, Eta
        errorsMeanCorrTrackReco[l]              =   TMath::Sqrt(pow(errorsMeanCorr[2][l],2)+    // TPCCluster
                                                                pow(errorsMeanCorr[3][l],2));   // Single pT

        errorsMeanCorrClusterProp[l]            =   TMath::Sqrt(pow(errorsMeanCorr[11][l],2)+    // ClusterMinEnergy
                                                                pow(errorsMeanCorr[14][l],2)+    // ClusterNCells
                                                                pow(errorsMeanCorr[10][l],2)+    // ClusterM02
                                                                pow(errorsMeanCorr[13][l],2));   // ClusterTiming

    }
    TGraphErrors* meanErrorsPID                 = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrPID ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsPhotonReco          = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrPhotonReco ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsSignalExtraction    = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrSignalExtraction ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsPileup              = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrPileup ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsTrackReco           = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrTrackReco ,ptBinsErr ,errorsMeanErrCorrSummed );
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
        Double_t heightLegend2      = 5*0.04;
        if (!spectrumName.Contains("Ratio"))
            heightLegend2           = 3*0.04;
        // create legend
        TLegend* legendSummedMeanNew= GetAndSetLegend2(minXLegend2, maxYLegend2-heightLegend2, minXLegend2+widthLegend2, maxYLegend2, 40, 2, "", 43, 0.1);
        Size_t markersizeSummed     = 1.3;
        // 0 "dEdxE"
        // 1 "dEdxPi"
        // 2 "TPCCluster"
        // 3 "SinglePt"
        // 4 "Chi2"
        // 5 "Qt"
        // 6 "CosPoint"
        // 7 "ToCloseV0s"
        // 8 "ConvPhi"
        // 9 "ClusterNonLinearity"
        // 10 "ClusterM02"
        // 11 "ClusterMinEnergy"
        // 12 "ClusterTrackMatching"
        // 13 "ClusterTiming"
        // 14 "ClusterNCells"
        // 15 "ClusterMaterialTRD"
        // 16 "OOBPileupGamma"
        // 17 "SPD"
        // 18 "BG"
        // 19 "Alpha"
        // 20 "Cocktail"
        // 21 "Material"
        // 22 "IntRange"
        // 23 "Efficiency"

        // Signal extraction error
        if (benable[22] || benable[18] || benable[19]){
            DrawGammaSetMarkerTGraphErr(meanErrorsSignalExtraction, GetMarkerStyleSystematics("IntRange"), markersizeSummed, GetColorSystematics("IntRange"),GetColorSystematics("IntRange"));
            meanErrorsSignalExtraction->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"signal ext. #pi^{0}","p");
        }
        if (benable[0] || benable[1]){
            DrawGammaSetMarkerTGraphErr(meanErrorsPID, GetMarkerStyleSystematics("dEdxE"), markersizeSummed, GetColorSystematics("dEdxE"),GetColorSystematics("dEdxE"));
            meanErrorsPID->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsPID,"electron PID","p");
        }
        if (benable[2] || benable[3]){
            DrawGammaSetMarkerTGraphErr(meanErrorsTrackReco, GetMarkerStyleSystematics("SinglePt"), markersizeSummed, GetColorSystematics("SinglePt"),GetColorSystematics("SinglePt"));
            meanErrorsTrackReco->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsTrackReco,"track reco.","p");
        }
        if (benable[4] || benable[5] || benable[6] || benable[7] || benable[8]){
            DrawGammaSetMarkerTGraphErr(meanErrorsPhotonReco, GetMarkerStyleSystematics("Qt"), markersizeSummed, GetColorSystematics("Qt"),GetColorSystematics("Qt"));
            meanErrorsPhotonReco->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsPhotonReco,"photon reco.","p");
        }
        if (benable[16] || benable[17] ){
            DrawGammaSetMarkerTGraphErr(meanErrorsPileup,  GetMarkerStyleSystematics("Pileup"), markersizeSummed, GetColorSystematics("Pileup"),GetColorSystematics("Pileup"));
            meanErrorsPileup->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsPileup,"pileup","p");
        }
        if (benable[10] || benable[11] || benable[13] || benable[14] ){
            DrawGammaSetMarkerTGraphErr(meanErrorsClusterProp, GetMarkerStyleSystematics("ClusterM02"), markersizeSummed, GetColorSystematics("ClusterM02"),GetColorSystematics("ClusterM02"));
            meanErrorsClusterProp->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsClusterProp,"cluster prop.","p");
        }
        if (benable[9]){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[9], markerStyle[9], markersizeSummed,color[9],color[9]);
            meanErrorsCorr[9]->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsCorr[9],"cl. energy scale","p");
        }
        if (benable[12]){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[12], markerStyle[12], markersizeSummed,color[12],color[12]);
            meanErrorsCorr[12]->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsCorr[12],"track match. cl.","p");
        }
        if (benable[23]){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[23], markerStyle[23], markersizeSummed,color[23],color[23]);
            meanErrorsCorr[23]->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsCorr[23],"efficiency","p");
        }
        if (benable[20]){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[20], markerStyle[20], markersizeSummed,color[20],color[20]);
            meanErrorsCorr[20]->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsCorr[20],"cocktail","p");
        }
        if (benable[15]){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[15], markerStyle[15], markersizeSummed,color[15],color[15]);
            meanErrorsCorr[15]->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsCorr[15],"outer material","p");
        }
        if (benable[21]){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[21], markerStyle[21], markersizeSummed,color[21],color[21]);
            meanErrorsCorr[21]->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsCorr[21],"inner material","p");
        }


        DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummed, 20, markersizeSummed,kBlack,kBlack);
        meanErrorsCorrSummed->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsCorrSummed,"qd. sum","p");
        legendSummedMeanNew->Draw();

        labelGamma->Draw();
        if (spectrumName.Contains("Ratio")) labelPi0->Draw();
        labelEnergy->Draw();
        labelSpectrum->Draw();

    canvasSummedErrMean->Update();
    canvasSummedErrMean->SaveAs(Form("GammaSystematicErrorsCalculated/SysErrorSummedVisuPCMEMC_%s_%s_%s.%s",spectrumName.Data(), energyForOutput.Data(), dateForOutput.Data(),suffix.Data()));
    delete canvasSummedErrMean;

    //+++++++++++++++++++++++++ SAVING DETAILED SYSTEMATICS ++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatnameMeanSinglePaperErr      = Form("GammaSystematicErrorsCalculated/SystematicErrorAveragedSinglePaperPCMEMC_%s_%s_%s.dat", spectrumName.Data(),energyForOutput.Data(), dateForOutput.Data());
    fstream SysErrDatAverSinglePaper;
    cout << SysErrDatnameMeanSinglePaperErr << endl;
    SysErrDatAverSinglePaper.open(SysErrDatnameMeanSinglePaperErr, ios::out);
    SysErrDatAverSinglePaper << "Pt bin\t" ;
    if (benable[22] || benable[18] || benable[19])                          SysErrDatAverSinglePaper << "SignalExtraction" << "\t";
    if (benable[0] || benable[1])                                           SysErrDatAverSinglePaper << "ElectronPID" << "\t";
    if (benable[2] || benable[3])                                           SysErrDatAverSinglePaper << "TrackReco" << "\t";
    if (benable[4] || benable[5] || benable[6] || benable[7] || benable[8]) SysErrDatAverSinglePaper << "PhotonReco" << "\t";
    if (benable[16] || benable[17] )                                        SysErrDatAverSinglePaper << "Pileup" << "\t";
    if (benable[10] || benable[11] || benable[13] || benable[14] )          SysErrDatAverSinglePaper << "ClusterProp" << "\t";
    if (benable[9])                                                         SysErrDatAverSinglePaper << "ClusterEnergyScale" << "\t";
    if (benable[12])                                                        SysErrDatAverSinglePaper << "ClusterTM" << "\t";
    if (benable[23])                                                        SysErrDatAverSinglePaper << "Efficiency" << "\t";
    if (benable[20])                                                        SysErrDatAverSinglePaper << "Cocktail" << "\t";
    if (benable[15])                                                        SysErrDatAverSinglePaper << "OuterMaterial" << "\t";
    if (benable[21])                                                        SysErrDatAverSinglePaper << "InnerMaterial" << "\t";
    SysErrDatAverSinglePaper << endl;

    for (Int_t l=0;l< nPtBinsActive;l++){
        SysErrDatAverSinglePaper << ptBins[l] << "\t";
        if (benable[22] || benable[18] || benable[19])                          SysErrDatAverSinglePaper << errorsMeanCorrSignalExtraction[l] << "\t";
        if (benable[0] || benable[1])                                           SysErrDatAverSinglePaper << errorsMeanCorrPID[l] << "\t";
        if (benable[2] || benable[3])                                           SysErrDatAverSinglePaper << errorsMeanCorrTrackReco[l] << "\t";
        if (benable[4] || benable[5] || benable[6] || benable[7] || benable[8]) SysErrDatAverSinglePaper << errorsMeanCorrPhotonReco[l] << "\t";
        if (benable[16] || benable[17] )                                        SysErrDatAverSinglePaper << errorsMeanCorrPileup[l] << "\t";
        if (benable[10] || benable[11] || benable[13] || benable[14] )          SysErrDatAverSinglePaper << errorsMeanCorrClusterProp[l] << "\t";
        if (benable[9])                                                         SysErrDatAverSinglePaper << errorsMeanCorr[9][l] << "\t";
        if (benable[12])                                                        SysErrDatAverSinglePaper << errorsMeanCorr[12][l] << "\t";
        if (benable[23])                                                        SysErrDatAverSinglePaper << errorsMeanCorr[23][l] << "\t";
        if (benable[20])                                                        SysErrDatAverSinglePaper << errorsMeanCorr[20][l] << "\t";
        if (benable[15])                                                        SysErrDatAverSinglePaper << errorsMeanCorr[15][l] << "\t";
        if (benable[21])                                                        SysErrDatAverSinglePaper << errorsMeanCorr[21][l] << "\t";
        SysErrDatAverSinglePaper << errorsMeanCorrSummed[l] << endl;
    }
    SysErrDatAverSinglePaper.close();
}

