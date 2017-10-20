// ------------------------------------------------
// This code finalizes the systematic errors for --
// the direct photon analysis. Gamma spectrum as --
// well as double ratio systematics are produced --
// ------------------------------------------------
// Code is maintained by:                        --
//          - Nicolas Schmidt                    --
//          - Lucas Altenkämper                  --
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
void FinaliseSystematicErrorsConv_Gammas_ppV2(  TString nameDataFileErrors      = "",
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
    TString detectionProcess                    = ReturnFullTextReconstructionProcess(mode);
    TString dateForOutput                       = ReturnDateStringForOutput();
    TString energyForOutput                     = ReturnCollisionEnergyOutputString(energy);

    TLatex *labelGamma                          = new TLatex(0.95,0.88,detectionProcess);
    SetStyleTLatex( labelGamma, 0.038,4);
    labelGamma->SetTextAlign(31);

    TLatex *labelEnergy                         = new TLatex(0.95,0.92,collisionSystem);
    SetStyleTLatex( labelEnergy, 0.038,4);
    labelEnergy->SetTextAlign(31);
    TLatex *labelSpectrum;
    if(!spectrumName.CompareTo("IncRatio"))
        labelSpectrum                           = new TLatex(0.95,0.84,"#gamma_{inc}/#pi^{0}");
    if(!spectrumName.CompareTo("DoubleRatio"))
        labelSpectrum                           = new TLatex(0.95,0.84,"R_{#gamma}");
    if(!spectrumName.CompareTo("Gamma"))
        labelSpectrum                           = new TLatex(0.95,0.84,"#gamma_{inc}");
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
            yRangesSysPlotting[1]               = 15.5;
        else
            yRangesSysPlotting[1]               = 28.5;
    }

    // Set names of cut variations for file input
    TString nameCutVariationSC[16]              = {"dEdxE", "TPCCluster", "SinglePt", "Chi2" , "Qt" ,
                                                    "DoubleCount" ,"BG" ,"Periods" ,"SPD" ,"OOBPileupGamma",
                                                    "CosPoint","dEdxPi", "Alpha" ,"Cocktail","RCut",
                                                    "IntRange" };

    // Set colors and markers
    Color_t color[20];
    Color_t markerStyle[20];
    // Set names of cut variations for legends
    TString nameCutVariation[16];

    if(!energy.CompareTo("2.76TeV"))
        nameCutVariationSC[7]                   = "OOBPileupPi0";
    if(!energy.CompareTo("8TeV")||!energy.CompareTo("7TeV"))
        nameCutVariationSC[14]                   = "OOBPileupPi0";
    if((!energy.CompareTo("8TeV")||!energy.CompareTo("7TeV"))&&!spectrumName.CompareTo("Gamma"))
        nameCutVariationSC[10]                   = "Cocktail";

    for (Int_t k =0; k<nCuts; k++ ){
        cout << "variation: " << nameCutVariationSC[k].Data() << endl;
        color[k]                                = GetColorSystematics( nameCutVariationSC[k] );
        markerStyle[k]                          = GetMarkerStyleSystematics( nameCutVariationSC[k] );
        nameCutVariation[k]                     = GetSystematicsName(nameCutVariationSC[k]);
        cout << "name for writing: " << nameCutVariation[k].Data() << "\t"<< color[k]  << "\t" << markerStyle[k] << endl;
    }


    // adapt some cut name for different energies
    if(!energy.CompareTo("7TeV"))
        nameCutVariationSC[7]                   = "7TeVPeriods";
    if(!energy.CompareTo("8TeV"))
        nameCutVariationSC[7]                   = "8TeVPeriods";
    if(!energy.CompareTo("900GeV"))
        nameCutVariationSC[7]                   = "BG";


    // Create output folder
    gSystem->Exec("mkdir -p GammaSystematicErrorsCalculated");

    // ***************************************************************************************************
    // ******************************** Booleans for enabling systematic errors **************************
    // ***************************************************************************************************
    Bool_t benable[16]                          = { 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0,  0 };
    Bool_t benableIncGamma900GeV[16]            = { 1, 1, 1, 1, 1,  1, 1, 0, 1, 1,
                                                    1, 1, 1, 0, 0,  0 };
    Bool_t benableIncRatio900GeV[16]            = { 1, 1, 1, 1, 1,  1, 1, 0, 1, 1,
                                                    1, 1, 1, 0, 0,  0 };
    Bool_t benableDR900GeV[16]                  = { 1, 1, 1, 1, 1,  1, 1, 0, 1, 1,
                                                    1, 1, 1, 1, 0,  0 };
    Bool_t benableIncGamma2760GeV[16]           = { 1, 1, 1, 1, 1,  0, 0, 0, 1, 1,
                                                    0, 1, 0, 0, 0,  0 };
    Bool_t benableIncRatio2760GeV[16]           = { 1, 1, 1, 1, 1,  0, 1, 1, 0, 1,
                                                    0, 1, 1, 0, 0,  1 };
    Bool_t benableDR2760GeV[16]                 = { 1, 1, 1, 1, 1,  0, 1, 1, 0, 1,
                                                    0, 1, 1, 1, 0,  1 };
    Bool_t benableIncGamma7TeV[16]              = { 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0,  0 };
    Bool_t benableIncRatio7TeV[16]              = { 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0,  0 };
    Bool_t benableDR7TeV[16]                    = { 1, 1, 1, 1, 1,  1, 1, 1, 1, 1,
                                                    1, 1, 1, 1, 0,  0 };
    Bool_t benableIncGamma8TeV[16]              = { 1, 1, 1, 1, 1,  1, 1, 1, 1, 1,
                                                    1, 1, 1, 0, 0,  0 };
    Bool_t benableIncRatio8TeV[16]              = { 1, 1, 1, 1, 1,  1, 1, 1, 1, 1,
                                                    1, 1, 1, 1, 1,  0 };
    Bool_t benableDR8TeV[16]                    = { 1, 1, 1, 1, 1,  1, 1, 1, 1, 1,
                                                    1, 1, 1, 1, 1,  1 };

    // ***************************************************************************************************
    // ******************************** Booleans for smoothing *******************************************
    // ***************************************************************************************************
    Bool_t bsmooth[16]                          = { 0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothIncGamma900GeV[16]            = { 0, 0, 0, 0, 1,  0, 0, 0, 1, 1,
                                                    0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothIncRatio900GeV[16]            = { 0, 0, 0, 0, 1,  0, 0, 0, 1, 1,
                                                    0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothDR900GeV[16]                  = { 0, 0, 0, 0, 1,  0, 0, 0, 1, 1,
                                                    0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothIncGamma2760GeV[16]           = { 1, 1, 1, 1, 1,  0, 0, 0, 1, 1,
                                                    0, 1, 0, 0, 0,  0 };
    Bool_t bsmoothIncRatio2760GeV[16]           = { 1, 1, 1, 1, 1,  0, 1, 0, 1, 1,
                                                    0, 1, 1, 0, 0,  0 };
    Bool_t bsmoothDR2760GeV[16]                 = { 1, 1, 1, 1, 1,  0, 1, 0, 1, 1,
                                                    0, 1, 1, 0, 0,  0 };
    Bool_t bsmoothIncGamma7TeV[16]              = { 0, 0, 0, 0, 1,  0, 0, 1, 1, 1,
                                                    0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothIncRatio7TeV[16]              = { 0, 0, 0, 0, 1,  0, 0, 1, 1, 1,
                                                    0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothDR7TeV[16]                    = { 0, 0, 1, 1, 1,  0, 0, 1, 1, 1,
                                                    0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothIncGamma8TeV[16]              = { 1, 0, 0, 1, 0,  0, 0, 1, 0, 1,
                                                    0, 0, 0, 0, 0,  0 };
    Bool_t bsmoothIncRatio8TeV[16]              = { 0, 0, 1, 1, 1,  1, 0, 1, 0, 0,
                                                    0, 0, 1, 0, 0,  0 };
    Bool_t bsmoothDR8TeV[16]                    = { 1, 0, 1, 1, 1,  1, 0, 1, 0, 1,
                                                    0, 0, 1, 0, 0,  1 };

    for (Int_t i = 0; i < numberCutStudies; i++){
        if (energy.CompareTo("900GeV") == 0){
            if(!spectrumName.CompareTo("IncRatio")){
                bsmooth[i]                      = bsmoothIncRatio900GeV[i];
                benable[i]                      = benableIncRatio900GeV[i];
            } else if(!spectrumName.CompareTo("DoubleRatio")){
                bsmooth[i]                      = bsmoothIncRatio900GeV[i];
                benable[i]                      = benableIncRatio900GeV[i];
            } else if(!spectrumName.CompareTo("Gamma")){
                bsmooth[i]                      = bsmoothIncGamma900GeV[i];
                benable[i]                      = benableIncGamma900GeV[i];
            }
        } else if (energy.CompareTo("2.76TeV") == 0){
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
        } else if (energy.CompareTo("7TeV") == 0){
            if(!spectrumName.CompareTo("IncRatio")){
                bsmooth[i]                      = bsmoothIncRatio7TeV[i];
                benable[i]                      = benableIncRatio7TeV[i];
            } else if(!spectrumName.CompareTo("DoubleRatio")){
                bsmooth[i]                      = bsmoothDR7TeV[i];
                benable[i]                      = benableDR7TeV[i];
            } else if(!spectrumName.CompareTo("Gamma")){
                bsmooth[i]                      = bsmoothIncGamma7TeV[i];
                benable[i]                      = benableIncGamma7TeV[i];
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
        if ( (nameCutVariationSC[i].CompareTo("SPD")==0 && !energy.CompareTo("2.76TeV")) ){
            TString nameGraphPos;
            TString nameGraphNeg;
            nameGraphPos                        = Form("%s_SystErrorRelPos_%s_pp",spectrumName.Data(),nameCutVariationSC[0].Data()  );
            nameGraphNeg                        = Form("%s_SystErrorRelNeg_%s_pp",spectrumName.Data(),nameCutVariationSC[0].Data()  );
            cout << "Cutstudies " << i << "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else if ( nameCutVariationSC[i].CompareTo("OOBPileupPi0")==0 ){
            TString nameGraphPos                = Form("Pi0_SystErrorRelPos_%s_pp",nameCutVariationSC[i].Data());
            TString nameGraphNeg                = Form("Pi0_SystErrorRelNeg_%s_pp",nameCutVariationSC[i].Data());
            cout << "Cutstudies " << i << "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
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
            TString nameGraphPos                = Form("%s_SystErrorRelPos_%s_pp",spectrumName.Data(),nameCutVariationSC[i].Data() );
            TString nameGraphNeg                = Form("%s_SystErrorRelNeg_%s_pp",spectrumName.Data(),nameCutVariationSC[i].Data() );
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        }
        if ( graphPosErrors == NULL ){
            cout << "systematic wasn't contained, setting it to 0" << endl;
            TString nameGraphPos            = Form("%s_SystErrorRelPos_%s_pp",spectrumName.Data(),nameCutVariationSC[0].Data() );
            TString nameGraphNeg            = Form("%s_SystErrorRelNeg_%s_pp",spectrumName.Data(),nameCutVariationSC[0].Data() );
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
            Double_t errorFixed                 = -1;
            Bool_t adjustPtDependent            = kFALSE;

            // fix dEdx e error #0
            if (!nameCutVariationSC[i].CompareTo("dEdxE")){
                if (spectrumName.Contains("Ratio")){
                    adjustPtDependent               = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                        if (!energy.CompareTo("2.76TeV")){
                            errorFixed              = 1.05+pow(ptBins[k],2.7)*0.045;
                        } else if (!energy.CompareTo("8TeV")){
                          if (ptBins[k] < 2.0)
                            errorFixed              = 0.7+pow(ptBins[k]-2.0,2)*0.75;
                          else 
                            errorFixed              = 0.7;
                        } else {
                            errorFixed              = 0.25+pow(ptBins[k]-2.0,2)*0.7;
                        }

                        if (errorFixed != -1){
                            errorsMean[i][k]        = errorFixed;
                            errorsMeanErr[i][k]     = errorFixed*0.01;
                            errorsMeanCorr[i][k]    = errorFixed;
                            errorsMeanErrCorr[i][k] = errorFixed*0.01;
                        }
                    }
                } else {
                    if (!energy.CompareTo("2.76TeV")){
                        errorFixed                  = 0.45;
                    } else if (!energy.CompareTo("8TeV")){
                        adjustPtDependent               = kTRUE;
                        for (Int_t k = 0; k < nPtBinsActive; k++){
                            if (ptBins[k] > 1.0)
                                errorFixed              = 0.5;
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    }
                }
            }

            // fix Single pt sys #1
            if (!nameCutVariationSC[i].CompareTo("TPCCluster")){
                if (spectrumName.Contains("Ratio")){
                    if (!energy.CompareTo("2.76TeV")){
                        adjustPtDependent               = kTRUE;
                        for (Int_t k = 0; k < nPtBinsActive; k++){
                            errorFixed          = 0.4+pow(ptBins[k],2.5)*0.045+0.2*ptBins[k];

                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } else {
                        errorFixed                  = 0.3;
                    }
                } else {
                    if (!energy.CompareTo("2.76TeV")){
                        errorFixed                  = 0.2;
                    }
                }
            }

            // fix Single pt sys #2
            if (!nameCutVariationSC[i].CompareTo("SinglePt")){
                if (spectrumName.Contains("Ratio")){
                    if (!energy.CompareTo("7TeV")||!energy.CompareTo("8TeV")){
                      errorFixed                    = 0.5;
                    } else if (!energy.CompareTo("2.76TeV")){
                        adjustPtDependent               = kTRUE;
                        for (Int_t k = 0; k < nPtBinsActive; k++){
                            if (ptBins[k] > 0.8)
                                errorFixed              = 0.5+pow(ptBins[k]+2.8,1.8)*0.05;

                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } else {
                        errorFixed                  = 0.3;
                    }
                } else {
                    if (!energy.CompareTo("2.76TeV")||!energy.CompareTo("7TeV")||!energy.CompareTo("8TeV")){
                        errorFixed                  = 0.25;
                    }
                }
            }

            // fix Chi2/psi pair uncertainties #3
            if (!nameCutVariationSC[i].CompareTo("Chi2")){
                if (spectrumName.Contains("Ratio")){
                    adjustPtDependent           = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                        if (!energy.CompareTo("8TeV")){
                          if (ptBins[k]<4)
                            errorFixed          = 1;
                            else
                             errorFixed          = 0.8+pow(ptBins[k]-4,2)*0.015+0.2*(ptBins[k]-3);
                        } else if (!energy.CompareTo("7TeV")){
                            if (ptBins[k]>4)
                                errorFixed          = 0.1+pow(ptBins[k]-4,2)*0.03+0.5*(ptBins[k]-3);
                        } else if (!energy.CompareTo("2.76TeV")){
                            if (ptBins[k]>0.8)
                                errorFixed          = 1.+pow(ptBins[k],2.2)*0.03+0.5*ptBins[k];
                        } else {
                            errorFixed              = 0.2+pow(ptBins[k],2)*0.025;
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
                        if (!energy.CompareTo("2.76TeV")||!energy.CompareTo("7TeV"))
                            errorFixed              = 0.2+pow(ptBins[k],1.9)*0.02+0.1*ptBins[k];
                        else if (!energy.CompareTo("8TeV"))
                            if (ptBins[k] > 1)
                                errorFixed              = 1.2;

                        if (errorFixed != -1){
                            errorsMean[i][k]        = errorFixed;
                            errorsMeanErr[i][k]     = errorFixed*0.01;
                            errorsMeanCorr[i][k]    = errorFixed;
                            errorsMeanErrCorr[i][k] = errorFixed*0.01;
                        }
                    }
                }
            }

            // fix Qt sys #4
            if (!nameCutVariationSC[i].CompareTo("Qt")){
                adjustPtDependent               = kTRUE;
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if (!energy.CompareTo("2.76TeV")){
                        if (spectrumName.Contains("Ratio")){
                            if (ptBins[k] > 0.6)
                                errorFixed      = 0.5+pow(ptBins[k],2.5)*0.04+0.5*ptBins[k];
                        } else {
                            errorFixed      = 0.15+pow(ptBins[k],2.0)*0.003+0.15*ptBins[k];
                        }
                    } else if (!energy.CompareTo("8TeV") && spectrumName.Contains("Double")){
                        if (ptBins[k] > 3.0)
                            errorFixed          = 0.1+pow(ptBins[k]+3,2)*0.02;
                    } else {
                        if (ptBins[k] > 3.0)
                            errorFixed          = 0.2+pow(ptBins[k]+3,2)*0.027;
                    }

                    if (errorFixed != -1){
                        errorsMean[i][k]        = errorFixed;
                        errorsMeanErr[i][k]     = errorFixed*0.01;
                        errorsMeanCorr[i][k]    = errorFixed;
                        errorsMeanErrCorr[i][k] = errorFixed*0.01;
                    }
                }
            }

            // fix double counting sys #5
            if (!nameCutVariationSC[i].CompareTo("DoubleCount")){
                adjustPtDependent               = kTRUE;
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if (!energy.CompareTo("8TeV")){
                        errorFixed              = 0.01+pow(ptBins[k]-0.5,2)*0.027;
                    }
                    if (errorFixed != -1){
                        errorsMean[i][k]        = errorFixed;
                        errorsMeanErr[i][k]     = errorFixed*0.01;
                        errorsMeanCorr[i][k]    = errorFixed;
                        errorsMeanErrCorr[i][k] = errorFixed*0.01;
                    }
                }
            }

            // fix BG sys #6
            if (!nameCutVariationSC[i].CompareTo("BG")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if(spectrumName.Contains("Ratio")){
                        errorFixed              = 2.2;
                    }else{
                        errorFixed              = 0;
                    }
                }
            }

            // fix periods sys #7
            if (!nameCutVariationSC[i].CompareTo("7TeVPeriods")||!nameCutVariationSC[i].CompareTo("8TeVPeriods")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if(spectrumName.CompareTo("Gamma")){
                        errorFixed              = 0;
                    }else{
                        errorFixed              = 1.5;
                    }
                }
            }

            // fix SPD pileup sys #8
            if (!nameCutVariationSC[i].CompareTo("SPD")){
                if (!energy.CompareTo("2.76TeV")){
                    errorFixed                  = 0.25;
                } else {
                    adjustPtDependent               = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                        if (ptBins[k] > 12.0)
                            errorFixed              = 1.0;
                        if (errorFixed != -1){
                            errorsMean[i][k]        = errorFixed;
                            errorsMeanErr[i][k]     = errorFixed*0.01;
                            errorsMeanCorr[i][k]    = errorFixed;
                            errorsMeanErrCorr[i][k] = errorFixed*0.01;
                        }
                    }
                }
            }
            // fix out-of-bunch sys #9
            if ( !nameCutVariationSC[i].CompareTo("OOBPileupGamma") ){
                if (!energy.CompareTo("2.76TeV")){
                    adjustPtDependent               = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                        if (ptBins[k] > 5)
                            errorFixed              = 1.5;
                        if (errorFixed != -1){
                            errorsMean[i][k]        = errorFixed;
                            errorsMeanErr[i][k]     = errorFixed*0.01;
                            errorsMeanCorr[i][k]    = errorFixed;
                            errorsMeanErrCorr[i][k] = errorFixed*0.01;
                        }
                    }
                } else if (!energy.CompareTo("8TeV")){
                    adjustPtDependent               = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                        if (ptBins[k] < 3)
                            errorFixed              = 2.2+pow(ptBins[k]-3,2)*0.25;
                        if (ptBins[k] > 3)
                            errorFixed              = 2.2+pow(ptBins[k]-3,2)*0.037;
                        if (errorFixed != -1){
                            errorsMean[i][k]        = errorFixed;
                            errorsMeanErr[i][k]     = errorFixed*0.01;
                            errorsMeanCorr[i][k]    = errorFixed;
                            errorsMeanErrCorr[i][k] = errorFixed*0.01;
                        }
                    }
                } else {
                    errorFixed                      = 1.0;
                }
            }

            // fix cosPoint sys #10

            // fix dEdxPi sys #11
            if (!nameCutVariationSC[i].CompareTo("dEdxPi")){
                if (spectrumName.Contains("Ratio")){
                    if (!energy.CompareTo("2.76TeV")){
                        adjustPtDependent               = kTRUE;
                        for (Int_t k = 0; k < nPtBinsActive; k++){
                                errorFixed              = 0.9+pow(ptBins[k],2.7)*0.04+ptBins[k]*0.1;

                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    }
                } else {
                    if (!energy.CompareTo("2.76TeV")){
                        adjustPtDependent               = kTRUE;
                        for (Int_t k = 0; k < nPtBinsActive; k++){
                            errorFixed                  = 0.2+pow(ptBins[k],2.5)*0.005+0.05*ptBins[k];

                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    }
                }
            }

            // fix alpha uncertainties #12
            if (!nameCutVariationSC[i].CompareTo("Alpha")){
                if (!spectrumName.CompareTo("IncRatio")){
                    if (!energy.CompareTo("2.76TeV"))
                        errorFixed              = 1.;
                    else if (!energy.CompareTo("8TeV")){
                        adjustPtDependent           = kTRUE;
                        for (Int_t k = 0; k < nPtBinsActive; k++){
                            errorFixed              = 0.01+pow(ptBins[k],2)*0.029;
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    }
                }
                if (!spectrumName.CompareTo("DoubleRatio")){
                    if (!energy.CompareTo("2.76TeV")){
                        errorFixed              = 1;
                    } else {
                        adjustPtDependent           = kTRUE;
                        for (Int_t k = 0; k < nPtBinsActive; k++){
                            errorFixed              = 0.01+pow(ptBins[k],2)*0.029;
                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    }
                }
            }

            // fix Cocktail sys #13
            if (!nameCutVariationSC[i].CompareTo("Cocktail")){
                if (!spectrumName.CompareTo("DoubleRatio")){
                    errorFixed                  = 1.0;
                }
            }

            // fix Rconv uncertainties #14
            if (!nameCutVariationSC[i].CompareTo("Rcut")){
                if (!spectrumName.CompareTo("DoubleRatio")){
                    adjustPtDependent           = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                        errorFixed              = 1.6+pow(ptBins[k]-0.5,2)*0.029;
                        if (errorFixed != -1){
                            errorsMean[i][k]        = errorFixed;
                            errorsMeanErr[i][k]     = errorFixed*0.01;
                            errorsMeanCorr[i][k]    = errorFixed;
                            errorsMeanErrCorr[i][k] = errorFixed*0.01;
                        }
                    }
                }
            }

            // fix IntRange sys #15
            if (!nameCutVariationSC[i].CompareTo("IntRange")){
                if (spectrumName.Contains("Ratio")){
                    adjustPtDependent           = kTRUE;
                    for (Int_t k = 0; k < nPtBinsActive; k++){
                        if (!energy.CompareTo("2.76TeV")){
                            if (ptBins[k] > 8.0)
                                errorFixed          = 3;
                        }
                        if (!energy.CompareTo("8TeV")){
                          if (ptBins[k] < 0.6)
                              errorFixed              = 4+pow(ptBins[k]-0.6,2)*90;
                          else if (ptBins[k] < 3 && ptBins[k] > 0.6)
                              errorFixed              =1.2+pow(ptBins[k]-3,2)*0.25;
                          else if (ptBins[k] > 3)
                              errorFixed              = 1.2+pow(ptBins[k]-3,2)*0.08;
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

    // Set the material budget error
    Double_t errorMaterial                      = 4.50;
    if (spectrumName.CompareTo("NoMatRatio") == 0)
        errorMaterial                           = 0.;

    for (Int_t l = 0; l < nPtBinsActive; l++){
        errorsPosSummed[l]                      = pow(errorsPosSummed[l],0.5);
        errorsMeanSummed[l]                     = pow(errorsMeanSummed[l],0.5);
        errorsPosErrSummed[l]                   = errorsPosSummed[l]*0.001;
        errorsMeanErrSummed[l]                  = errorsMeanSummed[l]*0.001;
        errorsNegSummed[l]                      = -pow(errorsNegSummed[l],0.5);
        errorsNegErrSummed[l]                   = errorsNegSummed[l]*0.001;
        errorsPosCorrMatSummed[l]               = pow(errorsPosCorrSummed[l]+ pow(errorMaterial ,2.),0.5);
        errorsMeanCorrMatSummed[l]              = pow(errorsMeanCorrSummed[l]+ pow(errorMaterial ,2.),0.5);
        errorsNegCorrMatSummed[l]               = -pow(errorsNegCorrSummed[l]+ pow(errorMaterial ,2.),0.5);

        errorsPosCorrSummed[l]                  = pow(errorsPosCorrSummed[l],0.5);
        errorsMeanCorrSummed[l]                 = pow(errorsMeanCorrSummed[l],0.5);
        errorsPosErrCorrSummed[l]               = errorsPosCorrSummed[l]*0.001;
        errorsMeanErrCorrSummed[l]              = errorsMeanCorrSummed[l]*0.001;
        errorsMeanErrCorrMatSummed[l]           = errorsMeanCorrMatSummed[l]*0.001;
        errorsNegCorrSummed[l]                  = -pow(errorsNegCorrSummed[l],0.5);
        errorsNegErrCorrSummed[l]               = errorsNegCorrSummed[l]*0.001;
    }

    Double_t errorsMat[nPtBinsActive];
    for (Int_t l = 0; l < nPtBinsActive; l++){
        errorsMat[l]                            = errorMaterial;
    }
    TGraphErrors* graphMaterialError            = new TGraphErrors(nPtBinsActive,ptBins ,errorsMat ,ptBinsErr ,errorsMeanErrSummed );

    negativeErrorsSummed                        = new TGraphErrors(nPtBinsActive,ptBins ,errorsNegSummed ,ptBinsErr ,errorsNegErrSummed );
    negativeErrorsCorrSummed                    = new TGraphErrors(nPtBinsActive,ptBins ,errorsNegCorrSummed ,ptBinsErr ,errorsNegErrCorrSummed );
    positiveErrorsSummed                        = new TGraphErrors(nPtBinsActive,ptBins ,errorsPosSummed ,ptBinsErr ,errorsPosErrSummed );
    positiveErrorsCorrSummed                    = new TGraphErrors(nPtBinsActive,ptBins ,errorsPosCorrSummed ,ptBinsErr ,errorsPosErrCorrSummed );
    meanErrorsSummed                            = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanSummed ,ptBinsErr ,errorsMeanErrSummed );
    meanErrorsCorrSummed                        = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrSummed ,ptBinsErr ,errorsMeanErrCorrSummed );
    meanErrorsCorrSummedIncMat                  = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrMatSummed ,ptBinsErr ,errorsMeanErrCorrMatSummed );

    //++++++++++++++++++++++++++++++ PLOTTING OF SYSMEAN +++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TCanvas* canvasNewSysErrMean       = new TCanvas("canvasNewSysErrMean", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasNewSysErrMean,  0.075, 0.01, 0.015, 0.095);
    canvasNewSysErrMean->SetLogx(1);

    Double_t minXLegend     = 0.12;
    Double_t maxYLegend     = 0.95;

    Double_t widthLegend    = 0.25;
    if ( nCutsActive+1 > 7)
        widthLegend         = 0.58;
    Double_t heightLegend   = 1.12* 0.035 * (nCutsActive+1);
    if ( nCutsActive+1 > 7)
        heightLegend        = 1.05* 0.035 * (nCutsActive/1.5+1);

    Double_t textSizeLabelsPixelmean        = 55;
    Double_t textSizeLabelsRelmean          = 55./1200;

        TH2D *histo2DSysErrMean                     = new TH2D("histo2DSysErrMean", "histo2DSysErrMean",1000, ptBins[0]-(ptBins[1]-ptBins[0])/2-0.07,1.2*(ptBins[nPtBinsActive-1]+ptBinsErr[nPtBinsActive-1]), 1000,
                                                               yRangesSysPlotting[0], yRangesSysPlotting[1]);
        SetStyleHistoTH2ForGraphs( histo2DSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean systematic Err %",
                                    0.85*textSizeLabelsRelmean, textSizeLabelsRelmean, 0.85*textSizeLabelsRelmean, textSizeLabelsRelmean, 0.9, 0.75);//(#times #epsilon_{pur})
        histo2DSysErrMean->Draw();

        TLegend* legendMean = GetAndSetLegend2(minXLegend,maxYLegend-heightLegend,minXLegend+widthLegend,maxYLegend, 40);
        legendMean->SetMargin(0.1);
        if ( nCutsActive+1 > 7) legendMean->SetNColumns(2);
        for(Int_t i = 0; i< numberCutStudies ; i++){
            if (benable[i]){
                DrawGammaSetMarkerTGraphErr(meanErrors[i], markerStyle[i], 1.,color[i],color[i]);
                meanErrors[i]->Draw("p,csame");
                legendMean->AddEntry(meanErrors[i],nameCutVariation[i].Data(),"p");
            }
        }

        DrawGammaSetMarkerTGraphErr(meanErrorsSummed, 20, 1.,1,1);
        meanErrorsSummed->Draw("p,csame");
        legendMean->AddEntry(meanErrorsSummed,"qd. sum w/o mat.","p");
        legendMean->Draw();
        labelEnergy->Draw();
        labelGamma->Draw();
        labelSpectrum->Draw();
    canvasNewSysErrMean->Update();
    canvasNewSysErrMean->SaveAs(Form("GammaSystematicErrorsCalculated/SysMean_%s_%s_%s.%s",spectrumName.Data(), energyForOutput.Data(),dateForOutput.Data(),suffix.Data()));

    canvasNewSysErrMean->cd();

        TH2F * histo2DNewSysErrMean;
        histo2DNewSysErrMean                = new TH2F("histo2DNewSysErrMean", "histo2DNewSysErrMean",1000, ptBins[0]-(ptBins[1]-ptBins[0])/2-0.07,1.2*(ptBins[nPtBinsActive-1]+ptBinsErr[nPtBinsActive-1]), 1000.,
                                                       yRangesSysPlotting[0], yRangesSysPlotting[1]);
        SetStyleHistoTH2ForGraphs( histo2DNewSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %",
                                0.85*textSizeLabelsRelmean, textSizeLabelsRelmean, 0.85*textSizeLabelsRelmean, textSizeLabelsRelmean, 0.9, 0.75);//(#times #epsilon_{pur})
        histo2DNewSysErrMean->GetYaxis()->SetLabelOffset(0.001);
        histo2DNewSysErrMean->GetXaxis()->SetLabelOffset(-0.01);
        histo2DNewSysErrMean->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DNewSysErrMean->DrawCopy();

        // create legend
        TLegend* legendMeanNew = GetAndSetLegend2(minXLegend,maxYLegend-heightLegend,minXLegend+widthLegend,maxYLegend, 40);
        legendMeanNew->SetMargin(0.1);
        if ( nCutsActive+1 > 7) legendMeanNew->SetNColumns(2);

        for(Int_t i = 0;i< nCuts ;i++){
            if (benable[i]){
                DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], markerStyle[i], 1.,color[i],color[i]);
                meanErrorsCorr[i]->Draw("pX0,csame");
                legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");

            }
        }

        DrawGammaSetMarkerTGraphErr(graphMaterialError, GetMarkerStyleSystematics( "Material" ), 1.,GetColorSystematics( "Material" ),GetColorSystematics( "Material" ));
        graphMaterialError->Draw("p,csame");
        legendMeanNew->AddEntry(graphMaterialError,"material","p");

        DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
        meanErrorsCorrSummedIncMat->Draw("p,csame");
        legendMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"qd. sum","p");
        legendMeanNew->Draw();

        // labeling
        labelGamma->Draw();
        labelEnergy->Draw();
        labelSpectrum->Draw();

    canvasNewSysErrMean->Update();
    canvasNewSysErrMean->SaveAs(Form("GammaSystematicErrorsCalculated/SysMeanNewWithMean_%s_%s_%s.%s",spectrumName.Data(), energyForOutput.Data(),dateForOutput.Data(),suffix.Data()));

    //+++++++++++++++++++++++++ SAVING SYSTEMATICS TO DAT FILE +++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatname                   = Form("GammaSystematicErrorsCalculated/SystematicError_%s_%s_%s.dat",spectrumName.Data(),energyForOutput.Data(),dateForOutput.Data());
    fstream SysErrDat;
    cout << SysErrDatname << endl;
    SysErrDat.open(SysErrDatname, ios::out);
    for (Int_t l=0; l< nPtBinsActive; l++){
        SysErrDat << ptBins[l] << "\t" <<errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
    }
    SysErrDat.close();

    //+++++++++++++++++++++++++ SAVING AVERAGE SYSTEMATICS TO DAT ++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatnameMean               = Form("GammaSystematicErrorsCalculated/SystematicErrorAveraged_%s_%s_%s.dat",spectrumName.Data(),energyForOutput.Data(),dateForOutput.Data());
    fstream SysErrDatAver;
    cout << SysErrDatnameMean << endl;
    SysErrDatAver.open(SysErrDatnameMean, ios::out);
    for (Int_t l=0; l< nPtBinsActive; l++){
        SysErrDatAver  << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
        // SysErrDatAver << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
    }
    SysErrDatAver.close();

    //+++++++++++++++++++++++++ SAVING DETAILED SYSTEMATICS ++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatnameMeanSingleErr      = Form("GammaSystematicErrorsCalculated/SystematicErrorAveragedSinglePCM_%s_%s_%s.dat",spectrumName.Data(),energyForOutput.Data(),dateForOutput.Data());
    fstream SysErrDatAverSingle;
    cout << SysErrDatnameMeanSingleErr << endl;
    SysErrDatAverSingle.open(SysErrDatnameMeanSingleErr, ios::out);
    SysErrDatAverSingle << "Pt bin\t" ;
    for (Int_t i= 0; i< numberCutStudies; i++){
        if (benable[i]) SysErrDatAverSingle << nameCutVariationSC[i] << "\t";
    }
    SysErrDatAverSingle << "Material";
    SysErrDatAverSingle << endl;
    for (Int_t l=0;l< nPtBinsActive;l++){
        SysErrDatAverSingle << ptBins[l] << "\t";
        for (Int_t i= 0; i< numberCutStudies; i++){
            if (benable[i]) SysErrDatAverSingle << errorsMeanCorr[i][l] << "\t";
        }
        SysErrDatAverSingle << 4.5 << "\t";
        SysErrDatAverSingle << errorsMeanCorrMatSummed[l] << endl;
    }
    SysErrDatAverSingle.close();

    Double_t errorsMeanCorrPID[nPtBinsActive];
    Double_t errorsMeanCorrSignalExtraction[nPtBinsActive];
    Double_t errorsMeanCorrPileup[nPtBinsActive];
    Double_t errorsMeanCorrTrackReco[nPtBinsActive];
    Double_t errorsMeanCorrPhotonReco[nPtBinsActive];

    for (Int_t l=0; l< nPtBinsActive; l++){
        //"0 dEdxE"
        //"1 TPCCluster"
        //"2 SinglePt"
        //"3 Chi2"
        //"4 Qt"
        //"5 DoubleCount"
        //"6 BG"
        //"7 Periods" // "OOBPileupPi0"
        //"8 SPD"
        //"9 OOBPileupGamma"
        //"10 CosPoint"
        //"11 dEdxPi"
        //"12 Alpha"
        //"13 Cocktail"
        //"14 Rcut"
        //"15 IntRange"

        // Signal extraction: Pi0Yield extraction, Alpha ,BG
        errorsMeanCorrSignalExtraction[l]       =   TMath::Sqrt(pow(errorsMeanCorr[12][l],2)+   // Alpha
                                                                pow(errorsMeanCorr[6][l],2)+    // BG
                                                                pow(errorsMeanCorr[15][l],2));  // Pi0 Yield Extraction

        // Pileup:  SPD, Out-of-bunch, Out-of-bunch pi0?
        errorsMeanCorrPileup[l]                 =   TMath::Sqrt(pow(errorsMeanCorr[7][l],2)+    // OOBPileupPi0
                                                                pow(errorsMeanCorr[8][l],2)+    // SPD
                                                                pow(errorsMeanCorr[9][l],2));   // OOBPileupGamma

        // PID: dEdxE, dEdxPi
        errorsMeanCorrPID[l]                    =   TMath::Sqrt(pow(errorsMeanCorr[0][l],2)+    // dEdxE
                                                                pow(errorsMeanCorr[11][l],2));  // dEdxPi

        // photon reco: Chi2+PsiPair, Qt, CosPoint, MinR
        errorsMeanCorrPhotonReco[l]             =   TMath::Sqrt(pow(errorsMeanCorr[3][l],2)+    // Chi2 PsiPair
                                                                pow(errorsMeanCorr[4][l],2)+    // Qt
                                                                pow(errorsMeanCorr[5][l],2)+    // DoubleCount
                                                                pow(errorsMeanCorr[10][l],2));  // CosPoint

        // track reconstruction: TPCCluster, Single pT, Eta
        errorsMeanCorrTrackReco[l]              =   TMath::Sqrt(pow(errorsMeanCorr[1][l],2)+    // TPCCluster
                                                                pow(errorsMeanCorr[2][l],2));   // Single pT
    }
    TGraphErrors* meanErrorsPID                 = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrPID ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsPhotonReco          = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrPhotonReco ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsSignalExtraction    = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrSignalExtraction ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsPileup              = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrPileup ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsTrackReco           = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrTrackReco ,ptBinsErr ,errorsMeanErrCorrSummed );

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
        histo2DSummedErrMean->GetXaxis()->SetLabelOffset(-0.01);
        histo2DSummedErrMean->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DSummedErrMean->DrawCopy();

        Double_t minXLegend2        = 0.13;
        Double_t maxYLegend2        = 0.95;
        Double_t widthLegend2       = 0.55;
        Double_t heightLegend2      = 4*0.04;
        if (!spectrumName.Contains("Ratio"))
            heightLegend2           = 3*0.04;


        // create legend
        TLegend* legendSummedMeanNew    = GetAndSetLegend2(minXLegend2, maxYLegend2-heightLegend2, minXLegend2+widthLegend2, maxYLegend2, 40, 2, "", 43, 0.1);
        Size_t markersizeSummed     = 1.3;

        // Signal extraction error
        if (benable[12] || benable[6] || benable[15]){
            DrawGammaSetMarkerTGraphErr(meanErrorsSignalExtraction, GetMarkerStyleSystematics("IntRange"), markersizeSummed, GetColorSystematics("IntRange"),GetColorSystematics("IntRange"));
            meanErrorsSignalExtraction->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"Signal Ext. #pi^{0}","p");
        }
        if (benable[0] || benable[11]){
            DrawGammaSetMarkerTGraphErr(meanErrorsPID, GetMarkerStyleSystematics("dEdxE"), markersizeSummed, GetColorSystematics("dEdxE"),GetColorSystematics("dEdxE"));
            meanErrorsPID->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsPID,"Electron PID","p");
        }
        if (benable[1] || benable[2]){
            DrawGammaSetMarkerTGraphErr(meanErrorsTrackReco, GetMarkerStyleSystematics("SinglePt"), markersizeSummed, GetColorSystematics("SinglePt"),GetColorSystematics("SinglePt"));
            meanErrorsTrackReco->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsTrackReco,"Track Reco.","p");
        }
        if (benable[3] || benable[4] || benable[5] || benable[10]){
            DrawGammaSetMarkerTGraphErr(meanErrorsPhotonReco, GetMarkerStyleSystematics("Qt"), markersizeSummed, GetColorSystematics("Qt"),GetColorSystematics("Qt"));
            meanErrorsPhotonReco->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsPhotonReco,"Photon Reco.","p");
        }
        if (benable[9] || benable[8] || benable[7]){
            DrawGammaSetMarkerTGraphErr(meanErrorsPileup,  GetMarkerStyleSystematics("Pileup"), markersizeSummed, GetColorSystematics("Pileup"),GetColorSystematics("Pileup"));
            meanErrorsPileup->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsPileup,"Pile-up","p");
        }
        if (benable[13]){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[13], markerStyle[13], markersizeSummed,color[13],color[13]);
            meanErrorsCorr[13]->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsCorr[13],"Cocktail","p");
        }
        DrawGammaSetMarkerTGraphErr(graphMaterialError,  GetMarkerStyleSystematics( "Material" ), 1.,GetColorSystematics( "Material" ),GetColorSystematics( "Material" ));
        graphMaterialError->Draw("p,csame");
        legendSummedMeanNew->AddEntry(graphMaterialError,"Material","p");

        DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, markersizeSummed,kBlack,kBlack);
        meanErrorsCorrSummedIncMat->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"qd. sum","p");
        legendSummedMeanNew->Draw();

        labelGamma->Draw();
        labelEnergy->Draw();
        labelSpectrum->Draw();

    canvasSummedErrMean->Update();
    canvasSummedErrMean->SaveAs(Form("GammaSystematicErrorsCalculated/SysErrorSummedVisu_%s_%s_%s.%s",spectrumName.Data(), energyForOutput.Data(),dateForOutput.Data(),suffix.Data()));

    delete canvasSummedErrMean;

    //+++++++++++++++++++++++++ SAVING SYSTEMATICS PAPER STYLE +++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatnameMeanPaper          = Form("GammaSystematicErrorsCalculated/SystematicErrorAveragedPaper_%s_%s_%s.dat",spectrumName.Data(),energyForOutput.Data(),dateForOutput.Data());
    fstream SysErrDatAverPaper;
    cout << SysErrDatnameMeanPaper << endl;
    SysErrDatAverPaper.open(SysErrDatnameMeanPaper, ios::out);
    SysErrDatAverPaper  << "p_{T}" << "\t Material \t Yield Extraction \t PID \t photon reco \t track recon \t summed" <<  endl;
    for (Int_t l=0; l< nPtBinsActive; l++){
        SysErrDatAverPaper << ptBins[l] <<"\t" << errorMaterial << "\t" << errorsMeanCorrSignalExtraction[l] << "\t" << errorsMeanCorrPID[l]<< "\t" << errorsMeanCorrPhotonReco[l]<< "\t" <<errorsMeanCorrTrackReco[l] <<"\t" << errorsMeanCorrMatSummed[l]<< endl;
    }
    SysErrDatAverPaper.close();
}
