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
void FinaliseSystematicErrorsConv_Gammas_PbPb(   TString nameDataFileErrors      = "",
                                                TString energy                  = "",
                                                TString cent                    = "0-10%",
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
    TString centForOutput                       = cent;
    centForOutput.ReplaceAll("%","");

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
    if(!spectrumName.CompareTo("Pi0"))
        labelSpectrum                           = new TLatex(0.95,0.84,"#pi^{0}");
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
    if (energy.CompareTo("PbPb_2.76TeV") == 0 ){
        if (spectrumName.Contains("Gamma"))
            yRangesSysPlotting[1]               = 20.5;
        else if (spectrumName.Contains("IncRatio"))
            yRangesSysPlotting[1]               = 20.5;
        else
            yRangesSysPlotting[1]               = 20.5;
    }

    // Set names of cut variations for file input
    TString nameCutVariationString[14] = { "dEdxE", "dEdxPi", "TOF", "TPCCluster", "SinglePt",
                                           "Chi2" , "Qt" , "CosPoint", "Asym" , "Cocktail",
                                           "IntRange", "Alpha", "BG", "ConvPhi"};

    // Set colors and markers
    Color_t color[20];
    Color_t markerStyle[20];
    // Set names of cut variations for legends
    TString nameCutVariation[14];

    for (Int_t k =0; k<nCuts; k++ ){
        cout << "variation: " << nameCutVariationString[k].Data() << endl;
        color[k]                                = GetColorSystematics( nameCutVariationString[k] );
        markerStyle[k]                          = GetMarkerStyleSystematics( nameCutVariationString[k] );
        nameCutVariation[k]                     = GetSystematicsName(nameCutVariationString[k]);
        cout << "name for writing: " << nameCutVariation[k].Data() << endl;
    }

    // Create output folder
	TString outputdir = Form("GammaSystematicErrorsCalculated_%s",dateForOutput.Data());
	TString outputsmooth = Form("%s/Smoothing",outputdir.Data());
	gSystem->Exec("mkdir -p "+outputsmooth);

    // ***************************************************************************************************
    // ******************************** Booleans for enabling systematic errors **************************
    // ***************************************************************************************************
    Bool_t benable[14]              = { 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0,
                                        0, 0, 0, 0};
    Bool_t benableIncGamma[14]      = { 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        0, 0, 0, 1};
    Bool_t benableIncRatio[14]      = { 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 0, 0, 1};
    Bool_t benableDoubleRatio[14]   = { 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 0, 0, 1};
    Bool_t benablePi0[14]           = { 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 0, 0, 1};

    // ***************************************************************************************************
    // ******************************** Booleans for smoothing *******************************************
    // ***************************************************************************************************
    Bool_t bsmooth[14]              = { 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0,
                                        0, 0, 0, 0};
    Bool_t bsmoothIncGamma[14]      = { 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        0, 0, 0, 1};
    Bool_t bsmoothIncRatio[14]      = { 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 0, 0, 1};
    Bool_t bsmoothDoubleRatio[14]   = { 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 0, 0, 1};
    Bool_t bsmoothPi0[14]           = { 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 0, 0, 1};
    Bool_t bsmoothvoid[14]          = { 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0,
                                        0, 0, 0, 0};

    Bool_t addedSig                 = kTRUE;
    for (Int_t i = 0; i < numberCutStudies; i++){
        if(!spectrumName.CompareTo("IncRatio")){
            bsmooth[i]                      = bsmoothIncRatio[i];
            benable[i]                      = benableIncRatio[i];
        } else if(!spectrumName.CompareTo("DoubleRatio")){
            bsmooth[i]                      = bsmoothDoubleRatio[i];
            benable[i]                      = benableDoubleRatio[i];
        } else if(!spectrumName.CompareTo("Gamma")){
            bsmooth[i]                      = bsmoothIncGamma[i];
            benable[i]                      = benableIncGamma[i];
        } else if(!spectrumName.CompareTo("Pi0")){
            bsmooth[i]                      = bsmoothPi0[i];
            benable[i]                      = benablePi0[i];
        }
    }

	Int_t classTypeCut[14]     = {0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0,
                                  0, 0, 0, 0};
	Int_t classGroupCut[14]    = {0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0,
                                  0, 0, 0, 0};


    // ***************************************************************************************************
    // ****************************** Initialize error vectors & graphs **********************************
    // ***************************************************************************************************
    Double_t* errorsNeg[nCuts];
    Double_t errorsNegCorr[nCuts][nPtBins];
    Double_t errorsNegSummed[nPtBins];
    Double_t errorsNegSummedTypeA[nPtBins];
    Double_t errorsNegSummedTypeB[nPtBins];
    Double_t errorsNegSummedTypeC[nPtBins];
    Double_t errorsNegCorrSummed[nPtBins];
    Double_t errorsNegCorrSummedTypeA[nPtBins];
    Double_t errorsNegCorrSummedTypeB[nPtBins];
    Double_t errorsNegCorrSummedTypeC[nPtBins];
    Double_t errorsNegCorrMatSummed[nPtBins];

    Double_t* errorsNegErr[nCuts];
    Double_t errorsNegErrCorr[nCuts][nPtBins];
    Double_t errorsNegErrSummed[nPtBins];
    Double_t errorsNegErrCorrSummed[nPtBins];

    Double_t errorsPosErrSummedTypeA[nCuts];
    Double_t errorsMeanErrSummedTypeA[nCuts];
    Double_t errorsNegErrSummedTypeA[nCuts];
    Double_t errorsPosCorrMatSummedTypeA[nCuts];
    Double_t errorsMeanCorrMatSummedTypeA[nCuts];
    Double_t errorsNegCorrMatSummedTypeA[nCuts];
    Double_t errorsPosErrCorrSummedTypeA[nCuts];
    Double_t errorsMeanErrCorrSummedTypeA[nCuts];
    Double_t errorsMeanErrCorrMatSummedTypeA[nCuts];
    Double_t errorsNegErrCorrSummedTypeA[nCuts];

    Double_t errorsPosErrSummedTypeB[nCuts];
    Double_t errorsMeanErrSummedTypeB[nCuts];
    Double_t errorsNegErrSummedTypeB[nCuts];
    Double_t errorsPosCorrMatSummedTypeB[nCuts];
    Double_t errorsMeanCorrMatSummedTypeB[nCuts];
    Double_t errorsNegCorrMatSummedTypeB[nCuts];
    Double_t errorsPosErrCorrSummedTypeB[nCuts];
    Double_t errorsMeanErrCorrSummedTypeB[nCuts];
    Double_t errorsMeanErrCorrMatSummedTypeB[nCuts];
    Double_t errorsNegErrCorrSummedTypeB[nCuts];

    Double_t errorsPosErrSummedTypeC[nCuts];
    Double_t errorsMeanErrSummedTypeC[nCuts];
    Double_t errorsNegErrSummedTypeC[nCuts];
    Double_t errorsPosCorrMatSummedTypeC[nCuts];
    Double_t errorsMeanCorrMatSummedTypeC[nCuts];
    Double_t errorsNegCorrMatSummedTypeC[nCuts];
    Double_t errorsPosErrCorrSummedTypeC[nCuts];
    Double_t errorsMeanErrCorrSummedTypeC[nCuts];
    Double_t errorsMeanErrCorrMatSummedTypeC[nCuts];
    Double_t errorsNegErrCorrSummedTypeC[nCuts];

    Double_t* errorsPos[nCuts];
    Double_t errorsPosCorr[nCuts][nPtBins];
    Double_t errorsPosSummed[nPtBins];
    Double_t errorsPosSummedTypeA[nPtBins];
    Double_t errorsPosSummedTypeB[nPtBins];
    Double_t errorsPosSummedTypeC[nPtBins];
    Double_t errorsPosCorrSummed[nPtBins];
    Double_t errorsPosCorrSummedTypeA[nPtBins];
    Double_t errorsPosCorrSummedTypeB[nPtBins];
    Double_t errorsPosCorrSummedTypeC[nPtBins];
    Double_t errorsPosCorrMatSummed[nPtBins];

    Double_t* errorsPosErr[nCuts];
    Double_t errorsPosErrSummed[nPtBins];
    Double_t errorsPosErrCorr[nCuts][nPtBins];
    Double_t errorsPosErrCorrSummed[nPtBins];

    Double_t errorsMean[nCuts][nPtBins];
    Double_t errorsMeanCorr[nCuts][nPtBins];
    Double_t errorsMeanSummed[nPtBins];
    Double_t errorsMeanSummedTypeA[nPtBins];
    Double_t errorsMeanSummedTypeB[nPtBins];
    Double_t errorsMeanSummedTypeC[nPtBins];
    Double_t errorsMeanCorrSummed[nPtBins];
    Double_t errorsMeanCorrSummedTypeA[nPtBins];
    Double_t errorsMeanCorrSummedTypeB[nPtBins];
    Double_t errorsMeanCorrSummedTypeC[nPtBins];
    Double_t errorsMeanCorrMatSummed[nPtBins];

    Double_t errorsMeanErr[nCuts][nPtBins];
    Double_t errorsMeanErrCorr[nCuts][nPtBins];
    Double_t errorsMeanErrSummed[nPtBins];
    Double_t errorsMeanErrCorrSummed[nPtBins];
    Double_t errorsMeanErrCorrMatSummed[nPtBins];

    Double_t errorsMeanNotSmooth[nCuts][nPtBins];
    Double_t errorsMeanErrNotSmooth[nCuts][nPtBins];
    Double_t errorsMeanCorrNotSmooth[nCuts][nPtBins];
    Double_t errorsMeanErrCorrNotSmooth[nCuts][nPtBins];

    TGraphErrors* negativeErrors[nCuts];
    TGraphErrors* negativeErrorsSummed;
    TGraphErrors* negativeErrorsSummedTypeA;
    TGraphErrors* negativeErrorsSummedTypeB;
    TGraphErrors* negativeErrorsSummedTypeC;
    TGraphErrors* positiveErrors[nCuts];
    TGraphErrors* positiveErrorsSummed;
    TGraphErrors* positiveErrorsSummedTypeA;
    TGraphErrors* positiveErrorsSummedTypeB;
    TGraphErrors* positiveErrorsSummedTypeC;
    TGraphErrors* negativeErrorsCorr[nCuts];
    TGraphErrors* negativeErrorsCorrSummed;
    TGraphErrors* negativeErrorsCorrSummedTypeA;
    TGraphErrors* negativeErrorsCorrSummedTypeB;
    TGraphErrors* negativeErrorsCorrSummedTypeC;
    TGraphErrors* positiveErrorsCorr[nCuts];
    TGraphErrors* positiveErrorsCorrSummed;
    TGraphErrors* positiveErrorsCorrSummedTypeA;
    TGraphErrors* positiveErrorsCorrSummedTypeB;
    TGraphErrors* positiveErrorsCorrSummedTypeC;
    TGraphErrors* meanErrors[nCuts];
    TGraphErrors* meanErrorsSummed;
    TGraphErrors* meanErrorsSummedTypeA;
    TGraphErrors* meanErrorsSummedTypeB;
    TGraphErrors* meanErrorsSummedTypeC;
    TGraphErrors* meanErrorsCorr[nCuts];
    TGraphErrors* meanErrorsCorrNotSmooth[nCuts];
    TGraphErrors* meanErrorsCorrSummed;
    TGraphErrors* meanErrorsCorrSummedTypeA;
    TGraphErrors* meanErrorsCorrSummedTypeB;
    TGraphErrors* meanErrorsCorrSummedTypeC;
    TGraphErrors* meanErrorsCorrSummedIncMat;
    TGraphErrors* meanErrorsCorrSummedTypeAIncMat;
    TGraphErrors* meanErrorsCorrSummedTypeBIncMat;
    TGraphErrors* meanErrorsCorrSummedTypeCIncMat;

    for (Int_t l = 0; l < nPtBins; l++){
        errorsPosSummed[l]                      = 0.;
        errorsNegSummed[l]                      = 0.;
        errorsMeanSummed[l]                     = 0.;
        errorsPosCorrSummed[l]                  = 0.;
        errorsNegCorrSummed[l]                  = 0.;
        errorsMeanCorrSummed[l]                 = 0.;

        errorsPosSummedTypeA[l]                      = 0.;
        errorsNegSummedTypeA[l]                      = 0.;
        errorsMeanSummedTypeA[l]                     = 0.;
        errorsPosCorrSummedTypeA[l]                  = 0.;
        errorsNegCorrSummedTypeA[l]                  = 0.;
        errorsMeanCorrSummedTypeA[l]                 = 0.;

        errorsPosSummedTypeB[l]                      = 0.;
        errorsNegSummedTypeB[l]                      = 0.;
        errorsMeanSummedTypeB[l]                     = 0.;
        errorsPosCorrSummedTypeB[l]                  = 0.;
        errorsNegCorrSummedTypeB[l]                  = 0.;
        errorsMeanCorrSummedTypeB[l]                 = 0.;

        errorsPosSummedTypeC[l]                      = 0.;
        errorsNegSummedTypeC[l]                      = 0.;
        errorsMeanSummedTypeC[l]                     = 0.;
        errorsPosCorrSummedTypeC[l]                  = 0.;
        errorsNegCorrSummedTypeC[l]                  = 0.;
        errorsMeanCorrSummedTypeC[l]                 = 0.;

    }

    for (Int_t i = 0; i < nCuts; i++){
        if (!benable[i]) {
            cout << "*****************************************************************" << endl;
            cout << "skipping: " << nameCutVariationString[i].Data() << endl;
            cout << "*****************************************************************" << endl;
            for (Int_t l = 0; l < nPtBins; l++){
                errorsMean[i][l]                = 0;
                errorsMeanErr[i][l]             = 0.0;
                errorsMeanCorr[i][l]            = 0;
                errorsMeanErrCorr[i][l]         = 0.0;

                errorsMeanNotSmooth[i][l]                = 0;
                errorsMeanErrNotSmooth[i][l]             = 0.0;
                errorsMeanCorrNotSmooth[i][l]            = 0;
                errorsMeanErrCorrNotSmooth[i][l]         = 0.0;

            }
            continue;
        }

        TGraphAsymmErrors* graphPosErrors       = NULL;
        TGraphAsymmErrors* graphNegErrors       = NULL;

        //Set currently undetermined uncertainties
        if ( nameCutVariationString[i].CompareTo("Pileup")==0 || nameCutVariationString[i].CompareTo("SPD")==0 || nameCutVariationString[i].CompareTo("Cocktail")==0 ){
            TString nameGraphPos = Form("%s_SystErrorRelPos_%s_%s",spectrumName.Data(),nameCutVariationString[i].Data(), cent.Data());
            TString nameGraphNeg = Form("%s_SystErrorRelNeg_%s_%s",spectrumName.Data(),nameCutVariationString[i].Data(), cent.Data());
            if(cent.CompareTo("20-50%")==0){
                nameGraphPos = Form("%s_SystErrorRelPos_%s_%s",spectrumName.Data(),nameCutVariationString[i].Data(), "20-40%");
                nameGraphNeg = Form("%s_SystErrorRelNeg_%s_%s",spectrumName.Data(),nameCutVariationString[i].Data(), "20-40%");
            }
            cout << "Cutstudies " << i << "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else {
            // Load input graphs from systematics file
            TString nameGraphPos                = Form("%s_SystErrorRelPos_%s_%s",spectrumName.Data(),nameCutVariationString[i].Data(), cent.Data() );
            TString nameGraphNeg                = Form("%s_SystErrorRelNeg_%s_%s",spectrumName.Data(),nameCutVariationString[i].Data(), cent.Data() );
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                      = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        }

        // Classify error category
        // A Errors
        if( nameCutVariationString[i].CompareTo("Chi2") == 0 || nameCutVariationString[i].CompareTo("dEdxE") == 0 || nameCutVariationString[i].CompareTo("PsiPair") == 0  ||
            nameCutVariationString[i].CompareTo("TPCCluster") == 0 || nameCutVariationString[i].CompareTo("IntRange") == 0 || /*nameCutVariationString[i].CompareTo("RCut") == 0 ||*/
            nameCutVariationString[i].CompareTo("TOF") == 0 || nameCutVariationString[i].CompareTo("SinglePt") == 0 || nameCutVariationString[i].CompareTo("ConvPhi") == 0 ){
            classTypeCut[i] = 1;
        }
        // B Errors
        if( nameCutVariationString[i].CompareTo("Qt") == 0 || nameCutVariationString[i].CompareTo("dEdxPi") == 0 || /*nameCutVariationString[i].CompareTo("IntRange") == 0 ||*/
            nameCutVariationString[i].CompareTo("Asym") == 0 || nameCutVariationString[i].CompareTo("Cocktail") == 0 /*nameCutVariationString[i].CompareTo("Generator") == 0 || nameCutVariationString[i].CompareTo("SharedTrackMinDist") == 0*/ ){
            classTypeCut[i] = 2;
        }
        // C Errors
        if(nameCutVariationString[i].CompareTo("Cocktail")){
            classTypeCut[i] = 3;
        }

        // Classify error category
        // track quality
        if( nameCutVariationString[i].CompareTo("TPCCluster") 	 || nameCutVariationString[i].CompareTo("SinglePt") ){
            classGroupCut[i] = 1;
        }
        // photon reco
        if( nameCutVariationString[i].CompareTo("Qt") == 0 || nameCutVariationString[i].CompareTo("Chi2") == 0 || nameCutVariationString[i].CompareTo("PsiPair") == 0 ||
            nameCutVariationString[i].CompareTo("Asym") == 0 /*|| nameCutVariationString[i].CompareTo("RCut") == 0 || nameCutVariationString[i].CompareTo("SharedTrackMinDist") == 0 || nameCutVariationString[i].CompareTo("Generator") == 0 */){
            classGroupCut[i] = 2;
        }
        // PID
        if( nameCutVariationString[i].CompareTo("dEdxPi") == 0 || nameCutVariationString[i].CompareTo("dEdxE") == 0 || nameCutVariationString[i].CompareTo("TOF") == 0 ){
            classGroupCut[i] = 3;
        }
        // cocktail
        if(nameCutVariationString[i].CompareTo("Cocktail") == 0){
            classGroupCut[i] = 4;
        }
        //Pi0 signal extraction, alpha
        if(nameCutVariationString[i].CompareTo("IntRange") == 0 || nameCutVariationString[i].CompareTo("Alpha") == 0){
            classGroupCut[i] = 5;
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
        CalculateMeanSysErr(            errorsMean[i],  errorsMeanErr[i],   errorsPos[i],       errorsNeg[i],           nPtBinsActive);
        CorrectSystematicErrorsWithMean(errorsPos[i],   errorsPosErr[i],    errorsPosCorr[i],   errorsPosErrCorr[i],    nPtBinsActive);
        CorrectSystematicErrorsWithMean(errorsNeg[i],   errorsNegErr[i],    errorsNegCorr[i],   errorsNegErrCorr[i],    nPtBinsActive);
        CorrectSystematicErrorsWithMean(errorsMean[i],  errorsMeanErr[i],   errorsMeanCorr[i],  errorsMeanErrCorr[i],   nPtBinsActive);

        // ***************************************************************************************************
        // ************************ Adjust errors if requested to fixed values *******************************
        // ***************************************************************************************************
        Double_t errorFixed                 = -1;
        Bool_t adjustPtDependent            = kFALSE;
        if (bsmooth[i] && !addedSig){

            for (Int_t k = 0; k < nPtBinsActive; k++){
                errorsMeanNotSmooth[i][k]        = errorsMean[i][k];
                errorsMeanErrNotSmooth[i][k]     = errorsMeanErr[i][k];
                errorsMeanCorrNotSmooth[i][k]    = errorsMeanCorr[i][k];
                errorsMeanErrCorrNotSmooth[i][k] = errorsMeanErrCorr[i][k];
            }

            // fix dEdxE error #0
            if (!nameCutVariationString[i].CompareTo("dEdxE")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if (!cent.CompareTo("0-10%")){
                        adjustPtDependent               = kTRUE;
                        if (spectrumName.Contains("Ratio")){
                            errorFixed      = 0.5+pow(ptBins[k]-5,2)*0.05;
                        } else if (spectrumName.Contains("Pi0")){
                            errorFixed      = 0.3+pow(ptBins[k]-6,2)*0.06;
                        } else {
                            errorFixed      = 0.2+pow(ptBins[k],1.8)*0.02;
                        }
                    } else {
                        adjustPtDependent               = kTRUE;
                        if (spectrumName.Contains("Ratio")){
                            errorFixed      = 0.2+pow(ptBins[k],1.7)*0.05;
                        } else if (spectrumName.Contains("Pi0")){
                            errorFixed      = 0.4+pow(ptBins[k],1.8)*0.03;
                        } else {
                            errorFixed      = 0.3+pow(ptBins[k],1.8)*0.02;
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
            if (!nameCutVariationString[i].CompareTo("dEdxPi")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if (!cent.CompareTo("0-10%")){
                        adjustPtDependent               = kTRUE;
                        if (spectrumName.Contains("Ratio")){
                            errorFixed                  = 0.35+0.05*ptBins[k]+pow(ptBins[k],2.2)*0.015;
                        } else if (spectrumName.Contains("Pi0")){
                            errorFixed                  = 0.35+0.05*ptBins[k]+pow(ptBins[k],2.2)*0.015;
                        } else {
                            errorFixed                  = 0.35+0.05*ptBins[k]+pow(ptBins[k],2.2)*0.01;
                        }
                    } else {
                        adjustPtDependent               = kTRUE;
                        if (spectrumName.Contains("Ratio")){
                            errorFixed                  = 0.2+0.05*ptBins[k]+pow(ptBins[k],2.2)*0.015;
                        } else if (spectrumName.Contains("Pi0")){
                            errorFixed                  = 0.35+0.05*ptBins[k]+pow(ptBins[k],2.2)*0.015;
                        } else {
                            errorFixed                  = 0.25+0.05*ptBins[k]+pow(ptBins[k],2.2)*0.01;
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

            // fix TOF sys #2
            if (!nameCutVariationString[i].CompareTo("TOF")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    adjustPtDependent               = kTRUE;
                    if (spectrumName.Contains("Ratio")){
                        errorFixed              = 0.15+0.05*ptBins[k]+pow(ptBins[k],2.7)*0.005;
                    } else if (spectrumName.Contains("Pi0")){
                        errorFixed              = 0.15+0.05*ptBins[k]+pow(ptBins[k],2.7)*0.002;
                    } else {
                        errorFixed                  = 0.05*ptBins[k]+pow(ptBins[k],2.2)*0.008;
                    }

                    if (errorFixed != -1){
                        errorsMean[i][k]        = errorFixed;
                        errorsMeanErr[i][k]     = errorFixed*0.01;
                        errorsMeanCorr[i][k]    = errorFixed;
                        errorsMeanErrCorr[i][k] = errorFixed*0.01;
                    }
                }
            }

            // fix TPC cls sys #3
            if (!nameCutVariationString[i].CompareTo("TPCCluster")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    adjustPtDependent               = kTRUE;
                    if (!cent.CompareTo("0-10%")){
                        if (spectrumName.Contains("Ratio") || spectrumName.Contains("Pi0")){
                            errorFixed      = 0.5+pow(ptBins[k]-5,2)*0.04;
                        } else {
                            errorFixed                  = 0.3;
                        }
                    } else {
                        if (spectrumName.Contains("Ratio") || spectrumName.Contains("Pi0")){
                            errorFixed      = 0.3+pow(ptBins[k]-5,2)*0.04;
                        } else {
                            errorFixed                  = 0.3;
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

            // fix Single pt sys #4
            if (!nameCutVariationString[i].CompareTo("SinglePt")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    adjustPtDependent               = kTRUE;
                    if (spectrumName.Contains("Ratio")){
                        errorFixed              = 0.1*ptBins[k]+pow(ptBins[k]-6,2)*0.06;
                    } else if (spectrumName.Contains("Pi0")){
                        errorFixed              = 0.1*ptBins[k]+pow(ptBins[k]-6,2)*0.08;
                    } else {
                        errorFixed                  = 0.15;
                    }

                    if (errorFixed != -1){
                        errorsMean[i][k]        = errorFixed;
                        errorsMeanErr[i][k]     = errorFixed*0.01;
                        errorsMeanCorr[i][k]    = errorFixed;
                        errorsMeanErrCorr[i][k] = errorFixed*0.01;
                    }
                }
            }

            // fix Chi2/psi sys #5
            if (!nameCutVariationString[i].CompareTo("Chi2")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if(!cent.CompareTo("0-10%")){
                        adjustPtDependent           = kTRUE;
                        if (spectrumName.Contains("Ratio")){
                            if (ptBins[k]>2.5)
                                errorFixed      = 0.5+pow(ptBins[k]-5,2)*0.05;
                            else
                                errorFixed      = -0.05+0.005*ptBins[k]+pow(ptBins[k]-6,2)*0.08;
                        } else if (spectrumName.Contains("Pi0")){
                            if (ptBins[k]<3.5)
                                errorFixed      = 0.3+pow(ptBins[k]-6,2)*0.07;
                            else
                                errorFixed      = 0.7+pow(ptBins[k],1.3)*0.02;
                        } else {
                            adjustPtDependent           = kFALSE;
                            errorFixed              = 0.7;
                        }
                    } else {
                        adjustPtDependent           = kTRUE;
                        if (spectrumName.Contains("Ratio")){
                            errorFixed      = 0.5+pow(ptBins[k]-5,2)*0.05;
                        } else if (spectrumName.Contains("Pi0")){
                            errorFixed      = 1.1+pow(ptBins[k],1.8)*0.03;
                        } else {
                            errorFixed      = 0.4+pow(ptBins[k],1.8)*0.02;
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

            // fix Qt sys #6
            if (!nameCutVariationString[i].CompareTo("Qt")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    adjustPtDependent               = kTRUE;
                    if (spectrumName.Contains("Ratio")){
                        errorFixed      = 0.15+0.05*ptBins[k]+pow(ptBins[k],2)*0.03;
                    } else if (spectrumName.Contains("Pi0")){
                        errorFixed      = 0.25+0.05*ptBins[k]+pow(ptBins[k],2)*0.025;
                    } else {
                        errorFixed      = 0.15+0.05*ptBins[k]+pow(ptBins[k],2.0)*0.007;
                    }

                    if (errorFixed != -1){
                        errorsMean[i][k]        = errorFixed;
                        errorsMeanErr[i][k]     = errorFixed*0.01;
                        errorsMeanCorr[i][k]    = errorFixed;
                        errorsMeanErrCorr[i][k] = errorFixed*0.01;
                    }
                }
            }

            // fix cosPoint sys #7
            if (!nameCutVariationString[i].CompareTo("CosPoint")){
                if (!cent.CompareTo("0-10%"))
                    errorFixed                  = 0.2;
                else
                    if (spectrumName.Contains("Ratio"))
                        errorFixed                  = 0.4;
                    else
                        errorFixed                  = 0.2;
            }

            // fix asym sys #8
            if (!nameCutVariationString[i].CompareTo("Asym")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if(!cent.CompareTo("0-10%")){
                        if (spectrumName.Contains("Ratio") || spectrumName.Contains("Pi0")){
                            adjustPtDependent           = kTRUE;
                            errorFixed              = 0.1+pow(ptBins[k],2)*0.005;
                        } else {
                            adjustPtDependent           = kFALSE;
                            errorFixed              = 0.05;
                        }
                    } else {
                        if (spectrumName.Contains("Ratio") || spectrumName.Contains("Pi0")){
                            adjustPtDependent           = kTRUE;
                            errorFixed              = 0.3+pow(ptBins[k],2)*0.005;
                        } else {
                            adjustPtDependent           = kFALSE;
                            errorFixed              = 0.05;
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

            // fix Cocktail sys #9
            if (!nameCutVariationString[i].CompareTo("Cocktail")){
                errorFixed      = 0.5;
            }

            // fix IntRange sys #10
            if (!nameCutVariationString[i].CompareTo("IntRange")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    adjustPtDependent           = kTRUE;
                    if (!cent.CompareTo("0-10%")){
                        if (spectrumName.Contains("Ratio") || spectrumName.Contains("Pi0") ){
                            if (ptBins[k]<6)
                                errorFixed      = 0.2+0.07*ptBins[k]+pow(ptBins[k]-5,2)*0.1;
                            else
                                errorFixed      = 0.35+pow(ptBins[k]-6,2)*0.05;
                        }
                    } else{
                        if (spectrumName.Contains("Ratio") || spectrumName.Contains("Pi0") ){
                            errorFixed      = 0.25+0.07*ptBins[k]+pow(ptBins[k]-5,2)*0.05;
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

            // fix alpha sys #11
            if (!nameCutVariationString[i].CompareTo("Alpha")){
            }

            // fix BG sys #12
            if (!nameCutVariationString[i].CompareTo("BG")){
            }

            // fix alpha sys #13
            if (!nameCutVariationString[i].CompareTo("ConvPhi")){
                adjustPtDependent           = kTRUE;
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if (!cent.CompareTo("0-10%")){
                        if (spectrumName.Contains("Ratio") || spectrumName.Contains("Pi0") )
                            errorFixed      = 0.2+pow(ptBins[k],1.7)*0.02;
                        else {
                            errorFixed      = -1;
                        }
                    } else {
                        if (spectrumName.Contains("Ratio") || spectrumName.Contains("Pi0") )
                            errorFixed      = 0.2+pow(ptBins[k],1.7)*0.05;
                        else
                            errorFixed      = 0.1+pow(ptBins[k],2)*0.02;
                    }

                    if (errorFixed != -1){
                        errorsMean[i][k]        = errorFixed;
                        errorsMeanErr[i][k]     = errorFixed*0.01;
                        errorsMeanCorr[i][k]    = errorFixed;
                        errorsMeanErrCorr[i][k] = errorFixed*0.01;
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
        } else  if (bsmooth[i] && addedSig){

            for (Int_t k = 0; k < nPtBinsActive; k++){
                errorsMeanNotSmooth[i][k]        = errorsMean[i][k];
                errorsMeanErrNotSmooth[i][k]     = errorsMeanErr[i][k];
                errorsMeanCorrNotSmooth[i][k]    = errorsMeanCorr[i][k];
                errorsMeanErrCorrNotSmooth[i][k] = errorsMeanErrCorr[i][k];
            }

            // fix dEdxE error #0
            if (!nameCutVariationString[i].CompareTo("dEdxE")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if (!cent.CompareTo("0-10%")){
                        adjustPtDependent               = kTRUE;
                        if (spectrumName.Contains("Ratio")){
                            if (ptBins[k]<=2.5)
                                errorFixed      = 0.6+pow(ptBins[k]-3,2)*0.4;
                            else
                                errorFixed      = 0.6+pow(ptBins[k],1.7)*0.008;
                        } else if (spectrumName.Contains("Pi0")){
                            if (ptBins[k]<=2.5)
                                errorFixed      = 0.8+pow(ptBins[k]-2.5,2)*0.7;
                            else
                                errorFixed      = 0.7;
                        } else {
                            errorFixed      = 0.2;
                        }
                    } else {
                        adjustPtDependent               = kTRUE;
                        if (spectrumName.Contains("Ratio")){
                            if (ptBins[k]<=2.5)
                                errorFixed      = 0.4+pow(ptBins[k]-2.2,2)*0.5;
                            else
                                errorFixed      = 0.4+pow(ptBins[k],1.7)*0.02;
                        } else if (spectrumName.Contains("Pi0")){
                            if (ptBins[k]<=2.5)
                                errorFixed      = 0.4+pow(ptBins[k]-2.2,2)*0.7;
                            else
                                errorFixed      = 0.5;
                        } else {
                            errorFixed      = 0.18+pow(ptBins[k],1.8)*0.02;
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
            if (!nameCutVariationString[i].CompareTo("dEdxPi")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if (!cent.CompareTo("0-10%")){
                        adjustPtDependent               = kTRUE;
                        if (spectrumName.Contains("Ratio")){
                            errorFixed                  = 0.65+0.05*ptBins[k]+pow(ptBins[k],2.2)*0.009;
                        } else if (spectrumName.Contains("Pi0")){
                            errorFixed                  = 0.3+0.05*ptBins[k]+pow(ptBins[k],2.2)*0.007;
                        } else {
                            errorFixed                  = 0.4+0.05*ptBins[k]+pow(ptBins[k],2.2)*0.009;
                        }
                    } else {
                        adjustPtDependent               = kTRUE;
                        if (spectrumName.Contains("Ratio")){
                            errorFixed                  = 0.2+0.05*ptBins[k]+pow(ptBins[k],2.2)*0.01;
                        } else if (spectrumName.Contains("Pi0")){
                            errorFixed                  = 0.2+0.05*ptBins[k]+pow(ptBins[k],2.2)*0.01;
                        } else {
                            errorFixed                  = 0.17+0.05*ptBins[k]+pow(ptBins[k],2.)*0.02;
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

            // fix TOF sys #2
            if (!nameCutVariationString[i].CompareTo("TOF")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    adjustPtDependent               = kTRUE;
                    if (spectrumName.Contains("Ratio") || spectrumName.Contains("Pi0")){
                        errorFixed              = 0.15+0.05*ptBins[k]+pow(ptBins[k],2.2)*0.002;
                    } else {
                        errorFixed                  = 0.05*ptBins[k]+pow(ptBins[k],2.2)*0.001;
                    }

                    if (errorFixed != -1){
                        errorsMean[i][k]        = errorFixed;
                        errorsMeanErr[i][k]     = errorFixed*0.01;
                        errorsMeanCorr[i][k]    = errorFixed;
                        errorsMeanErrCorr[i][k] = errorFixed*0.01;
                    }
                }
            }

            // fix TPC cls sys #3
            if (!nameCutVariationString[i].CompareTo("TPCCluster")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    adjustPtDependent               = kTRUE;
                    if (!cent.CompareTo("0-10%")){
                        if (spectrumName.Contains("Ratio")){
                            if (ptBins[k]<2)
                                errorFixed      = 0.6+pow(ptBins[k]-2.,2)*0.9;
                            else
                                errorFixed      = 0.5+pow(ptBins[k]-5,2)*0.02;
                        } else if (spectrumName.Contains("Pi0")){
                            if (ptBins[k]<2)
                                errorFixed      = 0.6+pow(ptBins[k]-2.5,2)*0.5;
                            else
                                errorFixed      = 0.5+pow(ptBins[k]-5,2)*0.03;
                        } else {
                            errorFixed                  = 0.3;
                        }
                    } else {
                        if (spectrumName.Contains("Ratio") || spectrumName.Contains("Pi0")){
                            errorFixed      = 0.2+pow(ptBins[k]-5,2)*0.025;
                        } else {
                            errorFixed                  = 0.2;
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

            // fix Single pt sys #4
            if (!nameCutVariationString[i].CompareTo("SinglePt")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    adjustPtDependent               = kTRUE;
                    if (!cent.CompareTo("0-10%")){
                        if (spectrumName.Contains("Ratio")){
                            if (ptBins[k]<3)
                                errorFixed      = .8+pow(ptBins[k]-2.5,2)*0.8;
                            else
                                errorFixed      = 0.55+0.05*ptBins[k]+pow(ptBins[k],2.2)*0.001;
                        } else if (spectrumName.Contains("Pi0")){
                            if (ptBins[k]<=2.5)
                                errorFixed      = 1.+pow(ptBins[k]-2.5,2)*0.8;
                            else
                                errorFixed      = 0.7;
                        } else {
                            errorFixed          = 0.15;
                        }
                    } else {
                        if (spectrumName.Contains("Ratio")){
                            if (ptBins[k]<2.5)
                                errorFixed      = 0.7+pow(ptBins[k]-2.5,2)*0.7;
                            else
                                errorFixed      = 0.5; //0.35+0.05*ptBins[k]+pow(ptBins[k],2.2)*0.001;
                        } else if (spectrumName.Contains("Pi0")){
                            if (ptBins[k]<2.5)
                                errorFixed      = 0.8+pow(ptBins[k]-2.5,2)*0.8;
                            else
                                errorFixed      = 0.45;
                        } else {
                            errorFixed          = 0.15;
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

            // fix Chi2/psi sys #5
            if (!nameCutVariationString[i].CompareTo("Chi2")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if(!cent.CompareTo("0-10%")){
                        adjustPtDependent           = kTRUE;
                        if (spectrumName.Contains("Ratio")){
//                             if (ptBins[k]>2.5)
                                errorFixed      = 1.5+pow(ptBins[k]-4,2)*0.03;
//                             else
//                                 errorFixed      = 0.005*ptBins[k]+pow(ptBins[k]-6,2)*0.08;
                        } else if (spectrumName.Contains("Pi0")){
//                             if (ptBins[k]<3.5)
//                                 errorFixed      = 0.7+pow(ptBins[k]-6,2)*0.07;
//                             else
//                                 errorFixed      = 1.2+pow(ptBins[k],1.3)*0.02;
                                errorFixed      = 1.5+pow(ptBins[k]-5,2)*0.03;
                        } else {
                            errorFixed      = 0.4+pow(ptBins[k],1.8)*0.015;
                        }
                    } else {
                        adjustPtDependent           = kTRUE;
                        if (spectrumName.Contains("Ratio")){
                            errorFixed      = 1.+pow(ptBins[k]-3,2)*0.02;
                        } else if (spectrumName.Contains("Pi0")){
//                             if (ptBins[k]<3.5)
//                                 errorFixed      = 0.2+pow(ptBins[k]-6,2)*0.07;
//                             else
                                errorFixed      = 1.2+pow(ptBins[k],1.3)*0.02;
                        } else {
                            errorFixed      = 0.4+pow(ptBins[k],1.8)*0.015;
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

            // fix Qt sys #6
            if (!nameCutVariationString[i].CompareTo("Qt")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    adjustPtDependent               = kTRUE;
                    if (!cent.CompareTo("0-10%")){
                        if (spectrumName.Contains("Ratio")){
                            errorFixed      = 0.3+0.05*ptBins[k]+pow(ptBins[k],2)*0.01;
                        } else if (spectrumName.Contains("Pi0")){
                            errorFixed      = 0.6+0.05*ptBins[k]+pow(ptBins[k],2)*0.01;
                        } else {
                            errorFixed      = 0.15+0.05*ptBins[k]+pow(ptBins[k],2.0)*0.007;
                        }

                    } else {

                        if (spectrumName.Contains("Ratio")){
                            errorFixed      = 0.15+0.05*ptBins[k]+pow(ptBins[k],2)*0.02;
                        } else if (spectrumName.Contains("Pi0")){
                            errorFixed      = 0.5+0.05*ptBins[k]+pow(ptBins[k],2)*0.012;
                        } else {
                            errorFixed      = 0.15+0.05*ptBins[k]+pow(ptBins[k],2.0)*0.007;
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

            // fix cosPoint sys #7
            if (!nameCutVariationString[i].CompareTo("CosPoint")){
                if (!cent.CompareTo("0-10%"))
                    errorFixed                  = 0.2;
                else
                    if (spectrumName.Contains("Ratio"))
                        errorFixed                  = 0.3;
                    else
                        errorFixed                  = 0.2;
            }

            // fix asym sys #8
            if (!nameCutVariationString[i].CompareTo("Asym")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if(!cent.CompareTo("0-10%")){
                        if (spectrumName.Contains("Ratio") || spectrumName.Contains("Pi0")){
                            adjustPtDependent           = kTRUE;
                            errorFixed              = 0.2+pow(ptBins[k],2)*0.01;
                        } else {
                            adjustPtDependent           = kTRUE;
                            errorFixed              = 0.05+pow(ptBins[k],2)*0.005;
                        }
                    } else {
                        if (spectrumName.Contains("Ratio")){
                            adjustPtDependent           = kTRUE;
                            errorFixed              = 0.1+pow(ptBins[k],2)*0.007;
                        } else if (spectrumName.Contains("Pi0")){
                            adjustPtDependent           = kTRUE;
                            errorFixed              = 0.13+pow(ptBins[k],2)*0.01;
                        } else {
                            adjustPtDependent           = kTRUE;
                            errorFixed              = 0.05+pow(ptBins[k],2)*0.005;
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

            // fix Cocktail sys #9
            if (!nameCutVariationString[i].CompareTo("Cocktail")){
                if (!cent.CompareTo("0-10%")){
                    if (spectrumName.Contains("Double")){
                        adjustPtDependent           = kTRUE;
                        for (Int_t k = 0; k < nPtBinsActive; k++){
                            errorFixed      = 1.5+pow(ptBins[k],1.7)*0.02;

                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } else
                        errorFixed                  = 0.04;
                } else {
                    if (spectrumName.Contains("Double")){
                        adjustPtDependent           = kTRUE;
                        for (Int_t k = 0; k < nPtBinsActive; k++){
                            errorFixed      = 1.+ptBins[k]*0.25;

                            if (errorFixed != -1){
                                errorsMean[i][k]        = errorFixed;
                                errorsMeanErr[i][k]     = errorFixed*0.01;
                                errorsMeanCorr[i][k]    = errorFixed;
                                errorsMeanErrCorr[i][k] = errorFixed*0.01;
                            }
                        }
                    } else
                        errorFixed                  = 0.05;
                }
            }

            // fix IntRange sys #10
            if (!nameCutVariationString[i].CompareTo("IntRange")){
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    adjustPtDependent           = kTRUE;
                    if (!cent.CompareTo("0-10%")){
                        if (spectrumName.Contains("Ratio") || spectrumName.Contains("Pi0") ){
                            errorFixed      = 1.8+pow(ptBins[k],1.7)*0.01;
                        }
                    } else{
                        if (spectrumName.Contains("Ratio") || spectrumName.Contains("Pi0") ){
                            errorFixed      = 1.+pow(ptBins[k],1.7)*0.01;
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

            // fix alpha sys #11
            if (!nameCutVariationString[i].CompareTo("Alpha")){
            }

            // fix BG sys #12
            if (!nameCutVariationString[i].CompareTo("BG")){
            }

            // fix conv phi sys #13
            if (!nameCutVariationString[i].CompareTo("ConvPhi")){
                adjustPtDependent           = kTRUE;
                for (Int_t k = 0; k < nPtBinsActive; k++){
                    if (!cent.CompareTo("0-10%")){
                        if (spectrumName.Contains("Ratio") || spectrumName.Contains("Pi0") )
                            errorFixed      = 0.2+pow(ptBins[k],1.7)*0.018;
                        else {
                            errorFixed      = -1;
                        }
                    } else {
                        if (spectrumName.Contains("Ratio") )
                            errorFixed      = 0.35+pow(ptBins[k],1.7)*0.018;
                        else if(spectrumName.Contains("Pi0"))
                            errorFixed      = 0.5+pow(ptBins[k],1.7)*0.02;
                        else
                            errorFixed      = 0.12+pow(ptBins[k],1.7)*0.012;
                    }

                    if (errorFixed != -1){
                        errorsMean[i][k]        = errorFixed;
                        errorsMeanErr[i][k]     = errorFixed*0.01;
                        errorsMeanCorr[i][k]    = errorFixed;
                        errorsMeanErrCorr[i][k] = errorFixed*0.01;
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
        for (Int_t l = 0; l < nPtBinsActive; l++){
            if( nameCutVariationString[i].CompareTo("Chi2") == 0 || nameCutVariationString[i].CompareTo("dEdxE") == 0 || nameCutVariationString[i].CompareTo("PsiPair") == 0  ||
            nameCutVariationString[i].CompareTo("TPCCluster") == 0 || nameCutVariationString[i].CompareTo("IntRange") == 0 ||
            nameCutVariationString[i].CompareTo("TOF") == 0 || nameCutVariationString[i].CompareTo("SinglePt") == 0 || nameCutVariationString[i].CompareTo("ConvPhi") == 0 ){
                errorsPosSummedTypeA[l]         = errorsPosSummedTypeA[l]+pow(errorsPos[i][l],2);
                errorsNegSummedTypeA[l]         = errorsNegSummedTypeA[l]+ pow(errorsNeg[i][l],2);
                errorsMeanSummedTypeA[l]        = errorsMeanSummedTypeA[l]+ pow(errorsMean[i][l],2);
                errorsPosCorrSummedTypeA[l]     = errorsPosCorrSummedTypeA[l]+pow(errorsPosCorr[i][l],2);
                errorsNegCorrSummedTypeA[l]     = errorsNegCorrSummedTypeA[l] +pow(errorsNegCorr[i][l],2);
                errorsMeanCorrSummedTypeA[l]    = errorsMeanCorrSummedTypeA[l]+ pow(errorsMeanCorr[i][l],2);

            }
            if( nameCutVariationString[i].CompareTo("Qt") == 0 || nameCutVariationString[i].CompareTo("dEdxPi") == 0 ||
            nameCutVariationString[i].CompareTo("Asym") == 0){
                errorsPosSummedTypeB[l]         = errorsPosSummedTypeB[l]+pow(errorsPos[i][l],2);
                errorsNegSummedTypeB[l]         = errorsNegSummedTypeB[l]+ pow(errorsNeg[i][l],2);
                errorsMeanSummedTypeB[l]        = errorsMeanSummedTypeB[l]+ pow(errorsMean[i][l],2);
                errorsPosCorrSummedTypeB[l]     = errorsPosCorrSummedTypeB[l]+pow(errorsPosCorr[i][l],2);
                errorsNegCorrSummedTypeB[l]     = errorsNegCorrSummedTypeB[l] +pow(errorsNegCorr[i][l],2);
                errorsMeanCorrSummedTypeB[l]    = errorsMeanCorrSummedTypeB[l]+ pow(errorsMeanCorr[i][l],2);

            }
            if(nameCutVariationString[i].CompareTo("Cocktail") == 0){
                errorsPosSummedTypeC[l]         = errorsPosSummedTypeC[l]+pow(errorsPos[i][l],2);
                errorsNegSummedTypeC[l]         = errorsNegSummedTypeC[l]+ pow(errorsNeg[i][l],2);
                errorsMeanSummedTypeC[l]        = errorsMeanSummedTypeC[l]+ pow(errorsMean[i][l],2);
                errorsPosCorrSummedTypeC[l]     = errorsPosCorrSummedTypeC[l]+pow(errorsPosCorr[i][l],2);
                errorsNegCorrSummedTypeC[l]     = errorsNegCorrSummedTypeC[l] +pow(errorsNegCorr[i][l],2);
                errorsMeanCorrSummedTypeC[l]    = errorsMeanCorrSummedTypeC[l]+ pow(errorsMeanCorr[i][l],2);

            }
        }
        negativeErrors[i]                       = new TGraphErrors(nPtBinsActive,ptBins ,errorsNeg[i] ,ptBinsErr ,errorsNegErr[i] );
        meanErrors[i]                           = new TGraphErrors(nPtBinsActive,ptBins ,errorsMean[i] ,ptBinsErr ,errorsMeanErr[i] );
        positiveErrors[i]                       = new TGraphErrors(nPtBinsActive,ptBins ,errorsPos[i] ,ptBinsErr ,errorsPosErr[i] );
        negativeErrorsCorr[i]                   = new TGraphErrors(nPtBinsActive,ptBins ,errorsNegCorr[i] ,ptBinsErr ,errorsNegErrCorr[i] );
        meanErrorsCorr[i]                       = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorr[i] ,ptBinsErr ,errorsMeanErrCorr[i] );
        positiveErrorsCorr[i]                   = new TGraphErrors(nPtBinsActive,ptBins ,errorsPosCorr[i] ,ptBinsErr ,errorsPosErrCorr[i] );
        meanErrorsCorrNotSmooth[i]              = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrNotSmooth[i] ,ptBinsErr ,errorsMeanErrCorrNotSmooth[i] );

    }

    // Set the material budget error
    Double_t errorMaterial                      = 4.50;
    if(spectrumName.CompareTo("Pi0")==0)
        errorMaterial                           = 9.;
    else if (spectrumName.CompareTo("NoMatRatio") == 0)
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

        errorsPosSummedTypeA[l]                      = pow(errorsPosSummedTypeA[l],0.5);
        errorsMeanSummedTypeA[l]                     = pow(errorsMeanSummedTypeA[l],0.5);
        errorsPosErrSummedTypeA[l]                   = errorsPosSummedTypeA[l]*0.001;
        errorsMeanErrSummedTypeA[l]                  = errorsMeanSummedTypeA[l]*0.001;
        errorsNegSummedTypeA[l]                      = -pow(errorsNegSummedTypeA[l],0.5);
        errorsNegErrSummedTypeA[l]                   = errorsNegSummedTypeA[l]*0.001;
        errorsPosCorrMatSummedTypeA[l]               = pow(errorsPosCorrSummedTypeA[l]+ pow(errorMaterial ,2.),0.5);
        errorsMeanCorrMatSummedTypeA[l]              = pow(errorsMeanCorrSummedTypeA[l]+ pow(errorMaterial ,2.),0.5);
        errorsNegCorrMatSummedTypeA[l]               = -pow(errorsNegCorrSummedTypeA[l]+ pow(errorMaterial ,2.),0.5);

        errorsPosCorrSummedTypeA[l]                  = pow(errorsPosCorrSummedTypeA[l],0.5);
        errorsMeanCorrSummedTypeA[l]                 = pow(errorsMeanCorrSummedTypeA[l],0.5);
        errorsPosErrCorrSummedTypeA[l]               = errorsPosCorrSummedTypeA[l]*0.001;
        errorsMeanErrCorrSummedTypeA[l]              = errorsMeanCorrSummedTypeA[l]*0.001;
        errorsMeanErrCorrMatSummedTypeA[l]           = errorsMeanCorrMatSummedTypeA[l]*0.001;
        errorsNegCorrSummedTypeA[l]                  = -pow(errorsNegCorrSummedTypeA[l],0.5);
        errorsNegErrCorrSummedTypeA[l]               = errorsNegCorrSummedTypeA[l]*0.001;

        errorsPosSummedTypeB[l]                      = pow(errorsPosSummedTypeB[l],0.5);
        errorsMeanSummedTypeB[l]                     = pow(errorsMeanSummedTypeB[l],0.5);
        errorsPosErrSummedTypeB[l]                   = errorsPosSummedTypeB[l]*0.001;
        errorsMeanErrSummedTypeB[l]                  = errorsMeanSummedTypeB[l]*0.001;
        errorsNegSummedTypeB[l]                      = -pow(errorsNegSummedTypeB[l],0.5);
        errorsNegErrSummedTypeB[l]                   = errorsNegSummedTypeB[l]*0.001;
        errorsPosCorrMatSummedTypeB[l]               = pow(errorsPosCorrSummedTypeB[l]+ pow(errorMaterial ,2.),0.5);
        errorsMeanCorrMatSummedTypeB[l]              = pow(errorsMeanCorrSummedTypeB[l]+ pow(errorMaterial ,2.),0.5);
        errorsNegCorrMatSummedTypeB[l]               = -pow(errorsNegCorrSummedTypeB[l]+ pow(errorMaterial ,2.),0.5);

        errorsPosCorrSummedTypeB[l]                  = pow(errorsPosCorrSummedTypeB[l],0.5);
        errorsMeanCorrSummedTypeB[l]                 = pow(errorsMeanCorrSummedTypeB[l],0.5);
        errorsPosErrCorrSummedTypeB[l]               = errorsPosCorrSummedTypeB[l]*0.001;
        errorsMeanErrCorrSummedTypeB[l]              = errorsMeanCorrSummedTypeB[l]*0.001;
        errorsMeanErrCorrMatSummedTypeB[l]           = errorsMeanCorrMatSummedTypeB[l]*0.001;
        errorsNegCorrSummedTypeB[l]                  = -pow(errorsNegCorrSummedTypeB[l],0.5);
        errorsNegErrCorrSummedTypeB[l]               = errorsNegCorrSummedTypeB[l]*0.001;

        errorsPosSummedTypeC[l]                      = pow(errorsPosSummedTypeC[l],0.5);
        errorsMeanSummedTypeC[l]                     = pow(errorsMeanSummedTypeC[l],0.5);
        errorsPosErrSummedTypeC[l]                   = errorsPosSummedTypeC[l]*0.001;
        errorsMeanErrSummedTypeC[l]                  = errorsMeanSummedTypeC[l]*0.001;
        errorsNegSummedTypeC[l]                      = -pow(errorsNegSummedTypeC[l],0.5);
        errorsNegErrSummedTypeC[l]                   = errorsNegSummedTypeC[l]*0.001;
        errorsPosCorrMatSummedTypeC[l]               = pow(errorsPosCorrSummedTypeC[l]+ pow(errorMaterial ,2.),0.5);
        errorsMeanCorrMatSummedTypeC[l]              = pow(errorsMeanCorrSummedTypeC[l]+ pow(errorMaterial ,2.),0.5);
        errorsNegCorrMatSummedTypeC[l]               = -pow(errorsNegCorrSummedTypeC[l]+ pow(errorMaterial ,2.),0.5);

        errorsPosCorrSummedTypeC[l]                  = pow(errorsPosCorrSummedTypeC[l],0.5);
        errorsMeanCorrSummedTypeC[l]                 = pow(errorsMeanCorrSummedTypeC[l],0.5);
        errorsPosErrCorrSummedTypeC[l]               = errorsPosCorrSummedTypeC[l]*0.001;
        errorsMeanErrCorrSummedTypeC[l]              = errorsMeanCorrSummedTypeC[l]*0.001;
        errorsMeanErrCorrMatSummedTypeC[l]           = errorsMeanCorrMatSummedTypeC[l]*0.001;
        errorsNegCorrSummedTypeC[l]                  = -pow(errorsNegCorrSummedTypeC[l],0.5);
        errorsNegErrCorrSummedTypeC[l]               = errorsNegCorrSummedTypeC[l]*0.001;


            cout  << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l]
                                        << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]
                                        << "\t" << "-"<< errorsMeanCorrMatSummedTypeA[l] << "\t" <<errorsMeanCorrMatSummedTypeA[l]
                                        << "\t"  << "-"<< errorsMeanCorrMatSummedTypeB[l] << "\t" <<errorsMeanCorrMatSummedTypeB[l]
                                        << "\t" << "-"<< errorsMeanCorrMatSummedTypeC[l] << "\t" <<errorsMeanCorrMatSummedTypeC[l] << "\t"  << endl;
    }

    Double_t errorsMat[nPtBinsActive];
    for (Int_t l = 0; l < nPtBinsActive; l++){
        errorsMat[l] = errorMaterial;
    }
    TGraphErrors* graphMaterialError            = new TGraphErrors(nPtBinsActive,ptBins ,errorsMat ,ptBinsErr ,errorsMeanErrSummed );
    negativeErrorsSummed                        = new TGraphErrors(nPtBinsActive,ptBins ,errorsNegSummed ,ptBinsErr ,errorsNegErrSummed );
    negativeErrorsCorrSummed                    = new TGraphErrors(nPtBinsActive,ptBins ,errorsNegCorrSummed ,ptBinsErr ,errorsNegErrCorrSummed );
    positiveErrorsSummed                        = new TGraphErrors(nPtBinsActive,ptBins ,errorsPosSummed ,ptBinsErr ,errorsPosErrSummed );
    positiveErrorsCorrSummed                    = new TGraphErrors(nPtBinsActive,ptBins ,errorsPosCorrSummed ,ptBinsErr ,errorsPosErrCorrSummed );
    meanErrorsSummed                            = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanSummed ,ptBinsErr ,errorsMeanErrSummed );
    meanErrorsCorrSummed                        = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrSummed ,ptBinsErr ,errorsMeanErrCorrSummed );
    meanErrorsCorrSummedIncMat                  = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrMatSummed ,ptBinsErr ,errorsMeanErrCorrMatSummed );

    negativeErrorsSummedTypeA                   = new TGraphErrors(nPtBinsActive,ptBins ,errorsNegSummedTypeA ,ptBinsErr ,errorsNegErrSummedTypeA );
    negativeErrorsCorrSummedTypeA               = new TGraphErrors(nPtBinsActive,ptBins ,errorsNegCorrSummedTypeA ,ptBinsErr ,errorsNegErrCorrSummedTypeA );
    positiveErrorsSummedTypeA                   = new TGraphErrors(nPtBinsActive,ptBins ,errorsPosSummedTypeA ,ptBinsErr ,errorsPosErrSummedTypeA );
    positiveErrorsCorrSummedTypeA               = new TGraphErrors(nPtBinsActive,ptBins ,errorsPosCorrSummedTypeA ,ptBinsErr ,errorsPosErrCorrSummedTypeA );
    meanErrorsSummedTypeA                       = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanSummedTypeA ,ptBinsErr ,errorsMeanErrSummedTypeA );
    meanErrorsCorrSummedTypeA                   = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrSummedTypeA ,ptBinsErr ,errorsMeanErrCorrSummedTypeA );
    meanErrorsCorrSummedTypeAIncMat             = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrMatSummedTypeA ,ptBinsErr ,errorsMeanErrCorrMatSummedTypeA );
    meanErrorsCorrSummedTypeAIncMat->Print();

    negativeErrorsSummedTypeB                   = new TGraphErrors(nPtBinsActive,ptBins ,errorsNegSummedTypeB ,ptBinsErr ,errorsNegErrSummedTypeB );
    negativeErrorsCorrSummedTypeB               = new TGraphErrors(nPtBinsActive,ptBins ,errorsNegCorrSummedTypeB ,ptBinsErr ,errorsNegErrCorrSummedTypeB );
    positiveErrorsSummedTypeB                   = new TGraphErrors(nPtBinsActive,ptBins ,errorsPosSummedTypeB ,ptBinsErr ,errorsPosErrSummedTypeB );
    positiveErrorsCorrSummedTypeB               = new TGraphErrors(nPtBinsActive,ptBins ,errorsPosCorrSummedTypeB ,ptBinsErr ,errorsPosErrCorrSummedTypeB );
    meanErrorsSummedTypeB                       = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanSummedTypeB ,ptBinsErr ,errorsMeanErrSummedTypeB );
    meanErrorsCorrSummedTypeB                   = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrSummedTypeB ,ptBinsErr ,errorsMeanErrCorrSummedTypeB );
    meanErrorsCorrSummedTypeBIncMat             = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrMatSummedTypeB ,ptBinsErr ,errorsMeanErrCorrMatSummedTypeB );

    negativeErrorsSummedTypeC                   = new TGraphErrors(nPtBinsActive,ptBins ,errorsNegSummedTypeC ,ptBinsErr ,errorsNegErrSummedTypeC );
    negativeErrorsCorrSummedTypeC               = new TGraphErrors(nPtBinsActive,ptBins ,errorsNegCorrSummedTypeC ,ptBinsErr ,errorsNegErrCorrSummedTypeC );
    positiveErrorsSummedTypeC                   = new TGraphErrors(nPtBinsActive,ptBins ,errorsPosSummedTypeC ,ptBinsErr ,errorsPosErrSummedTypeC );
    positiveErrorsCorrSummedTypeC               = new TGraphErrors(nPtBinsActive,ptBins ,errorsPosCorrSummedTypeC ,ptBinsErr ,errorsPosErrCorrSummedTypeC );
    meanErrorsSummedTypeC                       = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanSummedTypeC ,ptBinsErr ,errorsMeanErrSummedTypeC );
    meanErrorsCorrSummedTypeC                   = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrSummedTypeC ,ptBinsErr ,errorsMeanErrCorrSummedTypeC );
    meanErrorsCorrSummedTypeCIncMat             = new TGraphErrors(nPtBinsActive,ptBins ,errorsMeanCorrMatSummedTypeC ,ptBinsErr ,errorsMeanErrCorrMatSummedTypeC );


    //++++++++++++++++++++++++++++++ PLOTTING OF SYSMEAN +++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TCanvas* canvasNewSysErrMean       = new TCanvas("canvasNewSysErrMean", "", 200, 10, 1200, 1100);  // gives the page size
    DrawGammaCanvasSettings( canvasNewSysErrMean,  0.075, 0.01, 0.015, 0.095);
    canvasNewSysErrMean->SetLogx(1);

    Double_t minXLegend     = 0.11;
    Double_t maxYLegend     = 0.95;
    Double_t widthLegend    = 0.23;
    Double_t heightLegend   = 1.12* 0.035 * (nCutsActive+1);
    Int_t nColumnsLegend    = 1;
    if ( nCutsActive+1 > 7){
        widthLegend         = 0.55;
        heightLegend        = 0.035 * (nCutsActive/1.5+1);
        nColumnsLegend      = 2;
    }

    Double_t textSizeLabelsPixelmean        = 55;
    Double_t textSizeLabelsRelmean          = 55./1200;

        TH2D *histo2DSysErrMean             = new TH2D("histo2DSysErrMean", "histo2DSysErrMean",1000, ptBins[0]-(ptBins[1]-ptBins[0])/2-0.07,1.2*(ptBins[nPtBinsActive-1]+ptBinsErr[nPtBinsActive-1]), 1000,
                                                        yRangesSysPlotting[0], yRangesSysPlotting[1]);
        SetStyleHistoTH2ForGraphs( histo2DSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean systematic Err %",
                                    0.85*textSizeLabelsRelmean, textSizeLabelsRelmean, 0.85*textSizeLabelsRelmean, textSizeLabelsRelmean, 0.9, 0.75);//(#times #epsilon_{pur})
        histo2DSysErrMean->Draw();

        TLegend* legendMean                 = GetAndSetLegend2(minXLegend, maxYLegend-heightLegend, minXLegend+widthLegend, maxYLegend, 40, nColumnsLegend, "", 43, 0.1);
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
    canvasNewSysErrMean->SaveAs(Form("%s/SysMean_%s_%s%s_%s.%s",outputdir.Data(), spectrumName.Data(), energyForOutput.Data(), centForOutput.Data(), dateForOutput.Data(),suffix.Data()));

    canvasNewSysErrMean->cd();

        if(spectrumName.CompareTo("Pi0")==0) yRangesSysPlotting[1] = yRangesSysPlotting[1]+5;
        TH2F * histo2DNewSysErrMean;
        histo2DNewSysErrMean                = new TH2F("histo2DNewSysErrMean", "histo2DNewSysErrMean",1000, ptBins[0]-(ptBins[1]-ptBins[0])/2-0.07,1.2*(ptBins[nPtBinsActive-1]+ptBinsErr[nPtBinsActive-1]), 1000.,
                                                       yRangesSysPlotting[0], yRangesSysPlotting[1]);
        SetStyleHistoTH2ForGraphs( histo2DNewSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %",
                                0.85*textSizeLabelsRelmean, textSizeLabelsRelmean, 0.85*textSizeLabelsRelmean, textSizeLabelsRelmean, 0.9, 0.75);//(#times #epsilon_{pur})
        histo2DNewSysErrMean->GetYaxis()->SetLabelOffset(0.001);
        histo2DNewSysErrMean->GetXaxis()->SetLabelOffset(-0.01);
        histo2DNewSysErrMean->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DNewSysErrMean->DrawCopy();

        TLatex *labelThesis = new TLatex(minXLegend,maxYLegend-0.03,"This thesis");
        SetStyleTLatex( labelThesis, 0.85*textSizeLabelsRelmean,4);
        labelThesis->Draw();
        // create legend
        TLegend* legendMeanNew              = GetAndSetLegend2(minXLegend, maxYLegend-heightLegend-0.04, minXLegend+widthLegend, maxYLegend-0.05, 40, nColumnsLegend, "", 43, 0.1);
        for(Int_t i = 0;i< nCuts ;i++){
            if (benable[i]){
//                 if(!spectrumName.CompareTo("Gamma")||i!=7){
                    DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], markerStyle[i], 1.,color[i],color[i]);
                    meanErrorsCorr[i]->Draw("pX0,csame");
                    legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");
//                 }
            }
        }

        DrawGammaSetMarkerTGraphErr(graphMaterialError, GetMarkerStyleSystematics( "Material" ), 1.,GetColorSystematics( "Material" ),GetColorSystematics( "Material" ));
        graphMaterialError->Draw("px0,csame");
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
    canvasNewSysErrMean->SaveAs(Form("%s/SysMeanNewWithMean_%s_%s%s_%s.%s",outputdir.Data(), spectrumName.Data(), energyForOutput.Data(), centForOutput.Data(), dateForOutput.Data(),suffix.Data()));

    canvasNewSysErrMean->cd();
    histo2DNewSysErrMean->DrawCopy();

    for(Int_t i = 0;i< nCuts ;i++){
        if (benable[i]){
            if(!spectrumName.CompareTo("Gamma")||i!=7){
                DrawGammaSetMarkerTGraphErr(meanErrorsCorrNotSmooth[i], markerStyle[i], 1.,color[i],color[i]);
                meanErrorsCorrNotSmooth[i]->Draw("pX0,csame");
            }
        }
    }

    DrawGammaSetMarkerTGraphErr(graphMaterialError, GetMarkerStyleSystematics( "Material" ), 1.,GetColorSystematics( "Material" ),GetColorSystematics( "Material" ));
    graphMaterialError->Draw("p,csame");

    DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
    meanErrorsCorrSummedIncMat->Draw("p,csame");
    legendMeanNew->Draw();

    // labeling
    labelGamma->Draw();
    labelEnergy->Draw();
    labelSpectrum->Draw();

    canvasNewSysErrMean->Update();
    canvasNewSysErrMean->SaveAs(Form("%s/SysMeanNewWithMean_NotSmooth_%s_%s%s_%s.%s",outputdir.Data(), spectrumName.Data(), energyForOutput.Data(), centForOutput.Data(), dateForOutput.Data(),suffix.Data()));

    for(Int_t i = 0;i< nCuts ;i++){
        if (benable[i]){

            canvasNewSysErrMean->cd();
            TH2D *histo2DCheckSmooth = new TH2D("histo2DCheckSmooth", "histo2DCheckSmooth", 20,0.,ptBins[nPtBinsActive-1]+1,1000.,-0.5,20.);
                SetStyleHistoTH2ForGraphs( histo2DCheckSmooth, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %",
                                    0.85*55./1200, 55./1200, 0.85*55./1200, 55./1200, 0.9, 0.75);
                histo2DCheckSmooth->GetYaxis()->SetLabelOffset(0.001);
                histo2DCheckSmooth->GetXaxis()->SetLabelOffset(-0.01);
        //      histo2DCheckSmooth->GetXaxis()->SetMoreLogLabels(kTRUE);
                histo2DCheckSmooth->DrawCopy();

                DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], markerStyle[i], 1.,color[i],color[i]);
                meanErrorsCorr[i]->Draw("pX0,csame");

                DrawGammaSetMarkerTGraphErr(meanErrorsCorrNotSmooth[i], markerStyle[i], 1.,kGray+2,kGray+2);
                meanErrorsCorrNotSmooth[i]->Draw("pX0,csame");

            canvasNewSysErrMean->SaveAs(Form("%s/SysBeforeSmoothing_%s_%s_%s%s_%s.%s",outputsmooth.Data(), spectrumName.Data(), nameCutVariationString[i].Data(), energyForOutput.Data(), centForOutput.Data(), dateForOutput.Data(),suffix.Data()));
        }
    }


    //+++++++++++++++++++++++++ SAVING SYSTEMATICS TO DAT FILE +++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatname                   = Form("%s/SystematicError_%s_%s%s_%s.dat",outputdir.Data(), spectrumName.Data(), energyForOutput.Data(), centForOutput.Data(), dateForOutput.Data());
    fstream SysErrDat;
    cout << SysErrDatname << endl;
    SysErrDat.open(SysErrDatname, ios::out);
    for (Int_t l=0; l< nPtBinsActive; l++){
        SysErrDat << ptBins[l] << "\t" <<errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
    }
    SysErrDat.close();

    //+++++++++++++++++++++++++ SAVING AVERAGE SYSTEMATICS TO DAT ++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatnameMean               = Form("%s/SystematicErrorAveraged_%s_%s%s_%s.dat", outputdir.Data(), spectrumName.Data(), energyForOutput.Data(), centForOutput.Data(),  dateForOutput.Data());
    fstream SysErrDatAver;
    cout << SysErrDatnameMean << endl;
    SysErrDatAver.open(SysErrDatnameMean, ios::out);
    for (Int_t l=0; l< nPtBinsActive; l++){
        SysErrDatAver  << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
        cout  << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
    }
    SysErrDatAver.close();

    const char *SysErrDatnameMeanSepType               = Form("%s/SystematicErrorAveragedSepErrType_%s_%s%s_%s.dat", outputdir.Data(), spectrumName.Data(), energyForOutput.Data(), centForOutput.Data(),  dateForOutput.Data());
    fstream SysErrDatAverSepType;
    cout << SysErrDatnameMeanSepType << endl;
    SysErrDatAverSepType.open(SysErrDatnameMeanSepType, ios::out);
    for (Int_t l=0; l< nPtBinsActive; l++){
        if(!spectrumName.CompareTo("Pi0")){
            SysErrDatAverSepType  << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
            cout  << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l]
                                        << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
        } else {
            SysErrDatAverSepType  << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l]
                                        << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]
                                        << "\t" << "-"<< errorsMeanCorrMatSummedTypeA[l] << "\t" <<errorsMeanCorrMatSummedTypeA[l]
                                        << "\t"  << "-"<< errorsMeanCorrSummedTypeB[l] << "\t" <<errorsMeanCorrSummedTypeB[l]
                                        << "\t" << "-"<< errorsMeanCorrMatSummedTypeC[l] << "\t" <<errorsMeanCorrMatSummedTypeC[l] << "\t"  << endl;

            cout  << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l]
                                        << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]
                                        << "\t" << "-"<< errorsMeanCorrMatSummedTypeA[l] << "\t" <<errorsMeanCorrMatSummedTypeA[l]
                                        << "\t"  << "-"<< errorsMeanCorrMatSummedTypeB[l] << "\t" <<errorsMeanCorrMatSummedTypeB[l]
                                        << "\t" << "-"<< errorsMeanCorrMatSummedTypeC[l] << "\t" <<errorsMeanCorrMatSummedTypeC[l] << "\t"  << endl;
        }
    }
    SysErrDatAverSepType.close();

    //+++++++++++++++++++++++++ SAVING DETAILED SYSTEMATICS ++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatnameMeanSingleErr      = Form("%s/SystematicErrorAveragedSinglePCM_%s_%s%s_%s.dat",outputdir.Data(),  spectrumName.Data(),energyForOutput.Data(), centForOutput.Data(), dateForOutput.Data());
    fstream SysErrDatAverSingle;
    cout << SysErrDatnameMeanSingleErr << endl;
    SysErrDatAverSingle.open(SysErrDatnameMeanSingleErr, ios::out);
    SysErrDatAverSingle << "Pt bin\t" ;
    for (Int_t i= 0; i< numberCutStudies; i++){
        if (benable[i]) SysErrDatAverSingle << nameCutVariationString[i] << "\t";
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
        // 0 - "dEdxE", 1 - "dEdxPi", 2 - "TOF", 3 - "TPCCluster", 4 - "SinglePt",
        // 5 - "Chi2" , 6 - "Qt" , 7 - "CosPoint", 8 - "Asym" , 9 - "Cocktail"
        // 10 - "IntRange", 11 - "Alpha", 12 - "BG"

        // PID: dEdxE, dEdxPi, TOF
        errorsMeanCorrPID[l]                    =   TMath::Sqrt(pow(errorsMeanCorr[0][l],2)+    // dEdxE
                                                                pow(errorsMeanCorr[1][l],2)+  // dEdxPi
                                                                pow(errorsMeanCorr[2][l],2));  // TOF
        // track reconstruction: TPCCluster, Single pT, Eta
        errorsMeanCorrTrackReco[l]              =   TMath::Sqrt(pow(errorsMeanCorr[3][l],2)+    // TPCCluster
                                                                pow(errorsMeanCorr[4][l],2));   // Single pT

        // photon reco: Chi2+PsiPair, Qt, CosPoint, Asym
        errorsMeanCorrPhotonReco[l]             =   TMath::Sqrt(pow(errorsMeanCorr[5][l],2)+    // Chi2 PsiPair
                                                                pow(errorsMeanCorr[6][l],2)+    // Qt
                                                                pow(errorsMeanCorr[7][l],2)+    // CosPoint
                                                                pow(errorsMeanCorr[8][l],2));  // Asym

        // Signal extraction: Pi0Yield extraction, Alpha (,BG)
        errorsMeanCorrSignalExtraction[l]       =   TMath::Sqrt(//pow(errorsMeanCorr[11][l],2)+   // Alpha
                                                                //pow(errorsMeanCorr[12][l],2)+    // BG
                                                                pow(errorsMeanCorr[10][l],2));  // Pi0 Yield Extraction
        // Pileup:  SPD, Out-of-bunch, Out-of-bunch pi0?
//         errorsMeanCorrPileup[l]                 =   TMath::Sqrt(pow(errorsMeanCorr[13][l],2)+    // SPD
//                                                                 pow(errorsMeanCorr[14][l],2));   // Pileup

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
//     canvasSummedErrMean->SetLogx(1);

        TH2F * histo2DSummedErrMean;
        histo2DSummedErrMean            = new TH2F("histo2DSummedErrMean", "histo2DSummedErrMean",1000, ptBins[0]-(ptBins[1]-ptBins[0])/2-0.07, 1.2*(ptBins[nPtBinsActive-1]+ptBinsErr[nPtBinsActive-1]),
                                                   1000.,yRangesSysPlotting[0],yRangesSysPlotting[1] );
        SetStyleHistoTH2ForGraphs(  histo2DSummedErrMean, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %",
                                    0.85*textSizeLabelsRel, textSizeLabelsRel, 0.85*textSizeLabelsRel, textSizeLabelsRel, 0.9, 0.75);//(#times #epsilon_{pur})
        histo2DSummedErrMean->GetYaxis()->SetLabelOffset(0.001);
        histo2DSummedErrMean->GetXaxis()->SetLabelOffset(-0.01);
        histo2DSummedErrMean->GetXaxis()->SetMoreLogLabels(kTRUE);
        histo2DSummedErrMean->DrawCopy();

        Double_t minXLegend2        = 0.11;
        Double_t maxYLegend2        = 0.95;
        Double_t widthLegend2       = 0.5;
        Double_t heightLegend2      = 0.15;

        // create legend
        TLegend* legendSummedMeanNew    = GetAndSetLegend2(minXLegend2, maxYLegend2-heightLegend2, minXLegend2+widthLegend2, maxYLegend2, 40, 2, "", 43, 0.1);
        Size_t markersizeSummed     = 1.3;

        // PID
        if (benable[0] || benable[1]|| benable[2]){
            DrawGammaSetMarkerTGraphErr(meanErrorsPID, 21, markersizeSummed,color[1],color[1]);
            meanErrorsPID->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsPID,"Electron PID","p");
        }
        // Track reconstruction
        if (benable[3] || benable[4]){
            DrawGammaSetMarkerTGraphErr(meanErrorsTrackReco, 22, markersizeSummed,color[2],color[2]);
            meanErrorsTrackReco->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsTrackReco,"Track Reco.","p");
        }
        // Photon reconstruction
        if (benable[5] || benable[6] || benable[7] || benable[8] || benable[14]){
            DrawGammaSetMarkerTGraphErr(meanErrorsPhotonReco, 23, markersizeSummed,color[3],color[3]);
            meanErrorsPhotonReco->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsPhotonReco,"Photon Reco.","p");
        }
        // Signal extraction error
        if (benable[11] || benable[12] || benable[10]){
            DrawGammaSetMarkerTGraphErr(meanErrorsSignalExtraction, 20, markersizeSummed,color[0],color[0]);
            meanErrorsSignalExtraction->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"Signal Extraction #pi^{0}","p");
        }
        // Pile up
        if (benable[13]){
            DrawGammaSetMarkerTGraphErr(meanErrorsPileup, 25, markersizeSummed,color[5],color[5]);
            meanErrorsPileup->Draw("p,csame");
            legendSummedMeanNew->AddEntry(meanErrorsPileup,"Pileup","p");
        }
        // Cocktail
        if (benable[9]){
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
    canvasSummedErrMean->SaveAs(Form("%s/SysErrorSummedVisu_%s_%s%s_%s.%s",outputdir.Data(), spectrumName.Data(), energyForOutput.Data(), centForOutput.Data(), dateForOutput.Data(),suffix.Data()));

    delete canvasSummedErrMean;

    //+++++++++++++++++++++++++ SAVING SYSTEMATICS PAPER STYLE +++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const char *SysErrDatnameMeanPaper          = Form("%s/SystematicErrorAveragedPaper_%s_%s%s_%s.dat", outputdir.Data(), spectrumName.Data(),energyForOutput.Data(), centForOutput.Data(), dateForOutput.Data());
    fstream SysErrDatAverPaper;
    cout << SysErrDatnameMeanPaper << endl;
    SysErrDatAverPaper.open(SysErrDatnameMeanPaper, ios::out);
    SysErrDatAverPaper  << "p_{T}" << "\t Material \t Yield Extraction \t PID \t photon reco \t track recon \t summed" <<  endl;
    for (Int_t l=0; l< nPtBinsActive; l++){
        SysErrDatAverPaper << ptBins[l] <<"\t" << errorMaterial << "\t" << errorsMeanCorrSignalExtraction[l] << "\t" << errorsMeanCorrPID[l]<< "\t" << errorsMeanCorrPhotonReco[l]<< "\t" <<errorsMeanCorrTrackReco[l] <<"\t" << errorsMeanCorrMatSummed[l]<< endl;
    }
    SysErrDatAverPaper.close();
}
