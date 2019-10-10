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

void FinaliseSystematicErrorsMergedCalo_v2(
    TString nameDataFileErrors      = "",
    TString energy                  = "",
    TString meson                   = "",
    Int_t numberOfPtBins            = 1,
    Int_t numberCutStudies          = 1,
    Double_t startPtSys             = 0,
    TString additionalNameOutput    = "",
    TString suffix                  = "eps"
){

    if (numberCutStudies > 13) {
        cout << "ERROR: Too many variations set" << endl;
        return ;
    }

    // ***************************************************************************************************
    // ****************************** General style settings *********************************************
    // ***************************************************************************************************
    StyleSettingsThesis();
    SetPlotStyle();

    // ***************************************************************************************************
    // ****************************** Create output directory ********************************************
    // ***************************************************************************************************
    gSystem->Exec("mkdir -p SystematicErrorsCalculatedMergedCalo");

    // ***************************************************************************************************
    // ***************************** labeling and color settings *****************************************
    // ***************************************************************************************************
    TString date                            = ReturnDateString();
    TString dateForOutput                   = ReturnDateStringForOutput();
    TString collisionSystem                 = ReturnFullCollisionsSystem(energy);
    TString energyForOutput                 = energy;
    energyForOutput.ReplaceAll(".","_");

    TString textMeson                       = ReturnMesonString (meson);
    if(meson.Contains("Ratio"))
        textMeson = "#it{R}_{pA}^{#pi^{0}}";
    TString labelMesonReco                  = textMeson;
    TString labelTriggerText                = "";
    if (additionalNameOutput.Contains("INT7") ){
        labelTriggerText = "INT7";
    }else if (additionalNameOutput.Contains("EMC7") ){
        labelTriggerText = "EMC7";
    }else if (additionalNameOutput.Contains("EGA") ){
        labelTriggerText = "EGA";
    }else if (additionalNameOutput.Contains("EG2") ){
        labelTriggerText = "EG2";
    }else if (additionalNameOutput.Contains("EG1") ){
        labelTriggerText = "EG1";
    }
    TString labelDataSetText                = "";
    if(!energy.CompareTo("8TeV")){
        labelDataSetText = "LHC12[a-i]";
    }else if(!energy.CompareTo("pPb_8TeV")){
        labelDataSetText = "LHC16[r-s]";
    }
    // ***************************************************************************************************
    // ******************************* general variable definition  **************************************
    // ***************************************************************************************************
    Int_t   numberOfEntriesPos              = 0;
    Int_t   numberOfEntriesNeg              = 0;
    const Int_t nPtBins                     = numberOfPtBins;
    const Int_t nCuts                       = numberCutStudies;
    Double_t* ptBins;
    Double_t* ptBinsErr;
    // Double_t* ptBinsErr;
    TString nameCutVariation[13];
    TString nameCutVariationSC[13];
    TString nameCutVariationSCAll[13] = {
        "ClusterM02", "ClusterTrackMatchingCalo", "ClusterNonLinearity",  "CellMinE", "CellTiming",
        "ClusterMaterialTRD", "MesonResolution", "ClusterEnergyScale" , "Trigger", "Efficiency",
        "Secondary", "Purity", "DistanceBadChannel"
    };
    Color_t color[20];
    Color_t markerStyle[20];
    for (Int_t k = 0; k < 13; k++ ){
        color[k]        = GetColorSystematics( nameCutVariationSCAll[k], 10 );
        markerStyle[k]  = GetMarkerStyleSystematics( nameCutVariationSCAll[k], 10 );

        nameCutVariation[k]     = GetSystematicsName(nameCutVariationSCAll[k]);
        nameCutVariationSC[k]   = nameCutVariationSCAll[k];
    }

    // ***************************************************************************************************
    // ******************************** Booleans for smoothing *******************************************
    // ***************************************************************************************************
    Bool_t bsmooth[13]                      = { 0, 0, 0, 0, 0,      1, 0, 0, 0, 0,      0, 0, 0  };
    // pp 8TeV
    Bool_t bsmoothpp8TeVINT7Pi0[13]         = { 1, 1, 1, 1, 1,      1, 1, 1, 1, 1,      1, 1, 1  };
    Bool_t bsmoothpp8TeVEMC7Pi0[13]         = { 1, 1, 1, 1, 1,      1, 1, 1, 1, 1,      1, 1, 1  };
    Bool_t bsmoothpp8TeVEGAPi0[13]          = { 1, 1, 1, 1, 1,      1, 1, 1, 1, 1,      1, 1, 1  };
    // pPb 8.16TeV
    Bool_t bsmoothpPb8TeVINT7Pi0[13]        = { 1, 0, 1, 1, 1,      1, 1, 1, 1, 1,      1, 1, 1  };
    Bool_t bsmoothpPb8TeVEG2Pi0[13]         = { 1, 0, 1, 1, 1,      1, 1, 1, 1, 1,      1, 1, 1  };
    Bool_t bsmoothpPb8TeVEG1Pi0[13]         = { 1, 0, 1, 1, 1,      1, 1, 1, 1, 1,      1, 1, 1  };
    Bool_t bsmoothpPb8TeVINT7Pi0Ratio[13]   = { 1, 1, 1, 0, 1,      1, 1, 1, 1, 1,      1, 1, 1  };
    Bool_t bsmoothpPb8TeVEG2Pi0Ratio[13]    = { 1, 1, 1, 0, 1,      1, 1, 1, 1, 1,      1, 1, 1  };
    Bool_t bsmoothpPb8TeVEG1Pi0Ratio[13]    = { 1, 1, 1, 0, 1,      1, 1, 1, 1, 1,      1, 1, 1  };

    for (Int_t i = 0; i < numberCutStudies; i++){
        if(!energy.CompareTo("8TeV")){
            if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                bsmooth[i]      = bsmoothpp8TeVINT7Pi0[i];
            }else if (additionalNameOutput.Contains("EMC7") && meson.CompareTo("Pi0")==0){
                bsmooth[i]      = bsmoothpp8TeVEMC7Pi0[i];
            }else if (additionalNameOutput.Contains("EGA") && meson.CompareTo("Pi0")==0){
                bsmooth[i]      = bsmoothpp8TeVEGAPi0[i];
            }
        }
        else if(!energy.CompareTo("pPb_8TeV")){
            if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                bsmooth[i]      = bsmoothpPb8TeVINT7Pi0[i];
            }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0")==0){
                bsmooth[i]      = bsmoothpPb8TeVEG2Pi0[i];
            }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0")==0){
                bsmooth[i]      = bsmoothpPb8TeVEG1Pi0[i];
            }
            if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0Ratio")==0){
                bsmooth[i]      = bsmoothpPb8TeVINT7Pi0Ratio[i];
            }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0Ratio")==0){
                bsmooth[i]      = bsmoothpPb8TeVEG2Pi0Ratio[i];
            }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0Ratio")==0){
                bsmooth[i]      = bsmoothpPb8TeVEG1Pi0Ratio[i];
            }
        }
        else {
            cout << "smoothing not defined for " << energy.Data() << " returning!" << endl;
            return;
        }
    }


    // ***************************************************************************************************
    // ****************************** Initialize error vectors & graphs **********************************
    // ***************************************************************************************************

    Double_t* errorsNeg                     [nCuts];
    Double_t errorsNegCorr                  [nCuts][nPtBins];
    Double_t errorsNegSummed                [nPtBins];
    Double_t errorsNegCorrSummed            [nPtBins];
    Double_t errorsNegCorrMatSummed         [nPtBins];

    Double_t* errorsNegErr                  [nCuts];
    Double_t errorsNegErrCorr               [nCuts][nPtBins];
    Double_t errorsNegErrSummed             [nPtBins];
    Double_t errorsNegErrCorrSummed         [nPtBins];

    Double_t* errorsPos                     [nCuts];
    Double_t errorsPosCorr                  [nCuts][nPtBins];
    Double_t errorsPosSummed                [nPtBins];
    Double_t errorsPosCorrSummed            [nPtBins];
    Double_t errorsPosCorrMatSummed         [nPtBins];

    Double_t* errorsPosErr                  [nCuts];
    Double_t errorsPosErrSummed             [nPtBins];
    Double_t errorsPosErrCorr               [nCuts][nPtBins];
    Double_t errorsPosErrCorrSummed         [nPtBins];

    Double_t errorsMean                     [nCuts][nPtBins];
    Double_t errorsMeanCorr                 [nCuts][nPtBins];
    Double_t errorsMeanSummed               [nPtBins];
    Double_t errorsMeanCorrSummed           [nPtBins];
    Double_t errorsMeanCorrMatSummed        [nPtBins];

    Double_t errorsMeanErr                  [nCuts][nPtBins];
    Double_t errorsMeanErrCorr              [nCuts][nPtBins];
    Double_t errorsMeanErrSummed            [nPtBins];
    Double_t errorsMeanErrCorrSummed        [nPtBins];
    Double_t errorsMeanErrCorrMatSummed     [nPtBins];

    TGraphErrors* negativeErrors            [nCuts];
    TGraphErrors* positiveErrors            [nCuts];
    TGraphErrors* negativeErrorsCorr        [nCuts];
    TGraphErrors* positiveErrorsCorr        [nCuts];
    TGraphErrors* meanErrors                [nCuts];
    TGraphErrors* meanErrorsCorr            [nCuts];

    TGraphErrors* negativeErrorsSummed;
    TGraphErrors* positiveErrorsSummed;
    TGraphErrors* negativeErrorsCorrSummed;
    TGraphErrors* positiveErrorsCorrSummed;
    TGraphErrors* meanErrorsSummed;
    TGraphErrors* meanErrorsCorrSummed;
    TGraphErrors* meanErrorsCorrSummedIncMat;

    for (Int_t l = 0;l < nPtBins;l++){
        errorsPosSummed[l]              = 0.;
        errorsNegSummed[l]              = 0.;
        errorsMeanSummed[l]             = 0.;
        errorsPosCorrSummed[l]          = 0.;
        errorsNegCorrSummed[l]          = 0.;
        errorsMeanCorrSummed[l]         = 0.;
    }

    // ***************************************************************************************************
    // ****************************** Read & process data from file **************************************
    // ***************************************************************************************************
    TFile* fileErrorInput= new TFile(nameDataFileErrors.Data());

    for (Int_t i = 0;i < nCuts;i++){

        // read data
        TGraphAsymmErrors* graphPosErrors;
        TGraphAsymmErrors* graphNegErrors;
        TString nameGraphPos    = Form("%s_SystErrorRelPos_%s%s",meson.Data(),nameCutVariationSC[i].Data(),additionalNameOutput.Data()  );
        TString nameGraphNeg    = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),nameCutVariationSC[i].Data(),additionalNameOutput.Data()  );
        cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
        if ( nameCutVariationSC[i].CompareTo("ClusterEnergyScale") == 0 || nameCutVariationSC[i].CompareTo("Trigger") == 0 || 
            nameCutVariationSC[i].CompareTo("Efficiency") == 0 || nameCutVariationSC[i].CompareTo("MesonResolution") == 0 || 
            nameCutVariationSC[i].CompareTo("Secondary") == 0  || nameCutVariationSC[i].CompareTo("Purity") == 0
        ){
          nameGraphPos    = Form("%s_SystErrorRelPos_%s%s",meson.Data(),nameCutVariationSC[0].Data(),additionalNameOutput.Data()  );
          nameGraphNeg    = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),nameCutVariationSC[0].Data(),additionalNameOutput.Data()  );
        }

        graphPosErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
        graphNegErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());

        if (!graphPosErrors){
            nameGraphPos    = Form("%s_SystErrorRelPos_%s%s",meson.Data(),nameCutVariationSC[0].Data(),additionalNameOutput.Data()  );
            nameGraphNeg    = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),nameCutVariationSC[0].Data(),additionalNameOutput.Data()  );
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        }

        // take out offsets
        while (graphPosErrors->GetX()[0] < startPtSys){
            graphPosErrors->RemovePoint(0);
            graphNegErrors->RemovePoint(0);
        }
        cout << "****************************************************"<< endl;
        graphPosErrors->Print();
        cout << "****************************************************"<< endl;
        cout << "****************************************************"<< endl;
        cout << "****************************************************\n\n\n"<< endl;
        // Filling arrays
        if (i == 0) {
            ptBins      = graphNegErrors->GetX();
            ptBinsErr   = graphNegErrors->GetEXhigh();
        }
        errorsNeg[i]    = graphNegErrors->GetY();
        errorsNegErr[i] = graphNegErrors->GetEYhigh();
        errorsPos[i]    = graphPosErrors->GetY();
        errorsPosErr[i] = graphPosErrors->GetEYhigh();

        cout << nameCutVariationSC[i].Data() << endl;
        // Averaging of upper and lower errors
        CalculateMeanSysErr(errorsMean[i], errorsMeanErr[i], errorsPos[i], errorsNeg[i], nPtBins);
        // Automatic smoothing of 0 bins according to adjoining bins
        CorrectSystematicErrorsWithMean(errorsPos[i],errorsPosErr[i], errorsPosCorr[i], errorsPosErrCorr[i], nPtBins);
        CorrectSystematicErrorsWithMean(errorsNeg[i],errorsNegErr[i], errorsNegCorr[i], errorsNegErrCorr[i], nPtBins);
        CorrectSystematicErrorsWithMean(errorsMean[i], errorsMeanErr[i], errorsMeanCorr[i], errorsMeanErrCorr[i], nPtBins);

        // Routing for manual smoothing of systematic errors
        // ATTTENTION! you have to do this manually for each data set/trigger never trust the values mentioned here
        if (bsmooth[i]){
           // manual smoothing for cluster shape errors - variation 0
            if (nameCutVariationSC[i].CompareTo("ClusterM02")==0 ){
                cout << "Cluster M02 smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error          = 0;
                    if(!energy.CompareTo("8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 4.3;
                        }else if (additionalNameOutput.Contains("EMC7") && meson.CompareTo("Pi0")==0){
                            error   = 0.1+0.0493628*ptBins[k];
                        }else if (additionalNameOutput.Contains("EGA") && meson.CompareTo("Pi0")==0){
                            error   = 0.605507+0.0259021*ptBins[k]+-6.9671e-05*pow(ptBins[k],2);
                        }
                    }
                    else if(!energy.CompareTo("pPb_8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 2.07;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0")==0){
                            error   = -6.08732+0.919543*ptBins[k]+-0.0336618*pow(ptBins[k],2)+1.30913e-05*pow(ptBins[k],4);
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0")==0){
                            error   = 0.0915347+0.0368243*ptBins[k]+-0.000222092*pow(ptBins[k],2)+3.05061e-09*pow(ptBins[k],4);
                        }
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 3.0;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 3.70466+-0.076748*ptBins[k]+0.00226496*pow(ptBins[k],2);
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 1.02+0.0368243*ptBins[k]+-0.000222092*pow(ptBins[k],2)+3.05061e-09*pow(ptBins[k],4);
                        }
                    }
                    errorsMean[i][k]        = error;
                    errorsMeanErr[i][k]     = 0.01*error;
                    errorsMeanCorr[i][k]    = error;
                    errorsMeanErrCorr[i][k] = 0.01*error;
                }
            }

            // manual smoothing for cluster matching errors - variation 1
            if (nameCutVariationSC[i].CompareTo("ClusterTrackMatchingCalo")==0 ){
                cout << "Cluster track matching smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error          = 0;
                    if(!energy.CompareTo("8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = -1.32156+0.150567*ptBins[k];
                        }else if (additionalNameOutput.Contains("EMC7") && meson.CompareTo("Pi0")==0){
                            error   = -1.10329+0.0952144*ptBins[k];
                        }else if (additionalNameOutput.Contains("EGA") && meson.CompareTo("Pi0")==0){
                            error   = 51.3161+-50.5636/pow(1.00039,ptBins[k]);
                        }
                    }
                    else if(!energy.CompareTo("pPb_8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = -1.32156+0.150567*ptBins[k];
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0")==0){
                            error   = -1.10329+0.0952144*ptBins[k];
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0")==0){
                            error   = 51.3161+-50.5636/pow(1.00039,ptBins[k]);
                        }
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 51.3161+-50.5636/pow(1.00039,ptBins[k]);
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.5*(51.3161+-50.5636/pow(1.00039,ptBins[k]));
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.5*(51.3161+-50.5636/pow(1.00039,ptBins[k]));
                        }
                    }
                    errorsMean[i][k]        = error;
                    errorsMeanErr[i][k]     = 0.01*error;
                    errorsMeanCorr[i][k]    = error;
                    errorsMeanErrCorr[i][k] = 0.01*error;
                }
            }

            // manual smoothing for energy calibration errors - variation 2
            if (nameCutVariationSC[i].CompareTo("ClusterNonLinearity")==0 ){
                cout << "Cluster non linearity smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error          = 0;
                    if(!energy.CompareTo("8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 2.5;
                        }else if (additionalNameOutput.Contains("EMC7") && meson.CompareTo("Pi0")==0){
                            error   = 2.3;
                        }else if (additionalNameOutput.Contains("EGA") && meson.CompareTo("Pi0")==0){
                            error   = 1.26989+0.021*ptBins[k];
                        }
                    }
                    else if(!energy.CompareTo("pPb_8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 2.3;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0")==0){
                            error   = 2.3;
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0")==0){
                            error   = 2.3;
                        }
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.5;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.5;
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.5;
                        }
                    }
                    errorsMean[i][k]        = error;
                    errorsMeanErr[i][k]     = 0.01*error;
                    errorsMeanCorr[i][k]    = error;
                    errorsMeanErrCorr[i][k] = 0.01*error;
                }
            }

            // manual smoothing for cell aggregation - variation 3
            if (nameCutVariationSC[i].CompareTo("CellMinE")==0 ){
                cout << "cell Emin error smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error          = 0;
                    if(!energy.CompareTo("8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 1.7;
                        }else if (additionalNameOutput.Contains("EMC7") && meson.CompareTo("Pi0")==0){
                            error   = 1.7;
                        }else if (additionalNameOutput.Contains("EGA") && meson.CompareTo("Pi0")==0){
                            error   = 1.7;
                        }
                    }
                    else if(!energy.CompareTo("pPb_8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 2.5;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0")==0){
                            error   = 2.5;
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0")==0){
                            error   = 2.5;
                        }
                    }
                    errorsMean[i][k]        = error;
                    errorsMeanErr[i][k]     = 0.01*error;
                    errorsMeanCorr[i][k]    = error;
                    errorsMeanErrCorr[i][k] = 0.01*error;
                }
            }

            // manual smoothing for cell time uncertainties - variation 4
            if (nameCutVariationSC[i].CompareTo("CellTiming") == 0){
                cout << "Cell time smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error          = 0;
                    if(!energy.CompareTo("8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 3;
                        }else if (additionalNameOutput.Contains("EMC7") && meson.CompareTo("Pi0")==0){
                            error   = -3.10171+0.248002*ptBins[k]+-0.0031808*pow(ptBins[k],2);
                        }else if (additionalNameOutput.Contains("EGA") && meson.CompareTo("Pi0")==0){
                            error   = -0.838975+0.0990527*ptBins[k]+-0.000878754*pow(ptBins[k],2)+1.50378e-08*pow(ptBins[k],4);
                        }
                    }
                    else if(!energy.CompareTo("pPb_8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 1.06207+-399.636/pow(1.47445,ptBins[k]);
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0")==0){
                            error   = 0.410246+-4.1216/pow(1.16171,ptBins[k]);
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0")==0){
                            error   = 0.824983+-1.47897/pow(1.04395,ptBins[k]);
                        }
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 1.1;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.5;
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.4;
                        }
                    }
                    errorsMean[i][k]        = error;
                    errorsMeanErr[i][k]     = 0.01*error;
                    errorsMeanCorr[i][k]    = error;
                    errorsMeanErrCorr[i][k] = 0.01*error;
                }
            }

            // manual smoothing for Material infront of EMC - variation 5
            if (nameCutVariationSC[i].CompareTo("ClusterMaterialTRD")==0 ){
                cout << "Material smoothing" << endl;
                Double_t error                  = 4.24; //(3% for TRD mat, 3% for TOF mat added in quadrature)
                if(meson.Contains("Ratio"))
                    error = 0.0;
                for (Int_t k = 0;k < nPtBins;k++){
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }

            // manual smoothing for meson resolution errors - variation 6
            if (nameCutVariationSC[i].CompareTo("MesonResolution")==0 ){
                cout << "Meson resolution error smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error          = 0;
                    if(!energy.CompareTo("8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 7+0.02*(ptBins[k]-20);
                        }else if (additionalNameOutput.Contains("EMC7") && meson.CompareTo("Pi0")==0){
                            error   = 7+0.02*(ptBins[k]-20);
                        }else if (additionalNameOutput.Contains("EGA") && meson.CompareTo("Pi0")==0){
                            error   = 7+0.02*(ptBins[k]-20);
                        }
                    }
                    else if(!energy.CompareTo("pPb_8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 7+0.02*(ptBins[k]-20);
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0")==0){
                            error   = 7+0.02*(ptBins[k]-20);
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0")==0){
                            error   = 7+0.02*(ptBins[k]-20);
                        }
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.6;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.6;
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.6;
                        }
                    }
                    errorsMean[i][k]        = error;
                    errorsMeanErr[i][k]     = 0.01*error;
                    errorsMeanCorr[i][k]    = error;
                    errorsMeanErrCorr[i][k] = 0.01*error;
                }
            }

            // manual smoothing for energy scale errors (derived from mass difference MC & Data) - variation 7
            if (nameCutVariationSC[i].CompareTo("ClusterEnergyScale")==0 ){
                cout << "Cluster non linearity residual offset smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error          = 0;
                    if(!energy.CompareTo("8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 2.1;
                        }else if (additionalNameOutput.Contains("EMC7") && meson.CompareTo("Pi0")==0){
                            error   = 2.1;
                        }else if (additionalNameOutput.Contains("EGA") && meson.CompareTo("Pi0")==0){
                            error   = 2.1;
                        }
                    }
                    else if(!energy.CompareTo("pPb_8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 2.3;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0")==0){
                            error   = 2.3;
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0")==0){
                            error   = 2.3;
                        }
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.84;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.84;
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.84;
                        }
                    }
                    errorsMean[i][k]        = error;
                    errorsMeanErr[i][k]     = 0.01*error;
                    errorsMeanCorr[i][k]    = error;
                    errorsMeanErrCorr[i][k] = 0.01*error;
                }
            }

            // manual smoothing for Trigger normalization uncertainties - variation 8
            if (nameCutVariationSC[i].CompareTo("Trigger") == 0){
                cout << "Trigger rejection factor uncertainty smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error          = 0;
                    if(!energy.CompareTo("8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 0.;
                        }else if (additionalNameOutput.Contains("EMC7") && meson.CompareTo("Pi0")==0){
                            error   = TMath::Sqrt(1.22*1.22); // V2 clusterizer 3.21% trigger uncertainty, 2% pileup
                        }else if (additionalNameOutput.Contains("EGA") && meson.CompareTo("Pi0")==0){
                            error   = TMath::Sqrt(1.22*1.22+2.99*2.99); // V2 clusterizer 8.6, 2% pileup (trigger EG1/EG2: 5.99%)
                        }
                    }
                    else if(!energy.CompareTo("pPb_8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 0.;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0")==0){
                            error   = TMath::Sqrt(1.57*1.57); // V2 clusterizer 6.17% trigger uncertainty, 2% pileup (trigger: EG2/EMC7 5.27%)
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0")==0){
                            error   = TMath::Sqrt(1.57*1.57+1.05*1.05); // V2 clusterizer 8.6, 2% pileup (trigger EG1/EG2: 5.99%)
                        }
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0Ratio")==0){
                            error   = TMath::Sqrt(1.57*1.57+1.22*1.22); //
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0Ratio")==0){
                            error   = TMath::Sqrt(1.57*1.57+1.05*1.05+1.22*1.22+2.99*2.99); // 
                        }
                    }
                    errorsMean[i][k]        = error;
                    errorsMeanErr[i][k]     = 0.01*error;
                    errorsMeanCorr[i][k]    = error;
                    errorsMeanErrCorr[i][k] = 0.01*error;
                }
            }

            // manual smoothing for Efficiency uncertainties - variation 9
            if (nameCutVariationSC[i].CompareTo("Efficiency") == 0){
                // from comparison of power data and MC spectrum
                cout << "Efficiency smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error          = 0;
                    if(!energy.CompareTo("8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 3.;
                        }else if (additionalNameOutput.Contains("EMC7") && meson.CompareTo("Pi0")==0){
                            error   = 3.;
                        }else if (additionalNameOutput.Contains("EGA") && meson.CompareTo("Pi0")==0){
                            error   = 3.;
                        }
                    }
                    else if(!energy.CompareTo("pPb_8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 3.;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0")==0){
                            error   = 3.;
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0")==0){
                            error   = 3.;
                        }
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.9;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.9;
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.9;
                        }
                    }
                    errorsMean[i][k]        = error;
                    errorsMeanErr[i][k]     = 0.01*error;
                    errorsMeanCorr[i][k]    = error;
                    errorsMeanErrCorr[i][k] = 0.01*error;
                }
            }

            // manual smoothing for Secondary errors - variation 10
            if (nameCutVariationSC[i].CompareTo("Secondary")==0 ){
                cout << "Secondary Pi0 correction error smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error          = 0;
                    if(!energy.CompareTo("8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 3.8; // for now half of the correction
                        }else if (additionalNameOutput.Contains("EMC7") && meson.CompareTo("Pi0")==0){
                            error   = 3.8;
                        }else if (additionalNameOutput.Contains("EGA") && meson.CompareTo("Pi0")==0){
                            error   = 3.8;
                        }
                    }
                    else if(!energy.CompareTo("pPb_8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 4.2;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0")==0){
                            error   = 4.2;
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0")==0){
                            error   = 4.2;
                        }
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.75;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.75;
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.75;
                        }
                    }
                    errorsMean[i][k]        = error;
                    errorsMeanErr[i][k]     = 0.01*error;
                    errorsMeanCorr[i][k]    = error;
                    errorsMeanErrCorr[i][k] = 0.01*error;
                }
            }

            // manual smoothing for Purity correction errors - variation 11
            if (nameCutVariationSC[i].CompareTo("Purity")==0 ){
                cout << "Purity correction uncertainty error smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error          = 0;
                    if(!energy.CompareTo("8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 1.9;
                        }else if (additionalNameOutput.Contains("EMC7") && meson.CompareTo("Pi0")==0){
                            error   = 1.9;
                        }else if (additionalNameOutput.Contains("EGA") && meson.CompareTo("Pi0")==0){
                            error   = 1.9;
                        }
                    }
                    else if(!energy.CompareTo("pPb_8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 1.5;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0")==0){
                            error   = 1.5;
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0")==0){
                            error   = 1.5;
                        }
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.2;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.2;
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 0.2;
                        }
                    }
                    errorsMean[i][k]        = error;
                    errorsMeanErr[i][k]     = 0.01*error;
                    errorsMeanCorr[i][k]    = error;
                    errorsMeanErrCorr[i][k] = 0.01*error;
                }
            }

            // manual smoothing for Distance to Bad Channel correction errors - variation 12
            if (nameCutVariationSC[i].CompareTo("DistanceBadChannel")==0 ){
                cout << "Distance to Bad Channel uncertainty error smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error          = 0;
                    if(!energy.CompareTo("8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 2.0;
                        }else if (additionalNameOutput.Contains("EMC7") && meson.CompareTo("Pi0")==0){
                            error   = 1.5;
                        }else if (additionalNameOutput.Contains("EGA") && meson.CompareTo("Pi0")==0){
                            error   = 1.43573-1.01514/pow(1.00796,ptBins[k]);
                        }
                    }
                    else if(!energy.CompareTo("pPb_8TeV")){
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
                            error   = 1.3;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0")==0){
                            error   = 27.842+-28.4243/pow(1.00231,ptBins[k]);
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0")==0){
                            error   = 0.628326+-0.00693552*ptBins[k]+9.80999e-05*pow(ptBins[k],2);
                        }
                        if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 1.2;
                        }else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 1.2;
                        }else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0Ratio")==0){
                            error   = 1.2;
                        }
                    }
                    errorsMean[i][k]        = error;
                    errorsMeanErr[i][k]     = 0.01*error;
                    errorsMeanCorr[i][k]    = error;
                    errorsMeanErrCorr[i][k] = 0.01*error;
                }
            }

        } else {
            for (Int_t k = 0;k < nPtBins;k++){
                errorsMeanErr[i][k]         = 0.03;
                errorsMeanErrCorr[i][k]     = 0.03;
            }
        }
        // Quadratic sum of errors except material error infront of EMCal & inner material
        if (!(nameCutVariationSC[i].CompareTo("ClusterMaterialTRD")==0)){
            cout << "errors added quadratically" << endl;
            for (Int_t l = 0;l < nPtBins;l++){
                errorsPosSummed[l]      = errorsPosSummed[l]+pow(errorsPos[i][l],2);
                errorsNegSummed[l]      = errorsNegSummed[l]+ pow(errorsNeg[i][l],2);
                errorsMeanSummed[l]     = errorsMeanSummed[l]+ pow(errorsMean[i][l],2);
                errorsPosCorrSummed[l]  = errorsPosCorrSummed[l]+pow(errorsPosCorr[i][l],2);
                errorsNegCorrSummed[l]  = errorsNegCorrSummed[l] +pow(errorsNegCorr[i][l],2);
                errorsMeanCorrSummed[l] = errorsMeanCorrSummed[l]+ pow(errorsMeanCorr[i][l],2);
            }
        }
        // fill error graphs for plotting
        negativeErrors[i]       = new TGraphErrors(nPtBins,ptBins ,errorsNeg[i] ,ptBinsErr ,errorsNegErr[i] );
        meanErrors[i]           = new TGraphErrors(nPtBins,ptBins ,errorsMean[i] ,ptBinsErr ,errorsMeanErr[i] );
        positiveErrors[i]       = new TGraphErrors(nPtBins,ptBins ,errorsPos[i] ,ptBinsErr ,errorsPosErr[i] );
        negativeErrorsCorr[i]   = new TGraphErrors(nPtBins,ptBins ,errorsNegCorr[i] ,ptBinsErr ,errorsNegErrCorr[i] );
        meanErrorsCorr[i]       = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorr[i] ,ptBinsErr ,errorsMeanErrCorr[i] );
        positiveErrorsCorr[i]   = new TGraphErrors(nPtBins,ptBins ,errorsPosCorr[i] ,ptBinsErr ,errorsPosErrCorr[i] );

    }

    // Error for inner material budget
    Double_t errorMaterial = 0;

    // Calculate sqrt of summed errors for final errors, add material budget errors
    for (Int_t l = 0;l < nPtBins;l++){
        errorsPosSummed[l]              = pow(errorsPosSummed[l],0.5);
        errorsMeanSummed[l]             = pow(errorsMeanSummed[l],0.5);
        errorsPosErrSummed[l]           = errorsPosSummed[l]*0.001;
        errorsMeanErrSummed[l]          = errorsMeanSummed[l]*0.001;
        errorsNegSummed[l]              = -pow(errorsNegSummed[l],0.5);
        errorsNegErrSummed[l]           = errorsNegSummed[l]*0.001;

        // add EMCal material errors
        errorsPosCorrMatSummed[l]       = pow(errorsPosCorrSummed[l]+ pow(errorMaterial ,2.) + pow(errorsPosCorr[6][l],2) ,0.5);
        errorsMeanCorrMatSummed[l]      = pow(errorsMeanCorrSummed[l]+ pow(errorMaterial ,2.)+ pow(errorsMeanCorr[6][l],2),0.5);
        errorsNegCorrMatSummed[l]       = -pow(errorsNegCorrSummed[l]+ pow(errorMaterial ,2.)+ pow(errorsNegCorr[6][l],2),0.5);

        errorsPosCorrSummed[l]          = pow(errorsPosCorrSummed[l],0.5);
        errorsMeanCorrSummed[l]         = pow(errorsMeanCorrSummed[l],0.5);
        errorsPosErrCorrSummed[l]       = errorsPosCorrSummed[l]*0.001;
        errorsMeanErrCorrSummed[l]      = errorsMeanCorrSummed[l]*0.001;
        errorsMeanErrCorrMatSummed[l]   = errorsMeanCorrMatSummed[l]*0.001;
        errorsNegCorrSummed[l]          = -pow(errorsNegCorrSummed[l],0.5);
        errorsNegErrCorrSummed[l]       = errorsNegCorrSummed[l]*0.001;
    }

    // Create all other summed graphs
    cout << __LINE__ << endl;
    negativeErrorsSummed        = new TGraphErrors(nPtBins,ptBins ,errorsNegSummed ,ptBinsErr ,errorsNegErrSummed );
    negativeErrorsCorrSummed    = new TGraphErrors(nPtBins,ptBins ,errorsNegCorrSummed ,ptBinsErr ,errorsNegErrCorrSummed );
    positiveErrorsSummed        = new TGraphErrors(nPtBins,ptBins ,errorsPosSummed ,ptBinsErr ,errorsPosErrSummed );
    positiveErrorsCorrSummed    = new TGraphErrors(nPtBins,ptBins ,errorsPosCorrSummed ,ptBinsErr ,errorsPosErrCorrSummed );
    meanErrorsSummed            = new TGraphErrors(nPtBins,ptBins ,errorsMeanSummed ,ptBinsErr ,errorsMeanErrSummed );
    meanErrorsCorrSummed        = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSummed ,ptBinsErr ,errorsMeanErrCorrSummed );
    meanErrorsCorrSummedIncMat  = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrMatSummed ,ptBinsErr ,errorsMeanErrCorrMatSummed );

    cout << __LINE__ << endl;

    // Give legend position for plotting
    Double_t minXLegend     = 0.12;
    Double_t maxYLegend     = 0.95;
    Double_t widthLegend    = 0.25;
    if (numberCutStudies> 7)
        widthLegend         = 0.5;
    Double_t heightLegend   = 1.05* 0.035 * (numberCutStudies+3);
    if (numberCutStudies> 7)
        heightLegend        = 1.05* 0.035 * (numberCutStudies/2+1);

    // ***************************************************************************************************
    // ****************************** Plot all mean erros separately *************************************
    // ***************************************************************************************************

    TCanvas* canvasSysErrMean = new TCanvas("canvasSysErrMean","",200,10,1350,900);// gives the page size
    DrawGammaCanvasSettings( canvasSysErrMean, 0.08, 0.01, 0.015, 0.09);

        // create dummy histo
        TH2D *histo2DSysErrMean ;
        if (meson.Contains("Pi0") ){
            histo2DSysErrMean = new TH2D("histo2DSysErrMean", "", 20,startPtSys-2,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,0.,25.);
        } else {
            histo2DSysErrMean = new TH2D("histo2DSysErrMean", "", 20,startPtSys-2,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,0.,65.);
        }
        SetStyleHistoTH2ForGraphs( histo2DSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean systematic Err %", 0.03, 0.04, 0.03, 0.04,
                                1,0.9, 510, 510);
        histo2DSysErrMean->Draw();

        // create legend
        TLegend* legendMean = GetAndSetLegend2(minXLegend,maxYLegend-heightLegend,minXLegend+widthLegend,maxYLegend, 30);
        if (numberCutStudies> 7) legendMean->SetNColumns(2);

        for(Int_t i = 0;i< numberCutStudies ;i++){
            if(meson.Contains("Ratio") && i == 5){
                cout << "not drawing: " << nameCutVariation[i].Data() << endl;
                continue;
            }
            DrawGammaSetMarkerTGraphErr(meanErrors[i], markerStyle[i], 1.,color[i],color[i]);
            meanErrors[i]->Draw("pE0,csame");
            legendMean->AddEntry(meanErrors[i],nameCutVariation[i].Data(),"p");

        }
        legendMean->Draw();

        // plot labeling
        TLatex *labelMeson;
        labelMeson= new TLatex(0.95,0.885,labelMesonReco);
        SetStyleTLatex( labelMeson, 0.038,4,1,42,kTRUE,31);
        labelMeson->Draw();

        TLatex *labelCentrality = new TLatex(0.95,0.93,Form("%s",collisionSystem.Data() ));
        SetStyleTLatex( labelCentrality, 0.038,4,1,42,kTRUE,31);
        labelCentrality->Draw();

        TLatex *labelTrig= new TLatex(0.95,0.84,Form("%s %s",labelTriggerText.Data(),labelDataSetText.Data()));
        SetStyleTLatex( labelTrig, 0.038,4,1,42,kTRUE,31);
        labelTrig->Draw();

    canvasSysErrMean->Update();
    canvasSysErrMean->SaveAs(Form("SystematicErrorsCalculatedMergedCalo/SysMean_%s_%s%s_%s.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));

    delete canvasSysErrMean;

    // ***************************************************************************************************
    // ********************* Plot all mean erros separately after smoothing ******************************
    // ***************************************************************************************************    
    TCanvas* canvasNewSysErrMean = new TCanvas("canvasNewSysErrMean","",200,10,1350,900);// gives the page size
    DrawGammaCanvasSettings( canvasNewSysErrMean, 0.08, 0.01, 0.015, 0.09);

        // create dummy histo
        TH2D *histo2DNewSysErrMean ;
        if (meson.Contains("Pi0")){
            histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 20,startPtSys-2,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,25.);
        } else {
            histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 20,startPtSys-2,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,65.);
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
            if(meson.Contains("Ratio") && i == 5){
                cout << "not drawing: " << nameCutVariation[i].Data() << endl;
                continue;
            }
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], markerStyle[i], 1.,color[i],color[i]);
            meanErrorsCorr[i]->Draw("pX0,csame");
            legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");
        }
        DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
        meanErrorsCorrSummedIncMat->Draw("p,csame");
        legendMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
        legendMeanNew->Draw();

        // labeling
        labelMeson->Draw();
        labelCentrality->Draw();
        labelTrig->Draw();

    canvasNewSysErrMean->Update();
    canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedMergedCalo/SysMeanNewWithMean_%s_%s%s_%s.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
    
    // ***************************************************************************************************
    // ********************* Plot unsmoothed errors with fits ********************************************
    // ***************************************************************************************************    
    for (Int_t cut =0 ;cut < numberCutStudies;cut++ ){
        
        canvasNewSysErrMean->cd();
        histo2DNewSysErrMean->Draw();

            if (bsmooth[cut]) continue;
            cout <<endl << endl<<  "variation: " << cut << " \t"<< nameCutVariation[cut].Data() << endl;
            Double_t minPt = startPtSys;
            Double_t maxPt = ptBins[nPtBins-2]+1;

            TF1* pol0 = new TF1("pol0","[0]",minPt,maxPt);//
            TF1* pol1 = new TF1("pol1","[0]+[1]*x",minPt,maxPt);//
            TF1* pol2 = new TF1("pol2","[0]+[1]*x+[2]*x*x",minPt,maxPt);//
            TF1* pol4 = new TF1("pol4","[0]+[1]*x+[2]*x*x+[3]*x*x*x*x",minPt,maxPt);//
            TF1* bla  = new TF1("bla","[0]+[1]/pow([2],x)",minPt,maxPt);

            bla->SetParameter(0,1.5);
            bla->SetParameter(1,50);
            bla->SetParameter(2,5);
            pol4->SetParLimits(3,0,10);

            meanErrorsCorr[cut]->Fit(pol4,"QNRMEX0+","",minPt,maxPt);
            meanErrorsCorr[cut]->Fit(pol2,"QNRMEX0+","",minPt,maxPt);
            meanErrorsCorr[cut]->Fit(pol1,"QNRMEX0+","",minPt,maxPt);
            meanErrorsCorr[cut]->Fit(pol0,"QNRMEX0+","",minPt,maxPt);
            meanErrorsCorr[cut]->Fit(bla,"QNRMEX0+","",minPt,maxPt);
            pol4->SetLineColor(kRed+2);
            pol2->SetLineColor(kBlue+2);
            pol1->SetLineColor(kCyan+2);
            pol0->SetLineColor(kBlack);
            bla->SetLineColor(kMagenta+2);

            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[cut], 20+cut, 1.,color[cut],color[cut]);
            meanErrorsCorr[cut]->Draw("p,csame");
            pol4->Draw("same");
            pol2->Draw("same");
            pol1->Draw("same");
            pol0->Draw("same");
            bla->Draw("same");

            cout << "RED pol4:" << endl;
            cout <<"pol4:\t"<< pol4->GetParameter(0)<<"+"<<pol4->GetParameter(1)<<"*ptBins[k]+"<<pol4->GetParameter(2)<<"*pow(ptBins[k],2)+"<<pol4->GetParameter(3)<<"*pow(ptBins[k],4)" << endl;

            cout << "BLUE pol2:" << endl;
            cout <<"pol2:\t"<< pol2->GetParameter(0)<<"+"<<pol2->GetParameter(1)<<"*ptBins[k]+"<<pol2->GetParameter(2)<<"*pow(ptBins[k],2)" << endl;

            cout << "CYAN pol1:" << endl;
            cout <<"pol1:\t"<< pol1->GetParameter(0)<<"+"<<pol1->GetParameter(1)<<"*ptBins[k]" << endl;

            cout << "BLACK pol0:" << endl;
            cout <<"pol0:\t"<< pol0->GetParameter(0) << endl;

            cout << "MAGENTA 1/power:" << endl;
            cout <<"1/power:\t"<< bla->GetParameter(0) << "+" << bla->GetParameter(1) << "/pow(" << bla->GetParameter(2)<< ",ptBins[k])"<< endl;

        canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedMergedCalo/SysMeanNewWithMeanSingle_%s_%s%s_%s_Variation%d.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),cut,suffix.Data()));
    }

    // ***************************************************************************************************
    // ********************* Create output files with errors *********************************************
    // ***************************************************************************************************    
    const char *SysErrDatname = Form("SystematicErrorsCalculatedMergedCalo/SystematicErrorEMCmerged_%s_%s%s_%s.dat",meson.Data(),energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data());
    fstream SysErrDat;
    cout << SysErrDatname << endl;
    SysErrDat.open(SysErrDatname, ios::out);
    for (Int_t l=0;l< nPtBins;l++){
        SysErrDat << ptBins[l] << "\t" <<errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
    }

    SysErrDat.close();

    const char *SysErrDatnameMean = Form("SystematicErrorsCalculatedMergedCalo/SystematicErrorAveragedEMCmerged_%s_%s%s_%s.dat",meson.Data(),energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data());
    fstream SysErrDatAver;
    cout << SysErrDatnameMean << endl;
    SysErrDatAver.open(SysErrDatnameMean, ios::out);
    for (Int_t l=0;l< nPtBins;l++){
        SysErrDatAver << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
    }

    const char *SysErrDatnameMeanSingleErr = Form("SystematicErrorsCalculatedMergedCalo/SystematicErrorAveragedSingleEMCmerged_%s_%s%s_%s.dat",meson.Data(),energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data());
    fstream SysErrDatAverSingle;
    cout << SysErrDatnameMeanSingleErr << endl;
    SysErrDatAverSingle.open(SysErrDatnameMeanSingleErr, ios::out);
    SysErrDatAverSingle << "Pt bin\t" ; 
    for (Int_t i= 0; i< numberCutStudies; i++){
        SysErrDatAverSingle << nameCutVariationSC[i] << "\t";
    }
    SysErrDatAverSingle << endl; 
    for (Int_t l=0;l< nPtBins;l++){
        SysErrDatAverSingle << ptBins[l] << "\t";
        for (Int_t i= 0; i< numberCutStudies; i++){
            SysErrDatAverSingle << errorsMeanCorr[i][l] << "\t";
        }  
        
        SysErrDatAverSingle << errorsMeanCorrMatSummed[l] << endl;
    }
    
    
    SysErrDatAverSingle.close();
    
        
}