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

void FinaliseSystematicErrorsConvCalo_pp(   TString nameDataFileErrors    = "", 
                                            TString energy                = "", 
                                            TString meson                 = "", 
                                            Int_t numberOfPtBins          = 1, 
                                            Int_t numberCutStudies        = 1, 
                                            Float_t startPtSys            = 0, 
                                            TString additionalName        = "pp", 
                                            TString additionalNameOutput  = "", 
                                            TString suffix                = "eps"
                                        ){

    // ***************************************************************************************************
    // ****************************** General style settings *********************************************
    // ***************************************************************************************************
    StyleSettingsThesis();	
    SetPlotStyle();
    
    // ***************************************************************************************************
    // ****************************** Create output directory ********************************************
    // ***************************************************************************************************    
    gSystem->Exec("mkdir -p SystematicErrorsCalculatedConvCalo");
    
    // ***************************************************************************************************
    // ***************************** labeling and color settings *****************************************
    // ***************************************************************************************************
    TString date                            = ReturnDateString();
    TString dateForOutput                   = ReturnDateStringForOutput();
    TString collisionSystem                 = ReturnFullCollisionsSystem(energy);
    TString energyForOutput                 = energy;
    energyForOutput.ReplaceAll(".","_");
    
    Color_t color[20] = {860,894,807,880,418,403,802,923,634,432,404,435,420,407,416,830,404,608,kCyan-2,1};
    Color_t markerStyle[20] = {24,21,22,23,20,25,26,27,28,29,30,31,32,33,24,21,22,23,20,25};
    
    // ***************************************************************************************************
    // ******************************* general variable definition  **************************************
    // ***************************************************************************************************
    Int_t 	numberOfEntriesPos = 0;
    Int_t 	numberOfEntriesNeg = 0;
    const Int_t nPtBins = numberOfPtBins;
    const Int_t nCuts = numberCutStudies;
    Double_t* ptBins;
    Double_t* ptBinsErr;
    TString nameCutVariation[19];
    TString nameCutVariationSC[19];
    
    TString nameCutVariation2760GeV[19] = {"Yield extraction", "dE/dx e-line", "dE/dx #pi-line", "TPC cluster", "Single e^{#pm} #it{p}_{T}",
                                            "2D #chi^{2} #gamma, #psi_{pair} #gamma","2D q_{T}","#alpha meson", "#varphi_{conv}", "min E_{cluster}",
                                            "min # cells", "clu. energy calibration", "V0 tr. match. to cl.", "M_{02}", "cell time",
                                            "Mat. infront of EMCal", "Trigger normalization", "Efficiency", "Yield extraction #pi^{0}"};      
    TString nameCutVariationSC2760GeV[19] = {"YieldExtraction", "dEdxE", "dEdxPi", "TPCCluster", "SinglePt",
                                            "Chi2", "Qt", "Alpha", "ConvPhi", "ClusterMinEnergy",
                                            "ClusterNCells", "ClusterNonLinearity", "ClusterTrackMatching", "ClusterM02", "CellTiming",
                                            "ClusterMaterialTRD", "Trigger", "Efficiency", "YieldExtraction"};
    if (meson.CompareTo("EtaToPi0") == 0){
        nameCutVariation2760GeV[0]          = "Yield extraction #eta";
    }

    
    if (energy.CompareTo("2.76TeV") == 0) {
        for (Int_t i = 0; i < numberCutStudies; i++){
            nameCutVariation[i] = nameCutVariation2760GeV[i];
            nameCutVariationSC[i] = nameCutVariationSC2760GeV[i];
        }
    }
    
    // ***************************************************************************************************
    // ******************************** Booleans for smoothing *******************************************
    // ***************************************************************************************************
    Bool_t bsmooth[19]              = { 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0,
                                        0, 0, 0, 0 };
    // minimum bias trigger                      
    Bool_t bsmoothMBPi0[19]         = { 0, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    Bool_t bsmoothMBEta[19]         = { 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    Bool_t bsmoothMBEtaToPi0[19]    = { 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    // INT7 trigger
    Bool_t bsmoothINT7Pi0[19]       = { 0, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    Bool_t bsmoothINT7Eta[19]       = { 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    Bool_t bsmoothINT7EtaToPi0[19]  = { 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    // EMC1 trigger
    Bool_t bsmoothEMC1Pi0[19]       = { 0, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    Bool_t bsmoothEMC1Eta[19]       = { 0, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    Bool_t bsmoothEMC1EtaToPi0[19]  = { 0, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    // EMC7 trigger
    Bool_t bsmoothEMC7Pi0[19]       = { 0, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    Bool_t bsmoothEMC7Eta[19]       = { 0, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    Bool_t bsmoothEMC7EtaToPi0[19]  = { 0, 1, 1, 1, 1,
                                        1, 1, 1, 1, 0,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    // EG2 trigger
    Bool_t bsmoothEG2Pi0[19]        = { 0, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    Bool_t bsmoothEG2Eta[19]        = { 0, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    Bool_t bsmoothEG2EtaToPi0[19]   = { 0, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    // EG1 trigger
    Bool_t bsmoothEG1Pi0[19]        = { 0, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    Bool_t bsmoothEG1Eta[19]        = { 0, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
    Bool_t bsmoothEG1EtaToPi0[19]   = { 0, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1,
                                        1, 1, 1, 0 };
                          
    for (Int_t i = 0; i < numberCutStudies; i++){
        if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Pi0")==0){ 
            bsmooth[i] = bsmoothMBPi0[i];  
        } else if (additionalNameOutput.CompareTo("INT7") == 0 && meson.CompareTo("Pi0")==0){ 
            bsmooth[i] = bsmoothINT7Pi0[i];  
        } else if (additionalNameOutput.CompareTo("EMC1") == 0 && meson.CompareTo("Pi0")==0){ 
            bsmooth[i] = bsmoothEMC1Pi0[i];  
        } else if (additionalNameOutput.CompareTo("EMC7") == 0 && meson.CompareTo("Pi0")==0){ 
            bsmooth[i] = bsmoothEMC7Pi0[i];  
        } else if (additionalNameOutput.CompareTo("EG2") == 0 && meson.CompareTo("Pi0")==0){ 
            bsmooth[i] = bsmoothEG2Pi0[i];  
        } else if (additionalNameOutput.CompareTo("EG1") == 0 && meson.CompareTo("Pi0")==0){ 
            bsmooth[i] = bsmoothEG1Pi0[i];  
        } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Eta")==0){ 
            bsmooth[i] = bsmoothMBEta[i];  
        } else if (additionalNameOutput.CompareTo("INT7") == 0 && meson.CompareTo("Eta")==0){ 
            bsmooth[i] = bsmoothINT7Eta[i];  
        } else if (additionalNameOutput.CompareTo("EMC1") == 0 && meson.CompareTo("Eta")==0){ 
            bsmooth[i] = bsmoothEMC1Eta[i];  
        } else if (additionalNameOutput.CompareTo("EMC7") == 0 && meson.CompareTo("Eta")==0){ 
            bsmooth[i] = bsmoothEMC7Eta[i];  
        } else if (additionalNameOutput.CompareTo("EG2") == 0 && meson.CompareTo("Eta")==0){ 
            bsmooth[i] = bsmoothEG2Eta[i];  
        } else if (additionalNameOutput.CompareTo("EG1") == 0 && meson.CompareTo("Eta")==0){ 
            bsmooth[i] = bsmoothEG1Eta[i];  
        } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("EtaToPi0")==0){ 
            bsmooth[i] = bsmoothMBEtaToPi0[i];  
        } else if (additionalNameOutput.CompareTo("INT7") == 0 && meson.CompareTo("EtaToPi0")==0){ 
            bsmooth[i] = bsmoothINT7EtaToPi0[i];  
        } else if (additionalNameOutput.CompareTo("EMC1") == 0 && meson.CompareTo("EtaToPi0")==0){ 
            bsmooth[i] = bsmoothEMC1EtaToPi0[i];  
        } else if (additionalNameOutput.CompareTo("EMC7") == 0 && meson.CompareTo("EtaToPi0")==0){ 
            bsmooth[i] = bsmoothEMC7EtaToPi0[i];  
        } else if (additionalNameOutput.CompareTo("EG2") == 0 && meson.CompareTo("EtaToPi0")==0){ 
            bsmooth[i] = bsmoothEG2EtaToPi0[i];  
        } else if (additionalNameOutput.CompareTo("EG1") == 0 && meson.CompareTo("EtaToPi0")==0){ 
            bsmooth[i] = bsmoothEG1EtaToPi0[i];  
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
    
    for (Int_t l = 0; l < nPtBins; l++){
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
    
    for (Int_t i = 0; i < nCuts; i++){
        
        // read data
        TGraphAsymmErrors* graphPosErrors;
        TGraphAsymmErrors* graphNegErrors;
        // YieldExtraction - 0, Trigger - 16, Efficiency - 17
        if (i == 0 || i == 16 || i == 17){ // special treatment for Yield extraction error and calculated erros
            TString nameGraphPos    = "";
            TString nameGraphNeg    = "";
            if ( meson.CompareTo("EtaToPi0") != 0 ){
                nameGraphPos    = Form("%s_SystErrorRelPos_YieldExtraction_%s",meson.Data(),additionalName.Data() );
                nameGraphNeg    = Form("%s_SystErrorRelNeg_YieldExtraction_%s",meson.Data(),additionalName.Data() );                
            } else {
                nameGraphPos    = Form("Eta_SystErrorRelPos_YieldExtraction_%s",additionalName.Data() );
                nameGraphNeg    = Form("Eta_SystErrorRelNeg_YieldExtraction_%s",additionalName.Data() );
            }
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else if (i == 14 && meson.CompareTo("EtaToPi0") == 0 ){ // special treatment for Yield extraction error and calculated erros
            TString nameGraphPos    = Form("Eta_SystErrorRelPos_YieldExtraction_%s",additionalName.Data() );
            TString nameGraphNeg    = Form("Eta_SystErrorRelNeg_YieldExtraction_%s",additionalName.Data() );
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else if (i == 18){ // special treatment for eta to pi0 ratio    
            TString nameGraphPos    = Form("Pi0EtaBinning_SystErrorRelPos_YieldExtraction_%s",additionalName.Data() );
            TString nameGraphNeg    = Form("Pi0EtaBinning_SystErrorRelNeg_YieldExtraction_%s",additionalName.Data() );
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else { // read graphs from input file
            TString nameGraphPos    = Form("%s_SystErrorRelPos_%s%s",meson.Data(),nameCutVariationSC[i].Data(),additionalNameOutput.Data()  );
            TString nameGraphNeg    = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),nameCutVariationSC[i].Data(),additionalNameOutput.Data()  );
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        }
        
        // take out offsets
        while (graphPosErrors->GetX()[0] < startPtSys){
            graphPosErrors->RemovePoint(0);
            graphNegErrors->RemovePoint(0);
        }
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
            // manual smoothing for Yield extraction errors - variation 0
            if  (nameCutVariationSC[i].CompareTo("YieldExtraction") == 0){
                if ( meson.CompareTo("Pi0") == 0){
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (ptBins[k] < 1) continue;
                        Double_t error          = 2.;
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = error*0.01;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = error*0.01;
                    }
                } else {
                    for (Int_t k = 0; k < nPtBins; k++){
                        if (ptBins[k] < 2) continue;
                        Double_t error          = 8;
                        if (additionalNameOutput.CompareTo("")==0 || additionalNameOutput.CompareTo("INT7")==0){
                            error               = 10;
                        }  
                        
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = error*0.01;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = error*0.01;
                    }
                }
            }
            // manual smoothing for dEdx electron line errors - variation 1
            if (nameCutVariationSC[i].CompareTo("dEdxE")==0 ){
                if ( meson.CompareTo("Pi0")== 0){
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0.38;
                        if (additionalNameOutput.CompareTo("")==0 || additionalNameOutput.CompareTo("INT7")==0){
                            error = 0.38;
                        } else if (additionalNameOutput.CompareTo("EMC7")==0){
                            error = 0.7;    
                        } else if (additionalNameOutput.CompareTo("EMC1")==0 || additionalNameOutput.CompareTo("EG2")==0 ){
                            error = 1.5;
                        } else if (additionalNameOutput.CompareTo("EG1")==0){
                            error = 2.5;
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                } else if (meson.CompareTo("Eta") == 0) {
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0.;
                        if (additionalNameOutput.CompareTo("")==0 || additionalNameOutput.CompareTo("INT7")==0){
                            error = 2.;
                        } else if (additionalNameOutput.CompareTo("EMC7")==0){
                            error = 2.;    
                        } else if (additionalNameOutput.CompareTo("EMC1")==0 || additionalNameOutput.CompareTo("EG2")==0 ){
                            error = 2.;
                        } else if (additionalNameOutput.CompareTo("EG1")==0){
                            error = 2.;
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                } else {
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0.;
                        if (additionalNameOutput.CompareTo("")==0 || additionalNameOutput.CompareTo("INT7")==0){
                            error = TMath::Sqrt(2.*2.+0.38*0.38);
                        } else if (additionalNameOutput.CompareTo("EMC7")==0){
                            error = TMath::Sqrt(2.*2.+0.7*0.7);
                        } else if (additionalNameOutput.CompareTo("EMC1")==0 || additionalNameOutput.CompareTo("EG2")==0 ){
                            error = TMath::Sqrt(2.*2.+1.5*1.5);
                        } else if (additionalNameOutput.CompareTo("EG1")==0){
                            error = TMath::Sqrt(2.*2.+2.5*2.5);
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                }
            }
            // manual smoothing for dEdx pion line errors - variation 2
            if (nameCutVariationSC[i].CompareTo("dEdxPi")==0 ){
                if ( meson.CompareTo("Pi0") == 0){
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0.;
                        error   = 0.5+(-0.012)*ptBins[k]+(0.04)*ptBins[k]*ptBins[k]; // parametrisation
                        
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                } else if ( meson.CompareTo("Eta") == 0){
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 5;
                        error   = 1.5+(-0.012)*ptBins[k]+(0.04)*ptBins[k]*ptBins[k]; // parametrisation
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    } 
                } else {
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 5;
                        error   = 2.5+(-0.012)*ptBins[k]+(0.04)*ptBins[k]*ptBins[k]; // parametrisation
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    } 
                }
            }
            // manual smoothing for TPC cluster related errors - variation 3
            if (nameCutVariationSC[i].CompareTo("TPCCluster")==0 ){
                if (meson.CompareTo("Pi0") == 0){
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0.; // parametrisation
                        if (additionalNameOutput.CompareTo("")==0 || additionalNameOutput.CompareTo("INT7")==0){
                            error = 0.4;
                        } else if (additionalNameOutput.CompareTo("EMC1")==0 ||
                                   additionalNameOutput.CompareTo("EMC7")==0 ||
                                   additionalNameOutput.CompareTo("EG2")==0 ||
                                   additionalNameOutput.CompareTo("EG1")==0
                        ){
                            error = 0.8;
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = error*0.01;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = error*0.01;
                    }
                } else {
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 1.2;
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                }
            }
            // manual smoothing for single track momentum errors - variation 4
            if (nameCutVariationSC[i].CompareTo("SinglePt")==0){
                if(meson.CompareTo("Pi0") == 0){
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0;
                        if (additionalNameOutput.CompareTo("")==0 ||
                            additionalNameOutput.CompareTo("INT7")==0 ||
                            additionalNameOutput.CompareTo("EMC1")==0 ||
                            additionalNameOutput.CompareTo("EMC7")==0 ||
                            additionalNameOutput.CompareTo("EG2")==0  ||
                            additionalNameOutput.CompareTo("EG1")==0
                        ){
                            error   = 1.2+(-0.112)*ptBins[k]+(0.02)*ptBins[k]*ptBins[k]; // parametrisation
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                } else if(meson.CompareTo("Eta") == 0){
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0;
                        if (additionalNameOutput.CompareTo("")==0 ||
                            additionalNameOutput.CompareTo("INT7")==0 ||
                            additionalNameOutput.CompareTo("EMC1")==0 ||
                            additionalNameOutput.CompareTo("EMC7")==0 ||
                            additionalNameOutput.CompareTo("EG2")==0  ||
                            additionalNameOutput.CompareTo("EG1")==0
                        ){
                            error   = 3.2+(-0.112)*ptBins[k]+(0.02)*ptBins[k]*ptBins[k]; // parametrisation
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                } else {
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0;
                        if (additionalNameOutput.CompareTo("")==0 ||
                            additionalNameOutput.CompareTo("INT7")==0 ||
                            additionalNameOutput.CompareTo("EMC1")==0 ||
                            additionalNameOutput.CompareTo("EMC7")==0 ||
                            additionalNameOutput.CompareTo("EG2")==0  ||
                            additionalNameOutput.CompareTo("EG1")==0
                        ){
                            error   = 4.4+(-0.112)*ptBins[k]+(0.04)*ptBins[k]*ptBins[k]; // parametrisation
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                }
            }
            // manual smoothing for chi2/psi pair photon errors - variation 5
            if (nameCutVariationSC[i].CompareTo("Chi2")==0 ){
                if (meson.CompareTo("Pi0") == 0){
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0;
                        if (additionalNameOutput.CompareTo("")==0 ||
                            additionalNameOutput.CompareTo("INT7")==0 ||
                            additionalNameOutput.CompareTo("EMC1")==0 ||
                            additionalNameOutput.CompareTo("EMC7")==0 ||
                            additionalNameOutput.CompareTo("EG2")==0  ||
                            additionalNameOutput.CompareTo("EG1")==0
                        ){
                            error   = 1.1+(-0.112)*ptBins[k]+(0.035)*ptBins[k]*ptBins[k]; // parametrisation
//                             error   = 0.938+(-0.4177)*ptBins[k]+(0.1177)*ptBins[k]*ptBins[k]; // parametrisation
//                         } else if (additionalNameOutput.CompareTo("INT7")==0){
//                             error   = 0.868+(-0.112)*ptBins[k]+(0.0589)*ptBins[k]*ptBins[k]; // parametrisation
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = error*0.01;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = error*0.01;
                    }
                } else {
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 4;
                        if (additionalNameOutput.CompareTo("")==0 ||
                            additionalNameOutput.CompareTo("INT7")==0 ||
                            additionalNameOutput.CompareTo("EMC1")==0 ||
                            additionalNameOutput.CompareTo("EMC7")==0 ||
                            additionalNameOutput.CompareTo("EG2")==0  ||
                            additionalNameOutput.CompareTo("EG1")==0
                        ){
                            error   = 5.2+(-0.112)*ptBins[k]+(0.035)*ptBins[k]*ptBins[k]; // parametrisation
//                             error   = 0.938+(-0.4177)*ptBins[k]+(0.1177)*ptBins[k]*ptBins[k]; // parametrisation
//                         } else if (additionalNameOutput.CompareTo("INT7")==0){
//                             error   = 0.868+(-0.112)*ptBins[k]+(0.0589)*ptBins[k]*ptBins[k]; // parametrisation
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                }
            }
            // manual smoothing for Qt/alpha photon errors - variation 6
            if (nameCutVariationSC[i].CompareTo("Qt")==0 ){
                if (meson.CompareTo("Pi0") == 0){
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0;
                        if (additionalNameOutput.CompareTo("")==0 ||
                            additionalNameOutput.CompareTo("INT7")==0 ||
                            additionalNameOutput.CompareTo("EMC1")==0 ||
                            additionalNameOutput.CompareTo("EMC7")==0 ||
                            additionalNameOutput.CompareTo("EG2")==0  ||
                            additionalNameOutput.CompareTo("EG1")==0
                        ){
                            error = 0.6+(0.1)*ptBins[k]+(0.02)*ptBins[k]*ptBins[k]; // parametrisation
//                             error = 0.475+(-0.2895)*ptBins[k]+(0.0927)*ptBins[k]*ptBins[k]; // parametrisation
//                         } else if (additionalNameOutput.CompareTo("INT7")==0){
//                             error = 0.212+(0.051)*ptBins[k]+(0.062)*ptBins[k]*ptBins[k]; // parametrisation
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = error*0.01;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = error*0.01;
                    }
                } else if (meson.CompareTo("Eta") == 0){
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0.;
                        if (additionalNameOutput.CompareTo("")==0 ||
                            additionalNameOutput.CompareTo("INT7")==0 ||
                            additionalNameOutput.CompareTo("EMC1")==0 ||
                            additionalNameOutput.CompareTo("EMC7")==0 ||
                            additionalNameOutput.CompareTo("EG2")==0  ||
                            additionalNameOutput.CompareTo("EG1")==0
                        ){
                            error = 5.5+(0.1)*ptBins[k]+(0.02)*ptBins[k]*ptBins[k]; // parametrisation
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                } else {
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0.;
                        if (additionalNameOutput.CompareTo("")==0 ||
                            additionalNameOutput.CompareTo("INT7")==0 ||
                            additionalNameOutput.CompareTo("EMC1")==0 ||
                            additionalNameOutput.CompareTo("EMC7")==0 ||
                            additionalNameOutput.CompareTo("EG2")==0  ||
                            additionalNameOutput.CompareTo("EG1")==0
                        ){
                            error = 5.5+(0.1)*ptBins[k]+(0.02)*ptBins[k]*ptBins[k]; // parametrisation
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                } 
            }
            // manual smoothing for alpha meson errors - variation 7
            if (nameCutVariationSC[i].CompareTo("Alpha")==0 ){
                for (Int_t k = 0; k < nPtBins; k++){
                    if (ptBins[k] < 1.) continue;
                    Double_t error              = 0.;
                    if (additionalNameOutput.CompareTo("")==0 || additionalNameOutput.CompareTo("INT7")==0){
                        error   = 0.5; // parametrisation
                        if (meson.CompareTo("Eta") == 0) error  = 1.5;
                        if (meson.CompareTo("EtaToPi0") == 0) error  = 1.58;
                    } else if (additionalNameOutput.CompareTo("EMC1")==0 ||
                               additionalNameOutput.CompareTo("EMC7")==0 ||
                               additionalNameOutput.CompareTo("EG2")==0  ||
                               additionalNameOutput.CompareTo("EG1")==0
                    ){
                        error   = 0.; // parametrisation
                    }
                    
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = 0.01*error;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = 0.01*error;
                }
            }
            // manual smoothing for conversion acceptance cuts - variation 8
            if (nameCutVariationSC[i].CompareTo("ConvPhi")==0 ){
                if ( meson.CompareTo("Pi0") == 0){
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0.;
                        if (additionalNameOutput.CompareTo("")==0       ||
                            additionalNameOutput.CompareTo("INT7")==0   ||
                            additionalNameOutput.CompareTo("EMC7")==0
                        ){
                            error   = 0.5; // parametrisation
                        } else if (additionalNameOutput.CompareTo("EMC1")==0    ||
                                   additionalNameOutput.CompareTo("EG2")==0    ||
                                   additionalNameOutput.CompareTo("EG1")==0    
                        ){
                            error   = 1.8; // parametrisation
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                } else {
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0.;
                        if (additionalNameOutput.CompareTo("")==0       ||
                            additionalNameOutput.CompareTo("INT7")==0   
                        ){
                            error   = 1.; // parametrisation
                        } else if ( additionalNameOutput.CompareTo("EMC7")==0    ||
                                    additionalNameOutput.CompareTo("EMC1")==0    ||
                                    additionalNameOutput.CompareTo("EG2")==0     ||
                                    additionalNameOutput.CompareTo("EG1")==0    
                        ){
                            error   = 3.5; // parametrisation
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                }
            }
            // manual smoothing for minimum cluster energy errors - variation 9
            if (nameCutVariationSC[i].CompareTo("ClusterMinEnergy")==0  ){
                if(meson.CompareTo("EtaToPi0") == 0){
                    Double_t error              = TMath::Sqrt(2.2*2.2+1.13*1.13);                 
                    for (Int_t k = 0;k < nPtBins;k++){
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                } else if(meson.CompareTo("Pi0") == 0){
                    for (Int_t k = 0; k < nPtBins; k++){
                        errorsMean[i][k]        = 1.13;
                        errorsMeanErr[i][k]     = 0.011;
                        errorsMeanCorr[i][k]    = 1.13;
                        errorsMeanErrCorr[i][k] = 0.011;
                    }
                } else {
                    for (Int_t k = 0; k < nPtBins; k++){
                        errorsMean[i][k]        = 2.2;
                        errorsMeanErr[i][k]     = 0.02;
                        errorsMeanCorr[i][k]    = 2.2;
                        errorsMeanErrCorr[i][k] = 0.02;
                    }
                }
            }
            // manual smoothing for minimum number of cells in cluster errors - variation 10
            if (nameCutVariationSC[i].CompareTo("ClusterNCells")==0 ){ 
                for (Int_t k = 0; k < nPtBins; k++){
                    Double_t error              = 1.5;
                    if (meson.CompareTo("EtaToPi0") == 0 ) 
                        error                   = TMath::Sqrt(1.5*1.5+1.5*1.5);
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }
            // manual smoothing for energy calibration errors - variation 11
            if (nameCutVariationSC[i].CompareTo("ClusterNonLinearity")==0 ){ //&& meson.Contains("Pi0")
                for (Int_t k = 0; k < nPtBins; k++){
                    if ( ptBins[k] < 0.8 ) continue;
                    Double_t error              = 0;
                    if (additionalNameOutput.CompareTo("")==0){
                        error   = 2.+(0.01)*ptBins[k]+(0.01)*ptBins[k]*ptBins[k];
                    } else if (additionalNameOutput.CompareTo("INT7")==0){
                        error   = 2.+(0.01)*ptBins[k]+(0.01)*ptBins[k]*ptBins[k];
                    } else if (additionalNameOutput.CompareTo("EMC1")==0 || 
                               additionalNameOutput.CompareTo("EMC7")==0 ||
                               additionalNameOutput.CompareTo("EG2")==0  ||
                               additionalNameOutput.CompareTo("EG1")==0 
                    ){
                        error   = 2.+(0.01)*ptBins[k]+(0.01)*ptBins[k]*ptBins[k];
                    }
                    if (meson.CompareTo("Eta") == 0)
                        error   = 1.5* error;
                    if (meson.CompareTo("EtaToPi0") == 0) 
                        error   = 2*error;
    //              Double_t error = 0.79-0.21*ptBins[k]+0.17*ptBins[k]*ptBins[k]-0.00058*ptBins[k]*ptBins[k]*ptBins[k]*ptBins[k]; // parametrisation with No NonLinearity in 
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }
            // manual smoothing for cluster V0 matching errors - variation 12
            if (nameCutVariationSC[i].CompareTo("ClusterTrackMatching")==0 ){
                if ( meson.CompareTo("Pi0") == 0){
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0;
                        if (additionalNameOutput.CompareTo("")==0 ||
                            additionalNameOutput.CompareTo("INT7")==0 ||
                            additionalNameOutput.CompareTo("EMC1")==0 ||
                            additionalNameOutput.CompareTo("EMC7")==0 || 
                            additionalNameOutput.CompareTo("EG2")==0 ||
                            additionalNameOutput.CompareTo("EG1")==0
                        ){
                          error   = 0.5+(0.275)*ptBins[k]+(0.0045)*ptBins[k]*ptBins[k]; // parametrisation
//                             error   = 1.34+(-0.73)*ptBins[k]+(0.147)*ptBins[k]*ptBins[k]; // parametrisation
//                         } else if (additionalNameOutput.CompareTo("INT7")==0){
//                             error   = -0.0657+(0.275)*ptBins[k]+(0.0383)*ptBins[k]*ptBins[k]; // parametrisation
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                } else {
                    for (Int_t k = 0; k < nPtBins; k++){
                        Double_t error          = 0;
                        if (additionalNameOutput.CompareTo("")==0 ||
                            additionalNameOutput.CompareTo("INT7")==0 ||
                            additionalNameOutput.CompareTo("EMC1")==0 ||
                            additionalNameOutput.CompareTo("EMC7")==0 || 
                            additionalNameOutput.CompareTo("EG2")==0 ||
                            additionalNameOutput.CompareTo("EG1")==0
                        ){
                          error   = 4+(-0.2)*ptBins[k]+(0.04)*ptBins[k]*ptBins[k]+12/pow(2,ptBins[k]); // parametrisation
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                }
            }
           // manual smoothing for cluster shape errors - variation 13
            if (nameCutVariationSC[i].CompareTo("ClusterM02")==0 ){ //&& meson.Contains("Pi0")
                for (Int_t k = 0; k < nPtBins; k++){
                    Double_t error              =  1.2+(-0.01)*ptBins[k]+(0.015)*ptBins[k]*ptBins[k]; // parametrisation0.;
                    if ( meson.CompareTo("EtaToPi0") == 0 )
                        error = error*2;

                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }

            // manual smoothing for cell timing of EMC - variation 14
            if (nameCutVariationSC[i].CompareTo("CellTiming")==0 ){
                Double_t error                  = 0.95;
                if (meson.CompareTo("EtaToPi0") == 0)
                    error                       = 0;    // cancels fully for eta/pi0
                for (Int_t k = 0;k < nPtBins;k++){
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }
            
            // manual smoothing for Material infront of EMC - variation 15
            if (nameCutVariationSC[i].CompareTo("ClusterMaterialTRD")==0 ){
                Double_t error                  = 4.24; //(3% for TRD mat, 3% for TOF mat added in quadrature)
                if (meson.CompareTo("EtaToPi0") == 0)
                    error                       = 0;    // cancels fully for eta/pi0
                for (Int_t k = 0;k < nPtBins;k++){
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }
           
            // manual smoothing for Trigger normalization uncertainties - variation 16
            if (nameCutVariationSC[i].CompareTo("Trigger") == 0){
                for (Int_t k = 0; k < nPtBins; k++){
                    Double_t error              = 0.;
                    if (additionalNameOutput.CompareTo("")==0 ||
                        additionalNameOutput.CompareTo("INT7")==0 
                    ){
                        error   = 0.;
                    } else if (additionalNameOutput.CompareTo("EMC1")==0){
                        error   = TMath::Sqrt(4.755*4.755+2*2);
                    } else if (additionalNameOutput.CompareTo("EMC7")==0){
                        error   = TMath::Sqrt(4.208*4.208+2*2);
                    } else if (additionalNameOutput.CompareTo("EG2")==0){
                        error   = TMath::Sqrt(7.37*7.37+2*2);
                    } else if (additionalNameOutput.CompareTo("EG1")==0){
                        error   = TMath::Sqrt(10.59*10.59+2*2);
                    }
                    if (meson.CompareTo("EtaToPi0") == 0){
                        error   = 0.;
                    }
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
                
            }
           
            // manual smoothing for Efficiency uncertainties - variation 17
            if (nameCutVariationSC[i].CompareTo("Efficiency") == 0){
                for (Int_t k = 0; k < nPtBins; k++){
                    Double_t error              = 3.0;
                    if (meson.CompareTo("Pi0") == 0){
                        if (additionalNameOutput.CompareTo("")==0 ||
                            additionalNameOutput.CompareTo("INT7")==0 
                        ){
                            error   = 2.;
                        } else if (additionalNameOutput.CompareTo("EMC1")==0){
                            error   = 10000*pow(0.34,ptBins[k]+2.5)+2.;
                        } else if (additionalNameOutput.CompareTo("EMC7")==0){
                            error   = 5000*pow(0.5,ptBins[k]+6.2)+2.;
                        } else if (additionalNameOutput.CompareTo("EG2")==0){
                            error   = 10000*pow(0.34,ptBins[k]+2.5)+2.;
                        } else if (additionalNameOutput.CompareTo("EG1")==0){
                            error   = 200*pow(0.7,ptBins[k]+0.2)+2.;
                        }
                    } else if (meson.CompareTo("Eta") == 0){
                        if (additionalNameOutput.CompareTo("")==0 ||
                            additionalNameOutput.CompareTo("INT7")==0 
                        ){
                            error   = 5.;
                        } else if (additionalNameOutput.CompareTo("EMC1")==0){
                            error   = 10000*pow(0.34,ptBins[k]+2.5)+5.;
                        } else if (additionalNameOutput.CompareTo("EMC7")==0){
                            error   = 5000*pow(0.5,ptBins[k]+6.2)+5.;
                        } else if (additionalNameOutput.CompareTo("EG2")==0){
                            error   = 10000*pow(0.34,ptBins[k]+2.5)+5.;
                        } else if (additionalNameOutput.CompareTo("EG1")==0){
                            error   = 200*pow(0.7,ptBins[k]+0.2)+5.;
                        }
                    } else {
                        Double_t error1 = 0;
                        if (additionalNameOutput.CompareTo("")==0 ||
                            additionalNameOutput.CompareTo("INT7")==0 
                        ){
                            error1  = 2.;
                        } else if (additionalNameOutput.CompareTo("EMC1")==0){
                            error1  = 10000*pow(0.34,ptBins[k]+2.5)+2.;
                        } else if (additionalNameOutput.CompareTo("EMC7")==0){
                            error1  = 5000*pow(0.5,ptBins[k]+6.2)+2.;
                        } else if (additionalNameOutput.CompareTo("EG2")==0){
                            error1  = 10000*pow(0.34,ptBins[k]+2.5)+2.;
                        } else if (additionalNameOutput.CompareTo("EG1")==0){
                            error1  = 200*pow(0.7,ptBins[k]+0.2)+2.;
                        }
                        Double_t error2 = 0;
                        if (additionalNameOutput.CompareTo("")==0 ||
                            additionalNameOutput.CompareTo("INT7")==0 
                        ){
                            error2  = 5.;
                        } else if (additionalNameOutput.CompareTo("EMC1")==0){
                            error2  = 10000*pow(0.34,ptBins[k]+2.5)+5.;
                        } else if (additionalNameOutput.CompareTo("EMC7")==0){
                            error2  = 5000*pow(0.5,ptBins[k]+6.2)+5.;
                        } else if (additionalNameOutput.CompareTo("EG2")==0){
                            error2  = 10000*pow(0.34,ptBins[k]+2.5)+5.;
                        } else if (additionalNameOutput.CompareTo("EG1")==0){
                            error2  = 200*pow(0.7,ptBins[k]+0.2)+5.;
                        }
                        error = TMath::Sqrt(error1*error1+error2*error2);
                        cout << error << endl;                        
                    }    
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }

            
        } else {
            for (Int_t k = 0; k < nPtBins; k++){
                errorsMeanErr[i][k]         = 0.03;
                errorsMeanErrCorr[i][k]     = 0.03;
            }
        }
        // Quadratic sum of errors except material error infront of EMCal & inner material
        if (!nameCutVariationSC[i].CompareTo("ClusterMaterialTRD")==0){
            cout << "errors added quadratically" << endl;
            for (Int_t l = 0; l < nPtBins; l++){
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
    Double_t errorMaterial  = 4.50;
    if (meson.CompareTo("EtaToPi0") == 0) 
        errorMaterial       = 0.;
    
    // Calculate sqrt of summed errors for final errors, add material budget errors 
    for (Int_t l = 0; l < nPtBins; l++){
        errorsPosSummed[l]              = pow(errorsPosSummed[l],0.5);
        errorsMeanSummed[l]             = pow(errorsMeanSummed[l],0.5);
        errorsPosErrSummed[l]           = errorsPosSummed[l]*0.001;
        errorsMeanErrSummed[l]          = errorsMeanSummed[l]*0.001;
        errorsNegSummed[l]              = -pow(errorsNegSummed[l],0.5);
        errorsNegErrSummed[l]           = errorsNegSummed[l]*0.001;

        // add PCM & EMCal material errors
        errorsPosCorrMatSummed[l]       = pow(errorsPosCorrSummed[l]+ pow(errorMaterial ,2.) + pow(errorsPosCorr[14][l],2) ,0.5);
        errorsMeanCorrMatSummed[l]      = pow(errorsMeanCorrSummed[l]+ pow(errorMaterial ,2.)+ pow(errorsMeanCorr[14][l],2),0.5);
        errorsNegCorrMatSummed[l]       = -pow(errorsNegCorrSummed[l]+ pow(errorMaterial ,2.)+ pow(errorsNegCorr[14][l],2),0.5);
        
        errorsPosCorrSummed[l]          = pow(errorsPosCorrSummed[l],0.5);
        errorsMeanCorrSummed[l]         = pow(errorsMeanCorrSummed[l],0.5);
        errorsPosErrCorrSummed[l]       = errorsPosCorrSummed[l]*0.001;
        errorsMeanErrCorrSummed[l]      = errorsMeanCorrSummed[l]*0.001;
        errorsMeanErrCorrMatSummed[l]   = errorsMeanCorrMatSummed[l]*0.001;
        errorsNegCorrSummed[l]          = -pow(errorsNegCorrSummed[l],0.5);
        errorsNegErrCorrSummed[l]       = errorsNegCorrSummed[l]*0.001;
    }
    
    // Create material graph
    Double_t ptBinsMaterial [nPtBins];
    Double_t errorsMat      [nPtBins];
    for (Int_t l = 0; l < nPtBins; l++){
        errorsMat[l]            = errorMaterial;
        ptBinsMaterial[l]       = ptBins[l]+0.05;
        
    }
    TGraphErrors* graphMaterialError    = new TGraphErrors(nPtBins,ptBinsMaterial ,errorsMat ,ptBinsErr ,errorsMeanErrSummed );
    
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
    if (meson.CompareTo("Eta") == 0){
        minXLegend          = 0.23;
    }
    Double_t widthLegend    = 0.25;
    if (numberCutStudies> 9)  
        widthLegend         = 0.5;
    Double_t heightLegend   = 1.02 * 0.035 * (numberCutStudies+1);
    if (numberCutStudies> 9)  
        heightLegend        = 1.02 * 0.035 * (numberCutStudies/2+1);
    
    // ***************************************************************************************************
    // ****************************** Plot all mean erros separately *************************************
    // ***************************************************************************************************
    
    TCanvas* canvasSysErrMean = new TCanvas("canvasSysErrMean","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasSysErrMean, 0.08, 0.01, 0.015, 0.09);
 
        // create dummy histo
        TH2D *histo2DSysErrMean ;
        if ( meson.CompareTo("Pi0") == 0 ){
            histo2DSysErrMean = new TH2D("histo2DSysErrMean", "", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,30.);
        } else {
            histo2DSysErrMean = new TH2D("histo2DSysErrMean", "", 20,0.,ptBins[nPtBins-1]+1,1000.,0.,65.);
        }
        SetStyleHistoTH2ForGraphs( histo2DSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean systematic Err %", 0.03, 0.04, 0.03, 0.04,
                                1,0.9, 510, 510);
        histo2DSysErrMean->Draw();
        
        // create legend
        TLegend* legendMean = GetAndSetLegend2(minXLegend,maxYLegend-heightLegend,minXLegend+widthLegend,maxYLegend, 30);
        if (numberCutStudies> 9) legendMean->SetNColumns(2);
        
        for(Int_t i = 0; i< numberCutStudies ; i++){
            // ClusterMaterialTRD - 15, Trigger - 16
            if ( meson.CompareTo("EtaToPi0") == 0 && (i == 15 || i == 16 || i == 14) ){
                cout << "not drawing: " << nameCutVariation[i].Data() << endl;
                continue;
            }
            // Trigger - 16
            if ((additionalNameOutput.CompareTo("") == 0 || additionalNameOutput.CompareTo("INT7") == 0 ) && i == 16){
                cout << "not drawing: " << nameCutVariation[i].Data() << endl;
                continue;
            }
            // Alpha - 7
            if (!(additionalNameOutput.CompareTo("") == 0 || additionalNameOutput.CompareTo("INT7") == 0 ) ) {
                if ( i == 7 ){
                    cout << "not drawing: "<< nameCutVariation[i].Data() << endl;
                    continue;    
                }
            }

            
            DrawGammaSetMarkerTGraphErr(meanErrors[i], markerStyle[i], 1.,color[i],color[i]);
            meanErrors[i]->Draw("pE0,csame");
            legendMean->AddEntry(meanErrors[i],nameCutVariation[i].Data(),"p");
        }
        // PCM material error
        if ( meson.CompareTo("EtaToPi0") != 0 ){
            DrawGammaSetMarkerTGraphErr(graphMaterialError, 24, 1.,color[10],color[10]);
            graphMaterialError->Draw("pX0,csame");
            legendMean->AddEntry(graphMaterialError,"Inner Material","p");
        }
        legendMean->Draw();

        // plot labeling
        TLatex *labelMeson;
        if (meson.CompareTo("EtaToPi0") == 0){
            labelMeson= new TLatex(0.75,0.89,Form("#eta/#pi^{0} rec. #gamma_{conv}#gamma_{calo}"));
        } else if (meson.CompareTo("Pi0") == 0){
            labelMeson= new TLatex(0.75,0.89,Form("#pi^{0} #rightarrow #gamma_{conv}#gamma_{calo}"));
        } else {
            labelMeson= new TLatex(0.75,0.89,Form("#eta #rightarrow #gamma_{conv}#gamma_{calo}"));
        }
        SetStyleTLatex( labelMeson, 0.038,4);
        labelMeson->Draw();

        TLatex *labelCentrality = new TLatex(0.75,0.93,Form("%s",collisionSystem.Data() ));
        SetStyleTLatex( labelCentrality, 0.038,4);
        labelCentrality->Draw();

        TLatex *labelTrig;
        if (additionalNameOutput.CompareTo("")==0){
            labelTrig= new TLatex(0.75,0.84,Form("MB LHC11a"));
        } else if (additionalNameOutput.CompareTo("EMC1")==0){
            labelTrig= new TLatex(0.75,0.84,Form("EMC1 LHC11a"));
        } else if (additionalNameOutput.CompareTo("INT7")==0){
            labelTrig= new TLatex(0.75,0.84,Form("INT7 LHC13g"));
        } else if (additionalNameOutput.CompareTo("EMC7")==0){
            labelTrig= new TLatex(0.75,0.84,Form("EMC7 LHC13g"));
        } else if (additionalNameOutput.CompareTo("EG2")==0){
            labelTrig= new TLatex(0.75,0.84,Form("EG2 LHC13g"));
        } else if (additionalNameOutput.CompareTo("EG1")==0){
            labelTrig= new TLatex(0.75,0.84,Form("EG1 LHC13g"));
        }
        SetStyleTLatex( labelTrig, 0.038,4);
        labelTrig->Draw();
        
    canvasSysErrMean->Update();
    canvasSysErrMean->SaveAs(Form("SystematicErrorsCalculatedConvCalo/SysMean_%s_%s%s_%s.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
    
    delete canvasSysErrMean;
    
    
    // ***************************************************************************************************
    // ********************* Plot all mean erros separately after smoothing ******************************
    // ***************************************************************************************************    
    TCanvas* canvasNewSysErrMean = new TCanvas("canvasNewSysErrMean","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasNewSysErrMean, 0.08, 0.01, 0.015, 0.09);
    
        // create dummy histo
        TH2D *histo2DNewSysErrMean ;
        if ( meson.CompareTo("Pi0") == 0 ){
            histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,30.);
        } else { 
            histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,65.);
        }
        SetStyleHistoTH2ForGraphs( histo2DNewSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %", 0.03, 0.04, 0.03, 0.04,
                                1,0.9, 510, 510);
        histo2DNewSysErrMean->Draw();

        // create legend
        TLegend* legendMeanNew = GetAndSetLegend2(minXLegend,maxYLegend-heightLegend,minXLegend+widthLegend,maxYLegend, 30);
        legendMeanNew->SetMargin(0.1);
        if (numberCutStudies> 9) legendMeanNew->SetNColumns(2);

        for(Int_t i = 0; i< numberCutStudies ; i++){
            cout << i << "\t"<< additionalNameOutput.Data() << endl;
            // Trigger - 16
            if ((additionalNameOutput.CompareTo("") == 0 || additionalNameOutput.CompareTo("INT7") ==0 ) && i == 16){
                cout << "not drawing: " << nameCutVariation[i].Data() << endl;
                continue;
            }
            // Alpha - 7
            if (!(additionalNameOutput.CompareTo("") == 0 || additionalNameOutput.CompareTo("INT7") == 0) ) {
                if ( i == 7 ){
                    cout << "not drawing: "<< nameCutVariation[i].Data() << endl;
                    continue;    
                }
            }
            // ClusterMaterialTRD - 15, Trigger - 16
            if ( meson.CompareTo("EtaToPi0") == 0 && (i == 15 || i == 16 || i == 14) ){
                cout << "not drawing: " << nameCutVariation[i].Data() << endl;
                continue;
            }

            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], markerStyle[i], 1.,color[i],color[i]);
            meanErrorsCorr[i]->Draw("pX0,csame");
            legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");	
        }

        // PCM material
        if ( meson.CompareTo("EtaToPi0") != 0 ){
            DrawGammaSetMarkerTGraphErr(graphMaterialError, 24, 1.,color[10],color[10]);
            graphMaterialError->Draw("pX0,csame");
            legendMeanNew->AddEntry(graphMaterialError,"Inner Material","p");
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
    canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedConvCalo/SysMeanNewWithMean_%s_%s%s_%s.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
    
    // ***************************************************************************************************
    // ********************* Plot unsmoothed errors with fits ********************************************
    // ***************************************************************************************************    
    for (Int_t cut =0 ; cut < numberCutStudies; cut++ ){
        
        canvasNewSysErrMean->cd();
        histo2DNewSysErrMean->Draw();

            if (bsmooth[cut]) continue;
            cout <<endl << endl<<  "variation: " << cut << " \t"<< nameCutVariation[cut].Data() << endl;
            Double_t minPt = startPtSys;
//             if (additionalNameOutput.CompareTo("EMC1")==0)  minPt = 2.6;
            Double_t maxPt = ptBins[nPtBins-2]+1;
//             if (cut == 13) maxPt = 6;
//             if (cut == 12 || cut == 5) maxPt = 8;
//             if (cut == 6) maxPt = 8;

            TF1* pol0 = new TF1("pol0","[0]",minPt,maxPt); //
            TF1* pol1 = new TF1("pol1","[0]+[1]*x",minPt,maxPt); //
            TF1* pol2 = new TF1("pol2","[0]+[1]*x+[2]*x*x",minPt,maxPt); //
            TF1* pol4 = new TF1("pol4","[0]+[1]*x+[2]*x*x+[3]*x*x*x*x",minPt,maxPt); //
            pol4->SetParLimits(3,0,10);
            if (cut == 13) pol2->SetParLimits(2,0,0.1);
                
            meanErrorsCorr[cut]->Fit(pol4,"NRMEX0+","",minPt,maxPt);
            meanErrorsCorr[cut]->Fit(pol2,"NRMEX0+","",minPt,maxPt);
            meanErrorsCorr[cut]->Fit(pol1,"NRMEX0+","",minPt,maxPt);
            meanErrorsCorr[cut]->Fit(pol0,"NRMEX0+","",minPt,maxPt);
            pol4->SetLineColor(kRed+2);
            pol2->SetLineColor(kBlue+2);
            pol1->SetLineColor(kCyan+2);
            pol0->SetLineColor(kBlack);
            
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[cut], 20+cut, 1.,color[cut],color[cut]);
            meanErrorsCorr[cut]->Draw("p,csame");
            pol4->Draw("same");
            pol2->Draw("same");
            pol1->Draw("same");
            pol0->Draw("same");
        
        canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedConvCalo/SysMeanNewWithMeanSingle_%s_%s%s_%s_Variation%d.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),cut,suffix.Data()));
    }	
        
    // ***************************************************************************************************
    // ********************* Create output files with errors *********************************************
    // ***************************************************************************************************    
    const char *SysErrDatname = Form("SystematicErrorsCalculatedConvCalo/SystematicErrorPCMEMC_%s_%s%s_%s.dat",meson.Data(),energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data());
    fstream SysErrDat;
    cout << SysErrDatname << endl;
    SysErrDat.open(SysErrDatname, ios::out);
    for (Int_t l=0; l< nPtBins; l++){
        SysErrDat << ptBins[l] << "\t" <<errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
    }

    SysErrDat.close();

    const char *SysErrDatnameMean = Form("SystematicErrorsCalculatedConvCalo/SystematicErrorAveragedPCMEMC_%s_%s%s_%s.dat",meson.Data(),energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data());
    fstream SysErrDatAver;
    cout << SysErrDatnameMean << endl;
    SysErrDatAver.open(SysErrDatnameMean, ios::out);
    for (Int_t l=0; l< nPtBins; l++){
        SysErrDatAver << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
    }

    SysErrDatAver.close();
    
    const char *SysErrDatnameMeanSingleErr = Form("SystematicErrorsCalculatedConvCalo/SystematicErrorAveragedSinglePCMEMC_%s_%s%s_%s.dat",meson.Data(),energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data());
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

    
    // ***************************************************************************************************
    // ********************* Group errors according to topic *********************************************
    // ***************************************************************************************************    
    Double_t errorsMeanCorrPID[nPtBins];	
    Double_t errorsMeanCorrSignalExtraction[nPtBins];	
    Double_t errorsMeanCorrTrackReco[nPtBins];	
    Double_t errorsMeanCorrPhotonReco[nPtBins];	
    Double_t errorsMeanCorrClusterDescription[nPtBins];	
    
    for (Int_t l=0; l< nPtBins; l++){
        // "YieldExtraction"-0,"dEdxE"-1,"dEdxPi"-2, "TPCCluster"-3, "SinglePt"-4, "Chi2"-5, "Qt"-6, "Alpha"-7, "ConvPhi"-8, "ClusterMinEnergy"-9, "ClusterNCells"-10, 
        // "NonLinearity"-11, "ClusterTrackMatching" -12, "ClusterM02" -13, "CellTiming" -14,"ClusterMaterialTRD" -15
        // grouping:
        // Signal extraction: Yield extraction 0, Alpha 7 
        if (numberCutStudies>8){
            errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+errorsMeanCorr[7][l]*errorsMeanCorr[7][l]);	
        } else {
            errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]);
        }
        // PID: dEdxE 1, dEdxPi 2
        errorsMeanCorrPID[l] =TMath::Sqrt(errorsMeanCorr[1][l]*errorsMeanCorr[1][l]+ errorsMeanCorr[2][l]*errorsMeanCorr[2][l]);
        // photon reco: Chi2 5, Qt 6
        errorsMeanCorrPhotonReco[l] =TMath::Sqrt( errorsMeanCorr[5][l]* errorsMeanCorr[5][l]+errorsMeanCorr[6][l]*errorsMeanCorr[6][l]);
        // track reconstruction: TPCCluster 3, Single pt 4, ConvPhi 8
        errorsMeanCorrTrackReco[l] = TMath::Sqrt(errorsMeanCorr[3][l]*errorsMeanCorr[3][l]+errorsMeanCorr[4][l]*errorsMeanCorr[4][l]+errorsMeanCorr[8][l]*errorsMeanCorr[8][l]);				
        // cluster description in MC: ClusterMinEnergy 9, ClusterNCells 10, ClusterM02 13
        errorsMeanCorrClusterDescription[l] = TMath::Sqrt(errorsMeanCorr[9][l]*errorsMeanCorr[9][l]+errorsMeanCorr[10][l]*errorsMeanCorr[10][l]+errorsMeanCorr[13][l]*errorsMeanCorr[13][l]);
    }
        
        
    TGraphErrors* meanErrorsPID = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPID ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsPhotonReco = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPhotonReco ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsSignalExtraction = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSignalExtraction ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsTrackReco = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrTrackReco ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsClusterDescrip = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrClusterDescription ,ptBinsErr ,errorsMeanErrCorrSummed );
    
    cout << __LINE__ << endl;

    // ***************************************************************************************************
    // ********************* Plot grouped errors for better understanding ********************************
    // ***************************************************************************************************    
    Double_t minXLegend2 = 0.13;
    Double_t maxYLegend2 = 0.95;
    if (meson.CompareTo("Eta") == 0){
        minXLegend2 = 0.20;
    }
    Double_t widthLegend2 = 0.52;
    Double_t heightLegend2 = 0.25;
    
    TCanvas* canvasSummedErrMean = new TCanvas("canvasSummedErrMean","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasSummedErrMean, 0.08, 0.01, 0.015, 0.09);

        // create dummy histo
        TH2D *histo2DSummedErrMean ;
        if ( meson.CompareTo("Pi0") == 0 ){	
            histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,30.);
        } else {
            histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "", 20,0.,ptBins[nPtBins-1]+1,1000.,-0.5,65.);
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
        meanErrorsSignalExtraction->Draw("pX0,csame");
        // electron PID error
        DrawGammaSetMarkerTGraphErr(meanErrorsPID, 21, 1.,color[1],color[1]);
        meanErrorsPID->Draw("pX0,csame");
        // track reconstruction error
        DrawGammaSetMarkerTGraphErr(meanErrorsTrackReco, 22, 1.,color[2],color[2]);
        meanErrorsTrackReco->Draw("pX0,csame");
        // Photon reco error
        DrawGammaSetMarkerTGraphErr(meanErrorsPhotonReco, 23, 1.,color[3],color[3]);
        meanErrorsPhotonReco->Draw("pX0,csame");
        // Non linearity - 11
        DrawGammaSetMarkerTGraphErr(meanErrorsCorr[11], 25, 1.,color[5],color[5]);
        meanErrorsCorr[11]->Draw("pX0,csame");
        // Material infront of EMCAL -15
        if (meson.CompareTo("EtaToPi0") != 0){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[15], 21, 1.,color[7],color[7]);
            meanErrorsCorr[15]->Draw("pX0,csame");
        }
        // PCM Material budget 
        if (meson.CompareTo("EtaToPi0") != 0){
            DrawGammaSetMarkerTGraphErr(graphMaterialError, 24, 1.,color[4],color[4]);
            graphMaterialError->Draw("pX0,csame");
        }
        // Track matching V0 to EMCAL - 12
        DrawGammaSetMarkerTGraphErr(meanErrorsCorr[12], 20, 1.,color[18],color[18]);
        meanErrorsCorr[12]->Draw("pX0,csame");
        // Cluster description in MC
        DrawGammaSetMarkerTGraphErr(meanErrorsClusterDescrip, 22, 1.,color[8],color[8]);
        meanErrorsClusterDescrip->Draw("pX0,csame");
        // Efficiency - 17
        if (numberCutStudies>17){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[17], 20, 1.,color[12],color[12]);
            meanErrorsCorr[17]->Draw("pX0,csame");
        }
        
        if (meson.CompareTo("EtaToPi0") != 0){
            // Cell time - 14
            if (numberCutStudies>14){
                DrawGammaSetMarkerTGraphErr(meanErrorsCorr[14], 20, 1.,color[13],color[13]);
                meanErrorsCorr[14]->Draw("pX0,csame");
            }
            // trigger normalization -16
            if (numberCutStudies>16 && !(additionalNameOutput.CompareTo("") == 0 || additionalNameOutput.CompareTo("INT7") ==0 )){
                DrawGammaSetMarkerTGraphErr(meanErrorsCorr[16], 25, 1.,color[14],color[14]);
                meanErrorsCorr[16]->Draw("pX0,csame");
            }
            
        }
        
        legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"Signal Extraction","p");
        legendSummedMeanNew->AddEntry(meanErrorsPID,"Electron PID","p");
        legendSummedMeanNew->AddEntry(meanErrorsTrackReco,"Track Reconstruction","p");
        legendSummedMeanNew->AddEntry(meanErrorsPhotonReco,"Photon Reconstruction","p");
        if (meson.CompareTo("EtaToPi0") != 0){
            legendSummedMeanNew->AddEntry(graphMaterialError,"Inner Material","p");
        }
        legendSummedMeanNew->AddEntry(meanErrorsClusterDescrip,"Cluster Description","p");
        legendSummedMeanNew->AddEntry(meanErrorsCorr[11],"Cluster Energy Calibration","p");
        legendSummedMeanNew->AddEntry(meanErrorsCorr[12],"V0 tr. match. to cluster","p");
        if (meson.CompareTo("EtaToPi0") != 0){
            if (numberCutStudies>14) legendSummedMeanNew->AddEntry(meanErrorsCorr[14],"Cell time","p");
            legendSummedMeanNew->AddEntry(meanErrorsCorr[15],"Mat. infront of EMCal","p");
        }
        if (numberCutStudies>17) legendSummedMeanNew->AddEntry(meanErrorsCorr[17],"Efficiency","p");
        if (meson.CompareTo("EtaToPi0") != 0){
            if (numberCutStudies>16 && !(additionalNameOutput.CompareTo("") == 0 || additionalNameOutput.CompareTo("INT7") ==0 )) 
                legendSummedMeanNew->AddEntry(meanErrorsCorr[16],"Trigger normalization","p");
        }
        DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
        meanErrorsCorrSummedIncMat->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
        legendSummedMeanNew->Draw();
        
        labelMeson->Draw();
        labelCentrality->Draw();
        labelTrig->Draw();
    
    canvasSummedErrMean->Update();
    canvasSummedErrMean->SaveAs(Form("SystematicErrorsCalculatedConvCalo/SysErrorSummedVisu_%s_%s%s_%s.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
    
    delete canvasSummedErrMean;
// 
// 	
// 	
// 	const char *SysErrDatnameMeanPaper = Form("SystematicErrorsCalculatedConvCalo/SystematicErrorAveragedPaper_%s_%s%s_%s.dat",meson.Data(),energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data());
// 	fstream SysErrDatAverPaper;
// 	cout << SysErrDatnameMeanPaper << endl;
// 	SysErrDatAverPaper.open(SysErrDatnameMeanPaper, ios::out);
// 	SysErrDatAverPaper  << "#it{p}_{T}" << "\t Material \t Yield Extraction \t PID \t photon reco \t track recon \t summed" <<  endl;
// 	for (Int_t l=0; l< nPtBins; l++){
// 		SysErrDatAverPaper << ptBins[l] <<"\t" << errorMaterial*2 << "\t" << errorsMeanCorrSignalExtraction[l] << "\t" << errorsMeanCorrPID[l]<< "\t" << errorsMeanCorrPhotonReco[l]<< "\t" <<errorsMeanCorrTrackReco[l] <<"\t" << errorsMeanCorrMatSummed[l]<< endl;
// 	}
// 
// 	SysErrDatAverPaper.close();
	
}