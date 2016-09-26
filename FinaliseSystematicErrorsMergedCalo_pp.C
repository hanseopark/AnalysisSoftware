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

void FinaliseSystematicErrorsMergedCalo_pp( TString nameDataFileErrors      = "", 
                                            TString energy                  = "", 
                                            TString meson                   = "", 
                                            Int_t numberOfPtBins            = 1, 
                                            Int_t numberCutStudies          = 1, 
                                            Double_t startPtSys             = 0, 
                                            TString additionalNameOutput    = "", 
                                            TString suffix                  = "eps"
                                        ){
    
    if (numberCutStudies > 11) {
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
    TString labelMesonReco                  = textMeson;
    if (additionalNameOutput.Contains("NLM1") )
        labelMesonReco                      = Form("%s rec. w/ EMCal m. cl., %d lm", textMeson.Data(), 1);
    else if (additionalNameOutput.Contains("NLM2") )
        labelMesonReco                      = Form("%s rec. w/ EMCal m. cl., %d lm", textMeson.Data(), 2);
    
    // ***************************************************************************************************
    // ******************************* general variable definition  **************************************
    // ***************************************************************************************************
    Int_t   numberOfEntriesPos              = 0;    
    Int_t   numberOfEntriesNeg              = 0;
    const Int_t nPtBins                     = numberOfPtBins;
    const Int_t nCuts                       = numberCutStudies;
    Double_t* ptBins;
    Double_t* ptBinsErr;
    TString nameCutVariation[11];
    TString nameCutVariationSC[11];
    
    TString nameCutVariationSC2760GeV[11]   = { "ClusterTrackMatchingCalo", "ClusterM02", "ClusterNonLinearity",  "CellMinE", "CellTiming",  
                                                "ClusterMaterialTRD", "MesonResolution", "ClusterEnergyScale" , "Trigger", "Efficiency", "Secondary"};
    
    Color_t color[20];
    Color_t markerStyle[20];
    for (Int_t k = 0; k < 11; k++ ){
        color[k]        = GetColorSystematics( nameCutVariationSC2760GeV[k], 10 ); 
        markerStyle[k]  = GetMarkerStyleSystematics( nameCutVariationSC2760GeV[k], 10 );     
    }

    for (Int_t i = 0; i < numberCutStudies; i++){
        nameCutVariation[i]     = GetSystematicsName(nameCutVariationSC2760GeV[i]);
        nameCutVariationSC[i]   = nameCutVariationSC2760GeV[i];
    }
    
    // ***************************************************************************************************
    // ******************************** Booleans for smoothing *******************************************
    // ***************************************************************************************************
    Bool_t bsmooth[11]                      = { 0, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0,
                                                0
                                              };
    Bool_t bsmoothMBPi0[11]                 = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1
                                              };
    Bool_t bsmoothINT7Pi0[11]               = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1
                                              };
    Bool_t bsmoothEMC1Pi0[11]               = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1
                                              };
    Bool_t bsmoothEMC7Pi0[11]               = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1
                                              };
    Bool_t bsmoothEG2Pi0[11]                = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1
                                              };
    Bool_t bsmoothEG1Pi0[11]                = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1
                                              };
                          
    for (Int_t i = 0; i < numberCutStudies; i++){
        if ( (additionalNameOutput.Contains("INT1") || additionalNameOutput.Contains("MB") ) && meson.CompareTo("Pi0")==0){
            bsmooth[i]                      = bsmoothMBPi0[i];
        } else if (additionalNameOutput.Contains("INT7") && meson.CompareTo("Pi0")==0){
            bsmooth[i]                      = bsmoothINT7Pi0[i];
        } else if (additionalNameOutput.Contains("EMC1") && meson.CompareTo("Pi0")==0){
            bsmooth[i]                      = bsmoothEMC1Pi0[i];
        } else if (additionalNameOutput.Contains("EMC7") && meson.CompareTo("Pi0")==0){
            bsmooth[i]                      = bsmoothEMC7Pi0[i];
        } else if (additionalNameOutput.Contains("EG2") && meson.CompareTo("Pi0")==0){
            bsmooth[i]                      = bsmoothEG2Pi0[i];
        } else if (additionalNameOutput.Contains("EG1") && meson.CompareTo("Pi0")==0){
            bsmooth[i]                      = bsmoothEG1Pi0[i];
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
            nameCutVariationSC[i].CompareTo("Secondary") == 0
        ){
          nameGraphPos    = Form("%s_SystErrorRelPos_%s%s",meson.Data(),nameCutVariationSC[0].Data(),additionalNameOutput.Data()  );
          nameGraphNeg    = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),nameCutVariationSC[0].Data(),additionalNameOutput.Data()  );
        }  

        graphPosErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
        graphNegErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        
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

            // manual smoothing for cluster matching errors - variation 0
            if (nameCutVariationSC[i].CompareTo("ClusterTrackMatchingCalo")==0 ){
                cout << "Cluster track matching smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error          = 0;
                    if (additionalNameOutput.Contains("INT1") ||
                        additionalNameOutput.Contains("INT7") ||
                        additionalNameOutput.Contains("EMC1") ||
                        additionalNameOutput.Contains("EMC7") || 
                        additionalNameOutput.Contains("EG2") ||
                        additionalNameOutput.Contains("EG1")
                    ){
                        if (additionalNameOutput.Contains("NLM1"))
                            error   = 5+(-0.01)*ptBins[k]+(0.002)*ptBins[k]*ptBins[k];
//                             error   = 2+(-0.01)*ptBins[k]+(0.002)*ptBins[k]*ptBins[k];
//                            error   = 2.2+(-0.01)*ptBins[k]+(0.005)*ptBins[k]*ptBins[k]; // V1 clusterizer
                        else if (additionalNameOutput.Contains("NLM2"))
                            error   = 1.4+(-0.01)*ptBins[k]+(0.005)*ptBins[k]*ptBins[k];
                    }
//                     if (ptBins[k] > 12){
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
//                     }    
                }   
            }
              
           // manual smoothing for cluster shape errors - variation 1
            if (nameCutVariationSC[i].CompareTo("ClusterM02")==0 ){//&& meson.Contains("Pi0")
                cout << "Cluster M02 smoothing" << endl;
                if (additionalNameOutput.Contains("NLM1")){
                    for (Int_t k = 0;k < nPtBins;k++){
//                         if ( ptBins[k] > 12 ){
                            Double_t error              = 3.2+(-0.02)*ptBins[k]+(0.0027)*ptBins[k]*ptBins[k];                            
//                             Double_t error              = 1.6+(-0.02)*ptBins[k]+(0.0027)*ptBins[k]*ptBins[k];                            
//                             Double_t error              = 2.3+(-0.02)*ptBins[k]+(0.007)*ptBins[k]*ptBins[k]; // V1 clusterizer
                            errorsMean[i][k]            = error;
                            errorsMeanErr[i][k]         = error*0.01;
                            errorsMeanCorr[i][k]        = error;
                            errorsMeanErrCorr[i][k]     = error*0.01;
//                         }
                    }     
                } else if (additionalNameOutput.Contains("NLM2")) {
                    for (Int_t k = 0;k < nPtBins;k++){
                        if ( ptBins[k] > 8 ){                           
                            Double_t error              = 0.6+(-0.01)*ptBins[k]+(0.005)*ptBins[k]*ptBins[k];
                            errorsMean[i][k]            = error;
                            errorsMeanErr[i][k]         = error*0.01;
                            errorsMeanCorr[i][k]        = error;
                            errorsMeanErrCorr[i][k]     = error*0.01;
                        }    
                    }
                }     
            }
              
            // manual smoothing for energy calibration errors - variation 2
            if (nameCutVariationSC[i].CompareTo("ClusterNonLinearity")==0 ){//&& meson.Contains("Pi0")
                cout << "Cluster non linearity smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 0;
                    if ( (  additionalNameOutput.Contains("INT1")     || 
                            additionalNameOutput.Contains("EMC1") ) &&  additionalNameOutput.Contains("NLM2")
                    ){
                        
                        error   = 2.1+(0.02)*ptBins[k]+(0.005)*ptBins[k]*ptBins[k];
                    } else if ( (  additionalNameOutput.Contains("INT1")     || 
                            additionalNameOutput.Contains("EMC1") ) &&  additionalNameOutput.Contains("NLM1")
                    ){
                        error   = 3.0+(0.01)*ptBins[k]+(0.001)*ptBins[k]*ptBins[k];
//                         error   = 4.1+(0.02)*ptBins[k]+(0.001)*ptBins[k]*ptBins[k];
                        
                    } else if ( additionalNameOutput.Contains("INT7") ||
                                additionalNameOutput.Contains("EMC7") ||
                                additionalNameOutput.Contains("EG2")  ||
                                additionalNameOutput.Contains("EG1")  
                    ){
                       error   = 3.5+(0.01)*ptBins[k]+(0.0018)*ptBins[k]*ptBins[k]; 
//                        error   = 2.5+(0.01)*ptBins[k]+(0.004)*ptBins[k]*ptBins[k];  // V1 clusterizer
                    }  
                    
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }
            
            // manual smoothing for meson mass errors - variation 3
            if (nameCutVariationSC[i].CompareTo("MesonMass")==0 ){
                cout << "Mass error smoothing" << endl;
                if (additionalNameOutput.Contains("NLM2") ){
                    for (Int_t k = 0;k < nPtBins;k++){
                        if ( ptBins[k] > 7 ){
                            Double_t error              = 1.8;
                            if ( ptBins[k] > 27 ){
                                error                   = 11.2;
                            } else if ( ptBins[k] > 15 ){
                                error                   = 3.3;
                            }    
                            errorsMean[i][k]            = error;
                            errorsMeanErr[i][k]         = error*0.01;
                            errorsMeanCorr[i][k]        = error;
                            errorsMeanErrCorr[i][k]     = error*0.01;
                        }
                    }  
                } else if (additionalNameOutput.Contains("NLM1") ){
                    for (Int_t k = 0;k < nPtBins;k++){
                        if ( ptBins[k] > 10 ){
                            Double_t error              = 5;
//                             error                       = 2.+(0.001)*ptBins[k]+(0.004)*ptBins[k]*ptBins[k]; // V1 clusterizer
                            errorsMean[i][k]            = error;
                            errorsMeanErr[i][k]         = error*0.01;
                            errorsMeanCorr[i][k]        = error;
                            errorsMeanErrCorr[i][k]     = error*0.01;
                        }
                    }
                } 
            }    

            // manual smoothing for meson resolution errors 
            if (nameCutVariationSC[i].CompareTo("MesonResolution")==0 ){
                cout << "Meson resolution error smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 700*pow(0.7,ptBins[k]+0.2)+5.;
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }    
            // manual smoothing for Secondary errors 
            if (nameCutVariationSC[i].CompareTo("Secondary")==0 ){
                cout << "Meson resolution error smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 1.8;
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }    
            
            // manual smoothing for cell aggregation - variation 4
            if (nameCutVariationSC[i].CompareTo("CellMinE")==0 ){
                cout << "cell Emin error smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 0;
                    error                       = 2.;
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;

                }
            }    

            
            // manual smoothing for meson alpha errors - variation 4
            if (nameCutVariationSC[i].CompareTo("MesonAlpha")==0 ){
                cout << "Alpha error smoothing" << endl;
                if (additionalNameOutput.Contains("NLM2") ){
                    for (Int_t k = 0;k < nPtBins;k++){
                        if ( ptBins[k] > 7 ){
                            Double_t error              = 0;
                            error                       = 0.5+(0.001)*ptBins[k]+(0.002)*ptBins[k]*ptBins[k];
                            errorsMean[i][k]            = error;
                            errorsMeanErr[i][k]         = error*0.01;
                            errorsMeanCorr[i][k]        = error;
                            errorsMeanErrCorr[i][k]     = error*0.01;
                        }
                    }  
                } else if (additionalNameOutput.Contains("NLM1") ){
                    for (Int_t k = 0;k < nPtBins;k++){
//                         if ( ptBins[k] > 11 ){
                            Double_t error              = 0;
                            error                       = 20000*pow(0.5,ptBins[k])+0.75;
//                             error                       = 1.+(0.001)*ptBins[k]+(0.0015)*ptBins[k]*ptBins[k]; // V1 clusterizer
                            errorsMean[i][k]            = error;
                            errorsMeanErr[i][k]         = error*0.01;
                            errorsMeanCorr[i][k]        = error;
                            errorsMeanErrCorr[i][k]     = error*0.01;
//                         }
                    }
                } 
            }    

            // manual smoothing for cell time uncertainties - variation 5
            if (nameCutVariationSC[i].CompareTo("CellTiming") == 0){
                cout << "Cell time smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 0.;
                    if (additionalNameOutput.Contains("INT1") ||
                        additionalNameOutput.Contains("INT7") 
                    ){
                        error   = 1.5+(0.001)*ptBins[k]+(0.0015)*ptBins[k]*ptBins[k];
                    } else if (additionalNameOutput.Contains("EMC1")){
                        error   = 1.5+(0.001)*ptBins[k]+(0.0015)*ptBins[k]*ptBins[k];
                    } else if (additionalNameOutput.Contains("EMC7")){
                        error   = 1.5+(0.001)*ptBins[k]+(0.0015)*ptBins[k]*ptBins[k];
                    } else if (additionalNameOutput.Contains("EG2")){
                        error   = 1.5+(0.001)*ptBins[k]+(0.0015)*ptBins[k]*ptBins[k];
                    } else if (additionalNameOutput.Contains("EG1")){
                        error   = 1.5+(0.001)*ptBins[k]+(0.0015)*ptBins[k]*ptBins[k];
                    }    
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }                
            }    

            // manual smoothing for Material infront of EMC - variation 6
            if (nameCutVariationSC[i].CompareTo("ClusterMaterialTRD")==0 ){
                cout << "Material smoothing" << endl;
                Double_t error                  = 4.24; //(3% for TRD mat, 3% for TOF mat added in quadrature)
                for (Int_t k = 0;k < nPtBins;k++){
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }   
            }

            // manual smoothing for energy scale errors (derived from mass difference MC & Data) - variation 7
            if (nameCutVariationSC[i].CompareTo("ClusterEnergyScale")==0 ){//&& meson.Contains("Pi0")
                cout << "Cluster non linearity smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 0;
                    if (additionalNameOutput.Contains("INT1") || 
                        additionalNameOutput.Contains("INT7") ||
                        additionalNameOutput.Contains("EMC1") || 
                        additionalNameOutput.Contains("EMC7") ||
                        additionalNameOutput.Contains("EG2")  ||
                        additionalNameOutput.Contains("EG1") 
                    ){
                        error   = 0.5+7/pow(1.9,ptBins[k]);   
//                             1.5+(0.01)*ptBins[k]+(0.01)*ptBins[k]*ptBins[k]+1;
                    }    
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }

            // manual smoothing for Trigger normalization uncertainties - variation 8
            if (nameCutVariationSC[i].CompareTo("Trigger") == 0){
                cout << "Trigger smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 0.;
                    if (additionalNameOutput.Contains("INT1") ||
                        additionalNameOutput.Contains("INT7") 
                    ){
                        error   = 0.;
                    } else if (additionalNameOutput.Contains("EMC1")){
                        error   = TMath::Sqrt(10.5*10.5+2*2); // V2 clusterizer 10.5 % trigger uncertainty, 2% pileup
//                         error   = TMath::Sqrt(8.59*8.59+2*2); //V1 Clusterizer
                    } else if (additionalNameOutput.Contains("EMC7")){
                        error   = TMath::Sqrt(3.21*3.21+2*2); // V2 clusterizer 3.21% trigger uncertainty, 2% pileup
//                         error   = TMath::Sqrt(3.936*3.936+2*2);// V1 Clusterizer
                    } else if (additionalNameOutput.Contains("EG2")){
                        error   = TMath::Sqrt(6.17*6.17+2*2); // V2 clusterizer 6.17% trigger uncertainty, 2% pileup (trigger: EG2/EMC7 5.27%)
//                         error   = TMath::Sqrt(5.08*5.08+2*2); // V1 Clusterizer
                    } else if (additionalNameOutput.Contains("EG1")){
                        error   = TMath::Sqrt(8.6*8.6+2*2); // V2 clusterizer 8.6, 2% pileup (trigger EG1/EG2: 5.99%)
//                         error   = TMath::Sqrt(13.0*13.0+2*2); // V1 clusterizer
                    }    
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }                
            }    

            // manual smoothing for Efficiency uncertainties - variation 9
            if (nameCutVariationSC[i].CompareTo("Efficiency") == 0){
                cout << "Efficiency smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 5.0;
                    Double_t errorPi0           = 5.0;
                    if (additionalNameOutput.Contains("INT1") ||
                        additionalNameOutput.Contains("INT7") 
                    ){
                        errorPi0    = 5.;
                    } else if (additionalNameOutput.Contains("EMC1")){
                        errorPi0    = 10000*pow(0.34,ptBins[k]+2.5)+5.;
                    } else if (additionalNameOutput.Contains("EMC7")){
                        errorPi0    = 5000*pow(0.5,ptBins[k]+6.2)+5.;
                    } else if (additionalNameOutput.Contains("EG2")){
                        errorPi0    = 10000*pow(0.34,ptBins[k]+2.5)+5.;
                    } else if (additionalNameOutput.Contains("EG1")){
                        errorPi0    = 200*pow(0.7,ptBins[k]+0.2)+5.;
                    }    
                    error   = errorPi0; 
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }    

            
        } else {
            for (Int_t k = 0;k < nPtBins;k++){
                errorsMeanErr[i][k]         = 0.03;
                errorsMeanErrCorr[i][k]     = 0.03;
            }   
        }    
        // Quadratic sum of errors except material error infront of EMCal & inner material
        if (!nameCutVariationSC[i].CompareTo("ClusterMaterialTRD")==0){
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
            histo2DSysErrMean = new TH2D("histo2DSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,0.,30.);
        } else {
            histo2DSysErrMean = new TH2D("histo2DSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,0.,65.);
        }
        SetStyleHistoTH2ForGraphs( histo2DSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean systematic Err %", 0.03, 0.04, 0.03, 0.04,
                                1,0.9, 510, 510);
        histo2DSysErrMean->Draw();
        
        // create legend
        TLegend* legendMean = GetAndSetLegend2(minXLegend,maxYLegend-heightLegend,minXLegend+widthLegend,maxYLegend, 30);
        if (numberCutStudies> 7) legendMean->SetNColumns(2);
        
        for(Int_t i = 0;i< numberCutStudies ;i++){
//             if ((additionalNameOutput.CompareTo("INT1") == 0 || additionalNameOutput.CompareTo("INT7") ==0 ) && i == 10){
//                 cout << "not drawing: " << nameCutVariation[i].Data() << endl;
//                 continue;
//             }    
            
            DrawGammaSetMarkerTGraphErr(meanErrors[i], markerStyle[i], 1.,color[i],color[i]);
            meanErrors[i]->Draw("pE0,csame");           
            legendMean->AddEntry(meanErrors[i],nameCutVariation[i].Data(),"p");
        
        }
        legendMean->Draw();

        // plot labeling
        TLatex *labelMeson;
        labelMeson= new TLatex(0.68,0.885,labelMesonReco);
        SetStyleTLatex( labelMeson, 0.038,4);
        labelMeson->Draw();

        TLatex *labelCentrality = new TLatex(0.68,0.93,Form("%s",collisionSystem.Data() ));
        SetStyleTLatex( labelCentrality, 0.038,4);
        labelCentrality->Draw();

        TLatex *labelTrig;
        if (additionalNameOutput.Contains("INT1") || additionalNameOutput.Contains("MB") ){
            labelTrig= new TLatex(0.68,0.84,Form("MB LHC11a"));
        } else if (additionalNameOutput.Contains("EMC1")){
            labelTrig= new TLatex(0.68,0.84,Form("EMC1 LHC11a"));
        } else if (additionalNameOutput.Contains("INT7")){
            labelTrig= new TLatex(0.68,0.84,Form("INT7 LHC13g"));
        } else if (additionalNameOutput.Contains("EMC7")){
            labelTrig= new TLatex(0.68,0.84,Form("EMC7 LHC13g"));
        } else if (additionalNameOutput.Contains("EG2")){
            labelTrig= new TLatex(0.68,0.84,Form("EG2 LHC13g"));
        } else if (additionalNameOutput.Contains("EG1")){
            labelTrig= new TLatex(0.68,0.84,Form("EG1 LHC13g"));
        }
        SetStyleTLatex( labelTrig, 0.038,4);
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
            histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,30.);
        } else {
            histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,65.);
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
//             if ((additionalNameOutput.CompareTo("INT1") == 0 || additionalNameOutput.CompareTo("INT7") ==0 ) && i == 10){
//                 cout << "not drawing: " << nameCutVariation[i].Data() << endl;
//                 continue;
//             }    
            
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
            
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[cut], 20+cut, 1.,color[cut],color[cut]);
            meanErrorsCorr[cut]->Draw("p,csame");
            pol4->Draw("same");
            pol2->Draw("same");
            pol1->Draw("same");
            pol0->Draw("same");
            bla->Draw("same");
            
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