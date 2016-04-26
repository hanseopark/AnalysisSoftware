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

void FinalyseSystematicErrorsCalo_pPbV2(     const char* nameDataFileErrors  = "", 
                                            TString energy                  = "", 
                                            TString meson                   = "", 
                                            Int_t numberOfPtBins            = 1, 
                                            Int_t numberCutStudies          = 1, 
                                            Int_t offSetBeginning           = 0, 
                                            TString additionalName          = "pp", 
                                            TString additionalNameOutput    = "", 
                                            TString suffix                  = "eps",
                                            Int_t mode                      = 4
                                        ){
    
 // ***************************************************************************************************
    // ****************************** General style settings *********************************************
    // ***************************************************************************************************
    StyleSettingsThesis();
    SetPlotStyle();
    
    // ***************************************************************************************************
    // ****************************** Create output directory ********************************************
    // ***************************************************************************************************    
    gSystem->Exec("mkdir -p SystematicErrorsCalculatedCalo");
    
    // ***************************************************************************************************
    // ***************************** labeling and color settings *****************************************
    // ***************************************************************************************************
    TString date                            = ReturnDateString();
    TString dateForOutput                   = ReturnDateStringForOutput();
    TString collisionSystem                 = ReturnFullCollisionsSystem(energy);
    
    Color_t color[20]                       = { kBlue, kRed+1, kOrange+7, kPink+8, kGreen+2, 
                                                kYellow+2, kOrange+2, kBlue+2, kCyan-2, kViolet+1, 
                                                kAzure+1, 432, kPink+4, kOrange, 407, 
                                                416, 830, 404, kPink-6, 1};
    Color_t markerStyle[20]                 = { 24, 21, 22, 23, 20, 
                                                25, 26, 27, 28, 29,
                                                30, 31, 32, 33, 24, 
                                                21, 22, 23, 20, 25};
    
    // ***************************************************************************************************
    // ******************************* general variable definition  **************************************
    // ***************************************************************************************************
    Int_t   numberOfEntriesPos              = 0;
    Int_t   numberOfEntriesNeg              = 0;
    const Int_t nPtBins                     = numberOfPtBins;
    const Int_t nCuts                       = numberCutStudies;
    Double_t* ptBins;
    Double_t* ptBinsErr;
    TString nameCutVariation[12];
    TString nameCutVariationSC[12];
    
    TString nameCutVariation2760GeV[12]     = { "Yield extraction", "#theta_{#gamma#gamma}", "min E_{cluster}", "min # cells", "clu. energy calibration", 
                                                "tr. match. to cl.", "M_{02}", "Mat. infront of EMCal", "clu. energy scale", "Trigger normalization",
                                                "Efficiency", "Yield extraction #pi^{0}" };
    TString nameCutVariationSC2760GeV[12]   = { "YieldExtraction", "OpeningAngle", "ClusterMinEnergy", "ClusterNCells", "ClusterNonLinearity", 
                                                "ClusterTrackMatchingCalo", "ClusterM02","ClusterMaterialTRD", "ClusterEnergyScale" ,"Trigger", 
                                                "Efficiency", "YieldExtraction"};
    TString nameCutVariation5023GeV[12]     = { "Yield extraction", "#theta_{#gamma#gamma}", "min E_{cluster}", "min # cells", "clu. energy calibration", 
                                                "tr. match. to cl.", "M_{02}", "Mat. infront of EMCal", "Rapidity", "clu. energy scale",
						"Efficiency", "Yield extraction #pi^{0}"};
    TString nameCutVariationSC5023GeV[12]   = { "YieldExtraction", "OpeningAngle", "ClusterMinEnergy", "ClusterNCells", "ClusterNonLinearity", 
                                                "ClusterTrackMatchingCalo", "ClusterM02","ClusterMaterialTRD", "Rapidity", "ClusterEnergyScale",
						"Efficiency", "YieldExtraction"};
    if (meson.CompareTo("EtaToPi0") == 0){
        nameCutVariation2760GeV[0]          = "Yield extraction #eta";
	nameCutVariation5023GeV[0]          = "Yield extraction #eta";
    }
    if (energy.CompareTo("pPb_5.023TeV") == 0) {
        for (Int_t i = 0;i < numberCutStudies;i++){
            nameCutVariation[i]             = nameCutVariation5023GeV[i];
            nameCutVariationSC[i]           = nameCutVariationSC5023GeV[i];
        }
    }
    
    // ***************************************************************************************************
    // ******************************** Booleans for smoothing *******************************************
    // ***************************************************************************************************
    Bool_t bsmooth[12]                      = { 0, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0, 
                                                0, 0 };
    Bool_t bsmoothMBPi0[12]                 = { 0, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1, 0 };
    Bool_t bsmoothMBEta[12]                 = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1, 0 };
    Bool_t bsmoothMBEtaToPi0[12]            = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1, 0 };
                          
    for (Int_t i = 0; i < numberCutStudies; i++){
        if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Pi0")==0){
            bsmooth[i]                      = bsmoothMBPi0[i];
        } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Eta")==0){
            bsmooth[i]                      = bsmoothMBEta[i];
        } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("EtaToPi0")==0){
            bsmooth[i]                      = bsmoothMBEtaToPi0[i];
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
    TFile* fileErrorInput= new TFile(nameDataFileErrors);
    
    for (Int_t i = 0;i < nCuts;i++){
        
        // read data
        TGraphAsymmErrors* graphPosErrors;
        TGraphAsymmErrors* graphNegErrors;
        if (i == 0 || i == 9 || i==10){// || i == 8 || i == 9 || i == 10special treatment for Yield extraction error and calculated erros
            TString nameGraphPos    = "";
            TString nameGraphNeg    = "";
            if ( meson.CompareTo("EtaToPi0") != 0 ){
                nameGraphPos        = Form("%s_SystErrorRelPos_YieldExtraction_%s",meson.Data(),additionalName.Data() );
                nameGraphNeg        = Form("%s_SystErrorRelNeg_YieldExtraction_%s",meson.Data(),additionalName.Data() );
            } else {
                nameGraphPos        = Form("Eta_SystErrorRelPos_YieldExtraction_%s",additionalName.Data() );
                nameGraphNeg        = Form("Eta_SystErrorRelNeg_YieldExtraction_%s",additionalName.Data() );                
            }    
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        } else if ( i == 11) { // special treatment for eta to pi0 ratio    
            TString nameGraphPos    = Form("Pi0EtaBinning_SystErrorRelPos_YieldExtraction_%s",additionalName.Data() );
            TString nameGraphNeg    = Form("Pi0EtaBinning_SystErrorRelNeg_YieldExtraction_%s",additionalName.Data() );                
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
            
        } else {// read graphs from input file
            TString nameGraphPos    = Form("%s_SystErrorRelPos_%s%s",meson.Data(),nameCutVariationSC[i].Data(),additionalNameOutput.Data()  );
            TString nameGraphNeg    = Form("%s_SystErrorRelNeg_%s%s",meson.Data(),nameCutVariationSC[i].Data(),additionalNameOutput.Data()  );
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            graphPosErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        }
        
        // take out offsets
        for (Int_t j = 0;j < offSetBeginning;j++){
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
                cout << "Yield extraction smoothing" << endl;
                if (meson.CompareTo("Eta") == 0 || (meson.CompareTo("EtaToPi0") == 0 && i == 0) ){
                    for (Int_t k = 0;k < nPtBins;k++){
                        Double_t error;
                        if(ptBins[k]<7){
                          error = 5.0;
                        }else{
                          error = 7.17279 - 0.691497*ptBins[k] + 0.0640519*ptBins[k]*ptBins[k];
                        }
                        errorsMean[i][k]            = error;
                        errorsMeanErr[i][k]         = error*0.01;
                        errorsMeanCorr[i][k]        = error;
                        errorsMeanErrCorr[i][k]     = error*0.01;
                    }
                }
                
            }    
            // manual smoothing for Yield extraction errors - variation 1
            if  (nameCutVariationSC[i].CompareTo("OpeningAngle") == 0){
                cout << "Opening Angle smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                  Double_t error;
                  if(ptBins[k]<4){
                    error = 0.45;
                  }else{
                    error              = 0.753156-0.237041*ptBins[k]+0.0327857*ptBins[k]*ptBins[k];
                  }
                  if (meson.CompareTo("Eta") == 0 ) error = error*1.25;
                  if (meson.CompareTo("EtaToPi0") == 0 ) error = TMath::Sqrt(1+1.25*1.25)*error;
                  errorsMean[i][k]            = error;
                  errorsMeanErr[i][k]         = error*0.01;
                  errorsMeanCorr[i][k]        = error;
                  errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }    

            // manual smoothing for minimum cluster energy errors - variation 2
            if (nameCutVariationSC[i].CompareTo("ClusterMinEnergy")==0  ){
                cout << "Cluster minimum energy smoothing" << endl;
                Double_t error = 2;
                if (meson.CompareTo("Eta") == 0 ) error = error*2;
                if (meson.CompareTo("EtaToPi0") == 0 ) error = error*2;
                  for (Int_t k = 0;k < nPtBins;k++){
                      errorsMean[i][k]        = error;
                      errorsMeanErr[i][k]     = 0.01*error;
                      errorsMeanCorr[i][k]    = error;
                      errorsMeanErrCorr[i][k] = 0.01*error;
                  }
            }
            // manual smoothing for minimum number of cells in cluster errors - variation 3
            if (nameCutVariationSC[i].CompareTo("ClusterNCells")==0 ){
                cout << "Cluster NCells smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 1.5;
                    if (meson.CompareTo("EtaToPi0") == 0 ) error = TMath::Sqrt(2.)*1.5;
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }   
            }
            // manual smoothing for energy calibration errors - variation 4
            if (nameCutVariationSC[i].CompareTo("ClusterNonLinearity")==0 ){//&& meson.Contains("Pi0")
                cout << "Cluster non linearity smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 0.0104652+0.0860654*ptBins[k]+0.00316317*ptBins[k]*ptBins[k];
                    if (meson.CompareTo("EtaToPi0") == 0 ) error *= 2;
                    error = TMath::Sqrt(error*error+0.95*0.95);//adding 0.95% error for timing cut
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }
            // manual smoothing for cluster matching errors - variation 5
            if (nameCutVariationSC[i].CompareTo("ClusterTrackMatchingCalo")==0 ){
                cout << "Cluster track matching smoothing" << endl;
                    for (Int_t k = 0;k < nPtBins;k++){
                      Double_t error = 2.5+(+0.09)*ptBins[k];
                      if (meson.CompareTo("Eta") == 0){
                        error   = error*2;
                      }
                      if( meson.CompareTo("EtaToPi0") == 0 ) error = TMath::Sqrt(5)*error;
                      errorsMean[i][k]        = error;
                      errorsMeanErr[i][k]     = 0.01*error;
                      errorsMeanCorr[i][k]    = error;
                      errorsMeanErrCorr[i][k] = 0.01*error;
                    }   
            }
           // manual smoothing for cluster shape errors - variation 6
            if (nameCutVariationSC[i].CompareTo("ClusterM02")==0 ){//&& meson.Contains("Pi0")
                cout << "Cluster M02 smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                  Double_t error = 0.266401+0.191655*ptBins[k]+0.02*ptBins[k]*ptBins[k];
                  if( meson.CompareTo("EtaToPi0") == 0 ) error = TMath::Sqrt(2.)*error;
                  errorsMean[i][k]            = error;
                  errorsMeanErr[i][k]         = error*0.01;
                  errorsMeanCorr[i][k]        = error;
                  errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }
            // manual smoothing for Material infront of EMC - variation 7
            if (nameCutVariationSC[i].CompareTo("ClusterMaterialTRD")==0 ){
                cout << "Material smoothing" << endl;
                Double_t error                  = 5.83; //(3% for TRD mat, 5% for TOF mat added in quadrature)
                if (meson.CompareTo("EtaToPi0") == 0)
                    error                       = 0;    // cancels fully for eta/pi0
                for (Int_t k = 0;k < nPtBins;k++){
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }   
            }
           
            // manual smoothing for energy scale errors (derived from mass difference MC & Data) - variation 8
            if (nameCutVariationSC[i].CompareTo("ClusterEnergyScale")==0 ){//&& meson.Contains("Pi0")
                cout << "Cluster non linearity smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 0.5;
                    if (meson.CompareTo("Eta") == 0)
                        error   = 2* error;
                    if (meson.CompareTo("EtaToPi0") == 0)
                        error   = 2* error;
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }

           
            // manual smoothing for Trigger normalization uncertainties - variation 9
            if (nameCutVariationSC[i].CompareTo("Trigger") == 0){
                cout << "Trigger smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 0.;
                    if (additionalNameOutput.CompareTo("")==0 ){
                        error   = 0.;
                    }    
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
                
            }    
            // manual smoothing for Efficiency uncertainties - variation 10
            if (nameCutVariationSC[i].CompareTo("Efficiency") == 0){
                cout << "Efficiency smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 4.0;
                    if (meson.CompareTo("Eta") == 0) error   = 6.0;
                    if (meson.CompareTo("EtaToPi0")== 0) error = TMath::Sqrt(4*4+6*6);
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
                
            }
            
            // manual smoothing for Rapidity
            if (nameCutVariationSC[i].CompareTo("Rapidity")==0 ){
                cout << "Rapidity" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 1.5;
                    if (meson.CompareTo("Eta") == 0 || meson.CompareTo("EtaToPi0") == 0){
                      error *= 2;
                    }
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
    
    for (Int_t i = 0;i < nCuts;i++){
      for (Int_t k = 0;k < nPtBins;k++){
	cout << " error " << errorsMean[i][k] << endl;
      }
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
        errorsPosCorrMatSummed[l]       = pow(errorsPosCorrSummed[l]+ pow(errorMaterial ,2.) + pow(errorsPosCorr[7][l],2) ,0.5);
        errorsMeanCorrMatSummed[l]      = pow(errorsMeanCorrSummed[l]+ pow(errorMaterial ,2.)+ pow(errorsMeanCorr[7][l],2),0.5);
	cout << " quad sum " << errorsMeanCorrMatSummed[l] << endl;
        errorsNegCorrMatSummed[l]       = -pow(errorsNegCorrSummed[l]+ pow(errorMaterial ,2.)+ pow(errorsNegCorr[7][l],2),0.5);
        
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
    if (meson.CompareTo("Eta") == 0){
        minXLegend          = 0.23;
    }
    Double_t widthLegend    = 0.25;
    if (numberCutStudies> 7)  
        widthLegend         = 0.5;
    Double_t heightLegend   = 1.15 * 0.035 * (numberCutStudies+3);
    if (numberCutStudies> 7)  
        heightLegend        = 1.15 * 0.035 * (numberCutStudies/2+2);
    
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
            if ((additionalNameOutput.CompareTo("") == 0 || additionalNameOutput.CompareTo("INT7") ==0 ) && i == 12){
                cout << "not drawing: " << nameCutVariation[i].Data() << endl;
                continue;
            }    
//             if ( meson.CompareTo("Eta") == 0 && i == 1){
//                 cout << "not drawing: " << nameCutVariation[i].Data() << endl;
//                 continue;
//             }    
//             if ( meson.CompareTo("EtaToPi0") == 0 && i == 7){
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
        if (meson.CompareTo("EtaToPi0") == 0){
            labelMeson= new TLatex(0.75,0.89,Form("#eta/#pi^{0} rec. #gamma_{calo}"));
        } else if (meson.Contains("Pi0")){
            labelMeson= new TLatex(0.75,0.89,Form("#pi^{0} #rightarrow #gamma_{calo}#gamma_{calo}"));
        } else {
            labelMeson= new TLatex(0.75,0.89,Form("#eta #rightarrow #gamma_{calo}#gamma_{calo}"));
        }
        SetStyleTLatex( labelMeson, 0.038,4);
        labelMeson->Draw();

        TLatex *labelCentrality = new TLatex(0.75,0.93,Form("%s",collisionSystem.Data() ));
        SetStyleTLatex( labelCentrality, 0.038,4);
        labelCentrality->Draw();

        TLatex *labelTrig;
        if (additionalNameOutput.CompareTo("")==0){
            labelTrig= new TLatex(0.75,0.84,Form("MB LHC13bc"));
        } 
        SetStyleTLatex( labelTrig, 0.038,4);
        labelTrig->Draw();
        
    canvasSysErrMean->Update();
    canvasSysErrMean->SaveAs(Form("SystematicErrorsCalculatedCalo/SysMean_%s_%s%s_%s.%s",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
    
    delete canvasSysErrMean;
    
    
    // ***************************************************************************************************
    // ********************* Plot all mean erros separately after smoothing ******************************
    // ***************************************************************************************************    
    TCanvas* canvasNewSysErrMean = new TCanvas("canvasNewSysErrMean","",200,10,1350,900);// gives the page size
    DrawGammaCanvasSettings( canvasNewSysErrMean, 0.08, 0.01, 0.015, 0.09);
    
        // create dummy histo
        TH2D *histo2DNewSysErrMean ;
        if (meson.Contains("Pi0")){
            histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,40.);
        } else {
            histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,40.);
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
            if ((additionalNameOutput.CompareTo("") == 0 || additionalNameOutput.CompareTo("INT7") ==0 ) && i == 12){
                cout << "not drawing: " << nameCutVariation[i].Data() << endl;
                continue;
            }    
/*            if ( meson.CompareTo("Eta") == 0 && i == 1){
                cout << "not drawing: " << nameCutVariation[i].Data() << endl;
                continue;
            }   */ 
            if ( meson.CompareTo("EtaToPi0") == 0 && i == 7){
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
    canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedCalo/SysMeanNewWithMean_%s_%s%s_%s.%s",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
    
    // ***************************************************************************************************
    // ********************* Plot unsmoothed errors with fits ********************************************
    // ***************************************************************************************************    
    for (Int_t cut =0 ;cut < numberCutStudies;cut++ ){
        
        canvasNewSysErrMean->cd();
        histo2DNewSysErrMean->Draw();

            if (bsmooth[cut]) continue;
            cout <<endl << endl<<  "variation: " << cut << " \t"<< nameCutVariation[cut].Data() << endl;
            Double_t minPt = 0.6;
//             if (additionalNameOutput.CompareTo("EMC1")==0)  minPt = 2.6;
            Double_t maxPt = ptBins[nPtBins-2]+3;
	    if (meson.CompareTo("Eta") == 0)  minPt = 2.6;
	    if (meson.CompareTo("Eta") == 0)  maxPt = ptBins[nPtBins-2]+4;
//             if (cut == 5) maxPt = 8;
//             if (cut == 12 || cut == 5) maxPt = 8;
//             if (cut == 6) maxPt = 8;

            TF1* pol0 = new TF1("pol0","[0]",minPt,maxPt);//
            TF1* pol1 = new TF1("pol1","[0]+[1]*x",minPt,maxPt);//
            TF1* pol2 = new TF1("pol2","[0]+[1]*x+[2]*x*x",minPt,maxPt);//
            TF1* pol4 = new TF1("pol4","[0]+[1]*x+[2]*x*x+[3]*x*x*x*x",minPt,maxPt);//
            pol4->SetParLimits(3,0,10);
                
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
	    
	    TLegend* leg = new TLegend(0.6,0.6,0.8,0.9);
	    leg->AddEntry("pol0","pol0","l");
	    leg->AddEntry("pol1","pol1","l");
	    leg->AddEntry("pol2","pol2","l");
	    leg->AddEntry("pol4","pol4","l");
	    leg->Draw();
        
        canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedCalo/SysMeanNewWithMeanSingle_%s_%s%s_%s_Variation%d.%s",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),cut,suffix.Data()));
    }   
        
    // ***************************************************************************************************
    // ********************* Create output files with errors *********************************************
    // ***************************************************************************************************    
    const char *SysErrDatname = Form("SystematicErrorsCalculatedCalo/SystematicErrorEMCEMC_%s_%s%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
    fstream SysErrDat;
    cout << SysErrDatname << endl;
    SysErrDat.open(SysErrDatname, ios::out);
    for (Int_t l=0;l< nPtBins;l++){
        SysErrDat << ptBins[l] << "\t" <<errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
    }

    SysErrDat.close();

    const char *SysErrDatnameMean = Form("SystematicErrorsCalculatedCalo/SystematicErrorAveragedEMCEMC_%s_%s%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
    fstream SysErrDatAver;
    cout << SysErrDatnameMean << endl;
    SysErrDatAver.open(SysErrDatnameMean, ios::out);
    for (Int_t l=0;l< nPtBins;l++){
        SysErrDatAver << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
    }

    SysErrDatAver.close();
    
    // ***************************************************************************************************
    // ********************* Group errors according to topic *********************************************
    // ***************************************************************************************************
    Double_t errorsMeanCorrSignalExtraction[nPtBins];
    Double_t errorsMeanCorrClusterEnergy[nPtBins];
    Double_t errorsMeanCorrClusterDescription[nPtBins];
    
    for (Int_t l=0;l< nPtBins;l++){
//      "YieldExtraction"-0,"OpeningAngle"-1, "ClusterMinEnergy"-2, "ClusterNCells"-3, "NonLinearity"-4, "ClusterTrackMatchingCalo" -5, "ClusterM02" -6, "ClusterMaterialTRD" -7, "ClusterEnergyScale" -8
        // grouping:
        // Signal extraction: Yield extraction 0, Open-Angle 1 
        errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+errorsMeanCorr[1][l]*errorsMeanCorr[1][l]);
        // Signal extraction: NonLinearity 5 , Energy Scale 8 
        //errorsMeanCorrClusterEnergy[l]      = TMath::Sqrt(errorsMeanCorr[4][l]*errorsMeanCorr[4][l]+errorsMeanCorr[8][l]*errorsMeanCorr[8][l]);
        errorsMeanCorrClusterEnergy[l]      = TMath::Sqrt(errorsMeanCorr[4][l]*errorsMeanCorr[4][l]);//+errorsMeanCorr[8][l]*errorsMeanCorr[8][l]);
        // cluster description in MC: ClusterMinEnergy 2, ClusterNCells 3, ClusterM02 6
        errorsMeanCorrClusterDescription[l] = TMath::Sqrt(errorsMeanCorr[2][l]*errorsMeanCorr[2][l]+errorsMeanCorr[3][l]*errorsMeanCorr[3][l]+errorsMeanCorr[6][l]*errorsMeanCorr[6][l]);
    }
        
        
    TGraphErrors* meanErrorsSignalExtraction    = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSignalExtraction ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsClusterDescrip      = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrClusterDescription ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsClusterEnergy       = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrClusterEnergy ,ptBinsErr ,errorsMeanErrCorrSummed );
    
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
            histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,30.);
        } else {
            histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,65.);
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
        // Cluster description in MC
        DrawGammaSetMarkerTGraphErr(meanErrorsClusterDescrip, 22, 1.,color[1],color[1]);
        // Track matching to EMCAL
        DrawGammaSetMarkerTGraphErr(meanErrorsCorr[5], 25, 1.,color[2],color[2]);
        // Material infront of EMCAL
        DrawGammaSetMarkerTGraphErr(meanErrorsCorr[7], 21, 1.,kMagenta+2,kMagenta+2);
        // Cluster energy description
        DrawGammaSetMarkerTGraphErr(meanErrorsClusterEnergy, 20, 1.,color[4],color[4]);
        meanErrorsClusterEnergy->Draw("pX0,csame");
        meanErrorsCorr[7]->Draw("pX0,csame");
        meanErrorsCorr[5]->Draw("pX0,csame");
        meanErrorsClusterDescrip->Draw("pX0,csame");
        if (numberCutStudies>10){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[10], 23, 1.,color[8],color[8]);
            meanErrorsCorr[10]->Draw("pX0,csame");
        }
        if (numberCutStudies>9 && !(additionalNameOutput.CompareTo("") == 0 || additionalNameOutput.CompareTo("INT7") ==0 )){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[9], 25, 1.,color[6],color[6]);
            meanErrorsCorr[9]->Draw("pX0,csame");
        }

        legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"Signal Extraction","p");
        legendSummedMeanNew->AddEntry(meanErrorsClusterDescrip,"Cluster Description","p");
        legendSummedMeanNew->AddEntry(meanErrorsClusterEnergy,"Cluster Energy Description","p");
        legendSummedMeanNew->AddEntry(meanErrorsCorr[5],"track match. to cluster","p");
        legendSummedMeanNew->AddEntry(meanErrorsCorr[7],"Mat. infront of EMCal","p");
        if (numberCutStudies>10) legendSummedMeanNew->AddEntry(meanErrorsCorr[10],"Efficiency","p");
        if (numberCutStudies>9 && !(additionalNameOutput.CompareTo("") == 0 || additionalNameOutput.CompareTo("INT7") ==0 )) 
            legendSummedMeanNew->AddEntry(meanErrorsCorr[9],"Trigger normalization","p");
        DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
        meanErrorsCorrSummedIncMat->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
        legendSummedMeanNew->Draw();
        
        labelMeson->Draw();
        labelCentrality->Draw();
        labelTrig->Draw();
    
    canvasSummedErrMean->Update();
    canvasSummedErrMean->SaveAs(Form("SystematicErrorsCalculatedCalo/SysErrorSummedVisu_%s_%s%s_%s.%s",meson.Data(), energy.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
    
    delete canvasSummedErrMean;
// 
//  
//  
//  const char *SysErrDatnameMeanPaper = Form("SystematicErrorsCalculatedConvCalo/SystematicErrorAveragedPaper_%s_%s%s_%s.dat",meson.Data(),energy.Data(),additionalNameOutput.Data(),dateForOutput.Data());
//  fstream SysErrDatAverPaper;
//  cout << SysErrDatnameMeanPaper << endl;
//  SysErrDatAverPaper.open(SysErrDatnameMeanPaper, ios::out);
//  SysErrDatAverPaper  << "#it{p}_{T}" << "\t Material \t Yield Extraction \t PID \t photon reco \t track recon \t summed" <<  endl;
//  for (Int_t l=0;l< nPtBins;l++){
//      SysErrDatAverPaper << ptBins[l] <<"\t" << errorMaterial*2 << "\t" << errorsMeanCorrSignalExtraction[l] << "\t" << errorsMeanCorrPID[l]<< "\t" << errorsMeanCorrPhotonReco[l]<< "\t" <<errorsMeanCorrTrackReco[l] <<"\t" << errorsMeanCorrMatSummed[l]<< endl;
//  }
// 
//  SysErrDatAverPaper.close();
        
}