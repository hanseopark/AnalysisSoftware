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

void FinaliseSystematicErrorsCalo_pp8TeV(   const char* nameDataFileErrors  = "",
                                            TString energy                  = "", 
                                            TString meson                   = "", 
                                            Int_t numberOfPtBins            = 1, 
                                            Int_t numberCutStudies          = 1, 
                                            Double_t startPtSys             = 0, 
                                            TString additionalName          = "pp", 
                                            TString additionalNameOutput    = "", 
                                            TString suffix                  = "eps"
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
    TString energyForOutput                 = energy;
    energyForOutput.ReplaceAll(".","_");
    
    Color_t color[20]                       = { kBlue, kRed+1, kOrange+7, kPink+8, kGreen+2, 
                                                kYellow+2, kOrange+2, kBlue+2, kCyan-2, kViolet+1, 
                                                kAzure+1, kAzure+4, kPink+4, kOrange, 407, 
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
    TString nameCutVariation[16];
    TString nameCutVariationSC[16];
    
    TString nameCutVariation8TeV[16]     = { "yield extraction", "#theta_{#gamma#gamma}", "min E_{cluster}", "min # cells", "clu. energy calibration",
                                                "tr. match. to cl.", "M_{02}", "mat. infront of EMCal", "clu. energy scale", "cell time",
                                                "trigger normalization", "efficiency", "clu. timing","clusterization energy", "periods", "yield extraction #pi^{0}" };
    TString nameCutVariationSC8TeV[16]   = { "YieldExtraction", "OpeningAngle", "ClusterMinEnergy", "ClusterNCells", "ClusterNonLinearity",
                                                "ClusterTrackMatchingCalo", "ClusterM02","ClusterMaterialTRD", "ClusterEnergyScale" , "CellTiming", 
                                                "Trigger", "Efficiency", "ClusterTime", "ClusterizationEnergy", "Periods", "YieldExtraction"};
    if (meson.CompareTo("Pi0EtaBinning") == 0){
        nameCutVariation8TeV[0]          = "Yield extraction #eta";
    }
    
    if (energy.CompareTo("8TeV") == 0) {
        for (Int_t i = 0;i < numberCutStudies;i++){
            nameCutVariation[i]             = nameCutVariation8TeV[i];
            nameCutVariationSC[i]           = nameCutVariationSC8TeV[i];
        }
    }
    
    // ***************************************************************************************************
    // ******************************** Booleans for smoothing *******************************************
    // ***************************************************************************************************
    Bool_t bsmooth[16]                      = { 0, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1, 0};
    Bool_t bsmoothMBPi0[16]                 = { 0, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1, 0};
    Bool_t bsmoothMBEta[16]                 = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1, 0};
    Bool_t bsmoothMBPi0EtaBinning[16]       = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1, 1};
    Bool_t bsmoothEMC7Pi0[16]               = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1, 0};
    Bool_t bsmoothEMC7Eta[16]               = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1, 0};
    Bool_t bsmoothEMC7Pi0EtaBinning[16]     = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1, 1};
    Bool_t bsmoothEGAPi0[16]                = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1, 0};
    Bool_t bsmoothEGAEta[16]                = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1, 0};
    Bool_t bsmoothEGAPi0EtaBinning[16]      = { 1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1,
                                                1, 1, 1, 1, 1, 1};
                          
    for (Int_t i = 0; i < numberCutStudies; i++){
        if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Pi0")==0){
            bsmooth[i]                      = bsmoothMBPi0[i];
        } else if (additionalNameOutput.CompareTo("EMC7") == 0 && meson.CompareTo("Pi0")==0){
            bsmooth[i]                      = bsmoothEMC7Pi0[i];
        } else if (additionalNameOutput.CompareTo("EGA") == 0 && meson.CompareTo("Pi0")==0){
            bsmooth[i]                      = bsmoothEGAPi0[i];
        } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Eta")==0){
            bsmooth[i]                      = bsmoothMBEta[i];
        } else if (additionalNameOutput.CompareTo("EMC7") == 0 && meson.CompareTo("Eta")==0){
            bsmooth[i]                      = bsmoothEMC7Eta[i];
        } else if (additionalNameOutput.CompareTo("EGA") == 0 && meson.CompareTo("Eta")==0){
            bsmooth[i]                      = bsmoothEGAEta[i];
        } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Pi0EtaBinning")==0){
            bsmooth[i]                      = bsmoothMBPi0EtaBinning[i];
        } else if (additionalNameOutput.CompareTo("EMC7") == 0 && meson.CompareTo("Pi0EtaBinning")==0){
            bsmooth[i]                      = bsmoothEMC7Pi0EtaBinning[i];
        } else if (additionalNameOutput.CompareTo("EGA") == 0 && meson.CompareTo("Pi0EtaBinning")==0){
            bsmooth[i]                      = bsmoothEGAPi0EtaBinning[i];
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
    TString sFilePi0EtaBinning = "~/data/work/pcgGit/AnalysisSoftware/8TeV_CaloSystematics/CutStudies/8TeV/Eta_data_SystematicErrorCuts.root";
    if(additionalNameOutput.CompareTo("EMC7")==0) sFilePi0EtaBinning="~/data/work/pcgGit/AnalysisSoftware/8TeV_CaloSystematicsEMC7/CutStudies/8TeV/Eta_data_SystematicErrorCuts.root";
    else if(additionalNameOutput.CompareTo("EGA")==0) sFilePi0EtaBinning="~/data/work/pcgGit/AnalysisSoftware/8TeV_CaloSystematicsEGA/CutStudies/8TeV/Eta_data_SystematicErrorCuts.root";
    TFile* filePi0EtaBinning=new TFile(sFilePi0EtaBinning);
    
    for (Int_t i = 0;i < nCuts;i++){
        
        // read data
        TGraphAsymmErrors* graphPosErrors;
        TGraphAsymmErrors* graphNegErrors;
        if (i == 0 || i == 8 || i == 10 || i == 11 || (i == 13 && (additionalNameOutput.CompareTo("EMC7") == 0 || additionalNameOutput.CompareTo("EGA") == 0) ) || (i == 14 && (additionalNameOutput.CompareTo("EMC7") == 0 || additionalNameOutput.CompareTo("EGA") == 0) ) || (i == 9 && (additionalNameOutput.CompareTo("EMC7") == 0 || additionalNameOutput.CompareTo("EGA") == 0 || meson.CompareTo("Pi0EtaBinning") == 0))){// special treatment for Yield extraction error and calculated erros
            TString nameGraphPos    = "";
            TString nameGraphNeg    = "";
            if ( meson.CompareTo("Pi0EtaBinning") != 0 ){
                nameGraphPos        = Form("%s_SystErrorRelPos_YieldExtraction_%s",meson.Data(),additionalName.Data() );
                nameGraphNeg        = Form("%s_SystErrorRelNeg_YieldExtraction_%s",meson.Data(),additionalName.Data() );
            } else {
                nameGraphPos        = Form("Eta_SystErrorRelPos_YieldExtraction_%s",additionalName.Data() );
                nameGraphNeg        = Form("Eta_SystErrorRelNeg_YieldExtraction_%s",additionalName.Data() );                
            }    
            cout << "Cutstudies " << i<< "\t" <<nameGraphPos.Data() << "\t" << nameGraphNeg.Data()<<  endl;
            if ( meson.CompareTo("Pi0EtaBinning") == 0 ){
              graphPosErrors          = (TGraphAsymmErrors*)filePi0EtaBinning->Get(nameGraphPos.Data());
              graphNegErrors          = (TGraphAsymmErrors*)filePi0EtaBinning->Get(nameGraphNeg.Data());
            }else{
              graphPosErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
              graphNegErrors          = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
            }
        } else if ( i == 15) { // special treatment for eta to pi0 ratio
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
            // manual smoothing for Yield extraction errors - variation 0
            if  (nameCutVariationSC[i].CompareTo("YieldExtraction") == 0){
                cout << "Yield extraction smoothing" << endl;
                if ( (meson.CompareTo("Eta") == 0 || (meson.CompareTo("Pi0EtaBinning") == 0) && i == 0) ){
                    Double_t error = 7.5;
                    for (Int_t k = 0;k < nPtBins;k++){
                        error   = 7.5;
                        if(ptBins[k]>=8.) error = 12.65 + (-1.2)*ptBins[k] + 0.07*ptBins[k]*ptBins[k];
                        errorsMean[i][k]            = error;
                        errorsMeanErr[i][k]         = error*0.01;
                        errorsMeanCorr[i][k]        = error;
                        errorsMeanErrCorr[i][k]     = error*0.01;
                    }
                    
                    if (additionalNameOutput.CompareTo("EMC7")==0){
                        for (Int_t k = 0;k < nPtBins;k++){
                            error   = 6.14/pow(1.35,ptBins[k])+4.3;
                            if(ptBins[k]>=15.) error += (0.0125)*(pow((ptBins[k]-15.),2));

                            errorsMean[i][k]            = error;
                            errorsMeanErr[i][k]         = error*0.01;
                            errorsMeanCorr[i][k]        = error;
                            errorsMeanErrCorr[i][k]     = error*0.01;
                        }
                    } else if(additionalNameOutput.CompareTo("EGA")==0){
                      for (Int_t k = 0;k < nPtBins;k++){
                          error   = 8.;
                          if(ptBins[k]>=20.) error += (0.0125)*(pow((ptBins[k]-20.),2));

                          errorsMean[i][k]            = error;
                          errorsMeanErr[i][k]         = error*0.01;
                          errorsMeanCorr[i][k]        = error;
                          errorsMeanErrCorr[i][k]     = error*0.01;
                      }
                    }
                } else if (meson.Contains("Pi0")){
                    Double_t error              = 5.; // 0,75% at low pT, up to 2%
                    if ( additionalNameOutput.CompareTo("EMC7")==0 ){
                      for (Int_t k = 0;k < nPtBins;k++){
                            error              = 1.52;
                            if(ptBins[k]>=9.) error = 5.266 + (-0.855)*ptBins[k] + 0.049*ptBins[k]*ptBins[k];

                            errorsMean[i][k]            = error;
                            errorsMeanErr[i][k]         = error*0.01;
                            errorsMeanCorr[i][k]        = error;
                            errorsMeanErrCorr[i][k]     = error*0.01;
                        }
                    } else if ( additionalNameOutput.CompareTo("EGA")==0 ){
                      for (Int_t k = 0;k < nPtBins;k++){
                          error              = 4.;
                          error += 2*(ptBins[k]-15.);
                          errorsMean[i][k]            = error;
                          errorsMeanErr[i][k]         = error*0.01;
                          errorsMeanCorr[i][k]        = error;
                          errorsMeanErrCorr[i][k]     = error*0.01;
                      }
                    }
                }
            }    
              
              
            // manual smoothing for OpeningAngle errors - variation 1
            if  (nameCutVariationSC[i].CompareTo("OpeningAngle") == 0){
                cout << "Opening Angle smoothing" << endl;
                if ( additionalNameOutput.CompareTo("EMC7")==0 ){
                  for (Int_t k = 0;k < nPtBins;k++){
                      Double_t error              = 1.5;
                      if(ptBins[k]>=10.) error += (0.284)*(ptBins[k]-10.);
                      if( meson.CompareTo("Eta") == 0 ){
                          error   = 0.25;
                          if(ptBins[k]>=15.) error += (0.12)*(ptBins[k]-15.);
                      }else if( meson.CompareTo("Pi0EtaBinning") == 0){
                        Double_t error2 = 0.25;
                        if(ptBins[k]>=15.) error2 += (0.12)*(ptBins[k]-15.);
                        error = sqrt(pow(error,2)+pow(error2,2));
                      }
                      errorsMean[i][k]            = error;
                      errorsMeanErr[i][k]         = error*0.01;
                      errorsMeanCorr[i][k]        = error;
                      errorsMeanErrCorr[i][k]     = error*0.01;
                  }
                } else if ( additionalNameOutput.CompareTo("EGA")==0 ){
                  for (Int_t k = 0;k < nPtBins;k++){
                      Double_t error              = 2.;
                      if(ptBins[k]>=10.) error += (0.284)*(ptBins[k]-10.);
                      if( meson.CompareTo("Eta") == 0 ){
                          error   = 0.25;
                          if(ptBins[k]>=15.) error += (0.12)*(ptBins[k]-15.);
                      }else if( meson.CompareTo("Pi0EtaBinning") == 0){
                        Double_t error2 = 0.25;
                        if(ptBins[k]>=15.) error2 += (0.12)*(ptBins[k]-15.);
                        error = sqrt(pow(error,2)+pow(error2,2));
                      }
                      errorsMean[i][k]            = error;
                      errorsMeanErr[i][k]         = error*0.01;
                      errorsMeanCorr[i][k]        = error;
                      errorsMeanErrCorr[i][k]     = error*0.01;
                  }
                } else{
                  for (Int_t k = 0;k < nPtBins;k++){
                      Double_t error              = 0.5;
                      if( meson.CompareTo("Eta") == 0 ){
                          error   = 0.;
                      }
                      errorsMean[i][k]            = error;
                      errorsMeanErr[i][k]         = error*0.01;
                      errorsMeanCorr[i][k]        = error;
                      errorsMeanErrCorr[i][k]     = error*0.01;
                  }
                }
            }    

            // manual smoothing for minimum cluster energy errors - variation 2
            if (nameCutVariationSC[i].CompareTo("ClusterMinEnergy")==0  ){
                cout << "Cluster minimum energy smoothing" << endl;
                if(meson.CompareTo("Pi0EtaBinning") == 0){
                    Double_t error      = 0;
                    Double_t errorPi0   = 0;
                    Double_t errorEta   = 0;
                    
                    for (Int_t k = 0;k < nPtBins;k++){
                        errorPi0                = 1.2+6./pow(6.,ptBins[k]);
                        errorEta                = 1.2+6./pow(6.,ptBins[k]);
                        if( additionalNameOutput.CompareTo("EMC7") == 0 || additionalNameOutput.CompareTo("EGA") == 0 ) errorEta = 2.5;
                        error                   = TMath::Sqrt(errorPi0*errorPi0+errorEta*errorEta);
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                
                } else if(meson.Contains("Pi0")){
                    Double_t error = 0;
                    for (Int_t k = 0;k < nPtBins;k++){
                        error                   = 1.2+6./pow(6.,ptBins[k]);
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                    if ( additionalNameOutput.CompareTo("EMC7")==0 ){
                      for (Int_t k = 0;k < nPtBins;k++){
                          Double_t error              = 1.2;
                          errorsMean[i][k]            = error;
                          errorsMeanErr[i][k]         = error*0.01;
                          errorsMeanCorr[i][k]        = error;
                          errorsMeanErrCorr[i][k]     = error*0.01;
                      }
                    } else if ( additionalNameOutput.CompareTo("EGA")==0 ){
                      for (Int_t k = 0;k < nPtBins;k++){
                          Double_t error              = 1.2;
                          errorsMean[i][k]            = error;
                          errorsMeanErr[i][k]         = error*0.01;
                          errorsMeanCorr[i][k]        = error;
                          errorsMeanErrCorr[i][k]     = error*0.01;
                      }
                    }
                } else {
                    Double_t error = 0;
                    for (Int_t k = 0;k < nPtBins;k++){
                        error                   = 1.2+6./pow(6.,ptBins[k]);
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                    if ( additionalNameOutput.CompareTo("EMC7")==0 ){
                      for (Int_t k = 0;k < nPtBins;k++){
                          Double_t error              = 2.5;
                          errorsMean[i][k]            = error;
                          errorsMeanErr[i][k]         = error*0.01;
                          errorsMeanCorr[i][k]        = error;
                          errorsMeanErrCorr[i][k]     = error*0.01;
                      }
                    } else if ( additionalNameOutput.CompareTo("EGA")==0 ){
                      for (Int_t k = 0;k < nPtBins;k++){
                          Double_t error              = 2.5;
                          errorsMean[i][k]            = error;
                          errorsMeanErr[i][k]         = error*0.01;
                          errorsMeanCorr[i][k]        = error;
                          errorsMeanErrCorr[i][k]     = error*0.01;
                      }
                    }
                    
                }   
            }
            // manual smoothing for minimum number of cells in cluster errors - variation 3
            if (nameCutVariationSC[i].CompareTo("ClusterNCells")==0 ){
                cout << "Cluster NCells smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 0.5;
                    if (meson.CompareTo("Eta") == 0) error*= 2;
                    else if(meson.CompareTo("Pi0EtaBinning") == 0) error = error + error*2;
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
                    Double_t error              = 0;
                    if (additionalNameOutput.CompareTo("")==0  )
                         error = 1.0;
                    if ( additionalNameOutput.CompareTo("EMC7")==0  || 
                         additionalNameOutput.CompareTo("EGA")==0
                    )
                        error   = 1.25+error;
                      
                    if (meson.CompareTo("Eta") == 0)
                        error   = error*2;
                    if (meson.CompareTo("Pi0EtaBinning") == 0)
                        error   = error+error*2;

                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }
            // manual smoothing for cluster matching errors - variation 5
            if (nameCutVariationSC[i].CompareTo("ClusterTrackMatchingCalo")==0 ){
                cout << "Cluster track matching smoothing" << endl;
                
                if ( meson.Contains("Pi0") && meson.CompareTo("Pi0EtaBinning") != 0){
                    for (Int_t k = 0;k < nPtBins;k++){
                        Double_t error          = 0;
                        if (additionalNameOutput.CompareTo("")==0){
                          error = 2.0 + 16./pow(8.,ptBins[k]);
                        } else if( additionalNameOutput.CompareTo("EMC7")==0 ){
                          error = 3.3;
                          if(ptBins[k]>=13.) error += (0.285)*(ptBins[k]-13.);
                        } else if( additionalNameOutput.CompareTo("EGA")==0 ){
                          error = 3.3;
                          if(ptBins[k]>=13.) error += (0.285)*(ptBins[k]-13.);
                        }
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }   
                } else {
                    for (Int_t k = 0;k < nPtBins;k++){
                        Double_t error          = 0;
                        if (additionalNameOutput.CompareTo("")==0 ||
                            additionalNameOutput.CompareTo("EMC7")==0 || 
                            additionalNameOutput.CompareTo("EGA")==0
                        ){
                          if(meson.CompareTo("Pi0EtaBinning") == 0) {
                            error   = 4.0 + 16./pow(8.,ptBins[k]);
                            if(ptBins[k]>=20.) error += (0.1)*(ptBins[k]-20.);
                          }else {
                            error   = 4.0 + 16./pow(8.,ptBins[k]);
                            if(ptBins[k]>=20.) error += (0.1)*(ptBins[k]-20.);
                          }
                        }    
                        errorsMean[i][k]        = error;
                        errorsMeanErr[i][k]     = 0.01*error;
                        errorsMeanCorr[i][k]    = error;
                        errorsMeanErrCorr[i][k] = 0.01*error;
                    }
                }
            }
           // manual smoothing for cluster shape errors - variation 6
            if (nameCutVariationSC[i].CompareTo("ClusterM02")==0 ){//&& meson.Contains("Pi0")
                cout << "Cluster M02 smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 0.5+(-0.007)*ptBins[k]+(0.013)*ptBins[k]*ptBins[k];
                    if( additionalNameOutput.CompareTo("EMC7")==0 ){
                      error = 0.75;
                      if(ptBins[k]>=12.) error += (0.28)*(ptBins[k]-12.);
                    }else if( additionalNameOutput.CompareTo("EGA") == 0){
                      error = 1.;
                      if(ptBins[k]>=15.) error += (0.8)*(ptBins[k]-15.);
                    }
                    if (meson.CompareTo("Eta") == 0){
                        error   = 0.9+(-0.007)*ptBins[k]+(0.013)*ptBins[k]*ptBins[k];
                        if( additionalNameOutput.CompareTo("EMC7")==0 ){
                          error = 1.2;
                          if(ptBins[k]>=15.) error += (0.01)*(pow((ptBins[k]-15.),2));
                        }else if( additionalNameOutput.CompareTo("EGA") == 0){
                          error = 1.2;
                          if(ptBins[k]>=15.) error += (0.01)*(pow((ptBins[k]-15.),2));
                        }
                    } else if (meson.CompareTo("Pi0EtaBinning") == 0){
                        error   = 1.2+(-0.007)*ptBins[k]+(0.013)*ptBins[k]*ptBins[k];
                        if( additionalNameOutput.CompareTo("EMC7")==0 ){
                          error = 1.8;
                          if(ptBins[k]>=15.) error += (0.0125)*(pow((ptBins[k]-15.),2));
                        }else if( additionalNameOutput.CompareTo("EGA") == 0){
                          error = 1.8;
                          if(ptBins[k]>=15.) error += (0.0125)*(pow((ptBins[k]-15.),2));
                        }
                    }
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }
            // manual smoothing for Material infront of EMCal - variation 7
            if (nameCutVariationSC[i].CompareTo("ClusterMaterialTRD")==0 ){
                cout << "Material smoothing" << endl;
                Double_t error                  = TMath::Sqrt(3.*3.+3.*3.); //(3% for TRD mat, 3% for TOF mat added in quadrature)
                if (meson.CompareTo("Pi0EtaBinning") == 0)
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
                cout << "Cluster energy scale errors smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 0;
//                     if ( ptBins[k] < 1.2 ) {
//                         error = errorsMeanCorr[i][k] + 1;
//                     } else {
                        if (additionalNameOutput.CompareTo("")==0     || 
                            additionalNameOutput.CompareTo("EMC7")==0 ||
                            additionalNameOutput.CompareTo("EGA")==0){
                            error   = 1.5+25./pow(10,ptBins[k]);
                            //if(ptBins[k]>=6.) error += 0.125*(ptBins[k]-6.);
                        }    
//                     } 
                    if (meson.CompareTo("Eta") == 0 || meson.CompareTo("Pi0EtaBinning") == 0){
                      error   = 2.5+60./pow(10,ptBins[k]);
                      //if(ptBins[k]>=6.) error += 0.125*(ptBins[k]-6.);
                    }
//                    if (meson.CompareTo("Pi0EtaBinning") == 0)
//                        error   = 2* error;
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }

            // manual smoothing for cell time uncertainties - variation 9
            if (nameCutVariationSC[i].CompareTo("CellTiming") == 0){
                cout << "Cell time smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 0.;
                    if (additionalNameOutput.CompareTo("")==0 ){
                        error   = 1.5;
                    } else if (additionalNameOutput.CompareTo("EMC7")==0){
                        error   = 1.5;
                    } else if (additionalNameOutput.CompareTo("EGA")==0){
                        error   = 1.5;
                    }    
                    if (meson.CompareTo("Eta") == 0){
                        error   = error*1.5;
                    }
                    if (meson.CompareTo("Pi0EtaBinning") == 0){
                        error   = 0.;
                    }    
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
                
            }    
           
            // manual smoothing for Trigger normalization uncertainties - variation 10
            if (nameCutVariationSC[i].CompareTo("Trigger") == 0){
                cout << "Trigger smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 0.;
                    if (additionalNameOutput.CompareTo("")==0){
                        error   = 0.;
                    } else if (additionalNameOutput.CompareTo("EMC7")==0){
                        error   = TMath::Sqrt(2.0*2.0);
                    } else if (additionalNameOutput.CompareTo("EGA")==0){
                        error   = TMath::Sqrt(2.2*2.2+2.0*2.0);
                    }    
                    if (meson.CompareTo("Pi0EtaBinning") == 0){
                        error   = 0.;
                    }    
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
                
            }    
            // manual smoothing for Efficiency uncertainties - variation 11
            if (nameCutVariationSC[i].CompareTo("Efficiency") == 0){
                cout << "Efficiency smoothing" << endl;
                for (Int_t k = 0;k < nPtBins;k++){
                    Double_t error              = 1.;
                    Double_t errorPi0           = 1.;
                    Double_t errorEta           = 1.;

                    if (additionalNameOutput.CompareTo("")==0 ){
                        errorPi0    = 2.;
                        errorEta    = 4.;
                    } else if (additionalNameOutput.CompareTo("EMC7")==0){
                        errorPi0    = 0.5 + 25/pow(1.5,ptBins[k]) + 2.;
                        errorEta    = 0.5 + 25/pow(1.5,ptBins[k]) + 4.;
                    } else if (additionalNameOutput.CompareTo("EGA")==0){
                        errorPi0    = 5. + 2.;
                        errorEta    = 2. + 25/pow(1.2,ptBins[k]) + 4.;
                    }
                    if (meson.CompareTo("Pi0EtaBinning")==0){
                       error    = TMath::Sqrt(errorPi0*errorPi0+errorEta*errorEta);
                    } else if (meson.CompareTo("Pi0")==0) {
                        error   = errorPi0; 
                    } else {
                        error   = errorEta; 
                    }  
                   
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
                
            }    

            // manual smoothing for ClusterTime uncertainties - variation 12
            if (nameCutVariationSC[i].CompareTo("ClusterTime")==0 ){
                cout << "ClusterTime smoothing" << endl;
                Double_t error                  = 0.25;
                if(additionalNameOutput.CompareTo("EGA")) error = 0.75;
                if (meson.CompareTo("Eta") == 0){
                    error   = 0.5;
                    if(additionalNameOutput.CompareTo("EGA")) error = 1.5;
                }
                if (meson.CompareTo("Pi0EtaBinning") == 0){
                    error   = 0.;
                }
                for (Int_t k = 0;k < nPtBins;k++){
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }

            // manual smoothing for ClusterizationEnergy uncertainties - variation 13
            if (nameCutVariationSC[i].CompareTo("ClusterizationEnergy")==0 ){
                cout << "ClusterizationEnergy smoothing" << endl;
                Double_t error                  = 2.;
                for (Int_t k = 0;k < nPtBins;k++){
                    error = 3.0 + 20./pow(10,ptBins[k]);
//                    if (meson.CompareTo("Eta") == 0){
//                        error   = error*1.25;
//                    }
                    if (meson.CompareTo("Pi0EtaBinning") == 0){
                        error = 4.0 + 20./pow(10,ptBins[k]);
                    }
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }

            // manual smoothing for Periods uncertainties - variation 14
            if (nameCutVariationSC[i].CompareTo("Periods")==0 ){
                cout << "Periods smoothing" << endl;
                Double_t error                  = 3.;
                if (meson.CompareTo("Pi0EtaBinning") == 0){
                    error   = 0.;
                }
                for (Int_t k = 0;k < nPtBins;k++){
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
            }

            // manual smoothing for pi0 in eta binning - variation 15
            if (i==15 && nameCutVariationSC[i].CompareTo("YieldExtraction")==0 ){
                cout << "pi0etabinning smoothing" << endl;
                Double_t error                  = 3.;
                for (Int_t k = 0;k < nPtBins;k++){
                    error = 0.22 + 0.003*ptBins[k] + 0.014*ptBins[k]*ptBins[k];
                    errorsMean[i][k]            = error;
                    errorsMeanErr[i][k]         = error*0.01;
                    errorsMeanCorr[i][k]        = error;
                    errorsMeanErrCorr[i][k]     = error*0.01;
                }
                if ( additionalNameOutput.CompareTo("EMC7")==0 ){
                  for (Int_t k = 0;k < nPtBins;k++){
                        error              = 1.02;
                        if(ptBins[k]>=9.) error = 4.766 + (-0.855)*ptBins[k] + 0.049*ptBins[k]*ptBins[k];

                        errorsMean[i][k]            = error;
                        errorsMeanErr[i][k]         = error*0.01;
                        errorsMeanCorr[i][k]        = error;
                        errorsMeanErrCorr[i][k]     = error*0.01;
                    }
                } else if ( additionalNameOutput.CompareTo("EGA")==0 ){
                  for (Int_t k = 0;k < nPtBins;k++){
                      error              = 16.56 - 2.53*ptBins[k] + 0.1096*ptBins[k]*ptBins[k];
                      errorsMean[i][k]            = error;
                      errorsMeanErr[i][k]         = error*0.01;
                      errorsMeanCorr[i][k]        = error;
                      errorsMeanErrCorr[i][k]     = error*0.01;
                  }
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
        errorsPosCorrMatSummed[l]       = pow(errorsPosCorrSummed[l]+ pow(errorMaterial ,2.) + pow(errorsPosCorr[7][l],2) ,0.5);
        errorsMeanCorrMatSummed[l]      = pow(errorsMeanCorrSummed[l]+ pow(errorMaterial ,2.)+ pow(errorsMeanCorr[7][l],2),0.5);
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
    Double_t heightLegend   = 1.05* 0.035 * (numberCutStudies+3);
    if (numberCutStudies> 7)  
        heightLegend        = 1.05* 0.035 * (numberCutStudies/2+1);
    if (numberCutStudies> 7 && meson.CompareTo("Eta") == 0)
        heightLegend        = 1.05* 0.035 * (numberCutStudies/2+1);
    if (numberCutStudies> 7 && meson.CompareTo("Pi0EtaBinning") == 0)
        heightLegend        = 1.05* 0.035 * (numberCutStudies/2);
    
    // ***************************************************************************************************
    // ****************************** Plot all mean erros separately *************************************
    // ***************************************************************************************************
    
    TCanvas* canvasSysErrMean = new TCanvas("canvasSysErrMean","",200,10,1350,900);// gives the page size
    DrawGammaCanvasSettings( canvasSysErrMean, 0.08, 0.01, 0.015, 0.09);
 
        // create dummy histo
        TH2D *histo2DSysErrMean ;
        if (meson.Contains("Eta") ){
            Double_t max = 20;
            if(additionalNameOutput.CompareTo("EMC7")==0) max = 23;
            if(additionalNameOutput.CompareTo("EGA")==0) max = 26;
            if(meson.CompareTo("Pi0EtaBinning")==0){
              if(additionalNameOutput.CompareTo("EMC7")==0) max = 25;
              if(additionalNameOutput.CompareTo("EGA")==0) max = 25;
            }
            histo2DSysErrMean = new TH2D("histo2DSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,0.,max);
        } else {
            Double_t max = 17;
            if(additionalNameOutput.CompareTo("EGA")==0) max = 23;
            histo2DSysErrMean = new TH2D("histo2DSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,0.,max);
        }
        SetStyleHistoTH2ForGraphs( histo2DSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean systematic Err %", 0.03, 0.04, 0.03, 0.04,
                                1,0.9, 510, 510);
        histo2DSysErrMean->Draw();
        
        // create legend
        TLegend* legendMean = GetAndSetLegend2(minXLegend,maxYLegend-heightLegend,minXLegend+widthLegend,maxYLegend, 30);
        if (numberCutStudies> 7) legendMean->SetNColumns(2);
        
        for(Int_t i = 0;i< numberCutStudies ;i++){
            if ((additionalNameOutput.CompareTo("") == 0) && i == 10){
                cout << "not drawing: " << nameCutVariation[i].Data() << endl;
                continue;
            }    
            if ( meson.CompareTo("Eta") == 0 && i == 1){
                cout << "not drawing: " << nameCutVariation[i].Data() << endl;
                continue;
            }    
            if ( meson.CompareTo("Pi0EtaBinning") == 0 && (i == 7 || i == 10 || i == 9 || i == 12 || i == 14) ){
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
        if (meson.CompareTo("Pi0EtaBinning") == 0){
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
            labelTrig= new TLatex(0.75,0.84,Form("LHC12 - INT7"));
        } else if (additionalNameOutput.CompareTo("EMC7")==0){
            labelTrig= new TLatex(0.7,0.84,Form("LHC12 - EMC L0, INT7"));
        } else if (additionalNameOutput.CompareTo("EGA")==0){
            labelTrig= new TLatex(0.7,0.84,Form("LHC12 - EMC L1-GA, INT7"));
        }
        SetStyleTLatex( labelTrig, 0.038,4);
        labelTrig->Draw();
        
    canvasSysErrMean->Update();
    canvasSysErrMean->SaveAs(Form("SystematicErrorsCalculatedCalo/SysMean_%s_%s%s_%s.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
    
    delete canvasSysErrMean;
    
    
    // ***************************************************************************************************
    // ********************* Plot all mean erros separately after smoothing ******************************
    // ***************************************************************************************************    
    TCanvas* canvasNewSysErrMean = new TCanvas("canvasNewSysErrMean","",200,10,1350,900);// gives the page size
    DrawGammaCanvasSettings( canvasNewSysErrMean, 0.08, 0.01, 0.015, 0.09);
    
        // create dummy histo
        TH2D *histo2DNewSysErrMean ;
        if (meson.Contains("Eta")){
            Double_t max = 20;
            if(additionalNameOutput.CompareTo("EMC7")==0) max = 23;
            if(additionalNameOutput.CompareTo("EGA")==0) max = 26;
            if(meson.CompareTo("Pi0EtaBinning")==0){
              if(additionalNameOutput.CompareTo("EMC7")==0) max = 25;
              if(additionalNameOutput.CompareTo("EGA")==0) max = 25;
            }
            histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,max);
        } else {
            Double_t max = 17;
            if(additionalNameOutput.CompareTo("EGA")==0) max = 23;
            histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,max);
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
            if ((additionalNameOutput.CompareTo("") == 0 || additionalNameOutput.CompareTo("INT7") ==0 ) && i == 10){
                cout << "not drawing: " << nameCutVariation[i].Data() << endl;
                continue;
            }    
            if ( meson.CompareTo("Eta") == 0 && i == 1 && additionalNameOutput.CompareTo("")==0){
                cout << "not drawing: " << nameCutVariation[i].Data() << endl;
                continue;
            }    
            if ( meson.CompareTo("Pi0EtaBinning") == 0 && (i == 7 || i == 10 || i == 9 || i == 12 || i == 14) ){
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
    canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedCalo/SysMeanNewWithMean_%s_%s%s_%s.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
    
    // ***************************************************************************************************
    // ********************* Plot unsmoothed errors with fits ********************************************
    // ***************************************************************************************************    
    for (Int_t cut =0 ;cut < numberCutStudies;cut++ ){
        
        canvasNewSysErrMean->cd();
        histo2DNewSysErrMean->Draw();

            if (bsmooth[cut]) continue;
            cout <<endl << endl<<  "variation: " << cut << " \t"<< nameCutVariation[cut].Data() << endl;
            Double_t minPt = startPtSys;
//             if (additionalNameOutput.CompareTo("EMC1")==0)  minPt = 2.6;
            Double_t maxPt = ptBins[nPtBins-2]+1;
//             if (cut == 13) maxPt = 6;
//             if (cut == 12 || cut == 5) maxPt = 8;
             //if (cut == 6) minPt = 2.4;
             //if (cut == 6) maxPt = 12.5;

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
                
            cout << "red: " << endl;
            meanErrorsCorr[cut]->Fit(pol4,"NRMEX0+","",minPt,maxPt);
            cout << "blue: " << endl;
            meanErrorsCorr[cut]->Fit(pol2,"NRMEX0+","",minPt,maxPt);
            cout << "cyan: " << endl;
            meanErrorsCorr[cut]->Fit(pol1,"NRMEX0+","",minPt,maxPt);
            cout << "black: " << endl;
            meanErrorsCorr[cut]->Fit(pol0,"NRMEX0+","",minPt,maxPt);
            cout << "magenta: " << endl;
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
            
        canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedCalo/SysMeanNewWithMeanSingle_%s_%s%s_%s_Variation%d.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),cut,suffix.Data()));
    }   
        
    // ***************************************************************************************************
    // ********************* Create output files with errors *********************************************
    // ***************************************************************************************************    
    const char *SysErrDatname = Form("SystematicErrorsCalculatedCalo/SystematicErrorEMCEMC_%s_%s%s_%s.dat",meson.Data(),energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data());
    fstream SysErrDat;
    cout << SysErrDatname << endl;
    SysErrDat.open(SysErrDatname, ios::out);
    for (Int_t l=0;l< nPtBins;l++){
        SysErrDat << ptBins[l] << "\t" <<errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
    }

    SysErrDat.close();

    const char *SysErrDatnameMean = Form("SystematicErrorsCalculatedCalo/SystematicErrorAveragedEMCEMC_%s_%s%s_%s.dat",meson.Data(),energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data());
    fstream SysErrDatAver;
    cout << SysErrDatnameMean << endl;
    SysErrDatAver.open(SysErrDatnameMean, ios::out);
    for (Int_t l=0;l< nPtBins;l++){
        SysErrDatAver << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
    }

    const char *SysErrDatnameMeanSingleErr = Form("SystematicErrorsCalculatedCalo/SystematicErrorAveragedSingleEMCEMC_%s_%s%s_%s.dat",meson.Data(),energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data());
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
    Double_t errorsMeanCorrSignalExtraction[nPtBins];
    Double_t errorsMeanCorrClusterEnergy[nPtBins];
    Double_t errorsMeanCorrClusterDescription[nPtBins];
    Double_t errorsMeanCorrCellTiming[nPtBins];
    
    for (Int_t l=0;l< nPtBins;l++){
//      "YieldExtraction"-0,"OpeningAngle"-1, "ClusterMinEnergy"-2, "ClusterNCells"-3, "NonLinearity"-4, "ClusterTrackMatchingCalo" -5, "ClusterM02" -6, "ClusterMaterialTRD" -7, "ClusterEnergyScale" -8
//      "Cell Timing" -9 , "Trigger" - 10, "Efficiency" -11, "ClusterTime" -12, "ClusterizationEnergy" -13, "Periods" -14, "YieldExtraction Pi0EtaBinning" -15
        // grouping:
        // Signal extraction: Yield extraction 0, Open-Angle 1 
        errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+errorsMeanCorr[1][l]*errorsMeanCorr[1][l]);
        if (meson.CompareTo("Pi0EtaBinning") == 0)
            errorsMeanCorrSignalExtraction[l] = TMath::Sqrt(errorsMeanCorr[0][l]*errorsMeanCorr[0][l]+errorsMeanCorr[1][l]*errorsMeanCorr[1][l]+errorsMeanCorr[15][l]*errorsMeanCorr[15][l]);
        // Cluster energy extraction: NonLinearity 4 , Energy Scale 8
        errorsMeanCorrClusterEnergy[l]      = TMath::Sqrt(errorsMeanCorr[4][l]*errorsMeanCorr[4][l]+errorsMeanCorr[8][l]*errorsMeanCorr[8][l]);
        // cluster description in MC: ClusterMinEnergy 2, ClusterNCells 3, ClusterM02 6, Clusterization Energy-13
        errorsMeanCorrClusterDescription[l] = TMath::Sqrt(errorsMeanCorr[2][l]*errorsMeanCorr[2][l]+errorsMeanCorr[3][l]*errorsMeanCorr[3][l]+errorsMeanCorr[6][l]*errorsMeanCorr[6][l]+errorsMeanCorr[13][l]*errorsMeanCorr[13][l]);
        // cell timing: cell timing - -9, clusterTime - 12
        errorsMeanCorrCellTiming[l] = TMath::Sqrt(errorsMeanCorr[9][l]*errorsMeanCorr[9][l]+errorsMeanCorr[12][l]*errorsMeanCorr[12][l]);
    }
        
        
    TGraphErrors* meanErrorsSignalExtraction    = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSignalExtraction ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsClusterDescrip      = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrClusterDescription ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsClusterEnergy       = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrClusterEnergy ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsCellTiming          = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrCellTiming ,ptBinsErr ,errorsMeanErrCorrSummed );
    
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
        if (meson.Contains("Eta") ){
            Double_t max = 20.;
            if(additionalNameOutput.CompareTo("EMC7")==0) max = 23;
            if(additionalNameOutput.CompareTo("EGA")==0) max = 26;
            if(meson.CompareTo("Pi0EtaBinning")==0){
              if(additionalNameOutput.CompareTo("EMC7")==0) max = 25;
              if(additionalNameOutput.CompareTo("EGA")==0) max = 25;
            }
            histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,max);

        } else {
            Double_t max = 17;
            if(additionalNameOutput.CompareTo("EGA")==0) max = 23;
            histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,max);
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
        meanErrorsClusterDescrip->Draw("pX0,csame");
        // Track matching to EMCAL
        DrawGammaSetMarkerTGraphErr(meanErrorsCorr[5], 25, 1.,color[2],color[2]);
        meanErrorsCorr[5]->Draw("pX0,csame");
        // Material infront of EMCAL
        if (meson.CompareTo("Pi0EtaBinning") != 0){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[7], 21, 1.,kMagenta+2,kMagenta+2);
            meanErrorsCorr[7]->Draw("pX0,csame");
        }    
        // Cluster energy description
        DrawGammaSetMarkerTGraphErr(meanErrorsClusterEnergy, 20, 1.,color[4],color[4]);
        meanErrorsClusterEnergy->Draw("pX0,csame");

        // timing errors
        if (meson.CompareTo("Pi0EtaBinning") != 0){
            DrawGammaSetMarkerTGraphErr(meanErrorsCellTiming, 23, 1.,color[5],color[5]);
            meanErrorsCellTiming->Draw("pX0,csame");
        }
        
        // efficiency errors
        if (numberCutStudies>11){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[10], 23, 1.,color[8],color[8]);
            meanErrorsCorr[11]->Draw("pX0,csame");
        }
        // trigger errors
        if (numberCutStudies>10 && !(additionalNameOutput.CompareTo("") == 0 || additionalNameOutput.CompareTo("INT7") ==0 ) && meson.CompareTo("Pi0EtaBinning") != 0){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[10], 25, 1.,color[6],color[6]);
            meanErrorsCorr[10]->Draw("pX0,csame");
        }

        //period error
        if (meson.CompareTo("Pi0EtaBinning") != 0){
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[14], 26, 1.,color[6],color[6]);
            meanErrorsCorr[14]->Draw("pX0,csame");
        }

        legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"signal extraction","p");
        legendSummedMeanNew->AddEntry(meanErrorsClusterDescrip,"cluster description","p");
        legendSummedMeanNew->AddEntry(meanErrorsClusterEnergy,"cluster energy description","p");
        legendSummedMeanNew->AddEntry(meanErrorsCorr[5],"track match. to cluster","p");
        if (meson.CompareTo("Pi0EtaBinning") != 0)legendSummedMeanNew->AddEntry(meanErrorsCorr[7],"mat. infront of EMCal","p");
        if (meson.CompareTo("Pi0EtaBinning") != 0)legendSummedMeanNew->AddEntry(meanErrorsCellTiming,"cell timing","p");
        if (numberCutStudies>11) legendSummedMeanNew->AddEntry(meanErrorsCorr[11],"efficiency","p");
        if (numberCutStudies>10 && !(additionalNameOutput.CompareTo("") == 0) && meson.CompareTo("Pi0EtaBinning") != 0)
            legendSummedMeanNew->AddEntry(meanErrorsCorr[10],"trigger normalization","p");
        if (meson.CompareTo("Pi0EtaBinning") != 0)legendSummedMeanNew->AddEntry(meanErrorsCorr[14],"periods","p");
        DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
        meanErrorsCorrSummedIncMat->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
        legendSummedMeanNew->Draw();
        
        labelMeson->Draw();
        labelCentrality->Draw();
        labelTrig->Draw();
    
    canvasSummedErrMean->Update();
    canvasSummedErrMean->SaveAs(Form("SystematicErrorsCalculatedCalo/SysErrorSummedVisu_%s_%s%s_%s.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));
    
    delete canvasSummedErrMean;
// 
//  
//  
//  const char *SysErrDatnameMeanPaper = Form("SystematicErrorsCalculatedConvCalo/SystematicErrorAveragedPaper_%s_%s%s_%s.dat",meson.Data(),energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data());
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
