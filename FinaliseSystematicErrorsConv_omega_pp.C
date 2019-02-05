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

void FinaliseSystematicErrorsConv_omega_pp( TString nameDataFileErrors      = "",
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
    TString recoMethod                      = ReturnFullTextReconstructionProcess(40);
    energyForOutput.ReplaceAll(".","_");

    // ***************************************************************************************************
    // ******************************* general variable definition  **************************************
    // ***************************************************************************************************
    const Int_t nPtBins                     = numberOfPtBins;
    const Int_t nCuts                       = numberCutStudies;
    Double_t* ptBins                        = 0;
    Double_t* ptBinsErr                     = 0;
    TString nameCutVariationSC[19];
    TString nameCutVariationSCCurrent[19]   = { "RemovePileUp",                    // 0
                                                "Conversion_SinglePtCut",          // 1
                                                "Conversion_ClsTPCCut",            // 2
                                                "Conversion_TPCdEdxCutElectron",   // 3
                                                "Conversion_TPCdEdxCutPion",       // 4
                                                "Conversion_QtMaxCut",             // 5
                                                "Conversion_Chi2GammaCut",         // 6
                                                "Conversion_PsiPair",              // 7
                                                "Conversion_CosinePointingAngle",  // 8
                                                "ChargedPion_ClsTPCCut",           // 9
                                                "ChargedPion_DCACut",              // 10
                                                "ChargedPion_pTCut",               // 11
                                                "ChargedPion_TPCdEdxCutPion",      // 12
                                                "ChargedPion_MassCut",             // 13
                                                "NeutralPion_SelectionWindows",    // 14
                                                "Omega_BackgroundScheme",          // 15
                                                "Omega_RecoEfficiency",            // 16
                                                "Omega_BranchingRatio",            // 17
                                                "YieldExtraction"};                // 18
    TString nameCutVariation[19] = {"pileup","#gamma_{conv} track p_{T}","#gamma_{conv} track N_{cls,TPC}","#gamma_{conv} e^{#pm} PID","#gamma_{conv} #pi rej.","#gamma_{conv} q_{T}","#gamma_{conv} #chi^{2}","#gamma_{conv} #psi_{pair}","#gamma_{conv} cos #theta_{PA}",
                                   "#pi^{#pm} rec. N_{cls,TPC}","#pi^{#pm} rec. DCA","#pi^{#pm}  rec.min p_{T}","#pi^{#pm} rec. PID","M_{#pi^{+}#pi^{-}} cut","#pi^{0} rec. M_{#gamma#gamma} cut","background description",
                                    "reco. efficiency","branching ratio","yield extraction"};
    Color_t color[34] = {kTeal+3,
                         kTeal-5,800,kBlue-7,kCyan+4,kSpring+2,kTeal-7,kAzure-5,kRed-2,
                         kBlue,kRed-2,kBlue,kRed+1,kOrange+7,
                         kPink+8,kGreen+2,kYellow+2,
                         kOrange+2,kBlue+2,
                         kPink-3,
                        kCyan-2,kGray+2,kViolet+1,kAzure+1,kAzure+4,kPink+4,kOrange,404,kPink-6,860,kBlue-4,kRed-4,kGreen-4,kMagenta-2};
    Color_t markerStyle[23]={25,
                             24,28,33,29,24,30,21,23,
                             24,21,22,23,20,
                             26,26,27,
                             26,28,
                             30,
                             31,33,23};

    // Display names
    for (Int_t k =0; k<nCuts; k++ ){
        cout << "variation: " << nameCutVariationSCCurrent[k].Data() << endl;
        cout << "name for writing: " << nameCutVariationSCCurrent[k].Data() << endl;
        cout << "display name: " << nameCutVariation[k].Data() << endl;
    }

    for (Int_t i = 0; i < numberCutStudies; i++){
        nameCutVariationSC[i]               = nameCutVariationSCCurrent[i];
    }

    // ***************************************************************************************************
    // ******************************** Booleans for smoothing *******************************************
    // ***************************************************************************************************
    Bool_t bsmooth[17]                      = { 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,
                                                0,  0, 0, 0 ,0};
    Bool_t bsmoothMBOmega07TeV[19]             = { 1,                           // pileup
                                                   1, 1, 1, 1, 1, 1, 1, 0,      // conversion cuts
                                                   1, 1, 1, 1, 1,               // charged pion
                                                   1,                           // neutral pion
                                                   1, 1, 1,                       // omega
                                                   0};                          // yield extraction
    Bool_t bsmoothMBEta07TeV[19]             = { 1,                           // pileup
                                                   1, 1, 1, 1, 1, 1, 1, 1,      // conversion cuts
                                                   1, 1, 1, 1, 1,               // charged pion
                                                   1,                           // neutral pion
                                                   1, 1, 1,                     // omega
                                                   0};                          // yield extraction
    Bool_t bsmoothMBOmegaToPi07TeV[17]             = { 0, 0, 0, 0,  0, 0, 0, 0, 0,0,
                                                0, 0, 0,  0 ,0 , 0, 0}; // currently not used

    // get smoothing switches needed for meson
    for (Int_t i = 0; i < numberCutStudies; i++){
        if (energy.CompareTo("7TeV") == 0){
            if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Omega")==0){
                bsmooth[i]                      = bsmoothMBOmega07TeV[i];
            } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("Eta")==0){
                bsmooth[i]                      = bsmoothMBEta07TeV[i];
            } else if (additionalNameOutput.CompareTo("") == 0 && meson.CompareTo("OmegaToPi0")==0){
                bsmooth[i]                      = bsmoothMBOmegaToPi07TeV[i]; // not implemented yet
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

    Bool_t foundHisto[nCuts]; // keep track if histo is found or dummy is used
    for(Int_t i = 0; i < nCuts; i++){
       foundHisto[i] = kFALSE;
    }

    TFile* fileErrorInput                   = new TFile(nameDataFileErrors);
    for (Int_t i = 0; i < nCuts; i++){
        TGraphAsymmErrors* graphPosErrors   = NULL;
        TGraphAsymmErrors* graphNegErrors   = NULL;

        if (    nameCutVariationSC[i].CompareTo("RemovePileUp")==0
             || nameCutVariationSC[i].CompareTo("Conversion_SinglePtCut")==0
             || nameCutVariationSC[i].CompareTo("Conversion_ClsTPCCut")==0
             || nameCutVariationSC[i].CompareTo("Conversion_TPCdEdxCutElectron")==0
             || nameCutVariationSC[i].CompareTo("Conversion_TPCdEdxCutPion")==0
             || nameCutVariationSC[i].CompareTo("Conversion_QtMaxCut")==0
             || nameCutVariationSC[i].CompareTo("Conversion_Chi2GammaCut")==0
             || nameCutVariationSC[i].CompareTo("Conversion_PsiPair")==0
             || nameCutVariationSC[i].CompareTo("Conversion_CosinePointingAngle")==0
             || nameCutVariationSC[i].CompareTo("ChargedPionEtaCut")==0
             || nameCutVariationSC[i].CompareTo("ChargedPion_ClsTPCCut")==0
             || nameCutVariationSC[i].CompareTo("ChargedPion_DCACut")==0
             || nameCutVariationSC[i].CompareTo("ChargedPion_pTCut")==0
             || nameCutVariationSC[i].CompareTo("ChargedPion_TPCdEdxCutPion")==0
             || nameCutVariationSC[i].CompareTo("ChargedPion_MassCut")==0
             || nameCutVariationSC[i].CompareTo("NeutralPion_SelectionWindows")==0
             || nameCutVariationSC[i].CompareTo("Omega_BackgroundScheme")==0
             || nameCutVariationSC[i].CompareTo("Omega_NumberBckEvents")==0
             || nameCutVariationSC[i].CompareTo("Omega_RapidityCut")==0){

            TString nameGraphPos            = Form("%s_SystErrorRelPos_%spp",meson.Data(),nameCutVariationSC[i].Data()  );
            TString nameGraphNeg            = Form("%s_SystErrorRelNeg_%spp",meson.Data(),nameCutVariationSC[i].Data()  );
            cout << "name graph pos =" << nameGraphPos.Data() << endl;
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
            cout << "graphpos= " << graphPosErrors;
        } else if (nameCutVariationSC[i].CompareTo("YieldExtraction")==0){
            TString nameGraphPos            = Form("%s_SystErrorRelPos_%s_pp",meson.Data(),nameCutVariationSC[i].Data()  );
            TString nameGraphNeg            = Form("%s_SystErrorRelNeg_%s_pp",meson.Data(),nameCutVariationSC[i].Data()  );
            cout << "name graph pos =" << nameGraphPos.Data() << endl;
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
            cout << "graphpos= " << graphPosErrors;
        } else {
            // name is not in list
            printf("WARNING: trying to get %s, but no rule is known. Will try to get it anyways ... \n",nameCutVariationSC[i].Data());
            TString nameGraphPos            = Form("%s_SystErrorRelPos_%spp",meson.Data(),nameCutVariationSC[i].Data() );
            TString nameGraphNeg            = Form("%s_SystErrorRelNeg_%spp",meson.Data(),nameCutVariationSC[i].Data() );
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());
        }

        // check now if getting the histogram actually worked
        if(graphPosErrors==NULL){
            printf("WARNING: Didn't find %s graph in file %s!\n",nameCutVariationSC[i].Data(),fileErrorInput->GetName());
            printf("WARNING: Will use yield extraction histogram as dummy instead!");

            TString nameGraphPos            = Form("%s_SystErrorRelPos_%s_pp",meson.Data(),"YieldExtraction"  );
            TString nameGraphNeg            = Form("%s_SystErrorRelNeg_%s_pp",meson.Data(),"YieldExtraction"  );
            graphPosErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphPos.Data());
            graphNegErrors                  = (TGraphAsymmErrors*)fileErrorInput->Get(nameGraphNeg.Data());

            if(graphPosErrors==NULL){ // still didnt find histo
                printf("ERROR: Could't find yield extraction histo either, returning ...");
                return;
            }

            // Set everything to zero
            printf("INFO: Using YieldExtraction graph for %s and setting all values to zero ...\n",nameCutVariationSC[i].Data());
            for (Int_t k = 0; k< graphPosErrors->GetN(); k++){
                graphPosErrors->SetPoint(k, graphPosErrors->GetX()[k],0);
                graphPosErrors->SetPointEYhigh (k, 0);
                graphPosErrors->SetPointEYlow (k, 0);
                graphNegErrors->SetPoint(k, graphNegErrors->GetX()[k],0);
                graphNegErrors->SetPointEYhigh (k, 0);
                graphNegErrors->SetPointEYlow (k, 0);
            }
            graphPosErrors->SetName(Form("%s_SystErrorRelPos_%spp",meson.Data(),nameCutVariationSC[i].Data()));
            graphNegErrors->SetName(Form("%s_SystErrorRelNeg_%spp",meson.Data(),nameCutVariationSC[i].Data()));
        } else{
            printf("INFO: Loaded %s.graph succesfully!.. \n",nameCutVariationSC[i].Data());
            foundHisto[i]=kTRUE;
        }


        // remove points below minPt
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

            // pileup
            if (nameCutVariationSC[i].CompareTo("RemovePileUp")==0 ){
                minPt       = startPtSys;
                if (!energy.CompareTo("7TeV"))
                    errorReset  = 0.5; // 5% higher yield * 90% efficiency

                for (Int_t k = 0; k < nPtBins; k++){
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }

            // min pt of conversion track
            } else if (nameCutVariationSC[i].CompareTo("Conversion_SinglePtCut")==0 ){
                minPt       = startPtSys;

                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        errorReset = 1.8;
                        if(ptBins[k]>4)errorReset  = 0.5; // 1.85% with const smoothing * 0.68 for 1 sigma error
                    }
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }


            // min cls in TPC for conversion track
            } else if (nameCutVariationSC[i].CompareTo("Conversion_ClsTPCCut")==0 ){
                minPt       = startPtSys;
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        errorReset = 2.0;
                    }
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }

            // PID cut for conversion electron
            } else if (nameCutVariationSC[i].CompareTo("Conversion_TPCdEdxCutElectron")==0 ){
                minPt       = startPtSys;
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        errorReset = 1.0;
                    }
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }

            // rejection of pions when looking at conversion tracks
            } else if (nameCutVariationSC[i].CompareTo("Conversion_TPCdEdxCutPion")==0 ){
                minPt       = startPtSys;
                if (!energy.CompareTo("7TeV"))
                    errorReset  = 0.2; //

                for (Int_t k = 0; k < nPtBins; k++){
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }

            // conversion qt cut
            } else if (nameCutVariationSC[i].CompareTo("Conversion_QtMaxCut")==0 ){
                minPt       = startPtSys;
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        errorReset = 1.8;
                    }
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }

            // conversion chi2 cut
            } else if (nameCutVariationSC[i].CompareTo("Conversion_Chi2GammaCut")==0 ){
                minPt       =startPtSys;
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        errorReset = 1.5;
                    }
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }

            // conversion psipair
            } else if (nameCutVariationSC[i].CompareTo("Conversion_PsiPair")==0 ){
                minPt       = startPtSys;
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        errorReset = 1.5;
                    }
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }

            // conversion cosine pointing angle
            } else if (nameCutVariationSC[i].CompareTo("Conversion_CosinePointingAngle")==0 ){
                minPt       = startPtSys;
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        errorReset = 1.0;
                    }
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }

            // charged pion track cls in TPC
            } else if (nameCutVariationSC[i].CompareTo("ChargedPion_ClsTPCCut")==0 ){
                minPt       = startPtSys;
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        if(ptBins[k]<4.){
                            errorReset =  12.4986+-4.85793*ptBins[k]+0.584105*pow(ptBins[k],2);
                        } else{
                            errorReset = 2.3;
                        }
                    }
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }

                // charged pion PID
            } else if (nameCutVariationSC[i].CompareTo("ChargedPion_TPCdEdxCutPion")==0 ){
                minPt       = startPtSys;
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        errorReset = 3.0;
                    }
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }

            // charged pion PID
            } else if (nameCutVariationSC[i].CompareTo("ChargedPion_DCACut")==0 ){
                minPt       = startPtSys;
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        errorReset = 4.5;
                    }
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }
            // charged pion PID
            } else if (nameCutVariationSC[i].CompareTo("ChargedPion_pTCut")==0 ){
                minPt       = startPtSys;
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        errorReset = 4.;
                    }
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }
            // charged pion pair mass window
            } else if (nameCutVariationSC[i].CompareTo("ChargedPion_MassCut")==0 ){
                minPt       = startPtSys;
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        errorReset = 4.0;
                    }
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }

            // two gamma mass window
            } else if (nameCutVariationSC[i].CompareTo("NeutralPion_SelectionWindows")==0 ){
                minPt       = startPtSys;
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        errorReset = 9.;
                    }
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }

            // background description
            } else if (nameCutVariationSC[i].CompareTo("Omega_BackgroundScheme")==0 ){
                minPt       = startPtSys;
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        errorReset = 5.2;
                        if(ptBins[k]<5.0){
                            errorReset = 15.9952+-4.41681*ptBins[k]+0.450559*pow(ptBins[k],2);
                        }
                    }
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }
            // background description
            } else if (nameCutVariationSC[i].CompareTo("Omega_BranchingRatio")==0 ){
                minPt       = startPtSys;
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        errorReset = 0.8;
                    }
                    if (ptBins[k] > minPt){
                        errorsMean[i][k]            = errorReset;
                        errorsMeanErr[i][k]         = errorReset*0.01;
                        errorsMeanCorr[i][k]        = errorReset;
                        errorsMeanErrCorr[i][k]     = errorReset*0.01;
                    }
                }
            // Yield Extraction
            } else if (nameCutVariationSC[i].CompareTo("YieldExtraction")==0 ){
                minPt       = startPtSys;
                for (Int_t k = 0; k < nPtBins; k++){
                    if (!energy.CompareTo("7TeV")){
                        errorReset = 17.+0.2*pow(ptBins[k]-5.5,2);
                    }
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
    if (numberCutStudies> 11)
        widthLegend         = 0.75;
    Double_t heightLegend   = 1.05* 0.04 * (numberCutStudies+3);
    if (numberCutStudies> 7)
        heightLegend        = 1.05* 0.03 * (numberCutStudies/2+1);
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
            DrawGammaSetMarkerTGraphErr(meanErrors[i], markerStyle[i], 0.4,color[i],color[i]);
            meanErrors[i]->Draw("pE0,csame");
            legendMean->AddEntry(meanErrors[i],nameCutVariation[i].Data(),"p");
        }
        legendMean->SetMargin(0.05);
        legendMean->Draw();

        // plot labeling
        TLatex *labelMeson;
        if (meson.CompareTo("EtaToPi0") == 0){
            labelMeson= new TLatex(0.95,0.89,Form("#eta/#pi^{0} rec. #gamma_{conv}"));
        } else if (meson.Contains("Omega")){
            labelMeson= new TLatex(0.95,0.89,Form("#omega #rightarrow #pi^{+}#pi^{-}#pi^{0}"));
        } else {
            labelMeson= new TLatex(0.95,0.89,Form("#eta #rightarrow #pi^{+}#pi^{-}#pi^{0}"));
        }
        SetStyleTLatex( labelMeson, 0.038,4,1,42,kTRUE, 311);
        labelMeson->SetTextAlign(31);
        labelMeson->Draw();

        TLatex *labelCentrality = new TLatex(0.95,0.93,Form("%s",collisionSystem.Data() ));
        SetStyleTLatex( labelCentrality, 0.038,4,1,42,kTRUE, 11);
        labelCentrality->SetTextAlign(31);
        labelCentrality->Draw();

        TLatex *labelRecoMethod = new TLatex(0.95,0.84,Form("%s",recoMethod.Data() ));
        SetStyleTLatex( labelRecoMethod, 0.038,4,1,42,kTRUE, 11);
        labelRecoMethod->SetTextAlign(31);
        labelRecoMethod->Draw();

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
            histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,45.);
        } else if ( meson.Contains("Eta")){
            histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,70.);
        } else {
            histo2DNewSysErrMean = new TH2D("histo2DNewSysErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,45.);
        }
        SetStyleHistoTH2ForGraphs( histo2DNewSysErrMean, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %", 0.03, 0.04, 0.03, 0.04,
                                1,0.9, 510, 510);
        histo2DNewSysErrMean->Draw();

        // create legend
        TLegend* legendMeanNew = NULL;
        if(meson.Contains("Eta")){
            legendMeanNew = GetAndSetLegend2(minXLegend-0.08,maxYLegend-heightLegend,minXLegend+(0.8*widthLegend),maxYLegend, 30);
        } else{
            legendMeanNew = GetAndSetLegend2(minXLegend,maxYLegend-heightLegend,minXLegend+widthLegend,maxYLegend, 30);
        }
        if (numberCutStudies> 7) legendMeanNew->SetNColumns(2);
        if (numberCutStudies> 11) legendMeanNew->SetNColumns(3);

        for(Int_t i = 0;i< numberCutStudies ;i++){
            cout << i << "\t"<< additionalNameOutput.Data() << endl;
            DrawGammaSetMarkerTGraphErr(meanErrorsCorr[i], markerStyle[i], 1.,color[i],color[i]);
            meanErrorsCorr[i]->Draw("p,csame");
            legendMeanNew->AddEntry(meanErrorsCorr[i],nameCutVariation[i].Data(),"p");
        }
        if (meson.CompareTo("EtaToPi0")){
            DrawGammaSetMarkerTGraphErr(graphMaterialError, GetMarkerStyleSystematics( "Material"), 1., GetColorSystematics( "Material" ), GetColorSystematics( "Material"));
            graphMaterialError->Draw("p,csame");
            legendMeanNew->AddEntry(graphMaterialError,"material","p");
        }
        DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
        meanErrorsCorrSummedIncMat->Draw("p,csame");
        legendMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
        legendMeanNew->SetMargin(0.05);
        legendMeanNew->Draw();

        // labeling
        labelMeson->Draw();
        labelCentrality->Draw();
        labelRecoMethod->Draw();

    canvasNewSysErrMean->Update();
    canvasNewSysErrMean->SaveAs(Form("SystematicErrorsCalculatedConv/SysMeanNewWithMean_%s_%s%s_%s.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));

    // ***************************************************************************************************
    // ********************* Plot unsmoothed errors with fits ********************************************
    // ***************************************************************************************************
    for (Int_t cut =0 ;cut < numberCutStudies;cut++ ){

        canvasNewSysErrMean->cd();
        histo2DNewSysErrMean->GetYaxis()->SetRangeUser(0.,40.);
        histo2DNewSysErrMean->Draw();

            if (!showBeforeSmoothing) continue;

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

            Double_t fitminPt = minPt;
            Double_t fitmaxPt = maxPt;

            meanErrorsCorrBefore[cut]->Fit(pol4,"NRMEX0+","",fitminPt,fitmaxPt);
            meanErrorsCorrBefore[cut]->Fit(pol2,"NRMEX0+","",fitminPt,fitmaxPt);
            meanErrorsCorrBefore[cut]->Fit(pol1,"NRMEX0+","",fitminPt,fitmaxPt);
            meanErrorsCorrBefore[cut]->Fit(pol0,"NRMEX0+","",fitminPt,fitmaxPt);
            meanErrorsCorrBefore[cut]->Fit(bla,"NRMEX0+","",fitminPt,fitmaxPt);

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
    for (Int_t l=0; l< nPtBins; l++){
        SysErrDat << ptBins[l] <<"\t" << errorsNegCorrMatSummed[l] << "\t" <<errorsPosCorrMatSummed[l] << "\t"  <<errorsNegCorrSummed[l] << "\t" <<errorsPosCorrSummed[l]  << endl;
    }
    SysErrDat.close();

    const char *SysErrDatnameMean = Form("SystematicErrorsCalculatedConv/SystematicErrorAveragedPCM_%s_%s%s_%s.dat", meson.Data(), energyForOutput.Data(), additionalNameOutput.Data(), dateForOutput.Data());
    fstream SysErrDatAver;
    cout << SysErrDatnameMean << endl;
    SysErrDatAver.open(SysErrDatnameMean, ios::out);

    for (Int_t l=0; l< nPtBins; l++){
        SysErrDatAver  << ptBins[l] << "\t" << "-"<< errorsMeanCorrMatSummed[l] << "\t" <<errorsMeanCorrMatSummed[l] << "\t"  << "-"<< errorsMeanCorrSummed[l] << "\t" <<errorsMeanCorrSummed[l]  << endl;
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
    for (Int_t l=0;l< nPtBins;l++){
        SysErrDatAverSingle << ptBins[l] << "\t";
        for (Int_t i= 0; i< numberCutStudies; i++){
            if(!meson.CompareTo("EtaToPi0")&&i==1)
                continue;
            SysErrDatAverSingle << errorsMeanCorr[i][l] << "\t";
        }
        if(meson.CompareTo("EtaToPi0"))
            SysErrDatAverSingle << errorsMat[l] << "\t";
        SysErrDatAverSingle << errorsMeanCorrMatSummed[l] << endl;
    }
    SysErrDatAverSingle.close();

    // ***************************************************************************************************
    // ********************* Group errors according to topic *********************************************
    // ***************************************************************************************************
    Double_t errorsMeanCorrChargedPionReco[nPtBins];
    Double_t errorsMeanCorrSignalExtraction[nPtBins];
    Double_t errorsMeanCorrPi0Reco[nPtBins];
    Double_t errorsMeanCorrPhotonReco[nPtBins];
    Double_t errorsMeanCorrPileup[nPtBins];
    Double_t errorsMeanCorrOmegaReco[nPtBins];

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


        // Signal extraction: Yield extraction, Alpha ,BG, YieldExtraction, , "ChargedPion_MassCut"

        errorsMeanCorrSignalExtraction[l]   =   TMath::Sqrt(pow(errorsMeanCorr[13][l],2)+    // charged pion mass cut
                                                            pow(errorsMeanCorr[15][l],2)+    // background
                                                            pow(errorsMeanCorr[16][l],2)+    // reco efficiency
                                                            pow(errorsMeanCorr[17][l],2)+    // branching ratio
                                                            pow(errorsMeanCorr[18][l],2));   // YieldExtraction

        // charged pion cuts: "ChargedPion_ClsTPCCut","ChargedPion_TPCdEdxCutPion"
        errorsMeanCorrChargedPionReco[l]                =   TMath::Sqrt(pow(errorsMeanCorr[9][l],2)   +    // ClsTPC
                                                                        pow(errorsMeanCorr[10][l],2)  +    // DCA
                                                                        pow(errorsMeanCorr[11][l],2)  +    // pT
                                                                        pow(errorsMeanCorr[12][l],2));   // dEdxPi


        // photon reco: SinglePtCut, ClsTPCCut, TPCdEdxCutElectron, TPCdEdxCutPion, Conversion_QtMaxCut, Chi2GammaCut, PsiPair, CosinePointingAngle
        errorsMeanCorrPhotonReco[l]         =   TMath::Sqrt(pow(errorsMeanCorr[1][l],2)+    // SinglePt
                                                            pow(errorsMeanCorr[2][l],2)+    // ClsTPCCut
                                                            pow(errorsMeanCorr[3][l],2)+    // TPCdEdxElectron
                                                            pow(errorsMeanCorr[4][l],2)+    // TPCdEdxPion
                                                            pow(errorsMeanCorr[5][l],2)+    // Qt
                                                            pow(errorsMeanCorr[6][l],2)+    // Chi2
                                                            pow(errorsMeanCorr[7][l],2)+    // PsiPair
                                                            pow(errorsMeanCorr[8][l],2));   // CosineAngle

        // pi0 reco: pT, alpha, mass window
        errorsMeanCorrPi0Reco[l]           =   errorsMeanCorr[14][l];


        // pileup
        if(!energy.CompareTo("8TeV")||!energy.CompareTo("7TeV")){
            errorsMeanCorrPileup[l] = errorsMeanCorr[0][l];
        }
    }
    TGraphErrors* meanErrorsChargedPionReco                 = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrChargedPionReco ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsPhotonReco          = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPhotonReco ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsSignalExtraction    = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrSignalExtraction ,ptBinsErr ,errorsMeanErrCorrSummed );
    TGraphErrors* meanErrorsPi0Reco           = new TGraphErrors(nPtBins,ptBins ,errorsMeanCorrPi0Reco ,ptBinsErr ,errorsMeanErrCorrSummed );
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
            histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,45.);
        } else if (meson.Contains("Eta")){
            histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,60.);

        } else {
            histo2DSummedErrMean = new TH2D("histo2DSummedErrMean", "", 20,0.,ptBins[nPtBins-1]+ptBinsErr[nPtBins-1],1000.,-0.5,45.);
        }
        SetStyleHistoTH2ForGraphs( histo2DSummedErrMean, "#it{p}_{T} (GeV/#it{c})", "mean smoothed systematic Err %", 0.03, 0.04, 0.03, 0.04,
                                1,0.9, 510, 510);
        histo2DSummedErrMean->Draw();

        // create legend
        TLegend* legendSummedMeanNew = GetAndSetLegend2(minXLegend2,maxYLegend2-heightLegend2,minXLegend2+widthLegend2,maxYLegend2, 30);
        legendSummedMeanNew->SetNColumns(2);
        legendSummedMeanNew->SetMargin(0.1);

        // Signal extraction error
        DrawGammaSetMarkerTGraphErr(meanErrorsSignalExtraction, 20, 1.,kBlue,kBlue);
        meanErrorsSignalExtraction->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsSignalExtraction,"signal extraction","p");

        DrawGammaSetMarkerTGraphErr(meanErrorsChargedPionReco, 21, 1.,kRed-2,kRed-2);
        meanErrorsChargedPionReco->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsChargedPionReco,"#pi^{#pm} measurement","p");

        DrawGammaSetMarkerTGraphErr(meanErrorsPi0Reco, 22, 1.,kOrange+7,kOrange+7);
        meanErrorsPi0Reco->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsPi0Reco,"#pi^{0} reconstruction","p");

        DrawGammaSetMarkerTGraphErr(meanErrorsPhotonReco, 23, 1.,kGreen+2,kGreen+2);
        meanErrorsPhotonReco->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsPhotonReco,"#gamma reconstruction","p");

        DrawGammaSetMarkerTGraphErr(meanErrorsPileup, 25, 1.,kCyan-2,kCyan-2);
        meanErrorsPileup->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsPileup,"pileup","p");

        DrawGammaSetMarkerTGraphErr(graphMaterialError, 26, 1.,kViolet+1,kViolet+1);
        graphMaterialError->Draw("p,csame");
        legendSummedMeanNew->AddEntry(graphMaterialError,"material","p");

        DrawGammaSetMarkerTGraphErr(meanErrorsCorrSummedIncMat, 20, 1.,kBlack,kBlack);
        meanErrorsCorrSummedIncMat->Draw("p,csame");
        legendSummedMeanNew->AddEntry(meanErrorsCorrSummedIncMat,"quad. sum.","p");
        legendSummedMeanNew->Draw();

        labelMeson->Draw();
        labelCentrality->Draw();
        labelRecoMethod->Draw();

    canvasSummedErrMean->Update();
    canvasSummedErrMean->SaveAs(Form("SystematicErrorsCalculatedConv/SysErrorSummedVisu_%s_%s%s_%s.%s",meson.Data(), energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data(),suffix.Data()));

    delete canvasSummedErrMean;

    const char *SysErrDatnameMeanPaper = Form("SystematicErrorsCalculatedConv/SystematicErrorAveragedPCMPaper_%s_%s%s_%s.dat",meson.Data(),energyForOutput.Data(),additionalNameOutput.Data(),dateForOutput.Data());
    fstream SysErrDatAverPaper;
    cout << SysErrDatnameMeanPaper << endl;
    SysErrDatAverPaper.open(SysErrDatnameMeanPaper, ios::out);
    SysErrDatAverPaper  << "p_{T}" << "\t Material \t Signal Extraction \t charged pion reco \t photon reco \t pi0 reco  \t summed" <<  endl;
    for (Int_t l=0; l< nPtBins; l++){
        SysErrDatAverPaper << ptBins[l] <<"\t" << errorsMat[l] << "\t" << errorsMeanCorrSignalExtraction[l] << "\t" << errorsMeanCorrChargedPionReco[l]<< "\t" << errorsMeanCorrPhotonReco[l]<< "\t" <<errorsMeanCorrPi0Reco[l] <<"\t" << errorsMeanCorrMatSummed[l]<< endl;
    }
    SysErrDatAverPaper.close();

    printf("-----------------------------------SUMMARY----------------------------------\n");
    printf("%-40s \t %-20s \t %-15s \n","Name","Found Graph in File?","Did smoothing?");
    for(Int_t i = 0; i <nCuts;i++){
        if(meson.CompareTo("Eta")==0){
            printf("%-40s \t %-20d \t %-15d \n",nameCutVariationSC[i].Data(),foundHisto[i],bsmoothMBEta07TeV[i]);
        } else {
            printf("%-40s \t %-20d \t %-15d \n",nameCutVariationSC[i].Data(),foundHisto[i],bsmoothMBOmega07TeV[i]);
        }
    }

}
