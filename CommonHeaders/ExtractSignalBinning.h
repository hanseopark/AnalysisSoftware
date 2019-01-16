// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

#ifndef GAMMACONV_ExtractSignalBinning
#define GAMMACONV_ExtractSignalBinning

    #include <iostream>
    #include <stdio.h>
    using namespace std;

    #include "ConversionFunctionsBasicsAndLabeling.h"
    #include "ExtractSignalBinningpp900GeV.h"
    #include "ExtractSignalBinningpp2760GeV.h"
    #include "ExtractSignalBinningpp5TeV.h"
    #include "ExtractSignalBinningpp7TeV.h"
    #include "ExtractSignalBinningpp8TeV.h"
    #include "ExtractSignalBinningpp13TeV.h"
    #include "ExtractSignalBinningpPb5TeV.h"
    #include "ExtractSignalBinningpPb8TeV.h"
    #include "ExtractSignalBinningPbPb2760GeV.h"
    #include "ExtractSignalBinningPbPb5TeV.h"
    #include "ExtractSignalBinningXeXe5440GeV.h"

    Int_t fStartPtBin                               = 0;
    Int_t fColumn                                   = 0;
    Int_t fRow                                      = 0;
    Int_t fNBinsPt                                  = 0;
    Double_t *fBinsPt                               = NULL;
    Int_t* fNRebin                                  = NULL;
    Int_t fNBinsClusterPt                           = 0;
    Double_t *fBinsClusterPt                        = NULL;
    Int_t fNBinsPtDCAzDist                          = 0;
    Double_t *fBinsPtDCAzDist                       = NULL;
    Int_t fExampleBin                               = 0;
    Double_t fExampleBinScaleFac                    = 1.0;
    Int_t fNBinsPeakPt                              = 12;
    Int_t nIterBGFit                                = 10;
    TString optionBGSmoothingStandard               = "BackDecreasingWindow,BackSmoothing3";
    TString optionBGSmoothingVar1                   = "BackDecreasingWindow,BackSmoothing5";
    TString optionBGSmoothingVar2                   = "BackDecreasingWindow,BackSmoothing7";
    Double_t fMaxYFracBGOverIntHist                 = 30;
    Double_t fBGFitRange_SubPiZero[2]               = {0,0};
    Double_t fBGFitRange_FixedPzPiZero[2]           = {0,0};


    //****************************************************************************************************
    //****************** Pt binning for Inter/Extrapolations *********************************************
    //****************************************************************************************************
    Double_t fBinsInterAndExtrapolation[50]          = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                                         1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                                                         2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                                         4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0,10.0,11.0,
                                                        12.0,13.0,14.0,15.0,16.0,18.0,20.0,25.0,30.0};
    Double_t fBinsInterAndExtrapolationFine[62]      = { 0.0, 0.1,0.12,0.14,0.16,0.18, 0.2,0.25, 0.3,0.35,
                                                         0.4,0.45, 0.5,0.55, 0.6,0.65, 0.7,0.75, 0.8,0.85,
                                                         0.9,0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
                                                         1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
                                                         3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0,
                                                         9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,18.0,20.0,
                                                        25.0,30.0};

    //*************************************************************************************************
    //******************** CopyVectorToArray as helper function for Get/InitializeBinning *************
    //*************************************************************************************************
    // Function that COMPLETELY copies a vector containing bin definition to a C array
    // !!! Use return value to set binningMax, like so:
    //     binningMax = CopyVectorToArray( fBinsEtaPrime7TeVPt, binning );
    template<class TYPE_t>
    size_t CopyVectorToArray( vector<TYPE_t> vectorFrom, TYPE_t* arrayTo ) {
        for( size_t i=0; i<vectorFrom.size(); ++i ) arrayTo[i] = vectorFrom[i];
        return vectorFrom.size()-1;
    }
    // Function that copies a vector containing bin definition to a C array, UP TO A USER-DEFINED SIZE
    // !!! Use return value to set maxNBins, like so:
    //     maxNBins = CopyVectorToArray( fBinsEtaPrime7TeVPt, binning, 10 );
    template<class TYPE_t>
    size_t CopyVectorToArray( vector<TYPE_t> vectorFrom, TYPE_t* arrayTo, size_t maxNBins ) {
        if( maxNBins>=vectorFrom.size() ) {
            cout << "WARNING: Max bin too large (" << maxNBins << "), so set to maximum number of bins (" << vectorFrom.size()-1 << ")" << endl;
            maxNBins = vectorFrom.size()-1;
        }
        for( size_t i=0; i<=maxNBins; ++i ) arrayTo[i] = vectorFrom[i];
        return maxNBins;
    }
    // Function that copies a vector, with or without a user-defined maximum, and also sets binningMax
    // !!! Use return value to set binningMax, like so:
    //     maxNBins = CopyVectorToArray( binningMax, fBinsEtaPrime7TeVPt, binning );
    //     maxNBins = CopyVectorToArray( binningMax, fBinsEtaPrime7TeVPt, binning, 10 );
    template<class TYPE_t>
    size_t CopyVectorToArray( Int_t& binningMax, vector<TYPE_t> vectorFrom, TYPE_t* arrayTo ) {
        binningMax = CopyVectorToArray( vectorFrom, arrayTo );
        cout << "  binningMax = " << binningMax << endl;
        return binningMax;
    }
    template<class TYPE_t>
    size_t CopyVectorToArray( Int_t& binningMax, vector<TYPE_t> vectorFrom, TYPE_t* arrayTo, size_t maxNBins ) {
        binningMax = CopyVectorToArray( vectorFrom, arrayTo, maxNBins );
        cout << "  binningMax = " << binningMax << endl;
        cout << "  maxNBins   = " << maxNBins << endl;
        return binningMax;
    }


    //*************************************************************************************************
    //*********************  determine optimum number of rows and columns *****************************
    //*************************************************************************************************
    void GetOptimumNColumnsAndRows( Int_t totBins, Int_t startBin, Int_t &columns, Int_t &rows ) {
        if ( (totBins+1-startBin) < 5){
            columns     = 2;
            rows        = 2;
        } else if ( (totBins+1-startBin) < 7){
            columns     = 3;
            rows        = 2;
        } else if ( (totBins+1-startBin) < 9){
            columns     = 4;
            rows        = 2;
        } else if ( (totBins+1-startBin) < 10){
            columns     = 3;
            rows        = 3;
        } else if ( (totBins+1-startBin) < 13){
            columns     = 4;
            rows        = 3;
        } else if ( (totBins+1-startBin) < 16){
            columns     = 5;
            rows        = 3;
        } else if ( (totBins+1-startBin) < 21){
            columns     = 5;
            rows        = 4;
        } else if ( (totBins+1-startBin) < 25){
            columns     = 6;
            rows        = 4;
        } else if ( (totBins+1-startBin) < 31){
            columns     = 6;
            rows        = 5;
        } else if ( (totBins+1-startBin) < 36){
            columns     = 7;
            rows        = 5;
        } else if ( (totBins+2-startBin) < 41){
            columns     = 8;
            rows        = 5;
        } else if ( (totBins+2-startBin) < 46){
            columns     = 9;
            rows        = 5;
        } else if ( (totBins+2-startBin) < 61){
            columns     = 10;
            rows        = 6;
        } else if ( (totBins+2-startBin) < 67){
            columns     = 11;
            rows        = 6;
        } else if ( (totBins+2-startBin) < 73){
            columns     = 12;
            rows        = 6;
        } else if ( (totBins+2-startBin) < 78){
            columns     = 11;
            rows        = 7;
        } else if ( (totBins+2-startBin) < 84){
            columns     = 12;
            rows        = 7;
        } else if ( (totBins+2-startBin) < 92){
            columns     = 13;
            rows        = 7;
        } else if ( (totBins+2-startBin) < 97){
            columns     = 12;
            rows        = 8;
        } else if ( (totBins+2-startBin) < 109){
            columns     = 12;
            rows        = 9;
        } else if ( (totBins+2-startBin) < 121){
            columns     = 12;
            rows        = 10;
        } else if ( (totBins+2-startBin) < 133){
            columns     = 12;
            rows        = 11;
        } else {
            columns     = 12;
            rows        = 12;
        }
        cout << "nColumns: " << columns << "\t nRows: "  << rows << "\t nTotbins: " << (totBins+1-startBin) << endl;
    }

    //*************************************************************************************************
    //******************** Initialize Single bin for invariant mass plot ******************************
    //*************************************************************************************************
    Int_t ReturnSingleInvariantMassBinPlotting(
        TString meson,
        TString energy,
        Int_t mode,
        Int_t trigger,
        Double_t &scaleFac,
        Int_t triggerSet                    = -1,
        TString directPhotonRunningOption   = "",
        TString centrality                  = ""
    ) {

        // Heavy meson fix
        if( mode>=100 ) mode -= 100;

        if (triggerSet != -1){
            if (energy.CompareTo("2.76TeV") == 0){
                if (triggerSet == 1)
                    trigger     = 52;
                if (triggerSet == 2)
                    trigger     = 85;
                if (triggerSet == 3)
                    trigger     = 83;
                if (triggerSet == 4)
                    trigger     = 51;
                if (triggerSet == 5)
                    trigger     = 01;
            } else if (energy.CompareTo("8TeV") == 0 || energy.CompareTo("5TeV") == 0 || energy.CompareTo("13TeV") == 0 ){
                if (triggerSet == 1)
                    trigger     = 52;
                if (triggerSet == 2)
                    trigger     = 81;
                if (triggerSet == 3)
                    trigger     = 53;
                if (triggerSet == 4)
                    trigger     = 82;
            }
        }

        //***************************************************************************************
        //********************** Start setting pi0 example bins *********************************
        //***************************************************************************************
        if (meson.CompareTo("Pi0") == 0){
            if (energy.CompareTo("900GeV") == 0) {
                if (directPhotonRunningOption.CompareTo("directPhoton") != 0){
                    if (mode == 1)              // PCM-Dalitz
                        return 4;
                    else if (mode == 2 || mode == 13)
                        return 4;
                    else if (mode == 4 || mode == 12 )
                        return 4;
                    else
                        return 5;
                } else {
                    if (mode == 2)
                        return 4;
                    else if (mode == 4)
                        return 7;
                    else
                        return 5;
                }
            } else if (energy.CompareTo("2.76TeV") == 0) {
                if (mode == 0){             // PCM-PCM
                    return 7;
                } else if (mode == 1){      // PCM-Dalitz
                    return 3;
                } else if (mode == 2 || mode == 13) {     // PCM-EMC
                    switch (trigger){
                        case 0:             // INT1 13g
                        case 1:             // INT1 13g
                            return 7;
                            break;
                        case 3:             // INT1 11a
                            return 5;
                            break;
                        case 10:            // INT7 13g
                        case 11:            // INT8 13g
                            return 7;
                            break;
                        case 51:            // EMC1
                            return 21;
                            break;
                        case 52:            // EMC7
                            return 14;
                            break;
                        case 85:            // EG2
                            return 16;
                            break;
                        case 83:            // EG1
                            return 21;
                            break;
                        default:
                            return 7;
                            break;
                    }
                } else if ( mode == 4 || mode == 12  ){    // EMC-EMC
                    switch (trigger){
                        case 0:             // INT1 13g
                        case 1:             // INT1 13g
                        case 3:             // INT1 11a
                        case 10:            // INT7 13g
                        case 11:            // INT8 13g
                            return 9;
                            break;
                        case 51:            // EMC1
                            return 23;
                            break;
                        case 52:            // EMC7
                            return 15;
                            break;
                        case 85:            // EG2
                            return 18;
                            break;
                        case 83:            // EG1
                            return 24;
                            break;
                        default:
                            return 7;
                            break;
                    }
                } else if ( mode == 10 ){
                    switch (trigger){
                        case 0:             // INT1 13g
                        case 1:             // INT1 13g
                        case 3:             // INT1 11a
                        case 10:            // INT7 13g
                        case 11:            // INT8 13g
                            return 22;
                            break;
                        case 51:            // EMC1
                            return 24;
                            break;
                        case 52:            // EMC7
                            return 22;
                            break;
                        case 85:            // EG2
                            return 26;
                            break;
                        case 83:            // EG1
                            return 26;
                            break;
                        default:
                            return 20;
                            break;
                    }

                } else {
                    return 7;
                }
            } else if (energy.CompareTo("5TeV") == 0 || energy.CompareTo("5TeV2017") == 0) {
                if ( mode == 1 )
                    return 5;
                if ( mode == 2 )
                    return 9;
                if ( mode == 4 )
                    return 24;
                if ( mode == 12 )
                    return 12;
                if ( mode == 13 )
                    return 9;
                else
                  return 14;
            } else if (energy.CompareTo("7TeV") == 0) {
                if ( mode == 0 )
                    return 4;
                else if ( mode == 1 )
                    return 5;
                else if ( mode == 3 )
                    return 4;
                else
                    return 14;
            } else if (energy.CompareTo("8TeV") == 0) {
                if (mode == 0){             // PCM- PCM
                    switch (trigger){
                        case 0:
                        case 1:
                        case 10:
                            return 7;
                        case 11:
                            return 3;       // INT triggers
                            break;
                        case 52:
                            return 33;
                        case 53:
                            return 33;      // EMC triggers
                            break;
                        case 81:
                            return 34;      // EGA triggers
                            break;
                        case 82:
                            return 40;      // EGA triggers
                            break;
                        default:
                            return 3;
                            break;
                    }
                } else if (mode == 2 || mode == 13){      // PCM-EMC
                    if (directPhotonRunningOption.CompareTo("directPhoton") == 0){
                        return 5;
                    } else if (directPhotonRunningOption.CompareTo("directPhotonTagging") == 0){
                        return 1;
                    } else {
                        switch (trigger){
                            case 0:
                            case 1:
                            case 10:
                            case 11:
                                return 3;       // INT triggers
                                break;
                            case 52:
                            case 53:
                                return 31;      // EMC triggers
                                break;
                            case 81:
                            case 82:
                                return 41;      // EGA triggers
                                break;
                            default:
                                return 3;
                                break;
                        }
                    }
                } else if (mode == 4 || mode == 12 ){      // EMC-EMC
                    if (directPhotonRunningOption.CompareTo("directPhoton") != 0 ){
                        switch (trigger){
                            case 0:
                            case 1:
                            case 10:
                            case 11:
                                return 8;      // INT triggers
                                break;
                            case 52:
                            case 53:
                                return 31;      // EMC triggers
                                break;
                            case 81:
                            case 82:
                                return 35;      // EGA triggers
                                break;
                            default:
                                return 13;
                                break;
                        }
                    } else {
                        switch (trigger){
                            case 0:
                            case 1:
                            case 10:
                            case 11:
                                return 7;      // INT triggers
                                break;
                            case 52:
                            case 53:
                                return 25;      // EMC triggers
                                break;
                            case 81:
                            case 82:
                                return 36;      // EGA triggers
                                break;
                            default:
                                return 13;
                                break;
                        }

                    }
                } else if (mode == 10){      // EMC-EMC
                    switch (trigger){
                        case 0:
                        case 1:
                        case 10:
                        case 11:
                            return 38;      // INT triggers
                            break;
                        case 52:
                        case 53:
                            return 45;      // EMC triggers
                            break;
                        case 81:
                        case 82:
                            return 37;      // EGA triggers
                            break;
                        default:
                            return 38;
                            break;
                    }
                } else {                    // other modes
                    return 3;
                }
            } else if (energy.CompareTo("13TeV") == 0 || energy.CompareTo("13TeVRBins") == 0  ) {
                if (mode == 0){
                    return 5;
                } else if ( mode == 1 ){
                    return 5;
                } else if ( mode == 2 || mode == 13 ){
                    switch (trigger){
                        case 83:
                            return 20;
                        default:
                            return 5;
                    }
                } else {
                    return 10;
                }
            } else if (energy.CompareTo("13TeVLowB") == 0) {
                if (mode == 0){
                    scaleFac        = 4.;
                    return 2;
                } else {
                return 2;
                }
            } else if( energy.CompareTo("pPb_5.023TeVCent") == 0 ) {
                if (mode == 0){
                    return 7;
                } else if (mode == 2 || mode == 13){
                    return 14;
                } else if (mode == 3){
                    return 12;
                } else if (mode == 4){
                    return 15;
                }
            } else if( energy.CompareTo("pPb_5.023TeV") == 0 ) {
                if (mode == 0){
                    return 7;
                } else if (mode == 1){
                    return 5;
                } else if (mode == 2 || mode == 13){
                    if (directPhotonRunningOption.CompareTo("directPhotonTagging") == 0){
                        return 10;
                    } else if (directPhotonRunningOption.CompareTo("directPhotonA") == 0){
                        return 10;
                    } else {
                        switch (trigger){
                            case 0:
                            case 1:
                            case 10:
                            case 11:
                                return 7;      // INT triggers
                                break;
                            case 51:
                            case 52:
                            case 53:
                                return 10;      // EMC triggers
                                break;
                            case 85:
                                return 10;
                                break;
                            case 81:
                            case 82:
                            case 83:
                                return 28;      // EGA triggers
                                break;
                            default:
                                return 7;
                                break;
                        }
                    }
                } else if (mode == 3){
                    if (centrality.CompareTo("0-100%")){
                        return 7;
                    } else {
                        switch (trigger){
                            case 0:
                            case 1:
                            case 10:
                            case 11:
                                return 7;      // INT triggers
                                break;
                            case 62:
                                return 17;      // PHOS triggers
                                break;
                            default:
                                return 7;
                                break;
                        }
                    }
                } else if (mode == 4 || mode == 12 ){
                    if (directPhotonRunningOption.CompareTo("directPhoton") == 0)
                        return 13;
                    if (centrality.CompareTo("0-100%") != 0){
                        return 15;
                    } else {
                        switch (trigger){
                            case 0:
                            case 1:
                            case 10:
                            case 11:
                                return 15;      // INT triggers
                                break;
                            case 51:
                            case 52:
                            case 53:
                                return 14;      // EMC triggers
                                break;
                            case 85:
                                return 30;
                                break;
                            case 81:
                            case 82:
                            case 83:
                                return 15;      // EGA triggers
                                break;
                            default:
                                return 15;
                                break;
                        }
                    }
                } else if (mode == 5){
                    return 25;
                } else if (mode == 6){
                    return 7;
                } else if (mode == 7){
                    return 6;
                } else if (mode == 10){
                    return 6;
                } else {
                    return 7;
                }
            } else if( energy.CompareTo("pPb_5.023TeVRun2") == 0 ) {
                if (mode == 0){
                    return 7;
                } else if (mode == 1){
                    return 5;
                } else if (mode == 2 || mode == 13){
                    if (directPhotonRunningOption.CompareTo("directPhotonTagging") == 0){
                        return 10;
                    } else {
                        switch (trigger){
                            case 0:
                            case 1:
                            case 10:
                            case 11:
                                return 7;      // INT triggers
                                break;
                            case 51:
                            case 52:
                            case 53:
                                return 20;      // EMC triggers
                                break;
                            case 85:
                                return 26;
                                break;
                            case 81:
                            case 82:
                            case 83:
                                return 32;      // EGA triggers
                                break;
                            default:
                                return 7;
                                break;
                        }
                    }
                } else if (mode == 4 || mode == 12 ){
                    if (directPhotonRunningOption.CompareTo("directPhoton") == 0)
                        return 13;
                    else
                        return 16;
                } else if (mode == 5){
                    return 45;
                } else if (mode == 6){
                    return 7;
                } else if (mode == 7){
                    return 6;
                } else if (mode == 10){
                    return 6;
                } else {
                    return 7;
                }
            } else if( energy.CompareTo("pPb_8TeV") == 0 ) {
                if (mode == 0){
                    return 7;
                } else if (mode == 1){
                    return 5;
                } else if (mode == 5){
                    return 25;
                } else if (mode == 6){
                    return 7;
                } else if (mode == 7){
                    return 6;
                } else if (mode == 10){
                    switch (trigger){
                        case 10:            // INT7 13g
                            return 15;
                            break;
                        case 52:            // EMC7
                            return 15;
                            break;
                        case 85:            // EG2
                            return 18;
                            break;
                        case 83:            // EG1
                            return 24;
                            break;
                        default:
                            return 7;
                            break;
                    }
                } else {
                    return 7;
                }
            } else if( energy.CompareTo("PbPb_2.76TeV") == 0) {
                if (mode == 0){
                    scaleFac    = 10;
                    return 4;
                } else if (mode == 1){
                    scaleFac    = 20;
                    return 3;
                } else if (mode == 2 || mode == 13){
                    scaleFac    = 5;
                    return 8;
                } else if (mode == 4 || mode == 12 ){
                    scaleFac    = 1.5;
                    return 14;
                } else
                    return 4;
            } else if( energy.CompareTo("PbPb_5.02TeV") == 0 || energy.CompareTo("PbPb_5TeV") == 0) {
                if (mode == 0){
                    scaleFac    = 10;
                    return 4;
                } else if (mode == 1){
                    scaleFac    = 20;
                    return 3;
                } else if (mode == 2 || mode == 3 || mode == 13){
                    scaleFac    = 5;
                    return 20;
                } else if (mode == 4 || mode == 12 ){
                    scaleFac    = 1.5;
                    return 23;
                } else
                    return 4;
            } else if( energy.CompareTo("XeXe_5.44TeV") == 0) {
                if (mode == 0){
                    scaleFac    = 2;
                    return 5;
                } else if (mode == 2){
                    if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("20-40%") || !centrality.CompareTo("0-40%") || !centrality.CompareTo("40-80%")){
                        scaleFac        = 1;
                        return 7;
                    } else {
                        scaleFac        = 2;
                        return 12;
                    }
                } else if (mode == 3){
                    if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("20-40%") || !centrality.CompareTo("0-40%") || !centrality.CompareTo("40-80%")){
                        scaleFac        = 1;
                        return 6;
                    } else {
                        scaleFac        = 2;
                        return 10;
                    }
                } else if (mode == 4){
                    if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("20-40%") || !centrality.CompareTo("0-40%") || !centrality.CompareTo("40-80%")){
                        scaleFac        = 1;
                        return 10;
                    } else {
                        scaleFac        = 1;
                        return 15;
                    }
                } else if (mode == 5){
                    if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("20-40%") || !centrality.CompareTo("0-40%") || !centrality.CompareTo("40-80%")){
                        scaleFac        = 1;
                        return 10;
                    } else {
                        scaleFac        = 1;
                        return 12;
                    }

                } else {
                    scaleFac    = 2;
                    return 10;
                }
            }
        //***************************************************************************************
        //********************** Start setting eta example bins *********************************
        //***************************************************************************************
        } else if( meson.Contains("Eta") && !meson.EqualTo("EtaPrime") ) {
            if (energy.CompareTo("900GeV") == 0) {
                if (mode == 2 || mode == 13)
                    return 2;
                else
                    return 1;
            } else if (energy.CompareTo("2.76TeV") == 0) {
                if (mode == 0){             // PCM-PCM
                    return 4;
                } else if (mode == 2 || mode == 13) {     // PCM-EMC
                    switch (trigger){
                        case 0:             // INT1 13g
                        case 1:             // INT1 13g
                        case 10:            // INT7 13g
                        case 11:            // INT8 13g
                            return 4;
                            break;
                        case 3:             // INT1 11a
                            return 4;
                            break;
                        case 51:            // EMC1
                            return 8;
                            break;
                        case 52:            // EMC7
                            return 5;
                            break;
                        case 85:            // EG2
                            return 7;
                            break;
                        case 83:            // EG1
                            return 9;
                            break;
                        default:
                            return 7;
                            break;
                    }
                } else if ( mode == 4 || mode == 12  ){    // EMC-EMC
                    switch (trigger){
                        case 3:             // INT1 11a
                            return 6;
                            break;
                        case 0:             // INT1 13g
                        case 1:             // INT1 13g
                        case 10:            // INT7 13g
                        case 11:            // INT8 13g
                            return 7;
                            break;
                        case 52:            // EMC7
                            return 5;
                            break;
                        case 51:            // EMC1
                            return 8;
                            break;
                        case 85:            // EG2
                            return 7;
                            break;
                        case 83:            // EG1
                            return 10;
                            break;
                        default:
                            return 4;
                            break;
                    }
                } else {
                    return 4;
                }
            } else if (energy.CompareTo("5TeV") == 0 || energy.CompareTo("5TeV2017") == 0) {
                if (mode == 1 )
                    return 4;
                if (mode == 2 )
                    return 8;
                if (mode == 3 )
                    return 4;
                else if (mode == 4 )
                    return 9;
                else if (mode == 12 )
                    return 7;
                else
                    return 3;
            } else if (energy.CompareTo("7TeV") == 0) {
                if (mode == 1 )
                    return 4;
                if (mode == 3 ){
                    return 2;
                } else if(mode == 40){
                    scaleFac        = 1.;
                    return 5;
                } else if(mode == 41){
                    scaleFac        = 1.;
                    return 8;
                } else if(mode == 42){
                    scaleFac        = 1.;
                    return 6;
                } else if(mode == 44){
                    scaleFac        = 1.;
                    return 10;
                } else if(mode == 45){
                    scaleFac        = 1.;
                    return 7;
                } else
                    return 6;
            } else if (energy.CompareTo("8TeV") == 0) {
                if (mode == 0){             // PCM- PCM
                    switch (trigger){
                        case 0:
                        case 1:
                        case 10:
                            return 3;       // INT triggers
                        case 11:
                            return 6;       // INT triggers
                            break;
                        case 52:
                            return 14;
                        case 53:
                            return 11;      // EMC triggers
                            break;
                        case 81:
                        case 82:
                            return 13;      // EGA triggers
                            break;
                        default:
                            return 6;
                            break;
                    }
                } else if (mode == 2 || mode == 13){      // PCM-EMC
                    switch (trigger){
                        case 0:
                        case 1:
                        case 10:
                        case 11:
                            return 7;       // INT triggers
                            break;
                        case 52:
                        case 53:
                            return 14;      // EMC triggers
                            break;
                        case 81:
                        case 82:
                            return 19;      // EGA triggers
                            break;
                        default:
                            return 6;
                            break;
                    }

                } else if (mode == 4 || mode == 12 ){      // EMC-EMC
                    switch (trigger){
                        case 0:
                        case 1:
                        case 10:
                        case 11:
                            return 7;      // INT triggers
                            break;
                        case 52:
                        case 53:
                            if(meson.CompareTo("Pi0EtaBinning") == 0) return 15;
                            return 14;      // EMC triggers
                            break;
                        case 81:
                        case 82:
                            if(meson.CompareTo("Pi0EtaBinning") == 0) return 18;
                            return 20;      // EGA triggers
                            break;
                        default:
                            return 9;
                            break;
                    }
                } else {                    // other modes
                    return 6;
                }
            } else if (energy.CompareTo("13TeV") == 0 || energy.CompareTo("13TeVRBins") == 0) {
                if (mode == 0){
                    return 2;
                } else if (mode == 2){
                    switch (trigger){
                        case 83:
                            return 20;
                        default:
                            return 7;
                    }
                } else if(mode == 40){
                    scaleFac        = 4.;
                    return 2;
                } else if(mode == 60){
                    // scaleFac        = 4.;
                    return 2;
                } else{
                    return 7;
                }
            } else if (energy.CompareTo("13TeVLowB") == 0) {
                if (mode == 0){
                    scaleFac        = 30.;
                    return 2;
                } else if (mode == 4 || mode == 12){
                    scaleFac        = 4.;
                    return 14;
                } else if (mode == 2 || mode == 5){
                    scaleFac        = 4.;
                    return 13;
                } else {
                    return 2;
                }
            } else if( energy.CompareTo("pPb_5.023TeVCent") == 0 ) {
                if (mode == 0){
                    // scaleFac    = 2;
                    return 6;
                } else if (mode == 2 || mode == 13){
                    return 6;
                } else if (mode == 3){
                    return 7;
                } else if (mode == 4 || mode == 12 ){
                    return 10;
              }
            } else if( energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVRun2") == 0  ) {
                if (mode == 0){
                    // scaleFac    = 2;
                    return 6;
                } else if (mode == 1){
                    return 4;
                } else if (mode == 2 || mode == 13){
                    switch (trigger){
                        case 0:
                        case 1:
                        case 10:
                        case 11:
                            if (meson.CompareTo("Eta") == 0) scaleFac = 4.0;
                            return 8;      // INT triggers
                            break;
                        case 51:
                        case 52:
                        case 53:
                            return 5;      // EMC triggers
                            break;
                        case 85:
                            return 12;
                            break;
                        case 81:
                        case 82:
                        case 83:
                            return 13;      // EGA triggers
                            break;
                        default:
                            return 6;
                            break;
                    }
                } else if (mode == 3){
                    if (centrality.CompareTo("0-100%") == 0){
                        return 11;
                    } else {
                        return 7;
                    }
                } else if (mode == 4 || mode == 12 ){
                    switch (trigger){
                        case 0:
                        case 1:
                        case 10:
                        case 11:
                            return 10;      // INT triggers
                            break;
                        case 51:
                        case 52:
                        case 53:
                            return 8;      // EMC triggers
                            break;
                        case 85:
                            return 16;
                            break;
                        case 81:
                        case 82:
                        case 83:
                            return 15;      // EGA triggers
                            break;
                        default:
                            return 10;
                            break;
                    }
                } else if (mode == 5){
                    return 12;
                } else {
                    return 6;
                }
            } else if( energy.CompareTo("pPb_8TeV") == 0  ) {
                if (mode == 0){
                    // scaleFac    = 2;
                    return 6;
                } else if (mode == 1){
                    return 4;
                } else if (mode == 2 || mode == 13){
                    switch (trigger){
                        case 0:
                        case 1:
                        case 10:
                        case 11:
                            if (meson.CompareTo("Eta") == 0) scaleFac = 4.0;
                            return 8;      // INT triggers
                            break;
                        case 51:
                        case 52:
                        case 53:
                            return 12;      // EMC triggers
                            break;
                        case 85:
                            return 17;
                            break;
                        case 81:
                        case 82:
                        case 83:
                            return 20;      // EGA triggers
                            break;
                        default:
                            return 6;
                            break;
                    }
                } else if (mode == 3){
                    return 11;
                } else if (mode == 4 || mode == 12 ){
                    return 10;
                } else if (mode == 5){
                    return 10;
                } else {
                    return 6;
                }
            } else if( energy.CompareTo("PbPb_2.76TeV") == 0) {
                if (mode == 0){
                    scaleFac    = 40;
                    return 4;
                } else if (mode == 2 || mode == 13){
                    scaleFac    = 2;
                    return 7;
                } else if (mode == 4 || mode == 12 ){
                    scaleFac    = 1;
                    return 9;
                } else
                    return 4;
            } else if( energy.CompareTo("PbPb_5.02TeV") == 0) {
                if (mode == 2 || mode == 3 || mode == 13){
                    return 15;
                }else if (mode == 4 || mode == 12 ){
                    return 16;
		}else if (mode == 0){
		  if(meson.CompareTo("Pi0EtaBinning") == 0){
		    scaleFac = 1;
		    return 4;
		  } else {
		    scaleFac = 40;
		    return 2;
		  }
                }else
                    return 1;
            } else if( energy.CompareTo("XeXe_5.44TeV") == 0) {
                if (mode == 0)
                    return 4;
                else if (mode == 2)
                    return 6;
                else
                    return 6;
            }
        //***************************************************************************************
        //********************** Start setting eta prime example bins ***************************
        //***************************************************************************************
        } else if (meson.CompareTo("EtaPrime") == 0) {
            switch( mode ) {
                case 0: // PCM-PCM
                    switch( trigger ) {
                        case 10: scaleFac = 8.;  return 3; // INT7
                        case 52: scaleFac = 1.;  return 1; // L0
                        case 83: scaleFac = 1.;  return 2; // EG1
                        case 85: scaleFac = 1.;  return 2; // EG2
                    } break;
                case 2: // PCM-EMC
                    switch( trigger ) {
                        case 10: scaleFac = 10.; return 7; // INT7 (maybe bin 4?)
                        case 52: scaleFac = 1.;  return 1; // L0
                        case 83: scaleFac = 4.;  return 4; // EG1
                        case 85: scaleFac = 4.;  return 6; // EG2
                    } break;
                case 3: // PCM-PHOS
                    switch( trigger ) {
                        case 10: scaleFac = 4.;  return 3; // INT7
                        case 62: scaleFac = 8.;  return 2; // PHI7
                    } break;
                case 4: // EMC-EMC
                    switch( trigger ) {
                        case 10: scaleFac = 4.;  return 2; // INT7
                        case 52: scaleFac = 1.;  return 1; // L0
                        case 83: scaleFac = 4.;  return 3; // EG1
                        case 85: scaleFac = 2.;  return 3; // EG2
                    } break;
                case 5: // PHOS-PHOS
                    switch( trigger ) {
                        case 10: scaleFac = 1.;  return 2; // INT7
                        case 62: scaleFac = 1.;  return 6; // PHI7 (maybe 1 or 3 are nice?)
                    } break;
                default: scaleFac = 1.; return 1;
            }
        //***************************************************************************************
        //********************** Start setting omega example bins *******************************
        //***************************************************************************************
        } else if (meson.Contains("Omega")) {
            if(energy.CompareTo("13TeV") == 0) {
                if(mode == 40){
                    scaleFac        = 4.;
                    return 4;
                }
            } else {
                if(mode == 40){
                    scaleFac        = 1.;
                    return 9;
                } else if(mode == 41){
                    scaleFac        = 2.;
                    return 7;
                } else if(mode == 42){
                    scaleFac        = 2.;
                    return 9;
                } else if(mode == 44){
                    scaleFac        = 1.;
                    return 11;
                } else if(mode == 45){
                    scaleFac        = 5.;
                    return 5;
                } else{
                    scaleFac        = 2.;
                    return 2;
                }
            }

        } else {
            cout << "Single example bin for meson \"" << meson << "\" not defined" << endl;
            return 1;
        }
        return 0; // in case of switch fall through
    }

    //*************************************************************************************************
    //******************** GetStartBin for general combination ****************************************
    //*************************************************************************************************
    Int_t GetStartBin(
        TString   meson,
        TString   energy,
        Int_t     mode,
        Int_t     specialTrigg  =-1,
        TString   centrality    = "",
        TString   minECut       = ""
    ){

        // Heavy meson fix
        if( mode>=100 ) mode -= 100;

        Int_t startPtBin = 0;
        //*************************************************************************************************
        //******************** Determine startbin for Pi0  ************************************************
        //*************************************************************************************************
        if (meson.CompareTo("Pi0")==0){
            if (energy.CompareTo("2.76TeV") == 0){
                if ( mode == 0 ){
                    startPtBin     = 1;
                } else if ( mode == 1 ){
                    startPtBin     = 3;
                } else if ( mode == 2 || mode == 13 ){
                    startPtBin     = 3;
                } else if ( mode == 4 || mode == 12 ){
                    startPtBin     = 6;
                } else if ( mode == 5){
                    startPtBin     = 4;
                } else if ( mode == 10){
                    startPtBin     = 21;
                } else if (mode == 20){
                    startPtBin     = 1;
                }
            } else if (energy.CompareTo("5TeV") == 0 || energy.CompareTo("5TeV2017") == 0){
                if( energy.Contains("2017")){
                    if ( mode == 0){
                        startPtBin = 1;
                    } else if ( mode == 1){
                        startPtBin = 1;
                    } else if ( mode == 2 || mode == 13 ){
                        startPtBin = 3;
                    } else if ( mode == 3){
                        startPtBin = 1;
                    } else if ( mode == 4){
                      if (specialTrigg == 2) startPtBin = 49;
                      else startPtBin = 6;
                    } else if ( mode == 10){
                        startPtBin = 27;
                    } else if ( mode == 12){
                        startPtBin = 5;
                    } else
                      startPtBin     = 7;
                } else {
                    if ( mode == 2 ){
                        if (specialTrigg == 1) startPtBin = 7;
                        else if (specialTrigg == 2) startPtBin = 7;
                        else if (specialTrigg == 3) startPtBin = 7;
                        else startPtBin = 6;
                    }
                    else if( mode == 4 ){
                        if (specialTrigg == 1) startPtBin = 20;
                        else if (specialTrigg == 2) startPtBin = 20;
                        else if (specialTrigg == 3) startPtBin = 20;
                        else startPtBin = 8;
                    }
                    else if( mode == 12 || mode == 13 )
                      startPtBin = 2;
                    else if ( mode == 20 )
                      startPtBin     = 1;
                    else
                      startPtBin     = 7;
                }
            } else if (energy.CompareTo("8TeV") == 0){
                if ( mode == 0 ){
                    startPtBin     = 1;
                } else if ( mode == 1 ){
                    startPtBin     = 1;
                } else if ( mode == 2 || mode == 13 ){
                    startPtBin     = 2;
                } else if ( mode == 4 || mode == 12 ){
                    startPtBin     = 7;
                } else if ( mode == 5){
                    startPtBin     = 4;
                } else if ( mode == 10){
                    startPtBin     = 28;
                } else if (mode == 20){
                    startPtBin     = 1;
                }
            } else if (energy.CompareTo("13TeV") == 0 || energy.CompareTo("13TeVRBins") == 0 ){
                if ( mode == 0 ){
                    startPtBin     = 1;
                    if (specialTrigg == 1)  startPtBin = 7;
                } else if ( mode == 1 ){
                    startPtBin     = 1;
                } else if ( mode == 2 || mode == 13 ){
                    startPtBin     = 3;
                } else if ( mode == 3){
                    startPtBin     = 10;
                } else if ( mode == 4 || mode == 12 ){
                    startPtBin     = 1;
                } else if ( mode == 5){
                    startPtBin     = 1;
                } else if ( mode == 10){
                    startPtBin     = 28;
                } else if (mode == 20){
                    startPtBin     = 1;
                }
            } else if (energy.CompareTo("13TeVLowB") == 0){
                if ( mode == 0 || mode == 2 ){
                    startPtBin     = 1;
                } else if ( mode == 4  || mode == 12 ){
                    startPtBin     = 2;
                }

            } else if (energy.CompareTo("pPb_5.023TeVCent") == 0 ){
                if ( mode == 0 ){
                    startPtBin     = 1;
                } else if ( mode == 2 || mode == 13 ){
                    startPtBin     = 4;
                } else if ( mode == 3 ){
                    startPtBin     = 1;
                } else if ( mode == 4 || mode == 12 ){
                    startPtBin     = 7;
                } else if ( mode == -5 ){
                    startPtBin     = 1;
                }
            } else if (energy.CompareTo("pPb_5.023TeV") == 0 ){
                if ( mode == 0 ){
                    startPtBin     = 1;
                } else if ( mode == 1 ){
                    startPtBin     = 1;
                } else if ( mode == 2 || mode == 13 ){
                    if ((centrality.Contains("0-100%") && !centrality.Contains("60-100%")) && !(specialTrigg == 1 || specialTrigg == 2 || specialTrigg == 3))
                        startPtBin     = 6;
                    else if (!centrality.CompareTo("0-100%") && specialTrigg == 1)
                        startPtBin     = 1;
                    else if (!centrality.CompareTo("0-100%") && specialTrigg == 2)
                        startPtBin     = 4;
                    else if (!centrality.CompareTo("0-100%") && specialTrigg == 3)
                        startPtBin     = 5;
                    else
                        startPtBin     = 4;
                } else if ( mode == 3 ){
                    if ((centrality.Contains("0-100%") && !centrality.Contains("60-100%")) && specialTrigg != 4)
                        startPtBin     = 3;
                    else if (specialTrigg == 4)
                        startPtBin     = 10;
                    else
                        startPtBin     = 1;
                } else if ( mode == 4 || mode == 12 ){
                    if ((centrality.Contains("0-100%") && !centrality.Contains("60-100%")) && !(specialTrigg == 1 || specialTrigg == 2 || specialTrigg == 3))
                        startPtBin     = 9;
                    else if (!centrality.CompareTo("0-100%") && specialTrigg == 1)
                        startPtBin     = 3;
                    else if (!centrality.CompareTo("0-100%") && specialTrigg == 2)
                        startPtBin     = 8;
                    else if (!centrality.CompareTo("0-100%") && specialTrigg == 3)
                        startPtBin     = 6;
                    else
                        startPtBin     = 7;
                } else if ( mode == 5){
                    startPtBin     = 7;
                } else if ( mode == -5){
                    startPtBin     = 1;
                } else if ( mode == 10){
                    startPtBin     = 1;
                } else if (mode == 20){
                    startPtBin     = 1;
                } else {
                    startPtBin     = 1;
                }
            } else if (energy.CompareTo("pPb_5.023TeVRun2") == 0 ){
                if ( mode == 0 ){
                  if(centrality.Contains("0-5%") || centrality.Contains("0-1%"))
                    startPtBin     = 2;
                  else
                    startPtBin     = 1;
                } else if ( mode == 1 ){
                    startPtBin     = 1;
                } else if ( mode == 2 || mode == 13 ){
                    if (specialTrigg == 1)
                        startPtBin     = 14;
                    else if (specialTrigg == 2)
                        startPtBin     = 24;
                    else if (specialTrigg == 3)
                        startPtBin     = 29;
                    else if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-10%") || !centrality.CompareTo("60-100%"))
                        startPtBin     = 3;
                    else if (!centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") )
                        startPtBin     = 4;
                    else
                        startPtBin     = 2;
                } else if ( mode == 3 ){
                    startPtBin     = 5;
                } else if ( mode == 4 || mode == 12 ){
                    startPtBin     = 7;
                } else if ( mode == 5){
                    startPtBin     = 10;
                } else if ( mode == 10){
                    startPtBin     = 1;
                } else if (mode == 20){
                    startPtBin     = 1;
                } else {
                    startPtBin     = 1;
                }
            } else if (energy.CompareTo("pPb_8TeV") == 0 ){
                if ( mode == 0 ){
                    startPtBin     = 1;
                } else if (mode == 10){
                    startPtBin     = 28;
                } else if (mode == 2){
                    startPtBin     = 5;
                } else {
                    startPtBin     = 1;
                }
            } else if (energy.CompareTo("PbPb_2.76TeV") == 0){
                if ( mode == 0 ){
                    startPtBin      = 1;
                } else if ( mode == 1 ){
                    startPtBin      = 2;
                } else if ( mode == 2 || mode == 13 ){
                    startPtBin      = 4;
                } else if ( mode == 3 ){
                    startPtBin      = 3;
                } else if ( mode == 4 || mode == 12 ){
                    cout << minECut << endl;
                    if (minECut.Atoi() != 3)
                        startPtBin      = 12;
                    else
                        startPtBin      = 6;
                } else if ( mode == 5){
                    startPtBin      = 4;
                } else if (mode == 20){
                    startPtBin      = 1;
                } else {
                    startPtBin      = 1;
                }
            } else if (energy.CompareTo("PbPb_5.02TeV") == 0){
                if ( mode == 0 ){
                    startPtBin      = 1;
                } else if ( mode == 1 ){
                    startPtBin      = 2;
                } else if ( mode == 2 || mode == 13 ){
                    startPtBin      = 4;
                    if (centrality.CompareTo("0-20%") == 0 || centrality.CompareTo("0-10%") == 0 || centrality.CompareTo("0-5%") == 0 || centrality.CompareTo("5-10%") == 0){
                      startPtBin = 10;
                    }
                    if (centrality.CompareTo("10-20%") == 0){
                      startPtBin = 10;
                    }
                } else if ( mode == 3 ){
                    startPtBin      = 3;
                } else if ( mode == 4 || mode == 12 ){
                  //cout << "centrality.Data() = " << centrality.Data() << endl;
                    cout << minECut << endl;
                    if (minECut.Atoi() != 3)
                        startPtBin      = 6;
                    else
                        startPtBin      = 6;
                    if (centrality.CompareTo("0-20%") == 0 || centrality.CompareTo("0-10%") == 0 || centrality.CompareTo("0-5%") == 0 || centrality.CompareTo("5-10%") == 0){
                      startPtBin = 14;
                    }
                    if (centrality.CompareTo("10-20%") == 0){
                      startPtBin = 9;
                    }
                } else if ( mode == 5){
                    startPtBin      = 4;
                } else if (mode == 20){
                    startPtBin      = 1;
                } else {
                    startPtBin      = 1;
                }
            } else if (energy.CompareTo("XeXe_5.44TeV") == 0){
                if ( mode == 0 ){
                    if (!centrality.CompareTo("0-90%") || !centrality.CompareTo("0-80%") )
                        startPtBin     = 3;
                    else if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-40%") || !centrality.CompareTo("20-40%")  || !centrality.CompareTo("40-80%") )
                        startPtBin     = 1;
                    else
                        startPtBin     = 2;
                } else if ( mode == 1 ){
                    startPtBin     = 2;
                } else if ( mode == 2 || mode == 13 ){
                    if (!centrality.CompareTo("0-90%") || !centrality.CompareTo("0-80%") )
                        startPtBin     = 6;
                    else if (!centrality.CompareTo("0-20%") )
                        startPtBin     = 3;
                    else if ( !centrality.CompareTo("0-40%") )
                        startPtBin     = 2;
                    else if ( !centrality.CompareTo("20-40%")  || !centrality.CompareTo("40-80%") )
                        startPtBin     = 1;
                    else
                        startPtBin     = 7;
                } else if ( mode == 3 ){
                    if (!centrality.CompareTo("0-90%") || !centrality.CompareTo("0-80%") )
                        startPtBin     = 4;
                    else if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-40%") || !centrality.CompareTo("40-80%") )
                        startPtBin     = 2;
                    else if (!centrality.CompareTo("20-40%")  )
                        startPtBin     = 1;
                    else
                        startPtBin     = 6;
                } else if ( mode == 4 || mode == 12 ){
                    if (!centrality.CompareTo("0-90%") || !centrality.CompareTo("0-80%") )
                        startPtBin     = 8;
                    else if (!centrality.CompareTo("0-20%") )
                        startPtBin     = 3;
                    else if ( !centrality.CompareTo("0-40%") )
                        startPtBin     = 2;
                    else if ( !centrality.CompareTo("20-40%")  || !centrality.CompareTo("40-80%") )
                        startPtBin     = 1;
                    else
                        startPtBin     = 9;
                } else if ( mode == 5){
                    if (!centrality.CompareTo("0-90%") || !centrality.CompareTo("0-80%") )
                        startPtBin     = 6;
                    else if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-40%") || !centrality.CompareTo("20-40%")  || !centrality.CompareTo("40-80%") )
                        startPtBin     = 2;
                    else
                        startPtBin     = 8;
                } else if (mode == 20){
                    startPtBin     = 2;
                }
            }
        //*************************************************************************************************
        //******************** Determine startbin for Eta  ************************************************
        //*************************************************************************************************
        } else if( meson.Contains("Eta") && !meson.EqualTo("EtaPrime") ){
            if (energy.CompareTo("2.76TeV") == 0){
                if ( mode == 0 ){
                    startPtBin     = 1;
                } else if ( mode == 1 ){
                    startPtBin     = 1;
                } else if ( mode == 2 || mode == 13 ){
                    startPtBin     = 2;
                } else if ( mode == 3 ){
                    startPtBin     = 2;
                } else if ( mode == 4 || mode == 12 ){
                    startPtBin     = 4;
                } else if ( mode == 5){
                    startPtBin     = 3;
                } else if (mode == 20){
                    startPtBin     = 1;
                }
            } else if (energy.CompareTo("5TeV") == 0 || energy.CompareTo("5TeV2017") == 0){
              if( energy.Contains("2017")){
                  if ( mode == 0){
                      startPtBin = 1;
                 } else if ( mode == 1 ){
                      startPtBin = 1;
                  } else if ( mode == 2 ){
                      startPtBin = 2;
                  } else if ( mode == 3 ){
                      startPtBin = 1;
                  } else if ( mode == 4 ){
                    if (specialTrigg == 2) startPtBin = 16;
                    else startPtBin = 2;
                  } else if ( mode == 13 ){
                      startPtBin = 1;
                  } else if ( mode == 12 ){
                      startPtBin = 4;
                  } else
                    startPtBin     = 7;
              }else{
                if ( mode == 0 )
                    startPtBin     = 1;
                else if ( mode == 2 ){
                    startPtBin = 5;
                } else if ( mode == 4 ){
                    if (specialTrigg == 1) startPtBin = 11;
                    else if (specialTrigg == 2) startPtBin = 11;
                    else if (specialTrigg == 3) startPtBin = 11;
                    else startPtBin = 7;
                } else if ( mode == 20 ){
                  startPtBin     = 1;
                }
              }
            } else if (energy.CompareTo("7TeV") == 0){
                if (mode == 40){
                    startPtBin     = 4;
                } else if (mode == 41){
                    startPtBin     = 7;
                } else if (mode == 42){
                    startPtBin     = 5;
                } else if (mode == 44){
                    startPtBin     = 8;
                } else if (mode == 45){
                    startPtBin     = 4;
                }
            } else if (energy.CompareTo("8TeV") == 0){
                if ( mode == 0 ){
                    startPtBin     = 1;
                } else if ( mode == 1 ){
                    startPtBin     = 1;
                } else if ( mode == 2 || mode == 13 ){
                    startPtBin     = 3;
                } else if ( mode == 3 ){
                    startPtBin     = 2;
                } else if ( mode == 4 || mode == 12 ){
                    startPtBin     = 5;
                } else if ( mode == 5){
                    startPtBin     = 3;
                } else if (mode == 20){
                    startPtBin     = 1;
                }
            } else if (energy.CompareTo("13TeV") == 0 || energy.CompareTo("13TeVRBins") == 0){
                if( mode==0)        startPtBin = 1;
                else if( mode==2 || mode==13 )
                    if(specialTrigg == 2 )                        startPtBin = 11;
                    else if( specialTrigg==4 || specialTrigg==5 ) startPtBin = 3;
                    else                                          startPtBin = 1;
                else if( mode==3  ) startPtBin = 1;
                else if( mode==4 || mode==12 ) startPtBin = 1;
                else if( mode==5  ) startPtBin = 1;
                else if( mode==40 ) startPtBin = 2;
                else if( mode==41 ) startPtBin = 6;
                else if( mode==42 ) startPtBin = 4;
                else if( mode==44 ) startPtBin = 9;
                else if( mode==45 ) startPtBin = 7;
            } else if (energy.CompareTo("13TeVLowB") == 0){
                if ( mode == 0 ){
                    startPtBin     = 2;
                } else if ( mode == 4 || mode == 5 || mode == 12){
                    startPtBin     = 13;
                } else if ( mode == 2 ){
                    startPtBin     = 13;
                }
            } else if (energy.CompareTo("pPb_5.023TeVCent") == 0){
                if ( mode == 0 ){
                    startPtBin      = 3;
                } else if ( mode == 2 || mode == 13 ){
                    startPtBin      = 5;
                } else if ( mode == 3 ){
                    startPtBin      = 5;
                } else if ( mode == 4 || mode == 12 ){
                    startPtBin      = 6;
                }
            } else if (energy.CompareTo("pPb_5.023TeV") == 0 ){
                if ( mode == 0 ){
                    startPtBin      = 3;
                } else if ( mode == 1 ){
                    startPtBin      = 3;
                } else if ( mode == 2 || mode == 13 ){
                    if (!centrality.CompareTo("0-100%") && !(specialTrigg == 1 || specialTrigg == 2 || specialTrigg == 3))
                        startPtBin     = 5;
                    else if (!centrality.CompareTo("0-100%") && specialTrigg == 1 )
                        startPtBin     = 1;
                    else if (!centrality.CompareTo("0-100%") && specialTrigg == 2 )
                        startPtBin     = 5;
                    else if (!centrality.CompareTo("0-100%") && specialTrigg == 3 )
                        startPtBin     = 3;
                    else
                        startPtBin     = 5;
                } else if ( mode == 3 ){
                    if ((centrality.Contains("0-100%") && !centrality.Contains("60-100%")) && specialTrigg != 4)
                        startPtBin     = 4;
                    else if (specialTrigg == 4)
                        startPtBin     = 8;
                    else
                        startPtBin     = 5;
                } else if ( mode == 4 || mode == 12 ){
                    if (energy.CompareTo("pPb_5.023TeVCent") == 0)
                        startPtBin     = 6;
                    else if (!centrality.CompareTo("0-100%") && !(specialTrigg == 1 || specialTrigg == 2 || specialTrigg == 3))
                        startPtBin     = 7;
                    else if (!centrality.CompareTo("0-100%") && specialTrigg == 1 )
                        startPtBin     = 3;
                    else if (!centrality.CompareTo("0-100%") && specialTrigg == 2 )
                        startPtBin     = 8;
                    else if (!centrality.CompareTo("0-100%") && specialTrigg == 3 )
                        startPtBin     = 6;
                    else
                        startPtBin     = 6;
                } else if ( mode == 5){
                    startPtBin      = 11;
                } else if ( mode == -5){
                    startPtBin      = 3;
                } else if ( mode == 10){
                    startPtBin      = 11;
                } else if (mode == 20){
                    startPtBin      = 3;
                }
            } else if (energy.CompareTo("pPb_5.023TeVRun2") == 0){
                if ( mode == 0 ){
                    if (!centrality.CompareTo("0-5%")||!centrality.CompareTo("5-10%")||!centrality.CompareTo("0-1%"))
                        startPtBin     = 3;
                    else if (!centrality.CompareTo("0-1%"))
                        startPtBin     = 4;
                    else
                        startPtBin     = 2;
                } else if ( mode == 1 ){
                    startPtBin     = 3;
                } else if ( mode == 2 || mode == 13 ){
                    startPtBin     = 3;
                } else if ( mode == 3 ){
                    if (!centrality.CompareTo("0-20%") )
                        startPtBin     = 4;
                    else if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%"))
                        startPtBin     = 5;
                    else
                        startPtBin     = 3;
                } else if ( mode == 4 || mode == 12 ){
                    startPtBin     = 7;
                } else if ( mode == 5){
                    startPtBin     = 4;
                } else if (mode == 20){
                    startPtBin     = 3;
                }
            } else if (energy.CompareTo("pPb_8TeV") == 0){
                if ( mode == 0 ){
                    startPtBin     = 3;
                } else if ( mode == 1 ){
                    startPtBin     = 3;
                } else if ( mode == 2 || mode == 13 ){
                    startPtBin     = 5;
                } else if ( mode == 3 ){
                    startPtBin     = 4;
                } else if ( mode == 4 || mode == 12 ){
                    startPtBin     = 7;
                } else if ( mode == 5){
                    startPtBin     = 5;
                } else if (mode == 20){
                    startPtBin     = 3;
                }
            } else if (energy.CompareTo("XeXe_5.44TeV") == 0){
                if ( mode == 0 ){
                    if (centrality.CompareTo("0-80%") == 0)
                        startPtBin     = 1;
                    else
                        startPtBin     = 1;
                } else if ( mode == 1 ){
                    startPtBin     = 2;
                } else if ( mode == 2 || mode == 13 ){
                    startPtBin     = 3;
                } else if ( mode == 3 ){
                    startPtBin     = 3;
                } else if ( mode == 4 || mode == 12 ){
                    startPtBin     = 1;
                } else if ( mode == 5){
                    startPtBin     = 5;
                } else if (mode == 20){
                    startPtBin     = 1;
                }
            }
        //*************************************************************************************************
        //******************** Determine startbin for Eta Prime *******************************************
        //*************************************************************************************************
        } else if (meson.EqualTo("EtaPrime") ){
            if(      energy.EqualTo("7TeV")  ) startPtBin = 1;
            else if( energy.EqualTo("13TeV") ) startPtBin = 1;
        //*************************************************************************************************
        //******************** Determine startbin for Omega  **********************************************
        //*************************************************************************************************
        } else if (meson.CompareTo("Omega") == 0){
            if (energy.CompareTo("7TeV") == 0){
                if (mode == 40){
                    startPtBin     = 4;
                } else if (mode == 41){
                    startPtBin     = 6;
                } else if (mode == 42){
                    startPtBin     = 5;
                } else if (mode == 44){
                    startPtBin     = 8;
                } else if (mode == 45){
                    startPtBin     = 4;
                }
            } else if(energy.CompareTo("13TeV") == 0){
                if (mode == 40){
                    startPtBin     = 2;
                }
            }
        //*************************************************************************************************
        //******************** Determine startbin for direct Photon ***************************************
        //*************************************************************************************************
        } else if (meson.Contains("directPhoton") ) {
            if (energy.CompareTo("pPb_5.023TeV")==0 || energy.CompareTo("pPb_5.023TeVCent") == 0  || energy.CompareTo("pPb_5.023TeVRun2") == 0 ){
                if (mode == 0){
                    if (!centrality.CompareTo("0-1%") || !centrality.CompareTo("0-2%"))
                        startPtBin      = 2;
                    else
                        startPtBin      = 1;
                } else if (mode == 2 && meson.CompareTo("directPhotonA") == 0 ){
                    if (!centrality.CompareTo("0-100%"))
                        startPtBin      = 8;
                    else
                        startPtBin      = 6;
                } else if (mode == 2 && meson.CompareTo("directPhotonTagging") == 0 ){
                    startPtBin      = 1;
                } else if (mode == 3){
                    startPtBin      = 5;
                } else if (mode == 4){
                    startPtBin      = 8;
                } else {
                    startPtBin      = 1;
                }

            } else if (energy.CompareTo("5TeV") == 0 || energy.CompareTo("5TeV2017") == 0){
                if( energy.Contains("2017")){
                    if ( mode == 0){
                        startPtBin = 1;
                    }
                }
            }
        } else if ( meson.CompareTo("Rho") == 0 || meson.CompareTo("K0Star") == 0){
            startPtBin     = 1;
        } else if (meson.CompareTo("CPion") == 0){
            startPtBin     = 2;
        } else if (meson.CompareTo("Proton") == 0){
            startPtBin     = 3;
        } else if (meson.CompareTo("Phi") == 0 ){
            startPtBin     = 4;
        } else if (meson.CompareTo("Lambda") == 0){
            startPtBin     = 5;
        } else if (meson.CompareTo("CKaon") == 0){
            startPtBin     = 7;
        }
        return startPtBin;
    }

    //*************************************************************************************************
    //******************** GetBinning for general combination *****************************************
    //*************************************************************************************************
    Int_t GetBinning(
        Double_t*   binning,
        Int_t       &binningMax,
        TString     meson           = "Pi0",
        TString     energy          = "2.76TeV",
        Int_t       mode            = 2,
        Int_t       SpecialTrigger  = -1,
        Bool_t      DCAcase         = kFALSE,
        TString     centrality      = ""
    ){
        Int_t maxNBins      = 0;
        binningMax          = 0;
        //*************************************************************************************************
        //*************************** Set binning for pi0 spectra *****************************************
        //*************************************************************************************************
        if (meson.CompareTo("Pi0")==0){
            if (energy.CompareTo("900GeV") == 0){
                if ( mode == 0 ){
                    maxNBins = 11;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0900GeVPt[i];
                    }
                }else if (mode == 2 || mode == 13){
                    maxNBins = 11;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0900GeVPCMEMCPt[i];
                    }
                }else if (mode == 4 || mode == 12 ){
                    maxNBins = 11;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0900GeVEMCPt[i];
                    }
                }
            } else if (energy.CompareTo("2.76TeV") == 0){
                if ( mode == 2 || mode == 13 ){
                    maxNBins = 25;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi02760GeVPtTrigFullPCMEMC[i];
                    }
                } else if ( mode == 0 ){
                    maxNBins = 19;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi02760GeVPt[i];
                    }
                } else if ( mode == 4 || mode == 12 ){
                    maxNBins = 26;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi02760GeVPtTrig13g[i];
                    }
                } else if ( mode == 10){
                    maxNBins = 32;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi02760GeVPtmEMC[i];

                    }
                } else if (mode == 20){
                    maxNBins = 33;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi02760GeVFullHaitaomEMC[i];
                    }
                }
            } else if (energy.CompareTo("5TeV") == 0 || energy.CompareTo("5TeV2017") == 0){
              if (DCAcase) {
                  if ( mode == 0 && energy.Contains("2017")){
                    maxNBins = 18;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi05TeV2017PtDCA[i];
                    }
                  }else{
                    maxNBins = 15;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi05TeVPtDCA[i];
                    }
                  }
              } else {
                  if ( mode == 0 ){
                      if(energy.Contains("2017")){
                          if(fNBinsPt>15 && fNBinsPt<22){
                              maxNBins = 22; // binnign for PbPb
                              for(Int_t i = 0; i < maxNBins+1; i++)
                                  binning[i] = fBinsPi05TeV2017PCMforPbPbPt[i]; //fBinsPi05TeVPt[i];
                          } else if(fNBinsPt>22 && fNBinsPt<26){
                              maxNBins = 26; // binnign as Hikari
                              for(Int_t i = 0; i < maxNBins+1; i++)
                                  binning[i] = fBinsPi05TeVPt[i];
                          } else if(fNBinsPt>28 && fNBinsPt<32){
                              maxNBins = 30; // binning for combination
                              for(Int_t i = 0; i < maxNBins+1; i++)
                                  binning[i] = fBinsPi05TeV2017PCMCombinationPt[i];
                          } else if(fNBinsPt>39 && fNBinsPt<45){
                              maxNBins = 43; // binning standalone
                              for(Int_t i = 0; i < maxNBins+1; i++)
                                  binning[i] = fBinsPi05TeV2017Pt[i];
                          } else if(fNBinsPt>60 && fNBinsPt<70){
                              maxNBins = 65; // finer binning
                              for(Int_t i = 0; i < maxNBins+1; i++)
                                  binning[i] = fBinsPi05TeV2017ExtraFinePt[i];
                          } else {
                              maxNBins = 43;
                              for(Int_t i = 0; i < maxNBins+1; i++)
                                  binning[i] = fBinsPi05TeV2017Pt[i];
                          }
                      } else {
                          maxNBins = 26;
                          for(Int_t i = 0; i < maxNBins+1; i++){
                              binning[i] = fBinsPi05TeVPt[i];
                          }
                      }
                  } else if ( mode == 1 ) {
                      maxNBins = 29;
                      for(Int_t i = 0; i < maxNBins+1; i++){
                          binning[i] = fBinsPi05TeV2017DalitzPt[i];
                      }
                  } else if ( mode == 2 || mode == 20){
                      if(energy.Contains("2017")){
                          maxNBins = 84;
                          for(Int_t i = 0; i < maxNBins+1; i++){
                              binning[i] = fBinsPi05TeV2017PCMEMCPt[i];
                          }
                      }else{
                          maxNBins = 34;
                          for(Int_t i = 0; i < maxNBins+1; i++){
                              binning[i] = fBinsPi05TeVPCMEMCPt[i];
                          }
                      }
                  } else if ( mode == 3 ){
                    maxNBins = 33;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi05TeV2017PtCombination[i];
                    }
                  } else if ( mode == 4 ){
                    if(energy.Contains("2017")){
                        maxNBins = 33;
                        for(Int_t i = 0; i < maxNBins+1; i++){
                            binning[i] = fBinsPi05TeV2017PtCombination[i];
                        }
                    }else{
                      if(SpecialTrigger == 1 || SpecialTrigger == 2 || SpecialTrigger == 3){
                          maxNBins = 50;
                          for(Int_t i = 0; i < maxNBins+1; i++){
                              binning[i] = fBinsPi05TeVPtEMCTrigger1[i];
                          }
                      } else {
                          maxNBins = 34;
                          for(Int_t i = 0; i < maxNBins+1; i++){
                              binning[i] = fBinsPi05TeVPtEMC[i];
                          }
                      }
                    }
                  } else if ( mode == 10 ){
                      maxNBins = 55;
                      for(Int_t i = 0; i < maxNBins+1; i++){
                          binning[i] = fBinsPi05TeV2017PtmEMC[i];
                      }
                  } else if ( mode == 12 ){
                    if(energy.Contains("2017")){
                      maxNBins = 32;
                      for(Int_t i = 0; i < maxNBins+1; i++){
                          binning[i] = fBinsPi05TeV2017PtDMC[i];
                      }
                    } else {
                      maxNBins = 14;
                      for(Int_t i = 0; i < maxNBins+1; i++){
                          binning[i] = fBinsPi05TeVPtDCal[i];
                      }
                    }
                  } else if ( mode == 13 ){
                    if(energy.Contains("2017")){
                        maxNBins = 80;
                        for(Int_t i = 0; i < maxNBins+1; i++){
                            binning[i] = fBinsPi05TeV2017PtPCMDCal[i];
                        }
                    }else{
                        maxNBins = 24;
                        for(Int_t i = 0; i < maxNBins+1; i++){
                            binning[i] = fBinsPi05TeVPtPCMDCal[i];
                        }
                    }
                  } else {
                      maxNBins = 26;
                      for(Int_t i = 0; i < maxNBins+1; i++){
                          binning[i] = fBinsPi05TeVPt[i];
                      }
                  }
              }
            } else if (energy.CompareTo("7TeV") == 0){
                if ( mode == 2 ){
                    maxNBins = 38;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi07TeVPCMEMCPt[i];
                    }
                } else if ( mode == 4 || mode == 12 ){
                    maxNBins = 38;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi07TeVEMCPt[i];
                    }
                } else if ( mode == 0 ){
                    maxNBins = 38;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi07TeVPt[i];
                    }
                }
            } else if (energy.CompareTo("8TeV") == 0){
                if ( mode == 2 ){
                    maxNBins = 46;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0Comb8TeVPt[i];
                    }
                } else if ( mode == 4 || mode == 12 ){
                    maxNBins = 41;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0Comb8TeVPt[i];
                    }
                } else if ( mode == 0 ){
                    maxNBins = 33;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi08TeVPt[i];
                    }
                } else if ( mode == 10 ){
                    maxNBins = 59;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi08TeVPtmEMC[i];
                    }
                } else if ( mode == 11 ){
                    maxNBins = 64;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi08TeVPtmEMCComb[i];
                    }
                }
            } else if (energy.CompareTo("13TeV") == 0  || energy.CompareTo("13TeVRBins") == 0 ){
                if (DCAcase==kTRUE) CopyVectorToArray( binningMax, fBinsPi013TeVPCMTrigINT7PtDCA, binning, 27 );
                // Copy binning according to cases
                else switch(mode) {
                    case 0:
                        switch(SpecialTrigger) {
                            case -1: maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVPCMTrigCombPt, binning, 103 ); break;
                            case 0:
			      if(energy.Contains("RBins")) maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVPCMTrigINT7RBinsPt, binning, 19 );
			      else                         maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVPCMTrigINT7Pt, binning, 84 );
			      break;
                           case 4:
                            case 5:
                                if( energy.Contains("RBins")) maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVPCMTrigINT7RBinsPt, binning, 19 );
                                else                          maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVPCMTrigINT7Pt, binning, 84 );
                                break;
                            case 1: maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVPCMTrigEMC7Pt, binning, 64 ); break;
                            case 2: maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVPCMTrigEG1Pt, binning, 103 ); break;
                            case 3: maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVPCMTrigEG2Pt, binning, 94 ); break;
                        }
                        break;
                    case 1:
                        switch(SpecialTrigger) {
                            case 0: maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVDalitzPt, binning, 20 ); break;
                            case 4:
                            case 5: maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVDalitzPt, binning, 20 ); break;
                            case 2: maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVDalitzPt, binning, 20 ); break;
                            case 3: maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVDalitzPt, binning, 20 ); break;
                        }
                        break;

                    case 2:
                        switch(SpecialTrigger) {
                            case 0:
                            case 4:
                            case 5: maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVPCMEMCTrigINT7Pt, binning, 123 ); break;
                            case 2: maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVPCMEMCTrigEG1Pt, binning, 117 ); break;
                            case 3: maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVPCMEMCTrigEG2Pt, binning, 112 ); break;
                        }
                        break;
                    case 3: //PCM-PHOS
                        switch(SpecialTrigger) {
                            default:
                                maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVPCMPHOSTrigINT7Pt, binning, 82 );
                                break;
                        }
                        break;
                    case 4:
                        switch(SpecialTrigger) {
                            case 0:
                            case 4:
                            case 5:  maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVEMCTrigINT7Pt, binning ); break;
                            case 1:  maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVEMCTrigEMC7Pt, binning ); break;
                            case 2:  maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVEMCTrigEG1Pt, binning ); break;
                            case 3:  maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVEMCTrigEG2Pt, binning ); break;
                            default: maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVEMCTrigCombPt, binning, 201 ); break;
                        }
                        break;
                    case 10: maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVPtmEMC, binning, 59 ); break;
                    default: maxNBins = CopyVectorToArray( binningMax, fBinsPi013TeVPCMEMCTrigINT7Pt, binning, 99 ); break;
                }
                // Check max bins and array size
                CheckBinSize(maxNBins,binningMax,kFALSE);
                cout<<"Get Binning(), Pi0 13TeV, binningMax: "<<binningMax<<"; maxNBins: "<<maxNBins<<endl;
                for( Int_t i=0; i<binningMax+1; i++ ) cout << binning[i] << ", ";
                cout << endl;
            } else if (energy.CompareTo("13TeVLowB") == 0){
                if ( mode == 0 ){
                    if (DCAcase == kTRUE)
                        maxNBins = 23;
                    else
                        maxNBins = 64;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        if (DCAcase == kTRUE) {
                            binning[i] = fBinsPi013TeVLowBPtDCA[i];
                        } else
                        binning[i] = fBinsPi013TeVLowBPt[i];
                    }
                } else if (mode == 4 || mode == 12){
                    maxNBins = 40;
                    for(Int_t i = 0; i < maxNBins+1; i++)
                        binning[i] = fBinsPi013TeVLowBEMCPt[i];
                } else if (mode == 2){
                    maxNBins = 52;
                    for(Int_t i = 0; i < maxNBins+1; i++)
                        binning[i] = fBinsPi013TeVLowBPCMEMCPt[i];
                } else if ( mode == 5 ){
                    maxNBins = 52;
                    for(Int_t i = 0; i < maxNBins+1; i++)
                        binning[i] = fBinsPi013TeVLowBPHOSPt[i];
                }

            } else if (energy.CompareTo("pPb_5.023TeVCent") == 0 ){
                if (mode == 0 ){ // PCM
                    binningMax  = 25;
                    maxNBins    = 25;
                    if (DCAcase) binningMax  = 16;
                    for(Int_t i = 0; i < binningMax+1; i++){
                    if (DCAcase)
                        binning[i] = fBinsPi0pPb5TeVPtDCA[i];
                    else
                        binning[i] = fBinsPi0pPb5TeVPCMCentPt[i];
                    }
                } else if ( mode == 2 || mode == 13 ) {
                    maxNBins    = 24;
                    binningMax  = 24;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsPi0pPb5TeVEMCCentPt[i];
                    }
                } else if ( mode == 4 || mode == 12  ) {
                    maxNBins    = 24;
                    binningMax  = 24;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsPi0pPb5TeVEMCCentPt[i];
                    }
                } else if ( mode == 3 || mode == 5 ) {
                    maxNBins    = 24;
                    binningMax  = 24;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsPi0pPb5TeVEMCCentPt[i];
                    }
                } else if ( mode == -5 ) {
                    maxNBins    = 23;
                    binningMax  = 23;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsPi0pPb5TeVPHOSCentPt[i];
                    }
                } else if ( mode == 20 ) {
                    maxNBins    = 27;
                    binningMax  = 27;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsPi0pPb5TeVCombCentPt[i];
                    }
                }
            } else if (energy.CompareTo("pPb_5.023TeV") == 0 ){
                if (mode == 0 ){ // PCM
                    maxNBins    = 31;
                    binningMax  = 39;
                    if (DCAcase) binningMax  = 16;
                    else if ( !(centrality.Contains("0-100%") && !centrality.Contains("60-100%"))){
                        binningMax  = 25;
                        maxNBins    = 25;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if (DCAcase){
                            if ( !centrality.CompareTo("0-100%"))
                                binning[i] = fBinsPi0pPb5TeVPtDCA[i];
                            else
                                binning[i] = fBinsPi0pPb5TeVPtDCACent[i];
                        }else if ( !centrality.CompareTo("0-100%"))
                            binning[i] = fBinsPi0pPb5TeVPt[i];
                        else
                            binning[i] = fBinsPi0pPb5TeVPCMCentPt[i];
                    }
                } else if (mode == 1){ // Dalitz
                    maxNBins = 22;
                    binningMax  = 22;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsPi0pPb5TeVDalitzPt[i];
                    }
                } else if ( mode == 2 || mode == 13 ) {
                    maxNBins    = 32;
                    binningMax  = 36;
                    if (SpecialTrigger == 3 )
                        binningMax  = 36;
                    else if ( SpecialTrigger == 1 )
                        binningMax  = 17;
                    else if (SpecialTrigger == 2 )
                        binningMax  = 40;
                    else if ( !(centrality.Contains("0-100%") && !centrality.Contains("60-100%"))){
                        binningMax  = 24;
                        maxNBins    = 24;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if ( SpecialTrigger == 1 )
                            binning[i]  = fBinsPi0pPb5TeVPtEMCTrigEMC7[i];
                        else if ( SpecialTrigger == 2 )
                            binning[i]  = fBinsPi0pPb5TeVPtEMCTrigEG2[i];
                        else if ( SpecialTrigger == 3  )
                            binning[i]  = fBinsPi0pPb5TeVPtEMCTrigEG1[i];
                        else if ( (centrality.Contains("0-100%") && !centrality.Contains("60-100%")))
                            binning[i] = fBinsPi0pPb5TeVEMCPt[i];
                        else
                            binning[i] = fBinsPi0pPb5TeVEMCCentPt[i];
                    }
                } else if ( mode == 4 || mode == 12  ) {
                    maxNBins    = 32;
                    binningMax  = 36;
                    if (SpecialTrigger == 3 )
                        binningMax  = 36;
                    else if ( SpecialTrigger == 1 )
                        binningMax  = 17;
                    else if (SpecialTrigger == 2 )
                        binningMax  = 40;
                    else if ( !(centrality.Contains("0-100%") && !centrality.Contains("60-100%"))){
                        binningMax  = 24;
                        maxNBins    = 24;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if ( SpecialTrigger == 1 )
                            binning[i]  = fBinsPi0pPb5TeVPtEMCTrigEMC7[i];
                        else if ( SpecialTrigger == 2 )
                            binning[i]  = fBinsPi0pPb5TeVPtEMCTrigEG2[i];
                        else if ( SpecialTrigger == 3  )
                            binning[i]  = fBinsPi0pPb5TeVPtEMCTrigEG1[i];
                        else if ( (centrality.Contains("0-100%") && !centrality.Contains("60-100%")))
                            binning[i] = fBinsPi0pPb5TeVEMCPt[i];
                        else
                            binning[i] = fBinsPi0pPb5TeVEMCCentPt[i];
                    }
                } else if ( mode == 3 || mode == 5 ) {
                    maxNBins    = 30;
                    binningMax  = 36;
                    if ( !(centrality.Contains("0-100%") && !centrality.Contains("60-100%"))){
                        binningMax  = 24;
                        maxNBins    = 24;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if ( (centrality.Contains("0-100%") && !centrality.Contains("60-100%")))
                            binning[i] = fBinsPi0pPb5TeVPHOSPt[i];
                        else
                            binning[i] = fBinsPi0pPb5TeVEMCCentPt[i];
                    }
                } else if ( mode == -5 ) {
                    maxNBins    = 30;
                    binningMax  = 30;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsPi0pPb5TeVPHOSDmitriPt[i];
                    }
                } else if ( mode == 6 ) {
                  maxNBins    = 22;
                  binningMax  = 22;
                  for(Int_t i = 0; i < binningMax+1; i++){
                    binning[i] = fBinsPi0pPb5TeVEMCDalitzPt[i];
                  }
                } else if ( mode == 10 ){
                  maxNBins    = 31;
                  binningMax  = 31;
                  for(Int_t i = 0; i < binningMax+1; i++){
                    binning[i] = fBinsPi0pPb5TeVmEMCPt[i];
                  }
                } else if (mode == 20){ //combined
                    maxNBins = 32;
                    if ( !centrality.CompareTo("20-40%")) maxNBins = 26;
                    else if ( !centrality.CompareTo("40-60%")) maxNBins = 26;
                    else if ( !centrality.CompareTo("60-100%")) maxNBins = 25;
                    else if ( centrality.CompareTo("0-100%")) maxNBins = 27;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        if ( !centrality.CompareTo("0-100%"))
                            binning[i] = fBinsPi0pPb5TeVEMCPt[i];
                        else
                            binning[i] = fBinsPi0pPb5TeVCombCentPt[i];
                    }
                } else if (mode == 21){ //combined R1 updated with PHOS triggers
                    maxNBins = 35;
                    for(Int_t i = 0; i < maxNBins+1; i++)
                        binning[i] = fBinsPi0pPb5TeVCombPt[i];
                }
            } else if (energy.CompareTo("pPb_5.023TeVRun2") == 0){
                if (mode == 0 ){ // PCM
                    if ( !centrality.CompareTo("0-1%") ){
                        maxNBins    = 70;
                        binningMax  = 70;
                    } else {
                        maxNBins    = 88;
                        binningMax  = 88;
                    }
                    if (DCAcase) binningMax  = 16;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if (DCAcase)
                            binning[i] = fBinsPi0pPb5TeVPtDCA[i];
                        else if ( !centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") || !centrality.CompareTo("0-20%") || !centrality.CompareTo("20-40%") || centrality.Contains("40-60") || !centrality.CompareTo("60-100%") )
                            binning[i]  = fBinsPi0pPb5TeVPCMR2CentPt[i];
                        else if ( !centrality.CompareTo("0-1%") || !centrality.CompareTo("0-2%") )
                            binning[i]  = fBinsPi0pPb5TeVPCMR2Cent2Pt[i];
                        else
                            binning[i] = fBinsPi0pPb5TeVPCMR2Pt[i];
                    }
                } else if (mode == 1){ // Dalitz
                    maxNBins = 22;
                    binningMax  = 22;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsPi0pPb5TeVDalitzPt[i];
                    }
                } else if ( mode == 2 || mode == 13 ) {
                    if ( !centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") || !centrality.CompareTo("60-100%")){
                        maxNBins    = 81;
                        binningMax  = 81;
                    } else {
                        maxNBins    = 82;
                        binningMax  = 84;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if ( !centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") || !centrality.CompareTo("60-100%") )
                            binning[i]  = fBinsPi0pPb5TeVPCMEMCR2CentPt[i];
                        else
                            binning[i]  = fBinsPi0pPb5TeVPCMEMCR2Pt[i];
                    }
                } else if ( mode == 3 ) {
                    maxNBins    = 88;
                    binningMax  = 88;
                    if (!centrality.CompareTo("0-20%")){
                        maxNBins    = 87;
                    } else if (!centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") ){
                        binningMax  = 80;
                        maxNBins    = 80;
                    } else if (!centrality.CompareTo("60-100%")){
                        binningMax  = 79;
                        maxNBins    = 79;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if (!centrality.CompareTo("60-100%"))
                            binning[i] = fBinsPi0pPb5TeVPCMPHOSR2PerPt[i];
                        else if (!centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") )
                            binning[i] = fBinsPi0pPb5TeVPCMPHOSR2CentPt[i];
                        else
                            binning[i] = fBinsPi0pPb5TeVPCMPHOSR2Pt[i];
                    }
                } else if ( mode == 4 || mode == 12  ) {
                    maxNBins    = 33;
                    binningMax  = 33;
                    for(Int_t i = 0; i < binningMax+1; i++){
                      binning[i] = fBinsPi0pPb5TeVEMCR2CentPt[i];
                    }
                } else if ( mode == 5 ) {
                    if ( !centrality.CompareTo("0-20%") ) {
                        maxNBins    = 87;
                        binningMax  = 87;
                    } else if ( !centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") ) {
                        maxNBins    = 84;
                        binningMax  = 84;
                    } else if ( !centrality.CompareTo("60-100%") ) {
                        maxNBins    = 85;
                        binningMax  = 85;
                    }else {
                        maxNBins    = 90;
                        binningMax  = 90;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if ( !centrality.CompareTo("0-20%") ) {
                            binning[i] = fBinsPi0pPb5TeVPHOSR2CentPt[i];
                        } else if ( !centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") ) {
                            binning[i] = fBinsPi0pPb5TeVPHOSR2SCentPt[i];
                        } else if ( !centrality.CompareTo("60-100%") ) {
                            binning[i] = fBinsPi0pPb5TeVPHOSR2TCentPt[i];
                        } else {
                            binning[i] = fBinsPi0pPb5TeVPHOSR2Pt[i];
                        }
                    }
                } else if ( mode == 6 ) {
                  maxNBins    = 22;
                  binningMax  = 22;
                  for(Int_t i = 0; i < binningMax+1; i++){
                    binning[i] = fBinsPi0pPb5TeVEMCDalitzPt[i];
                  }
                } else if ( mode == 10 ){
                  maxNBins    = 31;
                  binningMax  = 31;
                  for(Int_t i = 0; i < binningMax+1; i++){
                    binning[i] = fBinsPi0pPb5TeVmEMCPt[i];
                  }
                } else if (mode == 20){ //combined
                    maxNBins = 92;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0pPb5TeVR2Pt[i];
                    }
                }
            } else if (energy.CompareTo("pPb_8TeV") == 0 ){
                if (mode == 0 ){ // PCM
                    maxNBins = 31;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0pPb8TeVPt[i];
                    }
                } else if (mode == 1){ // Dalitz
                    maxNBins = 22;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0pPb8TeVDalitzPt[i];
                    }
                } else if ( mode == 2 || mode == 13 ) {
                    maxNBins = 32;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0pPb8TeVEMCPt[i];
                    }
                } else if ( mode == 4 || mode == 12  ) {
                    maxNBins = 32;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0pPb8TeVEMCPt[i];
                    }
                } else if ( mode == 3 || mode == 5 ) {
                    maxNBins = 32;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0pPb8TeVEMCPt[i];
                    }
                } else if (mode == 20){ //combined
                    maxNBins = 32;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0pPb8TeVEMCPt[i];
                    }
                  } else if (mode == 10){ //combined
                    maxNBins = 59;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0pPb8TeVPtmEMC[i];
                    }
                }
            } else if (energy.CompareTo("PbPb_2.76TeV") == 0 ){
                if (mode == 0 ){ // PCM
                    maxNBins = 15;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0PbPb2760GeVPtLHC11h[i];
                    }
                } else if ( mode == 2 || mode == 13 ) {
                    maxNBins = 24;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0PbPb2760GeVPtLHC11h[i];
                    }
                } else if ( mode == 4 || mode == 12  ) {
                    maxNBins = 24;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0PbPb2760GeVPtLHC11h[i];
                    }
                } else if (mode == 20){ //combined
                    maxNBins = 15;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0PbPb2760GeVPtLHC11h[i];
                    }
                }
            } else if (energy.CompareTo("PbPb_5.02TeV") == 0 ){
                if (mode == 0 ){ // PCM
                    maxNBins = 28;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0PbPb5TeVPt[i];
                    }
                } else if ( mode == 2 || mode == 3 || mode == 13 ) {
                    maxNBins = 33;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0PbPb5TeVPCMEMCPt[i];
                    }
                } else if ( mode == 4 || mode == 12  ) {
                    maxNBins = 33;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0PbPb5TeVEMCPt[i];
                    }
                } else if (mode == 20){ //combined
                    maxNBins = 28;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsPi0PbPb5TeVPt[i];
                    }
                }
            } else if (energy.CompareTo("XeXe_5.44TeV") == 0 ){
                binningMax = 25;
                if (mode == 0 ) { // PCM
                    maxNBins = 22;
                    if ( !centrality.CompareTo("0-80%") || !centrality.CompareTo("0-90%") ){
                        maxNBins = 22;
                        for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i] = fBinsPi0XeXe5440GeVPt[i];
                        }
                    } else if ( !centrality.CompareTo("0-20%") || !centrality.CompareTo("0-40%")  || !centrality.CompareTo("20-40%") ){
                        maxNBins    = 20;
                        binningMax  = 24;
                        for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i] = fBinsPi0XeXe5440GeVCentPCMPt[i];
                        }
                    } else if ( !centrality.CompareTo("40-80%") ){
                        maxNBins    = 18;
                        binningMax  = 21;
                        for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i] = fBinsPi0XeXe5440GeVPerPCMPt[i];
                        }
                    }
                } else if ( mode == 2 || mode == 13 ) {
                    maxNBins = 23;
                    if ( !centrality.CompareTo("0-80%") || !centrality.CompareTo("0-90%") ){
                        maxNBins = 23;
                        for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i] = fBinsPi0XeXe5440GeVPt[i];
                        }
                    } else if ( !centrality.CompareTo("0-20%") || !centrality.CompareTo("0-40%")  || !centrality.CompareTo("20-40%") ){
                        maxNBins    = 17;
                        if (!centrality.CompareTo("0-40%") ) maxNBins    = 18;
                        binningMax  = 20;
                        for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i] = fBinsPi0XeXe5440GeVCentPCMEMCPt[i];
                        }
                    } else if ( !centrality.CompareTo("40-80%") ){
                        maxNBins    = 15;
                        binningMax  = 17;
                        for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i] = fBinsPi0XeXe5440GeVPerPCMEMCPt[i];
                        }
                    }
                } else if ( mode == 3 ) {
                    if ( !centrality.CompareTo("0-80%") || !centrality.CompareTo("0-90%") ){
                        maxNBins = 21;
                        for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i] = fBinsPi0XeXe5440GeVPt[i];
                        }
                    } else if ( !centrality.CompareTo("0-20%") || !centrality.CompareTo("0-40%")  || !centrality.CompareTo("20-40%") ){
                        maxNBins    = 16;
                        if (!centrality.CompareTo("20-40%")) maxNBins    = 15;
                        binningMax  = 20;
                        for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i] = fBinsPi0XeXe5440GeVCentPCMPHOSPt[i];
                        }
                    } else if ( !centrality.CompareTo("40-80%") ){
                        maxNBins    = 15;
                        binningMax  = 18;
                        for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i] = fBinsPi0XeXe5440GeVPerPCMPHOSPt[i];
                        }
                    }
                } else if ( mode == 4 || mode == 12  ) {
                    if ( !centrality.CompareTo("0-80%") || !centrality.CompareTo("0-90%") ){
                        maxNBins = 24;
                        for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i] = fBinsPi0XeXe5440GeVPt[i];
                        }
                    } else if ( !centrality.CompareTo("0-20%") || !centrality.CompareTo("0-40%")  || !centrality.CompareTo("20-40%") ){
                        maxNBins    = 16;
                        if (!centrality.CompareTo("0-40%")) maxNBins    = 17;
                        binningMax  = 17;
                        for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i] = fBinsPi0XeXe5440GeVCentEMCPt[i];
                        }
                    } else if ( !centrality.CompareTo("40-80%") ){
                        maxNBins    = 14;
                        binningMax  = 14;
                        for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i] = fBinsPi0XeXe5440GeVPerEMCPt[i];
                        }
                    }
                } else if ( mode == 5  ) {
                    if ( !centrality.CompareTo("0-80%") || !centrality.CompareTo("0-90%") ){
                        maxNBins = 22;
                        for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i] = fBinsPi0XeXe5440GeVPt[i];
                        }
                    } else if ( !centrality.CompareTo("0-20%") || !centrality.CompareTo("0-40%")  || !centrality.CompareTo("20-40%") ){
                        maxNBins    = 18;
                        if (!centrality.CompareTo("20-40%")) maxNBins    = 17;
                        binningMax  = 20;
                        for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i] = fBinsPi0XeXe5440GeVCentPHOSPt[i];
                        }
                    } else if ( !centrality.CompareTo("40-80%") ){
                        maxNBins    = 15;
                        binningMax  = 17;
                        for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i] = fBinsPi0XeXe5440GeVPerPHOSPt[i];
                        }
                    }
                } else if (mode == 20){ //combined
                    maxNBins = 24;
                    if (centrality.CompareTo("0-80%") != 0 ) maxNBins = 23;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsPi0XeXe5440GeVPt[i];
                    }
                }
            }
        //*************************************************************************************************
        //*************************** Set binning for eta spectra *****************************************
        //*************************************************************************************************
        } else if (meson.CompareTo("Eta") == 0){
            if (energy.CompareTo("2.76TeV") == 0){
                maxNBins = 12;
                for(Int_t i = 0; i < maxNBins+1; i++){
                    binning[i] = fBinsEta2760GeVPtTrig11a[i];
                }
            } else if (energy.CompareTo("5TeV") == 0 || energy.CompareTo("5TeV2017") == 0){
                if (DCAcase){
                    if ( mode == 0 && energy.Contains("2017")){
                        maxNBins = 18;
                        for(Int_t i = 0; i < maxNBins+1; i++){
                            binning[i] = fBinsEta5TeV2017PtDCA[i];
                        }
                    } else {
                        maxNBins = 8;
                        for(Int_t i = 0; i < maxNBins+1; i++){
                            binning[i] = fBinsEta5TeVPtDCA[i];
                        }
                    }
              } else {
                 if ( mode == 0 ){
                    if(energy.Contains("2017")){
                        if(fNBinsPt<6){
                            maxNBins = 5; // binnning for PbPb
                            for(Int_t i = 0; i < maxNBins+1; i++)
                                binning[i] = fBinsEta5TeV2017PCMforPbPbPt[i];
                        } else if(fNBinsPt>6 && fNBinsPt<9){
                            maxNBins = 8; // binning as Hikari
                            for(Int_t i = 0; i < maxNBins+1; i++)
                                binning[i] = fBinsEta5TeVPt[i];
                        } else if(fNBinsPt>=9 && fNBinsPt<10){
                            maxNBins = 9; // binning for combination
                            for(Int_t i = 0; i < maxNBins+1; i++)
                                binning[i] = fBinsEta5TeV2017PCMCombinationPt[i];
                        } else {
                            maxNBins = 20;
                            for(Int_t i = 0; i < maxNBins+1; i++)
                                binning[i] = fBinsEta5TeV2017Pt[i];
                        }
                    } else {
                        maxNBins = 13;
                        for(Int_t i = 0; i < maxNBins+1; i++){
                            binning[i] = fBinsEta5TeVPt[i];
                        }
                    }
                } else if ( mode == 1  ){
                    maxNBins = 10;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEta5TeV2017DalitzPt[i];
                    }
                } else if ( mode == 2 ){
                  if(energy.Contains("2017")){
                    maxNBins = 17;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEta5TeV2017PCMEMCPt[i];
                    }
                  } else {
                    maxNBins = 18;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEta5TeVPCMEMCPt[i];
                    }
                  }
                } else if ( mode == 3 ){
                        maxNBins = 15;
                        for(Int_t i = 0; i < maxNBins+1; i++){
                            binning[i] = fBinsEta5TeV2017PtCombination[i];
                        }
                } else if ( mode == 4 ){
                  if(energy.Contains("2017")){
                    maxNBins = 15;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEta5TeV2017PtCombination[i];
                    }
                  } else {
                    if(SpecialTrigger == 1 || SpecialTrigger == 2 || SpecialTrigger == 3){
                      maxNBins = 29;
                      for(Int_t i = 0; i < maxNBins+1; i++){
                          binning[i] = fBinsEta5TeVEMCPtTrigger1[i];
                      }
                    } else {
                      maxNBins = 22;
                      for(Int_t i = 0; i < maxNBins+1; i++){
                          binning[i] = fBinsEta5TeVEMCPt[i];
                      }
                    }
                  }
                } else if ( mode == 20 ){
                    maxNBins = 18;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEta5TeVPCMEMCPt[i];
                    }
                }else if ( mode == 13  ){
                  if(energy.Contains("2017")){
                    maxNBins = 29;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEta5TeV2017PCMDCalPt[i];
                    }
                  } else {
                    maxNBins = 18;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEta5TeVPCMEMCPt[i];
                    }
                  }
                }else if ( mode == 12  ){
                    maxNBins = 11;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEta5TeV2017DMCPt[i];
                    }
                }
              }
            } else if (energy.CompareTo("7TeV") == 0){
                if ( mode == 2 || mode == 13 || mode == 4 || mode == 12  ){
                    maxNBins = 18;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEta7TeVPCMEMCPt[i];
                    }
                } else if ( mode == 0 ){
                    maxNBins = 16;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEta7TeVPt[i];
                    }
                } else if(mode == 40){
                    maxNBins = 12;
                    binningMax  = 12;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEtaPiPlPiMiPiZero7TevPtPCM[i];
                    }
                } else if(mode == 41){
                    maxNBins = 9;
                    binningMax  = 9;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEtaPiPlPiMiPiZero7TevPtPCMEMC[i];
                    }
                } else if(mode == 42){
                    maxNBins = 12;
                    binningMax  = 12;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEtaPiPlPiMiPiZero7TevPtPCMPHOS[i];
                    }
                } else if(mode == 44){
                    maxNBins = 11;
                    binningMax  = 11;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEtaPiPlPiMiPiZero7TevPtEMC[i];
                    }
                } else if(mode == 45){
                    maxNBins = 10;
                    binningMax  = 10;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEtaPiPlPiMiPiZero7TevPtPHOS[i];
                    }
                }
            } else if (energy.CompareTo("13TeV") == 0 || energy.CompareTo("13TeVRBins") == 0){
                switch(mode) {
                    case 0:
                        if( DCAcase==kTRUE ) maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVPCMTrigINT7PtDCA, binning, 9 );
                        else switch(SpecialTrigger) {
                            case -1: maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVPCMTrigCombPt, binning, 24 ); break;
                            case 0:
                             if(energy.Contains("RBins")) maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVPCMTrigINT7RBinsPt, binning, 19 );
                              else                         maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVPCMTrigINT7Pt, binning, 40 );
                              break;

                            case 4:
                            case 5: maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVPCMTrigINT7Pt, binning, 40 ); break;
                            case 1: maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVPCMTrigEMC7Pt, binning, 22 ); break;
                            case 2: maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVPCMTrigEG1Pt, binning, 22 ); break;
                            case 3: maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVPCMTrigEG2Pt, binning, 22 ); break;
                        }
                        break;
                    case 2:
                    case 13:
                        switch(SpecialTrigger) {
                            case 0:
                            case 4:
                            case 5: maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVPCMEMCTrigINT7Pt, binning, 38 ); break;
                            case 3: maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVPCMEMCTrigEG2Pt, binning, 78 ); break;
                            case 2: maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVPCMEMCTrigEG1Pt, binning ); break;
                            default: maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVPCMEMCTrigINT7Pt, binning ); break;
                        }
                        break;
                    case 3: maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVPCMPHOSTrigINT7Pt, binning ); break;
                    case 4:
                    case 12:
                        switch(SpecialTrigger) {
                            case 0:
                            case 4:
                            case 5: maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVEMCTrigINT7Pt, binning ); break;
                            default: maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVEMCTrigINT7Pt, binning ); break;
                        }
                        break;
                    case 5: maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVPHOSTrigINT7Pt, binning ); break;
                    case 40:
                    case 41:
                    case 42:
                    case 44:
                    case 45: maxNBins = CopyVectorToArray( binningMax, fBinsEtaPiPlPiMiPiZero13TevPtPCM, binning, 17 ); break;
                    default:
                        switch(SpecialTrigger) {
                            case 2:  maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVEMCTrigEG1Pt, binning ); break;
                            case 3:  maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVEMCTrigEG2Pt, binning ); break;
                            default: maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVEMCTrigCombPt, binning, 155 ); break;
                        }
                        break;
                }
            } else if (energy.CompareTo("13TeVLowB") == 0) {
                if( mode==0 || mode == 2 || mode == 4 || mode == 5 || mode == 12) maxNBins = CopyVectorToArray( binningMax, fBinsEta13TeVLowBPt, binning );
            } else if (energy.CompareTo("8TeV") == 0){
                if ( mode == 2 || mode == 13 || mode == 4 || mode == 12  ){
                    maxNBins = 26;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEtaComb8TeVPt[i];
                    }
                } else if ( mode == 0 ){
                    maxNBins = 19;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEta8TeVPt[i];
                    }
                }
            } else if (energy.CompareTo("pPb_5.023TeVCent") == 0){
                if (mode == 0){ // PCM
                    binningMax  = 14;
                    maxNBins    = 12;
                    if (DCAcase) binningMax  = 16;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if (DCAcase)
                            binning[i] = fBinsEtapPb5TeVPtDCA[i];
                        else
                            binning[i] = fBinsEtapPb5TeVCentPt[i];
                    }
                } else if (mode == 2 || mode == 13 ){ // PCM-EMC, PCM-DMC
                    maxNBins    = 17;
                    binningMax  = 17;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsEtapPb5TeVCentPt[i];
                    }
                } else if (mode == 3 ){ // PCM-PHOS
                    maxNBins    = 11;
                    binningMax  = 17;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsEtapPb5TeVCentPt[i];
                    }
                } else if (mode == 4 || mode == 12 ){ // EMC, DMC
                    maxNBins    = 13;
                    binningMax  = 17;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsEtapPb5TeVCentPt[i];
                    }
                } else if (mode == -5 || mode == 5){ // PHOS
                    maxNBins    = 9;
                    binningMax  = 10;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsEtapPb5TeVPHOSCentPt[i];
                    }
                } else if (mode == 20 ){ // Comb
                    maxNBins    = 17;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEtapPb5TeVCentPt[i];
                    }
                }
            } else if (energy.CompareTo("pPb_5.023TeV") == 0){
                if (mode == 0){ // PCM
                    maxNBins    = 16;
                    binningMax  = 22;
                    if (DCAcase) binningMax  = 16;
                    else if (!(centrality.Contains("0-100%") && !centrality.Contains("60-100%"))){
                        binningMax  = 14;
                        maxNBins    = 12;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if (DCAcase)
                            binning[i] = fBinsEtapPb5TeVPtDCA[i];
                        else if ((centrality.Contains("0-100%") && !centrality.Contains("60-100%")))
                            binning[i] = fBinsEtapPb5TeVPt[i];
                        else
                            binning[i] = fBinsEtapPb5TeVCentPt[i];
                    }
                } else if (mode == 1){ // PCM-Dalitz
                    maxNBins    = 9;
                    binningMax  = 9;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i]      = fBinsEtapPb5TeVDalitzPt[i];
                    }
                } else if (mode == 2 || mode == 13 ){ // PCM-EMC, PCM-DMC
                    maxNBins    = 18;
                    binningMax  = 21;
                    if (SpecialTrigger == 3)
                        binningMax  = 20;
                    else if (SpecialTrigger == 2)
                        binningMax  = 22;
                    else if (SpecialTrigger == 1)
                        binningMax  = 11;
                    else if (!(centrality.Contains("0-100%") && !centrality.Contains("60-100%"))){
                        binningMax  = 14;
                        maxNBins    = 12;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if (SpecialTrigger == 1)
                            binning[i] = fBinsEtapPb5TeVPtEMCTrigEMC7[i];
                        else if (SpecialTrigger == 2)
                            binning[i] = fBinsEtapPb5TeVPtEMCTrigEG2[i];
                        else if (SpecialTrigger == 3)
                            binning[i] = fBinsEtapPb5TeVPtEMCTrigEG1[i];
                        else if ((centrality.Contains("0-100%") && !centrality.Contains("60-100%")))
                            binning[i] = fBinsEtapPb5TeVEMCPt[i];
                        else
                            binning[i] = fBinsEtapPb5TeVCentPt[i];
                    }
                } else if (mode == 3 ){ // PCM-PHOS
                    maxNBins    = 14;
                    binningMax  = 20;
                    if (!(centrality.Contains("0-100%") && !centrality.Contains("60-100%")) && SpecialTrigger != 4){
                        binningMax  = 14;
                        maxNBins    = 12;
                    } else  if (centrality.CompareTo("0-100%") && SpecialTrigger == 4){
                        binningMax  = 12;
                        maxNBins    = 12;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if ((centrality.Contains("0-100%") && !centrality.Contains("60-100%")))
                            binning[i] = fBinsEtapPb5TeVPCMPHOSPt[i];
                        else
                            binning[i] = fBinsEtapPb5TeVCentPt[i];
                    }
                } else if (mode == 4 || mode == 12 ){ // EMC, DMC
                    maxNBins    = 19;
                    binningMax  = 21;
                    if (SpecialTrigger == 3)
                        binningMax  = 20;
                    else if (SpecialTrigger == 2)
                        binningMax  = 22;
                    else if (SpecialTrigger == 1)
                        binningMax  = 11;
                    if (!(centrality.Contains("0-100%") && !centrality.Contains("60-100%"))){
                        binningMax  = 14;
                        maxNBins    = 13;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if (SpecialTrigger == 1)
                            binning[i] = fBinsEtapPb5TeVPtEMCTrigEMC7[i];
                        else if (SpecialTrigger == 2)
                            binning[i] = fBinsEtapPb5TeVPtEMCTrigEG2[i];
                        else if (SpecialTrigger == 3)
                            binning[i] = fBinsEtapPb5TeVPtEMCTrigEG1[i];
                        if ((centrality.Contains("0-100%") && !centrality.Contains("60-100%")))
                            binning[i] = fBinsEtapPb5TeVEMCPt[i];
                        else
                            binning[i] = fBinsEtapPb5TeVCentPt[i];
                    }
                } else if (mode == 5 ){ // PHOS
                    maxNBins    = 14;
                    binningMax  = 19;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsEtapPb5TeVPHOSPt[i];
                    }
                } else if (mode == -5 ){ // PHOS
                    maxNBins    = 12;
                    binningMax  = 12;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsEtapPb5TeVPHOSDmitriPt[i];
                    }
                } else if (mode == 20 ){
                    maxNBins = 19;
                    if (!centrality.CompareTo("20-40%") || !centrality.CompareTo("40-60%") || !centrality.CompareTo("0-20%"))
                        maxNBins    = 15;
                    else if (!centrality.CompareTo("60-100%") )
                        maxNBins    = 14;
                    else if (centrality.CompareTo("0-100%"))
                        maxNBins    = 16;

                    for(Int_t i = 0; i < maxNBins+1; i++){
                        if (!centrality.CompareTo("0-100%"))
                            binning[i] = fBinsEtapPb5TeVEMCPt[i];
                        else
                            binning[i] = fBinsEtapPb5TeVCentPt[i];

                    }
                } else if (mode == 21 ){
                    maxNBins = 21;
                    for(Int_t i = 0; i < maxNBins+1; i++)
                        binning[i] = fBinsEtapPb5TeVCombPt[i];
                }
            } else if (energy.CompareTo("pPb_5.023TeVRun2") == 0){
                if (mode == 0){ // PCM
                    maxNBins    = 35;
                    binningMax  = 35;
                    if (DCAcase)
                        binningMax  = 16;
                    else if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-10%") || !centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") || !centrality.CompareTo("0-20%") || !centrality.CompareTo("20-40%") || centrality.Contains("40-60") || !centrality.CompareTo("60-100%"))
                        binningMax  = 29;

                    for(Int_t i = 0; i < binningMax+1; i++){
                        if (DCAcase)
                            binning[i] = fBinsEtapPb5TeVPtDCA[i];
                        else if ( !centrality.CompareTo("0-10%") || !centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") || !centrality.CompareTo("0-20%") || !centrality.CompareTo("20-40%") || centrality.Contains("40-60") || !centrality.CompareTo("60-100%"))
                            binning[i] = fBinsEtapPb5TeVPCMR2CentPt[i];
                        else if (!centrality.CompareTo("0-1%") || !centrality.CompareTo("0-2%"))
                            binning[i] = fBinsEtapPb5TeVPCMR2Cent2Pt[i];
                        else
                            binning[i] = fBinsEtapPb5TeVPCMR2Pt[i];
                    }
                } else if (mode == 1){ // PCM-Dalitz
                    maxNBins    = 9;
                    binningMax  = 9;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i]      = fBinsEtapPb5TeVDalitzPt[i];
                    }
                } else if (mode == 2 || mode == 13 ){ // PCM-EMC, PCM-DMC

                    if (SpecialTrigger  > 0 ){
                        binningMax  = 26;
                        maxNBins    = 32;
                    // different cents
                    } else if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-10%")){
                        binningMax  = 29;
                        maxNBins    = 29;
                    } else if (!centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") || !centrality.CompareTo("60-100%")){
                        binningMax  = 28;
                        maxNBins    = 28;
                    } else{
                        binningMax  = 32;
                        maxNBins    = 32;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if (SpecialTrigger  > 0 )
                            binning[i] = fBinsEtapPb5TeVPtEMCTrig[i];
                        else if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-10%") )
                            binning[i] = fBinsEtapPb5TeVPCMEMCR2CentPt[i];
                        else if (!centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") || !centrality.CompareTo("60-100%"))
                            binning[i] = fBinsEtapPb5TeVPCMEMCR2SCentPt[i];
                        else
                            binning[i] = fBinsEtapPb5TeVPCMEMCR2Pt[i];
                    }
                } else if (mode == 3 ){ // PCM-PHOS
                    maxNBins    = 28;
                    binningMax  = 28;
                    if (!centrality.CompareTo("0-20%")){
                        maxNBins = 26;
                    } else if (!centrality.CompareTo("60-100%")){
                        maxNBins    = 26;
                        binningMax  = 26;
                    } else if (!centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%")) {
                        maxNBins    = 24;
                        binningMax  = 24;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if (!centrality.CompareTo("60-100%"))
                            binning[i] = fBinsEtapPb5TeVPCMPHOSR2PerPt[i];
                        else if (!centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%"))
                            binning[i] = fBinsEtapPb5TeVPCMPHOSR2SCentPt[i];
                        else
                            binning[i] = fBinsEtapPb5TeVPCMPHOSR2Pt[i];
                    }
                } else if (mode == 4 || mode == 12 ){ // EMC, DMC
                    maxNBins    = 21;
                    binningMax  = 21;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsEtapPb5TeVEMCR2CentPt[i];
                    }
                } else if (mode == 5 ){ // PHOS
                    if ( !centrality.CompareTo("0-20%") || !centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") || !centrality.CompareTo("60-100%") ) {
                        maxNBins    = 20;
                        binningMax  = 20;
                    } else {
                        maxNBins    = 30;
                        binningMax  = 30;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if ( !centrality.CompareTo("0-20%") || !centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") || !centrality.CompareTo("60-100%") ) {
                            binning[i] = fBinsEtapPb5TeVPHOSR2CentPt[i];
                        } else {
                            binning[i] = fBinsEtapPb5TeVPHOSR2Pt[i];
                        }
                    }
                } else if (mode == 20 ){
                    maxNBins = 39;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEtapPb5TeVR2Pt[i];
                    }
                }
            } else if (energy.CompareTo("pPb_8TeV") == 0){
                if (mode == 0){
                    maxNBins = 16;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEtapPb8TeVPt[i];
                    }
                } else if (mode == 2 || mode == 13 ){
                    maxNBins = 18;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEtapPb8TeVEMCPt[i];
                    }
                } else if (mode == 3 ){
                    maxNBins = 14;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEtapPb8TeVPCMPHOSPt[i];
                    }
                } else if (mode == 4 || mode == 12 ){
                    maxNBins = 19;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEtapPb8TeVEMCPt[i];
                    }
                } else if (mode == 20 ){
                    maxNBins = 19;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsEtapPb8TeVEMCPt[i];
                    }
                }
            } else if (energy.CompareTo("XeXe_5.44TeV") == 0 ){
                binningMax  = 7;
                maxNBins    = 7;
                for(Int_t i = 0; i < binningMax+1; i++){
                    binning[i] = fBinsEtaXeXe5440GeVPt[i];
                }
            }
        } else if (meson.CompareTo("EtaPrime") == 0) {
            if     (energy.EqualTo("7TeV")) {
                maxNBins = CopyVectorToArray( binningMax, fBinsEtaPrime7TeVPt, binning );
            } else if(energy.EqualTo("13TeV")) {
                switch(mode) {
                    case 0: // PCM-PCM
                        switch( SpecialTrigger ) {
                            case 0: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_PCM_INT7_Pt,binning); break; // 10
                            case 1: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_PCM_L0_Pt,  binning); break; // 52
                            case 2: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_PCM_EG1_Pt, binning); break; // 83
                            case 3: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_PCM_EG2_Pt, binning); break; // 85
                        } break;
                    case 2: // PCM-EMC
                        switch( SpecialTrigger ) {
                            case 0: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_PCMEMC_INT7_Pt,binning); break; // 10
                            case 1: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_PCMEMC_L0_Pt,  binning); break; // 52
                            case 2: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_PCMEMC_EG1_Pt, binning); break; // 83
                            case 3: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_PCMEMC_EG2_Pt, binning); break; // 85
                        } break;
                    case 3: // PCM-PHOS
                        switch( SpecialTrigger ) {
                            case 0: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_PCMPHOS_MinBias_Pt,binning); break; // 10
                            case 6: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_PCMPHOS_VZERO_Pt,  binning); break; // 62
                        } break;
                    case 4: // EMC-EMC
                        switch( SpecialTrigger ) {
                            case 0: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_EMC_INT7_Pt,binning); break; // 10
                            case 1: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_EMC_L0_Pt,  binning); break; // 52
                            case 2: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_EMC_EG1_Pt, binning); break; // 83
                            case 3: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_EMC_EG2_Pt, binning); break; // 85
                        } break;
                    case 5: // PHOS-PHOS
                        switch( SpecialTrigger ) {
                            case 0: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_PHOS_MinBias_Pt,binning); break; // 10
                            case 6: maxNBins = CopyVectorToArray(binningMax,fBinsEtaPrime13TeV_PHOS_VZERO_Pt,  binning); break; // 62
                        } break;
                    case 60: maxNBins = CopyVectorToArray( binningMax, fBinsEtaPrime13TeV_PCM_INT7_Pt,     binning ); break;
                    case 61: maxNBins = CopyVectorToArray( binningMax, fBinsEtaPrime13TeV_PCMEMC_INT7_Pt,  binning ); break;
                    case 63: maxNBins = CopyVectorToArray( binningMax, fBinsEtaPrime13TeV_PCMPHOS_MinBias_Pt, binning ); break;
                    case 64: maxNBins = CopyVectorToArray( binningMax, fBinsEtaPrime13TeV_EMC_INT7_Pt,     binning ); break;
                    case 65: maxNBins = CopyVectorToArray( binningMax, fBinsEtaPrime13TeV_PHOS_MinBias_Pt,    binning ); break;
                }
            }
        } else if (meson.Contains("Omega")){
            if (energy.CompareTo("7TeV") == 0){
                if(mode == 40){
                    maxNBins    = 13;
                    binningMax  = 13;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsOmegaPiPlPiMiPiZero7TevPtPCM[i];
                    }
                } else if(mode == 41){
                    maxNBins    = 11;
                    binningMax  = 11;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsOmegaPiPlPiMiPiZero7TevPtPCMEMC[i];
                    }
                } else if(mode == 42){
                    maxNBins    = 12;
                    binningMax  = 12;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsOmegaPiPlPiMiPiZero7TevPtPCMPHOS[i];
                    }
                } else if(mode == 44){
                    maxNBins    = 12;
                    binningMax  = 12;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsOmegaPiPlPiMiPiZero7TevPtEMC[i];
                    }
                } else if(mode == 45){
                    maxNBins    = 9;
                    binningMax  = 9;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i] = fBinsOmegaPiPlPiMiPiZero7TevPtPHOS[i];
                    }
                }
            }
        } else if (meson.CompareTo("Gamma") == 0){
            if (energy.CompareTo("2.76TeV") == 0){
                if (mode == 0 || mode == 2){
                    maxNBins = 18;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i]  = fBinsDirGamma2760GeVPt[i];
                    }
                } else if (mode == 4){
                    maxNBins = 19;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsDirGamma2760GeVPt[i];
                    }
                } else if (mode == 20) {
                    maxNBins = 19;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsDirGamma2760GeVPt[i];
                    }
                }
            } else if (energy.CompareTo("8TeV") == 0){
                if (mode == 0){
                    maxNBins = 23;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i]  = fBinsDirGamma8TeVPt[i];
                    }
                }
            } else if (energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0 || energy.CompareTo("pPb_5.023TeVRun2") == 0 ){
                if (mode == 0){
                    if(energy.CompareTo("pPb_5.023TeVRun2") == 0){
                        if (!centrality.CompareTo("5-10%") ){
                            maxNBins    = 88;
                            binningMax  = 88;
                        }else if (!centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") || !centrality.CompareTo("0-20%") || !centrality.CompareTo("20-40%") || !centrality.CompareTo("40-60%") || !centrality.CompareTo("60-100%")){
                            maxNBins    = 45;
                            binningMax  = 45;
                        } else if ( !centrality.CompareTo("0-2%") || !centrality.CompareTo("0-1%") ){
                            maxNBins    = 28;
                            binningMax  = 28;
                        } else {
                            maxNBins    = 88;
                            binningMax  = 88;
                        }
                        for(Int_t i = 0; i < binningMax+1; i++){
                            if ( !centrality.CompareTo("5-10%"))
                                binning[i]  = fBinsDirGammapPb5TeVR2CentPCMPt[i];
                            else if (!centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") || !centrality.CompareTo("0-20%") || !centrality.CompareTo("20-40%") || !centrality.CompareTo("40-60%") || !centrality.CompareTo("60-100%"))
                                binning[i]  = fBinsDirGammapPb5TeVR2Cent2PCMPt[i];
                            else if ( !centrality.CompareTo("0-1%") || !centrality.CompareTo("0-2%") )
                                binning[i]  = fBinsDirGammapPb5TeVR2Cent3PCMPt[i];
                            else
                                binning[i] = fBinsDirGammapPb5TeVR2PCMPt[i];
                            // if ( !centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") || !centrality.CompareTo("60-100%") )
                            //     binning[i]  = fBinsDirGammapPb5TeVR2CentPCMPt[i];
                            // else if ( !centrality.CompareTo("0-1%") || !centrality.CompareTo("0-2%") )
                            //     binning[i]  = fBinsDirGammapPb5TeVR2Cent2PCMPt[i];
                            // else
                            //     binning[i] = fBinsDirGammapPb5TeVR2PCMPt[i];
                        }
                    }else{
                        maxNBins    = 30;
                        binningMax  = 30;
                        if (centrality.CompareTo("0-100%")){ // applies if it is NOT min bias
                            maxNBins    = 25;
                            binningMax  = 25;
                        }
                        for(Int_t i = 0; i < binningMax+1; i++){
                            if (!centrality.CompareTo("0-100%"))
                                binning[i]  = fBinsDirGammapPb5TeVPCMPt[i];
                            else
                                binning[i]  = fBinsDirGammapPb5TeVCentPCMPt[i];
                        }
                    }
                } else if (mode == 2){
                    maxNBins    = 32;
                    binningMax  = 32;
                    if (centrality.CompareTo("0-100%") || energy.CompareTo("pPb_5.023TeVCent") == 0){ // applies if it is NOT min bias
                        maxNBins    = 25;
                        binningMax  = 25;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if (!centrality.CompareTo("0-100%") && !(energy.CompareTo("pPb_5.023TeVCent") == 0))
                            binning[i]  = fBinsDirGammapPb5TeVPCMEMCPt[i];
                        else
                            binning[i]  = fBinsDirGammapPb5TeVCentPCMPt[i];
                    }
                } else if (mode == 3){
                        maxNBins    = 25;
                        binningMax  = 25;
                    for(Int_t i = 0; i < binningMax+1; i++){
                            binning[i]  = fBinsDirGammapPb5TeVCentPCMPt[i];
                    }
                } else if (mode == 4){
                    maxNBins    = 33;
                    binningMax  = 33;
                    if (centrality.CompareTo("0-100%") || energy.CompareTo("pPb_5.023TeVCent") == 0){ // applies if it is NOT min bias
                        maxNBins    = 25;
                        binningMax  = 25;
                    }
                    for(Int_t i = 0; i < binningMax+1; i++){
                        if (!centrality.CompareTo("0-100%") && !(energy.CompareTo("pPb_5.023TeVCent") == 0))
                            binning[i]  = fBinsDirGammapPb5TeVEMCPt[i];
                        else
                            binning[i]  = fBinsDirGammapPb5TeVCentPCMPt[i];
                    }
                } else if (mode == 20) {
                    maxNBins    = 32;
                    binningMax  = 32;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i]  = fBinsDirGammapPb5TeVPCMEMCPt[i];
                    }
                } else if (mode == 21) {
                    maxNBins    = 41;
                    binningMax  = 41;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i]  = fBinsDirGammapPb5TeVPt[i];
                    }
                } else if (mode == 22) {
                  maxNBins    = 36;
                  binningMax  = 36;
                  for(Int_t i = 0; i < binningMax+1; i++){
                    binning[i]  = fBinsDirGammapPb5TeVAlterPt[i];
                  }
                } else if (mode == 23) {
                    maxNBins    = 37;
                    binningMax  = 37;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i]  = fBinsDirGammapPb5TeVAlter2Pt[i];
                    }
                } else if (mode == 24) {
                    maxNBins    = 28;
                    binningMax  = 28;
                    for(Int_t i = 0; i < binningMax+1; i++){
                        binning[i]  = fBinsDirGammapPb5TeVCombinationPt[i];
                    }
                }
            } else if (energy.CompareTo("pPb_8TeV") == 0  ){
                if (mode == 0){
                    maxNBins = 25;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i]  = fBinsDirGammapPb8TeVPt[i];
                    }
                } else if (mode == 2){
                    maxNBins = 28;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsDirGammapPb8TeVPCMEMCPt[i];
                    }
                } else if (mode == 4){
                    maxNBins = 25;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsDirGammapPb8TeVPt[i];
                    }
                } else if (mode == 20) {
                    maxNBins = 28;
                    for(Int_t i = 0; i < maxNBins+1; i++){
                        binning[i] = fBinsDirGammapPb8TeVPCMEMCPt[i];
                    }
                }
            } else if ( energy.CompareTo("5TeV") == 0 || energy.CompareTo("5TeV2017") == 0 ){
                if ( mode == 0 ){
                    if(energy.Contains("2017")){
                        maxNBins = 44;
                        for(Int_t i = 0; i < maxNBins+1; i++)
                            binning[i] = fBinsDirGamma5TeV2017PCMPt[i];
                    } else {
                        maxNBins = 17;
                        for(Int_t i = 0; i < maxNBins+1; i++)
                            binning[i] = fBinsDirGamma5TeVPt[i];
                    }
                }
            }
        } else if (meson.CompareTo("CKaon") == 0 || meson.CompareTo("CPion") == 0 ){
            maxNBins = 61;
            for(Int_t i = 0; i < maxNBins+1; i++){
                binning[i]  = fBinsInterAndExtrapolationFine[i];
            }
        } else if (meson.CompareTo("Lambda") == 0 || meson.CompareTo("Proton") == 0 || meson.CompareTo("Rho") == 0 || meson.CompareTo("K0Star") == 0 || meson.CompareTo("Phi") == 0){
            maxNBins = 48;
            if(meson.CompareTo("Phi") == 0)
                maxNBins = 36;
            for(Int_t i = 0; i < maxNBins+1; i++){
                binning[i]  = fBinsInterAndExtrapolation[i];
            }
        }
        cout << "maximum " << maxNBins << " pt bins for " << meson.Data() << endl;
        return maxNBins;
    }


    void InitializeClusterBinning( TString energy, Int_t modi ){

        // For heavy meson analysis
        Int_t modeHeavy = modi;
        if(modi>=100) modi -= 100;

        // Set fBinsClusterPt according to cases
        fBinsClusterPt          = new Double_t[400];
        if( energy.CompareTo("2.76TeV") == 0 || energy.CompareTo("PbPb_2.76TeV") == 0 || energy.CompareTo("PbPb_5.02TeV") == 0 ||  energy.CompareTo("5TeV") == 0 || energy.CompareTo("5TeV2017") == 0 || energy.CompareTo("pPb_8TeV") == 0 ){
            fNBinsClusterPt       = fNBinsCluster2760GeVPt;
            for(Int_t iPt=0;iPt<=fNBinsClusterPt;iPt++){
                fBinsClusterPt[iPt] = fBinsCluster2760GeVPt[iPt];
            }
        } else if( energy.CompareTo("pPb_5.023TeV") == 0  || energy.CompareTo("pPb_5.023TeVCent") == 0 || energy.CompareTo("pPb_5.023TeVRun2") == 0 ){
            fNBinsClusterPt       = fNBinsClusterpPb5TeVPt;
            for(Int_t iPt=0;iPt<=fNBinsClusterPt;iPt++){
                fBinsClusterPt[iPt] = fBinsClusterpPb5TeVPt[iPt];
            }

        } else if( energy.CompareTo("7TeV") == 0 || energy.CompareTo("900GeV") == 0){
            fNBinsClusterPt       = fNBinsCluster8TeVPt;
            for(Int_t iPt=0;iPt<=fNBinsClusterPt;iPt++){
                fBinsClusterPt[iPt] = fBinsCluster8TeVPt[iPt];
            }
        } else if( energy.CompareTo("8TeV") == 0 || energy.CompareTo("pPb_8TeV") == 0){
            if(modi == 2 || modi == 4){
                fNBinsClusterPt       = fNBinsCluster8TeVPt;
                for(Int_t iPt=0;iPt<=fNBinsClusterPt;iPt++){
                    fBinsClusterPt[iPt] = fBinsCluster8TeVPt[iPt];
                }
            }else{
                fNBinsClusterPt       = fNBinsCluster8TeVmEMCPt;
                for(Int_t iPt=0;iPt<=fNBinsClusterPt;iPt++){
                    fBinsClusterPt[iPt] = fBinsCluster8TeVmEMCPt[iPt];
                }
            }

        } else if( energy.EqualTo("13TeV") || energy.EqualTo("13TeVLowB") ) {
            if( modi!=0 && modeHeavy<100 ) {
                fNBinsClusterPt            = fNBinsCluster13TeVPt; // 335
                for(Int_t i=0; i<=fNBinsCluster13TeVPt; i++ ){
                    if (i < 1) fBinsCluster13TeVPt[i]          = 0.3*i;
                    else if(i<55) fBinsCluster13TeVPt[i]       = 0.3+0.05*(i-1);
                    else if(i<225) fBinsCluster13TeVPt[i]      = 3.+0.1*(i-55);
                    else if(i<265) fBinsCluster13TeVPt[i]      = 20.+0.25*(i-225);
                    else if(i<305) fBinsCluster13TeVPt[i]      = 30.+0.5*(i-265);
                    else if(i<325) fBinsCluster13TeVPt[i]      = 50.+1.0*(i-305);
                    else if(i<335) fBinsCluster13TeVPt[i]      = 70.+2.5*(i-325);
                    else fBinsCluster13TeVPt[i]                = 100;
                }
                for(Int_t iPt=0;iPt<=fNBinsClusterPt;iPt++){
                    fBinsClusterPt[iPt] = fBinsCluster13TeVPt[iPt];
                }
            } else if( modi==0 || modeHeavy>=100 ){
                fNBinsClusterPt        = fNBinsCluster13TeVPCMPt; // 301
                for(Int_t i=0; i<=fNBinsCluster13TeVPCMPt; i++ ){
                    if (i < 1) fBinsCluster13TeVPCMPt[i]          = 0.3*i;
                    else if(i<55) fBinsCluster13TeVPCMPt[i]       = 0.3+0.05*(i-1);
                    else if(i<125) fBinsCluster13TeVPCMPt[i]      = 3.+0.1*(i-55);
                    else if(i<155) fBinsCluster13TeVPCMPt[i]      = 10.+0.2*(i-125);
                    else if(i<211) fBinsCluster13TeVPCMPt[i]      = 16.+0.25*(i-155);
                    else if(i<251) fBinsCluster13TeVPCMPt[i]      = 30.+0.5*(i-211);
                    else if(i<301) fBinsCluster13TeVPCMPt[i]      = 50.+1.0*(i-251);
                    else fBinsCluster13TeVPCMPt[i]                = 100;
                }
                for(Int_t iPt=0;iPt<=fNBinsCluster13TeVPCMPt;iPt++){
                    fBinsClusterPt[iPt] = fBinsCluster13TeVPCMPt[iPt];
                }
            }
        } else if(  energy.CompareTo("XeXe_5.44TeV") == 0 ){
            fNBinsClusterPt       = fNBinsClusterXeXe5440GeVPt;
            for(Int_t iPt=0;iPt<=fNBinsClusterPt;iPt++){
                fBinsClusterPt[iPt] = fBinsClusterXeXe5440GeVPt[iPt];
            }
        } else {
            fNBinsClusterPt       = 0;
            fBinsClusterPt        = NULL;
        }
    }

    //*************************************************************************************************
    //******************** Determine special trigger set based on energy string & cutnumber ***********
    //*************************************************************************************************
    Int_t GetSpecialTriggerInt( TString energy, TString trigger ) {
        Int_t triggerSetTemp    = -1;
        if (energy.CompareTo("2.76TeV") == 0) {
            if (trigger.CompareTo("52") == 0){
                triggerSetTemp = 1;    // L0
            } else if ( trigger.CompareTo("85") == 0 ){
                triggerSetTemp = 2; //L1 G2 (lower threshold)
            } else if ( trigger.CompareTo("83") == 0    ){
                triggerSetTemp = 3; //L1 G1 (higher threshold)
            } else if ( trigger.CompareTo("51") == 0    ){
                triggerSetTemp = 4; //L0 LHC11a
            } else if ( trigger.CompareTo("01") == 0  || trigger.CompareTo("00") == 0   ){
                triggerSetTemp = 5; //INT7 LHC13g
            }
        } else if (energy.CompareTo("5TeV") == 0) {
            if (trigger.CompareTo("52") == 0){
                triggerSetTemp = 1;    // L0
            } else if ( trigger.CompareTo("85") == 0 ){
                triggerSetTemp = 2; //L1 G2 (lower threshold)
            } else if ( trigger.CompareTo("83") == 0    ){
                triggerSetTemp = 3; //L1 G2 (lower threshold)
            }
        } else if (energy.CompareTo("8TeV") == 0) {
            if (trigger.CompareTo("52") == 0){
                triggerSetTemp = 1; // L0 EMC7
            } else if ( trigger.CompareTo("81") == 0 ){
                triggerSetTemp = 2; //L1 INT7 EGA
            } else if ( trigger.CompareTo("53") == 0 ){
                triggerSetTemp = 3; // L0 EMC8
            } else if ( trigger.CompareTo("82") == 0 ) {
                triggerSetTemp = 4; // L1 INT8 EGA
            }
        } else if (energy.CompareTo("13TeV") == 0 || energy.CompareTo("13TeVRBins") == 0) {
            if     ( trigger.EqualTo("10") ) triggerSetTemp = 0; // MinBias
            else if( trigger.EqualTo("52") ) triggerSetTemp = 1; // L0 EMC7 3GeV
            else if( trigger.EqualTo("83") ) triggerSetTemp = 2; // EG1 8GeV
            else if( trigger.EqualTo("85") ) triggerSetTemp = 3; // EG2 4GeV
            else if( trigger.EqualTo("74") ) triggerSetTemp = 4; // VHM
            else if( trigger.EqualTo("76") ) triggerSetTemp = 5; // VHM+SPD2
            else if( trigger.EqualTo("62") ) triggerSetTemp = 6; // PHOS VZERO 4GeV
        } else if (energy.Contains("pPb_5.023TeV")) {
            if (trigger.CompareTo("52") == 0){
                triggerSetTemp = 1;    // L0
            } else if ( trigger.CompareTo("85") == 0 ){
                triggerSetTemp = 2; //L1 G2 (lower threshold)
            } else if ( trigger.CompareTo("83") == 0    ){
                triggerSetTemp = 3; //L1 G2 (lower threshold)
            } else if ( trigger.CompareTo("62") == 0    ){
                triggerSetTemp = 4; //PHOS PHI7
            } else {
                triggerSetTemp = 0;    // L0
            }
        } else if( energy.CompareTo("pPb_8TeV") == 0) {
            if (trigger.CompareTo("52") == 0){
                triggerSetTemp = 1;    // L0
            } else if ( trigger.CompareTo("85") == 0 ){
                triggerSetTemp = 2; //L1 G2 (lower threshold)
            } else if ( trigger.CompareTo("83") == 0    ){
                triggerSetTemp = 3; //L1 G2 (lower threshold)
            } else if ( trigger.CompareTo("62") == 0    ){
                triggerSetTemp = 4; //PHOS PHI7
            } else {
                triggerSetTemp = 0;    // L0
            }
        }
        return  triggerSetTemp;
    }

    //*************************************************************************************************
    //******************** Initialize binning for analysis stream  ************************************
    //*************************************************************************************************
    void InitializeBinning(
        TString setPi0,
        Int_t numberOfBins,
        TString energy,
        TString directPhoton,
        Int_t modi,
        TString eventCutSelection,
        TString clusterCutSelection,
        Int_t triggerSet = -1,
        Bool_t isDCA = kFALSE,
        TString centDCA = "",
        TString periodDCA = "",
        TString photonCutSelection = ""
    ) {


        //*************************************************************************************************
        //************************************ Binning for Cluster ****************************************
        //*************************************************************************************************

        // Heavy meson analysis
        Int_t modeHeavy = modi;
        if( modi>=100 ) modi -= 100;

        InitializeClusterBinning(energy, modeHeavy);

        //get centrality
        TString centrality      = GetCentralityString(eventCutSelection);
        // set trigger string
        TString trigger         = eventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
        Int_t specialTrigg      = 0;
        Int_t maxPtBinAvail     = 0;

        // Initialize bin for single invariant mass plot
        fExampleBin             = ReturnSingleInvariantMassBinPlotting (setPi0, energy, modi, trigger.Atoi(), fExampleBinScaleFac, triggerSet, directPhoton, centrality);
        cout << "Example pt bin: " <<  fExampleBin << endl;

        if (triggerSet == -1){
            specialTrigg        = GetSpecialTriggerInt(energy, trigger);
        } else {
            specialTrigg        = triggerSet;
        }

	cout<< "specialTrigg::"<< specialTrigg<<endl;
        //*************************************************************************************************
        //************************************ Binning for Pi0 ********************************************
        //*************************************************************************************************
        if (setPi0.CompareTo("Pi0") == 0){
            fNBinsPt                = numberOfBins;
            fBinsPt                 = new Double_t[200];
            fNRebin                 = new Int_t[199];
            //*********************************************************************************************
            //********************************** Pi0 for pp 0.9TeV*****************************************
            //*********************************************************************************************
            if (energy.CompareTo("900GeV") == 0) {
                if (directPhoton.CompareTo("directPhoton") == 0){
                    fStartPtBin     = 1;
                    if( modi == 2){
                      fStartPtBin = 3;
                    }else if(modi == 4){
                      fStartPtBin = 6;
                    }

                    if (fNBinsPt > 13) {
                        cout << "You have chosen Direct Photon Plots and more than 13 bins, this is not possible, it will be reduced to 14 bins." << endl;
                        fNBinsPt    = 13;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        fBinsPt[i]  = fBinsDirGamma900GeVPt[i];
                        if (i < fNBinsPt+1)
                            fNRebin[i] = fBinsDirGamma900GeVPtRebin[i];
                    }
                } else {
                    fStartPtBin     = 1;
                    if (fNBinsPt > 11) {
                        cout << "You have chosen to have more than 11 bins, this is not possible, it will be reduced to 11" << endl;
                        fNBinsPt    = 11;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                    if(modi != 2 && modi !=4){
                        fBinsPt[i]  = fBinsPi0900GeVPt[i];
                        if (i < fNBinsPt+1)
                            fNRebin[i] = fBinsPi0900GeVPtRebin[i];
                    } else if(modi == 2){
                        fBinsPt[i]  = fBinsPi0900GeVPCMEMCPt[i];
                        if (i < fNBinsPt+1)
                            fNRebin[i] = fBinsPi0900GeVPCMEMCPtRebin[i];
                    } else if(modi == 4){
                        fBinsPt[i]  = fBinsPi0900GeVEMCPt[i];
                        if (i < fNBinsPt+1)
                            fNRebin[i] = fBinsPi0900GeVEMCPtRebin[i];
                    }
                    }
                    nIterBGFit      = 11;
                }
            //*********************************************************************************************
            //********************************** Pi0 for pp 2.76TeV****************************************
            //*********************************************************************************************
            } else if (energy.CompareTo("2.76TeV") == 0) {
                if (directPhoton.CompareTo("directPhoton") == 0){
                    fStartPtBin     = 1;
                    if (modi == 2)
                        fStartPtBin = 3;

                    if (fNBinsPt > 14 && isDCA) {
                        cout << "You have chosen to have more than 14 bins, this is not possible, it will be reduced to 14" << endl;
                        fNBinsPt    = 14;
                    } else if (fNBinsPt > 21 && specialTrigg == 5 && modi!=0) {
                        cout << "You have chosen Direct Photon Plots and more than 21 bins, this is not possible, it will be reduced to 21 bins." << endl;
                        fNBinsPt    = 21;
                    } else if (fNBinsPt > 21 && modi ==0) {
                        cout << "You have chosen Direct Photon Plots and more than 21 bins, this is not possible, it will be reduced to 21 bins." << endl;
                        fNBinsPt    = 21;
                    } else if (fNBinsPt > 24 && modi!=0) {
                        cout << "You have chosen Direct Photon Plots and more than 24 bins, this is not possible, it will be reduced to 24 bins." << endl;
                        fNBinsPt    = 24;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        fBinsPt[i]  = fBinsDirGamma2760GeVPt[i];
                        if (i < fNBinsPt+1)
                            fNRebin[i] = fBinsDirGamma2760GeVPtRebin[i];
                    }
                } else {
                    fStartPtBin     = 1;
                    if (modi == 2 && specialTrigg == 0)
                        fStartPtBin = 3;
                    else if (modi == 2 && specialTrigg == 1)
                        fStartPtBin = 10;
                    else if (modi == 4 && specialTrigg == 1)
                        fStartPtBin = 13;
                    else if (modi == 2 && specialTrigg == 2)
                        fStartPtBin = 15;
                    else if (modi == 4 && specialTrigg == 2)
                        fStartPtBin = 16;
                    else if (modi == 2 && specialTrigg == 3)
                        fStartPtBin = 16;
                    else if (modi == 4 && specialTrigg == 3)
                        fStartPtBin = 18;
                    else if (modi == 2 && specialTrigg == 4)
                        fStartPtBin = 12;
                    else if (modi == 4 && specialTrigg == 4)
                        fStartPtBin = 15;
                    else if (modi == 10 && ( ReturnClusterNLM(clusterCutSelection) == 2 || ReturnClusterNLM(clusterCutSelection) == 0))
                        fStartPtBin = 17;
                    else if (modi == 10 && ReturnClusterNLM(clusterCutSelection) == 1 && (specialTrigg == 0 || specialTrigg == 5))
                        fStartPtBin = 21;
                    else if (modi == 10 && ReturnClusterNLM(clusterCutSelection) == 1)
                        fStartPtBin = 21; // 20 with V1 clusterizer
                    else if (modi == 4 )
                        fStartPtBin = 6;
                    else if (modi == 1 )
                        fStartPtBin = 3;

                    if (fNBinsPt > 14 && isDCA) {
                        cout << "You have chosen to have more than 14 bins, this is not possible, it will be reduced to 14" << endl;
                        fNBinsPt    = 14;
                    } else if (fNBinsPt > 19 && ( modi == 0 || modi == 1) && specialTrigg < 1) {
                        cout << "You have chosen to have more than 19 bins, this is not possible, it will be reduced to 19" << endl;
                        fNBinsPt    = 19;
                    } else if (fNBinsPt > 24 &&  modi == 0  && specialTrigg > 0) {
                        cout << "You have chosen to have more than 19 bins, this is not possible, it will be reduced to 19" << endl;
                        fNBinsPt    = 24;
                    } else if (fNBinsPt > 24 && (modi == 2 || modi == 3) && specialTrigg == 0){
                        cout << "You have chosen to have more than 24 bins, this is not possible, it will be reduced to 24" << endl;
                        fNBinsPt    = 24;
                    } else if (fNBinsPt > 29 && ( (modi == 2 && (specialTrigg == 1 || specialTrigg == 2 || specialTrigg == 3)) || modi ==4)){
                        cout << "You have chosen to have more than 29 bins, this is not possible, it will be reduced to 29" << endl;
                        fNBinsPt    = 29;
                    } else if (fNBinsPt > 28 &&  modi == 2 && specialTrigg == 3){
                        cout << "You have chosen to have more than 28 bins, this is not possible, it will be reduced to 28" << endl;
                        fNBinsPt    = 28;
                    } else if (fNBinsPt > 25 && ( (modi == 2  && (specialTrigg == 4 || specialTrigg == 1 || specialTrigg == 2 )) || (modi == 3 && specialTrigg == 4) )){
                        cout << "You have chosen to have more than 25 bins, this is not possible, it will be reduced to 25" << endl;
                        fNBinsPt    = 25;
                    } else if (fNBinsPt > 32 && (modi == 10)){
                        cout << "You have chosen to have more than 32 bins, this is not possible, it will be reduced to 32" << endl;
                        fNBinsPt    = 32;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        if (isDCA)
                            fBinsPt[i]          = fBinsPi02760GeVPtDCA[i];
                        else
                            fBinsPt[i]          = fBinsPi02760GeVPt[i];
                        if ((modi == 2) && specialTrigg == 0 ){
                            if (i < fNBinsPt+1)
                                fNRebin[i]  = fBinsPi02760GeVPCMEMCPtRebin[i];
                        } else if ( modi == 4 && (specialTrigg == 1 || specialTrigg == 2 || specialTrigg == 3 || specialTrigg == 4 ) ) {
                            if (i < fNBinsPt+1)
                                fNRebin[i]  = fBinsPi02760GeVEMCPtTrig13gRebin[i];
                            fBinsPt[i]      = fBinsPi02760GeVPtTrig13g[i];
                        } else if ( modi == 2 && specialTrigg == 3 ) {
                            if (i < fNBinsPt+1)
                                fNRebin[i]  = fBinsPi02760GeVPCMEMCPtTrig13gRebin[i];
                            fBinsPt[i]      = fBinsPi02760GeVPtTrig13gPCMEMC[i];
                        } else if ( modi == 2 && (specialTrigg == 4 || specialTrigg == 1 || specialTrigg == 2 ) ){
                            if (i < fNBinsPt+1)
                                fNRebin[i]  = fBinsPi02760GeVPCMEMCPtTrig11aRebin[i];
                            fBinsPt[i]      = fBinsPi02760GeVPtTrig11a[i];
                        } else if ( modi == 0 && (specialTrigg == 3 || specialTrigg == 1 || specialTrigg == 2 ) ){
                            if (i < fNBinsPt+1)
                                fNRebin[i]  = fBinsPi02760GeVPCMEMCPtTrig11aRebin[i];
                            fBinsPt[i]      = fBinsPi02760GeVPtTrig11a[i];

                        } else if ( modi == 10 ) {
                            if (i < fNBinsPt+1)
                                fNRebin[i]  = fBinsPi02760GeVPtmEMCRebin[i];
                            fBinsPt[i]      = fBinsPi02760GeVPtmEMC[i];
                        } else if ( modi == 1 ) {
                            if (i < fNBinsPt+1)
                                fNRebin[i]  = fBinsPi02760GeVDalitzPtRebin[i];
                            fBinsPt[i]      = fBinsPi02760GeVDalitzPt[i];
                        } else {
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsPi02760GeVPtRebin[i];
                        }
                    }
                    fMaxYFracBGOverIntHist              = 50;
                    nIterBGFit                          = 13;
                    optionBGSmoothingStandard           = "BackSmoothing9";
                    optionBGSmoothingVar1               = "BackSmoothing7";
                    optionBGSmoothingVar2               = "BackSmoothing11";
                }
            //*********************************************************************************************
            //********************************** Pi0 for pp 5TeV*******************************************
            //*********************************************************************************************
            } else if (energy.CompareTo("5TeV") == 0 || energy.CompareTo("5TeV2017") == 0) {
                if (directPhoton.CompareTo("directPhoton") == 0){
                    fStartPtBin                 = GetStartBin("Pi0", energy, modi, specialTrigg, centrality);
                    Int_t maxPtBinTheo          = GetBinning( fBinsPt, maxPtBinAvail, "Pi0", energy, modi, specialTrigg, isDCA, centrality );
                    if (fNBinsPt > maxPtBinTheo) {
                        cout << "**************************************************************************************************************************************" << endl;
                        cout << "********************** ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, **********************************" << endl;
                        cout << "You have chosen "<< fNBinsPt << " bins, which is more than the maximal " << maxPtBinTheo << " bins, this is not possible, it will be reduced to " << maxPtBinTheo << endl;
                        cout << "**************************************************************************************************************************************" << endl;
                        fNBinsPt    = maxPtBinTheo;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                    for (Int_t i = 0; i < fNBinsPt; i++) {
                        if ( modi == 0 ) {
                            if(energy.Contains("2017")){
                                fNRebin[i]  = fBinsDirGamma5TeV2017PCMPtRebin[i];
                            } else {
                                fNRebin[i]  = fBinsDirGamma5TeVPtRebin[i];
                            }
                        }
                    }

                    fNBinsPtDCAzDist    = 45;
                    fBinsPtDCAzDist     = new Double_t[fNBinsPtDCAzDist+1];
                    for (Int_t i = 0; i < fNBinsPtDCAzDist+1; i++) {
                        fBinsPtDCAzDist[i] = fBinsDirGamma5TeV2017PCMPt[i]; //fBinsDirGamma5TeV2017PCMPtDCA[i];
                    }

                } else {
                  fStartPtBin                 = GetStartBin("Pi0", energy, modi, specialTrigg, centrality);
                  Int_t maxPtBinTheo          = GetBinning( fBinsPt, maxPtBinAvail, "Pi0", energy, modi, specialTrigg, isDCA, centrality );
                  if (fNBinsPt > maxPtBinTheo) {
                      cout << "**************************************************************************************************************************************" << endl;
                      cout << "********************** ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, **********************************" << endl;
                      cout << "You have chosen "<< fNBinsPt << " bins, which is more than the maximal " << maxPtBinTheo << " bins, this is not possible, it will be reduced to " << maxPtBinTheo << endl;
                      cout << "**************************************************************************************************************************************" << endl;
                      fNBinsPt    = maxPtBinTheo;
                  }
                  GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                    for (Int_t i = 0; i < fNBinsPt; i++) {
                        if ( modi == 0 ) {
                            if(energy.Contains("2017")){
                              if(fNBinsPt>15 && fNBinsPt<22)
                                fNRebin[i]      = fBinsPi05TeV2017PCMforPbPbPtRebin[i]; //fBinsPi05TeVPtRebin[i];
                              else if(fNBinsPt>22 && fNBinsPt<26)
                                fNRebin[i]      = fBinsPi05TeVPtRebin[i];
                              else if(fNBinsPt>28 && fNBinsPt<32)
                                fNRebin[i]      = fBinsPi05TeV2017PCMCombinationPtRebin[i];
                              else if(fNBinsPt>39 && fNBinsPt<45)
                                fNRebin[i]      = fBinsPi05TeV2017PtRebin[i];
                              else if(fNBinsPt>60&&fNBinsPt<70)
                                fNRebin[i]      = fBinsPi05TeV2017ExtraFinePtRebin[i];
                            } else {
                                fNRebin[i]  = fBinsPi05TeVPtRebin[i];
                            }
                        } else if ( modi == 1 ) {
                          fNRebin[i]  = fBinsPi05TeV2017DalitzPtRebin[i];
                        } else if ( modi == 2 ) {
                          if(energy.Contains("2017"))
                            fNRebin[i] = fBinsPi05TeV2017PCMEMCPtRebin[i];
                          else
                            fNRebin[i] = fBinsPi05TeVPCMEMCPtRebin[i];
                        } else if ( modi == 3 ) {
                            fNRebin[i] = fBinsPi05TeV2017PtCombinationRebin[i];
                        } else if ( modi == 4 ) {
                          if(energy.Contains("2017")){
                            fNRebin[i] = fBinsPi05TeV2017PtCombinationRebin[i];
                          } else {
                            if(specialTrigg == 1 || specialTrigg == 2 || specialTrigg == 3)
                              fNRebin[i] = fBinsPi05TeVEMCPtRebinTrigger1[i];
                            else
                              fNRebin[i] = fBinsPi05TeVEMCPtRebin[i];
                          }
                        } else if ( modi == 10 ) {
                            fNRebin[i] = fBinsPi05TeV2017mEMCPtRebin[i];
                        } else if ( modi == 12 ) {
                          if(energy.Contains("2017"))
                            fNRebin[i] = fBinsPi05TeV2017DMCPtRebin[i];
                          else
                            fNRebin[i] = fBinsPi05TeVPtRebinDMC[i];
                        } else if ( modi == 13 ) {
                          if(energy.Contains("2017"))
                            fNRebin[i] = fBinsPi05TeV2017PtRebinPCMDCal[i];
                          else
                            fNRebin[i] = fBinsPi05TeVPtRebinPCMDCal[i];
                        } else  {
                          fNRebin[i]  = fBinsPi05TeVPtRebin[i];
                        }

                    }
                    if(modi == 0 && energy.Contains("2017")){
                        nIterBGFit                  = 8;
                        fMaxYFracBGOverIntHist      = 70;
                        optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                        optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                        optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing7";
                    } else {
                        nIterBGFit                  = 10;
                        fMaxYFracBGOverIntHist      = 60;
                        optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing5";
                        optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing7";
                        optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing3";
                    }
                }
            //*********************************************************************************************
            //********************************** Pi0 for pp 7TeV*******************************************
            //*********************************************************************************************
            } else if (energy.CompareTo("7TeV") == 0) {
                if (directPhoton.CompareTo("directPhoton") == 0){
                    fStartPtBin     = 1;
                    if (modi == 4)
                        fStartPtBin = 6;
                    else if(modi == 2)
                        fStartPtBin = 4;
                    else if(modi == 1)
                        fStartPtBin = 6;
                    if (fNBinsPt > 23) {
                        cout << "You have chosen Direct Photon Plots and more than 23 bins, this is not possible, it will be reduced to 23 bins." << endl;
                        fNBinsPt    = 23;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        fBinsPt[i]         = fBinsDirGamma7TeVPt[i];
                        if (i < fNBinsPt+1) {
                            if (modi == 4)
                                fNRebin[i] = fBinsDirGamma7TeVEMCPtRebin[i];
                            else
                                fNRebin[i] = fBinsDirGamma7TeVPtRebin[i];
                        }
                    }
                }else if (directPhoton.CompareTo("directPhotonTagging") == 0){
                    fStartPtBin     = 1;
                    if (fNBinsPt > 23) {
                        cout << "You have chosen Direct Photon Plots and more than 23 bins, this is not possible, it will be reduced to 23 bins." << endl;
                        fNBinsPt    = 23;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        fBinsPt[i]         = fBinsDirGamma7TeVPt[i];
                        if (i < fNBinsPt+1) {
                              fNRebin[i] = fBinsDirGamma7TeVPtRebin[i];
                        }
                    }
                } else {
                    fStartPtBin     = 1;
                    if (modi == 4)
                        fStartPtBin = 10;
                    else if(modi == 2)
                        fStartPtBin = 6;
                    if(modi == 5)
                        fStartPtBin = 2;

                    if (fNBinsPt > 27 && isDCA) {
                        cout << "You have chosen to have more than 27 bins, this is not possible, it will be reduced to 27" << endl;
                        fNBinsPt    = 27;
                    } else if (fNBinsPt > 40) {
                        cout << "You have chosen to have more than 40 bins, this is not possible, it will be reduced to 40" << endl;
                        fNBinsPt    = 40;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        if (isDCA){
                            fBinsPt[i]          = fBinsPi07TeVPtDCA[i];
                        } else {
                            if (modi == 2){
                                fBinsPt[i]      = fBinsPi07TeVPCMEMCPt[i];
                                if (i < fNBinsPt)
                                    fNRebin[i]  = fBinsPi07TeVPCMEMCPtRebin[i];
                            } else if (modi == 3){
                                fBinsPt[i]      = fBinsPi07TeVPCMPHOSPt[i];
                                if (i < fNBinsPt)
                                fNRebin[i]   = fBinsPi07TeVPCMPHOSPtRebin[i];
                            } else if (modi == 4){
                                fBinsPt[i]      = fBinsPi07TeVEMCPt[i];
                                if (i < fNBinsPt)
                                fNRebin[i]   = fBinsPi07TeVEMCPtRebin[i];
                            } else if (modi == 5){
                                fBinsPt[i]      = fBinsPi07TeVPCMPHOSPt[i];
                                if (i < fNBinsPt)
                                fNRebin[i]   = fBinsPi07TeVPCMPHOSPtRebin[i];
                            } else if (modi == 1){
                                fBinsPt[i]      = fBinsPi07TeVDalitzPt[i];
                                if (i < fNBinsPt)
                                fNRebin[i]   = fBinsPi07TeVDalitzPtRebin[i];
                            } else {
                                fBinsPt[i]      = fBinsPi07TeVPt[i];
                                if (i < fNBinsPt)
                                    fNRebin[i]  = fBinsPi07TeVPtRebin[i];
                            }
                        }
                    }
                    nIterBGFit                  = 9;
                    fMaxYFracBGOverIntHist      = 40;
                    optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing5";
                    optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing7";
                    optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing3";
                }
            //*********************************************************************************************
            //********************************** Pi0 for pp 8TeV*******************************************
            //*********************************************************************************************
            } else if (energy.CompareTo("8TeV") == 0) {
                if (directPhoton.CompareTo("directPhoton") == 0){
                    fStartPtBin     = 1;
                    if( modi == 2){
                      fStartPtBin = 4;
                    }else if(modi == 4){
                      fStartPtBin = 6;
                    }
                    if (modi == 4 && specialTrigg == 1){ fStartPtBin = 22; }
                    if (modi == 4 && specialTrigg == 2){ fStartPtBin = 33; }

                    if (fNBinsPt > 32 && (modi == 4 || modi == 2) ){
                      if( specialTrigg == 2 && fNBinsPt > 41){
                        cout << "You have chosen to have more than 41 bins, this is not possible, it will be reduced to 41" << endl;
                        fNBinsPt                        = 41;
                      } else if ( specialTrigg == 1 && fNBinsPt > 41){
                        cout << "You have chosen to have more than 41 bins, this is not possible, it will be reduced to 41" << endl;
                        fNBinsPt                        = 41;
                      } else if (specialTrigg!=1 && specialTrigg!=2){
                        cout << "You have chosen to have more than 23 bins, this is not possible, it will be reduced to 23" << endl;
                        fNBinsPt                        = 23;
                      }
                    } else if (fNBinsPt > 23 && (modi !=4 && modi !=2) ) {
                        cout << "You have chosen Direct Photon Plots and more than 23 bins, this is not possible, it will be reduced to 23 bins." << endl;
                        fNBinsPt    = 23;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        fBinsPt[i]         = fBinsDirGamma8TeVPt[i];
                        if (i < fNBinsPt+1) fNRebin[i] = fBinsDirGamma8TeVPtRebin[i];

                        if (modi == 4 ){
                            if( specialTrigg == 1 ){
                                fBinsPt[i]                 = fBinsDirGamma8TeVEMCalTriggerPt[i];
                            } else if ( specialTrigg == 2 ){
                                fBinsPt[i]                 = fBinsDirGamma8TeVEMCalTriggerPt[i];
                            } else
                                fBinsPt[i]                 = fBinsDirGamma8TeVPt[i];
                        }

                        if(modi == 4 ) {
                            if( specialTrigg == 1 ){
                                if (i < fNBinsPt+1) fNRebin[i] = fBinsDirGamma8TeVEMCalTriggerPtRebin[i];
                            } else if( specialTrigg == 2 ){
                                if (i < fNBinsPt+1) fNRebin[i] = fBinsDirGamma8TeVEMCalTriggerPtRebin[i];
                            } else{
                                if (i < fNBinsPt+1) fNRebin[i] = fBinsDirGamma8TeVEMCalPtRebin[i];
                            }
                        }
                    }
                } else if (directPhoton.CompareTo("directPhotonTagging") == 0){
                    fStartPtBin     = 1;
                    if( modi == 2){
                        fStartPtBin = 1;
                    }
                    if (fNBinsPt > 29 ) {
                        cout << "You have chosen Direct Photon Plots and more than 29 bins, this is not possible, it will be reduced to 29 bins." << endl;
                        fNBinsPt    = 29;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        fBinsPt[i]         = fBinsDirGammaTagging8TeVPt[i];
                        if (i < fNBinsPt+1) fNRebin[i] = fBinsDirGammaTagging8TeVPtRebin[i];
                    }
                } else {
                    fStartPtBin = 1;
                    if (modi == 4 ) fStartPtBin = 7;
                    if (modi == 2 ) fStartPtBin = 2;
                    if (modi == 0 && specialTrigg == 1) fStartPtBin = 21;
                    if (modi == 0 && specialTrigg == 2) fStartPtBin = 27;
                    if (modi == 2 && specialTrigg == 1) fStartPtBin = 21;
                    if (modi == 2 && specialTrigg == 2) fStartPtBin = 28;
                    if (modi == 4 && specialTrigg == 1) fStartPtBin = 22;
                    if (modi == 4 && specialTrigg == 2) fStartPtBin = 33;
                    if (modi == 10 && specialTrigg == 0) fStartPtBin = 28;
                    if (modi == 10 && specialTrigg == 1) fStartPtBin = 28;
                    if (modi == 10 && specialTrigg == 2) fStartPtBin = 28;

                    if (fNBinsPt > 21 && isDCA) {
                        cout << "You have chosen to have more than 21 bins, this is not possible, it will be reduced to 21" << endl;
                        fNBinsPt                        = 21;
                    } else if (fNBinsPt > 33 && modi != 2 && modi != 3 && modi != 4 && modi != 10) {
                        if ( specialTrigg == 1 && fNBinsPt > 42){
                            cout << "You have chosen to have more than 42 bins, this is not possible, it will be reduced to 42" << endl;
                            fNBinsPt                    = 42;
                        } else if ( specialTrigg == 2 && fNBinsPt > 41){
                            cout << "You have chosen to have more than 41 bins, this is not possible, it will be reduced to 41" << endl;
                            fNBinsPt                    = 41;
                        } else if (specialTrigg!=1 && specialTrigg!=2){
                            cout << "You have chosen to have more than 33 bins, this is not possible, it will be reduced to 32" << endl;
                            fNBinsPt                    = 33;
                        }
                    } else if (fNBinsPt > 32 && (modi ==4)){
                        if( specialTrigg == 2 && fNBinsPt > 42){
                            cout << "You have chosen to have more than 42 bins, this is not possible, it will be reduced to 42" << endl;
                            fNBinsPt                        = 42;
                        } else if ( specialTrigg == 1 && fNBinsPt > 42){
                            cout << "You have chosen to have more than 42 bins, this is not possible, it will be reduced to 42" << endl;
                            fNBinsPt                        = 42;
                        } else if (specialTrigg!=1 && specialTrigg!=2){
                            cout << "You have chosen to have more than 32 bins, this is not possible, it will be reduced to 32" << endl;
                            fNBinsPt                        = 32;
                        }
                    } else if (fNBinsPt > 29 && (modi == 2 || modi == 3)){
                        if( specialTrigg == 2 && fNBinsPt > 42){
                        cout << "You have chosen to have more than 42 bins, this is not possible, it will be reduced to 42" << endl;
                        fNBinsPt        = 42;
                        } else if ( specialTrigg == 1 && fNBinsPt > 42){
                        cout << "You have chosen to have more than 42 bins, this is not possible, it will be reduced to 42" << endl;
                        fNBinsPt = 42;
                        } else if(specialTrigg!=1 && specialTrigg!=2) {
                        cout << "You have chosen to have more than 31 bins, this is not possible, it will be reduced to 31" << endl;
                        fNBinsPt        = 29;
                        }
                    } else if (fNBinsPt > 59 && (modi == 10)){
                        cout << "You have chosen to have more than 59 bins, this is not possible, it will be reduced to 59" << endl;
                        fNBinsPt        = 59;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        if (modi == 4 ){
                            if( specialTrigg == 1 ){
                                fBinsPt[i]                 = fBinsPi08TeVEMCalTrigger1Pt[i];
                            } else if ( specialTrigg == 2 ){
                                fBinsPt[i]                 = fBinsPi08TeVEMCalTrigger2Pt[i];
                            } else
                                fBinsPt[i]                 = fBinsPi08TeVPtEMC[i];
                        } else if (modi == 2){
                            if( specialTrigg == 1 ){
                                fBinsPt[i]                 = fBinsPi08TeVPCMEMCalTrigger1Pt[i];
                            } else if ( specialTrigg == 2 ){
                                fBinsPt[i]                 = fBinsPi08TeVPCMEMCalTrigger2Pt[i];
                            } else
                                fBinsPt[i]                 = fBinsPi08TeVPtPCMEMC[i];
                        } else if (modi == 10){
                            fBinsPt[i]                     = fBinsPi08TeVPtmEMC[i];
                        } else {
                            if( specialTrigg == 1 ){
                                fBinsPt[i]                 = fBinsPi08TeVTrigger1Pt[i];
                            } else if ( specialTrigg == 2 ){
                                fBinsPt[i]                 = fBinsPi08TeVPCMTrigger2Pt[i];
                            } else if ( isDCA ) {
                                fBinsPt[i]                 = fBinsPi08TeVPtDCA[i];
                            } else
                                fBinsPt[i]                 = fBinsPi08TeVPt[i];
                        }

                        if (modi == 2){
                            if( specialTrigg == 1 ){
                                if (i < fNBinsPt+1) fNRebin[i] = fBinsPi08TeVPCMEMCTrigger1PtRebin[i];
                            } else if( specialTrigg == 2 ){
                                if (i < fNBinsPt+1) fNRebin[i] = fBinsPi08TeVPCMEMCTrigger2PtRebin[i];
                            } else{
                                if (i < fNBinsPt+1) fNRebin[i] = fBinsPi08TeVPCMEMCPtRebin[i];
                            }
                        } else if(modi == 4) {
                            if( specialTrigg == 1 ){
                                if (i < fNBinsPt+1) fNRebin[i] = fBinsPi08TeVEMCTrigger1PtRebin[i];
                            } else if( specialTrigg == 2 ){
                                if (i < fNBinsPt+1) fNRebin[i] = fBinsPi08TeVEMCTrigger2PtRebin[i];
                            } else{
                                if (i < fNBinsPt+1) fNRebin[i] = fBinsPi08TeVEMCPtRebin[i];
                            }
                        } else if(modi == 10) {
                            if (i < fNBinsPt+1) fNRebin[i]     = fBinsPi08TeVPtmEMCRebin[i];
                        } else {
                            if( specialTrigg == 1 ){
                                if (i < fNBinsPt+1) fNRebin[i] = fBinsPi08TeVPCMTrigger1PtRebin[i];
                            } else if( specialTrigg == 2 ){
                                if (i < fNBinsPt+1) fNRebin[i] = fBinsPi08TeVPCMTrigger2PtRebin[i];
                            } else{
                                if (i < fNBinsPt+1) fNRebin[i] = fBinsPi08TeVPtRebin[i];
                            }
                        }
                    }
                    nIterBGFit                  = 7;
                    fMaxYFracBGOverIntHist      = 60;
                    optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing5";
                    optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing7";
                    optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing3";

                }
            //*********************************************************************************************
            //********************************** Pi0 for pp 13TeV******************************************
            //*********************************************************************************************
            } else if (energy.CompareTo("13TeV") == 0 || energy.CompareTo("13TeVRBins") == 0) {
                if (directPhoton.CompareTo("directPhoton") == 0){
                    fStartPtBin                 = GetStartBin("Pi0", energy, modi, specialTrigg, centrality);
                    if (fNBinsPt > 24) {
                        cout << "You have chosen Direct Photon Plots and more than 24 bins, this is not possible, it will be reduced to 24 bins." << endl;
                        fNBinsPt    = 24;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        fBinsPt[i]         = fBinsDirGamma13TeVPt[i];
                        if (i < fNBinsPt+1) fNRebin[i] = fBinsDirGamma13TeVPtRebin[i];
                    }
                    fNBinsPtDCAzDist    = 15;
                    fBinsPtDCAzDist     = new Double_t[fNBinsPtDCAzDist+1];
                    for (Int_t i = 0; i < fNBinsPtDCAzDist+1; i++) {
                        fBinsPtDCAzDist[i] = fBinsDirGamma13TeVPtDCAzDist[i];
                    }
                }  else {//13TeV, not directPhoton
                    fStartPtBin                 = GetStartBin("Pi0", energy, modi, specialTrigg, centrality);
                    GetBinning( fBinsPt, maxPtBinAvail, "Pi0", energy, modi, specialTrigg, isDCA);
                    CheckBinSize(fNBinsPt,maxPtBinAvail,kTRUE);
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                    //Rebinning, because not implemented in getBinning
                    if (!isDCA) {
                        cout<<"ReBinning for Pi0, modi: "<<modi<<endl;
                        for (Int_t i = 0; i < fNBinsPt; i++) {
                            if (modi==0){
                                if (specialTrigg == 0 || specialTrigg == 4 || specialTrigg == 5){
                                    if (energy.Contains("RBins")){
                                        fNRebin[i]      = fBinsPi013TeVPCMTrigINT7RBinsPtRebin[i];
                                    }else {
                                        fNRebin[i]      = fBinsPi013TeVPCMTrigINT7PtRebin[i];
                                    }
                                } else if (specialTrigg==1){
                                    fNRebin[i]      = fBinsPi013TeVPCMTrigEMC7PtRebin[i];
                                } else if (specialTrigg==2){
                                    fNRebin[i]      = fBinsPi013TeVPCMTrigEG1PtRebin[i];
                                } else if (specialTrigg==3){
                                    fNRebin[i]      = fBinsPi013TeVPCMTrigEG2PtRebin[i];
                                }
                            } else if( modi == 2 ){
                                if (specialTrigg == 0 || specialTrigg == 4 || specialTrigg == 5){
                                    fNRebin[i]      = fBinsPi013TeVPCMEMCTrigINT7PtRebin[i];
                                } else if (specialTrigg==3){
                                    fNRebin[i]      = fBinsPi013TeVPCMEMCTrigEG2PtRebin[i];
                                } else if (specialTrigg==2){
                                    fNRebin[i]      = fBinsPi013TeVPCMEMCTrigEG1PtRebin[i];
                                }
                            } else if (modi == 3){
                                fNRebin[i]=fBinsPi013TeVPCMPHOSTrigINT7PtRebin[i];
                            } else if( modi == 4){
                                if (specialTrigg == 0 || specialTrigg == 4 || specialTrigg == 5){
                                    fNRebin[i]      = fBinsPi013TeVEMCTrigINT7PtRebin[i];
                                } else if (specialTrigg==1){
                                    fNRebin[i]      = fBinsPi013TeVEMCTrigEMC7PtRebin[i];
                                } else if (specialTrigg==3){
                                    fNRebin[i]      = fBinsPi013TeVEMCTrigEG2PtRebin[i];
                                } else if (specialTrigg==2){
                                    fNRebin[i]      = fBinsPi013TeVEMCTrigEG1PtRebin[i];
                                }
                            } else {
                                fNRebin[i]      = fBinsPi013TeVPCMEMCTrigINT7PtRebin[i];
                            }
                        }
                    }


                // 		    nIterBGFit                  = 7;
                //		    fMaxYFracBGOverIntHist      = 60;
                //		    optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                //		    optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing7";
                //		    optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing3";
                TString rBin = photonCutSelection(2,1);

                if (rBin.CompareTo("2") ==0){
                      nIterBGFit                  = 6;
                      fMaxYFracBGOverIntHist      = 100;
                      //                      optionBGSmoothingStandard   = "BackDecreasingWindow,NoSmoothing";
                      optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                      optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                      optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing7";
                    }else if( rBin.CompareTo("a") ==0){
                      nIterBGFit                  = 7;
                      fMaxYFracBGOverIntHist      = 100;
                      optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                      optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                      optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing7";
                    }else if( rBin.CompareTo("b") ==0){
                      nIterBGFit                  = 7;
                      fMaxYFracBGOverIntHist      = 100;
                      //                      optionBGSmoothingStandard   = "BackDecreasingWindow,NoSmoothing";
                      optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                      optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                      optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing7";
                    }else if( rBin.CompareTo("c") ==0){
                      nIterBGFit                  = 5;
                      fMaxYFracBGOverIntHist      = 110;
                      //optionBGSmoothingStandard   = "BackDecreasingWindow,NoSmoothing";
                      optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                      optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing3";
                      optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing3";
                      //fStartPtBin  = 1;
                    }else{
                      nIterBGFit                  = 6;
                      fMaxYFracBGOverIntHist      = 100;
                      //                      optionBGSmoothingStandard   = "BackDecreasingWindow,NoSmoothing";
                      optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                      optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                      optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing7";
                    }
                }
            //*********************************************************************************************
            //**************************************** Pi0 for 13TeV low B field **************************
            //*********************************************************************************************
            } else if (energy.CompareTo("13TeVLowB") == 0) {
                if (directPhoton.CompareTo("directPhoton") == 0){
                    fStartPtBin     = 1;
                    if (fNBinsPt > 24) {
                        cout << "You have chosen Direct Photon Plots and more than 24 bins, this is not possible, it will be reduced to 24 bins." << endl;
                        fNBinsPt    = 24;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        fBinsPt[i]         = fBinsDirGamma13TeVLowBPt[i];
                        if (i < fNBinsPt+1) fNRebin[i] = fBinsDirGamma13TeVLowBPtRebin[i];
                    }
                    fNBinsPtDCAzDist    = 15;
                    fBinsPtDCAzDist     = new Double_t[fNBinsPtDCAzDist+1];
                    for (Int_t i = 0; i < fNBinsPtDCAzDist+1; i++) {
                        fBinsPtDCAzDist[i] = fBinsDirGamma13TeVLowBPtDCAzDist[i];
                    }
                } else {
                    fStartPtBin                 = GetStartBin("Pi0", energy, modi, specialTrigg, centrality);
                    if (photonCutSelection.CompareTo("")) {
                        if (!((TString)photonCutSelection(GetPhotonSinglePtCutPosition(photonCutSelection),1)).CompareTo("0")){
                            cout << "Increase starting pT for higher minimum track pT cut (50 MeV)" << endl;
                            fStartPtBin += 0;
                        }
                    }
                    Int_t maxPtBinTheo          = GetBinning( fBinsPt, maxPtBinAvail, "Pi0", energy, modi, specialTrigg, isDCA, centrality );
                    if (fNBinsPt > maxPtBinTheo) {
                        cout << "**************************************************************************************************************************************" << endl;
                        cout << "********************** ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, **********************************" << endl;
                        cout << "You have chosen "<< fNBinsPt << " bins, which is more than the maximal " << maxPtBinTheo << " bins, this is not possible, it will be reduced to " << maxPtBinTheo << endl;
                        cout << "**************************************************************************************************************************************" << endl;
                        fNBinsPt    = maxPtBinTheo;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                    for (Int_t i = 0; i < fNBinsPt; i++) {
                        if ( modi == 0 ) {
                            fNRebin[i]  = fBinsPi013TeVLowBPtRebin[i];
                        } else if ( modi == 4 || modi == 12) {
                            fNRebin[i]  = fBinsPi013TeVLowBEMCPtRebin[i];
                        } else if ( modi == 5 ){
                            fNRebin[i]  = fBinsPi013TeVLowBPHOSPtRebin[i];
                        } else if ( modi == 2 )
                            fNRebin[i]  = fBinsPi013TeVLowBPCMEMCPtRebin[i];
                    }
                    nIterBGFit                  = 10;
                    fMaxYFracBGOverIntHist      = 60;
                    optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing5";
                    optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing7";
                    optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing3";
                }

            //*********************************************************************************************
            //********************************** Pi0 for pPb 5.023TeV**************************************
            //*********************************************************************************************
          } else if( energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0|| energy.CompareTo("pPb_5.023TeVRun2") == 0) {
                if (directPhoton.Contains("directPhoton") ){
                    fStartPtBin                 = GetStartBin(directPhoton, energy, modi, specialTrigg, centrality);
                    Int_t maxPtBinTheo          = GetBinning( fBinsPt, maxPtBinAvail, "Gamma", energy, modi, specialTrigg, isDCA, centrality );
                    if (fNBinsPt > maxPtBinAvail) {
                        cout << "**************************************************************************************************************************************" << endl;
                        cout << "********************** ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, **********************************" << endl;
                        cout << "You have chosen "<< fNBinsPt << " bins, which is more than the maximal " << maxPtBinAvail << " bins, this is not possible, it will be reduced to " << maxPtBinAvail << endl;
                        cout << "**************************************************************************************************************************************" << endl;
                        fNBinsPt    = maxPtBinAvail;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                    for (Int_t i = 0; i < fNBinsPt; i++) {
                        if (modi == 0){
                            if(energy.CompareTo("pPb_5.023TeVRun2") == 0){
                                if(!centrality.CompareTo("0-100%"))
                                    fNRebin[i]  = fBinsDirGammapPb5TeVPCMPtRebin[i];
                                else if(!centrality.CompareTo("0-1%") || !centrality.CompareTo("0-2%"))
                                    fNRebin[i]  = fBinsDirGammapPb5TeVR2Cent3PCMPtRebin[i];
                                else if (!centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") || !centrality.CompareTo("0-20%") || !centrality.CompareTo("20-40%") || !centrality.CompareTo("40-60%") || !centrality.CompareTo("60-100%"))
                                    fNRebin[i]  = fBinsPi0pPb5TeVPCMR2CentPtRebin[i];
                                else
                                    fNRebin[i]  = fBinsPi0pPb5TeVPCMR2CentPtRebin[i];
                            }else{
                                if(!centrality.CompareTo("0-100%"))
                                    fNRebin[i]  = fBinsDirGammapPb5TeVPCMPtRebin[i];
                                else
                                    fNRebin[i]  = fBinsDirGammapPb5TeVCentPCMPtRebin[i];
                            }
                        } else if (modi == 2 || modi == 3){
                            if(!centrality.CompareTo("0-100%"))
                                fNRebin[i]  = fBinsDirGammapPb5TeVPCMEMCPtRebin[i];
                            else
                                fNRebin[i]  = fBinsDirGammapPb5TeVCentPCMEMCPtRebin[i];
                        } else if (modi == 4){
                            fNRebin[i]  = fBinsDirGammapPb5TeVEMCPtRebin[i];
                        } else {
                            fNRebin[i]  = fBinsDirGammapPb5TeVPtRebin[i];
                        }
                    }
                    if(energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0){
                        optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                        optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                        optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing1";
                        nIterBGFit                  = 13;
                        fMaxYFracBGOverIntHist      = 50;
                    }else{
                        optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                        optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                        optionBGSmoothingVar2       = "noSmoothing";
                        nIterBGFit                  = 13;
                        fMaxYFracBGOverIntHist      = 50;
                    }
                } else {
                    fStartPtBin                 = GetStartBin("Pi0", energy, modi, specialTrigg, centrality);
                    Int_t maxPtBinTheo          = GetBinning( fBinsPt, maxPtBinAvail, "Pi0", energy, modi, specialTrigg, isDCA, centrality );
                    if (fNBinsPt > maxPtBinAvail) {
                      cout << "**************************************************************************************************************************************" << endl;
                      cout << "********************** ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, **********************************" << endl;
                      cout << "You have chosen "<< fNBinsPt << " bins, which is more than the maximal " << maxPtBinAvail << " bins, this is not possible, it will be reduced to " << maxPtBinAvail << endl;
                      cout << "**************************************************************************************************************************************" << endl;
                      fNBinsPt    = maxPtBinAvail;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                    for (Int_t i = 0; i < fNBinsPt; i++) {
                        if (modi == 0 ){
                          if (!energy.CompareTo("pPb_5.023TeVCent")){
                              fNRebin[i]          = fBinsPi0pPb5TeVCentPCMPtRebin[i];
                          } else if (!energy.CompareTo("pPb_5.023TeV")){
                            if((centrality.Contains("0-100%") && !centrality.Contains("60-100%")))                                                       // MB pi0 for PCM run 1
                              fNRebin[i]          = fBinsPi0pPb5TeVPtRebin[i];
                            else
                              fNRebin[i]          = fBinsPi0pPb5TeVCentPCMPtRebin[i];
                          } else if (!energy.CompareTo("pPb_5.023TeVRun2")){
                            if(centrality.Contains("0-5%") || centrality.Contains("5-10%") || centrality.Contains("60-100%"))
                              fNRebin[i]          = fBinsPi0pPb5TeVPCMR2CentPtRebin[i];
                            else if(centrality.Contains("0-1%") || centrality.Contains("0-2%"))
                              fNRebin[i]          = fBinsPi0pPb5TeVPCMR2Cent2PtRebin[i];
                            else
                              fNRebin[i]          = fBinsPi0pPb5TeVPCMR2PtRebin[i];                                          // MB pi0 for PCM run 2
                          }
                        } else if (modi == 1 )
                            fNRebin[i]          = fBinsPi0pPb5TeVDalitzPtRebin[i];
                        else if (   modi == 2 && energy.CompareTo("pPb_5.023TeV") == 0 && specialTrigg == 0)               // MB pi0 for PCM-EMC run 1
                            if((centrality.Contains("0-100%") && !centrality.Contains("60-100%")))                                                       // MB pi0 for PCM run 1
                              fNRebin[i]          = fBinsPi0pPb5TeVPCMEMCPtRebin[i];
                            else
                              fNRebin[i]          = fBinsPi0pPb5TeVCentPCMEMCPtRebin[i];
                        else if (modi == 2 && energy.CompareTo("pPb_5.023TeVCent") == 0 && specialTrigg == 0)                   // cent dependent pi0 for PCM-EMC run 1
                            fNRebin[i]          = fBinsPi0pPb5TeVCentPCMEMCPtRebin[i];
                        else if (   modi == 2 && energy.CompareTo("pPb_5.023TeVRun2") == 0 && specialTrigg == 0 &&          // cent dependent pi0 for PCM-EMC run 2: 0020, 0010
                                    (centrality.Contains("0-20%") || centrality.Contains("0-10%") ) )
                            fNRebin[i]          = fBinsPi0pPb5TeVPCMEMCR2CentPtRebin[i];
                        else if (   modi == 2 && energy.CompareTo("pPb_5.023TeVRun2") == 0 && specialTrigg == 0 &&          // cent dependent pi0 for PCM-EMC run 2: 0005, 0510, 60100
                            ( centrality.Contains("0-5%") || centrality.Contains("5-10%") || centrality.Contains("60-100%") ) )
                            fNRebin[i]          = fBinsPi0pPb5TeVPCMEMCR2SCentPtRebin[i];
                        else if ( (modi == 2 || modi == 13 ) && energy.CompareTo("pPb_5.023TeVRun2") == 0 && specialTrigg == 0)               // MB pi0 for PCM-EMC run 2
                            fNRebin[i]          = fBinsPi0pPb5TeVPCMEMCR2PtRebin[i];
                        else if (modi == 2 && specialTrigg == 1 )                                                           // MB pi0 for PCM-EMC run 1 - triggered EMC7
                            fNRebin[i]          = fBinsPi0pPb5TeVPCMEMCTrigEMC7PtRebin[i];
                        else if (modi == 2 && specialTrigg == 2 )                                                           // MB pi0 for PCM-EMC run 1 - triggered EG2
                            fNRebin[i]          = fBinsPi0pPb5TeVPCMEMCTrigEG2PtRebin[i];
                        else if (modi == 2 && specialTrigg == 3 )                                                           // MB pi0 for PCM-EMC run 1 - triggered EMC7, EG2
                            fNRebin[i]          = fBinsPi0pPb5TeVPCMEMCTrigEG1PtRebin[i];
                        else if (modi == 3 && energy.CompareTo("pPb_5.023TeV") == 0)       // MB pi0 for PCM-PHOS run 1
                            if((centrality.Contains("0-100%") && !centrality.Contains("60-100%")))                                                       // MB pi0 for PCM run 1
                              fNRebin[i]          = fBinsPi0pPb5TeVPCMPHOSPtRebin[i];
                            else
                              fNRebin[i]          = fBinsPi0pPb5TeVCentPCMPHOSPtRebin[i];
                        else if (modi == 3 && energy.CompareTo("pPb_5.023TeVCent") == 0 )                                       // cent dependent pi0 for PCM-PHOS run 1
                            fNRebin[i]          = fBinsPi0pPb5TeVCentPCMPHOSPtRebin[i];
                        else if (modi == 3 && energy.CompareTo("pPb_5.023TeVRun2") == 0 && centrality.Contains("0-100%"))   // MB pi0 for PCM-PHOS run 2
                            fNRebin[i]          = fBinsPi0pPb5TeVPCMPHOSR2PtRebin[i];
                        else if (modi == 3 && energy.CompareTo("pPb_5.023TeVRun2") == 0 && centrality.Contains("0-20%"))    // cent dep pi0 for PCM-PHOS run 2: 0020
                            fNRebin[i]          = fBinsPi0pPb5TeVPCMPHOSR2CentPtRebin[i];
                        else if (modi == 3 && energy.CompareTo("pPb_5.023TeVRun2") == 0 &&                                  // cent dep pi0 for PCM-PHOS run 2: 0005, 0510
                            ( centrality.Contains("0-5%") || centrality.Contains("5-10%")) )
                            fNRebin[i]          = fBinsPi0pPb5TeVPCMPHOSR2SCentPtRebin[i];
                        else if (modi == 3 && energy.CompareTo("pPb_5.023TeVRun2") == 0 && centrality.Contains("60-100%"))  // cent dep pi0 for PCM-PHOS run 2: 60100
                            fNRebin[i]          = fBinsPi0pPb5TeVPCMPHOSR2PerPtRebin[i];
                        else if (modi == 3 && energy.CompareTo("pPb_5.023TeVRun2") == 0 )                                   // cent dep pi0 for PCM-PHOS run 2
                            fNRebin[i]          = fBinsPi0pPb5TeVPCMPHOSR2PtRebin[i];
                        else if (modi == 4 && energy.CompareTo("pPb_5.023TeV") == 0 && specialTrigg == 1 )                  // EMC7 EMC run 1
                            fNRebin[i]          = fBinsPi0pPb5TeVEMCTrigEMC7PtRebin[i];
                        else if (modi == 4 && energy.CompareTo("pPb_5.023TeV") == 0 && specialTrigg == 2 )                  // EG2 EMC run 1
                            fNRebin[i]          = fBinsPi0pPb5TeVEMCTrigEG2PtRebin[i];
                        else if (modi == 4 && energy.CompareTo("pPb_5.023TeV") == 0 && specialTrigg == 3 )                  // EG1 EMC run 1
                            fNRebin[i]          = fBinsPi0pPb5TeVEMCTrigEG1PtRebin[i];
                        else if (modi == 4 && energy.CompareTo("pPb_5.023TeV") == 0 ) // MB pi0 for EMC run 1
                            if((centrality.Contains("0-100%") && !centrality.Contains("60-100%")))                                                       // MB pi0 for PCM run 1
                              fNRebin[i]          = fBinsPi0pPb5TeVEMCPtRebin[i];
                            else
                              fNRebin[i]          = fBinsPi0pPb5TeVCentEMCPtRebin[i];
                        else if (modi == 4 && energy.CompareTo("pPb_5.023TeVCent") == 0 )                                       // cent dependent pi0 for EMC run 1
                            fNRebin[i]          = fBinsPi0pPb5TeVCentEMCPtRebin[i];
                        else if (modi == 4 && energy.CompareTo("pPb_5.023TeVRun2") == 0 )                                   // MB pi0 for EMC run 2
                            fNRebin[i]          = fBinsPi0pPb5TeVEMCR2PtRebin[i];
                        else if (modi == 5 && energy.CompareTo("pPb_5.023TeV") == 0 && centrality.CompareTo("0-100%") == 0) // MB pi0 for PHOS run 1
                            fNRebin[i]          = fBinsPi0pPb5TeVPHOSPtRebin[i];
                        else if (modi == 5 && energy.CompareTo("pPb_5.023TeV") == 0 )                                       // cent dependent pi0  for PHOS run 1
                            fNRebin[i]          = fBinsPi0pPb5TeVCentPHOSPtRebin[i];
                        else if (modi == 5 && energy.CompareTo("pPb_5.023TeVRun2") == 0 && centrality.Contains("0-20%") )  // cent dependent pi0 for PHOS run 2: 0020
                            fNRebin[i]          = fBinsPi0pPb5TeVPHOSR2CentPtRebin[i];
                        else if (modi == 5 && energy.CompareTo("pPb_5.023TeVRun2") == 0 &&
                          (centrality.Contains("0-5%") || centrality.Contains("5-10%")))                                   // cent dependent pi0 for PHOS run 2: 0005, 0510
                            fNRebin[i]          = fBinsPi0pPb5TeVPHOSR2SCentPtRebin[i];
                        else if (modi == 5 && energy.CompareTo("pPb_5.023TeVRun2") == 0 && centrality.Contains("60-100%") ) // cent dependet pi0 for PHOS run 2: 60100
                            fNRebin[i]          = fBinsPi0pPb5TeVPHOSR2TCentPtRebin[i];
                        else if (modi == 5 && energy.CompareTo("pPb_5.023TeVRun2") == 0 )                                   // MB pi0 for PHOS run 2
                            fNRebin[i]          = fBinsPi0pPb5TeVPHOSR2PtRebin[i];
                        else if (modi == 6 )                                                                                // MB pi0 for Dalitz-PHOS
                            fNRebin[i]          = fBinsPi0pPb5TeVEMCDalitzPtRebin[i];
                        else if (modi == 7 )                                                                                // MB pi0 for Dalitz-PHOS
                            fNRebin[i]          = fBinsPi0pPb5TeVEMCDalitzPtRebin[i];
                        else if (modi == 10 )                                                                               // MB pi0 for mEMC
                            fNRebin[i]          = fBinsPi0pPb5TeVmEMCPtRebin[i];
                    }
                    optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                    optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                    optionBGSmoothingVar2       = "noSmoothing";
                    nIterBGFit                  = 11;
                    if (modi == 0){
                      if( !energy.CompareTo("pPb_5.023TeVRun2") ){
                        if(centrality.Contains("60-100%"))
                            nIterBGFit                  = 7;
                        else
                            nIterBGFit                  = 8;
                      }else if( !energy.CompareTo("pPb_5.023TeV") ){
                        if(centrality.Contains("0-20%") || centrality.Contains("40-60%"))
                          nIterBGFit                = 6;
                        else if(centrality.Contains("20-40%"))
                          nIterBGFit                = 7;
                        else if(centrality.Contains("60-100%"))
                          nIterBGFit                = 8;
                        else if(centrality.CompareTo(""))
                          nIterBGFit                = 7;
                      }
                    }
                    fMaxYFracBGOverIntHist      = 50;
                }
            //*********************************************************************************************
            //********************************** Pi0 for pPb 8TeV**************************************
            //*********************************************************************************************
           } else if( energy.CompareTo("pPb_8TeV") == 0) {
                if (directPhoton.Contains("directPhoton") ){
                    fStartPtBin     = 1;
                    if (modi == 2 && directPhoton.CompareTo("directPhoton") == 0){
                        fStartPtBin     = 1;
                    } else if (modi == 2 && directPhoton.CompareTo("directPhotonTagging") == 0){
                        fStartPtBin     = 8;
                    } else if (modi == 4){
                        fStartPtBin     = 8;
                    }

                    if (specialTrigg == 0 && ( modi == 2 || modi == 4 ) && fNBinsPt > 28 ){
                        cout << "You have chosen to have more than 28 bins, this is not possible, it will be reduced to 28 for calo analysis" << endl;
                        fNBinsPt    = 28;
                    } else if ( !( modi == 2 || modi == 4 ) && fNBinsPt > 25) {
                        cout << "You have chosen Direct Photon Plots and more than 25 bins, this is not possible, it will be reduced to 25 bins." << endl;
                        fNBinsPt    = 23;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        if (modi == 0){
                            fBinsPt[i]  = fBinsDirGammapPb8TeVPt[i];
                            if (i < fNBinsPt+1)
                                fNRebin[i]  = fBinsDirGammapPb8TeVPtRebin[i];
                        } else if (modi == 2){
                            fBinsPt[i]  = fBinsDirGammapPb8TeVPCMEMCPt[i];
                            if (i < fNBinsPt+1)
                                fNRebin[i]  = fBinsDirGammapPb8TeVPCMEMCPtRebin[i];
                        } else if (modi == 4){
                            fBinsPt[i]  = fBinsDirGammapPb8TeVPt[i];
                            if (i < fNBinsPt+1)
                                fNRebin[i]  = fBinsDirGammapPb8TeVPtRebin[i];
                        } else {
                            fBinsPt[i]  = fBinsDirGammapPb8TeVPt[i];
                            if (i < fNBinsPt+1)
                                fNRebin[i]  = fBinsDirGammapPb8TeVPtRebin[i];
                        }
                    }
                    optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                    optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                    optionBGSmoothingVar2       = "noSmoothing";
                    nIterBGFit                  = 13;
                    fMaxYFracBGOverIntHist      = 20;
                } else {
                    fStartPtBin                 = GetStartBin("Pi0", energy, modi, specialTrigg);
                    if (fNBinsPt > 16 && isDCA) {
                        cout << "You have chosen to have more than 16 bins, this is not possible, it will be reduced to 16" << endl;
                        fNBinsPt    = 16;
                    } else if (specialTrigg == 0 && ( modi == 2 || modi == 4 ) && fNBinsPt > 36 ){
                        cout << "You have chosen to have more than 36 bins, this is not possible, it will be reduced to 36 for calo analysis" << endl;
                        fNBinsPt    = 36;
                    } else if (fNBinsPt > 39 && specialTrigg == 0 && modi != 10){
                        cout << "You have chosen to have more than 39 bins, this is not possible, it will be reduced to 39 for conv analysis" << endl;
                        fNBinsPt    = 37;
                    } else if (fNBinsPt > 60 ){
                        cout << "You have chosen to have more than 60 bins, this is not possible, it will be reduced to 60 for calo analysis" << endl;
                        fNBinsPt    = 60;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        if ( modi == 10){
                            fBinsPt[i]  = fBinsPi0pPb8TeVPtmEMC[i];
                            if (i < fNBinsPt+1)
                                fNRebin[i]  = fBinsPi0pPb8TeVmEMCPtRebin[i];
                        } else {
                            if (isDCA)
                                fBinsPt[i]  = fBinsPi0pPb8TeVPtDCA[i];
                            else
                                fBinsPt[i]  = fBinsPi0pPb8TeVPt[i];
                            if (modi == 2 || modi == 4)
                                fBinsPt[i]  = fBinsPi0pPb8TeVEMCPt[i];
                            else if (modi == 3 || modi == 5)
                                fBinsPt[i]  = fBinsPi0pPb8TeVPHOSPt[i];
                            else if ( modi == 6 )
                                fBinsPt[i]  = fBinsPi0pPb8TeVEMCDalitzPt[i];
                            else if ( modi == 1 )
                                fBinsPt[i]  = fBinsPi0pPb8TeVDalitzPt[i];

                            if (i < fNBinsPt+1){
                                fNRebin[i]                         = fBinsPi0pPb8TeVPtRebin[i];
                                if (modi == 1 ) fNRebin[i]         = fBinsPi0pPb8TeVDalitzPtRebin[i];
                                if (modi == 2 ) fNRebin[i]         = fBinsPi0pPb8TeVPCMEMCPtRebin[i];
                                if (modi == 3 ) fNRebin[i]         = fBinsPi0pPb8TeVPCMPHOSPtRebin[i];
                                if (modi == 4 ) fNRebin[i]         = fBinsPi0pPb8TeVEMCPtRebin[i];
                                if (modi == 5 ) fNRebin[i]         = fBinsPi0pPb8TeVPHOSPtRebin[i];
                                if (modi == 6 ) fNRebin[i]         = fBinsPi0pPb8TeVEMCDalitzPtRebin[i];
                                if (modi == 7 ) fNRebin[i]         = fBinsPi0pPb8TeVEMCDalitzPtRebin[i];
                            }
                        }
                    }
                    optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                    optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                    optionBGSmoothingVar2       = "noSmoothing";
                    nIterBGFit                  = 11;
                    fMaxYFracBGOverIntHist      = 20;
                }
            //*********************************************************************************************
            //********************************** Pi0 for PbPb 2.76TeV**************************************
            //*********************************************************************************************
            } else if( energy.CompareTo("PbPb_2.76TeV") == 0) {
                if (directPhoton.CompareTo("directPhoton") == 0){
                    fStartPtBin     = 3;
                    if (fNBinsPt > 22) {
                        cout << "You have chosen Direct Photon Plots and more than 24 bins, this is not possible, it will be reduced to 24 bins." << endl;
                        fNBinsPt    = 22;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        if(fNBinsPt==22){
                            fBinsPt[i]  = fBinsDirGammaPbPb2760GeVPtLHC11h[i];
                            if (i < fNBinsPt+1){
                                if(centrality.CompareTo("20-40%")==0 || centrality.CompareTo("20-50%")==0) fNRebin[i] = fBinsDirGammaPbPb2760GeVPtLHC11hSemicRebin[i];
                                else fNRebin[i] = fBinsDirGammaPbPb2760GeVPtLHC11hRebin[i];
                            }
                        } else {
                            fBinsPt[i]  = fBinsDirGammaPbPb2760GeVPtLHC11hVar2[i];
                            if (i < fNBinsPt+1){
                                if(centrality.CompareTo("20-40%")==0 || centrality.CompareTo("20-50%")==0) fNRebin[i] = fBinsDirGammaPbPb2760GeVPtLHC11hSemicRebinVar2[i];
                                else fNRebin[i] = fBinsDirGammaPbPb2760GeVPtLHC11hRebinVar2[i];
                            }
                        }
                    }
                } else {
                    fStartPtBin     = GetStartBin("Pi0", energy, modi, specialTrigg, clusterCutSelection(GetClusterMinEnergyCutPosition(clusterCutSelection),1));
                    if (fNBinsPt > 15 && isDCA) {
                        cout << "You have chosen to have more than 15 bins, this is not possible, it will be reduced to 15" << endl;
                        fNBinsPt    = 15;
                    } else if (fNBinsPt > 26) {
                        cout << "You have chosen to have more than 26 bins, this is not possible, it will be reduced to 24" << endl;
                        fNBinsPt    = 26;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        if (isDCA) {
                            if (!centDCA.CompareTo("60-80%") || !centDCA.CompareTo("70-80%") || !centDCA.CompareTo("75-90%") ){
                                fBinsPt[i]          = fBinsPi0PbPb2760GeVPtDCAPer[i];
                            } else {
                                fBinsPt[i]          = fBinsPi0PbPb2760GeVPtDCA[i];
                            }
                        } else fBinsPt[i]  = fBinsPi0PbPb2760GeVPtLHC11h[i];
                        if (modi == 0){
                            if (i < fNBinsPt+1){
                                if(centrality.CompareTo("20-40%")==0 || centrality.CompareTo("20-50%")==0) fNRebin[i] = fBinsPi0PbPb2760GeVPtLHC11hSemicRebin[i];
                                else fNRebin[i] = fBinsPi0PbPb2760GeVPtLHC11hRebin[i];
                            }
                        } else {
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0PbPb2760GeVPtLHC11hPCMEMCRebin[i];
                        }
                    }
                    optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing5";
                    optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing3";
                    optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing7";
                    if (!centDCA.CompareTo("60-80%")){
                        nIterBGFit                  = 15;
                        fMaxYFracBGOverIntHist      = 15;
                    } else if ( (!centDCA.CompareTo("75-90%")) || (!centDCA.CompareTo("70-80%"))){
                        optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing9";
                        optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing7";
                        optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing11";
                        nIterBGFit                  = 14;
                        fMaxYFracBGOverIntHist      = 50;
                    } else if (!centDCA.CompareTo("70-80%")){
                        optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing9";
                        optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing7";
                        optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing11";
                        nIterBGFit                  = 14;
                        fMaxYFracBGOverIntHist      = 50;
                    } else if (!centDCA.CompareTo("60-70%")){
                        optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing7";
                        optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                        optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing9";
                        nIterBGFit                  = 17;
                        fMaxYFracBGOverIntHist      = 15;
                    } else if (!centDCA.CompareTo("50-60%")){
                        nIterBGFit                  = 14;
                        fMaxYFracBGOverIntHist      = 12;
                    } else if (!centDCA.CompareTo("40-60%")){
                        nIterBGFit                  = 17;
                        fMaxYFracBGOverIntHist      = 10;
                    } else if (!centDCA.CompareTo("40-50%")) {
                        nIterBGFit                  = 16;
                        fMaxYFracBGOverIntHist      = 12;
                    } else if ( (!centDCA.CompareTo("30-50%")) || (!centDCA.CompareTo("30-40%")) ) {
                        nIterBGFit                  = 18;
                        fMaxYFracBGOverIntHist      = 12;
                    } else if (!centDCA.CompareTo("20-50%")) {
                        nIterBGFit                  = 17;
                        fMaxYFracBGOverIntHist      = 12;
                    } else if (!centDCA.CompareTo("20-40%")){
                        nIterBGFit              = 17;
                        fMaxYFracBGOverIntHist  = 8;
                    } else if (!centDCA.CompareTo("20-30%")) {
                        nIterBGFit                  = 18;
                        fMaxYFracBGOverIntHist      = 12;
                    } else {
                        fMaxYFracBGOverIntHist      = 4;
                        nIterBGFit                  = 21;
                    }
                }
            //*********************************************************************************************
            //********************************** Pi0 for PbPb 5.02TeV**************************************
            //*********************************************************************************************
            } else if( energy.CompareTo("PbPb_5.02TeV") == 0) {
                if (directPhoton.CompareTo("directPhoton") == 0){
                    fStartPtBin     = 1;
                    if (fNBinsPt > 19) {
                        cout << "You have chosen Direct Photon Plots and more than 19 bins, this is not possible, it will be reduced to 19 bins." << endl;
                        fNBinsPt    = 19;
                    }
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        fBinsPt[i]  = fBinsDirGammaPbPb5TeVPt[i];
                        if (i < fNBinsPt+1) fNRebin[i] = fBinsDirGammaPbPb5TeVPtRebin[i];
                    }
                } else{
                    fStartPtBin     = GetStartBin("Pi0", energy, modi, specialTrigg, centrality);
                    if (fNBinsPt > 12 && isDCA) {
                        cout << "You have chosen to have more than 12 bins, this is not possible, it will be reduced to 12" << endl;
                        fNBinsPt        = 12;
                    }
                    cout << "calling GetOptimumNColumnsAndRows " <<  fNBinsPt << " , " << fStartPtBin << endl;
                    GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                    for (Int_t i = 0; i < fNBinsPt+1; i++) {
                        if (isDCA){
                            fBinsPt[i]  = fBinsPi0PbPb5TeVPtDCA[i];
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0PbPb5TeVPtRebin[i];
                        }else{
                          if (modi == 2 || modi == 3){
                            fBinsPt[i]  = fBinsPi0PbPb5TeVPCMEMCPt[i];
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0PbPb5TeVPCMEMCPtRebin[i];
                          }else if (modi == 4){
                            fBinsPt[i]  = fBinsPi0PbPb5TeVEMCPt[i];
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0PbPb5TeVEMCPtRebin[i];
                          }else{
                            fBinsPt[i]  = fBinsPi0PbPb5TeVPt[i];
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0PbPb5TeVPtRebin[i];
                          }
                        }
                    }
                }
            //*********************************************************************************************
            //********************************** Pi0 for XeXe 5.44TeV**************************************
            //*********************************************************************************************
            } else if( energy.CompareTo("XeXe_5.44TeV") == 0) {
                fStartPtBin     = GetStartBin("Pi0", energy, modi, -1, centrality);
                Int_t maxPtBinTheo          = GetBinning( fBinsPt, maxPtBinAvail, "Pi0", energy, modi, specialTrigg, isDCA, centrality );
                if (fNBinsPt > maxPtBinAvail) {
                    cout << "**************************************************************************************************************************************" << endl;
                    cout << "********************** ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, **********************************" << endl;
                    cout << "You have chosen "<< fNBinsPt << " bins, which is more than the maximal " << maxPtBinAvail << " bins, this is not possible, it will be reduced to " << maxPtBinAvail << endl;
                    cout << "**************************************************************************************************************************************" << endl;
                    fNBinsPt    = maxPtBinAvail;
                }
                GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                for (Int_t i = 0; i < fNBinsPt+1; i++) {
                    if (i < fNBinsPt+1){
                        fNRebin[i]            = fBinsPi0XeXe5440GeVPtRebin[i];
                        if (!centrality.CompareTo("0-90%") || !centrality.CompareTo("0-80%") ){
                            if (modi == 0)      fNRebin[i]        = fBinsPi0XeXe5440GeVPtPCMRebin[i];
                            else if (modi == 2) fNRebin[i]        = fBinsPi0XeXe5440GeVPtPCMEMCRebin[i];
                            else if (modi == 3) fNRebin[i]        = fBinsPi0XeXe5440GeVPtPCMPHOSRebin[i];
                            else if (modi == 4) fNRebin[i]        = fBinsPi0XeXe5440GeVPtEMCRebin[i];
                            else if (modi == 5) fNRebin[i]        = fBinsPi0XeXe5440GeVPtPHOSRebin[i];
                            else                fNRebin[i]        = fBinsPi0XeXe5440GeVPtRebin[i];
                        } else if (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-40%") || !centrality.CompareTo("20-40%")){
                            if (modi == 0)      fNRebin[i]        = fBinsPi0XeXe5440GeVPtCentPCMRebin[i];
                            else if (modi == 2) fNRebin[i]        = fBinsPi0XeXe5440GeVPtCentPCMEMCRebin[i];
                            else if (modi == 3) fNRebin[i]        = fBinsPi0XeXe5440GeVPtCentPCMPHOSRebin[i];
                            else if (modi == 4) fNRebin[i]        = fBinsPi0XeXe5440GeVPtCentEMCRebin[i];
                            else if (modi == 5) fNRebin[i]        = fBinsPi0XeXe5440GeVPtCentPHOSRebin[i];
                            else                fNRebin[i]        = fBinsPi0XeXe5440GeVPtRebinCent[i];
                        } else if ( !centrality.CompareTo("20-40%")){
                            if (modi == 0)      fNRebin[i]        = fBinsPi0XeXe5440GeVPtCentPCMRebin[i];
                            else if (modi == 2) fNRebin[i]        = fBinsPi0XeXe5440GeVPtCentPCMEMCRebin[i];
                            else if (modi == 3) fNRebin[i]        = fBinsPi0XeXe5440GeVPtSemiPCMPHOSRebin[i];
                            else if (modi == 4) fNRebin[i]        = fBinsPi0XeXe5440GeVPtCentEMCRebin[i];
                            else if (modi == 5) fNRebin[i]        = fBinsPi0XeXe5440GeVPtSemiPHOSRebin[i];
                            else                fNRebin[i]        = fBinsPi0XeXe5440GeVPtRebinCent[i];
                        } else if (!centrality.CompareTo("40-80%") ){
                            if (modi == 0)      fNRebin[i]        = fBinsPi0XeXe5440GeVPtPerPCMRebin[i];
                            else if (modi == 2) fNRebin[i]        = fBinsPi0XeXe5440GeVPtPerPCMEMCRebin[i];
                            else if (modi == 3) fNRebin[i]        = fBinsPi0XeXe5440GeVPtPerPCMPHOSRebin[i];
                            else if (modi == 4) fNRebin[i]        = fBinsPi0XeXe5440GeVPtPerEMCRebin[i];
                            else if (modi == 5) fNRebin[i]        = fBinsPi0XeXe5440GeVPtPerPHOSRebin[i];
                            else                fNRebin[i]        = fBinsPi0XeXe5440GeVPtRebinCent[i];
                        } else {
                            fNRebin[i]        = fBinsPi0XeXe5440GeVPtRebinCent[i];
                        }
                    }
                }

                optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing5";
                optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing3";
                optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing7";
                nIterBGFit                  = 15;
                fMaxYFracBGOverIntHist      = 15;
            }

        //*************************************************************************************************
        //********************************** Binning for Eta **********************************************
        //*************************************************************************************************
        } else if (setPi0.CompareTo("Eta") == 0 || setPi0.CompareTo("Pi0EtaBinning") == 0){
            fNBinsPt                = numberOfBins;
            fBinsPt                 = new Double_t[150];
            fNRebin                 = new Int_t[149];
            //*********************************************************************************************
            //********************************** Eta for pp 0.9 TeV****************************************
            //*********************************************************************************************
            if (energy.CompareTo("900GeV") == 0) {
                fStartPtBin         = 1;
                if(modi == 4) fStartPtBin         = 3;
                if (modi == 2){
                    if (fNBinsPt > 4) {
                        cout << "You have chosen to have more than 4 bins for Eta, this is not possible, it will be reduced to 4" << endl;
                        fNBinsPt        = 4;
                    }
                } else {
                    if (fNBinsPt > 3) {
                        cout << "You have chosen to have more than 3 bins for Eta, this is not possible, it will be reduced to 3" << endl;
                        fNBinsPt        = 3;
                    }
                }
                GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                for (Int_t i = 0; i < fNBinsPt+1; i++) {
                if( modi == 2){
                    fBinsPt[i]      = fBinsEta900GeVPCMEMCPt[i];
                    if (i < fNBinsPt+1)
                            fNRebin[i]  = fBinsEta900GeVPCMEMCPtRebin[i];
                    }else{
                        fBinsPt[i]      = fBinsEta900GeVPt[i];
                        if (i < fNBinsPt+1)
                            fNRebin[i]  = fBinsEta900GeVPtRebin[i];
                    }
                }
                nIterBGFit          = 13;
                fMaxYFracBGOverIntHist = 20;
            //*********************************************************************************************
            //********************************** Eta for pp 2.76TeV****************************************
            //*********************************************************************************************
            } else if (energy.CompareTo("2.76TeV") == 0) {
                fStartPtBin         = 1;
                if ( modi == 3)
                    fStartPtBin     = 2;
                else if (modi == 2 && specialTrigg == 0) // MB, PCM-EMC
                    fStartPtBin     = 2;
                else if (modi == 2 && specialTrigg == 5) { // INT7, PCM-EMC
                    fStartPtBin     = 4;
                } else if (modi == 2 && specialTrigg == 1) // EMC7, PCM-EMC
                    fStartPtBin     = 4;
                else if (modi == 2 && specialTrigg == 2) // L1 G2, PCM-EMC
                    fStartPtBin     = 6;
                else if (modi == 2 && specialTrigg == 3) // L1 G1, PCM-EMC
                    fStartPtBin     = 7;
                else if (modi == 2 && specialTrigg == 4) // EMC1, PCM-EMC
                    fStartPtBin     = 4;
                else if (modi == 4 && specialTrigg == 1)
                    fStartPtBin     = 5;
                else if (modi == 4 && specialTrigg == 2)
                    fStartPtBin     = 6;
                else if (modi == 4 && specialTrigg == 3)
                    fStartPtBin     = 7;
                else if (modi == 4 && specialTrigg == 4)
                    fStartPtBin     = 6;
                else if (modi == 4 )
                    fStartPtBin     = 4;
                else if (modi == 5 )
                    fStartPtBin     = 3;

                if (fNBinsPt > 14 && isDCA) {
                    cout << "You have chosen to have more than 15 DCA bins for Eta, this is not possible, it will be reduced to 15" << endl;
                    fNBinsPt        = 14;
                } else if (fNBinsPt > 7 && (modi == 0 || modi == 1) && specialTrigg < 1) {
                    cout << "You have chosen to have more than 7 bins for Eta, this is not possible, it will be reduced to 7" << endl;
                    fNBinsPt        = 7;
                } else if (fNBinsPt > 13 && (modi == 2 || modi == 3 || modi == 4 || modi == 0)){
                    cout << "You have chosen to have more than 13 bins for Eta, this is not possible, it will be reduced to 13" << endl;
                    fNBinsPt        = 13;
                }
                GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                for (Int_t i = 0; i < fNBinsPt+1; i++) {
                    if ( ( modi == 2 && specialTrigg == 0) ||
                        ( modi == 4 && specialTrigg == 0) ){
                        fBinsPt[i]  = fBinsEta2760GeVPt[i];
                        if (setPi0.CompareTo("Eta") == 0){
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsEta2760GeVPCMEMCPtRebin[i];
                        } else {
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0EtaBinning2760GeVPCMEMCPtRebin[i];
                        }
                    } else if ( (modi == 2 && specialTrigg == 5)||
                                (modi == 4 && specialTrigg == 5)){
                        fBinsPt[i]  = fBinsEta2760GeVPt[i];
                        if (setPi0.CompareTo("Eta") == 0){
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsEta2760GeVPCMEMCPtTrigINT7Rebin[i];
                        } else {
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0EtaBinning2760GeVPCMEMCPtRebin[i];
                        }
                    } else if ( (modi == 2 && specialTrigg == 4 ) ||
                                (modi == 4 && specialTrigg == 4 ) ||
                                (modi == 4 && specialTrigg == 2 ) ||
                                (modi == 4 && specialTrigg == 3 )
                            ){
                        fBinsPt[i]  = fBinsEta2760GeVPtTrig11a[i];
                        if (setPi0.CompareTo("Eta") == 0){
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsEta2760GeVPCMEMCPtTrig11aRebin[i];
                        } else {
                            if (i < fNBinsPt+1){
                                fNRebin[i] = fBinsPi0EtaBinning2760GeVPCMEMCPtTrig11aRebin[i];
                            }
                        }
                    } else if ( (modi == 2 && specialTrigg == 3 ) ||
                                (modi == 0 && specialTrigg > 0 )
                            ){
                        fBinsPt[i]  = fBinsEta2760GeVPtTrig11a[i];
                        if (setPi0.CompareTo("Eta") == 0){
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsEta2760GeVPCMEMCPtRebin[i];
                        } else {
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0EtaBinning2760GeVPCMEMCPtRebin[i];
                        }
                    } else if ( modi == 2 && specialTrigg == 2 ){
                        fBinsPt[i]  = fBinsEta2760GeVPtTrig11a[i];
                        if (setPi0.CompareTo("Eta") == 0){
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsEta2760GeVPCMEMCPtEG2Rebin[i];
                        } else {
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0EtaBinning2760GeVPCMEMCPtRebin[i];
                        }

                    } else {
                        if (isDCA) {
                            if (!setPi0.CompareTo("Pi0EtaBinning"))
                                fBinsPt[i]  = fBinsEta2760GeVPt[i];
                            else
                                fBinsPt[i]  = fBinsPi02760GeVPtDCA[i];
                        } else
                            fBinsPt[i]  = fBinsEta2760GeVPt[i];
                        if (setPi0.CompareTo("Eta") == 0){
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsEta2760GeVPtRebin[i];
                        } else {
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0EtaBinning2760GeVPtRebin[i];
                        }
                    }
                }
                if (!setPi0.CompareTo("Pi0EtaBinning"))
                    fMaxYFracBGOverIntHist      = 20;
                else
                    fMaxYFracBGOverIntHist      = 50;
                nIterBGFit                  = 13;
            //*********************************************************************************************
            //********************************** Eta for pp 5TeV*******************************************
            //*********************************************************************************************
            } else if (energy.CompareTo("5TeV") == 0 || energy.CompareTo("5TeV2017") == 0) {
              fStartPtBin                 = GetStartBin("Eta", energy, modi, specialTrigg, centrality);
              Int_t maxPtBinTheo          = GetBinning( fBinsPt, maxPtBinAvail, "Eta", energy, modi, specialTrigg, isDCA, centrality );
              if (fNBinsPt > maxPtBinTheo) {
                  cout << "**************************************************************************************************************************************" << endl;
                  cout << "********************** ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, **********************************" << endl;
                  cout << "You have chosen "<< fNBinsPt << " bins, which is more than the maximal " << maxPtBinTheo << " bins, this is not possible, it will be reduced to " << maxPtBinTheo << endl;
                  cout << "**************************************************************************************************************************************" << endl;
                  fNBinsPt    = maxPtBinTheo;
              }
              GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                for (Int_t i = 0; i < fNBinsPt; i++) {
                    if (setPi0.CompareTo("Eta") == 0){
                        if (i < fNBinsPt+1){
                            if ( modi == 0 ) {
                                if(energy.Contains("2017")){
                                    if(fNBinsPt<6)
                                        fNRebin[i]  = fBinsEta5TeV2017PCMforPbPbPtRebin[i];
                                    else if(fNBinsPt>6 && fNBinsPt<8)
                                        fNRebin[i]  = fBinsEta5TeVPtRebin[i];
                                    else if(fNBinsPt>8 && fNBinsPt<10)
                                        fNRebin[i]  = fBinsEta5TeV2017PCMCombinationPtRebin[i];
                                    else
                                        fNRebin[i]  = fBinsEta5TeV2017PtRebin[i];
                                } else {
                                    fNRebin[i]  = fBinsEta5TeVPtRebin[i];
                                }
                            } else if( (modi == 1) && (energy.Contains("2017"))){
                              fNRebin[i]  = fBinsEta5TeV2017DalitzPtRebin[i];
                            } else if( modi == 2 ){
                              if(energy.Contains("2017"))
                                fNRebin[i] = fBinsEta5TeV2017PCMEMCPtRebin[i];
                              else
                                fNRebin[i] = fBinsEta5TeVPCMEMCPtRebin[i];
                            } else if( modi == 3 ){
                              fNRebin[i]  = fBinsEta5TeV2017PtCombinationRebin[i];
                            } else if ( modi == 4 ){
                              if(energy.Contains("2017")){
                                fNRebin[i] = fBinsEta5TeV2017PtCombinationRebin[i];
                              } else {
                                if(specialTrigg == 1 || specialTrigg == 2 || specialTrigg == 3){
                                  fNRebin[i] = fBinsEta5TeVEMCPtRebinTrigger1[i];
                                }else{
                                  fNRebin[i] = fBinsEta5TeVEMCPtRebin[i];
                                }
                              }
                            } else if( modi == 13 ){
                              if(energy.Contains("2017"))
                                fNRebin[i] = fBinsEta5TeV2017PCMDCalPtRebin[i];
                              else
                                fNRebin[i] = fBinsEta5TeVPCMEMCPtRebin[i];
                            } else if( modi == 12 ){
                              if(energy.Contains("2017"))
                                fNRebin[i] = fBinsEta5TeV2017DMCPtRebin[i];
                            }
                      }
                    } else {
                        if (i < fNBinsPt+1) {
                            if( modi == 2 ){
                                fNRebin[i] = fBinsPi0EtaBinning5TeVPCMEMCPtRebin[i];
                            } else if( modi == 13 && energy.Contains("2017")){
                                fNRebin[i] = fBinsEta5TeV2017PCMDCalPtRebin[i];
                            } else if ( modi == 4 ){
                              if (specialTrigg == 1 || specialTrigg == 2 || specialTrigg == 3){
                                fNRebin[i]  = fBinsPi0EtaBinning5TeVPtRebinEMCTrigger1[i];
                              }else{
                                fNRebin[i]  = fBinsEta5TeV2017EMCPtRebin[i];
                              }
                            } else {
                              fNRebin[i]  = fBinsPi0EtaBinning5TeVPtRebin[i];
                            }
                        }
                    }
                }
                if(modi == 0 && energy.Contains("2017")){
                    nIterBGFit                  = 8;
                    fMaxYFracBGOverIntHist      = 70;
                    optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                    optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                    optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing7";
                } else {
                    nIterBGFit                  = 10;
                    fMaxYFracBGOverIntHist      = 60;
                    optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing5";
                    optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing6";
                    optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing7";
                }

            //*********************************************************************************************
            //********************************** Eta for pp 7TeV*******************************************
            //*********************************************************************************************
            } else if (energy.CompareTo("7TeV") == 0) {
                fStartPtBin         = 1;
                if(modi == 40 || modi == 41 || modi == 42 || modi == 44 || modi == 45) fStartPtBin     = GetStartBin("Eta","7TeV",modi);

                if (modi == 2 ) {
                    fStartPtBin     = 3;
                } else if (modi == 3 ) {
                    fStartPtBin     = 2;
                } else if (modi == 4) {
                    fStartPtBin     = 6;
                }

                if (fNBinsPt > 18) {
                    cout << "You have chosen to have more than 18 bins for Eta, this is not possible, it will be reduced to 18" << endl;
                    fNBinsPt        = 18;
                }
                GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                for (Int_t i = 0; i < fNBinsPt+2; i++) {
                    if (modi == 0){
                        fBinsPt[i]      = fBinsEta7TeVPt[i];
                        if (i < fNBinsPt+1){
                          if(!setPi0.CompareTo("Pi0EtaBinning"))
                            fNRebin[i]  = fBinsPi0EtaBinning7TeVPtRebin[i];
                          else
                            fNRebin[i]  = fBinsEta7TeVPtRebin[i];
                        }
                    } else if (modi == 1){
                        fBinsPt[i]      = fBinsEta7TeVDalitzPt[i];
                        if (i < fNBinsPt+1)
                            fNRebin[i]  = fBinsEta7TeVDalitzPtRebin[i];
                    } else if (modi == 2){
                        fBinsPt[i]      = fBinsEta7TeVPCMEMCPt[i];
                        if (i < fNBinsPt+1){
                          if(!setPi0.CompareTo("Pi0EtaBinning"))
                            fNRebin[i] = fBinsPi0EtaBinning7TeVPCMEMCPtRebin[i];
                          else
                            fNRebin[i] = fBinsEta7TeVPCMEMCPtRebin[i];
                        }
                    } else if (modi == 3){
                        fBinsPt[i]      = fBinsEta7TeVPCMPHOSPt[i];
                        if (i < fNBinsPt+1)
                            fNRebin[i]  = fBinsEta7TeVPCMPHOSPtRebin[i];
                    } else if(modi == 4){
                        fBinsPt[i]      = fBinsEta7TeVPCMEMCPt[i];
                        if (i < fNBinsPt+1){
                          if(!setPi0.CompareTo("Pi0EtaBinning"))
                            fNRebin[i] = fBinsPi0EtaBinning7TeVEMCPtRebin[i];
                          else
                            fNRebin[i] = fBinsEta7TeVEMCPtRebin[i];
                        }
                    } else if(modi == 5){
                        fBinsPt[i]      = fBinsEta7TeVPHOSPt[i];
                        if (i < fNBinsPt+1)
                            fNRebin[i]  = fBinsEta7TeVPHOSPtRebin[i];
                    }else if(modi == 40){
                        if (i < fNBinsPt+1) fBinsPt[i] = fBinsEtaPiPlPiMiPiZero7TevPtPCM[i];
                        if (i < fNBinsPt)   fNRebin[i] = fBinsEtaPiPlPiMiPiZero7TevPtRebinPCM[i];
                    } else if(modi == 41){
                        if (i < fNBinsPt+1) fBinsPt[i] = fBinsEtaPiPlPiMiPiZero7TevPtPCMEMC[i];
                        if (i < fNBinsPt)   fNRebin[i] = fBinsEtaPiPlPiMiPiZero7TevPtRebinPCMEMC[i];
                    } else if(modi == 42){
                        if (i < fNBinsPt+1) fBinsPt[i] = fBinsEtaPiPlPiMiPiZero7TevPtPCMPHOS[i];
                        if (i < fNBinsPt)   fNRebin[i] = fBinsEtaPiPlPiMiPiZero7TevPtRebinPCMPHOS[i];
                    } else if(modi == 44){
                        if (i < fNBinsPt+1) fBinsPt[i] = fBinsEtaPiPlPiMiPiZero7TevPtEMC[i];
                        if (i < fNBinsPt)   fNRebin[i] = fBinsEtaPiPlPiMiPiZero7TevPtRebinEMC[i];
                    } else if(modi == 45){
                        if (i < fNBinsPt+1) fBinsPt[i] = fBinsEtaPiPlPiMiPiZero7TevPtPHOS[i];
                        if (i < fNBinsPt)   fNRebin[i] = fBinsEtaPiPlPiMiPiZero7TevPtRebinPHOS[i];
                    } else {
                      fBinsPt[i]      = fBinsEta7TeVPt[i];
                      if (i < fNBinsPt+1)
                        fNRebin[i]  = fBinsEta7TeVPtRebin[i];
                    }

                }
                nIterBGFit                  = 12;
            //*********************************************************************************************
            //********************************** Eta for pp 8TeV*******************************************
            //*********************************************************************************************
            } else if (energy.CompareTo("8TeV") == 0) {

                fStartPtBin = 1;
                if (modi == 2 ) fStartPtBin = 2;
                if (modi == 4 ) fStartPtBin = 5;

                if (modi == 0 && specialTrigg == 1) fStartPtBin = 10;
                if (modi == 0 && specialTrigg == 2) fStartPtBin = 12;
                if (modi == 2 && specialTrigg == 1) fStartPtBin = 10;
                if (modi == 2 && specialTrigg == 2) fStartPtBin = 14;
                if (modi == 4 && specialTrigg == 1) fStartPtBin = 10;
                if (modi == 4 && specialTrigg == 2) fStartPtBin = 14;

                if ( fNBinsPt > 17 && isDCA ) {
                    cout << "You have chosen to have more than 17 bins for Eta, this is not possible, it will be reduced to 12" << endl;
                    fNBinsPt            = 17;
                } else if (fNBinsPt > 16 && modi != 2 && modi != 3 && modi != 4) {
                    if( specialTrigg == 2 && fNBinsPt > 23){
                        cout << "You have chosen to have more than 23 bins, this is not possible, it will be reduced to 23" << endl;
                        fNBinsPt        = 23;
                    } else if ( specialTrigg == 1 && fNBinsPt > 23){
                        cout << "You have chosen to have more than 23 bins, this is not possible, it will be reduced to 23" << endl;
                        fNBinsPt = 23;
                    } else if(specialTrigg!=1 && specialTrigg!=2 && fNBinsPt >21) {
                        cout << "You have chosen to have more than 21 bins for Eta, this is not possible, it will be reduced to 21" << endl;
                        fNBinsPt        = 21;
                    }
                } else if (fNBinsPt > 19 && (modi == 4)){
                    if( setPi0.CompareTo("Pi0EtaBinning") == 0 && ( specialTrigg ==1 && fNBinsPt > 24) ){
                        fNBinsPt        = 24;
                    } else if ( setPi0.CompareTo("Pi0EtaBinning") == 0 && ( specialTrigg == 2 && fNBinsPt > 19) ){
                        fNBinsPt        = 19;
                    } else if( specialTrigg == 2 && fNBinsPt > 23){
                        cout << "You have chosen to have more than 23 bins, this is not possible, it will be reduced to 23" << endl;
                        fNBinsPt        = 23;
                    } else if ( specialTrigg == 1 && fNBinsPt > 26){
                        cout << "You have chosen to have more than 26 bins, this is not possible, it will be reduced to 26" << endl;
                        fNBinsPt = 26;
                    } else if(specialTrigg!=1 && specialTrigg!=2 && fNBinsPt >21) {
                        cout << "You have chosen to have more than 21 bins for Eta, this is not possible, it will be reduced to 21" << endl;
                        fNBinsPt        = 21;
                    }
                } else if (fNBinsPt > 19 && (modi == 2 || modi == 3)){
                    if( specialTrigg == 2 && fNBinsPt > 23){
                        cout << "You have chosen to have more than 23 bins, this is not possible, it will be reduced to 23" << endl;
                        fNBinsPt        = 23;
                    } else if ( specialTrigg == 1 && fNBinsPt > 23){
                        cout << "You have chosen to have more than 23 bins, this is not possible, it will be reduced to 23" << endl;
                        fNBinsPt = 23;
                    } else if(specialTrigg!=1 && specialTrigg!=2) {
                        cout << "You have chosen to have more than 19 bins for Eta, this is not possible, it will be reduced to 19" << endl;
                        fNBinsPt        = 19;
                    }
                }
                GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                for (Int_t i = 0; i < fNBinsPt+1; i++) {
                    if( modi == 0 ){
                        if(specialTrigg == 1){
                            fBinsPt[i]      = fBinsEta8TeVPCMTrigger1Pt[i];
                        } else if (specialTrigg == 2){
                            fBinsPt[i]      = fBinsEta8TeVPCMTrigger2Pt[i];
                        } else
                            fBinsPt[i]      = fBinsEta8TeVPt[i];
                    } else if ( modi == 2 ){
                        if(specialTrigg == 1){
                            fBinsPt[i]      = fBinsEta8TeVTrigger1Pt[i];
                        } else if (specialTrigg == 2){
                            fBinsPt[i]      = fBinsEta8TeVTrigger2Pt[i];
                        } else
                            fBinsPt[i]      = fBinsEta8TeVPCMEMCPt[i];
                    } else if( modi == 4 ){
                        if(specialTrigg == 1){
                            fBinsPt[i]      = fBinsEta8TeVEMCTrigger1Pt[i];
                        } else if (specialTrigg == 2){
                            fBinsPt[i]      = fBinsEta8TeVTrigger2Pt[i];
                        } else
                            fBinsPt[i]      = fBinsEta8TeVEMCPt[i];
                    } else {
                        fBinsPt[i]      = fBinsEta8TeVPt[i];
                    }

                    if ( modi == 0 ){
                        if(specialTrigg == 1){
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsEta8TeVPCMTrigger1PtRebin[i];
                        } else if(specialTrigg == 2){
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsEta8TeVPCMTrigger2PtRebin[i];
                        } else if(!setPi0.CompareTo("Pi0EtaBinning")){
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0EtaBinning8TeVPtRebin[i];
                        } else {
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsEta8TeVPtRebin[i];
                        }
                    } else if ( modi == 2 ){
                        if(specialTrigg == 1){
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsEta8TeVPCMEMCTrigger1PtRebin[i];
                        } else if(specialTrigg == 2){
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsEta8TeVPCMEMCTrigger2PtRebin[i];
                        } else {
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsEta8TeVPCMEMCPtRebin[i];
                        }
                    } else if ( modi == 4 ) {
                        if(specialTrigg == 1){
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsEta8TeVEMCTrigger1PtRebin[i];
                        } else if(specialTrigg == 2){
                            if (i < fNBinsPt+1){
                                fNRebin[i] = fBinsEta8TeVEMCTrigger2PtRebin[i];
                                if(setPi0.CompareTo("Pi0EtaBinning") == 0 && fBinsPt[i]==18) fNRebin[i] = 16;
                            }
                        } else {
                            if (i < fNBinsPt+1) fNRebin[i] = fBinsEta8TeVEMCPtRebin[i];
                        }
                    } else {
                        if (i < fNBinsPt+1) fNRebin[i] = fBinsEta8TeVPtRebin[i];
                    }
                }
                nIterBGFit                  = 8;
                fMaxYFracBGOverIntHist      = 50;
                optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing5";
                optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing6";
                optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing7";

            //*********************************************************************************************
            //********************************** Eta for pp 13TeV******************************************
            //*********************************************************************************************
            } else if (energy.CompareTo("13TeV") == 0 || energy.CompareTo("13TeVRBins") == 0 ) {
                fStartPtBin                 = GetStartBin("Eta", energy, modi, specialTrigg);
                Int_t maxPtBinTheo          = GetBinning( fBinsPt, maxPtBinAvail, "Eta", energy, modi, specialTrigg, isDCA );
                if (fNBinsPt > maxPtBinTheo) {
                    cout << "**************************************************************************************************************************************" << endl;
                    cout << "********************** ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, **********************************" << endl;
                    cout << "You have chosen "<< fNBinsPt << " bins, which is more than the maximal " << maxPtBinAvail << " bins, this is not possible, it will be reduced to " << maxPtBinAvail << endl;
                    cout << "**************************************************************************************************************************************" << endl;
                    fNBinsPt    = maxPtBinTheo;
                }

                GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                if( setPi0.EqualTo("Eta") ) {
                    switch(modi) {
                        case 0:
                            switch(specialTrigg) {
                                case 0:
                                  if (energy.Contains("RBins")) CopyVectorToArray(fBinsEta13TeVPCMTrigINT7RBinsPtRebin,fNRebin); 
				  else                          CopyVectorToArray(fBinsEta13TeVPCMTrigINT7PtRebin,fNRebin); 
				  break;
                                case 4:
                                case 5: CopyVectorToArray(fBinsEta13TeVPCMTrigINT7PtRebin,fNRebin); break;
                                case 1: CopyVectorToArray(fBinsEta13TeVPCMTrigEMC7PtRebin,fNRebin); break;
                                case 2: CopyVectorToArray(fBinsEta13TeVPCMEMCTrigEG1PtRebin, fNRebin); break;
                                case 3: CopyVectorToArray(fBinsEta13TeVPCMEMCTrigEG2PtRebin, fNRebin); break;
                            }
                            break; // neeeded?
                        case 2:
                            switch(specialTrigg) {
                                case 0:
                                case 4:
                                case 5: CopyVectorToArray(fBinsEta13TeVPCMEMCTrigINT7PtRebin,fNRebin); break;
                                case 2: CopyVectorToArray(fBinsEta13TeVPCMEMCTrigEG1PtRebin, fNRebin); break;
                                case 3: CopyVectorToArray(fBinsEta13TeVPCMEMCTrigEG2PtRebin, fNRebin); break;
                            }
                            break;
                        case 3:
                            switch(specialTrigg) {
                                default:
                                    CopyVectorToArray(fBinsEta13TeVPCMPHOSTrigINT7PtRebin,fNRebin); break;
                            }
                            break;
                        case 4:
                            switch(specialTrigg) {
                                case 0:
                                case 4:
                                case 5: CopyVectorToArray(fBinsEta13TeVEMCTrigINT7PtRebin,fNRebin); break;
                                case 2: CopyVectorToArray(fBinsEta13TeVEMCTrigEG1PtRebin, fNRebin); break;
                                case 3: CopyVectorToArray(fBinsEta13TeVEMCTrigEG2PtRebin, fNRebin); break;
                            }
                            break;
                        case 5:  CopyVectorToArray(fBinsEta13TeVPHOSTrigINT7PtRebin,     fNRebin); break;
                        case 40: CopyVectorToArray(fBinsEtaPiPlPiMiPiZero13TevPtRebinPCM,fNRebin); break;
                        default: CopyVectorToArray(fBinsEta13TeVPCMEMCTrigINT7PtRebin,   fNRebin); break;
                    }
                } else { // Pi0Eta binning
                    switch(modi) {
                        case 0: CopyVectorToArray(fBinsPi0Eta13TeVPtPCMTrigINT7Rebin,fNRebin); break;
                        case 2:
                            switch(specialTrigg) {
                                case 0:
                                case 4:
                                case 5: CopyVectorToArray(fBinsPi0Eta13TeVPCMEMCTrigINT7PtRebin,fNRebin); break;
                                case 2: CopyVectorToArray(fBinsPi0Eta13TeVPCMEMCTrigEG1PtRebin, fNRebin); break;
                                case 3: CopyVectorToArray(fBinsPi0Eta13TeVPCMEMCTrigEG2PtRebin, fNRebin); break;
                            }
                            break;
                        case 3:
                            switch(specialTrigg) {
                                default:
                                    CopyVectorToArray(fBinsPi0EtaBinning13TeVPCMPHOSTrigINT7PtRebin,fNRebin); break;
                            }
                            break;
                        case 4:
                            switch(specialTrigg) {
                                case 0:
                                case 4:
                                case 5: CopyVectorToArray(fBinsPi0Eta13TeVPtEMCTrigINT7Rebin,fNRebin); break;
                                case 2: CopyVectorToArray(fBinsPi0Eta13TeVEMCTrigEG1PtRebin, fNRebin); break;
                                case 3: CopyVectorToArray(fBinsPi0Eta13TeVEMCTrigEG2PtRebin, fNRebin); break;
                            }
                            break;
                        default: CopyVectorToArray(fBinsPi0Eta13TeVPCMEMCTrigINT7PtRebin,fNRebin); break;
                    }
                }
                if ( setPi0.EqualTo("Pi0EtaBinning") ) nIterBGFit = 12;

                nIterBGFit                  = 7;
                fMaxYFracBGOverIntHist      = 60;
                optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing7";
                optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing3";

            //*********************************************************************************************
            // ********************************* Eta for 13TeV low B field ********************************
            //*********************************************************************************************
            } else if (energy.CompareTo("13TeVLowB") == 0) {
                fStartPtBin                 = GetStartBin("Eta", energy, modi, specialTrigg, centrality);
                Int_t maxPtBinTheo          = GetBinning( fBinsPt, maxPtBinAvail, "Eta", energy, modi, specialTrigg, isDCA, centrality );
                if (fNBinsPt > maxPtBinTheo) {
                    cout << "**************************************************************************************************************************************" << endl;
                    cout << "********************** ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, **********************************" << endl;
                    cout << "You have chosen "<< fNBinsPt << " bins, which is more than the maximal " << maxPtBinTheo << " bins, this is not possible, it will be reduced to " << maxPtBinTheo << endl;
                    cout << "**************************************************************************************************************************************" << endl;
                    fNBinsPt    = maxPtBinTheo;
                }
                GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                if( modi==0 || modi == 2  || modi == 4 || modi == 5 || modi == 12 ) CopyVectorToArray(fBinsEta13TeVLowBPtRebin,fNRebin);
                nIterBGFit                  = 8;
                fMaxYFracBGOverIntHist      = 70;
                optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing7";
            //*********************************************************************************************
            //********************************** Eta for pPb 5.023TeV**************************************
            //*********************************************************************************************
            } else if( energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0 || energy.CompareTo("pPb_5.023TeVRun2") == 0) {
                fStartPtBin                 = GetStartBin("Eta", energy, modi, specialTrigg, centrality);
                Int_t maxPtBinTheo          = GetBinning( fBinsPt, maxPtBinAvail, "Eta", energy, modi, specialTrigg, isDCA, centrality );
                if (fNBinsPt > maxPtBinAvail) {
                    cout << "**************************************************************************************************************************************" << endl;
                    cout << "********************** ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, **********************************" << endl;
                    cout << "You have chosen "<< fNBinsPt << " bins, which is more than the maximal " << maxPtBinAvail << " bins, this is not possible, it will be reduced to " << maxPtBinAvail << endl;
                    cout << "**************************************************************************************************************************************" << endl;
                    fNBinsPt    = maxPtBinAvail;
                }
                GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                for (Int_t i = 0; i < fNBinsPt; i++) {
                    if (!setPi0.CompareTo("Eta")){
                        if (modi == 0 && energy.CompareTo("pPb_5.023TeV") == 0 && (centrality.Contains("0-100%") && !centrality.Contains("60-100%")))                                       // MB eta for PCM run 1
                            fNRebin[i]  = fBinsEtapPb5TeVPtRebin[i];
                        else if (modi == 0 && (energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0))
                            fNRebin[i]  = fBinsEtapPb5TeVPCMCentPtRebin[i];
                        else if (modi == 0 && energy.CompareTo("pPb_5.023TeVRun2") == 0 &&
                            ( !centrality.CompareTo("0-10%") || !centrality.CompareTo("0-5%")  || !centrality.CompareTo("5-10%") || !centrality.CompareTo("0-20%") || !centrality.CompareTo("20-40%") || centrality.Contains("40-60") || !centrality.CompareTo("60-100%")) )
                            fNRebin[i]  = fBinsEtapPb5TeVPCMR2CentPtRebin[i];
                        else if (modi == 0 && energy.CompareTo("pPb_5.023TeVRun2") == 0 &&
                            (!centrality.CompareTo("0-1%") || !centrality.CompareTo("0-2%") ) )
                            fNRebin[i]  = fBinsEtapPb5TeVPCMR2Cent2PtRebin[i];
                        else if (modi == 0 && energy.CompareTo("pPb_5.023TeVRun2") == 0)                                                                // MB eta for PCM run 2
                            fNRebin[i]  = fBinsEtapPb5TeVPCMR2PtRebin[i];
                        else if (modi == 1)                                                                                                             // MB eta for Dalitz run 1
                            fNRebin[i]  = fBinsEtapPb5TeVDalitzPtRebin[i];
                        else if (modi == 2 && specialTrigg == 0  && energy.CompareTo("pPb_5.023TeV") == 0 && (centrality.Contains("0-100%") && !centrality.Contains("60-100%")))           // MB eta for PCM-EMC run 1
                            fNRebin[i]  = fBinsEtapPb5TeVPCMEMCPtRebin[i];
                        else if (modi == 2 && specialTrigg == 0  && (energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0) ) // cent dependent eta for PCM-EMC run 1
                            fNRebin[i]  = fBinsEtapPb5TeVPCMEMCCentPtRebin[i];
                        else if ((modi == 2 || modi == 13 ) && specialTrigg == 0  && energy.CompareTo("pPb_5.023TeVRun2") == 0  && !centrality.CompareTo("0-100%") )     // MB eta for PCM-EMC run 2
                            fNRebin[i]  = fBinsEtapPb5TeVPCMEMCR2PtRebin[i];
                        else if (   modi == 2 && specialTrigg == 0  && energy.CompareTo("pPb_5.023TeVRun2") == 0 &&                                     // cent dependent eta for PCM-EMC run 2: 0020, 0010
                                    (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-10%") ) )
                            fNRebin[i]  = fBinsEtapPb5TeVPCMEMCR2CentPtRebin[i];
                        else if (   modi == 2 && specialTrigg == 0  && energy.CompareTo("pPb_5.023TeVRun2") == 0 &&                                     // cent dependent eta for PCM-EMC run 2: 0510, 0005, 60100
                                    (!centrality.CompareTo("5-10%") || !centrality.CompareTo("0-5%") ||  !centrality.CompareTo("60-100%"))  )
                            fNRebin[i]  = fBinsEtapPb5TeVPCMEMCR2SCentPtRebin[i];
                        else if (modi == 2 && specialTrigg == 0  && energy.CompareTo("pPb_5.023TeVRun2") == 0 )                                         // cent dependent eta for PCM-EMC run 2: all but previously discussed
                            fNRebin[i]  = fBinsEtapPb5TeVCentPCMEMCR2PtRebin[i];
                        else if (modi == 2 && specialTrigg == 1)                                                                                        // EMC7 triggered MB eta for PCM-EMC run 1
                            fNRebin[i]  = fBinsEtapPb5TeVPCMEMCTrigEMC7PtRebin[i];
                        else if (modi == 2 && specialTrigg == 2)                                                                                        // EG2 triggered MB eta for PCM-EMC run 1
                            fNRebin[i]  = fBinsEtapPb5TeVPCMEMCTrigEG2PtRebin[i];
                        else if (modi == 2 && specialTrigg == 3)                                                                                        // EG1 triggered MB eta for PCM-EMC run 1
                            fNRebin[i]  = fBinsEtapPb5TeVPCMEMCTrigPtRebin[i];
                        else if (modi == 3 && energy.CompareTo("pPb_5.023TeV") == 0 && (centrality.Contains("0-100%") && !centrality.Contains("60-100%")))                                 // MB eta for PCM-PHOS run 1
                            fNRebin[i]  = fBinsEtapPb5TeVPCMPHOSPtRebin[i];
                        else if (modi == 3 && (energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0) )                 // cent dependent eta for PCM-PHOS run 1
                            fNRebin[i]  = fBinsEtapPb5TeVPCMPHOSCentPtRebin[i];
                        else if (modi == 3 && energy.CompareTo("pPb_5.023TeVRun2") == 0 && !centrality.CompareTo("0-100%"))                             // MB eta for PCM-PHOS run 2
                            fNRebin[i]  = fBinsEtapPb5TeVPCMPHOSR2PtRebin[i];
                        else if (modi == 3 && energy.CompareTo("pPb_5.023TeVRun2") == 0 && !centrality.CompareTo("0-20%"))                         // cent dependent eta for PCM-PHOS run 2: 0020
                            fNRebin[i]  = fBinsEtapPb5TeVPCMPHOSR2CentPtRebin[i];
                        else if (modi == 3 && energy.CompareTo("pPb_5.023TeVRun2") == 0 &&                                                         // cent dependent eta for PCM-PHOS run 2: 0005, 0510
                                (!centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%")) )
                            fNRebin[i]  = fBinsEtapPb5TeVPCMPHOSR2SCentPtRebin[i];
                        else if (modi == 3 && energy.CompareTo("pPb_5.023TeVRun2") == 0 && !centrality.CompareTo("60-100%"))                            // cent dependent eta for PCM-PHOS run 2: 60100
                            fNRebin[i]  = fBinsEtapPb5TeVPCMPHOSR2PerPtRebin[i];
                        else if (modi == 3 && energy.CompareTo("pPb_5.023TeVRun2") == 0)                                                                // cent dependent eta for PCM-PHOS run 2
                            fNRebin[i]  = fBinsEtapPb5TeVPCMPHOSR2PtRebin[i];
                        else if (modi == 4 && energy.CompareTo("pPb_5.023TeV") == 0 && specialTrigg == 1)                                               // EMC7 eta for EMC run 1
                            fNRebin[i]  = fBinsEtapPb5TeVEMCTrigEMC7PtRebin[i];
                        else if (modi == 4 && energy.CompareTo("pPb_5.023TeV") == 0 && specialTrigg == 2)                                               // EG2 eta for EMC run 1
                            fNRebin[i]  = fBinsEtapPb5TeVEMCTrigEG2PtRebin[i];
                        else if (modi == 4 && energy.CompareTo("pPb_5.023TeV") == 0 && specialTrigg == 3)                                               // EG1 eta for EMC run 1
                            fNRebin[i]  = fBinsEtapPb5TeVEMCTrigEG1PtRebin[i];
                        else if (modi == 4 && energy.CompareTo("pPb_5.023TeV") == 0 && (centrality.Contains("0-100%") && !centrality.Contains("60-100%")))                                 // MB eta for EMC run 1
                            fNRebin[i]  = fBinsEtapPb5TeVEMCPtRebin[i];
                        else if (modi == 4 && (energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0))                  // cent dependent eta for EMC run 1
                            fNRebin[i]  = fBinsEtapPb5TeVEMCCentPtRebin[i];
                        else if (modi == 4 && energy.CompareTo("pPb_5.023TeVRun2") == 0   )                           // MB eta for EMC run 2
                            fNRebin[i]  = fBinsEtapPb5TeVEMCR2CentPtRebin[i];
                        else if (modi == 5 && energy.CompareTo("pPb_5.023TeV") == 0 && !centrality.CompareTo("0-100%") )                                // MB eta for PHOS run 1
                            fNRebin[i]  = fBinsEtapPb5TeVPHOSPtRebin[i];
                        else if (modi == 5 && energy.CompareTo("pPb_5.023TeV") == 0)                                                                    // cent dependent eta for PHOS run 1
                            fNRebin[i]  = fBinsEtapPb5TeVPHOSCentPtRebin[i];
                        else if (modi == 5 && energy.CompareTo("pPb_5.023TeVRun2") == 0 && (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-5%") || !centrality.CompareTo("5-10%") || !centrality.CompareTo("60-100%") ))                                                                                          // cent dependent eta for PHOS run 2
                            fNRebin[i]  = fBinsEtapPb5TeVPHOSR2CentPtRebin[i];
                        else if (modi == 5 && energy.CompareTo("pPb_5.023TeVRun2") == 0)                                                                // MB eta for PHOS run 2
                            fNRebin[i]  = fBinsEtapPb5TeVPHOSR2PtRebin[i];

                    } else if (!setPi0.CompareTo("Pi0EtaBinning")){
                        if (modi == 0 && energy.CompareTo("pPb_5.023TeV") == 0 && (centrality.Contains("0-100%") && !centrality.Contains("60-100%")))                                      // MB pi0-eta for PCM run 1
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPtRebin[i];
                        else if (modi == 0 && (energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0  ))                                                                   // MB pi0-eta for PCM run 1
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPCMCentPtRebin[i];
                        else if (modi == 0 && energy.CompareTo("pPb_5.023TeVRun2") == 0 &&                                                              // cent dependent pi0-eta for PCM run 2: 0020, 0010, 0005, 0510, 60100
                            (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-10%") || !centrality.CompareTo("0-5%")  || !centrality.CompareTo("5-10%") ||  !centrality.CompareTo("60-100%") ) )
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPCMR2CentPtRebin[i];
                        else if (modi == 0 && energy.CompareTo("pPb_5.023TeVRun2") == 0)                                                                // MB pi0-eta for PCM run 2 + cent dependent
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPCMR2PtRebin[i];
                        else if (modi == 1)                                                                                                             // MB pi0-eta for PCM-Dalitz run 1
                            fNRebin[i]  = fBinsPi0EtapPb5TeVDalitzPtRebin[i];
                        else if (modi == 2 && specialTrigg == 0  && energy.CompareTo("pPb_5.023TeV") == 0 && (centrality.Contains("0-100%") && !centrality.Contains("60-100%")))           // MB pi0-eta for PCM-EMC run 1
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPCMEMCPtRebin[i];
                        else if (modi == 2 && specialTrigg == 0  && (energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0))                                             // cent dependent pi0-eta for PCM-EMC run 1
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPCMEMCCentPtRebin[i];
                        else if ((modi == 2 || modi == 13 )  && specialTrigg == 0  && energy.CompareTo("  pPb_5.023TeVRun2") == 0  && !centrality.CompareTo("0-100%") )   // MB pi0-eta for PCM-EMC run 2
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPCMEMCR2PtRebin[i];
                        else if (   modi == 2 && specialTrigg == 0  && energy.CompareTo("pPb_5.023TeVRun2") == 0  &&                                    // cent dependent pi0-eta for PCM-EMC run 2: 0020, 0010
                                    (!centrality.CompareTo("0-20%") || !centrality.CompareTo("0-10%") ))
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPCMEMCR2CentPtRebin[i];
                        else if (   modi == 2 && specialTrigg == 0  && energy.CompareTo("pPb_5.023TeVRun2") == 0  &&                                    // cent dependent pi0-eta for PCM-EMC run 2: 0000, 0510, 60100
                                    (!centrality.CompareTo("5-10%") || !centrality.CompareTo("0-5%") ||  !centrality.CompareTo("60-100%")) )
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPCMEMCR2SCentPtRebin[i];
                        else if (modi == 2 && specialTrigg == 1)                                                                                        // EMC7 triggered MB pi0-eta for PCM-EMC run 1
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPCMEMCTrigEMC7PtRebin[i];
                        else if (modi == 2 && specialTrigg == 2)                                                                                        // EG2 triggered MB pi0-eta for PCM-EMC run 1
                            fNRebin[i]  = fBinsEtapPb5TeVPCMEMCTrigEG2PtRebin[i];
                        else if (modi == 2 && specialTrigg == 3)                                                                                        // EG1 triggered MB pi0-eta for PCM-EMC run 1
                            fNRebin[i]  = fBinsEtapPb5TeVPCMEMCTrigPtRebin[i];
                        else if (modi == 3  && energy.CompareTo("pPb_5.023TeV") == 0 && (centrality.Contains("0-100%") && !centrality.Contains("60-100%")) )                               // MB pi0-eta for PCM-PHOS run 1
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPCMPHOSPtRebin[i];
                        else if (modi == 3  && (energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0))                                                                 // cent dependent pi0-eta for PCM-PHOS run 1
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPCMPHOSCentPtRebin[i];
                        else if (modi == 3  && energy.CompareTo("pPb_5.023TeVRun2") == 0 )                                                              // MB pi0-eta for PCM-PHOS run 2
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPCMPHOSR2PtRebin[i];
                        else if (modi == 4 && energy.CompareTo("pPb_5.023TeV") == 0 && (centrality.Contains("0-100%") && !centrality.Contains("60-100%")))                                 // MB pi0-eta for EMC run 1
                            fNRebin[i]  = fBinsPi0EtapPb5TeVEMCPtRebin[i];
                        else if (modi == 4 && (energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0))                                                                    // cent dependent for EMC run 1
                            fNRebin[i]  = fBinsPi0EtapPb5TeVEMCCentPtRebin[i];
                        else if (modi == 4 && energy.CompareTo("pPb_5.023TeVRun2") == 0)                                                                // MB pi0-eta for EMC run 2
                            fNRebin[i]  = fBinsPi0EtapPb5TeVEMCR2PtRebin[i];
                        else if (modi == 5 && energy.CompareTo("pPb_5.023TeV") == 0 && !centrality.CompareTo("0-100%"))                                 // MB pi0-eta for PHOS run 1
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPHOSPtRebin[i];
                        else if (modi == 5 && energy.CompareTo("pPb_5.023TeV") == 0)                                                                    // cent dependent PHOS run 1
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPHOSCentPtRebin[i];
                        else if (modi == 5 && energy.CompareTo("pPb_5.023TeVRun2") == 0)                                                                // MB pi0-eta for PHOS run 2
                            fNRebin[i]  = fBinsPi0EtapPb5TeVPHOSR2PtRebin[i];
                    }
                }

                if (!setPi0.CompareTo("Pi0EtaBinning")){
                    optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                    optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                    optionBGSmoothingVar2       = "noSmoothing";
                    nIterBGFit                  = 11;
                    if (modi == 0){
                      if( !energy.CompareTo("pPb_5.023TeVCent") )
                        nIterBGFit                  = 8;
                      else if( !energy.CompareTo("pPb_5.023TeVRun2") )
                        nIterBGFit                  = 8;
                      else if( !energy.CompareTo("pPb_5.023TeV") ){
                        if(centrality.Contains("0-20%"))
                          nIterBGFit                = 8;
                        else if(centrality.Contains("20-40%"))
                          nIterBGFit                = 8;
                        else if(centrality.Contains("40-60%"))
                          nIterBGFit                = 8;
                        else if(centrality.Contains("60-100%"))
                          nIterBGFit                = 8;
                      }
                    }
                } else {
                    optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                    optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                    optionBGSmoothingVar2       = "noSmoothing";
                    nIterBGFit                  = 11;
                    if (modi == 0){
                      if( !energy.CompareTo("pPb_5.023TeVRun2") )
                        nIterBGFit                  = 8;
                        if(centrality.Contains("0-100%"))
                            nIterBGFit                  = 7;
                      else if( !energy.CompareTo("pPb_5.023TeV") ){
                        if(centrality.CompareTo(""))
                          nIterBGFit                = 8;
                      }
                    }
                }
                fMaxYFracBGOverIntHist          = 20;
            //*********************************************************************************************
            //********************************** Eta for pPb 8TeV**************************************
            //*********************************************************************************************
            } else if( energy.CompareTo("pPb_8TeV") == 0 || energy.CompareTo("pPb_8TeVRun2") == 0) {
                fStartPtBin         = GetStartBin("Eta", energy, modi, specialTrigg);
                if (fNBinsPt > 16 && isDCA) {
                    cout << "You have chosen to have more than 16 DCA bins, this is not possible, it will be reduced to 16" << endl;
                    fNBinsPt        = 16;
                } else if (fNBinsPt > 20 && modi < 2) {
                    cout << "You have chosen to have more than 20 bins, this is not possible, it will be reduced to 20" << endl;
                    fNBinsPt        = 20;
                } else if (fNBinsPt > 21 && ( modi == 2 || modi == 4) && specialTrigg == 0 ) {
                    cout << "You have chosen to have more than 21 bins, this is not possible, it will be reduced to 21" << endl;
                    fNBinsPt        = 21;
                } else if (fNBinsPt > 22 && specialTrigg == 0){
                    cout << "You have chosen to have more than 22 bins, this is not possible, it will be reduced to 22" << endl;
                    fNBinsPt        = 22;
                } else if (fNBinsPt > 26 ){
                    cout << "You have chosen to have more than 26 bins, this is not possible, it will be reduced to 26" << endl;
                    fNBinsPt        = 26;
                }
                GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                for (Int_t i = 0; i < fNBinsPt+1; i++) {
                    if (modi == 2 && specialTrigg == 1){
                        fBinsPt[i]      = fBinsEtapPb8TeVPtEMCTrig[i];
                        if (i < fNBinsPt+1){
                            fNRebin[i]  = fBinsEtapPb8TeVPCMEMCTrigEMC7PtRebin[i];
                        }
                    } else if (modi == 2 && specialTrigg == 2){
                        fBinsPt[i]      = fBinsEtapPb8TeVPtEMCTrig[i];
                        if (i < fNBinsPt+1){
                            fNRebin[i]  = fBinsEtapPb8TeVPCMEMCTrigEG2PtRebin[i];
                        }
                    } else if (modi == 2 && specialTrigg == 3){
                        fBinsPt[i]      = fBinsEtapPb8TeVPtEMCTrig[i];
                        if (i < fNBinsPt+1){
                            fNRebin[i]  = fBinsEtapPb8TeVPCMEMCTrigPtRebin[i];
                        }
                    } else {
                        // PCM binning
                        if ( modi == 0){
                          if (isDCA )
                              fBinsPt[i]      = fBinsEtapPb8TeVPtDCA[i];
                          else
                              fBinsPt[i]      = fBinsEtapPb8TeVPt[i];
                        // Dalitz binning
                        } else if (modi == 1){
                            fBinsPt[i]      = fBinsEtapPb8TeVDalitzPt[i];
                        // EMC and PCM-EMC binning
                        } else if (modi == 2 || modi == 4){
                            fBinsPt[i]      = fBinsEtapPb8TeVEMCPt[i];
                        } else if (modi == 3){
                            fBinsPt[i]      = fBinsEtapPb8TeVPCMPHOSPt[i];
                        } else if (modi == 5){
                            fBinsPt[i]      = fBinsEtapPb8TeVPHOSPt[i];
                        } else {
                            fBinsPt[i]      = fBinsEtapPb8TeVPt[i];
                        }
                        // Rebin factors
                        if (i < fNBinsPt+1){
                            if (modi == 0 && !setPi0.CompareTo("Eta"))
                                fNRebin[i]  = fBinsEtapPb8TeVPtRebin[i];
                            else if (modi == 0 && !setPi0.CompareTo("Pi0EtaBinning"))
                                fNRebin[i]  = fBinsPi0EtapPb8TeVPtRebin[i];
                            else if (modi == 1)
                                fNRebin[i]  = fBinsEtapPb8TeVDalitzPtRebin[i];
                            else if (modi == 2 && !setPi0.CompareTo("Eta"))
                                fNRebin[i]  = fBinsEtapPb8TeVPCMEMCPtRebin[i];
                            else if (modi == 2 && !setPi0.CompareTo("Pi0EtaBinning"))
                                fNRebin[i]  = fBinsPi0EtapPb8TeVPCMEMCPtRebin[i];
                            else if (modi == 3 && !setPi0.CompareTo("Eta"))
                                fNRebin[i]  = fBinsEtapPb8TeVPCMPHOSPtRebin[i];
                            else if (modi == 3 && !setPi0.CompareTo("Pi0EtaBinning"))
                                fNRebin[i]  = fBinsPi0EtapPb8TeVPCMPHOSPtRebin[i];
                            else if (modi == 4 && !setPi0.CompareTo("Eta"))
                                fNRebin[i]  = fBinsEtapPb8TeVEMCPtRebin[i];
                            else if (modi == 4 && !setPi0.CompareTo("Pi0EtaBinning"))
                                fNRebin[i]  = fBinsPi0EtapPb8TeVEMCPtRebin[i];
                            else if (modi == 5 && !setPi0.CompareTo("Eta"))
                                fNRebin[i]  = fBinsEtapPb8TeVPHOSPtRebin[i];
                            else if (modi == 5 && !setPi0.CompareTo("Pi0EtaBinning"))
                                fNRebin[i]  = fBinsPi0EtapPb8TeVPHOSPtRebin[i];
                            else
                                fNRebin[i]  = fBinsEtapPb8TeVPtRebin[i];
                        }
                    }
                }

                if (!setPi0.CompareTo("Pi0EtaBinning")){
                    optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                    optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                    optionBGSmoothingVar2       = "noSmoothing";
                    nIterBGFit                  = 11;
                } else {
                    optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing3";
                    optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing5";
                    optionBGSmoothingVar2       = "noSmoothing";
                    nIterBGFit                  = 11;
                }
                fMaxYFracBGOverIntHist          = 20;
            //*********************************************************************************************
            //********************************** Eta for PbPb 2.76TeV**************************************
            //*********************************************************************************************
            } else if( energy.CompareTo("PbPb_2.76TeV") == 0) {
                fStartPtBin         = 2;
                if (modi == 4){
                    fStartPtBin     = 5;
                } else if (modi == 2) {
                    fStartPtBin     = 3;
                }
                if (isDCA){
                    if (!setPi0.CompareTo("Pi0EtaBinning"))
                        fStartPtBin     = 1; //otherwise usually 3
                    else
                        fStartPtBin     = 4; //otherwise usually 3
                }
                if (isDCA) {
                    if (!setPi0.CompareTo("Pi0EtaBinning")) {
                        if (fNBinsPt > 10) {
                            cout << "You have chosen to have more than 10 bins, this is not possible, it will be reduced to 10" << endl;
                            fNBinsPt            = 10;
                        }
                    } else {
                        if (fNBinsPt > 16) {
                            cout << "You have chosen to have more than 16 bins, this is not possible, it will be reduced to 16" << endl;
                            fNBinsPt                    = 16;
                        }
                    }
                } else if (modi != 4 && modi != 2 &&    fNBinsPt > 12) {
                    cout << "You have chosen to have more than 12 bins, this is not possible, it will be reduced to 12" << endl;
                    fNBinsPt        = 12;
                }
                if ((modi == 4 || modi == 2) &&   fNBinsPt > 14) {
                    cout << "You have chosen to have more than 14 bins, this is not possible, it will be reduced to 14" << endl;
                    fNBinsPt        = 14;
                }
                GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                for (Int_t i = 0; i < fNBinsPt+1; i++) {
                    if (isDCA){
                        if (!setPi0.CompareTo("Pi0EtaBinning")) {
                            fBinsPt[i]          = fBinsEtaPbPb2760GeVPt[i];
                        } else {
                            fBinsPt[i]          = fBinsEtaPbPb2760GeVPtDCA[i];
                        }
                    } else {
                        if (modi == 0)
                            fBinsPt[i]          = fBinsEtaPbPb2760GeVPtLHC11hLessBins[i]; //fBinsEtaPbPb2760GeVPtLHC11h[i];
                        else if (modi == 2 || modi == 4)
                            fBinsPt[i]          = fBinsEtaPbPb2760GeVPtLHC11hEMCBins[i];
                        else
                            fBinsPt[i]          = fBinsEtaPbPb2760GeVPtLHC11hLessBins[i];
                        if (i < fNBinsPt+1) fNRebin[i] = fBinsEtaPbPb2760GeVPtRebinLHC11hLessBins[i]; // fBinsEtaPbPb2760GeVPtRebinLHC11h[i]; //fBinsEtaPbPb2760GeVPtRebinLHC11hFinerBinning[i];
                    }
                }
                optionBGSmoothingStandard       = "BackDecreasingWindow,BackSmoothing5";
                optionBGSmoothingVar1           = "BackDecreasingWindow,BackSmoothing3";
                optionBGSmoothingVar2           = "BackDecreasingWindow,BackSmoothing7";
                fMaxYFracBGOverIntHist          = 8;
                if (!setPi0.CompareTo("Pi0EtaBinning")) {
                    optionBGSmoothingStandard       = "BackSmoothing9";
                    optionBGSmoothingVar1           = "BackSmoothing7";
                    optionBGSmoothingVar2           = "BackSmoothing11";
                    if (!centDCA.CompareTo("60-80%")){
                        nIterBGFit              = 15;
                    } else if (!centDCA.CompareTo("40-60%")){
                        nIterBGFit              = 17;
                    } else if (!centDCA.CompareTo("20-50%")){
                        nIterBGFit              = 19;
                    } else if (!centDCA.CompareTo("20-40%")){
                        nIterBGFit              = 19;
                    } else {
                        nIterBGFit              = 21;
                    }
                } else {
                    if (!centDCA.CompareTo("60-80%") || !centDCA.CompareTo("60-70%") || !centDCA.CompareTo("70-80%") || !centDCA.CompareTo("75-90%") ){
                        nIterBGFit                  = 15;
                        fMaxYFracBGOverIntHist      = 15;
                        optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing9";
                        optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing7";
                        optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing11";
                    } else if (!centDCA.CompareTo("50-60%")){
                        nIterBGFit                  = 14;
                        fMaxYFracBGOverIntHist      = 12;
                    } else if (!centDCA.CompareTo("40-60%") || !centDCA.CompareTo("40-50%") ){
                        nIterBGFit                  = 16;
                        fMaxYFracBGOverIntHist      = 10;
                    } else if (!centDCA.CompareTo("30-40%")){
                        nIterBGFit                  = 17;
                    } else if (!centDCA.CompareTo("20-50%")){
                        nIterBGFit                  = 16;
                    } else if (!centDCA.CompareTo("20-40%")){
                        nIterBGFit                  = 16;
                    } else if (!centDCA.CompareTo("20-30%")){
                        nIterBGFit                  = 18;
                    } else {
                        fMaxYFracBGOverIntHist      = 4;
                        nIterBGFit                  = 21;
                    }
                }

            //*********************************************************************************************
            //********************************** Eta for PbPb 5.02TeV**************************************
            //*********************************************************************************************
            } else if( energy.CompareTo("PbPb_5.02TeV") == 0) {
                fStartPtBin         = 1;
                if( modi == 2  || modi == 3){
                  fStartPtBin = 5;
                } else if( modi == 4 ){
                  fStartPtBin = 4;
                }

                if (fNBinsPt > 22) {
                    cout << "You have chosen to have more than 22 bins, this is not possible, it will be reduced to 22" << endl;
                    fNBinsPt        = 22;
                }
                GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
                for (Int_t i = 0; i < fNBinsPt+1; i++) {
                    if(modi == 2 || modi == 3){
                        fBinsPt[i]         = fBinsEtaPbPb5TeVPCMEMCPt[i];
                        if (i < fNBinsPt+1) fNRebin[i] = fBinsEtaPbPb5TeVPCMEMCPtRebin[i];
                    }else if(modi == 4){
                        fBinsPt[i]         = fBinsEtaPbPb5TeVEMCPt[i];
                        if (i < fNBinsPt+1) fNRebin[i] = fBinsEtaPbPb5TeVEMCPtRebin[i];
                    }else{
                        fBinsPt[i]         = fBinsEtaPbPb5TeVPt[i];
                        if (i < fNBinsPt+1) {
			  if (setPi0.CompareTo("Pi0EtaBinning") == 0) fNRebin[i] = fBinsPi0EtaBinningPbPb5TeVPtRebin[i];
			  else                                        fNRebin[i] = fBinsEtaPbPb5TeVPtRebin[i];
			}
                    }
                }

            //*********************************************************************************************
            //********************************** Eta for XeXe 5.44TeV**************************************
            //*********************************************************************************************
            } else if( energy.CompareTo("XeXe_5.44TeV") == 0) {
                fStartPtBin     = GetStartBin("Eta", energy, modi);
                Int_t maxPtBinTheo          = GetBinning( fBinsPt, maxPtBinAvail, "Eta", energy, modi, specialTrigg, isDCA );
                if (fNBinsPt > maxPtBinAvail) {
                    cout << "**************************************************************************************************************************************" << endl;
                    cout << "********************** ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, **********************************" << endl;
                    cout << "You have chosen "<< fNBinsPt << " bins, which is more than the maximal " << maxPtBinAvail << " bins, this is not possible, it will be reduced to " << maxPtBinAvail << endl;
                    cout << "**************************************************************************************************************************************" << endl;
                    fNBinsPt    = maxPtBinAvail;
                }

                GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                for (Int_t i = 0; i < fNBinsPt+1; i++) {
                    if (setPi0.CompareTo("Pi0EtaBinning")){
                        if (i < fNBinsPt+1) fNRebin[i] = fBinsEtaXeXe5440GeVPtRebin[i];
                    } else {
                        if (i < fNBinsPt+1) fNRebin[i] = fBinsPi0EtaBinningXeXe5440GeVPtRebin[i];
                    }
                }
                optionBGSmoothingStandard       = "BackDecreasingWindow,BackSmoothing5";
                optionBGSmoothingVar1           = "BackDecreasingWindow,BackSmoothing3";
                optionBGSmoothingVar2           = "BackDecreasingWindow,BackSmoothing7";
                fMaxYFracBGOverIntHist          = 8;
                if (setPi0.CompareTo("Pi0EtaBinning")) {
                    nIterBGFit                  = 15;
                    fMaxYFracBGOverIntHist      = 15;
                    optionBGSmoothingStandard   = "BackDecreasingWindow,BackSmoothing9";
                    optionBGSmoothingVar1       = "BackDecreasingWindow,BackSmoothing7";
                    optionBGSmoothingVar2       = "BackDecreasingWindow,BackSmoothing11";

                }
            }
        //*************************************************************************************************
        //********************************** Binning for Eta' *********************************************
        //*************************************************************************************************
        } else if (setPi0.CompareTo("EtaPrime") == 0){
            // Initialise members
            fNBinsPt           = numberOfBins;
            fBinsPt            = new Double_t[30];
            fNRebin            = new Int_t[29];
            fStartPtBin        = GetStartBin( "EtaPrime", energy, modi, specialTrigg );
            Int_t maxPtBinTheo = GetBinning( fBinsPt, maxPtBinAvail, "EtaPrime", energy, modi, specialTrigg);
            if (fNBinsPt > maxPtBinTheo) {
                cout << "**************************************************************************************************************************************" << endl;
                cout << "********************** ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, **********************************" << endl;
                cout << "You have chosen "<< fNBinsPt << " bins, which is more than the maximal " << maxPtBinAvail << " bins, this is not possible, it will be reduced to " << maxPtBinAvail << endl;
                cout << "**************************************************************************************************************************************" << endl;
                fNBinsPt    = maxPtBinTheo;
            }
            // Set binning according to energy/mode/trigger cases
            if     (energy.EqualTo("7TeV")) {
                CopyVectorToArray(fBinsEtaPrime7TeVPtRebin,fNRebin);
            } else if(energy.EqualTo("13TeV")) {
                switch(modi) {
                    case 0: // PCM-PCM
                        switch( specialTrigg ) {
                            case 0: CopyVectorToArray(fBinsEtaPrime13TeV_PCM_INT7_PtRebin,fNRebin); break; // 10
                            case 1: CopyVectorToArray(fBinsEtaPrime13TeV_PCM_L0_PtRebin,  fNRebin); break; // 52
                            case 2: CopyVectorToArray(fBinsEtaPrime13TeV_PCM_EG1_PtRebin, fNRebin); break; // 83
                            case 3: CopyVectorToArray(fBinsEtaPrime13TeV_PCM_EG2_PtRebin, fNRebin); break; // 85
                        } break;
                    case 2: // PCM-EMC
                        switch( specialTrigg ) {
                            case 0: CopyVectorToArray(fBinsEtaPrime13TeV_PCMEMC_INT7_PtRebin,fNRebin); break; // 10
                            case 1: CopyVectorToArray(fBinsEtaPrime13TeV_PCMEMC_L0_PtRebin,  fNRebin); break; // 52
                            case 2: CopyVectorToArray(fBinsEtaPrime13TeV_PCMEMC_EG1_PtRebin, fNRebin); break; // 83
                            case 3: CopyVectorToArray(fBinsEtaPrime13TeV_PCMEMC_EG2_PtRebin, fNRebin); break; // 85
                        } break;
                    case 3: // PCM-PHOS
                        switch( specialTrigg ) {
                            case 0: CopyVectorToArray(fBinsEtaPrime13TeV_PCMPHOS_MinBias_PtRebin,fNRebin); break; // 10
                            case 6: CopyVectorToArray(fBinsEtaPrime13TeV_PCMPHOS_VZERO_PtRebin,  fNRebin); break; // 62
                        } break;
                    case 4: // EMC-EMC
                        switch( specialTrigg ) {
                            case 0: CopyVectorToArray(fBinsEtaPrime13TeV_EMC_INT7_PtRebin,fNRebin); break; // 10
                            case 1: CopyVectorToArray(fBinsEtaPrime13TeV_EMC_L0_PtRebin,  fNRebin); break; // 52
                            case 2: CopyVectorToArray(fBinsEtaPrime13TeV_EMC_EG1_PtRebin, fNRebin); break; // 83
                            case 3: CopyVectorToArray(fBinsEtaPrime13TeV_EMC_EG2_PtRebin, fNRebin); break; // 85
                        } break;
                    case 5: // PHOS-PHOS
                        switch( specialTrigg ) {
                            case 0: CopyVectorToArray(fBinsEtaPrime13TeV_PHOS_MinBias_PtRebin,fNRebin); break; // 10
                            case 6: CopyVectorToArray(fBinsEtaPrime13TeV_PHOS_VZERO_PtRebin,  fNRebin); break; // 62
                        } break;
                    case 60: CopyVectorToArray( fBinsEtaPrime13TeV_PCM_INT7_PtRebin,        fNRebin ); break;
                    case 61: CopyVectorToArray( fBinsEtaPrime13TeV_EMC_INT7_PtRebin,        fNRebin ); break;
                    case 63: CopyVectorToArray( fBinsEtaPrime13TeV_PCMPHOS_MinBias_PtRebin, fNRebin ); break;
                    case 64: CopyVectorToArray( fBinsEtaPrime13TeV_EMC_INT7_PtRebin,        fNRebin ); break;
                    case 65: CopyVectorToArray( fBinsEtaPrime13TeV_PHOS_MinBias_PtRebin,    fNRebin ); break;
                }
            }
            // Set optimum columns and rows
            GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);
        //*************************************************************************************************
        //********************************** Binning for Omega ********************************************
        //*************************************************************************************************
        } else if (setPi0.CompareTo("Omega") == 0) {
            fNBinsPt        = numberOfBins;
            fBinsPt         = new Double_t[20];
            fNRebin         = new Int_t[19];
            fStartPtBin     = 0;

            if (energy.CompareTo("7TeV") == 0) {
                fStartPtBin                 = GetStartBin("Omega","7TeV",modi);
                Int_t maxPtBinTheo          = GetBinning( fBinsPt, maxPtBinAvail, "Omega", energy, modi);
                if (fNBinsPt > maxPtBinAvail) {
                    cout << "**************************************************************************************************************************************" << endl;
                    cout << "********************** ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, **********************************" << endl;
                    cout << "You have chosen "<< fNBinsPt << " bins, which is more than the maximal " << maxPtBinAvail << " bins, this is not possible, it will be reduced to " << maxPtBinAvail << endl;
                    cout << "**************************************************************************************************************************************" << endl;
                    fNBinsPt    = maxPtBinAvail;
                }
                GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                for (Int_t i = 0; i < fNBinsPt; i++) {
                    if(modi == 40)      fNRebin[i] = fBinsOmegaPiPlPiMiPiZero7TevPtRebinPCM[i];
                    else if(modi == 41) fNRebin[i] = fBinsOmegaPiPlPiMiPiZero7TevPtRebinPCMEMC[i];
                    else if(modi == 42) fNRebin[i] = fBinsOmegaPiPlPiMiPiZero7TevPtRebinPCMPHOS[i];
                    else if(modi == 44) fNRebin[i] = fBinsOmegaPiPlPiMiPiZero7TevPtRebinEMC[i];
                    else if(modi == 45) fNRebin[i] = fBinsOmegaPiPlPiMiPiZero7TevPtRebinPHOS[i];
                    else                fNRebin[i] = fBinsOmegaPiPlPiMiPiZero7TevPtRebinPCM[i];
                }
            } else if (energy.CompareTo("13TeV") == 0) {
                if (fNBinsPt > 20) {
                    cout << "You have chosen to have more than 15 bins for Omega, this is not possible, it will be reduced to 12" << endl;
                    fNBinsPt = 19;
                }
                fStartPtBin     = GetStartBin("Omega","13TeV",modi);
                GetOptimumNColumnsAndRows(fNBinsPt, fStartPtBin, fColumn, fRow);

                if(modi == 40){
                    for (Int_t i = 0; i < fNBinsPt+2; i++) {
                        fBinsPt[i] = fBinsOmegaPiPlPiMiPiZero13TevPtPCM[i];
                        if (i < fNBinsPt+1)
                            fNRebin[i] = fBinsOmegaPiPlPiMiPiZero13TevPtRebinPCM[i];
                    }
                } else{
                    for (Int_t i = 0; i < fNBinsPt+2; i++) {
                        fBinsPt[i] = fBinsOmegaPiPlPiMiPiZero13TevPtPCM[i];
                        if (i < fNBinsPt+1)
                            fNRebin[i] = fBinsOmegaPiPlPiMiPiZero13TevPtRebinPCM[i];
                    }
                }
            }
        }
        if (fExampleBin < fStartPtBin || fExampleBin > fNBinsPt){
            if (fExampleBin < fStartPtBin)
                fExampleBin = fStartPtBin;
            if (fExampleBin > fNBinsPt)
                fExampleBin = fNBinsPt-1;
            cout << "**************************************************************************************************************************************" << endl;
            cout << "********************** ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION, **********************************" << endl;
            cout << "You have chosen an incompatible Example bin it should lie between " << fStartPtBin << "\t" << fNBinsPt << endl;
            cout << "The example bin has been reset to: " << fExampleBin << ", please fix this in the code"<< endl;
            cout << "**************************************************************************************************************************************" << endl;
            fExampleBin = fStartPtBin;
        }
    }


#endif
