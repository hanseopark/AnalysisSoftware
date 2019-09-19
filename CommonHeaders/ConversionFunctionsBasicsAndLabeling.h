//**************************************************************************************************************************
//******************************** Conversion functions header for basics and plot labeling ********************************
//**************************************************************************************************************************
//*********** this header contains all basic functions for the analysis as well as the correct labeling of plots ***********
//**************************************************************************************************************************
#ifndef GAMMACONV_ConversionFunctionsBasics
#define GAMMACONV_ConversionFunctionsBasics

    #include <iostream>
    #include "TString.h"
    #include "TObjString.h"

    // ****************************************************************************************************************
    //********************* Enums *************************************************************************************
    // ****************************************************************************************************************
    enum Mesons {
        kPi0,
        kEta,
        kEtaPrime,
        kOmega,
        kCKaon,
        kCPion,
        kLambda,
        kRho0,
        kPhi,
        kK0Star,
        kProton
    };
    enum Modes {
        // Standard gamma-gamma analysis modes
        kModePCMPCM         = 0,
        kModePCMDalitz      = 1,
        kModePCMEMC         = 2,
        kModePCMPHOS        = 3,
        kModeEMCEMC         = 4,
        kModePHOSPHOS       = 5,
        kModeEMCDalitz      = 6,
        kModePHOSDalitz     = 7,
        kModePCMPCM_old     = 9,
        kMode_mEMC          = 10,
        kMode_mPHOS         = 11,
        kModeDCalDCal       = 12,
        kModePCMDCal        = 13,
        kModePCMEDC         = 14,
        kModeEDC            = 15,
        kModeComb           = 20,
        // Omega analysis
        kModePCMPCM_om      = 40,
        kModePCMEMC_om      = 41,
        kModePCMPHOS_om     = 42,
        kModePCMDCal_om     = 43,
        kModeEMCEMC_om      = 44,
        kModePHOSPHOS_om    = 45,
        kModeDCalDCal_om    = 46,
        kModePCMDalitz_om   = 47,
        kModeEMCDalitz_om   = 48,
        kModePHOSDalitz_om  = 49,
        kModeDCalDalitz_om  = 50,
        // Heavy meson analysis
        kModePCMPCM_Heavy   = 100,
        kModePCMEMC_Heavy   = 102,
        kModePCMPHOS_Heavy  = 103,
        kModeEMCEMC_Heavy   = 104,
        kModePHOSPHOS_Heavy = 105,
        kModePCMEDC_Heavy   = 114,
        kModeEDC_Heavy      = 115
    };

    enum SpecialTriggers {
        kSpecTrigINT7    = 10,
        kSpecTrigL0      = 52,
        kSpecTrigVHM     = 74,
        kSpecTrigVHMSPD2 = 76,
        kSpecTrigEG1     = 83,
        kSpecTrigEG2     = 85,
        kSpecTrigPHI7    = 62
    };

    // ****************************************************************************************************************
    //********************* Enum functions ****************************************************************************
    // ****************************************************************************************************************
    // ! backward compatibility fix
    Modes ReturnModeEnum( Int_t mode ) {
        switch( mode ) {
            // Standard gamma-gamma analysis modes
            case 0   : return kModePCMPCM;
            case 1   : return kModePCMDalitz;
            case 2   : return kModePCMEMC;
            case 3   : return kModePCMPHOS;
            case 4   : return kModeEMCEMC;
            case 5   : return kModePHOSPHOS;
            case 6   : return kModeEMCDalitz;
            case 7   : return kModePHOSDalitz;
            case 9   : return kModePCMPCM_old;
            case 10  : return kMode_mEMC;
            case 11  : return kMode_mPHOS;
            case 12  : return kModeDCalDCal;
            case 13  : return kModePCMDCal;
            case 14  : return kModePCMEDC;
            case 15  : return kModeEDC;
            case 20  : return kModeComb;
            // Omega analysis
            case 40  : return kModePCMPCM_om;
            case 41  : return kModePCMEMC_om;
            case 42  : return kModePCMPHOS_om;
            case 43  : return kModePCMDCal_om;
            case 44  : return kModeEMCEMC_om;
            case 45  : return kModePHOSPHOS_om;
            case 46  : return kModeDCalDCal_om;
            case 47  : return kModePCMDalitz_om;
            case 48  : return kModeEMCDalitz_om;
            case 49  : return kModePHOSDalitz_om;
            case 50  : return kModeDCalDalitz_om;
            // Heavy meson analysis
            case 100 : return kModePCMPCM_Heavy;
            case 102 : return kModePCMEMC_Heavy;
            case 103 : return kModePCMPHOS_Heavy;
            case 104 : return kModeEMCEMC_Heavy;
            case 105 : return kModePHOSPHOS_Heavy;
            case 114 : return kModePCMEDC_Heavy;
            case 115 : return kModeEDC_Heavy;
            default :
                std::cout << "ERROR in ReturnModeEnum: mode " << mode << " not defined as enum" << std::endl;
                return kModePCMPCM;
        }
    }
    Mesons ReturnMesonEnum( TString meson ) {
        if (meson.EqualTo("Pi0"))      return kPi0;
        if (meson.EqualTo("Eta"))      return kEta;
        if (meson.EqualTo("EtaPrim"))  return kEtaPrime;
        if (meson.EqualTo("EtaPrime")) return kEtaPrime;
        if (meson.EqualTo("Omega"))    return kOmega;
        if (meson.EqualTo("CKaon"))    return kCKaon;
        if (meson.EqualTo("CPion"))    return kCPion;
        if (meson.EqualTo("Lambda"))   return kLambda;
        if (meson.EqualTo("Rho0"))     return kRho0;
        if (meson.EqualTo("Phi"))      return kPhi;
        if (meson.EqualTo("K0Star"))   return kK0Star;
        if (meson.EqualTo("Proton"))   return kProton;
        std::cout << "Warning: particle \"" << meson << "\" not defined" << std:: endl;
        return kPi0;
    }


    // ****************************************************************************************************************
    //********************* global Variables **************************************************************************
    // ****************************************************************************************************************
    // PbPb 2.76TeV N_coll and T_AA
    // from https://arxiv.org/pdf/1301.4361.pdf
    // old values used for PbPb publications
    //                                             // 0-5,     5-10,   10-15,  15-20,  20-25,  25-30,  30-35,  35-40,  40-45,  45-50
    //     Double_t nCollPbPb2760GeVV0M5[20]       = {1684.4,  1316,   0,      0,      629.6,  483.7,  366.7,  273.4,  0,      0,
    //                                                 0,      0,      0,      0,      0,      0,      0,      0,      0,      0 };
    //     Double_t nCollPbPb2760GeVErrV0M5[20]    = { 190,    140,    0,      0,      62,     47,     35,     26,     0,      0,
    //                                                 0,      0,      0,      0,      0,      0,      0,      0,      0,      0 };
    //     Double_t nCollPbPb2760GeVV0M10[10]      = {1500, 921.2, (nCollPbPb2760GeVV0M5[4]+nCollPbPb2760GeVV0M5[5])/2, (nCollPbPb2760GeVV0M5[6]+nCollPbPb2760GeVV0M5[7])/2, 171.25,
    //                                                 84.28, 37.855, 15.575, 6.293, 0  };
    //     Double_t nCollPbPb2760GeVErrV0M10[10]   = { 165, 96, (nCollPbPb2760GeVErrV0M5[4]+nCollPbPb2760GeVΕρρV0M5[5])/2, (nCollPbPb2760GeVErrV0M5[6]+nCollPbPb2760GeVErrV0M5[7])/2, 16,
    //                                                 6.95, 2.85, 1.035, 0.325, 0  };
    //     Double_t nCollPbPb2760GeVV0M20[5]       = {1210.6, 438.4, 127.7, 26.71, 0};
    //     Double_t nCollPbPb2760GeVErrV0M20[5]    = {130.5, 42., 11, 2, 0};
    //
    //     Double_t nCollPbPb2760GeV2050           = 349.1;
    //     Double_t nCollPbPb2760GeVErr2050        = 51.;
    //     Double_t nCollPbPb2760GeV4080           = 77.205;
    //     Double_t nCollPbPb2760GeVErr4080        = 6.216;
    //     Double_t nCollPbPb2760GeV0040           = 740.;
    //     Double_t nCollPbPb2760GeV7590           = 8.219;
    //     Double_t nCollPbPb2760GeVErr7590        = 0.473;
    //     Double_t nCollPbPb2760GeV1030           = 11.6*64;
    //     Double_t nCollPbPb2760GeV30100          = 1.45*64;
    //
    //
    //     Double_t tAAPbPb2760GeVV0M5[20]         = { 26.32,  20.56,   0,      0,      0,      0,      0,      0,      0,      0,
    //                                                 0,      0,      0,      0,      0,      0,      0,      0,      0,      0 };
    //     Double_t tAAPbPb2760GeVErrV0M5[20]      = { 0.84224,    0.65792,    0,      0,      62,     47,     35,     26,     0,      0,
    //                                                 0,      0,      0,      0,      0,      0,      0,      0,      0,      0 };
    //     Double_t tAAPbPb2760GeVV0M10[10]        = { 23.44, 14.39, 0, 0, 0,   0, 0, 0, 0, 0  };
    //     Double_t tAAPbPb2760GeVErrV0M10[10]     = { 0.75008, 0.44609, 0, 0, 0,   0, 0, 0, 0, 0  };
    //     Double_t tAAPbPb2760GeVV0M20[10]        = { 18.915, 6.85, 1.996, 0.4174, 0};
    //     Double_t tAAPbPb2760GeVErrV0M20[10]     = { 0.5958225, 0.22605, 0.097804, 0.026296, 0};
    //
    //     Double_t tAAPbPb2760GeV2050             = 5.46;
    //     Double_t tAAPbPb2760GeVErr2050          = 0.195;

    // from https://cds.cern.ch/record/2636623/files/centrality%20determination%20note.pdf
    // updated 03.05.2019
                                            // 0-5,     5-10,   10-15,  15-20,  20-25,  25-30,  30-35,  35-40,  40-45,  45-50
    Double_t nCollPbPb2760GeVV0M5[20]       = { 1619,   1269,   1004,   791.1,  622.3,  485,    370.7,  279.4,  204.7,  148.4,
                                            // 50-55,   55-60,  60-65,  65-70,  70-75,  75-80,  80-85,  85-90,  90-95,  95-100
                                                104.4,  71.98,  48.15,  31.32,  19.69,  12.18,  7.244,  4.075,  2.173,  1.224 };
    Double_t nCollPbPb2760GeVErrV0M5[20]    = { 31,     27,     24,     19,     16,     14,     10,     9.2,    6.5,    5.7,
                                                4.2,    3.8,    2.9,    1.9,    0.97,   0.81,   0.48,   0.16,   0.11,   0.082 };
    Double_t nCollPbPb2760GeVV0M10[10]      = {1444, 897.7, 553.7, 325, 176.6,      88.21, 39.74, 15.96, 5.657, 1.709};
    Double_t nCollPbPb2760GeVErrV0M10[10]   = { 28, 21, 14, 9.7, 6,                 4, 2.4, 0.93, 0.31, 0.099};
    Double_t nCollPbPb2760GeVV0M20[5]       = {1171, 439.3, 132.4, 27.84, 3.682};
    Double_t nCollPbPb2760GeVErrV0M20[5]    = {24, 12, 4.9, 1.6, 0.16};

    // not updated
    Double_t nCollPbPb2760GeV2050           = 349.1;
    Double_t nCollPbPb2760GeVErr2050        = 51.;
    Double_t nCollPbPb2760GeV4080           = 77.205;
    Double_t nCollPbPb2760GeVErr4080        = 6.216;
    Double_t nCollPbPb2760GeV0040           = 740.;
    Double_t nCollPbPb2760GeV7590           = 8.219;
    Double_t nCollPbPb2760GeVErr7590        = 0.473;
    Double_t nCollPbPb2760GeV1030           = 11.6*64;
    Double_t nCollPbPb2760GeV30100          = 1.45*64;


    Double_t tAAPbPb2760GeVV0M5[20]         = { 26.2,  20.53,   16.25,  12.8,   10.07,  7.848,  5.998,  4.521,  3.313,  2.401,
                                                1.69,  1.165,   0.7792, 0.5067, 0.3186, 0.1972, 0.1172, 0.06594, 0.03517, 0.01981 };

    Double_t tAAPbPb2760GeVErrV0M5[20]      = { 0.27,  0.23,    0.23,   0.18,   0.16,   0.19,   0.15,   0.14,   0.088,  0.082,
                                                0.06,  0.056,   0.045,  0.031,  0.016,  0.013,  0.0078, 0.0028, 0.002,  0.0014 };
    Double_t tAAPbPb2760GeVV0M10[10]        = { 23.37, 14.53, 8.96, 5.259, 2.857,           1.427, 0.6431, 0.2582, 0.09153, 0.02765 };
    Double_t tAAPbPb2760GeVErrV0M10[10]     = { 0.2, 0.2, 0.17, 0.14, 0.084,                0.058, 0.038, 0.015, 0.0051, 0.0018  };
    Double_t tAAPbPb2760GeVV0M20[10]        = { 18.95, 7.109, 2.143, 0.4505, 0.05959};
    Double_t tAAPbPb2760GeVErrV0M20[10]     = { 0.19, 0.15, 0.07, 0.026, 0.0028};

    // not updated
    Double_t tAAPbPb2760GeV2050             = 5.46;
    Double_t tAAPbPb2760GeVErr2050          = 0.195;


    // PbPb 5.02TeV N_coll and T_AA
    //from https://cds.cern.ch/record/2636623/files/centrality%20determination%20note.pdf

    // NColl values for 10% V0M slices
    Double_t nCollPbPb5TeVV0M10[10]     = {1572, 973.4, 592.7, 343.8, 185.7,    91.41, 40.5, 16.12, 5.667, 1.708};
    Double_t nCollPbPb5TeVErrV0M10[10]  = {17.4, 11.3, 8.21, 5.76, 3.33,        2.11, 1.03, 0.341, 0.1, 0.0474};
    // NColl values for 20% V0M slices
    Double_t nCollPbPb5TeVV0M20[5]      = {1273, 468.2, 138.5, 28.31, 3.691};
    Double_t nCollPbPb5TeVErrV0M20[5]   = {14.1, 6.92, 2.7, 0.68, 0.0761};
    // NColl values for 20% V0M slices alternates
                                        // 10-30, 30-50, 50-70, 70-90
    Double_t nCollPbPb5TeVV0M20_2[5]      = {783.1, 264.8, 66.0, 10.9};
    Double_t nCollPbPb5TeVErrV0M20_2[5]   = {9.8, 4.51, 1.57, 0.221};
    // 0-5,     5-10,   10-15,  15-20,  20-25,  25-30,  30-35,  35-40,  40-45,  45-50
    Double_t nCollPbPb5TeVV0M5[20]      = { 1763,   1382,   1090,   857.3,  668.8,  516.6,  393,    294.6,  216.1,  155.2,
                                        // 50-55,   55-60,  60-65,  65-70,  70-75,  75-80,  80-85,  85-90,  90-95,  95-100
                                            108.6,  74.14,  49.18,  31.81,  20,     12.23,  7.198,  4.104,  2.172,  1.228 };
    Double_t nCollPbPb5TeVErrV0M5[20]   = { 19.4,   15.7,   12.3,   10.5,   9.21,   7.32,   6.51,   5.12,   3.83,   2.91,
                                            2.61,   1.74,   1.37,   0.763,  0.502,  0.253,  0.236,  0.0635, 0.0844, 0.0416 };
    Double_t nCollPbPb5TeVV0M5090       = 38.45;
    Double_t nCollPbPb5TeVErrV0M5090    = 0.9;

    // TAA values for 10% V0M slices
    Double_t tAAPbPb5TeVV0M10[10]       = {23.26, 14.4, 8.767, 5.086, 2.747,    1.352, 0.5992, 0.2385, 0.08383, 0.02527};
    Double_t tAAPbPb5TeVErrV0M10[10]    = {0.168, 0.126, 0.101, 0.0814, 0.0486, 0.0309, 0.0158, 0.00552, 0.00178, 0.000777};
    // TAA values for 20% V0M slices
                                          //0-20, 20-40, 40-60, 60-80, 80-100
    Double_t tAAPbPb5TeVV0M20[5]        = {18.83, 6.927, 2.049, 0.4188, 0.0546};
    Double_t tAAPbPb5TeVErrV0M20[5]     = {0.142, 0.0909, 0.0394, 0.0106, 0.00133};
    // TAA values alternate for 20% V0M slices
                                          // 10-30, 30-50, 50-70, 70-90 - not exact values: need to ask Alberica for precise numbers
    Double_t tAAPbPb5TeVV0M20_2[4]      = {11.58, 3.917, 0.9756, 0.1612};
    Double_t tAAPbPb5TeVErrV0M20_2[4]   = {0.1135, 0.0645, 0.02335, 0.00365};
    // 0-5,     5-10,   10-15,  15-20,  20-25,  25-30,  30-35,  35-40,  40-45,  45-50
    Double_t tAAPbPb5TeVV0M5[20]        = { 26.08,   20.44,   16.12,   12.68,  9.894,  7.641,  5.814,    4.358,  3.197,  2.296,
                                        // 50-55,   55-60,  60-65,  65-70,  70-75,  75-80,  80-85,  85-90,  90-95,  95-100
                                            1.607,  1.097,  0.7275,  0.4706,  0.2959,  0.1808,  0.1065,  0.06071,  0.03213,  0.01816 };
    Double_t tAAPbPb5TeVErrV0M5[20]     = { 0.176,  0.166,  0.135,   0.118,   0.11,    0.0942,  0.0917,  0.073,   0.0554,  0.0427,
                                            0.0378, 0.026,  0.0207, 0.012,   0.0079,  0.0041,  0.0037,  0.00126, 0.00132,  0.000673};
                                            // not exact values: need to ask Alberica for precise numbers
    Double_t tAAPbPb5TeVV0M5090        = 0.5684;
    Double_t tAAPbPb5TeVErrV0M5090     = 0.01234;

    // ******************************************************************************************
    // ************************ Set correct xSection for pp & pPb *******************************
    // ******************************************************************************************
    // general scale factors
    Double_t recalcBarn             = 1e12;         // NLO in pbarn!!!!
    Double_t factorToInel           = 1/1.12;       // this factor is multiplied with Raa and comes from trigger inelastic effiency for pp
    // pp 0.9 TeV
    Double_t xSection900GeVINEL     = 52.5*1e-3;
    Double_t xSection900GeV         = 47.78*1e-3;
    Double_t xSection900GeVV0AND    = 40.06*1e-3;
    Double_t xSection900GeVErrUp    = 2.39;
    Double_t xSection900GeVErrDown  = 1.86;
    // pp 2.76 TeV
    Double_t xSection2760GeV        = 55.416*1e-3;
    Double_t xSection2760GeVV0AND   = 47.73*1e-3;
    Double_t xSection2760GeVINEL    = 62.8*1e9;
    Double_t xSection2760GeVErr     = 3.9;
    // pp 5.023TeV
    Double_t xSection5023GeVV0AND   = 51.2*1e-3;    // from https://aliceinfo.cern.ch/ArtSubmission/sites/aliceinfo.cern.ch.ArtSubmission/files/draft/mgagliar/2016-Aug-22-paper_draft-vdmNote_5TeV.pdf
    Double_t xSection5023GeVV0ANDErr= 1.2;          // from https://aliceinfo.cern.ch/ArtSubmission/sites/aliceinfo.cern.ch.ArtSubmission/files/draft/mgagliar/2016-Aug-22-paper_draft-vdmNote_5TeV.pdf
    Double_t xSection5023GeVINEL    = 67.6*1e-3;    // from https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/stripath/2017-Jun-14-analysis_note-INEL_norm.pdf
    Double_t xSection5023GeVINELErr = 0.6;          // from https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/stripath/2017-Jun-14-analysis_note-INEL_norm.pdf
    // pp 7 TeV
    Double_t xSection7TeVINEL       = 73.2*1e-3;   // from https://arxiv.org/abs/1208.4968
    Double_t xSection7TeV           = 62.37*1e-3;  // from https://arxiv.org/abs/1208.4968
    Double_t xSection7TeVV0AND      = 54.31*1e-3;  // from https://arxiv.org/abs/1208.4968
    Double_t xSection7TeVErrUp      = 2.18;        // from https://arxiv.org/abs/1208.4968
    Double_t xSection7TeVErrDown    = 2.18;        // from https://arxiv.org/abs/1208.4968
    // pp 8 TeV
    Double_t xSection8TeVINEL       = 72.3*1e-3;   // from https://aliceinfo.cern.ch/Notes/node/665
    Double_t xSection8TeVINELErrUp  = 0.5;         // from https://aliceinfo.cern.ch/Notes/node/665
    Double_t xSection8TeVINELErrDown= 0.5;         // from https://aliceinfo.cern.ch/Notes/node/665
    Double_t xSection8TeVV0AND      = 55.8*1e-3;   // from https://aliceinfo.cern.ch/Notes/node/583
    Double_t xSection8TeVErrUp      = 1.45;        // from https://aliceinfo.cern.ch/Notes/node/583
    Double_t xSection8TeVErrDown    = 1.45;        // from https://aliceinfo.cern.ch/Notes/node/583
    Double_t xSection8TeVT0AND      = 25.5*1e-3;   // from https://aliceinfo.cern.ch/Notes/node/583
    Double_t xSection8TeVT0ErrUp    = 0.6;         // from https://aliceinfo.cern.ch/Notes/node/583
    Double_t xSection8TeVT0ErrDown  = 0.6;         // from https://aliceinfo.cern.ch/Notes/node/583
    // pp 13 TeV
    Double_t xSection13TeVV0AND     = 57.8*1e-3;   //
    Double_t xSection13TeVErrUp     = 1.27;         //
    Double_t xSection13TeVErrDown   = 1.27;         //
    Double_t xSection13TeVINEL      = 77.6*1e-3;   //
    Double_t xSection13TeVINELErr   = 1.0*1e-3;   //

    // pPb  5TeV
    // reference cent bins: https://twiki.cern.ch/twiki/bin/view/ALICE/PACentStudiesRun2 +
    // updated values for all cents 03.05.2019 according to: https://cds.cern.ch/record/2636623/files/centrality%20determination%20note.pdf

    Double_t xSection5023GeVINELpPb    = 67.6*1e-3;
    Double_t ncollpPb5023GeV           = 6.708;
    Double_t ncollErrpPb5023GeV        = 0.11;

    //   old values
    //                                                 // 0-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80, 80-90, 90-100
    //     Double_t nCollpPb5TeVBaseV0A10[10]      = { 13.285, 11.18, 10.5, 8.95, 7.38, 5.83, 4.42, 3.2, 2.23, 1.5};
    //     Double_t nCollpPbErr5TeVBaseV0A10[10]   = { 0.7178612504, 0.612, 3.9, 3.87, 3.72, 3.4, 2.94, 2.29, 1.59, 0.909};
    //     // 0-5, 5-10
    //     Double_t nCollpPb5TeVBaseV0A5[2]        = { 14.1, 12.47};
    //     Double_t nCollpPbErr5TeVBaseV0A5[2]     = { 0.747, 0.687};
    //     // 0-20, 20-40, 40-60, 60-80, 80-100
    //     Double_t nCollpPb5TeVBaseV0A20[5]       = { 12.2325, 9.089, 6.42, 3.884, 1.982};
    //     Double_t nCollpPbErr5TeVBaseV0A20[5]    = { 0.6653017029, 0.472, 0.272, 0.0974, 0.0331};
    //
    //     Double_t nCollpPb5TeVV0A60100           = 2.933;
    //     Double_t nCollpPbErr5TeVV0A60100        = 0.0612667663;
    //                                             // 0-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80, 80-90, 90-100
    //     Double_t nCollpPb5TeVBaseCL110[10]      = { 14.475, 11.87, 10.8, 9.03, 7.32, 5.67, 4.14, 2.93, 2.08, 1.45};
    //     Double_t nCollpPbErr5TeVBaseCL110[10]   = { 0.7494346839, 0.626, 3.53, 3.52, 3.38, 3.09, 2.62, 1.99, 1.4, 0.83};
    //                                             // 0-5, 5-10
    //     Double_t nCollpPb5TeVBaseCL15[2]        = { 15.51, 13.44};
    //     Double_t nCollpPbErr5TeVBaseCL15[2]     = { 0.804, 0.695};
    //                                             // 0-20, 20-40, 40-60, 60-80, 80-100
    //     Double_t nCollpPb5TeVBaseCL120[5]       = { 13.1725, 9.396, 6.169, 3.253, 1.6};
    //     Double_t nCollpPbErr5TeVBaseCL120[5]    = { 0.6883448722, 0.474, 0.251, 0.112, 0.0373};
    //     Double_t nCollpPb5TeVCL160100           = 2.4265;
    //     Double_t nCollpPbErr5TeVCL160100        = 0.0700557935;
    //     // 0-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80, 80-90, 90-100
    //     Double_t nCollpPb5TeVBaseHybrid10[10]   = { 13.2, 11.7, 10.8, 9.03, 7.32, 5.67, 4.14, 2.93, 2.08, 1.45};
    //     Double_t nCollpPbErr5TeVBaseHybrid10[10]= { 1.04095, 0.4446, 3.53, 3.52, 3.38, 3.09, 2.62, 1.99, 1.4, 0.83};
    //                                             // 0-5, 5-10
    //     Double_t nCollpPb5TeVBaseHybrid5[2]     = { 13.7, 12.7};
    //     Double_t nCollpPbErr5TeVBaseHybrid5[2]  = { 1.2056, 0.8763};
    //                                             // 0-20, 20-40, 40-60, 60-80, 80-100
    //     Double_t nCollpPb5TeVBaseHybrid20[5]    = { 12.45, 9.89, 6.94, 4.12, 2.12};
    //     Double_t nCollpPbErr5TeVBaseHybrid20[5] = { 0.7252125, 0.23736, 0.34006, 0.29664, 0.07844};
    //     Double_t nCollpPb5TeVHybrid60100        = 3.12;
    //     Double_t nCollpPbErr5TeVHybrid60100     = 0.17004;
    //
    //     Double_t tpPb5023GeV                    = 0.0983e3*(1/recalcBarn);
    //     Double_t tpPbErr5023GeV                 = 0.0017e3*(1/recalcBarn);

                                            // 0-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80, 80-90, 90-100
    Double_t nCollpPb5TeVBaseV0A10[10]      = { 13.29, 11.18, 9.759, 8.451, 7.075, 5.755, 4.437, 3.306, 2.365, 1.571};
    Double_t nCollpPbErr5TeVBaseV0A10[10]   = { 0.712, 0.612, 0.542, 0.391, 0.333, 0.222, 0.186, 0.0731, 0.0628, 0.0287};
                                            // 0-5, 5-10
    Double_t nCollpPb5TeVBaseV0A5[2]        = { 14.1, 12.47};
    Double_t nCollpPbErr5TeVBaseV0A5[2]     = { 0.747, 0.687};
                                            // 0-20, 20-40, 40-60, 60-80, 80-100
    Double_t nCollpPb5TeVBaseV0A20[5]       = { 12.23, 9.089, 6.42, 3.884, 1.982};
    Double_t nCollpPbErr5TeVBaseV0A20[5]    = { 0.673, 0.472, 0.272, 0.0974, 0.0331};

    Double_t nCollpPb5TeVV0A60100           = 2.93;
    Double_t nCollpPbErr5TeVV0A60100        = 0.1;

    // 0-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80, 80-90, 90-100
    Double_t nCollpPb5TeVBaseCL110[10]      = { 14.01, 11.62, 10.03, 8.54, 6.998, 5.476, 4.077, 2.914, 2.077, 1.441};
    Double_t nCollpPbErr5TeVBaseCL110[10]   = { 0.773, 0.607, 0.529, 0.42, 0.324, 0.241, 0.249, 0.208, 0.157, 0.105};
    // 0-5, 5-10
    Double_t nCollpPb5TeVBaseCL15[2]        = { 14.97, 13.05};
    Double_t nCollpPbErr5TeVBaseCL15[2]     = { 0.856, 0.743};
    // 0-20, 20-40, 40-60, 60-80, 80-100
    Double_t nCollpPb5TeVBaseCL120[5]       = { 12.85, 9.312, 6.273, 3.485, 1.744};
    Double_t nCollpPbErr5TeVBaseCL120[5]    = { 0.678, 0.458, 0.28, 0.187, 0.0863};
    Double_t nCollpPb5TeVCL160100           = 2.624;
    Double_t nCollpPbErr5TeVCL160100        = 0.128;

    // 0-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80, 80-90, 90-100 taking Ncoll mult
    Double_t nCollpPb5TeVBaseHybrid10[10]   = { 11.6, 10.7, 9.82, 8.76, 7.57, 6.17, 4.8, 3.55, 2.09, 1.58};
    Double_t nCollpPbErr5TeVBaseHybrid10[10]= { 0.8236, 0.4387, 0.1964, 0.15768, 0.3028, 0.38254, 0.36, 0.2556, 0.11495, 0.05688};
    // 0-5, 5-10
    Double_t nCollpPb5TeVBaseHybrid5[2]     = { 11.9, 11.3};
    Double_t nCollpPbErr5TeVBaseHybrid5[2]  = { 1.0234, 0.6215};
    // 0-20, 20-40, 40-60, 60-80, 80-100
    Double_t nCollpPb5TeVBaseHybrid20[5]    = { 11.1, 9.29, 6.87, 4.17, 2.03};
    Double_t nCollpPbErr5TeVBaseHybrid20[5] = { 0.6216, 0.22296, 0.35037, 0.29607, 0.06699};
    Double_t nCollpPb5TeVHybrid60100        = 3.1;
    Double_t nCollpPbErr5TeVHybrid60100     = 0.1674;

    Double_t tpPb5023GeV                    = 0.0983e3*(1/recalcBarn);
    Double_t tpPbErr5023GeV                 = 0.0017e3*(1/recalcBarn);

                                            // 0-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80, 80-90, 90-100
    Double_t tpPb5TeVBaseV0A10[10]          = { 0.1965, 0.1654, 0.1444, 0.125, 0.1047, 0.08513, 0.06564, 0.04891, 0.03499, 0.02324};
    Double_t tpPbErr5TeVBaseV0A10[10]       = { 0.0095, 0.00813, 0.00715, 0.00508, 0.00441, 0.0029, 0.00258, 0.00105, 0.000889, 0.00038};
                                            // 0-5, 5-10
    Double_t tpPb5TeVBaseV0A5[2]            = { 0.2085, 0.1845};
    Double_t tpPbErr5TeVBaseV0A5[2]         = { 0.00997, 0.00919};
                                            // 0-20, 20-40, 40-60, 60-80, 80-100
    Double_t tpPb5TeVBaseV0A20[5]           = { 0.181, 0.1345, 0.09496, 0.05745, 0.02932};
    Double_t tpPbErr5TeVBaseV0A20[5]        = { 0.00898, 0.00624, 0.00357, 0.00127, 0.000468};
    Double_t tpPb5TeVBaseV0A60100           = 0.04334;
    Double_t tpPbErr5TeVBaseV0A60100        = 0.00138;

    // 0-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80, 80-90, 90-100
    Double_t tpPb5TeVBaseCL110[10]          = { 0.2073, 0.1719, 0.1483, 0.1263, 0.1035, 0.081, 0.06031, 0.04311, 0.03073, 0.02132};
    Double_t tpPbErr5TeVBaseCL110[10]       = { 0.0106, 0.00815, 0.00707, 0.00551, 0.00426, 0.00329, 0.00353, 0.00291, 0.00224, 0.00154};
                                            // 0-5, 5-10
    Double_t tpPb5TeVBaseCL15[2]            = { 0.2214, 0.193};
    Double_t tpPbErr5TeVBaseCL15[2]         = { 0.0118, 0.0102};
                                            // 0-20, 20-40, 40-60, 60-80, 80-100
    Double_t tpPb5TeVBaseCL120[5]           = { 0.1901, 0.1377, 0.0928, 0.05156, 0.02579};
    Double_t tpPbErr5TeVBaseCL120[5]        = { 0.00913, 0.00599, 0.0037, 0.00255, 0.00116};
    Double_t tpPb5TeVBaseCL160100           = 0.03881;
    Double_t tpPbErr5TeVBaseCL160100        = 0.00181;


    // 0-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80, 80-90, 90-100 taking TAA mult
    Double_t tpPb5TeVBaseHybrid10[10]       = { 0.172, 0.158, 0.145, 0.13, 0.112, 0.0913, 0.071, 0.0525, 0.0309, 0.0234};
    Double_t tpPbErr5TeVBaseHybrid10[10]    = { 0.012212, 0.006478, 0.003045, 0.00247, 0.004592, 0.0056606, 0.005325, 0.00378, 0.0016995, 0.0008424};
                                            // 0-5, 5-10
    Double_t tpPb5TeVBaseHybrid5[2]         = { 0.176, 0.167};
    Double_t tpPbErr5TeVBaseHybrid5[2]      = { 0.015136, 0.009352};
                                            // 0-20, 20-40, 40-60, 60-80, 80-100
    Double_t tpPb5TeVBaseHybrid20[5]        = { 0.164, 0.137, 0.102, 0.0617, 0.03};
    Double_t tpPbErr5TeVBaseHybrid20[5]     = { 0.009348, 0.003288, 0.005202, 0.0043807, 0.00102};
    Double_t tpPb5TeVBaseHybrid60100        = 0.0459;
    Double_t tpPbErr5TeVBaseHybrid60100     = 0.0024786;

    // basic function to convert cutNumber to integere
    Int_t CutNumberToInteger(TString cutNumber){
      char tmpChar = cutNumber(0);
      Int_t tmpInt = tmpChar - '0' - 39;

      if(tmpInt >= 10 && tmpInt <= 35) return tmpInt;
      else return cutNumber.Atoi();
    }

    //************************************************************************************
    //********************* Separate cut numbers, old version ****************************
    //************************************************************************************
    void ReturnSeparatedCutNumber(  TString cutSel,
                                    TString& gammaCutNumber,
                                    TString& electronCutNumber,
                                    TString& mesonCutNumber,
                                    Bool_t kDalitz=kFALSE){

        TObjArray *arr;
        arr = cutSel.Tokenize("_");
        TObjString* objstrGamma;
        TObjString* objstrElectron;
        TObjString* objstrMeson;

        if (kDalitz){
            objstrGamma = (TObjString*)arr->At(0);
            objstrElectron = (TObjString*)arr->At(1);
            objstrMeson = (TObjString*)arr->At(2);

            gammaCutNumber= objstrGamma->GetString();
            electronCutNumber = objstrElectron->GetString();
            mesonCutNumber = objstrMeson->GetString();
        } else {
            objstrGamma = (TObjString*)arr->At(0);
            objstrMeson = (TObjString*)arr->At(1);

            gammaCutNumber= objstrGamma->GetString();
            mesonCutNumber = objstrMeson->GetString();
        }
        cout << cutSel.Data() << "\t" << gammaCutNumber.Data() << "\t" << electronCutNumber.Data() << "\t" << mesonCutNumber.Data() << endl;
        return;
    }

    //************************************************************************************
    //********************* Separate cut numbers, new version ****************************
    //** separates the possible combinations according to the mode into the             **
    //** respective eventCutNumber, gammaCutNumber (PCM), clusterCutNumber (EMCal/PHOS),**
    //** electronCutNumber (Dalitz) and mesonCutNumber                                  **
    //************************************************************************************
    void ReturnSeparatedCutNumberAdvanced(  TString cutSel,
                                            TString& eventCutNumber,
                                            TString& gammaCutNumber,
                                            TString& clusterCutNumber,
                                            TString& electronCutNumber,
                                            TString& mesonCutNumber,
                                            Int_t type=0){

        TObjArray *arr;
        arr = cutSel.Tokenize("_");
        TObjString* objstrEvent;
        TObjString* objstrGamma;
        TObjString* objstrCluster;
        TObjString* objstrElectron;
        TObjString* objstrMeson;

        // Switch for heavy meson modes
        Int_t pos = 0;
        if(type>=100) {
            ++pos;
            type -= 100;
        }
        // Read cut numbers based on given mode/type
        if (type == 0){ // PCM-PCM
            objstrEvent         = (TObjString*)arr->At(pos);
            objstrGamma         = (TObjString*)arr->At(1+pos);
            objstrMeson         = (TObjString*)arr->At(2+pos);

            eventCutNumber      = objstrEvent->GetString();
            gammaCutNumber      = objstrGamma->GetString();
            mesonCutNumber      = objstrMeson->GetString();

        } else if (type == 1){ //PCM dalitz
            objstrEvent         = (TObjString*)arr->At(pos);
            objstrGamma         = (TObjString*)arr->At(1+pos);
            objstrElectron      = (TObjString*)arr->At(2+pos);
            objstrMeson         = (TObjString*)arr->At(3+pos);

            eventCutNumber      = objstrEvent->GetString();
            gammaCutNumber      = objstrGamma->GetString();
            electronCutNumber   = objstrElectron->GetString();
            mesonCutNumber      = objstrMeson->GetString();

        } else if (type == 2 || type == 13 || type == 14){ //PCM-EMCal (PCM-DCal), PCM EMCal+DCal
            objstrEvent         = (TObjString*)arr->At(pos);
            objstrGamma         = (TObjString*)arr->At(1+pos);
            objstrCluster       = (TObjString*)arr->At(2+pos);
            objstrMeson         = (TObjString*)arr->At(3+pos);

            eventCutNumber      = objstrEvent->GetString();
            gammaCutNumber      = objstrGamma->GetString();
            clusterCutNumber    = objstrCluster->GetString();
            mesonCutNumber      = objstrMeson->GetString();
        } else if (type == 3){ //PCM-PHOS
            objstrEvent         = (TObjString*)arr->At(pos);
            objstrGamma         = (TObjString*)arr->At(1+pos);
            objstrCluster       = (TObjString*)arr->At(2+pos);
            objstrMeson         = (TObjString*)arr->At(3+pos);

            eventCutNumber      = objstrEvent->GetString();
            gammaCutNumber      = objstrGamma->GetString();
            clusterCutNumber    = objstrCluster->GetString();
            mesonCutNumber      = objstrMeson->GetString();
        } else if (type == 4 || type == 12 || type == 15 || type == 200){ //EMCal-EMCal (DCal-DCal), EMCal+DCal - EMCal+DCal, TriggerQA
            objstrEvent         = (TObjString*)arr->At(pos);
            objstrCluster       = (TObjString*)arr->At(1+pos);
            objstrMeson         = (TObjString*)arr->At(2+pos);

            eventCutNumber      = objstrEvent->GetString();
            clusterCutNumber    = objstrCluster->GetString();
            if (type != 200) mesonCutNumber      = objstrMeson->GetString();

            //temporary fix for train outputs from vAN20170525 to vAN20170601 for ConversionMesonCuts
            if(mesonCutNumber.Length()==17){
              TString tempMesonCutNumber = mesonCutNumber;
              tempMesonCutNumber.Remove(0,14);
              tempMesonCutNumber.Chop();
              Int_t tmpIntMesonCutNumber = tempMesonCutNumber.Atoi();
              char tmpChar = tmpIntMesonCutNumber+39+'0';
              mesonCutNumber.Replace(14,2,tmpChar);
            }
        } else if (type == 5 || type == -5){ //PHOS-PHOS
            objstrEvent         = (TObjString*)arr->At(pos);
            objstrCluster       = (TObjString*)arr->At(1+pos);
            objstrMeson         = (TObjString*)arr->At(2+pos);

            eventCutNumber      = objstrEvent->GetString();
            clusterCutNumber    = objstrCluster->GetString();
            mesonCutNumber      = objstrMeson->GetString();
        } else if (type == 6 ){ //Dalitz-EMCal
            objstrEvent         = (TObjString*)arr->At(pos);
            objstrGamma         = (TObjString*)arr->At(1+pos);
            objstrElectron      = (TObjString*)arr->At(2+pos);
            objstrCluster       = (TObjString*)arr->At(3+pos);
            objstrMeson         = (TObjString*)arr->At(4+pos);

            eventCutNumber      = objstrEvent->GetString();
            gammaCutNumber      = objstrGamma->GetString();
            clusterCutNumber    = objstrCluster->GetString();
            electronCutNumber   = objstrElectron->GetString();
            mesonCutNumber      = objstrMeson->GetString();
        } else if ( type == 7 ){ //Dalitz-PHOS
            objstrEvent         = (TObjString*)arr->At(pos);
            objstrGamma         = (TObjString*)arr->At(1+pos);
            objstrElectron      = (TObjString*)arr->At(2+pos);
            objstrCluster       = (TObjString*)arr->At(3+pos);
            objstrMeson         = (TObjString*)arr->At(4+pos);

            eventCutNumber      = objstrEvent->GetString();
            gammaCutNumber      = objstrGamma->GetString();
            clusterCutNumber    = objstrCluster->GetString();
            electronCutNumber   = objstrElectron->GetString();
            mesonCutNumber      = objstrMeson->GetString();
        } else if ( type == 10 ){ // EMCal merged
            objstrEvent         = (TObjString*)arr->At(pos);
            objstrGamma         = (TObjString*)arr->At(1+pos);
            objstrCluster       = (TObjString*)arr->At(2+pos);
            objstrMeson         = (TObjString*)arr->At(3+pos);

            eventCutNumber      = objstrEvent->GetString();
            gammaCutNumber      = objstrGamma->GetString();
            clusterCutNumber    = objstrCluster->GetString();
            mesonCutNumber      = objstrMeson->GetString();
        } else if ( type == 11 ){ // PHOS merged
            objstrEvent         = (TObjString*)arr->At(pos);
            objstrGamma         = (TObjString*)arr->At(1+pos);
            objstrCluster       = (TObjString*)arr->At(2+pos);
            objstrMeson         = (TObjString*)arr->At(3+pos);

            eventCutNumber      = objstrEvent->GetString();
            gammaCutNumber      = objstrGamma->GetString();
            clusterCutNumber    = objstrCluster->GetString();
            mesonCutNumber      = objstrMeson->GetString();
        // } else if (type == 12){ // flow
        //     objstrEvent         = (TObjString*)arr->At(0);
        //     objstrGamma         = (TObjString*)arr->At(1);
        //
        //     eventCutNumber      = objstrEvent->GetString();
        //     gammaCutNumber      = objstrGamma->GetString();
        }

        cout << cutSel.Data() << "\t" << eventCutNumber.Data() << "\t" << gammaCutNumber.Data() << "\t" <<  clusterCutNumber.Data() << "\t" <<electronCutNumber.Data() << "\t" << mesonCutNumber.Data() << endl;
        return;
    }

    //*****************************************************************************************
    //********************* Separate cut numbers, PiPlPiMiPiZero ******************************
    //** separates the possible combinations according to the mode into the respective       **
    //** typeCutNumber, eventCutNumber, gammaCutNumber (PCM), clusterCutNumber (EMCal/PHOS), **
    //** pionCutNumber, neutralPionCutNumber and mesonCutNumber                              **
    //*****************************************************************************************
    Int_t ReturnSeparatedCutNumberPiPlPiMiPiZero(TString cutSel,
                                                 TString& typeCutNumber,
                                                 TString& eventCutNumber,
                                                 TString& gammaCutNumber,
                                                 TString& clusterCutNumber,
                                                 TString& pionCutNumber,
                                                 TString& neutralPionCutNumber,
                                                 TString& mesonCutNumber,
                                                 Bool_t runningNewTask = kFALSE){

        TObjArray *arr;
        arr = cutSel.Tokenize("_");
        TObjString* objstrType;
        TObjString* objstrEvent;
        TObjString* objstrGamma;
        TObjString* objstrCluster;
        TObjString* objstrPion;
        TObjString* objstrNeutralPion;
        TObjString* objstrMeson;
        Int_t mode = -1;
        objstrType  = (TObjString*)arr->At(0);
        typeCutNumber  = objstrType->GetString();

        if (typeCutNumber.CompareTo("0") == 0){ // PCM-PCM
            objstrType  = (TObjString*)arr->At(0);
            objstrEvent = (TObjString*)arr->At(1);
            objstrGamma = (TObjString*)arr->At(2);
            objstrNeutralPion = (TObjString*)arr->At(3);
            objstrPion = (TObjString*)arr->At(4);
            objstrMeson = (TObjString*)arr->At(5);

            typeCutNumber  = objstrType->GetString();
            eventCutNumber= objstrEvent->GetString();
            gammaCutNumber= objstrGamma->GetString();
            pionCutNumber = objstrPion->GetString();
            neutralPionCutNumber = objstrNeutralPion->GetString();
            mesonCutNumber = objstrMeson->GetString();

            mode = 40;
            if(runningNewTask) mode = 60;
        } else if (typeCutNumber.CompareTo("1") == 0){ //PCM-calo
            objstrType  = (TObjString*)arr->At(0);
            objstrEvent = (TObjString*)arr->At(1);
            objstrGamma = (TObjString*)arr->At(2);
            objstrCluster = (TObjString*)arr->At(3);
            objstrNeutralPion = (TObjString*)arr->At(4);
            objstrPion = (TObjString*)arr->At(5);
            objstrMeson = (TObjString*)arr->At(6);

            typeCutNumber  = objstrType->GetString();
            eventCutNumber= objstrEvent->GetString();
            gammaCutNumber= objstrGamma->GetString();
            clusterCutNumber = objstrCluster->GetString();
            pionCutNumber = objstrPion->GetString();
            neutralPionCutNumber = objstrNeutralPion->GetString();
            mesonCutNumber = objstrMeson->GetString();

            // Check wich calo was used
            TString firstLetter(clusterCutNumber(0,1));
            if (firstLetter.CompareTo("1") == 0){        // EMCAL was used as calo
                mode = 41;
                if(runningNewTask) mode = 61;
            } else if (firstLetter.CompareTo("2") == 0){ // PHOS was used as calo
                mode = 42;
                if(runningNewTask) mode = 62;
            } else if (firstLetter.CompareTo("3") == 0){ // DCAL was used as calo
                mode = 43;
                if(runningNewTask) mode = 63;
            } else {
                mode = 41;
                if(runningNewTask) mode = 61;
            }
        }  else if (typeCutNumber.CompareTo("2") == 0){ // CALO-CALO
            objstrType  = (TObjString*)arr->At(0);
            objstrEvent = (TObjString*)arr->At(1);
            objstrCluster = (TObjString*)arr->At(2);
            objstrNeutralPion = (TObjString*)arr->At(3);
            objstrPion = (TObjString*)arr->At(4);
            objstrMeson = (TObjString*)arr->At(5);

            typeCutNumber  = objstrType->GetString();
            eventCutNumber= objstrEvent->GetString();
            clusterCutNumber = objstrCluster->GetString();
            pionCutNumber = objstrPion->GetString();
            neutralPionCutNumber = objstrNeutralPion->GetString();
            mesonCutNumber = objstrMeson->GetString();

            // Check wich calo was used
            TString firstLetter(clusterCutNumber(0,1));
            if (firstLetter.CompareTo("1") == 0){        // EMCAL was used as calo
                mode = 44;
                if(runningNewTask) mode = 64;
            } else if (firstLetter.CompareTo("2") == 0){ // PHOS was used as calo
                mode = 45;
                if(runningNewTask) mode = 65;
            } else if (firstLetter.CompareTo("3") == 0){ // DCAL was used as calo
                mode = 46;
                if(runningNewTask) mode = 66;
            } else {
                mode = 44;
                if(runningNewTask) mode = 64;
            }
        }

        cout << cutSel.Data() << "\t" << typeCutNumber.Data() << "\t" << eventCutNumber.Data() << "\t" << gammaCutNumber.Data() << "\t" <<  clusterCutNumber.Data() << "\t" <<pionCutNumber.Data() << "\t" <<neutralPionCutNumber.Data() << "\t" << mesonCutNumber.Data() << endl;
        return mode;
    }

    //************************************************************************************
    //********************* get number of events for PCM/calo analysis *******************
    //************************************************************************************
    Int_t GetNEvents (  TH1* histo,
                        Bool_t doCout=kTRUE){
        //cout<<"Debug, ConversionFunctionsBasicsAndLabeling.h, Line: "<<__LINE__<<"; histo->GetNbinsX(): "<<histo->GetNbinsX()<<endl;
        if (!histo) cout << "NO EVENT HISTO" << endl;
        if(histo->GetNbinsX()==11){
            if(histo->GetEntries()-histo->GetBinContent(5)-histo->GetBinContent(7)-histo->GetBinContent(4)==0) return 0;
            Int_t nEvents = histo->GetBinContent(1)+(histo->GetBinContent(1)/(histo->GetBinContent(1)+histo->GetBinContent(5)))*histo->GetBinContent(6);
            Int_t nEventsMB = histo->GetEntries()-histo->GetBinContent(4) -histo->GetBinContent(8)-histo->GetBinContent(9);
            for (Int_t i = 1; i<12; i++ ){
                if(doCout) cout << histo->GetBinContent(i) << "\t";
            }
            if(doCout) cout << nEventsMB  << endl;
            for (Int_t i = 1; i<12; i++ ){
                    if(doCout) cout << histo->GetXaxis()->GetBinLabel(i) << "\t" << histo->GetBinContent(i)/nEventsMB << "\n";
            }
            if(doCout) cout  << endl;
            // BinContent 1 - good events
            // BinContent 2 - centrality not selected
            // BinContent 3 - MC corrupt
            // BinContent 4 - no Trigger Bit
            // BinContent 5 - Zvertex-position,
            // BinContent 6 - no Contributors to vtx
            // BinContent 7 - PileUp
            // BinContent 8 - no SDD
            // BinContent 9 - no V0AND
            if(doCout)cout <<"nEvents new: "<< nEvents << "\t nEvents old: "<< histo->GetEntries()-histo->GetBinContent(5)-histo->GetBinContent(7)-histo->GetBinContent(4)<< endl;
            return nEvents;
        }else if(histo->GetNbinsX()==12){
            if(histo->GetEntries()-histo->GetBinContent(5)-histo->GetBinContent(7)-histo->GetBinContent(12)-histo->GetBinContent(4)==0) return 0;
            Int_t nEvents = histo->GetBinContent(1)+(histo->GetBinContent(1)/(histo->GetBinContent(1)+histo->GetBinContent(5)))*histo->GetBinContent(6);
            Int_t nEventsMB = histo->GetEntries()-histo->GetBinContent(4) -histo->GetBinContent(8)-histo->GetBinContent(9)-histo->GetBinContent(2);
            for (Int_t i = 1; i<13; i++ ){
                if(doCout) cout << histo->GetBinContent(i) << "\t";
            }
            if(doCout) cout << nEventsMB  << endl;
            for (Int_t i = 1; i<13; i++ ){
                    if(doCout) cout << histo->GetXaxis()->GetBinLabel(i) << "\t" << histo->GetBinContent(i)/nEventsMB << "\n";
            }
            if(doCout) cout << "accepted \t" << (Float_t)nEvents/nEventsMB << endl;
            if(doCout) cout << endl;
            // BinContent 1 - good events
            // BinContent 2 - centrality not selected
            // BinContent 3 - MC corrupt
            // BinContent 4 - no Trigger Bit
            // BinContent 5 - Zvertex-position,
            // BinContent 6 - no Contributors to vtx
            // BinContent 7 - PileUp
            // BinContent 8 - no SDD
            // BinContent 9 - no V0AND
            // BinContent 12 - SPD cluster vs tracklets
            if(doCout)cout <<"nEvents new: "<< nEvents <<  endl;
            return nEvents;
        }else if(histo->GetNbinsX()==13){
            if(histo->GetEntries()-histo->GetBinContent(5)-histo->GetBinContent(7)-histo->GetBinContent(12)-histo->GetBinContent(4)-histo->GetBinContent(13)==0) return 0;
            Int_t nEvents = histo->GetBinContent(1)+(histo->GetBinContent(1)/(histo->GetBinContent(1)+histo->GetBinContent(5)))*histo->GetBinContent(6);
            Int_t nEventsMB = histo->GetEntries()-histo->GetBinContent(4) -histo->GetBinContent(8)-histo->GetBinContent(9)-histo->GetBinContent(2);
            for (Int_t i = 1; i<14; i++ ){
                if(doCout) cout << histo->GetBinContent(i) << "\t";
            }
            if(doCout) cout << nEventsMB  << endl;
            for (Int_t i = 1; i<14; i++ ){
                    if(doCout) cout << histo->GetXaxis()->GetBinLabel(i) << "\t" << histo->GetBinContent(i)/nEventsMB << "\n";
            }
            if(doCout) cout << "accepted \t" << (Float_t)nEvents/nEventsMB << endl;
            if(doCout) cout << endl;
            // BinContent 1 - good events
            // BinContent 2 - centrality not selected
            // BinContent 3 - MC corrupt
            // BinContent 4 - no Trigger Bit
            // BinContent 5 - Zvertex-position,
            // BinContent 6 - no Contributors to vtx
            // BinContent 7 - PileUp
            // BinContent 8 - no SDD
            // BinContent 9 - no V0AND
            // BinContent 12 - SPD cluster vs tracklets
            // BinContent 13 - Out of bunch pileup
            if(doCout)cout <<"nEvents new: "<< nEvents <<  endl;
            return nEvents;
        }else if(histo->GetNbinsX()==14){
            if(histo->GetEntries()-histo->GetBinContent(5)-histo->GetBinContent(7)-histo->GetBinContent(12)-histo->GetBinContent(4)-histo->GetBinContent(13)-histo->GetBinContent(14)==0) return 0;
            Int_t nEvents = histo->GetBinContent(1)+(histo->GetBinContent(1)/(histo->GetBinContent(1)+histo->GetBinContent(5)))*histo->GetBinContent(6);
            Int_t nEventsMB = histo->GetEntries()-histo->GetBinContent(4) -histo->GetBinContent(8)-histo->GetBinContent(9)-histo->GetBinContent(2);
            for (Int_t i = 1; i<15; i++ ){
                if(doCout) cout << "("<<i<<")" << histo->GetBinContent(i) << "\t";
            }
            if(doCout) cout <<endl<< "nEvents: "<< nEvents <<"; nEventsMB:" << nEventsMB  << endl;
            for (Int_t i = 1; i<15; i++ ){
                    if(doCout) cout << "("<<i<<")"<< histo->GetXaxis()->GetBinLabel(i) << "\t" << histo->GetBinContent(i)/nEventsMB << "\n";
            }
            if(doCout) cout << "accepted \t" << (Float_t)nEvents/nEventsMB << endl;
            if(doCout) cout << endl;
            // BinContent 1 - good events
            // BinContent 2 - centrality not selected
            // BinContent 3 - MC corrupt
            // BinContent 4 - no Trigger Bit
            // BinContent 5 - Zvertex-position,
            // BinContent 6 - no Contributors to vtx
            // BinContent 7 - PileUp
            // BinContent 8 - no SDD
            // BinContent 9 - no V0AND
            // BinContent 12 - SPD cluster vs tracklets
            // BinContent 13 - Out of bunch pileup
            // BinContent 14 - Pileup V0M-TPCout
            if(doCout)cout <<"nEvents new: "<< nEvents <<  endl;
            return nEvents;
        }else if(histo->GetNbinsX()==15){
            if(histo->GetEntries()-histo->GetBinContent(5)-histo->GetBinContent(7)-histo->GetBinContent(12)-histo->GetBinContent(4)-histo->GetBinContent(13)-histo->GetBinContent(14)-histo->GetBinContent(15)==0) return 0;
            Int_t nEvents = histo->GetBinContent(1)+(histo->GetBinContent(1)/(histo->GetBinContent(1)+histo->GetBinContent(5)))*histo->GetBinContent(6);
            Int_t nEventsMB = histo->GetEntries()-histo->GetBinContent(4) -histo->GetBinContent(8)-histo->GetBinContent(9)-histo->GetBinContent(2);
            for (Int_t i = 1; i<16; i++ ){
                if(doCout) cout << histo->GetBinContent(i) << "\t";
            }
            if(doCout) cout << nEventsMB  << endl;
            for (Int_t i = 1; i<16; i++ ){
                    if(doCout) cout << histo->GetXaxis()->GetBinLabel(i) << "\t" << histo->GetBinContent(i)/nEventsMB << "\n";
            }
            if(doCout) cout << "accepted \t" << (Float_t)nEvents/nEventsMB << endl;
            if(doCout) cout << endl;
            // BinContent 1 - good events
            // BinContent 2 - centrality not selected
            // BinContent 3 - MC corrupt
            // BinContent 4 - no Trigger Bit
            // BinContent 5 - Zvertex-position,
            // BinContent 6 - no Contributors to vtx
            // BinContent 7 - PileUp
            // BinContent 8 - no SDD
            // BinContent 9 - no V0AND
            // BinContent 12 - SPD cluster vs tracklets
            // BinContent 13 - Out of bunch pileup
            // BinContent 14 - Pileup V0M-TPCout
            // BinContent 15 - Sphericity
            if(doCout)cout <<"nEvents new: "<< nEvents <<  endl;
            return nEvents;
        }else{
            cout << "ERROR: GetNEvents, dimension of histogram not known! Returning 0...!" << endl;
            return 0;
        }
    }

    //************************************************************************************
    //********************* get number of events for PCM/calo analysis *******************
    //************************************************************************************
    Double_t GetMissMCEventFrac (  TH1* histo ){
        if (!histo) cout << "NO EVENT HISTO" << endl;
        if(histo->GetNbinsX()==11){
            if(histo->GetEntries()-histo->GetBinContent(5)-histo->GetBinContent(7)-histo->GetBinContent(4)==0) return 0;
            Int_t nEvents = histo->GetBinContent(1)+(histo->GetBinContent(1)/(histo->GetBinContent(1)+histo->GetBinContent(5)))*histo->GetBinContent(6);
            Int_t nEventsMB = histo->GetEntries()-histo->GetBinContent(4) -histo->GetBinContent(8)-histo->GetBinContent(9);
            Double_t missEventFrac = 0;
            if (nEventsMB > 0)
                missEventFrac = histo->GetBinContent(3)/nEventsMB;
            return missEventFrac;
        }else if(histo->GetNbinsX()==12){
            if(histo->GetEntries()-histo->GetBinContent(5)-histo->GetBinContent(7)-histo->GetBinContent(12)-histo->GetBinContent(4)==0) return 0;
            Int_t nEvents = histo->GetBinContent(1)+(histo->GetBinContent(1)/(histo->GetBinContent(1)+histo->GetBinContent(5)))*histo->GetBinContent(6);
            Int_t nEventsMB = histo->GetEntries()-histo->GetBinContent(4) -histo->GetBinContent(8)-histo->GetBinContent(9);
            Double_t missEventFrac = 0;
            if (nEventsMB > 0)
                missEventFrac = histo->GetBinContent(3)/nEventsMB;
            return missEventFrac;
        }else if(histo->GetNbinsX()==13){
            if(histo->GetEntries()-histo->GetBinContent(5)-histo->GetBinContent(7)-histo->GetBinContent(12)-histo->GetBinContent(13)-histo->GetBinContent(4)==0) return 0;
            Int_t nEvents = histo->GetBinContent(1)+(histo->GetBinContent(1)/(histo->GetBinContent(1)+histo->GetBinContent(5)))*histo->GetBinContent(6);
            Int_t nEventsMB = histo->GetEntries()-histo->GetBinContent(4) -histo->GetBinContent(8)-histo->GetBinContent(9);
            Double_t missEventFrac = 0;
            if (nEventsMB > 0)
                missEventFrac = histo->GetBinContent(3)/nEventsMB;
            return missEventFrac;
        }else{
            cout << "ERROR: GetNEvents, dimension of histogram not known! Returning 0...!" << endl;
            return 0;
        }
    }

    //************************************************************************************
    //************************ return latex writing of meson name ************************
    //************************************************************************************
    const char* ReturnMesonPlainString ( Mesons meson ){
        switch(meson) {
            case kPi0      : return "Pi0";
            case kEta      : return "Eta";
            case kEtaPrime : return "EtaPrime";
            case kOmega    : return "Omega";
            case kCKaon    : return "CKaon";
            case kCPion    : return "CPion";
            case kLambda   : return "Lambda";
            case kRho0     : return "Rho0";
            case kPhi      : return "Phi";
            case kK0Star   : return "K0Star";
            case kProton   : return "Proton";
            default :
                std::cout << "ERROR: meson enum " << meson << " not defined" << std::endl;
                return "";
        }
    }

    //************************************************************************************
    //**** Decodes from the mode the respective reco process and return correct label ****
    //************************************************************************************
    TString ReturnTextReconstructionProcess(Int_t mode){
        if( mode>=100 ) mode -= 100;
        switch (mode){
            case 0:
                return "PCM";
            case 1:
                return "PCM-#gamma^{*}#gamma";
            case 2:
                return "PCM-EMC";
            case 3:
                return "PCM-PHOS";
            case 4:
                return "EMC";
            case 5: case -5:
                return "PHOS";
            case 6:
                return "EMC-#gamma^{*}#gamma";
            case 7:
                return "PHOS-#gamma^{*}#gamma";
            case 10:
                return "mEMC";
            case 11:
                return "mPHOS";
            case 12:
                return "DMC";
            case 13:
                return "PCM-DMC";
            case 14:
                return "PCM-EDC";
            case 15:
                return "EDC";
            case 20: case 21: case 22: case 23:
                return "Comb";
            // Cases added for omega analysis
            case 40: case 60:
                return "#pi^{0} rec w/ PCM";
            case 41: case 61:
                return "#pi^{0} rec w/ PCM-EMC";
            case 42: case 62:
                return "#pi^{0} rec w/ PCM, PHOS";
            case 43: case 63:
                return "#pi^{0} rec w/ PCM, DCAL";
            case 44: case 64:
                return "#pi^{0} rec w/ EMC";
            case 45: case 65:
                return "#pi^{0} rec w/ PHOS";
            case 46: case 66:
                return "#pi^{0} rec w/ DCAL, DCAL";
            case 47: case 67:
                return "#pi^{0} rec w/ PCM, DALITZ";
            case 48: case 68:
                return "#pi^{0} rec w/ EMCAL, DALITZ";
            case 49: case 69:
                return "#pi^{0} rec w/ PHOS, DALITZ";
            case 50: case 70:
                return "#pi^{0} rec w/ DCAL, DALITZ";
            default:
                return "not known";

        }
    }

    //************************************************************************************
    //**** Decodes from the mode the respective reco process and return correct label.****
    //**** This method is used for filename text.                                     ****
    //************************************************************************************
    TString ReturnTextReconstructionProcessWrite(Int_t mode){
        if( mode>=100 ) mode -= 100;
        switch (mode){
            case 0:
                return "PCM";
            case 1:
                return "PCM-Dal";
            case 2:
                return "PCM-EMC";
            case 3:
                return "PCM-PHOS";
            case 4:
                return "EMC";
            case 5: case -5:
                return "PHOS";
            case 6:
                return "EMC-Dal";
            case 7:
                return "PHOS-Dal";
            case 10:
                return "mEMC";
            case 11:
                return "mPHOS";
            case 14:
                return "PCM-EDC";
            case 15:
                return "EDC";
            case 20: case 21: case 22: case 23:
                return "Comb";
            // Cases added for omega analysis
            case 40: case 60:
                return "Pi0PCM";
            case 41: case 61:
                return "Pi0PCM-EMC";
            case 42: case 62:
                return "Pi0PCM-PHOS";
            case 43: case 63:
                return "Pi0PCM-DCAL";
            case 44: case 64:
                return "Pi0EMC";
            case 45: case 65:
                return "Pi0PHOS";
            case 46: case 66:
                return "Pi0DCAL";
            case 47: case 67:
                return "Pi0PCM-Dal";
            case 48: case 68:
                return "Pi0EMC-Dal";
            case 49: case 69:
                return "Pi0PHOS-Dal";
            case 50: case 70:
                return "Pi0DCAL-Dal";
            default:
                return "not known";
        }
    }

    //************************************************************************************
    //************************ return latex writing of meson name ************************
    //************************************************************************************
    TString ReturnMesonString ( TString mesonName){
        if(mesonName.CompareTo("Pi0") == 0 || mesonName.CompareTo("Pi0EtaBinning") == 0 ||  mesonName.CompareTo("Pi0OmegaBinning") == 0){
            return "#pi^{0}";
        } else if(mesonName.CompareTo("Eta") == 0)  {
            return "#eta";
        } else if(mesonName.CompareTo("EtaPrime") == 0 || mesonName.CompareTo("EtaPrim") == 0)  {
            return "#eta'";
        } else if(mesonName.CompareTo("Omega") == 0)  {
            return "#omega";
        } else if(mesonName.CompareTo("CKaon") == 0)  {
            return "K^{#pm}";
        } else if(mesonName.CompareTo("CPion") == 0)  {
            return "#pi^{#pm}";
        } else if(mesonName.CompareTo("Lambda") == 0)  {
            return "#Lambda";
        } else if(mesonName.CompareTo("Rho0") == 0)  {
            return "#rho^{0}";
        } else if(mesonName.CompareTo("Phi") == 0)  {
            return "#phi";
        } else if(mesonName.CompareTo("K0Star") == 0)  {
            return "K^{0*}";
        } else if(mesonName.CompareTo("Proton") == 0)  { // not a meson, but its easier having it in this function
            return "p#overline{p}";
        } else {
            cout << "No correct meson has been selected (" << mesonName << ")" << endl;
            return "";
        }
    }

    //************************************************************************************
    //************************ return bool for respective meson name *********************
    //************************************************************************************
    Bool_t ReturnMesonOption ( TString mesonName){
        if(mesonName.CompareTo("Pi0") == 0 || mesonName.CompareTo("Pi0EtaBinning") == 0){
            return kTRUE;
        } else if(mesonName.CompareTo("Eta") == 0)  {
            return kFALSE;
        } else {
            return kFALSE;
        }
    }


    //************************************************************************************
    //***** return detailed labels for different processes depending on energy ***********
    //************************************************************************************
    TString ReturnFullTextMeson(TString fEnergyFlagOpt,
                                TString textProcessOpt){

        if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
            return Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 7 TeV ",textProcessOpt.Data());
        } else if(fEnergyFlagOpt.CompareTo("8TeV") == 0){
            return Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 8 TeV ",textProcessOpt.Data());
        } else if(fEnergyFlagOpt.CompareTo("13TeV") == 0 || fEnergyFlagOpt.CompareTo("13TeVLowB") == 0){
            return Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 13 TeV ",textProcessOpt.Data());
        } else if(fEnergyFlagOpt.CompareTo("5TeV") == 0 || fEnergyFlagOpt.Contains("5TeV2017")  || fEnergyFlagOpt.CompareTo("5TeVSpecial") == 0 ){
            return Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 5.02 TeV ",textProcessOpt.Data());
        } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
            return  Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 900 GeV ",textProcessOpt.Data());
        } else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
            return Form("pp #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 2.76 TeV ",textProcessOpt.Data());
        } else if( fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) {
            return  Form("PbPb #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 2.76 TeV ",textProcessOpt.Data());
        } else if( fEnergyFlagOpt.CompareTo("PbPb_5.02TeV") == 0) {
            return  Form("PbPb #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 5.02 TeV ",textProcessOpt.Data());
        } else if( fEnergyFlagOpt.CompareTo("XeXe_5.44TeV") == 0) {
            return  Form("XeXe #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 5.44 TeV ",textProcessOpt.Data());
        } else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVCent") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVRun2") == 0 ) {
            return  Form("pPb #rightarrow %s (#rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}) + X @ 5.02 TeV ",textProcessOpt.Data());
        } else {
            cout << "No correct collision system specification, has been given" << endl;
            return "";
        }
    }

    //************************************************************************************
    //return detailed labels Pi+pi-gamma channel for different processes depending on energy
    //************************************************************************************
    TString ReturnFullTextMesonPiPlPiMiGamma(   TString fEnergyFlagOpt,
                                                TString textProcessOpt){

        if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
            return Form("pp #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 7 TeV ",textProcessOpt.Data());
        } else if( fEnergyFlagOpt.CompareTo("8TeV") == 0) {
            return  Form("pp #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 8 TeV ",textProcessOpt.Data());
        } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
            return  Form("pp #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 900 GeV ",textProcessOpt.Data());
        } else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
            return Form("pp #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 2.76 TeV ",textProcessOpt.Data());
        } else if( fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) {
            return  Form("PbPb #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 2.76 TeV ",textProcessOpt.Data());
        } else if( fEnergyFlagOpt.CompareTo("XeXe_5.44TeV") == 0) {
            return  Form("XeXe #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 5.44 TeV ",textProcessOpt.Data());
        } else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVCent") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVRun2") == 0 ) {
            return  Form("pPb #rightarrow %s (#rightarrow #pi^{+}#pi^{-}#gamma) + X @ 5.02 TeV ",textProcessOpt.Data());
        } else {
            cout << "No correct collision system specification, has been given" << endl;
            return "";
        }
    }

    //************************************************************************************
    //***************** return proper labeling for collision system **********************
    //************************************************************************************
    TString ReturnFullCollisionsSystem( TString fEnergyFlagOpt){
      if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
            return  "pp, #sqrt{#it{s}} = 7 TeV";
        } else if( fEnergyFlagOpt.CompareTo("7TeVSys") == 0) {
            return  "pp, #sqrt{#it{s}} = 7 TeV";
        } else if( fEnergyFlagOpt.CompareTo("8TeV") == 0) {
            return  "pp, #sqrt{#it{s}} = 8 TeV";
        } else if( fEnergyFlagOpt.CompareTo("13TeV") == 0) {
            return  "pp, #sqrt{#it{s}} = 13TeV";
        } else if( fEnergyFlagOpt.CompareTo("13TeVLowB") == 0) {
            return  "pp, #sqrt{#it{s}} = 13TeV (low B)";
        } else if( fEnergyFlagOpt.CompareTo("13TeVRBins") == 0) {
            return  "pp, #sqrt{#it{s}} = 13TeV (RBins)";
        } else if( fEnergyFlagOpt.CompareTo("5TeV") == 0 || fEnergyFlagOpt.CompareTo("5.023TeV") == 0 || fEnergyFlagOpt.CompareTo("5.02TeV") == 0 ) {
            return  "pp, #sqrt{#it{s}} = 5.02 TeV";
        } else if( fEnergyFlagOpt.Contains("5TeV2017") || fEnergyFlagOpt.CompareTo("5TeVSpecial") == 0) {
            return  "pp, #sqrt{#it{s}} = 5.02 TeV";
        } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
            return  "pp, #sqrt{#it{s}} = 900 GeV";
        } else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
            return  "pp, #sqrt{#it{s}} = 2.76 TeV";
        } else if( (fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) || (fEnergyFlagOpt.CompareTo("HI") == 0) ) {
            return "Pb-Pb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
        } else if( (fEnergyFlagOpt.CompareTo("PbPb_5.02TeV") == 0) ) {
            return "Pb-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
        } else if( (fEnergyFlagOpt.CompareTo("XeXe_5.44TeV") == 0) ) {
            return "Xe-Xe, #sqrt{#it{s}_{_{NN}}} = 5.44 TeV";
        } else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVCent") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.02TeV") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVRun2") == 0) {
            return "p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
        } else if( fEnergyFlagOpt.Contains("pPb_8TeV")) {
            return "p-Pb, #sqrt{#it{s}_{_{NN}}} = 8 TeV";
        } else {
            cout << "No correct collision system specification, has been given" << endl;
            return "";
        }
    }

    //************************************************************************************
    //***************** return proper cms-energy for collision system ********************
    //************************************************************************************
    Double_t ReturnCollisionEnergy( TString fEnergyFlagOpt){
        if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
            return  7000;
        } else if( fEnergyFlagOpt.CompareTo("8TeV") == 0) {
            return 8000;
        } else if( fEnergyFlagOpt.CompareTo("13TeV") == 0 || fEnergyFlagOpt.CompareTo("13TeVLowB") == 0 || fEnergyFlagOpt.CompareTo("13TeVRBins") == 0 ) {
            return 13000;
        } else if( fEnergyFlagOpt.CompareTo("5TeV") == 0 || fEnergyFlagOpt.Contains("5TeV2017") || fEnergyFlagOpt.CompareTo("5TeVSpecial") == 0  ) {
            return 5020;
        } else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
            return 2760;
        } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
            return 900;
        } else if( (fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) || (fEnergyFlagOpt.CompareTo("HI") == 0) ) {
            return 2760;
        } else if( (fEnergyFlagOpt.CompareTo("PbPb_5.02TeV") == 0) ) {
            return 5020;
        } else if( (fEnergyFlagOpt.CompareTo("XeXe_5.44TeV") == 0) ) {
            return 5444;
        } else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVCent") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.02TeV") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVRun2") == 0 )  {
            return 5023;
        } else if( fEnergyFlagOpt.Contains("pPb_8TeV") )  {
            return 8160;
        } else {
            cout << "No correct collision system energy specification, has been given" << endl;
            return 1;
        }
    }

    //************************************************************************************
    //***************** return proper cms-energy for collision system ********************
    //************************************************************************************
    TString ReturnCollisionEnergyOutputString( TString fEnergyFlagOpt){
        if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
            return  "pp7TeV";
        } else if( fEnergyFlagOpt.CompareTo("5TeV") == 0 || fEnergyFlagOpt.CompareTo("5.023TeV") == 0|| fEnergyFlagOpt.CompareTo("5.02TeV") == 0 || fEnergyFlagOpt.Contains("5TeV2017") || fEnergyFlagOpt.CompareTo("5TeVSpecial") == 0 ) {
            return  "pp5020GeV";
        } else if( fEnergyFlagOpt.CompareTo("8TeV") == 0) {
            return  "pp8TeV";
        } else if( fEnergyFlagOpt.CompareTo("13TeV") == 0 || fEnergyFlagOpt.CompareTo("13TeVLowB") == 0 || fEnergyFlagOpt.CompareTo("13TeVRBins") == 0  ) {
            return  "pp13TeV";
        } else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
            return  "pp2760GeV";
        } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
            return  "pp900GeV";
        } else if( (fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) ) {
            return  "PbPb2760GeV";
        } else if( (fEnergyFlagOpt.CompareTo("PbPb_5.02TeV") == 0) ) {
            return  "PbPb5TeV";
        } else if( (fEnergyFlagOpt.CompareTo("XeXe_5.44TeV") == 0) ) {
            return  "XeXe5440GeV";
        } else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVCent") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.02TeV") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVRun2") == 0 ) {
            return  "pPb5TeV";
        } else if( fEnergyFlagOpt.Contains("pPb_8TeV") ) {
            return  "pPb8TeV";
        } else {
            cout << "No correct collision system energy specification, has been given" << endl;
            return "";
        }
    }

    //************************************************************************************
    //***************** return proper cms-energy for collision system ********************
    //************************************************************************************
    TString ReturnCollisionEnergyStringForTheory( TString fEnergyFlagOpt){
        if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
            return  "7TeV";
        } else if( fEnergyFlagOpt.CompareTo("5TeV") == 0 || fEnergyFlagOpt.Contains("5TeV2017") || fEnergyFlagOpt.CompareTo("5TeVSpecial") == 0 ) {
            return  "5TeV";
        } else if( fEnergyFlagOpt.CompareTo("8TeV") == 0) {
            return  "8TeV";
        } else if( fEnergyFlagOpt.CompareTo("13TeV") == 0 || fEnergyFlagOpt.CompareTo("13TeVLowB") == 0 || fEnergyFlagOpt.CompareTo("13TeVRBins") == 0 ) {
            return  "13TeV";
        } else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
            return  "2760GeV";
        } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
            return  "900GeV";
        } else if( (fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) ) {
            return  "2760GeV";
        } else if( (fEnergyFlagOpt.CompareTo("PbPb_5.02TeV") == 0) ) {
            return  "5TeV";
        } else if( (fEnergyFlagOpt.CompareTo("XeXe_5.44TeV") == 0) ) {
            return  "5440GeV";
        } else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVCent") == 0 ||  fEnergyFlagOpt.CompareTo("pPb_5.023TeVRun2") == 0) {
            return  "5TeV";
        } else if( fEnergyFlagOpt.Contains("pPb_8TeV")) {
            return  "pPb8TeV";
        } else {
            cout << "No correct collision system energy specification, has been given" << endl;
            return "";
        }
    }

    //************************************************************************************
    //********************* EventCuts definition *****************************************
    //************************************************************************************
    Int_t GetEventSystemCutPosition ()                      {return 0;}
    Int_t GetEventCentralityMinCutPosition ()               {return 1;}
    Int_t GetEventCentralityMaxCutPosition ()               {return 2;}
    Int_t GetEventSelectSpecialTriggerCutPosition ()        {return 3;}
    Int_t GetEventSelectSpecialTriggerCutPositionHeavy ()   {return 5;}
    Int_t GetEventSelectSpecialSubTriggerCutPosition ()     {return 4;}
    Int_t GetEventRemovePileUpCutPosition ()                {return 5;}
    Int_t GetEventRejectExtraSignalsCutPosition ()          {return 6;}
    Int_t GetEventVertexCutPosition ()                      {return 7;}

    //************************************************************************************
    //************************** Conversion photonCut definition *************************
    //************************************************************************************
    Int_t GetPhotonV0FinderCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() > 23) return 0;
        else return 0;
    }
    Int_t GetPhotonEtaCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() > 23) return 1;
        else return 1;
    }
    Int_t GetPhotonMinRCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() > 23) return 2;
        else return 2;
    }
    Int_t GetPhotonEtaForPhiCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() > 23) return 3;
        else return -1;
    }
    Int_t GetPhotonMinPhiCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() > 23) return 4;
        else return -1;
    }
    Int_t GetPhotonMaxPhiCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() > 23) return 5;
        else return -1;
    }
    Int_t GetPhotonSinglePtCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() > 23) return 6;
        else return 3;
    }
    Int_t GetPhotonClsTPCCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() > 23) return 7;
        else return 4;
    }
    Int_t GetPhotonEDedxSigmaCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() > 23) return 8;
        else return 5;
    }
    Int_t GetPhotonPiDedxSigmaCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() > 23) return 9;
        else return 6;
    }
    Int_t GetPhotonPiMomDedxSigmaCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() > 23) return 10;
        else return 7;
    }
    Int_t GetPhotonPiMaxMomDedxSigmaCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() > 23) return 11;
        else return 8;
    }
    Int_t GetPhotonLowPRejectionSigmaCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() > 23) return 12;
        else return 9;
    }
    Int_t GetPhotonTOFelectronPIDCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() > 23) return 13;
        else return 10;
    }
    Int_t GetPhotonQtMaxCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() == 26) return 16;
        else if(gammaCutNumber.Length() == 24) return 14;
        else return 11;
    }
    Int_t GetPhotonChi2GammaCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() == 26) return 17;
        else if(gammaCutNumber.Length() == 24) return 15;
        else return 12;
    }
    Int_t GetPhotonPsiPairCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() == 26) return 18;
        else if(gammaCutNumber.Length() == 24) return 16;
        else return 13;
    }
    Int_t GetPhotonDoPhotonAsymmetryCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() == 26) return 19;
        else if(gammaCutNumber.Length() == 24) return 17;
        else return 14;
    }
    Int_t GetPhotonCosinePointingAngleCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() == 26) return 20;
        else if(gammaCutNumber.Length() == 24) return 18;
        else return 15;
    }
    Int_t GetPhotonSharedElectronCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() == 26) return 21;
        else if(gammaCutNumber.Length() == 24) return 19;
        else return 16;
    }
    Int_t GetPhotonRejectToCloseV0sCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() == 26) return 22;
        else if(gammaCutNumber.Length() == 24) return 20;
        else return 17;
    }
    Int_t GetPhotonDcaRPrimVtxCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() == 26) return 23;
        else if(gammaCutNumber.Length() == 24) return 21;
        else return 18;
    }
    Int_t GetPhotonDcaZPrimVtxCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() == 26) return 24;
        else if(gammaCutNumber.Length() == 24) return 22;
        else return 19;
    }
    Int_t GetPhotonEventPlaneCutPosition (TString gammaCutNumber){
        if(gammaCutNumber.Length() == 26) return 25;
        else if(gammaCutNumber.Length() == 24) return 23;
        else return 20;
    }

    //************************************************************************************
    //************************** ClusterCuts definition **********************************
    //************************************************************************************
    Int_t GetClusterTypeCutPosition ( TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 0;
        else return 0;
    }
    Int_t GetClusterEtaMinCutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 1;
        else return 1;
    }
    Int_t GetClusterEtaMaxCutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 2;
        else return 2;
    }
    Int_t GetClusterPhiMinCutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 3;
        else return 3;
    }
    Int_t GetClusterPhiMaxCutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 4;
        else return 4;
    }
    Int_t GetClusterNonLinearityCutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return -1;
        else return 5;
    }
    Int_t GetClusterDistanceToBadChannelCutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 5;
        else return 7;
    }
    Int_t GetClusterTimingCutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 6;
        else return 8;
    }
    Int_t GetClusterTrackMatchingCutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 7;
        else return 9;
    }
    Int_t GetClusterExoticCellCutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 8;
        else return 10;
    }
    Int_t GetClusterMinEnergyCutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 9;
        else return 11;
    }
    Int_t GetClusterMinNCellsCutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 10;
        else return 12;
    }
    Int_t GetClusterMinM02CutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 11;
        else return 13;
    }
    Int_t GetClusterMaxM02CutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 12;
        else return 14;
    }
    Int_t GetClusterMinM20CutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 13;
        else return 15;
    }
    Int_t GetClusterMaxM20CutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 14;
        else return 16;
    }
    Int_t GetClusterMaximumDispersionCutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 15;
        else return 17;
    }
    Int_t GetClusterNLMCutPosition (TString clusterCutNumber){
        if (clusterCutNumber.Length() == 17) 	return 16;
        else return 18;
    }

    //************************************************************************************
    //***************************** MesonCuts defintion **********************************
    //************************************************************************************
    Int_t GetMesonKindCutPosition ()                        {return 0;}
    Int_t GetMesonBGSchemeCutPosition ()                    {return 1;}
    Int_t GetMesonNumberOfBGEventsCutPosition ()            {return 2;}
    Int_t GetMesonDegreesForRotationMethodCutPosition ()    {return 3;}
    Int_t GetMesonRapidityCutPosition ()                    {return 4;}
    Int_t GetMesonPtCutPosition ()                           {return 5;}
    Int_t GetMesonAlphaCutPosition ()                       {return 6;}
    Int_t GetMesonSelectionWindowCutPosition ()             {return 7;}
    Int_t GetMesonSharedElectronCutPosition ()              {return 8;}
    Int_t GetMesonRejectToCloseV0sCutPosition ()            {return 9;}
    Int_t GetMesonUseMCPSmearingCutPosition ()              {return 10;}
    Int_t GetMesonDcaGammaGammaCutPosition ()               {return 11;}
    Int_t GetMesonDcaRCutPosition ()                        {return 12;}
    Int_t GetMesonDcaZCutPosition ()                        {return 13;}
    Int_t GetMesonOpeningAngleCutPosition ()                {return 14;}

    //************************************************************************************
    //***************************** Charged pion defintion **********************************
    //************************************************************************************
    Int_t GetPionEtaCut ()                        {return 0;}
    Int_t GetPionClsITSCut ()                     {return 1;}
    Int_t GetPionClsTPCCut ()                     {return 2;}
    Int_t GetPionDCACut ()                        {return 3;}
    Int_t GetPionPtCut ()                         {return 4;}
    Int_t GetPiondEdxSigmaITSCut ()               {return 5;}
    Int_t GetPiondEdxSigmaTPCCut ()               {return 6;}
    Int_t GetPiondEdxSigmaTOFCut ()               {return 7;}
    Int_t GetPionMassCut ()                       {return 8;}

    //************************************************************************************
    //*********************** date generation for plot labeling **************************
    //************************************************************************************
    TString ReturnDateString(Bool_t reversed = kFALSE){
        TDatime today;
        int iDate           = today.GetDate();
        int iYear           = iDate/10000;
        int iMonth          = (iDate%10000)/100;
        int iDay            = iDate%100;
        TString cMonth[12]  = {"Jan","Feb","Mar","Apr","May","Jun",
                            "Jul","Aug","Sep","Oct","Nov","Dec"};
        TString textDayth;
        if (iDay== 11){
            textDayth       = "th";
        } else if  (iDay== 12){
            textDayth       = "th";
        } else if  (iDay== 13){
            textDayth       = "th";
        } else if  (iDay%10 == 1){
            textDayth       = "st";
        } else if (iDay%10 == 2){
            textDayth       = "nd";
        } else if (iDay%10 == 3){
            textDayth       = "rd";
        } else {
            textDayth       = "th";
        }
        if(reversed) return Form("%s. %i%s, %i", cMonth[iMonth-1].Data(),iDay,textDayth.Data(), iYear);
        else return Form("%i^{%s} %s %i",iDay, textDayth.Data(),cMonth[iMonth-1].Data(), iYear);
    }

    //************************************************************************************
    //*********************** date generation for output names ***************************
    //************************************************************************************
    TString ReturnDateStringForOutput(){
        TDatime today;
        int iDate           = today.GetDate();
        int iYear           = iDate/10000;
        int iMonth          = (iDate%10000)/100;
        int iDay            = iDate%100;
        TString cMonth[12]  = {"Jan","Feb","Mar","Apr","May","Jun",
                            "Jul","Aug","Sep","Oct","Nov","Dec"};
        TString textDayth;
        if (iDay== 11){
            textDayth       = "th";
        } else if  (iDay== 12){
            textDayth       = "th";
        } else if  (iDay== 13){
            textDayth       = "th";
        } else if  (iDay%10 == 1){
            textDayth       = "st";
        } else if (iDay%10 == 2){
            textDayth       = "nd";
        } else if (iDay%10 == 3){
            textDayth       = "rd";
        } else {
            textDayth       = "th";
        }
        return Form("%i_%02d_%02d",iYear, iMonth, iDay);
    }

    //************************************************************************************
    //*********************** time generation for output names ***************************
    //************************************************************************************
    TString ReturnTimeStringForOutput(){
    TDatime 	today;
    int 		iHour = today.GetHour();
    int     iMinute = today.GetMinute();

    if(iHour < 10){
        if(iMinute < 10){
        return Form("_0%i0%i",iHour, iMinute);
        }
        else{
        return Form("_0%i%i",iHour, iMinute);
        }
    }
    else{
        if(iMinute < 10){
        return Form("_%i0%i",iHour, iMinute);
        }
        else{
        return Form("_%i%i",iHour, iMinute);
        }
    }
    }

    //************************************************************************************
    //** Analyzes meson rapidity cut, returns labeling string + double for normalization *
    //************************************************************************************
    Double_t ReturnRapidityStringAndDouble( TString cutSel,
                                            TString& rapidityRangeDummy){

        TString rapitdityCutNumberDummy     = cutSel(GetMesonRapidityCutPosition(),1);
        if (rapitdityCutNumberDummy.CompareTo("0") == 0){
            cout << "using rapidity of 0.9" << endl;
            rapidityRangeDummy              = "1.35";
            return 1.35*2;
        } else if (rapitdityCutNumberDummy.CompareTo("1") == 0){
            cout << "using rapidity of 0.8" << endl;
            rapidityRangeDummy              = "0.8";
            return 1.6;
        } else if (rapitdityCutNumberDummy.CompareTo("2") == 0){
            cout << "using rapidity of 0.7" << endl;
            rapidityRangeDummy              = "0.7";
            return 1.4;
        } else if (rapitdityCutNumberDummy.CompareTo("3") == 0){
            cout << "using rapidity of 0.6" << endl;
            rapidityRangeDummy              = "0.6";
            return 1.2;
        } else if (rapitdityCutNumberDummy.CompareTo("4") == 0){
            cout << "using rapidity of 0.5" << endl;
            rapidityRangeDummy              = "0.5";
            return 1.0;
        } else if (rapitdityCutNumberDummy.CompareTo("5") == 0){
            cout << "using rapidity of 0.85" << endl;
            rapidityRangeDummy              = "0.85";
            return 1.7;
        } else if (rapitdityCutNumberDummy.CompareTo("6") == 0){
            cout << "using rapidity of 0.75" << endl;
            rapidityRangeDummy              = "0.75";
            return 1.5;
        } else if (rapitdityCutNumberDummy.CompareTo("7") == 0){
            cout << "using rapidity of 0.3" << endl;
            rapidityRangeDummy              = "0.3";
            return 0.6;
        } else if (rapitdityCutNumberDummy.CompareTo("8") == 0){
            cout << "using rapidity of 0.35" << endl;
            rapidityRangeDummy              = "0.25";
            return 0.5;
        } else if (rapitdityCutNumberDummy.CompareTo("9") == 0){
            cout << "using rapidity of 0.4" << endl;
            rapidityRangeDummy              = "0.4";
            return 0.8;
        } else {
            cout <<  " no rapidity Range selected" << endl;
            return 1.;
        }
    }

    //************************************************************************************
    //************* Analyzes photon eta cut, returns double for normalization ************
    //************************************************************************************
    Double_t ReturnDeltaEta(TString gammaCutNumber){

        TString etaCutNumber(gammaCutNumber(GetPhotonEtaCutPosition(gammaCutNumber),1));
        if (etaCutNumber.CompareTo("0")==0){
            cout << "using eta for gammas of 0.9" << endl;
            return  1.8;
        } else if (etaCutNumber.CompareTo("1")==0){
            cout << "using eta for gammas of 0.6" << endl;
            return  1.2;
        } else if (etaCutNumber.CompareTo("2")==0){
            cout << "using eta for gammas of 1.4" << endl;
            return  2.8;
        } else if (etaCutNumber.CompareTo("3")==0){
            cout << "using eta for gammas of 0.65" << endl;
            return 0.65*2;
        } else if (etaCutNumber.CompareTo("4")==0){
            cout << "using eta for gammas of 0.75" << endl;
            return  1.5;
        } else if (etaCutNumber.CompareTo("5")==0){
            cout << "using eta for gammas of 0.5" << endl;
            return  1.0;
        } else if (etaCutNumber.CompareTo("6")==0){
            cout << "using eta for gammas of 5.0" << endl;
            return  10.0;
        } else if (etaCutNumber.CompareTo("7")==0){
            cout << "using eta for gammas of 0.3" << endl;
            return  0.6;
        } else if (etaCutNumber.CompareTo("8")==0){
            cout << "using eta for gammas of 0.4" << endl;
            return 0.8;
        } else if (etaCutNumber.CompareTo("9")==0){
            cout << "using eta for gammas of 10" << endl;
            return 20.;
        } else if (etaCutNumber.CompareTo("a")==0){
            cout << "using eta for gammas of 0.2-0.9" << endl;
            return 1.4;
        } else if (etaCutNumber.CompareTo("c")==0){
            cout << "using eta for gammas of 0.85" << endl;
            return 1.7;
        } else if (etaCutNumber.CompareTo("d")==0){
            cout << "using eta for gammas of 0.8" << endl;
            return 1.6;
        }
        cout << "Eta Value NOT found!!! using eta for gammas of 0.9" << endl;
        return 1.8;
    }


    //************************************************************************************
    //******************* Analyse minimum eta cut for clusters ***************************
    //************************************************************************************
    Double_t AnalyseClusterMinEtaCut (Int_t etaMin){
        switch (etaMin){
            case 0:
                return -10.;
            case 1:
                return -0.6687;
            case 2:
                return -0.5;
            case 3:
                return -2;
            case 4:
                return -0.13;
            case 5:
                return -0.7;
            case 6:
                return -0.3;
            case 7:
                return -0.4;
            case 8:
                return -0.66112;
            case 9:
                return -0.6687;
            default:
                return 0;
        }
    }

    //************************************************************************************
    //******************* Analyse maximum eta cut for clusters ***************************
    //************************************************************************************
    Double_t AnalyseClusterMaxEtaCut (Int_t etaMin){
        switch (etaMin){
            case 0:
                return 10.;
            case 1:
                return 0.66465;
            case 2:
                return 0.5;
            case 3:
                return 2;
            case 4:
                return 0.13;
            case 5:
                return 0.7;
            case 6:
                return 0.3;
            case 7:
                return 0.4;
            case 8:
                return 0.66112;
            case 9:
                return 0.66465;
            default:
                return 0;
        }
    }

    //************************************************************************************
    //******************* Analyse minimum phi cut for clusters ***************************
    //************************************************************************************
    Double_t AnalyseClusterMinPhiCut (Int_t etaMin){
        switch (etaMin){
            case 0:
                return -10000;
            case 1:
                return 1.39626;
            case 2:
                return 2.1;
            case 3:
                return 2.45;
            case 4:
                return 4.54;
            case 5:
                return 4.5572;
            case 6:
                return 4.36;
            case 7:
                return 1.39626;
            default:
                return 0;
        }
    }

    //************************************************************************************
    //******************* Analyse maximum phi cut for clusters ***************************
    //************************************************************************************
    Double_t AnalyseClusterMaxPhiCut (Int_t etaMin){
        switch (etaMin){
            case 0:
                return 10000;
            case 1:
                return 3.15;
            case 2:
                return 2.45;
            case 3:
                return 2.10;
            case 4:
                return 5.59;
            case 5:
                return 5.5658;
            case 6:
                return 5.59;
            case 7:
                return 5.7;
            case 9:
                return 5.7;
            default:
                return 0;
        }
    }

    //************************************************************************************
    //****** Analyzes photon (cluster) eta cut, returns double for normalization *********
    //************************************************************************************
    TString AnalyseEtaCalo(TString caloCutNumber){

        TString etaMinCutNumber(caloCutNumber(GetClusterEtaMinCutPosition(caloCutNumber),1));
        TString etaMaxCutNumber(caloCutNumber(GetClusterEtaMaxCutPosition(caloCutNumber),1));

        Float_t minEtaCut   = AnalyseClusterMinEtaCut(CutNumberToInteger(etaMinCutNumber));
        Float_t maxEtaCut   = AnalyseClusterMaxEtaCut(CutNumberToInteger(etaMaxCutNumber));


        return Form("%.2f < #gamma_{calo} < %.2f", minEtaCut, maxEtaCut);
    }


    //************************************************************************************
    //****** Analyzes photon (cluster) eta cut, returns double for normalization *********
    //************************************************************************************
    Double_t ReturnDeltaEtaCalo(TString caloCutNumber, Int_t mode = -1){

        TString etaMinCutNumber(caloCutNumber(GetClusterEtaMinCutPosition(caloCutNumber),1));
        TString etaMaxCutNumber(caloCutNumber(GetClusterEtaMaxCutPosition(caloCutNumber),1));

        Float_t minEtaCut   = AnalyseClusterMinEtaCut(CutNumberToInteger(etaMinCutNumber));
        Float_t maxEtaCut   = AnalyseClusterMaxEtaCut(CutNumberToInteger(etaMaxCutNumber));
        Float_t deltaEtaCut = TMath::Abs(minEtaCut) + TMath::Abs(maxEtaCut);
        if(mode == 12) deltaEtaCut = TMath::Abs(minEtaCut) + TMath::Abs(maxEtaCut) - 2*0.227579;

        return deltaEtaCut;
    }

    //************************************************************************************
    //****** Analyzes photon (cluster) phi cut, returns double for normalization *********
    //************************************************************************************
    Double_t ReturnDeltaPhiCalo(TString caloCutNumber, Int_t mode = -1){

        TString phiMinCutNumber(caloCutNumber(GetClusterPhiMinCutPosition(caloCutNumber),1));
        TString phiMaxCutNumber(caloCutNumber(GetClusterPhiMaxCutPosition(caloCutNumber),1));

        Float_t minPhiCut   = AnalyseClusterMinPhiCut(CutNumberToInteger(phiMinCutNumber));
        Float_t maxPhiCut   = AnalyseClusterMaxPhiCut(CutNumberToInteger(phiMaxCutNumber));
        Float_t deltaPhiCut = maxPhiCut - minPhiCut;
        if ( minPhiCut == -10000 && maxPhiCut == 10000 )
            deltaPhiCut     = 2*TMath::Pi();

        return deltaPhiCut;
    }


    //************************************************************************************
    //************* Analyzes meson BG mult cut, returns float for expected BG number *****
    //************************************************************************************
    Float_t ReturnBackgroundMult(TString cutSel){
        TString fBackgroundMultCutNumberDummy = cutSel(GetMesonNumberOfBGEventsCutPosition(),1);
        if (fBackgroundMultCutNumberDummy.CompareTo("0") == 0){
            cout << "using number of events for BG 5" << endl;
            return 5;
        } else if (fBackgroundMultCutNumberDummy.CompareTo("1") == 0){
            cout << "using number of events for BG 10" << endl;
            return 10;
        } else if (fBackgroundMultCutNumberDummy.CompareTo("2") == 0){
            cout << "using number of events for BG 15" << endl;
            return 15;
        } else if (fBackgroundMultCutNumberDummy.CompareTo("3") == 0){
            cout << "using number of events for BG 20" << endl;
            return 20;
        } else if (fBackgroundMultCutNumberDummy.CompareTo("4") == 0){
            cout << "using number of events for BG 2" << endl;
            return 2;
        } else if (fBackgroundMultCutNumberDummy.CompareTo("5") == 0){
            cout << "using number of events for BG 50" << endl;
            return 50;
        } else if (fBackgroundMultCutNumberDummy.CompareTo("6") == 0){
            cout << "using number of events for BG 80" << endl;
            return 80;
        } else if (fBackgroundMultCutNumberDummy.CompareTo("7") == 0){
            cout << "using number of events for BG 100" << endl;
            return 100;
        }
        return 0;
    }

    //************************************************************************************
    //****** Analyzes the eventCutnumber and returns the NColl for the respective cent ***
    //************************************************************************************
    Double_t GetNCollFromCutNumber (    TString cutNumber,
                                        TString energy      = "PbPb_2.76TeV"
                                ){
        TString systemCutNumber     = cutNumber(GetEventSystemCutPosition(),1);
        TString centralityCutNumber = cutNumber(GetEventCentralityMinCutPosition(),2);
        TString toTest              = "";
        if (energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0 || energy.CompareTo("pPb_5.023TeVRun2") == 0 || energy.CompareTo("pPb_5.02TeV") == 0){
            return ncollpPb5023GeV;
        } else if (energy.CompareTo("PbPb_2.76TeV") == 0){
            if (systemCutNumber.CompareTo("1") == 0 || systemCutNumber.CompareTo("2") == 0 || systemCutNumber.CompareTo("5") == 0){
                // find correct NColl err for 10% slices
                for (Int_t i = 0; i < 10; i++){
                    toTest          = Form("%d%d",i,i+1);
                    if (i == 9)
                        toTest          = Form("%d%d",i,0);
                    if (centralityCutNumber.CompareTo(toTest.Data()) == 0)
                        return nCollPbPb2760GeVV0M10[i];
                }
                // find correct NColl err for 20% slices
                for (Int_t i = 0; i < 5; i++){
                    toTest          = Form("%d%d",i*2,(i+1)*2);
                    if (i == 5)
                        toTest          = Form("%d%d",i*2,0);
                    if (centralityCutNumber.CompareTo(toTest.Data()) == 0)
                        return nCollPbPb2760GeVV0M20[i];
                }

                if (centralityCutNumber.CompareTo("13") == 0){ //10-30%
                    return nCollPbPb2760GeV1030;
                } else if (centralityCutNumber.CompareTo("25") == 0){ //20-50%
                    return nCollPbPb2760GeV2050;
                } else if (centralityCutNumber.CompareTo("48") == 0){ //40-80%
                    return nCollPbPb2760GeV4080;
                } else if (centralityCutNumber.CompareTo("30") == 0){ //10-30%
                    return nCollPbPb2760GeV30100;
                } else if (centralityCutNumber.CompareTo("04") == 0){ //0-40%
                    return nCollPbPb2760GeV0040;
                }
            } else if (systemCutNumber.CompareTo("3") == 0 || systemCutNumber.CompareTo("6") == 0){
                // find correct NColl err for 5% slices
                for (Int_t i = 0; i < 10; i++){
                    toTest          = Form("%d%d",i,i+1);
                    if (centralityCutNumber.CompareTo(toTest.Data()) == 0)
                        return nCollPbPb2760GeVV0M5[i];
                }
            } else if (systemCutNumber.CompareTo("4") == 0 || systemCutNumber.CompareTo("7") == 0){
                if (centralityCutNumber.CompareTo("69") == 0){ //75-90%
                    return nCollPbPb2760GeV7590;
                }
            }
        } else if (energy.CompareTo("XeXe_5.44TeV") == 0){
            return 1.;
        } else if (energy.CompareTo("PbPb_5.02TeV") == 0){
            if (systemCutNumber.CompareTo("1") == 0 || systemCutNumber.CompareTo("5") == 0){
                // find correct NColl for 10% slices
                for (Int_t i = 0; i < 10; i++){
                    toTest          = Form("%d%d",i,i+1);
                    if (i == 9)
                        toTest          = Form("%d%d",i,0);
                    if (centralityCutNumber.CompareTo(toTest.Data()) == 0)
                        return nCollPbPb5TeVV0M10[i];
                }
                // find correct NColl for 20% slices
                for (Int_t i = 0; i < 5; i++){
                    toTest          = Form("%d%d",i*2,(i+1)*2);
                    if (i == 5)
                        toTest          = Form("%d%d",i*2,0);
                    if (centralityCutNumber.CompareTo(toTest.Data()) == 0)
                      return nCollPbPb5TeVV0M20[i];
                }
                for (Int_t i = 0; i < 4; i++){
                  toTest          = Form("%d%d",i*2+1,(i+1)*2+1);
                  if (centralityCutNumber.CompareTo(toTest.Data()) == 0)
                    return nCollPbPb5TeVV0M20_2[i];
                }
                if  (centralityCutNumber.CompareTo("59") == 0)
                  return nCollPbPb5TeVV0M5090;

                cout << "ERROR: NColl values not implemented!" << endl;
            } else if (systemCutNumber.CompareTo("3") == 0 || systemCutNumber.CompareTo("6") == 0){
                // find correct NColl err for 5% slices
                for (Int_t i = 0; i < 10; i++){
                    toTest          = Form("%d%d",i,i+1);
                    if (centralityCutNumber.CompareTo(toTest.Data()) == 0)
                        return nCollPbPb5TeVV0M5[i];
                }

            } else
                cout << "ERROR: NColl values not implemented!" << endl;
        } else {
            return 1.;
        }
        return 1.;
    }

    //************************************************************************************
    //** Analyzes the eventCutnumber and returns the NColl error for the respective cent *
    //************************************************************************************
    Double_t GetNCollErrFromCutNumber ( TString cutNumber,
                                        TString energy      = "PbPb_2.76TeV"
                                    ){

        TString systemCutNumber     = cutNumber(GetEventSystemCutPosition(),1);
        TString centralityCutNumber = cutNumber(GetEventCentralityMinCutPosition(),2);
        TString toTest              = "";
        if (energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0 || energy.CompareTo("pPb_5.023TeVRun2") == 0 || energy.CompareTo("pPb_5.02TeV") == 0){
            return ncollErrpPb5023GeV;
        } else if (energy.CompareTo("PbPb_2.76TeV") == 0){
            if (systemCutNumber.CompareTo("1") == 0 || systemCutNumber.CompareTo("2") == 0 || systemCutNumber.CompareTo("5") == 0){
                // find correct NColl err for 10% slices
                for (Int_t i = 0; i < 10; i++){
                    toTest          = Form("%d%d",i,i+1);
                    if (i == 9)
                        toTest          = Form("%d%d",i,0);
                    if (centralityCutNumber.CompareTo(toTest.Data()) == 0)
                        return nCollPbPb2760GeVErrV0M10[i];
                }
                // find correct NColl err for 20% slices
                for (Int_t i = 0; i < 5; i++){
                    toTest          = Form("%d%d",i*2,(i+1)*2);
                    if (i == 5)
                        toTest          = Form("%d%d",i*2,0);
                    if (centralityCutNumber.CompareTo(toTest.Data()) == 0)
                        return nCollPbPb2760GeVErrV0M20[i];
                }
                if (centralityCutNumber.CompareTo("25") == 0){ //20-40%
                    return nCollPbPb2760GeVErr2050;
                } else if (centralityCutNumber.CompareTo("48") == 0){ //40-80%
                    return nCollPbPb2760GeVErr4080;
                }
            } else if (systemCutNumber.CompareTo("3") == 0 || systemCutNumber.CompareTo("6") == 0){
                for (Int_t i = 0; i < 10; i++){
                    if (centralityCutNumber.CompareTo(Form("%d%d",i,i+1)) == 0)
                        return nCollPbPb2760GeVV0M5[i];
                }
            } else if (systemCutNumber.CompareTo("4") == 0 || systemCutNumber.CompareTo("7") == 0){
                if (centralityCutNumber.CompareTo("69") == 0){ //75-90%
                    return nCollPbPb2760GeVErr7590;
                }
            }
            cout << "ERROR: NColl Err values not implemented!" << endl;
        } else if (energy.CompareTo("XeXe_5.44TeV") == 0){
            return 1.;
        } else if (energy.CompareTo("PbPb_5.02TeV") == 0){
            if (systemCutNumber.CompareTo("1") == 0 || systemCutNumber.CompareTo("5") == 0 ){
                // find correct NColl err for 10% slices
                for (Int_t i = 0; i < 10; i++){
                    toTest          = Form("%d%d",i,i+1);
                    if (i == 9)
                        toTest          = Form("%d%d",i,0);
                    if (centralityCutNumber.CompareTo(toTest.Data()) == 0)
                        return nCollPbPb5TeVErrV0M10[i];
                }
                // find correct NColl err for 20% slices
                for (Int_t i = 0; i < 5; i++){
                    toTest          = Form("%d%d",i*2,(i+1)*2);
                    if (i == 5)
                        toTest          = Form("%d%d",i*2,0);
                    if (centralityCutNumber.CompareTo(toTest.Data()) == 0)
                        return nCollPbPb5TeVErrV0M20[i];
                }
                for (Int_t i = 0; i < 4; i++){
                  toTest          = Form("%d%d",i*2+1,(i+1)*2+1);
                  if (centralityCutNumber.CompareTo(toTest.Data()) == 0)
                    return nCollPbPb5TeVErrV0M20_2[i];
                }
                if  (centralityCutNumber.CompareTo("59") == 0)
                  return nCollPbPb5TeVErrV0M5090;
                cout << "ERROR: NColl Err values not implemented!" << endl;

            } else if (systemCutNumber.CompareTo("3") == 0 || systemCutNumber.CompareTo("6") == 0){
                // find correct NColl err for 5% slices
                for (Int_t i = 0; i < 10; i++){
                    toTest          = Form("%d%d",i,i+1);
                    if (centralityCutNumber.CompareTo(toTest.Data()) == 0)
                        return nCollPbPb5TeVErrV0M5[i];
                }
                cout << "ERROR: NColl Err values not implemented!" << endl;
            } else
                cout << "ERROR: NColl Err values not implemented!" << endl;
        } else {
            return 1.;
        }
        return 1;
    }


    //************************************************************************************
    //***** Analyzes the name of the cent and return the NColl for the respective cent ***
    //************************************************************************************
    Double_t GetNCollFromName ( TString name,
                                TString energy  = "PbPb_2.76TeV"
                            ){
        TString toTest          = "";
        if (energy.CompareTo("pAu_0.2TeV") == 0){
            if (name.CompareTo("0005") == 0){ //0-5%
                return 9.59;
            } else if (name.CompareTo("00100") == 0){ //0-100%
                return 4.667;
            }
        } else if (energy.CompareTo("dAu_0.2TeV") == 0){
            if (name.CompareTo("00100") == 0) //0-100%
                return 7.59;
        } else if ( energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0 ||
                    energy.CompareTo("pPb_5.023TeVRun2") == 0 ||  energy.CompareTo("pPb_5.02TeV") == 0 || energy.CompareTo("pPb_5TeV") == 0){
            for (Int_t i = 0; i < 2; i++){
                toTest          = Form("%02d%02d_V0A",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPb5TeVBaseV0A5[i];
                toTest          = Form("%02d%02d_CL1",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPb5TeVBaseCL15[i];
                toTest          = Form("%02d%02d_ZNA",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPb5TeVBaseHybrid5[i];
            }
            for (Int_t i = 0; i < 5; i++){
                toTest          = Form("%02d%d_V0A",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPb5TeVBaseV0A20[i];
                toTest          = Form("%02d%d_CL1",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPb5TeVBaseCL120[i];
                toTest          = Form("%02d%d_ZNA",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPb5TeVBaseHybrid20[i];
            }
            for (Int_t i = 0; i < 10; i++){
                toTest          = Form("%02d%d_V0A",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPb5TeVBaseV0A10[i];
                toTest          = Form("%02d%d_CL1",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPb5TeVBaseCL110[i];
                toTest          = Form("%02d%d_ZNA",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPb5TeVBaseHybrid10[i];
            }

            if (name.CompareTo("60100_V0A") == 0){ //60-100%
                return nCollpPb5TeVV0A60100;
            } else if (name.CompareTo("60100_CL1") == 0){ //60-100%
                return nCollpPb5TeVCL160100;
            } else if (name.CompareTo("60100_ZNA") == 0){ //60-100%
                return nCollpPb5TeVHybrid60100;
            } else if (name.CompareTo("00100_CL1") == 0 || name.CompareTo("00100_V0A") == 0 || name.CompareTo("00100_ZNA") == 0){ //0-100%
                return ncollpPb5023GeV;
            } else {
                return ncollpPb5023GeV;
            }
        } else if (energy.CompareTo("CuCu_0.2TeV") == 0){
            //Centrality dNch/dη Ncoll Npart
            // 0%–40%   109.3±7.8   108.2±12.0  66.4±2.5
            // MB       51.7±3.6    51.8±5.6    34.6±1.2
            if (name.CompareTo("0040") == 0){ //0-40%
                return 108.2;
            } else if (name.CompareTo("0094") == 0){ //0-94%
                return 51.8;
            }
        } else if (energy.CompareTo("AuAu_0.2TeV") == 0){
            double dNdeta200AuAu[4]     = {5.1902e+02,  2.2543e+02,  8.5475e+01,  1.6362e+01};
            double edNdeta200AuAu[4]    = {2.6250e+01,  1.3175e+01,  8.0750e+00,  2.8137e+00};
//             0%–20% 770.6 ± 79.9 277.5 ± 6.5 735.2 ± 14.6 239 ± 25 ± 7
//             20%–40% 282.4 ± 28.4 135.6 ± 7.0 333.2 ± 10.7 260 ± 33 ± 8
//             40%–60% 82.6 ± 9.3 56.0 ± 5.3 126.6 ± 6.1 225 ± 28 ± 6
//             60%–92% 12.1 ± 3.1 12.5 ± 2.6 25.8 ± 4.0 238 ± 50 ± 6
//             0%–92% 251.1 ± 26.7 106.3 ± 5.0 268.8 ± 8.2 242 ± 28 ± 7
            if (name.CompareTo("0020") == 0){ //0-20%
                return 770.6;
            } else if (name.CompareTo("2040") == 0){ //20-40%
                return 282.4;
            } else if (name.CompareTo("4060") == 0){ //40-60%
                return 82.6;
            } else if (name.CompareTo("6092") == 0){ //40-60%
                return 12.1;
            } else if (name.CompareTo("0092") == 0){ //40-60%
                return 251.1;
            }
        } else if (energy.CompareTo("AuAu_0.2TeV_STAR") == 0){
            // 0-20% (766 ± 28)/42 mb
            // 20-40% (291 ± 30)/42 mb
            // 0-80% (292 ± 20)/42 mb
            // 40-60% (91 ± 20)/42 mb
            // 60-80% (22 ± 8)/42 mb
            if (name.CompareTo("0020") == 0){ //0-20%
                return 766;
            } else if (name.CompareTo("2040") == 0){ //20-40%
                return 291;
            } else if (name.CompareTo("4060") == 0){ //40-60%
                return 91;
            } else if (name.CompareTo("6080") == 0){ //60-80%
                return 22;
            } else if (name.CompareTo("0080") == 0){ //0-80%
                return 292;
            }
        } else if (energy.CompareTo("AuAu_62.4GeV") == 0){
            double dNdeta62AuAu[3]  = {3.412e+02, 1.518e+02, 1.315e+02};
            double edNdeta62AuAu[3] = {2.932e+01, 1.269e+01, 1.115e+01};
            if (name.CompareTo("0020") == 0){ //0-40%
                return 653;
            } else if (name.CompareTo("2040") == 0){ //0-94%
                return 240;
            } else if (name.CompareTo("0086") == 0){ //0-94%
                return 200;
            }
        } else if (energy.CompareTo("AuAu_39GeV") == 0){
            double dNdeta39AuAu     = 1.043e+02;
            double edNdeta39AuAu    = 8.882e+00;
            return 192;
        } else if (energy.CompareTo("PbPb_17.2GeV") == 0){
            return 660;
        } else if (energy.CompareTo("PbPb_2.76TeV_ATLAS") == 0){
            if (name.CompareTo("0010") == 0){ //0-10%
                return 1500.6;
            } else if (name.CompareTo("1020") == 0){ //10-20%
                return 923.3;
            } else if (name.CompareTo("0010") == 0){ //20-40%
                return 440.6;
            } else if (name.CompareTo("0010") == 0){ //40-80%
                return 77.8;
            }
        } else if (energy.CompareTo("PbPb_2.76TeV") == 0){
            // find correct NColl err for 5% slices
            for (Int_t i = 0; i < 20; i++){
                toTest          = Form("%02d%02d",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollPbPb2760GeVV0M5[i];
            }
            // find correct NColl err for 10% slices
            for (Int_t i = 0; i < 10; i++){
                toTest          = Form("%02d%d",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollPbPb2760GeVV0M10[i];
            }
            // find correct NColl err for 20% slices
            for (Int_t i = 0; i < 5; i++){
                toTest          = Form("%02d%d",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollPbPb2760GeVV0M20[i];
            }
            if (name.CompareTo("1030") == 0){ //10-30%
                return nCollPbPb2760GeV1030;
            } else if (name.CompareTo("2050") == 0){ //20-50%
                return nCollPbPb2760GeV2050;
            } else if (name.CompareTo("4080") == 0){ //40-80%
                return nCollPbPb2760GeV4080;
            } else if (name.CompareTo("7590") == 0){ //60-80%
                return nCollPbPb2760GeV7590;
            } else if (name.CompareTo("30100") == 0){ //30-100%
                return nCollPbPb2760GeV30100;
            }
        } else if (energy.CompareTo("PbPb_5.02TeV") == 0){
            // find correct NColl for 5% slices
            for (Int_t i = 0; i < 20; i++){
                toTest          = Form("%02d%02d",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollPbPb5TeVV0M5[i];
            }
            // find correct NColl for 10% slices
            for (Int_t i = 0; i < 10; i++){
                toTest          = Form("%02d%d",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollPbPb5TeVV0M10[i];
            }
            // find correct NColl for 20% slices
            for (Int_t i = 0; i < 5; i++){
                toTest          = Form("%02d%d",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollPbPb5TeVV0M20[i];
            }
            for (Int_t i = 0; i < 4; i++){
              toTest          = Form("%02d%d",i*20+10,(i+1)*20+10);
              if (name.CompareTo(toTest.Data()) == 0)
                return nCollPbPb5TeVV0M20_2[i];
            }
            if  (name.CompareTo("5090") == 0)
              return nCollPbPb5TeVV0M5090;
            cout << "ERROR: NColl values not implemented!" << endl;
        } else {
            return 1.;
        }
        return 1;
    }

    //************************************************************************************
    //* Analyzes the name of the cent and return the NColl error for the respective cent *
    //************************************************************************************
    Double_t GetNCollErrFromName (  TString name,
                                    TString energy  = "PbPb_2.76TeV"
                                ){
        TString toTest      = "";
        if (energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0 ||
            energy.CompareTo("pPb_5.023TeVRun2") == 0 || energy.CompareTo("pPb_5.02TeV") == 0 || energy.CompareTo("pPb_5TeV") == 0){
            for (Int_t i = 0; i < 2; i++){
                toTest          = Form("%02d%02d_V0A",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPbErr5TeVBaseV0A5[i];
                toTest          = Form("%02d%02d_CL1",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPbErr5TeVBaseCL15[i];
                toTest          = Form("%02d%02d_ZNA",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPbErr5TeVBaseHybrid5[i];
            }
            for (Int_t i = 0; i < 5; i++){
                toTest          = Form("%02d%d_V0A",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPbErr5TeVBaseV0A20[i];
                toTest          = Form("%02d%d_CL1",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPbErr5TeVBaseCL120[i];
                toTest          = Form("%02d%d_ZNA",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPbErr5TeVBaseHybrid20[i];
            }
            for (Int_t i = 0; i < 10; i++){
                toTest          = Form("%02d%d_V0A",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPbErr5TeVBaseV0A10[i];
                toTest          = Form("%02d%d_CL1",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPbErr5TeVBaseCL110[i];
                toTest          = Form("%02d%d_ZNA",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollpPbErr5TeVBaseHybrid10[i];
            }

            if (name.CompareTo("60100_V0A") == 0){ //60-100%
                return nCollpPbErr5TeVV0A60100;
            } else if (name.CompareTo("60100_CL1") == 0){ //60-100%
                return nCollpPbErr5TeVCL160100;
            } else if (name.CompareTo("60100_ZNA") == 0){ //60-100%
                return nCollpPbErr5TeVHybrid60100;
            } else if (name.CompareTo("00100_CL1") == 0 || name.CompareTo("00100_V0A") == 0 || name.CompareTo("00100_ZNA") == 0){ //0-100%
                return ncollErrpPb5023GeV;
            } else {
                return ncollErrpPb5023GeV;
            }
        } else if (energy.CompareTo("CuCu_0.2TeV") == 0){
            //Centrality dNch/dη Ncoll Npart
            // 0%–40%   109.3±7.8   108.2±12.0  66.4±2.5
            // MB       51.7±3.6    51.8±5.6    34.6±1.2
            if (name.CompareTo("0040") == 0){ //0-40%
                return 12.0;
            } else if (name.CompareTo("0094") == 0){ //0-94%
                return 5.6;
            }

        } else if (energy.CompareTo("AuAu_0.2TeV") == 0){
            //             0%–20% 770.6 ± 79.9 277.5 ± 6.5 735.2 ± 14.6 239 ± 25 ± 7
            //             20%–40% 282.4 ± 28.4 135.6 ± 7.0 333.2 ± 10.7 260 ± 33 ± 8
            //             40%–60% 82.6 ± 9.3 56.0 ± 5.3 126.6 ± 6.1 225 ± 28 ± 6
            //             60%–92% 12.1 ± 3.1 12.5 ± 2.6 25.8 ± 4.0 238 ± 50 ± 6
            //             0%–92% 251.1 ± 26.7 106.3 ± 5.0 268.8 ± 8.2 242 ± 28 ± 7
            if (name.CompareTo("0020") == 0){ //0-20%
                return 79.9;
            } else if (name.CompareTo("2040") == 0){ //20-40%
                return 28.4;
            } else if (name.CompareTo("4060") == 0){ //40-60%
                return 9.3;
            } else if (name.CompareTo("6092") == 0){ //40-60%
                return 3.1;
            } else if (name.CompareTo("0092") == 0){ //40-60%
                return 26.7;
            }
        } else if (energy.CompareTo("XeXe_5.44TeV") == 0){
            return 1.;
        } else if (energy.CompareTo("PbPb_2.76TeV") == 0){

            // find correct NColl err for 5% slices
            for (Int_t i = 0; i < 20; i++){
                toTest          = Form("%02d%02d",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollPbPb2760GeVErrV0M5[i];
            }
            // find correct NColl err for 10% slices
            for (Int_t i = 0; i < 10; i++){
                toTest          = Form("%02d%d",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollPbPb2760GeVErrV0M10[i];
            }
            // find correct NColl err for 20% slices
            for (Int_t i = 0; i < 5; i++){
                toTest          = Form("%02d%d",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollPbPb2760GeVErrV0M20[i];
            }

            if (name.CompareTo("2050") == 0){ //20-40%
                return nCollPbPb2760GeVErr2050;
            } else if (name.CompareTo("4080") == 0){ //40-80%
                return nCollPbPb2760GeVErr4080;
            } else if (name.CompareTo("7590") == 0){ //75-90%
                return nCollPbPb2760GeVErr7590;
            }

        } else if (energy.CompareTo("PbPb_5.02TeV") == 0){
            // find correct NColl for 5% slices
            for (Int_t i = 0; i < 20; i++){
                toTest          = Form("%02d%02d",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollPbPb5TeVErrV0M5[i];
            }
            // find correct NColl err for 10% slices
            for (Int_t i = 0; i < 10; i++){
                toTest          = Form("%02d%d",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollPbPb5TeVErrV0M10[i];
            }
            // find correct NColl err for 20% slices
            for (Int_t i = 0; i < 5; i++){
                toTest          = Form("%02d%d",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return nCollPbPb5TeVErrV0M20[i];
            }
            for (Int_t i = 0; i < 4; i++){
              toTest          = Form("%02d%d",i*20+10,(i+1)*20+10);
              if (name.CompareTo(toTest.Data()) == 0)
                return nCollPbPb5TeVErrV0M20_2[i];
            }
            if  (name.CompareTo("5090") == 0)
              return nCollPbPb5TeVErrV0M5090;
            cout << "ERROR: NColl Err values not implemented!" << endl;
        } else {
            return 1.;
        }
        return 1.;
    }


    //************************************************************************************
    //******* Analyzes the eventCutnumber and returns the secondary scaling factor for ***
    //******* K0s for the respective cent                                              ***
    //************************************************************************************
    Double_t GetScalingFactorSecCorrection(TString cutNumber){

        TString systemCutNumber = cutNumber(GetEventSystemCutPosition(),1);
        TString centralityCutNumber = cutNumber(GetEventCentralityMinCutPosition(),2);
        if (systemCutNumber.CompareTo("1") == 0 || systemCutNumber.CompareTo("2") == 0 || systemCutNumber.CompareTo("5") == 0){
            if (centralityCutNumber.CompareTo("02") == 0){ //0-20%
                return 1./0.302 -1.;
            } else if (centralityCutNumber.CompareTo("01") == 0){ //0-10%
                return 1./0.303 -1.;
            } else if (centralityCutNumber.CompareTo("12") == 0){ //10-20%
                return 1./0.2989 -1.;
            } else if (centralityCutNumber.CompareTo("24") == 0){ //20-40%
                return 1./0.308 -1.;
            } else if (centralityCutNumber.CompareTo("25") == 0){ //20-50%
                return 1./0.308 -1.;
            } else if (centralityCutNumber.CompareTo("46") == 0){ //40-60%
                return 1./0.344 -1.;
            } else if (centralityCutNumber.CompareTo("45") == 0){ //40-50%
                return 1./0.342 -1.;
            } else if (centralityCutNumber.CompareTo("56") == 0){ //50-60%
                return 1./0.346 -1.;
            } else if (centralityCutNumber.CompareTo("68") == 0){ //60-80%
                return 1./0.3979 -1.;
            } else if (centralityCutNumber.CompareTo("67") == 0){ //60-70%
                return 1./0.39 -1.;
            } else if (centralityCutNumber.CompareTo("78") == 0){ //70-80%
                return 1./0.41 -1.;
            } else if (centralityCutNumber.CompareTo("89") == 0){ //80-90%
                return 1./0.42 -1.;
            } else if (centralityCutNumber.CompareTo("48") == 0){ //40-80%
                return 1./0.36 -1.;
            } else if (centralityCutNumber.CompareTo("04") == 0){ //0-40%
                return 1./0.303 -1.;
            } else if (centralityCutNumber.CompareTo("08") == 0){ //0-80%
                return 1./0.3 -1.;
            }
        } else if (systemCutNumber.CompareTo("3") == 0 || systemCutNumber.CompareTo("6") == 0){
            if (centralityCutNumber.CompareTo("01") == 0){ //0-5%
                return 1./0.306  -1.;
            } else if (centralityCutNumber.CompareTo("12") == 0){ //0-5%
                return 1./0.301  -1.;
            }
        } else if (systemCutNumber.CompareTo("4") == 0 || systemCutNumber.CompareTo("7") == 0){
            if (centralityCutNumber.CompareTo("69") == 0){ //75-90%
            return 1./0.42 -1.;
            }
        } else if (systemCutNumber.CompareTo("8") == 0 || systemCutNumber.CompareTo("9") == 0){
            if (centralityCutNumber.CompareTo("00") == 0){ //MB
                return 1./0.62 -1.; //value for HIJING MC, for DPMJet: 1./0.52 -1.;
            } else if (centralityCutNumber.CompareTo("02") == 0){ //00-20%
                return 1./0.24 -1.;
            } else if (centralityCutNumber.CompareTo("24") == 0){ //20-40%
                return 1./0.39 -1.;
            } else if (centralityCutNumber.CompareTo("46") == 0){ //40-60%
                return 1./0.61 -1.;
            } else if (centralityCutNumber.CompareTo("68") == 0){ //60-80%
                return 1./1.11 -1.;
            } else if (centralityCutNumber.CompareTo("60") == 0){ //60-100%
                return 1./1.62 -1.;
            } else if (centralityCutNumber.CompareTo("80") == 0){ //80-100%
                return 1./3.00 -1.;
            }
        } else return 0.;
        return 0.;
    }

    //************************************************************************************
    //******* Analyzes the eventCutnumber and returns the secondary scaling factor for ***
    //******* K0s for the respective collision system                                  ***
    //************************************************************************************
    Double_t ReturnCorrectK0ScalingFactor(TString fEnergyFlagOpt, TString cutNr){
        if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
            return 1./0.75 -1.;
        } else if( fEnergyFlagOpt.CompareTo("8TeV") == 0) {
            return  1./0.75 -1.;
        } else if( fEnergyFlagOpt.CompareTo("13TeV") == 0 || fEnergyFlagOpt.CompareTo("13TeVLowB") == 0 || fEnergyFlagOpt.CompareTo("13TeVRBins") == 0  ) {
            cout << "Caution: no correct K0 Scaling factor for 13TeV available yet" << endl;
            return  1./1. -1.;
         } else if( fEnergyFlagOpt.CompareTo("5TeV") == 0 || fEnergyFlagOpt.Contains("5TeV2017")  || fEnergyFlagOpt.CompareTo("5TeVSpecial") == 0 ) {
            cout << "Same as 7 TeV K0 scaling factor" << endl;
            return  1./1. - 1; //1./0.75 -1.;
        } else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
            return  1./0.685 -1.;
        } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
            return 1./0.6 -1.;
        } else if( (fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) || (fEnergyFlagOpt.CompareTo("HI") == 0) ) {
            return GetScalingFactorSecCorrection(cutNr.Data());
        } else if( (fEnergyFlagOpt.CompareTo("PbPb_5.02TeV") == 0) ) {
            return GetScalingFactorSecCorrection(cutNr.Data());
        } else if( (fEnergyFlagOpt.CompareTo("XeXe_5.44TeV") == 0) ) {
            return GetScalingFactorSecCorrection(cutNr.Data());
        } else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVCent") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVRun2") == 0) {
            // return 0.;
            return  GetScalingFactorSecCorrection(cutNr.Data());
        } else {
            cout << "No correct collision system specification, has been given" << endl;
            return 0.;
        }
    }

    //************************************************************************************
    //***** Analyzes the name of the cent and return the TAA for the respective cent *****
    //************************************************************************************
    Double_t GetTAAFromName (   TString name,
                                TString energy  = "PbPb_2.76TeV"
                            ){
        TString toTest          = "";
        if (energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0 || energy.CompareTo("pPb_5.023TeVRun2") == 0 || energy.CompareTo("pPb_5.02TeV") == 0 ){

            for (Int_t i = 0; i < 2; i++){
                toTest          = Form("%02d%02d_V0A",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPb5TeVBaseV0A5[i]*1e3*(1/recalcBarn);
                toTest          = Form("%02d%02d_CL1",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPb5TeVBaseCL15[i]*1e3*(1/recalcBarn);
                toTest          = Form("%02d%02d_ZNA",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPb5TeVBaseHybrid5[i]*1e3*(1/recalcBarn);
            }
            for (Int_t i = 0; i < 5; i++){
                toTest          = Form("%02d%d_V0A",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPb5TeVBaseV0A20[i]*1e3*(1/recalcBarn);
                toTest          = Form("%02d%d_CL1",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPb5TeVBaseCL120[i]*1e3*(1/recalcBarn);
                toTest          = Form("%02d%d_ZNA",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPb5TeVBaseHybrid20[i]*1e3*(1/recalcBarn);
            }
            for (Int_t i = 0; i < 10; i++){
                toTest          = Form("%02d%d_V0A",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPb5TeVBaseV0A10[i]*1e3*(1/recalcBarn);
                toTest          = Form("%02d%d_CL1",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPb5TeVBaseCL110[i]*1e3*(1/recalcBarn);
                toTest          = Form("%02d%d_ZNA",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPb5TeVBaseHybrid10[i]*1e3*(1/recalcBarn);
            }

            if (name.CompareTo("60100_V0A") == 0){ //60-100%
                return tpPb5TeVBaseV0A60100*1e3*(1/recalcBarn);
            } else if (name.CompareTo("60100_CL1") == 0){ //60-100%
                return tpPb5TeVBaseCL160100*1e3*(1/recalcBarn);
            } else if (name.CompareTo("60100_ZNA") == 0){ //60-100%
                return tpPb5TeVBaseHybrid60100*1e3*(1/recalcBarn);
            } else if (name.CompareTo("00100_CL1") == 0 || name.CompareTo("00100_V0A") == 0 || name.CompareTo("00100_ZNA") == 0){ //0-100%
                return tpPb5023GeV;
            } else {
                return tpPb5023GeV;
            }
        } else if (energy.CompareTo("XeXe_5.44TeV") == 0){
            return 1.;
        } else if (energy.CompareTo("PbPb_2.76TeV") == 0){
            // find correct TAA for 5% slices
            for (Int_t i = 0; i < 20; i++){
                toTest          = Form("%02d%02d",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tAAPbPb2760GeVV0M5[i];
            }
            // find correct TAA for 10% slices
            for (Int_t i = 0; i < 10; i++){
                toTest          = Form("%02d%d",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tAAPbPb2760GeVV0M10[i];
            }
            // find correct TAA for 20% slices
            for (Int_t i = 0; i < 5; i++){
                toTest          = Form("%02d%d",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tAAPbPb2760GeVV0M20[i];
            }
            if (name.CompareTo("2050") == 0)        //20-50%
                return tAAPbPb2760GeVErr2050;
            cout << "ERROR: TAA values not implemented!" << endl;
        } else if (energy.CompareTo("PbPb_5.02TeV") == 0){
            // find correct TAA for 5% slices
            for (Int_t i = 0; i < 20; i++){
                toTest          = Form("%02d%02d",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tAAPbPb5TeVV0M5[i];
            }
            // find correct TAA for 10% slices
            for (Int_t i = 0; i < 10; i++){
                if (name.CompareTo(Form("%02d%d",i*10,(i+1)*10)) == 0)
                    return tAAPbPb5TeVV0M10[i];
            }
            // find correct TAA for 20% slices
            for (Int_t i = 0; i < 5; i++){
                if (name.CompareTo(Form("%02d%d",i*20,(i+1)*20)) == 0)
                    return tAAPbPb5TeVV0M20[i];
            }
            for (Int_t i = 0; i < 4; i++){
              toTest          = Form("%02d%d",i*20+10,(i+1)*20+10);
              if (name.CompareTo(toTest.Data()) == 0)
                return tAAPbPb5TeVV0M20_2[i];
            }
            if  (name.CompareTo("5090") == 0)
              return tAAPbPb5TeVV0M5090;


            cout << "ERROR: TAA values not implemented!" << endl;
        } else {
            return 1.;
        }
        return 1.;
    }

    //************************************************************************************
    //** Analyzes the name of the cent and return the TAA error for the respective cent **
    //************************************************************************************
    Double_t GetTAAErrFromName (TString name,
                                TString energy  = "PbPb_2.76TeV"
                                ){
        TString toTest  = "";
        if (energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0 || energy.CompareTo("pPb_5.023TeVRun2") == 0 || energy.CompareTo("pPb_5.02TeV") == 0 ){
            for (Int_t i = 0; i < 2; i++){
                toTest          = Form("%02d%02d_V0A",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPbErr5TeVBaseV0A5[i]*1e3*(1/recalcBarn);
                toTest          = Form("%02d%02d_CL1",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPbErr5TeVBaseCL15[i]*1e3*(1/recalcBarn);
                toTest          = Form("%02d%02d_ZNA",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPbErr5TeVBaseHybrid5[i]*1e3*(1/recalcBarn);
            }
            for (Int_t i = 0; i < 5; i++){
                toTest          = Form("%02d%d_V0A",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPbErr5TeVBaseV0A20[i]*1e3*(1/recalcBarn);
                toTest          = Form("%02d%d_CL1",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPbErr5TeVBaseCL120[i]*1e3*(1/recalcBarn);
                toTest          = Form("%02d%d_ZNA",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPbErr5TeVBaseHybrid20[i]*1e3*(1/recalcBarn);
            }
            for (Int_t i = 0; i < 10; i++){
                toTest          = Form("%02d%d_V0A",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPbErr5TeVBaseV0A10[i]*1e3*(1/recalcBarn);
                toTest          = Form("%02d%d_CL1",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPbErr5TeVBaseCL110[i]*1e3*(1/recalcBarn);
                toTest          = Form("%02d%d_ZNA",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tpPbErr5TeVBaseHybrid10[i]*1e3*(1/recalcBarn);
            }
            if (name.CompareTo("60100_V0A") == 0){ //60-100%
                return tpPbErr5TeVBaseV0A60100*1e3*(1/recalcBarn);
            } else if (name.CompareTo("60100_CL1") == 0){ //60-100%
                return tpPbErr5TeVBaseCL160100*1e3*(1/recalcBarn);
            } else if (name.CompareTo("60100_ZNA") == 0){ //60-100%
                return tpPbErr5TeVBaseHybrid60100*1e3*(1/recalcBarn);
            } else if (name.CompareTo("00100_CL1") == 0 || name.CompareTo("00100_V0A") == 0 || name.CompareTo("00100_ZNA") == 0){ //0-100%
                return tpPbErr5023GeV;
            } else {
                return tpPbErr5023GeV;
            }
        } else if (energy.CompareTo("XeXe_5.44TeV") == 0){
            return 1.;
        } else if (energy.CompareTo("PbPb_2.76TeV") == 0){
            // find correct TAA for 5% slices
            for (Int_t i = 0; i < 20; i++){
                toTest          = Form("%02d%02d",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tAAPbPb2760GeVErrV0M5[i];
            }
            // find correct TAA for 10% slices
            for (Int_t i = 0; i < 10; i++){
                toTest          = Form("%02d%d",i*10,(i+1)*10);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tAAPbPb2760GeVErrV0M10[i];
            }
            // find correct TAA for 20% slices
            for (Int_t i = 0; i < 5; i++){
                toTest          = Form("%02d%d",i*20,(i+1)*20);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tAAPbPb2760GeVErrV0M20[i];
            }
            if (name.CompareTo("2050") == 0) //20-50%
                return tAAPbPb2760GeVErr2050;
            cout << "ERROR: TAA values not implemented!" << endl;
        } else if (energy.CompareTo("PbPb_5.02TeV") == 0){
            // find correct TAA for 5% slices
            for (Int_t i = 0; i < 20; i++){
                toTest          = Form("%02d%02d",i*5,(i+1)*5);
                if (name.CompareTo(toTest.Data()) == 0)
                    return tAAPbPb5TeVErrV0M5[i];
            }
            // find correct TAA err for 10% slices
            for (Int_t i = 0; i < 10; i++){
                if (name.CompareTo(Form("%02d%d",i*10,(i+1)*10)) == 0)
                    return tAAPbPb5TeVErrV0M10[i];
            }
            // find correct TAA err for 20% slices
            for (Int_t i = 0; i < 5; i++){
                if (name.CompareTo(Form("%02d%d",i*20,(i+1)*20)) == 0)
                    return tAAPbPb5TeVErrV0M20[i];
            }
            for (Int_t i = 0; i < 4; i++){
              toTest          = Form("%02d%d",i*20+10,(i+1)*20+10);
              if (name.CompareTo(toTest.Data()) == 0)
                return tAAPbPb5TeVErrV0M20_2[i];
            }
            if  (name.CompareTo("5090") == 0)
              return tAAPbPb5TeVErrV0M5090;
            cout << "ERROR: TAA Err values not implemented!" << endl;
        } else {
            return 1.;
        }
        return 1.;
    }

    //************************************************************************************
    //** Analyzes the eventCutNumber and returns the staring parameters for the fit ******
    //************************************************************************************
    void ReturnParameterSetFittingPbPb(TString cutsel, Double_t* parameters){
        TString centralityCutNumber = cutsel(GetEventSystemCutPosition(),3);
        cout << "bla here" << endl;
        if ( centralityCutNumber.CompareTo("102") == 0 ||
            centralityCutNumber.CompareTo("101") == 0 ||
            centralityCutNumber.CompareTo("112") == 0 ||
            centralityCutNumber.CompareTo("502") == 0 ||
            centralityCutNumber.CompareTo("501") == 0 ||
            centralityCutNumber.CompareTo("512") == 0 ){ //0-20%
            parameters[0] = 100.;
            parameters[1] = 1500.;
            parameters[2] = 4.;
            parameters[3] = 40.;
            parameters[4] = 0.07;
            parameters[5] = 0.7;
            parameters[6] = 3.;
            parameters[7] = 4.5;
            parameters[8] = 2.;
            parameters[9] = 18.;
        } else if (centralityCutNumber.CompareTo("301") == 0 ||
                centralityCutNumber.CompareTo("601") == 0){ //0-5%
            parameters[0] = 50.;
            parameters[1] = 1500.;
            parameters[2] = 4.;
            parameters[3] = 40.;
            parameters[4] = 0.07;
            parameters[5] = 0.7;
            parameters[6] = 3.;
            parameters[7] = 6.;
            parameters[8] = 2.;
            parameters[9] = 18.;
        } else if (centralityCutNumber.CompareTo("312") == 0 ||
                centralityCutNumber.CompareTo("612") == 0 ){ //5-10%
            parameters[0] = 50.;
            parameters[1] = 1500.;
            parameters[2] = 4.;
            parameters[3] = 40.;
            parameters[4] = 0.07;
            parameters[5] = 0.7;
            parameters[6] = 3.;
            parameters[7] = 6.;
            parameters[8] = 2.;
            parameters[9] = 18.;
        } else if (centralityCutNumber.CompareTo("124") == 0 ||
                centralityCutNumber.CompareTo("524") == 0){ //20-40%
            parameters[0] = 50.;
            parameters[1] = 900.;
            parameters[2] = 3.;
            parameters[3] = 30.;
            parameters[4] = 0.06;
            parameters[5] = 0.5;
            parameters[6] = 3.;
            parameters[7] = 4.;
            parameters[8] = 2.;
            parameters[9] = 18.;
        } else if (centralityCutNumber.CompareTo("525") == 0){ //20-50%
            parameters[0] = 50.;
            parameters[1] = 900.;
            parameters[2] = 3.;
            parameters[3] = 30.;
            parameters[4] = 0.06;
            parameters[5] = 0.5;
            parameters[6] = 3.;
            parameters[7] = 4.;
            parameters[8] = 2.;
            parameters[9] = 18.;
        } else if (centralityCutNumber.CompareTo("146") == 0 ||
                centralityCutNumber.CompareTo("546") == 0){ //40-60%
            parameters[0] = 10.;
            parameters[1] = 200.;
            parameters[2] = 4.;
            parameters[3] = 40.;
            parameters[4] = 0.07;
            parameters[5] = 0.7;
            parameters[6] = 4.2;
            parameters[7] = 4.5;
            parameters[8] = 2.;
            parameters[9] = 18.;
        } else if (centralityCutNumber.CompareTo("168") == 0 ||
                centralityCutNumber.CompareTo("568") == 0){ //60-80%
            parameters[0] = 1.;
            parameters[1] = 80.;
            parameters[2] = 2.5;
            parameters[3] = 30.;
            parameters[4] = 0.05;
            parameters[5] = 0.5;
            parameters[6] = 2.;
            parameters[7] = 3.3;
            parameters[8] = 2.;
            parameters[9] = 18.;
        } else {
            parameters[0] = 1.;
            parameters[1] = 80.;
            parameters[2] = 2.5;
            parameters[3] = 30.;
            parameters[4] = 0.05;
            parameters[5] = 0.5;
            parameters[6] = 2.;
            parameters[7] = 3.3;
            parameters[8] = 2.;
            parameters[9] = 18.;
        }
        return;
    }

    //************************************************************************************
    //** Analyzes the name of the cent and returns the staring parameters for the fit ****
    //************************************************************************************
    void ReturnParameterSetFittingPbPbFromString(TString centralityCutNumber, Double_t* parameters){
        if (centralityCutNumber.CompareTo("0020") == 0 ||
            centralityCutNumber.CompareTo("0010") == 0 ||
            centralityCutNumber.CompareTo("1020") == 0  ){ //0-20%
            parameters[0] = 100.;
            parameters[1] = 1500.;
            parameters[2] = 4.;
            parameters[3] = 40.;
            parameters[4] = 0.07;
            parameters[5] = 0.7;
            parameters[6] = 3.;
            parameters[7] = 4.5;
            parameters[8] = 2.;
            parameters[9] = 18.;
            parameters[10] = 150.;
            parameters[11] = 11.5;
            parameters[12] = 0.135;
            parameters[13] = 4.3;
            parameters[14] = 7.6;
        } else if (centralityCutNumber.CompareTo("0005") == 0){ //0-5%
            parameters[0] = 50.;
            parameters[1] = 1500.;
            parameters[2] = 4.;
            parameters[3] = 40.;
            parameters[4] = 0.07;
            parameters[5] = 0.7;
            parameters[6] = 3.;
            parameters[7] = 6.;
            parameters[8] = 2.;
            parameters[9] = 18.;
            parameters[10] = 150.;
            parameters[11] = 11.5;
            parameters[12] = 0.135;
            parameters[13] = 4.3;
            parameters[14] = 7.6;
        } else if (centralityCutNumber.CompareTo("0510") == 0 ){ //5-10%
            parameters[0] = 50.;
            parameters[1] = 1500.;
            parameters[2] = 4.;
            parameters[3] = 40.;
            parameters[4] = 0.07;
            parameters[5] = 0.7;
            parameters[6] = 3.;
            parameters[7] = 6.;
            parameters[8] = 2.;
            parameters[9] = 18.;
            parameters[10] = 150.;
            parameters[11] = 11.5;
            parameters[12] = 0.135;
            parameters[13] = 4.3;
            parameters[14] = 7.6;
        } else if (centralityCutNumber.CompareTo("2040") == 0){ //20-40%
            parameters[0] = 50.;
            parameters[1] = 900.;
            parameters[2] = 3.;
            parameters[3] = 30.;
            parameters[4] = 0.06;
            parameters[5] = 0.5;
            parameters[6] = 3.;
            parameters[7] = 4.;
            parameters[8] = 2.;
            parameters[9] = 18.;
            parameters[10] = 90.;
            parameters[11] = 10.;
            parameters[12] = 0.135;
            parameters[13] = 3.5;
            parameters[14] = 7.7;
        } else if (centralityCutNumber.CompareTo("4060") == 0){ //40-60%
            parameters[0] = 10.;
            parameters[1] = 200.;
            parameters[2] = 4.;
            parameters[3] = 40.;
            parameters[4] = 0.07;
            parameters[5] = 0.7;
            parameters[6] = 4.2;
            parameters[7] = 4.5;
            parameters[8] = 2.;
            parameters[9] = 18.;
            parameters[10] = 27.;
            parameters[11] = 10.;
            parameters[12] = 0.135;
            parameters[13] = 4.4;
            parameters[14] = 7.5;
        } else if (centralityCutNumber.CompareTo("6080") == 0){ //60-80%
            parameters[0] = 1.;
            parameters[1] = 80.;
            parameters[2] = 2.5;
            parameters[3] = 30.;
            parameters[4] = 0.05;
            parameters[5] = 0.5;
            parameters[6] = 2.;
            parameters[7] = 3.3;
            parameters[8] = 2.;
            parameters[9] = 18.;
            parameters[10] = 10.;
            parameters[11] = 10.;
            parameters[12] = 0.135;
            parameters[13] = 3.;
            parameters[14] = 7.;
            cout << "bla here" << endl;
        } 	else {
        // 		Double_t parameter6080[10] = {10.,80.,2.5,30.,0.05,.5,3.,3.3,2.,18.};
            parameters[0] = 1.;
            parameters[1] = 80.;
            parameters[2] = 2.5;
            parameters[3] = 30.;
            parameters[4] = 0.05;
            parameters[5] = 0.5;
            parameters[6] = 2.;
            parameters[7] = 3.3;
            parameters[8] = 2.;
            parameters[9] = 18.;
            parameters[10] = 10.;
            parameters[11] = 10.;
            parameters[12] = 0.135;
            parameters[13] = 3.;
            parameters[14] = 7.;
        }
        return;
    }


    //************************************************************************************
    //** Analyzes the eventCutNumber and returns the centrality estimator name for files *
    //************************************************************************************
    TString GetCentralityEstimatorString(TString cutNumber){
        TString ppCutNumber                 = cutNumber(GetEventSystemCutPosition(),1);
        if (ppCutNumber.CompareTo("0") ==0){
            return "";
        } else if ( ppCutNumber.CompareTo("1") ==0 || ppCutNumber.CompareTo("3") ==0 || ppCutNumber.CompareTo("4") ==0 || ppCutNumber.CompareTo("5") ==0 || ppCutNumber.CompareTo("6") ==0 || ppCutNumber.CompareTo("7") ==0){
            return "_V0M";
        } else if ( ppCutNumber.CompareTo("8") ==0 || ppCutNumber.CompareTo("a") ==0 || ppCutNumber.CompareTo("c") ==0){
            return "_V0A";
        } else if ( ppCutNumber.CompareTo("2") ==0 || ppCutNumber.CompareTo("9") ==0 || ppCutNumber.CompareTo("b") ==0 || ppCutNumber.CompareTo("d") ==0 || ppCutNumber.CompareTo("g") ==0){
            return "_CL1";
        } else if ( ppCutNumber.CompareTo("e") ==0 || ppCutNumber.CompareTo("f") ==0){
            return "_ZNA";
        } else return "";
    }
    //************************************************************************************
    //** Analyzes the eventCutNumber and returns the centrality name with % for plotting *
    //************************************************************************************
    TString GetCentralityString(TString cutNumber){
        TString centralityCutNumberStart    = cutNumber(GetEventCentralityMinCutPosition(),1);
        TString centralityCutNumberEnd      = cutNumber(GetEventCentralityMaxCutPosition(),1);
        TString ppCutNumber                 = cutNumber(GetEventSystemCutPosition(),1);
        TString triggerNumber               = cutNumber(3,2);
        if (ppCutNumber.CompareTo("0") ==0 || ppCutNumber.CompareTo("h") ==0 || ppCutNumber.CompareTo("i") ==0 || ppCutNumber.CompareTo("j") ==0){
          if( triggerNumber.CompareTo("a0") ==0 ) { // INT7 CALOFAST
            return "calofastINT7";
          } else if ( triggerNumber.CompareTo("a1") ==0 ){ // EMC7 CALOFAST
            return "calofastEMC7";
          } else if ( triggerNumber.CompareTo("a2") ==0 ){ // EG2 CALOFAST
            return "calofastEG2";
          } else if ( triggerNumber.CompareTo("a3") ==0 ){ // EG1 CALOFAST
            return "calofastEG1";
          } else if ( triggerNumber.CompareTo("ap") ==0 ){ // PHI7 CALOFAST
            return "calofastPHI7";
          } else {
            return "pp";
          }
        } else if ( ppCutNumber.CompareTo("1") ==0 || ppCutNumber.CompareTo("5") ==0 || ppCutNumber.CompareTo("8") ==0 || ppCutNumber.CompareTo("4") ==0){
            if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
                return "0-100%";
            } else if( centralityCutNumberEnd.CompareTo("0")!=0){
                return Form("%i-%i%s", CutNumberToInteger(centralityCutNumberStart)*10,CutNumberToInteger(centralityCutNumberEnd)*10,"%");
            } else {
                if (centralityCutNumberEnd.CompareTo("0") == 0){
                    return Form("%i-100%s", CutNumberToInteger(centralityCutNumberStart)*10,"%");
                } else {
                    return Form("%i-%i%s", CutNumberToInteger(centralityCutNumberStart)*10,CutNumberToInteger(centralityCutNumberEnd)*10,"%");
                }
            }
        } else if ( ppCutNumber.CompareTo("2") ==0 || ppCutNumber.CompareTo("9") ==0 ){
            if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
                return "CL1 0-100%";
            } else if( centralityCutNumberEnd.CompareTo("0")!=0){
                return Form("CL1 %i-%i%s", CutNumberToInteger(centralityCutNumberStart)*10,CutNumberToInteger(centralityCutNumberEnd)*10,"%");
            } else {
                if (centralityCutNumberEnd.CompareTo("0") == 0){
                    return Form("CL1 %i-100%s", CutNumberToInteger(centralityCutNumberStart)*10,"%");
                } else {
                    return Form("CL1 %i-%i%s", CutNumberToInteger(centralityCutNumberStart)*10,CutNumberToInteger(centralityCutNumberEnd)*10,"%");
                }
            }
        } else if (ppCutNumber.CompareTo("3") ==0 || ppCutNumber.CompareTo("6") ==0 || ppCutNumber.CompareTo("7") ==0 || ppCutNumber.CompareTo("a") ==0){
            if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
                return "0-100%";
            } else {
                return Form("%i-%i%s", CutNumberToInteger(centralityCutNumberStart)*5,CutNumberToInteger(centralityCutNumberEnd)*5,"%");
            }
        } else if (ppCutNumber.CompareTo("c") ==0){
            if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
                return "0-100%";
            } else {
                return Form("%i-%i%s", CutNumberToInteger(centralityCutNumberStart),CutNumberToInteger(centralityCutNumberEnd),"%");
            }
        } else if ( ppCutNumber.CompareTo("e") ==0 || ppCutNumber.CompareTo("f") ==0){
            if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
                return "ZNA 0-100%";
            } else if( centralityCutNumberEnd.CompareTo("0")!=0){
                return Form("ZNA %i-%i%s", CutNumberToInteger(centralityCutNumberStart)*10,CutNumberToInteger(centralityCutNumberEnd)*10,"%");
            } else {
                if (centralityCutNumberEnd.CompareTo("0") == 0){
                    return Form("ZNA %i-100%s", CutNumberToInteger(centralityCutNumberStart)*10,"%");
                } else {
                    return Form("ZNA %i-%i%s", CutNumberToInteger(centralityCutNumberStart)*10,CutNumberToInteger(centralityCutNumberEnd)*10,"%");
                }
            }
        } else return "";
    }

    //************************************************************************************
    //** Analyzes the eventCutNumber and returns the centrality name for labeling w/o % **
    //************************************************************************************
    TString GetCentralityStringWoPer(TString cutNumber){
        TString centralityCutNumberStart    = cutNumber(GetEventCentralityMinCutPosition(),1);
        TString centralityCutNumberEnd      = cutNumber(GetEventCentralityMaxCutPosition(),1);
        TString ppCutNumber                 = cutNumber(GetEventSystemCutPosition(),1);
        if (ppCutNumber.CompareTo("0") ==0 || ppCutNumber.CompareTo("h") ==0 || ppCutNumber.CompareTo("i") ==0 || ppCutNumber.CompareTo("j") ==0){
            return "pp";
        } else if ( ppCutNumber.CompareTo("1") ==0 || ppCutNumber.CompareTo("5") ==0 || ppCutNumber.CompareTo("8") ==0 || ppCutNumber.CompareTo("4") ==0 ){
            if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
                return "0-100";
            } else if( centralityCutNumberEnd.CompareTo("0")!=0){
                return Form("%i-%i", CutNumberToInteger(centralityCutNumberStart)*10,CutNumberToInteger(centralityCutNumberEnd)*10);
            } else {
                if (centralityCutNumberEnd.CompareTo("0") == 0){
                    return Form("%i-100", CutNumberToInteger(centralityCutNumberStart)*10);
                } else {
                    return Form("%i-%i", CutNumberToInteger(centralityCutNumberStart)*10,CutNumberToInteger(centralityCutNumberEnd)*10);
                }
            }
        } else if ( ppCutNumber.CompareTo("2") ==0 || ppCutNumber.CompareTo("9") ==0){
            if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
                return "CL1_0-100";
            } else {
                return Form("CL1_%i-%i", CutNumberToInteger(centralityCutNumberStart)*10,CutNumberToInteger(centralityCutNumberEnd)*10);
            }
        } else if (ppCutNumber.CompareTo("3") ==0 || ppCutNumber.CompareTo("6") ==0 || ppCutNumber.CompareTo("7") ==0 || ppCutNumber.CompareTo("a") ==0){
            if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
                return "0-100";
            } else {
                return Form("%i-%i", CutNumberToInteger(centralityCutNumberStart)*5,CutNumberToInteger(centralityCutNumberEnd)*5);
            }
        } else if (ppCutNumber.CompareTo("c") ==0 ){
            if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
                return "0-100";
            } else {
                return Form("%i-%i", CutNumberToInteger(centralityCutNumberStart),CutNumberToInteger(centralityCutNumberEnd));
            }
        } else if ( ppCutNumber.CompareTo("e") ==0 || ppCutNumber.CompareTo("f") ==0){
            if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
                return "ZNA 0-100";
            } else {
                return Form("ZNA %i-%i", CutNumberToInteger(centralityCutNumberStart),CutNumberToInteger(centralityCutNumberEnd));
            }
        } else return "";
    }

    //************************************************************************************
    //** Analyzes the eventCutNumber and returns the centrality name as joing numbers ****
    //************************************************************************************
    TString GetCentralityStringOutput(TString cutNumber){
        TString centralityCutNumberStart    = cutNumber(GetEventCentralityMinCutPosition(),1);
        TString centralityCutNumberEnd      = cutNumber(GetEventCentralityMaxCutPosition(),1);
        TString ppCutNumber                 = cutNumber(GetEventSystemCutPosition(),1);
        if (ppCutNumber.CompareTo("0") ==0 || ppCutNumber.CompareTo("h") ==0 || ppCutNumber.CompareTo("i") ==0 || ppCutNumber.CompareTo("j") ==0){
            return "pp";
        } else if ( ppCutNumber.CompareTo("1") ==0 || ppCutNumber.CompareTo("5") ==0 || ppCutNumber.CompareTo("8") ==0 || ppCutNumber.CompareTo("4") ==0 ) {
            if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
                return "00100";
            } else {
                if (centralityCutNumberStart.CompareTo("0") == 0){
                    return Form("00%i", CutNumberToInteger(centralityCutNumberEnd)*10);
                } else if (centralityCutNumberEnd.CompareTo("0") == 0){
                    return Form("%i100", CutNumberToInteger(centralityCutNumberStart)*10);
                } else {
                    return Form("%i%i", CutNumberToInteger(centralityCutNumberStart)*10,CutNumberToInteger(centralityCutNumberEnd)*10);
                }
            }
        } else if ( ppCutNumber.CompareTo("2") ==0 || ppCutNumber.CompareTo("9") ==0){
            if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
                return "CL1_00100";
            } else {
                if (centralityCutNumberStart.CompareTo("0") == 0){
                    return Form("CL1_00%i", CutNumberToInteger(centralityCutNumberEnd)*10);
                } else if (centralityCutNumberEnd.CompareTo("0") == 0){
                    return Form("CL1_%i100", CutNumberToInteger(centralityCutNumberStart)*10);
                } else {
                    return Form("CL1_%i%i", CutNumberToInteger(centralityCutNumberStart)*10,CutNumberToInteger(centralityCutNumberEnd)*10);
                }
            }
        } else if (ppCutNumber.CompareTo("3") ==0 || ppCutNumber.CompareTo("6") ==0 || ppCutNumber.CompareTo("7") ==0 || ppCutNumber.CompareTo("a") ==0){
            if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
                return "00100";
            } else {
                if (centralityCutNumberStart.CompareTo("0") == 0){
                    return Form("00%i", CutNumberToInteger(centralityCutNumberEnd)*5);
                } else {
                    return Form("%i%i", CutNumberToInteger(centralityCutNumberStart)*5,CutNumberToInteger(centralityCutNumberEnd)*5);
                }
            }
        } else if (ppCutNumber.CompareTo("c") ==0 ){
            if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
                return "00100";
            } else {
                if (centralityCutNumberStart.CompareTo("0") == 0){
                    return Form("00%i", CutNumberToInteger(centralityCutNumberEnd));
                } else {
                    return Form("%i%i", CutNumberToInteger(centralityCutNumberStart),CutNumberToInteger(centralityCutNumberEnd));
                }
            }
        } else if ( ppCutNumber.CompareTo("e") ==0 || ppCutNumber.CompareTo("f") ==0){
            if (centralityCutNumberStart.CompareTo("0") == 0 && centralityCutNumberEnd.CompareTo("0") == 0  ){
                return "ZNA_00100";
            } else {
                if (centralityCutNumberStart.CompareTo("0") == 0){
                    return Form("ZNA_00%i", CutNumberToInteger(centralityCutNumberEnd)*10);
                } else if (centralityCutNumberEnd.CompareTo("0") == 0){
                    return Form("ZNA_%i100", CutNumberToInteger(centralityCutNumberStart)*10);
                } else {
                    return Form("ZNA_%i%i", CutNumberToInteger(centralityCutNumberStart)*10,CutNumberToInteger(centralityCutNumberEnd)*10);
                }
            }
        } else return "";
    }
    //************************************************************************************
    //** Analyzes the TPC dEdx E electron cut, return correct cut label ******************
    //************************************************************************************
    TString AnalyseTPCdEdxCutElectronLineTPC(Int_t ededxSigmaCut){
        switch(ededxSigmaCut){
            case 0: // -10,10  kTrue
                return "-10 #sigma < TPC dE/dx - #LTdE/dx#GT < 10 #sigma";
            case 1: // -10,10 kFalse
                return "-10 #sigma < TPC dE/dx - #LTdE/dx#GT < 10 #sigma";
            case 2:// -6,7
                return "-6 #sigma < TPC dE/dx - #LTdE/dx#GT < 7 #sigma";
            case 3: // -5,5
                return "-5 #sigma < TPC dE/dx - #LTdE/dx#GT < 5 #sigma";
            case 4: // -4,5
                return "-4 #sigma < TPC dE/dx - #LTdE/dx#GT < 5 #sigma";
            case 5: // -3,5
                return "-3 #sigma < TPC dE/dx - #LTdE/dx#GT < 5 #sigma";
            case 6: // -4,4
                return "-4 #sigma < TPC dE/dx - #LTdE/dx#GT < 4 #sigma";
            case 7: // -2.5,4
                return "-2.5 #sigma < TPC dE/dx - #LTdE/dx#GT < 4 #sigma";
            case 8: // -2,3.5
                return "-2 #sigma < TPC dE/dx - #LTdE/dx#GT 3.5 < #sigma";
            case 9: // -3,4
                return "-3 #sigma < TPC dE/dx - #LTdE/dx#GT < 4 #sigma";
            default:
                return "no dEdx cut defined";
        }
        return kTRUE;
    }

    //************************************************************************************
    //** Analyzes the M_{e^+e^-}, electron cut, return correct cut label ******************
    //************************************************************************************
    TString AnalyseMassEPCut(Int_t ReconstructionMassEPCut){
        switch(ReconstructionMassEPCut){
            case 0: // -999, 999-999
                return "No cut";
            case 1: // -999, 0.135
                return "M_{e^{+}e^{-}}<0.0135";
            case 2: // -999, 0.100
                return "M_{e^{+}e^{-}}<0.1";
            case 3: // 1.0, 0.015-0.032
                return "p_{T}(M_{e^{+}e^{-}})<1.0 and M_{e^{+}e^{-}}<0.015 or p_{T}(M_{e^{+}e^{-}})>1.0 and M_{e^{+}e^{-}}<0.032";
            case 4: // -999, 0.050
                return "M_{e^{+}e^{-}}<0.05";
            case 5: // -999, 0.035
                return "M_{e^{+}e^{-}}<0.035";
            case 6: // -999, 0.015
                return "M_{e^{+}e^{-}}<0.015";
            case 7: // 1.0, 0.015-0.035
                return "p_{T}(M_{e^{+}e^{-}})<1.0 and M_{e^{+}e^{-}}<0.015 or p_{T}(M_{e^{+}e^{-}})>1.0 and M_{e^{+}e^{-}}<0.035";
            case 8: // 1.0, 0.015-0.050
                return "p_{T}(M_{e^{+}e^{-}})<1.0 and M_{e^{+}e^{-}}<0.015 or p_{T}(M_{e^{+}e^{-}})>1.0 and M_{e^{+}e^{-}}<0.05";
            case 9: // 1.0, 0.025-0.035
                return "p_{T}(M_{e^{+}e^{-}})<1.0 and M_{e^{+}e^{-}}<0.025 or p_{T}(M_{e^{+}e^{-}})>1.0 and M_{e^{+}e^{-}}<0.035";
            case 10: // 1.0, 0.02-0.03
                return "p_{T}(M_{e^{+}e^{-}})<1.0 and M_{e^{+}e^{-}}<0.02 or p_{T}(M_{e^{+}e^{-}})>1.0 and M_{e^{+}e^{-}}<0.03";
            case 11: // 1.0, 0.027-0.054
                return "p_{T}(M_{e^{+}e^{-}})<1.0 and M_{e^{+}e^{-}}<0.027 or p_{T}(M_{e^{+}e^{-}})>1.0 and M_{e^{+}e^{-}}<0.054";
            case 12: // 1.0, 0.02-0.02
                return "M_{e^{+}e^{-}}<0.02";
            default:
                return "no Mass cut defined for Electron and Positron";
        }
        return kTRUE;
    }
    //************************************************************************************
    //** Analyzes the TPC cluster cut, return correct cut label **************************
    //************************************************************************************
    TString AnalyseTPCClusterCut(Int_t clsTPCCut){
        switch(clsTPCCut){
            case 0: // 0
                return "min. TPC cl.: 0";
            case 1:  // 60
                return "min. TPC cl.r: 60";
            case 2:  // 80
                return "min. TPC cl.: 80";
            case 3:  // 100
                return "min. TPC cl.: 100";
            case 4:  // 95% of findable clusters
                return "TPC cl./uncorr find. cl.: 0.95";
            case 5:  // 0% of findable clusters
                return "TPC cl./corr find. cl.: 0.";
            case 6:  // 80% of findable clusters
                return "TPC cl./corr find. cl.: 0.7";
            case 7:  // 0% of findable clusters
                return "TPC cl./uncorr find. cl.: 0.35";
            case 8:
                return "TPC cl./corr find. cl.: 0.35";
            case 9:
                return "TPC cl./corr find. cl.: 0.6";
            default:
                return "no cluster cut defined";
        }
    }

    //************************************************************************************
    //** Analyzes the TPC dEdx electron cut, return correct cut label ********************
    //************************************************************************************
    TString AnalyseTPCdEdxCutElectronLine(Int_t ededxSigmaCut){
        switch(ededxSigmaCut){
            case 0: // -10,10
                return "-10 < #sigma_{e} < 10";
            case 1: // -5,5
                return "-5 < #sigma_{e} < 5";
            case 2: // -3,5
                return "-3 < #sigma_{e} < 5";
            case 3: // -4,5
                return "-4 < #sigma_{e} < 5";
            case 4: // -6,7
                return "-6 < #sigma_{e} < 7";
            case 5: // -4,4
                return "-4 < #sigma_{e} < 4";
            case 6: // -2.5,4
                return "-2.5 < #sigma_{e} < 4";
            case 7: // -2,3.5
                return "-2 < #sigma_{e} < 3.5";
            case 8: // -2.5,3
                return "-2.5 < #sigma_{e} < 3";
            case 9: // -2.5,5
                return "-2.5 < #sigma_{e} < 5";
            default:
                return "no dEdx cut defined";
        }
        return kTRUE;
    }

    //************************************************************************************
    //** Analyzes the TPC dEdx pion cuts, return correct cut label ***********************
    //************************************************************************************
    TString AnalyseTPCdEdxCutPionLine(TString sPionCut){
        cout << sPionCut << endl;
        TString sPidedxSigmaCut                 = sPionCut(0,1);
        Int_t pidedxSigmaCut                    = CutNumberToInteger(sPidedxSigmaCut);
        cout << "pidedxSigmaCut: " << pidedxSigmaCut << endl;
        Double_t fPIDnSigmaAbovePionLine        = 0;
        Double_t fPIDnSigmaAbovePionLineHighPt  = 0;
        switch(pidedxSigmaCut){
            case 0:  // -10
                fPIDnSigmaAbovePionLine         = -10;
                fPIDnSigmaAbovePionLineHighPt   = -10;
                break;
            case 1:   // 0
                fPIDnSigmaAbovePionLine         = 0;
                fPIDnSigmaAbovePionLineHighPt   = -10;
                break;
            case 2:  // 1
                fPIDnSigmaAbovePionLine         = 1;
                fPIDnSigmaAbovePionLineHighPt   = -10;
                break;
            case 3:  // 1
                fPIDnSigmaAbovePionLine         = 2.5;
                fPIDnSigmaAbovePionLineHighPt   = -10;
                break;
            case 4:  // 1
                fPIDnSigmaAbovePionLine         = 3.;
                fPIDnSigmaAbovePionLineHighPt   = 1.;
                break;
            case 5:  // 1
                fPIDnSigmaAbovePionLine         = 2.;
                fPIDnSigmaAbovePionLineHighPt   = -10;
                break;
            case 6:  // 1
                fPIDnSigmaAbovePionLine         = 2.;
                fPIDnSigmaAbovePionLineHighPt   = 0.5;
                break;
            case 7:  // 1
                fPIDnSigmaAbovePionLine         = 3.5;
                fPIDnSigmaAbovePionLineHighPt   = -10;
                break;
            case 8:  // 1
                fPIDnSigmaAbovePionLine         = 2.;
                fPIDnSigmaAbovePionLineHighPt   = 1.;
                break;
            case 9:
                fPIDnSigmaAbovePionLine         = 3.0; // We need a bit less tight cut on dE/dx
                fPIDnSigmaAbovePionLineHighPt   = -10;
                break;
            default:
                cout << "pion line cut unknown" << endl;
        }

        TString sPiMomdedxSigmaCut              = sPionCut(1,1);
        Int_t piMomdedxSigmaCut                 = CutNumberToInteger(sPiMomdedxSigmaCut);
        Double_t fPIDMinPnSigmaAbovePionLine    = 0;
        cout << "piMomdedxSigmaCut: " << piMomdedxSigmaCut << endl;
        switch(piMomdedxSigmaCut){
            case 0:  // 0.5 GeV
                fPIDMinPnSigmaAbovePionLine     = 0.5;
                break;
            case 1:  // 1. GeV
                fPIDMinPnSigmaAbovePionLine     = 1.;
                break;
            case 2:  // 1.5 GeV
                fPIDMinPnSigmaAbovePionLine     = 1.5;
                break;
            case 3:  // 20.0 GeV
                fPIDMinPnSigmaAbovePionLine     = 20.;
                break;
            case 4:  // 50.0 GeV
                fPIDMinPnSigmaAbovePionLine     = 50.;
                break;
            case 5:  // 0.3 GeV
                fPIDMinPnSigmaAbovePionLine     = 0.3;
                break;
            case 6:  // 0.25 GeV
                fPIDMinPnSigmaAbovePionLine     = 0.25;
                break;
            case 7:  // 0.4 GeV
                fPIDMinPnSigmaAbovePionLine     = 0.4;
                break;
            case 8:  // 0.2 GeV
                fPIDMinPnSigmaAbovePionLine     = 0.2;
                break;
            default:
                cout << "pion line minimum pt cut unknown" << endl;
        }

        TString sPiMaxMomdedxSigmaCut           = sPionCut(2,1);
        Int_t piMaxMomdedxSigmaCut              = CutNumberToInteger(sPiMaxMomdedxSigmaCut);
        Double_t fPIDMaxPnSigmaAbovePionLine    = 0;
        cout << "piMaxMomdedxSigmaCut: " << piMaxMomdedxSigmaCut << endl;
        switch(piMaxMomdedxSigmaCut){
            case 0:  // 100. GeV
                fPIDMaxPnSigmaAbovePionLine     = 100.;
                break;
            case 1:  // 5. GeV
                fPIDMaxPnSigmaAbovePionLine     = 5.;
                break;
            case 2:  // 4. GeV
                fPIDMaxPnSigmaAbovePionLine     = 4.;
                break;
            case 3:  // 3.5 GeV
                fPIDMaxPnSigmaAbovePionLine     = 3.5;
                break;
            case 4:  // 3. GeV
                fPIDMaxPnSigmaAbovePionLine     = 3.;
                break;
            case 5:  // 7. GeV
                fPIDMaxPnSigmaAbovePionLine     = 7.;
                break;
            case 6:  // 2. GeV
                fPIDMaxPnSigmaAbovePionLine     = 2.;
                break;
            default:
                cout << "pion line minimum pt cut unknown" << endl;
        }
        return Form("rejected #pi for #sigma_{#pi} < %.2f (%.2f GeV/c < p_{T}_{#pi} < %.2f GeV/c), #sigma_{#pi} < %.2f (p_{T}_{#pi} > %.2f GeV/c)",fPIDnSigmaAbovePionLine, fPIDMinPnSigmaAbovePionLine,fPIDMaxPnSigmaAbovePionLine,fPIDnSigmaAbovePionLineHighPt,fPIDMaxPnSigmaAbovePionLine);
    }


    //************************************************************************************
    //** Analyzes the TOF electron PID cut, return correct cut label *********************
    //************************************************************************************
    TString AnalyseTOFelectronPIDCut(Int_t TOFpidCut){
        switch(TOFpidCut){
            case 0: // 0
                return "-100 < #sigma^{TOF}_{e} < 100";
                break;
            case 1:
                return "-7 < #sigma^{TOF}_{e} < 7";
                break;
            case 2:
                return "-5 < #sigma^{TOF}_{e} < 5";
                break;
            case 3:
                return "-3 < #sigma^{TOF}_{e} < 5";
                break;
            case 4:
                return "-2 < #sigma^{TOF}_{e} < 3";
                break;
            case 5:
                return "-3 < #sigma^{TOF}_{e} < 3";
                break;
            default:
                return "no TOF cut defined";
        }
    }


    //************************************************************************************
    //********* Analyzes the chi2 gamma cut, return correct cut label ********************
    //************************************************************************************
    TString AnalyseChi2GammaCut(    Int_t chi2GammaCut,
                                    Int_t psiPairCut){   // Set Cut

        TString psiPairCutString            = "";
        Bool_t k2DPsiPairChi2               = kFALSE;
        switch(psiPairCut) {
            case 0:
                psiPairCutString = "|#Psi_{Pair}| < 10000";
                break;
            case 1:
                psiPairCutString = "|#Psi_{Pair}| < 0.1";
                break;
            case 2:
                psiPairCutString = "|#Psi_{Pair}| < 0.05";
                break;
            case 3:
                psiPairCutString = "|#Psi_{Pair}| < 0.035";
                break;
            case 4:
                psiPairCutString = "|#Psi_{Pair}| < 0.2";
                break;
            case 5:
                psiPairCutString = "|#Psi_{Pair}| < 0.1";
                k2DPsiPairChi2= kTRUE;
                break;
            case 6:
                psiPairCutString = "|#Psi_{Pair}| < 0.05";
                k2DPsiPairChi2= kTRUE;
                break;
            case 7:
                psiPairCutString = "|#Psi_{Pair}| < 0.035";
                k2DPsiPairChi2= kTRUE;
                break;
            case 8:
                psiPairCutString =  "|#Psi_{Pair}| < 0.2";
                k2DPsiPairChi2= kTRUE;
                break;
            case 9:
                psiPairCutString =  "|#Psi_{Pair}| < 0.5";
                break;
            default:
                psiPairCutString =  "#Psi_{Pair} cut not defined";
                break;
        }

        TString chi2CutString               = "";
        switch(chi2GammaCut){
        case 0: // 100
            chi2CutString = "#chi_{#gamma}^{2} < 100";
            break;
        case 1:  // 50
            chi2CutString = "#chi_{#gamma}^{2} < 50";
            break;
        case 2:  // 30
            chi2CutString = "#chi_{#gamma}^{2} < 30";
            break;
        case 3:
            chi2CutString = "#chi_{#gamma}^{2} < 200";
            break;
        case 4:
            chi2CutString = "#chi_{#gamma}^{2} < 500";
            break;
        case 5:
            chi2CutString = "#chi_{#gamma}^{2} < 1000000";
            break;
        case 6:
            chi2CutString = "#chi_{#gamma}^{2} < 5";
            break;
        case 7:
            chi2CutString = "#chi_{#gamma}^{2} < 10";
            break;
        case 8:
            chi2CutString = "#chi_{#gamma}^{2} < 20";
            break;
        case 9:
            chi2CutString = "#chi_{#gamma}^{2} < 15";
            break;
        case 10:
            chi2CutString = "#chi_{#gamma}^{2} < 25";
            break;
        default:
            chi2CutString = "#chi_{#gamma}^{2} cut unknown";
            break;
        }
        if (k2DPsiPairChi2){
            return Form("2D cut: %s, %s", chi2CutString.Data(), psiPairCutString.Data());
        } else {
            return Form("1D cut: %s, 1D cut %s", chi2CutString.Data(), psiPairCutString.Data());
        }
        return "";
    }

    //************************************************************************************
    //********* Analyzes the qt gamma cut, return correct cut label **********************
    //************************************************************************************
    TString AnalyseQtMaxCut(Int_t QtMaxCut){

        switch(QtMaxCut){
            case 0: //
                return "no q_{T}_{#gamma} cut applied";
            case 1:
                return "q_{T}_{#gamma} < 0.1 GeV/c";
            case 2:
                return "2D ellipse q_{T}_{#gamma} < 0.06 GeV/c, #alpha < 0.95";
            case 3:
                return "q_{T}_{#gamma} < 0.05 GeV/c";
            case 4:
                return "q_{T}_{#gamma} < 0.03 GeV/c";
            case 5:
                return "q_{T}_{#gamma} < 0.02 GeV/c";
            case 6:
                return "2D ellipse q_{T}_{#gamma} < 0.02 GeV/c, #alpha < 0.95";
            case 7:
                return "q_{T}_{#gamma} < 0.15 GeV/c";
            case 8:
                return "2D ellipse q_{T}_{#gamma} < 0.05 GeV/c, #alpha < 0.95";
            case 9:
                return "2D ellipse q_{T}_{#gamma} < 0.03 GeV/c, #alpha < 0.95";
            default:
                return "no q_{T}_{#gamma} cut defined";
        }
    }

    //************************************************************************************
    //********* Analyzes the single leg pt cut, return correct cut label *****************
    //************************************************************************************
    TString AnalyseSinglePtCut(Int_t singlePtCut, Bool_t printGammaPt = kFALSE){

        if(printGammaPt){
            switch(singlePtCut){
            case 0: // 0.050 GeV + min gamma pT cut of 20 MeV
                return "p_{T}_{e^{#pm}} > 0.050, p_{T}_{#gamma} > 0.020 GeV/c";
            case 1:  // 0.100 GeV + min gamma pT cut of 20 MeV
                return "p_{T}_{e^{#pm}} > 0.100, p_{T}_{#gamma} > 0.020 GeV/c";
            case 2:  // 0.150 GeV + min gamma pT cut of 20 MeV
                return "p_{T}_{e^{#pm}} > 0.150, p_{T}_{#gamma} > 0.020 GeV/c";
            case 3:  // 0.200 GeV + min gamma pT cut of 20 MeV
                return "p_{T}_{e^{#pm}} > 0.200, p_{T}_{#gamma} > 0.020 GeV/c";
            case 4:  // 0.075 GeV + min gamma pT cut of 20 MeV
                return "p_{T}_{e^{#pm}} > 0.075, p_{T}_{#gamma} > 0.020 GeV/c";
            case 5:  // 0.125 GeV + min gamma pT cut of 20 MeV
                return "p_{T}_{e^{#pm}} > 0.125, p_{T}_{#gamma} > 0.020 GeV/c";
            case 6:  // 0.04 GeV  + min gamma pT cut of 10 MeV
                return "p_{T}_{e^{#pm}} > 0.040, p_{T}_{#gamma} > 0.010 GeV/c";
            case 7:  // 0.0 GeV  + min gamma pT cut of 0 MeV
                return "p_{T}_{e^{#pm}} > 0.0, p_{T}_{#gamma} > 0.0 GeV/c";
            case 8:  // 0.02 GeV + min gamma pT cut of 10 MeV
                return "p_{T}_{e^{#pm}} > 0.020, p_{T}_{#gamma} > 0.010 GeV/c";
            case 9: // 0.050 GeV + min gamma pT cut of 100 MeV
                return "p_{T}_{e^{#pm}} > 0.050, p_{T}_{#gamma} > 0.100 GeV/c";
            default:
                return "p_{T} cut not defined";
            }
        }else{
            switch(singlePtCut){
            case 0: // 0.050 GeV
                return "p_{T}_{e^{#pm}} > 0.050 GeV/c";
            case 1:  // 0.100 GeV
                return "p_{T}_{e^{#pm}} > 0.100 GeV/c";
            case 2:  // 0.150 GeV
                return "p_{T}_{e^{#pm}} > 0.150 GeV/c";
            case 3:  // 0.200 GeV
                return "p_{T}_{e^{#pm}} > 0.200 GeV/c";
            case 4:  // 0.075 GeV
                return "p_{T}_{e^{#pm}} > 0.075 GeV/c";
            case 5:  // 0.125 GeV
                return "p_{T}_{e^{#pm}} > 0.125 GeV/c";
            case 6:  // 0.04 GeV
                return "p_{T}_{e^{#pm}} > 0.040 GeV/c";
            case 7:  // 0.0 GeV
                return "p_{T}_{e^{#pm}} > 0.0 GeV/c";
            case 8:  // 0.02 GeV
                return "p_{T}_{e^{#pm}} > 0.020 GeV/c";
            default:
                return "p_{T}_{e^{#pm}} cut not defined";
            }
        }
    }

    //************************************************************************************
    //********* Analyzes the dca z photon cu, return correct cut label *******************
    //************************************************************************************
    TString AnalyseDCAZPhotonCut(Int_t dcaZPhoton){
        switch(dcaZPhoton){
            case 0:  //
                return "|dca_{Z}| < 1000 cm";
            case 1:  //
                return "|dca_{Z}| < 10 cm";
            case 2:  //
                return "|dca_{Z}| < 5 cm";
            case 3:  //
                return "|dca_{Z}| < 4 cm";
            case 4:  //
                return "|dca_{Z}| < 3 cm";
            case 5:  //
                return "|dca_{Z}| < 2.5 cm";
            case 6:  //
                return "|dca_{Z}| < 2 cm";
            case 7:  //
                return "|dca_{Z}| < 1.5 cm";
            case 8:  //
                return "|dca_{Z}| < 1 cm";
            case 9:  //
                return "|dca_{Z}| < 0.5 cm";
            default:
                return "|dca_{Z}| cut not defined";
        }
    }

    //************************************************************************************
    //********* Analyzes the dca z gamma cut, return max dca value ***********************
    //************************************************************************************
    Double_t AnalyseDCAZPhotonCutValue(Int_t dcaZPhoton){
        switch(dcaZPhoton){
            case 0:  //
                return 1000.;
            case 1:  //
                return 10.;
            case 2:  //
                return 5.;
            case 3:  //
                return 4.;
            case 4:  //
                return 3.;
            case 5:  //
                return 2.5;
            case 6:  //
                return 2.;
            case 7:  //
                return 1.5;
            case 8:  //
                return 1.;
            case 9:  //
                return 0.5;
            default:
                return 1000;
        }
    }

    //************************************************************************************
    //********* Analyzes the cos(theta_point) gamma cut, return correct cut label ********
    //************************************************************************************
    TString AnalyseCosPointCut(Int_t cosPoint){
        switch(cosPoint){
            case 0:  //
                return "cos(#Theta_{point}) > -1";
            case 1:  //
                return "cos(#Theta_{point}) > 0";
            case 2:  //
                return "cos(#Theta_{point}) > 0.5";
            case 3:  //
                return "cos(#Theta_{point}) > 0.75";
            case 4:  //
                return "cos(#Theta_{point}) > 0.85";
            case 5:  //
                return "cos(#Theta_{point}) > 0.88";
            case 6:  //
                return "cos(#Theta_{point}) > 0.9";
            case 7:  //
                return "cos(#Theta_{point}) > 0.95";
            default:
                return "cos(#Theta_{point}) cut not defined";
        }
    }

    //************************************************************************************
    //********* Analyzes the eta gamma/electron cut, return correct cut label ************
    //************************************************************************************
    TString AnalyseEtaCut(Int_t etaCut){

        switch(etaCut){
            case 0: // 0.9
                return "#eta_{#gamma,e^{#pm}} < 0.9";
            case 1:  // 0.6
                return "#eta_{#gamma,e^{#pm}} < 0.6";
            case 2:  // 1.4
                return "#eta_{#gamma,e^{#pm}} < 1.4";
            case 3: // 0.65
                return "#eta_{#gamma,e^{#pm}} < 0.65";
            case 4: // 0.75
                return "#eta_{#gamma,e^{#pm}} < 0.75";
            case 5: // 0.5
                return "#eta_{#gamma,e^{#pm}} < 0.5";
            case 6: // 5.
                return "#eta_{#gamma,e^{#pm}} < 5";
            case 7: // 0.1 - 0.8
                return "#eta_{#gamma,e^{#pm}} < 0.7";
            case 8: // 0.1 - 0.8
                return "#eta_{#gamma,e^{#pm}} < 0.4";
            case 9: // 10
                return "#eta_{#gamma,e^{#pm}} < 10";
            default:
                return "no #eta_{#gamma,e^{#pm}} cut defined";
        }
    }

    //************************************************************************************
    //****** Analyzes the eta gamma/electron cut for pPb, return correct cut label *******
    //************************************************************************************
    TString AnalyseEtaCutpPb(Int_t etaCut){

        switch(etaCut){
            case 0: // 0.9
                return "#eta_{#gamma,e^{#pm}} < 0.9";
            case 1:  // 1.2
                return "#eta_{#gamma,e^{#pm}} < 0.6";
            case 2:  // 1.4
                return "#eta_{#gamma,e^{#pm}} < 1.4";
            case 3: // 0.8
                return "#eta_{#gamma,e^{#pm}} < 0.8";
            case 4: // 0.75
                return "#eta_{#gamma,e^{#pm}} < 0.75";
            case 5: // 0.9 - 1.4
                return "#eta_{#gamma,e^{#pm}} < 0.5";
            case 6: // 5.
                return "#eta_{#gamma,e^{#pm}} < 5";
            case 7: // 0.1 - 0.8
                return "#eta_{#gamma,e^{#pm}} < 0.3";
            case 8: // 0.1 - 0.8
                return "#eta_{#gamma,e^{#pm}} < 0.4";
            case 9: // 10
                return "#eta_{#gamma,e^{#pm}} < 10";
            default:
                return "no #eta_{#gamma,e^{#pm}} cut defined";
        }
    }

    //************************************************************************************
    //**************** Analyzes the R gamma cut, return correct cut label ****************
    //************************************************************************************
    TString AnalyseRCut(Int_t RCut){
        // Set Cut
        switch(RCut){
            case 0:
                return "0 cm < R_{conv, #gamma} < 180 cm";
            case 1:
                return "2.8 cm < R_{conv, #gamma} < 180 cm";
            case 2:
                return "5 cm < R_{conv, #gamma} < 180 cm";
            case 3:
                return "10 cm < R_{conv, #gamma} < 70 cm";
            case 4:
                return "5 cm < R_{conv, #gamma} < 70 cm";
            case 5:
                return "10 cm < R_{conv, #gamma} < 180 cm";
            case 6:
                return "20 cm < R_{conv, #gamma} < 180 cm";
            case 7:
                return "35 cm < R_{conv, #gamma} < 180 cm"; // before 26 (19.4.2014)qdel
            case 8:
                return "12.5 cm < R_{conv, #gamma} < 180 cm";
            case 9:
                return "7.5 cm < R_{conv, #gamma} < 180 cm";
            default:
                return "R cut not defined";
        }
    }

    //************************************************************************************
    // Analyzes the R gamma cut, together with photonQuality cut, return correct cut label
    //************************************************************************************
    TString AnalyseRCutAndQuality(Int_t RCut, Int_t photonQualitCut){
        // Set Cut
        TString stringRCut = "";
        switch(RCut){
            case 0:
                stringRCut= "0 cm < R_{conv, #gamma} < 180 cm";
                break;
            case 1:
                stringRCut= "2.8 cm < R_{conv, #gamma} < 180 cm";
                break;
            case 2:
                stringRCut=  "5 cm < R_{conv, #gamma} < 180 cm";
                break;
            case 3:
                stringRCut= "10 cm < R_{conv, #gamma} < 70 cm";
                break;
            case 4:
                stringRCut= "5 cm < R_{conv, #gamma} < 70 cm";
                break;
            case 5:
                stringRCut= "10 cm < R_{conv, #gamma} < 180 cm";
                break;
            case 6:
                stringRCut= "20 cm < R_{conv, #gamma} < 180 cm";
                break;
            case 7:
                stringRCut= "35 cm < R_{conv, #gamma} < 180 cm"; // before 26 (19.4.2014)qdel
                break;
            case 8:
                stringRCut= "12.5 cm < R_{conv, #gamma} < 180 cm";
                break;
            case 9:
                stringRCut= "7.5 cm < R_{conv, #gamma} < 180 cm";
                break;
            default:
                stringRCut= "R cut not defined";
                break;
        }

        TString stringPhotonQuality = "";
        switch(photonQualitCut) {
            case 0:
                stringPhotonQuality =  "photon Quality: 1,2,3";
                break;
            case 2:
                stringPhotonQuality =  "photon Quality: 1";
                break;
            case 3:
                stringPhotonQuality =  "photon Quality: 2";
                break;
            case 4:
                stringPhotonQuality =  "photon Quality: 3";
                break;
            default:
                stringPhotonQuality =  "photon Quality cut not defined";
                break;
        }
        return Form("%s, %s", stringRCut.Data(), stringPhotonQuality.Data());

    }

    //************************************************************************************
    //************* Analyzes the psi_pair gamma cut, return correct cut label ************
    //************************************************************************************
    TString AnalysePsiPair(Int_t PsiPairCut, Int_t chi2GammaCut){

        TString psiPairCutString = "";
        Bool_t k2DPsiPairChi2 = kFALSE;
        switch(PsiPairCut) {
            case 0:
                psiPairCutString = "|#Psi_{Pair}| < 10000";
                break;
            case 1:
                psiPairCutString = "|#Psi_{Pair}| < 0.1";
                break;
            case 2:
                psiPairCutString = "|#Psi_{Pair}| < 0.05";
                break;
            case 3:
                psiPairCutString = "|#Psi_{Pair}| < 0.035";
                break;
            case 4:
                psiPairCutString = "|#Psi_{Pair}| < 0.2";
                break;
            case 5:
                psiPairCutString = "|#Psi_{Pair}| < 0.1";
                k2DPsiPairChi2= kTRUE;
                break;
            case 6:
                psiPairCutString = "|#Psi_{Pair}| < 0.05";
                k2DPsiPairChi2= kTRUE;
                break;
            case 7:
                psiPairCutString = "|#Psi_{Pair}| < 0.035";
                k2DPsiPairChi2= kTRUE;
                break;
            case 8:
                psiPairCutString =  "|#Psi_{Pair}| < 0.2";
                k2DPsiPairChi2= kTRUE;
                break;
            case 9:
                psiPairCutString =  "|#Psi_{Pair}| < 0.5";
                break;
            default:
                psiPairCutString =  "#Psi_{Pair} cut not defined";
                break;
        }

        TString chi2CutString = "";
        switch(chi2GammaCut){
            case 0: // 100
                chi2CutString = "#chi_{#gamma}^{2} < 100";
                break;
            case 1:  // 50
                chi2CutString = "#chi_{#gamma}^{2} < 50";
                break;
            case 2:  // 30
                chi2CutString = "#chi_{#gamma}^{2} < 30";
                break;
            case 3:
                chi2CutString = "#chi_{#gamma}^{2} < 200";
                break;
            case 4:
                chi2CutString = "#chi_{#gamma}^{2} < 500";
                break;
            case 5:
                chi2CutString = "#chi_{#gamma}^{2} < 1000000";
                break;
            case 6:
                chi2CutString = "#chi_{#gamma}^{2} < 5";
                break;
            case 7:
                chi2CutString = "#chi_{#gamma}^{2} < 10";
                break;
            case 8:
                chi2CutString = "#chi_{#gamma}^{2} < 20";
                break;
            case 9:
                chi2CutString = "#chi_{#gamma}^{2} < 15";
                break;
            default:
                chi2CutString = "#chi_{#gamma}^{2} cut unknown";
                break;
        }
        if (k2DPsiPairChi2){
            return Form("2D cut: %s, %s", psiPairCutString.Data(), chi2CutString.Data());
        } else {
            return Form("1D cut: %s", psiPairCutString.Data());
        }
        return "";
    }

    //************************************************************************************
    //******** Analyzes the psi_pair and R gamma cut, return correct cut label ***********
    //************************************************************************************
    TString AnalysePsiPairAndR(Int_t PsiPairCut, Int_t RCut ){
        TString psiPairCut = "";
        switch(PsiPairCut) {
            case 0:
                psiPairCut = "|#Psi_{Pair}| < 10000";
                break;
            case 1:
                psiPairCut = "|#Psi_{Pair}| < 0.1";
                break;
            case 2:
                psiPairCut = "|#Psi_{Pair}| < 0.05";
                break;
            case 3:
                psiPairCut = "|#Psi_{Pair}| < 0.035";
                break;
            case 4:
                psiPairCut = "|#Psi_{Pair}| < 0.15";
                break;
            case 5:
                psiPairCut = "#Psi_{Pair} < (0.1 - 0.1/1 * #Delta #Phi )";
                break;
            case 6:
                psiPairCut = "#Psi_{Pair} < (0.05 - 0.05/1 * #Delta #Phi )";
                break;
            case 7:
                psiPairCut = "#Psi_{Pair} < (0.035 - 0.035/1 * #Delta #Phi )";
                break;
            case 8:
                psiPairCut = "#Psi_{Pair} < (0.2 - 0.2/1 * #Delta #Phi )";
                break;
            case 9:
                psiPairCut = "#Psi_{Pair} < 0.5";
                break;
            default:
                psiPairCut = "#Psi_{Pair} cut not defined";
                break;
        }
        TString RCutString = "";
        switch(RCut){
            case 0:
                RCutString =  "0 cm < R_{conv, #gamma} < 180 cm";
                break;
            case 1:
                RCutString =  "2.8 cm < R_{conv, #gamma} < 180 cm";
                break;
            case 2:
                RCutString =  "5 cm < R_{conv, #gamma} < 180 cm";
                break;
            case 3:
                RCutString =  "10 cm < R_{conv, #gamma} < 70 cm";
                break;
            case 4:
                RCutString =  "5 cm < R_{conv, #gamma} < 70 cm";
                break;
            case 5:
                RCutString =  "10 cm < R_{conv, #gamma} < 180 cm";
                break;
            case 6:
                RCutString =  "20 cm < R_{conv, #gamma} < 180 cm";
                break;
            case 7:
                RCutString =  "26 cm < R_{conv, #gamma} < 180 cm";
                break;
            case 8:
                RCutString =  "35 cm < R_{conv, #gamma} < 180 cm";
                break;
            case 9:
                RCutString =  "5 cm < R_{conv, #gamma} < 35 cm";
                break;
            default:
                RCutString =  "R cut not defined";
                break;
        }
        return Form("%s, %s", psiPairCut.Data(), RCutString.Data());
    }

    //************************************************************************************
    //********** Analyzes the V0reader gamma cut, return correct cut label ***************
    //************************************************************************************
    TString AnalyseV0ReaderCut(Int_t V0ReaderCut){
        // Set Cut
        switch(V0ReaderCut){
            case 0:  //
                return "Onfly V0finder";
            case 1:  //
                return "Offline V0finder";
            default:
                return "V0finder cut not defined";
        }
    }


    //************************************************************************************
    //***** Analyzes the 2D psi pair & chi2 gamma cut, return correct cut label **********
    //************************************************************************************
    TString AnalyseChi2PsiPair(Int_t PsiPairCut ){
        switch(PsiPairCut) {
            case 26:
                return "2D: #chi_{#gamma}^{2} = 30 |#Psi_{Pair}| < 0.05";
            case 22:
                return "1D: #chi_{#gamma}^{2} = 30 |#Psi_{Pair}| < 0.05";
            case 16:
                return "2D: #chi_{#gamma}^{2} = 50 |#Psi_{Pair}| < 0.05";
            case 86:
                return "2D: #chi_{#gamma}^{2} = 20|#Psi_{Pair}| < 0.05";
            case 25:
                return "2D: #chi_{#gamma}^{2} = 30 |#Psi_{Pair}| < 0.1";
            case 27:
                return "2D: #chi_{#gamma}^{2} = 30 |#Psi_{Pair}| < 0.035";
            default:
                return " #chi_{#gamma}^{2} #Psi_{Pair} cut not defined";
        }
    }

    //************************************************************************************
    //***** Analyzes the photon quality gamma cut, return correct cut label **************
    //************************************************************************************
    TString AnalysePhotonQuality(Int_t photonQualitCut ){
        switch(photonQualitCut) {
            case 0:
                return "photon Quality: 1,2,3";
            case 2:
                return "photon Quality: 1 (TPC only photons)";
            case 3:
                return "photon Quality: 2 (1 leg with #geq 2 ITS hits)";
            case 4:
                return "photon Quality: 3 (both legs with #geq 2 ITS hits)";
            default:
                return "photon Quality cut not defined";
        }
    }

    //************************************************************************************
    //********* Analyzes the photon asymmetry cut, return correct cut label **************
    //************************************************************************************
    TString AnalysePhotonAsymmetry(Int_t photonQualitCut ){
        switch(photonQualitCut) {
            case 0:
                return "no cut";
            case 1:
                return "for p_{T,track} > 3.5,  A_{gamma} < 0.04";
            case 2:
                return "for p_{T,track} > 3.5,  A_{gamma} < 0.06";
            case 3:
                return "for p_{T,track} > 0.0,  A_{gamma} < 0.05";
            case 4:
                return "p-dependent asymmetry cut ";
            default:
                return "photon asymmetry cut not defined";
        }
    }

    //************************************************************************************
    //********* Analyzes the Sphericity cut, return correct cut label       **************
    //************************************************************************************
    TString AnalyseSphericity(TString string ){
        if(string.CompareTo("h0a") == 0){
          return "0 < S_{T} < 1";
        } else if(string.CompareTo("h05") == 0){
          return "S_{T} < 0.5";
        } else if(string.CompareTo("h5a") == 0){
          return "S_{T} > 0.5";
        } else if(string.CompareTo("h03") == 0){
          return "S_{T} < 0.3";
        } else if(string.CompareTo("h7a") == 0){
          return "S_{T} > 0.7";
        } else {
          return "Sphericity cut not defined";
        }
    }
    //************************************************************************************
    //********* Analyzes the Multiplicity cut, return correct cut label       ************
    //************************************************************************************
    TString AnalyseMultiplicity(TString string){
      if(string.CompareTo("m01") == 0){
        return "0-1% V0M";
      } else if(string.CompareTo("m02") == 0){
        return "0-2% V0M";
      } else if(string.CompareTo("m15") == 0){
        return "1-5% V0M";
      } else if(string.CompareTo("m05") == 0){
        return "0-5% V0M";
      } else if(string.CompareTo("m5k") == 0){
        return "5-20% V0M";
      } else if(string.CompareTo("n24") == 0){
        return "20-40% V0M";
      } else if(string.CompareTo("n26") == 0){
        return "20-60% V0M";
      } else if(string.CompareTo("n47") == 0){
        return "40-70% V0M";
      } else if(string.CompareTo("n6a") == 0){
        return "60-100% V0M";
      } else if(string.CompareTo("n7a") == 0){
        return "70-100% V0M";
      } else if(string.CompareTo("o01") == 0){
        return "0-1% SPD";
      } else if(string.CompareTo("o02") == 0){
        return "0-2% SPD";
      } else if(string.CompareTo("o05") == 0){
        return "0-5% SPD";
      } else if(string.CompareTo("o5k") == 0){
        return "5-20% SPD";
      } else if(string.CompareTo("p26") == 0){
        return "20-60% SPD";
      } else if(string.CompareTo("p6a") == 0){
        return "60-100% SPD";
      } else if(string.CompareTo("000") == 0){
        return "0-100%";
      } else {
        return "Mult cut not defined";
      }
    }

    //************************************************************************************
    //******** Analyzes the phi exclusion photon cuts, return correct cut label **********
    //************************************************************************************
    TString AnalyseConvPhiExclusionCut(TString gammaCutNumber){

        Int_t minCutNumberPhi   = GetPhotonMinPhiCutPosition(gammaCutNumber);
        if (minCutNumberPhi == -1) return "Full TPC acceptance";

        TString minPhiCutNumber(gammaCutNumber(minCutNumberPhi,1));
        TString maxPhiCutNumber(gammaCutNumber(GetPhotonMaxPhiCutPosition(gammaCutNumber),1));
        TString etaExclusion(gammaCutNumber(GetPhotonEtaForPhiCutPosition(gammaCutNumber),1));

        TString etaAcc          = "";
        if (etaExclusion.CompareTo("1") == 0) etaAcc = "#eta < 0" ;
            else if (etaExclusion.CompareTo("2") == 0) etaAcc = "#eta > 0" ;

        if ( etaExclusion.CompareTo("0") == 0 && minPhiCutNumber.CompareTo("0")==0 && maxPhiCutNumber.CompareTo("0")==0 ){
            return "Full TPC acceptance";
        } else if (minPhiCutNumber.CompareTo("0")==0 && maxPhiCutNumber.CompareTo("0")==0){
            return Form("Full $varphi acceptance, %s",etaAcc.Data()) ;
        }

        Double_t fMinPhiCut     = 0;
        Double_t fMaxPhiCut     = 2* TMath::Pi();
        switch(CutNumberToInteger(minPhiCutNumber)) {
            case 0:
                fMinPhiCut = 0; //no cut on phi
                break;
            case 1:
                fMinPhiCut = 1.7; //OROC C08
                break;
            case 2:
                fMinPhiCut = 4.4; //EMCal
                break;
            case 3:
                fMinPhiCut = 1.0; //PHOS
                break;
            case 4:
                fMinPhiCut = 3.4; //EMCal tight
                break;
            case 5:
                fMinPhiCut = 2.0; //OROC C08 medium cut
                break;
            case 6:
                fMinPhiCut = 2.2; //OROC C08 small cut
                break;
            case 7:
                fMinPhiCut = 2.4; //OROC C08 tightest cut
                break;

            default:
                fMinPhiCut = 0.;
                break;
        }

        switch (CutNumberToInteger(maxPhiCutNumber)) {
            case 0:
                fMaxPhiCut = 2* TMath::Pi(); //no cut
                break;
            case 1:
                fMaxPhiCut = 4.3; //OROC C08
                break;
            case 2:
                fMaxPhiCut = 5.8; //EMCal
                break;
            case 3:
                fMaxPhiCut = 3.0; //PHOS
                break;
            case 4:
                fMaxPhiCut = 1.; //EMCal
                break;
            case 5:
                fMaxPhiCut = 4.0; //OROC C08 medium cut
                break;
            case 6:
                fMaxPhiCut = 3.8; //OROC C08 small cut
                break;
            case 7:
                fMaxPhiCut = 3.6; //OROC C08 tighest cut
                break;

            default:
                fMaxPhiCut = 2* TMath::Pi(); //no cut
                break;
        }

        if (fMinPhiCut < fMaxPhiCut){
            if (etaAcc.CompareTo("") == 0)
                return Form("exclude #gamma_{conv} in %1.3f < #varphi < %1.3f", fMinPhiCut, fMaxPhiCut);
            else
                return Form("exclude #gamma_{conv} in %1.3f < #varphi < %1.3f, %s", fMinPhiCut, fMaxPhiCut, etaAcc.Data());
        } else {
            return Form("restrict #gamma_{conv} to %1.3f < #varphi < %1.3f",  fMaxPhiCut, fMinPhiCut );
        }

        return "";
    }


    //************************************************************************************
    //******** Analyzes the trigger event cut, return correct cut label ******************
    //************************************************************************************
    TString AnalyseSpecialTriggerCut(Int_t SpecialTrigger, TString periodName = "No"){
        // Set Cut
        switch(SpecialTrigger){
            case 0:
                if (periodName.CompareTo("LHC11a")==0||periodName.CompareTo("LHC11a_p4_woSDD")==0)
                    return "without SDD, V0OR";
                else if (periodName.Contains("LHC12") || periodName.Contains("LHC13") )
                    return "V0AND";
                else
                    return "V0OR";
            case 1:
                if (periodName.CompareTo("LHC11a")==0||periodName.CompareTo("LHC11a_p4_woSDD")==0)
                    return "without SDD, V0OR";
                else
                    return "V0OR";
            case 3:
                if (periodName.CompareTo("LHC11a")==0||periodName.CompareTo("LHC11a_p4_wSDD")==0)
                    return "with SDD, V0OR";
                else
                    return "V0OR";
            case 10:
                return "V0AND";
            case 11:
                return "T0AND";
            case 12:
                if (periodName.CompareTo("LHC11a")==0||periodName.CompareTo("LHC11a_p4_wSDD")==0)
                    return "V0AND";
                return ""; // To prevent fall through
            case 13:
                if (periodName.CompareTo("LHC11a")==0||periodName.CompareTo("LHC11a_p4_wSDD")==0)
                    return "with SDD, V0AND";
                else
                    return "V0AND";
            case 310 :
                return "INT7 & 0V0M";
            case 311 :
                return "INT7 & 0V0L";
            case 312 :
                return "INT7 & 0VHM";
            case 313 :
                return "INT7 & V0L7";
            case 51:
                return "EMC L0, INT1";
            case 52:
                return "EMC L0, INT7";
            case 53:
                return "EMC L0, INT8";
            case 54:
                return "DMC L0, INT1";
            case 55:
                return "DMC L0, INT7";
            case 56:
                return "DMC L0, INT8";
            case 61:
                return "PHOS L0, INT1";
            case 62:
                return "PHOS L0, INT7";
            case 63:
                return "PHOS L0, INT8";
            case 71:
                return "SPD HM, INT1";
            case 72:
                return "SPD HM, INT7";
            case 73:
                return "SPD HM, INT8";
            case 74:
                return "VZERO HM";
            case 75:
                return "VZERO & SPD HM";
            case 76:
                return "VZERO HM, SPD pile-up veto";
            case 81:
                return "EMC L1-GA, INT7";
            case 82:
                return "EMC L1-GA, INT8";
            case 83:
                return "EMC L1-G1, INT7";
            case 84:
                return "EMC L1-G1, INT8";
            case 85:
                return "EMC L1-G2, INT7";
            case 86:
                return "EMC L1-G2, INT8";
            case 87:
                return "DMC L1-GA, INT7";
            case 88:
                return "DMC L1-GA, INT8";
            case 89:
                return "DMC L1-G1, INT7";
            case 810:
                return "DMC L1-G1, INT8";
            case 811:
                return "DMC L1-G2, INT7";
            case 812:
                return "DMC L1-G2, INT8";
            case 91:
                return "EMC L1-JE, INT7";
            case 92:
                return "EMC L1-JE, INT8";
            case 93:
                return "EMC L1-J1, INT7";
            case 94:
                return "EMC L1-J1, INT8";
            case 95:
                return "EMC L1-J2, INT7";
            case 96:
                return "EMC L1-J2, INT8";
            case 97:
                return "DMC L1-J1, INT7";
            case 98:
                return "DMC L1-J1, INT8";
            case 99:
                return "DMC L1-J2, INT7";
            case 910:
                return "DMC L1-J2, INT8";
            default:
                return "special Trigger cut not defined";

        }
    }

    //************************************************************************************
    //******** Analyzes the multiplicity cut in pp events ********************************
    //************************************************************************************
    TString AnalysePPMultiplicityCut(Int_t minMult, Int_t maxMult){
        // Set Cut
        Int_t multBins[9] =  {0,   2,   5,    10,   15,  30,  50,  100,  1000 };

        TString multiplicityString = "";
        if ( (minMult == 0 && maxMult == 0) || maxMult<minMult){
            multiplicityString = "No selection";
        } else {
            multiplicityString =  Form("%d #leq TPC track mult < %d", multBins[minMult], multBins[maxMult]);
        }
        return multiplicityString;
    }

    //************************************************************************************
    //***** Analyzes the photon quality gamma cut, return correct cut label **************
    //************************************************************************************
    TString AnalyseHeaderSelection(Int_t eventHeader ){
        switch(eventHeader) {
            case 0:
                return "All headers";
            case 1:
                return "only MB header";
            case 2:
                return "specified added signal header";
            case 3:
                return "only MB header for gamma, specified headers for rest";
            default:
                return "header cut not defined";
        }
    }

    //************************************************************************************
    //******** Analyzes the cluster track matching cuts, return correct cut label ********
    //************************************************************************************
    TString AnalyseTrackMatchingCut(Int_t trackmatching, Int_t clusterType ){
        if (clusterType == 1){
            switch(trackmatching) {
                case 0:
                    return "TM disabled";
                case 1:
                    return "TM #scale[0.5]{#Delta#eta < 0.008, -0.03 < #Delta#varphi_{+,-} < 0.03} for V0s";
                case 2:
                    return "TM #scale[0.5]{#Delta#eta < 0.012, -0.05(-0.04) < #Delta#varphi_{+(-)} < 0.04(0.05)} for V0s";
                case 3:
                    return "TM #scale[0.5]{#Delta#eta < 0.016, -0.09(-0.06) < #Delta#varphi_{+(-)} < 0.06(0.09)} for V0s";
                case 4:
                    return "TM #scale[0.5]{#Delta#eta < 0.018, -0.11(-0.07) < #Delta#varphi_{+(-)} < 0.07(0.11)} for V0s";
                case 5:
                    return "TM #scale[0.5]{#Delta#eta < 0.020, -0.13(-0.08) < #Delta#varphi_{+(-)} < 0.08(0.13)} for V0s";
                case 6:
                    return "TM #scale[0.5]{#Delta#eta #leq 0.010 + (#frac{1}{#it{p}_{T} + (1/(0.03-0.010))^{1/2.5}})^{2.5}, #Delta#phi #leq 0.015 + (#frac{1}{#it{p}_{T} + (1/(0.08-0.015))^{1/2}})^{2}} for V0s";
                case 7:
                    return "TM #scale[0.5]{#Delta#eta #leq 0.010 + (#frac{1}{#it{p}_{T} + (1/(0.04-0.010))^{1/2.5}})^{2.5}, #Delta#phi #leq 0.015 + (#frac{1}{#it{p}_{T} + (1/(0.09-0.015))^{1/2}})^{2}} for V0s";
                case 8:
                    return "TM #scale[0.5]{#Delta#eta #leq 0.010 + (#frac{1}{#it{p}_{T} + (1/(0.05-0.010))^{1/2.5}})^{2.5}, #Delta#phi #leq 0.015 + (#frac{1}{#it{p}_{T} + (1/(0.10-0.015))^{1/1.75}})^{1.75}}   for V0s";
                case 9:
                    return "TM #scale[0.5]{#Delta#eta #leq 0.015 + (#frac{1}{#it{p}_{T} + (1/(0.06-0.015))^{1/2.5}})^{2.5}, #Delta#phi #leq 0.020 + (#frac{1}{#it{p}_{T} + (1/(0.12-0.020))^{1/1.75}})^{1.75}} for V0s";
                default:
                    return "track missmatch not defined ";
            }
        } else if (clusterType == 2){
            switch(trackmatching) {
                case 0:
                    return "TM disabled";
                case 1:
                    return "TrMatch excl. #Delta#eta < 0.005, -0.03 < #Delta#varphi_{+} < 0.03 for V0s";
                case 2:
                    return "TrMatch excl. #Delta#eta < 0.010, -0.09(-0.07) < #Delta#varphi_{+} < 0.07(0.09) for V0s";
                case 3:
                    return "TrMatch excl. #Delta#eta < 0.015, -0.15(-0.11) < #Delta#varphi_{+} < 0.11(0.15) for V0s";
                default:
                    return "track missmatch not defined ";
            }
        }  else {
            cout << "AnalyseTrackMatchingCut: clusterType " << clusterType << " not recognized!" << endl;
            return "track missmatch not defined ";
        }
    }

    //************************************************************************************
    //**** Analyzes the cluster track matching cuts for calo, return correct cut label ***
    //************************************************************************************
    TString AnalyseTrackMatchingCaloCut(Int_t trackmatching, Int_t clusterType ){
        if (clusterType == 1){
            switch(trackmatching) {
                case 0:
                    return "TM disabled";
                case 1:
                    return "TM #scale[0.5]{#Delta#eta < 0.008, -0.03 < #Delta#varphi_{+,-} < 0.03}";
                case 2:
                    return "TM #scale[0.5]{#Delta#eta < 0.012, -0.05(-0.04) < #Delta#varphi_{+(-)} < 0.04(0.05)}";
                case 3:
                    return "TM #scale[0.5]{#Delta#eta < 0.016, -0.09(-0.06) < #Delta#varphi_{+(-)} < 0.06(0.09)}";
                case 4:
                    return "TM #scale[0.5]{#Delta#eta < 0.018, -0.11(-0.07) < #Delta#varphi_{+(-)} < 0.07(0.11)}";
                case 5:
                    return "TM #scale[0.5]{#Delta#eta < 0.020, -0.13(-0.08) < #Delta#varphi_{+(-)} < 0.08(0.13)}";
                case 6:
                    return "TM #scale[0.5]{#Delta#eta #leq 0.010 + (#frac{1}{#it{p}_{T} + (1/(0.03-0.010))^{1/2.5}})^{2.5}, #Delta#phi #leq 0.015 + (#frac{1}{#it{p}_{T} + (1/(0.08-0.015)})^{1/2})^{2}}";
                case 7:
                    return "TM #scale[0.5]{#Delta#eta #leq 0.010 + (#frac{1}{#it{p}_{T} + (1/(0.04-0.010))^{1/2.5}})^{2.5}, #Delta#phi #leq 0.015 + (#frac{1}{#it{p}_{T} + (1/(0.09-0.015)})^{1/2})^{2}}";
                case 8:
                    return "TM #scale[0.5]{#Delta#eta #leq 0.010 + (#frac{1}{#it{p}_{T} + (1/(0.05-0.010))^{1/2.5}})^{2.5}, #Delta#phi #leq 0.015 + (#frac{1}{#it{p}_{T} + (1/(0.10-0.015)})^{1/1.75})^{1.75}}";
                case 9:
                    return "TM #scale[0.5]{#Delta#eta #leq 0.015 + (#frac{1}{#it{p}_{T} + (1/(0.06-0.015))^{1/2.5}})^{2.5}, #Delta#phi #leq 0.020 + (#frac{1}{#it{p}_{T} + (1/(0.12-0.020)})^{1/1.75})^{1.75}}";
                case 18:
                    return "TM #scale[0.5]{MIP subtraction: #Delta#eta < 0.010, -0.011 < #Delta#varphi_{+,-} < 0.011}";
                case 19:
                    return "TM #scale[0.5]{MIP subtraction: #Delta#eta < 0.015, -0.02 < #Delta#varphi_{+,-} < 0.02}";
                case 20:
                    return "TM #scale[0.5]{MIP subtraction: #Delta#eta #leq 0.010 + (#frac{1}{#it{p}_{T} + (1/(0.04-0.010))^{1/2.5}})^{2.5}, #Delta#phi #leq 0.015 + (#frac{1}{#it{p}_{T} + (1/(0.09-0.015)})^{1/2})^{2}}";
                default:
                    return "track missmatch not defined ";
            }
        } else if (clusterType == 2){
            switch(trackmatching) {
                case 0:
                    return "TM disabled";
                case 1:
                    return "TrMatch excl. #Delta#eta < 0.005, -0.03 < #Delta#varphi_{+} < 0.03";
                case 2:
                    return "TrMatch excl. #Delta#eta < 0.010, -0.09(-0.07) < #Delta#varphi_{+} < 0.07(0.09)";
                case 3:
                    return "TrMatch excl. #Delta#eta < 0.015, -0.15(-0.11) < #Delta#varphi_{+} < 0.11(0.15)";
                default:
                    return "track missmatch not defined ";
            }
        } else {
            cout << "AnalyseTrackMatchingCaloCut: clusterType " << clusterType << " not recognized!" << endl;
            return "track missmatch not defined ";
        }
    }

    //************************************************************************************
    //***************** Analyzes the min energy, return correct cut value ****************
    //************************************************************************************
    Double_t ReturnMinClusterEnergy(TString clusterCutSelection){
        TString minEnergyCutNumber(clusterCutSelection(GetClusterMinEnergyCutPosition(clusterCutSelection),1));
        Int_t minEnergyCut          = CutNumberToInteger(minEnergyCutNumber);
        switch(minEnergyCut){
            case 0:
                return 0;
                break;
            case 1:
                return 0.5;
                break;
            case 2:
                return 0.6;
                break;
            case 3:
                return 0.7;
                break;
            case 4:
                return 0.8;
                break;
            case 5:
                return 0.9;
                break;
            case 6:
                return 4.5;
                break;
            case 7:
                return 5.0;
                break;
            case 8:
                return 5.5;
                break;
            case 9:
                return 6.0;
                break;
            default:
                return 0;
        }
    }

    //************************************************************************************
    //***************** Analyzes the min energy, return correct cut value ****************
    //************************************************************************************
    Double_t ReturnMinNCells(TString clusterCutSelection){
        TString nCellsCutNumber(clusterCutSelection(GetClusterMinNCellsCutPosition(clusterCutSelection),1));
        return CutNumberToInteger(nCellsCutNumber);
    }

    //************************************************************************************
    //***************** Analyzes the N cells cut, return correct cut label ****************
    //************************************************************************************
    TString AnalyseNCellsCut (Int_t nCellsCut){
        switch(nCellsCut){
            case 0:
                return "No Ncells cut";
            case 1:
                return "n Cells #geq 1";
            case 2:
                return "n Cells #geq 2";
            case 3:
                return "n Cells #geq 3";
            case 4:
                return "n Cells #geq 4";
            case 5:
                return "n Cells #geq 5";
            case 6:
                return "n Cells #geq 6";
            default:
                return "n Cells cut not defined";
        }
    }

    //************************************************************************************
    //********************* Analyzes the M02 cut, return correct cut label ***************
    //************************************************************************************
    TString AnalyseM02Cut(Int_t minM02, Int_t maxM02){
        Double_t fMinM02 = 0.;
        cout << minM02 << "\t" << maxM02 << endl;
        switch(minM02){
            case 0:
                fMinM02=0;
                break;
            case 1:
                fMinM02=0.002;
                break;
            case 2:
                fMinM02=0.1;
                break;
            case 3:
                fMinM02=0.2;
                break;
            default:
                fMinM02 = -10;
                break;

        }

        Double_t fMaxM02 = 1000.;
        Bool_t isPtDepM02 = kFALSE;
        TString pTdepM02 = "";
        switch(maxM02){
            case 0:
                fMaxM02=1000;
                break;
            case 1:
                fMaxM02=1.;
                break;
            case 2:
                fMaxM02=0.7;
                break;
            case 3:
                fMaxM02=0.5;
                break;
            case 4:
                fMaxM02=0.4;
                break;
            case 5:
                fMaxM02=0.3;
                break;
            case 6:
                fMaxM02=0.27;
                break;
            case 7:
                fMaxM02=1.3;
                break;
            case 8:
                fMaxM02=2.5;
                break;
            case 9:
                fMaxM02=0.35;
                break;
            case 10:
                fMaxM02=0.33;
                break;
            case 11:
                fMaxM02=0.28;
                break;
            case 12:
                fMaxM02=0.32;
                break;
            case 13: // d
                fMaxM02=0.4;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 4.2: (0.27+0.0072#it{E}^{2}, #it{E} > 4.2: 0.4)";
                break;
            case 14: // e
                fMaxM02=0.5;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 5.1: (0.31+0.0072#it{E}^{2}, #it{E} > 5.1: 0.5)";
                break;
            case 15: // f
                fMaxM02=0.7;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 6.9: (0.36+0.0072#it{E}^{2}, #it{E} > 6.9: 0.7)";
                break;
            case 16: // g
                fMaxM02=0.7;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 6.8: (0.37+0.0072#it{E}^{2}, #it{E} > 6.8: 0.7)";
                break;
            case 17: // h
                fMaxM02=0.5;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 5.3: (0.30+0.0072#it{E}^{2}, #it{E} > 5.3: 0.5)";
                break;
            case 18: // i
                fMaxM02=0.7;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 7.0: (0.35+0.0072#it{E}^{2}, #it{E} > 7.0: 0.7)";
                break;
            case 19: // j
                fMaxM02=0.39;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 4.4: (0.25+0.0072#it{E}^{2}, #it{E} > 4.4: 0.39)";
                break;
            case 20: // k
                fMaxM02=0.5;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 5.0: (0.27+0.0092#it{E}^{2}, #it{E} > 5.0: 0.5)";
                break;
            case 21: // l
                fMaxM02=0.5;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 5.0: (0.32+0.0072#it{E}^{2}, #it{E} > 5.0: 0.5)";
                break;
            case 22: // m
                fMaxM02=0.5;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 3.4: (0.32+0.0152#it{E}^{2}, #it{E} > 3.4: 0.5)";
                break;
            case 23: // n
                fMaxM02=0.5;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 2.8: (0.32+0.0238#it{E}^{2}, #it{E} > 2.8: 0.5)";
                break;
            case 24: // o
                fMaxM02=0.7;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 6.8: (0.27+0.0092#it{E}^{2}, #it{E} > 6.8: 0.7)";
                break;
            case 25: // p
                fMaxM02=0.7;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 7.3: (0.32+0.0072#it{E}^{2}, #it{E} > 7.3: 0.7)";
                break;
            case 26: // q
                fMaxM02=0.7;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 7.1: (0.34+0.0072#it{E}^{2}, #it{E} > 7.1: 0.7)";
                break;
            case 27: // r
                fMaxM02=0.7;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 5.9: (0.25+0.0072#it{E}^{2}, #it{E} > 5.9: 0.5)";
                break;
            case 28: // s
                fMaxM02=0.7;
                isPtDepM02=kTRUE;
                pTdepM02="#it{E} #leq 4.0: (0.32+0.0238#it{E}^{2}, #it{E} > 4.0: 0.7)";
                break;
            default:
                fMaxM02 = -10;
                break;
        }
        if(isPtDepM02) return Form("#scale[0.7]{%1.1f < #sigma_{long}^{2} < for %s}", fMinM02, pTdepM02.Data());
        else return Form("%1.1f < #sigma_{long}^{2} < %3.2f", fMinM02, fMaxM02);

    }

    //************************************************************************************
    //********* Analyzes the M02 cut for merged ana, return correct cut label ************
    //************************************************************************************
    TString AnalyseMergedM02Cut(Int_t minM02, Int_t maxM02){
        Double_t fMinM02 = 0.;
        cout << minM02 << "\t" << maxM02 << endl;
        switch(minM02){
            case 0:
                fMinM02=0.1;
                break;
            case 1:
                fMinM02=-10;
                break;
            case 2:
                fMinM02=-10;
                break;
            case 3:
                fMinM02=-10;
                break;
            case 4:
                fMinM02=-10;
                break;
            case 5:
                fMinM02=-10;
                break;
            case 6:
                fMinM02=0.30;
                break;
            case 7:
                fMinM02=0.27;
                break;
            case 8:
                fMinM02=0.25;
                break;
            default:
                fMinM02=-10;
                break;

        }

        Double_t fMaxM02 = 1000.;
        switch(maxM02){
            case 0:
                fMaxM02=1000;
                break;
            case 1:
                fMaxM02=-10;
                break;
            case 2:
                fMaxM02=-10;
                break;
            case 3:
                fMaxM02=-10;
                break;
            default:
                fMaxM02 = -10;
                break;
        }
        if (fMinM02 != -10 && fMaxM02 != -10)
        return Form("%1.3f < M_{02} < %3.3f", fMinM02, fMaxM02);
        else if (fMinM02 != -10)
        return Form("M_{02} > %1.3f, #it{p}_{T} dependent upper limit %d ", fMinM02, maxM02);
        else if (fMinM02 != -10)
        return Form("M_{02} < %1.3f, #it{p}_{T} dependent lower limit %d ", fMaxM02, minM02);
        else
        return Form("M_{02}: #it{p}_{T} dependent lower (%d) & upper (%d) limit  ", minM02, maxM02);
    }


    //************************************************************************************
    //********************* Analyzes the phi cut, return correct cut value ***************
    //************************************************************************************
    TString AnalyseAcceptanceCutPhiCluster(Int_t minPhi, Int_t maxPhi){
        if (minPhi == 0 && maxPhi == 0){
            return "Full acceptance";
        } else if (minPhi == 1 && maxPhi ==1){
            return "restricted to EMCal acceptance";
        } else if (minPhi == 2 && maxPhi ==1){
            return "EMCal 2012/13 geometry with TRD (2.1 < #varphi < 3.15)";
        } else if (minPhi == 1 && maxPhi ==3){
            return "EMCal 2012/13 geometry no TRD (1.39626 < #varphi < 2.1)";
        } else if (minPhi == 3 && maxPhi ==1){
            return "EMCal 2011 geometry with TRD (2.45 < #varphi < 3.15)";
        } else if (minPhi == 1 && maxPhi ==2){
            return "EMCal 2011 geometry no TRD (1.39636 < #varphi < 2.45)";
        }
        return "";
    }

    //************************************************************************************
    //***************** Analyzes the min energy cut, return correct cut label ************
    //************************************************************************************
    TString AnalyseMinEnergyCut(Int_t minEnergyCut, Int_t clusterType=1){
        Double_t fMinEnergy = 0.;
        switch(minEnergyCut){
            case 0:
                fMinEnergy = 0;
                break;
            case 1:
                if(clusterType!=2){
                    fMinEnergy = 0.5;
                } else{ //PHOS
                    fMinEnergy = 0.3;
                }
                break;
            case 2:
                if(clusterType!=2){
                    fMinEnergy = 0.6;
                } else{ //PHOS
                    fMinEnergy = 0.5;
                }
                break;
            case 3:
                if(clusterType!=2){
                    fMinEnergy = 0.7;
                } else{ //PHOS
                    fMinEnergy = 0.6;
                }
                break;
            case 4:
                if(clusterType!=2){
                    fMinEnergy = 0.8;
                } else{ //PHOS
                    fMinEnergy = 0.7;
                }
                break;
            case 5:
                if(clusterType!=2){
                    fMinEnergy = 0.9;
                } else{ //PHOS
                    fMinEnergy = 0.8;
                }
                break;
                break;
            case 6:
                if(clusterType!=2){
                    fMinEnergy = 4.5;
                } else{ //PHOS
                    fMinEnergy = 0.9;
                }
                break;
            case 7:
                if(clusterType!=2){
                    fMinEnergy = 5.0;
                } else{ //PHOS
                    fMinEnergy = 0.2;
                }
                break;
            case 8:
                if(clusterType!=2){
                    fMinEnergy = 5.5;
                } else{ //PHOS
                    fMinEnergy = 0.4;
                }
                break;
            case 9:
                if(clusterType!=2){
                    fMinEnergy = 6.0;
                } else{ //PHOS
                    fMinEnergy = -1; // not defined
                }
                break;
            default:
                fMinEnergy = 0;
                break;
        }

        return Form("E_{clus} > %3.3f GeV/c", fMinEnergy);
    }

    //************************************************************************************
    //***************** Analyzes the min energy cut, return correct cut label ************
    //************************************************************************************
    TString AnalyseMinEnergyCutPHOS(Int_t minEnergyCut){
        Double_t fMinEnergy = 0.;
        switch(minEnergyCut){
          case 0:
            fMinEnergy=0.1;
            break;
          case 1:
            fMinEnergy=0.3;
            break;
          case 2:
            fMinEnergy=0.5;
            break;
          case 3:
            fMinEnergy=0.6;
            break;
          case 4:
            fMinEnergy=0.7;
            break;
          case 5:
            fMinEnergy=0.8;
            break;
          case 6:
            fMinEnergy=0.9;
            break;
          case 7:
            fMinEnergy=0.2;
            break;
          case 8:
            fMinEnergy=0.4;
            break;
          default:
            fMinEnergy = 0;
            break;
          }
          return Form("E_{clus} > %3.3f GeV/c", fMinEnergy);
    }



    //************************************************************************************
    //***************** Analyzes the min energy cut, return correct cut label ************
    //************************************************************************************
    TString AnalyseMinEnergyMergedCut(Int_t minEnergyCut){
        Double_t fMinEnergy = 0.;
        switch(minEnergyCut){
            case 0:
                fMinEnergy = 0.1;
                break;
            case 1:
                fMinEnergy = 4.0;
                break;
            case 2:
                fMinEnergy = 5.0;
                break;
            case 3:
                fMinEnergy = 6.0;
                break;
            case 4:
                fMinEnergy = 7.0;
                break;
            case 5:
                fMinEnergy = 7.5;
                break;
            case 6:
                fMinEnergy = 8.0;
                break;
            case 7:
                fMinEnergy = 8.5;
                break;
            case 8:
                fMinEnergy = 9.0;
                break;
            case 9:
                fMinEnergy = 9.5;
                break;
            default:
                fMinEnergy = 0;
                break;
        }

        return Form("E_{clus} > %3.3f GeV/c", fMinEnergy);
    }


    //************************************************************************************
    //** Analyzes the cluster timing cut, return correct cut label **************************
    //************************************************************************************
    TString AnalyseClusterTimingCut(Int_t timing){
        switch(timing){
            case 0:
                return "-500 s < t_{clus} < 500 s";
            case 1:
                return "-1 #mu s < t_{clus} < 1 #mu s";
            case 2:
                return "-500 ns < t_{clus} < 500 ns";
            case 3:
                return "-200 ns < t_{clus} < 200 ns";
            case 4:
                return "-100 ns < t_{clus} < 100 ns";
            case 5:
                return "-50 ns < t_{clus} < 50 ns";
            case 6:
                return "-30 ns < t_{clus} < 35 ns";
            case 7:
                return "-30 ns < t_{clus} < 30 ns";
            case 8:
                return "-20 ns < t_{clus} < 30 ns";
            case 9:
                return "-20 ns < t_{clus} < 25 ns";
            case 10:
                return "-12.5 ns < t_{clus} < 13 ns";
            case 11:
                return "-130 ns < t_{clus} < 130 ns";
            case 12:
                return "-110 ns < t_{clus} < 110 ns";
            case 13:
                return "-120 ns < t_{clus} < 120 ns";
            case 14:
                return "-90 ns < t_{clus} < 90 ns";
            case 15:
                return "-80 ns < t_{clus} < 80 ns";
            case 16:
                return "-30 ns < t_{clus} < 30 ns - w. TimingCutEfficiency";
            case 17:
                return "-50 ns < t_{clus} < 50 ns - w. TimingCutEfficiency";
            case 18:
                return "-12.5 ns < t_{clus} < 13 ns - w. TimingCutEfficiency";
            case 19:
                return "-30 ns < t_{clus} < 30 ns - w. TimingCutEfficiency";
            case 20:
                return "-50 ns < t_{clus} < 50 ns - w. TimingCutEfficiency";
            case 21:
                return "-12.5 ns < t_{clus} < 13 ns - w. TimingCutEfficiency";
            default:
                return "no timing cut defined";
        }
    }

    //************************************************************************************
    //** Analyzes the cluster timing cut, return correct cut label **************************
    //************************************************************************************
    TString AnalyseClusterNonLinearityCut(Int_t nonlinearity){


    //   return Form("non linearity cut %d",nonlinearity);
    switch(nonlinearity){
        case 0:
            return "none";
        case 1:
            return "SDM (Jason)";
        case 2:
            return "test-beam";
        case 11:
            return "CCRF";
        case 12:
            return "CRF";
        case 13:
            return "TBv3 w/ CCRF";
        case 14:
            return "TBv3 w/ CRF";
        case 15:
            return "ConvCalo with data shifted";
        case 16:
            return "SDM w/ data corr (new)";
        case 21:
            return "CCMF";
        case 22:
            return "CMF";
        case 31:
            return "ConvCalo, open timing";
        case 32:
            return "SDM w/o data corr (new), open timing";
        case 41:
            return "CCRF";
        case 42:
            return "CRF";
        case 51:
            return "CCMF";
        case 52:
            return "CMF";

        default:
            return Form("non linearity cut %d",nonlinearity);
    }
    //  switch(nonlinearity){
    //    case 0:
    //      return "from tender?";
    //    default:
    //      return "no nonlinearity defined";
    //  }
    }

    //************************************************************************************
    //***************************** Returns the cluster NLM ******************************
    //************************************************************************************
    Int_t ReturnClusterNLM( TString clusterCutString){
        TString nCellsCutNumber(clusterCutString(GetClusterNLMCutPosition(clusterCutString),1));
        return CutNumberToInteger(nCellsCutNumber);
    }

    //************************************************************************************
    //** Analyzes the alpha meson cut, return correct cut label **************************
    //************************************************************************************
    TString AnalyseAlphaMesonCut(Int_t alphaMesonCut){
        switch(alphaMesonCut){
            case 0:	// 0- 0.7
                return "|#alpha_{meson}| < 0.7";
            case 1:	// 0-0.5
                return "0.65*tanh(1.8*x)"; //"|#alpha_{meson}| < 0.5";
            case 2:	// 0.5-1
                return "0.5 < |#alpha_{meson}| < 1";
            case 3:	// 0.0-1
                return "|#alpha_{meson}| < 1";
            case 4:	// 0-0.65
                return "|#alpha_{meson}| < 0.65";
            case 5:	// 0-0.75
                return "|#alpha_{meson}| < 0.75";
            case 6:	// 0-0.8
                return "|#alpha_{meson}| < 0.8";
            case 7:	// 0.0-0.85
                return "|#alpha_{meson}| < 0.85";
            case 8:	// 0.0-0.6
                return "|#alpha_{meson}| < 0.6";
            case 9:	// 0.0-0.3
                return "|#alpha_{meson}| < 0.3";

    	    case 10:	// 0.0-0.2
                return "|#alpha_{meson}| < 0.2";

            case 13:	// 0.0-0.1
                return "|#alpha_{meson}| < 0.1";

            case 14:	// 0.0-0.3
                return "|#alpha_{meson}| < 0.3";

            case 15:	// 0.0-0.4
                return "|#alpha_{meson}| < 0.4";

	    case 16:	// 0.0-0.5
                return "|#alpha_{meson}| < 0.5";
            default:
                return "no alpha cut defined";
        }
    }

    //************************************************************************************
    //** Analyzes the MC smearing cut for mesons, return correct cut label ***************
    //************************************************************************************
    TString AnalyseMCSmearingCut(Int_t mcSmearingCut){
        switch(mcSmearingCut){
            case 0:
                return "no additional smearing" ;
            case 1:
                return "#sqrt{0.011^{2} + 0.007^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
            case 2:
                return "#sqrt{0.022^{2} + 0.007^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
            case 3:
                return "#sqrt{0.044^{2} + 0.007^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
            case 7:
                return "#sqrt{0.011^{2} + 0.014^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
            case 8:
                return "#sqrt{0.011^{2} + 0.0035^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
            case 9:
                return "#sqrt{0.011^{2} + 0.028^{2} #times P_{#gamma}^{2}} #times Rand.Gaus(0,1)" ;
            default:
                return "smearing cut not defined";
        }
    }

    //************************************************************************************
    //********* Analyzes the meson BG scheme, return correct cut label *******************
    //************************************************************************************
    TString AnalyseBackgroundScheme(TString sBackgroundScheme){
        TString sBackgroundSchemeB      = sBackgroundScheme(0,1);
        Int_t BackgroundScheme          = CutNumberToInteger(sBackgroundSchemeB);
        TString bGScheme                = "";
        cout << "BackgroundScheme: " << BackgroundScheme << endl;

        switch(BackgroundScheme){
            case 0: //Rotation
                bGScheme="Rotation";
                break;
            case 1: // mixed event with V0 multiplicity
                bGScheme="mixed event, V0 mult";
                break;
            case 2: // mixed event with track multiplicity
                bGScheme="mixed event, track mult";
                break;
            case 3: //Rotation
                bGScheme="Rotation, with prob.";
                break;
            case 4: //No BG calculation
                bGScheme="no BG calculated";
                break;
            case 5: //Rotation
                bGScheme="Rotation, new BG handler";
                break;
            case 6: // mixed event with V0 multiplicity
                bGScheme="mixed event, new BG handler, V0 mult";
                break;
            case 7: // mixed event with track multiplicity
                bGScheme="mixed event, new BG handler, track mult";
                break;
            case 8: //Rotation
                bGScheme="Rotation, new BG handler, with prob. ";
                break;
            case 9: // mixed event with PtMax method
                bGScheme="mixed event, Pt max method ";
                break;
            case 10: // likesign mixing
                bGScheme="likesign mixing";
                break;
            case 11: // sideband mixing right side
                bGScheme="sideband mixing (180-220 MeV)";
                break;
            case 12: // sideband mixing left side
                bGScheme="sideband mixing (10-50 MeV)";
                break;
            case 13: // sideband mixing both sides
                bGScheme="sideband mixing both sides";
                break;
            default:
                bGScheme="no BG method selected";

        }

        TString sNumberOfBGEvents       = sBackgroundScheme(1,1);
        Int_t NumberOfBGEvents          = CutNumberToInteger(sNumberOfBGEvents);
        Int_t fNumberOfBGEvents         = 0;
        cout << "NumberOfBGEvents: " << NumberOfBGEvents << endl;

        switch(NumberOfBGEvents){
            case 0:
                fNumberOfBGEvents = 5;
                break;
            case 1:
                fNumberOfBGEvents = 10;
                break;
            case 2:
                fNumberOfBGEvents = 15;
                break;
            case 3:
                fNumberOfBGEvents = 20;
                break;
            case 4:
                fNumberOfBGEvents = 2;
                break;
            case 5:
                fNumberOfBGEvents = 50;
                break;
            case 6:
                fNumberOfBGEvents = 80;
                break;
            case 7:
                fNumberOfBGEvents = 100;
                break;
            default:
                cout<<"Warning: NumberOfBGEvents not defined "<<NumberOfBGEvents<<endl;
        }

        TString sDegreesForRotationMethod   = sBackgroundScheme(1,2);
        Int_t DegreesForRotationMethod      = CutNumberToInteger(sDegreesForRotationMethod);
        Int_t fnDegreeRotationPMForBG       = 0;
        cout << "DegreesForRotationMethod: " << DegreesForRotationMethod << endl;

        switch(DegreesForRotationMethod){
            case 0:
                fnDegreeRotationPMForBG = 5;
                break;
            case 1:
                fnDegreeRotationPMForBG = 10;
                break;
            case 2:
                fnDegreeRotationPMForBG = 15;
                break;
            case 3:
                fnDegreeRotationPMForBG = 20;
                break;
            default:
                cout<<"Warning: DegreesForRotationMethod not defined "<<DegreesForRotationMethod<<endl;
        }

        if (BackgroundScheme == 0 || BackgroundScheme == 3 || BackgroundScheme == 5 || BackgroundScheme == 8){
            return Form("%s within #pm %i°, %i photons per pool",bGScheme.Data(), fnDegreeRotationPMForBG, fNumberOfBGEvents)   ;
        } else {
            return Form("%s with %i photons per pool",bGScheme.Data(), fNumberOfBGEvents)   ;
        }
        return "";
    }

    //************************************************************************************
    //*********** Analyzes the rapidity meson cut, return correct cut label **************
    //************************************************************************************
    TString AnalyseRapidityMesonCut(Int_t RapidityMesonCut){
        // Set Cut
        switch(RapidityMesonCut){
            case 0:  //
                return "|y_{meson}| < 1.35";
            case 1:  //
                return "|y_{meson}| < 0.8";
            case 2:  //
                return "|y_{meson}| < 0.7";
            case 3:  //
                return "|y_{meson}| < 0.6";
            case 4:  //
                return "|y_{meson}| < 0.5";
            case 5:  //
                return "|y_{meson}| < 0.85";
            case 6:  //
                return "|y_{meson}| < 0.75";
            case 7:  //
                return "|y_{meson}| < 0.3";
            case 8:  //
                return "|y_{meson}| < 0.25";
            case 9:  //
                return "|y_{meson}| < 0.4";
            default:
                return "rapidity cut not defined";
        }
    }

    //************************************************************************************
    //*********** Analyzes the rapidity meson cut for pPv, return correct cut label ******
    //************************************************************************************
    TString AnalyseRapidityMesonCutpPb(Int_t RapidityMesonCut){
        // Set Cut
        switch(RapidityMesonCut){
            case 0:  //
                return "|y_{meson}| < 1.35";
            case 1:  //
                return "|y_{meson}| < 0.8";
            case 2:  //
                return "|y_{meson}| < 0.7";
            case 3:  //
                return "|y_{meson}| < 0.6";
            case 4:  //
                return "|y_{meson}| < 0.5";
            case 5:  //
                return "|y_{meson}| < 0.85";
            case 6:  //
                return "|y_{meson}| < 0.75";
            case 7:  //
                return "|y_{meson}| < 0.3";
            case 8:  //
                return "|y_{meson}| < 0.25";
            case 9:  //
                return "|y_{meson}| < 0.4";
            default:
                return "rapidity cut not defined";
        }
    }

    //************************************************************************************
    //******** Analyzes the chi2 meson cut, return correct cut label *********************
    //************************************************************************************
    TString AnalyseChi2MesonCut(Int_t chi2GammaCut){   // Set Cut
        switch(chi2GammaCut){
            case 0: // 100
                return "#chi_{meson}^{2} < 100";
            case 1:  // 50
                return "#chi_{meson}^{2} < 50";
            case 2:  // 30
                return "#chi_{meson}^{2} < 30";
            case 3:
                return "#chi_{meson}^{2} < 200";
            case 4:
                return "#chi_{meson}^{2} < 500";
            case 5:
                return "#chi_{meson}^{2} < 1000";
            default:
                return "#chi_{meson}^{2} cut unknown";
        }
    }

    //************************************************************************************
    //******** Analyzes the openign angle meson cut, return correct cut label *********************
    //************************************************************************************
    TString AnalyseMesonOpeningAngleCut(Int_t OpeningAngleCut){   // Set Cut
    switch(OpeningAngleCut){
        case 0:
        return "#theta_{meson} > 0";
        case 1:
        return "#theta_{meson} > 0.005";
        case 2:
        return "#theta_{meson} > f(#it{p}_{T})";
        case 3:
        return "#theta_{meson} > 0.01";
        case 4:
        return "#theta_{meson} > 0.0152";
        case 5:
        return "#theta_{meson} > 0.0202";
        case 6:
        return "#theta_{meson} > 0.017";
        case 7:
        return "#theta_{meson} > 0.016";
        case 8:
        return "#theta_{meson} > 0.018";
        case 9:
        return "#theta_{meson} > 0.019";
        case 10:
        return "#theta_{meson} > 0. + 1 cell dist";
        case 11:
        return "#theta_{meson} > 0.0152 + 1 cell dist";
        case 12:
        return "#theta_{meson} > 0.016 + 1 cell dist";
        case 13:
        return "#theta_{meson} > 0.017 + 1 cell dist";
        case 14:
        return "#theta_{meson} > 0.018 + 1 cell dist";
        case 15:
        return "#theta_{meson} > 0.019 + 1 cell dist";
        case 16:
        return "#theta_{meson} > 0.020 + 1 cell dist";
        default:
        return "#theta_{meson} cut unknown";
    }
    }

    //************************************************************************************
    //******** Analyze neutral meson mass selection window cut *********************
    //************************************************************************************
    TString AnalyseMesonSelectionWindowCut(Int_t SelectionWindowCut){   // Set Cut
        switch(SelectionWindowCut){
            case 0:
        return "0 MeV < M_{#gamma#gamma} < 4000 MeV";
        case 1:
        return "100 MeV < M_{#gamma#gamma} < 145 MeV";
        case 2:
        return "110 MeV < M_{#gamma#gamma} < 145 MeV";
        case 3:
        return "120 MeV < M_{#gamma#gamma} < 145 MeV";
        case 4:
        return "100 MeV < M_{#gamma#gamma} < 150 MeV";
        case 5:
        return "110 MeV < M_{#gamma#gamma} < 150 MeV";
        case 6:
        return "120 MeV < M_{#gamma#gamma} < 150 MeV";
        case 7:
        return "100 MeV < M_{#gamma#gamma} < 155 MeV";
        case 8:
        return "125 MeV < M_{#gamma#gamma} < 145 MeV";
        case 9:
        return "110 MeV < M_{#gamma#gamma} < 155 MeV";
        case 10:
        return "80 MeV < M_{#gamma#gamma} < 145 MeV";
        default:
        return "M_{#gamma#gamma} cut unknown";
    }
    }

    //************************************************************************************
    //** Analyzes the RBins cut, return correct cut label **************************
    //************************************************************************************
    TString AnalyseRBinCut(Int_t RBinCut){
        switch(RBinCut){
            case 0:	// 5-33.5
                return "0 cm < R < 180. cm ";
            case 2:	// 5- 180
                return " 5 cm < R < 180 cm ";
            case 10:	//a 5-33.5
                return "5 cm < R < 33.5 cm ";
            case 11:	//b 33.5-72.
                return "33.5 cm < R < 72 cm ";
            case 12:	//c 72-180
                return "72 cm < R < 180 cm ";
            case 16:	//g 95-180
                return "95 cm < R < 180 cm ";
            case 17:	//h 5-13
                return "5 cm < R < 13. cm ";
            case 18:	//i 13-33.5
                return "13 cm < R < 33.5 cm ";
            case 19:	//i 33.5-55
                return "33.5 cm < R < 55 cm ";
            case 20:	//i 55-72
                return "55 cm < R < 72 cm ";
            case 21:	//i 72-95
                return "72 cm < R < 95 cm ";
            case 22:	//i 5-180, except 55-72
                return "5 cm < R < 55 cm OR 72 cm < R < 180 cm ";
            default:
                return "RBins cut unknown";
        }
    }
    //******** Analyze neutral meson pt cut *********************
    //************************************************************************************
    TString AnalyseMesonPtCut(Int_t pTCut){   // Set Cut
        switch(pTCut){
        case 0:
        return "no cut";
        case 1:
        return "p_{T} > 0.4 GeV/c";
        case 2:
        return "p_{T} > 0.7 GeV/c";
        case 3:
        return "p_{T} > 0.9 GeV/c";
        case 4:
        return "p_{T} > 1.0 GeV/c";
        case 5:
        return "p_{T} > 1.2 GeV/c";
        case 6:
        return "p_{T} > 1.5 GeV/c";
        case 7:
        return "p_{T} > 0.5 GeV/c";

        default:
        return "cut unknown";
    }
    }
    //************************************************************************************
    //******** Analyze Charged Pion Cut ClsTPC *********************
    //************************************************************************************
    TString AnalysePionClsTPCCut(Int_t clsTPCCut){   // Set Cut
        switch(clsTPCCut){
            case 0:
            return "N_{Cls} > 0";
            case 1:
            return "N_{Cls} > 70";
            case 2:
            return "N_{Cls} > 80";
            case 3:
            return "N_{Cls} > 100";
            case 4:
            return "N_{Cls} > 70 and 0 percent of findable clusters";
            case 5:
            return "N_{Cls} > 70 and 35 percent of findable clusters";
            case 6:
            return "N_{Cls} > 70 and 60 percent of findable clusters";
            case 7:
            return "N_{Cls} > 70 and 35 percent of findable clusters";
            case 8:
            return "N_{Cls} > 0 and 35 percent of findable clusters";
            case 9:
            return "N_{Cls} > 70 and 35 percent of findable clusters (use corrected info)";
            default:
            return "cut unknown";
         }
    }

    //************************************************************************************
    //******** Analyze Charged Pion DCA Cut                          *********************
    //************************************************************************************
    TString AnalysePionDCACut(Int_t DCACut){   // Set Cut
        switch(DCACut){
            case 0:
            return "open";
            case 1:
            return "DCA_{XY} < 0.0182+0.0350/pt^1.01";
            case 2:
            return "DCA_{XY} < 1 cm and max DCA_{Z} < 2 cm";
            case 3:
            return "DCA_{XY} < 0.0182+0.0350/pt^1.01 cm and max DCA_{Z} < 3 cm";
            case 4:
            return "DCA_{XY} < 0.5 cm and max DCA_{Z} < 3 cm";
            default:
            return "cut unknown";
         }
    }

    //************************************************************************************
    //******** Analyze Charged Pion Min Pt ***********************************************
    //************************************************************************************
    TString AnalysePionPtCut(Int_t minPtCut){   // Set Cut
        switch(minPtCut){
            case 0:
            return "p_{T} > 0.075 GeV/c";
            case 1:
            return "p_{T} > 0.1 GeV/c";
            case 2:
            return "p_{T} > 0.125 GeV/c";
            case 3:
            return "p_{T} > 0.150 GeV/c";
            case 4:
            return "p_{T} > 0.400 GeV/c";
            default:
            return "cut unknown";
         }
    }

    //************************************************************************************
    //******** Analyze Charged Pion Cut dEdx TPC******************************************
    //************************************************************************************
    TString AnalysePiondEdxTPCCut(Int_t dEdxTPCCut){   // Set Cut
        switch(dEdxTPCCut){
            case 0:
            return "no cut";
            case 1:
            return "-10 < n #sigma < 10";
            case 2:
            return "-6 < n #sigma < 7";
            case 3:
            return "-5 < n #sigma < 5";
            case 4:
            return "-4 < n #sigma < 5";
            case 5:
            return "-4 < n #sigma < 4";
            case 6:
            return "-3 < n #sigma < 4";
            case 7:
            return "-3 < n #sigma < 3";
            case 8:
            return "-2 < n #sigma < 3";
            default:
            return "cut unknown";
         }
    }

    //************************************************************************************
    //******** Analyze Charged Pion Cut dEdx TPC******************************************
    //************************************************************************************
    TString AnalysePionMassCut(Int_t PionMassCut){   // Set Cut
        switch(PionMassCut){
            case 0:
            return "no cut";
            case 1:
            return "N_{#pi #pi} < 1 GeV/c^{2}";
            case 2:
            return "N_{#pi #pi} < 0.75 GeV/c^{2}";
            case 3:
            return "N_{#pi #pi} < 0.6 GeV/c^{2}";
            case 4:
            return "N_{#pi #pi} < 0.547853 GeV/c^{2}";
            case 5:
            return "N_{#pi #pi} < 0.5 GeV/c^{2}";
            case 6:
            return "N_{#pi #pi} < 0.65 GeV/c^{2}";
            case 7:
            return "N_{#pi #pi} < 0.7 GeV/c^{2}";
            case 8:
            return "N_{#pi #pi} < 0.85 GeV/c^{2}";
            case 9:
            return "N_{#pi #pi} < 1.5 GeV/c^{2}";
            default:
            return "cut unknown";
         }
    }
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    // ****************************************************************************************************************
    Double_t ReturnCorrectXSection (    TString energy,
                                        Int_t selTrig
                                   ){

        Double_t xSectionInt = 0;
        if(energy.CompareTo("900GeV") == 0){
            if (selTrig == 1){
                xSectionInt = xSection900GeVV0AND;
                cout << "V0AND xSection taken: \t" << xSectionInt << endl;
            } else if (selTrig == 3){
                xSectionInt = xSection900GeVINEL;
                cout << "INEL xSection taken: \t" << xSectionInt << endl;
            } else {
                xSectionInt = xSection900GeV;
                cout << "V0OR xSection taken: \t" << xSectionInt << endl;
            }
        } else if(energy.CompareTo("2.76TeV") == 0){
            if (selTrig == 1){
                xSectionInt = xSection2760GeVV0AND;
                cout << "V0AND xSection taken: \t" << xSectionInt << endl;
            } else if (selTrig == 3){
                xSectionInt = xSection2760GeVINEL*1e-12;
                cout << "INEL xSection taken: \t" << xSectionInt << endl;
            } else {
                xSectionInt = xSection2760GeV;
                cout << "V0OR xSection taken: \t" << xSectionInt << endl;
            }
        } else if( energy.CompareTo("5TeV") == 0 || energy.CompareTo("5.02TeV") == 0 || energy.CompareTo("5.023TeV") == 0 || energy.CompareTo("5023GeV") == 0 || energy.Contains("5TeV2017") || energy.CompareTo("5TeVSpecial") == 0){
            if (selTrig == 1){
                xSectionInt = xSection5023GeVV0AND;
                cout << "V0AND xSection taken: \t" << xSectionInt << endl;
            } else if (selTrig == 3){
                xSectionInt = xSection5023GeVINEL;
                cout << "INEL xSection taken: \t" << xSectionInt << endl;
            }
        } else if(energy.CompareTo("7TeV") == 0){
            if (selTrig == 1){
                xSectionInt = xSection7TeVV0AND;
                cout << "V0AND xSection taken: \t" << xSectionInt << endl;
            } else if (selTrig == 3){
                xSectionInt = xSection7TeVINEL;
                cout << "INEL xSection taken: \t" << xSectionInt << endl;
            } else {
                xSectionInt = xSection7TeV;
                cout << "V0OR xSection taken: \t" << xSectionInt << endl;
            }
        } else if(energy.CompareTo("8TeV") == 0){
            if (selTrig == 1){
                xSectionInt = xSection8TeVV0AND;
                cout << "V0AND xSection taken: \t" << xSectionInt << endl;
            } else if (selTrig == 2){
                xSectionInt = xSection8TeVT0AND;
                cout << "T0AND xSection taken: \t" << xSectionInt << endl;
            } else if (selTrig == 3){
                xSectionInt = xSection8TeVINEL;
                cout << "INEL xSection taken: \t" << xSectionInt << endl;
            } else {
                xSectionInt = 0.;
                cout << "ERROR: V0OR xSection not determined, set to \t" << xSectionInt << endl;
            }
        } else if(energy.CompareTo("13TeV") == 0 || energy.CompareTo("13TeVLowB") == 0 || energy.CompareTo("13TeVRBins") == 0 ){
            if (selTrig == 1){
                xSectionInt = xSection13TeVV0AND;
                cout << "V0AND xSection taken: \t" << xSectionInt << endl;
            } else if (selTrig == 3){
                xSectionInt = xSection13TeVINEL;
                cout << "INEL xSection taken: \t" << xSectionInt << endl;
            } else {
                cout << "ERROR: V0OR xSection not deterimined, set to \t" << xSectionInt << endl;
            }

        } else if( energy.CompareTo("pPb_5TeV") == 0  || energy.CompareTo("pPb_5.023TeV") == 0 || energy.CompareTo("pPb_5.023TeVCent") == 0 || energy.CompareTo("pPb_5.023TeVRun2") == 0 || energy.CompareTo("pPb_5.02TeV") == 0 ){
            if (selTrig == 3){
                xSectionInt = xSection5023GeVINELpPb;
                cout << "INEL xSection taken: \t" << xSectionInt << endl;
            }
        } else {
            cout << "ERROR: energy not deterimined, xsection set to \t" << xSectionInt << endl;
        }
        return xSectionInt;
    }


    //************************************************************************************
    //* Decodes from the mode the respective reco process and return correct label + details
    //************************************************************************************
    TString ReturnFullTextReconstructionProcess( Int_t mode, Int_t separate = 0, TString meson = "", TString clusterCutNumber = "" ){
        if( mode >= 100 ) mode -= 100;
        if (separate == 0){
            switch (mode){
                case 0: case 40: case 60:
                    return "#gamma's rec. with PCM";
                case 1: case 47: case 67:
                    return "#gamma's rec. with PCM, Dalitz";
                case 2: case 41: case 61:
                    return "#gamma's rec. with PCM, EMCal";
                case 3: case 42: case 62:
                    return "#gamma's rec. with PCM, PHOS";
                case 4: case 44: case 64:
                    return "#gamma's rec. with EMCal";
                case 5: case 45: case -5: case 65:
                    return "#gamma's rec. with PHOS";
                case 6: case 48: case 68:
                    return "#gamma's rec. with EMCal, Dalitz";
                case 7: case 49: case 69:
                    return "#gamma's rec. with PHOS, Dalitz";
                case 10:
                    if (clusterCutNumber.CompareTo("") != 0){
                        Int_t nlm = ReturnClusterNLM(clusterCutNumber);
                        if (meson.CompareTo("") == 0) return Form("rec. w/ mEMCal, %d lm", nlm );
                        else return Form("%s rec. mEMCal, %d lm", meson.Data(), nlm);
                    } else {
                        if (meson.CompareTo("") == 0) return "rec. w/ mEMCal";
                        else return Form("%s rec. mEMCal", meson.Data());
                    }
                case 11:
                    if (clusterCutNumber.CompareTo("") != 0){
                        Int_t nlm = ReturnClusterNLM(clusterCutNumber);
                        if (meson.CompareTo("") == 0) return Form("rec. w/ mPHOS, %d lm", nlm );
                        else return Form("%s rec. w/ mPHOS, %d lm", meson.Data(), nlm);
                    } else {
                        if (meson.CompareTo("") == 0) return "rec. w/ mPHOS";
                        else return Form("%s rec. w/ mPHOS", meson.Data());
                    }
                case 12: case 46:
                    return "#gamma's rec. with DCal";
                case 13: case 43:
                    return "#gamma's rec. with PCM, DCal";
                case 14:
                    return "#gamma's rec. with PCM, EMCal + DCal";
                case 15:
                    return "#gamma's rec. with DCal + EMCal";
                case 20: case 21: case 22: case 23:
                    return "combined";
                default:
                    return "not known";
            }
        } else if (separate == 1){
            switch (mode){
                case 0:
                case 1:
                case 2:
                case 3:
                    return "#gamma's rec. with PCM";
                case 6:
                case 7:
                    return "#gamma*'s rec. with Dalitz";
                case 12:
                    return "#gamma's rec. with DCAL";
                case 13:
                    return "#gamma's rec. with PCM, DCal";
                case 14:
                    return "#gamma's rec. with PCM, EMCAL/DCAL";
                case 15:
                    return "#gamma's rec. with EMCAL/DCAL";
                default:
                    return "not known";
            }
        } else if (separate == 2){
            switch (mode){
                case 1:
                    return "#gamma*'s rec. with Dalitz";
                case 0:
                case 2:
                case 4:
                case 6:
                    return "#gamma's rec. with EMCal";
                case -5:
                case 3:
                case 5:
                case 7:
                    return "#gamma's rec. with PHOS";
                case 12:
                    return "#gamma's rec. with DCAL";
                case 13:
                    return "#gamma's rec. with PCM, DCal";
                case 14:
                    return "#gamma's rec. with PCM, EMCAL/DCAL";
                case 15:
                    return "#gamma's rec. with EMCAL/DCAL";
                default:
                    return "not known";
            }
        } else {
            cout << "ReturnFullTextReconstructionProcess: separate" << separate << " not recognized!" << endl;
            return "not known";
        }
    }

    //************************************************************************************
    //******* Return correct generator name from MC name *********************************
    //************************************************************************************
    TString ReturnGeneratorNameFromMCName(TString MCname){
        if (MCname.CompareTo("LHC12f1a") == 0){
            return "Pythia 8.1, Tune 4C";
        } else if (MCname.CompareTo("LHC12f1b") == 0){
            return "Phojet";
        } else if (MCname.Contains("LHC14j4")){
            return "Pythia 6, Perugia-2011";
        } else if (MCname.CompareTo("LHC15g2") == 0){
            return "Pythia 8.2, Tune 4C";
        } else if (MCname.CompareTo("LHC15g1a") == 0 || MCname.CompareTo("LHC15a3a") == 0 ){
            return "Pythia 6, Perugia 0, p_{T,hard}";
        } else if (MCname.CompareTo("LHC15h") == 0 ){
            return "Pythia 8.2";
        } else if (MCname.CompareTo("LHC16c2") == 0 ){
            return "Pythia 8.2, p_{T,hard}";
        } else if (MCname.CompareTo("LHC13b2_efix") == 0 ){
            return "DPMJet";
        } else if (MCname.CompareTo("LHC13e7") == 0 ){
            return "HIJING";
        } else if (MCname.CompareTo("LHC17f2a") == 0 ){
            return "EPOSLHC";
        } else if (MCname.CompareTo("LHC17f2b") == 0 || MCname.CompareTo("LHC18f3") == 0 ){
            return "DPMJet";
        } else if (MCname.Contains("LHC173l") == 0 ){
            return "Pythia 8.2";
        } else if (MCname.CompareTo("LHC16P1Pyt8") == 0 || MCname.CompareTo("LHC17P1Pyt8") == 0 || MCname.CompareTo("LHC18P1Pyt8") == 0 ){
            return "Pythia 8.2, Monash2013";
        } else {
            return "undefined";
        }
        return "";
    }

    //************************************************************************************
    //** Return correct trigger rejection factor based on trigger cutnumber and energy ***
    //************************************************************************************
    Double_t ReturnTriggerRejectionFactor(TString energy, Int_t trigger, TString strTrigger = ""){
        Double_t triggerRejec   = 1;
        if (energy.CompareTo("2.76TeV") == 0){
            cout << "Trigger used: " << trigger << endl;
            if (trigger == 51){         // EMC1
                triggerRejec    = 1228;
            } else if (trigger == 52){  // EMC7
                triggerRejec    = 125;
            } else if (trigger == 85){  // EG2
                triggerRejec    = 1909;
            } else if (trigger == 83){  // EG1
                triggerRejec    = 7217;
            }
        } else  if (energy.Contains("5TeV2017") ){
            cout << "Trigger used: " << strTrigger.Data() << endl;
            if (!strTrigger.CompareTo("a1")){  // EMC7
                triggerRejec    = 2464;
            } else if (!strTrigger.CompareTo("a2")){  // EG2
                triggerRejec    = 909;
            }
        } else  if (energy.CompareTo("8TeV") == 0){
            cout << "Trigger used: " << trigger << endl;
            if (trigger == 52){  // EMC7
                triggerRejec    = 67.3;
            } else if (trigger == 81){  // EGA
                triggerRejec    = 15075;
            }
        }else  if (energy.Contains("pPb_8TeV") ){
        cout << "Trigger used: " << trigger << endl;
            if (!strTrigger.CompareTo("8e")){  // EG2
                triggerRejec    = 276;
            } else if (!strTrigger.CompareTo("8d")){  // EG1
                triggerRejec    = 874;
            } else if (!strTrigger.CompareTo("9c")){  // EJ2
                triggerRejec    = 1;
            } else if (!strTrigger.CompareTo("9b")){  // EJ1
                triggerRejec    = 1;
            }
    }
        return triggerRejec;
    }

    //************************************************************************************
    //** Return mean reconstructed R of meson reconstruction for given mode **************
    //************************************************************************************
    Double_t ReturnMeanR(Int_t mode){
        Double_t meanR   = 1.;
        if(mode == 0 || (mode == 2 || mode == 13) || mode == 3  || mode == 14){
        meanR = 60.;
        }else if(mode == 4 || mode == 12 || mode == 10  || mode == 15){
        meanR = 428.;
        }else if(mode == 5 || mode == 11 || mode == -5){
        meanR = 460.;
        } else meanR = 60.;

        return meanR;
    }

    //************************************************************************************
    //************ Return cocktail normalization factor **********************************
    //************************************************************************************
    Double_t ReturnCocktailNormalization(TString energy, TString eventCutString) {

        // cocktail is normalized per INEL event, except for: PbPb and pPb 5TeV

        Int_t       selTrig             = 0;    // 0 = V0OR, 1 = V0AND, 2 = T0OR, 3 = INEL
        TString     trigger             = eventCutString(GetEventSelectSpecialTriggerCutPosition(),2);
        if (trigger.Atoi() == 10 || trigger.Atoi() == 52 || trigger.Atoi() == 83  || trigger.Atoi() == 85 || trigger.Atoi() == 81 ) {
            selTrig                     = 1;
        }

        Double_t    xSec                = ReturnCorrectXSection(energy, selTrig);
        Double_t    xSecINEL            = ReturnCorrectXSection(energy, 3);
        Double_t    triggerRejection    = ReturnTriggerRejectionFactor(energy, trigger.Atoi(),trigger);
        Double_t    scaleFactor         = 1.;
        if (energy.BeginsWith("PbPb") || energy.BeginsWith("XeXe")) {
            scaleFactor                 = 1.;
        } else if (energy.BeginsWith("pPb")) {
            scaleFactor                 = 1.;
        } else {
            if (xSec && xSecINEL)
                scaleFactor         = xSecINEL/xSec;
            else
                cout << "ERROR: xSec not found for " << energy.Data() << ", " << selTrig << endl;
        }

        if (triggerRejection != 1.) {
            cout << "Trigger rejection factor is " << triggerRejection << endl;
            scaleFactor                 = scaleFactor * triggerRejection;
        }

        cout << "Additional normalization factor for spectra from cocktail simulation: " << scaleFactor << endl;

        return scaleFactor;
    }

    //************************************************************************************
    //************ Return correct trigger name based on trigger cutnumber ****************
    //************************************************************************************
    TString ReturnTriggerName(Int_t trigger, TString energy = "", TString strTrigger = ""){
        cout << trigger << endl;
        if (trigger == 1 || trigger == 3 || trigger == 0){  // INT1
            cout << "energy: " << energy << endl;
            if(energy.Contains("p-Pb")) return "MB";
            return "INT1";
        } else if (trigger == 10 ){                         // INT7
            return "INT7";
        } else if (trigger == 11){                          // INT7
            return "INT8";
        } else if (trigger == 51){                          // EMC1
            return "EMC1";
        } else if (trigger == 52){                          // EMC7
            return "EMC7";
        } else if (trigger == 62){                          // PHOS L0, INT7
            return "PHI7";
        } else if (trigger == 71){                          // SHM1
            return "SHM1";
        } else if (trigger == 72){                          // SHM1
            return "SHM7";
        } else if (trigger == 73){                          // SHM1
            return "SHM8";
        } else if (trigger == 74){                          // SHM1
            return "VHM";
        } else if (trigger == 75){                          // SHM1
            return "VHMSH2";
        } else if (trigger == 76){                          // SHM1
            return "VHM-SPD2";
        } else if (trigger == 80){                          // EGA
            return "EGA";
        } else if (trigger == 81){                          // EGA
            return "EGA";
        } else if (trigger == 85 || !strTrigger.CompareTo("8e")){                          // EG2
            return "EG2";
        } else if (trigger == 83 || !strTrigger.CompareTo("8d")){                          // EG1
            return "EG1";
        } else if (trigger == 95 || !strTrigger.CompareTo("9c")){                          // EG2
            return "EG2";
        } else if (trigger == 93 || !strTrigger.CompareTo("9b")){                          // EG1
            return "EG1";
        }
        return "";
    }

    //************************************************************************************
    //******************* Analyse kappa cut **********************************************
    //************************************************************************************
    Bool_t AnalyseKappaCut(Int_t cutValue, Double_t kappaMin, Double_t kappaMax){
        switch (cutValue) {
            case 0: // completely open
                kappaMin = -200;
                kappaMax = 200;
                break;
            case 1: // mainly pi pi
                kappaMin = -20;
                kappaMax = -13;
                break;
            case 2: // mainly pi e
                kappaMin = -11;
                kappaMax = -6;
                break;
            case 3: // signal
                kappaMin = -3;
                kappaMax = 5;
                break;
            case 4: // remaining
                kappaMin = 11;
                kappaMax = 20;
                break;
            case 5: // -5-10 full signal peak(including background)
                kappaMin = -5;
                kappaMax = 10;
                break;
            default:
                cout << "wrong kappa cut value" << endl;
                return kFALSE;
        }

        return kTRUE;
    }
    void CheckBinSize(Int_t &iActualBinSize, Int_t iMaximalBinSize, Bool_t giveCout=kTRUE){
        if (iActualBinSize > iMaximalBinSize) {
            if (giveCout){
                cout <<"You have chosen to have more than "<<iMaximalBinSize<<" bins, this is not possible, it will be reduced to "<<iMaximalBinSize<<"" << endl;}
            iActualBinSize    = iMaximalBinSize;}
    }

    template<Int_t ArraySize>
    Bool_t BinClamper(Double_t (&Array)[ArraySize], Int_t& PreviouslyUsedBins, Bool_t giveCout=kTRUE) {
      if (ArraySize < PreviouslyUsedBins){
          if (giveCout){
              cout <<"You have chosen to have more than "<<ArraySize<<" bins ("<<PreviouslyUsedBins<<"), this is not possible, it will be reduced to "<<ArraySize<<"" << endl;
          }
          PreviouslyUsedBins = ArraySize;
          //std::terminate(1);
          return kFALSE;
      }
      return kTRUE;
    }

    //************************************************************************************
    // Compute single cut number as int
    //************************************************************************************
    Int_t ReturnSingleAlphaNumericCutAsInt (TString cut){
        cut.ToLower();
        const char *cutChar = cut.Data();

        Int_t cutInt    = ((int)cutChar[0]>=(int)'a') ? cutChar[0]-'a'+10 : cutChar[0]-'0';
        return cutInt;
    }

#endif
