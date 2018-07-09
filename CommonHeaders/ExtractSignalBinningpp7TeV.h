#ifndef GAMMACONV_ExtractSignalBinningpp7TeV
#define GAMMACONV_ExtractSignalBinningpp7TeV

    //****************************************************************************************************
    //******************** Pt binning for pp, 7 TeV ******************************************************
    //****************************************************************************************************
    Double_t fBinsPi07TeVPt[39]                     = { 0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
                                                        1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0,
                                                        2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8,
                                                        4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0,
                                                        10.0,12.0,16.0,20.0,25.0};
    Double_t fBinsPi07TeVPCMPHOSPt[44]              = { 0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
                                                        1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1,
                                                        2.2, 2.3, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8,
                                                        4.0, 4.3, 4.6, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0,
                                                        12.0, 16.0, 20.0, 25.0};
    Double_t fBinsPi07TeVPtDCA[28]                  = { 0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
                                                        1.2, 1.3, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6,
                                                        4.0, 5.0, 6.0, 8.0, 12.0, 16.0, 20.0, 25.0};
    Double_t fBinsPi07TeVDalitzPt[23]               =  {0, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
                                                        2.2, 2.4, 2.6, 2.8, 3.2, 3.6, 4.0, 4.5, 5.0, 6.0,
                                                        8.0, 10.0, 15.};
    Int_t fBinsPi07TeVDalitzPtRebin[22]             =  {5, 5, 5, 5, 4, 4, 4, 4, 4, 4,
                                                        4, 5, 5, 5, 5,  5, 5, 5, 8, 8,
                                                        8, 10};
    Int_t fBinsPi07TeVPtRebin[38]                   = { 3, 2, 1, 1, 1, 1, 1, 1, 1,
                                                        1, 1, 1, 1, 1, 1, 1, 1, 1,
                                                        1, 1, 1, 1, 1, 1, 1, 1, 1,
                                                        1, 2, 3, 3, 4, 4, 4,
                                                        4, 4, 4, 5};
    Int_t fBinsPi07TeVPCMPHOSPtRebin[43]            = { 3, 2, 2, 2, 2, 2, 2, 2,
                                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                                        1, 2, 2, 3, 3, 4, 4, 4, 5, 5,
                                                        5, 5, 1};
    Double_t fBinsPi07TeVPCMEMCPt[39]               = { 0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
                                                        1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0,
                                                        2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8,
                                                        4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0,
                                                        10.0, 12.0, 16.0, 20.0};
    Int_t fBinsPi07TeVPCMEMCPtRebin[38]             = { 2, 2, 2, 2, 2, 2, 8, 5, 4, 4,
                                                        2, 2, 2, 2, 2, 2, 2, 2,
                                                        2, 2, 2, 2, 2, 2, 2, 2,
                                                        2, 2, 2, 2, 2, 2, 4, 4, 3,
                                                        5, 8, 1};
    Double_t fBinsPi07TeVEMCPt[39]                  = { 0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
                                                        1.2, 1.4, 1.6, 1.8, 2.0,
                                                        2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8,
                                                        4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0,
                                                        10.0, 11.0, 12.0, 14.0, 16.0, 20.0, 25.0};
    Int_t fBinsPi07TeVEMCPtRebin[38]                = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                                        4, 4, 4, 2, 2, 2, 2,
                                                        2, 2, 2, 2, 2, 2, 2,
                                                        2, 2, 2, 2, 2, 2, 1, 5,
                                                        5, 5, 5, 8, 2, 2};

    Double_t fBinsEta7TeVPt[18]                     = { 0.0, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0,
                                                        3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0};
    Double_t fBinsEta7TeVPCMPHOSPt[19]              = { 0.0, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0,
                                                        3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0};
    Double_t fBinsEta7TeVPHOSPt[18]                 = { 0.0, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0,
                                                        3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0};
    Double_t fBinsEta7TeVPCMEMCPt[19]               = { 0.0, 0.4, 0.6, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 3.5,
                                                        4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0, 20, 25.0};
    Int_t fBinsEta7TeVPtRebin[17]                   = { 8, 7, 7, 4, 4, 4, 4, 4, 5, 5,
                                                        5, 5, 5, 5, 6, 8, 8};
    Int_t fBinsPi0EtaBinning7TeVPtRebin[17]         = { 8, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                                        2, 2, 2, 2, 2, 2, 4};
    Int_t fBinsEta7TeVPCMPHOSPtRebin[18]            = { 8, 8, 8, 8, 5, 5, 5, 5, 5, 5,
                                                        10, 10, 8, 8, 8, 8, 8};
    Int_t fBinsEta7TeVPHOSPtRebin[17]               = { 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                                        8, 8, 8, 8, 8, 8, 8};

    Int_t fBinsEta7TeVPCMEMCPtRebin[18]             = { 2, 2, 2, 12, 10, 8, 8, 8, 6, 8,
                                                        8, 10, 10, 16, 16, 16, 16, 20};
    Int_t fBinsPi0EtaBinning7TeVPCMEMCPtRebin[18]   = { 2, 2, 2, 4, 4, 2, 2, 2, 2, 2,
                                                        2, 4, 4, 4, 4, 4, 16, 20};
    Int_t fBinsEta7TeVEMCPtRebin[18]                = { 2, 2, 2, 2, 2, 2, 10, 10, 10, 8,
                                                        8, 8, 8, 10, 12, 20, 2, 2};
    Int_t fBinsPi0EtaBinning7TeVEMCPtRebin[18]      = { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                                        2, 2, 2, 2, 2, 4, 2, 2};
    Double_t fBinsEta7TeVDalitzPt[10]               = { 0., 0.6, 1.0, 1.4, 1.8, 2.2, 2.8, 4.4, 6., 10.};
    Int_t fBinsEta7TeVDalitzPtRebin[9]              = { 10, 10, 10, 10, 10, 10, 10, 10, 10};

    Int_t fBinsPi0EtaBinning7TeVDalitzPtRebin[9]    = { 8, 2, 2, 2, 2, 2, 4, 4, 4};

    Double_t fBinsEtaPrim7TeVPt[8]                  = { 0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 10.0};
    Int_t fBinsEtaPrim7TeVPtRebin[7]                = { 8, 2, 2, 2, 2, 2, 2};

    Double_t fBinsDirGamma7TeVPt[25]                = { 0.0, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,
                                                        1.8, 2.0, 2.2, 2.4, 2.7, 3.0, 3.5, 4.0,
                                                        4.5, 5.0, 6.0, 7.0, 9.0, 12., 16., 20.};
    Int_t fBinsDirGamma7TeVPtRebin[24]              = { 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
                                                        2, 2, 2, 3, 3, 4, 4, 4, 5, 5,
                                                        5, 5, 5};
    Int_t fBinsDirGamma7TeVEMCPtRebin[23]           = { 5, 5, 5, 5, 5, 4, 4, 3, 2, 2,
                                                        2, 2, 2, 3, 3, 4, 4, 4, 5, 5,
                                                        5, 5, 5};

    // Eta->pi+pi-pi0
    Double_t fBinsEtaPiPlPiMiPiZero7TevPtPCM[14]         = {0,1,1.4,1.6,1.8,2,2.5,3,3.5,4,5,6.,8.,10.};
    Int_t fBinsEtaPiPlPiMiPiZero7TevPtRebinPCM[13]       = {3,6,6,6,3,3,3,3,3,3,3,3,4};

    Double_t fBinsEtaPiPlPiMiPiZero7TevPtPCMEMC[11]         = {0,1,1.5,2,2.5,3,3.5,4,5,6.,10.};
    Int_t fBinsEtaPiPlPiMiPiZero7TevPtRebinPCMEMC[10]       = {4,4,4,4,4,4,4,6,6,6};

    Double_t fBinsEtaPiPlPiMiPiZero7TevPtEMC[16]         = {0,1,1.5,2,2.5,3,3.5,4,5,6.,8.,10.,12,14,16,20.};
    Int_t fBinsEtaPiPlPiMiPiZero7TevPtRebinEMC[15]       = {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};

    Double_t fBinsEtaPiPlPiMiPiZero7TevPtPCMPHOS[12]         = {0,1,1.5,2,2.5,3,3.5,4,5,6.,8.,10.};
    Int_t fBinsEtaPiPlPiMiPiZero7TevPtRebinPCMPHOS[11]       = {4,4,4,4,4,4,4,4,4,4,4};

    Double_t fBinsEtaPiPlPiMiPiZero7TevPtPHOS[11]         = {0,1,1.5,2.,4.,6.,8.,10.,12.,14.,17.};
    Int_t fBinsEtaPiPlPiMiPiZero7TevPtRebinPHOS[10]       = {4,4,4,4,4,4,4,4,4,4};

    // omega->pi+pi-pi0
    Double_t fBinsOmegaPiPlPiMiPiZero7TevPtPCM[14]         = {0,1,1.4,1.6,1.8,2,2.5,3,3.5,4,5,6.,8.,10.};
    Int_t fBinsOmegaPiPlPiMiPiZero7TevPtRebinPCM[13]       = {10,10,7,7,7,5,5,5,5,5,5,10,13};

    Double_t fBinsOmegaPiPlPiMiPiZero7TevPtPCMEMC[11]         = {0,1,1.5,2,2.5,3,3.5,4,5,6.,10.};
    Int_t fBinsOmegaPiPlPiMiPiZero7TevPtRebinPCMEMC[10]       = {10,10,10,10,10,10,8,8,7,7};

    Double_t fBinsOmegaPiPlPiMiPiZero7TevPtEMC[16]         = {0,1,1.5,2,2.5,3,3.5,4,5,6.,8.,10.,12,14,16,20.};
    Int_t fBinsOmegaPiPlPiMiPiZero7TevPtRebinEMC[15]       = {10,10,10,10,10,10,8,8,10,10,10,10,8,10,10};

    Double_t fBinsOmegaPiPlPiMiPiZero7TevPtPCMPHOS[12]         = {0,1,1.5,2,2.5,3,3.5,4,5,6.,8.,10.};
    Int_t fBinsOmegaPiPlPiMiPiZero7TevPtRebinPCMPHOS[11]       = {5,5,5,8,8,8, 8,8,10,16,16};

    Double_t fBinsOmegaPiPlPiMiPiZero7TevPtPHOS[11]         = {0,1,1.5,2,4.,6.,8.,10.,12.,14.,17.};
    Int_t fBinsOmegaPiPlPiMiPiZero7TevPtRebinPHOS[10]       = {12,12,12,12,12,12,12,12,12,12};

#endif