#ifndef GAMMACONV_ExtractSignalBinningXeXe5440GeV
#define GAMMACONV_ExtractSignalBinningXeXe5440GeV

    //****************************************************************************************************
    //***************************** Pt binning for XeXe, 5.44 TeV ***********************************
    //****************************************************************************************************
    Double_t fBinsPi0XeXe5440GeVPt[25]              = { 0.0, 0.3, 0.4, 0.5, 0.6,    0.8, 1.0, 1.2, 1.4, 1.6,
                                                        1.8, 2.0, 2.2, 2.4, 2.6,    3.0, 3.5, 4.0, 5.0, 6.0,
                                                        8.0, 10.0, 12.0, 14.0, 20.0};
    Double_t fBinsPi0XeXe5440GeVCentPt[25]          = { 0.0, 0.3, 0.4, 0.5, 0.6,    0.8, 1.0, 1.2, 1.4, 1.6,
                                                        1.8, 2.0, 2.2, 2.4, 2.6,    3.0, 3.5, 4.0, 5.0, 6.0,
                                                        8.0, 10.0, 12.0, 14.0, 20.0};
    Double_t fBinsPi0XeXe5440GeVCentPCMPt[25]       = { 0.0, 0.3, 0.4, 0.5, 0.6,    0.8, 1.0, 1.2, 1.4, 1.6,
                                                        1.8, 2.0, 2.2, 2.4, 2.6,    3.0, 3.5, 4.0, 5.0, 6.0,
                                                        8.0, 10.0, 12.0, 14.0, 20.0};
    Double_t fBinsPi0XeXe5440GeVCentPCMEMCPt[21]    = { 0.0, 0.8, 1.0, 1.2, 1.4,    1.6, 1.8, 2.0, 2.2, 2.4,
                                                        2.6, 3.0, 3.5, 4.0, 5.0,    6.0, 8.0, 10.0, 12.0, 14.0,
                                                        20.0};
    Double_t fBinsPi0XeXe5440GeVCentPCMPHOSPt[21]   = { 0.0, 0.6, 0.8, 1.0, 1.2,    1.4, 1.6, 1.8, 2.0, 2.2,
                                                        2.4, 2.6, 3.0, 3.5, 4.0,    6.0, 8.0, 10.0, 12.0, 14.0,
                                                        20.0};
    Double_t fBinsPi0XeXe5440GeVCentEMCPt[18]       = { 0.0, 1.4, 1.6, 1.8, 2.0,    2.2, 2.4, 2.6, 3.0, 3.5,
                                                        4.0, 5.0, 6.0, 8.0, 10.0,   12.0, 14.0, 20.0};
    Double_t fBinsPi0XeXe5440GeVCentPHOSPt[21]      = { 0.0, 0.8, 1.0, 1.2, 1.4,    1.6, 1.8, 2.0, 2.2, 2.4,
                                                        2.6, 3.0, 3.5, 4.0, 6.0,    8.0, 10.0, 12.0, 14.0,
                                                        20.0};
    Double_t fBinsPi0XeXe5440GeVPerPt[22]           = { 0.0, 0.3, 0.4, 0.5, 0.6,    0.8, 1.0, 1.2, 1.4, 1.6,
                                                        1.8, 2.0, 2.2, 2.4, 2.6,    3.0, 3.5, 4.0, 6.0, 10.0,
                                                        14.0, 20.0};
    Double_t fBinsPi0XeXe5440GeVPerPCMPt[22]        = { 0.0, 0.3, 0.4, 0.5, 0.6,    0.8, 1.0, 1.2, 1.4, 1.6,
                                                        1.8, 2.0, 2.2, 2.4, 2.6,    3.0, 3.5, 4.0, 6.0, 10.0,
                                                        14.0, 20.0};

    Double_t fBinsPi0XeXe5440GeVPerPCMEMCPt[18]     = { 0.0, 0.8, 1.0, 1.2, 1.4,    1.6, 1.8, 2.0, 2.2, 2.4,
                                                        2.6, 3.0, 3.5, 4.0, 6.0,    10.0, 14.0, 20.0};
    Double_t fBinsPi0XeXe5440GeVPerPCMPHOSPt[19]    = { 0.0, 0.6, 0.8, 1.0, 1.2,    1.4, 1.6, 1.8, 2.0, 2.2,
                                                        2.4, 2.6, 3.0, 3.5, 4.0,    6.0, 10.0, 14.0, 20.0};
    Double_t fBinsPi0XeXe5440GeVPerEMCPt[15]        = { 0.0, 1.4, 1.6, 1.8, 2.0,    2.2, 2.4, 2.6, 3.0, 3.5,
                                                        4.0, 6.0, 10.0, 14.0, 20.0};
    Double_t fBinsPi0XeXe5440GeVPerPHOSPt[18]       = { 0.0, 0.8, 1.0, 1.2, 1.4,    1.6, 1.8, 2.0, 2.2, 2.4,
                                                        2.6, 3.0, 3.5, 4.0, 6.0,    10.0, 14.0, 20.0};
    Int_t fBinsPi0XeXe5440GeVPtRebin[24]            = { 10, 8, 5, 5, 4,     4, 4, 4, 4, 2,
                                                        2, 2, 2, 4, 4,      4, 4, 4, 4, 4,
                                                        10, 12, 12, 12};
    Int_t fBinsPi0XeXe5440GeVPtPCMRebin[24]         = { 10, 8, 8, 4, 4,     4, 2, 2, 2, 2,
                                                        2, 2, 2, 4, 4,      4, 4, 5, 5, 8,
                                                        10, 12, 12, 12};
    Int_t fBinsPi0XeXe5440GeVPtPCMEMCRebin[25]      = { 10, 8, 8, 5, 5,     5, 5, 5, 4, 4,
                                                        4, 4, 4, 4, 4,      4, 4, 4, 4, 5,
                                                        8, 12, 20, 20,      20};
    Int_t fBinsPi0XeXe5440GeVPtPCMPHOSRebin[25]     = { 10, 8, 8, 8, 8,     8, 5, 4, 4, 4,
                                                        4, 4, 4, 4, 4,      4, 4, 5, 8, 8,
                                                        10, 12, 12, 12, 12};
    Int_t fBinsPi0XeXe5440GeVPtEMCRebin[24]         = { 10, 8, 8, 4, 4,     4, 4, 4, 4, 4,
                                                        4, 4, 4, 4, 4,      4, 4, 4, 4, 5,
                                                        5, 5, 10, 12};
    Int_t fBinsPi0XeXe5440GeVPtPHOSRebin[24]        = { 10, 8, 5, 5, 4,     4, 4, 4, 4, 4,
                                                        4, 4, 4, 4, 4,      4, 4, 4, 4, 4,
                                                        5, 8, 10, 12};
    Int_t fBinsPi0XeXe5440GeVPtCentPCMRebin[24]     = { 10, 8, 8, 8, 5,     5, 4, 4, 4, 4,
                                                        4, 4, 4, 5, 5,      5, 5, 5, 8, 8,
                                                        10, 10, 12, 12};
    Int_t fBinsPi0XeXe5440GeVPtCentPCMEMCRebin[20]  = { 10, 10, 8, 5, 5,    5, 5, 5, 5, 5,
                                                        5, 5, 5, 8, 8,      8, 10, 10, 12, 12};
    Int_t fBinsPi0XeXe5440GeVPtCentPCMPHOSRebin[20] = { 10, 8, 8, 8, 8,     8, 8, 8, 8, 8,
                                                        8, 8, 8, 8, 8,      8, 8, 10, 10, 10};
    Int_t fBinsPi0XeXe5440GeVPtCentEMCRebin[17]     = { 10, 8, 8, 8, 5,     5, 5, 5, 5, 8,
                                                        8, 8, 8, 10, 10,    12, 12};
    Int_t fBinsPi0XeXe5440GeVPtCentPHOSRebin[20]    = { 10, 8, 5, 5, 4,     4, 4, 4, 4, 4,
                                                        4, 4, 5, 5, 5,      8, 8, 10, 10, 12};
    Int_t fBinsPi0XeXe5440GeVPtSemiPCMPHOSRebin[20] = { 10, 8, 8, 8, 8,     8, 8, 8, 8, 8,
                                                        8, 8, 8, 8, 10,     10, 10, 10, 10, 10};
    Int_t fBinsPi0XeXe5440GeVPtSemiPHOSRebin[20]    = { 10, 8, 8, 8, 5,     5, 5, 5, 5, 5,
                                                        5, 5, 5, 5, 5,      8, 8, 10, 12};
    Int_t fBinsPi0XeXe5440GeVPtPerPCMRebin[21]      = { 10, 8, 8, 8, 4,     4, 4, 4, 4, 4,
                                                        4, 4, 4, 5, 5,      5, 5, 5, 8, 8,
                                                        12};
    Int_t fBinsPi0XeXe5440GeVPtPerPCMEMCRebin[17]   = { 10, 8, 4, 4, 4,     4, 4, 4, 5, 5,
                                                        5, 5, 5, 8, 8,      8, 12};
    Int_t fBinsPi0XeXe5440GeVPtPerPCMPHOSRebin[18]  = { 10, 8, 8, 8, 8,     8, 8, 8, 8, 8,
                                                        8, 8, 8, 8, 8,      10, 10, 12};
    Int_t fBinsPi0XeXe5440GeVPtPerEMCRebin[14]      = { 10, 8, 4, 4, 4,     4, 5, 5, 5, 5,
                                                        5, 8, 8, 12};
    Int_t fBinsPi0XeXe5440GeVPtPerPHOSRebin[17]     = { 10, 8, 4, 4, 4,     4, 4, 4, 5, 5,
                                                        5, 5, 5, 8, 8,      8, 12};
    Int_t fBinsPi0XeXe5440GeVPtRebinCent[24]        = { 10, 8, 8, 8, 5, 5, 5, 4, 4, 4,
                                                        4, 4, 4, 4, 4, 4, 5, 5, 8, 8,
                                                        8, 12, 12, 12};

    Double_t fBinsEtaXeXe5440GeVPt[8]               = { 0.0, 1.0, 2.0, 3.0, 4.0,    6.0, 8.0, 10.0};
    Int_t fBinsEtaXeXe5440GeVPtRebin[7]             = { 10, 10, 10, 10, 10,   10, 10};

    Int_t fBinsPi0EtaBinningXeXe5440GeVPtRebin[9]   = { 10, 4, 2, 2, 2, 2, 4, 4, 8};

    Double_t fBinsDirGammaXeXe5440GeVPt[20]         = { 0.0, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3,
                                                        2.7, 3.1, 3.5, 4.0, 4.5, 5.5, 6.5, 8.0, 11.0, 14.0};
    Int_t fBinsDirGammaXeXe5440GeVPtRebin[19]       = { 4, 4, 2, 2, 2, 2, 2, 2, 2, 2,
                                                        2, 2, 2, 2, 2, 2, 4, 4, 4};
    Int_t fNBinsClusterXeXe5440GeVPt                =  55;
    Double_t fBinsClusterXeXe5440GeVPt[56]          = { 0.0, 0.7, 0.8, 0.9, 1.0,    1.1, 1.2, 1.3, 1.4, 1.5,
                                                        1.6, 1.7, 1.8, 1.9, 2.0,    2.2, 2.4, 2.6, 2.8, 3.0,
                                                        3.2, 3.4, 3.6, 3.8, 4.0,    4.2, 4.4, 4.6, 4.8, 5.0,
                                                        5.2, 5.4, 5.6, 5.8, 6.0,    6.2, 6.4, 6.6, 6.8, 7.0,
                                                        7.4, 7.8, 8.2, 8.6, 9.0,    9.5, 10,  12,  14,  16,
                                                        18, 20,  25, 30, 35,        40 };
#endif