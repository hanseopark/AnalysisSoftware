#ifndef GAMMACONV_ExtractSignalBinningPbPb2760GeV
#define GAMMACONV_ExtractSignalBinningPbPb2760GeV

    //****************************************************************************************************
    //***************************** Pt binning for PbPb 2010, 2.76 TeV ***********************************
    //****************************************************************************************************
    Double_t fBinsPi0PbPb2760GeVPt[25]              = { 0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
                                                        2.2, 2.4, 2.6, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0,
                                                        12.0, 14.0,16.0, 20.,25.};
    Double_t fBinsPi0PbPb2760GeVPtNew[18]           = { 0.0, 0.5, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2,
                                                        2.4, 2.6, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0};
    Double_t fBinsPi0PbPb2760GeVPeripheralPt[16]    = { 0.0, 0.5, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5,
                                                        3.0, 4.0, 6.0, 8.0, 10.0, 12.0 };
    Double_t fBinsPi0PbPb2760GeVPtDCA[16]           = { 0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
                                                        2.25, 2.5,3.0, 4.0, 6.0, 12.};
    Double_t fBinsPi0PbPb2760GeVPtDCAPer[12]        = { 0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 2.0, 2.5, 3.,
                                                        6., 10.};
    Double_t fBinsEtaPbPb2760GeVPtDCA[14]           = { 0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0,
                                                        4.0, 6.0, 10., 12.};
    Int_t fBinsPi0PbPb2760GeVPtRebin[24]            = { 10, 8, 2, 2, 2, 2, 2, 2, 2, 2,
                                                        2, 2, 2, 2, 2, 2, 4, 4, 4, 4,
                                                        4, 8, 8, 8};
    Int_t fBinsPi0PbPb2760GeVPtRebinNew[17]         = { 10, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                                                        4, 4, 4, 4, 5, 5, 8};
    Int_t fBinsPi0PbPb2760GeVPeripheralPtRebin[15]  = { 10, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                                                        4, 4, 8, 8, 8};
    Double_t fBinsPi0PbPb2760GeVDalitzPt[10]        =  {0, 0.6, 1., 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 15.};
    Int_t fBinsPi0PbPb2760GeVDalitzPtRebin[9]       =  {5, 5, 5, 5, 5, 5, 5, 5, 5};

    Double_t fBinsEtaPbPb2760GeVDalitzPt[5]         =  {0.0, 1.5, 2.0, 4.0, 7.0};
    Int_t fBinsEtaPbPb2760GeVDalitzPtRebin[4]       =  {10, 5, 5, 5};
    Double_t fBinsEtaPbPb2760GeVPt[5]               = { 0.0, 1.5, 2.0, 4.0, 7.0};
    Int_t fBinsEtaPbPb2760GeVPtRebin[4]             = { 10, 8, 5, 5};

    Int_t fBinsPi0EtaBinningPbPb2760GeVPtRebin[4]   = { 10, 2, 2, 2};
    Int_t fBinsPi0EtaBinningPbPb2760GeVDalitzPtRebin[4]=  {10, 2, 2, 2};

    Double_t fBinsDirGammaPbPb2760GeVPt[20]         = { 0.0, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3,
                                                        2.7, 3.1, 3.5, 4.0, 4.5, 5.5, 6.5, 8.0, 11.0, 14.0};
    Int_t fBinsDirGammaPbPb2760GeVPtRebin[19]       = { 4, 4, 2, 2, 2, 2, 2, 2, 2, 2,
                                                        2, 2, 2, 2, 2, 2, 4, 4, 4};

    //****************************************************************************************************
    //***************************** Pt binning for PbPb 2011, 2.76 TeV ***********************************
    //****************************************************************************************************
    //same as 10h binning but for the last bins = EMCal bins {4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 20.0, 30.0};
    Double_t fBinsPi0PbPb2760GeVPtLHC11h[27]        = { 0.,  0.4, 0.6, 0.8, 1.,  1.2, 1.4, 1.6, 1.8, 2.,
                                                        2.2, 2.4, 2.6, 3.,  3.5, 4.,  5.,  6.,  8.,  10.,
                                                        12., 14., 16., 18., 20., 25., 30.};
    Int_t fBinsPi0PbPb2760GeVPtLHC11hRebin[26]      = {    10,   8,   2,   2,   2,   2,   2,   2,   2,   2,
                                                            2,   2,   2,   2,   2,   2,   4,   4,   4,   4,
                                                            5,   8,   8,   8,   10,  10 };
    Int_t fBinsPi0PbPb2760GeVPtLHC11hSemicRebin[26] = {    10,   5,   2,   2,   2,   2,   2,   2,   2,   2,
                                                            2,   2,   2,   2,   2,   2,   4,   4,   4,   5,
                                                            5,   8,   8,   8,   10,  10 };
    Double_t fBinsEtaPbPb2760GeVPtLHC11h[17]        = { 0.0, 0.6, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0 ,3.5, 4.0,
                                                        5.0, 6.0, 8.0, 10,  12., 15., 19.};
    Int_t fBinsEtaPbPb2760GeVPtRebinLHC11h[16]      = {     10,  8,   8,   4,   4,   4,   5,  5,    8,   8,
                                                            8,   8,   10,  10,  10, 10};
    Double_t fBinsEtaPbPb2760GeVPtLHC11hLessBins[13]= { 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.,
                                                        12.0, 15.0, 19.0};
    Int_t fBinsEtaPbPb2760GeVPtRebinLHC11hLessBins[14] = {    10,   8,   5,   5,   5,   5,   5,   8, 10, 10,
                                                        10, 10, 10, 10};

    Double_t fBinsDirGammaPbPb2760GeVPtLHC11h[23]   = { 0.0, 0.4, 0.9, 1.1,  1.3,  1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 3.,
                                                        3.3, 3.7, 4.1, 4.6,  5.4,  6.2,  7.,  8., 11., 14.};
    Int_t fBinsDirGammaPbPb2760GeVPtLHC11hRebin[22] = {2,2, 2,  2, 2, 2,    2,    2,   2,  2,  2,    2,   2,   2,  2,
                                                           2,   2,   4,   4,   4,   4,   4};
    Int_t fBinsDirGammaPbPb2760GeVPtLHC11hSemicRebin[22] = {2, 2, 2, 2, 2, 2,    2,    2,   2,  2,  2,    2,   2,   2,  2,
                                                           2,   2,   4,   4,   4,   4,   5};
    Double_t fBinsDirGammaPbPb2760GeVPtLHC11hVar2[19]   = {0.0, 0.4, 0.8, 1.,  1.2,  1.4, 1.6, 1.8, 2., 2.3, 2.6, 3.,
                                                            3.5, 4.,  5.,  6.,  8.,  10., 14.};
    Int_t fBinsDirGammaPbPb2760GeVPtLHC11hRebinVar2[18] = {2, 2, 2, 2,   2,   2,  2,  2,    2,   2,   2,  2,
                                                           2,   2,   4,   4,   4,   4};
    Int_t fBinsDirGammaPbPb2760GeVPtLHC11hSemicRebinVar2[18] = {2, 2, 2, 2,    2,   2,  2,  2,    2,   2,   2,  2,
                                                           2,   2,   4,   4,   4,   5};

    Int_t fBinsPi0PbPb2760GeVPtLHC11hPCMEMCRebin[26]= { 10, 4, 2, 2, 2, 2, 2, 2, 2, 2,
                                                        2, 2, 2, 2, 2, 2, 4, 4, 4, 5,
                                                        5, 10, 10, 10, 10, 10 };
    Double_t fBinsEtaPbPb2760GeVPtLHC11hEMCBins[15] = { 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.,
                                                        12., 15., 18., 24., 30.};
#endif