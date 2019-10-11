#ifndef GAMMACONV_ExtractSignalBinningpPb8TeV
#define GAMMACONV_ExtractSignalBinningpPb8TeV

    //****************************************************************************************************
    //****************** Pt binning for pPb, 8 TeV *******************************************************
    //****************************************************************************************************
    Double_t fBinsPi0CombpPb8TeVPt[59] = {
        0.0, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,       1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6,
        3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0,       8.5, 9.0, 10., 11., 12., 13., 14., 15., 16., 17.,
        18., 20., 22., 26., 30., 35., 40., 45., 50., 55.,       60., 65., 70., 80., 100, 125, 150, 175, 200 };
    Double_t fBinsPi0pPb8TeVPt[34] = {
        0.0, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,       1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6,
        3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10., 12.,       16., 20., 25., 30. };
    Double_t fBinsPi0pPb8TeVPtDCA[24] = {
        0.0, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,       2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 8.0, 10.,
        12., 16., 20., 25.};
    Double_t fBinsPi0pPb8TeVPtPCMEMC[31] = {
        0.0, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2,       2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5,
        5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 10., 12., 16.,       25.};
    Double_t fBinsPi0pPb8TeVPCMEMCalTrigger1Pt[44] = {
        0.0, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,       1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6,
        3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0,       8.5, 9.0, 10., 11., 12., 13., 14., 15., 16., 17.,
        18., 20., 22., 26.};
    Double_t fBinsPi0pPb8TeVPCMEMCalTrigger2Pt[44] = {
        0.0, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,       1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6,
        3.8, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 8.5, 9.0, 10.,       11., 12., 13., 14., 15., 16., 17., 18., 20., 22.,
        26., 30., 35., 40.};
    Double_t fBinsPi0pPb8TeVPtEMC[33]= {
        0.0, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,       1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6,
        3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.,       12., 16., 20.};
    Double_t fBinsPi0pPb8TeVEMCalTrigger1Pt[45] = {
        0.0, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,       1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6,
        3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0,       8.5, 9.0, 10., 11., 12., 13., 14., 15., 16., 17.,
        18., 20., 25., 30., 40.};
    Double_t fBinsPi0pPb8TeVEMCalTrigger2Pt[43] = {
        0.0, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,       1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6,
        3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10., 11.,       12., 13., 14., 15., 16., 17., 18., 20., 22., 26.,
        30., 35., 40.};

    Double_t fBinsPi0pPb8TeVPtmEMC[59] = {
        0.0, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,       1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6,
        3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0,       8.5, 9.0, 10., 11., 12., 13., 14., 15., 16., 17.,
        18., 20., 22., 26., 30., 35., 40. ,45., 50., 55.,       60., 65., 70., 80., 100, 125, 150, 175, 200};
    Double_t fBinsPi0pPb8TeVPHOSPt[39] = {
        0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4,       1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
        3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10.,       12., 16., 20., 25., 30., 35., 40., 45., 50.0};
    Double_t fBinsPi0pPb8TeVDalitzPt[23] = {
        0.0, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,       2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.6, 4.0, 5.0, 6.0,
        8.0, 10., 15.};
    Int_t fBinsPi0pPb8TeVDalitzPtRebin[22] = {
        5, 5, 5, 5, 4, 4, 4, 4, 4, 4,       4, 4, 5, 5, 5, 5, 5, 5, 5, 5,
        8, 8};
    Int_t fBinsPi0pPb8TeVPtRebin[33] = {
        2, 8, 4, 2, 1, 1, 1, 1, 1, 1,       1, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 4, 4, 4,       5, 5, 5};
    Int_t fBinsPi0pPb8TeVEMCPtRebin[32] = {
        2, 2, 2, 2, 2, 2, 2, 4, 2, 2,       2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 4, 4, 4, 4, 4, 8,       2, 2};

    Int_t fBinsPi0pPb8TeVEMCTrigger1PtRebin[44] = {
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,       2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 4, 4, 5, 4, 4, 4, 5, 4,       4, 4, 4, 4, 4, 5, 5, 5, 8,10,
        8, 2, 2, 2};
    Int_t fBinsPi0pPb8TeVEMCTrigger2PtRebin[42] = {
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,       2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 4, 5,       5, 4, 5, 8, 8, 8, 8, 5, 2, 2,
        2, 2};
    Int_t fBinsPi0pPb8TeVPCMEMCPtRebin[30] = {
        2, 2, 4, 2, 2, 2, 2, 2, 2, 2,       2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        4, 4, 4, 4, 5, 5, 5, 8,16, 2};
    Int_t fBinsPi0pPb8TeVPCMEMCTrigger1PtRebin[43] = {
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,       2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 5, 5, 5, 4, 4, 5, 5, 5, 5,       5, 5, 5, 5, 5, 5, 5, 5,10,10,
        10,10,10};
    Int_t fBinsPi0pPb8TeVPCMEMCTrigger2PtRebin[43] = {
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,       2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2,10,16,      10,10, 8,10, 8,16, 8, 8,10,16,
        10,10,10};
    Int_t fBinsPi0pPb8TeVPtmEMCRebin[59] = {
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,       1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 2, 2, 2, 2, 2, 2, 2, 2,       2, 2, 2, 2, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 5, 5,       5, 5, 5, 5, 5,10,10,10,10};
    Int_t fBinsPi0pPb8TeVPCMPHOSPtRebin[36] = {
        7, 7, 7, 7, 4, 3, 3, 3, 3, 3,       3, 3, 3, 3, 3, 3, 3, 3, 4, 4,
        4, 4, 4, 4, 4, 5, 7, 7, 7, 7,       7, 7, 7, 7, 10, 10};
    Int_t fBinsPi0pPb8TeVPHOSPtRebin[38] = {
        5, 4, 4, 4, 4, 4, 2, 2, 2, 2,       2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,       4, 4, 5, 5, 8, 8, 8, 8};
    Int_t fBinsPi0pPb8TeVEMCDalitzPtRebin[22] = {
        5, 5, 8, 5, 4, 4, 4, 4, 4, 4,       4, 4, 5, 5, 5, 5, 5, 5, 5, 5,
        8, 8};
    Double_t fBinsEtaCombpPb8TeVPt[26] = {
        0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,       6.0, 7.0, 8.0, 9.0, 10., 12., 14., 16., 18., 20.,
        25., 30., 35., 40., 50., 60.};
    Double_t fBinsEtapPb8TeVPt[19] = {
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0,       6.0, 7.0, 8.0, 10., 12., 14., 16., 18., 20.};
    Double_t fBinsEtapPb8TeVPtDCA[17] = {
        0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4,       1.8, 2.4, 3.5, 5.0, 7.0, 10., 14.};
    Double_t fBinsEtapPb8TeVEMCPt[20] = {
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,       5.0, 6.0, 7.0, 8.0, 10., 12., 14., 16., 18., 20.};
    Double_t fBinsEtapPb8TeVPCMEMCPt[21] = {
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,       5.0, 6.0, 7.0, 8.0, 9.0, 10., 12., 14., 16., 18.,
        20.};
    Double_t fBinsEtapPb8TeVTrigger1Pt[23] = {
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,       5.0, 6.0, 7.0, 8.0, 9.0, 10., 12., 14., 16., 18.,
        20., 25., 30.};
    Double_t fBinsEtapPb8TeVEMCTrigger1Pt[24] = {
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,       5.0, 6.0, 7.0, 8.0, 9.0, 10., 12., 14., 16., 18.,
        20., 25., 30., 40.};
    Double_t fBinsEtapPb8TeVTrigger2Pt[24] = {
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0,       6.0, 8.0, 10., 12., 14., 16., 18., 20., 25., 30.,
        35., 40., 50., 60.};
    Double_t fBinsEtapPb8TeVPCMTrigger1Pt[26] = {
        0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6,       4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10., 12.,
        14., 16., 20., 25., 30., 40.};
    Double_t fBinsEtapPb8TeVPCMTrigger2Pt[24] = {
        0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6,       4.0, 5.0, 6.0, 8.0, 10., 12., 14., 16., 18., 20.,
        25., 30., 35., 40.};
    Double_t fBinsEtapPb8TeVPHOSPt[24] = {
        0., 0.3, 0.5, 0.7, 0.9, 1.1, 1.4, 1.8, 2.2, 3.0,
        4.,  5.,  6., 8.,  10,  12., 16., 20., 25., 30., 35, 40, 45, 50};
    Double_t fBinsEtapPb8TeVPCMPHOSPt[21] = {
        0., 0.3, 0.5, 0.7, 1.1, 1.4, 1.8, 2.2, 2.6, 3.0,
        3.5, 4.,  5.,  6., 8.,  10,  12., 16., 20., 25.,
        30.};
    Int_t fBinsEtapPb8TeVPtRebin[22] = {
        10, 10,  5,  5,  5, 5,  8,  8,  8,  8,
        4,   5,  8,  8,  8, 8,   10, 10, 10, 10,
        10,  10};
    Int_t fBinsEtapPb8TeVEMCPtRebin[21] = {
        4, 4, 4,16,10,16, 8, 8,10, 8,       8, 8, 8, 8,10,16,16,16,20, 4,
        4};
    Int_t fBinsEtapPb8TeVPCMEMCPtRebin[21] = {
        10, 10,  10 /*1st*/,  10,  10, 10,   10,  10,  10,  10,
        10,   10,  10,  10,  16, 20,  20, 20, 20, 20,
        20};
    Int_t fBinsEtapPb8TeVPCMTrigger1PtRebin[26] = {
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5,       5, 5, 5, 6, 6, 6, 6, 6, 8, 8,
        8, 8, 10, 10, 16, 2};
    Int_t fBinsEtapPb8TeVPCMTrigger2PtRebin[23] = {
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5,       5, 8, 8, 8, 8, 8, 8, 8,10,20,
        5};
    Int_t fBinsEtapPb8TeVEMCTrigger1PtRebin[25] = {
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5,       5, 8, 8, 8, 6, 8, 5, 8, 8,10,
        10, 20, 2};
    Int_t fBinsEtapPb8TeVEMCTrigger2PtRebin[25] = {
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5,       5, 5, 5,16,16,10,10, 8, 8, 8,
        10, 10, 10, 10,10};
    Int_t fBinsEtapPb8TeVPCMEMCTrigger1PtRebin[24] = {
        5, 5, 5, 5, 5, 5, 5, 5,10,10,       8, 8,10,10,12, 5,10,10,20, 5,
        5, 5};
    Int_t fBinsEtapPb8TeVPCMEMCTrigger2PtRebin[23] = {
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5,       5,20,20,20,20,16,20,16,20,20,
        2};
    Int_t fBinsPi0EtaBinningpPb8TeVPtRebin[19]  = {
        8, 1, 1, 1, 1, 1, 1, 2, 2, 2,       2, 4, 4, 4, 4, 4, 4, 4, 4};
    Int_t fBinsEtapPb8TeVPHOSPtRebin[23] = {
        19, 10, 10, 10, 10, 10, 10, 10, 10, 10,     10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
    Int_t fBinsEtapPb8TeVPCMPHOSPtRebin[20]  = {
        16, 16, 16, 16, 16, 10, 10, 10, 10, 10,     10, 16, 16, 16, 16, 21, 21, 21, 21, 21};
    // DIR GAMMA
    Double_t fBinsDirGammapPb8TeVPt[26] = {
        0.0, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6,
        1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2, 3.6, 4.0, 4.8,
        5.6, 6.4, 7.2, 8.0, 10.0, 14.0};
    Int_t fBinsDirGammapPb8TeVPtRebin[25] = {
        4, 2, 1, 1, 1, 1, 1, 1, 1, 1,       1, 1, 1, 1, 1, 1, 1, 1, 2, 2,
        2, 2, 4, 4, 4};
    Double_t fBinsDirGammapPb8TeVPCMEMCPt[29] = {
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0,
        1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2,
        3.6, 4.0, 4.8, 5.6, 6.4, 7.2, 8.0, 10.0, 14.0};
    Int_t fBinsDirGammapPb8TeVPCMEMCPtRebin[28] = {
        10, 4, 4, 2, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
        2, 4, 4, 4, 4, 5, 5, 5 };
#endif
