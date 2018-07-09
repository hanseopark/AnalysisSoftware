#ifndef GAMMACONV_ExtractSignalBinningPbPb2760GeV
#define GAMMACONV_ExtractSignalBinningPbPb2760GeV
    //****************************************************************************************************
    //****************** Pt binning for PbPb, 5.02 TeV ***************************************************
    //****************************************************************************************************
    Double_t fBinsPi0PbPb5TeVPt[16]                 = { 0.0, 1.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.0,
                                                        3.5, 4.0, 5.0, 7.0, 9.0, 12.0};
    Double_t fBinsPi0PbPb5TeVEMCPt[25]              = { 0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
                                                        2.2, 2.4, 2.6, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0,
                                                        12.0, 14.0,16.0, 20.,25.};
    Int_t fBinsPi0PbPb5TeVEMCPtRebin[24]            = { 10, 8, 2, 2, 2, 2, 2, 2, 2, 2,
                                                        2, 2, 2, 2, 2, 2, 4, 4, 4, 4,
                                                        4, 8, 8, 8};
    Double_t fBinsPi0PbPb5TeVPCMEMCPt[25]           = { 0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
                                                        2.2, 2.4, 2.6, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0,
                                                        12.0, 14.0,16.0, 20.,25.};
    Int_t fBinsPi0PbPb5TeVPCMEMCPtRebin[24]         = { 10, 8, 2, 2, 2, 2, 2, 2, 2, 2,
                                                        2, 2, 2, 2, 2, 2, 4, 4, 4, 4,
                                                        4, 8, 8, 8};
    Double_t fBinsPi0PbPb5TeVPtDCA[14]              = { 0.0, 1.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.0,
                                                        3.5, 4.0, 5.0, 7.0};
    Int_t fBinsPi0PbPb5TeVPtRebin[15]               = { 10, 4, 2, 2, 2, 2, 2, 2, 2, 2,
                                                        4, 4, 4, 4, 4};

    Double_t fBinsEtaPbPb5TeVPt[4]                  = { 0.0, 1.0, 3.0, 6.0};
    Int_t fBinsEtaPbPb5TeVPtRebin[3]                = { 10, 8, 8};
    Double_t fBinsEtaPbPb5TeVPtDCA[4]               = { 0.0, 1.0, 3.0, 6.0};
    Double_t fBinsEtaPbPb5TeVEMCPt[23]              = { 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6,
                                                        4.0, 5.0, 6.0, 8., 10., 12., 14., 16., 18., 20., 25., 30., 35.};
    Int_t fBinsEtaPbPb5TeVEMCPtRebin[22]            = { 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                                        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};
    Double_t fBinsEtaPbPb5TeVPCMEMCPt[23]           = { 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6,
                                                        4.0, 5.0, 6.0, 8., 10., 12., 14., 16., 18., 20., 25., 30., 35.};
    Int_t fBinsEtaPbPb5TeVPCMEMCPtRebin[22]         = { 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                                        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};

    Int_t fBinsPi0EtaBinningPbPb5TeVPtRebin[3]      = { 10, 2, 2};

    Double_t fBinsDirGammaPbPb5TeVPt[20]            = { 0.0, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3,
                                                        2.7, 3.1, 3.5, 4.0, 4.5, 5.5, 6.5, 8.0, 11.0, 14.0};
    Int_t fBinsDirGammaPbPb5TeVPtRebin[19]          = { 4, 4, 2, 2, 2, 2, 2, 2, 2, 2,
                                                        2, 2, 2, 2, 2, 2, 4, 4, 4};

#endif