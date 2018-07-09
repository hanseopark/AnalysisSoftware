#ifndef GAMMACONV_ExtractSignalBinningpp900GeV
#define GAMMACONV_ExtractSignalBinningpp900GeV

    //****************************************************************************************************
    //******************** Pt binning for pp, 0.9 TeV ****************************************************
    //****************************************************************************************************
    Double_t fBinsPi0900GeVPt[12]                   = { 0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5,
                                                        3.5, 4.5};
    Int_t fBinsPi0900GeVPtRebin[11]                 = { 4, 4, 2, 2, 2, 2, 2, 2, 2, 4,
                                                        4};

    Double_t fBinsPi0900GeVPCMEMCPt[12]             = { 0.0, 0.8, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0,
                                                        7.0, 10.0};
    Int_t fBinsPi0900GeVPCMEMCPtRebin[11]           = { 2, 5, 4, 4, 4, 4, 4, 4, 8, 8, 8};

    Double_t fBinsPi0900GeVEMCPt[13]                = { 0.0, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0,
                                                        10.0, 16.0};
    Int_t fBinsPi0900GeVEMCPtRebin[12]              = { 2, 8, 5, 4, 4, 5, 8, 8, 8, 16, 2,
                                                        2};

    Double_t fBinsEta900GeVPt[4]                    =  {0., 0.9, 1.8, 3.0};
    Int_t fBinsEta900GeVPtRebin[3]                  =  {8, 5, 5};
    Int_t fBinsPi0EtaBinning900GeVPtRebin[3]        =  {8, 4, 4};

    Double_t fBinsEta900GeVPCMEMCPt[6]              =  {0., 0.9, 1.8, 3.0, 5.0, 7.0};
    Int_t fBinsEta900GeVPCMEMCPtRebin[5]            =  {5, 16, 10, 20, 5};

    Double_t fBinsDirGamma900GeVPt[14]              = { 0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5,
                                                        3.5, 4.5, 6.5, 8.5};
    Int_t fBinsDirGamma900GeVPtRebin[13]            = { 4, 4, 2, 2, 2, 2, 2, 2, 2, 4,
                                                        4, 4, 4};
#endif