#ifndef GAMMACONV_ExtractSignalBinningPbPb5TeV
#define GAMMACONV_ExtractSignalBinningPbPb5TeV

//****************************************************************************************************
//****************** Pt binning for PbPb, 5.02 TeV ***************************************************
//****************************************************************************************************

//****************************************************************************************************
//****************** Pi0 *****************************************************************************
//****************************************************************************************************

// PCM
Double_t fBinsPi0PbPb5TeVPt[29]                 = { 0.0, 0.4, 0.6, 0.8, 1.,
                                                    1.2, 1.4, 1.6, 1.8, 2.,
                                                    2.2, 2.4, 2.6, 2.8, 3.,
                                                    3.2, 3.4, 3.6, 3.8, 4.,
                                                    4.5, 5., 5.5, 6., 7.,
                                                    8., 10., 12., 20.};

Int_t fBinsPi0PbPb5TeVPtRebin[28]               = { 2, 2, 2, 2, 2,
                                                    2, 2, 2, 2, 2,
                                                    2, 2, 2, 2, 2,
                                                    2, 2, 2, 2, 2,
                                                    2, 2, 4, 4, 4,
                                                    4, 8, 8};

Double_t fBinsPi0PbPb5TeVPtDCA[14]              = { 0.0, 1.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.0,
						    3.5, 4.0, 5.0, 7.0};

// EMCal:
Double_t fBinsPi0PbPb5TeVEMCPt[34]                  =  {  0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
                                                          2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                                          4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10., 12., 14.,
                                                          16., 20., 25., 30.};
    Int_t fBinsPi0PbPb5TeVEMCPtRebin[33]            =  {  8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                                          8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                                          8, 8, 8, 8, 8, 8, 16, 16, 16, 16,
                                                          16, 16, 16};
// PCM-EMCal:
    Double_t fBinsPi0PbPb5TeVPCMEMCPt[34]           =  {  0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
                                                          2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                                          4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10., 12., 14.,
                                                          16., 20., 25., 30.};
    Int_t fBinsPi0PbPb5TeVPCMEMCPtRebin[33]         =  {  10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                                                          10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                                                          10, 10, 10, 10, 10, 10, 10, 10, 16, 16,
                                                          16, 16, 16};

// PCM binning for AOD validation
// Double_t fBinsPi0PbPb5TeVPt[11]                  = { 0.0, 0.4, 5.0, 6.0, 7.0, 8.0, 10.0,
//                                                    12.0, 16.0, 20., 25.};
// Int_t fBinsPi0PbPb5TeVPtRebin[10]               = { 10, 4, 4, 4, 4, 4,
//                                                    4, 4, 4, 8};


//****************************************************************************************************
//****************** Eta *****************************************************************************
//****************************************************************************************************
// 400 inv mass bins

// PCM
Double_t fBinsEtaPbPb5TeVPt[11]                  = { 0., 1., 2., 2.5, 3.,3.5, 4., 6., 8., 10., 20.};

Int_t fBinsEtaPbPb5TeVPtRebin[10]                = { 10, 4, 4, 8, 8, 8, 8 , 8, 8, 8};

Double_t fBinsEtaPbPb5TeVPtDCA[4]               = { 0.0, 1.0, 3.0, 6.0};

// EMCal
Double_t fBinsEtaPbPb5TeVEMCPt[13]              = { 0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 25.0, 30.0};
Int_t fBinsEtaPbPb5TeVEMCPtRebin[12]            = { 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
Int_t fBinsPi0EtaBinningPbPb5TeVEMCPtRebin[12]  = { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
// PCM-EMCal:
Double_t fBinsEtaPbPb5TeVPCMEMCPt[23]           = { 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6,
                                                    4.0, 5.0, 6.0, 8., 10., 12., 14., 16., 18., 20., 25., 30., 35.};
Int_t fBinsEtaPbPb5TeVPCMEMCPtRebin[22]         = { 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                                                    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};


//****************************************************************************************************
//****************** Pi0EtaBinning *******************************************************************
//****************************************************************************************************

// PCM
Int_t fBinsPi0EtaBinningPbPb5TeVPtRebin[10]     = { 10, 2, 2, 2, 2, 2, 2 , 2, 2, 2};

//****************************************************************************************************
//****************** Gamma ***************************************************************************
//****************************************************************************************************

Double_t fBinsDirGammaPbPb5TeVPt[20]            = { 0.0, 0.4, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3,
						    2.7, 3.1, 3.5, 4.0, 4.5, 5.5, 6.5, 8.0, 11.0, 14.0};

Int_t fBinsDirGammaPbPb5TeVPtRebin[19]          = { 4, 4, 2, 2, 2, 2, 2, 2, 2, 2,
						    2, 2, 2, 2, 2, 2, 4, 4, 4};


#endif
