// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion
Int_t fStartPtBin                   = 0;
Int_t fColumn                       = 0; 
Int_t fRow                          = 0;	
Int_t fNBinsPt                      = 0;	
Double_t *fBinsPt                   = NULL; 
Int_t *fNRebin                      = NULL;

//Double_t fBinsPi07TeVPt[22]       =  {0, 0.6, 0.8, 1.0, 1.2,
//                                         1.4, 1.6, 1.8, 2.0, 2.2,
//                                         2.4, 2.6, 2.8, 3.2, 3.6,
//                                         4.0, 4.5, 5.0, 6.0, 8.0,
//                                         10.0, 15.};						
//Int_t    fBinsPi07TeVPtRebin[21]  =  {5, 5, 5, 4, 4,
//                                         4, 4, 4, 4, 4, 
//                                         5, 5, 5, 5, 5,
//                                         5, 5, 8, 8, 8,
//                                         10};
											
Double_t fBinsPi07TeVPt[23]         =  {0, 0.6, 0.7, 0.8, 1.0,
                                        1.2, 1.4, 1.6, 1.8, 2.0,
                                        2.2, 2.4, 2.6, 2.8, 3.2,
                                        3.6, 4.0, 4.5, 5.0, 6.0,
                                        8.0, 10.0, 15.};						
Int_t fBinsPi07TeVPtRebin[22]       =  {5, 5, 5, 5, 4,
                                        4, 4, 4, 4, 4, 
                                        4, 5, 5, 5, 5, 
                                        5, 5, 5, 8, 8,
                                        8, 10};

Double_t fBinsEta7TeVPt[10]         =  {0., 0.6, 1.0, 1.4, 1.8,
                                        2.2, 2.8, 4.4, 6., 10.};
Int_t fBinsEta7TeVPtRebin[9]        =  {10, 10, 10, 10, 10, 
                                        10, 10, 10, 10};
Int_t fBinsPi0EtaBinning7TeVPtRebin[9]    =  {8, 2, 2, 2, 2,
                                              2, 4, 4, 4};


Double_t fBinsPi05023GeVPt[23]      =  {0, 0.6, 0.7, 0.8, 1.0, 
                                        1.2, 1.4, 1.6, 1.8, 2.0,
                                        2.2, 2.4, 2.6, 2.8, 3.0,
                                        3.2, 3.6, 4.0, 5.0, 6.0,
                                        8.0, 10., 15.};
Int_t fBinsPi05023GeVPtRebin[22]    =  {5, 5, 5, 5, 4, 
                                        4, 4, 4, 4, 4,
                                        4, 4, 5, 5, 5,
                                        5, 5, 5, 5, 5,  
                                        8, 8};

Double_t fBinsPi0EMCALDalitz5023GeVPt[23]       =  {0, 0.6, 0.7, 0.8, 1.0, 
                                                    1.2, 1.4, 1.6, 1.8, 2.0, 
                                                    2.2, 2.4, 2.6, 2.8, 3.0, 
                                                    3.2, 3.6, 4.0, 5.0, 6.0,
                                                    8.0, 10., 15.};
Int_t fBinsPi0EMCALDalitz5023GeVPtRebin[22]     =  {5, 5, 8, 5, 4,
                                                    4, 4, 4, 4, 4,
                                                    4, 4, 5, 5, 5,
                                                    5, 5, 5, 5, 5,
                                                    8, 8};

Double_t fBinsEta5023GeVPt[10]      =  {0., 0.6, 1.0, 1.4, 1.8,
                                        2.2, 2.8, 4.4, 6., 10.};
Int_t fBinsEta5023GeVPtRebin[9]     =  {10, 10, 10, 10, 10,
                                        10, 10, 10, 10};
Int_t fBinsPi0EtaBinning5023GeVPtRebin[9]   =  {8, 2, 2, 2, 2, 
                                                2, 4, 4, 4};

Double_t fBinsEtaPrim7TeVPt[8]      =  {0., 0.5, 1., 2., 3.,
                                        4., 6., 10.}; 
Int_t fBinsEtaPrim7TeVPtRebin[7]    =  {8, 2, 2, 2, 2,
                                        2, 2};

Double_t fBinsPi0900GeVPt[12] 		=  {0, 0.4, 0.6, 0.8, 1.0,
                                        1.2, 1.4, 1.6, 2., 2.5,
                                        3.5,4.5};
Int_t fBinsPi0900GeVPtRebin[11] 	=  {4, 2, 2, 2, 2,
                                        2, 2, 2, 2, 4,
                                        4};
Double_t fBinsEta900GeVPt[4] 	    =  {0., 0.9, 1.8, 3.0};
Int_t fBinsEta900GeVPtRebin[4] 	    =  {8, 5, 5};
Int_t fBinsPi0EtaBinning900GeVPtRebin[4]   =  {8, 4, 4};

Double_t fBinsPi02760GeVPt[8]       =  {0.0, 0.6, 1.0, 1.4, 2.0, 
                                        3.0, 5.0, 10.0};
Int_t fBinsPi02760GeVPtRebin[7]     =  {5, 5, 4, 4, 4,
                                        5, 5};


Double_t fBinsEta2760GeVPt[8]       =  {0., 0.5, 1.0, 1.5, 2.0, 
                                        2.5, 4., 6.};
Int_t fBinsEta2760GeVPtRebin[7]     =  {8, 8, 5, 5, 5, 
                                        5, 8};
Int_t fBinsPi0EtaBinning2760GeVPtRebin[7]   =  {8, 2, 2, 2, 2, 
                                                4, 4};

//Double_t fBinsPi0HIPt[14] 	    =  {0, 0.6, 0.8, 1.2, 1.6,
//                                         2.0, 2.4, 3.0, 3.5, 4.0, 
//                                         5.0, 6.0, 10.0, 15.};
//Int_t fBinsPi0HIPtRebin[13] 	    =  {5, 5, 5, 4, 4,
//                                         4, 5, 5, 5, 5,
//                                         8, 8, 10};


Double_t fBinsPi0HIPt[10]           =  {0, 0.6, 1., 2.0, 3.0,
                                        4.0, 5.0, 6.0, 10.0, 15.};
Int_t fBinsPi0HIPtRebin[9]          =  {5, 5, 5, 5, 5,
                                        5, 5, 5, 5};

Double_t fBinsEtaHIPt[5]            =  {0.0, 1.5, 2.0, 4.0, 7.0};
Int_t fBinsEtaHIPtRebin[4]          =  {10, 5, 5, 5};
Int_t fBinsPi0EtaBinningHIPtRebin[4]=  {10, 2, 2, 2};


Int_t fExampleBin                   = 0;
