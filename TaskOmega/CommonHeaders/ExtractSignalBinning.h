// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion
Int_t fStartPtBin = 		0;
Int_t fColumn = 			0; 
Int_t fRow = 				0;	
Int_t fNBinsPt = 			0;	
Double_t *fBinsPt = 		NULL; 
Int_t *fNRebin = 			NULL;
Int_t fExampleBin =     	0;
Int_t fNBinsPeakPt =		12;
Double_t fBinsPeakPt[13] = 						{0.0, 0.4, 0.6, 0.8, 1.0,
												1.2, 1.4, 1.6, 2.0, 3.0,
												4.0, 5.0, 7.0};
Int_t fBinsPeakPtRebin[12] = 					{4, 4, 2, 2, 2, 
												2, 2, 2, 2, 4, 
												4, 8};
Double_t fBinsPeakPtHalf[13] = 					{0.0, 0.2, 0.3, 0.4, 0.5,
												0.6, 0.7, 0.8, 1.0, 1.5,
												2.0, 2.5, 3.5};

//******************** Pt binning for pp, 7 TeV ***************************************************
Double_t fBinsPi07TeVPt[33] =					{0.0, 0.3, 0.4, 0.5, 0.6, 
												0.8, 1.0, 1.2, 1.4, 1.6, 
												1.8, 2.0, 2.2, 2.4, 2.6,
												2.8, 3.0, 3.2, 3.4, 3.6,
												3.8, 4.0, 4.5, 5.0, 5.5,
												6.0, 7.0, 8.0, 10.0, 12.0,
												16.0, 20.0, 25.0};
Double_t fBinsPi07TeVPtDCA[22] =			 	{0.0, 0.3, 0.4, 0.6, 0.8,
												1.0, 1.2, 1.4, 1.6, 2.0,
												2.4, 2.8, 3.2, 3.6, 4.0,
												5.0, 6.0, 8.0, 12.0, 16.0,
												20.0, 25.0};
Int_t fBinsPi07TeVPtRebin[32] =					{2, 2, 1, 1, 1, 
												1, 1, 1, 1, 1, 
												1, 1, 1, 1, 1,
												1, 1, 1, 1, 1,
												1, 1, 2, 2, 2,
												2, 4, 4, 4, 5,
												5, 5};
Double_t fBinsEta7TeVPt[9] = 					{0.1, 2.5,3,3.5,

												 4.0,4.5,5.0, 6.0, 10};

Int_t fBinsEta7TeVPtRebin[8] = 					{2, 10,  10,10,
												10,10, 10,10};

Double_t fBinsOmega7TevPt[12] = 					{0.1, 2.5,3,3.5,

                                                 4.0,4.5,5.0, 6.0, 7.0,8.0,10.,15.};

Int_t fBinsOmega7TevPtRebin[11] = 					{2, 10,  10,10,
                                                10,10, 10,10,10,10,10};

Int_t fBinsPi0EtaBinning7TeVPtRebin[12] = 		{8, 1, 1, 1, 1, 
												1, 1, 2, 2, 2, 
												20, 4};
Double_t fBinsEtaPrim7TeVPt[8] =				{0.0, 0.5, 1.0, 2.0, 3.0,
												4.0, 6.0, 10.0}; 
Int_t fBinsEtaPrim7TeVPtRebin[7] =	 			{8, 2, 2, 2, 2, 
												2, 2};
Double_t fBinsPi07TeVDirectPhotonPt[24] = 		{0.0, 0.3, 0.6, 0.8, 1.0,
												1.2, 1.4, 1.6, 1.8, 2.0,
												2.4, 2.8, 3.2, 3.6, 4.2,
												4.8, 5.8, 7.0, 8.5, 10.0,
												12.0, 16.0, 20.0, 25.0};
Int_t fBinsPi07TeVDirectPhotonPtRebin[23] = 	{2, 2, 1, 1, 1, 
												1, 1, 1, 1, 1, 
												1, 1, 1, 2, 2,
												2, 2, 4, 4, 4,
												5, 5, 5};

//******************** Pt binning for pp, 8 TeV ***************************************************
Double_t fBinsPi08TeVPt[33] =					{0.0, 0.3, 0.4, 0.5, 0.6, 
												0.8, 1.0, 1.2, 1.4, 1.6, 
												1.8, 2.0, 2.2, 2.4, 2.6,
												2.8, 3.0, 3.2, 3.4, 3.6,
												3.8, 4.0, 4.5, 5.0, 5.5,
												6.0, 7.0, 8.0, 10.0, 12.0,
												16.0, 20.0, 25.0};
Double_t fBinsPi08TeVPtDCA[22] =			 	{0.0, 0.3, 0.4, 0.6, 0.8,
												1.0, 1.2, 1.4, 1.6, 2.0,
												2.4, 2.8, 3.2, 3.6, 4.0,
												5.0, 6.0, 8.0, 12.0, 16.0,
												20.0, 25.0};
Int_t fBinsPi08TeVPtRebin[32] =					{2, 2, 1, 1, 1, 
												1, 1, 1, 1, 1, 
												1, 1, 1, 1, 1,
												1, 1, 1, 1, 1,
												1, 1, 2, 2, 2,
												2, 4, 4, 4, 5,
												5, 5};
Double_t fBinsEta8TeVPt[15] = 					{0.0, 0.4, 0.7, 1.0, 1.4,
												1.8, 2.2, 2.6, 3.0, 3.5, 
												4.0, 6.0, 8.0, 12.0, 16.0}; 
Int_t fBinsEta8TeVPtRebin[14] = 				{8, 5, 5, 4, 4,
												4, 4, 5, 5, 5, 
												5, 5, 5, 5};
Int_t fBinsPi0EtaBinning8TeVPtRebin[12] = 		{8, 1, 1, 1, 1, 
												1, 1, 2, 2, 2, 
												2, 4};
												
//******************** Pt binning for pp, 0.9 TeV ***************************************************
Double_t fBinsPi0900GeVPt[12] = 				{0.0, 0.4, 0.6, 0.8, 1.0,
												1.2, 1.4, 1.6, 2.0, 2.5,
												3.5, 4.5};
Int_t fBinsPi0900GeVPtRebin[11] = 				{4, 2, 2, 2, 2, 
												2, 2, 2, 2, 4, 
												4};
Double_t fBinsEta900GeVPt[4] = 					{0., 0.9, 1.8, 3.0};
Int_t fBinsEta900GeVPtRebin[4] = 				{8, 5, 5};
Int_t fBinsPi0EtaBinning900GeVPtRebin[4] =	 	{8, 4, 4};
Double_t fBinsPi0900GeVDirectPhotonPt[8] = 		{0.0, 0.6, 0.8, 1.0, 1.3, 
												2.0, 3.0, 4.5};
Int_t fBinsPi0900GeVDirectPhotonPtRebin[7] = 	{4, 2, 2, 2, 2, 
												2, 4};

//******************** Pt binning for pp, 2.76 TeV ***************************************************
Double_t fBinsPi02760GeVPt[20] = 				{0.0, 0.4, 0.6, 0.8, 1.0,
												1.2, 1.4, 1.6, 1.8, 2.0, 
												2.2, 2.4, 2.6, 3.0, 3.5, 
												4.0, 5.0, 6.0, 8.0, 10.0};
Double_t fBinsPi02760GeVPtDCA[15] =				{0.0, 0.4, 0.6, 0.8, 1.0, 
												1.2, 1.4, 1.6, 2.0, 2.4, 
												3.0, 4.0, 6.0, 8.0, 10.0};
Int_t fBinsPi02760GeVPtRebin[19] = 				{4, 4, 2, 2, 2, 
												2, 2, 2, 2, 2, 
												2, 2, 2, 2, 2, 
												2, 4, 4, 4};
//Double_t  fBinsPi02760GeVPt[13] = 			{0.0, 0.4, 0.7, 1.0, 1.5, 
// 												2.0, 2.5, 3.0, 4.0, 5.0, 
// 												7.0, 9.0, 12.0};
//Int_t fBinsPi02760GeVPtRebin[12] = 			{2, 2, 2, 2, 2, 
// 												2, 2, 2, 2, 4, 
// 												4, 4};
Double_t fBinsEta2760GeVPt[8] = 				{0., 0.5, 1.0, 1.5, 2.0, 
												2.5, 4.0, 6.0};
Double_t fBinsEta2760GeVPtDCA[16] =				{0., 0.4, 0.6, 0.8, 1.0,
												1.125, 1.5, 1.75, 2.0, 2.5,
												3.0, 4.0, 6.0, 10.0, 15.0, 
												20.0 };
Int_t fBinsEta2760GeVPtRebin[7] = 				{8, 8, 5, 5, 5, 
												5, 8};
Int_t fBinsPi0EtaBinning2760GeVPtRebin[9] = 	{8, 2, 2, 2, 2, 
												2, 2, 4, 4};
Double_t fBinsPi02760GeVDirectPhotonPt[16] = 	{0.0, 0.3, 0.6, 0.8, 1.0, 
												1.2, 1.5, 1.8, 2.3, 2.8, 
												3.3, 3.9, 4.5, 5.5, 7.0, 
												10.0};
Int_t fBinsPi02760GeVDirectPhotonPtRebin[15] = 	{4,2,2,2,2,2,2,2,2,2,2,2,4,4,4};

// ***************** Pt binning for PbPb, 2.76 TeV *************************************
Double_t fBinsPi0HIPt[25] =						{0.0, 0.4, 0.6, 0.8, 1.0,
												1.2, 1.4, 1.6, 1.8, 2.0,
												2.2, 2.4, 2.6, 3.0, 3.5,
												4.0, 5.0, 6.0, 8.0, 10.0, 
												12.0, 14.0,16.0, 20.,25.};
Double_t fBinsPi0HIPtNew[18] =  				{0.0, 0.5, 0.8, 1.0, 1.2,
												1.4, 1.6, 1.8, 2.0, 2.2, 
												2.4, 2.6, 3.0, 4.0, 6.0, 
												8.0, 10.0, 12.0};
Double_t fBinsPi0HIPeripheralPt[16] =			{0.0, 0.5, 0.8, 1.0, 1.2,
												1.4, 1.6, 1.8, 2.0, 2.5,
												3.0, 4.0, 6.0, 8.0, 10.0,
												12.0 };
Double_t fBinsPi0HIPtDCA[16] = 					{0.0, 0.4, 0.6, 0.8, 1.0,
												1.2, 1.4, 1.6, 1.8, 2.0,
												2.25, 2.5,3.0, 4.0, 6.0,
												12.};
Double_t fBinsPi0HIPtDCAPer[12] =				{0.0, 0.4, 0.6, 0.8, 1.0,
												1.2, 1.4, 2.0, 2.5, 3., 
												6., 10.};  
Double_t fBinsEtaHIPtDCA[14] =					{0.0, 0.4, 0.6, 0.8, 1.0,
												1.2, 1.5, 2.0, 2.5, 3.0,
												4.0, 6.0, 10., 25.};
Int_t fBinsPi0HIPtRebin[24] = 					{10, 8, 8, 5, 4,
												4, 4, 4, 4, 4,
												4, 4, 4, 4, 4,
												4, 5, 5, 5, 5,
												5, 5, 8, 8};
Int_t fBinsPi0HIPtRebinNew[17] =				{10, 4, 4, 4, 4,
												4, 4, 4, 4, 4,
												4, 4, 4, 4, 5,
												5, 8};
Int_t fBinsPi0HIPeripheralPtRebin[15] = 	  	{10, 4, 4, 4, 4,
												4, 4, 4, 4, 4,
												4, 4, 8, 8, 8};
Double_t fBinsEtaHIPt[5] = 						{0.0, 1.5, 2.0, 4.0, 7.0};
Double_t fBinsEtaHIPtLHC11h[13] =				{0.0, 0.8, 1.0, 1.4, 1.8,
												2.2, 2.6, 3.0 ,3.5, 4.,
												6.0, 8.0, 25.};
// Double_t fBinsEtaHIPtLHC11h[11] =			{0.0, 1.0, 1.5, 2.0, 2.5,
//												3.0, 4.0, 6.0, 8.0, 10.0,
// 												15};

Int_t fBinsEtaHIPtRebin[4] = 					{10, 8, 5, 5};
Int_t fBinsEtaHIPtRebinLHC11h[12] =      		{10, 8, 4, 4, 4,
												4, 5, 5, 5, 5,
												8, 10};
Int_t fBinsEtaHIPtRebinLHC11hFinerBinning[12] =	{10, 8, 6, 5, 5,
												5, 5, 6, 6, 6,
												8, 10};

// Int_t fBinsEtaHIPtRebinLHC11h[10] =			{10, 8, 8, 8, 8,
//												8, 8, 8, 8, 8};
Int_t fBinsPi0EtaBinningHIPtRebin[4] = 			{10, 2, 2, 2};

// Double_t fBinsPi0HIDirectPhotonPt[21] = 		{0.0, 0.4, 0.8, 1.0, 1.2, 
// 												1.4, 1.6, 1.8, 2.0, 2.3, 
// 												2.6, 2.9, 3.2, 3.5, 4.0, 
// 												4.5, 5.5, 6.5, 8.0, 11.0,
// 												14.0}; 
// Int_t fBinsPi0HIDirectPhotonPtRebin[20] = 	{4, 4, 2, 2, 2, 
// 												2, 2, 2, 2, 2, 
// 												2, 2, 2, 2, 2,
// 												2, 2, 4, 4, 4};
//  Double_t fBinsPi0HIDirectPhotonPt[21] = 	{0.0, 0.4, 0.8, 1.0, 1.2, 
// 												1.4, 1.6, 1.8, 2.1, 2.4, 
// 												2.7, 3.0, 3.3, 3.6, 4.0, 
// 												4.5, 5.5, 6.5, 8.0, 11.0,
// 												14.0}; 
//  Int_t fBinsPi0HIDirectPhotonPtRebin[20] = 	{4, 4, 2, 2, 2, 
// 												2, 2, 2, 2, 2, 
// 												2, 2, 2, 2, 2, 
// 												2, 2, 4, 4, 4};

Double_t fBinsPi0HIDirectPhotonPt[20] = 		{0.0, 0.4, 0.8, 1.0, 1.2, 
												1.4, 1.6, 1.8, 2.0, 2.3, 
												2.7, 3.1, 3.5, 4.0, 4.5,
												5.5, 6.5, 8.0, 11.0, 14.0};
Int_t fBinsPi0HIDirectPhotonPtRebin[19] = 		{4, 4, 2, 2, 2, 
												2, 2, 2, 2, 2, 
												2, 2, 2, 2, 2, 
												2, 4, 4, 4};

// ***************** Pt binning for pPb, 5.023 TeV *************************************
Double_t fBinsPi0pPbPt[32] =					{0.0, 0.3, 0.4, 0.5, 0.6,
												0.7, 0.8, 1.0, 1.2, 1.4,
												1.6, 1.8, 2.0, 2.2, 2.4,
												2.6, 2.8, 3.0, 3.2, 3.4,
												3.6, 3.8, 4.0, 4.5, 5.0,
												5.5, 6.0, 7.0, 8.0, 10.0,
												12.0, 14.0};
//Double_t fBinsPi0pPbPt[10] = 					{0.0, 0.5, 1.0, 1.5, 2.0,
// 												3.0, 4.0, 6.0, 8.0, 14.0};
Double_t fBinsPi0pPbPtDCA[13] =					{0.0, 0.4, 0.5, 0.6, 0.7,
												0.8, 1.0, 1.4, 1.8, 2.4,
												4.0, 6.0, 14.0};
Double_t fBinsPi0pPbPt_Cent[25] =				{0.0, 0.4, 0.5, 0.6, 0.7,
												0.8, 1.0, 1.2, 1.4, 1.6,
												1.8, 2.0, 2.2, 2.4, 2.6,
												2.8, 3.0, 3.5, 4.0, 5.0,
												6.0, 8.0, 10.0, 12.0, 14.0};
Int_t fBinsPi0pPbPtRebin[31] =					{10, 8, 4, 2, 2,
												1, 1, 1, 1, 1,
												1, 1, 1, 1, 1,
												1, 1, 1, 1, 1,
												1, 1, 1, 2, 2,
												2, 4, 4, 5, 10,
												10};
//Int_t fBinsPi0pPbPtRebin[9] = 				{10, 4, 4, 4, 4, 
// 												4, 4, 4, 10};
Int_t fBinsPi0pPbPtRebin_Cent[24] =				{10, 4, 4, 4, 4,
												4, 4, 4, 4, 4, 
												4, 4, 4, 4, 4, 
												4, 4, 4, 4, 4, 
												4, 5, 10, 10}; 

Double_t fBinsEtapPbPt[14] = 					{0., 0.3, 0.5, 0.7, 0.9, 
												1.1, 1.4, 1.8, 2.2,
												3.0, 4., 5.,
												7.,   10.};
Double_t fBinsOmegapPbPt[8] = 					{0., 1., 2.0, 3.0, 4.0, 6.0, 8.0, 12.0};
Double_t fBinsEtapPbPt_Cent[15] = 				{0.,  0.4,  0.6,  0.8,  1.0,
												1.2, 1.4,  1.6,  2.0,  2.5,
												3.0, 4.,   6.,   8.,   25.};
Double_t fBinsOmegapPbPt_Cent[15] = 			{0.,  0.4,  0.6,  0.8,  1.0,
												1.2, 1.4,  1.6,  2.0,  2.5,
												3.0, 4.,   6.,   8.,   25.};
										
//Double_t fBinsEtapPbPt[9] = 					{0., 0.5, 1.0, 1.5, 2.0, 
// 												2.5, 4., 6., 8.};
Double_t fBinsEtapPbPtDCA[6] = 					{0., 0.5, 1.0, 2.0, 4., 8.};
Int_t fBinsEtapPbPtRebin[13] =	 				{10, 8, 8, 8, 
												6,  6,  4,  4,
												6,  8,  8, 8, 10};
Int_t fBinsOmegapPbPtRebin[7] =	 				{10, 8, 8, 6, 6, 8, 10};

Int_t fBinsEtapPbPtRebin_Cent[14] = 			{10, 10, 8, 8, 8, 
												5,  5,  4,  4, 4,
												4,  5,  8,  
												10};   
Int_t fBinsPi0EtaBinningpPbPtRebin[17] = 		{8, 2, 1, 1, 1,
												1, 1, 1, 1, 1,
												1, 1, 2, 2, 2,
												4, 4};
Double_t fBinsEtapPbPt3Body[15] = 				{0., 0.4, 0.6, 0.8, 1.0,
												1.2, 1.4, 1.6, 2.0, 2.5,
												3.0, 4.0, 6.0, 8.0, 25.};
//Double_t fBinsEtapPbPt[9] = 					{0., 0.5, 1.0, 1.5, 2.0,
// 												2.5, 4.0, 6.0, 8.0};
Int_t fBinsEtapPbPt3BodyRebin[14] = 			{5, 5, 5, 5, 5, 
												5, 5, 5, 5, 5,
												5, 5, 5, 5};
Double_t fBinsPi0pPbDirectPhotonPt[20] = 		{0.0, 0.4, 0.8, 1.0, 1.2, 
												1.4, 1.6, 1.8, 2.0, 2.3, 
												2.7, 3.1, 3.5, 4.0, 4.5,
												5.5, 6.5, 8.0, 11.0, 14.0};
Int_t fBinsPi0pPbDirectPhotonPtRebin[19] = 		{4, 4, 2, 2, 2, 
												2, 2, 2, 2, 2, 
												2, 2, 2, 2, 2, 
												2, 4, 4, 4};

