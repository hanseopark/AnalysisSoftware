TString     fdate                 = "";
TString     fEnergyFlag           = "";
TString     fCollisionSystem      = "";
TString     fTextMeasurement      = "";
TString     fDetectionProcess     = "";
TString     fCutSelection         = "";
TString     fCutSelectionRead     = "";
TString     fEventCutSelection    = "";
TString     fGammaCutSelection    = "";
TString     fClusterCutSelection  = "";
TString     fElectronCutSelection = "";
TString     fMesonCutSelection    = "";
TString     nameMainDir           = "GammaConvMaterial";

Double_t markerSize               = 1.2;
Color_t colorData                 = kBlack;
Color_t colorMC                   = kRed;
Color_t colorTrueMC               = kMagenta+2;
Color_t colorTrueCombMC           = kGreen+2;
Color_t colorComparisonMC         = kBlue-4;

Double_t fMinPt                   = 0.;
Double_t fMinPt1                  = 0.8;
Double_t fMaxPt                   = 20.;

Int_t rebinRPlots                 = 1;
Int_t rebinZPlots                 = 4;
Int_t rebinPtPlots                = 2;
Int_t rebinPhiPlots               = 3;

Double_t rMinGas                  = 95.;
Double_t rMaxGas                  = 145.;

const int nBinsPt=69;
Double_t arrayPtBins[nBinsPt]; 
const int nBinsPtTwo=31;
Double_t arrayPtBinsTwo[nBinsPtTwo]; 

//const int nBinsPtFine             = 15;
//Double_t projPtBinsFine[nBinsPtFine]     = {0.05, 0.1, 0.15, 0.2, 0.25,
//                                            0.3,  0.35, 0.4, 0.45, 0.5,
//                                            0.55, 0.6, 0.65, 0.7, 0.75};

Double_t projRBins[4]      = {0., 5., 35., 180.};
TH1F *histoPtinRBinData[3] = {NULL, NULL, NULL};
TH1F *histoPtinRBinMC[3]  = {NULL, NULL, NULL};
//AM . Changed from 0.15 to 0.1 to be able to use the secondary subtracted
Double_t projPtBins[6]     = {0.1, 0.3, 0.4, 0.5, 0.6, 0.7};
TH1F *histoRinPtBinData[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
TH1F *histoRinPtBinMC[6] = {NULL, NULL, NULL, NULL, NULL, NULL};

const int nBinsPtFine             = 8;
Double_t projPtBinsFine[nBinsPtFine]     = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};


TH1F *histoRinPtBinDataFine[nBinsPtFine] = {NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL}; 
                                           
TH1F *histoRinPtBinMCFine[nBinsPtFine]   = {NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL};
Double_t nconvInRangeDataFine[nBinsPtFine];
Double_t dataStatErrorGasFine[nBinsPtFine];
Double_t dataStatRelErrorGasFine[nBinsPtFine];

TH1D * histoIntegralGasDataFine[nBinsPtFine]   = {NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL};
TH1F * histoRinPtBinDataRebinFine[nBinsPtFine] = {NULL, NULL, NULL, NULL, NULL,NULL, NULL, NULL};
TH1F * histoRinPtBinDataScaledToGasRebinFine[nBinsPtFine] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

Double_t nconvInRangeMCFine[nBinsPtFine];
Double_t dataStatErrorGasMCFine[nBinsPtFine];
Double_t dataStatRelErrorGasMCFine[nBinsPtFine];
TH1D * histoIntegralGasMCFine[nBinsPtFine]= {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
TH1F * histoRinPtBinMCRebinFine[nBinsPtFine]= {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
TH1F * histoRinPtBinMCScaledToGasRebinFine[nBinsPtFine]= {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
TH1F * histoDataMCRatioRinPtBinScaledToGasFine[nBinsPtFine]= {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};


const int nBinsR                  = 12;
Double_t arrayRBins[13]           = {0., 1.5, 5., 8.5, 13., 21., 33.5, 41., 55., 72., 95., 145., 180};
TString arrayRangesRBins[12]      = { "0 cm < R #leq 1.5 cm",     //0
                                      "1.5 cm < R #leq 5. cm",    //1
                                      "5. cm < R #leq 8.5 cm",    //2
                                      "8.5 cm < R #leq 13 cm",    //3
                                      "13 cm < R #leq 21 cm",     //4
                                      "21 cm < R #leq 33.5 cm",   //5
                                      "33.5 cm < R #leq 41 cm",   //6
                                      "41 cm < R #leq 55 cm",     //7
                                      "55 cm < R #leq 72 cm",     //8
                                      "72 cm < R #leq 95 cm",     //9
                                      "95 cm < R #leq 145 cm",    //10
                                      "145 cm < R #leq 180 cm"    //11
                                    };

TString arrayNamesRBins[12]       = { "Vertex",                                       //0
                                      "BeamPipe+SPD 1",                               //1
                                      "SPD 2",                                        //2
                                      "Thermal shield/Support between SPD/SDD",       //3
                                      "SDD 1 +Thermal shield",                        //4
                                      "SDD 2 +Thermal shield",                        //5
                                      "SSD 1",                                        //6
                                      "SSD 2",                                        //7
                                      "Air + TPC in. cont. vessel + CO_{2}",          //8
                                      "CO_{2} + TPC in. field cage vessel+TPC rods ", //9
                                      "Ne: CO_{2}: N_{2}",                            //10
                                      "Ne: CO_{2}: N_{2}"                             //11
                                    };

Double_t fMaterialWeightsForSecEffCor[12]={1.0377,1.19547,1.24405,1.22941,1.07585,1.12396, 1.09077, 1.03798,0.948079,0.981036, 1.,0.964843};
//1.0313, 1.18782, 1.23548, 1.22107 ,1.07171, 1.11886, 1.08748, 1.0367, 0.950803, 0.982077, 1., 0.966804};
//{1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};

TH1F * histoMCSecGammaPtBySource[4]={NULL,NULL,NULL,NULL};
TH1F * histoMCSecConvGammaPtBySource[4]={NULL,NULL,NULL,NULL};
TH1F * histoMCSecConvGammaPtEachRBinBySource[4][nBinsR];

TH1F * histoMCTrueRecSecGammaPtEachRBinBySource[4][nBinsR];
TH1F * histoTrueEffSecEachRBinBySource[4][nBinsR];
TH1F * histoTrueRecEffSecEachRBinBySource[4][nBinsR];
TH1F * histoConvProbSecEachRBinBySource[4][nBinsR];
TH1F * histoConvProbSecBySource[4];
TH1F * histoFracSecPtEachRBinBySource[4][nBinsR];
TH1F * histoPtEachRBinDataSecYieldFromSecFracBySourceRaw[4][nBinsR];
TH1F * histoPtEachRBinDataSecSubtractedUsingCocktail[nBinsR];
TH1F * histoDataFromCocktailSecConvGammaPtEachRBinBySourceRaw[4][nBinsR];

Int_t       nBins                                           = 0;
Double_t    xMin                                            = 0.;
Double_t    xMax                                            = 20;


// pT distribution of data reconstructed in each R Bin
TH1F *histoPtEachRBinData[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};


// pT distribution of MC reconstructed in each R Bin, validated true and validated true primary
TH1F *histoPtEachRBinMC[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
TH1F *histoPtTrueMCEachRBin[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};    
TH1F *histoPtTruePrimMCEachRBin[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

//  pT distribution of secondary from K0S in each R bin
TH1F * histoPtTrueSecEachRBinFromK0S[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
TH1F * histoPtTrueSecEachRBinAllSources[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
TH1F * histoPurityPtEachRBin[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}; 
TH1F * histoPurityPrimPtEachRBin[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}; 
TH1F * histoFracSecPtEachRBin[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}; 
TH1F * histoEfficiencySecEachRBinAllSources[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}; 

TH1F *histoPtEachRBinDataSecSubtracted[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
TH1F *histoPtEachRBinMCSecSubtracted[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

TH1F *histoPtEachRBinDataSecYieldFromSecFrac[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
TH1F *histoPtEachRBinMCSecYieldFromSecFrac[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};


const int nBinsPtMin= 8;
Double_t  arrayBinsPtMin[nBinsPtMin+1];
TH1F * histoWeightsEachRPtMin[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
TH1F * histoWeightsEachRPtMinSecSub[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
TH1F * histoWeightsEachRPtMinSecSubUsingCocktail[nBinsR] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

Double_t nConvInRangeFromPtMinSecSubtractedData[nBinsR][nBinsPtFine];
Double_t nConvInRangeFromPtMinSecSubtractedDataRelErr[nBinsR][nBinsPtFine];

Double_t nConvInRangeFromPtMinSecSubtractedDataUsingCocktail[nBinsR][nBinsPtFine];
Double_t nConvInRangeFromPtMinSecSubtractedDataUsingCocktailRelErr[nBinsR][nBinsPtFine];

Double_t nConvInRangeFromPtMinSecSubtractedMC[nBinsR][nBinsPtFine];
Double_t nConvInRangeFromPtMinSecSubtractedMCRelErr[nBinsR][nBinsPtFine];
Double_t nConvInRangeFromPtMinSecSubtractedDataToGas[nBinsR][nBinsPtFine];
Double_t nConvInRangeFromPtMinSecSubtractedDataUsingCocktailToGas[nBinsR][nBinsPtFine];
Double_t nConvInRangeFromPtMinSecSubtractedMCToGas[nBinsR][nBinsPtFine];
Double_t weightInRangeFromPtMinSecSubtracted[nBinsR][nBinsPtFine];
Double_t weightInRangeFromPtMinSecSubtractedUsingCocktail[nBinsR][nBinsPtFine];

Double_t nConvInRangeFromPtMinSecSubtractedDataToGasRelErr[nBinsR][nBinsPtFine];
Double_t nConvInRangeFromPtMinSecSubtractedDataUsingCocktailToGasRelErr[nBinsR][nBinsPtFine];
Double_t nConvInRangeFromPtMinSecSubtractedMCToGasRelErr[nBinsR][nBinsPtFine];
Double_t weightInRangeFromPtMinSecSubtractedRelErr[nBinsR][nBinsPtFine];
Double_t weightInRangeFromPtMinSecSubtractedUsingCocktailRelErr[nBinsR][nBinsPtFine];

TH1F *histoPtinRBinTrueMC[3]     = {NULL, NULL, NULL};
TH1F *histoPtinRBinTruePrimMC[3] = {NULL, NULL, NULL};
TH1F *histoPtinRBinTrueSecMC[3]  = {NULL, NULL, NULL};

TH1F *histoRinPtBinTrueMC[6]     = {NULL, NULL, NULL, NULL, NULL, NULL};
TH1F *histoRinPtBinTruePrimMC[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
TH1F *histoRinPtBinTrueSecMC[6]  = {NULL, NULL, NULL, NULL, NULL, NULL};


const int nBinsZ                  = 12;
Double_t arrayZBins[13]           = {10.,15.,20.,40.,40.,60.,60.,80.,120.,120.,200.,200.,200};

Double_t eps = 0.001;
Size_t sizeTextNameBins              = 0.05;

TString nameHistoRatioPhiInRMC;
TString nameHistoRatioZInRMC;

TH1F*   histoRDataScaledToGas;
TH1F*   histoRFullPtDataScaledToGas;
TH1F*   histoRMCScaledToGas;
TH1F*   histoRFullPtMCScaledToGas;
TH1F*   histoRDataScaledToGasRebin;
TH1F*   histoRMCScaledToGasRebin;
TH1F*   histoDataMCRatioRScaledToGas;
TH1F*   histoDataMCRatioRScaledToGasSecSub;
TH1F*   histoRFullPtDataScaledToGasRebin;
TH1F*   histoRFullPtMCScaledToGasRebin;
TH1F*   histoDataMCRatioRFullPtScaledToGas;

TH1F*   histoRinPtBinDataScaledToGasPtBin1;
TH1F*   histoRinPtBinMCScaledToGasPtBin1;
TH1F*   histoRinPtBinDataScaledToGasPtBin1Rebin;
TH1F*   histoRinPtBinMCScaledToGasPtBin1Rebin;
TH1F*   histoDataMCRatioRinPtBinScaledToGasPtBin1;
TH1F*   histoDataMCRatioRScaledToGasSecSubPtBin1;

TH1F*   histoRinPtBinDataScaledToGasPtBin2;
TH1F*   histoRinPtBinMCScaledToGasPtBin2;
TH1F*   histoRinPtBinDataScaledToGasPtBin2Rebin;
TH1F*   histoRinPtBinMCScaledToGasPtBin2Rebin;
TH1F*   histoDataMCRatioRinPtBinScaledToGasPtBin2;
TH1F*   histoDataMCRatioRScaledToGasSecSubPtBin2;

TH1F*   histoRinPtBinDataScaledToGasPtBin3;
TH1F*   histoRinPtBinMCScaledToGasPtBin3;
TH1F*   histoRinPtBinDataScaledToGasPtBin3Rebin;
TH1F*   histoRinPtBinMCScaledToGasPtBin3Rebin;
TH1F*   histoDataMCRatioRinPtBinScaledToGasPtBin3;
TH1F*   histoDataMCRatioRScaledToGasSecSubPtBin3;

TH1F*	histoDataMCRatioR;
TH1F*	histoDataMCRatioRRebin;
TH1F*	histoMidPtDataMCRatioR;
TH1F*	histoPurityR;
TH1F*	histoPurityPrimR;
TH1F*	histoPurityPrimRPtBin3;
TH1F*	histoPuritySecR;
TH1F*	histoPuritySecRPtBin3;
TH1F*	histoPurityPt;
TH1F*	histoPurityPt5cm;
TH1F*	histoPurityPrimPt5cm;
TH1F*	histoPuritySecPt5cm;
TH1F*	histoEffiR;
TH1F*	histoEffiPt;
TH1F*	histoEffiPhi;
TH1F*	histoEffiEta;

TH1D*   histoPhiInRData[nBinsR];
TH1D*   histoPhiInRMC[nBinsR];
TH1D*   histoPhiInRTrueMC[nBinsR];
TH1D*	histoDataMCRatioPhiInR[nBinsR];

TH1D*   histoIntegralGasStatErrorDataFullRange[nBinsR];
TH1D*   histoIntegralGasStatErrorDataPtBin1[nBinsR];
TH1D*   histoIntegralGasStatErrorDataPtBin2[nBinsR];
TH1D*   histoIntegralGasStatErrorDataPtBin3[nBinsR];
TH1D*   histoIntegralGasStatErrorMC[nBinsR];
TH1D*   histoIntegralGasStatErrorMCPtBin1[nBinsR];
TH1D*   histoIntegralGasStatErrorMCPtBin2[nBinsR];
TH1D*   histoIntegralGasStatErrorMCPtBin3[nBinsR];

TH1D*   histoIntegralGasDataFullRange;
TH1D*   histoIntegralGasDataPtBin1;
TH1D*   histoIntegralGasDataPtBin2;
TH1D*   histoIntegralGasDataPtBin3;
TH1D*   histoIntegralGasMCFullRange;
TH1D*   histoIntegralGasMCPtBin1;
TH1D*   histoIntegralGasMCPtBin2;
TH1D*   histoIntegralGasMCPtBin3;

TH1F*   histoRDataRebin;
TH1F*   histoRinPtBinDataPtBin1Rebin;
TH1F*   histoRinPtBinDataPtBin2Rebin;
TH1F*   histoRinPtBinDataPtBin3Rebin;

TH1F*   histoRMCRebin;
TH1F*   histoRinPtBinMCPtBin1Rebin;
TH1F*   histoRinPtBinMCPtBin2Rebin;
TH1F*   histoRinPtBinMCPtBin3Rebin;

TH1D*   histoZInRData[nBinsR];
TH1D*   histoZInRMC[nBinsR];
TH1D*   histoZInRTrueMC[nBinsR];
TH1D*	histoDataMCRatioZInR[nBinsR];

Double_t arrayX[2];
Double_t arrayY[3];
Double_t relX[3];
Double_t relY[3];
Int_t textSizeLabels = 30;
Double_t textsizeLabelsDown = 0;
Double_t textsizeFacDown = 0;
Double_t textsizeLabelsUp = 0;
Double_t textsizeFacUp = 0;

Double_t arrayXRdistrib[3];
Double_t arrayYRdistrib[3];
Double_t relXRdistrib[3];
Double_t relYRdistrib[3];

Double_t arrayXRConv2Pad[2];
Double_t arrayYRConv2Pad[3];
Double_t relXRConv2Pad[3];
Double_t relYRConv2Pad[3];
Double_t textsizeLabelsDownRConv2Pad = 0;
Double_t textsizeFacDownRConv2Pad = 0;
Double_t textsizeLabelsUpRConv2Pad = 0;
Double_t textsizeFacUpRConv2Pad = 0;

Double_t arrayXpurity[2];
Double_t arrayYpurity[3];
Double_t relXpurity[3];
Double_t relYpurity[3];
Double_t textsizeLabelsDownpurity = 0;
Double_t textsizeFacDownpurity = 0;
Double_t textsizeLabelsUppurity = 0;
Double_t textsizeFacUppurity = 0;

Double_t arrayXPhotonChar[3];
Double_t arrayYPhotonChar[3];
Double_t relXPhotonChar[3];
Double_t relYPhotonChar[3];

Double_t arrayXPhiInRbins[2];
Double_t arrayYPhiInRbins[3];
Double_t relXPhiInRbins[3];
Double_t relYPhiInRbins[3];

Double_t arrayXZInRbins[2];
Double_t arrayYZInRbins[3];
Double_t relXZInRbins[3];
Double_t relYZInRbins[3];

Double_t EtaRange[2] = {-1.5,1.5};
Double_t PhiRange[2] = {0.,6.28};
Double_t PtRange[2]  = {0.,15.};
Double_t RRange[2]   = {0.,180.};
Double_t minYRatio = 0.9;
Double_t maxYRatio = 1.5;

Color_t colorOnFly = kRed+1;
Color_t colorOffline = kBlue+1;

Color_t color[12] = { kBlack, kAzure, kGreen+2, kOrange+2, kCyan+2,
                      kYellow+2, kViolet-3, kSpring+10, kRed+1, kMagenta-8,
                      kGray, kGray+3};

//TH1F* histoMCSecGammaPtBySource[4];
//TH1F* histoMCSecConvGammaPtBySource[4][nBinsR];
TH1F* histoMCTrueRecSecGammaPtBySource[4][nBinsR];
//TH1F * histoTrueRecEffSecEachRBinBySource[4][nBinsR];
//TH1F * histoConvProbSecEachRBinBySource[4][nBinsR];

//TH1F* histoMCAllGammaPt = NULL;
//TH2F * histoConvRPtMC = { NULL, NULL, NULL, NULL };


TH1F * histoRecEffPrimaryEachRBin[nBinsR];
TH1F * histoConvProbPrimaryEachRBin[nBinsR];



// Taken from ExtractGammaSignalV2.h
TFile*      fFileCocktailInput                                          = NULL;
TH1F*       fHistoSecGammaCocktailFromXPt[4]                                              = { NULL, NULL, NULL, NULL };
TH1F*       fHistoSecGammaCocktailFromXPtOrBin[4]                                         = { NULL, NULL, NULL, NULL };
TString     fSecondaries[4]                                         = {"K0s", "K0l", "Lambda", "Rest"};

Bool_t   LoadSecondariesFromCocktailFile    (   TString,
						TString                                         );



void PlotRDistribRunwise(   TH1D** fHistoData,
                            TH1D** fHistoMC,
                            TString namePlot,
                            TString nameCanvas,
                            TString namePad,
                            Int_t fRowPlot,
                            Int_t fColumnPlot,
                            TString *runNumberName,
                            Int_t startNrun,
                            Int_t endNrun
                        ){
    TGaxis::SetMaxDigits(3);

    TCanvas * canvasDataSpectra     = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    canvasDataSpectra->SetTopMargin(0.0);
    canvasDataSpectra->SetBottomMargin(0.0);
    canvasDataSpectra->SetRightMargin(0.0);
    canvasDataSpectra->SetLeftMargin(0.0);
    canvasDataSpectra->SetLogy(1);

    TPad * padDataSpectra           = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataSpectra->SetFillColor(0);
    padDataSpectra->GetFrame()->SetFillColor(0);
    padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(1);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();

    Int_t place = 0;
    for(Int_t i = startNrun; i<endNrun; i++){

        place                       = place + 1; //give the right place in the page
        if (place == fColumnPlot){
            i--;
            padDataSpectra->cd(place);

            TString textAlice       = "";

            Double_t nPixels        = 13;
            Double_t textHeight     = 0.08;
            if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
                textHeight          = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
            } else {
                textHeight          = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
            }

            Double_t startTextX     = 0.10;
            Double_t startTextY     = 0.9;
            Double_t differenceText = textHeight*1.25;

            TLatex *alice           = new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
//             TLatex *energy          = new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);

            alice->SetNDC();
            alice->SetTextColor(1);
            alice->SetTextSize(textHeight*1.3);
            alice->Draw();

//             energy->SetNDC();
//             energy->SetTextColor(1);
//             energy->SetTextSize(textHeight);
//             energy->Draw();

            TLegend* legendData     = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+2.)*differenceText);
            legendData->SetTextSize(textHeight);
            legendData->SetTextFont(62);
            legendData->SetFillColor(0);
            legendData->SetFillStyle(0);
            legendData->SetLineWidth(0);
            legendData->SetLineColor(0);
            legendData->SetMargin(0.15);
            Size_t markersize       = fHistoData[i]->GetMarkerSize();
            fHistoData[i]->SetMarkerSize(3*markersize);
            Size_t linesize         = fHistoData[i]->GetLineWidth();
            legendData->AddEntry(fHistoData[i],"Data","l");
            legendData->AddEntry(fHistoMC[i],"MC","l");
            legendData->Draw();

        } else {

            padDataSpectra->cd(place);
            padDataSpectra->cd(place)->SetLogy(1);
            padDataSpectra->cd(place)->SetTopMargin(0.08);
            padDataSpectra->cd(place)->SetBottomMargin(0.08);
            padDataSpectra->cd(place)->SetRightMargin(0.015);
            int remaining           = (place-1)%fColumnPlot;
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.08);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);
            TString titlePt = Form("Run number %s",runNumberName[i].Data());
//             cout << titlePt << endl;

            TH2F * histoDummyR = new TH2F("histoDummyR","histoDummyR",1000,0.,180.,1000,1.e-6,0.1);
            SetStyleHistoTH2ForGraphs(histoDummyR, "R (cm)",Form("%s Counts",titlePt.Data()), 0.035,0.04,0.035,0.04,1.,1.);
    		histoDummyR->GetYaxis()->SetRangeUser(1.e-6, 1.e-2);
            histoDummyR->DrawCopy();

            DrawGammaSetMarker(fHistoData[i], 20, 1.5, colorData, colorData);
            DrawGammaSetMarker(fHistoMC[i], 20, 1.5, colorMC, colorMC);
            fHistoMC[i]->Draw("same,hist");
            fHistoData[i]->Draw("same,hist");

        }
    }
    canvasDataSpectra->Print(namePlot.Data());
    delete padDataSpectra;
    delete canvasDataSpectra;
}

//*********************************************************************************************************
//*********** Convert secondary spectrum from cocktail to raw (conversion method)            **************
// - multiply with secondary conv. prob.
// - use unfolding for resolution correction if response matrix given
// - multiply with secondary reco. eff.
// - scale with nEvts for same event norm as data spectrum
//*********************************************************************************************************
Bool_t ConvertCocktailSecondaryToRaw(TH1F*       histoGammaSec,
                                     TH1F*       histoConvProb,
                                     TH1F*       histoRecoEff,
                                     TH2F*       responseMatrix,
                                     Double_t    nEvt,
                                     Bool_t      useResponseMatrix,
                                     Int_t       nIterationsUnfolding = 5
                                     ){
    // multiply with conv. prob.
    if (!histoConvProb) return kFALSE;
    //histoGammaSec->Multiply(histoConvProb);
    // multiply with reco. eff. (in MC pT if unfolding is used)
    if (!histoRecoEff) return kFALSE;
    histoGammaSec->Multiply(histoRecoEff);
    // correct for resolution, if response matrix given (otherwise included in reco. eff.)
    //    if (useResponseMatrix) {
    //        if (!responseMatrix) return kFALSE;
    //        cout << "will use unfolding for " << histoGammaSec->GetName() << endl;
    //        RooUnfoldResponse response(0,0,responseMatrix);
    //        RooUnfoldBayes unfold_SpectrumCocktail (&response, histoGammaSec, nIterationsUnfolding);
    //        histoGammaSec = (TH1D*)unfold_SpectrumCocktail.Hreco();
    //    }
    // scale with nEvt from data

    histoGammaSec->Scale(nEvt);
    
    return kTRUE;
}
