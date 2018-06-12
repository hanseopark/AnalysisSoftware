TString     fdate                                                       = "";
TString     fEnergyFlag                                                 = "";
TString     fCollisionSystem                                            = "";
TString     fTextMeasurement                                            = "";
TString     fDetectionProcess                                           = "";
TString     fCutSelection                                               = "";
TString     fCutSelectionRead                                           = "";
TString     fEventCutSelection                                          = "";
TString     fGammaCutSelection                                          = "";
TString     fClusterCutSelection                                        = "";
TString     fElectronCutSelection                                       = "";
TString     fMesonCutSelection                                          = "";

Double_t fMinPt     = 0.;
Double_t fMinPt1    = 0.8;
Double_t fMaxPt     = 20.;

Int_t rebinRPlots   = 2;
Int_t rebinZPlots   = 4;
Int_t rebinPtPlots  = 4;

Double_t rMinGas    = 95.;
Double_t rMaxGas    = 145.;

const int nBinsR = 12;
// Double_t arrayRBins[14] = {0., 1.5, 3.5, 5.7, 8.6, 13., 21.,  33.5, 41., 55., 72., 90.,  150., 180};
Double_t arrayRBins[13] = {0., 1.5, 5., 8.6, 13., 21., 33.5, 41., 55., 72., 95., 145., 180};
//                      = {0., 1.5, 5., 8.6, 13., 21., 33.5, 41., 55., 72., 90., 150., 180};

TString arrayRangesRBins[12] =
                       {"0 cm < R < 1.5 cm",     //0
                        "1.5 cm < R < 5. cm",    //1
                        "5. cm < R < 8.6 cm",    //2
                        "8.6 cm < R < 13 cm",    //3
                        "13 cm < R < 21 cm",     //4
                        "21 cm < R < 33.5 cm",   //5
                        "33.5 cm < R < 41 cm",   //6
                        "41 cm < R < 55 cm",     //7
                        "55 cm < R < 72 cm",     //8
                        "72 cm < R < 90 cm",     //9
                        "90 cm < R < 150 cm",    //10
                        "150 cm < R < 180 cm"    //11
                       };

TString arrayNamesRBins[12] =
                       {"Vertex",                                       //0
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

// const int 	nBinsR = 				10;
// Float_t         arrayRBins[11] =                {0.,2.,3.5,13.,21.,33.5,55.,72.,90.,150.,180};
// TString         arrayNamesRangesRBins[11]=      {"0 cm < R < 2 cm",       //0
// 						     "2. cm < R < 3.5 cm",    //1
// 						     "3.5 cm < R < 13 cm",    //2
// 						     "13 cm < R < 21. cm",   //3
// 						     "21 cm < R < 33.5 cm",   //4
// 						     "33.5 cm < R < 55 cm",   //5
// 						     "55 cm < R < 72 cm",     //6
// 						     "72 cm < R < 90 cm",     //7
// 						     "90 cm < R < 150 cm",    //8
// 						     "150 cm < R < 180 cm"};  //9

Double_t arrayZBins[13] = {10.,15.,20.,40.,40.,60.,60.,80.,120.,120.,200.,200.,200};

TLine *lineRLimits[13];
TLine *lineRLimitsPurity[13];

Double_t doubleLatexNamingBinsX = 0.16;
Double_t doubleLatexNamingBinsRatioX = 0.11;
Double_t doubleLatexNamingBinsX2 = 0.24;
Double_t doubleLatexNamingBinsY = 0.9;
Double_t doubleLatexNamingBinsY2 = 0.86;
Double_t doubleLatexNamingBinsY3 = 0.82;
Size_t sizeTextNameBins = 0.05;

Double_t doubleLatexNamingCutX=0.16;
Double_t doubleLatexNamingCutY=0.20;

TString 	nameHistoRatioPhiInRMC;
TString 	nameHistoRatioZInRMC;

TH1F*   histoRDataScaledToGas;
TH1F*   histoRFullPtDataScaledToGas;
TH1F*   histoRMCScaledToGas;
TH1F*   histoRFullPtMCScaledToGas;
TH1F*   histoRDataScaledToGasRebin;
TH1F*   histoRMCScaledToGasRebin;
TH1F*   histoDataMCRatioRScaledToGas;
TH1F*   histoRFullPtDataScaledToGasRebin;
TH1F*   histoRFullPtMCScaledToGasRebin;
TH1F*   histoDataMCRatioRFullPtScaledToGas;

TH1F*   histoRinPtBinDataScaledToGas01;
TH1F*   histoRinPtBinMCScaledToGas01;
TH1F*   histoRinPtBinDataScaledToGasRebin01;
TH1F*   histoRinPtBinMCScaledToGasRebin01;
TH1F*   histoDataMCRatioRinPtBinScaledToGas01; 

TH1F*   histoRinPtBinDataScaledToGas02;
TH1F*   histoRinPtBinMCScaledToGas02;
TH1F*   histoRinPtBinDataScaledToGasRebin02;
TH1F*   histoRinPtBinMCScaledToGasRebin02;
TH1F*   histoDataMCRatioRinPtBinScaledToGas02; 

TH1F*   histoRinPtBinDataScaledToGas03;
TH1F*   histoRinPtBinMCScaledToGas03;
TH1F*   histoRinPtBinDataScaledToGasRebin03;
TH1F*   histoRinPtBinMCScaledToGasRebin03;
TH1F*   histoDataMCRatioRinPtBinScaledToGas03; 

TH1F*	histoDataMCRatioR;
TH1F*	histoMidPtDataMCRatioR;
TH1F*	histoPurityR;
TH1F*	histoPurityPt;
TH1F*	histoPurityPt5cm;

TH1D*   histoPhiInRData[nBinsR];
TH1D*   histoPhiInRMC[nBinsR];
TH1D*   histoPhiInRTrueMC[nBinsR];
TH1D*	histoDataMCRatioPhiInR[nBinsR];

TH1D*   histoIntegralGasStatErrorData[nBinsR];
TH1D*   histoIntegralGasStatErrorData01[nBinsR];
TH1D*   histoIntegralGasStatErrorData02[nBinsR];
TH1D*   histoIntegralGasStatErrorData03[nBinsR];
TH1D*   histoIntegralGasStatErrorMC[nBinsR];
TH1D*   histoIntegralGasStatErrorMC01[nBinsR];
TH1D*   histoIntegralGasStatErrorMC02[nBinsR];
TH1D*   histoIntegralGasStatErrorMC03[nBinsR];

TH1D*   histoIntegralGasData;
TH1D*   histoIntegralGasData01;
TH1D*   histoIntegralGasData02;
TH1D*   histoIntegralGasData03;
TH1D*   histoIntegralGasMC;
TH1D*   histoIntegralGasMC01;
TH1D*   histoIntegralGasMC02;
TH1D*   histoIntegralGasMC03;

TH1F*   histoRDataRebin;
TH1F*   histoRinPtBinDataRebin01;
TH1F*   histoRinPtBinDataRebin02;
TH1F*   histoRinPtBinDataRebin03;

TH1F*   histoRinPtBinMCRebin01;
TH1F*   histoRinPtBinMCRebin02;
TH1F*   histoRinPtBinMCRebin03;

TH1D*   histoZInRData[nBinsR];
TH1D*   histoZInRMC[nBinsR];
TH1D*   histoZInRTrueMC[nBinsR];
TH1D*	histoDataMCRatioZInR[nBinsR];

Int_t rebinPhiPlots=3;

Double_t markerSize =1.05;

Color_t colorData = kBlack;
Color_t colorMC = kRed;
Color_t colorTrueMC = kMagenta+2;
Color_t colorTrueCombMC = kGreen+2;

Color_t colorComparisonMC = kBlue-4;



//__________________________________________ Plotting all Invariant Mass bins _______________________________________________
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
