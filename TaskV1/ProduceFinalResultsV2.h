// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion
TString 		collisionSystem;
TString			detSystem;
const char* 	fileNameSysErrEta;
Bool_t 			thesis 		=	kFALSE;
Bool_t 			conference 	=	kFALSE;
Double_t		minPtForFits=	0.2;
Double_t		minPtForFitsEta = 		0.4;
TString 		rapitdityCutNumber;
Double_t 		deltaEta;
TString 		rapidityRange;
Bool_t 			plotPythiaPhojetInRatio = kFALSE;
Double_t    	xSection                = 0; 
TString 		prefix2;
Bool_t	 		kLevy;
Bool_t 			kHag;
Double_t 		scaling = 					1./(2.*TMath::Pi());
Double_t 		mesonMassExpectPi0 = 			0;
Double_t 		mesonMassExpectEta = 			0;

TString 		fileNamePi0ch;
TFile* 		filePi0;
TH1D*		histoCorrectedYieldPi0;
TH1D*		histoUncorrectedYieldPi0;
TH1D*		histoFWHMMesonPi0;
TH1D*		histoMassMesonPi0;
TH1D*		histoAccPi0;
TH1D*		histoTrueEffPtPi0;
TH1D*		histoTrueFWHMMesonPi0;
TH1D*		histoTrueMassMesonPi0;
TH1D*    	histoMCInputPi0;
TH1F*		histoEventQualtityPi0;
TH1D* 		histoMassMesonPi0MinusExp;
TH1D* 		histoTrueMassMesonPi0MinusExp;
TH1D* 		histoMesonSignalFullPtInvMass;
TH1D* 		histoFWHMMesonPi0MeV;
TH1D*		histoTrueFWHMMesonPi0MeV;

TLatex * 		textInvMass;
TString 		date;

Float_t 		nEvt;
TFile* 		fileEta;
TH1D *		histoCorrectedYieldEta;
TH1D *		histoUnCorrectedYieldEta;
TH1D *		histoFWHMMesonEta;
TH1D *		histoMassMesonEta;
TH1D *		histoAccEta;
TH1D *		histoTrueEffPtEta;
TH1D *		histoTrueFWHMMesonEta;
TH1D *		histoTrueMassMesonEta;
TH1D*       histoMCInputEta;
TH1D * 		histoMassMesonEtaMinusExp;
TH1D * 		histoTrueMassMesonEtaMinusExp;
TH1D* 		histoFWHMMesonEtaMeV;
TH1D*		histoTrueFWHMMesonEtaMeV;

TH1D *		histoEtaToPi0Phojet;
TH1D *		histoEtaToPi0Pythia;
TH1D * 		histoRatioEtaPi0;

//Ratio as TGraphAsymError with +20 - 40
Double_t 		ratioXValue[50];
Double_t 		ratioYValue[50];
Double_t 		ratioXError[50];
Double_t 		ratioSysUpError[50];
Double_t 		ratioSysDownError[50];

const char *	nameFinalResDat;
fstream 		fileFinalResults;

//************ Fitting variables *****************
TF1 *		fitPtLevy;
TF1 *		fitPtHagedorn;
TF1 *		fitPtTCM;
TF1 *		fitPtBoltzmann;
TF1 *		fitPtExp;
TF1 *		fitPtPowerlaw;
TF1 *		fitPtModPowerlaw;
TH1D *		histoRatioFitLevy;
TH1D *		histoRatioFitHag;
TH1D *		histoRatioFitTCM;
TH1D *		histoRatioFitBoltz;
TH1D *		histoRatioFitExp;
TH1D *		histoRatioFitPow;
TH1D *		histoRatioFitModPow;
//*************** Variables Systematic errors *******************
TGraphAsymmErrors* 	graphCorrectedYieldPi0SysErr;
TGraphAsymmErrors* 	graphCorrectedYieldPi0SysErrA;
TGraphAsymmErrors* 	graphCorrectedYieldPi0StatPlusSys;

TH1D*			histoInvCrossSectionPi0;
TGraphAsymmErrors* 	graphInvCrossSectionSysPi0;
TGraphAsymmErrors* 	graphInvCrossSectionSysAPi0;

TGraphAsymmErrors* 	graphCorrectedYieldEtaSysErr;
TGraphAsymmErrors* 	graphCorrectedYieldEtaSysErrA;
TGraphAsymmErrors* 	graphCorrectedYieldEtaStatPlusSys;
TH1D	*			histoInvCrossSectionEta;
TGraphAsymmErrors* 	graphInvCrossSectionSysEta;
TGraphAsymmErrors* 	graphInvCrossSectionSysAEta;

TGraphAsymmErrors* 	graphRatioFitLevySysErr;
TGraphAsymmErrors* 	graphRatioFitTCMSysErr;
TGraphAsymmErrors* 	graphRatioFitPowSysErr;
TGraphAsymmErrors* 	graphRatioFitModPowSysErr;

ifstream 		fileSysErrPi0;
Int_t 		nPointsPi0 = 			0;
Double_t 		relSystErrorPi0Up[50];
Double_t 		relSystErrorPi0Down[50];
Double_t 		relSystErrorWOMaterialPi0Up[50];
Double_t 		relSystErrorWOMaterialPi0Down[50];
Double_t 		systErrorPi0Up[50];
Double_t 		systErrorPi0Down[50];
Double_t 		xValueCorrYieldPi0[50];
Double_t 		yValueCorrYieldPi0[50];
Double_t 		yErrorCorrYieldPi0[50];
Double_t 		xErrorCorrYieldPi0[50];
Double_t 		yErrorCombinedPi0Up[50];
Double_t 		yErrorCombinedPi0Down[50];
Double_t 		sysErrPi0RatioFitConvUp[50];
Double_t 		sysErrPi0RatioFitConvDown[50];
Double_t 		yValuePi0RatioFitConv[50];
ifstream 		fileSysErrEta;
Int_t 			nPointsEta = 0;
Double_t 		relSystErrorEtaUp[50];
Double_t 		relSystErrorEtaDown[50];
Double_t 		relSystErrorWOMaterialEtaUp[50];
Double_t 		relSystErrorWOMaterialEtaDown[50];
Double_t 		systErrorEtaUp[50];
Double_t 		systErrorEtaDown[50];
Double_t 		xValueCorrYieldEta[50];
Double_t 		yValueCorrYieldEta[50];
Double_t 		yErrorCorrYieldEta[50];
Double_t 		xErrorCorrYieldEta[50];
Double_t 		yErrorCombinedEtaUp[50];
Double_t 		yErrorCombinedEtaDown[50];
Double_t 		sysErrEtaRatioFitConvUp[50];
Double_t 		sysErrEtaRatioFitConvDown[50];
Double_t 		yValueEtaRatioFitConv[50];

TF1 *		fitPtLevySysErr;
TF1 *		fitPtHagedornSysErr;
TF1 *		fitPtBoltzmannSysErr;
TF1 *		fitPtTCMSysErr;
TF1 *		fitPtPowerlawSysErr;
TF1 *		fitPtModPowerlawSysErr;

TString 		forOutput;
TH1D* 		histoNumberOfEvents;

TH1D*		histoPi0ToChargedPhojet;
TH1D*		histoPi0ToChargedPythia;
TGraphAsymmErrors* graphSystErrRatio;

Double_t fPileUpCorrectionConv7TeV = 1-0.0105;

Int_t colorScheme;
Bool_t use7TeVPytPho=kFALSE;
Int_t fExampleBinPi0     = 7;
Int_t fExampleBinEta     = 6;

void SelectExampleBin(TString optionEnergy, TString useSameBinningPi0Eta){
	
    if(optionEnergy.CompareTo("7TeV") == 0){
        fExampleBinPi0     = 7;
        fExampleBinEta     = 6;
        if (!useSameBinningPi0Eta.CompareTo("")==0)fExampleBinPi0 = fExampleBinEta;
    } else if( optionEnergy.CompareTo("8TeV") == 0) {
        fExampleBinPi0     = 7;
        fExampleBinEta     = 6;
        if (!useSameBinningPi0Eta.CompareTo("")==0)fExampleBinPi0 = fExampleBinEta;
    } else if( optionEnergy.CompareTo("2.76TeV") == 0) {
        fExampleBinPi0     = 7;
        fExampleBinEta     = 4;
        if (!useSameBinningPi0Eta.CompareTo("")==0)fExampleBinPi0 = fExampleBinEta;
    } else if( optionEnergy.CompareTo("900GeV") == 0) {
        fExampleBinPi0     = 7;
        fExampleBinEta     = 2;
        if (!useSameBinningPi0Eta.CompareTo("")==0)fExampleBinPi0 = fExampleBinEta;
    } else if( optionEnergy.CompareTo("pPb_5.023TeV") == 0 ) {
        fExampleBinPi0     = 7;
        fExampleBinEta     = 4;
        if (!useSameBinningPi0Eta.CompareTo("")==0)fExampleBinPi0 = fExampleBinEta;
    } else if( optionEnergy.CompareTo("PbPb_2.76TeV") == 0 ) {
        fExampleBinPi0     = 7;
        fExampleBinEta     = 6;
    } else if( optionEnergy.CompareTo("13TeV") == 0 ) {
        fExampleBinPi0     = 2;
        fExampleBinEta     = 2;
        if (!useSameBinningPi0Eta.CompareTo("")==0)fExampleBinPi0 = fExampleBinEta;
    }
}
// ****** FUNCTION TO PLOT MASS AND FWHM TOGETHER ********
void PlotFWHMMass(TString plotName,TH1D* histoData,TH1D* histoMC,TH1D* histoData2,TH1D* histoMC2,Double_t yMin,Double_t yMax,TString yAxisLabel,TString meson,TString optionEnergy,Int_t mode,TString outputDir,TString prefix2,TString cutSelection,TString suffix,TString period = "",Bool_t plotBoth = kFALSE,TGraphAsymmErrors* graphSystErr=NULL,TGraphAsymmErrors* graphSystErr2=NULL){
	
	TString mesonString		="#pi^{0}";
	if(plotName.Contains("Eta"))mesonString		="#eta";
	
	Double_t maxpT      = 1.2*histoData->GetXaxis()->GetBinCenter(histoData->GetNbinsX())+(histoData->GetXaxis()->GetBinWidth(histoData->GetNbinsX()))/2;
	
	histoData->Scale(1000);
	histoData2->Scale(1000);
	histoMC->Scale(1000);
	histoMC2->Scale(1000);
	
	Double_t minYMass;
	Double_t maxYMass;
	if(plotName.Contains("Pi0")){
		minYMass=131.1;
		maxYMass=143.9;
	}else{
		minYMass=535.1;
		maxYMass=579.9;
	}
	Double_t maxYFWHM=2.8*histoData->GetBinContent(histoData->GetMaximumBin());
	
	Double_t arrayBoundariesX1_4[2];
	Double_t arrayBoundariesY1_4[3];
	Double_t relativeMarginsX[3];
	Double_t relativeMarginsY[3];
	Size_t textSizeLabelsPixel             = 50;
	ReturnCorrectValuesForCanvasScaling(1350,1250, 1, 2,0.09, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

	TCanvas* canvasMassWidthPi0     = new TCanvas("canvasMassWidthPi0","",0,0,1350,1250);  // gives the page size
	DrawGammaCanvasSettings( canvasMassWidthPi0,  0.13, 0.02, 0.03, 0.06);

	TPad* padWidthPi0               = new TPad("padWidthPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
	DrawGammaPadSettings( padWidthPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
	padWidthPi0->Draw();

	TPad* padMassPi0                = new TPad("padMassPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
	DrawGammaPadSettings( padMassPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
	padMassPi0->Draw();

	TPad* padMassLegend1            = new TPad("padMassLegend1", "", 0.13, 0.36, 0.52, 0.52,-1, -1, -2);
	DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
	padMassLegend1->SetFillStyle(0);
	padMassLegend1->Draw();

	padWidthPi0->cd();
	padWidthPi0->SetLogx(); 
	
	Double_t margin                 = relativeMarginsX[0]*2.7*1350;
	Double_t textsizeLabelsWidth    = 0;
	Double_t textsizeFacWidth       = 0;
	if (padWidthPi0->XtoPixel(padWidthPi0->GetX2()) < padWidthPi0->YtoPixel(padWidthPi0->GetY1())){
		textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi0->XtoPixel(padWidthPi0->GetX2()) ;
		textsizeFacWidth            = (Double_t)1./padWidthPi0->XtoPixel(padWidthPi0->GetX2()) ;
	} else {
		textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi0->YtoPixel(padWidthPi0->GetY1());
		textsizeFacWidth            = (Double_t)1./padWidthPi0->YtoPixel(padWidthPi0->GetY1());
	}

	TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, 0.23, maxpT ,1000., -30, 40);
	SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
							  0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
	histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-1.,maxYFWHM);
	histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
	histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
	histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
	histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
	histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
	histo2DAllPi0FWHM->DrawCopy(); 
	
	Color_t colorData   = kBlack;//, kGray+1, kRed+2, kBlue+2, kGreen+3, kCyan+2, kViolet, kMagenta+2,  kRed-2, kBlue-2};
    Color_t colorMC    	= kRed+2;//, kRed+2, kBlue+2, kGreen+3, kCyan+2, kViolet, kMagenta+2,  kRed-2, kBlue-2}; 
    if(colorScheme==1){
		colorData   = GetColorDefaultColor(optionEnergy,"","");
	}
	if(colorScheme==2){
		colorData   = GetDefaultColorDiffDetectors( detSystem,0,kFALSE);
	}
    Marker_t marker    	= 20;//, 20, 21, 34, 29, 33, 21, 27, 28, 30 };21
    Marker_t markerMC  	= 24;//, 24, 25, 28, 30, 27, 25, 27, 28, 30 };25
	Size_t markerSize 	= 2.5;
	
	DrawGammaSetMarker(histoData, marker, markerSize, colorData , colorData);
	histoData->Draw("p,same,z");
	DrawGammaSetMarker(histoMC, markerMC, markerSize, colorMC , colorMC);
	histoMC->Draw("p,same,z");


	TLatex *labelLegendAMass    = new TLatex(0.13,0.06,"a)");
	SetStyleTLatex( labelLegendAMass, textSizeLabelsPixel,4);
	labelLegendAMass->SetTextFont(43);
	labelLegendAMass->Draw();

	TLatex *labelMassPerf       = new TLatex(0.13,0.87,"ALICE work-in-progress");
	SetStyleTLatex( labelMassPerf, textSizeLabelsPixel,4);
	labelMassPerf->SetTextFont(43);
	labelMassPerf->Draw();        
	TLatex *labelMassEnergy     = new TLatex(0.13,0.78,collisionSystem.Data());
	SetStyleTLatex( labelMassEnergy, textSizeLabelsPixel,4);
	labelMassEnergy->SetTextFont(43);
	labelMassEnergy->Draw();
	TLatex *labelMassPi0        = new TLatex(0.13,0.69,Form("%s #rightarrow #gamma#gamma",mesonString.Data()));
	SetStyleTLatex( labelMassPi0, textSizeLabelsPixel,4);
	labelMassPi0->SetTextFont(43);
	labelMassPi0->Draw();  
	      
	padMassPi0->cd();
	padMassPi0->SetLogx();

	Double_t textsizeLabelsMass         = 0;
	Double_t textsizeFacMass            = 0;
	if (padMassPi0->XtoPixel(padMassPi0->GetX2()) <padMassPi0->YtoPixel(padMassPi0->GetY1()) ){
		textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
		textsizeFacMass                 = (Double_t)1./padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
	} else {
		textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi0->YtoPixel(padMassPi0->GetY1());
		textsizeFacMass                 = (Double_t)1./padMassPi0->YtoPixel(padMassPi0->GetY1());
	}

	TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, 0.23, maxpT, 1000., minYMass, maxYMass);
	SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass,textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
	histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
	histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
	histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
	histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
	histo2DAllPi0Mass->GetXaxis()->SetLabelOffset(-0.015);
	histo2DAllPi0Mass->DrawCopy(); 
	
	if(plotName.Contains("Pi0")){
		DrawGammaLines(0.23, 25. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);
	}else{
		DrawGammaLines(0.23, 25. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,0.1, kGray);
	}
	
	DrawGammaSetMarker(histoData2, marker, markerSize, colorData , colorData);
	histoData2->Draw("p,same,z");
	DrawGammaSetMarker(histoMC2, markerMC, markerSize, colorMC , colorMC);
	histoMC2->Draw("p,same,z");
	
	TLatex *labelLegendBMass            = new TLatex(0.13,0.22,"b)");
	SetStyleTLatex( labelLegendBMass, textSizeLabelsPixel,4);
	labelLegendBMass->SetTextFont(43);
	labelLegendBMass->Draw();
	//********************************** Defintion of the Legend **************************************************    
	Double_t columnsLegendMass2[3]      = {0.,0.57,0.84};
	Double_t rowsLegendMass2[4]         = {0.75,0.5,0.25,0.01};
	//******************* Offsets ***********************
	Double_t offsetMarkerXMass2         = 0.1;
	Double_t offsetMarkerYMass2         = 0.1;
	//****************** Scale factors ******************
	Double_t scaleMarkerMass2           = 1.2;
	padMassLegend1->cd();
	//****************** first Column **************************************************
	TLatex *textMassPCM                 = new TLatex(columnsLegendMass2[0],rowsLegendMass2[1],detSystem);
	SetStyleTLatex( textMassPCM, textSizeLabelsPixel,4);
	textMassPCM->SetTextFont(43);
	textMassPCM->Draw();
	//****************** second Column *************************************************
	TLatex *textMassData                = new TLatex(columnsLegendMass2[1],rowsLegendMass2[0] ,"Data");
	SetStyleTLatex( textMassData, textSizeLabelsPixel,4);
	textMassData->SetTextFont(43);
	textMassData->Draw();
	TLatex *textMassMC                  = new TLatex(columnsLegendMass2[2] ,rowsLegendMass2[0],"MC");
	SetStyleTLatex( textMassMC, textSizeLabelsPixel,4);
	textMassMC->SetTextFont(43);
	textMassMC->Draw();

	TMarker* markerPCMPi0Mass        = CreateMarkerFromHisto(histoData,columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
	markerPCMPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2);

	TMarker* markerPCMPi0MassMC      = CreateMarkerFromHisto(histoMC,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
	markerPCMPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[1]+ offsetMarkerYMass2);

    canvasMassWidthPi0->Update();
    canvasMassWidthPi0->SaveAs(Form("%s/%s_%s.%s",outputDir.Data(),prefix2.Data(),plotName.Data(),suffix.Data()));
}

// ****** FUNCTION TO PLOT THE ETA/PI0 RATIO *************
void calculateEtaPi0Ratio(TH1D* histoCorrectedYieldPi0,TH1D* histoCorrectedYieldEta, Double_t maxPt,TString optionEnergy,TGraphAsymmErrors* graphSystErrRatio){
	Size_t textSizeEtaToPi0 = 0.04;            
	DrawGammaLines(0., maxPt , 0.45, 0.45, 0.3, kGray+2);
	DrawGammaLines(0., maxPt , 0.4, 0.4, 0.3, kGray, 7);
	DrawGammaLines(0., maxPt , 0.5, 0.5, 0.3, kGray, 7);
	// create legend: case1->with pythia and phojet spectra, case2->without
	TLegend* legendEtaToPi0; 
	if(plotPythiaPhojetInRatio){
		legendEtaToPi0 	= GetAndSetLegend2(0.14, 0.83-4*0.85*0.04, 0.47,0.83,28);
	}else{
		legendEtaToPi0	= GetAndSetLegend2(0.14, 0.83-2*0.85*0.04, 0.47,0.83,28);
	}
	histoRatioEtaPi0 = (TH1D*) histoCorrectedYieldEta->Clone();
    histoRatioEtaPi0->Divide(histoRatioEtaPi0,histoCorrectedYieldPi0,1.,1.,"");
	// plot pythia and phojet spectra
	if(plotPythiaPhojetInRatio){
		DrawGammaSetMarker(histoEtaToPi0Phojet, 24, 1, kRed+2, kRed+2);
		histoEtaToPi0Phojet->Draw("same");
		DrawGammaSetMarker(histoEtaToPi0Pythia, 24, 1, kBlue+2, kBlue+2);
		histoEtaToPi0Pythia->Draw("same");
    }
    // set color scheme for data
    Color_t colorDataEPR = kBlack;
    Color_t colorSysEPR = kBlack;
    if(colorScheme==1){
		colorDataEPR   = GetColorDefaultColor(optionEnergy,"","");
		colorSysEPR    = GetColorDefaultColor(optionEnergy,"","",kTRUE);
	}
    if(colorScheme==2){
		colorDataEPR   = GetDefaultColorDiffDetectors( detSystem,0,kFALSE);
		colorSysEPR    = GetDefaultColorDiffDetectors( detSystem,0,kTRUE);
	}
	// draw systematic errors and data
    DrawGammaSetMarkerTGraphAsym(graphSystErrRatio,26,0,colorSysEPR,colorSysEPR,2,kTRUE);
	graphSystErrRatio->Draw("p,2,same");
	DrawGammaSetMarker(histoRatioEtaPi0, 24, 1.5, colorDataEPR, colorDataEPR);
	histoRatioEtaPi0->DrawCopy("same,x0,e1"); 
	
	// fill legend with entries
	legendEtaToPi0->AddEntry(histoRatioEtaPi0,"data","p");
	legendEtaToPi0->AddEntry(graphSystErrRatio,"systematic uncertainty","f");
	if(plotPythiaPhojetInRatio){
		if(use7TeVPytPho){
			legendEtaToPi0->AddEntry(histoEtaToPi0Phojet,"Phojet 7TeV","p");
			legendEtaToPi0->AddEntry(histoEtaToPi0Pythia,"Pythia 7TeV","p");
		}else{
			legendEtaToPi0->AddEntry(histoEtaToPi0Phojet,"Phojet","p");
			legendEtaToPi0->AddEntry(histoEtaToPi0Pythia,"Pythia","p");
		}
    }
	legendEtaToPi0->Draw();
}
// ****** FUNCTION TO FIT THE FINAL SPECTRA **************
void makeFinalFits(TH1D* histoCorrYield,Double_t maxPt,TString useSameBinningPi0Eta,TString strMeson,Bool_t plotRatios=kFALSE,TGraphAsymmErrors* graphSystErrRatio=NULL){
	TH1D* histoFitting = (TH1D*) histoCorrYield->Clone();
	
	// Fit histCorr with Levy
    fitPtLevy = FitObject("l","fitPtLevy",strMeson,histoFitting,minPtForFits,maxPt,NULL,"QNRME+");
    DrawGammaSetMarkerTF1(fitPtLevy, 1, 1.5, kBlue);
    
    // Fit histCorr with TCM
    Double_t paramTCM[5] = {histoCorrYield->GetBinContent(0),0.3,histoCorrYield->GetBinContent(0)/1000,0.8,3};
    fitPtTCM = FitObject("tcm","fitPtTCM",strMeson,histoFitting,minPtForFits,maxPt,paramTCM,"QNRMEX0+");
    DrawGammaSetMarkerTF1( fitPtTCM, 1, 1.5, kGreen+2);
    fitPtTCM->SetLineStyle(7);
    
    // Fit histCorr with Powerlaw
    fitPtPowerlaw = FitObject("p","fitPtPowerlaw",strMeson,histoFitting,1.5,maxPt,NULL,"QNRME+");
    DrawGammaSetMarkerTF1( fitPtPowerlaw, 1, 1.5, kTeal);
    
    // Fit histCorr with ModPowerlaw
    fitPtModPowerlaw = FitObject("m","fitPtModPowerlaw",strMeson,histoFitting,0.3,maxPt,NULL,"QNRME+");
    DrawGammaSetMarkerTF1( fitPtModPowerlaw, 1, 1.5, kMagenta+2);

    if(plotRatios){
		// Prepare histograms for ratios 
		histoRatioFitLevy 	= (TH1D*) histoCorrYield->Clone();
		histoRatioFitTCM 	= (TH1D*) histoCorrYield->Clone();
		histoRatioFitPow 	= (TH1D*) histoCorrYield->Clone();
		histoRatioFitModPow = (TH1D*) histoCorrYield->Clone();
		// Calculating ratios 
		histoRatioFitLevy 	= CalculateHistoRatioToFit (histoRatioFitLevy, fitPtLevy); 
		histoRatioFitTCM 	= CalculateHistoRatioToFit (histoRatioFitTCM, fitPtTCM);
		histoRatioFitPow 	= CalculateHistoRatioToFit (histoRatioFitPow, fitPtPowerlaw);
		histoRatioFitModPow = CalculateHistoRatioToFit (histoRatioFitModPow, fitPtModPowerlaw);
		// Prepare histograms for ratios 
		graphRatioFitLevySysErr 	= (TGraphAsymmErrors*) graphSystErrRatio->Clone();
		graphRatioFitTCMSysErr 		= (TGraphAsymmErrors*) graphSystErrRatio->Clone();
		graphRatioFitPowSysErr 		= (TGraphAsymmErrors*) graphSystErrRatio->Clone();
		graphRatioFitModPowSysErr 	= (TGraphAsymmErrors*) graphSystErrRatio->Clone();
		cout << __LINE__ << endl;
		// Calculating ratios 
		graphRatioFitLevySysErr 	= CalculateGraphErrRatioToFit (graphRatioFitLevySysErr, fitPtLevy); 
		graphRatioFitTCMSysErr 		= CalculateGraphErrRatioToFit (graphRatioFitTCMSysErr, fitPtTCM);
		graphRatioFitPowSysErr 		= CalculateGraphErrRatioToFit (graphRatioFitPowSysErr, fitPtPowerlaw);
		graphRatioFitModPowSysErr 	= CalculateGraphErrRatioToFit (graphRatioFitModPowSysErr, fitPtModPowerlaw);
		cout << __LINE__ << endl;
		// Draw ratios
		DrawGammaSetMarkerTGraphAsym(graphRatioFitLevySysErr,26,0,kBlue-8,kBlue-8,2,kTRUE);
		graphRatioFitLevySysErr->Draw("p,2,same");
		DrawGammaSetMarkerTGraphAsym(graphRatioFitTCMSysErr,26,0,kGreen-8,kGreen-8,2,kTRUE);
		graphRatioFitTCMSysErr->Draw("p,2,same");
		DrawGammaSetMarkerTGraphAsym(graphRatioFitPowSysErr,26,0,kTeal-8,kTeal-8,2,kTRUE);
		graphRatioFitPowSysErr->Draw("p,2,same");
		DrawGammaSetMarkerTGraphAsym(graphRatioFitModPowSysErr,26,0,kMagenta-8,kMagenta-8,2,kTRUE);
		graphRatioFitModPowSysErr->Draw("p,2,same");
		cout << __LINE__ << endl;
		
		DrawGammaSetMarker(histoRatioFitLevy, 20, 1.5, kBlue, kBlue);
		histoRatioFitLevy->Draw("e1,x0,same");
		DrawGammaSetMarker(histoRatioFitTCM, 24, 1.5, kGreen+2, kGreen+2);
		histoRatioFitTCM->Draw("e1,x0,same");
		DrawGammaSetMarker(histoRatioFitPow, 21, 1.5, kTeal, kTeal);
		histoRatioFitPow->Draw("e1,x0,same");
		DrawGammaSetMarker(histoRatioFitModPow, 25, 1.5, kMagenta+2, kMagenta+2);
		histoRatioFitModPow->Draw("e1,x0,same");
		// Draw +-10% lines around unity
		DrawGammaLines(0., maxPt , 1.0, 1.0,0.1, kGray+2 );
        DrawGammaLines(0., maxPt , 1.1, 1.1,0.1, kGray, 7);
        DrawGammaLines(0., maxPt , 0.9, 0.9,0.1, kGray, 7);
        // Add legend with fit entries
        TLegend* legendFit = GetAndSetLegend2(0.22, 0.825, 0.55,0.825+4*0.84*0.04,28);
		legendFit->AddEntry(histoRatioFitLevy,"Levy fit");
		legendFit->AddEntry(histoRatioFitTCM,"TCM fit");
		legendFit->AddEntry(histoRatioFitPow,"Powerlaw fit");
		legendFit->AddEntry(histoRatioFitModPow,"ModPowerlaw fit");
        legendFit->Draw();
	}else{
		// Draw fits
		fitPtLevy->Draw("same");
		fitPtTCM->Draw("same");
		fitPtPowerlaw->Draw("same");
		fitPtModPowerlaw->Draw("same");

		// Add legend with fit entries
		TLegend* legendFit = GetAndSetLegend2(0.22, 0.125, 0.55,0.125+3*0.95*0.04,28);
		legendFit->AddEntry(fitPtLevy,"Levy fit");
		legendFit->AddEntry(fitPtTCM,"TCM fit");
		legendFit->AddEntry(fitPtPowerlaw,"Powerlaw fit");
		legendFit->AddEntry(fitPtModPowerlaw,"ModPowerlaw fit");
		legendFit->Draw();
	}
}

// ****** MAIN PLOTTING FUNCTION *************************
void PlotFinalOutput(TString plotName,TH1D* histoData,TH1D* histoMC,TH1D* histoData2,TH1D* histoMC2,Double_t yMin,Double_t yMax,TString yAxisLabel,TString meson,TString optionEnergy,Int_t mode,TString outputDir,TString prefix2,TString useSameBinningPi0Eta,TString cutSelection,TString suffix,TString period = "",Bool_t plotBoth = kFALSE,TGraphAsymmErrors* graphSystErr=NULL,TGraphAsymmErrors* graphSystErr2=NULL){
	// Define strings for legends
	TString strEtaBoth		="";
	TString strPi0Both		="";
	TString dataStr			="data";
	TString dataStr2		="data";
	TString MCStr			="MC";
	TString MCStr2			="MC";
	if(plotName.Contains("Comb")&&!plotName.Contains("FWHM")){dataStr="#pi^{0}";dataStr2="#eta";}
	TString mesonString		="Pi0";
	if(plotName.Contains("Eta"))mesonString		="Eta";
	if(plotBoth)strEtaBoth	="#eta";
	if(plotBoth)strPi0Both	="#pi^{0}";
	Size_t textSizeSpectra	= 0.04;
	TString collisionSystem = ReturnFullCollisionsSystem(optionEnergy);
    TString detectionProcess= ReturnFullTextReconstructionProcess(mode);
    
    // Create Canvas with three possible sizes depending on plot
    Int_t canvasSize		=900;
    if(plotName.Contains("Raw")||plotName.Contains("Corr")&&!plotName.Contains("Ratio"))canvasSize=1350;  
    if(plotName.Contains("EtaToPi0"))canvasSize=915;  
	TCanvas* canvasPlotFinal     = new TCanvas("canvasPlotFinal","",0,0,1000,canvasSize);
	
	// Set margins and title offsets for the canvas
	Double_t yTitleOffset = 1.1;
	Double_t xTitleOffset = 0.85;
	if(plotName.Contains("FWHM")){
		DrawGammaCanvasSettings( canvasPlotFinal, 0.09, 0.015, 0.035, 0.08);
		yTitleOffset = 1.;
	}else if(plotName.Contains("Mass")){
		DrawGammaCanvasSettings( canvasPlotFinal, 0.11, 0.015, 0.015, 0.08);
		yTitleOffset = 1.4;
	}else if(plotName.Contains("EtaToPi0")){
		DrawGammaCanvasSettings( canvasPlotFinal, 0.09, 0.015, 0.02, 0.1);
		yTitleOffset = 1.;
	}else if(plotName.Contains("Eff")){
		DrawGammaCanvasSettings( canvasPlotFinal, 0.09, 0.015, 0.015, 0.08);
		canvasPlotFinal->SetLogy(1);
	}else if(plotName.Contains("Acc")){
		DrawGammaCanvasSettings( canvasPlotFinal, 0.1, 0.015, 0.015, 0.08);
		yTitleOffset = 1.25;
	}else if(plotName.Contains("Raw")||plotName.Contains("Corr")){
		if(plotName.Contains("Ratio")){
			DrawGammaCanvasSettings( canvasPlotFinal, 0.09, 0.015, 0.015, 0.08);
			canvasPlotFinal->SetLogy(0);
			yTitleOffset = 1.2;
		}else{
			DrawGammaCanvasSettings( canvasPlotFinal, 0.16, 0.02, 0.015, 0.07);
			canvasPlotFinal->SetLogy(1);
			yTitleOffset = 1.7;
			xTitleOffset = 0.77;
		}
	}else{
		DrawGammaCanvasSettings( canvasPlotFinal, 0.09, 0.015, 0.035, 0.08);
		yTitleOffset = 1.;
	}
	// Get the maximum pT value from the histogram (last bin+0.5*error)
	Double_t maxpT      = histoData->GetXaxis()->GetBinCenter(histoData->GetNbinsX())+(histoData->GetXaxis()->GetBinWidth(histoData->GetNbinsX()))/2;
	// Create dummy histogram for plotting
    TH2F * histoDummy  = new TH2F("histoDummy","histoDummy",1000,0., maxpT,10000,yMin, yMax);
    SetStyleHistoTH2ForGraphs(histoDummy, "#it{p}_{T} (GeV/#it{c})",yAxisLabel,0.85*textSizeSpectra,textSizeSpectra, 0.85*textSizeSpectra,textSizeSpectra, xTitleOffset,yTitleOffset);
    // Set legend and labels for the plots
    TLegend* 	legendPlot;TLatex* 	labelEnergyPlot;TLatex* labelPi0Plot;TLatex* labelDetProcPlot;
    if(plotName.Contains("Acc")||plotName.Contains("Eff")){
		legendPlot 		= GetAndSetLegend2(0.62, 0.125, 0.95,0.125+0.85*textSizeSpectra,28);
		labelEnergyPlot = new TLatex(0.62, 0.13+3*0.85*textSizeSpectra,collisionSystem.Data());
		labelPi0Plot    = new TLatex(0.62, 0.13+2*0.85*textSizeSpectra,Form("%s #rightarrow #gamma#gamma",meson.Data()));
		labelDetProcPlot= new TLatex(0.62, 0.13+0.85*textSizeSpectra,detectionProcess.Data());
	}else if(plotName.Contains("Raw")||plotName.Contains("Corr")){
		if(plotName.Contains("SysErr")&&!plotName.Contains("Ratio")){
			legendPlot 		= GetAndSetLegend2(0.62, 0.825-0.85*textSizeSpectra, 0.95,0.825+0.85*textSizeSpectra,28);
		}else{
			legendPlot 		= GetAndSetLegend2(0.62, 0.825, 0.95,0.825+0.85*textSizeSpectra,28);
		}
		labelEnergyPlot = new TLatex(0.62, 0.83+3*0.85*textSizeSpectra,collisionSystem.Data());
		labelPi0Plot    = new TLatex(0.62, 0.83+2*0.85*textSizeSpectra,Form("%s #rightarrow #gamma#gamma",meson.Data()));
		labelDetProcPlot= new TLatex(0.62, 0.83+0.85*textSizeSpectra,detectionProcess.Data());
	}else{
		legendPlot	= GetAndSetLegend2(0.52, 0.87, 0.95, 0.87+(1.05*4/2*0.85*textSizeSpectra),28);
		labelEnergyPlot = new TLatex(0.14, 0.84+(1.02*2*textSizeSpectra*0.85),collisionSystem.Data());
		labelPi0Plot    = new TLatex(0.14, 0.84+0.99*textSizeSpectra*0.85,Form("%s #rightarrow #gamma#gamma",meson.Data()));
		labelDetProcPlot= new TLatex(0.14, 0.84,detectionProcess.Data());
	}

    canvasPlotFinal->cd();
    histoDummy->DrawCopy(); 
    legendPlot->SetNColumns(2);
    
    // set colors defaul and depending on the colorscheme set in ProduceFinalResultsV2.C
    Color_t colorData= kBlack;
    Color_t colorData2= kAzure+2;
    Color_t colorMC= kRed+2;
    Color_t colorMC2= kRed+2;
    Color_t colorSys= kGray+1;
    Color_t colorSys2= kAzure-8;
    
    if(colorScheme==1&&!plotName.Contains("Comb")){
		colorData   = GetColorDefaultColor(optionEnergy,"","");
		colorData2  = kAzure+2;
		colorMC    	= kRed+2;
		colorMC2   	= kRed+2;
		colorSys	= GetColorDefaultColor(optionEnergy,"","",kTRUE);
		colorSys2	= kAzure-8;
	}else if(colorScheme==2){
		colorData   = GetDefaultColorDiffDetectors( detSystem,0,kFALSE);
		colorData2  = kAzure+2;
		colorMC    	= kRed+2;
		colorMC2   	= kRed+2;
		colorSys	= GetDefaultColorDiffDetectors( detSystem,0,kTRUE);
		colorSys2	= kAzure-8;
	}
    Marker_t marker    	= 20;
    Marker_t marker2    = 21;
    Marker_t markerMC  	= 24;
    Marker_t markerMC2  = 25;

    Size_t pointSize    = 1.5;
    
    if(plotName.Contains("SysErr")){
		if(graphSystErr){
			DrawGammaSetMarkerTGraphAsym(graphSystErr,26,0,colorSys,colorSys,2,kTRUE);
			graphSystErr->Draw("p,2,same");
			legendPlot->SetNColumns(1);
		}
		if(graphSystErr2){
			DrawGammaSetMarkerTGraphAsym(graphSystErr2,26,0,colorSys2,colorSys2,2,kTRUE);
			graphSystErr2->Draw("p,2,same");
		}
	}
	DrawGammaSetMarker(histoData, marker, pointSize, colorData, colorData);
	if(!plotName.Contains("Ratio")){
		if(histoData){
			if(histoMC && plotName.Contains("AccEff")){histoData->Multiply(histoMC);histoMC=NULL;}
			histoData->DrawCopy("e1,same"); 
			if(!plotName.Contains("Acc")&&!plotName.Contains("Eff")||plotName.Contains("Comb"))legendPlot->AddEntry(histoData, Form("%s %s %s",period.Data(),strPi0Both.Data(),dataStr.Data()), "p");
		}
		if(histoData2){
			if(histoMC2 && plotName.Contains("AccEff")){histoData2->Multiply(histoMC2);histoMC2=NULL;}
			DrawGammaSetMarker(histoData2, marker2, pointSize, colorData2, colorData2);
			histoData2->DrawCopy("e1,same"); 
			legendPlot->AddEntry(histoData2, Form("%s %s %s",period.Data(),strEtaBoth.Data(),dataStr2.Data()), "p"); 
		}
		if(histoMC){
			DrawGammaSetMarker(histoMC, markerMC, pointSize, colorMC, colorMC);
			histoMC->DrawCopy("e1,same"); 
			legendPlot->AddEntry(histoMC, Form("%s %s %s", period.Data(),strPi0Both.Data(),MCStr.Data()), "p"); 
		}
		if(histoMC2){
			DrawGammaSetMarker(histoMC2, markerMC2, pointSize, colorMC2, colorMC2);
			histoMC2->DrawCopy("e1,same"); 
			legendPlot->AddEntry(histoMC2, Form("%s %s %s", period.Data(),strEtaBoth.Data(),MCStr2.Data()), "p"); 
		}
		
	}
	if(graphSystErr&&plotName.Contains("SysErr"))legendPlot->AddEntry(graphSystErr, "systematic uncertainty", "f"); 
	if(plotName.Contains("Fitted")){
		if(plotName.Contains("Ratio")){
			makeFinalFits(histoData,maxpT,useSameBinningPi0Eta,mesonString,kTRUE,graphSystErr);
		}else{
			makeFinalFits(histoData,maxpT,useSameBinningPi0Eta,mesonString);
		}
	}
	if(plotName.Contains("EtaToPi0"))calculateEtaPi0Ratio(histoData,histoData2,maxpT,optionEnergy,graphSystErr);
	
    legendPlot->Draw();
    SetStyleTLatex( labelEnergyPlot, 0.85*textSizeSpectra,4);
    labelEnergyPlot->Draw();

    SetStyleTLatex( labelPi0Plot, 0.85*textSizeSpectra,4);
    labelPi0Plot->Draw();
    
    SetStyleTLatex( labelDetProcPlot, 0.85*textSizeSpectra,4);
    labelDetProcPlot->Draw();

    canvasPlotFinal->Update();
    canvasPlotFinal->SaveAs(Form("%s/%s_%s.%s",outputDir.Data(),prefix2.Data(),plotName.Data(),suffix.Data()));
}
