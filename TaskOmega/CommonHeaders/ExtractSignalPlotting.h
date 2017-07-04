// provided by Gamma Conversion Group, $ALICE_ROOT/PWG4/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

// //______________________________ Declaration of the functions in this header file_____________________________________________
// void PlotInvMassSinglePtBin(TH1D* , TH1D* , TString , TString , Double_t* , TString);
// void PlotExampleInvMassBins(TH1D*, TH1D*, TH1D* , TF1* , Int_t , TString , TString , Double_t* , Float_t* , Double_t , TString , TString , TString , TString , TString , Double_t*, TString);
// void PlotInvMassInPtBins(TH1D** , TH1D** , TString , TString , TString , Double_t* , TString , TString , Int_t , Int_t , Int_t , Int_t , Double_t* , TString , Bool_t, TString , TString);
// void PlotInvMassSecondaryInPtBins(TH1D** , TH1D** , TH1D** ,TString , TString , TString , Double_t* , TString , TString ,  Int_t , Int_t , Int_t , Int_t , Double_t* , TString , Bool_t, TString , TString);
// void PlotInvMassRatioInPtBins(TH1D** , TString , TString , TString ,  Double_t* , TString , TString ,  Int_t , Int_t , Int_t , Int_t , Double_t* , TString , Bool_t , TString, TString);
// void PlotWithFitSubtractedInvMassSinglePtBin(TH1D * , TH1D** , TF1 * , TString , TString ,  Double_t* , Bool_t, TString);
// void PlotWithFitSubtractedInvMassSinglePtBin2(TH1D * , TF1 * , TString , TString ,  Double_t* , Bool_t, TString );
// void PlotWithFitSubtractedInvMassInPtBins(TH1D ** , TH1D** , TF1 ** , TString , TString , TString , Double_t* , TString , TString ,  Int_t , Int_t , Int_t , Int_t , Double_t* , TString , Bool_t, TString, TString);
// void PlotWithFitPeakPosInvMassInPtBins(TH1D ** ,TH1D **, TF1 ** , TString , TString , TString , Double_t* , TString , TString ,  Int_t , Int_t , Int_t , Int_t , Double_t* , TString , Bool_t, TString , TString);
// 

//MY ADDITION!!!
//________________________________ Plotting Single 1D Histogram ____________________________________________
void PlotSingle1DHistogram(TH1D* histo, TString namePlot, TString nameCanvas, TString XaxisName, TString YaxisName, TString suffix, Int_t XaxisScale, Int_t YaxisScale, Double_t* RangeX, Double_t* RangeY,TString addText, Double_t xPosLine1 = 0., Double_t xPosLine2 = 0.){
    TCanvas * canvas = new TCanvas(nameCanvas.Data(),"",1400,800);  // gives the page size
	canvas->SetTopMargin(0.08);
	canvas->SetBottomMargin(0.13);
	canvas->SetRightMargin(0.04);
	canvas->SetLeftMargin(0.12);

	
	canvas->SetLogx(XaxisScale);
	canvas->SetLogy(YaxisScale);
	
	histo->GetXaxis()->SetRangeUser(RangeX[0], RangeX[1]);
	histo->GetYaxis()->SetRangeUser(RangeY[0], RangeY[1]);
	
	DrawSingleHisto(histo, nameCanvas, XaxisName, YaxisName, kFALSE);
	TLine *line1,*line2;
	line1 = new TLine(xPosLine1,RangeY[0],xPosLine1,RangeY[1]);
	line2 = new TLine(xPosLine2,RangeY[0],xPosLine2,RangeY[1]);
	
	if (xPosLine1!=0. && xPosLine2!=0.)
	{
		line1->SetLineWidth(2);
		line1->SetLineColor(2);
		line1->Draw();
		
		line2->SetLineWidth(2);
		line2->SetLineColor(2);
		line2->Draw();
	}
	
    DrawAliceLogoOmega(0.7,0.85,0.1,0.05,1400,800,histo->GetEntries(),addText);

	canvas->Update();

	canvas->Print(namePlot.Data(),suffix.Data());
	delete canvas;
	delete line1;
	delete line2;
}

void PlotSingle2DHistogram(TH2D* histo, TString namePlot, TString nameCanvas, TString XaxisName, TString YaxisName, TString suffix, Int_t XaxisScale, Int_t YaxisScale, Int_t ZaxisScale, Double_t* RangeX, Double_t* RangeY,TString addText, Double_t xPosLine1 = 0., Double_t xPosLine2 = 0.){
    TCanvas * canvas = new TCanvas(nameCanvas.Data(),"",1600,800);  // gives the page size
	canvas->SetTopMargin(0.08);
	canvas->SetBottomMargin(0.13);
	canvas->SetRightMargin(0.1);
	canvas->SetLeftMargin(0.12);

	
	canvas->SetLogx(XaxisScale);
	canvas->SetLogy(YaxisScale);
	canvas->SetLogz(ZaxisScale);
	
	histo->GetXaxis()->SetRangeUser(RangeX[0], RangeX[1]);
	histo->GetYaxis()->SetRangeUser(RangeY[0], RangeY[1]);
	
	DrawSingleHisto(histo, nameCanvas, XaxisName, YaxisName, kTRUE);
	TLine *line1,*line2;
	line1 = new TLine(xPosLine1,RangeY[0],xPosLine1,RangeY[1]);
	line2 = new TLine(xPosLine2,RangeY[0],xPosLine2,RangeY[1]);
	
	if (xPosLine1!=0. && xPosLine2!=0.)
	{
		line1->SetLineWidth(2);
		line1->SetLineColor(2);
		line1->Draw();
		
		line2->SetLineWidth(2);
		line2->SetLineColor(2);
		line2->Draw();
	}
	
    DrawAliceLogoOmega(0.2,0.85,0.1,0.05,1400,800,histo->GetEntries(),addText);

	canvas->Update();

	canvas->Print(namePlot.Data(),suffix.Data());
	delete canvas;
}

//__________________________________________ Plotting all Invariant Mass bins _______________________________________________
void PlotMultipleSlicesOf2DHisto(TH2D* histo2D, TString namePlot, TString nameCanvas, TString namePad, Double_t* fPlottingRangeMeson, TString dateDummy, TString fMesonType, Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, Double_t* fRangeBinsPt, Double_t* RebinFactorsPtSlice, TString fDecayChannel, TString decayChannel,  TString fDetectionChannel, TString fEnergy){
	
	TH1D *histo[fNumberPtBins];
	TH1D *emptyhisto = new TH1D("empty", "empty", 1, fPlottingRangeMeson[0], fPlottingRangeMeson[1]);
	emptyhisto->SetBinContent(1,0.);
	TGaxis::SetMaxDigits(3);

    TCanvas * canvasDataSpectra = new TCanvas(nameCanvas.Data(),"",1400,800);  // gives the page size
	canvasDataSpectra->SetTopMargin(0.02);
	canvasDataSpectra->SetBottomMargin(0.02);
	canvasDataSpectra->SetRightMargin(0.02);
	canvasDataSpectra->SetLeftMargin(0.02);

	TPad * padDataSpectra = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
	padDataSpectra->SetFillColor(0);
	padDataSpectra->GetFrame()->SetFillColor(0);
	padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(0);
	padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
	padDataSpectra->Draw();

	cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
	cout << fMesonType.Data() << endl;
	Int_t place = 0;
	for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
		cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
		Double_t startPt = fRangeBinsPt[iPt];
		Double_t endPt = fRangeBinsPt[iPt+1];
		
		TString Name = Form("histo_%i",iPt);
		Double_t startBinIn2DHisto = histo2D->GetYaxis()->FindBin(startPt);
		Double_t endBinIn2DHisto = histo2D->GetYaxis()->FindBin(endPt);
		
		histo[iPt] = histo2D->ProjectionX(Name,startBinIn2DHisto,endBinIn2DHisto,"");
		histo[iPt]->Rebin(RebinFactorsPtSlice[iPt]);
		
		place = place + 1;						//give the right place in the page
		if (place == fColumnPlot){
			iPt--;
			padDataSpectra->cd(place);

			TString textAlice = "ALICE performance";
			TString textEvents;
			//if(fMonteCarloInfo){textEvents = "MC";} else
			textEvents = "Data";
			Double_t nPixels = 13;
			Double_t textHeight = 0.08;
			if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
				textHeight = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
			} else {
				textHeight = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
			}	
			
			Double_t startTextX = 0.15;
			Double_t startTextY = 0.8;
			Double_t differenceText = textHeight*1.25;

			TLatex *alice = 		new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
			TLatex *latexDate = 	new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
			TLatex *energy = 		new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
			TLatex *process = 		new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
			TLatex *detprocess = 	new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
			TLatex *events = 		new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight*1.3);
			alice->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
	
			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();

			process->SetNDC(); 
			process->SetTextColor(1);
			process->SetTextSize(textHeight);
			process->Draw();

			detprocess->SetNDC(); 
			detprocess->SetTextColor(1);
			detprocess->SetTextSize(textHeight);
			detprocess->Draw();

			/*events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();*/
		} else {

			padDataSpectra->cd(place);
			padDataSpectra->cd(place)->SetTopMargin(0.12);
			padDataSpectra->cd(place)->SetBottomMargin(0.15);
			padDataSpectra->cd(place)->SetRightMargin(0.02);
			padDataSpectra->cd(place)->SetLeftMargin(0.15);

			DrawGammaHisto( histo[iPt],
						Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
							Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
							fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
			DrawGammaHisto( emptyhisto,
				Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
				Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
				fPlottingRangeMeson[0],fPlottingRangeMeson[1],3);
		}
	}
	canvasDataSpectra->Print(namePlot.Data());
	delete padDataSpectra;
	delete canvasDataSpectra;
	delete emptyhisto;
	for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++) {delete histo[iPt];}
}

//END OF MY ADDITION!!!

//________________________________ Plotting Invariant mass in a single p_t bin ____________________________________________
void PlotInvMassSinglePtBin(TH1D* fHistoMappingGGInvMassPtBinPlot, TH1D* fHistoMappingBackNormInvMassPtBinPlot, TString namePlot, TString nameCanvas, Double_t* fPlottingRangeMeson, TString decayChannel){
    TCanvas * canvasDataSpectra = new TCanvas(nameCanvas.Data(),"",1400,800);  // gives the page size
	canvasDataSpectra->SetTopMargin(0.12);
	canvasDataSpectra->SetBottomMargin(0.15);
	canvasDataSpectra->SetRightMargin(0.05);
	canvasDataSpectra->SetLeftMargin(0.15);


	fHistoMappingGGInvMassPtBinPlot->GetXaxis()->SetRangeUser(10, 10000);
	DrawGammaHisto( fHistoMappingGGInvMassPtBinPlot,
				 "2.4 GeV/c < p_{T} < 2.6 GeV/c",
				 Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
				 fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);

	DrawGammaHisto( fHistoMappingBackNormInvMassPtBinPlot,
				"2.4 GeV/c < p_{T} < 2.6 GeV/c",
				Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
				fPlottingRangeMeson[0],fPlottingRangeMeson[1],1);

	canvasDataSpectra->Print(namePlot.Data());
	delete canvasDataSpectra;
}

//_______________________ Plotting Invariant mass with BG and subtraction in a single p_t bin __________________________________
void PlotExampleInvMassBins(TH1D* histoInvMassSignalWithBG, TH1D* histoInvMassSubstracted, TH1D* histoInvMassBG, TF1* fitSignal, Int_t exampleBin, TString outputDir, TString suffix, Double_t* fPlottingRangeMeson, Float_t* pictDrawingCoordinatesDummy, Double_t fNumberOfEvents, TString dateDummy, TString fMesonType, TString fSimulation, TString fPlottingType, TString fCollisionSystemDummy, Double_t* fRangeBinsPt, TString decayChannel){

	TCanvas* canvasLonelyBin = new TCanvas("canvasLonelyBin","",200,10,700,500);  // gives the page size
	DrawGammaCanvasSettings( canvasLonelyBin,  0.1, 0.03, 0.07, 0.1);

	Double_t startPt = fRangeBinsPt[exampleBin];
	Double_t endPt = fRangeBinsPt[exampleBin+1];

	Bool_t kDalitz = kFALSE;
	if(decayChannel.CompareTo("e^{+}e^{-}#gamma")==0){
		kDalitz=kTRUE;
	}
	
	canvasLonelyBin->cd();
	DrawGammaHistoBigger( histoInvMassSignalWithBG,
                 Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
				 Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
				 fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);

	DrawGammaHistoBigger( histoInvMassBG,
                  Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
				  Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
				  fPlottingRangeMeson[0],fPlottingRangeMeson[1],1);


    TLegend* legendExamplePlotInvMass = new TLegend(0.6,0.2,0.9,0.3);
	legendExamplePlotInvMass->SetTextSize(0.04);
	legendExamplePlotInvMass->SetFillColor(0);
	legendExamplePlotInvMass->AddEntry(histoInvMassSignalWithBG,"Data","ep");
	legendExamplePlotInvMass->AddEntry(histoInvMassBG,"combinatorial BG","l");
	legendExamplePlotInvMass->Draw();

	Bool_t lable;
	Bool_t mcFile;
	if(fMesonType.CompareTo("Pi0") == 0 || fMesonType.CompareTo("Pi0EtaBinning") == 0){ lable = kTRUE;
		} else { lable = kFALSE; }
	
	if(fSimulation.CompareTo("data")==0){ mcFile = kFALSE;
		} else {mcFile = kTRUE;}
	
	if (!(fMesonType.CompareTo("Omega") == 0))//ADDED BY ME
	if (!fPlottingType.CompareTo("thesis") == 0)DrawAliceLogoPi0PerformanceExtract(pictDrawingCoordinatesDummy[0], pictDrawingCoordinatesDummy[1], pictDrawingCoordinatesDummy[2], pictDrawingCoordinatesDummy[3], pictDrawingCoordinatesDummy[4], pictDrawingCoordinatesDummy[5], pictDrawingCoordinatesDummy[6], pictDrawingCoordinatesDummy[7], fNumberOfEvents ,fCollisionSystemDummy ,mcFile, kFALSE, lable ,700 ,500 ,dateDummy,kDalitz);

//THIS WAS ADDED BY ME!!!!
	if ((fMesonType.CompareTo("Omega") == 0) && (!fPlottingType.CompareTo("thesis") == 0))
		DrawAliceLogoOmegaPerformance(pictDrawingCoordinatesDummy[0], pictDrawingCoordinatesDummy[1], pictDrawingCoordinatesDummy[2], pictDrawingCoordinatesDummy[3], pictDrawingCoordinatesDummy[4], pictDrawingCoordinatesDummy[5], pictDrawingCoordinatesDummy[6], pictDrawingCoordinatesDummy[7], fNumberOfEvents ,fCollisionSystemDummy ,mcFile, kFALSE, lable ,700 ,500 ,dateDummy,kDalitz);
//END

	canvasLonelyBin->Update();
	canvasLonelyBin->Print(Form("%s/%s_%s_InvMassDistributionWithBG_PtBin_%i.%s",outputDir.Data(),fMesonType.Data(),fSimulation.Data(),exampleBin,suffix.Data()));

	TCanvas* canvasLonelyBin2 = new TCanvas("canvasLonelyBin2","",200,10,700,500);  // gives the page size
	DrawGammaCanvasSettings( canvasLonelyBin2, 0.1, 0.03, 0.07, 0.1);

	canvasLonelyBin2->cd();
	DrawGammaHistoBigger( histoInvMassSubstracted,
				 Form("       %3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
				 Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
				 fPlottingRangeMeson[0],fPlottingRangeMeson[1],0	);
	DrawGammaHistoBigger( histoInvMassSubstracted,
				Form("       %3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
				Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
				fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);

	if (fitSignal!=0x00){
		fitSignal->SetLineColor(kCyan+3);
		fitSignal->SetLineWidth(2);
		fitSignal->DrawCopy("same");
	}
	TLegend* legendExamplePlotInvMass2;
	if(fMesonType.CompareTo("Pi0") == 0 || fMesonType.CompareTo("Pi0EtaBinning") == 0){
		legendExamplePlotInvMass2 = new TLegend(0.6,0.4,0.95,0.5);
	} else {
        legendExamplePlotInvMass2 = new TLegend(0.6,0.82,0.95,0.92);
	}
    if(fMesonType.CompareTo("Omega") == 0 || fMesonType.CompareTo("Eta") == 0) legendExamplePlotInvMass2 = new TLegend(0.15,0.55,0.5,0.65);
	
	legendExamplePlotInvMass2->SetTextSize(0.033);	legendExamplePlotInvMass2->SetFillColor(0);
	legendExamplePlotInvMass2->AddEntry(histoInvMassSubstracted,"Data (BG subtracted)","ep");
	legendExamplePlotInvMass2->AddEntry(fitSignal,"Fit","l");
	legendExamplePlotInvMass2->Draw();

	if (!(fMesonType.CompareTo("Omega") == 0))//ADDED BY ME
	if (!fPlottingType.CompareTo("thesis") == 0)  DrawAliceLogoPi0PerformanceExtract(pictDrawingCoordinatesDummy[0], pictDrawingCoordinatesDummy[1], pictDrawingCoordinatesDummy[2], pictDrawingCoordinatesDummy[3], pictDrawingCoordinatesDummy[4], pictDrawingCoordinatesDummy[5], pictDrawingCoordinatesDummy[6], pictDrawingCoordinatesDummy[7], fNumberOfEvents ,fCollisionSystemDummy, mcFile, kFALSE, lable, 700, 500, dateDummy,kDalitz);

//THIS WAS ADDED BY ME!!!!
	if ((fMesonType.CompareTo("Omega") == 0) && (!fPlottingType.CompareTo("thesis") == 0))
		DrawAliceLogoOmegaPerformance(pictDrawingCoordinatesDummy[0], pictDrawingCoordinatesDummy[1], pictDrawingCoordinatesDummy[2], pictDrawingCoordinatesDummy[3], pictDrawingCoordinatesDummy[4], pictDrawingCoordinatesDummy[5], pictDrawingCoordinatesDummy[6], pictDrawingCoordinatesDummy[7], fNumberOfEvents ,fCollisionSystemDummy ,mcFile, kFALSE, lable ,700 ,500 ,dateDummy,kDalitz);
//END

	canvasLonelyBin2->Update();
	canvasLonelyBin2->Print(Form("%s/%s_%s_InvMassDistributionSubstracted_PtBin_%i.%s",outputDir.Data(),fMesonType.Data(),fSimulation.Data(),exampleBin,suffix.Data()));
}

//__________________________________________ Plotting all Invariant Mass bins _______________________________________________
void PlotInvMassInPtBins(TH1D** fHistoMappingGGInvMassPtBinPlot, TH1D** fHistoMappingBackNormInvMassPtBinPlot, TString namePlot, TString nameCanvas, TString namePad, Double_t* fPlottingRangeMeson, TString dateDummy, TString fMesonType, Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, Double_t* fRangeBinsPt, TString fDecayChannel, Bool_t fMonteCarloInfo, TString decayChannel,  TString fDetectionChannel, TString fEnergy, TString bckLegend = "mixed evt."){
	TGaxis::SetMaxDigits(3);

    TCanvas * canvasDataSpectra = new TCanvas(nameCanvas.Data(),"",1400,800);  // gives the page size
	canvasDataSpectra->SetTopMargin(0.0);
	canvasDataSpectra->SetBottomMargin(0.0);
	canvasDataSpectra->SetRightMargin(0.0);
    canvasDataSpectra->SetLeftMargin(0.0);

    TPad * padDataSpectra = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
	padDataSpectra->SetFillColor(0);
	padDataSpectra->GetFrame()->SetFillColor(0);
	padDataSpectra->SetBorderMode(0);
    padDataSpectra->SetLogy(0);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
	padDataSpectra->Draw();

	cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
	cout << fMesonType.Data() << endl;
	Int_t place = 0;
	for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
		cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
		Double_t startPt = fRangeBinsPt[iPt];
		Double_t endPt = fRangeBinsPt[iPt+1];

		place = place + 1;						//give the right place in the page
		if (place == fColumnPlot){
			iPt--;
            padDataSpectra->cd(place);


			TString textAlice = "ALICE performance";
			TString textEvents;
			if(fMonteCarloInfo){textEvents = "MC";}
			else {textEvents = "Data";}
			Double_t nPixels = 13;
			Double_t textHeight = 0.08;
			if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
				textHeight = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
			} else {
				textHeight = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
			}	
			
			Double_t startTextX = 0.10;
			Double_t startTextY = 0.9;
			Double_t differenceText = textHeight*1.25;

			TLatex *alice = 		new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
			TLatex *latexDate = 	new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
			TLatex *energy = 		new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
			TLatex *process = 		new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
			TLatex *detprocess = 	new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
			TLatex *events = 		new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight*1.3);
			alice->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
	
			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();

			process->SetNDC(); 
			process->SetTextColor(1);
			process->SetTextSize(textHeight);
			process->Draw();

			detprocess->SetNDC(); 
			detprocess->SetTextColor(1);
			detprocess->SetTextSize(textHeight);
			detprocess->Draw();

			events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();

			TLegend* legendData = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+2.)*differenceText);
			legendData->SetTextSize(textHeight);
			legendData->SetTextFont(62);
			legendData->SetFillColor(0);
			legendData->SetFillStyle(0);
			legendData->SetLineWidth(0);
			legendData->SetLineColor(0);
			legendData->SetMargin(0.15);
			Size_t markersize = fHistoMappingGGInvMassPtBinPlot[iPt-1]->GetMarkerSize();
			fHistoMappingGGInvMassPtBinPlot[iPt-1]->SetMarkerSize(3*markersize);
			legendData->AddEntry(fHistoMappingGGInvMassPtBinPlot[iPt-1],Form("same evt. M_{%s} (BG+Signal)",decayChannel.Data()),"ep");
			Size_t linesize = fHistoMappingBackNormInvMassPtBinPlot[iPt-1]->GetLineWidth();
			fHistoMappingBackNormInvMassPtBinPlot[iPt-1]->SetLineWidth(5*linesize);
			legendData->AddEntry(fHistoMappingBackNormInvMassPtBinPlot[iPt-1],Form("%s M_{%s}", bckLegend.Data(),decayChannel.Data()),"l");
			legendData->Draw();				

		} else {

            padDataSpectra->cd(place);
			padDataSpectra->cd(place)->SetTopMargin(0.12);
			padDataSpectra->cd(place)->SetBottomMargin(0.15);
			padDataSpectra->cd(place)->SetRightMargin(0.02);
            padDataSpectra->cd(place)->SetLeftMargin(0.15);
            //gPad->SetLogy();


            DrawGammaHisto( fHistoMappingGGInvMassPtBinPlot[iPt],
						Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
							Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
							fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);

			DrawGammaHisto( fHistoMappingBackNormInvMassPtBinPlot[iPt],
							Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
					Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
					fPlottingRangeMeson[0],fPlottingRangeMeson[1],1);
// 			if(fMonteCarloInfo){
// 				DrawGammaHistoColored( fHistoMappingTrueAllBckInvMassPtBins[iPt],
// 									Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
// 									Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
// 									fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,5);
// 				DrawGammaHistoColored( fHistoMappingTrueGGBckInvMassPtBins[iPt],
// 								Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
// 								Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
// 								fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,3);
// 				DrawGammaHistoColored( fHistoMappingTrueContBckInvMassPtBins[iPt],
// 									Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
// 									Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
// 									fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,6);
// 		}
		}
	}
	canvasDataSpectra->Print(namePlot.Data());
	delete padDataSpectra;
	delete canvasDataSpectra;
}

//________________________________________ Plot Invariant Mass Bin With Secondary Contribution _______________________________
void PlotInvMassSecondaryInPtBins(TH1D** fHistoMappingGGInvMassPtBinPlot, TH1D** fHistoMappingSecondaryTotalInvMassPtBinPlot, TH1D** fHistoMappingSecondaryK0sInvMassPtBinPlot,TString namePlot, TString nameCanvas, TString namePad, Double_t* fPlottingRangeMeson, TString dateDummy, TString fMesonType,  Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, Double_t* fRangeBinsPt, TString fDecayChannel, Bool_t fMonteCarloInfo, TString decayChannel,  TString fDetectionChannel, TString fEnergy){
	TGaxis::SetMaxDigits(3);

    TCanvas * canvasDataSpectra = new TCanvas(nameCanvas.Data(),"",1400,800);  // gives the page size
	canvasDataSpectra->SetTopMargin(0.0);
	canvasDataSpectra->SetBottomMargin(0.0);
	canvasDataSpectra->SetRightMargin(0.0);
	canvasDataSpectra->SetLeftMargin(0.0);

	TPad * padDataSpectra = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
	padDataSpectra->SetFillColor(0);
	padDataSpectra->GetFrame()->SetFillColor(0);
	padDataSpectra->SetBorderMode(0);
	padDataSpectra->SetLogy(0);
	padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
	padDataSpectra->Draw();

	cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
	cout << fMesonType.Data() << endl;
	Int_t place = 0;
	for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
		cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
		Double_t startPt = fRangeBinsPt[iPt];
		Double_t endPt = fRangeBinsPt[iPt+1];

		place = place + 1;						//give the right place in the page
		if (place == fColumnPlot){
			iPt--;
			padDataSpectra->cd(place);

			TString textAlice = "ALICE performance";
			TString textEvents;
			if(fMonteCarloInfo){textEvents = "MC";}
			else {textEvents = "Data";}
			Double_t nPixels = 13;
			Double_t textHeight = 0.08;
			if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
				textHeight = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
			} else {
				textHeight = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
			}	
			Double_t startTextX = 0.10;
			Double_t startTextY = 0.8;
			Double_t differenceText = textHeight*1.25;

			TLatex *alice = 		new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
			TLatex *latexDate = 	new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
			TLatex *energy = 		new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
			TLatex *process = 		new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
			TLatex *detprocess = 	new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
			TLatex *events = 		new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight*1.3);
			alice->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
	
			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();

			process->SetNDC(); 
			process->SetTextColor(1);
			process->SetTextSize(textHeight);
			process->Draw();

			detprocess->SetNDC(); 
			detprocess->SetTextColor(1);
			detprocess->SetTextSize(textHeight);
			detprocess->Draw();

			events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();
		} else {

			padDataSpectra->cd(place);
			padDataSpectra->cd(place)->SetTopMargin(0.12);
			padDataSpectra->cd(place)->SetBottomMargin(0.15);
			padDataSpectra->cd(place)->SetRightMargin(0.02);
            padDataSpectra->cd(place)->SetLeftMargin(0.15);

			DrawGammaHisto( fHistoMappingGGInvMassPtBinPlot[iPt],
						Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
							Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
							fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);

			DrawGammaHisto( fHistoMappingSecondaryTotalInvMassPtBinPlot[iPt],
							Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
					Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
					fPlottingRangeMeson[0],fPlottingRangeMeson[1],1);
			
			DrawGammaHistoColored( fHistoMappingSecondaryK0sInvMassPtBinPlot[iPt],
									Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
									Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
									fPlottingRangeMeson[0],fPlottingRangeMeson[1],0,2);
				
			
			
			}
	}
	canvasDataSpectra->Print(namePlot.Data());
	delete padDataSpectra;
	delete canvasDataSpectra;

}

//________________________________________ Plot Invariant Mass Bin With GG contamination for Dalitz _______________________________
void PlotInvMassBckGGInPtBins(TH1D** fHistoMappingGGInvMassPtBinPlot, TH1D** fHistoMappingSecondaryTotalInvMassPtBinPlot, TString namePlot, TString nameCanvas, TString namePad, Double_t* fPlottingRangeMeson, TString dateDummy, TString fMesonType,  Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, Double_t* fRangeBinsPt, TString fDecayChannel, Bool_t fMonteCarloInfo, TString decayChannel, TString fDetectionChannel, TString fEnergy){
	TGaxis::SetMaxDigits(3);

    TCanvas * canvasDataSpectra = new TCanvas(nameCanvas.Data(),"",1400,800);  // gives the page size
	canvasDataSpectra->SetTopMargin(0.00);
	canvasDataSpectra->SetBottomMargin(0.00);
	canvasDataSpectra->SetRightMargin(0.0);
	canvasDataSpectra->SetLeftMargin(0.00);

	TPad * padDataSpectra = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
	padDataSpectra->SetFillColor(0);
	padDataSpectra->GetFrame()->SetFillColor(0);
	padDataSpectra->SetBorderMode(0);
	padDataSpectra->SetLogy(0);
	padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
	padDataSpectra->Draw();

	cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
	cout << fMesonType.Data() << endl;
	
	Int_t place = 0;
	for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
		cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
		Double_t startPt = fRangeBinsPt[iPt];
		Double_t endPt = fRangeBinsPt[iPt+1];

		place = place + 1;						//give the right place in the page
		if (place == fColumnPlot){
			iPt--;
			padDataSpectra->cd(place);

			TString textAlice = "ALICE performance";
			TString textEvents;
			if(fMonteCarloInfo){textEvents = "MC";}
			else {textEvents = "Data";}
			Double_t nPixels = 13;
			Double_t textHeight = 0.08;
			if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
				textHeight = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
			} else {
				textHeight = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
			}	
			Double_t startTextX = 0.10;
			Double_t startTextY = 0.8;
			Double_t differenceText = textHeight*1.25;

			TLatex *alice = 		new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
			TLatex *latexDate = 	new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
			TLatex *energy = 		new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
			TLatex *process = 		new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
			TLatex *detprocess = 	new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
			TLatex *events = 		new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight*1.3);
			alice->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
	
			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();

			process->SetNDC(); 
			process->SetTextColor(1);
			process->SetTextSize(textHeight);
			process->Draw();

			detprocess->SetNDC(); 
			detprocess->SetTextColor(1);
			detprocess->SetTextSize(textHeight);
			detprocess->Draw();

			events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();

		} else {

			padDataSpectra->cd(place);
			padDataSpectra->cd(place)->SetTopMargin(0.12);
			padDataSpectra->cd(place)->SetBottomMargin(0.15);
			padDataSpectra->cd(place)->SetRightMargin(0.02);
			padDataSpectra->cd(place)->SetLeftMargin(0.15);

			DrawGammaHisto( fHistoMappingGGInvMassPtBinPlot[iPt],
						Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
							Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
							fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);

			DrawGammaHisto( fHistoMappingSecondaryTotalInvMassPtBinPlot[iPt],
							Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
					Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
					fPlottingRangeMeson[0],fPlottingRangeMeson[1],1);
			}
	}
	canvasDataSpectra->Print(namePlot.Data());
	delete padDataSpectra;
	delete canvasDataSpectra;

}


//__________________________________ Plotting Invariant Mass as Ratio _________________________________________________________
void PlotInvMassRatioInPtBins(TH1D** fHistoMappingGGInvMassPtBinPlot, TString namePlot, TString nameCanvas, TString namePad,  Double_t* fPlottingRangeMeson, TString dateDummy, TString fMesonType,  Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, Double_t* fRangeBinsPt, TString fDecayChannel, Bool_t fMonteCarloInfo, TString decayChannel,  TString fDetectionChannel, TString fEnergy){
    TCanvas * canvasDataSpectra = new TCanvas(nameCanvas.Data(),"",1400,800);  // gives the page size
	canvasDataSpectra->SetTopMargin(0.00);
	canvasDataSpectra->SetBottomMargin(0.00);
	canvasDataSpectra->SetRightMargin(0.00);
	canvasDataSpectra->SetLeftMargin(0.00);

	TPad * padDataSpectra = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
	padDataSpectra->SetFillColor(0);
	padDataSpectra->GetFrame()->SetFillColor(0);
	padDataSpectra->SetBorderMode(0);
   	padDataSpectra->SetLogy(1);
	padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
	padDataSpectra->Draw();

	cout<<"fColumnPlot: "<<fColumnPlot<<" fRowPlot: "<<fRowPlot<<endl;
	cout << fMesonType.Data() << endl;
	Int_t place = 0;

	for(Int_t iPt=fStartBinPtRange;iPt<fNumberPtBins;iPt++){
		cout<<"Pt: "<<iPt<<" of "<<fNumberPtBins<<endl;
		Double_t startPt = fRangeBinsPt[iPt];
		Double_t endPt = fRangeBinsPt[iPt+1];

		place = place + 1;						//give the right place in the page

		if (place == fColumnPlot){
			iPt--;
			padDataSpectra->cd(place);

			TString textAlice = "ALICE performance";
			TString textEvents;
			if(fMonteCarloInfo){textEvents = "MC";}
			else {textEvents = "Data";}
			Double_t nPixels = 13;
			Double_t textHeight = 0.08;
			if (padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) < padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1())){
				textHeight = (Double_t)nPixels/padDataSpectra->cd(place)->XtoPixel(padDataSpectra->cd(place)->GetX2()) ;
			} else {
				textHeight = (Double_t)nPixels/padDataSpectra->cd(place)->YtoPixel(padDataSpectra->cd(place)->GetY1());
			}	
			Double_t startTextX = 0.10;
			Double_t startTextY = 0.8;
			Double_t differenceText = textHeight*1.25;

			TLatex *alice = 		new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
			TLatex *latexDate = 	new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
			TLatex *energy = 		new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
			TLatex *process = 		new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
			TLatex *detprocess = 	new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
			TLatex *events = 		new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight*1.3);
			alice->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
	
			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();

			process->SetNDC(); 
			process->SetTextColor(1);
			process->SetTextSize(textHeight);
			process->Draw();

			detprocess->SetNDC(); 
			detprocess->SetTextColor(1);
			detprocess->SetTextSize(textHeight);
			detprocess->Draw();

			events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();

		} else {

			padDataSpectra->cd(place);
			padDataSpectra->cd(place)->SetTopMargin(0.12);
			padDataSpectra->cd(place)->SetBottomMargin(0.15);
			padDataSpectra->cd(place)->SetRightMargin(0.02);
			padDataSpectra->cd(place)->SetLeftMargin(0.15);

			DrawGammaHisto( fHistoMappingGGInvMassPtBinPlot[iPt],
						 Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
						 Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
						 fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);

			DrawGammaHisto( fHistoMappingGGInvMassPtBinPlot[iPt],
						Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
						Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
						fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);
		}
	}
	canvasDataSpectra->Print(namePlot.Data());
	delete padDataSpectra;
	delete canvasDataSpectra;
}

//____________________________ Plotting Invariant Mass with Subtraction for Single Bin ______________________________________
void PlotWithFitSubtractedInvMassSinglePtBin(TH1D * fHistoMappingSignalInvMassPtBinPlot, TH1D** fHistoMappingTrueMesonInvMassPtBinsPlot, TF1 * fFitSignalInvMassPtBinPlot, TString namePlot, TString nameCanvas,  Double_t* fPlottingRangeMeson, Bool_t fMonteCarloInfo, TString decayChannel){
    TCanvas *canvasDataFit = new TCanvas(nameCanvas.Data(),"",1400,800);  // gives the page size
	canvasDataFit->SetTopMargin(0.12);
	canvasDataFit->SetBottomMargin(0.15);
	canvasDataFit->SetRightMargin(0.05);
	canvasDataFit->SetLeftMargin(0.15);

	fHistoMappingSignalInvMassPtBinPlot->SetAxisRange(fPlottingRangeMeson[0],fPlottingRangeMeson[1]);
// 	cout<<"Maximum::"<<fHistoMappingSignalInvMassPtBinPlot->GetMaximum()<<endl;
	DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot,
					"2.4 GeV/c < p_{T} < 2.6 GeV/c",
					Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
					fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
	if(fMonteCarloInfo){
		fHistoMappingTrueMesonInvMassPtBinsPlot[14]->SetLineColor(kRed);
		fHistoMappingTrueMesonInvMassPtBinsPlot[14]->SetLineWidth(0.5);
		fHistoMappingTrueMesonInvMassPtBinsPlot[14]->DrawCopy("same");
	}
	if(fMonteCarloInfo) fHistoMappingTrueMesonInvMassPtBinsPlot[14]->DrawCopy("same");
	if (fFitSignalInvMassPtBinPlot!=0x00){
		fFitSignalInvMassPtBinPlot->SetLineColor(kCyan+3);
		fFitSignalInvMassPtBinPlot->SetLineWidth(0.7);
		fFitSignalInvMassPtBinPlot->DrawCopy("same");
	} 
	canvasDataFit->Print(namePlot.Data());
	delete canvasDataFit;
}

//____________________________ Plotting Invariant Mass with Subtraction for Single Bin ______________________________________
void PlotWithFitSubtractedInvMassSinglePtBin2(TH1D * fHistoMappingSignalInvMassPtBinPlot, TH1D * fHistoMappingSignalInvMassPtBinPlot2, TH1D * fHistoMappingSignalInvMassPtBinPlot3, TF1 * fFitSignalInvMassPtBinPlot, TF1 * fFitBGInvMassPtBinPlot,TString namePlot, TString nameCanvas,  Double_t* fPlottingRangeMeson, Bool_t fMonteCarloInfo, TString decayChannel) {
    TCanvas *canvasDataFit = new TCanvas(nameCanvas.Data(),"",1400,800);  // gives the page size
	canvasDataFit->SetTopMargin(0.12);
	canvasDataFit->SetBottomMargin(0.15);
	canvasDataFit->SetRightMargin(0.05);
	canvasDataFit->SetLeftMargin(0.15);
	fHistoMappingSignalInvMassPtBinPlot->SetAxisRange(fPlottingRangeMeson[0],fPlottingRangeMeson[1]);
// 	cout<<"Maximum::"<<fHistoMappingSignalInvMassPtBinPlot->GetMaximum()<<endl;
	DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot,
					"2.4 GeV/c < p_{T} < 2.6 GeV/c",
					Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
					fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);

	fHistoMappingSignalInvMassPtBinPlot2->SetMarkerColor(kRed+2);	
	fHistoMappingSignalInvMassPtBinPlot2->SetMarkerSize(0.5);	
	fHistoMappingSignalInvMassPtBinPlot2->SetLineColor(kRed+2);
	fHistoMappingSignalInvMassPtBinPlot2->DrawCopy("same,p,e1");

	fHistoMappingSignalInvMassPtBinPlot3->SetMarkerColor(kCyan+1);	
	fHistoMappingSignalInvMassPtBinPlot3->SetMarkerSize(0.5);	
	fHistoMappingSignalInvMassPtBinPlot3->SetLineColor(kCyan+1);
	fHistoMappingSignalInvMassPtBinPlot3->DrawCopy("same,p,e1");
	
	if (fFitBGInvMassPtBinPlot!=0x00){
		fFitBGInvMassPtBinPlot->SetLineColor(kCyan+3);
		fFitBGInvMassPtBinPlot->SetLineWidth(0.7);
		fFitBGInvMassPtBinPlot->DrawCopy("same");
	} 
	if (fFitSignalInvMassPtBinPlot!=0x00){
		fFitSignalInvMassPtBinPlot->SetLineColor(kBlack);
		fFitSignalInvMassPtBinPlot->SetLineWidth(0.7);
		fFitSignalInvMassPtBinPlot->DrawCopy("same");
	} 

	canvasDataFit->Print(namePlot.Data());
	delete canvasDataFit;

        if(fMonteCarloInfo)fMonteCarloInfo=kFALSE; 
}

//____________________________ Plotting Invariant Mass Subtracted for all bins ________________________________________________________
void PlotWithFitSubtractedInvMassInPtBins(TH1D ** fHistoMappingSignalInvMassPtBinPlot, TH1D** fHistoMappingTrueMesonInvMassPtBinsPlot, TF1 ** fFitSignalInvMassPtBinPlot, TString namePlot, TString nameCanvas, TString namePad, Double_t* fPlottingRangeMeson, TString dateDummy, TString fMesonType,  Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, Double_t* fRangeBinsPt, TString fDecayChannel, Bool_t fMonteCarloInfo, TString decayChannel, TString fDetectionChannel, TString fEnergy, TString fTextMCvalidated ="", Bool_t labelData = kTRUE, TString fTextFit = "Fit", TString fTextMGammaGamma ="mixed evt. subtr. M_{#gamma#gamma}"){
    TCanvas *canvasDataFit = new TCanvas(nameCanvas.Data(),"",1400,800);  // gives the page size
	canvasDataFit->SetTopMargin(0.00);
	canvasDataFit->SetBottomMargin(0.00);
	canvasDataFit->SetRightMargin(0.00);
	canvasDataFit->SetLeftMargin(0.00);

	//fTextMGammaGamma = Form("mixed evt. subtr. M_{%s}",fDecayChannel.Data());

	TPad * padDataFit = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
	padDataFit->SetFillColor(0);
	padDataFit->GetFrame()->SetFillColor(0);
	padDataFit->SetBorderMode(0);
	padDataFit->Divide(fColumnPlot,fRowPlot,0.0,0.0);
	padDataFit->SetLeftMargin(0.);
	padDataFit->SetRightMargin(0.);
	padDataFit->SetTopMargin(0.);
	padDataFit->SetBottomMargin(0.);
	padDataFit->Draw();

	cout << fMesonType.Data() << endl;
	Int_t place = 0;
    for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
		Double_t startPt = fRangeBinsPt[iPt];
		Double_t endPt = fRangeBinsPt[iPt+1];

		place = place + 1;						//give the right place in the page
		if (place == fColumnPlot){
			iPt--;
			padDataFit->cd(place);

			TString textAlice = "ALICE performance";
			TString textEvents;
			if(fMonteCarloInfo){textEvents = "MC";}
			else {textEvents = "Data";}
			Double_t nPixels = 13;
			Double_t textHeight = 0.08;
			if (padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) < padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1())){
				textHeight = (Double_t)nPixels/padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) ;
			} else {
				textHeight = (Double_t)nPixels/padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1());
			}	
			Double_t startTextX = 0.10;
			Double_t startTextY = 0.9;
			Double_t differenceText = textHeight*1.25;

			TLatex *alice = 		new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
			TLatex *latexDate = 	new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
			TLatex *energy = 		new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
			TLatex *process = 		new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
			TLatex *detprocess = 	new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
			TLatex *events = 		new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight*1.3);
			alice->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
	
			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();

			process->SetNDC(); 
			process->SetTextColor(1);
			process->SetTextSize(textHeight);
			process->Draw();

			detprocess->SetNDC(); 
			detprocess->SetTextColor(1);
			detprocess->SetTextSize(textHeight);
			detprocess->Draw();

			events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();

			if (fMonteCarloInfo) {
				Double_t totalheightLeg = 2.;
				if (fTextMCvalidated.CompareTo("") != 0 && labelData){
					totalheightLeg = 3.;
				}
				TLegend* legendMC = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+totalheightLeg)*differenceText);
				legendMC->SetTextSize(textHeight);
				legendMC->SetTextFont(62);
				legendMC->SetLineColor(0);
				legendMC->SetLineWidth(0);
				legendMC->SetFillStyle(0);
				legendMC->SetFillColor(0);
				legendMC->SetMargin(0.15);
				if (labelData){
					Size_t markersize = fHistoMappingSignalInvMassPtBinPlot[iPt-1]->GetMarkerSize();
					fHistoMappingSignalInvMassPtBinPlot[iPt-1]->SetMarkerSize(3*markersize);
					legendMC->AddEntry(fHistoMappingSignalInvMassPtBinPlot[iPt-1],fTextMGammaGamma.Data(),"ep");
				}
				if (fTextMCvalidated.CompareTo("") != 0){
					Size_t markersize2 = fHistoMappingTrueMesonInvMassPtBinsPlot[iPt-1]->GetMarkerSize();
					if (labelData)fHistoMappingTrueMesonInvMassPtBinsPlot[iPt-1]->SetMarkerSize(3*markersize2);
					legendMC->AddEntry(fHistoMappingTrueMesonInvMassPtBinsPlot[iPt-1],fTextMCvalidated.Data(),"ep");
				}	
				if (fFitSignalInvMassPtBinPlot[iPt-1]!=0x00){
					Size_t linewidth = fFitSignalInvMassPtBinPlot[iPt-1]->GetLineWidth();
					fFitSignalInvMassPtBinPlot[iPt-1]->SetLineWidth(5*linewidth);
					legendMC->AddEntry(fFitSignalInvMassPtBinPlot[iPt-1],fTextFit.Data(),"l");
				}
				legendMC->Draw();				
			}else {
				TLegend* legendData = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+2.)*differenceText);
				legendData->SetTextSize(textHeight);
				legendData->SetTextFont(62);
				legendData->SetFillColor(0);
				legendData->SetFillStyle(0);
				legendData->SetLineWidth(0);
				legendData->SetLineColor(0);
				legendData->SetMargin(0.15);
				Size_t markersize = fHistoMappingSignalInvMassPtBinPlot[iPt-1]->GetMarkerSize();
				fHistoMappingSignalInvMassPtBinPlot[iPt-1]->SetMarkerSize(3*markersize);
				legendData->AddEntry(fHistoMappingSignalInvMassPtBinPlot[iPt-1],fTextMGammaGamma.Data(),"ep");
				if (fFitSignalInvMassPtBinPlot[iPt-1]!=0x00){
					Size_t linewidth = fFitSignalInvMassPtBinPlot[iPt-1]->GetLineWidth();
					fFitSignalInvMassPtBinPlot[iPt-1]->SetLineWidth(5*linewidth);
					legendData->AddEntry(fFitSignalInvMassPtBinPlot[iPt-1],fTextFit.Data(),"l");
				}	
				legendData->Draw();				
			}	

		} else {

			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.02);
			padDataFit->cd(place)->SetLeftMargin(0.15);

			cout << startPt << "-" << endPt <<endl;

			if (labelData) {
				DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot[iPt],
								Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
								Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
								fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
				DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot[iPt],
							Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
							Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
							fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);

				
				
				if(fMonteCarloInfo){
					fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetMarkerColor(kRed+2);	
					fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineColor(kRed+2);
					fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineWidth(0.99);
					fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->DrawCopy("same");
				}
			} else {
				DrawGammaHisto( fHistoMappingTrueMesonInvMassPtBinsPlot[iPt],
								Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
								Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
								fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
				DrawGammaHisto( fHistoMappingTrueMesonInvMassPtBinsPlot[iPt],
							Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
							Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
							fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);

				fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetMarkerColor(kRed+2);	
				fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetMarkerSize(0.5);	
				fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineColor(kRed+2);
// 				fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineWidth(0.99);
				fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->DrawCopy("same,p,e1");
				
			}	
			if (fFitSignalInvMassPtBinPlot[iPt]!=0x00){
				fFitSignalInvMassPtBinPlot[iPt]->SetLineColor(kCyan+3);
				fFitSignalInvMassPtBinPlot[iPt]->SetLineWidth(0.99);
				fFitSignalInvMassPtBinPlot[iPt]->DrawCopy("same");
				TLine * l1 = new TLine (fFitSignalInvMassPtBinPlot[iPt]->GetParameter(1),fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),fFitSignalInvMassPtBinPlot[iPt]->GetParameter(1),fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum());
				l1->SetLineColor(4);
				l1->SetLineWidth(0.7);
				//l1->Draw("same");
			} 

		}
	}
	canvasDataFit->Print(namePlot.Data());
	delete padDataFit;
	delete canvasDataFit;
}

//____________________________ Plotting Invariant Mass Subtracted for all bins ________________________________________________________
void PlotWithBGFitSubtractedInvMassInPtBins(TH1D ** fHistoMappingSignalPlusBGInvMassPtBinPlot, TH1D ** fHistoMappingBG, TH1D** fHistoMappingSignal, TF1 ** fFitBGInvMassPtBinPlot, TString namePlot, TString nameCanvas, TString namePad, Double_t* fPlottingRangeMeson, TString dateDummy, TString fMesonType,  Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, Double_t* fRangeBinsPt, TString fDecayChannel, Bool_t fMonteCarloInfo, TString decayChannel, TString fDetectionChannel, TString fEnergy){
    TCanvas *canvasDataFit = new TCanvas(nameCanvas.Data(),"",1400,800);  // gives the page size
	canvasDataFit->SetTopMargin(0.00);
	canvasDataFit->SetBottomMargin(0.00);
	canvasDataFit->SetRightMargin(0.00);
	canvasDataFit->SetLeftMargin(0.00);

	TPad * padDataFit = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
	padDataFit->SetFillColor(0);
	padDataFit->GetFrame()->SetFillColor(0);
	padDataFit->SetBorderMode(0);
	padDataFit->Divide(fColumnPlot,fRowPlot,0.0,0.0);
	padDataFit->SetLeftMargin(0.);
	padDataFit->SetRightMargin(0.);
	padDataFit->SetTopMargin(0.);
	padDataFit->SetBottomMargin(0.);
	padDataFit->Draw();

	cout << fMesonType.Data() << endl;
	Int_t place = 0;
	for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
		Double_t startPt = fRangeBinsPt[iPt];
		Double_t endPt = fRangeBinsPt[iPt+1];

		place = place + 1;						//give the right place in the page
		if (place == fColumnPlot){
			iPt--;
			padDataFit->cd(place);

			TString textAlice = "ALICE performance";
			TString textEvents;
			if(fMonteCarloInfo){textEvents = "MC";}
			else {textEvents = "Data";}
			Double_t nPixels = 13;
			Double_t textHeight = 0.08;
			if (padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) < padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1())){
				textHeight = (Double_t)nPixels/padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) ;
			} else {
				textHeight = (Double_t)nPixels/padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1());
			}	
			Double_t startTextX = 0.10;
			Double_t startTextY = 0.9;
			Double_t differenceText = textHeight*1.25;

			TLatex *alice = 		new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
			TLatex *latexDate = 	new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
			TLatex *energy = 		new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
			TLatex *process = 		new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
			TLatex *detprocess = 	new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
			TLatex *events = 		new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight*1.3);
			alice->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
	
			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();

			process->SetNDC(); 
			process->SetTextColor(1);
			process->SetTextSize(textHeight);
			process->Draw();

			detprocess->SetNDC(); 
			detprocess->SetTextColor(1);
			detprocess->SetTextSize(textHeight);
			detprocess->Draw();

			events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();

			TLegend* legendData = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+4.)*differenceText);
			legendData->SetTextSize(textHeight);
			legendData->SetTextFont(62);
			legendData->SetFillColor(0);
			legendData->SetFillStyle(0);
			legendData->SetLineWidth(0);
			legendData->SetLineColor(0);
			legendData->SetMargin(0.15);
			Size_t markersize = fHistoMappingSignalPlusBGInvMassPtBinPlot[iPt-1]->GetMarkerSize();
			fHistoMappingSignalPlusBGInvMassPtBinPlot[iPt-1]->SetMarkerSize(3*markersize);
			legendData->AddEntry(fHistoMappingSignalPlusBGInvMassPtBinPlot[iPt-1],Form("mixed evt. subtr. M_{%s}",decayChannel.Data()),"ep");
			fHistoMappingBG[iPt-1]->SetMarkerSize(3*markersize);
			legendData->AddEntry(fHistoMappingBG[iPt-1],"fitted add. BG","ep");
			fHistoMappingSignal[iPt-1]->SetMarkerSize(5*markersize);
			legendData->AddEntry(fHistoMappingSignal[iPt-1],"signal","ep");
			if (fFitBGInvMassPtBinPlot[iPt-1]!=0x00){
				Size_t linewidth = fFitBGInvMassPtBinPlot[iPt-1]->GetLineWidth();
				fFitBGInvMassPtBinPlot[iPt-1]->SetLineWidth(5*linewidth);
				legendData->AddEntry(fFitBGInvMassPtBinPlot[iPt-1],"BG fit","l");
			}	
			legendData->Draw();				
	
		} else {

			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.02);
			padDataFit->cd(place)->SetLeftMargin(0.15);

			cout << startPt << "-" << endPt <<endl;

			DrawGammaHisto( fHistoMappingSignalPlusBGInvMassPtBinPlot[iPt],
							Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
							Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
							fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
			DrawGammaHisto( fHistoMappingSignalPlusBGInvMassPtBinPlot[iPt],
						Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
						Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
						fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);

			fHistoMappingBG[iPt]->SetMarkerSize(0.2);
			fHistoMappingBG[iPt]->SetMarkerColor(kCyan+1);	
			fHistoMappingBG[iPt]->SetLineColor(kCyan+1);
			fHistoMappingBG[iPt]->SetLineWidth(0.99);
			fHistoMappingBG[iPt]->DrawCopy("same,pe1");				
			
			fHistoMappingSignal[iPt]->SetMarkerSize(0.2);
			fHistoMappingSignal[iPt]->SetMarkerColor(kRed+2);	
			fHistoMappingSignal[iPt]->SetLineColor(kRed+2);
			fHistoMappingSignal[iPt]->SetLineWidth(0.99);
			fHistoMappingSignal[iPt]->DrawCopy("same,pe1");

			if (fFitBGInvMassPtBinPlot[iPt]!=0x00){
				fFitBGInvMassPtBinPlot[iPt]->SetLineColor(kCyan+3);
				fFitBGInvMassPtBinPlot[iPt]->SetLineWidth(0.99);
				fFitBGInvMassPtBinPlot[iPt]->DrawCopy("same");
				TLine * l1 = new TLine (fFitBGInvMassPtBinPlot[iPt]->GetParameter(1),fHistoMappingSignalPlusBGInvMassPtBinPlot[iPt]->GetMinimum(),fFitBGInvMassPtBinPlot[iPt]->GetParameter(1),fHistoMappingSignalPlusBGInvMassPtBinPlot[iPt]->GetMaximum());
				l1->SetLineColor(4);
				l1->SetLineWidth(0.7);
				//l1->Draw("same");
			} 

		}
	}
	canvasDataFit->Print(namePlot.Data());
	delete padDataFit;
	delete canvasDataFit;
}



//______________________________ Invariant Mass with Fit to Peak Position for all Pt bins __________________________________
void PlotWithFitPeakPosInvMassInPtBins(TH1D ** fHistoMappingSignalInvMassPtBinPlot,TH1D **fHistoMappingTrueMesonInvMassPtBinsPlot, TF1 ** fFitSignalInvMassPtBinPlot, TString namePlot, TString nameCanvas, TString namePad, Double_t* fPlottingRangeMeson, TString dateDummy, TString fMesonType,  Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, Double_t* fRangeBinsPt, TString fDecayChannel, Bool_t fMonteCarloInfo, TString decayChannel, TString fDetectionChannel, TString fEnergy){
        TCanvas *canvasDataFit = new TCanvas(nameCanvas.Data(),"",1400,800);  // gives the page size
		canvasDataFit->SetTopMargin(0.02);
		canvasDataFit->SetBottomMargin(0.02);
		canvasDataFit->SetRightMargin(0.02);
		canvasDataFit->SetLeftMargin(0.02);

		TPad * padDataFit = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
		padDataFit->SetFillColor(0);
		padDataFit->GetFrame()->SetFillColor(0);
		padDataFit->SetBorderMode(0);
		padDataFit->Divide(fColumnPlot,fRowPlot,0.0,0.0);
		padDataFit->SetLeftMargin(0.2);
		padDataFit->Draw();

		cout << fMesonType.Data() << endl;
		Int_t place = 0;
		for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
			Double_t startPt = fRangeBinsPt[iPt];
			Double_t endPt = fRangeBinsPt[iPt+1];

			place = place + 1;						//give the right place in the page
			if (place == fColumnPlot){
				iPt--;
				padDataFit->cd(place);

				TString textAlice = "ALICE performance";
				TString textEvents;
				if(fMonteCarloInfo){textEvents = "MC";}
				else {textEvents = "Data";}
				Double_t nPixels = 13;
				Double_t textHeight = 0.08;
				if (padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) < padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1())){
					textHeight = (Double_t)nPixels/padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) ;
				} else {
					textHeight = (Double_t)nPixels/padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1());
				}	
				Double_t startTextX = 0.10;
				Double_t startTextY = 0.8;
				Double_t differenceText = textHeight*1.25;

				TLatex *alice = 		new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
				TLatex *latexDate = 	new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
				TLatex *energy = 		new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
				TLatex *process = 		new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
				TLatex *detprocess = 	new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
				TLatex *events = 		new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

				alice->SetNDC();
				alice->SetTextColor(1);
				alice->SetTextSize(textHeight*1.3);
				alice->Draw();

				latexDate->SetNDC();
				latexDate->SetTextColor(1);
				latexDate->SetTextSize(textHeight);
				latexDate->Draw();
		
				energy->SetNDC();
				energy->SetTextColor(1);
				energy->SetTextSize(textHeight);
				energy->Draw();

				process->SetNDC(); 
				process->SetTextColor(1);
				process->SetTextSize(textHeight);
				process->Draw();

				detprocess->SetNDC(); 
				detprocess->SetTextColor(1);
				detprocess->SetTextSize(textHeight);
				detprocess->Draw();

				events->SetNDC();
				events->SetTextColor(1);
				events->SetTextSize(textHeight);
				events->Draw();
			} else {

				padDataFit->cd(place);
				padDataFit->cd(place)->SetTopMargin(0.12);
				padDataFit->cd(place)->SetBottomMargin(0.15);
				padDataFit->cd(place)->SetRightMargin(0.02);
				padDataFit->cd(place)->SetLeftMargin(0.15);

				cout << startPt << "-" << endPt <<endl;

				DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot[iPt],
							 Form("%3.2f GeV/c < E < %3.2f GeV/c",startPt,endPt),
							 Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
							 fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
				DrawGammaHisto( fHistoMappingSignalInvMassPtBinPlot[iPt],
							Form("%3.2f GeV/c < E < %3.2f GeV/c",startPt,endPt),
							Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
							fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);

				if(fMonteCarloInfo){
					fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineColor(kRed);
					fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->SetLineWidth(0.99);
					fHistoMappingTrueMesonInvMassPtBinsPlot[iPt]->DrawCopy("same");
				}
				if (fFitSignalInvMassPtBinPlot[iPt]!=0x00){
					fFitSignalInvMassPtBinPlot[iPt]->SetLineColor(kCyan+3);
					fFitSignalInvMassPtBinPlot[iPt]->SetLineWidth(0.99);
					fFitSignalInvMassPtBinPlot[iPt]->DrawCopy("same");
					TLine * l1 = new TLine (fFitSignalInvMassPtBinPlot[iPt]->GetParameter(1),fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMinimum(),fFitSignalInvMassPtBinPlot[iPt]->GetParameter(1),fHistoMappingSignalInvMassPtBinPlot[iPt]->GetMaximum());
					l1->SetLineColor(4);
					l1->SetLineWidth(0.7);
					//l1->Draw("same");
				} 

		}
	}
	canvasDataFit->Print(namePlot.Data());
	delete padDataFit;
	delete canvasDataFit;
}



//____________________________ Plotting Invariant Mass validated splitted in different categories for all bins ________________________________________________________
void PlotTrueInvMassSplittedInPhotonAndElectronInPtBins(TH1D ** fHistoTrueSignal, TH1D** fHistoTrueSignalPhotons, TH1D** fHistoTrueSignalElectrons, TH1D** fHistoTrueSignalConvPhotons, TH1D** fHistoTrueSignalMixed, TString namePlot, TString nameCanvas, TString namePad, Double_t* fPlottingRangeMeson, TString dateDummy, TString fMesonType,  Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, Double_t* fRangeBinsPt, TString fDecayChannel, Bool_t fMonteCarloInfo, TString decayChannel, TString fDetectionChannel, TString fEnergy, Int_t mode){
    TCanvas *canvasDataFit = new TCanvas(nameCanvas.Data(),"",1400,800);  // gives the page size
	canvasDataFit->SetTopMargin(0.00);
	canvasDataFit->SetBottomMargin(0.00);
	canvasDataFit->SetRightMargin(0.00);
	canvasDataFit->SetLeftMargin(0.00);

	TPad * padDataFit = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
	padDataFit->SetFillColor(0);
	padDataFit->GetFrame()->SetFillColor(0);
	padDataFit->SetBorderMode(0);
	padDataFit->Divide(fColumnPlot,fRowPlot,0.0,0.0);
	padDataFit->SetLeftMargin(0.);
	padDataFit->SetRightMargin(0.);
	padDataFit->SetTopMargin(0.);
	padDataFit->SetBottomMargin(0.);
	padDataFit->Draw();

	cout << fMesonType.Data() << endl;
	Int_t place = 0;
	for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
		Double_t startPt = fRangeBinsPt[iPt];
		Double_t endPt = fRangeBinsPt[iPt+1];

		place = place + 1;						//give the right place in the page
		if (place == fColumnPlot){
			iPt--;
			padDataFit->cd(place);

			TString textAlice = "ALICE performance";
			TString textEvents;
			if(fMonteCarloInfo){textEvents = "MC";}
			else {textEvents = "Data";}
			Double_t nPixels = 10;
			Double_t textHeight = 0.08;
			if (padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) < padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1())){
				textHeight = (Double_t)nPixels/padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) ;
			} else {
				textHeight = (Double_t)nPixels/padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1());
			}	
			Double_t startTextX = 0.10;
			Double_t startTextY = 0.9;
			Double_t differenceText = textHeight*1.25;

			TLatex *alice = 		new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
			TLatex *latexDate = 	new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
			TLatex *energy = 		new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
			TLatex *process = 		new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
			TLatex *detprocess = 	new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
			TLatex *events = 		new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight*1.3);
			alice->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
	
			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();

			process->SetNDC(); 
			process->SetTextColor(1);
			process->SetTextSize(textHeight);
			process->Draw();

			detprocess->SetNDC(); 
			detprocess->SetTextColor(1);
			detprocess->SetTextSize(textHeight);
			detprocess->Draw();

			events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();
			Double_t nSignals = 4;
			if (fHistoTrueSignalMixed != NULL) nSignals=5;
			TLegend* legendMC = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+nSignals)*differenceText);
			legendMC->SetTextSize(textHeight);
			legendMC->SetTextFont(62);
			legendMC->SetLineColor(0);
			legendMC->SetLineWidth(0);
			legendMC->SetFillStyle(0);
			legendMC->SetFillColor(0);
			legendMC->SetMargin(0.15);
			Size_t markersize2 = fHistoTrueSignal[iPt-1]->GetMarkerSize();
			fHistoTrueSignal[iPt-1]->SetMarkerSize(3*markersize2);
			fHistoTrueSignalPhotons[iPt-1]->SetMarkerSize(3*markersize2);
			fHistoTrueSignalElectrons[iPt-1]->SetMarkerSize(3*markersize2);
			fHistoTrueSignalConvPhotons[iPt-1]->SetMarkerSize(3*markersize2);
			if (fHistoTrueSignalMixed != NULL)fHistoTrueSignalMixed[iPt-1]->SetMarkerSize(3*markersize2);
			legendMC->AddEntry(fHistoTrueSignal[iPt-1],"validated meson","ep");
			if (mode == 4 || mode == 5){
				legendMC->AddEntry(fHistoTrueSignalPhotons[iPt-1],"val. #gamma#gamma","ep");
				legendMC->AddEntry(fHistoTrueSignalElectrons[iPt-1],"val. e^{#pm}e^{#pm}","ep");
				legendMC->AddEntry(fHistoTrueSignalConvPhotons[iPt-1],"val. #gamma_{conv}#gamma_{conv}","ep");
				if (fHistoTrueSignalMixed != NULL) legendMC->AddEntry(fHistoTrueSignalMixed[iPt-1],"val. #gamma#gamma_{conv}","ep");
			} else if (mode == 2 || mode == 3) {
				legendMC->AddEntry(fHistoTrueSignalPhotons[iPt-1],"val. #gamma_{conv}#gamma","ep");
				legendMC->AddEntry(fHistoTrueSignalElectrons[iPt-1],"val. #gamma_{conv}e^{#pm}","ep");
				legendMC->AddEntry(fHistoTrueSignalConvPhotons[iPt-1],"val. #gamma_{conv}#gamma_{conv}","ep");
			}	
			legendMC->Draw();				
			
		} else {

			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.02);
			padDataFit->cd(place)->SetLeftMargin(0.15);

			cout << startPt << "-" << endPt <<endl;

			DrawGammaHisto( fHistoTrueSignal[iPt],
							Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
							Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
							fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
			DrawGammaHisto( fHistoTrueSignal[iPt],
						Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
						Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
						fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);

			fHistoTrueSignalPhotons[iPt]->SetMarkerColor(kRed+2);
			fHistoTrueSignalPhotons[iPt]->SetMarkerStyle(21);	
			fHistoTrueSignalPhotons[iPt]->SetMarkerSize(0.3);	
			fHistoTrueSignalPhotons[iPt]->SetLineColor(kRed+2);
			fHistoTrueSignalPhotons[iPt]->DrawCopy("same,p,e1");
				
			fHistoTrueSignalElectrons[iPt]->SetMarkerColor(kBlue+2);
			fHistoTrueSignalElectrons[iPt]->SetMarkerStyle(33);	
			fHistoTrueSignalElectrons[iPt]->SetMarkerSize(0.5);	
			fHistoTrueSignalElectrons[iPt]->SetLineColor(kBlue+2);
			fHistoTrueSignalElectrons[iPt]->DrawCopy("same,p,e1");

			fHistoTrueSignalConvPhotons[iPt]->SetMarkerColor(kCyan+2);
			fHistoTrueSignalConvPhotons[iPt]->SetMarkerStyle(27);	
			fHistoTrueSignalConvPhotons[iPt]->SetMarkerSize(0.3);	
			fHistoTrueSignalConvPhotons[iPt]->SetLineColor(kCyan+2);
			fHistoTrueSignalConvPhotons[iPt]->DrawCopy("same,p,e1");
			
			if (fHistoTrueSignalMixed != NULL){
				fHistoTrueSignalMixed[iPt]->SetMarkerColor(kViolet+2);
				fHistoTrueSignalMixed[iPt]->SetMarkerStyle(24);	
				fHistoTrueSignalMixed[iPt]->SetMarkerSize(0.3);	
				fHistoTrueSignalMixed[iPt]->SetLineColor(kViolet+2);
				fHistoTrueSignalMixed[iPt]->DrawCopy("same,p,e1");
			}	
			
		}
	}
	canvasDataFit->Print(namePlot.Data());
	delete padDataFit;
	delete canvasDataFit;
}

//____________________________ Plotting Invariant Mass validated splitted in different categories for all bins ________________________________________________________
void PlotTrueInvMassSplittedInMergedInPtBins(TH1D ** fHistoTrueSignal, TH1D** fHistoTrueSignalMerged, TH1D** fHistoTrueSignalMergedPartConv, TString namePlot, TString nameCanvas, TString namePad, Double_t* fPlottingRangeMeson, TString dateDummy, TString fMesonType,  Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, Double_t* fRangeBinsPt, TString fDecayChannel, Bool_t fMonteCarloInfo, TString decayChannel, TString fDetectionChannel, TString fEnergy, Int_t mode){
    TCanvas *canvasDataFit = new TCanvas(nameCanvas.Data(),"",1400,800);  // gives the page size
	canvasDataFit->SetTopMargin(0.00);
	canvasDataFit->SetBottomMargin(0.00);
	canvasDataFit->SetRightMargin(0.00);
	canvasDataFit->SetLeftMargin(0.00);

	TPad * padDataFit = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
	padDataFit->SetFillColor(0);
	padDataFit->GetFrame()->SetFillColor(0);
	padDataFit->SetBorderMode(0);
	padDataFit->Divide(fColumnPlot,fRowPlot,0.0,0.0);
	padDataFit->SetLeftMargin(0.);
	padDataFit->SetRightMargin(0.);
	padDataFit->SetTopMargin(0.);
	padDataFit->SetBottomMargin(0.);
	padDataFit->Draw();

	cout << fMesonType.Data() << endl;
	Int_t place = 0;
	for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
		Double_t startPt = fRangeBinsPt[iPt];
		Double_t endPt = fRangeBinsPt[iPt+1];

		place = place + 1;						//give the right place in the page
		if (place == fColumnPlot){
			iPt--;
			padDataFit->cd(place);

			TString textAlice = "ALICE performance";
			TString textEvents;
			if(fMonteCarloInfo){textEvents = "MC";}
			else {textEvents = "Data";}
			Double_t nPixels = 10;
			Double_t textHeight = 0.08;
			if (padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) < padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1())){
				textHeight = (Double_t)nPixels/padDataFit->cd(place)->XtoPixel(padDataFit->cd(place)->GetX2()) ;
			} else {
				textHeight = (Double_t)nPixels/padDataFit->cd(place)->YtoPixel(padDataFit->cd(place)->GetY1());
			}	
			Double_t startTextX = 0.10;
			Double_t startTextY = 0.9;
			Double_t differenceText = textHeight*1.25;

			TLatex *alice = 		new TLatex(startTextX, startTextY, Form("%s",textAlice.Data()));
			TLatex *latexDate = 	new TLatex(startTextX, (startTextY-1.25*differenceText), dateDummy.Data());
			TLatex *energy = 		new TLatex(startTextX, (startTextY-2.25*differenceText), fEnergy);
			TLatex *process = 		new TLatex(startTextX, (startTextY-3.25*differenceText), fDecayChannel);
			TLatex *detprocess = 	new TLatex(startTextX, (startTextY-4.25*differenceText), fDetectionChannel);
			TLatex *events = 		new TLatex(startTextX, (startTextY-5.25*differenceText), Form("%s: %2.1e events",textEvents.Data(), fNEvents));

			alice->SetNDC();
			alice->SetTextColor(1);
			alice->SetTextSize(textHeight*1.3);
			alice->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
	
			energy->SetNDC();
			energy->SetTextColor(1);
			energy->SetTextSize(textHeight);
			energy->Draw();

			process->SetNDC(); 
			process->SetTextColor(1);
			process->SetTextSize(textHeight);
			process->Draw();

			detprocess->SetNDC(); 
			detprocess->SetTextColor(1);
			detprocess->SetTextSize(textHeight);
			detprocess->Draw();

			events->SetNDC();
			events->SetTextColor(1);
			events->SetTextSize(textHeight);
			events->Draw();
			Double_t nSignals = 3;
			TLegend* legendMC = new TLegend(startTextX,startTextY-5.75*differenceText,1,startTextY-(5.75+nSignals)*differenceText);
			legendMC->SetTextSize(textHeight);
			legendMC->SetTextFont(62);
			legendMC->SetLineColor(0);
			legendMC->SetLineWidth(0);
			legendMC->SetFillStyle(0);
			legendMC->SetFillColor(0);
			legendMC->SetMargin(0.15);
			Size_t markersize2 = fHistoTrueSignal[iPt-1]->GetMarkerSize();
			fHistoTrueSignal[iPt-1]->SetMarkerSize(3*markersize2);
			fHistoTrueSignalMerged[iPt-1]->SetMarkerSize(3*markersize2);
			fHistoTrueSignalMergedPartConv[iPt-1]->SetMarkerSize(3*markersize2);
			legendMC->AddEntry(fHistoTrueSignal[iPt-1],"validated meson","ep");
			if (mode == 4 || mode == 5){
				legendMC->AddEntry(fHistoTrueSignalMerged[iPt-1],"at least 1 EM_{merged}","ep");
				legendMC->AddEntry(fHistoTrueSignalMergedPartConv[iPt-1],"at least 1 EM_{merged}, part conv","ep");
			} else {
				legendMC->AddEntry(fHistoTrueSignalMerged[iPt-1],"val. #gamma_{conv}EM_{merged}","ep");
				legendMC->AddEntry(fHistoTrueSignalMergedPartConv[iPt-1],"val. #gamma_{conv}EM_{merged}, part conv","ep");
			}	
			legendMC->Draw();				
			
		} else {

			padDataFit->cd(place);
			padDataFit->cd(place)->SetTopMargin(0.12);
			padDataFit->cd(place)->SetBottomMargin(0.15);
			padDataFit->cd(place)->SetRightMargin(0.02);
			padDataFit->cd(place)->SetLeftMargin(0.15);

			cout << startPt << "-" << endPt <<endl;

			DrawGammaHisto( fHistoTrueSignal[iPt],
							Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
							Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
							fPlottingRangeMeson[0],fPlottingRangeMeson[1],0);
			DrawGammaHisto( fHistoTrueSignal[iPt],
						Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt),
						Form("M_{%s} (GeV/c^{2})",decayChannel.Data()), Form("dN_{%s}/dM_{%s}",decayChannel.Data(), decayChannel.Data()),
						fPlottingRangeMeson[0],fPlottingRangeMeson[1],2);

			fHistoTrueSignalMerged[iPt]->SetMarkerColor(kRed+2);
			fHistoTrueSignalMerged[iPt]->SetMarkerStyle(21);	
			fHistoTrueSignalMerged[iPt]->SetMarkerSize(0.3);	
			fHistoTrueSignalMerged[iPt]->SetLineColor(kRed+2);
			fHistoTrueSignalMerged[iPt]->DrawCopy("same,p,e1");
				
			fHistoTrueSignalMergedPartConv[iPt]->SetMarkerColor(kBlue+2);
			fHistoTrueSignalMergedPartConv[iPt]->SetMarkerStyle(33);	
			fHistoTrueSignalMergedPartConv[iPt]->SetMarkerSize(0.5);	
			fHistoTrueSignalMergedPartConv[iPt]->SetLineColor(kBlue+2);
			fHistoTrueSignalMergedPartConv[iPt]->DrawCopy("same,p,e1");
			
		}
	}
	canvasDataFit->Print(namePlot.Data());
	delete padDataFit;
	delete canvasDataFit;
}
